/* Bitonic sort cost bench for distinguished-point records (spec section 9). Standalone: it does not
 * include step.c and shares no code with it.
 *
 * Usage: sortnet <inputs> <rounds> <sortnet|sortscan> [seed [n]]
 *   sortnet:  <rounds> ascending bitonic sorts of n records (fingerprint key, 128 bits).
 *   sortscan: the same, with the adjacent-equality scan after each sort.
 *   seed (default 1) seeds the xorshift64 fill; n (default 4096) is a power of two.
 * The inputs file is accepted for the organizer harness (measure.sh) and is not read. Any other
 * mode (the harness's default walk and dump runs) returns at once. Counts go to stderr. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>

typedef struct { uint64_t v[4]; } u256;
typedef struct { uint64_t fp[2]; u256 a, b; } dp_t;

/* Constant-instruction comparison of 128-bit fingerprints (fp[0] is the low word): the fp_lt_mask
 * limb loop over two limbs. Returns all ones when x < y. */
static inline uint64_t fp128_lt_mask(const uint64_t x[2], const uint64_t y[2]) {
    uint64_t lt = 0;
    for (int i = 0; i < 2; i++) {
        uint64_t ltw = (uint64_t)(x[i] < y[i]);
        uint64_t ne = ltw | (uint64_t)(x[i] > y[i]);
        uint64_t m = 0 - ne;
        lt = (ltw & m) | (lt & ~m);
    }
    return 0 - lt;
}

/* Exchange x and y when m is all ones: XOR masking over all ten words of the record. */
static inline void cswap_rec(dp_t *x, dp_t *y, uint64_t m) {
    for (int i = 0; i < 2; i++) {
        uint64_t t = (x->fp[i] ^ y->fp[i]) & m;
        x->fp[i] ^= t;
        y->fp[i] ^= t;
    }
    for (int i = 0; i < 4; i++) {
        uint64_t t = (x->a.v[i] ^ y->a.v[i]) & m;
        x->a.v[i] ^= t;
        y->a.v[i] ^= t;
    }
    for (int i = 0; i < 4; i++) {
        uint64_t t = (x->b.v[i] ^ y->b.v[i]) & m;
        x->b.v[i] ^= t;
        y->b.v[i] ^= t;
    }
}

/* Bitonic sorting network, ascending on the fingerprint. Pairs and directions depend on the indices
 * only. The exchange decision is (fp[i] < fp[l]) xor asc: an ascending comparator exchanges when
 * fp[i] >= fp[l], so equal keys may be exchanged, which leaves the key sequence unchanged. */
__attribute__((noinline)) static void bitonic_sort(dp_t *R, uint32_t n) {
    for (uint32_t k = 2; k <= n; k <<= 1) {
        for (uint32_t j = k >> 1; j > 0; j >>= 1) {
            for (uint32_t i = 0; i < n; i++) {
                uint32_t l = i ^ j;
                if (l > i) {
                    uint64_t asc = 0 - (uint64_t)((i & k) == 0);
                    uint64_t m = fp128_lt_mask(R[i].fp, R[l].fp) ^ asc;
                    cswap_rec(&R[i], &R[l], m);
                }
            }
        }
    }
}

/* Number of adjacent record pairs with equal fingerprints, branch-free. */
__attribute__((noinline)) static uint64_t fp_adjacent_equal(const dp_t *R, uint32_t n) {
    uint64_t acc = 0;
    for (uint32_t i = 0; i + 1 < n; i++) {
        uint64_t d = (R[i].fp[0] ^ R[i + 1].fp[0]) | (R[i].fp[1] ^ R[i + 1].fp[1]);
        acc += ((d | (0 - d)) >> 63) ^ 1;
    }
    return acc;
}

static uint64_t xs64(uint64_t *s) {
    uint64_t x = *s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *s = x;
    return x;
}

/* The records come from xorshift64 seeded by `seed`; they are filled before collection starts. */
__attribute__((noinline)) static void sort_mode(uint64_t rounds, uint32_t n, uint64_t seed, int with_scan) {
    dp_t *R = malloc(sizeof(dp_t) * n);
    uint64_t s = (seed + 1) * 0x9E3779B97F4A7C15ULL;
    if (s == 0) s = 1;
    for (uint32_t i = 0; i < n; i++) {
        dp_t *r = &R[i];
        r->fp[0] = xs64(&s);
        r->fp[1] = xs64(&s);
        for (int j = 0; j < 4; j++) r->a.v[j] = xs64(&s);
        for (int j = 0; j < 4; j++) r->b.v[j] = xs64(&s);
    }
    volatile uint64_t sink = 0;
#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    uint64_t acc = 0;
    if (with_scan) {
        for (uint64_t r = 0; r < rounds; r++) {
            bitonic_sort(R, n);
            acc += fp_adjacent_equal(R, n);
        }
    } else {
        for (uint64_t r = 0; r < rounds; r++) bitonic_sort(R, n);
    }
    acc ^= R[0].fp[0] ^ R[n - 1].fp[1];
    sink = acc;
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif
    free(R);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: sortnet <inputs> <rounds> <sortnet|sortscan> [seed [n]]\n");
        return 2;
    }
    uint64_t rounds = strtoull(argv[2], NULL, 10);
    const char *mode = argc > 3 ? argv[3] : "";
    int with_scan = strcmp(mode, "sortscan") == 0;
    if (!with_scan && strcmp(mode, "sortnet") != 0) return 0;
    uint64_t seed = argc > 4 ? strtoull(argv[4], NULL, 10) : 1;
    uint64_t n = argc > 5 ? strtoull(argv[5], NULL, 10) : 4096;
    if (n < 2 || (n & (n - 1)) != 0 || n > (1ULL << 30)) {
        fprintf(stderr, "n must be a power of two in [2, 2^30]\n");
        return 2;
    }
    sort_mode(rounds, (uint32_t)n, seed, with_scan);
    fprintf(stderr, "%s=%llu n=%llu seed=%llu\n", mode, (unsigned long long)rounds,
            (unsigned long long)n, (unsigned long long)seed);
    return 0;
}
