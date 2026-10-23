/* 64-bit-index bench for the distinguished-point sort (spec section 9). Standalone: it does not
 * include step.c or sortnet.c.
 *
 * Usage: sortnet64 <inputs> <rounds> <mode> [seed [n [verify]]]
 *   sortnet64:  <rounds> ascending bitonic sorts of n records by the 128-bit key (fp[1] is the high word).
 *   sortscan64: the same, with the equal-neighbour scan after each sort.
 *   scatter64:  <rounds> bucket scatter passes. Record bucket b = fp[1] >> 27 (37 bits, the production
 *               bucket). The measurement keeps 2^16 cursors, indexed by b >> 21 (the top 16 bits of b),
 *               because 2^37 cursors cannot be allocated here. The per-record work (one cursor load and
 *               store, one 80-byte record copy) does not depend on the number of buckets. Each round
 *               resets the cursors (one store per cursor, counted) and scatters the n records into
 *               slots of cap = 2n / 2^16 rounded up to a power of two; cap overflows are counted.
 *   seed (default 1) seeds xorshift64; n (default 4096) is a power of two, at most 2^30.
 *   verify:     sort modes only. After the rounds (outside collection) counts the adjacent pairs out of
 *               order by the 128-bit key, and checks that the output is a permutation of the input with
 *               an order-independent checksum.
 * The inputs file is accepted for the organizer harness and is not read. Any other mode (the harness's
 * default walk and dump runs) returns at once. Counts go to stderr. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>

typedef struct { uint64_t v[4]; } u256;
typedef struct { uint64_t fp[2]; u256 a, b; } dp_t;

#define SC_BUCKET_BITS 16
#define SC_NB ((uint64_t)1 << SC_BUCKET_BITS)

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

/* Bitonic sorting network, ascending on the 128-bit fingerprint, 64-bit indices. Pairs and directions
 * depend on the indices only. The exchange decision is (fp[i] < fp[l]) xor asc: an ascending comparator
 * exchanges when fp[i] >= fp[l], so equal keys may be exchanged, which leaves the key sequence unchanged. */
__attribute__((noinline)) static void bitonic_sort64(dp_t *R, uint64_t n) {
    for (uint64_t k = 2; k <= n; k <<= 1) {
        for (uint64_t j = k >> 1; j > 0; j >>= 1) {
            for (uint64_t i = 0; i < n; i++) {
                uint64_t l = i ^ j;
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
__attribute__((noinline)) static uint64_t fp_adjacent_equal64(const dp_t *R, uint64_t n) {
    uint64_t acc = 0;
    for (uint64_t i = 0; i + 1 < n; i++) {
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

static uint64_t mix64(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

/* Order-independent record hash: the sum over records of this is unchanged by a permutation. */
static uint64_t rec_hash(const dp_t *r) {
    uint64_t h = mix64(r->fp[0] ^ 0x1) + mix64(r->fp[1] ^ 0x2);
    for (int j = 0; j < 4; j++) h += mix64(r->a.v[j] ^ (0x10 + (uint64_t)j)) + mix64(r->b.v[j] ^ (0x20 + (uint64_t)j));
    return h;
}

static uint64_t checksum(const dp_t *R, uint64_t n) {
    uint64_t c = 0;
    for (uint64_t i = 0; i < n; i++) c += rec_hash(&R[i]);
    return c;
}

static uint64_t unsorted_pairs(const dp_t *R, uint64_t n) {
    uint64_t bad = 0;
    for (uint64_t i = 0; i + 1 < n; i++) {
        uint64_t gt_hi = (uint64_t)(R[i].fp[1] > R[i + 1].fp[1]);
        uint64_t eq_hi = (uint64_t)(R[i].fp[1] == R[i + 1].fp[1]);
        uint64_t gt_lo = (uint64_t)(R[i].fp[0] > R[i + 1].fp[0]);
        bad += gt_hi | (eq_hi & gt_lo);
    }
    return bad;
}

/* The records come from xorshift64 seeded by `seed`; they are filled before collection starts. */
static dp_t *fill_records(uint64_t n, uint64_t seed) {
    dp_t *R = malloc(sizeof(dp_t) * n);
    if (!R) {
        fprintf(stderr, "allocation of %llu records failed\n", (unsigned long long)n);
        exit(2);
    }
    uint64_t s = (seed + 1) * 0x9E3779B97F4A7C15ULL;
    if (s == 0) s = 1;
    for (uint64_t i = 0; i < n; i++) {
        dp_t *r = &R[i];
        r->fp[0] = xs64(&s);
        r->fp[1] = xs64(&s);
        for (int j = 0; j < 4; j++) r->a.v[j] = xs64(&s);
        for (int j = 0; j < 4; j++) r->b.v[j] = xs64(&s);
    }
    return R;
}

__attribute__((noinline)) static void sort_mode64(uint64_t rounds, uint64_t n, uint64_t seed, int with_scan, int verify) {
    dp_t *R = fill_records(n, seed);
    uint64_t before = verify ? checksum(R, n) : 0;
    volatile uint64_t sink = 0;
#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    uint64_t acc = 0;
    if (with_scan) {
        for (uint64_t r = 0; r < rounds; r++) {
            bitonic_sort64(R, n);
            acc += fp_adjacent_equal64(R, n);
        }
    } else {
        for (uint64_t r = 0; r < rounds; r++) bitonic_sort64(R, n);
    }
    acc ^= R[0].fp[0] ^ R[n - 1].fp[1];
    sink = acc;
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif
    if (verify) {
        uint64_t bad = unsorted_pairs(R, n);
        uint64_t after = checksum(R, n);
        fprintf(stderr, "verify n=%llu rounds=%llu unsorted_pairs=%llu permutation=%s\n",
                (unsigned long long)n, (unsigned long long)rounds, (unsigned long long)bad,
                after == before ? "ok" : "BROKEN");
    }
    free(R);
}

/* One scatter pass: record i goes to bucket b = fp[1] >> 27, measured as the cursor bk = b >> 21. The slot
 * is bk * cap + (cursor mod cap); cursor >= cap counts an overflow (the slot is then overwritten; the
 * count is what is reported). Branch-free on the data. */
__attribute__((noinline)) static uint64_t scatter_pass(const dp_t *R, uint64_t n, dp_t *out, uint64_t *cursor, uint64_t cap) {
    uint64_t ovf = 0;
    for (uint64_t i = 0; i < n; i++) {
        uint64_t b = R[i].fp[1] >> 27;
        uint64_t bk = b >> 21;
        uint64_t c = cursor[bk];
        out[bk * cap + (c & (cap - 1))] = R[i];
        ovf += (uint64_t)(c >= cap);
        cursor[bk] = c + 1;
    }
    return ovf;
}

__attribute__((noinline)) static void scatter_mode64(uint64_t rounds, uint64_t n, uint64_t seed) {
    dp_t *R = fill_records(n, seed);
    uint64_t cap = 1;
    while (cap * SC_NB < 2 * n) cap <<= 1;
    dp_t *out = malloc(sizeof(dp_t) * SC_NB * cap);
    uint64_t *cursor = malloc(sizeof(uint64_t) * SC_NB);
    if (!out || !cursor) {
        fprintf(stderr, "allocation failed for the scatter buckets\n");
        exit(2);
    }
    volatile uint64_t sink = 0;
    uint64_t ovf_total = 0;
#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    for (uint64_t r = 0; r < rounds; r++) {
        for (uint64_t b = 0; b < SC_NB; b++) cursor[b] = 0;
        ovf_total += scatter_pass(R, n, out, cursor, cap);
    }
    sink = ovf_total ^ out[0].fp[0] ^ out[SC_NB * cap - 1].fp[1];
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif
    fprintf(stderr, "scatter64 overflows=%llu cap=%llu\n", (unsigned long long)ovf_total, (unsigned long long)cap);
    free(out);
    free(cursor);
    free(R);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: sortnet64 <inputs> <rounds> <sortnet64|sortscan64|scatter64> [seed [n [verify]]]\n");
        return 2;
    }
    uint64_t rounds = strtoull(argv[2], NULL, 10);
    const char *mode = argc > 3 ? argv[3] : "";
    int kind = 0;
    if (strcmp(mode, "sortnet64") == 0) kind = 1;
    else if (strcmp(mode, "sortscan64") == 0) kind = 2;
    else if (strcmp(mode, "scatter64") == 0) kind = 3;
    else return 0;
    uint64_t seed = argc > 4 ? strtoull(argv[4], NULL, 10) : 1;
    uint64_t n = argc > 5 ? strtoull(argv[5], NULL, 10) : 4096;
    int verify = argc > 6 && strcmp(argv[6], "verify") == 0;
    if (n < 2 || (n & (n - 1)) != 0 || n > (1ULL << 30)) {
        fprintf(stderr, "n must be a power of two in [2, 2^30]\n");
        return 2;
    }
    if (kind == 3) {
        scatter_mode64(rounds, n, seed);
    } else {
        sort_mode64(rounds, n, seed, kind == 2, verify);
        fprintf(stderr, "%s=%llu n=%llu seed=%llu\n", mode, (unsigned long long)rounds,
                (unsigned long long)n, (unsigned long long)seed);
    }
    return 0;
}
