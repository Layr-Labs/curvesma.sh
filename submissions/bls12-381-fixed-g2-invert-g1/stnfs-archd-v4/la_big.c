/* stnfs-archd-v2: full-size (Schirokauer-map) matrix-entry MAC for block Wiedemann
 * over F_r (r = BLS12-381 group order, 255 bits, 4 limbs).
 *
 * The raw (unfiltered) relation matrix carries, per row, w_v SMALL valuation entries
 * (|c| < 256; priced by the shipped v1 la_mac kernel at 33 retired instr/entry) plus
 * S_COUNT = 41 FULL-SIZE entries: the Schirokauer-map values sigma_j(relation) in F_r
 * (dgp2021 section 6.5: "the 7 dense columns of Schirokauer maps"; here 35 f-side +
 * 6 g-side, the unit ranks of the two absolute fields).  This kernel prices one
 * full-size entry: acc += a * b with a, b in [0, r).
 *
 * Inner loop: column-ordered schoolbook 4x4 -> 8-limb product added into a per-row
 * 10-limb lazy accumulator (W <= 2^6 entries x 2^510 < 2^516 fits).  Row finish:
 * one Barrett reduction of the 10-limb accumulator (mu = floor(2^516 / r),
 * qhat = (acc >> 258) * mu >> 258, error <= 10 r handled by compare-subtract).
 *
 * Build: clang -O2 -o la_big la_big.c
 *   ./la_big check          prints "i W result-hex" lines (verify with la_big_check.py)
 *   ./la_big bench R W N    deterministic benchmark, N repetitions (N=0: setup only)
 * Retired instr/entry = (instrs(N=8) - instrs(N=0)) / (8 * R * W), /usr/bin/time -l.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct { uint64_t v[4]; } fe_t;
typedef struct { uint64_t v[10]; } acc10_t;

static const uint64_t RMOD[4] = {0xffffffff00000001ull, 0x53bda402fffe5bfeull,
                                 0x3339d80809a1d805ull, 0x73eda753299d7d48ull};
static const uint64_t MU[5] = {0x2737a020c0d63931ull, 0x5043eb4be4bad714ull,
                               0x8b5dcb707e08ed36ull, 0x55094edfede377c3ull, 0x23ull};

static uint64_t rngs = 88172645463325252ull;
static uint64_t xs64(void) {
    rngs ^= rngs << 13; rngs ^= rngs >> 7; rngs ^= rngs << 17;
    return rngs;
}
static void fe_rand(fe_t *e) {
    for (;;) {
        uint64_t w[4] = {xs64(), xs64(), xs64(), xs64() & 0x7fffffffffffffffull};
        int ge = 0;
        for (int i = 3; i >= 0; i--) {
            if (w[i] != RMOD[i]) { ge = w[i] > RMOD[i]; break; }
        }
        if (!ge) { memcpy(e->v, w, 32); return; }
    }
}

static inline void big_mac(acc10_t *acc, const fe_t *a, const fe_t *b) {
    uint64_t pc[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            __uint128_t p = (__uint128_t)a->v[i] * b->v[j];
            uint64_t lo = (uint64_t)p, hi = (uint64_t)(p >> 64);
            uint64_t t1 = pc[i + j] + lo;
            uint64_t c1 = (t1 < lo);
            uint64_t t2 = t1 + carry;
            uint64_t c2 = (t2 < t1);
            pc[i + j] = t2;
            carry = hi + c1 + c2;
        }
        pc[i + 4] += carry;
    }
    uint64_t carry = 0;
    for (int k = 0; k < 8; k++) {
        uint64_t t1 = acc->v[k] + pc[k];
        uint64_t c1 = (t1 < pc[k]);
        uint64_t t2 = t1 + carry;
        uint64_t c2 = (t2 < t1);
        acc->v[k] = t2;
        carry = c1 + c2;
    }
    for (int k = 8; k < 10 && carry; k++) {
        uint64_t t = acc->v[k] + carry;
        carry = (t == 0);
        acc->v[k] = t;
    }
}

static int ge_r(const uint64_t x[10]) {
    for (int i = 9; i >= 4; i--)
        if (x[i]) return 1;
    for (int i = 3; i >= 0; i--)
        if (x[i] != RMOD[i]) return x[i] > RMOD[i];
    return 1;   /* equal */
}
static void sub_r(uint64_t x[10]) {
    uint64_t borrow = 0;
    for (int k = 0; k < 4; k++) {
        uint64_t t = x[k] - RMOD[k] - borrow;
        borrow = (x[k] < RMOD[k]) | ((x[k] == RMOD[k]) & (borrow != 0));
        x[k] = t;
    }
    for (int k = 4; k < 10 && borrow; k++) {
        uint64_t t = x[k] - 1;
        borrow = (x[k] == 0);
        x[k] = t;
    }
}

static void row_finish(const acc10_t *acc, fe_t *out) {
    /* t = acc >> 258 (258 = 4*64 + 2): 5 limbs */
    uint64_t t[5];
    for (int k = 0; k < 5; k++) {
        uint64_t a = acc->v[4 + k];
        uint64_t b = (5 + k < 10) ? acc->v[5 + k] : 0;
        t[k] = (a >> 2) | (b << 62);
    }
    /* prod = t * MU (5x5 -> 10 limbs), q = prod >> 258 (5 limbs) */
    uint64_t prod[10] = {0};
    for (int i = 0; i < 5; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 5; j++) {
            __uint128_t p = (__uint128_t)t[i] * MU[j];
            uint64_t s1 = prod[i + j] + carry;
            uint64_t c1 = (s1 < carry);
            uint64_t s2 = s1 + (uint64_t)p;
            uint64_t c2 = (s2 < (uint64_t)p);
            prod[i + j] = s2;
            carry = (uint64_t)(p >> 64) + c1 + c2;
        }
        prod[i + 5] += carry;
    }
    uint64_t q[5];
    for (int k = 0; k < 5; k++) {
        uint64_t a = prod[4 + k];
        uint64_t b = (5 + k < 10) ? prod[5 + k] : 0;
        q[k] = (a >> 2) | (b << 62);
    }
    /* qr = q * r (5x4 -> 9 limbs); rem = acc - qr (>= 0 by Barrett) */
    uint64_t qr[10] = {0};
    for (int i = 0; i < 5; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; j++) {
            __uint128_t p = (__uint128_t)q[i] * RMOD[j];
            uint64_t s1 = qr[i + j] + carry;
            uint64_t c1 = (s1 < carry);
            uint64_t s2 = s1 + (uint64_t)p;
            uint64_t c2 = (s2 < (uint64_t)p);
            qr[i + j] = s2;
            carry = (uint64_t)(p >> 64) + c1 + c2;
        }
        qr[i + 4] += carry;
    }
    uint64_t rem[10];
    uint64_t borrow = 0;
    for (int k = 0; k < 10; k++) {
        uint64_t sub = (k < 9) ? qr[k] : 0;
        uint64_t t3 = acc->v[k] - sub - borrow;
        borrow = (acc->v[k] < sub) | ((acc->v[k] == sub) & (borrow != 0));
        rem[k] = t3;
    }
    /* Barrett slack: qhat can undershoot by up to ~10; loop the compare-subtract */
    while (ge_r(rem)) sub_r(rem);
    memcpy(out->v, rem, 32);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: la_big check | la_big bench rows W reps\n"); return 1; }
    int ischeck = !strcmp(argv[1], "check");
    int rows = ischeck ? 64 : atoi(argv[2]);
    int W = ischeck ? 48 : atoi(argv[3]);
    int reps = ischeck ? 1 : atoi(argv[4]);
    fe_t *ents = malloc(sizeof(fe_t) * (size_t)rows * W);
    fe_t *xv = malloc(sizeof(fe_t) * W);
    acc10_t *accs = malloc(sizeof(acc10_t) * rows);
    fe_t *outs = malloc(sizeof(fe_t) * rows);
    if (!ents || !xv || !accs || !outs) return 1;
    for (int i = 0; i < rows * W; i++) fe_rand(&ents[i]);
    for (int j = 0; j < W; j++) fe_rand(&xv[j]);
    for (int r = 0; r < reps; r++) {
        memset(accs, 0, sizeof(acc10_t) * rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < W; j++)
                big_mac(&accs[i], &ents[(size_t)i * W + j], &xv[j]);
            row_finish(&accs[i], &outs[i]);
        }
    }
    if (reps == 0) {   /* setup-only baseline still touches all data deterministically */
        for (int i = 0; i < rows; i++) {
            memset(&accs[i], 0, sizeof(acc10_t));
            for (int j = 0; j < W; j++)
                big_mac(&accs[i], &ents[(size_t)i * W + j], &xv[j]);
            row_finish(&accs[i], &outs[i]);
            memset(&accs[i], 0, sizeof(acc10_t));
        }
    }
    if (ischeck) {
        for (int i = 0; i < rows; i++)
            printf("%d %d %016llx%016llx%016llx%016llx\n", i, W,
                   (unsigned long long)outs[i].v[3], (unsigned long long)outs[i].v[2],
                   (unsigned long long)outs[i].v[1], (unsigned long long)outs[i].v[0]);
        if (getenv("LA_BIG_RAW")) {
            for (int i = 0; i < rows; i++)
                printf("acc %d %016llx%016llx%016llx%016llx%016llx%016llx%016llx%016llx%016llx%016llx\n", i,
                       (unsigned long long)accs[i].v[9], (unsigned long long)accs[i].v[8],
                       (unsigned long long)accs[i].v[7], (unsigned long long)accs[i].v[6],
                       (unsigned long long)accs[i].v[5], (unsigned long long)accs[i].v[4],
                       (unsigned long long)accs[i].v[3], (unsigned long long)accs[i].v[2],
                       (unsigned long long)accs[i].v[1], (unsigned long long)accs[i].v[0]);
        }
    } else {
        uint64_t sink = 0;
        for (int i = 0; i < rows; i++) sink ^= outs[i].v[0];
        printf("bench rows=%d W=%d reps=%d sink=%llu\n", rows, W, reps,
               (unsigned long long)sink);
    }
    return 0;
}
