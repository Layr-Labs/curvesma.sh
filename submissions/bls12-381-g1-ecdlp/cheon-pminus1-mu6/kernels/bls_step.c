/* Per-step instruction counts on BLS12-381 G1 with pinned blst (Montgomery-form Fp).
   Program "rho" (argv[3]=0): class-rho walk step, same structure as secp_step.c.
     mode 0: batched incomplete affine addition only; mode 2: mode 0 without the per-batch inversion;
     mode 1: full step (beta*x, -x-beta*x, min of three on Montgomery limbs, even-limb y, index,
             look-ahead and DP tests, state store).
   Program "cheon" (argv[3]=1): one BSGS element of the Cheon attack.
     mode 0: two batched additions T0[w0]+T1[w1], then +T2[w2], indices from a counter;
     mode 2: mode 0 without the per-batch inversions;
     mode 1: baby element = fr Montgomery multiply, 85-bit window split, 3 table loads, two batched
             additions, class key, open-addressing insert (marginal load 0.25..0.5);
     mode 3: giant element = same but unsuccessful lookup in a table prefilled to load 0.5.
   Iterations = steps/elements; batch B = 256. Field values are random (constant-time blst ops). */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"
#define B 256
#define R 2048
#define TW 1024
#define HSLOTS 16384
typedef struct { blst_fp x, y; } aff;
static aff X[B], TR[R], T0[TW], T1[TW], T2[TW], S[B];
static blst_fp D[B], ACC[B], BETA;
static uint32_t H[B];
static blst_fr C[B], XI;
typedef struct { uint64_t k[6]; uint64_t j[4]; } slot;
static slot *HT;
static uint64_t seed = 88172645463325252ULL, ctr = 0, hits = 0;
static uint64_t rnd(void) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; return seed; }
static void rfp(blst_fp *a) { uint64_t v[6]; for (int i = 0; i < 6; i++) v[i] = rnd(); v[5] &= 0x0fffffffffffffffULL; blst_fp_from_uint64(a, v); }
static void rfr(blst_fr *a) { uint64_t v[4]; for (int i = 0; i < 4; i++) v[i] = rnd(); v[3] &= 0x0fffffffffffffffULL; blst_fr_from_uint64(a, v); }
static int lt(const blst_fp *a, const blst_fp *b) {
    for (int i = 5; i >= 0; i--) { if (a->l[i] != b->l[i]) return a->l[i] < b->l[i]; }
    return 0;
}

/* out[k] = P[k] + Q[k] for k < B, batched (Montgomery trick); inv toggles the one inversion. */
static void batch_add(aff *out, aff *const *P, aff *const *Q, int noinv) {
    blst_fp acc, inv;
    for (int k = 0; k < B; k++) {
        blst_fp_sub(&D[k], &Q[k]->x, &P[k]->x);
        if (k == 0) acc = D[0]; else blst_fp_mul(&acc, &ACC[k - 1], &D[k]);
        ACC[k] = acc;
    }
    if (noinv) inv = ACC[B - 1]; else blst_fp_inverse(&inv, &ACC[B - 1]);
    for (int k = B - 1; k >= 0; k--) {
        blst_fp ik, lam, t, x3, y3;
        if (k) blst_fp_mul(&ik, &ACC[k - 1], &inv); else ik = inv;
        blst_fp_mul(&inv, &inv, &D[k]);
        blst_fp_sub(&t, &Q[k]->y, &P[k]->y);
        blst_fp_mul(&lam, &t, &ik);
        blst_fp_sqr(&x3, &lam);
        blst_fp_sub(&x3, &x3, &P[k]->x);
        blst_fp_sub(&x3, &x3, &Q[k]->x);
        blst_fp_sub(&t, &P[k]->x, &x3);
        blst_fp_mul(&y3, &lam, &t);
        blst_fp_sub(&y3, &y3, &P[k]->y);
        out[k].x = x3; out[k].y = y3;
    }
}

static void canon_key(blst_fp *best, const blst_fp *x) {
    blst_fp bx, b2x;
    blst_fp_mul(&bx, x, &BETA);
    blst_fp_add(&b2x, x, &bx);
    blst_fp_cneg(&b2x, &b2x, 1);
    *best = *x;
    if (lt(&bx, best)) *best = bx;
    if (lt(&b2x, best)) *best = b2x;
}

__attribute__((noinline)) static void rho_round(int mode, uint64_t *dps) {
    static aff *P[B], *Q[B];
    static aff OUT[B];
    for (int k = 0; k < B; k++) { P[k] = &X[k]; Q[k] = &TR[H[k]]; }
    batch_add(OUT, P, Q, mode == 2);
    for (int k = 0; k < B; k++) {
        if (mode != 1) { X[k] = OUT[k]; H[k] = (uint32_t)(OUT[k].x.l[0] & (R - 1)); continue; }
        blst_fp best;
        canon_key(&best, &OUT[k].x);
        blst_fp_cneg(&OUT[k].y, &OUT[k].y, (int)(OUT[k].y.l[0] & 1));
        uint32_t nh = (uint32_t)(best.l[0] & (R - 1));
        if (nh == H[k]) nh = (nh + 1) & (R - 1);
        if (((best.l[0] >> 11) & 0xffffffffffULL) == 0) (*dps)++;
        X[k].x = best; X[k].y = OUT[k].y; H[k] = nh;
    }
}

static void insert(const blst_fp *key, uint64_t j) {
    uint64_t i = key->l[0] & (HSLOTS - 1);
    while (HT[i].k[0] | HT[i].k[1]) i = (i + 1) & (HSLOTS - 1);
    memcpy(HT[i].k, key->l, 48);
    HT[i].j[0] = j; HT[i].j[1] = HT[i].j[2] = HT[i].j[3] = 0;
}
static int lookup(const blst_fp *key) {
    uint64_t i = key->l[0] & (HSLOTS - 1);
    while (HT[i].k[0] | HT[i].k[1]) {
        if (HT[i].k[0] == key->l[0] && memcmp(HT[i].k, key->l, 48) == 0) return 1;
        i = (i + 1) & (HSLOTS - 1);
    }
    return 0;
}

__attribute__((noinline)) static void cheon_round(int mode) {
    static aff *P[B], *Q[B];
    static uint32_t w0[B], w1[B], w2[B];
    for (int k = 0; k < B; k++) {
        if (mode == 1 || mode == 3) {
            blst_fr_mul(&C[k], &C[k], &XI);
            const uint64_t *l = C[k].l;
            uint64_t a0 = l[0], a1 = l[1] & ((1ULL << 21) - 1);
            uint64_t b0 = (l[1] >> 21) | (l[2] << 43), b1 = (l[2] >> 21) & ((1ULL << 21) - 1);
            uint64_t c0 = (l[2] >> 42) | (l[3] << 22), c1 = l[3] >> 42;
            /* real tables are indexed by the full 85-bit windows; here reduced mod TW for memory */
            w0[k] = (uint32_t)((a0 ^ a1) & (TW - 1));
            w1[k] = (uint32_t)((b0 ^ b1) & (TW - 1));
            w2[k] = (uint32_t)((c0 ^ c1) & (TW - 1));
        } else {
            w0[k] = (uint32_t)(ctr++ & (TW - 1)); w1[k] = (uint32_t)((ctr * 7) & (TW - 1)); w2[k] = (uint32_t)((ctr * 13) & (TW - 1));
        }
        P[k] = &T0[w0[k]]; Q[k] = &T1[w1[k]];
    }
    batch_add(S, P, Q, mode == 2);
    for (int k = 0; k < B; k++) { P[k] = &S[k]; Q[k] = &T2[w2[k]]; }
    static aff OUT[B];
    batch_add(OUT, P, Q, mode == 2);
    if (mode == 0 || mode == 2) { X[0] = OUT[B - 1]; return; }
    for (int k = 0; k < B; k++) {
        blst_fp best;
        canon_key(&best, &OUT[k].x);
        if (mode == 1) insert(&best, ctr++); else hits += (uint64_t)lookup(&best);
    }
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 4096, mode = argc > 2 ? atoi(argv[2]) : 1;
    int prog = argc > 3 ? atoi(argv[3]) : 0;
    uint64_t dps = 0;
    HT = calloc(HSLOTS, sizeof(slot));
    rfp(&BETA); rfr(&XI);
    for (int i = 0; i < R; i++) { rfp(&TR[i].x); rfp(&TR[i].y); }
    for (int i = 0; i < TW; i++) { rfp(&T0[i].x); rfp(&T0[i].y); rfp(&T1[i].x); rfp(&T1[i].y); rfp(&T2[i].x); rfp(&T2[i].y); }
    for (int k = 0; k < B; k++) { rfp(&X[k].x); rfp(&X[k].y); H[k] = (uint32_t)(rnd() & (R - 1)); rfr(&C[k]); }
    if (prog == 1 && mode == 3) for (int i = 0; i < HSLOTS / 2; i++) { blst_fp k; rfp(&k); insert(&k, i); }
    if (prog == 0) rho_round(mode, &dps); else cheon_round(mode == 1 ? 0 : mode);  /* warm-up */
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for (int i = 0; i < n / B; i++) { if (prog == 0) rho_round(mode, &dps); else cheon_round(mode); }
    CALLGRIND_TOGGLE_COLLECT;
    printf("%lu %lu %lu\n", (unsigned long)X[0].x.l[0], (unsigned long)dps, (unsigned long)hits);
    return 0;
}
