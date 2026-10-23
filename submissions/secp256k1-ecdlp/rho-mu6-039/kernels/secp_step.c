/* Per-step instruction count of the class-rho walk on secp256k1 (libsecp256k1 field code).
   mode 0: batched incomplete affine addition X + R_h only (Montgomery trick, no normalisation).
   mode 1: full step = mode 0 + class canonicalisation (beta*x, -x-beta*x, normalise, min of three,
           even y), index h, look-ahead test, DP test, state store.
   mode 2: mode 0 without the single per-batch inversion (isolates its amortised cost).
   Overhead per step = slope(mode1) - slope(mode0); inversion share = slope(mode0) - slope(mode2). Iterations = walker steps; batch B = 256. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <valgrind/callgrind.h>
#define SECP256K1_WIDEMUL_INT128 1
#define SECP256K1_INT128_NATIVE 1
#include "/secp/src/int128_impl.h"
#include "/secp/src/field_impl.h"
#include "/secp/src/modinv64_impl.h"
#define B 256
#define R 2048
static secp256k1_fe X[B], Y[B], TX[R], TY[R], D[B], ACC[B];
static uint32_t H[B];
static secp256k1_fe BETA;
static uint64_t seed = 88172645463325252ULL;
static uint64_t rnd(void) { seed ^= seed << 13; seed ^= seed >> 7; seed ^= seed << 17; return seed; }
static void rfe(secp256k1_fe *a) {
    unsigned char b[32];
    for (int i = 0; i < 32; i++) b[i] = (unsigned char)rnd();
    b[0] &= 0x7f;
    secp256k1_fe_set_b32_mod(a, b);
    secp256k1_fe_normalize(a);
}
__attribute__((noinline)) static void round_(int mode, uint64_t *dps) {
    secp256k1_fe acc, inv, t, nx;
    for (int k = 0; k < B; k++) {
        uint32_t h = H[k];
        secp256k1_fe_negate(&nx, &X[k], 1);
        D[k] = TX[h];
        secp256k1_fe_add(&D[k], &nx);
        if (k == 0) acc = D[0]; else secp256k1_fe_mul(&acc, &ACC[k - 1], &D[k]);
        ACC[k] = acc;
    }
    if (mode == 2) inv = ACC[B - 1]; else secp256k1_fe_inv_var(&inv, &ACC[B - 1]);
    for (int k = B - 1; k >= 0; k--) {
        uint32_t h = H[k];
        secp256k1_fe ik, lam, l2, x3, y3, ny;
        if (k) secp256k1_fe_mul(&ik, &ACC[k - 1], &inv); else ik = inv;
        secp256k1_fe_mul(&inv, &inv, &D[k]);
        secp256k1_fe_negate(&ny, &Y[k], 1);
        t = TY[h];
        secp256k1_fe_add(&t, &ny);
        secp256k1_fe_mul(&lam, &t, &ik);
        secp256k1_fe_sqr(&l2, &lam);
        x3 = l2;
        secp256k1_fe_negate(&nx, &X[k], 1);
        secp256k1_fe_add(&x3, &nx);
        secp256k1_fe_negate(&nx, &TX[h], 1);
        secp256k1_fe_add(&x3, &nx);
        secp256k1_fe_negate(&nx, &x3, 5);
        t = X[k];
        secp256k1_fe_add(&t, &nx);
        secp256k1_fe_mul(&y3, &lam, &t);
        secp256k1_fe_add(&y3, &ny);
        if (mode != 1) {
            X[k] = x3; Y[k] = y3;
            H[k] = (uint32_t)(x3.n[0] & (R - 1));
            continue;
        }
        secp256k1_fe bx, b2x, best;
        secp256k1_fe_normalize_var(&x3);
        secp256k1_fe_mul(&bx, &x3, &BETA);
        secp256k1_fe_normalize_var(&bx);
        secp256k1_fe_negate(&b2x, &x3, 1);
        secp256k1_fe_negate(&nx, &bx, 1);
        secp256k1_fe_add(&b2x, &nx);
        secp256k1_fe_normalize_var(&b2x);
        best = x3;
        if (secp256k1_fe_cmp_var(&bx, &best) < 0) best = bx;
        if (secp256k1_fe_cmp_var(&b2x, &best) < 0) best = b2x;
        secp256k1_fe_normalize_var(&y3);
        if (secp256k1_fe_is_odd(&y3)) { secp256k1_fe_negate(&y3, &y3, 1); secp256k1_fe_normalize_var(&y3); }
        uint32_t nh = (uint32_t)(best.n[0] & (R - 1));
        if (nh == h) nh = (nh + 1) & (R - 1);          /* look-ahead branch (taken ~1/R) */
        if (((best.n[0] >> 11) & 0xffffffffffULL) == 0) (*dps)++;  /* DP test */
        X[k] = best; Y[k] = y3; H[k] = nh;
    }
}
int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 4096, mode = argc > 2 ? atoi(argv[2]) : 1;
    uint64_t dps = 0;
    static const unsigned char beta_b[32] = {
        0x7a,0xe9,0x6a,0x2b,0x65,0x7c,0x07,0x10,0x6e,0x64,0x47,0x9e,0xac,0x34,0x34,0xe9,
        0x9c,0xf0,0x49,0x75,0x12,0xf5,0x89,0x95,0xc1,0x39,0x6c,0x28,0x71,0x95,0x01,0xee};
    secp256k1_fe_set_b32_mod(&BETA, beta_b);
    for (int i = 0; i < R; i++) { rfe(&TX[i]); rfe(&TY[i]); }
    for (int k = 0; k < B; k++) { rfe(&X[k]); rfe(&Y[k]); H[k] = (uint32_t)(rnd() & (R - 1)); }
    round_(mode, &dps);
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for (int i = 0; i < n / B; i++) round_(mode, &dps);
    CALLGRIND_TOGGLE_COLLECT;
    printf("%lu %lu\n", (unsigned long)X[0].n[0], (unsigned long)dps);
    return 0;
}
