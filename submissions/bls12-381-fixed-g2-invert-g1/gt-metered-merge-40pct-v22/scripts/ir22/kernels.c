/* kernels.c -- counted AArch64 blocks of the v22 ledger (symbol boundaries disjoint).
 *
 *  key_six_plusminus : twelve class keys from one Fp12 element: two 6-limb
 *                      Montgomery products by zeta + ten modular add/sub.
 *  min_rotation_booth: argmin over the twelve keys (6-limb compare, no early
 *                      exit).  Equal to the minimal cyclic rotation when the
 *                      minimum is unique; ties take the counted full-rotation
 *                      fallback `min_rotation_full` at the measured tie rate.
 *  min_rotation_full : full quadratic rotation compare (fallback symbol).
 *  cyclotomic_sqr    : y -> y^2 in the shipped tower (Fp2=u^2+1, Fp6=v^3-(1+u),
 *                      Fp12=w^2-v), formula (1+2 v c1^2, 2 c0 c1), Karatsuba
 *                      Fp6, specialised Fp2 squaring.
 * P6/ZETA6 are extern (consts.c) so no priced table can be constant-folded.
 */
#include <stdint.h>
extern const uint64_t P6[6], NP0, ZETA6[6];
typedef uint64_t u64;
typedef unsigned __int128 u128;
#define AI __attribute__((always_inline)) static inline

AI void fpcopy(u64 *r, const u64 *a) {
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) r[i] = a[i];
}

/* 6-limb Montgomery product r = a*b*2^-384 mod p (CIOS) */
AI void fpmul(u64 r[6], const u64 a[6], const u64 b[6]) {
    u64 t[7] = {0,0,0,0,0,0,0};
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) {
        u64 carry = 0;
#pragma clang loop unroll(full)
        for (int j = 0; j < 6; j++) {
            u128 s = (u128)a[j] * b[i] + t[j] + carry;
            t[j] = (u64)s; carry = (u64)(s >> 64);
        }
        u128 s = (u128)t[6] + carry; t[6] = (u64)s;
        u64 m = t[0] * NP0;
        u128 s2 = (u128)m * P6[0] + t[0];
        carry = (u64)(s2 >> 64);
#pragma clang loop unroll(full)
        for (int j = 1; j < 6; j++) {
            u128 s3 = (u128)m * P6[j] + t[j] + carry;
            t[j - 1] = (u64)s3; carry = (u64)(s3 >> 64);
        }
        u128 s4 = (u128)t[6] + carry; t[5] = (u64)s4; t[6] = (u64)(s4 >> 64);
    }
    u64 d[6]; u64 borrow = 0;
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) { u128 s = (u128)t[i] - P6[i] - borrow; d[i] = (u64)s; borrow = (u64)((s >> 64) & 1); }
    u64 mask = borrow - 1;
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) r[i] = (t[i] & ~mask) | (d[i] & mask);
}

AI void fpadd(u64 r[6], const u64 a[6], const u64 b[6]) {
    u64 c = 0; u64 t[6];
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) { u128 s = (u128)a[i] + b[i] + c; t[i] = (u64)s; c = (u64)(s >> 64); }
    u64 d[6], br = 0;
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) { u128 s = (u128)t[i] - P6[i] - br; d[i] = (u64)s; br = (u64)((s >> 64) & 1); }
    u64 m = br - 1;
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) r[i] = (t[i] & ~m) | (d[i] & m);
}

AI void fpsub(u64 r[6], const u64 a[6], const u64 b[6]) {
    u64 br = 0; u64 t[6];
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) { u128 s = (u128)a[i] - b[i] - br; t[i] = (u64)s; br = (u64)((s >> 64) & 1); }
    u64 m = 0 - br; u64 c = 0;
#pragma clang loop unroll(full)
    for (int i = 0; i < 6; i++) { u128 s = (u128)t[i] + (P6[i] & m) + c; r[i] = (u64)s; c = (u64)(s >> 64); }
}

AI void fpneg(u64 r[6], const u64 a[6]) {
    const u64 z[6] = {0,0,0,0,0,0};
    fpsub(r, z, a);
}

void key_six_plusminus(u64 out[12][6], const u64 A[6], const u64 B[6]) {
    u64 Az[6], Bz[6], t[6], nA[6], nB[6], nAz[6], nBz[6];
    fpmul(Az, A, ZETA6);
    fpmul(Bz, B, ZETA6);
    fpneg(nA, A); fpneg(nB, B); fpneg(nAz, Az); fpneg(nBz, Bz);
    fpsub(t, Az, A); fpcopy(out[4], t);
    fpsub(t, Bz, B); fpcopy(out[5], t);
    fpsub(t, A, Az); fpcopy(out[10], t);
    fpsub(t, B, Bz); fpcopy(out[11], t);
    fpcopy(out[0], A); fpcopy(out[1], B); fpcopy(out[2], Az); fpcopy(out[3], Bz);
    fpcopy(out[6], nA); fpcopy(out[7], nB); fpcopy(out[8], nAz); fpcopy(out[9], nBz);
}

AI int keycmp(const u64 a[6], const u64 b[6]) {
#pragma clang loop unroll(full)
    for (int i = 5; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

int min_rotation_booth(const u64 kv[12][6]) {
    int best = 0;
#pragma clang loop unroll(full)
    for (int i = 1; i < 12; i++)
        if (keycmp(kv[i], kv[best]) < 0) best = i;
    return best;
}

int min_rotation_full(const u64 kv[12][6]) {
    int best = 0;
#pragma clang loop unroll(full)
    for (int i = 1; i < 12; i++) {
        int c = 0;
#pragma clang loop unroll(full)
        for (int j = 0; j < 12 && c == 0; j++)
            c = keycmp(kv[(i + j) % 12], kv[(best + j) % 12]);
        if (c < 0) best = i;
    }
    return best;
}

typedef struct { u64 c[2][6]; } fp2;
typedef struct { fp2 c[3]; } fp6;
typedef struct { fp6 c[2]; } fp12;

AI void f2mul(fp2 *r, const fp2 *a, const fp2 *b) {
    u64 t0[6], t1[6], s0[6], s1[6], rr[6], ii[6];
    fpmul(t0, a->c[0], b->c[0]); fpmul(t1, a->c[1], b->c[1]);
    fpadd(s0, a->c[0], a->c[1]); fpadd(s1, b->c[0], b->c[1]);
    fpmul(rr, s0, s1); fpsub(rr, rr, t0); fpsub(rr, rr, t1);
    fpadd(ii, t0, t1);
    fpcopy(r->c[0], rr); fpcopy(r->c[1], ii);
}
AI void f2sqr(fp2 *r, const fp2 *a) {          /* (a0+a1)(a0-a1), 2 a0 a1 */
    u64 s[6], d[6], t0[6], t1[6];
    fpadd(s, a->c[0], a->c[1]); fpsub(d, a->c[0], a->c[1]);
    fpmul(t0, s, d); fpmul(t1, a->c[0], a->c[1]);
    fpcopy(r->c[0], t0); fpadd(r->c[1], t1, t1);
}
AI void f2add(fp2 *r, const fp2 *a, const fp2 *b) {
    fpadd(r->c[0], a->c[0], b->c[0]); fpadd(r->c[1], a->c[1], b->c[1]);
}
AI void f2sub(fp2 *r, const fp2 *a, const fp2 *b) {
    fpsub(r->c[0], a->c[0], b->c[0]); fpsub(r->c[1], a->c[1], b->c[1]);
}
AI void f2mulxi(fp2 *r, const fp2 *a) {        /* xi = 1+u : (a0-a1, a0+a1) */
    u64 t0[6], t1[6];
    fpsub(t0, a->c[0], a->c[1]); fpadd(t1, a->c[0], a->c[1]);
    fpcopy(r->c[0], t0); fpcopy(r->c[1], t1);
}
/* Karatsuba Fp6, 6 Fp2 mults: v^3 = xi */
AI void f6mul(fp6 *r, const fp6 *a, const fp6 *b) {
    fp2 t0, t1, t2, u0, u1, u2, k;
    f2mul(&t0, &a->c[0], &b->c[0]);
    f2mul(&t1, &a->c[1], &b->c[1]);
    f2mul(&t2, &a->c[2], &b->c[2]);
    f2add(&u0, &a->c[1], &a->c[2]); f2add(&k, &b->c[1], &b->c[2]); f2mul(&u0, &u0, &k);
    f2add(&u1, &a->c[0], &a->c[2]); f2add(&k, &b->c[0], &b->c[2]); f2mul(&u1, &u1, &k);
    f2add(&u2, &a->c[0], &a->c[1]); f2add(&k, &b->c[0], &b->c[1]); f2mul(&u2, &u2, &k);
    f2sub(&u0, &u0, &t1); f2sub(&u0, &u0, &t2);
    f2sub(&u1, &u1, &t0); f2sub(&u1, &u1, &t2);
    f2sub(&u2, &u2, &t0); f2sub(&u2, &u2, &t1);
    f2mulxi(&k, &u0); f2add(&r->c[0], &t0, &k);
    f2mulxi(&k, &t2); f2add(&r->c[1], &u1, &k);
    f2add(&r->c[2], &u2, &t1);
}
AI void f6sqr(fp6 *r, const fp6 *a) { f6mul(r, a, a); }
AI void f6add(fp6 *r, const fp6 *a, const fp6 *b) {
#pragma clang loop unroll(full)
    for (int i = 0; i < 3; i++) f2add(&r->c[i], &a->c[i], &b->c[i]);
}
AI void f6sub(fp6 *r, const fp6 *a, const fp6 *b) {
#pragma clang loop unroll(full)
    for (int i = 0; i < 3; i++) f2sub(&r->c[i], &a->c[i], &b->c[i]);
}
AI void f6mulv(fp6 *r, const fp6 *a) {          /* (a0,a1,a2)*v = (xi a2, a0, a1) */
    fp2 t0;
    f2mulxi(&t0, &a->c[2]);
    r->c[0] = t0; r->c[1] = a->c[0]; r->c[2] = a->c[1];
}
void cyclotomic_sqr(fp12 *r, const fp12 *a) {
    fp6 t, vt, s, one, x, two_vt;
#pragma clang loop unroll(full)
    for (int i = 0; i < 3; i++) {
        one.c[i].c[0][0] = 1;
#pragma clang loop unroll(full)
        for (int k = 1; k < 6; k++) one.c[i].c[0][k] = 0;
#pragma clang loop unroll(full)
        for (int k = 0; k < 6; k++) one.c[i].c[1][k] = 0;
    }
    f6sqr(&t, &a->c[1]);
    f6mulv(&vt, &t);
    f6add(&x, &a->c[0], &a->c[1]);
    f6sqr(&s, &x);
    f6add(&two_vt, &vt, &vt);
    f6add(&r->c[0], &two_vt, &one);
    f6sub(&x, &s, &vt);
    f6sub(&x, &x, &one);
    f6sub(&r->c[1], &x, &t);
}
