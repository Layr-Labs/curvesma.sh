/* secp256k1 batched affine addition for the six-class rho walk (freestanding C).
 *
 * Field elements are 4 little-endian 64-bit limbs. p = 2^256 - C0, C0 = 2^32 + 977,
 * so 2^256 = C0 (mod p). "weak" = any value in [0, 2^256); "canonical" = [0, p).
 *   fe_mul, fe_sqr : any weak inputs -> weak output congruent to the product
 *   fe_sub(a, b)   : a weak, b canonical -> weak output congruent to a - b
 *   fe_canon       : weak -> canonical
 * Walk points (x, y) are always canonical.
 */
#include <stddef.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;
#define C0 0x1000003D1ULL
#define AI __attribute__((always_inline)) static inline

typedef struct { u64 v[4]; } fe;
typedef struct { u64 x[4]; u64 y[4]; } pt;

AI u64 mhi(u64 x, u64 y) { return (u64)(((u128)x * y) >> 64); }

AI void fe_mul(u64 r[4], const u64 a[4], const u64 b[4]) {
  u64 z[8], c;
  u64 b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
  {
    u64 a0 = a[0];
    u64 l0 = a0 * b0, h0 = mhi(a0, b0);
    u64 l1 = a0 * b1, h1 = mhi(a0, b1);
    u64 l2 = a0 * b2, h2 = mhi(a0, b2);
    u64 l3 = a0 * b3, h3 = mhi(a0, b3);
    z[0] = l0;
    z[1] = __builtin_addcll(l1, h0, 0, &c);
    z[2] = __builtin_addcll(l2, h1, c, &c);
    z[3] = __builtin_addcll(l3, h2, c, &c);
    z[4] = h3 + c;
  }
  for (int i = 1; i < 4; i++) {
    u64 ai = a[i];
    u64 l0 = ai * b0, h0 = mhi(ai, b0);
    u64 l1 = ai * b1, h1 = mhi(ai, b1);
    u64 l2 = ai * b2, h2 = mhi(ai, b2);
    u64 l3 = ai * b3, h3 = mhi(ai, b3);
    u64 p1 = __builtin_addcll(l1, h0, 0, &c);
    u64 p2 = __builtin_addcll(l2, h1, c, &c);
    u64 p3 = __builtin_addcll(l3, h2, c, &c);
    u64 p4 = h3 + c;
    z[i] = __builtin_addcll(z[i], l0, 0, &c);
    z[i + 1] = __builtin_addcll(z[i + 1], p1, c, &c);
    z[i + 2] = __builtin_addcll(z[i + 2], p2, c, &c);
    z[i + 3] = __builtin_addcll(z[i + 3], p3, c, &c);
    z[i + 4] = p4 + c;
  }
  /* z = L + H*2^256 = L + H*C0 (mod p); s = L + H*C0 < 2^256 * (1 + C0) */
  u64 g0 = z[4] * C0, k0 = mhi(z[4], C0);
  u64 g1 = z[5] * C0, k1 = mhi(z[5], C0);
  u64 g2 = z[6] * C0, k2 = mhi(z[6], C0);
  u64 g3 = z[7] * C0, k3 = mhi(z[7], C0);
  u64 q1 = __builtin_addcll(g1, k0, 0, &c);
  u64 q2 = __builtin_addcll(g2, k1, c, &c);
  u64 q3 = __builtin_addcll(g3, k2, c, &c);
  u64 q4 = k3 + c;
  u64 s0 = __builtin_addcll(z[0], g0, 0, &c);
  u64 s1 = __builtin_addcll(z[1], q1, c, &c);
  u64 s2 = __builtin_addcll(z[2], q2, c, &c);
  u64 s3 = __builtin_addcll(z[3], q3, c, &c);
  u64 s4 = q4 + c;                       /* s4 <= C0 */
  u64 w0 = s4 * C0, w1 = mhi(s4, C0);    /* s4*C0 < 2^65, w1 <= 1 */
  u64 t0 = __builtin_addcll(s0, w0, 0, &c);
  u64 t1 = __builtin_addcll(s1, w1, c, &c);
  u64 t2 = __builtin_addcll(s2, 0, c, &c);
  u64 t3 = __builtin_addcll(s3, 0, c, &c);
  u64 m = c ? C0 : 0;                    /* carry out: t < 2^65, so t + C0 < 2^256 */
  r[0] = __builtin_addcll(t0, m, 0, &c);
  r[1] = __builtin_addcll(t1, 0, c, &c);
  r[2] = __builtin_addcll(t2, 0, c, &c);
  r[3] = t3 + c;
}

AI void fe_sqr(u64 r[4], const u64 a[4]) {
  u64 a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
  u64 c;
  u64 l01 = a0 * a1, h01 = mhi(a0, a1);
  u64 l02 = a0 * a2, h02 = mhi(a0, a2);
  u64 l03 = a0 * a3, h03 = mhi(a0, a3);
  u64 l12 = a1 * a2, h12 = mhi(a1, a2);
  u64 l13 = a1 * a3, h13 = mhi(a1, a3);
  u64 l23 = a2 * a3, h23 = mhi(a2, a3);
  u64 o1 = l01;
  u64 o2 = __builtin_addcll(h01, l02, 0, &c);
  u64 o3 = __builtin_addcll(h02, l03, c, &c);
  u64 o4 = h03 + c;
  u64 p3 = l12;
  u64 p4 = __builtin_addcll(h12, l13, 0, &c);
  u64 p5 = h13 + c;
  o3 = __builtin_addcll(o3, p3, 0, &c);
  o4 = __builtin_addcll(o4, p4, c, &c);
  u64 o5 = __builtin_addcll(p5, l23, c, &c);
  u64 o6 = h23 + c;
  u64 d7 = o6 >> 63;
  u64 d6 = (o6 << 1) | (o5 >> 63);
  u64 d5 = (o5 << 1) | (o4 >> 63);
  u64 d4 = (o4 << 1) | (o3 >> 63);
  u64 d3 = (o3 << 1) | (o2 >> 63);
  u64 d2 = (o2 << 1) | (o1 >> 63);
  u64 d1 = o1 << 1;
  u64 z0 = a0 * a0, g0 = mhi(a0, a0);
  u64 y2 = a1 * a1, g1 = mhi(a1, a1);
  u64 y4 = a2 * a2, g2 = mhi(a2, a2);
  u64 y6 = a3 * a3, g3 = mhi(a3, a3);
  u64 z1 = __builtin_addcll(d1, g0, 0, &c);
  u64 z2 = __builtin_addcll(d2, y2, c, &c);
  u64 z3 = __builtin_addcll(d3, g1, c, &c);
  u64 z4 = __builtin_addcll(d4, y4, c, &c);
  u64 z5 = __builtin_addcll(d5, g2, c, &c);
  u64 z6 = __builtin_addcll(d6, y6, c, &c);
  u64 z7 = __builtin_addcll(d7, g3, c, &c);
  u64 e0 = z4 * C0, k0 = mhi(z4, C0);
  u64 e1 = z5 * C0, k1 = mhi(z5, C0);
  u64 e2 = z6 * C0, k2 = mhi(z6, C0);
  u64 e3 = z7 * C0, k3 = mhi(z7, C0);
  u64 q1 = __builtin_addcll(e1, k0, 0, &c);
  u64 q2 = __builtin_addcll(e2, k1, c, &c);
  u64 q3 = __builtin_addcll(e3, k2, c, &c);
  u64 q4 = k3 + c;
  u64 s0 = __builtin_addcll(z0, e0, 0, &c);
  u64 s1 = __builtin_addcll(z1, q1, c, &c);
  u64 s2 = __builtin_addcll(z2, q2, c, &c);
  u64 s3 = __builtin_addcll(z3, q3, c, &c);
  u64 s4 = q4 + c;
  u64 w0 = s4 * C0, w1 = mhi(s4, C0);
  u64 t0 = __builtin_addcll(s0, w0, 0, &c);
  u64 t1 = __builtin_addcll(s1, w1, c, &c);
  u64 t2 = __builtin_addcll(s2, 0, c, &c);
  u64 t3 = __builtin_addcll(s3, 0, c, &c);
  u64 m = c ? C0 : 0;
  r[0] = __builtin_addcll(t0, m, 0, &c);
  r[1] = __builtin_addcll(t1, 0, c, &c);
  r[2] = __builtin_addcll(t2, 0, c, &c);
  r[3] = t3 + c;
}

/* r = a - b mod p, b canonical. A borrow out means a < b < p, so the wrapped value
 * a - b + 2^256 lies in (C0, 2^256) and subtracting C0 (= adding p mod 2^256) cannot borrow. */
AI void fe_sub(u64 r[4], const u64 a[4], const u64 b[4]) {
  u64 bw;
  u64 d0 = __builtin_subcll(a[0], b[0], 0, &bw);
  u64 d1 = __builtin_subcll(a[1], b[1], bw, &bw);
  u64 d2 = __builtin_subcll(a[2], b[2], bw, &bw);
  u64 d3 = __builtin_subcll(a[3], b[3], bw, &bw);
  u64 m = bw ? C0 : 0;
  r[0] = __builtin_subcll(d0, m, 0, &bw);
  r[1] = __builtin_subcll(d1, 0, bw, &bw);
  r[2] = __builtin_subcll(d2, 0, bw, &bw);
  r[3] = __builtin_subcll(d3, 0, bw, &bw);
}

/* r < 2^256: r + C0 carries out of 2^256 exactly when r >= p, and then r - p = (r + C0) mod 2^256 */
AI void fe_canon(u64 r[4]) {
  u64 c;
  u64 t0 = __builtin_addcll(r[0], C0, 0, &c);
  u64 t1 = __builtin_addcll(r[1], 0, c, &c);
  u64 t2 = __builtin_addcll(r[2], 0, c, &c);
  u64 t3 = __builtin_addcll(r[3], 0, c, &c);
  u64 r0 = c ? t0 : r[0], r1 = c ? t1 : r[1], r2 = c ? t2 : r[2], r3 = c ? t3 : r[3];
  r[0] = r0; r[1] = r1; r[2] = r2; r[3] = r3;
}

AI void ld(u64 d[4], const u64 s[4]) { d[0] = s[0]; d[1] = s[1]; d[2] = s[2]; d[3] = s[3]; }

/* Pass 1 (forward), n >= 2 walks W_i and table points T_i = *tp[i], all canonical,
 * x(T_i) != x(W_i): acc[i] = dx_0 * ... * dx_i with dx_i = x(T_i) - x(W_i). */
void pass1(size_t n, const pt *w, const pt *const *tp, fe *acc) {
  u64 dx[4];
  fe_sub(dx, tp[0]->x, w[0].x);
  ld(acc[0].v, dx);
  for (size_t i = 1; i < n; i++) {
    fe_sub(dx, tp[i]->x, w[i].x);
    fe_mul(acc[i].v, acc[i - 1].v, dx);
  }
}

/* W <- W + T given inv = 1/(x(T) - x(W)) */
AI void finish(pt *w, const pt *t, const u64 inv[4]) {
  u64 dy[4], lam[4], l2[4], x3[4], tt[4], uu[4], y3[4];
  fe_sub(dy, t->y, w->y);
  fe_mul(lam, dy, inv);
  fe_sqr(l2, lam);
  fe_sub(x3, l2, w->x);
  fe_sub(x3, x3, t->x);
  fe_canon(x3);
  fe_sub(tt, w->x, x3);
  fe_mul(uu, lam, tt);
  fe_sub(y3, uu, w->y);
  fe_canon(y3);
  ld(w->x, x3);
  ld(w->y, y3);
}

/* Pass 2 (backward). Rinv = 1/acc[n-1]. Replaces every W_i by W_i + T_i. */
void pass2(size_t n, pt *w, const pt *const *tp, const fe *acc, const u64 Rinv[4]) {
  u64 R[4], inv[4], dx[4];
  ld(R, Rinv);
  for (size_t i = n - 1; i > 0; i--) {
    const pt *t = tp[i];
    fe_mul(inv, R, acc[i - 1].v);
    fe_sub(dx, t->x, w[i].x);
    fe_mul(R, R, dx);
    finish(&w[i], t, inv);
  }
  finish(&w[0], tp[0], R);
}

/* Callable forms of the field operations, used by the inversion and by the tests. */
__attribute__((noinline)) void fe_mul_f(u64 r[4], const u64 a[4], const u64 b[4]) { fe_mul(r, a, b); }
__attribute__((noinline)) void fe_sqr_f(u64 r[4], const u64 a[4]) { fe_sqr(r, a); }

static void sqn(u64 r[4], const u64 a[4], int n) {
  fe_sqr_f(r, a);
  for (int i = 1; i < n; i++) fe_sqr_f(r, r);
}

/* a^(p-2): the libsecp256k1 addition chain, 255 squarings and 15 multiplications */
void fe_inv(u64 r[4], const u64 a[4]) {
  u64 x2[4], x3[4], x6[4], x9[4], x11[4], x22[4], x44[4], x88[4], x176[4], x220[4], x223[4], t[4];
  fe_sqr_f(x2, a); fe_mul_f(x2, x2, a);
  fe_sqr_f(x3, x2); fe_mul_f(x3, x3, a);
  sqn(x6, x3, 3); fe_mul_f(x6, x6, x3);
  sqn(x9, x6, 3); fe_mul_f(x9, x9, x3);
  sqn(x11, x9, 2); fe_mul_f(x11, x11, x2);
  sqn(x22, x11, 11); fe_mul_f(x22, x22, x11);
  sqn(x44, x22, 22); fe_mul_f(x44, x44, x22);
  sqn(x88, x44, 44); fe_mul_f(x88, x88, x44);
  sqn(x176, x88, 88); fe_mul_f(x176, x176, x88);
  sqn(x220, x176, 44); fe_mul_f(x220, x220, x44);
  sqn(x223, x220, 3); fe_mul_f(x223, x223, x3);
  sqn(t, x223, 23); fe_mul_f(t, t, x22);
  sqn(t, t, 5); fe_mul_f(t, t, a);
  sqn(t, t, 3); fe_mul_f(t, t, x2);
  sqn(t, t, 2); fe_mul_f(r, t, a);
}

/* Test driver (not counted): one full round of batched additions. */
void batch_add(size_t n, pt *w, const pt *const *tp, fe *acc) {
  u64 Rinv[4];
  pass1(n, w, tp, acc);
  fe_inv(Rinv, acc[n - 1].v);
  pass2(n, w, tp, acc, Rinv);
}
