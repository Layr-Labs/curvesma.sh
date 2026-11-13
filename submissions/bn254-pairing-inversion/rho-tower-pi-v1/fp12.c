#include "fp12.h"
#include <string.h>

/* BN254 base field p = 36u^4+36u^3+24u^2+6u+1, u = 4965661367192848881 */
static const uint64_t P[4] = {
  0x3c208c16d87cfd47ULL, 0x97816a916871ca8dULL,
  0xb85045b68181585dULL, 0x30644e72e131a029ULL };
static const uint64_t R2[4] = {
  0xf32cfc5b538afa89ULL, 0xb5e71911d44501fbULL,
  0x47ab1eff0a417ff6ULL, 0x06d89f71cab8351fULL };
static const uint64_t N0 = 0x87d20782e4866389ULL;
const fp FP_P = {{0x3c208c16d87cfd47ULL,0x97816a916871ca8dULL,0xb85045b68181585dULL,0x30644e72e131a029ULL}};
const fp FP_ZERO = {{0,0,0,0}};
const fp FP_ONE_MONT = {{0xd35d438dc58f0d9dULL,0x0a78eb28f5c70b3dULL,0x666ea36f7879462cULL,0x0e0a77c19a07df2fULL}};

static inline int lt4(const uint64_t a[4], const uint64_t b[4]) {
  for (int i = 3; i >= 0; i--) { if (a[i] != b[i]) return a[i] < b[i]; }
  return 0;
}
static inline void sub4(uint64_t z[4], const uint64_t a[4], const uint64_t b[4]) {
  uint64_t borrow = 0;
  for (int i = 0; i < 4; i++) {
    __uint128_t d = (__uint128_t)a[i] - b[i] - borrow;
    z[i] = (uint64_t)d; borrow = (uint64_t)(d >> 64) & 1ULL;
  }
}
static inline void add4(uint64_t z[4], const uint64_t a[4], const uint64_t b[4]) {
  uint64_t carry = 0;
  for (int i = 0; i < 4; i++) {
    __uint128_t s = (__uint128_t)a[i] + b[i] + carry;
    z[i] = (uint64_t)s; carry = (uint64_t)(s >> 64);
  }
}

/* Montgomery multiplication, CIOS with a 9-word accumulator. */
void fp_mul(fp *z, const fp *x, const fp *y) {
  uint64_t t[6] = {0,0,0,0,0,0};
  for (int i = 0; i < 4; i++) {
    uint64_t C = 0;
    for (int j = 0; j < 4; j++) {
      __uint128_t uv = (__uint128_t)x->v[j] * y->v[i] + t[j] + C;
      t[j] = (uint64_t)uv; C = (uint64_t)(uv >> 64);
    }
    { __uint128_t uv = (__uint128_t)t[4] + C; t[4] = (uint64_t)uv; t[5] = (uint64_t)(uv >> 64); }
    uint64_t m = t[0] * N0;
    { __uint128_t uv = (__uint128_t)m * P[0] + t[0]; C = (uint64_t)(uv >> 64); }
    for (int j = 1; j < 4; j++) {
      __uint128_t uv = (__uint128_t)m * P[j] + t[j] + C;
      t[j-1] = (uint64_t)uv; C = (uint64_t)(uv >> 64);
    }
    { __uint128_t uv = (__uint128_t)t[4] + C; t[3] = (uint64_t)uv; t[4] = t[5] + (uint64_t)(uv >> 64); }
  }
  uint64_t d[4]; sub4(d, t, P);
  if (t[4] || !lt4(t, P)) memcpy(z->v, d, sizeof(d));
  else memcpy(z->v, t, sizeof(uint64_t)*4);
}

void fp_add(fp *z, const fp *x, const fp *y) {
  uint64_t s[4]; add4(s, x->v, y->v);
  uint64_t d[4]; sub4(d, s, P);
  if (s[3] >> 63 || !lt4(s, P)) memcpy(z->v, d, sizeof(d));
  else memcpy(z->v, s, sizeof(uint64_t)*4);
}
void fp_sub(fp *z, const fp *x, const fp *y) {
  uint64_t d[4]; sub4(d, x->v, y->v);
  if (lt4(x->v, y->v)) add4(d, d, P);
  memcpy(z->v, d, sizeof(d));
}
int fp_eq(const fp *x, const fp *y) {
  uint64_t d = 0; for (int i=0;i<4;i++) d |= x->v[i]^y->v[i]; return d==0;
}

/* Fp2 = Fp[u]/(u^2+1) */
void fp2_mul(fp2 *z, const fp2 *x, const fp2 *y) {
  fp t0,t1,t2,a,b;
  fp_mul(&t0,&x->a,&y->a);
  fp_mul(&t1,&x->b,&y->b);
  fp_add(&a,&x->a,&x->b); fp_add(&b,&y->a,&y->b);
  fp_mul(&t2,&a,&b);
  fp_sub(&z->a,&t0,&t1);
  fp_sub(&t2,&t2,&t0); fp_sub(&z->b,&t2,&t1);
}
void fp2_add(fp2 *z, const fp2 *x, const fp2 *y){ fp_add(&z->a,&x->a,&y->a); fp_add(&z->b,&x->b,&y->b); }
void fp2_sub(fp2 *z, const fp2 *x, const fp2 *y){ fp_sub(&z->a,&x->a,&y->a); fp_sub(&z->b,&x->b,&y->b); }
static void fp2_mul_xi(fp2 *z, const fp2 *x) { /* (9+u)(a+bu) = (9a-b) + (a+9b)u */
  fp a8,a9,b8,b9,t;
  fp_add(&a8,&x->a,&x->a); fp_add(&a9,&a8,&x->a); fp_add(&a9,&a9,&a9); /* 2a, 4a, 8a? */
  /* recompute: 9a = 8a + a */
  fp_add(&a8,&x->a,&x->a); fp_add(&t,&a8,&a8); fp_add(&a8,&t,&t); /* 8a */
  fp_add(&a9,&a8,&x->a);
  fp_add(&b8,&x->b,&x->b); fp_add(&t,&b8,&b8); fp_add(&b8,&t,&t); /* 8b */
  fp_add(&b9,&b8,&x->b);
  fp_sub(&z->a,&a9,&x->b);
  fp_add(&z->b,&x->a,&b9);
}

/* Fp6 = Fp2[v]/(v^3-xi), xi = 9+u */
void fp6_mul(fp6 *z, const fp6 *x, const fp6 *y) {
  fp2 t0,t1,t2,s0,s1,s2,xi0,xi1;
  fp2_mul(&t0,&x->c[0],&y->c[0]);
  fp2_mul(&t1,&x->c[1],&y->c[1]);
  fp2_mul(&t2,&x->c[2],&y->c[2]);
  /* c0 = t0 + xi*(x1*y2 + x2*y1) */
  fp2_mul(&s0,&x->c[1],&y->c[2]);
  fp2_mul(&s1,&x->c[2],&y->c[1]);
  fp2_add(&s0,&s0,&s1); fp2_mul_xi(&xi0,&s0); fp2_add(&s0,&t0,&xi0);
  /* c1 = x0*y1+x1*y0 + xi*t2 */
  fp2_mul(&s1,&x->c[0],&y->c[1]);
  fp2_mul(&s2,&x->c[1],&y->c[0]);
  fp2_add(&s1,&s1,&s2); fp2_mul_xi(&xi1,&t2); fp2_add(&s1,&s1,&xi1);
  /* c2 = x0*y2+x1*y1+x2*y0 */
  fp2_mul(&s2,&x->c[0],&y->c[2]);
  fp2_add(&s2,&s2,&t1);
  fp2_mul(&xi0,&x->c[2],&y->c[0]); fp2_add(&s2,&s2,&xi0);
  z->c[0]=s0; z->c[1]=s1; z->c[2]=s2;
}
void fp6_add(fp6 *z, const fp6 *x, const fp6 *y){ for(int i=0;i<3;i++) fp2_add(&z->c[i],&x->c[i],&y->c[i]); }
void fp6_sub(fp6 *z, const fp6 *x, const fp6 *y){ for(int i=0;i<3;i++) fp2_sub(&z->c[i],&x->c[i],&y->c[i]); }
static void fp6_mul_v(fp6 *z, const fp6 *x) { /* v*(c0+c1 v+c2 v^2) = xi*c2 + c0 v + c1 v^2 */
  fp2 t; fp2_mul_xi(&t,&x->c[2]);
  z->c[2]=x->c[1]; z->c[1]=x->c[0]; z->c[0]=t;
}

/* Fp12 = Fp6[w]/(w^2-v) */
void fp12_mul(fp12 *z, const fp12 *x, const fp12 *y) {
  fp6 t0,t1,t2,a,b,v1;
  fp6_mul(&t0,&x->c[0],&y->c[0]);
  fp6_mul(&t1,&x->c[1],&y->c[1]);
  fp6_add(&a,&x->c[0],&x->c[1]); fp6_add(&b,&y->c[0],&y->c[1]);
  fp6_mul(&t2,&a,&b);
  fp6_mul_v(&v1,&t1);
  fp6_add(&z->c[0],&t0,&v1);
  fp6_sub(&t2,&t2,&t0); fp6_sub(&z->c[1],&t2,&t1);
}
void fp12_add(fp12 *z, const fp12 *x, const fp12 *y){ fp6_add(&z->c[0],&x->c[0],&y->c[0]); fp6_add(&z->c[1],&x->c[1],&y->c[1]); }
void fp12_sub(fp12 *z, const fp12 *x, const fp12 *y){ fp6_sub(&z->c[0],&x->c[0],&y->c[0]); fp6_sub(&z->c[1],&x->c[1],&y->c[1]); }
int fp12_eq(const fp12 *x, const fp12 *y) {
  int e=1; for(int i=0;i<2;i++) for(int j=0;j<3;j++){ e &= fp_eq(&x->c[i].c[j].a,&y->c[i].c[j].a); e &= fp_eq(&x->c[i].c[j].b,&y->c[i].c[j].b);} return e;
}
void fp12_set_one(fp12 *z) {
  memset(z,0,sizeof(*z)); z->c[0].c[0].a=FP_ONE_MONT;
}
void fp12_pow(fp12 *z, const fp12 *x, const uint64_t e[4]) {
  fp12 r; fp12_set_one(&r);
  for (int i=253;i>=0;i--) {
    fp12_mul(&r,&r,&r);
    if ((e[i>>6]>>(i&63))&1) fp12_mul(&r,&r,x);
  }
  *z=r;
}
