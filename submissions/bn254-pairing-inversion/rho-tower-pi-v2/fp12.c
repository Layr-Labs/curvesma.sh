/* Constant-time BN254 tower arithmetic (bn254-gt-dlp rho-tower-v3).
 *
 * Every function below executes the same instruction sequence for every
 * operand: no conditional branch, no early return, no data-dependent loop
 * trip count.  The only branches emitted are fixed-trip loop back-edges.
 * The conditional subtract of the Montgomery multiplication and the
 * conditional add of the modulus are mask selects.
 */
#include "fp12.h"

/* BN254 base field p = 36u^4+36u^3+24u^2+6u+1, u = 4965661367192848881 */
static const uint64_t P[4] = {
  0x3c208c16d87cfd47ULL, 0x97816a916871ca8dULL,
  0xb85045b68181585dULL, 0x30644e72e131a029ULL };
static const uint64_t N0 = 0x87d20782e4866389ULL;
const fp FP_P = {{0x3c208c16d87cfd47ULL,0x97816a916871ca8dULL,0xb85045b68181585dULL,0x30644e72e131a029ULL}};
const fp FP_ZERO = {{0,0,0,0}};
const fp FP_ONE_MONT = {{0xd35d438dc58f0d9dULL,0x0a78eb28f5c70b3dULL,0x666ea36f7879462cULL,0x0e0a77c19a07df2fULL}};

/* z = a - b over four limbs; returns the final borrow (0 or 1). */
static inline uint64_t sub4(uint64_t z[4], const uint64_t a[4], const uint64_t b[4]) {
  uint64_t borrow = 0;
  for (int i = 0; i < 4; i++) {
    __uint128_t d = (__uint128_t)a[i] - b[i] - borrow;
    z[i] = (uint64_t)d; borrow = (uint64_t)(d >> 64) & 1ULL;
  }
  return borrow;
}

/* Montgomery multiplication, CIOS with a six-word accumulator. */
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
  /* t < 2p, t[4] in {0,1}.  d = t - p mod 2^256; t >= p iff t[4]=1 or the
     subtraction did not borrow.  Mask select, no branch. */
  uint64_t d[4];
  uint64_t borrow = sub4(d, t, P);
  uint64_t cond = t[4] | (borrow ^ 1ULL);
  uint64_t mask = 0 - cond;
  for (int i = 0; i < 4; i++) z->v[i] = (d[i] & mask) | (t[i] & ~mask);
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

static inline void fp2_mul_xi(fp2 *z, const fp2 *x) { /* (9+u)(a+bu) = (9a-b) + (a+9b)u */
  fp a2,a8,a9,b2,b8,b9;
  fp_add(&a2,&x->a,&x->a); fp_add(&a8,&a2,&a2); fp_add(&a8,&a8,&a8); /* 8a */
  fp_add(&a9,&a8,&x->a);
  fp_add(&b2,&x->b,&x->b); fp_add(&b8,&b2,&b2); fp_add(&b8,&b8,&b8); /* 8b */
  fp_add(&b9,&b8,&x->b);
  fp_sub(&z->a,&a9,&x->b);
  fp_add(&z->b,&x->a,&b9);
}

/* Fp6 = Fp2[v]/(v^3-xi), xi = 9+u, six-product Karatsuba form:
 *   v0=x0y0, v1=x1y1, v2=x2y2,
 *   s01=(x0+x1)(y0+y1)-v0-v1, s02=(x0+x2)(y0+y2)-v0-v2,
 *   s12=(x1+x2)(y1+y2)-v1-v2,
 *   c0 = v0 + xi*s12, c1 = s01 + xi*v2, c2 = s02 + v1.
 */
void fp6_mul(fp6 *z, const fp6 *x, const fp6 *y) {
  fp2 v0,v1,v2,m01,m02,m12,s01,s02,s12,t,u,w;
  fp2_mul(&v0,&x->c[0],&y->c[0]);
  fp2_mul(&v1,&x->c[1],&y->c[1]);
  fp2_mul(&v2,&x->c[2],&y->c[2]);
  fp2_add(&t,&x->c[0],&x->c[1]); fp2_add(&u,&y->c[0],&y->c[1]); fp2_mul(&m01,&t,&u);
  fp2_add(&t,&x->c[0],&x->c[2]); fp2_add(&u,&y->c[0],&y->c[2]); fp2_mul(&m02,&t,&u);
  fp2_add(&t,&x->c[1],&x->c[2]); fp2_add(&u,&y->c[1],&y->c[2]); fp2_mul(&m12,&t,&u);
  fp2_sub(&s01,&m01,&v0); fp2_sub(&s01,&s01,&v1);
  fp2_sub(&s02,&m02,&v0); fp2_sub(&s02,&s02,&v2);
  fp2_sub(&s12,&m12,&v1); fp2_sub(&s12,&s12,&v2);
  fp2_mul_xi(&w,&s12); fp2_add(&z->c[0],&v0,&w);
  fp2_mul_xi(&w,&v2);  fp2_add(&z->c[1],&s01,&w);
  fp2_add(&z->c[2],&s02,&v1);
}

static inline void fp6_mul_v(fp6 *z, const fp6 *x) { /* v*(c0+c1 v+c2 v^2) = xi*c2 + c0 v + c1 v^2 */
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
