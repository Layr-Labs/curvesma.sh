#ifndef FP12_H
#define FP12_H
#include <stdint.h>
#include <string.h>
typedef struct { uint64_t v[4]; } fp;
typedef struct { fp a, b; } fp2;
typedef struct { fp2 c[3]; } fp6;
typedef struct { fp6 c[2]; } fp12;
extern const fp FP_P, FP_ZERO, FP_ONE_MONT;
void fp_mul(fp *z, const fp *x, const fp *y);

/* Branchless canonical addition and subtraction; inputs and outputs are
 * canonical (< p) four-limb values. Inline because the tower multiplication
 * makes far more calls to these than to fp_mul. */
static inline void fp_add(fp *z, const fp *x, const fp *y) {
  uint64_t s[4], d[4];
  __uint128_t acc;
  uint64_t carry = 0;
  for (int i = 0; i < 4; i++) {
    acc = (__uint128_t)x->v[i] + y->v[i] + carry;
    s[i] = (uint64_t)acc; carry = (uint64_t)(acc >> 64);
  }
  uint64_t borrow = 0;
  for (int i = 0; i < 4; i++) {
    acc = (__uint128_t)s[i] - FP_P.v[i] - borrow;
    d[i] = (uint64_t)acc; borrow = (uint64_t)(acc >> 64) & 1ULL;
  }
  uint64_t cond = carry | (borrow ^ 1ULL);
  uint64_t mask = 0 - cond;
  for (int i = 0; i < 4; i++) z->v[i] = (d[i] & mask) | (s[i] & ~mask);
}

static inline void fp_sub(fp *z, const fp *x, const fp *y) {
  uint64_t d[4];
  __uint128_t acc;
  uint64_t borrow = 0;
  for (int i = 0; i < 4; i++) {
    acc = (__uint128_t)x->v[i] - y->v[i] - borrow;
    d[i] = (uint64_t)acc; borrow = (uint64_t)(acc >> 64) & 1ULL;
  }
  uint64_t mask = 0 - borrow;
  uint64_t carry = 0;
  for (int i = 0; i < 4; i++) {
    acc = (__uint128_t)d[i] + (FP_P.v[i] & mask) + carry;
    z->v[i] = (uint64_t)acc; carry = (uint64_t)(acc >> 64);
  }
}

static inline int fp_eq(const fp *x, const fp *y) {
  uint64_t d = 0;
  for (int i = 0; i < 4; i++) d |= x->v[i] ^ y->v[i];
  return d == 0;
}
static inline void fp_neg(fp *z, const fp *x) {
  fp zero; memset(&zero, 0, sizeof zero);
  fp_sub(z, &zero, x);
}

void fp2_mul(fp2 *z, const fp2 *x, const fp2 *y);
static inline void fp2_add(fp2 *z, const fp2 *x, const fp2 *y){ fp_add(&z->a,&x->a,&y->a); fp_add(&z->b,&x->b,&y->b); }
static inline void fp2_sub(fp2 *z, const fp2 *x, const fp2 *y){ fp_sub(&z->a,&x->a,&y->a); fp_sub(&z->b,&x->b,&y->b); }
void fp6_mul(fp6 *z, const fp6 *x, const fp6 *y);
/* Written out explicitly so the compiler cannot keep a loop back-edge. */
static inline void fp6_add(fp6 *z, const fp6 *x, const fp6 *y){
  fp2_add(&z->c[0],&x->c[0],&y->c[0]);
  fp2_add(&z->c[1],&x->c[1],&y->c[1]);
  fp2_add(&z->c[2],&x->c[2],&y->c[2]);
}
static inline void fp6_sub(fp6 *z, const fp6 *x, const fp6 *y){
  fp2_sub(&z->c[0],&x->c[0],&y->c[0]);
  fp2_sub(&z->c[1],&x->c[1],&y->c[1]);
  fp2_sub(&z->c[2],&x->c[2],&y->c[2]);
}
void fp12_mul(fp12 *z, const fp12 *x, const fp12 *y);
static inline void fp12_add(fp12 *z, const fp12 *x, const fp12 *y){ fp6_add(&z->c[0],&x->c[0],&y->c[0]); fp6_add(&z->c[1],&x->c[1],&y->c[1]); }
static inline void fp12_sub(fp12 *z, const fp12 *x, const fp12 *y){ fp6_sub(&z->c[0],&x->c[0],&y->c[0]); fp6_sub(&z->c[1],&x->c[1],&y->c[1]); }
int  fp12_eq(const fp12 *x, const fp12 *y);
void fp12_pow(fp12 *z, const fp12 *x, const uint64_t e[4]);
void fp12_set_one(fp12 *z);
#endif
