#ifndef FP12_H
#define FP12_H
#include <stdint.h>
typedef struct { uint64_t v[4]; } fp;
typedef struct { fp a, b; } fp2;
typedef struct { fp2 c[3]; } fp6;
typedef struct { fp6 c[2]; } fp12;
extern const fp FP_P, FP_ZERO, FP_ONE_MONT;
void fp_mul(fp *z, const fp *x, const fp *y);
void fp_add(fp *z, const fp *x, const fp *y);
void fp_sub(fp *z, const fp *x, const fp *y);
int  fp_eq(const fp *x, const fp *y);
void fp2_mul(fp2 *z, const fp2 *x, const fp2 *y);
void fp2_add(fp2 *z, const fp2 *x, const fp2 *y);
void fp2_sub(fp2 *z, const fp2 *x, const fp2 *y);
void fp6_mul(fp6 *z, const fp6 *x, const fp6 *y);
void fp6_add(fp6 *z, const fp6 *x, const fp6 *y);
void fp6_sub(fp6 *z, const fp6 *x, const fp6 *y);
void fp12_mul(fp12 *z, const fp12 *x, const fp12 *y);
void fp12_add(fp12 *z, const fp12 *x, const fp12 *y);
void fp12_sub(fp12 *z, const fp12 *x, const fp12 *y);
int  fp12_eq(const fp12 *x, const fp12 *y);
void fp12_pow(fp12 *z, const fp12 *x, const uint64_t e[4]);
void fp12_set_one(fp12 *z);
#endif
