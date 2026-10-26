#include <stdint.h>
#include "la_consts.h"
typedef struct { uint64_t v[4]; } fe_t;
typedef struct { uint64_t v[5]; } acc_t;
/* reduce a 5-limb accumulator v < 2^275 into [0, r): t = v>>255 (< 2^20),
 * qhat = (t*MBAR)>>20 with MBAR = floor(2^275/r); v -= qhat*r; then <= 2 conditional subtracts. */
static void red(uint64_t o[4], const acc_t *a) {
    uint64_t v0=a->v[0],v1=a->v[1],v2=a->v[2],v3=a->v[3],v4=a->v[4];
    uint64_t t = (v4 << 1) | (v3 >> 63);
    uint64_t q = (uint64_t)(((__uint128_t)t * MBAR) >> 20);
    __uint128_t c, b; uint64_t r4[5];
    c = (__uint128_t)q*R0; r4[0]=(uint64_t)c; c>>=64;
    c += (__uint128_t)q*R1; r4[1]=(uint64_t)c; c>>=64;
    c += (__uint128_t)q*R2; r4[2]=(uint64_t)c; c>>=64;
    c += (__uint128_t)q*R3; r4[3]=(uint64_t)c; r4[4]=(uint64_t)(c>>64);
    b = (__uint128_t)v0 - r4[0]; v0=(uint64_t)b; uint64_t br=(uint64_t)(b>>64)&1;
    b = (__uint128_t)v1 - r4[1] - br; v1=(uint64_t)b; br=(uint64_t)(b>>64)&1;
    b = (__uint128_t)v2 - r4[2] - br; v2=(uint64_t)b; br=(uint64_t)(b>>64)&1;
    b = (__uint128_t)v3 - r4[3] - br; v3=(uint64_t)b; br=(uint64_t)(b>>64)&1;
    v4 = v4 - r4[4] - br;
    for (int k = 0; k < 3; k++) {      /* v in [0, 3r) here */
        b = (__uint128_t)v0 - R0; uint64_t w0=(uint64_t)b; br=(uint64_t)(b>>64)&1;
        b = (__uint128_t)v1 - R1 - br; uint64_t w1=(uint64_t)b; br=(uint64_t)(b>>64)&1;
        b = (__uint128_t)v2 - R2 - br; uint64_t w2=(uint64_t)b; br=(uint64_t)(b>>64)&1;
        b = (__uint128_t)v3 - R3 - br; uint64_t w3=(uint64_t)b; br=(uint64_t)(b>>64)&1;
        uint64_t w4 = v4 - br;
        uint64_t keep = w4 >> 63;          /* borrow out: v < r, keep v */
        uint64_t mk = keep - 1;            /* all-ones if we take w */
        v0 = (w0&mk)|(v0&~mk); v1=(w1&mk)|(v1&~mk); v2=(w2&mk)|(v2&~mk); v3=(w3&mk)|(v3&~mk); v4=(w4&mk)|(v4&~mk);
    }
    o[0]=v0;o[1]=v1;o[2]=v2;o[3]=v3;
}
/* y = (pos - neg) mod r */
void row_finish(fe_t *y, const acc_t *pos, const acc_t *neg) {
    uint64_t p[4], n[4];
    red(p, pos); red(n, neg);
    __uint128_t b; uint64_t br, d[4];
    b=(__uint128_t)p[0]-n[0]; d[0]=(uint64_t)b; br=(uint64_t)(b>>64)&1;
    b=(__uint128_t)p[1]-n[1]-br; d[1]=(uint64_t)b; br=(uint64_t)(b>>64)&1;
    b=(__uint128_t)p[2]-n[2]-br; d[2]=(uint64_t)b; br=(uint64_t)(b>>64)&1;
    b=(__uint128_t)p[3]-n[3]-br; d[3]=(uint64_t)b; br=(uint64_t)(b>>64)&1;
    uint64_t mk = 0 - br;   /* add r back if negative */
    __uint128_t c;
    c=(__uint128_t)d[0]+(R0&mk); y->v[0]=(uint64_t)c; c>>=64;
    c+=(__uint128_t)d[1]+(R1&mk); y->v[1]=(uint64_t)c; c>>=64;
    c+=(__uint128_t)d[2]+(R2&mk); y->v[2]=(uint64_t)c; c>>=64;
    c+=(__uint128_t)d[3]+(R3&mk); y->v[3]=(uint64_t)c;
}
