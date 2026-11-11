/* Scalar 64-bit reference implementation of the priced walk arithmetic.
 * Every function below is straight-line (no loops, no data-dependent
 * branches).  Leaf field operations are separate calls; the Fp2/Fp6/Fp12
 * tower functions are fixed sequences of calls.  The executed instruction
 * count is therefore the static count of each body plus its callees; see
 * machine.md and development-evidence.json.  __uint128_t expands to
 * MUL/UMULH and ADDS/ADCS pairs on 64-bit registers; no SIMD is used.
 */
#include <stdint.h>
typedef uint64_t word;
typedef __uint128_t wide;
typedef struct { word l[4]; } fe;          /* Montgomery residue < m      */
typedef struct { fe c0, c1; } fp2;          /* c0 + c1*u, u^2 = -1         */
typedef struct { fp2 c0, c1, c2; } fp6;     /* over Fp2, v^3 = xi = 9+u    */
typedef struct { fp6 c0, c1; } fp12;        /* over Fp6, w^2 = v           */
#define NOINLINE __attribute__((noinline))
#define INLINE static inline __attribute__((always_inline))

static const word P[4] = {0x3c208c16d87cfd47ULL, 0x97816a916871ca8dULL,
                          0xb85045b68181585dULL, 0x30644e72e131a029ULL};
static const word RR[4] = {0x43e1f593f0000001ULL, 0x2833e84879b97091ULL,
                           0xb85045b68181585dULL, 0x30644e72e131a029ULL};
#define P_N0 0x87d20782e4866389ULL          /* -p^{-1} mod 2^64 */

/* out = a+b mod m for canonical a,b; 2m < 2^255 so no carry out. */
INLINE void add_m(fe *o, const fe *a, const fe *b, const word *m) {
    unsigned long long c, s0,s1,s2,s3,d0,d1,d2,d3;
    s0=__builtin_addcll(a->l[0],b->l[0],0,&c); s1=__builtin_addcll(a->l[1],b->l[1],c,&c);
    s2=__builtin_addcll(a->l[2],b->l[2],c,&c); s3=__builtin_addcll(a->l[3],b->l[3],c,&c);
    d0=__builtin_subcll(s0,m[0],0,&c); d1=__builtin_subcll(s1,m[1],c,&c);
    d2=__builtin_subcll(s2,m[2],c,&c); d3=__builtin_subcll(s3,m[3],c,&c);
    o->l[0]=c?s0:d0; o->l[1]=c?s1:d1; o->l[2]=c?s2:d2; o->l[3]=c?s3:d3;
}
NOINLINE void fp_add(fe *o, const fe *a, const fe *b) { add_m(o,a,b,P); }
NOINLINE void fr_add(fe *o, const fe *a, const fe *b) { add_m(o,a,b,RR); }
/* out = a-b mod p */
NOINLINE void fp_sub(fe *o, const fe *a, const fe *b) {
    unsigned long long c, d0,d1,d2,d3,k;
    d0=__builtin_subcll(a->l[0],b->l[0],0,&c); d1=__builtin_subcll(a->l[1],b->l[1],c,&c);
    d2=__builtin_subcll(a->l[2],b->l[2],c,&c); d3=__builtin_subcll(a->l[3],b->l[3],c,&c);
    k=0-c;
    o->l[0]=__builtin_addcll(d0,P[0]&k,0,&c); o->l[1]=__builtin_addcll(d1,P[1]&k,c,&c);
    o->l[2]=__builtin_addcll(d2,P[2]&k,c,&c); o->l[3]=__builtin_addcll(d3,P[3]&k,c,&c);
}
/* Montgomery product a*b*2^-256 mod p, coarsely integrated operand scanning,
 * fully unrolled.  Inputs canonical (<p<2^254); result canonical. */
#define MULROW(bi) do { wide v; word c, q; \
    v=(wide)a->l[0]*(bi)+t0; t0=(word)v; c=(word)(v>>64); \
    v=(wide)a->l[1]*(bi)+t1+c; t1=(word)v; c=(word)(v>>64); \
    v=(wide)a->l[2]*(bi)+t2+c; t2=(word)v; c=(word)(v>>64); \
    v=(wide)a->l[3]*(bi)+t3+c; t3=(word)v; c=(word)(v>>64); \
    v=(wide)t4+c; t4=(word)v; t5=(word)(v>>64); \
    q=t0*P_N0; \
    v=(wide)q*P[0]+t0; c=(word)(v>>64); \
    v=(wide)q*P[1]+t1+c; t0=(word)v; c=(word)(v>>64); \
    v=(wide)q*P[2]+t2+c; t1=(word)v; c=(word)(v>>64); \
    v=(wide)q*P[3]+t3+c; t2=(word)v; c=(word)(v>>64); \
    v=(wide)t4+c; t3=(word)v; t4=t5+(word)(v>>64); } while (0)
NOINLINE void fp_mul(fe *o, const fe *a, const fe *b) {
    word t0=0,t1=0,t2=0,t3=0,t4=0,t5=0; unsigned long long d0,d1,d2,d3,br;
    MULROW(b->l[0]); MULROW(b->l[1]); MULROW(b->l[2]); MULROW(b->l[3]);
    d0=__builtin_subcll(t0,P[0],0,&br); d1=__builtin_subcll(t1,P[1],br,&br);
    d2=__builtin_subcll(t2,P[2],br,&br); d3=__builtin_subcll(t3,P[3],br,&br);
    o->l[0]=br?t0:d0; o->l[1]=br?t1:d1; o->l[2]=br?t2:d2; o->l[3]=br?t3:d3;
}

/* ---- Fp2: 3 mul + 5 add/sub (Karatsuba) ---- */
NOINLINE void fp2_mul(fp2 *o, const fp2 *a, const fp2 *b) {
    fe s, t, v0, v1, m;
    fp_add(&s,&a->c0,&a->c1); fp_add(&t,&b->c0,&b->c1);
    fp_mul(&v0,&a->c0,&b->c0); fp_mul(&v1,&a->c1,&b->c1); fp_mul(&m,&s,&t);
    fp_sub(&o->c0,&v0,&v1);
    fp_sub(&m,&m,&v0); fp_sub(&o->c1,&m,&v1);
}
NOINLINE void fp2_add(fp2 *o, const fp2 *a, const fp2 *b) {
    fp_add(&o->c0,&a->c0,&b->c0); fp_add(&o->c1,&a->c1,&b->c1); }
NOINLINE void fp2_sub(fp2 *o, const fp2 *a, const fp2 *b) {
    fp_sub(&o->c0,&a->c0,&b->c0); fp_sub(&o->c1,&a->c1,&b->c1); }
/* o = (9+u)*a = (9x-y) + (9y+x)u : 10 Fp add/sub.  o must not alias a. */
NOINLINE void fp2_mul_xi(fp2 *o, const fp2 *a) {
    fe x8, y8;
    fp_add(&x8,&a->c0,&a->c0); fp_add(&x8,&x8,&x8); fp_add(&x8,&x8,&x8);
    fp_add(&y8,&a->c1,&a->c1); fp_add(&y8,&y8,&y8); fp_add(&y8,&y8,&y8);
    fp_add(&x8,&x8,&a->c0); fp_add(&y8,&y8,&a->c1);   /* 9x, 9y */
    fp_sub(&o->c0,&x8,&a->c1); fp_add(&o->c1,&y8,&a->c0);
}
/* ---- Fp6: 6 Fp2 mul, 15 Fp2 add/sub, 2 xi (Karatsuba) ---- */
NOINLINE void fp6_mul(fp6 *o, const fp6 *a, const fp6 *b) {
    fp2 v0,v1,v2,s,t,m,x;
    fp2_mul(&v0,&a->c0,&b->c0); fp2_mul(&v1,&a->c1,&b->c1); fp2_mul(&v2,&a->c2,&b->c2);
    fp2_add(&s,&a->c1,&a->c2); fp2_add(&t,&b->c1,&b->c2); fp2_mul(&m,&s,&t);
    fp2_sub(&m,&m,&v1); fp2_sub(&m,&m,&v2); fp2_mul_xi(&x,&m);
    fp2_add(&s,&a->c0,&a->c1); fp2_add(&t,&b->c0,&b->c1); fp2_mul(&m,&s,&t);
    fp2_add(&s,&a->c0,&a->c2); fp2_add(&t,&b->c0,&b->c2);
    fp2_add(&o->c0,&v0,&x);   /* all reads of a,b precede this first write: o may alias */
    fp2_sub(&m,&m,&v0); fp2_sub(&m,&m,&v1); fp2_mul_xi(&x,&v2);
    fp2_add(&o->c1,&m,&x);
    fp2_mul(&m,&s,&t);
    fp2_sub(&m,&m,&v0); fp2_sub(&m,&m,&v2); fp2_add(&o->c2,&m,&v1);
}
NOINLINE void fp6_add(fp6 *o, const fp6 *a, const fp6 *b) {
    fp2_add(&o->c0,&a->c0,&b->c0); fp2_add(&o->c1,&a->c1,&b->c1); fp2_add(&o->c2,&a->c2,&b->c2); }
NOINLINE void fp6_sub(fp6 *o, const fp6 *a, const fp6 *b) {
    fp2_sub(&o->c0,&a->c0,&b->c0); fp2_sub(&o->c1,&a->c1,&b->c1); fp2_sub(&o->c2,&a->c2,&b->c2); }
/* ---- Fp12: 3 Fp6 mul, 4 Fp6 add/sub, 3 Fp2 add, 1 xi.  All reads of a,b
 * precede the first write to o. ---- */
NOINLINE void fp12_mul(fp12 *o, const fp12 *a, const fp12 *b) {
    fp6 v0,v1,s,t; fp2 x;
    fp6_mul(&v0,&a->c0,&b->c0); fp6_mul(&v1,&a->c1,&b->c1);
    fp6_add(&s,&a->c0,&a->c1); fp6_add(&t,&b->c0,&b->c1);
    fp6_mul(&o->c1,&s,&t);
    fp6_sub(&o->c1,&o->c1,&v0); fp6_sub(&o->c1,&o->c1,&v1);
    /* c0 = v0 + v*v1, v*(c0,c1,c2) = (xi*c2, c0, c1) */
    fp2_mul_xi(&x,&v1.c2);
    fp2_add(&o->c0.c0,&v0.c0,&x); fp2_add(&o->c0.c1,&v0.c1,&v1.c0); fp2_add(&o->c0.c2,&v0.c2,&v1.c1);
}

/* ---- one walk step ---- */
typedef struct { fp12 m; fe c, d; } record;            /* 448 bytes */
typedef struct { fp12 y; fe a, b; } state;
#define K1 0x9E3779B97F4A7C15ULL
#define K2 0xD6E8FEB86659FD93ULL
INLINE word mix(const fp12 *y, word salt) {
    word h = (y->c0.c0.c0.l[0] ^ salt) * K1;
    h ^= h >> 32; h += y->c1.c1.c1.l[1]; h *= K2; h ^= h >> 32;
    return h;
}
/* returns 1 when the new state is distinguished (theta = 2^-100) */
NOINLINE int walk_step(state *out, const state *in, const record *table, word salt) {
    const record *e = &table[mix(&in->y, salt) >> 54];      /* 1024 records */
    fp12_mul(&out->y, &in->y, &e->m);
    fr_add(&out->a, &in->a, &e->c); fr_add(&out->b, &in->b, &e->d);
    return (out->y.c0.c1.c0.l[0] | (out->y.c1.c0.c1.l[0] & 0xFFFFFFFFFULL)) == 0;
}
