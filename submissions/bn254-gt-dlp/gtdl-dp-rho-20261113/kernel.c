/* BN254 GT multiplication kernel (CIOS), scalar 64-bit limbs only.
 * Counted like the challenge reference costs: gcc -O3 -fno-tree-vectorize
 * -fno-builtin, valgrind callgrind event Ir, loops of 4096 and 8192.
 * "t" mode prints the 96 words of a fixed product for external checking.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
typedef struct { uint64_t v[4]; } fp;

static const uint64_t P[4] = {0x3c208c16d87cfd47ULL, 0x97816a916871ca8dULL,
                              0xb85045b68181585dULL, 0x30644e72e131a029ULL};
static uint64_t N0;

/* Derived from tracks/bn254-gt-dlp/target-profile.json:
   GCOEF[i] = coefficient i of generator_g, in Montgomery form (x*2^256 mod p);
   REXP = subgroup order r; RMONT = 2^256 mod p; NINEMONT = 9*2^256 mod p. */
static const uint64_t GCOEF[12][4] = {
  {0x53e1d9fc3a8329efULL,0x9254a1949ff5465fULL,0x3d01af561fad5084ULL,0x2ae045056b3b7c1eULL},
  {0xd514f44343aae165ULL,0x6d330380ed14f311ULL,0xe6490e6130d72457ULL,0x0e220885bfce1f4dULL},
  {0x7f138852111c85f6ULL,0x6e42d4f25e64fff5ULL,0x1e2c29f438c5faf0ULL,0x05a226123d0a333dULL},
  {0x0221eab5ce4d20adULL,0x69b3e161979120e2ULL,0xea0115f3d02e5cd4ULL,0x18d86b9b7ecf1950ULL},
  {0x5869976aa182adffULL,0xb87c7bb0c0a87668ULL,0xf7d1d1a0ddced101ULL,0x07918bfca0015959ULL},
  {0x0be63e0e6309bae7ULL,0x2efbc923df47f05fULL,0xf5f5b66cc5d7e7f8ULL,0x253dac9cb214fa73ULL},
  {0x2e02a64acbd60549ULL,0xd618018ea58e4addULL,0x14d585f1a45ba647ULL,0x1832226987c434fcULL},
  {0x9458abcb56d24998ULL,0xb17540bd2a9e5adbULL,0x9a9983c82e401a9fULL,0x1614817a84c16291ULL},
  {0x172d1f257a4d598eULL,0xddf5bc7b7ffb5ac0ULL,0xae0b22c0bbb0f602ULL,0x1b158f3c2fae9b18ULL},
  {0x7002907c28ebfe11ULL,0x7b0591d3d080da67ULL,0xde7e5aa2181f138eULL,0x210e437dfc43d951ULL},
  {0xf16c96d081754cdbULL,0xce0394312bceeb55ULL,0x644e4dcf1f01ff0aULL,0x0cbea85ee0b236ccULL},
  {0xd34bab373157aa84ULL,0x3511ed44fd0d8598ULL,0x67e42a0bc2ced972ULL,0x2b8f1d5dfd20c55bULL},
};
static const uint64_t REXP[4] = {0x43e1f593f0000001ULL,0x2833e84879b97091ULL,0xb85045b68181585dULL,0x30644e72e131a029ULL };
static const uint64_t RMONT[4] = {0xd35d438dc58f0d9dULL,0x0a78eb28f5c70b3dULL,0x666ea36f7879462cULL,0x0e0a77c19a07df2fULL };
static const uint64_t NINEMONT[4] = {0xf60647ce410d7ff7ULL,0x2f3d6f4dd31bd011ULL,0x2943337e3940c6d1ULL,0x1d9598e8a7e39857ULL };


static inline void fadd(fp *r, const fp *a, const fp *b) {
    uint64_t t0, t1, t2, t3, c;
    u128 w;
    w = (u128)a->v[0] + b->v[0]; t0 = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)a->v[1] + b->v[1] + c; t1 = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)a->v[2] + b->v[2] + c; t2 = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)a->v[3] + b->v[3] + c; t3 = (uint64_t)w;
    uint64_t d = 0, s0, s1, s2, s3;
    w = (u128)t0 - P[0]; s0 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t1 - P[1] - d; s1 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t2 - P[2] - d; s2 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t3 - P[3] - d; s3 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    uint64_t m = (uint64_t)0 - (uint64_t)(d == 0);
    r->v[0] = (s0 & m) | (t0 & ~m); r->v[1] = (s1 & m) | (t1 & ~m);
    r->v[2] = (s2 & m) | (t2 & ~m); r->v[3] = (s3 & m) | (t3 & ~m);
}

static inline void fsub(fp *r, const fp *a, const fp *b) {
    uint64_t t0, t1, t2, t3, d;
    u128 w;
    w = (u128)a->v[0] - b->v[0]; t0 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)a->v[1] - b->v[1] - d; t1 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)a->v[2] - b->v[2] - d; t2 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)a->v[3] - b->v[3] - d; t3 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    uint64_t m = (uint64_t)0 - d, c;
    w = (u128)t0 + (P[0] & m); r->v[0] = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)t1 + (P[1] & m) + c; r->v[1] = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)t2 + (P[2] & m) + c; r->v[2] = (uint64_t)w; c = (uint64_t)(w >> 64);
    w = (u128)t3 + (P[3] & m) + c; r->v[3] = (uint64_t)w;
}

/* Montgomery CIOS: r = a*b*R^{-1} mod P, a,b < P, fully register-resident */
static inline void fmul(fp *r, const fp *a, const fp *b) {
    const uint64_t *A = a->v;
    uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t bi = b->v[i], c, m;
        u128 w;
        w = (u128)A[0] * bi + t0;        t0 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)A[1] * bi + t1 + c;    t1 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)A[2] * bi + t2 + c;    t2 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)A[3] * bi + t3 + c;    t3 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)t4 + c;                t4 = (uint64_t)w;
        m = t0 * N0;
        w = (u128)m * P[0] + t0;         c = (uint64_t)(w >> 64);
        w = (u128)m * P[1] + t1 + c;     t0 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)m * P[2] + t2 + c;     t1 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)m * P[3] + t3 + c;     t2 = (uint64_t)w; c = (uint64_t)(w >> 64);
        w = (u128)t4 + c;                t3 = (uint64_t)w; t4 = (uint64_t)(w >> 64);
    }
    uint64_t d, s0, s1, s2, s3;
    u128 w;
    w = (u128)t0 - P[0];         s0 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t1 - P[1] - d;     s1 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t2 - P[2] - d;     s2 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    w = (u128)t3 - P[3] - d;     s3 = (uint64_t)w; d = (uint64_t)((w >> 64) & 1);
    uint64_t m = (uint64_t)0 - (uint64_t)(d == 0);
    r->v[0] = (s0 & m) | (t0 & ~m); r->v[1] = (s1 & m) | (t1 & ~m);
    r->v[2] = (s2 & m) | (t2 & ~m); r->v[3] = (s3 & m) | (t3 & ~m);
}

typedef struct { fp c0, c1; } fp2;
typedef struct { fp2 c[3]; } fp6;
typedef struct { fp6 c0, c1; } fp12;

static inline void f2add(fp2 *r, const fp2 *a, const fp2 *b) { fadd(&r->c0, &a->c0, &b->c0); fadd(&r->c1, &a->c1, &b->c1); }
static inline void f2sub(fp2 *r, const fp2 *a, const fp2 *b) { fsub(&r->c0, &a->c0, &b->c0); fsub(&r->c1, &a->c1, &b->c1); }

static inline void f2mul(fp2 *r, const fp2 *a, const fp2 *b) {
    fp t0, t1, s0, s1;
    fmul(&t0, &a->c0, &b->c0);
    fmul(&t1, &a->c1, &b->c1);
    fadd(&s0, &a->c0, &a->c1);
    fadd(&s1, &b->c0, &b->c1);
    fmul(&s0, &s0, &s1);
    fsub(&s0, &s0, &t0);
    fsub(&s0, &s0, &t1);
    fsub(&r->c0, &t0, &t1);
    r->c1 = s0;
}

static inline void mul_by_xi(fp2 *r, const fp2 *a) { /* (a0+a1u)(9+u) = (9a0-a1) + (9a1+a0)u */
    fp2 d, t, au;
    f2add(&d, a, a); f2add(&d, &d, &d); f2add(&d, &d, &d);   /* 8a */
    f2add(&t, &d, a);                                        /* 9a */
    au.c0.v[0] = au.c0.v[1] = au.c0.v[2] = au.c0.v[3] = 0;
    au.c0.v[0] = 0;
    fsub(&au.c0, &au.c0, &a->c1);
    au.c1 = a->c0;
    f2add(r, &t, &au);
}

static inline void f6add(fp6 *r, const fp6 *a, const fp6 *b) { for (int i = 0; i < 3; i++) f2add(&r->c[i], &a->c[i], &b->c[i]); }
static inline void f6sub(fp6 *r, const fp6 *a, const fp6 *b) { for (int i = 0; i < 3; i++) f2sub(&r->c[i], &a->c[i], &b->c[i]); }

static inline void f6mul(fp6 *r, const fp6 *a, const fp6 *b) {
    fp2 v0, v1, v2, t0, t1, t2, s0, s1, x;
    f2mul(&v0, &a->c[0], &b->c[0]);
    f2mul(&v1, &a->c[1], &b->c[1]);
    f2mul(&v2, &a->c[2], &b->c[2]);
    f2add(&s0, &a->c[1], &a->c[2]); f2add(&s1, &b->c[1], &b->c[2]);
    f2mul(&t0, &s0, &s1); f2sub(&t0, &t0, &v1); f2sub(&t0, &t0, &v2);
    f2add(&s0, &a->c[0], &a->c[1]); f2add(&s1, &b->c[0], &b->c[1]);
    f2mul(&t1, &s0, &s1); f2sub(&t1, &t1, &v0); f2sub(&t1, &t1, &v1);
    f2add(&s0, &a->c[0], &a->c[2]); f2add(&s1, &b->c[0], &b->c[2]);
    f2mul(&t2, &s0, &s1); f2sub(&t2, &t2, &v0); f2sub(&t2, &t2, &v2);
    mul_by_xi(&x, &t0); f2add(&r->c[0], &v0, &x);
    mul_by_xi(&x, &v2); f2add(&r->c[1], &t1, &x);
    f2add(&r->c[2], &t2, &v1);
}

static inline void mul_by_v(fp6 *r, const fp6 *a) {
    fp2 x;
    mul_by_xi(&x, &a->c[2]);
    r->c[0] = x; r->c[1] = a->c[0]; r->c[2] = a->c[1];
}

static inline void f12mul(fp12 *r, const fp12 *a, const fp12 *b) {
    fp6 ac, bd, s0, s1, t, vbd;
    f6mul(&ac, &a->c0, &b->c0);
    f6mul(&bd, &a->c1, &b->c1);
    f6add(&s0, &a->c0, &a->c1);
    f6add(&s1, &b->c0, &b->c1);
    f6mul(&t, &s0, &s1);
    f6sub(&t, &t, &ac);
    f6sub(&t, &t, &bd);
    mul_by_v(&vbd, &bd);
    f6add(&r->c0, &ac, &vbd);
    r->c1 = t;
}

static uint64_t seed_state = 0x243f6a8885a308d3ULL;
static inline uint64_t xr(void) { seed_state ^= seed_state << 13; seed_state ^= seed_state >> 7; seed_state ^= seed_state << 17; return seed_state; }

/* w-basis (profile) -> tower: for m=0..5, j=m&1, k=m>>1,
   a = (c_m + 9 c_{m+6}) + c_{m+6} u  (u = w^6-9). */
static void from_profile(fp12 *z, const uint64_t c[12][4]) {
    for (int i = 0; i < 48; i++) ((uint64_t *)z)[i] = 0;
    for (int m = 0; m < 6; m++) {
        int j = m & 1, k = m >> 1;
        fp nine, cbeta, cm, ninec, re;
        for (int t = 0; t < 4; t++) nine.v[t] = NINEMONT[t];
        for (int t = 0; t < 4; t++) cbeta.v[t] = c[m + 6][t];
        for (int t = 0; t < 4; t++) cm.v[t] = c[m][t];
        fmul(&ninec, &nine, &cbeta);
        fadd(&re, &cm, &ninec);
        fp2 a; a.c0 = re; a.c1 = cbeta;
        if (j == 0) z->c0.c[k] = a; else z->c1.c[k] = a;
    }
}
static int check_g(void) {
    fp12 base, one, acc;
    from_profile(&base, GCOEF);
    for (int i = 0; i < 48; i++) ((uint64_t *)&one)[i] = 0;
    for (int t = 0; t < 4; t++) ((uint64_t *)&one)[t] = RMONT[t];
    acc = one;
    int top = 0;
    for (int i = 255; i >= 0; i--) if ((REXP[i / 64] >> (i % 64)) & 1) { top = i; break; }
    for (int i = top; i >= 0; i--) {
        fp12 t; f12mul(&t, &acc, &acc); acc = t;
        if ((REXP[i / 64] >> (i % 64)) & 1) { fp12 u; f12mul(&u, &acc, &base); acc = u; }
    }
    int eq = 1, nz = 0;
    for (int i = 0; i < 48; i++) {
        if (((uint64_t *)&acc)[i] != ((uint64_t *)&one)[i]) eq = 0;
        if (((uint64_t *)&base)[i] != ((uint64_t *)&one)[i]) nz = 1;
    }
    printf("g^r == 1: %s ;  g != 1: %s\n", eq ? "YES" : "NO", nz ? "YES" : "NO");
    return eq && nz;
}

int main(int argc, char **argv) {
    uint64_t n0 = 1;
    for (int i = 0; i < 6; i++) n0 *= 2 - P[0] * n0;
    N0 = (uint64_t)0 - n0;


    if (argc > 1 && argv[1][0] == 'g') return check_g() ? 0 : 1;
    if (argc > 1 && argv[1][0] == 't') {
        fp12 x, y, z;
        uint64_t *px = (uint64_t *)&x, *py = (uint64_t *)&y;
        for (int i = 0; i < 48; i++) {
            px[i] = (uint64_t)((i * 0x9e3779b97f4a7c15ULL + 12345) & 0x0fffffffffffffffULL);
            py[i] = (uint64_t)((i * 0xc2b2ae3d27d4eb4fULL + 7) & 0x0fffffffffffffffULL);
        }
        f12mul(&z, &x, &y);
        uint64_t *pz = (uint64_t *)&z;
        for (int i = 0; i < 48; i++) printf("%016llx\n", (unsigned long long)pz[i]);
        return 0;
    }
    long iters = argc > 1 ? atol(argv[1]) : 4096;
    fp12 a, b;
    uint64_t *pa = (uint64_t *)&a;
    for (int i = 0; i < 48; i++) pa[i] = xr() & 0x0fffffffffffffffULL;
    b = a; b.c1.c[2].c0.v[0] ^= 1ULL;
    volatile uint64_t sink = 0;
    for (long i = 0; i < iters; i++) {
        fp12 c; f12mul(&c, &a, &b);
        b = c;
        sink ^= c.c0.c[0].c0.v[0];
    }
    printf("%llu\n", (unsigned long long)sink);
    return 0;
}
