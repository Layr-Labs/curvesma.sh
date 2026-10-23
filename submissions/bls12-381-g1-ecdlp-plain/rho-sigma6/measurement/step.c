/* Batched order-6 class rho step on BLS12-381 G1 (sigma walk), for Callgrind instruction counts.
 *
 * Uses the organizer-pinned blst (commit f7a88067...), linked as in the organizer's
 * harness: server.c plus the upstream AArch64 assembly. Field elements stay in blst's
 * Montgomery form M(v) = v*2^384 mod p (six 64-bit limbs, always fully reduced).
 *
 * The six multiplications of the batched incomplete affine addition (3M batch, 3M
 * formula) go through ref_mul/ref_sqr so their instructions can be separated from
 * everything else; every other instruction of the step is charged.
 *
 * Classes and walk (see walk_ref.py, which replays this exactly):
 *  - class key: key(W) = M(y) if M(y) is even, else p - M(y)  (= M(+-y), even).
 *    y determines the class: y^2 = x^3 + 4 fixes x up to beta.
 *  - canonical point: Wc = (beta^k x, s y) with M(beta^k x) least of the three and
 *    s = +-1 making M(s y) even.
 *  - the walk stores the actual point W = sigma(Wc), sigma = (x -> beta^-k x, y -> s y),
 *    and steps W <- W + sigma(R_i), i = (key >> 1) mod 2^ADDER_BITS. Since
 *    W + sigma(R_i) = sigma(Wc + R_i), the sequence of classes is that of the walk
 *    Wc <- canon(Wc + R_i). sigma(R_i) is read from a table of the six images of R_i.
 *  - k is found without computing beta*x: the 64-bit fractions 2^64*M(x)/p and
 *    2^64*M(beta x)/p are approximated from M(x) with a few multiplications (error in
 *    [0,13)); 2^64*M(beta^2 x)/p follows from M(x)+M(beta x)+M(beta^2 x) in {p,2p}.
 *    When the error intervals do not separate the least value, k is computed exactly.
 *
 * Usage: step <inputs.txt> <batch-steps> [dump]
 *        step <inputs.txt> <count> single   (count unbatched f(W) evaluations)
 * Inputs: 2^ADDER_BITS adders then B start points, as written by walk_ref.py gen.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#include "/blst/bindings/blst.h"

typedef uint64_t limb_t;
typedef limb_t vec384[6];
/* blst AArch64 assembly entry points (build/elf/*-armv8.S). */
void mul_mont_384(vec384 ret, const vec384 a, const vec384 b, const vec384 p, limb_t n0);
void sqr_mont_384(vec384 ret, const vec384 a, const vec384 p, limb_t n0);
void add_mod_384(vec384 ret, const vec384 a, const vec384 b, const vec384 p);
void sub_mod_384(vec384 ret, const vec384 a, const vec384 b, const vec384 p);

#ifndef ADDER_BITS
#define ADDER_BITS 16
#endif
#define NA (1 << ADDER_BITS)
#ifndef B
#define B 16384
#endif
#define ANCHOR_PERIOD 64
#define MAX_CYCLE 64
#define DP_BITS 60
#define KEY_SHIFT 1 /* key is even: index uses bits 1..ADDER_BITS */

typedef struct { limb_t l[6]; } fe;
typedef struct { fe x, y; } pt; /* affine, Montgomery form */

static const vec384 P = {0xb9feffffffffaaab, 0x1eabfffeb153ffff, 0x6730d2a0f6b0f624,
                         0x64774b84f38512bf, 0x4b1ba7b6434bacd7, 0x1a0111ea397fe69a};
static const limb_t N0 = 0x89f3fffcfffcfffd;
/* beta = 0x5f19672f...fffefffe (cube root of unity in Fp), Montgomery form. */
static const fe BETA = {{0x30f1361b798a64e8, 0xf3b8ddab7ece5a2a, 0x16a8ca3ac61577f7,
                         0xc26a2ff874fd029b, 0x3636b76660701c6e, 0x051ba4ab241b6160}};
static const fe ZERO = {{0}};
/* C0 = floor(2^512/p) (limbs 1,2 used), C1 = floor(2^512*beta/p) (limbs 1..7 used). */
/* Not static: the compiler must load these rather than rebuild 64-bit immediates. */
limb_t FRAC_K[10] = {0x8d4792687dbc18bb, 0x32492fe34d42e1f4, 0xb35ac25664a8d148,
                     0x093138b09cce4e83, 0x51363cf7e4d25b79, 0x10b32afd5bc3c517,
                     0xa8360e5f04680424, 0x3, 0xd835d2f3cc9e45ce, 0x9};
#define C1 FRAC_K
#define C0_1 FRAC_K[8]
#define C0_2 FRAC_K[9]
#define FRAC_ERR 13

static pt T[NA];        /* adders R_i as read */
static pt T6[NA][3][2]; /* T6[i][k][s] = (beta^-k X_i, (-1)^s Y_i) */
static pt WA[B], WB[B];
static pt *W = WA, *WN = WB; /* current and next states, swapped every batch step */
static const pt *SEL[B];     /* sigma(R_i) for each walk's next step */
static fe D[B], PRE[B];      /* PRE[w] = D[0]*...*D[w-1] */
static uint64_t curkey[B], prevkey[B], anchorkey[B]; /* low key limbs: W, its predecessor, anchor */
static uint64_t two_cycles, long_cycles, dps, exceptional, unbatched, frac_fallbacks;

/* Reference-operation multiplications. The barrier after the call keeps it from
 * becoming a tail call, so --toggle-collect brackets exactly the multiplication. */
__attribute__((noinline)) static void ref_mul(fe *r, const fe *a, const fe *b) {
    blst_fp_mul((blst_fp *)r, (const blst_fp *)a, (const blst_fp *)b);
    __asm__ volatile("" ::: "memory");
}
__attribute__((noinline)) static void ref_sqr(fe *r, const fe *a) {
    blst_fp_sqr((blst_fp *)r, (const blst_fp *)a);
    __asm__ volatile("" ::: "memory");
}

static inline void fmul(fe *r, const fe *a, const fe *b) { mul_mont_384(r->l, a->l, b->l, P, N0); }
static inline void fsqr(fe *r, const fe *a) { sqr_mont_384(r->l, a->l, P, N0); }
static inline void fadd(fe *r, const fe *a, const fe *b) { add_mod_384(r->l, a->l, b->l, P); }

/* r = a - b mod p for reduced a, b: subtract, then add p back only on borrow. */
static inline void fsub(fe *r, const fe *a, const fe *b) {
    uint64_t t0, t1, t2, t3, t4, t5, u0, u1, u2, u3, u4, u5;
    __asm__ volatile(
        "ldp %[t0], %[t1], [%[a]]\n\t"
        "ldp %[u0], %[u1], [%[b]]\n\t"
        "ldp %[t2], %[t3], [%[a], #16]\n\t"
        "ldp %[u2], %[u3], [%[b], #16]\n\t"
        "ldp %[t4], %[t5], [%[a], #32]\n\t"
        "ldp %[u4], %[u5], [%[b], #32]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        "sbcs %[t4], %[t4], %[u4]\n\t"
        "sbcs %[t5], %[t5], %[u5]\n\t"
        "b.cs 1f\n\t"
        "ldp %[u0], %[u1], [%[p]]\n\t"
        "ldp %[u2], %[u3], [%[p], #16]\n\t"
        "ldp %[u4], %[u5], [%[p], #32]\n\t"
        "adds %[t0], %[t0], %[u0]\n\t"
        "adcs %[t1], %[t1], %[u1]\n\t"
        "adcs %[t2], %[t2], %[u2]\n\t"
        "adcs %[t3], %[t3], %[u3]\n\t"
        "adcs %[t4], %[t4], %[u4]\n\t"
        "adc %[t5], %[t5], %[u5]\n"
        "1:\n\t"
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        "stp %[t4], %[t5], [%[r], #32]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3), [t4] "=&r"(t4), [t5] "=&r"(t5),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3), [u4] "=&r"(u4), [u5] "=&r"(u5)
        : [r] "r"(r->l), [a] "r"(a->l), [b] "r"(b->l), [p] "r"(P)
        : "cc", "memory");
}

/* r = a - b - c mod p for reduced a, b, c. */
static inline void fsub2(fe *r, const fe *a, const fe *b, const fe *c) {
    uint64_t t0, t1, t2, t3, t4, t5, u0, u1, u2, u3, u4, u5;
    __asm__ volatile(
        "ldp %[t0], %[t1], [%[a]]\n\t"
        "ldp %[u0], %[u1], [%[b]]\n\t"
        "ldp %[t2], %[t3], [%[a], #16]\n\t"
        "ldp %[u2], %[u3], [%[b], #16]\n\t"
        "ldp %[t4], %[t5], [%[a], #32]\n\t"
        "ldp %[u4], %[u5], [%[b], #32]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        "sbcs %[t4], %[t4], %[u4]\n\t"
        "sbcs %[t5], %[t5], %[u5]\n\t"
        "b.cs 1f\n\t"
        "ldp %[u0], %[u1], [%[p]]\n\t"
        "ldp %[u2], %[u3], [%[p], #16]\n\t"
        "ldp %[u4], %[u5], [%[p], #32]\n\t"
        "adds %[t0], %[t0], %[u0]\n\t"
        "adcs %[t1], %[t1], %[u1]\n\t"
        "adcs %[t2], %[t2], %[u2]\n\t"
        "adcs %[t3], %[t3], %[u3]\n\t"
        "adcs %[t4], %[t4], %[u4]\n\t"
        "adc %[t5], %[t5], %[u5]\n"
        "1:\n\t"
        "ldp %[u0], %[u1], [%[c]]\n\t"
        "ldp %[u2], %[u3], [%[c], #16]\n\t"
        "ldp %[u4], %[u5], [%[c], #32]\n\t"
        "subs %[t0], %[t0], %[u0]\n\t"
        "sbcs %[t1], %[t1], %[u1]\n\t"
        "sbcs %[t2], %[t2], %[u2]\n\t"
        "sbcs %[t3], %[t3], %[u3]\n\t"
        "sbcs %[t4], %[t4], %[u4]\n\t"
        "sbcs %[t5], %[t5], %[u5]\n\t"
        "b.cs 2f\n\t"
        "ldp %[u0], %[u1], [%[p]]\n\t"
        "ldp %[u2], %[u3], [%[p], #16]\n\t"
        "ldp %[u4], %[u5], [%[p], #32]\n\t"
        "adds %[t0], %[t0], %[u0]\n\t"
        "adcs %[t1], %[t1], %[u1]\n\t"
        "adcs %[t2], %[t2], %[u2]\n\t"
        "adcs %[t3], %[t3], %[u3]\n\t"
        "adcs %[t4], %[t4], %[u4]\n\t"
        "adc %[t5], %[t5], %[u5]\n"
        "2:\n\t"
        "stp %[t0], %[t1], [%[r]]\n\t"
        "stp %[t2], %[t3], [%[r], #16]\n\t"
        "stp %[t4], %[t5], [%[r], #32]\n\t"
        : [t0] "=&r"(t0), [t1] "=&r"(t1), [t2] "=&r"(t2), [t3] "=&r"(t3), [t4] "=&r"(t4), [t5] "=&r"(t5),
          [u0] "=&r"(u0), [u1] "=&r"(u1), [u2] "=&r"(u2), [u3] "=&r"(u3), [u4] "=&r"(u4), [u5] "=&r"(u5)
        : [r] "r"(r->l), [a] "r"(a->l), [b] "r"(b->l), [c] "r"(c->l), [p] "r"(P)
        : "cc", "memory");
}

/* Integer comparison of Montgomery representatives: <0, 0, >0. */
static inline int fcmp(const fe *a, const fe *b) {
    for (int i = 5; i >= 0; i--) {
        if (a->l[i] != b->l[i]) return a->l[i] < b->l[i] ? -1 : 1;
    }
    return 0;
}
static inline int feq(const fe *a, const fe *b) { return memcmp(a->l, b->l, sizeof a->l) == 0; }
static inline int fis_zero(const fe *a) {
    return (a->l[0] | a->l[1] | a->l[2] | a->l[3] | a->l[4] | a->l[5]) == 0;
}
static inline void fneg(fe *r, const fe *a) { fsub(r, &ZERO, a); }

/* Exact canonical exponent k: index of the least of M(x), M(beta x), M(beta^2 x). */
static int canon_k_exact(const fe *x) {
    fe bx, s, b2x;
    fmul(&bx, x, &BETA);
    fadd(&s, x, &bx);
    fneg(&b2x, &s);
    int k = 0;
    const fe *m = x;
    if (fcmp(&bx, m) < 0) { k = 1; m = &bx; }
    if (fcmp(&b2x, m) < 0) k = 2;
    return k;
}

#define HI(a, b) ((uint64_t)(((unsigned __int128)(a) * (b)) >> 64))

/* k from 64-bit fraction approximations; exact fallback when not separated. */
static inline int canon_k(const fe *x) {
    const uint64_t a0 = x->l[0], a1 = x->l[1], a2 = x->l[2], a3 = x->l[3], a4 = x->l[4], a5 = x->l[5];
    /* g0, g1 in [f - 12, f], f = floor(2^64 * M(beta^j x) / p), barring wraparound. */
    uint64_t g0 = a5 * C0_2 + HI(a4, C0_2) + HI(a5, C0_1);
    uint64_t g1 = a0 * C1[7] + a1 * C1[6] + a2 * C1[5] + a3 * C1[4] + a4 * C1[3] + a5 * C1[2] +
                  HI(a0, C1[6]) + HI(a1, C1[5]) + HI(a2, C1[4]) + HI(a3, C1[3]) + HI(a4, C1[2]) +
                  HI(a5, C1[1]);
    /* 2^64*(M(x)+M(beta x)+M(beta^2 x))/p is 2^64 or 2^65: the third value lies in (g2, g2 + 2E]. */
    uint64_t g2 = 0 - g0 - g1 - 2 * FRAC_ERR;
    const uint64_t LIM = ~0ULL - 2 * FRAC_ERR - 1;
    if (g0 > LIM || g1 > LIM || g2 > LIM) goto exact;
    if (g0 < g1) {
        if (g0 < g2) {
            if (g1 - g0 > 2 * FRAC_ERR && g2 - g0 > 2 * FRAC_ERR) return 0;
        } else if (g0 - g2 > 2 * FRAC_ERR && g1 - g2 > 2 * FRAC_ERR) return 2;
    } else {
        if (g1 < g2) {
            if (g0 - g1 > 2 * FRAC_ERR && g2 - g1 > 2 * FRAC_ERR) return 1;
        } else if (g0 - g2 > 2 * FRAC_ERR && g1 - g2 > 2 * FRAC_ERR) return 2;
    }
exact:
    frac_fallbacks++;
    return canon_k_exact(x);
}

/* Low limb of the class key and the sign bit s (1: y is negated in the canonical point). */
static inline uint64_t key_lo(const pt *p, int *s) {
    uint64_t y0 = p->y.l[0];
    *s = (int)(y0 & 1);
    return (y0 & 1) ? P[0] - y0 : y0;
}
static void key_full(fe *k, const pt *p) {
    if (p->y.l[0] & 1) fneg(k, &p->y); else *k = p->y;
}

/* Table entry sigma(R_i) for the walk's next step from W. */
static inline const pt *select_adder(const pt *p, uint64_t *klo) {
    int s;
    uint64_t kl = key_lo(p, &s);
    int k = canon_k(&p->x);
    *klo = kl;
    return &T6[(kl >> KEY_SHIFT) & (NA - 1)][k][s];
}

static int pt_eq(const pt *a, const pt *b) { return feq(&a->x, &b->x) && feq(&a->y, &b->y); }

/* Unbatched affine addition/doubling with its own inversion (rare paths). */
static int add_single(pt *r, const pt *a, const pt *b) {
    fe num, den, lam, t, x3, y3;
    unbatched++;
    if (feq(&a->x, &b->x)) {
        if (!feq(&a->y, &b->y) || fis_zero(&a->y)) return 0;
        fsqr(&num, &a->x);
        fadd(&t, &num, &num);
        fadd(&num, &num, &t);
        fadd(&den, &a->y, &a->y);
    } else {
        fsub(&num, &b->y, &a->y);
        fsub(&den, &b->x, &a->x);
    }
    blst_fp_inverse((blst_fp *)&den, (const blst_fp *)&den);
    fmul(&lam, &num, &den);
    fsqr(&x3, &lam);
    fsub(&x3, &x3, &a->x);
    fsub(&x3, &x3, &b->x);
    fsub(&t, &a->x, &x3);
    fmul(&y3, &lam, &t);
    fsub(&y3, &y3, &a->y);
    r->x = x3;
    r->y = y3;
    return 1;
}

/* One walk step from the actual point a, unbatched. */
static int f_single(pt *r, const pt *a) {
    uint64_t kl;
    return add_single(r, a, select_adder(a, &kl));
}

/* Canonical point Wc = (beta^k x, s y), exact. */
static void canon_point(pt *c, const pt *a) {
    fe s;
    int k = canon_k_exact(&a->x);
    c->x = a->x;
    if (k >= 1) fmul(&c->x, &a->x, &BETA);
    if (k == 2) { fadd(&s, &a->x, &c->x); fneg(&c->x, &s); }
    key_full(&c->y, a);
}

/* Deterministic exit from a fruitless cycle: double the canonical point of the
 * cycle's least class key. */
static int escape(pt *r, const pt *cyc, int len) {
    int m = 0;
    fe km, ki;
    key_full(&km, &cyc[0]);
    for (int i = 1; i < len; i++) {
        key_full(&ki, &cyc[i]);
        if (fcmp(&ki, &km) < 0) { m = i; km = ki; }
    }
    pt c;
    canon_point(&c, &cyc[m]);
    return add_single(r, &c, &c);
}

static int same_class(const pt *a, const pt *b) {
    fe ka, kb;
    key_full(&ka, a);
    key_full(&kb, b);
    return feq(&ka, &kb);
}

static int is_dp_rest(const pt *p) {
    /* the low DP bits (key bits 12..63) were zero; check key bits 64..71. */
    fe k;
    key_full(&k, p);
    return (k.l[1] & ((1ULL << (DP_BITS + ADDER_BITS + KEY_SHIFT - 64)) - 1)) == 0;
}

static void batch_step(uint64_t t) {
    fe acc, inv, iw, dy, lam, l2, tt, y3;
    pt *Wc = W, *Wn = WN;
    const int anchor_now = t % ANCHOR_PERIOD == 0;
    fsub(&D[0], &SEL[0]->x, &Wc[0].x);
    PRE[1] = D[0];
    for (int w = 1; w < B - 1; w++) {
        fsub(&D[w], &SEL[w]->x, &Wc[w].x);
        ref_mul(&PRE[w + 1], &PRE[w], &D[w]);
    }
    fsub(&D[B - 1], &SEL[B - 1]->x, &Wc[B - 1].x);
    ref_mul(&acc, &PRE[B - 1], &D[B - 1]);
    if (fis_zero(&acc)) { exceptional++; abort(); }
    blst_fp_inverse((blst_fp *)&inv, (const blst_fp *)&acc);
    for (int w = B - 1; w >= 0; w--) {
        const pt *R = SEL[w];
        pt *nw = &Wn[w];
        if (w > 0) {
            ref_mul(&iw, &inv, &PRE[w]);
            ref_mul(&inv, &inv, &D[w]);
        } else {
            iw = inv;
        }
        fsub(&dy, &R->y, &Wc[w].y);
        ref_mul(&lam, &dy, &iw);
        ref_sqr(&l2, &lam);
        fsub2(&nw->x, &l2, &Wc[w].x, &R->x);
        fsub(&tt, &Wc[w].x, &nw->x);
        ref_mul(&y3, &lam, &tt);
        fsub(&nw->y, &y3, &Wc[w].y);
        uint64_t kl;
        const pt *sel = select_adder(nw, &kl);
        if (__builtin_expect(kl == prevkey[w], 0)) {
            pt g, cyc[2];
            if (f_single(&g, nw) && same_class(&g, &Wc[w])) {
                cyc[0] = Wc[w];
                cyc[1] = *nw;
                if (!escape(nw, cyc, 2)) { exceptional++; abort(); }
                two_cycles++;
                sel = select_adder(nw, &kl);
            }
        } else if (__builtin_expect(kl == anchorkey[w], 0)) {
            pt cyc[MAX_CYCLE + 1], s;
            int len = 1;
            cyc[0] = *nw;
            if (!f_single(&s, nw)) { exceptional++; abort(); }
            while (!same_class(&s, nw) && len <= MAX_CYCLE) {
                cyc[len++] = s;
                if (!f_single(&s, &s)) { exceptional++; abort(); }
            }
            if (same_class(&s, nw)) {
                if (!escape(nw, cyc, len)) { exceptional++; abort(); }
                long_cycles++;
                sel = select_adder(nw, &kl);
            }
        }
        prevkey[w] = curkey[w];
        curkey[w] = kl;
        SEL[w] = sel;
        if (anchor_now) anchorkey[w] = kl;
        if (__builtin_expect((kl >> (ADDER_BITS + KEY_SHIFT)) == 0, 0) && is_dp_rest(nw)) dps++;
    }
    W = Wn;
    WN = Wc;
}

static void read_fe(FILE *f, fe *r) {
    char hex[97];
    unsigned char b[48];
    if (fscanf(f, "%96s", hex) != 1 || strlen(hex) != 96) { fprintf(stderr, "bad input\n"); exit(2); }
    for (int i = 0; i < 48; i++) {
        unsigned v;
        sscanf(hex + 2 * i, "%2x", &v);
        b[i] = (unsigned char)v;
    }
    blst_fp_from_bendian((blst_fp *)r, b);
}

static void print_fe(const fe *a) {
    unsigned char b[48];
    blst_bendian_from_fp(b, (const blst_fp *)a);
    for (int i = 0; i < 48; i++) printf("%02x", b[i]);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: step <inputs> <batch-steps> [dump|single]\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) return 2;
    for (int i = 0; i < NA; i++) { read_fe(f, &T[i].x); read_fe(f, &T[i].y); }
    for (int w = 0; w < B; w++) { read_fe(f, &W[w].x); read_fe(f, &W[w].y); }
    fclose(f);
    for (int i = 0; i < NA; i++) {
        fe bx, b2x, s, ny;
        fmul(&bx, &T[i].x, &BETA);
        fadd(&s, &T[i].x, &bx);
        fneg(&b2x, &s);
        fneg(&ny, &T[i].y);
        /* beta^-k X: k=0 -> X, k=1 -> beta^2 X, k=2 -> beta X */
        const fe *xs[3] = {&T[i].x, &b2x, &bx};
        for (int k = 0; k < 3; k++) {
            T6[i][k][0].x = T6[i][k][1].x = *xs[k];
            T6[i][k][0].y = T[i].y;
            T6[i][k][1].y = ny;
        }
    }
    for (int w = 0; w < B; w++) {
        uint64_t kl;
        SEL[w] = select_adder(&W[w], &kl);
        curkey[w] = kl;
        prevkey[w] = anchorkey[w] = ~0ULL;
    }
    frac_fallbacks = 0;
    int n = atoi(argv[2]);
    if (argc > 3 && strcmp(argv[3], "single") == 0) {
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
        for (int i = 0; i < n; i++) {
            if (!f_single(&W[i % B], &W[i % B])) { exceptional++; abort(); }
        }
        CALLGRIND_TOGGLE_COLLECT;
        fprintf(stderr, "unbatched=%llu\n", (unsigned long long)unbatched);
        return 0;
    }
#ifndef REF_ONLY
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
#endif
    for (int t = 1; t <= n; t++) batch_step((uint64_t)t);
#ifndef REF_ONLY
    CALLGRIND_TOGGLE_COLLECT;
#endif
    fprintf(stderr, "walk_steps=%llu two_cycles=%llu long_cycles=%llu dps=%llu exceptional=%llu unbatched=%llu frac_fallbacks=%llu\n",
            (unsigned long long)n * B, (unsigned long long)two_cycles, (unsigned long long)long_cycles,
            (unsigned long long)dps, (unsigned long long)exceptional, (unsigned long long)unbatched,
            (unsigned long long)frac_fallbacks);
    if (argc > 3) {
        for (int w = 0; w < B; w++) { print_fe(&W[w].x); printf(" "); print_fe(&W[w].y); printf("\n"); }
    }
    return 0;
}
