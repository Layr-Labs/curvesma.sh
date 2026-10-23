/* Batched order-6 class rho step on secp256k1 (sigma walk), every instruction charged.
 *
 * "Instruction route": the step does not use the 6M reference operation. All field
 * arithmetic is secp_field.h (4 x 64-bit limbs, fully reduced, GCC inline AArch64
 * assembly), and every instruction of the step is charged at 1/C.
 *
 * Classes and walk (as in the reference-operation package; see walk_ref.py):
 *  - class key: key(W) = y if y is even, else p - y (even); y determines the class.
 *  - canonical point Wc = (beta^k x, s y), beta^k x least, s y even.
 *  - the walk stores W = sigma(Wc) and steps W <- W + sigma(R_i), i = (key >> 1) mod NA,
 *    with sigma(R_i) read from a table of the six images of R_i.
 *  - k from 64-bit fractions: with C1 = floor(2^384 beta/p) (limbs 1..5),
 *    g1 = sum_{i+j=5} lo(a_i C1_j) + sum_{i+j=4} hi(a_i C1_j) approximates
 *    2^64 (beta x mod p)/p with error in [0,9); 2^64 x/p is approximated by a_3
 *    (p is within 2^33 of 2^256) with error in [0,2); the third value follows from
 *    x + beta x + beta^2 x in {p, 2p}. Ambiguous cases are computed exactly.
 *  - trails: a walk that reaches a distinguished point restarts from a fresh start
 *    point, so each trail reports one distinguished point once; a trail that crosses
 *    more than 2^DP_BITS anchor boundaries is abandoned and restarted the same way.
 *
 * Inputs: 2^ADDER_BITS adders, B start points, SPARES further start points, then 2^ADDER_BITS escape points.
 * Usage: step <inputs.txt> <batch-steps> [dump]
 *        step <inputs.txt> <count> single   (count unbatched f(W) evaluations)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#include "secp_field.h"
#include "secp_sub.h"

#ifndef ADDER_BITS
#define ADDER_BITS 16
#endif
#define NA (1 << ADDER_BITS)
#ifndef B
#define B 65536
#endif
#define ANCHOR_PERIOD 64
#define MAX_CYCLE 64
#ifndef DP_BITS
#define DP_BITS 60
#endif
#define SPARES 1024
/* A trail is abandoned once it crosses more than 2^DP_BITS anchor boundaries without a
 * distinguished point, i.e. after between CAP - 64 and CAP = 64/theta steps. */
#ifndef CAP_PERIODS
#define CAP_PERIODS (1ULL << DP_BITS)
#endif
/* The predecessor and anchor filters compare the low 64 key bits; a test build can narrow
 * them to force false filter matches. Exact checks decide every escape. */
#ifndef FILTER_MASK
#define FILTER_MASK (~0ULL)
#endif
#define P0 0xFFFFFFFEFFFFFC2FULL

typedef struct { fe x, y; } pt;

static const fe ZERO = {{0}};
static const fe BETA = {{0xc1396c28719501eeULL, 0x9cf0497512f58995ULL, 0x6e64479eac3434e9ULL, 0x7ae96a2b657c0710ULL}};
/* Not static: loaded rather than rebuilt from immediates. C1 limbs 1..5. */
uint64_t FRAC_K[5] = {0x7a4a36aebcbb3d53, 0xc1396c28ec7e6dee, 0x9cf0497512f58995,
                      0x6e64479eac3434e9, 0x7ae96a2b657c0710};
#define FRAC_ERR 11

static pt T[NA];
static pt T6[NA][3][2];
static pt ESC[NA];          /* escape table, independent of the walk adders */
static pt WA[B], WB[B];
static pt *W = WA, *WN = WB;
static const pt *SEL[B];
static fe D[B], PRE[B];
static uint64_t curkey[B], prevkey[B], anchorkey[B];
static uint64_t tper[B];     /* anchor boundaries crossed by the current trail of each walk */
static pt SPARE[SPARES];
static uint64_t two_cycles, long_cycles, dps, exceptional, unbatched, frac_fallbacks;
static uint64_t abandoned, next_spare, ahead_dps;

static inline int feq(const fe *a, const fe *b) { return memcmp(a->v, b->v, sizeof a->v) == 0; }
static inline int fis_zero(const fe *a) { return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0; }
static inline int fcmp(const fe *a, const fe *b) {
    for (int i = 3; i >= 0; i--)
        if (a->v[i] != b->v[i]) return a->v[i] < b->v[i] ? -1 : 1;
    return 0;
}

static inline void sqr_n(fe *r, const fe *a, int n) {
    sf_sqr(r, a);
    while (--n) sf_sqr(r, r);
}

/* a^(p-2) with libsecp256k1's addition chain: 255 squarings, 15 multiplications. */
static __attribute__((noinline)) void finv(fe *r, const fe *a) {
    fe x2, x3, x6, x9, x11, x22, x44, x88, x176, x220, x223, t;
    sf_sqr(&x2, a);       sf_mul(&x2, &x2, a);
    sf_sqr(&x3, &x2);     sf_mul(&x3, &x3, a);
    sqr_n(&x6, &x3, 3);   sf_mul(&x6, &x6, &x3);
    sqr_n(&x9, &x6, 3);   sf_mul(&x9, &x9, &x3);
    sqr_n(&x11, &x9, 2);  sf_mul(&x11, &x11, &x2);
    sqr_n(&x22, &x11, 11); sf_mul(&x22, &x22, &x11);
    sqr_n(&x44, &x22, 22); sf_mul(&x44, &x44, &x22);
    sqr_n(&x88, &x44, 44); sf_mul(&x88, &x88, &x44);
    sqr_n(&x176, &x88, 88); sf_mul(&x176, &x176, &x88);
    sqr_n(&x220, &x176, 44); sf_mul(&x220, &x220, &x44);
    sqr_n(&x223, &x220, 3); sf_mul(&x223, &x223, &x3);
    sqr_n(&t, &x223, 23); sf_mul(&t, &t, &x22);
    sqr_n(&t, &t, 5);     sf_mul(&t, &t, a);
    sqr_n(&t, &t, 3);     sf_mul(&t, &t, &x2);
    sqr_n(&t, &t, 2);     sf_mul(r, &t, a);
}

static __attribute__((noinline)) int canon_k_exact(const fe *x) {
    fe bx, s, b2x;
    sf_mul(&bx, x, &BETA);
    sf_add(&s, x, &bx);
    sf_neg(&b2x, &s);
    int k = 0;
    const fe *m = x;
    if (fcmp(&bx, m) < 0) { k = 1; m = &bx; }
    if (fcmp(&b2x, m) < 0) k = 2;
    return k;
}

#define HI(a, b) ((uint64_t)(((unsigned __int128)(a) * (b)) >> 64))

/* Branch-free rotation decision. For values g0, g1, g2 (g2 = -g0-g1-2E mod 2^64) it returns
 * the index of the least value and sets *ok when all three lie below LIM and the least is more
 * than 2E below the other two; the caller falls back to the exact computation otherwise. */
static inline int canon_pick(uint64_t g0, uint64_t g1, uint64_t g2, int *ok) {
    const uint64_t E2 = 2 * FRAC_ERR;
    const uint64_t LIM = ~0ULL - E2 - 1;
    uint64_t m = g0;
    int k = 0;
    int lt = g1 < m;
    m = lt ? g1 : m;
    k = lt ? 1 : k;
    lt = g2 < m;
    m = lt ? g2 : m;
    k = lt ? 2 : k;
    int sep = (g0 - m > E2) + (g1 - m > E2) + (g2 - m > E2);
    *ok = (g0 <= LIM) & (g1 <= LIM) & (g2 <= LIM) & (sep == 2);
    return k;
}

static inline int canon_k(const fe *x) {
    const uint64_t a0 = x->v[0], a1 = x->v[1], a2 = x->v[2], a3 = x->v[3];
    uint64_t g0 = a3;
    uint64_t g1 = a1 * FRAC_K[3] + a2 * FRAC_K[2] + a3 * FRAC_K[1] + a0 * FRAC_K[4] +
                  HI(a0, FRAC_K[3]) + HI(a1, FRAC_K[2]) + HI(a2, FRAC_K[1]) + HI(a3, FRAC_K[0]);
    uint64_t g2 = 0 - g0 - g1 - 2 * FRAC_ERR;
    int ok;
    int k = canon_pick(g0, g1, g2, &ok);
    if (__builtin_expect(!ok, 0)) {
        frac_fallbacks++;
        return canon_k_exact(x);
    }
    return k;
}

static inline uint64_t key_lo(const pt *p, int *s) {
    uint64_t y0 = p->y.v[0];
    *s = (int)(y0 & 1);
    const uint64_t m = 0 - (y0 & 1);
    return (y0 & ~m) | ((P0 - y0) & m);
}
static void key_full(fe *k, const pt *p) {
    if (p->y.v[0] & 1) sf_neg(k, &p->y); else *k = p->y;
}
/* Folds each distinguished point's key into a digest, which the replica checks: it shows
 * which point every trail records, look-ahead points included. */
static uint64_t dp_digest;
static void record_dp(const pt *p) {
    fe k;
    key_full(&k, p);
    for (int i = 0; i < 4; i++) dp_digest = dp_digest * 1099511628211ULL + k.v[i];
}
static inline const pt *select_adder(const pt *p, uint64_t *klo) {
    int s;
    uint64_t kl = key_lo(p, &s);
    *klo = kl;
    return &T6[(kl >> 1) & (NA - 1)][canon_k(&p->x)][s];
}
static int same_class(const pt *a, const pt *b) {
    fe ka, kb;
    key_full(&ka, a);
    key_full(&kb, b);
    return feq(&ka, &kb);
}

static int add_single(pt *r, const pt *a, const pt *b) {
    fe num, den, lam, t, x3, y3;
    unbatched++;
    if (feq(&a->x, &b->x)) {
        if (!feq(&a->y, &b->y) || fis_zero(&a->y)) return 0;
        sf_sqr(&num, &a->x);
        sf_add(&t, &num, &num);
        sf_add(&num, &num, &t);
        sf_add(&den, &a->y, &a->y);
    } else {
        sf_sub(&num, &b->y, &a->y);
        sf_sub(&den, &b->x, &a->x);
    }
    finv(&den, &den);
    sf_mul(&lam, &num, &den);
    sf_sqr(&x3, &lam);
    sf_sub(&x3, &x3, &a->x);
    sf_sub(&x3, &x3, &b->x);
    sf_sub(&t, &a->x, &x3);
    sf_mul(&y3, &lam, &t);
    sf_sub(&y3, &y3, &a->y);
    r->x = x3;
    r->y = y3;
    return 1;
}
static int f_single(pt *r, const pt *a) {
    uint64_t kl;
    return add_single(r, a, select_adder(a, &kl));
}
static void canon_point(pt *c, const pt *a) {
    fe s;
    int k = canon_k_exact(&a->x);
    c->x = a->x;
    if (k >= 1) sf_mul(&c->x, &a->x, &BETA);
    if (k == 2) { sf_add(&s, &a->x, &c->x); sf_neg(&c->x, &s); }
    key_full(&c->y, a);
}
/* Adder index for an escape: key bits ADDER_BITS+1+DP_BITS and up, which the walk index
 * (bits 1..ADDER_BITS) and the distinguished-point test (the DP_BITS bits below) do not use. */
static inline unsigned esc_index(const fe *k) {
    const int shift = ADDER_BITS + 1 + DP_BITS;
    unsigned __int128 w = ((unsigned __int128)k->v[1] << 64) | k->v[0];
    return (unsigned)((w >> shift) & (NA - 1));
}
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
    /* The escape adds a point of the separate table ESC, chosen by the least member's key from
     * the bits above the distinguished-point test. ESC is independent of the walk adders, so the
     * escape point is fresh given the cycle and cannot lie in the cycle. */
    return add_single(r, &c, &ESC[esc_index(&km)]);
}
#if DP_BITS + ADDER_BITS + 1 <= 64
#define IS_DP(kl, p) ((((kl) >> (ADDER_BITS + 1)) & ((1ULL << DP_BITS) - 1)) == 0)
#else
static int is_dp_rest(const pt *p) {
    fe k;
    key_full(&k, p);
    return (k.v[1] & ((1ULL << (DP_BITS + ADDER_BITS + 1 - 64)) - 1)) == 0;
}
#define IS_DP(kl, p) (((kl) >> (ADDER_BITS + 1)) == 0 && is_dp_rest(p))
#endif

/* Distinguished-point test for a look-ahead point. */
static int point_is_dp(const pt *p) {
    uint64_t k = 0;
    int sg;
    k = key_lo(p, &sg);
    (void)sg;
    return IS_DP(k, p);
}

/* Start a new trail in walk w. In production the finished trail's distinguished key and
 * coins go to the table and fresh coins give the new start point (charged separately);
 * this harness takes the start points from the spare points of the inputs. */
static __attribute__((noinline)) void new_trail(int w, pt *nw, uint64_t *kl, const pt **sel) {
    *nw = SPARE[next_spare++ % SPARES];
    *sel = select_adder(nw, kl);
    prevkey[w] = anchorkey[w] = ~0ULL;
    tper[w] = 0;
}

/* A distinguished point ends its trail: record it, then start a new trail. */
static __attribute__((noinline)) void dp_restart(int w, pt *nw, uint64_t *kl, const pt **sel) {
    record_dp(nw);
    new_trail(w, nw, kl, sel);
}

/* Rare paths, out of line: an exact cycle check and its escape. The common step only
 * compares the filters and calls these when a filter matches. */
static __attribute__((noinline)) void rare_2cycle(int w, pt *nw, const pt *Wc, uint64_t *kl, const pt **sel,
                                                  int *dp, int *esc) {
    pt g, cyc[2];
    if (f_single(&g, nw) && same_class(&g, &Wc[w])) {
        cyc[0] = Wc[w];
        cyc[1] = *nw;
        if (!escape(nw, cyc, 2)) { exceptional++; abort(); }
        two_cycles++;
        *esc = 1;
        *sel = select_adder(nw, kl);
        *dp = IS_DP(*kl, nw);
    }
}

static __attribute__((noinline)) void rare_anchor(int w, pt *nw, uint64_t *kl, const pt **sel, int *dp) {
    pt cyc[MAX_CYCLE + 1], s;
    int len = 1, sg;
    cyc[0] = *nw;
    if (!f_single(&s, nw)) { exceptional++; abort(); }
    while (!same_class(&s, nw) && len <= MAX_CYCLE) {
        if (point_is_dp(&s)) {
            /* The walk reaches s next, so s is the trail's distinguished point and its record. */
            *nw = s;
            *kl = key_lo(nw, &sg);
            *dp = 1;
            ahead_dps++;
            break;
        }
        cyc[len++] = s;
        if (!f_single(&s, &s)) { exceptional++; abort(); }
    }
    if (!*dp && same_class(&s, nw)) {
        if (!escape(nw, cyc, len)) { exceptional++; abort(); }
        long_cycles++;
        *sel = select_adder(nw, kl);
        *dp = IS_DP(*kl, nw);
    }
}

static __attribute__((noinline)) void batch_step(uint64_t t) {
    fe acc, inv, iw, dy, lam, l2, tt, y3;
    pt *Wc = W, *Wn = WN;
    const int anchor_now = t % ANCHOR_PERIOD == 0;
    s_subb_bl(&D[0], &SEL[0]->x, &Wc[0].x);
    PRE[1] = D[0];
    for (int w = 1; w < B - 1; w++) {
        s_subb_bl(&D[w], &SEL[w]->x, &Wc[w].x);
        sf_mul(&PRE[w + 1], &PRE[w], &D[w]);
    }
    s_subb_bl(&D[B - 1], &SEL[B - 1]->x, &Wc[B - 1].x);
    sf_mul(&acc, &PRE[B - 1], &D[B - 1]);
    if (fis_zero(&acc)) { exceptional++; abort(); }
    finv(&inv, &acc);
    for (int w = B - 1; w >= 0; w--) {
        const pt *R = SEL[w];
        pt *nw = &Wn[w];
        if (w > 0) {
            sf_mul(&iw, &inv, &PRE[w]);
            sf_mul(&inv, &inv, &D[w]);
        } else {
            iw = inv;
        }
        s_subb_bl(&dy, &R->y, &Wc[w].y);
        sf_mul(&lam, &dy, &iw);
        sf_sqr(&l2, &lam);
        s_sub2_bl(&nw->x, &l2, &Wc[w].x, &R->x);
        s_subb_bl(&tt, &Wc[w].x, &nw->x);
        sf_mul(&y3, &lam, &tt);
        s_subb_bl(&nw->y, &y3, &Wc[w].y);
        uint64_t kl;
        const pt *sel = select_adder(nw, &kl);
        /* The distinguished-point test comes first, on every point the walk produces
         * (the new point, look-ahead points, escape points), so an escape never skips one;
         * the filters only decide when an exact cycle check runs. */
        int dp = IS_DP(kl, nw), esc = 0;
        if (!dp && __builtin_expect((kl & FILTER_MASK) == prevkey[w], 0))
            rare_2cycle(w, nw, Wc, &kl, &sel, &dp, &esc);
        /* Not "else": a false predecessor match must not suppress the anchor check, or a walk
         * whose anchor shares filter bits with the class two steps back could never escape. */
        if (!dp && !esc && __builtin_expect((kl & FILTER_MASK) == anchorkey[w], 0))
            rare_anchor(w, nw, &kl, &sel, &dp);
        if (__builtin_expect(dp, 0)) {
            dps++;
            dp_restart(w, nw, &kl, &sel);
        } else if (anchor_now) {
            anchorkey[w] = kl & FILTER_MASK;
            if (__builtin_expect(++tper[w] > CAP_PERIODS, 0)) {
                abandoned++;
                new_trail(w, nw, &kl, &sel);
            } else {
                prevkey[w] = curkey[w];
            }
        } else {
            prevkey[w] = curkey[w];
        }
        curkey[w] = kl & FILTER_MASK;
        SEL[w] = sel;
    }
    W = Wn;
    WN = Wc;
}

static void read_fe(FILE *f, fe *r) {
    char hex[65];
    if (fscanf(f, "%64s", hex) != 1 || strlen(hex) != 64) { fprintf(stderr, "bad input\n"); exit(2); }
    for (int i = 0; i < 4; i++) {
        char limb[17];
        memcpy(limb, hex + 16 * (3 - i), 16);
        limb[16] = 0;
        r->v[i] = strtoull(limb, NULL, 16);
    }
}
static void print_fe(const fe *a) {
    for (int i = 3; i >= 0; i--) printf("%016llx", (unsigned long long)a->v[i]);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: step <inputs> <batch-steps> [dump|single]\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) return 2;
    for (int i = 0; i < NA; i++) { read_fe(f, &T[i].x); read_fe(f, &T[i].y); }
    for (int w = 0; w < B; w++) { read_fe(f, &W[w].x); read_fe(f, &W[w].y); }
    for (int i = 0; i < SPARES; i++) { read_fe(f, &SPARE[i].x); read_fe(f, &SPARE[i].y); }
    for (int i = 0; i < NA; i++) { read_fe(f, &ESC[i].x); read_fe(f, &ESC[i].y); }
    fclose(f);
    for (int i = 0; i < NA; i++) {
        fe bx, b2x, s, ny;
        sf_mul(&bx, &T[i].x, &BETA);
        sf_add(&s, &T[i].x, &bx);
        sf_neg(&b2x, &s);
        sf_neg(&ny, &T[i].y);
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
        curkey[w] = kl & FILTER_MASK;
        prevkey[w] = anchorkey[w] = ~0ULL;
    }
    frac_fallbacks = 0;
    int n = atoi(argv[2]);
    if (argc > 3 && strcmp(argv[3], "single") == 0) {
        CALLGRIND_ZERO_STATS;
        CALLGRIND_TOGGLE_COLLECT;
        for (int i = 0; i < n; i++)
            if (!f_single(&W[i % B], &W[i % B])) { exceptional++; abort(); }
        CALLGRIND_TOGGLE_COLLECT;
        fprintf(stderr, "unbatched=%llu\n", (unsigned long long)unbatched);
        return 0;
    }
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for (int t = 1; t <= n; t++) batch_step((uint64_t)t);
    CALLGRIND_TOGGLE_COLLECT;
    fprintf(stderr, "walk_steps=%llu two_cycles=%llu long_cycles=%llu dps=%llu ahead_dps=%llu abandoned=%llu exceptional=%llu unbatched=%llu frac_fallbacks=%llu dp_digest=%016llx\n",
            (unsigned long long)n * B, (unsigned long long)two_cycles, (unsigned long long)long_cycles,
            (unsigned long long)dps, (unsigned long long)ahead_dps, (unsigned long long)abandoned, (unsigned long long)exceptional,
            (unsigned long long)unbatched, (unsigned long long)frac_fallbacks, (unsigned long long)dp_digest);
    if (argc > 3)
        for (int w = 0; w < B; w++) { print_fe(&W[w].x); printf(" "); print_fe(&W[w].y); printf("\n"); }
    return 0;
}
