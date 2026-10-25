/* Batched negation-class rho step on P-256 (sigma walk) for Callgrind counts; every instruction is
 * charged. Field: p256_field.h, Montgomery M(v) = v*2^256 mod p, four 64-bit limbs, GCC inline
 * AArch64 asm; batch_step uses p256_fused.h. Walk (walk_ref.py replays it): class {W,-W}, key(W) =
 * M(x); canonical Wc = (x, s y), s = +-1 making M(s y) even; W <- W + sigma(R_i), i = key mod
 * 2^ADDER_BITS, from a table of R_i and -R_i. A step is one batched affine addition: 3M for the
 * simultaneous inversion, lambda (M), lambda^2 (S), y3 (M). A trail at a distinguished point (a
 * look-ahead one too) records it once in dp_digest and restarts from a fresh start point, as does
 * one crossing over 2^DP_BITS anchor boundaries (abandoned). loops counts repeats among the last
 * RING_K escape keys. Usage: step <inputs> <batch-steps> [dump]; step <inputs> <count> single.
 * Inputs (walk_ref.py gen): 2^ADDER_BITS adders, B start points, SPARES spare starts. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <valgrind/callgrind.h>
#include "p256_field.h"
#include "p256_sub.h"
#include "p256_fused.h"

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
 * them to force false matches. Exact checks decide every escape. */
#ifndef FILTER_MASK
#define FILTER_MASK (~0ULL)
#endif
#if DP_BITS + ADDER_BITS <= 64
#define IS_DP(kl, p) ((((kl) >> ADDER_BITS) & ((1ULL << DP_BITS) - 1)) == 0)
#define IS_DP1 IS_DP
#else
#define IS_DP(kl, p) (((kl) >> ADDER_BITS) == 0 && \
                      ((p)->x.v[1] & ((1ULL << (DP_BITS + ADDER_BITS - 64)) - 1)) == 0)
/* IS_DP as one compare (batch_step): both limbs are read at every step. */
#define IS_DP1(kl, p) ((((kl) >> ADDER_BITS) | \
                        ((p)->x.v[1] & ((1ULL << (DP_BITS + ADDER_BITS - 64)) - 1))) == 0)
#endif

typedef struct { fe x, y; } pt; /* affine, Montgomery form */

static const fe P = {{0xffffffffffffffffULL, 0x00000000ffffffffULL, 0, 0xffffffff00000001ULL}};
static const fe ZERO = {{0}};
/* p - 2, for Fermat inversion */
static const uint64_t PM2[4] = {0xfffffffffffffffdULL, 0x00000000ffffffffULL, 0, 0xffffffff00000001ULL};

static pt T[NA]; /* adders R_i as read */
static pt T2[NA][2]; /* T2[i][s] = (X_i, (-1)^s Y_i) */
static pt WA[B], WB[B];
static pt *W = WA, *WN = WB; /* current and next states, swapped every batch step */
static const pt *SEL[B]; /* sigma(R_i) for each walk's next step */
static fe D[B], PRE[B]; /* PRE[w] = D[0]*...*D[w-1] mod p, below 2^256 */
static uint64_t curkey[B], prevkey[B], anchorkey[B]; /* low key limbs: W, its predecessor, anchor */
static uint64_t tper[B]; /* anchor boundaries crossed by the current trail of each walk */
static pt SPARE[SPARES];
#define RING_K 4
static uint64_t ring[B][RING_K];
static int ring_n[B], ring_pos[B];
static uint64_t loops;
static uint64_t two_cycles, long_cycles, dps, exceptional, unbatched;
static uint64_t abandoned, next_spare, ahead_dps;

static inline int feq(const fe *a, const fe *b) { return memcmp(a->v, b->v, sizeof a->v) == 0; }
static inline int fis_zero(const fe *a) { return (a->v[0] | a->v[1] | a->v[2] | a->v[3]) == 0; }
static inline int fcmp(const fe *a, const fe *b) {
    for (int i = 3; i >= 0; i--) {
        if (a->v[i] != b->v[i]) return a->v[i] < b->v[i] ? -1 : 1;
    }
    return 0;
}

/* r = a^(p-2) = a^-1 (Fermat), left-to-right square and multiply. */
static __attribute__((noinline)) void finv(fe *r, const fe *a) {
    fe acc = *a;
    for (int i = 254; i >= 0; i--) {
        k_sqr(&acc, &acc);
        if ((PM2[i >> 6] >> (i & 63)) & 1) k_mul(&acc, &acc, a);
    }
    *r = acc;
}

static inline uint64_t key_lo(const pt *p) { return p->x.v[0]; }

/* Table entry sigma(R_i) for the walk's next step from W. */
static inline const pt *select_adder(const pt *p) {
    return &T2[p->x.v[0] & (NA - 1)][p->y.v[0] & 1];
}

/* select_adder for batch_step; the empty asm keeps the mask its own instruction at every A. */
static inline const pt *select_adder1(const pt *p) {
    uint64_t i = p->x.v[0] & (NA - 1);
    __asm__("" : "+r"(i));
    return &T2[i][p->y.v[0] & 1];
}

/* Distinguished-point test for a look-ahead point. */
static int point_is_dp(const pt *p) {
    uint64_t k = key_lo(p);
    return IS_DP(k, p);
}

/* Folds each distinguished point's class key (the four limbs of M(x), least significant
 * first) into a digest the replica checks: it shows which point every trail records. */
static uint64_t dp_digest;
static void record_dp(const pt *p) {
    for (int i = 0; i < 4; i++) dp_digest = dp_digest * 1099511628211ULL + p->x.v[i];
}

/* Start a new trail in walk w. In production the finished trail's key and coins go to the
 * table and fresh coins give the new start point (charged separately); this harness takes
 * the spare points of the inputs. */
static __attribute__((noinline)) void new_trail(int w, pt *nw, uint64_t *kl, const pt **sel) {
    *nw = SPARE[next_spare++ % SPARES];
    *kl = key_lo(nw);
    *sel = select_adder(nw);
    prevkey[w] = anchorkey[w] = ~0ULL;
    tper[w] = 0;
    ring_n[w] = ring_pos[w] = 0;
}

/* Diagnostic (loops); no control flow depends on it. */
static inline void esc_seen(int w, uint64_t kl) {
    int n = ring_n[w] < RING_K ? ring_n[w] : RING_K;
    for (int i = 0; i < n; i++)
        if (ring[w][i] == kl) { loops++; return; }
    ring[w][ring_pos[w]] = kl;
    ring_pos[w] = (ring_pos[w] + 1) % RING_K;
    ring_n[w]++;
}

/* Record the distinguished point nw and start a new trail; out of line. */
static __attribute__((noinline)) void dp_restart(int w, pt *nw, uint64_t *kl, const pt **sel) {
    record_dp(nw);
    new_trail(w, nw, kl, sel);
}

static int pt_eq(const pt *a, const pt *b) { return feq(&a->x, &b->x) && feq(&a->y, &b->y); }

/* Unbatched affine addition/doubling with its own inversion (rare paths). */
static int add_single(pt *r, const pt *a, const pt *b) {
    fe num, den, lam, t, x3, y3;
    unbatched++;
    if (feq(&a->x, &b->x)) {
        if (!feq(&a->y, &b->y) || fis_zero(&a->y)) return 0;
        fe three = {{0}};
        k_sqr(&num, &a->x);
        k_add(&t, &num, &num);
        k_add(&num, &num, &t);
        /* a = -3: num = 3x^2 - 3 */
        three.v[0] = 3;
        k_to(&three, &three);
        k_sub(&num, &num, &three);
        k_add(&den, &a->y, &a->y);
    } else {
        k_sub(&num, &b->y, &a->y);
        k_sub(&den, &b->x, &a->x);
    }
    finv(&den, &den);
    k_mul(&lam, &num, &den);
    k_sqr(&x3, &lam);
    k_sub(&x3, &x3, &a->x);
    k_sub(&x3, &x3, &b->x);
    k_sub(&t, &a->x, &x3);
    k_mul(&y3, &lam, &t);
    k_sub(&y3, &y3, &a->y);
    r->x = x3;
    r->y = y3;
    return 1;
}

/* One walk step from the actual point a, unbatched. */
static int f_single(pt *r, const pt *a) { return add_single(r, a, select_adder(a)); }

/* Deterministic exit from a fruitless cycle: double the canonical point of the cycle's
 * least class key. */
static int escape(pt *r, const pt *cyc, int len) {
    int m = 0;
    for (int i = 1; i < len; i++)
        if (fcmp(&cyc[i].x, &cyc[m].x) < 0) m = i;
    pt c = cyc[m];
    if (c.y.v[0] & 1) k_sub(&c.y, &ZERO, &c.y);
    return add_single(r, &c, &c);
}

static inline int same_class(const pt *a, const pt *b) { return feq(&a->x, &b->x); }

/* Rare paths, out of line: an exact cycle check and its escape, called when a filter matches. */
static __attribute__((noinline)) void rare_2cycle(int w, pt *nw, const pt *Wc, uint64_t *kl, const pt **sel,
                                                  int *dp, int *esc) {
    pt g, cyc[2];
    if (f_single(&g, nw) && same_class(&g, &Wc[w])) {
        cyc[0] = Wc[w];
        cyc[1] = *nw;
        if (!escape(nw, cyc, 2)) { exceptional++; abort(); }
        two_cycles++;
        *esc = 1;
        *kl = key_lo(nw);
        *sel = select_adder(nw);
        *dp = IS_DP(*kl, nw);
        if (!*dp) esc_seen(w, *kl);
    }
}

static __attribute__((noinline)) void rare_anchor(int w, pt *nw, uint64_t *kl, const pt **sel, int *dp) {
    pt cyc[MAX_CYCLE + 1], s;
    int len = 1;
    cyc[0] = *nw;
    if (!f_single(&s, nw)) { exceptional++; abort(); }
    while (!same_class(&s, nw) && len <= MAX_CYCLE) {
        if (point_is_dp(&s)) {
            /* The walk reaches s next, so s is the trail's distinguished point and its record. */
            *nw = s;
            *kl = key_lo(nw);
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
        *kl = key_lo(nw);
        *sel = select_adder(nw);
        *dp = IS_DP(*kl, nw);
        if (!*dp) esc_seen(w, *kl);
    }
}

/* A step's tail when its point is distinguished, matches a filter, or would cross the cap; out of line. */
static __attribute__((noinline)) void step_rare(int w, pt *nw, const pt *Wc, int anchor_now) {
    uint64_t kl = key_lo(nw);
    const pt *sel = select_adder(nw);
    /* The distinguished-point test comes first, on every point the walk produces
     * (new, look-ahead and escape points), so an escape never skips one. */
    int dp = IS_DP1(kl, nw), esc = 0;
    if (!dp && (kl & FILTER_MASK) == prevkey[w])
        rare_2cycle(w, nw, Wc, &kl, &sel, &dp, &esc);
    /* Not "else": a false predecessor match must not suppress the anchor check, or a walk
     * whose anchor matches the class two steps back could never escape. */
    if (!dp && !esc && (kl & FILTER_MASK) == anchorkey[w])
        rare_anchor(w, nw, &kl, &sel, &dp);
    if (dp) {
        dps++;
        dp_restart(w, nw, &kl, &sel);
    } else if (anchor_now) {
        anchorkey[w] = kl & FILTER_MASK;
        if (++tper[w] > CAP_PERIODS) {
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

static __attribute__((noinline)) void batch_step(uint64_t t) {
    fe acc, S[3]; /* inv, iw, lam */
    pt *Wc = W, *Wn = WN;
    const int anchor_now = t % ANCHOR_PERIOD == 0;
    k_subb_bl(&D[0], &SEL[0]->x, &Wc[0].x);
    PRE[1] = D[0];
    uintptr_t d = (uintptr_t)&D[1], r = (uintptr_t)&PRE[2], c = (uintptr_t)&Wc[1];
    for (int w = 1; w < B - 1; w++)
        k_submul(&d, &r, &SEL[w]->x, &c);
    k_subb_bl(&D[B - 1], &SEL[B - 1]->x, &Wc[B - 1].x);
    k_subm(&acc, &SEL[B - 1]->x, &Wc[B - 1].x, &PRE[B - 1]);
    if (fis_zero(&acc)) { exceptional++; abort(); }
    finv(&S[0], &acc);
    uintptr_t m = (uintptr_t)&PRE[B], n = (uintptr_t)&D[B];
    for (int w = B - 1; w >= 0; w--) {
        const pt *R = SEL[w];
        pt *nw = &Wn[w];
        if (w > 0)
            k_mul2(S, &m, &n); /* iw = inv*PRE[w], inv = inv*D[w]: below 2^256, not canonical */
        else
            S[1] = S[0];
        k_lxy(S, &R->x, &Wc[w].x, &nw->x);
        uint64_t kl = key_lo(nw), f = kl & FILTER_MASK;
        const pt *sel = select_adder1(nw);
        /* step_rare runs the whole tail; on any other step none of its rare branches applies */
        if (__builtin_expect(IS_DP1(kl, nw) | (f == prevkey[w]) | (f == anchorkey[w]) |
                             (anchor_now && tper[w] >= CAP_PERIODS), 0)) {
            step_rare(w, nw, Wc, anchor_now);
            continue;
        }
        if (anchor_now) {
            anchorkey[w] = f;
            tper[w]++;
        }
        prevkey[w] = curkey[w];
        curkey[w] = f;
        SEL[w] = sel;
    }
    W = Wn;
    WN = Wc;
}

static void read_fe(FILE *f, fe *r) {
    char hex[65];
    if (fscanf(f, "%64s", hex) != 1 || strlen(hex) != 64) { fprintf(stderr, "bad input\n"); exit(2); }
    fe a;
    for (int i = 0; i < 4; i++) {
        char limb[17];
        memcpy(limb, hex + 16 * (3 - i), 16);
        limb[16] = 0;
        a.v[i] = strtoull(limb, NULL, 16);
    }
    k_to(r, &a);
}

static void print_fe(const fe *a) {
    fe t;
    k_from(&t, a);
    for (int i = 3; i >= 0; i--) printf("%016llx", (unsigned long long)t.v[i]);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: step <inputs> <batch-steps> [dump|single]\n"); return 2; }
    FILE *f = fopen(argv[1], "r");
    if (!f) return 2;
    for (int i = 0; i < NA; i++) { read_fe(f, &T[i].x); read_fe(f, &T[i].y); }
    for (int w = 0; w < B; w++) { read_fe(f, &W[w].x); read_fe(f, &W[w].y); }
    for (int i = 0; i < SPARES; i++) { read_fe(f, &SPARE[i].x); read_fe(f, &SPARE[i].y); }
    fclose(f);
    for (int i = 0; i < NA; i++) {
        T2[i][0] = T2[i][1] = T[i];
        k_sub(&T2[i][1].y, &ZERO, &T[i].y);
    }
    for (int w = 0; w < B; w++) {
        SEL[w] = select_adder(&W[w]);
        curkey[w] = key_lo(&W[w]) & FILTER_MASK;
        prevkey[w] = anchorkey[w] = ~0ULL;
    }
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
    CALLGRIND_ZERO_STATS;
    CALLGRIND_TOGGLE_COLLECT;
    for (int t = 1; t <= n; t++) batch_step((uint64_t)t);
    CALLGRIND_TOGGLE_COLLECT;
    fprintf(stderr, "walk_steps=%llu two_cycles=%llu long_cycles=%llu dps=%llu ahead_dps=%llu abandoned=%llu exceptional=%llu unbatched=%llu dp_digest=%016llx loops=%llu\n",
            (unsigned long long)n * B, (unsigned long long)two_cycles, (unsigned long long)long_cycles,
            (unsigned long long)dps, (unsigned long long)ahead_dps, (unsigned long long)abandoned, (unsigned long long)exceptional,
            (unsigned long long)unbatched, (unsigned long long)dp_digest, (unsigned long long)loops);
    if (argc > 3) {
        for (int w = 0; w < B; w++) { print_fe(&W[w].x); printf(" "); print_fe(&W[w].y); printf("\n"); }
    }
    return 0;
}
