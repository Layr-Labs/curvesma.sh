/* ARCH-D v3 hit-cost microbenchmark.
 *
 * v3 of submissions/bls12-381-fixed-g2-invert-g1/stnfs-archd-v2/archd_hit2.c.
 * The enumeration (Schnorr-Euchner DFS over a 12-dim lattice coset inside a
 * 12-dim ellipsoid, ambient state carried one level at a time) is UNCHANGED.
 *
 * v3 repair (review objection 1): v2 packed with radix 1497 (A = 2^9.55 = 748)
 * while the operating region has A = 2^9.95 and, because the fitted ellipsoid is
 * not clipped to the reference box, coordinates up to 1.66*A.  Radix 1497
 * aliased, e.g. (749,0,0,0,0,1) and (-748,1,0,0,0,1).  The production slab id
 * therefore uses balanced radix 4097 (digits in [-2048,2048]), which covers the
 * measured production bound |v_i| <= 1700 at A = 2^10.00 (evidence.py part [2c]).
 *
 *   - the slab id is the FULL signed-coordinate id: 9 of the 12 ambient
 *     coordinates packed in balanced radix B = 4097, held in a pair of uint64
 *     words (a-part, 5 digits; b-part, 4 digits).  9 * log2(4097) = 108.0
 *     significant bits.  No bit of the id is masked before hashing or
 *     comparison.  The packing is injective whenever every digit is in
 *     [-(B-1)/2, (B-1)/2]; packcheck verifies that bound and tests for
 *     collisions, including the rejected v2 counterexample.
 *
 * The hot hit path does not multiply by the radix constant: the packed words
 * are carried incrementally by the walk (setup_dlt), and a hit only negates
 * them for canonicalisation.  Changing 1497 -> 4097 therefore leaves the
 * measured per-hit instruction count unchanged; it changes the setup constants
 * and widens the 128-bit key comparison (still one 64-bit compare per word).
 *   - the point is canonicalised BEFORE the id is computed: a-coordinates are
 *     negated if a5 < 0, b-coordinates if b5 < 0 (the X(-a) = C - X(a) flip is
 *     exact in the balanced-radix packing: flipping every digit negates the
 *     packed word).
 *   - every node owns a descriptor hash table keyed by the full 128-bit id
 *     pair, open addressing with linear probing, exact key comparison on both
 *     words.  The table is prepopulated from a pre-pass over the very same
 *     point set, so every measured probe HITS.  A probe miss in the correctness
 *     or measured path is fatal.
 *   - routing: the destination node index in [0, 2^20) is taken from a
 *     multiplicative hash of BOTH id words; a route table entry maps the node
 *     to one of 8 destination queues.  The update is appended to the
 *     per-destination ring buffer; the owner drains it and applies: probe the
 *     descriptor table with the full key, recompute the 64-bit cell from the
 *     in-slab coordinates, add logp.
 *
 * Coordinate split (9 slab coords + 3 in-slab coords):
 *   a-word digits : v0, v1, v2, v3, v5     (5 digits, 60.0 bits)
 *   b-word digits : v6, v7, v9, v11        (4 digits, 48.0 bits)
 *   in-slab       : v4, v8, v10
 *   sign keys     : a5 = v5, b5 = v11
 *
 * In production the three in-slab coordinates use EXACT strides derived from
 * the slab's 3-dimensional bounding box: cell = off + (v4-lo4) + wa*(v8-lo8)
 * + wab*(v10-lo10), wa = hi4-lo4+1, wab = wa*(hi8-lo8+1).  Together with the
 * injective 9-digit id this makes the full 12-coordinate address injective.
 * The benchmark's desc_payload uses synthetic fixed-width strides so that the
 * 2^25-byte bucket array stays in cache; the addressing instruction sequence
 * is identical (one multiply-add per stride).  packcheck tests the exact-stride
 * cell formula separately.
 *
 * Modes:
 *   ./archd_hit2 stcheck            carried state vs direct v at every hit
 *   ./archd_hit2 refcheck [reps]    reference (recompute from ambient) vs
 *                                   single-memory vs distributed: bit-identical
 *   ./archd_hit2 stats              table load, probe depths, destination use
 *   ./archd_hit2 packcheck          v3: id injectivity on the production range,
 *                                   the rejected v2 alias, exact in-slab strides
 *   ./archd_hit2 icount reps hmode [apply]
 *       20 walk only, state kept live by an asm barrier (no hit work)
 *       21 walk only, v1's dead xor sink (state removed by clang -O2)
 *       23 full hit, single shared memory (descriptor probe + byte update)
 *       24 full hit, send to per-destination buffer (+ apply if arg 4 = 1)
 *       26 all descriptor keys made absent (unsuccessful-probe variant)
 *       27 miss-tolerant probe with the true key (isolates the miss cost)
 *   ./archd_hit2 count              hit count only
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "archd_lattice_d12.h"

#define D LAT_D

/* ------------------------------------------------------------- geometry --- */
/* Production region at A = 2^10.00: max |coordinate| over the unclipped
 * ellipsoid is 1700 (evidence.py part [2c] computes 1699.1 from the embedded
 * covariance).  Digits must therefore cover [-1700,1700]; balanced radix 4097
 * gives [-2048,2048], i.e. ~20% headroom.  AMAX below is the DIGIT bound used
 * by the range check, not the region bound. */
#define AMAX 2048
#define RADIX 4097LL
#define IDTAG (1ULL << 63)          /* occupancy tag in the high bit of k1 */

/* p_j positions: a-word digits v0,v1,v2,v3,v5 ; b-word digits v6,v7,v9,v11 */
typedef struct {
    uint64_t k0, k1;                /* k1 |= IDTAG when occupied */
    uint64_t off, wa, wab;          /* slab descriptor payload */
} ent_t;

#define DLOG 22
#define DN ((size_t)1 << DLOG)
static ent_t *dtab;
static uint64_t dmask;

#define NODELOG 20
#define NNODES ((uint32_t)1 << NODELOG)
#define NDEST 8
static uint8_t *route;              /* node -> destination queue */

#define BLOG 25
#define BN_BYTES ((size_t)1 << BLOG)
static unsigned char *bucket8;
static uint64_t b8mask;
static unsigned short logp = 137;

#define BUFCAP ((uint32_t)1 << 19)
typedef struct { uint64_t k0, k1; int64_t t0, t1, t2; uint16_t lp; } msg_t;
typedef struct { msg_t *m; uint32_t n; } queue_t;
static queue_t msgbuf[NDEST];

static double MU[D][D], BN[D], CST[D];
static int64_t DLT[D][7];
static long long g_hits;

static uint64_t now_ns(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (uint64_t)t.tv_sec * 1000000000ull + (uint64_t)t.tv_nsec;
}

/* ------------------------------------------------------------------ setup */
static void setup(void) {
    double G[D][D], BS[D][D], w[D];
    int i, j, a, b;
    for (i = 0; i < D; i++)
        for (j = 0; j < D; j++) {
            double s = 0.0;
            for (a = 0; a < D; a++) {
                double t = 0.0;
                for (b = 0; b < D; b++) t += LAT_Q[a][b] * (double)LAT_BRED[j][b];
                s += (double)LAT_BRED[i][a] * t;
            }
            G[i][j] = s;
        }
    for (i = 0; i < D; i++) {
        for (j = 0; j < i; j++) {
            double s = G[i][j];
            for (a = 0; a < j; a++) s -= MU[i][a] * MU[j][a] * BN[a];
            MU[i][j] = s / BN[j];
        }
        double s = G[i][i];
        for (a = 0; a < i; a++) s -= MU[i][a] * MU[i][a] * BN[a];
        BN[i] = s;
    }
    for (i = 0; i < D; i++) {
        for (a = 0; a < D; a++) {
            double s = (double)LAT_BRED[i][a];
            for (j = 0; j < i; j++) s -= MU[i][j] * BS[j][a];
            BS[i][a] = s;
        }
    }
    for (a = 0; a < D; a++) {
        double s = 0.0;
        for (b = 0; b < D; b++) s += LAT_Q[a][b] * LAT_CENTER[b];
        w[a] = s;
    }
    for (i = 0; i < D; i++) {
        double s = 0.0;
        for (a = 0; a < D; a++) s += w[a] * BS[i][a];
        CST[i] = s / BN[i];
    }
    for (i = 0; i < D; i++)
        if (!(BN[i] > 0)) { fprintf(stderr, "bad GSO bn[%d]=%g\n", i, BN[i]); exit(1); }
}

static void setup_dlt(void) {
    for (int k = 0; k < D; k++) {
        const int *b = &LAT_BRED[k][0];
        DLT[k][0] = (int64_t)b[0] + RADIX * ((int64_t)b[1] + RADIX * ((int64_t)b[2]
                  + RADIX * ((int64_t)b[3] + RADIX * (int64_t)b[5])));
        DLT[k][1] = (int64_t)b[6] + RADIX * ((int64_t)b[7] + RADIX * ((int64_t)b[9]
                  + RADIX * (int64_t)b[11]));
        DLT[k][2] = b[5];           /* a5 */
        DLT[k][3] = b[11];          /* b5 */
        DLT[k][4] = b[4];           /* in-slab a4  */
        DLT[k][5] = b[8];           /* in-slab b2  */
        DLT[k][6] = b[10];          /* in-slab b4  */
    }
}

static inline uint64_t mixid(uint64_t a, uint64_t b) {
    uint64_t h = a * 0x9E3779B97F4A7C15ULL;
    h ^= (b + 0x165667B19E3779F9ULL) * 0xC2B2AE3D27D4EB4FULL;
    h ^= h >> 31;
    return h;
}

static void fatal_miss(uint64_t k0, uint64_t k1) __attribute__((noreturn));
static void fatal_miss(uint64_t k0, uint64_t k1) {
    fprintf(stderr, "FATAL: slab descriptor probe miss for id (%016llx,%016llx)\n",
            (unsigned long long)k0, (unsigned long long)k1);
    exit(1);
}

/* Exact full-width descriptor probe.  Miss (empty slot) is fatal. */
static inline ent_t *dlookup(uint64_t k0, uint64_t k1) {
    uint64_t i = mixid(k0, k1) & dmask;
    for (;;) {
        ent_t *e = &dtab[i];
        if (e->k1 == (k1 | IDTAG)) { if (e->k0 == k0) return e; }
        else if (e->k1 == 0) fatal_miss(k0, k1);
        i = (i + 1) & dmask;
    }
}

/* Same probe, tolerating absence: returns a fixed fallback descriptor when the
 * key is absent.  Used only by the miss-cost variant and the statistics pass;
 * it never runs in the correctness or reported single/distributed paths. */
static ent_t miss_entry = {0, 0, 0, 1, 1};
static inline ent_t *dlookup_opt(uint64_t k0, uint64_t k1) {
    uint64_t i = mixid(k0, k1) & dmask;
    for (;;) {
        ent_t *e = &dtab[i];
        if (e->k1 == (k1 | IDTAG)) { if (e->k0 == k0) return e; }
        else if (e->k1 == 0) return &miss_entry;
        i = (i + 1) & dmask;
    }
}

static void desc_payload(uint64_t k0, uint64_t k1, ent_t *e) {
    uint64_t h = mixid(k0, k1);
    e->off = h & 0xFFFFFFFULL;                        /* 28-bit slab offset */
    e->wa  = (1ULL << 20) + ((h >> 28) & 0xFFFFFULL); /* ~2^20 stride       */
    e->wab = e->wa * ((1ULL << 20) + ((h >> 48) & 0xFFFFULL)); /* ~2^40     */
}

/* ------------------------------------------------- direct (reference) id --- */
/* Recompute v = sum_i z_i b_i, canonicalise, build the full id from the 9
 * signed coordinates.  Used by pre-pass, reference implementation, stats. */
static inline void direct_id(const int *z, uint64_t *k0, uint64_t *k1,
                             int64_t *t0, int64_t *t1, int64_t *t2,
                             int64_t *vout, int64_t *maxabs) {
    int64_t v[D];
    for (int a = 0; a < D; a++) v[a] = 0;
    for (int i = 0; i < D; i++) {
        int zi = z[i];
        if (!zi) continue;
        const int *b = &LAT_BRED[i][0];
        for (int a = 0; a < D; a++) v[a] += (int64_t)zi * b[a];
    }
    if (vout) for (int a = 0; a < D; a++) vout[a] = v[a];
    if (maxabs) for (int a = 0; a < D; a++) {
        int64_t m = v[a] < 0 ? -v[a] : v[a];
        if (m > *maxabs) *maxabs = m;
    }
    int64_t fa = v[5] >> 63, fb = v[11] >> 63;
    for (int a = 0; a < 6; a++) v[a] = (v[a] ^ fa) - fa;
    for (int a = 6; a < D; a++) v[a] = (v[a] ^ fb) - fb;
    *k0 = (uint64_t)(v[0] + RADIX * (v[1] + RADIX * (v[2] + RADIX * (v[3]
        + RADIX * v[5]))));
    *k1 = (uint64_t)(v[6] + RADIX * (v[7] + RADIX * (v[9] + RADIX * v[11])));
    *t0 = v[4]; *t1 = v[8]; *t2 = v[10];
}

/* ------------------------------------------ v3: slab-id injectivity test --- */
#define PACK_TRIALS (1u << 20)
static uint64_t prng_s = 0x243F6A8885A308D3ULL;
static inline uint64_t prng(void) {
    uint64_t z = (prng_s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/* Pack the nine slab coordinates (v0,v1,v2,v3,v5 | v6,v7,v9,v11). */
static inline void pack9(const int64_t *v, int64_t radix, uint64_t *k0, uint64_t *k1) {
    *k0 = (uint64_t)(v[0] + radix * (v[1] + radix * (v[2] + radix * (v[3]
        + radix * v[5]))));
    *k1 = (uint64_t)(v[6] + radix * (v[7] + radix * (v[9] + radix * v[11])));
}

typedef struct { uint64_t k0, k1; uint64_t t; unsigned char used; } pkey_t;
static pkey_t *ptab;
static uint64_t ptabmask;

static int pput(uint64_t k0, uint64_t k1) {
    uint64_t i = mixid(k0, k1) & ptabmask;
    for (;;) {
        pkey_t *e = &ptab[i];
        if (!e->used) { e->used = 1; e->k0 = k0; e->k1 = k1; return 1; }
        if (e->k0 == k0 && e->k1 == k1) return 0;   /* repeated key = collision */
        i = (i + 1) & ptabmask;
    }
}

static int run_packcheck(void) {
    int64_t v1[12] = {749, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
    int64_t v2[12] = {-748, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1};
    uint64_t a0, a1, b0, b1;
    pack9(v1, 1497, &a0, &a1);
    pack9(v2, 1497, &b0, &b1);
    int alias_old = (a0 == b0 && a1 == b1);
    pack9(v1, RADIX, &a0, &a1);
    pack9(v2, RADIX, &b0, &b1);
    int distinct_new = !(a0 == b0 && a1 == b1);

    /* Production bound at A = 2^10.00 from the embedded covariance: 1699.1. */
    const int64_t B = 1700;
    int ok_range = (B <= (RADIX - 1) / 2);

    uint64_t nkeys = 1ULL << 21;
    ptab = calloc(nkeys, sizeof *ptab);
    if (!ptab) { perror("calloc ptab"); return 2; }
    ptabmask = nkeys - 1;

    uint64_t coll = 0;
    int64_t v[12];
    for (uint64_t it = 0; it < PACK_TRIALS; it++) {
        for (int j = 0; j < 12; j++)
            v[j] = (int64_t)(prng() % (uint64_t)(2 * B + 1)) - B;
        uint64_t k0, k1;
        pack9(v, RADIX, &k0, &k1);
        /* A repeated full tuple cannot occur in 2^20 draws; any repeated key is
         * therefore a genuine packing collision. */
        if (!pput(k0, k1)) coll++;
    }
    /* boundary tuples: each of the nine coordinates at +-B, others zero */
    for (int j = 0; j < 9; j++) {
        for (int s = -1; s <= 1; s += 2) {
            int idx[9] = {0, 1, 2, 3, 5, 6, 7, 9, 11};
            int64_t w[12] = {0};
            w[idx[j]] = s * B;
            uint64_t k0, k1;
            pack9(w, RADIX, &k0, &k1);
            if (!pput(k0, k1)) coll++;
        }
    }
    free(ptab);

    /* Exact in-slab cell with strides wa = 2B+1, wab = wa^2.  The map is a
     * balanced mixed-radix representation, hence injective; test it anyway. */
    int64_t wa = 2 * B + 1, wab = wa * wa;
    uint64_t ncell = 1ULL << 21, cellcoll = 0;
    pkey_t *ctab = calloc(ncell, sizeof *ctab);
    if (!ctab) { perror("calloc ctab"); return 2; }
    uint64_t cmask = ncell - 1;
    for (uint64_t it = 0; it < PACK_TRIALS; it++) {
        int64_t t0 = (int64_t)(prng() % (uint64_t)(2 * B + 1)) - B;
        int64_t t1 = (int64_t)(prng() % (uint64_t)(2 * B + 1)) - B;
        int64_t t2 = (int64_t)(prng() % (uint64_t)(2 * B + 1)) - B;
        int64_t cell = t0 + wa * t1 + wab * t2;
        uint64_t i = mixid((uint64_t)cell, 0x5BD1E995ULL) & cmask;
        for (;;) {
            pkey_t *e = &ctab[i];
            if (!e->used) {
                e->used = 1; e->k0 = (uint64_t)cell;
                e->t = (uint64_t)t0; e->k1 = (uint64_t)t1;
                break;
            }
            if (e->k0 == (uint64_t)cell) {
                if (e->t != (uint64_t)t0 || e->k1 != (uint64_t)t1) cellcoll++;
                break;
            }
            i = (i + 1) & cmask;
        }
    }
    free(ctab);

    printf("packcheck: production bound |v_i|<=%lld, radix %lld covers +/-%lld %s\n",
           (long long)B, (long long)RADIX, (long long)(RADIX - 1) / 2,
           ok_range ? "ok" : "FAIL");
    printf("  v2 alias with radix 1497 (must be 1): %d\n", alias_old);
    printf("  v3 ids distinct for the same pair (must be 1): %d\n", distinct_new);
    printf("  random/boundary id collisions: %llu of %llu samples\n",
           (unsigned long long)coll, (unsigned long long)PACK_TRIALS);
    printf("  exact-stride in-slab collisions: %llu (wa=%lld, wab=%lld)\n",
           (unsigned long long)cellcoll, (long long)wa, (long long)wab);
    return (ok_range && alias_old && distinct_new && coll == 0 && cellcoll == 0) ? 0 : 1;
}

static void dinsert(uint64_t k0, uint64_t k1) {
    uint64_t i = mixid(k0, k1) & dmask;
    for (;;) {
        ent_t *e = &dtab[i];
        if (e->k1 == 0) {
            e->k0 = k0; e->k1 = k1 | IDTAG;
            desc_payload(k0, k1, e);
            return;
        }
        if (e->k1 == (k1 | IDTAG) && e->k0 == k0) return;
        i = (i + 1) & dmask;
    }
}

static void build_route(void) {
    route = malloc(NNODES);
    if (!route) { perror("malloc route"); exit(1); }
    for (uint64_t i = 0; i < NNODES; i++) {
        uint64_t x = i + 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        x ^= x >> 31;
        route[i] = (uint8_t)(x & (NDEST - 1));
    }
}

/* ------------------------------------------------------- per-hit work ---- */
/* The macros below run with the carried ambient state s0..s6 in scope:
 *   s0 = packed a-word (v0,v1,v2,v3,v5), s1 = packed b-word (v6,v7,v9,v11),
 *   s2 = v5 (a5 sign key), s3 = v11 (b5 sign key),
 *   s4 = v4, s5 = v8, s6 = v10 (in-slab coordinates). */
#define CANON { \
    int64_t fa = s2 >> 63, fb = s3 >> 63;                                    \
    uint64_t k0 = (uint64_t)((s0 ^ fa) - fa);                                \
    uint64_t k1 = (uint64_t)((s1 ^ fb) - fb);                                \
    int64_t t0 = (s4 ^ fa) - fa, t1 = (s5 ^ fb) - fb, t2 = (s6 ^ fb) - fb;

#define HIT_SINGLE                                                           \
    CANON                                                                    \
    ent_t *e = dlookup(k0, k1);                                              \
    uint64_t cell = e->off + (uint64_t)t0 + e->wa * (uint64_t)t1             \
                  + e->wab * (uint64_t)t2;                                   \
    bucket8[cell & b8mask] += (unsigned char)logp; }

#define HIT_SEND                                                             \
    CANON                                                                    \
    uint64_t h = mixid(k0, k1);                                              \
    uint32_t node = (uint32_t)((h >> 24) & (NNODES - 1));                    \
    uint32_t dst = route[node];                                              \
    queue_t *q = &msgbuf[dst];                                               \
    msg_t *m = &q->m[q->n];                                                  \
    m->k0 = k0; m->k1 = k1; m->t0 = t0; m->t1 = t1; m->t2 = t2; m->lp = logp;\
    q->n++; }

/* Miss-cost variant: every key is made absent by flipping an unused high bit of
 * the b-word (the id uses only ~43 of its 64 bits), so the probe is always
 * unsuccessful and ends at an empty slot.  Branch-free corruption. */
#define HIT_PROBE_OPT                                                        \
    CANON                                                                    \
    ent_t *e = dlookup_opt(k0, k1);                                          \
    uint64_t cell = e->off + (uint64_t)t0 + e->wa * (uint64_t)t1             \
                  + e->wab * (uint64_t)t2;                                   \
    bucket8[cell & b8mask] += (unsigned char)logp; }

#define HIT_MISS                                                             \
    CANON                                                                    \
    ent_t *e = dlookup_opt(k0, k1 ^ (1ULL << 45));                           \
    uint64_t cell = e->off + (uint64_t)t0 + e->wa * (uint64_t)t1             \
                  + e->wab * (uint64_t)t2;                                   \
    bucket8[cell & b8mask] += (unsigned char)logp; }

/* Apply at the owner: consume the full id, probe the descriptor table,
 * recompute the 64-bit cell, update the bucket. */
static void apply_all(void) {
    for (int d = 0; d < NDEST; d++) {
        queue_t *q = &msgbuf[d];
        for (uint32_t i = 0; i < q->n; i++) {
            const msg_t *m = &q->m[i];
            ent_t *e = dlookup(m->k0, m->k1);
            uint64_t cell = e->off + (uint64_t)m->t0
                          + e->wa * (uint64_t)m->t1 + e->wab * (uint64_t)m->t2;
            bucket8[cell & b8mask] += (unsigned char)m->lp;
        }
        q->n = 0;
    }
}

/* ------------------------------------------------------------ enumeration */
#define BND_EPS 1e-9

/* Generic DFS (unchanged from v1), used by the reference path and pre-pass. */
static long long g_nodes;
static long long enum_direct(double R2, int do_insert) {
    int z[D], lo[D], hi[D];
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;
    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    g_nodes = 0;
    g_hits = 0;

    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            g_nodes++;
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;
            if (k == 0) {
                g_hits++;
                uint64_t k0, k1; int64_t t0, t1, t2;
                direct_id(z, &k0, &k1, &t0, &t1, &t2, NULL, NULL);
                if (do_insert) dinsert(k0, k1);
                else {
                    ent_t *e = dlookup(k0, k1);
                    uint64_t cell = e->off + (uint64_t)t0
                                  + e->wa * (uint64_t)t1 + e->wab * (uint64_t)t2;
                    bucket8[cell & b8mask] += (unsigned char)logp;
                }
            } else {
                rem[k - 1] = r;
                double zz = (double)z[k];
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
                z[k] = lo[k] - 1;
            }
        } else {
            if (k == D - 1) break;
            k++;
        }
    }
    return g_hits;
}

/* Optimised DFS: generic levels >= 1; level 0 is a tight register loop over
 * its exact interval, carrying the 7 ambient state words.  Identical walk to
 * v1's FAST_BODY. */
#define FAST_BODY(HITSTMT)                                                                  \
    int z[D], lo[D], hi[D];                                                                 \
    double C[D][D], rem[D];                                                                 \
    int64_t ST[D][7];                                                                       \
    int k = D - 1, i;                                                                       \
    double rr;                                                                              \
    long long hits = 0;                                                                     \
    for (i = 0; i < D; i++) C[k][i] = CST[i];                                               \
    rem[k] = R2;                                                                            \
    rr = sqrt(rem[k] / BN[k]);                                                              \
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);                                              \
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);                                             \
    z[k] = lo[k] - 1;                                                                       \
    for (int q = 0; q < 7; q++) ST[k][q] = (int64_t)(lo[k] - 1) * DLT[k][q];                \
    const int64_t d0 = DLT[0][0], d1 = DLT[0][1], d2 = DLT[0][2], d3 = DLT[0][3],           \
                  d4 = DLT[0][4], d5 = DLT[0][5], d6 = DLT[0][6];                           \
    for (;;) {                                                                              \
        if (z[k] < hi[k]) {                                                                 \
            z[k]++;                                                                         \
            for (int q = 0; q < 7; q++) ST[k][q] += DLT[k][q];                              \
            double t = C[k][k] - (double)z[k];                                              \
            double r = rem[k] - t * t * BN[k];                                              \
            if (r < 0.0) r = 0.0;                                                           \
            double zz = (double)z[k];                                                       \
            if (k == 1) {                                                                   \
                double c0 = C[1][0] - zz * MU[1][0];                                        \
                double r0 = sqrt(r / BN[0]);                                                \
                int l0 = (int)ceil(c0 - r0 - BND_EPS), h0 = (int)floor(c0 + r0 + BND_EPS);  \
                int64_t s0 = ST[1][0] + l0 * d0, s1 = ST[1][1] + l0 * d1,                   \
                        s2 = ST[1][2] + l0 * d2, s3 = ST[1][3] + l0 * d3,                   \
                        s4 = ST[1][4] + l0 * d4, s5 = ST[1][5] + l0 * d5,                   \
                        s6 = ST[1][6] + l0 * d6;                                            \
                for (int j = l0; j <= h0; j++) {                                            \
                    HITSTMT;                                                                \
                    s0 += d0; s1 += d1; s2 += d2; s3 += d3; s4 += d4; s5 += d5; s6 += d6;   \
                }                                                                           \
                if (h0 >= l0) hits += h0 - l0 + 1;                                          \
            } else {                                                                        \
                rem[k - 1] = r;                                                             \
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];              \
                k--;                                                                        \
                rr = sqrt(rem[k] / BN[k]);                                                  \
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);                                  \
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);                                 \
                z[k] = lo[k] - 1;                                                           \
                for (int q = 0; q < 7; q++)                                                 \
                    ST[k][q] = ST[k + 1][q] + (int64_t)(lo[k] - 1) * DLT[k][q];             \
            }                                                                               \
        } else {                                                                            \
            if (k == D - 1) break;                                                          \
            k++;                                                                            \
        }                                                                                   \
    }                                                                                       \
    return hits;

/* True "walk only, ambient state carried, no hit work": an empty asm with the
 * 7 state words as register inputs forces them live at each hit without adding
 * any instruction (Google-Benchmark DoNotOptimize idiom). */
static long long enum_fast_walk_live(double R2) {
    FAST_BODY(__asm__ __volatile__("" :: "r"(s0), "r"(s1), "r"(s2), "r"(s3),
                                   "r"(s4), "r"(s5), "r"(s6)))
}
/* v1's published methodology: a dead xor expression.  clang -O2 removes the
 * whole state chain, so this number is the plain enumeration floor. */
static long long enum_fast_walk_dead(double R2) {
    FAST_BODY((void)(s0 ^ s1 ^ s2 ^ s3 ^ s4 ^ s5 ^ s6))
}
static long long enum_fast_single(double R2) { FAST_BODY(HIT_SINGLE) }
static long long enum_fast_send(double R2)   { FAST_BODY(HIT_SEND) }
static long long enum_fast_probe_opt(double R2) { FAST_BODY(HIT_PROBE_OPT) }
static long long enum_fast_miss(double R2)   { FAST_BODY(HIT_MISS) }

/* ------------------------------------------------ carried-state checking -- */
static long long g_stbad;
static int64_t g_maxabs;

static void check_state(const int *z, const int64_t *st) {
    int64_t r[D];
    for (int a = 0; a < D; a++) r[a] = 0;
    for (int i = 0; i < D; i++) {
        int zi = z[i];
        if (!zi) continue;
        const int *b = &LAT_BRED[i][0];
        for (int a = 0; a < D; a++) r[a] += (int64_t)zi * b[a];
    }
    for (int a = 0; a < D; a++) {
        int64_t m = r[a] < 0 ? -r[a] : r[a];
        if (m > g_maxabs) g_maxabs = m;
    }
    int64_t pa = r[0] + RADIX * (r[1] + RADIX * (r[2] + RADIX * (r[3]
        + RADIX * r[5])));
    int64_t pb = r[6] + RADIX * (r[7] + RADIX * (r[9] + RADIX * r[11]));
    if (st[0] != pa || st[1] != pb || st[2] != r[5] || st[3] != r[11]
        || st[4] != r[4] || st[5] != r[8] || st[6] != r[10]) g_stbad++;
}

/* Generic DFS with carried state, checking the state at every hit. */
static long long enum_check(double R2) {
    int z[D], lo[D], hi[D];
    int64_t ST[D + 1][7];
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;
    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    for (int q = 0; q < 7; q++) ST[k][q] = (int64_t)(lo[k] - 1) * DLT[k][q];
    g_hits = 0;
    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            for (int q = 0; q < 7; q++) ST[k][q] += DLT[k][q];
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;
            if (k == 0) { g_hits++; check_state(z, ST[0]); }
            else {
                rem[k - 1] = r;
                double zz = (double)z[k];
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
                z[k] = lo[k] - 1;
                for (int q = 0; q < 7; q++)
                    ST[k][q] = ST[k + 1][q] + (int64_t)(lo[k] - 1) * DLT[k][q];
            }
        } else {
            if (k == D - 1) break;
            k++;
        }
    }
    return g_hits;
}

/* ------------------------------------------------------------ statistics -- */
static uint64_t g_depth_hist[128], g_depth_max, g_depth_sum, g_nprobe;
static uint64_t g_miss_hist[128], g_miss_max, g_miss_sum, g_nmissprobe, g_miss_found;
static uint64_t g_dest_hist[NDEST], g_node_used[NNODES];

static uint64_t probe_depth(uint64_t k0, uint64_t k1, int *found) {
    uint64_t i = mixid(k0, k1) & dmask, d = 1;
    for (;;) {
        ent_t *e = &dtab[i];
        if (e->k1 == (k1 | IDTAG)) { if (e->k0 == k0) { *found = 1; return d; } }
        else if (e->k1 == 0) { *found = 0; return d; }
        i = (i + 1) & dmask;
        d++;
    }
}

static long long enum_stats(double R2) {
    int z[D], lo[D], hi[D];
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;
    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    g_hits = 0;
    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;
            if (k == 0) {
                g_hits++;
                uint64_t k0, k1; int64_t t0, t1, t2;
                direct_id(z, &k0, &k1, &t0, &t1, &t2, NULL, NULL);
                int found;
                uint64_t d = probe_depth(k0, k1, &found);
                g_depth_sum += d; g_nprobe++;
                g_depth_hist[d < 128 ? d : 127]++;
                if (d > g_depth_max) g_depth_max = d;
                {
                    uint64_t cm = k1 ^ (1ULL << 45);   /* absent by construction */
                    uint64_t dm = probe_depth(k0, cm, &found);
                    g_miss_sum += dm; g_nmissprobe++;
                    if (found) g_miss_found++;
                    g_miss_hist[dm < 128 ? dm : 127]++;
                    if (dm > g_miss_max) g_miss_max = dm;
                }
                uint64_t h = mixid(k0, k1);
                uint32_t node = (uint32_t)((h >> 24) & (NNODES - 1));
                g_node_used[node] = 1;
                g_dest_hist[route[node]]++;
            } else {
                rem[k - 1] = r;
                double zz = (double)z[k];
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
                z[k] = lo[k] - 1;
            }
        } else {
            if (k == D - 1) break;
            k++;
        }
    }
    return g_hits;
}

/* ------------------------------------------------------------------ misc -- */
static uint64_t digest(void) {
    unsigned long long s = 1469598103934665603ULL;
    for (size_t i = 0; i < BN_BYTES; i++) s = (s ^ bucket8[i]) * 1099511628211ULL;
    return s;
}
static uint64_t digest_of(const unsigned char *p) {
    unsigned long long s = 1469598103934665603ULL;
    for (size_t i = 0; i < BN_BYTES; i++) s = (s ^ p[i]) * 1099511628211ULL;
    return s;
}

static long long run_mode(int hm, double R2, int apply) {
    long long h = 0;
    for (int d = 0; d < NDEST; d++) msgbuf[d].n = 0;
    if (hm == 20) h = enum_fast_walk_live(R2);
    else if (hm == 21) h = enum_fast_walk_dead(R2);
    else if (hm == 23) h = enum_fast_single(R2);
    else if (hm == 24) { h = enum_fast_send(R2); if (apply) apply_all(); }
    else if (hm == 26) h = enum_fast_miss(R2);
    else if (hm == 27) h = enum_fast_probe_opt(R2);
    else if (hm == 3) h = enum_direct(R2, 0);
    else { fprintf(stderr, "bad hmode %d\n", hm); exit(2); }
    return h;
}

static void alloc_tables(void) {
    dtab = calloc(DN, sizeof *dtab);
    if (!dtab) { perror("calloc dtab"); exit(1); }
    dmask = DN - 1;
    bucket8 = calloc(BN_BYTES, 1);
    if (!bucket8) { perror("calloc bucket"); exit(1); }
    b8mask = BN_BYTES - 1;
    for (int d = 0; d < NDEST; d++) {
        msgbuf[d].m = malloc((size_t)BUFCAP * sizeof(msg_t));
        if (!msgbuf[d].m) { perror("malloc msgbuf"); exit(1); }
        /* Pre-touch every page once, in every invocation including the 0-rep
         * baseline: otherwise the first repetition pays ~43M instructions of
         * page-fault work and the 0-vs-8 delta is inflated by ~2.5/pt.  The
         * stores are volatile so the compiler cannot fold the touch into the
         * allocation (a plain malloc+memset becomes calloc, which leaves the
         * pages unbacked and defeats the pre-touch). */
        {
            volatile unsigned char *p = (volatile unsigned char *)msgbuf[d].m;
            for (size_t off = 0; off < (size_t)BUFCAP * sizeof(msg_t); off += 4096)
                p[off] = 0;
        }
        msgbuf[d].n = 0;
    }
    build_route();
}

int main(int argc, char **argv) {
    const char *mode = argc > 1 ? argv[1] : "stcheck";
    int reps = argc > 2 ? atoi(argv[2]) : 1;
    double R2 = LAT_R2;
    for (int i = 2; i < argc; i++)
        if (!strncmp(argv[i], "R2=", 3)) R2 = atof(argv[i] + 3);
    if (!strcmp(mode, "packcheck")) return run_packcheck();   /* v3 */
    setup();
    setup_dlt();
    alloc_tables();

    if (!strcmp(mode, "stcheck")) {
        long long h = enum_check(R2);
        printf("stcheck hits=%lld state mismatches=%lld max|v_i|=%lld (id bound %d)\n",
               h, g_stbad, (long long)g_maxabs, AMAX);
        if (g_maxabs > AMAX) printf("WARNING: coordinate out of id range\n");
        return g_stbad != 0 || g_maxabs > AMAX;
    }

    if (!strcmp(mode, "build")) {   /* build descriptor table (pre-pass) */
        g_hits = 0;
        long long h = enum_direct(R2, 1);
        size_t used = 0, maxprobe = 0;
        (void)maxprobe;
        for (size_t i = 0; i < DN; i++) if (dtab[i].k1) used++;
        printf("prepass hits=%lld distinct_ids=%zu table=%zu load=%.4f\n",
               h, used, DN, (double)used / (double)DN);
        return 0;
    }

    if (!strcmp(mode, "stats")) {
        g_hits = 0;
        long long h = enum_direct(R2, 1);          /* build table */
        size_t used = 0;
        for (size_t i = 0; i < DN; i++) if (dtab[i].k1) used++;
        memset(g_depth_hist, 0, sizeof g_depth_hist);
        memset(g_miss_hist, 0, sizeof g_miss_hist);
        memset(g_dest_hist, 0, sizeof g_dest_hist);
        memset(g_node_used, 0, sizeof g_node_used);
        g_depth_max = g_depth_sum = g_nprobe = 0;
        g_miss_max = g_miss_sum = g_nmissprobe = g_miss_found = 0;
        g_hits = 0;
        long long h2 = enum_stats(R2);
        uint64_t nodes_used = 0;
        for (uint32_t i = 0; i < NNODES; i++) nodes_used += g_node_used[i];
        printf("hits=%lld (pass1 %lld) distinct_ids=%zu load=%.4f\n",
               h2, h, used, (double)used / (double)DN);
        printf("descriptor probe (present keys): mean=%.4f max=%llu\n",
               (double)g_depth_sum / (double)g_nprobe, (unsigned long long)g_depth_max);
        for (int i = 1; i < 16; i++)
            if (g_depth_hist[i]) printf("  depth %2d: %llu\n", i,
                                        (unsigned long long)g_depth_hist[i]);
        printf("absent-key probes (bit45 flip, all hits): n=%llu mean=%.4f max=%llu "
               "accidental-key-found=%llu\n",
               (unsigned long long)g_nmissprobe,
               (double)g_miss_sum / (double)g_nmissprobe,
               (unsigned long long)g_miss_max, (unsigned long long)g_miss_found);
        for (int i = 1; i < 16; i++)
            if (g_miss_hist[i]) printf("  depth %2d: %llu\n", i,
                                       (unsigned long long)g_miss_hist[i]);
        printf("route table: %llu nodes touched of %u, dest histogram:",
               (unsigned long long)nodes_used, (unsigned)NNODES);
        for (int d = 0; d < NDEST; d++) printf(" %llu", (unsigned long long)g_dest_hist[d]);
        printf("\n");
        uint64_t mx = 0;
        for (int d = 0; d < NDEST; d++) if (g_dest_hist[d] > mx) mx = g_dest_hist[d];
        printf("max messages per destination (one rep) = %llu, buffer cap = %u  %s\n",
               (unsigned long long)mx, (unsigned)BUFCAP, mx <= BUFCAP ? "OK" : "OVERFLOW");
        return mx <= BUFCAP ? 0 : 1;
    }

    if (!strcmp(mode, "refcheck")) {
        unsigned char *bref = calloc(BN_BYTES, 1);
        unsigned char *bsin = calloc(BN_BYTES, 1);
        unsigned char *bdis = calloc(BN_BYTES, 1);
        if (!bref || !bsin || !bdis) { perror("calloc"); exit(1); }
        g_hits = 0;
        long long hpre = enum_direct(R2, 1);
        printf("prepass hits=%lld\n", hpre);

        long long hr = 0, hs = 0, hd = 0;
        for (int r = 0; r < reps; r++) {
            bucket8 = bref;
            memset(bref, 0, BN_BYTES);
            hr += enum_direct(R2, 0);

            bucket8 = bsin;
            memset(bsin, 0, BN_BYTES);
            hs += enum_fast_single(R2);

            bucket8 = bdis;
            memset(bdis, 0, BN_BYTES);
            int of = 0;
            for (int d = 0; d < NDEST; d++) if (msgbuf[d].n) of = 1;
            (void)of;
            for (int d = 0; d < NDEST; d++) msgbuf[d].n = 0;
            hd += enum_fast_send(R2);
            uint64_t mx = 0;
            for (int d = 0; d < NDEST; d++) if (msgbuf[d].n > mx) mx = msgbuf[d].n;
            if (mx > BUFCAP) { printf("FATAL buffer overflow %llu\n", (unsigned long long)mx); return 1; }
            apply_all();
        }
        int ok1 = memcmp(bref, bsin, BN_BYTES) == 0;
        int ok2 = memcmp(bref, bdis, BN_BYTES) == 0;
        printf("reps=%d hits ref=%lld single=%lld dist=%lld\n", reps, hr, hs, hd);
        printf("digest ref=%016llx single=%016llx dist=%016llx\n",
               (unsigned long long)digest_of(bref),
               (unsigned long long)digest_of(bsin),
               (unsigned long long)digest_of(bdis));
        printf("bucket arrays: ref==single %s, ref==dist %s\n",
               ok1 ? "IDENTICAL" : "DIFFER", ok2 ? "IDENTICAL" : "DIFFER");
        free(bref); free(bsin); free(bdis);
        bucket8 = calloc(BN_BYTES, 1);
        return (ok1 && ok2 && hr == hs && hs == hd) ? 0 : 1;
    }

    if (!strcmp(mode, "icount")) {
        int hm = argc > 3 ? atoi(argv[3]) : 20;
        int apply = argc > 4 ? atoi(argv[4]) : 0;
        g_hits = 0;
        enum_direct(R2, 1);                 /* pre-pass: fill descriptor table */
        memset(bucket8, 0, BN_BYTES);
        long long h = 0;
        for (int r = 0; r < reps; r++) {
            for (int d = 0; d < NDEST; d++) msgbuf[d].n = 0;
            h += run_mode(hm, R2, apply);
        }
        uint64_t dg = (hm == 20 || (hm == 24 && !apply)) ? digest() : digest();
        printf("reps=%d hits=%lld digest=%016llx\n", reps, h, (unsigned long long)dg);
        return 0;
    }

    if (!strcmp(mode, "bench")) {
        long long h = enum_direct(R2, 0);
        printf("reference hits=%lld\n", h);
        uint64_t t0 = now_ns();
        long long h2 = enum_fast_single(R2);
        uint64_t t1 = now_ns();
        printf("single hits=%lld  %.2f ns/pt\n", h2,
               (double)(t1 - t0) / (double)h2);
        return 0;
    }

    if (!strcmp(mode, "count")) {
        g_hits = 0;
        long long h = enum_direct(R2, 1);
        printf("%lld\nnodes=%lld nodes/pt=%.3f\n", h, g_nodes,
               (double)g_nodes / (double)h);
        return 0;
    }

    fprintf(stderr, "unknown mode %s\n", mode);
    return 2;
}
