/* ARCH-D hit-cost microbenchmark.
 *
 * Measures the amortised AArch64 cost per enumerated point of a
 * Schnorr-Euchner-style DFS enumeration of a random coset of a 12-dimensional
 * integer lattice (covolume 2^40) inside a 12-dimensional ellipsoid, including
 * the per-hit work: a slab id built from 9 of the 12 coordinates and one
 * bucket[sid] += logp store.
 *
 * The enumeration is the standard GSO/Cholesky DFS with incremental partial
 * sums: the centre vector is carried down one level at a time
 * (C[k-1][i] = C[k][i] - z_k * mu[k][i]) and the remaining squared radius is
 * decremented by bn[k] * (C[k][k] - z_k)^2.  A node at level k is skipped
 * unless (C[k][k]-z_k)^2 bn[k] <= rem[k], so the tree is pruned by the radius
 * exactly as in Schnorr-Euchner / Fincke-Pohst enumeration.
 *
 * Build (see make.sh):
 *   clang -O2 -o enum12 enum.c -lm
 *   clang -O2 -DDIM6 ... (via -DCDATA_H='"cdata_lat_d6.h"')
 *
 * Modes:
 *   ./enum12 bench [reps]      time enum-only / enum+hit / hit-only
 *   ./enum12 hist              node-count histogram per level (reported)
 *   ./enum12 check             count, duplicates, direct-form containment
 *   ./enum12 count             print hit count only (fpylll cross-check)
 *   ./enum12 dump N out.txt    dump N sampled points (reduced-basis coords)
 *   ./enum12 verify M          small d only: full set vs brute-force box scan
 *   ./enum12 brute M           brute-force count in box [-M,M]^d
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef CDATA_H
#define CDATA_H "archd_lattice_d12.h"
#endif
#include CDATA_H

#define D LAT_D

/* ---------------------------------------------------------------- globals */
static double MU[D][D], BN[D], CST[D];
static unsigned short *bucket;
static uint32_t bmask;
static unsigned short logp = 137;
static long long g_nodes, g_hits;
static long long g_hist[D];

static const int *BRED = &LAT_BRED[0][0];

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
    /* Gram-Schmidt vectors (real coordinates) for the coset centre */
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

/* ------------------------------------------------------- per-hit work ---- */
/* Slab id from 9 of the 12 coordinates: 9 offset+mask ops, a couple of
 * multiplies, one bucket store.  (Variant used for the central estimate.) */
static inline void hit9(const int *z) {
    uint32_t g0 = (((uint32_t)(z[0] + 512)) & 511u) + (((uint32_t)(z[1] + 512)) & 511u)
                + (((uint32_t)(z[2] + 512)) & 511u);
    uint32_t g1 = (((uint32_t)(z[3] + 512)) & 511u) + (((uint32_t)(z[4] + 512)) & 511u)
                + (((uint32_t)(z[5] + 512)) & 511u);
    uint32_t g2 = (((uint32_t)(z[6] + 512)) & 511u) + (((uint32_t)(z[7] + 512)) & 511u)
                + (((uint32_t)(z[8] + 512)) & 511u);
    uint32_t sid = ((g0 * 1000003u + g1) * 1000033u + g2) & bmask;
    bucket[sid] += logp;
}
/* Leaner 3-coordinate variant, for the reported range. */
static inline void hit3(const int *z) {
    uint32_t sid = (((((uint32_t)(z[0] + 512)) & 511u) * 61u
                   + (((uint32_t)(z[1] + 512)) & 511u)) * 61u
                   + (((uint32_t)(z[2] + 512)) & 511u)) & bmask;
    bucket[sid] += logp;
}


/* ------------------------------------------- realistic ARCH-D hit (added) ---
 * The hit must be addressed by the AMBIENT point v = sum_i z_i BRED[i] (the (a,b)
 * coefficient vector), not by the lattice coordinates z.  Per enumeration level we
 * carry 7 integers, all affine in v, updated with one add each per node visit:
 *   SA, SB : dense slab index split into its a-part (coords 0,1,5) and b-part (6..11)
 *   A5, B5 : the sign keys a5 = v[5], b5 = v[11]
 *   T2..T4 : the 3 in-slab coordinates v[2], v[3], v[4]
 * At the hit, the point is canonicalised to its orbit representative under
 * (a,b) -> (-a,-b) and the automorphism (a,b) -> (a,-b): flip the a-part if a5 < 0 and
 * the b-part if b5 < 0 (affine: X(-a) = CX - X(a)).  The canonical slab descriptor
 * {offset, w0, w01} is loaded and the cell is updated with a byte read-modify-write. */
#define NST 7
typedef struct { int64_t off; int64_t w0, w01; } slabdesc_t;
static slabdesc_t *slabtab;
static uint64_t smask;
static int64_t DLT[D][NST];             /* per-basis-vector deltas of the 7 state words */
static const int64_t CSA = 1000003, CSB = 7000001;
static unsigned char *bucket8;
static uint64_t b8mask;
/* distributed variant: each hit is appended to an outgoing buffer instead */
typedef struct { uint64_t dst; uint32_t cell; uint32_t lp; } msg_t;
static msg_t *outbuf; static uint64_t outcap, outn;

static void setup_real(void) {
    static const int64_t R[12] = {1, 2311, 0, 0, 0, 5340721, 1, 2311, 5340721, 12342406231LL,
                                  28523000000000LL, 65916600000000000LL};
    for (int k = 0; k < D; k++) {
        const int *b = &LAT_BRED[k][0];
        DLT[k][0] = b[0] * R[0] + b[1] * R[1] + b[5] * R[5];
        DLT[k][1] = b[6] * R[6] + b[7] * R[7] + b[8] * R[8] + b[9] * R[9] + b[10] * R[10] + b[11] * R[11];
        DLT[k][2] = b[5]; DLT[k][3] = b[11];
        DLT[k][4] = b[2]; DLT[k][5] = b[3]; DLT[k][6] = b[4];
    }
    size_t ns = (size_t)1 << 16;
    slabtab = malloc(ns * sizeof *slabtab); smask = ns - 1;
    for (size_t i = 0; i < ns; i++) { slabtab[i].off = (int64_t)(i * 2654435761u) & 0xffffff;
        slabtab[i].w0 = 1500 + (i & 255); slabtab[i].w01 = slabtab[i].w0 * (1400 + (i >> 8 & 255)); }
    size_t nb = (size_t)1 << 24; bucket8 = calloc(nb, 1); b8mask = nb - 1;
    outcap = (size_t)1 << 23; outbuf = malloc(outcap * sizeof *outbuf);
}

static inline void hit_real(const int64_t *st) {
    int64_t fa = st[2] >> 63, fb = st[3] >> 63;          /* -1 if negative */
    int64_t sa = fa ? CSA - st[0] : st[0];
    int64_t sb = fb ? CSB - st[1] : st[1];
    const slabdesc_t *d = &slabtab[(uint64_t)(sa + sb) & smask];
    int64_t t2 = (st[4] ^ fa) - fa, t3 = (st[5] ^ fa) - fa, t4 = (st[6] ^ fa) - fa;
    uint64_t cell = (uint64_t)(d->off + t2 + d->w0 * t3 + d->w01 * t4);
    bucket8[cell & b8mask] += (unsigned char)logp;     /* real: no mask (descriptor bounds) */
}
static inline void hit_send(const int64_t *st) {
    int64_t fa = st[2] >> 63, fb = st[3] >> 63;
    int64_t sa = fa ? CSA - st[0] : st[0];
    int64_t sb = fb ? CSB - st[1] : st[1];
    const slabdesc_t *d = &slabtab[(uint64_t)(sa + sb) & smask];
    int64_t t2 = (st[4] ^ fa) - fa, t3 = (st[5] ^ fa) - fa, t4 = (st[6] ^ fa) - fa;
    msg_t *m = &outbuf[outn++ & (outcap - 1)];
    m->dst = (uint64_t)(sa + sb); m->cell = (uint32_t)(d->off + t2 + d->w0 * t3 + d->w01 * t4); m->lp = logp;
}
/* receiver side of the distributed variant: apply a delivered buffer */
static void apply_buf(const msg_t *m, uint64_t n) {
    for (uint64_t i = 0; i < n; i++) bucket8[m[i].cell & b8mask] += (unsigned char)m[i].lp;
}

/* ------------------------------------------------------------ enumeration */
/* tolerance on the interval bounds: lattice points exactly on the
 * ellipsoid boundary are otherwise lost to rounding (see FINDINGS). */
#define BND_EPS 1e-9

#define HITMODE_NONE 0
#define HITMODE_HIT9 1
#define HITMODE_HIT3 2
#define HITMODE_REAL 3
#define HITMODE_SEND 4

static long long g_stbad;
static void check_state(const int *z, const int64_t *st) {
    static const int64_t R[12] = {1, 2311, 0, 0, 0, 5340721, 1, 2311, 5340721, 12342406231LL,
                                  28523000000000LL, 65916600000000000LL};
    int64_t v[D] = {0};
    for (int i = 0; i < D; i++) for (int a = 0; a < D; a++) v[a] += (int64_t)z[i] * LAT_BRED[i][a];
    int64_t sa = v[0] * R[0] + v[1] * R[1] + v[5] * R[5];
    int64_t sb = v[6] * R[6] + v[7] * R[7] + v[8] * R[8] + v[9] * R[9] + v[10] * R[10] + v[11] * R[11];
    if (st[0] != sa || st[1] != sb || st[2] != v[5] || st[3] != v[11] || st[4] != v[2] || st[5] != v[3] || st[6] != v[4]) g_stbad++;
}
static long long enumerate(double R2, int hitmode, int do_hist) {
    int z[D], lo[D], hi[D];
    int64_t ST[D + 1][NST];
    int track = hitmode >= HITMODE_REAL;
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;

    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    if (track) for (int s = 0; s < NST; s++) ST[k][s] = (int64_t)(lo[k] - 1) * DLT[k][s];
    if (do_hist) memset(g_hist, 0, sizeof g_hist);
    g_nodes = g_hits = 0;

    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            g_nodes++;
            if (do_hist) g_hist[k]++;
            if (track) for (int s = 0; s < NST; s++) ST[k][s] += DLT[k][s];
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;          /* boundary rounding only */
            if (k == 0) {
                g_hits++;
                if (hitmode == HITMODE_HIT9) hit9(z);
                else if (hitmode == HITMODE_HIT3) hit3(z);
                else if (hitmode == HITMODE_REAL) hit_real(ST[0]);
                else if (hitmode == HITMODE_SEND) hit_send(ST[0]);
                if (do_hist && hitmode >= HITMODE_REAL) check_state(z, ST[0]);
            } else {
                rem[k - 1] = r;
                double zz = (double)z[k];
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
                z[k] = lo[k] - 1;
                if (track) for (int s = 0; s < NST; s++) ST[k][s] = ST[k + 1][s] + (int64_t)(lo[k] - 1) * DLT[k][s];
            }
        } else {
            if (k == D - 1) break;
            k++;
        }
    }
    return g_hits;
}

static long long enum_track_only(double R2) {
    const int do_hist = 0;
    int z[D], lo[D], hi[D];
    int64_t ST[D + 1][NST];
    
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;

    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    for (int s = 0; s < NST; s++) ST[k][s] = (int64_t)(lo[k] - 1) * DLT[k][s];
    if (do_hist) memset(g_hist, 0, sizeof g_hist);
    g_nodes = g_hits = 0;

    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            g_nodes++;
            if (do_hist) g_hist[k]++;
            for (int s = 0; s < NST; s++) ST[k][s] += DLT[k][s];
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;          /* boundary rounding only */
            if (k == 0) {
                g_hits++;
                /* state carried, no hit work */ (void)0;
            } else {
                rem[k - 1] = r;
                double zz = (double)z[k];
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
                z[k] = lo[k] - 1;
                for (int s = 0; s < NST; s++) ST[k][s] = ST[k + 1][s] + (int64_t)(lo[k] - 1) * DLT[k][s];
            }
        } else {
            if (k == D - 1) break;
            k++;
        }
    }
    return g_hits;
}

static long long enum_hit_real(double R2) {
    const int do_hist = 0;
    int z[D], lo[D], hi[D];
    int64_t ST[D + 1][NST];
    
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;

    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    for (int s = 0; s < NST; s++) ST[k][s] = (int64_t)(lo[k] - 1) * DLT[k][s];
    if (do_hist) memset(g_hist, 0, sizeof g_hist);
    g_nodes = g_hits = 0;

    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            g_nodes++;
            if (do_hist) g_hist[k]++;
            for (int s = 0; s < NST; s++) ST[k][s] += DLT[k][s];
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;          /* boundary rounding only */
            if (k == 0) {
                g_hits++;
                hit_real(ST[0]);
            } else {
                rem[k - 1] = r;
                double zz = (double)z[k];
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
                z[k] = lo[k] - 1;
                for (int s = 0; s < NST; s++) ST[k][s] = ST[k + 1][s] + (int64_t)(lo[k] - 1) * DLT[k][s];
            }
        } else {
            if (k == D - 1) break;
            k++;
        }
    }
    return g_hits;
}

static long long enum_hit_send(double R2) {
    const int do_hist = 0;
    int z[D], lo[D], hi[D];
    int64_t ST[D + 1][NST];
    
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;

    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    for (int s = 0; s < NST; s++) ST[k][s] = (int64_t)(lo[k] - 1) * DLT[k][s];
    if (do_hist) memset(g_hist, 0, sizeof g_hist);
    g_nodes = g_hits = 0;

    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            g_nodes++;
            if (do_hist) g_hist[k]++;
            for (int s = 0; s < NST; s++) ST[k][s] += DLT[k][s];
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;          /* boundary rounding only */
            if (k == 0) {
                g_hits++;
                hit_send(ST[0]);
            } else {
                rem[k - 1] = r;
                double zz = (double)z[k];
                for (i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
                hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
                z[k] = lo[k] - 1;
                for (int s = 0; s < NST; s++) ST[k][s] = ST[k + 1][s] + (int64_t)(lo[k] - 1) * DLT[k][s];
            }
        } else {
            if (k == D - 1) break;
            k++;
        }
    }
    return g_hits;
}


/* Optimised realistic kernel: generic DFS for levels >= 1; the last level is a tight
 * register loop over its exact interval [lo0, hi0] (every z0 in it is a hit, exactly as
 * in enumerate()), carrying the 7 ambient state words in registers. */
#define FAST_BODY(HITSTMT)                                                                  \
    int z[D], lo[D], hi[D];                                                                 \
    double C[D][D], rem[D];                                                                 \
    int64_t ST[D][NST];                                                                     \
    int k = D - 1, i;                                                                       \
    double rr;                                                                              \
    long long hits = 0;                                                                     \
    for (i = 0; i < D; i++) C[k][i] = CST[i];                                               \
    rem[k] = R2;                                                                            \
    rr = sqrt(rem[k] / BN[k]);                                                              \
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);                                              \
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);                                             \
    z[k] = lo[k] - 1;                                                                       \
    for (int q = 0; q < NST; q++) ST[k][q] = (int64_t)(lo[k] - 1) * DLT[k][q];              \
    const int64_t d0 = DLT[0][0], d1 = DLT[0][1], d2 = DLT[0][2], d3 = DLT[0][3],           \
                  d4 = DLT[0][4], d5 = DLT[0][5], d6 = DLT[0][6];                           \
    for (;;) {                                                                              \
        if (z[k] < hi[k]) {                                                                 \
            z[k]++;                                                                         \
            for (int q = 0; q < NST; q++) ST[k][q] += DLT[k][q];                            \
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
                for (int q = 0; q < NST; q++)                                               \
                    ST[k][q] = ST[k + 1][q] + (int64_t)(lo[k] - 1) * DLT[k][q];             \
            }                                                                               \
        } else {                                                                            \
            if (k == D - 1) break;                                                          \
            k++;                                                                            \
        }                                                                                   \
    }                                                                                       \
    return hits;

#define HIT_REAL_REGS {                                                                     \
        int64_t fa = s2 >> 63, fb = s3 >> 63;                                               \
        int64_t sa = fa ? CSA - s0 : s0;                                                    \
        int64_t sb = fb ? CSB - s1 : s1;                                                    \
        const slabdesc_t *dd = &slabtab[(uint64_t)(sa + sb) & smask];                       \
        int64_t t2 = (s4 ^ fa) - fa, t3 = (s5 ^ fa) - fa, t4 = (s6 ^ fa) - fa;              \
        uint64_t cell = (uint64_t)(dd->off + t2 + dd->w0 * t3 + dd->w01 * t4);              \
        bucket8[cell & b8mask] += (unsigned char)logp; }
#define HIT_SEND_REGS {                                                                     \
        int64_t fa = s2 >> 63, fb = s3 >> 63;                                               \
        int64_t sa = fa ? CSA - s0 : s0;                                                    \
        int64_t sb = fb ? CSB - s1 : s1;                                                    \
        const slabdesc_t *dd = &slabtab[(uint64_t)(sa + sb) & smask];                       \
        int64_t t2 = (s4 ^ fa) - fa, t3 = (s5 ^ fa) - fa, t4 = (s6 ^ fa) - fa;              \
        msg_t *mm = &outbuf[outn++ & (outcap - 1)];                                         \
        mm->dst = (uint64_t)(sa + sb);                                                      \
        mm->cell = (uint32_t)(dd->off + t2 + dd->w0 * t3 + dd->w01 * t4); mm->lp = logp; }

static long long enum_fast_walk(double R2) { FAST_BODY((void)(s0 ^ s1 ^ s2 ^ s3 ^ s4 ^ s5 ^ s6)) }
static long long enum_fast_real(double R2) { FAST_BODY(HIT_REAL_REGS) }
static long long enum_fast_send(double R2) { FAST_BODY(HIT_SEND_REGS) }

/* collect the points (not timed) */
static int (*pts)[D];
static long long npts, ptcap;

static long long enumerate_collect(double R2) {
    int z[D], lo[D], hi[D];
    double C[D][D], rem[D];
    int k = D - 1, i;
    double rr;
    npts = 0;
    for (i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    rr = sqrt(rem[k] / BN[k]);
    lo[k] = (int)ceil(C[k][k] - rr - BND_EPS);
    hi[k] = (int)floor(C[k][k] + rr + BND_EPS);
    z[k] = lo[k] - 1;
    g_nodes = g_hits = 0;
    for (;;) {
        if (z[k] < hi[k]) {
            z[k]++;
            g_nodes++;
            double t = C[k][k] - (double)z[k];
            double r = rem[k] - t * t * BN[k];
            if (r < 0.0) r = 0.0;
            if (k == 0) {
                g_hits++;
                if (npts < ptcap) memcpy(pts[npts++], z, sizeof z);
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

/* --------------------------------------------------- direct-form helpers */
/* norm^2_Q of (sum z_i B_i - centre), evaluated with the full quadratic form
 * (independent of the GSO machinery used by the enumeration). */
static double direct_norm2(const int *z) {
    double x[D], s = 0.0;
    int a, b, i;
    for (a = 0; a < D; a++) {
        double t = -LAT_CENTER[a];
        for (i = 0; i < D; i++) t += (double)LAT_BRED[i][a] * (double)z[i];
        x[a] = t;
    }
    for (a = 0; a < D; a++) {
        double t = 0.0;
        for (b = 0; b < D; b++) t += LAT_Q[a][b] * x[b];
        s += t * x[a];
    }
    return s;
}

/* brute force: scan z in [-M,M]^d, direct form (no pruning).  Fills out[] if
 * non-NULL, returns the count, sets *nboundary to the number of hits on the
 * box boundary (must be 0 for the box to have been large enough). */
static long long brute(double R2, int M, int (*out)[D], long long cap,
                       long long *nboundary) {
    int z[D];
    long long cnt = 0;
    int i, a;
    *nboundary = 0;
    for (i = 0; i < D; i++) z[i] = -M;
    for (;;) {
        double n2 = direct_norm2(z);
        if (n2 <= R2) {
            cnt++;
            for (a = 0; a < D; a++)
                if (z[a] == M || z[a] == -M) { (*nboundary)++; break; }
            if (out && cnt <= cap) memcpy(out[cnt - 1], z, sizeof z);
        }
        /* odometer */
        for (i = 0; i < D; i++) {
            if (z[i] < M) { z[i]++; break; }
            z[i] = -M;
        }
        if (i == D) break;
    }
    return cnt;
}

static int cmp_pts(const void *pa, const void *pb) {
    const int *a = pa, *b = pb;
    for (int i = 0; i < D; i++)
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    return 0;
}

/* ------------------------------------------------------------------ main */
static void alloc_bucket(void) {
    size_t nb = (size_t)1 << 22;             /* 4M * 2 B = 8 MB */
    bucket = calloc(nb, sizeof(unsigned short));
    if (!bucket) { perror("calloc"); exit(1); }
    bmask = (uint32_t)(nb - 1);
}

int main(int argc, char **argv) {
    const char *mode = argc > 1 ? argv[1] : "bench";
    int reps = argc > 2 ? atoi(argv[2]) : 5;
    double R2 = LAT_R2;                    /* any arg "R2=<value>" overrides */
    for (int i = 2; i < argc; i++)
        if (!strncmp(argv[i], "R2=", 3)) R2 = atof(argv[i] + 3);
    setup();
    alloc_bucket();

    if (!strcmp(mode, "bench")) {
        long long hits;
        uint64_t t0, t1;
        double ns_enum = 0, ns_full = 0, ns_hit = 0, ns_hit3 = 0, ns_full3 = 0;
        /* warm-up */
        hits = enumerate(R2, HITMODE_HIT9, 0);
        printf("hits=%lld  nodes=%lld  nodes/pt=%.3f  internal/pt=%.3f\n",
               hits, g_nodes, (double)g_nodes / hits,
               (double)(g_nodes - hits) / hits);
        for (int r = 0; r < reps; r++) {
            t0 = now_ns();
            enumerate(R2, HITMODE_NONE, 0);
            t1 = now_ns();
            ns_enum += (double)(t1 - t0) / hits;
            memset(bucket, 0, ((size_t)1 << 22) * sizeof(unsigned short));
            t0 = now_ns();
            enumerate(R2, HITMODE_HIT9, 0);
            t1 = now_ns();
            ns_full += (double)(t1 - t0) / hits;
            memset(bucket, 0, ((size_t)1 << 22) * sizeof(unsigned short));
            t0 = now_ns();
            enumerate(R2, HITMODE_HIT3, 0);
            t1 = now_ns();
            ns_full3 += (double)(t1 - t0) / hits;
        }
        ns_enum /= reps; ns_full /= reps; ns_full3 /= reps;
        /* hit-only pass over collected points */
        ptcap = hits + (hits >> 4) + 16;
        pts = malloc((size_t)ptcap * sizeof *pts);
        if (!pts) { perror("malloc"); exit(1); }
        long long tot = enumerate_collect(R2);
        if (tot != hits) printf("WARN collect hits %lld != %lld\n", tot, hits);
        for (int r = 0; r < reps; r++) {
            memset(bucket, 0, ((size_t)1 << 22) * sizeof(unsigned short));
            t0 = now_ns();
            for (long long q = 0; q < npts; q++) hit9(pts[q]);
            t1 = now_ns();
            ns_hit += (double)(t1 - t0) / npts;
            memset(bucket, 0, ((size_t)1 << 22) * sizeof(unsigned short));
            t0 = now_ns();
            for (long long q = 0; q < npts; q++) hit3(pts[q]);
            t1 = now_ns();
            ns_hit3 += (double)(t1 - t0) / npts;
        }
        ns_hit /= reps; ns_hit3 /= reps;
        printf("points=%lld reps=%d\n", hits, reps);
        printf("enum-only      %.3f ns/pt  (%.1f eff.instr @3.5GHz)\n",
               ns_enum, ns_enum * 3.5);
        printf("enum+hit9      %.3f ns/pt  (%.1f eff.instr @3.5GHz)\n",
               ns_full, ns_full * 3.5);
        printf("enum+hit3      %.3f ns/pt  (%.1f eff.instr @3.5GHz)\n",
               ns_full3, ns_full3 * 3.5);
        printf("hit9-only      %.3f ns/pt  (%.1f eff.instr @3.5GHz)\n",
               ns_hit, ns_hit * 3.5);
        printf("hit3-only      %.3f ns/pt  (%.1f eff.instr @3.5GHz)\n",
               ns_hit3, ns_hit3 * 3.5);
        return 0;
    }
    if (!strcmp(mode, "stcheck")) {
        setup_real();
        long long h = enumerate(R2, HITMODE_REAL, 1);
        printf("stcheck hits=%lld state mismatches=%lld\n", h, g_stbad);
        return g_stbad != 0;
    }
    if (!strcmp(mode, "icount")) {
        int hm = argc > 3 ? atoi(argv[3]) : 0;
        setup_real();
        long long h = 0;
        for (int r = 0; r < reps; r++) {
            outn = 0;
            if (hm == 13) h += enum_track_only(R2);
            else if (hm == 3) h += enum_hit_real(R2);
            else if (hm == 4) h += enum_hit_send(R2);
            else if (hm == 20) h += enum_fast_walk(R2);
            else if (hm == 23) h += enum_fast_real(R2);
            else if (hm == 24) h += enum_fast_send(R2);
            else h += enumerate(R2, hm, 0);
            if ((hm == HITMODE_SEND || hm == 24) && argc > 4) apply_buf(outbuf, outn < outcap ? outn : outcap);
        }
        unsigned long long s = 1469598103934665603ULL;
        for (size_t i = 0; i < ((size_t)1 << 24); i++) s = (s ^ bucket8[i]) * 1099511628211ULL;
        printf("reps=%d hits=%lld digest=%016llx\n", reps, h, s);
        return 0;
    }
    if (!strcmp(mode, "hist")) {
        long long hits = enumerate(R2, HITMODE_NONE, 1);
        printf("hits=%lld nodes=%lld\n", hits, g_nodes);
        for (int i = D - 1; i >= 0; i--)
            printf("  level %2d: %8lld node visits\n", i, g_hist[i]);
        return 0;
    }
    if (!strcmp(mode, "count")) {
        long long hits = enumerate(R2, HITMODE_NONE, 0);
        printf("%lld\n", hits);
        printf("nodes=%lld nodes/pt=%.3f\n", g_nodes, (double)g_nodes / hits);
        return 0;
    }
    if (!strcmp(mode, "check") || !strcmp(mode, "dump") || !strcmp(mode, "verify")) {
        ptcap = (long long)(LAT_EXPECTED * 1.1) + 1000;
        pts = malloc((size_t)ptcap * sizeof *pts);
        if (!pts) { perror("malloc"); exit(1); }
        long long hits = enumerate_collect(R2);
        printf("hits=%lld expected(gh)=%.1f ratio=%.4f\n", hits, LAT_EXPECTED,
               (double)hits / LAT_EXPECTED);
        printf("nodes=%lld nodes/pt=%.3f internal/pt=%.3f\n", g_nodes,
               (double)g_nodes / hits, (double)(g_nodes - hits) / hits);
        /* containment via the direct form */
        double maxr = 0.0;
        long long bad = 0;
        for (long long q = 0; q < npts; q++) {
            double n2 = direct_norm2(pts[q]);
            if (n2 > maxr) maxr = n2;
            if (n2 > R2 * (1.0 + 1e-12)) bad++;
        }
        printf("direct-form: max norm^2 / R2 = %.12f   outside=%lld\n",
               maxr / R2, bad);
        /* duplicates */
        qsort(pts, (size_t)npts, sizeof *pts, cmp_pts);
        long long dup = 0;
        for (long long q = 1; q < npts; q++)
            if (!memcmp(pts[q], pts[q - 1], sizeof pts[q])) dup++;
        printf("duplicates=%lld\n", dup);
        if (!strcmp(mode, "dump")) {
            const char *fn = argc > 2 ? argv[2] : "points.txt";
            int want = argc > 3 ? atoi(argv[3]) : 1000;
            FILE *f = fopen(fn, "w");
            if (!f) { perror(fn); exit(1); }
            for (int q = 0; q < want; q++) {
                long long idx = (long long)q * npts / want;   /* sample evenly */
                for (int i = 0; i < D; i++)
                    fprintf(f, "%d%c", pts[idx][i], i + 1 < D ? ' ' : '\n');
            }
            fclose(f);
            printf("dumped %d points to %s\n", want, fn);
        }
        if (!strcmp(mode, "verify")) {
            int M = argc > 2 ? atoi(argv[2]) : 8;
            int (*bp)[D] = malloc((size_t)(hits + 1000) * sizeof *bp);
            long long nb = 0;
            long long bc = brute(R2, M, bp, hits + 1000, &nb);
            printf("brute box M=%d: count=%lld boundary_hits=%lld  %s\n", M, bc,
                   nb, bc == hits ? "COUNT MATCH" : "COUNT MISMATCH");
            qsort(bp, (size_t)bc, sizeof *bp, cmp_pts);
            printf("sets: %s\n", (bc == hits && !memcmp(bp, pts, (size_t)hits * sizeof *bp))
                   ? "IDENTICAL (sorted)" : "DIFFER");
            /* box adequacy: no hit may sit on the box boundary */
            return nb == 0 && bc == hits ? 0 : 1;
        }
        return 0;
    }
    fprintf(stderr, "unknown mode %s\n", mode);
    return 2;
}
