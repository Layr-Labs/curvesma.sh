/* slice12: ARCH-E a-slice sieve kernel with PACKED ROWS (stnfs-aslice-v12).
 *
 * One a-slice is the 6-dimensional half-ellipsoid S = {b in Z^6 : Q_b(b) <= R^2, b5 >= 1}.
 *
 * Layout (one cell per point of S, no bounding-box bytes):
 *   - LINES k = b5 = 1..K; SEGMENTS of G = 2^gs consecutive lines, segment index (k-1) >> gs.
 *   - A ROW is the set of b0 with (b0, b1..b4, k) in S for a fixed suffix (b1..b4, k).  Rows are packed one
 *     after the other inside their segment; cells are segment-local indices (int64, < 2^62 / 4 per segment).
 *   - Per line, a dense METADATA array over the 4-dimensional suffix box prod_{i=1..4} [-r_i, r_i] (the
 *     projection of the whole slice ellipsoid, the same for every line):  meta[k][idx4] = (cell of b0 = 0 in
 *     that row) = row start - lo, so a point's cell is meta[k][idx4] + b0.  idx4 = sum_{i=1..4} (b_i + r_i) T_i is
 *     affine in the walk coordinates; the fold b -> -b maps idx4 to (L4 - 1) - idx4 and b0 to -b0.
 *   - Every cell holds two uint16 sums (side f and side g, separate arrays), in units of 1/16 bit.  A sum never
 *     overflows: only ideals containing the point add, so the sum is <= 16 log2 |N| <= 16 * 1302 < 65535.
 *   Rows are built with the radius R^2 (1 + 2^-30) + 2^-20, so every point the walk emits (Q <= R^2 up to the
 *   walk's 1e-9 coordinate guard) lies in its row; the check mode compares every cell with brute force.
 *
 * Paths (cosets recentred before enumeration, so walk coordinates are small for every p < 2^63):
 *   LINE path (p < 2^20): per line k and each of the two cosets (ideal at +k; ideal at -k folded = conjugate's coset
 *     at +k), a 5-dimensional Schnorr-Euchner walk of the coset of L5(p,s) inside the line ellipsoid; direct update
 *     arr[meta[k][idx4] + b0] += logp.
 *   RECORD path (p >= 2^20): one 6-dimensional walk of the coset of L(p,s) over the full ellipsoid (both signs of b5)
 *     with affine states (idx4, b0, b5); a hit with b5 < 0 is folded, b5 = 0 is skipped, and the record
 *     (segment-local cell uint64, global ideal id uint64) is appended to its segment's RETAINED record store
 *     (64-bit counters; a full store is grown, never wrapped).  The APPLY pass of a prime batch (all ideals of one
 *     side with the same scaled log) adds the batch log to the cells of the records appended since the last apply.
 *   SCAN: per row, per threshold block (thresholds come from the certified norm-bound blocks, nbound12.c), a cell is
 *     a (preliminary) survivor iff f >= Tf and g >= Tg; survivors are listed and their f sum set to the sentinel
 *     0xFFFF (unreachable by a true sum).
 *   RECOVER: one pass over each segment's retained records; a record whose cell holds the sentinel is copied to the
 *     survivor list, giving every survivor ALL its sieved ideals >= 2^20 with their global ids.
 *
 * Arithmetic: modular setups use the Barrett reduction bar_red (validated by mode `barrett`); there is no % in any
 *   sieving path.  Recentring: q_i = llround(t_i) with t = W o in binary64; ANY integer q gives a valid coset
 *   representative o - sum q_i B_i, which is computed exactly in wrapping uint64 arithmetic (the exact result is small,
 *   so the wrapped value equals it).  Floating error only moves the representative by a few basis vectors; mode
 *   `recentre` measures the deviation from exact rounding (computed with int128) up to p ~ 2^62.
 *
 * Modes:
 *   check log2pts gs                   line + record paths (both sides, several primes), apply, scan, recover vs brute force
 *   barrett n                          n random Barrett reductions per bit length 3..63 vs the C remainder
 *   recentre n                         deviation of the floating rounding from exact rounding, p up to 2^62
 *   bench <what> log2pts gs p s c reps one timed function (callgrind --toggle-collect=<function>):
 *        what = line (line_path), rec (record_path), apply (apply_batch), setup6 (setup6), setupl (setup_line),
 *               scan (scan_segment), recover (recover_segment), build (build_layout), zero (zero_cells)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "slice_mb.h"
#define D 6
#define BND_EPS 1e-9
#define NOINLINE __attribute__((noinline))
#define POISON INT64_MIN
#define SENTINEL 0xFFFFu
#define AMAX 4096ULL

/* ---------------- Barrett reduction, valid for every odd p in [3, 2^63) that is not a power of two ---------------- */
typedef struct { uint64_t p, mu; unsigned k; } bar_t;
static bar_t bar_init(uint64_t p) {            /* factor-base construction (once per prime) */
    bar_t b; b.p = p; b.k = 64 - (unsigned)__builtin_clzll(p);
    b.mu = (uint64_t)((((__uint128_t)1) << (2 * b.k)) / p);   /* < 2^(k+1) <= 2^64 since 2^(k-1) < p < 2^k */
    return b;
}
/* x < 2^(2k): q = floor(floor(x / 2^(k-1)) mu / 2^(k+1)) satisfies q_true - 2 <= q <= q_true (HAC 14.42) */
static inline uint64_t bar_red(const bar_t *b, __uint128_t x) {
    uint64_t xs = (uint64_t)(x >> (b->k - 1));
    uint64_t q = (uint64_t)(((__uint128_t)xs * b->mu) >> (b->k + 1));
    __uint128_t r = x - (__uint128_t)q * b->p;
    if (r >= b->p) r -= b->p;
    if (r >= b->p) r -= b->p;
    return (uint64_t)r;
}

/* ---------------- lattice tools (dimension n <= 6), integer bases, Gram-Schmidt in a metric ---------------- */
typedef struct { int n; long long B[D][D]; double MU[D][D], BN[D]; } lat_t;
static void gso(lat_t *L, const double *Mx) {
    int n = L->n; double Bst[D][D];
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) Bst[i][k] = (double)L->B[i][k];
        for (int j = 0; j < i; j++) {
            double num = 0;
            for (int k = 0; k < n; k++) { double r = 0; for (int l = 0; l < n; l++) r += Mx[k * n + l] * Bst[j][l]; num += (double)L->B[i][k] * r; }
            L->MU[i][j] = num / L->BN[j];
            for (int k = 0; k < n; k++) Bst[i][k] -= L->MU[i][j] * Bst[j][k];
        }
        double nn = 0; for (int k = 0; k < n; k++) { double r = 0; for (int l = 0; l < n; l++) r += Mx[k * n + l] * Bst[i][l]; nn += Bst[i][k] * r; }
        L->BN[i] = nn;
    }
}
static void lll(lat_t *L, const double *Mx) {
    int n = L->n; gso(L, Mx); int k = 1;
    while (k < n) {
        for (int j = k - 1; j >= 0; j--) { long long q = llround(L->MU[k][j]); if (q) { for (int i = 0; i < n; i++) L->B[k][i] -= q * L->B[j][i]; gso(L, Mx); } }
        if (L->BN[k] >= (0.99 - L->MU[k][k - 1] * L->MU[k][k - 1]) * L->BN[k - 1]) k++;
        else { long long t[D]; memcpy(t, L->B[k], sizeof t); memcpy(L->B[k], L->B[k - 1], sizeof t); memcpy(L->B[k - 1], t, sizeof t); gso(L, Mx); k = k > 1 ? k - 1 : 1; }
    }
}
static void hyperplane_lattice(lat_t *L, int n, uint64_t p, uint64_t s) {   /* {x in Z^n : sum x_i s^i = 0 mod p} */
    L->n = n; memset(L->B, 0, sizeof L->B); L->B[0][0] = (long long)p; uint64_t si = 1;
    for (int i = 1; i < n; i++) { si = (uint64_t)((__uint128_t)si * s % p); L->B[i][0] = (long long)((p - si) % p); L->B[i][i] = 1; }
}
static void basis_inverse(const lat_t *L, double W[D][D]) {   /* W[i][j]: coordinates t_i = sum_j W[i][j] v_j */
    int n = L->n; double A[D][2 * D];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { A[i][j] = (double)L->B[j][i]; A[i][n + j] = i == j; }
    for (int c = 0; c < n; c++) { int pr = c; for (int r = c + 1; r < n; r++) if (fabs(A[r][c]) > fabs(A[pr][c])) pr = r;
        for (int k = 0; k < 2 * n; k++) { double t = A[c][k]; A[c][k] = A[pr][k]; A[pr][k] = t; }
        double d = A[c][c]; for (int k = 0; k < 2 * n; k++) A[c][k] /= d;
        for (int r = 0; r < n; r++) if (r != c) { double f = A[r][c]; for (int k = 0; k < 2 * n; k++) A[r][k] -= f * A[c][k]; } }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) W[i][j] = A[i][n + j];
}
/* origin (c, 0, ..., 0): q_i = round(W[i][0] c); out = o - sum q_i B_i exactly (wrapping uint64; the true value is small) */
static inline void recentre(const lat_t *L, const double W[D][D], uint64_t c, long long *out) {
    int n = L->n; long long q[D]; double cd = (double)c;
    for (int i = 0; i < n; i++) q[i] = llround(W[i][0] * cd);
    for (int k = 0; k < n; k++) {
        uint64_t v = k == 0 ? c : 0;
        for (int i = 0; i < n; i++) v -= (uint64_t)q[i] * (uint64_t)L->B[i][k];
        out[k] = (long long)v;
    }
}
static inline void gs_centre(const lat_t *L, const double W[D][D], const double *ctr, double *cst) {
    int n = L->n; double t[D];
    for (int i = 0; i < n; i++) { t[i] = 0; for (int j = 0; j < n; j++) t[i] += W[i][j] * ctr[j]; }
    for (int k = 0; k < n; k++) { double v = t[k]; for (int i = k + 1; i < n; i++) v += t[i] * L->MU[i][k]; cst[k] = v; }
}

/* ---------------- region and packed layout ---------------- */
static double M6[D][D], M5[5][5], Minv5[5][5], mvec[5], shift5[5], dconst, R2, R2B;
static long long r6[D], T4[5], L4, L4m1, K;   /* box radii; strides of b1..b4 (T4[1..4]); suffix box size; lines */
static int gs; static long long nseg;
typedef struct { long long start; int lo, len; long long idx4; long long k; } row_t;
typedef struct {
    uint16_t *f, *g; long long ncell;
    row_t *rows; long long nrows, rowcap;
    uint64_t *rcell, *rid; uint64_t rn, rcap, applied;
} seg_t;
static seg_t *SEG; static long long **meta;   /* meta[k] for k = 1..K */
static void inv_n(const double *A0, int n, double *out) {
    double A[D][2 * D]; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { A[i][j] = A0[i * n + j]; A[i][n + j] = i == j; }
    for (int c = 0; c < n; c++) { double d = A[c][c]; for (int k = 0; k < 2 * n; k++) A[c][k] /= d; for (int r = 0; r < n; r++) if (r != c) { double f = A[r][c]; for (int k = 0; k < 2 * n; k++) A[r][k] -= f * A[c][k]; } }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) out[i * n + j] = A[i][n + j];
}
static double det_n(const double *A0, int n) {
    double A[D][D]; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = A0[i * n + j];
    double d = 1; for (int c = 0; c < n; c++) { d *= A[c][c]; for (int r = c + 1; r < n; r++) { double f = A[r][c] / A[c][c]; for (int k = c; k < n; k++) A[r][k] -= f * A[c][k]; } } return d;
}
static void setup_region(int log2pts, int gshift) {
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) M6[i][j] = MB0[i][j];
    for (int i = 0; i < 5; i++) { for (int j = 0; j < 5; j++) M5[i][j] = M6[i][j]; mvec[i] = M6[i][5]; }
    inv_n((double *)M5, 5, (double *)Minv5);
    for (int i = 0; i < 5; i++) { shift5[i] = 0; for (int j = 0; j < 5; j++) shift5[i] += Minv5[i][j] * mvec[j]; }
    dconst = M6[5][5]; for (int i = 0; i < 5; i++) dconst -= mvec[i] * shift5[i];
    double vol = ldexp(1.0, log2pts + 1), det6 = det_n((double *)M6, 6);     /* full ellipsoid ~ 2^(log2pts+1) points */
    R2 = pow(vol * sqrt(det6) * 6 / (M_PI * M_PI * M_PI), 1.0 / 3);
    R2B = R2 * (1 + ldexp(1.0, -30)) + ldexp(1.0, -20);                      /* row-building radius (margin) */
    double Minv6[D][D]; inv_n((double *)M6, 6, (double *)Minv6);
    for (int i = 0; i < D; i++) r6[i] = (long long)floor(sqrt(R2B * Minv6[i][i]));
    K = r6[5]; gs = gshift; nseg = ((K - 1) >> gs) + 1;
    L4 = 1; for (int i = 1; i <= 4; i++) { T4[i] = L4; L4 *= 2 * r6[i] + 1; } L4m1 = L4 - 1;
}
static void push_row(seg_t *sg, row_t r) {
    if (sg->nrows == sg->rowcap) { sg->rowcap = sg->rowcap ? 2 * sg->rowcap : 1024; sg->rows = realloc(sg->rows, sizeof(row_t) * sg->rowcap); }
    sg->rows[sg->nrows++] = r;
}
/* rows of every line: the b0 interval of each suffix of the 4-dim box (production: once per slice, charged per cell) */
static NOINLINE void build_layout(void) {
    for (long long k = 1; k <= K; k++) {
        seg_t *sg = &SEG[(k - 1) >> gs]; long long *ml = meta[k];
        long long b[D]; b[5] = k; long long idx = 0;
        for (b[4] = -r6[4]; b[4] <= r6[4]; b[4]++) for (b[3] = -r6[3]; b[3] <= r6[3]; b[3]++)
        for (b[2] = -r6[2]; b[2] <= r6[2]; b[2]++) {
            /* partial sums over b2..b5, then the b1 loop */
            double beta0 = 0, gam0 = 0, beta1 = 0;
            for (int j = 2; j <= 5; j++) { beta0 += M6[0][j] * b[j]; beta1 += M6[1][j] * b[j]; for (int i = 2; i <= 5; i++) gam0 += M6[i][j] * b[i] * b[j]; }
            for (b[1] = -r6[1]; b[1] <= r6[1]; b[1]++, idx++) {
                double b1 = (double)b[1];
                double beta = beta0 + M6[0][1] * b1, gam = gam0 + 2 * b1 * beta1 + M6[1][1] * b1 * b1;
                double disc = beta * beta - M6[0][0] * (gam - R2B);
                if (disc < 0) { ml[idx] = POISON; continue; }
                double w = sqrt(disc) / M6[0][0], c = -beta / M6[0][0];
                long long lo = (long long)ceil(c - w), hi = (long long)floor(c + w);
                if (lo > hi) { ml[idx] = POISON; continue; }
                ml[idx] = sg->ncell - lo;
                row_t r = {sg->ncell, (int)lo, (int)(hi - lo + 1), idx, k}; push_row(sg, r);
                sg->ncell += hi - lo + 1;
            }
        }
    }
}
/* zero both sum arrays with 64-bit general-register stores (no library memset, no vector registers) */
static NOINLINE void zero_cells(seg_t *sg) {
    for (int side = 0; side < 2; side++) {
        uint16_t *a = side ? sg->g : sg->f; long long n = sg->ncell + 1, i = 0;
        for (; i < n && ((uintptr_t)(a + i) & 7); i++) a[i] = 0;
        uint64_t *w = (uint64_t *)(a + i); long long nw = (n - i) >> 2, j = 0;   /* built with -fno-tree-loop-distribute-patterns */
        for (; j + 4 <= nw; j += 4) { w[j] = 0; w[j + 1] = 0; w[j + 2] = 0; w[j + 3] = 0; }
        for (; j < nw; j++) w[j] = 0;
        for (i += nw << 2; i < n; i++) a[i] = 0;
    }
}
static void alloc_slice(void) {
    SEG = calloc(nseg, sizeof(seg_t)); meta = calloc(K + 1, sizeof(long long *));
    for (long long k = 1; k <= K; k++) meta[k] = malloc(sizeof(long long) * L4);
    build_layout();
    for (long long s = 0; s < nseg; s++) {
        SEG[s].f = malloc(sizeof(uint16_t) * (SEG[s].ncell + 8)); SEG[s].g = malloc(sizeof(uint16_t) * (SEG[s].ncell + 8));
        zero_cells(&SEG[s]);
        SEG[s].rcap = 1 << 16; SEG[s].rcell = malloc(8 * SEG[s].rcap); SEG[s].rid = malloc(8 * SEG[s].rcap);
    }
}
static long long total_cells(void) { long long t = 0; for (long long s = 0; s < nseg; s++) t += SEG[s].ncell; return t; }
static NOINLINE void grow_records(seg_t *sg) {   /* production: spill the full buffer to the segment's record store */
    sg->rcap *= 2; sg->rcell = realloc(sg->rcell, 8 * sg->rcap); sg->rid = realloc(sg->rid, 8 * sg->rcap);
}

/* ---------------- Schnorr-Euchner walk with three affine 64-bit states ---------------- */
typedef struct { long long base; long long dl[D]; } aff_t;
#define WALK3(L, cst, rad2, A1, A2, A3, HITSTMT)                                                                  \
    do {                                                                                                          \
        int n = (L)->n; long long z[D], zlo[D], zhi[D]; double C[D][D], rem[D];                                   \
        long long SA[D], SB[D], SC[D];                                                                            \
        int k = n - 1; for (int i = 0; i < n; i++) C[k][i] = (cst)[i];                                            \
        rem[k] = (rad2); double rr = sqrt(rem[k] / (L)->BN[k]);                                                   \
        zlo[k] = (long long)ceil(C[k][k] - rr - BND_EPS); zhi[k] = (long long)floor(C[k][k] + rr + BND_EPS); z[k] = zlo[k] - 1; \
        SA[k] = (A1).base + (zlo[k] - 1) * (A1).dl[k]; SB[k] = (A2).base + (zlo[k] - 1) * (A2).dl[k];            \
        SC[k] = (A3).base + (zlo[k] - 1) * (A3).dl[k];                                                            \
        const long long da0 = (A1).dl[0], db0 = (A2).dl[0], dc0 = (A3).dl[0];                                     \
        for (;;) {                                                                                                \
            if (z[k] < zhi[k]) {                                                                                  \
                z[k]++; SA[k] += (A1).dl[k]; SB[k] += (A2).dl[k]; SC[k] += (A3).dl[k];                            \
                double tt = C[k][k] - (double)z[k]; double r_ = rem[k] - tt * tt * (L)->BN[k]; if (r_ < 0) r_ = 0;   \
                double zz = (double)z[k];                                                                         \
                if (k == 1) {                                                                                     \
                    double c0 = C[1][0] - zz * (L)->MU[1][0]; double r0 = sqrt(r_ / (L)->BN[0]);                  \
                    long long l0 = (long long)ceil(c0 - r0 - BND_EPS), h0 = (long long)floor(c0 + r0 + BND_EPS); \
                    long long a_ = SA[1] + l0 * da0, b_ = SB[1] + l0 * db0, k_ = SC[1] + l0 * dc0;               \
                    for (long long j = l0; j <= h0; j++) { HITSTMT; a_ += da0; b_ += db0; k_ += dc0; }            \
                    if (h0 >= l0) hits += h0 - l0 + 1;                                                            \
                } else {                                                                                          \
                    rem[k - 1] = r_; for (int i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * (L)->MU[k][i];      \
                    k--; rr = sqrt(rem[k] / (L)->BN[k]);                                                          \
                    zlo[k] = (long long)ceil(C[k][k] - rr - BND_EPS); zhi[k] = (long long)floor(C[k][k] + rr + BND_EPS); z[k] = zlo[k] - 1; \
                    SA[k] = SA[k + 1] + (zlo[k] - 1) * (A1).dl[k]; SB[k] = SB[k + 1] + (zlo[k] - 1) * (A2).dl[k]; \
                    SC[k] = SC[k + 1] + (zlo[k] - 1) * (A3).dl[k];                                                \
                }                                                                                                 \
            } else { if (k == n - 1) break; k++; }                                                                \
        }                                                                                                         \
    } while (0)

/* ---------------- an ideal pair: lattices are factor-base data (built once per (p, s), charged there) ---------------- */
typedef struct {
    uint64_t p, s, s5, rinv, aoff; bar_t bar; uint64_t id; int side; uint16_t lp;
    lat_t L5, L6; double W5[D][D], W6[D][D];
} pair_t;
static void prepare_pair(pair_t *P, uint64_t p, uint64_t s, uint64_t rinv, uint64_t id, int side) {
    P->p = p; P->s = s; P->rinv = rinv; P->id = id; P->side = side; P->bar = bar_init(p);
    P->aoff = p * ((AMAX + p - 1) / p);              /* multiple of p >= AMAX, so a_i + aoff >= 0 */
    P->lp = (uint16_t)floor(16 * log2((double)p));   /* floor of 16 log2 p (the batch log) */
    uint64_t s5 = 1; for (int i = 0; i < 5; i++) s5 = (uint64_t)((__uint128_t)s5 * s % p); P->s5 = s5;
    hyperplane_lattice(&P->L5, 5, p, s); lll(&P->L5, (double *)M5); basis_inverse(&P->L5, P->W5);
    hyperplane_lattice(&P->L6, 6, p, s); lll(&P->L6, (double *)M6); basis_inverse(&P->L6, P->W6);
}

/* per-(slice, pair) setup of both paths: c = a(s) r^-1 mod p by Horner with Barrett steps.
   |a_i| <= AMAX = 2^12 (the a-ellipsoid radius is < 1600); a_i + aoff < 2^13 + p < 2^(2k) for p >= 2^7. */
static NOINLINE uint64_t setup6(const pair_t *P, const long long *a) {
    const bar_t *br = &P->bar; uint64_t v = 0;
    for (int i = 5; i >= 0; i--) {
        v = bar_red(br, (__uint128_t)v * P->s);
        v += bar_red(br, (__uint128_t)((uint64_t)a[i] + P->aoff)); if (v >= P->p) v -= P->p;
    }
    return bar_red(br, (__uint128_t)v * P->rinv);
}

/* LINE path: per (pair, line, coset) setup + 5-dim walk with direct updates */
typedef struct { aff_t A1, A2; double cst[D]; double rad2; } lset_t;
static NOINLINE int setup_line(const pair_t *P, uint64_t c, long long k, int sign, lset_t *o) {
    double rad2 = R2 - (double)k * k * dconst; if (rad2 <= 0) return 0;
    const bar_t *br = &P->bar; uint64_t p = P->p;
    uint64_t ks5 = bar_red(br, (__uint128_t)(uint64_t)k * P->s5);
    uint64_t ck;
    if (!sign) ck = c >= ks5 ? c - ks5 : c + p - ks5;                          /* ideal at +k: b'(s) = c - k s^5 */
    else { uint64_t t = c + ks5; if (t >= p) t -= p; ck = t ? p - t : 0; }    /* conjugate at +k: b'(s) = -c - k s^5 */
    long long oo[D]; recentre(&P->L5, P->W5, ck, oo);
    double ctr[5]; for (int i = 0; i < 5; i++) ctr[i] = -(double)k * shift5[i] - (double)oo[i];
    gs_centre(&P->L5, P->W5, ctr, o->cst);
    o->A1.base = 0; for (int i = 1; i <= 4; i++) o->A1.base += (oo[i] + r6[i]) * T4[i];
    o->A2.base = oo[0];
    for (int j = 0; j < 5; j++) { o->A1.dl[j] = 0; for (int i = 1; i <= 4; i++) o->A1.dl[j] += P->L5.B[j][i] * T4[i]; o->A2.dl[j] = P->L5.B[j][0]; }
    o->rad2 = rad2; return 1;
}
static NOINLINE long long line_walk(const pair_t *P, const lset_t *o, uint16_t *arr, const long long *ml) {
    long long hits = 0; const uint16_t lp = P->lp; aff_t A3 = {0, {0}};
    WALK3(&P->L5, o->cst, o->rad2, o->A1, o->A2, A3, { arr[ml[a_] + b_] += lp; (void)k_; });
    return hits;
}
static NOINLINE long long line_path(const pair_t *P, uint64_t c) {
    long long hits = 0; lset_t o;
    for (long long k = 1; k <= K; k++) {
        seg_t *sg = &SEG[(k - 1) >> gs]; uint16_t *arr = P->side ? sg->g : sg->f;
        for (int sign = 0; sign < 2; sign++) if (setup_line(P, c, k, sign, &o)) hits += line_walk(P, &o, arr, meta[k]);
    }
    return hits;
}

/* RECORD path: 6-dim walk over the full ellipsoid, fold, append to the segment's retained record store */
static NOINLINE long long record_path(const pair_t *P, uint64_t c) {
    long long hits = 0; long long oo[D]; recentre(&P->L6, P->W6, c, oo);
    double ctr[D], cst[D]; for (int i = 0; i < D; i++) ctr[i] = -(double)oo[i];
    gs_centre(&P->L6, P->W6, ctr, cst);
    aff_t A1, A2, A3;                                   /* idx4, b0, b5 */
    A1.base = 0; for (int i = 1; i <= 4; i++) A1.base += (oo[i] + r6[i]) * T4[i];
    A2.base = oo[0]; A3.base = oo[5];
    for (int j = 0; j < D; j++) { A1.dl[j] = 0; for (int i = 1; i <= 4; i++) A1.dl[j] += P->L6.B[j][i] * T4[i]; A2.dl[j] = P->L6.B[j][0]; A3.dl[j] = P->L6.B[j][5]; }
    const long long Lm = L4m1; const int g = gs; const uint64_t id = P->id; seg_t *const S = SEG; long long **const MT = meta;
    WALK3(&P->L6, cst, R2, A1, A2, A3, {
        long long m_ = k_ >> 63;
        long long kf = (k_ ^ m_) - m_;                            /* |b5| */
        if (kf) {                                                 /* b5 = 0 is outside the half-slice */
            long long i4 = ((Lm - a_) & m_) | (a_ & ~m_);         /* fold of the suffix index */
            long long b0 = (b_ ^ m_) - m_;                        /* fold of b0 */
            seg_t *sg = &S[(kf - 1) >> g]; uint64_t n_ = sg->rn;
            if (n_ == sg->rcap) grow_records(sg);
            sg->rcell[n_] = (uint64_t)(MT[kf][i4] + b0); sg->rid[n_] = id; sg->rn = n_ + 1;
        }
    });
    return hits;
}
/* APPLY: the records appended since the last apply all belong to one prime batch (same side, same scaled log) */
static NOINLINE uint64_t apply_batch(seg_t *sg, int side, uint16_t lp) {
    uint16_t *arr = side ? sg->g : sg->f; const uint64_t *rc = sg->rcell; uint64_t n = sg->rn;
    for (uint64_t i = sg->applied; i < n; i++) arr[rc[i]] += lp;
    uint64_t done = n - sg->applied; sg->applied = n; return done;
}

/* SCAN: thresholds per block (Tf, Tg); survivors listed and marked */
typedef struct { long long cell; long long seg; } surv_t;
static surv_t *SURV; static long long nsurv, survcap;
static NOINLINE long long scan_segment(long long si, const long long *blk_start, const uint16_t *tf, const uint16_t *tg, long long nblk) {
    seg_t *sg = &SEG[si]; uint16_t *f = sg->f; const uint16_t *g = sg->g; long long found = 0;
    for (long long b = 0; b < nblk; b++) {
        const uint16_t Tf = tf[b], Tg = tg[b]; const long long e = blk_start[b + 1];
        for (long long c = blk_start[b]; c < e; c++) {
            if ((f[c] >= Tf) & (g[c] >= Tg)) {
                if (nsurv == survcap) { survcap = survcap ? 2 * survcap : 1024; SURV = realloc(SURV, sizeof(surv_t) * survcap); }
                SURV[nsurv].cell = c; SURV[nsurv].seg = si; nsurv++; f[c] = SENTINEL; found++;
            }
        }
    }
    return found;
}
/* RECOVER: one pass over the segment's retained records */
typedef struct { uint64_t cell, id; } rec_t;
static rec_t *RECOV; static uint64_t nrecov, recovcap;
static NOINLINE uint64_t recover_segment(seg_t *sg) {
    const uint16_t *f = sg->f; const uint64_t *rc = sg->rcell, *ri = sg->rid; uint64_t n = sg->rn, found = 0;
    for (uint64_t i = 0; i < n; i++) {
        if (f[rc[i]] == SENTINEL) {
            if (nrecov == recovcap) { recovcap = recovcap ? 2 * recovcap : 1024; RECOV = realloc(RECOV, sizeof(rec_t) * recovcap); }
            RECOV[nrecov].cell = rc[i]; RECOV[nrecov].id = ri[i]; nrecov++; found++;
        }
    }
    return found;
}

#include "slice12_check.c"   /* brute-force reference, check/barrett/recentre/bench modes and main: one translation unit */
