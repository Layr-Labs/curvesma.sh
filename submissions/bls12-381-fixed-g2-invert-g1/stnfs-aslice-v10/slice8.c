/* slice8: ARCH-E a-slice sieve kernel, production-width version (stnfs-aslice-v8).
 *
 * One a-slice is the 6-dimensional half-ellipsoid S = {b in Z^6 : Q_b(b) <= R^2, b5 >= 1}.  It is organised in LINES
 * (b5 = k, k = 1..K) and SEGMENTS (ranges of consecutive lines).  Every line has the same 5-dimensional bounding box
 * BOX5 = prod_{i<5} [-r_i, r_i] (the projection of the full ellipsoid; symmetric), so a point (b', k) has the address
 *     pos = (k - 1) * L5 + sum_{i<5} (b_i + r_i) * S_i,   L5 = |BOX5|,   0 <= pos < K * L5,
 * and a segment of G lines owns the byte array pos in [(k0-1) L5, (k0-1+G) L5): at most 2^62 bytes in production,
 * 64-bit addressable; the segment index of a point is (k - 1) / G.  The fold b -> -b of the conjugate ideal keeps
 * the line (b5 -> -b5 -> +k) and maps the local part to K5 - local, K5 = sum_{i<5} 2 r_i S_i.
 *
 * Two sieving paths (both with exactly recentred cosets, so every enumeration coordinate is small):
 * Primes are unsigned 64-bit (every p < 2^64 is representable).  Record positions are SEGMENT-LOCAL
 * ((k-1) mod G) L5 + local < G L5 <= 2^62 and records are a structure of two arrays (uint64 position, uint32 prime
 * index: 12 bytes per record, two stores), the segment index is (k-1)/G (a shift in production, G a power of two) and
 * the apply pass indexes the segment's own bucket array with the local position; no global offset is ever formed.
 *   LINE path (p < P_MID = 2^20): for each line k and each coset of the pair (the ideal at +k, and the ideal at -k which is
 *     folded onto +k), the 5-dimensional coset {b' : b'(s) = c_k mod p}, c_k = c - k s^5, of the lattice
 *     L5(p,s) = {b' : sum_{i<5} b_i s^i = 0 mod p} (reduced once per (p,s) in the metric of the line ellipsoid) is
 *     enumerated inside the line's 5-dim ellipsoid (centre -k M5^{-1} m, radius^2 R^2 - k^2 d) by a Schnorr-Euchner
 *     walk carrying the 64-bit local address; bytes are updated directly.  Setup per (pair, line, coset): c_k, the
 *     origin (c_k,0,0,0,0) recentred by exact Babai rounding, the walk centre.
 *   RECORD path (p >= P_MID = 2^20): one 6-dimensional coset enumeration per pair over the FULL ellipsoid (both signs of b5),
 *     carrying the line index k and the local address as affine state; each hit is folded if k < 0 and appended as a
 *     record {segment-local pos (8 bytes), prime index (4 bytes)} to the buffers of its segment; an APPLY pass per segment adds the
 *     prime's log to bucket[pos].  Survivors read their segment's records to recover their large sieved primes.
 *
 * Correctness: `check` runs both paths for a set of pairs below and above P_MID and compares the complete bucket
 * array of the half-slice with a brute-force reference (every point of BOX5 x [1..K] with b'(s) + k s^5 = c mod p
 * for the ideal, or = -c for its conjugate, tested against Q_b directly): bit-identical digests.
 * All coordinates and addresses are 64-bit integers; coset origins are recentred before any floating-point use, so
 * the Schnorr-Euchner centre coordinates are O(1) regardless of p.
 *
 * Modes:
 *   slice8 check  log2pts K p s c                 both paths vs brute force on the half-slice with K lines
 *   slice8 icount5 log2pts K p s c reps           line path: K lines x 2 cosets per rep (retired instructions)
 *   slice8 icount6 log2pts K p s c reps           record path: walk + records (no apply)
 *   slice8 icount6a log2pts K p s c reps          record path: walk + records + apply
 *   slice8 setup  p s trials                      per-(pair, line, coset) setup of the line path
 *   slice8 setup6 p s trials                      per-(slice, pair) setup of the record path (a(s) mod p, c, recentring, centre)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "slice_mb.h"
#define D 6
#define BND_EPS 1e-9
static double M6[D][D], M5[5][5], Minv5[5][5], mvec[5], shift5[5], dconst;   /* Q_b = b'^T M5 b' + 2 k m^T b' + k^2 M55 */
static double R2;
static long long r5[5], S5[5], L5, K5, KL;   /* box radii, strides, line size, fold constant, number of lines */
static unsigned char *bucket; static size_t NB;
static unsigned char logp = 61;
/* ---------------- generic lattice tools (dimension n <= 6), integer bases, GS in a metric ---------------- */
typedef struct { int n; long long B[D][D]; double MU[D][D], BN[D]; } lat_t;
static double ipn(const lat_t *L, const double *Mx, const long long *x, const long long *y) {
    int n = L->n; double t = 0;
    for (int i = 0; i < n; i++) { double r = 0; for (int j = 0; j < n; j++) r += Mx[i * n + j] * y[j]; t += x[i] * r; }
    return t;
}
static void gso(lat_t *L, const double *Mx) {
    int n = L->n; double Bst[D][D];
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) Bst[i][k] = L->B[i][k];
        for (int j = 0; j < i; j++) {
            double num = 0;
            for (int k = 0; k < n; k++) { double r = 0; for (int l = 0; l < n; l++) r += Mx[k * n + l] * Bst[j][l]; num += L->B[i][k] * r; }
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
/* lattice {x in Z^n : sum x_i s^i = 0 mod p} */
static void hyperplane_lattice(lat_t *L, int n, uint64_t p, uint64_t s) {
    L->n = n; memset(L->B, 0, sizeof L->B); L->B[0][0] = (long long)p; uint64_t si = 1;
    for (int i = 1; i < n; i++) { si = (uint64_t)((__uint128_t)si * s % p); L->B[i][0] = (long long)((p - si) % p); L->B[i][i] = 1; }
}
/* inverse of the basis (rows) as a matrix of doubles: Winv[i][j] such that t_i = sum_j Winv[i][j] v_j are the coordinates of v */
static void basis_inverse(const lat_t *L, double W[D][D]) {
    int n = L->n; double A[D][2 * D];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { A[i][j] = L->B[j][i]; A[i][n + j] = i == j; }
    for (int c = 0; c < n; c++) { int pr = c; for (int r = c + 1; r < n; r++) if (fabs(A[r][c]) > fabs(A[pr][c])) pr = r;
        for (int k = 0; k < 2 * n; k++) { double t = A[c][k]; A[c][k] = A[pr][k]; A[pr][k] = t; }
        double d = A[c][c]; for (int k = 0; k < 2 * n; k++) A[c][k] /= d;
        for (int r = 0; r < n; r++) if (r != c) { double f = A[r][c]; for (int k = 0; k < 2 * n; k++) A[r][k] -= f * A[c][k]; } }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) W[i][j] = A[i][n + j];
}
/* exact recentring: given origin o (integer), return o' = o - sum round(t_i) B_i (integer), small */
static void recentre(const lat_t *L, const double W[D][D], const long long *o, long long *out) {
    int n = L->n; long long q[D];
    for (int i = 0; i < n; i++) { double t = 0; for (int j = 0; j < n; j++) t += W[i][j] * (double)o[j]; q[i] = llround(t); }
    for (int k = 0; k < n; k++) { long long v = o[k]; for (int i = 0; i < n; i++) v -= q[i] * L->B[i][k]; out[k] = v; }
}
/* ---------------- region ---------------- */
static double det_n(const double *A0, int n) {
    double A[D][D]; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = A0[i * n + j];
    double d = 1; for (int c = 0; c < n; c++) { d *= A[c][c]; for (int r = c + 1; r < n; r++) { double f = A[r][c] / A[c][c]; for (int k = c; k < n; k++) A[r][k] -= f * A[c][k]; } } return d;
}
static void inv_n(const double *A0, int n, double *out) {
    double A[D][2 * D]; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { A[i][j] = A0[i * n + j]; A[i][n + j] = i == j; }
    for (int c = 0; c < n; c++) { double d = A[c][c]; for (int k = 0; k < 2 * n; k++) A[c][k] /= d; for (int r = 0; r < n; r++) if (r != c) { double f = A[r][c]; for (int k = 0; k < 2 * n; k++) A[r][k] -= f * A[c][k]; } }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) out[i * n + j] = A[i][n + j];
}
static void setup_region(int log2pts, long long K) {
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) M6[i][j] = MB0[i][j];
    for (int i = 0; i < 5; i++) { for (int j = 0; j < 5; j++) M5[i][j] = M6[i][j]; mvec[i] = M6[i][5]; }
    inv_n((double *)M5, 5, (double *)Minv5);
    for (int i = 0; i < 5; i++) { shift5[i] = 0; for (int j = 0; j < 5; j++) shift5[i] += Minv5[i][j] * mvec[j]; }   /* centre of line k is -k*shift5 */
    dconst = M6[5][5]; for (int i = 0; i < 5; i++) dconst -= mvec[i] * shift5[i];                              /* radius^2 of line k is R2 - k^2 dconst */
    /* choose R2 so that the full 6-dim ellipsoid holds ~2^(log2pts+1) points: vol = vball6 R^6 / sqrt det */
    double vol = ldexp(1.0, log2pts + 1); double det6 = det_n((double *)M6, 6);
    R2 = pow(vol * sqrt(det6) * 6 / (M_PI * M_PI * M_PI), 1.0 / 3);
    double Minv6[D][D]; inv_n((double *)M6, 6, (double *)Minv6);
    for (int i = 0; i < 5; i++) r5[i] = (long long)floor(sqrt(R2 * Minv6[i][i]) + 1e-9);
    long long kmax = (long long)floor(sqrt(R2 * Minv6[5][5]) + 1e-9);
    KL = K > 0 ? K : kmax; if (KL > kmax) KL = kmax;
    L5 = 1; for (int i = 0; i < 5; i++) { S5[i] = L5; L5 *= 2 * r5[i] + 1; }
    K5 = 0; for (int i = 0; i < 5; i++) K5 += 2 * r5[i] * S5[i];
    NB = (size_t)(KL * L5); bucket = calloc(NB, 1);
    fprintf(stderr, "region: R2=%.3f box5=", R2); for (int i = 0; i < 5; i++) fprintf(stderr, "[%lld,%lld]", -r5[i], r5[i]);
    fprintf(stderr, " lines=%lld (of %lld) L5=%lld bytes=%zu\n", KL, kmax, L5, NB);
}
static inline double q6(const long long *b) { double t = 0; for (int i = 0; i < D; i++) { double r = 0; for (int j = 0; j < D; j++) r += M6[i][j] * b[j]; t += b[i] * r; } return t; }
static inline long long local5(const long long *b) { long long a = 0; for (int i = 0; i < 5; i++) a += (b[i] + r5[i]) * S5[i]; return a; }
static uint64_t digest(void) { uint64_t h = 1469598103934665603ULL; for (size_t i = 0; i < NB; i++) { h ^= bucket[i]; h *= 1099511628211ULL; } return h; }
/* ---------------- brute-force reference: ideal (c) and its conjugate (-c) over the half-slice ---------------- */
static long long ref_enum(uint64_t p, uint64_t s, uint64_t c) {
    long long b[D], hits = 0; uint64_t spow[D]; spow[0] = 1; for (int i = 1; i < D; i++) spow[i] = (uint64_t)((__uint128_t)spow[i - 1] * s % p);
    for (b[5] = 1; b[5] <= KL; b[5]++) for (b[4] = -r5[4]; b[4] <= r5[4]; b[4]++) for (b[3] = -r5[3]; b[3] <= r5[3]; b[3]++)
    for (b[2] = -r5[2]; b[2] <= r5[2]; b[2]++) for (b[1] = -r5[1]; b[1] <= r5[1]; b[1]++) {
        uint64_t rest = 0; for (int i = 1; i < D; i++) { uint64_t bi = (uint64_t)(((__int128)b[i] % (__int128)p + p) % p); rest = (uint64_t)(((__uint128_t)rest + (__uint128_t)bi * spow[i]) % p); }
        for (int sign = 0; sign < 2; sign++) {
            uint64_t target = sign ? (p - c) % p : c;            /* conjugate ideal: b(s) = -c */
            uint64_t r0 = (target + p - rest) % p;
            long long b0 = -r5[0] + (long long)((r0 + (uint64_t)(((__int128)r5[0] % (__int128)p + p) % p)) % p);
            for (b[0] = b0; b[0] <= r5[0]; b[0] += (long long)p) {
                if (q6(b) > R2 + 1e-9) continue;
                hits++; bucket[(b[5] - 1) * L5 + local5(b)] += logp;
            }
        }
    }
    return hits;
}
/* ---------------- Schnorr-Euchner walk (generic n) with affine 64-bit state ---------------- */
typedef struct { long long base; long long dl[D]; } aff_t;     /* value = base + sum z_i dl[i] */
/* enumerate z with ||sum z_i B_i - ctr||_M^2 <= rad2 ; ctr given in GS coordinates (cst); on each point call HIT with the
   affine states (addr) and (line) ; implemented as a macro over a body to keep the hot loop simple */
#define WALK(L, cst, rad2, A1, A2, HITSTMT)                                                                       \
    do {                                                                                                          \
        int n = (L)->n; long long z[D], zlo[D], zhi[D]; double C[D][D], rem[D]; long long SA[D], SB[D];           \
        int k = n - 1; for (int i = 0; i < n; i++) C[k][i] = (cst)[i];                                            \
        rem[k] = (rad2); double rr = sqrt(rem[k] / (L)->BN[k]);                                                   \
        zlo[k] = (long long)ceil(C[k][k] - rr - BND_EPS); zhi[k] = (long long)floor(C[k][k] + rr + BND_EPS); z[k] = zlo[k] - 1; \
        SA[k] = (A1).base + (zlo[k] - 1) * (A1).dl[k]; SB[k] = (A2).base + (zlo[k] - 1) * (A2).dl[k];            \
        const long long da0 = (A1).dl[0], db0 = (A2).dl[0];                                                       \
        for (;;) {                                                                                                \
            if (z[k] < zhi[k]) {                                                                                  \
                z[k]++; SA[k] += (A1).dl[k]; SB[k] += (A2).dl[k];                                                 \
                double tt = C[k][k] - (double)z[k]; double r_ = rem[k] - tt * tt * (L)->BN[k]; if (r_ < 0) r_ = 0;   \
                double zz = (double)z[k];                                                                         \
                if (k == 1) {                                                                                     \
                    double c0 = C[1][0] - zz * (L)->MU[1][0]; double r0 = sqrt(r_ / (L)->BN[0]);                  \
                    long long l0 = (long long)ceil(c0 - r0 - BND_EPS), h0 = (long long)floor(c0 + r0 + BND_EPS); \
                    long long a_ = SA[1] + l0 * da0, b_ = SB[1] + l0 * db0;                                       \
                    for (long long j = l0; j <= h0; j++) { HITSTMT; a_ += da0; b_ += db0; }                       \
                    if (h0 >= l0) hits += h0 - l0 + 1;                                                            \
                } else {                                                                                          \
                    rem[k - 1] = r_; for (int i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * (L)->MU[k][i];      \
                    k--; rr = sqrt(rem[k] / (L)->BN[k]);                                                          \
                    zlo[k] = (long long)ceil(C[k][k] - rr - BND_EPS); zhi[k] = (long long)floor(C[k][k] + rr + BND_EPS); z[k] = zlo[k] - 1; \
                    SA[k] = SA[k + 1] + (zlo[k] - 1) * (A1).dl[k]; SB[k] = SB[k + 1] + (zlo[k] - 1) * (A2).dl[k]; \
                }                                                                                                 \
            } else { if (k == n - 1) break; k++; }                                                                \
        }                                                                                                         \
    } while (0)
/* GS centre coordinates of a real vector ctr: solve B^T t = ctr then c*_k = t_k + sum_{i>k} t_i mu_ik */
static void gs_centre(const lat_t *L, const double W[D][D], const double *ctr, double *cst) {
    int n = L->n; double t[D];
    for (int i = 0; i < n; i++) { t[i] = 0; for (int j = 0; j < n; j++) t[i] += W[i][j] * ctr[j]; }
    for (int k = 0; k < n; k++) { double v = t[k]; for (int i = k + 1; i < n; i++) v += t[i] * L->MU[i][k]; cst[k] = v; }
}
/* ---------------- LINE path ---------------- */
static lat_t LAT5; static double W5[D][D];
static void prepare_line_lattice(uint64_t p, uint64_t s) {   /* once per (p, s): factor-base construction in production */
    hyperplane_lattice(&LAT5, 5, p, s); lll(&LAT5, (double *)M5); basis_inverse(&LAT5, W5);
}
static long long line_path(uint64_t p, uint64_t s, uint64_t c, int apply_bytes) {
    long long hits = 0;
    uint64_t s5 = 1; for (int i = 0; i < 5; i++) s5 = (uint64_t)((__uint128_t)s5 * s % p);
    aff_t A1, A2; for (int i = 0; i < 5; i++) { A1.dl[i] = 0; for (int k = 0; k < 5; k++) A1.dl[i] += LAT5.B[i][k] * S5[k]; A2.dl[i] = 0; }
    unsigned char *bk = bucket; const unsigned char lp = logp;
    for (long long k = 1; k <= KL; k++) {
        double rad2 = R2 - (double)k * k * dconst; if (rad2 <= 0) continue;
        for (int sign = 0; sign < 2; sign++) {
            /* coset of the ideal at line +k (sign 0) or of the conjugate at +k (sign 1: b'(s) = -c - k s^5) */
            uint64_t ks5 = (uint64_t)((__uint128_t)(uint64_t)k * s5 % p);
            uint64_t ck = sign ? (p - (c + ks5) % p) % p : (c + p - ks5) % p;
            long long o[D] = {(long long)ck, 0, 0, 0, 0, 0}, oo[D];
            recentre(&LAT5, W5, o, oo);                               /* exact integer recentring */
            double ctr[5], cst[D]; for (int i = 0; i < 5; i++) ctr[i] = -(double)k * shift5[i] - (double)oo[i];
            gs_centre(&LAT5, W5, ctr, cst);
            A1.base = (k - 1) * L5; for (int i = 0; i < 5; i++) A1.base += (oo[i] + r5[i]) * S5[i];
            A2.base = 0;
            if (apply_bytes) { WALK(&LAT5, cst, rad2, A1, A2, { bk[a_] += lp; (void)b_; }); }
            else { WALK(&LAT5, cst, rad2, A1, A2, { __asm__ __volatile__("" :: "r"(a_), "r"(b_)); }); }
        }
    }
    return hits;
}
/* ---------------- RECORD path ---------------- */
#define NSEG 64
static uint64_t *segpos[NSEG]; static uint32_t *segpidx[NSEG]; static uint32_t segn[NSEG]; static long long seg_lines; static uint32_t segcap;
static unsigned char *segbase[NSEG];   /* each segment's own bucket array (here: slices of one allocation) */
static lat_t LAT6; static double W6[D][D];
static void prepare_slice_lattice(uint64_t p, uint64_t s) {  /* once per (p, s): factor-base construction in production */
    hyperplane_lattice(&LAT6, 6, p, s); lll(&LAT6, (double *)M6); basis_inverse(&LAT6, W6);
}
static long long record_path(uint64_t p, uint64_t s, uint64_t c, uint32_t pidx, int emit) {
    long long hits = 0;
    long long o[D] = {(long long)c, 0, 0, 0, 0, 0}, oo[D]; recentre(&LAT6, W6, o, oo);
    double ctr[D], cst[D]; for (int i = 0; i < D; i++) ctr[i] = -(double)oo[i];
    gs_centre(&LAT6, W6, ctr, cst);
    aff_t A1, A2;                                   /* A1 = local 5-dim address, A2 = line index k = b5 */
    for (int i = 0; i < D; i++) { A1.dl[i] = 0; for (int kk = 0; kk < 5; kk++) A1.dl[i] += LAT6.B[i][kk] * S5[kk]; A2.dl[i] = LAT6.B[i][5]; }
    A1.base = 0; for (int kk = 0; kk < 5; kk++) A1.base += (oo[kk] + r5[kk]) * S5[kk]; A2.base = oo[5];
    const long long KK5 = K5, LL5 = L5, SL = seg_lines; const uint32_t cap = segcap;
    if (emit) {
        WALK(&LAT6, cst, R2, A1, A2, {
            long long m_ = b_ >> 63;                              /* b5 < 0: fold */
            long long kf = (b_ ^ m_) - m_;                        /* |b5| */
            long long lf = ((KK5 - a_) & m_) | (a_ & ~m_);        /* local address, folded */
            long long sg = (kf - 1) / SL;                         /* segment index (production: a shift, G a power of two) */
            if (kf >= 1 && sg < NSEG) {                           /* b5 = 0 points are outside the half-slice */
                uint32_t n_ = segn[sg];
                segpos[sg][n_] = (uint64_t)(((kf - 1) - sg * SL) * LL5 + lf);   /* SEGMENT-LOCAL position < G L5 */
                segpidx[sg][n_] = pidx;
                segn[sg] = (n_ + 1) % cap;                        /* benchmark ring; production flushes a full buffer */
            }
        });
    } else {
        WALK(&LAT6, cst, R2, A1, A2, { __asm__ __volatile__("" :: "r"(a_), "r"(b_)); });
    }
    return hits;
}
static void apply_records(void) {
    const unsigned char lp = logp;
    for (int sg = 0; sg < NSEG; sg++) {
        unsigned char *bk = segbase[sg]; const uint64_t *ps = segpos[sg]; uint32_t n = segn[sg];
        for (uint32_t i = 0; i < n; i++) bk[ps[i]] += lp;          /* the prime index is read by survivors only */
        segn[sg] = 0;
    }
}
static void alloc_segments(void) {
    seg_lines = (KL + NSEG - 1) / NSEG; if (seg_lines < 1) seg_lines = 1;
    segcap = 1u << 22;
    for (int i = 0; i < NSEG; i++) { segpos[i] = malloc(sizeof(uint64_t) * segcap); segpidx[i] = malloc(sizeof(uint32_t) * segcap); segn[i] = 0;
        segbase[i] = bucket + (size_t)i * (size_t)seg_lines * (size_t)L5; }
}
int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 2; }
    const char *mode = argv[1];
    if (!strcmp(mode, "check")) {
        int log2pts = atoi(argv[2]); long long K = atoll(argv[3]); uint64_t p = strtoull(argv[4], 0, 10), s = strtoull(argv[5], 0, 10), c = strtoull(argv[6], 0, 10);
        setup_region(log2pts, K); alloc_segments(); prepare_line_lattice(p, s); prepare_slice_lattice(p, s);
        long long h1 = line_path(p, s, c, 1); uint64_t d1 = digest();
        memset(bucket, 0, NB); long long h2 = record_path(p, s, c, 7, 1); apply_records(); uint64_t d2 = digest();
        memset(bucket, 0, NB); long long h3 = ref_enum(p, s, c); uint64_t d3 = digest();
        /* the record path walks the full ellipsoid (b5 = 0 points excluded at emission): its hit count includes b5 = 0 */
        printf("line path: hits=%lld digest=%016llx | record path: walked=%lld digest=%016llx | reference: hits=%lld digest=%016llx | %s\n",
               h1, (unsigned long long)d1, h2, (unsigned long long)d2, h3, (unsigned long long)d3, (d1 == d3 && d2 == d3) ? "MATCH" : "MISMATCH");
        return (d1 == d3 && d2 == d3) ? 0 : 1;
    }
    if (!strcmp(mode, "icount5") || !strcmp(mode, "icount6") || !strcmp(mode, "icount6a")) {
        int log2pts = atoi(argv[2]); long long K = atoll(argv[3]); uint64_t p = strtoull(argv[4], 0, 10), s = strtoull(argv[5], 0, 10), c = strtoull(argv[6], 0, 10); int reps = atoi(argv[7]);
        setup_region(log2pts, K); alloc_segments(); long long tot = 0; prepare_line_lattice(p, s); prepare_slice_lattice(p, s);
        for (int r = 0; r < reps; r++) {
            if (mode[6] == '5') tot += line_path(p, s, c + r, 1);
            else { tot += record_path(p, s, c + r, 7, 1); if (mode[7] == 'a') apply_records(); else for (int i = 0; i < NSEG; i++) segn[i] = 0; }
        }
        printf("reps=%d hits_per_rep=%lld digest=%016llx\n", reps, reps ? tot / reps : 0, (unsigned long long)digest());
        return 0;
    }
    if (!strcmp(mode, "setup6")) {  /* per-(slice, pair) setup of the record path: a(s) mod p (Horner), c = a(s) r^-1, recentring, centre */
        uint64_t p = strtoull(argv[2], 0, 10), s = strtoull(argv[3], 0, 10); long long trials = atoll(argv[4]);
        setup_region(20, 0); prepare_slice_lattice(p, s);
        const uint64_t P = p; const uint64_t mu = (uint64_t)((((__uint128_t)1) << 126) / P);
        #define REDC6(x) ({ __uint128_t _x = (x); uint64_t _q = (uint64_t)((((__uint128_t)(uint64_t)(_x >> 62)) * mu) >> 64); \
                            uint64_t _r = (uint64_t)(_x - (__uint128_t)_q * P); _r = _r >= P ? _r - P : _r; _r = _r >= P ? _r - P : _r; _r; })
        uint64_t rinv = 123456789 % P; long long acc = 0; long long a[6]; for (int i = 0; i < 6; i++) a[i] = (long long)(i * 7919 % 2000) - 1000;
        for (long long t = 0; t < trials; t++) {
            a[0] += 1; if (a[0] > 1000) a[0] = -1000;
            uint64_t as = 0; for (int i = 5; i >= 0; i--) as = REDC6((__uint128_t)as * s + (uint64_t)(a[i] + 1000));
            uint64_t c = REDC6((__uint128_t)as * rinv);
            long long o[D] = {(long long)c, 0, 0, 0, 0, 0}, oo[D]; recentre(&LAT6, W6, o, oo);
            double ctr[D], cst[D]; for (int i = 0; i < D; i++) ctr[i] = -(double)oo[i];
            gs_centre(&LAT6, W6, ctr, cst);
            acc += oo[0] + (long long)cst[0];
        }
        printf("trials=%lld acc=%lld\n", trials, acc);
        return 0;
    }
    if (!strcmp(mode, "setup")) {   /* per-(pair, line, coset) setup of the line path: c_k, recentring, centre */
        uint64_t p = strtoull(argv[2], 0, 10), s = strtoull(argv[3], 0, 10); long long trials = atoll(argv[4]);
        setup_region(20, 0); prepare_line_lattice(p, s);
        uint64_t s5 = 1; for (int i = 0; i < 5; i++) s5 = (uint64_t)((__uint128_t)s5 * s % p);
        const uint64_t P = p; const uint64_t mu = (uint64_t)((((__uint128_t)1) << 126) / P);
        #define REDC(x) ({ __uint128_t _x = (x); uint64_t _q = (uint64_t)((((__uint128_t)(uint64_t)(_x >> 62)) * mu) >> 64); \
                           uint64_t _r = (uint64_t)(_x - (__uint128_t)_q * P); _r = _r >= P ? _r - P : _r; _r = _r >= P ? _r - P : _r; _r; })
        long long acc = 0; uint64_t c = 123456789 % P;
        for (long long t = 0; t < trials; t++) {
            long long k = 1 + (t % 1000);
            uint64_t ck = REDC((__uint128_t)(c + P - REDC((__uint128_t)(uint64_t)k * s5)));
            long long o[D] = {(long long)ck, 0, 0, 0, 0, 0}, oo[D]; recentre(&LAT5, W5, o, oo);
            double ctr[5], cst[D]; for (int i = 0; i < 5; i++) ctr[i] = -(double)k * shift5[i] - (double)oo[i];
            gs_centre(&LAT5, W5, ctr, cst);
            acc += oo[0] + (long long)cst[0];
        }
        printf("trials=%lld acc=%lld\n", trials, acc);
        return 0;
    }
    fprintf(stderr, "unknown mode\n"); return 2;
}
