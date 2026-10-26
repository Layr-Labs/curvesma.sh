/* slice6: ARCH-E a-slice sieve kernel (6-dim b-space), retired-instruction benchmark.
 *
 * Region: half-ellipsoid {b in Z^6 : b^T M b <= R2, b5 >= 0} with M = MB0 (fixed shape) scaled; bucket array
 * over the bounding box of the half-region with fixed-cube addressing addr = sum_i (b_i - lo_i) * S_i (injective).
 * Ideal coset: L(p,s) = {b : sum b_i s^i = 0 mod p}, coset c*e0 + L (all degree-1 ideals above (p, y-s)
 * share L; c = a(s)/r mod p is the only ideal-and-slice-dependent datum).
 * The sieve walks the coset over the FULL ellipsoid (both signs of b5) and FOLDS b -> -b when b5 < 0
 * (addr(-b) = K - addr(b)), which supplies the conjugate ideal's updates (x -> -x automorphism).
 *
 * Modes:
 *   slice6 check  log2pts p s c          fast enumeration vs brute-force reference: bit-identical bucket digest
 *   slice6 icount log2pts p s c reps     bucket updates, `reps` repetitions (count with /usr/bin/time -l)
 *   slice6 babai  log2pts p s trials     large-p regime: Babai nearest-coset-point test vs exact enumeration
 *   slice6 setup  p s trials             per-(slice,pair) setup: c = a(s) r^-1 mod p and Babai fractional coords
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include "slice_mb.h"
#define D 6
static double M[D][D];           /* metric */
static double R2;
static long long B[D][D];        /* lattice basis rows (integer, LLL-reduced in metric M) */
static double MU[D][D], BN[D];   /* Gram-Schmidt in metric M: mu[i][j], ||b*_i||^2 */
static long long lo[D], hi[D], S[D]; static long long K; static size_t NB;
static unsigned char *bucket; static unsigned char logp = 61;
static long long p_, s_;

static double ip(const long long *x, const long long *y) {   /* <x,y>_M */
    double t = 0; for (int i = 0; i < D; i++) { double r = 0; for (int j = 0; j < D; j++) r += M[i][j] * y[j]; t += x[i] * r; } return t;
}
static void gso(void) {
    double Bst[D][D];
    for (int i = 0; i < D; i++) {
        for (int k = 0; k < D; k++) Bst[i][k] = B[i][k];
        for (int j = 0; j < i; j++) {
            double num = 0, den = BN[j];
            for (int k = 0; k < D; k++) { double r = 0; for (int l = 0; l < D; l++) r += M[k][l] * Bst[j][l]; num += B[i][k] * r; }
            MU[i][j] = num / den;
            for (int k = 0; k < D; k++) Bst[i][k] -= MU[i][j] * Bst[j][k];
        }
        double n = 0; for (int k = 0; k < D; k++) { double r = 0; for (int l = 0; l < D; l++) r += M[k][l] * Bst[i][l]; n += Bst[i][k] * r; }
        BN[i] = n;
    }
}
static void lll(void) {   /* textbook LLL (delta=0.99) in metric M, integer basis */
    gso();
    int k = 1;
    while (k < D) {
        for (int j = k - 1; j >= 0; j--) {
            long long q = llround(MU[k][j]);
            if (q) { for (int i = 0; i < D; i++) B[k][i] -= q * B[j][i]; gso(); }
        }
        if (BN[k] >= (0.99 - MU[k][k - 1] * MU[k][k - 1]) * BN[k - 1]) k++;
        else { long long t[D]; memcpy(t, B[k], sizeof t); memcpy(B[k], B[k - 1], sizeof t); memcpy(B[k - 1], t, sizeof t); gso(); k = k > 1 ? k - 1 : 1; }
    }
}
static void setup_lattice(long long p, long long s) {
    /* basis of {b : sum b_i s^i = 0 mod p}: (p,0,..), (-s^i mod p, e_i) */
    memset(B, 0, sizeof B);
    B[0][0] = p;
    long long si = 1;
    for (int i = 1; i < D; i++) { si = (__int128)si * s % p; B[i][0] = (p - si) % p; B[i][i] = 1; }
    lll();
}
static void setup_region(int log2pts) {
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) M[i][j] = MB0[i][j];
    /* scale R2 so the half-ellipsoid holds ~2^log2pts points: vol = pi^3/6 * R^6 / sqrt(det M) / 2 */
    double det = 1; { double A[D][D]; memcpy(A, M, sizeof A);
        for (int c = 0; c < D; c++) { det *= A[c][c]; for (int r = c + 1; r < D; r++) { double f = A[r][c] / A[c][c]; for (int k = c; k < D; k++) A[r][k] -= f * A[c][k]; } } }
    double vol = ldexp(1.0, log2pts) * 2;
    double R6 = vol * sqrt(det) * 6 / (M_PI * M_PI * M_PI);
    R2 = pow(R6, 1.0 / 3);
    /* bounding box: |b_i| <= sqrt(R2 * (M^-1)_ii) */
    double Minv[D][D]; { double A[D][2 * D]; for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) { A[i][j] = M[i][j]; A[i][D + j] = i == j; }
        for (int c = 0; c < D; c++) { double d = A[c][c]; for (int k = 0; k < 2 * D; k++) A[c][k] /= d; for (int r = 0; r < D; r++) if (r != c) { double f = A[r][c]; for (int k = 0; k < 2 * D; k++) A[r][k] -= f * A[c][k]; } }
        for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) Minv[i][j] = A[i][D + j]; }
    for (int i = 0; i < D; i++) { long long e = (long long)floor(sqrt(R2 * Minv[i][i]) + 1e-9); lo[i] = -e; hi[i] = e; }
    lo[D - 1] = 0;
    NB = 1; for (int i = 0; i < D; i++) { S[i] = (long long)NB; NB *= (size_t)(hi[i] - lo[i] + 1); }
    K = 0; for (int i = 0; i < D; i++) K += (-lo[i] - lo[i]) * S[i];   /* addr(-b) = K - addr(b): sum_i (-b_i - lo_i) S_i = sum_i(-2 lo_i) S_i - addr */
    bucket = calloc(NB, 1);
    fprintf(stderr, "region: R2=%.3f bbox=", R2); for (int i = 0; i < D; i++) fprintf(stderr, "[%lld,%lld]", lo[i], hi[i]); fprintf(stderr, " bytes=%zu K=%lld\n", NB, K);
}
static inline long long addr_of(const long long *b) { long long a = 0; for (int i = 0; i < D; i++) a += (b[i] - lo[i]) * S[i]; return a; }
static inline double quad(const long long *b) { return ip(b, b); }

/* ---- reference: brute force over the bbox (both signs of b5) ---- */
static long long ref_enum(long long p, long long s, long long c) {
    long long b[D], hits = 0; long long spow[D]; spow[0] = 1; for (int i = 1; i < D; i++) spow[i] = (__int128)spow[i - 1] * s % p;
    long long full_lo5 = -hi[D - 1];
    for (b[5] = full_lo5; b[5] <= hi[5]; b[5]++) for (b[4] = lo[4]; b[4] <= hi[4]; b[4]++) for (b[3] = lo[3]; b[3] <= hi[3]; b[3]++)
    for (b[2] = lo[2]; b[2] <= hi[2]; b[2]++) for (b[1] = lo[1]; b[1] <= hi[1]; b[1]++) {
        /* solve for b0: b0 = c - sum_{i>=1} b_i s^i mod p; enumerate b0 in [lo0,hi0] in that class */
        long long r = c; for (int i = 1; i < D; i++) r = ((r - (__int128)b[i] * spow[i]) % p + p) % p;
        long long b0 = lo[0] + (((r - lo[0]) % p) + p) % p;
        for (b[0] = b0; b[0] <= hi[0]; b[0] += p) {
            if (quad(b) > R2 + 1e-9) continue;
            hits++;
            long long a = addr_of(b);
            if (b[5] < 0) a = K - a;
            bucket[a] += logp;
        }
    }
    return hits;
}

/* ---- fast: Schnorr-Euchner over the coset (center shift), carrying addr and b5 ---- */
static double CST[D];            /* coset offset coordinates in the basis (negated: center of the shifted ellipsoid) */
static long long DA[D], D5[D];   /* per-level deltas of addr and b5 */
static long long ADDR0, B50;
static double CST[D]; static long long DA[D], D5[D]; static long long ADDR0, B50;
static double WCOL[D];           /* column 0 of the basis inverse: t = WCOL * c are the basis coordinates of (c,0,..,0) */
static void setup_coset_tables(void) {
    double A[D][2 * D]; for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) { A[i][j] = B[j][i]; A[i][D + j] = i == j; }
    for (int c2 = 0; c2 < D; c2++) { int pr = c2; for (int r = c2 + 1; r < D; r++) if (fabs(A[r][c2]) > fabs(A[pr][c2])) pr = r;
        for (int k = 0; k < 2 * D; k++) { double t = A[c2][k]; A[c2][k] = A[pr][k]; A[pr][k] = t; }
        double d = A[c2][c2]; for (int k = 0; k < 2 * D; k++) A[c2][k] /= d;
        for (int r = 0; r < D; r++) if (r != c2) { double f = A[r][c2]; for (int k = 0; k < 2 * D; k++) A[r][k] -= f * A[c2][k]; } }
    for (int i = 0; i < D; i++) WCOL[i] = A[i][D + 0];
    for (int i = 0; i < D; i++) { DA[i] = 0; for (int k = 0; k < D; k++) DA[i] += B[i][k] * S[k]; D5[i] = B[i][5]; }
}
/* per-coset work (production hot path): t = -c*WCOL (6 mul), GS centre (15 fma), addr/b5 origin */
static inline void setup_coset(long long c) {
    double t[D]; for (int i = 0; i < D; i++) t[i] = -WCOL[i] * (double)c;
    for (int k = 0; k < D; k++) { double v = t[k]; for (int i = k + 1; i < D; i++) v += t[i] * MU[i][k]; CST[k] = v; }
    ADDR0 = (c - lo[0]) * S[0]; for (int k = 1; k < D; k++) ADDR0 += (0 - lo[k]) * S[k]; B50 = 0;
}
#define BND_EPS 1e-9
static long long fast_enum(void) {
    int z[D], zlo[D], zhi[D]; double C[D][D], rem[D]; long long SA[D], S5[D];
    int k = D - 1; long long hits = 0;
    for (int i = 0; i < D; i++) C[k][i] = CST[i];
    rem[k] = R2;
    double rr = sqrt(rem[k] / BN[k]);
    zlo[k] = (int)ceil(C[k][k] - rr - BND_EPS); zhi[k] = (int)floor(C[k][k] + rr + BND_EPS); z[k] = zlo[k] - 1;
    SA[k] = ADDR0 + (long long)(zlo[k] - 1) * DA[k]; S5[k] = B50 + (long long)(zlo[k] - 1) * D5[k];
    const long long da0 = DA[0], d50 = D5[0], KK = K; unsigned char *bk = bucket; const unsigned char lp = logp;
    for (;;) {
        if (z[k] < zhi[k]) {
            z[k]++; SA[k] += DA[k]; S5[k] += D5[k];
            double tt = C[k][k] - (double)z[k];
            double r = rem[k] - tt * tt * BN[k]; if (r < 0) r = 0;
            double zz = (double)z[k];
            if (k == 1) {
                double c0 = C[1][0] - zz * MU[1][0];
                double r0 = sqrt(r / BN[0]);
                int l0 = (int)ceil(c0 - r0 - BND_EPS), h0 = (int)floor(c0 + r0 + BND_EPS);
                long long a = SA[1] + (long long)l0 * da0, b5 = S5[1] + (long long)l0 * d50;
                for (int j = l0; j <= h0; j++) {
                    long long m = b5 >> 63;                 /* all-ones if b5 < 0 */
                    long long af = ((KK - a) & m) | (a & ~m);
                    bk[af] += lp;
                    a += da0; b5 += d50;
                }
                if (h0 >= l0) hits += h0 - l0 + 1;
            } else {
                rem[k - 1] = r;
                for (int i = 0; i < k; i++) C[k - 1][i] = C[k][i] - zz * MU[k][i];
                k--;
                rr = sqrt(rem[k] / BN[k]);
                zlo[k] = (int)ceil(C[k][k] - rr - BND_EPS); zhi[k] = (int)floor(C[k][k] + rr + BND_EPS); z[k] = zlo[k] - 1;
                SA[k] = SA[k + 1] + (long long)(zlo[k] - 1) * DA[k]; S5[k] = S5[k + 1] + (long long)(zlo[k] - 1) * D5[k];
            }
        } else { if (k == D - 1) break; k++; }
    }
    return hits;
}
static uint64_t digest(void) { uint64_t h = 1469598103934665603ULL; for (size_t i = 0; i < NB; i++) { h ^= bucket[i]; h *= 1099511628211ULL; } return h; }

/* ---- Babai regime: nearest coset point by rounding in the reduced basis ---- */
static double Winv[D][D];  /* inverse of basis (rows) for fractional coordinates */
static double G[D][D];     /* Gram of basis in metric M */
static void setup_babai(void) {
    double A[D][2 * D]; for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) { A[i][j] = B[j][i]; A[i][D + j] = i == j; }   /* columns = basis vectors */
    for (int c = 0; c < D; c++) { int pr = c; for (int r = c + 1; r < D; r++) if (fabs(A[r][c]) > fabs(A[pr][c])) pr = r;
        for (int k = 0; k < 2 * D; k++) { double t = A[c][k]; A[c][k] = A[pr][k]; A[pr][k] = t; }
        double d = A[c][c]; for (int k = 0; k < 2 * D; k++) A[c][k] /= d;
        for (int r = 0; r < D; r++) if (r != c) { double f = A[r][c]; for (int k = 0; k < 2 * D; k++) A[r][k] -= f * A[c][k]; } }
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) Winv[i][j] = A[i][D + j];   /* t = Winv * off  gives off = sum t_j B_j */
    for (int i = 0; i < D; i++) for (int j = 0; j < D; j++) G[i][j] = ip(B[i], B[j]);
}
/* returns 1 and the point if the rounded coset point is inside; cost-relevant path */
static inline int babai_point(long long c, long long *out, double *q) {
    double t[D]; for (int i = 0; i < D; i++) t[i] = Winv[i][0] * (double)c;   /* off = (c,0..0): only column 0 */
    double f[D]; for (int i = 0; i < D; i++) f[i] = t[i] - nearbyint(t[i]);
    double qq = 0; for (int i = 0; i < D; i++) { double r = 0; for (int j = 0; j < D; j++) r += G[i][j] * f[j]; qq += f[i] * r; }
    *q = qq;
    if (qq > R2) return 0;
    for (int k = 0; k < D; k++) { double v = 0; for (int i = 0; i < D; i++) v += f[i] * B[i][k]; out[k] = llround(v); }
    return 1;
}
static uint64_t rng_s = 0x9E3779B97F4A7C15ULL;
static inline uint64_t rng(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 7; rng_s ^= rng_s << 17; return rng_s; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage\n"); return 2; }
    const char *mode = argv[1];
    if (!strcmp(mode, "check") || !strcmp(mode, "icount")) {
        int log2pts = atoi(argv[2]); long long p = atoll(argv[3]), s = atoll(argv[4]), c = atoll(argv[5]);
        setup_region(log2pts); setup_lattice(p, s); setup_coset_tables(); setup_coset(c);
        if (!strcmp(mode, "check")) {
            long long h1 = fast_enum(); uint64_t d1 = digest();
            memset(bucket, 0, NB); long long h2 = ref_enum(p, s, c); uint64_t d2 = digest();
            printf("fast hits=%lld digest=%016llx | ref hits=%lld digest=%016llx | %s\n", h1, (unsigned long long)d1, h2, (unsigned long long)d2, d1 == d2 && h1 == h2 ? "MATCH" : "MISMATCH");
            return d1 == d2 ? 0 : 1;
        }
        int reps = atoi(argv[6]); long long tot = 0;
        for (int r = 0; r < reps; r++) tot += fast_enum();
        printf("reps=%d hits_per_rep=%lld digest=%016llx\n", reps, reps ? tot / reps : 0, (unsigned long long)digest());
        return 0;
    }
    if (!strcmp(mode, "icountc")) {   /* icountc log2pts p s ncos : enumerate cosets c=1..ncos (per-coset setup + walk) */
        int log2pts = atoi(argv[2]); long long p = atoll(argv[3]), s = atoll(argv[4]); long long ncos = atoll(argv[5]);
        setup_region(log2pts); setup_lattice(p, s); setup_coset_tables();
        long long tot = 0;
        for (long long c = 1; c <= ncos; c++) { setup_coset(c % p); tot += fast_enum(); }
        printf("ncos=%lld hits=%lld hits/coset=%.4f digest=%016llx\n", ncos, tot, (double)tot / ncos, (unsigned long long)digest());
        return 0;
    }
    if (!strcmp(mode, "babai")) {
        int log2pts = atoi(argv[2]); long long p = atoll(argv[3]), s = atoll(argv[4]); int trials = atoi(argv[5]);
        setup_region(log2pts); setup_lattice(p, s); setup_coset_tables(); setup_babai();
        int found = 0, exist = 0, agree = 0, inside_found = 0;
        for (int t = 0; t < trials; t++) {
            long long c = (long long)(rng() % (uint64_t)p);
            long long pt[D]; double q;
            int f = babai_point(c, pt, &q);
            /* exact: enumerate the coset in the full ellipsoid (fast_enum counts hits) */
            setup_coset(c); memset(bucket, 0, NB); long long h = fast_enum();
            if (h > 0) exist++;
            if (f) { found++; if (quad(pt) <= R2 + 1e-9) inside_found++; }
            if ((h > 0) == f) agree++;
        }
        printf("p=%lld trials=%d exact-exists=%d babai-found=%d (inside-verified %d) agree=%d  miss-rate=%.3f false-pos=%.3f\n", p, trials, exist, found, inside_found, agree,
               exist ? 1.0 - (double)found / exist : 0.0, found > exist ? (double)(found - exist) / trials : 0.0);
        return 0;
    }
    if (!strcmp(mode, "setup")) {
        /* Barrett: mu = floor(2^126 / p); for x < p^2 < 2^124: q = ((x >> 62) * mu) >> 64 ; x -= q p ; <= 2 corrections */
        /* per-(slice, pair) setup cost: a(s) mod p by Horner (5 steps), c = a(s) * rinv mod p, Babai fractional coords + Gram test */
        long long p = atoll(argv[2]), s = atoll(argv[3]); long long trials = atoll(argv[4]);
        setup_region(20); setup_lattice(p, s); setup_coset_tables(); setup_babai();
        uint64_t rinv = 123456789 % p; long long acc = 0; double qsum = 0;
        long long a[6]; for (int i = 0; i < 6; i++) a[i] = (long long)(rng() % 2000) - 1000;
        const uint64_t P = (uint64_t)p; const uint64_t mu = (uint64_t)((((__uint128_t)1) << 126) / P);
        #define REDC(x) ({ __uint128_t _x = (x); uint64_t _q = (uint64_t)((((__uint128_t)(uint64_t)(_x >> 62)) * mu) >> 64); \
                           uint64_t _r = (uint64_t)(_x - (__uint128_t)_q * P); _r = _r >= P ? _r - P : _r; _r = _r >= P ? _r - P : _r; _r; })
        for (long long t = 0; t < trials; t++) {
            a[0] += 1; if (a[0] > 1000) a[0] = -1000;          /* sweep a */
            uint64_t as = 0;
            for (int i = 5; i >= 0; i--) as = REDC((__uint128_t)as * (uint64_t)s + (uint64_t)(a[i] + 1000));
            long long c = (long long)REDC((__uint128_t)as * rinv);
            long long pt[D]; double q; int f = babai_point(c, pt, &q);
            acc += f; qsum += q;
        }
        printf("trials=%lld found=%lld qsum=%.3f\n", trials, acc, qsum);
        return 0;
    }
    fprintf(stderr, "unknown mode\n"); return 2;
}
