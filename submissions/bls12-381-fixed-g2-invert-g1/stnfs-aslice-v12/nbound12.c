/* nbound12: certified per-block bounds on 16 log2 |N_f| and 16 log2 |N_g| for the scan thresholds (stnfs-aslice-v12).
 *
 * Operation-for-operation C port of embedding_bounds.NormBounds (block, row, adaptive_row with singleton_policy='reject'):
 * the four places of K_h (two real roots, two complex roots of h with positive imaginary part) are certified intervals
 * (rational Rouche disks, rounded outward); every interval operation is evaluated in binary64 round-to-nearest and
 * widened by one ulp on each side (IEEE nextafter, implemented on the bit pattern), so every enclosure is outward.
 * For a block b0 in [lo, hi] of a row (b1..b5 fixed) and each place v: bv = b(v), w = bv^2, F_v = sum_k c_k(v) w^k (the
 * x-homogenised f evaluated at a(v), v), G_v = a(v)^2 + |u| v w, and 16 log2 |N| is bounded by summing per-place
 * integer bounds of 16 log2 |F_v| (|F_v|^2 for complex places) from a 17-entry table of 2^(j/16) bounds.
 * A block is ACCEPTED when both sides' bounds are finite and at most 32 (two bits) wide; otherwise it is halved
 * (floor midpoint) down to single points, and an unresolved single point is REJECTED (dropped, never a survivor).
 * Each row is cut into aligned 16-point blocks (b0 // 16), so a block costs at most 31 evaluations.
 *
 * The thresholds of an accepted leaf are T = lower_bound - 16 slack (preliminary survivor test, see proof.md);
 * a rejected point gets T = 0xFFFF on the f side, which no sum reaches (sums are <= 16 * 1302 < 65535).
 *
 * Modes:
 *   check <vectors>             replay the Python reference vectors (gen_nbound_vectors.py): every leaf and every bound
 *                               must be identical; prints the weighted evaluation density
 *   explore <vectors> L0 W [reps]   aligned L0-blocks of the full rows around each pool sample, acceptance width W
 *   dump <vectors> L0 W         explore, plus one line per pool row (evaluations, leaves, sample rejected, sample bounds)
 *   bench <vectors> reps        the same blocks, reps times (whole-program callgrind Ir at two values of reps)
 *   benchsetup <vectors> reps   only the per-a and per-row setups (the difference with bench prices the blocks)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#define NOINLINE __attribute__((noinline))

typedef struct { double lo, hi; } iv;
typedef struct { iv re, im; } civ;
static const double TINY = 0x1p-1074;
/* IEEE nextafter(x, -inf) / nextafter(x, +inf) for finite x, on the bit pattern (one ulp; +-0 -> -+2^-1074) */
static inline double dn(double x) {
    uint64_t u; memcpy(&u, &x, 8); uint64_t s = (uint64_t)((int64_t)u >> 63);
    uint64_t r = u - 1 + (s & 2); r = (u << 1) == 0 ? 0x8000000000000001ULL : r; memcpy(&x, &r, 8); return x;
}
static inline double up(double x) {
    uint64_t u; memcpy(&u, &x, 8); uint64_t s = (uint64_t)((int64_t)u >> 63);
    uint64_t r = u + 1 - (s & 2); r = (u << 1) == 0 ? 1ULL : r; memcpy(&x, &r, 8); return x;
}
/* No intermediate can overflow: with |a_i|, |b_i| < 2^12 (guarded in setup_a / row_setup) and |roots| < 1.4, every place
   value is < 2^215 (|F_v|^2 < 2^430), far below 2^1024; log_bounds still rejects any non-finite or empty magnitude. */
static int nonfinite;   /* set when an input is outside the guarded range */
static inline iv mk(double lo, double hi) { iv r = {lo, hi}; return r; }
static inline iv ivrat(double c) { return mk(dn(c), up(c)); }                  /* Interval.rational of an exact double */
static inline iv ivadd(iv a, iv b) { return mk(dn(a.lo + b.lo), up(a.hi + b.hi)); }
static inline iv ivneg(iv a) { iv r = {-a.hi, -a.lo}; return r; }
static inline iv ivsub(iv a, iv b) { return ivadd(a, ivneg(b)); }
static inline double mn2(double a, double b) { return b < a ? b : a; }
static inline double mx2(double a, double b) { return b > a ? b : a; }
static inline iv ivmul(iv a, iv b) {
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    return mk(dn(mn2(mn2(p1, p2), mn2(p3, p4))), up(mx2(mx2(p1, p2), mx2(p3, p4))));
}
static inline iv ivabs(iv a) {
    double l = fabs(a.lo), h = fabs(a.hi); iv r;
    r.lo = (a.lo <= 0 && 0 <= a.hi) ? 0.0 : mn2(l, h); r.hi = mx2(l, h); return r;
}
static inline iv ivsq(iv a) { iv b = ivabs(a); return mk(mx2(0.0, dn(b.lo * b.lo)), up(b.hi * b.hi)); }
static inline civ cmk(iv re, iv im) { civ r = {re, im}; return r; }
static inline civ cconst(double c) { iv z = {0, 0}; return cmk(ivrat(c), z); }
static inline civ cadd(civ a, civ b) { return cmk(ivadd(a.re, b.re), ivadd(a.im, b.im)); }
static inline civ cmul(civ a, civ b) {
    iv A = ivmul(a.re, b.re), B = ivmul(a.im, b.im), re = ivsub(A, B);
    iv C = ivmul(a.re, b.im), E = ivmul(a.im, b.re), im = ivadd(C, E);
    return cmk(re, im);
}
static inline iv cabs2(civ a) { iv v = ivadd(ivsq(a.re), ivsq(a.im)); v.lo = mx2(0.0, v.lo); return v; }

/* 16 log2 bounds: lower = 16 (e - 1) + #{j : LU[j] <= 2m} - 1, upper = 16 (e - 1) + #{j : LL[j] < 2m}, x = m 2^e */
static double LL[17], LU[17];
static inline int count_le(const double *T, double x) { int n = 0; for (int k = 16; k; k >>= 1) if (n + k <= 17 && T[n + k - 1] <= x) n += k; return n; }
static inline int count_lt(const double *T, double x) { int n = 0; for (int k = 16; k; k >>= 1) if (n + k <= 17 && T[n + k - 1] < x) n += k; return n; }
static inline void frexp2(double x, double *m2, int *e) {          /* 2m in [1,2) and frexp exponent e */
    uint64_t u; memcpy(&u, &x, 8); int be = (int)((u >> 52) & 0x7FF);
    if (be == 0) { int ee; double m = frexp(x, &ee); *m2 = 2 * m; *e = ee; return; }   /* subnormal: library path */
    u = (u & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL; memcpy(m2, &u, 8); *e = be - 1022;
}
static inline int log_bounds(iv m, int *lo, int *hi) {
    if (!(m.lo > 0) || !isfinite(m.lo) || !isfinite(m.hi) || m.hi < m.lo) return 0;
    double x; int e;
    frexp2(m.lo, &x, &e); *lo = 16 * (e - 1) + count_le(LU, x) - 1;
    frexp2(m.hi, &x, &e); *hi = 16 * (e - 1) + count_lt(LL, x);
    return 1;
}

/* ---------------- per-a context (once per slice) ---------------- */
#define NPL 4
typedef struct { int real; iv rr; civ cr; iv one_r; civ one_c; iv a2r, uyr, cr_[7]; civ a2c, uyc, cc[7]; } place_t;
static place_t PL[NPL]; static int npl;
static iv ROOT_RE[NPL], ROOT_IM[NPL]; static int ROOT_REAL[NPL];
static const double UABS = 15132376222941642752.0;   /* |u| rounds exactly to this double; ivrat widens it */
static void setup_a(const long long *a) {
    npl = NPL; for (int i = 0; i < 6; i++) if (a[i] <= -4096 || a[i] >= 4096) nonfinite = 1;
    for (int v = 0; v < NPL; v++) {
        place_t *P = &PL[v]; P->real = ROOT_REAL[v];
        if (P->real) {
            iv root = ROOT_RE[v], one = {1, 1}, av = ivmul(one, ivrat(0.0)); P->rr = root; P->one_r = one;
            for (int i = 5; i >= 0; i--) av = ivadd(ivmul(av, root), ivrat((double)a[i]));
            iv ap[13], yp[7]; ap[0] = one; yp[0] = one;
            for (int k = 1; k <= 12; k++) ap[k] = ivmul(ap[k - 1], av);
            for (int k = 1; k <= 6; k++) yp[k] = ivmul(yp[k - 1], root);
            P->cr_[0] = ap[12]; P->cr_[1] = ivmul(ivmul(root, ivrat(-2.0)), ap[10]); P->cr_[2] = ivmul(one, ivrat(0.0));
            P->cr_[3] = ivmul(ivmul(yp[3], ivrat(2.0)), ap[6]); P->cr_[4] = ivmul(one, ivrat(0.0));
            P->cr_[5] = ivmul(yp[5], ap[2]); P->cr_[6] = yp[6];
            P->a2r = ap[2]; P->uyr = ivmul(root, ivrat(UABS));
        } else {
            civ root = cmk(ROOT_RE[v], ROOT_IM[v]), one = cconst(1.0), av = cmul(one, cconst(0.0)); P->cr = root; P->one_c = one;
            for (int i = 5; i >= 0; i--) av = cadd(cmul(av, root), cconst((double)a[i]));
            civ ap[13], yp[7]; ap[0] = one; yp[0] = one;
            for (int k = 1; k <= 12; k++) ap[k] = cmul(ap[k - 1], av);
            for (int k = 1; k <= 6; k++) yp[k] = cmul(yp[k - 1], root);
            P->cc[0] = ap[12]; P->cc[1] = cmul(cmul(root, cconst(-2.0)), ap[10]); P->cc[2] = cmul(one, cconst(0.0));
            P->cc[3] = cmul(cmul(yp[3], cconst(2.0)), ap[6]); P->cc[4] = cmul(one, cconst(0.0));
            P->cc[5] = cmul(yp[5], ap[2]); P->cc[6] = yp[6];
            P->a2c = ap[2]; P->uyc = cmul(root, cconst(UABS));
        }
    }
}
/* ---------------- per-row context (once per row of ~1300-1850 cells) ---------------- */
typedef struct { iv br[NPL]; civ bc[NPL]; } rowctx_t;
static NOINLINE void row_setup(const long long *suffix, rowctx_t *R) {
    for (int i = 0; i < 5; i++) if (suffix[i] <= -4096 || suffix[i] >= 4096) nonfinite = 1;
    for (int v = 0; v < NPL; v++) {
        place_t *P = &PL[v];
        if (P->real) { iv bv = ivmul(P->one_r, ivrat(0.0)); for (int i = 4; i >= 0; i--) bv = ivadd(ivmul(bv, P->rr), ivrat((double)suffix[i])); R->br[v] = ivmul(bv, P->rr); }
        else { civ bv = cmul(P->one_c, cconst(0.0)); for (int i = 4; i >= 0; i--) bv = cadd(cmul(bv, P->cr), cconst((double)suffix[i])); R->bc[v] = cmul(bv, P->cr); }
    }
}
/* ---------------- one block evaluation ---------------- */
typedef struct { int ok[2], lo[2], hi[2]; } logs_t;
static long long EVALS;
static NOINLINE void block_eval(const rowctx_t *R, long long lo, long long hi, logs_t *L) {
    EVALS++;
    L->ok[0] = L->ok[1] = 1; L->lo[0] = L->lo[1] = L->hi[0] = L->hi[1] = 0;
    iv b0 = {(double)lo, (double)hi}, z = {0, 0};
    for (int v = 0; v < NPL; v++) {
        const place_t *P = &PL[v]; iv mag[2];
        if (P->real) {
            iv bv = ivadd(R->br[v], b0), w = ivmul(bv, bv), fv = P->cr_[6];
            for (int k = 5; k >= 0; k--) { fv = ivmul(fv, w); if (k != 2 && k != 4) fv = ivadd(fv, P->cr_[k]); }
            iv gv = ivadd(P->a2r, ivmul(P->uyr, w));
            mag[0] = ivabs(fv); mag[1] = ivabs(gv);
        } else {
            civ bv = cadd(R->bc[v], cmk(b0, z)), w = cmul(bv, bv), fv = P->cc[6];
            for (int k = 5; k >= 0; k--) { fv = cmul(fv, w); if (k != 2 && k != 4) fv = cadd(fv, P->cc[k]); }
            civ gv = cadd(P->a2c, cmul(P->uyc, w));
            mag[0] = cabs2(fv); mag[1] = cabs2(gv);
        }
        for (int s = 0; s < 2; s++) {
            int l, h;
            if (!log_bounds(mag[s], &l, &h)) L->ok[s] = 0;
            else if (L->ok[s]) { L->lo[s] += l; L->hi[s] += h; }
        }
    }
}
/* ---------------- adaptive aligned block (<= 16 points), reject policy; emits leaves ---------------- */
typedef struct { long long lo, hi; int rejected; logs_t L; } leaf_t;
static int WIDTH = 32;   /* acceptance width in 1/16 bit (the v11 reference uses 32) */
static NOINLINE int adaptive_block(const rowctx_t *R, long long lo, long long hi, leaf_t *out) {
    long long st[64][2]; int sp = 0, n = 0; st[sp][0] = lo; st[sp][1] = hi; sp++;
    const int W = WIDTH;
    while (sp) {
        sp--; long long l = st[sp][0], h = st[sp][1]; logs_t L; block_eval(R, l, h, &L);
        if (L.ok[0] && L.ok[1] && L.hi[0] - L.lo[0] <= W && L.hi[1] - L.lo[1] <= W) { out[n].lo = l; out[n].hi = h; out[n].rejected = 0; out[n].L = L; n++; }
        else if (l == h) { out[n].lo = l; out[n].hi = h; out[n].rejected = 1; n++; }
        else { long long m = (l + h) >> 1; st[sp][0] = m + 1; st[sp][1] = h; sp++; st[sp][0] = l; st[sp][1] = m; sp++; }
    }
    return n;
}

/* ---------------- vector replay ---------------- */
typedef struct { long long a[6], suf[5], lo, hi, rlo, rhi, b0; double wt; int evals, nleaf; leaf_t *leaf; } case_t;
static case_t *CS; static int ncase; static double ref_density;
static void load(const char *path) {
    FILE *f = fopen(path, "r"); if (!f) { perror(path); exit(2); }
    char tok[64]; int cap = 4096; CS = malloc(sizeof(case_t) * cap); int pv = 0;
    while (fscanf(f, "%63s", tok) == 1) {
        if (!strcmp(tok, "places")) { int n; if (fscanf(f, "%d", &n) != 1 || n != NPL) exit(3); }
        else if (!strcmp(tok, "place")) {
            char s[4][64]; int real; if (fscanf(f, "%d %63s %63s %63s %63s", &real, s[0], s[1], s[2], s[3]) != 5) exit(3);
            ROOT_REAL[pv] = real; ROOT_RE[pv].lo = strtod(s[0], 0); ROOT_RE[pv].hi = strtod(s[1], 0);
            ROOT_IM[pv].lo = strtod(s[2], 0); ROOT_IM[pv].hi = strtod(s[3], 0); pv++;
        }
        else if (!strcmp(tok, "logtab")) { char s[64]; for (int j = 0; j < 17; j++) { if (fscanf(f, "%63s", s) != 1) exit(3); LL[j] = strtod(s, 0); } for (int j = 0; j < 17; j++) { if (fscanf(f, "%63s", s) != 1) exit(3); LU[j] = strtod(s, 0); } }
        else if (!strcmp(tok, "row")) { case_t *c = &CS[ncase]; int sd, ix;
            if (fscanf(f, "%lld %lld %lld %lf %d %d", &c->rlo, &c->rhi, &c->b0, &c->wt, &sd, &ix) != 6) exit(3); }
        else if (!strcmp(tok, "case")) {
            case_t *c = &CS[ncase++]; for (int i = 0; i < 6; i++) if (fscanf(f, "%lld", &c->a[i]) != 1) exit(3);
            for (int i = 0; i < 5; i++) if (fscanf(f, "%lld", &c->suf[i]) != 1) exit(3);
            if (fscanf(f, "%lld %lld %d %d", &c->lo, &c->hi, &c->evals, &c->nleaf) != 4) exit(3);
            c->leaf = malloc(sizeof(leaf_t) * c->nleaf);
            for (int j = 0; j < c->nleaf; j++) {
                char t2[16], r[64]; leaf_t *L = &c->leaf[j];
                if (fscanf(f, "%15s %lld %lld %63s", t2, &L->lo, &L->hi, r) != 4) exit(3);
                if (!strcmp(r, "R")) L->rejected = 1;
                else { L->rejected = 0; L->L.lo[0] = atoi(r); if (fscanf(f, "%d %d %d", &L->L.hi[0], &L->L.lo[1], &L->L.hi[1]) != 3) exit(3); }
            }
        }
        else if (!strcmp(tok, "end")) { char s[64]; if (fscanf(f, "%63s %lf", s, &ref_density) != 2) exit(3); }
    }
    fclose(f);
}
/* explore L0 W: the aligned L0-block of the full row containing the sample's b0, acceptance width W; slice-weighted */
static int DUMP;
static int explore(int L0, int W, int benchreps) {
    static leaf_t out[8192]; WIDTH = W; double sw = 0, se = 0, sl = 0, srej = 0, swidth = 0; long long pts = 0, rejpts = 0;
    for (int r = 0; r < (benchreps ? benchreps : 1); r++)
    for (int i = 0; i < ncase; i++) {
        case_t *c = &CS[i]; rowctx_t R; setup_a(c->a); row_setup(c->suf, &R);
        long long blo = c->b0 >= 0 ? c->b0 / L0 * L0 : -((-c->b0 + L0 - 1) / L0) * L0;   /* floor(b0 / L0) * L0 */
        long long lo = blo > c->rlo ? blo : c->rlo, hi = blo + L0 - 1 < c->rhi ? blo + L0 - 1 : c->rhi;
        long long e0 = EVALS; int n = adaptive_block(&R, lo, hi, out); double np = (double)(hi - lo + 1);
        if (r) continue;
        sw += c->wt; se += c->wt * (EVALS - e0) / np; sl += c->wt * n / np; pts += hi - lo + 1;
        if (DUMP) {   /* per pool row: evaluations, leaves, sample rejected, bounds at the sample (yield mask input) */
            int rej = 0; logs_t at = {{0, 0}, {0, 0}, {0, 0}};
            for (int j = 0; j < n; j++) if (out[j].lo <= c->b0 && c->b0 <= out[j].hi) { rej = out[j].rejected; at = out[j].L; }
            printf("row %d %lld %lld %lld %d %d %d %d %d %d %d\n", i, lo, hi, EVALS - e0, n, rej, at.lo[0], at.hi[0], at.lo[1], at.hi[1], nonfinite);
        }
        for (int j = 0; j < n; j++) {
            if (out[j].rejected) { rejpts += out[j].hi - out[j].lo + 1; if (out[j].lo <= c->b0 && c->b0 <= out[j].hi) srej += c->wt; }
            else if (out[j].lo <= c->b0 && c->b0 <= out[j].hi) swidth += c->wt * (out[j].L.hi[0] - out[j].L.lo[0] + out[j].L.hi[1] - out[j].L.lo[1]) / 16.0;
        }
    }
    printf("evaluations_total=%lld ", EVALS);
    printf("L0=%d W=%d weighted: evaluations/point=%.5f leaves/point=%.5f sampled-point rejected weight=%.6f mean (f+g) bound width at the sample=%.3f bits | points=%lld rejected points=%lld nonfinite=%d\n",
           L0, W, se / sw, sl / sw, srej / sw, swidth / sw, pts, rejpts, nonfinite);
    return 0;
}
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: nbound12 check|bench|explore vectors [reps | L0 W [benchreps]]\n"); return 2; }
    load(argv[2]);
    if (!strcmp(argv[1], "dump")) { DUMP = 1; return explore(atoi(argv[3]), atoi(argv[4]), 0); }
    if (!strcmp(argv[1], "explore")) return explore(atoi(argv[3]), atoi(argv[4]), argc > 5 ? atoi(argv[5]) : 0);
    if (!strcmp(argv[1], "benchrow")) {   /* per-a setup once per case, then row_setup reps times */
        int reps = atoi(argv[3]); long long rows = 0; rowctx_t R;
        for (int i = 0; i < ncase; i++) { setup_a(CS[i].a); for (int r = 0; r < reps; r++) { row_setup(CS[i].suf, &R); rows++; __asm__ __volatile__("" :: "r"(&R) : "memory"); } }
        printf("rows=%lld\n", rows); return 0;
    }
    int bench = !strncmp(argv[1], "bench", 5); int reps = bench && argc > 3 ? atoi(argv[3]) : 1;
    int setups_only = !strcmp(argv[1], "benchsetup");   /* per-a and per-row setups only (difference = adaptive blocks) */
    leaf_t out[64]; long long bad = 0, points = 0, leaves = 0, rejected = 0, rows = 0;
    for (int r = 0; r < reps; r++)
        for (int i = 0; i < ncase; i++) {
            case_t *c = &CS[i]; rowctx_t R; nonfinite = 0;
            setup_a(c->a); row_setup(c->suf, &R); rows++;
            if (setups_only) { points += c->hi - c->lo + 1; __asm__ __volatile__("" :: "r"(&R) : "memory"); continue; }
            long long e0 = EVALS; int n = adaptive_block(&R, c->lo, c->hi, out);
            points += c->hi - c->lo + 1; leaves += n;
            if (r) continue;
            int ok = n == c->nleaf && EVALS - e0 == c->evals && !nonfinite;
            for (int j = 0; ok && j < n; j++) {
                leaf_t *x = &out[j], *y = &c->leaf[j];
                ok = x->lo == y->lo && x->hi == y->hi && x->rejected == y->rejected;
                if (ok && !x->rejected) ok = x->L.lo[0] == y->L.lo[0] && x->L.hi[0] == y->L.hi[0] && x->L.lo[1] == y->L.lo[1] && x->L.hi[1] == y->L.hi[1];
                rejected += x->rejected;
            }
            if (!ok) { if (bad < 5) fprintf(stderr, "mismatch case %d\n", i); bad++; }
        }
    printf("cases=%d reps=%d rows=%lld points=%lld evaluations=%lld leaves=%lld rejected=%lld mismatches=%lld ref_weighted_density=%.12f\n",
           ncase, reps, rows, points, EVALS, leaves, rejected, bad, ref_density);
    return bad != 0;
}
