/* eval17: certified per-block norm-bound evaluation, g-first side split + interval-Horner
 * (stnfs-sub123-v1, relation-collection eval lever; replaces the shipped nbound12 evaluator).
 *
 * WHAT IS CERTIFIED (identical to the shipped scheme, embedding_bounds.py / nbound12.c):
 *   For a candidate relation (a, b) in O_h^2 (a, b in Z[y], deg <= 5; K_h = Q[y]/(h),
 *   h = y^6+y^3+y^2-y-1) the two norm forms are
 *       N_f(a,b) = N_{K_h/Q}( b(y)^12 f(a(y)/b(y); y) ),   f = x^12 - 2y x^10 + 2y^3 x^6
 *                                                          + y^5 x^2 + (1 + y - y^2 - y^3),
 *       N_g(a,b) = N_{K_h/Q}( a(y)^2 + |u| y b(y)^2 ),      |u| = 0xd201000000010000.
 *   f is even in x, so per place v of K_h the f-form value is the sextic in w = b(v)^2
 *       F_v = c6 w^6 + c5 w^5 + c3 w^3 + c1 w + c0,
 *       c0 = a(v)^12, c1 = -2 v a(v)^10, c3 = 2 v^3 a(v)^6, c5 = v^5 a(v)^2, c6 = v^6,
 *   and the g-form value is linear in w: G_v = a(v)^2 + |u| v w.  Over the six embeddings
 *   these are products of 72 (f) resp. 12 (g) linear factors a(eta) - xi b(eta); the evaluation
 *   uses the four certified places (two real roots of h, two complex conjugate pairs counted
 *   through |.|^2) of the Rouché root disks, exactly as shipped.
 *   A block is an aligned interval of b_0 cells (production: 64 cells, b_1..b_5 fixed by the
 *   row); the block box is b_0 in [lo, hi] mapped through each place.  Every operation of the
 *   block evaluation is an outward-rounded binary64 enclosure, so per side
 *       16 log2 |N_side| in [LO, HI]   (integer sixteenth-bit bounds; floor/ceil to 2^(j/16))
 *   rigorously for every cell of the block.  A block is ACCEPTED as a leaf iff both sides are
 *   finite and at most WIDTH sixteenth-bits wide (production WIDTH = 128 = 8 bits: the width
 *   bounds survivors' cofactors by slack + 8 bits); otherwise it is split at the floor midpoint;
 *   an unresolved SINGLE POINT is REJECTED (never a survivor).  Accepted leaves carry the
 *   thresholds T_side = LO_side - 16*slack used by the sieve (proof.md section 3): every
 *   B-smooth cell has sieve mass >= T_side, so no smooth cell is lost (the conservative end of
 *   the certified interval - a lower bound on 16 log2|N| - is what sets the threshold; the
 *   evaluator never drops a cell because a norm UPPER bound is large).
 *
 * WHAT CHANGES (the lever; survivor sets stay bit-identical - see eval17_validate.py):
 *   1. g-first side split with early aborts.  A leaf needs BOTH sides certified, so a block
 *      whose g-side already fails certification (unboundable magnitude, or certified width
 *      > WIDTH) can only be split: the f-side evaluation is skipped entirely.  The width sums
 *      are monotone (each place contributes a nonnegative number of sixteenth-bits), so
 *      failure is decided as early as per place.  No cell is ever dropped by the abort: the
 *      block is split and re-resolved; the only dropped cells are unresolvable singletons,
 *      exactly the cells the shipped scheme drops.  Measured on all 2,401 pool rows
 *      (production 64/128): 15.3% of block evaluations abort on the g-side.
 *   2. Interval-Horner with explicit ulp-widening.  The block evaluation keeps the shipped
 *      Horner association EXACTLY (fv = fv*w; then +c_k, k = 5..0, the zero c2/c4 skipped -
 *      reassociating the sparse sextic changes the enclosure at interval-width level, not the
 *      ulp level, and is NOT bit-identical; a first design that tried it failed the check and
 *      was reverted).  What changes is only the outward-rounding discipline of the block
 *      evaluation: every operation rounds outward by at least one ulp via  x -+ |x| * 2^-51
 *      (wdn/wup, proved below) instead of the nextafter bit twiddle.  The two disciplines
 *      differ by at most ~2 ulps per operation, and the 16 log2 binning grid is 2^(1/16)-coarse,
 *      so the difference cannot move a bin in practice; bit-identity of every leaf and bound is
 *      verified on the full pool.
 *   The per-slice (setup_a) and per-row (row_setup) contexts keep the shipped nextafter
 *   discipline verbatim (they are priced separately as C_NB_ROW / C_NB_SLICE).
 *
 * RIGOR OF wdn/wup (explicit ulp-widening; IEEE-754 binary64, round-to-nearest, no FENV access):
 *   For normal x with ulp(x) = 2^(e-52) and x in [2^(e-1), 2^e), |x| * 2^-51 lies in
 *   [ulp(x), 2 ulp(x)), so fl(|x| * 2^-51) >= ulp(x) * (1 - 2^-53) > ulp(x)/2 + ulp(x)/2^54,
 *   and  x - fl(|x|*2^-51)  is strictly below x - ulp(x)/2, hence rounds to <= x - ulp(x)
 *   = nextafter(x, -inf); likewise  x + fl(|x|*2^-51)  rounds to >= nextafter(x, +inf).
 *   Every enclosure is therefore widened by at least one ulp on each side (at most two), which
 *   also covers the <= 0.5 ulp rounding error of the enclosing operation itself.  At +-0 (or
 *   subnormal endpoints, which cannot arise before the magnitudes'  m.lo > 0  gate fails) the
 *   widening may be zero, and a zero endpoint fails the same  m.lo > 0  test in both
 *   disciplines, so acceptance decisions agree.  (2^-51 rather than 2^-52: at exact powers of
 *   two |x|*2^-52 is exactly half an ulp and ties-to-even can leave x unwidened - not rigorous.)
 *
 * Modes (vector file = gen_nbound_vectors.py format, plus an optional `width <n>` token):
 *   check <vectors>          replay every reference leaf and bound bit-identically; exit 1 on any
 *                            difference; prints weighted evaluation density and abort counters
 *   dump <vectors> L0 W      aligned L0-blocks of every row (production: L0=64, W=128): one row
 *                            line (evals, leaves, sample rejected, sample bounds) + leaf lines
 *   explore <vectors> L0 W [reps]  summary only (weighted evaluation density, abort counters)
 *   The acceptance width comes from the vector file (`width` token, default 32); dump/explore
 *   override it with W.  No instruction counts are claimed here: eval17_bench.c is the
 *   measurement harness for the later callgrind run (linux/arm64).
 */
#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#define NOINLINE __attribute__((noinline))

typedef struct { double lo, hi; } iv;
typedef struct { iv re, im; } civ;

/* ---------- outward rounding: two disciplines ---------- */
/* shipped discipline (verbatim nbound12.c): nextafter(x, -/+inf) on the bit pattern;
   used by setup_a / row_setup only (per-slice and per-row contexts, priced separately). */
static inline double dn(double x) {
    uint64_t u; memcpy(&u, &x, 8); uint64_t s = (uint64_t)((int64_t)u >> 63);
    uint64_t r = u - 1 + (s & 2); r = (u << 1) == 0 ? 0x8000000000000001ULL : r; memcpy(&x, &r, 8); return x;
}
static inline double up(double x) {
    uint64_t u; memcpy(&u, &x, 8); uint64_t s = (uint64_t)((int64_t)u >> 63);
    uint64_t r = u + 1 - (s & 2); r = (u << 1) == 0 ? 1ULL : r; memcpy(&x, &r, 8); return x;
}
/* block-evaluation discipline: explicit ulp-widening x -+ |x|*2^-51 (header proof: always
   >= 1 ulp, at most ~2 ulps - an ulp-level perturbation of the shipped nextafter widening,
   verified bit-identical on the full pool; 2^-52 would tie-round to zero widening at exact
   powers of two).  EPS = 2^-51 is exact in binary64; |x|*EPS is one rounded multiply. */
static const double EPS = 0x1p-51;
static inline double wdn(double x) { return x - fabs(x) * EPS; }
static inline double wup(double x) { return x + fabs(x) * EPS; }

/* ---------- interval arithmetic (block-evaluation discipline unless noted) ---------- */
static inline iv mk(double lo, double hi) { iv r = {lo, hi}; return r; }
static inline double mn2(double a, double b) { return b < a ? b : a; }
static inline double mx2(double a, double b) { return b > a ? b : a; }
static inline iv ivadd(iv a, iv b) { return mk(wdn(a.lo + b.lo), wup(a.hi + b.hi)); }
static inline iv ivsub(iv a, iv b) { return mk(wdn(a.lo - b.hi), wup(a.hi - b.lo)); }
static inline iv ivmul(iv a, iv b) {
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    return mk(wdn(mn2(mn2(p1, p2), mn2(p3, p4))), wup(mx2(mx2(p1, p2), mx2(p3, p4))));
}
static inline iv ivabs(iv a) {
    double l = fabs(a.lo), h = fabs(a.hi); iv r;
    r.lo = (a.lo <= 0 && 0 <= a.hi) ? 0.0 : mn2(l, h); r.hi = mx2(l, h); return r;
}
static inline iv ivsq(iv a) { iv b = ivabs(a); return mk(mx2(0.0, wdn(b.lo * b.lo)), wup(b.hi * b.hi)); }
static inline civ cmk(iv re, iv im) { civ r = {re, im}; return r; }
static inline civ cadd(civ a, civ b) { return cmk(ivadd(a.re, b.re), ivadd(a.im, b.im)); }
static inline civ cmul(civ a, civ b) {
    iv A = ivmul(a.re, b.re), B = ivmul(a.im, b.im), re = ivsub(A, B);
    iv C = ivmul(a.re, b.im), E = ivmul(a.im, b.re), im = ivadd(C, E);
    return cmk(re, im);
}
static inline iv cabs2(civ a) { iv v = ivadd(ivsq(a.re), ivsq(a.im)); v.lo = mx2(0.0, v.lo); return v; }
/* shipped-discipline helpers for the per-slice/per-row contexts (verbatim nbound12.c) */
static inline iv nmk(double lo, double hi) { iv r = {lo, hi}; return r; }
static inline iv nadd(iv a, iv b) { return nmk(dn(a.lo + b.lo), up(a.hi + b.hi)); }
static inline iv nmul(iv a, iv b) {
    double p1 = a.lo * b.lo, p2 = a.lo * b.hi, p3 = a.hi * b.lo, p4 = a.hi * b.hi;
    return nmk(dn(mn2(mn2(p1, p2), mn2(p3, p4))), up(mx2(mx2(p1, p2), mx2(p3, p4))));
}
static inline iv nrat(double c) { return nmk(dn(c), up(c)); }
static inline civ ncmk(iv re, iv im) { civ r = {re, im}; return r; }
static inline civ nconst(double c) { iv z = {0, 0}; return ncmk(nrat(c), z); }
static inline civ ncadd(civ a, civ b) { return ncmk(nadd(a.re, b.re), nadd(a.im, b.im)); }
static inline civ ncmul(civ a, civ b) {
    iv A = nmul(a.re, b.re), B = nmul(a.im, b.im), re = nadd(A, nmk(-B.hi, -B.lo));
    iv C = nmul(a.re, b.im), E = nmul(a.im, b.re), im = nadd(C, E);
    return ncmk(re, im);
}

/* ---------- 16 log2 bounds (verbatim nbound12.c; integer-exact) ---------- */
static double LL[17], LU[17];
static inline int count_le(const double *T, double x) { int n = 0; for (int k = 16; k; k >>= 1) if (n + k <= 17 && T[n + k - 1] <= x) n += k; return n; }
static inline int count_lt(const double *T, double x) { int n = 0; for (int k = 16; k; k >>= 1) if (n + k <= 17 && T[n + k - 1] < x) n += k; return n; }
static inline void frexp2(double x, double *m2, int *e) {
    uint64_t u; memcpy(&u, &x, 8); int be = (int)((u >> 52) & 0x7FF);
    if (be == 0) { int ee; double m = frexp(x, &ee); *m2 = 2 * m; *e = ee; return; }
    u = (u & 0x000FFFFFFFFFFFFFULL) | 0x3FF0000000000000ULL; memcpy(m2, &u, 8); *e = be - 1022;
}
static inline int log_bounds(iv m, int *lo, int *hi) {
    if (!(m.lo > 0) || !isfinite(m.lo) || !isfinite(m.hi) || m.hi < m.lo) return 0;
    double x; int e;
    frexp2(m.lo, &x, &e); *lo = 16 * (e - 1) + count_le(LU, x) - 1;
    frexp2(m.hi, &x, &e); *hi = 16 * (e - 1) + count_lt(LL, x);
    return 1;
}

/* ---------- per-a context (once per slice; shipped verbatim, nextafter discipline) ---------- */
#define NPL 4
typedef struct { int real; iv rr; civ cr; iv one_r; civ one_c; iv a2r, uyr, cr_[7]; civ a2c, uyc, cc[7]; } place_t;
static place_t PL[NPL];
static iv ROOT_RE[NPL], ROOT_IM[NPL]; static int ROOT_REAL[NPL];
static const double UABS = 15132376222941642752.0;   /* |u| rounds exactly to this double */
static int nonfinite;   /* set when an input is outside the guarded range (|coeff| < 2^12) */
static void setup_a(const long long *a) {
    for (int i = 0; i < 6; i++) if (a[i] <= -4096 || a[i] >= 4096) nonfinite = 1;
    for (int v = 0; v < NPL; v++) {
        place_t *P = &PL[v]; P->real = ROOT_REAL[v];
        if (P->real) {
            iv root = ROOT_RE[v], one = {1, 1}, av = nmul(one, nrat(0.0)); P->rr = root; P->one_r = one;
            for (int i = 5; i >= 0; i--) av = nadd(nmul(av, root), nrat((double)a[i]));
            iv ap[13], yp[7]; ap[0] = one; yp[0] = one;
            for (int k = 1; k <= 12; k++) ap[k] = nmul(ap[k - 1], av);
            for (int k = 1; k <= 6; k++) yp[k] = nmul(yp[k - 1], root);
            P->cr_[0] = ap[12]; P->cr_[1] = nmul(nmul(root, nrat(-2.0)), ap[10]); P->cr_[2] = nmul(one, nrat(0.0));
            P->cr_[3] = nmul(nmul(yp[3], nrat(2.0)), ap[6]); P->cr_[4] = nmul(one, nrat(0.0));
            P->cr_[5] = nmul(yp[5], ap[2]); P->cr_[6] = yp[6];
            P->a2r = ap[2]; P->uyr = nmul(root, nrat(UABS));
        } else {
            civ root = ncmk(ROOT_RE[v], ROOT_IM[v]), one = nconst(1.0), av = ncmul(one, nconst(0.0)); P->cr = root; P->one_c = one;
            for (int i = 5; i >= 0; i--) av = ncadd(ncmul(av, root), nconst((double)a[i]));
            civ ap[13], yp[7]; ap[0] = one; yp[0] = one;
            for (int k = 1; k <= 12; k++) ap[k] = ncmul(ap[k - 1], av);
            for (int k = 1; k <= 6; k++) yp[k] = ncmul(yp[k - 1], root);
            P->cc[0] = ap[12]; P->cc[1] = ncmul(ncmul(root, nconst(-2.0)), ap[10]); P->cc[2] = ncmul(one, nconst(0.0));
            P->cc[3] = ncmul(ncmul(yp[3], nconst(2.0)), ap[6]); P->cc[4] = ncmul(one, nconst(0.0));
            P->cc[5] = ncmul(yp[5], ap[2]); P->cc[6] = yp[6];
            P->a2c = ap[2]; P->uyc = ncmul(root, nconst(UABS));
        }
    }
}
/* ---------- per-row context (once per row; shipped verbatim, nextafter discipline) ---------- */
typedef struct { iv br[NPL]; civ bc[NPL]; } rowctx_t;
static NOINLINE void row_setup(const long long *suffix, rowctx_t *R) {
    for (int i = 0; i < 5; i++) if (suffix[i] <= -4096 || suffix[i] >= 4096) nonfinite = 1;
    for (int v = 0; v < NPL; v++) {
        place_t *P = &PL[v];
        if (P->real) { iv bv = nmul(P->one_r, nrat(0.0)); for (int i = 4; i >= 0; i--) bv = nadd(nmul(bv, P->rr), nrat((double)suffix[i])); R->br[v] = nmul(bv, P->rr); }
        else { civ bv = ncmul(P->one_c, nconst(0.0)); for (int i = 4; i >= 0; i--) bv = ncadd(ncmul(bv, P->cr), nconst((double)suffix[i])); R->bc[v] = ncmul(bv, P->cr); }
    }
}

/* ---------- one block evaluation: g-first side split (the lever) ----------
 * pass 1 computes the g-side bounds place by place and returns as soon as the g-side
 * certification is certain to fail (unboundable magnitude, or partial sixteenth-bit width
 * already > WIDTH; widths only grow, since each place adds a nonnegative amount);
 * pass 2 (f-side) runs only for g-passing blocks, with the same early-abort rules.
 * The emitted decision and bounds are identical to the shipped both-sides evaluation:
 * a block is accepted iff BOTH sides are finite and <= WIDTH wide, and every emitted bound
 * is the same per-side place-order sum. */
typedef struct { int ok, lo, hi; } side_t;
typedef struct { side_t f, g; } logs_t;
static long long EVALS;                    /* block evaluations (the priced quantity) */
static long long GABORT_OK[NPL], GABORT_WD[NPL], FABORT_OK[NPL], FABORT_WD[NPL];
static long long G_ABORTS, F_ABORTS, PASSES;
static int ABORTMODE;                      /* 0 = g-first with early aborts (design); 1 = no aborts (A/B); 2 = g-side only (pricing) */
/* acceptance width in sixteenth-bits (production 128 = 8 bits); set from the vector
   file / CLI before any block evaluation */
static int WIDTH_RUN = 32;
static NOINLINE void block_eval(const rowctx_t *R, long long lo, long long hi, logs_t *L) {
    EVALS++;
    L->g.ok = 1; L->g.lo = L->g.hi = 0; L->f.ok = 1; L->f.lo = L->f.hi = 0;
    iv b0 = {(double)lo, (double)hi}, w[NPL], z = {0, 0};
    civ cw[NPL];
    const int W = WIDTH_RUN;
    for (int v = 0; v < NPL; v++) {                 /* pass 1: g-side */
        const place_t *P = &PL[v]; iv mag; int l, h;
        if (P->real) {
            iv bv = ivadd(R->br[v], b0); w[v] = ivmul(bv, bv);
            iv gv = ivadd(P->a2r, ivmul(P->uyr, w[v]));
            mag = ivabs(gv);
        } else {
            civ bv = cadd(R->bc[v], cmk(b0, z)); cw[v] = cmul(bv, bv);
            civ gv = cadd(P->a2c, cmul(P->uyc, cw[v]));
            mag = cabs2(gv);
        }
        if (!log_bounds(mag, &l, &h)) { L->g.ok = 0; if (ABORTMODE == 0) { GABORT_OK[v]++; G_ABORTS++; return; } }
        else if (L->g.ok) { L->g.lo += l; L->g.hi += h;
            if (ABORTMODE == 0 && L->g.hi - L->g.lo > W) { GABORT_WD[v]++; G_ABORTS++; return; } }
    }
    if (ABORTMODE == 2) { PASSES++; return; }       /* pricing probe: g-side only */
    for (int v = 0; v < NPL; v++) {                 /* pass 2: f-side, the shipped Horner association
                                                       verbatim (fv = fv*w then +c_k, k=5..0, c2/c4
                                                       zero and skipped): reassociating the sparse
                                                       sextic changes the enclosure at interval-width
                                                       level and is NOT bit-identical - do not. */
        const place_t *P = &PL[v]; iv mag; int l, h;
        if (P->real) {
            iv fv = P->cr_[6];
            for (int k = 5; k >= 0; k--) { fv = ivmul(fv, w[v]); if (k != 2 && k != 4) fv = ivadd(fv, P->cr_[k]); }
            mag = ivabs(fv);
        } else {
            civ fv = P->cc[6];
            for (int k = 5; k >= 0; k--) { fv = cmul(fv, cw[v]); if (k != 2 && k != 4) fv = cadd(fv, P->cc[k]); }
            mag = cabs2(fv);
        }
        if (!log_bounds(mag, &l, &h)) { L->f.ok = 0; if (ABORTMODE == 0) { FABORT_OK[v]++; F_ABORTS++; return; } }
        else if (L->f.ok) { L->f.lo += l; L->f.hi += h;
            if (ABORTMODE == 0 && L->f.hi - L->f.lo > W) { FABORT_WD[v]++; F_ABORTS++; return; } }
    }
    PASSES++;
}

/* ---------- adaptive aligned block, singleton policy 'reject' (shipped recursion) ---------- */
typedef struct { long long lo, hi; int rejected; logs_t L; } leaf_t;
static NOINLINE int adaptive_block(const rowctx_t *R, long long lo, long long hi, leaf_t *out) {
    long long st[64][2]; int sp = 0, n = 0; st[sp][0] = lo; st[sp][1] = hi; sp++;
    while (sp) {
        sp--; long long l = st[sp][0], h = st[sp][1]; logs_t L; block_eval(R, l, h, &L);
        if (L.f.ok && L.g.ok && L.f.hi - L.f.lo <= WIDTH_RUN && L.g.hi - L.g.lo <= WIDTH_RUN) {
            out[n].lo = l; out[n].hi = h; out[n].rejected = 0; out[n].L = L; n++;
        } else if (l == h) { out[n].lo = l; out[n].hi = h; out[n].rejected = 1; n++; }
        else { long long m = (l + h) >> 1; st[sp][0] = m + 1; st[sp][1] = h; sp++; st[sp][0] = l; st[sp][1] = m; sp++; }
    }
    return n;
}

/* ---------- vector file (gen_nbound_vectors.py format + optional `width` token) ---------- */
typedef struct { long long a[6], suf[5], lo, hi, rlo, rhi, b0; double wt; int evals, nleaf; leaf_t *leaf; } case_t;
static case_t *CS; static int ncase; static double ref_density; static int file_width = 32;
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
        else if (!strcmp(tok, "width")) { if (fscanf(f, "%d", &file_width) != 1) exit(3); }
        else if (!strcmp(tok, "logtab")) { char s[64]; for (int j = 0; j < 17; j++) { if (fscanf(f, "%63s", s) != 1) exit(3); LL[j] = strtod(s, 0); } for (int j = 0; j < 17; j++) { if (fscanf(f, "%63s", s) != 1) exit(3); LU[j] = strtod(s, 0); } }
        else if (!strcmp(tok, "row")) { case_t *c = &CS[ncase]; int sd, ix;
            if (fscanf(f, "%lld %lld %lld %lf %d %d", &c->rlo, &c->rhi, &c->b0, &c->wt, &sd, &ix) != 6) exit(3); }
        else if (!strcmp(tok, "case")) {
            case_t *c = &CS[ncase++]; if (ncase > cap) { cap *= 2; CS = realloc(CS, sizeof(case_t) * cap); c = &CS[ncase - 1]; }
            for (int i = 0; i < 6; i++) if (fscanf(f, "%lld", &c->a[i]) != 1) exit(3);
            for (int i = 0; i < 5; i++) if (fscanf(f, "%lld", &c->suf[i]) != 1) exit(3);
            if (fscanf(f, "%lld %lld %d %d", &c->lo, &c->hi, &c->evals, &c->nleaf) != 4) exit(3);
            c->leaf = malloc(sizeof(leaf_t) * (c->nleaf ? c->nleaf : 1));
            for (int j = 0; j < c->nleaf; j++) {
                char t2[16], r[64]; leaf_t *L = &c->leaf[j];
                if (fscanf(f, "%15s %lld %lld %63s", t2, &L->lo, &L->hi, r) != 4) exit(3);
                if (!strcmp(r, "R")) L->rejected = 1;
                else { L->rejected = 0; L->L.f.lo = atoi(r); if (fscanf(f, "%d %d %d", &L->L.f.hi, &L->L.g.lo, &L->L.g.hi) != 3) exit(3); }
            }
        }
        else if (!strcmp(tok, "end")) { char s[64]; if (fscanf(f, "%63s %lf", s, &ref_density) != 2) exit(3); }
    }
    fclose(f);
}
static long long LEAVES_ACC, LEAVES_REJ, PTS_ACC, PTS_REJ;
static void tally_block(const leaf_t *out, int n) {
    for (int j = 0; j < n; j++) {
        LEAVES_ACC++;
        if (out[j].rejected) { LEAVES_REJ++; PTS_REJ += out[j].hi - out[j].lo + 1; }
        else PTS_ACC += out[j].hi - out[j].lo + 1;
    }
}
static void print_counters(void) {
    long long e = EVALS;
    printf("evaluations=%lld (passes=%lld g-aborts=%lld f-aborts=%lld) leaves=%lld rejected-leaves=%lld points accepted=%lld rejected=%lld",
           e, PASSES, G_ABORTS, F_ABORTS, LEAVES_ACC, LEAVES_REJ, PTS_ACC, PTS_REJ);
    printf(" g-abort by place: ok[0..3]=%lld,%lld,%lld,%lld width=%lld,%lld,%lld,%lld",
           GABORT_OK[0], GABORT_OK[1], GABORT_OK[2], GABORT_OK[3],
           GABORT_WD[0], GABORT_WD[1], GABORT_WD[2], GABORT_WD[3]);
    printf(" f-abort by place: ok[0..3]=%lld,%lld,%lld,%lld width=%lld,%lld,%lld,%lld\n",
           FABORT_OK[0], FABORT_OK[1], FABORT_OK[2], FABORT_OK[3],
           FABORT_WD[0], FABORT_WD[1], FABORT_WD[2], FABORT_WD[3]);
}
static leaf_t OUT[8192];
/* explore/dump: aligned L0-block of the row containing the sample b0, clipped to [rlo, rhi] */
static int DUMPMODE;
static int __attribute__((unused)) explore(int L0, int W, int reps) {
    WIDTH_RUN = W; double sw = 0, se = 0, sl = 0, srej = 0; long long pts = 0;
    for (int r = 0; r < (reps ? reps : 1); r++)
    for (int i = 0; i < ncase; i++) {
        case_t *c = &CS[i]; rowctx_t R; nonfinite = 0;
        setup_a(c->a); row_setup(c->suf, &R);
        long long blo = c->b0 >= 0 ? c->b0 / L0 * L0 : -((-c->b0 + L0 - 1) / L0) * L0;   /* floor(b0/L0)*L0 */
        long long lo = blo > c->rlo ? blo : c->rlo, hi = blo + L0 - 1 < c->rhi ? blo + L0 - 1 : c->rhi;
        long long e0 = EVALS; int n = adaptive_block(&R, lo, hi, OUT); double np = (double)(hi - lo + 1);
        if (r) continue;
        sw += c->wt; se += c->wt * (EVALS - e0) / np; sl += c->wt * n / np; pts += hi - lo + 1;
        tally_block(OUT, n);
        int rej = -1; logs_t at; memset(&at, 0, sizeof(at));
        for (int j = 0; j < n; j++) if (OUT[j].lo <= c->b0 && c->b0 <= OUT[j].hi) { rej = OUT[j].rejected; at = OUT[j].L; }
        if (DUMPMODE) {
            printf("row %d %lld %lld %lld %d %d %d %d %d %d %d\n", i, lo, hi, EVALS - e0, n, rej,
                   at.f.lo, at.f.hi, at.g.lo, at.g.hi, nonfinite);
            for (int j = 0; j < n; j++)
                if (OUT[j].rejected) printf(" leaf %lld %lld R\n", OUT[j].lo, OUT[j].hi);
                else printf(" leaf %lld %lld %d %d %d %d\n", OUT[j].lo, OUT[j].hi, OUT[j].L.f.lo, OUT[j].L.f.hi, OUT[j].L.g.lo, OUT[j].L.g.hi);
        }
        if (rej > 0) srej += c->wt;
    }
    printf("L0=%d W=%d weighted: evaluations/point=%.5f leaves/point=%.5f sampled-point rejected weight=%.6f points=%lld nonfinite=%d\n",
           L0, W, se / sw, sl / sw, srej / sw, pts, nonfinite);
    print_counters();
    printf("reference_weighted_density=%.12f\n", ref_density);
    return 0;
}
static int __attribute__((unused)) check(const char *path) {
    (void)path;
    WIDTH_RUN = file_width;
    long long bad = 0, points = 0, rows = 0;
    for (int i = 0; i < ncase; i++) {
        case_t *c = &CS[i]; rowctx_t R; nonfinite = 0;
        setup_a(c->a); row_setup(c->suf, &R); rows++;
        long long e0 = EVALS; int n = adaptive_block(&R, c->lo, c->hi, OUT);
        points += c->hi - c->lo + 1; tally_block(OUT, n);
        int ok = n == c->nleaf && EVALS - e0 == c->evals && !nonfinite;
        for (int j = 0; ok && j < n; j++) {
            leaf_t *x = &OUT[j], *y = &c->leaf[j];
            ok = x->lo == y->lo && x->hi == y->hi && x->rejected == y->rejected;
            if (ok && !x->rejected)
                ok = x->L.f.lo == y->L.f.lo && x->L.f.hi == y->L.f.hi && x->L.g.lo == y->L.g.lo && x->L.g.hi == y->L.g.hi;
        }
        if (!ok) { if (bad < 5) fprintf(stderr, "MISMATCH case %d\n", i); bad++; }
    }
    printf("check: cases=%d rows=%lld points=%lld mismatches=%lld width=%d ref_weighted_density=%.12f\n",
           ncase, rows, points, bad, file_width, ref_density);
    print_counters();
    return bad != 0;
}
#ifndef EVAL17_NO_MAIN
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: eval17 check|dump|explore <vectors> [L0 W [reps]]\n"); return 2; }
    load(argv[2]);
    if (!strcmp(argv[1], "check")) return check(argv[2]);
    if (argc < 5) { fprintf(stderr, "usage: eval17 dump|explore <vectors> L0 W [reps]\n"); return 2; }
    int L0 = atoi(argv[3]), W = atoi(argv[4]);
    if (!strcmp(argv[1], "dump")) { DUMPMODE = 1; return explore(L0, W, 0); }
    if (!strcmp(argv[1], "explore")) return explore(L0, W, argc > 5 ? atoi(argv[5]) : 0);
    fprintf(stderr, "unknown mode %s\n", argv[1]); return 2;
}
#endif
