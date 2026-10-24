/* census_v10o.c — FACTORIAL attribution census for the p256 guarded-transition plateau.
 *
 * v10n shipped ONE guarded arm (prev-xor + ring4 + escape) against ONE control compiled
 * with -DUNGUARDED, whose index was mix64(x ^ SALT) with FIXED SALT.  census_v10f41.c,
 * the source of the nine pre-guard plateau cells, uses the PARITY-SELECTED salts
 * SA=0xa5a5.. / SB=0x5a5a.. (widx(x,par)).  The verdict names that mismatch, so the
 * v10n control is not the pre-guard walk and its 0.142 "collapse" is unattributed
 * between (a) the salt rule, (b) the carried prev row, (c) the ring/escape.
 *
 * This file makes the three mechanisms INDEPENDENT compile-time switches so a full
 * 2x2x2 attribution can be run on one binary build grammar and one seed policy:
 *   ARM_SALT_PARITY : index salt = (par ? SB : SA)   [census_v10f41.c rule]
 *                     else         = fixed SALT       [v10n -DUNGUARDED rule]
 *   ARM_PREV        : row = (mix64(x^salt) ^ prev) & 2047, prev carried
 *                     else  row = mix64(x^salt) & 2047          [class-only]
 *   ARM_RING        : 4-deep x ring, x(i)==x(i-4) => out-of-line escape (row 1, resalt)
 *   RING_INIT_ZERO  : ring seeded with the start x (v10n behaviour); default is a
 *                     sentinel so the first four comparisons are not self-matches.
 *
 * Everything else — fixture primes, the 2^40 class count, the two 2048-row tables, the
 * DP test mix64(x^DPSALT)&mask, the queue/merge/recovery solve — is byte-identical to
 * census_guarded_v10n.c / census_v10f41.c.
 *
 * Usage: census_v10o <dp_bits> <trials> <seed> <lam_num> <lam_den>
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;

#ifndef FP
#define FP 1099511629627ULL
#endif
#ifndef FQ
#define FQ 1099513640503ULL
#endif
#ifndef NCLASS
#define NCLASS 549756820252ULL
#endif
#ifndef CB
#define CB 2ULL
#endif
#ifndef GX
#define GX 2ULL
#endif
#ifndef GY
#define GY 185782196093ULL
#endif
#define SLOTS 16
#define NB 2048
#define DPSALT  0xd1b54a32d192ed03ULL
#define FIXSALT 0xd1b54a32d192ed03ULL
#define SA 0xa5a5a5a5a5a5a5a5ULL
#define SB 0x5a5a5a5a5a5a5a5aULL
#define RESALT 0x5a5a5a5a5a5a5a5aULL

static inline uint64_t mm(uint64_t a, uint64_t b) { return (uint64_t)((u128)a * b % FP); }
static inline uint64_t madd(uint64_t a, uint64_t b) { uint64_t s = a + b; return s >= FP ? s - FP : s; }
static inline uint64_t msub(uint64_t a, uint64_t b) { return a >= b ? a - b : a + FP - b; }
static inline uint64_t rmul(uint64_t a, uint64_t b) { return (uint64_t)((u128)a * b % FQ); }
static inline uint64_t radd(uint64_t a, uint64_t b) { uint64_t s = a + b; return s >= FQ ? s - FQ : s; }
static inline uint64_t rsub(uint64_t a, uint64_t b) { return a >= b ? a - b : a + FQ - b; }

static uint64_t mpow(uint64_t b, uint64_t e, uint64_t m)
{ uint64_t r = 1 % m; b %= m; while (e) { if (e & 1) r = (uint64_t)((u128)r * b % m); b = (uint64_t)((u128)b * b % m); e >>= 1; } return r; }
static inline uint64_t minv(uint64_t a) { return mpow(a, FP - 2, FP); }
static inline uint64_t rinv(uint64_t a) { return mpow(a, FQ - 2, FQ); }

/* curve a-coefficient: 0 = the v10o fixture curve y^2=x^3+b; FP-3 = the a=-3
 * doubling numerator 3x^2+a (the only place `a` enters the walk transition). */
static uint64_t CAV = 0;

typedef struct { uint64_t x, y; int inf; } pt;

static pt padd(pt A, pt B)
{
    if (A.inf) return B; if (B.inf) return A;
    if (A.x == B.x) {
        if (A.y != B.y || A.y == 0) { pt O = {0, 0, 1}; return O; }
        uint64_t l = mm(madd(mm(3, mm(A.x, A.x)), CAV), minv(madd(A.y, A.y)));
        uint64_t x3 = msub(mm(l, l), madd(A.x, A.x));
        uint64_t y3 = msub(mm(l, msub(A.x, x3)), A.y);
        pt R = {x3, y3, 0}; return R;
    }
    uint64_t l = mm(msub(B.y, A.y), minv(msub(B.x, A.x)));
    uint64_t x3 = msub(msub(mm(l, l), A.x), B.x);
    uint64_t y3 = mm(l, msub(A.x, x3)); y3 = msub(y3, A.y);
    pt R = {x3, y3, 0}; return R;
}

static pt smul(pt P, uint64_t k)
{
    pt R = {0, 0, 1};
    while (k) { if (k & 1) R = padd(R, P); P = padd(P, P); k >>= 1; }
    return R;
}

static inline uint64_t mix64(uint64_t v)
{ v *= 0x9e3779b97f4a7c15ULL; v ^= v >> 29; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 32; return v; }

static inline uint64_t splitmix64(uint64_t *s)
{ uint64_t z = (*s += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }

static pt     TBL[NB], TBL2[NB];
static int64_t TLOG[NB], TLOG2[NB];
static inline const pt *wrow(uint64_t j, int par) { return par ? &TBL2[j] : &TBL[j]; }
static inline int64_t rowlog(uint64_t j, int par) { return par ? TLOG2[j] : TLOG[j]; }

static int is_prime_u64(uint64_t n)
{ if (n < 2) return 0; if (n % 2 == 0) return n == 2;
  for (uint64_t d = 3; d * d <= n; d += 2) if (n % d == 0) return 0; return 1; }
static int on_curve(pt A)
{ if (A.inf) return 1; return mm(A.y, A.y) == madd(madd(mm(mm(A.x, A.x), A.x), mm(CAV, A.x)), CB); }

typedef struct { uint64_t key; uint64_t a, c; int owner; } rec;
static rec *HT;
static unsigned HBITS = 13, HSIZE;
static rec Q_[4096];
static int qn = 0;
static inline uint64_t hsh(uint64_t k, unsigned shift)
{ return (k * 0x9e3779b97f4a7c15ULL) >> (64 - shift); }

typedef struct { uint64_t x, y, a, c; uint64_t prev; uint64_t ring[4]; uint64_t ridx; } slot;

#ifndef ARM_NAME
#define ARM_NAME "guarded"
#endif
#ifndef PAIR_REPS
#define PAIR_REPS 1
#endif

static long g_ring_hits = 0, g_ring_esc = 0, g_tblx_esc = 0, g_dbl = 0;

int main(int argc, char **argv)
{
    if (argc < 6) { fprintf(stderr, "usage: %s dp_bits trials seed lam_num lam_den\n", argv[0]); return 2; }
    long cav_env = 0;
    if (getenv("CURVE_A")) cav_env = atol(getenv("CURVE_A"));
    /* applied AFTER fixture validation: the a=-3 probe keeps the a=0 curve and changes only the
     * doubling numerator, which is the sole place `a` enters the transition. */
    int b = atoi(argv[1]);
    long trials = atol(argv[2]);
    uint64_t seed0 = strtoull(argv[3], 0, 0);
    uint64_t lam_n = strtoull(argv[4], 0, 0), lam_d = strtoull(argv[5], 0, 0);

    int okp = is_prime_u64(FP), okr = is_prime_u64(FQ);
    pt G = {GX, GY, 0};
    int okorder = smul(G, FQ).inf;
    int okcurve = on_curve(G);
    if (cav_env) CAV = (cav_env < 0) ? (uint64_t)(FP - (uint64_t)(-cav_env) % FP) % FP : (uint64_t)cav_env % FP;
    if (!okp || !okr || !okorder || !okcurve) {
        fprintf(stderr, "FIXTURE INVALID p%d r%d ord%d curve%d\n", okp, okr, okorder, okcurve);
        return 3;
    }

    uint64_t s = 0xabcdef0123456789ULL;
    for (int j = 0; j < NB; j++) {
        uint64_t t; int dup;
        do { t = splitmix64(&s) % (FQ - 1) + 1; dup = 0;
             for (int k = 0; k < j; k++) if ((uint64_t)(TLOG[k] < 0 ? FQ - (uint64_t)(-TLOG[k]) : (uint64_t)TLOG[k]) == t) { dup = 1; break; }
        } while (dup);
        pt T = smul(G, t);
        if (T.y & 1) { T.y = FP - T.y; t = FQ - t; }
        TBL[j] = T; TLOG[j] = (int64_t)t;
    }
    {
        uint64_t s2 = 0x13579bdf2468ace0ULL;
        for (int j = 0; j < NB; j++) {
            uint64_t t; int dup;
            do { t = splitmix64(&s2) % (FQ - 1) + 1; dup = 0;
                 for (int k = 0; k < NB; k++) {
                     uint64_t u = (uint64_t)(TLOG[k] < 0 ? FQ - (uint64_t)(-TLOG[k]) : (uint64_t)TLOG[k]);
                     if (u == t || FQ - u == t) { dup = 1; break; }
                 }
            } while (dup);
            pt T = smul(G, t);
            if (T.y & 1) { T.y = FP - T.y; t = FQ - t; }
            TBL2[j] = T; TLOG2[j] = (int64_t)t;
        }
    }

    double lam = (double)lam_n / (double)lam_d;
    long Ttot = (long)(0.5 + sqrt(2.0 * (double)NCLASS * lam));
    long per = Ttot / SLOTS; if (per < 1) per = 1;
    if (getenv("CENSUS_STEPS")) per = atol(getenv("CENSUS_STEPS"));
    uint64_t mask = (b >= 63) ? ~0ULL : ((1ULL << b) - 1);

    double recmax = (double)per * SLOTS / (double)(1ULL << (b >= 40 ? 40 : b));
    HBITS = 13;
    while (HBITS < 21 && (double)(1u << HBITS) < 1.3 * recmax) HBITS++;
    HSIZE = 1u << HBITS;
    HT = calloc(HSIZE, sizeof(rec));
    if (!HT) { fprintf(stderr, "oom\n"); return 4; }

    uint64_t rng = seed0 ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1));
    long acc_dc0 = 0, acc_fail = 0, acc_rec = 0, acc_mer = 0, acc_self = 0, acc_ver = 0, acc_suc = 0;
    long acc_suc1 = 0, acc_run2 = 0;

    for (long tr = 0; tr < trials; tr++) {
        uint64_t x = splitmix64(&rng) % (FQ - 1) + 1;
        pt Qp = smul(G, x);
        if (Qp.inf) continue;
        int got = 0, got1 = 0;
        long nrec = 0, nmer = 0, nself = 0, nver = 0;
        for (int rep = 0; rep < PAIR_REPS; rep++) {
        uint64_t rngc = (PAIR_REPS > 1)
            ? (seed0 * 0x9e3779b97f4a7c15ULL + (uint64_t)tr * 0xbf58476d1ce4e5b9ULL
               + (uint64_t)rep * 0x94d049bb133111ebULL + 0x51ed270bULL)
            : rng;
        slot S[SLOTS];
        for (int i = 0; i < SLOTS; i++) {
            uint64_t ci = (uint64_t)(i + 1), ai; pt St;
            for (;;) {
                ai = splitmix64(&rngc) % (FQ - 1) + 1;
                uint64_t k0 = radd(ai, rmul(ci, x));
                St = smul(G, k0);
                if (!St.inf) break;
            }
            if (St.y & 1) { St.y = FP - St.y; ai = rsub(0, ai); ci = rsub(0, ci); }
            S[i].x = St.x; S[i].y = St.y; S[i].a = ai; S[i].c = ci;
            S[i].prev = 0;
#ifdef RING_INIT_ZERO
            for (int k = 0; k < 4; k++) S[i].ring[k] = St.x;
#else
            for (int k = 0; k < 4; k++) S[i].ring[k] = (uint64_t)-1;
#endif
            S[i].ridx = 0;
        }
        memset(HT, 0, (size_t)HSIZE * sizeof(rec));
        qn = 0;
        if (PAIR_REPS > 1 && got) break;   /* run 2 only if run 1 failed */

        for (long st = 0; st < per; st++) {
            uint64_t num[SLOTS], den[SLOTS], pre[SLOTS + 1];
            int esc[SLOTS], dbl[SLOTS];
            uint64_t J[SLOTS];
            for (int i = 0; i < SLOTS; i++) {
                int par = (int)(st & 1);
#ifdef ARM_SALT_PARITY
                uint64_t c = mix64(S[i].x ^ (par ? SB : SA));
#else
                uint64_t c = mix64(S[i].x ^ FIXSALT);
#endif
#ifdef ARM_PREV
                uint64_t j = (c ^ S[i].prev) & 2047;
                S[i].prev = j;
#else
                uint64_t j = c & 2047;
                (void)par;
#endif
                J[i] = j;
            }
            for (int i = 0; i < SLOTS; i++) {
                const pt *Tj = wrow(J[i], (int)(st & 1));
                dbl[i] = 0; esc[i] = 0;
                if (S[i].x == Tj->x) {
                    dbl[i] = 1; g_dbl++;
                    if (S[i].y != Tj->y) { esc[i] = 1; g_tblx_esc++; }
                    num[i] = madd(mm(3, mm(S[i].x, S[i].x)), 0);
                    den[i] = madd(S[i].y, S[i].y);
                } else {
                    num[i] = msub(Tj->y, S[i].y);
                    den[i] = msub(Tj->x, S[i].x);
                }
                if (den[i] == 0) den[i] = 1;
            }
            pre[0] = 1;
            for (int i = 0; i < SLOTS; i++) pre[i + 1] = mm(pre[i], den[i]);
            uint64_t inv = minv(pre[SLOTS]);
            uint64_t lam_[SLOTS];
            for (int i = SLOTS - 1; i >= 0; i--) {
                lam_[i] = mm(num[i], mm(inv, pre[i]));
                inv = mm(inv, den[i]);
            }
            for (int i = 0; i < SLOTS; i++) {
                const pt *Tj = wrow(J[i], (int)(st & 1));
                uint64_t x1 = S[i].x, y1 = S[i].y, na = S[i].a, nc = S[i].c;
                uint64_t x3, y3;
                if (esc[i]) {
                    pt Sp = padd((pt){x1, y1, 0}, (pt){x1, y1, 0});
                    x3 = Sp.x; y3 = Sp.y; na = radd(na, na); nc = radd(nc, nc);
                } else {
                    uint64_t l = lam_[i];
                    x3 = msub(msub(mm(l, l), x1), Tj->x);
                    y3 = msub(mm(l, msub(x1, x3)), y1);
                    if (dbl[i]) { na = radd(na, na); nc = radd(nc, nc); }
                    else        { na = radd(na, (uint64_t)rowlog(J[i], (int)(st & 1))); }
                }
                if (y3 & 1) { y3 = FP - y3; na = rsub(0, na); nc = rsub(0, nc); }
                int do_escape = 0;
#ifdef ARM_RING
                {   uint64_t k = S[i].ridx & 3;
                    if (x3 == S[i].ring[k]) { do_escape = 1; g_ring_hits++; }
                    S[i].ring[k] = x3;
                    S[i].ridx++;
                }
#endif
                S[i].x = x3; S[i].y = y3; S[i].a = na; S[i].c = nc;
                if (do_escape) {
                    g_ring_esc++;
                    pt E = TBL[1];
                    uint64_t ae = S[i].a, ce = S[i].c;
                    pt Sp;
                    if (S[i].x == E.x) {
                        Sp = padd((pt){S[i].x, S[i].y, 0}, (pt){S[i].x, S[i].y, 0});
                    } else {
                        uint64_t l = mm(msub(E.y, S[i].y), minv(msub(E.x, S[i].x)));
                        uint64_t xe = msub(msub(mm(l, l), S[i].x), E.x);
                        uint64_t ye = msub(mm(l, msub(S[i].x, xe)), S[i].y);
                        Sp.x = xe; Sp.y = ye; Sp.inf = 0;
                    }
                    ae = radd(ae, (uint64_t)TLOG[1]);
                    if (Sp.y & 1) { Sp.y = FP - Sp.y; ae = rsub(0, ae); ce = rsub(0, ce); }
                    S[i].x = Sp.x; S[i].y = Sp.y; S[i].a = ae; S[i].c = ce;
                    { uint64_t k = (S[i].ridx - 1) & 3; S[i].ring[k] = Sp.x; }
                    S[i].prev ^= RESALT;
                }
                if (((uint64_t)mix64(x3 ^ DPSALT) & mask) == 0) {
                    if (qn < 4096) { Q_[qn].key = x3; Q_[qn].a = na; Q_[qn].c = nc; Q_[qn].owner = i; qn++; }
                    nrec++;
                }
            }
            if (qn >= 32 || (st + 1 == per && qn)) {
                for (int q = 0; q < qn; q++) {
                    rec r = Q_[q];
                    uint64_t k = r.key + 1;
                    unsigned h = (unsigned)hsh(k, HBITS);
                    for (unsigned probe = 0; probe < HSIZE; probe++) {
                        unsigned idx = (h + probe) & (HSIZE - 1);
                        if (HT[idx].key == 0) { HT[idx].key = k; HT[idx].a = r.a; HT[idx].c = r.c; HT[idx].owner = r.owner; break; }
                        if (HT[idx].key == k) {
                            if (HT[idx].owner != r.owner) {
                                nmer++;
                                uint64_t dc = rsub(HT[idx].c, r.c);
                                if (dc) {
                                    uint64_t xs = rmul(rsub(r.a, HT[idx].a), rinv(dc));
                                    pt chk = smul(G, xs);
                                    if (!chk.inf && chk.x == Qp.x && chk.y == Qp.y) { nver++; got = 1; }
                                    else acc_fail++;
                                } else acc_dc0++;
                            } else nself++;
                            break;
                        }
                    }
                }
                qn = 0;
            }
        }
        if (rep == 0 && got) got1 = 1;
        if (rep == 1) acc_run2 += got;
        }  /* rep: second run uses FRESH walker coins on the SAME challenge x */
        acc_rec += nrec; acc_mer += nmer; acc_self += nself; acc_ver += nver;
        acc_suc += got; acc_suc1 += got1;
    }

    double ltheta = (double)per / (double)(1ULL << b);
    printf("{\"arm\":\"%s\",\"dp_bits\":%d,\"trials\":%ld,\"successes\":%ld,\"p_hat\":%.6f,"
           "\"lambda\":%.6f,\"ltheta\":%.6f,\"per_walker_steps\":%ld,"
           "\"records\":%ld,\"meets\":%ld,\"self\":%ld,\"verified\":%ld,"
           "\"dc0\":%ld,\"fail\":%ld,\"ring_hits\":%ld,\"ring_escapes\":%ld,\"tblx_escapes\":%ld,\"doubles\":%ld,"
           "\"fixture_classes_log2\":%.6f,\"run1_solved\":%ld,\"run2_solved\":%ld,\"pair_reps\":%d,\"curve_a\":%llu}\n",
           ARM_NAME, b, trials, acc_suc, (double)acc_suc / (double)trials,
           lam, ltheta, per, acc_rec, acc_mer, acc_self, acc_ver,
           acc_dc0, acc_fail, g_ring_hits, g_ring_esc, g_tblx_esc, g_dbl,
           log2((double)NCLASS) + 1.0, acc_suc1, acc_run2, PAIR_REPS, (unsigned long long)CAV);
    return 0;
}