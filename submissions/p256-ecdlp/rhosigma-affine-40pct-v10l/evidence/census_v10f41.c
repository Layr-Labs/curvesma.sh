/* census_v10f41.c — v10f part 4: ATTRACTOR-FREE recovery census on the 2^40 fixture.
 *
 * Defect fixed here (see HOLD-part4.md): with one shared table and the index taken
 * from the state, the canon sign-flip makes a repeated index an INVOLUTION --
 * p -> +-(p+T_j), and when idx repeats (birthday over 2048 rows, ~2^11 steps) the
 * second step cancels the first: p -> q -> p, a 2-cycle.  Measured exactly (bigint):
 * 100% of starts, median entry 1,600 steps, 40% of a 5,000-step trail is a repeat,
 * and it is INDEPENDENT of the index rule (x&2047, mix64(x)&2047, per-walker coset
 * offset and xor-fold all identical) -- which is why rescaling could never help.
 *
 * Fix (option O3): TWO shared tables, alternated by step parity, each indexed by a
 * salted mix of the canonical x.  A repeat then requires T^A_j == T^B_j' exactly
 * (a linear coincidence, ~2^-40) instead of an index collision (~2^-11).
 * Basin probe: 0/500 starts enter a cycle within 5,000 steps and 0/1.5M trail
 * steps repeat (baseline 62% / 40%).  Canon retained => the proven (a,c) solve
 * is untouched.  Storage: +2048 rows (2x table) = 2.6e-11 of W.
 *
 * census_v10f40.c — production-shaped recovery census on a 2^40-class prime-order fixture.
 *
 * RESCALE of census_v10f.c.  The 2^24-class fixture put 2,048 table rows against
 * 2^24 classes and the walk self-captured: verified recovery fell with DP
 * density (0.32/0.26/0.21/0.13) which is the cycle-domination signature.  This
 * fixture is the same construction with q ~ 2^40 (N' = (q+1)/2 ~ 2^39):
 *
 *   p = 1099511629627 (prime, 2^40.0),  E: y^2 = x^3 + 2
 *   #E = 1099513640503 = q  (prime; cofactor 1 -- CM, D=-3, 4p = t^2+3v^2)
 *   G = (2, 185782196093), [q]G = O, checked at startup
 *   N' = (q+1)/2 = 549756820252 ~ 2^39.0 classes
 *
 * Walk budget: T = round(sqrt(2 N' lambda)) total steps over 16 slots, i.e.
 * S_t = T/16 ~ 2^15.5 per walker, while a random map on N' classes first
 * revisits at ~sqrt(N') ~ 2^19.5: a ~16x per-walker margin (4x budget = the
 * deliberately-captured control arm, argv[6]).
 *
 * Everything else is the same mechanism as census_v10f.c: even-y canon with
 * (a,c) negated with it, records by mix64(x^DPSALT) & mask, queue drained by
 * hash insert, merges solved and VERIFIED by [x]P == Q.
 *
 * Arithmetic: field and scalar ops use compile-time constants so the modular
 * reductions become multiply-high sequences, and the 16 slot additions share
 * ONE inversion per step (the Montgomery batch trick -- the same shared
 * inversion the priced kernel performs over the walk slots).
 *
 * Usage: census_v10f40 <dp_bits> <trials> <seed> <lam_num> <lam_den> [budget_mult]
 */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;

#define FP 1099511629627ULL
#define FQ 1099513640503ULL
#define NCLASS 549756820252ULL
#define CB 2ULL
#define GX 2ULL
#define GY 185782196093ULL
#define SLOTS 16
#define NB 2048

static uint64_t SALT = 0xd1b54a32d192ed03ULL;

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

typedef struct { uint64_t x, y; int inf; } pt;

static pt padd(pt A, pt B)              /* full addition, infinity-safe */
{
    if (A.inf) return B;
    if (B.inf) return A;
    if (A.x == B.x) {
        if (A.y != B.y) { pt O = {0, 0, 1}; return O; }
        uint64_t num = mm(3, mm(A.x, A.x));
        uint64_t l = mm(num, minv(mm(2, A.y)));
        uint64_t x3 = msub(msub(mm(l, l), A.x), A.x);
        return (pt){x3, msub(mm(l, msub(A.x, x3)), A.y), 0};
    }
    uint64_t l = mm(msub(B.y, A.y), minv(msub(B.x, A.x)));
    uint64_t x3 = msub(msub(mm(l, l), A.x), B.x);
    return (pt){x3, msub(mm(l, msub(A.x, x3)), A.y), 0};
}

static pt smul(pt B, uint64_t k)
{
    pt Acc = {0, 0, 1};
    int top = 63;
    while (top > 0 && ((k >> top) & 1) == 0) top--;
    for (int i = top; i >= 0; i--) {
        Acc = padd(Acc, Acc);
        if ((k >> i) & 1) Acc = padd(Acc, B);
    }
    return Acc;
}

static inline uint64_t mix64(uint64_t v)
{ v *= 0x9e3779b97f4a7c15ULL; v ^= v >> 29; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 32; return v; }

static int IMODE = 1;   /* 1 = O3 two-table salted mix (default, attractor-free); 0 = part-3 baseline */
static const uint64_t SA = 0xa5a5a5a5a5a5a5a5ULL, SB = 0x5a5a5a5a5a5a5a5aULL;
static inline uint64_t widx(uint64_t x, int par)
{
    if (!IMODE) return x & (NB - 1);
    return par ? (mix64(x ^ SB) & (NB - 1)) : (mix64(x ^ SA) & (NB - 1));
}


static inline uint64_t splitmix64(uint64_t *s)
{ uint64_t z = (*s += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }

static pt    TBL[NB], TBL2[NB];
static int64_t TLOG[NB], TLOG2[NB];
static inline const pt *wrow(uint64_t j, int par) { return par ? &TBL2[j] : &TBL[j]; }

static int is_prime_u64(uint64_t n)
{ if (n < 2) return 0; if (n % 2 == 0) return n == 2;
  for (uint64_t d = 3; d * d <= n; d += 2) if (n % d == 0) return 0; return 1; }
static int on_curve(pt A)
{ if (A.inf) return 1; return mm(A.y, A.y) == (mm(mm(A.x, A.x), A.x) + CB) % FP; }

/* ---- per-trial store ---- */
typedef struct { uint64_t key; uint64_t a, c; int owner; } rec;
static rec *HT;
static unsigned HBITS = 13, HSIZE;
static rec Q_[4096];
static int qn = 0;

static inline uint64_t hsh(uint64_t k, unsigned shift)
{ return (k * 0x9e3779b97f4a7c15ULL) >> (64 - shift); }

typedef struct { uint64_t x, y; uint64_t a, c; } slot;

int main(int argc, char **argv)
{
    if (argc < 6) { fprintf(stderr, "usage: %s dp_bits trials seed lam_num lam_den [budget_mult]\n", argv[0]); return 2; }
    int b = atoi(argv[1]);
    long trials = atol(argv[2]);
    uint64_t seed0 = strtoull(argv[3], 0, 0);
    uint64_t lam_n = strtoull(argv[4], 0, 0), lam_d = strtoull(argv[5], 0, 0);
    long mult = argc > 6 ? atol(argv[6]) : 1;
    if (mult < 1) mult = 1;

    int okp = is_prime_u64(FP), okr = is_prime_u64(FQ);
    pt G = {GX, GY, 0};
    int okorder = smul(G, FQ).inf;
    int okcurve = on_curve(G);
    if (!okp || !okr || !okorder || !okcurve) {
        fprintf(stderr, "FIXTURE INVALID p%d r%d ord%d curve%d\n", okp, okr, okorder, okcurve);
        return 3;
    }

    /* deterministic known-log jump table, canonical (even-y) rows */
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
    {   /* O3: second table, independent stream, cross-table dup-checked */
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

    if (getenv("CENSUS_INVARIANT"))
        for (int j = 0; j < NB; j++) {
            if (TBL2[j].y & 1) { fprintf(stderr, "TABLE2 ROW %d ODD\n", j); return 8; }
            { pt v2 = smul(G, (uint64_t)(TLOG2[j] < 0 ? FQ - (uint64_t)(-TLOG2[j]) : (uint64_t)TLOG2[j]));
              if (v2.x != TBL2[j].x || v2.y != TBL2[j].y) { fprintf(stderr, "TABLE2 ROW %d MISMATCH\n", j); return 8; } }
            pt v = smul(G, (uint64_t)(TLOG[j] < 0 ? FQ - (uint64_t)(-TLOG[j]) : (uint64_t)TLOG[j]));
            int64_t logged = TLOG[j];
            pt vv = smul(G, (uint64_t)(logged < 0 ? FQ - (uint64_t)(-logged) : (uint64_t)logged));
            if (v.x != TBL[j].x || v.y != TBL[j].y) { fprintf(stderr, "TABLE ROW %d MISMATCH\n", j); return 8; }
            if (TBL[j].y & 1) { fprintf(stderr, "TABLE ROW %d ODD\n", j); return 8; }
            (void)v;
        }

    double lam = (double)lam_n / (double)lam_d;
    long Ttot = (long)(0.5 + sqrt(2.0 * (double)NCLASS * lam)) * mult;
    long per = Ttot / SLOTS; if (per < 1) per = 1;
    if (getenv("CENSUS_STEPS")) per = atol(getenv("CENSUS_STEPS"));
    uint64_t mask = (b >= 63) ? ~0ULL : ((1ULL << b) - 1);
    { const char *im = getenv("CENSUS_INDEX"); if (im && !strcmp(im, "low")) IMODE = 0; }

    if (getenv("CENSUS_ARITH")) {              /* arithmetic cross-check rows for bigint audit */
        long n = atol(getenv("CENSUS_ARITH"));
        uint64_t ss = 0x1234567890abcdefULL;
        for (long i = 0; i < n; i++) {
            pt A = smul(G, splitmix64(&ss) % (FQ - 1) + 1);
            pt B = smul(G, splitmix64(&ss) % (FQ - 1) + 1);
            if (A.y & 1) A.y = FP - A.y;
            if (B.y & 1) B.y = FP - B.y;
            pt C = padd(A, B);
            printf("%llu %llu %llu %llu %llu %llu\n", (unsigned long long)A.x, (unsigned long long)A.y,
                   (unsigned long long)B.x, (unsigned long long)B.y,
                   (unsigned long long)C.x, (unsigned long long)C.y);
        }
        return 0;
    }

    /* hash slots ~1.3x the largest possible record count for this density */
    double recmax = (double)per * SLOTS / (double)(1ULL << (b >= 40 ? 40 : b));
    HBITS = 13;
    while (HBITS < 21 && (double)(1u << HBITS) < 1.3 * recmax) HBITS++;
    HSIZE = 1u << HBITS;
    HT = calloc(HSIZE, sizeof(rec));
    if (!HT) { fprintf(stderr, "oom\n"); return 4; }

    uint64_t rng = seed0 ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(b + 1));
    long acc_dc0 = 0, acc_fail = 0, acc_rec = 0, acc_mer = 0, acc_self = 0, acc_ver = 0, acc_suc = 0, acc_esc = 0, acc_dbl = 0;

    for (long tr = 0; tr < trials; tr++) {
        uint64_t x = splitmix64(&rng) % (FQ - 1) + 1;
        pt Qp = smul(G, x);
        if (Qp.inf) continue;
        slot S[SLOTS];
        for (int i = 0; i < SLOTS; i++) {
            uint64_t ci = (uint64_t)(i + 1), ai; pt St;
            for (;;) {
                ai = splitmix64(&rng) % (FQ - 1) + 1;
                uint64_t k0 = radd(ai, rmul(ci, x));
                St = smul(G, k0);
                if (!St.inf) break;
            }
            if (St.y & 1) { St.y = FP - St.y; ai = rsub(0, ai); ci = rsub(0, ci); }
            S[i].x = St.x; S[i].y = St.y; S[i].a = ai; S[i].c = ci;
        }
        memset(HT, 0, (size_t)HSIZE * sizeof(rec));
        qn = 0;
        int got = 0;
        long nrec = 0, nmer = 0, nself = 0, nver = 0, nesc = 0, ndbl = 0, ndc0 = 0, nfail = 0;

        for (long st = 0; st < per; st++) {
            uint64_t num[SLOTS], den[SLOTS], pre[SLOTS + 1];
            int dbl[SLOTS], esc[SLOTS];
            if (getenv("CENSUS_SCALAR")) {                 /* reference path: one inversion per slot */
                for (int i = 0; i < SLOTS; i++) {
                    uint64_t j = widx(S[i].x, (int)(st & 1));
                    pt Tj = *wrow(j, (int)(st & 1));
                    uint64_t na = S[i].a, nc = S[i].c;
                    int d = 0;
                    if (S[i].x == Tj.x) { d = 1; ndbl++; if (S[i].y != Tj.y) nesc++; }
                    pt Sp = padd((pt){S[i].x, S[i].y, 0}, Tj);
                    if (Sp.inf) { ndbl++; nesc++; Sp = padd((pt){S[i].x, S[i].y, 0}, (pt){S[i].x, S[i].y, 0}); }
                    if (d) { na = radd(na, na); nc = radd(nc, nc); }
                    else   { na = radd(na, (uint64_t)((st & 1) ? TLOG2[j] : TLOG[j])); }
                    if (Sp.y & 1) { Sp.y = FP - Sp.y; na = rsub(0, na); nc = rsub(0, nc); }
                    S[i].x = Sp.x; S[i].y = Sp.y; S[i].a = na; S[i].c = nc;
                    if (getenv("CENSUS_INVARIANT")) {
                        pt v = smul(G, radd(na, rmul(nc, x)));
                        if (v.inf || v.x != Sp.x || v.y != Sp.y) { fprintf(stderr, "INVARIANT FAIL tr%ld slot%d\n", tr, i); return 9; }
                    }
                    if (getenv("CENSUS_DUMP") && st < 6) fprintf(stderr, "D %ld %d %llu %llu\n", st, i, (unsigned long long)Sp.x, (unsigned long long)Sp.y);
                    if (((uint64_t)mix64(Sp.x ^ SALT) & mask) == 0) {
                        if (qn < 4096) { Q_[qn].key = Sp.x; Q_[qn].a = na; Q_[qn].c = nc; Q_[qn].owner = i; qn++; }
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
                                        else { if (nfail < 2) fprintf(stderr, "FAILDBG tr%ld\n", tr); nfail++; }
                                    } else ndc0++;
                                } else nself++;
                                break;
                            }
                        }
                    }
                    qn = 0;
                }
                continue;
            }
            for (int i = 0; i < SLOTS; i++) {
                uint64_t j = widx(S[i].x, (int)(st & 1));
                const pt *row = wrow(j, (int)(st & 1));
                pt Tj = *row;
                dbl[i] = 0; esc[i] = 0;
                if (S[i].x == Tj.x) {
                    dbl[i] = 1; ndbl++;
                    if (S[i].y != Tj.y) { esc[i] = 1; nesc++; }   /* zero denominator */
                    num[i] = madd(mm(3, mm(S[i].x, S[i].x)), 0);  /* a = 0 */
                    den[i] = madd(S[i].y, S[i].y);
                } else {
                    num[i] = msub(Tj.y, S[i].y);
                    den[i] = msub(Tj.x, S[i].x);
                }
                if (den[i] == 0) den[i] = 1;                      /* never on this fixture */
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
                uint64_t j = widx(S[i].x, (int)(st & 1));
                pt Tj = *wrow(j, (int)(st & 1));
                uint64_t x1 = S[i].x, y1 = S[i].y, na = S[i].a, nc = S[i].c;
                uint64_t x3, y3;
                if (esc[i]) {                    /* A = -Tj: old code doubles; keep the counter */
                    pt Sp = padd((pt){x1, y1, 0}, (pt){x1, y1, 0});
                    x3 = Sp.x; y3 = Sp.y; na = radd(na, na); nc = radd(nc, nc);
                } else {
                    uint64_t l = lam_[i];
                    x3 = msub(msub(mm(l, l), x1), Tj.x);
                    y3 = msub(mm(l, msub(x1, x3)), y1);
                    if (dbl[i]) { na = radd(na, na); nc = radd(nc, nc); }
                    else        { na = radd(na, (uint64_t)((st & 1) ? TLOG2[j] : TLOG[j])); }
                }
                if (y3 & 1) { y3 = FP - y3; na = rsub(0, na); nc = rsub(0, nc); }
                S[i].x = x3; S[i].y = y3; S[i].a = na; S[i].c = nc;
                if (getenv("CENSUS_DUMP") && st < 6) fprintf(stderr, "D %ld %d %llu %llu\n", st, i, (unsigned long long)x3, (unsigned long long)y3);
                if (getenv("CENSUS_INVARIANT")) {
                    pt v = smul(G, radd(na, rmul(nc, x)));
                    if (v.inf || v.x != x3 || v.y != y3) { fprintf(stderr, "INVARIANT FAIL tr%ld slot%d\n", tr, i); return 9; }
                }
                if (((uint64_t)mix64(x3 ^ SALT) & mask) == 0) {
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
                        if (HT[idx].key == 0) {
                            HT[idx].key = k; HT[idx].a = r.a; HT[idx].c = r.c; HT[idx].owner = r.owner;
                            break;
                        }
                        if (HT[idx].key == k) {
                            if (HT[idx].owner != r.owner) {
                                nmer++;
                                uint64_t dc = rsub(HT[idx].c, r.c);
                                if (dc) {
                                    uint64_t xs = rmul(rsub(r.a, HT[idx].a), rinv(dc));
                                    pt chk = smul(G, xs);
                                    if (!chk.inf && chk.x == Qp.x && chk.y == Qp.y) { nver++; got = 1; }
                                    else { if (nfail < 2) fprintf(stderr, "FAILDBG tr%ld\n", tr); nfail++; }
                                } else ndc0++;
                            } else nself++;
                            break;
                        }
                    }
                }
                qn = 0;
            }
        }
        acc_dc0 += ndc0; acc_fail += nfail; acc_rec += nrec; acc_mer += nmer; acc_self += nself;
        acc_ver += nver; acc_suc += got; acc_esc += nesc; acc_dbl += ndbl;
    }

    printf("{\"dp_bits\":%d,\"seed\":%llu,\"trials\":%ld,\"slots\":%d,\"per_walker\":%ld,"
           "\"budget_mult\":%ld,\"steps\":%ld,\"records\":%ld,\"merges\":%ld,\"self_hits\":%ld,"
           "\"escaping\":%ld,\"doubling_path\":%ld,\"verified_recoveries\":%ld,\"successes\":%ld,"
           "\"dc0\":%ld,\"check_fail\":%ld,"
           "\"fixture\":{\"p\":%llu,\"a\":%d,\"b\":%llu,\"order\":%llu,\"classes\":%llu,"
           "\"table_rows\":%d,\"tables\":%d,\"index_mode\":\"%s\",\"lambda_num\":%llu,\"lambda_den\":%llu}}\n",
           b, (unsigned long long)seed0, trials, SLOTS, per, mult, (long)Ttot * trials, acc_rec,
           acc_mer, acc_self, acc_esc, acc_dbl, acc_ver, acc_suc, acc_dc0, acc_fail,
           (unsigned long long)FP, 0, (unsigned long long)CB, (unsigned long long)FQ,
           (unsigned long long)NCLASS, NB, IMODE ? 2 : 1, IMODE ? "o3-two-table-salted-mix" : "low11",
           (unsigned long long)lam_n, (unsigned long long)lam_d);
    return 0;
}
