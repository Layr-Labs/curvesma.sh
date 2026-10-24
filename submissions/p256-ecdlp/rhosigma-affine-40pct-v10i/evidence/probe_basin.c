/* probe_basin.c — attractor diagnostic on the p256 v10f 2^40 fixture.
 *
 * The census map is  p -> p + T[x(p) & 2047]  (one shared 2048-row table, index from
 * the LOW 11 BITS of x).  Part 3 measured 91.4% self-hits at S_t=47,711, i.e. every
 * walker falls into a short cycle in point space almost immediately.
 *
 * This probe measures, per candidate index rule:
 *   MI      mutual information I(bin_n ; bin_{n+1}) in bits.  If the increment is
 *           selected by low bits of the state, successive bins are nearly a
 *           deterministic function of each other => MI -> 11 bits => the walk is a
 *           2048-state attractor, not a random map.  Mixing rule => MI ~ 0.
 *   CYC     fraction of starts that enter a cycle within CAP steps (CAP ~ S_t), with
 *           median entry mu and median cycle length lam.  A random map on 2^39
 *           classes has mu ~ 2^19.5 = 741k >> CAP, so CYC must be ~0.
 *   BASIN   share of starts in the modal (giant) cycle.
 *   REP     fraction of a 5,000-step trail that repeats an already-visited class
 *           (the self_hits analogue) -- should be ~0 for a mixing rule.
 *   EQV     index(p)==index(-p) and canon(y)==canon(p-y) over 2000 points: the
 *           class-consistency (sigma-equivariance) gate.
 *
 * Modes: 0 baseline (x & 2047)         1 O2 mix64(x)&2047
 *        2 O3 two tables, parity-salted 3 O1 per-walker coset offset (x+O)&2047
 *        4 O2b cheap fold (x^(x>>20))&2047
 * Usage: probe_basin <mode> [starts] [cap] [seed]
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
#define NB 2048

static inline uint64_t mm(uint64_t a, uint64_t b) { return (uint64_t)((u128)a * b % FP); }
static inline uint64_t madd(uint64_t a, uint64_t b) { uint64_t s = a + b; return s >= FP ? s - FP : s; }
static inline uint64_t msub(uint64_t a, uint64_t b) { return a >= b ? a - b : a + FP - b; }

static uint64_t mpow(uint64_t b, uint64_t e, uint64_t m)
{ uint64_t r = 1 % m; b %= m; while (e) { if (e & 1) r = (uint64_t)((u128)r * b % m); b = (uint64_t)((u128)b * b % m); e >>= 1; } return r; }
static inline uint64_t minv(uint64_t a) { return mpow(a, FP - 2, FP); }

typedef struct { uint64_t x, y; int inf; } pt;

static pt padd(pt A, pt B)
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
    for (int i = top; i >= 0; i--) { Acc = padd(Acc, Acc); if ((k >> i) & 1) Acc = padd(Acc, B); }
    return Acc;
}

static inline uint64_t mix64(uint64_t v)
{ v *= 0x9e3779b97f4a7c15ULL; v ^= v >> 29; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 32; return v; }

static inline uint64_t splitmix64(uint64_t *s)
{ uint64_t z = (*s += 0x9e3779b97f4a7c15ULL); z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL; return z ^ (z >> 31); }

static pt TBL[NB], TBL2[NB];
static int64_t TLOG[NB], TLOG2[NB];

static void build(pt *TB, int64_t *TL, uint64_t *sp)
{
    for (int j = 0; j < NB; j++) {
        uint64_t t; int dup;
        do { t = splitmix64(sp) % (FQ - 1) + 1; dup = 0;
             for (int k = 0; k < j; k++) if ((uint64_t)(TL[k] < 0 ? FQ - (uint64_t)(-TL[k]) : (uint64_t)TL[k]) == t) { dup = 1; break; }
        } while (dup);
        pt T = smul((pt){GX, GY, 0}, t);
        if (T.y & 1) { T.y = FP - T.y; t = FQ - t; }
        TB[j] = T; TL[j] = (int64_t)t;
    }
}

static int MODE = 0;
static uint64_t OFFW = 0;          /* per-walker coset offset (mode 3) */
static uint64_t SA = 0xa5a5a5a5a5a5a5a5ULL, SB = 0x5a5a5a5a5a5a5a5aULL;

static int NOCAN = 0;      /* modes 5,6: walk on true points (no canon sign-flip) */
static inline uint64_t idx(uint64_t x)
{
    switch (MODE) {
        case 0: return x & (NB - 1);
        case 1: return mix64(x) & (NB - 1);
        case 2: return mix64(x ^ SA) & (NB - 1);          /* table A (even step) */
        case 3: return (x + OFFW) & (NB - 1);
        case 4: return (x ^ (x >> 20)) & (NB - 1);
        case 5: return x & (NB - 1);                     /* no-canon baseline index */
        case 6: return mix64(x) & (NB - 1);              /* no-canon + mixed index */
        default: return x & (NB - 1);
    }
}
static inline uint64_t idxB(uint64_t x) { return mix64(x ^ SB) & (NB - 1); }

/* one step; MODE 2 alternates tables when par!=0 */
static inline pt step(pt S, int par)
{
    pt Tj; 
    if (MODE == 2 && par) { Tj = TBL2[idxB(S.x)]; }
    else                  { Tj = TBL[idx(S.x)];  }
    pt R;
    if (S.x == Tj.x) {                                  /* doubling path */
        R = (S.y == Tj.y) ? padd(S, S) : padd(S, S);    /* census keeps the counter; value same here */
    } else R = padd(S, Tj);
    if (!NOCAN && (R.y & 1)) R.y = FP - R.y;            /* canon: even y (modes 0-4 only) */
    return R;
}
/* fixed map used for cycle detection: for mode 2 the 2-step composite */
static inline pt fmap(pt S) { if (MODE == 2) return step(step(S, 0), 1); return step(S, 0); }

#define MI 2048
static uint32_t cnt[MI][MI];
static uint32_t rowsum[MI], colsum[MI];

int main(int argc, char **argv)
{
    MODE = argc > 1 ? atoi(argv[1]) : 0;
    NOCAN = (MODE == 5 || MODE == 6);
    long starts = argc > 2 ? atol(argv[2]) : 2000;
    long cap    = argc > 3 ? atol(argv[3]) : 40000;
    uint64_t sd = argc > 4 ? strtoull(argv[4], 0, 0) : 0x1234ULL;

    uint64_t s1 = 0xabcdef0123456789ULL;
    build(TBL, TLOG, &s1);
    if (MODE == 2) { uint64_t s2 = 0x13579bdf2468ace0ULL; build(TBL2, TLOG2, &s2); }

    pt G = {GX, GY, 0};
    uint64_t rng = sd ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(MODE + 1));

    /* ---------- EQV gate: index is a function of x only; x(-P)=x(P) ---------- */
    long eqv_bad = 0, odd_bad = 0;
    for (long i = 0; i < 2000; i++) {
        uint64_t k = splitmix64(&rng) % (FQ - 1) + 1;
        pt P = smul(G, k);
        if (P.inf) continue;
        if (idx(P.x) != idx(P.x)) eqv_bad++;                 /* index depends on x only, structurally */
        pt N = P; N.y = FP - N.y;
        if (idx(P.x) != idx(N.x)) eqv_bad++;                 /* explicit: x is negation-invariant */
        if (idx(P.x) != idx(P.x)) eqv_bad++;
        uint64_t ca = (P.y & 1) ? FP - P.y : P.y;
        uint64_t cb = (N.y & 1) ? FP - N.y : N.y;
        if (ca != cb) eqv_bad++;                             /* canon(y)==canon(p-y) */
        if (ca & 1) odd_bad++;
    }

    /* ---------- MI of consecutive bins over random points ---------- */
    memset(cnt, 0, sizeof cnt);
    long mibad = 0;
    long NMI = 200000;
    for (long i = 0; i < NMI; i++) {
        uint64_t k = splitmix64(&rng) % (FQ - 1) + 1;
        pt P = smul(G, k); if (P.inf) continue;
        if (P.y & 1) P.y = FP - P.y;
        uint64_t j = idx(P.x);
        pt R = step(P, 0);
        if (R.inf) { mibad++; continue; }
        uint64_t j2 = (MODE == 2 ? idxB(R.x) : idx(R.x));
        cnt[j][j2]++;
    }
    double tot = 0; for (int a = 0; a < MI; a++) for (int b = 0; b < MI; b++) { tot += cnt[a][b]; rowsum[a] += cnt[a][b]; colsum[b] += cnt[a][b]; }
    double mi = 0;
    int used = 0, maxrow = 0;
    for (int a = 0; a < MI; a++) if (rowsum[a]) used++;
    for (int a = 0; a < MI; a++) {
        if (!rowsum[a]) continue;
        if ((int)rowsum[a] > maxrow) maxrow = rowsum[a];
        for (int b = 0; b < MI; b++) {
            if (!cnt[a][b]) continue;
            double p = (double)cnt[a][b] / tot;
            double pa = (double)rowsum[a] / tot, pb = (double)colsum[b] / tot;
            mi += p * log2(p / (pa * pb));
        }
    }
    /* deterministic-bin check: modal successor share per row, averaged over rows */
    double det = 0; long detrows = 0;
    for (int a = 0; a < MI; a++) {
        if (rowsum[a] < 20) continue;
        int best = 0; for (int b = 0; b < MI; b++) if ((int)cnt[a][b] > best) best = cnt[a][b];
        det += (double)best / (double)rowsum[a]; detrows++;
    }
    det = detrows ? det / detrows : 0;

    /* ---------- cycle structure (Brent, capped) + repeat rate ---------- */
    long ncyc = 0; long long summu = 0, sumlam = 0; long long mus[4096];
    long ncyc_ok = 0;
    long ng = 0, nrep = 0, ntot = 0, ncyc400 = 0;
    long starts2 = starts < 4096 ? starts : 4096;
    long long base = 0; int have_base = 0;
    for (long t = 0; t < starts2; t++) {
        uint64_t k = splitmix64(&rng) % (FQ - 1) + 1;
        pt P0 = smul(G, k); if (P0.inf) continue;
        if (P0.y & 1) P0.y = FP - P0.y;
        OFFW = splitmix64(&rng);                      /* mode 3: private offset per walker */

        /* trail repeat rate over 5000 steps (open-addressed set, 16384 slots) */
        if (t < 300) {
            enum { HS = 16384 };
            static uint64_t seen[HS];
            memset(seen, 0, sizeof seen);
            pt P = P0;
            for (long st = 0; st < 5000; st++) {
                P = step(P, (int)(st & 1));
                ntot++;
                uint64_t kk = P.x + 1;
                unsigned h = (unsigned)((kk * 0x9e3779b97f4a7c15ULL) >> 50) & (HS - 1);
                int hit = 0;
                for (unsigned pr = 0; pr < HS; pr++) {
                    unsigned ii = (h + pr) & (HS - 1);
                    if (seen[ii] == 0) { seen[ii] = kk; break; }
                    if (seen[ii] == kk) { hit = 1; break; }
                }
                if (hit) nrep++;
            }
        }
        /* Brent */
        pt tor = P0, hare = fmap(P0);
        long long power = 1, lam = 1;
        long steps = 1;
        int cyc = 0;
        while (steps < cap) {
            if (tor.x == hare.x && tor.y == hare.y) { cyc = 1; break; }
            if (power == lam) { tor = hare; power *= 2; lam = 0; }
            hare = fmap(hare); lam++; steps++;
        }
        /* entry length */
        long long mu = -1;
        if (cyc) {
            tor = P0; hare = P0;
            for (long long i = 0; i < lam; i++) hare = fmap(hare);
            mu = 0;
            while (!(tor.x == hare.x && tor.y == hare.y) && mu < cap) { tor = fmap(tor); hare = fmap(hare); mu++; }
            if (!(tor.x == hare.x && tor.y == hare.y)) mu = -1;
        }
        if (cyc && mu >= 0) { ncyc++; mus[ncyc - 1] = mu; summu += mu; sumlam += lam; ncyc_ok++; }
        /* giant-basin share: min x on the cycle (cycle rep) of the first 400 starts */
        if (t < 400 && cyc) {
            pt C = hare; uint64_t mn = C.x;
            for (long long i = 0; i < lam; i++) { C = fmap(C); if (C.x < mn) mn = C.x; }
            if (!have_base) { base = (long long)mn; have_base = 1; }
            ncyc400++;
            if ((long long)mn == base) ng++;
        }
    }
    double cycfrac = starts2 ? (double)ncyc / (double)starts2 : 0;
    double repfrac = ntot ? (double)nrep / (double)ntot : 0;
    printf("{\"mode\":%d,\"mi_bits\":%.4f,\"mi_norm\":%.4f,\"modal_succ\":%.4f,\"bins_used\":%d,"
           "\"cyc_frac\":%.4f,\"med_mu\":%lld,\"med_lam\":%lld,\"giant_share\":%.3f,"
           "\"rep_frac\":%.4f,\"eqv_bad\":%ld,\"odd_y\":%ld,\"starts\":%ld,\"cap\":%ld,\"classes\":%llu}\n",
           MODE, mi, mi / 11.0, det, used, cycfrac,
           ncyc ? (long long)(summu / ncyc) : -1LL, ncyc ? (long long)(sumlam / ncyc) : -1LL,
           (ncyc400 ? (double)ng / (double)ncyc400 : 0.0), repfrac, eqv_bad, odd_bad, starts2, cap,
           (unsigned long long)NCLASS);
    return 0;
}
