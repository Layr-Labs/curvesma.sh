/* walkmap_census.c — direct census of the p256-ecdlp sigma walk map, guarded and
 * unguarded, at fixture scale.
 *
 * WHY THIS EXISTS (judge c0376361 item 2): the v10k map selected the jump row from the
 * current class alone, so an A,B,A,B class orbit with canonicalisation signs +,-,+,-
 * returns H to itself exactly: H -> H+A -> -(H+A+B) -> -H-B -> H.  No v10k evidence
 * censused that family, and the 50,764-step census could not have seen it anyway.
 *
 * THIS PROGRAM censuses the map, not the curve: the map is
 *     par  = len & 1
 *     c    = mix64(x.low ^ SALT[par])          (the class)
 *     j    = c & 2047                          UNGUARDED (v10k)
 *     j    = (c ^ prev) & 2047, prev = j_prev   GUARDED   (v10l)
 *     P    = canon(P + T_par[j])               canon: y <- -y if y odd; log sign flips
 * and it is scale free: it consumes only the low limb of x and the parity of y.  The
 * curve, its field and its group order therefore do not enter the trap structure; the
 * fixture is a real short-Weierstrass curve over a 2^25 prime with the same 2048-row,
 * two-table, canonicalising map and a 2^25-point group.
 *
 * OUTPUT: for each map, per-start exact cycle detection (a repeated walk STATE is a
 * permanent trap: the map is deterministic and stateless apart from the guard register),
 * period histogram, basin mass, and entry latency.  A period-4 A,B,A,B trap is the family
 * the judge named; every period <= 6 is reported.
 *
 * build: clang -O2 -o walkmap_census walkmap_census.c
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef uint64_t u64;
static u64 P, CA, CB;

/* ---------- modular arithmetic (P < 2^26, products fit in u64) ---------- */
static u64 sub_(u64 a, u64 b) { return a >= b ? a - b : a + P - b; }
static u64 add_(u64 a, u64 b) { u64 s = a + b; return s >= P ? s - P : s; }
static u64 mul_(u64 a, u64 b) { return (a * b) % P; }
static u64 inv_(u64 a)
{
    int64_t t = 0, newt = 1;
    int64_t r = (int64_t)P, newr = (int64_t)a;
    while (newr) {
        int64_t q = r / newr, tmp;
        tmp = t - q * newt; t = newt; newt = tmp;
        tmp = r - q * newr; r = newr; newr = tmp;
    }
    if (t < 0) t += (int64_t)P;
    return (u64)t;
}

typedef struct { u64 x, y; int inf; } pt;

static pt padd(pt a, pt b)
{
    pt r; r.inf = 0;
    if (a.inf) return b;
    if (b.inf) return a;
    u64 lam;
    if (a.x == b.x) {
        if (add_(a.y, b.y) == 0) { r.inf = 1; r.x = r.y = 0; return r; }   /* P + (-P) */
        lam = mul_(add_(mul_(3, mul_(a.x, a.x)), CA), inv_(add_(a.y, a.y)));
    } else {
        lam = mul_(sub_(b.y, a.y), inv_(sub_(b.x, a.x)));
    }
    r.x = sub_(sub_(mul_(lam, lam), a.x), b.x);
    r.y = sub_(mul_(lam, sub_(a.x, r.x)), a.y);
    return r;
}
static pt pmul(pt g, u64 k)
{
    pt acc; acc.inf = 1; acc.x = acc.y = 0;
    for (int i = 63; i >= 0; i--) {
        acc = padd(acc, acc);
        if ((k >> i) & 1) acc = padd(acc, g);
    }
    return acc;
}

/* ---------- the map ---------- */
static u64 mix64(u64 v)
{ v *= 0x9e3779b97f4a7c15ULL; v ^= v >> 29; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 32; return v; }
#define SALT0 0xa5a5a5a5a5a5a5a5ULL
#define SALT1 0x5a5a5a5a5a5a5a5aULL
#define ROWS 2048
static pt T0[ROWS], T1[ROWS], T0n[ROWS], T1n[ROWS];   /* canonical copies */
static u64 T0odd[ROWS], T1odd[ROWS];                 /* stored sign flip */

static pt canon(pt p, int *sign)
{
    if (!p.inf && (p.y & 1)) { p.y = P - p.y; *sign ^= 1; }
    return p;
}
static pt step_unguarded(pt p, u64 len, int *sign)
{
    u64 par = len & 1;
    u64 c = mix64(p.x ^ (par ? SALT1 : SALT0));
    int j = (int)(c & (ROWS - 1));
    pt t = par ? T1[j] : T0[j];
    int s2 = par ? (int)T1odd[j] : (int)T0odd[j];
    pt r = padd(p, t);
    *sign ^= s2;
    return canon(r, sign);
}
static pt step_guarded(pt p, u64 len, u64 *prev, int *sign)
{
    u64 par = len & 1;
    u64 c = mix64(p.x ^ (par ? SALT1 : SALT0));
    int j = (int)((c ^ *prev) & (ROWS - 1));
    *prev = (u64)j;
    pt t = par ? T1[j] : T0[j];
    int s2 = par ? (int)T1odd[j] : (int)T0odd[j];
    pt r = padd(p, t);
    *sign ^= s2;
    return canon(r, sign);
}

/* ---------- visited-state hash set (generation stamped) ---------- */
#define HS 22
static u64 *hk, *hstamp;
static u64 gen;
static u64 *hval;
static void hclear(void) { gen++; }
static int hprobe(u64 key, u64 at)
{
    u64 i = (key * 0x9e3779b97f4a7c15ULL) >> (64 - HS);
    for (;;) {
        if (hstamp[i] != gen) { hstamp[i] = gen; hk[i] = key; hval[i] = at; return -1; }
        if (hk[i] == key) return (int)1;
        i = (i + 1) & ((1u << HS) - 1);
    }
}
static u64 hget(u64 key)
{
    u64 i = (key * 0x9e3779b97f4a7c15ULL) >> (64 - HS);
    for (;;) {
        if (hstamp[i] != gen) return 0;
        if (hk[i] == key) return hval[i];
        i = (i + 1) & ((1u << HS) - 1);
    }
}

/* ---------- modular square root (Tonelli-Shanks, works for any odd prime) ---------- */
static u64 powm_(u64 b, u64 e) { u64 a = 1; while (e) { if (e & 1) a = mul_(a, b); b = mul_(b, b); e >>= 1; } return a; }
static int modsqrt_(u64 n, u64 *out)
{
    if (n == 0) { *out = 0; return 1; }
    if (powm_(n, (P - 1) / 2) != 1) return 0;
    if (P % 4 == 3) { *out = powm_(n, (P + 1) / 4); return mul_(*out, *out) == n; }
    u64 q = P - 1, S = 0;
    while ((q & 1) == 0) { q >>= 1; S++; }
    u64 z = 2;
    while (powm_(z, (P - 1) / 2) != P - 1) z++;
    u64 M = S, c = powm_(z, q), t = powm_(n, q), R = powm_(n, (q + 1) / 2);
    while (t != 1) {
        u64 i = 0, t2 = t;
        while (t2 != 1 && i < M) { t2 = mul_(t2, t2); i++; }
        if (i == M) return 0;
        u64 b = c, j;
        for (j = 0; j + i < M - 1; j++) b = mul_(b, b);
        M = i; c = mul_(b, b); t = mul_(t, c); R = mul_(R, b);
    }
    *out = R; return mul_(R, R) == n;
}
static pt FIND_G;

/* ---------- curve construction: pick G, fill the two jump tables ---------- */
static void build_curve(void)
{
    P = 33554393ULL;                 /* 2^25 - 39, prime */
    CA = P - 3; CB = 12345 % P;
    int found = 0;
    for (u64 x = 1; x < 100000 && !found; x++) {
        u64 rhs = add_(add_(mul_(mul_(x, x), x), mul_(CA, x)), CB);
        u64 y;
        if (modsqrt_(rhs, &y)) { FIND_G.x = x; FIND_G.y = y; FIND_G.inf = 0; found = 1; }
    }
    if (!found) { fprintf(stderr, "no base point\n"); exit(2); }
    u64 s1 = 0xabcdef0123456789ULL, s2 = 0x13579bdf2468ace0ULL;
    for (int j = 0; j < ROWS; j++) {
        u64 t = (mix64(s1 += 0x9e3779b97f4a7c15ULL) | 1) % P;
        pt A = pmul(FIND_G, t);
        int sg = 0; A = canon(A, &sg); T0[j] = A; T0odd[j] = (u64)sg;
        t = (mix64(s2 += 0x9e3779b97f4a7c15ULL) | 1) % P;
        pt B = pmul(FIND_G, t);
        sg = 0; B = canon(B, &sg); T1[j] = B; T1odd[j] = (u64)sg;
    }
}

int main(int argc, char **argv)
{
    u64 starts = argc > 1 ? strtoull(argv[1], 0, 10) : 2000;
    u64 L = argc > 2 ? strtoull(argv[2], 0, 10) : 50000;
    build_curve();
    hk = calloc(1u << HS, 8); hstamp = calloc(1u << HS, 8); hval = calloc(1u << HS, 8);
    pt G = FIND_G;
    printf("{\"field_prime\":%llu,\"rows\":%d,\"starts\":%llu,\"max_steps\":%llu,"
           "\"base\":[%llu,%llu],\"maps\":{", (unsigned long long)P, ROWS,
           (unsigned long long)starts, (unsigned long long)L,
           (unsigned long long)G.x, (unsigned long long)G.y);

    for (int guarded = 0; guarded < 2; guarded++) {
        u64 per[13]; memset(per, 0, sizeof per);
        u64 ncyc = 0, entrysum = 0, entered_short = 0;
        u64 return4 = 0, classlag2 = 0, return2 = 0, stepsrun = 0;   /* direct family counters */
        u64 period4 = 0, period2 = 0, period6 = 0;
        u64 rng = 0x1234567890abcdefULL;
        for (u64 s = 0; s < starts; s++) {
            u64 k = (mix64(rng += 0x9e3779b97f4a7c15ULL) % P) | 1;
            pt p = pmul(G, k);
            u64 prev = 0; int sign = 0;
            hclear();
            u64 i;
            int cyc = 0, per_ = 0, ent = 0;
            u64 rx[4]; int nrx = 0; u64 cls[2]; int ncls = 0;
            for (i = 0; i < L; i++) {
                u64 key = (p.x << (26 + 11)) | (p.y << 11) | (guarded ? (prev & 0x7f) : 0) | ((u64)(i & 1) << 62);
                /* DIRECT family census: x == x_{i-4} is the A,B,A,B period-4 closure
                   (x determines the point: y is canonicalised even), and a class
                   lag-2 repeat is its necessary precursor.  Counted independently of
                   the first-repeat detector, which a finite group dominates. */
                stepsrun++;
                if (nrx >= 4 && rx[(nrx) & 3] == p.x) return4++;
                if (nrx >= 2 && rx[(nrx + 2) & 3] == p.x) return2++;
                { u64 c = mix64(p.x ^ (i & 1 ? SALT1 : SALT0)) & 2047;
                  if (ncls >= 2 && cls[ncls & 1] == c) classlag2++;
                  cls[ncls & 1] = c; ncls++; }
                rx[nrx & 3] = p.x; nrx++;
                if (hprobe(key, i) > 0) { cyc = 1; ent = (int)hget(key); per_ = (int)(i - ent); break; }
                if (guarded) p = step_guarded(p, i, &prev, &sign);
                else p = step_unguarded(p, i, &sign);
            }
            if (cyc) {
                ncyc++; entrysum += ent;
                if (per_ <= 12) per[per_]++;
                if (per_ == 2) period2++;
                if (per_ == 4) period4++;
                if (per_ == 6) period6++;
                if (per_ <= 6) entered_short++;
            }
        }
        printf("%s\"%s\":{\"cycles\":%llu,\"basin_mass\":%.6f,\"mean_entry_step\":%.1f,"
               "\"periods\":{", guarded ? "," : "", guarded ? "guarded" : "unguarded",
               (unsigned long long)ncyc, (double)ncyc / (double)starts,
               ncyc ? (double)entrysum / (double)ncyc : 0.0);
        int first = 1;
        for (int q = 1; q <= 12; q++) if (per[q]) {
            printf("%s\"%d\":%llu", first ? "" : ",", q, (unsigned long long)per[q]); first = 0; }
        printf("},\"period2\":%llu,\"period4\":%llu,\"period6\":%llu,\"cycles_le6\":%llu,"
               "\"steps_run\":%llu,\"return4\":%llu,\"return4_rate\":%.3e,"
               "\"return2\":%llu,\"class_lag2\":%llu,\"class_lag2_rate\":%.6f}",
               (unsigned long long)period2, (unsigned long long)period4,
               (unsigned long long)period6, (unsigned long long)entered_short,
               (unsigned long long)stepsrun, (unsigned long long)return4,
               stepsrun ? (double)return4 / (double)stepsrun : 0.0,
               (unsigned long long)return2, (unsigned long long)classlag2,
               stepsrun ? (double)classlag2 / (double)stepsrun : 0.0);
    }
    printf("}}\n");
    return 0;
}
