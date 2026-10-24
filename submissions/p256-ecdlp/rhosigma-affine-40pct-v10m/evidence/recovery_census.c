/* recovery_census.c — measure DP-record COALESCENCE on the shipped walk map, guarded and
 * unguarded, at identical fixture parameters.
 *
 * WHY (judge d739ebd2 / cbc03d46 item 3): the nine recovery cells in invariant_fit.json were
 * walked on the PRE-GUARD map.  The shipped map selects the row from (class XOR previous
 * row), so two walkers that reach the SAME point with DIFFERENT guard registers generally
 * take different next jumps and need not follow a common tail.  Whether that costs
 * coalescence is a measurable question, and this program measures it: same curve, same
 * classes, same DP mask, same walker count, same length -- only the row-selection rule
 * differs between the two arms.
 *
 * METRIC: a "meet" is a DP record whose key x is already held by a DIFFERENT walker.  x
 * determines the point (y is canonicalised even), so a meet is a genuine full-point
 * collision between distinct walkers -- exactly the recovery event the success law counts.
 * meets/records is the per-record coalescence rate; the ratio guarded/unguarded is the
 * factor the guarded transition costs, and it multiplies the fitted p directly.
 *
 * build: clang -O2 -o recovery_census recovery_census.c ; ./recovery_census [walkers] [steps] [dpbits]
 */
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

typedef uint64_t u64;
static u64 P = 33554393ULL, CA, CB;

static u64 sub_(u64 a, u64 b) { return a >= b ? a - b : a + P - b; }
static u64 add_(u64 a, u64 b) { u64 s = a + b; return s >= P ? s - P : s; }
static u64 mul_(u64 a, u64 b) { return (a * b) % P; }
static u64 inv_(u64 a)
{
    int64_t t = 0, nt = 1, r = (int64_t)P, nr = (int64_t)a;
    while (nr) { int64_t q = r / nr, tmp; tmp = t - q * nt; t = nt; nt = tmp; tmp = r - q * nr; r = nr; nr = tmp; }
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
        if (add_(a.y, b.y) == 0) { r.inf = 1; r.x = r.y = 0; return r; }
        lam = mul_(add_(mul_(3, mul_(a.x, a.x)), CA), inv_(add_(a.y, a.y)));
    } else lam = mul_(sub_(b.y, a.y), inv_(sub_(b.x, a.x)));
    r.x = sub_(sub_(mul_(lam, lam), a.x), b.x);
    r.y = sub_(mul_(lam, sub_(a.x, r.x)), a.y);
    return r;
}
static pt pmul(pt g, u64 k)
{
    pt acc; acc.inf = 1; acc.x = acc.y = 0;
    for (int i = 63; i >= 0; i--) { acc = padd(acc, acc); if ((k >> i) & 1) acc = padd(acc, g); }
    return acc;
}
static u64 mix64(u64 v)
{ v *= 0x9e3779b97f4a7c15ULL; v ^= v >> 29; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 32; return v; }
#define SALT0 0xa5a5a5a5a5a5a5a5ULL
#define SALT1 0x5a5a5a5a5a5a5a5aULL
#define ROWS 2048
static pt T0[ROWS], T1[ROWS];
static u64 T0o[ROWS], T1o[ROWS];
static u64 powm_(u64 b, u64 e) { u64 a = 1; while (e) { if (e & 1) a = mul_(a, b); b = mul_(b, b); e >>= 1; } return a; }
static int modsqrt_(u64 n, u64 *out)
{
    if (n == 0) { *out = 0; return 1; }
    if (powm_(n, (P - 1) / 2) != 1) return 0;
    if (P % 4 == 3) { *out = powm_(n, (P + 1) / 4); return mul_(*out, *out) == n; }
    u64 q = P - 1, S = 0;
    while ((q & 1) == 0) { q >>= 1; S++; }
    u64 z = 2; while (powm_(z, (P - 1) / 2) != P - 1) z++;
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
static pt G;
static void build(void)
{
    CA = P - 3; CB = 12345 % P;
    for (u64 x = 1; x < 100000; x++) {
        u64 rhs = add_(add_(mul_(mul_(x, x), x), mul_(CA, x)), CB), y;
        if (modsqrt_(rhs, &y)) { G.x = x; G.y = y; G.inf = 0; break; }
    }
    u64 s1 = 0xabcdef0123456789ULL, s2 = 0x13579bdf2468ace0ULL;
    for (int j = 0; j < ROWS; j++) {
        u64 t = (mix64(s1 += 0x9e3779b97f4a7c15ULL) | 1) % P;
        pt A = pmul(G, t); int sg = 0;
        if (!A.inf && (A.y & 1)) { A.y = P - A.y; sg = 1; }
        T0[j] = A; T0o[j] = (u64)sg;
        t = (mix64(s2 += 0x9e3779b97f4a7c15ULL) | 1) % P;
        pt B = pmul(G, t); sg = 0;
        if (!B.inf && (B.y & 1)) { B.y = P - B.y; sg = 1; }
        T1[j] = B; T1o[j] = (u64)sg;
    }
}

/* DP record table: key x -> owner */
#define HS (1u << 18)
static u64 *hk, *ho, *hst, *hfd; static u64 gen;
static void hclear(void) { gen++; }
/* returns 1 if a meet with a different owner; inserts otherwise */
static int dp_insert(u64 key, u64 owner, u64 *records, u64 *selfhits)
{
    u64 i = (key * 0x9e3779b97f4a7c15ULL) >> (64 - 18);
    for (;;) {
        if (hst[i] != gen) { hst[i] = gen; hk[i] = key; ho[i] = owner; (*records)++; return 0; }
        if (hk[i] == key) {
            if (ho[i] == owner) { (*selfhits)++; return 0; }
            if (!hfd[i]) { hfd[i] = 1; return 1; }   /* count each collision ONCE */
            return 0;
        }
        i = (i + 1) & (HS - 1);
    }
}

int main(int argc, char **argv)
{
    u64 NW = argc > 1 ? strtoull(argv[1], 0, 10) : 16;
    u64 L  = argc > 2 ? strtoull(argv[2], 0, 10) : 200000;
    int dpb = argc > 3 ? atoi(argv[3]) : 10;
    build();
    hk = calloc(HS, 8); ho = calloc(HS, 8); hst = calloc(HS, 8); hfd = calloc(HS, 8);
    printf("{\"field_prime\":%llu,\"rows\":%d,\"walkers\":%llu,\"steps\":%llu,\"dp_bits\":%d,\"base\":[%llu,%llu],\"arms\":{",
           (unsigned long long)P, ROWS, (unsigned long long)NW, (unsigned long long)L, dpb,
           (unsigned long long)G.x, (unsigned long long)G.y);
    for (int guarded = 0; guarded < 2; guarded++) {
        u64 *prev = calloc(NW, 8), *sign = calloc(NW, 8);
        pt *wl = calloc(NW, sizeof(pt));
        u64 rng = 0xfeedfacecafebeefULL, records = 0, meets = 0, selfhits = 0, dprecs = 0;
        for (u64 w = 0; w < NW; w++) {
            u64 k = (mix64(rng += 0x9e3779b97f4a7c15ULL) % P) | 1;
            wl[w] = pmul(G, k); prev[w] = 0; sign[w] = 0;
        }
        hclear();
        for (u64 step = 0; step < L; step++) {
            for (u64 w = 0; w < NW; w++) {
                u64 par = step & 1;
                u64 c = mix64(wl[w].x ^ (par ? SALT1 : SALT0));
                int j = (int)((guarded ? (c ^ prev[w]) : c) & (ROWS - 1));
                if (guarded) prev[w] = (u64)j;
                pt t = par ? T1[j] : T0[j];
                u64 sg = par ? T1o[j] : T0o[j];
                wl[w] = padd(wl[w], t);
                sign[w] ^= sg;
                if (!wl[w].inf && (wl[w].y & 1)) { wl[w].y = P - wl[w].y; sign[w] ^= 1; }
                if ((mix64(wl[w].x ^ 0xd1b54a32d192ed03ULL) & ((1ULL << dpb) - 1)) == 0) {
                    dprecs++;
                    meets += (u64)dp_insert(wl[w].x, w, &records, &selfhits);
                }
            }
        }
        printf("%s\"%s\":{\"dp_records\":%llu,\"meets\":%llu,\"self_hits\":%llu,\"walker_steps\":%llu,"
               "\"coalescence_per_record\":%.6e}",
               guarded ? "," : "", guarded ? "guarded" : "unguarded",
               (unsigned long long)records, (unsigned long long)meets, (unsigned long long)selfhits,
               (unsigned long long)(NW * L), records ? (double)meets / (double)records : 0.0);
        free(prev); free(sign); free(wl);
    }
    printf("}}\n");
    return 0;
}
