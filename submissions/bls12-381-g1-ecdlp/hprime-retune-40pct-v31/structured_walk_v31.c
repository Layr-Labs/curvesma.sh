// v30 CORRECTED ARM. Fixes verdict-0ad3f80 flaw 1 (wrong group object) and flaw 3 (kappa on the
// wrong statistic), and implements the gap CAP of the cap/abort policy (flaw 4).
//
// ONE object, ONE counter path for every arm: records(side) = { a_i * w_k mod p } over a SHARED
// walk list W of consecutive H-primes (H = index-8 subgroup of order M_H), hit test
//     hit <=> exists u in mu6 : v*u*b_j*w_l == a_i*w_k
// The ONLY knob that differs between arms is draw_base():
//   span_bits = 0  -> a = g^k with k uniform in H's exponent range (full support; the crown arm)
//   span_bits = B  -> a = g^k with k uniform in [1,2^B)   (PRODUCTION POOL LAW: an H-element whose
//                     DISCRETE LOG is small, i.e. exactly "g^k with small discrete log")
// v29 built this arm as pow8(x) = x^8 mod p, i.e. the image of an INTEGER under the 8th-power map.
// That is not g^k for small k: log_g(x^8) is (8 * log of x) and is not small, and for k=8 v29's own
// wrap_free test passed on a value a = 2^24 whose products DO exceed p. Here the arm is built by
// g^k with k < 2^B, so log_g(a_i) = k < 2^B is small BY CONSTRUCTION and span+walk_band < 41 makes
// a_i*w_k an exact (no-mod-p-wrap) product of the exponents as the production 97+77=174-bit key is.
//
// COUNTERS (identical in every arm; this is the "charged arm == reported arm" fix):
//   * residual collisions: r = N - #distinct residues              (the v29 statistic, kept for continuity)
//   * KEY collisions: sum_b C(m_b,2) over the TRUNCATED key the store actually indexes
//     (key = residue & (2^KEY_BITS - 1)); this is the statistic the failure debit is priced on.
// Both are counted in ONE census object and emitted in ONE JSON line -> ledger reads that file.
//
// GAP CAP: the produced walk list is scanned for its gap law (g_k = w_{k+1}-w_k). With gap_cap > 0
// a gap above the cap ABORTS the trail (the cap/abort policy of the note); aborts are counted, so
// the abort rate is a MEASURED number and not a story.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

// v31 CAP POLICY constant, byte-identical to the ledger's spec text and the census JSON it emits:
// "K_max = 2^32 = 4294967296".  The derive script greps THIS literal (three-way identity assert).
#define K_MAX_DEFAULT 4294967296ULL

typedef unsigned __int128 u128;
static uint64_t P, M_H, G, GH, MU6[6], SPAN_BITS, KEY_BITS, GAP_CAP;

static inline uint64_t mulmod(uint64_t a, uint64_t b) { return (uint64_t)((u128)a * b % P); }
static uint64_t powmod(uint64_t b, uint64_t e) {
    uint64_t r = 1; b %= P;
    while (e) { if (e & 1) r = mulmod(r, b); b = mulmod(b, b); e >>= 1; }
    return r;
}
static inline uint64_t mulmod_n(uint64_t a, uint64_t b, uint64_t n) { return (uint64_t)((u128)a * b % n); }
static uint64_t powmod_n(uint64_t b, uint64_t e, uint64_t n) {
    uint64_t r = 1; b %= n;
    while (e) { if (e & 1) r = mulmod_n(r, b, n); b = mulmod_n(b, b, n); e >>= 1; }
    return r;
}
static int is_prime_u64(uint64_t n) {
    if (n < 2) return 0;
    static const uint64_t sp[12] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) { if (n % sp[i] == 0) return n == sp[i]; }
    uint64_t d = n - 1; int s = 0;
    while (!(d & 1)) { d >>= 1; s++; }
    for (int i = 0; i < 12; i++) {
        if (sp[i] >= n) continue;
        uint64_t x = powmod_n(sp[i], d, n);
        if (x == 1 || x == n - 1) continue;
        int ok = 0;
        for (int r = 1; r < s; r++) { x = mulmod_n(x, x, n); if (x == n - 1) { ok = 1; break; } }
        if (!ok) return 0;
    }
    return 1;
}
static inline uint64_t rnd(uint64_t *s) { uint64_t x = *s; x ^= x << 13; x ^= x >> 7; x ^= x << 17; *s = x; return x; }

// PER-TRAIL WINDOW: each trail starts at its own offset in the H-prime list, so its realised gap
// support (and its maximum gap) is a genuine per-trail random variable -- the object the cap and
// the support distribution are about. Production walks start at a fresh point per challenge too.
#define GAP_SPREAD 8192

static uint64_t *W_LIST; static uint64_t W_LEN, W_MIN, W_MAX;
static uint64_t *GAPS; static uint64_t NGAP, GAP_MAXV, GAP_SUM;

static void build_walk(uint64_t len, uint64_t lo) {
    W_LIST = malloc(sizeof(uint64_t) * len);
    GAPS = malloc(sizeof(uint64_t) * len);
    uint64_t q = lo | 1, got = 0, prev = 0;
    while (got < len) {
        uint64_t qq = q % P;
        if (qq > 1 && is_prime_u64(q) && powmod(qq, M_H) == 1) {
            if (got) GAPS[NGAP++] = qq > prev ? qq - prev : prev - qq;
            prev = qq; W_LIST[got++] = qq;
        }
        q += 2;
    }
    W_LIST = realloc(W_LIST, sizeof(uint64_t) * got);
    GAPS = realloc(GAPS, sizeof(uint64_t) * got);
    W_LEN = got;
    W_MIN = W_LIST[0]; W_MAX = W_LIST[0];
    GAP_SUM = 0; GAP_MAXV = 0;
    for (uint64_t i = 1; i < got; i++) { if (W_LIST[i] < W_MIN) W_MIN = W_LIST[i]; if (W_LIST[i] > W_MAX) W_MAX = W_LIST[i]; }
    for (uint64_t i = 0; i < NGAP; i++) { GAP_SUM += GAPS[i]; if (GAPS[i] > GAP_MAXV) GAP_MAXV = GAPS[i]; }
}
static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

// hash set with multiplicity: tracks exact colliding PAIRS = sum_b C(m_b,2) incrementally
typedef struct { uint64_t *key; uint32_t *tag, *cnt; uint32_t epoch; uint64_t mask; uint64_t pairs; } HS;

static void hs_init(HS *h, uint64_t n) {
    uint64_t b = 1; while ((1ULL << b) < 2 * n) b++;
    h->mask = (1ULL << b) - 1;
    h->key = calloc(h->mask + 1, 8); h->tag = calloc(h->mask + 1, 4); h->cnt = calloc(h->mask + 1, 4);
    h->epoch = 0; h->pairs = 0;
}
static void hs_free(HS *h) { free(h->key); free(h->tag); free(h->cnt); }
static inline void hs_clear(HS *h) { h->epoch++; h->pairs = 0; if (h->epoch == 0) { memset(h->tag, 0, (h->mask + 1) * 4); h->epoch = 1; } }
// returns 1 if the key was new; maintains exact pair count
static inline int hs_put(HS *h, uint64_t k) {
    uint64_t i = (k * 0x9E3779B97F4A7C15ULL) & h->mask;
    while (h->tag[i] == h->epoch) {
        if (h->key[i] == k) { h->pairs += h->cnt[i]; h->cnt[i]++; return 0; }
        i = (i + 1) & h->mask;
    }
    h->tag[i] = h->epoch; h->key[i] = k; h->cnt[i] = 1; return 1;
}
static inline int hs_has(HS *h, uint64_t k) {
    uint64_t i = (k * 0x9E3779B97F4A7C15ULL) & h->mask;
    while (h->tag[i] == h->epoch) { if (h->key[i] == k) return 1; i = (i + 1) & h->mask; }
    return 0;
}

typedef struct {
    int wid; uint64_t nb, ell, trials, seed;
    uint64_t hits, distinct_sum, maxdup, keypairs_sum, respairs_sum, max_keypairs, max_mixpairs;
    uint64_t aborts, mixpairs_sum;
    uint64_t maxgap;
    uint64_t *kpair_walk, *rpair_walk, *mg_walk, *mpair_walk;
    double seconds;
} Job;

// PRODUCTION PROPOOL LAW: base = g^k with k < 2^SPAN_BITS (small DISCRETE LOG, group membership
// automatic because every g^k lies in H). v29 returned pow8(x): the 8th power of a small integer,
// whose discrete log is NOT small -- a different object.
static inline uint64_t draw_base(uint64_t *s) {
    if (SPAN_BITS == 0) {                      // full-support arm: any H-element = GH^k, k uniform mod M_H
        uint64_t k = rnd(s) % M_H; if (!k) k = 1;
        return powmod(GH, k);
    }
    uint64_t span = 1ULL << SPAN_BITS;
    uint64_t k = 1 + rnd(s) % (span - 1);      // k IS the discrete log: small by construction
    return powmod(GH, k);
}
static inline uint64_t key_of(uint64_t r) {          // v30 low-bit key (retained for the diagnostic)
    return KEY_BITS ? (r & ((1ULL << KEY_BITS) - 1)) : r;
}
// v31 MIXED-DIGEST key: 2-round xorshift-multiply finaliser (murmur3 fmix64 shape) applied to the
// FULL residue BEFORE truncation, so the truncated key inherits every residue bit. 6 ALU ops
// (3 xor-shift, 2 multiply, 1 xor-shift) -- charged as +3 over the low-bit pack in the ledger.
static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}
static inline uint64_t key_mix(uint64_t r) {
    return KEY_BITS ? (mix64(r) & ((1ULL << KEY_BITS) - 1)) : mix64(r);
}

static void *worker(void *arg) {
    Job *j = (Job *)arg;
    uint64_t N = j->nb * j->ell;
    uint64_t *A = malloc(sizeof(uint64_t) * N), *B = malloc(sizeof(uint64_t) * N);
    uint64_t *GHY = malloc(sizeof(uint64_t) * (j->ell ? j->ell : 1));
    HS h, hd, hk, hm;
    hs_init(&h, 2 * N); hs_init(&hd, 2 * N); hs_init(&hk, 2 * N); hs_init(&hm, 2 * N);
    uint64_t s = j->seed * 0x2545F4914F6CDD1DULL + 12345;
    double t0 = (double)clock() / CLOCKS_PER_SEC;
    for (uint64_t t = 0; t < j->trials; t++) {
        // PER-TRAIL WINDOW + CAP/ABORT POLICY: the trail uses its own offset into the H-prime list;
        // if any realised gap in that window exceeds GAP_CAP the trail ABORTS (charged, counted).
        uint64_t off = rnd(&s) % GAP_SPREAD;
        uint64_t mg = 0;
        for (uint64_t k = 0; k < j->ell; k++) if (GAPS[off + k] > mg) mg = GAPS[off + k];
        j->mg_walk[t] = mg;
        if (GAP_CAP && mg > GAP_CAP) { j->aborts++; j->kpair_walk[t] = 0; j->rpair_walk[t] = 0; continue; }
        if (mg > j->maxgap) j->maxgap = mg;
        j->mg_walk[t] = mg;
        hs_clear(&h); hs_clear(&hd); hs_clear(&hk); hs_clear(&hm);
        uint64_t *W = W_LIST + off;
        // v31 RECORD (unchanged object, v30's): A[i][k] = a_i * W_k mod p, a_i = GH^k_i with k_i <
        // 2^SPAN, W_k the k-th H-prime of the shared list.  v31 RETRACTS v30's "exact integer
        // product / no-wrap" claim instead of redefining the object: the field product reduces mod
        // p BY DEFINITION (that is the group law, exact), and the exponent sum k_i + log(W_k) is
        // never asserted small.  What v31 DOES assert and machine-checks is the exponent-domain
        // statement about the POOL: k_i + h_k < 2^97 + 2^77 < M_H, so no exponent reduction occurs
        // (see the wrap row in derive_v31.py), and distinctness of the records is the measured
        // residue-collision statistic, not a no-wrap assumption.
        for (uint64_t i = 0; i < j->nb; i++) {
            uint64_t a = draw_base(&s), b = draw_base(&s);
            for (uint64_t k = 0; k < j->ell; k++) {
                A[i * j->ell + k] = mulmod(a, W[k]);
                B[i * j->ell + k] = mulmod(b, W[k]);
            }
        }
        uint64_t nd = 0;
        for (uint64_t x = 0; x < N; x++) { if (hs_put(&hd, A[x])) nd++;
            hs_put(&h, A[x]); hs_put(&hk, key_of(A[x])); hs_put(&hm, key_mix(A[x])); }
        uint64_t v = powmod(GH, 1 + rnd(&s) % (M_H - 1));  // target law: v uniform in H (production)
        int hit = 0;
        for (int u = 0; u < 6 && !hit; u++) {
            uint64_t vu = mulmod(v, MU6[u]);
            for (uint64_t x = 0; x < N; x++) if (hs_has(&h, mulmod(vu, B[x]))) { hit = 1; break; }
        }
        j->hits += hit; j->distinct_sum += nd;
        uint64_t col = N - nd;
        if (col > j->maxdup) j->maxdup = col;
        j->keypairs_sum += hk.pairs; j->respairs_sum += hd.pairs;
        if (hk.pairs > j->max_keypairs) j->max_keypairs = hk.pairs;
        if (hm.pairs > j->max_mixpairs) j->max_mixpairs = hm.pairs;
        j->kpair_walk[t] = hk.pairs; j->rpair_walk[t] = hd.pairs; j->mpair_walk[t] = hm.pairs;
        j->mixpairs_sum += hm.pairs;
    }
    j->seconds = (double)clock() / CLOCKS_PER_SEC - t0;
    free(A); free(B); free(GHY); hs_free(&h); hs_free(&hd); hs_free(&hk); hs_free(&hm);
    return NULL;
}

static double pct(uint64_t *v, uint64_t n, double p) {
    if (!n) return 0.0;
    uint64_t idx = (uint64_t)ceil(p * (double)n) - 1;
    if (idx >= n) idx = n - 1;
    return (double)v[idx];
}

int main(int argc, char **argv) {
    uint64_t nb = strtoull(argv[1], 0, 10), ell = strtoull(argv[2], 0, 10);
    uint64_t trials = strtoull(argv[3], 0, 10), workers = strtoull(argv[4], 0, 10);
    uint64_t seed0 = argc > 5 ? strtoull(argv[5], 0, 10) : 1;
    SPAN_BITS = argc > 6 ? strtoull(argv[6], 0, 10) : 0;
    KEY_BITS  = argc > 7 ? strtoull(argv[7], 0, 10) : 0;
    GAP_CAP   = argc > 8 ? strtoull(argv[8], 0, 10) : K_MAX_DEFAULT;
    uint64_t p = (1ULL << 41) + 1;
    while (!(is_prime_u64(p) && (p - 1) % 48 == 0)) p += 2;
    P = p; M_H = (p - 1) / 8;
    for (uint64_t g = 2; g < 100000; g++) {
        if (powmod(g, (p - 1) / 2) != 1 && powmod(g, (p - 1) / 3) != 1) { G = g; break; }
    }
    GH = powmod(G, 8);                        // order-M_H generator of the index-8 subgroup H
    for (int j = 0; j < 6; j++) MU6[j] = powmod(G, (p - 1) / 6 * j);
    build_walk(ell + GAP_SPREAD, 1ULL << 20);

    double N = (double)nb * (double)ell;
    Job *jobs = calloc(workers, sizeof(Job));
    pthread_t *th = calloc(workers, sizeof(pthread_t));
    uint64_t base = trials / workers;
    double w0 = time(NULL);
    for (uint64_t i = 0; i < workers; i++) {
        jobs[i].wid = i; jobs[i].nb = nb; jobs[i].ell = ell;
        jobs[i].trials = base + (i < trials - base * workers ? 1 : 0);
        jobs[i].seed = seed0 * 1000 + i;
        jobs[i].kpair_walk = calloc(jobs[i].trials ? jobs[i].trials : 1, 8);
        jobs[i].rpair_walk = calloc(jobs[i].trials ? jobs[i].trials : 1, 8);
        jobs[i].mg_walk = calloc(jobs[i].trials ? jobs[i].trials : 1, 8);
        jobs[i].mpair_walk = calloc(jobs[i].trials ? jobs[i].trials : 1, 8);
        pthread_create(&th[i], NULL, worker, &jobs[i]);
    }
    uint64_t hits = 0, dsum = 0, maxdup = 0, ksum = 0, rsum = 0, kmax = 0, aborts = 0;
    uint64_t msum = 0, mmax = 0;
    for (uint64_t i = 0; i < workers; i++) {
        pthread_join(th[i], NULL);
        hits += jobs[i].hits; dsum += jobs[i].distinct_sum; aborts += jobs[i].aborts;
        ksum += jobs[i].keypairs_sum; rsum += jobs[i].respairs_sum;
        if (jobs[i].maxdup > maxdup) maxdup = jobs[i].maxdup;
        if (jobs[i].max_keypairs > kmax) kmax = jobs[i].max_keypairs;
        msum += jobs[i].mixpairs_sum;
        if (jobs[i].max_mixpairs > mmax) mmax = jobs[i].max_mixpairs;
    }
    // merge per-walk pair arrays for exact percentiles
    uint64_t tot = 0; for (uint64_t i = 0; i < workers; i++) tot += jobs[i].trials;
    uint64_t *kp = malloc(8 * (tot ? tot : 1)), *rp = malloc(8 * (tot ? tot : 1)), o = 0;
    uint64_t *mp = malloc(8 * (tot ? tot : 1));
    for (uint64_t i = 0; i < workers; i++) {
        memcpy(kp + o, jobs[i].kpair_walk, 8 * jobs[i].trials);
        memcpy(rp + o, jobs[i].rpair_walk, 8 * jobs[i].trials);
        memcpy(mp + o, jobs[i].mpair_walk, 8 * jobs[i].trials);
        o += jobs[i].trials;
    }
    qsort(kp, tot, 8, cmp_u64); qsort(rp, tot, 8, cmp_u64); qsort(mp, tot, 8, cmp_u64);
    qsort(GAPS, NGAP, 8, cmp_u64);

    double done = (double)(tot - aborts);
    double NN = (double)(nb * ell);
    double kexp = NN * (NN - 1.0) / 2.0 / (double)(KEY_BITS ? (1ULL << KEY_BITS) : M_H);
    double rexp = NN * (NN - 1.0) / (2.0 * (double)M_H);
    double gmean = (double)GAP_SUM / (double)(NGAP ? NGAP : 1);
    double mu_fix = 8.0 * log((double)W_MAX);   // fixture density law: d=8, q ~ W_MAX
    uint64_t *mg = malloc(8 * (tot ? tot : 1)), mo = 0;
    for (uint64_t i = 0; i < workers; i++) { memcpy(mg + mo, jobs[i].mg_walk, 8 * jobs[i].trials); mo += jobs[i].trials; }
    qsort(mg, tot, 8, cmp_u64);
    printf("{\"p\":%llu,\"m\":%llu,\"log2_m\":%.6f,\"nb\":%llu,\"ell\":%llu,\"span_bits\":%llu,"
           "\"key_bits\":%llu,\"gap_cap\":%llu,\"walk_len\":%llu,\"walk_band_log2\":%.4f,"
           "\"records_per_side\":%.0f,\"K\":%.9f,\"trials\":%llu,\"workers\":%llu,\"hits\":%llu,\"q\":%.9f,"
           "\"aborts\":%llu,\"abort_rate\":%.6f,"
           "\"gap_mean\":%.4f,\"gap_p50\":%.0f,\"gap_p99\":%.0f,\"gap_p999\":%.0f,\"gap_max\":%llu,"
           "\"gap_mean_over_mu\":%.6f,"
           "\"trail_maxgap_p50\":%.0f,\"trail_maxgap_p95\":%.0f,\"trail_maxgap_p99\":%.0f,"
           "\"trail_maxgap_p999\":%.0f,\"trail_maxgap_max\":%.0f,\"trail_maxgap_p50_over_mu\":%.6f,"
           "\"key_pair_expect_per_walk\":%.6f,\"key_pairs_mean\":%.6f,"
           "\"kappa_key_mean\":%.6f,\"kappa_key_p50\":%.6f,\"kappa_key_p95\":%.6f,\"kappa_key_p99\":%.6f,"
           "\"kappa_key_max\":%.6f,\"key_pairs_max\":%llu,"
           "\"key_pairs_p99\":%.0f,\"key_pairs_p999\":%.0f,"
           "\"res_pair_expect_per_walk\":%.6f,\"res_pairs_mean\":%.6f,"
           "\"kappa_res_mean\":%.6f,\"kappa_res_max\":%.6f,\"max_dup_records\":%llu,"
           "\"mix_pair_expect_per_walk\":%.6f,\"mix_pairs_mean\":%.6f,\"mix_pairs_max\":%llu,"
           "\"kappa_mix_mean\":%.6f,\"kappa_mix_p50\":%.6f,\"kappa_mix_p95\":%.6f,"
           "\"kappa_mix_p99\":%.6f,\"kappa_mix_max\":%.6f,"
           "\"span_plus_band\":%llu,\"exp_no_reduce\":%s,\"log2_m2\":%.6f,\"wall_seconds\":%.2f}\n",
           (unsigned long long)P, (unsigned long long)M_H, log2((double)M_H),
           (unsigned long long)nb, (unsigned long long)ell, (unsigned long long)SPAN_BITS,
           (unsigned long long)KEY_BITS, (unsigned long long)GAP_CAP, (unsigned long long)W_LEN,
           log2((double)W_MAX), N, 6.0 * N * N / (double)M_H,
           (unsigned long long)tot, (unsigned long long)workers, (unsigned long long)hits,
           done > 0 ? (double)hits / done : 0.0,
           (unsigned long long)aborts, tot ? (double)aborts / (double)tot : 0.0,
           gmean, pct(GAPS, NGAP, 0.50), pct(GAPS, NGAP, 0.99), pct(GAPS, NGAP, 0.999),
           (unsigned long long)GAP_MAXV, gmean / mu_fix,
           pct(mg, tot, 0.50), pct(mg, tot, 0.95), pct(mg, tot, 0.99), pct(mg, tot, 0.999),
           pct(mg, tot, 1.0), pct(mg, tot, 0.50) / mu_fix,
           kexp, done > 0 ? (double)ksum / done : 0.0,
           kexp > 0 ? (double)ksum / done / kexp : 0.0, kexp > 0 ? pct(kp, tot, 0.50) / kexp : 0.0,
           kexp > 0 ? pct(kp, tot, 0.95) / kexp : 0.0, kexp > 0 ? pct(kp, tot, 0.99) / kexp : 0.0,
           kexp > 0 ? (double)kmax / kexp : 0.0, (unsigned long long)kmax,
           pct(kp, tot, 0.99), pct(kp, tot, 0.999),
           rexp, done > 0 ? (double)rsum / done : 0.0,
           rexp > 0 ? (double)rsum / done / rexp : 0.0, rexp > 0 ? (double)maxdup / rexp : 0.0,
           (unsigned long long)maxdup,
           kexp, done > 0 ? (double)msum / done : 0.0, (unsigned long long)mmax,
           kexp > 0 ? (double)msum / done / kexp : 0.0, kexp > 0 ? pct(mp, tot, 0.50) / kexp : 0.0,
           kexp > 0 ? pct(mp, tot, 0.95) / kexp : 0.0, kexp > 0 ? pct(mp, tot, 0.99) / kexp : 0.0,
           kexp > 0 ? (double)mmax / kexp : 0.0,
           (unsigned long long)(SPAN_BITS + (uint64_t)ceil(log2((double)W_MAX))),
           ((double)(SPAN_BITS) + ceil(log2((double)W_MAX)) < log2((double)M_H)) ? "true" : "false",
           log2((double)M_H),
           (double)time(NULL) - w0);
    for (uint64_t i = 0; i < workers; i++) { free(jobs[i].kpair_walk); free(jobs[i].rpair_walk);
        free(jobs[i].mg_walk); free(jobs[i].mpair_walk); }
    free(kp); free(rp); free(mp); free(mg); free(jobs); free(th);
    return 0;
}
