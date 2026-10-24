// v28 STRUCTURED-WALK TRANSFER ARM
// Samples the PRODUCTION object, not the ratio layer:
//   records(side) = { a_i * w_k mod p }  with a SHARED walk list W (same for both
//   sides AND across all bases), W = consecutive H-primes (index-8 subgroup,
//   prime integers, genuine gaps), a_i,b_j uniform in H, target v uniform in H.
//   hit  <=>  exists u in mu6, i,j,k,l : v*u*b_j*w_l == a_i*w_k  (mod p)
// Density matched to production: K = 6*N^2/m, m = |H| = (p-1)/8.
// (nb,ell) LADDER: production is nb=2^66 bases x ell=2^49.65 steps (nb/ell=2^16.35).
// Reports hits, distinct records (set-law), and duplicate max (kappa table feed).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <pthread.h>
#include <time.h>

typedef unsigned __int128 u128;
static uint64_t P, M_H, G, MU6[6];

static inline uint64_t mulmod(uint64_t a, uint64_t b) { return (uint64_t)((u128)a * b % P); }
static uint64_t powmod(uint64_t b, uint64_t e) {
    uint64_t r = 1; b %= P;
    while (e) { if (e & 1) r = mulmod(r, b); b = mulmod(b, b); e >>= 1; }
    return r;
}
static inline uint64_t pow8(uint64_t x) { uint64_t y = mulmod(x % P, x % P); y = mulmod(y, y); return mulmod(y, y); }
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

static uint64_t *W_LIST; static uint64_t W_LEN;

static void build_walk(uint64_t len, uint64_t lo) {
    W_LIST = malloc(sizeof(uint64_t) * len);
    uint64_t q = lo | 1, got = 0;
    while (got < len) {
        uint64_t qq = q % P;
        if (qq > 1 && is_prime_u64(q) && powmod(qq, M_H) == 1) W_LIST[got++] = qq;
        q += 2;
    }
    W_LEN = got;
}

typedef struct { uint64_t *key; uint32_t *tag; uint32_t epoch; uint64_t mask; } HS;

static void hs_init(HS *h, uint64_t n) {
    uint64_t b = 1; while ((1ULL << b) < 2 * n) b++;
    h->mask = (1ULL << b) - 1;
    h->key = calloc(h->mask + 1, 8); h->tag = calloc(h->mask + 1, 4); h->epoch = 0;
}
static void hs_free(HS *h) { free(h->key); free(h->tag); }
static inline void hs_clear(HS *h) { h->epoch++; if (h->epoch == 0) { memset(h->tag, 0, (h->mask + 1) * 4); h->epoch = 1; } }
static inline void hs_put(HS *h, uint64_t k) {
    uint64_t i = (k * 0x9E3779B97F4A7C15ULL) & h->mask;
    while (h->tag[i] == h->epoch) { if (h->key[i] == k) return; i = (i + 1) & h->mask; }
    h->tag[i] = h->epoch; h->key[i] = k;
}
static inline int hs_ins_new(HS *h, uint64_t k) {   // 1 if newly inserted
    uint64_t i = (k * 0x9E3779B97F4A7C15ULL) & h->mask;
    while (h->tag[i] == h->epoch) { if (h->key[i] == k) return 0; i = (i + 1) & h->mask; }
    h->tag[i] = h->epoch; h->key[i] = k; return 1;
}
static inline int hs_has(HS *h, uint64_t k) {
    uint64_t i = (k * 0x9E3779B97F4A7C15ULL) & h->mask;
    while (h->tag[i] == h->epoch) { if (h->key[i] == k) return 1; i = (i + 1) & h->mask; }
    return 0;
}

typedef struct {
    int wid; uint64_t nb, ell, trials, seed;
    uint64_t hits, distinct_sum, maxdup; double seconds;
} Job;

static void *worker(void *arg) {
    Job *j = (Job *)arg;
    uint64_t N = j->nb * j->ell;
    uint64_t *A = malloc(sizeof(uint64_t) * N), *B = malloc(sizeof(uint64_t) * N);
    HS h, hd; hs_init(&h, 2 * N); hs_init(&hd, 2 * N);
    uint64_t s = j->seed * 0x2545F4914F6CDD1DULL + 12345;
    double t0 = (double)clock() / CLOCKS_PER_SEC;
    for (uint64_t t = 0; t < j->trials; t++) {
        hs_clear(&h); hs_clear(&hd);
        for (uint64_t i = 0; i < j->nb; i++) {
            uint64_t a = pow8(rnd(&s) | 1), b = pow8(rnd(&s) | 1);
            for (uint64_t k = 0; k < j->ell; k++) {
                A[i * j->ell + k] = mulmod(a, W_LIST[k]);
                B[i * j->ell + k] = mulmod(b, W_LIST[k]);
            }
        }
        uint64_t nd = 0;
        for (uint64_t x = 0; x < N; x++) { if (hs_ins_new(&hd, A[x])) nd++; hs_put(&h, A[x]); }
        uint64_t v = pow8(rnd(&s) | 1);
        int hit = 0;
        for (int u = 0; u < 6 && !hit; u++) {
            uint64_t vu = mulmod(v, MU6[u]);
            for (uint64_t x = 0; x < N; x++) if (hs_has(&h, mulmod(vu, B[x]))) { hit = 1; break; }
        }
        j->hits += hit; j->distinct_sum += nd;
        if (N - nd > j->maxdup) j->maxdup = N - nd;
    }
    j->seconds = (double)clock() / CLOCKS_PER_SEC - t0;
    free(A); free(B); hs_free(&h); hs_free(&hd);
    return NULL;
}

int main(int argc, char **argv) {
    uint64_t nb = strtoull(argv[1], 0, 10), ell = strtoull(argv[2], 0, 10);
    uint64_t trials = strtoull(argv[3], 0, 10), workers = strtoull(argv[4], 0, 10);
    uint64_t seed0 = argc > 5 ? strtoull(argv[5], 0, 10) : 1;
    uint64_t p = (1ULL << 41) + 1;
    while (!(is_prime_u64(p) && (p - 1) % 48 == 0)) p += 2;
    P = p; M_H = (p - 1) / 8;
    for (uint64_t g = 2; g < 100000; g++) {
        if (powmod(g, (p - 1) / 2) != 1 && powmod(g, (p - 1) / 3) != 1) { G = g; break; }
    }
    for (int j = 0; j < 6; j++) MU6[j] = powmod(G, (p - 1) / 6 * j);
    build_walk(ell, 1ULL << 20);
    double N = (double)nb * (double)ell;
    Job *jobs = calloc(workers, sizeof(Job));
    pthread_t *th = calloc(workers, sizeof(pthread_t));
    uint64_t base = trials / workers;
    double t0 = (double)clock() / CLOCKS_PER_SEC, w0 = time(NULL);
    for (uint64_t i = 0; i < workers; i++) {
        jobs[i].wid = i; jobs[i].nb = nb; jobs[i].ell = ell;
        jobs[i].trials = base + (i < trials - base * workers ? 1 : 0);
        jobs[i].seed = seed0 * 1000 + i;
        pthread_create(&th[i], NULL, worker, &jobs[i]);
    }
    uint64_t hits = 0, dsum = 0, maxdup = 0;
    double secs = 0;
    for (uint64_t i = 0; i < workers; i++) {
        pthread_join(th[i], NULL);
        hits += jobs[i].hits; dsum += jobs[i].distinct_sum; secs += jobs[i].seconds;
        if (jobs[i].maxdup > maxdup) maxdup = jobs[i].maxdup;
    }
    printf("{\"p\":%llu,\"m\":%llu,\"log2_m\":%.6f,\"nb\":%llu,\"ell\":%llu,\"records_per_side\":%.0f,"
           "\"K\":%.9f,\"trials\":%llu,\"workers\":%llu,\"walk_len\":%llu,\"hits\":%llu,\"q\":%.9f,"
           "\"distinct_mean\":%.3f,\"max_duplicate_records\":%llu,\"core_seconds\":%.2f,\"wall_seconds\":%.2f}\n",
           (unsigned long long)P, (unsigned long long)M_H, log2((double)M_H),
           (unsigned long long)nb, (unsigned long long)ell, N, 6.0 * N * N / (double)M_H,
           (unsigned long long)trials, (unsigned long long)workers, (unsigned long long)W_LEN,
           (unsigned long long)hits, (double)hits / (double)trials,
           (double)dsum / (double)trials, (unsigned long long)maxdup,
           secs, (double)time(NULL) - w0);
    (void)t0;
    return 0;
}
