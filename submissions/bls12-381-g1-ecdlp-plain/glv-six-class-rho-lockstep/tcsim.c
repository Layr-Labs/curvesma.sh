// Offline toy simulator: lockstep six-class rho (production rules of proof.md
// section 2 at toy scale) measuring T_c = class arrivals of all trails until the
// first arrival at a class already visited by another trail.
// usage: tcsim q r h R dp_bits M k trials seed
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
static u64 Q, RR, HH, BETA, LAM;
static u64 mulm(u64 a, u64 b) { return (u64)((u128)a * b % Q); }
static u64 addm(u64 a, u64 b) { u64 s = a + b; return s >= Q ? s - Q : s; }
static u64 subm(u64 a, u64 b) { return a >= b ? a - b : a + Q - b; }
static u64 powm(u64 a, u64 e, u64 m) { u128 r = 1, b = a % m; while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; } return (u64)r; }
static u64 invm(u64 a) { // extended Euclid
    int64_t t = 0, nt = 1; u64 r = Q, nr = a;
    while (nr) { u64 qq = r / nr; int64_t tt = t - (int64_t)qq * nt; t = nt; nt = tt; u64 rr = r - qq * nr; r = nr; nr = rr; }
    return t < 0 ? (u64)(t + (int64_t)Q) : (u64)t;
}
typedef struct { u64 x, y; int inf; } Pt;
static Pt padd(Pt P, Pt R2) {
    if (P.inf) return R2; if (R2.inf) return P;
    u64 l;
    if (P.x == R2.x) {
        if (addm(P.y, R2.y) == 0) { Pt o = {0, 0, 1}; return o; }
        l = mulm(mulm(3, mulm(P.x, P.x)), invm(addm(P.y, P.y)));
    } else l = mulm(subm(R2.y, P.y), invm(subm(R2.x, P.x)));
    Pt o; o.inf = 0; o.x = subm(subm(mulm(l, l), P.x), R2.x); o.y = subm(mulm(l, subm(P.x, o.x)), P.y); return o;
}
static Pt pmul(u64 k, Pt P) { Pt Rr = {0, 0, 1}; while (k) { if (k & 1) Rr = padd(Rr, P); P = padd(P, P); k >>= 1; } return Rr; }
static u64 rs = 88172645463325252ULL;
static u64 rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static u64 rndr(u64 n) { return rnd() % n; }
static Pt FB[8][256]; // fixed base [j*256^i]G
static Pt fmul(u64 k) { Pt Rr = {0, 0, 1}; int i = 0; k %= RR; while (k) { Rr = padd(Rr, FB[i][k & 255]); k >>= 8; i++; } return Rr; }
static u64 sqrtm(u64 a) { // Tonelli-Shanks
    if (a == 0) return 0; if (powm(a, (Q - 1) / 2, Q) != 1) return (u64)-1;
    u64 s = 0, t = Q - 1; while (!(t & 1)) { t >>= 1; s++; }
    u64 z = 2; while (powm(z, (Q - 1) / 2, Q) != Q - 1) z++;
    u64 m = s, c = powm(z, t, Q), tt = powm(a, t, Q), x = powm(a, (t + 1) / 2, Q);
    while (tt != 1) { u64 i = 1, t2 = mulm(tt, tt); while (t2 != 1) { t2 = mulm(t2, t2); i++; } u64 b = c; for (u64 j = 0; j < m - i - 1; j++) b = mulm(b, b); m = i; c = mulm(b, b); tt = mulm(tt, c); x = mulm(x, b); }
    return x;
}
// class arrival hash: key K (even, < Q) -> trail id
static u64 *hk; static u64 *hv; static u64 hmask;
static u64 hset_get(u64 K, u64 tid) { // returns previous tid or 0 (inserted)
    u64 i = (K * 0x9E3779B97F4A7C15ULL) >> 20 & hmask;
    while (hk[i]) { if (hk[i] == K + 1) return hv[i]; i = (i + 1) & hmask; }
    hk[i] = K + 1; hv[i] = tid; return 0;
}
static int RB; static u64 Rn, DPM;
typedef struct { Pt P; u64 ts; u64 kcp; int kv; u64 tid; } W;
static Pt *TAB; // R*6
static void canon(Pt P, u64 *K, int *t) {
    u64 bx = mulm(BETA, P.x), b2x = subm(subm(0, P.x), bx);
    int j = 0; u64 cx = P.x; if (bx < cx) { cx = bx; j = 1; } if (b2x < cx) { cx = b2x; j = 2; }
    if (P.y & 1) { *K = Q - P.y; *t = 3 + j; } else { *K = P.y; *t = j; }
}
int main(int argc, char **argv) {
    if (argc < 10) { fprintf(stderr, "usage\n"); return 1; }
    Q = strtoull(argv[1], 0, 10); u64 n = strtoull(argv[2], 0, 10); HH = strtoull(argv[3], 0, 10); RR = n / HH;
    Rn = strtoull(argv[4], 0, 10); int dpb = atoi(argv[5]); int M = atoi(argv[6]); int k = atoi(argv[7]);
    long trials = atol(argv[8]); rs ^= strtoull(argv[9], 0, 10) * 0x9E3779B97F4A7C15ULL; for (int i = 0; i < 10; i++) rnd();
    double kappa = argc > 10 ? atof(argv[10]) : 1.005;
    RB = 0; while ((1ULL << RB) < Rn) RB++; DPM = (1ULL << dpb) - 1;
    u64 Pref = 64 * (u64)k, Lmax = 20ULL << dpb;
    // generator
    Pt G;
    for (;;) { u64 x = rndr(Q); u64 y = sqrtm(addm(mulm(mulm(x, x), x), 4)); if (y == (u64)-1) continue; Pt P = {x, y, 0}; G = pmul(HH, P); if (!G.inf && pmul(RR, G).inf) break; }
    for (int i = 0; i < 8; i++) { Pt B = i == 0 ? G : padd(FB[i - 1][255], FB[i - 1][1]); FB[i][0].inf = 1; for (int j = 1; j < 256; j++) FB[i][j] = padd(FB[i][j - 1], B); }
    // beta / lambda matching
    u64 bc[2], lc[2]; { u64 g = 2; for (;; g++) { u64 c = powm(g, (Q - 1) / 3, Q); if (c != 1) { bc[0] = c; bc[1] = mulm(c, c); break; } } }
    { u64 g = 2; for (;; g++) { u64 c = powm(g, (RR - 1) / 3, RR); if (c != 1) { lc[0] = c; lc[1] = (u64)((u128)c * c % RR); break; } } }
    for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) { Pt L = fmul(lc[b]); if (L.x == mulm(G.x, bc[a]) && L.y == G.y) { BETA = bc[a]; LAM = lc[b]; } }
    u64 binv = mulm(BETA, BETA);
    double N = (double)(RR - 1) / 6.0, Np = N * Rn / (Rn - 1.0);
    double S0 = sqrt(2 * Np * log(1 / 0.6097)), S1 = kappa * S0;
    u64 hs = 1; while (hs < (u64)(16 * S1 + 1024)) hs <<= 1; hmask = hs - 1;
    hk = calloc(hs, 8); hv = calloc(hs, 8); TAB = malloc(sizeof(Pt) * Rn * 6);
    W *ws = malloc(sizeof(W) * M);
    double sum_tc = 0, sum_tc2 = 0; long le1 = 0, le0 = 0, done = 0; long long det = 0, ops_all = 0;
    for (long tr = 0; tr < trials; tr++) {
        for (u64 h = 0; h < Rn; h++) {
            Pt Mh; do Mh = fmul(1 + rndr(RR - 1)); while (Mh.inf);
            for (int t = 0; t < 6; t++) { int s = t / 3, j = t % 3; Pt Mt = Mh; u64 f = j == 0 ? 1 : (j == 1 ? binv : BETA); Mt.x = mulm(Mh.x, f); if (s) Mt.y = subm(0, Mt.y); TAB[h * 6 + t] = Mt; }
        }
        memset(hk, 0, hs * 8);
        u64 arrivals = 0, tc = 0, tidc = 0;
#define ARRIVE(KK, TID) do { arrivals++; u64 pv = hset_get((KK), (TID)); if (pv && pv != (TID) && !tc) tc = arrivals; } while (0)
        for (int i = 0; i < M; i++) { Pt P; do P = fmul(rndr(RR)); while (P.inf); ws[i].P = P; ws[i].ts = 0; ws[i].kv = 0; ws[i].tid = ++tidc; u64 K; int t; canon(P, &K, &t); ARRIVE(K, ws[i].tid); }
        for (u64 bt = 0; !tc; bt++) {
            for (int i = 0; i < M && !tc; i++) {
                W *w = &ws[i]; u64 K; int t; Pt tgt;
                for (;;) {
                    canon(w->P, &K, &t);
                    if (((K >> (1 + RB)) & DPM) == 0) { // distinguished: restart
                        Pt P; do P = fmul(rndr(RR)); while (P.inf); w->P = P; w->ts = bt; w->kv = 0; w->tid = ++tidc; canon(P, &K, &t); ARRIVE(K, w->tid); if (tc) break; continue; }
                    if (bt % k == 0 && w->kv && w->kcp == K) {
                        det++;
                        u64 bestK = K; Pt best = w->P, Qp = w->P; u64 cnt = 0; int ok = 1;
                        for (;;) { u64 KQ; int tq; canon(Qp, &KQ, &tq); Pt T2 = TAB[((KQ >> 1) & (Rn - 1)) * 6 + tq]; if (Qp.x == T2.x) { ok = 0; break; } Qp = padd(Qp, T2); canon(Qp, &KQ, &tq); ARRIVE(KQ, w->tid); if (KQ == K) break; if (KQ < bestK) { bestK = KQ; best = Qp; } if (++cnt > Pref) { ok = 0; break; } }
                        if (tc) break;
                        if (!ok) { Pt P; do P = fmul(rndr(RR)); while (P.inf); w->P = P; w->ts = bt; w->kv = 0; w->tid = ++tidc; canon(P, &K, &t); ARRIVE(K, w->tid); if (tc) break; continue; }
                        w->P = padd(best, best); w->kv = 0; canon(w->P, &K, &t); ARRIVE(K, w->tid); if (tc) break; continue;
                    }
                    if (bt % Pref == 0) {
                        if (bt - w->ts >= Lmax) { Pt P; do P = fmul(rndr(RR)); while (P.inf); w->P = P; w->ts = bt; w->kv = 0; w->tid = ++tidc; canon(P, &K, &t); ARRIVE(K, w->tid); if (tc) break; continue; }
                        w->kcp = K; w->kv = 1;
                    }
                    tgt = TAB[((K >> 1) & (Rn - 1)) * 6 + t];
                    break;
                }
                if (tc) break;
                if (w->P.x == tgt.x) { Pt P; do P = fmul(rndr(RR)); while (P.inf); w->P = P; w->ts = bt + 1; w->kv = 0; w->tid = ++tidc; canon(P, &K, &t); ARRIVE(K, w->tid); continue; }
                w->P = padd(w->P, tgt); ops_all++; canon(w->P, &K, &t); ARRIVE(K, w->tid);
            }
        }
        double x = tc / sqrt(Np); sum_tc += x; sum_tc2 += x * x; done++;
        if (tc <= S1) le1++; if (tc <= S0) le0++;
    }
    double mean = sum_tc / done, sd = sqrt((sum_tc2 / done - mean * mean) / done);
    double model = sqrt(M_PI / 2);
    printf("{\"q\":%llu,\"R\":%llu,\"dp\":%d,\"M\":%d,\"k\":%d,\"trials\":%ld,\"trail_over_sqrtN\":%.6f,\"mean_tc_over_sqrtNp\":%.5f,\"se\":%.5f,\"model\":%.5f,\"slowdown_pct\":%.3f,\"slowdown_se_pct\":%.3f,\"frac_le_S1\":%.5f,\"frac_le_S0\":%.5f,\"model_S1\":%.5f,\"det_per_op\":%.3g}\n",
           (unsigned long long)Q, (unsigned long long)Rn, dpb, M, k, done, (double)(1ULL << dpb) / sqrt(N), mean, sd, model, 100 * (mean / model - 1), 100 * sd / model,
           (double)le1 / done, (double)le0 / done, 1 - exp(-S1 * S1 / (2 * Np)), (double)det / (ops_all ? ops_all : 1));
    return 0;
}
