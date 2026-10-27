/* ecm12: the modular multiplication behind the cofactorisation price (stnfs-aslice-v12).
 *
 * Survivor cofactors are at most 2^300 (f slack <= 292 bits incl. rounding; g smaller), so ECM arithmetic uses
 * L = 5 64-bit limbs.  mont_mul is CIOS Montgomery multiplication (a b 2^-320 mod n, inputs and output < n, n odd,
 * n < 2^319 so the final conditional subtraction suffices), 64-bit general-register arithmetic with 128-bit products.
 * ECM curve operations are counted in these multiplications (S = M); proof.md section 6 gives the curve counts.
 *
 * Modes:
 *   check n          n random products compared with a slow exact reference (shift-subtract reduction)
 *   bench reps       reps chained multiplications (whole-program Ir at two values of reps)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define L 5
typedef unsigned __int128 u128;
static void mont_mul(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *n, uint64_t ninv) {
    uint64_t t[L + 2]; memset(t, 0, sizeof t);
    for (int i = 0; i < L; i++) {
        u128 c = 0;
        for (int j = 0; j < L; j++) { c += (u128)a[j] * b[i] + t[j]; t[j] = (uint64_t)c; c >>= 64; }
        c += t[L]; t[L] = (uint64_t)c; t[L + 1] = (uint64_t)(c >> 64);
        uint64_t m = t[0] * ninv; c = (u128)m * n[0] + t[0]; c >>= 64;
        for (int j = 1; j < L; j++) { c += (u128)m * n[j] + t[j]; t[j - 1] = (uint64_t)c; c >>= 64; }
        c += t[L]; t[L - 1] = (uint64_t)c; t[L] = t[L + 1] + (uint64_t)(c >> 64);
    }
    /* t < 2n: subtract n once if t >= n */
    uint64_t d[L], br = 0;
    for (int j = 0; j < L; j++) { u128 x = (u128)t[j] - n[j] - br; d[j] = (uint64_t)x; br = (uint64_t)(x >> 64) & 1; }
    int ge = t[L] || !br;
    for (int j = 0; j < L; j++) r[j] = ge ? d[j] : t[j];
}
static uint64_t inv64(uint64_t n0) { uint64_t x = 1; for (int i = 0; i < 7; i++) x *= 2 - n0 * x; return -x; }   /* -n^-1 mod 2^64 */
/* slow exact reference: (a b) mod n by shift-subtract on 2L+1 limbs, then multiply by 2^-320 via repeated halving */
static int geq(const uint64_t *x, const uint64_t *y, int k) { for (int i = k - 1; i >= 0; i--) { if (x[i] != y[i]) return x[i] > y[i]; } return 1; }
static void sub(uint64_t *x, const uint64_t *y, int k) { uint64_t br = 0; for (int i = 0; i < k; i++) { u128 v = (u128)x[i] - y[i] - br; x[i] = (uint64_t)v; br = (uint64_t)(v >> 64) & 1; } }
static void ref(uint64_t *r, const uint64_t *a, const uint64_t *b, const uint64_t *n) {
    uint64_t p[2 * L + 1]; memset(p, 0, sizeof p);
    for (int i = 0; i < L; i++) { u128 c = 0; for (int j = 0; j < L; j++) { c += (u128)a[i] * b[j] + p[i + j]; p[i + j] = (uint64_t)c; c >>= 64; } p[i + L] = (uint64_t)c; }
    /* reduce p mod n bit by bit */
    uint64_t rem[L + 1]; memset(rem, 0, sizeof rem); uint64_t nn[L + 1]; memcpy(nn, n, 8 * L); nn[L] = 0;
    for (int bit = 64 * 2 * L - 1; bit >= 0; bit--) {
        for (int i = L; i > 0; i--) rem[i] = (rem[i] << 1) | (rem[i - 1] >> 63); rem[0] = (rem[0] << 1) | ((p[bit / 64] >> (bit % 64)) & 1);
        if (geq(rem, nn, L + 1)) sub(rem, nn, L + 1);
    }
    /* divide by 2^320 mod n: halve 320 times (add n when odd) */
    for (int k = 0; k < 64 * L; k++) {
        uint64_t c = 0; if (rem[0] & 1) { u128 s = 0; for (int i = 0; i <= L; i++) { s += (u128)rem[i] + nn[i]; rem[i] = (uint64_t)s; s >>= 64; } }
        for (int i = 0; i <= L; i++) { uint64_t nc = rem[i] & 1; rem[i] = (rem[i] >> 1) | (i < L ? rem[i + 1] << 63 : 0); (void)nc; } (void)c;
    }
    memcpy(r, rem, 8 * L);
}
static uint64_t st = 0x2545F4914F6CDD1DULL;
static uint64_t rnd(void) { st ^= st << 13; st ^= st >> 7; st ^= st << 17; return st; }
static void rand_mod(uint64_t *x, const uint64_t *n) { do { for (int i = 0; i < L; i++) x[i] = rnd(); x[L - 1] &= n[L - 1] | (n[L - 1] >> 1) | 0x7FFFFFFFFFFFFFFFULL; } while (geq(x, n, L)); }
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: ecm12 check n | bench reps\n"); return 2; }
    uint64_t n[L]; for (int i = 0; i < L; i++) n[i] = rnd(); n[0] |= 1; n[L - 1] = (n[L - 1] >> 2) | (1ULL << 61);   /* odd, < 2^319 */
    uint64_t ninv = inv64(n[0]);
    if (!strcmp(argv[1], "check")) {
        long long cnt = atoll(argv[2]), bad = 0;
        for (long long t = 0; t < cnt; t++) {
            if (t % 100 == 0) { for (int i = 0; i < L; i++) n[i] = rnd(); n[0] |= 1; n[L - 1] = (n[L - 1] >> (2 + t % 7)) | (1ULL << (61 - t % 7)); ninv = inv64(n[0]); }
            uint64_t a[L], b[L], r1[L], r2[L]; rand_mod(a, n); rand_mod(b, n);
            if (t % 10 == 1) { memcpy(a, n, sizeof a); a[0] -= 1; memcpy(b, a, sizeof b); }   /* (n-1)^2 */
            mont_mul(r1, a, b, n, ninv); ref(r2, a, b, n); bad += memcmp(r1, r2, sizeof r1) != 0;
        }
        printf("mont_mul (5 limbs): %lld products, %lld mismatches\n", cnt, bad); return bad != 0;
    }
    long long reps = atoll(argv[2]); uint64_t x[L], y[L]; rand_mod(x, n); rand_mod(y, n);
    for (long long t = 0; t < reps; t++) { mont_mul(x, x, y, n, ninv); __asm__ __volatile__("" :: "r"(x) : "memory"); }
    printf("reps=%lld x0=%llu\n", reps, (unsigned long long)x[0]); return 0;
}
