/* Projection through X = U*D, field-correct Montgomery arithmetic mod r.
   D holds the diagonal in Montgomery form (d*2^256 mod r).  Per column element: one
   Montgomery multiply v*d*2^-256 mod r (4x4 schoolbook + REDC), then one signed add into each
   of the column's 3 accumulator rows with a conditional subtraction of r.  U: 3 rows per
   column, uniform signs, branch-free two's complement. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define M 28
static const uint64_t R[4] = {0xffffffff00000001ULL, 0x53bda402fffe5bfeULL, 0x3339d80809a1d805ULL, 0x73eda753299d7d48ULL};
static const uint64_t R0INV = 0xfffffffeffffffffULL;     /* -r^-1 mod 2^64 */
static const uint64_t R2[4] = {0xc999e990f3f29c6dULL, 0x2b6cedcb87925c23ULL, 0x05d314967254398fULL, 0x0748d9d99f59ff11ULL};
static uint64_t V[1<<20][4], ACC[1<<12][4], D[1<<20][4];
static uint32_t ROW[1<<20][3]; static int8_t SG[1<<20];

static inline void mul4(uint64_t out[8], const uint64_t a[4], const uint64_t b[4]) {
  for (int i = 0; i < 8; i++) out[i] = 0;
  for (int i = 0; i < 4; i++) {
    unsigned __int128 carry = 0;
    for (int j = 0; j < 4; j++) {
      unsigned __int128 z = (unsigned __int128)a[j] * b[i] + out[i + j] + carry;
      out[i + j] = (uint64_t)z; carry = z >> 64;
    }
    int k = i + 4; uint64_t cy = (uint64_t)carry;
    while (cy) { unsigned __int128 w = (unsigned __int128)out[k] + cy; out[k] = (uint64_t)w; cy = (uint64_t)(w >> 64); k++; }
  }
}
/* REDC of an 8-limb value: returns the 4-limb result (may be >= r by < r once) */
static inline void redc(uint64_t t[8], uint64_t out[4]) {
  unsigned __int128 c = 0;
  for (int i = 0; i < 4; i++) {
    uint64_t m = t[i] * R0INV;
    c = 0;
    for (int j = 0; j < 4; j++) {
      unsigned __int128 z = (unsigned __int128)m * R[j] + t[i + j] + c;
      t[i + j] = (uint64_t)z; c = z >> 64;
    }
    unsigned __int128 z = (unsigned __int128)t[i + 4] + c; t[i + 4] = (uint64_t)z;
    int k = i + 5; uint64_t cy = (uint64_t)(z >> 64);
    while (cy && k < 8) { unsigned __int128 w = (unsigned __int128)t[k] + cy; t[k] = (uint64_t)w; cy = (uint64_t)(w >> 64); k++; }
  }
  for (int i = 0; i < 4; i++) out[i] = t[i + 4];
}
static inline void condsub(uint64_t a[4]) {
  int ge = 1;
  for (int i = 3; i >= 0; i--) { if (a[i] != R[i]) { ge = a[i] > R[i]; break; } }
  if (!ge) return;
  unsigned char b = 0;
  for (int i = 0; i < 4; i++) { unsigned __int128 z = (unsigned __int128)a[i] - R[i] - b; a[i] = (uint64_t)z; b = (unsigned char)((z >> 64) ? 1 : 0); }
}
static inline void sub4(uint64_t out[4], const uint64_t x[4]) {   /* r - x, x != 0 */
  unsigned char b = 0;
  for (int i = 0; i < 4; i++) { unsigned __int128 z = (unsigned __int128)R[i] - x[i] - b; out[i] = (uint64_t)z; b = (unsigned char)((z >> 64) ? 1 : 0); }
}
static inline void acc4(uint64_t a[4], const uint64_t p[4], uint64_t mask) {
  uint64_t q[4];
  if (mask) { sub4(q, p); if (!(q[0]|q[1]|q[2]|q[3])) { /* p = 0 */ } }
  else { q[0]=p[0]; q[1]=p[1]; q[2]=p[2]; q[3]=p[3]; }
  unsigned char c = 0;
  for (int i = 0; i < 4; i++) { unsigned __int128 z = (unsigned __int128)a[i] + q[i] + c; a[i] = (uint64_t)z; c = (unsigned char)(z >> 64); }
  condsub(a);
}
int main(int argc, char **argv) {
  int reps = (argc > 1 && argv[1][0] >= '0' && argv[1][0] <= '9') ? atoi(argv[1]) : (argc > 2 ? atoi(argv[2]) : 0);
  size_t N = 1 << 20;
  if (argc > 1 && !strcmp(argv[1], "selftest")) {
    uint64_t p[8], o[4], t[4], am[4], bm[4], m[4];
    uint64_t a[4] = {3,0,0,0}, b[4] = {5,0,0,0}, one[4] = {1,0,0,0};
    mul4(p, a, R2); redc(p, am);           /* mont(3) */
    mul4(p, b, R2); redc(p, bm);           /* mont(5) */
    mul4(p, am, bm); redc(p, m);           /* mont(15) */
    mul4(p, m, one); redc(p, t);           /* unmont -> 15 */
    printf("mont(3)*mont(5) unmonted = %llu (expect 15)\n", (unsigned long long)t[0]);
    uint64_t big[4] = {R[0]-1,R[1],R[2],R[3]}, bm2[4], r2m[4], x[4] = {0,0,0,0}, y[4];
    mul4(p, big, R2); redc(p, bm2); mul4(p, R2, R2); redc(p, r2m);
    y[0]=bm2[0]; y[1]=bm2[1]; y[2]=bm2[2]; y[3]=bm2[3]; acc4(x, y, 0); acc4(x, y, 0);   /* (r-1)+(r-1) in mont form */
    mul4(p, x, one); redc(p, t);
    printf("(r-1)+(r-1) mod r = %llu (expect r-2 = %llu)\n", (unsigned long long)t[0], (unsigned long long)(R[0]-2));
    uint64_t z0[4] = {0,0,0,0}, oo[4] = {1,0,0,0};
    acc4(z0, oo, ~0ULL);
    printf("0 - 1 = %llu %llu %llu %llu (expect r-1 = %llu ...)\n", (unsigned long long)z0[0], (unsigned long long)z0[1],
      (unsigned long long)z0[2], (unsigned long long)z0[3], (unsigned long long)(R[0]-1));
    return 0; }
  for (size_t i = 0; i < N; i++) {
    { uint64_t vv[4] = {i * 2654435761u, i * 7 + 1, i * 13 + 2, i * 29 + 3}, pp[8];
      mul4(pp, vv, R2); redc(pp, V[i]); }
    /* random 128-bit d, converted to Montgomery form as d*R2*2^-256 = d*2^256 mod r */
    { uint64_t dd[4] = {i * 40503u + 7, i * 2246822519u + 3, 0, 0}, pp[8];
      mul4(pp, dd, R2); redc(pp, D[i]); }
    for (int k = 0; k < 3; k++) ROW[i][k] = (uint32_t)((i * 2654435761u + k * 97) & (M - 1));
    SG[i] = (int8_t)((i & 1) ? 1 : -1);
  }
  for (int r = 0; r < reps; r++)
    for (size_t j = 0; j < N; j++) {
      uint64_t p[8], o[4]; mul4(p, V[j], D[j]); redc(p, o);
      uint64_t mask = (uint64_t)((int64_t)SG[j] >> 63);
      for (int k = 0; k < 3; k++) acc4(ACC[ROW[j][k]], o, mask);
    }
  printf("%llu\n", (unsigned long long)(ACC[0][0] + ACC[5][1]));
  return 0;
}
