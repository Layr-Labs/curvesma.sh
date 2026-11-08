/* X = U*D projection: acc += sign * (v * d) with v a 4-limb, d a 2-limb (128-bit) diagonal entry,
   accumulator 6-limb, full carries.  U: 3 rows per column, uniform signs, branch-free two's
   complement.  Per column element this is the real arithmetic: 4x2 limb products with carry
   propagation, then a 6-limb signed accumulate. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define M 28
static uint64_t V[1<<20][4], ACC[1<<12][6], D[1<<20][2]; static uint32_t ROW[1<<20][3]; static int8_t SG[1<<20];
static inline void mul48(const uint64_t v[4], const uint64_t d[2], uint64_t p[6]) {
  unsigned __int128 c = 0;
  for (int i = 0; i < 6; i++) {
    unsigned __int128 s = c;
    for (int j = (i < 4 ? 0 : i - 3); j <= i && j < 4; j++) s += (unsigned __int128)v[j] * d[i - j];
    p[i] = (uint64_t)s; c = s >> 64;
  }
}
static inline void acc6(uint64_t a[6], const uint64_t p[6], uint64_t mask) {
  unsigned char c = (unsigned char)(mask & 1);
  for (int i = 0; i < 6; i++) { unsigned __int128 zz = (unsigned __int128)a[i] + (p[i] ^ mask) + c; a[i] = (uint64_t)zz; c = (unsigned char)(zz >> 64); }
}
int main(int argc, char **argv) {
  int reps = (argc > 1 && argv[1][0] >= '0' && argv[1][0] <= '9') ? atoi(argv[1]) : (argc > 2 ? atoi(argv[2]) : 0);
  size_t N = 1 << 20;
  if (argc > 1 && !strcmp(argv[1], "selftest")) {   /* v=1, d=1 -> 1; v=0 -> 0 */
    uint64_t v[4] = {1,0,0,0}, d[2] = {1,0}, p[6], a[6] = {0,0,0,0,0,0};
    mul48(v, d, p); acc6(a, p, 0); printf("1*1 -> %llu %llu\n", (unsigned long long)a[0], (unsigned long long)a[1]);
    uint64_t z[6] = {0,0,0,0,0,0}, v0[4] = {0,0,0,0}; mul48(v0, d, p);
    printf("0*1 -> %d  (p all zero: %d)\n", (int)((p[0]|p[1]|p[2]|p[3]|p[4]|p[5])==0), (int)((z[0]|z[1])==0));
    return 0; }
  for (size_t i = 0; i < N; i++) { for (int q = 0; q < 4; q++) V[i][q] = i * 2654435761u + q * 7;
    D[i][0] = i * 40503u + 7; D[i][1] = i * 2246822519u + 3;
    for (int k = 0; k < 3; k++) ROW[i][k] = (uint32_t)((i * 2654435761u + k * 97) & (M - 1));
    SG[i] = (int8_t)((i & 1) ? 1 : -1); }
  for (int r = 0; r < reps; r++)
    for (size_t j = 0; j < N; j++) {
      uint64_t p[6]; mul48(V[j], D[j], p);
      uint64_t mask = (uint64_t)((int64_t)SG[j] >> 63);
      for (int k = 0; k < 3; k++) acc6(ACC[ROW[j][k]], p, mask); }
  printf("%llu\n", (unsigned long long)(ACC[0][0] + ACC[5][1]));
  return 0;
}
