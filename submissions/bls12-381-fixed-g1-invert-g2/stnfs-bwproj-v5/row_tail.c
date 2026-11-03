typedef unsigned long uint64_t;
typedef unsigned __int128 u128;
/* Per-row tail of one product (STORAGE.md). A (pos, neg) pair is two six-limb sums.
   A (block, sign) group holds distinct columns of one 2^56-element block, so its sum is
   < 2^56 x 255 x r < 2^319; a row has < 2^60 entries, so every combined sum is < 2^323. */
uint64_t *send_pair(uint64_t *__restrict out, const uint64_t *__restrict loc) {
  for (int i = 0; i < 12; i++) out[i] = loc[i];
  return out + 12;
}
const uint64_t *combine(uint64_t *loc, const uint64_t *in) {
  for (int s = 0; s < 2; s++) {
    u128 c = 0;
    for (int i = 0; i < 6; i++) { c += (u128)loc[6*s+i] + in[6*s+i]; loc[6*s+i] = (uint64_t)c; c >>= 64; }
  }
  return in + 12;
}
/* reduce6: v < 2^384 (six limbs) -> v mod r (four limbs), r < 2^255.
   hi = floor(v / 2^192) < 2^192, MU = floor(2^448 / r) < 2^194 (four limbs),
   qh = floor(hi MU / 2^256) is q or q - 1 (q = floor(v / r) < 2^130), so v - qh r < 2r. */
__attribute__((noinline)) void reduce6(uint64_t *y, uint64_t *v, const uint64_t *MU, const uint64_t *r) {
  uint64_t hi[3] = {v[3], v[4], v[5]}, p[7] = {0};
  for (int i = 0; i < 3; i++) {
    u128 c = 0;
    for (int j = 0; j < 4; j++) { c += (u128)hi[i] * MU[j] + p[i+j]; p[i+j] = (uint64_t)c; c >>= 64; }
    p[i+4] = (uint64_t)c;
  }
  uint64_t q[3] = {p[4], p[5], p[6]};          /* qh < 2^131 */
  uint64_t t[6] = {0};
  for (int i = 0; i < 3; i++) {
    u128 c = 0;
    for (int j = 0; j < 4; j++) { c += (u128)q[i] * r[j] + t[i+j]; t[i+j] = (uint64_t)c; c >>= 64; }
    if (i + 4 < 6) t[i+4] = (uint64_t)c;
  }
  uint64_t b = 0;
  for (int i = 0; i < 6; i++) { u128 s = (u128)v[i] - t[i] - b; v[i] = (uint64_t)s; b = (uint64_t)(s >> 64) & 1; }
  /* now v < 2r < 2^256: v[4] = v[5] = 0 */
  uint64_t w[4]; b = 0;
  for (int i = 0; i < 4; i++) { u128 s = (u128)v[i] - r[i] - b; w[i] = (uint64_t)s; b = (uint64_t)(s >> 64) & 1; }
  uint64_t keep = (uint64_t)0 - b;               /* all ones when v < r */
  for (int i = 0; i < 4; i++) y[i] = (v[i] & keep) | (w[i] & ~keep);
}
/* fold9: x < 2^576 (nine limbs, the dense bw_full accumulator) -> six limbs, value
   x0..3 + sum_{k=4..8} x_k R64[k] < 2^256 + 5 x 2^64 x 2^255 < 2^323, with R64[k] = 2^(64k) mod r. */
__attribute__((noinline)) void fold9(uint64_t *a, const uint64_t *x, const uint64_t (*R64)[4]) {
  a[0] = x[0]; a[1] = x[1]; a[2] = x[2]; a[3] = x[3]; a[4] = 0; a[5] = 0;
  for (int k = 4; k < 9; k++) {
    u128 c = 0; uint64_t t = x[k];
    for (int i = 0; i < 4; i++) { c += (u128)R64[k][i] * t + a[i]; a[i] = (uint64_t)c; c >>= 64; }
    c += a[4]; a[4] = (uint64_t)c; a[5] += (uint64_t)(c >> 64);
  }
}
/* finish at the row owner: y = (pos - neg + dense) mod r, with K = k r in [2^323, 2^324):
   v = pos + (K - neg) + fold9(dense) < 2^323 + 2^324 + 2^323 = 2^325. */
void finish(uint64_t *y, const uint64_t *pair, const uint64_t *dense, const uint64_t *K,
            const uint64_t (*R64)[4], const uint64_t *MU, const uint64_t *r) {
  uint64_t v[6], a[6]; u128 c = 0; uint64_t b = 0;
  fold9(a, dense, R64);
  for (int i = 0; i < 6; i++) {
    u128 t = (u128)K[i] - pair[6+i] - b; b = (uint64_t)(t >> 64) & 1;
    c += (u128)pair[i] + (uint64_t)t + a[i]; v[i] = (uint64_t)c; c >>= 64;
  }
  reduce6(y, v, MU, r);
}
/* reduce9: a dense accumulator alone (projection entries, Horner elements) -> mod r. */
void reduce9(uint64_t *y, const uint64_t *dense, const uint64_t (*R64)[4], const uint64_t *MU, const uint64_t *r) {
  uint64_t a[6];
  fold9(a, dense, R64);
  reduce6(y, a, MU, r);
}
