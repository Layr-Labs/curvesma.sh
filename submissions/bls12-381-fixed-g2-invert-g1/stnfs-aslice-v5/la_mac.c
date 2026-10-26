#include <stdint.h>
/* Block-Wiedemann mat-vec inner loop over F_r (r = BLS12-381 255-bit, 4 limbs).
 * Matrix entries after filtering are small nonzero integers, |c| <= 255.  Each row is
 * stored sign-split: its positive entries, then its negative entries (magnitudes as
 * uint8), with the split point npos.  Both halves run the SAME unsigned
 * multiply-by-small + lazy 5-limb accumulate loop, into two accumulators; the row
 * result is (acc_pos - acc_neg) mod r, one reduction per row (row_finish).
 * x[] holds fully reduced representatives in [0, r). */
typedef struct { uint64_t v[4]; } fe_t;
typedef struct { uint64_t v[5]; } acc_t;

void mac_run(acc_t *restrict acc, const fe_t *restrict x,
             const uint32_t *restrict idx, const uint8_t *restrict val, long w) {
    uint64_t a0 = acc->v[0], a1 = acc->v[1], a2 = acc->v[2], a3 = acc->v[3], a4 = acc->v[4];
    for (long j = 0; j < w; j++) {
        const uint64_t *xe = x[idx[j]].v;
        uint64_t k = val[j];
        __uint128_t p0 = (__uint128_t)xe[0] * k, p1 = (__uint128_t)xe[1] * k;
        __uint128_t p2 = (__uint128_t)xe[2] * k, p3 = (__uint128_t)xe[3] * k;
        /* acc += p0 + p1<<64 + p2<<128 + p3<<192 */
        __uint128_t s = (__uint128_t)a0 + (uint64_t)p0;
        a0 = (uint64_t)s;
        s = (s >> 64) + a1 + (uint64_t)(p0 >> 64) + (uint64_t)p1;
        a1 = (uint64_t)s;
        s = (s >> 64) + a2 + (uint64_t)(p1 >> 64) + (uint64_t)p2;
        a2 = (uint64_t)s;
        s = (s >> 64) + a3 + (uint64_t)(p2 >> 64) + (uint64_t)p3;
        a3 = (uint64_t)s;
        a4 += (uint64_t)(s >> 64) + (uint64_t)(p3 >> 64);
    }
    acc->v[0] = a0; acc->v[1] = a1; acc->v[2] = a2; acc->v[3] = a3; acc->v[4] = a4;
}

void row_dot(acc_t *pos, acc_t *neg, const fe_t *x, const uint32_t *idx,
             const uint8_t *val, long npos, long w) {
    mac_run(pos, x, idx, val, npos);
    mac_run(neg, x, idx + npos, val + npos, w - npos);
}
