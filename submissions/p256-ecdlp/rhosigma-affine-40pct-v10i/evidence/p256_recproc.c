/* recproc.c — PRODUCTION record processor: DP-triggered sparse record build,
 * hash insert, and one Montgomery batch-inversion round over n=2^16 walkers. */
#include <stdint.h>
#include "p256_fe_cc.h"
#define NHASH (1u<<18)
typedef struct { uint64_t key; uint64_t a[4], b[4]; } srec;   /* 8+32+32 = 72 B */
extern srec HT[NHASH];
extern uint64_t QREC[4][6];
extern uint32_t QN, QCAP;
static inline uint64_t mix64(uint64_t v)
{ v *= 0x9e3779b97f4a7c15ULL; v ^= v >> 29; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 32; return v; }
/* one report: detect was done in the step; here = build record + insert */
void recproc(uint32_t i)
{
    uint64_t h = mix64(QREC[i][0]) & (NHASH - 1);
    srec *p = &HT[h];
    if (p->key == 0) {
        p->key = QREC[i][0] + 1;
        for (int k = 0; k < 4; k++) { p->a[k] = QREC[i][1 + k]; p->b[k] = QREC[i][5 - 0 + k - 4 + 4]; }
    }
}
/* Montgomery batch inversion of n denominators: 3(n-1) muls + 1 inv per n steps */
void batch_inv(uint64_t d[][4], uint64_t out[][4], int n)
{
    uint64_t pre[64][4], inv[4];
    pre[0][0] = 1; pre[0][1] = pre[0][2] = pre[0][3] = 0;
    for (int i = 0; i < n - 1; i++) fe_mul(pre[i + 1], pre[i], d[i]);
    fe_inv(inv, pre[n - 1]);
    for (int i = n - 1; i > 0; i--) {
        fe_mul(out[i], inv, pre[i - 1]);
        fe_mul(inv, inv, d[i]);
    }
    fe_mul(out[0], inv, pre[0]);
}
