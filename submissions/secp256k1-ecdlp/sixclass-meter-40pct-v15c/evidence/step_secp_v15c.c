/* step_secp_v15c.c — asm 6-mul affine walk step at secp256k1 prices.
 *
 * Levers vs v15 (1046): lazy Solinas mul (no final conditional subtract),
 * fe_sub (masked subtract of C, 10 Ir instead of 27), canonicalisation
 * only where the construction needs the canonical representative (table
 * index, DP test, parity canon of the stored point).
 */
#include <stdint.h>
typedef unsigned __int128 u128;
#include "secp_fe_v15c.h"

#ifndef V15C_TABLES_DEFINED
typedef struct { uint64_t x[4], y[4]; } tblent;
typedef struct { uint64_t x[4], y[4]; uint64_t len, esc; } awalker;
extern tblent TBL[2048];
extern awalker A[16];
extern uint64_t PROD[16][4], SUF[16][4];
#endif

__attribute__((noinline)) void step_asm(int n)
{
    for (int i = 0; i < n; i++) {
        u64 x1[4], y1[4];
        for (int k = 0; k < 4; k++) { x1[k] = A[i].x[k]; y1[k] = A[i].y[k]; }
        u64 j = x1[0] & 2047;                 /* stored point is canonical (< p) */
        u64 d[4], num[4], inv[4], lam[4], x3[4], y3[4], t[4];
        fe_sub(d, TBL[j].x, x1);
        fe_sub(num, TBL[j].y, y1);
        fe_mul(PROD[i], PROD[i], d);      /* 1 prefix   */
        fe_mul(SUF[i], SUF[i], d);        /* 2 suffix   */
        fe_mul(inv, PROD[i], SUF[i]);     /* 3 d^-1     */
        fe_mul(lam, num, inv);            /* 4 lambda   */
        fe_mul(t, lam, lam);              /* 5 lambda^2 */
        fe_sub(x3, t, x1); fe_sub(x3, x3, TBL[j].x);
        fe_sub(t, x1, x3);
        fe_mul(y3, lam, t);               /* 6          */
        fe_sub(y3, y3, y1);
        /* canon_pick: representative with even y (negation map) */
        u64 neg[4], e;
        neg[0] = __builtin_subcll(P0, y3[0], 0, &e);
        neg[1] = __builtin_subcll(P1, y3[1], e, &e);
        neg[2] = __builtin_subcll(P1, y3[2], e, &e);
        neg[3] = __builtin_subcll(P3, y3[3], e, &e);
        u64 odd = y3[0] & 1;                 /* y3 canonical: mul output minus canonical */
        u64 m = (u64)0 - odd;
        for (int k = 0; k < 4; k++) y3[k] = (y3[k] & ~m) | (neg[k] & m);
        /* distinguished-point test: 32 low bits of canonical x3 */
        u64 dp = x3[0] & 0xffffffffULL;      /* x3 canonical (mul output) */
        u64 len = A[i].len + 1, esc = A[i].esc | (dp == 0);
        for (int k = 0; k < 4; k++) { A[i].x[k] = x3[k]; A[i].y[k] = y3[k]; }
        A[i].len = len; A[i].esc = esc;
    }
}
