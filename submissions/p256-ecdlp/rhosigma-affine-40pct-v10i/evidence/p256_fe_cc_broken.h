/* MUTATION CONTROL: the v10f header VERBATIM (P1 = 0, P2 = 0x00000000ffffffff,
 * limbs 1 and 2 of p transposed).  Only an include guard was added so this file
 * can be included by the same driver TU; no arithmetic line was touched.
 * The kernel test MUST fail against this header or it certifies nothing.
 */
/* p256_fe_cc.h — P-256 field kernels as carry chains (v9f).
 *
 * Every carry is an adc/sbb chain (__builtin_addcll / __builtin_subcll), the
 * 16 limb products of the Montgomery multiply are accumulated in CIOS order,
 * and each reduction round is ONE carry chain over the four limbs it touches.
 *
 * Two structural facts about p are used, and nothing else:
 *   p0 = 2^64 - 1  =>  -p^-1 = 1 (NP0 = 1), so a round's quotient digit is the
 *                      low limb itself, and  q*p0 + t0 = q*2^64 exactly: the
 *                      round's low limb is zero and its carry-in is q.
 *   p2 = 0         =>  the third product of a round is q*0 = 0, i.e. one add,
 *                      not a mul/umulh pair.
 * Contract: canonical inputs a,b < p (every producer in the walk is
 * fe_add/fe_sub/fe_mul on canonical operands, and each returns canonical).
 * Invariant used by the reduction: the accumulator stays < 2p, so its top limb
 * is 0 or 1 and no sixth limb exists.
 */
#ifndef P256_FE_CC_H
#define P256_FE_CC_H

#include <stdint.h>

typedef uint64_t u64;
typedef unsigned __int128 u128;

#define P0 0xffffffffffffffffULL
#define P1 0x0000000000000000ULL
#define P2 0x00000000ffffffffULL
#define P3 0xffffffff00000001ULL

/* r = a + b mod p : one adc chain, one sbb chain, 4 csel */
static __attribute__((always_inline)) inline void fe_add(u64 r[4], const u64 a[4], const u64 b[4])
{
    u64 c, d;
    u64 t0 = __builtin_addcll(a[0], b[0], 0, &c);
    u64 t1 = __builtin_addcll(a[1], b[1], c, &c);
    u64 t2 = __builtin_addcll(a[2], b[2], c, &c);
    u64 t3 = __builtin_addcll(a[3], b[3], c, &c);
    u64 cy = c;                                  /* carry out of limb 3 */
    u64 s0 = __builtin_subcll(t0, P0, 0, &d);
    u64 s1 = __builtin_subcll(t1, P1, d, &d);
    u64 s2 = __builtin_subcll(t2, P2, d, &d);
    u64 s3 = __builtin_subcll(t3, P3, d, &d);
    u64 take = cy | (d ^ 1);      /* t >= p  <=>  no borrow out, or carry out */
    r[0] = take ? s0 : t0;
    r[1] = take ? s1 : t1;
    r[2] = take ? s2 : t2;
    r[3] = take ? s3 : t3;
}

/* r = a - b mod p : one sbb chain, one adc chain, 4 csel on the borrow */
static __attribute__((always_inline)) inline void fe_sub(u64 r[4], const u64 a[4], const u64 b[4])
{
    u64 d;
    u64 t0 = __builtin_subcll(a[0], b[0], 0, &d);
    u64 t1 = __builtin_subcll(a[1], b[1], d, &d);
    u64 t2 = __builtin_subcll(a[2], b[2], d, &d);
    u64 t3 = __builtin_subcll(a[3], b[3], d, &d);
    u64 m = 0 - d;                               /* borrow => add p back */
    u64 c;
    u64 s0 = __builtin_addcll(t0, P0 & m, 0, &c);
    u64 s1 = __builtin_addcll(t1, P1 & m, c, &c);
    u64 s2 = __builtin_addcll(t2, P2, c, &c);
    u64 s3 = __builtin_addcll(t3, P3 & m, c, &c);
    r[0] = d ? s0 : t0;
    r[1] = d ? s1 : t1;
    r[2] = d ? s2 : t2;
    r[3] = d ? s3 : t3;
}

/* ---- Montgomery multiply, CIOS with 5-limb accumulator ---- */

typedef struct { u64 t0, t1, t2, t3, t4; } acc5;

/* T += a * bi : four mul/umulh pairs, every carry in an adc chain */
__attribute__((always_inline)) static inline void cios_row(acc5 *T, const u64 a[4], u64 bi)
{
    u64 cc, hi, lo;
    u128 m = (u128)a[0] * bi;
    lo = (u64)m; hi = (u64)(m >> 64);
    u64 n0 = __builtin_addcll(lo, T->t0, 0, &cc);
    hi = __builtin_addcll(hi, 0, cc, &cc);

    m = (u128)a[1] * bi + hi;
    lo = (u64)m; hi = (u64)(m >> 64);
    u64 n1 = __builtin_addcll(lo, T->t1, 0, &cc);
    hi = __builtin_addcll(hi, 0, cc, &cc);

    m = (u128)a[2] * bi + hi;
    lo = (u64)m; hi = (u64)(m >> 64);
    u64 n2 = __builtin_addcll(lo, T->t2, 0, &cc);
    hi = __builtin_addcll(hi, 0, cc, &cc);

    m = (u128)a[3] * bi + hi;
    lo = (u64)m; hi = (u64)(m >> 64);
    u64 n3 = __builtin_addcll(lo, T->t3, 0, &cc);
    hi = __builtin_addcll(hi, 0, cc, &cc);

    T->t0 = n0; T->t1 = n1; T->t2 = n2; T->t3 = n3;
    T->t4 = __builtin_addcll(T->t4, hi, 0, &cc);
}

/* T <- (T + q*p) / 2^64 with q = T->t0 : one adc chain over three products */
__attribute__((always_inline)) static inline void cios_red(acc5 *T)
{
    u64 q = T->t0, cc, cc2;
    u128 m = (u128)q * P1 + T->t1 + q;            /* carry-in is q, not 0 */
    u64 n0 = (u64)m, hi = (u64)(m >> 64);
    u64 n1 = __builtin_addcll(T->t2, hi, 0, &cc); /* p2 = 0: add, no product */
    m = (u128)q * P3 + T->t3 + cc;
    u64 n2 = (u64)m;
    u64 n3 = __builtin_addcll(T->t4, (u64)(m >> 64), 0, &cc2);
    T->t0 = n0; T->t1 = n1; T->t2 = n2; T->t3 = n3; T->t4 = cc2;
}

__attribute__((always_inline)) static inline void fe_mul(u64 r[4], const u64 a[4], const u64 b[4])
{
    acc5 T = {0, 0, 0, 0, 0};
    cios_row(&T, a, b[0]); cios_red(&T);
    cios_row(&T, a, b[1]); cios_red(&T);
    cios_row(&T, a, b[2]); cios_red(&T);
    cios_row(&T, a, b[3]); cios_red(&T);
    u64 d;
    u64 s0 = __builtin_subcll(T.t0, P0, 0, &d);
    u64 s1 = __builtin_subcll(T.t1, P1, d, &d);
    u64 s2 = __builtin_subcll(T.t2, 0, d, &d);
    u64 s3 = __builtin_subcll(T.t3, P3, d, &d);
    u64 take = T.t4 | (d ^ 1);
    r[0] = take ? s0 : T.t0;
    r[1] = take ? s1 : T.t1;
    r[2] = take ? s2 : T.t2;
    r[3] = take ? s3 : T.t3;
}

static void fe_inv(u64 r[4], const u64 a[4])
{
    static const u64 E[4] = {0xfffffffffffffffdULL, 0x00000000ffffffffULL,
                             0, 0xffffffff00000001ULL};
    u64 acc[4] = {1, 0, 0, 0}, base[4];
    for (int i = 0; i < 4; i++) base[i] = a[i];
    for (int i = 255; i >= 0; i--) {
        fe_mul(acc, acc, acc);
        if ((E[i >> 6] >> (i & 63)) & 1) fe_mul(acc, acc, base);
    }
    for (int i = 0; i < 4; i++) r[i] = acc[i];
}

#endif /* P256_FE_CC_H */
