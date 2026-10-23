/* secp_fe_v15c.h — HAND-SCHEDULED aarch64 secp256k1 field kernels (v15c).
 *
 * p = 2^256 - C, C = 2^32 + 977.  2^256 == C (mod p).
 *
 * fe_mul: 4x4 schoolbook, 8-limb accumulator resident in GPRs, then the
 * Solinas fold as a *multiply* by C (not a shift/add pair):
 *     L + H*2^256 == L + H*C (mod p),  C < 2^33  =>  H*C is one 4x4-by-constant
 *     mul/umulh chain, 5 limbs.
 * Two further folds (w4*C, then carry*C) close the reduction in 4 limbs, then
 * ONE branchless conditional subtract of p written as "add C, take on carry"
 * (exact: 2^256 == C, so y + C mod 2^256 == y - p when y >= p; y < 2^256 is
 * invariant). Output is CANONICAL, r < p.
 *
 * fe_sub: 4 sbb + csel(C on borrow) + 4 sbb = 9 Ir, exact for a,b < 2^256
 *         whenever |a-b| <= p (v15b gate: the 4-limb lazy values keep this).
 *
 * fe_canon_low32: canonical low 32 bits without materialising the value.
 */
#ifndef SECP_FE_V15C_H
#define SECP_FE_V15C_H
#include <stdint.h>

typedef uint64_t u64;

#define P0 0xFFFFFFFEFFFFFC2FULL
#define P1 0xFFFFFFFFFFFFFFFFULL
#define P3 0xFFFFFFFFFFFFFFFFULL
#define CSOL 0x1000003D1ULL

/* r = a*b mod p, canonical (r < p). val */
__attribute__((always_inline)) static inline void fe_mul(u64 r[4], const u64 a[4], const u64 b[4])
{
    static const u64 KC = CSOL;

    u64 o0, o1, o2, o3, q0, q1, q2, q3, s0, s1, s2, s3, s4, s5, s6;
    u64 a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
    u64 b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
    (void)q0; (void)q1; (void)q2; (void)q3;
    (void)s0; (void)s1; (void)s2; (void)s3; (void)s4; (void)s5; (void)s6;
    asm volatile(
    "	/* a,b values in registers */\n\t"
    "	mul	%[o0], %[a0], %[b0]\n\t"
    "	umulh	%[s1], %[a0], %[b0]\n\t"
    "	mul	%[s5], %[a1], %[b0]\n\t"
    "	umulh	%[s6], %[a1], %[b0]\n\t"
    "	adds	%[o1], %[s5], %[s1]\n\t"
    "	adc	%[s2], %[s6], xzr\n\t"
    "	mul	%[s5], %[a2], %[b0]\n\t"
    "	umulh	%[s6], %[a2], %[b0]\n\t"
    "	adds	%[o2], %[s5], %[s2]\n\t"
    "	adc	%[s3], %[s6], xzr\n\t"
    "	mul	%[s5], %[a3], %[b0]\n\t"
    "	umulh	%[s6], %[a3], %[b0]\n\t"
    "	adds	%[o3], %[s5], %[s3]\n\t"
    "	adc	%[q0], %[s6], xzr\n\t"
    "	mul	%[s0], %[a0], %[b1]\n\t"
    "	umulh	%[s1], %[a0], %[b1]\n\t"
    "	mul	%[s5], %[a1], %[b1]\n\t"
    "	umulh	%[s6], %[a1], %[b1]\n\t"
    "	adds	%[s1], %[s1], %[s5]\n\t"
    "	adc	%[s2], %[s6], xzr\n\t"
    "	mul	%[s5], %[a2], %[b1]\n\t"
    "	umulh	%[s6], %[a2], %[b1]\n\t"
    "	adds	%[s2], %[s2], %[s5]\n\t"
    "	adc	%[s3], %[s6], xzr\n\t"
    "	mul	%[s5], %[a3], %[b1]\n\t"
    "	umulh	%[s6], %[a3], %[b1]\n\t"
    "	adds	%[s3], %[s3], %[s5]\n\t"
    "	adc	%[s4], %[s6], xzr\n\t"
    "	adds	%[o1], %[o1], %[s0]\n\t"
    "	adcs	%[o2], %[o2], %[s1]\n\t"
    "	adcs	%[o3], %[o3], %[s2]\n\t"
    "	adcs	%[q0], %[q0], %[s3]\n\t"
    "	adc	%[q1], %[s4], xzr\n\t"
    "	mul	%[s0], %[a0], %[b2]\n\t"
    "	umulh	%[s1], %[a0], %[b2]\n\t"
    "	mul	%[s5], %[a1], %[b2]\n\t"
    "	umulh	%[s6], %[a1], %[b2]\n\t"
    "	adds	%[s1], %[s1], %[s5]\n\t"
    "	adc	%[s2], %[s6], xzr\n\t"
    "	mul	%[s5], %[a2], %[b2]\n\t"
    "	umulh	%[s6], %[a2], %[b2]\n\t"
    "	adds	%[s2], %[s2], %[s5]\n\t"
    "	adc	%[s3], %[s6], xzr\n\t"
    "	mul	%[s5], %[a3], %[b2]\n\t"
    "	umulh	%[s6], %[a3], %[b2]\n\t"
    "	adds	%[s3], %[s3], %[s5]\n\t"
    "	adc	%[s4], %[s6], xzr\n\t"
    "	adds	%[o2], %[o2], %[s0]\n\t"
    "	adcs	%[o3], %[o3], %[s1]\n\t"
    "	adcs	%[q0], %[q0], %[s2]\n\t"
    "	adcs	%[q1], %[q1], %[s3]\n\t"
    "	adc	%[q2], %[s4], xzr\n\t"
    "	mul	%[s0], %[a0], %[b3]\n\t"
    "	umulh	%[s1], %[a0], %[b3]\n\t"
    "	mul	%[s5], %[a1], %[b3]\n\t"
    "	umulh	%[s6], %[a1], %[b3]\n\t"
    "	adds	%[s1], %[s1], %[s5]\n\t"
    "	adc	%[s2], %[s6], xzr\n\t"
    "	mul	%[s5], %[a2], %[b3]\n\t"
    "	umulh	%[s6], %[a2], %[b3]\n\t"
    "	adds	%[s2], %[s2], %[s5]\n\t"
    "	adc	%[s3], %[s6], xzr\n\t"
    "	mul	%[s5], %[a3], %[b3]\n\t"
    "	umulh	%[s6], %[a3], %[b3]\n\t"
    "	adds	%[s3], %[s3], %[s5]\n\t"
    "	adc	%[s4], %[s6], xzr\n\t"
    "	adds	%[o3], %[o3], %[s0]\n\t"
    "	adcs	%[q0], %[q0], %[s1]\n\t"
    "	adcs	%[q1], %[q1], %[s2]\n\t"
    "	adcs	%[q2], %[q2], %[s3]\n\t"
    "	adc	%[q3], %[s4], xzr\n\t"
    "	mul	%[s0], %[q0], %[C]\n\t"
    "	umulh	%[s1], %[q0], %[C]\n\t"
    "	mul	%[s5], %[q1], %[C]\n\t"
    "	umulh	%[s6], %[q1], %[C]\n\t"
    "	adds	%[s1], %[s1], %[s5]\n\t"
    "	adc	%[s2], %[s6], xzr\n\t"
    "	mul	%[s5], %[q2], %[C]\n\t"
    "	umulh	%[s6], %[q2], %[C]\n\t"
    "	adds	%[s2], %[s2], %[s5]\n\t"
    "	adc	%[s3], %[s6], xzr\n\t"
    "	mul	%[s5], %[q3], %[C]\n\t"
    "	umulh	%[s6], %[q3], %[C]\n\t"
    "	adds	%[s3], %[s3], %[s5]\n\t"
    "	adc	%[s4], %[s6], xzr\n\t"
    "	adds	%[o0], %[o0], %[s0]\n\t"
    "	adcs	%[o1], %[o1], %[s1]\n\t"
    "	adcs	%[o2], %[o2], %[s2]\n\t"
    "	adcs	%[o3], %[o3], %[s3]\n\t"
    "	adc	%[q0], %[s4], xzr\n\t"
    "	mul	%[s0], %[q0], %[C]\n\t"
    "	umulh	%[s1], %[q0], %[C]\n\t"
    "	adds	%[o0], %[o0], %[s0]\n\t"
    "	adcs	%[o1], %[o1], %[s1]\n\t"
    "	adcs	%[o2], %[o2], xzr\n\t"
    "	adcs	%[o3], %[o3], xzr\n\t"
    "	adc	%[s2], xzr, xzr\n\t"
    "	adds	%[s0], %[o0], %[C]\n\t"
    "	adcs	%[s1], %[o1], xzr\n\t"
    "	adcs	%[s3], %[o2], xzr\n\t"
    "	adcs	%[s4], %[o3], xzr\n\t"
    "	adc	%[s5], xzr, xzr\n\t"
    "	adds	%[s2], %[s2], %[s5]\n\t"
    "	csel	%[o0], %[s0], %[o0], ne\n\t"
    "	csel	%[o1], %[s1], %[o1], ne\n\t"
    "	csel	%[o2], %[s3], %[o2], ne\n\t"
    "	csel	%[o3], %[s4], %[o3], ne\n\t"
    : [o0] "=&r"(o0), [o1] "=&r"(o1), [o2] "=&r"(o2), [o3] "=&r"(o3),
      [q0] "=&r"(q0), [q1] "=&r"(q1), [q2] "=&r"(q2), [q3] "=&r"(q3),
      [s0] "=&r"(s0), [s1] "=&r"(s1), [s2] "=&r"(s2), [s3] "=&r"(s3),
      [s4] "=&r"(s4), [s5] "=&r"(s5), [s6] "=&r"(s6)
    : [a0] "r"(a0), [a1] "r"(a1), [a2] "r"(a2), [a3] "r"(a3),
      [b0] "r"(b0), [b1] "r"(b1), [b2] "r"(b2), [b3] "r"(b3),
      [C] "r"(KC)
    : "cc");
    r[0] = o0; r[1] = o1; r[2] = o2; r[3] = o3;

}

/* r = a-b mod p ; exact for 4-limb a,b with |a-b| <= p */
__attribute__((always_inline)) static inline void fe_sub(u64 r[4], const u64 a[4], const u64 b[4])
{
    static const u64 KC = CSOL;
    u64 a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
    u64 b0 = b[0], b1 = b[1], b2 = b[2], b3 = b[3];
    u64 o0, o1, o2, o3, m;
    asm volatile(
        "subs\t%[o0], %[a0], %[b0]\n\t"
        "sbcs\t%[o1], %[a1], %[b1]\n\t"
        "sbcs\t%[o2], %[a2], %[b2]\n\t"
        "sbcs\t%[o3], %[a3], %[b3]\n\t"
        "csel\t%[m], %[C], xzr, cc\n\t"
        "subs\t%[o0], %[o0], %[m]\n\t"
        "sbcs\t%[o1], %[o1], xzr\n\t"
        "sbcs\t%[o2], %[o2], xzr\n\t"
        "sbcs\t%[o3], %[o3], xzr"
        : [o0] "=&r"(o0), [o1] "=&r"(o1), [o2] "=&r"(o2), [o3] "=&r"(o3), [m] "=&r"(m)
        : [a0] "r"(a0), [a1] "r"(a1), [a2] "r"(a2), [a3] "r"(a3),
          [b0] "r"(b0), [b1] "r"(b1), [b2] "r"(b2), [b3] "r"(b3), [C] "r"(KC)
        : "cc");
    r[0] = o0; r[1] = o1; r[2] = o2; r[3] = o3;
}

/* canonicalise a lazy value: subtract p once if v >= p */
__attribute__((always_inline)) static inline void fe_canon(u64 r[4], const u64 v[4])
{
    u64 d;
    u64 q0 = __builtin_subcll(v[0], P0, 0, &d);
    u64 q1 = __builtin_subcll(v[1], P1, d, &d);
    u64 q2 = __builtin_subcll(v[2], P1, d, &d);
    u64 q3 = __builtin_subcll(v[3], P3, d, &d);
    u64 take = d ^ 1;
    r[0] = take ? q0 : v[0];
    r[1] = take ? q1 : v[1];
    r[2] = take ? q2 : v[2];
    r[3] = take ? q3 : v[3];
}

/* canonical low 32 bits of a lazy 4-limb value; v >= p iff v + C carries out
 * of 2^256, i.e. limb0 + C carries and limbs 1..3 are all ones. */
__attribute__((always_inline)) static inline u64 canon_low32(const u64 v[4])
{
    u64 s = v[0] + CSOL;
    u64 c0 = (s < CSOL);
    u64 all1 = v[1] & v[2] & v[3];
    u64 take = c0 & (u64)(all1 == ~0ULL);
    u64 lo = take ? s : v[0];
    return lo & 0xffffffffULL;
}
#endif
