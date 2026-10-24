/* Counted routine for the label renormalisation multiplier of the g1-invert-g2
   class walk:  mulmod_r(a, c) = a * c mod r,  r = BLS12-381 scalar order.

   The multiplier c = p^k mod r is one of exactly TWELVE fixed constants
   (k = kappa in 0..11), so the package pre-scales them once at table-construction
   time:  cs[k] = p^k * 2^256 mod r.  The per-step work is then a single 4-limb
   CIOS Montgomery multiplication REDC(a * cs[k]) = a * p^k mod r.  No branches,
   no division, no conversion of a required.

   Called twice per step (labels a and b), so the per-step ledger row is
   2 x (counted instructions of this function). */
#include <stdint.h>

#define N0INV 0xfffffffeffffffffULL   /* -r^{-1} mod 2^64 */
#define R0 0xffffffff00000001ULL
#define R1 0x53bda402fffe5bfeULL
#define R2 0x3339d80809a1d805ULL
#define R3 0x73eda753299d7d48ULL

void mulmod_r(uint64_t out[4], const uint64_t a[4], const uint64_t cs[4])
{
    uint64_t t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0;
    for (int i = 0; i < 4; i++) {
        uint64_t ai = a[i], m, c;
        unsigned __int128 p;
        /* t += ai * cs */
        p = (unsigned __int128)ai * cs[0] + t0;        t0 = (uint64_t)p; c = (uint64_t)(p >> 64);
        p = (unsigned __int128)ai * cs[1] + t1 + c;    t1 = (uint64_t)p; c = (uint64_t)(p >> 64);
        p = (unsigned __int128)ai * cs[2] + t2 + c;    t2 = (uint64_t)p; c = (uint64_t)(p >> 64);
        p = (unsigned __int128)ai * cs[3] + t3 + c;    t3 = (uint64_t)p; c = (uint64_t)(p >> 64);
        p = (unsigned __int128)t4 + c;                 t4 = (uint64_t)p; t5 = (uint64_t)(p >> 64);
        /* m = t0 * n0inv mod 2^64 ; t = (t + m*r) / 2^64 */
        m = t0 * N0INV;
        p = (unsigned __int128)m * R0 + t0;            c = (uint64_t)(p >> 64);
        p = (unsigned __int128)m * R1 + t1 + c;        t0 = (uint64_t)p; c = (uint64_t)(p >> 64);
        p = (unsigned __int128)m * R2 + t2 + c;        t1 = (uint64_t)p; c = (uint64_t)(p >> 64);
        p = (unsigned __int128)m * R3 + t3 + c;        t2 = (uint64_t)p; c = (uint64_t)(p >> 64);
        p = (unsigned __int128)t4 + c;                 t3 = (uint64_t)p; c = (uint64_t)(p >> 64);
        t4 = t5 + c;
        t5 = 0;
    }
    /* branchless final conditional subtract of r */
    uint64_t d0, d1, d2, d3, borrow = 0;
    unsigned __int128 d;
    d = (unsigned __int128)t0 - R0 - borrow; d0 = (uint64_t)d; borrow = (uint64_t)((d >> 64) & 1);
    d = (unsigned __int128)t1 - R1 - borrow; d1 = (uint64_t)d; borrow = (uint64_t)((d >> 64) & 1);
    d = (unsigned __int128)t2 - R2 - borrow; d2 = (uint64_t)d; borrow = (uint64_t)((d >> 64) & 1);
    d = (unsigned __int128)t3 - R3 - borrow; d3 = (uint64_t)d; borrow = (uint64_t)((d >> 64) & 1);
    uint64_t mask = borrow - 1;    /* all ones iff no borrow, i.e. t >= r */
    out[0] = (d0 & mask) | (t0 & ~mask);
    out[1] = (d1 & mask) | (t1 & ~mask);
    out[2] = (d2 & mask) | (t2 & ~mask);
    out[3] = (d3 & mask) | (t3 & ~mask);
    (void)t4; (void)t5;
}
