/* cios4.c -- 4x64-bit Montgomery multiplication (CIOS), the kernel the ledger
 * prices.  Build, count and self-test:
 *
 *     clang -O2 -c cios4.c -o cios4.o            # static span (row_extract.py)
 *     clang -O2 -DSELFTEST cios4.c -o cios4_test  # correctness vs a Python reference
 *     python3 scripts/row_extract.py
 *
 * v17a made two mistakes a reviewer caught in verdict 40b720e2 point 3, both
 * fixed here:
 *   (1) the reduction chain accumulated m*n[j] into t[j] instead of into t[j-1]
 *       -- i.e. it never performed the CIOS shift-down, so the routine was not
 *       the multiplication it claimed to be.  The accumulator now has the extra
 *       limb, the shift is explicit, and the last carry is propagated.
 *   (2) the COUNT was the static span.  One call executes the outer loop four
 *       times, so the honest number is the EXECUTED count (row_extract.py reports
 *       both: static span and executed rows = once-rows + 4 x loop-body-rows).
 *
 * n0inv = -n[0]^{-1} mod 2^64 is computed by the caller; cios4_selftest below
 * computes it with Newton iteration so the test needs no helper.
 */
#include <stdint.h>
#include <stdio.h>

typedef struct { uint64_t l[4]; } u256;

/* t += a*b + c ; returns the new carry (mac over a 128-bit product). */
static inline uint64_t mac(uint64_t *t, uint64_t a, uint64_t b, uint64_t c)
{
    unsigned __int128 s = (unsigned __int128)a * b + *t + c;
    *t = (uint64_t)s;
    return (uint64_t)(s >> 64);
}

__attribute__((noinline, visibility("default")))
void cios4(u256 *out, const u256 *a, const u256 *b, const u256 *n, uint64_t n0inv)
{
    uint64_t t[5] = {0, 0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
        uint64_t c = 0;
        for (int j = 0; j < 4; j++)             /* t += a[i]*b */
            c = mac(&t[j], a->l[i], b->l[j], c);
        t[4] = c;

        uint64_t m = t[0] * n0inv;              /* cancels t[0] mod 2^64 */
        c = mac(&t[0], m, n->l[0], 0);          /* t[0] becomes 0            */
        for (int j = 1; j < 4; j++)
            c = mac(&t[j], m, n->l[j], c);
        t[4] += c;

        for (int j = 0; j < 4; j++)             /* shift the accumulator down */
            t[j] = t[j + 1];
        t[4] = 0;
    }
    uint64_t borrow = 0, r[4];
    for (int j = 0; j < 4; j++) {               /* final conditional subtract */
        unsigned __int128 d = (unsigned __int128)t[j] - n->l[j] - borrow;
        r[j] = (uint64_t)d;
        borrow = (uint64_t)(d >> 64) & 1;
    }
    uint64_t mask = ~(uint64_t)0 + borrow;      /* all ones when t >= n */
    for (int j = 0; j < 4; j++)
        out->l[j] = (r[j] & mask) | (t[j] & ~mask);
}

/* span fence: the extractor counts from _cios4 up to this label */
__attribute__((noinline, visibility("default")))
void cios4_end(void) { }

#ifdef SELFTEST
static uint64_t n0inv64(uint64_t n0)            /* -n0^{-1} mod 2^64 */
{
    uint64_t x = 1;
    for (int k = 0; k < 6; k++)
        x *= 2 - n0 * x;                        /* Newton, doubles the bits */
    return (uint64_t)0 - x;
}

int main(void)
{
    for (int k = 0; k < 6; k++) {
        u256 a, b, n, o;
        uint64_t s = 0x9E3779B97F4A7C15ULL + (uint64_t)k * 0xD1B54A32D192ED03ULL;
        for (int j = 0; j < 4; j++) {
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            a.l[j] = s;
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            b.l[j] = s;
            s = s * 6364136223846793005ULL + 1442695040888963407ULL;
            n.l[j] = s;
        }
        n.l[0] |= 1ULL;                         /* odd modulus */
        n.l[3] |= 1ULL << 63;                   /* top bit set  */
        /* the routine's domain is a, b < n (its bound is T < R*n).  Clearing the
         * top limb puts a, b < 2^192, and n >= 2^255 here, so the contract holds
         * without a reduction routine the package would have to price. */
        a.l[3] = 0;
        b.l[3] = 0;
        cios4(&o, &a, &b, &n, n0inv64(n.l[0]));
        printf("%016llx%016llx%016llx%016llx %016llx%016llx%016llx%016llx "
               "%016llx%016llx%016llx%016llx %016llx%016llx%016llx%016llx\n",
               (unsigned long long)n.l[3], (unsigned long long)n.l[2],
               (unsigned long long)n.l[1], (unsigned long long)n.l[0],
               (unsigned long long)a.l[3], (unsigned long long)a.l[2],
               (unsigned long long)a.l[1], (unsigned long long)a.l[0],
               (unsigned long long)b.l[3], (unsigned long long)b.l[2],
               (unsigned long long)b.l[1], (unsigned long long)b.l[0],
               (unsigned long long)o.l[3], (unsigned long long)o.l[2],
               (unsigned long long)o.l[1], (unsigned long long)o.l[0]);
    }
    return 0;
}
#endif
