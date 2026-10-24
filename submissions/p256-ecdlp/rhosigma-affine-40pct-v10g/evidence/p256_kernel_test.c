/* p256_kernel_test.c — exercises the SHIPPED kernel at the REAL P-256 modulus.
 *
 * Reviewer b8b00688 item 1: the 2,000/2,000 bigint comparisons in v10f were made
 * on the separate 40-bit fixture arithmetic (census_v10f41.c), never on the
 * shipped evidence/p256_fe_cc.h.  This driver closes that gap: it includes the
 * shipped header AND the shipped counted step source (for canon_y), prints its
 * raw outputs, and p256_kernel_test.py recomputes every row in Python big
 * integers and diffs.  The same driver is rebuilt against
 * evidence/p256_fe_cc_broken.h (the v10f header, limbs 1/2 transposed) as the
 * mutation control: if the control does not FAIL, the test proves nothing.
 *
 * Output records (all decimal uint64, little-endian limb order):
 *   V a0 a1 a2 a3 b0 b1 b2 b3 m0..m3 p0..p3 s0..s3 t0..t3
 *        m = fe_mul(a,b)  p = fe_add(a,b)  s = fe_sub(a,b)  t = fe_sub(b,a)
 *   I a0 a1 a2 a3 r0 r1 r2 r3        r = fe_inv(a)
 *   C y0 y1 y2 y3 c0 c1 c2 c3        c = canon_y(y)   (shipped step source)
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#include "p256_fe_cc.h"
/* the counted step's own canonicaliser, so the shipped expression is tested */
#include "step_v10f_o3.c"

/* The counted step TU references the production walk's storage; here the symbols
 * only need to exist so that TU links.  They are never read or written: this
 * driver calls canon_y (and the field kernels) only. */
apt TBL[2048], TBL2[2048];
awalker A[16];
uint64_t PROD[16][4], SUF[16][4];

/* Independent reference limbs of p (NOT the header's macros): used only to
 * generate canonical inputs, so the fixed and control builds see identical
 * vectors and the comparison target never moves. */
static const uint64_t PR[4] = {0xffffffffffffffffULL, 0x00000000ffffffffULL,
                               0x0000000000000000ULL, 0xffffffff00000001ULL};

static uint64_t S = 0x243f6a8885a308d3ULL;
static uint64_t nxt(void)
{
    S ^= S << 13;
    S ^= S >> 7;
    S ^= S << 17;
    return S;
}

/* a < 2^256 < 2p, so one conditional subtraction of p makes a canonical */
static void reduce_p(uint64_t a[4])
{
    uint64_t d, t0, t1, t2, t3;
    t0 = __builtin_subcll(a[0], PR[0], 0, &d);
    t1 = __builtin_subcll(a[1], PR[1], d, &d);
    t2 = __builtin_subcll(a[2], PR[2], d, &d);
    t3 = __builtin_subcll(a[3], PR[3], d, &d);
    if (!d) { a[0] = t0; a[1] = t1; a[2] = t2; a[3] = t3; }
}

static void rand_fe(uint64_t a[4])
{
    a[0] = nxt(); a[1] = nxt(); a[2] = nxt(); a[3] = nxt();
    reduce_p(a);
}

/* edge values, little-endian limbs: 0,1,2,p-1,p-2,(p-1)/2,2^63,2^64,2^128,
 * 2^192, 2^224, 2^255, 2^255+1, 2^256-2^224, p-1-2^96 */
#define NEDGE 15
static const uint64_t EDGE[NEDGE][4] = {
    {0ULL, 0ULL, 0ULL, 0ULL},
    {1ULL, 0ULL, 0ULL, 0ULL},
    {2ULL, 0ULL, 0ULL, 0ULL},
    {0xfffffffffffffffeULL, 0x00000000ffffffffULL, 0ULL, 0xffffffff00000001ULL}, /* p-1 */
    {0xfffffffffffffffdULL, 0x00000000ffffffffULL, 0ULL, 0xffffffff00000001ULL}, /* p-2 */
    {0xffffffffffffffffULL, 0x000000007fffffffULL, 0x8000000000000000ULL, 0x7fffffff80000000ULL}, /* (p-1)/2 */
    {0x8000000000000000ULL, 0ULL, 0ULL, 0ULL},                                   /* 2^63 */
    {0ULL, 1ULL, 0ULL, 0ULL},                                                    /* 2^64 */
    {0ULL, 0ULL, 1ULL, 0ULL},                                                    /* 2^128 */
    {0ULL, 0ULL, 0ULL, 1ULL},                                                    /* 2^192 */
    {0ULL, 0ULL, 0ULL, 0x0000000100000000ULL},                                   /* 2^224 */
    {0ULL, 0ULL, 0ULL, 0x8000000000000000ULL},                                   /* 2^255 */
    {1ULL, 0ULL, 0ULL, 0x8000000000000000ULL},                                   /* 2^255+1 */
    {0ULL, 0ULL, 0ULL, 0xffffffff00000000ULL},                                   /* 2^256-2^224 */
    {0xfffffffffffffffeULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffff00000000ULL}, /* p-1-2^96 */
};

#define NRAND 10000

static void emit(const uint64_t a[4], const uint64_t b[4])
{
    uint64_t m[4], s[4], d[4], u[4], v[4];
    fe_mul(m, a, b);
    fe_add(s, a, b);
    fe_sub(d, a, b);
    fe_sub(u, b, a);
    (void)v;
    printf("V");
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, a[i]);
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, b[i]);
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, m[i]);
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, s[i]);
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, d[i]);
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, u[i]);
    printf("\n");
}

static void emit_inv(const uint64_t a[4])
{
    uint64_t r[4];
    fe_inv(r, a);
    printf("I");
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, a[i]);
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, r[i]);
    printf("\n");
}

static void emit_canon(const uint64_t y[4])
{
    uint64_t c[4];
    for (int i = 0; i < 4; i++) c[i] = y[i];
    canon_y(c);
    printf("C");
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, y[i]);
    for (int i = 0; i < 4; i++) printf(" %" PRIu64, c[i]);
    printf("\n");
}

int main(void)
{
    uint64_t a[4], b[4];
    for (int t = 0; t < NRAND; t++) {
        rand_fe(a);
        rand_fe(b);
        emit(a, b);
    }
    for (int i = 0; i < NEDGE; i++) {
        emit(EDGE[i], EDGE[i]);                       /* squares / self-ops */
        for (int j = 0; j < NEDGE; j++) emit(EDGE[i], EDGE[j]);
    }
    for (int t = 0; t < 200; t++) {
        rand_fe(a);
        if (t == 0) { for (int i = 0; i < 4; i++) a[i] = EDGE[0][i]; }   /* 0, not invertible */
        emit_inv(a);
    }
    for (int i = 0; i < NEDGE; i++) emit_inv(EDGE[i]);
    for (int t = 0; t < 300; t++) {
        rand_fe(a);
        emit_canon(a);
    }
    for (int i = 0; i < NEDGE; i++) emit_canon(EDGE[i]);
    return 0;
}
