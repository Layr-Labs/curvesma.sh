/* cios4.c -- 4x64-bit Montgomery multiplication (CIOS), the kernel the ledger's
 * M_insn row prices.  Build and count:
 *
 *     clang -O2 -c cios4.c -o cios4.o
 *     python3 scripts/row_extract.py          # disassembles, counts the span,
 *                                             # asserts rows == ledger.M_insn
 *
 * The point of this file is that the instruction count is READ from the
 * disassembly of the compiled routine (otool -tv), not asserted in prose.
 * Nothing here is executed by the submission; it is the counting artifact.
 */
#include <stdint.h>

typedef struct { uint64_t l[4]; } u256;

/* t[j] += a*b + c ; returns the new carry. */
static inline uint64_t mac(uint64_t *t, uint64_t a, uint64_t b, uint64_t c)
{
    unsigned __int128 s = (unsigned __int128)a * b + *t + c;
    *t = (uint64_t)s;
    return (uint64_t)(s >> 64);
}

/* out = a * b * R^{-1} mod n,  R = 2^256,  n0inv = -n[0]^{-1} mod 2^64. */
__attribute__((noinline, visibility("default")))
void cios4(u256 *out, const u256 *a, const u256 *b, const u256 *n, uint64_t n0inv)
{
    uint64_t t[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
        uint64_t c = 0;
        for (int j = 0; j < 4; j++)
            c = mac(&t[j], a->l[i], b->l[j], c);
        uint64_t s = t[4] + c;
        t[4] = s;
        t[5] = (s < c);

        uint64_t m = t[0] * n0inv;
        c = mac(&t[1], m, n->l[1], mac(&t[0], m, n->l[0], 0));
        c = mac(&t[2], m, n->l[2], c);
        c = mac(&t[3], m, n->l[3], c);
        s = t[4] + c;
        t[3] = s;
        t[4] = t[5] + (s < c);
        t[5] = 0;
    }
    uint64_t borrow = 0;
    uint64_t r[4];
    for (int j = 0; j < 4; j++) {              /* final conditional subtract */
        unsigned __int128 d = (unsigned __int128)t[j] - n->l[j] - borrow;
        r[j] = (uint64_t)d;
        borrow = (uint64_t)(d >> 64) & 1;
    }
    uint64_t mask = ~(uint64_t)0 + borrow;     /* all ones when t >= n */
    for (int j = 0; j < 4; j++)
        out->l[j] = (r[j] & mask) | (t[j] & ~mask);
}

/* span fence: the extractor counts from _cios4 up to this label */
__attribute__((noinline, visibility("default")))
void cios4_end(void) { }
