/* mulmod33.c -- the 255-bit modular product the coefficient row charges.
 *
 * v32 charged 64 Ir for "one 255-bit mulmod" without an implementation or a
 * derivation (verdict 3cec3e8 flaw 3).  This file IS the implementation the
 * ledger prices: four 64-bit limbs, Montgomery reduction with a 64-bit
 * modulus inverse, schoolbook 4x4 product (16 mul) + 4x4 reduction (16 mul).
 * scripts/keys33.py disassembles the shipped object of this function and
 * counts the instructions in its body; the ledger charges that count.
 */
#include <stdint.h>

typedef struct { uint64_t l[4]; } u255;

static inline void rd(u255 *r, const u255 *a, const u255 *b, const u255 *n,
                      uint64_t n0inv) {
    uint64_t t[6] = {0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 4; i++) {
        uint64_t c = 0;
        for (int j = 0; j < 4; j++) {
            unsigned __int128 p = (unsigned __int128)a->l[i] * b->l[j]
                                + t[i + j] + c;
            t[i + j] = (uint64_t)p;
            c = (uint64_t)(p >> 64);
        }
        uint64_t s = t[i + 4] + c;
        t[i + 4] = s;
        t[i + 5] += (s < c);
    }
    for (int i = 0; i < 4; i++) {
        uint64_t m = t[i] * n0inv;
        uint64_t c = 0;
        for (int j = 0; j < 4; j++) {
            unsigned __int128 p = (unsigned __int128)m * n->l[j]
                                + t[i + j] + c;
            t[i + j] = (uint64_t)p;
            c = (uint64_t)(p >> 64);
        }
        uint64_t s = t[i + 4] + c;
        t[i + 4] = s;
        t[i + 5] += (s < c);
        t[i] = 0;
    }
    for (int j = 0; j < 4; j++) r->l[j] = t[4 + j];
    /* one conditional subtraction, branch-free: r -= n if r >= n */
    uint64_t borrow = 0, d[4];
    for (int j = 0; j < 4; j++) {
        unsigned __int128 p = (unsigned __int128)r->l[j] - n->l[j] - borrow;
        d[j] = (uint64_t)p;
        borrow = (uint64_t)(p >> 64) & 1;
    }
    uint64_t mask = borrow - 1;             /* all ones iff no borrow */
    for (int j = 0; j < 4; j++) r->l[j] = (d[j] & mask) | (r->l[j] & ~mask);
}

void mulmod255(u255 *r, const u255 *a, const u255 *b, const u255 *n,
               uint64_t n0inv) {
    rd(r, a, b, n, n0inv);
}
