/* step_v10f_o3.c — v10f step with the O3 two-table attractor fix (counted).

Identical to step_v10f.c except the jump index.  Part 4 measured why the single
index rule cannot work: with one shared table, a repeated index is an INVOLUTION.
The step is S -> canon(S + T_j); when idx(S) repeats, the second step adds the
same T_j to the negated representative and returns to the start, so every walker
enters a 2-cycle at ~2^11 steps (100% of starts, median entry 1,600; measured in
bigint and in C).  That is why x&2047, mix64(x)&2047, a per-walker coset offset
and an xor-fold index all measure identically -- the index rule is irrelevant.

Fix: two independent shared tables, selected by the walker's step parity, each
indexed by a salted mix of the canonical x.  A repeat now needs
T^A_j == T^B_j' as group elements (~2^-40), not an 11-bit index collision.

Cost over v10f: parity extract (1), salt select (2 csel), salted mix64 (6),
table-base select (1).
*/
#include <stdint.h>
typedef unsigned __int128 u128;
#include "p256_fe_cc.h"

static inline uint64_t mix64(uint64_t v)
{ v *= 0x9e3779b97f4a7c15ULL; v ^= v >> 29; v *= 0xbf58476d1ce4e5b9ULL; v ^= v >> 32; return v; }

typedef struct { uint64_t x[4], y[4]; } apt;
extern apt TBL[2048], TBL2[2048];
typedef struct { uint64_t x[4], y[4]; uint64_t len, esc; } awalker;
extern awalker A[16];
extern uint64_t PROD[16][4], SUF[16][4];

static const uint64_t SA = 0xa5a5a5a5a5a5a5a5ULL, SB = 0x5a5a5a5a5a5a5a5aULL;

/* y <- even-y representative of {y, p-y}: mvn + 3 sbc + tst + 4 csel = 9 */
__attribute__((always_inline)) static inline void canon_y(uint64_t y[4])
{
    uint64_t c, n0, n1, n2, n3;
    n0 = ~y[0];
    n1 = __builtin_subcll(P1, y[1], 0, &c);
    n2 = __builtin_subcll(P2, y[2], c, &c);
    n3 = __builtin_subcll(P3, y[3], c, &c);
    __asm__("tst %4, #1\n\t"
            "csel %0, %0, %5, eq\n\t"
            "csel %1, %1, %6, eq\n\t"
            "csel %2, %2, %7, eq\n\t"
            "csel %3, %3, %8, eq"
            : "+&r"(y[0]), "+&r"(y[1]), "+&r"(y[2]), "+&r"(y[3])
            : "r"(y[0]), "r"(n0), "r"(n1), "r"(n2), "r"(n3));
}

__attribute__((noinline)) void step_v10f_o3(int n)
{
    for (int i = 0; i < n; i++) {
        uint64_t x1[4], y1[4];
        for (int k = 0; k < 4; k++) { x1[k] = A[i].x[k]; y1[k] = A[i].y[k]; }
        uint64_t par = A[i].len & 1;
        const apt *tab = par ? TBL2 : TBL;
        uint64_t j = mix64(x1[0] ^ (par ? SB : SA)) & 2047;
        uint64_t d[4], num[4], inv[4], lam[4], x3[4], y3[4], t[4];
        fe_sub(d, tab[j].x, x1);
        fe_sub(num, tab[j].y, y1);
        fe_mul(PROD[i], PROD[i], d);
        fe_mul(SUF[i], SUF[i], d);
        fe_mul(inv, PROD[i], SUF[i]);
        fe_mul(lam, num, inv);
        fe_mul(t, lam, lam);
        fe_sub(x3, t, x1); fe_sub(x3, x3, tab[j].x);
        fe_sub(t, x1, x3);
        fe_mul(y3, lam, t);
        fe_sub(y3, y3, y1);
        canon_y(y3);
        uint64_t dp = mix64(x3[0] ^ 0xd1b54a32d192ed03ULL);
        uint64_t len = A[i].len + 1, esc = A[i].esc | ((dp & ((1ULL << 61) - 1)) == 0);
        for (int k = 0; k < 4; k++) { A[i].x[k] = x3[k]; A[i].y[k] = y3[k]; }
        A[i].len = len; A[i].esc = esc;
    }
}
