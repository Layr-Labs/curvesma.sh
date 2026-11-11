/* Priced tag-tracing step loop for the BN254 GT walk (see machine.md).
 *
 * State of a window: y = the 48 Montgomery limbs of the window base Y0
 * (12 Fp coordinates, 4 limbs each); e = current table entry, standing for
 * the group element Y0*T_e.  Each entry stores, for every limb (k,i) of y,
 * the 128-bit fraction F = floor(2^128 * u/p), u = v_k*2^(64i-256) mod p,
 * where v is the row with  coord0(Y*T_e) = sum_k y_k v_k (mod p).
 * One iteration computes A = sum (umulh(y,Flo) + y*Fhi) mod 2^64, which
 * satisfies  2^64*t/p in [A, A+96)  for the exact coordinate t (proof.md),
 * reads the 3 selector bits from the top of A when they are certified, and
 * follows one stored pointer.  No field multiplication is performed.
 *
 * The loop is the only loop; its three exits are not-taken branches on the
 * priced path.  The trip budget is a 64-bit count supplied by the driver,
 * which keeps the 128-bit global and per-walk counters (machine.md).
 */
#include <stdint.h>
typedef uint64_t word;
typedef __uint128_t wide;
typedef struct entry {
    word F[96];              /* (Flo,Fhi) for limb j at F[2j], F[2j+1]   */
    struct entry *next[8];   /* child for selector s; 0 at window depth   */
    word T[48];              /* T_e, Montgomery limbs, for the re-base    */
    word alpha[4], beta[4];  /* label increments mod r                    */
} entry;                     /* 1280 bytes                                */
typedef struct { const entry *e; word A; word budget; } cursor;
enum { RUN_BUDGET = 0, RUN_EXACT = 1, RUN_REBASE = 2 };
#define ERR 96u
#define MASK61 0x1fffffffffffffffULL
#define TERM(j) (word)(((wide)y[j] * F[2*(j)]) >> 64) + y[j] * F[2*(j)+1]
#define T4(j) (TERM(j) + TERM(j+1) + TERM(j+2) + TERM(j+3))
#define T16(j) (T4(j) + T4(j+4) + T4(j+8) + T4(j+12))

__attribute__((noinline))
int tt_run(cursor *c, const word *y) {
    const entry *e = c->e;
    word budget = c->budget, A;
    for (;;) {
        const word *F = e->F;
        const entry *n;
        A = T16(0) + T16(16) + T16(32);
        if (__builtin_expect(((A + ERR) & MASK61) <= ERR, 0)) { c->e = e; c->A = A; c->budget = budget; return RUN_EXACT; }
        n = e->next[A >> 61];
        if (__builtin_expect(n == 0, 0)) { c->e = e; c->A = A; c->budget = budget; return RUN_REBASE; }
        e = n;
        if (__builtin_expect(--budget == 0, 0)) { c->e = e; c->A = A; c->budget = 0; return RUN_BUDGET; }
    }
}
