/* rec_layout.h -- THE canonical record layout, v15p.
 *
 * ONE definition, included by the shipped module (walk_v15o.c), by the
 * layout assertion TU (layout_assert.c) and checked against the assembler
 * constants of evidence/step_secp_v15g.s by scripts/check_layout.py.
 *
 * The record is 88 bytes, NOT 80: x[4] + y[4] + idx + len + esc = 11 u64.
 * The y coordinate is needed mid-step (the point step adds the table point to
 * it), so the honest fix is the 88-byte record and the listing's stride and
 * field offsets are corrected to it -- never a forward-only claim about 80.
 */
#ifndef REC_LAYOUT_V15P_H
#define REC_LAYOUT_V15P_H

typedef unsigned long long u64;
typedef unsigned __int128 u128;

typedef struct {
    u64 x[4];      /*   0 .. 31  canonical class key (min of the 3 phi-images) */
    u64 y[4];      /*  32 .. 63  point y, after the alignment's sign choice    */
    u64 idx;       /*  64        table row, DERIVED FROM THE KEY ONLY           */
    u64 len;       /*  72        steps since the last report                    */
    u64 esc;       /*  80        escape/flag word                               */
} rec_t;           /*  88 bytes                                                 */

#define REC_SIZE      88u
#define REC_OFF_X      0u
#define REC_OFF_Y     32u
#define REC_OFF_IDX   64u
#define REC_OFF_LEN   72u
#define REC_OFF_ESC   80u

/* Compile-time law: the struct and the stride the asm walks are ONE number.
 * Any drift is a build failure, not a review finding. */
_Static_assert(sizeof(rec_t) == REC_SIZE, "rec_t must be exactly 88 bytes");
_Static_assert(sizeof(rec_t) % 8 == 0, "rec_t must be a whole number of u64");
_Static_assert(__builtin_offsetof(rec_t, x) == REC_OFF_X,   "x offset");
_Static_assert(__builtin_offsetof(rec_t, y) == REC_OFF_Y,   "y offset");
_Static_assert(__builtin_offsetof(rec_t, idx) == REC_OFF_IDX, "idx offset");
_Static_assert(__builtin_offsetof(rec_t, len) == REC_OFF_LEN, "len offset");
_Static_assert(__builtin_offsetof(rec_t, esc) == REC_OFF_ESC, "esc offset");
_Static_assert(REC_OFF_ESC + 8 == REC_SIZE, "esc must be the last u64");


/* ---- worker state (verdict flaw 4, v15s): the LIVE per-worker struct. ------
 * The v15r ledger charged 88 B/worker, i.e. rec_t alone, and omitted the live
 * label coefficient pair, the detector history ring and the counters.  ONE
 * definition here; sizeof/offsetof are machine-derived by layout_assert.c.
 * The batch-inversion prefix array (BATCH x 64 B) is a separate scratch block
 * charged in the peak-storage row, not part of the struct. */
typedef struct {
    rec_t rec;              /*   0 ..  87  point + table row + step + flags   */
    u64 a[4];               /*  88 .. 119  label coefficient a               */
    u64 c[4];               /* 120 .. 151  label coefficient c               */
    u64 prev0[4];           /* 152 .. 183  detector history (previous key)   */
    u64 prev1[4];           /* 184 .. 215  detector history (older key)      */
    u64 minx[4];            /* 216 .. 247  detector window minimum key       */
    u64 rep, win;           /* 248 .. 263  detector counters                 */
} worker_t;                 /* 264 bytes                                     */

#define WORKER_SIZE      264u
#define WORKER_BATCH     67u          /* PROD/SUF steps batched per inversion   */
#define PREFIX_BYTES     (WORKER_BATCH * 64u)   /* Montgomery batch prefix array */

_Static_assert(sizeof(worker_t) == WORKER_SIZE, "worker_t must be exactly 264 bytes");
_Static_assert(__builtin_offsetof(worker_t, a) == 88, "worker a offset");
_Static_assert(__builtin_offsetof(worker_t, c) == 120, "worker c offset");
_Static_assert(__builtin_offsetof(worker_t, prev0) == 152, "worker prev0 offset");
_Static_assert(__builtin_offsetof(worker_t, prev1) == 184, "worker prev1 offset");
_Static_assert(__builtin_offsetof(worker_t, minx) == 216, "worker minx offset");
_Static_assert(__builtin_offsetof(worker_t, rep) == 248, "worker rep offset");
_Static_assert(__builtin_offsetof(worker_t, win) == 256, "worker win offset");
_Static_assert(WORKER_SIZE == REC_SIZE + 8u * 22u, "worker = rec + 22 u64 of live state");

#endif
