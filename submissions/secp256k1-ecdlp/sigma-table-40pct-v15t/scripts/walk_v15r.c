/* walk_v15o.c -- SHIPPED per-step bookkeeping of the six-class sigma walk, v15p.
 *
 * Repairs every executable finding of the v15o verdict (evidence/verdict-v15o.md):
 *
 *  (1) LAYOUT.  One record type, defined ONCE in scripts/rec_layout.h with
 *      compile-time asserts: rec_t = {x[4], y[4], idx, len, esc} = 88 bytes,
 *      offsets 0/32/64/72/80.  The listing evidence/step_secp_v15g.s walks
 *      REC_SIZE (88) and stores len at 72 and esc at 80 -- previously it strode
 *      80 and used offset 64 as BOTH the index and the length, which is the
 *      stride AND field-semantics defect the verdict names.  The y coordinate
 *      is genuinely needed mid-step (the point step adds the table point to
 *      it), so 88 is the minimal honest layout, not 80.
 *
 *  (2) LABELS.  v15o_label now
 *        - uses R mod r as Montgomery ONE (integer 1 multiplied under Mont
 *          reduction gives a*R^-1, not a -- the old slot-0 row was wrong);
 *        - negates with an UNCONDITIONAL r - m for the multiplier m, which is
 *          a unit (lambda^k), so the negation is exact for every m != 0.  The
 *          old code stored r-a only when the subtraction borrowed, which never
 *          happens for a reduced a, so it never negated;
 *        - ADDS the selected jump's exponent pair (ju, jv) mod r before the
 *          scaling: the invariant is A = [a]G + [c]Q with A' = s*phi^-slot
 *          (A + J), so the label must carry the jump, not just the scaling;
 *        - consumes the SIGN chosen by v15o_align (even-y canonicalisation
 *          negates the point, so the multiplier is negated with it).  The
 *          point step no longer forces even y on its own: alignment is the ONE
 *          canonicalisation site, and it reports the sign it used.
 *
 *  (3) INDEX.  The table row is derived from the CANONICAL KEY ONLY
 *      (idx = mix(key) & (R-1)); the argmin slot is no longer part of the
 *      address.  The slot changes between orbit partners by construction, so
 *      an address that contains it is not a class function; the slot is still
 *      needed, but only to select lambda^-slot in the LABEL, where its effect
 *      is proven by the shipped partner walk in scripts/label_invariant.py.
 *
 *  (4) DETECTOR.  Two-step full-key comparison plus a running-minimum
 *      checkpoint over a WINDOW counter, matching the escape the ledger prices.
 *
 *  (5) REPORT POOL.  Probe count is capped at PROBE_MAX (a deterministic
 *      bound, not a load-factor premise) and a full bucket replaces its last
 *      probed node instead of refusing.
 */
#include "rec_layout.h"

/* v15q: the field constants are the ONE thing a fixture build may override, so
 * the same shipped source is compiled for secp256k1 (harness) and for the small
 * census fixture.  Every other line is identical in both builds. */
#ifndef V15Q_P0
#define V15Q_P0 0xFFFFFFFEFFFFFC2FULL
#endif
#ifndef V15Q_C
#define V15Q_C 0x1000003D1ULL
#endif

extern const u64 BETA[4];       /* cube root of unity mod p    (driver TU) */
extern u64 NODE_HEAD[];         /* pool: 12 limbs per node, 96 B (driver TU) */
extern u64 NODE_MASK;           /* nodes - 1                   (driver TU) */
extern u64 R_CONST[4];          /* group order r               (driver TU) */
extern u64 R_N0;                /* -r^-1 mod 2^64              (driver TU) */
extern u64 LAMBDA_R[4];         /* lambda   in Mont form (driver TU) */
extern u64 LAMBDA2_R[4];        /* lambda^2 in Mont form (driver TU) */
extern u64 R_MONT_ONE[4];       /* R mod r  == Montgomery one  (driver TU) */

#define C_LIMB V15Q_C
static const u64 P0 = V15Q_P0;
#define ROWS 2048u
/* table WIDTH is the one resource the census fixture scales (2048 -> 128);
 * the production build leaves it at ROWS.  The BUDGET formula always uses the
 * production R in both, so the normalization row compares like with like. */
#ifndef V15R_ROWS
#define V15R_ROWS ROWS
#endif
/* Mutation switches for scripts/label_invariant.py.  They are compiled-in only
 * by the harness; the shipped object is built without any of them.  Each one
 * re-introduces exactly one v15o defect so the invariant test has a control
 * that must FAIL. */
#ifdef MUT_MONT_ONE
#define M_ONE_IS_INTEGER 1
#endif
#ifdef MUT_NEG_BORROW
#define M_NEG_ON_BORROW 1
#endif
#ifdef MUT_NO_JUMP
#define M_DROP_JUMP 1
#endif
#ifdef MUT_NO_SIGN
#define M_DROP_SIGN 1
#endif
#ifdef MUT_IDX_SLOT
#define M_IDX_SLOT 1
#endif
#ifdef MUT_IDX_RAW
#define M_IDX_RAW 1
#endif
#define PROBE_MAX 4u
#define WINDOW 67u

static inline int fcmp(const u64 a[4], const u64 b[4]) {
    #pragma clang loop unroll(full)
    for (int i = 3; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

static inline int eq4(const u64 a[4], const u64 b[4]) {
    return ((a[0]^b[0]) | (a[1]^b[1]) | (a[2]^b[2]) | (a[3]^b[3])) ? 0 : 1;
}

#ifndef V15Q_P0
#define V15Q_P0 0xFFFFFFFEFFFFFC2FULL
#endif
#ifndef V15Q_C
#define V15Q_C 0x1000003D1ULL
#endif

/* one field product mod p = 2^256 - C, C = 0x1000003D1 (secp256k1 fold).
 * the fold overflowed), caught by scripts/label_invariant.py against Python big
 * integers.  This version folds limb by limb from the top: 2^{64i} = C*2^{64(i-4)}
 * for i >= 4, and the carry of hi*C can never reach limb i-1, so one descending
 * pass leaves v[4..7] == 0. */
static inline void fmul(u64 r[4], const u64 a[4], const u64 b[4]) {
#ifdef V15Q_SMALL_FIELD
    /* fixture build: p < 2^32, so one u128 product + one remainder is exact.
     * Same function contract, same call sites; only the width differs. */
    u128 w = (u128)a[0] * (u128)b[0];
    u64 z = (u64)(w % (u128)P0);
    r[0] = z; r[1] = 0; r[2] = 0; r[3] = 0;
    (void)C_LIMB;
    return;
#else
    u64 v[8] = {0};
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) {
        u64 cy = 0;
        #pragma clang loop unroll(full)
        for (int j = 0; j < 4; j++) {
            u128 s = (u128)a[i] * b[j] + v[i + j] + cy;
            v[i + j] = (u64)s;
            cy = (u64)(s >> 64);
        }
        #pragma clang loop unroll(full)
        for (int k = i + 4; k < 8 && cy; k++) {
            u128 s = (u128)v[k] + cy;
            v[k] = (u64)s;
            cy = (u64)(s >> 64);
        }
    }
    #pragma clang loop unroll(full)
    for (int i = 7; i >= 4; i--) {
        u64 hi = v[i];
        v[i] = 0;
        u64 cy = 0;
        u128 s = (u128)hi * C_LIMB + v[i - 4];
        v[i - 4] = (u64)s; cy = (u64)(s >> 64);
        s = (u128)v[i - 3] + cy; v[i - 3] = (u64)s; cy = (u64)(s >> 64);
        s = (u128)v[i - 2] + cy; v[i - 2] = (u64)s;   /* cy == 0 here */
    }
    u64 z[4] = {v[0], v[1], v[2], v[3]};
    #pragma clang loop unroll(full)
    for (int pass = 0; pass < 2; pass++) {
        u64 d[4]; u64 br = 0;
        #pragma clang loop unroll(full)
        for (int i = 0; i < 4; i++) {
            u128 w = (u128)z[i] - (i == 0 ? P0 : ~0ULL) - br;
            d[i] = (u64)w; br = (u64)(w >> 64) & 1;
        }
        if (!br) { z[0]=d[0]; z[1]=d[1]; z[2]=d[2]; z[3]=d[3]; }
    }
    r[0]=z[0]; r[1]=z[1]; r[2]=z[2]; r[3]=z[3];
#endif
}

/* ---- v15o_align: canonical class key, CLASS-DEFINED row, sign, slot ----
 * key = min(x, BETA*x, BETA^2*x) is invariant on the six-element orbit and
 * (key, even y) is one representative per class.  idx = mix(key) & (ROWS-1)
 * uses canonical-key bits ONLY, so two orbit partners select the same row by
 * construction (the partner walk in scripts/label_invariant.py measures it).
 * The argmin slot is returned for the LABEL only (it selects lambda^-slot);
 * it is never part of an address.  Returns the sign of y in *sgn after
 * replacing y by its even representative.
 */
__attribute__((noinline))
u64 v15o_align(rec_t *rec, int *sgn) {
    u64 b1[4], b2[4];
    fmul(b1, rec->x, BETA);
    fmul(b2, b1, BETA);
    u64 key[4]; u64 slot = 0;
    const u64 *cand[3] = {rec->x, b1, b2};
    key[0]=cand[0][0]; key[1]=cand[0][1]; key[2]=cand[0][2]; key[3]=cand[0][3];
    #pragma clang loop unroll(full)
    for (int i = 1; i < 3; i++) {
        if (fcmp(cand[i], key) < 0) { key[0]=cand[i][0]; key[1]=cand[i][1]; key[2]=cand[i][2]; key[3]=cand[i][3]; slot = (u64)i; }
    }
#if defined(M_IDX_SLOT) || defined(M_IDX_RAW)
    u64 raw0 = rec->x[0];
#endif
    rec->x[0]=key[0]; rec->x[1]=key[1]; rec->x[2]=key[2]; rec->x[3]=key[3];
#if defined(M_IDX_RAW)
    rec->idx = 3 * (raw0 & 2047) + slot;        /* input-dependent address */
#elif defined(M_IDX_SLOT)
    /* v15p DEFECT, kept as a mutation control: the argmin slot of the PRE-
     * canonical point is folded into the address, so two representatives of one
     * class store the same (x,y) but different slots -> different next jump. */
    rec->idx = 3 * (((key[0] ^ (key[1] * 0x9E3779B97F4A7C15ULL)) >> 3) & (V15R_ROWS - 1)) + slot;
#else
    /* v15q: the row is a function of the CANONICAL KEY ALONE.  slot is returned
     * for the LABEL only (it names the phi power the canonicalisation applied);
     * it is never an address bit.  With slot out of the address, the stored
     * (x,y) determines idx exactly: two orbit partners store the SAME point and
     * therefore select the SAME row, so the transition is a map on classes -
     * which is what the partner harness now measures on the STORED points. */
    rec->idx = (((key[0] ^ (key[1] * 0x9E3779B97F4A7C15ULL)) >> 3) & (V15R_ROWS - 1));
#endif
    *sgn = (int)(rec->y[0] & 1);              /* 1 => negate for the even rep */
    if (*sgn) {                                /* y <- p - y */
#ifdef V15Q_SMALL_FIELD
        rec->y[0] = P0 - rec->y[0];
        rec->y[1] = 0; rec->y[2] = 0; rec->y[3] = 0;
#else
        u64 d[4]; u64 br = 0; u128 w;
        w = (u128)P0 - rec->y[0];        d[0]=(u64)w; br=(u64)(w>>64)&1;
        w = (u128)~0ULL - rec->y[1] - br; d[1]=(u64)w; br=(u64)(w>>64)&1;
        w = (u128)~0ULL - rec->y[2] - br; d[2]=(u64)w; br=(u64)(w>>64)&1;
        w = (u128)~0ULL - rec->y[3];      d[3]=(u64)w;
        rec->y[0]=d[0]; rec->y[1]=d[1]; rec->y[2]=d[2]; rec->y[3]=d[3];
#endif
    }
    return slot;
}

/* Montgomery product mod r: out = a*b*R^-1 mod r.
 * Textbook CIOS with a two-limb top accumulator.  The v15o body carried the
 * t[4] overflow into t[3] and then cleared t[4], so it silently dropped the
 * top carry (measured 115/200 mismatches against Python big integers). */
static inline void rmont(u64 out[4], const u64 a[4], const u64 b[4]) {
    u64 t[6] = {0,0,0,0,0,0};
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) {
        u64 cy = 0;
        #pragma clang loop unroll(full)
        for (int j = 0; j < 4; j++) {
            u128 s = (u128)a[j] * b[i] + t[j] + cy;
            t[j] = (u64)s; cy = (u64)(s >> 64);
        }
        u128 s = (u128)t[4] + cy;
        t[4] = (u64)s; t[5] = (u64)(s >> 64);
        u64 m = t[0] * R_N0;
        s = (u128)m * R_CONST[0] + t[0];
        cy = (u64)(s >> 64);
        #pragma clang loop unroll(full)
        for (int j = 1; j < 4; j++) {
            s = (u128)m * R_CONST[j] + t[j] + cy;
            t[j - 1] = (u64)s; cy = (u64)(s >> 64);
        }
        s = (u128)t[4] + cy;
        t[3] = (u64)s;
        t[4] = t[5] + (u64)(s >> 64);
    }
    /* CIOS leaves a value < 2r, which does NOT fit four limbs (r ~ 2^256):
     * t[4] is part of the result, so the final subtraction is 5-limb. */
    u64 z[5] = {t[0], t[1], t[2], t[3], t[4]};
    u64 d[5]; u64 br = 0;
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) {
        u128 w = (u128)z[i] - R_CONST[i] - br;
        d[i] = (u64)w; br = (u64)(w >> 64) & 1;
    }
    { u128 w = (u128)z[4] - br; d[4] = (u64)w; br = (u64)(w >> 64) & 1; }
    if (!br) { z[0]=d[0]; z[1]=d[1]; z[2]=d[2]; z[3]=d[3]; }
    out[0]=z[0]; out[1]=z[1]; out[2]=z[2]; out[3]=z[3];
    (void)z[4];
}

/* x += y mod r for x,y < r: 4-limb add, then ONE conditional subtract of r. */
static inline void rmod_add(u64 x[4], const u64 y[4]) {
    u64 s[4]; u64 cy = 0;
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) { u128 w = (u128)x[i] + y[i] + cy; s[i] = (u64)w; cy = (u64)(w >> 64); }
    u64 d[4]; u64 br = 0;
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) { u128 w = (u128)s[i] - R_CONST[i] - br; d[i] = (u64)w; br = (u64)(w >> 64) & 1; }
    if (cy || !br) { x[0]=d[0]; x[1]=d[1]; x[2]=d[2]; x[3]=d[3]; }
    else           { x[0]=s[0]; x[1]=s[1]; x[2]=s[2]; x[3]=s[3]; }
}

/* m <- r - m, UNCONDITIONAL.  m is a unit (lambda^k), so m != 0 and r - m is
 * the exact additive inverse in [1, r-1]; there is no borrow condition to test
 * and no reduced nonzero value escapes the negation. */
static inline void rneg_unit(u64 m[4]) {
    u64 d[4]; u64 br = 0;
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) { u128 w = (u128)R_CONST[i] - m[i] - br; d[i] = (u64)w; br = (u64)(w >> 64) & 1; }
#ifdef M_NEG_ON_BORROW
    if (br) { m[0]=d[0]; m[1]=d[1]; m[2]=d[2]; m[3]=d[3]; }   /* stores only on borrow */
#else
    m[0]=d[0]; m[1]=d[1]; m[2]=d[2]; m[3]=d[3];
#endif
}

/* ---- v15o_label: (a,c) <- s * lambda^-slot * ((a,c) + (ju,jv)) mod r ----
 * The step is A' = canon(A + J) with canon = s * phi^-slot, phi = [lambda].
 * So the label MUST add the jump's exponents and then scale by s*lambda^-slot.
 * slot selects lambda^-slot = {1, lambda^2, lambda} in Montgomery form (the
 * multiplier is a table constant, not a product); the sign multiplies it. */
__attribute__((noinline))
void v15o_label(u64 a[4], u64 c[4], u64 slot, int sgn, const u64 ju[4], const u64 jv[4]) {
    u64 m[4];
    if (slot == 0)      {
#ifdef M_ONE_IS_INTEGER
        m[0]=1; m[1]=0; m[2]=0; m[3]=0;         /* v15o: integer one, not R mod r */
#else
        m[0]=R_MONT_ONE[0]; m[1]=R_MONT_ONE[1]; m[2]=R_MONT_ONE[2]; m[3]=R_MONT_ONE[3];
#endif
    }
    /* align WRITES x' = beta^slot * x, and phi(P) = [lambda]P, so the point it
     * produced is [lambda^+slot]P -- the multiplier is lambda^+slot, NOT
     * lambda^-slot (v15o and the first v15p build had the direction inverted;
     * the verdict names this as wrong even with repaired multiplication). */
    else if (slot == 1) { m[0]=LAMBDA_R[0];   m[1]=LAMBDA_R[1];   m[2]=LAMBDA_R[2];   m[3]=LAMBDA_R[3]; }
    else                { m[0]=LAMBDA2_R[0];  m[1]=LAMBDA2_R[1];  m[2]=LAMBDA2_R[2];  m[3]=LAMBDA2_R[3]; }
#ifndef M_DROP_SIGN
    if (sgn) rneg_unit(m);
#endif
#if !defined(M_DROP_JUMP)
    rmod_add(a, ju);
    rmod_add(c, jv);
#endif
    rmont(a, a, m);
    rmont(c, c, m);
}

/* ---- v15o_detect: dp1 predicate, TWO-step key comparison, running-min
 * checkpoint over a WINDOW counter -- the detector the ledger's escape row
 * prices.  State: prev0/prev1 (two previous canonical keys), minx (running
 * minimum), win (steps since the minimum last moved), rep (repeat count). */
__attribute__((noinline))
u64 v15o_detect(const rec_t *rec, u64 *prev0, u64 *prev1, u64 *rep, u64 *minx, u64 *win) {
    u64 dp = 1 - (rec->x[0] & 1);                  /* 1 => report (density 2^-1) */
    if (eq4(rec->x, prev0) || eq4(rec->x, prev1)) *rep += 1;
    if (fcmp(rec->x, minx) < 0) {
        minx[0]=rec->x[0]; minx[1]=rec->x[1]; minx[2]=rec->x[2]; minx[3]=rec->x[3];
        *win = 0;
    }
    *win += 1;
    if (*win > WINDOW) *win = WINDOW;              /* escape armed, priced in ledger */
    prev1[0]=prev0[0]; prev1[1]=prev0[1]; prev1[2]=prev0[2]; prev1[3]=prev0[3];
    prev0[0]=rec->x[0]; prev0[1]=rec->x[1]; prev0[2]=rec->x[2]; prev0[3]=rec->x[3];
    return dp;
}

/* ---- v15o_report: retained-record pool, BOUNDED probing ----
 * Node = 12 limbs: key 0..3, a 4..7, c 8..11 -- 96 bytes, addressed as
 * NODE_HEAD + 12*node (limb stride).  At most PROBE_MAX nodes are touched per
 * call, so the per-report instruction row is a deterministic bound and needs
 * no load-factor premise; a full bucket REPLACES its last probed node (retention
 * policy: newest-in-bucket), so the routine never refuses and never scans the
 * whole table.  Returns 1 on a key match (pair delivered), 0 otherwise. */
__attribute__((noinline))
u64 v15o_report(u64 *pool, u64 mask, const rec_t *rec, const u64 a[4], const u64 c[4],
                u64 out_a[4], u64 out_c[4]) {
    u64 h = (rec->x[0] ^ (rec->x[1] * 0x9E3779B97F4A7C15ULL)) & mask;
    u64 *node = pool + 12 * h;
    #pragma clang loop unroll(full)
    for (u64 probe = 0; probe < PROBE_MAX; probe++) {
        node = pool + 12 * ((h + probe) & mask);
        u64 k0 = node[0], k1 = node[1], k2 = node[2], k3 = node[3];
        if ((k0 | k1 | k2 | k3) == 0) {            /* empty: store and return */
            node[0]=rec->x[0]; node[1]=rec->x[1]; node[2]=rec->x[2]; node[3]=rec->x[3];
            node[4]=a[0]; node[5]=a[1]; node[6]=a[2]; node[7]=a[3];
            node[8]=c[0]; node[9]=c[1]; node[10]=c[2]; node[11]=c[3];
            return 0;
        }
        if (k0==rec->x[0] && k1==rec->x[1] && k2==rec->x[2] && k3==rec->x[3]) {
            out_a[0]=node[4]; out_a[1]=node[5]; out_a[2]=node[6]; out_a[3]=node[7];
            out_c[0]=node[8]; out_c[1]=node[9]; out_c[2]=node[10]; out_c[3]=node[11];
            return 1;
        }
    }
    /* bucket full and no match: replace the last probed node */
    node[0]=rec->x[0]; node[1]=rec->x[1]; node[2]=rec->x[2]; node[3]=rec->x[3];
    node[4]=a[0]; node[5]=a[1]; node[6]=a[2]; node[7]=a[3];
    node[8]=c[0]; node[9]=c[1]; node[10]=c[2]; node[11]=c[3];
    return 0;
}

/* ================= v15r additions (verdict 1b2ee652, flaws 1 and 2) ========
 *
 * FLAW 1 -- PRECOMPUTATION BAKED TO Q0.  The v15q table was built once, from a
 * FIXED Q0 = [K0]P0, while the census drew a fresh challenge Q for every trial.
 * The jump (u_k, v_k) then advanced the point by [u_k]P0 + [v_k]Q0 but the label
 * added (u_k, v_k) to the coefficients of (P0, Q) -- a wrong pair for every
 * trial but the first.  v15r_init_for_Q is the SHIPPED per-challenge entry
 * point: it rebuilds the table as
 *
 *     J_k = [u_k]P0 + [v_k]Q        v_k = 1 + (k mod 15)
 *
 * so the label's (u_k, v_k) is the true decomposition in the (P0, Q) basis for
 * ANY challenge Q.  The [u_k]P0 half is challenge-independent and lives in a
 * one-time static table (UBASE, charged once); the per-challenge work is the 15
 * small multiples of Q plus one point addition per row -- priced as its own
 * full one-time accounting row (INIT_FOR_Q), charged once PER CHALLENGE.
 *
 * FLAW 2 -- THE DETECTOR COUNTERS NEVER ACTED.  v15o_detect maintained *rep and
 * *win, and nothing branched on them.  v15r_action turns the counters into the
 * mechanism the escape row prices, and v15r_reset_hist performs the state reset
 * (the caller then respawns the walk from fresh coins).  MUT_DETACH_ACTION
 * keeps the counter and removes the branch, so the reset test has a control
 * that must FAIL.
 */
extern u64 UBASE[];              /* 8*ROWS limbs: static [u_k]P0 (driver TU) */
extern u64 SMALL[];              /* 8*16 limbs: per-challenge [v]Q (this TU) */
extern int v15q_add(u64*, u64*, const u64*, const u64*, const u64*, const u64*);

#ifndef V15R_REP_THRESH
#define V15R_REP_THRESH 16u
#endif
#define REP_THRESH V15R_REP_THRESH

__attribute__((noinline))
u64 v15r_init_for_Q(u64 *tbl, const u64 qx[4], const u64 qy[4]) {
    u64 fails = 0;
    /* SMALL[v] = [v]Q for v = 1..15: SMALL[v] = SMALL[v-1] + Q */
    for (int i = 0; i < 4; i++) { SMALL[8 + i] = qx[i]; SMALL[12 + i] = qy[i]; }
    for (u64 v = 2; v <= 15; v++) {
        u64 *dst = SMALL + 8 * v, *src = SMALL + 8 * (v - 1);
        u64 rx[4], ry[4];
        if (v15q_add(rx, ry, src, src + 4, SMALL + 8, SMALL + 12) != 0) { fails++; continue; }
        for (int i = 0; i < 4; i++) { dst[i] = rx[i]; dst[4 + i] = ry[i]; }
    }
    for (u64 k = 0; k < V15R_ROWS; k++) {
        u64 v = 1 + (k % 15);
        u64 rx[4], ry[4];
        const u64 *a = UBASE + 8 * k;
        const u64 *b = SMALL + 8 * v;
        if (v15q_add(rx, ry, a, a + 4, b, b + 4) != 0) { fails++; rx[0]=ry[0]=0; }
        u64 *o = tbl + 8 * k;
        for (int i = 0; i < 4; i++) { o[i] = rx[i]; o[4 + i] = ry[i]; }
    }
    return fails;
}

__attribute__((noinline))
u64 v15r_action(const u64 *rep, const u64 *win) {
#ifdef MUT_DETACH_ACTION
    /* v15q behaviour kept as the control: the counter moves, nothing acts. */
    (void)rep; (void)win;
    return 0;
#else
    return (u64)((*rep >= REP_THRESH) || (*win >= WINDOW));
#endif
}

__attribute__((noinline))
void v15r_reset_hist(u64 *prev0, u64 *prev1, u64 *rep, u64 *minx, u64 *win) {
    for (int i = 0; i < 4; i++) { prev0[i] = 0; prev1[i] = 0; minx[i] = ~0ULL; }
    *rep = 0; *win = 0;
}

/* ---- v15r_escape: the restart routine the detector's counters invoke -----
 * Fresh coins are a random 16-bit offset on each label plus the matching point
 * action, so the escaped trail lands on a point it has not visited and its
 * history is genuinely cleared (v15r_reset_hist).  The caller then continues
 * from the returned state with the updated (a, c).
 */
static inline int small_mul(u64 ox[4], u64 oy[4], u64 k, const u64 *px, const u64 *py) {
    u64 rx[4], ry[4], ax[4], ay[4];
    int have = 0;
    for (int i = 15; i >= 0; i--) {
        if (have) {
            if (v15q_add(ax, ay, rx, ry, rx, ry) != 0) return -1;
            for (int j = 0; j < 4; j++) { rx[j]=ax[j]; ry[j]=ay[j]; }
        }
        if ((k >> i) & 1) {
            if (!have) { for (int j = 0; j < 4; j++) { rx[j]=px[j]; ry[j]=py[j]; } have = 1; }
            else {
                if (v15q_add(ax, ay, rx, ry, px, py) != 0) return -1;
                for (int j = 0; j < 4; j++) { rx[j]=ax[j]; ry[j]=ay[j]; }
            }
        }
    }
    if (!have) { for (int j = 0; j < 4; j++) { ox[j]=0; oy[j]=0; } return 1; }
    for (int j = 0; j < 4; j++) { ox[j]=rx[j]; oy[j]=ry[j]; }
    return 0;
}

__attribute__((noinline))
int v15r_escape(u64 x[4], u64 y[4], u64 a[4], u64 c[4],
                const u64 *p0x, const u64 *p0y, const u64 *qx, const u64 *qy,
                u64 s_off, u64 t_off,
                u64 *prev0, u64 *prev1, u64 *rep, u64 *minx, u64 *win) {
    u64 sx[4], sy[4], tx[4], ty[4], nx[4], ny[4], tmp[4];
    if (small_mul(sx, sy, s_off, p0x, p0y) < 0) return -1;
    if (small_mul(tx, ty, t_off, qx, qy) < 0) return -1;
    u64 mx[4], my[4];
    if (v15q_add(nx, ny, x, y, sx, sy) != 0) return -1;
    /* v15q_add must not alias its output with an input (x1 is read after x3) */
    if (v15q_add(mx, my, nx, ny, tx, ty) != 0) return -1;
    for (int i = 0; i < 4; i++) { nx[i] = mx[i]; ny[i] = my[i]; }
    for (int i = 0; i < 4; i++) { x[i] = nx[i]; y[i] = ny[i]; }
    for (int i = 0; i < 4; i++) tmp[i] = (i == 0) ? s_off : 0;
    rmod_add(a, tmp);
    for (int i = 0; i < 4; i++) tmp[i] = (i == 0) ? t_off : 0;
    rmod_add(c, tmp);
    v15r_reset_hist(prev0, prev1, rep, minx, win);
    return 0;
}
