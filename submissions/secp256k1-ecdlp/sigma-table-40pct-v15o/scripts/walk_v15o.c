/* walk_v15o.c -- the SHIPPED per-step bookkeeping of the six-class sigma walk,
 * v15o.  ONE object carries every row the composition table prices:
 *
 *     v15o_align      canonical (key, idx) selection of the output x  [called 1/step]
 *     v15o_label      exponent (a,c) update under the canonicalisation [called 1/step]
 *     v15o_detect     dp1 predicate + repeat/checkpoint detector       [called 1/step]
 *     v15o_report     open-addressed retained-record pool              [called per report]
 *
 * Disjointness is by SYMBOL: no call of one of these functions contains code of
 * another (all four are leaf functions; the two field products inside v15o_align
 * are the inlined fmul below, counted once, and v15o_label's products are mod r,
 * a different modulus, so no instruction is shared).  The record layout is the
 * one step_secp_v15g.s walks: 80 bytes, x at 0..31, y at 32..63, idx at 64,
 * len at 72 -- so the step's stride and this file's struct are the same object,
 * not two descriptions of one.
 *
 * Curve: secp256k1, p = 2^256 - C, C = 0x1000003D1.  Four 64-bit limbs, LE.
 */
typedef unsigned long long u64;
typedef unsigned __int128 u128;

extern const u64 BETA[4];       /* cube root of unity mod p    (driver TU) */
extern u64 NODE_HEAD[];         /* pool: 12 limbs per node, 96 B (driver TU) */
extern u64 NODE_MASK;           /* nodes - 1                   (driver TU) */
extern u64 R_CONST[4];          /* group order r               (driver TU) */
extern u64 R_N0;                /* -r^-1 mod 2^64              (driver TU) */
extern u64 LAMBDA_R[4];         /* (lambda, lambda^2) in Mont form (driver TU) */

/* the record the point listing walks; 80 bytes, identical in both TUs */
typedef struct { u64 x[4], y[4], idx, len, esc; } rec_t;

#define C_LIMB 0x1000003D1ULL
static const u64 P0 = 0xFFFFFFFEFFFFFC2FULL;

static inline int fcmp(const u64 a[4], const u64 b[4]) {
    for (int i = 3; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

/* one field product mod p = 2^256 - C (fold reduction, secp256k1 form) */
static inline void fmul(u64 r[4], const u64 a[4], const u64 b[4]) {
    u64 t[8] = {0};
    for (int i = 0; i < 4; i++) {
        u128 c = 0;
        for (int j = 0; j < 4; j++) {
            u128 s = (u128)a[i] * b[j] + t[i + j] + c;
            t[i + j] = (u64)s;
            c = s >> 64;
        }
        t[i + 4] = (u64)c;
    }
    /* fold the high half in twice: t = lo + hi*C */
    for (int pass = 0; pass < 2; pass++) {
        u128 c = 0;
        u64 hi[4] = {t[4], t[5], t[6], t[7]};
        for (int i = 0; i < 4; i++) {
            u128 s = (u128)hi[i] * C_LIMB + t[i] + c;
            t[i] = (u64)s;
            c = s >> 64;
        }
        u64 extra = (u64)c;
        t[4] = extra * C_LIMB;          /* extra < 2^34, fits */
        t[5] = t[6] = t[7] = 0;
        if (pass == 0) continue;
        {   /* second pass already folded `extra`; add it once more */
            u128 s = (u128)t[0] + t[4];
            t[0] = (u64)s; u64 cy = (u64)(s >> 64);
            for (int i = 1; i < 4; i++) { u128 q = (u128)t[i] + cy; t[i] = (u64)q; cy = (u64)(q >> 64); }
        }
    }
    u64 z[4] = {t[0], t[1], t[2], t[3]};
    /* conditional subtract p once (result < 2p) */
    u64 d[4]; u64 br = 0;
    u128 w;
    w = (u128)z[0] - P0;        d[0] = (u64)w; br = (u64)(w >> 64) & 1;
    w = (u128)z[1] - ~0ULL - br; d[1] = (u64)w; br = (u64)(w >> 64) & 1;
    w = (u128)z[2] - ~0ULL - br; d[2] = (u64)w; br = (u64)(w >> 64) & 1;
    w = (u128)z[3] - ~0ULL - br; d[3] = (u64)w; br = (u64)(w >> 64) & 1;
    if (!br) { r[0]=d[0]; r[1]=d[1]; r[2]=d[2]; r[3]=d[3]; }
    else     { r[0]=z[0]; r[1]=z[1]; r[2]=z[2]; r[3]=z[3]; }
}

/* ---- v15o_align: canonical class key and table index of the output x ----
 * key = min(x, BETA*x, BETA^2*x) is invariant on the six-element orbit; idx is
 * the argmin slot folded into the table address 3*(key & 2047) + slot, so the
 * NEXT step selects the stored three-slot row and nothing is left to a
 * separate adapter.  Returns the sign bit of y in *sgn. */
__attribute__((noinline))
u64 v15o_align(rec_t *rec, int *sgn) {
    u64 b1[4], b2[4];
    fmul(b1, rec->x, BETA);
    fmul(b2, b1, BETA);
    u64 key[4]; u64 slot = 0;
    const u64 *cand[3] = {rec->x, b1, b2};
    for (int i = 0; i < 3; i++) {
        if (fcmp(cand[i], key) < 0 || i == 0) { key[0]=cand[i][0]; key[1]=cand[i][1]; key[2]=cand[i][2]; key[3]=cand[i][3]; slot = (u64)i; }
    }
    rec->x[0]=key[0]; rec->x[1]=key[1]; rec->x[2]=key[2]; rec->x[3]=key[3];
    rec->idx = 3 * (key[0] & 2047) + slot;
    *sgn = (int)(rec->y[0] & 1);              /* 1 => negate for the even rep */
    if (*sgn) {                                /* y <- p - y */
        u64 d[4]; u64 br = 0; u128 w;
        w = (u128)P0 - rec->y[0];        d[0]=(u64)w; br=(u64)(w>>64)&1;
        w = (u128)~0ULL - rec->y[1] - br; d[1]=(u64)w; br=(u64)(w>>64)&1;
        w = (u128)~0ULL - rec->y[2] - br; d[2]=(u64)w; br=(u64)(w>>64)&1;
        w = (u128)~0ULL - rec->y[3] - br; d[3]=(u64)w;
        rec->y[0]=d[0]; rec->y[1]=d[1]; rec->y[2]=d[2]; rec->y[3]=d[3];
    }
    return slot;
}

/* ---- v15o_label: (a,c) <- s * lambda^slot * (a,c) mod r ----
 * The canonicalisation is a group automorphism [s*lambda^slot]; the exponents
 * that make A = [a]G + [c]Q must be transformed by the SAME scalar or the
 * recovery equation is wrong.  Montgomery multiplication mod r, two products. */
static inline void rmont(u64 out[4], const u64 a[4], const u64 b[4]) {
    u64 t[6] = {0,0,0,0,0,0};
    for (int i = 0; i < 4; i++) {
        u64 cy = 0;
        for (int j = 0; j < 4; j++) {
            u128 s = (u128)a[i] * b[j] + t[j] + cy;
            t[j] = (u64)s; cy = (u64)(s >> 64);
        }
        u64 s2 = t[4] + cy; t[4] = s2;
        u64 m = t[0] * R_N0;
        u128 s3 = (u128)m * R_CONST[0] + t[0];
        u64 c2 = (u64)(s3 >> 64);
        for (int j = 1; j < 4; j++) {
            u128 s4 = (u128)m * R_CONST[j] + t[j] + c2;
            t[j-1] = (u64)s4; c2 = (u64)(s4 >> 64);
        }
        u64 s5 = t[4] + c2; t[3] = s5; t[4] = 0;
    }
    u64 z[4] = {t[0],t[1],t[2],t[3]};
    u64 d[4]; u64 br = 0; u128 w;
    for (int i = 0; i < 4; i++) {
        w = (u128)z[i] - R_CONST[i] - br; d[i] = (u64)w; br = (u64)(w >> 64) & 1;
    }
    if (!br) { out[0]=d[0]; out[1]=d[1]; out[2]=d[2]; out[3]=d[3]; }
    else     { out[0]=z[0]; out[1]=z[1]; out[2]=z[2]; out[3]=z[3]; }
}

__attribute__((noinline))
void v15o_label(u64 a[4], u64 c[4], u64 slot, int sgn) {
    u64 m[4];
    if (slot == 0) { m[0]=1; m[1]=0; m[2]=0; m[3]=0; }   /* Montgomery 1 */
    else if (slot == 1) { m[0]=LAMBDA_R[0]; m[1]=LAMBDA_R[1]; m[2]=LAMBDA_R[2]; m[3]=LAMBDA_R[3]; }
    else { rmont(m, (const u64[4]){LAMBDA_R[0],LAMBDA_R[1],LAMBDA_R[2],LAMBDA_R[3]},
                    (const u64[4]){LAMBDA_R[0],LAMBDA_R[1],LAMBDA_R[2],LAMBDA_R[3]}); }
    rmont(a, a, m);
    rmont(c, c, m);
    if (sgn) {                                            /* negate mod r */
        u64 d[4]; u64 br = 0; u128 w;
        for (int i = 0; i < 4; i++) {
            w = (u128)R_CONST[i] - a[i] - br; d[i] = (u64)w; br = (u64)(w >> 64) & 1;
        }
        if (br) { a[0]=d[0]; a[1]=d[1]; a[2]=d[2]; a[3]=d[3]; }
        for (int i = 0; i < 4; i++) {
            w = (u128)R_CONST[i] - c[i] - (i==0?0:0) - br; d[i] = (u64)w; br = (u64)(w >> 64) & 1;
        }
        if (br) { c[0]=d[0]; c[1]=d[1]; c[2]=d[2]; c[3]=d[3]; }
    }
}

/* ---- v15o_detect: dp1 predicate on the canonical x, plus the repeat detector.
 * dp1 (declared density 2^-1) is the LOW BIT of the canonical x -- not 32 zero
 * bits.  The detector keeps the previous class key and a repeat counter so a
 * two-cycle is seen, not assumed away. */
__attribute__((noinline))
u64 v15o_detect(const rec_t *rec, u64 *prev, u64 *rep) {
    u64 dp = 1 - (rec->x[0] & 1);                  /* 1 => report */
    if (rec->x[0] == prev[0] && rec->x[1] == prev[1] &&
        rec->x[2] == prev[2] && rec->x[3] == prev[3]) *rep += 1;
    prev[0]=rec->x[0]; prev[1]=rec->x[1]; prev[2]=rec->x[2]; prev[3]=rec->x[3];
    return dp;
}

/* ---- v15o_report: retained-record pool, open addressing ----
 * Node = 12 limbs: key 0..3, a 4..7, c 8..11 -- 96 bytes, addressed as
 * NODE_HEAD + 12*node (LIMB stride, so neighbouring nodes never overlap and
 * the top index stays inside the array; the v15n one-limb stride is gone).
 * On a key match the stored exponents are returned: that is the collision the
 * success argument uses.  Returns 1 on a match (pair delivered), 0 on insert. */
__attribute__((noinline))
u64 v15o_report(u64 *pool, u64 mask, const rec_t *rec, const u64 a[4], const u64 c[4],
                u64 out_a[4], u64 out_c[4]) {
    u64 h = (rec->x[0] ^ (rec->x[1] * 0x9E3779B97F4A7C15ULL)) & mask;
    for (u64 probe = 0; probe <= mask; probe++) {
        u64 *node = pool + 12 * ((h + probe) & mask);
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
    return 2;                                      /* pool full: refusal, charged */
}
