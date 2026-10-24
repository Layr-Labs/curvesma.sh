/* align_step_v15n.c -- the SHIPPED alignment routine of the six-class rho step,
 * plus the SHIPPED per-report processor.  Both are compiled with
 * `clang -O2 -arch arm64 -S` and counted from the emitted listing
 * (scripts/count_iter.py, evidence/counts-step-v15n.json); the equivariance and
 * attribution probes execute THIS object through ctypes (scripts/probe_v15n.py).
 *
 * Why the alignment is separate from the 748-Ir point listing: the base listing
 * (evidence/step_secp_v15g.s) is unchanged and priced at 748 Ir.  The v15m defect
 * was an alignment *block appended inside that loop* behind `tbz w0,#0,.Ltail`,
 * which (a) skipped alignment on even-y outputs, (b) copied the computed y into
 * the registers that .Ltail stores as y, and (c) used a 32-byte stride against a
 * 64-byte affine table.  Here the alignment is a function of the OUTPUT x only
 * and is called once per step on EVERY output -- it cannot be bypassed, it never
 * sees y as a multiplicand, and the table row it selects is the same 64-byte
 * affine entry the point listing loads.
 *
 * Curve: secp256k1, p = 2^256 - C, C = 0x1000003D1.  Four 64-bit limbs, LE.
 */
typedef unsigned long long u64;
typedef unsigned __int128 u128;

extern const u64 BETA[4];      /* cube root of unity mod p   (driver TU) */

typedef struct { u64 x[4], y[4], slot, len, esc, pad; } rec_t;   /* 96 B */
extern rec_t A[];              /* walk records               (driver TU) */
extern u64 NREC;               /* records in this sweep      (driver TU) */

#define C_LIMB 0x1000003D1ULL
static const u64 P0 = 0xFFFFFFFEFFFFFC2FULL;
static const u64 P1 = 0xFFFFFFFFFFFFFFFFULL;

static inline int fcmp(const u64 a[4], const u64 b[4]) {
    for (int i = 3; i >= 0; i--) {
        if (a[i] < b[i]) return -1;
        if (a[i] > b[i]) return 1;
    }
    return 0;
}

static inline void fsub(u64 r[4], const u64 a[4], const u64 b[4]) {
    u64 t[4]; u128 br = 0;
    br = (u128)a[0] - b[0]; t[0] = (u64)br; br = (br >> 64) & 1;
    br = (u128)a[1] - b[1] - br; t[1] = (u64)br; br = (br >> 64) & 1;
    br = (u128)a[2] - b[2] - br; t[2] = (u64)br; br = (br >> 64) & 1;
    br = (u128)a[3] - b[3] - br; t[3] = (u64)br;
    if (br >> 64) {
        u128 c = 0;
        c += (u128)t[0] + P0; r[0] = (u64)c; c >>= 64;
        c += (u128)t[1] + P1; r[1] = (u64)c; c >>= 64;
        c += (u128)t[2] + P1; r[2] = (u64)c; c >>= 64;
        c += (u128)t[3] + P1; r[3] = (u64)c;
    } else { r[0]=t[0]; r[1]=t[1]; r[2]=t[2]; r[3]=t[3]; }
}

static inline void fmul(u64 r[4], const u64 a[4], const u64 b[4]) {
    u64 t[8] = {0,0,0,0,0,0,0,0};
    for (int i = 0; i < 4; i++) {
        u64 cy = 0;
        for (int j = 0; j < 4; j++) {
            u128 acc = (u128)a[i] * b[j] + t[i+j] + cy;
            t[i+j] = (u64)acc; cy = (u64)(acc >> 64);
        }
        t[i+4] = cy;
    }
    u64 s[6] = {t[0], t[1], t[2], t[3], 0, 0};
    for (int k = 7; k >= 4; k--) {
        u128 m = (u128)t[k] * C_LIMB;
        int q0 = k - 4;
        u128 c = (u128)s[q0] + (u64)m; s[q0] = (u64)c; c >>= 64;
        for (int q = q0 + 1; q < 6; q++) {
            c += (u128)s[q] + (q == q0 + 1 ? (u64)(m >> 64) : 0);
            s[q] = (u64)c; c >>= 64;
        }
        if (c) { u64 z[4] = {0,0,0,0}; (void)z; s[5] += (u64)c; }
    }
    while (s[4] | s[5]) {
        u128 m = (u128)s[4] * C_LIMB;
        u64 hi = s[5];
        s[4] = 0; s[5] = 0;
        u128 c = (u128)s[0] + (u64)m; s[0] = (u64)c; c >>= 64;
        for (int q = 1; q < 4; q++) { c += s[q]; s[q] = (u64)c; c >>= 64; }
        c += (u128)hi * C_LIMB + (u64)(m >> 64);
        u128 c2 = (u128)s[0] + (u64)c; s[0] = (u64)c2; c2 >>= 64;
        for (int q = 1; q < 4; q++) { c2 += s[q]; s[q] = (u64)c2; c2 >>= 64; }
    }
    r[0]=s[0]; r[1]=s[1]; r[2]=s[2]; r[3]=s[3];
    for (int k = 0; k < 2; k++) {
        u64 q[4]; u128 br = 0;
        br = (u128)r[0] - P0; q[0] = (u64)br; br = (br >> 64) & 1;
        br = (u128)r[1] - P1 - br; q[1] = (u64)br; br = (br >> 64) & 1;
        br = (u128)r[2] - P1 - br; q[2] = (u64)br; br = (br >> 64) & 1;
        br = (u128)r[3] - P1 - br; q[3] = (u64)br;
        if (!(br >> 64)) { r[0]=q[0]; r[1]=q[1]; r[2]=q[2]; r[3]=q[3]; }
    }
}

/* (key, slot): key = min(x, BETA*x, BETA^2*x); slot = argmin in {0,1,2}.
 * Orbit-invariant by construction; called on EVERY output x of the step. */
void class_slot(const u64 x[4], u64 key[4], u64 *slot) {
    u64 bx[4], b2x[4], m[4];
    u64 im = 0;
    fmul(bx, x, BETA);
    fmul(b2x, bx, BETA);
    for (int i = 0; i < 4; i++) m[i] = x[i];
    if (fcmp(bx, m) < 0)  { for (int i=0;i<4;i++) m[i]=bx[i];  im = 1; }
    if (fcmp(b2x, m) < 0) { for (int i=0;i<4;i++) m[i]=b2x[i]; im = 2; }
    for (int i = 0; i < 4; i++) key[i] = m[i];
    *slot = im;
}

/* even-y sign normalisation, applied on the same every-output path */
void even_y(u64 y[4]) {
    u64 z[4] = {0,0,0,0}, ny[4];
    if (y[0] & 1ULL) { fsub(ny, z, y); for (int i = 0; i < 4; i++) y[i] = ny[i]; }
}

/* align_sweep: the priced alignment pass, one call per walk step over NREC
 * records.  Writes the 64-byte table row index 3*(key[0]&2047)+slot and the
 * sign-normalised y back into the state, and accumulates len. */
void align_sweep(void) {
    u64 key[4], slot;
    for (u64 i = 0; i < NREC; i++) {
        rec_t *r = &A[i];
        class_slot(r->x, key, &slot);                  /* ALWAYS: no parity test */
        r->pad = 3ULL * (key[0] & 2047ULL) + slot;     /* 64-byte affine row     */
        even_y(r->y);                                  /* sign fold, all outputs */
        for (int q = 0; q < 4; q++) r->x[q] = key[q];
        r->slot = slot;
        r->len += 1;
    }
}

/* ---- per-report processor: hash / route / compare / store -------------- */
/* One call per emitted report.  ~S_ops/2 reports per run at the dp1 band, so
 * this row is charged at (reports) x (this count) / 1000 units. */
extern u64 NODE_HEAD[65536];      /* 2^16 routing nodes, 4-limb key slot */
extern u64 REPORT_STORE[8];       /* the node's stored 256-bit key       */
extern u64 REPORT_REC[4];         /* payload: class key of the report    */

void report_proc(u64 node_out[1]) {
    u64 h[4], key[4], slot;
    for (int i = 0; i < 4; i++) key[i] = REPORT_REC[i];
    fmul(h, key, BETA);                       /* universal-ish mix: 2 products */
    fmul(h, h, key);
    u64 node = h[0] & 0xFFFFULL;              /* route to a 2^16 partition    */
    node_out[0] = node;
    u64 *head = &NODE_HEAD[node];
    if (fcmp(key, head) == 0) {               /* compare against the stored key */
        REPORT_STORE[0] = 1;                  /* duplicate: mark, no store    */
        return;
    }
    for (int i = 0; i < 4; i++) head[i] = key[i];   /* store the new key      */
    REPORT_STORE[0] = 0;
}

/* exported wrappers so the probes exercise the SAME translation unit */
void api_fmul(u64 r[4], const u64 a[4], const u64 b[4]) { fmul(r, a, b); }
void api_fsub(u64 r[4], const u64 a[4], const u64 b[4]) { fsub(r, a, b); }
void api_class(const u64 x[4], u64 key[4], u64 *slot) { class_slot(x, key, slot); }
