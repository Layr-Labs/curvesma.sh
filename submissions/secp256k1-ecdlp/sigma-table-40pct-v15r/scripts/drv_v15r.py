#!/usr/bin/env python3
"""drv_v15q.py -- build the SHIPPED bookkeeping (scripts/walk_v15r.c) into a
dylib for a chosen curve, plus the one function the C source does not contain:
the affine point step.

Both the real-curve harness and the census fixture compile THE SAME source with
different CONSTANTS only (p, beta, order r, lambda, curve b).  No line of the
canonicalisation, the index rule, the label, the detector or the pool is
re-implemented in Python: Python supplies the table points, the walk order and
the big-integer ORACLE, nothing else.

The point step is the priced 6-product affine step (evidence/step_secp_v15c.c
body: d = x_J - x, lam = (y_J - y)/d, x3 = lam^2 - x - x_J, y3 = lam(x-x3) - y)
with the inversion by Fermat.  The listing's PROD/SUF prefix-suffix batch-
inversion schedule is priced as its own ledger row; it changes how d^-1 is
obtained, not the transition this harness certifies.
"""
import ctypes
import os
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = "/Users/andrewgordon/.hermes/cache/scratch/v15q-build"
os.makedirs(BUILD, exist_ok=True)


def limbs(v):
    return [(v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(4)]


def c4(v):
    return "{" + ",".join("0x%016xULL" % l for l in limbs(v)) + "}"


DRIVER = r"""
typedef unsigned long long u64;
typedef unsigned __int128 u128;
#define V15Q_P0 %(P0)s
#define V15Q_C  %(C)s
extern const u64 BETA[4];
u64 NODE_HEAD[12 * %(NODES)d];
u64 NODE_MASK = %(NODES)d - 1;
u64 R_CONST[4] = %(R)s;
u64 R_N0 = 0;
u64 LAMBDA_R[4] = %(LAM)s;
u64 LAMBDA2_R[4] = %(LAM2)s;
u64 R_MONT_ONE[4] = %(MONT1)s;
u64 R2_R[4] = %(R2)s;
const u64 BETA[4] = %(BETA)s;
const u64 P_EXP[4] = %(PEXP)s;          /* p - 2 */
const u64 R_EXP[4] = %(REXP)s;          /* r - 2 */
u64 TBL[8 * 3 * 2048];
u64 UBASE[8 * 2048];
u64 SMALL[8 * 16];

#include "walk_v15r.c"
"""

# The secp fsub needs the real borrow chain; the fixture one is a single limb.
FSUB_SMALL = r"""
static inline void fsub(u64 o[4], const u64 a[4], const u64 b[4]) {
    o[0] = (a[0] + V15Q_P0 - (b[0] % V15Q_P0)) % V15Q_P0;
    o[1] = o[2] = o[3] = 0;
}
"""

FSUB_SECP = r"""
static inline void fsub(u64 o[4], const u64 a[4], const u64 b[4]) {
    u64 d[4]; u64 br = 0;
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) { u128 w = (u128)a[i] - b[i] - br; d[i] = (u64)w; br = (u64)(w >> 64) & 1; }
    u64 s[4]; u64 cy = 0;
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) { u128 w = (u128)d[i] + (i == 0 ? V15Q_P0 : ~0ULL) + cy; s[i] = (u64)w; cy = (u64)(w >> 64) & 1; }
    if (br) { o[0]=s[0]; o[1]=s[1]; o[2]=s[2]; o[3]=s[3]; }
    else    { o[0]=d[0]; o[1]=d[1]; o[2]=d[2]; o[3]=d[3]; }
}
"""


def build(name, p, r, beta, lam, nodes, defines=()):
    """Compile walk_v15r.c + driver for the given curve.  Returns (lib, Rec)."""
    small = p.bit_length() <= 32
    src = DRIVER % {
        "P0": ("0x%016xULL" % (p & 0xFFFFFFFFFFFFFFFF)) if not small else ("0x%xULL" % p),
        "C": "0x1ULL" if small else "0x1000003D1ULL",
        "NODES": nodes,
        "R": c4(r), "LAM": c4(lam * (1 << 256) % r), "LAM2": c4(lam * lam % r * (1 << 256) % r),
        "MONT1": c4((1 << 256) % r), "R2": c4(pow(2, 512, r)),
        "BETA": c4(beta), "PEXP": c4(p - 2), "REXP": c4(r - 2),
    }
    body = r"""
/* exponentiation in the field (Fermat inverse) and in Z_r */
static inline void fpow(u64 o[4], const u64 a[4], const u64 *e,
                        void (*mul)(u64*, const u64*, const u64*)) {
    u64 t[4] = {1,0,0,0};
    int top = 255;
    while (top >= 0 && !((e[top >> 6] >> (top & 63)) & 1)) top--;
    for (int i = top; i >= 0; i--) {
        u64 s[4]; mul(s, t, t); t[0]=s[0]; t[1]=s[1]; t[2]=s[2]; t[3]=s[3];
        if ((e[i >> 6] >> (i & 63)) & 1) {
            u64 q[4]; mul(q, t, a); t[0]=q[0]; t[1]=q[1]; t[2]=q[2]; t[3]=q[3];
        }
    }
    o[0]=t[0]; o[1]=t[1]; o[2]=t[2]; o[3]=t[3];
}
static inline void fmul_w(u64 o[4], const u64 a[4], const u64 b[4]) { fmul(o, a, b); }
static inline void rmont_w(u64 o[4], const u64 a[4], const u64 b[4]) { rmont(o, a, b); }
static inline void finv(u64 o[4], const u64 a[4]) { fpow(o, a, P_EXP, fmul_w); }
/* Montgomery-domain power: start from R mod r, every step rmont, canonicalise out.
 * v15q started from the INTEGER 1 under rmont, which returns a^k*R^-(k-1), so its
 * inverse -- and therefore its recovery solve -- was silently wrong (measured:
 * 20/20 wrong x with rc==0). */
static inline void fpow_mont(u64 o[4], const u64 a[4], const u64 *e) {
    u64 t[4]; for (int i = 0; i < 4; i++) t[i] = R_MONT_ONE[i];
    int top = 255;
    while (top >= 0 && !((e[top >> 6] >> (top & 63)) & 1)) top--;
    for (int i = top; i >= 0; i--) {
        u64 s[4]; rmont(s, t, t);
        for (int j = 0; j < 4; j++) t[j] = s[j];
        if ((e[i >> 6] >> (i & 63)) & 1) {
            u64 q[4]; rmont(q, t, a);
            for (int j = 0; j < 4; j++) t[j] = q[j];
        }
    }
    for (int i = 0; i < 4; i++) o[i] = t[i];
}
/* canonical product xy mod r via two Montgomery products */
static inline void rmul(u64 o[4], const u64 a[4], const u64 b[4]) {
    u64 t[4]; rmont(t, a, b); rmont(o, t, R2_R);
}
static inline void rinv(u64 o[4], const u64 a[4]) {
    /* the exponentiation runs on the MONTGOMERY REPRESENTATION aR: rmont(t,a)
     * with canonical a leaves the domain (a^e R * a * R^-1 = a^(e+1)), which is
     * what made v15q's inverse wrong for every input. */
    u64 aR[4]; rmont(aR, a, R2_R);
    u64 t[4]; fpow_mont(t, aR, R_EXP);
    u64 one[4] = {1,0,0,0};
    rmont(o, t, one);
}

/* ---- the priced affine point step --------------------------------------- */
__attribute__((noinline))
int v15q_step(u64 x[4], u64 y[4], u64 idx) {
    const u64 *jx = TBL + 8 * (idx % (3 * 2048));
    const u64 *jy = jx + 4;
    u64 d[4], num[4], inv[4], lam[4], t[4], x3[4], y3[4], e[4];
    fsub(d, jx, x);
    fsub(num, jy, y);
    if ((d[0] | d[1] | d[2] | d[3]) == 0) return -1;
    finv(inv, d);
    fmul(lam, num, inv);
    fmul(t, lam, lam);
    fsub(x3, t, x); fsub(x3, x3, jx);
    fsub(e, x, x3);
    fmul(y3, lam, e); fsub(y3, y3, y);
    x[0]=x3[0]; x[1]=x3[1]; x[2]=x3[2]; x[3]=x3[3];
    y[0]=y3[0]; y[1]=y3[1]; y[2]=y3[2]; y[3]=y3[3];
    return 0;
}

/* general affine add, used by the oracle-side scalar multiplication only */
__attribute__((noinline))
int v15q_add(u64 x3[4], u64 y3[4], const u64 x1[4], const u64 y1[4],
             const u64 x2[4], const u64 y2[4]) {
    static const u64 TWO[4] = {2,0,0,0};
    static const u64 THREE[4] = {3,0,0,0};
    u64 d[4], num[4], inv[4], lam[4], t[4], e[4];
    int same = (x1[0]^x2[0] | x1[1]^x2[1] | x1[2]^x2[2] | x1[3]^x2[3]) == 0;
    if (same) {
        if ((y1[0]^y2[0] | y1[1]^y2[1] | y1[2]^y2[2] | y1[3]^y2[3]) != 0) return -1;
        fmul(t, x1, x1); fmul(num, t, THREE);
        fmul(d, y1, TWO);
    } else {
        fsub(d, x2, x1); fsub(num, y2, y1);
    }
    if ((d[0] | d[1] | d[2] | d[3]) == 0) return -1;
    finv(inv, d);
    fmul(lam, num, inv);
    fmul(t, lam, lam);
    fsub(x3, t, x1); fsub(x3, x3, x2);
    fsub(e, x1, x3);
    fmul(y3, lam, e); fsub(y3, y3, y1);
    return 0;
}

/* solve x from two records of one class: x = -(a1-a2)/(c1-c2) mod r */
static inline void rsub(u64 o[4], const u64 a[4], const u64 b[4]) {
    u64 br = 0;
    #pragma clang loop unroll(full)
    for (int i = 0; i < 4; i++) { u128 w = (u128)a[i] - b[i] - br; o[i] = (u64)w; br = (u64)(w >> 64) & 1; }
    if (br) {
        u64 cy = 0;
        #pragma clang loop unroll(full)
        for (int i = 0; i < 4; i++) { u128 w = (u128)o[i] + R_CONST[i] + cy; o[i] = (u64)w; cy = (u64)(w >> 64) & 1; }
    }
}

__attribute__((noinline))
int v15q_solve(const u64 a1[4], const u64 c1[4], const u64 a2[4], const u64 c2[4], u64 out[4]) {
    u64 da[4], dc[4], inv[4], t[4];
    rsub(da, a1, a2);
    rsub(dc, c1, c2);
    if ((dc[0] | dc[1] | dc[2] | dc[3]) == 0) return -1;
    rinv(inv, dc);
    rmul(t, da, inv);
    if ((t[0] | t[1] | t[2] | t[3]) == 0) return -1;
    rneg_unit(t);
    out[0]=t[0]; out[1]=t[1]; out[2]=t[2]; out[3]=t[3];
    return 0;
}

__attribute__((noinline))
void v15r_dbg(u64 out[4], const u64 a[4], const u64 b[4]) { rmul(out, a, b); }
__attribute__((noinline))
void v15r_dbg_inv(u64 out[4], const u64 a[4]) { rinv(out, a); }
/* [k]P by double-and-add (oracle side only) */
__attribute__((noinline))
int v15q_mul(u64 x[4], u64 y[4], const u64 k[4], const u64 px[4], const u64 py[4]) {
    u64 rx[4]={0,0,0,0}, ry[4]={0,0,0,0}, ax[4], ay[4];
    int have = 0;
    for (int i = 255; i >= 0; i--) {
        if (have) {
            if (v15q_add(ax, ay, rx, ry, rx, ry) != 0) return -1;
            for (int j = 0; j < 4; j++) { rx[j]=ax[j]; ry[j]=ay[j]; }
        }
        if ((k[i >> 6] >> (i & 63)) & 1) {
            if (!have) {
                for (int j = 0; j < 4; j++) { rx[j]=px[j]; ry[j]=py[j]; }
                have = 1;
            } else {
                if (v15q_add(ax, ay, rx, ry, px, py) != 0) return -1;
                for (int j = 0; j < 4; j++) { rx[j]=ax[j]; ry[j]=ay[j]; }
            }
        }
    }
    if (!have) return -1;
    x[0]=rx[0];x[1]=rx[1];x[2]=rx[2];x[3]=rx[3];
    y[0]=ry[0];y[1]=ry[1];y[2]=ry[2];y[3]=ry[3];
    return 0;
}
"""
    src = src + (FSUB_SMALL if small else FSUB_SECP) + body
    cfile = os.path.join(BUILD, "drv_%s.c" % name)
    open(cfile, "w").write(src)
    dylib = os.path.join(BUILD, "lib%s.dylib" % name)
    cmd = ["clang", "-O2", "-arch", "arm64", "-dynamiclib", "-I", HERE, cfile, "-o", dylib]
    if small:
        cmd.insert(1, "-DV15Q_SMALL_FIELD")
    for d in defines:
        cmd.insert(1, "-D" + d)
    subprocess.run(cmd, check=True)
    lib = ctypes.CDLL(dylib)
    U64x4 = ctypes.c_uint64 * 4

    class Rec(ctypes.Structure):
        _fields_ = [("x", U64x4), ("y", U64x4), ("idx", ctypes.c_uint64),
                    ("len", ctypes.c_uint64), ("esc", ctypes.c_uint64)]

    for nm in ("v15o_align", "v15o_label", "v15o_detect", "v15o_report",
               "v15q_step", "v15q_add", "v15q_solve", "v15q_mul"):
        getattr(lib, nm).restype = ctypes.c_int if nm.startswith("v15q") else ctypes.c_uint64
    lib.v15r_init_for_Q.restype = ctypes.c_uint64
    lib.v15r_action.restype = ctypes.c_uint64
    lib.v15r_action.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(ctypes.c_uint64)]
    lib.v15r_reset_hist.argtypes = [ctypes.c_void_p] * 5
    lib.v15r_dbg.argtypes = [U64x4, U64x4, U64x4]
    lib.v15r_dbg_inv.argtypes = [U64x4, U64x4]
    lib.v15r_escape.restype = ctypes.c_int
    lib.v15r_escape.argtypes = [U64x4, U64x4, U64x4, U64x4, U64x4, U64x4, U64x4, U64x4,
                                ctypes.c_uint64, ctypes.c_uint64] + [ctypes.c_void_p] * 5
    lib.v15r_init_for_Q.argtypes = [ctypes.POINTER(ctypes.c_uint64), U64x4, U64x4]
    lib.v15o_align.argtypes = [ctypes.POINTER(Rec), ctypes.POINTER(ctypes.c_int)]
    lib.v15o_label.argtypes = [U64x4, U64x4, ctypes.c_uint64, ctypes.c_int, U64x4, U64x4]
    lib.v15o_detect.argtypes = [ctypes.POINTER(Rec)] + [ctypes.c_void_p] * 5
    lib.v15o_report.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.c_uint64,
                                ctypes.POINTER(Rec), U64x4, U64x4, U64x4, U64x4]
    lib.v15q_step.argtypes = [U64x4, U64x4, ctypes.c_uint64]
    lib.v15q_add.argtypes = [U64x4, U64x4, U64x4, U64x4, U64x4, U64x4]
    lib.v15q_solve.argtypes = [U64x4, U64x4, U64x4, U64x4, U64x4]
    lib.v15q_mul.argtypes = [U64x4, U64x4, U64x4, U64x4, U64x4]
    ctypes.c_uint64.in_dll(lib, "R_N0").value = (-pow(r, -1, 2 ** 64)) % 2 ** 64
    return lib, Rec


def set_tbl(lib, pts):
    """TBL is u64[8*3*2048]; fill k -> (x,y) for k < 3*2048 (mod ROWS)."""
    T = (ctypes.c_uint64 * (8 * 3 * 2048)).in_dll(lib, "TBL")
    for k in range(3 * 2048):
        P = pts[k % len(pts)]
        for i in range(4):
            T[8 * k + i] = (P[0] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
            T[8 * k + 4 + i] = (P[1] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
