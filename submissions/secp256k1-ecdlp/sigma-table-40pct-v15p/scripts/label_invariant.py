#!/usr/bin/env python3
"""label_invariant.py -- execute the SHIPPED bookkeeping and check its invariant.

Builds scripts/walk_v15o.c (plus a driver TU defining the field constants) into
a dylib and drives v15o_align / v15o_label through ctypes on the real
secp256k1 curve, with independent big-integer arithmetic as the oracle.

  INVARIANT  2000-step walk: at every step the canonical key the C wrote equals
             min(beta^k * Araw.x) computed independently and the stored y is the
             even representative; at the end  A_end == [a]G + [c]Q  for the
             (a,c) the C carried.  This is the recovery invariant the v15o
             verdict says the label bookkeeping fails.

  PARTNER    the step map must be class-equivariant: for 1000 random points P
             and each orbit partner (phi P, phi^2 P, -P), step(P) and
             step(partner) must land in the SAME class.  3000 comparisons.

  CONTROLS   six mutations of the same source, each re-introducing exactly one
             v15o defect: the four label defects must FAIL the invariant, and
             the input-dependent-address mutation must split the partner walk
             (that is the reviewer's argmin-slot finding, reproduced).

No production execution is claimed: 2000 steps, one target, secp256k1 itself.
"""
import ctypes
import json
import os
import random
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
BUILD = "/Users/andrewgordon/.hermes/cache/scratch/v15p-build"
os.makedirs(BUILD, exist_ok=True)

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE
LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
N = 2000
NPTS = 1000
RAD = 2 ** 256 % r
LAM_INV = [1, LAM * LAM % r, LAM]          # lambda^-u for u = 0,1,2


def limbs(v):
    return [(v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF for i in range(4)]


def mont(v):
    return v * RAD % r


def inv(a, m=p):
    return pow(a, m - 2, m)


def padd(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    lam = (3 * x1 * x1) * inv(2 * y1) % p if P == Q else (y2 - y1) * inv(x2 - x1) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def pneg(P):
    return None if P is None else (P[0], (-P[1]) % p)


def pmul(k, P):
    R, Q = None, P
    while k:
        if k & 1:
            R = padd(R, Q)
        Q = padd(Q, Q)
        k >>= 1
    return R


def phi(P, j=1):
    return None if P is None else (pow(BETA, j, p) * P[0] % p, P[1])


G = (GX, GY)
assert pmul(r, G) is None and phi(G) == pmul(LAM, G)

DRIVER = r"""
typedef unsigned long long u64;
extern const u64 BETA[4];
u64 NODE_HEAD[12 * 1024];
u64 NODE_MASK = 1023;
u64 R_CONST[4], R_N0, LAMBDA_R[4], LAMBDA2_R[4], R_MONT_ONE[4];
const u64 BETA[4] = {%s};
""" % ",".join("0x%016xULL" % v for v in limbs(BETA))


def build(name, defines=()):
    dylib = os.path.join(BUILD, "lib%s.dylib" % name)
    cfile = os.path.join(BUILD, "driver%s.c" % name)
    open(cfile, "w").write(DRIVER)
    cmd = ["clang", "-O2", "-arch", "arm64", "-dynamiclib", "-I", HERE, cfile,
           os.path.join(HERE, "walk_v15o.c"), "-o", dylib]
    for d in defines:
        cmd.insert(1, "-D" + d)
    subprocess.run(cmd, check=True)
    lib = ctypes.CDLL(dylib)
    U64x4 = ctypes.c_uint64 * 4

    class Rec(ctypes.Structure):
        _fields_ = [("x", U64x4), ("y", U64x4), ("idx", ctypes.c_uint64),
                    ("len", ctypes.c_uint64), ("esc", ctypes.c_uint64)]

    lib.v15o_align.restype = ctypes.c_uint64
    lib.v15o_align.argtypes = [ctypes.POINTER(Rec), ctypes.POINTER(ctypes.c_int)]
    for nm, val in (("R_CONST", r), ("LAMBDA_R", mont(LAM)),
                    ("LAMBDA2_R", mont(LAM * LAM % r)), ("R_MONT_ONE", RAD)):
        arr = U64x4.in_dll(lib, nm)
        for i, l in enumerate(limbs(val)):
            arr[i] = l
    ctypes.c_uint64.in_dll(lib, "R_N0").value = (-pow(r, -1, 2 ** 64)) % 2 ** 64
    return lib, Rec


def align(lib, Rec, Pt):
    rec = Rec()
    for i, l in enumerate(limbs(Pt[0])):
        rec.x[i] = l
    for i, l in enumerate(limbs(Pt[1])):
        rec.y[i] = l
    sgn = ctypes.c_int(0)
    slot = lib.v15o_align(ctypes.byref(rec), ctypes.byref(sgn))
    gx = sum(rec.x[i] << (64 * i) for i in range(4))
    gy = sum(rec.y[i] << (64 * i) for i in range(4))
    return gx, gy, slot, sgn.value, rec.idx


def label(lib, a, c, slot, sgn, u, v):
    U64x4 = ctypes.c_uint64 * 4
    A, C, JU, JV = U64x4(), U64x4(), U64x4(), U64x4()
    for i in range(4):
        A[i] = (a >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
        C[i] = (c >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
        JU[i] = (u >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
        JV[i] = (v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
    lib.v15o_label(A, C, ctypes.c_uint64(slot), ctypes.c_int(sgn), JU, JV)
    return (sum(A[i] << (64 * i) for i in range(4)),
            sum(C[i] << (64 * i) for i in range(4)))


# orbit-generated rows: row[k][u] = phi^-u J_k, coefficients scaled by lambda^-u
K0 = 0x5EC0FFEE1234567890ABCDEF0123456789ABCDEF0123456789ABCDEF01234567 % r
Qs = pmul(K0, G)
ROWS = 2048
ROW = []
for k in range(ROWS):
    uk = (0x1234567 * k + 0xDEADBEEF) % r
    vk = (0x89ABCDE * k + 0xFEEDFACE) % r
    Jk = padd(pmul(uk, G), pmul(vk, Qs))
    for u in range(3):
        sc = LAM_INV[u]
        ROW.append((sc * uk % r, sc * vk % r, phi(Jk, (3 - u) % 3)))
assert all(len(x) == 3 for x in ROW) and len(ROW) == 3 * ROWS


def row_of(idx):
    return ROW[idx % (3 * ROWS)]


def step(P, idx):
    """one shipped-coded step: A' = canon(A + row(idx))"""
    u, v, J = row_of(idx)
    return padd(P, J), u, v


out = {"kind": "label-invariant-v15p", "curve": "secp256k1", "steps": N,
       "partner_points": NPTS, "rows": []}

lib, Rec = build("ship")
rnd = random.Random(20261024)


def invariant(lib, Rec, n=N):
    start = pmul(0x1234567890ABCDEF % r, G)
    gx, gy, slot, sgn, idx = align(lib, Rec, start)
    P, a, c = (gx, gy), 0, 0
    bad = 0
    for _ in range(n):
        u, v, J = row_of(idx)
        Araw = padd(P, J)
        gx, gy, slot, sgn, idx = align(lib, Rec, Araw)
        cand = [pow(BETA, k, p) * Araw[0] % p for k in range(3)]
        if gx != min(cand):
            bad += 1
        ey = Araw[1] if Araw[1] % 2 == 0 else (-Araw[1]) % p
        if gy != ey:
            bad += 1
        a, c = label(lib, a, c, slot, sgn, u, v)
        P = (gx, gy)
    ok = P == padd(pmul(a, G), pmul(c, Qs))
    return bad, bool(ok), P, a, c


bad, ok, P, a, c = invariant(lib, Rec)
out["rows"].append({"row": "fixed_shipped", "invariant_violations": bad,
                    "end_point_equals_aG_plus_cQ": ok,
                    "verdict": "PASS" if bad == 0 and ok else "FAIL"})


def partner(lib, Rec, npts=NPTS):
    mism = 0
    tot = 0
    for i in range(npts):
        P = pmul(rnd.randrange(1, r), G)
        _, _, _, _, i0 = align(lib, Rec, P)
        S0, _, _ = step(P, i0)
        k0, _, _, _, _ = align(lib, Rec, S0)
        for Pt in (phi(P, 1), phi(P, 2), pneg(P)):
            _, _, _, _, i1 = align(lib, Rec, Pt)
            S1, _, _ = step(Pt, i1)
            k1, _, _, _, _ = align(lib, Rec, S1)
            tot += 1
            if k0 != k1:
                mism += 1
    return mism, tot


m, t = partner(lib, Rec)
out["partner_rows"] = [{"row": "fixed_shipped", "comparisons": t, "class_mismatches": m,
                        "verdict": "PASS" if m == 0 else "FAIL"}]

MUTS = [("control_mont_one_is_integer", ["MUT_MONT_ONE"], "invariant"),
        ("control_negate_on_borrow", ["MUT_NEG_BORROW"], "invariant"),
        ("control_jump_not_added", ["MUT_NO_JUMP"], "invariant"),
        ("control_sign_not_consumed", ["MUT_NO_SIGN"], "invariant"),
        ("control_index_from_raw_x", ["MUT_IDX_RAW"], "partner"),
        ("control_index_v15o_slot_only", ["MUT_IDX_SLOT"], "partner")]
for name, defs, kind in MUTS:
    lb, Rc = build(name, defs)
    if kind == "invariant":
        b, o, _, _, _ = invariant(lb, Rc, n=300)
        out["rows"].append({"row": name, "kind": "invariant",
                            "invariant_violations": b, "end_point_equals_aG_plus_cQ": o,
                            "steps": 300,
                            "verdict": "PASS" if b == 0 and o else "FAIL"})
    else:
        mm, tt = partner(lb, Rc, npts=200)
        out["partner_rows"].append({"row": name, "comparisons": tt,
                                    "class_mismatches": mm,
                                    "verdict": "PASS" if mm == 0 else "FAIL"})

out["verdict"] = ("PASS" if out["rows"][0]["verdict"] == "PASS"
                  and out["partner_rows"][0]["verdict"] == "PASS"
                  and all(x["verdict"] == "FAIL" for x in out["rows"][1:])
                  and all(x["verdict"] == "FAIL" for x in out["partner_rows"][1:])
                  else "FAIL")
os.makedirs(os.path.join(PKG, "evidence"), exist_ok=True)
json.dump(out, open(os.path.join(PKG, "evidence", "label-invariant-v15p.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
