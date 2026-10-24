#!/usr/bin/env python3
"""harness_v15q.py -- the SHIPPED bookkeeping executed, on secp256k1 itself.

Every operation below is a call into the C object built from scripts/walk_v15o.c
(v15q): v15o_align (canonical key, sign, slot, index), v15q_step (the priced
affine point step), v15o_label, v15o_detect, v15o_report.  Python supplies the
table points, the walk order, and the big-integer ORACLE -- it never performs a
group operation that the walk is priced for.

ROWS
  invariant   24 walks x 500 steps.  EVERY step is the C step.  Per step the
              stored key must equal min(beta^j * x_raw) (oracle) and the stored
              y must be even; every 25 steps and at the end the stored point
              must equal [a]G + [c]Q for the (a,c) the C label carried, with
              [a]G+cQ computed by the C scalar multiplication.
  partner     120 classes x the SIX orbit representatives (phi, phi^2 and the
              THREE NEGATIVES -P, -phi P, -phi^2 P).  For each representative
              the SHIPPED align is applied first; the test then requires the
              six STORED points to be equal and the six stored indices to be
              equal -- i.e. the state, hence the whole future walk, is a
              function of the class.  Each of the six stored states is then
              advanced 3 steps by the C step and the six class keys must agree
              at every step.  Jumps are added to the STORED (post-align) point,
              never to the raw representative.
  controls    six mutations of the same source.  The four label mutations must
              FAIL the invariant; the slot-in-index mutation (the v15p defect)
              and the raw-x index must FAIL the partner row.  A control that
              cannot fail certifies nothing, so all six are shipped.
"""
import ctypes
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import drv_v15q as D

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE
LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
ROWS = 2048
NK = 500
NWALKS = 24
NPTS = 120
U4 = ctypes.c_uint64 * 4
G = (GX, GY)
Qs = None


def u4(v):
    a = U4()
    for i in range(4):
        a[i] = (v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
    return a


def val(a):
    return sum(a[i] << (64 * i) for i in range(4))


def padd(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    lam = (3 * x1 * x1) * pow(2 * y1, p - 2, p) % p if P == Q else (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def pmul(k, P=G):
    R, Q = None, P
    while k:
        if k & 1:
            R = padd(R, Q)
        Q = padd(Q, Q)
        k >>= 1
    return R


def phi(P, j=1):
    return None if P is None else (pow(BETA, j, p) * P[0] % p, P[1])


# orbit-generated row family: row k has jump [u_k]G + [v_k]Q and labels (u_k,v_k)
K0 = 0x5EC0FFEE1234567890ABCDEF0123456789ABCDEF0123456789ABCDEF01234567 % r
Qs = pmul(K0)
UV = [((0x1234567 * k + 0xDEADBEEF) % r, (0x89ABCDE * k + 0xFEEDFACE) % r) for k in range(ROWS)]
JUMPS = [pmul((UV[k][0] + UV[k][1] * K0) % r) for k in range(ROWS)]


def mk(name, defines=()):
    lib, Rec = D.build(name, p, r, BETA, LAM, 1024, defines)
    D.set_tbl(lib, JUMPS)
    return lib, Rec


def align(lib, Rec, P):
    rec = Rec()
    for i in range(4):
        rec.x[i] = (P[0] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
        rec.y[i] = (P[1] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
    sg = ctypes.c_int(0)
    slot = lib.v15o_align(ctypes.byref(rec), ctypes.byref(sg))
    return (val(rec.x), val(rec.y)), slot, sg.value, rec.idx


def step(lib, P, idx):
    x, y = u4(P[0]), u4(P[1])
    rc = lib.v15q_step(x, y, ctypes.c_uint64(idx))
    return rc, (val(x), val(y))


def label(lib, a, c, slot, sgn, k):
    A, C, JU, JV = u4(a), u4(c), u4(UV[k][0]), u4(UV[k][1])
    lib.v15o_label(A, C, ctypes.c_uint64(slot), ctypes.c_int(sgn), JU, JV)
    return val(A), val(C)


def mulC(lib, k, P):
    if k % r == 0:
        return None
    x, y, K = u4(0), u4(0), u4(k)
    rc = lib.v15q_mul(x, y, K, u4(P[0]), u4(P[1]))
    assert rc == 0
    return (val(x), val(y))


def invariant(lib, Rec, nwalks=NWALKS, nk=NK):
    bad = keybad = ybad = endbad = 0
    checks = 0
    rnd = random.Random(424242)
    for w in range(nwalks):
        k0 = rnd.randrange(1, r)
        P0 = pmul(k0)
        key, slot0, sgn0, idx = align(lib, Rec, P0)
        # labels describe the canonical START point, not a= c= 0
        mult = 1 if slot0 == 0 else (LAM if slot0 == 1 else LAM * LAM % r)
        a = (mult * k0 % r) * (-1 if sgn0 else 1) % r
        c = 0
        P = key
        for s in range(nk):
            kk = idx                     # the row the step will ADD; the label
            rc, Praw = step(lib, P, idx)  # must carry that row's (u,v), not the new one
            if rc != 0:
                bad += 1
                break
            key, slot, sgn, idx = align(lib, Rec, Praw)
            cand = min(pow(BETA, j, p) * Praw[0] % p for j in range(3))
            if key[0] != cand:
                keybad += 1
            if key[1] % 2:
                ybad += 1
            a, c = label(lib, a, c, slot, sgn, kk % ROWS)
            P = key
            if (s + 1) % 25 == 0 or s + 1 == nk:
                checks += 1
                if mulC(lib, a, G) != padd(mulC(lib, c, Qs), P):
                    endbad += 1
    return {"row": "shipped", "walks": nwalks, "steps_per_walk": nk,
            "step_failures": bad, "key_violations": keybad, "even_y_violations": ybad,
            "endpoint_checks": checks, "endpoint_violations": endbad,
            "verdict": "PASS" if bad == keybad == ybad == endbad == 0 else "FAIL"}


def partner(lib, Rec, npts=NPTS):
    mism_state = mism_idx = mism_key = 0
    tot = 0
    rnd = random.Random(99)
    for _ in range(npts):
        P = pmul(rnd.randrange(1, r))
        reps = [P, phi(P, 1), phi(P, 2), (P[0], (-P[1]) % p),
                (phi(P, 1)[0], (-P[1]) % p), (phi(P, 2)[0], (-P[1]) % p)]
        stored = []
        for Q in reps:
            key, slot, sgn, idx = align(lib, Rec, Q)
            stored.append((key, idx))
        tot += 1
        if len({s[0] for s in stored}) != 1:
            mism_state += 1
        if len({s[1] for s in stored}) != 1:
            mism_idx += 1
        # advance every STORED state 3 steps; the six class keys must agree
        cur = [s[0] for s in stored]
        idc = [s[1] for s in stored]
        for _s in range(3):
            nxt = []
            for j in range(6):
                rc, Praw = step(lib, cur[j], idc[j])
                if rc != 0:
                    mism_key += 1
                    nxt.append(cur[j])
                    continue
                key, slot, sgn, idx = align(lib, Rec, Praw)
                nxt.append(key)
                idc[j] = idx
            if len(set(nxt)) != 1:
                mism_key += 1
            cur = nxt
    return {"row": "shipped", "classes": tot, "representatives": 6,
            "stored_point_mismatches": mism_state, "stored_index_mismatches": mism_idx,
            "next_key_mismatches": mism_key,
            "verdict": "PASS" if mism_state == mism_idx == mism_key == 0 else "FAIL"}


out = {"kind": "harness-v15q", "curve": "secp256k1", "rows": [], "partner_rows": [],
       "note": ("stored = post-align point; jumps are added to the stored point; "
                "the six representatives are P, phi P, phi^2 P and their negatives")}

lib, Rec = mk("ship")
out["rows"].append(invariant(lib, Rec))
out["partner_rows"].append(partner(lib, Rec))

MUTS = [("control_mont_one_is_integer", ["MUT_MONT_ONE"], "inv"),
        ("control_negate_on_borrow", ["MUT_NEG_BORROW"], "inv"),
        ("control_jump_not_added", ["MUT_NO_JUMP"], "inv"),
        ("control_sign_not_consumed", ["MUT_NO_SIGN"], "inv"),
        ("control_slot_in_index_v15p", ["MUT_IDX_SLOT"], "part"),
        ("control_index_from_raw_x", ["MUT_IDX_RAW"], "part")]
for name, defs, kind in MUTS:
    lb, Rc = mk(name, defs)
    if kind == "inv":
        row = invariant(lb, Rc, nwalks=4, nk=120)
        row["row"] = name
        row["kind"] = "invariant"
        out["rows"].append(row)
    else:
        row = partner(lb, Rc, npts=20)
        row["row"] = name
        out["partner_rows"].append(row)

ok = (out["rows"][0]["verdict"] == "PASS"
      and out["partner_rows"][0]["verdict"] == "PASS"
      and all(r["verdict"] == "FAIL" for r in out["rows"][1:])
      and all(r["verdict"] == "FAIL" for r in out["partner_rows"][1:]))
out["verdict"] = "PASS" if ok else "FAIL"
out["wrong_direction_for_oracle_note"] = (
    "the oracle uses Python big integers for the canonical key and the C scalar "
    "multiplication for [a]G+[c]Q; the walk itself is C throughout")
os.makedirs(os.path.join(PKG, "evidence"), exist_ok=True)
json.dump(out, open(os.path.join(PKG, "evidence", "harness-v15q.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
sys.exit(0 if ok else 1)
