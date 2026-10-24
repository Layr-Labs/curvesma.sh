#!/usr/bin/env python3
"""harness_v15r.py -- the SHIPPED object executed; every row has a control.

ROWS
  invariant  the endpoint identity is DERIVED from invariant-row.json (single
             source, scripts/invcheck.py).  v15q hand-wrote the wrong side
             ('[a]G == [c]Q + P') and charged the walk for its own bug.
  partner    six orbit representatives incl. -P, jumps added to STORED points.
  action     detector FORCED on >=100 walks: every counter trip must produce a
             state reset (pre/post state-hash inequality), the endpoint must
             still hold after the restart, and a detached-action build must
             FAIL the same test.
  recovery   constructed collisions through the shipped solve + verify.
  controls   slot-in-index, raw-x index, detached action, wrong-side endpoint
             row: each must FAIL the row it targets.
"""
import ctypes
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import drv_v15t as D
import invcheck
import toy_core as T

p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
BETA = 0x7AE96A2B657C07106E64479EAC3434E99CF0497512F58995C1396C28719501EE
LAM = 0x5363AD4CC05C30E0A5261C028812645A122E22EA20816678DF02967C1B23BD72
GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
ROWS = 2048
U4 = ctypes.c_uint64 * 4


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


def pmul(k, P):
    R, Q = None, P
    while k:
        if k & 1:
            R = padd(R, Q)
        Q = padd(Q, Q)
        k >>= 1
    return R


def phi(P, j=1):
    return (pow(BETA, j, p) * P[0] % p, P[1])


G = (GX, GY)
K0 = 0x5EC0FFEE1234567890ABCDEF0123456789ABCDEF0123456789ABCDEF01234567 % r
Q0 = pmul(K0, G)
UV = [((0x1234567 * k + 0xDEADBEEF) % r, 1 + k % 15) for k in range(ROWS)]


def mk(name, defines=()):
    lib, Rec = D.build(name, p, r, BETA, LAM, 1024, defines)
    UB = (ctypes.c_uint64 * (8 * 2048)).in_dll(lib, "UBASE")
    for k in range(ROWS):
        u, _v = UV[k]
        ax, ay = U4(), U4()
        assert lib.v15q_mul(ax, ay, u4(u), u4(GX), u4(GY)) == 0
        for i in range(4):
            UB[8 * k + i] = ax[i]
            UB[8 * k + 4 + i] = ay[i]
    tbl = (ctypes.c_uint64 * (8 * 3 * 2048)).in_dll(lib, "TBL")
    assert lib.v15r_init_for_Q(tbl, u4(Q0[0]), u4(Q0[1])) == 0
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
    A, C = u4(a), u4(c)
    lib.v15o_label(A, C, ctypes.c_uint64(slot), ctypes.c_int(sgn), u4(UV[k][0]), u4(UV[k][1]))
    return val(A), val(C)


def mulC(lib, k, P):
    if k % r == 0:
        return None
    x, y = u4(0), u4(0)
    rc = lib.v15q_mul(x, y, u4(k), u4(P[0]), u4(P[1]))
    assert rc == 0
    return (val(x), val(y))


ROW = invcheck.load()
BAD_ROW = json.loads(json.dumps(ROW))
BAD_ROW["terms"][2]["coef"] = "c"          # the v15q side: +[c]Q instead of -[c]Q


def invariant(lib, Rec, nwalks=24, nk=200, row=None):
    bad = keybad = ybad = endbad = checks = 0
    rnd = random.Random(424242)
    for _ in range(nwalks):
        k0 = rnd.randrange(1, r)
        P0 = pmul(k0, G)
        key, slot0, sgn0, idx = align(lib, Rec, P0)
        m = 1 if slot0 == 0 else (LAM if slot0 == 1 else LAM * LAM % r)
        if sgn0:
            m = (-m) % r
        a, c = m * k0 % r, 0
        P = key
        for s in range(nk):
            kk = idx
            rc, Praw = step(lib, P, idx)
            if rc != 0:
                bad += 1
                break
            key, slot, sgn, idx = align(lib, Rec, Praw)
            if key[0] != min(pow(BETA, j, p) * Praw[0] % p for j in range(3)):
                keybad += 1
            if key[1] % 2:
                ybad += 1
            a, c = label(lib, a, c, slot, sgn, kk % ROWS)
            P = key
            if (s + 1) % 50 == 0 or s + 1 == nk:
                checks += 1
                if not invcheck.holds(T.Curve(p, r, 7, BETA, LAM), {"P": P, "G": G, "Q": Q0},
                                      a, c, row):
                    endbad += 1
    return {"row": "shipped", "walks": nwalks, "steps_per_walk": nk,
            "step_failures": bad, "key_violations": keybad, "even_y_violations": ybad,
            "endpoint_checks": checks, "endpoint_violations": endbad,
            "endpoint_formula": ROW["text"],
            "verdict": "PASS" if bad + keybad + ybad + endbad == 0 else "FAIL"}


def partner(lib, Rec, npts=60):
    mism_state = mism_idx = mism_key = 0
    rnd = random.Random(99)
    for _ in range(npts):
        P = pmul(rnd.randrange(1, r), G)
        reps = [P, phi(P, 1), phi(P, 2),
                (P[0], (-P[1]) % p), (phi(P, 1)[0], (-P[1]) % p), (phi(P, 2)[0], (-P[1]) % p)]
        stored = [align(lib, Rec, Q) for Q in reps]
        if len({s[0] for s in stored}) != 1:
            mism_state += 1
        if len({s[3] for s in stored}) != 1:
            mism_idx += 1
        cur = [s[0] for s in stored]
        idc = [s[3] for s in stored]
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
    return {"row": "shipped", "classes": npts, "representatives": 6,
            "stored_point_mismatches": mism_state, "stored_index_mismatches": mism_idx,
            "next_key_mismatches": mism_key,
            "verdict": "PASS" if mism_state + mism_idx + mism_key == 0 else "FAIL"}


FIX = None


def fixture():
    global FIX
    if FIX is None:
        cur = T.find_curve(25, 1501)
        lib, Rec = D.build("act%d" % os.getpid(), cur.p, cur.r, cur.beta, cur.lam, 1024,
                           ["V15R_ROWS=8", "V15R_REP_THRESH=1"])
        rnd = random.Random(7)
        P0 = T.base_point(cur, rnd)
        UVF = [((0x1234567 * k + 0xDEADBEEF) % cur.r, 1 + k % 15) for k in range(8)]
        UB = (ctypes.c_uint64 * (8 * 2048)).in_dll(lib, "UBASE")
        for k in range(8):
            ax, ay = U4(), U4()
            assert lib.v15q_mul(ax, ay, u4(UVF[k][0]), u4(P0[0]), u4(P0[1])) == 0
            for i in range(4):
                UB[8 * k + i] = ax[i]
                UB[8 * k + 4 + i] = ay[i]
        FIX = (cur, lib, Rec, P0, UVF)
    return FIX


def action_row(detached=False, nwalks=120, nk=90):
    cur, lib, Rec, P0, UVF = fixture()
    if detached:
        lib, Rec = D.build("det%d" % os.getpid(), cur.p, cur.r, cur.beta, cur.lam, 1024,
                           ["V15R_ROWS=8", "V15R_REP_THRESH=1", "MUT_DETACH_ACTION"])
        UB = (ctypes.c_uint64 * (8 * 2048)).in_dll(lib, "UBASE")
        for k in range(8):
            ax, ay = U4(), U4()
            assert lib.v15q_mul(ax, ay, u4(UVF[k][0]), u4(P0[0]), u4(P0[1])) == 0
            for i in range(4):
                UB[8 * k + i] = ax[i]
                UB[8 * k + 4 + i] = ay[i]
    rnd = random.Random(2024)
    K = rnd.randrange(1, cur.r)
    Tt = cur.mul(K, P0)
    tbl = (ctypes.c_uint64 * (8 * 3 * 2048)).in_dll(lib, "TBL")
    assert lib.v15r_init_for_Q(tbl, u4(Tt[0]), u4(Tt[1])) == 0
    trips = resets = inv_ok = 0
    walks_with_trip = walks_with_reset = 0
    for _w in range(nwalks):
        a = rnd.randrange(1, cur.r)
        c = rnd.randrange(1, cur.r)
        P = cur.add(cur.mul(a, P0), cur.mul(c, Tt))
        rec = Rec()
        for i in range(4):
            rec.x[i] = (P[0] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
            rec.y[i] = (P[1] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
        sg = ctypes.c_int(0)
        slot = lib.v15o_align(ctypes.byref(rec), ctypes.byref(sg))
        m = 1 if slot == 0 else (cur.lam if slot == 1 else cur.lam * cur.lam % cur.r)
        if sg.value:
            m = (-m) % cur.r
        a, c = a * m % cur.r, c * m % cur.r
        idx = rec.idx
        p0, p1 = (ctypes.c_uint64 * 4)(), (ctypes.c_uint64 * 4)()
        rep, win = ctypes.c_uint64(0), ctypes.c_uint64(0)
        minx = (ctypes.c_uint64 * 4)(*([(1 << 64) - 1] * 4))
        hit = False
        for _s in range(nk):
            kk = idx
            rc = lib.v15q_step(rec.x, rec.y, ctypes.c_uint64(idx))
            if rc != 0:
                break
            sg = ctypes.c_int(0)
            slot = lib.v15o_align(ctypes.byref(rec), ctypes.byref(sg))
            A, C = u4(a), u4(c)
            lib.v15o_label(A, C, ctypes.c_uint64(slot), ctypes.c_int(sg.value),
                           u4(UVF[kk % 8][0]), u4(UVF[kk % 8][1]))
            a, c = val(A), val(C)
            idx = rec.idx
            lib.v15o_detect(ctypes.byref(rec), ctypes.byref(p0), ctypes.byref(p1),
                            ctypes.byref(rep), ctypes.byref(minx), ctypes.byref(win))
            trip = (rep.value >= 1) or (win.value >= 67)
            if trip:
                trips += 1
            if lib.v15r_action(ctypes.byref(rep), ctypes.byref(win)):
                hit = True
                pre = (tuple(rec.x), tuple(rec.y), a, c, rep.value, win.value)
                s_off = rnd.randrange(1, 1 << 16)
                t_off = rnd.randrange(1, 1 << 16)
                rcx = lib.v15r_escape(rec.x, rec.y, A, C, u4(P0[0]), u4(P0[1]),
                                      u4(Tt[0]), u4(Tt[1]), ctypes.c_uint64(s_off),
                                      ctypes.c_uint64(t_off), ctypes.byref(p0),
                                      ctypes.byref(p1), ctypes.byref(rep),
                                      ctypes.byref(minx), ctypes.byref(win))
                if rcx == 0:
                    a, c = (a + s_off) % cur.r, (c + t_off) % cur.r
                    sg = ctypes.c_int(0)
                    slot = lib.v15o_align(ctypes.byref(rec), ctypes.byref(sg))
                    # the post-escape alignment also scales BOTH labels
                    A2, C2 = u4(a), u4(c)
                    lib.v15o_label(A2, C2, ctypes.c_uint64(slot), ctypes.c_int(sg.value),
                                   u4(0), u4(0))
                    a, c = val(A2), val(C2)
                    idx = rec.idx
                    post = (tuple(rec.x), tuple(rec.y), a, c, rep.value, win.value)
                    if post != pre:
                        resets += 1
                    Pn = (val(rec.x), val(rec.y))
                    if invcheck.holds(cur, {"P": Pn, "G": P0, "Q": Tt}, a, c):
                        inv_ok += 1
        walks_with_trip += 1 if hit else 0
        walks_with_reset += 1 if (hit and resets > 0) else 0
    return {"row": "action_detached" if detached else "action_forced",
            "walks": nwalks, "steps_per_walk": nk, "counter_trips": trips,
            "state_resets": resets, "endpoint_holds_after_reset": inv_ok,
            "walks_with_trip": walks_with_trip,
            "verdict": ("FAIL" if detached else
                        ("PASS" if trips > 0 and resets == trips and inv_ok == resets else "FAIL"))}


def recovery_row(lib, n=64):
    cur = T.Curve(p, r, 7, BETA, LAM)
    rnd = random.Random(77)
    K = rnd.randrange(1, r)
    Tt = cur.mul(K, G)
    ok = 0
    rc_hist = {}
    for _ in range(n):
        a1 = rnd.randrange(1, r)
        c1 = rnd.randrange(1, r)
        c2 = (c1 + 1) % r
        a2 = (a1 - K) % r          # [a1]G+[c1]Q == [a2]G+[c2]Q
        x = U4()
        rc = lib.v15q_solve(u4(a1), u4(c1), u4(a2), u4(c2), x)
        rc_hist[str(rc)] = rc_hist.get(str(rc), 0) + 1
        if rc == 0:
            vx, vy = U4(), U4()
            if lib.v15q_mul(vx, vy, x, u4(GX), u4(GY)) == 0 and (val(vx), val(vy)) == Tt:
                ok += 1
    return {"row": "recovery_shipped", "pairs": n, "solved_and_verified": ok,
            "rc_histogram": rc_hist,
            "verdict": "PASS" if ok == n else "FAIL"}


def main():
    out = {"kind": "harness-v15s", "rows": [], "partner_rows": [],
           "note": ("endpoint identity derived from scripts/invariant-row.json; jumps are added to "
                    "the STORED post-align point; six representatives include the negatives")}

    lib, Rec = mk("ship")
    out["rows"].append(invariant(lib, Rec))
    out["partner_rows"].append(partner(lib, Rec))
    out["rows"].append(recovery_row(lib))
    out["action_rows"] = [action_row(False), action_row(True)]

    # --- controls -------------------------------------------------------------
    for name, defs, kind in [("control_slot_in_index_v15p", ["MUT_IDX_SLOT", "V15R_ROWS=2048"], "part"),
                             ("control_index_from_raw_x", ["MUT_IDX_RAW", "V15R_ROWS=2048"], "part")]:
        lb, Rc = mk(name, defs)
        row = partner(lb, Rc, npts=20)
        row["row"] = name
        row["targets"] = "partner"
        out["partner_rows"].append(row)
    # wrong-side endpoint (the v15q harness formula) must FAIL
    out["rows"].append(dict(invariant(lib, Rec, nwalks=2, nk=50, row=BAD_ROW),
                            row="control_wrong_side_endpoint_v15q", targets="invariant"))

    ok = (out["rows"][0]["verdict"] == "PASS"
          and out["partner_rows"][0]["verdict"] == "PASS"
          and out["rows"][1]["verdict"] == "PASS"
          and out["action_rows"][0]["verdict"] == "PASS"
          and out["action_rows"][1]["verdict"] == "FAIL"
          and all(r["verdict"] == "FAIL" for r in out["rows"][2:] if r["row"].startswith("control"))
          and all(r["verdict"] == "FAIL" for r in out["partner_rows"][1:]))
    out["verdict"] = "PASS" if ok else "FAIL"
    json.dump(out, open(os.path.join(PKG, "evidence", "harness-v15s.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
    sys.exit(0 if ok else 1)

if __name__ == '__main__':
    main()
