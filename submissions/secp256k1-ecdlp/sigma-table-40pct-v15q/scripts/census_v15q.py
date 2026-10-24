#!/usr/bin/env python3
"""census_v15q.py -- the success census run THROUGH THE SHIPPED OBJECT.

Every mechanism that the claim depends on is the compiled object, not a Python
stand-in:

  * class key, canonical representative, sign, argmin slot and the row index --
    v15o_align (scripts/walk_v15o.c, compiled for the fixture curve);
  * the affine point step -- v15q_step (same TU);
  * the exponent-pair bookkeeping -- v15o_label;
  * the DP predicate, the two-step repeat counter and the running-minimum
    checkpoint -- v15o_detect;
  * RETENTION -- v15o_report, the production four-probe pool over the
    2^node_bits-node array, insert/lookup VERBATIM.  The v15p objection (a
    dictionary retains differently from the priced pool) cannot apply: there is
    no dictionary.  The fixture's pool is the same code at the same probe bound;
    only the node count and the field width are scaled.

STATISTIC.  A trial draws a target T=[K]G and runs w=2 SEQUENTIAL, INDEPENDENT
runs of that production pool (fresh coins, fresh pool, each run charged its own
budget).  A trial is a success iff one of the two runs returns a two-record pair
whose solved exponent verifies, [x]G == T.  The reported interval is the exact
one-sided Clopper-Pearson lower bound on THAT directly measured fraction -- no
1-(1-p)^w expression is used anywhere, so no independence-of-p assumption is
needed.  mu is measured from the same runs (class visits / distinct classes).
"""
import ctypes
import json
import math
import os
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import drv_v15q as D
import toy_core as T

getcontext().prec = 50
SEED = 1501
BITS = 25
R_TABLE = 2048
NTRAILS = 16
NODE_BITS = 16
TRIALS = int(os.environ.get("V15Q_TRIALS", "512"))
MU = F(1107, 1000)
FLOOR = F(39, 100)
W_TRIALS = 2


def cp95_lower(k, n, alpha=Decimal("0.05")):
    """one-sided Clopper-Pearson lower bound, exact binomial tail, bisection."""
    k, n = int(k), int(n)
    if k == 0:
        return Decimal(0)
    lnC = sum(Decimal(i).ln() for i in range(n - k + 1, n + 1)) - \
          sum(Decimal(i).ln() for i in range(1, k + 1))

    def tail(p):
        p = Decimal(p)
        q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal("1e-45") * s:
                break
        return s

    lo, hi = Decimal("0.0"), Decimal("0.999")
    for _ in range(90):
        mid = (lo + hi) / 2
        if tail(mid) > alpha:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-14"):
            break
    return lo


def s_ops(N6, mu=MU, floor=FLOOR):
    q = Decimal(1000 - int(floor * 1000)) / Decimal(1000)
    L = F(int((-q.ln()) * Decimal(10) ** 70) + 1, 10 ** 70)
    X2 = 2 * mu * F(N6 * R_TABLE, R_TABLE - 1) * L
    s = math.isqrt(X2.numerator // X2.denominator)
    if F(s * s) < X2:
        s += 1
    return s


CURVE = None
LIB = None
REC = None
U4 = ctypes.c_uint64 * 4


def u4(v):
    a = U4()
    for i in range(4):
        a[i] = (v >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
    return a


def val(a):
    return sum(a[i] << (64 * i) for i in range(4))


def init_worker():
    """One fixture, ONE shared jump table (production uses one table too)."""
    global CURVE, LIB, REC, P0, UVS
    cur = T.find_curve(BITS, SEED)
    CURVE = cur
    LIB, REC = D.build("cen%d" % os.getpid(), cur.p, cur.r, cur.beta, cur.lam, 1 << NODE_BITS)
    rnd = random.Random(SEED + 1)
    P0 = T.base_point(cur, rnd)
    UVS, pts, _K0 = build_table(cur, P0, rnd)
    Tarr = (ctypes.c_uint64 * (8 * 3 * 2048)).in_dll(LIB, "TBL")
    for k in range(3 * 2048):
        Q = pts[k % R_TABLE]
        for i in range(4):
            Tarr[8 * k + i] = (Q[0] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
            Tarr[8 * k + 4 + i] = (Q[1] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF


def build_table(cur, P0, rnd):
    """row k: [u_k]G + [v_k]Q with Q=[K0]G, deterministic in k."""
    K0 = 0x5EC0FFEE1234567890ABCDEF0123456789ABCDEF0123456789ABCDEF01234567 % cur.r
    Q = cur.mul(K0, P0)
    uvs, pts = [], []
    for k in range(R_TABLE):
        u = (0x1234567 * k + 0xDEADBEEF) % cur.r
        v = (0x89ABCDE * k + 0xFEEDFACE) % cur.r
        uvs.append((u, v))
        pts.append(cur.add(cur.mul(u, P0), cur.mul(v, Q)))
    return uvs, pts, K0


def one_run(cur, lib, Rec, P0, Tt, uvs, S, rnd, measure_mu):
    """One production run: NTRAILS walkers, shared four-probe pool, budget S."""
    pool = (ctypes.c_uint64 * (12 * (1 << NODE_BITS))).in_dll(lib, "NODE_HEAD")
    ctypes.memset(pool, 0, ctypes.sizeof(pool))
    mask = (1 << NODE_BITS) - 1
    st = []
    for _ in range(NTRAILS):
        a = rnd.randrange(1, cur.r)
        c = rnd.randrange(1, cur.r)
        P = cur.add(cur.mul(a, P0), cur.mul(c, Tt))
        rec = Rec()
        for i in range(4):
            rec.x[i] = (P[0] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
            rec.y[i] = (P[1] >> (64 * i)) & 0xFFFFFFFFFFFFFFFF
        sg = ctypes.c_int(0)
        slot = lib.v15o_align(ctypes.byref(rec), ctypes.byref(sg))
        mult = 1 if slot == 0 else (cur.lam if slot == 1 else cur.lam * cur.lam % cur.r)
        a = a * mult % cur.r
        if sg.value:
            a = (-a) % cur.r
        st.append(dict(rec=rec, a=a, c=c, idx=rec.idx,
                       prev0=(ctypes.c_uint64 * 4)(), prev1=(ctypes.c_uint64 * 4)(),
                       rep=ctypes.c_uint64(0), minx=(ctypes.c_uint64 * 4)(),
                       win=ctypes.c_uint64(0), seen=set()))
        for k in range(4):
            st[-1]["minx"][k] = (1 << 64) - 1
    visits = 0
    seen_all = set()
    matches = 0
    for step in range(S):
        for tr in st:
            rec = tr["rec"]
            rc = lib.v15q_step(rec.x, rec.y, ctypes.c_uint64(tr["idx"]))
            if rc != 0:
                continue
            sg = ctypes.c_int(0)
            slot = lib.v15o_align(ctypes.byref(rec), ctypes.byref(sg))
            u, v = uvs[tr["idx"] % R_TABLE]
            A, C = u4(tr["a"]), u4(tr["c"])
            lib.v15o_label(A, C, ctypes.c_uint64(slot), ctypes.c_int(sg.value),
                           u4(u), u4(v))
            tr["a"], tr["c"] = val(A), val(C)
            tr["idx"] = rec.idx
            visits += 1
            if measure_mu:
                seen_all.add(val(rec.x))
            dp = lib.v15o_detect(ctypes.byref(rec), ctypes.byref(tr["prev0"]),
                                 ctypes.byref(tr["prev1"]), ctypes.byref(tr["rep"]),
                                 ctypes.byref(tr["minx"]), ctypes.byref(tr["win"]))
            if dp:
                oa, oc = U4(), U4()
                hit = lib.v15o_report(pool, ctypes.c_uint64(mask), ctypes.byref(rec),
                                      A, C, oa, oc)
                if hit:
                    matches += 1
                    x = U4()
                    rc2 = lib.v15q_solve(A, C, oa, oc, x)
                    if rc2 == 0:
                        vx, vy = U4(), U4()
                        rc3 = lib.v15q_mul(vx, vy, x, u4(P0[0]), u4(P0[1]))
                        if rc3 == 0 and (val(vx), val(vy)) == tuple(Tt):
                            return dict(ok=True, steps=step + 1, matches=matches,
                                        visits=visits, distinct=len(seen_all))
        # early exit: nothing left to do once every class is stored
    return dict(ok=False, steps=S, matches=matches, visits=visits,
                distinct=len(seen_all))


def trial(args):
    idx, = args
    cur, uvs = CURVE, UVS
    rnd = random.Random(90000 + idx)
    N6 = (cur.r - 1) // 6
    S = s_ops(N6)
    K = rnd.randrange(1, cur.r)
    Tt = cur.mul(K, P0)
    res = []
    for w in range(W_TRIALS):
        r = one_run(cur, LIB, REC, P0, Tt, uvs, S, random.Random(70000 + 17 * idx + w), w == 0)
        r["K_ok"] = None
        res.append(r)
        if r["ok"]:
            break
    return dict(idx=idx, ok=any(r["ok"] for r in res), runs=res,
                visits=res[0]["visits"], distinct=res[0]["distinct"])


def main():
    cur = T.find_curve(BITS, SEED)
    N6 = (cur.r - 1) // 6
    S = s_ops(N6)
    out = dict(kind="census-v15q", curve=dict(name="q%d" % BITS, p_bits=BITS,
                                              p=str(cur.p), order=str(cur.r),
                                              beta=str(cur.beta), lam=str(cur.lam),
                                              R=R_TABLE, ntrails=NTRAILS),
               node_bits=NODE_BITS, nodes=1 << NODE_BITS, trials=TRIALS, w=W_TRIALS,
               S_ops=str(S), S_ops_log2=math.log2(S),
               method=("w=2 sequential independent runs of the SHIPPED four-probe pool "
                       "(v15o_report) with the SHIPPED align/label/detect/step; success "
                       "measured directly on the 2-run set"))
    with Pool(8, initializer=init_worker) as pool:
        rows = pool.map(trial, [(i,) for i in range(TRIALS)])
    k2 = sum(1 for r in rows if r["ok"])
    k1 = sum(1 for r in rows if r["runs"][0]["ok"])
    vis = sum(r["visits"] for r in rows)
    dis = sum(r["distinct"] for r in rows)
    out.update(trials=TRIALS, success_2run=k2, success_1run=k1,
               p_1run=float(F(k1, TRIALS)), p_2run=float(F(k2, TRIALS)),
               cp95_2run=str(cp95_lower(k2, TRIALS)),
               cp95_1run=str(cp95_lower(k1, TRIALS)),
               mu_measured=float(F(vis, dis)), mu_measured_exact=str(F(vis, dis)),
               visits=vis, distinct_classes=dis,
               matches=sum(sum(r["matches"] for r in x["runs"]) for x in rows))
    out["floor_met"] = bool(Decimal(out["cp95_2run"]) >= Decimal("0.39"))
    json.dump(out, open(os.path.join(PKG, "evidence", "census-v15q.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "curve"}, indent=1))


if __name__ == "__main__":
    main()
