#!/usr/bin/env python3
"""fullstore_census.py -- p256 v13: INSERT-TIME FULL-STORE detection census.

WHY THIS FILE EXISTS.  v12's detector stored a class only when the report gate
(DP density 2^-d) fired, so measured success depended on d: the ladder read
-1.326 bits of success odds per bit of density, and extrapolating 55 bits to the
production density 2^-60 was the objection that killed the ticket.  v13 removes
the gate: EVERY visited class is inserted and looked up, and a repeat is seen at
INSERT time.  Success then depends only on x = visits/sqrt(N), a scale-free
quantity, and not on any density parameter at all.

THE WALK IS UNCHANGED.  Same qwalk.step, same orbit-row map, same lag-1/lag-3
escape, same fixed start block.  Only the DETECTOR differs.  Because the whole
trajectory is retained, the same traces are also read back through the v12
gate (d = 5, 8, 11) as a MATCHED ABLATION: identical trajectories, gate on vs
gate off.  That row shows the v12 slope is a property of the gate, not of the
walk, and that deleting the gate deletes the dimension.

RECOVERY.  Every step carries the invariant a + c*x == class (mod q).  Two steps
that land on the same class therefore give x = (a1-a2)/(c2-c1) mod q; the
census recomputes x and compares it with the trial's own challenge.  A detected
repeat is only counted as a success when that comparison passes.
"""
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402

p = int('ffffffff00000001000000000000000000000000ffffffffffffffffffffffff', 16)
r = int('ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551', 16)
N_CLASSES = (r - 1) // 2
R_TAB = 2 ** 20
A_H1 = F(1015, 1000)
MODEL_SUCCESS = F(3903, 10000)
getcontext().prec = 80
_z0 = Decimal(1) - Decimal(MODEL_SUCCESS.numerator) / Decimal(MODEL_SUCCESS.denominator)
_ln0 = -_z0.ln()


def frac_of_dec(dd, digits=70):
    return F(int(dd * Decimal(10) ** digits) + (1 if dd >= 0 else -1), 10 ** digits)


N_CELLS = 2 * R_TAB
Nprime = F(N_CLASSES * N_CELLS, N_CELLS - 1)
_x2 = 2 * Nprime * frac_of_dec(_ln0)
S0 = math.isqrt(_x2.numerator // _x2.denominator)
if F(S0 * S0) < _x2:
    S0 += 1
_S = A_H1 * A_H1 * F(S0) * F(S0)
S_OPS = math.isqrt(_S.numerator // _S.denominator)
if F(S_OPS * S_OPS) < _S:
    S_OPS += 1
# The census must run at the budget the ledger CHARGES: x = 1.6002 visits/sqrt(N)
# (v12's executed shape, census-design.json).  S_OPS/sqrt(N) above is the
# birthday shape at the 0.39 floor; the ledger charges A_H1 * x_meas, so the
# fixture is run at x_meas and the ledger over-charges it by A_H1.
SHAPE = float(os.environ.get("SHAPE", repr(1.6002022719673725)))

ALPHA = Decimal("0.05")


def cp95_lower(k, n):
    """One-sided Clopper-Pearson lower bound, bracket (0,1) (clip bug fixed)."""
    k, n = int(k), int(n)
    if n == 0 or k == 0:
        return Decimal(0)
    if k >= n:
        return ALPHA ** (Decimal(1) / Decimal(n))
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p_):
        p_ = Decimal(p_)
        q_ = Decimal(1) - p_
        t = (lnC + k * p_.ln() + (n - k) * q_.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p_ / ((j + 1) * q_)
            s += t
            if t < Decimal("1e-40") * s:
                break
        return s

    lo, hi = Decimal("0.000001"), Decimal("0.999999")
    for _ in range(300):
        mid = (lo + hi) / 2
        if tail(mid) > ALPHA:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-14"):
            break
    return lo


def make_start(U, bt_max, w, seed):
    q = U.q
    rng = np.random.default_rng(seed & 0x7FFFFFFF)
    xt = rng.integers(1, q - 1, size=bt_max).astype(np.uint64)
    xv = np.repeat(xt, w).astype(np.uint64)
    lanes = bt_max * w
    a = rng.integers(0, q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, q, size=lanes).astype(np.uint64)
    t0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                   for ai, ci, xi in zip(a.tolist(), cc.tolist(), xv.tolist())],
                  dtype=np.uint64)
    canon, e, _h = U.canon(t0)
    return xt, xv, (e * a) % q, (e * cc) % q, canon


def walk_full(U, canon0, a0, cc0, xv, per, esc, lag, w, bt, coins_seed=0xC0FFEE):
    """Shipped walk; returns every step's (class, a, c) -- the full store."""
    lanes = canon0.size
    canon = canon0.copy(); a = a0.copy(); cc = cc0.copy()
    hist = [canon.copy(), canon.copy(), canon.copy()]
    esc_run = np.zeros(lanes, dtype=np.int64)
    coins = np.random.default_rng(coins_seed)
    CLS = np.empty((per, lanes), dtype=np.uint32)
    A = np.empty((per, lanes), dtype=np.uint64)
    C = np.empty((per, lanes), dtype=np.uint64)
    fires = fires1 = fires3 = restarts = 0
    inv_bad = 0
    kcap = qwalk.K_ESC if esc else None
    for t in range(per):
        cs, as_, cs_ = canon, a, cc
        cn, a2, cc2, _h, _pv, fired1 = qwalk.step(
            U, cs, as_, cs_, xv, 0, False, hist[0], esc, 0, t,
            None, None, None, coins)
        fired1 = fired1 if esc else np.zeros(lanes, dtype=bool)
        fired3 = np.zeros(lanes, dtype=bool)
        if lag >= 3:
            extra = (cn == hist[1]) & ~fired1
            if extra.any():
                idx = np.flatnonzero(extra)
                m = (U.col(cs[idx]) + 1) % U.R
                cn2, a3, cc3, _e = qwalk._apply(U, cs[idx], as_[idx], cs_[idx], xv[idx], m)
                cn = cn.copy(); a2 = a2.copy(); cc2 = cc2.copy()
                cn[idx] = cn2; a2[idx] = a3; cc2[idx] = cc3
                fired3 = extra
        pv = cs.copy()
        rst = np.zeros(lanes, dtype=bool)
        if kcap is not None:
            esc_run[:] = np.where(fired1 | fired3, esc_run + 1, 0)
            over = esc_run > kcap
            if over.any():
                idx = np.flatnonzero(over)
                cn3, a4, cc4, _a0f, _c0f = qwalk.fresh_state(U, xv[idx], coins, idx.size)
                cn = cn.copy(); a2 = a2.copy(); cc2 = cc2.copy()
                cn[idx] = cn3; a2[idx] = a4; cc2[idx] = cc4
                esc_run[idx] = 0
                rst[idx] = True
        a, cc = a2, cc2
        canon = cn
        CLS[t] = canon.astype(np.uint32)
        A[t] = a
        C[t] = cc
        fires1 += int(fired1.sum()); fires3 += int(fired3.sum())
        fires += int((fired1 | fired3).sum()); restarts += int(rst.sum())
        hist = [pv, hist[0], hist[1]]
        if t % 16 == 0:
            sl = slice(0, min(96, lanes))
            inv_bad += sum(1 for ai, ci, xi, cni in zip(a[sl].tolist(), cc[sl].tolist(),
                                                       xv[sl].tolist(), canon[sl].tolist())
                           if (ai + ci * xi) % U.q != cni)
    return {"CLS": CLS, "A": A, "C": C, "fires": fires, "fires_lag1": fires1,
            "fires_lag3": fires3, "restarts": restarts, "inv_bad": inv_bad}


def first_repeat(U, cls, a, c, xt_b, w, per, dp_d=0):
    """Earliest insert-time repeat (dp_d=0: every step) -> (found, solved, step).

    dp_d>0 is the v12 gate replayed on the SAME trace: only gated steps are
    inserted, so the repeat can be missed.  This is the matched ablation.
    """
    q = U.q
    if dp_d > 0:
        keep = np.flatnonzero(U.dp(cls, dp_d))
        if keep.size < 2:
            return False, False, -1
        ordr = keep[np.argsort(cls[keep], kind="stable")]
    else:
        ordr = np.argsort(cls, kind="stable")
    sv = cls[ordr]
    eq = np.flatnonzero(sv[1:] == sv[:-1])
    if eq.size == 0:
        return False, False, -1
    # first detection = smallest SECOND index among the equal pairs
    later = ordr[eq + 1]
    best = np.lexsort((ordr[eq], later))[0]
    i1, i2 = int(ordr[eq[best]]), int(ordr[eq[best] + 1])
    if i2 < i1:
        i1, i2 = i2, i1
    den = (int(c[i2]) - int(c[i1])) % q
    if den == 0:
        return True, False, i2 // w
    cand = ((int(a[i1]) - int(a[i2])) * pow(den, -1, q)) % q
    return True, cand == int(xt_b), i2 // w


def run_cell(qmin, R, W, bt, seed, scale, out, chunk=200):
    U = qwalk.Universe(qmin, R, seed=seed)
    assert U.q < (1 << 32), ("uint64 canon wrap", U.q)
    per = int(round(float(SHAPE) * scale * math.sqrt(U.classes) / W))
    shape_real = per * W / math.sqrt(U.classes)
    assert abs(shape_real - float(SHAPE) * scale) <= 5e-3, (shape_real, float(SHAPE) * scale)
    xt, xv, a_all, cc_all, canon_all = make_start(U, bt, W, seed)
    acc = {"succ": 0, "det": 0, "meet": [], "fires": 0, "restarts": 0, "inv_bad": 0,
           "dp": {}}
    for dd in (5, 8, 11):
        acc["dp"][dd] = {"succ": 0, "det": 0}
    for c0 in range(0, bt, chunk):
        b = min(chunk, bt - c0)
        sl = slice(c0 * W, (c0 + b) * W)
        res = walk_full(U, canon_all[sl], a_all[sl], cc_all[sl], xv[sl],
                        per, True, 3, W, b)
        acc["fires"] += res["fires"]; acc["restarts"] += res["restarts"]
        acc["inv_bad"] += res["inv_bad"]
        for j in range(b):
            cls = res["CLS"][:, j * W:(j + 1) * W].ravel()
            av = res["A"][:, j * W:(j + 1) * W].ravel()
            cv = res["C"][:, j * W:(j + 1) * W].ravel()
            found, solved, step = first_repeat(U, cls, av, cv, xt[c0 + j], W, per, 0)
            acc["det"] += int(found); acc["succ"] += int(solved)
            if found and step >= 0:
                acc["meet"].append(step)
            for dd in (5, 8, 11):
                f2, s2, _ = first_repeat(U, cls, av, cv, xt[c0 + j], W, per, dd)
                acc["dp"][dd]["det"] += int(f2); acc["dp"][dd]["succ"] += int(s2)
    n = bt
    rec = {"qmin": qmin, "q": int(U.q), "classes": int(U.classes),
           "classes_log2": math.log2(U.classes), "R": R, "walkers": W,
           "trials": n, "per_walker_steps": per, "visits_per_trial": per * W,
           "x_visits_over_sqrtN": shape_real, "shape_scale": scale,
           "fullstore": {"succ": acc["succ"], "trials": n,
                         "success_rate": acc["succ"] / n,
                         "cp95_success_f": float(cp95_lower(acc["succ"], n)),
                         "cp95_success": str(cp95_lower(acc["succ"], n)),
                         "detected": acc["det"],
                         "detection_given_repeat": (acc["det"] / n),
                         "model_birthday": 1 - math.exp(-shape_real ** 2 / 2),
                         "meet_step_median": (sorted(acc["meet"])[len(acc["meet"]) // 2]
                                              if acc["meet"] else None),
                         "meet_step_max": (max(acc["meet"]) if acc["meet"] else None)},
           "gate_replay": {}, "escape_fires": acc["fires"],
           "restarts": acc["restarts"], "invariant_violations_sampled": acc["inv_bad"]}
    for dd in (5, 8, 11):
        s_ = acc["dp"][dd]["succ"]
        rec["gate_replay"]["d%d" % dd] = {
            "succ": s_, "trials": n, "success_rate": s_ / n,
            "cp95_success_f": float(cp95_lower(s_, n)),
            "missing_vs_fullstore": acc["succ"] - s_}
    out.append(rec)
    print("[cell] qmin=%d classes=2^%.2f x=%.4f W=%d per=%d  fullstore %d/%d "
          "(model %.4f)  d5 %d d8 %d d11 %d  %.1fs"
          % (qmin, math.log2(U.classes), shape_real, W, per, acc["succ"], n,
             rec["fullstore"]["model_birthday"], acc["dp"][5]["succ"],
             acc["dp"][8]["succ"], acc["dp"][11]["succ"], 0.0), flush=True)
    return rec


def main():
    t_all = time.time()
    BT = int(os.environ.get("BT", "800"))
    W = int(os.environ.get("W", "16"))
    R = int(os.environ.get("RCOLS", "16384"))
    seed = int(os.environ.get("SEED", "20250"))
    qmin = int(os.environ.get("QMIN", "4290000000"))
    qmin2 = qmin // 255
    cells = []
    cells_sel = os.environ.get("CELLS", "design,low").split(",")
    if "design" in cells_sel:            # x = 1.6002, the budget the ledger charges
        run_cell(qmin, R, W, BT, seed, 1.0, cells)
    if "xlow" in cells_sel:              # model row p(x) = 1 - exp(-x^2/2)
        run_cell(qmin, R, W, BT, seed + 1, 0.80, cells)
    if "xhi" in cells_sel:
        run_cell(qmin, R, W, BT, seed + 2, 1.25, cells)
    if "low" in cells_sel:               # second order: 255x fewer classes, same x
        run_cell(qmin2, R, W, BT, seed + 3, 1.0, cells)
    out = {"test": "fullstore_census_v13",
           "protocol": ("shipped qwalk.step, escape armed kcap=8 (lag-1+lag-3); "
                        "EVERY visited class inserted and looked up (no report "
                        "gate on the success path); repeat at insert -> "
                        "x = (a1-a2)/(c2-c1) mod q recomputed and compared with "
                        "the trial challenge; gate_replay reads the SAME traces "
                        "through the v12 density gate as a matched ablation"),
           "ledger_shape": float(SHAPE), "cells": cells,
           "elapsed_s": round(time.time() - t_all, 1)}
    os.makedirs(os.path.join(D, "evidence"), exist_ok=True)
    with open(os.path.join(D, "evidence",
                           os.environ.get("OUT", "fullstore-census.json")), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("[done]", out["elapsed_s"], "s", flush=True)
    for c in cells:
        print("x=%.3f classes=2^%.1f fullstore %.4f [%.4f] model %.4f | d5 %.4f "
              "d8 %.4f d11 %.4f"
              % (c["x_visits_over_sqrtN"], c["classes_log2"],
                 c["fullstore"]["success_rate"], c["fullstore"]["cp95_success_f"],
                 c["fullstore"]["model_birthday"],
                 c["gate_replay"]["d5"]["success_rate"],
                 c["gate_replay"]["d8"]["success_rate"],
                 c["gate_replay"]["d11"]["success_rate"]), flush=True)


if __name__ == "__main__":
    main()
