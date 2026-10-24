#!/usr/bin/env python3
"""walker_slope.py -- v20: WALKER / TABLE / DENSITY arms at FIXED ORDER.

verdict-v18 point 1: the ladder moved only the order.  The reviewer asked for the
OTHER parameters moved -- selector, R, walker count, DP density -- because those
differ between the toy (R=64, W=16 walkers, d=5) and production (R=2^20, W=2^20
walkers, d=60).  This script moves R, W and d at ONE order (q = 2^28) with the
dimensionless shape x = visits/sqrt(N_classes) held at the shipped value, and
measures the same same-challenge w=2 statistic.

WHY THE SHAPE IS INVARIANT UNDER W.  per = x*sqrt(N)/W, so total visits per set
W*per = x*sqrt(N) is FIXED; changing W only changes how those visits are split
among walkers.  The work per set is therefore constant across the W arms and the
arms are directly comparable at equal cost.

CRASHED-CELL INSTRUMENT (v20).  A cell that writes p_set = 0 after a wall time far
below the numpy step rate is NOT a measurement: an instant exit leaves the success
mask empty and the census returns a silent zero.  Every arm here asserts

    per_walker_steps > 0, lanes > 0, walks_executed = sets*W*per > 0,
    elapsed_s >= walks_executed / 2e7        (numpy floor ~5e6 lane-steps/s),
    0 < p_set < 1, trials == sets

so a zero-try or instant-exit cell RAISES instead of writing a zero.  The five
cells v19 logged as "9600 sets in 4-6 s at p_set exactly 0" violate the time law
by 10-40x and are re-run here.

Arms (q = 2^28, N = 4.47e7 classes, d = 5 unless stated, sets as listed):
    W4    W=4    per=1688
    W16   W=16   per=422    (the shipped arm, reproduction)
    W64   W=64   per=106
    W256  W=256  per=26
    R512  R=512  W=16
    R1024 R=1024 W=16
    d7    W=16   d=7
    d9    W=16   d=9

Output: evidence/walker-arms.json with every cell plus the fitted logit slope in
log2(W) and in d, each with a parametric-bootstrap 95% CI.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import transfer_arm as T      # noqa: E402
import census_v20n as C       # noqa: E402

OUT = os.path.join(D, "evidence", "walker-arms.json")
Q28 = 1 << 28
ARMS = [
    ("W4",    dict(q_min=Q28, sets=2400, W=4,   d=5, R=64)),
    ("W16",   dict(q_min=Q28, sets=2400, W=16,  d=5, R=64)),
    ("W64",   dict(q_min=Q28, sets=2400, W=64,  d=5, R=64)),
    ("W256",  dict(q_min=Q28, sets=2400, W=256, d=5, R=64)),
    ("R512",  dict(q_min=Q28, sets=1600, W=16,  d=5, R=512)),
    ("R1024", dict(q_min=Q28, sets=1600, W=16,  d=5, R=1024)),
    ("d7",    dict(q_min=Q28, sets=1600, W=16,  d=7, R=64)),
    ("d9",    dict(q_min=Q28, sets=1600, W=16,  d=9, R=64)),
]


def check(out, spec):
    """Instrument: a silent zero / instant exit must DIE here, not reach evidence."""
    walks = out["sets"] * out["walks_per_set"] * out["per_walker_steps"]
    assert out["per_walker_steps"] > 0, ("zero-length trail", out["arm"])
    assert out["walks_per_set"] > 0 and out["sets"] > 0, ("empty cell", out["arm"])
    assert walks > 0, ("zero tries", out["arm"])
    assert out["runs"]["sets"] == spec["sets"], ("trial mismatch", out["arm"])
    assert out["elapsed_s"] >= walks / 2e7, (
        "instant exit: %s ran %d lane-steps in %.1fs" % (out["arm"], walks, out["elapsed_s"]))
    assert 0.0 < out["p_set_measured"] < 1.0, (
        "degenerate p_set %r in %s" % (out["p_set_measured"], out["arm"]))
    out["walks_executed"] = walks
    out["lane_steps_per_s"] = round(walks / max(out["elapsed_s"], 1e-9), 1)
    return out


def logit(p):
    return math.log(p / (1.0 - p))


def fit(xs, ks, ns, rng, draws=4000):
    """Weighted logit fit + parametric bootstrap CI (binomial resampling)."""
    def slope(kk):
        x = np.asarray(xs, dtype=float)
        y = np.array([logit(k / n) for k, n in zip(kk, ns)])
        w = np.array([k * (1 - k / n) for k, n in zip(kk, ns)], dtype=float)  # binomial info
        A = np.vstack([np.ones_like(x), x]).T
        Wm = np.diag(w)
        beta = np.linalg.solve(A.T @ Wm @ A, A.T @ Wm @ y)
        return float(beta[1])

    b0 = slope(ks)
    boot = []
    for _ in range(draws):
        kk = [int(rng.binomial(n, k / n)) for k, n in zip(ks, ns)]
        try:
            boot.append(slope(kk))
        except np.linalg.LinAlgError:
            continue
    lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"slope": b0, "ci95": [float(lo), float(hi)], "draws": len(boot),
            "method": "weighted logit OLS, binomial parametric bootstrap"}


def main():
    t0 = time.time()
    rows = []
    for name, spec in ARMS:
        out = T.run_arm(name, spec["q_min"], spec["sets"], R=spec["R"], d=spec["d"],
                        W=spec["W"], seed=20250)
        rows.append(check(out, spec))
    # walker slope at fixed d=5, R=64
    wsel = [(r["walks_per_set"], r["runs"]["set_solved"], r["sets"])
            for r in rows if r["dp_bits"] == 5 and r["R"] == 64]
    wsel.sort()
    rng = np.random.default_rng(7)
    wfit = fit([math.log2(w) for w, _k, _n in wsel],
               [k for _w, k, _n in wsel], [n for _w, _k, n in wsel], rng)
    # density slope at fixed W=16, R=64
    dsel = sorted((r["dp_bits"], r["runs"]["set_solved"], r["sets"])
                  for r in rows if r["walks_per_set"] == 16 and r["R"] == 64)
    dfit = fit([d for d, _k, _n in dsel], [k for _d, k, _n in dsel],
               [n for _d, _k, n in dsel], rng)
    # table slope (R 64 -> 1024) at W=16, d=5
    rsel = sorted((r["R"], r["runs"]["set_solved"], r["sets"])
                  for r in rows if r["walks_per_set"] == 16 and r["dp_bits"] == 5)
    rfit = fit([math.log2(x) for x, _k, _n in rsel], [k for _x, k, _n in rsel],
               [n for _x, _k, n in rsel], rng)
    out = {
        "test": "walker_table_density_arms_fixed_order",
        "q_log2": 28, "order_q": None, "sets": {r["arm"]: r["sets"] for r in rows},
        "shape_visits_over_sqrtN": {r["arm"]: r["visits_over_sqrtN"] for r in rows},
        "cells": rows,
        "slope_walker_logit_per_log2W": wfit,
        "slope_density_logit_per_bit": dfit,
        "slope_table_logit_per_log2R": rfit,
        "walker_slope_negative": bool(wfit["slope"] < 0),
        "walker_slope_sig_negative": bool(wfit["ci95"][1] < 0),
        "charge_log2W_16_to_20": 16.0,
        "elapsed_s": round(time.time() - t0, 1),
        "note": ("all cells pass the crash instrument (nonzero tries, elapsed >= "
                 "lane-steps/2e7); p_set is the same-challenge w=2 pairing statistic"),
    }
    out["order_q"] = rows[0]["q"]
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    for r in rows:
        print("%-6s R=%-5d W=%-4d d=%d per=%-5d sets=%-5d p_set %.5f cp95 %.5f "
              "%.1fs %s" % (r["arm"], r["R"], r["walks_per_set"], r["dp_bits"],
                            r["per_walker_steps"], r["sets"], r["p_set_measured"],
                            r["cp95_set_lower_f"], r["elapsed_s"],
                            "OK" if r["floor_met_by_cp95"] else "BELOW"), flush=True)
    print("slope_walker %.4f CI %s" % (wfit["slope"], wfit["ci95"]))
    print("slope_density %.4f CI %s" % (dfit["slope"], dfit["ci95"]))
    print("slope_table %.4f CI %s" % (rfit["slope"], rfit["ci95"]))


if __name__ == "__main__":
    main()
