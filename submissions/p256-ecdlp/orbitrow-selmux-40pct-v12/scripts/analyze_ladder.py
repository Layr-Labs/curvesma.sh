#!/usr/bin/env python3
"""analyze_ladder.py -- slopes + the production-extrapolation budget table.

Method: secp v19/v20 (sign-honest charge).  For each swept parameter fit
log2(odds) = a + slope * log2(param) on the three measured cells (OLS, SE with
1 dof, +/-1.96 SE reported), then extrapolate to the production value:

    charge_bits = |slope| * dlog2   if slope < 0   (UNFAVOURABLE: charged in full)
                = 0                if slope >= 0  (favourable: charged ZERO)

Total charge is compared against the floor slack log2(success/0.39) -- the only
success-odds budget the claim has.  If the charge exceeds the slack the ladder says
the transfer is not supported and the ledger must say so.
"""
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
EV = os.path.join(D, "evidence")

PROD = {"W": 2 ** 20, "R": 2 ** 20, "dp_bits": 60}


def odds(p):
    p = min(max(p, 1e-9), 1 - 1e-9)
    return p / (1 - p)


def fit(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx
    a = my - b * mx
    if n > 2:
        ss = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
        se = math.sqrt(ss / (n - 2) / sxx)
    else:
        se = float("nan")
    return a, b, se


def main():
    merged = os.path.join(EV, "ladder-cells.json")
    if os.path.exists(merged):
        cells = json.load(open(merged))["cells"]
    else:
        cells = {}
        for f in sorted(glob.glob(os.path.join(EV, "ladder-*.json"))):
            d = json.load(open(f))
            cells[d["cell"]] = d
    need = ["base", "W64", "W256", "R4096", "R1024", "d8", "d11", "selval"]
    missing = [c for c in need if c not in cells]
    assert not missing, ("ladder cells missing: %s (a crashed cell is not evidence)" % missing)
    base = cells["base"]
    for c in cells.values():
        assert c["runtime_guard_ok"], ("runtime guard failed: %s" % c["cell"])
        assert c["invariant_violations"] == 0, c
    rows = []
    budget = {}
    for dim, keys, field in (("walkers", ["base", "W64", "W256"], "W"),
                             ("columns", ["base", "R4096", "R1024"], "R"),
                             ("DP density", ["base", "d8", "d11"], "dp_bits")):
        xs = [math.log2(cells[k][field]) for k in keys]
        ys = [math.log2(odds(cells[k]["success_rate"])) for k in keys]
        a, b, se = fit(xs, ys)
        dlog = math.log2(PROD[field]) - xs[0]          # base -> production
        charge = max(0.0, -b) * dlog
        ci = 1.96 * se if se == se else 0.0
        budget[dim] = {"param": field, "base": cells[keys[0]][field], "production": PROD[field],
                       "delta_log2": dlog, "slope_bits_per_log2": b,
                       "slope_ci95": ci, "sign": "favourable" if b > 0 else "unfavourable",
                       "charge_bits": charge,
                       "charge_bits_worst_ci": (max(0.0, -(b - ci)) * dlog)}
        rows.append({"dim": dim, "points": [{"cell": k, "param": cells[k][field],
                                            "p": cells[k]["success_rate"],
                                            "cp95": cells[k]["cp95_success_f"]} for k in keys],
                     "slope": b, "slope_ci95": ci, "slope_se": se, "n": len(keys),
                     "charge_bits": charge})
    total = sum(v["charge_bits"] for v in budget.values())
    total_worst = sum(v["charge_bits_worst_ci"] for v in budget.values())
    sel = cells["selval"]
    swap = {"cell": "selector-swap", "uniform_key": {"p": base["success_rate"],
                                                     "fires_per_lane_step": base["fires_per_lane_step"]},
            "class_value": {"p": sel["success_rate"],
                            "fires_per_lane_step": sel["fires_per_lane_step"]},
            "delta_p": base["success_rate"] - sel["success_rate"],
            "fire_ratio": sel["fires_per_lane_step"] / base["fires_per_lane_step"],
            "note": ("class-value truncation is the measured WRONG model (v20g): the "
                     "canonical class is concentrated near 0, so equal-width buckets "
                     "collide far more often.  The shipped selector uses a uniform key "
                     "and the row is the justification, not a claim input.")}
    slack = math.log2(base["success_rate"] / 0.39)
    out = {"test": "sensitivity_ladder_v12", "method": "secp v19/v20 sign-honest ladder",
           "sets_per_cell": base["sets"], "shape_x": base["visits_over_sqrtN"],
           "base_success": base["success_rate"], "floor": 0.39,
           "floor_slack_bits": slack, "dims": rows, "budget": budget,
           "total_charge_bits": total, "total_charge_bits_worst_ci": total_worst,
           "transfer_supported": bool(total <= slack),
           "selector_swap": swap,
           "cells": {k: {kk: v[kk] for kk in ("W", "R", "dp_bits", "successes", "sets",
                                             "success_rate", "cp95_success_f",
                                             "fires_per_lane_step", "elapsed_s")}
                     for k, v in cells.items()}}
    with open(os.path.join(EV, "ladder-analysis.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    for r in rows:
        print("%-10s slope %+7.4f +/- %.4f bits/log2  -> charge %.4f bits"
              % (r["dim"], r["slope"], r["slope_ci95"], r["charge_bits"]))
    print("selector swap: delta_p %+.4f  fire ratio %.2fx (class-value is worse)"
          % (swap["delta_p"], swap["fire_ratio"]))
    print("total charge %.4f bits (worst CI %.4f) vs floor slack %.4f bits -> %s"
          % (total, total_worst, slack, "SUPPORTED" if total <= slack else "NOT SUPPORTED"))


if __name__ == "__main__":
    main()
