"""Planner pair-fraction scaling evidence for the H4 pair-fraction claim.

Runs the shipped evidence.py planner (plan/class_lists/toy_prime, unmodified)
on a predeclared ladder of (d, X2, seed) specs, prints the measured
roots-per-element ladder, fits the window scaling, and checks the production
extrapolation margins for f=0.985 (direct) and f=0.99 (conservative half-fit).

Standard library only. Runtime is a few minutes (dominated by the three
largest specs). Not part of the hosted experiment suite; run on demand:

    python3 planner_f99.py

Production parameters for comparison: d = 2^20, X2 = 2^76, class lists of
about 2^56 elements, planner window 2 GMAX/d = 2^25 members, span/(2 GMAX)
about 2^31 windows. Every toy below is far smaller; the claim is an explicit
heuristic extrapolation, not a proof (see proof.md section 1).
"""

import json
import math
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import evidence as E

# (d, logX2, seed, classes). GMAX = 16 d^2 and the toy prime follow
# evidence.planner_mode exactly. List covers windows 256..4096 members and
# spans 4..512 windows. Specs with span below 8 windows (list shorter than a
# few windows, opposite of the production regime) are measured for context but
# excluded from every fit and bound (see PATHOLOGICAL below).
SPECS = [
    (8, 13, 0, 8),
    (8, 14, 0, 8),
    (8, 17, 0, 4),
    (8, 17, 1, 4),
    (16, 15, 0, 8),
    (16, 15, 1, 8),
    (16, 15, 2, 8),
    (16, 18, 0, 4),
    (16, 20, 0, 4),
    (16, 22, 0, 2),
    (32, 21, 0, 4),
    (32, 21, 1, 4),
    (64, 20, 0, 4),
    (64, 23, 0, 2),
    (128, 23, 0, 2),
]

# Production planner window in class-list members: 2*GMAX/d at d = 2^20.
W_PROD = 1 << 25


def run_spec(dd, logX2, seed, classes):
    rng = random.Random(seed)
    rr = E.toy_prime(28 + int(math.log2(dd)) + 4, dd, rng)
    X2 = 1 << logX2
    cl = E.class_lists(rr, dd, X2)
    GMAX = 16 * dd * dd
    window = 2 * GMAX // dd
    span = round(X2 / (2 * GMAX), 1)
    rows = []
    for c in range(min(classes, dd)):
        L = cl[c]
        if len(L) < 50:
            continue
        _, st = E.plan(L, GMAX, max_tree=1 << 20)
        rows.append(st)
    ns = [st["n"] for st in rows]
    pfs = [st["pair_fraction"] for st in rows]
    return {
        "d": dd, "X2": X2, "GMAX": GMAX, "window_members": window,
        "span_over_2GMAX": span, "r_bits": rr.bit_length(),
        "classes": len(rows), "n_total": sum(ns),
        "n_min": min(ns), "n_max": max(ns),
        "min_pair_fraction": min(pfs), "mean_pair_fraction": sum(pfs) / len(pfs),
        "mean_roots_per_element": sum(st["roots"] for st in rows) / sum(ns),
        "roots_per_class": [st["roots"] for st in rows],
        "max_cycle_rejections": max(st["cycle_rejections"] for st in rows),
        "size_rejections": sum(st["size_rejections"] for st in rows),
        "depth_rejections": sum(st.get("depth_rejections", 0) for st in rows),
        "demoted": sum(st["demoted"] for st in rows),
        "max_depth_seen": max(st["max_depth"] for st in rows),
        "max_tree_size_seen": max(st["max_tree_size"] for st in rows),
    }


def main():
    t0 = time.time()
    specs = [run_spec(*s) for s in SPECS]
    # Non-pathological rows: span of at least 8 windows (production: 2^31).
    good = [r for r in specs if r["span_over_2GMAX"] >= 8]
    # Fit ladder: largest span per window, span >= 16, d >= 16.
    ladder = {}
    for r in good:
        w = r["window_members"]
        if r["d"] >= 16 and r["span_over_2GMAX"] >= 16:
            if w not in ladder or r["X2"] > ladder[w]["X2"]:
                ladder[w] = r
    ws = sorted(ladder)
    xs = [math.log(w) for w in ws]
    ys = [math.log(ladder[w]["mean_roots_per_element"]) for w in ws]
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    beta = -sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
    beta_c = beta / 2  # conservative: half the fitted exponent
    ref = ladder[max(ws)]
    r_ref = ref["mean_roots_per_element"]
    w_ref = ref["window_members"]
    r_prod = r_ref * (W_PROD / w_ref) ** (-beta_c)
    f_prod = 1 - r_prod
    mono = all(
        ladder[ws[i]]["mean_roots_per_element"]
        > ladder[ws[i + 1]]["mean_roots_per_element"]
        for i in range(len(ws) - 1)
    )
    checks = {
        # (a) sanity: every non-pathological row is well behaved.
        "all_good_rows_roots_le_2pt1pct": all(
            r["mean_roots_per_element"] <= 0.021 for r in good
        ),
        "good_max_roots": max(
            r["mean_roots_per_element"] for r in good
        ),
        # (b) direct support for f = 0.985 at d >= 16, span >= 16.
        "ladder_rows_roots_le_1pt5pct": all(
            ladder[w]["mean_roots_per_element"] <= 0.015 for w in ws
        ),
        "ladder_max_roots": max(
            ladder[w]["mean_roots_per_element"] for w in ws
        ),
        # (c) conservative half-fit extrapolation supports f = 0.99.
        "beta_fit": beta,
        "beta_conservative": beta_c,
        "reference": {"window": w_ref, "roots_per_element": r_ref},
        "extrapolated_production_roots": r_prod,
        "extrapolated_production_f": f_prod,
        "extrapolation_margin_to_1pct": 0.01 / r_prod,
        "extrapolated_roots_le_1pct": r_prod <= 0.01,
        # (d) strictly decreasing roots/n in window on the ladder.
        "window_monotone_decreasing": mono,
        # (e) planner capacity: no binding constraint in any toy.
        "size_rejections_total": sum(r["size_rejections"] for r in specs),
        "depth_rejections_total": sum(r["depth_rejections"] for r in specs),
        "max_depth_any": max(r["max_depth_seen"] for r in specs),
        "max_tree_any": max(r["max_tree_size_seen"] for r in specs),
    }
    assert checks["all_good_rows_roots_le_2pt1pct"]
    assert checks["ladder_rows_roots_le_1pt5pct"]
    assert checks["extrapolated_roots_le_1pct"], checks
    assert checks["window_monotone_decreasing"]
    assert checks["size_rejections_total"] == 0
    assert checks["depth_rejections_total"] == 0
    assert checks["max_depth_any"] <= 1 << 14
    assert checks["max_tree_any"] <= 1 << 20
    out = {"specs": specs, "checks": checks,
           "elapsed_seconds": round(time.time() - t0, 1)}
    print(json.dumps(out, indent=1))
    print("---- ladder (window members, span, mean roots/element) ----")
    for r in specs:
        print(
            "d=%d X2=2^%d span=%.1f n=%d window=%d min_f=%.4f mean_f=%.4f "
            "roots/n=%.4f rejections(c/s/d)=%d/%d/%d demoted=%d depth=%d "
            "tree=%d" % (
                r["d"], int(math.log2(r["X2"])), r["span_over_2GMAX"],
                r["n_total"], r["window_members"], r["min_pair_fraction"],
                r["mean_pair_fraction"], r["mean_roots_per_element"],
                r["max_cycle_rejections"], r["size_rejections"],
                r["depth_rejections"], r["demoted"], r["max_depth_seen"],
                r["max_tree_size_seen"]))
    print("fit: beta=%.4f beta_c=%.4f ref_w=%d ref_r=%.4f r_prod=%.4f "
          "f_prod=%.4f margin_to_1pct=%.2f" % (
              beta, beta_c, w_ref, r_ref, r_prod, f_prod,
              0.01 / r_prod))


if __name__ == "__main__":
    main()
