#!/usr/bin/env python3
"""budget2.py -- v20k two-phase budget: the DECLARED success row, per arm.

THE DEFECT THIS REPAIRS (reviewer point 2 of bb309448).  v20j's headline success
was the drain-inclusive rate while the ledger priced only the walk, so the
declared success and the priced budget were two different experiments.  v20k
declares the budget TWO-PHASE up front

    B_total = B_walk + B_drain,        both phases charged in one ledger,

and reports, for every arm, the success measured at BOTH budgets:

    success_walk_only = recoveries completed inside B_walk
    success_total     = recoveries completed inside B_total  <- THE DESIGN ROW

so the row that gates is the row measured at the work the ledger charges.  The
walk-only row is published as the disclosed sensitivity (it is what the same
experiment would report if the drain were removed; it is NOT the production
protocol, because a production walker's drain is charged work that it executes).

CP-95 is the one-sided Clopper-Pearson lower bound, computed by the Decimal
bisection copied verbatim from
curvesmash/submissions/secp256k1-ecdlp/sigma-table-40pct-v15k/scripts/gen_cert_v15k.py
(`cp95_lower`): the same helper whose 0.398285 row certified on that track.  The
bisection bracket is widened from [0.30,0.45] to [0.05,0.95] because two arms sit
below 0.30; the result is ASSERTED to lie strictly inside the bracket so a
clipped root can never be reported as a bound.
"""
import json
import os
import sys
from decimal import Decimal, getcontext

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
SRC = os.path.join(D, "evidence", "prod-sweep.json")
OUT = os.path.join(D, "evidence", "budget2.json")

getcontext().prec = 60
alpha = Decimal("0.05")
LO, HI = Decimal("0.001"), Decimal("0.95")


def cp95_lower(k, n):
    """one-sided Clopper-Pearson lower bound: smallest p with P[X>=k | p,n] >= alpha.
    P[X>=k] is strictly INCREASING in p, so bisection converges on the lower end.
    (body verbatim from gen_cert_v15k.py; bracket widened, clipping asserted)"""
    k = int(k)
    n = int(n)
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p):
        p = Decimal(p)
        q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()   # term_k
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal("1e-40") * s:
                break
        return s

    lo, hi = LO, HI
    for _ in range(200):
        mid = (lo + hi) / 2
        if tail(mid) > alpha:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-13"):
            break
    assert LO < lo < HI, ("cp95 bracket clipped", lo, k, n)
    return lo


def main():
    d = json.load(open(SRC))
    rows = []
    for r in d["rows"]:
        n = int(r["trials"])
        k_walk = int(r["successes_budget_raw"])
        k_tot = int(r["successes_raw"])
        assert k_walk + int(r["successes_drain_raw"]) == k_tot, r["dp_bits"]
        rows.append({
            "dp_bits": int(r["dp_bits"]),
            "trials": n,
            "walks": int(r["walks"]),
            "budget_steps": int(r["budget_steps"]),
            "drain_cap_steps": int(r["drain_cap_steps"]),
            "drain_window_over_budget": r["drain_window_over_budget"],
            "spacings_per_walker": r["budget_spacings_per_walker"],
            "successes_walk_only": k_walk,
            "successes_total": k_tot,
            "success_walk_only": k_walk / n,
            "success_total": k_tot / n,
            "cp95_walk_only": str(cp95_lower(k_walk, n)),
            "cp95_total": str(cp95_lower(k_tot, n)),
            "cp95_walk_only_f": float(cp95_lower(k_walk, n)),
            "cp95_total_f": float(cp95_lower(k_tot, n)),
        })
    # DESIGN ARM.  Production runs 2^46 DP spacings inside one walker budget; the
    # arm with the most spacings per walker is the closest to that regime, and it
    # is also the arm with the worst CP-95 on BOTH rows, so declaring it design
    # cannot hide a worse arm.  The other arms are reported in full.
    design = max(rows, key=lambda r: r["spacings_per_walker"])
    inproc = [r for r in rows if r["dp_bits"] in d["gate"]["gated_densities"]]
    worst_tot = min(inproc, key=lambda r: r["cp95_total_f"])
    worst_walk = min(inproc, key=lambda r: r["cp95_walk_only_f"])
    assert design is worst_tot, (
        "design arm must be the worst in-protocol arm on the DESIGN row",
        design["dp_bits"], worst_tot["dp_bits"])
    assert len(inproc) == 3, [r["dp_bits"] for r in inproc]
    FLOOR = 0.39
    bud = {
        "test": "two_phase_budget",
        "source": "evidence/prod-sweep.json",
        "design_row": "success_total",
        "sensitivity_row": "success_walk_only",
        "why_design_row": ("the ledger charges B_total = B_walk + B_drain as one "
                           "budget, so the success statistic must be the one "
                           "measured over B_total; success_walk_only is the same "
                           "experiment with the charged drain removed and is "
                           "published as the disclosed sensitivity"),
        "design_arm_dp_bits": design["dp_bits"],
        "design_arm_why": ("most DP spacings per walker (closest to the 2^46 "
                           "production regime) and the worst CP-95 on both rows"),
        "floor": FLOOR,
        "cp95_total_design": design["cp95_total"],
        "cp95_walk_only_design": design["cp95_walk_only"],
        "floor_met_total": design["cp95_total_f"] >= FLOOR,
        "floor_met_walk_only": design["cp95_walk_only_f"] >= FLOOR,
        "arms": rows,
    }
    with open(OUT, "w") as f:
        json.dump(bud, f, sort_keys=True, indent=1)
    for r in rows:
        print("d=%-2d  walk_only %.6f (CP95 %.6f)   total %.6f (CP95 %.6f)"
              "   drain %.3f of walk"
              % (r["dp_bits"], r["success_walk_only"], r["cp95_walk_only_f"],
                 r["success_total"], r["cp95_total_f"],
                 r["drain_window_over_budget"]))
    print("design arm d=%d  CP95_total %s >= 0.39: %s"
          % (design["dp_bits"], design["cp95_total"], bud["floor_met_total"]))
    print("                CP95_walk_only %s >= 0.39: %s"
          % (design["cp95_walk_only"], bud["floor_met_walk_only"]))


if __name__ == "__main__":
    main()
