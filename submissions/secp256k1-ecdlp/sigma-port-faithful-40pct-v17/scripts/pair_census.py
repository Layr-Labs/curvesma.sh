#!/usr/bin/env python3
"""pair_census.py -- MEASURED w=2 amplification (two sequential fresh runs per set).

Why this file exists.  The single-run success of the shipped walk is ~0.283, below
the track's 0.39 floor.  The bls lineage's valid ticket rode a *measured* w=2
amplification instead of a model: run the walk twice, independently, on the SAME
challenge, and count the pairing as solved when EITHER run solves.  This script
measures that on secp256k1's census machinery, with no birthday-model premise:

    p_set = P(run A solves  OR  run B solves)       with A, B independent

so p_set = 1 - (1 - p)^2 by independence, and the SHIPPED number is the measured
set fraction's one-sided 95% Clopper-Pearson lower bound -- not the formula.  The
formula is reported beside it and the agreement row is |measured - 1-(1-p)^2|.

Independence is constructed, not assumed: run B draws a DIFFERENT fixed challenge
block and a different walker coin stream (walk_arm coins_seed), so no state, no
start and no coin is shared with run A; the census' own fixed-seed block design
is preserved inside each run.

The charge for this is +1 bit to W (two runs), disclosed in ledger.py's W2 row.

CP helper calibration (mandatory): cp95_lower is Monte-Carlo checked here at the
observed n and p -- the fraction of simulated Binomial(n, p) draws whose CP lower
bound is <= p must be ~0.95 (one-sided coverage), plus the k == n closed form.
"""
import json
import math
import os
import sys
from decimal import Decimal

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk               # noqa: E402
import census_v20n as C    # noqa: E402


def calibrate(n, p_hat, sims=2000, seed=7):
    """Monte-Carlo coverage of C.cp95_lower + the k==n closed form."""
    rng = np.random.default_rng(seed)
    ks = rng.binomial(n, p_hat, size=sims)
    cache = {}
    covered = 0
    for k in ks.tolist():
        if k not in cache:
            cache[k] = float(C.cp95_lower(k, n))
        covered += int(cache[k] <= p_hat)
    closed = float(Decimal("0.05") ** (Decimal(1) / Decimal(n)))
    exact_n = float(C.cp95_lower(n, n))
    return {"sims": sims, "n": n, "p": p_hat,
            "coverage": covered / sims,
            "k_equals_n_closed_form": closed,
            "k_equals_n_helper": exact_n,
            "k_equals_n_abs_err": abs(closed - exact_n),
            "k_zero_helper": float(C.cp95_lower(0, n))}


def main():
    arm = os.environ.get("ARM", "lag3")
    lag = int(arm[3:])
    BT = int(os.environ.get("BT", "3200"))
    SEED = int(os.environ.get("SEED", "20250"))
    W, d = 16, 5
    U = qwalk.Universe(4_290_000_000, 64, seed=SEED)
    per = int(round(float(C.SHAPE) * math.sqrt(U.classes) / W))
    out = {"test": "pair_census_w2", "arm": arm, "sets": BT, "walks_per_set": W,
           "per_walker_steps": per, "dp_bits": d, "universe_q": int(U.q),
           "classes": int(U.classes), "seed": SEED,
           "interpreter": sys.executable, "numpy": np.__version__,
           "protocol": ("two sequential independent runs per set: different fixed "
                        "challenge block (seed*1000+r) and different walker coin "
                        "stream (0xC0FFEE+r); a set is solved when either run "
                        "solves; escape armed, kcap=8")}
    solved = []
    fires3 = inv_bad = 0
    for r in (0, 1):
        xt, xv, a_all, cc_all, canon_all = C.make_start(U, BT, W, SEED * 1000 + r + 1)
        res = C.walk_arm(U, canon_all, a_all, cc_all, xv, d, per, True, lag, W, BT,
                         coins_seed=0xC0FFEE + r)
        solved.append(res["succ_step"] >= 0)
        fires3 += res["fires_lag3"]
        inv_bad += res["invariant_violations"]
    A, B = solved
    s_set = (A | B)
    n = BT
    kA, kB, kS = int(A.sum()), int(B.sum()), int(s_set.sum())
    pA, pB, pS = kA / n, kB / n, kS / n
    formula = 1 - (1 - pA) * (1 - pB)
    out["runs"] = {"A_solved": kA, "B_solved": kB,
                   "both": int((A & B).sum()), "neither": int((~A & ~B).sum()),
                   "set_solved": kS, "sets": n}
    out["p_run"] = {"A": pA, "B": pB, "pooled": (kA + kB) / (2 * n)}
    out["p_set_measured"] = pS
    out["p_set_formula"] = formula
    out["formula_minus_measured"] = formula - pS
    out["cp95_set_lower"] = str(C.cp95_lower(kS, n))
    out["cp95_set_lower_f"] = float(C.cp95_lower(kS, n))
    out["cp95_formula_lower"] = str(C.cp95_lower(int(round(formula * n)), n))
    out["cp95_run_A_lower"] = str(C.cp95_lower(kA, n))
    out["floor"] = 0.39
    out["floor_met_by_cp95"] = bool(float(C.cp95_lower(kS, n)) >= 0.39)
    out["lag3_fires"] = fires3
    out["invariant_violations_sampled"] = inv_bad
    out["calibration"] = calibrate(n, pS)
    # does the *single-run* lower bound also clear the floor?  (it should not)
    out["single_run_floor_met"] = bool(float(C.cp95_lower(kA, n)) >= 0.39)
    os.makedirs(os.path.join(D, "evidence"), exist_ok=True)
    with open(os.path.join(D, "evidence", os.environ.get("OUT", "pair-census.json")), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
