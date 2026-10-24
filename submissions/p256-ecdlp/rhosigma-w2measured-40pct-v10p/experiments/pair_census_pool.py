#!/usr/bin/env python3
"""pair_census_pool.py -- pool the w=2 pair-census shards, CP95 lower, calibration.

Protocol per set: ONE challenge x, run 1 with fresh walker coins; if run 1 fails, run 2
with FRESH walker coins on the SAME x; the set is solved iff either run solves.  A shard
prints one JSON line {trials, successes (sets), run1_solved, run2_solved, ...}.

CP one-sided 95% lower bound: p_L = Beta^{-1}(0.05; k, n-k+1).  The helper is CALIBRATED
before use (Monte-Carlo coverage at the observed n and p, plus the k==n closed form
0.05^(1/n)), as the fleet law requires.
"""
import glob
import json
import math
import os
import sys
from decimal import Decimal

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def betainc_inv_05(a, b):
    """Solve P(X<=x)=0.05 for X~Beta(a,b): the one-sided 95 pct LOWER Clopper-Pearson limit.

    The bound must sit BELOW the point estimate, so the target is the 5th percentile of
    Beta(k, n-k+1), not the 95th (solving for 0.95 returns the UPPER limit - that inversion
    bug was caught by asserting cp95_lower <= k/n on every row).
    """
    from math import lgamma

    def betacf(x, a, b, MAXIT=300, EPS=3e-16, FPMIN=1e-300):
        qab, qap, qam = a + b, a + 1.0, a - 1.0
        c = 1.0
        d = 1.0 - qab * x / qap
        if abs(d) < FPMIN:
            d = FPMIN
        d = 1.0 / d
        h = d
        for m in range(1, MAXIT + 1):
            m2 = 2 * m
            aa = m * (b - m) * x / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d
            if abs(d) < FPMIN:
                d = FPMIN
            c = 1.0 + aa / c
            if abs(c) < FPMIN:
                c = FPMIN
            d = 1.0 / d
            h *= d * c
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d
            if abs(d) < FPMIN:
                d = FPMIN
            c = 1.0 + aa / c
            if abs(c) < FPMIN:
                c = FPMIN
            d = 1.0 / d
            de = d * c
            h *= de
            if abs(de - 1.0) < EPS:
                break
        return h

    def betai(x, a, b):
        if x <= 0.0:
            return 0.0
        if x >= 1.0:
            return 1.0
        lbeta = lgamma(a + b) - lgamma(a) - lgamma(b)
        front = math.exp(lbeta + a * math.log(x) + b * math.log1p(-x))
        if x < (a + 1.0) / (a + b + 2.0):
            return front * betacf(x, a, b) / a
        return 1.0 - math.exp(lbeta + b * math.log1p(-x) + a * math.log(x)) * betacf(1.0 - x, b, a) / b

    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if betai(mid, a, b) < 0.05:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def cp95_lower(k, n):
    if n <= 0:
        return 0.0
    if k == 0:
        return 0.0
    if k == n:
        return float(Decimal("0.05") ** (Decimal(1) / Decimal(n)))
    v = betainc_inv_05(k, n - k + 1)
    assert v <= k / n, (k, n, v)          # a lower bound that exceeds p_hat is a bug, not a result
    return v


def calibrate(n, p, sims=4000, seed=11):
    rng = np.random.default_rng(seed)
    ks = rng.binomial(n, p, size=sims)
    cache, covered = {}, 0
    for k in ks.tolist():
        if k not in cache:
            cache[k] = cp95_lower(k, n)
        covered += int(cache[k] <= p)
    closed = float(Decimal("0.05") ** (Decimal(1) / Decimal(n)))
    return {"sims": sims, "n": n, "p": p, "coverage_one_sided": covered / sims,
            "k_equals_n_closed_form": closed,
            "k_equals_n_helper": cp95_lower(n, n),
            "k_equals_n_abs_err": abs(closed - cp95_lower(n, n)),
            "k_zero": cp95_lower(0, n), "k_half": cp95_lower(n // 2, n)}


def pool(pattern):
    trials = suc = r1 = r2 = recs = meets = ver = dbl = ring = 0
    shards = []
    for f in sorted(glob.glob(pattern)):
        try:
            d = json.load(open(f))
        except Exception:
            continue
        shards.append(os.path.basename(f))
        trials += d["trials"]; suc += d["successes"]
        r1 += d["run1_solved"]; r2 += d["run2_solved"]
        recs += d["records"]; meets += d["meets"]; ver += d["verified"]
        dbl += d["doubles"]; ring += d["ring_hits"]
    return {"shards": len(shards), "trials": trials, "sets_solved": suc,
            "run1_solved": r1, "run2_solved": r2, "records": recs,
            "meets": meets, "verified": ver, "doubles": dbl, "ring_hits": ring,
            "shard_files": shards}


if __name__ == "__main__":
    root = sys.argv[1] if len(sys.argv) > 1 else "/tmp/pv10p/out"
    out = {}
    for arm in ("shipped_guarded", "control_pre_parity"):
        for av in ("0", "-3"):
            key = f"{arm}|a={av}"
            agg = pool(os.path.join(root, f"{'guarded' if arm=='shipped_guarded' else 'control'}_a{av}_s*.json"))
            n, k = agg["trials"], agg["sets_solved"]
            k1 = agg["run1_solved"]
            p2 = k / n if n else float("nan")
            p1 = k1 / n if n else float("nan")
            formula = 1 - (1 - p1) ** 2 if n else float("nan")
            agg.update({"p1_single_run": p1, "p2_set": p2,
                        "p2_cp95_lower": cp95_lower(k, n),
                        "p1_cp95_lower": cp95_lower(k1, n),
                        "formula_1_minus_1mp1_sq": formula,
                        "measured_minus_formula": p2 - formula,
                        "run2_given_run1_failed": (k - k1) / (n - k1) if n > k1 else None,
                        "floor": 0.39, "floor_met": cp95_lower(k, n) >= 0.39,
                        "w2_threshold": 1 - math.sqrt(0.61)})
            out[key] = agg
    cal_n_big = out["shipped_guarded|a=0"]["trials"]
    out["cp_calibration"] = {
        "shipped_arm": calibrate(out["shipped_guarded|a=0"]["trials"],
                                 out["shipped_guarded|a=0"]["p2_set"]),
        "control_arm": calibrate(out["control_pre_parity|a=0"]["trials"],
                                 out["control_pre_parity|a=0"]["p2_set"]),
        "n_sets_shipped": cal_n_big,
        "helper": "CP one-sided 95% lower = Beta^-1(0.05; k, n-k+1) (bisection on I_x(a,b))",
    }
    p = os.path.join(HERE, "pair_census_v10p.json")
    json.dump(out, open(p, "w"), indent=1, sort_keys=True)
    for k, v in out.items():
        if k == "cp_calibration":
            continue
        print(k, "n=%d k=%d p1=%.5f p2=%.5f cplow=%.5f floor_met=%s" % (
            v["trials"], v["sets_solved"], v["p1_single_run"], v["p2_set"],
            v["p2_cp95_lower"], v["floor_met"]))
    print("wrote", p)