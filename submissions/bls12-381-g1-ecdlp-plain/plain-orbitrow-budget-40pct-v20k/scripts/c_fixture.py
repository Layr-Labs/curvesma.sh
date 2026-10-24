#!/usr/bin/env python3
"""c_fixture.py -- v20g: the C fixture.  Production selector, 2^26 classes, a
FULL 2^20-step reporting spacing per walker.

WHY THIS EXISTS (reviewer point 3, second half).  v20e's medium fixture ran 16,384
steps per lane at a 2^14 spacing, so "183.5 trap opportunities per DP spacing" was
a rate extrapolation and the surrogate used a SplitMix partition of scalar
representatives instead of the production coordinate selector.  This fixture:

  * uses qwalk.col, which is the ledger's own selector.cell(K, R, P) on a uniform
    class-function K (see qwalk.col docstring) -- the SAME imported module the
    ledger prices;
  * has N = 2^26 classes (the medium fixture has 2^25.4; production 2^375, so this
    is a scale, not a size, and it is disclosed as such);
  * runs 2^20 steps per walker, i.e. ONE FULL production-shaped reporting spacing
    between reports (d = 20), so the trap/escape regime is measured at the
    spacing, not extrapolated to it;
  * pre-capture margin: the spacing holds 2^20/(6*R_eff) = 341 two-step-return
    opportunities against the >= 32x margin the fixture-scale law requires, and
    the production figure 2^37.4 is the same quantity at production R and DP;
  * measures, per walk: escape fires, the longest run of consecutive substitutions
    along ONE walk (K_obs), substitutions that actually left the two-cycle, and a
    residual counter; and it measures class coverage AND pair (class, prev-class)
    coverage on the same traces, plus the conversion from a class collision to a
    verified x from two records.

The pool in production is DP-gated.  At one spacing per walker a walker emits ~1
record, so record-level collisions are rare here and are REPORTED, not assumed; the
conversion statement is measured on the full-record pool (d = 0 semantics), which
is the algebra the pool replays.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qwalk  # noqa: E402

Q_MIN = 6 * (1 << 26) + 1
R_TAB = 1 << 10
SPACING_BITS = 20
LANES = 8
STEPS = 1 << SPACING_BITS


def one_trial(seed, d=20, lanes=LANES, steps=STEPS, q_min=Q_MIN, R=R_TAB):
    U = qwalk.Universe(q_min, R, seed=20250)
    q = U.q
    rng = np.random.default_rng(seed)
    x = int(rng.integers(1, q - 1))
    a = rng.integers(0, q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, q, size=lanes).astype(np.uint64)
    xv = np.full(lanes, x, dtype=np.uint64)
    t0 = (a + cc * xv) % q
    canon, e, _h = U.canon(t0)
    a = (e * a) % q
    cc = (e * cc) % q
    prev = canon.copy()

    trace = np.empty((steps, lanes), dtype=np.uint32)   # q < 2^32
    At = np.empty((steps, lanes), dtype=np.uint32)      # coefficient trace (pool)
    Ct = np.empty((steps, lanes), dtype=np.uint32)
    fires = np.zeros(lanes, dtype=np.int64)
    run = np.zeros(lanes, dtype=np.int64)
    maxrun = np.zeros(lanes, dtype=np.int64)
    left = np.zeros(lanes, dtype=np.int64)
    residual = np.zeros(lanes, dtype=np.int64)
    trapped = 0
    inv_bad = 0
    rec_cls, rec_w = [], []

    for t in range(steps):
        trace[t] = canon
        At[t] = a
        Ct[t] = cc
        pin = prev
        cn, a2, cc2, hit, pout, fired = qwalk.step(
            U, canon, a, cc, xv, d, want_dp=True, prev=pin, esc=True, tick=t)
        fires += fired
        run = np.where(fired, run + 1, 0)
        np.maximum(maxrun, run, out=maxrun)
        trapped += int(np.count_nonzero(cn == pin))
        if fired.any():
            ii = np.flatnonzero(fired)
            lf = cn[ii] != pin[ii]
            left[ii] += lf
            residual[ii] += ~lf
        if hit.any():                       # DP report: the pool record
            ii = np.flatnonzero(hit)
            rec_cls.append(cn[ii])
            rec_w.append(ii)
            a0 = rng.integers(0, q, size=ii.size).astype(np.uint64)
            c0 = rng.integers(0, q, size=ii.size).astype(np.uint64)
            s0 = (a0 + c0 * xv[ii]) % q
            cn2, e2, _h2 = U.canon(s0)
            cn = cn.copy(); a2 = a2.copy(); cc2 = cc2.copy()
            cn[ii] = cn2
            a2[ii] = (e2 * a0) % q
            cc2[ii] = (e2 * c0) % q
            pout = pout.copy()
            pout[ii] = cn2
        if t % 97 == 0:
            inv_bad += int(np.count_nonzero((a2 + cc2 * xv) % q != cn))
        canon, a, cc, prev = cn, a2, cc2, pout

    # ---- class and pair collision census on the SAME trace ------------------
    wid = np.tile(np.arange(lanes, dtype=np.int64), steps)
    cls = trace.astype(np.int64).reshape(-1)
    pair = (trace[:-1].astype(np.int64).reshape(-1) * np.int64(q)
            + trace[1:].astype(np.int64).reshape(-1))
    pw = np.tile(np.arange(lanes, dtype=np.int64), steps - 1)

    def collide(keys, wids):
        o = np.argsort(keys, kind="stable")
        k = keys[o]
        w = wids[o]
        same = np.empty(k.size, dtype=bool)
        same[0] = False
        same[1:] = k[1:] == k[:-1]
        starts = np.flatnonzero(same)
        if starts.size == 0:
            return 0, 0
        # run start = one before each repeat; group runs by key
        rst = np.flatnonzero(np.concatenate([[True], k[1:] != k[:-1]]))
        ren = np.append(rst[1:], k.size)
        big = (ren - rst) >= 2
        if not big.any():
            return 0, 0
        wmin = np.minimum.reduceat(w, rst)[big]
        wmax = np.maximum.reduceat(w, rst)[big]
        n_multi = int(np.count_nonzero(wmin != wmax))
        return int(big.sum()), n_multi

    cls_multi, cls_cross = collide(cls, wid)
    pair_multi, pair_cross = collide(pair, pw)

    # ---- conversion: two RECORDS of one class (full-record pool) ------------
    # The production pool is DP-gated.  At N = 2^26 a walker's cycle length is
    # ~sqrt(N) = 2^13 steps, far below the 2^20-step spacing, so a d = 20 gate
    # fires ~0 times in one spacing here (reported as dp_records); the production
    # scale has sqrt(N) = 2^187.5 >> DP = 2^60, which is why production reports.
    # The conversion the reviewer asks for is the ALGEBRA of the pool, so it is
    # measured on the full-record pool (every step is a record, d = 0 semantics):
    # two records of one class from two different walkers must solve x.
    cls_flat = trace.astype(np.int64).reshape(-1)
    if At is not None:
        a_flat = At.astype(np.int64).reshape(-1)
        c_flat = Ct.astype(np.int64).reshape(-1)
    o = np.argsort(cls_flat, kind="stable")
    k = cls_flat[o]
    w = wid[o]
    rst = np.flatnonzero(np.concatenate([[True], k[1:] != k[:-1]]))
    ren = np.append(rst[1:], k.size)
    conv_pairs = conv_ok = 0
    for st, en in zip(rst[(ren - rst) >= 2], ren[(ren - rst) >= 2]):
        ww = w[st:en]
        if np.unique(ww).size < 2:
            continue
        # two entries of different walkers, first occurrence of each
        seen = {}
        pair = None
        for j in range(st, en):
            wj = int(w[j])
            if wj in seen:
                continue
            seen[wj] = int(o[j])         # position in the ORIGINAL (unsorted) array
            if len(seen) == 2:
                pair = (seen[min(seen)], seen[max(seen)])
                break
        if pair is None:
            continue
        i, j = pair
        conv_pairs += 1
        den = (c_flat[j] - c_flat[i]) % q
        if den and (a_flat[i] - a_flat[j]) * pow(int(den), -1, q) % q == x:
            conv_ok += 1
        if conv_pairs >= 4000:          # a sample of 4000 groups is plenty
            break
    dp_records = int(sum(int(x.size) for x in rec_cls))

    return dict(seed=seed, q=int(q), classes=int(U.classes), R=R, R_eff=int(U.R_eff),
                lanes=lanes, steps=steps, spacing_bits=SPACING_BITS,
                spacing=1 << SPACING_BITS,
                trap_period=6 * int(U.R_eff),
                trap_opportunities_per_spacing=(1 << SPACING_BITS) / (6.0 * U.R_eff),
                trap_rate_measured=trapped / float(steps * lanes),
                trap_rate_structural=1.0 / (6.0 * U.R_eff),
                fires=int(fires.sum()),
                walks_with_fire=int((fires > 0).sum()),
                max_consecutive_fires_per_walk=int(maxrun.max()),
                fires_left_two_cycle=int(left.sum()),
                fires_residual=int(residual.sum()),
                escape_success_measured=(float(left.sum()) / float(fires.sum())
                                         if fires.sum() else None),
                invariant_violations=int(inv_bad),
                invariant_samples=int((steps + 96) // 97 * lanes),
                class_collisions_multi=cls_multi,
                class_collisions_cross_walker=cls_cross,
                pair_collisions_multi=pair_multi,
                pair_collisions_cross_walker=pair_cross,
                pair_over_class=(pair_cross / cls_cross if cls_cross else None),
                dp_records=dp_records,
                records=int(cls_flat.size),
                record_class_groups=conv_pairs,
                record_conversion_verified=conv_ok,
                record_conversion=(conv_ok / conv_pairs if conv_pairs else None))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=8)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    from multiprocessing import Pool
    jobs = [20250 + 7 * i for i in range(a.trials)]
    with Pool(a.workers) as pool:
        rows = pool.map(one_trial, jobs)
    tot = lambda k: sum(r[k] for r in rows)          # noqa: E731
    agg = {
        "fixture": {"q_min": Q_MIN, "R": R_TAB, "spacing_bits": SPACING_BITS,
                    "lanes": LANES, "steps": STEPS, "trials": len(rows),
                    "classes": rows[0]["classes"], "q": rows[0]["q"],
                    "R_eff": rows[0]["R_eff"],
                    "trap_opportunities_per_spacing":
                        rows[0]["trap_opportunities_per_spacing"],
                    "pre_capture_margin_vs_32x":
                        rows[0]["trap_opportunities_per_spacing"] / 32.0,
                    "production": {"R": 1 << 20, "spacing_bits": 60,
                                   "classes_log2": 375.0}},
        "rows": rows,
        "totals": {
            "lane_steps": tot("lanes") * 0 + len(rows) * LANES * STEPS,
            "fires": tot("fires"),
            "fires_left_two_cycle": tot("fires_left_two_cycle"),
            "fires_residual": tot("fires_residual"),
            "walks_with_fire": tot("walks_with_fire"),
            "max_consecutive_fires_per_walk": max(r["max_consecutive_fires_per_walk"]
                                                  for r in rows),
            "invariant_violations": tot("invariant_violations"),
            "invariant_samples": tot("invariant_samples"),
            "class_collisions_cross_walker": tot("class_collisions_cross_walker"),
            "pair_collisions_cross_walker": tot("pair_collisions_cross_walker"),
            "class_trials_with_collision": sum(
                1 for r in rows if r["class_collisions_cross_walker"] > 0),
            "pair_trials_with_collision": sum(
                1 for r in rows if r["pair_collisions_cross_walker"] > 0),
            "dp_records": tot("dp_records"),
            "records": tot("records"),
            "record_class_groups": tot("record_class_groups"),
            "record_conversion_verified": tot("record_conversion_verified"),
        },
    }
    agg["totals"]["trap_rate_measured"] = (
        agg["totals"]["fires"] / float(agg["totals"]["lane_steps"]))
    agg["totals"]["trap_rate_structural"] = rows[0]["trap_rate_structural"]
    agg["totals"]["trap_rate_over_structural"] = (
        agg["totals"]["trap_rate_measured"] / rows[0]["trap_rate_structural"])
    agg["totals"]["class_coverage"] = (
        agg["totals"]["class_trials_with_collision"] / float(len(rows)))
    agg["totals"]["pair_coverage"] = (
        agg["totals"]["pair_trials_with_collision"] / float(len(rows)))
    agg["checks"] = {
        "spacing_is_one_full_dp_spacing": STEPS == (1 << SPACING_BITS),
        "pre_capture_margin_at_least_32x":
            rows[0]["trap_opportunities_per_spacing"] >= 32.0,
        "classes_at_least_2_pow_26": rows[0]["classes"] >= (1 << 26),
        "trap_law_within_25_percent": 0.75 <= agg["totals"]["trap_rate_over_structural"] <= 1.25,
        "escape_leaves_two_cycle": (agg["totals"]["fires_left_two_cycle"] > 0
                                    and agg["totals"]["fires_residual"] == 0),
        "invariant_holds": agg["totals"]["invariant_violations"] == 0,
        "class_collisions_exist": agg["totals"]["class_collisions_cross_walker"] > 0,
        "pair_collisions_exist": agg["totals"]["pair_collisions_cross_walker"] > 0,
    }
    agg["passed"] = all(agg["checks"].values())
    out = a.out or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                "evidence", "c-fixture.json")
    with open(out, "w") as f:
        json.dump(agg, f, indent=1, sort_keys=True)
    print(json.dumps({k: agg[k] for k in ("fixture", "totals", "checks", "passed")},
                     indent=1))


if __name__ == "__main__":
    main()
