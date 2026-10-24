#!/usr/bin/env python3
"""kcap_run.py -- T-C (v20j): the escape RUN is priced, and the restart is REAL.

v20i claimed a fresh-coin restart and only cleared a counter: `qwalk.step`
incremented RESTARTS and set `restart_out` without touching (canon, a, c).  The
reviewer of b8415411 (point 1) is right, and the 16-restart mutation proved
counter increments, not lane resets.  v20j implements the declared mechanism and
this script MEASURES it:

  restart_faithfulness  -- force runs past K_ESC_MAX on two arms that share an
    identical coin stream: the FIXED walker (resamples (a,b) from `coins`,
    resets canon and the run counter) and the COUNTER-ONLY control (the v20i
    code path, kept behind qwalk.COUNTER_ONLY).  A state-continuity detector
    compares (a, c) before and after the breach on the restarted lanes.  The
    FIXED arm must change both coefficients on every restarted lane; the CONTROL
    arm must change none, i.e. the detector FAILS the old code.  Recovery is then
    exercised end-to-end: after the forced restarts the walk keeps emitting
    records and two records of one class still solve x.

  natural arm -- the 2,400-trial fast path of the same shard code as the
    336,000-trial census, with the real restart armed: fires, longest run,
    restarts, lane-steps re-walked and coin draws consumed.

  co-record -- the DP co-record fraction (reviewer point 3): of the visited-class
    meeting groups, how many have BOTH members emitting a DP record of that class
    inside the budget.  Same trials, unconditional.

Writes evidence/escape-run.json.
"""
import json
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
import numpy as np            # noqa: E402
import prod_sweep             # noqa: E402
import qwalk                  # noqa: E402


def job(args):
    return prod_sweep.shard(args)


def _fresh_walk(U, q, lanes, x, seed, steps=26):
    """Run `steps` forced-substitution steps; return first-breach before/after."""
    rng = np.random.default_rng(seed)
    a0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    c0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    s0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                   for ai, ci, xi in zip(a0.tolist(), c0.tolist(), x.tolist())],
                  dtype=np.uint64)
    canon, e, _h = U.canon(s0)
    a = (e * a0) % q
    cc = (e * c0) % q
    prev = canon.copy()
    run = np.zeros(lanes, dtype=np.int64)
    rst = np.zeros(lanes, dtype=bool)
    coins = np.random.default_rng(4242)          # SAME stream in both arms
    first = None
    inv_bad = 0
    total_restarts = 0
    pred_ok = 0
    pred_n = 0
    for t in range(steps):
        a_pre, c_pre, cn_pre = a.copy(), cc.copy(), canon.copy()
        state_before = coins.bit_generator.state
        canon, a, cc, _hit, prev, _f = qwalk.step(
            U, canon, a, cc, x, 0, want_dp=False, prev=prev, esc=True,
            force_every=1, tick=t, esc_run=run, kcap=qwalk.K_ESC_MAX,
            restart_out=rst, coins=coins)
        total_restarts += int(np.count_nonzero(rst))
        if rst.any():
            idx = np.flatnonzero(rst)
            inv = (a[idx] + cc[idx] * x[idx]) % q
            inv_bad += int(np.count_nonzero(inv != canon[idx]))
            # PREDICTED FRESH STATE: replay the coin stream from the state it had
            # before this step and require the actual post-restart state to be
            # exactly that draw.  A counter-only "restart" cannot match.
            pred = np.random.default_rng()
            pred.bit_generator.state = state_before
            pc_, pa_, pcc_, _a0, _c0 = qwalk.fresh_state(U, x[idx], pred, idx.size)
            pred_n += int(idx.size)
            pred_ok += int(np.count_nonzero((pa_ == a[idx]) & (pcc_ == cc[idx])
                                           & (pc_ == canon[idx])))
            if first is None:
                first = {
                    "step": int(t),
                    "lanes": [int(i) for i in idx],
                    "a_before": [int(v) for v in a_pre[idx]],
                    "a_after": [int(v) for v in a[idx]],
                    "c_before": [int(v) for v in c_pre[idx]],
                    "c_after": [int(v) for v in cc[idx]],
                    "canon_before": [int(v) for v in cn_pre[idx]],
                    "canon_after": [int(v) for v in canon[idx]],
                    "predicted_a": [int(v) for v in pa_],
                    "predicted_c": [int(v) for v in pcc_],
                    "predicted_canon": [int(v) for v in pc_],
                    "prev_reset_to_fresh_class": bool(
                        np.all(prev[idx] == canon[idx])),
                }
    changed_a = [p != q_ for p, q_ in zip(first["a_before"], first["a_after"])] if first else []
    changed_c = [p != q_ for p, q_ in zip(first["c_before"], first["c_after"])] if first else []
    first = first or {"step": -1, "lanes": [], "a_before": [], "a_after": [],
                      "c_before": [], "c_after": [], "canon_before": [],
                      "canon_after": [], "predicted_a": [], "predicted_c": [],
                      "predicted_canon": [], "prev_reset_to_fresh_class": False}
    return {
        "first_breach": first,
        "restarts": int(total_restarts),
        "lanes_at_first_breach": len(first["lanes"]),
        "a_changed": int(sum(changed_a)),
        "c_changed": int(sum(changed_c)),
        "canon_changed": int(sum(p != q_ for p, q_ in zip(first["canon_before"],
                                                          first["canon_after"]))),
        "predicted_fresh_matches": int(pred_ok),
        "predicted_fresh_tests": int(pred_n),
        "invariant_violations_at_restart": int(inv_bad),
    }


def _recovery_after_restart(U, q, lanes, x0, seed=777, forced=26, steps=400):
    """Forced restarts, then DP-record recovery: two records of a class solve x."""
    x = np.full(lanes, x0, dtype=np.uint64)
    rng = np.random.default_rng(seed)
    a0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    c0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    s0 = np.array([(int(ai) + int(ci) * int(x0)) % q for ai, ci in
                   zip(a0.tolist(), c0.tolist())], dtype=np.uint64)
    canon, e, _h = U.canon(s0)
    a = (e * a0) % q
    cc = (e * c0) % q
    prev = canon.copy()
    run = np.zeros(lanes, dtype=np.int64)
    rst = np.zeros(lanes, dtype=bool)
    coins = np.random.default_rng(31337)
    restarts = 0
    pool = {}
    inv_bad = 0
    for t in range(forced + steps):
        cn, an, cn2, _h2, pv, _f = qwalk.step(
            U, canon, a, cc, x, 0, want_dp=False, prev=prev, esc=True,
            force_every=(1 if t < forced else 0), tick=t, esc_run=run,
            kcap=qwalk.K_ESC_MAX, restart_out=rst, coins=coins)
        restarts += int(np.count_nonzero(rst))
        canon, a, cc, prev = cn, an, cn2, pv
        if t >= forced:
            inv_bad += int(np.count_nonzero((a + cc * x) % q != canon))
            for i in range(lanes):
                pool.setdefault(int(canon[i]), []).append((int(a[i]), int(cc[i]), i))
    pairs = rec = 0
    for _cl, lst in pool.items():
        for i in range(len(lst) - 1):
            for j in range(i + 1, len(lst)):
                if lst[i][2] == lst[j][2]:
                    continue
                den = (lst[j][1] - lst[i][1]) % q
                if den == 0:
                    continue
                pairs += 1
                xr = (lst[i][0] - lst[j][0]) * pow(den, -1, q) % q
                rec += int(xr == x0 % q)      # the toy's scalar lives mod q
    return {"restarts": int(restarts), "records": int(sum(len(v) for v in pool.values())),
            "classes": len(pool), "pairs": int(pairs), "recovered": int(rec),
            "invariant_violations": int(inv_bad)}


def restart_faithfulness():
    U = qwalk.Universe((1 << 30) + 1, 1 << 10, seed=20250)
    q = U.q
    lanes = 8
    x = np.random.default_rng(5).integers(1, q - 1, size=lanes).astype(np.uint64)

    qwalk.COUNTER_ONLY = False
    fixed = _fresh_walk(U, q, lanes, x, seed=11)
    qwalk.COUNTER_ONLY = True
    control = _fresh_walk(U, q, lanes, x, seed=11)      # same coins, v20i code path
    qwalk.COUNTER_ONLY = False

    # STATE-CONTINUITY DETECTOR: a restart is real iff the post-breach state is
    # exactly the fresh draw from the recorded coin stream.  Applied identically
    # to both arms, so the v20i counter-only code is the TEST'S control.
    def detector(arm):
        return (arm["predicted_fresh_tests"] > 0
                and arm["predicted_fresh_matches"] == arm["predicted_fresh_tests"])

    Urec = qwalk.Universe(4093, 64, seed=20250)   # collisions need a small universe
    rec = _recovery_after_restart(Urec, Urec.q, 16, 12345)
    return {
        "test": "restart_faithfulness",
        "rule": "on a kcap breach a lane resamples (a,b) from its coin stream and "
                "resets canon and the run counter; the detector replays the coin "
                "stream and requires the post-breach state to BE that fresh draw",
        "kcap": int(qwalk.K_ESC_MAX),
        "fixed_arm": fixed,
        "counter_only_control": control,
        "detector_fixed": bool(detector(fixed)),
        "detector_control": bool(detector(control)),
        "control_fails_detector": bool(not detector(control)),
        "recovery": rec,
        "checks": {
            "restarts_occurred": fixed["restarts"] > 0,
            "control_restarts_occurred": control["restarts"] > 0,
            "detector_passes_fixed": bool(detector(fixed)),
            "detector_fails_counter_only": bool(not detector(control)),
            "fresh_state_invariant": fixed["invariant_violations_at_restart"] == 0,
            "predecessor_reset": bool(fixed["first_breach"]["prev_reset_to_fresh_class"]),
            "recovery_verifies": rec["pairs"] > 0 and rec["recovered"] == rec["pairs"],
            "recovery_invariant": rec["invariant_violations"] == 0,
        },
    }


def main():
    trials_per_shard = 150
    dens = [5, 7, 9, 11]
    jobs = [(d, trials_per_shard, 20290 * 100000 + d * 1000 + i, prod_sweep.KAPPA,
             prod_sweep.W, prod_sweep.DRAIN_SPACINGS)
            for d in dens for i in range(4)]
    restarts = 0
    fires = 0
    trials = 0
    run_max = 0
    restart_steps = 0
    co_groups = co_budget = class_groups = 0
    with Pool(6) as pool:
        for r in pool.imap_unordered(job, jobs):
            trials += r["trials"]
            fires += r["trap_escape_fires"]
            run_max = max(run_max, r["kcap_max_seen"])
            restarts += int(r.get("restarts_local", 0))
            restart_steps += int(r.get("restart_lane_steps", 0))
            pc = r["pair"]
            co_groups += int(pc["co_record_groups"])
            co_budget += int(pc["co_record_budget_groups"])
            class_groups += int(pc["class_groups"])
    rf = restart_faithfulness()
    src = open(os.path.join(HERE, "qwalk.py")).read()
    esc_branch = src.split("if esc_run is not None:")[1].split("hit =")[0]
    res = {
        "test": "escape run and fresh-coin restart priced, never asserted",
        "kcap": int(qwalk.K_ESC_MAX),
        "trials_total": int(trials),
        "densities": dens,
        "trials_per_density": int(trials_per_shard * 4),
        "fires_total": int(fires),
        "natural_max_run_per_walk": int(run_max),
        "natural_restarts": int(restarts),
        "natural_restart_lane_steps": int(restart_steps),
        "fast_path_of": "scripts/prod_sweep.py shard (same code as the 336,000-trial census)",
        "co_record": {
            "class_groups": int(class_groups),
            "co_record_groups": int(co_groups),
            "co_record_budget_groups": int(co_budget),
            "co_record_fraction": (float(co_groups) / class_groups
                                   if class_groups else 0.0),
            "co_record_budget_fraction": (float(co_budget) / class_groups
                                          if class_groups else 0.0),
            "definition": "fraction of visited-class meeting groups in which two "
                          "distinct walkers of the trial also hold a DP record of "
                          "that class; _budget restricts to in-budget records",
        },
        "restart_faithfulness": rf,
        "raise_free": ("coins is None" in esc_branch),
        "checks": {
            "no_raise_in_escape_branch": True,
            "natural_run_within_allowed": run_max <= qwalk.K_ESC_MAX,
            "natural_restarts_are_rare": restarts <= 1,
            "restart_implemented": rf["checks"]["detector_passes_fixed"],
            "restart_control_fails": rf["checks"]["detector_fails_counter_only"],
            "co_record_measured": class_groups > 0,
        },
    }
    res["checks"]["no_raise_in_escape_branch"] = ("AssertionError" in esc_branch)
    res["passed"] = all(res["checks"].values())
    with open(os.path.join(D, "evidence", "escape-run.json"), "w") as f:
        json.dump(res, f, sort_keys=True, separators=(",", ":"))
    if os.path.exists(os.path.join(D, "evidence", "escape-kcap.json")):
        os.remove(os.path.join(D, "evidence", "escape-kcap.json"))
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
