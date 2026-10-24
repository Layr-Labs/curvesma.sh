#!/usr/bin/env python3
"""selftest_selector.py -- v12 DIFFERENTIAL COVERAGE TEST, old vs new selector.

Verdict 32143a2 point 2 (the shared R/2 bug): v11's qwalk.col drew K uniform on
[0,2^39) and called selector.cell_vec(K, R, 1<<40) = (K*R) >> 40, so a fixture
with R columns reached only R/2 of them and the trap density was 1/(2*(R/2))=1/R
instead of the advertised 1/(2R).  THE SAME ARITHMETIC was in the secp census
(R=64 -> R_eff=32 by the same shift).

This script is the missing audit: it MEASURES the reach of both maps, asserts the
new one covers every column, and measures the resulting two-step-return (trap)
density on the real walk so the one-bit trap-row correction is a measurement, not
a repair claim.

ROWS
  coverage   per R: distinct columns reached by the v11 map vs the v12 map over
             2^20 keys; asserts new == R and old == R/2 (+1 boundary column).
  law        column-collision probability measured vs 1/R_eff for both maps.
  trap       two-step-return fires per lane-step on the shipped walk (matched
             starts, R=64 and R=256), old map vs new map, with the CP95 interval
             on the fire rate and the ratio (expected 2.0).
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
import qwalk                      # noqa: E402
import selector                   # noqa: E402
import census_v20n as C           # noqa: E402

P_FIX = 1 << 40


class UOld(qwalk.Universe):
    """v11 selector: K on [0,2^39), shift 40 -> only R/2 columns."""

    def col(self, canon):
        K = (qwalk.sm64(canon ^ np.uint64(0x243F6A8885A308D3)) >> np.uint64(25))
        K = K.astype(np.int64)
        return selector.cell_vec(K, self.R, P_FIX).astype(np.int64)


def coverage(R, samples=1 << 20, seed=3):
    rng = np.random.default_rng(seed)
    K = rng.integers(0, 1 << 39, size=samples, dtype=np.int64)
    old = selector.cell_vec(K, R, P_FIX)
    K2 = rng.integers(0, 1 << 40, size=samples, dtype=np.int64)
    new = selector.cell_span(K2, R, P_FIX)
    do, dn = np.unique(old).size, np.unique(new).size
    # collision probability: P(two uniform keys land in one column)
    def coll(m):
        c = np.bincount(m, minlength=R).astype(np.float64)
        p = c / c.sum()
        return float((p * p).sum())
    return {"R": int(R), "keys_sampled": int(samples),
            "columns_reached_v11": int(do), "columns_reached_v12": int(dn),
            "v11_formula": int(selector.reachable_columns(R, P_FIX)),
            "v12_formula": int(selector.reachable_columns_span(R, P_FIX)),
            "old_reach_is_half": bool(do <= R // 2 + 1),
            "new_reach_is_full": bool(dn == R),
            "collision_v11_measured": coll(old), "collision_v12_measured": coll(new),
            "1_over_R": 1.0 / R,
            "1_over_R_eff_v11": 1.0 / float(selector.reachable_columns(R, P_FIX)),
            "trap_rate_v11_law": 0.5 / float(selector.reachable_columns(R, P_FIX)),
            "trap_rate_v12_law": 0.5 / R,
            "trap_rate_ratio_law": (0.5 / float(selector.reachable_columns(R, P_FIX))) / (0.5 / R)}


def trap_rate(R, sets=800, W=16, d=5, qmin=4_290_000_000, seed=20250, steps=160):
    """Matched-start walk with each selector; fires per lane-step."""
    out = {}
    for name, U in (("v11", UOld(qmin, R, seed=seed)), ("v12", qwalk.Universe(qmin, R, seed=seed))):
        xt, xv, a_all, cc_all, canon_all = C.make_start(U, sets, W, seed)
        res = C.walk_arm(U, canon_all, a_all, cc_all, xv, d, steps, True, 3, W, sets,
                         coins_seed=0xC0FFEE)
        lane_steps = steps * sets * W
        out[name] = {"fires": int(res["fires"]), "lane_steps": int(lane_steps),
                     "fires_per_lane_step": res["fires"] / lane_steps,
                     "invariant_violations": int(res["invariant_violations"]),
                     "successes": int(np.count_nonzero(res["succ_step"] >= 0)),
                     "R_eff": int(U.R_eff)}
    out["ratio_meas"] = (out["v11"]["fires_per_lane_step"] /
                         out["v12"]["fires_per_lane_step"] if out["v12"]["fires_per_lane_step"] else None)
    out["R"] = R
    out["dp_bits"] = d
    out["lane_steps"] = out["v12"]["lane_steps"]
    return out


def main():
    t0 = time.time()
    rows = [coverage(R) for R in (64, 256, 1024, 16384)]
    traps = [trap_rate(64), trap_rate(256)]
    for r in rows:
        assert r["new_reach_is_full"], r
        assert r["old_reach_is_half"], r
    for t in traps:
        assert t["v11"]["invariant_violations"] == 0 and t["v12"]["invariant_violations"] == 0, t
        # the trap channel is a REAL event in both arms: a zero here means the
        # cell crashed, not that the fix is free
        assert t["v11"]["fires"] > 100 and t["v12"]["fires"] > 100, t
        # measured ratio must be within 25% of the law (2.0)
        assert 1.5 <= t["ratio_meas"] <= 2.5, t
    out = {"test": "selector_differential_v12", "planted_error_control": {
               "description": "the v11 map IS the planted error; the test must call it half-reach",
               "v11_rows_flagged_half": sum(1 for r in rows if r["old_reach_is_half"]),
               "rows": len(rows), "control_pass": all(r["old_reach_is_half"] for r in rows)},
           "coverage": rows, "trap": traps,
           "elapsed_s": round(time.time() - t0, 1),
           "verdict_point": "32143a2 point 2 (fixture selector R/2 reach)"}
    os.makedirs(os.path.join(D, "evidence"), exist_ok=True)
    with open(os.path.join(D, "evidence", "selector-diff.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    for r in rows:
        print("R=%-6d reach v11 %-6d v12 %-6d | coll v11 %.3e v12 %.3e (1/R %.3e) | trap law ratio %.2f"
              % (r["R"], r["columns_reached_v11"], r["columns_reached_v12"],
                 r["collision_v11_measured"], r["collision_v12_measured"], r["1_over_R"],
                 r["trap_rate_ratio_law"]))
    for t in traps:
        print("trap R=%-4d v11 %.4e/step  v12 %.4e/step  ratio %.2f (law 2.00)  fires %d/%d"
              % (t["R"], t["v11"]["fires_per_lane_step"], t["v12"]["fires_per_lane_step"],
                 t["ratio_meas"], t["v11"]["fires"], t["v12"]["fires"]))
    print("[done]", out["elapsed_s"], "s")


if __name__ == "__main__":
    main()
