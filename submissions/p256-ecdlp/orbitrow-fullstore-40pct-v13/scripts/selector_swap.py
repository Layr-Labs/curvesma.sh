#!/usr/bin/env python3
"""selector_swap.py -- RIDER 2: the two selector arms at MATCHED range/config.

v12's swap was range-confounded (verdict 2620b615 point 3): the class-value arm
computed (canon*R)>>40 while canon <= 2.1e9, so it reached 32 of 16384 columns
and its 501x trap ratio was a 512x range artefact.  Here BOTH arms

  * take the SAME start class, drawn uniformly on [0, q/2) -- which at lam = 1
    is exactly the production class law, K = min(y, p-y) uniform on [0, p/2];
  * use the SAME R, the SAME u-table, the SAME escape, the SAME step budget;
  * reach ALL R columns (asserted per arm by direct count, not by formula);
  * differ ONLY in how the class value is turned into a row index:
      arm EC   : m = (canon * R) >> (bitlen(q) - 1)      (real class value)
      arm HASH : m = cell_span(sm64(canon) >> 24, R, 2^40)  (the shipped surrogate)

The measured rows are the column-reach count, the two-step-return (trap) rate
with a 95% interval per arm, and the difference with a 95% interval.  A
degenerate arm (fewer than R/2 columns) is refused rather than reported.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402
import selector  # noqa: E402

R = int(os.environ.get("RCOLS", "256"))
SETS = int(os.environ.get("SETS", "800"))
W = int(os.environ.get("W", "16"))
STEPS = int(os.environ.get("STEPS", "400"))
SEED = int(os.environ.get("SEED", "20250"))


class Arm(qwalk.Universe):
    """Universe whose row map is one of the two selector arms."""

    def __init__(self, q_min, R_, seed, arm):
        super().__init__(q_min, R_, seed=seed)
        self.arm = arm
        self.reached = set()

    def col(self, canon):
        canon = np.asarray(canon, dtype=np.uint64)
        if self.arm == "EC":
            # real class value on [0, q/2): shift bitlen(q)-1 spans exactly R
            m = (canon.astype(object) * self.R >> (self.q.bit_length() - 1))
            m = np.array([int(v) for v in m], dtype=np.int64)
        elif self.arm == "HASH":
            K2 = qwalk.sm64(canon ^ np.uint64(0x243F6A8885A308D3)) >> np.uint64(24)
            m = selector.cell_span(K2.astype(np.int64), self.R, 1 << 40).astype(np.int64)
        else:
            raise ValueError(self.arm)
        assert int(m.min(initial=0)) >= 0 and int(m.max(initial=0)) < self.R, (self.arm, m.max())
        self.reached.update(m.tolist())
        return m


def cp95_wilson(k, n):
    if n == 0:
        return (0.0, 0.0)
    z = 1.959963984540054
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def run_arm(qmin, arm, starts, xv, a0, cc0):
    U = Arm(qmin, R, SEED, arm)
    lanes = starts.size
    canon = starts.copy(); a = a0.copy(); cc = cc0.copy()
    hist0 = canon.copy()
    hist1 = canon.copy()
    traps = 0
    tot = 0
    for t in range(STEPS):
        cs, as_, cs_ = canon, a, cc
        cn, a2, cc2, _h, _pv, _f = qwalk.step(
            U, cs, as_, cs_, xv, 0, False, hist0, False, 0, t, None, None, None, None)
        traps += int(np.count_nonzero(cn == hist1))     # class two steps back
        tot += lanes
        hist1 = hist0
        hist0 = cs
        a, cc, canon = a2, cc2, cn
    return {"arm": arm, "traps": traps, "lane_steps": tot, "reach": len(U.reached),
            "rate": traps / tot}


def main():
    qmin = int(os.environ.get("QMIN", "4290000000"))
    U0 = qwalk.Universe(qmin, R, seed=SEED)
    q = U0.q
    rng = np.random.default_rng(SEED)
    # start classes uniform on [0, q/2): the production class law at lam = 1
    starts = rng.integers(0, q // 2, size=SETS * W).astype(np.uint64)
    xt = rng.integers(1, q - 1, size=SETS).astype(np.uint64)
    xv = np.repeat(xt, W).astype(np.uint64)
    a0 = rng.integers(0, q, size=SETS * W).astype(np.uint64)
    cc0 = rng.integers(0, q, size=SETS * W).astype(np.uint64)
    # rescale the coefficient pair so the invariant a + c*x == class holds
    a0 = (a0 - cc0 * xv) % q
    out = {"test": "selector_swap_matched", "q": int(q), "R": R, "sets": SETS,
           "walkers": W, "steps": STEPS, "seed": SEED, "arms": {}}
    for arm in ("EC", "HASH"):
        rec = run_arm(qmin, arm, starts, xv, a0, cc0)
        assert rec["reach"] >= R // 2, ("degenerate arm", arm, rec["reach"])
        lo, hi = cp95_wilson(rec["traps"], rec["lane_steps"])
        rec.update({"rate_ci95": [lo, hi], "rate_per_step": rec["rate"]})
        out["arms"][arm] = rec
    e, h = out["arms"]["EC"], out["arms"]["HASH"]
    dp = h["rate"] - e["rate"]
    se = math.sqrt(e["rate"] * (1 - e["rate"]) / e["lane_steps"]
                   + h["rate"] * (1 - h["rate"]) / h["lane_steps"])
    out["delta_hash_minus_ec"] = dp
    out["delta_ci95"] = [dp - 1.959963984540054 * se, dp + 1.959963984540054 * se]
    out["range_matched"] = {"EC_reach": e["reach"], "HASH_reach": h["reach"],
                            "R": R, "equal_reach": e["reach"] == h["reach"] == R}
    out["law_1_over_2R"] = 1 / (2 * R)
    out["both_arms_in_ci_of_law"] = all(
        a["rate_ci95"][0] <= out["law_1_over_2R"] <= a["rate_ci95"][1]
        for a in (e, h))
    out["verdict"] = ("selector surrogate is FAITHFUL (delta within CI)" if
                      out["delta_ci95"][0] <= 0 <= out["delta_ci95"][1] else
                      "selector arms differ beyond CI -- surrogate NOT faithful")
    with open(os.path.join(D, "evidence", "selector-swap.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps({k: out[k] for k in ("delta_hash_minus_ec", "delta_ci95",
                                          "range_matched", "law_1_over_2R",
                                          "both_arms_in_ci_of_law", "verdict")},
                     indent=1))
    for a in (e, h):
        print(a["arm"], "reach", a["reach"], "traps", a["traps"], "/", a["lane_steps"],
              "rate %.6f" % a["rate"], "CI", ["%.6f" % v for v in a["rate_ci95"]])


if __name__ == "__main__":
    main()
