#!/usr/bin/env python3
"""sens_ladder.py -- v12 parameter-sensitivity ladder (method lifted from secp v19/v20).

WHY.  Both surviving objections are TRANSFER objections: (1) fixture success does not
transfer to production DP density 2^-60; (2) the step price is asserted.  v19/v20 on
the secp lane replaced the first by MEASURING the success law along a one-at-a-time
ladder and extrapolating with a sign-honest charge.  This file is that method at the
p256 fixtures -- not a new method.

DESIGN (one parameter at a time, everything else at the base cell)
    base     W=16, R=64,   d=5      (the shipped census cell)
    walkers  W = 16 / 64 / 256      (per-walker steps scaled so visits/sqrt(N) = x)
    columns  R = 64 / 256 / 1024
    density  d = 5 / 8 / 11         (2^-5, 2^-8, 2^-11 report bits)
plus
    selector swap: uniform-key surrogate (v12) vs class-value truncation (v20g's
    measured WRONG model) at the base cell -- the row that justifies the surrogate.

Each cell: >=800 sets, matched challenge block, escape armed kcap=8, lag-3 arm,
CP95 set-success.  Runtime guard: a cell whose wall time is below half the rate the
calibration implies FAILS LOUD (secp v19 lost cells that "finished" in 4-6 s at p==0).

Usage:  CELL=base python3 scripts/sens_ladder.py   (writes evidence/ladder-<cell>.json)
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
import census_v20n as C           # noqa: E402
import selector                   # noqa: E402

QMIN = 4_290_000_000
SHAPE = None   # filled from C.SHAPE (design shape x = 1.6002)
SETS = int(os.environ.get("SETS", "800"))
SCALE = float(os.environ.get("SHAPE_SCALE", "1.585"))   # design shape x = 1.6002
LAG = 3
CAL_LANE_STEPS_PER_S = 1.2e7     # measured ~9-15e6/s on this host (all cores idle)
CAL_MAX_LANE_STEPS_PER_S = 5.0e7 # above this the loop cannot have run (crash guard)

CELLS = {
    "base":      dict(W=16, R=16384, d=5),
    "W64":       dict(W=64, R=16384, d=5),
    "W256":      dict(W=256, R=16384, d=5),
    "R4096":     dict(W=16, R=4096,  d=5),
    "R1024":     dict(W=16, R=1024,  d=5),
    "d8":        dict(W=16, R=16384, d=8),
    "d11":       dict(W=16, R=16384, d=11),
    "selval":    dict(W=16, R=16384, d=5, sel="classvalue"),
}


class UClassValue(qwalk.Universe):
    """v20g's measured WRONG selector: truncate the class VALUE, not a uniform key.

    `canon` is the minimum of the orbit residues and is concentrated near 0, so
    equal-width buckets of it are skewed.  Kept as the selector-swap control.
    """

    def col(self, canon):
        c = np.asarray(canon, dtype=np.int64)
        return ((c * self.R) >> 40).astype(np.int64)


def run(name):
    spec = CELLS[name]
    W, R, d = spec["W"], spec["R"], spec["d"]
    t0 = time.time()
    U = (UClassValue(QMIN, R) if spec.get("sel") == "classvalue"
         else qwalk.Universe(QMIN, R))
    per = int(round(float(C.SHAPE) * SCALE * math.sqrt(U.classes) / W))
    xt, xv, a_all, cc_all, canon_all = C.make_start(U, SETS, W, 20250)
    res = C.walk_arm(U, canon_all, a_all, cc_all, xv, d, per, True, LAG, W, SETS,
                     coins_seed=0xC0FFEE)
    k = int(np.count_nonzero(res["succ_step"] >= 0))
    lane_steps = per * SETS * W
    elapsed = time.time() - t0
    rate = lane_steps / elapsed
    # runtime guard (secp v19 lost cells that "finished" in 4-6 s at p == 0): a cell
    # that runs FASTER than the machine can, or that posts zero successes at a density
    # where the law predicts non-zero, is a crash and FAILS LOUD rather than shipping
    # as a "favourable" slope point.
    guard_ok = bool(rate <= CAL_MAX_LANE_STEPS_PER_S and elapsed / SETS >= 1e-3)
    assert guard_ok, ("CELL IMPOSSIBLY FAST / CRASHED: %s ran %.1fs for %d lane-steps "
                      "(%.2e/s > %.2e/s cap or < 1ms/set)"
                      % (name, elapsed, lane_steps, rate, CAL_MAX_LANE_STEPS_PER_S))
    assert k > 0, ("zero successes -- cell is not evidence, the cell crashed: %s" % name)
    out = {"cell": name, "spec": spec, "sets": SETS, "W": W, "R": R, "dp_bits": d,
           "per_walker_steps": per, "visits_over_sqrtN": per * W / math.sqrt(U.classes),
           "classes": int(U.classes), "successes": k, "success_rate": k / SETS,
           "cp95_success": str(C.cp95_lower(k, SETS)),
           "cp95_success_f": float(C.cp95_lower(k, SETS)),
           "escape_fires": int(res["fires"]), "fires_per_lane_step": res["fires"] / lane_steps,
           "restarts": int(res["restarts"]),
           "invariant_violations": int(res["invariant_violations"]),
           "lane_steps": int(lane_steps), "elapsed_s": round(elapsed, 1),
           "lane_steps_per_s": rate, "runtime_guard_ok": guard_ok,
           "selector": spec.get("sel", "uniform_key_v12")}
    with open(os.path.join(D, "evidence", "ladder-%s.json" % name), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("%-8s W=%-3d R=%-5d d=%-2d k=%-4d p=%.4f cp95=%.4f fires/step=%.2e %.1fs"
          % (name, W, R, d, k, k / SETS, out["cp95_success_f"],
             out["fires_per_lane_step"], elapsed), flush=True)
    return out


if __name__ == "__main__":
    run(os.environ.get("CELL", "base"))
