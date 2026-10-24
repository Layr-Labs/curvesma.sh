#!/usr/bin/env python3
"""sensitivity_arm.py -- v19: the PARAMETER-SENSITIVITY ladder (verdict 46b1b85b).

The v18 verdict accepts the four-order ladder as relevant evidence and rejects the
TRANSFER, naming exactly why: the surrogate differs from production in R (64 vs
2^20), walkers (16 vs 2^20), DP density (2^-5 vs 2^-60) and the selector type, and
"short-cycle absorption, history-dependent escape and distinguished-point detection
depend on these other parameters".  Varying ORDER alone cannot see them.

This script varies ONE parameter at a time at FIXED ORDER (q_min = 2^28,
N_classes ~ 2^25.42 -- the cheap arm whose v18 baseline is p_set = 0.5400) with
800 sets per cell, holding the dimensionless shape x = visits/sqrt(N) fixed, so
TOTAL WORK PER SET is constant across every cell and the only thing that changes
is how that work is laid out / which parameter the map is built with.

CELLS (each 800 sets, same-challenge two-run pairing, same CP helper):
    WALKERS  W = 16, 64, 256          (per-walker steps scale as 1/W: work fixed)
    COLUMNS  R = 64, 256, 1024        (R_eff = reachable_columns(R, p))
    DP GAP   d = 5, 8, 11             (report gate 2^-d per visit)
    CAP      K_ESC = 8 (shipped) vs 10**9 (no cap)   -- direct absorption-loss row

Every cell reports p_set, its CP-95 lower end, per-walker steps, the DP-record count
per walker (per/2^d) and its own restart count on its OWN denominator (800 sets x
2 runs = 1600 armed run-trials per cell; NOT pooled across cells -- verdict
46b1b85b's restart objection).

Sensitivity is the slope of ln p_set against ln(parameter), least-squares on the
three log-spaced points, with the standard error of that slope.  The transfer
budget multiplies |slope| by the log-distance from the surrogate value to the
production value, per parameter, and is reported as a table (ledger.py reads
evidence/sensitivity.json; nothing is typed into the note).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import transfer_arm as T           # noqa: E402  (run_arm: same statistic)
import qwalk                       # noqa: E402
import census_v20n as C            # noqa: E402

OUTDIR = os.path.join(D, "evidence")
QMIN = 1 << 28                     # FIXED ORDER for every cell (v18 baseline arm L28)
SETS = 9600
SEED = 20250

CELLS = {
    # name          R      d    W     kcap
    "S_w16":  dict(R=64,   d=5,  W=16,  kcap=8),      # v18 L28 baseline, 800 sets
    "S_w64":  dict(R=64,   d=5,  W=64,  kcap=8),
    "S_w256": dict(R=64,   d=5,  W=256, kcap=8),
    "S_R256": dict(R=256,  d=5,  W=16,  kcap=8),
    "S_R1024": dict(R=1024, d=5, W=16,  kcap=8),
    "S_d8":   dict(R=64,   d=8,  W=16,  kcap=8),
    "S_d11":  dict(R=64,   d=11, W=16,  kcap=8),
    "S_nocap": dict(R=64,  d=5,  W=16,  kcap=10 ** 9),  # cap OFF: absorption-loss row
    # d = 0 (report gate OFF): every visit is a record, so p_set measures PURE
    # COALESCENCE.  The d=5 walker row above is confounded with detection (at
    # W=256 the per-walker budget is 26 steps < 2^5), so the walker and column
    # slopes are taken on THESE cells and the DP row is the measured ratio.
    # d = 0 is NOT usable: walk_arm sets `hit = U.dp(cn, d) if d > 0 else None`, so
    # a gate-off arm produces no records at all and reads p_set = 0.  Measured and
    # discarded (the five S0_* cells in this file's history), not shipped as data.
    #
    # PLATEAU cells: f = per/2^d is the detection variable.  d = 4 and d = 3 raise f
    # above the d=5 point by 2x and 4x, so the SLOPE OF THE SATURATED PLATEAU is
    # measured rather than assumed -- production sits 2^33.5 further out.
    "S_d4":  dict(R=64, d=4, W=16, kcap=8),
    "S_d3":  dict(R=64, d=3, W=16, kcap=8),
    # INVARIANCE cells: f = per/2^d = 0.83 at three DIFFERENT (W, d) pairs.  If
    # p_set collapses onto f, the walker and DP dimensions are one variable, and
    # neither is a free extrapolation.
    "S_w64d7":  dict(R=64, d=7,  W=64,  kcap=8),   # per 106 / 2^7  = 0.83
    "S_w128d6": dict(R=64, d=6,  W=128, kcap=8),   # per  53 / 2^6  = 0.83
}


def run_cell(name, R, d, W, kcap):
    saved = qwalk.K_ESC
    qwalk.K_ESC = int(kcap)
    try:
        out = T.run_arm(name, QMIN, SETS, R=R, d=d, W=W, seed=SEED)
    finally:
        qwalk.K_ESC = saved
    out["kcap"] = int(kcap)
    out["armed_run_trials"] = 2 * SETS          # this cell's OWN denominator
    out["cp_denominator_note"] = ("800 sets x 2 runs; NOT pooled with any other cell")
    out["dp_records_per_walker"] = (out["per_walker_steps"] / (1 << d)) if d < 64 else 0.0
    with open(os.path.join(OUTDIR, "transfer-arm-%s.json" % name), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    return out


def main():
    names = sys.argv[1:] or list(CELLS)
    os.makedirs(OUTDIR, exist_ok=True)
    for n in names:
        run_cell(n, **CELLS[n])


if __name__ == "__main__":
    main()
