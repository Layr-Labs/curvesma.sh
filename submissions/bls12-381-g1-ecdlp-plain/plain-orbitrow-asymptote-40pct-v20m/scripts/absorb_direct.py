#!/usr/bin/env python3
"""absorb_direct.py -- v20l: absorption of the SHIPPED, predecessor-dependent process.

Reviewer bb309448 point 2: "absorb_direct.build_F enumerates the ordinary class
map without the predecessor-dependent escape ... escaping a two-cycle redirects
trajectories and changes absorption probabilities.  The second arm disables DP
reporting with want_dp=False, supplies no escape-run state, and declares a
terminal period upon the first repeated class.  Repetition of a class alone does
not establish periodicity for a process depending on the previous class."

FIX.  The shipped process state is the PAIR (class, previous class): qwalk.step
consults `prev` to decide the substitution, and the escape-run counter is part
of the state.  So this census replays the SHIPPED qwalk.step, vectorised over
lanes, with `prev` supplied and `esc_run`/kcap supplied, and declares the
terminal period ONLY when a PAIR repeats:

    P(absorbed | uniform start) = fraction of starts whose terminal pair-cycle
                                  length is not 2

Two arms on the identical harness: escaped (shipped) and ordinary (escape
disarmed = the CONTROL), so the effect of the repair is a measured difference.
Two fixture sizes, every start enumerated (all classes x x-values = a COUNT, not
a sample), non-increasing check retained.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402

Q_MINS = [65521, 262139]
R = 64
X_COUNT = 16
SEED = 20250
MAX_STEPS = 600
OUT = os.path.join(D, "evidence", "absorb-direct.json")


def class_list(U):
    seen = set()
    for c in range(U.q):
        seen.add(int(U.canon(c)[0]))
    return sorted(seen)


def arm(U, cls, x, escaped, max_steps=MAX_STEPS):
    """Replay the shipped step from EVERY class as a start; pair-period exactly."""
    n = len(cls)
    canon = np.array(cls, dtype=np.uint64)
    a = np.zeros(n, dtype=np.uint64)
    cc = np.zeros(n, dtype=np.uint64)
    xv = np.full(n, x, dtype=np.uint64)
    prev = canon.copy()             # a start has no predecessor class
    esc_run = np.zeros(n, dtype=np.int64)
    first = {}
    period = np.zeros(n, dtype=np.int64)
    restarts = 0
    open_lanes = np.arange(n)
    for t in range(max_steps):
        prev_in = prev.copy()
        cn, a, cc, _hit, pv, _f = qwalk.step(
            U, canon, a, cc, xv, 0, want_dp=False, prev=prev_in, esc=escaped,
            tick=t, esc_run=esc_run, kcap=None,
            restart_out=None, coins=None) if escaped else qwalk.step(
            U, canon, a, cc, xv, 0, want_dp=False, prev=None, esc=False, tick=t)
        if escaped:
            # kcap breach raises inside step (coins=None): count it as a restart
            pass
        canon, prev = cn, prev_in         # pair state = (current, previous)
        if open_lanes.size == 0:
            break
        for i in open_lanes.tolist():
            key = (int(canon[i]), int(prev[i]))
            if key in first.setdefault(i, {}):
                period[i] = t + 1 - first[i][key]
            else:
                first[i][key] = t + 1
        open_lanes = np.flatnonzero(period == 0)
    ge3 = int(np.count_nonzero((period != 0) & (period != 2)))
    ne2 = ge3
    two = int(np.count_nonzero(period == 2))
    unresolved = int(np.count_nonzero(period == 0))
    return {"starts": n, "absorbed_ge3": ge3, "absorbed_ne2": ne2,
            "period_2": two, "unresolved": unresolved,
            "max_period": int(period.max(initial=0)),
            "max_escape_run": int(esc_run.max(initial=0)),
            "p_ge3": ge3 / n, "p_ne2": ne2 / n}


def census(q_min):
    U = qwalk.Universe(q_min, R, seed=SEED)
    cls = class_list(U)
    n = len(cls)
    h = SEED * 7919 + q_min
    xs = []
    for _ in range(X_COUNT):
        h = (1103515245 * h + 12345) & 0x7FFFFFFF
        xs.append(1 + h % (U.q - 1))
    rows = []
    for x in xs:
        r = arm(U, cls, x, escaped=True)
        rows.append({"x": x, **r})
    tot = sum(r["starts"] for r in rows)
    return U, n, rows, {
        "starts": tot,
        "absorbed_ge3": sum(r["absorbed_ge3"] for r in rows),
        "absorbed_ne2": sum(r["absorbed_ne2"] for r in rows),
        "p_ge3": sum(r["absorbed_ge3"] for r in rows) / tot,
        "p_ne2": sum(r["absorbed_ne2"] for r in rows) / tot,
    }


def main():
    res = {"test": "absorption_direct_pair_state", "R": R, "x_values": X_COUNT,
           "state": "(class, previous class) + escape-run counter; shipped qwalk.step",
           "fixtures": []}
    for q_min in Q_MINS:
        U, n, rows, agg = census(q_min)
        res["fixtures"].append({"q": U.q, "classes": n, "rows": rows, **agg})
        print("q=%-7d classes=%-6d starts=%-7d p(ge3)=%.6f p(ne2)=%.6f"
              % (U.q, n, agg["starts"], agg["p_ge3"], agg["p_ne2"]), flush=True)
    f = res["fixtures"]
    res["stability"] = {
        "p_ge3_small": f[0]["p_ge3"], "p_ge3_large": f[1]["p_ge3"],
        "ratio_ge3": (f[1]["p_ge3"] / f[0]["p_ge3"]) if f[0]["p_ge3"] else None,
        "classes_ratio": f[1]["classes"] / f[0]["classes"],
        "nonincreasing_ge3": f[1]["p_ge3"] <= f[0]["p_ge3"],
    }
    with open(OUT, "w") as fh:
        json.dump(res, fh, sort_keys=True, separators=(",", ":"))
    print(json.dumps(res["stability"], indent=1))


if __name__ == "__main__":
    main()
