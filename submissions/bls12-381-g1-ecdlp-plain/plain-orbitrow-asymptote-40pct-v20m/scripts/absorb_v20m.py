#!/usr/bin/env python3
"""absorb_v20m.py -- absorption of the SHIPPED predecessor-dependent process (v20m).

Reviewer 2d98d2a8 point 3, fixed exactly:
  (a) `canon, prev = cn, prev_in` pinned prev to the INITIAL class because the
      step's RETURNED predecessor `pv` was ignored.  Here prev <- pv every step
      (real sliding two-step test), and a restart resets prev as step() does.
  (b) kcap=None is replaced by the production kcap = qwalk.K_ESC = 8, with a
      coin Generator supplied so an over-long run really resamples state.
  (c) the repetition key includes the escape-run counter: key = (class, prev,
      esc_run), so a period is declared only on the FULL shipped state.

Mutation control: the same harness with prev RE-PINNED to the initial class
(the v20l defect) must change the numbers; it is run and reported.
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
X_COUNT = 8
CLASSES_PER_X = 2048
SEED = 20250
MAX_STEPS = 900
OUT = os.path.join(D, "evidence", "absorb-direct.json")


def class_list(U):
    seen = set()
    for c in range(U.q):
        seen.add(int(U.canon(c)[0]))
    return np.array(sorted(seen), dtype=np.uint64)


def arm(U, cls, xv, escaped, pin_prev=False, max_steps=MAX_STEPS, seed=SEED):
    """Replay the shipped step from every start; period on the FULL pair state."""
    n = cls.size
    canon = cls.copy()
    a = np.zeros(n, dtype=np.uint64)
    cc = np.zeros(n, dtype=np.uint64)
    prev = canon.copy()                 # a start has no predecessor class
    esc_run = np.zeros(n, dtype=np.int64)
    coins = np.random.default_rng((seed ^ 0xA5A5A5A5) & 0x7FFFFFFF)
    init_prev = canon.copy()
    keys = np.empty((max_steps, n), dtype=np.uint64)
    restarts = 0
    fired_total = 0
    for t in range(max_steps):
        # key BEFORE the step: (class, prev, esc_run) with class,prev < 2^18
        keys[t] = (canon.astype(np.uint64) << np.uint64(26)) | \
                  (prev.astype(np.uint64) << np.uint64(8)) | \
                  (esc_run.astype(np.uint64) & np.uint64(0xFF))
        rst = np.zeros(n, dtype=bool)
        cn, a, cc, _hit, pv, fired = qwalk.step(
            U, canon, a, cc, xv, 0, want_dp=False, prev=prev, esc=escaped,
            tick=t, esc_run=esc_run, kcap=(qwalk.K_ESC if escaped else None),
            restart_out=rst, coins=coins)
        restarts += int(rst.sum())
        fired_total += int(fired.sum())
        canon = cn
        prev = init_prev if pin_prev else pv      # THE FIX (and its control)
    # first repeated full-state key per lane (stable sort: order[start+1] is the
    # second-smallest original index within an equal-value group)
    period = np.zeros(n, dtype=np.int64)
    cyc = np.zeros(n, dtype=np.int64)
    for i in range(n):
        col = keys[:, i]
        order = np.argsort(col, kind="stable")
        sv = col[order]
        brk = np.flatnonzero(np.r_[True, sv[1:] != sv[:-1]])
        ends = np.r_[brk[1:], sv.size]
        g = np.flatnonzero(ends - brk >= 2)
        if g.size == 0:
            continue
        best = None
        for gi in g.tolist():
            i1 = int(order[brk[gi]])
            i2 = int(order[brk[gi] + 1])
            if best is None or i2 < best[1]:
                best = (i1, i2)
        period[i] = best[1]
        cyc[i] = best[1] - best[0]
    open_ = int(np.count_nonzero(period == 0))
    ne2 = int(np.count_nonzero((cyc != 0) & (cyc != 2)))
    ge3 = int(np.count_nonzero(cyc >= 3))
    two = int(np.count_nonzero(cyc == 2))
    return {"starts": int(n), "absorbed_ne2": ne2, "absorbed_ge3": ge3,
            "period_2": two, "open_no_repeat": open_,
            "max_cycle": int(cyc.max(initial=0)),
            "max_escape_run": int(esc_run.max(initial=0)),
            "restarts": restarts, "escape_fires": fired_total,
            "p_ne2": ne2 / n, "p_ge3": ge3 / n}


def census(q_min, pin_prev=False, tag=""):
    U = qwalk.Universe(q_min, R, seed=SEED)
    cls = class_list(U)
    nc = cls.size
    rng = np.random.default_rng((SEED + q_min) & 0x7FFFFFFF)
    xvals = (1 + rng.integers(0, U.q - 1, size=X_COUNT)).astype(np.uint64)
    cols = []
    xv = []
    for x in xvals.tolist():
        pick = rng.integers(0, nc, size=min(CLASSES_PER_X, nc))
        cols.append(cls[pick])
        xv.append(np.full(cols[-1].size, x, dtype=np.uint64))
    cls_all = np.concatenate(cols)
    xv_all = np.concatenate(xv)
    r = arm(U, cls_all, xv_all, escaped=True, pin_prev=pin_prev)
    r.update({"q": int(U.q), "classes": int(nc), "tag": tag})
    return r


def main():
    res = {"test": "absorption_direct_full_pair_state",
           "R": R, "x_values": X_COUNT, "classes_per_x": CLASSES_PER_X,
           "state": "(class, previous class, escape-run counter); shipped qwalk.step",
           "production_kcap": qwalk.K_ESC, "max_steps": MAX_STEPS,
           "fixtures": [], "mutation_control": []}
    for q_min in Q_MINS:
        r = census(q_min, tag="escape_armed")
        res["fixtures"].append(r)
        print(json.dumps(r), flush=True)
        c = census(q_min, pin_prev=True, tag="MUTATION_prev_pinned")
        res["mutation_control"].append(c)
        print(json.dumps(c), flush=True)
    f = res["fixtures"]
    m = res["mutation_control"]
    res["stability"] = {
        "p_ge3_small": f[0]["p_ge3"], "p_ge3_large": f[1]["p_ge3"],
        "ratio_ge3": (f[1]["p_ge3"] / f[0]["p_ge3"]) if f[0]["p_ge3"] else None,
        "classes_ratio": f[1]["classes"] / f[0]["classes"],
        "nonincreasing_ge3": f[1]["p_ge3"] <= f[0]["p_ge3"],
        "total_starts": sum(x["starts"] for x in f),
        "restarts_total": sum(x["restarts"] for x in f),
    }
    res["mutation"] = {
        "p_ge3_pinned_small": m[0]["p_ge3"], "p_ge3_pinned_large": m[1]["p_ge3"],
        "changed_small": m[0]["p_ge3"] != f[0]["p_ge3"],
        "changed_large": m[1]["p_ge3"] != f[1]["p_ge3"],
        "abs_delta_small": abs(m[0]["p_ge3"] - f[0]["p_ge3"]),
        "abs_delta_large": abs(m[1]["p_ge3"] - f[1]["p_ge3"]),
    }
    with open(OUT, "w") as fh:
        json.dump(res, fh, sort_keys=True, separators=(",", ":"))
    print(json.dumps(res["stability"], indent=1))
    print(json.dumps(res["mutation"], indent=1))


if __name__ == "__main__":
    main()
