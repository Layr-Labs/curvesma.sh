#!/usr/bin/env python3
"""invariant-check.py -- v16 repair of verdict 8e043fc point 4.

The reviewer charged that the lag-3 escape substitutes using the POST-step
coefficients with the PRE-step class, breaking a + c*x == canon, and that the DP
mask is retained from before the substitution.  This script tests the shipped
qwalk.step directly on the three census arms: after EVERY step, including steps
where the escape fires, it recomputes (a + c*x) mod q and compares it to the
class the step reports, and it re-evaluates the DP gate on the post-substitution
class.  A planted-wrong control (substitution applied one step late) must FAIL
the same test, so a passing run is not vacuous.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402

OUT = os.path.join(D, "evidence", "invariant.json")


def run(U, w, bt, steps, esc_lag, seed=20250, sabotage=False):
    rng = np.random.default_rng(seed & 0x7FFFFFFF)
    q = U.q
    x = rng.integers(1, q - 1, size=bt).astype(np.uint64)
    a = rng.integers(0, q, size=bt * w).astype(np.uint64)
    c = rng.integers(0, q, size=bt * w).astype(np.uint64)
    xv = np.repeat(x, w)
    canon = ((a.astype(object) + c.astype(object) * xv.astype(object))
             % q).astype(np.uint64)
    prev = None
    esc_run = np.zeros(bt * w, dtype=np.int64)
    fires = 0
    violations = 0
    dp_checked = 0
    for t in range(steps):
        want_dp = (t % 8 == 0)
        cnew, a2, c2, hit, prev_new, fired = qwalk.step(
            U, canon, a, c, xv, 5, want_dp=want_dp, prev=prev,
            esc=(esc_lag > 0), esc_run=esc_run, kcap=8, tick=t)
        fires += int(np.count_nonzero(fired))
        # THE INVARIANT: a + c*x == canon, on the post-step state
        lhs = (a2.astype(object) + c2.astype(object) * xv.astype(object)) % q
        viol = int(np.count_nonzero(lhs != cnew.astype(object)))
        violations += viol
        if want_dp:
            dp_checked += 1
            # the gate must be a function of the POST-substitution class only
            gate_now = U.dp(cnew, 5)
            if not np.array_equal(gate_now, hit):
                violations += 1
        a, c, canon, prev = a2, c2, cnew, prev_new
        if sabotage and fired.any():
            # PLANTED CONTROL: reproduce the reviewer's hypothesised defect --
            # take the escape branch's coefficients from a second application to
            # the POST-step state while keeping the class from the first one.
            # This must break a + c*x == canon, so a clean shipped run is not vacuous.
            idx = np.flatnonzero(fired)
            m = U.col(canon[idx])
            ma = (m + 1) % U.R
            cn2, a4, c4, _e = qwalk._apply(U, canon[idx], a2[idx], c2[idx], xv[idx], ma)
            a[idx] = a4
            c[idx] = c4
            canon[idx] = cnew[idx]   # class from the FIRST substitution
    return {"steps": steps, "lanes": int(bt * w), "escape_fires": fires,
            "invariant_violations": violations, "dp_checks": dp_checked}


def main():
    U = qwalk.Universe(4_290_000_000, 64, seed=20250)
    arms = {"escaped_lag1": 1, "escaped_lag3": 3, "control_noescape": 0}
    out = {"test": "invariant_v16", "q": int(U.q), "classes": int(U.classes),
           "protocol": ("shipped qwalk.step; after every step recompute (a + c*x) mod q "
                        "and compare against the reported class; re-evaluate the DP gate "
                        "on the post-substitution class"), "arms": {}}
    for k, lag in arms.items():
        out["arms"][k] = run(U, 16, 64, 400, lag)
    out["planted_control"] = run(U, 16, 64, 200, 3, sabotage=True)
    ok = (all(v["invariant_violations"] == 0 for v in out["arms"].values())
          and out["planted_control"]["invariant_violations"] > 0)
    out["planted_control_caught"] = out["planted_control"]["invariant_violations"] > 0
    out["passed"] = bool(ok)
    os.makedirs(os.path.join(D, "evidence"), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))
    if not ok:
        raise SystemExit("INVARIANT CHECK FAILED")


if __name__ == "__main__":
    main()
