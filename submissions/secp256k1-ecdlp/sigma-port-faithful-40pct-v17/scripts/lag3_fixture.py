#!/usr/bin/env python3
"""lag3_fixture.py -- PLANTED TEST: the lag-3 trigger fires on a real 3-cycle.

The v16 package was rejected because `esc_lag` was only ever read as a Boolean:
the "lag-1" and "lag-3" arms ran the same qwalk.step path, so nothing in the
package demonstrated that a lag-3 match was detected at all.  This script is the
missing demonstration.  It does not touch the production universe; it builds a
SMALL one (q ~ 2e4, R = 8, classes ~ 3e3) where a 3-cycle is a common event,
finds starts whose lag-1 trace really is eventually periodic with period 3, and
then runs the SHIPPED census arm (census_v20n.walk_arm) on the same starts with
lag = 1 and lag = 3.

Asserted (all read from the run, none typed):
  * planted_period3_starts > 0                  -- the fixture really contains 3-cycles
  * lag3_fires >= planted_period3_starts        -- the lag-3 arm fires on them
  * lag1_fires_lag3_channel == 0                -- the lag-1 arm cannot fire it
  * invariant_violations == 0 in BOTH arms      -- every substituted step keeps
                                                   a + c*x == class (the v16 defect)
  * for every planted start the lag-3 trajectory is NOT periodic with period 3
    over its last window              -- the substitution actually leaves the cycle
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk          # noqa: E402
import census_v20n as C  # noqa: E402

PER = 240
W = 1          # one lane per start: we are testing the trigger, not the meeting


def period_of(trace, lane=0):
    """Smallest p in 1..60 with trace[-p:] == trace[-2p:-p]; 0 if none."""
    v = trace[:, lane]
    for p in range(1, 61):
        if v.size < 2 * p:
            break
        if np.array_equal(v[-p:], v[-2 * p:-p]):
            return p
    return 0


def main():
    U = qwalk.Universe(20_000, 8, seed=99)
    assert U.classes > 1000, U.classes
    rng = np.random.default_rng(4242)
    n_starts = 96
    starts = rng.integers(1, U.q, size=n_starts).astype(np.uint64)
    a0 = rng.integers(0, U.q, size=n_starts).astype(np.uint64)
    c0 = rng.integers(0, U.q, size=n_starts).astype(np.uint64)
    t0 = (a0 + c0 * starts) % U.q
    canon0, e, _h = U.canon(t0)
    a0 = (e * a0) % U.q
    c0 = (e * c0) % U.q
    xv = starts.astype(np.uint64)

    out = {"test": "lag3_planted_fixture", "q": int(U.q), "R": U.R,
           "classes": int(U.classes), "starts": n_starts, "per": PER}
    arms = {}
    for lag, esc in ((1, True), (3, True)):
        r = C.walk_arm(U, canon0, a0, c0, xv, 0, PER, esc, lag, W, n_starts)
        arms["lag%d" % lag] = r
    # planted set: starts whose LAG-1 trace is eventually periodic with period 3
    planted = [b for b in range(n_starts)
               if period_of(arms["lag1"]["trace"], b) == 3]
    escaped = [b for b in planted
               if period_of(arms["lag3"]["trace"], b) != 3]
    per_ok = sum(1 for b in planted
                 if period_of(arms["lag1"]["trace"], b) == 3)
    out["planted_period3_starts"] = len(planted)
    out["period3_starts"] = per_ok
    out["lag1_fires"] = arms["lag1"]["fires"]
    out["lag1_fires_lag3_channel"] = arms["lag1"]["fires_lag3"]
    out["lag3_fires_total"] = arms["lag3"]["fires"]
    out["lag3_fires_lag3_channel"] = arms["lag3"]["fires_lag3"]
    out["lag3_fires_lag1_channel"] = arms["lag3"]["fires_lag1"]
    out["lag3_planted_starts_escaping_period3"] = len(escaped)
    out["invariant_violations"] = {
        "lag1": arms["lag1"]["invariant_violations"],
        "lag3": arms["lag3"]["invariant_violations"]}
    out["per_walker_steps"] = PER
    out["trace_periods_lag1"] = [period_of(arms["lag1"]["trace"], b) for b in planted[:24]]
    out["trace_periods_lag3"] = [period_of(arms["lag3"]["trace"], b) for b in planted[:24]]
    assert len(planted) > 0, "fixture planted no 3-cycle start"
    assert out["lag1_fires_lag3_channel"] == 0, out
    assert out["lag3_fires_lag3_channel"] >= len(planted), out
    assert out["invariant_violations"]["lag1"] == 0, out
    assert out["invariant_violations"]["lag3"] == 0, out
    assert out["lag3_planted_starts_escaping_period3"] == len(planted), out
    out["asserts_passed"] = True
    os.makedirs(os.path.join(D, "evidence"), exist_ok=True)
    with open(os.path.join(D, "evidence", "lag3-firing.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps({k: v for k, v in out.items() if k != "trace_periods_lag1"
                      and k != "trace_periods_lag3"}, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
