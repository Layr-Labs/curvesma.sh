#!/usr/bin/env python3
"""cycle3.py -- T-D: the THREE-CYCLE channel, censused like the two-cycle.

The v20g reviewer's point 4: the two-step-return test cannot see a cycle of
length three.  For lambda^3 = 1 the affine map f(z) = lambda*(z+v) satisfies
f^3(z) = z because 1 + lambda + lambda^2 = 0 and lambda^3 = 1, so a walk that
takes THAT map three times in a row is periodic and the two-step test
(successor == previous class) never fires on it.  v20g neither ruled this out
nor bounded it.

Census, on the same trap-regime universe as the medium fixture (q ~ 2^28,
R = 2^10, 256 lanes, 2^15 steps), escape ARMED -- and the v20e-era NO-ESCAPE
configuration as the control that shows its own cycles:

  one/two/three-cycle visits   c_{t+1} == c_t / c_{t+1} == c_{t-1} /
                               c_{t+1} == c_{t-2} (and neither of the shorter ones)
  STRUCTURAL three-cycles      the last AND the same column on all three steps
                               (m_t == m_{t-1} == m_{t-2}): the reviewer's channel,
                               where equal u lets the translation cancel.  A cycle
                               with unequal columns is not the lambda*(z+v) map.
  settled cycle lengths        smallest L in 1..64 with c_t == c_{t-L} over the last
                               64 steps, per lane: which cycles the control falls into.

The escape is armed through qwalk.step with the K_ESC run-length array, so this
census also re-exercises the enforced cap.

Writes evidence/cycle3.json.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
import qwalk               # noqa: E402

Q_MIN = (1 << 28) + 1
R_TAB = 1 << 10
LANES = 256
STEPS = 1 << 15
SPACING_BITS = 14
TAIL = 128


def census(esc, lanes=LANES, steps=STEPS, seed=4242):
    U = qwalk.Universe(Q_MIN, R_TAB, seed=20250)
    q = U.q
    rng = np.random.default_rng(seed)
    x = rng.integers(1, q - 1, size=lanes).astype(np.uint64)
    a0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    c0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    s0 = (a0 + c0 * x) % q
    c1, e, _ = U.canon(s0)          # c_t
    a = (e * a0) % q
    cc = (e * c0) % q
    c2 = c1.copy()                   # c_{t-1}
    c3 = c1.copy()                   # c_{t-2}
    c4 = c1.copy()                   # c_{t-3}
    run = np.zeros(lanes, dtype=np.int64)
    one = two = three = three_struct = three_diff = three_twoeq = three_allne = 0
    fires = 0
    kcap_max = 0
    tail = np.empty((TAIL, lanes), dtype=np.uint64)
    for t in range(steps):
        tail[t % TAIL] = c1
        m1 = U.col(c1)
        cn, a, cc, _hit, pv, fired = qwalk.step(
            U, c1, a, cc, x, SPACING_BITS, want_dp=False, prev=c2,
            esc=esc, esc_run=run, kcap=(qwalk.K_ESC if esc else None))
        kcap_max = max(kcap_max, int(run.max(initial=0)))
        fires += int(np.count_nonzero(fired))
        one += int(np.count_nonzero(cn == c1))
        two += int(np.count_nonzero(cn == c2))
        hit3 = (cn == c3) & (cn != c2) & (cn != c1)
        n3 = int(np.count_nonzero(hit3))
        if n3:
            three += n3
            m2 = U.col(c2)
            m3 = U.col(c3)
            same = hit3 & (m1 == m2) & (m2 == m3)
            two_eq = hit3 & (((m1 == m2) & (m2 != m3))
                             | ((m1 == m3) & (m1 != m2))
                             | ((m2 == m3) & (m1 != m2)))
            three_struct += int(np.count_nonzero(same))
            three_twoeq += int(np.count_nonzero(two_eq))
            three_allne += n3 - int(np.count_nonzero(same)) - int(np.count_nonzero(two_eq))
        c4, c3, c2, c1 = c3, c2, c1, cn
    # settled cycle length from the tail window (chronological order restored)
    idx = np.arange(steps - TAIL, steps) % TAIL
    T = tail[idx]                      # tail window, oldest -> newest
    lens = np.zeros(lanes, dtype=np.int64)
    seg = T[TAIL - 64:]
    for L in range(1, 65):
        ref = T[TAIL - 64 - L:TAIL - L]
        hitL = np.all(seg == ref, axis=0)
        lens[(lens == 0) & hitL] = L
    lens[lens == 0] = 999              # no period <= 64 in the window
    vals, cnts = np.unique(lens, return_counts=True)
    lane_steps = lanes * steps
    return {
        "esc": bool(esc), "q": int(q), "classes": int(U.classes),
        "R": int(R_TAB), "R_eff": int(U.R_eff), "lanes": int(lanes),
        "steps": int(steps), "lane_steps": int(lane_steps),
        "fires": int(fires),
        "kcap_enforced": int(qwalk.K_ESC) if esc else 0,
        "kcap_max_seen": int(kcap_max),
        "one_cycle_visits": int(one), "two_cycle_visits": int(two),
        "three_cycle_visits": int(three),
        "three_cycle_structural": int(three_struct),
        "three_cycle_two_columns_equal": int(three_twoeq),
        "three_cycle_all_columns_distinct": int(three_allne),
        "two_cycle_rate_per_step": two / float(lane_steps),
        "three_cycle_rate_per_step": three / float(lane_steps),
        "settled_cycle_lengths": {int(v): int(c) for v, c in
                                  zip(vals.tolist(), cnts.tolist())},
    }


def main():
    esc = census(True)
    ctl = census(False)
    res = {
        "test": "three_cycle_census",
        "fixture": {"q_min": Q_MIN, "R": R_TAB, "lanes": LANES, "steps": STEPS,
                    "spacing_bits": SPACING_BITS,
                    "desc": "same trap regime as evidence/medium-fixture.json"},
        "escape_armed": esc,
        "control_no_escape": ctl,
        "checks": {
            "control_settles_into_cycles": any(
                int(k) <= 64 for k in ctl["settled_cycle_lengths"]),
            "control_has_cycles_longer_than_two": any(
                int(k) >= 3 for k in ctl["settled_cycle_lengths"] if int(k) <= 64),
            "structural_three_cycle_absent_with_escape":
                esc["three_cycle_structural"] == 0,
            "control_structural_three_cycles":
                ctl["three_cycle_structural"],
            "escape_kcap_never_tripped": esc["kcap_max_seen"] <= qwalk.K_ESC,
            "escape_two_cycle_collapses": esc["two_cycle_rate_per_step"]
                < 0.01 * max(ctl["two_cycle_rate_per_step"], 1e-12),
        },
    }
    res["passed"] = all([res["checks"][k] for k in (
        "control_settles_into_cycles", "control_has_cycles_longer_than_two",
        "structural_three_cycle_absent_with_escape", "escape_kcap_never_tripped",
        "escape_two_cycle_collapses")])
    out = os.path.join(D, "evidence", "cycle3.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
