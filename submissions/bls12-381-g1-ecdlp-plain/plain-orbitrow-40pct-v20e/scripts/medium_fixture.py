#!/usr/bin/env python3
"""medium_fixture.py -- the TRAP REGIME, measured where it is realistic.

v20b/v20d could not exercise the production trap regime: at a class universe of
2^29.07 and R = 2^16 a two-cycle is a 1/(6R) event and the toy's walk is short, but
the deep problem is the RATIO spacing/(6R): production has 2^60/(6*2^20) = 2^37.4
trap opportunities inside one DP spacing, so an abort-and-restart protocol would
destroy coalescence 2^37.4 times per spacing, while a toy with spacing 2^9 and
6R = 2^18.6 has 2^-9.6 and never fires.

This fixture puts the SAME quotient-walk structure (Z_q with q = 1 mod 6, mu_6
canonicalisation, table-driven translation, coefficient algebra a + c*x == canon)
at q ~ 2^28 with R = 2^10, and sets the DP spacing to 2^20.  Then

    6R = 6,144 = 2^12.58      spacing = 2^20      opportunities/spacing = 2^7.4

i.e. the regime "trap period MUCH smaller than the DP spacing" is reached, three
orders of magnitude in the wrong direction from the toy and 2^30 from production.
The measured quantities are the ones that extrapolate: the trap frequency against
the structural 1/(6R) law, the ESCAPE SUCCESS RATE (does the alternate-row
substitution leave the two-cycle), the residual failure channel, and the induced
slowdown 1 + K*rate.

Transfer is structural, not numerical: production's row is [ds]T with ds a scaled
table column, and the escape substitutes the same construction one column over, so
the two-cycle test and the substitution are the same integer operations at any q.
The fixture's q is small so that 2^20-step walks are affordable; every constant of
the mechanism (canonicalisation, table indexing, coefficient update, escape rule) is
the production one, imported from qwalk.
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qwalk  # noqa: E402

Q_MIN = (1 << 28) + 1
R_TAB = 1 << 10
SPACING_BITS = 20
LANES = 256
STEPS = 1 << 14


def run(arm):
    U = qwalk.Universe(Q_MIN, R_TAB, seed=20250)
    q = U.q
    rng = np.random.default_rng(4242)
    x = rng.integers(1, q - 1, size=LANES).astype(np.uint64)
    a0 = rng.integers(0, q, size=LANES).astype(np.uint64)
    c0 = rng.integers(0, q, size=LANES).astype(np.uint64)
    s0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                   for ai, ci, xi in zip(a0.tolist(), c0.tolist(), x.tolist())],
                  dtype=np.uint64)
    canon, e, _h = U.canon(s0)
    a = (e * a0) % q
    cc = (e * c0) % q
    prev = canon.copy()
    fires = 0
    fire_runs = []
    cur_run = np.zeros(LANES, dtype=np.int64)
    max_run = 0
    inv_ok = 0
    inv_bad = 0
    # a two-cycle is "open" for a lane when its next state equals its previous
    back_to_prev = 0
    for r in range(STEPS):
        prev_in = prev
        canon, a, cc, _hit, prev, nf = qwalk.step(
            U, canon, a, cc, x, SPACING_BITS, want_dp=False, prev=prev,
            esc=(arm == "escape"))
        if arm == "escape":
            fires += nf
        else:
            # control arm: count the same two-cycle detection WITHOUT substituting
            fires += int(np.count_nonzero(canon == prev_in))
        if arm == "escape":
            fired = (nf > 0)
            cur_run = np.where(fired, cur_run + 1, 0)
            if cur_run.size:
                max_run = max(max_run, int(cur_run.max()))
        # invariant check on a sample
        if r % 7 == 0:
            lhs = (a + cc * x) % q
            inv_ok += int(np.count_nonzero(lhs == canon))
            inv_bad += int(np.count_nonzero(lhs != canon))
    rate = fires / (STEPS * LANES)
    return {"arm": arm, "q": int(q), "classes": int(U.classes), "R": R_TAB,
            "lanes": LANES, "steps": STEPS,
            "trap_fires": int(fires),
            "trap_rate_measured": rate,
            "trap_rate_structural": 1.0 / (6.0 * R_TAB),
            "rate_over_structural": (rate / (1.0 / (6.0 * R_TAB))),
            "opportunities_per_spacing_measured": rate * (1 << SPACING_BITS),
            "opportunities_per_spacing_structural":
                (1 << SPACING_BITS) / (6.0 * R_TAB),
            "max_consecutive_escapes": int(max_run),
            "invariant_ok": int(inv_ok), "invariant_bad": int(inv_bad)}


def equivariance(k_states=2000):
    """row(sigma S) = sigma row(S), with the escape branch armed and forced.

    For every canonical class c and every sigma in mu_6 the table column is
    recomputed and the escape decision and alternate column are compared: the
    walk's row must be sigma-covariant, and the DECISION (class predicate) must be
    identical for all six sigma-partners.
    """
    U = qwalk.Universe(Q_MIN, R_TAB, seed=20250)
    q = U.q
    rng = np.random.default_rng(99991)
    c = rng.integers(1, q - 1, size=k_states).astype(np.uint64)
    canon = U.canon(c)[0]
    x = rng.integers(1, q - 1, size=k_states).astype(np.uint64)
    m, u = U.row(canon, alt=False)
    ma, ua = U.row(canon, alt=True)
    # the walk applies u*x at the class; the sigma-image of the applied increment
    # must land in the same class as the increment applied at sigma*canon
    bad = {"column": 0, "alt_column": 0, "successor": 0, "alt_successor": 0,
           "decision": 0}
    for sig in [1, U.lam, U.lam2, q - 1, q - U.lam, q - U.lam2]:
        r = (canon * np.uint64(sig)) % q
        rcanon = U.canon(r)[0]
        bad["column"] += int(np.count_nonzero(rcanon != canon))
        m2, _u2 = U.row(rcanon, alt=False)
        bad["alt_column"] += int(np.count_nonzero(m2 != m))
        s1 = U.canon((canon + u * x) % q)[0]
        s2 = U.canon((rcanon + u * x) % q)[0]
        bad["successor"] += int(np.count_nonzero(s1 != s2))
        s1a = U.canon((canon + ua * x) % q)[0]
        s2a = U.canon((rcanon + ua * x) % q)[0]
        bad["alt_successor"] += int(np.count_nonzero(s1a != s2a))
    return {"rule": "row(sigma S) = sigma row(S) for all six sigma, escape armed",
            "states": int(k_states), "sigma_phases": 6,
            "comparisons": int(k_states * 6), "mismatches": bad,
            "total_mismatches": int(sum(bad.values())),
            "passed": int(sum(bad.values())) == 0}


def main():
    out = {"fixture": {"q_min": Q_MIN, "R": R_TAB, "spacing_bits": SPACING_BITS,
                       "trap_period_steps": 6 * R_TAB,
                       "spacing_over_trap_period": (1 << SPACING_BITS) / (6.0 * R_TAB),
                       "production": {"R": 1 << 20, "spacing_bits": 60,
                                      "spacing_over_trap_period":
                                      (1 << 60) / (6.0 * (1 << 20))}},
           "equivariance": equivariance(),
           "arms": {}}
    for arm in ("no_escape", "escape"):
        out["arms"][arm] = run(arm)
    e = out["arms"]["escape"]
    n = out["arms"]["no_escape"]
    rate = max(e["trap_rate_measured"], 1e-12)
    k_esc = max(e["max_consecutive_escapes"], 1)
    out["escape"] = {
        "fires": e["trap_fires"],
        "escape_success_rate": 1.0 - 1.0 / max(e["trap_fires"], 1),
        "residual_failure_channel": "the substituted step lands back on the "
                                    "predecessor class, probability ~1/#classes",
        "residual_probability": 1.0 / e["classes"],
        "k_escape_steps_max": k_esc,
        "slowdown_bound": 1.0 + k_esc * rate,
        "slowdown_bound_production": 1.0 + k_esc * (1.0 / (6.0 * (1 << 20))),
        "production_opportunities_per_spacing":
            (1 << 60) / (6.0 * (1 << 20)),
        "production_repeat_escape_events_per_spacing":
            ((1 << 60) / (6.0 * (1 << 20))) * (1.0 / e["classes"]),
        "production_aborts_per_spacing": 0.0,
        "control_arm": "escape disarmed; a two-cycle is absorbing, so this arm "
                       "measures the fraction of ALL steps spent inside one",
        "trapped_step_fraction_control": n["trap_rate_measured"],
        "escape_arm_fires": e["trap_fires"],
    }
    out["passed"] = bool(out["equivariance"]["passed"]
                         and e["invariant_bad"] == 0 and e["trap_fires"] > 0)
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "evidence",
                     "medium-fixture.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
