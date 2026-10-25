"""v21 return/cycle measurement -- CORRECTED after ticket 90e73945.

The review is right and v20 was wrong.  R is a power of two dividing 2^64, so
    prod = (S XOR (t*2654435761 mod 2^64) XOR trail*7919) mod R
satisfies prod(t + R) = prod(t): the EFFECTIVE phase is t mod R even though the
stored counter grows without bound.  A repeated true point at the same effective
phase therefore repeats the whole subsequent point dynamics, and v20's claim that
"the pair (S,t) never recurs" was false -- it tested the raw counter, not the phase.

What this script measures instead:

* state repeats keyed by (TRUE point T, t mod R) -- the correct recurrence unit;
* each repeat classified: USEFUL (coefficient difference nonzero, so the loop
  recovers x and the trial has already succeeded) or FATAL (coefficient difference
  (0,0), the walker cycles with no new information);
* a HIGH-LOAD cell where phase-matched repeats are frequent, so the fatal fraction
  is measured rather than asserted.  The load lambda = V^2/(2*classes) is large
  there by construction: the expected repeat count is computed analytically BEFORE
  the run and printed, which is the non-vacuity check the review demanded;
* the production charge: repeats per trail at production lambda, times the measured
  fatal fraction, compared with the CP-95 upper bound on the per-trail restart rate
  that the package charges (the larger of the two is charged).
"""
import json, math, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")
sys.path.insert(0, HERE)
import walk_v20 as W20


def cp95_upper(k, n):
    return W20.cp95_upper(k, n)


def walk_tracked(q, R, STEPS, x, tab, c, d, S, trail):
    """Declared walk with phase-aware exact tracking on the TRUE point."""
    tmul = 2654435761
    seen = {}
    reps = 0
    fatal = 0
    useful = 0
    ident = 0
    for t in range(STEPS):
        prod = (S ^ ((t * tmul) & 0xFFFFFFFFFFFFFFFF) ^ (trail * 7919)) % R
        v, a, b = tab[prod]
        T = S + v
        if T >= q:
            T -= q
        if T == 0:
            ident += 1
            break
        c2 = c + a
        if c2 >= q:
            c2 -= q
        d2 = d + b
        if d2 >= q:
            d2 -= q
        S = T
        c = c2
        d = d2
        key = (S, t % R)                     # TRUE point + EFFECTIVE phase
        prev = seen.get(key)
        if prev is None:
            seen[key] = (c, d)
        else:
            reps += 1
            pc, pd = prev
            da = (c - pc) % q
            db = (d - pd) % q
            if da == 0 and db == 0:
                fatal += 1
            else:
                cand = (-da * pow(db, -1, q)) % q if db else None
                if cand == x:
                    useful += 1
                else:
                    fatal += 1               # unresolvable repeat: treat as a loss
    return {"repeats": reps, "fatal": fatal, "useful": useful,
            "identity_breaks": ident, "visits": STEPS}


def main():
    # high-load cell: expected phase-matched repeats >> 1 by construction
    q0, R0, W0, S0, n0 = 1048573, 64, 1, 8192, 200      # q ~ 2^20, R = 64
    visits = W0 * S0
    expected = visits * (visits - 1) / (2.0 * q0) / R0   # phase-matched, analytic
    tot = {"repeats": 0, "fatal": 0, "useful": 0, "identity_breaks": 0, "visits": 0}
    trails = 0
    for ti in range(n0):
        rng = random.Random((20251024 << 26) ^ (ti * W20.GOLDEN) ^ (R0 << 3) ^ (W0 << 50))
        x = rng.randrange(2, q0 - 1)
        tab = W20.make_table(q0, x, R0)
        for w in range(W0):
            c = rng.randrange(1, q0)
            d = rng.randrange(1, q0)
            Sv = (c + d * x) % q0
            r = walk_tracked(q0, R0, S0, x, tab, c, d, Sv, (ti * W0 + w) & 0xFFFFFF)
            trails += 1
            for k in tot:
                tot[k] += r[k]
    fatal_frac = tot["fatal"] / tot["repeats"] if tot["repeats"] else None
    fatal_ub = cp95_upper(tot["fatal"], tot["repeats"]) if tot["repeats"] else 1.0
    assert tot["repeats"] > expected * 0.3, (tot["repeats"], expected)

    # production: phase-matched repeats per trail at lambda = 0.634's own geometry
    LOG2_V, WALKERS = 127.67, 1 << 20
    S_prod = 2 ** (LOG2_V - 20)                 # steps per walker
    q_prod = 2 ** 256
    R_prod = 1 << 20
    reps_per_trail_prod = S_prod * (S_prod - 1) / (2.0 * q_prod * R_prod)
    charge_analytic = reps_per_trail_prod * fatal_ub

    prev = json.load(open(os.path.join(OUT, "returns_v20.json")))
    charge_measured_restart = prev["identity_break_rate_cp95_upper"]
    charge = max(charge_analytic, charge_measured_restart)

    out = {
        "correction": {
            "retracted": "v20 asserted that (S,t) never recurs because t strictly increases. "
                         "That is FALSE: R is a power of two dividing 2^64, so the selector is "
                         "R-periodic in t and the effective phase is t mod R.",
            "therefore": "a repeated TRUE point at the same effective phase repeats the "
                         "subsequent point dynamics, exactly as the review says",
            "fixed_here": "tracking keys on (true point, t mod R); every repeat is classified "
                          "USEFUL (nonzero coefficient difference: the loop recovers x) or "
                          "FATAL ((0,0) difference, or an unresolvable relation)",
        },
        "high_load_cell": {
            "why_non_vacuous": "expected phase-matched repeats computed analytically BEFORE the "
                               "run and printed below; the cell is chosen so the expectation is "
                               "large",
            "q": q0, "R": R0, "walkers": W0, "steps_per_walker": S0, "trials": n0,
            "expected_repeats_analytic": expected,
            "measured": tot, "trails": trails,
            "fatal_fraction": fatal_frac,
            "fatal_fraction_cp95_upper": fatal_ub,
            "detects_repeats": tot["repeats"] > 0,
        },
        "production_charge": {
            "repeats_per_trail_analytic": reps_per_trail_prod,
            "formula": "S(S-1)/(2*q*R), S = 2^107.67 steps, q = 2^256, R = 2^20",
            "fatal_fraction_cp95_upper": fatal_ub,
            "analytic_cycle_loss_charge": charge_analytic,
            "measured_restart_charge": charge_measured_restart,
            "charged_relative_loss": charge,
            "basis": "the LARGER of the analytic cycle charge and the per-trail CP-95 upper "
                     "bound on identity restarts from the v20 stability runs; the loss model is "
                     "per trail (a dead walker loses its own work), never divided by 2^20",
        },
        "v20_error": prev["pricing"]["v19_error"],
        "note": "this row is a correctness fix, not a price reduction: the charged loss is the "
                "same order as v20's, and the justification is now the one the review asked for",
    }
    json.dump(out, open(os.path.join(OUT, "returns_v21.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({k: out[k] for k in ("high_load_cell", "production_charge")}, indent=1)[:1600])


if __name__ == "__main__":
    main()
