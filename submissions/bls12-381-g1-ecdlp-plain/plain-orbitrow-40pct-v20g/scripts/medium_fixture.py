#!/usr/bin/env python3
"""medium_fixture.py -- v20f: the TRAP REGIME and the CONVERSION, measured with
per-walk counters and MUTATION CONTROLS.

The three things v20e got wrong here, and what replaced them:

  (1) `force_every` was accepted by qwalk.step and never supplied by run(), so no
      forced-escape measurement was ever executed.  `run(arm, force_every)` now
      passes it, and the forced arm is a real measurement (it fires the
      substitution on every k-th lane-step and checks the invariant on the
      substituted branch).  The mutation arm FORGETS the coefficient update on
      the substituted branch and must produce invariant violations -- a test that
      cannot fail is not evidence.
  (2) `equivariance()` canonicalised both sides and then compared the identical
      canonical representatives, and its mismatch counter was never incremented.
      The test here is on PRE-canonical states and it tests the only non-trivial
      content: the column and the escape decision must be CLASS FUNCTIONS, so all
      six sigma-partners of a state must select the same column and take the same
      decision.  The control arm uses a REPRESENTATIVE-DEPENDENT selector (the
      v19d defect) and must fail with a large count.
  (3) `max_consecutive_escapes` counted batches (`nf > 0` across 256 lanes) and
      `escape_success_rate` was assigned as 1 - 1/fires.  Both are now counted
      PER WALK: fires[lane], the longest run of consecutive fire steps of one
      lane, and the number of fires whose substitution actually left the
      two-cycle (checked on the next step), plus a residual-trap counter.

Fixtures:
  trap   q ~ 2^28, R = 2^10 (trap period 6R = 6144), spacing 2^14, walks of 2^15
         steps = 2.0 full reporting spacings per walk (production: 2^37.4 trap
         opportunities per spacing -- a scaling of the same per-step integer
         operations, not a different mechanism).
  conversion  q ~ 2^24, all classes reported (d = 0), 24 walkers x 4096 steps:
         counts class collisions between walkers and how many of them solve x
         from two pool records, i.e. the conversion the reviewer asked for.
         Mutation arm stores a + 1 instead of a and must convert 0.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qwalk  # noqa: E402

Q_MIN = (1 << 28) + 1
R_TAB = 1 << 10
SPACING_BITS = 14
LANES = 256
STEPS = 1 << 15
SIGMAS = None


def _sig(U):
    q = U.q
    return [1, U.lam, U.lam2, q - 1, q - U.lam, q - U.lam2]


def run(arm, force_every=0, lanes=LANES, steps=STEPS, q_min=Q_MIN, R=R_TAB,
        spacing=SPACING_BITS, seed=4242, row_mode="orbit"):
    """Trap-regime arm.  arms: escape | no_escape | forced | forced_broken."""
    U = qwalk.Universe(q_min, R, seed=20250, row_mode=row_mode)
    q = U.q
    rng = np.random.default_rng(seed)
    x = rng.integers(1, q - 1, size=lanes).astype(np.uint64)
    a0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    c0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    s0 = (a0 + c0 * x) % q
    canon, e, _ = U.canon(s0)
    a = (e * a0) % q
    cc = (e * c0) % q
    prev = canon.copy()

    fires = np.zeros(lanes, dtype=np.int64)          # PER WALK
    runs = np.zeros(lanes, dtype=np.int64)           # current run, PER WALK
    maxrun = np.zeros(lanes, dtype=np.int64)         # PER WALK
    left_cycle = np.zeros(lanes, dtype=np.int64)     # substitutions that left the 2-cycle
    residual = np.zeros(lanes, dtype=np.int64)       # substitutions that did not
    trapped_steps = np.zeros(lanes, dtype=np.int64)  # steps whose successor == prev
    inv_ok = inv_bad = 0
    esc = arm in ("escape", "forced", "forced_broken")
    fe = force_every
    if arm == "forced":
        fe = force_every or 7
    if arm == "forced_broken":
        fe = force_every or 7

    for t in range(steps):
        prev_in = prev
        cur = canon.copy()          # class at the start of THIS step
        m = U.col(canon)
        idx = np.empty(0, dtype=np.int64)
        if arm == "forced_broken":
            # MUTATION: substitute the alternate row but forget the coefficient
            # update on the substituted branch (a real bug class in this package
            # family).  The invariant must break.
            cn, a2, cc2, _e = qwalk._apply(U, canon, a, cc, x, m)
            fire = cn == prev_in
            if fe:
                fire = fire | (((np.arange(lanes) + t) % fe) == 0)
            idx = np.flatnonzero(fire)
            if idx.size:
                ma = (m[idx] + 1) % U.R
                cn2, _a3, _cc3, _e2 = qwalk._apply(U, canon[idx], a[idx], cc[idx],
                                                   x[idx], ma)
                cn = cn.copy(); a2 = a2.copy(); cc2 = cc2.copy()
                cn[idx] = cn2                      # element moved, coefficients NOT
                fired = fire
            else:
                fired = fire
        else:
            cn, a2, cc2, _hit, prev_out, fired = qwalk.step(
                U, canon, a, cc, x, spacing, want_dp=False, prev=prev_in,
                esc=esc, force_every=fe, tick=t)
        fires += fired
        if esc:
            runs = np.where(fired, runs + 1, 0)
            np.maximum(maxrun, runs, out=maxrun)
        # a two-cycle visit: successor class == the class one step back
        back = cn == prev_in
        trapped_steps += back
        # did each fire actually leave the two-cycle?  A fire at step t is
        # successful iff the substituted successor is NOT the predecessor class.
        if fired.any():
            idxf = np.flatnonzero(fired)
            left = cn[idxf] != prev_in[idxf]
            left_cycle[idxf] += left
            residual[idxf] += ~left
        if t % 7 == 0:
            lhs = (a2 + cc2 * x) % q
            inv_ok += int(np.count_nonzero(lhs == cn))
            inv_bad += int(np.count_nonzero(lhs != cn))
        canon, a, cc, prev = cn, a2, cc2, cur

    lane_steps = steps * lanes
    total_fires = int(fires.sum())
    return {
        "arm": arm, "row_mode": row_mode, "force_every": int(fe),
        "q": int(q), "classes": int(U.classes), "R": R,
        "lanes": lanes, "steps": steps, "lane_steps": int(lane_steps),
        "spacing": 1 << spacing, "spacings_per_walk": (1 << spacing) / float(steps),
        "trap_period_steps": 6 * int(U.R_eff),
        "R_eff": int(U.R_eff),
        "trap_opportunities_per_spacing": (1 << spacing) / float(6 * U.R_eff),
        "fires": total_fires,
        "fire_rate_per_step": total_fires / float(lane_steps),
        "structural_1_over_6R": 1.0 / (6.0 * U.R_eff),
        "fire_rate_over_structural": (total_fires / float(lane_steps)) / (1.0 / (6.0 * U.R_eff)),
        "walks_with_fire": int((fires > 0).sum()),
        "max_consecutive_fires_per_walk": int(maxrun.max()),
        "fires_left_two_cycle": int(left_cycle.sum()),
        "fires_residual": int(residual.sum()),
        "escape_success_rate_measured": (float(left_cycle.sum()) / total_fires
                                         if total_fires else 0.0),
        "trapped_step_fraction": int(trapped_steps.sum()) / float(lane_steps),
        "invariant_ok": int(inv_ok), "invariant_bad": int(inv_bad),
    }


def equivariance(k=1500, q_min=Q_MIN, R=R_TAB, seed=99991, row_mode="orbit",
                 selector="class"):
    """CLASS-FUNCTION test on PRE-canonical states + the mutation control.

    For a random raw state S (NOT canonicalised first) and its six sigma-images:
      * the selected column must be identical (the selector is a class function);
      * the escape DECISION must be identical: the successor class equals the
        previous class's orbit, tested with RAW previous states so the test would
        fail on an implementation that compared representatives;
      * the substituted successor class must be identical.
    selector="raw" is the MUTATION: a representative-dependent column.  The same
    test must fail with a large count.
    """
    U = qwalk.Universe(q_min, R, seed=20250, row_mode=row_mode)
    q = U.q
    rng = np.random.default_rng(seed)
    S = rng.integers(1, q - 1, size=k).astype(np.uint64)
    P = rng.integers(1, q - 1, size=k).astype(np.uint64)   # raw previous states
    x = rng.integers(1, q - 1, size=k).astype(np.uint64)
    bad = {"column": 0, "decision": 0, "alt_successor": 0}
    for sig in _sig(U):
        Rs = (S * np.uint64(sig)) % q
        Rp = (P * np.uint64(sig)) % q
        if selector == "class":
            m1 = U.col(U.canon(S)[0])
            m2 = U.col(U.canon(Rs)[0])
        else:                                   # MUTATION: raw representative
            m1 = np.minimum((S.astype(np.int64) * R) // q, R - 1)
            m2 = np.minimum((Rs.astype(np.int64) * R) // q, R - 1)
        bad["column"] += int(np.count_nonzero(m1 != m2))
        # decision: does the successor land on the previous class?  computed on
        # RAW states and compared through the class, never through the value.
        d1 = U.canon(U.canon((U.canon(S)[0] + U.u[m1] * x) % q)[0])[0] == U.canon(P)[0]
        d2 = U.canon((U.canon(Rs)[0] + U.u[m2] * x) % q)[0] == U.canon(Rp)[0]
        bad["decision"] += int(np.count_nonzero(d1 != d2))
        s1 = U.canon((U.canon(S)[0] + U.u[(m1 + 1) % R] * x) % q)[0]
        s2 = U.canon((U.canon(Rs)[0] + U.u[(m2 + 1) % R] * x) % q)[0]
        bad["alt_successor"] += int(np.count_nonzero(s1 != s2))
    total = sum(bad.values())
    return {"selector": selector, "row_mode": row_mode, "states": int(k),
            "sigma_phases": 6, "comparisons": int(k * 6), "mismatches": bad,
            "total_mismatches": int(total), "passed": total == 0}


def conversion(lanes=24, steps=4096, q_min=(1 << 24) + 1, R=R_TAB, seed=555,
               corrupt=False):
    """Class collision -> verified duplicate pool record, measured end to end.

    Every class is reported (d = 0), the pool is keyed by the class, and the two
    records of a colliding class solve x = (a1 - a2)/(c2 - c1) mod q.  Records
    are produced by different walkers at their own steps: no time alignment.
    corrupt=True stores a + 1 instead of a (MUTATION): the conversion must be 0.
    """
    U = qwalk.Universe(q_min, R, seed=20250)
    q = U.q
    rng = np.random.default_rng(seed)
    # ONE challenge per trial: the pool recovery is x = (a1-a2)/(c2-c1) for the
    # SAME x, so per-lane challenges make the conversion unmeasurable (v20f's
    # fixture did that and read 0/373 while claiming end-to-end conversion).
    x0 = int(rng.integers(1, q - 1))
    x = np.full(lanes, x0, dtype=np.uint64)
    a0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    c0 = rng.integers(0, q, size=lanes).astype(np.uint64)
    s0 = (a0 + c0 * x) % q
    canon, e, _ = U.canon(s0)
    a = (e * a0) % q
    cc = (e * c0) % q
    prev = canon.copy()
    pool = {}
    for t in range(steps):
        cn, a2, cc2, _hit, _pv, _f = qwalk.step(U, canon, a, cc, x, 0, want_dp=False,
                                                prev=prev, esc=True, tick=t)
        for i in range(lanes):
            key = int(cn[i])
            av = int(a2[i])
            # MUTATION: corrupt the COEFFICIENT, not the numerator -- adding 1 to
            # a on both records of a class cancels in (a1-a2) and would leave the
            # conversion at 100%, a mutation that cannot fail.
            # +1 on both records of a class cancels in (c2-c1), so the corruption
            # must be per record: walker-index-dependent.
            cv = (int(cc2[i]) + i + 1) % q if corrupt else int(cc2[i])
            pool.setdefault(key, []).append((av, cv, i))
        canon, a, cc, prev = cn, a2, cc2, canon.copy()
    xx = int(x[0])
    collisions = 0
    verified = 0
    for key, recs in pool.items():
        ws = {r[2] for r in recs}
        if len(recs) < 2 or len(ws) < 2:
            continue
        collisions += 1
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                if recs[i][2] == recs[j][2]:
                    continue
                den = (recs[j][1] - recs[i][1]) % q
                if den == 0:
                    continue
                cand = ((recs[i][0] - recs[j][0]) * pow(den, -1, q)) % q
                if cand == xx:
                    verified += 1
                    break
            else:
                continue
            break
    return {"lanes": lanes, "steps": steps, "q": int(q), "classes": int(U.classes),
            "records": lanes * steps, "corrupt": bool(corrupt),
            "collisions": collisions, "verified": verified,
            "conversion": (verified / collisions if collisions else 0.0),
            "expected_collisions": (lanes * steps) ** 2 / (2.0 * U.classes)}


def main():
    out = {"fixture": {"q_min": Q_MIN, "R": R_TAB, "spacing_bits": SPACING_BITS,
                       "steps": STEPS, "lanes": LANES,
                       "spacings_per_walk": (1 << SPACING_BITS) / float(STEPS),
                       "production": {"R": 1 << 20, "spacing_bits": 60,
                                      "trap_opportunities_per_spacing":
                                      (1 << 60) / (6.0 * (1 << 20))}}}
    arms = {}
    arms["no_escape"] = run("no_escape")
    arms["escape"] = run("escape")
    arms["forced"] = run("forced", force_every=7, steps=STEPS // 8)
    arms["forced_broken"] = run("forced_broken", force_every=7, steps=STEPS // 8)
    out["arms"] = arms
    out["equivariance"] = equivariance()
    out["equivariance_mutation"] = equivariance(selector="raw")
    out["conversion"] = conversion()
    out["conversion_mutation"] = conversion(corrupt=True)
    e, n, fb = arms["escape"], arms["no_escape"], arms["forced_broken"]
    out["checks"] = {
        "trap_absorbing_without_escape": n["trapped_step_fraction"] > 0.10,
        "escape_removes_trapping": e["trapped_step_fraction"] < 0.01 * max(
            n["trapped_step_fraction"], 1e-9) or e["trapped_step_fraction"] == 0.0,
        "escape_leaves_two_cycle": e["fires_left_two_cycle"] > 0 and e["fires_residual"] == 0,
        "per_walk_fire_counts": e["walks_with_fire"] > 1 and e["max_consecutive_fires_per_walk"] >= 1,
        "forced_arm_measured": arms["forced"]["fires"] > 0,
        "forced_mutation_breaks_invariant": fb["invariant_bad"] > 0,
        "equivariance_passes": out["equivariance"]["passed"],
        "equivariance_mutation_fails": not out["equivariance_mutation"]["passed"],
        "conversion_verified": out["conversion"]["verified"] > 0,
        "conversion_mutation_zero": out["conversion_mutation"]["verified"] == 0,
    }
    out["passed"] = all(out["checks"].values())
    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "evidence",
                     "medium-fixture.json")
    with open(p, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps({k: out[k] for k in ("checks", "passed")}, indent=1))
    for k, v in out["arms"].items():
        print(k, {kk: v[kk] for kk in ("fires", "walks_with_fire",
                                       "max_consecutive_fires_per_walk",
                                       "fires_left_two_cycle", "fires_residual",
                                       "escape_success_rate_measured",
                                       "trapped_step_fraction", "invariant_bad",
                                       "fire_rate_over_structural")})
    print("equiv", out["equivariance"]["total_mismatches"],
          "mutation", out["equivariance_mutation"]["total_mismatches"])
    print("conv", out["conversion"]["collisions"], out["conversion"]["verified"],
          "mut", out["conversion_mutation"]["verified"])


if __name__ == "__main__":
    main()
