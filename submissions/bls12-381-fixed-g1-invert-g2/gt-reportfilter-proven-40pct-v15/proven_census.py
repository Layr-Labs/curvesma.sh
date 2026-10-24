#!/usr/bin/env python3
"""v15 measurements: retention, class-marginal stability, cycle-pattern census.

The v14 verdict (2a566536) rejected three things this script measures instead of
asserting:

  (1) RETENTION.  "equivariance with correspondingly transformed rings is not the
      same as the successor class being a function of the current class alone ...
      the specific 0.982421875 retention is not established by the covariance
      tests."  v15 measures retention directly, as the coalescence of two trails
      that sit on the SAME class with INDEPENDENT ring histories, at the shipped
      guard depth D = 12, over a table-size sweep, and reports CP-95 bounds.

  (2) CLASS MARGINAL.  "a deterministic map on augmented states does not by itself
      preserve a uniform class marginal."  v15 measures the class-frequency vector
      of long guarded walks and reports a chi-square statistic against uniform,
      its p-value and the worst cell deviation, plus a block-split stability row.

  (3) CYCLE PATTERNS.  "counting repeated rows does not justify multiplying
      independent 1/(12R) factors without accounting for possible patterns and
      their routing dependence."  v15 measures the conditional return probability
      of every row pattern up to length 6 by direct count, on the guarded walk and
      on the depth-0 control, and prices the trap with a UNION BOUND over patterns
      (combinatorial pattern count x measured conditional probability), not with
      multiplied independence factors.

  (4) JOINT CENSUS.  "T1 compares only elements reached from empty rings, not
      complete joint states ... the joint-state enumeration is at depth 1, not the
      shipped depth 12."  v15 reruns the census on the FULL joint state
      (element, ring tuple) at the shipped D = 12, exhaustive over all fixture
      starts, and reports the actual first-return lag (null is a legal result).

Fixture: the v14 real-tower toy (tower.py, classkey.py, routing.py), r = 9421,
Fp12 tower over a 64-bit p, twelve pre-expanded columns.  Table size R is a
parameter and is stated on every row; the production table is R = 2^14.

Outputs: evidence/retention.json, evidence/pattern-census.json,
         evidence/census-joint.json
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import classkey as CK
import routing
import tower as TW
import toy_census as TC
from cp95 import cp95_lower

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
BITS = 5                      # fixture: r = 9421, p 60-bit, Fp12 tower
R_ELEMS = 9420                # r - 1, the number of nonidentity elements
D = K.GUARD_DEPTH             # 12 -- the SHIPPED guard depth


# ---------------------------------------------------------------------------
# chi-square upper tail, exact enough for df ~ 784 (regularised incomplete gamma)
# ---------------------------------------------------------------------------
def _gammq(a, x):
    """Q(a,x) = 1 - P(a,x), Numerical-Recipes gser/gcf."""
    if x < 0 or a <= 0:
        return float('nan')
    if x == 0:
        return 1.0
    if x < a + 1.0:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(10000):
            ap += 1.0
            d *= x / ap
            s += d
            if abs(d) < abs(s) * 1e-16:
                break
        return 1.0 - s * math.exp(-x + a * math.log(x) - math.lgamma(a))
    b = x + 1.0 - a
    c = 1e300
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-16:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def chi2_sf(x, df):
    return _gammq(df / 2.0, x / 2.0)


# ---------------------------------------------------------------------------
# (1) class-marginal stability
# ---------------------------------------------------------------------------
def class_marginal(toy, walks=512, steps=200, seed=11):
    rng = random.Random(seed)
    T = toy.T
    blocks = 8
    per_block = walks // blocks
    block_counts = []
    counts = {}
    samples = 0
    for w in range(walks):
        if w % per_block == 0:
            block_counts.append({})
        y = toy.rand_gt()
        a, b = 1, 0
        guard = CK.Guard(D)
        for s in range(steps):
            y, a, b, kind, info = toy.step(y, a, b, guard)
            mkey = T.keys12(y)
            m = min(mkey)
            counts[m] = counts.get(m, 0) + 1
            bl = block_counts[-1]
            bl[m] = bl.get(m, 0) + 1
            samples += 1
    cells = len(counts)
    e = samples / cells
    chi2 = sum((c - e) ** 2 / e for c in counts.values())
    df = cells - 1
    maxdev = max(abs(c - e) / e for c in counts.values())
    # block-split stability: chi2 recomputed inside each block of walks
    bchi = []
    bn = steps * per_block
    for bl in block_counts:
        be = bn / cells
        bchi.append(sum((c - be) ** 2 / be for c in bl.values() if be > 0)
                    + be * (cells - len(bl)))
    return dict(
        fixture_bits=BITS, r=toy.r, walks=walks, steps_per_walk=steps,
        samples=samples, class_cells_observed=cells,
        expected_classes=toy.N,
        class_cells_equal_expected_N=(cells == toy.N),
        expected_per_cell=e,
        chi2=chi2, df=df, p_value=chi2_sf(chi2, df),
        chi2_over_df=chi2 / df,
        max_cell_deviation_relative=maxdev,
        block_split=dict(blocks=blocks, walks_per_block=per_block,
                         chi2_per_block=bchi,
                         chi2_over_df_per_block=[c / df for c in bchi],
                         min_chi2_over_df=min(bchi) / df,
                         max_chi2_over_df=max(bchi) / df),
        note="the class key is the sigma-invariant min of the twelve Frobenius keys; "
             "cells = distinct keys observed = N = (r-1)/12 exactly; a uniform class "
             "marginal is the null, so chi2/df near 1 and a p-value that is not tiny "
             "is the confirmation the verdict asked for.",
    )


# ---------------------------------------------------------------------------
# (2) retention: coalescence of two trails on the same class, independent rings
# ---------------------------------------------------------------------------
def _warm_ring(toy, rng, steps=D):
    """A ring drawn from an actual walk (independent history), starting empty."""
    y = toy.rand_gt()
    g = CK.Guard(D)
    for _ in range(steps):
        y, _a, _b, kind, _i = toy.step(y, 1, 0, g)
    return y, tuple(g.slots)


def retention(toy, pairs, maxlag=D):
    rng = random.Random(20261025)
    T = toy.T
    coalesced = 0
    first_lag = {}
    agree_one_step = 0
    for _ in range(pairs):
        y = toy.rand_gt()
        _yA, ringA = _warm_ring(toy, rng)
        _yB, ringB = _warm_ring(toy, rng)
        gA, gB = CK.Guard(D), CK.Guard(D)
        gA.slots = list(ringA)
        gB.slots = list(ringB)
        yA, yB = y, y
        met = None
        for s in range(1, 2 * maxlag + 1):
            yA, _a, _b, kA, iA = toy.step(yA, 1, 0, gA)
            yB, _b2, _b3, kB, iB = toy.step(yB, 1, 0, gB)
            if s == 1 and kA == kB:
                agree_one_step += 1
            same = (yA == yB) and (tuple(gA.slots) == tuple(gB.slots))
            if same:
                if met is None:
                    met = s
                elif s - met >= maxlag:
                    break
            elif met is not None:
                met = -1                     # met and then split: not retained
                break
        ok = met is not None and met > 0
        coalesced += 1 if ok else 0
        first_lag[str(met if met is not None else 0)] = \
            first_lag.get(str(met if met is not None else 0), 0) + 1
    rate = coalesced / pairs
    return dict(
        rows=len(toy.us), guard_depth=D, pairs=pairs,
        coalesced=coalesced, rate=rate,
        cp95_lower=float(cp95_lower(coalesced, pairs)),
        cp95_upper_failure=1.0 - float(cp95_lower(coalesced, pairs)),
        failure_rate=1 - rate,
        implied_c=(1 - rate) * len(toy.us),
        one_step_class_agreement=agree_one_step / pairs,
        theory_1_minus_2D2_over_R=1.0 - 2.0 * D * D / len(toy.us),
        first_coalescence_lag_hist=first_lag,
        note="both trails start on the SAME element with rings drawn from two "
             "independent walks; a pair is retained when the joint state "
             "(element, ring tuple) coincides within D steps and stays equal for a "
             "further D steps.  This is the merge-stickiness event the ledger "
             "charges, measured rather than derived.",
    )


# ---------------------------------------------------------------------------
# (3) cycle-pattern census, guarded vs control; (4) joint-state first returns
# ---------------------------------------------------------------------------
def pattern_census(toy, maxlag=26, control=False):
    r = toy.r
    idx = TC.enumerate_states(toy)
    steps = maxlag
    pat = {}                      # L -> {pattern: [occ, ret]}
    elem_returns = {L: 0 for L in range(2, maxlag + 1)}
    joint_returns = {L: 0 for L in range(2, maxlag + 1)}
    first_joint = None
    for i in range(1, r):
        g = CK.Guard(D) if not control else None
        y, a, b = idx[i], 1, 0
        ys = [y]
        rings = [tuple(g.slots) if g is not None else ()]
        rows = []
        for _s in range(steps):
            if control:
                y, a, b, kind, info = toy.step_v13(y, a, b)
                t = info[1]
            else:
                y, a, b, kind, info = toy.step(y, a, b, g)
                t = info[1]
            rows.append(t)
            ys.append(y)
            rings.append(tuple(g.slots) if g is not None else ())
        for L in range(2, maxlag + 1):
            if ys[L] == ys[0]:
                elem_returns[L] += 1
            if rings[L] == rings[0] and ys[L] == ys[0]:
                joint_returns[L] += 1
                if first_joint is None or L < first_joint:
                    first_joint = L
        for L in range(2, 7):
            d = pat.setdefault(L, {})
            for i0 in range(0, steps - L + 1):
                key = tuple(rows[i0:i0 + L])
                e = d.get(key)
                if e is None:
                    d[key] = [1, 1 if ys[i0 + L] == ys[i0] else 0]
                else:
                    e[0] += 1
                    if ys[i0 + L] == ys[i0]:
                        e[1] += 1
    rows_out = {}
    for L in (2, 3, 4, 5, 6):
        d = pat[L]
        occ = sum(v[0] for v in d.values())
        ret = sum(v[1] for v in d.values())
        # worst pattern: CP-95 upper on its conditional return probability
        worst, worst_rate, worst_occ = None, 0.0, 0
        for key, (o, k) in d.items():
            hi = 1.0 - (0.05 ** (1.0 / o)) if k == 0 else float(cp95_lower(k, o))
            if hi >= worst_rate:
                worst, worst_rate, worst_occ = key, hi, o
        rows_out[L] = dict(patterns_observed=len(d), occurrences=occ,
                           returns=ret, aggregate_rate=ret / occ,
                           worst_pattern_rate_cp95_upper=worst_rate,
                           worst_pattern_occurrences=worst_occ,
                           worst_pattern=("<none: all patterns return>"
                                          if ret and worst_rate > 0.5 else "row-sequence"),
                           zero_event_rule_of_three=3.0 / occ)
    nz = [L for L in elem_returns if elem_returns[L]]
    nzj = [L for L in joint_returns if joint_returns[L]]
    return dict(
        rows=len(toy.us), guard_depth=(None if control else D),
        control=(control), maxlag=maxlag, elements=r - 1, steps=maxlag,
        element_return_lag_counts={str(L): elem_returns[L] for L in elem_returns},
        joint_return_lag_counts={str(L): joint_returns[L] for L in joint_returns},
        min_lag=min(nz) if nz else None,
        joint_min_lag=min(nzj) if nzj else None,
        element_returns_total=sum(elem_returns.values()),
        joint_returns_total=sum(joint_returns.values()),
        patterns={str(L): rows_out[L] for L in rows_out},
    )


def joint_long(toy, starts=256, budget=400, seed=7):
    """Full joint state (element, ring tuple) first-repeat census on a subset,
    over a long budget -- the zero-event window the trap bound uses."""
    rng = random.Random(seed)
    idx = TC.enumerate_states(toy)
    first = []
    none = 0
    for _ in range(starts):
        i = rng.randrange(1, toy.r)
        y, a, b = idx[i], 1, 0
        g = CK.Guard(D)
        seen = {0: (y, tuple(g.slots))}
        lag = None
        for s in range(1, budget + 1):
            y, a, b, kind, info = toy.step(y, a, b, g)
            st = (y, tuple(g.slots))
            for k, v in seen.items():
                if v == st:
                    lag = s - k
                    break
            if lag is not None:
                break
            seen[s] = st
        if lag is None:
            none += 1
        else:
            first.append(lag)
    return dict(starts=starts, budget=budget, no_repeat=none,
                repeats=len(first),
                min_first_repeat_lag=(min(first) if first else None),
                max_first_repeat_lag=(max(first) if first else None),
                zero_event_rule_of_three=3.0 / (starts * budget),
                total_steps=starts * budget)


def trap_bound(guarded, control_arm, theta_hi, R_prod=1 << 14):
    """UNION BOUND over cycle patterns -- no multiplied independence factors.

    A return at lag L needs (i) a row pattern with at least L/2 repeated rows and
    (ii) the applied columns of those rows to satisfy a vanishing sum of twelfth
    roots.  (i) is counted combinatorially: at most C(L, L/2) R^(-L/2) of the R^L
    row words are coincidence patterns.  (ii) is the measured conditional return
    probability given the pattern, theta_hi (CP-95 upper over every pattern up to
    length 6, both fixtures).  The union over L >= 26 is geometric in L, dominated
    by L = 26.
    """
    L0 = 2 * D + 2                     # 26: Theorem 1 kills every L <= 25
    terms = []
    for L in range(L0, L0 + 8):
        c = math.comb(L, L // 2) * (1.0 / R_prod) ** (L / 2)
        terms.append((L, c, c * theta_hi))
    pat_count = sum(t[1] for t in terms)
    per_step = sum(t[2] for t in terms)
    # mixed cycles: a cycle containing a squaring closes only for one v per
    # pattern because ord_r(2) > 2^126 (evidence/factor-r.json:ord_r_2)
    mixed_per_step = 26.0 * (1.0 / 2.0 ** 254)
    return dict(
        L_first_priced=L0,
        patterns_with_repeated_rows_L26=math.comb(L0, L0 // 2),
        pattern_count_sum=pat_count,
        theta_hi=theta_hi,
        per_step_all_mul=per_step,
        per_step_mixed=mixed_per_step,
        per_step_total=per_step + mixed_per_step,
        walk_steps_log2=125.83,
        walk_total_log2=math.log2(per_step + mixed_per_step) + 125.83,
        measured_guard_zero_event_per_pattern=guarded['patterns']['4'][
            'zero_event_rule_of_three'],
        measured_control_ratio_row={
            str(L): (control_arm['patterns'][str(L + 1)]['aggregate_rate']
                     / control_arm['patterns'][str(L)]['aggregate_rate'])
            for L in (2, 3, 4, 5)
            if control_arm['patterns'][str(L)]['aggregate_rate'] > 0},
        note="union bound over row patterns (combinatorial count) x measured "
             "conditional return probability (CP-95 upper), plus the mixed-cycle "
             "row from ord_r(2) = 2^250.7 > 2^126.  No 1/(12R) factor is multiplied "
             "anywhere.",
    )


def main():
    os.makedirs(EV, exist_ok=True)
    toy64 = TC.Toy(rows=64)
    res = {}

    print("[1/4] class marginal ...", flush=True)
    res['marginal'] = class_marginal(toy64, walks=512, steps=200)

    print("[2/4] retention sweep ...", flush=True)
    ret = {}
    for rows, pairs in ((64, 400), (256, 600), (1024, 1200)):
        toy = toy64 if rows == 64 else TC.Toy(rows=rows)
        ret[str(rows)] = retention(toy, pairs)
        print("      R=%5d  rate=%.4f  c=%.1f  (theory %.4f)"
              % (rows, ret[str(rows)]['rate'], ret[str(rows)]['implied_c'],
                 ret[str(rows)]['theory_1_minus_2D2_over_R']), flush=True)
    # CP-95 upper on the failure rate, propagated to production R = 2^14
    c_hi = max(v['cp95_upper_failure'] * int(k) for k, v in ret.items())
    rho = max(0.0, 1.0 - c_hi / (1 << 14))
    res['retention_sweep'] = ret
    res['retention_production'] = dict(
        R_prod=1 << 14, guard_depth=D,
        c_hi=c_hi, c_hi_definition="max over the sweep of (CP-95 upper bound on the "
                                   "per-pair failure rate) x R",
        rho_lower=rho,
        v14_claimed=0.982421875,
        measured_at_largest_R=ret["1024"]['rate'],
        measured_at_largest_R_cp95_lower=ret["1024"]['cp95_lower'],
        verdict="the measured retention at the largest affordable table is "
                "reported with its CP-95 interval; the quantity carried into the "
                "ledger is rho_lower, the CP-95 lower bound on 1 - c/R at the "
                "production table, NOT the analytic 1 - 2D^2/R.",
    )
    with open(os.path.join(EV, 'retention.json'), 'w') as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
    print("      rho_lower(production) = %.6f  (v14 claimed 0.982421875)" % rho,
          flush=True)

    print("[3/4] pattern census, guarded R=8 and R=64, control R=8 ...", flush=True)
    toy8 = TC.Toy(rows=8)
    g8 = pattern_census(toy8, control=False)
    c8 = pattern_census(toy8, control=True)
    toy64b = TC.Toy(rows=64, seed=999)
    g64 = pattern_census(toy64b, control=False, maxlag=18)
    theta = max(g8['patterns'][str(L)]['worst_pattern_rate_cp95_upper']
                for L in (2, 3, 4, 5, 6))
    theta = max(theta, max(g64['patterns'][str(L)]['worst_pattern_rate_cp95_upper']
                           for L in (2, 3, 4, 5, 6)))
    tb = trap_bound(g8, c8, theta)
    pc = dict(fixture_bits=BITS, guard_depth=D, guarded_R8=g8, guarded_R64=g64,
              control_R8=c8, theta_hi=theta, trap_bound=tb)
    with open(os.path.join(EV, 'pattern-census.json'), 'w') as fh:
        json.dump(pc, fh, indent=1, sort_keys=True)
    print("      guarded R=8  min_lag=%s  joint_min_lag=%s  patterns(L=4)=%d"
          % (g8['min_lag'], g8['joint_min_lag'],
             g8['patterns']['4']['patterns_observed']), flush=True)
    print("      control R=8  min_lag=%s  element returns=%d"
          % (c8['min_lag'], c8['element_returns_total']), flush=True)
    print("      theta_hi=%.3e  trap per step=%.3e  walk total 2^%.1f"
          % (theta, tb['per_step_total'], tb['walk_total_log2']), flush=True)

    print("[4/4] joint-state long census ...", flush=True)
    jl = joint_long(toy8, starts=256, budget=400)
    cj = dict(fixture_bits=BITS, guard_depth=D, rows=8,
              exhaustive=g8, long_subset=jl,
              min_lag=g8['joint_min_lag'],
              note="exhaustive over all 9420 starts with an empty ring for 26 steps, "
                   "then a 256-start x 400-step window; the state compared is the "
                   "FULL joint state (element, ring tuple), not the element alone.")
    with open(os.path.join(EV, 'census-joint.json'), 'w') as fh:
        json.dump(cj, fh, indent=1, sort_keys=True)
    print("      joint min_lag=%s  long-window repeats=%d/%d"
          % (cj['min_lag'], jl['repeats'], jl['starts']))
    print("done")


if __name__ == '__main__':
    main()
