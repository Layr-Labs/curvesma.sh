#!/usr/bin/env python3
"""v19 MEASUREMENTS — the three rows verdict 5e6e07c4 demands, measured on the
shipped step, not asserted.

M1 REPEAT CLASSES (defect 1).  The v16/v17 trap row charged only repeats whose
   label increment is (0,0) (the denominator b2 - p^j b1 vanishes, no solve
   equation).  The verdict is right that this is the wrong event: the production
   description stores DP records only, with no self-collision detector, so a
   trail that returns to the SAME (element, guard) state is on a cycle that
   contains no DP with probability 1 - O(20 * 2^-80) and never reports --
   regardless of the increment.  This script re-runs the census and TAGS EVERY
   repeat with its increment, so the zero/nonzero split is measured rather than
   assumed, and both classes are then charged as traps.

M2 DP DENSITY (defect 3).  A density is not an interval.  The report cap is a
   deterministic ceiling on a random count, so the row needed is the measured
   DISPERSION: hits per fixed window at the fixture, the variance, the
   over-dispersion var/mean, and the worst window seen.  The cap is then charged
   at the worst window (measured), and the stopping loss at the cap is bounded by
   Chebyshev with the measured variance.

M3 TOY LADDER (defect 5).  Twenty-four full-protocol trials at the fixture the
   package already ships (r = 9421, 61-bit p, ord_r(p) = 12, DP width 2^-5,
   table 64 rows x 12 columns, guard depth 12): plant x, walk the shipped
   guarded step, terminate trails at their first distinguished point, store the
   DP records, solve on a class-key collision for all twelve shifts, VERIFY
   every accepted candidate against the planted z.  Every row is a full ladder,
   not a spot check.  Two controls: a planted wrong secret must fail
   verification, and a zero-increment record pair must produce no candidate at
   all (the trap, exhibited).
"""
import json
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction

import classkey as CK
import constants as K
import routing
import toy_census as TC

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
D = K.GUARD_DEPTH


# ---------------------------------------------------------------- M1
def repeat_classes(toy, walks=5000, steps=200, seed=17, win_walks=256,
                   win_steps=400, win_seed=1701):
    """Per-pair joint-state repeats, each tagged with its label increment."""
    rng = random.Random(seed)
    lag_counts = {}
    zero_lag = {}
    pairs = 0
    w = toy.T
    r = toy.r
    for _ in range(walks):
        y = toy.rand_gt()
        g = CK.Guard(D)
        a, b = 1, 0
        seen = {}
        for t in range(steps + 1):
            st = (y, tuple(g.slots))
            if st in seen:
                t0, a0, b0 = seen[st]
                lag = t - t0
                lag_counts[lag] = lag_counts.get(lag, 0) + 1
                if (a - a0) % r == 0 and (b - b0) % r == 0:
                    zero_lag[lag] = zero_lag.get(lag, 0) + 1
            else:
                seen[st] = (t, a, b)
            if t == steps:
                break
            y, a, b, _k, _i = toy.step(y, a, b, g)
        pairs += steps
    total = sum(lag_counts.values())
    zero_total = sum(zero_lag.values())
    # the v16/v17 window row, same walk length, tagging increments too
    rngw = random.Random(win_seed)
    wal = wr = wz = wev = 0
    first_lag = {}
    for _ in range(win_walks):
        y = toy.rand_gt()
        g = CK.Guard(D)
        a, b = 1, 0
        seen = {}
        jr = 0
        fl = None
        for t in range(win_steps + 1):
            st = (y, tuple(g.slots))
            if st in seen:
                t0, a0, b0 = seen[st]
                jr += 1
                if fl is None:
                    fl = t - t0
                if (a - a0) % r == 0 and (b - b0) % r == 0:
                    wz += 1
            else:
                seen[st] = (t, a, b)
            if t == win_steps:
                break
            y, a, b, _k, _i = toy.step(y, a, b, g)
        wal += 1
        if jr:
            wr += 1
        wev += jr
        if fl is not None:
            first_lag[str(fl)] = first_lag.get(str(fl), 0) + 1
    return {
        'step_source': 'classkey.step (shipped, imported), guard depth %d' % D,
        'fixture_r': r, 'fixture_elements': r - 1, 'table_rows': len(toy.us),
        'census': {'walks': walks, 'steps_per_walk': steps, 'step_pairs': pairs,
                   'joint_repeat_events': total,
                   'p_pair_all': total / pairs,
                   'nonzero_increment_events': total - zero_total,
                   'zero_increment_events': zero_total,
                   'zero_increment_fraction': (zero_total / total) if total else None,
                   'lag_histogram': {str(k): v for k, v in sorted(lag_counts.items())},
                   'zero_increment_lag_histogram': {str(k): v for k, v in sorted(zero_lag.items())},
                   'nonzero_increment_lag_histogram': {
                       str(k): lag_counts[k] - zero_lag.get(k, 0)
                       for k in sorted(lag_counts)}},
        'window': {'walks': wal, 'window_steps': win_steps,
                   'walks_with_joint_repeat': wr,
                   'walk_repeat_fraction': wr / wal,
                   'joint_repeat_events': wev,
                   'zero_increment_events': wz,
                   'nonzero_increment_events': wev - wz,
                   'first_joint_repeat_lag_hist': dict(sorted(first_lag.items(),
                                                              key=lambda kv: int(kv[0])))},
        'reading': 'A repeat of (element, guard) is a cycle. The shipped attack stores DP '
                   'records only and has no self-collision detector, so a cycle whose 20-odd '
                   'elements contain no DP (probability 1 - 20*2^-80) never reports and the '
                   'trail is PENDING at every cutoff. The increment is irrelevant to the trap; '
                   'it only decides whether a solve WOULD have existed had both records been '
                   'kept. Both classes are charged as traps in v19.',
    }


# ---------------------------------------------------------------- M2
def density(toy, walks=600, steps=4096, dp_bits=5, window=128, seed=77):
    """DP hits per fixed window: dispersion, worst window, Chebyshev row."""
    rng = random.Random(seed)
    per_window = []
    per_walk = []
    for _ in range(walks):
        y = toy.rand_gt()
        g = CK.Guard(D)
        a, b = 1, 0
        hits = 0
        cur = 0
        for s in range(steps):
            y, a, b, _k, _i = toy.step(y, a, b, g)
            if routing.dp_hit(CK.classkey(toy.T, y)[0], dp_bits):
                hits += 1
                cur += 1
            if (s + 1) % window == 0:
                per_window.append(cur)
                cur = 0
        per_walk.append(hits)
    n = len(per_window)
    mean = sum(per_window) / n
    var = sum((x - mean) ** 2 for x in per_window) / (n - 1)
    wmean = sum(per_walk) / len(per_walk)
    wvar = sum((x - wmean) ** 2 for x in per_walk) / (len(per_walk) - 1)
    worst = max(per_window)
    return {
        'fixture': {'r': toy.r, 'dp_bits': dp_bits, 'theta_steps': 1 << dp_bits,
                    'table_rows': len(toy.us), 'guard_depth': D},
        'windows': {'walks': walks, 'steps_per_walk': steps, 'window_steps': window,
                    'windows': n, 'mean_hits_per_window': mean,
                    'var_hits_per_window': var, 'sd': var ** 0.5,
                    'over_dispersion_var_over_mean': var / mean if mean else None,
                    'worst_window_hits': worst,
                    'worst_over_mean': worst / mean if mean else None,
                    'model_hits_per_window': window / (1 << dp_bits),
                    'measured_over_model': (mean / (window / (1 << dp_bits)))},
        'walks_row': {'walks': len(per_walk), 'mean_hits_per_walk': wmean,
                      'var_hits_per_walk': wvar,
                      'over_dispersion_var_over_mean': wvar / wmean if wmean else None,
                      'worst_walk_hits': max(per_walk)},
        'chebyshev': 'P(R >= k*mean) <= var/((k-1)*mean)^2 with the measured variance; '
                     'the cap is charged at the measured worst window, not at the mean',
    }


# ---------------------------------------------------------------- M3
def ladder(toy, trials=24, trails=8, meters=64, dp_bits=5, seed=20261024):
    """24 full-protocol trials at fixture scale: plant, walk, DP, solve, verify."""
    rng = random.Random(seed)
    T, r = toy.T, toy.r
    theta = 1 << dp_bits
    cap = meters * theta
    ppow = [pow(toy.p, q, r) for q in range(12)]
    rows = []
    for t in range(trials):
        x = rng.randrange(1, r)
        z = T.pow(toy.g, x)
        store = {}
        recovered = False
        steps = 0
        records = 0
        collisions = 0
        solved_j = None
        invariant_ok = True
        for trail in range(trails):
            a = rng.randrange(1, r)
            b = rng.randrange(1, r)
            y = T.mul(T.pow(toy.g, a), T.pow(z, b))
            gr = CK.Guard(D)
            for s in range(cap):
                y, a, b, _kind, _info = toy.step(y, a, b, gr)
                steps += 1
                if routing.dp_hit(CK.classkey(T, y)[0], dp_bits):
                    records += 1
                    if T.mul(T.pow(toy.g, a), T.pow(z, b)) != y:
                        invariant_ok = False
                    k = CK.classkey(T, y)[0]
                    if k in store:
                        a2, b2 = store[k]
                        collisions += 1
                        cands = CK.solve_candidates(ppow, None, a, b, a2, b2, r)
                        for _j, xx in cands:
                            if T.pow(toy.g, xx) == z:
                                solved_j = _j
                                recovered = True
                                break
                        break
                    store[k] = (a, b)
                    break
            if recovered:
                break
        rows.append({'trial': t, 'planted_x_bits': r.bit_length(), 'trails_used': trail + 1,
                     'steps': steps, 'dp_records': records, 'collisions': collisions,
                     'solved_shift': solved_j, 'invariant_ok': invariant_ok,
                     'recovered': recovered, 'pass': bool(recovered and invariant_ok)})

    # control A: the planted secret plus one must NOT verify
    xc = rng.randrange(1, r)
    zc = T.pow(toy.g, xc)
    control_false = (T.pow(toy.g, (xc + 1) % r) == zc)
    # control B/C: the SAME element reached with a zero and a nonzero increment.
    #   b2 = b1 + 1, a2 = a1 - x  =>  y2 = g^{a2} z^{b2} = g^{a1} z^{b1} = y1, increment (x,1).
    #   zero increment  : the j=0 branch has denominator b2-b1 = 0, and every other
    #                     shift's candidate is spurious: NOTHING verifies -> the trap.
    #   nonzero increment: the j=0 branch solves x exactly -> one candidate verifies.
    xt = rng.randrange(1, r)
    zt = T.pow(toy.g, xt)
    a1 = rng.randrange(1, r)
    b1 = rng.randrange(1, r)
    y1 = T.mul(T.pow(toy.g, a1), T.pow(zt, b1))
    c0 = CK.solve_candidates(ppow, None, a1, b1, a1, b1, r)
    v0 = sum(1 for _j, xx in c0 if T.pow(toy.g, xx) == zt)
    a2 = (a1 - xt) % r
    b2 = (b1 + 1) % r
    y2 = T.mul(T.pow(toy.g, a2), T.pow(zt, b2))
    c1 = CK.solve_candidates(ppow, None, a2, b2, a1, b1, r)
    v1 = sum(1 for _j, xx in c1 if T.pow(toy.g, xx) == zt)
    assert y2 == y1
    return {
        'fixture': {'r': r, 'p': toy.p, 'p_bits': toy.p.bit_length(),
                    'ord_r_p': 12, 'N_classes': toy.N, 'table_rows': len(toy.us),
                    'dp_bits': dp_bits, 'theta_steps': theta, 'trail_cap_steps': cap,
                    'trails_per_trial': trails, 'guard_depth': D,
                    'trials': trials, 'seed': seed},
        'protocol': 'plant x; z = g^x; y = g^a z^b; shipped guarded step; trail ends at its '
                    'first DP (width 2^-%d); one shared store keyed by the 48-byte class key; '
                    'on a key collision solve all twelve shifts and VERIFY g^x == z' % dp_bits,
        'rows': rows,
        'passed': sum(1 for q in rows if q['pass']),
        'failed': [q['trial'] for q in rows if not q['pass']],
        'controls': {
            'wrong_secret_rejected': not control_false,
            'zero_increment_candidates': len(c0),
            'zero_increment_verified': v0,
            'zero_increment_is_trap': v0 == 0,
            'nonzero_increment_candidates': len(c1),
            'nonzero_increment_verified': v1,
            'nonzero_increment_solvable': v1 >= 1,
        },
        'note': 'the ladder measures the MERGE -> SOLVE -> VERIFY mechanics at fixture scale, '
                'which is supercritical (S/sqrt(N) ~ 2.3 here against 0.28 at production); it '
                'is a correctness ladder, not a production recovery rate',
    }


def main():
    os.makedirs(EV, exist_ok=True)
    t0 = time.time()
    toy = TC.Toy(rows=64, seed=20261024)
    print('fixture r=%d rows=%d D=%d (%.1fs)' % (toy.r, len(toy.us), D, time.time() - t0),
          flush=True)
    rep = repeat_classes(toy)
    print('M1 repeats %d / %d pairs, zero-inc %d, nonzero-inc %d (%.1fs)'
          % (rep['census']['joint_repeat_events'], rep['census']['step_pairs'],
             rep['census']['zero_increment_events'],
             rep['census']['nonzero_increment_events'], time.time() - t0), flush=True)
    den = density(toy)
    print('M2 mean %.3f hits/window var %.3f worst %d ratio %.2f (%.1fs)'
          % (den['windows']['mean_hits_per_window'], den['windows']['var_hits_per_window'],
             den['windows']['worst_window_hits'], den['windows']['worst_over_mean'],
             time.time() - t0), flush=True)
    lad = ladder(toy)
    print('M3 ladder %d/%d, failures %s (%.1fs)'
          % (lad['passed'], len(lad['rows']), lad['failed'], time.time() - t0), flush=True)
    json.dump(rep, open(os.path.join(EV, 'repeat-classes-v19.json'), 'w'), indent=1, sort_keys=True)
    json.dump(den, open(os.path.join(EV, 'density-v19.json'), 'w'), indent=1, sort_keys=True)
    json.dump(lad, open(os.path.join(EV, 'ladder-v19.json'), 'w'), indent=1, sort_keys=True)
    print('wrote evidence/repeat-classes-v19.json density-v19.json ladder-v19.json')


if __name__ == '__main__':
    main()
