#!/usr/bin/env python3
"""v27 measurement layer: the SQUARED-STEP walk, with the judge's three defects fixed.

F1  THE DP PREDICATE IS EXECUTED AT CENSUS TIME.  v26 stored and checked EVERY
    visited class key and called any repeat a success.  Production stores a
    report only when routing.dp_hit(key) fires (density 2^-80) and restarts the
    trail at each report.  Here success is DP-hit AND agree, and the count row is
    kept: every-key store -> DP hits -> reports -> collisions -> agreed.

F1b THE MODEL CARRIES (1-f^2)/mu INSIDE THE EXPONENT: 1-exp(-S^2(1-f^2)/(2 mu N)).
    At S=32, N=785, f=3/8, mu=1.1 that is 0.3992, not 0.4791.

F3  THE RETURN EXPONENT IS PER-j: sigma^j(y)=y^(p^j), so the two-step return is
    y^(4-p^j)=(m1^4 m2^2)^-1, not y^3=sigma^j/(m1^4 m2^2).  e_j=(4-p^j) mod r is
    inverted on the order-r subgroup (gcd on r, NOT on p-1).

F3b full_map_census keys states by the FULL twelve-value class signature, so
    sigma-partners are distinct states; the v26 collapse rate is reported.
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import classkey as CK
import ledger as L
import routing
import tower as TW
from cp95 import cp95_lower

TRAILS = 4
WORKERS = 24
TRIALS = 60000


def step26(toy, y, a, b):
    return CK.step_v26(toy.T, y, a, b, toy.TBL, toy.U, toy.r, toy.tbits, toy.rows)


def make_toy(bits=5, rows=64, seed=20261024):
    import toy_census as TC
    return TC.Toy(bits=bits, seed=seed, rows=rows)


def model_v27(S, N, f, mu):
    """The DISCLOSED production model: 1-exp(-S^2 (1-f^2) / (2 mu N))."""
    lam = (float(S) ** 2) * (1.0 - float(f) ** 2) / (2.0 * float(mu) * float(N))
    return 1.0 - math.exp(-lam), lam


def model_v26_code(S, N):
    """What v26_measure.census actually computed (no factor): the defect."""
    return 1.0 - math.exp(-(float(S) ** 2) / (2.0 * float(N)))


def model_reconciliation(S=32, N=785, f=0.375, mu=1.1):
    m27, lam = model_v27(S, N, f, mu)
    return dict(S=S, N=N, f=f, mu=mu, lambda_v27=lam, model_v27=m27,
                model_v26_code=model_v26_code(S, N), judge_recomputed=0.399,
                delta_vs_v26=model_v26_code(S, N) - m27,
                factor=(1 - f * f) / mu,
                agrees_with_judge=abs(m27 - 0.399) < 5e-4)


def _dp_batch(args):
    bits, seed, trials, dp_bits, ratio = args
    toy = make_toy(bits=bits, seed=seed)
    rng = random.Random(seed ^ 0xB07)
    s_toy = K.step_budget(n=toy.N)
    budget = max(2, int(round(ratio * s_toy)))
    theta = 1 << dp_bits
    cap = K.TRAIL_METERS * theta
    rec = dict(trials=trials, dp_bits=dp_bits, theta=theta, cap=cap, budget=budget,
               steps=0, trials_used=0, trail_starts=0, reports=0, dp_hits=0,
               mul_steps=0, sq_steps=0, stored_all=0, collisions=0, recovered=0,
               zeroden=0, wrong_j=0, invariant_fail=0, first_report_steps=0)
    for _ in range(trials):
        x = rng.randrange(1, toy.r)
        z = toy.T.pow(toy.g, x)
        store = {}
        allkeys = set()
        used = 0
        got = False
        while used < budget:
            a = rng.randrange(1, toy.r)
            b = rng.randrange(1, toy.r)
            y = toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b))
            rec['trail_starts'] += 1
            for _s in range(cap):
                if used >= budget:
                    break
                y, a, b, kind, _i = step26(toy, y, a, b)
                used += 1
                rec['steps'] += 1
                if kind == 'mul':
                    rec['mul_steps'] += 1
                else:
                    rec['sq_steps'] += 1
                k = toy.key(y)[0]
                allkeys.add(k)
                if not routing.dp_hit(k, dp_bits):
                    continue
                rec['dp_hits'] += 1
                rec['reports'] += 1
                if toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b)) != y:
                    rec['invariant_fail'] += 1
                rec['first_report_steps'] += _s + 1
                if k in store:
                    rec['collisions'] += 1
                    a2, b2 = store[k]
                    cands = CK.solve_candidates([pow(toy.p, j, toy.r) for j in range(12)],
                                                None, a, b, a2, b2, toy.r)
                    if not cands:
                        rec['zeroden'] += 1
                        got = False
                    else:
                        got = False
                        for _j, xx in cands:
                            if toy.T.pow(toy.g, xx) == z:
                                got = True
                                break
                    if got:
                        rec['recovered'] += 1
                    else:
                        rec['wrong_j'] += 1
                    break
                store[k] = (a, b)
                break
            if got:
                break
        rec['stored_all'] += len(allkeys)
        rec['trials_used'] += 1
    return rec


def _merge(parts):
    out = dict(parts[0])
    for p in parts[1:]:
        for k, v in p.items():
            if isinstance(v, int) and k not in ('cap', 'theta', 'budget', 'dp_bits'):
                out[k] += v
    return out


def dp_census(trials=TRIALS, bits=5, ratio=1.0, dp_bits=5, workers=WORKERS):
    import multiprocessing as mp
    nw = max(1, min(workers, trials))
    per = [trials // nw] * nw
    for i in range(trials - sum(per)):
        per[i] += 1
    jobs = [(bits, 4242 + 17 * i, n, dp_bits, ratio) for i, n in enumerate(per) if n]
    with mp.Pool(len(jobs)) as pool:
        rec = _merge(pool.map(_dp_batch, jobs))
    n = rec['trials']
    rec['dp_density'] = 2.0 ** (-dp_bits)
    rec['reports_per_trial'] = rec['reports'] / n
    rec['stored_per_trial'] = rec['stored_all'] / n
    rec['p_collide_per_trial'] = rec['collisions'] / n
    rec['p_recover'] = rec['recovered'] / n
    rec['cp95_lower'] = float(cp95_lower(rec['recovered'], n))
    rec['cp95_upper'] = 1.0 - float(cp95_lower(n - rec['recovered'], n))
    rec['p_agree_given_collision'] = (rec['recovered'] / rec['collisions']
                                      if rec['collisions'] else None)
    rec['model_success'] = model_v27(rec['budget'], toyN(bits), K.F_SQR, K.MU_WALK)[0]
    rec['model_success_v26_code'] = model_v26_code(rec['budget'], toyN(bits))
    rec['sqr_fraction'] = (rec['sq_steps'] / rec['steps']) if rec['steps'] else None
    return rec


def toyN(bits):
    import toy_census as TC
    return TC.FIX[bits]['N']


def return_exponent(p, r, j):
    e = (4 - pow(p, j, r)) % r
    return e, math.gcd(e, r)


def recover_y_fast(T, r, m1, m2, j, p):
    """y with ((y m1)^2 m2)^2 = sigma^j(y), i.e. y^(4-p^j) = (m1^4 m2^2)^-1."""
    e, g = return_exponent(p, r, j)
    if g != 1:
        return None, e, g
    prod = T.mul(T.pow(m1, 4), T.pow(m2, 2))
    inv = T.pow(prod, r - 1)
    return T.pow(inv, pow(e, -1, r)), e, g


def recover_y_old(T, r, m1, m2, j, p):
    """The v26 formula: y^3 = sigma^j/(m1^4 m2^2).  Correct only at j=0."""
    prod = T.mul(T.pow(m1, 4), T.pow(m2, 2))
    inv = T.pow(prod, r - 1)
    return T.pow(inv, pow(3, -1, r)), 3, math.gcd(3, r)


def recovery_proof(bits=5, seed=77021, j_list=None, msgs=50):
    toy = make_toy(bits=bits, seed=seed)
    T, r, p = toy.T, toy.r, toy.p
    rng = random.Random(seed ^ 0xABCD)
    j_list = [0, 1, 2, 3] if j_list is None else j_list
    out = dict(samples=0, misses=0, subgroup_fail=0, gcd_fail=0, per_j={},
               j0_regression_ok=0, j0_regression_total=0,
               old_formula_hits_j0=0, old_formula_hits_jne0=0, jne0_total=0,
               miss_examples=[])
    for j in j_list:
        st = dict(samples=0, misses=0, old_hits=0, differ=0)
        for _ in range(200):
            t = rng.randrange(toy.rows)
            kap = rng.randrange(12)
            t2 = rng.randrange(toy.rows)
            kap2 = rng.randrange(12)
            m1 = toy.TBL[t][kap]
            m2 = toy.TBL[t2][kap2]
            for _m in range(msgs):
                y, e, g = recover_y_fast(T, r, m1, m2, j, p)
                out['samples'] += 1
                st['samples'] += 1
                if y is None:
                    out['gcd_fail'] += 1
                    continue
                if T.pow(y, r) != T.one12():
                    out['subgroup_fail'] += 1
                y1 = T.pow(T.mul(y, m1), 2)
                y2 = T.pow(T.mul(y1, m2), 2)
                if y2 != toy.frob_pow(y, j):
                    out['misses'] += 1
                    st['misses'] += 1
                    if len(out['miss_examples']) < 5:
                        out['miss_examples'].append([j, t, kap, t2, kap2])
                yo, _e3, _g3 = recover_y_old(T, r, m1, m2, j, p)
                if j == 0:
                    out['j0_regression_total'] += 1
                    if yo == y:
                        out['j0_regression_ok'] += 1
                    if (T.pow(T.mul(T.pow(T.mul(yo, m1), 2), m2), 2) == yo):
                        out['old_formula_hits_j0'] += 1
                else:
                    out['jne0_total'] += 1
                    if (T.pow(T.mul(T.pow(T.mul(yo, m1), 2), m2), 2)
                            == toy.frob_pow(yo, j)):
                        out['old_formula_hits_jne0'] += 1
                        st['old_hits'] += 1
                if yo != y:
                    st['differ'] += 1
        out['per_j'][str(j)] = st
    out['zero_misses'] = out['misses'] == 0
    out['j0_regression_kept'] = (out['j0_regression_ok'] == out['j0_regression_total']
                                 and out['j0_regression_total'] > 0)
    out['old_formula_wrong_for_j_ne_0'] = out['old_formula_hits_jne0'] == 0
    return out


def full_map_census(toy, which='v26'):
    T, r = toy.T, toy.r
    elems = [None] * r
    y = toy.g
    for i in range(1, r):
        elems[i] = y
        y = T.mul(y, toy.g)
    exact = {}
    classkey = {}
    for i in range(1, r):
        exact[tuple(toy.T.keys12(elems[i]))] = i
        classkey.setdefault(toy.key(elems[i])[0], i)
    nxt = [0] * r
    for i in range(1, r):
        if which == 'v26':
            y2, _a, _b, kind, _info = step26(toy, elems[i], 1, 0)
        else:
            y2, _a, _b, kind, _info = CK.step_v13(T, elems[i], 1, 0, toy.TBL,
                                                  toy.U, r, toy.tbits)
        nxt[i] = exact[tuple(toy.T.keys12(y2))]
    state = [0] * r
    hist = {}
    two_cycles = 0
    for i in range(1, r):
        if state[i]:
            continue
        path, j = [], i
        while j and not state[j]:
            state[j] = 1
            path.append(j)
            j = nxt[j]
        if j and state[j] == 1:
            cyc = path[path.index(j):]
            hist[len(cyc)] = hist.get(len(cyc), 0) + 1
            if len(cyc) == 2:
                two_cycles += 1
        for q in path:
            state[q] = 2
    return dict(rule=which, elements=r - 1, distinct_exact_states=len(exact),
                distinct_class_keys=len(classkey),
                collapse_factor=(r - 1) / max(1, len(classkey)),
                cycles=sum(hist.values()),
                cycle_lengths={str(k): v for k, v in sorted(hist.items())},
                two_cycles=two_cycles,
                predicted_two_cycles_expectation=(r - 1) / (12.0 * len(toy.us)))


def theorem_scan(toy, plant=None):
    T = toy.T
    TBL = [row[:] for row in toy.TBL]
    planted = None
    if plant is not None:
        t0, k0 = plant
        TBL[t0][k0] = T.mul(TBL[t0][k0], toy.tbl[(t0 + 1) % toy.rows])
        planted = (t0, k0)
    n_rows, n_cols = len(TBL), len(TBL[0])
    pairs, pair_fail, orbit_ok, orbit_fail = 0, [], 0, []
    for t in range(n_rows):
        for k in range(6):
            if T.mul(TBL[t][k], TBL[t][k + 6]) == T.one12():
                pairs += 1
            else:
                pair_fail.append([t, k, k + 6])
        for k in range(n_cols):
            if toy.frob_pow(TBL[t][0], (-k) % 12) == TBL[t][k]:
                orbit_ok += 1
            else:
                orbit_fail.append([t, k])
    inv = all(toy.frob_pow(y, 6) == T.pow(y, toy.r - 1)
              for y in toy.tbl[::max(1, toy.rows // 8)])
    return dict(rows=n_rows, cols=n_cols, sigma6_is_inversion=inv,
                inverse_pair_columns=pairs,
                inverse_pair_columns_total=n_rows * 6,
                inverse_pair_failures=pair_fail[:20], orbit_columns=orbit_ok,
                orbit_columns_total=n_rows * n_cols, orbit_failures=orbit_fail[:20],
                plant=planted,
                clean=(pairs == n_rows * 6 and orbit_ok == n_rows * n_cols
                       and not pair_fail and not orbit_fail and inv))


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'all'
    out = {}
    os.makedirs('evidence', exist_ok=True)
    if cmd in ('all', 'model'):
        out['model'] = model_reconciliation()
        print('MODEL v27 %.6f  v26code %.6f  judge 0.399  agrees %s'
              % (out['model']['model_v27'], out['model']['model_v26_code'],
                 out['model']['agrees_with_judge']))
    if cmd in ('all', 'recovery'):
        out['recovery'] = recovery_proof()
        rr = out['recovery']
        print('RECOVERY samples %d misses %d subgroup_fail %d gcd_fail %d j0_regression %s'
              % (rr['samples'], rr['misses'], rr['subgroup_fail'], rr['gcd_fail'],
                 rr['j0_regression_kept']))
        for j, st in rr['per_j'].items():
            print('  j=%s samples %d misses %d old_hits %d differ %d'
                  % (j, st['samples'], st['misses'], st['old_hits'], st['differ']))
    if cmd in ('all', 'cycle'):
        toy = make_toy()
        out['cycle'] = dict(v26=full_map_census(toy, 'v26'), v13=full_map_census(toy, 'v13'))
        c = out['cycle']['v26']
        print('CYCLE v26 elements %d distinct_exact %d distinct_class %d collapse %.2f '
              'two_cycles %d cycles %d'
              % (c['elements'], c['distinct_exact_states'], c['distinct_class_keys'],
                 c['collapse_factor'], c['two_cycles'], c['cycles']))
    if cmd in ('all', 'theorem'):
        toy = make_toy()
        clean = theorem_scan(toy)
        planted = theorem_scan(toy, plant=(7, 3))
        out['theorem'] = dict(clean=clean, planted=planted,
                              planted_flagged=(not planted['clean']))
        print('THEOREM pairs %d/%d orbit %d/%d sigma6 %s planted_flagged %s'
              % (clean['inverse_pair_columns'], clean['inverse_pair_columns_total'],
                 clean['orbit_columns'], clean['orbit_columns_total'],
                 clean['sigma6_is_inversion'], out['theorem']['planted_flagged']))
    with open('evidence/v27-measure-%s.json' % cmd, 'w') as f:
        json.dump(out, f, sort_keys=True, default=str, separators=(',', ':'))
