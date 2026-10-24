#!/usr/bin/env python3
"""v28 measurement layer: the judge's three flaws measured, not declared.

F1  MERGE-AGE LAW.  v27 charged P(miss|merge) = (1-e^-128)/128 from a UNIFORM-age
    premise over the 128*theta cap.  First-DP restarts make the renewal age
    nonuniform (a restart biases young ages).  Here every trial records, for each
    trail START, its executed length L (steps until its first report, or until the
    budget/cap), and for every MERGE the follower's age A (steps since its own
    restart) plus the stored elder's age.  From those we emit
      - the age histogram in theta units, conditional on a merge,
      - the empirical survival P(A > k*theta) for k = 1..cap/theta,
      - the retention integral computed on the EMPIRICAL age law, evaluated at the
        PRODUCTION cap 128*theta:  delta_meas = mean_i exp(-(128 - A_i/theta)).
    No uniform assumption anywhere; the uniform value is computed only as the
    DISCLOSED-REJECTED method (retention_uniform).

F2  LENGTH MOMENTS.  Conditioning on a collision weights lengths by L^2, so a
    uniform prefix of the conditioned trail has mean E[L^3]/(2 E[L^2]) -- not
    E[L^2]/(2 E[L]).  The sums of L, L^2, L^3 are accumulated as EXACT INTEGERS and
    the ratio is recomputed from them (moments()), so a perturbation of one sample
    moves the ratio.

F3  STOPPING-RULE SIM.  Each trial records the cumulative step at which it
    recovered (or None).  The halt-point distribution over budget fractions is
    emitted so the success loss of "stop launching at W - reserve" can be read off
    the measured CDF, and the production value is bounded by the S-truncation
    arithmetic (2*reserve/S).
"""
import json
import math
import os
import random
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import classkey as CK
import routing
import toy_census as TC

WORKERS = 12
PLAN = {1: 20000, 2: 20000, 3: 12000, 4: 12000}


def _batch(args):
    bits, seed, trials, dp_bits = args
    toy = TC.Toy(bits=bits, seed=seed)
    rng = random.Random(seed ^ 0xA9E5)
    theta = 1 << dp_bits
    cap = K.TRAIL_METERS * theta            # fixture cap 16*theta (budget is binding)
    S_toy = K.step_budget(n=toy.N)
    rec = dict(trials=trials, dp_bits=dp_bits, theta=theta, cap=cap, budget=S_toy,
               n_trails=0, sum_l=0, sum_l2=0, sum_l3=0, max_l=0, censored=0,
               merges=0, age_hist={}, age_sum=0, age_sum2=0, age_max=0,
               elder_sum=0, elder_max=0,
               ret_meas_sum=0.0, ret_meas_max=0.0, ret_cap16_sum=0.0,
               recoveries=0, rec_step_sum=0, rec_frac_hist={}, lost_at=0,
               age_over_theta_hist={})
    one = toy.T.one12()
    for _ in range(trials):
        x = rng.randrange(1, toy.r)
        z = toy.T.pow(toy.g, x)
        store = {}
        used = 0
        got = False
        while used < S_toy:
            a = rng.randrange(1, toy.r)
            b = rng.randrange(1, toy.r)
            y = toy.T.mul(toy.T.pow(toy.g, a), toy.T.pow(z, b))
            L = 0
            end = 'budget'
            for _s in range(cap):
                if used >= S_toy:
                    end = 'budget'
                    break
                y, a, b, kind, _i = CK.step_v26(toy.T, y, a, b, toy.TBL, toy.U,
                                                toy.r, toy.tbits, toy.rows)
                used += 1
                L += 1
                k = toy.key(y)[0]
                if not routing.dp_hit(k, dp_bits):
                    continue
                end = 'report'
                if k in store:
                    a2, b2, L2 = store[k]
                    rec['merges'] += 1
                    rec['age_hist'][str(L)] = rec['age_hist'].get(str(L), 0) + 1
                    rec['age_over_theta_hist'][str(int(L // theta))] = \
                        rec['age_over_theta_hist'].get(str(int(L // theta)), 0) + 1
                    rec['age_sum'] += L
                    rec['age_sum2'] += L * L
                    rec['age_max'] = max(rec['age_max'], L)
                    rec['elder_sum'] += L2
                    rec['elder_max'] = max(rec['elder_max'], L2)
                    r16 = math.exp(-(cap - L) / theta)
                    rec['ret_cap16_sum'] += r16
                    rec['ret_meas_max'] = max(rec['ret_meas_max'], r16)
                    d = math.exp(-(128.0 - L / float(theta)))
                    rec['ret_meas_sum'] += d
                    got = True
                    break
                store[k] = (a, b, L)
                break
            rec['n_trails'] += 1
            rec['sum_l'] += L
            rec['sum_l2'] += L * L
            rec['sum_l3'] += L ** 3
            rec['max_l'] = max(rec['max_l'], L)
            if end == 'budget':
                rec['censored'] += 1
            if got:
                break
        if got:
            rec['recoveries'] += 1
            rec['rec_step_sum'] += used
            b = min(9, int(10 * used / float(S_toy)))
            rec['rec_frac_hist'][str(b)] = rec['rec_frac_hist'].get(str(b), 0) + 1
    return rec


def _merge(parts):
    out = dict(parts[0])
    for p in parts[1:]:
        for k, v in p.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    out[k][kk] = out[k].get(kk, 0) + vv
            elif isinstance(v, float):
                out[k] += v
            elif isinstance(v, int) and k in ('age_max', 'max_l', 'elder_max'):
                out[k] = max(out[k], v)
            elif isinstance(v, int) and k not in ('dp_bits', 'theta', 'cap', 'budget',
                                                  'trials'):
                out[k] += v
            elif k in ('trials',):
                out[k] += v
    return out


def run_traces(plan=None, bits=5, workers=WORKERS):
    import multiprocessing as mp
    plan = PLAN if plan is None else plan
    jobs = []
    for dp_bits, trials in sorted(plan.items()):
        nw = max(1, min(workers, trials))
        per = [trials // nw] * nw
        for i in range(trials - sum(per)):
            per[i] += 1
        for i, n in enumerate(per):
            if n:
                jobs.append((bits, 9001 + 37 * dp_bits + 11 * i, n, dp_bits))
    with mp.Pool(min(len(jobs), workers)) as pool:
        recs = pool.map(_batch, jobs)
    out = {}
    for dp_bits in plan:
        part = [r for r in recs if r['dp_bits'] == dp_bits]
        m = _merge(part)
        m['age_mean'] = m['age_sum'] / m['merges'] if m['merges'] else None
        m['age_mean_over_theta'] = (m['age_mean'] / m['theta']) if m['merges'] else None
        m['elder_mean_over_theta'] = (m['elder_sum'] / m['merges'] / m['theta']
                                      if m['merges'] else None)
        m['delta_meas_prodcap'] = (m['ret_meas_sum'] / m['merges']
                                   if m['merges'] else None)
        m['delta_meas_prodcap_log2'] = (math.log2(m['delta_meas_prodcap'])
                                        if m['merges'] and m['delta_meas_prodcap'] > 0
                                        else None)
        m['delta_meas_fixturecap'] = (m['ret_cap16_sum'] / m['merges']
                                      if m['merges'] else None)
        m['p_recover'] = m['recoveries'] / m['trials']
        m['censored_frac'] = m['censored'] / m['n_trails'] if m['n_trails'] else None
        m['moments'] = moments(m)
        out[str(dp_bits)] = m
    return out


def moments(rec):
    """Exact moment sums -> conditional wasted-work prefix mean.

    Collision conditioning weights a trail of length L by L^2.  The mean distance
    from the collision point back to a UNIFORM prefix position of that trail is
    L/2, so averaged over the L^2-weighted law the prefix mean is
        E[L^3] / (2 E[L^2]),
    which for an untruncated exponential law is 3 theta / 2 (NOT <= theta).
    """
    n = rec['n_trails']
    s1, s2, s3 = rec['sum_l'], rec['sum_l2'], rec['sum_l3']
    r3 = Fraction(s3, 2 * s2) if s2 else None
    return dict(n_trails=n, sum_l=s1, sum_l2=s2, sum_l3=s3,
                e_l=Fraction(s1, n) if n else None,
                e_l2=Fraction(s2, n) if n else None,
                e_l3=Fraction(s3, n) if n else None,
                ratio_e3_over_2e2=r3,
                ratio_over_theta=(float(r3) / rec['theta']) if r3 else None,
                max_l=rec['max_l'], censored=rec['censored'],
                censored_frac=(rec['censored'] / n) if n else None,
                sum_str=[str(s1), str(s2), str(s3)])


def retention_uniform(cap_mult=128):
    """The v27 method, now DISCLOSED-REJECTED: age uniform on [0, cap]."""
    return (1.0 - math.exp(-cap_mult)) / cap_mult


def retention_analytic(cap_mult=128):
    """Same age law as the trail-length law already declared: A ~ size-biased
    min(Exp(theta), cap).  Density of L is (1/theta)e^{-L/theta} on [0,cap), and
    conditioning on a collision weights by L^2, so
        E[e^{A/theta}] = int_0^cap (a/theta^2) da = cap^2/(2 theta^2),
        delta = e^{-cap/theta} * cap^2/(2 theta^2) = e^-128 * 8192 = 2^-171.66."""
    return math.exp(-cap_mult) * (cap_mult ** 2) / 2.0


def stopping_model(S, cap_mult=128, theta_bits=K.THETA_BITS, c_stress=None):
    """Stopping rule arithmetic: launches stop at W - reserve, each in-flight
    trail charges at most reserve, so total <= W.  Success loss = the S^2 share
    removed by the reserve."""
    theta = 1 << int(theta_bits)
    reserve_steps = cap_mult * theta
    return dict(reserve_steps=str(reserve_steps),
                reserve_frac_of_S=reserve_steps / float(S),
                loss_rel_S2=2.0 * reserve_steps / float(S),
                loss_bits=-math.log2(1.0 - 2.0 * reserve_steps / float(S))
                if 2.0 * reserve_steps < S else None)


if __name__ == '__main__':
    os.makedirs('evidence', exist_ok=True)
    tr = run_traces()
    age_law = dict(rule='age at merge = steps since the FOLLOWER trail start, '
                        'conditional on the merge event',
                   fixture_cap_mult=K.TRAIL_METERS,
                   production_cap_mult=128,
                   uniform_rejected=retention_uniform(),
                   arithmetic_under_declared_law=retention_analytic(),
                   arithmetic_under_declared_law_log2=math.log2(retention_analytic()),
                   arms={})
    mom = dict(rule='E[L^3]/(2 E[L^2]) from exact integer sums, per arm',
               arms={})
    stop = dict(rule='halt launches when cumulative charged work >= W - reserve; '
                     'reserve = one full trail charge = 128*theta*c_mul',
                arms={})
    for k, m in tr.items():
        age_law['arms'][k] = dict(
            theta=m['theta'], cap=m['cap'], budget=m['budget'], trials=m['trials'],
            merges=m['merges'], age_mean_over_theta=m['age_mean_over_theta'],
            elder_mean_over_theta=m['elder_mean_over_theta'],
            age_max=m['age_max'], age_max_over_theta=m['age_max'] / m['theta'],
            age_hist=m['age_hist'], age_over_theta_hist=m['age_over_theta_hist'],
            delta_measured_production_cap=m['delta_meas_prodcap'],
            delta_measured_production_cap_log2=m['delta_meas_prodcap_log2'],
            delta_measured_fixture_cap=m['delta_meas_fixturecap'],
            survival_gt_k_theta={str(kk): sum(v for kk2, v in m['age_hist'].items()
                                              if int(kk2) > kk * m['theta'])
                                 for kk in range(1, 17)},
            p_recover=m['p_recover'])
        mom['arms'][k] = m['moments']
        s_toy = m['budget']
        stop['arms'][k] = dict(
            S_toy=s_toy, trials=m['trials'], recoveries=m['recoveries'],
            rec_frac_hist=m['rec_frac_hist'],
            rec_step_mean_over_S=(m['rec_step_sum'] / m['recoveries'] / s_toy
                                  if m['recoveries'] else None),
            halt_experiments={str(r): dict(
                halt_step=s_toy - r,
                lost_trials=sum(v for kk, v in m['rec_frac_hist'].items()
                                if int(kk) >= min(9, int(10 * (s_toy - r) / s_toy))),
                lost_frac=sum(v for kk, v in m['rec_frac_hist'].items()
                              if int(kk) >= min(9, int(10 * (s_toy - r) / s_toy)))
                / m['trials'])
                for r in (1, 2, 4, 8)})
    with open('evidence/v28-age-law.json', 'w') as f:
        json.dump(age_law, f, indent=1, sort_keys=True, default=str)
    with open('evidence/v28-moments.json', 'w') as f:
        json.dump(mom, f, indent=1, sort_keys=True, default=str)
    with open('evidence/v28-stopping.json', 'w') as f:
        json.dump(stop, f, indent=1, sort_keys=True, default=str)
    print('AGE LAW uniform(rejected) %.9f   analytic(declared law) 2^%.2f'
          % (retention_uniform(), math.log2(retention_analytic())))
    for k, a in sorted(age_law['arms'].items()):
        print('  dp=%s theta=%d merges=%d age/theta mean %.3f max %d  delta_meas 2^%s'
              % (k, a['theta'], a['merges'], a['age_mean_over_theta'] or -1,
                 a['age_max'], ('%.2f' % a['delta_measured_production_cap_log2'])
                 if a['delta_measured_production_cap_log2'] is not None else 'na'))
    for k, m in sorted(mom['arms'].items()):
        print('  moments dp=%s n=%d SL=%s SL2=%s SL3=%s  E[L3]/2E[L2]/theta=%.4f'
              % (k, m['n_trails'], m['sum_l'], m['sum_l2'], m['sum_l3'],
                 m['ratio_over_theta']))
    for k, s in sorted(stop['arms'].items()):
        print('  stopping dp=%s rec_mean/S=%.3f frac_hist=%s'
              % (k, s['rec_step_mean_over_S'] or -1, s['rec_frac_hist']))
