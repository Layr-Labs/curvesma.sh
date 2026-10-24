#!/usr/bin/env python3
"""v29 measurement layer.  Four judge-named defects measured, not declared.

F1  AGE AT THE MERGE, NOT AT THE REPORT.  v28 stored only distinguished-point
    (report) keys, so the recorded age was conditioned on the follower reaching a
    report key that was already stored.  The retention functional needs the age at
    the MERGE EVENT (the repeat of ANY visited class key), because after that
    repeat the follower's class-key sequence is the elder's, so its next report is
    the elder's report at a deterministic remaining distance D.  Here every visited
    class key is stored, every repeat is a merge event, and for each merge we
    record (a) the follower age A in exact integer steps (1-step histogram, no
    coarse bins), (b) the remaining distance D to the follower's next report, and
    (c) whether that report is credited (solve + verify gate).  The endpoint-only
    age law is kept for comparison.

F2  L^2 SIZE BIAS, EXACT SUMS.  Unchanged from v28: sums of L, L^2, L^3 are exact
    integers and E[L^3]/(2 E[L^2]) is recomputed from them.

F3  STOPPING: UNIT-EXACT RECOVERY STEPS.  Each trial records the exact step at
    which a verified recovery happened; the halt experiments are exact counts of
    trials with rec_step > S_toy - r, never decile-histogram arithmetic.

F4  CONCURRENCY.  The launch schedule is the one the ledger already declares: one
    trail start per theta executed steps (measured in the census as trail_starts
    ~ steps/theta).  A trail lasts at most cap = 128 theta steps, so at most
    ceil(cap/theta) = 128 trails are in flight at any completed-step count.  The
    simulator below checks that bound at the maximum launch rate.
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
CAP_MULT_PROD = 128


def _solve_verify(toy, a, b, a2, b2, z, shifts=12):
    """TWELVE-SHIFT recovery gate (v30).  A class collision is y1 = sigma^j y2
    for one j in Z/12, i.e. a + b x = p^j (a2 + b2 x), NOT the j=0 equation.
    Try every shift, solve x = (p^j a2 - a)/(b - p^j b2), and VERIFY g^x = z.
    v29 tried only j=0 and rejected 11/12 of legitimate recoveries."""
    for j in range(shifts):
        pj = toy.ppow[j]
        den = (b - pj * b2) % toy.r
        if den == 0:
            continue
        x = ((pj * a2 - a) % toy.r) * pow(den, -1, toy.r) % toy.r
        if toy.T.pow(toy.g, x) == z:
            return True, j
    return False, None


def _batch(args):
    bits, seed, trials, dp_bits = args
    toy = TC.Toy(bits=bits, seed=seed)
    rng = random.Random(seed ^ 0xA9E5)
    theta = 1 << dp_bits
    cap = K.TRAIL_METERS * theta            # fixture cap 16*theta
    S_toy = K.step_budget(n=toy.N)
    rec = dict(trials=trials, dp_bits=dp_bits, theta=theta, cap=cap, budget=S_toy,
               n_trails=0, sum_l=0, sum_l2=0, sum_l3=0, max_l=0, censored=0,
               merges_all=0, merges_report=0, verified=0, unverified=0,
               age_hist_all={}, age_hist_report={}, age_over_theta_hist={},
               age_sum=0, age_max=0, elder_sum=0,
               ret_all_sum=0.0, ret_rep_sum=0.0,
               d_hist={}, d_censored=0, report_age_hist={},
               recoveries=0, rec_step_hist={}, rec_frac_hist={},
               shift_hist={}, len_hist={})
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
            pending = []                     # merge ages awaiting this trail's report
            for _s in range(cap):
                if used >= S_toy:
                    end = 'budget'
                    break
                y, a, b, kind, _i = CK.step_v26(toy.T, y, a, b, toy.TBL, toy.U,
                                                toy.r, toy.tbits, toy.rows)
                used += 1
                L += 1
                k = toy.key(y)[0]
                is_dp = routing.dp_hit(k, dp_bits)
                if k in store:
                    a2, b2, L2 = store[k]
                    rec['merges_all'] += 1
                    rec['age_hist_all'][str(L)] = rec['age_hist_all'].get(str(L), 0) + 1
                    rec['age_over_theta_hist'][str(int(L // theta))] = \
                        rec['age_over_theta_hist'].get(str(int(L // theta)), 0) + 1
                    rec['age_sum'] += L
                    rec['age_max'] = max(rec['age_max'], L)
                    rec['elder_sum'] += L2
                    rec['ret_all_sum'] += math.exp(-(CAP_MULT_PROD - L / float(theta)))
                    ok, sj = _solve_verify(toy, a, b, a2, b2, z)
                    if ok:
                        rec['verified'] += 1
                        rec['shift_hist'][str(sj)] = rec['shift_hist'].get(str(sj), 0) + 1
                    else:
                        rec['unverified'] += 1
                    if not is_dp:
                        pending.append(L)     # internal merge; trail continues
                        continue
                    rec['merges_report'] += 1
                    rec['age_hist_report'][str(L)] = \
                        rec['age_hist_report'].get(str(L), 0) + 1
                    rec['ret_rep_sum'] += math.exp(-(CAP_MULT_PROD - L / float(theta)))
                    for A in pending:
                        rec['d_hist'][str(L - A)] = rec['d_hist'].get(str(L - A), 0) + 1
                    pending = []
                    if ok and not got:
                        got = True
                        rec['recoveries'] += 1
                        rec['rec_step_hist'][str(used)] = \
                            rec['rec_step_hist'].get(str(used), 0) + 1
                        b10 = min(9, int(10 * used / float(S_toy)))
                        rec['rec_frac_hist'][str(b10)] = \
                            rec['rec_frac_hist'].get(str(b10), 0) + 1
                    end = 'report'
                    break
                store[k] = (a, b, L)
                if is_dp:
                    # this trail reports here, at age L; every pending merge's
                    # remaining distance D = L - A
                    rec['report_age_hist'][str(L)] = \
                        rec['report_age_hist'].get(str(L), 0) + 1
                    for A in pending:
                        rec['d_hist'][str(L - A)] = rec['d_hist'].get(str(L - A), 0) + 1
                    pending = []
                    end = 'report'
                    break
            # merges still pending when the trail is censored: D > cap - A
            rec['d_censored'] += len(pending)
            rec['n_trails'] += 1
            rec['sum_l'] += L
            rec['sum_l2'] += L * L
            rec['sum_l3'] += L ** 3
            rec['len_hist'][str(L)] = rec['len_hist'].get(str(L), 0) + 1
            rec['max_l'] = max(rec['max_l'], L)
            if end == 'budget':
                rec['censored'] += 1
            if got:
                break
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
            elif isinstance(v, int) and k in ('age_max', 'max_l'):
                out[k] = max(out[k], v)
            elif isinstance(v, int) and k not in ('dp_bits', 'theta', 'cap', 'budget'):
                out[k] += v
    return out


def moments(rec):
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


def halt_row(rec, r):
    """Unit-exact: trials whose verified recovery step exceeds S_toy - r."""
    s = rec['budget']
    halt = s - r
    lost = sum(v for k, v in rec['rec_step_hist'].items() if int(k) > halt)
    return dict(halt_step=halt, lost_trials=lost,
                lost_frac=lost / float(rec['trials']),
                recoveries=rec['recoveries'], trials=rec['trials'])


def max_inflight(lengths, interval, horizon=None):
    """Launch one trail every `interval` steps; return the max concurrent count."""
    if not lengths:
        return 0
    horizon = horizon or (len(lengths) * interval + max(lengths) + 1)
    launches = list(range(0, horizon, interval))
    events = []
    for i, t in enumerate(launches):
        ln = lengths[i % len(lengths)]
        events.append((t, 1))
        events.append((t + ln, -1))
    events.sort(key=lambda e: (e[0], e[1]))     # finishes before launches at a tie
    cur = best = 0
    for _t, d in events:
        cur += d
        best = max(best, cur)
    return best


def scheduler_trace(launch_period, trail_len, rounds, enforce=None):
    """Discrete synchronous scheduler on the DECLARED clock: a launch is allowed
    when (completed aggregate steps since the last launch) >= launch_period; each
    in-flight trail executes one own step per round, and a trail is capped at
    trail_len own steps.  `enforce` is an atomic in-flight cap (None = the bare
    schedule).  Returns the observed maximum in-flight count."""
    inflight = []
    done = 0
    last = 0
    best = 0
    launched = 0
    for _r in range(rounds):
        first = (launched == 0 and not inflight)
        if (first or done - last >= launch_period) and (enforce is None or len(inflight) < enforce):
            inflight.append(trail_len)
            last = done
            launched += 1
        done += len(inflight)              # aggregate completed steps this round
        inflight = [x - 1 for x in inflight if x - 1 > 0]
        best = max(best, len(inflight))
    return dict(max_inflight=best, launched=launched, completed_steps=done,
                rounds=rounds, launch_period=launch_period, trail_len=trail_len,
                enforced_cap=enforce)


def concurrency_control(theta, cap_mult=CAP_MULT_PROD):
    """v30: the bare schedule is simulated, not assumed.  ceil(cap/theta) is NOT
    a bound of the bare schedule (the launch clock is aggregate, the cap is own
    steps: the measured max exceeds it), so v30 enforces an atomic in-flight cap
    and ships both trace maxima."""
    cap = cap_mult * theta
    derived = -(-cap // theta)
    scaled_theta, scaled_cap = 8, 8 * cap_mult        # same ratio, tractable clock
    uncapped = scheduler_trace(scaled_theta, scaled_cap, rounds=4 * scaled_cap)
    capped = scheduler_trace(scaled_theta, scaled_cap, rounds=4 * scaled_cap,
                             enforce=derived)
    toosmall = scheduler_trace(scaled_theta, scaled_cap, rounds=4 * scaled_cap,
                               enforce=derived - 1)
    return dict(launch_interval_steps=theta, trail_cap_steps=cap,
                derived_ratio=derived, launch_cadence_clock='aggregate completed steps',
                trail_cap_clock='own executed steps',
                bare_schedule_max_inflight=uncapped['max_inflight'],
                bare_schedule_bound_attained=(uncapped['max_inflight'] == derived),
                enforced_cap=derived,
                enforced_trace_max_inflight=capped['max_inflight'],
                enforce_holds=(capped['max_inflight'] <= derived),
                toosmall_cap_trace_max_inflight=toosmall['max_inflight'],
                toosmall_bound_fails=(toosmall['max_inflight'] >= derived - 1))


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
        m['age_mean'] = m['age_sum'] / m['merges_all'] if m['merges_all'] else None
        m['age_mean_over_theta'] = (m['age_mean'] / m['theta']) if m['merges_all'] else None
        m['elder_mean_over_theta'] = (m['elder_sum'] / m['merges_all'] / m['theta']
                                      if m['merges_all'] else None)
        m['delta_residual_all'] = (m['ret_all_sum'] / m['merges_all']
                                   if m['merges_all'] else None)
        m['delta_residual_report'] = (m['ret_rep_sum'] / m['merges_report']
                                      if m['merges_report'] else None)
        m['delta_direct_miss'] = (m['d_censored'] / m['merges_all']
                                  if m['merges_all'] else None)
        m['p_recover'] = m['recoveries'] / m['trials']
        m['moments'] = moments(m)
        m['halts'] = {str(r): halt_row(m, r) for r in (1, 2, 4, 8, 16)}
        m['concurrency'] = concurrency_control(m['theta'])
        out[str(dp_bits)] = m
    return out


if __name__ == '__main__':
    os.makedirs('evidence', exist_ok=True)
    tr = run_traces()
    age = dict(rule='age A at the MERGE EVENT (repeat of any visited class key), '
                    'in exact integer steps since the follower trail start; '
                    'D = remaining distance from the merge to the follower report',
               fixture_cap_mult=K.TRAIL_METERS, production_cap_mult=CAP_MULT_PROD,
               uniform_rejected=(1.0 - math.exp(-CAP_MULT_PROD)) / CAP_MULT_PROD,
               arms={})
    mom = dict(rule='E[L^3]/(2 E[L^2]) from exact integer sums, per arm', arms={})
    stop = dict(rule='halt launches when completed steps >= S - reserve; '
                     'reserve = N_inflight * cap, N_inflight = ceil(cap/theta) = 128',
                arms={})
    for k, m in tr.items():
        age['arms'][k] = dict(
            theta=m['theta'], cap=m['cap'], budget=m['budget'], trials=m['trials'],
            n_trails=m['n_trails'], merges_all=m['merges_all'],
            merges_report=m['merges_report'], verified=m['verified'],
            unverified=m['unverified'],
            age_mean_over_theta=m['age_mean_over_theta'],
            elder_mean_over_theta=m['elder_mean_over_theta'],
            age_max=m['age_max'], age_max_over_theta=m['age_max'] / m['theta'],
            age_hist_all=m['age_hist_all'], age_hist_report=m['age_hist_report'],
            age_over_theta_hist=m['age_over_theta_hist'],
            delta_residual_all=m['delta_residual_all'],
            delta_residual_all_log2=(math.log2(m['delta_residual_all'])
                                     if m['delta_residual_all'] else None),
            delta_residual_report=m['delta_residual_report'],
            delta_residual_report_log2=(math.log2(m['delta_residual_report'])
                                        if m['delta_residual_report'] else None),
            delta_direct_miss=m['delta_direct_miss'], d_censored=m['d_censored'],
            d_hist=m['d_hist'], report_age_hist=m['report_age_hist'],
            survival_gt_k_theta={str(kk): sum(v for kk2, v in m['age_hist_all'].items()
                                              if int(kk2) > kk * m['theta'])
                                 for kk in range(1, 17)},
            p_recover=m['p_recover'],
            gate_p_per_merge=(m['verified'] / m['merges_all'] if m['merges_all'] else None),
            shift_hist=m['shift_hist'],
            rec_step_hist=m['rec_step_hist'],
            rec_frac_hist=m['rec_frac_hist'], len_hist=m['len_hist'],
            censored=m['censored'])
        mom['arms'][k] = m['moments']
        rec = dict(m['halts'])
        stop['arms'][k] = dict(S_toy=m['budget'], trials=m['trials'],
                               recoveries=m['recoveries'],
                               rec_step_hist=m['rec_step_hist'],
                               rec_frac_hist=m['rec_frac_hist'],
                               halts=rec, concurrency=m['concurrency'])
    with open('evidence/v30-age-law.json', 'w') as f:
        json.dump(age, f, indent=1, sort_keys=True, default=str)
    with open('evidence/v30-moments.json', 'w') as f:
        json.dump(mom, f, indent=1, sort_keys=True, default=str)
    with open('evidence/v30-stopping.json', 'w') as f:
        json.dump(stop, f, indent=1, sort_keys=True, default=str)
    print('AGE LAW uniform(rejected) %.9f' % age['uniform_rejected'])
    for k, a in sorted(age['arms'].items()):
        print('  dp=%s theta=%d merges_all=%d merges_rep=%d verified=%d age/theta %.3f '
              'delta_all 2^%s delta_rep 2^%s miss %s'
              % (k, a['theta'], a['merges_all'], a['merges_report'], a['verified'],
                 a['age_mean_over_theta'] or -1,
                 ('%.2f' % a['delta_residual_all_log2'])
                 if a['delta_residual_all_log2'] is not None else 'na',
                 ('%.2f' % a['delta_residual_report_log2'])
                 if a['delta_residual_report_log2'] is not None else 'na',
                 a['delta_direct_miss']))
    for k, m in sorted(mom['arms'].items()):
        print('  moments dp=%s n=%d SL=%s SL2=%s SL3=%s ratio/theta=%.4f'
              % (k, m['n_trails'], m['sum_l'], m['sum_l2'], m['sum_l3'],
                 m['ratio_over_theta']))
    for k, s in sorted(stop['arms'].items()):
        print('  stopping dp=%s halts=%s bare=%d enforced=%d toosmall=%d'
              % (k, {r: v['lost_frac'] for r, v in s['halts'].items()},
                 s['concurrency']['bare_schedule_max_inflight'],
                 s['concurrency']['enforced_trace_max_inflight'],
                 s['concurrency']['toosmall_cap_trace_max_inflight']))
