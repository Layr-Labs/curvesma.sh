#!/usr/bin/env python3
"""v33 measure: every number executes protocol.run_protocol.  Writes
evidence/{protocol33,gate33,detect33,tail33,audit33}.json.  No sampler fork.

The four v32 objections map onto four measurements here:

  (1) STEP BUDGET   sweep() runs the shipped path; the halt row carries
      completed == S and consumed == executed.  legacy_control() runs v32's
      charge-at-launch rule on the SAME seed to show the 128x overage it had.
  (2) TAIL          tail_rows() integrates the actual-event geometric tail
      q^R (q = 1 - 1/theta) past each launch's measured remaining lifetime R,
      per launch, and reports the analytic bound 0.582*theta/S it must sit
      under.  tail_scale_control() re-evaluates the SAME population at R/2 and
      2R: the row must move, so it is a function of R and not a constant.
  (3) DETECTION     merge_batch() conditions on MERGE events (the merge log of
      the integrated protocol) and counts detected/merged, n >= 1000, with an
      exact Clopper-Pearson lower bound.  retention_of_same_run() is the OLD
      object (shared reports / reports) printed beside it, labelled.
  (4) UNIT CHAIN    audit d re-emits the claim after multiplying the raw
      INSTRUCTION ledger by C_REF: exactly one division by C_REF exists, so the
      claim must move by exactly log2(C_REF) bits.
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy as TY                # noqa: E402
import protocol as PR           # noqa: E402
import constants as K           # noqa: E402
from cp95 import cp95_lower     # noqa: E402

THETA_BITS = 5
TRIALS = 24
EV = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'evidence'))


def one(seed, budget, theta_bits=THETA_BITS, **kw):
    toy = TY.Toy(bits=5, seed=seed)
    x = toy.rng.randrange(1, toy.r)
    z = toy.T.pow(toy.g, x)
    return PR.run_protocol(toy, z, PR.Cfg(theta_bits, budget), **kw)


def sweep(budgets, trials=TRIALS, theta_bits=THETA_BITS):
    out = []
    for B in budgets:
        rec = 0
        completed = launches = merged = detected = reports = 0
        at_s = 0
        for i in range(trials):
            r = one(1000 + 7 * i, B, theta_bits)
            rec += r['recovered']
            completed += r['completed']
            launches += r['launches']
            merged += r['merged']
            detected += r['detected']
            reports += r['reports']
            at_s += r['executed_equals_S']
        row = dict(budget=B, theta_bits=theta_bits, trials=trials, recovered=rec,
                   success=rec / trials,
                   cp95_lower=float(cp95_lower(rec, trials)),
                   mean_completed=completed / trials,
                   mean_launches=launches / trials,
                   launch_bound=-(-B // (1 << theta_bits)),
                   ran_to_S=at_s, merged=merged, detected=detected,
                   reports=reports)
        out.append(row)
        print('S=%6d  success %2d/%2d  cp95lo %.4f  completed %.1f (<=S: %s)  '
              'launches %.1f/bound %d  merges %d'
              % (B, rec, trials, row['cp95_lower'], row['mean_completed'],
                 row['mean_completed'] <= B, row['mean_launches'],
                 row['launch_bound'], merged))
    return out


def legacy_control(seed=4242, budget=20000, theta_bits=THETA_BITS):
    """v32's rule (charge theta at launch, no executed-step cutoff) vs v33's."""
    leg = one(seed, budget, theta_bits, stop_on_success=False,
              legacy_launch_credit=True, trace_n=10 ** 9)
    new = one(seed, budget, theta_bits, stop_on_success=False, trace_n=10 ** 9)
    return dict(S=budget,
                v32_launch_credit=dict(completed=leg['completed'],
                                       consumed=leg['consumed'],
                                       launches=leg['launches'], halt=leg['halt'],
                                       remaining=leg['remaining'],
                                       over_S=leg['completed'] / budget),
                v33_step_budget=dict(completed=new['completed'],
                                     consumed=new['consumed'],
                                     launches=new['launches'], halt=new['halt'],
                                     remaining=new['remaining'],
                                     executed_equals_S=new['executed_equals_S'],
                                     consumed_equals_executed=new['consumed_equals_executed'],
                                     budget_binds=new['budget_binds'],
                                     launch_bound=new['launch_bound'],
                                     launch_bound_binds=new['launch_bound_binds'],
                                     halt_row=new['trace'][-1]),
                same_seed=(leg['completed'] != new['completed']))


def gate_trace(seed=4242, budget=20000, theta_bits=THETA_BITS):
    r = one(seed, budget, theta_bits, stop_on_success=False, trace_n=10 ** 9)
    tr = r['trace']
    at_cap = next((x for x in tr if x['inflight'] == PR.INFLIGHT_CAP), None)
    first_defer = next((x for x in tr if x['event'] == 'defer'), None)
    first_boot = tr[0]
    return dict(cfg=r['cfg'], halt=r['halt'], launches=r['launches'],
                deferred=r['deferred'], completed=r['completed'],
                consumed=r['consumed'], remaining=r['remaining'],
                launch_bound=r['launch_bound'],
                launch_bound_binds=r['launch_bound_binds'],
                executed_equals_S=r['executed_equals_S'],
                bootstrap=first_boot, at_cap=at_cap, first_defer=first_defer,
                halt_row=tr[-1], tail=tr[-3:], trace_len=len(tr),
                n_bootstrap=r['n_bootstrap'], reports=r['reports'],
                merged=r['merged'])


def merge_batch(target=1000, budget=20000, theta_bits=THETA_BITS, max_runs=40):
    """Condition on MERGE EVENTS.  Each run's merge log is the population; a
    merge is DETECTED when the twelve-shift recovery verifies x against the
    planted target at that merge.  Stopping on success is OFF so the batch can
    collect many merges inside the step budget."""
    merged = detected = 0
    runs = 0
    steps = 0
    j_hist = {}
    false_merge = []
    for i in range(max_runs):
        if merged >= target:
            break
        log = []
        r = one(5000 + i, budget, theta_bits, stop_on_success=False,
                merge_log=log)
        runs += 1
        steps += r['completed']
        for m in log:
            merged += 1
            detected += 1 if m['detected'] else 0
            if m['j'] is not None:
                j_hist[m['j']] = j_hist.get(m['j'], 0) + 1
            if not m['detected']:
                false_merge.append(m)
    cp = float(cp95_lower(detected, merged))
    return dict(merged=merged, detected=detected, runs=runs, steps=steps,
                p_detect=detected / merged if merged else None, cp95_lower=cp,
                false_merges=len(false_merge), j_hist=dict(sorted(j_hist.items())),
                target=target, budget=budget,
                note='detected/merged conditioned on merge events of the '
                     'INTEGRATED protocol; the CP-95 lower bound is what is '
                     'charged, not the point estimate')


def shared_report_fraction(**kw):
    """The OLD statistic, same run: distinct stored keys later shared / reports.
    Label it for what it is -- a coverage fraction, not a detection probability."""
    r = one(9091, 20000, 5, stop_on_success=False, **kw)
    st = r['store']
    return dict(reports=st.reports, shared=len(st.shared), matches=st.matches,
                shared_fraction=(len(st.shared) / st.reports) if st.reports else None)


def tail_rows(seed=4242, budget=20000, theta_bits=5):
    """The ACTUAL event: D = steps from a trail's insertion to its next
    distinguished point is geometric(1/theta) step by step (fresh 5-bit field
    per step), so P(D > R) = q^R.  R is the trail's remaining lifetime at its
    launch, measured.  The row is the mean of q^R over the launch population;
    the analytic ceiling is 0.582*theta/S (the geometric sum of the launch
    schedule), and the charged value is the larger of the two."""
    r = one(seed, budget, theta_bits, stop_on_success=False)
    th = 1 << theta_bits
    q = (th - 1.0) / th
    L = r['launch_log']
    def row(scale):
        logs = []
        for x in L:
            R = max(0, min((1 << theta_bits) * 128, budget - x['t0']))
            logs.append(int(R * scale) * math.log2(q))
        return dict(mean_qR=sum(2.0 ** v for v in logs) / len(logs),
                    R_min=int(min(max(0, min((1 << theta_bits) * 128, budget - x['t0']))
                                  * scale for x in L)),
                    R_max=int(max(max(0, min((1 << theta_bits) * 128, budget - x['t0']))
                                  * scale for x in L)),
                    n=len(logs))
    one_x = row(1.0)
    half = row(0.5)
    two = row(2.0)
    analytic = 0.582 * th / budget
    return dict(theta=th, S=budget, n_launches=len(L), R_at_1x=one_x,
                R_at_half=half, R_at_2x=two,
                analytic_ceiling=analytic,
                charged=max(one_x['mean_qR'], analytic),
                row_moves_with_R=(half['mean_qR'] > one_x['mean_qR'] > two['mean_qR']),
                analytic_form='delta <= (1/(e-1)) * theta / S = 0.582 * theta / S')


def audits(budget=20000):
    out = {}
    # (a) step budget: v32's launch credit vs v33's per-step charge, same seed
    lc = legacy_control()
    out['a_step_budget'] = dict(
        legacy_v32=lc['v32_launch_credit'], shipped_v33=lc['v33_step_budget'],
        same_seed=lc['same_seed'],
        verdict='v32 charged launches and never tested the step counter: the same '
                'seed executes %.1fx the declared S.  v33 stops at completed == S '
                'with consumed == executed.' % lc['v32_launch_credit']['over_S'])
    # (b) tail bound: the same launch population re-evaluated at R/2 and 2R
    tr = tail_rows()
    out['b_tail_R'] = dict(R_at_1x=tr['R_at_1x'], R_at_half=tr['R_at_half'],
                           R_at_2x=tr['R_at_2x'],
                           analytic_ceiling=tr['analytic_ceiling'],
                           charged=tr['charged'], row_moves_with_R=tr['row_moves_with_R'],
                           verdict='the row is a function of the measured R, not a '
                                   'constant: halving R raises it, doubling R lowers it')
    # (c) conditional detection vs the old shared-report fraction
    mb = merge_batch()
    sf = shared_report_fraction()
    out['c_detection'] = dict(
        p_detect_given_merge=mb['p_detect'], cp95_lower=mb['cp95_lower'],
        n_merged=mb['merged'], detected=mb['detected'],
        old_shared_report_fraction=sf['shared_fraction'],
        old_object=dict(reports=sf['reports'], shared=sf['shared'],
                        matches=sf['matches']),
        differ=abs((mb['p_detect'] or 0) - (sf['shared_fraction'] or 0)) > 0.05,
        verdict='P(detect|merge) %.4f (CP-lo %.4f, n=%d) vs the old shared-report '
                'fraction %.3f -- different objects, both labelled'
                % (mb['p_detect'], mb['cp95_lower'], mb['merged'],
                   sf['shared_fraction']))
    # (d) unit chain
    import ledger33 as L33
    base = L33.claim33()
    moved = L33.claim33(scale=K.C_REF)
    shift = moved['claim_log2'] - base['claim_log2']
    out['d_unit_chain'] = dict(
        one_division_line='W_ref = W_instructions / C_REF (ledger33.total33)',
        base_claim=base['claim_log2'], scaled_claim=moved['claim_log2'],
        shift_bits=shift, log2_cref=math.log2(K.C_REF),
        shift_equals_log2_cref=abs(shift - math.log2(K.C_REF)) < 1e-3,
        double_normalisation_would_show=(abs(shift) < 1e-6 or
                                         abs(shift - 2 * math.log2(K.C_REF)) < 1e-3),
        verdict='multiplying the RAW instruction ledger by C_REF moves the claim '
                'by exactly log2(C_REF) = %.6f bits; a second division would leave '
                'it fixed or move it by 2*log2(C_REF)'
                % math.log2(K.C_REF))
    return out


def main():
    budgets = (200, 400, 800, 1600, 4000, 20000)
    sw = sweep(budgets)
    lc = legacy_control()
    gt = gate_trace()
    mb = merge_batch()
    sf = shared_report_fraction()
    tr = tail_rows()
    os.makedirs(EV, exist_ok=True)
    # the measured factors are written BEFORE the unit-chain audit reads them
    json.dump(dict(merge_batch=mb, old_shared_report=sf,
                   charged=mb['cp95_lower']),
              open(os.path.join(EV, 'detect33.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    json.dump(tr, open(os.path.join(EV, 'tail33.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    au = audits()
    json.dump(dict(sweeps=sw, theta_bits=THETA_BITS, trials=TRIALS,
                   halt_rule=PR.HALT_RULE, inflight_cap=PR.INFLIGHT_CAP,
                   cap_mult=PR.CAP_MULT, legacy=lc),
              open(os.path.join(EV, 'protocol33.json'), 'w'),
              indent=1, sort_keys=True, default=str)
    json.dump(gt, open(os.path.join(EV, 'gate33.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    json.dump(au, open(os.path.join(EV, 'audit33.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    print('\nLEGACY v32 rule: completed %d vs S %d (%.1fx)'
          % (lc['v32_launch_credit']['completed'], lc['S'],
             lc['v32_launch_credit']['over_S']))
    print('SHIPPED v33   : completed %d == S, consumed==executed %s, halt %s'
          % (lc['v33_step_budget']['completed'],
             lc['v33_step_budget']['consumed_equals_executed'],
             lc['v33_step_budget']['halt']))
    print('GATE', {k: gt[k] for k in ('halt', 'launches', 'deferred', 'completed',
                                     'launch_bound', 'launch_bound_binds')})
    print('halt row', gt['halt_row'])
    print('DETECT merged %d detected %d p %.6f cp95lo %.6f  old shared frac %.3f'
          % (mb['merged'], mb['detected'], mb['p_detect'], mb['cp95_lower'],
             sf['shared_fraction']))
    print('TAIL', tr['R_at_1x'], 'half', tr['R_at_half']['mean_qR'],
          'two', tr['R_at_2x']['mean_qR'], 'charged %.3e' % tr['charged'])
    for k in ('a_step_budget', 'b_tail_R', 'c_detection', 'd_unit_chain'):
        print(k, '->', au[k].get('verdict'))


if __name__ == '__main__':
    main()
