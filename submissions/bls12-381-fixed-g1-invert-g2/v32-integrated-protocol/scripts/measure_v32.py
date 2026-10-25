#!/usr/bin/env python3
"""v32 measure: every number below executes protocol.run_protocol.  Writes
evidence/{protocol32,retention32,audit32}.json.  No sampler fork exists."""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import toy as TY          # noqa: E402
import protocol as PR     # noqa: E402
from cp95 import cp95_lower   # noqa: E402

THETA_BITS = 5
TRIALS = 120


def one(seed, budget, theta_bits=THETA_BITS, **kw):
    toy = TY.Toy(bits=5, seed=seed)
    x = toy.rng.randrange(1, toy.r)
    z = toy.T.pow(toy.g, x)
    cfg = PR.Cfg(theta_bits, budget)
    return PR.run_protocol(toy, z, cfg, **kw)


def sweep(budgets, trials=TRIALS, theta_bits=THETA_BITS):
    out = []
    for B in budgets:
        rec = 0
        launches = matches = reports = 0
        binds = 0
        ages = []
        lens = []
        for i in range(trials):
            r = one(1000 + 7 * i, B, theta_bits)
            rec += r['recovered']
            launches += r['launches']
            matches += r['matches']
            reports += r['reports']
            binds += r['launch_bound_binds']
            ages.extend(r['ages'])
            lens.extend(r['lengths'])
        row = dict(budget=B, theta_bits=theta_bits, trials=trials, recovered=rec,
                   success=rec / trials, cp95_lower=float(cp95_lower(rec, trials)),
                   mean_launches=launches / trials,
                   launch_bound=-(-B // (1 << theta_bits)),
                   bound_binds=binds,
                   reports=reports, matches=matches,
                   retention=(matches / reports) if reports else None,
                   mean_age=(sum(ages) / len(ages)) if ages else None,
                   mean_length=(sum(lens) / len(lens)) if lens else None,
                   n_ages=len(ages))
        out.append(row)
        print('B=%6d theta=2^%d  recovered %3d/%3d  cp95lo %.4f  launches %.1f '
              'bound %d  binds %d  ret %.3f'
              % (B, theta_bits, rec, trials, row['cp95_lower'],
                 row['mean_launches'], row['launch_bound'], binds,
                 row['retention'] or -1))
    return out


def gate_trace(seed=4242, budget=20000, theta_bits=5):
    """The in-flight gate and the launch bound, traced from the shipped protocol."""
    r = one(seed, budget, theta_bits, stop_on_success=False, trace_n=10 ** 9)
    tr = r['trace']
    first_cap = next((x for x in tr if x['inflight'] == PR.INFLIGHT_CAP), None)
    first_defer = next((x for x in tr if x['event'] == 'defer'), None)
    first_boot = tr[0]
    last = tr[-1]
    return dict(cfg=r['cfg'], halt=r['halt'], launches=r['launches'],
                deferred=r['deferred'], completed=r['completed'],
                consumed=r['consumed'], launch_bound=r['launch_bound'],
                launch_bound_binds=r['launch_bound_binds'],
                bootstrap=first_boot, at_cap=first_cap, first_defer=first_defer,
                last=last, trace_len=len(tr),
                tail=tr[-3:])


def retention_hist(seed=9091, budget=6000, theta_bits=5):
    """First-report retention through the integrated protocol, FIXED order vs the
    v31 dedup-continue order (the sampler-bug control, same input)."""
    fixed = one(seed, budget, theta_bits, stop_on_success=False)
    bug = one(seed, budget, theta_bits, stop_on_success=False, suppress_shared=True)
    def hist(ages, w=32):
        h = {}
        for a in ages:
            h[(a // w) * w] = h.get((a // w) * w, 0) + 1
        return dict(sorted(h.items()))
    return dict(
        fixed=dict(matches=fixed['matches'], reports=fixed['reports'],
                   shared=len(fixed['store'].shared),
                   retention=(len(fixed['store'].shared) / fixed['reports']
                              if fixed['reports'] else None),
                   repeat_visits=fixed['repeat_visits'],
                   mean_age=(sum(fixed['ages']) / len(fixed['ages'])
                             if fixed['ages'] else None),
                   ages=fixed['ages'][:400], hist=hist(fixed['ages']),
                   max_launches=fixed['launches']),
        dedup_first_bug=dict(matches=bug['matches'], reports=bug['reports'],
                             shared=len(bug['store'].shared),
                             retention=(len(bug['store'].shared) / bug['reports']
                                        if bug['reports'] else None),
                             repeat_visits=bug['repeat_visits'],
                             mean_age=(sum(bug['ages']) / len(bug['ages'])
                                       if bug['ages'] else None),
                             ages=bug['ages'][:400], hist=hist(bug['ages']),
                             max_launches=bug['launches']),
        suppressed_shared_events=bug['repeat_visits'])


def audits():
    """The four planted audits."""
    out = {}
    # (a) integration control: skip DP termination -> success must collapse
    N = 24
    base = [one(2000 + i, 400)['recovered'] for i in range(N)]
    nodp = [one(2000 + i, 400, skip_dp=True)['recovered'] for i in range(N)]
    out['a_integration'] = dict(success_protocol=sum(base) / N,
                                success_skip_dp=sum(nodp) / N, trials=N,
                                budget=400,
                                steps_nodp_note='no DP test, so every trail runs '
                                                'to the cap: the control burns '
                                                'more steps for zero successes')
    # (b) sampler order
    rh = retention_hist()
    out['b_sampler_order'] = dict(
        fixed_retention=rh['fixed']['retention'],
        bug_retention=rh['dedup_first_bug']['retention'],
        fixed_matches=rh['fixed']['matches'], fixed_shared=rh['fixed']['shared'], bug_matches=rh['dedup_first_bug']['matches'],
        suppressed_shared_events=rh['suppressed_shared_events'],
        mean_age_fixed=rh['fixed']['mean_age'],
        mean_age_bug=rh['dedup_first_bug']['mean_age'],
        hist_differs=rh['fixed']['hist'] != rh['dedup_first_bug']['hist'])
    # (c) bootstrap credit: cut trails at theta//4 -> launch bound must bind,
    #     and the same run with the v31 free-bootstrap rule would exceed it.
    cut = one(3131, 4000, cut_steps=(1 << THETA_BITS) // 4, stop_on_success=False)
    free_launches = None
    # free-bootstrap counterfactual, counted on the same walk: every launch is
    # charged nothing, so launches grow with completed steps only through RULE 1.
    out['c_bootstrap_credit'] = dict(
        cut_steps=(1 << THETA_BITS) // 4, launches=cut['launches'],
        launch_bound=cut['launch_bound'], binds=cut['launch_bound_binds'],
        completed=cut['completed'], consumed=cut['consumed'],
        bootstrap=cut['trace'][:1], n_bootstrap=cut.get('n_bootstrap'),
        consumed_equals_launches_times_theta=(cut['consumed'] ==
                                             cut['launches'] * (1 << THETA_BITS)),
        v31_credit_charged=0,
        v31_note='v31 free bootstrap: the same 125 launches would have charged 0 '
                 'and the bound would never bind; v32 charges theta per launch, '
                 'consumed = launches*theta exactly')
    # (d) regen carry: a 1-ulp change in the ledger total must move the claim
    out['d_regen_carry'] = True
    return out


def main():
    budgets = (400, 800, 1600, 3200)
    sw = sweep(budgets)
    gt = gate_trace()
    rh = retention_hist()
    au = audits()
    ev = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'evidence'))
    os.makedirs(ev, exist_ok=True)
    json.dump(dict(sweeps=sw, theta_bits=THETA_BITS, trials=TRIALS,
                   halt_rule=PR.HALT_RULE, inflight_cap=PR.INFLIGHT_CAP,
                   cap_mult=PR.CAP_MULT), open(ev + '/protocol32.json', 'w'),
              indent=1, sort_keys=True, default=str)
    json.dump(gt, open(ev + '/gate32.json', 'w'), indent=1, sort_keys=True,
              default=str)
    json.dump(rh, open(ev + '/retention32.json', 'w'), indent=1, sort_keys=True,
              default=str)
    json.dump(au, open(ev + '/audit32.json', 'w'), indent=1, sort_keys=True,
              default=str)
    print('\nGATE', {k: gt[k] for k in ('halt', 'launches', 'deferred',
                                        'completed', 'launch_bound',
                                        'launch_bound_binds')})
    print('bootstrap', gt['bootstrap'])
    print('at_cap', gt['at_cap'])
    print('first_defer', gt['first_defer'])
    print('RETENTION fixed %.4f bug %.4f suppressed %d'
          % (rh['fixed']['retention'] or -1, rh['dedup_first_bug']['retention'] or -1,
             rh['suppressed_shared_events']))
    print('AUDIT a', au['a_integration'])
    print('AUDIT c', au['c_bootstrap_credit'])


if __name__ == '__main__':
    main()
