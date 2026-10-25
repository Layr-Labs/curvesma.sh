#!/usr/bin/env python3
"""v34 finalize: turns the per-run rows into CLUSTER-LEVEL charges.

WHY THIS FILE EXISTS.  Verdict-v33 objection 2: the 7-run pool had no
independence justification.  The class map is deterministic, so events inside a
run are NOT exchangeable Bernoulli draws, and the run is the natural cluster.
Two things follow, and both are done here:

  delta   P(T > 128*theta) is charged at the MAXIMUM over runs of its CP-95 upper
          bound (a penalty, so the conservative direction is UP and the cluster
          is the run, not the trail).
  p_det   two factorizations of the same measured loss are computed and BOTH are
          written down:
            frame A  p_all = the merge-all rate, cluster-lower-bounded
                     (delta = 0: the whole loss lives in the merge event)
            frame B  (1 - delta_max) * p_det_cond, where p_det_cond is the
                     detection rate CONDITIONED on a report, cluster-lower-bounded
          The ledger charges whichever frame yields the larger S (more
          conservative), and argument.md states which one that was.

The cluster interval is the normal-approximation interval on the mean of the 24
per-run rates (mean - 1.96 * sd / sqrt(n_runs)), floored at 0.  It is stated as
such: an interval on the RUN-LEVEL MEAN, which is the quantity that transfers to a
production instance (the production z is one draw from the same seed law, so the
mean, not the best or worst run, is the estimand).
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cp95 import cp95_lower       # noqa: E402
from measure_v34 import cp95_upper  # noqa: E402

EV = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                  '..', 'evidence'))
Z = 1.96


def cluster(vals):
    n = len(vals)
    m = sum(vals) / n
    var = sum((v - m) ** 2 for v in vals) / (n - 1) if n > 1 else 0.0
    sd = math.sqrt(var)
    se = sd / math.sqrt(n)
    return dict(n=n, mean=m, sd=sd, se=se, lower=max(0.0, m - Z * se),
                upper=min(1.0, m + Z * se), z=Z,
                interval='normal approximation on the mean of the per-run rates')


def main():
    tl = json.load(open(os.path.join(EV, 'tail34.json')))
    mg = json.load(open(os.path.join(EV, 'merge34.json')))
    never_rate, p_cond = [], []
    for r in tl['per_run']:
        n = r['n_trails']
        nr = n - r['n_dp']
        never_rate.append(dict(seed=r['seed'], n=n, never=nr, rate=nr / n,
                               cp95_upper=cp95_upper(nr, n)))
        if r['merged']:
            p_cond.append(dict(seed=r['seed'], merged=r['merged'],
                               detected=r['detected'],
                               rate=r['detected'] / r['merged'],
                               cp95_lower=float(cp95_lower(r['detected'], r['merged']))))
    delta_max = max(x['cp95_upper'] for x in never_rate)
    delta_mean = cluster([x['rate'] for x in never_rate])
    cl_never = cluster([x['rate'] for x in never_rate])
    cl_never['cp_upper_of_mean'] = cp95_upper(
        int(sum(x['never'] for x in never_rate)),
        int(sum(x['n'] for x in never_rate)))
    cl_cond = cluster([x['rate'] for x in p_cond])
    cl_all = cluster([x['p'] for x in mg['per_run'] if x['p'] is not None])
    frames = dict(
        A=dict(p_all=cl_all['lower'], delta=0.0,
               note='merge-all rate, cluster lower bound, delta deleted'),
        B=dict(p_cond=cl_cond['lower'], delta=delta_max,
               note='(1 - delta_max) * conditioned detection'))
    frames['A']['effective'] = frames['A']['p_all']
    frames['B']['effective'] = (1 - frames['B']['delta']) * frames['B']['p_cond']
    tl['delta_rows'] = never_rate
    tl['delta_cluster'] = dict(max_cp95_upper=delta_max, mean_rate=cl_never['mean'],
                               mean_rate_interval=cl_never,
                               cp_upper_of_pooled=cl_never['cp_upper_of_mean'])
    tl['charged'] = delta_max
    tl['charged_kind'] = ('MAX over the {n} runs of the CP-95 UPPER bound on the '
                          'per-run P(T > 128*theta); a penalty is charged UP, and '
                          'the RUN is the cluster (the class map is deterministic, '
                          'so trails inside a run are not independent)'.format(n=len(never_rate)))
    mg['p_det_conditioned'] = dict(per_run=p_cond, cluster=cl_cond,
                                   charged=cl_cond['lower'],
                                   min_per_run_cp95_lower=min(x['cp95_lower'] for x in p_cond),
                                   pooled=float(cp95_lower(
                                       sum(x['detected'] for x in p_cond),
                                       sum(x['merged'] for x in p_cond))))
    mg['merge_all_cluster'] = cl_all
    mg['charged'] = cl_all['lower']
    mg['charged_kind'] = ('cluster lower bound on the mean of the {n} per-run '
                          'merge-all rates (normal approximation, z={z}); the '
                          'pooled binomial and the minimum per-run CP are both '
                          'reported beside it as the non-exchangeable diagnostics'
                          .format(n=cl_all['n'], z=Z))
    mg['frames'] = frames
    json.dump(tl, open(os.path.join(EV, 'tail34.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    json.dump(mg, open(os.path.join(EV, 'merge34.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    print('delta  max CP-upper %.6f   mean %.6f  interval [%.4f, %.4f]  pooled CP-up %.6f'
          % (delta_max, cl_never['mean'], cl_never['lower'], cl_never['upper'],
             cl_never['cp_upper_of_mean']))
    print('p_all  cluster mean %.6f  lower %.6f   (min per-run CP-lo %.6f, pooled %.6f)'
          % (cl_all['mean'], cl_all['lower'],
             min(x['cp95_lower'] for x in mg['per_run'] if x['cp95_lower'] is not None),
             mg['pooled_cp95_lower']))
    print('p_cond cluster mean %.6f  lower %.6f  min per-run CP-lo %.6f'
          % (cl_cond['mean'], cl_cond['lower'],
             min(x['cp95_lower'] for x in p_cond)))
    print('frames A effective %.6f   B effective %.6f' %
          (frames['A']['effective'], frames['B']['effective']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
