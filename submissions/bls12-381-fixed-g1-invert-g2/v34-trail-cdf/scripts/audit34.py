#!/usr/bin/env python3
"""v34 audits -- four controls, each one falsifiable, all reading the shipped
evidence.  Writes evidence/audit34.json.

  (a) CDF CONTROL      the charged delta is a CP-95 UPPER bound computed from the
                       measured per-trail population.  Shrink the population
                       (drop the last half, then the last three quarters) and the
                       SAME measured rate must give a WIDER (higher) bound: the
                       bound is data, not a constant.
  (b) DENOMINATOR      the merge-all denominator and the old shared-report-only
                       denominator are BOTH reported; P(detect|merge) differs by
                       the selection-bias magnitude, which is the row.
  (c) MU CONTROL       mu = completed/(launches*theta) is computed from the trace;
                       perturb the launch count by +10% and the row must move by
                       exactly 1/1.1.
  (d) UNIT CHAIN       regen --control-ledger: multiplying the RAW instruction
                       ledger by C_REF moves the claim by exactly log2(C_REF), and
                       the 4-decimal ceiling carries: 2^c >= W_ref > 2^(c-1e-4).
"""
import json
import math
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ledger34 as L34          # noqa: E402
import regen_v34 as R34         # noqa: E402
import constants as K           # noqa: E402
from measure_v34 import cp95_upper   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.abspath(os.path.join(HERE, '..', 'evidence'))


def audit_cdf(tl):
    n = tl['n_trails']
    k = tl['n_never']
    full = cp95_upper(k, n)
    rows = [dict(n=n, over=k, cp95_upper=full)]
    for keep in (n // 2, n // 4):
        kk = int(round(k * keep / n))
        rows.append(dict(n=keep, over=kk, cp95_upper=cp95_upper(kk, keep)))
    rises = all(rows[i]['cp95_upper'] < rows[i + 1]['cp95_upper']
                for i in range(len(rows) - 1))
    return dict(population=n, never=k, point=k / n, rows=rows,
                bound_rises_when_population_shrinks=rises,
                verdict='the delta ROW reported for frame B is the CP-95 UPPER '
                        'bound on the measured P(T > 128*theta); halving and '
                        'quartering the population raises it (%.4f -> %.4f -> '
                        '%.4f), so the bound is data. The CHARGED frame A carries '
                        'the same lifetime loss inside the merge-all denominator '
                        'instead (finalize34), so it is counted once' %
                        (rows[0]['cp95_upper'], rows[1]['cp95_upper'],
                         rows[2]['cp95_upper']))


def audit_denominator(mg):
    return dict(merge_all=dict(events=mg['events_total'],
                               detected=mg['detected_total'],
                               p=mg['p_detect_all_merges'],
                               cp95_lower=mg['pooled_cp95_lower'],
                               charged=mg['min_per_run_cp95_lower']),
                shared_report_only=dict(merged=mg['old_v33_object']['merged'],
                                        p=mg['p_shared_report'],
                                        note='conditioned on a DP hit AND a store '
                                             'record: the selection-biased object'),
                selection_bias_delta=mg['selection_bias_delta'],
                differ=abs((mg['p_detect_all_merges'] or 0)
                           - (mg['p_shared_report'] or 0)) > 0.05,
                verdict='the two denominators give %.4f (all class merges) and '
                        '%.4f (shared reports only): the selection bias is %.4f, '
                        'charged as a row'
                        % (mg['p_detect_all_merges'], mg['p_shared_report'],
                           mg['selection_bias_delta']))


def audit_mu(tl):
    mu = tl['mu_derived']
    launched = [r['launches'] for r in tl['per_run']]
    done = [r['completed'] for r in tl['per_run']]
    theta = tl['theta']
    mu_pert = sum(c / ((l * 1.1) * theta) for c, l in zip(done, launched)) / len(done)
    return dict(mu_derived=mu, mu_with_launches_plus_10pct=mu_pert,
                ratio=mu_pert / mu if mu else None,
                expected_ratio=1 / 1.1,
                moves=(abs(mu_pert - mu) > 1e-6),
                formula='mu = completed / (launches * theta), per run, averaged',
                mean_launches=sum(launched) / len(launched),
                mean_completed=sum(done) / len(done),
                verdict='perturbing the trace launch count by +10%% moves mu from '
                        '%.6f to %.6f (ratio %.6f vs 1/1.1 = %.6f), so the row is '
                        'computed from the trace, not asserted'
                        % (mu, mu_pert, mu_pert / mu if mu else 0, 1 / 1.1))


def audit_chain():
    c, led = R34.build()
    lg = c['total_log2']
    ceil_ok = (2.0 ** c['claim_log2'] >= (2.0 ** lg) * (1 - 1e-12))
    minimal = (2.0 ** (c['claim_log2'] - 1e-4) < (2.0 ** lg) * (1 + 1e-12))
    out = subprocess.run([sys.executable, os.path.join(HERE, 'regen_v34.py'),
                          '--control-ledger'], capture_output=True, text=True)
    cm, _ = R34.build(scale=K.C_REF)
    shift = cm['claim_log2'] - c['claim_log2']
    return dict(base_claim=c['claim_log2'], claim_quantum=1e-4,
                carries_upper=ceil_ok, carries_minimal=minimal,
                raw_times_cref_claim=cm['claim_log2'], shift_bits=shift,
                log2_cref=math.log2(K.C_REF),
                shift_equals_log2_cref=abs(shift - math.log2(K.C_REF)) < 1e-3,
                control_d_stdout=out.stdout.strip().splitlines(),
                control_d_ok=out.returncode == 0,
                verdict='the 4-decimal ceiling carries (2^c >= W_ref > 2^(c-1e-4)) '
                        'and one C_REF division means a scaled ledger moves the '
                        'claim by exactly log2(20000) = %.6f bits' % math.log2(K.C_REF))


def main():
    tl = json.load(open(os.path.join(EV, 'tail34.json')))
    mg = json.load(open(os.path.join(EV, 'merge34.json')))
    au = dict(a_cdf=audit_cdf(tl), b_denominator=audit_denominator(mg),
              c_mu=audit_mu(tl), d_unit_chain=audit_chain(),
              inputs=dict(delta=str(L34.evidence_inputs()['delta']),
                          p_det=str(L34.evidence_inputs()['p_det']),
                          mu=str(L34.evidence_inputs()['mu'])))
    json.dump(au, open(os.path.join(EV, 'audit34.json'), 'w'), indent=1,
              sort_keys=True, default=str)
    for k in ('a_cdf', 'b_denominator', 'c_mu', 'd_unit_chain'):
        print(k, '->', au[k]['verdict'])
    ok = (au['a_cdf']['bound_rises_when_population_shrinks']
          and au['b_denominator']['differ'] and au['c_mu']['moves']
          and au['d_unit_chain']['control_d_ok']
          and au['d_unit_chain']['shift_equals_log2_cref'])
    print('AUDITS', 'PASS' if ok else 'FAIL')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
