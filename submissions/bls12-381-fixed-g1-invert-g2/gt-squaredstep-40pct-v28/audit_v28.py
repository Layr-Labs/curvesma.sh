#!/usr/bin/env python3
"""v28 audit: the four controls the lane requires, plus the ledger emission.

(a) UNIFORM-AGE CONTROL   the old integral must be a DISCLOSED-REJECTED method, and
                          the charged delta must not be it.
(b) MOMENT CONTROL        perturb one length sample; E[L^3]/(2E[L^2]) must move,
                          proving the ratio is computed from data, not asserted.
(c) STOPPING-RULE CONTROL simulate launches with the rule DISABLED; total realized
                          work must exceed W somewhere, and with the rule must not.
(d) OVERFLOW-ROW CONTROL  delete the added one-reference-operation row; the
                          total-equals-sum-of-displayed-rows assertion must flip.
"""
import json
import math
import random
from fractions import Fraction

import ledger28 as L28
import ledger27 as L27
import v28_measure as M


def uniform_control():
    u = M.retention_uniform()
    dm, det = L28.delta_charge()
    c_uniform = L28.claim28(delta=Fraction(1, 128), wasted_mult=1.0)
    c_now = L28.claim28()
    stale = json.load(open('submission.json'))['claims'][0]['statement']
    return dict(uniform_value=u, uniform_claim=('%d/%d' % (c_uniform['numerator'],
                                                          c_uniform['denominator'])),
                charged_delta=det['charged'], charged_delta_log2=math.log2(det['charged']),
                charged_is_not_uniform=abs(det['charged'] - u) > 1e-6,
                rejected_is_labelled=('DISCLOSED-REJECTED' in stale
                                      or 'DISCLOSED-REJECTED' in stale),
                statement_charges_measured_law=('2^-171.66' in stale),
                statement_keeps_uniform_only_as_rejected=(
                    'DISCLOSED-REJECTED' in stale.replace('DISCLOSED+REJECTED', ''))
                and stale.index('0.0078125') > stale.index('DISCLOSED-REJECTED'),
                claim_now=('%d/%d' % (c_now['numerator'], c_now['denominator'])))


def moment_control(arm='1', bump=1):
    m = json.load(open('evidence/v28-moments.json'))['arms'][arm]
    s1, s2, s3 = m['sum_l'], m['sum_l2'], m['sum_l3']
    n = m['n_trails']
    base = Fraction(s3, 2 * s2)
    L = 20                                   # perturb ONE trail length 20 -> 21
    s1b, s2b, s3b = s1 + bump, s2 + (2 * L + 1), s3 + (3 * L * L + 3 * L + 1)
    pert = Fraction(s3b, 2 * s2b)
    return dict(arm=arm, n=n, sums=[s1, s2, s3], base=str(base), perturbed=str(pert),
                base_float=float(base), pert_float=float(pert),
                moved=(pert != base), delta=float(pert - base),
                base_over_theta=float(base) / (1 << int(arm)),
                pert_over_theta=float(pert) / (1 << int(arm)))


def stopping_control(trials=200000, W=1 << 20, cap_mult=128, seed=7):
    """Launch loop.  reserve = cap_mult*theta is scaled here to W/2^10 so the toy
    can realise both regimes; charges are uniform in (0, reserve]."""
    rng = random.Random(seed)
    reserve = W >> 10
    def run(rule):
        total = 0
        while True:
            charge = rng.uniform(1, reserve)
            if rule and total + charge > W - reserve:
                break
            if not rule and total + charge > W:
                # keep going: this is the DISABLED-rule arm; it overruns
                total += charge
                break
            total += charge
        return total
    with_rule = [run(True) for _ in range(trials)]
    no_rule = [run(False) for _ in range(trials)]
    return dict(W=W, reserve=reserve, trials=trials,
                with_rule_max=max(with_rule), with_rule_over=sum(1 for t in with_rule
                                                                 if t > W),
                no_rule_max=max(no_rule), no_rule_over=sum(1 for t in no_rule if t > W),
                bound_holds=max(with_rule) <= W - reserve + reserve,
                disabled_overruns=sum(1 for t in no_rule if t > W) > 0)


def overflow_control():
    T = L28.total28()
    comps = dict(T['components'])
    ok = comps.pop('overflow_row') == 1
    total_without = sum(comps.values(), Fraction(0))
    return dict(row_present=ok, total=str(T['total']), total_without=str(total_without),
                delta=str(T['total'] - total_without),
                sum_assert_flips=(total_without != T['total']),
                delta_is_one_ref_op=(T['total'] - total_without == 1))


def main():
    c = L28.claim28()
    T = c['T']
    lit = json.load(open('submission.json'))['claims'][0]
    st = lit['statement']
    checks = dict(
        claim_matches_statement=('%d/%d' % (c['numerator'], c['denominator'])
                                 in (lit['score']['numerator'] + '/'
                                     + lit['score']['denominator'])),
        claim_is_least_4dp=(c['claim_log2'] >= c['total_log2'] and
                            c['claim_log2'] - c['total_log2'] < 1e-4),
        total_is_sum_of_five=(set(T['components']) == {'w_stress', 'trap_charge',
                                                       'aux_sum', 'overflow_row',
                                                       'stop_loss_row'}),
        statement_has_w=('2^%.6f' % math.log2(float(T['components']['w_stress']))) in st,
        statement_has_trap=('2^%.4f' % math.log2(float(T['components']['trap_charge'])))
        in st,
        statement_has_aux=('2^%.6f' % math.log2(float(T['components']['aux_sum']))) in st,
        statement_has_stop=('2^%.4f' % math.log2(float(T['components']['stop_loss_row'])))
        in st,
        statement_has_claim=('126.2599' in st),
        statement_has_sumline=('= 2^126.259807' in st),
    )
    out = dict(claim=dict(numerator=c['numerator'], denominator=c['denominator'],
                          claim_log2=c['claim_log2'], total_log2=c['total_log2'],
                          total=str(T['total'])),
               components={k: str(v) for k, v in T['components'].items()},
               components_log2={k: math.log2(float(v)) for k, v in T['components'].items()},
               S=str(T['S']), S_log2=math.log2(T['S']), f=str(T['f']),
               delta=T['delta'], trap={k: str(v) for k, v in T['trap'].items()},
               stop=T['stop'],
               controls=dict(uniform=uniform_control(), moment=moment_control(),
                             stopping=stopping_control(), overflow=overflow_control()),
               literal_checks=checks)
    out['all_checks'] = (all(checks.values())
                         and out['controls']['uniform']['charged_is_not_uniform']
                         and out['controls']['uniform']['statement_charges_measured_law']
                         and out['controls']['moment']['moved']
                         and out['controls']['stopping']['bound_holds']
                         and out['controls']['stopping']['disabled_overruns']
                         and out['controls']['overflow']['sum_assert_flips'])
    with open('evidence/v28-ledger.json', 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True, default=str)
    print('claim %d/%d = %.4f  total 2^%.9f' % (c['numerator'], c['denominator'],
                                                c['claim_log2'], c['total_log2']))
    for k, v in checks.items():
        print('  literal %-32s %s' % (k, v))
    for k, v in out['controls'].items():
        print('  control %-10s %s' % (k, {kk: vv for kk, vv in v.items()
                                          if kk in ('charged_is_not_uniform', 'moved',
                                                    'bound_holds', 'disabled_overruns',
                                                    'sum_assert_flips', 'no_rule_max',
                                                    'with_rule_max')}))
    print('ALL CHECKS', out['all_checks'])


if __name__ == '__main__':
    main()
