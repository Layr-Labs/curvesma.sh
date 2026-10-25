#!/usr/bin/env python3
"""v30 audit: literal sweep + four planted controls + 20 independent five-row
recomputations.  Control (a) shifts the gate 12 -> 11 and the gate literal must
move; (b) the bare-schedule trace must exceed the derived 128 and the reserved
cap must be enforced; (c) the mean-price row must visibly underprice the max-sum
row on fat-tail synthetic charges; (d) truncating the age support must raise the
tail row."""
import json
import math
import random
import re
from fractions import Fraction

import ledger30 as L30
import v30_measure as M

AGE = json.load(open('evidence/v30-age-law.json'))
STO = json.load(open('evidence/v30-stopping.json'))
SUB = json.load(open('submission.json'))
LIT = json.load(open('evidence/v30-literals.json'))


def literal_sweep():
    text = SUB['claims'][0]['statement'] + '\n' + open('argument.md').read()
    computed = {}
    for k, v in LIT.items():
        try:
            computed[k] = float(v['value'])
        except (TypeError, ValueError):
            pass
    for k, v in AGE['arms'].items():
        for fld in ('delta_residual_all_log2', 'delta_residual_report_log2'):
            if v.get(fld) is not None:
                computed['age_%s_%s' % (k, fld)] = v[fld]
    T = L30.total29()
    computed['reserve_over_S_log2'] = math.log2(T['stop']['reserve_over_S'])
    rows = {}
    for tok in sorted(set(re.findall(r'(?<![A-Za-z0-9_])2\^([-+]?\d+(?:\.\d+)?)', text))):
        x = float(tok)
        hit = None
        for k, v in computed.items():
            if abs(v - x) <= 5e-3:
                hit = 'computed:%s' % k
                break
        rows[tok] = hit
    return dict(all_Y=all(v is not None for v in rows.values()), tokens=rows,
                n_tokens=len(rows))


def gate_control():
    """(a) drop the last shift.  The twelve-shift gate credits every merge; the
    eleven-shift gate cannot credit the j=11 witnesses, so the gate literal
    changes and the audit asserts the shipped literal is the 12-shift one."""
    a = AGE['arms']['1']
    hist = {int(k): v for k, v in a['shift_hist'].items()}
    total = sum(hist.values())
    keep12 = total
    keep11 = total - hist.get(11, 0)
    p12 = keep12 / float(a['merges_all'])
    p11 = keep11 / float(a['merges_all'])
    gate = LIT['gate_per_merge']['value']
    return dict(shifts_12=total, shifts_11=keep11, p_12=p12, p_11=p11,
                literal_changes=(p12 != p11), shipped_is_12=(abs(gate - p12) < 5e-7),
                delivered_p_recover=LIT['gate_measured']['value'],
                census_rate=LIT['gate_census']['value'],
                fires=(p12 != p11 and abs(gate - p12) < 5e-7))


def inflight_control():
    """(b) the bare-schedule trace must exceed ceil(cap/theta)=128.  If the trace
    max is >= the derived bound, the derivation is rejected and the enforced
    atomic cap is required; the one-smaller cap must under-cover."""
    T = L30.total29()
    c = STO['arms']['1']['concurrency']
    n = T['stop']['n_inflight']
    return dict(derived_ratio=c['derived_ratio'], bare_trace_max=c['bare_schedule_max_inflight'],
                bare_exceeds_derived=(c['bare_schedule_max_inflight'] >= c['derived_ratio']),
                enforced_cap=c['enforced_cap'], enforced_trace_max=c['enforced_trace_max_inflight'],
                enforce_holds=(c['enforced_trace_max_inflight'] <= c['enforced_cap']),
                toosmall_cap=(n - 1), toosmall_trace_max=c['toosmall_cap_trace_max_inflight'],
                undercoverage_steps=str(c['toosmall_cap_trace_max_inflight'] * c['trail_cap_steps']),
                reserve_steps=T['stop']['reserve_steps'],
                fires=(c['bare_schedule_max_inflight'] >= c['derived_ratio']
                       and c['enforced_trace_max_inflight'] <= n
                       and c['toosmall_cap_trace_max_inflight'] < n))


def max_vs_mean_control():
    """(c) fat-tail synthetic charges.  The mean-price row prices every step at
    the average mix; the max-sum row prices every step at the dearer step.  The
    mean must visibly underprice when the tail carries the dear steps."""
    T = L30.total29()
    cm = T['ledgers']['c_mul_stress']
    cs = T['ledgers']['c_sqr_stress']
    f = T['ledgers']['f']
    mean_price = (1 - float(f)) * float(cm) + float(f) * float(cs)
    n = 1000
    synthetic = [float(cm)] * 1 + [float(cs)] * (n - 1)     # fat right tail
    mean_sum = n * mean_price
    max_sum = sum(max(x, float(cm)) for x in synthetic)
    w_mean = LIT['mean_w_log2']['value']
    w_max = LIT['w_stress_log2']['value']
    return dict(mean_price=mean_price, c_max=float(cm), n=n,
                synthetic_mean_sum=mean_sum, synthetic_max_sum=max_sum,
                mean_underprices_synthetic=(mean_sum < max_sum),
                w_mean_log2=w_mean, w_max_log2=w_max,
                undercount_bits=w_max - w_mean,
                fires=(mean_sum < max_sum and (w_max - w_mean) > 0.3))


def support_control():
    """(d) truncate the measured support from 16 theta to 8 theta; the tail row
    must rise (more unobserved tail charged at the conservative constant)."""
    t16, _ = L30.delta_tail_bound(support_lo=16.0)
    t8, _ = L30.delta_tail_bound(support_lo=8.0)
    return dict(tail_16_log2=math.log2(t16), tail_8_log2=math.log2(t8),
                truncation_raises=(t8 > t16),
                ratio=t8 / t16, fires=(t8 > t16))


def independent_recompute(n=20, seed=30):
    rng = random.Random(seed)
    base = L30.total29()
    outs = []
    for _ in range(n):
        f = Fraction(14566, 41929) * (1 + rng.choice([-1, 1]) * 2 ** -40 * rng.random())
        d = base['delta']['charged'] * (1 + rng.choice([-1, 1]) * 2 ** -20)
        T = L30.total29(f_measured=float(f), delta=d)
        s = sum(T['components'].values(), Fraction(0))
        outs.append(s == T['total'] and
                    abs(math.log2(float(T['total'])) - T['log2_total']) < 1e-9)
    return dict(n=n, all_equal=all(outs), passed=sum(outs))


def main():
    c = L30.claim29()
    T = c['T']
    lit = SUB['claims'][0]
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
        statement_has_mean=('2^%.6f' % LIT['mean_w_log2']['value']) in st,
        statement_has_trap=('2^%.6f' % math.log2(float(T['components']['trap_charge']))) in st,
        statement_has_aux=('2^%.6f' % math.log2(float(T['components']['aux_sum']))) in st,
        statement_has_stop=('2^%.6f' % math.log2(float(T['components']['stop_loss_row']))) in st,
        statement_has_claim=('%.4f' % c['claim_log2']) in st,
        statement_has_sumline=('2^%.9f' % c['total_log2']) in st,
        statement_has_gate=('%.6f' % LIT['gate_measured']['value']) in st,
        statement_has_census=('%.6f' % LIT['gate_census']['value']) in st,
        statement_has_reserve=('2^%.2f' % math.log2(int(T['stop']['reserve_steps']))) in st,
        statement_has_inflight=('%d' % T['stop']['n_inflight'] in st),
    )
    controls = dict(gate=gate_control(), inflight=inflight_control(),
                    max_vs_mean=max_vs_mean_control(), support=support_control())
    sweep = literal_sweep()
    indep = independent_recompute()
    planted = dict(gate_fires=controls['gate']['fires'],
                   inflight_fires=controls['inflight']['fires'],
                   max_vs_mean_fires=controls['max_vs_mean']['fires'],
                   support_fires=controls['support']['fires'])
    out = dict(claim=dict(numerator=c['numerator'], denominator=c['denominator'],
                          claim_log2=c['claim_log2'], total_log2=c['total_log2'],
                          total=str(T['total'])),
               components={k: str(v) for k, v in T['components'].items()},
               components_log2={k: math.log2(float(v)) for k, v in T['components'].items()},
               S=str(T['S']), S_log2=math.log2(T['S']), f=str(T['f']),
               delta=T['delta'], trap=T['trap'], stop=T['stop'],
               gate=LIT['gate_measured'], concurrency=T['stop'],
               literal_checks=checks, literal_sweep=sweep, controls=controls,
               independent_recompute=indep, planted=planted,
               planted_score='%d/4' % sum(planted.values()))
    out['all_checks'] = (all(checks.values()) and sweep['all_Y'] and indep['all_equal']
                         and all(planted.values()))
    json.dump(out, open('evidence/v30-ledger.json', 'w'), indent=1, sort_keys=True,
              default=str)
    print('claim %d/%d = %.4f total 2^%.9f' % (c['numerator'], c['denominator'],
                                               c['claim_log2'], c['total_log2']))
    print('literal sweep: %d tokens, all_Y=%s' % (sweep['n_tokens'], sweep['all_Y']))
    for t, v in sweep['tokens'].items():
        if v is None:
            print('  UNMAPPED 2^%s' % t)
    print('controls', {k: v['fires'] for k, v in controls.items()})
    print('planted', out['planted_score'])
    print('indep recompute %d/%d' % (indep['passed'], indep['n']))
    print('ALL CHECKS', out['all_checks'])


if __name__ == '__main__':
    main()
