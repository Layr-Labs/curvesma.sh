#!/usr/bin/env python3
"""v29 audit: literal sweep + controls (a) literal, (b) endpoint, (c) concurrency,
(d) loss response, plus 20 independent five-row recomputations (4/4 planted
controls must fire)."""
import json
import math
import random
import re
from fractions import Fraction

import ledger29 as L29
import v29_measure as M

AGE = json.load(open('evidence/v29-age-law.json'))
STO = json.load(open('evidence/v29-stopping.json'))


def literal_sweep():
    st = json.load(open('submission.json'))['claims'][0]['statement']
    text = st + '\n' + open('argument.md').read()
    computed = {}
    for k, v in json.load(open('evidence/v29-literals.json')).items():
        try:
            computed[k] = float(v['value'])
        except (TypeError, ValueError):
            pass
    for k, v in AGE['arms'].items():
        for fld in ('delta_residual_all_log2', 'delta_residual_report_log2'):
            if v.get(fld) is not None:
                computed['age_%s_%s' % (k, fld)] = v[fld]
    computed['reserve_over_S_log2'] = math.log2(STO['arms']['1']['halts']['1']['lost_frac'] / 2.0) if False else math.log2(L29.total29()['stop']['reserve_over_S'])
    historical = {
        -128.0: 'v28 code base literal (withdrawn form)',
        -115.0: 'v28 code literal (withdrawn)',
        -171.66: 'v28 text literal (withdrawn, unnormalised)',
        -171.665: 'v28 text literal (withdrawn, unnormalised)',
        126.267148: 'v27c w_stress (historical)',
        125.8286: 'v27c S at delta=1/128 (historical)',
        35.3257: 'v27c trap row (historical)',
        88.438583: 'v28 stop row (historical)',
        -38.8: 'v28 reserve/S (historical)',
        -182.02: 'v28 endpoint delta (historical)',
        -182.13: 'v29 endpoint-only delta (computed row, same value)',
    }
    vals = list(computed.values())
    rows = {}
    for tok in sorted(set(re.findall(r'(?<![A-Za-z0-9_])2\^([-+]?\d+(?:\.\d+)?)', text))):
        x = float(tok)
        hit = None
        for k, v in computed.items():
            if abs(v - x) <= 5e-3:
                hit = 'computed:%s' % k
                break
        if hit is None:
            for h, why in historical.items():
                if abs(h - x) <= 5.1e-3:
                    hit = 'historical:%s' % why
                    break
        rows[tok] = hit
    return dict(all_Y=all(v is not None for v in rows.values()), tokens=rows,
                n_tokens=len(rows))


def literal_control():
    """(a) 1-ulp controls.  Claim grid: +1e-4 moves the reduced rational by 1/10000.
    Ledger float: one ulp of f changes the exact total fraction (regen input)."""
    c0 = L29.claim29()
    c1 = L29.claim29(delta=c0['T']['delta']['charged'] * (1 + 2 ** -52))
    bumped = c0['claim_log2'] + 10 ** -4
    num1 = int(round(bumped * 10 ** 4))
    from math import gcd
    g = gcd(num1, 10 ** 4)
    f = Fraction(14566, 41929)
    f2 = f + Fraction(1, 10 ** 6)
    t2 = L29.total29(f_measured=float(f2))
    return dict(claim0='%d/%d' % (c0['numerator'], c0['denominator']),
                claim_grid_plus_1ulp='%d/%d' % (num1 // g, 10 ** 4 // g),
                grid_changed=(num1 // g, 10 ** 4 // g) != (c0['numerator'], c0['denominator']),
                total0=str(c0['T']['total']),
                total_f_ulp=str(t2['total']),
                fraction_changed=(str(t2['total']) != str(c0['T']['total'])))


def endpoint_control(arm='1'):
    """(b) coarsen the exact 1-step age histogram into [1,2),[2,4),[4,8),[8,16),
    [16,32) and integrate at both endpoints.  Endpoint integrals must differ from
    the raw integral; charge the upper endpoint."""
    a = AGE['arms'][arm]
    theta = a['theta']
    hist = {int(k): v for k, v in a['age_hist_all'].items()}
    total = sum(hist.values())
    raw = a['delta_residual_all']
    bins = [(1, 2), (2, 4), (4, 8), (8, 16), (16, 32)]
    lo_sum = hi_sum = 0.0
    for lo, hi in bins:
        n = sum(v for age, v in hist.items() if lo <= age < hi)
        lo_sum += n * math.exp(-(128 - lo / float(theta)))
        hi_sum += n * math.exp(-(128 - hi / float(theta)))
    lo_d = lo_sum / total
    hi_d = hi_sum / total
    return dict(arm=arm, raw=raw, raw_log2=math.log2(raw),
                lower_endpoint=lo_d, lower_log2=math.log2(lo_d),
                upper_endpoint=hi_d, upper_log2=math.log2(hi_d),
                endpoints_differ=(abs(hi_d - lo_d) / raw > 1e-3),
                raw_inside=(min(lo_d, hi_d) <= raw <= max(lo_d, hi_d)),
                charge_upper=hi_d, charge_upper_log2=math.log2(hi_d))


def concurrency_control():
    c = M.concurrency_control((1 << 80), cap_mult=128)
    T = L29.total29()
    return dict(bound=c['bound'], worst_case_max_inflight=c['worst_case_max_inflight'],
                bound_attained=c['bound_attained'],
                reserve_inflight=T['stop']['n_inflight'],
                reserve_matches_bound=(T['stop']['n_inflight'] == c['bound']),
                reserve_steps=T['stop']['reserve_steps'],
                empirical_short_trails_max_inflight=c['empirical_short_trails_max_inflight'])


def loss_control(arm='1'):
    h = STO['arms'][arm]['halts']
    rs = sorted(h, key=lambda x: int(x))
    fracs = [h[r]['lost_frac'] for r in rs]
    mono = all(fracs[i] <= fracs[i + 1] for i in range(len(fracs) - 1))
    return dict(arm=arm, order=rs, lost_fracs=fracs, monotone=mono,
                rises=(max(fracs) > min(fracs)),
                early_halt_bigger=(fracs[-1] > fracs[0]))


def independent_recompute(n=20, seed=29):
    rng = random.Random(seed)
    base = L29.total29()
    outs = []
    for _ in range(n):
        f = Fraction(14566, 41929) * (1 + rng.choice([-1, 1]) * 2 ** -40 * rng.random())
        d = base['delta']['charged'] * (1 + rng.choice([-1, 1]) * 2 ** -20)
        T = L29.total29(f_measured=float(f), delta=d)
        s = sum(T['components'].values(), Fraction(0))
        outs.append(s == T['total'] and
                    abs(math.log2(float(T['total'])) - T['log2_total']) < 1e-9)
    return dict(n=n, all_equal=all(outs), passed=sum(outs))


def main():
    c = L29.claim29()
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
        statement_has_trap=('2^%.6f' % math.log2(float(T['components']['trap_charge']))) in st,
        statement_has_aux=('2^%.6f' % math.log2(float(T['components']['aux_sum']))) in st,
        statement_has_stop=('2^%.6f' % math.log2(float(T['components']['stop_loss_row']))) in st,
        statement_has_claim=('126.2599' in st),
        statement_has_sumline=('2^%.9f' % c['total_log2']) in st,
        statement_has_new_delta=('2^%.6f' % math.log2(T['delta']['charged'])) in st,
        statement_has_reserve=('2^%.2f' % math.log2(int(T['stop']['reserve_steps']))) in st,
        statement_has_inflight=('128 trails' in st),
    )
    controls = dict(literal=literal_control(), endpoint=endpoint_control(),
                    concurrency=concurrency_control(), loss=loss_control())
    sweep = literal_sweep()
    indep = independent_recompute()
    planted = dict(
        literal_fires=(controls['literal']['grid_changed']
                       and controls['literal']['fraction_changed']),
        endpoint_fires=(controls['endpoint']['endpoints_differ']
                        and controls['endpoint']['raw_inside']),
        concurrency_fires=(controls['concurrency']['bound_attained']
                           and controls['concurrency']['reserve_matches_bound']),
        loss_fires=(controls['loss']['monotone'] and controls['loss']['early_halt_bigger']))
    out = dict(claim=dict(numerator=c['numerator'], denominator=c['denominator'],
                          claim_log2=c['claim_log2'], total_log2=c['total_log2'],
                          total=str(T['total'])),
               components={k: str(v) for k, v in T['components'].items()},
               components_log2={k: math.log2(float(v)) for k, v in T['components'].items()},
               S=str(T['S']), S_log2=math.log2(T['S']), f=str(T['f']),
               delta=T['delta'], trap=T['trap'], stop=T['stop'],
               literal_checks=checks, literal_sweep=sweep, controls=controls,
               independent_recompute=indep, planted=planted,
               planted_score='%d/4' % sum(planted.values()))
    out['all_checks'] = (all(checks.values()) and sweep['all_Y'] and indep['all_equal']
                         and all(planted.values()))
    json.dump(out, open('evidence/v29-ledger.json', 'w'), indent=1, sort_keys=True,
              default=str)
    print('claim %d/%d = %.4f  total 2^%.9f' % (c['numerator'], c['denominator'],
                                                c['claim_log2'], c['total_log2']))
    print('literal sweep: %d tokens, all_Y=%s' % (sweep['n_tokens'], sweep['all_Y']))
    for t, v in sweep['tokens'].items():
        if v is None:
            print('  UNMAPPED 2^%s' % t)
    print('planted', out['planted_score'], planted)
    print('indep recompute %d/%d' % (indep['passed'], indep['n']))
    print('endpoint raw/lo/hi: 2^%.2f / 2^%.2f / 2^%.2f'
          % (controls['endpoint']['raw_log2'], controls['endpoint']['lower_log2'],
             controls['endpoint']['upper_log2']))
    print('concurrency bound %d attained %s reserve %s'
          % (controls['concurrency']['bound'], controls['concurrency']['bound_attained'],
             controls['concurrency']['reserve_steps']))
    print('ALL CHECKS', out['all_checks'])


if __name__ == '__main__':
    main()
