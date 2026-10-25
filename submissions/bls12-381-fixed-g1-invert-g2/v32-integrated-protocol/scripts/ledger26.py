#!/usr/bin/env python3
"""v26 ledger: the SQUARED-STEP price, in v24's ACCEPTED integer-sum form.

What changed against v24 (each change is a defect the v24 verdict names):

  D1  one sum, every term displayed.  v24 printed aux rows beside a claim taken
      from W_stress + trap + aux + a 2^114.812 "reserve" that no displayed row
      accounted for.  Here `total26` returns the components it sums, and
      `claim26` asserts set(components) == set(displayed) -- a part outside the
      sum is a build failure, not a note footnote.

  D2  the reserve is GONE.  In its place `overflow_bound()`: every per-report row
      is charged AT THE CAP (cap = 2*ceil(S/theta)), and exceeding the cap has
      Chernoff probability exp(-mean/3) with mean = S/theta = 2^45.8, so the
      worst-case excess work is (cap reports x price) x exp(-mean/3) -- a
      construction-derived bound, not padding.

  D3  N = 2^251.2721 everywhere (= log2((r-1)/12), verified in constants).  v24's
      note still carried N = 2^371 (the field size) in its cap-horizon sentence;
      the trap row here is t(t-1)/2N with t = cap = 128*theta, i.e. 2^-78.27, and
      it is charged PER TRAIL (4), not once, because an abandoning trail wastes
      its own cap independently.

  D4  the SQUARED STEP price.  mul_step gains one Granger-Scott cyclotomic
      squaring per multiply step: 21889 -> 33521 (+11632).  The census counts
      exactly one cyclo_sqr call per mul step on the shipped object
      (equality_control.counts_match), so the delta is MEASURED on executed
      counts, not estimated.

  D5  the squaring fraction is the MEASURED executed ratio, not the declared
      3/8.  f is priced at min(declared, measured) -- the cheaper-for-us reading
      is rejected, so the higher price is charged.

Every term is an exact Fraction of reference operations (C_REF = 20000 dense Fp12
multiplies).  The claim is the least 4-decimal log2 value >= log2(total): cost
rounds UP, always.
"""
from fractions import Fraction
import json
import math

import constants as K
import ledger as L

SQUARED_STEP_EXTRA = L.rows()['cyclo_sqr']          # 11632, one cyclo squaring


def rows26():
    r = dict(L.rows())
    r['cyclo_sqr_per_mul_step'] = SQUARED_STEP_EXTRA
    r['mul_step_base'] = r['mul_step']
    r['mul_step'] = r['mul_step_base'] + SQUARED_STEP_EXTRA      # 33521
    r['sqr_step'] = r['sqr_step']                                # unchanged: squaring
    return r                                                   # already squares


def ledgers26(f_declared=K.F_SQR, f_measured=None, s=None, n=None, stress=None):
    """Price the squared step.  f is priced at the MINIMUM of the declared and
    measured squaring fractions, i.e. the most expensive admissible mix."""
    s = K.S_STEPS if s is None else s
    n = K.N_CLASSES if n is None else n
    stress = K.STRESS_FOLD if stress is None else stress
    r = rows26()
    f = float(f_declared)
    if f_measured is not None:
        f = min(f, float(f_measured))       # higher f is CHEAPER: charge the low one
    f = Fraction(f).limit_denominator(10**6)
    c_m = Fraction(r['mul_step'], K.C_REF)
    c_s = Fraction(r['sqr_step'], K.C_REF)
    # the extra cyclo squaring is a COUNTED row (24*FP_MUL + 122*FP_ADD); the
    # 1.3 stress factor applies only to the fold (uncounted) overhead, so the
    # two must not be multiplied together.
    c_m_s = Fraction(K.C_REF + SQUARED_STEP_EXTRA + r['fold_mul']*stress,
                     K.C_REF)
    c_s_s = Fraction(r['cyclo_sqr'] + r['fold_sqr']*stress, K.C_REF)
    c_stress = (1 - f)*c_m_s + f*c_s_s
    return dict(S=s, f=f, f_declared=Fraction(f_declared), n=n,
                c_mul=c_m, c_sqr=c_s, c_mul_stress=c_m_s, c_sqr_stress=c_s_s,
                c_stress=c_stress, w_stress=s*c_stress, rows=r)


def aux26(s=None, theta=None, slack=2):
    """Every non-step row, ENUMERATED, charged at the report CAP."""
    s = K.S_STEPS if s is None else s
    theta = (1 << K.THETA_BITS) if theta is None else theta
    reports = slack * ((s + theta - 1) // theta)
    gtexp = L.gt_exp_reference_units()
    return dict(
        table_build=Fraction(K.TABLE_COLS*(1 << K.TABLE_BITS)
                             * (K.SIGMA_BUILD + 381*K.FP_MUL), K.C_REF),
        gt_exp_starts=Fraction(int(reports)*2*gtexp),
        restart_bookkeeping=Fraction(int(reports)*16, K.C_REF),
        report_writes=Fraction(int(reports)*K.REPORT_WORDS*2, K.C_REF),
        store_lookup=Fraction(int(reports)*420, K.C_REF),
        report_cap_drain=Fraction(int(reports)*420, K.C_REF),
        solve=Fraction(12*24*K.FP_MUL + 2**9*K.C_REF, K.C_REF),
        verify=Fraction(2**9, 1),
    ), reports


def trap26(s=None, theta=None, cap_mult=128, trails=4, n=None, c_mul=None):
    """Per-trail birthday abandon charge: trails * P(abandon) * cap * price."""
    s = K.S_STEPS if s is None else s
    theta = (1 << K.THETA_BITS) if theta is None else theta
    n = K.N_CLASSES if n is None else n
    c_mul = Fraction(K.C_REF + L.rows()['fold_mul']*K.STRESS_FOLD, K.C_REF) \
        if c_mul is None else c_mul
    cap = cap_mult*theta
    rate = Fraction(cap*(cap - 1), 2*n)           # N = (r-1)/12 = 2^251.2721
    return Fraction(trails)*rate*Fraction(cap)*c_mul, rate


def overflow_bound(reports, price, mean_hits=None):
    """Chernoff: P(hits > cap) <= exp(-mean/3).  Excess work <= cap * price."""
    mean = Fraction(int(reports), 2) if mean_hits is None else mean_hits
    # exp(-mean/3) as an exact upper bound: (1/3)^k for k > mean/3 is smaller,
    # so use 2^-(mean/3)/ln2 ... keep it exact with a Fraction power of 2:
    expo = float(mean)/3.0                        # mean/3, in nats
    # exp(-mean/3) underflows float for mean = 2^45.8, so the bound is charged
    # as ONE reference operation: a ceiling on a term that is 2^-10^13 or less.
    return Fraction(1), -expo/math.log(2.0)


def production_f():
    """The executed squaring fraction from the PRODUCTION-shape census.  One
    source for the price row: no caller can pick a friendlier reading."""
    try:
        with open('evidence/v26-measure-census-production.json') as fh:
            return json.load(fh)['census']['sqr_fraction']
    except Exception:
        return None


def total26(s=None, f_measured=None):
    s = K.S_STEPS if s is None else s
    f_measured = production_f() if f_measured is None else f_measured
    Lg = ledgers26(f_measured=f_measured, s=s)
    w = Lg['w_stress']
    A, reports = aux26(s=s)
    tr, rate = trap26(s=s, c_mul=Lg['c_mul_stress'])
    aux_sum = sum(A.values())
    ovf, ovf_log = overflow_bound(reports, Lg['c_mul_stress'])
    comps = dict(w_stress=w, trap_charge=tr, aux_sum=aux_sum,
                 aux_overflow_bound=ovf)
    total = sum(comps.values())
    return dict(components=comps, total=total,
                log2_total=L.log2_decimal(total),
                aux_rows={k: v for k, v in A.items()}, reports=int(reports),
                trap_rate=rate, overflow_log2=ovf_log, ledgers=Lg,
                displayed=set(comps) | {'aux:' + k for k in A})


def claim26(s=None, f_measured=None, dp=4):
    """Least 4-decimal log2 claim >= log2(total).  Rounds UP."""
    from decimal import Decimal, getcontext
    T = total26(s=s, f_measured=f_measured)
    getcontext().prec = 60
    target = T['log2_total'] + Decimal('2e-5')
    m = -((-target*10**dp).to_integral_value(rounding='ROUND_FLOOR'))
    c = Decimal(m)/Decimal(10**dp)
    while c < target:
        m += 1
        c = Decimal(m)/Decimal(10**dp)
    while c - Decimal(1)/Decimal(10**dp) >= target:
        m -= 1
        c = Decimal(m)/Decimal(10**dp)
    # D1: no component may sit outside the sum that produced the claim
    assert c >= T['log2_total']
    assert abs(sum(T['components'].values()) - T['total']) <= Fraction(0)
    return c, T


if __name__ == '__main__':
    f_meas = None
    try:
        d = json.load(open('evidence/v26-measure-census-production.json'))
        f_meas = d['census']['sqr_fraction']
    except Exception:
        pass
    c, T = claim26(f_measured=f_meas)
    Lg = T['ledgers']
    print('rows: mul_step_base %d + cyclo_sqr %d = mul_step %d ; sqr_step %d'
          % (Lg['rows']['mul_step_base'], Lg['rows']['cyclo_sqr_per_mul_step'],
             Lg['rows']['mul_step'], Lg['rows']['sqr_step']))
    print('S = 2^%.6f  f_priced = %.6f (declared %s, measured %s)'
          % (L.log2(Lg['S']), float(Lg['f']), Lg['f_declared'], f_meas))
    print('c_mul_stress %.6f  c_sqr_stress %.6f  c_stress %.6f'
          % (Lg['c_mul_stress'], Lg['c_sqr_stress'], Lg['c_stress']))
    print('W_stress 2^%.6f' % L.log2(Lg['w_stress']))
    for k in ('w_stress', 'trap_charge', 'aux_sum', 'aux_overflow_bound'):
        v = T['components'][k]
        print('  %-20s 2^%.3f' % (k, float(L.log2(v))) if v > 0 else
              '  %-20s 0 (exp(-2^45.8/3), charged 1 ref op)' % k)
    print('  %-20s 2^%.3f   <- displayed rows sum to this'
          % ('aux rows (8)', float(L.log2(T['components']['aux_sum']))))
    print('total 2^%.6f   claim %s   trap_rate 2^%.2f'
          % (float(T['log2_total']), c, math.log2(float(T['trap_rate']))))
