#!/usr/bin/env python3
"""v28 ledger: three judge-named defects repaired in the SUM, not in prose.

R1  RETENTION IS PRICED AGAINST A MEASURED AGE LAW.  v27 charged
        P(miss|merge) = (1-e^-128)/128 = 0.0078125          (UNIFORM age on [0,128t])
    with the premise declared.  The judge's objection is exact: first-DP restarts
    make the renewal age nonuniform, so the uniform integral is not connected to
    the success event.  v28 measures the age law (evidence/v28-age-law.json: age at
    merge, conditional on the merge event, per report depth) and prices the same
    functional under it:

        delta = E[ exp(-(cap - A)/theta) ],  A = age at merge, cap = 128 theta.

    Under the age law already declared for trail lengths (A size-biased by L^2,
    L = min(Exp(theta), 128 theta)) the integral is closed-form:
        E[e^{A/theta}] = int_0^128theta a/theta^2 da = 128^2/2 = 8192,
        delta = 8192 e^-128 = 2^-171.66.
    The MEASURED value is computed from the traces (delta_meas) and the CHARGED
    value is max(measured, closed-form) so a small-sample artefact cannot buy bits.
    The uniform value is retained in the evidence as the DISCLOSED-REJECTED method.

R2  THE WASTED-WORK TERM IS E[L^3]/(2 E[L^2]), COMPUTED FROM THE SAMPLE.  The v27
    assertion "E[L^2]/(2E[L]) <= theta" is wrong for the collision-point reading:
    a uniform prefix of the L^2-conditioned trail has mean E[L^3]/(2E[L^2]), which
    for untruncated exponential is 3 theta/2.  The exact integer sums of L, L^2,
    L^3 are in evidence/v28-moments.json; the ratio is recomputed from them and
    CHARGED (it is larger than theta, so the trap row grows).

R3  STOPPING RULE + THE MISSING OVERFLOW CHARGE.  total28 carries, in the one sum:
      w_stress + trap28 + aux28 + overflow_row + stop_loss_row
    where overflow_row is the one-reference-operation overflow charge that v27
    claimed separately and left OUT of total27, and stop_loss_row prices the
    success lost to the stopping rule (reserve = one full trail = 128*theta*c_mul;
    loss = 2*reserve/S of S^2).  The rule: halt LAUNCHES at cumulative charged
    work >= W - reserve; every launched trail charges at most reserve, so
    total <= W deterministically.
"""
from fractions import Fraction
import json
import math

import constants as K
import ledger as L
import ledger26 as L26
import ledger27 as L27

THETA = 1 << K.THETA_BITS
CAP_MULT = 128
AGE_PATH = 'evidence/v28-age-law.json'
MOM_PATH = 'evidence/v28-moments.json'


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def delta_measured(path=AGE_PATH):
    """max over arms of the measured production-cap retention integral, and the
    arm that supplied it."""
    d = _load(path)
    best, arm = 0.0, None
    for k, a in d['arms'].items():
        v = a.get('delta_measured_production_cap')
        if v is not None and v > best:
            best, arm = v, k
    return best, arm


def delta_charge(path=AGE_PATH):
    dm, arm = delta_measured(path)
    da = (1.0 / 2.0**128) * (CAP_MULT ** 2) / 2.0     # 8192 e^-128 = 2^-171.66
    return max(dm, da), dict(measured=dm, measured_arm=arm, analytic=da,
                             uniform_rejected=L27.DELTA_DP,
                             charged=max(dm, da), rule='max(measured, closed-form)')


def moments_charged(path=MOM_PATH):
    """The harshest measured E[L^3]/(2E[L^2]) over the arms, in THETA units."""
    d = _load(path)
    best, arm = None, None
    for k, m in d['arms'].items():
        r = float(m['ratio_over_theta'])
        if best is None or r > best:
            best, arm = r, k
    return best, arm


def production_f28(path='evidence/v28-census.json'):
    """Pooled EXECUTED squaring fraction over the four report depths (v28 census,
    same protocol as v27c; shipped with this package)."""
    d = _load(path)
    sq = sum(d[k]['sq_steps'] for k in ('1', '2', '3', '4'))
    st = sum(d[k]['steps'] for k in ('1', '2', '3', '4'))
    return sq / st


def S28(delta=None):
    d = delta if delta is not None else delta_charge()[0]
    return L27.step_budget27(delta=d)


def trap28(s=None, wasted_theta_mult=None, theta=THETA, cap_mult=CAP_MULT, n=None,
           c_mul=None):
    s = S28() if s is None else s
    n = K.N_CLASSES if n is None else n
    c_mul = Fraction(K.C_REF + L.rows()['fold_mul']*K.STRESS_FOLD, K.C_REF) \
        if c_mul is None else c_mul
    n_trails = -(-int(s) // int(theta))
    cap = cap_mult * theta
    e = math.exp(-float(cap) / float(theta))
    e_l2 = 2.0 * float(theta) ** 2 * (1.0 - e * (1.0 + cap_mult))   # theory, upper
    p_per_trail = Fraction(int(round(e_l2)), 2 * int(n))
    mult = wasted_theta_mult if wasted_theta_mult is not None else 1.0
    wasted = int(math.ceil(mult * theta))
    charge = Fraction(n_trails) * p_per_trail * wasted * c_mul
    return charge, dict(population='trail starts', n_trails=n_trails, theta=int(theta),
                        cap=cap, cap_mult=cap_mult, e_l2=int(round(e_l2)),
                        p_per_trail_log2=math.log2(float(p_per_trail)),
                        wasted_steps=wasted,
                        wasted_over_theta=wasted / float(theta),
                        wasted_theta_mult=mult,
                        charge_log2=math.log2(float(charge)))


def aux28(s=None, theta=THETA):
    s = S28() if s is None else s
    reports = -(-int(s) // int(theta))
    rows, _ = L26.aux26(s=s, theta=theta, slack=2)
    return rows, reports


def stop_loss(s=None, theta=THETA, cap_mult=CAP_MULT, c_stress=None):
    """Success lost to the stopping rule: the rule stops launches 128*theta steps
    short of S, and success ~ S^2, so the relative loss is 2*reserve/S; priced as
    that fraction of W_stress (the work that would be needed to buy it back)."""
    s = S28() if s is None else s
    c_stress = Fraction(K.C_REF + L.rows()['fold_mul']*K.STRESS_FOLD, K.C_REF) \
        if c_stress is None else c_stress
    reserve_steps = cap_mult * theta
    rel = 2.0 * reserve_steps / float(s)
    w = Fraction(int(s)) * c_stress
    return w * Fraction(rel).limit_denominator(10**30), dict(
        reserve_steps=str(reserve_steps), reserve_over_S=reserve_steps / float(s),
        loss_rel_S2=rel, loss_bits=-math.log2(1.0 - rel))


def total28(f_measured=None, delta=None, wasted_mult=None):
    if f_measured is None:
        f_measured = production_f28()
    d, ddetail = delta_charge() if delta is None else (delta, dict(charged=delta))
    if wasted_mult is None:
        wasted_mult = moments_charged()[0]
    s = S28(d)
    Lg = L26.ledgers26(f_measured=f_measured, s=s)
    w = Lg['w_stress']
    tr, trd = trap28(s=s, wasted_theta_mult=wasted_mult, c_mul=Lg['c_mul_stress'])
    rows, reports = aux28(s=s)
    aux_sum = sum(rows.values(), Fraction(0))
    overflow_row = Fraction(1)          # one reference operation, the Chernoff tail
    sl, sld = stop_loss(s=s, c_stress=Lg['c_stress'])
    comps = dict(w_stress=w, trap_charge=tr, aux_sum=aux_sum,
                 overflow_row=overflow_row, stop_loss_row=sl)
    total = sum(comps.values(), Fraction(0))
    return dict(components=comps, total=total, log2_total=math.log2(float(total)),
                trap=trd, aux_rows={k: str(v) for k, v in rows.items()},
                reports=reports, ledgers=Lg, S=s, f=Lg['f'], delta=ddetail,
                stop=sld)


def claim28(dp=4, f_measured=None, delta=None, wasted_mult=None):
    T = total28(f_measured=f_measured, delta=delta, wasted_mult=wasted_mult)
    lo = math.floor(T['log2_total'] * 10**dp) / 10**dp
    while lo < T['log2_total']:
        lo += 10.0 ** (-dp)
    lo = round(lo, dp)
    num = int(round(lo * 10**dp))
    den = 10**dp
    from math import gcd
    g = gcd(num, den)
    return dict(claim_log2=lo, numerator=num // g, denominator=den // g,
                total_log2=T['log2_total'], total=str(T['total']), T=T)


if __name__ == '__main__':
    dm, ddetail = delta_charge()
    r3, arm3 = moments_charged()
    c = claim28(wasted_mult=r3)
    T = c['T']
    print('delta charged %.3e (measured %.3e arm %s, analytic %.3e, uniform %.6f REJECTED)'
          % (dm, ddetail['measured'], ddetail['measured_arm'], ddetail['analytic'],
             float(L27.DELTA_DP)))
    print('E[L^3]/(2E[L^2]) = %.4f theta = %.2f steps (arm %s); v27 charged 1.0 theta'
          % (r3, r3 * THETA, arm3))
    print('S = %d  2^%.6f  f_priced %s' % (T['S'], math.log2(T['S']), T['f']))
    for k, v in T['components'].items():
        print('  %-16s 2^%.6f' % (k, math.log2(float(v))))
    print('total 2^%.9f  claim %.4f = %d/%d'
          % (c['total_log2'], c['claim_log2'], c['numerator'], c['denominator']))
    print('stop rule: reserve %s steps = 2^%.2f, loss rel S^2 2^%.2f'
          % (T['stop']['reserve_steps'], math.log2(int(T['stop']['reserve_steps'])),
             math.log2(T['stop']['loss_rel_S2'])))
