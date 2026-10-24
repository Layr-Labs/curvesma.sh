#!/usr/bin/env python3
"""v29 ledger.  Four judge-named defects repaired in the arithmetic, not in prose.

FIX 1  LITERAL/CONVENTION.  delta_charge now computes the quantity the text has
       always described, with the L^2 normalisation and the cap atom written out:

           L ~ min(Exp(theta), 128 theta), collision-weighted by L^2
           E[L^2]        = 2 theta^2 (1 - e^-128 (1+128))
           E[e^(A/theta)] = (128^3/3 + 128^2) / (2 (1 - 129 e^-128))
           delta_closed   = e^-128 * E[e^(A/theta)] = 2^-166.216...

       The v28 code computed 8192*2^-128 = 2^-115 while the text said
       8192*e^-128 = 2^-171.66: the code's 2^-115 was a typo-class literal, and
       the text's 2^-171.66 also dropped the normalisation.  Both are replaced by
       ONE computed value, and regen_claim.py emits every literal in the claim
       from this module's output (no prose-typed exponent).

FIX 2  AGE AT THE MERGE.  The charged measured value comes from v29_measure's
       merge-age law (every visited class key stored, every repeat counted),
       not from endpoint report ages.  Both are emitted; the raw-merge residual
       integral is the one charged against the closed form.

FIX 3  CONCURRENCY.  Launch schedule: one trail start per theta executed steps
       (the schedule ledger27 already declares, measured in the census as
       trail_starts ~ S/theta).  A trail lasts at most cap = 128 theta steps, so
       in-flight trails <= ceil(cap/theta) = 128.  reserve = 128 * 128 theta
       = 2^94 step-equivalents.  The rule is stated on COMPLETED STEPS, not on
       realized work at an averaged price: halt launches when completed steps
       >= S - reserve.  Every launched trail then finishes inside the reserve.

FIX 4  STOPPING LOSS.  The loss is the arithmetic share 2*reserve/S of success
       (success ~ S^2), charged as stop_loss_row = W_stress * 2*reserve/S;
       it is NOT asserted to be zero.  The unit-exact stop row from the simulator
       (v29-stopping.json) is shipped as the measured response curve, and the
       audit's control (d) forces an early halt and sees the measured loss rise.
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
N_INFLIGHT = 128                     # ceil(cap/theta) under one launch per theta steps
AGE_PATH = 'evidence/v29-age-law.json'
MOM_PATH = 'evidence/v29-moments.json'


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def delta_measured(path=AGE_PATH):
    """Harshest measured residual integral over the RAW merge ages, and the arm."""
    d = _load(path)
    best, arm = 0.0, None
    for k, a in d['arms'].items():
        v = a.get('delta_residual_all')
        if v is not None and v > best:
            best, arm = v, k
    return best, arm


def delta_closed_form(cap_mult=CAP_MULT):
    """e^-m * (m^3/3 + m^2) / (2 (1 - e^-m (1+m))), the normalised L^2-size-biased
    law WITH the cap atom.  (v28's 8192 e^-128 dropped both the normalisation and
    the atom; it is intentionally not reproduced here.)"""
    m = float(cap_mult)
    e = math.exp(-m)
    e_l2 = 2.0 * (1.0 - e * (1.0 + m))          # in theta^2 units
    return e * (m ** 3 / 3.0 + m ** 2) / e_l2


def delta_charge(path=AGE_PATH):
    dm, arm = delta_measured(path)
    da = delta_closed_form()
    return max(dm, da), dict(measured=dm, measured_arm=arm, analytic=da,
                             analytic_log2=math.log2(da),
                             measured_log2=math.log2(dm) if dm > 0 else None,
                             charged=max(dm, da), rule='max(measured raw-merge, '
                             'normalised closed form with cap atom)',
                             uniform_rejected=L27.DELTA_DP,
                             v28_text_value=math.exp(-128.0) * 8192.0)


def moments_charged(path=MOM_PATH):
    d = _load(path)
    best, arm = None, None
    for k, m in d['arms'].items():
        r = float(m['ratio_over_theta'])
        if best is None or r > best:
            best, arm = r, k
    return best, arm


def production_f29(path='evidence/v29-census.json'):
    d = _load(path)
    sq = sum(d[k]['sq_steps'] for k in ('1', '2', '3', '4'))
    st = sum(d[k]['steps'] for k in ('1', '2', '3', '4'))
    return sq / st


def S29(delta=None):
    d = delta if delta is not None else delta_charge()[0]
    return L27.step_budget27(delta=d)


def trap29(s=None, wasted_theta_mult=None, theta=THETA, cap_mult=CAP_MULT, n=None,
           c_mul=None):
    s = S29() if s is None else s
    n = K.N_CLASSES if n is None else n
    c_mul = Fraction(K.C_REF + L.rows()['fold_mul'] * K.STRESS_FOLD, K.C_REF) \
        if c_mul is None else c_mul
    n_trails = -(-int(s) // int(theta))
    cap = cap_mult * theta
    e = math.exp(-float(cap) / float(theta))
    e_l2 = 2.0 * float(theta) ** 2 * (1.0 - e * (1.0 + cap_mult))
    p_per_trail = Fraction(int(round(e_l2)), 2 * int(n))
    mult = wasted_theta_mult if wasted_theta_mult is not None else 1.0
    wasted = int(math.ceil(mult * theta))
    charge = Fraction(n_trails) * p_per_trail * wasted * c_mul
    return charge, dict(population='trail starts', n_trails=n_trails, theta=int(theta),
                        cap=cap, cap_mult=cap_mult, e_l2=int(round(e_l2)),
                        p_per_trail_log2=math.log2(float(p_per_trail)),
                        wasted_steps=wasted, wasted_over_theta=wasted / float(theta),
                        wasted_theta_mult=mult, charge_log2=math.log2(float(charge)))


def aux29(s=None, theta=THETA):
    s = S29() if s is None else s
    reports = -(-int(s) // int(theta))
    rows, _ = L26.aux26(s=s, theta=theta, slack=2)
    return rows, reports


def stop_loss29(s=None, theta=THETA, cap_mult=CAP_MULT, c_stress=None,
                n_inflight=N_INFLIGHT):
    """reserve = N_inflight * cap, in STEPS.  Loss = 2*reserve/S of success priced
    at W_stress.  N_inflight is derived, not asserted: launch interval = theta
    steps, trail ceiling = cap steps => at most cap/theta trails in flight."""
    s = S29() if s is None else s
    c_stress = Fraction(K.C_REF + L.rows()['fold_mul'] * K.STRESS_FOLD, K.C_REF) \
        if c_stress is None else c_stress
    cap = cap_mult * theta
    n_inflight = min(int(n_inflight), -(-cap // theta))
    reserve_steps = n_inflight * cap
    rel = 2.0 * reserve_steps / float(s)
    w = Fraction(int(s)) * c_stress
    return w * Fraction(rel).limit_denominator(10 ** 30), dict(
        n_inflight=n_inflight, launch_interval_steps=theta, trail_cap_steps=cap,
        reserve_steps=str(reserve_steps), reserve_over_S=reserve_steps / float(s),
        loss_rel_S2=rel, loss_bits=-math.log2(1.0 - rel),
        rule='halt launches when completed steps >= S - reserve')


def total29(f_measured=None, delta=None, wasted_mult=None):
    if f_measured is None:
        f_measured = production_f29()
    d, ddetail = delta_charge() if delta is None else (delta, dict(charged=delta))
    if wasted_mult is None:
        wasted_mult = moments_charged()[0]
    s = S29(d)
    Lg = L26.ledgers26(f_measured=f_measured, s=s)
    w = Lg['w_stress']
    tr, trd = trap29(s=s, wasted_theta_mult=wasted_mult, c_mul=Lg['c_mul_stress'])
    rows, reports = aux29(s=s)
    aux_sum = sum(rows.values(), Fraction(0))
    overflow_row = Fraction(1)
    sl, sld = stop_loss29(s=s, c_stress=Lg['c_stress'])
    comps = dict(w_stress=w, trap_charge=tr, aux_sum=aux_sum,
                 overflow_row=overflow_row, stop_loss_row=sl)
    total = sum(comps.values(), Fraction(0))
    return dict(components=comps, total=total, log2_total=math.log2(float(total)),
                trap=trd, aux_rows={k: str(v) for k, v in rows.items()},
                reports=reports, ledgers=Lg, S=s, f=Lg['f'], delta=ddetail,
                stop=sld)


def claim29(dp=4, f_measured=None, delta=None, wasted_mult=None):
    T = total29(f_measured=f_measured, delta=delta, wasted_mult=wasted_mult)
    lo = math.floor(T['log2_total'] * 10 ** dp) / 10 ** dp
    while lo < T['log2_total']:
        lo += 10.0 ** (-dp)
    lo = round(lo, dp)
    num = int(round(lo * 10 ** dp))
    den = 10 ** dp
    from math import gcd
    g = gcd(num, den)
    return dict(claim_log2=lo, numerator=num // g, denominator=den // g,
                total_log2=T['log2_total'], total=str(T['total']), T=T)


if __name__ == '__main__':
    dm, ddetail = delta_charge()
    r3, arm3 = moments_charged()
    c = claim29(wasted_mult=r3)
    T = c['T']
    print('delta charged %.3e = 2^%.6f (measured %.3e 2^%.2f arm %s; closed %.3e '
          '2^%.6f; v28 text 2^%.2f)'
          % (dm, math.log2(dm), ddetail['measured'], ddetail['measured_log2'] or 0,
             ddetail['measured_arm'], ddetail['analytic'], ddetail['analytic_log2'],
             math.log2(ddetail['v28_text_value'])))
    print('E[L^3]/(2E[L^2]) = %.4f theta (arm %s)' % (r3, arm3))
    print('S = %d  2^%.6f  f_priced %s' % (T['S'], math.log2(T['S']), T['f']))
    for k, v in T['components'].items():
        print('  %-16s 2^%.6f' % (k, math.log2(float(v))))
    print('total 2^%.9f  claim %.4f = %d/%d'
          % (c['total_log2'], c['claim_log2'], c['numerator'], c['denominator']))
    print('stop: N_inflight %d reserve %s steps 2^%.2f loss rel 2^%.2f'
          % (T['stop']['n_inflight'], T['stop']['reserve_steps'],
             math.log2(int(T['stop']['reserve_steps'])),
             math.log2(T['stop']['loss_rel_S2'])))
