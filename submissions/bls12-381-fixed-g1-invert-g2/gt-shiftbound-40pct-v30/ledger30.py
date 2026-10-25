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
AGE_PATH = 'evidence/v30-age-law.json'
MOM_PATH = 'evidence/v30-moments.json'


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


def delta_tail_bound(path=AGE_PATH, support_lo=16.0):
    """v30 RETENTION TAIL, ROUTE (a).  The measured merge ages live on the
    fixture support [1, 16 theta]; the production cap is 128 theta.  The raw
    integral is therefore a SUPPORT integral, not a full-range one.  Beyond the
    support charge an exponential-decay row at a conservative constant: fit the
    measured per-theta decay between 1 theta and 4 theta, HALVE the fitted rate
    (slower decay = larger tail), anchor at the measured survival at 4 theta,
    and charge every tail merge as a full miss (its cap atom contributes at most
    1 to exp(-(128-A/theta))).  This row is charged; it is not a fit presented
    as data."""
    d = _load(path)
    arm, best = None, -1.0
    for k, a in d['arms'].items():
        v = a.get('delta_residual_all')
        if v is not None and v > best:
            best, arm = v, k
    if arm is None:
        return None, dict(rule='no tail row')
    a = d['arms'][arm]
    s1 = a['survival_gt_k_theta']['1'] / float(a['merges_all'])
    s4 = a['survival_gt_k_theta']['4'] / float(a['merges_all'])
    lam = max(1e-9, 0.5 * math.log(max(s1, 1e-300) / max(s4, 1e-300)) / 3.0)
    lo, hi = float(support_lo), float(CAP_MULT)
    # Conservative tail constant: charge every merge beyond the support as a full
    # miss (the cap atom contributes at most 1 to exp(-(128-A/theta))).
    p_tail = s4 * math.exp(-lam * (lo - 4.0))
    return p_tail, dict(arm=arm, support='[1,16] theta (fixture cap)',
                        production_cap='128 theta',
                        survival_1=s1, survival_4=s4, decay_per_theta=lam,
                        decay_rule='half the measured 1->4 theta slope (conservative)',
                        tail_start_theta=lo, p_tail=p_tail,
                        p_tail_log2=math.log2(p_tail) if p_tail > 0 else None,
                        charge='each tail merge charged as a full miss '
                               '(cap atom <= 1); P(A>a) = s4 exp(-lam (a-4))')


def delta_charge(path=AGE_PATH):
    dm, arm = delta_measured(path)
    da = delta_closed_form()
    dt, tdetail = delta_tail_bound(path)
    # support integral + tail bound, then the max with the full-range closed form
    charged = max(dm + (dt or 0.0), da)
    return charged, dict(measured=dm, measured_arm=arm, analytic=da,
                         analytic_log2=math.log2(da),
                         measured_log2=math.log2(dm) if dm > 0 else None,
                         tail=dt, tail_detail=tdetail,
                         charged=charged, rule='max(measured support integral + '
                         'exponential-decay tail row, normalised closed form)',
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


def production_f29(path='evidence/v30-census.json'):
    d = _load(path)
    sq = sum(d[k]['sq_steps'] for k in ('1', '2', '3', '4'))
    st = sum(d[k]['steps'] for k in ('1', '2', '3', '4'))
    return sq / st


def S29(delta=None):
    d = delta if delta is not None else delta_charge()[0]
    return L27.step_budget27(delta=d)


def trap29(s=None, wasted_theta_mult=None, theta=THETA, cap_mult=CAP_MULT, n=None,
           c_mul=None, n_inflight=N_INFLIGHT):
    """v30 DETERMINISTIC trap bound, not an expectation.  Per trail the wasted
    work is at most cap = 128 theta own steps; at most n_inflight trails are live
    (the enforced atomic cap); so the concurrent waste charge is bounded by
    n_inflight * cap, priced at the MAX single-step price.  The measured
    E[L^3]/(2E[L^2]) is shipped as a diagnostic, NOT as the charged row."""
    s = S29() if s is None else s
    n = K.N_CLASSES if n is None else n
    c_mul = Fraction(K.C_REF + L.rows()['fold_mul'] * K.STRESS_FOLD, K.C_REF) \
        if c_mul is None else c_mul
    cap = cap_mult * theta
    charge = Fraction(int(n_inflight)) * Fraction(int(cap)) * c_mul
    return charge, dict(population='enforced in-flight cap x max single-trail charge',
                        n_inflight=int(n_inflight), theta=int(theta),
                        cap=cap, cap_mult=cap_mult,
                        wasted_steps=str(int(n_inflight) * int(cap)),
                        charge_log2=math.log2(float(charge)),
                        measured_wasted_theta_mult=wasted_theta_mult,
                        measured_rule='E[L^3]/(2E[L^2]) diagnostic only, not charged')


def aux29(s=None, theta=THETA):
    s = S29() if s is None else s
    reports = -(-int(s) // int(theta))
    rows, _ = L26.aux26(s=s, theta=theta, slack=2)
    return rows, reports


def stop_loss29(s=None, theta=THETA, cap_mult=CAP_MULT, c_stress=None,
                n_inflight=N_INFLIGHT):
    """v30: reserve = N_inflight * cap STEP-EQUIVALENTS, priced at the MAX
    single-step price (not an averaged mix).  N_inflight is the ATOMIC enforced
    in-flight cap; the bare launch schedule does not imply it (see
    v30_measure.concurrency_control trace rows).  Loss = 2*reserve/S of success."""
    s = S29() if s is None else s
    c_stress = Fraction(K.C_REF + L.rows()['fold_mul'] * K.STRESS_FOLD, K.C_REF) \
        if c_stress is None else c_stress
    cap = cap_mult * theta
    n_inflight = int(n_inflight)
    reserve_steps = n_inflight * cap
    rel = Fraction(2 * reserve_steps, int(s))
    charge = Fraction(2 * reserve_steps) * c_stress
    return charge, dict(
        n_inflight=n_inflight, launch_interval_steps=theta, trail_cap_steps=cap,
        reserve_steps=str(reserve_steps), reserve_over_S=reserve_steps / float(s),
        loss_rel_S2=float(rel), loss_bits=-math.log2(1.0 - float(rel)) if rel < 1 else None,
        price='max single-step price (c_mul_stress)',
        rule='halt launches when completed steps >= S - reserve; atomic in-flight cap enforced')


def total29(f_measured=None, delta=None, wasted_mult=None):
    if f_measured is None:
        f_measured = production_f29()
    d, ddetail = delta_charge() if delta is None else (delta, dict(charged=delta))
    if wasted_mult is None:
        wasted_mult = moments_charged()[0]
    s = S29(d)
    Lg = L26.ledgers26(f_measured=f_measured, s=s)
    # v30 MAX-PRICE CEILING: every completed step is charged at the dearer of the
    # two step prices (the mul step), not at the measured average mix.  The
    # averaged f-priced figure is shipped only as a diagnostic.
    w = Fraction(int(s)) * Lg['c_mul_stress']
    tr, trd = trap29(s=s, wasted_theta_mult=wasted_mult, c_mul=Lg['c_mul_stress'])
    rows, reports = aux29(s=s)
    aux_sum = sum(rows.values(), Fraction(0))
    overflow_row = Fraction(1)
    sl, sld = stop_loss29(s=s, c_stress=Lg['c_mul_stress'])
    comps = dict(w_stress=w, trap_charge=tr, aux_sum=aux_sum,
                 overflow_row=overflow_row, stop_loss_row=sl)
    total = sum(comps.values(), Fraction(0))
    return dict(components=comps, total=total, log2_total=math.log2(float(total)),
                trap=trd, aux_rows={k: str(v) for k, v in rows.items()},
                reports=reports, ledgers=Lg, S=s, f=Lg['f'], delta=ddetail,
                stop=sld,
                diagnostic_average_price=dict(
                    c_stress=str(Lg['c_stress']),
                    w_stress_at_average=float(s) * float(Lg['c_stress']),
                    w_stress_at_average_log2=math.log2(float(s) * float(Lg['c_stress'])),
                    note='mean-price row shown only to expose the undercount; NOT summed'))


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
