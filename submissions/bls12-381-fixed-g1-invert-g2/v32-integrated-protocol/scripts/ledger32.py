#!/usr/bin/env python3
"""v32 ledger.  Three repairs against v31, one row per verdict defect.

(a) THE OUTPUT PATH IS Q = [x]Q0, priced by instruction count (verdict defect 3).
    The GT discrete-log scalar x is what the walk returns.  Q0 is any fixed G2
    point; z0 = e(P,Q0) is computed once per instance, and then
        z = z0^x  =>  e(P, [x]Q0) = e(P,Q0)^x = z,
    so the returned point is [x]Q0.  There is NO affine square-root chain and no
    "unchanged return algebra": the row is
        Miller(e(P,Q0))            1533 Fp12 mul  (64 DBL x 9 + 87 ADD x 11)
        final exp (easy+hard)       300 Fp12 mul  (1 inv + ~10 mul + 3x64-bit exp)
        G2 reconstruction [x]Q0     226 Fp12 mul  (254 DBL x 24 + 127 ADD x 48 Fp mul,
                                                  /54 Fp mul per Fp12 mul)
        verify e(P,Q) = z           1833 Fp12 mul  (both sides variable)
    = 3892 Fp12 mul = 3892 * C_REF instruction equivalents.

(b) THE LABELS ARE 255-BIT, NOT 80-BIT/TWO-WORD (verdict defect 3).  step_v26
    updates (a,b) modulo r, so each label is four u64 limbs; the add is the
    5-limb add-then-conditional-subtract chain priced as 4 FP_ADD + 1 compare.

(c) THE STEP COUNTER IS 87 BITS (verdict defect 3): the trail cap is 128*theta =
    2^87 at production, so the own-step counter is 16 B (two u64), not 8 B.

The claim is emitted by regen_v32.py from this ledger; no claim literal is typed.
"""
import math
from fractions import Fraction

import constants as K
import ledger as L
import ledger26 as L26
import ledger27 as L27

THETA = 1 << K.THETA_BITS
CAP_MULT = 128
CAP = CAP_MULT * THETA

# ---- (a) output construction, instruction-counted ---------------------------
MEASURED_RETENTION = Fraction(3, 5)   # evidence/retention32.json: 6 of 10 reports shared

FP12_FROM_KARATSUBA = 54          # Fp12 mul = 18 Fp2 mul = 54 Fp mul
MILLER_DBL, MILLER_ADD = 9, 11    # Fp12 mul-equivalents per Miller DBL / ADD
MILLER_ITERS, MILLER_ADDS = 64, 87
FINAL_EXP = 300                   # Fp12 mul-equivalents
G2_DBL, G2_ADD = 24, 48           # Fp mul per Jacobian G2 DBL / mixed-ADD
SCALAR_BITS = K.R_ORDER.bit_length()   # 255


def output_rows():
    miller = MILLER_ITERS * MILLER_DBL + MILLER_ADDS * MILLER_ADD
    g2 = (SCALAR_BITS - 1) * G2_DBL + (SCALAR_BITS // 2) * G2_ADD
    g2_fp12 = g2 / float(FP12_FROM_KARATSUBA)
    fixed_p = miller + FINAL_EXP            # z0 = e(P,Q0)
    verify = miller + FINAL_EXP             # e(P,Q) = z
    total = fixed_p + FINAL_EXP * 0 + g2_fp12 + verify
    instr = total * K.C_REF
    return dict(miller_fp12=miller, final_exp_fp12=FINAL_EXP,
                g2_fp_mul=g2, g2_fp12mul=g2_fp12,
                fixed_p_fp12=fixed_p, verify_fp12=verify,
                pairing_fp12mul=total, pairing_instructions=instr,
                pairing_log2=math.log2(instr),
                pairing_row=Fraction(int(round(instr)), K.C_REF))


def width_rows():
    """(b)/(c): the widths the shipped step actually touches."""
    return dict(
        label_bits=SCALAR_BITS, label_limbs=4, label_bytes=32,
        label_add_limbs=5,                       # 4-limb add + carry limb
        label_add_cost=5 * K.FP_ADD + K.LOOP,
        step_counter_bytes=16,                   # 2^87 cap needs 87 bits
        step_counter_bits=87,
        key_bytes=48, key_words=6,
        record_bytes=K.RECORD_BYTES,
        live_gt_states=CAP_MULT, live_bytes=CAP_MULT * 656,
        v31_label_claim='80-bit / two-word', v32_label='255-bit / four u64 limbs',
        v31_counter_claim='8 B', v32_counter='16 B')


def rows32(delta=None, retention=None):
    """delta = P(miss | merge); retention = the MEASURED first-report retention.

    v31 charged max(fixture estimate, 1/64).  v32 charges the measured value in
    the slot the sizing law actually divides by (v27's `retention`), so the
    fixture number enters the claim instead of a declared floor."""
    d = Fraction(1, 128) if delta is None else Fraction(delta).limit_denominator(10 ** 6)
    if retention is None:
        retention = MEASURED_RETENTION
    ret = Fraction(retention).limit_denominator(10 ** 6)
    S = L27.step_budget27(delta=d, retention=ret)
    Lg = L26.ledgers26(s=S)
    c_max = Lg['c_mul_stress']
    w_eff = S * c_max
    launches = -(-S // THETA)
    e = math.exp(-CAP / float(THETA))
    e_l2 = 2.0 * float(THETA) ** 2 * (1.0 - e * (1.0 + CAP_MULT))
    p_trail = Fraction(int(round(e_l2)), 2 * int(K.N_CLASSES))
    abandon = Fraction(launches) * p_trail * THETA * c_max
    in_flight_row = CAP_MULT * CAP * c_max
    aux_rows, reports = L27.aux27(s=S)
    aux_sum = sum(aux_rows.values(), Fraction(0))
    orow = output_rows()
    return dict(S=S, theta=THETA, cap=CAP, launches=launches, delta=d,
                retention=ret,
                c_max=c_max, w_effective=w_eff, c_mean=Lg['c_stress'],
                abandon_row=abandon, in_flight_row=in_flight_row,
                aux_sum=aux_sum, reports=reports,
                pairing_row=orow['pairing_row'],
                overflow_row=Fraction(0), stop_loss_row=Fraction(0),
                pairing_detail=orow, width=width_rows(), aux_rows=aux_rows,
                report_mu=K.MU_WALK)


def total32(delta=None, retention=None):
    r = rows32(delta, retention)
    # stop_loss: the cumulative abandoned-work row plus the in-flight row are the
    # pre-halt population that a stop-loss schedule must be able to pay for.
    r['stop_loss_row'] = r['abandon_row'] + r['in_flight_row']
    comps = dict(w_effective=r['w_effective'], abandon_row=r['abandon_row'],
                 aux_sum=r['aux_sum'], pairing_row=r['pairing_row'],
                 overflow_row=r['overflow_row'], stop_loss_row=r['stop_loss_row'])
    tot = sum(comps.values(), Fraction(0))
    return dict(components=comps, total=tot, log2_total=math.log2(float(tot)),
                rows=r)


def claim32(delta=None, retention=None):
    T = total32(delta, retention)
    lg = T['log2_total']
    lo = math.floor(lg * 10 ** 4) / 10.0 ** 4
    while lo < lg:
        lo += 1e-4
    lo = round(lo, 4)
    num, den = int(round(lo * 10 ** 4)), 10 ** 4
    g = math.gcd(num, den)
    return dict(claim_log2=lo, numerator=num // g, denominator=den // g,
                total_log2=lg, total=str(T['total']), T=T)


if __name__ == '__main__':
    c = claim32()
    T = c['T']
    print('S=%d delta=%s c_max=%.6f' % (T['rows']['S'], T['rows']['delta'],
                                        float(T['rows']['c_max'])))
    for k, v in T['components'].items():
        print('  %-14s 2^%.6f' % (k, math.log2(float(v)) if v else -math.inf))
    print('pairing', T['rows']['pairing_detail'])
    print('TOTAL 2^%.6f  CLAIM %s = %d/%d'
          % (c['total_log2'], c['claim_log2'], c['numerator'], c['denominator']))
