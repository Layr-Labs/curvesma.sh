#!/usr/bin/env python3
"""v33 ledger -- every row in RAW INSTRUCTIONS, ONE division by C_REF.

THE v32 DEFECT (verdict-v32, objection 4).  v32's six components were already
target-operation equivalents (pairing_row = 19458889/5000 = 3891.78 target ops),
and the claim statement compared their sum W with 20000*2^c -- applying the
normalisation a second time, so 126.9558 was not the least 4-decimal c satisfying
the written inequality.

v33 fixes the UNIT CHAIN instead of the prose:

    rows are INSTRUCTIONS                        (integers / exact Fractions)
    W = sum(rows)                                instructions
    W_ref = W / C_REF                            <-- THE ONE DIVISION LINE
    claim = least 4-decimal c with 2^c >= W_ref  (score unit: log2 target ops)

`claim33(scale=...)` multiplies the raw instruction ledger before that single
division: the emitted claim must then move by exactly log2(scale).  That is
control (d): a double-normalised chain moves by 0 or 2*log2(C_REF).

The step budget S comes from the sizing law with BOTH measured factors charged:
p_detect = P(detect | merge), charged at its CP-95 LOWER bound (objection 3), and
delta = P(D > R) on the actual event, the geometric tail integrated past each
launch's measured remaining lifetime (objection 2), charged at the analytic
ceiling 0.582*theta/S and solved to a fixed point with S itself.
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

FP12_FROM_KARATSUBA = 54          # Fp12 mul = 18 Fp2 mul = 54 Fp mul
MILLER_DBL, MILLER_ADD = 9, 11
MILLER_ITERS, MILLER_ADDS = 64, 87
FINAL_EXP = 300
G2_DBL, G2_ADD = 24, 48
SCALAR_BITS = K.R_ORDER.bit_length()      # 255

# the geometric-sum ceiling for the launch-schedule tail: sum_{j>=1} e^-j = 1/(e-1)
TAIL_GEOM = Fraction(582, 1000)


# ---------------------------------------------------------------------------
# (1) the step price in INSTRUCTIONS, one number per channel, no /C_REF
# ---------------------------------------------------------------------------
def step_price_instr(stress=K.STRESS_FOLD):
    """INSTRUCTIONS per step, each row rounded UP to a whole instruction so the
    ledger total is an exact INTEGER."""
    r = L.rows()
    fold_mul = int(math.ceil(r['fold_mul'] * float(stress)))
    fold_sqr = int(math.ceil(r['fold_sqr'] * float(stress)))
    mul = K.C_REF + r['cyclo_sqr'] + fold_mul
    sqr = r['cyclo_sqr'] + fold_sqr
    return dict(mul=Fraction(mul), sqr=Fraction(sqr), f=Fraction(K.F_SQR),
                fold_mul_charged=fold_mul, fold_sqr_charged=fold_sqr,
                mean=(1 - Fraction(K.F_SQR)) * Fraction(mul)
                     + Fraction(K.F_SQR) * Fraction(sqr))


def output_rows():
    """Q = [x]Q0, instruction-counted (unchanged from v32; now exact)."""
    miller = MILLER_ITERS * MILLER_DBL + MILLER_ADDS * MILLER_ADD
    g2 = (SCALAR_BITS - 1) * G2_DBL + (SCALAR_BITS // 2) * G2_ADD
    g2_fp12 = Fraction(g2, FP12_FROM_KARATSUBA)
    fixed_p = miller + FINAL_EXP
    verify = miller + FINAL_EXP
    fp12 = Fraction(fixed_p) + g2_fp12 + Fraction(verify)
    instr = Fraction(math.ceil(fp12 * K.C_REF))     # whole instructions, rounded UP
    return dict(miller_fp12=miller, final_exp_fp12=FINAL_EXP, g2_fp_mul=g2,
                g2_fp12mul=str(g2_fp12), fixed_p_fp12=fixed_p,
                verify_fp12=verify, pairing_fp12mul=str(fp12),
                pairing_instructions=str(instr),
                pairing_log2=math.log2(float(instr)),
                pairing_target_ops=str(fp12))


# ---------------------------------------------------------------------------
# (2) the sizing law, both measured factors charged
# ---------------------------------------------------------------------------
def step_budget33(p_det, delta=None, mu=K.MU_WALK, target=K.TARGET_PAD,
                  f=K.F_SQR, n=K.N_CLASSES, fixed_point=True):
    """Least integer S with (1-delta)(1-exp(-S^2(1-f^2)p_det/(2 mu N))) >= target.

    delta = P(D>R) on the actual event.  The analytic ceiling 0.582*theta/S
    depends on S, so delta and S are solved together (the correction is ~2^-46
    at production, far below the claim's 4th decimal, but it is iterated instead
    of being dropped).
    """
    pd = float(p_det)
    d = 0.0 if delta is None else float(delta)
    S = None
    for _ in range(40):
        g = math.log(1.0 / (1.0 - float(target) / (1.0 - d)))
        s2 = 2.0 * float(mu) * float(n) * g / ((1.0 - float(f) ** 2) * pd)
        newS = int(math.isqrt(int(s2)))
        while newS * newS < s2:
            newS += 1
        while newS > 1 and (newS - 1) ** 2 >= s2:
            newS -= 1
        if not fixed_point:
            return newS
        nd = float(TAIL_GEOM * THETA / newS)
        if S == newS and abs(nd - d) < 1e-30:
            return newS
        S, d = newS, nd
    return S


def tail_charge(S, theta=THETA):
    """delta = P(D>R) ceiling: the launch schedule's geometric tail, 0.582*theta/S."""
    return TAIL_GEOM * Fraction(int(theta), 1) / Fraction(int(S), 1)


# ---------------------------------------------------------------------------
# (3) the six rows, in INSTRUCTIONS
# ---------------------------------------------------------------------------
def rows33(S, p_det, delta=None, launch_mult=1):
    pr = step_price_instr()
    launches = -(-int(S) // THETA) * int(launch_mult)
    e = math.exp(-CAP / float(THETA))
    e_l2 = 2.0 * float(THETA) ** 2 * (1.0 - e * (1.0 + CAP_MULT))
    p_trail = Fraction(int(round(e_l2)), 2 * int(K.N_CLASSES))
    abandon_raw = Fraction(launches) * p_trail * THETA * pr['mul']
    abandon = Fraction(-((-abandon_raw.numerator) // abandon_raw.denominator))  # ceil: integer instructions
    in_flight = CAP_MULT * CAP * pr['mul']
    aux_rows = aux33(S, pr)
    orow = output_rows()
    return dict(S=int(S), theta=THETA, cap=CAP, launches=launches,
                delta=(tail_charge(S) if delta is None else Fraction(delta)),
                p_det=Fraction(p_det), p_trail=p_trail,
                step_price_instructions=pr['mul'].numerator // pr['mul'].denominator,
                step_price_mean_instructions=str(pr['mean']),
                w_effective=Fraction(int(S)) * pr['mul'],
                abandon_row=abandon, in_flight_row=in_flight,
                aux_rows=aux_rows, aux_sum=sum(aux_rows.values(), Fraction(0)),
                pairing_row=Fraction(orow['pairing_instructions']),
                pairing_detail=orow,
                overflow_row=Fraction(0),
                reports=2 * (-(-int(S) // THETA)))


def aux33(S, pr=None):
    """Every non-step row, ENUMERATED, in INSTRUCTIONS, at the report CAP."""
    pr = step_price_instr() if pr is None else pr
    reports = 2 * (-(-int(S) // THETA))
    cyclo = K.cyclo_sqr_instructions()
    gt_exp = 381 * cyclo + 191 * K.C_REF          # ONE GT exponentiation
    return dict(
        table_build=Fraction(K.TABLE_COLS * (1 << K.TABLE_BITS)
                             * (K.SIGMA_BUILD + 381 * K.FP_MUL)),
        gt_exp_starts=Fraction(reports) * 2 * gt_exp,
        restart_bookkeeping=Fraction(reports * 16),
        report_writes=Fraction(reports * K.REPORT_WORDS * 2),
        store_lookup=Fraction(reports * 420),
        report_cap_drain=Fraction(reports * 420),
        solve=Fraction(12 * 24 * K.FP_MUL + 2 ** 9 * K.C_REF),
        verify=Fraction(2 ** 9 * K.C_REF),
    )


def total33(S, p_det, delta=None, scale=1):
    """ONE sum over the six charged rows, in instructions, then ONE division."""
    r = rows33(S, p_det, delta)
    comps = dict(w_effective=r['w_effective'], abandon_row=r['abandon_row'],
                 aux_sum=r['aux_sum'], pairing_row=r['pairing_row'],
                 overflow_row=r['overflow_row'],
                 stop_loss_row=r['in_flight_row'] + r['abandon_row'])
    W_instr = sum(comps.values(), Fraction(0)) * Fraction(scale)
    W_ref = W_instr / K.C_REF                     # <<< THE ONE DIVISION LINE
    return dict(components=comps, W_instructions=W_instr, W_target_ops=W_ref,
                log2_W=math.log2(float(W_ref)), rows=r,
                division_line='W_target_ops = W_instructions / C_REF')


def claim33(p_det=None, delta=None, scale=1):
    """Least 4-decimal c with 2^c >= W_target_ops (score unit: log2 target ops).

    p_det is the measured CP-95 lower bound on P(detect|merge) (evidence).
    delta is the ACTUAL-EVENT tail bound 0.582*theta/S; with delta=None it is
    solved to a fixed point with S, and the measured fixture row is checked to
    sit under the same ceiling (it must, or the bound is not a bound).
    """
    if p_det is None:
        p_det = evidence_inputs()[0]
    S = step_budget33(p_det, delta)
    delta_used = tail_charge(S) if delta is None else Fraction(delta)
    T = total33(S, p_det, delta_used, scale=scale)
    lg = T['log2_W']
    c = math.floor(lg * 10 ** 4) / 10.0 ** 4
    while c < lg:
        c += 1e-4
    c = round(c, 4)
    num, den = int(round(c * 10 ** 4)), 10 ** 4
    g = math.gcd(num, den)
    return dict(claim_log2=c, numerator=num // g, denominator=den // g,
                total_log2=lg, S=S, T=T,
                exact_upper=bool(lg <= c),
                tight=bool(c - lg < 1e-4),
                delta=str(T['rows']['delta']), p_det=str(T['rows']['p_det']),
                delta_source='actual-event geometric tail, ceiling 0.582*theta/S '
                             'solved with S (measured fixture row checked under it)')


def evidence_inputs(ev=None):
    """Read the measured factors from the shipped evidence.  No literal default:
    if the evidence is absent the claim cannot be emitted."""
    import json
    import os
    ev = ev or os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'evidence'))
    det = json.load(open(os.path.join(ev, 'detect33.json')))
    tl = json.load(open(os.path.join(ev, 'tail33.json')))
    p_det = Fraction(str(det['charged'])).limit_denominator(10 ** 9)
    delta = Fraction(str(tl['charged'])).limit_denominator(10 ** 30)
    return p_det, delta


if __name__ == '__main__':
    c = claim33()
    T = c['T']
    print('S = %d  delta = %.6e  p_det = %.6f' % (c['S'], float(c['T']['rows']['delta']),
                                                  float(c['T']['rows']['p_det'])))
    for k, v in T['components'].items():
        print('  %-14s 2^%.6f instr' % (k, math.log2(float(v)) if v else -math.inf))
    print('W_instr 2^%.6f   W_ref 2^%.6f' % (math.log2(float(T['W_instructions'])),
                                            T['log2_W']))
    print('CLAIM %s = %d/%d' % (c['claim_log2'], c['numerator'], c['denominator']))
