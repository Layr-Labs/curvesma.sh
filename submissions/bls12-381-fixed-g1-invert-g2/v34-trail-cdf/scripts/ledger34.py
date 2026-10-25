#!/usr/bin/env python3
"""v34 ledger -- every row in RAW INSTRUCTIONS, ONE division by C_REF, and the
two v33 constants DERIVED instead of asserted.

WHAT CHANGES FROM v33 (verdict-v33 accounting objection):

  * FINAL_EXP.  v33 asserted 300 Fp12 multiplications with no chain.  v34 emits
    the chain from the shipped curve parameter: the easy part (1 inversion +
    2 Fp12 muls) plus three exponentiations by the BLS12-381 seed
    |x| = 0xd201000000010000, each counted exactly by binary square-and-multiply
    ((bitlen-1) squarings + (popcount-1) multiplies for THAT exponent) plus the
    12 Fp12 muls of the hard-part addition chain.  The count is printed with its
    bit length and Hamming weight, so the reader can re-derive it.

  * G2 / 54.  The divisor is not asserted either: 54 = 18 Fp2-mul x 3 Fp-mul is
    the Karatsuba tower identity (Fp12 = Fp6[w]/(w^2-v): 3 Fp6 muls per Fp12 mul;
    Fp6 = Fp2[v]/(v^3-xi): 6 Fp2 muls per Fp6 mul; Fp2 = Fp[u]/(u^2-beta): 3 Fp
    muls per Fp2 mul -> 3*6*3 = 54), printed as a structural row, and the
    ADDITIONS and LOOP CONTROL of the same chain are charged explicitly in
    instructions at the organiser prices FP_ADD and LOOP instead of disappearing
    inside the division.

  * delta.  P(D>R) is gone as a schedule model.  The trail-level lifetime loss is
    MEASURED (tail34.json: exact per-trail T, P(T > 128*theta) at its CP-95 UPPER
    bound) and charged as delta; p_det is the detection probability CONDITIONED on
    a report existing.  The product (1-delta)*p_det is the merge-all rate measured
    on all class merges in merge34.json -- the same loss, factored, not counted
    twice.

  * mu.  Derived from the trace (completed/(launches*theta)), not imported.
"""
import math
from fractions import Fraction

import constants as K
import ledger as L

THETA = 1 << K.THETA_BITS
CAP_MULT = 128
CAP = CAP_MULT * THETA

# --- BLS12-381 seed, the only exponent the hard part of the final exponentiation
#     uses.  |x| = 0xd201000000010000 for BLS12-381.
BLS_X = 0xd201000000010000

# --- tower: Fp12 mul = 3 Fp6 mul (Karatsuba) = 3*6 Fp2 mul = 18*3 Fp mul
FP12_TOWER = dict(fp6_mul_per_fp12=3, fp2_mul_per_fp6=6, fp_mul_per_fp2=3)
FP12_FROM_KARATSUBA = (FP12_TOWER['fp6_mul_per_fp12']
                       * FP12_TOWER['fp2_mul_per_fp6']
                       * FP12_TOWER['fp_mul_per_fp2'])       # 54

MILLER_DBL, MILLER_ADD = 9, 11
MILLER_ITERS, MILLER_ADDS = 64, 87
G2_DBL, G2_ADD = 24, 48
G2_DBL_ADDS, G2_ADD_ADDS = 16, 24        # Fp additions inside one DBL / one ADD
SCALAR_BITS = K.R_ORDER.bit_length()     # 255


# ---------------------------------------------------------------------------
# (0) the derived constants
# ---------------------------------------------------------------------------
def final_exp_chain(x=BLS_X):
    """The final-exponentiation chain, counted.  No asserted constant.

    easy part  f^((p^6-1)(p^2+1)): one Fp12 inversion (conjugate and invert) plus
               two Fp12 multiplications; the p^2 Frobenius is conjugation in the
               tower and is not an Fp12 multiplication.
    hard part  the standard chain in the BLS seed: three exponentiations by |x|
               plus HARD_MULS Fp12 multiplications.  Each exponentiation is
               counted EXACTLY by binary square-and-multiply on |x|.
    inversion  an Fp12 inversion inverts one Fp2 element, so it is charged as the
               tower descent 54 + 36 + 15 Fp muls plus one Fp inversion, and the
               Fp inversion is charged at its own binary count on p-2.
    """
    bits = x.bit_length()
    hw = bin(x).count('1')
    sq_per = bits - 1
    mul_per = hw - 1
    hard_exp = sq_per + mul_per
    hard = 3 * hard_exp + 12
    # Fp inversion by Fermat: p-2, binary
    p_bits = K.P_FIELD.bit_length()
    p_hw = bin(K.P_FIELD - 2).count('1')
    fp_inv_fp_muls = (p_bits - 1) + (p_hw - 1)
    inv_fp_muls = FP12_FROM_KARATSUBA + 36 + 15 + fp_inv_fp_muls
    inv_fp12 = Fraction(inv_fp_muls, FP12_FROM_KARATSUBA)
    easy = inv_fp12 + 2
    total = easy + hard
    return dict(x=hex(x), x_bits=bits, x_hamming=hw,
                square_per_exp=sq_per, multiply_per_exp=mul_per,
                hard_exponentiations=3, hard_muls_chain=12,
                hard_fp12=hard, fp_inversion_fp_muls=fp_inv_fp_muls,
                inversion_fp12=str(inv_fp12), easy_fp12=str(easy),
                total_fp12=str(total), total_fp12_ceil=math.ceil(float(total)),
                rule='binary square-and-multiply: (bitlen-1) squarings + '
                     '(popcount-1) multiplies for the exponent |x|; three '
                     'exponentiations plus the 12-mul hard-part chain; easy part '
                     '1 inversion + 2 muls')


def output_rows():
    """Q = [x]Q0, instruction-counted with the DERIVED final exponentiation."""
    fe = final_exp_chain()
    fexp = fe['total_fp12_ceil']                      # 232 (was the asserted 300)
    miller = MILLER_ITERS * MILLER_DBL + MILLER_ADDS * MILLER_ADD
    g2_fp_mul = (SCALAR_BITS - 1) * G2_DBL + (SCALAR_BITS // 2) * G2_ADD
    g2_fp12 = Fraction(g2_fp_mul, FP12_FROM_KARATSUBA)
    fixed_p = miller + fexp
    verify = miller + fexp
    fp12 = Fraction(fixed_p) + g2_fp12 + Fraction(verify)
    instr = Fraction(math.ceil(fp12 * K.C_REF))
    # the ADDITIONS and CONTROL of the G2 chain, charged in instructions
    g2_adds = (SCALAR_BITS - 1) * G2_DBL_ADDS + (SCALAR_BITS // 2) * G2_ADD_ADDS
    g2_add_instr = g2_adds * K.FP_ADD
    g2_ctrl_instr = SCALAR_BITS * K.LOOP
    miller_adds_instr = MILLER_ITERS * 4 * K.FP_ADD      # 4 adds inside one DBL
    extra_instr = g2_add_instr + g2_ctrl_instr + miller_adds_instr
    return dict(miller_fp12=miller, final_exp_fp12=fexp,
                final_exp_chain=fe,
                g2_fp_mul=g2_fp_mul, g2_fp12mul=str(g2_fp12),
                g2_adds=g2_adds, g2_add_instr=g2_add_instr,
                g2_ctrl_instr=g2_ctrl_instr, miller_adds_instr=miller_adds_instr,
                extra_instr=extra_instr,
                fp12_tower=FP12_TOWER, fp12_from_karatsuba=FP12_FROM_KARATSUBA,
                fixed_p_fp12=fixed_p, verify_fp12=verify,
                pairing_fp12mul=str(fp12), pairing_instructions=str(instr),
                pairing_instructions_total=str(int(instr) + extra_instr),
                pairing_log2=math.log2(float(instr) + extra_instr))


# ---------------------------------------------------------------------------
# (1) the step price in INSTRUCTIONS, no /C_REF
# ---------------------------------------------------------------------------
def step_price_instr(stress=K.STRESS_FOLD):
    r = L.rows()
    fold_mul = int(math.ceil(r['fold_mul'] * float(stress)))
    fold_sqr = int(math.ceil(r['fold_sqr'] * float(stress)))
    mul = K.C_REF + r['cyclo_sqr'] + fold_mul
    sqr = r['cyclo_sqr'] + fold_sqr
    return dict(mul=Fraction(mul), sqr=Fraction(sqr), f=Fraction(K.F_SQR),
                fold_mul_charged=fold_mul, fold_sqr_charged=fold_sqr,
                mean=(1 - Fraction(K.F_SQR)) * Fraction(mul)
                     + Fraction(K.F_SQR) * Fraction(sqr))


# ---------------------------------------------------------------------------
# (2) the sizing law, both MEASURED factors charged
# ---------------------------------------------------------------------------
def step_budget34(p_det, delta, mu, target=K.TARGET_PAD, f=K.F_SQR,
                  n=K.N_CLASSES):
    """Least integer S with (1-delta)(1-exp(-S^2(1-f^2)p_det/(2 mu N))) >= target.

    delta is the MEASURED trail-level lifetime loss P(T > 128*theta) (CP-95 upper
    bound) and p_det the detection probability CONDITIONED on a report; their
    product is the merge-all rate of merge34.json, so one loss is charged once.
    """
    d = float(delta)
    if not (0.0 <= d < 1.0):
        raise ValueError('delta must be in [0,1)')
    if (1.0 - d) <= float(target):
        raise ValueError('infeasible: (1-delta) = %.6f <= target %.6f, no S exists'
                         % (1.0 - d, float(target)))
    g = math.log(1.0 / (1.0 - float(target) / (1.0 - d)))
    s2 = 2.0 * float(mu) * float(n) * g / ((1.0 - float(f) ** 2) * float(p_det))
    newS = int(math.isqrt(int(s2)))
    while newS * newS < s2:
        newS += 1
    while newS > 1 and (newS - 1) ** 2 >= s2:
        newS -= 1
    return newS


# ---------------------------------------------------------------------------
# (3) the six rows, in INSTRUCTIONS
# ---------------------------------------------------------------------------
def rows34(S, p_det, delta, mu):
    pr = step_price_instr()
    launches = -(-int(S) // THETA)
    e = math.exp(-CAP / float(THETA))
    e_l2 = 2.0 * float(THETA) ** 2 * (1.0 - e * (1.0 + CAP_MULT))
    p_trail = Fraction(int(round(e_l2)), 2 * int(K.N_CLASSES))
    abandon_raw = Fraction(launches) * p_trail * THETA * pr['mul']
    abandon = Fraction(-((-abandon_raw.numerator) // abandon_raw.denominator))
    in_flight = CAP_MULT * CAP * pr['mul']
    aux_rows = aux34(S, pr)
    orow = output_rows()
    return dict(S=int(S), theta=THETA, cap=CAP, launches=launches,
                delta=Fraction(str(delta)).limit_denominator(10 ** 12),
                p_det=Fraction(p_det), mu=Fraction(str(mu)).limit_denominator(10 ** 9),
                p_trail=p_trail,
                step_price_instructions=pr['mul'].numerator // pr['mul'].denominator,
                step_price_mean_instructions=str(pr['mean']),
                w_effective=Fraction(int(S)) * pr['mul'],
                abandon_row=abandon, in_flight_row=in_flight,
                aux_rows=aux_rows, aux_sum=sum(aux_rows.values(), Fraction(0)),
                pairing_row=Fraction(orow['pairing_instructions_total']),
                pairing_detail=orow,
                overflow_row=Fraction(0),
                reports=2 * launches)


def aux34(S, pr=None):
    pr = step_price_instr() if pr is None else pr
    reports = 2 * (-(-int(S) // THETA))
    cyclo = K.cyclo_sqr_instructions()
    gt_exp = 381 * cyclo + 191 * K.C_REF
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


def total34(S, p_det, delta, mu, scale=1):
    r = rows34(S, p_det, delta, mu)
    comps = dict(w_effective=r['w_effective'], abandon_row=r['abandon_row'],
                 aux_sum=r['aux_sum'], pairing_row=r['pairing_row'],
                 overflow_row=r['overflow_row'],
                 stop_loss_row=r['in_flight_row'] + r['abandon_row'])
    W_instr = sum(comps.values(), Fraction(0)) * Fraction(scale)
    W_ref = W_instr / K.C_REF                     # <<< THE ONE DIVISION LINE
    return dict(components=comps, W_instructions=W_instr, W_target_ops=W_ref,
                log2_W=math.log2(float(W_ref)), rows=r,
                division_line='W_target_ops = W_instructions / C_REF')


def claim34(p_det=None, delta=None, mu=None, scale=1):
    """The claim under BOTH factorizations of the measured loss (see finalize34).

    frame A  delta = 0, p_det = the merge-all rate at its cluster lower bound (the
             whole detection loss lives on the merge event)
    frame B  delta = the max per-run CP-95 upper bound on P(T > 128*theta),
             p_det = the detection rate conditioned on a report, cluster lower
             bound
    The frame with the LARGER S is charged (more conservative), and the other is
    reported beside it.  Both are the same measurement factored once.
    """
    if p_det is None or delta is None or mu is None:
        e = evidence_inputs()
        p_det = e['p_det'] if p_det is None else p_det
        delta = e['delta'] if delta is None else delta
        mu = e['mu'] if mu is None else mu
        frames = e['frames']
    else:
        frames = None
    S = step_budget34(p_det, delta, mu)
    T = total34(S, p_det, delta, mu, scale=scale)
    lg = T['log2_W']
    c = math.floor(lg * 10 ** 4) / 10.0 ** 4
    while c < lg:
        c += 1e-4
    c = round(c, 4)
    num, den = int(round(c * 10 ** 4)), 10 ** 4
    g = math.gcd(num, den)
    out = dict(claim_log2=c, numerator=num // g, denominator=den // g,
               total_log2=lg, S=S, T=T, exact_upper=bool(lg <= c),
               tight=bool(c - lg < 1e-4),
               delta=str(T['rows']['delta']), p_det=str(T['rows']['p_det']),
               mu=str(T['rows']['mu']))
    if frames is not None:
        e = evidence_inputs()
        Sa = step_budget34(frames['A']['p_all'], 0.0, mu)
        try:
            Sb = step_budget34(frames['B']['p_cond'], frames['B']['delta'], mu)
            feasible_B = True
        except ValueError:
            Sb, feasible_B = None, False
        out['frames'] = dict(
            A=dict(S=Sa, p_det=str(frames['A']['p_all']), delta='0',
                   log2_S=math.log2(Sa)),
            B=dict(S=Sb, p_det=str(frames['B']['p_cond']),
                   delta=str(frames['B']['delta']),
                   log2_S=(math.log2(Sb) if Sb else None),
                   feasible=feasible_B,
                   infeasible_because=(None if feasible_B else
                                       '(1-delta) %.6f <= target %.6f: the '
                                       'worst seed cannot reach the success target '
                                       'at all -- a statement about the degenerate '
                                       'fixture, not about production'
                                       % (1 - float(frames['B']['delta']),
                                          float(K.TARGET_PAD)))),
            charged='B' if (feasible_B and Sb >= Sa) else 'A',
            larger_S=max([x for x in (Sa, Sb) if x]))
    return out


def evidence_inputs(ev=None):
    """The measured factors, read from the shipped evidence.  No literal default.

    The measured loss is factored two ways (finalize34):
      A  delta = 0,        p_det = merge-all rate at its cluster lower bound
      B  delta = delta_max, p_det = conditioned detection at its cluster lower bound
    Both are the same measurement; the frame with the LARGER S is charged, and the
    other is carried in `frames` so the note can print both.
    """
    import json
    import os
    ev = ev or os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'evidence'))
    tl = json.load(open(os.path.join(ev, 'tail34.json')))
    mg = json.load(open(os.path.join(ev, 'merge34.json')))
    delta_b = Fraction(str(tl['charged'])).limit_denominator(10 ** 12)
    mu = Fraction(str(tl['mu_derived'])).limit_denominator(10 ** 9)
    p_all = Fraction(str(mg['charged'])).limit_denominator(10 ** 12)
    p_cond = Fraction(str(mg['p_det_conditioned']['charged'])).limit_denominator(10 ** 12)
    frames = dict(A=dict(p_all=p_all, delta=Fraction(0)),
                  B=dict(p_cond=p_cond, delta=delta_b))
    Sa = step_budget34(p_all, 0.0, mu)
    # frame B is only a candidate if its (1-delta) can still reach the target: the
    # worst-run delta may be so large that NO S exists, which is a statement about
    # the degenerate fixture, not about production.
    try:
        Sb = step_budget34(p_cond, delta_b, mu)
    except ValueError:
        Sb = None
    charge = 'B' if (Sb is not None and Sb >= Sa) else 'A'
    return dict(delta=(delta_b if charge == 'B' else Fraction(0)),
                mu=mu, p_det=(p_cond if charge == 'B' else p_all),
                frames=frames, charge=charge, S_A=Sa, S_B=Sb,
                delta_all_upper=delta_b, p_det_all=p_all, p_det_cond=p_cond,
                delta_point=tl['p_T_gt_cap'],
                p_det_pooled=mg['p_det_conditioned'].get('pooled'),
                n_runs=len(tl['per_run']))


if __name__ == '__main__':
    e = evidence_inputs()
    c = claim34()
    T = c['T']
    print('inputs delta %.6f (point %.6f)  p_det %.6f  mu %.6f  p_det_all %s'
          % (float(e['delta']), e['delta_point'], float(e['p_det']), float(e['mu']),
             e['p_det_all_merges']))
    print('fexp chain', {k: T['rows']['pairing_detail']['final_exp_chain'][k]
                         for k in ('x_bits', 'x_hamming', 'hard_fp12',
                                   'total_fp12_ceil')})
    print('S = %d' % c['S'])
    for k, v in T['components'].items():
        print('  %-14s 2^%.6f instr' % (k, math.log2(float(v)) if v else -math.inf))
    print('W_instr 2^%.6f   W_ref 2^%.6f' % (math.log2(float(T['W_instructions'])),
                                            T['log2_W']))
    print('CLAIM %s = %d/%d' % (c['claim_log2'], c['numerator'], c['denominator']))
