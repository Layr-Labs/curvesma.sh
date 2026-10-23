#!/usr/bin/env python3
"""Exact itemised ledger for the fixed-G1 -> G2 class walk.

Every row is an instruction count at the pinned kernel price (constants.FP_MUL for
a base-field multiply, 383 measured instructions) or an organiser reference
operation (constants.C_REF).  Nothing is padded: the two ledgers below are
  W_raw    : each row at its derived bound,
  W_stress : every judgement-bounded row (the fold) carried at 13/10.
Both are exact Fractions; log2 is only for display.
"""
from fractions import Fraction
from math import log2
import constants as K


def rows():
    """Per-step instruction rows.  Returns dict of named rows.

    A step is one dense Fp12 multiply (reference operation) or one cyclotomic
    squaring (non-reference: charged at its instruction count), plus the fold
    that turns the result into a class key.  The fold is paid on every step
    because the class key of the current state decides the next step type, so
    it cannot be deferred past a step.
    """
    selector_ab   = 2*K.FP_MUL       # A = Re(d0), B = Re(gamma d0)
    selector_zeta = 4*K.FP_MUL       # A*zeta, A*zeta^2, B*zeta, B*zeta^2
    fold_adds     = 4*K.FP_ADD       # B accumulation, absolute-value negations
    fold_cmp      = 11*K.CMP6        # six |.| tests + five-way maximum, with cmov
    hashing       = 3*K.HASH_W       # table index word + two distinguished words
    dp_test       = 8                # two-comparison distinguished test
    detector      = 2*K.CMP6 + 16    # two-cycle rule: key vs previous, min slot
    hist_slot     = 32               # 48-byte previous-key copy in/out
    loop          = 32               # table addressing, counters, control
    fold = (selector_ab + selector_zeta + fold_adds + fold_cmp + hashing
            + dp_test + detector + hist_slot + loop)
    # Granger-Scott cyclotomic squaring, exact operation count of tower.cyclo_sqr_fast
    # (instrumented run, part_a prints the same number): two Fp6 squarings of six Fp2
    # squarings each (12 FpM + 37 FpA), the v-shift multiply, four Fp6 add/sub.
    cyc = 24*K.FP_MUL + 122*K.FP_ADD
    return dict(
        dense_reference=K.C_REF,
        selector_ab=selector_ab, selector_zeta=selector_zeta,
        fold_adds=fold_adds, fold_cmp=fold_cmp, hashing=hashing,
        dp_test=dp_test, detector=detector, hist_slot=hist_slot, loop=loop,
        fold=fold, cyclo_sqr=cyc,
        mul_step=K.C_REF + fold, sqr_step=cyc + fold,
    )


def ledgers(f=K.FRAC_SQR, mu=K.MU_WALK, target=K.TARGET, stress=K.STRESS_FOLD,
            s=K.S_STEPS, n=K.N_CLASSES):
    r = rows()
    c_m = Fraction(r['mul_step'], K.C_REF)
    c_s = Fraction(r['sqr_step'], K.C_REF)
    c_raw = (1 - f)*c_m + f*c_s
    # stressed: judgement-bounded fold rows at 13/10, field-arithmetic rows unchanged
    c_m_s = Fraction(K.C_REF + r['fold']*stress, K.C_REF)
    c_s_s = Fraction(r['cyclo_sqr'] + r['fold']*stress, K.C_REF)
    c_stress = (1 - f)*c_m_s + f*c_s_s
    w_raw = s*c_raw
    w_str = s*c_stress
    return dict(S=s, c_mul=c_m, c_sqr=c_s, c_raw=c_raw, w_raw=w_raw,
                c_stress=c_stress, w_stress=w_str,
                delta=w_str/w_raw - 1, rows=r)


def l2(x):
    return log2(float(x))


if __name__ == "__main__":
    L = ledgers()
    r = L['rows']
    print("--- per-step rows (instructions)")
    for k in ('dense_reference','selector_ab','selector_zeta','fold_adds','fold_cmp',
              'hashing','dp_test','detector','hist_slot','loop','fold','cyclo_sqr',
              'mul_step','sqr_step'):
        print("  %-16s %7d" % (k, r[k]))
    print("fold sum check:", r['fold'] == sum(r[k] for k in
          ('selector_ab','selector_zeta','fold_adds','fold_cmp','hashing','dp_test',
           'detector','hist_slot','loop')))
    print("f = %s   mu = %s   target = %s" % (K.FRAC_SQR, K.MU_WALK, K.TARGET))
    print("S = %d  (log2 %.5f)" % (L['S'], log2(L['S'])))
    print("c_mul = %.6f  c_sqr = %.6f  c_raw = %.6f" %
          (L['c_mul'], L['c_sqr'], L['c_raw']))
    print("W_raw    = 2^%.5f" % l2(L['w_raw']))
    print("W_stress = 2^%.5f   delta = %.5f" % (l2(L['w_stress']), l2(L['w_stress'])-l2(L['w_raw'])))
    print("claim    = %s = 2^%.5f" % (K.CLAIM, l2(K.CLAIM)))
    print("raw pad   claim - W_raw    = %.5f  (rule >= 0.10)" % (float(K.CLAIM) - l2(L['w_raw'])))
    print("stress pad claim - W_stress = %.5f  (rule 0.005 .. 0.06)" % (float(K.CLAIM) - l2(L['w_stress'])))
    print("gate: claim < 125.90 :", float(K.CLAIM) < 125.90)
    print("success budget: 1-exp(-S^2(1-f^2)/(2 mu N)) = %.4f" %
          (1 - pow(2.718281828459045,
                   -L['S']*L['S']*(1-float(K.FRAC_SQR)**2)/(2*float(K.MU_WALK)*K.N_CLASSES))))
