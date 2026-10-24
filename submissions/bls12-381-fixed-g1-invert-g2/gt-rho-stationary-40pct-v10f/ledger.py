#!/usr/bin/env python3
"""Exact itemised ledger for the fixed-G1 -> G2 class walk, WINDOWLESS build.

Every row is an instruction count at the pinned kernel price (constants.FP_MUL for
a base-field multiply, 383 measured instructions) or an organiser reference
operation (constants.C_REF = 20000, the dense Fp12 multiply of COST_MODEL.md).

The only difference from the v10b ledger is that the 130-instruction look-back
window row is GONE: b68a086c voided the window (its A/B was confounded and a
per-trail history window breaks coalescence, since two trails at a common class
carry different windows and therefore diverge).  The fold is nine rows summing to
1813 = 403 + 766 + 180 + 264 + 72 + 8 + 32 + 32 + 56.  Nothing is padded: the
two ledgers below are

  W_raw    : each row at its derived bound,
  W_stress : every judgement-bounded row (the fold) carried at 13/10.

Both are exact Fractions; log2 is only for display.  The auxiliary rows (table
construction, trail starts, DP report writes, the solve and the verification) are
itemised separately and each is proved below 2^-60 of W_stress.
"""
from fractions import Fraction
from math import log2
import constants as K


def rows():
    """Per-step instruction rows.  Returns dict of named rows.

    A step is one dense Fp12 multiply (reference operation) or one cyclotomic
    squaring (non-reference: charged at its instruction count), plus the fold
    that turns the result into a class key.  The fold is paid on every step
    because the class key of the current state decides the next step type, so it
    cannot be deferred past a step.  The step is a stationary function of the
    class key: no history, no window, no escape, no abort.
    """
    # v10: the three-multiplication selector.  In this tower gamma = g0*(1-u)
    # exactly, so d = gamma0*a + gamma1*b = g0*(a-b) is ONE base-field multiply,
    # and the twelve keys are a, d, zeta*a, zeta*d and their zeta^2 forms
    # (zeta^2 x = zeta x - x, free).  Verified on the printed constants in
    # part_a.py: g1 = -g0, g0^6 = -1/8, zeta = 2 g0^2, and the three-multiply
    # form returns the general form's value on 200 random (a,b).
    selector_ab   = 1*K.FP_MUL + 1*K.FP_ADD   # a-b, then g0*(a-b)
    selector_zeta = 2*K.FP_MUL       # zeta*a, zeta*d  (zeta^2 forms free)
    fold_adds     = 9*K.FP_ADD       # a-b, two zeta^2 subtractions, six sign flips
    fold_cmp      = 11*K.CMP6        # six |.| tests + five-way maximum, with cmov
    hashing       = 3*K.HASH_W       # table index word + two distinguished words
    dp_test       = 8                # two-comparison distinguished test
    hist_slot     = 32               # 48-byte class-key traffic: two keys in/out per step
    loop          = 32               # table addressing, (twist) index add, counters
    witness_add   = 56               # a += u, b += v mod r: two 4-limb adds + conditional subtract
    fold = (selector_ab + selector_zeta + fold_adds + fold_cmp + hashing
            + dp_test + hist_slot + loop + witness_add)
    # Granger-Scott cyclotomic squaring, exact operation count of tower.cyclo_sqr_fast
    # (instrumented run, part_a prints the same number): two Fp6 squarings of six Fp2
    # squarings each (12 FpM + 37 FpA), the v-shift multiply, four Fp6 add/sub.
    cyc = 24*K.FP_MUL + 122*K.FP_ADD
    return dict(
        dense_reference=K.C_REF,
        selector_ab=selector_ab, selector_zeta=selector_zeta,
        fold_adds=fold_adds, fold_cmp=fold_cmp, hashing=hashing,
        dp_test=dp_test, hist_slot=hist_slot, loop=loop,
        witness_add=witness_add,
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


def aux(s=K.S_STEPS, theta=1 << K.THETA_BITS, table_entries=1 << 14, twists=12,
        record_bytes=96, trails=4, verify_muls=1 << 9):
    """Every non-step row of the attack, in reference-operation units, exactly.

    - table: 2^14 multipliers x 12 Frobenius forms, one dense Fp12 multiply each.
    - trail starts: 4 fresh starts (g^a z^b), one dense multiply each.
    - reports: S/theta expected reports; a 96-byte record write costs 12
      instructions per byte-group (12 stores of 8 bytes) / C.
    - solve: the exponent solve is 2x2 linear algebra in (a,b) mod r, charged at
      2^9 dense-multiply equivalents (it is a few dozen 4-limb ops in truth).
    - verify: one scalar multiplication and one comparison per recovery, 2^9.
    """
    reports = Fraction(s, theta)
    rec = Fraction(int(reports) * (record_bytes//8) * 12, K.C_REF)
    a = dict(
        table_build=Fraction(table_entries*twists, 1),
        trail_starts=Fraction(trails, 1),
        report_writes=rec,
        solve=Fraction(2**9, 1),
        verify=Fraction(verify_muls, 1),
    )
    a['total'] = sum(a.values())
    return a


def l2(x):
    return log2(float(x))


if __name__ == "__main__":
    L = ledgers()
    r = L['rows']
    print("--- per-step rows (instructions), windowless")
    for k in ('dense_reference','selector_ab','selector_zeta','fold_adds','fold_cmp',
              'hashing','dp_test','hist_slot','loop','witness_add','fold','cyclo_sqr',
              'mul_step','sqr_step'):
        print("  %-16s %7d" % (k, r[k]))
    print("fold sum check:", r['fold'] == sum(r[k] for k in
          ('selector_ab','selector_zeta','fold_adds','fold_cmp','hashing','dp_test',
           'hist_slot','loop','witness_add')))
    print("window/probe rows present:", 'detector' in r, "(must be False)")
    print("f = %s   mu = %s   target = %s" % (K.FRAC_SQR, K.MU_WALK, K.TARGET))
    print("S = %d  (log2 %.5f)" % (L['S'], log2(L['S'])))
    print("c_mul = %.6f  c_sqr = %.6f  c_raw = %.6f" %
          (L['c_mul'], L['c_sqr'], L['c_raw']))
    print("W_raw    = 2^%.5f" % l2(L['w_raw']))
    print("W_stress = 2^%.5f   delta = %.5f" % (l2(L['w_stress']), l2(L['w_stress'])-l2(L['w_raw'])))
    print("claim    = %s = 2^%.5f" % (K.CLAIM, l2(K.CLAIM)))
    print("raw pad   claim - W_raw    = %.5f  (rule >= 0.10)" % (float(K.CLAIM) - l2(L['w_raw'])))
    print("stress pad claim - W_stress = %.5f  (rule 0.005 .. 0.06)" % (float(K.CLAIM) - l2(L['w_stress'])))
    print("--- auxiliary rows (reference ops, exact)")
    A = aux()
    for k in ('table_build','trail_starts','report_writes','solve','verify','total'):
        print("  %-16s 2^%.3f   fraction of W_stress 2^%.2f"
              % (k, l2(A[k]), l2(A[k]) - l2(L['w_stress'])))
    print("all aux below 2^-60 of W_stress:",
          all(l2(v) - l2(L['w_stress']) < -60 for k, v in A.items() if k != 'total'))
    print("model success: 1-exp(-S^2/(2N)) = %.5f"
          % (1 - pow(2.718281828459045, -float(Fraction(L['S']*L['S'], 2*K.N_CLASSES)))))
