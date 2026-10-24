#!/usr/bin/env python3
"""Exact itemised ledger for the fixed-G1 -> G2 class walk, CANONICALISING build.

Every row is an instruction count at the pinned kernel price (constants.FP_MUL for
a base-field multiply, 383 measured instructions) or an organiser reference
operation (constants.C_REF = 20000, the dense Fp12 multiply of COST_MODEL.md).

Differences from the v10f ledger (27e3ab98):

  * the RENORM row is new and is paid on BOTH step types.  Canonicalising the
    state every step costs the label renormalisation (a,b) <- ((a+delta) p^k,
    b p^k) mod r: two 255-bit modular multiplications, counted at
    constants.MODMUL_INSTR = 143 each (evidence/modmul-count.json, clang -O2
    -arch arm64, 0/255 differential mismatches), so RENORM = 286.
  * the nine-row fold is unchanged at 1813 and stays windowless: no look-back,
    no checkpoint field, no escape.

The two ledgers below are

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
    that turns the result into a class key, plus the renormalisation that returns
    the state to the canonical representative.  The fold and the renormalisation
    are paid on every step because the next step's type is decided by the class
    key of the current state and the state must be canonical before the next
    table multiplication.
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
    loop          = 32               # table addressing, column index add, counters
    witness_add   = 56               # a += u, b += v mod r: two 4-limb adds + conditional subtract
    fold = (selector_ab + selector_zeta + fold_adds + fold_cmp + hashing
            + dp_test + hist_slot + loop + witness_add)
    # Granger-Scott cyclotomic squaring, exact operation count of tower.cyclo_sqr_fast
    # (instrumented run, part_a prints the same number): two Fp6 squarings of six Fp2
    # squarings each (12 FpM + 37 FpA), the v-shift multiply, four Fp6 add/sub.
    cyc = 24*K.FP_MUL + 122*K.FP_ADD
    renorm = K.RENORM
    return dict(
        dense_reference=K.C_REF,
        selector_ab=selector_ab, selector_zeta=selector_zeta,
        fold_adds=fold_adds, fold_cmp=fold_cmp, hashing=hashing,
        dp_test=dp_test, hist_slot=hist_slot, loop=loop,
        witness_add=witness_add,
        fold=fold, cyclo_sqr=cyc, renorm=renorm,
        mul_step=K.C_REF + fold + renorm,
        sqr_step=cyc + fold + renorm,
    )


def ledgers(f=K.FRAC_SQR, mu=K.MU_WALK, target=K.TARGET, stress=K.STRESS_FOLD,
            s=K.S_STEPS, n=K.N_CLASSES):
    r = rows()
    c_m = Fraction(r['mul_step'], K.C_REF)
    c_s = Fraction(r['sqr_step'], K.C_REF)
    c_raw = (1 - f)*c_m + f*c_s
    # stressed: judgement-bounded fold rows at 13/10, field-arithmetic and the
    # counted renormalisation row unchanged
    c_m_s = Fraction(K.C_REF + r['fold']*stress + r['renorm'], K.C_REF)
    c_s_s = Fraction(r['cyclo_sqr'] + r['fold']*stress + r['renorm'], K.C_REF)
    c_stress = (1 - f)*c_m_s + f*c_s_s
    w_raw = s*c_raw
    w_str = s*c_stress
    return dict(S=s, c_mul=c_m, c_sqr=c_s, c_raw=c_raw, w_raw=w_raw,
                c_stress=c_stress, w_stress=w_str,
                delta=w_str/w_raw - 1, rows=r)


def aux(s=K.S_STEPS, theta=1 << K.THETA_BITS, table_entries=1 << 14, twists=12,
        record_bytes=96, trails=4, verify_muls=1 << 9):
    """Every non-step row of the attack, in reference-operation units, exactly.

    - table: 2^14 rows x 12 independent columns, one dense Fp12 multiply each.
    - trail starts: 4 fresh starts (g^a z^b), one dense multiply each.
    - reports: S/theta expected reports; a 96-byte record write costs 12
      instructions per byte-group (12 stores of 8 bytes) / C.
    - solve: the exponent solve is now a single 2x2 cancellation in (a,b) mod r,
      x = (a1-a2)/(b2-b1), with NO Frobenius search: charged at 2^9 dense-multiply
      equivalents (it is a few dozen 4-limb ops plus one inversion).
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


if __name__ == '__main__':
    L = ledgers()
    r = L['rows']
    print("--- per-step rows (instructions), canonicalising")
    for k in ('dense_reference','selector_ab','selector_zeta','fold_adds','fold_cmp',
              'hashing','dp_test','hist_slot','loop','witness_add','fold','cyclo_sqr',
              'renorm','mul_step','sqr_step'):
        print("  %-16s %7d" % (k, r[k]))
    print("fold sum check:", r['fold'] == sum(r[k] for k in
          ('selector_ab','selector_zeta','fold_adds','fold_cmp','hashing','dp_test',
           'hist_slot','loop','witness_add')))
    print("window/probe rows present:", 'detector' in r, "(must be False)")
    print("f = %s   mu = %s   target = %s" % (K.FRAC_SQR, K.MU_WALK, K.TARGET))
    print("S = %d  (log2 %.6f)" % (L['S'], log2(L['S'])))
    print("c_mul = %.6f  c_sqr = %.6f  c_raw = %.6f" %
          (float(L['c_mul']), float(L['c_sqr']), float(L['c_raw'])))
    print("W_raw    = 2^%.6f" % l2(L['w_raw']))
    print("W_stress = 2^%.6f   delta = %.6f" %
          (l2(L['w_stress']), l2(L['w_stress'])-l2(L['w_raw'])))
    print("claim    = %s = 2^%.6f" % (K.CLAIM, l2(K.CLAIM)))
    print("stress pad claim - W_stress = %.6f  (rule >= 0.005)" %
          (float(K.CLAIM) - l2(L['w_stress'])))
    print("--- auxiliary rows (reference ops, exact)")
    A = aux()
    for k in ('table_build','trail_starts','report_writes','solve','verify','total'):
        print("  %-16s 2^%.3f   fraction of W_stress 2^%.2f"
              % (k, l2(A[k]), l2(A[k]) - l2(L['w_stress'])))
    print("all aux below 2^-60 of W_stress:",
          all(l2(v) - l2(L['w_stress']) < -60 for k, v in A.items() if k != 'total'))
    print("model success: 1-exp(-S^2/(2N)) = %.6f"
          % (1 - pow(2.718281828459045, -float(Fraction(L['S']*L['S'], 2*K.N_CLASSES)))))
