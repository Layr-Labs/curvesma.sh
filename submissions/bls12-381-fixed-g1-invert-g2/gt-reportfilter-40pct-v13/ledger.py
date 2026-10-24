#!/usr/bin/env python3
"""Exact itemised ledger for the v13 EQUIVARIANT mu6-CLASS walk.

Every row is an instruction count at the pinned kernel price (constants.FP_MUL for
a base-field multiply, 383 measured instructions) or an organiser reference
operation (constants.C_REF = 20000, the dense Fp12 multiply of COST_MODEL.md).

What v13 charges that v12d (HOLD 126.7101) did, and what it no longer charges:

  * SELECTOR = 1613 on EVERY step.  Unchanged and unavoidable: the twelve class
    keys give BOTH the sigma-invariant routing/record key and the rotation index
    that selects the multiplier column.  It is v12b's counted row.
  * COLUMN_ADDR = 16 on the multiply channel: index = 12*t + kappa into the flat
    2^14 x 12 table.  Charged, not assumed free.
  * SIGMA_ON_STATE = 0.  v12d's 14 028 row is GONE because the step applies no
    sigma power to the state; the sigma powers are pre-applied to the twelve
    table columns and charged once in aux().  This is the answer to v12b's
    finding 2 that is a deletion, not a re-estimate.
  * RENORM = 0.  Nothing canonicalises, so no label is renormalised: the label
    update is a += u_t p^{-kappa}, one 255-bit modular add.
  * ESCAPE = 0.  No branch can adopt an inverted product, because no branch forms
    one.  f_sqr is the routing density 3/8 exactly, measured in T1.

The two ledgers below are

  W_raw    : each row at its derived bound,
  W_stress : every judgement-bounded row (the routing/selector fold) at 13/10;
             field-arithmetic rows (the cyclotomic squaring) and the counted
             selector stay at their counted values.

Both are exact Fractions; log2 is only for display.
"""
from decimal import Decimal
from fractions import Fraction
from math import log2
import constants as K


def rows():
    fold_mul = (K.KEY_HASH + K.DP_TEST + K.HIST_SLOT + K.LOOP + K.TABLE_FETCH
                + K.LABEL_ADD + K.SELECTOR + K.COLUMN_ADDR)
    fold_sqr = (K.KEY_HASH + K.DP_TEST + K.HIST_SLOT + K.LOOP + K.SELECTOR
                + 2*K.LABEL_ADD)                       # label doubling, no fetch
    cyc = 24*K.FP_MUL + 122*K.FP_ADD                   # Granger-Scott cyclotomic sqr
    return dict(
        dense_reference=K.C_REF, key_hash=K.KEY_HASH,
        dp_test=K.DP_TEST, hist_slot=K.HIST_SLOT, loop=K.LOOP,
        table_fetch=K.TABLE_FETCH, label_add=K.LABEL_ADD,
        selector=K.SELECTOR, column_addr=K.COLUMN_ADDR,
        fold_mul=fold_mul, fold_sqr=fold_sqr,
        cyclo_sqr=cyc, sigma_on_state=0, renorm=0, escape=0,
        mul_step=K.C_REF + fold_mul,
        sqr_step=cyc + fold_sqr,
    )


def ledgers(f=K.F_SQR, mu=K.MU_WALK, target=K.TARGET, stress=K.STRESS_FOLD,
            s=K.S_STEPS, n=K.N_CLASSES):
    r = rows()
    c_m = Fraction(r['mul_step'], K.C_REF)
    c_s = Fraction(r['sqr_step'], K.C_REF)
    c_raw = (1 - f)*c_m + f*c_s
    c_m_s = Fraction(K.C_REF + r['fold_mul']*stress, K.C_REF)
    c_s_s = Fraction(r['cyclo_sqr'] + r['fold_sqr']*stress, K.C_REF)
    c_stress = (1 - f)*c_m_s + f*c_s_s
    w_raw = s*c_raw
    w_str = s*c_stress
    return dict(S=s, c_mul=c_m, c_sqr=c_s, c_raw=c_raw,
                c_mul_stress=c_m_s, c_sqr_stress=c_s_s, c_stress=c_stress,
                w_raw=w_raw, w_stress=w_str, delta=w_str/w_raw - 1,
                model=model_success(s, n, f), f=f, mu=mu, target=target, rows=r)


def aux(s=K.S_STEPS, theta=1 << K.THETA_BITS, exp_ref=381, trails=4,
        record_bytes=96, solve_sigma=12):
    """Every non-step row of the attack, in reference-operation units, exactly."""
    entries = (1 << K.TABLE_BITS)*K.TABLE_COLS
    # each column entry = one sigma application to the base entry + its exponent
    table_build = Fraction(K.TABLE_COLS*(1 << K.TABLE_BITS)*(K.SIGMA_BUILD + exp_ref*K.FP_MUL),
                           K.C_REF)
    starts = Fraction(trails*2*exp_ref, 1)
    reports = s // theta
    rec = Fraction(int(reports)*(record_bytes//8)*12, K.C_REF)
    a = dict(
        table_build=table_build,
        trail_starts=starts,
        report_writes=rec,
        solve=Fraction(solve_sigma*24*K.FP_MUL + 2**9*K.C_REF, K.C_REF),
        verify=Fraction(2**9, 1),
    )
    a['total'] = sum(a.values())
    return a


def model_success(s, n, f):
    """1 - exp(-Lambda) with Lambda = S^2 (1-f^2) / (2 mu N), exact rational form."""
    lam = Fraction(s*s, 1) * (1 - f*f) / (2*K.MU_WALK*n)
    pw, term = Fraction(1, 1), Fraction(1, 1)
    for k in range(1, 60):
        term = term*lam/k
        pw += term
    return lam, pw


def log2_decimal(x, prec=60):
    """log2 of an exact positive Fraction, to `prec` significant Decimal digits."""
    from decimal import Decimal, getcontext
    getcontext().prec = prec
    if isinstance(x, Fraction):
        v = Decimal(x.numerator).ln() - Decimal(x.denominator).ln()
    else:
        v = Decimal(x).ln()
    return v / Decimal(2).ln()


def claim_log2(w_stress, pad=Decimal('2e-5'), dp=4, prec=60):
    """Least d-decimal claim c (a LOG2 value) with c >= log2(W_stress) + pad.

    The submitted score IS a log2 value, so the ceiling is taken in the log
    domain; the pad of 2e-5 bits is the family's rounding allowance.  Both
    directions are asserted, and the assertion is checked to sit away from the
    4-decimal boundary by more than 1e-30 bits so the Decimal log cannot
    misround.
    """
    from decimal import Decimal, getcontext
    getcontext().prec = prec
    logw = log2_decimal(w_stress, prec)
    target = logw + pad
    m = -((-target * 10**dp).to_integral_value(rounding='ROUND_FLOOR'))
    c = Decimal(m) / Decimal(10**dp)
    while c < target:
        m += 1
        c = Decimal(m) / Decimal(10**dp)
    while c - Decimal(1) / Decimal(10**dp) >= target:
        m -= 1
        c = Decimal(m) / Decimal(10**dp)
    assert c >= target and (c - Decimal(1) / Decimal(10**dp)) < target
    assert abs(c - target) > Decimal('1e-30')
    return c


def claim_floor_value(w_stress, pad=Fraction(1, 10**5), dp=4):
    """(reference only) least d-decimal VALUE >= w_stress + pad."""
    target = w_stress + pad
    num, den = target.numerator, target.denominator
    m = -((-num * 10**dp) // den)
    c = Fraction(m, 10**dp)
    while c < target:
        m += 1
        c = Fraction(m, 10**dp)
    return c


def l2(x):
    return log2(float(x))


if __name__ == '__main__':
    L = ledgers()
    r = L['rows']
    print("--- per-step rows (instructions)")
    for k in ('dense_reference', 'key_hash', 'dp_test', 'hist_slot', 'loop',
              'table_fetch', 'label_add', 'selector', 'column_addr',
              'fold_mul', 'fold_sqr', 'cyclo_sqr', 'sigma_on_state', 'renorm',
              'escape', 'mul_step', 'sqr_step'):
        print("  %-16s %7d" % (k, r[k]))
    print("fold_mul sum check:", r['fold_mul'] == sum(r[k] for k in
          ('key_hash', 'dp_test', 'hist_slot', 'loop', 'table_fetch',
           'label_add', 'selector', 'column_addr')))
    print("fold_sqr sum check:", r['fold_sqr'] == sum(r[k] for k in
          ('key_hash', 'dp_test', 'hist_slot', 'loop', 'selector')) + 2*K.LABEL_ADD)
    print("f_sqr = %s (declared %s, measured %s)  mu = %s  target = %s  S = %d (2^%.6f)"
          % (L['f'], K.F_SQR_DECLARED, K.F_SQR_MEASURED, L['mu'], L['target'],
             L['S'], log2(L['S'])))
    print("c_mul %.6f  c_sqr %.6f  c_raw %.6f" % (L['c_mul'], L['c_sqr'], L['c_raw']))
    print("W_raw    = 2^%.6f" % l2(L['w_raw']))
    print("W_stress = 2^%.6f   delta = %.6f" % (l2(L['w_stress']), L['delta']))
    C = claim_log2(L['w_stress'])
    lam, pw = L['model']
    print("claim    = %s (log2, pad 2e-5 bits)" % C)
    print("           above W_stress by %s bits, above W_raw by %s bits"
          % (C - log2_decimal(L['w_stress']), C - log2_decimal(L['w_raw'])))
    print("Lambda   = %.10f  model success = %.7f  (floor 0.39: %s)"
          % (float(lam), 1 - 1/float(pw), 1 - 1/float(pw) >= 0.39))
    A = aux()
    for k in ('table_build', 'trail_starts', 'report_writes', 'solve', 'verify', 'total'):
        print("  %-14s 2^%.3f  frac of W_stress 2^%.2f"
              % (k, l2(A[k]), l2(A[k]) - l2(L['w_stress'])))
    print("all aux below 2^-60 of W_stress:",
          all(l2(v) - l2(L['w_stress']) < -60 for k, v in A.items() if k != 'total'))
    print("legality: frontier 125.84, legality line 125.8345, claim %.4f" % float(C))
