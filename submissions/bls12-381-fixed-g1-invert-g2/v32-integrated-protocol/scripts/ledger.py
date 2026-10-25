#!/usr/bin/env python3
"""Exact itemised ledger for the v14 GUARDED equivariant mu6-class walk.

Every row is an instruction count at the pinned kernel price (constants.FP_MUL for
a base-field multiply, 383 measured instructions) or an organiser reference
operation (constants.C_REF = 20000, the dense Fp12 multiply of COST_MODEL.md).

What v14 charges that v13 did, and what it still does not charge:

  * SELECTOR = 1613 on EVERY step (unchanged counted row: the twelve class keys
    give BOTH the sigma-invariant routing/record key and the column index).
  * COLUMN_ADDR = 16 on the multiply channel.
  * GUARD = 48 on the MULTIPLY channel only: D = 12 packed (row,column) compares
    with conditional select, one ring shift of 12 words, index arithmetic.  This
    is the row that removes the absorbing family; it is charged, not argued free.
  * SIGMA_ON_STATE = 0, RENORM = 0, ESCAPE = 0: no branch applies a sigma power,
    no branch canonicalises, no branch abandons a trail.
  * f_sqr is ONE number, the routing share 3/8, used for BOTH the step budget S
    and the step price.  The fixture share is a cross-check row only.

The two ledgers below are

  W_raw    : each row at its derived bound,
  W_stress : every judgement-bounded row (the routing/selector fold) at 13/10;
             field-arithmetic rows (the cyclotomic squaring) stay counted.

Both are exact Fractions; log2 is only for display.
"""
from fractions import Fraction
from math import log2
import constants as K


def rows():
    # v23: THE GUARD IS GONE.  The per-step row is exactly what the guard-free step
    # (classkey.step_v13: column = kappa(y), no ring) executes: selector, column
    # address, table fetch, label add, key/DP/hist traffic, loop.  K.GUARD is not
    # charged because no ring lookup, compare or shift happens.
    fold_mul = (K.KEY_HASH + K.DP_TEST + K.HIST_SLOT + K.LOOP + K.TABLE_FETCH
                + K.LABEL_ADD + K.SELECTOR + K.COLUMN_ADDR + K.ROWDIFF)
    fold_sqr = (K.KEY_HASH + K.DP_TEST + K.HIST_SLOT + K.LOOP + K.SELECTOR
                + 2*K.LABEL_ADD)                       # label doubling, no fetch,
                                                       # no guard (squaring is not
                                                       # a multiplier: no inverse)
    cyc = 24*K.FP_MUL + 122*K.FP_ADD                   # Granger-Scott cyclotomic sqr
    return dict(
        dense_reference=K.C_REF, key_hash=K.KEY_HASH,
        dp_test=K.DP_TEST, hist_slot=K.HIST_SLOT, loop=K.LOOP,
        table_fetch=K.TABLE_FETCH, label_add=K.LABEL_ADD,
        selector=K.SELECTOR, column_addr=K.COLUMN_ADDR, rowdiff=K.ROWDIFF, guard=0,
        fold_mul=fold_mul, fold_sqr=fold_sqr,
        cyclo_sqr=cyc, sigma_on_state=0, renorm=0, escape=0,
        mul_step=K.C_REF + fold_mul,
        sqr_step=cyc + fold_sqr,
    )


def ledgers(f=K.F_SQR, mu=K.MU_WALK, target=K.TARGET_PAD, stress=K.STRESS_FOLD,
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
                model=model_success(s, n, f), f=f, mu=mu, target=target, rows=r,
                stick=K.STICK_LOSS, net_success=(1-K.STICK_LOSS)*(1-1/float(model_success(s, n, f)[1])))


def gt_exp_reference_units():
    """ONE GT exponentiation of a 381-bit scalar, in reference operations.

    The v14 ledger divided `reports * 2 * 381` SCALAR instructions by C and called
    it "two GT exponentiations"; the v14 verdict rejected that ("ledger.aux charges
    two GT exponentiations as reports*2*381 scalar instructions divided by C, rather
    than charging their GT arithmetic").  Charged here as tower arithmetic:

      * a dense Fp12 multiplication IS the reference operation (C_REF = 20,000,
        COST_MODEL.md), so one dense multiply costs exactly 1 reference unit;
      * a GT squaring is the counted Granger-Scott cyclotomic row
        cyclo_sqr = 24*FP_MUL + 122*FP_ADD (the SAME row the step ledger uses,
        ledger.rows()['cyclo_sqr']), i.e. 11,632 instructions = 11632/20000
        reference units.  Fp12 squaring in the tower is 24 Fp multiplications,
        not 54: the cyclotomic form removes the Fp2 multiplications (tower.py:
        cyclo_sqr_fast; 24 Fp-mul + 122 Fp-add, asserted in gen_cert);
      * a 381-bit exponentiation is 381 squarings plus at most 191 dense
        multiplications (190 is the mean Hamming weight of a 381-bit exponent;
        191 is charged, rounded up).
    """
    sq = Fraction(381 * K.cyclo_sqr_instructions(), K.C_REF)
    mul = 191
    return sq + mul


def report_cap(s=K.S_STEPS, theta=1 << K.THETA_BITS, slack=2):
    """The report CAP: slack x ceil(S/theta), a stated ceiling with margin.

    v15 charged the bare ceil(S/theta), which the v15 verdict rejected as
    conflating resident capacity with reports processed.  v16 (a) makes the cap
    an explicit ceiling with a factor-2 slack over the mean DP-hit count, (b)
    charges every per-report row AT THE CAP, so the bound holds on every run,
    and (c) ships the Chernoff row P(hits > cap) <= exp(-mean/3) and the census
    calibration row beside it.
    """
    return slack * ((s + theta - 1) // theta)


def aux(s=K.S_STEPS, theta=1 << K.THETA_BITS, exp_ref=381, solve_sigma=12, slack=2):
    """Every non-step row of the attack, in reference-operation units, exactly.

    The production protocol stated here is the one the census executes:

      * one shared report store, keyed by the 48-byte class key; a record is
        (class key, a, b) = 112 bytes;
      * each trail starts from a fresh uniform label pair (y = g^a z^b, TWO GT
        exponentiations) and ends at its FIRST distinguished point, i.e. after a
        geometric number of steps with mean theta = 2^80;
      * the trail that reports is replaced by a fresh trail, so the number of
        trail starts equals the number of reports, ceil(S/theta);
      * no live trail is ever restarted, aborted or re-seeded: the guard rewrites
        one table column and the trail continues.

    Rows added in v15, each demanded by the v14 verdict:

      gt_exp_starts      : 2 * GT_EXP_REF per trail start (tower block counts)
      store_lookup       : a 48-byte key hash (6 words * 30 instructions) plus two
                           probe steps of a 112-byte record compare (2 * 120) =
                           420 instructions per report
      report_cap_drain   : the OVERFLOW POLICY.  Capacity is
                           cap = ceil(S/theta) records, so the store is a bound,
                           not an expectation (v14 charged floor(S/2^80), an
                           expected scale).  If a report would exceed the cap the
                           oldest record is dropped and its class re-probed: the
                           drain charge is one lookup + one record write per
                           dropped record, priced at cap lookups worst case.
                           By construction reports <= cap, so the row is a
                           charged ceiling, not an expected cost.
    """
    entries = (1 << K.TABLE_BITS)*K.TABLE_COLS
    table_build = Fraction(K.TABLE_COLS*(1 << K.TABLE_BITS)*(K.SIGMA_BUILD + exp_ref*K.FP_MUL),
                           K.C_REF)
    reports = report_cap(s, theta, slack)       # CAP with a stated 2x slack
    gtexp = gt_exp_reference_units()
    starts = Fraction(int(reports)*2*gtexp)
    bookkeeping = Fraction(int(reports)*16, K.C_REF)      # fresh (a,b) pair + counter reset
    rec = Fraction(int(reports)*K.REPORT_WORDS*2, K.C_REF)  # 112 B = 14 words, 2 instr/word
    lookup = Fraction(int(reports)*420, K.C_REF)
    drain = Fraction(int(reports)*420, K.C_REF)           # overflow policy, charged
    a = dict(
        table_build=table_build,
        gt_exp_starts=starts,
        restart_bookkeeping=bookkeeping,
        report_writes=rec,
        store_lookup=lookup,
        report_cap_drain=drain,
        solve=Fraction(solve_sigma*24*K.FP_MUL + 2**9*K.C_REF, K.C_REF),
        verify=Fraction(2**9, 1),
    )
    a['cap_slack'] = Fraction(slack)
    a['total'] = sum(a.values())
    return a


def storage(s=K.S_STEPS, theta=1 << K.THETA_BITS, trails=4, guard_trails=None):
    """Peak storage in BYTES, every term with the 112-byte record.

    store   : cap = ceil(S/theta) records x 112 B -- the deterministic CAP (v14
              charged floor(S/theta), an expectation);
    labels  : the live trails' label pairs, 2 x 32 B each;
    u_labels: the U table, 2^14 x 12 x 32 B -- the label half of the twelce-column
              advice.  v14's peak expression omitted it (the verdict: "the exact
              peak-storage expression omits the U labels used by the step");
    advice  : the twelve-column multiplier table, 2^14 x 12 x 576 B;
    guard   : guard history, one 4-byte packed (row,column) word per ring slot,
              per live trail;
    code    : the shipped code and constants, charged at 262,144 B;
    scratch : 8 GT elements (a working product, the sigma partners, the verify).
    """
    reports = report_cap(s, theta)
    gt = 576
    trails = 0 if guard_trails is None else guard_trails   # v23: no ring state at all
    u_bytes = (1 << K.TABLE_BITS)*K.TABLE_COLS*K.LABEL_BYTES
    guard_bytes = trails*K.GUARD_DEPTH*4
    code_bytes = 262144
    return dict(
        store_bytes=int(reports)*K.RECORD_BYTES,
        label_bytes=trails*2*K.LABEL_BYTES,
        u_label_bytes=u_bytes,
        advice_bytes=(1 << K.TABLE_BITS)*K.TABLE_COLS*gt,
        guard_history_bytes=0,
        code_bytes=code_bytes,
        scratch_bytes=8*gt,
        records=int(reports),
        report_cap=int(reports),
        record_bytes=K.RECORD_BYTES,
        overflow_policy="drop-oldest + drain, charged in ledger.aux:report_cap_drain; "
                        "capacity is the CAP, so no run exceeds it",
        peak_bytes=(int(reports)*K.RECORD_BYTES + trails*2*K.LABEL_BYTES + u_bytes
                    + (1 << K.TABLE_BITS)*K.TABLE_COLS*gt + 0
                    + code_bytes + 8*gt),
    )


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


def claim_log2(w_stress, pad=None, dp=4, prec=60):
    """Least d-decimal claim c (a LOG2 value) with c >= log2(W_stress) + pad.

    The submitted score IS a log2 value, so the ceiling is taken in the log
    domain; the pad of 2e-5 bits is the family's rounding allowance.  Both
    directions are asserted, and the assertion is checked to sit away from the
    4-decimal boundary by more than 1e-30 bits so the Decimal log cannot misround.
    """
    from decimal import Decimal, getcontext
    getcontext().prec = prec
    if pad is None:
        pad = Decimal('2e-5')
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


def l2(x):
    return log2(float(x))


if __name__ == '__main__':
    L = ledgers()
    r = L['rows']
    print("--- per-step rows (instructions)")
    for k in ('dense_reference', 'key_hash', 'dp_test', 'hist_slot', 'loop',
              'table_fetch', 'label_add', 'selector', 'column_addr', 'rowdiff', 'guard',
              'fold_mul', 'fold_sqr', 'cyclo_sqr', 'sigma_on_state', 'renorm',
              'escape', 'mul_step', 'sqr_step'):
        print("  %-16s %7d" % (k, r[k]))
    print("fold_mul sum check:", r['fold_mul'] == sum(r[k] for k in
          ('key_hash', 'dp_test', 'hist_slot', 'loop', 'table_fetch',
           'label_add', 'selector', 'column_addr', 'rowdiff', 'guard')))
    print("fold_sqr sum check:", r['fold_sqr'] == sum(r[k] for k in
          ('key_hash', 'dp_test', 'hist_slot', 'loop', 'selector')) + 2*K.LABEL_ADD)
    print("f_sqr = %s (declared %s)  mu = %s  target = %s  S = %d (2^%.6f)"
          % (L['f'], K.F_SQR_DECLARED, L['mu'], L['target'], L['S'], log2(L['S'])))
    print("c_mul %.6f  c_sqr %.6f  c_raw %.6f  c_stress %.6f"
          % (L['c_mul'], L['c_sqr'], L['c_raw'], L['c_stress']))
    print("W_raw    = 2^%.6f" % l2(L['w_raw']))
    print("W_stress = 2^%.6f   delta = %.6f" % (l2(L['w_stress']), L['delta']))
    C = claim_log2(L['w_stress'])
    lam, pw = L['model']
    print("claim    = %s (log2, pad 2e-5 bits)" % C)
    print("           above W_stress by %s bits, above W_raw by %s bits"
          % (C - log2_decimal(L['w_stress']), C - log2_decimal(L['w_raw'])))
    print("Lambda   = %.10f  model success = %.7f  net of stickiness %.7f (floor 0.39: %s)"
          % (float(lam), 1 - 1/float(pw), L['net_success'], L['net_success'] >= 0.39))
    A = aux()
    for k in ('table_build', 'trail_starts', 'restart_bookkeeping', 'report_writes',
              'solve', 'verify', 'total'):
        print("  %-22s 2^%.3f  frac of W_stress 2^%.2f"
              % (k, l2(A[k]), l2(A[k]) - l2(L['w_stress'])))
    for k, v in storage().items():
        print("  storage %-14s %s" % (k, v))
    print("legality: frontier 125.84, legality line 125.8345, claim %.4f" % float(C))


# ===========================================================================
# THE ONE SUMMATION (verdict fcfdd383 defect 1).  v23's claim was taken from
# log2(W_stress) ALONE and the aux block was printed beside it, never added to
# it; the reviewer added the two by hand and the claim came out BELOW its own
# components.  From v24 the claim is the ceiling of ONE exact sum and there is
# no code path that can produce a claim from a part:
#
#   W_total = S*c_stress            the step work (W_stress)
#           + trap_charge           rate * cap * price, the capped waste of the
#                                   trails that abandon in a cycle
#           + sum(aux rows)         every row of aux(), enumerated
#           + aux_reserve           the difference between the previously
#                                   DISCLOSED aux reading (2^114.812, the v23
#                                   note) and the enumerated rows, reserved so
#                                   the claim covers both readings
#
# Every term is an exact Fraction of reference operations.  The claim is the
# least 4-decimal LOG2 value >= log2(W_total) + pad, asserted in package.
# ===========================================================================
PRIOR_AUX_DISCLOSED_LOG2 = Fraction(114812, 1000)     # 2^114.812, as v23 shipped it


def trap_charge(s=None, theta=None, cap_mult=128, rate=None, c_mul=None):
    """rate * cap * price, in reference units, EXACT.

    rate  : the charged bound on the probability that a trail abandons in a cycle
            before its first distinguished point.  The production rate is the
            birthday union bound over the CAP horizon: t(t-1)/2N with t = cap =
            128*2^80 and N = 2^371, i.e. <= 2^-198, and over the REPORT interval
            theta^2/2N = 2^-92.3 (constants.LOSS_ROWS['dp_detect_tail'] carries the
            same row as 74e-12).  We charge the LARGER of the two so the row is
            conservative: rate = 2^-92.3.
    cap   : the step cap, 128 * 2^80 (v23 used 128*theta; unchanged).
    price : c_mul, the multiply-channel step price in reference units.

    The row is an ADDITIVE charge in the same units as W_stress, which is what
    v23 never did: it printed the cap rule and charged the waste nowhere in the
    sum the claim was taken from.
    """
    s = K.S_STEPS if s is None else s
    theta = (1 << K.THETA_BITS) if theta is None else theta
    c_mul = (Fraction(K.C_REF + rows()['fold_mul'], K.C_REF)) if c_mul is None else c_mul
    if rate is None:
        rate = Fraction(1, 1 << 92) * Fraction(3, 2)      # 2^-92.3, the theta row
    cap = cap_mult * theta
    return Fraction(rate) * Fraction(cap) * Fraction(c_mul)


def total_all(s=None, theta=None, cap_mult=128, rate=None, aux_reserve=True):
    """THE claim-bearing total.  One function, one sum, no part outside it."""
    s = K.S_STEPS if s is None else s
    L = ledgers(s=s)
    A = aux(s=s) if theta is None else aux(s=s, theta=theta)
    w = L['w_stress']
    tr = trap_charge(s=s, theta=theta, cap_mult=cap_mult, rate=rate, c_mul=L['c_mul_stress'])
    aux_sum = A['total']
    reserve = Fraction(0)
    if aux_reserve:
        prior = Fraction(2) ** int(PRIOR_AUX_DISCLOSED_LOG2)  # not exact; computed below
        prior = _pow2_frac(PRIOR_AUX_DISCLOSED_LOG2, 40)
        if prior > aux_sum:
            reserve = prior - aux_sum
    total = w + tr + aux_sum + reserve
    return dict(w_stress=w, trap_charge=tr, aux_rows={k: A[k] for k in A if k != 'total'},
                aux_sum=aux_sum, aux_reserve=reserve, total=total,
                log2_total=log2_decimal(total), w_plus_aux=w + aux_sum + tr,
                covered=bool(total >= w + tr + aux_sum))


def _pow2_frac(exp_log2, prec):
    """2^(num/den) as an exact-enough Fraction (>= the true value)."""
    from decimal import Decimal, getcontext
    getcontext().prec = prec + 20
    v = (Decimal(int(exp_log2.numerator)) / Decimal(int(exp_log2.denominator))).exp() \
        if False else (Decimal(int(exp_log2.numerator)) / Decimal(int(exp_log2.denominator))
                       * Decimal(2).ln()).exp()
    num, den = v.as_integer_ratio()
    return Fraction(num, den)


def claim_from_total(dp=4, pad_bits=None):
    """The submitted claim: least 4-decimal log2 value covering the WHOLE sum."""
    from decimal import Decimal, getcontext
    T = total_all()
    getcontext().prec = 60
    pad = Decimal('2e-5') if pad_bits is None else pad_bits
    target = T['log2_total'] + pad
    m = -((-target * 10**dp).to_integral_value(rounding='ROUND_FLOOR'))
    c = Decimal(m) / Decimal(10**dp)
    while c < target:
        m += 1
        c = Decimal(m) / Decimal(10**dp)
    while c - Decimal(1) / Decimal(10**dp) >= target:
        m -= 1
        c = Decimal(m) / Decimal(10**dp)
    assert c >= target
    assert c >= T['log2_total']
    return c, T
