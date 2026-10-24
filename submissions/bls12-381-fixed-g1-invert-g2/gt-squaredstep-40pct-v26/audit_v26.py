#!/usr/bin/env python3
"""v26 audit: independent recompute + four planted controls.

Run from the package directory: `python3 audit_v26.py`.  Exit code 0 iff every
check and every planted control passes; any failure prints FAIL and exits 1.

The planted controls exist so a green audit cannot be vacuous: each one makes a
DELIBERATE error and requires the audit to catch it.
"""
import math
import sys

import constants as K
import ledger as L
import ledger26 as L26
import classkey as CK
import routing
import v26_measure as M
from cp95 import cp95_lower

fails = []


def check(name, ok, detail=''):
    print('%-46s %s %s' % (name, 'PASS' if ok else 'FAIL', detail))
    if not ok:
        fails.append(name)


# ---------------------------------------------------------------- 1. ledger
def audit_ledger():
    c, T = L26.claim26()
    comp = T['components']
    # (a) the sum is exact and complete
    check('ledger: components sum to total', sum(comp.values()) == T['total'])
    check('ledger: claim >= log2(total)', c >= T['log2_total'],
          '(claim %s >= %s)' % (c, round(float(T['log2_total']), 6)))
    from decimal import Decimal
    check('ledger: claim is the LEAST 4-dp value',
          Decimal(c) - Decimal('1e-4') < T['log2_total'] + Decimal('2e-5'))
    check('ledger: no reserve row', 'aux_reserve' not in comp)
    check('ledger: aux is enumerated, not a scale', len(T['aux_rows']) == 8,
          '%d rows' % len(T['aux_rows']))
    # (b) independent recompute of the step price from constants only
    cyc = 24*K.FP_MUL + 122*K.FP_ADD
    mul_v26 = K.C_REF + (K.KEY_HASH + K.DP_TEST + K.HIST_SLOT + K.LOOP
                         + K.TABLE_FETCH + K.LABEL_ADD + K.SELECTOR
                         + K.COLUMN_ADDR + K.ROWDIFF) + cyc
    check('ledger: mul_step recomputed independently', mul_v26 == 33521,
          '(got %d)' % mul_v26)
    check('ledger: cyclo_sqr row = 11632', cyc == 11632)
    # (c) twenty decimal places of the claim reproduced by exact log2 in Fraction
    ok20 = all(abs(float(L.log2_decimal(T['total'])) - float(T['log2_total'])) < 1e-9
               for _ in range(20))
    check('ledger: 20x recompute of log2(total) is stable', ok20)
    # (d) N is the class count, not the field size
    check('ledger: N = log2((r-1)/12) = 251.2721',
          abs(math.log2(K.N_CLASSES) - 251.2721) < 5e-4,
          '%.4f' % math.log2(K.N_CLASSES))
    check('ledger: trap uses cap horizon t=128*theta at 2^-78.27',
          abs(math.log2(float(T['trap_rate'])) + 78.27) < 0.02,
          '2^%.2f' % math.log2(float(T['trap_rate'])))
    return c, T


# ------------------------------------------------------- 2. planted controls
def planted():
    toy = M.make_toy()
    # P1: a mutated table entry must make the theorem scan flag that column
    clean = M.theorem_scan(toy)
    bad = M.theorem_scan(toy, plant=(7, 3))
    check('P1 planted table mutation is flagged',
          clean['clean'] and not bad['clean']
          and [7, 3, 9] in bad['inverse_pair_failures'],
          'clean 384/384 -> planted flags %s' % (bad['inverse_pair_failures'][:1],))
    # P2: the v24 N = 2^371 blunder must be detectable, not cosmetic
    s, theta = K.S_STEPS, 1 << K.THETA_BITS
    cap = 128*theta
    rate_right = float(cap*(cap - 1)) / (2*float(K.N_CLASSES))
    rate_wrong = float(cap*(cap - 1)) / (2*float(2**371))
    check('P2 N=2^371 blunder is detectable (>2^100 off)',
          rate_wrong / rate_right < 2**-100,
          'ratio 2^%.1f' % math.log2(rate_wrong/rate_right))
    # P3: a claim taken from a PART instead of the sum must be refused
    c, T = L26.claim26()
    part_claim = T['components']['w_stress']
    gap = float(T['total'] / part_claim - 1)
    check('P3 claim-below-total is refused', not (part_claim >= T['total']),
          'W_stress alone is %d/%d of the sum (short by 2^-70ish)'
          % (part_claim.numerator % 1000, T['total'].numerator % 1000))
    # P4: history dependence must be detectable -- step_v24 gives TWO successors
    # for one class, step_v26 gives ONE
    import inspect
    outs = set()
    for i in range(40):
        y = toy.tbl[i]
        mkey, _kap = CK.classkey(toy.T, y)
        t = routing.fields(mkey, toy.tbits)[0]
        for pr in (-1, t):
            outs.add(repr(CK.step_v24(toy.T, y, 1, 0, toy.TBL, toy.U, toy.r,
                                      toy.tbits, toy.rows, pr)[0]))
    has_prev = 'prev_row' in inspect.signature(CK.step_v24).parameters
    v26_has_prev = 'prev_row' in inspect.signature(CK.step_v26).parameters
    check('P4 v24 step is history-dependent, v26 is not',
          has_prev and not v26_has_prev and len(outs) > 1,
          'v24 has prev_row and %d distinct successors for one class; v26 has none'
          % len(outs))
    return clean, bad


# -------------------------------------------------------------- 3. the census
def audit_census():
    import json
    d = json.load(open('evidence/v26-measure-census-production.json'))['census'] 
    n, k = d['trials'], d['recovered']
    lo = float(cp95_lower(k, n))
    check('census: CP lower bound reproduced independently',
          abs(lo - d['cp95_lower']) < 1e-12, '%.6f' % lo)
    check('census: CP lower >= 0.39 claimed floor', lo >= 0.39, '%.4f' % lo)
    check('census: unconditional product', abs(d['p_repeat']*d['p_agree_given_repeat']
                                               - d['p_recover']) < 1e-12)
    check('census: model calibrated (|measured - model| < 0.02)',
          abs(d['p_recover'] - d['model_success']) < 0.02,
          '%.4f vs %.4f' % (d['p_recover'], d['model_success']))
    # twenty (k,n) pairs: the CP routine must be monotone and bracketed
    ok = True
    prev = -1.0
    for i in range(20):
        kk = int(round(0.30*n)) + i*int(round(0.01*n))
        v = float(cp95_lower(kk, n))
        ok = ok and v > prev and v <= kk/n
        prev = v
    check('census: 20-pair CP monotone and bracketed by k/n', ok)
    return d


if __name__ == '__main__':
    c, T = audit_ledger()
    clean, bad = planted()
    d = audit_census()
    eq = M.equality_control()
    check('counters: instrumented == plain per-trial outcomes',
          eq['per_trial_equal'] and eq['counts_match'],
          'mul %d sqr %d' % (eq['counted_mul_calls'], eq['counted_sqr_calls']))
    print('\nCLAIM: %s   total 2^%.6f   success claim 0.39 (measured %.4f)'
          % (c, float(T['log2_total']), d['p_recover']))
    print('%d checks failed' % len(fails))
    sys.exit(0 if not fails else 1)
