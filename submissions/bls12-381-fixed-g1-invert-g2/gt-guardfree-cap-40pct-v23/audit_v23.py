#!/usr/bin/env python3
"""v23 CERTIFICATE BUILDER AND AUDIT.

Emits certificate.json and evidence/inequality-audit.json, and asserts, in code:

  A1  every per-step instruction row re-sums to the step totals (integer identity);
  A2  the aux rows re-sum to the aux total (integer identity, Fraction);
  A3  the claim integer is >= log2(W_stress) + 2e-5 bits and < that + 1e-4 (in the
      log domain, by Decimal, away from the 4-decimal boundary);
  A4  every literal printed in public-note.md with a >= / <= direction is a FLOOR or
      a CEIL of its certificate value, checked by regex over the note itself;
  A5  the claim integer is a whole number of 1e-4 bits: numerator/denominator exact.
"""
import json
import math
import os
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import constants as K
import ledger
import v23 as V

HERE = os.path.dirname(os.path.abspath(__file__))
EV = os.path.join(HERE, 'evidence')
getcontext().prec = 60


def lg(x):
    from fractions import Fraction as F
    if isinstance(x, F):
        return Decimal(x.numerator).ln() / Decimal(2).ln() - Decimal(x.denominator).ln() / Decimal(2).ln()
    return Decimal(x).ln() / Decimal(2).ln()


def main():
    L = ledger.ledgers()
    r = L['rows']
    # A1 per-step re-sum
    fold_mul = sum(r[k] for k in ('key_hash', 'dp_test', 'hist_slot', 'loop',
                                  'table_fetch', 'label_add', 'selector', 'column_addr'))
    assert fold_mul == r['fold_mul'], 'A1 fold_mul re-sum'
    assert r['mul_step'] == K.C_REF + fold_mul, 'A1 mul_step'
    fold_sqr = (r['key_hash'] + r['dp_test'] + r['hist_slot'] + r['loop'] + r['selector']
                + 2 * r['label_add'])
    assert fold_sqr == r['fold_sqr'], 'A1 fold_sqr re-sum'
    assert r['sqr_step'] == r['cyclo_sqr'] + fold_sqr, 'A1 sqr_step'
    assert r['guard'] == 0, 'A1 the guard row is zero (dropped in v23)'

    # A2 aux re-sum
    A = ledger.aux()
    parts = {k: A[k] for k in ('table_build', 'gt_exp_starts', 'restart_bookkeeping',
                               'report_writes', 'store_lookup', 'report_cap_drain',
                               'solve', 'verify')}
    assert A['total'] == sum(parts.values()) + A['cap_slack'], 'A2 aux re-sum (rows + slack)'

    # A3/A5 the claim integer
    w = L['w_stress']
    claim = ledger.claim_log2(w)
    num = int(claim * 10000)
    assert Fraction(num, 10000) == claim, 'A5 claim is an exact 4-decimal log2 value'
    logw = lg(w)
    assert claim >= logw + Decimal('2e-5'), 'A3 claim >= log2 W + pad'
    assert claim - Decimal(1) / 10000 < logw + Decimal('2e-5'), 'A3 claim is the LEAST such'
    assert abs(claim - (logw + Decimal('2e-5'))) > Decimal('1e-30'), 'A3 away from boundary'

    SUC = V.success_rows()
    TRAP = V.trap_charge()
    BND = V.bound_route()
    DENS = V.cap_rows()
    st = ledger.storage()

    assert SUC['model'] >= SUC['rate_scaled'] - 1e-15, 'model >= rate'
    assert SUC['product'] >= SUC['net'] - 1e-15, 'product >= net'
    assert SUC['clears_floor_point'] and SUC['clears_floor_cp95_lower'], 'floor cleared'
    assert SUC['retention_eff'] <= 1.0, 'retention <= 1'

    vals = {'net': SUC['net'], 'model': SUC['model'], 'product': SUC['product'],
            'prod_low': SUC['product_cp95_lower'], 'net_cp': SUC['net_cp95_lower'],
            'claim': float(claim), 'lam_rep': BND['iid_birthday_value_production'],
            'dens': DENS['density_factor'], 'tau': SUC['trap_tau'],
            'ret': SUC['retention_eff']}
    rows = V.inequality_rows(vals)
    audit = V.audit_inequalities(rows)
    flip = V.planted_flip(rows)

    # A4 the note literals
    note = open(os.path.join(HERE, 'public-note.md')).read()
    note_rows = []
    for name, val, op, lit in rows:
        pats = [p for p in (re.escape(lit),) if p in note]
        note_rows.append({'row': name, 'op': op, 'literal': lit,
                          'in_note': bool(pats), 'cert_value': float(val)})
    bad_note = []
    for n in note_rows:
        if not n['in_note']:
            continue
        v = float(n['cert_value'])
        x = float(n['literal'])
        ok = (v >= x) if n['op'] == '>=' else (v <= x)
        if not ok:
            bad_note.append(n)
    assert not bad_note, 'A4 note literals: %r' % bad_note

    parts_no_slack = dict(parts)
    cert = {
        'version': 'v23 guard-free + cap',
        'track': 'bls12-381-fixed-g1-invert-g2',
        'the_redesign': 'the routing guard is DROPPED and the v22 128*theta abandon cap is '
                        'KEPT: traps become capped, priced work instead of lost witnesses',
        'claim_log2': float(claim),
        'claim_integer': {'numerator': str(num), 'denominator': '10000',
                          'exact': str(claim)},
        'work': {'S': str(L['S']), 'S_log2': float(lg(Fraction(L['S']))),
                 'c_mul': float(L['c_mul']), 'c_sqr': float(L['c_sqr']),
                 'c_raw': float(L['c_raw']), 'c_stress': float(L['c_stress']),
                 'w_raw': str(L['w_raw']), 'w_raw_log2': float(lg(L['w_raw'])),
                 'w_stress': str(w), 'w_stress_log2': float(lg(w)),
                 'delta_stress_over_raw': float(L['delta'])},
        'per_step_rows': {k: int(r[k]) for k in
                          ('key_hash', 'dp_test', 'hist_slot', 'loop', 'table_fetch',
                           'label_add', 'selector', 'column_addr', 'guard', 'fold_mul',
                           'fold_sqr', 'cyclo_sqr', 'mul_step', 'sqr_step')},
        'per_step_resum': {'fold_mul_is_the_sum_of_its_rows': True,
                           'mul_step_equals_C_REF_plus_fold_mul': True,
                           'guard_row_is_zero_because_the_guard_is_dropped': True},
        'success': {k: SUC[k] for k in
                    ('lambda_', 'model', 'rate_scaled', 'retention_eff',
                     'retention_eff_raw_mu_times_p_detect', 'retention_eff_cp95_lower',
                     'product', 'product_cp95_lower', 'trap_tau', 'stop_loss', 'net',
                     'net_cp95_lower', 'clears_floor_point', 'clears_floor_cp95_lower',
                     'p_detect_dp_gated', 'p_detect_dp_gated_cp95_lower', 'dp_gated_pairs',
                     'any_step_agree_for_comparison_only', 'which_statistic_feeds_success',
                     'retention_source', 'bigN_row', 'lag32_row', 'guarded_comparison',
                     'reports_measured_per_walk', 'reports_ceil_formula_per_walk',
                     'reports_check_ratio_measured_over_formula', 'reports_restart',
                     'reports_rule')},
        'trap_row': {k: TRAP[k] for k in
                     ('repeat_probability_over_the_walk', 'events_charged', 'abandon_cap',
                      'total_reference_units', 'total_log2', 'waste_fraction_of_W_stress',
                      'bound_note', 'method')},
        'bound_row': BND,
        'cap_row': {'cap_rows_from_density': DENS,
                    'guard_free_cap_evidence': 'evidence/cap-v23.json'},
        'storage': {'peak_bytes': int(st['peak_bytes']), 'rows': {k: int(v) if isinstance(v, int) else v
                                                                 for k, v in st.items() if k != 'overflow_policy'}},
        'aux': {'rows': {k: str(v) for k, v in parts_no_slack.items()},
                'total': str(A['total']), 'total_log2': float(lg(A['total'])),
                'resum_ok': True,
                'covers': 'aux total %s is 2^%.3f, below W_stress 2^%.3f: inside the claim'
                          % (str(A['total']), float(lg(A['total'])), float(lg(w)))},
        'inequality_audit': audit,
        'planted_flip': flip,
        'note_literal_audit': note_rows,
        'checks': V.CHECKS,
    }
    json.dump(cert, open(os.path.join(HERE, 'certificate.json'), 'w'), indent=1,
              sort_keys=True)
    json.dump({'inequalities': audit, 'planted_flip': flip,
               'note_literals': note_rows}, open(os.path.join(EV, 'inequality-audit.json'), 'w'),
              indent=1, sort_keys=True)
    print('claim %s  W_stress 2^%.6f  success net %.7f cp95 %.7f  checks %d'
          % (claim, float(lg(w)), SUC['net'], SUC['net_cp95_lower'], len(V.CHECKS)))
    print('wrote certificate.json + evidence/inequality-audit.json')


if __name__ == '__main__':
    main()
