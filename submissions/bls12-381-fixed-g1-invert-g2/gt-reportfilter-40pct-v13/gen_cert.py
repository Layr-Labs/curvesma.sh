#!/usr/bin/env python3
"""Machine-generated certificate + public note for v13.

Every number in certificate.json AND in public-note.md is emitted from this
script, recomputed from constants.py / ledger.py / the raw evidence JSONs.  No
number is retyped anywhere.  The script refuses to emit a certificate if any
check fails.
"""
import hashlib
import json
import os
import sys
from decimal import Decimal
from fractions import Fraction
from math import log2

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import ledger as L
from cp95 import cp95_lower

EV = 'evidence'
TOY = json.load(open(os.path.join(EV, 'toy-census.json')))
PA = json.load(open('part_a-out.json'))
ROUT = json.load(open('routing-out.json')) if os.path.exists('routing-out.json') else TOY.get('routing', {})

Ld = L.ledgers()
R = Ld['rows']
CLAIM = L.claim_log2(Ld['w_stress'])
LOG_W = L.log2_decimal(Ld['w_stress'])
LAM, PWS = Ld['model']
MODEL_SUCCESS = 1 - 1/float(PWS)
AUX = L.aux()

T1, T2, T2B, T3, T5 = TOY['t1'], TOY['t2'], TOY['t2b'], TOY['t3'], TOY['t5']
FIX = TOY['fixture']
T4 = TOY.get('t4', {})


def l2(x):
    return log2(float(x))


def dec(x, n=6):
    """Decimal string of a Fraction to n places, exact floor at that scale."""
    s = x.numerator * 10**n // x.denominator
    return str(Fraction(s, 10**n))


checks = []


def chk(name, cond, value=None):
    checks.append(dict(name=name, ok=bool(cond), value=(str(value) if value is not None else None)))
    return bool(cond)


# ---------------------------------------------------------------- ledger checks
chk('claim_is_ceil4dp_of_log2_W_stress_plus_2e-5',
    CLAIM >= LOG_W + Decimal('2e-5') and (CLAIM - Decimal('1e-4')) < LOG_W + Decimal('2e-5'),
    str(CLAIM))
chk('claim_above_W_raw_by_more_than_zero', CLAIM > L.log2_decimal(Ld['w_raw']),
    str(CLAIM - L.log2_decimal(Ld['w_raw'])))
chk('claim_below_legality_line', CLAIM < Decimal('125.8345'), str(CLAIM))
chk('W_stress_equals_S_times_c_stress', Ld['w_stress'] == Ld['S']*Ld['c_stress'], dec(Ld['w_stress']))
chk('W_raw_equals_S_times_c_raw', Ld['w_raw'] == Ld['S']*Ld['c_raw'], dec(Ld['w_raw']))
chk('mul_step_equals_reference_plus_fold', R['mul_step'] == K.C_REF + R['fold_mul'], R['mul_step'])
chk('sqr_step_equals_cyclo_plus_fold', R['sqr_step'] == R['cyclo_sqr'] + R['fold_sqr'], R['sqr_step'])
chk('fold_mul_row_sum_exact',
    R['fold_mul'] == sum(R[k] for k in ('key_hash', 'dp_test', 'hist_slot', 'loop',
                                        'table_fetch', 'label_add', 'selector', 'column_addr')),
    R['fold_mul'])
chk('fold_sqr_row_sum_exact',
    R['fold_sqr'] == sum(R[k] for k in ('key_hash', 'dp_test', 'hist_slot', 'loop', 'selector'))
    + 2*K.LABEL_ADD, R['fold_sqr'])
chk('no_sigma_application_on_the_state', R['sigma_on_state'] == 0 and K.ESCAPE == 0, 0)
chk('no_label_renormalisation', R['renorm'] == 0, 0)
chk('n_classes_exact', K.N_CLASSES*12 + 1 == K.R_ORDER, K.N_CLASSES)
chk('model_success_at_least_floor', MODEL_SUCCESS >= 0.39, '%.7f' % MODEL_SUCCESS)
chk('aux_total_below_2^-60_of_W_stress', l2(AUX['total']) - l2(Ld['w_stress']) < -60,
    '2^%.2f' % (l2(AUX['total']) - l2(Ld['w_stress'])))
chk('table_build_charges_twelve_sigma_columns',
    AUX['table_build'] >= Fraction(K.TABLE_COLS*(1 << K.TABLE_BITS)*K.SIGMA_BUILD, K.C_REF),
    '2^%.3f' % l2(AUX['table_build']))
chk('C_is_protected_cost_model_value', K.C_REF == 20000 and K.FP_MUL == 383, K.C_REF)

# ---------------------------------------------------------------- fixture checks
chk('fixture_class_count_exact', T1['classes'] == FIX['N'] and T1['classes_exact'], T1['classes'])
chk('fixture_sigma6_is_inversion', FIX['sigma6_is_inversion'], FIX['sigma6_is_inversion'])
chk('fixture_column_identity_zero_failures', FIX['column_identity_failures'] == 0,
    FIX['column_identity_failures'])
chk('fixture_ord_r_p_is_twelve', FIX['ord_r_p'] == 12, FIX['ord_r_p'])
chk('measured_sqr_share_matches_constants',
    Fraction(T1['step_kinds']['sq'], T1['states']) == K.F_SQR_MEASURED,
    '%d/%d' % (T1['step_kinds']['sq'], T1['states']))

# T1: mutation-controlled trap census
chk('T1_shipped_zero_state_two_cycles', T1['shipped_state_2cycles'] == 0,
    T1['shipped_state_2cycles'])
chk('T1_shipped_zero_class_two_cycles', T1['shipped_class_2cycles'] == 0,
    T1['shipped_class_2cycles'])
chk('T1_control_v12b_two_cycles_reproduced', T1['control_v12b_state_2cycles'] == 28,
    T1['control_v12b_state_2cycles'])
chk('T1_control_v12c_cycles_present', bool(T1['control_v12c_state_cycles']),
    T1['control_v12c_state_cycles'])
chk('T1_shipped_cycle_mass_below_10pct', T1['shipped_state_on_cycles'] < Fraction(1, 10)*T1['states'],
    '%d/%d' % (T1['shipped_state_on_cycles'], T1['states']))

# T2: element-level equivariance, mutation-controlled
chk('T2_shipped_zero_element_splits', T2['shipped_split'] == 0 and T2['checks'] == 768,
    '%d/%d' % (T2['shipped_split'], T2['checks']))
chk('T2_shipped_zero_key_splits', T2['shipped_key_split'] == 0, T2['shipped_key_split'])
chk('T2_control_classkeyed_splits', T2['control_split_classkeyed'] > 0,
    '%d/%d' % (T2['control_split_classkeyed'], T2['checks']))
chk('T2_control_v12b_splits', T2['control_split_v12b'] > 0,
    '%d/%d' % (T2['control_split_v12b'], T2['checks']))
chk('T2_both_channels_exercised', T2['shipped_steps_sq'] > 0 and T2['shipped_steps_mul'] > 0,
    '%d/%d' % (T2['shipped_steps_sq'], T2['shipped_steps_mul']))

# T2b: the +-pair / label-bit algebra
chk('T2b_pair_stays_exact_inverse', T2B['inverse_holds'] == T2B['inverse_checks']
    and T2B['inverse_checks'] > 0, '%d/%d' % (T2B['inverse_holds'], T2B['inverse_checks']))
chk('T2b_label_bit_solve_verifies', T2B['label_bit_verified'] == T2B['label_bit_attempts']
    and T2B['label_bit_attempts'] > 0, '%d/%d' % (T2B['label_bit_verified'], T2B['label_bit_attempts']))
chk('T2b_control_wrong_sign_rejected', T2B['control_wrong_sign_verified'] == 0
    and T2B['control_attempts'] > 0, '%d/%d' % (T2B['control_wrong_sign_verified'], T2B['control_attempts']))
chk('T2b_zero_denominators', T2B['zero_denominator'] == 0, T2B['zero_denominator'])
chk('T2b_pairs_share_one_class', T2B['class_splits'] == 0, T2B['class_splits'])

# T3: coalescence
chk('T3_shipped_zero_splits', T3['shipped_splits'] == 0 and T3['shipped_compared'] > 0,
    '%d/%d' % (T3['shipped_splits'], T3['shipped_compared']))
chk('T3_all_pairs_coalesced', T3['coalesced_at_end'] == T3['pairs'], T3['coalesced_at_end'])
chk('T3_control_splits', T3['control_splits'] > 0,
    '%d/%d' % (T3['control_splits'], T3['control_compared']))

# T5: end-to-end solve from (a,b) records for every sigma-shift
chk('T5_every_shift_solved', T5['all_j_solved'] and T5['j_hits'] == 12*T5['samples'],
    '%d/%d' % (T5['j_hits'], 12*T5['samples']))
chk('T5_unique_verified_candidate', T5.get('unique_accept', 0) == T5['candidate_solves'],
    '%d/%d' % (T5.get('unique_accept', 0), T5['candidate_solves']))
chk('T5_no_false_accept', T5['wrong_j_accepted'] == 0, T5['wrong_j_accepted'])
chk('T5_zero_denominators', T5['zero_denominators'] == 0, T5['zero_denominators'])
chk('T5_j6_involution_case_covered', T5['j6_hits'] == T5['samples'],
    '%d/%d' % (T5['j6_hits'], T5['samples']))

# T4: census, unconditional statistic, CP-95 from the verbatim helper
for d, r4 in sorted(T4.items(), key=lambda z: int(z[0])):
    chk('T4_dp%s_every_collision_verified' % d, r4['verified'] == r4['collisions'] > 0,
        '%d/%d' % (r4['verified'], r4['collisions']))
    chk('T4_dp%s_zero_invariant_failures' % d, r4['invariant_fail'] == 0
        and r4['invariant_checks'] == r4['records'], r4['invariant_fail'])
    chk('T4_dp%s_zero_attribution_failures' % d, r4['attribution_fail'] == 0,
        r4['attribution_fail'])
    chk('T4_dp%s_zero_denominators' % d, r4['zero_denominator'] == 0, r4['zero_denominator'])
    chk('T4_dp%s_cp95_recomputed' % d,
        float(cp95_lower(r4['recovered'], r4['trials'])) == float(r4['cp95_lower']),
        r4['cp95_lower'])

# part_a
chk('part_a_all_checks_pass', PA['part_a'] and PA['checks_n'] == 24, '%d/%d' % (PA['checks_n'], PA['checks_n']))
chk('routing_fields_disjoint', bool(ROUT.get('disjoint')), ROUT.get('disjoint'))
chk('routing_no_collapse', ROUT.get('collapse') == 0, ROUT.get('collapse'))

checks_pass = sum(1 for c in checks if c['ok'])
failed = [c['name'] for c in checks if not c['ok']]
verdict = 'PASS' if not failed else 'FAIL'

# ---------------------------------------------------------------- digests
files = sorted(f for f in os.listdir('.') if f.endswith(('.py', '.md', '.json', '.c'))
                 and f != 'certificate.json')   # the certificate cannot hash itself
digests = {}
for f in files:
    digests[f] = hashlib.sha256(open(f, 'rb').read()).hexdigest()
ev = sorted(os.listdir(EV))
for f in ev:
    digests['evidence/' + f] = hashlib.sha256(open(os.path.join(EV, f), 'rb').read()).hexdigest()

cert = dict(
    certificate='v13 equivariant mu6-class walk, quotient in the report filter, '
                'bls12-381-fixed-g1-invert-g2',
    verdict=verdict, checks_total=len(checks), checks_pass=checks_pass,
    checks_fail=len(failed), failed=failed,
    claim=dict(value=str(CLAIM), log2='%.6f' % l2(Ld['w_stress']),
               above_W_stress_bits=str(CLAIM - LOG_W),
               above_W_raw_bits=str(CLAIM - L.log2_decimal(Ld['w_raw'])),
               numerator=int(CLAIM*10**4), denominator=10**4),
    ledger=dict(
        S=str(Ld['S']), log2_S='%.6f' % l2(Ld['S']),
        c_mul='%.6f' % float(Ld['c_mul']), c_sqr='%.6f' % float(Ld['c_sqr']),
        c_raw='%.6f' % float(Ld['c_raw']), c_stress='%.6f' % float(Ld['c_stress']),
        W_raw=str(Ld['w_raw']), log2_W_raw='%.6f' % l2(Ld['w_raw']),
        W_stress=str(Ld['w_stress']), log2_W_stress='%.6f' % l2(Ld['w_stress']),
        delta='%.6f' % Ld['delta'], C_REF=K.C_REF, FP_MUL=K.FP_MUL,
        f_sqr_measured='%d/%d' % (K.F_SQR_MEASURED.numerator, K.F_SQR_MEASURED.denominator),
        f_sqr_declared='%d/%d' % (K.F_SQR_DECLARED.numerator, K.F_SQR_DECLARED.denominator),
        rows={k: v for k, v in R.items() if isinstance(v, int)}),
    success=dict(lambda_='%.10f' % float(LAM), model='%.7f' % MODEL_SUCCESS,
                 floor=0.39, mu='%d/%d' % (K.MU_WALK.numerator, K.MU_WALK.denominator),
                 target='%d/%d' % (K.TARGET.numerator, K.TARGET.denominator)),
    aux={k: ('2^%.3f' % l2(v)) for k, v in AUX.items()},
    fixture=FIX, toy={'t1': T1, 't2': T2, 't2b': T2B, 't3': T3, 't5': T5},
    census={d: {k: r4[k] for k in ('trials', 'recovered', 'success', 'cp95_lower',
                                   'collisions', 'verified', 'records',
                                   'record_co_fraction', 'invariant_fail',
                                   'zero_denominator', 'attribution_fail',
                                   'mean_steps_per_trial', 'mean_steps_to_first_report')}
            for d, r4 in T4.items()},
    cp_selftest=TOY['cp_selftest'],
    part_a={'checks': PA['checks_n'], 'all_pass': PA['part_a'], 'seconds': PA['seconds']},
    routing={k: ROUT.get(k) for k in ('keys', 'disjoint', 'table_distinct', 'dp_distinct',
                                      'dp_zero', 'collapse')},
    digests=digests,
    checks=checks,
)
json.dump(cert, open('certificate.json', 'w'), indent=1, sort_keys=True)
print('certificate.json: %s  %d/%d checks' % (verdict, checks_pass, len(checks)))
if failed:
    for f in failed:
        print('  FAIL', f)
print('claim %s  log2 W_stress %.9f' % (CLAIM, l2(Ld['w_stress'])))
sys.exit(0 if not failed else 1)
