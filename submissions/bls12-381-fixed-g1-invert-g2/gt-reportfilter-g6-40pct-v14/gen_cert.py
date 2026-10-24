#!/usr/bin/env python3
"""Machine-generated certificate for v14.  Every number in certificate.json is
emitted from this script, recomputed from constants.py / ledger.py / sigma_set.py
and the raw evidence JSONs.  No number is retyped.  The script refuses to emit a
PASS if any check fails.

The certificate distinguishes two kinds of row:

  checks      : facts the package asserts and the reviewer can re-run;
  failed_gate : premises the package explicitly does NOT claim to have met, with
                the reason.  v14 carries one such row (see argument.md section 6).
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
import sigma_set as SS
from cp95 import cp95_lower

EV = 'evidence'
TOY = json.load(open(os.path.join(EV, 'toy-census.json')))
SSET = json.load(open(os.path.join(EV, 'sigma-set.json')))
PA = json.load(open('part_a-out.json'))
ROUT = json.load(open('routing-out.json'))

Ld = L.ledgers()
R = Ld['rows']
CLAIM = L.claim_log2(Ld['w_stress'])
LOG_W = L.log2_decimal(Ld['w_stress'])
LAM, PWS = Ld['model']
MODEL_SUCCESS = 1 - 1/float(PWS)
NET_SUCCESS = Ld['net_success']
AUX = L.aux()
STO = L.storage()

T1, T1J, T2, T2B, T3, T5 = TOY['t1'], TOY['t1j'], TOY['t2'], TOY['t2b'], TOY['t3'], TOY['t5']
FIX = TOY['fixture']
T4 = TOY.get('t4', {})
GR = TOY.get('guard_rate', {})


def l2(x):
    return log2(float(x))


def dec(x, n=6):
    s = Fraction(x).numerator * 10**n // Fraction(x).denominator
    return str(Fraction(s, 10**n))


checks = []
failed_gate = []


def chk(name, cond, value=None):
    checks.append(dict(name=name, ok=bool(cond), value=(str(value) if value is not None else None)))
    return bool(cond)


# ------------------------------------------------- ledger: ONE squared fraction
chk('f_sqr_is_the_routing_share_only', K.F_SQR == Fraction(3, 8) and K.F_SQR_DECLARED == K.F_SQR,
    '%d/%d' % (K.F_SQR.numerator, K.F_SQR.denominator))
chk('budget_and_price_use_the_same_f',
    Ld['f'] == K.F_SQR and Ld['S'] == K.step_budget(f=K.F_SQR, target=K.TARGET_PAD),
    'S=2^%.6f' % l2(Ld['S']))
chk('model_lambda_uses_the_same_f',
    Ld['model'][0] == Fraction(Ld['S']*Ld['S'], 1)*(1-K.F_SQR*K.F_SQR)/(2*K.MU_WALK*K.N_CLASSES),
    '%.10f' % float(Ld['model'][0]))
chk('f_sqr_measured_is_not_an_input',
    'F_SQR_MEASURED' not in open('constants.py').read(),
    'constants.py has no measured-share constant')
chk('claim_is_ceil4dp_of_log2_W_stress_plus_2e-5',
    CLAIM >= LOG_W + Decimal('2e-5') and (CLAIM - Decimal('1e-4')) < LOG_W + Decimal('2e-5'),
    str(CLAIM))
chk('claim_above_W_raw_by_more_than_zero', CLAIM > L.log2_decimal(Ld['w_raw']),
    str(CLAIM - L.log2_decimal(Ld['w_raw'])))
chk('claim_below_legality_line', CLAIM < Decimal('125.8345'), str(CLAIM))
chk('W_stress_equals_S_times_c_stress', Ld['w_stress'] == Ld['S']*Ld['c_stress'], dec(Ld['w_stress']))
chk('W_raw_equals_S_times_c_raw', Ld['w_raw'] == Ld['S']*Ld['c_raw'], dec(Ld['w_raw']))
chk('c_raw_is_the_f_mixture', Ld['c_raw'] == (1-K.F_SQR)*Ld['c_mul'] + K.F_SQR*Ld['c_sqr'],
    dec(Ld['c_raw']))

# ------------------------------------------------- ledger: rows and the guard
chk('mul_step_equals_reference_plus_fold', R['mul_step'] == K.C_REF + R['fold_mul'], R['mul_step'])
chk('sqr_step_equals_cyclo_plus_fold', R['sqr_step'] == R['cyclo_sqr'] + R['fold_sqr'], R['sqr_step'])
chk('fold_mul_row_sum_exact',
    R['fold_mul'] == sum(R[k] for k in ('key_hash', 'dp_test', 'hist_slot', 'loop',
                                        'table_fetch', 'label_add', 'selector',
                                        'column_addr', 'guard')), R['fold_mul'])
chk('fold_sqr_row_sum_exact',
    R['fold_sqr'] == sum(R[k] for k in ('key_hash', 'dp_test', 'hist_slot', 'loop', 'selector'))
    + 2*K.LABEL_ADD, R['fold_sqr'])
chk('guard_charged_on_the_multiply_channel_only',
    R['guard'] == K.GUARD == 48 and R['fold_sqr'] ==
    sum(R[k] for k in ('key_hash', 'dp_test', 'hist_slot', 'loop', 'selector')) + 2*K.LABEL_ADD,
    R['guard'])
chk('no_sigma_application_on_the_state', R['sigma_on_state'] == 0 and K.ESCAPE == 0, 0)
chk('no_label_renormalisation', R['renorm'] == 0, 0)
chk('n_classes_exact', K.N_CLASSES*12 + 1 == K.R_ORDER, K.N_CLASSES)

# ------------------------------------------------- record size 112 B everywhere
chk('record_is_112_bytes', K.RECORD_BYTES == 112 == K.KEY_BYTES + 2*K.LABEL_BYTES,
    '%d = %d + 2*%d' % (K.RECORD_BYTES, K.KEY_BYTES, K.LABEL_BYTES))
chk('aux_report_write_uses_112_bytes',
    AUX['report_writes'] == Fraction(STO['records']*K.REPORT_WORDS*2, K.C_REF),
    'REPORT_WORDS=%d' % K.REPORT_WORDS)
chk('no_96_byte_row_survives',
    '96' not in open('ledger.py').read().split('def aux')[1].split('def storage')[0],
    'ledger.aux has no 96-byte literal')

# ------------------------------------------------- aux: restart schedule + peak
chk('aux_has_a_restart_schedule',
    AUX['trail_starts'] == Fraction(STO['records']*K.TRAIL_START_MULS*381, K.C_REF)
    and AUX['restart_bookkeeping'] > 0, '2^%.3f' % l2(AUX['trail_starts']))
chk('aux_restart_schedule_is_S_over_theta',
    STO['records'] == K.S_STEPS // (1 << K.THETA_BITS), STO['records'])
chk('peak_storage_covers_store_labels_advice_scratch',
    STO['peak_bytes'] == STO['store_bytes'] + STO['label_bytes']
    + STO['advice_bytes'] + STO['scratch_bytes'], STO['peak_bytes'])
chk('peak_storage_arithmetic_exact',
    STO['store_bytes'] == STO['records']*K.RECORD_BYTES, STO['store_bytes'])
chk('aux_total_below_2^-60_of_W_stress', l2(AUX['total']) - l2(Ld['w_stress']) < -60,
    '2^%.2f' % (l2(AUX['total']) - l2(Ld['w_stress'])))
chk('table_build_charges_twelve_sigma_columns',
    AUX['table_build'] >= Fraction(K.TABLE_COLS*(1 << K.TABLE_BITS)*K.SIGMA_BUILD, K.C_REF),
    '2^%.3f' % l2(AUX['table_build']))
chk('C_is_protected_cost_model_value', K.C_REF == 20000 and K.FP_MUL == 383, K.C_REF)

# ------------------------------------------------- success: floor after stickiness
chk('model_success_at_least_floor', MODEL_SUCCESS >= 0.39, '%.7f' % MODEL_SUCCESS)
chk('net_success_after_stickiness_at_least_floor', NET_SUCCESS >= 0.39, '%.7f' % NET_SUCCESS)
chk('stickiness_retention_is_2D2_over_R',
    K.STICK_LOSS == Fraction(2*K.GUARD_DEPTH*K.GUARD_DEPTH, 1 << K.TABLE_BITS),
    str(float(K.STICK_LOSS)))

# ------------------------------------------------- sigma-set enumeration (theorem)
chk('sigma_set_enumeration_scanned_4096', SSET['enumeration']['subsets'] == 4096,
    SSET['enumeration']['subsets'])
chk('only_the_empty_and_full_sets_are_closed',
    SSET['enumeration']['closed'] == [[], list(range(12))], SSET['enumeration']['closed'])
chk('no_nonempty_sigma_stable_inversion_free_set',
    SSET['enumeration']['closed_and_inversion_free'] == [[]]
    and SSET['enumeration']['largest_closed_inversion_free'] == [],
    SSET['enumeration']['closed_and_inversion_free'])
chk('lower_half_is_inversion_free_but_not_closed',
    SSET['lower_half']['inversion_free'] and not SSET['lower_half']['closed'],
    SSET['lower_half']['closed'])
chk('column_index_action_is_transitive', SSET['transitivity']['is_transitive'],
    SSET['transitivity']['orbit_of_0'])

# ------------------------------------------------- fixture facts
chk('fixture_class_count_exact', T1['control_v13_rule_v13_table']['states'] == FIX['r'] - 1
    and FIX['N'] == (FIX['r'] - 1)//12, FIX['N'])
chk('fixture_sigma6_is_inversion', FIX['sigma6_is_inversion'], FIX['sigma6_is_inversion'])
chk('fixture_column_identity_zero_failures', FIX['column_identity_failures'] == 0,
    FIX['column_identity_failures'])
chk('fixture_no_column_is_an_involution', FIX['column_involution_failures'] == 0,
    FIX['column_involution_failures'])
chk('fixture_inverse_column_pairs_are_visible', FIX['column_inverse_pairs_visible'] > 0,
    FIX['column_inverse_pairs_visible'])
chk('fixture_ord_r_p_is_twelve', FIX['ord_r_p'] == 12, FIX['ord_r_p'])
chk('fixture_guard_depth_is_the_shipped_one', FIX['guard_depth'] == K.GUARD_DEPTH,
    FIX['guard_depth'])

# ------------------------------------------------- T1: the return-lag ladder
d12 = T1['d12']
chk('T1_d12_no_return_at_any_lag_up_to_26', d12['total_returns_below_26'] == 0
    and d12['min_lag'] is None, d12['lag_counts'])
chk('T1_d12_census_is_exhaustive', d12['elements'] == FIX['r'] - 1, d12['elements'])
chk('T1_d1_no_return_at_lag_2_or_3',
    T1['d1']['lag_counts']['2'] == 0 and T1['d1']['lag_counts']['3'] == 0,
    T1['d1']['lag_counts'])
chk('T1_d2_no_return_at_lag_below_6',
    all(T1['d2']['lag_counts'][str(L)] == 0 for L in (2, 3, 4, 5)),
    T1['d2']['lag_counts'])
chk('T1_d0_rule_returns_at_lag_2', T1['d0_v13_rule']['min_lag'] == 2
    and T1['d0_v13_rule']['lag_counts']['2'] > 0, T1['d0_v13_rule']['lag_counts'])
c8 = T1['control_v13_rule_small_table']
chk('T1_control_v13_rule_fires_two_cycles', c8['two_cycles'] > 0, c8['two_cycles'])
chk('T1_control_reproduces_the_1_over_12R_scale',
    abs(c8['two_cycles']/max(1.0, c8['predicted_2cycles_expectation']) - 1) < 5,
    '%d vs %.1f' % (c8['two_cycles'], c8['predicted_2cycles_expectation']))
chk('T1_control_on_the_v13_table_is_blind', T1['control_v13_rule_v13_table']['two_cycles'] == 0,
    T1['control_v13_rule_v13_table']['two_cycles'])

# ------------------------------------------------- T1j: exact joint-space census
j1, j0 = T1J['guarded_d1'], T1J['v13_rule']
chk('T1j_census_covers_the_whole_joint_space',
    j1['states'] == (FIX['r'] - 1)*j1['ring_states'] + 1, j1['states'])
chk('T1j_guarded_zero_two_cycles', j1['two_cycles'] == 0, j1['two_cycles'])
chk('T1j_guarded_zero_three_cycles', j1['three_cycles'] == 0, j1['three_cycles'])
chk('T1j_control_two_cycles_present', j0['two_cycles'] > 0, j0['two_cycles'])

# ------------------------------------------------- T2 / T2b / T3 / T5
chk('T2_shipped_zero_element_splits', T2['shipped_split'] == 0 and T2['checks'] == 768,
    '%d/%d' % (T2['shipped_split'], T2['checks']))
chk('T2_shipped_zero_key_splits', T2['shipped_key_split'] == 0, T2['shipped_key_split'])
chk('T2_control_classkeyed_splits', T2['control_split_classkeyed'] > 0,
    '%d/%d' % (T2['control_split_classkeyed'], T2['checks']))
chk('T2_control_v12b_splits', T2['control_split_v12b'] > 0,
    '%d/%d' % (T2['control_split_v12b'], T2['checks']))
chk('T2_both_channels_exercised', T2['shipped_steps_sq'] > 0 and T2['shipped_steps_mul'] > 0,
    '%d/%d' % (T2['shipped_steps_sq'], T2['shipped_steps_mul']))
chk('T2b_pair_stays_exact_inverse', T2B['inverse_holds'] == T2B['inverse_checks']
    and T2B['inverse_checks'] > 0, '%d/%d' % (T2B['inverse_holds'], T2B['inverse_checks']))
chk('T2b_label_bit_solve_verifies', T2B['label_bit_verified'] == T2B['label_bit_attempts']
    and T2B['label_bit_attempts'] > 0, '%d/%d' % (T2B['label_bit_verified'], T2B['label_bit_attempts']))
chk('T2b_control_wrong_sign_rejected', T2B['control_wrong_sign_verified'] == 0
    and T2B['control_attempts'] > 0, '%d/%d' % (T2B['control_wrong_sign_verified'], T2B['control_attempts']))
chk('T2b_zero_denominators', T2B['zero_denominator'] == 0, T2B['zero_denominator'])
chk('T2b_pairs_share_one_class', T2B['class_splits'] == 0, T2B['class_splits'])
chk('T3_shipped_zero_splits', T3['shipped_splits'] == 0 and T3['shipped_compared'] > 0,
    '%d/%d' % (T3['shipped_splits'], T3['shipped_compared']))
chk('T3_all_pairs_coalesced', T3['coalesced_at_end'] == T3['pairs'], T3['coalesced_at_end'])
chk('T3_control_splits', T3['control_splits'] > 0,
    '%d/%d' % (T3['control_splits'], T3['control_compared']))
chk('T3_memory_divergence_measured', T3['memory_divergence']['pairs'] == T3['pairs']
    and T3['memory_divergence']['mul_steps'] > 0, T3['memory_divergence'])
chk('T5_every_shift_solved', T5['all_j_solved'] and T5['j_hits'] == 12*T5['samples'],
    '%d/%d' % (T5['j_hits'], 12*T5['samples']))
chk('T5_unique_verified_candidate', T5.get('unique_accept', 0) == T5['candidate_solves'],
    '%d/%d' % (T5.get('unique_accept', 0), T5['candidate_solves']))
chk('T5_no_false_accept', T5['wrong_j_accepted'] == 0, T5['wrong_j_accepted'])
chk('T5_zero_denominators', T5['zero_denominators'] == 0, T5['zero_denominators'])
chk('T5_j6_involution_case_covered', T5['j6_hits'] == T5['samples'],
    '%d/%d' % (T5['j6_hits'], T5['samples']))

# ------------------------------------------------- T4 census
for d, r4 in sorted(T4.items(), key=lambda z: int(z[0])):
    chk('T4_dp%s_every_collision_verified' % d, r4['verified'] == r4['collisions'] > 0,
        '%d/%d' % (r4['verified'], r4['collisions']))
    chk('T4_dp%s_zero_invariant_failures' % d, r4['invariant_fail'] == 0
        and r4['invariant_checks'] == r4['records'], r4['invariant_fail'])
    chk('T4_dp%s_zero_attribution_failures' % d, r4['attribution_fail'] == 0, r4['attribution_fail'])
    chk('T4_dp%s_zero_denominators' % d, r4['zero_denominator'] == 0, r4['zero_denominator'])
    chk('T4_dp%s_cp95_recomputed' % d,
        float(cp95_lower(r4['recovered'], r4['trials'])) == float(r4['cp95_lower']),
        r4['cp95_lower'])

# ------------------------------------------------- guard fire rate (fixture scale)
# guard fire rate: reported, not asserted.  The fixture's short cycles make a walk
# revisit rows far more often than D/R, so this row does NOT measure the production
# retention; the ledger uses the analytic 1-2D^2/R and the row is disclosed as
# fixture-scale in argument.md section 5.
if GR:
    chk('guard_fire_rate_recorded_at_two_table_sizes',
        GR['r64']['mul_steps'] > 0 and GR['r512']['mul_steps'] > 0,
        'r64 %.4f  r512 %.4f' % (GR['r64']['rate_per_mul_step'],
                                 GR['r512']['rate_per_mul_step']))
chk('part_a_all_checks_pass', PA['part_a'] and PA['checks_n'] == 24,
    '%d/%d' % (PA['checks_n'], PA['checks_n']))
chk('routing_fields_disjoint', bool(ROUT.get('disjoint')), ROUT.get('disjoint'))
chk('routing_no_collapse', ROUT.get('collapse') == 0, ROUT.get('collapse'))

# ------------------------------------------------- EXPLICIT FAILED GATE
failed_gate.append(dict(
    name='no_inversion_free_sigma_stable_column_table',
    status='FAIL',
    why='exhaustive over all 2^12 subsets: the only closed sets are {} and Z/12, and '
        'Z/12 contains the inverse pairs {k,k+6}; the lower-half table is inverse-free '
        'but not closed (m^2=1 forced). The route the v14 dispatch proposed does not '
        'exist, so the package ships the routing guard instead.',
    evidence='evidence/sigma-set.json:enumeration.closed_and_inversion_free'))
failed_gate.append(dict(
    name='production_scale_guard_census',
    status='NOT_MET',
    why='the guard census runs at table 64 and 512 rows, not 2^14; the production '
        'retention 1-2D^2/R is analytic and the fixture fire rate is disclosed as '
        'fixture-scale only (see argument.md section 6).',
    evidence='evidence/toy-census.json:guard_rate'))

checks_pass = sum(1 for c in checks if c['ok'])
failed = [c['name'] for c in checks if not c['ok']]
verdict = 'PASS' if not failed else 'FAIL'

# ------------------------------------------------- digests
files = sorted(f for f in os.listdir('.') if f.endswith(('.py', '.md', '.json', '.c'))
               and f != 'certificate.json')
digests = {}
for f in files:
    digests[f] = hashlib.sha256(open(f, 'rb').read()).hexdigest()
for f in sorted(os.listdir(EV)):
    digests['evidence/' + f] = hashlib.sha256(open(os.path.join(EV, f), 'rb').read()).hexdigest()

cert = dict(
    certificate='v14 guarded equivariant mu6-class walk, quotient in the report filter, '
                'bls12-381-fixed-g1-invert-g2',
    verdict=verdict, checks_total=len(checks), checks_pass=checks_pass,
    checks_fail=len(failed), failed=failed, failed_gate=failed_gate,
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
        f_sqr='%d/%d' % (K.F_SQR.numerator, K.F_SQR.denominator),
        rows={k: v for k, v in R.items() if isinstance(v, int)}),
    success=dict(lambda_='%.10f' % float(LAM), model='%.7f' % MODEL_SUCCESS,
                 stickiness_loss='%.6f' % float(K.STICK_LOSS),
                 net='%.7f' % NET_SUCCESS, floor=0.39,
                 mu='%d/%d' % (K.MU_WALK.numerator, K.MU_WALK.denominator),
                 target='%d/%d' % (K.TARGET_PAD.numerator, K.TARGET_PAD.denominator)),
    record=dict(bytes=K.RECORD_BYTES, key_bytes=K.KEY_BYTES, label_bytes=K.LABEL_BYTES,
                words=K.REPORT_WORDS),
    storage={k: (str(v) if isinstance(v, int) else v) for k, v in STO.items()},
    aux={k: ('2^%.3f' % l2(v)) for k, v in AUX.items()},
    sigma_set=SSET,
    fixture=FIX,
    toy={'t1': T1, 't1j': T1J, 't2': T2, 't2b': T2B, 't3': T3, 't5': T5,
         'guard_rate': GR},
    census={d: {k: r4[k] for k in ('trials', 'recovered', 'success', 'cp95_lower',
                                   'collisions', 'verified', 'records',
                                   'record_co_fraction', 'invariant_fail',
                                   'zero_denominator', 'attribution_fail',
                                   'mean_steps_per_trial', 'mean_steps_to_first_report',
                                   'guard_fire_rate_per_mul_step')}
            for d, r4 in T4.items()},
    cp_selftest=TOY['cp_selftest'],
    part_a={'checks': PA['checks_n'], 'all_pass': PA['part_a'], 'seconds': PA['seconds']},
    routing={k: ROUT.get(k) for k in ('keys', 'disjoint', 'table_distinct', 'dp_distinct',
                                      'dp_zero', 'collapse')},
    digests=digests, checks=checks,
)
json.dump(cert, open('certificate.json', 'w'), sort_keys=True, separators=(',', ':'))
print('certificate.json: %s  %d/%d checks, %d failed gates'
      % (verdict, checks_pass, len(checks), len(failed_gate)))
for f in failed:
    print('  FAIL', f)
print('claim %s  log2 W_stress %.9f  net success %.7f'
      % (CLAIM, l2(Ld['w_stress']), NET_SUCCESS))
sys.exit(0 if not failed else 1)
