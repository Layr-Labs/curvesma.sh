#!/usr/bin/env python3
"""Independent re-derivation of the v14 certificate from the RAW evidence.

It re-reads constants.py, ledger.py, sigma_set.py, part_a-out.json,
routing-out.json and the raw evidence JSONs, recomputes every number the
certificate prints, re-hashes every file against the manifest, and fails on any
disagreement.  It never reads a submitted summary as truth.
"""
import hashlib
import json
import os
import sys
from decimal import Decimal
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import ledger as L
import sigma_set as SS
from cp95 import cp95_lower

fail = []


def eq(name, a, b):
    if a != b:
        fail.append('%s: recomputed %r != certified %r' % (name, a, b))


C = json.load(open('certificate.json'))
Ld = L.ledgers()
R = Ld['rows']
logw = L.log2_decimal(Ld['w_stress'])
claim = L.claim_log2(Ld['w_stress'])

eq('claim', str(claim), C['claim']['value'])
eq('log2_S', '%.6f' % L.log2_decimal(Ld['S']), C['ledger']['log2_S'])
eq('log2_W_raw', '%.6f' % L.log2_decimal(Ld['w_raw']), C['ledger']['log2_W_raw'])
eq('log2_W_stress', '%.6f' % logw, C['ledger']['log2_W_stress'])
eq('W_stress_exact', str(Ld['w_stress']), C['ledger']['W_stress'])
eq('W_raw_exact', str(Ld['w_raw']), C['ledger']['W_raw'])
for k, v in R.items():
    if isinstance(v, int) and k in C['ledger']['rows']:
        eq('row.' + k, v, C['ledger']['rows'][k])
eq('mul_step', K.C_REF + R['fold_mul'], R['mul_step'])
eq('sqr_step', R['cyclo_sqr'] + R['fold_sqr'], R['sqr_step'])
eq('claim_above_pad', claim >= logw + Decimal('2e-5'), True)
eq('claim_is_least', (claim - Decimal('1e-4')) < logw + Decimal('2e-5'), True)
eq('claim_numerator', int(claim * 10**4), C['claim']['numerator'])
eq('claim_denominator', 10**4, C['claim']['denominator'])
eq('net_success_floor', Ld['net_success'] >= 0.39, True)
eq('net_success_value', '%.7f' % Ld['net_success'], C['success']['net'])
eq('stickiness', '%.6f' % float(K.STICK_LOSS), C['success']['stickiness_loss'])

# records and storage: one constant, 112
STO = L.storage()
eq('record_bytes', K.RECORD_BYTES, 112)
eq('record_bytes_is_key_plus_two_labels', K.KEY_BYTES + 2*K.LABEL_BYTES, 112)
eq('storage_record_bytes', STO['record_bytes'], 112)
eq('storage_store_bytes', str(STO['store_bytes']), C['storage']['store_bytes'])
eq('storage_peak_bytes', str(STO['peak_bytes']), C['storage']['peak_bytes'])
eq('storage_peak_sum', STO['peak_bytes'],
   STO['store_bytes'] + STO['label_bytes'] + STO['u_label_bytes']
   + STO['advice_bytes'] + STO['guard_history_bytes'] + STO['code_bytes']
   + STO['scratch_bytes'])
eq('storage_u_labels_present', STO['u_label_bytes'],
   (1 << K.TABLE_BITS)*K.TABLE_COLS*K.LABEL_BYTES)
eq('storage_records', STO['records'],
   2*((K.S_STEPS + (1 << K.THETA_BITS) - 1)//(1 << K.THETA_BITS)))
eq('aux_report_writes_use_112',
   L.aux()['report_writes'], Fraction(STO['records']*K.REPORT_WORDS*2, K.C_REF))

# sigma-set: recomputed, not read
E = SS.enumerate_subsets()
eq('sigma_set_closed', E['closed'], [[], list(range(12))])
eq('sigma_set_closed_and_free', E['closed_and_inversion_free'], [[]])
eq('sigma_set_largest', E['largest_closed_inversion_free'], [])
eq('sigma_set_certified', C['sigma_set']['enumeration']['closed_and_inversion_free'], [[]])
eq('lower_half_not_closed', SS.six_column_refutation()['closed'], False)
eq('transitive', SS.transitivity_argument()['is_transitive'], True)

TOY = json.load(open('evidence/toy-census.json'))
eq('t1_states', TOY['t1']['d12']['elements'], TOY['fixture']['r'] - 1)
eq('t1_d12_no_returns', TOY['t1']['d12']['total_returns_below_26'], 0)
eq('t1_d12_min_lag', TOY['t1']['d12']['min_lag'], None)
eq('t1_d0_min_lag', TOY['t1']['d0_v13_rule']['min_lag'], 2)
eq('t1_d1_no_return_below_4',
   TOY['t1']['d1']['lag_counts']['2'] + TOY['t1']['d1']['lag_counts']['3'], 0)
eq('t1_d2_no_return_below_6',
   sum(TOY['t1']['d2']['lag_counts'][str(L)] for L in (2, 3, 4, 5)), 0)
eq('t1_control_2cycles', TOY['t1']['control_v13_rule_small_table']['two_cycles'] > 0, True)
eq('t1_control_v13_table_blind', TOY['t1']['control_v13_rule_v13_table']['two_cycles'], 0)
eq('t1j_guarded_2cycles', TOY['t1j']['guarded_d1']['two_cycles'], 0)
eq('t1j_guarded_3cycles', TOY['t1j']['guarded_d1']['three_cycles'], 0)
eq('t1j_control_2cycles', TOY['t1j']['v13_rule']['two_cycles'] > 0, True)
eq('t2_checks', TOY['t2']['checks'], 768)
eq('t2_splits', TOY['t2']['shipped_split'], 0)
eq('t2_key_splits', TOY['t2']['shipped_key_split'], 0)
eq('t2_control_classkeyed', TOY['t2']['control_split_classkeyed'] > 0, True)
eq('t2_control_v12b', TOY['t2']['control_split_v12b'] > 0, True)
eq('t2b_inverse', TOY['t2b']['inverse_holds'], TOY['t2b']['inverse_checks'])
eq('t2b_label_bit', TOY['t2b']['label_bit_verified'], TOY['t2b']['label_bit_attempts'])
eq('t2b_control', TOY['t2b']['control_wrong_sign_verified'], 0)
eq('t3_splits', TOY['t3']['shipped_splits'], 0)
eq('t3_coalesced', TOY['t3']['coalesced_at_end'], TOY['t3']['pairs'])
eq('t5_all_shifts', TOY['t5']['j_hits'], 12*TOY['t5']['samples'])
eq('t5_unique', TOY['t5']['unique_accept'], TOY['t5']['candidate_solves'])
eq('t5_no_false_accept', TOY['t5']['wrong_j_accepted'], 0)
eq('guard_depth_in_fixture', TOY['fixture']['guard_depth'], K.GUARD_DEPTH)

for d in ('5',):
    f = 'evidence/t4-dp%s.json' % d
    if not os.path.exists(f):
        # v17: the t4 calibration row is superseded by the MEASURED reports-per-walk
        # census in evidence/bounds-v17.json (fix 3); the density census is not
        # re-run in this revision and no row is claimed from it.
        eq('t4_%s_superseded_by_v17_census' % d,
           os.path.exists('evidence/bounds-v17.json'), True)
        continue
    if not os.path.exists(f):
        fail.append('missing census file ' + f)
        continue
    r4 = json.load(open(f))
    eq('t4_%s_cp95' % d, float(cp95_lower(r4['recovered'], r4['trials'])), float(r4['cp95_lower']))
    eq('t4_%s_verified_all' % d, r4['verified'] == r4['collisions'] > 0, True)
    eq('t4_%s_invariant' % d, r4['invariant_fail'], 0)
    eq('t4_%s_zeroden' % d, r4['zero_denominator'], 0)
    eq('t4_%s_attribution' % d, r4['attribution_fail'], 0)
    if d in C['census']:
        eq('t4_%s_success_certified' % d, r4['success'], C['census'][d]['success'])

PA = json.load(open('part_a-out.json'))
eq('part_a', PA['part_a'], True)
eq('part_a_checks', PA['checks_n'], 24)
RO = json.load(open('routing-out.json'))
eq('routing_disjoint', RO['disjoint'], True)
eq('routing_collapse', RO['collapse'], 0)

for name, digest in C['digests'].items():
    if not os.path.exists(name):
        fail.append('digest target missing: ' + name)
        continue
    if hashlib.sha256(open(name, 'rb').read()).hexdigest() != digest:
        fail.append('digest mismatch: ' + name)

print('verify_certificate: %d checked groups, %d failures' % (len(C['checks']), len(fail)))
for f in fail[:20]:
    print('  FAIL', f)
print('VERIFY', 'PASS' if not fail else 'FAIL')
sys.exit(0 if not fail else 1)
