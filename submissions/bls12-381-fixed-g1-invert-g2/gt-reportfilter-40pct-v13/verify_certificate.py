#!/usr/bin/env python3
"""Independent re-derivation of the v13 certificate from the RAW evidence.

It re-reads constants.py, ledger.py, part_a-out.json, routing-out.json and the
toy census JSONs, recomputes every number the certificate prints, re-hashes every
file against the manifest, and fails on any disagreement.  It never reads a
submitted summary as truth and it never takes the certificate's word for a count.
"""
import hashlib
import json
import os
import sys
from decimal import Decimal

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import constants as K
import ledger as L
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
eq('claim_above_floor_check', claim >= logw + Decimal('2e-5'), True)
eq('model_floor', 1 - 1/float(Ld['model'][1]) >= 0.39, True)
eq('claim_numerator', int(claim * 10**4), C['claim']['numerator'])
eq('claim_denominator', 10**4, C['claim']['denominator'])

TOY = json.load(open('evidence/toy-census.json'))
eq('t1_states', TOY['t1']['states'], C['toy']['t1']['states'])
eq('t1_shipped_state_2cycles', TOY['t1']['shipped_state_2cycles'], 0)
eq('t1_shipped_class_2cycles', TOY['t1']['shipped_class_2cycles'], 0)
eq('t1_control_v12b_2cycles', TOY['t1']['control_v12b_state_2cycles'], 28)
from fractions import Fraction as _F
eq('measured_f_sqr',
   _F(TOY['t1']['step_kinds']['sq'], TOY['t1']['states']),
   _F(*[int(v) for v in C['ledger']['f_sqr_measured'].split('/')]))
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

for d in ('5', '7', '9'):
    f = 'evidence/t4-dp%s.json' % d
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
