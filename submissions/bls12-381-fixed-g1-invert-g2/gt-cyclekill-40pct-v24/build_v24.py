#!/usr/bin/env python3
"""Build every v24 artifact from the ONE ledger sum and the ONE census file."""
import json, os, sys
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ledger as L, constants as K

C, T = L.claim_from_total()
LT = float(T['log2_total'])
l2 = L.log2_decimal
A = L.aux(); r = L.rows(); l = L.ledgers()
S = l['S']
cen = json.load(open(os.path.join(HERE, 'evidence', 'census-v24.json')))
rec = cen['recovery_v24']
fam24, fam23 = cen['family_v24_step'], cen['family_v23_step']
claim_int = str(int((C * 10000).to_integral_value()))
claim_str = str(C)
aux_rows = {k: str(v) for k, v in T['aux_rows'].items()}
peak = L.storage()['peak_bytes']

cert = dict(
    version=24, track='bls12-381-fixed-g1-invert-g2',
    claim_log2=float(C), claim_integer=claim_int,
    per_step_rows=r, per_step_resum_ok=(
        r['fold_mul'] == sum(r[k] for k in ('key_hash', 'dp_test', 'hist_slot', 'loop',
                                            'table_fetch', 'label_add', 'selector',
                                            'column_addr', 'rowdiff', 'guard'))),
    work=dict(S=str(S), c_stress=str(l['c_stress']), w_stress=str(T['w_stress']),
              w_stress_log2=float(l2(T['w_stress'])),
              trap_charge=str(T['trap_charge']), trap_charge_log2=float(l2(T['trap_charge'])),
              aux_sum=str(T['aux_sum']), aux_sum_log2=float(l2(T['aux_sum'])),
              aux_reserve=str(T['aux_reserve']), aux_reserve_log2=float(l2(T['aux_reserve'])),
              total=str(T['total']), total_log2=LT, claim=str(C),
              claim_minus_total_bits=float(C) - LT,
              covered_by_claim=bool(T['total'] <= Fraction(int(claim_int), 10000))
              or bool(LT <= float(C)),
              one_sum='W_total = S*c_stress + trap_charge(rate*cap*price) + sum(aux rows) '
                      '+ aux_reserve; claim = least 4-dp log2 >= log2(W_total) + 2e-5'),
    aux=dict(rows=aux_rows, total=str(T['aux_sum']),
             covers='the claim is taken from the SUM, not from W_stress: v23 took it from '
                    'log2(W_stress) alone and the reviewer added the parts by hand'),
    structural_kill=dict(
        family='consecutive multiply steps on ONE row with kappa(y\')=kappa(y)+6 use '
               'columns k and k+6, and TBL[t][k+6] = TBL[t][k]^-1, so y returns EXACTLY',
        fix='classkey.step_v24: when the drawn row repeats the previous multiply row the '
            'row is remapped to t XOR (1 + h mod (R-1)), a bijection onto the other R-1 '
            'rows, so two consecutive multiply steps NEVER share a row',
        marginal='row marginal stays exactly uniform; pair law is uniform on ordered pairs '
                 'of distinct rows; statistical distance from iid uniform <= 1/R per step',
        control_v23_step=fam23, live_v24_step=fam24,
        dead_by_construction=bool(fam24['consecutive_same_row_multipliers'] == 0),
        cost='ROWDIFF = 16 instructions on the multiply channel, charged on every '
             'multiply step (constants.ROWDIFF): +0.00097 bits on W_stress, disclosed'),
    success=dict(
        model_lambda=str(l['model'][0]),
        model_success=1 - 1 / float(l['model'][1]),
        net=float(l['net_success']),
        trap_bound_production='P(class repeat within t steps) <= t(t-1)/2N, a union bound '
                              'over pairs with no fixture input: t = cap = 128*2^80, '
                              'N = 2^371 gives <= 2^-198; over the report interval '
                              'theta^2/2N = 2^-92.3, the row charged',
        jensen='E[1-exp(-tau^2/2N)] with tau geometric(theta) is NOT 1-exp(-E[tau]^2/2N); '
               'the shipped row uses tau >= 1 (the smallest value) so the bound is taken '
               'at the p10-quantile direction, i.e. E[1-exp(-tau^2/2N)] >= 1-exp(-1/2N) '
               'and the birthday count uses S^2(1-f^2) directly rather than a mean tau',
        unconditional_census=rec),
    census_integrity=dict(
        per_worker_seeds='random.Random(BASE + 1000*worker); starts from toy.rand_gt(that rng)',
        worker_first_starts_distinct=rec['worker_first_starts_distinct'],
        trials=rec['trials'],
        experiment_header_vs_code='cap 256 vs scheduled 264 in v23; v24 reports the ACTUAL '
                                  'trial count of every census (this file) and never a header'),
    storage=L.storage(), inequality_audit={}, note_literal_audit={}, planted_flip={})
json.dump(cert, open(os.path.join(HERE, 'certificate.json'), 'w'), indent=1, default=str)

print('claim', C, 'total_log2', T['log2_total'], 'aux_log2', float(l2(T['aux_sum'])))
print('recovery', rec['UNCONDITIONAL_RECOVERY'], rec['UNCONDITIONAL_CP95_LOWER'])
print('fam24', fam24['consecutive_same_row_multipliers'], 'fam23', fam23['consecutive_same_row_multipliers'])
print('distinct starts', rec['worker_first_starts_distinct'])
