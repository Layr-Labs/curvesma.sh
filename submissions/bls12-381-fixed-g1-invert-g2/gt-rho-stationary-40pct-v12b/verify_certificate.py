#!/usr/bin/env python3
"""Regenerate certificate.json from the raw evidence JSONs and check every number.

Nothing here retypes a number from the note: the ledger is recomputed from
constants.py, the toy rows are re-summed from the raw a-*.json / d-*.json /
e-equivariance.json, the modular-multiply row is re-read from
evidence/modmul-count.json, the Clopper-Pearson floor is recomputed from the
pooled sweep counts, and the claim is tested in exact integer arithmetic.

Exit status 1 on any failed check.  Run with PYTHONDONTWRITEBYTECODE=1 (intake
refuses .pyc files).
"""
import json
import os
import sys
from decimal import Decimal
from fractions import Fraction
from math import log2

D = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, D)
import constants as K
import ledger as L

CHECKS = []


def chk(name, ok, detail=''):
    CHECKS.append(dict(name=name, ok=bool(ok), detail=str(detail)))
    return bool(ok)


def ev(name):
    p = os.path.join(D, 'evidence', name)
    if not os.path.exists(p):
        p = os.path.join(D, name)
    return json.load(open(p))


def cp95_lower(k, n, alpha=Decimal('0.05')):
    """One-sided Clopper-Pearson lower bound: smallest p with P[X>=k|p,n] >= alpha."""
    k = int(k); n = int(n)
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p):
        p = Decimal(p); q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal('1e-40') * s:
                break
        return s

    lo, hi = Decimal('0.30'), Decimal('0.50')
    for _ in range(80):
        mid = (lo + hi) / 2
        if tail(mid) > alpha:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal('1e-14'):
            break
    return lo


def main():
    L1 = L.ledgers()
    r = L1['rows']
    # ---------- 1. ledger rows: exact cancellation, no hidden row
    chk('fold_is_sum_of_nine_rows', r['fold'] == sum(r[k] for k in (
        'selector_ab', 'selector_zeta', 'fold_adds', 'fold_cmp', 'hashing',
        'dp_test', 'hist_slot', 'loop', 'witness_add')), r['fold'])
    chk('fold_is_1813', r['fold'] == 1813, r['fold'])
    chk('cyclo_sqr_is_11632', r['cyclo_sqr'] == 11632, r['cyclo_sqr'])
    chk('renorm_row_present_and_counted', r['renorm'] == K.RENORM == 2 * K.MODMUL_INSTR,
        r['renorm'])
    chk('no_window_or_detector_row', 'window' not in r and 'detector' not in r, 'absent')
    chk('mul_step_identity', r['mul_step'] == K.C_REF + r['fold'] + r['renorm'],
        r['mul_step'])
    chk('sqr_step_identity', r['sqr_step'] == r['cyclo_sqr'] + r['fold'] + r['renorm'],
        r['sqr_step'])
    # ---------- 2. the counted renormalisation routine
    mm = ev('modmul-count.json')
    chk('modmul_counted_at_143', mm['instructions_O2'] == K.MODMUL_INSTR,
        mm['instructions_O2'])
    chk('modmul_differential_zero_mismatch', mm['differential_mismatches'] == 0,
        '%d cases' % mm['differential_cases'])
    chk('renorm_equals_two_modmul', K.RENORM == 2 * mm['instructions_O2'], K.RENORM)
    # ---------- 3. S, prices, W, claim -- exact integer arithmetic only
    chk('c_raw_is_exact_94805_over_100000', L1['c_raw'] == Fraction(94805, 100000),
        str(L1['c_raw']))
    chk('w_raw_is_S_times_c_raw', L1['w_raw'] == L1['S'] * L1['c_raw'], 'exact')
    pad = float(K.CLAIM) - log2(float(L1['w_stress']))
    num, den = K.CLAIM.numerator, K.CLAIM.denominator
    chk('claim_ge_W_stress_exact_integer', Fraction(2) ** num >= L1['w_stress'] ** den,
        '%.6f bits of pad' % pad)
    chk('claim_reserve_at_least_0p005', pad >= 0.005, '%.6f' % pad)
    chk('claim_clears_gate_minus_0p0065', float(K.CLAIM) <= 125.84 - 0.0065, 125.84)
    # ---------- 4. auxiliary rows
    A = L.aux()
    worst = min(log2(float(v)) - log2(float(L1['w_stress']))
                for k, v in A.items() if k != 'total')
    chk('aux_all_below_2^-60_of_W_stress', worst < -60, '%.2f' % worst)
    # ---------- 5. toy audits
    for b in (5, 7, 9, 11):
        a = ev('a-%d.json' % b)
        chk('audit_A_zero_fail_%d' % b, a['audit_A_fail'] == 0 and a['audit_A_pass'] > 0,
            a['audit_A_pass'])
        chk('audit_B_equivariance_zero_fail_%d' % b,
            a['audit_B_fail'] == 0 and a['audit_B_pass'] > 0, a['audit_B_pass'])
        chk('audit_E_inverse_free_zero_fail_%d' % b,
            a['audit_E_fail'] == 0 and a['audit_E_pass'] > 0, a['audit_E_pass'])
        chk('table_inverse_free_%d' % b, a['table_inverse_free'] is True, 'true')
        chk('table_all_lower_half_%d' % b, a['table_all_lower_half'] is True, 'true')
        chk('frob_pow_differential_zero_%d' % b,
            a['frob_pow_differential']['fail'] == 0, a['frob_pow_differential'])
        chk('control_old_rule_splits_nonzero_%d' % b,
            a['control_C_old_rule_splits'] > 0, a['control_C_old_rule_splits'])
        chk('kappa_no_degenerate_tie_%d' % b, a['kappa_degenerate'] == 0, 0)
    # ---------- 6. sweeps: verified recovery, invariant, orbit census
    sweeps = {}
    tot_s = tot_n = 0
    present = [b for b in (5, 7, 9, 11) if os.path.exists(os.path.join(D, 'd-%d.json' % b))]
    chk('binding_densities_present', 7 in present and 9 in present, present)
    for b in present:
        d = ev('d-%d.json' % b)
        sweeps[b] = d
        chk('sweep_trials_equal_request_%d' % b, d['trials'] == d['trials_requested'],
            d['trials'])
        chk('sweep_invariant_zero_fail_%d' % b, d['inv_failures'] == 0, d['inv_checks'])
        chk('sweep_no_false_confirm_%d' % b, d['false_confirm'] == 0, 0)
        chk('sweep_every_merge_verified_%d' % b, d['solve_success'] == d['merges'],
            '%d/%d' % (d['solve_success'], d['merges']))
        chk('sweep_table_inverse_free_%d' % b, d['table_inverse_free'] is True, 'true')
        chk('sweep_no_short_orbit_%d' % b, d['orbit_censor'] == 0, d['orbit_censor'])
        tot_s += d['success']; tot_n += d['trials']
    chk('pooled_trials_at_least_6400_at_2^-7_and_2^-9',
        sweeps[7]['trials'] + sweeps[9]['trials'] >= 6400,
        sweeps[7]['trials'] + sweeps[9]['trials'])
    cp = cp95_lower(tot_s, tot_n)
    chk('pooled_cp95_lower_at_least_0p39', cp >= Decimal('0.39'), str(cp))
    # ---------- 7. equivariance: 0 splits, and a control that does split
    e = ev('e-equivariance.json')
    chk('equivariance_zero_splits', e['totals']['splits'] == 0, e['totals']['splits'])
    chk('equivariance_d1_zero_inverse_closure',
        e['totals']['d1_hits'] == 0, e['totals']['d1_hits'])
    chk('equivariance_control_reproduces_the_defect',
        e['totals']['control_d1_hits'] == sum(f['control_v10f_d1_checks'] for f in e['per_fixture'])
        and sum(f['control_v10f_d1_checks'] for f in e['per_fixture']) > 0,
        '%d/%d ladder multipliers have their inverse in the ladder'
        % (e['totals']['control_d1_hits'],
           sum(f['control_v10f_d1_checks'] for f in e['per_fixture'])))
    chk('p6_is_minus_one_on_every_fixture', e['all_fixtures_p6_is_minus_one'] is True,
        'true')
    # ---------- 8. routing and part A
    ro = ev('routing-out.json')
    chk('routing_fields_bit_disjoint', ro['disjoint'] is True, 'true')
    chk('routing_no_field_collapse', ro['collapse'] == 0, ro['collapse'])
    chk('routing_no_checkpoint_field', ro['has_ckpt_field'] is False, 'false')
    chk('routing_window_is_zero', ro['window'] == 0, ro['window'])
    pa = ev('part_a-out.json')
    chk('part_a_all_pass', pa.get('fail', 0) == 0 and pa.get('failures', 0) == 0,
        json.dumps({k: v for k, v in pa.items() if 'fail' in k.lower()}))

    failed = [c for c in CHECKS if not c['ok']]
    res = dict(
        certificate='v12b canonicalising class walk, bls12-381-fixed-g1-invert-g2',
        checks_total=len(CHECKS), checks_pass=len(CHECKS) - len(failed),
        checks_fail=len(failed), verdict='PASS' if not failed else 'FAIL',
        failed=[c['name'] for c in failed],
        ledger=dict(rows={k: v for k, v in r.items() if isinstance(v, int)},
                    S=L1['S'], log2_S=log2(L1['S']),
                    c_raw=str(L1['c_raw']), c_stress=str(L1['c_stress']),
                    log2_W_raw=log2(float(L1['w_raw'])),
                    log2_W_stress=log2(float(L1['w_stress'])),
                    claim=str(K.CLAIM), claim_minus_W_stress=pad),
        modmul=dict(instructions=mm['instructions_O2'],
                    cases=mm['differential_cases'],
                    mismatches=mm['differential_mismatches'], renorm=K.RENORM),
        toy=dict(pooled_recovery='%d/%d' % (tot_s, tot_n),
                 pooled_cp95_lower=str(cp),
                 per_density={str(b): dict(recovery=d['recovery'],
                                           invariants=d['inv_checks'],
                                           inv_failures=d['inv_failures'],
                                           exact_state_repeats=d['exact_state_repeats'],
                                           orbit_censor=d['orbit_censor'])
                              for b, d in sweeps.items()}),
        equivariance=dict(checks=e['totals']['checks'], splits=e['totals']['splits'],
                          d1_checks=e['totals']['d1_checks'],
                          d1_hits=e['totals']['d1_hits'],
                          control_splits=e['totals']['control_splits'],
                          control_d1_hits=e['totals']['control_d1_hits']),
        checks=CHECKS,
    )
    json.dump(res, open(os.path.join(D, 'certificate.json'), 'w'), indent=1)
    print('CERTIFICATE', res['verdict'], '%d/%d' % (res['checks_pass'],
                                                    res['checks_total']))
    for c in CHECKS:
        print('  %-52s %-5s %s' % (c['name'], 'PASS' if c['ok'] else 'FAIL', c['detail']))
    return 0 if not failed else 1


if __name__ == '__main__':
    sys.exit(main())
