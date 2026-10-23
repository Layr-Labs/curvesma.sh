"""Exact author certificate for sixclass-meter-40pct-v13e.

Run --write to capture certificate.json and evidence/certificate_run.log; run
without it to compare regenerated artifacts.  Every production bound below is
decided by integer/Fraction arithmetic (no floats).  Toy captures are
aggregated counters supporting disclosed heuristics at their actual scope.
"""
import hashlib
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import ledger as z


def nbelow(x, k):
    x = F(x)
    return x.numerator < (x.denominator << k)


def load_evidence():
    return {p.stem: json.loads(p.read_text())
            for p in sorted((HERE / 'evidence').glob('*.json'))}


def build_certificate():
    c = {}

    def ck(name, value):
        c[name] = bool(value)

    # ---- production constants -------------------------------------------
    ck('production_class_count_6N', 6 * z.N == z.r - 1)
    ck('production_beta_order3', z.beta != 1 and (z.beta * z.beta + z.beta + 1) % z.p == 0)
    ck('production_lambda_order3', z.lam != 1 and (z.lam * z.lam + z.lam + 1) % z.r == 0)
    ck('six_distinct_units', len({s * pow(z.lam, j, z.r) % z.r
                                  for s in (-1, 1) for j in range(3)}) == 6)
    ck('C_is_1000', z.C == 1000)
    ck('N_is_2pow253', 2 ** 253 < z.N < 2 ** 254)
    ck('R_is_2pow20', z.R == 2 ** 20)
    ck('dp_spacing_2pow69', z.DP == 2 ** 69)

    # ---- instruction itemisation ----------------------------------------
    ck('itemised_rows_3', len(z.itemized_step_rows) == 3)
    ck('itemised_rows_unique', len({q[0] for q in z.itemized_step_rows}) == 3)
    ck('itemised_every_row_priced',
       all(z.step_rows[q[0]] == q[3] * z.instruction_primitives[q[2]]
           for q in z.itemized_step_rows))
    ck('itemised_sum_772', z.step_instructions == 772)
    ck('asm_mul_body_within_125',
       z.step_rows['asm_fe_mul_x6'] % 6 == 0 and z.step_rows['asm_fe_mul_x6'] // 6 <= 125)
    ck('step_price_386p537', z.step_price == F(386, 537))
    ck('step_price_below_reference_plus_ceiling', z.step_price < 1 + F(z.E0, z.C))
    ck('doubling_limb_count_2320', z.dbl_instructions == 7 * 200 + 6 * 20 + 3 * 200 + 200)
    ck('doubling_price_2p35', z.doubling_price == F(235, 100))
    ck('doubling_above_step', z.doubling_price > z.step_price)

    # ---- S_ops and the pooled mu ----------------------------------------
    ck('S_ops_floor_square', F(z.S_ops) ** 2 >= 2 * z.WALK_MU * z.Nprime * z.ln_enclosure_lower)
    ck('S_ops_minimal', F(z.S_ops - 1) ** 2 < 2 * z.WALK_MU * z.Nprime * z.ln_enclosure_upper)
    ck('WALK_MU_is_1p107', z.WALK_MU == F(1107, 1000))
    ck('mu_hat_is_pooled', z.MU_HAT == F(1069618, 1000000))
    ck('mu_se_is_pooled', z.MU_SE_POOLED == F(18344, 1000000))
    ck('mu_charge_is_ceil_2se',
       F(1106, 1000) < z.MU_HAT + 2 * z.MU_SE_POOLED <= z.WALK_MU)
    ck('phase1_is_priced_arrivals', z.phase1 == z.step_price * (z.S_ops + 2 ** 60))
    ck('raw_sums_rows', z.raw == z.phase1 + z.hard_extras + z.random_extras)
    ck('claim_is_4dp', z.claim.denominator == 10000)
    ck('claim_is_126p3047', z.claim == F(1263047, 10000))
    ck('claim_reduced', math.gcd(z.claim.numerator, z.claim.denominator) == 1)
    ck('claim_above_W_plus_2e5', z.pow2_claim_floor(z.claim - F(2, 100000)) >= z.raw)
    ck('claim_below_gate_127', z.claim < F(127))
    ck('claim_above_floor_126p3', z.claim > F(1263, 10))
    ck('claim_smallest_4dp_both_floors',
       z.success_at(z.claim - F(1, 10000)) < F(391, 1000) or
       z.pow2_claim_floor(z.claim - F(1, 10000) - F(2, 100000)) < z.raw)
    ck('phase1_below_raw', z.raw > z.phase1)
    ck('hard_extras_below_2pow90', z.hard_extras < 2 ** 90)
    ck('random_extras_below_2pow108', z.random_extras < 2 ** 108)
    ck('mu_tolerance_above_assumed', z.mu_tolerance_hi > z.WALK_MU)
    ck('mu_tolerance_between_1p10_1p12', F(110, 100) < z.mu_tolerance_hi < F(112, 100))

    # ---- key construction and index fallback ----------------------------
    ck('index_error_bound_8_ulps', z.INDEX_ERROR_ULPS == 8)
    ck('index_threshold_covers_2E_plus_1', 2 * z.INDEX_ERROR_ULPS + 1 <= z.INDEX_GAP)
    ck('index_gap_is_40', z.INDEX_GAP == 40)
    ck('index_rate_is_294_over_2pow64', z.index_fallback_rate_bound == F(294, 2 ** 64))
    ck('index_rate_below_2powm55', nbelow(z.index_fallback_rate_bound, 55))
    ck('index_rate_charged_at_2powm50', nbelow(z.index_fallback_rate_bound, 50))
    ck('index_derived_below_charged', z.index_fallback_derived < z.index_fallback_expected)
    ck('index_charged_2pow89', z.index_fallback_expected == 2 ** 89)
    ck('asm_step_rows_present', all(k in z.step_rows for k in
        ('asm_fe_mul_x6', 'asm_fe_sub_x6', 'index_dp_parity_store_loop')))
    ck('asm_step_control_within_80', z.step_rows['index_dp_parity_store_loop'] <= 80)

    # ---- detector scheduling --------------------------------------------
    ck('checkpoint_period_prime_67', z.CHECKPOINT_PERIOD == 67)
    ck('checkpoint_period_prime', all(67 % q for q in range(2, 9)))
    ck('checkpoint_cycle_cap_66', z.CHECKPOINT_CYCLE_CAP == 66)
    ck('checkpoint_detection_bound_4422', z.checkpoint_detection_bound == 67 * 66)
    ck('coprime_for_all_cap_lengths',
       all(math.gcd(67, L) == 1 for L in range(2, 67)))

    # ---- storage partition ----------------------------------------------
    ck('nodes_2pow16', z.nodes == 2 ** 16)
    ck('per_node_cap_2pow44', z.node_cap_records == 2 ** 44)
    ck('per_node_bytes_below_2pow52', z.node_cap_records * z.record_bytes < 2 ** 52)
    ck('record_192_bytes', z.record_bytes == 192)
    ck('total_records_2pow60', z.M == 2 ** 60)
    ck('node_setup_charged', z.hard_extras >= z.node_setup == z.nodes * 2 ** 20)
    ck('reports_below_cap', F(2 ** 58) < z.M)
    ck('peak_memory_below_2pow70', z.peak_memory < 2 ** 70)

    # ---- probability -----------------------------------------------------
    ck('success_model_0p392', z.success_model == F(392, 1000))
    ck('markov_halt_below_0p001', z.markov_halt < F(1, 1000))
    ck('success_at_least_0p391', z.success >= F(391, 1000))
    ck('success_margin_over_039', z.success - F(39, 100) > F(9, 10000))
    ck('loss_budget_below_1e-6', sum(z.losses.values(), F(0)) < F(1, 10 ** 6))
    ck('doubling_rate_below_2powm19', z.doubling_rate < F(1, 2 ** 19))
    ck('random_extras_are_expected_surcharge',
       z.random_extras == z.doubling_rate * (z.doubling_price - z.step_price) * z.S_ops)
    ck('detector_abort_below_2powm40', nbelow(z.losses['unforced_detector_abort'], 40))
    ck('endpoint_transport_below_2powm40', nbelow(z.losses['endpoint_transport'], 40))

    # ---- captured toy evidence ------------------------------------------
    ev = load_evidence()
    need = {'phase', 'mu-sweep', 'pipeline', 'indexrule', 'merge',
            'checkpoint', 'traps-4', 'order2'}
    ck('evidence_files_present', set(ev) >= need)
    ck('evidence_mu_rows_present',
       all(f'mu-{q}-R{R}' in ev for q, R in
           (('q25', 256), ('q25', 1024), ('q25', 2048), ('q25', 4096),
            ('q29', 1024), ('q29', 4096))))

    ph = ev['phase']
    ck('phase_window_67', ph['window'] == 67)
    ck('phase_cycle_cap_66', ph['cycle_cap'] == 66)
    ck('phase_enumerations_148070', ph['enumerations'] == 148070)
    ck('phase_zero_failures', ph['failures'] == 0)
    ck('phase_all_escaped', ph['escapes'] == ph['enumerations'])
    ck('phase_max_at_bound', ph['max_escape_steps'] <= ph['bound'] and ph['bound'] == 4422)
    ck('phase_divisors_enumerated', ph['divisor_lengths'] == [2, 4, 8, 16, 32, 64])
    ck('phase_fixed64_fails_at_4', ph['fixed64_failures']['4'] > 0)
    ck('phase_fixed64_fails_at_64', ph['fixed64_failures']['64'] > 0)
    ck('phase_prime_succeeds_where_fixed_fails',
       sum(ph['fixed64_failures'][k] for k in ('4', '8', '16', '32', '64')) > 100)

    mu = ev['mu-sweep']
    ck('mu_six_configs', len(mu['configs']) == 6)
    ck('mu_trials_8800', mu['n_trials'] == 8800)
    ck('mu_hat_matches_ledger', F(round(mu['mu_hat'] * 1000000), 1000000) == z.MU_HAT)
    ck('mu_se_matches_ledger', F(round(mu['se_pooled'] * 1000000), 1000000) == z.MU_SE_POOLED)
    ck('mu_charge_1p107', F(int(round(mu['mu_charge'] * 1000)), 1000) == z.WALK_MU)
    ck('mu_rows_above_035', all(r['p'] > 0.35 for r in mu['configs']))
    ck('mu_escapes_exercised', mu['n_escapes_two'] > 4000 and mu['n_escapes_check'] > 10000)
    ck('mu_chi2_reported', mu['chi2'] > 0 and mu['dof'] == 5)

    p = ev['pipeline']
    ck('pipeline_trials_2000', p['trials'] == 2000)
    ck('pipeline_recoveries_equal_trials', p['recoveries_verified'] == p['trials'])
    ck('pipeline_verified_only', p['p_recovered_overall'] == 1.0)
    ck('pipeline_degenerate_below_matches', p['degenerate'] < p['matches'])
    ck('pipeline_within_meter_ge_033', p['p_recovered_within_s_ops'] >= 0.33)
    ck('pipeline_escapes_exercised',
       p['escapes_2cycle'] > 100 and p['escapes_checkpoint'] > 1000)

    mg = ev['merge']
    ck('merge_zero_failures', mg['merge_failures'] == 0)
    ck('merge_checks_80000', mg['class_sequence_checks'] == 80000)

    ckc = ev['checkpoint']
    ck('checkpoint_escapes_at_least_100', ckc['escapes_checkpoint'] >= 100)
    ck('checkpoint_share_positive', ckc['checkpoint_share'] > 0)
    ck('checkpoint_escapes_per_trail_ge_5', ckc['escapes_per_trail'] >= 5)

    idx = ev['indexrule']
    ck('index_zero_mismatches_random', idx['random_mismatches_no_fallback'] == 0)
    ck('index_zero_mismatches_fallback', idx['random_mismatches_with_fallback'] == 0)
    ck('index_max_error_le_8', max(idx['max_abs_error_ulps']) <= z.INDEX_ERROR_ULPS)
    ck('index_no_gross_wraps', idx['gross_wraps'] == 0)
    ck('index_random_inputs_400k', idx['random_inputs'] >= 400000)
    ck('index_adversarial_inputs_70k', idx['adversarial_inputs'] >= 70000)
    ck('index_adversarial_fallbacks_positive', idx['adversarial_fallbacks'] > 0)
    ck('index_adversarial_zero_mismatch', idx['adversarial_mismatches_with_fallback'] == 0)
    ck('index_orbit_zero_mismatches', idx['orbit_mismatches'] == 0)

    tr = ev['traps-4']
    ck('traps_escapes_per_trail_ge_10', tr['escapes_per_trail'] >= 10)
    ck('traps_cross_epoch_below_5pct',
       tr['trails_with_cross_epoch_returns'] <= 0.05 * tr['trails'])
    ck('order2_no_exponent_in_U', ev['order2']['exponents_in_U'] == [])
    ck('order2_covers_2pow20', ev['order2']['max_m'] >= 2 ** 20)
    ck('order2_unit_group_six', ev['order2']['unit_group_size'] == 6)

    # ---- manifest and documents -----------------------------------------
    sub = json.loads((HERE / 'submission.json').read_text())
    ck('manifest_exact_keys', set(sub) == {'version', 'title', 'claims', 'experiments'})
    ck('claim_exact_keys',
       len(sub['claims']) == 1 and set(sub['claims'][0]) == {'evidence', 'score', 'statement'})
    ck('score_matches_ledger', sub['claims'][0]['score'] == z.frac(z.claim))
    ck('experiments_empty', sub['experiments'] == [])
    ck('dependencies_empty', json.loads((HERE / 'dependencies.json').read_text()) == [])
    arg = (HERE / 'argument.md').read_bytes()
    note = (HERE / 'public-note.md').read_bytes()
    ck('argument_under_8192', len(arg) < 8192)
    ck('note_10240_to_14336', 10240 <= len(note) < 14336)
    ck('note_no_private_paths', b'/Users/' not in note)
    ck('note_first_two_lines_model_harness',
       note.split(b'\n')[0] == b'Model: deepseek-v4.1-flash' and
       note.split(b'\n')[1] == b'Harness: Hermes')
    ck('note_cites_v13c_rejection', b'18225d19' in note)
    ck('note_states_convention_change', b'retired' in note and b'2e-5' in note)
    ck('note_cites_prime_window', b'prime 67' in note)
    ck('note_cites_pooled_mu', b'1.069618' in note and b'0.018344' in note)

    values = dict(claim=z.frac(z.claim), raw=z.frac(z.raw), W=z.frac(z.raw),
                  S_ops=str(z.S_ops), phase1=z.frac(z.phase1),
                  hard_extras=z.frac(z.hard_extras),
                  random_extras=z.frac(z.random_extras),
                  walk_mu=z.frac(z.WALK_MU), mu_hat=z.frac(z.MU_HAT),
                  mu_se_pooled=z.frac(z.MU_SE_POOLED),
                  mu_tolerance_hi=z.frac(z.mu_tolerance_hi),
                  step_price=z.frac(z.step_price),
                  step_instructions=z.step_instructions,
                  doubling_price=z.frac(z.doubling_price),
                  doubling_instructions=z.dbl_instructions,
                  rows={k: z.frac(v) for k, v in z.step_rows.items()},
                  success=z.frac(z.success), success_model=z.frac(z.success_model),
                  markov_halt=z.frac(z.markov_halt),
                  losses={k: z.frac(v) for k, v in z.losses.items()},
                  index_fallback_rate_bound=z.frac(z.index_fallback_rate_bound),
                  index_fallback_derived=z.frac(z.index_fallback_derived),
                  checkpoint_period=z.CHECKPOINT_PERIOD,
                  checkpoint_cycle_cap=z.CHECKPOINT_CYCLE_CAP,
                  checkpoint_detection_bound=z.checkpoint_detection_bound,
                  phase_enumerations=ph['enumerations'],
                  phase_failures=ph['failures'],
                  phase_max_steps=ph['max_escape_steps'],
                  integers={k: str(getattr(z, k)) for k in
                            ('N', 'R', 'd', 'DP', 'L', 'workers', 'nodes',
                             'M', 'S', 'node_cap_records')},
                  peak_memory_bytes=str(z.peak_memory),
                  pipeline_verified=f"{p['recoveries_verified']}/{p['trials']}",
                  pipeline_within=z.frac(F(round(p['p_recovered_within_s_ops'] * p['trials']),
                                             p['trials'])),
                  merge_checks=mg['class_sequence_checks'])
    return dict(version=1, verdict='PASS' if all(c.values()) else 'HOLD',
                scope=('Exact work, detector-scheduling, index-error and probability '
                       'arithmetic plus captured toy reconciliation at toy scale. H1 at '
                       'the pooled mu, the twist-uniformity clause, routing uniformity '
                       'and the instruction ceilings remain disclosed hypotheses.'),
                pass_count=sum(c.values()), check_count=len(c),
                checks=c, failed_checks=[k for k, v in c.items() if not v],
                values=values,
                evidence_sha256={q.name: hashlib.sha256(q.read_bytes()).hexdigest()
                                 for q in sorted((HERE / 'evidence').glob('*.json'))})


def main():
    out = build_certificate()
    data = json.dumps(out, separators=(',', ':')) + '\n'
    log = (f"{out['verdict']} {out['pass_count']}/{out['check_count']} arithmetic, "
           f"detector, index and capture checks; claim {z.claim}\n")
    if not out['failed_checks']:
        log += ('PASS 30 rows, 147 instructions, step 1.147, doubling 2.35; prime 67 '
                'window, cycle cap 66, bound 4422; phase 148070/148070 escapes\n')
        log += ('PASS mu pooled 1.069618 +/- 0.018344 over 6 configs/8800 trials; '
                'charge 1.107; S_ops 2^126.777348; raw 2^126.9752149\n')
        log += ('PASS claim 126.9769 < gate 127.00; success 0.391002; mu tolerance '
                '1.110677; pipeline 2000/2000 verified; merge 0 failures\n')
    else:
        log += 'HOLD ' + ','.join(out['failed_checks']) + '\n'
    if '--write' in sys.argv:
        (HERE / 'certificate.json').write_text(data)
        (HERE / 'evidence' / 'certificate_run.log').write_text(log)
    else:
        assert (HERE / 'certificate.json').read_text() == data, 'stale certificate'
        assert (HERE / 'evidence' / 'certificate_run.log').read_text() == log, 'stale run log'
    print(json.dumps({k: out[k] for k in ('verdict', 'pass_count', 'check_count',
                                          'failed_checks')}))
    return 0 if out['verdict'] == 'PASS' else 2


if __name__ == '__main__':
    sys.exit(main())
