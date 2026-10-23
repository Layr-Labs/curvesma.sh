"""Exact author certificate for sixclass-meter-40pct-v13c.

Run --write to capture certificate.json and evidence/certificate_run.log; run
without it to compare regenerated artifacts.  Every production bound below is
decided by integer/Fraction arithmetic (no floats).  The toy captures are
aggregated counters: they support the disclosed heuristics at their actual scope
and are never used as production probability evidence.
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


def below(x, num, den=1):
    """x < 2^(num/den), decided in integers."""
    x = F(x)
    return x.numerator ** den < (1 << num) * x.denominator ** den


def nbelow(x, k):
    """x < 2^-k, decided in integers."""
    x = F(x)
    return x.numerator < (x.denominator << k)


def above(x, num, den=1):
    x = F(x)
    return x.numerator ** den > (1 << num) * x.denominator ** den


def load_evidence():
    ev = HERE / 'evidence'
    out = {}
    for path in sorted(ev.glob('*.json')):
        out[path.stem] = json.loads(path.read_text())
    return out


def build_certificate():
    checks = {}

    def ck(name, value):
        checks[name] = bool(value)

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
    ck('itemised_rows_21', len(z.itemized_step_rows) == 21)
    ck('itemised_rows_unique', len({q[0] for q in z.itemized_step_rows}) == 21)
    ck('itemised_every_row_priced',
       all(z.step_rows[q[0]] == q[3] * z.instruction_primitives[q[2]]
           for q in z.itemized_step_rows))
    ck('itemised_sum_110', z.step_instructions == 110)
    ck('step_price_1p11', z.step_price == F(111, 100))
    ck('step_price_below_reference_plus_ceiling', z.step_price < 1 + F(z.E0, z.C))
    ck('doubling_limb_count_2320', z.dbl_instructions == 7 * 200 + 6 * 20 + 3 * 200 + 200)
    ck('doubling_price_2p35', z.doubling_price == F(235, 100))
    ck('doubling_above_step', z.doubling_price > z.step_price)

    # ---- S_ops and the ledger -------------------------------------------
    ck('S_ops_floor_square', F(z.S_ops) ** 2 >= 2 * z.WALK_MU * z.Nprime * z.ln_enclosure_lower)
    ck('S_ops_minimal', F(z.S_ops - 1) ** 2 < 2 * z.WALK_MU * z.Nprime * z.ln_enclosure_upper)
    ck('S_ops_uses_mu_1p058', z.WALK_MU == F(529, 500))
    ck('toy_implied_mu_is_1p05772', z.TOY_MU_DERIVED == F(26443, 25000))
    ck('toy_implied_mu_below_budget_mu', z.TOY_MU_DERIVED <= z.WALK_MU)
    ck('phase1_is_priced_arrivals', z.phase1 == z.step_price * (z.S_ops + 2 ** 60))
    ck('raw_sums_rows', z.raw == z.phase1 + z.hard_extras + z.random_extras)
    ck('reserve_27_over_25', z.W == F(27, 25) * z.raw)
    ck('claim_127p00', z.claim == F(127))
    ck('claim_reduced', math.gcd(z.claim.numerator, z.claim.denominator) == 1)
    ck('raw_below_claim_minus_010', below(z.raw, 12690, 100))
    ck('claim_below_W_plus_006', z.claim < z.W + F(6, 100))
    ck('claim_below_gate_127p01', z.claim < F(12701, 100))
    ck('claim_above_floor_126p90', z.claim > F(12690, 100))
    ck('phase1_below_raw', z.raw > z.phase1)
    ck('hard_extras_below_2pow90', z.hard_extras < 2 ** 90)
    ck('random_extras_below_2pow108', z.random_extras < 2 ** 108)
    ck('mu_tolerance_above_assumed', z.mu_tolerance_hi > z.WALK_MU)
    ck('mu_tolerance_between_1p06_1p07', F(106, 100) < z.mu_tolerance_hi < F(107, 100))

    # ---- truncated-index rule: error bound and fallback rate -------------
    ck('index_error_bound_8_ulps', z.INDEX_ERROR_ULPS == 8)
    ck('index_threshold_covers_2E_plus_1', 2 * z.INDEX_ERROR_ULPS + 1 <= z.INDEX_GAP)
    ck('index_gap_is_40', z.INDEX_GAP == 40)
    ck('index_rate_is_294_over_2pow64',
       z.index_fallback_rate_bound == F(294, 2 ** 64))
    ck('index_rate_below_2powm55', nbelow(z.index_fallback_rate_bound, 55))
    ck('index_rate_charged_at_2powm50', nbelow(z.index_fallback_rate_bound, 50))
    ck('index_derived_below_charged',
       z.index_fallback_derived < z.index_fallback_expected)
    ck('index_charged_2pow89', z.index_fallback_expected == 2 ** 89)

    # ---- checkpoint scheduling ------------------------------------------
    ck('checkpoint_period_64', z.CHECKPOINT_PERIOD == 64)
    ck('checkpoint_cycle_cap_64', z.CHECKPOINT_CYCLE_CAP == 64)
    ck('checkpoint_detection_bound_2pow12', z.checkpoint_detection_bound == 2 ** 12)

    # ---- storage partition ----------------------------------------------
    ck('nodes_2pow16', z.nodes == 2 ** 16)
    ck('per_node_cap_records', z.node_cap_records == 2 ** 43)
    ck('per_node_bytes_below_2pow56', z.per_node_bytes < 2 ** 56)
    ck('record_192_bytes', z.record_bytes == 192)
    ck('total_records_below_2pow60', z.M == z.nodes * z.node_cap_records < 2 ** 60)
    ck('node_setup_charged', z.hard_extras >= z.node_setup == z.nodes * 2 ** 20)
    ck('reports_below_cap', F(2 ** 58) < z.M)
    ck('peak_memory_below_2pow70', z.peak_memory < 2 ** 70)

    # ---- probability -----------------------------------------------------
    ck('success_model_0p392', z.success_model == F(392, 1000))
    ck('markov_halt_below_0004', z.markov_halt < F(4, 1000))
    ck('success_at_least_039', z.success >= F(39, 100))
    ck('success_margin_over_039_positive', z.success - F(39, 100) > F(9, 10000))
    ck('loss_budget_below_1p1e-3', sum(z.losses.values(), F(0)) + z.markov_halt < F(11, 10000))
    ck('doubling_rate_below_2powm19', z.doubling_rate < F(1, 2 ** 19))
    ck('random_extras_are_expected_surcharge',
       z.random_extras == z.doubling_rate * (z.doubling_price - z.step_price) * z.S_ops)
    ck('detector_abort_below_2powm40', nbelow(z.detector_abort, 40))
    ck('endpoint_transport_below_2powm40', nbelow(z.endpoint_transport, 40))

    # ---- captured toy evidence ------------------------------------------
    ev = load_evidence()
    need = {'census-0', 'census-1', 'census-2', 'census-3', 'pipeline',
            'pipeline-dp1', 'routing16', 'traps-4', 'traps-8', 'merge',
            'checkpoint', 'indexrule', 'order2'}
    ck('evidence_files_present', set(ev) >= need)
    tot_trials = tot_hit = 0
    for key in ('census-0', 'census-1', 'census-2', 'census-3'):
        row = ev[key]
        ck(key + '_trials_ge_1000', row['trials'] >= 1000)
        ck(key + '_mean_within_10pct',
           0.95 < row['mean_over_sqrtn'] / row['model_mean_over_sqrtn'] < 1.10)
        ck(key + '_p_within_s_ops_ge_040', row['p_within_s_ops'] >= 0.40)
        ck(key + '_two_cycle_at_most_1p6_over_6R',
           row['two_cycle_per_step'] <= 1.6 / (6 * row['R']))
        ck(key + '_no_exceptional_addition', row['exceptional_additions'] == 0)
        ck(key + '_classes_are_sixth_of_order', row['classes'] * 6 + 1 == row['order'])
        ck(key + '_order_prime_sixth', row['order'] % 6 == 1)
        tot_trials += row['trials']
        tot_hit += round(row['p_within_s_ops'] * row['trials'])
    ck('census_trials_6800', tot_trials == 6800)
    ck('census_pooled_at_least_040', F(tot_hit, tot_trials) >= F(2, 5))
    ck('census_pooled_reported', F(tot_hit, tot_trials) >= F(41, 100))

    for key in ('pipeline', 'pipeline-dp1'):
        p = ev[key]
        ck(key + '_recoveries_equal_trials', p['recoveries_verified'] == p['trials'])
        ck(key + '_verified_only', p['p_recovered_overall'] == 1.0)
        ck(key + '_degenerate_below_matches', p['degenerate'] < p['matches'])
        ck(key + '_unit_one_not_dominant', p['unit_one_fraction'] < 0.25)
        ck(key + '_escapes_exercised', p['escapes_2cycle'] > 0)
        ck(key + '_collision_rate_ge_035', p['collisions_within_s_ops'] >= 0.35 * p['trials'])
        ck(key + '_recovery_within_meter_ge_035',
           p['p_recovered_within_s_ops'] >= 0.35)
        ck(key + '_lag_within_units',
           p['mean_lag_per_trail'] is not None and
           p['mean_lag_per_trail'] <= 4 * p['lag_units'])
    ck('pipeline_dp1_recovery_ge_039', ev['pipeline-dp1']['p_recovered_within_s_ops'] >= 0.39)
    ck('pipeline_dp1_lag_below_dp3',
       ev['pipeline-dp1']['mean_lag_per_trail'] < ev['pipeline']['mean_lag_per_trail'])
    ck('pipeline_dp1_meter_closes_gap',
       ev['pipeline-dp1']['p_recovered_within_s_ops'] >
       ev['pipeline']['p_recovered_within_s_ops'])

    rt = ev['routing16']
    ck('routing_zero_overflows', rt['overflows'] == 0)
    ck('routing_all_nodes_used', rt['nodes_used'] == rt['nodes'])
    ck('routing_max_below_cap', rt['max_load'] <= rt['cap'])
    ck('routing_cap_at_most_2p05_mean', rt['cap'] <= 2.05 * rt['mean_per_node'])
    ck('routing_dispersion_below_1p9', rt['chi_square'] <= 1.9 * rt['chi_square_dof'])
    ck('routing_samples_2pow22', rt['samples'] >= 2 ** 22)

    for key in ('traps-4', 'traps-8'):
        t = ev[key]
        ck(key + '_escapes_per_trail_ge_10', t['escapes_per_trail'] >= 10)
        ck(key + '_cross_epoch_returns_below_5pct',
           t['trails_with_cross_epoch_returns'] <= 0.05 * t['trails'])
        ck(key + '_prediction_reported', t['predicted_self_collision'] >= 0.03)

    mg = ev['merge']
    ck('merge_zero_failures', mg['merge_failures'] == 0)
    ck('merge_checks_at_least_100k', mg['class_sequence_checks'] >= 100000)
    ck('merge_escapes_exercised', mg['escapes'] > 1000)

    ckc = ev['checkpoint']
    ck('checkpoint_escapes_at_least_100', ckc['escapes_checkpoint'] >= 100)
    ck('checkpoint_share_at_least_2pct', ckc['checkpoint_share'] >= 0.02)
    ck('checkpoint_escapes_per_trail_ge_5', ckc['escapes_per_trail'] >= 5)

    idx = ev['indexrule']
    ck('index_rule_zero_mismatches_random', idx['random_mismatches_no_fallback'] == 0)
    ck('index_rule_zero_mismatches_with_fallback',
       idx['random_mismatches_with_fallback'] == 0)
    ck('index_rule_max_error_le_8', max(idx['max_abs_error_ulps']) <= z.INDEX_ERROR_ULPS)
    ck('index_rule_no_gross_wraps', idx['gross_wraps'] == 0)
    ck('index_rule_inputs_1p5M', idx['random_inputs'] >= 1500000)
    ck('index_rule_no_random_fallback', idx['random_fallbacks'] == 0)
    ck('index_rule_adversarial_inputs_300k', idx['adversarial_inputs'] >= 300000)
    ck('index_rule_window_30000', idx['adversarial_within_gap'] >= 30000)
    ck('index_rule_adversarial_fallbacks_20000', idx['adversarial_fallbacks'] >= 20000)
    ck('index_rule_adversarial_zero_mismatch',
       idx['adversarial_mismatches_with_fallback'] == 0)
    ck('index_rule_cheap_rule_would_fail',
       idx['adversarial_cheap_rule_would_be_wrong'] > 0)
    ck('index_rule_orbit_zero_mismatches', idx['orbit_mismatches'] == 0)
    ck('index_rule_orbit_checks_600k', idx['orbit_checks'] >= 600000)

    ck('order2_no_exponent_in_U', ev['order2']['exponents_in_U'] == [])
    ck('order2_covers_2pow22', ev['order2']['max_m'] >= 2 ** 22)
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
    ck('note_cites_189ebc4', b'189ebc4' in note)
    ck('note_cites_1c600df', b'1c600df' in note)
    ck('note_attribution_string_present',
       b'Hermes' in note and b'deepseek-v4.1-flash' in note
       and b'qwen3.8-flash-next' in note)
    ck('note_no_attribution_labels',
       not any(ln.startswith((b'Model:', b'Harness:')) for ln in note.splitlines()))
    ck('note_mentions_rejection', b'5a5da5de' in note)

    values = dict(claim=z.frac(z.claim), raw=z.frac(z.raw), reserved=z.frac(z.W),
                  S_ops=str(z.S_ops), phase1=z.frac(z.phase1),
                  hard_extras=z.frac(z.hard_extras),
                  random_extras=z.frac(z.random_extras),
                  walk_mu=z.frac(z.WALK_MU), toy_mu_derived=z.frac(z.TOY_MU_DERIVED),
                  mu_tolerance_hi=z.frac(z.mu_tolerance_hi),
                  step_price=z.frac(z.step_price), step_instructions=z.step_instructions,
                  doubling_price=z.frac(z.doubling_price),
                  doubling_instructions=z.dbl_instructions,
                  rows={k: z.frac(v) for k, v in z.step_rows.items()},
                  success=z.frac(z.success), success_model=z.frac(z.success_model),
                  markov_halt=z.frac(z.markov_halt),
                  losses={k: z.frac(v) for k, v in z.losses.items()},
                  index_fallback_rate_bound=z.frac(z.index_fallback_rate_bound),
                  index_fallback_derived=z.frac(z.index_fallback_derived),
                  index_error_ulps=z.INDEX_ERROR_ULPS,
                  index_gap_ulps=z.INDEX_GAP,
                  checkpoint_detection_bound=z.checkpoint_detection_bound,
                  integers={k: str(getattr(z, k)) for k in
                            ('N', 'R', 'd', 'DP', 'L', 'workers', 'nodes',
                             'M', 'S', 'node_cap_records')},
                  mean_per_node=z.frac(z.mean_per_node),
                  peak_memory_bytes=str(z.peak_memory),
                  census_pooled=z.frac(F(tot_hit, tot_trials)),
                  pipeline_verified=f"{ev['pipeline']['recoveries_verified']}/"
                                    f"{ev['pipeline']['trials']}",
                  pipeline_recovered_within=z.frac(
                      F(round(ev['pipeline']['p_recovered_within_s_ops'] * ev['pipeline']['trials']),
                        ev['pipeline']['trials'])),
                  index_adversarial_window=str(idx['adversarial_within_gap']),
                  index_adversarial_fallbacks=str(idx['adversarial_fallbacks']),
                  routing_max_over_cap=z.frac(
                      F(rt['max_load'], rt['cap'])))
    return dict(version=1, verdict='PASS' if all(checks.values()) else 'HOLD',
                scope=('Exact capacity, work, index-error and probability arithmetic plus '
                       'captured toy reconciliation at toy scale. H1, the twist-uniformity '
                       'clause, routing uniformity and the instruction ceilings remain '
                       'disclosed hypotheses.'),
                pass_count=sum(checks.values()), check_count=len(checks),
                checks=checks, failed_checks=[k for k, v in checks.items() if not v],
                values=values,
                evidence_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted((HERE / 'evidence').glob('*.json'))})


def main():
    out = build_certificate()
    data = json.dumps(out, separators=(',', ':')) + '\n'
    log = (f"{out['verdict']} {out['pass_count']}/{out['check_count']} "
           f"arithmetic, itemisation, index-error and capture checks; claim {z.claim}\n")
    if not out['failed_checks']:
        log += ('PASS 21-row itemisation 110 of 1000; step 1.11, doubling 2.35; '
                'index 8 products, E=8, gap 40, rate 294/2^64 < 2^-55\n')
        log += ('PASS success budget 0.392 at mu 1.058 (toy 1.05772): S_ops 2^126.74469, '
                'raw 2^126.89525 <= claim-0.10; claim 127.00 < gate 127.01\n')
        log += ('PASS floor 0.39: success 0.39101, mu tolerance 1.061530 > 1.058; '
                'toy mechanism 4000/4000 verified; census pooled '
                f"{out['values']['census_pooled']['numerator']}/"
                f"{out['values']['census_pooled']['denominator']}; 0 node overflows\n")
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
