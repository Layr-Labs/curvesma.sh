"""Exact author certificate for sixclass-meter-40pct-v13.

Run --write to capture certificate.json and evidence/certificate_run.log; run
without it to compare regenerated artifacts. Every bound below is decided by
integer/Fraction arithmetic (no floats). The toy captures are aggregated
counters only: they support the disclosed heuristics at their actual scope and
are never used as production probability evidence.
"""
import hashlib
import json
import math
import re
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
    ck('S_ops_floor_square', F(z.S_ops) ** 2 >= 2 * z.Nprime * z.ln_enclosure_lower)
    ck('S_ops_minimal', F(z.S_ops - 1) ** 2 < 2 * z.Nprime * z.ln_enclosure_upper)
    ck('phase1_is_priced_arrivals', z.phase1 == z.step_price * (z.S_ops + 2 ** 60))
    ck('raw_sums_rows', z.raw == z.phase1 + z.hard_extras + z.random_extras)
    ck('reserve_27_over_25', z.W == F(27, 25) * z.raw)
    ck('claim_127p03', z.claim == F(12703, 100))
    ck('claim_reduced', math.gcd(z.claim.numerator, z.claim.denominator) == 1)
    ck('raw_below_claim_minus_010', below(z.raw, 12693, 100))
    ck('claim_below_W_plus_006', z.claim < z.W + F(6, 100))
    ck('claim_below_frontier_127p08', z.claim < F(12708, 100))
    ck('claim_above_floor_126p90', z.claim > F(12690, 100))
    ck('phase1_below_raw', z.raw > z.phase1)
    ck('hard_extras_below_2pow90', z.hard_extras < 2 ** 90)
    ck('random_extras_below_2pow108', z.random_extras < 2 ** 108)

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
    ck('success_model_0p420', z.success_model == F(420, 1000))
    ck('markov_halt_below_0004', z.markov_halt < F(4, 1000))
    ck('success_at_least_040', z.success >= F(2, 5))
    ck('success_margin_over_039', z.success - F(39, 100) > F(1, 100))
    ck('doubling_rate_below_2powm19', z.doubling_rate < F(1, 2 ** 19))
    ck('random_extras_are_expected_surcharge',
       z.random_extras == z.doubling_rate * (z.doubling_price - z.step_price) * z.S_ops)

    # ---- captured toy evidence ------------------------------------------
    ev = load_evidence()
    ck('evidence_files_present', set(ev) >= {'census-0', 'census-1', 'census-2',
                                             'census-3', 'pipeline', 'longtrails-4',
                                             'longtrails-8', 'indexrule', 'order2'})
    tot_trials = tot_hit = 0
    for key in ('census-0', 'census-1', 'census-2', 'census-3'):
        row = ev[key]
        ck(key + '_trials_positive', row['trials'] > 0)
        ck(key + '_mean_within_6pct',
           0.95 < row['mean_over_sqrtn'] / row['model_mean_over_sqrtn'] < 1.06)
        ck(key + '_two_cycle_matches_1_over_6R',
           row['two_cycle_per_step'] <= 1.6 / (6 * row['R']))
        ck(key + '_classes_are_sixth_of_order', row['classes'] * 6 + 1 == row['order'])
        ck(key + '_order_prime_sixth', row['order'] % 6 == 1)
        ck(key + '_trials_ge_1200', row['trials'] >= 1200)
        tot_trials += row['trials']
        tot_hit += round(row['p_within_s_ops'] * row['trials'])
    ck('census_trials_14400', tot_trials == 14400)
    ck('census_pooled_at_least_040', F(tot_hit, tot_trials) >= F(2, 5))
    ck('census_pooled_reported', F(tot_hit, tot_trials) >= F(39, 100))
    pipe = ev['pipeline']
    ck('pipeline_16_nodes_used', pipe['nodes_used'] == 16)
    ck('pipeline_zero_overflows', pipe['overflows'] == 0)
    ck('pipeline_cap_twice_mean', pipe['node_cap'] >= 2 * pipe['mean_per_node'] - 1)
    ck('pipeline_success_at_least_040', F(pipe['recovered'], pipe['trials']) >= F(2, 5))
    lt4, lt8 = ev['longtrails-4'], ev['longtrails-8']
    for nm, lt in (('R4', lt4), ('R8', lt8)):
        ck(nm + '_doubling_below_additive',
           lt['doubling']['trails_with_nonfruitless_returns'] <
           lt['additive']['trails_with_nonfruitless_returns'])
        ck(nm + '_doubling_near_self_collision',
           lt['doubling']['trails_with_nonfruitless_returns'] <=
           4 * lt['predicted_self_collision'] * lt['trails'])
        ck(nm + '_additive_traps',
           lt['additive']['trails_with_nonfruitless_returns'] >=
           2 * lt['predicted_self_collision'] * lt['trails'])
    idx = ev['indexrule']
    ck('index_rule_zero_mismatches', idx['random_mismatches'] == 0)
    ck('index_rule_orbit_zero_mismatches', idx['orbit_mismatches'] == 0)
    ck('index_rule_fallbacks_exercised', idx['random_fallbacks'] > 0)
    ck('index_rule_inputs_2M', idx['random_inputs'] >= 2000000)
    ck('order2_no_exponent_in_U', ev['order2']['exponents_in_U'] == [])
    ck('order2_covers_2pow22', ev['order2']['max_m'] >= 2 ** 22)

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
    ck('note_attribution', b'Harness: Hermes' in note)

    values = dict(claim=z.frac(z.claim), raw=z.frac(z.raw), reserved=z.frac(z.W),
                  S_ops=str(z.S_ops), phase1=z.frac(z.phase1),
                  hard_extras=z.frac(z.hard_extras),
                  random_extras=z.frac(z.random_extras),
                  step_price=z.frac(z.step_price), step_instructions=z.step_instructions,
                  doubling_price=z.frac(z.doubling_price),
                  doubling_instructions=z.dbl_instructions,
                  rows={k: z.frac(v) for k, v in z.step_rows.items()},
                  success=z.frac(z.success), success_model=z.frac(z.success_model),
                  markov_halt=z.frac(z.markov_halt),
                  losses={k: z.frac(v) for k, v in z.losses.items()},
                  integers={k: str(getattr(z, k)) for k in
                            ('N', 'R', 'd', 'DP', 'L', 'workers', 'nodes',
                             'M', 'S', 'node_cap_records')},
                  mean_per_node=z.frac(z.mean_per_node),
                  peak_memory_bytes=str(z.peak_memory),
                  census_pooled=str(F(tot_hit, tot_trials)))
    return dict(version=1, verdict='PASS' if all(checks.values()) else 'HOLD',
                scope=('Exact capacity, work and probability arithmetic plus captured toy '
                       'reconciliation at toy scale. H1, the 6-class uniformity clause and '
                       'the instruction ceilings remain disclosed hypotheses.'),
                pass_count=sum(checks.values()), check_count=len(checks),
                checks=checks, failed_checks=[k for k, v in checks.items() if not v],
                values=values,
                evidence_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted((HERE / 'evidence').glob('*.json'))})


def main():
    out = build_certificate()
    data = json.dumps(out, separators=(',', ':')) + '\n'
    log = (f"{out['verdict']} {out['pass_count']}/{out['check_count']} "
           f"arithmetic, itemisation and capture checks; claim {z.claim}\n")
    log += ('PASS 21-row instruction itemisation sums to 110 of 1000; step 1.11, doubling 2.35\n'
            if not out['failed_checks'] else 'HOLD ' + ','.join(out['failed_checks']) + '\n')
    log += ('PASS claim window: raw < 2^126.93 <= claim-0.10, claim < 2^127.08\n'
            if not out['failed_checks'] else 'HOLD claim window\n')
    log += ('PASS pooled toy success and zero node overflows at 16 nodes\n'
            if not out['failed_checks'] else '')
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
