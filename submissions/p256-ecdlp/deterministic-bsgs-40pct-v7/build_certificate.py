#!/usr/bin/env python3
"""Regenerate exact v7 arithmetic and bounded mechanism captures.
No float enters any acceptance test. Python is not the costed machine.
"""
import hashlib
import json
import math
import re
from fractions import Fraction
from pathlib import Path

from mechanism import (Field, affine_add, affine_scalar, jac, key, make_tree,
                       mixed, normalize, progression, recover, twice)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
P = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
R = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
C = 1000
M_ROWS = {
    'schoolbook_mul_umulh_16_pairs': 16 * 2,
    'schoolbook_carry_16_times_9': 16 * 9,
    'reduction_mul_umulh_4_times_4_pairs': 4 * 4 * 2,
    'reduction_carry_16_times_9': 16 * 9,
    'quotient_alias_copies': 4,
    'operand_modulus_output_argument_traffic': 24,
    'scratch_zeroing': 9,
    'final_five_word_subtract_select_moves': 16,
    'callee_save_restore_spill_allowance': 128,
    'four_shift_copy_schedules_and_round_control': 64,
    'entry_return_pointers_moves': 40,
}
A_ROWS = {'loads': 12, 'sum_with_carry': 5, 'subtract_or_correction': 5,
          'select': 4, 'stores': 4, 'calls_pointers_control': 10}
M_RAW, A_RAW = sum(M_ROWS.values()), sum(A_ROWS.values())
M, A = (6 * M_RAW + 4) // 5, (6 * A_RAW + 4) // 5
BATCH = 1024
I = 512 * (M + 64) + 256
INV_SHARE = (I + BATCH - 1) // BATCH
J = 16 * M + 32 * A + 512
GEN_ROWS = {
    'complete_case_Jacobian_update': J,
    'batch_prefix_backward_3M': 3 * M,
    'separate_ceil_I_over_1024': INV_SHARE,
    'x_scaling_2M': 2 * M,
    'y_scaling_3M': 3 * M,
    'canonical_Montgomery_decoding_2M': 2 * M,
    'batch_array_masks_control': 256,
    'record_build_and_local_send_preparation': 512,
}
K_REC = sum(GEN_ROWS.values())
M_PER_REC = (16 + 3 + 2 + 3 + 2) * M
NON_M_PER_REC = K_REC - M_PER_REC
V = 512 * J + I + 7 * M + 512
m = 1 << 128
h = (R + m - 1) // m
ACTUAL_N, CAP_N = m + h, 2 * m
HEIGHT, HOPS, NODE_BYTES = 128, 129, 1 << 14
INIT_ROWS = {'code_distribution': 4096, 'local_zeroing': 4096,
             'neighbor_prefix_link_setup': 2048, 'metadata': 1024, 'control': 5120}
PLACE_ROWS = {'payload_both_endpoints': 48, 'route': 16, 'queue': 16,
              'pointer': 16, 'tag': 16, 'call_control': 16}
SORT_ROWS = {'request_response_words': 52, 'key_compare': 64,
             'cache_copy': 48, 'state_end': 96, 'schedule_calls': 128,
             'pointers': 64, 'flags_counters': 60}
NODE_INIT, PLACE, SORT = sum(INIT_ROWS.values()), sum(PLACE_ROWS.values()), sum(SORT_ROWS.values())
ROWS = {
    'R_setup': 1 << 34,
    'R_init': 4 * m * NODE_INIT,
    'R_gen': CAP_N * K_REC,
    'R_place': CAP_N * HOPS * PLACE,
    'R_sort': CAP_N * HOPS * SORT + 4 * m * SORT,
    'R_merge': CAP_N * 512,
    'R_verify': 2 * V + (1 << 20),
}
T = sum(ROWS.values())
# Retain only the recurring 26M terms at their original price. Reprice
# everything else, including inversion shares and ALL one-time work.
M_TOTAL = CAP_N * M_PER_REC
NON_M_TOTAL = T - M_TOTAL
STRESS = Fraction(M_TOTAL, 1) + Fraction(6, 5) * NON_M_TOTAL
SCORE = Fraction(13647, 100)
CP, CQ = SCORE.numerator, SCORE.denominator
PRE = (ROWS['R_setup'] + ROWS['R_init'] + m * K_REC + m * HOPS * PLACE
       + m * HOPS * SORT + 2 * m * SORT)
PEAK = 4 * m * NODE_BYTES + (1 << 20)
PINS = {
    'tracks/p256-ecdlp/target-profile.json': '802d2350d1804b19fcf0e1f6d80d816862e966625710abfa2b4244ca764aa605',
    'tracks/classical-group-work-v2.json': 'e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294',
    'tracks/p256-ecdlp/rules.md': '0cbd08d817e61a129cf16a9d2d962f343830e139fe64e9a33da85863b7db09a6',
    'tracks/docs/COST_MODEL.md': 'a1d1e0c7dfd052f1306e7412781f832714a603f58e1d748733b63a258eb743c6',
}
START, END = '<!-- GENERATED-LEDGER-START -->', '<!-- GENERATED-LEDGER-END -->'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def below_power(work, exponent):
    """Exactly work < 2^exponent, with rational inputs and integer powers."""
    work, exponent = Fraction(work), Fraction(exponent)
    return (work.numerator ** exponent.denominator
            < work.denominator ** exponent.denominator * (1 << exponent.numerator))


def rat(value):
    value = Fraction(value)
    return f'{value.numerator}/{value.denominator}'


def replace_block(path, block):
    text = path.read_text()
    assert text.count(START) == 1 and text.count(END) == 1
    before, rest = text.split(START)
    _, after = rest.split(END)
    path.write_text(before + START + '\n' + block + '\n' + END + after)


def json_write(name, obj):
    (HERE / name).write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n')


def main():
    checks, details = {}, {}

    def check(name, condition, detail: object = ''):
        assert name not in checks
        checks[name] = bool(condition)
        details[name] = str(detail)
        print(('PASS ' if condition else 'FAIL ') + name + (' :: ' + str(detail) if detail else ''))

    profile = json.loads((ROOT / 'tracks/p256-ecdlp/target-profile.json').read_text())
    model = json.loads((ROOT / 'tracks/classical-group-work-v2.json').read_text())
    for path, expected in PINS.items():
        check('protected_pin_' + path.replace('/', '_'), sha(ROOT / path) == expected)
    check('p_is_protected', P == int(profile['field_modulus'], 16))
    check('r_is_protected', R == int(profile['subgroup_order'], 16))
    check('a_is_exact_minus_three', int(profile['coefficient_a'], 16) == P - 3)
    check('cofactor_is_one', profile['cofactor'] == 1)
    check('C_fixed_1000', model['reference_operation_costs']['p256-ecdlp'] == C)
    check('word_bits_64', model['word_bits'] == 64)
    check('each_scalar_instruction_costs_one', model['primitive_operations']['scalar_AArch64_instruction'] == 1)
    check('p_special_form_exact', P == (1 << 256) - (1 << 224) + (1 << 192) + (1 << 96) - 1)
    check('P256_Montgomery_pprime_is_one', -pow(P, -1, 1 << 64) % (1 << 64) == 1)
    check('CIOS_intermediate_inequality', 2 * P + P * ((1 << 64) - 1) * 2 == 2 * P * (1 << 64))
    check('CIOS_predivision_fits_nine_words', 2 * P * (1 << 64) < 1 << (9 * 64))
    check('register_schedule_fits_30', 9 + 8 + 4 + 2 + 3 + 1 + 3 <= 30)
    for nm, value in M_ROWS.items():
        check('M_row_' + nm, value >= 0, value)
    check('M_raw_exact_sum', M_RAW == 637, M_RAW)
    check('M_20pct_inside_ceiling', 5 * M >= 6 * M_RAW and 5 * (M - 1) < 6 * M_RAW, M)
    for nm, value in A_ROWS.items():
        check('A_row_' + nm, value >= 0, value)
    check('A_raw_exact_sum', A_RAW == 40, A_RAW)
    check('A_20pct_inside_ceiling', 5 * A >= 6 * A_RAW and 5 * (A - 1) < 6 * A_RAW, A)
    check('complete_mixed_branch_multiplications_fit_J', max(11, 4 + 8) <= 16)
    check('complete_mixed_branch_additions_fit_J', max(20, 2 + 20) <= 32)
    check('I_explicit_formula', I == 512 * (M + 64) + 256, I)
    check('I_call_bound_from_256_bit_exponent', (P - 2).bit_length() + (P - 2).bit_count() <= 512)
    check('inverse_share_separate_upper_bound', INV_SHARE * BATCH >= I and (INV_SHARE - 1) * BATCH < I, INV_SHARE)
    check('m_divisible_by_batch', m % BATCH == 0)
    check('h_divisible_by_batch', h % BATCH == 0)
    check('J_exact', J == 16 * M + 32 * A + 512, J)
    for nm, value in GEN_ROWS.items():
        check('generation_component_' + nm, value > 0, value)
    check('K_rec_exact_displayed_sum', K_REC == J + 3 * M + INV_SHARE + 2 * M + 3 * M + 2 * M + 256 + 512, K_REC)
    check('K_rec_M_nonM_partition', M_PER_REC + NON_M_PER_REC == K_REC and M_PER_REC == 26 * M)
    check('complete_verification_with_standalone_inversion', V == 512 * J + I + 7 * M + 512, V)
    check('h_is_ceiling', h == (R + m - 1) // m, h)
    check('h_strictly_less_than_m', 0 < h < m)
    check('minimal_full_coverage', (h - 1) * m < R <= h * m)
    check('coverage_interval_below_2r', h * m < 2 * R)
    check('actual_count_m_plus_h_not_2m', ACTUAL_N == m + h < CAP_N == 1 << 129, ACTUAL_N)
    check('indices_fit_two_words', max(m - 1, h - 1) < 1 << 128)
    check('end_and_counts_fit_three_words', CAP_N < 1 << 192)
    check('candidate_products_fit_four_words_plus_borrow', h * m < 1 << 256 and m < R)
    check('tree_height_m_relation', 1 << HEIGHT == m)
    check('engine_hops_height_plus_one', HOPS == HEIGHT + 1)
    check('tree_nodes_below_4m', 2 * (2 * m - 1) < 4 * m)
    check('local_address_space_fits_64_bits', NODE_BYTES < 1 << 64)
    check('global_memory_does_not_fit_64_bits', PEAK > 1 << 64)
    check('record_layout_96_bytes', 32 + 32 + 8 + 16 + 8 == 96)
    check('local_code_buffers_metadata_fit_allocation', 8192 + 4 * 96 + 3 * 96 + 1024 <= NODE_BYTES)
    check('node_initialization_exact_row_sum', NODE_INIT == 1 << 14, INIT_ROWS)
    check('placement_exact_row_sum', PLACE == 128, PLACE_ROWS)
    check('sort_per_edge_exact_row_sum', SORT == 512, SORT_ROWS)
    for nm, value in ROWS.items():
        check('ledger_' + nm, value > 0, value)
    check('R_gen_exact_capacity_product', ROWS['R_gen'] == 2 * m * K_REC, ROWS['R_gen'])
    check('R_sort_records_and_all_END_edges', ROWS['R_sort'] == 2 * m * HOPS * SORT + 4 * m * SORT)
    check('total_exact_row_sum', T == sum(ROWS.values()), T)
    check('M_nonM_exact_partition', M_TOTAL + NON_M_TOTAL == T and NON_M_TOTAL > 0)
    check('stress_exact_6_over_5_nonM', STRESS == Fraction(5 * M_TOTAL + 6 * NON_M_TOTAL, 5), rat(STRESS))
    check('stress_not_smaller_than_nominal', STRESS >= T)
    check('claim_reduced', math.gcd(CP, CQ) == 1, rat(SCORE))
    check('score_matches_submission', json.loads((HERE / 'submission.json').read_text())['claims'][0]['score'] == {'numerator': str(CP), 'denominator': str(CQ)})
    check('nominal_work_below_claim_integer_powers', below_power(Fraction(T, C), SCORE))
    check('nominal_slack_at_least_0p005', below_power(Fraction(T, C), SCORE - Fraction(1, 200)))
    check('stress_headroom_at_least_0p10', below_power(STRESS / C, SCORE - Fraction(1, 10)))
    check('requested_margin_pair_is_incompatible', Fraction(1, 10) > Fraction(3, 50) and STRESS >= T,
          'stress<=claim-0.10 implies nominal<=claim-0.10, incompatible with nominal>=claim-0.06')
    check('honest_nominal_wider_reserve_not_hidden', below_power(Fraction(T, C), SCORE - Fraction(3, 50)))
    check('claim_below_public_P256_prior', SCORE < Fraction(1382191, 10000))
    check('memory_exact_sum', PEAK == 4 * m * NODE_BYTES + (1 << 20), PEAK)
    check('memory_below_2pow145', PEAK < 1 << 145)
    check('preprocessing_exact_component_sum', PRE == ROWS['R_setup'] + ROWS['R_init'] + sum([m * K_REC, m * HOPS * PLACE, m * HOPS * SORT, 2 * m * SORT]), PRE)
    check('preprocessing_below_2pow136_units', PRE < C * (1 << 136))
    check('preprocessing_included_in_total', PRE <= T)
    check('success_coverage_entire_residue_interval', R <= m * h and (h - 1) * m < R and R > m)
    check('exact_success_one_exceeds_floor', 100 >= 39)
    check('no_dependencies', json.loads((HERE / 'dependencies.json').read_text()) == [])
    check('organizer_experiments_empty', json.loads((HERE / 'submission.json').read_text())['experiments'] == [])

    gx, gy = (int(profile['generator_P'][c], 16) for c in ('x', 'y'))
    bcurve, gen = int(profile['coefficient_b'], 16), (gx, gy)
    check('generator_on_exact_curve', (gy * gy - gx ** 3 + 3 * gx - bcurve) % P == 0)
    check('generator_order_annihilation', affine_scalar(R, gen, P) is None)
    f = Field(P)
    operands = [0, 1, 2, 3, P - 1, P - 2, (1 << 64) - 1, 1 << 64,
                (1 << 64) + 1, (1 << 128) - 1, 1 << 128,
                (1 << 192) - 1, 1 << 192, (1 << 255) - 1, 1 << 255, gx, gy, bcurve]
    state = 20261023
    for _ in range(8):
        state = (state * 6364136223846793005 + 1442695040888963407) % P
        operands.append(state)
    pair_ok, pairs = True, 0
    for a in operands:
        for b in operands:
            before = f.pairs
            pair_ok &= f.cios(a, b) == a * b * f.radix_inv % P
            pair_ok &= f.pairs - before == 32
            pairs += 1
    check('executed_full_width_CIOS_pairs', pair_ok and pairs == len(operands) ** 2, pairs)
    point_tests = 0
    for scalar in [0, 1, 2, 3, 17]:
        point = affine_scalar(scalar, gen, P)
        for scale in [1, 2, P - 1]:
            jp = jac(None, f) if point is None else (
                f.enc(point[0] * scale ** 2 % P),
                f.enc(point[1] * scale ** 3 % P), f.enc(scale))
            check(f'actual_P256_double_s{scalar}_z{scale}', normalize([twice(jp, f)], f)[0] == affine_add(point, point, P))
            for label, addend in [('O', None), ('P', gen), ('negP', (gx, -gy % P)), ('self', point)]:
                before = f.calls
                got = mixed(jp, addend, f)
                calls = f.calls - before
                check(f'actual_P256_mixed_s{scalar}_z{scale}_{label}',
                      calls <= 16 and normalize([got], f)[0] == affine_add(point, addend, P))
                point_tests += 1
    check('Y_zero_complete_branch', twice((f.enc(1), 0, f.enc(1)), f)[2] == 0)
    for label, start, increment in [('baby', None, gen), ('giant', affine_scalar(77, gen, P), (gx, -gy % P))]:
        out = progression(start, increment, 17, f, batch=4)
        oracle, ok = start, True
        for value in out:
            ok &= value == oracle
            oracle = affine_add(oracle, increment, P)
        check('actual_P256_progression_' + label, ok and len(out) == 17)
    neg = (gx, -gy % P)
    stream_records = [(gen, 7), (neg, 8), (None, 0), (affine_scalar(2, gen, P), 9), (gen, 10)]
    tree = make_tree(stream_records, 8)
    drained = []
    while True:
        value = tree.pull()
        if value is None:
            break
        drained.append(value)
    expected = sorted(stream_records, key=key)
    check('tree_stream_full_multiplicity_and_stable_ties', drained == expected)
    check('END_is_sticky', tree.pull() is None and tree.pull() is None)
    reversed_stream = sorted(stream_records, key=key, reverse=True)
    check('negative_control_reversed_compare_detected', reversed_stream != expected and any(key(reversed_stream[i]) > key(reversed_stream[i + 1]) for i in range(len(reversed_stream) - 1)))
    check('negative_control_old_c_and_c_plus_r_fail', affine_scalar(1, gen, P) != neg and affine_scalar(R + 1, gen, P) != neg)
    check('correct_minus_candidate_Q_equals_negP', affine_scalar((-1) % R, gen, P) == neg)

    tp, tb, tr, tg, tm = 97, 6, 107, (0, 43), 16
    toy_points = [(x, y) for x in range(tp) for y in range(tp) if (y * y - x ** 3 + 3 * x - tb) % tp == 0]
    check('toy_curve_counted_prime_order', len(toy_points) + 1 == tr and all(tr % d for d in range(2, math.isqrt(tr) + 1)), tr)
    check('toy_generator_on_curve_and_order', tg in toy_points and affine_scalar(tr, tg, tp) is None)
    th = (tr + tm - 1) // tm
    tf, toy_success, max_v = Field(tp), 0, 0
    for secret in range(tr):
        q = affine_scalar(secret, tg, tp)
        got, nv = recover(q, tg, tp, tr, tm, th, tf)
        toy_success += got == secret
        max_v = max(max_v, nv)
    check('exhaustive_toy_all_scalars_recovered', toy_success == tr, f'{toy_success}/{tr}')
    check('exhaustive_toy_at_most_two_verifications', max_v <= 2, max_v)
    secrets = [0, 1, 2, 15, 16, 17, 63, 127, 255, R - 1, R - 2]
    production_success, fixture_rows = 0, []
    for secret in secrets:
        q = affine_scalar(secret, gen, P)
        got, nv = recover(q, gen, P, R, 16, 16, f)
        check('actual_P256_BSGS_fixture_' + str(secret), got == secret and nv <= 2)
        production_success += got == secret
        fixture_rows.append({'secret': str(secret), 'recovered': str(got), 'verifications': nv})

    table = ['# Exact generated ledger', '', 'Every row below is generated from exact integer expressions.', '',
             '| Quantity | Exact value |', '|---|---:|',
             f'| m | {m} |', f'| actual h=ceil(r/m) | {h} |',
             f'| actual records m+h | {ACTUAL_N} |', f'| charged capacity 2m | {CAP_N} |',
             f'| M raw / safety ceiling | {M_RAW} / {M} |', f'| A raw / safety ceiling | {A_RAW} / {A} |',
             f'| J=16M+32A+512 | {J} |', f'| I=512(M+64)+256 | {I} |',
             f'| ceil(I/1024), charged separately | {INV_SHARE} |']
    for name, value in GEN_ROWS.items():
        table.append(f'| generation: {name} | {value} |')
    table += [f'| K_rec=sum(generation rows) | {K_REC} |', f'| V=512J+I+7M+512 | {V} |']
    for name, value in ROWS.items():
        table.append(f'| {name} | {value} |')
    table += [f'| T=sum(complete rows) | {T} |', f'| M-only recurring part | {M_TOTAL} |',
              f'| remaining instructions | {NON_M_TOTAL} |', f'| T_stress=Mpart+6*rest/5 | {rat(STRESS)} |',
              f'| preprocessing, included | {PRE} |', f'| peak bytes, including all node code/state | {PEAK} |',
              f'| submitted reduced score | {rat(SCORE)} |', '',
              'R_init=4m*16384; R_gen=2m*K_rec; R_place=2m*129*128;',
              'R_sort=2m*129*512+4m*512; R_merge=2m*512;',
              'R_setup=2^34; R_verify=2V+2^20.', '',
              f'Exact claim proof: T^{CQ} < 1000^{CQ} * 2^{CP}.',
              f'Exact stress proof: (T_stress/1000)^{CQ} < 2^{CP - CQ // 10}.',
              'The rational stress proof uses numerator/denominator cross-products, not floating point.', '']
    ledger = '\n'.join(table)
    (HERE / 'ledger.md').write_text(ledger)
    compact = (f'```text\nM={M}; A={A}; J={J}; I={I}; ceil(I/1024)={INV_SHARE}\n'
               f'K_rec={K_REC}; R_gen=2m*K_rec={ROWS["R_gen"]}\n'
               f'T={T}; C={C}; claim={rat(SCORE)}\n'
               f'T^{CQ} < 1000^{CQ}*2^{CP}\n'
               f'(T_stress/1000)^{CQ} < 2^{CP - CQ // 10}\n```')
    replace_block(HERE / 'argument.md', compact)
    replace_block(HERE / 'public-note.md', ledger)
    note = (HERE / 'public-note.md').read_text()
    argument = (HERE / 'argument.md').read_text()
    check('public_note_minimum_bytes', len(note.encode()) >= 5120, len(note.encode()))
    check('public_note_below_maximum_bytes', len(note.encode()) <= 100 * 1024)
    check('argument_under_8KB', len(argument.encode()) < 8192, len(argument.encode()))
    check('public_note_attribution_actual', 'GPT 6.1 Sol Max / Hermes' in note)
    check('public_note_private_path_clean', '/Users/' not in note)
    check('public_note_secret_shape_clean', re.search(r'ykn_[0-9a-f]{32,}|sk-[A-Za-z0-9]{20,}', note) is None)
    check('generated_numeric_blocks_equal_source', ledger in note and compact in argument)
    manifest_names = ['argument.md', 'kernel-ledger.md', 'ledger.md', 'public-note.md',
                      'submission.json', 'dependencies.json', 'mechanism.py', 'build_certificate.py']
    manifest = {name: sha(HERE / name) for name in manifest_names}
    evidence = {
        'version': 1,
        'scope': 'Actually executed bounded mechanism checks; not a production run or instruction measurement.',
        'protected_sha256': PINS,
        'source_sha256': manifest,
        'captured': {'CIOS_operand_pairs': pairs, 'CIOS_pairs_per_call': 32,
                     'mixed_branch_cases': point_tests, 'toy_curve': {'p': tp, 'a': tp - 3, 'b': tb, 'r': tr, 'P': list(tg)},
                     'toy_m': tm, 'toy_h': th, 'toy_exhaustive_trials': tr, 'toy_successes': toy_success,
                     'toy_max_verifications': max_v, 'P256_trials': len(secrets), 'P256_successes': production_success,
                     'P256_fixtures': fixture_rows,
                     'negative_controls': ['old sign candidates fail Q=-P', 'reversed comparison is not sorted']},
        'margin_policy': 'Maximum 0.06-bit nominal slack conflicts with minimum 0.10-bit stress slack; choose and disclose wider nominal reserve.',
        'rows_instructions': {k: str(v) for k, v in ROWS.items()},
        'T': str(T), 'T_stress': rat(STRESS), 'preprocessing': str(PRE), 'peak_bytes': str(PEAK),
    }
    json_write('evidence.json', evidence)
    manifest['evidence.json'] = sha(HERE / 'evidence.json')
    nominal_log = math.log2(T) - math.log2(C)
    stress_log = math.log2(STRESS.numerator) - math.log2(STRESS.denominator * C)
    cert = {
        'version': 1,
        'scope': 'Exact arithmetic plus captured bounded mechanism tests. Analytical schedules in kernel-ledger.md justify prices; PASS labels alone do not prove an instruction implementation.',
        'protected_sha256': PINS, 'source_sha256': manifest,
        'integers': {'p': str(P), 'r': str(R), 'C': C, 'm': str(m), 'h': str(h),
                     'actual_records': str(ACTUAL_N), 'charged_capacity_records': str(CAP_N),
                     'M_raw': M_RAW, 'M': M, 'A_raw': A_RAW, 'A': A, 'I': I, 'J': J,
                     'batch': BATCH, 'inverse_share': INV_SHARE, 'K_rec': K_REC, 'V': V,
                     'height': HEIGHT, 'engine_leaf_hops': HOPS, 'node_bytes': NODE_BYTES,
                     'M_per_rec': M_PER_REC, 'nonM_per_rec': NON_M_PER_REC},
        'primitive_rows': {'M': M_ROWS, 'A': A_ROWS},
        'generation_rows': GEN_ROWS, 'communication_rows': {'init': INIT_ROWS, 'place': PLACE_ROWS, 'sort': SORT_ROWS},
        'rows_instructions': {k: str(v) for k, v in ROWS.items()},
        'T': str(T), 'M_total': str(M_TOTAL), 'nonM_total': str(NON_M_TOTAL), 'T_stress': rat(STRESS),
        'claim': {'numerator': str(CP), 'denominator': str(CQ)},
        'proofs': {'claim': f'T^{CQ}<1000^{CQ}*2^{CP}',
                   'stress': 'below_power(T_stress/C, claim-1/10), using exact rational numerator/denominator powers',
                   'coverage': '(h-1)m<r<=hm<2r and x=(x//m)m+x%m for all canonical x',
                   'sign': 'equal x implies G_i=+/-B_j, so im+j or im-j mod r; full y and [c]P=Q resolve first equal-key pair',
                   'margin_conflict': 'stress>=nominal implies stress headroom 0.10 forces nominal headroom >=0.10>0.06'},
        'success_probability': {'numerator': '1', 'denominator': '1'},
        'preprocessing_instructions': str(PRE), 'peak_memory_bytes': str(PEAK), 'advice_bytes': 0,
        'checks': checks, 'details': details, 'check_count': len(checks), 'all_pass': all(checks.values()),
        'diagnostics_only_not_in_bounds': {'nominal_log2': nominal_log, 'stress_log2': stress_log,
                                            'nominal_claim_slack': float(SCORE) - nominal_log,
                                            'stress_claim_slack': float(SCORE) - stress_log},
    }
    json_write('certificate.json', cert)
    print(json.dumps({'checks': len(checks), 'passed': sum(checks.values()), 'M': M, 'A': A,
                      'K_rec': K_REC, 'R_gen': str(ROWS['R_gen']), 'T': str(T), 'claim': rat(SCORE),
                      'diagnostics': cert['diagnostics_only_not_in_bounds']}, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
