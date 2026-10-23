"""Exact v8 instruction ledger. Python is a verifier, not the priced machine."""
from fractions import Fraction

P = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
R = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
C = 1000
SCORE = Fraction(2721, 20)  # 136.05; retain the predecessor's stress reserve.
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
M_PER_REC = 26 * M
V = 512 * J + I + 7 * M + 512
SINGLE_NORM = I + 7 * M + 512
m = 1 << 127
L = 2 * m  # 4^64; the signed 64-digit lattice is all odd integers |b|<L.
h = (R + L - 1) // L
ACTUAL_N, CAP_N = m + h, m + L
BABY_HEIGHT, GIANT_HEIGHT = 127, 128
BABY_HOPS, GIANT_HOPS = BABY_HEIGHT + 1, GIANT_HEIGHT + 1
NODES = 2 * CAP_N  # Strict upper bound on nodes/END edges in both trees.
EDGE_RECORDS = m * BABY_HOPS + L * GIANT_HOPS
NODE_BYTES = 1 << 14
INIT_ROWS = {'code_distribution': 4096, 'local_zeroing': 4096,
             'neighbor_prefix_link_setup': 2048, 'metadata': 1024, 'control': 5120}
PLACE_ROWS = {'payload_both_endpoints': 48, 'route': 16, 'queue': 16,
              'pointer': 16, 'tag': 16, 'call_control': 16}
SORT_ROWS = {'request_response_words': 52, 'key_compare': 64,
            'cache_copy': 48, 'state_end': 96, 'schedule_calls': 128,
            'pointers': 64, 'flags_counters': 60}
NODE_INIT, PLACE, SORT = sum(INIT_ROWS.values()), sum(PLACE_ROWS.values()), sum(SORT_ROWS.values())
BABY_EXTRA_ROWS = {'increment_table_loads': 8, 'two_word_index_odd_label': 16,
                   'two_entry_selection_pointers': 16, 'loop_and_control': 24}
BABY_EXTRA = sum(BABY_EXTRA_ROWS.values())
# Three charged scalar actions: [1]P, [3]P and P0=[m]P.
# Four complete updates: derive 2P=3P-P, 2P0, 4P0, and Q-2P0.
# Four standalone normalizations (no batch amortization), negations/control.
SMALL_EXPANDED = 3 * V + 4 * J + 4 * SINGLE_NORM + (1 << 20)
SMALL_CHARGED = 4 * V + 4 * J + (1 << 20)
ROWS = {
    'R_setup': 1 << 34,
    'R_small_table_and_seeds': SMALL_CHARGED,
    'R_init': NODES * NODE_INIT,
    'R_gen': CAP_N * K_REC,
    'R_baby_lookup_labels': m * BABY_EXTRA,
    'R_place': EDGE_RECORDS * PLACE,
    'R_sort': EDGE_RECORDS * SORT + NODES * SORT,
    'R_merge': CAP_N * 512,
    'R_verify': 2 * V + (1 << 20),
}
T = sum(ROWS.values())
# Conservatively stress ALL one-time work, including its M calls, too.
M_TOTAL = CAP_N * M_PER_REC
NON_M_TOTAL = T - M_TOTAL
STRESS = Fraction(5 * M_TOTAL + 6 * NON_M_TOTAL, 5)
PRE = (ROWS['R_setup'] + ROWS['R_small_table_and_seeds'] + ROWS['R_init']
       + m * (K_REC + BABY_EXTRA + BABY_HOPS * (PLACE + SORT)) + 2 * m * SORT)
PEAK = NODES * NODE_BYTES + (1 << 20)
PINS = {
    'tracks/p256-ecdlp/target-profile.json': '802d2350d1804b19fcf0e1f6d80d816862e966625710abfa2b4244ca764aa605',
    'tracks/classical-group-work-v2.json': 'e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294',
    'tracks/p256-ecdlp/rules.md': '0cbd08d817e61a129cf16a9d2d962f343830e139fe64e9a33da85863b7db09a6',
    'tracks/docs/COST_MODEL.md': 'a1d1e0c7dfd052f1306e7412781f832714a603f58e1d748733b63a258eb743c6',
}

def below_power(work, exponent):
    work, exponent = Fraction(work), Fraction(exponent)
    return (work.numerator ** exponent.denominator
            < work.denominator ** exponent.denominator * (1 << exponent.numerator))

def rational(value):
    value = Fraction(value)
    return f'{value.numerator}/{value.denominator}'

def ledger_text():
    vals = {'positive odd babies m': m, 'L=2m=4^64': L, 'actual h=ceil(r/L)': h,
            'actual records m+h': ACTUAL_N, 'charged records m+L=3*2^127': CAP_N,
            'M raw': M_RAW, 'M ceiling': M, 'A raw': A_RAW, 'A ceiling': A,
            'J=16M+32A+512': J, 'I=512(M+64)+256': I,
            'ceil(I/1024), separate': INV_SHARE, **{'generation: '+k:v for k,v in GEN_ROWS.items()},
            'K_rec exact sum': K_REC, 'baby extra lookup/label': BABY_EXTRA,
            'V=512J+I+7M+512': V, 'expanded table/seed cost': SMALL_EXPANDED,
            'charged table/seed envelope': SMALL_CHARGED,
            'baby / giant hops': f'{BABY_HOPS} / {GIANT_HOPS}',
            'charged nodes / END edges': NODES, **ROWS,
            'T exact instruction sum': T, 'M-only recurring part': M_TOTAL,
            'remaining instructions': NON_M_TOTAL, 'T_stress=Mpart+6*rest/5': rational(STRESS),
            'preprocessing included': PRE, 'peak bytes including code': PEAK,
            'submitted reduced score': rational(SCORE)}
    lines = ['# Exact generated v8 ledger', '', '| Quantity | Exact value |', '|---|---:|']
    lines += [f'| {k} | {v} |' for k,v in vals.items()]
    lines += ['', 'H=0; W=T/1000. R_gen=(m+L)*23121; baby extra=64m.',
              'R_init=2(m+L)*16384; E=128m+129L; R_place=128E;',
              'R_sort=512E+2(m+L)*512; R_merge=(m+L)*512.',
              'R_setup=2^34; R_small_table_and_seeds=4V+4J+2^20; R_verify=2V+2^20.',
              'Exact claim: T^20 < 1000^20*2^2721.',
              'Exact stress: (T_stress/1000)^20 < 2^2719 (claim minus 0.10).',
              'No floating logarithm enters an acceptance comparison.', '']
    return '\n'.join(lines)
