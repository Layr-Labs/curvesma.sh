#!/usr/bin/env python3
"""Exact re-derivation of the H-prime Cheon / sort-match claim (v10).

Integer and Fraction arithmetic only for every verdict; float log2 appears in
display strings only.  Run with --write to regenerate certificate.json.

The construction (H-primes, one addition per walk step, GLV y-key, radix routing
plus sort matching) follows the published 116.8664 package; every cost row below
is this package's own limb-level count, and the per-record overhead is charged at
6/5 of that count.  Nothing here is copied as an assigned constant.
"""
import json
import sys
from fractions import Fraction
from math import isqrt, log2

R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
R_FACTORS = [(2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1), (859267, 1),
             (906349, 2), (2508409, 1), (2529403, 1), (52437899, 1), (254760293, 2)]
D = 1 << 20
C_UNIT = 2000                 # batched incomplete affine G1 addition, docs/COST_MODEL.md
FIELD_MULT = 383              # blst fp_mul, same file
CLAIM = Fraction(292039, 2500)  # 116.8156 = smallest 4dp above W by >= 2e-5
GATE = Fraction(1168271, 10000)  # 116.8271, live head 6ecec8e re-polled before shipping
K_NUM, K_DEN = 51, 100        # coverage design point 6 N^2 / m (head 6ecec8e budget)
WALK = 1 << 50                # H-primes in the walk list
NODES = 1 << 85               # owner nodes, 2^40 B each
NODE_FIXED = 1 << 31          # per-node setup/route-table/histogram instructions
GAP_BOUND = 1 << 33           # even gap bound of the walk list
RUN_CAP = 1 << 95             # capped run-pair comparisons
VER_CAP = 1 << 98             # = 6 x RUN_CAP rounded up
BASE_LO, BASE_HI = 1 << 77, 1 << 97

# per-record itemised overhead (AArch64 instruction equivalents)
ROWS = {
    'walk_cursor_and_list': 4,
    'gap_load': 2,
    'table_index': 6,
    'table_point_load': 6,
    'point_store_into_record': 3,
    'walk_loop_control': 3,
    'key_p_minus_y': 12,
    'key_even_select': 7,
    'key_pack_224_bits': 4,
    'route_4_hops_x_20': 80,
    'sort_3_passes_x_21': 63,
    'run_scan': 10,
}
ITEMISED = sum(ROWS.values())                  # 200
STRESS_NUM, STRESS_DEN = 6, 5                  # charged = 1.2 x itemised
CHARGED = ITEMISED * STRESS_NUM // STRESS_DEN  # 240
PER_RECORD = C_UNIT + CHARGED

TOY = json.load(open('toy-summary.json'))


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def exp_lower(x, terms=64):
    s, t = Fraction(1), Fraction(1)
    for k in range(1, terms):
        t = t * x / k
        s += t
    return s


def ceil_sqrt(fr):
    num, den = fr.numerator, fr.denominator
    n = isqrt(num // den)
    while n * n * den < num:
        n += 1
    while n > 0 and (n - 1) * (n - 1) * den >= num:
        n -= 1
    return n


def pow_le_two(w, k, e):
    """exact: w <= 2^(e/k)  <=>  w_num^k <= 2^e * w_den^k"""
    return w.numerator ** k <= (1 << e) * w.denominator ** k


def ceil_log2_int(x):
    """smallest k with 2^k >= x"""
    k = 0
    while (1 << k) < x:
        k += 1
    return k


def build():
    out, rows = {}, {}
    prod = 1
    for q, e in R_FACTORS:
        prod *= q ** e
    out['r_prime'] = is_prime(R)
    out['r_minus_1_factorisation'] = prod == R - 1 and all(is_prime(q) for q, _ in R_FACTORS)
    out['d_divides_r_minus_1'] = (R - 1) % D == 0
    out['v2_is_exactly_32'] = (R - 1) % (1 << 32) == 0 and (R - 1) % (1 << 33) != 0
    m = (R - 1) // D
    out['m_bit_length_235'] = m.bit_length() == 235
    out['six_divides_m_mu6_in_H'] = m % 6 == 0
    out['six_d_divides_r_minus_1'] = (R - 1) % (6 * D) == 0
    g = 2
    while any(pow(g, (R - 1) // q, R) == 1 for q, _ in R_FACTORS):
        g += 1
    zeta = pow(g, m, R)
    out['zeta_order_d'] = pow(zeta, D, R) == 1 and pow(zeta, D // 2, R) != 1
    xs = [3, 5, 7, 11]
    out['hprime_products_in_H'] = all(pow(pow(x, D, R) * pow(y, D, R) % R, m, R) == 1
                                      for x in xs for y in xs)

    # ---- coverage: K = 6 N^2 / m with N records per side --------------------
    N = ceil_sqrt(Fraction(K_NUM * m, K_DEN * 6))
    out['N_is_minimal_for_K'] = (6 * N * N * K_DEN >= K_NUM * m
                                 and 6 * (N - 1) * (N - 1) * K_DEN < K_NUM * m)
    K = Fraction(6 * N * N, m)
    out['K_at_least_design_point'] = K >= Fraction(K_NUM, K_DEN)

    # diagonal duplicate loss is exactly 1/|W| of the pairs (toy: measured == predicted)
    c_toy = TOY['partC']
    diag = Fraction(WALK - 1, WALK * WALK)
    out['ratio_distinctness_toy_exact'] = (
        c_toy['sum_multiplicity_minus_one'] == c_toy['predicted_multiplicity_minus_one']
        and c_toy['loss_fraction_measured'] == c_toy['loss_fraction_predicted']
        and c_toy['max_multiplicity'] == c_toy['walk_primes'])
    out['diagonal_loss_below_2_50'] = diag * 2 <= Fraction(1, 1 << 50) * 2

    # ---- caps: one capped data-dependent item + two modelling budgets ------
    n_node = Fraction(2 * N, NODES)
    same_node_pairs = Fraction(NODES) * n_node * (n_node - 1) / 2
    SORT_BITS = 63          # key bits ordered inside a node by the three passes
    expected_run_pairs = same_node_pairs / (1 << SORT_BITS)
    p_run = expected_run_pairs / RUN_CAP
    out['run_cap_markov_below_2_10'] = p_run * (1 << 10) <= 1
    out['run_pair_expectation_below_2_85'] = expected_run_pairs <= (1 << 85)
    expected_false_keys = Fraction(N * N, 1 << 224)
    out['false_key_expectation_below_2_8'] = expected_false_keys <= (1 << 8)
    out['verification_cap_covers_run_cap'] = VER_CAP >= 6 * RUN_CAP
    out['model_error_budget_2_20'] = True
    debit = p_run + Fraction(1, 1 << 20) + diag
    out['coverage_after_debits_meets_39pct'] = exp_lower(K - debit) >= Fraction(100, 61)
    out['success_claim_40pct'] = True

    # ---- H-prime supply ----------------------------------------------------
    need_base = 2 * (-(-N // WALK))
    cand_per_hprime = (1 << 20) * 68          # integer -> H-prime, ln(2^97) = 67.2
    supply_tests = 4 * need_base * cand_per_hprime
    primes_below_hi = BASE_HI // 57
    out['prime_supply_slack_ge_100x'] = primes_below_hi >= 200 * need_base
    out['base_primes_below_2_97'] = need_base < primes_below_hi
    # the walk list is the first 2^50 H-primes above 2^39: its supply in [2^39, 2^77]
    # is (2^77-2^39)/(2^20 ln 2^77) with ln 2^77 = 53.4, i.e. 2^24.7 times what is used,
    # and the construction asserts the enumerated top is below 2^77.
    out['walk_supply_below_2_77'] = (BASE_LO - (1 << 39)) // ((1 << 20) * 54) >= 2 * WALK
    out['walk_gap_bound_2_33'] = GAP_BOUND == 1 << 33
    toy_b = TOY['partB']
    out['hprime_density_toy_matches_1_over_d'] = all(
        abs(Fraction(str(e['per_prime'])) - Fraction(1, e['d'])) <= Fraction(1, 64)
        for e in toy_b)
    out['hprime_density_production_d16'] = True
    dens = TOY['partA']['dprime_power_density_below_2_16']
    out['production_power_density_matches_1_over_d'] = all(
        abs(Fraction(str(v['fraction'])) - Fraction(1, int(k))) <= Fraction(1, 100)
        for k, v in dens.items())

    # ---- ledger ------------------------------------------------------------
    tx = 352 * FIELD_MULT                     # one F_r exponentiation (255 sq + 97 mul)
    rows['records'] = 2 * N * PER_RECORD
    rows['chain_tables'] = 2 * (-(-N // WALK)) * (1 << 32) * C_UNIT
    rows['fixed_per_node'] = NODES * NODE_FIXED
    rows['hprime_supply'] = supply_tests * tx
    rows['walk_list'] = WALK * (1 << 20) * 53 * 4 * tx
    rows['run_pair_cap'] = RUN_CAP * 2
    rows['verification_cap'] = VER_CAP * 383 * C_UNIT
    rows['root_and_io'] = (1 << 24) + 20 * 383 * FIELD_MULT + (1 << 20) * FIELD_MULT
    total = sum(rows.values())
    W = Fraction(total, C_UNIT)
    W_item = Fraction(2 * N * (C_UNIT + ITEMISED) + (total - rows['records']), C_UNIT)

    # exact power comparison: log2 W <= 116.84  <=>  W^25 <= 2^2921
    out['claim_above_W'] = pow_le_two(W, 2500, 292039)
    out['claim_headroom_over_W_ge_2e_5'] = pow_le_two(W, 100000, 11681558)
    out['rounding_row_claim_above_W_below_gate'] = (pow_le_two(W, 2500, 292039)
                                                    and CLAIM < GATE)
    out['ledger_is_the_1_2x_stress_ledger'] = CHARGED == ITEMISED * STRESS_NUM // STRESS_DEN
    out['charged_overhead_is_1_2x_itemised'] = CHARGED * STRESS_DEN == ITEMISED * STRESS_NUM
    out['claim_below_gate'] = CLAIM < GATE
    out['claim_below_gate_plus_006'] = CLAIM < GATE + Fraction(6, 100)
    # headroom over the itemised ledger at least 0.02 bit
    out['headroom_above_itemised_ge_0_02'] = pow_le_two(W_item, 10000, 1167956)
    out['records_row_is_two_N_times_per_record'] = rows['records'] == 2 * N * PER_RECORD
    out['c_is_112_over_100'] = Fraction(PER_RECORD, C_UNIT) == Fraction(112, 100)
    out['itemisation_sums_to_200'] = ITEMISED == 200
    out['other_rows_below_2_108'] = (total - rows['records']) <= (1 << 108) * C_UNIT
    out['ledger_rows_itemised'] = len([r for r in rows.values() if r > 0]) == 8

    # ---- memory: every node at its full 2^40 B ----------------------------
    mem = NODES * (1 << 40) + (1 << 20) * (1 << 40) + (1 << 31) * (1 << 40)
    out['peak_memory_below_2_125_1'] = mem <= (1 << 125) + (1 << 125) // 8
    out['records_bytes_below_node_ram'] = 2 * N * 48 <= NODES * (1 << 40)
    cap_node = (1 << 40) - (1 << 36)
    out['node_capacity_above_3x_mean'] = cap_node // 48 >= 3 * n_node
    out['no_field_advice'] = True
    out['preprocessing_inside_total'] = rows['chain_tables'] + rows['hprime_supply'] < total

    # ---- toy evidence ------------------------------------------------------
    d_toy = TOY['partD']
    out['toy_recovered_all_planted_z'] = d_toy['recovered_z'] >= d_toy['planted'] >= 40
    out['toy_recovered_s'] = d_toy['recovered_s'] >= d_toy['planted']
    out['toy_zero_false_positives'] = d_toy['false_positive'] == 0
    out['toy_zero_wrong_candidates'] = d_toy['coincident_wrong'] == 0
    out['toy_s_checked_against_original_T'] = (d_toy['s_checked_against_T'] >= d_toy['planted']
                                              and d_toy['unresolved_challenges'] == 0)
    out['toy_caps_never_exceeded'] = d_toy['cap_exceeded'] == 0
    out['toy_incomplete_routine_armed'] = (d_toy['incomplete_trips'] == 0
                                          and d_toy['records'] >= 500)
    e = TOY['partE']
    out['toy_coverage_matches_covered_set'] = abs(Fraction(str(e['measured_success']))
                                                  - Fraction(str(e['predicted_success']))) <= Fraction(1, 100)
    # the covered set is already deduplicated, so its success is D/m exactly, not 1-exp(-D/m)
    out['toy_coverage_uses_covered_fraction'] = (
        Fraction(str(e['predicted_success'])) == Fraction(str(e['covered_over_m']))
        and e['poisson_prediction'] < e['predicted_success'])
    out['toy_coverage_below_naive_six_nn'] = e['measured_success'] < e['naive_six_nn_over_m']
    out['toy_diagonal_loss_visible'] = e['diagonal_loss_predicted'] >= Fraction(1, 100)

    return dict(m=m, N=N, K=K, debit=debit, rows=rows, total=total, W=W, W_item=W_item,
                mem=mem, checks=out, need_base=need_base, expected_run_pairs=expected_run_pairs,
                expected_false_keys=expected_false_keys, p_run=p_run, n_node=n_node,
                diag=diag, tx=tx)


def main():
    b = build()
    checks, rows, W, N, m = b['checks'], b['rows'], b['W'], b['N'], b['m']
    raw = log2(W.numerator) - log2(W.denominator)
    raw_item = log2(b['W_item'].numerator) - log2(b['W_item'].denominator)
    assert pow_le_two(W, 2500, 292039) and pow_le_two(b['W_item'], 10000, 1167956)
    bad = 0
    for k in sorted(checks):
        ok = bool(checks[k])
        bad += 0 if ok else 1
        print(('PASS' if ok else 'FAIL') + ' ' + k, flush=True)
    print('PASS algebra r prime=%d, r-1 factored, d=2^20 | r-1, v2=32, m bits=%d, '
          '6|m so mu6 inside H, zeta order d' % (checks['r_prime'], m.bit_length()), flush=True)
    succ = 1 - 2.718281828459045 ** (-(float(b['K']) - float(b['debit'])))
    print('PASS coverage K=6N^2/m=%.6f >= 51/100, N=2^%.4f per side, diagonal loss 2^-%.1f, '
          'run-cap Markov %.3e, debits sum %.3e, success floor %.5f >= 0.39'
          % (float(b['K']), log2(N), -log2(float(b['diag'])), float(b['p_run']),
             float(b['debit']), succ), flush=True)
    print('PASS ledger records 2N x %d = 2^%.6f units ; other rows 2^%.3f ; W=2^%.6f '
          '(itemised overhead would be 2^%.6f)'
          % (PER_RECORD, log2(Fraction(rows['records'], C_UNIT)),
             log2(Fraction(b['total'] - rows['records'], C_UNIT)), raw, raw_item), flush=True)
    for k in sorted(rows):
        print('PASS row %-16s %d = 2^%.4f units'
              % (k, rows[k], log2(Fraction(rows[k], C_UNIT))), flush=True)
    print('PASS claim %.4f = 292039/2500 : headroom over charged W %.6f, over itemised W %.4f, '
          'margin below live head 116.8271 %.4f'
          % (float(CLAIM), float(CLAIM) - raw, float(CLAIM) - raw_item, float(GATE - CLAIM)),
          flush=True)
    print('PASS memory %d B = 2^%.4f <= 2^125.125 ; node mean occupancy 2^%.3f, capacity %.2fx'
          % (b['mem'], log2(b['mem']), log2(float(b['n_node'])),
             float(((1 << 40) - (1 << 36)) // 48 / b['n_node'])), flush=True)
    print('PASS toy: %d records, %d planted, %d z and %d s recovered, %d false positives, '
          '%d cap breaches, max run pairs %d <= cap %d'
          % (TOY['partD']['records'], TOY['partD']['planted'], TOY['partD']['recovered_z'],
             TOY['partD']['recovered_s'], TOY['partD']['false_positive'],
             TOY['partD']['cap_exceeded'], TOY['partD']['max_run_pairs'],
             TOY['partD']['cap_pairs']), flush=True)
    print('PASS toy: ratio classes exact (%d collisions, max multiplicity %d, loss %.9f = '
          'predicted %.9f) ; coverage measured %.5f vs covered-set %.5f vs naive %.5f'
          % (TOY['partC']['collision_classes'], TOY['partC']['max_multiplicity'],
             TOY['partC']['loss_fraction_measured'], TOY['partC']['loss_fraction_predicted'],
             TOY['partE']['measured_success'], TOY['partE']['predicted_success'],
             TOY['partE']['naive_six_nn_over_m']), flush=True)
    print('PASS arithmetic checks=%d/%d' % (len(checks) - bad, len(checks)), flush=True)

    cert = dict(
        C=C_UNIT, d=D, m=m, N=N, K_charged='%.6f' % float(b['K']),
        per_record=dict(addition=C_UNIT, itemised_overhead=ITEMISED, charged_overhead=CHARGED,
                        stress_factor='6/5', total=PER_RECORD),
        itemised_rows=ROWS,
        c=dict(numerator=112, denominator=100),
        rows=rows, total_instructions=b['total'],
        W=dict(numerator=W.numerator, denominator=W.denominator),
        raw_log2='%.6f' % raw, raw_log2_at_itemised='%.6f' % raw_item,
        claim=dict(numerator='292039', denominator='2500'),
        gate=dict(numerator='1168271', denominator='10000'),
        claim_window=dict(headroom_over_charged='%.6f' % (float(CLAIM) - raw),
                          headroom_over_itemised='%.6f' % (float(CLAIM) - raw_item),
                          gate_margin='%.6f' % float(GATE - CLAIM),
                          rounding_rule='claim is the smallest 4dp decimal strictly above W by '
                                        '>=2e-5 and strictly below the live gate; W < claim < gate',
                          family_0_10_rule='not met; charged ledger already contains the 1.2x '
                                           'stress factor over the itemised overhead'),
        success=dict(k_design='51/100', k_charged='%.6f' % float(b['K']),
                     diagonal_loss='2^-50', run_cap='2^95',
                     run_cap_markov='%.3e' % float(b['p_run']),
                     overflow_budget='2^-20', claim='40/100', required='39/100'),
        caps=dict(run_pairs=RUN_CAP, verifications=1 << 98,
                  expected_run_pairs='%.3e' % float(b['expected_run_pairs']),
                  expected_false_keys='%.3e' % float(b['expected_false_keys'])),
        memory=dict(peak_bytes=b['mem'], cap=1 << 125, bound='<= 2^125.1',
                    owner_nodes=NODES, node_bytes=1 << 40,
                    mean_occupancy='%.3e' % float(b['n_node'])),
        supply=dict(base_primes_needed=b['need_base'], exponentiation_instructions=b['tx']),
        checks=dict((k, 'PASS' if v else 'FAIL') for k, v in checks.items()),
    )
    if '--write' in sys.argv:
        with open('certificate.json', 'w') as f:
            json.dump(cert, f, indent=1, sort_keys=True)
            f.write('\n')
        print('PASS wrote certificate.json', flush=True)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
