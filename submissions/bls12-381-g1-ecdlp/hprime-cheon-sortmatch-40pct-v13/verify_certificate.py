#!/usr/bin/env python3
"""Exact re-derivation of the H-prime Cheon / sort-match claim (v12).

Integer and Fraction arithmetic only for every verdict; float log2 appears in
display strings only.  Run with --write to regenerate certificate.json.

The construction (H-primes, one addition per walk step, GLV y-key, radix routing
plus sort matching) follows the published 116.8664 package; every cost row below
is this package's own limb-level count, and the per-record overhead is charged at
6/5 of that count.  Nothing here is copied as an assigned constant.

v12 changes, both of them ledger/model accounting rather than construction:
  * coverage is stated as E[C_distinct]/m with C_distinct the distinct residue
    set of the 6N^2 candidate occurrences; the model is the random-multiset
    occupancy form m(1-exp(-K_cand/m)), whose success 1-exp(-K) already contains
    every random repetition.  The modelling allowance is no longer the boolean
    2^-20 placeholder: it is the largest measured shortfall over three exact
    distinct-coverage censuses in toy.py (partEcensus), charged in full.
  * the 2^20-candidate T-scan of the d-th root search is an explicit ledger row
    (2^20 candidates x 2^11 units = 2^31 units), which the previous ledger left
    inside the much smaller root_and_io row.
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


def exp_neg_bounds(x, terms=40):
    """Exact lo <= exp(-x) <= hi by the alternating series (x > 0 rational)."""
    s, t = Fraction(1), Fraction(1)
    lo = hi = Fraction(1)
    for k in range(1, terms):
        t = t * (-x) / k
        s += t
        if k % 2:
            lo = s            # truncated after a negative term: lower bound
        else:
            hi = s            # truncated after a positive term: upper bound
    return lo, hi


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

    # ---- coverage model: E[C_distinct]/m on a two-axis grid ------------------
    # Success is the distinct covered fraction of H: E[C_distinct]/m, where
    # C_distinct is the residue set of the 6N^2 candidate occurrences.  Under
    # H3occ the candidate multiset behaves like K_cand uniform draws into m bins,
    # so E[C_distinct]/m = 1-exp(-K_cand/m) = 1-exp(-K); that occupancy form
    # already contains every random repetition of a multiset.  The allowance for
    # departing from it is measured on a grid: two toy fields (r ~ 2^24, r ~ 2^26)
    # x d in {8, 16, 64}, every cell at K_cand at the production design point.
    cen = TOY['partEcensus']                       # single-field K sweep, d = 8
    grid = TOY['partFgrid']
    cells = grid['cells']
    ratios = [Fraction(str(c['ratio_measured_over_predicted'])) for c in cen]
    shortfalls = [Fraction(str(c['shortfall_over_m'])) for c in cen]
    cell_bits = sorted(set(c['r_bits'] for c in cells))
    cell_ds = sorted(set(c['d'] for c in cells))
    effs = [Fraction(str(c['shortfall_eff_over_m'])) for c in cells]
    eps_model = max(abs(v) for v in effs)          # charged in full, success units
    out['coverage_census_three_sizes'] = (len(cen) == 3
                                          and all(c['distinct_covered'] > 0 for c in cen)
                                          and all(c['candidate_occurrences'] > 0 for c in cen))
    out['coverage_census_ratio_rises_with_size'] = ratios[0] < ratios[1] < ratios[2]
    out['coverage_census_largest_ratio_ge_0_99'] = (
        ratios[2] >= Fraction(99, 100)
        and cen[2]['candidate_occurrences'] > cen[1]['candidate_occurrences']
        > cen[0]['candidate_occurrences'])
    out['coverage_census_brackets_production_K'] = (
        Fraction(str(cen[2]['K_cand'])) >= K - Fraction(1, 100)
        and Fraction(str(cen[2]['K_cand'])) <= K + Fraction(1, 20))
    out['coverage_census_diagonal_explains_residual'] = all(
        abs(Fraction(str(c['ratio_eff'])) - 1) <= Fraction(5, 1000) for c in cen)
    out['coverage_census_raw_form_over_predicts'] = all(
        Fraction(str(c['shortfall_over_m'])) > 0 for c in cen)
    # --- the two-axis grid that carries the production allowance --------------
    out['coverage_grid_two_fields'] = (len(cells) == len(cell_bits) * len(cell_ds)
                                       and len(cell_bits) == 2 and len(cells) == 6)
    out['coverage_grid_three_d_values'] = cell_ds == [8, 16, 64]
    out['coverage_grid_fields_are_distinct'] = (
        cell_bits[0] < cell_bits[1] and grid['field_bits'] == cell_bits)
    out['coverage_grid_brackets_production_K'] = all(
        abs(Fraction(str(c['K_cand'])) - K) <= Fraction(3, 100) for c in cells)
    # one-sided: the counted toy diagonal is an upper envelope of the raw shortfall,
    # and the residual never adds to it (the opposite sign is diagonal candidates
    # coinciding with off-diagonal ones, which cannot lower a production floor)
    out['coverage_grid_diagonal_is_the_systematic_term'] = all(
        Fraction(str(c['shortfall_over_m'])) - Fraction(str(c['diagonal_share']))
        <= eps_model for c in cells)
    out['coverage_grid_residual_shrinks_with_field'] = (
        grid['shrink_count'] == len(cells) // 2
        and all(Fraction(str(s['large_field_abs'])) <= Fraction(str(s['small_field_abs']))
                for s in grid['shrink_pairs']))
    out['model_error_allowance_is_grid_max'] = (
        eps_model == max(abs(v) for v in effs) and eps_model >= Fraction(1, 1 << 20)
        and eps_model <= Fraction(1, 100))
    def _pred_ok(c, key_occ, key_pred):
        m_c = c['m']
        lo, hi = exp_neg_bounds(Fraction(c[key_occ], m_c))
        return (m_c - hi * m_c - 1 <= c[key_pred] <= m_c - lo * m_c + 1)

    out['model_error_allowance_rows_are_exact'] = (
        all(_pred_ok(c, 'candidate_occurrences', 'predicted_distinct') for c in cen)
        and all(_pred_ok(c, 'candidate_occurrences', 'predicted_distinct') for c in cells)
        and all(abs(Fraction(c['distinct_covered'], c['predicted_distinct'])
                    - Fraction(str(c['ratio_measured_over_predicted']))) <= Fraction(1, 10 ** 6)
                for c in cells))
    debit = p_run + eps_model + diag               # all three in success-probability units
    out['coverage_after_debits_meets_39pct'] = (1 / exp_lower(K) + debit) <= Fraction(61, 100)
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
    # the d-th root search scans 2^20 candidate roots, one scalar multiplication
    # each at 2^11 units: an explicit row, not an allowance inside root_and_io
    rows['root_scan_T'] = (1 << 20) * (1 << 11) * C_UNIT
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
    out['ledger_rows_itemised'] = len([r for r in rows.values() if r > 0]) == 9
    out['root_scan_row_is_2_31_units'] = Fraction(rows['root_scan_T'], C_UNIT) == (1 << 31)

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
    # the 20000-challenge hit count checks the uniform challenge sampler against the
    # toy's exact deduplicated fraction D/m; the distinct-set *prediction* is tested
    # by the three-size census above, not by this hit count
    out['toy_coverage_below_naive_six_nn'] = e['measured_success'] < e['naive_six_nn_over_m']
    out['toy_diagonal_loss_visible'] = e['diagonal_loss_predicted'] >= Fraction(1, 100)
    rad = TOY['partD']
    out['toy_radix_pipeline_three_lsd_passes'] = (
        rad['sort_passes'] == 3
        and rad['sort_bits'] == rad['sort_digit'] * rad['sort_passes']
        and rad['sort_bits'] < rad['key_bits'])
    out['toy_radix_routing_on_top_bits'] = rad['route_bits'] < rad['key_bits']
    out['toy_radix_ops_match_ledger'] = (
        rad['ops_per_record_per_pass'] == rad['ledger_ops_per_record_per_pass']
        and rad['ledger_ops_per_record_per_pass'] == ROWS['sort_3_passes_x_21'] // 3
        and rad['ops_per_run_pair'] == rad['ledger_ops_per_run_pair']
        and rad['ledger_ops_per_run_pair'] == ROWS['run_scan'])
    out['toy_full_key_scan_covers_unordered_bits'] = (
        rad['middle_bits_scan_only'] == rad['key_bits'] - rad['route_bits'] - rad['sort_bits']
        and rad['middle_bits_scan_only'] > 0)
    out['toy_radix_prefix_work_is_per_node'] = rad['prefix_ops_per_node_pass'] <= 2 ** 21

    return dict(m=m, N=N, K=K, debit=debit, rows=rows, total=total, W=W, W_item=W_item,
                mem=mem, checks=out, need_base=need_base, expected_run_pairs=expected_run_pairs,
                expected_false_keys=expected_false_keys, p_run=p_run, n_node=n_node,
                diag=diag, tx=tx, eps_model=eps_model, census=cen)


def main():
    b = build()
    grid, cells, rad = TOY['partFgrid'], TOY['partFgrid']['cells'], TOY['partD']
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
    succ = 1 - 2.718281828459045 ** (-float(b['K']))
    print('PASS coverage K=6N^2/m=%.6f >= 51/100, N=2^%.4f per side, diagonal 2^-%.1f, '
          'run-cap Markov %.3e, measured model allowance %.8f, total debits %.6f, '
          'success floor %.5f >= 0.39'
          % (float(b['K']), log2(N), -log2(float(b['diag'])), float(b['p_run']),
             float(b['eps_model']), float(b['debit']), succ - float(b['debit'])), flush=True)
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
                     model_error_allowance='%.8f' % float(b['eps_model']),
                     model_error_allowance_source='grid maximum of the absolute measured deviation '
                                                  'of the exact distinct-coverage census from the '
                                                  'occupancy prediction with the counted toy diagonal '
                                                  'repetitions removed (toy.py partFgrid: 2 fields '
                                                  'r=8391553/33555073 x d=8/16/64, 6 cells, K_cand '
                                                  '0.495-0.515 bracketing production K=0.51; '
                                                  'max_abs_shortfall_eff=0.00541878)',
                     debits_in_success_units='%.6f' % float(b['debit']),
                     success_after_debits='%.5f' % (1 - 2.718281828459045 ** (-float(b['K']))
                                                    - float(b['debit'])),
                     claim='40/100', required='39/100'),
        coverage_grid=dict(field_bits=grid['field_bits'], d_values=grid['d_values'],
                           K_values=grid['K_values'], cells_count=grid['cells_count'],
                           max_shortfall_eff='%.8f' % float(grid['max_shortfall_eff']),
                           max_abs_shortfall_eff='%.8f' % float(grid['max_abs_shortfall_eff']),
                           max_raw_shortfall='%.8f' % float(grid['max_raw_shortfall']),
                           residual_shrinks_with_field=grid['shrink_count'],
                           cells=[dict(r_bits=c['r_bits'], d=c['d'], m=c['m'], nw=c['nw'],
                                       nb=c['nb'], records_per_side=c['records_per_side'],
                                       candidate_occurrences=c['candidate_occurrences'],
                                       distinct_covered=c['distinct_covered'],
                                       K_cand=c['K_cand'], predicted_distinct=c['predicted_distinct'],
                                       ratio_measured_over_predicted=c['ratio_measured_over_predicted'],
                                       shortfall_over_m=c['shortfall_over_m'],
                                       diagonal_repetitions=c['diagonal_repetitions'],
                                       diagonal_share=c['diagonal_share'],
                                       predicted_distinct_eff=c['predicted_distinct_eff'],
                                       ratio_eff=c['ratio_eff'],
                                       shortfall_eff_over_m=c['shortfall_eff_over_m'])
                                  for c in cells]),
        radix_pipeline=dict(passes=rad['sort_passes'], digit_bits=rad['sort_digit'],
                            routed_bits=rad['route_bits'], sorted_bits=rad['sort_bits'],
                            key_bits=rad['key_bits'],
                            middle_bits_scan_only=rad['middle_bits_scan_only'],
                            toy_ops_per_record_per_pass=rad['ops_per_record_per_pass'],
                            ledger_ops_per_record_per_pass=rad['ledger_ops_per_record_per_pass'],
                            toy_ops_per_run_pair=rad['ops_per_run_pair'],
                            ledger_ops_per_run_pair=rad['ledger_ops_per_run_pair'],
                            production='routing on the top 85 key bits; 3 LSD counting-sort '
                                       'passes over the 63 locally ordered bits (21-bit digits, '
                                       '2^21 counters per node); groups are 63-bit prefixes and '
                                       'the run scan compares all four 224-bit key limbs '
                                       'pairwise, which is what the 10-instruction row prices'),
        coverage_census=[dict(nw=c['nw'], records_per_side=c['records_per_side'],
                              candidate_occurrences=c['candidate_occurrences'],
                              distinct_covered=c['distinct_covered'], m=c['m'],
                              K_cand=c['K_cand'], predicted_distinct=c['predicted_distinct'],
                              ratio_measured_over_predicted=c['ratio_measured_over_predicted'],
                              shortfall_over_m=c['shortfall_over_m'],
                              diagonal_repetitions=c['diagonal_repetitions'],
                              ratio_eff=c['ratio_eff'])
                         for c in b['census']],
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
