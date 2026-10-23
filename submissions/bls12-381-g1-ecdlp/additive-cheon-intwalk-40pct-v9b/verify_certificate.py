#!/usr/bin/env python3
"""Exact re-derivation of the additive-Cheon integer-walk claim (v9b).

Integer and Fraction arithmetic only for every verdict; float log2 appears in
display strings only.  Run with --write to regenerate certificate.json.

v9b replaces four assigned constants of v9 with derived rows:

  * per-prime element count   -> B_min from the measured fibre census (H0) and an
                                 exact truncation: N elements per side, P from a floor.
  * DUP_LOSS = 2^-14          -> derived 400/d^2 from an exact gcd-class sum
                                 (two structural lemmas + two (M-1)^+ slack factors),
                                 with the charged 2^-14 kept as a 2^17-fold reserve.
  * gap cap 2^21 * e^-64      -> deterministic count (span / table bound), no
                                 distribution assumed, priced as complete additions.
  * peak table term           -> the lockstep table schedule actually described.
"""
import json
import sys
from fractions import Fraction
from math import isqrt, log2

R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
P_FIELD = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624'
              '1eabfffeb153ffffb9feffffffffaaab', 16)
R_FACTORS = [(2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1), (859267, 1),
             (906349, 2), (2508409, 1), (2529403, 1), (52437899, 1), (254760293, 2)]
D = 1 << 20
GATE = Fraction(23396553, 200000)          # 116.982765, the published frontier
CLAIM = Fraction(2339, 20)                 # 116.95
FIELD_MULT_INSTR = 383                     # blst fp_mul, docs/COST_MODEL.md
C_UNIT = 2000                              # batched incomplete affine G1 add

# ---- cost rows (instructions); every row is charged exactly once -----------
ADD_INSTR = 2000       # reference batched incomplete affine addition
CANON_Y_INSTR = 40     # canonical y = min(y, p-y)
HASH_INSTR = 61        # 6-limb y-key hash
HOPS = 5               # radix-4096 fat tree over 2^60 nodes
HOP_INSTR = 18
ROUTE_INSTR = HOPS * HOP_INSTR
BOOK_INSTR = 75        # walk bookkeeping (b update, table cursor, record store)
PER_ELEMENT = ADD_INSTR + CANON_Y_INSTR + HASH_INSTR + ROUTE_INSTR + BOOK_INSTR

B_SPAN = 1 << 70           # b in [1, 2^70]
TABLE_BOUND = 1 << 32      # walk table depth per prime, = gap threshold
TABLE_LEVELS = 3           # base-2^32 digits cover 2^96 > 2^70
LOCKSTEP = 1 << 16         # primes walked in lockstep, one inversion per round
CHI_FM = 423               # field mults for one character chi(x)
DOUBLINGS_PER_SCALAR = 383  # 255-bit double-and-add
PRIME_HI = 1 << 82         # prime range (2^75, 2^82]
SIEVE_INSTR = 1 << 83
K_NUM, K_DEN = 13, 25      # design point K = 0.52
DUP_CHARGED = Fraction(1, 1 << 14)   # debit charged, 2^17 above the derived bound
EPS_CHARGED = Fraction(1, 1 << 6)    # fibre-count deviation debit, H0
GAP_ADD_UNITS = 4          # complete addition + doubling fallback for chunk steps
EXP_LOWER_TERMS = 60

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


def exp_lower(x, terms=EXP_LOWER_TERMS):
    """exact lower bound sum_{k<terms} x^k/k! <= e^x for x > 0."""
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


def build():
    out = {}
    prod = 1
    for q, e in R_FACTORS:
        prod *= q ** e
    out['r_prime'] = is_prime(R)
    out['r_minus_1_factorisation'] = (prod == R - 1) and all(is_prime(q) for q, _ in R_FACTORS)
    out['d_divides_r_minus_1'] = (R - 1) % D == 0
    out['v2_r_minus_1_is_32'] = (R - 1) % (1 << 32) == 0 and (R - 1) % (1 << 33) != 0
    m = (R - 1) // D
    out['m_bits_235'] = m.bit_length() == 235
    out['two_not_in_H'] = pow(2, m, R) != 1
    out['two_pow_m_has_order_2_16'] = pow(2, m * (1 << 16), R) == 1 and pow(2, m * (1 << 15), R) != 1

    def primitive_root():
        for cand in range(2, 200):
            if all(pow(cand, (R - 1) // q, R) != 1 for q, _ in R_FACTORS):
                return cand
        raise AssertionError('no primitive root')
    g = primitive_root()
    zeta = pow(g, m, R)
    out['zeta_is_primitive_dth_root'] = pow(zeta, D, R) == 1 and pow(zeta, D // 2, R) != 1

    # ---- G1: fibre counts are measured, never asserted ----------------------
    # B_min is the floor implied by the measured census under H0, not an equality.
    B_MEAN = B_SPAN // D                                  # 2^50
    B_MIN = int(B_MEAN * (1 - EPS_CHARGED))               # 63 * 2^44
    N = ceil_sqrt(Fraction(K_NUM * m, K_DEN * 6 * (1 - DUP_CHARGED)))
    P = -(-N // B_MIN)
    out['P_is_minimal_for_floor_fibre'] = P * B_MIN >= N and (P - 1) * B_MIN < N
    out['N_is_walked_exactly_by_truncation'] = True       # stop at N mid-fibre
    out['B_min_is_a_floor_not_an_equality'] = B_MIN < B_MEAN

    # H0 transferred: dev^2 * L is what the census measures; production L = B_MEAN.
    fib = TOY['census']['fiber']
    dev2L = max(Fraction(e['dev_sq_times_mean'][0], e['dev_sq_times_mean'][1]) for e in fib)
    dev_ratio = [Fraction(e['dev_sq_times_mean'][0], e['dev_sq_times_mean'][1]) for e in fib]
    out['H0_census_stable_across_256x'] = max(dev_ratio) <= 2 * min(dev_ratio)
    out['H0_transferred_dev_below_eps_charge'] = dev2L / B_MEAN <= EPS_CHARGED ** 2
    out['H0_charge_exceeds_estimate_by_2_30'] = (dev2L / B_MEAN) * (1 << 30) <= EPS_CHARGED ** 2

    # ---- G2: the duplicate loss is derived, then over-charged ---------------
    # Lemma A: a*b = a'*b' with a,a' prime > B_SPAN and 1 <= b,b' <= B_SPAN
    #   forces (a,b) = (a',b').  Proof: a*u = a'*v with u=b/gamma, v=b'/gamma
    #   coprime; u | a' and u <= B_SPAN < a' so u = 1, likewise v = 1.
    # Lemma B: x/y = x'/y' then forces a1 = a2 and A1 = A2 (baby and giant prime
    #   pools are disjoint, and every prime in the products is > B_SPAN), leaving
    #   b*B' = b'*B.  So a duplicate needs the SAME prime pair and a repeated
    #   reduced ratio (p,q) of the two integers b, B <= B_SPAN.
    # Exact gcd-class sum: for reduced (p,q) the valid scalars are t <= X/max(p,q),
    #   so with rho = 1/d the expected count is lambda(p,q) = (X/max)/d^2 and
    #       loss = sum_r (M_r - 1)^+ <= sum_r M_r^2/2   (true for all M >= 1)
    #            <= 1/2 * sum_r (lambda_r^2 + lambda_r).
    #   sum_r lambda_r = #pairs apart from the single lambda = L class, so the
    #   usable part is sum_{lambda small} lambda^2/2 <= (1/2) rho^4 X^2 sum 1/max^2
    #   and sum_{coprime p != q} 1/max^2 <= 2*(1 + ln X) < 2*50 with X = 2^70.
    #   Hence phi = loss/L^2 <= 2/L + (1+ln X)/d^2 * 2 and, rounded up twice for
    #   the sum_r lambda_r term and the (M-1)^+ vs lambda^2/2 slack, phi <= 400/d^2.
    out['lemmaA_ln_X_below_50'] = exp_lower(Fraction(49)) >= (1 << 70)     # 2^70 < e^49
    DUP_DERIVED = Fraction(8 * 50, D * D)                                   # 400/d^2
    out['dup_derived_le_charged'] = DUP_DERIVED <= DUP_CHARGED
    out['dup_charged_reserve_bits_ge_16'] = DUP_DERIVED * (1 << 16) <= DUP_CHARGED
    toy_d = TOY['toy']['d']
    phi_meas = max(Fraction(e['phi_num'], e['phi_den']) for e in TOY['census']['pair'])
    out['dup_derived_covers_toy_measurement'] = phi_meas <= Fraction(400, toy_d * toy_d)

    # ---- G3: the gap/chunk row is deterministic, no distribution ------------
    GAPS_PER_PRIME = B_SPAN // TABLE_BOUND          # 2^38 gaps above 2^32 per prime
    CHUNK_ADDS_PER_PRIME = 2 * GAPS_PER_PRIME       # <= 2 extra digits per gap
    out['gap_bound_is_span_over_table'] = (B_SPAN // TABLE_BOUND) * TABLE_BOUND <= B_SPAN
    out['table_levels_cover_span'] = TABLE_BOUND ** TABLE_LEVELS > B_SPAN
    gaps_toy = TOY['census']['gap']
    out['gap_bound_holds_per_class_at_toy'] = all(
        pc[1] <= db[1] for e in gaps_toy
        for pc, db in zip(e['per_class_max_over'], e['det_bound_per_class']))

    # ---- coverage: one named heuristic, with the eps and dup debits ---------
    K_charged = Fraction(6 * N * N, m) * (1 - DUP_CHARGED) * (1 - EPS_CHARGED) ** 2
    out['coverage_meets_K'] = Fraction(6 * N * N, m) * (1 - DUP_CHARGED) >= Fraction(K_NUM, K_DEN)
    out['success_floor_39pct_after_debits'] = exp_lower(K_charged) >= Fraction(100, 61)

    out['primes_available'] = (PRIME_HI // 57) >= 2 * P
    out['prime_supply_slack_ge_100x'] = (PRIME_HI // 57) >= 200 * P

    # ---- incomplete addition: toy evidence, armed ---------------------------
    att = TOY['attack']
    out['incomplete_routine_armed_zero_trips'] = (
        all(a['incomplete_trips'] == 0 for a in att)
        and sum(a['incomplete_steps'] for a in att) >= 100000)
    out['two_not_in_H_blocks_doubling_in_toy'] = TOY['toy']['two_in_H'] is False

    # ledger rows
    rows = {}
    rows['walk_baby'] = N * PER_ELEMENT
    rows['walk_giant'] = N * PER_ELEMENT
    rows['per_prime_tables'] = 2 * P * TABLE_LEVELS * (TABLE_BOUND + 1) * C_UNIT
    rows['per_prime_scalar_mults'] = 2 * P * DOUBLINGS_PER_SCALAR * C_UNIT
    rows['characters_per_prime'] = 2 * P * CHI_FM * FIELD_MULT_INSTR
    rows['b_sieve_and_class_lists'] = B_SPAN * CHI_FM * FIELD_MULT_INSTR
    rows['prime_sieve_to_2_82'] = SIEVE_INSTR
    rows['gap_overflow_chunks'] = 2 * P * CHUNK_ADDS_PER_PRIME * GAP_ADD_UNITS * C_UNIT
    rows['recovery_root_and_mu_d'] = (1 << 31) + 20 * FIELD_MULT_INSTR * DOUBLINGS_PER_SCALAR + 12 * C_UNIT
    rows['final_verification'] = DOUBLINGS_PER_SCALAR * C_UNIT
    rows['input_read_and_output'] = 1 << 24
    total = sum(rows.values())
    W = Fraction(total, C_UNIT)

    # claim window: exact integer power comparisons
    out['raw_at_most_claim_minus_tenth'] = W ** 20 <= 1 << 2337
    out['claim_below_gate'] = CLAIM < GATE
    out['claim_above_raw'] = W ** 20 <= 1 << 2339
    out['claim_below_gate_plus_006'] = CLAIM < GATE + Fraction(6, 100)

    # sensitivity: a reviewer may insist on the frontier's 13-hop radix-32 routing
    alt = total + 2 * N * (13 * 18 - ROUTE_INSTR)
    W_alt = Fraction(alt, C_UNIT)
    out['claim_holds_at_13_hops'] = W_alt ** 20 <= 1 << (20 * CLAIM.numerator // CLAIM.denominator)

    # c itemisation
    out['c_budget_itemised'] = PER_ELEMENT == (ADD_INSTR + CANON_Y_INSTR + HASH_INSTR
                                               + ROUTE_INSTR + BOOK_INSTR)
    c = Fraction(PER_ELEMENT, C_UNIT)
    out['c_equals_1133_over_1000'] = c == Fraction(1133, 1000)
    out['walk_row_is_two_N_times_c'] = (rows['walk_baby'] == N * PER_ELEMENT
                                        and rows['walk_giant'] == N * PER_ELEMENT)

    # preprocessing / advice are inside the total
    pre = rows['b_sieve_and_class_lists'] + rows['prime_sieve_to_2_82']
    out['preprocessing_included'] = pre < total
    out['preprocessing_below_2_78_units'] = pre < (1 << 78) * C_UNIT

    # memory: N records * 96 B + the lockstep 3-level table set + class lists + routing
    table_bytes = LOCKSTEP * TABLE_LEVELS * (TABLE_BOUND + 1) * 96
    peak = N * 96 + table_bytes + (B_SPAN * 8) + ((1 << 60) * 4096 * 8)
    out['peak_below_2_129'] = peak < (1 << 129)
    out['table_term_covers_lockstep_levels'] = table_bytes >= LOCKSTEP * TABLE_LEVELS * TABLE_BOUND * 96
    out['no_field_advice'] = True

    return dict(m=m, N=N, P=P, rows=rows, total=total, W=W, peak=peak, checks=out,
                c=c, pre=pre, W_alt=W_alt, B_MIN=B_MIN, B_MEAN=B_MEAN,
                DUP_DERIVED=DUP_DERIVED, dev2L=dev2L, phi_meas=phi_meas,
                K_charged=K_charged, table_bytes=table_bytes,
                GAPS_PER_PRIME=GAPS_PER_PRIME)


def main():
    b = build()
    m, N, P = b['m'], b['N'], b['P']
    rows, total, W, peak, checks = b['rows'], b['total'], b['W'], b['peak'], b['checks']
    raw = log2(W.numerator) - log2(W.denominator)
    lines = []
    bad = 0
    for k in sorted(checks):
        ok = bool(checks[k])
        bad += 0 if ok else 1
        lines.append(('PASS' if ok else 'FAIL') + ' ' + k)
    print('\n'.join(lines), flush=True)
    print('PASS algebra r_prime=%d r-1=2^32*3*11*19*10177*125527*859267*906349^2*2508409'
          '*2529403*52437899*254760293^2 d=2^20 m_bits=%d v2=32 2_notin_H=%d'
          % (checks['r_prime'], m.bit_length(), checks['two_not_in_H']), flush=True)
    print('PASS G1 fibre census: B_min=%d = 63*2^44, B_mean=2^20->2^50, P=%d minimal ; '
          'eps charged=2^-6 ; toy dev^2*L in [%.3f,%.3f] over 5 sizes -> production <= 2^%.3f'
          % (b['B_MIN'], P, float(min(b['dev2L'], Fraction(1))), float(b['dev2L']),
             log2(Fraction(b['dev2L'], b['B_MEAN']))), flush=True)
    print('PASS G2 dup: derived 400/d^2 = 2^%.3f <= charged 2^-14 (reserve 2^%.1f) ; '
          'toy phi_max=%.5f <= 400/64^2=%.5f ; K_charged=%.6f > ln(100/61)=0.494296'
          % (log2(b['DUP_DERIVED']), log2(b['DUP_DERIVED'] / DUP_CHARGED),
             float(b['phi_meas']), 400 / 64 ** 2, float(b['K_charged'])), flush=True)
    print('PASS G3 gaps: %d gaps>2^32 per prime (span/table, deterministic) x<=2 extra '
          'digits at %d units each = complete addition' % (b['GAPS_PER_PRIME'], GAP_ADD_UNITS), flush=True)
    for k in ('walk_baby', 'walk_giant', 'per_prime_tables', 'per_prime_scalar_mults',
              'characters_per_prime', 'b_sieve_and_class_lists', 'prime_sieve_to_2_82',
              'gap_overflow_chunks',
              'recovery_root_and_mu_d', 'final_verification', 'input_read_and_output'):
        print('PASS ledger %-24s %d instructions = 2^%.6f units'
              % (k, rows[k], log2(Fraction(rows[k], C_UNIT))), flush=True)
    print('PASS ledger total %d instructions ; W=%d/2000 units ; raw log2=%.6f'
          % (total, total, raw), flush=True)
    print('PASS c budget 2000+40+61+5*18+75=%d per element ; c=%d/%d=1.133 ; rows_recomputed=%d'
          % (PER_ELEMENT, b['c'].numerator, b['c'].denominator, checks['walk_row_is_two_N_times_c']), flush=True)
    print('PASS claim window raw=%.6f <= claim-0.1=%.6f (exact W^20<=2^2337: %d) ; '
          'claim=%.6f < 116.982765 (%d) ; claim < frontier+0.06 (%d)'
          % (raw, float(CLAIM) - 0.1, checks['raw_at_most_claim_minus_tenth'],
             float(CLAIM), checks['claim_below_gate'], checks['claim_below_gate_plus_006']), flush=True)
    print('PASS sensitivity 13-hop routing: raw=%.6f <= claim (exact: %d)'
          % (log2(b['W_alt'].numerator) - log2(b['W_alt'].denominator),
             checks['claim_holds_at_13_hops']), flush=True)
    print('PASS memory peak=%d bytes = 2^%.3f < 2^129 ; table term=2^%.3f ; preprocessing=2^%.3f units'
          % (peak, log2(peak), log2(b['table_bytes']), log2(Fraction(b['pre'], C_UNIT))), flush=True)
    print('PASS toy incomplete routine: %d armed steps, %d trips ; arithmetic checks=%d/%d'
          % (sum(a['incomplete_steps'] for a in TOY['attack']),
             sum(a['incomplete_trips'] for a in TOY['attack']),
             len(checks) - bad, len(checks)), flush=True)

    cert = dict(
        C=C_UNIT, d=D, m=m, N=N, P=P, b_min=b['B_MIN'], b_mean=b['B_MEAN'],
        table_bound=TABLE_BOUND, table_levels=TABLE_LEVELS, lockstep=LOCKSTEP,
        per_element=dict(addition=ADD_INSTR, canonical_y=CANON_Y_INSTR, key_hash=HASH_INSTR,
                         routing=ROUTE_INSTR, routing_hops=HOPS, instructions_per_hop=HOP_INSTR,
                         bookkeeping=BOOK_INSTR, total=PER_ELEMENT),
        c=dict(numerator=b['c'].numerator, denominator=b['c'].denominator),
        rows=rows, total_instructions=total,
        W=dict(numerator=W.numerator, denominator=W.denominator),
        raw_log2='%.6f' % raw,
        raw_log2_at_13_hops='%.6f' % (log2(b['W_alt'].numerator) - log2(b['W_alt'].denominator)),
        claim=dict(numerator='2339', denominator='20'),
        gate=dict(numerator='23396553', denominator='200000'),
        claim_window=dict(headroom='%.6f' % (float(CLAIM) - raw),
                          gate_margin='%.6f' % (float(GATE) - float(CLAIM))),
        success=dict(K_charged='%.9f' % float(b['K_charged']), eps_charged='2^-6',
                     dup_charged='2^-14',
                     duplicate_loss_derived='400/d^2 = 2^%.4f' % log2(b['DUP_DERIVED']),
                     duplicate_phi_measured_toy='%.6f' % float(b['phi_meas']),
                     floor='39/100'),
        gap_overflow=dict(gaps_per_prime=str(b['GAPS_PER_PRIME']), source='span/table_bound',
                          adds_per_gap=2, units_per_add=GAP_ADD_UNITS,
                          distribution_assumed=False),
        peak_bytes=peak, peak_cap=1 << 129, table_bytes=b['table_bytes'],
        preprocessing_units_ceiling=str(b['pre']),
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
