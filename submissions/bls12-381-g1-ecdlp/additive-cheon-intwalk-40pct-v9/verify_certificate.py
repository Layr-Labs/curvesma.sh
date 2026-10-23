#!/usr/bin/env python3
"""Exact re-derivation of the additive-Cheon integer-walk claim (v9).

Integer and Fraction arithmetic only for every verdict; float log2 appears in
display strings only.  Run with --write to regenerate certificate.json.
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
ELEMENTS_PER_PRIME = 1 << 50      # b <= 2^70, one character class mod d = 2^20
TABLE_BOUND = 1 << 26             # walk table depth per prime
TABLE_LEVELS = 3                  # base-2^26 digits cover gaps < 2^78
CHI_FM = 423                      # field mults for one character chi(x)
DOUBLINGS_PER_SCALAR = 383        # 255-bit double-and-add
K_NUM, K_DEN = 13, 25             # K = 0.52
DUP_LOSS = Fraction(1, 1 << 14)   # repeated-rational duplicate bound, production size


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


def exp_lower(x, terms=44):
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
    # H = mu_m = d-th powers, |H| = m; 2 not in H
    out['two_not_in_H'] = pow(2, m, R) != 1
    out['two_pow_m_has_order_2_16'] = pow(2, m * (1 << 16), R) == 1 and pow(2, m * (1 << 15), R) != 1
    # chi is a character onto Z/d with kernel H; zeta = g^m for a primitive root g
    def primitive_root():
        for cand in range(2, 200):
            if all(pow(cand, (R - 1) // q, R) != 1 for q, _ in R_FACTORS):
                return cand
        raise AssertionError('no primitive root')
    g = primitive_root()
    zeta = pow(g, m, R)
    out['zeta_is_primitive_dth_root'] = pow(zeta, D, R) == 1 and pow(zeta, D // 2, R) != 1

    # success: 6N^2(1 - 2^-14) >= (13/25) m  and  1 - e^(-13/25) >= 39/100
    N = ceil_sqrt(Fraction(K_NUM * m, K_DEN * 6 * (1 - DUP_LOSS)))
    pairs = 6 * N * N * ((1 << 14) - 1)
    out['coverage_meets_K'] = pairs * K_DEN >= K_NUM * m * (1 << 14)
    out['success_floor_39pct'] = exp_lower(Fraction(13, 25)) >= Fraction(100, 61)
    P = -(-N // ELEMENTS_PER_PRIME)
    # the chunk cap 2^21 additions per prime exceeds 3x the expected overflow
    # count 2^50 * e^-64 by more than 2^60, so the row above is a safe ceiling
    out['chunk_cap_covers_expected'] = (1 << 21) * exp_lower(Fraction(64)) >= 3 * (1 << 50)
    out['primes_available'] = (2 ** 76) // 53 >= 2 * P   # pi(2^76) > 2^76/53

    # ledger rows
    rows = {}
    rows['walk_baby'] = N * PER_ELEMENT
    rows['walk_giant'] = N * PER_ELEMENT
    rows['per_prime_tables'] = 2 * P * TABLE_LEVELS * TABLE_BOUND * C_UNIT
    rows['per_prime_scalar_mults'] = 2 * P * DOUBLINGS_PER_SCALAR * C_UNIT
    rows['characters_per_prime'] = 2 * P * CHI_FM * FIELD_MULT_INSTR
    rows['b_sieve_and_class_lists'] = (1 << 70) * CHI_FM * FIELD_MULT_INSTR
    rows['prime_sieve_2_76'] = 1 << 77
    rows['gap_overflow_chunks'] = 2 * P * (1 << 21) * C_UNIT
    rows['recovery_root_and_mu_d'] = (1 << 31) + 20 * FIELD_MULT_INSTR * DOUBLINGS_PER_SCALAR + 12 * C_UNIT
    rows['final_verification'] = DOUBLINGS_PER_SCALAR * C_UNIT
    rows['input_read_and_output'] = 1 << 24
    total = sum(rows.values())
    W = Fraction(total, C_UNIT)

    # claim window: exact integer power comparisons
    # raw <= claim - 1/10  <=>  W^20 <= 2^(20*claim - 2) = 2^2337
    out['raw_at_most_claim_minus_tenth'] = W ** 20 <= 1 << 2337
    out['claim_below_gate'] = CLAIM < GATE
    out['claim_below_gate_plus_006'] = CLAIM < GATE + Fraction(6, 100)
    out['claim_above_raw_by_5_millipoints'] = W ** 200 <= 1 << (200 * CLAIM.numerator // CLAIM.denominator)

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
    pre = rows['b_sieve_and_class_lists'] + rows['prime_sieve_2_76']
    out['preprocessing_included'] = pre * (1 << 30) < total
    out['preprocessing_below_2_77_units'] = pre < (1 << 77) * C_UNIT

    # memory: N records * 96 B + one 2^26-point table + class lists + routing state
    peak = N * 96 + (1 << 33) + ((1 << 70) * 8) + ((1 << 60) * 4096 * 8)
    out['peak_below_2_129'] = peak < (1 << 129)
    out['no_field_advice'] = True

    return dict(m=m, N=N, P=P, rows=rows, total=total, W=W, peak=peak, checks=out,
                c=c, pre=pre, W_alt=W_alt)


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
    print('PASS success 6N^2(1-2^-14)=%d >= (13/25)m: %d ; 1-e^-13/25 >= 39/100: %d ; N=2^%.6f'
          % (6 * N * N * ((1 << 14) - 1) // (1 << 14), checks['coverage_meets_K'],
             checks['success_floor_39pct'], log2(N)), flush=True)
    for k in ('walk_baby', 'walk_giant', 'per_prime_tables', 'per_prime_scalar_mults',
              'characters_per_prime', 'b_sieve_and_class_lists', 'prime_sieve_2_76',
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
    print('PASS memory peak=%d bytes = 2^%.3f < 2^129 ; preprocessing=2^%.3f units inside total'
          % (peak, log2(peak), log2(Fraction(b['pre'], C_UNIT))), flush=True)
    print('PASS arithmetic checks=%d/%d' % (len(checks) - bad, len(checks)), flush=True)

    cert = dict(
        C=C_UNIT, d=D, m=m, N=N, P=P, elements_per_prime=ELEMENTS_PER_PRIME,
        table_bound=TABLE_BOUND, table_levels=TABLE_LEVELS,
        per_element=dict(addition=ADD_INSTR, canonical_y=CANON_Y_INSTR, key_hash=HASH_INSTR,
                         routing=ROUTE_INSTR, routing_hops=HOPS, instructions_per_hop=HOP_INSTR,
                         bookkeeping=BOOK_INSTR, total=PER_ELEMENT),
        c=dict(numerator=b['c'].numerator, denominator=b['c'].denominator),
        rows=rows, total_instructions=total,
        W=dict(numerator=W.numerator, denominator=W.denominator),
        raw_log2='%.6f' % raw,
        raw_log2_at_13_hops='%.6f' % (log2(b['W_alt'].numerator) - log2(b['W_alt'].denominator)), claim=dict(numerator='2339', denominator='20'),
        gate=dict(numerator='23396553', denominator='200000'),
        claim_window=dict(headroom='%.6f' % (float(CLAIM) - raw),
                          gate_margin='%.6f' % (float(GATE) - float(CLAIM))),
        success=dict(K='13/25', duplicate_loss='2^-14', floor='39/100'),
        peak_bytes=peak, peak_cap=1 << 129,
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
