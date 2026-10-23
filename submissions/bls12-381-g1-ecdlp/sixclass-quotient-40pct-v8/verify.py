#!/usr/bin/env python3
"""Exact integer/Fraction verifier; no floating-point arithmetic or timing prices.
Normal run writes certificate.json and inventory.json after validating all files.
--arithmetic-only is a bootstrap: no complete package certificate is written.
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import gcd, isqrt, prod
from pathlib import Path
from toy import Curve, DIGIT_SHIFTS, DIGIT_TABLE

HERE = Path(__file__).resolve().parent
R_ORDER = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
FACTORS = [(2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1),
           (859267, 1), (906349, 2), (2508409, 1), (2529403, 1),
           (52437899, 1), (254760293, 2)]


def prime(n):
    return n == 2 or (n > 2 and n % 2 != 0 and all(n % j for j in range(3, isqrt(n) + 1, 2)))


def under_power(work, exponent):
    work, exponent = F(work), F(exponent)
    a, b = work.numerator ** exponent.denominator, work.denominator ** exponent.denominator
    if exponent.numerator >= 0:
        b <<= exponent.numerator
    else:
        a <<= -exponent.numerator
    return a < b


def ln_unit(x):
    assert 1 <= x <= 2
    u = (x - 1) / (x + 1)
    v, term, acc, count = u * u, u, F(0), 36
    for j in range(count):
        acc += term / (2 * j + 1)
        term *= v
    lower = 2 * acc
    return lower, lower + 2 * term / ((2 * count + 1) * (1 - v))


@lru_cache(maxsize=None)
def log2_bounds(work):
    work = F(work)
    assert work > 0
    k = work.numerator.bit_length() - work.denominator.bit_length()
    unit = work / (2 ** k if k >= 0 else F(1, 2 ** (-k)))
    if unit < 1:
        unit *= 2
        k -= 1
    assert 1 <= unit < 2
    low, high = ln_unit(unit)
    lo2, hi2 = ln_unit(F(2))
    lower, upper = k + low / hi2, k + high / lo2
    assert upper - lower < F(1, 10 ** 25)
    return lower, upper


def fixed(n, places=12):
    sign, n = ('-' if n < 0 else ''), abs(n)
    return f'{sign}{n // 10 ** places}.{n % 10 ** places:0{places}d}'


def display(bounds):
    lo, hi = bounds
    scale = 10 ** 12
    return [fixed(lo.numerator * scale // lo.denominator),
            fixed(-(-hi.numerator * scale // hi.denominator))]


def rational(x):
    x = F(x)
    return dict(numerator=str(x.numerator), denominator=str(x.denominator))


def certificate(extra_per_record=8192, claim=F(2649, 20)):
    r, d, g, C = R_ORDER, 2 ** 20, 7, 2000
    p = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
    P = (int('17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb', 16),
         int('08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1', 16))
    beta = int('1a0111ea397fe699ec02408663d4de85aa0d857d89759ad4897d29650fb85f9b409427eb4f49fffd8bfd00000000aaac', 16)
    q = (r - 1) // d
    Q, K = q // 6, -(-(39 * (q // 6)) // 100)
    N = isqrt(4 * K)
    N += N * N < 4 * K
    a, b = N // 2, N - N // 2
    L = 1 << (N - 1).bit_length()
    zeta, lam = pow(g, d, r), pow(g, (r - 1) // 3, r)
    curve = Curve(p, r, P, beta, lam)
    multipliers = [pow(zeta, n * Q, r) for n in DIGIT_SHIFTS]
    E = 510 * (64 * 2048 + 256 * 64 + 4096) + 1024 * 2048 + 2 ** 18
    route = (3 * 104 + 8) * 129 * 512
    canonical_rows = dict(two_fixed_constant_Fp_products=2 * 2048, one_y_negation=64,
                          five_full_lex_comparisons=5 * 256, five_selections=5 * 128,
                          loads_scratch=512, digit_side_packing=128, reserve=1472)
    rows = dict(grid_scalar_actions=N * E, routed_sorting=N * route,
                metadata_and_initialization=N * 65536,
                orbit_record_maintenance=N * extra_per_record,
                dense_root_and_final_verification=21 * E,
                cold_setup=2 ** 80, input_scan=2 ** 40)
    R = sum(rows.values())
    W, stress = F(R, C), F(6 * R, 5 * C)
    raw, stressed = log2_bounds(W), log2_bounds(stress)
    floor = tuple(x + F(1, 10) for x in stressed)
    roots = [dict(t=t, e=d >> (t + 1), coefficient=(d >> (t + 1)) * q * 2 ** t) for t in range(20)]
    endpoints = [[k, k // a, k % a] for k in [0, a - 1, a, a * b - 1]]
    lifts = [[t * Q, t * Q + a * b] for t in range(6)]
    cold = dict(factor_trial_discovery_and_verification=2 ** 49,
                known_modular_arithmetic=2 ** 40, extra_actions=2 ** 14 * E,
                counters_negations_code_output=2 ** 40)
    checks = dict(
        factors_product=prod(ell ** e for ell, e in FACTORS) == r - 1,
        factors_prime=all(prime(ell) for ell, _ in FACTORS),
        lucas_r_prime=pow(g, r - 1, r) == 1 and all(gcd(pow(g, (r - 1) // ell, r) - 1, r) == 1 for ell, _ in FACTORS),
        d_is_supplied=d == 2 ** 20 and d * q == r - 1,
        q_divisible_by_six=q % 6 == 0 and Q * 6 == q,
        q_composite_not_EC_order=q % 4096 == 0 and q != r,
        primitive_zeta=pow(zeta, q, r) == 1 and all(pow(zeta, q // ell, r) != 1 for ell, _ in FACTORS),
        beta_order_three=0 < beta < p and beta != 1 and pow(beta, 3, p) == 1,
        lambda_order_three=lam != 1 and pow(lam, 3, r) == 1,
        lambda_zeta_alignment=pow(zeta, 2 * Q, r) == lam,
        negative_scalar=pow(zeta, 3 * Q, r) == r - 1,
        generator_on_curve=(P[1] ** 2 - P[0] ** 3 - 4) % p == 0,
        generator_order=curve.mul(r, P) is None and P is not None,
        production_endomorphism_alignment=curve.mul(lam, P) == curve.tau(P),
        all_six_multipliers_distinct=len(set(multipliers)) == 6 and set(multipliers) == {pow(zeta, t * Q, r) for t in range(6)},
        no_nonidentity_stabilizers=all(gcd(m - 1, r) == 1 for m in multipliers[1:]),
        digit_enumeration=multipliers == [1, lam, lam * lam % r, r - 1, -lam % r, -(lam * lam) % r],
        digit_table_36=all(pow(zeta, DIGIT_TABLE[u][v] * Q, r) == multipliers[u] * pow(multipliers[v], -1, r) % r for u in range(6) for v in range(6)),
        digit_correct_ring=all(pow(m, 6, r) == 1 for m in multipliers),
        rectangle_meets_target=a * b >= K,
        rectangle_minimal=(N - 1) // 2 * ((N - 1) - (N - 1) // 2) < K,
        rectangle_balanced=a + b == N and 0 <= b - a <= 1,
        no_quotient_wrap=0 < a * b < Q,
        baby_distinct_bound=a - 1 < Q,
        giant_distinct_bound=a * (b - 1) < Q,
        quotient_tiling=all(0 <= i < b and 0 <= j < a and a * i + j == k for k, i, j in endpoints)
                       and (b - 1) * a + (a - 1) + 1 == a * b,
        six_lifts_disjoint=all(0 <= lo < hi <= q for lo, hi in lifts) and all(lifts[t][1] < lifts[t + 1][0] for t in range(5)),
        integer_coverage_identity=6 * (a * b) == sum(hi - lo for lo, hi in lifts),
        exact_success_floor=100 * 6 * a * b >= 39 * q,
        dense_exponents_supplied=all(1 <= row['e'] <= d for row in roots),
        dense_binary_identity=all(2 * row['coefficient'] == r - 1 for row in roots),
        dense_sign_nontrivial=pow(g, (r - 1) // 2, r) == r - 1,
        full_scalar_action_price=E == 79650816,
        routing_row=route == 21135360,
        orbit_extra_exact=extra_per_record == 8192,
        record_layout=48 + 48 + 8 + 16 + 8 == 128,
        canonicalization_ceiling=sum(canonical_rows.values()) == extra_per_record,
        three_bit_digit_packed=max(range(6)) < 2 ** 3 and (1 | (5 << 1)) < 2 ** 64,
        counters_and_addresses=N < 2 ** 128 and 2 * L * 128 < 2 ** 192,
        routing_height=L.bit_length() - 1 <= 129,
        capacity=N <= L < 2 * N,
        peak_memory=1536 * L + 2 ** 64 < 2 ** 129,
        metadata_bound=8192 * L + 8192 * N < 65536 * N,
        cold_bound=sum(cold.values()) < 2 ** 50 < 2 ** 80,
        ledger_sum=R == N * (E + route + 65536 + extra_per_record) + 21 * E + 2 ** 80 + 2 ** 40,
        source_truth_operating_point=a == b == 57012617742321962895928165416718183,
        source_truth_R=R == 11500574304558579830758679329268362223869952,
        claim_reduced=gcd(claim.numerator, claim.denominator) == 1 and claim == F(2649, 20),
        base_claim=under_power(W, claim),
        all_row_reprice_headroom=under_power(stress, claim - F(1, 10)),
        hundredth_132_44_fails=not under_power(stress, F(3311, 25) - F(1, 10)),
        claim_window_min=claim - floor[1] >= F(5, 1000),
        claim_window_max=claim - floor[0] <= F(6, 100),
        log_interval_width=all(hi - lo < F(1, 10 ** 25) for lo, hi in [raw, stressed]),
        n_range_input_memory=2 ** 32 * 512 < 2 ** 64,
    )
    assert all(checks.values()), [k for k, v in checks.items() if not v]
    return dict(checks=checks, claim=rational(claim), r=r, p=p, P=list(P), beta=beta, lam=lam, g=g, d=d, q=q, Q=Q,
                K=K, baby=a, giant=b, grid_records=N, C=C, H=0, E=E, route=route,
                orbit_extra_per_record=extra_per_record, canonical_rows=canonical_rows,
                digit_shifts=DIGIT_SHIFTS, digit_table=DIGIT_TABLE, multipliers=multipliers,
                rows=rows, R=R, W=rational(W), stress_work=rational(stress), success=rational(F(6 * a * b, q)),
                exact_covered=6 * a * b, coverage_intervals=lifts,
                exact_coverage_identity='disjoint union t=0..5 [tQ,tQ+ab); cardinality=6ab',
                raw_log2_interval=display(raw), stressed_log2_interval=display(stressed),
                policy_floor_interval=display(floor), claim_slack_interval=display((claim - floor[1], claim - floor[0])),
                repriced_headroom_interval=display((claim - stressed[1], claim - stressed[0])),
                endpoints=endpoints, factors=[list(x) for x in FACTORS], dense_root_identities=roots,
                cold_subceilings=cold, L=L, peak_bytes=1536 * L + 2 ** 64,
                input_n_range=[2 ** 20, 2 ** 32], record_bytes=128, key_bytes=104)


def package_checks(cert):
    from orbit_audit import audit_orbits
    sub = json.loads((HERE / 'submission.json').read_text())
    note = (HERE / 'public-note.md').read_bytes()
    argument = (HERE / 'argument.md').read_bytes()
    rows = sum([json.loads((HERE / name).read_text()) for name in ['toy-grid-a.json', 'toy-grid-b.json']], [])
    summary = json.loads((HERE / 'toy-summary.json').read_text())
    curve, inputs = Curve(), [None]
    for _ in range(1, curve.r):
        inputs.append(curve.add(inputs[-1], curve.P))
    audit = audit_orbits(curve, inputs, DIGIT_SHIFTS, DIGIT_TABLE)
    squares = {y * y % curve.p for y in range(1, (curve.p + 1) // 2)}
    actual_order = 1 + sum(1 if (v := (x * x * x + 4) % curve.p) == 0 else 2 if v in squares else 0 for x in range(curve.p))
    expected = 0
    for row in rows:
        a, b, q = row['a'], row['b'], row['q']
        Q = q // 6
        keys = {a * i + j + t * Q for i in range(b) for j in range(a) for t in range(6)}
        assert 6 * Q == q and 0 < a * b < Q and len(keys) == 6 * a * b
        assert row['random_trials'] >= 2000
        assert row['exhaustive_k2_trials'] == q and row['exhaustive_k2_hits'] == len(keys)
        assert row['wrong_outputs'] == row['coverage_prediction_errors'] == row['duplicate_keys'] == 0
        assert row['exact_coverage'] == [6 * a * b, q] and row['integer_tiling'] == 'PASS'
        assert sum(row['k2_histogram']) == row['random_trials']
        assert row['random_hits'] == sum(row['k2_histogram'][k] for k in keys) == sum(row['hit_k2_histogram'])
        assert all(row['hit_k2_histogram'][k] == (row['k2_histogram'][k] if k in keys else 0) for k in range(q))
        assert len(row['digit_pair_histogram']) == 36 and sum(row['digit_pair_histogram']) == row['random_hits']
        assert all(sum(x) == row['random_hits'] for x in row['dense_bit_histogram'])
        actions = {str(a + b): row['random_misses'], str(a + b + row['d'].bit_length()): row['random_hits']}
        assert row['scalar_action_histogram'] == actions
        expected += row['random_trials']
    statement = sub['claims'][0]['statement']
    checks = dict(manifest_schema=set(sub) == {'version', 'title', 'claims', 'experiments'} and sub['version'] == 1,
                  manifest_claim=sub['claims'][0]['score'] == cert['claim'] and set(sub['claims'][0]) == {'statement', 'score', 'evidence'},
                  statement_self_contained=all(s not in statement.lower() for s in ['improves', '-v7', '037d0cbc', 'frontier', 'accepted']),
                  no_declared_experiments=sub['experiments'] == [],
                  dependencies_empty=json.loads((HERE / 'dependencies.json').read_text()) == [],
                  note_size=10240 <= len(note) < 14336,
                  argument_size=len(argument) < 8192,
                  note_public=b'/Users/' not in note and b'ykn_' not in note,
                  note_lineage_attribution=b'qwen3.8-flash-next / Hermes' in note and b'gpt-6.1-sol / Hermes' in note,
                  evidence_matrix=len(rows) >= 7 and expected == summary['random_trials'],
                  evidence_summaries=summary['configuration_count'] == len(rows)
                      and summary['random_hits'] == sum(row['random_hits'] for row in rows)
                      and summary['exhaustive_k2_trials'] == sum(row['q'] for row in rows),
                  exhaustive_all_secrets=summary['exhaustive_all_secrets']['trials'] == 12288
                      and summary['exhaustive_all_secrets']['hits'] == 6144
                      and summary['exhaustive_all_secrets']['wrong_outputs'] == 0,
                  toy_p_r_prime=prime(curve.p) and prime(curve.r),
                  toy_curve_actual_order=actual_order == 258069 == 21 * curve.r,
                  orbit_exhaustive_replayed=audit == summary['orbit_audit'],
                  orbit_no_exceptions=audit['orbit_size_exceptions'] == sum(audit['fixed_nonidentity_counts']) == 0,
                  digit_36_exhaustive=audit['digit_pair_tests'] == 73728 and audit['wrong_digit_relations'] == 0,
                  evidence_capture_cap=sum((HERE / name).stat().st_size for name in ['toy-grid-a.json', 'toy-grid-b.json', 'toy-summary.json']) <= 16384,
                  evidence_files=all(name == 'certificate.json' or (HERE / name).is_file() for name in sub['claims'][0]['evidence']),
                  bounded_files=all(p.stat().st_size <= 16384 for p in HERE.rglob('*') if p.is_file()),
                  no_bytecode=not any(HERE.rglob('*.pyc')),
                  file_count=len([p for p in HERE.rglob('*') if p.is_file()]) <= 32)
    assert all(checks.values()), [k for k, v in checks.items() if not v]
    cert['checks'].update(checks)
    cert['toy_summary_sha256'] = sha256((HERE / 'toy-summary.json').read_bytes()).hexdigest()
    cert['embedded_source_sha256'] = {name: sha256((HERE / name).read_bytes()).hexdigest() for name in ['verify.py', 'toy.py', 'orbit_audit.py']}
    cert['toy_configurations'], cert['random_trials'] = len(rows), expected
    cert['orbit_audit'] = audit


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--arithmetic-only', action='store_true')
    args = parser.parse_args()
    cert = certificate()
    if args.arithmetic_only:
        print('PASS arithmetic checks=' + str(len(cert['checks'])) + '; DOCUMENTATION PENDING; no package certificate written', flush=True)
        return
    package_checks(cert)
    (HERE / 'certificate.json').write_text(json.dumps(cert, indent=2, sort_keys=True) + '\n')
    inventory = {str(p.relative_to(HERE)): sha256(p.read_bytes()).hexdigest()
                 for p in sorted(HERE.rglob('*')) if p.is_file() and p.name != 'inventory.json'}
    assert sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file()) < 262144
    (HERE / 'inventory.json').write_text(json.dumps(inventory, indent=2, sort_keys=True) + '\n')
    print('PASS exact certificate checks=' + str(len(cert['checks'])) + ' claim=2649/20 C=2000; coverage=6ab/q', flush=True)
    print('PASS raw=' + str(cert['raw_log2_interval']) + ' policy_floor=' + str(cert['policy_floor_interval']) + '; 132.44 FAILS same policy', flush=True)
    print('PASS inventory files=' + str(len(inventory)) + ' note_bytes=' + str((HERE / 'public-note.md').stat().st_size), flush=True)


if __name__ == '__main__':
    main()
