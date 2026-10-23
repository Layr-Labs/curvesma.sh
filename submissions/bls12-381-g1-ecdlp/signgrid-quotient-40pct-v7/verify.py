#!/usr/bin/env python3
"""Exact integer/Fraction verifier; no floating-point arithmetic or timing prices.
Normal run writes certificate.json and inventory.json after validating all files.
--arithmetic-only is a bootstrap: no complete package certificate is written.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import gcd, isqrt, prod
from pathlib import Path

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


def certificate(extra_per_record=256, claim=F(3331, 25)):
    r, d, g, C = R_ORDER, 2 ** 20, 7, 2000
    q = (r - 1) // d
    Q, K = q // 2, -(-(39 * (q // 2)) // 100)
    N = isqrt(4 * K)
    N += N * N < 4 * K
    a, b = N // 2, N - N // 2
    L = 1 << (N - 1).bit_length()
    zeta = pow(g, d, r)
    E = 510 * (64 * 2048 + 256 * 64 + 4096) + 1024 * 2048 + 2 ** 18
    route = (3 * 104 + 8) * 129 * 512
    rows = dict(grid_scalar_actions=N * E, routed_sorting=N * route,
                metadata_and_initialization=N * 65536,
                sign_record_maintenance=N * extra_per_record,
                dense_root_and_final_verification=21 * E,
                cold_setup=2 ** 80, input_scan=2 ** 40)
    R = sum(rows.values())
    W, stress = F(R, C), F(6 * R, 5 * C)
    raw, stressed = log2_bounds(W), log2_bounds(stress)
    floor = tuple(x + F(1, 10) for x in stressed)
    roots = [dict(t=t, e=d >> (t + 1), coefficient=(d >> (t + 1)) * q * 2 ** t) for t in range(20)]
    endpoints = [[k, k // a, k % a] for k in [0, a - 1, a, a * b - 1]]
    cold = dict(factor_trial_discovery_and_verification=2 ** 49,
                known_modular_arithmetic=2 ** 40, extra_actions=2 ** 14 * E,
                counters_negations_code_output=2 ** 40)
    checks = dict(
        factors_product=prod(ell ** e for ell, e in FACTORS) == r - 1,
        factors_prime=all(prime(ell) for ell, _ in FACTORS),
        lucas_r_prime=pow(g, r - 1, r) == 1 and all(gcd(pow(g, (r - 1) // ell, r) - 1, r) == 1 for ell, _ in FACTORS),
        d_is_supplied=d == 2 ** 20 and d * q == r - 1,
        q_even=q % 2 == 0 and Q * 2 == q,
        q_not_assumed_prime=q % 4096 == 0,
        sign_scalar=pow(zeta, Q, r) == r - 1,
        primitive_zeta=pow(zeta, q, r) == 1 and all(pow(zeta, q // ell, r) != 1 for ell, _ in FACTORS),
        rectangle_meets_target=a * b >= K,
        rectangle_minimal=(N - 1) // 2 * ((N - 1) - (N - 1) // 2) < K,
        rectangle_balanced=a + b == N and 0 <= b - a <= 1,
        no_quotient_wrap=0 < a * b < Q,
        baby_distinct_bound=a - 1 < Q,
        giant_distinct_bound=a * (b - 1) < Q,
        quotient_tiling=all(0 <= i < b and 0 <= j < a and a * i + j == k for k, i, j in endpoints)
                       and (b - 1) * a + (a - 1) + 1 == a * b,
        two_lifts_disjoint=0 <= a * b <= Q and Q + a * b <= q,
        integer_coverage_identity=2 * (a * b) == a * b + ((Q + a * b) - Q),
        exact_success_floor=100 * 2 * a * b >= 39 * q,
        orientation_reconstruction=all(((u - v) * Q) % q == ((u ^ v) * Q) % q for u in [0, 1] for v in [0, 1]),
        dense_exponents_supplied=all(1 <= row['e'] <= d for row in roots),
        dense_binary_identity=all(2 * row['coefficient'] == r - 1 for row in roots),
        dense_sign_nontrivial=pow(g, (r - 1) // 2, r) == r - 1,
        full_scalar_action_price=E == 79650816,
        routing_row=route == 21135360,
        sign_extra_exact=extra_per_record == 256,
        record_layout=48 + 48 + 8 + 16 + 8 == 128,
        canonicalization_ceiling=sum([48, 32, 64, 48, 32, 32]) == extra_per_record,
        counters_and_addresses=N < 2 ** 128 and 2 * L * 128 < 2 ** 192,
        routing_height=L.bit_length() - 1 <= 129,
        capacity=N <= L < 2 * N,
        peak_memory=1536 * L + 2 ** 64 < 2 ** 129,
        metadata_bound=8192 * L + 8192 * N < 65536 * N,
        cold_bound=sum(cold.values()) < 2 ** 50 < 2 ** 80,
        ledger_sum=R == N * (E + route + 65536 + extra_per_record) + 21 * E + 2 ** 80 + 2 ** 40,
        source_truth_operating_point=a == 98748750602204457565883422358598783 and b == 98748750602204457565883422358598784,
        source_truth_R=R == 19918011671547009368992591426495508435762944,
        claim_reduced=gcd(claim.numerator, claim.denominator) == 1 and claim == F(3331, 25),
        base_claim=under_power(W, claim),
        all_row_reprice_headroom=under_power(stress, claim - F(1, 10)),
        hundredth_133_23_fails=not under_power(stress, F(13323, 100) - F(1, 10)),
        claim_window_min=claim - floor[1] >= F(5, 1000),
        claim_window_max=claim - floor[0] <= F(6, 100),
        log_interval_width=all(hi - lo < F(1, 10 ** 25) for lo, hi in [raw, stressed]),
        n_range_input_memory=2 ** 32 * 512 < 2 ** 64,
    )
    assert all(checks.values()), [k for k, v in checks.items() if not v]
    return dict(checks=checks, claim=rational(claim), r=r, g=g, d=d, q=q, Q=Q,
                K=K, baby=a, giant=b, grid_records=N, C=C, H=0, E=E, route=route,
                sign_extra_per_record=extra_per_record, rows=rows, R=R, W=rational(W),
                stress_work=rational(stress), success=rational(F(2 * a * b, q)),
                exact_covered=2 * a * b, exact_coverage_identity='[0,ab) disjoint-union [Q,Q+ab); cardinality=2ab',
                raw_log2_interval=display(raw), stressed_log2_interval=display(stressed),
                policy_floor_interval=display(floor), claim_slack_interval=display((claim - floor[1], claim - floor[0])),
                repriced_headroom_interval=display((claim - stressed[1], claim - stressed[0])),
                endpoints=endpoints, factors=[list(x) for x in FACTORS], dense_root_identities=roots,
                cold_subceilings=cold, L=L, peak_bytes=1536 * L + 2 ** 64,
                input_n_range=[2 ** 20, 2 ** 32], record_bytes=128, key_bytes=104)


def package_checks(cert):
    sub = json.loads((HERE / 'submission.json').read_text())
    note = (HERE / 'public-note.md').read_bytes()
    argument = (HERE / 'argument.md').read_bytes()
    rows = sum([json.loads((HERE / name).read_text()) for name in ['toy-grid-a.json', 'toy-grid-b.json']], [])
    summary = json.loads((HERE / 'toy-summary.json').read_text())
    expected = 0
    for row in rows:
        a, b, q = row['a'], row['b'], row['q']
        keys = {a * i + j + t * (q // 2) for i in range(b) for j in range(a) for t in range(2)}
        assert len(keys) == 2 * a * b
        assert row['random_trials'] >= 2000
        assert row['exhaustive_k2_trials'] == q and row['exhaustive_k2_hits'] == len(keys)
        assert row['wrong_outputs'] == row['coverage_prediction_errors'] == row['duplicate_keys'] == 0
        assert sum(row['k2_histogram']) == row['random_trials']
        assert row['random_hits'] == sum(row['k2_histogram'][k] for k in keys) == sum(row['hit_k2_histogram'])
        assert all(row['hit_k2_histogram'][k] == (row['k2_histogram'][k] if k in keys else 0) for k in range(q))
        assert sum(row['orientation_pair_histogram']) == row['random_hits']
        assert all(sum(x) == row['random_hits'] for x in row['dense_bit_histogram'])
        assert sum(row['scalar_action_histogram'].values()) == row['random_trials']
        expected += row['random_trials']
    checks = dict(manifest_schema=set(sub) == {'version', 'title', 'claims', 'experiments'} and sub['version'] == 1,
                  manifest_claim=sub['claims'][0]['score'] == cert['claim'] and set(sub['claims'][0]) == {'statement', 'score', 'evidence'},
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
                      and summary['exhaustive_all_secrets']['hits'] == 5376
                      and summary['exhaustive_all_secrets']['wrong_outputs'] == 0,
                  evidence_files=all(name == 'certificate.json' or (HERE / name).is_file()
                                     for name in sub['claims'][0]['evidence']),
                  bounded_files=all(p.stat().st_size <= 16384 for p in HERE.rglob('*') if p.is_file()),
                  no_bytecode=not any(HERE.rglob('*.pyc')),
                  file_count=len([p for p in HERE.rglob('*') if p.is_file()]) <= 32)
    assert all(checks.values()), [k for k, v in checks.items() if not v]
    cert['checks'].update(checks)
    cert['toy_summary_sha256'] = sha256((HERE / 'toy-summary.json').read_bytes()).hexdigest()
    cert['toy_configurations'] = len(rows)
    cert['random_trials'] = expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--arithmetic-only', action='store_true')
    args = parser.parse_args()
    cert = certificate()
    if args.arithmetic_only:
        print('PASS arithmetic checks=' + str(len(cert['checks'])) + '; DOCUMENTATION PENDING; no package certificate written')
        return
    package_checks(cert)
    (HERE / 'certificate.json').write_text(json.dumps(cert, indent=2, sort_keys=True) + '\n')
    inventory = {str(p.relative_to(HERE)): sha256(p.read_bytes()).hexdigest()
                 for p in sorted(HERE.rglob('*')) if p.is_file() and p.name != 'inventory.json'}
    assert sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file()) < 262144
    (HERE / 'inventory.json').write_text(json.dumps(inventory, indent=2, sort_keys=True) + '\n')
    print('PASS exact certificate checks=' + str(len(cert['checks'])) + ' claim=3331/25 C=2000; coverage=2ab/q')
    print('PASS raw=' + str(cert['raw_log2_interval']) + ' policy_floor=' + str(cert['policy_floor_interval']) + '; 133.23 FAILS same policy')
    print('PASS inventory files=' + str(len(inventory)) + ' note_bytes=' + str((HERE / 'public-note.md').stat().st_size))


if __name__ == '__main__':
    main()
