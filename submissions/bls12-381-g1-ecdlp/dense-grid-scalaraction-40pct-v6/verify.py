#!/usr/bin/env python3
"""Exact v6 certificate. All asserted bounds use integers, never runtime prices."""
import json
from math import isqrt, prod

R_ORDER = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
FACTORS = [(2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1),
           (859267, 1), (906349, 2), (2508409, 1), (2529403, 1),
           (52437899, 1), (254760293, 2)]


def prime(n):
    return n >= 2 and all(n % j for j in range(2, isqrt(n) + 1))


def ceil_sqrt(n):
    z = isqrt(n)
    return z + (z * z < n)


def certificate():
    r, d, C, g = R_ORDER, 2 ** 20, 2000, 7
    q = (r - 1) // d
    m = ceil_sqrt((81 * q + 199) // 200)
    N = 2 * m
    L = 1 << (N - 1).bit_length()
    # This is the passing v4 scalar-action row, unchanged.
    E = 510 * (64 * 2048 + 256 * 64 + 4096) + 1024 * 2048 + 2 ** 18
    route = (3 * 104 + 8) * 129 * 512
    rows = dict(grid_scalar_actions=N * E,
                routing_sorting=N * route,
                metadata_and_initialization=N * 65536,
                dense_root_and_verification_actions=21 * E,
                cold_setup=2 ** 80,
                input_scan=2 ** 40)
    R = sum(rows.values())
    peak = 1536 * L + 2 ** 64
    # Endpoints and quotient/remainder certify the m-by-m interval identity.
    endpoints = [0, m - 1, m, m * m - 1]
    endpoint_decomposition = [[k, k // m, k % m] for k in endpoints]
    dense_identities = []
    for t in range(20):
        e = d >> (t + 1)
        dense_identities.append(dict(t=t, e=e, coefficient=q * (1 << t) * e))
    # Cold-setup sub-ceilings come from ledger.md, not measured Python time.
    cold_subceilings = dict(factor_prime_testing=2 ** 24,
                           bounded_modular_arithmetic=2 ** 40,
                           extra_cold_scalar_actions=2 ** 14 * E,
                           counters_code_negation_output=2 ** 40)
    checks = dict(
        factorization=prod(l ** k for l, k in FACTORS) == r - 1,
        factors_prime=all(prime(l) for l, _ in FACTORS),
        primitive_root=pow(g, r - 1, r) == 1 and all(pow(g, (r - 1) // l, r) != 1 for l, _ in FACTORS),
        divisor=d * q == r - 1,
        minimal_grid=(m - 1) ** 2 * 200 < 81 * q <= m * m * 200,
        partial_grid_no_wrap=m * m < q,
        exact_coverage_floor=200 * m * m >= 81 * q,
        interval_tiling=all(0 <= i < m and 0 <= j < m and i * m + j == k
                            for k, i, j in endpoint_decomposition)
                        and (m - 1) * m + (m - 1) == m * m - 1,
        dense_supplied_exponents=all(1 <= x['e'] <= d for x in dense_identities),
        dense_binary_identity=all(x['coefficient'] * 2 == r - 1 for x in dense_identities),
        dense_sign_nontrivial=pow(g, (r - 1) // 2, r) == r - 1,
        scalar_row_unchanged=E == 79650816,
        counters_three_limbs=N < 2 ** 128 and L <= 2 ** 128 and 2 * L < 2 ** 192,
        scalar_limbs=r < 2 ** 255 and (r - 1) ** 2 < 2 ** 510,
        routing_height=L.bit_length() - 1 <= 129,
        fixed_capacity=N <= L < 2 * N,
        memory=peak < 2 ** 129,
        metadata_ceiling=8192 * L + 8192 * N < 65536 * N,
        cold_analytic_ceiling=sum(cold_subceilings.values()) < 2 ** 44 < 2 ** 80,
        ledger_sum=R == sum(rows.values()),
        claim=R ** 5 < C ** 5 * 2 ** 669,
        repriced_point_one_bit_headroom=(6 * R) ** 10 < (5 * C) ** 10 * 2 ** 1337,
        rank_strike=669 < 134 * 5,
    )
    assert all(checks.values()), checks
    return dict(checks=checks, claim=[669, 5], C=C, H=0, R=R, rows=rows,
                r=r, g=g, d=d, q=q, m=m, grid_records=N, E=E,
                route=route, route_transfers_per_record=320,
                key_bytes=104, record_bytes=128, L=L,
                peak_bytes=peak, hard_cap_bytes=2 ** 129,
                input_n_range=[2 ** 20, 2 ** 32],
                success=[m * m, q], success_floor=[81, 200],
                stress_work=[6 * R, 5 * C],
                factors=[list(x) for x in FACTORS],
                cold_subceilings=cold_subceilings,
                endpoint_decomposition=endpoint_decomposition,
                dense_exponent_identities=dense_identities)


if __name__ == '__main__':
    print(json.dumps(certificate(), sort_keys=True, indent=2))
