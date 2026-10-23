"""Faithful scalar-action grid plus dense supplied-power root lifting.
Run python3 toy.py --trials 2000 --seed 20261023; stdout is captured evidence.
Input generation knows s; neither scalar_grid nor dense_lift receives it.
"""
import argparse
import hashlib
import json
from math import isqrt
import random


def prime(n):
    return n >= 2 and all(n % j for j in range(2, isqrt(n) + 1))


def ceil_sqrt(n):
    z = isqrt(n)
    return z + (z * z < n)


def curve():
    p, r, P = 221201, 12289, (169338, 64967)

    def neg(A):
        return None if A is None else (A[0], -A[1] % p)

    def add(A, B):
        if A is None:
            return B
        if B is None:
            return A
        x, y = A
        u, v = B
        if x == u and (y + v) % p == 0:
            return None
        slope = ((3 * x * x) * pow(2 * y, -1, p) if A == B
                 else (v - y) * pow(u - x, -1, p)) % p
        z = (slope * slope - x - u) % p
        return z, (slope * (x - z) - y) % p

    def mul(k, A):
        assert k >= 0
        B = None
        while k:
            if k & 1:
                B = add(B, A)
            A = add(A, A)
            k >>= 1
        return B

    def key(A):
        width = (p.bit_length() + 7) // 8
        if A is None:
            return bytes(2 * width) + (1).to_bytes(8, 'big')
        return A[0].to_bytes(width, 'big') + A[1].to_bytes(width, 'big') + bytes(8)

    return p, r, P, mul, neg, key


def run_toy(trials, seed):
    p, r, P, mul, neg, key = curve()
    assert prime(p) and prime(r)
    assert P[1] ** 2 % p == (P[0] ** 3 + 4) % p
    assert mul(r, P) is None and P is not None
    assert mul(0, P) is None and mul(r + 1, P) == P
    assert key(P) != key(neg(P)) and neg(neg(P)) == P
    g, d = 11, 64
    q = (r - 1) // d
    assert all(pow(g, (r - 1) // l, r) != 1 for l in (2, 3))
    m = ceil_sqrt((81 * q + 199) // 200)
    assert m * m < q
    rng = random.Random(seed)
    hits, wrong, prediction_errors = 0, 0, 0
    histogram, hit_histogram = [0] * q, [0] * q
    bit_histogram = [[0, 0] for _ in range(d.bit_length() - 1)]
    baby_histogram, giant_histogram = [0] * m, [0] * m
    work_histogram, transcript = {}, hashlib.sha256()
    for _ in range(trials):
        # Uniform nonzero s via the bijective mixed-radix exponent pair.
        k1, k2 = rng.randrange(d), rng.randrange(q)
        s = pow(g, q * k1 + k2, r)
        powers = [mul(pow(s, i, r), P) for i in range(1, d + 1)]
        T, Ad = powers[0], powers[d - 1]
        recovered_k2 = scalar_grid(P, Ad, pow(g, d, r), r, m, mul, key)
        predicted = k2 < m * m
        prediction_errors += (recovered_k2 is not None) != predicted
        histogram[k2] += 1
        recovered = None
        actions = 2 * m
        if recovered_k2 is not None:
            recovered, bits = dense_lift(P, T, powers, recovered_k2, g, d, r, mul, neg)
            actions += (d.bit_length() - 1) + 1
            wrong += recovered != s or recovered_k2 != k2
            hits += 1
            hit_histogram[k2] += 1
            baby_histogram[k2 % m] += 1
            giant_histogram[k2 // m] += 1
            for t, bit in enumerate(bits):
                bit_histogram[t][bit] += 1
        work_histogram[str(actions)] = work_histogram.get(str(actions), 0) + 1
        row = [k1, k2, recovered_k2, recovered, actions]
        transcript.update((json.dumps(row, separators=(',', ':')) + '\n').encode())
    assert wrong == 0 and prediction_errors == 0
    assert sum(histogram) == trials and sum(hit_histogram) == hits
    return dict(p=p, r=r, P=list(P), curve_b=4, g=g, d=d, q=q, m=m,
                trials=trials, seed=seed, successes=hits, failures=trials - hits,
                wrong_outputs=wrong, coverage_prediction_errors=prediction_errors,
                exact_coverage=[m * m, q], grid_action_cap=2 * m,
                success_action_cap=2 * m + d.bit_length(),
                input_power_points_per_trial=d, full_xy_negation_distinct=True,
                k2_histogram=histogram, success_k2_histogram=hit_histogram,
                baby_index_histogram=baby_histogram, giant_index_histogram=giant_histogram,
                dense_bit_histogram=bit_histogram, scalar_action_histogram=work_histogram,
                transcript_sha256=transcript.hexdigest(),
                scope='Real tiny EC arithmetic and full-key radix sort; no production execution or cost calibration.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--trials', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=20261023)
    args = parser.parse_args()
    print(json.dumps(run_toy(args.trials, args.seed), sort_keys=True, indent=2))


def scalar_grid(P, A, base, r, m, mul, key):
    """Two finite geometric point chains and full-key stable LSD radix sort."""
    records = []
    B, G = P, A
    giant_step = pow(base, -m, r)
    for j in range(m):
        records.append((key(B), 0, j))
        B = mul(base, B)
    for i in range(m):
        records.append((key(G), 1, i))
        G = mul(giant_step, G)
    # Production charges a routed 104-byte-key variant; no hash lookup is used.
    for pos in range(len(records[0][0]) - 1, -1, -1):
        buckets = [[] for _ in range(256)]
        for rec in records:
            buckets[rec[0][pos]].append(rec)
        records = [rec for bucket in buckets for rec in bucket]
    for a, b in zip(records, records[1:]):
        if a[0] == b[0] and a[1] != b[1]:
            baby, giant = (a, b) if a[1] == 0 else (b, a)
            return giant[2] * m + baby[2]
    return None


def dense_lift(P, T, powers, k2, g, d, r, mul, neg):
    """Use supplied powers only. Never exponentiate an unknown encoded scalar."""
    assert d > 0 and d & (d - 1) == 0 and (r - 1) % d == 0
    q = (r - 1) // d
    a, bits = 0, []
    for t in range(d.bit_length() - 1):
        e = d >> (t + 1)
        candidate = pow(g, e * (q * a + k2), r)
        C = mul(candidate, P)
        A = powers[e - 1]
        if C == A:
            bit = 0
        elif neg(C) == A:
            bit = 1
        else:
            raise ValueError('honest supplied-power identity violated')
        bits.append(bit)
        a += bit << t
    s = pow(g, q * a + k2, r)
    if mul(s, P) != T:
        raise ValueError('final scalar verification failed')
    return s, bits


if __name__ == '__main__':
    main()
