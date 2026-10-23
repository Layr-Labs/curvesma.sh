#!/usr/bin/env python3
"""Own toy arithmetic evidence. Exponents are visible in the toy, NOT in the attack.
No production BLS cost, partition distribution, or production probability is tested.
"""
import argparse
import json
import hashlib
from math import isqrt
import random


def pipeline(trials, seed):
    r, R, processors, theta_bits = 65521, 64, 2, 1
    assert all(r % d for d in range(2, isqrt(r) + 1))
    lam = next(pow(g, (r - 1) // 12, r) for g in range(2, r)
               if len({pow(pow(g, (r - 1) // 12, r), k, r) for k in range(12)}) == 12)
    lp = [pow(lam, k, r) for k in range(12)]
    N = (r - 1) // 12 + 1
    S = isqrt((6 * N + 4) // 5) + 1
    rng = random.Random(seed)
    totals = dict(verified_successes=0, tortoise_iterations=0, evaluations=0,
                  dense_multiplications=0, starts=0, reports=0, floyd_restarts=0,
                  timeout_restarts=0, false_fingerprint_pairs=0,
                  incorrect_verified_outputs=0, class_equivariance_failures=0)

    def canonical(y):
        k = min(range(12), key=lambda j: lp[j] * y % r)
        return lp[k] * y % r, k

    def h(domain, c, bits):
        return int.from_bytes(hashlib.sha256(domain + c.to_bytes(4, 'big')).digest(), 'big') & ((1 << bits) - 1)

    for _ in range(trials):
        x = rng.randrange(1, r)  # Nonidentity branch; the identity target is exact separately.
        table = [(rng.randrange(r), rng.randrange(r)) for j in range(R)]

        def fresh():
            totals['starts'] += 1
            a, b = rng.randrange(r), rng.randrange(r)
            return ((a + b * x) % r, a, b)

        def step(state, count=True):
            y, a, b = state
            c, k = canonical(y)
            j = h(b'partition-v7', c, 6)
            aa, bb = table[j]
            inv = lp[-k % 12]
            out = ((2 * y + inv * (aa + bb * x)) % r,
                   (2 * a + inv * aa) % r, (2 * b + inv * bb) % r)
            assert out[0] == (out[1] + out[2] * x) % r
            if count:
                totals['evaluations'] += 1
                totals['dense_multiplications'] += 2
            return out

        # Check the same table induces the same successor class on all 12 representatives.
        probe = fresh()
        baseline = canonical(step(probe, False)[0])[0]
        for ell in lp:
            shifted = tuple(ell * v % r for v in probe)
            totals['class_equivariance_failures'] += canonical(step(shifted, False)[0])[0] != baseline

        walkers = []
        for i in range(processors):
            state = fresh()
            walkers.append([state, state, 0])
        reports = []
        success = False

        def normalized(state):
            y, a, b = state
            c, k = canonical(y)
            return c, lp[k] * a % r, lp[k] * b % r

        def recover(one, two):
            den = (one[2] - two[2]) % r
            if not den:
                return False
            guessed = (two[1] - one[1]) * pow(den, -1, r) % r
            # This is g^guess=z in the toy's explicitly visible exponent model.
            return guessed == x

        for t in range(S):
            i = t % processors
            tort, hare, age = walkers[i]
            tort, hare = step(tort), step(step(hare))
            totals['tortoise_iterations'] += 1
            age += 1
            norm = normalized(tort)
            distinguished = h(b'dp-v7', norm[0], theta_bits) == 0
            if distinguished:
                reports.append((h(b'fingerprint-v7', norm[0], 8), norm))
                totals['reports'] += 1
                state = fresh()
                walkers[i] = [state, state, 0]
            elif canonical(tort[0])[0] == canonical(hare[0])[0]:
                success |= recover(norm, normalized(hare))
                totals['floyd_restarts'] += 1
                state = fresh()
                walkers[i] = [state, state, 0]
            elif age >= 128 * (1 << theta_bits):
                totals['timeout_restarts'] += 1
                state = fresh()
                walkers[i] = [state, state, 0]
            else:
                walkers[i] = [tort, hare, age]
        reports.sort(key=lambda q: q[0])
        for i, (fp, one) in enumerate(reports):
            j = i + 1
            while j < len(reports) and reports[j][0] == fp:
                two = reports[j][1]
                if one[0] != two[0]:
                    totals['false_fingerprint_pairs'] += 1
                success |= recover(one, two)
                j += 1
        totals['verified_successes'] += bool(success)
    return {'trials': trials, 'seed': seed, 'prime': r, 'lambda': lam,
            'nonidentity_classes': (r - 1) // 12, 'classes_with_identity': N,
            'table_entries': R, 'processors': processors, 'theta_bits': theta_bits,
            'fingerprint_bits': 8, 'budget_tortoise_per_trial': S, **totals,
            'scope': 'Visible-exponent toy, hashed whole canonical exponent, not production four-product key. Supports mechanism and model only, not BLS costs or production probabilities.'}


def floyd_iterations(mapping, start):
    tortoise, hare = start, start
    for iterations in range(1, len(mapping) + 1):
        tortoise = mapping[tortoise]
        hare = mapping[mapping[hare]]
        if tortoise == hare:
            return iterations
    raise AssertionError('finite functional graph did not meet')


def alignment_fixture():
    r, lam = 13, 2  # lam has exact order 12 modulo this toy prime.
    powers = [pow(lam, k, r) for k in range(12)]
    assert len(set(powers)) == 12
    rows = []
    excluded = 0
    for k, ell in enumerate(powers):
        row = {'k': k, 'usable': 0, 'aligned_correct': 0, 'normalized_correct': 0,
               'unaligned_correct': 0, 'unaligned_wrong': 0, 'unaligned_undefined': 0}
        for x in range(r):
            for y in range(1, r):
                kt = min(range(12), key=lambda i: powers[i] * y % r)
                for b in range(r):
                    a = (y - b * x) % r
                    for bp in range(r):
                        ap = (pow(ell, -1, r) * y - bp * x) % r
                        den = (b - ell * bp) % r
                        if not den:
                            excluded += 1
                            continue
                        row['usable'] += 1
                        aligned = (ell * ap - a) * pow(den, -1, r) % r
                        row['aligned_correct'] += aligned == x
                        yp = (ap + bp * x) % r
                        kh = min(range(12), key=lambda i: powers[i] * yp % r)
                        A, B = powers[kt] * a % r, powers[kt] * b % r
                        Ap, Bp = powers[kh] * ap % r, powers[kh] * bp % r
                        normalized = (Ap - A) * pow((B - Bp) % r, -1, r) % r
                        row['normalized_correct'] += normalized == x
                        if b == bp:
                            row['unaligned_undefined'] += 1
                        else:
                            raw = (ap - a) * pow((b - bp) % r, -1, r) % r
                            row['unaligned_correct'] += raw == x
                            row['unaligned_wrong'] += raw != x
        rows.append(row)
    return {'prime': r, 'lambda': lam, 'x_values': r, 'orbit_elements': r - 1,
            'rows': rows, 'usable_pairs': sum(q['usable'] for q in rows),
            'aligned_correct': sum(q['aligned_correct'] for q in rows),
            'normalized_correct': sum(q['normalized_correct'] for q in rows),
            'unaligned_wrong': sum(q['unaligned_wrong'] for q in rows),
            'unaligned_undefined': sum(q['unaligned_undefined'] for q in rows),
            'zero_aligned_denominator_excluded': excluded,
            'scope': 'Exhaustive representative-alignment identity on the specified toy prime only.'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--alignment', action='store_true')
    parser.add_argument('--trials', type=int, default=1000)
    parser.add_argument('--seed', type=int, default=71)
    args = parser.parse_args()
    print(json.dumps(alignment_fixture() if args.alignment else pipeline(args.trials, args.seed), indent=2))
