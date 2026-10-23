#!/usr/bin/env python3
"""Own toy arithmetic evidence. Exponents are visible in the toy, NOT in the attack.
No production BLS cost, partition distribution, or production probability is tested.
"""
import argparse
import json
from math import isqrt
import random


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
    if not args.alignment:
        parser.error('pipeline evidence moved to mechanism_toy.py; use --alignment here')
    print(json.dumps(alignment_fixture(), indent=2))
