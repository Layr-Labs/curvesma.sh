"""Toy check of Frobenius-class cardinality and an adding walk.

The production attack is analytic. This script only exercises a small
additive cyclic group of prime order r = 1 (mod 12), modelling GT with
the map x |-> p x that stands in for Frobenius. It does not measure
BLS12-381 costs or production collision probabilities.
"""

from __future__ import annotations

import json
import math
import sys


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def find_primitive_root(r: int) -> int:
    order = r - 1
    factors = []
    n = order
    p = 2
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        factors.append(n)
    for g in range(2, r):
        if all(pow(g, order // q, r) != 1 for q in factors):
            return g
    raise RuntimeError("no primitive root")


def class_rep(x: int, pmod: int, r: int) -> int:
    best = x
    y = x
    for _ in range(11):
        y = (y * pmod) % r
        if y < best:
            best = y
    return best


def run(order: int) -> dict:
    r = int(order)
    if r % 12 != 1 or not is_prime(r):
        raise ValueError("order must be a prime congruent to 1 mod 12")
    gen = find_primitive_root(r)
    pmod = pow(gen, (r - 1) // 12, r)
    if pow(pmod, 12, r) != 1 or pow(pmod, 6, r) == 1:
        raise AssertionError("p-map does not have order 12")

    seen = {}
    n_classes = 0
    for x in range(1, r):
        orbit = set()
        y = x
        for _ in range(12):
            orbit.add(y)
            y = (y * pmod) % r
        if len(orbit) != 12:
            raise AssertionError("orbit size %d for %d" % (len(orbit), x))
        if 0 in orbit:
            raise AssertionError("orbit contains 0")
        rep = min(orbit)
        if rep not in seen:
            seen[rep] = True
            n_classes += 1
    expected = (r - 1) // 12
    if n_classes != expected:
        raise AssertionError("class count %d != %d" % (n_classes, expected))

    multipliers = [((i * 17 + 3) % (r - 1)) + 1 for i in range(16)]
    x = 1
    visited = {}
    cap = int(8 * math.isqrt(expected) + 16)
    collided = False
    steps = 0
    for step in range(1, cap + 1):
        rep = class_rep(x, pmod, r)
        if rep in visited:
            collided = True
            steps = step
            break
        visited[rep] = step
        h = rep % 16
        x = (rep + multipliers[h]) % r
        if x == 0:
            x = 1
    if not collided:
        raise AssertionError("no class collision within cap %d" % cap)
    return {
        "order": r,
        "classes": n_classes,
        "expected_classes": expected,
        "orbit_ok": True,
        "collided": collided,
        "steps": steps,
        "cap": cap,
    }


if __name__ == "__main__":
    spec = json.load(sys.stdin)
    print(json.dumps(run(spec["order"]), sort_keys=True))
