"""Toy Pollard-rho check with orbit reduction; not evidence of production costs.

Runs a tiny Pollard rho (with negation-map canonical representatives, the
order-2 subgroup of the secp256k1 automorphism group) on the additive
cyclic group Z_order, recovering every secret. This validates the walk,
coefficient tracking, orbit reduction and collision-solving logic at toy
scale only. It does not establish production curve behavior, the GLV
automorphism, memory bounds, or operation costs.
"""
import json
import math
import random
import sys


def canonical(v, n):
    # Orbit {v, -v}: canonical = min(v, n - v); multiplier mu in {1, -1}
    # with canonical = mu * v mod n.
    if v % n <= (-v) % n:
        return v % n, 1
    return (-v) % n, n - 1


def pollard_rho(order, secret, seed):
    if secret % order == 0:
        return 0
    rng = random.Random(seed)
    r = 32
    # Precompute steps: T_j = c_j * 1 + d_j * secret  (additive notation:
    # group is Z_order, generator is 1, target is `secret`).
    steps = []
    for _ in range(r):
        c = rng.randrange(order)
        d = rng.randrange(order)
        t = (c + d * secret) % order
        steps.append((t, c, d))

    # Walk state: (rep, a, b) with rep = a + b*secret (mod order).
    a = rng.randrange(order)
    b = rng.randrange(order)
    rep, mu = canonical((a + b * secret) % order, order)
    a = (mu * a) % order
    b = (mu * b) % order
    table = {}
    # Budget: generous multiple of the heuristic expectation.
    # Toy scale: store every visited representative (no distinguished filter);
    # the production analysis uses distinguished points for memory efficiency.
    budget = int(8 * math.sqrt(math.pi * order / 4)) + 100
    for _ in range(budget):
        j = rep % r
        t, c, d = steps[j]
        a = (a + c) % order
        b = (b + d) % order
        raw = (rep + t) % order
        rep, mu = canonical(raw, order)
        a = (mu * a) % order
        b = (mu * b) % order
        if rep in table:
            a1, b1 = table[rep]
            if (a1, b1) != (a, b) and (b - b1) % order != 0:
                inv = pow((b - b1) % order, -1, order)
                x = (-(a - a1) * inv) % order
                return x
        else:
            table[rep] = (a, b)
    return None


def main():
    n = json.load(sys.stdin)["order"]
    assert 2 <= n <= 10000
    checked = 0
    for secret in range(n):
        x = None
        # Retry with different seeds on fruitless cycles (toy-scale artifact;
        # production uses cycle escape; the retry just makes the toy robust).
        for trial in range(10):
            x = pollard_rho(n, secret, seed=1234 + secret * 100 + trial)
            if x is not None:
                break
        assert x == secret % n, (n, secret, x)
        checked += 1
    print(json.dumps({
        "group": "additive cyclic group Z_order with negation orbit reduction",
        "order": n,
        "checked": checked,
        "all_recovered": True,
        "scope": "toy group only; does not validate production cost bounds",
    }))


if __name__ == "__main__":
    main()
