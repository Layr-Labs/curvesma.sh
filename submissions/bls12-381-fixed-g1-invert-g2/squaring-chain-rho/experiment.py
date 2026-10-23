"""Toy check of the block update used in proof.md.

The production walk lives in GT inside Fp12. This experiment keeps only the
combinatorial content: a cyclic group of prime order `order` embedded in Z_p*,
blocks of one multiplication by an adder followed by three squarings, a
distinguished predicate on the low bits of the element, and the discrete log
recovered from the linear congruence in the exponents.

State (y, alpha, beta) satisfies y = g^alpha * z^beta. A block multiplies by
g^(u_s) * z and then squares three times, so it emits four elements.

Scope: small orders only. It measures neither Fp12 cost, nor memory, nor the
255-bit random-mapping tail, and it cannot distinguish the seven-ratio
accounting of proof.md from a naive per-element sample count.
"""

import json
import sys
from math import gcd, isqrt

BRANCHES = 20
CHAIN = 3  # squarings per block, so a block emits CHAIN + 1 elements


def is_prime(m):
    if m < 2:
        return False
    if m % 2 == 0:
        return m == 2
    d = 3
    while d * d <= m:
        if m % d == 0:
            return False
        d += 2
    return True


def prime_one_mod(order):
    k = 1
    while True:
        p = k * order + 1
        if is_prime(p):
            return p
        k += 1


def order_generator(order, p):
    h = 2
    while h < p:
        g = pow(h, (p - 1) // order, p)
        if g != 1:
            return g
        h += 1
    raise RuntimeError("no generator")


def walk(order, x, g, p, u, adders, budget):
    """Return x when the walk's distinguished-point collision recovers it."""
    z = pow(g, x, p)
    y = 1
    alpha = 0
    beta = 0
    store = {}
    for _ in range(budget):
        s = y % BRANCHES
        y = y * adders[s] % p * z % p
        alpha = (alpha + u[s]) % order
        beta = (beta + 1) % order
        for e in range(CHAIN + 1):
            if e:
                y = y * y % p
                alpha = 2 * alpha % order
                beta = 2 * beta % order
            if y % 8:
                continue
            previous = store.get(y)
            if previous is None:
                store[y] = (alpha, beta)
                continue
            a0, b0 = previous
            d = (beta - b0) % order
            if d and gcd(d, order) == 1:
                cand = (a0 - alpha) * pow(d, -1, order) % order
                if pow(g, cand, p) == z:
                    return cand
    return None


def run(order):
    if order < 17 or not is_prime(order):
        return {"order": order, "tested": 0, "recovered": 0,
                "note": "order must be a prime at least 17"}
    p = prime_one_mod(order)
    g = order_generator(order, p)
    u = [pow(2, 12 * s, order) + 1 for s in range(BRANCHES)]
    adders = [pow(g, u[s], p) for s in range(BRANCHES)]
    budget = 12 * isqrt(order) + 128
    if order <= 257:
        targets = list(range(1, order))
    else:
        step = max(1, order // 64)
        targets = list(range(1, order, step))[:64]
    recovered = 0
    for x in targets:
        if walk(order, x, g, p, u, adders, budget) == x:
            recovered += 1
    return {"order": order, "modulus": p, "generator": g,
            "tested": len(targets), "recovered": recovered,
            "budget_blocks": budget,
            "success": round(recovered / len(targets), 6)}


def main():
    order = 101
    if len(sys.argv) > 1:
        order = int(sys.argv[1])
    print(json.dumps(run(order), sort_keys=True))


if __name__ == "__main__":
    main()
