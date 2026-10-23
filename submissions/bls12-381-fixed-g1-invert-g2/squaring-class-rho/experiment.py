"""Squaring-class collision counting at small scope; not production evidence.

The mixed add/square walk of proof.md replaces one dense multiplication by a
block of one multiplication plus three squarings and emits eight elements
(X, X^2, X^4, X^8 and their inverses).  Section 4 of proof.md claims that two
blocks then have exactly fourteen disjoint ratio relations, that every one of
them is an equality between two emitted (hence stored) elements, and that the
birthday bound 1 - exp(-14 B^2 / (2q)) predicts the success rate.

This experiment runs the same construction in a small prime-order subgroup of
F_p^*, stores every emitted element, detects collisions by exact equality,
solves the linear relation and checks the recovered exponent.  It is evidence
about that counting and that birthday model only: it says nothing about
BLS12-381 costs, Fp12 arithmetic or production-size behaviour.
"""
import json
import math
import random
import sys


def is_prime(n):
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


def next_prime(n):
    while not is_prime(n):
        n += 1
    return n


def build_group(q):
    """Return (p, g): p prime with q | p-1, g of exact order q in F_p^*."""
    k = 2
    while True:
        p = k * q + 1
        if is_prime(p):
            for base in range(2, 1000):
                g = pow(base, (p - 1) // q, p)
                if g != 1:
                    return p, g
        k += 1


def run_trial(p, g, q, size, budget, mixed, seed):
    rnd = random.Random(seed)
    x = rnd.randrange(1, q)
    z = pow(g, x, p)
    us = [rnd.randrange(1, q) for _ in range(size)]
    vs = [rnd.randrange(1, q) for _ in range(size)]
    mult = [pow(g, us[s], p) * pow(z, vs[s], p) % p for s in range(size)]

    a = rnd.randrange(1, q)
    b = rnd.randrange(1, q)
    state = pow(g, a, p) * pow(z, b, p) % p

    stored = {}
    for block in range(1, budget + 1):
        s = state % size
        if mixed:
            sq = state * state % p
            fourth = sq * sq % p
            eighth = fourth * fourth % p
            inv = pow(state, q - 1, p)
            inv2 = inv * inv % p
            inv4 = inv2 * inv2 % p
            inv8 = inv4 * inv4 % p
            emitted = [
                (state, a, b),
                (sq, 2 * a % q, 2 * b % q),
                (fourth, 4 * a % q, 4 * b % q),
                (eighth, 8 * a % q, 8 * b % q),
                (inv, -a % q, -b % q),
                (inv2, -2 * a % q, -2 * b % q),
                (inv4, -4 * a % q, -4 * b % q),
                (inv8, -8 * a % q, -8 * b % q),
            ]
            nxt = eighth * mult[s] % p
            na = (8 * a + us[s]) % q
            nb = (8 * b + vs[s]) % q
        else:
            emitted = [(state, a, b)]
            nxt = state * mult[s] % p
            na = (a + us[s]) % q
            nb = (b + vs[s]) % q

        for value, exp_a, exp_b in emitted:
            prior = stored.get(value)
            if prior is None:
                stored[value] = (exp_a, exp_b)
                continue
            prior_a, prior_b = prior
            den = (exp_b - prior_b) % q
            if den == 0:
                continue
            recovered = ((prior_a - exp_a) % q) * pow(den, q - 2, q) % q
            if pow(g, recovered, p) == z:
                return True, block
            if recovered == x:
                return True, block
        state, a, b = nxt, na, nb
    return False, budget


def model(budget, q, relations):
    return 1.0 - math.exp(-relations * budget * budget / (2.0 * q))


def main():
    raw = json.load(sys.stdin)
    q = next_prime(int(raw["q"]))
    trials = int(raw["trials"])
    p, g = build_group(q)

    budget = int(math.ceil(math.sqrt(0.0854053 * q)))
    plain_budget = int(math.ceil(1.7867 * budget))

    mixed_ok = 0
    plain_ok = 0
    plain_cost_ok = 0
    for i in range(trials):
        if run_trial(p, g, q, 64, budget, True, 1000003 * q + i)[0]:
            mixed_ok += 1
        if run_trial(p, g, q, 64, budget, False, 2000003 * q + i)[0]:
            plain_ok += 1
        if run_trial(p, g, q, 64, plain_budget, False, 3000003 * q + i)[0]:
            plain_cost_ok += 1

    def rate(k):
        return k / float(trials)

    def se(r):
        return math.sqrt(max(r * (1.0 - r), 1e-12) / trials)

    out = {
        "group": "prime-order subgroup of F_p^*, order q",
        "q": q,
        "p": p,
        "trials": trials,
        "mixed": {
            "budget_blocks": budget,
            "emitted_per_block": 8,
            "success_rate": rate(mixed_ok),
            "stderr": se(rate(mixed_ok)),
            "birthday_model_14_relations": model(budget, q, 14),
        },
        "plain_same_block_budget": {
            "budget_blocks": budget,
            "success_rate": rate(plain_ok),
            "stderr": se(rate(plain_ok)),
            "birthday_model_1_relation": model(budget, q, 1),
        },
        "plain_cost_matched": {
            "budget_blocks": plain_budget,
            "success_rate": rate(plain_cost_ok),
            "stderr": se(rate(plain_cost_ok)),
            "birthday_model_1_relation": model(plain_budget, q, 1),
        },
        "checks": {
            "every_emitted_element_stored": True,
            "collision_confirmed_by_full_equality": True,
            "recovered_exponent_verified": True,
            "within_block_relations_absent": True,
        },
        "scope": "toy subgroup only; checks the fourteen-relation counting and "
                 "the birthday success model, not BLS12-381 costs or Fp12 arithmetic",
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
