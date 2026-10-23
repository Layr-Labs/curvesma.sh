#!/usr/bin/env python3
"""Toy model of the Frobenius-class block birthday walk.

The scored attack works in GT = mu_r with the cheap maps

    sq : X -> X^2                 cyclotomic squaring
    fr : X -> X^p                 Frobenius, exponent multiplier p
    cj : X -> X^{p^6} = X^{-1}    free conjugation (p^6 = -1 mod r)

and emits, for each walk block, the set { X^{eps p^k 2^d} } with
d = 0..depth, k = 0..5, eps = +-1, which is 12 (depth + 1) group elements per
block at the cost of one dense multiplication plus `depth` squarings plus
5 (depth + 1) Frobenius maps.

This experiment replaces GT by the additive cyclic group Z_q of prime order q.
A walk element is its own exponent a + b x mod q; squaring is x -> 2x, Frobenius
is x -> rho x for an element rho of order 12 (so rho^6 = -1 plays the role of
p^6 = -1), and conjugation is x -> -x. Element equality, the collision equation
and the birthday count are then exact, so the experiment measures the
combinatorial model at its own scope. It says nothing about field-arithmetic
costs, production-size behaviour or the production bound.

It also checks, on the real BLS12-381 profile constants, that p^6 = -1 (mod r)
and that the 384 emitted multipliers are distinct mod r.
"""

import json
import math
import random
import sys

P_MOD = int(
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
    "1eabfffeb153ffffb9feffffffffaaab",
    16,
)
R_ORDER = int(
    "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    16,
)


def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        y = pow(a, d, n)
        if y == 1 or y == n - 1:
            continue
        for _ in range(s - 1):
            y = y * y % n
            if y == n - 1:
                break
        else:
            return False
    return True


def find_group(bits, rng):
    """Prime q near 2^bits with q = 1 (mod 12), plus an element rho of order 12."""
    q = (1 << bits) + 1
    while True:
        if q % 12 == 1 and is_prime(q):
            break
        q += 12
    while True:
        a = rng.randrange(2, q - 1)
        rho = pow(a, (q - 1) // 12, q)
        if rho != 1 and pow(rho, 6, q) == q - 1 and pow(rho, 12, q) == 1:
            return q, rho


def check_real_constants(depth):
    """Exact checks on the real profile constants used by the scored attack."""
    conj_ok = pow(P_MOD, 6, R_ORDER) == R_ORDER - 1
    pow2 = [pow(2, d, R_ORDER) for d in range(depth + 1)]
    powp = [pow(P_MOD, k, R_ORDER) for k in range(6)]
    mult = set()
    for d in range(depth + 1):
        for k in range(6):
            for s in (1, -1):
                mult.add(s * powp[k] % R_ORDER * pow2[d] % R_ORDER)
    return {
        "p6_equals_minus_one_mod_r": conj_ok,
        "squaring_multipliers_distinct": len(set(pow2)) == depth + 1,
        "emitted_multipliers": len(mult),
        "emitted_multipliers_expected": 12 * (depth + 1),
        "emitted_multipliers_all_distinct": len(mult) == 12 * (depth + 1),
    }


def run_trial(q, rho, depth, blocks, rng):
    """One full attempt. Returns (success, collisions, emitted, predicted_lambda)."""
    x_true = rng.randrange(1, q)

    steps = [(rng.randrange(1, 1 << 16), rng.randrange(1, 1 << 16))
             for _ in range(32)]
    powp = [pow(rho, k, q) for k in range(6)]
    pow2 = [pow(2, d, q) for d in range(depth + 1)]

    table = {}
    a, b = 0, 1                      # X_0 = h = g^x, so the exponent is 1 * x
    emitted = 0
    collisions = 0
    solved = None

    for i in range(1, blocks + 1):
        c, d = steps[rng.randrange(32)]
        a += c
        b += d
        xval = (a + b * x_true) % q   # X_i = g^{a_i} h^{b_i}
        for dd in range(depth + 1):
            base = pow2[dd] * xval % q        # X_i^{2^d}
            for k in range(6):
                for s in (1, -1):
                    m = s * powp[k] % q * pow2[dd] % q
                    val = m * base % q
                    emitted += 1
                    hit = table.get(val)
                    if hit is None:
                        table[val] = (i, m, a % q, b % q)
                        continue
                    j, mv, aj, bj = hit
                    if j == i:
                        continue              # deterministic within-block relation
                    collisions += 1
                    if solved is not None:
                        continue
                    den = (m * (b % q) - mv * bj) % q
                    if den == 0:
                        continue
                    cand = (mv * aj - m * (a % q)) % q * pow(den, q - 2, q) % q
                    if cand == x_true:
                        solved = cand

    lam = emitted * emitted / (2.0 * q)
    return solved is not None, collisions, emitted, lam


def main(argv):
    cfg = json.loads(argv[1]) if len(argv) > 1 else {}
    q_bits = int(cfg.get("q_bits", 30))
    trials = int(cfg.get("trials", 24))
    depth = int(cfg.get("depth", 31))
    q_bits = max(20, min(34, q_bits))
    trials = max(4, min(40, trials))
    depth = max(1, min(31, depth))

    rng = random.Random(0x5EED_1234 + 1000 * q_bits + trials + depth)
    q, rho = find_group(q_bits, rng)

    per_block = 12 * (depth + 1)
    target = math.sqrt(2.0 * 0.59784 * q)      # lambda = 0.598 at success 0.45
    blocks = max(1, int(round(target / per_block)))

    results = []
    for _ in range(trials):
        ok, coll, emitted, lam = run_trial(q, rho, depth, blocks, rng)
        results.append((ok, coll, emitted, lam))

    succ = sum(1 for r in results if r[0])
    colls = sum(r[1] for r in results)
    emitted = results[0][2]
    lam = results[0][3]
    predicted_success = 1.0 - math.exp(-lam)

    out = {
        "toy_group": {"prime": q, "bits": q_bits, "rho_order": 12,
                      "rho_pow6_is_minus_one": pow(rho, 6, q) == q - 1},
        "depth": depth,
        "elements_per_block": per_block,
        "blocks_per_attempt": blocks,
        "emitted_per_attempt": emitted,
        "predicted_collisions_per_attempt": round(lam, 4),
        "predicted_success": round(predicted_success, 4),
        "measured_collisions_total": colls,
        "measured_collisions_per_attempt": round(colls / trials, 4),
        "trials": trials,
        "successes": succ,
        "measured_success": round(succ / trials, 4),
        "real_constants": check_real_constants(depth),
    }
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
