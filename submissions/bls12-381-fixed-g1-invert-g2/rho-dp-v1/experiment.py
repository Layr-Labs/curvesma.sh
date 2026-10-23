"""Toy check of the exponent bookkeeping and distinguished-point collision logic.

This runs in an additive group of prime order, not on a curve and not in Fp12.
Walks iterate until a distinguished endpoint, exactly as in the production
algorithm, and collisions are resolved by the same linear equation. It is
evidence only for the algebra and walk structure, not for production costs.
"""
import json
import random
import sys


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def step(x, u, v, target, order):
    chi = (((x * 0x9E3779B97F4A7C15) ^ (x >> 7)) >> 13) & 0xFF
    kind = chi % 3
    if kind == 0:
        c = (chi >> 1) + 1
        return (x + c) % order, (u + c) % order, v
    if kind == 1:
        return (x + target) % order, u, (v + 1) % order
    return (2 * x) % order, (2 * u) % order, (2 * v) % order


def solve(order, target, dp_bits, rng, abort_limit=1 << 19):
    """Recover a with target = a in Z_order under addition (generator g = 1)."""
    if target % order == 0:
        return 0, 0, 0
    mask = (1 << dp_bits) - 1
    table = {}
    degenerate = 0
    aborts = 0
    while True:
        start_s = rng.randrange(order)
        u = start_s
        v = 1
        x = (target + start_s) % order
        steps = 0
        while x & mask:
            x, u, v = step(x, u, v, target, order)
            steps += 1
            if steps > abort_limit:
                break
        if x & mask:
            aborts += 1
            continue
        if x in table:
            u1, v1 = table[x]
            if v1 == v:
                degenerate += 1
                continue
            a = ((u1 - u) * pow(v - v1, -1, order)) % order
            if (u1 + v1 * a - x) % order != 0:
                raise AssertionError("collision equation failed")
            return a, degenerate, aborts
        table[x] = (u, v)


def main():
    cfg = json.load(sys.stdin)
    order = int(cfg["order"])
    instances = int(cfg["instances"])
    dp_bits = int(cfg.get("dp_bits", 4))
    if not is_prime(order):
        raise SystemExit("order must be prime")
    if not (2 <= dp_bits <= 16):
        raise SystemExit("dp_bits out of range")
    rng = random.Random(20261022)
    ok = 0
    total_deg = 0
    total_abort = 0
    for _ in range(instances):
        a = rng.randrange(order)
        recovered, deg, ab = solve(order, a, dp_bits, rng)
        if recovered != a:
            raise AssertionError("wrong discrete logarithm")
        ok += 1
        total_deg += deg
        total_abort += ab
    print(json.dumps({
        "group": "integers modulo order under addition (toy, g=1)",
        "order": order,
        "instances": instances,
        "dp_bits": dp_bits,
        "all_recovered": True,
        "degenerate_collisions": total_deg,
        "aborted_walks": total_abort,
        "scope": "toy group only; validates walk algebra and collision recovery, not production cost",
        "note": "degenerate collisions are retried with a fresh walk; the toy's short walks inflate this rate relative to production",
    }))


if __name__ == "__main__":
    main()
