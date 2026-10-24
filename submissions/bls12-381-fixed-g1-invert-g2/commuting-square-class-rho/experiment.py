#!/usr/bin/env python3
"""Toy witness for the commuting-square Frobenius-class walk.

The group is the additive model of a prime-order cyclic group whose order is
1 mod 12. The secret discrete log is x. A state (a, b) stands for the
exponent a + x*b. Frobenius is multiplication of the exponent by a fixed
omega of order 12. The canonical representative is the minimum residue in
the orbit. The branch bit is a fixed hash of that canonical residue, which depends
only on the orbit. The low bit selects multiply versus square.

Multiply steps add a gamma-jump or a z-jump and then re-canonicalise, which
multiplies both exponents by omega^k. Square steps double and re-canonicalise.
Squaring commutes with the orbit map; the re-canonicalisation in this toy is
the mathematical orbit function, not the production instruction schedule.
The production schedule defers the full-element Frobenius off the square
steps. This file does not measure Fp12 instruction counts.
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
    p = 3
    while p * p <= n:
        if n % p == 0:
            return False
        p += 2
    return True


def find_omega(q: int) -> int:
    target = (q - 1) // 12
    for g in range(2, q):
        w = pow(g, target, q)
        if (
            pow(w, 12, q) == 1
            and pow(w, 6, q) != 1
            and pow(w, 4, q) != 1
            and pow(w, 3, q) != 1
        ):
            return w
    raise RuntimeError("no element of order 12")


def canon(exp: int, omega: int, q: int) -> tuple[int, int]:
    best = exp % q
    twist = 0
    cur = best
    for i in range(1, 12):
        cur = (cur * omega) % q
        if cur < best:
            best = cur
            twist = i
    return best, twist


def solve(a1: int, b1: int, a2: int, b2: int, q: int) -> int | None:
    den = (b1 - b2) % q
    if den == 0:
        return None
    return ((a2 - a1) % q) * pow(den, q - 2, q) % q


def walk_once(q: int, omega: int, secret: int, steps: int) -> str:
    # Element z = gamma^secret, represented as a=0, b=1, then canonicalised.
    a, b = 0, 1
    e = secret % q
    crep, twist = canon(e, omega, q)
    if twist:
        w = pow(omega, twist, q)
        a = (a * w) % q
        b = (b * w) % q
        e = crep
    seen: dict[int, tuple[int, int]] = {e: (a, b)}
    for _ in range(steps):
        # Stand-in for the production invariant limb: a fixed odd hash of the
        # canonical residue. Low bit 0 multiplies; the next bit selects the jump.
        h = (e * 123456789) % q
        if h & 1 == 0:
            if h & 2:
                e = (e + secret) % q
                b = (b + 1) % q
            else:
                e = (e + 1) % q
                a = (a + 1) % q
        else:
            e = (2 * e) % q
            a = (2 * a) % q
            b = (2 * b) % q
        crep, twist = canon(e, omega, q)
        if twist:
            w = pow(omega, twist, q)
            a = (a * w) % q
            b = (b * w) % q
            e = crep
        prev = seen.get(e)
        if prev is not None:
            if prev == (a, b):
                continue
            got = solve(a, b, prev[0], prev[1], q)
            if got is None:
                return "singular"
            if got != secret:
                return "wrong"
            return "hit"
        seen[e] = (a, b)
    return "miss"


def main() -> None:
    spec = json.load(sys.stdin)
    q = int(spec["order"])
    trials = int(spec["trials"])
    if q < 13 or q > 20000 or trials < 1 or trials > 400:
        raise SystemExit("input outside declared toy bounds")
    if not is_prime(q) or (q - 1) % 12 != 0:
        print(json.dumps({"skipped": "need a prime order congruent to 1 mod 12", "order": q}))
        return
    omega = find_omega(q)
    n_orbits = (q - 1) // 12
    steps = math.ceil(math.sqrt(math.pi * n_orbits / 2))
    counts = {"hit": 0, "miss": 0, "singular": 0, "wrong": 0}
    for s in range(trials):
        secret = (3 * s + 5) % q
        if secret == 0:
            secret = 1
        counts[walk_once(q, omega, secret, steps)] += 1
    if counts["wrong"]:
        raise SystemExit("collision produced a wrong logarithm")
    fraction = counts["hit"] / trials
    print(
        json.dumps(
            {
                "order": q,
                "orbits": n_orbits,
                "omega": omega,
                "trials": trials,
                "steps": steps,
                "hits": counts["hit"],
                "misses": counts["miss"],
                "singular": counts["singular"],
                "wrong": counts["wrong"],
                "fraction": fraction,
                "scope": "toy prime-order orbit walk only; not a production cost proof",
            }
        )
    )


if __name__ == "__main__":
    main()
