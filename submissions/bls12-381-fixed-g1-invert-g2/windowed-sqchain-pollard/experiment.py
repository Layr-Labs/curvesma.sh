"""Toy witness for the windowed squaring-backbone Pollard schedule.

The group is the additive model of a prime-order cyclic group: the discrete
log is a secret x mod n, the walk state (a, b) stands for a + x*b, and a
repeated state value is the collision. Steps whose index is divisible by 3
double the element (cyclotomic square). Other steps add a jump from a table
built by a shared doubling chain (the squaring backbone). This file checks
recovery correctness and the fraction of walks that collide by the classical
Pollard step bound. It says nothing about BLS12-381 field costs.
"""

from __future__ import annotations

import json
import math
import sys


def jumps(n: int) -> list[tuple[int, int]]:
    # Shared doubling backbone: powers 2^i of the two bases, then the sums
    # of one power from each. (u, v) means multiply by γ^u * τ^v.
    backbone = [1, 2, 4, 8]
    table = []
    for u in backbone:
        for v in backbone:
            table.append((u % n, v % n))
    return table


def solve(a1: int, b1: int, a2: int, b2: int, n: int) -> int | None:
    num = (a1 - a2) % n
    den = (b2 - b1) % n
    if den == 0:
        return None
    return (num * pow(den, n - 2, n)) % n


def walk_once(n: int, secret: int, seed: int, steps: int) -> str:
    table = jumps(n)
    # Deterministic start from the seed: γ^{s+1} * τ^{1}, secret unknown to the update.
    a = (seed + 1) % n
    b = 1 % n
    seen: dict[int, tuple[int, int]] = {}
    for t in range(steps):
        e = (a + secret * b) % n
        prev = seen.get(e)
        if prev is not None:
            got = solve(a, b, prev[0], prev[1], n)
            if got is None:
                return "singular"
            if got != secret:
                return "wrong"
            return "hit"
        seen[e] = (a, b)
        if t % 3 == 0:
            a = (2 * a) % n
            b = (2 * b) % n
        else:
            u, v = table[e & 15]
            a = (a + u) % n
            b = (b + v) % n
    return "miss"


def main() -> None:
    spec = json.load(sys.stdin)
    n = int(spec["order"])
    trials = int(spec["trials"])
    if n < 3 or n > 20000 or trials < 1 or trials > 400:
        raise SystemExit("input outside declared toy bounds")
    if any(n % p == 0 for p in range(2, int(n**0.5) + 1)):
        print(json.dumps({"skipped": "composite order", "order": n}))
        return
    # Same 28/25 inflation as the production step bound.
    steps = math.ceil((28 / 25) * math.sqrt(math.pi * n / 2))
    counts = {"hit": 0, "miss": 0, "singular": 0, "wrong": 0}
    for s in range(trials):
        secret = (3 * s + 5) % n
        counts[walk_once(n, secret, s, steps)] += 1
    if counts["wrong"]:
        raise SystemExit("collision produced a wrong logarithm")
    ok = counts["hit"]
    fraction = ok / trials
    print(
        json.dumps(
            {
                "order": n,
                "trials": trials,
                "steps": steps,
                "hits": ok,
                "misses": counts["miss"],
                "singular": counts["singular"],
                "wrong": counts["wrong"],
                "fraction": fraction,
                "scope": "toy prime-order walk only; not a production cost proof",
            }
        )
    )


if __name__ == "__main__":
    main()
