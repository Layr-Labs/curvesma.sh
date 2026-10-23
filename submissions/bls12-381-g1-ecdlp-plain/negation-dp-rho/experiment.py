#!/usr/bin/env python3
"""Score ceiling and toy recovery for negation-canonical distinguished-point rho.

Production integers are recomputed from the protected subgroup order.
The walk runs only on Z/nZ at the supplied odd order.
"""

from __future__ import annotations

import json
import sys
from decimal import Decimal, ROUND_CEILING, getcontext

P_BLS = int(
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
    "1eabfffeb153ffffb9feffffffffaaab",
    16,
)
R_BLS = int(
    "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    16,
)

FMUL = 383
C = 2000
PARTITIONS = 2048
BATCH = 4096
# Fermat p-2 is 380 squares + 228 multiplications = 608. Montgomery's
# trick costs 1 inversion + 3(B-1) multiplications, i.e. 3 + 605/B per
# inverse. The reference 6M already counts 3. The difference is the excess.
EXCESS = (605 * FMUL + BATCH - 1) // BATCH
E = 249
SIDE_UNITS = 2**100
SCORE_NUM = 127424633
SCORE_DEN = 1000000
TOY_PARTS = 64
TOY_WALKS = 24


def isqrt_ceil(n: int) -> int:
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            break
        x = y
    if x * x < n:
        x += 1
    return x


def production_steps() -> tuple[int, int, int]:
    if (R_BLS - 1) % 2 != 0:
        raise SystemExit("r-1 is odd")
    n_classes = (R_BLS - 1) // 2
    # pi < 355/113, so pi*N'/2 < 355*N'/226.
    t2 = (355 * n_classes + 225) // 226
    mean = isqrt_ceil(t2)
    num = (PARTITIONS + 1) * (2 * PARTITIONS + 1)
    den = PARTITIONS * (2 * PARTITIONS)
    steps = (mean * num + den - 1) // den
    return n_classes, mean, steps


def score_ceiling(steps: int) -> int:
    getcontext().prec = 80
    wnum = steps * (C + E) + SIDE_UNITS * C
    log = (Decimal(wnum) / Decimal(C)).ln() / Decimal(2).ln()
    scaled = (log * Decimal(SCORE_DEN)).to_integral_value(rounding=ROUND_CEILING)
    return int(scaled)


def check_ledger() -> dict:
    if EXCESS != 57:
        raise SystemExit("batch excess drifted")
    if E != 192 + EXCESS:
        raise SystemExit("instruction ledger drifted")
    if P_BLS.bit_length() != 381 or R_BLS.bit_length() != 255:
        raise SystemExit("parameter width")
    exp = P_BLS - 2
    squares = exp.bit_length() - 1
    extra = bin(exp).count("1") - 1
    if squares != 380 or extra != 228 or squares + extra != 608:
        raise SystemExit("Fermat chain drifted")
    n_classes, mean, steps = production_steps()
    got = score_ceiling(steps)
    if got != SCORE_NUM:
        raise SystemExit("score ceiling %s != %s" % (got, SCORE_NUM))
    # 2048 jumps, two 255-bit scalar multiplications, <= 508 adds each.
    pre = 2048 * 2 * 508 + 508
    if pre > 2**21:
        raise SystemExit("preprocessing budget")
    return {
        "classes": str(n_classes),
        "mean": str(mean),
        "steps": str(steps),
        "score_numerator": got,
        "score_denominator": SCORE_DEN,
        "excess_instructions": EXCESS,
        "extra_instructions": E,
        "preprocessing_adds": pre,
    }


def splitmix(x: int) -> int:
    mask = (1 << 64) - 1
    z = (x + 0x9E3779B97F4A7C15) & mask
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & mask
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & mask
    return (z ^ (z >> 31)) & mask


def canon_neg(e: int, a: int, b: int, n: int) -> tuple[int, int, int]:
    """Representative in {0} union the lower half, folding the sign into (a, b)."""
    e %= n
    a %= n
    b %= n
    flipped = (n - e) % n
    if e != 0 and flipped < e:
        return flipped, (-a) % n, (-b) % n
    return e, a, b


def curve_parity() -> dict:
    """On an odd prime field exactly one of y and p-y is even, and negation stays on the curve."""
    p = 12253
    point = None
    for x in range(1, 80):
        rhs = (pow(x, 3, p) + 4) % p
        for y in range(p):
            if (y * y) % p == rhs:
                point = (x, y)
                break
        if point:
            break
    if point is None:
        raise SystemExit("no toy point")
    x, y = point
    y2 = (p - y) % p
    if (y & 1) == (y2 & 1):
        raise SystemExit("parity rule failed")
    if (y * y) % p != (y2 * y2) % p:
        raise SystemExit("negation left the curve")
    return {"p": p, "x": x, "even_y": y if (y & 1) == 0 else y2}


def recover_once(n: int, secret: int, seed: int) -> tuple[int, int]:
    js = []
    for j in range(TOY_PARTS):
        c = splitmix(seed + 2 * j + 1) % (n - 1) + 1
        d = splitmix(seed + 2 * j + 2) % n
        js.append((c, d))
    walks = []
    for i in range(TOY_WALKS):
        a = splitmix(seed + 10007 + 2 * i) % n
        b = splitmix(seed + 20011 + 2 * i) % n
        if b == 0 and a == 0:
            b = 1
        e, a, b = canon_neg((a + b * secret) % n, a, b, n)
        walks.append([a, b, e])
    table: dict[int, tuple[int, int]] = {}
    cap = 25 * n
    steps = 0
    while steps < cap:
        for w in walks:
            a, b, e = w
            j = e % TOY_PARTS
            for attempt in range(2):
                c, d = js[j]
                e2, a2, b2 = canon_neg(
                    (e + c + d * secret) % n,
                    (a + c) % n,
                    (b + d) % n,
                    n,
                )
                if e2 % TOY_PARTS != j or attempt == 1:
                    a, b, e = a2, b2, e2
                    break
                j = (j + 1) % TOY_PARTS
            if e != (a + b * secret) % n:
                raise SystemExit("canonical state left the coefficient line")
            w[0], w[1], w[2] = a, b, e
            steps += 1
            # Partition is bits 0..5 (TOY_PARTS = 64). The distinguished
            # bit is bit 6, disjoint from the partition, as in production.
            if e == 0 or ((e >> 6) & 1) != 0:
                continue
            prev = table.get(e)
            if prev is None:
                table[e] = (a, b)
                continue
            a1, b1 = prev
            den = (b1 - b) % n
            if den == 0:
                continue
            got = ((a - a1) % n) * pow(den, -1, n) % n
            if got == secret:
                return steps, 1
    return steps, 0


def recover(n: int, secret: int, seed: int) -> tuple[int, int]:
    total = 0
    for attempt in range(6):
        steps, ok = recover_once(n, secret, seed + 104729 * attempt)
        total += steps
        if ok:
            return total, 1
    return total, 0


def toy(order: int, trials: int, seed: int) -> dict:
    if order < 3:
        raise SystemExit("order must be an integer >= 3")
    n_classes = (order - 1) // 2
    mean = (3.141592653589793 * n_classes / 2) ** 0.5
    ratios = []
    for t in range(trials):
        secret = 1 + splitmix(seed + 17 * t + 3) % (order - 1)
        steps, ok = recover(order, secret, seed + 1000 * (t + 1))
        if not ok:
            raise SystemExit("missed secret %s" % secret)
        ratios.append(steps / mean)
    ratios.sort()
    # Index of the smallest observation that is at least the 0.39 quantile.
    q_index = min(len(ratios) - 1, int((len(ratios) - 1) * 0.39))
    return {
        "order": order,
        "trials": trials,
        "recovered": trials,
        "quantile_0_39_over_mean": ratios[q_index],
        "median_over_mean": ratios[len(ratios) // 2],
        "max_over_mean": ratios[-1],
    }


def main() -> None:
    spec = json.load(sys.stdin)
    ledger = check_ledger()
    parity = curve_parity()
    walked = toy(int(spec["order"]), int(spec["trials"]), int(spec["seed"]))
    json.dump(
        {"ledger": ledger, "parity": parity, "walk": walked},
        sys.stdout,
        sort_keys=True,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
