"""Toy negation-quotient kangaroo check; not evidence of production cost."""
import json
import math
import sys


def quotient(k, n):
    """Canonical representative of the pair {k, -k} in Z_n."""
    k %= n
    neg = (-k) % n
    return k if k <= neg else neg


def recover(n, secret, seed, attempts=16):
    """Two-animal walk on the negation quotient of the additive group Z_n."""
    width = n // 2
    trail_len = math.isqrt(width) + 2
    for attempt in range(attempts):
        state = 1 + ((seed + 7919 * attempt) % 104729)
        count = max(2, min(trail_len, 64))
        jumps = []
        for _ in range(count):
            state = (state * 1103515245 + 12345) % 2147483648
            jumps.append(1 + state % (2 * trail_len))
        start = n // 2
        trail = {}
        pos = start % n
        for _ in range(trail_len):
            trail[quotient(pos, n)] = pos
            pos = (pos + jumps[quotient(pos, n) % count]) % n
        pos, offset = secret % n, 0
        for _ in range(4 * trail_len + 8):
            key = quotient(pos, n)
            if key in trail:
                tame = trail[key]
                for candidate in ((tame - offset) % n, (-tame - offset) % n):
                    if candidate == secret:
                        return candidate
            step = jumps[key % count]
            pos = (pos + step) % n
            offset += step
    return None


def main():
    order = json.load(sys.stdin)["order"]
    failures = [x for x in range(order) if recover(order, x, 12345) != x]
    if failures:
        raise AssertionError("unrecovered scalars: %r" % failures[:8])
    print(json.dumps({
        "group": "integers modulo order under addition, quotient by negation",
        "order": order,
        "checked": order,
        "all_recovered": True,
        "scope": "toy group only; does not validate production costs, memory or success bounds",
    }))


if __name__ == "__main__":
    main()
