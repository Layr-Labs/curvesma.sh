"""Toy exhaustive check of negation-canonicalized Pollard rho with fruitless-cycle
escape; models <P> as Z_n under addition. Does not validate production curve
arithmetic, instruction costs, or full-size timing -- see proof.md."""
import json
import sys
import random


def run(n, w, unit):
    half = w // 2
    assert pow(unit, half, n) == 1
    powers = [pow(unit, i, n) for i in range(half)]
    group = sorted(set(powers) | set((n - g) % n for g in powers))
    assert len(group) == w, (n, w, unit, group)

    L = 8
    rng = random.Random(1234567)
    branch_c = [rng.randrange(1, n) for _ in range(L)]
    branch_d = [rng.randrange(1, n) for _ in range(L)]
    escape_c, escape_d = rng.randrange(1, n), rng.randrange(1, n)

    def canon(value, a, b):
        best = None
        for g in group:
            cand = ((g * value) % n, (g * a) % n, (g * b) % n)
            if best is None or cand[0] < best[0]:
                best = cand
        return best

    def step(value, a, b, x):
        h = value % L
        return canon((value + branch_c[h] + branch_d[h] * x) % n,
                      (a + branch_c[h]) % n, (b + branch_d[h]) % n)

    def escape(value, a, b, x):
        return canon((value + escape_c + escape_d * x) % n,
                      (a + escape_c) % n, (b + escape_d) % n)

    ok_count = 0
    for x in range(n):
        value, a, b = canon(1, 1, 0)
        table = {value: (a, b)}
        recovered = None
        for _ in range(6 * n + 20):
            value, a, b = step(value, a, b, x)
            if value in table:
                a0, b0 = table[value]
                denom = (b0 - b) % n
                if denom == 0:
                    value, a, b = escape(value, a, b, x)
                    table[value] = (a, b)
                    continue
                recovered = ((a - a0) * pow(denom, n - 2, n)) % n
                break
            table[value] = (a, b)
        if recovered == x % n:
            ok_count += 1
    return ok_count, n


if __name__ == "__main__":
    params = json.load(sys.stdin)
    n, w, unit = params["order"], params["w"], params["unit"]
    ok_count, total = run(n, w, unit)
    print(json.dumps({
        "group": "integers modulo n under addition, with multiply-by-unit automorphism",
        "order": n, "w": w, "recovered": ok_count, "total": total,
        "all_recovered": ok_count == total,
        "scope": "toy group only; tests the collision/canonicalization/escape and "
                 "linear-solve logic, not production curve arithmetic or cost bounds",
    }))
