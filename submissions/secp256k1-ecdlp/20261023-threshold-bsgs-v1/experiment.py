"""Bounded specification checks; execute only via the shared verifier runner."""
import json
import math
import sys

P = 2**256 - 2**32 - 977
R = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
BASE = (
    int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
    int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16),
)


def add(u, v, p):
    if u is None:
        return v
    if v is None:
        return u
    x, y = u
    a, b = v
    if x == a:
        if (y + b) % p == 0:
            return None
        slope = 3 * x * x * pow(2 * y, -1, p) % p
    else:
        slope = (b - y) * pow((a - x) % p, -1, p) % p
    z = (slope * slope - x - a) % p
    return z, (slope * (x - z) - y) % p


def mul(n, u, p):
    out = None
    while n:
        if n & 1:
            out = add(out, u, p)
        u = add(u, u, p)
        n >>= 1
    return out


def neg(u, p):
    return None if u is None else (u[0], -u[1] % p)


def batch_step(states, step, active, p):
    denominators = [
        (step[0] - u[0]) % p
        if used and u is not None and u[0] != step[0] else 1
        for u, used in zip(states, active)
    ]
    prefix = [1]
    for d in denominators:
        prefix.append(prefix[-1] * d % p)
    inv = pow(prefix[-1], -1, p)
    inverses = [0] * len(states)
    for t in reversed(range(len(states))):
        inverses[t] = inv * prefix[t] % p
        inv = inv * denominators[t] % p
    out = []
    for u, used, di in zip(states, active, inverses):
        if not used:
            out.append(None)
        elif u is None or u[0] == step[0]:
            out.append(add(u, step, p))
        else:
            slope = (step[1] - u[1]) * di % p
            x = (slope * slope - u[0] - step[0]) % p
            out.append((x, (slope * (u[0] - x) - u[1]) % p))
    return out


def key(u):
    return 0 if u is None else ((2 + (u[1] & 1)) << 256) + u[0]


def stream(start, step, length, chains, p):
    layers = (length + chains - 1) // chains
    states = [add(start, mul(t * layers, step, p), p) for t in range(chains)]
    records = []
    for k in range(layers):
        active = [t * layers + k < length for t in range(chains)]
        for t, used in enumerate(active):
            if used:
                records.append((key(states[t]), t * layers + k))
        states = batch_step(states, step, active, p)
    assert sorted(i for _, i in records) == list(range(length))
    return records


def radix(records, digit_bits=4):
    mask = (1 << digit_bits) - 1
    for shift in range(0, 258, digit_bits):
        counts = [0] * (mask + 1)
        for value, _ in records:
            counts[(value >> shift) & mask] += 1
        cursor, total = [], 0
        for n in counts:
            cursor.append(total)
            total += n
        assert total == len(records)
        dest = [None] * total
        for record in records:
            d = (record[0] >> shift) & mask
            dest[cursor[d]] = record
            cursor[d] += 1
        records = dest
    return records


def solve(q, base, p, order, m, chains):
    babies = radix(stream(None, base, m, chains, p))
    d = neg(mul(m, base, p), p)
    giants = radix(stream(q, d, m, chains, p))
    a = b = 0
    while a < m and b < m:
        bk, j = babies[a]
        gk, i = giants[b]
        if bk == gk:
            x = i * m + j
            assert 0 <= x < order and mul(x, base, p) == q
            return x
        if bk < gk:
            a += 1
        else:
            b += 1
    return None


def fold_mul(a, b):
    z = a * b
    for _ in range(3):
        z = (z & (2**256 - 1)) + (z >> 256) * (2**32 + 977)
    assert z < 2 * P
    return z - P if z >= P else z


def increment_three_limbs(value):
    limbs = [(value >> (64 * i)) & (2**64 - 1) for i in range(3)]
    carry = 1
    for i in range(3):
        v = limbs[i] + carry
        limbs[i], carry = v & (2**64 - 1), v >> 64
    assert not carry
    return sum(v << (64 * i) for i, v in enumerate(limbs))


def main():
    assert json.load(sys.stdin) == {}
    coverage = (39 * R + 99) // 100
    m = math.isqrt(coverage - 1) + 1
    b = 2**16
    layers = (m + b - 1) // b
    n = 2 * m
    f, g, inv, buckets = 8192, 78037248, 2130176, 2**43
    rows = {
        "padded_slots": 2 * b * layers * (6 * 4096 + 8 * 256 + 1024),
        "batch_inversions": 2 * layers * inv,
        "seeds_and_D": b * 1025 * g + 512 * g,
        "record_writes": n * (f + 1024),
        "radix": n * 6 * (3 * f + 512),
        "buckets": 12 * buckets * 192,
        "merge": n * (f + 1024),
        "exceptions": 16 * g,
        "verification": 4096 * g,
        "router_initialization": 2**73 * 8192,
        "setup_output": 2**44,
    }
    total = sum(rows.values())
    peak = 4 * m * 56 + 2**73 * 8192 + buckets * 64 + 2**28
    assert total == 83574680044327118202874855850367101169472512
    assert total**4 < 1000**4 * 2**545
    assert peak < 2**136
    assert (m - 1)**2 < coverage <= m**2 < R
    assert 100 * m**2 >= 39 * R
    assert 2 * b * layers >= n
    assert (b - 1) * layers < m
    assert 56 * 2**58 + 8192 < 2**64
    assert ((4 + 11) * 4 + 2 * 16) * 72 + 512 < f

    toy_p, toy_r, toy_base, toy_m = 211, 199, (3, 33), 9
    count = 1 + sum((y*y - x*x*x - 7) % toy_p == 0
                    for x in range(toy_p) for y in range(toy_p))
    assert count == toy_r and mul(toy_r, toy_base, toy_p) is None
    successes = 0
    for x in range(toy_r):
        result = solve(mul(x, toy_base, toy_p), toy_base, toy_p, toy_r, toy_m, 4)
        assert result == (x if x < toy_m**2 else None), (x, result)
        successes += result is not None
    assert successes == 81

    assert (BASE[1]**2 - BASE[0]**3 - 7) % P == 0 and mul(R, BASE, P) is None
    production_m = 17
    fixtures = [0, 1, 16, 17, 18, 127, 271, 272, 287, 288, 289, 290, R - 1]
    for x in fixtures:
        result = solve(mul(x, BASE, P), BASE, P, R, production_m, 4)
        assert result == (x if x < production_m**2 else None), (x, result)
    batch = [None, BASE, neg(BASE, P), mul(2, BASE, P), mul(12345, BASE, P)]
    assert batch_step(batch, BASE, [True] * 5, P) == [add(u, BASE, P) for u in batch]
    assert batch_step(batch, BASE, [False] * 5, P) == [None] * 5
    values = [0, 1, 2, 2**32 + 977, 2**64 - 1, 2**128 - 1,
              2**255, P - 2, P - 1]
    for a in values:
        for c in values:
            assert fold_mul(a, c) == a * c % P
    for v in [0, 2**64 - 1, 2**128 - 1, 2**128]:
        assert increment_three_limbs(v) == v + 1
    tied = [(7, i) for i in range(20)]
    assert radix(tied) == tied
    assert key(None) != key(BASE) != key(neg(BASE, P))
    print(json.dumps({
        "verdict": "PASS", "m": m, "covered_scalars": m*m,
        "coverage_margin_100": 100*m*m - 39*R,
        "rows_instructions": rows, "total_instructions": total, "peak_bytes": peak,
        "toy_trials": toy_r, "toy_successes": successes,
        "toy_expected_failures": toy_r - successes,
        "secp256k1_small_grid_fixtures": len(fixtures),
        "field_product_fixtures": len(values)**2, "counter_fixtures": 4,
        "scope": "bounded mechanism and exact ledger checks, not production execution or instruction measurements",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
