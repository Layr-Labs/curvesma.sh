"""Toy DETERMINISTIC baby-step/giant-step by SORTING (no hash table).

Confirms that sort-based BSGS recovers the discrete logarithm for EVERY instance
(success probability exactly 1), using at most 2*ceil(sqrt(l)) group
multiplications.  The collision between baby and giant lists is found by sorting
the combined list on the exact group element and scanning adjacent entries -- a
deterministic, comparison-based procedure with no hash table, no probing, and no
randomness.  (Production sorts on a 512-bit fingerprint with full-element
confirmation on equal-fingerprint neighbours; see proof.md.)

Group: the prime-order-l subgroup of F_q* (cyclic, order l), standing in for GT.
"""
import json
import math
import random
import sys


def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def setup(bits, rng):
    while True:
        l = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not is_prime(l):
            continue
        for c in range(2, 4000, 2):
            q = c * l + 1
            if is_prime(q):
                g = 1
                while g == 1:
                    g = pow(rng.randrange(2, q - 1), (q - 1) // l, q)
                return l, q, g


def bsgs_sort(l, q, g, target):
    """Deterministic sort-based BSGS. Returns (a, ops) with g^a = target."""
    m = int(math.isqrt(l - 1)) + 1
    ops = 0
    baby = []                                  # (element, ('B', j))
    cur = 1
    for j in range(m):
        baby.append((cur, ('B', j)))
        if j + 1 < m:
            cur = cur * g % q
            ops += 1
    factor = pow(g, (q - 1) - (m % (q - 1)), q)  # g^{-m}
    h = (l + m - 1) // m
    giant = []                                 # (element, ('G', i))
    gamma = target
    for i in range(h):
        giant.append((gamma, ('G', i)))
        gamma = gamma * factor % q
        ops += 1
    combined = baby + giant
    combined.sort(key=lambda t: t[0])          # deterministic sort on exact element
    for idx in range(len(combined) - 1):       # scan adjacent for equal element, opposite tags
        if combined[idx][0] == combined[idx + 1][0]:
            (t1, k1), (t2, k2) = combined[idx][1], combined[idx + 1][1]
            tags = {t1: k1, t2: k2}
            if 'B' in tags and 'G' in tags:
                a = (tags['G'] * m + tags['B']) % l
                if pow(g, a, q) == target:
                    return a, ops
    raise AssertionError("no collision - impossible for a generator")


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 24))
    trials = int(args.get("trials", 200))
    seed = int(args.get("seed", 1))
    rng = random.Random(seed)
    l, q, g = setup(bits, rng)
    m = int(math.isqrt(l - 1)) + 1
    ok = 0
    max_ops = 0
    for _ in range(trials):
        x = rng.randrange(0, l)
        target = pow(g, x, q)
        a, ops = bsgs_sort(l, q, g, target)
        if a == x % l:
            ok += 1
        max_ops = max(max_ops, ops)
    print(json.dumps({
        "l": l, "trials": trials, "m_ceil_sqrt_l": m,
        "success_probability": ok / trials,
        "all_recovered": ok == trials,
        "max_ops_observed": max_ops,
        "two_sqrt_l_bound": 2 * m,
        "max_ops_le_2sqrt_l": bool(max_ops <= 2 * m),
        "method": "sort-based (no hash table)",
    }, indent=1))


if __name__ == "__main__":
    main()
