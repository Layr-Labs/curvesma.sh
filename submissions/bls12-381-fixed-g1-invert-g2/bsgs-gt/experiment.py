"""Toy DETERMINISTIC baby-step/giant-step for discrete log in a cyclic group.

Confirms that BSGS recovers the discrete logarithm for EVERY instance (success
probability exactly 1), and reports the operation count m + (giant steps), which
is at most 2*ceil(sqrt(l)).  No randomness, no walk, no cycles, no truncation, no
success heuristic.

Group: the prime-order-l subgroup of F_q* (cyclic, order l), standing in for GT
(cyclic, prime order r).  The production attack uses the same algorithm on GT with
g = e(P,Q0); see proof.md.
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


def bsgs(l, q, g, target):
    """Return (a, ops) with g^a = target, 0<=a<l. Deterministic; always succeeds."""
    m = int(math.isqrt(l - 1)) + 1          # m = ceil(sqrt(l))
    baby = {}
    cur = 1
    ops = 0
    for j in range(m):                      # baby steps g^0..g^{m-1}
        if cur not in baby:
            baby[cur] = j
        if j + 1 < m:
            cur = cur * g % q               # one multiplication
            ops += 1
    factor = pow(g, (q - 1) - (m % (q - 1)), q)   # g^{-m} mod q (g has order l | q-1)
    gamma = target
    h = (l + m - 1) // m                     # ceil(l/m) giant steps suffice
    for i in range(h):                       # giant steps target*g^{-i m}
        if gamma in baby:
            a = (i * m + baby[gamma]) % l
            if pow(g, a, q) == target:       # verify (rejects any false coincidence)
                return a, ops
        gamma = gamma * factor % q
        ops += 1
    raise AssertionError("BSGS failed - impossible for a generator")


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
        a, ops = bsgs(l, q, g, target)
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
    }, indent=1))


if __name__ == "__main__":
    main()
