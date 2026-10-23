"""Toy O(1)-memory Floyd (tortoise-hare) rho for discrete log.

A fixed-budget Pollard rho with FLOYD cycle detection: two pointers (tortoise
one step, hare two steps) walk a deterministic R-adding map f(x)=x*M_{H(x)} on a
cyclic prime-order group; when they meet, g^{a_t}z^{b_t}=g^{a_h}z^{b_h} solves the
log. Memory is O(1) (two elements and their exponent pairs) -- natively
addressable, no table, no distinguished points, no sorting, no large/segmented
memory. The step budget is a hard deterministic cap.

Group: the prime-order-l subgroup of F_q* (cyclic, order l), standing in for GT.
Success = a usable collision found within the hard budget of 1.25*k0 tortoise
steps (3 group multiplications each), a Monte-Carlo estimate of the bounded
success probability (proof.md lower-bounds it at >= 0.39). Toy scale only.
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


def trial(l, q, g, rng, partitions, budget_tort, counters):
    """One Floyd rho; return True iff a usable collision is found within budget."""
    x = rng.randrange(1, l)
    h = pow(g, x, q)
    table = []
    for _ in range(partitions):
        ca, cb = rng.randrange(l), rng.randrange(l)
        table.append((pow(g, ca, q) * pow(h, cb, q) % q, ca, cb))

    def step(u, a, b):
        j = (u * 0x9E3779B1 >> 7) % partitions
        m, ma, mb = table[j]
        counters['mults'] += 1
        return u * m % q, (a + ma) % l, (b + mb) % l

    # start tortoise and hare at the same random point
    a0, b0 = rng.randrange(l), rng.randrange(l)
    s = pow(g, a0, q) * pow(h, b0, q) % q
    tx, ta, tb = s, a0, b0
    hx, ha, hb = s, a0, b0
    for _ in range(budget_tort):
        tx, ta, tb = step(tx, ta, tb)                 # tortoise +1
        hx, ha, hb = step(hx, ha, hb)
        hx, ha, hb = step(hx, ha, hb)                 # hare +2
        if tx == hx:                                  # Floyd meeting
            if (tb - hb) % l:
                xr = (ha - ta) * pow(tb - hb, l - 2, l) % l
                if pow(g, xr, q) == h:
                    return True
            return False                              # unsolvable (b_t==b_h); rare
    return False


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 24))
    trials = int(args.get("trials", 150))
    seed = int(args.get("seed", 1))
    partitions = int(args.get("partitions", 256))
    rng = random.Random(seed)
    l, q, g = setup(bits, rng)
    k0 = math.sqrt(2 * l * math.log(1 / 0.61))
    budget_tort = int(math.ceil(float(args.get("budget_factor",1.1)) * k0))           # hard cap on tortoise steps
    ok = 0
    for _ in range(trials):
        c = {'mults': 0}
        if trial(l, q, g, rng, partitions, budget_tort, c):
            ok += 1
    print(json.dumps({
        "l": l, "trials": trials, "partitions_R": partitions,
        "budget_tortoise_steps": budget_tort,
        "budget_tortoise_over_sqrt_l": budget_tort / math.sqrt(l),
        "success_within_budget": ok / trials,
        "memory": "O(1): two elements + two exponent pairs",
        "random_mapping_rho_le_1.25k0": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
