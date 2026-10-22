"""Toy check of parallel distinguished-point rho on inverse classes {x, 1/x}.

Works in the prime-order-l subgroup of F_q^* (a cyclic group where, as in GT,
inversion is an automorphism).  Counts every walk step, including steps wasted
in fruitless cycles and after the last distinguished point, and reports the
empirical success rate within the claimed budget
    K = 1.25 * sqrt(2 * (l/2) * ln(1/0.61)).
Only evidence at toy scale; it does not establish production behaviour.
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
        for c in range(2, 2000, 2):
            q = c * l + 1
            if is_prime(q):
                while True:
                    g = pow(rng.randrange(2, q - 1), (q - 1) // l, q)
                    if g != 1:
                        return l, q, g


def trial(l, q, g, rng, partitions, walks, theta_bits):
    x = rng.randrange(1, l)
    h = pow(g, x, q)
    table = []
    for _ in range(partitions):
        a, b = rng.randrange(l), rng.randrange(l)
        table.append((pow(g, a, q) * pow(h, b, q) % q, a, b))
    mask = (1 << theta_bits) - 1
    stuck = 20 << theta_bits
    seen = {}
    steps = 0

    def rep(y, a, b):
        yi = pow(y, q - 2, q)
        if yi < y:
            return yi, (-a) % l, (-b) % l
        return y, a, b

    def fresh():
        a, b = rng.randrange(l), rng.randrange(l)
        y, a, b = rep(pow(g, a, q) * pow(h, b, q) % q, a, b)
        return [y, a, b, 0]

    state = [fresh() for _ in range(walks)]
    while True:
        for w in state:
            y, a, b, n = w
            j = (y * 0x9E3779B1 >> 7) % partitions
            # Wiener-Zuccherato look-ahead: avoid the dominant 2-cycles.
            for t in range(partitions):
                jj = (j + t) % partitions
                m, ma, mb = table[jj]
                ny, na, nb = rep(y * m % q, (a + ma) % l, (b + mb) % l)
                if (ny * 0x9E3779B1 >> 7) % partitions != jj:
                    break
            steps += 1
            n += 1
            if (ny * 0x2545F491 >> 3) & mask == 0:
                if ny in seen:
                    pa, pb = seen[ny]
                    if (pb - nb) % l:
                        xr = (na - pa) * pow(pb - nb, l - 2, l) % l
                        if pow(g, xr, q) == h:
                            return steps
                else:
                    seen[ny] = (na, nb)
                w[:] = fresh()
                continue
            if n > stuck:
                w[:] = fresh()
                continue
            w[:] = [ny, na, nb, n]


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 26))
    trials = int(args.get("trials", 100))
    seed = int(args.get("seed", 1))
    walks = int(args.get("walks", 16))
    theta_bits = int(args.get("theta_bits", 4))
    rng = random.Random(seed)
    l, q, g = setup(bits, rng)
    budget = 1.25 * math.sqrt(2 * (l / 2) * math.log(1 / 0.61))
    expected = math.sqrt(math.pi * l / 4)
    results = [trial(l, q, g, rng, 64, walks, theta_bits) for _ in range(trials)]
    ok = sum(1 for s in results if s <= budget)
    print(json.dumps({
        "l": l, "q": q, "trials": trials, "walks": walks,
        "distinguished_fraction": 2.0 ** -theta_bits,
        "approx_dp_lag_steps": walks * 2 ** theta_bits,
        "mean_steps": sum(results) / trials,
        "mean_over_sqrt_pi_l_over_4": sum(results) / trials / expected,
        "budget_steps": budget,
        "success_within_budget": ok / trials,
        "random_mapping_success_at_k0": 0.39,
        "random_mapping_success_at_budget": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
