"""Toy plain distinguished-point rho for discrete log (no automorphism).

This is the clean baseline attack: a standard R-adding distinguished-point rho on
a cyclic group of prime order, with van Oorschot-Wiener parallel collision search
and cheap incremental restarts.  It uses NO automorphism/orbit quotient, hence
has NO fruitless cycles -- the only assumption is the textbook random-mapping
model for the walk's collision time.

Group: the prime-order-l subgroup of F_q* (a cyclic group of order l), standing
in for GT (cyclic, prime order r).  Every group multiplication of the walk,
look-ahead (none here) and restarts is counted in 'mults'.

Reports mean mults / sqrt(pi*l/2), and the success fraction within the budget
1.25*k0, k0 = sqrt(2*l*ln(1/0.61)) (the random-mapping 0.39 point).  Toy scale.
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


def trial(l, q, g, rng, partitions, walks, theta_bits, counters):
    x = rng.randrange(1, l)
    h = pow(g, x, q)
    table = []
    for _ in range(partitions):
        a, b = rng.randrange(l), rng.randrange(l)
        table.append((pow(g, a, q) * pow(h, b, q) % q, a, b))
    mask = (1 << theta_bits) - 1
    seen = {}

    def fresh(prev):
        if prev is None:
            a, b = rng.randrange(l), rng.randrange(l)
            u = pow(g, a, q) * pow(h, b, q) % q
        else:
            pu, pa, pb = prev
            jj = rng.randrange(partitions)
            m, ma, mb = table[jj]
            u = pu * m % q                # one multiplication (counted)
            a, b = (pa + ma) % l, (pb + mb) % l
            counters['mults'] += 1
        return [u, a, b]

    state = [fresh(None) for _ in range(walks)]
    while True:
        for w in state:
            u, a, b = w
            j = (u * 0x9E3779B1 >> 7) % partitions
            m, ma, mb = table[j]
            nu, na, nb = u * m % q, (a + ma) % l, (b + mb) % l
            counters['mults'] += 1
            if (nu * 0x2545F491 >> 3) & mask == 0:        # distinguished point
                if nu in seen:
                    pa, pb = seen[nu]
                    if (pb - nb) % l:
                        xr = (na - pa) * pow(pb - nb, l - 2, l) % l
                        if pow(g, xr, q) == h:
                            return
                else:
                    seen[nu] = (na, nb)
                w[:] = fresh((nu, na, nb))
                continue
            w[:] = [nu, na, nb]


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 26))
    trials = int(args.get("trials", 120))
    seed = int(args.get("seed", 1))
    walks = int(args.get("walks", 8))
    theta_bits = int(args.get("theta_bits", 4))
    partitions = int(args.get("partitions", 256))
    rng = random.Random(seed)
    l, q, g = setup(bits, rng)
    budget = 1.25 * math.sqrt(2 * l * math.log(1 / 0.61))
    expected = math.sqrt(math.pi * l / 2)
    totals = []
    for _ in range(trials):
        c = {'mults': 0}
        trial(l, q, g, rng, partitions, walks, theta_bits, c)
        totals.append(c['mults'])
    ok = sum(1 for t in totals if t <= budget)
    print(json.dumps({
        "l": l, "trials": trials, "partitions_R": partitions, "walks": walks,
        "distinguished_fraction": 2.0 ** -theta_bits,
        "mean_mults": sum(totals) / trials,
        "mean_mults_over_sqrt_pi_l_over_2": sum(totals) / trials / expected,
        "budget_mults": budget,
        "success_within_budget": ok / trials,
        "random_mapping_success_at_budget": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
