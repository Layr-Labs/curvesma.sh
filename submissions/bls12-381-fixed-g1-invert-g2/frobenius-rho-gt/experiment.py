"""Toy check of parallel distinguished-point rho on orbits of an order-12 automorphism.

The group is the prime-order-l subgroup of F_q^*, with 12 | l-1.  The
automorphism phi(x) = x^lam, where lam has order 12 mod l, plays the role of the
p-power Frobenius on GT (on GT, Frobenius is x -> x^(p mod r), and p has
order 12 mod r).  The walk uses the twisted-state method from proof.md: the
state u is never canonicalised.  Instead k = argmin_k key(phi^k(u)) is found and
the next state is u * phi^-k(M_j), using a pre-twisted table T[j][k] = phi^-k(M_j).
Every step is counted, including steps of abandoned (stuck) walks and the tail
after the last distinguished point.  The output reports mean steps relative to
sqrt(pi*N/2), N = (l-1)/12 + 1, and the success rate within the claimed budget
1.25*sqrt(2*N*ln(1/0.61)).  This is evidence at toy scale only.
"""
import json
import math
import random
import sys

ORBIT = 12


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
        l = (rng.getrandbits(bits) | (1 << (bits - 1))) // 12 * 12 + 1
        if not is_prime(l):
            continue
        for c in range(2, 4000, 2):
            q = c * l + 1
            if not is_prime(q):
                continue
            g = 1
            while g == 1:
                g = pow(rng.randrange(2, q - 1), (q - 1) // l, q)
            while True:
                lam = pow(rng.randrange(2, l), (l - 1) // 12, l)
                if all(pow(lam, d, l) != 1 for d in (1, 2, 3, 4, 6)):
                    return l, q, g, lam


def trial(l, q, g, lam, rng, partitions, walks, theta_bits):
    x = rng.randrange(1, l)
    h = pow(g, x, q)
    lam_pows = [pow(lam, k, l) for k in range(ORBIT)]
    lam_inv = [pow(lam, (ORBIT - k) % ORBIT, l) for k in range(ORBIT)]
    table = []
    for _ in range(partitions):
        a, b = rng.randrange(l), rng.randrange(l)
        row = []
        for k in range(ORBIT):
            ak, bk = a * lam_inv[k] % l, b * lam_inv[k] % l
            row.append((pow(g, ak, q) * pow(h, bk, q) % q, ak, bk))
        table.append(row)
    mask = (1 << theta_bits) - 1
    stuck = 20 << theta_bits
    seen = {}
    steps = 0

    def canon(u):
        # key(phi^k(u)) = u^(lam^k); return the minimal key and its k.
        best, bk = None, 0
        for k in range(ORBIT):
            v = pow(u, lam_pows[k], q)
            if best is None or v < best:
                best, bk = v, k
        return best, bk

    def fresh():
        a, b = rng.randrange(l), rng.randrange(l)
        u = pow(g, a, q) * pow(h, b, q) % q
        key, k = canon(u)
        return [u, a, b, key, k, 0]

    state = [fresh() for _ in range(walks)]
    while True:
        for w in state:
            u, a, b, key, k, n = w
            j = (key * 0x9E3779B1 >> 7) % partitions
            # Wiener-Zuccherato look-ahead on the class-invariant index.
            for t in range(partitions):
                jj = (j + t) % partitions
                m, ma, mb = table[jj][k]
                nu, na, nb = u * m % q, (a + ma) % l, (b + mb) % l
                nkey, nk = canon(nu)
                if (nkey * 0x9E3779B1 >> 7) % partitions != jj:
                    break
            steps += 1
            n += 1
            if (nkey * 0x2545F491 >> 3) & mask == 0:
                # Store exponents of the canonical representative phi^nk(nu).
                ra, rb = na * lam_pows[nk] % l, nb * lam_pows[nk] % l
                if nkey in seen:
                    pa, pb = seen[nkey]
                    if (pb - rb) % l:
                        xr = (ra - pa) * pow(pb - rb, l - 2, l) % l
                        if pow(g, xr, q) == h:
                            return steps
                else:
                    seen[nkey] = (ra, rb)
                w[:] = fresh()
                continue
            if n > stuck:
                w[:] = fresh()
                continue
            w[:] = [nu, na, nb, nkey, nk, n]


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 30))
    trials = int(args.get("trials", 100))
    seed = int(args.get("seed", 1))
    walks = int(args.get("walks", 8))
    theta_bits = int(args.get("theta_bits", 4))
    rng = random.Random(seed)
    l, q, g, lam = setup(bits, rng)
    classes = (l - 1) // ORBIT + 1
    budget = 1.25 * math.sqrt(2 * classes * math.log(1 / 0.61))
    expected = math.sqrt(math.pi * classes / 2)
    results = [trial(l, q, g, lam, rng, 64, walks, theta_bits) for _ in range(trials)]
    ok = sum(1 for s in results if s <= budget)
    print(json.dumps({
        "l": l, "q": q, "lam": lam, "orbit": ORBIT, "classes": classes,
        "trials": trials, "walks": walks,
        "distinguished_fraction": 2.0 ** -theta_bits,
        "approx_dp_lag_steps": walks * 2 ** theta_bits,
        "mean_steps": sum(results) / trials,
        "mean_over_sqrt_pi_N_over_2": sum(results) / trials / expected,
        "budget_steps": budget,
        "success_within_budget": ok / trials,
        "random_mapping_success_at_k0": 0.39,
        "random_mapping_success_at_budget": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
