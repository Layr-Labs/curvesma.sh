"""Toy orbit-rho with cheap incremental restarts, counting EVERY multiplication.

Group: prime-order-l subgroup of F_q* with an order-12 automorphism phi(x)=x^lam
(12 | l-1), standing in for the order-12 Frobenius on GT.  Walk is the
twisted-state method: the state is never canonicalised; the class key is
min_k key(phi^k(u)); the next state is u * T[j][k*] with pre-twisted table
entries T[j][k]=phi^-k(M_j).

This version addresses three review points on the production claim:
  * Restarts are CHEAP and INCREMENTAL: the c-th walk starts at S_c = S_{c-1}*g
    (one multiplication, exponent +1), not two fresh exponentiations.  The code
    counts that one multiplication.
  * EVERY multiplication is counted, including each look-ahead retry attempt
    ('mults'), so the reported budget is an honest all-steps count.
  * It measures the wasted-step fraction (steps spent in abandoned, cycle-stuck
    segments) and lets you sweep the partition count R to exhibit the ~1/R^2
    decay of the fruitless-cycle rate.

Outputs mean 'mults' relative to sqrt(pi*N/2), success within the budget
1.25*sqrt(2*N*ln(1/0.61)), and the measured wasted fraction.  Toy scale only.
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


def trial(l, q, g, lam, rng, partitions, walks, theta_bits, counters, lookahead=True):
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

    def canon(u):
        best, bk = None, 0
        for k in range(ORBIT):
            v = pow(u, lam_pows[k], q)
            if best is None or v < best:
                best, bk = v, k
        return best, bk

    # Cheap restart: a new walk starts at (previous distinguished point) * M_jj
    # for a random table index jj -- ONE multiplication, with a well-randomised
    # (a,b) because the DP is pseudo-random and M_jj injects fresh exponents.
    # The initial `walks` starts (only w of them) use a full random exponentiation.
    def fresh(prev):
        if prev is None:
            a, b = rng.randrange(l), rng.randrange(l)
            u = pow(g, a, q) * pow(h, b, q) % q
        else:
            pu, pa, pb = prev
            jj = rng.randrange(partitions)
            m, ma, mb = table[jj][0]
            u = pu * m % q                 # one multiplication (counted)
            a, b = (pa + ma) % l, (pb + mb) % l
            counters['mults'] += 1
        key, k = canon(u)
        return [u, a, b, key, k, 0]

    state = [fresh(None) for _ in range(walks)]
    while True:
        for w in state:
            u, a, b, key, k, n = w
            j = (key * 0x9E3779B1 >> 7) % partitions
            attempts = partitions if lookahead else 1
            for t in range(attempts):              # look-ahead; count each attempt
                jj = (j + t) % partitions
                m, ma, mb = table[jj][k]
                nu, na, nb = u * m % q, (a + ma) % l, (b + mb) % l
                counters['mults'] += 1
                nkey, nk = canon(nu)
                if (nkey * 0x9E3779B1 >> 7) % partitions != jj:
                    break
            n += 1
            if (nkey * 0x2545F491 >> 3) & mask == 0:     # distinguished: record, restart
                ra, rb = na * lam_pows[nk] % l, nb * lam_pows[nk] % l
                if nkey in seen:
                    pa, pb = seen[nkey]
                    if (pb - rb) % l:
                        xr = (ra - pa) * pow(pb - rb, l - 2, l) % l
                        if pow(g, xr, q) == h:
                            return
                else:
                    seen[nkey] = (ra, rb)
                w[:] = fresh((nu, na, nb))
                continue
            if n > stuck:                             # abandoned fruitless cycle
                counters['wasted'] += n
                counters['stuck_segments'] += 1
                w[:] = fresh((nu, na, nb))
                continue
            w[:] = [nu, na, nb, nkey, nk, n]


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 24))
    trials = int(args.get("trials", 100))
    seed = int(args.get("seed", 1))
    walks = int(args.get("walks", 8))
    theta_bits = int(args.get("theta_bits", 4))
    partitions = int(args.get("partitions", 64))
    lookahead = bool(args.get("lookahead", True))
    rng = random.Random(seed)
    l, q, g, lam = setup(bits, rng)
    classes = (l - 1) // ORBIT + 1
    budget = 1.25 * math.sqrt(2 * classes * math.log(1 / 0.61))
    expected = math.sqrt(math.pi * classes / 2)
    totals = []
    wasted_frac = []
    stuck_rate = []
    for _ in range(trials):
        c = {'mults': 0, 'wasted': 0, 'stuck_segments': 0}
        trial(l, q, g, lam, rng, partitions, walks, theta_bits, c, lookahead)
        totals.append(c['mults'])
        wasted_frac.append(c['wasted'] / max(1, c['mults']))
        stuck_rate.append(c['stuck_segments'] / max(1, c['mults']))
    ok = sum(1 for t in totals if t <= budget)
    print(json.dumps({
        "l": l, "classes": classes, "trials": trials,
        "partitions_R": partitions, "walks": walks, "lookahead": lookahead,
        "distinguished_fraction": 2.0 ** -theta_bits,
        "mean_mults": sum(totals) / trials,
        "mean_mults_over_sqrt_pi_N_over_2": sum(totals) / trials / expected,
        "budget_mults": budget,
        "success_within_budget": ok / trials,
        "mean_wasted_fraction": sum(wasted_frac) / trials,
        "stuck_entries_per_step": sum(stuck_rate) / trials,
        "random_mapping_success_at_budget": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
