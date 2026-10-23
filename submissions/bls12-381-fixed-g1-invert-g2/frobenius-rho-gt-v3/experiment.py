"""Toy orbit-rho with Brent cycle detection; measures fruitless-cycle entries
UNCENSORED and the walk's multiplication count.

Group: prime-order-l subgroup of F_q* with an order-12 automorphism phi(x)=x^lam
(12 | l-1), standing in for the order-12 Frobenius on GT.  Twisted-state walk:
the state is never canonicalised; the class key is min_k key(phi^k(u)); the next
state is u*T[j][k*] with pre-twisted entries T[j][k]=phi^-k(M_j).

This version answers the reviewer's remaining points:
 * Fruitless cycles are found by BRENT cycle detection on the class key within a
   segment, so a stuck walk is abandoned within O(cycle length) steps -- NOT only
   after a 20/theta timeout, and NOT censored by another walk ending the trial.
   Every cycle entry is counted ('cycle_entries'), whether or not the trial ends.
 * The partition index is H(class key); it depends only on H being a good hash of
   the class, not on which canonical key is used, so the toy full-element key and
   the production one-coordinate key induce the same partition law (demonstrated
   by the look-ahead-OFF entry rate ~ 1/R below).
 * The counter 'walk_mults' counts exactly the group multiplications of the walk,
   look-ahead retries and restarts -- i.e. the WALK LENGTH in steps.  It does NOT
   include the toy's canonicalisation exponentiations (production replaces these
   with the 8-Fp-mult one-coordinate key, accounted in proof.md), nor one-time
   table/pairing setup or final verification.  Scope is walk-step work only.

Reports mean walk_mults / sqrt(pi*N/2), success within 1.25*k0, the measured
cycle-entry rate per step, and the wasted-step fraction.  Toy scale only.
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
    hard_cap = 50 << theta_bits            # backstop only; Brent should fire first
    seen = {}

    def canon(u):
        best, bk = None, 0
        for k in range(ORBIT):
            v = pow(u, lam_pows[k], q)
            if best is None or v < best:
                best, bk = v, k
        return best, bk

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
            counters['walk_mults'] += 1
        key, k = canon(u)
        # state: u,a,b,key,k,n, brent_saved_key, brent_next_pow, brent_since
        return [u, a, b, key, k, 0, key, 1, 0]

    state = [fresh(None) for _ in range(walks)]
    while True:
        for w in state:
            u, a, b, key, k, n, bsaved, bpow, bsince = w
            j = (key * 0x9E3779B1 >> 7) % partitions
            attempts = partitions if lookahead else 1
            for t in range(attempts):
                jj = (j + t) % partitions
                m, ma, mb = table[jj][k]
                nu, na, nb = u * m % q, (a + ma) % l, (b + mb) % l
                counters['walk_mults'] += 1
                nkey, nk = canon(nu)
                if (nkey * 0x9E3779B1 >> 7) % partitions != jj:
                    break
            n += 1
            if (nkey * 0x2545F491 >> 3) & mask == 0:        # distinguished -> record, restart
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
            # Brent cycle detection on the class key (uncensored fruitless-cycle test)
            bsince += 1
            if nkey == bsaved:                               # cycle: class repeated
                counters['cycle_entries'] += 1
                counters['wasted'] += n
                w[:] = fresh((nu, na, nb))
                continue
            if bsince == bpow:                               # advance Brent checkpoint
                bsaved, bpow, bsince = nkey, bpow * 2, 0
            if n > hard_cap:                                 # backstop
                counters['wasted'] += n
                w[:] = fresh((nu, na, nb))
                continue
            w[:] = [nu, na, nb, nkey, nk, n, bsaved, bpow, bsince]


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 24))
    trials = int(args.get("trials", 120))
    seed = int(args.get("seed", 41))
    walks = int(args.get("walks", 8))
    theta_bits = int(args.get("theta_bits", 4))
    partitions = int(args.get("partitions", 256))
    lookahead = bool(args.get("lookahead", True))
    rng = random.Random(seed)
    l, q, g, lam = setup(bits, rng)
    classes = (l - 1) // ORBIT + 1
    budget = 1.25 * math.sqrt(2 * classes * math.log(1 / 0.61))
    expected = math.sqrt(math.pi * classes / 2)
    totals, waste_f, entry_rate = [], [], []
    for _ in range(trials):
        c = {'walk_mults': 0, 'wasted': 0, 'cycle_entries': 0}
        trial(l, q, g, lam, rng, partitions, walks, theta_bits, c, lookahead)
        totals.append(c['walk_mults'])
        waste_f.append(c['wasted'] / max(1, c['walk_mults']))
        entry_rate.append(c['cycle_entries'] / max(1, c['walk_mults']))
    ok = sum(1 for t in totals if t <= budget)
    print(json.dumps({
        "l": l, "classes": classes, "trials": trials,
        "partitions_R": partitions, "walks": walks, "lookahead": lookahead,
        "distinguished_fraction": 2.0 ** -theta_bits,
        "mean_walk_mults": sum(totals) / trials,
        "mean_walk_mults_over_sqrt_pi_N_over_2": sum(totals) / trials / expected,
        "budget_walk_mults": budget,
        "success_within_budget": ok / trials,
        "mean_wasted_fraction": sum(waste_f) / trials,
        "cycle_entries_per_step": sum(entry_rate) / trials,
        "one_over_R": 1.0 / partitions,
        "random_mapping_success_at_budget": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
