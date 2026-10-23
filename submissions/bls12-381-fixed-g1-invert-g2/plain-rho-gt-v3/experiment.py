"""Toy BOUNDED distinguished-point rho for discrete log (hard work cutoff).

A fixed-budget distinguished-point Pollard rho on a cyclic prime-order group
(standing in for GT).  No automorphism, so the only cycles are rho tails; a
trajectory may rarely enter a cycle containing no distinguished point, which the
bounded implementation handles explicitly:

  * HARD GLOBAL CUTOFF: the attack performs at most W_max = ceil(1.25*k0) total
    group multiplications (summed over all walks) and then stops, success or not.
    This W_max is the charged work; every restart, truncated and exceptional-path
    multiplication is inside it (all counted in 'mults', which never exceeds
    W_max).
  * PER-WALK TRUNCATION: a walk that runs seg_cap/theta steps without producing a
    distinguished point (an exceptional distinguished-point-free cycle) is
    abandoned and restarted.  A healthy walk produces a distinguished point in
    ~1/theta steps, so it is truncated only with probability ~e^{-seg_cap},
    negligible; the rare truncated walks waste at most seg_cap/theta steps, all
    inside W_max.

Success is "collision found within W_max": a Monte-Carlo estimate of the bounded
algorithm's success probability, which the proof lower-bounds at >= 0.39.  Toy
scale only.
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


def trial(l, q, g, rng, partitions, walks, theta_bits, seg_cap, max_steps, counters):
    """Return True iff a usable collision is found within max_steps mults."""
    x = rng.randrange(1, l)
    h = pow(g, x, q)
    table = []
    for _ in range(partitions):
        a, b = rng.randrange(l), rng.randrange(l)
        table.append((pow(g, a, q) * pow(h, b, q) % q, a, b))
    mask = (1 << theta_bits) - 1
    seg_limit = seg_cap << theta_bits           # per-walk truncation length
    seen = {}

    def fresh(prev):
        if prev is None:
            a, b = rng.randrange(l), rng.randrange(l)
            u = pow(g, a, q) * pow(h, b, q) % q
        else:
            pu, pa, pb = prev
            jj = rng.randrange(partitions)
            m, ma, mb = table[jj]
            u = pu * m % q                      # one multiplication (counted)
            a, b = (pa + ma) % l, (pb + mb) % l
            counters['mults'] += 1
        return [u, a, b, 0]                      # u, A, B, steps-in-segment

    state = [fresh(None) for _ in range(walks)]
    while counters['mults'] < max_steps:
        for w in state:
            if counters['mults'] >= max_steps:
                return False
            u, a, b, n = w
            j = (u * 0x9E3779B1 >> 7) % partitions
            m, ma, mb = table[j]
            nu, na, nb = u * m % q, (a + ma) % l, (b + mb) % l
            counters['mults'] += 1
            n += 1
            if (nu * 0x2545F491 >> 3) & mask == 0:        # distinguished point
                if nu in seen:
                    pa, pb = seen[nu]
                    if (pb - nb) % l:
                        xr = (na - pa) * pow(pb - nb, l - 2, l) % l
                        if pow(g, xr, q) == h:
                            return True
                else:
                    seen[nu] = (na, nb)
                w[:] = fresh((nu, na, nb))
                continue
            if n >= seg_limit:                            # exceptional DP-free cycle: truncate
                counters['truncations'] += 1
                w[:] = fresh((nu, na, nb))
                continue
            w[:] = [nu, na, nb, n]
    return False


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 26))
    trials = int(args.get("trials", 120))
    seed = int(args.get("seed", 1))
    walks = int(args.get("walks", 8))
    theta_bits = int(args.get("theta_bits", 4))
    partitions = int(args.get("partitions", 256))
    seg_cap = int(args.get("seg_cap", 64))
    rng = random.Random(seed)
    l, q, g = setup(bits, rng)
    k0 = math.sqrt(2 * l * math.log(1 / 0.61))
    max_steps = int(math.ceil(1.25 * k0))                 # hard global charged-work cutoff
    successes = 0
    total_trunc = 0
    for _ in range(trials):
        c = {'mults': 0, 'truncations': 0}
        if trial(l, q, g, rng, partitions, walks, theta_bits, seg_cap, max_steps, c):
            successes += 1
        total_trunc += c['truncations']
    print(json.dumps({
        "l": l, "trials": trials, "partitions_R": partitions, "walks": walks,
        "distinguished_fraction": 2.0 ** -theta_bits,
        "seg_cap_in_1_over_theta": seg_cap,
        "hard_cutoff_W_max": max_steps,
        "W_max_over_sqrt_pi_l_over_2": max_steps / math.sqrt(math.pi * l / 2),
        "success_within_hard_cutoff": successes / trials,
        "mean_truncations_per_trial": total_trunc / trials,
        "random_mapping_success_at_1.25k0": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
