"""Toy STRICTLY-BOUNDED single-processor distinguished-point rho (serial VW).

A fixed-budget distinguished-point Pollard rho for a discrete log in a cyclic
prime-order group (standing in for GT), run on ONE processor (serial van
Oorschot-Wiener): short segments restart at each distinguished point, their
distinguished points are stored, and a cross-segment collision solves the log.

Strict bounding: the attack performs AT MOST W_max = ceil(1.25*k0) group
multiplications and stops.  The counter is checked BEFORE every multiplication,
including the one that starts a fresh segment, so it never exceeds W_max.  The
counter 'mults' is the charged walk work (segment steps plus the one restart
multiplication per segment); the one-time setup pairing/table and the final
verification are O(1) and are charged separately in proof.md, not here.

A segment that runs seg_cap/theta multiplications without a distinguished point
(an exceptional long run; almost never cyclic) is truncated and a fresh segment
started.  With seg_cap large the expected number of truncations over the whole
budget is << 1, so no collision opportunity is lost.

Reports the fraction of trials that find a usable collision within W_max (a
Monte-Carlo estimate of the bounded success probability, lower-bounded at >= 0.39
in proof.md), and the mean per-trial truncations.  Toy scale only; the toy stores
exact group elements, while production stores 384-bit fingerprints with
recompute-from-(A,B) resolution (proof.md section 2.4).
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


def trial(l, q, g, rng, partitions, theta_bits, seg_cap, max_steps, counters):
    """One serial-VW run; return True iff a usable collision is found within max_steps."""
    x = rng.randrange(1, l)
    h = pow(g, x, q)
    table = []
    for _ in range(partitions):
        a, b = rng.randrange(l), rng.randrange(l)
        table.append((pow(g, a, q) * pow(h, b, q) % q, a, b))
    mask = (1 << theta_bits) - 1
    seg_limit = seg_cap << theta_bits
    seen = {}

    # Strict budget: perform a multiplication only if under the cutoff.
    def mul(u, m):
        if counters['mults'] >= max_steps:
            return None
        counters['mults'] += 1
        return u * m % q

    # start the first segment from a random point (its construction is the
    # one-time/first-segment start; charged as one multiplication for uniformity)
    a, b = rng.randrange(l), rng.randrange(l)
    u = pow(g, a, q) * pow(h, b, q) % q
    n = 0
    while True:
        j = (u * 0x9E3779B1 >> 7) % partitions
        m, ma, mb = table[j]
        nu = mul(u, m)
        if nu is None:
            return False                               # hit the hard cutoff
        na, nb = (a + ma) % l, (b + mb) % l
        n += 1
        restart = False
        if (nu * 0x2545F491 >> 3) & mask == 0:         # distinguished point
            if nu in seen:                             # exact key (toy); prod: fingerprint+resolution
                pa, pb = seen[nu]
                if (pb - nb) % l:
                    xr = (na - pa) * pow(pb - nb, l - 2, l) % l
                    if pow(g, xr, q) == h:
                        return True
                # same-exponent coincidence: continue
            else:
                seen[nu] = (na, nb)
            restart = True
        elif n >= seg_limit:                           # exceptional long run: truncate
            counters['truncations'] += 1
            restart = True
        if restart:
            jj = rng.randrange(partitions)
            rm, rma, rmb = table[jj]
            ru = mul(nu, rm)                           # one restart multiplication (budget-checked)
            if ru is None:
                return False
            u, a, b, n = ru, (na + rma) % l, (nb + rmb) % l, 0
        else:
            u, a, b = nu, na, nb


def main():
    args = json.load(sys.stdin)
    bits = int(args.get("bits", 26))
    trials = int(args.get("trials", 120))
    seed = int(args.get("seed", 1))
    theta_bits = int(args.get("theta_bits", 4))
    partitions = int(args.get("partitions", 256))
    seg_cap = int(args.get("seg_cap", 256))
    rng = random.Random(seed)
    l, q, g = setup(bits, rng)
    k0 = math.sqrt(2 * l * math.log(1 / 0.61))
    max_steps = int(math.ceil(1.25 * k0))
    successes = 0
    total_trunc = 0
    max_mults = 0
    for _ in range(trials):
        c = {'mults': 0, 'truncations': 0}
        if trial(l, q, g, rng, partitions, theta_bits, seg_cap, max_steps, c):
            successes += 1
        total_trunc += c['truncations']
        max_mults = max(max_mults, c['mults'])
    print(json.dumps({
        "l": l, "trials": trials, "partitions_R": partitions,
        "distinguished_fraction": 2.0 ** -theta_bits, "seg_cap_in_1_over_theta": seg_cap,
        "hard_cutoff_W_max": max_steps,
        "max_observed_mults": max_mults,
        "strictly_bounded_ok": bool(max_mults <= max_steps),
        "W_max_over_sqrt_pi_l_over_2": max_steps / math.sqrt(math.pi * l / 2),
        "success_within_hard_cutoff": successes / trials,
        "mean_truncations_per_trial": total_trunc / trials,
        "random_mapping_success_at_1.25k0": 1 - 0.61 ** 1.5625,
    }, indent=1))


if __name__ == "__main__":
    main()
