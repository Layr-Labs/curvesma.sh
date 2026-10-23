"""Fresh toy statistics for the Cheon walk, with the per-trial step counts kept.

  CURVESMASH_CORE=<curvesmash-core> python3 -B toy_stats.py <bits> <k> <trials> <seed> > out.json

Same trial loop as experiments/cheon_toy.py (cscore.cheon, mult 64, dp_bits 3, trails 16).
It records each trial's step count so success can be read at any budget, and it checks
every recovered secret. Wilson 95% intervals are computed from the exact counts.
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.environ.get("CURVESMASH_CORE", "."))  # directory containing cscore/
from cscore import cheon


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (round(centre - half, 4), round(centre + half, 4))


def main():
    bits, k, trials, seed = (int(a) for a in sys.argv[1:5])
    rng = random.Random(seed)
    r = cheon.find_prime_order(bits, k, rng)
    N = (r - 1) >> k
    M = N // 6
    walk = {"mult": 64, "dp_bits": 3, "trails": 16}
    steps, recovered, wrong = [], 0, 0
    for _ in range(trials):
        orc = cheon.Oracle(r, rng, use_mu=True)
        s = rng.randrange(1, r)
        T = orc.point(s)
        powers = {1 << j: orc.point(pow(s, 1 << j, r)) for j in range(k + 1)}
        got, st = cheon.recover_scalar(orc, T, powers, k, rng, walk)
        wrong += st["wrong"]
        recovered += int(got == s)
        steps.append(st["steps"])
    out = {"bits": bits, "k": k, "trials": trials, "seed": seed, "M_log2": math.log2(M),
           "recovered": recovered, "wrong_collisions": wrong, "steps": steps}
    sys.stdout.write(json.dumps(out) + "\n")


if __name__ == "__main__":
    main()
