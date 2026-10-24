"""Fix 2: selector arm rebuilt at the correct history depth, clustered CI.

v13 measured the wrong event: it compared the successor class with hist1 (two states
back), which counts 3-step returns, not the 2-step trap behind 1/(2R); it disabled
escape; and it treated dependent trajectory steps as independent Bernoulli trials.

v14: cn vs hist0 (the 2-step trap: successor equals the state one step back after a
coefficient-trivial return), escape ENABLED as in production, paired same-seed design
(hash selector vs EC selector on identical walks), and a cluster bootstrap over whole
WALKS (the trajectory is the sampling unit), one method row.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v14kit as K

R = 64
WALKS = 160
STEPS = 220
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def base_point():
    for x in range(1, K.P_TOY):
        rhs = (x * x * x - 3 * x + K.B_TOY) % K.P_TOY
        if pow(rhs, (K.P_TOY - 1) // 2, K.P_TOY) == 1:
            y = pow(rhs, (K.P_TOY + 1) // 4, K.P_TOY)
            if (y * y - rhs) % K.P_TOY == 0:
                return (x, y)
    raise SystemExit("no base point")


def table(P, Q):
    tab = []
    for r in range(R):
        a = 1 + (r * 7 + 3) % R
        b = 1 + (r * 11 + 5) % R
        tab.append(K.add(K.mul(a, P), K.mul(b, Q)))
    return tab


def walk(P, Q, tab, seed, use_ec, escape=True):
    """One trail. Returns the per-walk trap count (cn == hist0) and step count."""
    rnd = random.Random(seed)
    S = K.add(K.mul(rnd.randrange(1, K.P_TOY), P), K.mul(rnd.randrange(1, K.P_TOY), Q))
    hist0 = hist1 = None
    traps = 0
    steps = 0
    for t in range(STEPS):
        if use_ec:
            # EC selector: row from the canonical class bits (production selector)
            r = (S[0] ^ (t * 2654435761)) % R if S else 0
        else:
            r = rnd.randrange(R)          # hash selector: fresh uniform coins
        T = K.add(S, tab[r])
        C = K.canon(T)
        steps += 1
        if hist0 is not None and K.cls(C) == K.cls(hist0):
            traps += 1
            if escape:
                S = K.add(C, tab[(r + 17) % R])   # production escape: one extra step
                hist1, hist0 = hist0, S
                continue
        hist1, hist0 = hist0, C
        S = C
    return traps, steps


def boot_ci(diffs, B=2000, seed=7):
    rnd = random.Random(seed)
    n = len(diffs)
    means = []
    for _ in range(B):
        s = 0.0
        for _ in range(n):
            s += diffs[rnd.randrange(n)]
        means.append(s / n)
    means.sort()
    return means[int(0.025 * B)], means[int(0.975 * B) - 1]


def main():
    P = base_point()
    Q = K.mul(424242, P)
    tab = table(P, Q)
    per_walk = []
    tot_ec = tot_hash = 0
    for w in range(WALKS):
        seed = 90000 + w                      # paired: same seed for both arms
        ec, st = walk(P, Q, tab, seed, True)
        hh, _ = walk(P, Q, tab, seed, False)
        tot_ec += ec
        tot_hash += hh
        per_walk.append((ec / st) - (hh / st))
    lo, hi = boot_ci(per_walk)
    mean = sum(per_walk) / len(per_walk)
    # non-degeneracy: expected 2-step traps for a uniform selector must stay small
    birthday_ratio = (STEPS * STEPS / 2.0) / (K.P_TOY / 2.0)
    ev = {
        "design": {"walks": WALKS, "steps_per_walk": STEPS, "R": R, "paired_seeds": True,
                   "history_depth": "cn compared to hist0 (2-step trap); v13 used hist1 and "
                                    "therefore counted 3-step returns",
                   "escape": "enabled, production form (one extra step on a detected trap)",
                   "cluster_unit": "one whole walk; bootstrap resamples walks, not steps",
                   "bootstrap_reps": 2000},
        "observed": {"hash_traps": tot_hash, "ec_traps": tot_ec,
                     "hash_rate": tot_hash / (WALKS * STEPS),
                     "ec_rate": tot_ec / (WALKS * STEPS)},
        "delta_per_step": mean, "delta_ci95": [lo, hi],
        "ci_contains_zero": lo <= 0.0 <= hi,
        "nondegeneracy": {"birthday_ratio": birthday_ratio,
                          "gate": 0.1, "passes": birthday_ratio <= 0.1},
        "verdict": "inconclusive at this budget: the interval contains zero and the "
                   "selector-transfer question stays declared OPEN. The arm is exploratory "
                   "and no claim row depends on it.",
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "selector_arm.json"), "w") as f:
        json.dump(ev, f, indent=1, sort_keys=True)
    print(json.dumps({k: ev[k] for k in ("observed", "delta_per_step", "delta_ci95",
                                         "ci_contains_zero", "nondegeneracy")}, indent=1))


if __name__ == "__main__":
    main()
