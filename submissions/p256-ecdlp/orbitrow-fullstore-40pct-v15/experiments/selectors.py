"""v15 selector arm: the TWO-step history test the production schedule actually keys.

v14 compared the successor class with a variable named hist0 that the update order
had already set equal to the CURRENT state, so it measured a one-step return.  Here
the successor S_{t+1} is compared with S_{t-1}: a return after two steps, which is
the trap the 1/(2R) schedule term defends.  Escape enabled, paired seeds, cluster
bootstrap over walks.  Exploratory only: no claim row depends on it.
"""
import json, os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v15kit as K

R = 64; WALKS = 160; STEPS = 220
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "evidence")


def table(P, Q):
    return [t[2] for t in K.row_table(P, Q, R)]


def walk(P, Q, tab, seed, use_ec, escape=True):
    rnd = random.Random(seed)
    S = K.add(K.mul(rnd.randrange(1, K.P_TOY), P), K.mul(rnd.randrange(1, K.P_TOY), Q))
    prev2 = prev1 = None
    traps = steps = 0
    for t in range(STEPS):
        r = ((S[0] ^ (t * 2654435761)) % R) if (use_ec and S) else rnd.randrange(R)
        C = K.canon(K.add(S, tab[r]))
        steps += 1
        if prev2 is not None and K.cls(C) == K.cls(prev2):   # S_{t+1} == S_{t-1}
            traps += 1
            if escape:
                S = K.add(C, tab[(r + 17) % R])
                prev2, prev1 = prev1, S
                continue
        prev2, prev1 = prev1, C
        S = C
    return traps, steps


def main():
    P = K.base_point(); Q = K.mul(424242, P); tab = table(P, Q)
    per = []
    te = th = 0
    for w in range(WALKS):
        seed = 90000 + w
        e, st = walk(P, Q, tab, seed, True)
        h, _ = walk(P, Q, tab, seed, False)
        te += e; th += h
        per.append((e - h) / st)
    rnd = random.Random(7); n = len(per); means = []
    for _ in range(2000):
        means.append(sum(per[rnd.randrange(n)] for _ in range(n)) / n)
    means.sort()
    lo, hi = means[50], means[1949]
    ev = {"design": {"walks": WALKS, "steps": STEPS, "history_depth":
                     "successor vs S_{t-1} (two-step return); v14 compared with the current state",
                     "escape": "enabled", "cluster": "walk", "reps": 2000},
          "observed": {"ec_traps": te, "hash_traps": th},
          "delta_per_step": sum(per) / n, "delta_ci95": [lo, hi],
          "ci_contains_zero": lo <= 0 <= hi,
          "verdict": "exploratory; selector transfer stays OPEN and no claim row depends on it"}
    json.dump(ev, open(os.path.join(OUT, "selector_arm_v15.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"delta": ev["delta_per_step"], "ci": ev["delta_ci95"], "traps": [te, th]}, indent=1))


if __name__ == "__main__":
    main()
