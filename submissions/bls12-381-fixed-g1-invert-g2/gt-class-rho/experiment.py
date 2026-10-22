"""Toy end-to-end run of the class-rho method: group Z_n (n prime, n = 1 mod 12) with the order-12 automorphism
u -> lam*u standing in for z -> z^p. Parallel walks, distinguished points, and the class-defined cycle rule: a walk
that revisits a class label within the last P steps has found a cycle of the class map; it takes the cycle member with the
smallest label (a function of the cycle only), discards the lap, and escapes from that member by squaring its canonical
element (u -> 2u in the toy, the analogue of point doubling). Control rules: an escape through an ordinary table step with
index j+1 (the earlier, flawed rule: the next step can undo it and the walker loops to the cap), and a walker-history escape.
Toy only; it validates the logic (class canonicalisation, exponent recovery, success probability, and that walks which
share a class stay identical afterwards, also through cycles), not production costs."""
import json
import math
import random
import sys

P = 24                                           # cycle window (production: 24)


def prime(n):
    return n > 1 and all(n % q for q in range(2, math.isqrt(n) + 1))


class Toy:
    def __init__(self, n, R, seed):
        self.n, self.R = n, R
        self.events = self.repeats = 0
        self.lam = next(pow(g, (n - 1) // 12, n) for g in range(2, n)
                        if pow(pow(g, (n - 1) // 12, n), 6, n) != 1 and pow(pow(g, (n - 1) // 12, n), 4, n) != 1)
        self.pw = [pow(self.lam, i, n) for i in range(12)]
        rnd = random.Random(seed)
        self.table = [rnd.randrange(1, n) for _ in range(R)]       # step j adds table[j]: u += table[j]

    def canon(self, a, b, x):
        n = self.n
        i = min(range(12), key=lambda s: self.pw[s] * ((a + b * x) % n) % n)
        a, b = self.pw[i] * a % n, self.pw[i] * b % n
        return (a + b * x) % n, a, b                                # label, exponents of the canonical element

    def step(self, st, x, shift):
        v, a, b = st
        j = ((v * 2654435761 >> 8) + shift) % self.R               # index is a function of the class label
        return self.canon((a + self.table[j]) % self.n, b, x)

    def walk(self, a, b, x, dp, cut, rule="class"):
        """Returns (path, executed_steps, status). path = logical class path (labels); status DP or CUT."""
        path = [self.canon(a, b, x)]
        used = 0
        hist = []                                                   # control rule only
        seen = set()                                                # cycle members this walk already escaped from
        while True:
            cur = path[-1]
            if rule in ("class", "jplus1"):
                new = self.step(cur, x, 0)
                used += 1
                if new[0] % dp == 0:
                    path.append(new)
                    return path, used, "DP"
                w = next((w for w in range(max(0, len(path) - P), len(path)) if path[w][0] == new[0]), None)
                if w is not None:                                   # revisit within P: cycle path[w:]
                    m = min(range(w, len(path)), key=lambda t: path[t][0])
                    del path[m + 1:]                                # discard the lap past the minimum
                    _, ea, eb = path[-1]
                    self.events += 1
                    self.repeats += path[-1][0] in seen             # same minimum again: the escape fell back into the cycle
                    seen.add(path[-1][0])
                    if rule == "class":
                        new = self.canon(2 * ea % self.n, 2 * eb % self.n, x)   # escape: square the minimum member
                    else:
                        new = self.step(path[-1], x, 1)             # control: ordinary table step j+1
                    used += 1
                    if new[0] % dp == 0:
                        path.append(new)
                        return path, used, "DP"
                path.append(new)
            else:                                                   # control: history of this walker's own labels
                new = self.step(cur, x, 1 if cur[0] in hist[-3:] else 0)
                hist = (hist + [cur[0]])[-4:]
                used += 1
                path.append(new)
                if new[0] % dp == 0:
                    return path, used, "DP"
            if used > cut:
                return path, used, "CUT"


def merge_test(n, R, rule, walkers=160, cut=900, dp=150, seed=7):
    """Walks from random starts on one table. For each pair, find the first class one shares with the other; the
    claim is that from then on the two logical class sequences coincide up to the earlier end."""
    toy = Toy(n, R, seed)
    rnd = random.Random(seed + 1)
    x = rnd.randrange(1, n)
    paths = []
    cyc = 0
    for _ in range(walkers):
        p, used, st = toy.walk(rnd.randrange(n), rnd.randrange(n), x, dp, cut, rule)
        paths.append([s[0] for s in p])
        cyc += used - (len(p) - 1)
    pos = [{l: t for t, l in enumerate(p)} for p in paths]
    merged = bad = 0
    for i in range(len(paths)):
        for j in range(i + 1, len(paths)):
            hit = next((t for t, l in enumerate(paths[i]) if l in pos[j]), None)
            if hit is None:
                continue
            merged += 1
            u = pos[j][paths[i][hit]]
            m = min(len(paths[i]) - hit, len(paths[j]) - u)
            if paths[i][hit:hit + m] != paths[j][u:u + m]:
                bad += 1
    return {"merged_pairs": merged, "diverged_pairs": bad, "wasted_steps_fraction": cyc / max(1, sum(map(len, paths)))}


def trap_test(n, R, rule, walkers=1500, cut=900, dp=10 ** 9, seed=3):
    """Direct test of the escape: walks of `cut` steps with no distinguished points on a group with many classes, so
    cycles of f (2-cycles at 1/(12R) per step) are frequent and every escape event is counted. A repeat is an escape
    from a member already escaped from in the same walk, i.e. the walker fell back into the cycle (the trap)."""
    toy = Toy(n, R, seed)
    rnd = random.Random(seed + 1)
    x = rnd.randrange(1, n)
    for _ in range(walkers):
        toy.walk(rnd.randrange(n), rnd.randrange(n), x, dp, cut, rule)
    return {"escape_events": toy.events, "escapes_that_fell_back": toy.repeats,
            "fall_back_rate": round(toy.repeats / max(1, toy.events), 4)}


def run(n, trials):
    toy0 = Toy(n, 64, 0)
    classes = (n - 1) // 12 + 1
    k = math.sqrt(-2 * math.log(0.61) * classes)
    theta, walks = 4.0 / math.sqrt(classes), 6
    dpm = int(1 / theta)
    ok = 0
    steps_used = 0
    for t in range(trials):
        rnd = random.Random(1000 + t)
        x = rnd.randrange(1, n)
        toy = Toy(n, 64, 5000 + t)
        dps = {}
        budget = 1.1 * k + walks * 1.0 / theta
        used = 0
        found = None
        while used < budget and found is None:
            for w in range(walks):
                path, u, st = toy.walk(rnd.randrange(n), rnd.randrange(n), x, dpm, 20 * dpm)
                used += u
                if st != "DP":
                    continue
                v, a, b = path[-1]
                if v in dps and dps[v][1] != b:
                    a2, b2 = dps[v]
                    found = (a - a2) * pow(b2 - b, -1, n) % n
                    break
                dps[v] = (a, b)
        ok += found == x
        steps_used += used
    return {"order": n, "classes": classes, "trials": trials, "recovered": ok, "success_rate": ok / trials,
            "mean_steps_over_sqrt_classes": steps_used / trials / math.sqrt(classes)}


if __name__ == "__main__":
    n = json.load(sys.stdin)["order"]
    n = next(m for m in range(n, 2 * n) if m % 12 == 1 and prime(m))
    res = run(n, 400)
    assert res["success_rate"] >= 0.39, res
    # cycle stress: table of 6 entries makes 2-cycles about 1.4% per step, so most segments meet several cycles
    res["cycle_stress_class_rule"] = merge_test(n, 6, "class")
    res["cycle_stress_walker_history_control"] = merge_test(n, 6, "history")
    res["cycle_stress_table_step_escape_control"] = merge_test(n, 6, "jplus1")
    assert res["cycle_stress_class_rule"]["diverged_pairs"] == 0, res
    big = next(m for m in range(10 ** 8, 2 * 10 ** 8) if m % 12 == 1 and prime(m))
    res["escape_fall_back_squaring"] = trap_test(big, 6, "class")
    res["escape_fall_back_table_step_control"] = trap_test(big, 6, "jplus1")
    assert res["escape_fall_back_squaring"]["escapes_that_fell_back"] == 0, res
    res["scope"] = "toy group only; does not validate production cost bounds"
    print(json.dumps(res))
