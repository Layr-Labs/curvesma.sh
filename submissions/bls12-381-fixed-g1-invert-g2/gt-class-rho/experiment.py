"""Toy run of the class-rho walk with the production walk's own rules, on Z_n (n prime, n = 1 mod 12) with the order-12
automorphism u -> lam*u standing in for z -> z^p (lam^6 = -1 as in GT).
Same rules as the claim. The canonical coordinate of a class is the minimum of the 12 coordinate values of its elements, which
come in pairs (v, M - v) as in GT (v_{i+6} = p - v_i); here M = 2^40 and the coordinate of an element is a fixed pseudorandom
40-bit function of it (the production heuristic is exactly that the coordinates of a GT element behave so; an additive toy
group has no such coordinate of its own). The distinguished-point test is "the low d bits of the canonical coordinate are zero"
and the table index is the next log2(R) bits (disjoint low bit fields). A walker keeps the last 24 labels; a revisit of a
label within 24 steps is a cycle of the class map: it takes the cycle member with the smallest label, discards the lap and
escapes by squaring that member's canonical element (u -> 2u in the toy, the analogue of point doubling). Trails are cut after
20 * 2^d steps. The additive table (u += table[j]) plays the role of the multiplicative table of GT.
Experiments: (1) steps to the first collision between two different trails, for fresh instances, fresh tables and fresh
start coins, as a ratio to the random-draw benchmark sqrt(pi N'/2) + 2/3 for N' classes, with a 95% interval, the fraction of
runs that collide within the production step budget, the exponent recovery at each collision and the measured cycle-event
rate against the 2-cycle rate 1/(12R); (2) the same with the smallest table Teske covers (R = 16); (3) controls for the
escape rule. Toy scale only: it does not validate production costs."""
import json
import math
import random
import sys

P = 24                                           # cycle window (production: 24)
MU = 1.20                                        # production extra step margin (proof.md)
BITS = 40                                        # coordinate width of the toy (production: 381 bits)
M = 1 << BITS


def prime(n):
    return n > 1 and all(n % q for q in range(2, math.isqrt(n) + 1))


class Group:
    """Z_n with the order-12 automorphism and a pseudorandom coordinate; the canonical-coordinate table is built once."""

    def __init__(self, n, tabulate):
        self.n = n
        lam = next(pow(g, (n - 1) // 12, n) for g in range(2, n)
                   if pow(pow(g, (n - 1) // 12, n), 6, n) != 1 and pow(pow(g, (n - 1) // 12, n), 4, n) != 1)
        self.pw = [pow(lam, i, n) for i in range(6)]       # the other six conjugates are the negatives
        rnd = random.Random(n)
        self.H = [rnd.getrandbits(BITS) for _ in range(n // 2 + 1)] if tabulate else None
        self.tab = None
        if tabulate:
            self.tab = [self.canon_u(u) for u in range(n)]

    def hv(self, s):
        if self.H is not None:
            return self.H[s]
        z = (s + 0x9E3779B97F4A7C15) & (2 ** 64 - 1)
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9 & (2 ** 64 - 1)
        z = (z ^ (z >> 27)) * 0x94D049BB133111EB & (2 ** 64 - 1)
        return (z ^ (z >> 31)) & (M - 1)

    def canon_u(self, u):
        """(label, m): the canonical coordinate of the class of u and the multiplier m with m*u the canonical element."""
        if self.tab is not None:
            return self.tab[u]
        n, best = self.n, None
        for q in self.pw:
            t = u * q % n
            s, flip = (t, False) if t <= n - t else (n - t, True)   # s: the element of the pair {t, -t} below n/2
            h = self.hv(s)
            if h <= M - h:                                         # the pair's smaller value sits on s
                lab, m = h, (n - q if flip else q) % n
            else:                                                  # ... or on -s
                lab, m = M - h, (q if flip else n - q) % n
            if best is None or lab < best[0]:
                best = (lab, m)
        return best


class Toy:
    def __init__(self, grp, R, d, ishift=None):
        self.g, self.n, self.R, self.d, self.mask = grp, grp.n, R, d, (1 << d) - 1
        self.ishift = d if ishift is None else ishift          # first bit of the index field (default: just above the DP bits)
        self.events = self.repeats = self.two_cycles = 0
        self.table = None

    def set_table(self, rnd):
        self.table = [rnd.randrange(1, self.n) for _ in range(self.R)]   # step j adds table[j]

    def canon(self, a, b, x):
        lab, m = self.g.canon_u((a + b * x) % self.n)
        return lab, m * a % self.n, m * b % self.n          # label, exponents of the canonical element

    def step(self, st, x, shift):
        v, a, b = st
        j = ((v >> self.ishift) + shift) % self.R           # index: bits of the label, disjoint from the DP bits
        return self.canon((a + self.table[j]) % self.n, b, x)

    def walk(self, a, b, x, cut, rule="class"):
        """Returns (path, executed_steps, status). path = logical class path; status DP or CUT."""
        path = [self.canon(a, b, x)]
        pos = {path[0][0]: 0}
        used = 0
        hist = []                                           # control rule only
        seen = set()                                        # cycle members this walk already escaped from
        mask, n = self.mask, self.n
        while True:
            cur = path[-1]
            if rule in ("class", "jplus1"):
                new = self.step(cur, x, 0)
                used += 1
                if new[0] & mask == 0:
                    path.append(new)
                    return path, used, "DP"
                w = pos.get(new[0])
                if w is not None and len(path) - w <= P:    # revisit within P: cycle path[w:]
                    m = min(range(w, len(path)), key=lambda t: path[t][0])
                    for t in range(m + 1, len(path)):
                        if pos.get(path[t][0]) == t:
                            del pos[path[t][0]]
                    self.events += 1
                    self.two_cycles += len(path) - w == 2
                    del path[m + 1:]                        # discard the lap past the minimum
                    _, ea, eb = path[-1]
                    self.repeats += path[-1][0] in seen     # same minimum again: the escape fell back into the cycle
                    seen.add(path[-1][0])
                    if rule == "class":
                        new = self.canon(2 * ea % n, 2 * eb % n, x)        # escape: square the minimum member
                    else:
                        new = self.step(path[-1], x, 1)                    # control: ordinary table step j+1
                    used += 1
                    if new[0] & mask == 0:
                        path.append(new)
                        return path, used, "DP"
                pos[new[0]] = len(path)
                path.append(new)
            else:                                           # control: history of this walker's own labels
                new = self.step(cur, x, 1 if cur[0] in hist[-3:] else 0)
                hist = (hist + [cur[0]])[-4:]
                used += 1
                path.append(new)
                if new[0] & mask == 0:
                    return path, used, "DP"
            if used > cut:
                return path, used, "CUT"


def collision_runs(grp, R, d, trials, seed):
    """Fresh instance, fresh table and fresh start coins per trial. Draws are the start classes and the retained steps of
    successive trails; a collision is a label already drawn by a different trail (a repeat within the same trail is a
    self-collision, counted and not used)."""
    n = grp.n
    classes = (n - 1) // 12 + 1
    bench = math.sqrt(math.pi * classes / 2) + 2 / 3
    budget = MU * math.sqrt(2 * math.log(1 / 0.61) * classes * (1 + 1.0 / R))
    cut = 20 * (1 << d)
    toy = Toy(grp, R, d)
    S, rec, within, selfc, logical, executed = [], 0, 0, 0, 0, 0
    for t in range(trials):
        rnd = random.Random(seed + t)
        x = rnd.randrange(1, n)
        toy.set_table(rnd)
        seen, steps, tid, found = {}, 0, 0, None
        while found is None:
            tid += 1
            path, used, st = toy.walk(rnd.randrange(n), rnd.randrange(n), x, cut)
            executed += used
            logical += len(path) - 1
            for lab, pa, pb in path:
                steps += 1
                hit = seen.get(lab)
                if hit is not None:
                    if hit[0] == tid:
                        selfc += 1
                        break
                    found = (steps, hit, (pa, pb))
                    break
                seen[lab] = (tid, pa, pb)
        S.append(found[0])
        (_, ha, hb), (pa, pb) = found[1], found[2]
        if hb != pb:
            rec += (ha - pa) * pow(pb - hb, -1, n) % n == x
        within += found[0] <= budget
    mean = sum(S) / trials
    sd = math.sqrt(sum((s - mean) ** 2 for s in S) / (trials - 1))
    half = 1.96 * sd / math.sqrt(trials) / bench
    p = within / trials
    ph = 1.96 * math.sqrt(p * (1 - p) / trials)
    return {"order": n, "classes": classes, "R": R, "dp_bits": d, "trials": trials,
            "mean_steps_ratio": round(mean / bench, 4),
            "ratio_ci95": [round(mean / bench - half, 4), round(mean / bench + half, 4)],
            "within_production_budget": round(p, 4), "within_ci95": [round(p - ph, 4), round(p + ph, 4)],
            "exponent_recovered": rec, "recovered_of": trials, "self_collisions": selfc,
            "executed_over_retained_minus_1": round(executed / logical - 1, 6),
            "two_cycles_per_executed_step": round(toy.two_cycles / executed, 6),
            "two_cycle_rate_1_over_12R": round(1 / (12 * R), 6),
            "all_cycle_events_per_executed_step": round(toy.events / executed, 6)}


def merge_test(grp, R, d, rule, walkers=90, cut=600, seed=7):
    """Walks from random starts on one table. For each pair, find the first class one shares with the other; the
    claim is that from then on the two logical class sequences coincide up to the earlier end."""
    toy = Toy(grp, R, d)
    rnd = random.Random(seed + 1)
    toy.set_table(rnd)
    x = rnd.randrange(1, grp.n)
    paths = []
    for _ in range(walkers):
        p, _, _ = toy.walk(rnd.randrange(grp.n), rnd.randrange(grp.n), x, cut, rule)
        labels = []
        for s in p:                                         # a trail that repeats a label beyond the window is trapped in
            if s[0] in labels:                              # a long cycle (random-mapping event); keep it up to there
                break
            labels.append(s[0])
        paths.append(labels)
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
    return {"merged_pairs": merged, "diverged_pairs": bad}


def trap_test(grp, R, rule, walkers=50, cut=250, seed=3):
    """Direct test of the escape: walks of `cut` steps with no distinguished points (d = 40 never fires) on a group with
    many classes, so cycles of f are frequent and every escape event is counted. A repeat is an escape from a member
    already escaped from in the same walk, i.e. the walker fell back into the cycle (the trap)."""
    toy = Toy(grp, R, BITS, ishift=0)
    rnd = random.Random(seed + 1)
    toy.set_table(rnd)
    x = rnd.randrange(1, grp.n)
    for _ in range(walkers):
        toy.walk(rnd.randrange(grp.n), rnd.randrange(grp.n), x, cut, rule)
    return {"escape_events": toy.events, "escapes_that_fell_back": toy.repeats,
            "fall_back_rate": round(toy.repeats / max(1, toy.events), 4)}


if __name__ == "__main__":
    inp = json.load(sys.stdin)
    n = next(m for m in range(inp["order"], 2 * inp["order"]) if m % 12 == 1 and prime(m))
    trials = inp["trials"]
    grp = Group(n, True)
    out = {"scope": "toy group only; does not validate production cost bounds"}
    for R, d, share in ((64, 2, 0.6), (16, 2, 0.4)):
        r = collision_runs(grp, R, d, max(50, int(trials * share)), 1000 * R)
        out["walk_R%d_d%d" % (R, d)] = r
        assert r["exponent_recovered"] == r["recovered_of"], r
        assert r["within_ci95"][0] >= 0.39, r
    out["cycle_stress_class_rule"] = merge_test(grp, 6, 6, "class")
    out["cycle_stress_walker_history_control"] = merge_test(grp, 6, 6, "history")
    out["cycle_stress_table_step_escape_control"] = merge_test(grp, 6, 6, "jplus1")
    assert out["cycle_stress_class_rule"]["diverged_pairs"] == 0, out
    big = Group(next(m for m in range(10 ** 8, 2 * 10 ** 8) if m % 12 == 1 and prime(m)), False)
    out["escape_fall_back_squaring"] = trap_test(big, 6, "class")
    out["escape_fall_back_table_step_control"] = trap_test(big, 6, "jplus1")
    assert out["escape_fall_back_squaring"]["escapes_that_fell_back"] == 0, out
    print(json.dumps(out))
