"""Toy end-to-end run of the negation-map distinguished-point rho of proof.md.

Group: y^2 = x^3 - 3x + 4 over F_65537 (a = -3, as for P-256), prime order n = 65851, so the
nonzero points split into N' = (n+1)/2 = 32926 classes {X, -X}.  The canonical representative of a
class is the one with the smaller y as an integer, exactly as proof.md picks it.

The walk is the production walk with the constants scaled down and the rule unchanged.  A table of
T_j = c_j P + d_j Q with fresh nonzero coefficients c_j, d_j uniform below n; the table index is
read from the top bits of the canonical x; a point is distinguished when the low 6 bits of its
canonical x are zero (theta = 2^-6, toy value); a step adds T_j and re-canonicalises.  Each walker
keeps signed coefficients (A, B) and a sign s with point = s (A P + B Q), reduced modulo n only at
a distinguished point.  Every 16 steps the walker stores the full 17-bit canonical x as its
checkpoint and every new point is compared with all of those bits (the toy analogue of the 256-bit
compare); a match is a return to that exact class, so the walker is on a cycle of the walk.  It then
walks the cycle once, keeps the member with the smallest canonical x with its complete state,
restores it and adds one escape point from an independent escape table, indexed by its own bits of
that member's x.  The escape point is a function of the cycle only.  Table index, escape index and
distinguished-point test use disjoint bit fields, as in proof.md.

Two table sizes are used, because the two questions need different regimes.

  * The rule tests (merge, escape re-entry) use 64 table entries, so that a two-cycle arises about
    every 128 steps and the cycle rule is exercised on almost every walk.
  * The collision-ratio measurement uses 1024 table entries.  With 64 entries the cycle overhead
    alone is about 2^-7.7 per step, four orders of magnitude above the production value of
    2^-16.87 (proof.md section 2), and it would dominate the measured ratio.  With 1024 entries it
    is below 2^-13 per step, still far above production, so the measured ratio is an upper bound on
    the production walk's overhead, not a measurement of it.

Reported, all at toy scope:
  recovery            end-to-end scalar recovery by distinguished-point collisions, k of n;
  merge               two walkers at the same class, one keeping its checkpoint history and one
                      with an empty checkpoint (possibly starting inside a cycle):
                        same distinguished point      k1 of n,
                        both looping in one cycle     k2 of n (the escape landed back in the basin
                                                      of the same cycle: the toy basin fraction is
                                                      about 2^-5, the production bound 2^-135.4),
                        genuinely diverged           k3 of n (a non-zero k3 would refute the
                                                      merge claim);
  escapes/reentries   escape events, and how many land back in the basin of the cycle just left;
  collision ratio     total walk steps across 8 walkers until two different walkers share a class,
                      divided by sqrt(pi N'/2), with a 95% interval for the mean ratio.

Nothing here measures production work, production walk quality or the production memory bound; the
toy exercises the rule and reports the one ratio it can measure, as an upper bound.

Usage: python3 experiment.py [--inputs JSON]   (the default inputs run in under two seconds)
"""
import json
import math
import random
import sys

PFIELD = 65537          # field modulus, prime
ACOEF = -3              # curve coefficient a
BCOEF = 4               # curve coefficient b
NORDER = 65851          # prime order of the group (cofactor 1)
NCLASS = (NORDER + 1) // 2
DPBITS = 6              # theta = 2^-6
CKPT = 16               # checkpoint period (toy value; the production period is 64)
MAXL = 2 * CKPT         # a detected cycle has at most CKPT classes
GENERATOR = None
CFG = {"entries": 64, "main_shift": 10, "esc_shift": 7, "esc_entries": 8}


def is_prime(v):
    if v < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if v % q == 0:
            return v == q
    d, s = v - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, v)
        if x in (1, v - 1):
            continue
        for _ in range(s - 1):
            x = x * x % v
            if x == v - 1:
                break
        else:
            return False
    return True


def padd(p1, p2):
    """Affine addition on y^2 = x^3 + a x + b over F_PFIELD; None is the identity."""
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    (x1, y1), (x2, y2) = p1, p2
    if x1 == x2:
        if (y1 + y2) % PFIELD == 0:
            return None
        slope = (3 * x1 * x1 + ACOEF) * pow(2 * y1, -1, PFIELD) % PFIELD
    else:
        slope = (y2 - y1) * pow((x2 - x1) % PFIELD, -1, PFIELD) % PFIELD
    x3 = (slope * slope - x1 - x2) % PFIELD
    return (x3, (slope * (x1 - x3) - y1) % PFIELD)


def pmul(k, p1):
    out = None
    while k:
        if k & 1:
            out = padd(out, p1)
        p1 = padd(p1, p1)
        k >>= 1
    return out


def canon(p1):
    """Canonical representative of the class {X, -X}: the one with the smaller y."""
    if p1 is None:
        return None, 1
    x, y = p1
    neg = (-y) % PFIELD
    if y < neg:
        return (x, y), 1
    return (x, neg), -1


def setup():
    """A generator of the prime-order group, found through a field square-root table."""
    roots = {}
    for y in range(1, PFIELD):
        roots.setdefault(y * y % PFIELD, y)
    for x in range(1, PFIELD):
        rhs = (x * x * x + ACOEF * x + BCOEF) % PFIELD
        if rhs in roots:
            pt = (x, roots[rhs])
            if pmul(NORDER, pt) is None:
                return pt
    raise SystemExit("no point of the stated prime order")


def make_table(rng, q, count):
    """T_j = [c_j]P + [d_j]Q with fresh nonzero coefficients and distinct pairs."""
    table, pairs = [], set()
    while len(table) < count:
        c, d = rng.randrange(1, NORDER), rng.randrange(1, NORDER)
        if (c, d) in pairs:
            continue
        pairs.add((c, d))
        table.append((padd(pmul(c, GENERATOR), pmul(d, q)), c, d))
    return table


class Walker(object):
    __slots__ = ("pt", "a", "b", "s", "cnt", "chk", "last_min", "escaped", "fell_back")

    def __init__(self, pt, a=0, b=0, s=1, cnt=0, chk=None):
        self.pt, self.a, self.b, self.s = pt, a, b, s
        self.cnt = cnt
        self.chk = pt[0] if chk is None else chk
        self.last_min = None
        self.escaped = False
        self.fell_back = False

    def copy(self):
        return Walker(self.pt, self.a, self.b, self.s, self.cnt, self.chk)

    def coeffs(self):
        return (self.s * self.a) % NORDER, (self.s * self.b) % NORDER


def step(w, table, stats):
    """One walk step. When X = -T_j the affine sum is the identity (the zero-denominator case of
    proof.md, here probability about 1/65851 per step): the step doubles X instead, a function of
    the class, like the slow path the proof budgets."""
    tj, cj, dj = table[(w.pt[0] >> CFG["main_shift"]) % CFG["entries"]]
    raw = padd(w.pt, tj)
    if raw is None:
        raw = padd(w.pt, w.pt)
        w.a *= 2
        w.b *= 2
        stats["zero_denominator_steps"] += 1
    else:
        w.a += w.s * cj
        w.b += w.s * dj
    nxt, sign = canon(raw)
    w.s *= sign
    w.pt = nxt


def cycle_event(w, table, esc_table, stats):
    """Walk the cycle once, restore its smallest-x member, add the escape point of that member."""
    start = w.pt
    best = (start[0], (start, w.a, w.b, w.s))
    laps = 0
    while True:
        step(w, table, stats)
        laps += 1
        if w.pt == start or laps > MAXL:
            break
        if w.pt[0] < best[0]:
            best = (w.pt[0], (w.pt, w.a, w.b, w.s))
    w.pt, w.a, w.b, w.s = best[1]
    te, ce, de = esc_table[(w.pt[0] >> CFG["esc_shift"]) % CFG["esc_entries"]]
    nxt, sign = canon(padd(w.pt, te))
    w.a += w.s * ce
    w.b += w.s * de
    w.s *= sign
    w.pt = nxt
    w.chk, w.cnt = w.pt[0], 0
    stats["escapes"] += 1
    if not w.escaped:
        w.escaped = True
        stats["escaping_walks"] += 1
    if w.last_min == best[0] and not w.fell_back:
        w.fell_back = True
        stats["reentries"] += 1
    w.last_min = best[0]


def advance(w, table, esc_table, stats):
    """One step, then the distinguished-point test, then the full-width cycle compare."""
    step(w, table, stats)
    if (w.pt[0] & ((1 << DPBITS) - 1)) == 0:
        return "dp"
    if w.pt[0] == w.chk:
        cycle_event(w, table, esc_table, stats)
        return "cycle"
    w.cnt += 1
    if w.cnt >= CKPT:
        w.chk, w.cnt = w.pt[0], 0
    return "walk"


def start_walker(rng, q):
    a0, b0 = rng.randrange(1, NORDER), rng.randrange(1, NORDER)
    pt, sign = canon(padd(pmul(a0, GENERATOR), pmul(b0, q)))
    return Walker(pt, a0 * sign, b0 * sign, 1)


def solve(r1, r2, q):
    """x from a1 P + b1 Q = e (a2 P + b2 Q), e = +-1, checked by [x]P = Q."""
    a1, b1 = r1
    a2, b2 = r2
    for e in (1, -1):
        den = (b1 - e * b2) % NORDER
        if den:
            x = (e * a2 - a1) * pow(den, -1, NORDER) % NORDER
            if pmul(x, GENERATOR) == q:
                return x
    return None


def recovery_trial(rng, stats):
    x = rng.randrange(1, NORDER)
    q = pmul(x, GENERATOR)
    table, esc = make_table(rng, q, CFG["entries"]), make_table(rng, q, CFG["esc_entries"])
    walkers = [start_walker(rng, q) for _ in range(8)]
    seen = {}
    for _ in range(400):
        for w in walkers:
            if advance(w, table, esc, stats) != "dp":
                continue
            got = w.coeffs()
            if w.pt[0] in seen:
                if seen[w.pt[0]] != got:
                    val = solve(seen[w.pt[0]], got, q)
                    if val is not None:
                        stats["recovered"] += 1
                        return True
            else:
                seen[w.pt[0]] = got
    return False


def same_orbit(w_a, w_b, table, esc, stats, limit=8 * MAXL):
    """True when the two walkers' orbits meet (both loop in one cycle)."""
    seen = set()
    for _ in range(limit):
        advance(w_a, table, esc, stats)
        seen.add(w_a.pt[0])
    for _ in range(limit):
        advance(w_b, table, esc, stats)
        if w_b.pt[0] in seen:
            return True
    return False


def merge_pair(rng, stats, q, table, esc):
    """Two walkers at one class: one keeps the checkpoint history, one starts with an empty one."""
    w1 = start_walker(rng, q)
    for _ in range(rng.randrange(20, 60)):
        advance(w1, table, esc, stats)
    wa, wb = w1.copy(), Walker(w1.pt)
    dpa = dpb = None
    for _ in range(400):
        if dpa is None and advance(wa, table, esc, stats) == "dp":
            dpa = wa.pt
        if dpb is None and advance(wb, table, esc, stats) == "dp":
            dpb = wb.pt
        if dpa is not None and dpb is not None:
            break
    if dpa is not None and dpb is not None:
        if dpa == dpb:
            stats["merge_same_dp"] += 1
            return "same_dp"
        stats["merge_diverged"] += 1
        return "diverged"
    if same_orbit(wa, wb, table, esc, stats):
        stats["merge_trapped"] += 1
        return "trapped"
    stats["merge_diverged"] += 1
    return "diverged"


def collision_steps(rng, stats, q, table, esc):
    """Total walk steps across 8 walkers until two different walkers share a class."""
    walkers = [start_walker(rng, q) for _ in range(8)]
    where, steps = {}, 0
    for _ in range(400):
        for i, w in enumerate(walkers):
            advance(w, table, esc, stats)
            steps += 1
            key = w.pt[0]
            if key in where and where[key] != i:
                return steps
            where[key] = i
    return None


def new_world(rng):
    x = rng.randrange(1, NORDER)
    q = pmul(x, GENERATOR)
    return q, make_table(rng, q, CFG["entries"]), make_table(rng, q, CFG["esc_entries"])


def run(inputs):
    rng = random.Random(inputs["seed"])
    trials = int(inputs["trials"])
    stats = {"escapes": 0, "reentries": 0, "escaping_walks": 0, "recovered": 0,
             "zero_denominator_steps": 0, "merge_same_dp": 0, "merge_trapped": 0,
             "merge_diverged": 0}
    rule_entries = CFG["entries"]

    for _ in range(trials):
        recovery_trial(rng, stats)

    CFG.update({"entries": rule_entries, "main_shift": 10, "esc_shift": 7, "esc_entries": 8})
    merge_tables, merge_pairs = 12, 4 * trials
    for i in range(merge_tables):
        q, table, esc = new_world(rng)
        for _ in range(merge_pairs // merge_tables + (1 if i < merge_pairs % merge_tables else 0)):
            merge_pair(rng, stats, q, table, esc)
    escape_worlds = trials + trials // 3
    for _ in range(escape_worlds):
        q, table, esc = new_world(rng)
        for w in [start_walker(rng, q) for _ in range(8)]:
            for _ in range(80):
                advance(w, table, esc, stats)

    CFG.update({"entries": 1024, "main_shift": 7, "esc_shift": 6, "esc_entries": 2})
    collision_tables, collision_samples = 4, 8 * trials
    ratios = []
    for i in range(collision_tables):
        q, table, esc = new_world(rng)
        draws = collision_samples // collision_tables + (
            1 if i < collision_samples % collision_tables else 0)
        for _ in range(draws):
            got = collision_steps(rng, stats, q, table, esc)
            if got:
                ratios.append(got / math.sqrt(math.pi * NCLASS / 2.0))
    mean = sum(ratios) / len(ratios)
    var = sum((r - mean) ** 2 for r in ratios) / (len(ratios) - 1)
    half = 1.96 * math.sqrt(var / len(ratios))
    return {
        "group": {"p": PFIELD, "a": ACOEF, "b": BCOEF, "order": NORDER, "classes": NCLASS,
                  "theta_bits": DPBITS, "checkpoint_period": CKPT, "checkpoint_bits": 17,
                  "generator": list(GENERATOR)},
        "recovery": {"recovered": stats["recovered"], "trials": trials,
                     "table_entries": rule_entries},
        "merge": {"same_dp": stats["merge_same_dp"], "both_in_one_cycle": stats["merge_trapped"],
                  "diverged": stats["merge_diverged"], "pairs": merge_pairs,
                  "tables": merge_tables, "table_entries": rule_entries},
        "escapes": {"events": stats["escapes"], "escaping_walks": stats["escaping_walks"],
                    "walkers_that_fell_back": stats["reentries"],
                    "zero_denominator_steps": stats["zero_denominator_steps"],
                    "fresh_tables": escape_worlds + trials + merge_tables},
        "collision": {"samples": len(ratios), "tables": collision_tables, "table_entries": 1024,
                      "reference_sqrt_pi_N_over_2": round(math.sqrt(math.pi * NCLASS / 2.0), 3),
                      "ratio_mean": round(mean, 4),
                      "ratio_ci95": [round(mean - half, 4), round(mean + half, 4)]},
    }


INPUTS = {"seed": 20261023, "trials": 30}


def main(argv):
    global GENERATOR
    inputs = dict(INPUTS)
    if "--inputs" in argv:
        inputs.update(json.loads(argv[argv.index("--inputs") + 1]))
    if not is_prime(PFIELD) or not is_prime(NORDER):
        raise SystemExit("toy parameters are not prime")
    GENERATOR = setup()
    print(json.dumps(run(inputs), indent=2, sort_keys=True))


if __name__ == "__main__":
    main(sys.argv[1:])
