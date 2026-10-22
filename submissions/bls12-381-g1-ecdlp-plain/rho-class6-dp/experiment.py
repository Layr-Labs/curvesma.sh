"""Toy end-to-end run of the proof's method: parallel rho with distinguished points on classes under
negation and the order-3 endomorphism, r-adding walk, canonical class representative, fixed-period cycle detection with a cycle-determined escape (double the smallest class of the cycle).
Toy curve y^2=x^3+2 over F_30091 with prime order 29881. Not evidence of production costs."""
import json
import math
import random
import sys

p, b, n = 30091, 2, 29881
beta = next(pow(g, (p - 1) // 3, p) for g in range(2, p) if pow(g, (p - 1) // 3, p) != 1)


def add(A, B):
    if A is None:
        return B
    if B is None:
        return A
    (x1, y1), (x2, y2) = A, B
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def mul(k, A):
    R = None
    while k:
        if k & 1:
            R = add(R, A)
        A = add(A, A)
        k >>= 1
    return R


def point():
    for x in range(1, p):
        v = (x ** 3 + b) % p
        if pow(v, (p - 1) // 2, p) == 1:
            for y in range(p):
                if y * y % p == v:
                    return (x, y)


P = point()
assert mul(n, P) is None
LAM = next(l for l in range(2, n) if (l * l + l + 1) % n == 0 and mul(l, P) == (beta * P[0] % p, P[1]))
MU = [(sg * pow(LAM, k, n)) % n for k in range(3) for sg in (1, -1)]  # scalar of each automorphism


def image(Q, s):  # apply automorphism index s=(k,sign)
    k, sg = divmod(s, 2)
    return (pow(beta, k, p) * Q[0] % p, Q[1] if sg == 0 else (-Q[1]) % p)


def canon(Q):  # smallest x of the three images, then y < p/2; returns (class point, index s with s.Q = class point)
    best = min(range(3), key=lambda k: pow(beta, k, p) * Q[0] % p)
    s = 2 * best
    C = image(Q, s)
    if 2 * C[1] > p:
        s += 1
        C = image(Q, s)
    return C, s


ALPHA = 16  # toy checkpoint period (production: 64)
INV = {}
for _s in range(6):
    _k, _sg = divmod(_s, 2)
    INV[_s] = (-_k % 3) * 2 + _sg


def make_table(rng, T, tab):
    return [(a, c, add(mul(a, P), mul(c, T))) for a, c in ((rng.randrange(n), rng.randrange(n)) for _ in range(tab))]


def step(state, entries):
    """One class step. state=(R,a,c). Returns (class point C, new state). The next class depends only on the class."""
    R, a, c = state
    C, s = canon(R)
    h = (C[0] * 2654435761 >> 5) % len(entries)
    ea, ec, ET = entries[h]
    i = INV[s]
    return C, (add(R, image(ET, i)), (a + MU[i] * ea) % n, (c + MU[i] * ec) % n)


def escape(state, entries):
    """state is on a cycle (its class just returned to the checkpoint). Walk the cycle once, find the member with the
    smallest class x (a function of the cycle), restore the state held there and double it."""
    C0, _ = canon(state[0])
    best, cur, L = None, state, 0
    while True:
        C, nxt = step(cur, entries)
        if best is None or C[0] < best[0][0]:
            best = (C, cur)
        cur, L = nxt, L + 1
        if canon(cur[0])[0] == C0:
            break
    R, a, c = best[1]
    return (add(R, R), 2 * a % n, 2 * c % n), L, best[0], canon(add(R, R))[0]


def solve(x, rng, tab=64, theta_bits=3, walkers=8, stats=None):
    T = mul(x, P)
    if T is None:
        return 0, 0
    entries = make_table(rng, T, tab)
    seen, steps = {}, 0

    def fresh():
        a, c = rng.randrange(n), rng.randrange(n)
        return [(add(mul(a, P), mul(c, T)), a, c), None, 0]  # state, saved checkpoint S, step counter t

    ws = [fresh() for _ in range(walkers)]
    while True:
        for w in ws:
            state = w[0]
            if state[0] is None:
                w[:] = fresh()
                continue
            C, nxt = step(state, entries)
            steps += 1
            if C[0] % (1 << theta_bits) == 0:  # distinguished class (checked first)
                _, s = canon(state[0])
                ka, kc = MU[s] * state[1] % n, MU[s] * state[2] % n
                if C in seen and seen[C] != (ka, kc):
                    a2, c2 = seen[C]
                    if (c2 - kc) % n:
                        return (ka - a2) * pow(c2 - kc, -1, n) % n, steps
                seen[C] = (ka, kc)
                w[:] = fresh()
                continue
            w[2] += 1
            if w[1] == C:  # (i) back at the checkpoint: C is on a cycle of the class walk
                esc, L, member, cls_after = escape(state, entries)
                if stats is not None:
                    stats["events"] += 1
                    stats["cycle_len"][L] = stats["cycle_len"].get(L, 0) + 1
                    assert L <= ALPHA
                    steps += 2 * L  # cycle walk (and a second pass for the escape check) are real work in the toy count
                w[:] = [esc, None, 0]
                continue
            if w[2] % ALPHA == 0:  # (ii) compare first, then save
                w[1] = C
            w[0] = nxt


def cycle_test(rng, trials=200, tab=64, horizon=3000):
    """Walker-independence of the escape: for every detected cycle, the escape computed from each member of the cycle
    (each possible detection point) gives the same class, namely the double of the smallest member."""
    T = mul(rng.randrange(1, n), P)
    entries = make_table(rng, T, tab)
    events = checked = 0
    for _ in range(trials):
        a, c = rng.randrange(n), rng.randrange(n)
        state = (add(mul(a, P), mul(c, T)), a, c)
        saved, t, seen_cls = None, 0, []
        for _ in range(horizon):
            if state[0] is None:
                break
            C, nxt = step(state, entries)
            t += 1
            if saved == C:
                L = t - saved_t
                events += 1
                ends, cur = set(), state
                for _ in range(L):  # detect from every member of the cycle
                    ends.add(escape(cur, entries)[3])
                    cur = step(cur, entries)[1]
                assert len(ends) == 1, "escape class depends on the entry point"
                checked += 1
                break
            if t % ALPHA == 0:
                saved, saved_t = C, t
            state = nxt
    return events, checked


def main():
    cfg = json.load(sys.stdin)
    trials = cfg["trials"]
    rng = random.Random(7)
    tot = 0
    stats = {"events": 0, "cycle_len": {}}
    for _ in range(trials):
        x = rng.randrange(n)
        got, steps = solve(x, rng, stats=stats)
        assert got == x, (x, got)
        tot += steps
    ev, ok = cycle_test(rng)
    print(json.dumps({"group": "y^2=x^3+2 over F_30091, prime order 29881", "classes": (n - 1) // 6,
                      "trials": trials, "all_recovered": True,
                      "mean_steps": tot / trials, "sqrt_class_count": math.sqrt((n - 1) / 6),
                      "cycle_events_in_solves": stats["events"], "cycle_lengths": stats["cycle_len"],
                      "cycle_test_events": ev, "cycle_escape_entry_independent": ok == ev,
                      "scope": "toy group only; does not validate production cost bounds"}))


main()
