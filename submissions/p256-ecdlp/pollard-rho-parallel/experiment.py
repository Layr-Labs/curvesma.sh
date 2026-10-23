"""Toy parallel Pollard rho on group points; not evidence of full-size costs.

Runs the same construction as proof.md end-to-end on a small prime-order
Weierstrass curve: an R-adding walk with high-bit partition, look-ahead,
distinguished points, and collision solving by re-derived coefficients.
Negation map intentionally absent. Evidence only at this scope: it checks
that the measured steps-to-solve track the sqrt(pi r / 2) idealization and
that every detected collision carries a solvable relation.
"""
import json
import math
import random
import sys

# Prime-order Weierstrass curve y^2 = x^3 + b over Fp, #E = r prime.
p = 1999
a = 0
b = 3
r = 2089
G = (1, 2)

RPART = 16        # partition count R
DPMASK = 7        # distinguished when x mod 8 == 0, theta = 1/8
WALKS = 8         # concurrent walks
LOOKAHEAD = True
STEP_CAP = 200000


def finv(z):
    return pow(z % p, p - 2, p)


def ginv(z):
    return pow(z % r, r - 2, r)


def add(Pt, Qt):
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    if Pt[0] == Qt[0] and (Pt[1] + Qt[1]) % p == 0:
        return None
    if Pt == Qt:
        if Pt[1] % p == 0:
            return None
        lam = (3 * Pt[0] * Pt[0] + a) * finv(2 * Pt[1]) % p
    else:
        lam = (Qt[1] - Pt[1]) * finv(Qt[0] - Pt[0]) % p
    x3 = (lam * lam - Pt[0] - Qt[0]) % p
    return (x3, (lam * (Pt[0] - x3) - Pt[1]) % p)


def smul(k, Pt):
    k %= r
    acc = None
    while k:
        if k & 1:
            acc = add(acc, Pt)
        Pt = add(Pt, Pt)
        k >>= 1
    return acc


def h(Pt):
    return Pt[0] % RPART


def check_curve():
    assert (G[1] * G[1] - (G[0] * G[0] * G[0] + a * G[0] + b)) % p == 0
    assert smul(r, G) is None
    assert smul(1, G) == G
    assert r > p  # prime order confirmed by construction search


def make_table(rng, Q):
    tab = []
    for _ in range(RPART):
        c = rng.randrange(r)
        d = rng.randrange(r)
        m = add(smul(c, G), smul(d, Q))
        if m is None:
            c, d = 1, 0
            m = G
        tab.append((c, d, m))
    return tab


def step(state, tab):
    """One look-ahead adding step. state = (point, alpha, beta)."""
    (X, al, be) = state
    i = h(X)
    for t in range(RPART):
        j = (i + t) % RPART
        c, d, m = tab[j]
        Y = add(X, m)
        if Y is None:
            continue
        if LOOKAHEAD and h(Y) == j:
            continue
        return (Y, (al + c) % r, (be + d) % r)
    return None


def solve_relation(a1, b1, a2, b2, s):
    # X1 = s * X2 gives (a1 - s*a2) + (b1 - s*b2) x = 0 mod r.
    den = (b1 - s * b2) % r
    if den == 0:
        return None
    return ((s * a2 - a1) * ginv(den)) % r


def recover_one(secret, rng):
    Q = smul(secret, G)
    if Q is None:
        return 0, 0, 0
    if Q == G:
        return 1, 0, 0
    tab = make_table(rng, Q)
    table = {}
    states = []
    for _ in range(WALKS):
        u, v = rng.randrange(r), rng.randrange(r)
        states.append([add(smul(u, G), smul(v, Q)), u, v])
    steps = 0
    fruitless = 0
    while steps < STEP_CAP:
        for st in states:
            nxt = step(tuple(st), tab)
            if nxt is None:
                u, v = rng.randrange(r), rng.randrange(r)
                st[0], st[1], st[2] = add(smul(u, G), smul(v, Q)), u, v
                continue
            st[0], st[1], st[2] = nxt
            steps += 1
            X = st[0]
            if X is None or (X[0] & DPMASK) != 0:
                continue
            key = X[0]
            if key in table:
                X2, a2, b2 = table[key]
                # Same x means X2 is X or -X; try both signs, as in proof.md.
                cand = None
                for s in (1, -1):
                    y2 = X2[1] if s == 1 else (-X2[1]) % p
                    if (X[0], y2 % p) == X:
                        cand = solve_relation(st[1], st[2], a2, b2, s)
                        break
                if cand is not None and smul(cand, G) == Q:
                    return cand, steps, fruitless
                fruitless += 1
                u, v = rng.randrange(r), rng.randrange(r)
                st[0], st[1], st[2] = add(smul(u, G), smul(v, Q)), u, v
            else:
                table[key] = (X, st[1], st[2])
    return None, steps, fruitless


def recover(trials):
    check_curve()
    rng = random.Random(0xC0FFEE)
    out = []
    total_steps = 0
    total_fruitless = 0
    recovered = 0
    for t in range(trials):
        secret = rng.randrange(r)
        x, steps, fl = recover_one(secret, rng)
        ok = x is not None and smul(x, G) == smul(secret, G)
        recovered += 1 if ok else 0
        total_steps += steps
        total_fruitless += fl
        out.append({"trial": t, "ok": ok, "steps": steps,
                    "fruitless_collisions": fl})
    ideal = math.sqrt(math.pi * r / 2.0)
    return {
        "curve": {"p": p, "a": a, "b": b, "r": r, "G": list(G)},
        "params": {"R": RPART, "theta": "1/8", "walks": WALKS,
                   "step_cap": STEP_CAP},
        "trials": trials,
        "recovered": recovered,
        "mean_steps_to_solve": total_steps / max(1, trials),
        "idealized_sqrt_pi_r_over_2": ideal,
        "steps_ratio_vs_ideal": (total_steps / max(1, trials)) / ideal,
        "fruitless_collisions": total_fruitless,
        "results": out,
    }


inp = json.load(sys.stdin)
trials = int(inp["trials"])
print(json.dumps(recover(min(trials, 64))))
