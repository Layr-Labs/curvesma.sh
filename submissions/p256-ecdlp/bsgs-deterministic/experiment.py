"""Toy deterministic BSGS; not evidence of full-size costs.

Same recurrence as proof.md on a prime-order Weierstrass curve of order 61:
baby points [i]P for i < m = ceil(sqrt(r)), giant steps by [m]P, solve
x = i + j*m. Lookup is exhaustive equality, which is exact on this group;
production cost uses the radix-sort bound in proof.md and is not measured here.
"""
import json
import math
import sys

p = 53
a = 1
b = 8
r = 61
G = (1, 13)


def finv(z):
    return pow(z % p, p - 2, p)


def add(Pt, Qt):
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    x1, y1 = Pt
    x2, y2 = Qt
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a) * finv(2 * y1) % p
    else:
        lam = (y2 - y1) * finv(x2 - x1) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def smul(k, Pt):
    k %= r
    acc = None
    base = Pt
    while k:
        if k & 1:
            acc = add(acc, base)
        base = add(base, base)
        k >>= 1
    return acc


def check_curve():
    assert (G[1] * G[1] - (G[0] ** 3 + a * G[0] + b)) % p == 0
    assert smul(r, G) is None
    assert smul(1, G) == G
    order = 1
    for x in range(p):
        rhs = (x * x * x + a * x + b) % p
        if rhs == 0:
            order += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            order += 2
    assert order == r


def recover_one(secret):
    Q = smul(secret, G)
    if Q is None:
        return 0
    m = math.isqrt(r)
    if m * m < r:
        m += 1
    table = {}
    cur = None
    for i in range(m):
        if i == 1:
            cur = G
        elif i >= 2:
            cur = add(cur, G)
        table[cur] = i
    xs = [pt[0] for pt in table if pt is not None]
    if len(xs) != len(set(xs)):
        return None
    step = smul(m, G)
    neg = None if step is None else (step[0], (-step[1]) % p)
    S = Q
    for j in range(m):
        if S in table:
            return table[S] + j * m
        S = add(S, neg)
    return None


def recover(trials):
    check_curve()
    # Full group, not a sample. `trials` is accepted so the runner schema
    # stays the usual object; the sweep does not grow with it.
    _ = trials
    results = []
    recovered = 0
    for secret in range(r):
        x = recover_one(secret)
        ok = x is not None and smul(x, G) == smul(secret, G) and x == secret
        recovered += 1 if ok else 0
        results.append({"secret": secret, "x": x, "ok": ok})
    return {
        "curve": {"p": p, "a": a, "b": b, "r": r, "G": list(G)},
        "m": math.isqrt(r) + (0 if math.isqrt(r) ** 2 == r else 1),
        "trials_requested": trials,
        "scalars_checked": r,
        "recovered": recovered,
        "all_ok": recovered == r,
        "results": results,
    }


inp = json.load(sys.stdin)
trials = int(inp["trials"])
if not 1 <= trials <= 64:
    raise SystemExit(2)
print(json.dumps(recover(trials)))
