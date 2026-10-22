"""Toy negation-map Pollard rho; not evidence of full-size costs or security."""
import json
import random
import sys

# Prime-order Weierstrass curve y^2 = x^3 + a x + b over Fp.
p = 53
a = 1
b = 8
r = 61
Px, Py = 1, 13
G = (Px, Py)


def finv(z):
    return pow(z % p, p - 2, p)


def ginv(z):
    return pow(z % r, r - 2, r)


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + a) * finv(2 * y1) % p
    else:
        lam = (y2 - y1) * finv(x2 - x1) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def smul(k, P):
    k %= r
    acc = None
    base = P
    while k:
        if k & 1:
            acc = add(acc, base)
        base = add(base, base)
        k >>= 1
    return acc


def canon(P):
    if P is None:
        return None
    x, y = P
    return (x, min(y, p - y))


def h(P, npart):
    if P is None:
        return 0
    return P[0] % npart


def check_curve():
    assert (Py * Py - (Px * Px * Px + a * Px + b)) % p == 0
    assert smul(r, G) is None
    assert smul(1, G) == G


def solve_collision(a1, b1, a2, b2, Q):
    den = (b2 - b1) % r
    if den != 0:
        x = ((a1 - a2) * ginv(den)) % r
        if smul(x, G) == Q:
            return x
    den = (b2 + b1) % r
    if den != 0:
        x = ((-a1 - a2) * ginv(den)) % r
        if smul(x, G) == Q:
            return x
    return None


def recover_one(secret, rng):
    Q = smul(secret, G)
    if Q is None:
        return secret == 0
    npart = 8
    table = []
    for _ in range(npart):
        c = rng.randrange(1, r)
        d = rng.randrange(1, r)
        M = add(smul(c, G), smul(d, Q))
        table.append((c, d, M))
    dps = {}
    for _walk in range(32):
        u = rng.randrange(r)
        v = rng.randrange(r)
        raw = add(smul(u, G), smul(v, Q))
        X = canon(raw)
        if X is None:
            continue
        alpha, beta = u, v
        if X[1] != raw[1]:
            alpha, beta = (-u) % r, (-v) % r
        for _step in range(128):
            i = h(X, npart)
            c, d, M = table[i]
            Y = add(X, M)
            a2, b2 = (alpha + c) % r, (beta + d) % r
            tries = 0
            while Y is not None and h(canon(Y), npart) == i and tries < npart:
                i = (i + 1) % npart
                c, d, M = table[i]
                Y = add(X, M)
                a2, b2 = (alpha + c) % r, (beta + d) % r
                tries += 1
            if Y is None:
                break
            Xn = canon(Y)
            if Xn[1] != Y[1]:
                a2, b2 = (-a2) % r, (-b2) % r
            X, alpha, beta = Xn, a2, b2
            if X[0] % 4 != 0:
                continue
            key = X[0]
            if key in dps:
                a1, b1 = dps[key]
                if solve_collision(a1, b1, alpha, beta, Q) is not None:
                    return True
            else:
                dps[key] = (alpha, beta)
    return False


def recover(trials):
    check_curve()
    rng = random.Random(1)
    recovered = 0
    for _ in range(trials):
        secret = rng.randrange(r)
        if recover_one(secret, rng):
            recovered += 1
    return recovered


inp = json.load(sys.stdin)
trials = int(inp["trials"])
ok = recover(trials)
print(json.dumps({
    "curve": "toy Weierstrass p=53 a=1 b=8 r=61 P=(1,13)",
    "trials": trials,
    "recovered": ok,
    "all_recovered": ok == trials,
    "scope": "toy curve only; does not validate production cost bounds",
}))
