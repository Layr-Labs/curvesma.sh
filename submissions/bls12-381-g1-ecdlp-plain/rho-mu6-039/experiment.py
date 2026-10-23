"""Toy-scale run of the walk in proof.md on small j=0 curves of prime order.

Parallel r-adding walks on classes of six points (negation and the GLV map),
canonical representative (min{x, beta*x, -x-beta*x}, even y), 2-cycle
look-ahead, restart on fruitless cycles. Every visited class is stored
(distinguished-point delay is accounted analytically in proof.md). Stops at the
claimed budget T = ceil(SLACK * sqrt(2 ln(1/0.61) * n/6)) total additions,
including look-ahead and restart work, and reports the fraction of secrets
recovered. Toy scope only; not a production cost proof.
"""
import json
import math
import random
import sys
import time

inp = json.load(sys.stdin)
BITS, TRIALS, SEED, R = inp["bits"], inp["trials"], inp["seed"], inp["r"]
WALKS = inp["walks"]
SLACK = 1.04
rng = random.Random(SEED)
t0 = time.time()


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    dd, s = n - 1, 0
    while dd % 2 == 0:
        dd //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, dd, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def sqrt_mod(a, p):
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, x = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, x = i, b * b % p, t * b * b % p, x * b % p
    return x


def add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0:
            return None
        lam = 3 * P[0] * P[0] * pow(2 * P[1], -1, p) % p
    else:
        lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x3 = (lam * lam - P[0] - Q[0]) % p
    return (x3, (lam * (P[0] - x3) - P[1]) % p)


def smul(k, P):
    out = None
    while k:
        if k & 1:
            out = add(out, P)
        P = add(P, P)
        k >>= 1
    return out


while True:
    p = rng.randrange(1 << (BITS - 1), 1 << BITS)
    if p % 3 != 1 or not is_prime(p):
        continue
    b = rng.randrange(1, p)
    x = rng.randrange(p)
    rhs = (x ** 3 + b) % p
    if pow(rhs, (p - 1) // 2, p) != 1:
        continue
    G = (x, sqrt_mod(rhs, p))
    # j=0 curves: #E = p+1-t with 4p = t^2+3v^2; test the six candidate orders.
    found = None
    for v in range(1, math.isqrt(4 * p // 3) + 1):
        t2 = 4 * p - 3 * v * v
        tt = math.isqrt(t2)
        if tt * tt == t2:
            for T in {tt, -tt, (tt + 3 * v) // 2, -(tt + 3 * v) // 2, (tt - 3 * v) // 2, -(tt - 3 * v) // 2}:
                n = p + 1 - T
                if is_prime(n) and smul(n, G) is None:
                    found = n
            break
    if found:
        n = found
        break

beta = next(c for c in (pow(a, (p - 1) // 3, p) for a in range(2, 50)) if c != 1)
lamb = next(l for l in (pow(g, (n - 1) // 3, n) for g in range(2, 50)) if l != 1
            and smul(l, G) == (beta * G[0] % p, G[1]))
lamb2 = lamb * lamb % n
assert smul(lamb2, G) == (beta * beta * G[0] % p, G[1])


def canon(P, a, c):
    """Canonical class representative and its (a, c) coefficients: P = aG + cQ."""
    x, y = P
    bx = beta * x % p
    b2x = (-x - bx) % p
    if x <= bx and x <= b2x:
        nx, mu = x, 1
    elif bx <= b2x:
        nx, mu = bx, lamb
    else:
        nx, mu = b2x, lamb2
    if y & 1:
        y, mu = p - y, n - mu
    return (nx, y), a * mu % n, c * mu % n


def h(P):
    return (P[0] * 0x9E3779B1 >> 7) % R


budget = math.ceil(SLACK * math.sqrt(2 * math.log(1 / 0.61) * n / 6))
succ = done = adds_total = lookahead_adds = restarts = useless = 0
for trial in range(TRIALS):
    if time.time() - t0 > 15:
        break
    secret = rng.randrange(1, n)
    Q = smul(secret, G)
    table = []
    for _ in range(R):
        ca, cc = rng.randrange(n), rng.randrange(n)
        table.append((add(smul(ca, G), smul(cc, Q)), ca, cc))

    def fresh():
        a0, c0 = rng.randrange(n), rng.randrange(n)
        return canon(add(smul(a0, G), smul(c0, Q)), a0, c0)

    walks = [list(fresh()) + [k] for k in range(WALKS)]
    seen = {}
    steps, got = 0, None
    while steps < budget and got is None:
        for w in walks:
            if steps >= budget:
                break
            P, a, c, wid = w
            i = h(P)
            # 2-cycle look-ahead: skip an index that would map straight back.
            for k in range(R):
                j = (i + k) % R
                Y = add(P, table[j][0])
                steps += 1
                if k:
                    lookahead_adds += 1
                if Y is None:
                    break
                Yc, ya, yc = canon(Y, a + table[j][1], c + table[j][2])
                if h(Yc) != j:
                    break
            if Y is None:
                w[:3] = fresh()
                restarts += 1
                continue
            key = Yc[0]
            if key in seen:
                sa, sc = seen[key]
                if (yc - sc) % n:
                    got = (sa - ya) * pow(yc - sc, -1, n) % n
                    break
                useless += 1
                w[:3] = fresh()
                restarts += 1
                continue
            seen[key] = (ya, yc)
            w[0], w[1], w[2] = Yc, ya, yc
    adds_total += steps
    done += 1
    if got is not None and got == secret:
        succ += 1
    elif got is not None:
        raise AssertionError("wrong logarithm")

print(json.dumps({
    "params": {"p": p, "b": b, "n": n, "r": R, "walks": WALKS, "slack": SLACK},
    "budget_additions_per_trial": budget,
    "random_mapping_steps_for_0.39": math.sqrt(2 * math.log(1 / 0.61) * n / 6),
    "trials_run": done,
    "recovered": succ,
    "success_rate_within_budget": succ / done if done else None,
    "binomial_sd_at_0.39": math.sqrt(0.39 * 0.61 / done) if done else None,
    "lookahead_extra_additions_fraction": lookahead_adds / max(adds_total, 1),
    "fruitless_or_useless_restarts": restarts,
    "useless_self_collisions": useless,
    "seconds": round(time.time() - t0, 2),
    "scope": "toy curve only; checks the class walk and its success at the claimed budget, not production cost",
}))
