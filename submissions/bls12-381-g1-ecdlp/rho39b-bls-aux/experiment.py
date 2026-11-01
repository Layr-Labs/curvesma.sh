"""Toy-scale validation of the parallel distinguished-points rho attack.

Runs P = 8 walks with a shared distinguished-point table (r-adding walk
partitioned by x-coordinate, collision solve, independent verification) on
a small prime-order elliptic curve. Measures the empirical success rate
against the birthday prediction 1 - exp(-T(T-1)/2n) used in proof.md. This
validates algorithm logic and heuristic fit at small scale only; it does
not execute or prove production cost bounds.
"""
import json
import math
import random
import sys

R_PARTITIONS = 16
WALKS = 8
DISTINGUISHED_MOD = 4  # toy density 1/4 (production uses 2^-32)
TAIL_FACTOR = 8
TRIALS = 2000
MIN_RATE = 0.34  # nominal prediction ~0.40; margin for noise

# Fixed toy curve: order -> (field prime, a, b). Group order is prime.
CURVES = {
    967: (1009, 0, 11),
    9883: (9973, 1, 21),
}

INF = None


def padd(P, Q, a, p):
    if P is INF:
        return Q
    if Q is INF:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return INF
        s = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        s = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (s * s - x1 - x2) % p
    return (x3, (s * (x1 - x3) - y1) % p)


def pmul(k, P, a, p):
    R = INF
    while k:
        if k & 1:
            R = padd(R, P, a, p)
        P = padd(P, P, a, p)
        k >>= 1
    return R


def group_order(a, b, p):
    n = 1
    for x in range(p):
        r = (pow(x, 3, p) + a * x + b) % p
        if r == 0:
            n += 1
        elif pow(r, (p - 1) // 2, p) == 1:
            n += 2
    return n


def find_point(a, b, p):
    for x in range(p):
        r = (pow(x, 3, p) + a * x + b) % p
        if r == 0:
            return (x, 0)
        if pow(r, (p - 1) // 2, p) == 1:
            return (x, tonelli(r, p))
    raise AssertionError("no point found")


def tonelli(r, p):
    # p prime, r a quadratic residue. Brute force is fine at toy scale.
    for y in range(1, p):
        if (y * y) % p == r:
            return y
    raise AssertionError("no sqrt")


def walk_length(n):
    # Same formula as production: smallest T with T(T-1) >= 2*n*0.5113.
    target = 2 * n * 0.5113
    t = int(math.sqrt(target)) + 1
    while t * (t - 1) < target:
        t += 1
    return t


def partition(Y):
    if Y is INF:
        return 0
    return Y[0] % R_PARTITIONS


def distinguished(Y):
    return Y is not INF and Y[0] % DISTINGUISHED_MOD == 0


def rho_attack(n, P, Q, a, b_curve, p, secret, rng):
    t_walk = walk_length(n)
    per_walk = t_walk // WALKS + 1
    tail = TAIL_FACTOR * DISTINGUISHED_MOD
    table = {}
    for _ in range(WALKS):
        u = [rng.randrange(n) for _ in range(R_PARTITIONS)]
        v = [rng.randrange(n) for _ in range(R_PARTITIONS)]
        M = [padd(pmul(u[j], P, a, p), pmul(v[j], Q, a, p), a, p)
             for j in range(R_PARTITIONS)]
        ca, cb = rng.randrange(n), rng.randrange(n)
        Y = padd(pmul(ca, P, a, p), pmul(cb, Q, a, p), a, p)
        # Production keys by (x, y-parity); full-point keys are equivalent
        # here (same point iff same (x,y)) and simpler at toy scale.
        for _ in range(per_walk + tail):
            f = partition(Y)
            ca = (ca + u[f]) % n
            cb = (cb + v[f]) % n
            Y = padd(Y, M[f], a, p)
            if distinguished(Y):
                if Y in table:
                    a0, b0 = table[Y]
                    da = (ca - a0) % n
                    db = (cb - b0) % n
                    if db != 0:
                        rec = (-da * pow(db, -1, n)) % n
                        if pmul(rec, P, a, p) == Q:
                            return rec
                        raise AssertionError("wrong scalar from solve")
                else:
                    table[Y] = (ca, cb)
    return None


def main():
    order = json.load(sys.stdin)["order"]
    assert isinstance(order, int) and 2 <= order <= 10000
    assert order in CURVES, f"unsupported toy order {order}"
    p, ca, cb = CURVES[order]
    n = group_order(ca, cb, p)
    assert n == order, "curve order mismatch"
    P = find_point(ca, cb, p)
    assert pmul(n, P, ca, p) is INF, "generator order wrong"
    rng = random.Random(0xC0FFEE + order)
    t_walk = walk_length(n)
    pred = 1.0 - math.exp(-t_walk * (t_walk - 1) / (2 * n))
    ok = 0
    for _ in range(TRIALS):
        secret = rng.randrange(n)
        Q = pmul(secret, P, ca, p)
        if rho_attack(n, P, Q, ca, cb, p, secret, rng) is not None:
            ok += 1
    rate = ok / TRIALS
    assert rate >= MIN_RATE, f"success rate {rate:.3f} below {MIN_RATE}"
    print(json.dumps({
        "curve": f"y^2=x^3+{ca}x+{cb} over F_{p}",
        "order": n,
        "walks": WALKS,
        "trials": TRIALS,
        "successes": ok,
        "success_rate": round(rate, 4),
        "predicted_walk_collision": round(pred, 4),
        "t_walk": t_walk,
        "all_recovered_correct": True,
        "scope": "toy curve only; does not validate production cost bounds",
    }))


main()
