"""Toy end-to-end check of the stored six-leaf birthday construction (v2).

Mirrors the production structure that review found untested in v1:
  * the Y axis (and its scalars) is drawn once and reused for every X,
  * leaves live in a chained bucket table indexed by the low bits of x, with no
    overwrite, so every record stays reachable,
  * a collision is resolved by the two linear equations and verified with
    [x]P = Q.

Curve: y^2 = x^3 + 7 over F_127 (prime order 127, p = 1 mod 3).
Input: JSON on stdin, {"seed": int, "trials": int}. Output: JSON on stdout.
Toy scale only; it supports no production cost or success claim.
"""
import json
import random
import sys

P_MOD = 127
A_CURVE, B_CURVE = 0, 7
BUCKETS = 64            # toy bucket count; production uses 2^129


def inv(a, m):
    return pow(a % m, -1, m)


def padd(P, Q, m=P_MOD):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % m == 0:
        return None
    lam = ((3 * x1 * x1 + A_CURVE) * inv(2 * y1, m) if P == Q
           else (y2 - y1) * inv(x2 - x1, m)) % m
    x3 = (lam * lam - x1 - x2) % m
    return (x3, (lam * (x1 - x3) - y1) % m)


def pmul(k, P, m=P_MOD):
    R, Q = None, P
    while k:
        if k & 1:
            R = padd(R, Q, m)
        Q = padd(Q, Q, m)
        k >>= 1
    return R


def point_order(P, m=P_MOD):
    o, Q = 1, P
    while Q is not None and o <= 3 * m:
        Q = padd(Q, P, m)
        o += 1
    return o


def curve_order(m=P_MOD):
    n = 1
    for x in range(m):
        rhs = (x * x * x + A_CURVE * x + B_CURVE) % m
        for y in range(m):
            if y * y % m == rhs:
                n += 1
    return n


def find_generator(m=P_MOD):
    order = curve_order(m)
    for x in range(1, m):
        rhs = (x * x * x + A_CURVE * x + B_CURVE) % m
        for y in range(m):
            if y * y % m == rhs and point_order((x, y), m) == order:
                return (x, y), order
    raise AssertionError("no full-order point")


def cube_root_unity(m=P_MOD):
    for b in range(2, m):
        if (b * b + b + 1) % m == 0:
            return b
    raise AssertionError


def main():
    cfg = json.load(sys.stdin)
    rng = random.Random(int(cfg.get("seed", 1)))
    trials = int(cfg.get("trials", 3))

    G, N = find_generator()
    beta = cube_root_unity()
    phi = lambda P: None if P is None else ((beta * P[0]) % P_MOD, P[1])
    lam = next(k for k in range(2, N) if pmul(k, G) == phi(G))
    bt = [1, beta, beta * beta % P_MOD]

    recovered = 0
    for _ in range(trials):
        x = rng.randrange(1, N)
        Q = pmul(x, G)

        n = 3
        A = [rng.randrange(1, N) for _ in range(n)]
        B = [rng.randrange(1, N) for _ in range(n)]
        C = [rng.randrange(1, N) for _ in range(n)]
        Y = [pmul(C[l], G) for l in range(n)]

        buckets = [[] for _ in range(BUCKETS)]
        found = None
        for k in range(n):
            X = padd(pmul(A[k], G), pmul(B[k], Q))
            for l in range(n):
                xX, yX = X
                xY, yY = Y[l]
                d = [(xY - bt[t] * xX) % P_MOD for t in range(3)]
                if 0 in d:
                    continue
                D = d[0] * d[1] % P_MOD * d[2] % P_MOD
                Di = inv(D, P_MOD)
                u = [(yY - yX) % P_MOD, (yY + yX) % P_MOD]
                for t in range(3):
                    w = Di * d[(t + 1) % 3] % P_MOD * d[(t + 2) % 3] % P_MOD
                    for s in range(2):
                        xx = (u[s] * u[s] % P_MOD * w % P_MOD * w - xY - bt[t] * xX) % P_MOD
                        b = xx % BUCKETS
                        for other in buckets[b]:
                            if other["x"] == xx:
                                found = (other, {"x": xx, "k": k, "l": l, "t": t, "s": s})
                                break
                        if found:
                            break
                        buckets[b].append({"x": xx, "k": k, "l": l, "t": t, "s": s})
                    if found:
                        break
                if found:
                    break
            if found:
                break
        assert found, "toy run found no collision"

        # resolve: m_i = lam^t a +- c, q_i = lam^t b
        def coeffs(rec):
            sgn = 1 if rec["s"] == 0 else -1
            return ((pow(lam, rec["t"], N) * A[rec["k"]] + sgn * C[rec["l"]]) % N,
                    (pow(lam, rec["t"], N) * B[rec["k"]]) % N)

        r1, r2 = found
        m1, q1 = coeffs(r1)
        m2, q2 = coeffs(r2)
        got = None
        for num, den in ((m2 - m1, q1 - q2), (-(m1 + m2), q1 + q2)):
            den %= N
            if den == 0:
                continue
            cand = (num % N) * inv(den, N) % N
            if pmul(cand, G) == Q:
                got = cand
                break
        assert got == x, ("toy recovery failed", got, x)
        recovered += 1

    print(json.dumps({
        "curve": "y^2 = x^3 + 7 over F_%d" % P_MOD,
        "group_order": N,
        "beta": beta,
        "lambda": lam,
        "buckets": BUCKETS,
        "recovered": recovered,
        "result": "PASS",
    }, indent=1))


if __name__ == "__main__":
    main()
