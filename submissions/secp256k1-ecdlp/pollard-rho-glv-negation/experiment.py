"""Toy Pollard rho with negation map and GLV quotient walk.

Implements the exact walk logic of proof.md on the additive cyclic group
Z/qZ (prime q, q = 1 mod 3): the endomorphism is multiplication by a cube
root of unity lambda mod q (found by search), classes are {+-lambda^j * R},
walks use a 32-entry r-adding table, distinguished points, cycle detection
with restart, and collision resolution. This checks algorithm logic only;
it does not validate production cost bounds, memory bounds, or secp256k1
arithmetic.
"""
import json
import math
import random
import sys


def find_lambda(q):
    # Smallest lambda > 1 with lambda^3 = 1 mod q, lambda != 1.
    for lam in range(2, q):
        if pow(lam, 3, q) == 1:
            return lam
    raise AssertionError("no cube root of unity mod q")


def canonical(R, a, b, lam, q):
    # Among {+-lam^j * R}, pick the representative with the smallest
    # value in [0, q); transform (a, b) by the same factor.
    best = None
    for j in range(3):
        f = pow(lam, j, q)
        for s in (1, q - 1):
            cand = (s * f * R) % q
            key = (cand, 0 if s == 1 else 1, j)
            if best is None or key < best[0]:
                best = (key, cand, (s * f * a) % q, (s * f * b) % q)
    _, Rc, ac, bc = best
    return Rc, ac, bc


def rho(q, secret, rng, max_steps=500000, max_restarts=200):
    lam = find_lambda(q)
    # Class-structure check: orbits of size 6 except 0; lam^3 = 1.
    assert pow(lam, 3, q) == 1 and lam != 1
    for R in range(1, min(q, 500)):
        orbit = {(s * pow(lam, j, q) * R) % q
                 for j in range(3) for s in (1, q - 1)}
        assert len(orbit) == 6
    # 32-entry r-adding table: T_k = ta_k + tb_k * secret (mod q).
    table = []
    for _ in range(32):
        ta, tb = rng.randrange(q), rng.randrange(q)
        table.append((ta, tb, (ta + tb * secret) % q))
    distinguished = {}  # R -> (a, b)
    steps = 0
    restarts = 0
    while steps < max_steps:
        # Fresh random walk start: R = a + b * secret.
        a, b = rng.randrange(q), rng.randrange(q)
        R = (a + b * secret) % q
        seen = set()
        while steps < max_steps:
            Rc, a, b = canonical(R, a, b, lam, q)
            if Rc in seen:
                # Fruitless cycle: abort this walk and restart.
                restarts += 1
                if restarts > max_restarts:
                    raise AssertionError("too many restarts")
                break
            seen.add(Rc)
            k = Rc & 31  # low 5 bits: partition index
            ta, tb, Tk = table[k]
            R = (Rc + Tk) % q
            a = (a + ta) % q
            b = (b + tb) % q
            steps += 1
            if R % 16 == 0:  # distinguished point (toy threshold 2^4)
                if R in distinguished:
                    a1, b1 = distinguished[R]
                    if b != b1:
                        inv = pow((b - b1) % q, -1, q)
                        x = ((a1 - a) * inv) % q
                        assert (x == secret), "recovered wrong secret"
                        return {"q": q, "lambda": lam, "secret": secret,
                                "recovered": x, "steps": steps,
                                "restarts": restarts,
                                "distinguished_stored": len(distinguished)}
                    # Trivial collision (same coefficients): keep walking.
                else:
                    distinguished[R] = (a, b)
                break  # walk ends at a distinguished point; restart
    raise AssertionError("no collision within step budget")


inp = json.load(sys.stdin)
q, secret = inp["order"], inp["secret"]
assert 2 <= q and (q - 1) % 3 == 0 and 1 <= secret < q
rng = random.Random(12345)
result = rho(q, secret, rng)
print(json.dumps({"group": "integers modulo prime order under addition",
                  "rho": result,
                  "scope": "toy group only; does not validate production cost bounds"}))
