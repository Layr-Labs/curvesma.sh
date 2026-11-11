"""Toy Pollard rho with the negation map (no endomorphism).

Implements the exact walk logic of proof.md on the additive cyclic group
Z/qZ (prime q): classes are {R, -R}, canonicalised to the smaller of the
two representatives; walks use a 32-entry r-adding table, distinguished
points, cycle detection with restart, and collision resolution. This
checks algorithm logic only; it does not validate production cost
bounds, memory bounds, or P-256 arithmetic.
"""
import json
import random
import sys


def canonical(R, a, b, q):
    # Among {R, -R}, pick the representative with the smaller value in
    # [0, q); transform (a, b) by the same sign.
    if R <= q - R:
        return R, a, b
    return q - R, (-a) % q, (-b) % q


def rho(q, secret, rng, max_steps=500000, max_restarts=200):
    # Class-structure check: orbits of size 2 except 0.
    for R in range(1, min(q, 500)):
        assert len({R, (-R) % q}) == 2
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
            Rc, a, b = canonical(R, a, b, q)
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
                        assert x == secret, "recovered wrong secret"
                        return {"q": q, "secret": secret,
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
assert 2 <= q and 1 <= secret < q
rng = random.Random(12345)
result = rho(q, secret, rng)
print(json.dumps({"group": "integers modulo prime order under addition",
                  "rho": result,
                  "scope": "toy group only; does not validate production cost bounds"}))
