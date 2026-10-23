"""Toy 20-adding walk. Not evidence of production cost or mixing."""
import json
import sys

n = json.load(sys.stdin)["order"]
# Branches v_s = 1, u_s = 2^(s mod 8) + 1, as in the production adder shape.
u = [(1 << (s % 8)) + 1 for s in range(20)]


def recover(secret):
    # y = alpha + beta * secret  in Z/nZ, start at 0.
    seen = {}
    alpha = 0
    beta = 0
    y = 0
    for step in range(1, 8 * n + 3):
        s = y % 20
        alpha = (alpha + u[s]) % n
        beta = (beta + 1) % n
        y = (alpha + beta * secret) % n
        if y % 4 == 0:
            if y in seen:
                a0, b0 = seen[y]
                db = (beta - b0) % n
                da = (a0 - alpha) % n
                if db == 0:
                    continue
                got = (da * pow(db, -1, n)) % n
                return got == secret
            seen[y] = (alpha, beta)
    return False


checked = sum(1 for secret in range(n) if recover(secret))
if checked != n:
    raise AssertionError("adding walk missed a toy discrete log")
print(json.dumps({
    "group": "integers modulo order under addition",
    "order": n,
    "checked": checked,
    "all_recovered": True,
    "scope": "toy group only; does not validate production cost bounds",
}))
