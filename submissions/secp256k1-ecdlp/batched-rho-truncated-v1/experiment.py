"""Scalar model of the double-or-add walk. Not a production cost or curve proof."""
import json
import sys


def invertible(value, modulus):
    try:
        return pow(value, -1, modulus)
    except ValueError:
        return None


n = json.load(sys.stdin)["order"]
steps = 1
while steps * (steps - 1) < 2 * n:
    steps += 1

recovered = 0
collisions = 0
algebra_ok = True
for secret in range(n):
    a = 0
    b = 1
    element = secret % n
    seen = {element: (a, b)}
    for _ in range(steps):
        branch = element % 32
        if branch < 16:
            a = (2 * a) % n
            b = (2 * b) % n
        else:
            a = (a + (branch - 15)) % n
            b = (b + 1 + (branch % 5)) % n
        element = (a + b * secret) % n
        if element not in seen:
            seen[element] = (a, b)
            continue
        collisions += 1
        a0, b0 = seen[element]
        delta_b = (b0 - b) % n
        inverse = invertible(delta_b, n)
        if inverse is None:
            break
        solved = ((a - a0) % n) * inverse % n
        if solved != secret % n:
            algebra_ok = False
        else:
            recovered += 1
        break

assert algebra_ok, "collision algebra failed"
assert steps * (steps - 1) >= 2 * n
if n == 1009:
    assert recovered / n >= 0.39, "birthday sample below the acceptance floor"

print(json.dumps({
    "group": "integers modulo order under addition",
    "order": n,
    "steps": steps,
    "secrets": n,
    "collisions": collisions,
    "recovered": recovered,
    "recovery_rate": recovered / n,
    "algebra_ok": algebra_ok,
    "exponent_rule": "steps*(steps-1) >= 2*order",
    "scope": "toy scalar walk only; does not validate production curve cost",
}))
