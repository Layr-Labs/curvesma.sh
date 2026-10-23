"""Finite truncated-coverage BSGS check; not evidence of full-size costs or security."""
import json
import math
import sys

n = json.load(sys.stdin)["order"]
m = max(1, math.isqrt((2 * n) // 5) + 1)
H = (2 * n + 5 * m - 1) // (5 * m)
babies = {j % n: j for j in range(m)}


def giant_walk(secret):
    value = secret
    for i in range(H):
        if value in babies:
            return (i * m + babies[value]) % n
        value = (value - m) % n
    return None


covered = min(H * m, n)
for secret in range(covered):
    assert giant_walk(secret) == secret, "covered scalar not recovered"
for secret in range(covered, n):
    assert giant_walk(secret) is None, "uncovered scalar unexpectedly recovered"
assert 5 * H * m >= 2 * n, "coverage below 2/5"
print(json.dumps({"group": "integers modulo order under addition", "order": n,
                  "baby_steps": m, "giant_steps": H,
                  "covered": covered, "uncovered": n - covered,
                  "coverage": covered / n, "truncated": H * m < n,
                  "all_covered_recovered": True,
                  "all_uncovered_exhaust_loop": True,
                  "scope": "toy group only; does not validate production cost bounds"}))
