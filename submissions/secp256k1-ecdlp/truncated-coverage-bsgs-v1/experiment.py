"""Finite check of the truncated BSGS coverage; not evidence of full-size costs."""
import json
import math
import sys

COVERAGE_NUM = 2
COVERAGE_DEN = 5

data = json.load(sys.stdin)
n = data["order"]
m = math.isqrt(n) + 1
limit = (COVERAGE_NUM * n + COVERAGE_DEN - 1) // COVERAGE_DEN
steps = (limit + m - 1) // m
babies = {j % n: j for j in range(m)}
recovered = 0
for secret in range(n):
    value = secret
    for i in range(steps):
        if value in babies:
            a = (i * m + babies[value]) % n
            assert a == secret
            recovered += 1
            break
        value = (value - m) % n
fraction = recovered / n
assert steps * m >= limit
assert fraction + 1e-12 >= COVERAGE_NUM / COVERAGE_DEN
print(json.dumps({"group": "integers modulo order under addition", "order": n,
                  "checked": n, "baby_steps": m, "giant_steps": steps,
                  "covered_prefix": steps * m,
                  "recovered": recovered, "recovered_fraction": fraction,
                  "declared_coverage": COVERAGE_NUM / COVERAGE_DEN,
                  "scope": "toy group only; does not validate production cost bounds"}))
