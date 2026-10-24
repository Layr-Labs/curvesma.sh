"""Finite BSGS sanity check; not evidence of full-size costs or security."""
import json
import math
import sys

n = json.load(sys.stdin)["order"]
m = math.isqrt(n) + 1
babies = {j % n: j for j in range(m)}
for secret in range(n):
    value = secret
    for i in range((n + m - 1) // m):
        if value in babies:
            recovered = (i * m + babies[value]) % n
            assert recovered == secret
            break
        value = (value - m) % n
    else:
        raise AssertionError("Discrete logarithm not recovered")
print(json.dumps({"group": "integers modulo order under addition", "order": n,
                  "checked": n, "all_recovered": True,
                  "scope": "toy group only; does not validate production cost bounds"}))
