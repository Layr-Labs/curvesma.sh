"""Finite signed-window BSGS sanity check; not evidence of full-size costs."""
import json
import math
import sys

n = json.load(sys.stdin)["order"]
m = math.isqrt(n // 2) or 1
stride = 2 * m - 1
h = (n + m - 1 + stride - 1) // stride
babies = {}
for j in range(m):
    key = min(j % n, (-j) % n)
    babies[key] = j
for secret in range(n):
    recovered = None
    gamma = (secret + (m - 1)) % n
    for i in range(h):
        t = gamma if gamma <= n // 2 else n - gamma
        if t in babies:
            j = babies[t]
            c = (i * stride - (m - 1)) % n
            for a in ((c + j) % n, (c - j) % n):
                if a == secret:
                    recovered = a
                    break
            if recovered is not None:
                break
        gamma = (gamma - stride) % n
    if recovered != secret:
        raise AssertionError("Discrete logarithm not recovered")
print(json.dumps({
    "group": "integers modulo order under addition",
    "order": n,
    "m": m,
    "stride": stride,
    "h": h,
    "checked": n,
    "all_recovered": True,
    "scope": "toy group only; does not validate production cost bounds"
}))
