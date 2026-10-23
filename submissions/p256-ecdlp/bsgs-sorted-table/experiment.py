"""Finite sorted-table BSGS sanity check; not evidence of full-size costs."""
import json
import math
import sys

n = json.load(sys.stdin)["order"]
m = math.isqrt(n) + 1
h = (n + m - 1) // m

records = sorted((j % n, j) for j in range(m))
keys = [record[0] for record in records]
assert len(set(keys)) == len(keys), "duplicate key in the baby table"


def lookup(key):
    lo, hi = 0, len(keys)
    while lo < hi:
        mid = (lo + hi) // 2
        if keys[mid] < key:
            lo = mid + 1
        else:
            hi = mid
    if lo < len(keys) and keys[lo] == key:
        return records[lo][1]
    return None


checked = 0
for secret in range(n):
    value = secret % n
    recovered = None
    for i in range(h):
        j = lookup(value)
        if j is not None:
            recovered = (i * m + j) % n
            break
        value = (value - m) % n
    assert recovered == secret % n, (secret, recovered)
    checked += 1

print(json.dumps({
    "group": "integers modulo order under addition",
    "order": n,
    "m": m,
    "giant_steps": h,
    "records": len(records),
    "identity_recovered": True,
    "checked": checked,
    "all_recovered": True,
    "scope": "toy group only; does not validate production cost bounds",
}))
