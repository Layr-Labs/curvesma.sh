"""Truncated negation-class BSGS sanity check on a toy cyclic group.

The production construction runs on secp256k1; here we run the identical
indexing, window coverage and hash-probe accounting on the additive group
Z/nZ, where the key function k(v) = min(v mod n, (-v) mod n) plays the role
of the x-coordinate (two scalars share a key iff they are equal or opposite).
This is evidence only at toy scale: it checks that every covered secret is
recovered, that uncovered secrets contribute the declared coverage fraction,
and that the simulated linear-probing hash table stays inside the charged
probe budget. It does not validate production costs.
"""
import json
import math
import sys

n = json.load(sys.stdin)["order"]
if not (3 <= n <= 100000):
    raise ValueError("order out of range")
m = math.isqrt(n // 5) + 1
stride = 2 * m - 1
h = (2 * n + 5 * stride - 1) // (5 * stride)   # ceil(0.4 n / stride)
cap = ((m * 10) // 7) + 1                      # load factor 0.7

def key(v):
    v %= n
    return min(v, n - v)

def hsh(k):
    # multiplicative fold on a fixed-width key; mirrors the production
    # double-odd multiply-and-shift hash on the 256-bit x-coordinate
    return (k * 11400714819323198485) % cap

# Baby table: open addressing with linear probing; slots hold (key, j)
# or None. j = 0 is not stored: t = 0 is handled by the gamma=0 branch.
table = [None] * cap
insert_probes = 0
insert_ops = 0
worst_insert = 0
for j in range(1, m):
    k = j % n
    p = hsh(k)
    probes = 1
    while table[p] is not None and table[p][0] != k:
        p = (p + 1) % cap
        probes += 1
    table[p] = (k, j)
    insert_probes += probes
    insert_ops += 1
    worst_insert = max(worst_insert, probes)

# Covered secrets are x + (m - 1) mod n lying in the h disjoint giant
# windows {in - (m-1) + t : |t| <= m-1}, i.e. x mod n in (in - n, in].
covered = set()
for i in range(h):
    for u in range(i * stride - stride + 1, i * stride + 1):
        covered.add(u % n)

lookup_probes = 0
lookup_ops = 0
worst_lookup = 0
hits = 0
for secret in range(n):
    gamma = (secret + m - 1) % n
    recovered = None
    for i in range(h):
        if gamma == 0:
            recovered = (i * stride - (m - 1)) % n
            break
        k = key(gamma)
        p = hsh(k)
        probes = 1
        while table[p] is not None:
            if table[p][0] == k:
                j = table[p][1]
                c = (i * stride - (m - 1)) % n
                for a in ((c + j) % n, (c - j) % n):
                    if a == secret:
                        recovered = a
                        break
                break
            p = (p + 1) % cap
            probes += 1
        lookup_probes += probes
        lookup_ops += 1
        if recovered is not None:
            hits += 1
            break
        gamma = (gamma - stride) % n
    expect = secret in covered
    if (recovered is not None) != expect:
        raise AssertionError(
            "coverage contract broken at %d" % secret)

mean_insert = insert_probes / max(1, insert_ops)
mean_lookup = lookup_probes / max(1, lookup_ops)
# Charged budget assumes <= 6.06 mean probes per operation at load 0.7
# (textbook unsuccessful linear-probing mean); keep a 1.5x safety factor.
if mean_insert > 6.06 * 1.5 or mean_lookup > 6.06 * 1.5:
    raise AssertionError("probe budget exceeded")
if h * stride < (2 * n) // 5:
    raise AssertionError("coverage below 0.4 n")
print(json.dumps({
    "group": "integers modulo order under addition",
    "order": n,
    "m": m,
    "stride": stride,
    "giant_steps": h,
    "covered_scalars": len(covered),
    "coverage_fraction": len(covered) / n,
    "required_fraction": 0.4,
    "hits": hits,
    "hits_match_covered": hits == len(covered),
    "mean_insert_probes": mean_insert,
    "mean_lookup_probes": mean_lookup,
    "worst_insert_probes": worst_insert,
    "worst_lookup_probes": worst_lookup,
    "scope": "toy group only; does not validate production cost bounds"
}))
