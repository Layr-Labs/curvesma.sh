"""Exhaustive coverage check for negation-expanded exact-interval BSGS.

In an additive cyclic group the compressed key is the scalar. Field negation is
the group negative, so the intersection is the same one the production encoding
implements by flipping the 0x20 sign flag. For every admitted order the
parameter pair meets h*(2m-1) >= ceil(39*order/100) and h*(2m-1) <= order, and
the scan recovers x exactly on the claimed residue set.

This is evidence about the coverage identity only. It says nothing about
production point encodings, the 2^128-record sorts, instruction counts or the
storage bound, which are accounted analytically in proof.md.
"""
import json
import sys

def _is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1 if d == 2 else 2
    return True


order = json.load(sys.stdin)["order"]
if not (2 <= order <= 10000):
    raise SystemExit("order out of range")

TARGET = (39 * order + 99) // 100
BUCKET = 12 * (1 << 32) * 6


def ledger(m, h):
    records = (2 * m - 1) + h
    group_ops = (m - 1) + (h - 1) + 2 + 512
    instructions = 920 * records + 16 * (m - 1) + 2 * BUCKET + 4096
    return 2000 * group_ops + instructions


best = None
for m in range(1, order):
    s = 2 * m - 1
    if s >= order or 2 * (m - 1) >= order:
        break
    h = (TARGET + s - 1) // s
    if h < 1 or h * s > order:
        continue
    num = ledger(m, h)
    if best is None or num < best[0]:
        best = (num, m, h, s)
if best is None:
    raise SystemExit("no feasible parameters")

_, m, h, s = best
babies = {}
for j in range(m):
    babies.setdefault(j % order, []).append(j)
    if j:
        babies.setdefault((-j) % order, []).append(-j)

claimed = set()
reps = {}
for i in range(h):
    for e in range(-(m - 1), m):
        residue = (i * s + e) % order
        if residue in reps:
            raise AssertionError("duplicate representation at %d" % residue)
        reps[residue] = (i, e)
        claimed.add(residue)

if len(claimed) != h * s:
    raise AssertionError("cardinality %d != %d" % (len(claimed), h * s))

recovered = 0
for x in range(order):
    found = None
    matches = 0
    for i in range(h):
        labels = babies.get((x - i * s) % order)
        if not labels:
            continue
        for e in labels:
            matches += 1
            cand = (i * s + e) % order
            if cand == x:
                found = cand
    if (found is not None) != (x in claimed):
        raise AssertionError("coverage mismatch at x=%d" % x)
    if found is not None:
        recovered += 1
    if _is_prime(order) and matches > 1:
        raise AssertionError("prime order emitted %d matches at x=%d" % (matches, x))

if recovered != len(claimed):
    raise AssertionError("recovered %d claimed %d" % (recovered, len(claimed)))

print(json.dumps({
    "group": "additive cyclic group of the given order",
    "order": order,
    "m": m,
    "h": h,
    "stride": s,
    "coverage": h * s,
    "checked": order,
    "recovered": recovered,
    "exact_signed_interval_coverage": True,
    "scope": "toy group only; does not validate production encodings or cost bounds",
}))
