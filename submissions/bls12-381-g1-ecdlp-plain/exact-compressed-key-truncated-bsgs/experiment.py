"""Exhaustive coverage check for the truncated baby-step/giant-step rectangle.

The scored algorithm runs in G1 with the 48-byte compressed encoding.  In an
additive cyclic group of order n the compressed key of a scalar is the scalar
itself, so the intersection logic is identical and can be checked exhaustively:
with m the least integer satisfying m^2 >= ceil(39n/100), the tables
{j : 0 <= j < m} and {(x - i*m) mod n : 0 <= i < m} intersect exactly when
x < m^2, and every emitted pair (i, j) satisfies i*m + j == x (mod n).

This is evidence about the coverage identity only.  It says nothing about
production point encodings, the 2^127.7-record sort, instruction counts or the
storage bound, which are accounted analytically in proof.md.
"""
import json
import math
import sys

order = json.load(sys.stdin)["order"]
if not (2 <= order <= 10000):
    raise SystemExit("order out of range")

target = -(-39 * order // 100)          # ceil(39*order/100)
m = math.isqrt(target)
if m * m < target:
    m += 1

recovered = 0
for x in range(order):
    giant = {(x - i * m) % order: i for i in range(m)}
    hit = None
    for j in range(m):
        if j in giant:
            hit = (giant[j], j)
            break
    expect_hit = x < m * m
    if (hit is not None) != expect_hit:
        raise AssertionError(
            "coverage mismatch at x=%d (hit=%r, expected=%r)" % (x, hit, expect_hit)
        )
    if hit is not None:
        i, j = hit
        if (i * m + j) % order != x:
            raise AssertionError("wrong candidate at x=%d: %r" % (x, hit))
        recovered += 1

assert recovered == min(order, m * m), (recovered, m * m, order)
print(json.dumps({
    "group": "additive cyclic group of the given order",
    "order": order,
    "m": m,
    "m_squared": m * m,
    "checked": order,
    "recovered": recovered,
    "exact_interval_coverage": True,
    "scope": "toy group only; does not validate production encodings or cost bounds",
}))
