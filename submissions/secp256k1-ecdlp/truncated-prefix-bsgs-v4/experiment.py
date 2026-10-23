"""Production-arithmetic witness for the truncated-prefix BSGS ledger.

The experiment runs the identical secp256k1 affine code path used by the
construction - the same field modulus, the same generator, the same point
formulas, the same sort key and the same binary-search lookup - at reduced
table size, and it:

  * counts the field operations performed by every point addition and checks
    the observed maximum against the ledger's structural coefficients (one
    inversion, sixteen multiplications, thirty-two modular additions);
  * checks that the sorted baby keys are strictly increasing, which is what
    makes the binary search exact and collision free;
  * recovers scalars inside the covered prefix, including the endpoints, and
    confirms each by an independent [a]P = Q equality;
  * shows that a scalar just outside the prefix is not found within h giant
    steps, exhibiting the truncation boundary;
  * re-checks the exact full-size coverage and work ledger in integer
    arithmetic.

Scope: the reduced-size runs are evidence about the arithmetic and about the
list-generation and matching logic only. They do not execute production table
sizes and do not establish the production cost bound by measurement.
"""

import json
import sys
from bisect import bisect_left
from decimal import Decimal, ROUND_CEILING, getcontext
from math import isqrt

P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
G = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)

LEDGER_ADD = {"A": 32, "M": 16, "I": 1}

FIELD = {"A": 0, "M": 0, "I": 0}
ADD_MAX = {"A": 0, "M": 0, "I": 0}
ADD_CALLS = 0
MUL_MAX = 0


def fadd(a, b):
    FIELD["A"] += 1
    return (a + b) % P


def fsub(a, b):
    FIELD["A"] += 1
    return (a - b) % P


def fmul(a, b):
    FIELD["M"] += 1
    return (a * b) % P


def finv(a):
    FIELD["I"] += 1
    return pow(a, P - 2, P)


def neg(pt):
    return None if pt is None else (pt[0], (-pt[1]) % P)


def padd(p, q):
    global ADD_CALLS
    before = dict(FIELD)
    out = _padd(p, q)
    ADD_CALLS += 1
    for k in FIELD:
        d = FIELD[k] - before[k]
        if d > ADD_MAX[k]:
            ADD_MAX[k] = d
    return out


def _padd(p, q):
    if p is None:
        return q
    if q is None:
        return p
    x1, y1 = p
    x2, y2 = q
    if x1 == x2:
        if fadd(y1, y2) == 0:
            return None
        num = fmul(3, fmul(x1, x1))
        den = fmul(2, y1)
    else:
        num = fsub(y2, y1)
        den = fsub(x2, x1)
    lam = fmul(num, finv(den))
    x3 = fsub(fsub(fmul(lam, lam), x1), x2)
    y3 = fsub(fmul(lam, fsub(x1, x3)), y1)
    return (x3, y3)


def pmul(k, pt):
    global MUL_MAX
    before = ADD_CALLS
    out = None
    for bit in bin(k)[2:]:
        out = padd(out, out)
        if bit == "1":
            out = padd(out, pt)
    used = ADD_CALLS - before
    if used > MUL_MAX:
        MUL_MAX = used
    return out


def rec_key(rec):
    pt, _ = rec
    return (1, 0, 0) if pt is None else (0, pt[0], pt[1])


def build_baby(m):
    table = [(None, 0)]
    cur = None
    for j in range(1, m):
        cur = padd(cur, G)
        table.append((cur, j))
    return table, padd(cur, G)


def solve(q, keys, ordered, m, h, mp):
    step = neg(mp)
    gamma = q
    for i in range(h):
        key = rec_key((gamma, 0))
        pos = bisect_left(keys, key)
        if pos < len(keys) and keys[pos] == key:
            return (i * m + ordered[pos][1]) % N
        gamma = padd(gamma, step)
    return None


def log2_ceil(x):
    return (Decimal(x).ln() / Decimal(2).ln()).quantize(
        Decimal("1.000001"), rounding=ROUND_CEILING)


def full_ledger():
    getcontext().prec = 80
    s39 = -(-(39 * N) // 100)
    m = isqrt(s39 - 1) + 1
    h = -(-s39 // m)
    assert m * m >= s39 > (m - 1) * (m - 1)
    assert m * (m - 1) < s39
    assert h == m
    assert 100 * m * h >= 39 * N
    assert m * h < N

    a_k = 128 * 2
    m_k = 256 * (2 * a_k + 64) + 128 * 2
    i_k = 2 * 256 * (m_k + 64) + 128 * 2
    c_k = i_k + 16 * m_k + 32 * a_k + 4096
    assert (a_k, m_k, i_k, c_k) == (256, 147712, 75661568, 78037248)

    b, l, k, w = 256, 128, 3, 4
    b_k = 16 * k + 64
    d_k = 16 * w + 64
    v_k = (l + 1) * (b_k + 64) + d_k + 256

    setup = 1 << 24
    baby = m * (c_k + d_k)
    sort = m * l * (b_k + 2 * d_k + 128) + 4096 * l
    giant = h * (c_k + v_k)
    ver = 4 * b * c_k + (1 << 24)
    total = setup + baby + sort + giant + ver
    pre = setup + baby + sort
    mem = 2 * m * w * 32 + (1 << 20)

    time_log2 = log2_ceil(Decimal(total) / Decimal(c_k))
    pre_log2 = log2_ceil(Decimal(pre) / Decimal(c_k))
    mem_log2 = log2_ceil(mem)
    assert time_log2 == Decimal("128.321575")
    assert pre_log2 == Decimal("127.321949")
    assert mem_log2 == Decimal("135.320774")
    return {
        "m": m,
        "h": h,
        "covered_prefix_log2": str(log2_ceil(m * h)),
        "success_lower_bound": "39/100",
        "total_word_ops": str(total),
        "preprocessing_word_ops": str(pre),
        "peak_bytes": str(mem),
        "time_log2": str(time_log2),
        "preprocessing_log2": str(pre_log2),
        "memory_log2_bytes": str(mem_log2),
    }


def main():
    cfg = json.load(sys.stdin)
    m = int(cfg["m"])
    samples = int(cfg["samples"])
    assert 2 <= m <= 512
    assert 1 <= samples <= 32
    h = m
    span = m * h

    table, mp = build_baby(m)
    ordered = sorted(table, key=rec_key)
    keys = [rec_key(rec) for rec in ordered]
    for left, right in zip(keys, keys[1:]):
        assert left < right, "duplicate baby key"

    targets = [0, 1, 2, m - 1, m, m + 1, span - 1, span // 2]
    state = 0x9E3779B97F4A7C15
    while len(targets) < samples + 8:
        state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        targets.append(state % span)
    targets = sorted({t for t in targets if 0 <= t < span})

    recovered = 0
    for x in targets:
        q = pmul(x, G)
        a = solve(q, keys, ordered, m, h, mp)
        assert a is not None, "scalar inside the prefix was not recovered"
        assert a == x, "recovered scalar differs from the planted scalar"
        assert pmul(a, G) == q, "independent verification failed"
        recovered += 1

    outside = span
    assert outside < N
    missed = solve(pmul(outside, G), keys, ordered, m, h, mp) is None
    assert missed, "scalar outside the prefix was recovered"

    for op, bound in LEDGER_ADD.items():
        assert ADD_MAX[op] <= bound, "point addition exceeded the C coefficient"
    assert MUL_MAX <= 2 * 256, "scalar multiplication exceeded 2*b point additions"

    out = {
        "curve": "secp256k1",
        "arithmetic": "affine, production field modulus and generator",
        "m": m,
        "h": h,
        "covered_prefix": span,
        "planted_scalars": len(targets),
        "recovered": recovered,
        "outside_prefix_not_found": missed,
        "point_additions": ADD_CALLS,
        "max_field_ops_per_point_add": dict(ADD_MAX),
        "ledger_coefficients_per_point_add": dict(LEDGER_ADD),
        "max_point_additions_per_scalar_mul": MUL_MAX,
        "scalar_mul_bound": 2 * 256,
        "sorted_keys_strictly_increasing": True,
        "full_size_ledger": full_ledger(),
        "scope": "reduced table size; arithmetic and list logic only, not the production cost bound",
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
