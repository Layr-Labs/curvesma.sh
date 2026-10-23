"""Reference implementation: deterministic truncated baby-step/giant-step on
BLS12-381 G1 with a sorted-array intersection.

Track: bls12-381-g1-ecdlp-plain.  Instance: T = x*P0 with 0 <= x < r, P0 the
fixed generator of the prime-order group G1 (order r, cofactor irrelevant).

The algorithm is deterministic.  Both the baby array and the giant array are
sorted by affine x-coordinate and then merge-intersected; the recoverable
scalars are exactly the rectangle [0, m) x [0, h] inside Z/r, so the success
probability is the exact rational m*h/r and no record is dropped, evicted or
left unindexed.  Every reported scalar is confirmed by a full scalar
multiplication, so a false x-coordinate collision can never produce a wrong
answer.
"""

from bisect import bisect_left
from fractions import Fraction

P_MOD = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R_ORDER = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
CURVE_B = 4

X_BYTES = 48          # a 381-bit Fp coefficient
INDEX_BYTES = 16      # a step index below 2^128
RECORD_BYTES = X_BYTES + INDEX_BYTES


def _inv(a):
    return pow(a % P_MOD, P_MOD - 2, P_MOD)


def point_add(a, b):
    """Affine addition on y^2 = x^3 + 4 over Fp; None is the identity."""
    if a is None:
        return b
    if b is None:
        return a
    x1, y1 = a
    x2, y2 = b
    if x1 == x2:
        if (y1 + y2) % P_MOD == 0:
            return None
        lam = (3 * x1 * x1 % P_MOD) * _inv(2 * y1) % P_MOD
    else:
        lam = ((y2 - y1) % P_MOD) * _inv(x2 - x1) % P_MOD
    x3 = (lam * lam - x1 - x2) % P_MOD
    y3 = (lam * (x1 - x3) - y1) % P_MOD
    return (x3, y3)


def point_neg(a):
    if a is None:
        return None
    return (a[0], (-a[1]) % P_MOD)


def point_mul(k, point):
    acc = None
    while k > 0:
        if k & 1:
            acc = point_add(acc, point)
        point = point_add(point, point)
        k >>= 1
    return acc


def _isqrt(n):
    if n < 0:
        raise ValueError("isqrt of negative")
    if n == 0:
        return 0
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


def parameters(success_num=39, success_den=100, r=R_ORDER):
    """Smallest balanced (m, h) with m*h >= ceil(r * success_num/success_den).

    Balanced means m = ceil(sqrt(target)) and h = ceil(target/m), which
    minimises m + h for a fixed product.  The product is taken exactly.
    """
    target = -(-(r * success_num) // success_den)
    m = _isqrt(target)
    while m * m < target:
        m += 1
    h = -(-target // m)
    return m, h


def baby_records(P, m):
    """Sorted (x, j) records for j in [1, m); j = 0 is the identity."""
    out = []
    acc = None
    for j in range(1, m):
        acc = point_add(acc, P)
        out.append((acc[0], j))
    out.sort()
    return out


def giant_records(Q, P, m, h):
    """Sorted (x, i) records for i in [0, h) with G_i = Q - i*[m]P non-identity.

    The identity giants are returned separately as the list of indices i for
    which x == i*m exactly.
    """
    step = point_mul(m, P)
    out = []
    identities = []
    gamma = Q
    for i in range(h):
        if gamma is None:
            identities.append(i)
        else:
            out.append((gamma[0], i))
        gamma = point_add(gamma, point_neg(step))
    out.sort()
    return out, identities


def intersect(baby, giant):
    """Yield (i, j) pairs whose x-coordinates agree, in sorted order."""
    a = b = 0
    nb, ng = len(baby), len(giant)
    while a < nb and b < ng:
        xa, ja = baby[a]
        xg, ib = giant[b]
        if xa == xg:
            # Advance over every record sharing this x-coordinate on both
            # sides; duplicates are rare but must not be skipped.
            a_end = a
            while a_end < nb and baby[a_end][0] == xa:
                a_end += 1
            b_end = b
            while b_end < ng and giant[b_end][0] == xg:
                b_end += 1
            for k in range(a, a_end):
                for l in range(b, b_end):
                    yield giant[l][1], baby[k][1]
            a, b = a_end, b_end
        elif xa < xg:
            a += 1
        else:
            b += 1


def recover(Q, P, m, h):
    """Return x with Q = x*P, or None when x is outside the covered rectangle."""
    if Q is None:
        return 0
    baby = baby_records(P, m)
    giant, identities = giant_records(Q, P, m, h)
    for i in identities:
        cand = (i * m) % R_ORDER
        if point_mul(cand, P) == Q:
            return cand
    for i, j in intersect(baby, giant):
        for cand in ((i * m + j) % R_ORDER, (i * m - j) % R_ORDER):
            if point_mul(cand, P) == Q:
                return cand
    return None


def claim(m=None, h=None, r=R_ORDER):
    """Exact claim numbers for the balanced parameters at the 0.39 floor."""
    if m is None or h is None:
        m, h = parameters(r=r)
    covered = min(r, m * h)
    return {
        "m": m,
        "h": h,
        "covered": covered,
        "success": Fraction(covered, r),
        "memory_bytes": RECORD_BYTES * (m + h + max(m, h)) + 2 ** 35 + 2 ** 20,
        "baby_group_ops": m,
        "giant_group_ops": h,
    }


if __name__ == "__main__":
    c = claim()
    print("m =", c["m"])
    print("h =", c["h"])
    print("success =", c["success"])
    print("memory_bytes =", c["memory_bytes"])
