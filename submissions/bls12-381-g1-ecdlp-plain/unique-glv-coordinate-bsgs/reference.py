"""Unique GLV-coordinate rectangle on plain BLS12-381 G1.

Baby points are a*P0 for a in [0, m). Giant points are T - b*[lambda]P0 for
b in [0, h). Each point is stored once, under the minimum compressed encoding
of its order-6 orbit. The scored set is the four-unit rectangle

    S = { b*lambda + sigma*a  :  0 <= b < h, 0 <= a < m,
          sigma in {1, -1, 1+lambda, -(1+lambda)},
          sigma = 1 whenever a = 0,
          and the -(1+lambda) family is taken only for b < a }

with h <= m-1. Its cardinality is exactly h*(4*m - h - 2). The parameters
below are a feasible pair inside the disjointness constraints
2*m - 1 <= lambda and (m-1)+(h-1) <= lambda-1, with |S| >= ceil(39*r/100).
proof.md is the argument. This module only pins the integers.
"""

from __future__ import annotations

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
LAMBDA = 0xAC45A4010001A40200000000FFFFFFFF
BETA = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAC
BETA2 = 0x5F19672FDF76CE51BA69C6076A0F77EADDB3A93BE6F89688DE17D813620A00022E01FFFFFFFEFFFE

# h/m = 2/3 band. Both are 126-bit, h <= m-1, 2*m-1 <= LAMBDA.
M_PROD = 95929641369258993122728749523035932376
H_PROD = 63953094246172662083740702189702033217

PER_POINT = 1750
PER_IDENTITY = 920
VERIFY_OPS = 2048
BUCKET_INSTR = 12 * (1 << 32) * 6
FIXED_INSTR = 8192
C = 2000


def coverage(m: int, h: int) -> int:
    """|S| = h*(4*m - h - 2) for 1 <= h <= m-1."""
    if m < 2 or h < 1 or h > m - 1:
        raise ValueError(f"rectangle needs 1 <= h <= m-1, got m={m} h={h}")
    return h * (4 * m - h - 2)


def work_fraction(m: int = M_PROD, h: int = H_PROD) -> tuple[int, int]:
    """W = H + R/C as a reduced positive fraction (numerator, denominator)."""
    group_ops = (m - 1) + (h - 1) + VERIFY_OPS
    instructions = (
        PER_IDENTITY
        + PER_POINT * (m - 1)
        + PER_POINT * h
        + 2 * BUCKET_INSTR
        + FIXED_INSTR
    )
    num = C * group_ops + instructions
    den = C
    g = _gcd(num, den)
    return num // g, den // g


def memory_bytes(m: int = M_PROD, h: int = H_PROD) -> int:
    """64*N + 16*max(m,h) + 2^36 counters + 2^20 fixed."""
    n = m + h
    return 64 * n + 16 * max(m, h) + (1 << 36) + (1 << 20)


def preprocessing_fraction(m: int = M_PROD) -> tuple[int, int]:
    """Baby chain, baby canonicalisation, baby sort, one bucket pass. In W."""
    per_nonidentity = 766 + 64 + 32 + 840  # 1702; merge is online
    instructions = 872 + per_nonidentity * (m - 1) + BUCKET_INSTR
    num = C * (m - 1) + instructions
    den = C
    g = _gcd(num, den)
    return num // g, den // g


def check_production() -> None:
    lam = LAMBDA
    if lam * lam + lam + 1 != R:
        raise AssertionError("lambda^2 + lambda + 1 != r")
    m, h = M_PROD, H_PROD
    if not (1 <= h <= m - 1):
        raise AssertionError("h not in 1..m-1")
    if 2 * m - 1 > lam:
        raise AssertionError("signed A-range meets")
    if (m - 1) + (h - 1) > lam - 1:
        raise AssertionError("coefficient sum wraps")
    target = (39 * R + 99) // 100
    got = coverage(m, h)
    if got < target:
        raise AssertionError("coverage below 39/100")
    num, den = work_fraction(m, h)
    if num != 74945032319733588378032555490423231295467 or den != 250:
        raise AssertionError(f"unexpected work fraction {num}/{den}")


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return a


if __name__ == "__main__":
    check_production()
