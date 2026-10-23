"""Exact arithmetic for the submitted bound.  Every printed integer is a
quantity that appears in proof.md and submission.json.

Run with ``python3 budget.py``; nothing here touches the network or the
filesystem.
"""

import math
from fractions import Fraction

R_ORDER = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
C = 2000                      # BLS12-381 G1 reference-operation price
SUCCESS_NUM, SUCCESS_DEN = 39, 100

RECORD_BYTES = 64             # 48-byte x-coordinate + 16-byte step index
SORT_INSTR_PER_RECORD = 240   # 12 radix passes x 20 instructions
MERGE_INSTR_PER_ELEMENT = 20
FIXED_GROUP_OPS = 1056        # [m]P0 (254) + two verifications (800) + slack
FIXED_INSTRUCTIONS = 2 ** 36 + 1024


def parameters():
    """Smallest balanced (m, h) with m*h >= ceil(39r/100)."""
    target = -(-(R_ORDER * SUCCESS_NUM) // SUCCESS_DEN)
    m = math.isqrt(target)
    while m * m < target:
        m += 1
    h = -(-target // m)
    return target, m, h


def ledger(m, h):
    H = m + h + FIXED_GROUP_OPS
    R = (SORT_INSTR_PER_RECORD + MERGE_INSTR_PER_ELEMENT) * (m + h) + FIXED_INSTRUCTIONS
    W = Fraction(C * H + R, C)
    return H, R, W


def memory_bytes(m, h):
    return RECORD_BYTES * (m + h + max(m, h)) + 2 ** 35 + 2 ** 20


def preprocessing(m):
    return Fraction(C * m + SORT_INSTR_PER_RECORD * m, C)


def main():
    target, m, h = parameters()
    H, R, W = ledger(m, h)
    log2_w = math.log2(W.numerator) - math.log2(W.denominator)
    print("ceil(39r/100) =", target)
    print("m             =", m)
    print("h             =", h)
    print("m*h           =", m * h)
    print("m*h >= target =", m * h >= target, " m*h < r =", m * h < R_ORDER)
    print("H             =", H)
    print("R             =", R)
    print("W             =", W.numerator, "/", W.denominator)
    print("log2(W)       =", repr(log2_w))
    print("memory bytes  =", memory_bytes(m, h))
    print("log2 memory   =", repr(math.log2(memory_bytes(m, h))))
    pre = preprocessing(m)
    print("preprocessing =", pre.numerator, "/", pre.denominator)
    print("log2 preproc  =", repr(math.log2(pre.numerator) - math.log2(pre.denominator)))


if __name__ == "__main__":
    main()
