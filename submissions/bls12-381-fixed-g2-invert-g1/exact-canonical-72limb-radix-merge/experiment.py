"""Exact 72-limb LSD radix merge on canonical Fp12 residues.

The production ledger lives in `ledger`. The Docker entry reads a small JSON
object and checks that same join on real Fp12 products. It does not price the
full attack and it does not use an additive toy group.
"""

from __future__ import annotations

import json
import sys
from array import array

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 20000
D = 1 << 16
PASSES = 288  # 72 limbs * 4 digits of 16 bits
N_LIMBS = 72

# Each touched u64 costs WORD instructions, which covers a two-limb address
# and the load or store. Ceilings sit above the schedule in proof.md.
WORD = 16
CLEAR = (2 * WORD + 8) * D + 64
PREFIX = (4 * WORD + 16) * D + 64
PASS_SETUP = 64
HIST = 7 * WORD + 32
SCAT = 9 * WORD + 40
PASS_PER_N = HIST + SCAT
MERGE_PER_CMP = (144 * WORD) + (72 * 4) + 32
BUILD_PER_RECORD = (2 * WORD) + 16
MISC = 1 << 20
ALLOWANCE = 1 << 19  # preprocessing half and online half
INVERSE_MULS = 512


def window_m() -> int:
    threshold = (391 * R + 999) // 1000
    m = isqrt(threshold)
    if m * m < threshold:
        m += 1
    if m * m < threshold or (m - 1) ** 2 >= threshold:
        raise AssertionError("window")
    if m * m > R:
        raise AssertionError("window exceeds r")
    return m


def isqrt(n: int) -> int:
    if n < 0:
        raise ValueError("negative")
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


def ledger() -> dict[str, int]:
    m = window_m()
    one_sort = PASSES * (CLEAR + PREFIX + PASS_SETUP) + PASSES * PASS_PER_N * m
    pre_h = m + INVERSE_MULS + ALLOWANCE
    online_h = (m - 1) + ALLOWANCE
    pre_r = one_sort + BUILD_PER_RECORD * m + (MISC // 2)
    online_r = one_sort + (2 * m) * MERGE_PER_CMP + BUILD_PER_RECORD * m + (MISC // 2)
    total_h = pre_h + online_h
    total_r = pre_r + online_r
    # bytes: two 72-limb keys, three 2-word index vectors, 128-bit counters, 1 MiB
    peak = 2 * m * (N_LIMBS * 8) + 3 * m * 16 + D * 16 + (1 << 20)
    return {
        "m": m,
        "threshold": (391 * R + 999) // 1000,
        "covered": m * m,
        "pre_h": pre_h,
        "online_h": online_h,
        "total_h": total_h,
        "pre_r": pre_r,
        "online_r": online_r,
        "total_r": total_r,
        "pre_num": pre_h * C + pre_r,
        "total_num": total_h * C + total_r,
        "den": C,
        "peak_bytes": peak,
        "advice_bytes": 0,
    }


def limbs_from_fp(x: int) -> list[int]:
    if not 0 <= x < P:
        raise ValueError("residue")
    return [(x >> (64 * i)) & ((1 << 64) - 1) for i in range(6)]


def pack_fp12(coeffs: list[int]) -> tuple[int, ...]:
    if len(coeffs) != 12:
        raise ValueError("degree")
    out: list[int] = []
    for c in coeffs:
        out.extend(limbs_from_fp(c % P))
    return tuple(out)


def digits_lsd(limbs: tuple[int, ...], bits: int) -> list[int]:
    """Digit 0 is least significant for lexicographic order on the limb tuple.

    Within a limb the low bits are less significant. Limbs are stored with
    limb 0 as the most significant tuple entry, matching Python's tuple order.
    Production pricing uses bits=16. The runner uses bits=8 so the same
    recurrence finishes inside the 20-second experiment limit.
    """
    if 64 % bits:
        raise ValueError("digit width")
    mask = (1 << bits) - 1
    digits: list[int] = []
    for limb in reversed(limbs):
        for shift in range(0, 64, bits):
            digits.append((limb >> shift) & mask)
    return digits


def radix_lsd(records: list[tuple[tuple[int, ...], int]], bits: int) -> list[tuple[tuple[int, ...], int]]:
    width = 1 << bits
    keyed = [(digits_lsd(limbs, bits), limbs, idx) for limbs, idx in records]
    n = len(keyed)
    passes = N_LIMBS * (64 // bits)
    # 64-bit counters are exact for the scaled n. Production n needs 128-bit
    # counters; the digit recurrence and the stable scatter are unchanged.
    zeros = bytes(width * 8)
    for d in range(passes):
        count = array("Q")
        count.frombytes(zeros)
        for digits, _, _ in keyed:
            count[digits[d]] += 1
        total = 0
        for i in range(width):
            c = count[i]
            count[i] = total
            total += c
        out: list = [None] * n
        for item in keyed:
            slot = count[item[0][d]]
            out[slot] = item
            count[item[0][d]] = slot + 1
        keyed = out
    ordered = [(limbs, idx) for _, limbs, idx in keyed]
    for i in range(1, n):
        if ordered[i - 1][0] > ordered[i][0]:
            raise AssertionError("radix order")
    return ordered


def merge_join(baby, giant, bits: int = 8):
    a = radix_lsd(baby, bits)
    b = radix_lsd(giant, bits)
    i = j = 0
    while i < len(a) and j < len(b):
        if a[i][0] == b[j][0]:
            return a[i][1], b[j][1]
        if a[i][0] < b[j][0]:
            i += 1
        else:
            j += 1
    return None


def fp_add(a: int, b: int) -> int:
    c = a + b
    return c - P if c >= P else c


def fp_sub(a: int, b: int) -> int:
    return a - b if a >= b else a + P - b


def fp_mul(a: int, b: int) -> int:
    return (a * b) % P


def fp2_add(a, b):
    return (fp_add(a[0], b[0]), fp_add(a[1], b[1]))


def fp2_sub(a, b):
    return (fp_sub(a[0], b[0]), fp_sub(a[1], b[1]))


def fp2_mul(a, b):
    a0b0 = fp_mul(a[0], b[0])
    a1b1 = fp_mul(a[1], b[1])
    cross = fp_mul(fp_add(a[0], a[1]), fp_add(b[0], b[1]))
    return (fp_sub(a0b0, a1b1), fp_sub(fp_sub(cross, a0b0), a1b1))


XI = (1, 1)


def fp6_mul(a, b):
    z0 = fp2_mul(a[0], b[0])
    z1 = fp2_add(fp2_mul(a[0], b[1]), fp2_mul(a[1], b[0]))
    z2 = fp2_add(fp2_add(fp2_mul(a[0], b[2]), fp2_mul(a[2], b[0])), fp2_mul(a[1], b[1]))
    z3 = fp2_add(fp2_mul(a[1], b[2]), fp2_mul(a[2], b[1]))
    z4 = fp2_mul(a[2], b[2])
    return (
        fp2_add(z0, fp2_mul(z3, XI)),
        fp2_add(z1, fp2_mul(z4, XI)),
        z2,
    )


def fp6_add(a, b):
    return tuple(fp2_add(a[i], b[i]) for i in range(3))


def fp6_sub(a, b):
    return tuple(fp2_sub(a[i], b[i]) for i in range(3))


def fp12_mul(a, b):
    a0b0 = fp6_mul(a[0], b[0])
    a1b1 = fp6_mul(a[1], b[1])
    cross = fp6_mul(fp6_add(a[0], a[1]), fp6_add(b[0], b[1]))
    c0, c1, c2 = a1b1
    v_a1b1 = (fp2_mul(c2, XI), c0, c1)
    low = fp6_add(a0b0, v_a1b1)
    high = fp6_sub(fp6_sub(cross, a0b0), a1b1)
    return (low, high)


def zero2():
    return (0, 0)


def fp12_one():
    z = zero2()
    return ((1, 0), z, z), (z, z, z)


def fp12_coeffs(e) -> list[int]:
    out = []
    for fp6 in e:
        for fp2 in fp6:
            out.append(fp2[0])
            out.append(fp2[1])
    return out


def embed(k: int):
    """β = k + u in Fp2, sitting in Fp12. u^2 = -1."""
    z = zero2()
    return ((k % P, 1), z, z), (z, z, z)


def embed_inv(k: int):
    den = fp_add(fp_mul(k % P, k % P), 1)
    den_inv = pow(den, P - 2, P)
    z = zero2()
    return ((fp_mul(k % P, den_inv), fp_sub(0, den_inv)), z, z), (z, z, z)


def tower_bsgs(span: int, base_k: int) -> dict:
    if not 2 <= span <= 24:
        raise ValueError("span")
    base = embed(base_k)
    invb = embed_inv(base_k)
    if pack_fp12(fp12_coeffs(fp12_mul(base, invb))) != pack_fp12(fp12_coeffs(fp12_one())):
        raise AssertionError("inverse")
    baby_elems = [fp12_one()]
    acc = fp12_one()
    for _ in range(1, span):
        acc = fp12_mul(acc, base)
        baby_elems.append(acc)
    baby = [(pack_fp12(fp12_coeffs(e)), j) for j, e in enumerate(baby_elems)]
    if len({limbs for limbs, _ in baby}) != span:
        raise AssertionError("baby collision")
    gstep_inv = fp12_one()
    for _ in range(span):
        gstep_inv = fp12_mul(gstep_inv, invb)
    checked = 0
    corners = ((0, 0), (span - 1, span - 1))
    for i, j in corners:
        secret = i * span + j
        target = fp12_one()
        for _ in range(secret):
            target = fp12_mul(target, base)
        giant = []
        acc = target
        giant.append((pack_fp12(fp12_coeffs(acc)), 0))
        for t in range(1, span):
            acc = fp12_mul(acc, gstep_inv)
            giant.append((pack_fp12(fp12_coeffs(acc)), t))
        found = merge_join(baby, giant)
        if found != (j, i):
            raise AssertionError(f"join {found} != {(j, i)} secret {secret}")
        if found[1] * span + found[0] != secret:
            raise AssertionError("rebuild")
        checked += 1
    # Planted full-width keys, including a difference only in limb 0 and only in limb 71.
    hi = pack_fp12([1] + [0] * 11)
    lo_coeffs = [0] * 12
    lo = list(pack_fp12(lo_coeffs))
    lo[71] = 1
    lo_t = tuple(lo)
    mid = pack_fp12([2] + [0] * 11)
    planted = merge_join([(mid, 4), (hi, 1), (lo_t, 2)], [(pack_fp12([3] + [0] * 11), 7), (lo_t, 9)])
    if planted != (2, 9):
        raise AssertionError(f"planted {planted}")
    if merge_join([(hi, 1)], [(mid, 2)]) is not None:
        raise AssertionError("false match")
    return {"span": span, "checked": checked, "base_k": base_k}


def main() -> None:
    spec = json.load(sys.stdin)
    report = tower_bsgs(int(spec["span"]), int(spec["base_k"]))
    numbers = ledger()
    # Keep stdout one JSON object. Huge integers are exact decimals.
    out = {
        "group": "Fp12 multiplication, canonical 72-limb keys, LSD radix merge",
        "checked_pairs": report["checked"],
        "span": report["span"],
        "base_k": report["base_k"],
        "false_match": False,
        "scope": "scaled Fp12 products and the exact join; not a production-size cost measurement",
        "m": str(numbers["m"]),
        "covered_over_r": [str(numbers["covered"]), str(R)],
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
