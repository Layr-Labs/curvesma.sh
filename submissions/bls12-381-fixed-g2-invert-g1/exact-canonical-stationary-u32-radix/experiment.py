"""Exact canonical join: stationary keys, 32-bit production digits.

Production prices 144 LSD passes of 32-bit digits (72 limbs x 2). The runner
executes the same counting-sort recurrence at 8-bit digits. A 2^32 counter
vector does not fit the 256 MiB runner, and a 2^16 histogram does not finish
inside the 20-second runner; the proof charges the 32-bit vector in full.
Keys stay put. Only the index permutation moves. Equality is 72-limb identity.
"""

from __future__ import annotations

import json
import sys
from array import array

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
C = 20000
WORD = 16
# Production digit. The runner never allocates this histogram.
PROD_BITS = 32
D = 1 << PROD_BITS
PASSES = 72 * (64 // PROD_BITS)  # 144
N_LIMBS = 72
CLEAR = (2 * WORD + 8) * D + 64
PREFIX = (4 * WORD + 16) * D + 64
PASS_SETUP = 64
HIST = 7 * WORD + 32
SCAT = 9 * WORD + 40
PASS_PER_N = HIST + SCAT
MERGE_PER_CMP = (144 * WORD) + (72 * 4) + 32
BUILD_PER_RECORD = (2 * WORD) + 16
MISC = 1 << 20
ALLOWANCE = 1 << 19
INVERSE_MULS = 512


def isqrt(n: int) -> int:
    if n < 0:
        raise ValueError("negative")
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


def window_m() -> int:
    threshold = (391 * R + 999) // 1000
    m = isqrt(threshold)
    if m * m < threshold:
        m += 1
    if m * m < threshold or (m - 1) ** 2 >= threshold or m * m > R:
        raise AssertionError("window")
    return m


def ledger() -> dict[str, int]:
    m = window_m()
    one_sort = PASSES * (CLEAR + PREFIX + PASS_SETUP) + PASSES * PASS_PER_N * m
    pre_h = m + INVERSE_MULS + ALLOWANCE
    online_h = (m - 1) + ALLOWANCE
    pre_r = one_sort + BUILD_PER_RECORD * m + (MISC // 2)
    online_r = one_sort + (2 * m) * MERGE_PER_CMP + BUILD_PER_RECORD * m + (MISC // 2)
    total_h = pre_h + online_h
    total_r = pre_r + online_r
    peak = 2 * m * (N_LIMBS * 8) + 3 * m * 16 + D * 16 + (1 << 20)
    return {
        "m": m,
        "threshold": (391 * R + 999) // 1000,
        "covered": m * m,
        "passes": PASSES,
        "digit_bits": PROD_BITS,
        "one_sort": one_sort,
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


def radix_order(keys: list[tuple[int, ...]], bits: int) -> list[int]:
    """LSD counting sort. `keys` is stationary. Returns a permutation of indices.

    Digit 0 is the low `bits` of limb 71. The last digit is the high `bits` of
    limb 0. That is lexicographic order on the 72-tuple. Production uses
    bits=32; this function is the same recurrence for any bits dividing 64.
    """
    if bits not in (8, 16, 32) or 64 % bits:
        raise ValueError("digit width")
    n = len(keys)
    order = list(range(n))
    width = 1 << bits
    mask = width - 1
    per_limb = 64 // bits
    passes = N_LIMBS * per_limb
    zeros = bytes(width * 8)
    for d in range(passes):
        limb = (N_LIMBS - 1) - (d // per_limb)
        shift = (d % per_limb) * bits
        count = array("Q")
        count.frombytes(zeros)
        digs = array("I")
        for idx in order:
            digit = (keys[idx][limb] >> shift) & mask
            digs.append(digit)
            count[digit] += 1
        total = 0
        for i in range(width):
            c = count[i]
            count[i] = total
            total += c
        nxt = [0] * n
        for pos, idx in enumerate(order):
            digit = digs[pos]
            slot = count[digit]
            count[digit] = slot + 1
            nxt[slot] = idx
        order = nxt
    for a, b in zip(order, order[1:]):
        if keys[a] > keys[b]:
            raise AssertionError("radix order")
    return order


def merge_join(
    baby_keys: list[tuple[int, ...]],
    baby_idx: list[int],
    giant_keys: list[tuple[int, ...]],
    giant_idx: list[int],
    bits: int,
) -> tuple[int, int] | None:
    bo = radix_order(baby_keys, bits)
    go = radix_order(giant_keys, bits)
    i = j = 0
    while i < len(bo) and j < len(go):
        bk = baby_keys[bo[i]]
        gk = giant_keys[go[j]]
        if bk == gk:
            return baby_idx[bo[i]], giant_idx[go[j]]
        if bk < gk:
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
    """beta = k + u in Fp2, sitting in Fp12. u^2 = -1."""
    z = zero2()
    return ((k % P, 1), z, z), (z, z, z)


def embed_inv(k: int):
    den = fp_add(fp_mul(k % P, k % P), 1)
    den_inv = pow(den, P - 2, P)
    z = zero2()
    return ((fp_mul(k % P, den_inv), fp_sub(0, den_inv)), z, z), (z, z, z)


def tower_bsgs(span: int, base_k: int, bits: int) -> dict:
    if not 2 <= span <= 8:
        raise ValueError("span")
    if bits != 8:
        raise ValueError("runner digit width")
    base = embed(base_k)
    invb = embed_inv(base_k)
    if pack_fp12(fp12_coeffs(fp12_mul(base, invb))) != pack_fp12(fp12_coeffs(fp12_one())):
        raise AssertionError("inverse")
    baby_elems = [fp12_one()]
    acc = fp12_one()
    for _ in range(1, span):
        acc = fp12_mul(acc, base)
        baby_elems.append(acc)
    baby_keys = [pack_fp12(fp12_coeffs(e)) for e in baby_elems]
    baby_idx = list(range(span))
    if len(set(baby_keys)) != span:
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
        giant_keys = []
        giant_idx = []
        acc = target
        giant_keys.append(pack_fp12(fp12_coeffs(acc)))
        giant_idx.append(0)
        for t in range(1, span):
            acc = fp12_mul(acc, gstep_inv)
            giant_keys.append(pack_fp12(fp12_coeffs(acc)))
            giant_idx.append(t)
        found = merge_join(baby_keys, baby_idx, giant_keys, giant_idx, bits)
        if found != (j, i):
            raise AssertionError(f"join {found} != {(j, i)} secret {secret}")
        if found[1] * span + found[0] != secret:
            raise AssertionError("rebuild")
        checked += 1
    hi = pack_fp12([1] + [0] * 11)
    lo = list(pack_fp12([0] * 12))
    lo[71] = 1
    lo_t = tuple(lo)
    mid = pack_fp12([2] + [0] * 11)
    other = pack_fp12([3] + [0] * 11)
    planted = merge_join(
        [mid, hi, lo_t],
        [4, 1, 2],
        [other, lo_t],
        [7, 9],
        bits,
    )
    if planted != (2, 9):
        raise AssertionError(f"planted {planted}")
    if merge_join([hi], [1], [mid], [2], bits) is not None:
        raise AssertionError("false match")
    return {"span": span, "checked": checked, "base_k": base_k, "bits": bits}


def main() -> None:
    spec = json.load(sys.stdin)
    report = tower_bsgs(int(spec["span"]), int(spec["base_k"]), int(spec["bits"]))
    numbers = ledger()
    out = {
        "group": "Fp12 multiplication, canonical 72-limb keys, stationary LSD radix merge",
        "checked_pairs": report["checked"],
        "span": report["span"],
        "base_k": report["base_k"],
        "runner_bits": report["bits"],
        "priced_digit_bits": numbers["digit_bits"],
        "priced_passes": numbers["passes"],
        "false_match": False,
        "scope": "scaled Fp12 products and the exact join; not a production-size cost measurement",
        "m": str(numbers["m"]),
        "covered_over_r": [str(numbers["covered"]), str(R)],
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
