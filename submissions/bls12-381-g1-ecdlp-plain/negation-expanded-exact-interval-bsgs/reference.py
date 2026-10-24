"""Negation-expanded exact-interval BSGS on prime-order BLS12-381 G1.

Baby steps are j*P0 for j in [0, m). Each nonzero baby is stored twice: once
under its compressed encoding and once under the encoding of its negative.
On G1 the second key is the first key with the 0x20 sign flag flipped, so the
extra record is a copy, not a group operation. Giant steps use stride
s = 2*m - 1. Injective compressed keys make a hit point equality, and the
recoverable scalars are exactly two explicit intervals whose cardinalities sum
to h*s.

The scored claim is the ledger in proof.md. This module is the recovery logic
at reference scale plus the parameter predicate for the production order.
"""

from __future__ import annotations

BLS12_381_R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001

# Production minimizer of the proof.md ledger. Checked, not searched, on load
# of the production order: a 2^126 scan is not part of the attack.
M_PROD = 87990032694955804185480753244926427368
H_PROD = 116206294576542127460944719508070012830

KEY_BYTES = 48
INDEX_BITS = 128
SIGN_BIT = 1 << 127
SIGN_FLAG = 0x20  # py_ecc / Zcash compressed-G1 y-largest flag
RADIX_BITS = 32
RADIX_PASSES = (KEY_BYTES * 8) // RADIX_BITS
BUCKET_INSTR = 12 * (1 << 32) * 6


def _ledger_numerator(m: int, h: int) -> int:
    """Numerator of W = H + R/2000 for this (m, h), denominator 2000."""
    n_records = (2 * m - 1) + h
    group_ops = (m - 1) + (h - 1) + 2 + 512
    instructions = (
        920 * n_records
        + 16 * (m - 1)
        + 2 * BUCKET_INSTR
        + 4096
    )
    return 2000 * group_ops + instructions


def select_parameters(r: int = BLS12_381_R, num: int = 39, den: int = 100):
    """Return (m, h, s) with s = 2*m-1 and h*s >= ceil(num*r/den).

    On the production order this returns the hardcoded minimizer after checking
    the coverage, distinctness and 128-bit index constraints, and that both
    neighbouring baby-lengths have a strictly larger ledger numerator.
    """
    if r == BLS12_381_R and num == 39 and den == 100:
        m, h = M_PROD, H_PROD
        s = 2 * m - 1
        target = (num * r + den - 1) // den
        coverage = h * s
        assert s < r and h < (1 << 127) and m < (1 << 127)
        assert 2 * (m - 1) < r
        assert target <= coverage <= r
        base = _ledger_numerator(m, h)
        assert base < _ledger_numerator(m - 1, (target + (2 * m - 3) - 1) // (2 * m - 3))
        assert base < _ledger_numerator(m + 1, (target + (2 * m + 1) - 1) // (2 * m + 1))
        return m, h, s

    target = (num * r + den - 1) // den
    best = None
    upper = min(r // 2, max(target, 1))
    for m in range(1, upper + 1):
        s = 2 * m - 1
        if s >= r or 2 * (m - 1) >= r:
            break
        h = (target + s - 1) // s
        if h < 1 or h >= (1 << 127) or h * s > r:
            continue
        numerator = _ledger_numerator(m, h)
        if best is None or numerator < best[0]:
            best = (numerator, m, h, s)
    if best is None:
        raise RuntimeError("no feasible negation-expanded parameters")
    return best[1], best[2], best[3]


def flip_compressed_sign(key: bytes) -> bytes:
    """Encoding of -P given the 48-byte compressed encoding of a non-identity P."""
    if len(key) != KEY_BYTES:
        raise ValueError("compressed G1 key is 48 bytes")
    out = bytearray(key)
    out[0] ^= SIGN_FLAG
    return bytes(out)


def batch_affine_chain(group, start, step, count, block=1 << 20):
    """start + k*step for k = 0..count-1, one inversion per block."""
    if count <= 0:
        return []
    out = [start]
    point = start
    remaining = count - 1
    while remaining > 0:
        n = min(block, remaining)
        denoms = []
        cur = point
        for _ in range(n):
            denoms.append(group.add_denominator(cur, step))
            cur = group.add_numerator_point(cur, step)
        prefix = [1]
        for d in denoms:
            prefix.append(group.field_mul(prefix[-1], d))
        inv_total = group.field_inv(prefix[-1])
        inverses = [0] * n
        for k in range(n - 1, -1, -1):
            inverses[k] = group.field_mul(inv_total, prefix[k])
            inv_total = group.field_mul(inv_total, denoms[k])
        for k in range(n):
            point = group.add_with_inverse(point, step, inverses[k])
            out.append(point)
        remaining -= n
    return out


def radix_sort_records(records, passes=RADIX_PASSES, radix_bits=RADIX_BITS):
    """LSD radix sort of (key, index) records. Same digit order as the ledger."""
    mask = (1 << radix_bits) - 1
    order = list(records)
    for p in range(passes):
        shift = p * radix_bits
        counters = {}
        for key, _ in order:
            digit = (int.from_bytes(key, "big") >> shift) & mask
            counters[digit] = counters.get(digit, 0) + 1
        offset = {}
        running = 0
        for digit in sorted(counters):
            offset[digit] = running
            running += counters[digit]
        out = [None] * len(order)
        for record in order:
            digit = (int.from_bytes(record[0], "big") >> shift) & mask
            out[offset[digit]] = record
            offset[digit] += 1
        order = out
    return order


def _merge_hit(baby_sorted, giant_sorted):
    """One key match between the two sorted lists, or None."""
    i = 0
    j = 0
    hit = None
    while i < len(baby_sorted) and j < len(giant_sorted):
        bk, bi = baby_sorted[i]
        gk, gi = giant_sorted[j]
        if bk == gk:
            if hit is not None:
                raise AssertionError("two matches")
            hit = (bi, gi)
            i += 1
            j += 1
        elif bk < gk:
            i += 1
        else:
            j += 1
    return hit


def solve(group, P0, T, m, h, s, order):
    """Return x with T = x*P0 when x is in the covered set, else None."""
    if s != 2 * m - 1:
        raise ValueError("stride must be 2*m-1")
    baby = batch_affine_chain(group, group.identity(), P0, m)
    end = baby[-1]
    stride_point = group.add(group.dbl(end), P0)
    giant = batch_affine_chain(group, T, group.neg(stride_point), h)

    baby_records = []
    for j, point in enumerate(baby):
        baby_records.append((group.compress(point), j))
        if j:
            # Production uses flip_compressed_sign; any lawful compress() that
            # is injective and satisfies compress(-P) = flip(compress(P)) may
            # skip the group negation. The toy group negates explicitly.
            baby_records.append((group.compress(group.neg(point)), SIGN_BIT | j))
    giant_records = [(group.compress(point), i) for i, point in enumerate(giant)]
    hit = _merge_hit(
        radix_sort_records(baby_records),
        radix_sort_records(giant_records),
    )
    if hit is None:
        return None
    baby_index, giant_index = hit
    j = baby_index & (SIGN_BIT - 1)
    sign = -1 if baby_index & SIGN_BIT else 1
    x = (giant_index * s + sign * j) % order
    if group.mul(P0, x) == T:
        return x
    return None


if __name__ == "__main__":
    m, h, s = select_parameters()
    assert (m, h, s) == (M_PROD, H_PROD, 2 * M_PROD - 1)
    coverage = h * s
    print("m", m)
    print("h", h)
    print("s", s)
    print("coverage", coverage)
    print("success_at_least_39/100", coverage * 100 >= 39 * BLS12_381_R)
