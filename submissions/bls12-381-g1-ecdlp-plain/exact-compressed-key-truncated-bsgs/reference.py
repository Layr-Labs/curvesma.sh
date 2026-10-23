"""Reference implementation of exact compressed-key truncated BSGS on BLS12-381 G1.

The scored claim lives in proof.md; this file is a small-scale reference for the
recovery logic, the compressed-key intersection and the batched affine chains.
It is not a production implementation and does not exercise the production sort
widths (2^32 buckets, 128-bit counters, 2^127.7 records), which are accounted
analytically in proof.md.

Exactness in one line: the 48-byte compressed encoding of a G1 point stores the
x-coordinate plus a sign bit and is injective on G1, so a key match is point
equality, so a cross-table match at (i, j) is exactly T = (i*m + j)*P0.  With
m*h < r the map (i, j) -> i*m + j is injective on [0, m) x [0, m-1] and its
image is the integer interval [0, m^2); the recoverable set is exactly that
interval, the success probability is exactly m^2/r, and at most one candidate is
ever emitted.
"""

from __future__ import annotations

# Order of the BLS12-381 G1 group.
BLS12_381_R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001

KEY_BYTES = 48          # canonical compressed G1 encoding
INDEX_BYTES = 16        # 128-bit record index; production N = 2m needs 128 bits
RADIX_BITS = 32         # production digit width; 48 bytes = 12 passes
RADIX_PASSES = (KEY_BYTES * 8) // RADIX_BITS


def select_parameters(r: int = BLS12_381_R, num: int = 39, den: int = 100):
    """Least m with m^2 >= ceil(num*r/den), and h = m - 1.

    The exact success probability m^2/r is then at least num/den, while
    m*h = m^2 - m < r keeps the rectangle injective.
    """
    target = -(-num * r // den)
    lo, hi = 1, 1 << ((target.bit_length() + 2) // 2 + 1)
    while lo < hi:
        mid = (lo + hi) // 2
        if mid * mid >= target:
            hi = mid
        else:
            lo = mid + 1
    m = lo
    assert m * m >= target and m * (m - 1) < r
    return m, m - 1


def batch_affine_chain(group, start, step, count, block=1 << 20):
    """start + k*step for k = 0..count-1 via batched affine additions.

    Montgomery's trick: accumulate the denominators of one block, invert their
    product once, then back-substitute.  The production schedule uses
    block = 2^20, which amortises one ~218,693-instruction inversion over 2^20
    additions (0.209 instructions per addition) inside the 3M batch allowance of
    the 6M reference price.  The reference keeps the same construction with a
    small block.
    """
    out = [start]
    point = start
    remaining = count - 1
    while remaining > 0:
        n = min(block, remaining)
        # Denominators of the n chained additions point + k*step.
        denoms = []
        cur = point
        for _ in range(n):
            denoms.append(group.add_denominator(cur, step))
            cur = group.add_numerator_point(cur, step)
        # Prefix products, one inversion, back-substitution.
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
    """LSD radix sort of (48-byte key, index) records on the key.

    Production uses 2^32 buckets with 128-bit counters and 128-bit scatter
    offsets (see the ledger in proof.md).  At reference scale the counters are
    sparse Python integers, which is the same algorithm with the same ordering.
    """
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


def solve(group, P0, T, m, r=BLS12_381_R):
    """Return x with T = x*P0 when x < m^2, else None.

    Both tables are keyed on the compressed encoding, so no sign branch exists.
    At most one cross-table pair is emitted, and it is verified once.
    """
    assert m * (m - 1) < r
    step = group.mul(P0, m)

    baby = batch_affine_chain(group, group.identity(), P0, m)
    giant = batch_affine_chain(group, T, group.neg(step), m)

    records = [(group.compress(p), (0 << 128) | j) for j, p in enumerate(baby)]
    records += [(group.compress(p), (1 << 128) | i) for i, p in enumerate(giant)]
    records = radix_sort_records(records)

    candidates = []
    for left, right in zip(records, records[1:]):
        if left[0] == right[0] and (left[1] >> 128) != (right[1] >> 128):
            candidates.append((left[1], right[1]))
    assert len(candidates) <= 1, "rectangle injectivity forbids two matches"

    for a, b in candidates:
        i = a & ((1 << 128) - 1) if (a >> 128) == 1 else b & ((1 << 128) - 1)
        j = b & ((1 << 128) - 1) if (a >> 128) == 0 else a & ((1 << 128) - 1)
        x = i * m + j
        if x < r and group.mul(P0, x) == T:
            return x
    return None


if __name__ == "__main__":
    m, h = select_parameters()
    print("r =", BLS12_381_R)
    print("m =", m)
    print("h =", h)
    print("m^2 =", m * m)
    print("success = m^2/r, at least 39/100:", m * m * 100 >= 39 * BLS12_381_R)
