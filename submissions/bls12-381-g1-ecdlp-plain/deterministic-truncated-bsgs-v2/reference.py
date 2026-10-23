"""Scope-limited analogue of the truncated BSGS of proof.md.

The production algorithm works in G1 = <P0>, a cyclic group of prime order r,
and keys its baby table on the compressed x-coordinate.  On an elliptic curve
x(P) = x(-P), so the key is exactly the quotient of the group by negation:
the pair {v, -v}.  The key space is far larger than the slot space, so the
128-bit reduction of the key is what makes slots collide.  This module
reproduces the identical table construction, the seeded multiply-shift probe
family, the eight-probe window with its drop rule, the sign branch, the
interval coverage and the failure path in the additive cyclic group Z_n, with a
separate slot count S standing in for 2^128.  It is evidence only at the scale
it is run: it validates the recovery logic, the exact indexed-set coverage and
the drop path, not production costs, the BLS12-381 arithmetic or the memory
bound.
"""

M64 = (1 << 64) - 1


def key_of(v, n):
    """Key and sign of the negation class of v in Z_n (n odd).

    The key is a large label of the class representative, standing in for the
    381-bit compressed x-coordinate; the sign is 0 when v is the representative
    and 1 when -v is.  Key 0 is the identity, the only self-opposite element.
    """
    w = (-v) % n
    if v < w:
        rep, sign = v, 0
    elif w < v:
        rep, sign = w, 1
    else:
        rep, sign = 0, 0
    if rep == 0:
        return 0, 0
    return _mix(0x5EED5EED5EED5EED, rep, 0x9E3779B97F4A7C15) | 1, sign


def _mix(seed, t, salt):
    """Deterministic stand-in for one word of the setup seed."""
    x = (seed * 0x9E3779B97F4A7C15 + t * 0xBF58476D1CE4E5B9 + salt) & M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)


def probe_family(seed, cap, slots):
    """cap independent multiply-shift probes s_t(k) = (a_t*k + b_t) mod slots.

    a_t is odd, exactly as in the production family, so the probe is a
    bijection of the reduced key and collisions come only from the reduction
    of the key modulo the slot count.
    """
    family = []
    for t in range(cap):
        a = 1 + 2 * (_mix(seed, t, 0xA5A5A5A5A5A5A5A5) % (slots // 2))
        b = _mix(seed, t, 0x5A5A5A5A5A5A5A5A) % slots
        family.append((a, b))
    return family


def build_table(n, m, family, slots):
    """Index i = 0..m-1 in the first free slot of its capped probe window.

    Keys are distinct for m < n/2 (Lemma 1), so no slot ever holds two records
    with the same key.  A record whose whole window is occupied is left
    unindexed.
    """
    index = {}
    sign = [0] * m
    placed = []
    for i in range(m):
        k, s = key_of(i % n, n)
        sign[i] = s
        kk = k % slots
        for a, b in family:
            slot = (a * kk + b) % slots
            if slot not in index:
                index[slot] = (k, i)
                placed.append(i)
                break
    return index, sign, placed


def solve(n, m, h, target, index, sign, family, slots):
    """Return x with x = target in Z_n, or None when target is not covered."""
    if target % n == 0:
        return 0
    gamma = target % n
    step = (-m) % n
    for j in range(h):
        k, s = key_of(gamma, n)
        kk = k % slots
        for a, b in family:
            slot = (a * kk + b) % slots
            e = index.get(slot)
            if e is None:
                break                       # this j misses
            ek, i = e
            if ek == k:
                if i == 0:
                    cand = (j * m) % n
                elif sign[i] == s:
                    cand = (j * m + i) % n
                else:
                    cand = (j * m - i) % n
                if cand == target % n:
                    return cand
                break                       # unreachable, defensive
        gamma = (gamma + step) % n
    return None


def covered_set(n, m, h, placed):
    """Residues on which solve must succeed, from Theorem 1."""
    indexed = set(placed)
    seen = {0}                              # the identity pre-check
    for i in indexed:
        for j in range(h):
            seen.add(j * m + i)
        if i:
            seen.add(n - i)
    return seen


def selfcheck(n, m, h, slots, cap, seed):
    """Exhaustive check: exact coverage of the indexed set and no wrong answer."""
    assert n % 2 == 1 and 0 < m < n // 2 and m * h + m < n and 2 * m <= slots
    family = probe_family(seed, cap, slots)
    index, sign, placed = build_table(n, m, family, slots)
    seen = set()
    for x in range(n):
        got = solve(n, m, h, x, index, sign, family, slots)
        if got is None:
            continue
        assert got % n == x, (n, m, h, x, got)
        seen.add(x)
    expected = covered_set(n, m, h, placed)
    assert expected <= seen, (n, m, h, len(seen), len(expected))
    a = len(placed)
    if a == m:
        assert seen == expected, (n, m, h, len(seen), len(expected))
    assert len(seen) >= a * h + a - 1, (n, m, h, len(seen), a)
    return len(seen), a, m - a


def main():
    cases = (
        (101, 7, 7, 16, 8),
        (1009, 20, 20, 64, 8),
        (7919, 61, 61, 128, 8),
        (1009, 200, 4, 1024, 2),    # short window: exercises the drop path
    )
    for n, m, h, slots, cap in cases:
        got, a, dropped = selfcheck(n, m, h, slots, cap, 0x1234567890ABCDEF)
        print("n=%d m=%d h=%d slots=%d cap=%d load=%.3f covered=%d indexed=%d dropped=%d"
              % (n, m, h, slots, cap, m / slots, got, a, dropped))
    print("selfcheck ok")


if __name__ == "__main__":
    main()
