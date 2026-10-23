"""Scope-limited analogue of the deterministic truncated BSGS of proof.md.

The production algorithm works in G1 = <P0>, a cyclic group of prime order r,
and keys its baby table on the compressed x-coordinate.  On an elliptic curve
x(P) = x(-P), so the key is exactly the quotient of the group by negation:
the pair {v, -v}.  This module reproduces the identical table construction,
open-addressing probe, sign branch, interval coverage and failure path in the
additive cyclic group Z_n, where the same quotient is {v, n-v}.  It is
evidence only at the scale it is run: it validates the recovery logic, the
exact covered set and the exact success count, not production costs, the
BLS12-381 arithmetic or the memory bound.
"""


def key_of(v, n):
    """Key and sign of the negation class of v in Z_n (n odd).

    Returns (representative, sign) where representative = min(v, n - v) and
    sign is 0 when v is the representative and 1 when -v is.
    """
    w = (-v) % n
    if v < w:
        return v, 0
    if w < v:
        return w, 1
    return 0, 0          # v == 0: the identity, the only self-opposite element


def build_table(n, m):
    """Records i = 0..m-1 and a hash index keyed by negation class.

    Keys are distinct for m < n/2 (Lemma 1 of proof.md), so the index is
    one-to-one and no probe ever compares two different records.
    """
    index = {}
    sign = [0] * m
    for i in range(m):
        k, s = key_of(i % n, n)
        index[k] = i
        sign[i] = s
    return index, sign


def solve(n, m, h, target, index=None, sign=None):
    """Return x with x = target in Z_n, or None when target is not covered."""
    if target % n == 0:
        return 0
    if index is None:
        index, sign = build_table(n, m)
    gamma = target % n
    step = (-m) % n
    for j in range(h):
        k, s = key_of(gamma, n)
        i = index.get(k)
        if i is not None:
            if i == 0:
                a = (j * m) % n
            elif sign[i] == s:
                a = (j * m + i) % n
            else:
                a = (j * m - i) % n
            if a == target % n:
                return a
        gamma = (gamma + step) % n
    return None


def covered_set(n, m, h):
    """Set of x in [0, n) on which solve must succeed, by Theorem 1."""
    return set(range(min(m * h, n))) | set(range(n - m + 1, n))


def selfcheck(n, m, h):
    """Exhaustive check: exact covered set and exact success count."""
    assert n % 2 == 1 and 0 < m < n // 2 and m * h < n - m
    index, sign = build_table(n, m)
    seen = set()
    for x in range(n):
        got = solve(n, m, h, x, index, sign)
        if got is None:
            continue
        assert got % n == x, (n, m, h, x, got)
        seen.add(x)
    expected = covered_set(n, m, h)
    assert seen == expected, (n, m, h, len(seen), len(expected))
    assert len(seen) == m * h + m - 1
    return len(seen), len(expected)


def main():
    for n, m, h in ((101, 7, 7), (1009, 20, 20), (7919, 61, 61)):
        got, want = selfcheck(n, m, h)
        print("n=%d m=%d h=%d covered=%d expected=%d" % (n, m, h, got, want))
    print("selfcheck ok")


if __name__ == "__main__":
    main()
