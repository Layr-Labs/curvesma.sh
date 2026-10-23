"""selector.py -- v17: the ONE implementation of the cell map and the restart draw.

Production and toy share this module; the ledger prices exactly this function.

Spec (argument.md section 2):

    K      = min(y, P - y)                       # class key of a six-class orbit
    cell   = (K * R) >> (P.bit_length() - 1)     # R = 2048 table rows

This is a POWER-OF-TWO quotient: it divides K*R by 2^(bitlen P - 1), not by P.
It is deliberately NOT floor(K*R/P) -- the denominator is the power of two that
bounds the class key, which is what makes the map a shift (4 instructions) rather
than a 6-limb reciprocal multiply (~134).  The two maps differ: for the 20-bit
fixture with K = 2045, R = 4096 the shift gives 15 while floor(K*R/P) gives 13.

What makes the shift map sound:

  * K <= P//2 < 2^(bitlen P - 1), so cell < R always: no mask is needed and no
    column can overflow.
  * the reachable K-range [0, P//2] is a prefix of [0, 2^(bitlen P - 1)), which
    the shift cuts into R equal cells of 2^(bitlen P - 1)/R keys each.  Every
    reachable column therefore has exactly that many preimages except the last,
    partial one, so the column law is uniform up to that single cut-off, and the
    number of reachable columns is printed by `reachable_columns` and used as R in
    the self-map correction N' = N*R/(R-1).  No column is counted twice.

`reachable_columns(R, P)` is the number of table rows that can actually be
reached; v17 reports it (1665 of 2048 for BLS12-381 G1) and uses it as R in the
self-map correction N' = N*R/(R-1), so no row is counted twice.

Restarts: `RestartStream` is the declared RNG for infinity restarts -- a SHA-256
counter stream from a fixed seed, so every restart draws a FRESH pair (a, c)
rather than incrementing the old one, and the draw is reproducible from the seed
and the draw index.
"""
import hashlib


def cell_exact(K, R, P):
    """The spec, verbatim: (K*R) >> (bitlen(P) - 1), exact integers.

    Independent reference for the shared binding below; the certificate re-derives
    the same expression inline so a broken alias cannot certify itself.
    """
    if not (0 <= K <= P // 2):
        raise ValueError("class key out of range (must be in [0, P//2])")
    return (K * R) >> (P.bit_length() - 1)


def cell_reciprocal(K, R, P):
    """The same value by exact integer division instead of a shift.

    floor(K*R / 2^(bitlen P - 1)) is what the shift computes; this evaluates the
    quotient with // so a broken shift cannot hide.
    """
    return (K * R) // (1 << (P.bit_length() - 1))


# The production and the toy share this binding; nothing else may select a cell.
cell = cell_exact


def reachable_columns(R, P):
    """Number of rows a class key can reach: floor((P//2)*R / 2^(bitlen P-1)) + 1."""
    return ((P // 2) * R >> (P.bit_length() - 1)) + 1


def columns_are_equally_weighted(K_lo, K_hi, R, P):
    """Exact uniformity check on a K-range: count preimages column by column."""
    counts = {}
    for K in range(K_lo, K_hi):
        c = cell(K, R, P)
        counts[c] = counts.get(c, 0) + 1
    inner = [v for k, v in counts.items() if k not in
             (min(counts), max(counts))]
    return (not inner) or (max(inner) - min(inner) <= 1)


class RestartStream:
    """Declared RNG stream for infinity restarts (fresh pairs, logged)."""

    def __init__(self, seed, Q):
        if isinstance(seed, str):
            seed = seed.encode()
        self.seed = seed
        self.Q = Q
        self.chain = hashlib.sha256(b"plain-v17-restart:" + seed).digest()
        self.draws = 0
        self.log = []

    def next_pair(self):
        out = []
        for _ in range(2):
            self.chain = hashlib.sha256(self.chain).digest()
            out.append(int.from_bytes(self.chain[:16], "big") % self.Q)
        self.draws += 1
        a, c = out
        if self.draws <= 4:
            self.log.append([a, c])
        return a, c

    def digest(self):
        return hashlib.sha256(
            b"plain-v17-restart:" + self.seed + b":%d" % self.draws).hexdigest()[:16]
