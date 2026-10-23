#!/usr/bin/env python3
"""p256-ecdlp v9c SHARED PRODUCTION PROTOCOL (one place, imported by toy and ledger).

Everything that the toy and the certificate both depend on lives here and nowhere
else:

  KAPPA / FAIL_MODEL / S_ops()   the walk-budget formula (ONE constant, 1.10)
  dp predicate                   reads ONLY the low min(64, half) bits of the key
  router                         reads ONLY bits half..half+63 of the key
  node cap                       fixed a priori from the declared report mean

PRODUCTION SPEC (key width 256 bits, i.e. the x-coordinate of the class
representative):

  key K = x  in F_p  (the class {P,-P} shares x, so K is a function of the class)
  half  = 128
  DP predicate at width b:  mix64((K mod 2^64) ^ DPSALT) mod 2^b == 0
                            -> reads bits 0..63 of K only
  route:                    rw = mix64(((K >> 128) mod 2^64) ^ RSALT)
                            node = rw >> 48  in [0, 2^16)
                            -> reads bits 128..191 of K only

The two input fields {0..63} and {128..191} are DISJOINT, so the DP event never
constrains the routing field: for a key that is uniform on its field, the route
is exactly uniform over 2^16 conditional on the DP event.  mix64 is a bijection
on 64-bit words, so at b = 64 the predicate pins exactly one value of bits
0..63 and leaves bits 128..191 free.  (The rejected v9b spec fed the SAME
low-64-bit field to both the predicate and the router, so at b = 64 the route
was a constant: the collapse the reviewer identified.)
"""
from fractions import Fraction

C = 1000                          # instructions per reference unit (borrowed)
KAPPA = Fraction(11, 10)          # 1.10 walk-quality allowance, measured in the toy
FAIL_MODEL = Fraction(60, 100)    # 0.60 modelled per-walk failure probability
DP_BITS_PROD = 64                 # production DP predicate width
ROUTE_BITS = 16                  # production routing field width
NODES_PROD = 1 << 16              # 2^16 storage nodes
CAP_RECORDS_PER_NODE = 1 << 50    # declared fixed per-node cap
RECORD_BYTES = 64
MASK64 = (1 << 64) - 1
DPSALT = 0xC2B2AE3D27D4EB4F
RSALT = 0x5BF036354E1B2C3D
MUL1 = 0xBF58476D1CE4E5B9
MUL2 = 0x94D049BB133111EB


# ------------------------------------------------------------------ mixing
def mix64(z):
    z &= MASK64
    z = ((z ^ (z >> 30)) * MUL1) & MASK64
    z = ((z ^ (z >> 27)) * MUL2) & MASK64
    return (z ^ (z >> 31)) & MASK64


def _unxor(y, shift):
    x = y
    s = shift
    while s < 64:
        x = y ^ (x >> shift)
        s += shift
    return x & MASK64


def inv_mix64(z):
    """exact inverse permutation of mix64 (used only to build a DP-conditional
    census sample; the walk itself never calls it)."""
    z &= MASK64
    z = _unxor(z, 31)
    z = (z * pow(MUL2, -1, 1 << 64)) & MASK64
    z = _unxor(z, 27)
    z = (z * pow(MUL1, -1, 1 << 64)) & MASK64
    return _unxor(z, 30)


# ------------------------------------------------------ key field geometry
def half_of(key_bits):
    return key_bits >> 1


def dp_field_mask(key_bits):
    """bits of the key that the DP predicate is allowed to read: the low
    min(64, half) bits of the low half."""
    return (1 << min(64, half_of(key_bits))) - 1


def route_field_mask(key_bits):
    """bits of the key that the router is allowed to read: bits half..half+63."""
    return MASK64 << half_of(key_bits)


def fields_disjoint(key_bits):
    return dp_field_mask(key_bits) & route_field_mask(key_bits) == 0


def dp_word(key, key_bits):
    return mix64((key & dp_field_mask(key_bits)) ^ DPSALT)


def is_dp(key, key_bits, dp_bits):
    return (dp_word(key, key_bits) & ((1 << dp_bits) - 1)) == 0


def route_field(key, key_bits):
    return (key & route_field_mask(key_bits)) >> half_of(key_bits)


def route_of(key, key_bits):
    return (mix64(route_field(key, key_bits) ^ RSALT) >> (64 - ROUTE_BITS)) \
        & ((1 << ROUTE_BITS) - 1)


# ----------------------------------------------- rejected v9b spec (control)
def legacy_dp_word(key):
    return mix64(key ^ DPSALT)


def legacy_is_dp(key, dp_bits):
    return (legacy_dp_word(key) & ((1 << dp_bits) - 1)) == 0


def legacy_route_of(key):
    return (mix64(key ^ RSALT) >> (64 - ROUTE_BITS)) & ((1 << ROUTE_BITS) - 1)


# ------------------------------------------------------------ fixed cap rule
def cap_over_mean(reports_prod, nodes_prod=NODES_PROD):
    """declared production cap / declared production mean per node (exact)."""
    return Fraction(CAP_RECORDS_PER_NODE) / Fraction(reports_prod, nodes_prod)


def node_cap(expected_reports, nodes, ratio):
    """FIXED per-node cap, computed before the run from the DECLARED expected
    report count (never from the observed count)."""
    v = Fraction(ratio) * Fraction(expected_reports, nodes)
    return max(1, -(-v.numerator // v.denominator))


# ------------------------------------------------------- exact ln / sqrt / S_ops
def ln_bounds(num, den, terms=200):
    """exact rational bracket [lo, hi] for ln(num/den), lo <= ln <= hi."""
    z = Fraction(num - den, num + den)
    z2 = z * z
    s, t = Fraction(0), z
    for k in range(terms):
        s += t / (2 * k + 1)
        t *= z2
    lo = 2 * s
    tail = 2 * t / ((2 * terms + 1) * (1 - z2))
    return lo, lo + tail


def ceil_sqrt(fr):
    s = 1
    while Fraction(s * s) < fr:
        s *= 2
    lo, hi = s // 2, s
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if Fraction(mid * mid) >= fr:
            hi = mid
        else:
            lo = mid
    return hi


def S_ops(n_prime):
    """S_ops = ceil(KAPPA*sqrt(2*n_prime*ln(1/0.60))), exact rational bracket.
    This is the ONLY definition of the walk budget: the toy calls it with the
    toy N', the certificate calls it with production N'."""
    ln_lo, ln_hi = ln_bounds(100, 60)
    return ceil_sqrt(KAPPA * KAPPA * Fraction(2) * Fraction(n_prime) * ln_hi)
