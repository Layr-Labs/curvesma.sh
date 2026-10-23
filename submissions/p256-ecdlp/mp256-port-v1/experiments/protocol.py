#!/usr/bin/env python3
"""p256-ecdlp mp256-port-v1 -- SHARED PRODUCTION PROTOCOL.

Independent reimplementation of the negation-class distinguished-point rho
construction published in 5fa78d62 (Meganpark980320) with the compiled
batch-round step price of 2ec1a0b0 (jungjipdo).  One place for every constant
that the toy and the ledger both read:

  KAPPA = 21/20 (1.05), FAIL = 0.605      walk budget / success model
  step price 1.35 units, doubling 1.80, traversal 1.45, inversion 2^17 per
  round of at least 2^15 walkers, C = 1000 instructions per unit
  DP predicate width 64 on the low key field, 16-bit route on a DISJOINT field
  2^16 storage nodes, per-node record < 2^56.4 bytes, R = 2048, L_max = 20*2^64
  W_max = 2^128.06, success floor 0.3925

Key geometry (production): key K = x of the class {P,-P}; the DP predicate
reads only key bits 0..63 through mix64, the router reads only bits 128..191,
so the two fields are disjoint and mix64 is a bijection on 64-bit words.
"""
from fractions import Fraction
from math import isqrt

C = 1000                             # instructions per reference unit (borrowed)
KAPPA = Fraction(21, 20)             # 1.05 walk-quality constant (one place)
FAIL_NUM, FAIL_DEN = 605, 1000       # modelled per-walk failure probability 0.605
DP_BITS_PROD = 64
ROUTE_BITS = 16
NODES_PROD = 1 << 16
RECORD_BITS = 26 * 8                 # route16 + dp-fp64 + trail-id64 + step64
CAP_BYTES_PER_NODE = 1 << 56.4 if False else 102211695333519936   # < 2^56.4 B
W_MAX_EXP = Fraction(12806, 100)     # W_max = 2^(12806/100)
SUCCESS_FLOOR = Fraction(3925, 10000)
R_PROD = 2048
L_MAX = 20 * (1 << 64)
MASK64 = (1 << 64) - 1
DPSALT = 0xC2B2AE3D27D4EB4F
RSALT = 0x5BF036354E1B2C3D
MUL1 = 0xBF58476D1CE4E5B9
MUL2 = 0x94D049BB133111EB

# prices, in units (1 unit = C = 1000 AArch64 instructions)
STEP_PRICE = Fraction(117, 100)      # 1.17   <= 1,050 counted + 11.4% margin
DBL_PRICE = Fraction(180, 100)       # 1.80   <= 1,635 counted + 10.1% margin
TRAV_PRICE = Fraction(145, 100)      # 1.45   <= 1,308 counted + 10.9% margin
INV_PRICE = 1 << 17                  # 131,072 instructions, shared by >= 2^15 walkers


def mix64(z):
    z &= MASK64
    z = ((z ^ (z >> 30)) * MUL1) & MASK64
    z = ((z ^ (z >> 27)) * MUL2) & MASK64
    return (z ^ (z >> 31)) & MASK64


def half_of(key_bits):
    return key_bits >> 1


def dp_field_mask(key_bits):
    return (1 << min(64, half_of(key_bits))) - 1


def route_field_mask(key_bits):
    return MASK64 << half_of(key_bits)


def fields_disjoint(key_bits):
    return dp_field_mask(key_bits) & route_field_mask(key_bits) == 0


def dp_word(key, key_bits):
    return mix64((key & dp_field_mask(key_bits)) ^ DPSALT)


def is_dp(key, key_bits, dp_bits):
    return (dp_word(key, key_bits) & ((1 << dp_bits) - 1)) == 0


def class_key(x, y_canon, key_bits=128):
    """Toy class key: low 64-bit field mix64(x) (DP predicate), high field
    mix64(y_canon) (router) -- same geometry as the 256-bit production key."""
    return ((mix64(y_canon) << 64) | mix64(x)) & ((1 << key_bits) - 1)


def route_of(key, key_bits=128, route_bits=ROUTE_BITS):
    f = (key & route_field_mask(key_bits)) >> half_of(key_bits)
    return (mix64(f ^ RSALT) >> (64 - route_bits)) & ((1 << route_bits) - 1)


# --------------------------------------------------- exact decimal-free bounds
def ln_bounds(num, den, terms=48):
    """exact rational bracket [lo,hi] for ln(num/den); atanh series, z small."""
    z = Fraction(num - den, num + den)
    z2 = z * z
    s, t = Fraction(0), z
    for k in range(terms):
        s += t / (2 * k + 1)
        t *= z2
    lo = 2 * s
    tail = 2 * t / ((2 * terms + 1) * (1 - z2))
    return lo, lo + tail


def ceil_sqrt_frac(fr):
    """ceil(sqrt(fr)) exactly, for a positive Fraction, without denominator blowup."""
    num, den = fr.numerator, fr.denominator
    v = isqrt(num // den)                      # sqrt of the integer part
    while v * v * den < num:
        v += 1
    while v > 0 and (v - 1) * (v - 1) * den >= num:
        v -= 1
    return v


def N_prime(r, R=R_PROD):
    """N' = N*R/(R-1) with N = (r-1)/2 negation classes."""
    return Fraction(((r - 1) // 2) * R, R - 1)


def S_ops(n_prime, kappa=KAPPA):
    """S_ops = ceil(kappa*sqrt(2*n_prime*ln(1/0.605))) -- the ONLY walk budget."""
    ln_lo, ln_hi = ln_bounds(FAIL_DEN, FAIL_NUM)
    return ceil_sqrt_frac(kappa * kappa * 2 * Fraction(n_prime) * ln_hi)


def node_cap(expected_reports, nodes, ratio):
    """FIXED per-node cap from the DECLARED expected report count."""
    v = Fraction(ratio) * Fraction(expected_reports, nodes)
    return max(1, -(-v.numerator // v.denominator))
