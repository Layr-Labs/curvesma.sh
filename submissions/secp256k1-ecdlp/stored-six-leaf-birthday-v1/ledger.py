"""Exact ledger for the stored six-leaf birthday submission.

Computes, in exact integers, the total instruction count T of the attack, the
success bound, and the 15-decimal claim: the least multiple of 10^-15 that is at
least log2(T/1000). Prints every number quoted in proof.md.

No randomness is used; re-running prints identical values.
"""
import math
from decimal import Decimal, getcontext, ROUND_CEILING
from fractions import Fraction

getcontext().prec = 90

# ---- secp256k1 parameters (tracks/secp256k1-ecdlp/target-profile.json) ----
P = 2 ** 256 - 2 ** 32 - 977
R = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

# ---- cost model (tracks/docs/COST_MODEL.md) ------------------------------
C = 1000       # scalar AArch64 instructions per reference group operation
MUL = 179      # measured 5x52 base-field multiplication
SQR = 179      # squaring charged as a multiplication (S = M convention)
ADD = 30       # modular add/sub with conditional reduction, conservative
INS = 15       # one table insert: index, load, 4-word compare, store
REC_FIX = 24   # per-record pointer/loop work
SV = 64        # one scalar coin draw (RNDR + range test), charged per scalar

# ---- success target ------------------------------------------------------
# Poisson: P(no colliding pair) <= exp(-mu); require 1 - exp(-mu) >= 0.39.
TARGET = Fraction(494296, 10 ** 6)          # -ln(0.61) = 0.494296


def mu_of(s: int) -> Fraction:
    """Exact mean number of unordered leaf pairs with equal x: C(s,2) * 2/r."""
    return Fraction(s * (s - 1), 2) * Fraction(2, R)


s_target = math.isqrt(int(TARGET * R))       # smallest s with s(s-1) >= TARGET*r
while Fraction(s_target * (s_target - 1)) < TARGET * R:
    s_target += 1
# rectangular layout: n X-points and n Y-points give n^2 records and 6 n^2 leaves
n_axis = math.isqrt(s_target // 6 + 1)
while 6 * n_axis * n_axis < s_target:
    n_axis += 1
records = n_axis * n_axis
s = 6 * records
assert mu_of(s) >= TARGET
GEN_X = 645410                               # [a]P + [b]Q, two scalar mults
GEN_Y = 642930                               # [c]P
VERIFY = 2 ** 21                             # resolve + [x]P = Q + output
ZERO = 2 ** 133 // 64                        # bulk-zero the 2^128-slot table

# ---- per-record arithmetic ----------------------------------------------
# d_t = x_Y - beta^t x_X                 3 A
# D   = d_0 d_1 d_2                      2 M
# batched inverse of D (amortised)       3 M
# w_t = D^-1 d_j d_k                     6 M
# u_+/-^2                                2 S
# w_t^2                                  3 S
# x_{t,+/-} = u^2 w^2 - x_Y - x_t        6 M + 12 A
REC_M, REC_S, REC_A = 2 + 3 + 6 + 6, 2 + 3, 3 + 2 + 12
REC_ARITH = REC_M * MUL + REC_S * SQR + REC_A * ADD + REC_FIX
REC = REC_ARITH + 6 * INS

T = records * REC + (2 * n_axis) * (GEN_X + GEN_Y) + (2 * n_axis) * 2 * SV + VERIFY + ZERO
W = Fraction(T, C)

# ---- claim: least 15-decimal c with 2^c >= W -----------------------------


def pow2(c: Fraction, bits: int = 240):
    """Interval [lo, hi] for 2^c from a rational c, via 2^floor * 2^frac."""
    whole = c.numerator // c.denominator
    frac = c - whole
    # 2^frac by Taylor series on exp(frac*ln2), with ln2 bracketed
    ln2_lo = Decimal(2).ln()
    x = Decimal(frac.numerator) / Decimal(frac.denominator) * ln2_lo
    term = Decimal(1)
    acc = Decimal(1)
    for k in range(1, bits):
        term = term * x / k
        acc += term
    scale = Decimal(2) ** whole
    v = acc * scale
    return v


def claim_of(w: Fraction) -> Fraction:
    start = (Decimal(w.numerator).ln() - Decimal(w.denominator).ln()) / Decimal(2).ln() * 10 ** 15
    c = Fraction(int(start.to_integral_value(rounding=ROUND_CEILING)), 10 ** 15)
    for _ in range(200):
        if pow2(c) >= Decimal(w.numerator) / Decimal(w.denominator) and \
           pow2(c - Fraction(1, 10 ** 15)) < Decimal(w.numerator) / Decimal(w.denominator):
            return c
        c += Fraction(1, 10 ** 15) if pow2(c) < Decimal(w.numerator) / Decimal(w.denominator) else -Fraction(1, 10 ** 15)
    raise AssertionError("claim not bracketed")


CLAIM = claim_of(W)
log2W = (Decimal(W.numerator).ln() - Decimal(W.denominator).ln()) / Decimal(2).ln()

if __name__ == "__main__":
    print("leaves S                %d  (2^%.4f)" % (s, math.log2(s)))
    print("records                 %d  (2^%.4f)" % (records, math.log2(records)))
    print("mean colliding pairs    %s = %.6f" % (mu_of(s), float(mu_of(s))))
    print("success 1-exp(-mu)      %.6f" % (1 - math.exp(-float(mu_of(s)))))
    print("per-record arithmetic   %d M %d S %d A + %d = %d" % (REC_M, REC_S, REC_A, REC_FIX, REC_ARITH))
    print("per-record total        %d  (6 inserts)" % REC)
    print("per-leaf instructions   %.3f" % (REC / 6))
    print("generators              %d" % ((2 * n_axis) * (GEN_X + GEN_Y)))
    print("table zeroing           %d" % ZERO)
    print("verification/output     %d" % VERIFY)
    print("total T                 %d" % T)
    print("log2 T                  %.9f" % math.log2(T))
    print("log2(T/1000)            %.12f" % float(log2W))
    print("claim                   %s / 10^15 = %.12f" % (CLAIM.numerator, float(CLAIM)))
    print("2^claim >= W            %s" % (pow2(CLAIM) >= Decimal(T) / C))
    print("2^(claim-1e-15) < W     %s" % (pow2(CLAIM - Fraction(1, 10 ** 15)) < Decimal(T) / C))

