"""Exact ledger for stored-six-leaf-birthday-v2.

Fixes the four rejected points of v1:
  1. the claim is emitted as an exact reduced rational equal to the least
     15-decimal bound on log2(T/1000);
  2. success is bounded with the collision-pairing (each cross-index collision
     forces a partner), so the requirement is P(no collision) <= exp(-mu/2);
  3. the table never overwrites: 2^129 bucket heads chain every record;
  4. storage, records, axis scalars and buffers are itemised and charged.

Run: python3 ledger.py
"""
import math
from decimal import Decimal, getcontext
from fractions import Fraction

getcontext().prec = 90

P = 2 ** 256 - 2 ** 32 - 977
R = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

C = 1000
MUL = 179
SQR = 179
ADD = 30
INS = 24          # one chained insert: index, head load, chain compare, record store, link
REC_FIX = 24
GEN_X = 645410
GEN_Y = 642930
COIN = 64
VERIFY = 2 ** 21
HEADS_BITS = 129  # bucket count 2^129

# success: mu/2 >= ln(1/0.61) with margin -> mu >= 1.0
MU_MIN = Fraction(1, 1)


def mu_of(s):
    return Fraction(s * (s - 1), R)   # C(s,2) * 2/r


S0 = math.isqrt(int(MU_MIN * R)) + 1          # S(S-1) >= mu*R  =>  S ~ sqrt(R)
n = math.isqrt(S0 // 6) + 1              # S = 6 n^2
while 6 * n * n * (6 * n * n - 1) < MU_MIN * R:
    n += 1
while n > 1 and 6 * (n - 1) ** 2 * (6 * (n - 1) ** 2 - 1) >= MU_MIN * R:
    n -= 1
records = n * n
S = 6 * records
MU = mu_of(S)
SUCCESS = 1 - math.exp(-float(MU) / 2)

REC_M, REC_S, REC_A = 17, 5, 17
REC_ARITH = REC_M * MUL + REC_S * SQR + REC_A * ADD + REC_FIX
REC = REC_ARITH + 6 * INS

ZERO = 2 ** HEADS_BITS * 8 // 64
axis_stores = (2 * n) * 2 + (3 * n)          # X,Y points (2 STP each) and 3n scalars
T = records * REC + 2 * n * (GEN_X + GEN_Y) + 3 * n * COIN + axis_stores + ZERO + VERIFY
W = Fraction(T, C)


def pow2(c: Fraction):
    whole = c.numerator // c.denominator
    frac = c - whole
    x = Decimal(frac.numerator) / Decimal(frac.denominator) * Decimal(2).ln()
    term, acc = Decimal(1), Decimal(1)
    for k in range(1, 220):
        term = term * x / k
        acc += term
    return acc * (Decimal(2) ** whole)


def claim_of(w: Fraction) -> Fraction:
    lo, hi = 0, 1 << 10
    # bracket by bisection on the 10^-15 grid
    est = (Decimal(w.numerator).ln() - Decimal(w.denominator).ln()) / Decimal(2).ln()
    c = int((est * 10 ** 15).to_integral_value(rounding="ROUND_CEILING"))
    cand = Fraction(c, 10 ** 15)
    for _ in range(200):
        if pow2(cand) >= Decimal(w.numerator) / Decimal(w.denominator) and \
           pow2(cand - Fraction(1, 10 ** 15)) < Decimal(w.numerator) / Decimal(w.denominator):
            return cand
        cand += Fraction(1, 10 ** 15) if pow2(cand) < Decimal(w.numerator) / Decimal(w.denominator) else -Fraction(1, 10 ** 15)
    raise AssertionError


CLAIM = claim_of(W)

if __name__ == "__main__":
    print("n                    %d  (2^%.4f)" % (n, math.log2(n)))
    print("records n^2          %d  (2^%.4f)" % (records, math.log2(records)))
    print("leaves S=6n^2        %d  (2^%.4f)" % (S, math.log2(S)))
    print("mu = S(S-1)/r        %.9f" % float(MU))
    print("success lb 1-e^-mu/2 %.9f" % SUCCESS)
    print("per record           %d = %d arith + %d inserts" % (REC, REC_ARITH, 6 * INS))
    print("records row          %d" % (records * REC))
    print("axis stores          %d" % axis_stores)
    print("zeroing heads        %d" % ZERO)
    print("total T              %d" % T)
    print("log2 T               %.12f" % math.log2(T))
    print("log2(T/1000)         %.15f" % (math.log2(T) - math.log2(1000)))
    print("claim                %d/%d" % (CLAIM.numerator, CLAIM.denominator))
    print("claim as float       %.15f" % float(CLAIM))
    print("2^claim >= W         %s" % (pow2(CLAIM) >= Decimal(T) / C))
    print("2^(claim-1e-15) < W  %s" % (pow2(CLAIM - Fraction(1, 10 ** 15)) < Decimal(T) / C))
    print("memory bytes         %d  (2^%.2f)" % (2 ** HEADS_BITS * 8 + S * 64 + 3 * n * 32 + 2 * n * 64,
                                                  math.log2(2 ** HEADS_BITS * 8 + S * 64 + 3 * n * 32 + 2 * n * 64)))
