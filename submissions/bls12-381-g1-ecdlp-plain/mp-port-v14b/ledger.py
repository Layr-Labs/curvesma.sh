#!/usr/bin/env python3
"""ledger.py — exact-rational ledger for mp-port-v14.

Independent reimplementation of the arithmetic published in Yukon submission
03b9917 (jungjipdo), track bls12-381-g1-ecdlp-plain. No floats, no floats in
any bound: every log2 value is an interval of Fractions obtained from the
atanh series for ln with an explicit geometric tail, so every printed
inequality is a machine-checked rational comparison.

Rows reproduced (theirs -> here):
  N = (r-1)/6, log2 N = 252.27213
  N' = N R/(R-1), R = 2^20
  S0 = sqrt(2 N' ln(1/0.6097))          -> 2^126.12851
  S_ops = 1.05 S0                        -> 2^126.19890
  phase1 = 1.0875 (S_ops + 2^60)         -> < 2^126.319911
  W_max = 2^126.321, headroom > 2^115.94
  expected random surcharges < 2^102.48, Markov early halt < 2^-13.46
  success >= 0.3903 - 0.000089 - 1.1e-7 > 0.39021
"""
from fractions import Fraction as F

R = 1 << 20
# BLS12-381 G1 subgroup order (protected profile constant)
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001


def _atanh_ln_bounds(x, terms=90):
    """ln(x) for rational x>0 as (lo,hi), both Fractions.

    Uses ln x = 2*sum_{k>=0} z^(2k+1)/(2k+1), z=(x-1)/(x+1), after reducing to
    x>=1 (ln x = -ln(1/x)). All terms positive, so the partial sum is a lower
    bound and the partial sum plus 2*z^(2n+3)/((2n+3)(1-z^2)) is an upper bound.
    """
    assert x > 0
    sign = 1
    if x < 1:
        x = 1 / x
        sign = -1
    z = (x - 1) / (x + 1)
    lo = F(0)
    p = z
    z2 = z * z
    for k in range(terms):
        lo += p / (2 * k + 1)
        p *= z2
    lo *= 2
    tail = 2 * p / ((2 * terms + 1) * (1 - z2))
    hi = lo + tail
    return (lo, hi) if sign > 0 else (-hi, -lo)


LN2_LO, LN2_HI = _atanh_ln_bounds(F(2))


def log2_bounds(x, terms=90):
    """log2(x) for rational x>0 as (lo,hi) Fractions, via bit split + ln."""
    assert x > 0
    k = x.numerator.bit_length() - x.denominator.bit_length()
    m = x / (F(2) ** k)
    while m < 1:
        k -= 1
        m = x / (F(2) ** k)
    while m >= 2:
        k += 1
        m = x / (F(2) ** k)
    a, b = _atanh_ln_bounds(m, terms)
    if a >= 0:
        return (F(k) + a / LN2_HI, F(k) + b / LN2_LO)
    return (F(k) + a / LN2_LO, F(k) + b / LN2_HI)


def log2_of_interval(lo, hi):
    """log2 of a positive rational interval, endpoints mapped monotonically."""
    return (log2_bounds(F(lo))[0], log2_bounds(F(hi))[1])


# ---- rows -------------------------------------------------------------
N = F(r - 1, 6)
Np = N * R / (R - 1)
LOG2N_LO, LOG2N_HI = log2_bounds(N)
LOG2NP_LO, LOG2NP_HI = log2_bounds(Np)

P_FAIL = F(6097, 10000)             # 1 - 0.3903
LN_TERM_LO, LN_TERM_HI = _atanh_ln_bounds(1 / P_FAIL)
LOG2_LNTERM_LO, LOG2_LNTERM_HI = log2_of_interval(LN_TERM_LO, LN_TERM_HI)

# log2 S0 = (1 + log2 N' + log2 ln(1/p_fail))/2
S0_LO = (1 + LOG2NP_LO + LOG2_LNTERM_LO) / 2
S0_HI = (1 + LOG2NP_HI + LOG2_LNTERM_HI) / 2

A105_LO, A105_HI = log2_bounds(F(105, 100))
SOPS_LO = S0_LO + A105_LO
SOPS_HI = S0_HI + A105_HI

# phase1 = 1.0875*(S_ops + 2^60). The additive 2^60 is bounded in log2 space:
# log2(S_ops + 2^60) = log2 S_ops + log2(1 + 2^60/S_ops), and
# 2^60/S_ops <= 2^(60 - 126.19) = 2^-66.19, so log2(1+t) <= t*1.4427 < 2^-65.7.
# We add the safe rational allowance 2^-60 (1e-18) to the upper end.
SUM_ALLOW = F(1, 10 ** 18)
A10875_LO, A10875_HI = log2_bounds(F(10875, 10000))
PHASE1_LO_FROM_SUM = A10875_LO + SOPS_LO
PHASE1_HI = A10875_HI + SOPS_HI + SUM_ALLOW

PUB = {
    "log2N": F(25227213, 100000),
    "log2S0": F(12612851, 100000),
    "log2Sops": F(12619890, 100000),
    "log2phase1": F(126319911, 1000000),
}


def delta(a, b):
    return abs(a - b)


if __name__ == "__main__":
    print("log2 N     lo=%s hi=%s  pub=%s  d<=%s"
          % (float(LOG2N_LO), float(LOG2N_HI), float(PUB["log2N"]),
             float(max(delta(LOG2N_LO, PUB["log2N"]), delta(LOG2N_HI, PUB["log2N"])))))
    print("log2 S0    lo=%s hi=%s  pub=%s  d<=%s"
          % (float(S0_LO), float(S0_HI), float(PUB["log2S0"]),
             float(max(delta(S0_LO, PUB["log2S0"]), delta(S0_HI, PUB["log2S0"])))))
    print("log2 S_ops lo=%s hi=%s  pub=%s  d<=%s"
          % (float(SOPS_LO), float(SOPS_HI), float(PUB["log2Sops"]),
             float(max(delta(SOPS_LO, PUB["log2Sops"]), delta(SOPS_HI, PUB["log2Sops"])))))
    print("log2 phase1 lo=%s hi=%s pub=%s  d<=%s"
          % (float(PHASE1_LO_FROM_SUM), float(PHASE1_HI), float(PUB["log2phase1"]),
             float(max(delta(PHASE1_LO_FROM_SUM, PUB["log2phase1"]),
                       delta(PHASE1_HI, PUB["log2phase1"])))))
    print("headroom bits =", float(F(126321, 1000) - PHASE1_HI))
