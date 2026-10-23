#!/usr/bin/env python3
"""Exact ledger for gt-squaring-branch-40pct-v14.  Integers/Fractions only.

Unit rule (docs/COST_MODEL.md): cost = H + R/C, C = 20,000 instructions per dense
Fp12 multiplication (the reference operation), H = number of reference operations.
Non-reference instruction prices, measured blst convention: Fp multiplication 383,
modular add/sub 45 (both as called routines).
"""
import json
from fractions import Fraction as F
from math import isqrt

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
RS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 20000
RTAB = 1024

# instruction rows.  The Fp-level counts of the two arithmetic blocks are the
# numbers counted by running scripts/gt_prod.py (evidence/prod_run.json).
prod = json.load(open(__file__.rsplit("/", 2)[0] + "/evidence/prod_run.json"))
SQR_ARITH = prod["sqr_rows"]["instr"]          # cyclotomic squaring, counted
KEY_ARITH = prod["key_rows"]["instr"]          # 6-key +-pair key computation, counted

COMMON = (("trace_re_plus_branch_bits", 60), ("exponent_pair_update", 100),
          ("cycle_compares_loop_margin", 140))
MUL_EXTRA = (("key_six_plusminus", KEY_ARITH), ("min_six_argmax_index", 180),
             ("table_address_and_load", 150), ("exponent_pair_add", 60))
SQR_EXTRA = (("cyclotomic_squaring", SQR_ARITH), ("exponent_pair_double", 40))
I_COMMON = sum(v for _, v in COMMON)
I_MUL = I_COMMON + sum(v for _, v in MUL_EXTRA)
I_SQR = I_COMMON + sum(v for _, v in SQR_EXTRA)

S_NUM, S_DEN = 3, 8                              # branch probability, K mod 8 < 3
C_MUL = 1 + F(I_MUL, C)                          # one dense mul + overhead
C_SQR = F(I_SQR, C)                              # no dense mul
C_AVG = F(S_DEN - S_NUM, S_DEN) * C_MUL + F(S_NUM, S_DEN) * C_SQR

# collision law: two classes merge only if their branch labels differ (both
# squaring cannot merge, squaring is injective) and, for two multiply steps,
# only if the multiplier index differs.  Merge probability per pair
#   rho/N with rho = (1 - s^2)(1 - 1/R)
MIX_INFLATION = 1 / (1 - F(S_NUM, S_DEN) ** 2)   # 64/55
TAB_INFLATION = F(RTAB, RTAB - 1)

N = (RS - 1) // 12
assert (RS - 1) == 12 * N
NP = F(N) * MIX_INFLATION * TAB_INFLATION

P_SUCCESS = F(39, 100)
Q = 1 - P_SUCCESS                                 # 61/100


def ln_enclosure(q=Q, terms=40):
    """[lo,hi] enclosure of ln(1/q) from the atanh series, ALL terms added (no
    alternating-sign bug: z>0, every term positive)."""
    z = (1 - q) / (1 + q)
    z2 = z * z
    lo = hi = F(0)
    pw = z
    for k in range(1, 2 * terms, 2):
        t = pw / k
        lo += t
        hi += t
        pw *= z2
    tail = pw / (2 * terms + 1)
    return 2 * (lo), 2 * (hi + tail)


LN_LO, LN_HI = ln_enclosure()
assert LN_LO < F(49430, 100000) and LN_HI > F(49429, 100000)   # ln(100/61) = 0.4942963...


def sqrt_ceil(x):
    x = F(x)
    n_, d_ = x.numerator, x.denominator
    s = isqrt(n_ // d_)
    while s * s * d_ < n_:
        s += 1
    while s > 0 and (s - 1) * (s - 1) * d_ >= n_:
        s -= 1
    return s


S_STEPS = sqrt_ceil(2 * NP * LN_HI)               # steps for success 0.39
E_STEPS = 1 << 85                                 # detection window, steps
SETUP = 1 << 65                                   # tables, start points, pairings, output

W = F(S_STEPS) * C_AVG + F(E_STEPS) * C_AVG + SETUP


def cmp_pow2(x, k, scale=1000):
    x = F(x)
    return x.numerator ** scale < (1 << k) * x.denominator ** scale


def log2_interval(x, scale=1000):
    lo = (x.numerator.bit_length() - x.denominator.bit_length() - 1) * scale
    hi = lo + scale
    while not cmp_pow2(x, hi, scale):
        hi += scale
    while cmp_pow2(x, lo, scale):
        lo -= scale
    while hi - lo > 1:
        m = (lo + hi) // 2
        if cmp_pow2(x, m, scale):
            hi = m
        else:
            lo = m
    return lo, hi


def fmt(x):
    lo, hi = log2_interval(x)
    return "%d.%03d (%d.%03d)" % (lo // 1000, lo % 1000, hi // 1000, hi % 1000)


if __name__ == "__main__":
    print("rows: common %d  mul %d (%.5f u)  sqr %d (%.5f u)  sqr_arith %d  key %d"
          % (I_COMMON, I_MUL, C_MUL, I_SQR, C_SQR, SQR_ARITH, KEY_ARITH))
    print("c_mul %.5f c_sqr %.5f c_avg %.5f mix_inflation %.5f" % (C_MUL, C_SQR, C_AVG, MIX_INFLATION))
    print("N 2^%s ; N' 2^%s" % (fmt(F(N)), fmt(NP)))
    print("S(P=.39) = 2^%s steps" % fmt(F(S_STEPS)))
    print("W = 2^%s" % fmt(W))
    lo, hi = log2_interval(W)
    print("W log2 in (%d.%03d, %d.%03d); claim 125.73; gate 125.74" % (lo // 1000, lo % 1000, hi // 1000, hi % 1000))
    json.dump({"I_COMMON": I_COMMON, "I_MUL": I_MUL, "I_SQR": I_SQR,
               "SQR_ARITH": SQR_ARITH, "KEY_ARITH": KEY_ARITH,
               "c_mul": str(C_MUL), "c_sqr": str(C_SQR), "c_avg": str(C_AVG),
               "mix_inflation": str(MIX_INFLATION), "N": N, "Np_num": NP.numerator, "Np_den": NP.denominator,
               "S_steps": S_STEPS, "E_steps": E_STEPS, "setup": SETUP,
               "W_num": W.numerator, "W_den": W.denominator,
               "W_log2_lo": lo, "W_log2_hi": hi,
               "ln_lo": str(LN_LO), "ln_hi": str(LN_HI)},
              open(__file__.rsplit("/", 2)[0] + "/evidence/ledger.json", "w"), indent=1)
