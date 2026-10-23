#!/usr/bin/env python3
"""Single source of every production constant used by this package.

No module may restate a numerical constant of the walk: the step mixture, the
walk-quality allowance, the success target, the prices and the claim all live
here and are imported.  verify_certificate.py asserts this import discipline.
"""
from fractions import Fraction

# --- curve constants (BLS12-381, tracks/bls12-381-fixed-g1-invert-g2/target-profile.json)
P_FIELD = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORDER = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
N_CLASSES = (R_ORDER - 1)//12          # Frobenius classes of GT, ord_r(p) = 12

# --- organiser prices (tracks/docs/COST_MODEL.md)
C_REF   = 20000     # reference dense Fp12 multiply, "1 target operation"
FP_MUL  = 383       # blst base-field multiply, measured 383 instructions
FP_ADD  = 20        # 6-limb add/sub with conditional reduction (6 ADC + 6 SBC + mask)
CMP6    = 24        # 6-limb compare with conditional move
HASH_W  = 24        # one hashing word folded from two key limbs

# --- walk constants
FRAC_SQR = Fraction(3, 8)     # fraction of squaring steps: low three bits of h below 3
MU_WALK  = Fraction(11, 10)   # walk-quality allowance (random-map inefficiency)
TARGET   = Fraction(392, 1000)  # main-event success budgeted before losses
LOSSES   = Fraction(3, 10000)   # post-event losses (cycles, cut trails, degeneracy)
SUCCESS  = Fraction(39, 100)    # required success probability

# --- claim
CLAIM = Fraction(12589, 100)   # 125.89 < live head 125.90 (strictly below)
STRESS_FOLD = Fraction(13, 10) # hostility factor carried on every judgement-bounded row


def step_budget(mu=MU_WALK, target=TARGET, f=FRAC_SQR, n=N_CLASSES):
    """Least integer S with S^2 (1-f^2) >= 2 mu N ln(1/(1-target)).

    Exact integer search; the logarithm is evaluated once in floating point and
    then corrected by integer comparisons, so the returned S is exact for the
    stated inequality (checked again in verify_certificate.py).
    """
    from math import log, isqrt
    s2 = 2*float(mu)*n*log(1/(1-float(target)))/(1 - float(f)*float(f))
    s = isqrt(int(s2))
    while s*s < s2:
        s += 1
    while s > 1 and (s-1)*(s-1) >= s2:
        s -= 1
    return s


S_STEPS = step_budget()
