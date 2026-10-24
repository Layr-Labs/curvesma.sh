#!/usr/bin/env python3
"""Single source of every production constant used by this package.

No module may restate a numerical constant of the walk: the step mixture, the
walk-quality allowance, the success target, the prices and the claim all live
here and are imported.  verify_certificate.py asserts this import discipline.

v10f is the WINDOWLESS build.  There is no look-back window, no checkpoint field
and no escape of any kind: the step is a stationary function of the class key
alone.  That is the merge-safe shape (equal class -> equal next class, for every
trail, forever) and it removes the b68a086c defect class outright.
"""
from fractions import Fraction

# --- curve constants (BLS12-381, tracks/bls12-381-fixed-g1-invert-g2/target-profile.json)
P_FIELD = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORDER = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
N_CLASSES = (R_ORDER - 1)//12          # Frobenius classes of GT, ord_r(p) = 12

# --- organiser prices (tracks/docs/COST_MODEL.md, row "BLS12-381 fixed-G1 pairing inversion into G2")
C_REF   = 20000     # reference dense Fp12 multiply, "1 target operation"
FP_MUL  = 383       # blst base-field multiply, measured 383 instructions
FP_ADD  = 20        # 6-limb add/sub with conditional reduction (6 ADC + 6 SBC + mask)
CMP6    = 24        # 6-limb compare with conditional move
HASH_W  = 24        # one hashing word folded from two key limbs

# --- walk constants
FRAC_SQR = Fraction(3, 8)     # fraction of squaring steps: routing step word below 3
MU_WALK  = Fraction(11, 10)   # walk-quality allowance (random-map inefficiency)
TARGET   = Fraction(392, 1000)  # step-budget target; the model success at S is 0.4711
SUCCESS  = Fraction(39, 100)    # required success probability (COST_MODEL.md line 45)

# --- trail cutoff and report cap (both are constructions, not tail claims)
TRAIL_METERS = 16             # trail aborts after 16 expected DP distances = 16/theta steps
CAP_SLACK    = 1 << 25        # report cap = (11/10) S theta + 2^25   (Chernoff)
THETA_BITS   = 80             # production DP density 2^-80

# --- no escape / no history (v10f): declared, so the certificate can assert it
WINDOW       = 0              # class keys remembered: none
PROBES       = 0              # deterministic probes: none

# --- loss rows: every component of the post-event budget, each derived in argument.md
LOSS_ROWS = (
    ("cutoff_waste",    Fraction(5, 10**7)),    # 2 exp(-TRAIL_METERS), twice, on S
    ("cap_overflow",    Fraction(1, 10**20)),   # Chernoff at 10 per cent over the mean
    ("fp_false_confirm", Fraction(1, 10**25)),  # 2^125.8 compares x 2^-62, each rejected
    ("degenerate_solve", Fraction(1, 10**70)),  # den = 0 mod r, <= 12/r
    ("twist_bias",      Fraction(1, 10**6)),    # 32-bit twist field mod 12, 12 x 2^-32
    ("model_reserve",   Fraction(297, 10**6)),  # undisclosed random-map deviation
)
LOSSES = sum(x for _, x in LOSS_ROWS)           # = 2.985e-4 by construction

# --- claim
# W_raw = 2^125.71174, W_stress = 2^125.75316 (windowless ledger: nine fold rows
# summing to 1813, C = 20000, FP_MUL = 383, cyclo_sqr = 11632).  Claim = 125.79 is
# 0.03684 above W_stress and 0.07826 above W_raw, 0.020 under the 125.83 tie cap
# and 0.030 under the live lowest claim 125.84 (kaankolcu a9c23016).
CLAIM = Fraction(12579, 100)
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
