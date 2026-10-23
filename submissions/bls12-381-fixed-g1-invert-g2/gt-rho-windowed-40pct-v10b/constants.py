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
FRAC_SQR = Fraction(3, 8)     # fraction of squaring steps: routing step word below 3
SQR_THRESHOLD = 3             # st < 3 of the eight routing step values, i.e. f = 3/8
MU_WALK  = Fraction(11, 10)   # walk-quality allowance (random-map inefficiency)
TARGET   = Fraction(392, 1000)  # main-event success budgeted before losses
SUCCESS  = Fraction(39, 100)    # required success probability

# --- trail cutoff and report cap (both are constructions, not tail claims)
TRAIL_METERS = 16             # trail aborts after 16 expected DP distances = 16/theta steps
CAP_SLACK    = 1 << 25        # report cap = (11/10) S theta + 2^25   (Chernoff)
THETA_BITS   = 80             # production DP density 2^-80

# --- escape (ported: look-back window, deterministic probing)
WINDOW       = 24             # class keys remembered
PROBES       = 8              # deterministic probes per blocked step

# --- loss rows: every component of the post-event budget, each derived in argument.md
LOSS_ROWS = (
    ("cutoff_waste",    Fraction(5, 10**7)),    # 2 exp(-TRAIL_METERS), twice, on S
    ("cap_overflow",    Fraction(1, 10**20)),   # Chernoff at 10 per cent over the mean
    ("escape_aborts",   Fraction(1, 10**30)),   # <= 4 tricks/trail cap, Markov on 2^-124
    ("fp_false_confirm", Fraction(1, 10**25)),  # 2^125.8 compares x 2^-62, each rejected
    ("degenerate_solve", Fraction(1, 10**70)),  # den = 0 mod r, <= 12/r
    ("twist_bias",      Fraction(1, 10**6)),    # 32-bit twist field mod 12, 12 x 2^-32
    ("model_reserve",   Fraction(297, 10**6)),  # undisclosed random-map deviation
)
LOSSES = sum(x for _, x in LOSS_ROWS)           # = 3.0e-4 by construction

# --- claim
# W_raw = 2^125.72175, W_stress = 2^125.76579 (v10b ledger: three-multiply selector,
# look-back-24 window).  W_combined = W_stress / (1 - L_hi) = 2^125.79678 is the
# ledger that also carries the measured one-sided loss debit L_hi = 0.02125.
# Claim = 125.81 is 0.01322 above W_combined, 0.04421 above W_stress, 0.08825 above
# W_raw, 0.020 under the 125.83 cap and 0.030 under the live lowest claim 125.84
# (kaankolcu a9c23016, polled 2026-10-23 23:1xZ).
CLAIM = Fraction(12581, 100)
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
