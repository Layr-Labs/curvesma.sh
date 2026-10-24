#!/usr/bin/env python3
"""Single source of every production constant used by this package.

No module may restate a numerical constant of the walk: the step mixture, the
walk-quality allowance, the success target, the prices and the claim all live
here and are imported.  verify_certificate.py asserts this import discipline.

v12b is the CANONICALISING build.  Two changes over v10f (27e3ab98), both forced
by its rejection (D1: the p^6 = -1 stationary two-cycle):

  1. The state is held at the canonical Frobenius representative every step, so
     the multiplier applied to the element is always the table entry itself and
     never sigma^(-k) of one.  Holding the representative costs the label
     renormalisation (a,b) <- ((a+delta) p^k, b p^k) mod r: two 255-bit modular
     multiplications per step, a new counted ledger row (RENORM below).
  2. The table is built from INDEPENDENT lower-half exponents 1 <= u <= (r-1)/2,
     so no entry's inverse is an entry.  With (1) the whole (m, m^-1) two-cycle
     family -- the v10f defect -- is unreachable by construction, not merely
     improbable.

There is still no window, no checkpoint field and no escape of any kind: the
step is a stationary function of the class key alone (the merge-safe shape).
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

# --- the canonicalising step: counted renormalisation row
# One label's multiply is evidence/modmul-count.json: clang -O2 -arch arm64 emits
# 143 instructions for mulmod_r (4-limb CIOS Montgomery multiply, branchless final
# subtract), 0/255 differential mismatches against bigint arithmetic across all
# twelve p^k and the edge labels 0,1,2,r-2,r-1.
MODMUL_INSTR = 143            # counted, evidence/modmul-count.json
RENORM = 2*MODMUL_INSTR       # two multiplications per step: label a and label b

# --- trail cutoff and report cap (both are constructions, not tail claims)
TRAIL_METERS = 16             # trail aborts after 16 expected DP distances = 16/theta steps
CAP_SLACK    = 1 << 25        # report cap = (11/10) S theta + 2^25   (Chernoff)
THETA_BITS   = 80             # production DP density 2^-80

# --- no escape / no history (v12b): declared, so the certificate can assert it
WINDOW       = 0              # class keys remembered: none
PROBES       = 0              # deterministic probes: none

# --- loss rows: every component of the post-event budget, each derived in argument.md
LOSS_ROWS = (
    ("cutoff_waste",    Fraction(5, 10**7)),    # 2 exp(-TRAIL_METERS), twice, on S
    ("cap_overflow",    Fraction(1, 10**20)),   # Chernoff at 10 per cent over the mean
    ("fp_false_confirm", Fraction(1, 10**25)),  # 2^125.8 compares x 2^-62, each rejected
    ("degenerate_solve", Fraction(1, 10**70)),  # den = 0 mod r, <= 12/r
    ("column_bias",     Fraction(1, 10**6)),    # 32-bit column field mod 12, 12 x 2^-32
    ("model_reserve",   Fraction(297, 10**6)),  # undisclosed random-map deviation
)
LOSSES = sum(x for _, x in LOSS_ROWS)           # = 2.985e-4 by construction

# --- claim
# W_raw = 2^125.733665, W_stress = 2^125.774477 (canonicalising ledger: the
# windowless nine-row fold 1813 plus the counted RENORM = 286 row on both step
# types, C = 20000, FP_MUL = 383, cyclo_sqr = 11632).  Claim is the least
# 3-decimal value at least W_stress + 0.005: it clears the live 125.84 gate
# (kaankolcu a9c23016, re-polled at submit) and the 125.8335 tie cap.
CLAIM = Fraction(12578, 100)
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
