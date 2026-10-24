#!/usr/bin/env python3
"""v13 constants: the QUOTIENT-IN-THE-REPORT-FILTER build.

What changes over v12d (force-squaring HOLD 126.7101) is one thing only, and it
is a structural deletion, not a re-pricing:

  * THE TRANSITION NEVER CANONICALISES.  The state is carried as whatever
    representative the walk reached; the multiplier absorbs the representative
    instead of the state absorbing a canonicalisation:

        m = sigma^{-kappa(y)} M_t          (twelve pre-expanded table columns)
        y <- y * m,   a <- a + u_t p^{-kappa},   b <- b

    For a sigma-partner y' = sigma^j y the key set is the same SET, so t is the
    same, kappa(y') = kappa(y) - j, and m' = sigma^j m; hence y'*m' =
    sigma^j(y*m).  The successor CLASS is a function of the class: merge-safe by
    algebra, on every representative, with one selector per step.

  * THEREFORE no sigma^k is applied to the Fp12 state and no label
    renormalisation exists.  v12b's finding 2 (sigma^k on the state uncharged)
    and the 14 028-instruction row it forced in v12d are both answered by
    DELETION: the sigma powers live in the twelve table columns, built once and
    charged in aux().  The v12b involution trap cannot exist because the state
    is never inverted: sigma^6 = inversion is never applied to a state.

  * THE QUOTIENT LIVES IN THE REPORT FILTER.  The distinguished-point predicate
    and the record key are the sigma-INVARIANT class key min{keys12(y)}; the
    record stores the full label pair (a,b) and the twelve-way solve recovers
    the unknown sigma-shift, including j = 6, which IS the involution (the
    +-pair / label-bit case: x = -(a1+a2)/(b1+b2)).

Nothing is restated outside this module.
"""
from fractions import Fraction

# --- curve constants (tracks/bls12-381-fixed-g1-invert-g2/target-profile.json)
P_FIELD = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORDER = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
N_CLASSES = (R_ORDER - 1)//12          # Frobenius classes of GT, ord_r(p) = 12

# --- organiser prices (tracks/docs/COST_MODEL.md)
C_REF   = 20000     # reference dense Fp12 multiply, "1 target operation"
FP_MUL  = 383       # blst base-field multiply, measured 383 instructions
FP_ADD  = 20        # 6-limb add/sub with conditional reduction
CMP6    = 24        # 6-limb compare with conditional move
HASH_W  = 24        # one hashing word folded from two key limbs

KAPPA_BITS = P_FIELD.bit_length()      # 381; the class key is an Fp element

# --- walk constants
FRAC_SQR = Fraction(3, 8)     # routing step word < 3 -> squaring channel (3 of 8 classes)
MU_WALK  = Fraction(11, 10)   # walk-quality allowance (random-map inefficiency)
TARGET   = Fraction(392, 1000)  # step-budget target (floor 0.39 + 0.002 reserve)
SUCCESS  = Fraction(39, 100)    # required success probability (COST_MODEL.md line 45)

# --- the twelve-key selector: ONE per step.  It produces BOTH the sigma-invariant
#     routing/record key min{keys12(y)} AND the rotation index kappa(y) that
#     selects the multiplier column.  Rows are v12b's counted selector rows
#     (SELECTOR_AB + SELECTOR_ZETA + FOLD_ADDS + FOLD_CMP); v12d charged the SAME
#     row on every step for its canonicalisation, so this is a counted row, not a
#     new estimate.
SELECTOR_AB   = 1*FP_MUL + 1*FP_ADD   # a-b, then g0*(a-b)
SELECTOR_ZETA = 2*FP_MUL              # zeta*a, zeta*d (zeta^2 forms free)
FOLD_ADDS     = 9*FP_ADD              # a-b, two zeta^2 subtractions, six sign flips
FOLD_CMP      = 11*CMP6               # six |.| tests + five-way maximum with cmov
SELECTOR      = SELECTOR_AB + SELECTOR_ZETA + FOLD_ADDS + FOLD_CMP   # = 1613

# --- twelve-column table addressing: index = 12*t + kappa, one multiply-add over
#     a 2^14 x 12 flat array plus the 576-byte fetch.  Charged, not assumed.
COLUMN_ADDR  = 16

# --- the sigma powers moved OFF the per-step path and ONTO the table build.
#     Each of the 12 columns is sigma^{-kappa}(M_t): nine Fp2 constant products
#     = 36 Fp muls (odd-k generic Fp2 product, the worst case) + six conjugations
#     + six accumulations.  Same row as v12d's sigma_apply, applied 12*2^14 times
#     once instead of once per step.
SIGMA_MULS    = 36
SIGMA_CONJ    = 6
SIGMA_BUILD   = SIGMA_MULS*FP_MUL + (SIGMA_CONJ + 6)*FP_ADD   # = 14 028

# --- NO label renormalisation exists in v13: a += u_t p^{-kappa} is one add.
RENORM       = 0

# --- per-step rows: imported by ledger.rows(), restated nowhere
KEY_HASH     = 3*HASH_W       # one 256-bit mix -> table word + DP word + step word
DP_TEST      = 8
HIST_SLOT    = 32             # class-key traffic, 48 bytes
LOOP         = 32             # table addressing, counters
TABLE_FETCH  = 72             # 576-byte table entry: 72 eight-byte loads
LABEL_ADD    = 28             # a += u mod r

# --- trail cutoff and report cap (constructions, not tail claims)
TRAIL_METERS = 16
CAP_SLACK    = 1 << 25
THETA_BITS   = 80
TABLE_BITS   = 14
TABLE_COLS   = 12

# --- declared, so the certificate can assert it: no window, no checkpoint, no
#     counter, no escape.  The step is a stationary map on classes.
WINDOW       = 0
PROBES       = 0
ESCAPE       = 0

# --- loss rows (declared reserves; the work claim uses W_stress only)
LOSS_ROWS = (
    ("cutoff_waste",     Fraction(5, 10**7)),
    ("cap_overflow",     Fraction(1, 10**20)),
    ("fp_false_confirm", Fraction(1, 10**25)),
    ("degenerate_solve", Fraction(1, 10**70)),
    ("dp_detect_tail",   Fraction(74, 10**12)),
    ("model_reserve",    Fraction(297, 10**6)),
)
LOSSES = sum(x for _, x in LOSS_ROWS)

STRESS_FOLD = Fraction(13, 10)

F_ESC = Fraction(0)                 # v13 has no escape branch at all
# Squaring share: DECLARED 3/8 by the routing word, MEASURED 3444/9420 = 0.365605
# over all r-1 fixture states (evidence/toy-census.json t1.sqr_fraction).  The
# smaller value is carried because the squaring channel is the CHEAPER branch, so
# a smaller share is the more expensive price.  Both are reported.
F_SQR_DECLARED = FRAC_SQR
F_SQR_MEASURED = Fraction(3444, 9420)
F_SQR = min(F_SQR_DECLARED, F_SQR_MEASURED)


def step_budget(mu=MU_WALK, target=TARGET, f=F_SQR, n=N_CLASSES):
    """Least integer S with S^2 (1-f^2) >= 2 mu N ln(1/(1-target))."""
    from math import log, isqrt
    s2 = 2*float(mu)*n*log(1/(1-float(target)))/(1 - float(f)*float(f))
    s = isqrt(int(s2))
    while s*s < s2:
        s += 1
    while s > 1 and (s-1)*(s-1) >= s2:
        s -= 1
    return s


S_STEPS = step_budget()
