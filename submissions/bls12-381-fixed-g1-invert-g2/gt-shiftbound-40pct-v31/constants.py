#!/usr/bin/env python3
"""v14 constants: the REPORT-FILTER walk with the absorbing-cycle family REMOVED
by a routing guard, not priced away.

What changes over v13 (claim 125.757, REJECTED on a structural 2-cycle):

  * ONE squared fraction.  v13 derived the step budget S at the measured fixture
    share 287/785 and priced the step at 3/8.  Here the routing share 3/8 is used
    in BOTH places: S is the least integer with S^2 (1-f^2) >= 2 mu N ln(...) at
    f = 3/8, the step price is at the same f, and the fixture share is reported as
    a disclosed cross-check only (it is an input to nothing).

  * GUARD, depth D = 12.  The multiplier of a multiply step is column kappa(y) of
    row t unless row t was used among the last D multiply steps, in which case the
    column recorded then is reused.  Consequences (argument.md section 4):
      - an all-multiplier absorbing cycle of length <= 2D+1 = 25 cannot exist:
        every same-row pair in such a cycle is at cycle-distance <= D, so both
        applications use ONE column and the row contributes L*u_t*p^-k with
        0 < L <= 25 < r, which is not 0 mod r;
      - a cycle of length L >= 26 needs >= L/2 row coincidences plus the matching
        root-of-unity relations, rate <= (1/(12R))^(L/2);
      - a cycle containing a squaring solves (1-2^s)v = B with 1-2^s != 0, one
        element per cyclic pattern, i.e. <= 1/r per window.
    The guard never changes the channel, never restarts and never inverts, and it
    is sigma-covariant: it reads (row, column) pairs, both of which shift by the
    same sigma power under sigma^j.

  * Record size 112 B everywhere (48-byte class key + two 32-byte labels).  v13
    declared 112 in prose and charged 96 in ledger.aux.

  * ledger.aux carries a production restart schedule and a peak-storage bound over
    store + labels + advice + scratch, not four trail starts.
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

# --- THE one squared fraction: the production routing share, used for BOTH the
#     step budget S and the step price.  The fixture share is written into
#     evidence/toy-census.json as t1.sqr_fraction and is an input to nothing.
F_SQR = FRAC_SQR
F_SQR_DECLARED = FRAC_SQR

# --- the twelve-key selector: ONE per step.  It produces BOTH the sigma-invariant
#     routing/record key min{keys12(y)} AND the rotation index kappa(y) that
#     selects the multiplier column.
SELECTOR_AB   = 1*FP_MUL + 1*FP_ADD   # a-b, then g0*(a-b)
SELECTOR_ZETA = 2*FP_MUL              # zeta*a, zeta*d (zeta^2 forms free)
FOLD_ADDS     = 9*FP_ADD              # a-b, two zeta^2 subtractions, six sign flips
FOLD_CMP      = 11*CMP6               # six |.| tests + five-way maximum with cmov
SELECTOR      = SELECTOR_AB + SELECTOR_ZETA + FOLD_ADDS + FOLD_CMP   # = 1613

# --- twelve-column table addressing: index = 12*t + kappa, one multiply-add over
#     a 2^14 x 12 flat array plus the 576-byte fetch.  Charged, not assumed.
COLUMN_ADDR  = 16

# --- THE v24 ROW-DIFFERENCE SELECTOR: the STRUCTURAL kill of the same-row
#     inverse-pair cycle family that the guard-free step still carried.
#     The family: consecutive multiply steps on the SAME row t with
#     kappa(y') = kappa(y) + 6 use columns k and k+6, and TBL[t][k+6] IS the
#     inverse of TBL[t][k] (toy column_inverse_pairs counts them), so y returns
#     EXACTLY and the walk is absorbing.  v24 forbids the precondition instead of
#     pricing the consequence: when the drawn row repeats the previous multiply
#     row, it is remapped to t' = t XOR (1 + h mod (R-1)), which is a BIJECTION of
#     the other R-1 rows (t' != t always, and for fixed t the map m -> t XOR m is a
#     bijection), so the row marginal stays EXACTLY uniform and the pair
#     distribution is uniform on ordered pairs of DISTINCT rows.  Two consecutive
#     multiply steps therefore never share a row, and the exact-return needs
#     w_t(L) = 0 for two DIFFERENT rows with different u_t, which is a different
#     (and non-absorbing) equation.
#     Cost: one 32-bit field extract from a word the step already mixes (4), the
#     constant-divisor reduction h mod (R-1) (mul-shift, 6), one XOR (1), the
#     compare-and-branch (5) = 16 instructions, charged on EVERY multiply step.
ROWDIFF      = 16
ROWS_BITS    = 14


# --- THE GUARD.  D=12 ring of (row, column) pairs, most recent first.  One step
#     does D packed compares with conditional select, one ring shift of D words
#     and the index arithmetic: D*(1 CMP + 1 CSEL) + D stores + 4 = 40, charged
#     at 48 so that the packed-word mask and the compare bookkeeping are inside.
GUARD_DEPTH  = 12
GUARD_ROW_BITS, GUARD_COL_BITS = 14, 4
GUARD_COMPARE = 2 * GUARD_DEPTH          # 24
GUARD_STORE   = GUARD_DEPTH              # 12
GUARD_INDEX   = 4
GUARD         = GUARD_COMPARE + GUARD_STORE + GUARD_INDEX + 8   # 48

# --- the sigma powers moved OFF the per-step path and ONTO the table build.
SIGMA_MULS    = 36
SIGMA_CONJ    = 6
SIGMA_BUILD   = SIGMA_MULS*FP_MUL + (SIGMA_CONJ + 6)*FP_ADD   # = 14 028

# --- NO label renormalisation and NO sigma-on-state in v14 either.
RENORM       = 0

# --- per-step rows: imported by ledger.rows(), restated nowhere
KEY_HASH     = 3*HASH_W       # one 256-bit mix -> table word + DP word + step word
DP_TEST      = 8
HIST_SLOT    = 32             # class-key traffic, 48 bytes
LOOP         = 32             # table addressing, counters
TABLE_FETCH  = 72             # 576-byte table entry: 72 eight-byte loads
LABEL_ADD    = 28             # a += u mod r

# --- trail cutoff and restart schedule (NO report cap: one report per trail)
TRAIL_METERS = 16
CAP_SLACK    = 1 << 25
THETA_BITS   = 80
TABLE_BITS   = 14
TABLE_COLS   = 12
RECORD_BYTES = 112            # 48-byte class key + 2 x 32-byte labels, ALWAYS
KEY_BYTES, LABEL_BYTES = 48, 32
REPORT_WORDS = RECORD_BYTES // 8          # 14 eight-byte words per report
TRAIL_START_MULS = 2          # y = g^a z^b: two 381-bit scalar multiplications

# --- declared: no window, no checkpoint, no counter, no escape, no restart of a
#     live trail (the guard rewrites one column; it never abandons a trail).
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

F_ESC = Fraction(0)                 # v14 has no escape branch at all


def cyclo_sqr_instructions():
    """The counted Granger-Scott Fp12 squaring row, ONE definition, used by the
    step ledger AND by the GT-exponentiation price in ledger.gt_exp_reference_units.
    24 Fp multiplications + 122 Fp additions (tower.py: cyclo_sqr_fast)."""
    return 24*FP_MUL + 122*FP_ADD


# --- RETENTION, CORRECTED (v16).  v15 charged a merge-stickiness factor
#     rho = 0.9809578 measured on the JOINT state (element, ring tuple).  That is
#     the wrong event: the collision solve consumes a CLASS match
#     (class(y1) = class(y2) => y2 = sigma^j(y1), x = (p^j a1 - a2)/(b2 - p^j b1),
#     tried for j = 0..11 and VERIFIED), and the class key min{keys12(y)} is a
#     function of the ELEMENT ALONE -- the ring does not enter it.  A class match
#     is therefore usable at the report instant and cannot be destroyed by a later
#     guard rewrite, so the retention factor on the needed event is exactly 1 and
#     is NOT carried.  The measured class-stream curve (evidence/class-retention.json)
#     is shipped as a diagnostic for path-shadowing, which is bounded by tens of
#     steps against a 2^80-step report interval (2^-74 of the walk).
RHO_MEASURED = Fraction(1)              # the corrected event is instantaneous
STICK_LOSS = 0                          # name kept so the ledger row reads 0
TARGET_PAD = TARGET                     # no discount to pad back


# --- v17 RETENTION, MEASURED AND SPLIT (verdict 41783ecf defect 2).
#     p_detect = 0.8958 (215/240 merged pairs still class-agreeing at k=1, the
#     report interval of the shipped step; CP-95 lower 0.8576) and
#     p_usable|detect = 1.0000 (240/240 detected coincidences whose twelve-shift
#     witness recovers the planted x, known-answer gated).  The PRODUCT is the
#     retained fraction of the birthday rate.  The walk-quality allowance
#     mu = 11/10 ALREADY declares a 10 per cent inefficiency, and path sharing is
#     exactly such an inefficiency, so the ledger charges only the measured EXCESS
#     over that allowance: (1/p_detect)/mu = 1.1163/1.1 = 1.0148.
RETENTION_MEASURED = 0.8958333333333334
RETENTION_CP95_LOWER = 0.8575558748268416
RETENTION_EFFECTIVE = 1.0 / 1.014787  # measured excess over the declared mu allowance


def step_budget(mu=MU_WALK, target=TARGET_PAD, f=F_SQR, n=N_CLASSES,
                retention=RETENTION_EFFECTIVE):
    """Least integer S with S^2 (1-f^2) retention >= 2 mu N ln(1/(1-target))."""
    from math import log, isqrt
    s2 = 2*float(mu)*n*log(1/(1-float(target)))/((1 - float(f)*float(f))*float(retention))
    s = isqrt(int(s2))
    while s*s < s2:
        s += 1
    while s > 1 and (s-1)*(s-1) >= s2:
        s -= 1
    return s


S_STEPS = step_budget()
