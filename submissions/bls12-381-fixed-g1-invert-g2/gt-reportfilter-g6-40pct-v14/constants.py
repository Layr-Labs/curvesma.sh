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

# --- trail cutoff, report cap, restart schedule
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


# --- MERGE-STICKINESS RETENTION.  The ring is per trail, so two trails that meet
#     at a common class can disagree about the guard for up to D steps, after
#     which their rings are identical (they recorded the same D rows and the same
#     applied columns) and the merge is exact for ever.  Per step the disagreement
#     probability is at most 2D/R (each ring holds at most D rows out of R), so
#     over the D-step reconciliation window the merge survives with probability
#        1 - 2*D*D/R  =  1 - 288/2^14  =  0.982421875
#     and that factor is charged against the success target below.
STICK_LOSS = Fraction(2*GUARD_DEPTH*GUARD_DEPTH, 1 << TABLE_BITS)
TARGET_PAD = TARGET / (1 - STICK_LOSS)


def step_budget(mu=MU_WALK, target=TARGET_PAD, f=F_SQR, n=N_CLASSES):
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
