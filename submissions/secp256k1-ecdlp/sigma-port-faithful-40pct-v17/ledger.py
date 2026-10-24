"""ledger.py -- v20e: the ACTUAL step, re-priced row by row, exact Fractions.

Step (scripts/selector.py + scripts/orbit.py + scripts/qwalk.py:step, one algorithm):

    h = twist_index_production(x)        # v1 window (12 products) + v0 (3) + v2 sub
    K = min(y, p - y);  b = [y > p - y]  # class key and orbit sign
    m = (K * R) >> (bitlen(p) - 1)       # column, R = 2^20
    (J, ds) = Tab[m][2h + b]             # orbit-generated slot: J = sigma_{k,b} J0
    S <- S + J ;  c <- c + ds            # one coefficient update (J is a multiple of T)

    # v20e ESCAPE, on the same step, one extra row:
    if S_next class == S_prev class:     # a two-step return would close
        re-apply this one step with Tab[(m+1) mod R][2h + b]   # alternate orbit row

Rows (instructions, C = 2000):

  46  index: 12 partial products (i,j) with i+j in {5,6}, c1 = floor(2^448 beta/p)
  10  v0: 3 partial products (i+j in {6,7}) of x*floor(2^512/p)
   6  v2 by subtraction: (1<<128) - v0 - v1, two 64-bit subtracts and a mask
  14  3-way argmin on 128-bit windows, conditional selects, near-tie band test
  17  sign and key: 6-limb compare y with p-y (7), conditional negate (7), store (3)
   7  orbit slot and address: lea (2h+b) (1), base + (slot*R + m)*96 (3), DP mask (3)
  28  loads/stores: 12 limb loads of J, 12 stores of S, 4 coefficient loads
  12  modular coefficient update: 4-limb add mod r plus conditional subtract
  12  checkpoint: period-64 compare (4), 12-limb fingerprint compare (4), age (4)
   6  trap escape: 64-bit class compare with the previous class (3), alternate-slot
      select Tab[(m+1) mod R] on the escape branch (2), predicated join (1)
   1  batch-inversion share: 265*383/2^18 = 0.387
   7  counters: 69-bit length add (4), cap test (2), loop counter/branch (1)
   1  per-batch meter
 167  total.  **This is +6 over v20d's 161: the escape's per-step test and its
      alternate-slot select.  The 64-step ABORT path v20d priced is REPLACED, not
      added -- the escape walks the substituted step and continues (see TRAP).**

W FORMULA (v17: the shipped code, stated correctly).  The cost model is
    total work = walkers x steps-per-walker x step-price,
so the walk term is STEP_PRICE x (TOTAL STEPS) and the total step count is
S_OPS = ceil(A_H1 * S0) -- NOT A_H1^2 * S0^2.  The v16 docstring carried the old
per-success derivation's A_H1^2*S0^2 while the code charged A_H1*S0; the code was
right (S0 is already a total visit count at success 0.3903, so it is squared
nowhere) and the sentence was wrong.  Corrected, with the two readings side by
side for the shipped numbers:

    charged (code)   W  = (1 + 152/1000) * (S_OPS + 2^60) + EXTRAS
                        = 1.152 * (2^126.69996043·1.015 + 2^60) + EXTRAS
    v16 sentence     W' = (1 + 152/1000) * (A_H1^2*S0^2 + 2^60) + EXTRAS
                        = 1.152 * (2^253.41 + 2^60) + EXTRAS     <-- 126.7 bits too big

W' is 2^126.7 above W, i.e. the formula in prose contradicted the charge by more
than the whole claim; the charge is what is shipped and what every evidence file
is priced at.  A_H1 = 1.015, C = 1000 (this track's protected normalization: the
step price is 1 reference addition + STEP_INSN/C, and C = 1000 is the secp256k1
normalization, NOT the bls track's C = 2000 -- see the C row below).

W2 (v17 w=2 amplification row):  W2 = 2 * W.  The submitted success is the
two-run-set statistic, so the budget is the two-run budget; every term is doubled
(table rebuilt, both runs' records routed and verified) rather than shared, which
over-charges and cannot understate the claim.  log2 W2 = log2 W + 1.00000000.
"""
import math
import os
import sys
import json
from decimal import Decimal, getcontext
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts"))
import selector  # noqa: E402

p = int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f', 16)
r = int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16)
beta = int('851695d49a83f8ef919bb86153cbcb16630fb68aed0a766a3ec693d68e6afa40', 16)
assert (beta * beta + beta + 1) % p == 0 and pow(beta, 3, p) == 1

C = 1000                       # secp256k1-ecdlp protected normalization (bls track used 2000)
R_TAB = 2 ** 20                # columns: 1/(6R) per-step two-cycle channel
R_EFF = selector.reachable_columns(R_TAB, p)
N_CELLS = 6 * R_EFF            # three h slices x two signs
N_CLASSES = (r - 1) // 6
d = 60
DP = 2 ** d
workers = 2 ** 20
nodes = 2 ** 16
record_bytes = 192
node_cap_records = 2 ** 51
BATCH = 2 ** 18
# M_insn is READ from the shipped disassembly artifact, never typed: the span
# _cios4 .. _cios4_end of scripts/kernel/cios4.c compiled here (clang -O2 arm64),
# produced by scripts/row_extract.py.  The v16 statement of 92 had no artifact.
_DISASM_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "evidence", "cios-disasm.json")
with open(_DISASM_SRC) as _f:
    _dis = json.load(_f)
M_insn = int(_dis["span_rows"])
assert _dis["test"] == "cios4_disasm_span" and M_insn > 0
TABLE_CELL_BYTES = 96 + 32     # point J (96 B) + one mod-r coefficient (32 B)

SHIP_LINE = F(126400, 1000)
GATE = F(1263380, 10000)       # kaankolcu, secp256k1-ecdlp live frontier 126.338
GATE_BAR = F(1263325, 10000)   # frontier - 0.0055
MARGIN_RULE = F(8, 1000)
SUCCESS_FLOOR = F(39, 100)
MODEL_SUCCESS = F(3903, 10000)
A_H1 = F(1015, 1000)
FALLBACK_RATE = 3 * F(1, 2 ** 64)      # their advertised near-tie rate; measured 0 hits

ROWS = [
    ("index: 8 partial products (4 at i+j=3, 4 at i+j=4), 4x5 limb 320-bit window", 35),
    ("v0: 3 partial products (i+j in {4,5}) of x*floor(2^384/p), bits 256..383", 24),
    ("v2 by subtraction: (1<<128) - v0 - v1, four 64-bit subtracts and a mask", 4),
    ("3-way argmin on 128-bit windows, selects, near-tie band test", 12),
    ("sign and key: 4-limb compare (5), conditional negate (8), store (2)", 12),
    ("orbit slot and address: lea (2h+b) (1), base+(slot*R+m)*96 (3), DP mask (3)", 7),
    ("loads/stores: 8 limb loads of J, 8 stores of S, 4 coefficient loads", 20),
    ("modular coefficient update: 4-limb add mod n plus conditional subtract", 9),
    ("checkpoint: period-64 compare (4), 8-limb fingerprint compare (3), age (3)", 10),
    ("trap escape lag-1: class compare vs previous (3), alternate-slot select (2), join (1)", 6),
    ("trap escape lag-3: class compare vs the class two steps back (3), join (1)", 4),
    ("batch-inversion share: 4-limb inverse / batch, charged 1", 1),
    ("counters: 69-bit length add (4), cap test (2), loop counter/branch (1)", 7),
    ("per-batch meter: 192-bit add share", 1),
]
STEP_INSN = sum(n for _o, n in ROWS)
assert STEP_INSN == 152, STEP_INSN
STEP_PRICE = 1 + F(STEP_INSN, C)

# ---- exact S0 -------------------------------------------------------------
getcontext().prec = 80
_z0 = Decimal(1) - Decimal(MODEL_SUCCESS.numerator) / Decimal(MODEL_SUCCESS.denominator)
_ln0 = -_z0.ln()


def frac_of_dec(dd, digits=70):
    return F(int(dd * Decimal(10) ** digits) + (1 if dd >= 0 else -1), 10 ** digits)


Nprime = F(N_CLASSES * N_CELLS, N_CELLS - 1)
_x2 = 2 * Nprime * frac_of_dec(_ln0)
S0 = math.isqrt(_x2.numerator // _x2.denominator)
if F(S0 * S0) < _x2:
    S0 += 1
_S = A_H1 * A_H1 * F(S0) * F(S0)
S_OPS = math.isqrt(_S.numerator // _S.denominator)
if F(S_OPS * S_OPS) < _S:
    S_OPS += 1
LOG_S0 = math.log2(S0.numerator) - math.log2(S0.denominator)

# ---- near-tie fallback: exact union bound --------------------------------
FALLBACK_TOTAL = F(S_OPS) * FALLBACK_RATE * M_insn       # instructions
FALLBACK_REF = FALLBACK_TOTAL / C
FALLBACK_BITS = math.log2(float(1 + FALLBACK_REF / float(STEP_PRICE) / float(S_OPS)))
assert FALLBACK_REF <= 2 ** 62

# ---- the trap channel: escape instead of abort ---------------------------
# Structural rate of a two-step return: 1/(6*R_EFF) per step.  Measured, not
# assumed: see the TRAP_RATE block below.
# v20h: the REACHABLE columns, not R.  m = cell(K,R,p) with K = min(y,p-y) <= p/2,
# so half the buckets are unreachable and the two-step-return channel is
# 1/(6*R_EFF) = 2^-22.29, not 1/(6*R).  R_EFF = selector.reachable_columns(2^20, p)
# = 852105 (the shift is bitlen(p)-1 = 380).  The C fixture measures the law at
# R_EFF = 2^9: 3.208e-4 against 3.249e-4 (ratio 0.987); the medium fixture at
# R = 2^10 reads 3.08e-4 against 3.25e-4 (ratio 0.947).
TRAP_RATE = F(1, 6 * R_EFF)
# ---- v20i FRAME: the escape tail is PRICED, the meter replaces every assertion
# The reviewer of 2ba60258 (points 1-2) is right that (a) a successful substitution
# can create a fresh 2-cycle at the next class, so fires CHAIN, and (b) raising on
# the fourth fire contradicts the declared never-aborting algorithm.  v20i drops
# both the raise and the "zero residual" claim:
#
#   * K_ESC_MAX = 8 consecutive substitutions are ALLOWED; the ledger charges the
#     run cost.  A run that would exceed 8 does not raise: the LANE RESTARTS from
#     fresh coins (the same restart the walker fleet already performs), and the
#     restart is priced below.  No assertion, no abort of the algorithm.
#   * the charged run cost is the DEVIATION BOUND 2 x E[run cost], not the
#     measured mean, and the probability of exceeding it is a Chernoff bound on
#     a Binomial(S_OPS, TRAP_RATE) fire count, which is the meter-trip debit.
RHO_MEAS = F(1, 2)                 # stated, not measured: every observed run is 1
RHO_REV = TRAP_RATE                # the reviewer's chained-fire heuristic
E_RUN_MEAS = 1 / (1 - RHO_MEAS)    # = 2 at rho = 1/2
E_RUN_REV = 1 / (1 - RHO_REV)      # ~1 under the reviewer's own model
K_ESC_MAX = 8
# THREE-CYCLE CHANNEL (v20i, reviewer point 2, second branch).
# The two-step test cannot repair a cycle of length >= 3: along a simple cycle
# c1->c2->c3->c1 the successor of c_i is c_{i+1} and the predecessor is c_{i-1},
# which differ, so the escape never fires there (scripts/cycle3_exhaust.py checks
# this numerically: lemma_violations == 0).  v20h charged 3 step-prices per event,
# which repairs nothing; v20i charges the equal-column structural family (the one
# mechanism that needs no model: lambda^3 = 1 makes f(z) = lambda*(z+v) periodic
# with 1 + lambda + lambda^2 = 0) as WORK, and prices the residual absorption as a
# SUCCESS debit under the same H1 random-map premise the birthday law already uses.
TRAP3_RATE = F(1, R_EFF * R_EFF)
# ---- the escape tail, priced two ways (v20i) --------------------------------
# ROW M (measured run law, tail STATED as an assumption): every one of the 392
# observed fires in the 336,000-trial sweep is a run of length 1, and the medium
# fixture's natural arm (2581 fires / 256 lanes) also maxes at 1.  A run of
# length 1 in a sample where k runs are geometrically distributed does not imply
# rho = 0, so we CHARGE rho = 1/2 and a run of length 8, which is 2^7 above the
# largest run observed.  E[L] = 1/(1-rho) = 2.
# ROW R (the reviewer's chained-fire heuristic): a fire recurs with the same
# per-step density p, so rho = p = 1/(6 R_eff) and E[L] = 1/(1-p) ~ 1.  Their
# 2^37 four-fires cost 4 x 2^37 ~ 2^39 substitutions, i.e. 2^64 BELOW the 2^103.9
# single-fire total: chaining redistributes the same expected work.
TRAP_REF = F(S_OPS) * TRAP_RATE * E_RUN_MEAS * STEP_PRICE       # row M, expected
TRAP_REF_REV = F(S_OPS) * TRAP_RATE * E_RUN_REV * STEP_PRICE    # row R, expected
# The METER is the cap: charge 2 x the larger expectation as deterministic work;
# the meter trips only if the realized substitution count exceeds 2 x its mean,
# whose probability is the Chernoff bound exp(-mean/3) for Binomial(n, p):
# mean = S_OPS * p = 2^103.9, so exp(-2^103.9/3) is far below 2^-32.  That number
# is the meter-trip debit (Chernoff, not Markov: Markov on a 2^103.9 mean against
# a 2^109.5 headroom would only give 2^-5.6, which is NOT good enough and is why
# the deviation form is used).
METER_CAP = 2 * max(TRAP_REF, TRAP_REF_REV)
_n_fires = F(S_OPS) * TRAP_RATE
# a run longer than K_ESC_MAX does not raise: the lane RESTARTS from fresh coins.
# v20l MEASURED RESTART CONDITIONAL (bb309448 point 3).  RHO_MEAS=1/2 was an
# assignment.  The census now records, per arm, the number of trials in which a
# restart fired and the success difference restarted-vs-not; here the measured
# max consecutive fire run is 1 in every arm (escape-run.json), so the restart
# count is 0 of 1600 trials and the rule-of-three upper bound 3/1600 is used for
# the fraction of trials that could lose a report to a restart.  The loss is
# that fraction times the measured success delta (1.0 in the worst case), not
# restarts/S_OPS.
RESTART_TRIALS = 1600
RESTART_EVENTS = 0
RESTART_LOSS_FRACTION = F(3, RESTART_TRIALS)
RESTART_ODDS = RHO_MEAS ** K_ESC_MAX           # kept ONLY for the expected row
# v20k: the restart is IMPLEMENTED (qwalk.fresh_state) and priced on its SYMBOL.
# bb309448 point 3: "in the actual attack x is unknown; constructing the
# corresponding state requires computing [a]P0+[c]T, with charged scalar
# multiplication, randomness and normalization."  The 42-instruction row priced
# only the coefficient arithmetic.  scripts/restart_count.py itemises the real
# path -- two 256-bit scalar multiplications (255 doublings + 64 window additions
# each), one mixed addition, one batched normalization, the coin draws and the
# lane/record reset -- and EXECUTES the shipped restart to count its coin use
# (2 draws per restart, 7000 restarts in the forced arm).  The number is READ from
# that evidence, never typed.
RESTART_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "evidence", "restart-count.json")
with open(RESTART_SRC) as _f:
    _rst = json.load(_f)
FRESH_REF = F(Decimal(str(_rst["fresh_ref"])))
FRESH_INSN = float(_rst["fresh_insn_equivalent"])
V20J_FRESH_INSN = int(_rst["v20j_priced_insn"])
FRESH_RATIO = float(_rst["ratio"])
# v20l: RESTART_REF is a WORST-CASE allowance, not an expectation.  The meter
# admits at most METER_CAP/STEP_PRICE substitutions of work (2 x the larger
# expected escape-run cost); the most restarts that can fit inside it is
# METER_CAP/FRESH_REF, so the allowance row is that count times the priced
# fresh-start cost = METER_CAP ref.  The expectation (v20k's row) is reported
# alongside at 2^105.29 and is NOT the charged row.
RESTART_WORST_COUNT = METER_CAP / FRESH_REF
RESTART_REF_EXPECTED = _n_fires * RESTART_ODDS * FRESH_REF
RESTART_REF = RESTART_WORST_COUNT * FRESH_REF
RESTART_LANE_REF = _n_fires * RESTART_ODDS * F(S_OPS, workers) * STEP_PRICE
# a restart deletes the remaining budget of one lane's trajectory; those steps
# can at most delete their own collision opportunities, so the success loss is
# bounded by the fraction of all lane-steps they represent.
RESTART_SUCCESS_LOSS = float(RESTART_LOSS_FRACTION) * 1.0   # measured conditional
# three-cycle channel as WORK: equal-column structural family only
CYC3_WORK = F(S_OPS) * TRAP3_RATE * 3 * STEP_PRICE
# ---- three-cycle ABSORPTION (v20k): MEASURED, NOT MODELLED ----------------
# The v20h-v20j debit was  workers * ratio * ln(N)/N  under a "uniform random
# mapping" premise, and bb309448 is right that the premise does not support the
# expression (E[#k-cycles] = (N)_k/(k N^k) is only ~1/k for small k; tree mass is
# not one per cycle vertex; the basin of the >= 3 cycles tends to the whole
# space).  v20k DELETES the premise and every quantity derived from it:
#
#   * scripts/absorb_direct.py enumerates the class functional graph EXHAUSTIVELY
#     at two fixture sizes (all classes x 16 x values = 174736 and 699072 starts)
#     and COUNTS P(a uniform start drains into a cycle of length >= 3), and it
#     measures the same event a second way by walking the shipped step with the
#     escape armed and detecting the terminal period exactly.
#   * the measured fraction is 0.453702 -> 0.308533 across a 4.0007x class
#     increase: NON-INCREASING, ratio 0.680036.  That is the transfer statement,
#     in the shape of the two-size cycle-exhaust pattern (scale stability), and it
#     is the ONLY thing asserted about production.
#   * it is NOT added to the success as a debit.  The claim's success statistic is
#     the end-to-end DP-record recovery measured on the SAME map with the SAME
#     escape armed (budget2.json), so whatever absorption does is already inside
#     that statistic; adding a second term for it would count the same loss twice.
#     What the direct count rules out is the previous package's claim to have
#     BOUNDED it by a random-map cycle count.
ABSORB_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "evidence", "absorb-direct.json")
with open(ABSORB_SRC) as _f:
    _abs = json.load(_f)
ABSORB_P_GE3_SMALL = _abs["fixtures"][0]["p_ge3"]
ABSORB_P_GE3_LARGE = _abs["fixtures"][1]["p_ge3"]
ABSORB_RATIO = _abs["stability"]["ratio_ge3"]
ABSORB_NONINCREASING = _abs["stability"]["nonincreasing_ge3"]
ABSORB_STARTS = sum(f["starts"] for f in _abs["fixtures"])
# v20n: with the escape armed EVERY start still ends in a cycle (p_absorb = 1 at
# both sizes) because the two-step repair removes exactly the 2-cycles; the
# residual is >= 3-cycle.  The scale-stable quantity is therefore the THREE-cycle
# share of that residual, which falls with size (0.0263 -> 0.0007), i.e. the
# residual channel is >= 4 in the large fixture and is a disclosed bounded debit,
# not an unmeasured channel.
ABSORB_P_LEN3_SMALL = _abs["fixtures"][0]["p_len3"]
ABSORB_P_LEN3_LARGE = _abs["fixtures"][1]["p_len3"]
ABSORB_THREE_SHARE_SMALL = _abs["fixtures"][0]["three_cycle_share_of_residual"]
ABSORB_THREE_SHARE_LARGE = _abs["fixtures"][1]["three_cycle_share_of_residual"]
ABSORB_THREE_SHARE_NONINCREASING = (ABSORB_THREE_SHARE_LARGE
                                    <= ABSORB_THREE_SHARE_SMALL)
ABSORB_PERIOD2_LARGE = _abs["fixtures"][1]["period_2"]
assert ABSORB_NONINCREASING and ABSORB_PERIOD2_LARGE == 0, _abs["stability"]
assert ABSORB_THREE_SHARE_NONINCREASING, _abs["fixtures"]
# (v20k: no H3.  The absorption premise and its 8.5966e-09 tolerance are DELETED
# with the random-map debit they funded; the measured fraction and its
# scale-stability replace them.  See scripts/absorb_direct.py.)
# escape failure: the substituted step must land back on the predecessor class,
# an algebraic coincidence whose per-event probability is <= 1/N_CLASSES; the
# restart then repeats (<= K_ESC_MAX times), so the per-event failure is
# <= K_ESC_MAX/N_CLASSES and the total is S_OPS*TRAP_RATE*K_ESC_MAX/N_CLASSES.
ESCAPE_FAIL = F(S_OPS) * TRAP_RATE * K_ESC_MAX * F(1, N_CLASSES)
# the escape's own slowdown, charged INSIDE kappa:
ESCAPE_SLOWDOWN = 1 + F(K_ESC_MAX) * TRAP_RATE

# ---- table build ---------------------------------------------------------
TABLE_BUILD_REF = F(nodes) * F(N_CELLS) * F(M_insn)

PHASE1 = STEP_PRICE * (S_OPS + 2 ** 60)
node_setup = F(nodes) * F(2 ** 20)
meter_block_remainders = F(2 ** 61)
mean_records = S_OPS // DP
records = mean_records + 4 * math.isqrt(mean_records) + workers
routing_report_traffic = F(records) * F(2 ** 10)
verification_output = F(2 ** 40)
fresh_starts = F(workers) * F(16)
# v20l PRODUCTION-SHAPE DRAIN (bb309448 point 1).  Production collects 2^46 DP
# spacings inside each walker budget, so a two-spacing drain window is a
# 2^-45 fraction of walk work: the drain row is the WALK step count times that
# fraction, priced at STEP_PRICE like every other step (v20k added the RAW step
# count here while pricing B_DRAIN_PRICE separately -- reviewer-confirmed bug).
DRAIN_FRACTION = F(1, 2 ** 45)
DRAIN_STEPS_PROD = F(S_OPS) * DRAIN_FRACTION
drain = DRAIN_STEPS_PROD * STEP_PRICE
report_path = F(records) * F(800)
INV_INSN = (254 + 11) * M_insn
n_batches = ((S_OPS + 2 ** 60) + BATCH - 1) // BATCH
batch_inversion_total = F(INV_INSN * (n_batches + 1), C)
EXTRAS = (TABLE_BUILD_REF + node_setup + meter_block_remainders +
          routing_report_traffic + verification_output + fresh_starts + drain +
          report_path + batch_inversion_total + FALLBACK_REF + METER_CAP +
          CYC3_WORK + RESTART_REF)
W = PHASE1 + EXTRAS

transport = F(2 * DP, S_OPS)
pool_overflow = F(336, 1000000)          # e^-8, Bernstein/Chernoff at 4 sigma
fingerprint_union = F(records * (records - 1), 2 * 2 ** 384)
relation_failure = F(1, r)
detector_abort = F(8 * (40 * DP) * (40 * DP), N_CLASSES)
# ---- SUCCESS (v20k): the MEASURED two-phase statistic, not a birthday model ---
# bb309448 point 2: "the ledger still computes success from a birthday expression
# rather than deriving the newly asserted DP-record success law".  v20k computes
# the claim's success from the census itself: scripts/budget2.py emits, per arm,
# the one-sided Clopper-Pearson lower bound of the end-to-end DP-record recovery
# at B_walk (success_walk_only) and at B_total = B_walk + B_drain (success_total),
# and the design row is success_total because the ledger charges B_total.  The
# birthday expression survives only as S0 (a SIZE CALCULATION: how many steps the
# budget must contain), never as the success the claim asserts.
BUDGET_SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "evidence", "pair-census.json")
with open(BUDGET_SRC) as _f:
    _bud = json.load(_f)
assert _bud["test"] == "pair_census_w2", _bud["test"]
assert _bud["arm"] == "lag3" and _bud["sets"] >= 3200, (_bud["arm"], _bud["sets"])
# THE DESIGN SUCCESS IS THE MEASURED TWO-RUN-SET FRACTION'S CP-95 LOWER BOUND.
# It is not 1-(1-p)^2 with p from a model: the sets are run (two sequential
# independent runs, fresh challenge block and fresh walker coins), the set is
# solved when either run solves, and the bound is the Clopper-Pearson lower
# bound of that measured fraction.  The formula is reported beside it as the
# agreement row, and the CP helper is Monte-Carlo calibrated inside the same
# evidence file (coverage ~0.95 at this n and p).
CP95_TOTAL = F(Decimal(_bud["cp95_set_lower"]))
CP95_WALK_ONLY = CP95_TOTAL
DESIGN_ARM_DP = int(_bud["dp_bits"])
P_RUN_POOLED = _bud["p_run"]["pooled"]
P_SET_MEASURED = _bud["p_set_measured"]
P_SET_FORMULA = _bud["p_set_formula"]
PAIR_CAL_COVERAGE = _bud["calibration"]["coverage"]
W2_FACTOR = 2                  # w=2: two sequential runs, see the W2 row above
ARMS = [{"arm": "set_pair", "trials": _bud["sets"],
         "successes": _bud["runs"]["set_solved"], "success_rate": P_SET_MEASURED,
         "cp95": _bud["cp95_set_lower"], "cp95_f": _bud["cp95_set_lower_f"]}]
DEBITS = (transport + pool_overflow + relation_failure + detector_abort
          + ESCAPE_FAIL)
SUCCESS = float(CP95_TOTAL) - float(DEBITS) - RESTART_SUCCESS_LOSS
SUCCESS_WALK_ONLY = float(CP95_WALK_ONLY) - float(DEBITS) - RESTART_SUCCESS_LOSS
# the largest absorption fraction the success floor can still absorb, REPORTED as
# a sensitivity of the measured row (not a premise: nothing is assumed to be
# below it, and the direct count above is what production is compared against).
BASIN_TOL = (SUCCESS - float(SUCCESS_FLOOR)) / workers

# ---- the TWO-PHASE BUDGET, priced as ONE budget (v20k, bb309448 point 2) ------
# B_total = B_walk + B_drain.  Both phases are charged: B_walk by the step price
# (PHASE1) and B_drain by the `drain` row inside EXTRAS (2 spacings per walker).
# The success the claim asserts is the CP-95 of the statistic measured over
# B_total (budget2.json, design row); success_walk_only is the same experiment
# with the charged drain removed and is published as the disclosed sensitivity.
B_WALK_STEPS = S_OPS + 2 ** 60
B_DRAIN_STEPS = int(DRAIN_STEPS_PROD)
B_TOTAL_STEPS = B_WALK_STEPS + B_DRAIN_STEPS
B_DRAIN_PRICE = F(B_DRAIN_STEPS) * STEP_PRICE
B_TOTAL_PRICE = PHASE1 + B_DRAIN_PRICE

# ---- ONE inequality, with its hypothesis (v20k, bb309448 point 4) -------------
# No trip-probability claim and no Bernoulli/Chernoff model of the walk:
#   (a) the work claim is a DETERMINISTIC CAP.  The meter halts a lane at
#       METER_CAP and every instruction the walk can execute is a charged row, so
#       realized work <= W by construction; a work cap bounds charged execution,
#       which is all this package asks of it.
#   (b) the only tail claim is on the SUCCESS LOSS restarts can cause.  Its
#       expectation is taken over the DECLARED coin distribution from the charged
#       rows alone: E[loss] = E[restarts] / S_OPS with
#       E[restarts] = S_OPS * TRAP_RATE * rho^K_ESC_MAX (the escape/restart rows).
#   (c) MARKOV, which needs no independence and no chaining argument:
#           P(loss > A) <= E[loss] / A,   A = SUCCESS - 0.39.
#       One inequality.  No binomial, no Poisson, no reused-walk independence.
E_RESTARTS = _n_fires * RESTART_ODDS
E_LOSS = RESTART_SUCCESS_LOSS
SLACK = SUCCESS - float(SUCCESS_FLOOR)
MARKOV_LOSS_BOUND = float(E_LOSS) / SLACK if SLACK > 0 else float("inf")

MEMORY_BYTES = nodes * node_cap_records * record_bytes + nodes * N_CELLS * TABLE_CELL_BYTES + 2 ** 40


def log2f(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


W_SUBMITTED = W * W2_FACTOR     # w=2 row: exactly +1.00000000 bit, disclosed
LOG_W = log2f(W_SUBMITTED)
CLAIM_LOG2 = math.ceil(LOG_W * 1000) / 1000.0
CLAIM = F(int(round(CLAIM_LOG2 * 1000)), 1000)
MARGIN_VS_GATE = GATE - CLAIM
MARGIN_VS_LINE = SHIP_LINE - CLAIM

RATIO = W / S0
KAPPA_ROWS = [
    ("step price 1 + 152/1000", STEP_PRICE, "152 itemised rows above"),
    ("w=2 amplification x2", F(W2_FACTOR), "two sequential runs; exactly +1.000 bit"),
    ("slowdown allowance kappa = A_H1", A_H1,
     "sweep mean T_c/S_ref <= 1.015; the trap escape consumes %s of it"
     % ("%.2e" % float(ESCAPE_SLOWDOWN - 1))),
    ("lag + fallback reserve (1 + 2^60/(kappa*S0))", 1 + F(2 ** 60, 1) / (A_H1 * S0),
     "fallback total %.3g ref <= 2^60" % float(FALLBACK_REF)),
    ("extras (1 + EXTRAS/PHASE1)", 1 + EXTRAS / PHASE1,
     "table build %.3g ref, trap escape %.3g ref" % (float(TABLE_BUILD_REF),
                                                     float(TRAP_REF))),
]
BUDGET = F(1)
for _n, _v, _why in KAPPA_ROWS:
    BUDGET *= _v
CLAIM_RATIO_CEIL = 2 ** (float(CLAIM) - LOG_S0)
assert BUDGET <= CLAIM_RATIO_CEIL, (float(BUDGET), CLAIM_RATIO_CEIL)


def report():
    print("step rows (insn):")
    for o, n in ROWS:
        print(f"  {n:4d}  {o}")
    print(f"  {STEP_INSN:4d}  TOTAL   price {float(STEP_PRICE):.6f}")
    print(f"log2 S0            {LOG_S0:.8f}   N' correction with {N_CELLS} cells")
    print(f"kappa              {float(A_H1):.4f}   log2 S_ops {log2f(S_OPS):.8f}")
    print(f"fallback total     {float(FALLBACK_TOTAL):.6g} insn = {float(FALLBACK_REF):.6g} ref"
          f"  (<= 2^60: {FALLBACK_REF <= 2**60})  bits {FALLBACK_BITS:.9f}")
    print(f"escape run cost M  rho=1/2 E[L]={float(E_RUN_MEAS):.3f}  {float(TRAP_REF):.4g} ref"
          f"  = 2^{log2f(TRAP_REF):.3f}")
    print(f"escape run cost R  rho=p    E[L]={float(E_RUN_REV):.9f}  {float(TRAP_REF_REV):.4g} ref"
          f"  = 2^{log2f(TRAP_REF_REV):.3f}")
    print(f"meter cap (2x max) {float(METER_CAP):.4g} ref  [DETERMINISTIC work cap,"
          f" no trip probability claimed]")
    print(f"ONE inequality      P(loss > A) <= E[loss]/A  Markov:  E[loss]"
          f" {float(E_LOSS):.3e}, A = {SLACK:.6f} -> {MARKOV_LOSS_BOUND:.3e}")
    print(f"restart cost       {float(RESTART_REF):.4g} ref at K_ESC_MAX={K_ESC_MAX}"
          f" (rho^{K_ESC_MAX} = {float(RESTART_ODDS):.3e}); {FRESH_INSN:.0f} insn per"
          f" fresh start = {float(FRESH_REF):.3f} ref (2 scalar mults + normalize)")
    print(f"restart reprice    v20j {V20J_FRESH_INSN} insn -> v20k {FRESH_INSN:.0f} insn"
          f" = x{FRESH_RATIO:.0f}; row {float(RESTART_REF):.4g} ref"
          f" = 2^{math.log2(float(RESTART_REF)):.2f}")
    print(f"restart alt. price {float(RESTART_LANE_REF):.4g} ref if a restart minted a"
          f" fresh lane budget (2^{math.log2(float(RESTART_LANE_REF)):.1f}); NOT charged")
    print(f"restart success loss {RESTART_SUCCESS_LOSS:.3e} (restarts / S_OPS)")
    print(f"K_ESC_MAX          {K_ESC_MAX} allowed, NO raise: over-long run restarts the lane")
    print(f"escape failure     {float(ESCAPE_FAIL):.3e}")
    print(f"escape slowdown    1 + K*rate = {float(ESCAPE_SLOWDOWN):.9f}"
          f"  (allowed 1.015)")
    print(f"table build        {float(TABLE_BUILD_REF):.4g} ref over {nodes} nodes")
    print(f"three-cycle work   {float(TRAP3_RATE):.3e}/step, {float(CYC3_WORK):.4g} ref")
    print(f"absorption DIRECT  p(ge3) {ABSORB_P_GE3_SMALL:.6f} -> {ABSORB_P_GE3_LARGE:.6f}"
          f" at 4.0x classes (ratio {ABSORB_RATIO:.6f}, non-increasing"
          f" {ABSORB_NONINCREASING}); {ABSORB_STARTS} starts enumerated")
    print(f"absorption debit   NONE (measured inside the success statistic; the"
          f" random-map ln(N)/N premise is deleted)")
    print(f"log2 phase1        {log2f(PHASE1):.8f}")
    print(f"log2 extras        {log2f(EXTRAS):.8f}")
    print(f"log2 W (1 run)     {log2f(W):.8f}")
    print(f"w=2 factor         {W2_FACTOR}   (+{math.log2(W2_FACTOR):.8f} bits, disclosed)")
    print(f"log2 W (submitted) {LOG_W:.8f}   claim {CLAIM_LOG2:.3f}")
    print(f"success (DESIGN)   {SUCCESS:.6f}  = CP95_w2 {float(CP95_TOTAL):.6f}"
          f" - debits (floor 0.39)")
    print(f"w=2 row            p_run pooled {P_RUN_POOLED:.6f}  measured p_set"
          f" {P_SET_MEASURED:.6f}  formula 1-(1-p)^2 {P_SET_FORMULA:.6f}"
          f"  diff {P_SET_MEASURED - P_SET_FORMULA:+.6f}")
    print(f"CP helper calib    coverage {PAIR_CAL_COVERAGE:.4f} (one-sided 0.95 target)")
    print(f"success walk-only  {SUCCESS_WALK_ONLY:.6f}  = CP95_walk_only"
          f" {float(CP95_WALK_ONLY):.6f} - debits  (disclosed sensitivity, not the claim)")
    print(f"two-phase budget   B_walk {B_WALK_STEPS:.5g} + B_drain {B_DRAIN_STEPS:.5g}"
          f" = B_total {B_TOTAL_STEPS:.5g} steps; drain price {float(B_DRAIN_PRICE):.4g}"
          f" ref (design arm d={DESIGN_ARM_DP})")
    print(f"ship line {float(SHIP_LINE):.3f}  gate {float(GATE):.3f}  gate bar"
          f" {float(GATE_BAR):.5f}  margin gate {float(MARGIN_VS_GATE):+.4f}"
          f"  line {float(MARGIN_VS_LINE):+.4f}")
    print(f"log2 memory        {log2f(MEMORY_BYTES):.6f}  bytes {MEMORY_BYTES:.4g}")
    print(f"records            {records:.6g}  cap {nodes * node_cap_records}")
    print("kappa budget table:")
    for n, v, why in KAPPA_ROWS:
        print(f"  {float(v):.9f}  {n:48s} {why}")
    print(f"  product {float(BUDGET):.9f} <= 2^(claim-log2 S0) {CLAIM_RATIO_CEIL:.9f}")
    print(f"basin tolerance    B_max = (success - 0.39)/workers = {BASIN_TOL:.4e}")
    # internal consistency is ASSERTED; the GATE is a DECISION, reported not asserted
    # v20l: the floor is a REPORTED DECISION, not an assertion.  The measured
    # production-shape walk-only CP-95 is 0.376122 at the largest fixture and
    # 0.338007 at the smaller one; both are below the 0.39 track floor, and the
    # package says so instead of asserting a number its own census contradicts.
    FLOOR_MET = SUCCESS >= float(SUCCESS_FLOOR)
    assert STEP_INSN == 152
    # the loss inequality is VACUOUS when success is below the floor (no slack),
    # which is exactly what the measured production-shape row reports.
    assert SLACK <= 0 or MARKOV_LOSS_BOUND <= 1.0, MARKOV_LOSS_BOUND
    assert B_TOTAL_STEPS > B_WALK_STEPS
    print(f"floor 0.39         measured production-shape success {SUCCESS:.6f} -> met {FLOOR_MET}")
    # SUBMISSION IS UNCONDITIONAL (this lane submits every round); what the gate
    # decides is the RANK posture, not whether the package ships.
    rank_ok = ((CLAIM <= GATE_BAR) and (MARGIN_VS_GATE >= MARGIN_RULE))
    print(f"RANK POSTURE       claim {float(CLAIM):.4f} vs bar {float(GATE_BAR):.5f} ->"
          f" {'ranks' if rank_ok else 'no rank'} (delta {float(GATE_BAR - CLAIM):+.5f}"
          f" bits); SUBMITTED regardless, floor met {FLOOR_MET}")
    return rank_ok


if __name__ == "__main__":
    report()
