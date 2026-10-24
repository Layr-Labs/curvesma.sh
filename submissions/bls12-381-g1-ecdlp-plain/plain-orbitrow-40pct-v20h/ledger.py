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

W = (1 + 167/2000) * (A_H1^2 * S0^2 + 2^60) + EXTRAS,  A_H1 = 1.015, C = 2000.
"""
import math
import os
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts"))
import selector  # noqa: E402

p = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624'
        '1eabfffeb153ffffb9feffffffffaaab', 16)
r = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
beta = int('5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a0002'
           '2e01fffffffefffe', 16)
assert (beta * beta + beta + 1) % p == 0 and pow(beta, 3, p) == 1

C = 2000                       # COST_MODEL.md: BLS12-381 G1 batched affine addition
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
M_insn = 383                   # blst_fp_mul, COST_MODEL.md
TABLE_CELL_BYTES = 96 + 32     # point J (96 B) + one mod-r coefficient (32 B)

SHIP_LINE = F(126266, 1000)
GATE = F(126274, 1000)         # mpjunior92 7b34f78d, live at this build (re-polled)
GATE_BAR = F(126268, 1000) + F(5, 10000)   # lowest - 0.004 - 0.0015 = 126.26850
MARGIN_RULE = F(8, 1000)
SUCCESS_FLOOR = F(39, 100)
MODEL_SUCCESS = F(3903, 10000)
A_H1 = F(1015, 1000)
FALLBACK_RATE = 3 * F(1, 2 ** 64)      # their advertised near-tie rate; measured 0 hits

ROWS = [
    ("index: 12 partial products, (i,j) with i+j in {5,6}, c1 = floor(2^448 beta/p)", 46),
    ("v0: 3 partial products (i+j in {6,7}) of x*floor(2^512/p); exact to <1 unit", 10),
    ("v2 by subtraction: (1<<128) - v0 - v1, two 64-bit subtracts and a mask", 6),
    ("3-way argmin on 128-bit windows, selects, near-tie band test", 14),
    ("sign and key: 6-limb compare (7), conditional negate (7), store (3)", 17),
    ("orbit slot and address: lea (2h+b) (1), base+(slot*R+m)*96 (3), DP mask (3)", 7),
    ("loads/stores: 12 limb loads of J, 12 stores of S, 4 coefficient loads", 28),
    ("modular coefficient update: 4-limb add mod r plus conditional subtract", 12),
    ("checkpoint: period-64 compare (4), 12-limb fingerprint compare (4), age (4)", 12),
    ("trap escape: class compare vs previous (3), alternate-slot select (2), join (1)", 6),
    ("batch-inversion share: 265*383/2^18 = 0.387, charged 1", 1),
    ("counters: 69-bit length add (4), cap test (2), loop counter/branch (1)", 7),
    ("per-batch meter: 192-bit add share", 1),
]
STEP_INSN = sum(n for _o, n in ROWS)
assert STEP_INSN == 167, STEP_INSN
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
# THREE-CYCLE CHANNEL (v20h, reviewer point 4).  The two-step test is blind to a
# cycle of length 3: for lambda^3 = 1 the affine map f(z) = lambda*(z+v) has
# f^3(z) = z because 1 + lambda + lambda^2 = 0.  That map needs the SAME column on
# all three steps (equal u, so the translation cancels), an event of density
# <= 1/R_EFF^2 per step; it is charged here.  Cycles with unequal columns are not
# that mechanism: they are rho-closures of the finite map, whose rate is bounded
# by walk^2/N (2^120/2^252.27 = 2^-132.27 per walk at production), so they cannot
# recur inside a production budget and are excluded structurally.  Measured:
# scripts/cycle3.py, 0 structural three-cycles in 16,777,216 lane-steps.
TRAP3_RATE = F(1, R_EFF * R_EFF)
K_ESC = 3
import qwalk                                    # noqa: E402
assert qwalk.K_ESC == K_ESC, (qwalk.K_ESC, K_ESC)   # the enforced cap is the charged one
# measured, not chosen: the natural escape arm's longest run along ONE walk is 1
# (medium: 2581 fires / 256 lanes; C fixture: 5382 fires at 2^26 classes and one
# full 2^20-step spacing).  The 1-in-7 forced stress arm reaches 3, and the ledger
# charges that stress maximum.  qwalk.step ENFORCES it: the walk raises if a lane
# exceeds K_ESC consecutive substitutions (v20h).
TRAP_REF = F(S_OPS) * TRAP_RATE * K_ESC * STEP_PRICE   # reference units, worst case
TRAP3_REF = F(S_OPS) * TRAP3_RATE * K_ESC * STEP_PRICE  # three-cycle channel, charged
TRAP_BITS = math.log2(float(1 + TRAP_REF / float(STEP_PRICE) / float(S_OPS)))
# escape failure: the substituted step must land back on the predecessor class,
# an algebraic coincidence whose per-event probability is <= 1/N_CLASSES; the
# escape then simply repeats (<= K_ESC times), so the per-event failure is
# <= K_ESC/N_CLASSES and the total is S_OPS*TRAP_RATE*K_ESC/N_CLASSES.
ESCAPE_FAIL = F(S_OPS) * TRAP_RATE * K_ESC * F(1, N_CLASSES)
# the escape's own slowdown, charged INSIDE kappa:
ESCAPE_SLOWDOWN = 1 + F(K_ESC) * TRAP_RATE

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
drain = F(workers) * F(2 * DP)
report_path = F(records) * F(800)
INV_INSN = (254 + 11) * M_insn
n_batches = ((S_OPS + 2 ** 60) + BATCH - 1) // BATCH
batch_inversion_total = F(INV_INSN * (n_batches + 1), C)
EXTRAS = (TABLE_BUILD_REF + node_setup + meter_block_remainders +
          routing_report_traffic + verification_output + fresh_starts + drain +
          report_path + batch_inversion_total + FALLBACK_REF + TRAP_REF + TRAP3_REF)
W = PHASE1 + EXTRAS

transport = F(2 * DP, S_OPS)
pool_overflow = F(336, 1000000)          # e^-8, Bernstein/Chernoff at 4 sigma
fingerprint_union = F(records * (records - 1), 2 * 2 ** 384)
relation_failure = F(1, r)
detector_abort = F(8 * (40 * DP) * (40 * DP), N_CLASSES)
_rho = 1 - math.exp(-float(A_H1) ** 2 * float(_ln0))
SUCCESS = _rho - float(transport + pool_overflow + relation_failure + detector_abort
                       + ESCAPE_FAIL)

MEMORY_BYTES = nodes * node_cap_records * record_bytes + nodes * N_CELLS * TABLE_CELL_BYTES + 2 ** 40


def log2f(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


LOG_W = log2f(W)
CLAIM_LOG2 = math.ceil(LOG_W * 1000) / 1000.0
CLAIM = F(int(round(CLAIM_LOG2 * 1000)), 1000)
MARGIN_VS_GATE = GATE - CLAIM
MARGIN_VS_LINE = SHIP_LINE - CLAIM

RATIO = W / S0
KAPPA_ROWS = [
    ("step price 1 + 167/2000", STEP_PRICE, "167 itemised rows above"),
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
    print(f"trap escape        rate {float(TRAP_RATE):.3e}/step, K_esc {K_ESC},"
          f" {float(TRAP_REF):.4g} ref, bits {TRAP_BITS:.9f},"
          f" failure {float(ESCAPE_FAIL):.3e}")
    print(f"escape slowdown    1 + K*rate = {float(ESCAPE_SLOWDOWN):.9f}"
          f"  (allowed 1.015)")
    print(f"table build        {float(TABLE_BUILD_REF):.4g} ref over {nodes} nodes")
    print(f"three-cycle charge {float(TRAP3_RATE):.3e}/step, {float(TRAP3_REF):.4g} ref")
    print(f"log2 phase1        {log2f(PHASE1):.8f}")
    print(f"log2 extras        {log2f(EXTRAS):.8f}")
    print(f"log2 W             {LOG_W:.8f}   claim {CLAIM_LOG2:.3f}")
    print(f"success            {SUCCESS:.6f}   (floor 0.39)")
    print(f"ship line {float(SHIP_LINE):.3f}  gate {float(GATE):.3f}  gate bar"
          f" {float(GATE_BAR):.5f}  margin gate {float(MARGIN_VS_GATE):+.4f}"
          f"  line {float(MARGIN_VS_LINE):+.4f}")
    print(f"log2 memory        {log2f(MEMORY_BYTES):.6f}  bytes {MEMORY_BYTES:.4g}")
    print(f"records            {records:.6g}  cap {nodes * node_cap_records}")
    print("kappa budget table:")
    for n, v, why in KAPPA_ROWS:
        print(f"  {float(v):.9f}  {n:48s} {why}")
    print(f"  product {float(BUDGET):.9f} <= 2^(claim-log2 S0) {CLAIM_RATIO_CEIL:.9f}")
    assert SUCCESS >= float(SUCCESS_FLOOR)
    assert MARGIN_VS_GATE >= MARGIN_RULE
    assert CLAIM <= GATE_BAR
    assert STEP_INSN == 167


if __name__ == "__main__":
    report()
