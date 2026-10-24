"""ledger_v11.py -- p256-ecdlp v11: the plain-lineage orbit-row walk PORTED to P-256.

WHAT IS PORTED.  This package is a retarget, not a new family.  It carries the
judge-validated plain-track construction (bls12-381-g1-ecdlp-plain ticket 5ecbf15,
verdict-valid at 126.2364) through the secp256k1 retarget (sigma-port-transfer-v18,
verdict cd18d785: every architecture point ACCEPTED -- orbit-row walk, refresh
store, lag-3 escape, metered ledger, selfaudit -- with only the transfer premise, a
stale restart denominator and a withdrawn kernel price held open) to NIST P-256.

THE ONE STRUCTURAL CHANGE, AND WHY IT MUST BE MADE.  The lineage's class map is the
mu_6 orbit quotient: it exists because BLS12-381 G1 and secp256k1 are j = 0 curves
(a = 0), so x -> beta*x is a curve automorphism and six points share one class key.
NIST P-256 is a = -3, b != 0 (j = 1728): x -> beta*x sends y^2 = x^3 - 3x + b to
y^2 = x^3 - 3*beta^2*x + b, which is P-256 only for beta = +-1, and (x, y) ->
(-x, i*y) needs b = 0.  P-256's automorphism group is therefore {+1, -1} and the
retarget's class map is canon(c) = min(c, q - c) -- the SAME code with lam = 1, so
the orbit min/max degenerates to the signed pair and e in {+1, -1}.  Nothing in the
step, the escape, the recovery or the ledger form changes; the class space becomes
N = (r-1)/2 instead of (r-1)/6, which is the +0.79-bit cost of the retarget and is
charged in S0, not argued.

W FORM (unchanged from the accepted lineage).  The cost model is
    total work = walkers x steps-per-walker x step-price,
so the walk term is STEP_PRICE x (TOTAL STEPS):

    W = (1 + STEP_INSN/C) * (S_OPS + 2^60) + EXTRAS

with C = 1000 (this track's protected normalization) and STEP_INSN itemised below.
No amplification factor: the success asserted is a MEASURED single-run Clopper-
Pearson lower bound, so W is the one-run budget (no +1.000 bit).

S_OPS IS THE BUDGET THE CENSUS ACTUALLY RAN.  evidence/census-design.json records
visits_over_sqrt(N) = x_meas for the arm whose CP-95 is claimed, so

    S_OPS = ceil(A_H1 * x_meas * sqrt(N'))    ,  N' = N*N_CELLS/(N_CELLS-1)

i.e. the executed visits times the kappa slowdown allowance.  Nothing is charged
below what was executed.
"""
import json
import math
import os
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "scripts"))
import selector  # noqa: E402

p = int('ffffffff00000001000000000000000000000000ffffffffffffffffffffffff', 16)
r = int('ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551', 16)
CURVE_A = int('ffffffff00000001000000000000000000000000fffffffffffffffffffffffc', 16)
CURVE_B = int('5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b', 16)
GEN = (int('6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296', 16),
       int('4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5', 16))
assert CURVE_A == p - 3 and 0 < CURVE_B < p
_gx, _gy = GEN
assert (_gy * _gy - (_gx ** 3 + CURVE_A * _gx + CURVE_B)) % p == 0
# NO NONTRIVIAL AUTOMORPHISM -- checked, not asserted in prose.  A nontrivial cube
# root of unity beta in F_p (p = 1 mod 3, so they exist) would preserve the curve
# only if a*(beta^2 - 1) = 0; with a = -3 that needs beta = +-1.
_BETAS = []
for _g in range(2, 400):
    _v = pow(_g, (p - 1) // 3, p)
    if _v != 1 and _v not in _BETAS:
        _BETAS.append(_v)
assert _BETAS and any((_b * _b * _b) % p == 1 and _b != 1 for _b in _BETAS)
assert all((CURVE_A * (b * b - 1)) % p != 0 for b in _BETAS)
assert (r - 1) % 2 == 0
# j-INVARIANT (verdict 27d28f66 point 4): the earlier text called P-256 j = 1728; with
# a != 0 AND b != 0 the correct statement is j NOT IN {0, 1728}.  j = 1728*4a^3/(4a^3+27b^2).
_J_NUM = 1728 * 4 * pow(CURVE_A, 3, p) % p
_J_DEN = (4 * pow(CURVE_A, 3, p) + 27 * pow(CURVE_B, 2, p)) % p
assert _J_DEN % p != 0
J_INV = _J_NUM * pow(_J_DEN, -1, p) % p
assert J_INV not in (0, 1728) and (CURVE_A % p) != 0 and (CURVE_B % p) != 0

C = 1000                       # p256-ecdlp protected normalization
R_TAB = 2 ** 20
R_EFF = selector.reachable_columns(R_TAB, p)
N_CLASSES = (r - 1) // 2       # orbit size 2: the retarget's cost, charged in S0
N_CELLS = 2 * R_EFF
d = 60
DP = 2 ** d
workers = 2 ** 20
nodes = 2 ** 16
record_bytes = 192
node_cap_records = 2 ** 53   # v11r: 2^69 slots >= the 2^68.20 reports the budget can emit
BATCH = 2 ** 18
M_insn = 97                    # STRUCTURAL count (no disassembly claimed anywhere)
TABLE_CELL_BYTES = 96 + 32

GATE = F(127237678, 1000000)   # aburan28 6d5e2adb, p256-ecdlp live frontier 127.237678
SUCCESS_FLOOR = F(39, 100)
MODEL_SUCCESS = F(3903, 10000)
A_H1 = F(1015, 1000)

# ---- p256 STEP RECOUNT ----------------------------------------------------
# The walk step IS one reference group operation (S <- S + J, affine addition: six
# field multiplications with a batched-inversion share) and that operation is the
# leading 1 in STEP_PRICE.  Policy: "Do not double-count a reference operation and
# its implementation", so the addition's field muls are not re-charged -- EXCEPT the
# P-256 reduction premium, because the borrowed reference unit is a secp256k1
# addition price and this track's policy says to account actual P-256 work in
# instructions.  Solinas reduction of a 512-bit product costs 8 mul-adds plus
# add/carry work against secp256k1's 4 mul-add pseudo-Mersenne fold: +30 insn per
# field multiplication, six of them = 180.
ROWS = [
    ("class key: 4-limb compare x vs p-x (4), masked 4-limb subtract (8), select (2)", 14),
    ("column and orbit slot: 64-bit shift by constant (1), table base + index (2)", 3),
    ("sign bit: y low-limb test and key pack (3)", 3),
    ("P-256 REDUCTION PREMIUM over the borrowed reference addition: 6 field muls x "
     "+30 insn (Solinas: 8 mul-adds + 20 add/carry, vs secp256k1's 4 mul-add fold)", 180),
    ("loads/stores: 4 limb loads of J, 4 stores of S, 4 coefficient loads (20)", 20),
    ("coefficient update mod n: 4-limb add plus conditional subtract", 9),
    ("checkpoint: period-64 compare (4), 4-limb fingerprint compare (3), age (3)", 10),
    ("trap escape lag-1: class compare vs previous (3), alternate-slot select (2), join (1)", 6),
    ("trap escape lag-3: class compare vs the class two steps back (3), join (1)", 4),
    ("DP fingerprint: 64-bit mix (4) and threshold compare (2)", 6),
    ("counters: 64-bit length add (4), cap test (2), loop counter/branch (1)", 7),
    ("per-batch meter: 192-bit add share", 1),
]
STEP_INSN = sum(n for _o, n in ROWS)
assert STEP_INSN == 263, STEP_INSN
STEP_PRICE = 1 + F(STEP_INSN, C)

# ---- S0 and the executed budget -------------------------------------------
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
SHAPE_BIRTHDAY = math.sqrt(2 * float(frac_of_dec(_ln0)))     # x at p = 0.3903


def load(rel):
    with open(os.path.join(HERE, rel)) as f:
        return json.load(f)


DESIGN = load("evidence/census-design.json")
LOW = load("evidence/census-low.json")
PAIR = load("evidence/pair-census.json")
ABS = load("evidence/absorb-direct.json")

FIXTURE_R = 16384              # RCOLS of the shipped census command (scripts/census_v20n.py)
X_MEAS = float(DESIGN["visits_over_sqrtN"])
SQRT_NP = math.sqrt(float(Nprime))
S_OPS = int(math.ceil(float(A_H1) * X_MEAS * SQRT_NP))
LOG_S0 = math.log2(S0.numerator) - math.log2(S0.denominator)

ARM_DESIGN = DESIGN["arms"]["escaped_lag3"]
CP95_TOTAL = F(Decimal(ARM_DESIGN["cp95_success"]))
CP95_LOW = F(Decimal(LOW["arms"]["escaped_lag3"]["cp95_success"]))
ORDER_SPAN = float(DESIGN["classes"]) / float(LOW["classes"])
DESIGN_ARM_DP = int(DESIGN["dp_bits"])
PAIR_CP = F(Decimal(PAIR["cp95_set_lower"]))
PAIR_SETS = int(PAIR["sets"])

# ---- trap channel, escape cost, restart ----------------------------------
TRAP_RATE = F(1, 2 * R_EFF)          # two-step return: P(e = -1) = 1/2, times 1/R_eff
RHO_MEAS = F(1, 2)
E_RUN_MEAS = 1 / (1 - RHO_MEAS)
K_ESC_MAX = 8
TRAP_REF = F(S_OPS) * TRAP_RATE * E_RUN_MEAS * STEP_PRICE
METER_CAP = 2 * TRAP_REF
TRAP3_RATE = F(1, R_EFF * R_EFF)
CYC3_WORK = F(S_OPS) * TRAP3_RATE * 3 * STEP_PRICE

# fresh start, priced structurally IN THIS FILE (2 scalar mults + normalize):
# 255 doublings at 2M+5S and 64 window additions at 11M, 4-limb, 4 insn/product
DBL_INSN = 7 * 4 * 24 + 20      # 7 field muls x 16 products x 4 insn, + adds
ADD_INSN = 11 * 4 * 24 + 20
SCALAR_MULT_INSN = 255 * DBL_INSN + 64 * ADD_INSN
NORMALIZE_INSN = 3 * 152
FRESH_INSN = 2 * SCALAR_MULT_INSN + ADD_INSN + NORMALIZE_INSN + 120
FRESH_REF = F(FRESH_INSN, C)
# v11r (verdict 27d28f66 point 1, restart pooling): the lag-1 and lag-3 arms are
# MATCHED pairs sharing starts and the lag-3 trigger fires 0 times on this orbit-2
# walk, so the two blocks are not 6,400 independent restart trials.  The denominator
# is the lag-1 armed block alone.
RESTART_TRIALS = int(DESIGN["arms"]["escaped_lag1"]["trials"])
RESTART_EVENTS = int(DESIGN["arms"]["escaped_lag1"]["restarts"])
assert RESTART_EVENTS == 0, RESTART_EVENTS
RESTART_LOSS_FRACTION = F(3, RESTART_TRIALS)
RESTART_ODDS = RHO_MEAS ** K_ESC_MAX
RESTART_REF = (METER_CAP / FRESH_REF) * FRESH_REF     # worst-case allowance
RESTART_SUCCESS_LOSS = float(RESTART_LOSS_FRACTION)

# ---- table build and the rest of EXTRAS ----------------------------------
TABLE_BUILD_REF = F(nodes) * F(N_CELLS) * F(M_insn)
node_setup = F(nodes) * F(2 ** 20)
meter_block_remainders = F(2 ** 61)
mean_records = S_OPS // DP
records = mean_records + 4 * math.isqrt(mean_records) + workers
routing_report_traffic = F(records) * F(2 ** 10)
verification_output = F(2 ** 40)
fresh_starts = F(workers) * F(16)
DRAIN_FRACTION = F(1, 2 ** 45)
drain = F(S_OPS) * DRAIN_FRACTION * STEP_PRICE
report_path = F(records) * F(800)
INV_INSN = (254 + 11) * M_insn
n_batches = ((S_OPS + 2 ** 60) + BATCH - 1) // BATCH
batch_inversion_total = F(INV_INSN * (n_batches + 1), C)
FALLBACK_RATE = 3 * F(1, 2 ** 64)
FALLBACK_REF = F(S_OPS) * FALLBACK_RATE * M_insn / C
EXTRAS = (TABLE_BUILD_REF + node_setup + meter_block_remainders +
          routing_report_traffic + verification_output + fresh_starts + drain +
          report_path + batch_inversion_total + FALLBACK_REF + METER_CAP +
          CYC3_WORK + RESTART_REF)
PHASE1 = STEP_PRICE * (S_OPS + 2 ** 60)
W = PHASE1 + EXTRAS

# ---- success: the measured single-run CP-95, minus priced debits ----------
transport = F(2 * DP, S_OPS)
pool_overflow = F(336, 1000000)
fingerprint_union = F(records * (records - 1), 2 * 2 ** 384)
relation_failure = F(1, r)
detector_abort = F(8 * (40 * DP) * (40 * DP), N_CLASSES)
ESCAPE_FAIL = F(S_OPS) * TRAP_RATE * K_ESC_MAX * F(1, N_CLASSES)
DEBITS = transport + pool_overflow + relation_failure + detector_abort + ESCAPE_FAIL
SUCCESS = float(CP95_TOTAL) - float(DEBITS) - RESTART_SUCCESS_LOSS
FLOOR_MET = SUCCESS >= float(SUCCESS_FLOOR)
SLACK = SUCCESS - float(SUCCESS_FLOOR)
E_LOSS = RESTART_SUCCESS_LOSS
MARKOV_LOSS_BOUND = (E_LOSS / SLACK) if SLACK > 0 else float("inf")

MEMORY_BYTES = nodes * node_cap_records * record_bytes + nodes * N_CELLS * TABLE_CELL_BYTES + 2 ** 40


def log2f(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


# SENSITIVITY (declared, not claimed): if a P-256 addition must be charged as TWO
# reference units (C=1000 is a normalization, not an implementation bound -- verdict
# point 2), the step price becomes 2 + 263/1000 and the claim moves by
# log2(2.263/1.263) bits.
STEP_PRICE_2REF = 2 + F(STEP_INSN, C)
REPORT_SLOTS = nodes * node_cap_records
LOG_W = log2f(W)
CLAIM_2REF = math.ceil((LOG_W + math.log2(float(STEP_PRICE_2REF / STEP_PRICE))) * 1000) / 1000.0
CLAIM_LOG2 = math.ceil(LOG_W * 1000) / 1000.0
CLAIM = F(int(round(CLAIM_LOG2 * 1000)), 1000)
MARGIN_VS_GATE = GATE - CLAIM

KAPPA_ROWS = [
    ("step price 1 + 263/1000", STEP_PRICE, "263 itemised rows above"),
    ("slowdown allowance kappa = A_H1", A_H1,
     "escape consumes %.2e of it" % float(K_ESC_MAX * TRAP_RATE)),
    ("lag + fallback reserve", 1 + F(2 ** 60, 1) / (A_H1 * S0),
     "fallback total %.3g ref" % float(FALLBACK_REF)),
    ("extras (1 + EXTRAS/PHASE1)", 1 + EXTRAS / PHASE1,
     "table build %.3g ref, escape %.3g ref" % (float(TABLE_BUILD_REF), float(TRAP_REF))),
]
BUDGET = F(1)
for _n, _v, _why in KAPPA_ROWS:
    BUDGET *= _v
CLAIM_RATIO_CEIL = 2 ** (float(CLAIM) - LOG_S0)
assert BUDGET <= CLAIM_RATIO_CEIL, (float(BUDGET), CLAIM_RATIO_CEIL)


def report():
    print("P-256 step rows (insn):")
    for o, n in ROWS:
        print("  %4d  %s" % (n, o))
    print("  %4d  TOTAL   price %.6f" % (STEP_INSN, float(STEP_PRICE)))
    print("class space        N = (r-1)/2 = 2^%.6f  (orbit size 2: no P-256 automorphism)"
          % math.log2(N_CLASSES))
    print("log2 S0            %.8f  (birthday x = %.6f)" % (LOG_S0, SHAPE_BIRTHDAY))
    print("executed shape     x_meas = %.6f, S_OPS %.6g log2 %.8f"
          % (X_MEAS, S_OPS, math.log2(S_OPS)))
    print("log2 phase1        %.8f" % log2f(PHASE1))
    print("log2 extras        %.8f  (table build 2^%.2f, meter cap 2^%.2f, restart 2^%.2f)"
          % (log2f(EXTRAS), log2f(TABLE_BUILD_REF), log2f(METER_CAP), log2f(RESTART_REF)))
    print("log2 W             %.8f   claim %.3f" % (LOG_W, CLAIM_LOG2))
    print("success DESIGN     %.6f = CP95 %.6f - debits %.3e - restart %.3e   floor met %s"
          % (SUCCESS, float(CP95_TOTAL), float(DEBITS), RESTART_SUCCESS_LOSS, FLOOR_MET))
    print("success LOW order  %.6f at 1/%.1f the class count -> floor met %s"
          % (float(CP95_LOW), ORDER_SPAN,
             float(CP95_LOW) - float(DEBITS) - RESTART_SUCCESS_LOSS >= float(SUCCESS_FLOOR)))
    print("w=2 comparison     pair census CP95 %.6f over %d sets (NOT the claim)"
          % (float(PAIR_CP), PAIR_SETS))
    print("absorption         p(ge3) %.6f -> %.6f at %.2fx classes"
          % (ABS["fixtures"][0]["p_ge3"], ABS["fixtures"][1]["p_ge3"],
             ABS["fixtures"][1]["classes"] / ABS["fixtures"][0]["classes"]))
    print("gate               %.6f  claim %.3f  delta %+.4f bits -> %s"
          % (float(GATE), float(CLAIM), float(MARGIN_VS_GATE),
             "RANKS" if CLAIM < GATE else "NO RANK (honest port)"))
    print("log2 memory        %.6f" % log2f(MEMORY_BYTES))
    print("records/reports    log2 records %.4f <= slots %.4f (nodes x cap)"
          % (math.log2(records), math.log2(REPORT_SLOTS)))
    print("SENSITIVITY        if one P-256 addition is charged as 2 reference units:"
          " claim %.3f (+%.3f bits), still no rank"
          % (CLAIM_2REF, CLAIM_2REF - CLAIM_LOG2))
    print("kappa table:")
    for n, v, why in KAPPA_ROWS:
        print("  %.9f  %-40s %s" % (float(v), n, why))
    print("  product %.9f <= 2^(claim-log2 S0) %.9f" % (float(BUDGET), CLAIM_RATIO_CEIL))
    out = {
        "test": "ledger_p256_v11", "C": C, "step_rows": ROWS, "step_insn": STEP_INSN,
        "step_price": float(STEP_PRICE), "N_classes_log2": math.log2(N_CLASSES),
        "log2_S0": LOG_S0, "shape_measured": X_MEAS, "S_OPS_log2": math.log2(S_OPS),
        "log2_phase1": log2f(PHASE1), "log2_extras": log2f(EXTRAS),
        "log2_W": LOG_W, "claim": CLAIM_LOG2, "gate": float(GATE),
        "margin_vs_gate": float(MARGIN_VS_GATE), "ranks": bool(CLAIM < GATE),
        "success": SUCCESS, "cp95_design": float(CP95_TOTAL),
        "cp95_low_order": float(CP95_LOW), "order_span": ORDER_SPAN,
        "floor": float(SUCCESS_FLOOR), "floor_met": bool(FLOOR_MET),
        "debits": float(DEBITS), "restart_loss": RESTART_SUCCESS_LOSS,
        "markov_loss_bound": MARKOV_LOSS_BOUND,
        "pair_cp95_w2": float(PAIR_CP), "pair_sets": PAIR_SETS,
        "memory_log2_bytes": log2f(MEMORY_BYTES),
        "fresh_insn": FRESH_INSN, "fresh_ref": float(FRESH_REF),
        "restart_ref_log2": log2f(RESTART_REF),
        "kappa_product": float(BUDGET), "claims_ratio_ceiling": CLAIM_RATIO_CEIL,
        "dp_bits": DESIGN_ARM_DP, "fixture_R": FIXTURE_R,
        "j_invariant": J_INV, "claim_at_2_reference_units": CLAIM_2REF,
        "claim_2ref_delta": CLAIM_2REF - CLAIM_LOG2,
        "report_slots_log2": math.log2(REPORT_SLOTS),
        "records_log2": math.log2(records),
        "fixture_trap_rate_log2": math.log2(1.0 / (2 * FIXTURE_R)),
        "trap_rate_log2": log2f(TRAP_RATE),
        "preprocessing_log2": log2f(TABLE_BUILD_REF),
    }
    with open(os.path.join(HERE, "evidence", "ledger-out.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    return out


if __name__ == "__main__":
    report()
