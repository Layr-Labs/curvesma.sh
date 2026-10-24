"""ledger_v13.py -- p256-ecdlp v13: FULL-STORE detection, no report gate.

THE ONE STRUCTURAL CHANGE.  v12 stored a class only when the DP gate fired, so
measured success depended on the gate density (ladder slope -1.326 bits of odds
per bit of density) and the production density 2^-60 was an extrapolation the
reviewer refused.  v13 deletes the gate from the success path: every visited
class is inserted and looked up, and a repeat is seen AT INSERT.  Success is
then a function of x = visits/sqrt(N) only.  The traces are unchanged (shipped
qwalk.step, same escape) -- only the detector differs, which is why the same
census reads both ways (evidence/census-design.json: gate-off d0 row = 0.685
shot at 3200 sets; evidence/fullstore-census.json re-reads its own traces
through both detectors as a matched ablation).

PRICE.  Insert+lookup is charged per visited step, itemised in ROWS_INSERT.
The DP-fingerprint row of v12 is retained (over-charge, harmless).  The step
price is 1 + STEP_INSN/C with C = 1000, the track's protected normalization.

MEMORY IS THE PRICE OF DELETING THE GATE, AND IT IS STATED, NOT HIDDEN:
records = S_OPS (one per visit), entry 40 B, capacity 2x for load factor 0.5.
The bound is a TIME bound; rules.md makes memory supplementary.  The
conservative alternative charging one full reference unit per random table
access is emitted as a SENSITIVITY row, not claimed.
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
C = 1000
R_TAB = 2 ** 20
R_EFF = selector.reachable_columns_span(R_TAB, p)
assert R_EFF == R_TAB
N_CLASSES = (r - 1) // 2
N_CELLS = 2 * R_EFF
workers = 2 ** 20
A_H1 = F(1015, 1000)
MODEL_SUCCESS = F(3903, 10000)
GATE = F(127237678, 1000000)
SUCCESS_FLOOR = F(39, 100)
ENTRY_BYTES = 40
LOAD_DEN = 2                  # capacity = 2 x arrivals (load factor 0.5)

# ---- step rows ------------------------------------------------------------
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
# v13 insert rows: the FULL-STORE detector charged on every visited step
ROWS_INSERT = [
    ("64-bit class mix for the bucket index (4)", 4),
    ("bucket address: base + index*40 wide multiply-add (3)", 3),
    ("bucket load of the 32-byte stored class and compare (5)", 5),
    ("probe advance on a miss, average 1.5 probes (6)", 6),
    ("occupancy store: 32-byte class + 8-byte (trail, step) tag (7)", 7),
    ("table load-factor counter add and threshold reload (2)", 2),
]
STEP_INSN = sum(n for _o, n in ROWS)
INSERT_INSN = sum(n for _o, n in ROWS_INSERT)
assert STEP_INSN == 263, STEP_INSN
assert INSERT_INSN == 27, INSERT_INSN
TOTAL_INSN = STEP_INSN + INSERT_INSN        # spec constant, three-way asserted
STEP_PRICE = 1 + F(TOTAL_INSN, C)
STEP_PRICE_1REF = 2 + F(TOTAL_INSN, C)      # sensitivity: 1 ref unit per table access

# ---- S0 / executed budget -------------------------------------------------
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
SHAPE_BIRTHDAY = math.sqrt(2 * float(frac_of_dec(_ln0)))


def load(rel):
    with open(os.path.join(HERE, rel)) as f:
        return json.load(f)


DESIGN = load("evidence/census-design.json")
FS = load("evidence/fullstore-census.json")

ARM = DESIGN["arms"]["escaped_lag3"]
N_TRIALS = int(ARM["trials"])
# FULL-STORE success on the 3200-set shipped census: the d0 (gate-off) row is
# exactly insert-time detection -- every coalescence is seen, none is gated out.
FS_SUCCESSES = int(round(ARM["d0_implied_success"] * N_TRIALS))
assert FS_SUCCESSES == 2192, FS_SUCCESSES


ALPHA = Decimal("0.05")


def cp95_lower(k, n):
    """Exact one-sided Clopper-Pearson lower bound, bracket (0,1) (no clip)."""
    k, n = int(k), int(n)
    if k <= 0:
        return 0.0
    if k >= n:
        return float(ALPHA ** (Decimal(1) / Decimal(n)))
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p_):
        p_ = Decimal(p_); q_ = Decimal(1) - p_
        t = (lnC + k * p_.ln() + (n - k) * q_.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p_ / ((j + 1) * q_)
            s += t
            if t < Decimal("1e-40") * s:
                break
        return s

    lo, hi = Decimal("0.000001"), Decimal("0.999999")
    for _ in range(300):
        mid = (lo + hi) / 2
        if tail(mid) > ALPHA:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-14"):
            break
    return float(lo)


CP95_FS = cp95_lower(FS_SUCCESSES, N_TRIALS)
# corroborating full-store cell re-run at the same charged budget (v13 census)
FS_CELL = [c for c in FS["cells"] if abs(c["x_visits_over_sqrtN"] - 1.6002022719673725) < 5e-3][0]
CP95_FS_CELL = float(FS_CELL["fullstore"]["cp95_success_f"])
assert CP95_FS_CELL > 0.0
X_MEAS = float(DESIGN["visits_over_sqrtN"])
SQRT_NP = math.sqrt(float(Nprime))
S_OPS = int(math.ceil(float(A_H1) * X_MEAS * SQRT_NP))
LOG_S0 = math.log2(S0.numerator) - math.log2(S0.denominator)

# ---- cap / support row (g1-ecdlp accepted form) ---------------------------
# support = records a trail contributes = its step count, a DESIGN CONSTANT of
# this detector (a full-store trail has no DP truncation), so the cap row is a
# deterministic function of the cap and carries no distributional premise:
K_MAX = int(DESIGN["per_walker_steps"])          # = 4632 at the fixture
SUPPORT_MEASURED_MAX = int(DESIGN["arms"]["escaped_lag3"]["solve_step_max"])
assert SUPPORT_MEASURED_MAX <= K_MAX, (SUPPORT_MEASURED_MAX, K_MAX)
SUPPORT_PER_TRAIL_PROD = S_OPS // workers
assert SUPPORT_PER_TRAIL_PROD * workers <= S_OPS

# ---- table / memory -------------------------------------------------------
RECORDS = S_OPS
BUCKETS = LOAD_DEN * RECORDS
TABLE_INIT_INSN = 4 * BUCKETS                   # 4 insn per bucket, one time
TABLE_INIT_REF = F(TABLE_INIT_INSN, C)
TABLE_BYTES = BUCKETS * ENTRY_BYTES

# ---- other extras (unchanged from v12, re-derived here) -------------------
nodes = 2 ** 16
M_insn = 97
node_setup = F(nodes) * F(2 ** 20)
meter_block_remainders = F(2 ** 61)
verification_output = F(2 ** 40)
fresh_starts = F(workers) * F(16)
BATCH = 2 ** 18
INV_INSN = (254 + 11) * M_insn
n_batches = ((S_OPS + 2 ** 60) + BATCH - 1) // BATCH
batch_inversion_total = F(INV_INSN * (n_batches + 1), C)
TRAP_RATE = F(1, 2 * R_EFF)
E_RUN_MEAS = 2
TRAP_REF = F(S_OPS) * TRAP_RATE * E_RUN_MEAS * STEP_PRICE
METER_CAP = 2 * TRAP_REF
CYC3_WORK = F(S_OPS) * F(1, R_EFF * R_EFF) * 3 * STEP_PRICE
FALLBACK_RATE = 3 * F(1, 2 ** 64)
FALLBACK_REF = F(S_OPS) * FALLBACK_RATE * M_insn / C
EXTRAS = (TABLE_INIT_REF + node_setup + meter_block_remainders + verification_output +
          fresh_starts + batch_inversion_total + METER_CAP + CYC3_WORK + FALLBACK_REF)
PHASE1 = STEP_PRICE * (S_OPS + 2 ** 60)
W = PHASE1 + EXTRAS

# ---- success --------------------------------------------------------------
relation_failure = F(1, r)
pool_overflow = F(336, 1000000)
table_abort = F(1, 2 ** 40)          # load-factor abort, one-sided bound
DEBITS = pool_overflow + relation_failure + table_abort
RESTART_SUCCESS_LOSS = 0.0           # no restarts: escape kcap never exceeded (measured 0)
SUCCESS = CP95_FS - float(DEBITS) - RESTART_SUCCESS_LOSS
SLACK = SUCCESS - float(SUCCESS_FLOOR)


def log2f(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


LOG_W = log2f(W)
CLAIM_LOG2 = math.ceil(LOG_W * 1000) / 1000.0
CLAIM = F(int(round(CLAIM_LOG2 * 1000)), 1000)
CLAIM_1REF = math.ceil((LOG_W + math.log2(float(STEP_PRICE_1REF / STEP_PRICE))) * 1000) / 1000.0

# integer-sum audit: the charged rows must be the sum of the itemised rows
assert sum(n for _o, n in ROWS) + sum(n for _o, n in ROWS_INSERT) == TOTAL_INSN
assert TOTAL_INSN == 290

out = {
    "test": "ledger_p256_v13_fullstore",
    "C": C, "step_rows": ROWS, "insert_rows": ROWS_INSERT,
    "step_insn": STEP_INSN, "insert_insn": INSERT_INSN, "total_insn": TOTAL_INSN,
    "step_price": float(STEP_PRICE), "step_price_1ref_access": float(STEP_PRICE_1REF),
    "insn_row_sum": sum(n for _o, n in ROWS) + sum(n for _o, n in ROWS_INSERT),
    "N_classes_log2": math.log2(N_CLASSES), "log2_S0": LOG_S0,
    "shape_measured": X_MEAS, "S_OPS_log2": math.log2(S_OPS),
    "log2_phase1": log2f(PHASE1), "log2_extras": log2f(EXTRAS),
    "log2_W": LOG_W, "claim": CLAIM_LOG2, "claim_at_1ref_access": CLAIM_1REF,
    "claim_delta_1ref": CLAIM_1REF - CLAIM_LOG2,
    "gate": float(GATE), "margin_vs_gate": float(GATE - CLAIM),
    "ranks": bool(CLAIM < GATE),
    "success": SUCCESS, "cp95_fullstore": CP95_FS,
    "fullstore_successes": FS_SUCCESSES, "fullstore_trials": N_TRIALS,
    "fullstore_rate": FS_SUCCESSES / N_TRIALS,
    "debits": float(DEBITS), "floor": float(SUCCESS_FLOOR),
    "floor_met": bool(SUCCESS >= float(SUCCESS_FLOOR)), "floor_slack": SLACK,
    "dp_gate_on_success_path": False,
    "records_log2": math.log2(RECORDS), "entry_bytes": ENTRY_BYTES,
    "buckets_log2": math.log2(BUCKETS), "table_bytes_log2": math.log2(TABLE_BYTES),
    "table_bytes": float(TABLE_BYTES), "table_init_insn": TABLE_INIT_INSN,
    "table_init_ref_log2": log2f(TABLE_INIT_REF),
    "memory_bytes_log2": math.log2(float(TABLE_BYTES) + float(nodes * 2 ** 20 * 192)),
    "k_max_per_trail": K_MAX, "support_measured_max": SUPPORT_MEASURED_MAX,
    "support_per_trail_production": SUPPORT_PER_TRAIL_PROD,
    "support_is_design_constant": True, "p_abort_support": 0.0,
    "abort_bound_one_sided": 3.0 / N_TRIALS,
    "restart_success_loss": RESTART_SUCCESS_LOSS,
    "cp95_fullstore_cell_1600": CP95_FS_CELL,
    "fullstore_cell_rate": FS_CELL["fullstore"]["success_rate"],
    "gate_ablation": {k: {"drop_vs_fullstore": v["missing_vs_fullstore"] / v["trials"],
                          "success_rate": v["success_rate"]}
                      for k, v in FS_CELL["gate_replay"].items()},
    "gate_ablation_cell_x": FS_CELL["x_visits_over_sqrtN"],
    "gate_ablation_cell_trials": FS_CELL["trials"],
    "fullstore_census_cells": [
        {"classes_log2": c["classes_log2"], "x": c["x_visits_over_sqrtN"],
         "trials": c["trials"], "success_rate": c["fullstore"]["success_rate"],
         "cp95": c["fullstore"]["cp95_success_f"],
         "model_birthday": c["fullstore"]["model_birthday"]} for c in FS["cells"]],
    "trap_rate_log2": log2f(TRAP_RATE),
    "preprocessing_log2": log2f(TABLE_INIT_REF),
    "worker_count_log2": math.log2(workers),
    "fresh_ref": float(F(492300, C)),
}
with open(os.path.join(HERE, "evidence", "ledger-out.json"), "w") as f:
    json.dump(out, f, indent=1, sort_keys=True)

if __name__ == "__main__":
    print("P-256 v13 full-store ledger")
    print("  step rows %d + insert rows %d = %d insn -> price %.6f"
          % (STEP_INSN, INSERT_INSN, TOTAL_INSN, float(STEP_PRICE)))
    print("  log2 S0 %.8f  x %.6f  log2 S_OPS %.8f" % (LOG_S0, X_MEAS, math.log2(S_OPS)))
    print("  log2 W %.8f  claim %.3f  (1-ref-access sensitivity %.3f)"
          % (LOG_W, CLAIM_LOG2, CLAIM_1REF))
    print("  success %.6f = CP95 %.6f - debits %.3e  floor met %s (slack %.4f)"
          % (SUCCESS, CP95_FS, float(DEBITS), SUCCESS >= float(SUCCESS_FLOOR), SLACK))
    print("  records 2^%.4f x %d B, buckets 2^%.4f, table 2^%.2f B"
          % (math.log2(RECORDS), ENTRY_BYTES, math.log2(BUCKETS), math.log2(TABLE_BYTES)))
    print("  K_max per trail %d (design constant), measured max %d"
          % (K_MAX, SUPPORT_MEASURED_MAX))
    print("  gate ablation at x=%.3f: %s"
          % (out["gate_ablation_cell_x"], json.dumps(out["gate_ablation"])))
    print("  gate %.6f  delta %+.4f bits" % (float(GATE), float(GATE - CLAIM)))
