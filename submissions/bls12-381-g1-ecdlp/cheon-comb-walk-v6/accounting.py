"""Exact accounting for the claim in proof.md (stdlib only; prints JSON).

The attack runs under a hard work meter. Every item is listed in ITEMS with its H (reference
units: batched incomplete affine G1 additions, one unit each, the 6M convention) and R (other
scalar instructions, loads, stores and control included). The meter halts the attack when its
total reaches W_MAX, so time_log2 = log2(W_MAX) bounds the work by construction.

Measured per-item costs come from evidence/measure-v4/summary.json. Charges derived from them are
worst-case: each measured mean is increased by the bound on its data-dependent branches, stated
in the comments. Analytic bounds are marked. All integers and rationals are exact; the score is
rounded up with exact.log2_up.
"""
import json
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exact  # noqa: E402  (vendored cscore/exact.py)

with open(os.path.join(HERE, "evidence", "measure-v6exact", "summary.json")) as fh:
    MEAS = json.load(fh)

r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 2000
K_EXP = 20
d = 1 << K_EXP
N = (r - 1) // d
assert (r - 1) % d == 0 and N % 6 == 0
M = N // 6
assert r.bit_length() == 255

# Walk parameters.
TAU = Fraction(13, 20)                 # model success 1 - exp(-TAU)
WBITS = 20                             # comb window width
WINDOWS = 13                           # ceil(256 / 20) signed windows
ADDS_PER_POINT = WINDOWS - 1           # exact worst case: the first nonzero window is a copy
RECORDS = 4096                         # round-robin records (active trails)
DP_BITS = 90
MU_DP = 1 << DP_BITS
CAP = 1 << 100
DETECT_K = 10
DETECT_SPAN = RECORDS * MU_DP * (1 << DETECT_K)          # 2^112 points
S_MODEL = exact.iroot_ceil(-(-(4 * M * TAU.numerator) // TAU.denominator), 2)
B_MIN = S_MODEL + DETECT_SPAN
RESTARTS_MAX = B_MIN // CAP
TRAILS = RECORDS + RESTARTS_MAX

# Bounded distinguished-point table.
N_BUCKETS = 1 << 20
BUCKET_SLOTS = 256
SLOT_BYTES = 72
DP_CAP = 1 << 27
DP_ACTION_BOUND = 16000

# Analytic bounds.
TBLN = 1 << 19                         # entries per comb window and base (magnitudes 1..2^19)
SAMPLE = 4000                          # per random element (PRNG and rejection)
FR_BOUND = 1000                        # per mod-r Montgomery multiplication
CHAIN_CHECKS = 2 * K_EXP + 1
TONELLI_MULS = 2 * 32 ** 2
SOLVE_BOUND = 1 << 18
CONST_BOUND = 1 << 16
MULTIPLIERS = 64
RESTART_CHARGE = Fraction(MEAS["restart_other_per_restart_exact"])   # forced path: 12 unbatched additions (no reference core)

# Measured costs (summary.json; one round = 4096 points, differenced over two rounds).
R_POINT_MEAS = Fraction(MEAS["natural_other_per_point_Ir"]).limit_denominator(1 << 20)   # natural mean, reported only (float)
R_POINT_FORCED = Fraction(MEAS["forced_other_per_point_exact"])   # worst case, exact (reference core by exact call counts)
BUILD_STAGES = (1 << 19) - 2                              # table stages j = 2 .. 2^19 - 1
R_BUILD_MEAS = Fraction(MEAS["build_fixed_other_exact"]) + BUILD_STAGES * Fraction(MEAS["build_other_per_stage_exact"])
BUILD_ADDS = int(MEAS["build_adds_per_stage"]) * BUILD_STAGES   # 26 chains (2 bases x 13 windows) per stage, exact count
INIT_ADDS = RECORDS * ADDS_PER_POINT   # initial comb evaluations (worst case)

# Worst-case per-point other work: the maximal path, measured with every digit nonzero and all 15 stages
# occupied (forced run). Every other data-dependent path executes less work (proof.md section 5).
R_POINT_CHARGE = R_POINT_FORCED

BUCKET_BYTES = 4 + BUCKET_SLOTS * SLOT_BYTES
DP_TABLE_BYTES = N_BUCKETS * BUCKET_BYTES
TABLE_BYTES = 2 * WINDOWS * TBLN * 96   # two bases, affine points of 2 x 48 bytes (1.31 GB)
RECORD_BYTES = RECORDS * 256
CODE_BYTES = 1 << 24

ITEMS = []


def item(name, phase, H=0, R=0):
    ITEMS.append((name, phase, Fraction(H), Fraction(R)))


# Preprocessing.
item("comb tables: batched affine additions (reference, 26 chains per stage)", "preprocessing", H=BUILD_ADDS)
item("comb tables: other work, measured per stage and fixed part (inversions, doublings, bookkeeping)", "preprocessing",
     R=R_BUILD_MEAS)
item("initial records: batched comb evaluations (12 additions each, worst case)", "preprocessing", H=INIT_ADDS)
item("initial records: other work, bound (one point's worst-case other work each)", "preprocessing",
     R=RECORDS * R_POINT_CHARGE)
item("initial records: fresh start exponents, sampling and 20 squarings (bound)", "preprocessing",
     R=RECORDS * (SAMPLE + K_EXP * FR_BOUND))
item("multipliers: sampling and 20 squarings (bound)", "preprocessing",
     R=MULTIPLIERS * (SAMPLE + K_EXP * FR_BOUND))
item("DP table zeroing (one store per 8-byte word)", "preprocessing",
     R=N_BUCKETS * (4 + BUCKET_SLOTS * SLOT_BYTES) // 8)
item("constants and input parse (bound)", "preprocessing", R=CONST_BOUND)

# Walk: points, restarts and worst-case distinguished-point work.
item("walk points: 12 batched additions each, B_min points", "walk", H=ADDS_PER_POINT * B_MIN)
item("walk points: other work, worst case per point", "walk", R=R_POINT_CHARGE * B_MIN)
item("restarts: fresh start exponent, sampling and 20 squarings (bound), at most floor(B/CAP)", "walk",
     R=RESTARTS_MAX * (SAMPLE + K_EXP * FR_BOUND))
item("restarts: unbatched comb evaluation, forced worst case per restart (12 additions)", "walk",
     R=RESTARTS_MAX * RESTART_CHARGE)
item("distinguished points: worst case, DP_CAP x bound", "walk", R=DP_CAP * DP_ACTION_BOUND)

# Finish: the chain's comb checks are unbatched (single exponents), so they are other work.
item("chain: 41 unbatched comb evaluations (bound)", "finish", R=CHAIN_CHECKS * RESTART_CHARGE)
item("chain: 20 square roots in Fr (bound)", "finish", R=K_EXP * TONELLI_MULS * FR_BOUND)
item("collision solve (bound)", "finish", R=SOLVE_BOUND)

H_TOTAL = sum(H for _, _, H, _ in ITEMS)
R_TOTAL = sum(R for _, _, _, R in ITEMS)
W_MAX = H_TOTAL + R_TOTAL / C
time_log2 = exact.log2_up(W_MAX)

pre_H = sum(H for _, ph, H, _ in ITEMS if ph == "preprocessing")
pre_R = sum(R for _, ph, _, R in ITEMS if ph == "preprocessing")
preprocessing_log2 = exact.log2_up(pre_H + pre_R / C)

memory_bytes = DP_TABLE_BYTES + TABLE_BYTES + RECORD_BYTES + CODE_BYTES
memory_log2 = exact.log2_up(Fraction(memory_bytes))

# Success losses (upper bounds; proof.md section 4).
LOSS_CAP = 2 * (MU_DP / CAP)
LOSS_DETECT = 2 * math.exp(-DETECT_SPAN / (RECORDS * MU_DP))
LOSS_NODP = TRAILS * math.exp(-CAP / MU_DP)
LOSS_CYCLE = float(Fraction(TRAILS) * Fraction(CAP * CAP, 2 * M))
lam = DP_CAP / N_BUCKETS
k = BUCKET_SLOTS + 1
log_tail = -lam + k * (1 + math.log(lam) - math.log(k))
LOSS_OVERFLOW = N_BUCKETS * math.exp(log_tail)
LOSS_LEMMA_LOG2 = -100.0                # window-12 coincidence: <= 2^117 points x 2^21 / 2^240 (proof.md section 2)
LOSS_TOTAL = LOSS_CAP + LOSS_DETECT + LOSS_NODP + LOSS_CYCLE + LOSS_OVERFLOW + 2.0 ** LOSS_LEMMA_LOG2

STAT_LOWER = Fraction(4227, 10000)
CLAIM_SUCCESS = Fraction(40, 100)
assert float(STAT_LOWER) * (1 - LOSS_TOTAL) >= float(CLAIM_SUCCESS)

out = {
    "tau_model": str(TAU),
    "adds_per_point": ADDS_PER_POINT,
    "R_point_natural_mean": float(R_POINT_MEAS),
    "R_point_charge_forced_max": float(R_POINT_CHARGE),
    "S_model_log2": float(exact.log2(S_MODEL)),
    "detect_span_log2": float(exact.log2(DETECT_SPAN)),
    "B_min_log2": float(exact.log2(B_MIN)),
    "restarts_max": RESTARTS_MAX,
    "trails": TRAILS,
    "dp_expected_log2": float(exact.log2(B_MIN / 2 ** DP_BITS)),
    "loss_cap": LOSS_CAP,
    "loss_detect": LOSS_DETECT,
    "loss_no_dp_trail": LOSS_NODP,
    "loss_cycle": LOSS_CYCLE,
    "loss_overflow": LOSS_OVERFLOW,
    "loss_lemma_log2": LOSS_LEMMA_LOG2,
    "loss_total": LOSS_TOTAL,
    "success_model": 1 - math.exp(-float(TAU)),
    "success_stat_lower": str(STAT_LOWER),
    "success_claim": str(CLAIM_SUCCESS),
    "H_total": str(H_TOTAL),
    "H_log2": float(exact.log2(H_TOTAL)),
    "R_total_log2": float(exact.log2(R_TOTAL)),
    "R_over_C_log2": float(exact.log2(R_TOTAL / C)),
    "W_max_log2": float(exact.log2(W_MAX)),
    "time_log2": str(time_log2),
    "time_log2_float": float(exact.log2(W_MAX)),
    "score_rational": {"numerator": str(time_log2.numerator), "denominator": str(time_log2.denominator)},
    "preprocessing_log2": str(preprocessing_log2),
    "memory_log2_bytes": str(memory_log2),
    "memory_bytes": memory_bytes,
    "advice_bytes": 0,
    "items": [{"name": n, "phase": ph, "H_units": str(H), "R_instructions": str(R)}
              for n, ph, H, R in ITEMS],
}

if __name__ == "__main__":
    print(json.dumps(out, indent=2))
