"""Exact accounting for the claim in proof.md (stdlib only; prints JSON).

The attack is run under a hard work meter. Every item of the attack is listed in ITEMS with its
H (reference units: batched incomplete affine G1 additions and doublings) and R (scalar AArch64
instructions, loads, stores and control included). The meter stops the attack when its total
reaches W_MAX, so time_log2 = log2(W_MAX) is an upper bound on the work by construction.
Measured per-item costs come from evidence/measure-v3/summary.json; analytic bounds are marked.
All integers and rationals are exact; the score is rounded up with exact.log2_up.
"""
import json
import math
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import exact  # noqa: E402  (vendored cscore/exact.py)

with open(os.path.join(HERE, "evidence", "measure-v3", "summary.json")) as fh:
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
RECORDS = 256                          # round-robin walk records (active trails)
DP_BITS = 90                           # distinguished with probability 2^-90
MU_DP = 1 << DP_BITS                   # mean distance between distinguished points
CAP = 1 << 100                         # trail cap, steps per trail
DETECT_K = 10                          # detection margin in units of RECORDS * MU_DP
DETECT_SPAN = RECORDS * MU_DP * (1 << DETECT_K)     # 2^108 points
S_MODEL = exact.iroot_ceil(-(-(4 * M * TAU.numerator) // TAU.denominator), 2)
B_MIN = S_MODEL + DETECT_SPAN          # points the meter guarantees to allow
RESTARTS_MAX = B_MIN // CAP            # restarts never exceed floor(points / CAP)
TRAILS = RECORDS + RESTARTS_MAX        # trails started in total

# Distinguished-point table (bounded): 2^20 buckets of 256 slots, abort on overflow.
N_BUCKETS = 1 << 20
BUCKET_SLOTS = 256
SLOT_BYTES = 72                        # 37-byte key prefix, 32-byte exponent, family, padding
DP_CAP = 1 << 27                       # abort above this many stored points (load 1/2)
DP_ACTION_BOUND = 16000                # per distinguished point: hash, scan of <= 256 slots, store

# Constants of the cost model and analytic bounds.
WINDOWS = 16
TBLN = 1 << 15                         # entries per comb window, both bases
ADDS_PER_POINT = WINDOWS               # reference units charged per walk point
CHAIN_CHECKS = 2 * K_EXP + 1           # comb checks in the square-root chain
SAMPLE = 4000                          # analytic bound per random element (PRNG and rejection)
FR_BOUND = 1000                        # analytic bound per mod-r Montgomery multiplication
TONELLI_MULS = 2 * 32 ** 2             # per square root in Fr, 2-adicity 32
SOLVE_BOUND = 1 << 18                  # inversion in Fr for the collision solve
CONST_BOUND = 1 << 16                  # parsing T and the powers used, beta and mu6
MULTIPLIERS = 64                       # walk multipliers, sampled once
RESTART_SAMPLING = True               # every restart and initial record uses a fresh random start

# Measured per-item costs (evidence/measure-v3/summary.json, differenced over 256 and 512 steps).
R_POINT = Fraction(MEAS["R_point_instructions"]).limit_denominator(1 << 20)
R_RESTART = Fraction(MEAS["restart_R_instructions_per_restart"]).limit_denominator(1 << 20)
NORM_PER_ENTRY_BOUND = 10000          # analytic: batch normalisation of one entry (6 mults, 1 squaring)
R_BUILD = 2 * WINDOWS * TBLN * NORM_PER_ENTRY_BOUND   # see evidence/measure-v3/summary.json (note_build)
H_BUILD_ADDS = MEAS["build_proj_calls"]                 # measured table additions (exact)
BUILD_DOUBLINGS = 2 * WINDOWS * 16                      # base shifts, one unit each
INIT_RECORDS = RECORDS

# Bucket memory: a header of 4 bytes and 256 slots of 72 bytes each.
BUCKET_BYTES = 4 + BUCKET_SLOTS * SLOT_BYTES
DP_TABLE_BYTES = N_BUCKETS * BUCKET_BYTES
TABLE_BYTES = 2 * WINDOWS * TBLN * 96  # two bases, affine points of 2 x 48 bytes
RECORD_BYTES = RECORDS * 256
CODE_BYTES = 1 << 24

ITEMS = []


def item(name, phase, H=0, R=0):
    ITEMS.append((name, phase, Fraction(H), Fraction(R)))


# Preprocessing.
item("comb tables: projective table additions (measured count)", "preprocessing", H=H_BUILD_ADDS)
item("comb tables: base shifts (16 doublings per window and base)", "preprocessing", H=BUILD_DOUBLINGS)
item("comb tables: batch normalisation of every entry (analytic bound)", "preprocessing", R=R_BUILD)
item("walk records: initial comb evaluations (bound: one step per record)", "preprocessing",
     H=INIT_RECORDS * ADDS_PER_POINT, R=INIT_RECORDS * R_POINT)
item("multipliers: sampling and 20 squarings (bound)", "preprocessing",
     R=MULTIPLIERS * (SAMPLE + K_EXP * FR_BOUND))
item("initial records: fresh start exponent, sampling and 20 squarings (bound)", "preprocessing",
     R=RECORDS * (SAMPLE + K_EXP * FR_BOUND))
item("DP table zeroing (one store per 8-byte word)", "preprocessing",
     R=N_BUCKETS * (4 + BUCKET_SLOTS * SLOT_BYTES) // 8)
item("constants and input parse (bound)", "preprocessing", R=CONST_BOUND)

# Walk: every point, restarts and the worst-case distinguished-point work.
item("walk points: 16 units each", "walk", H=ADDS_PER_POINT * B_MIN)
item("walk points: measured other instructions", "walk", R=R_POINT * B_MIN)
item("restarts: 16 units each (at most floor(points / CAP))", "walk",
     H=ADDS_PER_POINT * RESTARTS_MAX)
item("restarts: measured other instructions", "walk", R=R_RESTART * RESTARTS_MAX)
item("restarts: fresh start exponent, sampling and 20 squarings (bound)", "walk",
     R=RESTARTS_MAX * (SAMPLE + K_EXP * FR_BOUND))
item("distinguished points: worst case, DP_CAP x bound", "walk", R=DP_CAP * DP_ACTION_BOUND)

# Finish.
item("chain: comb checks (41)", "finish", H=ADDS_PER_POINT * CHAIN_CHECKS)
item("chain: one walk point's other work per check", "finish", R=CHAIN_CHECKS * R_POINT)
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

# Success losses (each an upper bound; proof.md section 3).
LOSS_CAP = 2 * (MU_DP / CAP)                       # two trails, cap truncation of a merged path
LOSS_DETECT = 2 * math.exp(-DETECT_SPAN / (RECORDS * MU_DP))
LOSS_NODP = TRAILS * math.exp(-CAP / MU_DP)        # a whole trail without any distinguished point
LOSS_CYCLE = float(Fraction(TRAILS) * Fraction(CAP * CAP, 2 * M))   # in-trail cycles, birthday bound
# Bucket overflow: n_b x P(Poisson(lambda) >= BUCKET_SLOTS + 1) with lambda = DP_CAP / n_b.
lam = DP_CAP / N_BUCKETS
k = BUCKET_SLOTS + 1
log_tail = -lam + k * (1 + math.log(lam) - math.log(k))   # Chernoff: e^-lam (e lam / k)^k
LOSS_OVERFLOW = N_BUCKETS * math.exp(log_tail)
LOSS_TOTAL = LOSS_CAP + LOSS_DETECT + LOSS_NODP + LOSS_CYCLE + LOSS_OVERFLOW

# Statistics: the fresh toy runs give the Wilson 95% lower bound 0.4227 at the claim budget.
STAT_LOWER = Fraction(4227, 10000)
CLAIM_SUCCESS = Fraction(40, 100)
assert float(STAT_LOWER) * (1 - LOSS_TOTAL) >= float(CLAIM_SUCCESS)

out = {
    "tau_model": str(TAU),
    "S_model_log2": float(exact.log2(S_MODEL)),
    "detect_span_log2": float(exact.log2(DETECT_SPAN)),
    "B_min_log2": float(exact.log2(B_MIN)),
    "restarts_max": RESTARTS_MAX,
    "trails": TRAILS,
    "dp_cap_log2": float(exact.log2(DP_CAP)),
    "dp_expected_log2": float(exact.log2(B_MIN / 2 ** DP_BITS)),
    "loss_cap": LOSS_CAP,
    "loss_detect": LOSS_DETECT,
    "loss_no_dp_trail": LOSS_NODP,
    "loss_cycle": LOSS_CYCLE,
    "loss_overflow": LOSS_OVERFLOW,
    "loss_total": LOSS_TOTAL,
    "success_model": 1 - math.exp(-float(TAU)),
    "success_stat_lower": str(STAT_LOWER),
    "success_claim": str(CLAIM_SUCCESS),
    "H_total": str(H_TOTAL),
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
