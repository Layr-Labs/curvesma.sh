"""Exact accounting for the claim in proof.md (stdlib only; prints JSON).

W = H + R/C for n = 2^20. H counts reference operations (batched incomplete affine G1
additions and doublings, one unit each). R counts other scalar AArch64 instructions, including
loads, stores and control. C = 2000 is the BLS12-381 G1 price in tracks/docs/COST_MODEL.md.
Every item of the attack appears once in ITEMS. All integers and rationals are exact; the
score is rounded up at six decimals by exact.log2_up.

Measured inputs (evidence/measure-full-step/summary.json):
  R_STEP   instructions per walk step, the full step as run, both families
  REF_STEP instructions per step inside the reference functions ref_add, ref_fp_mul,
           ref_fr_mul; each call of ref_fr_mul is at most REF_STEP / 2 (two calls per step)
Analytic bounds (proof.md section 5): SAMPLE, DP_ACTION, the chain and the solve.
"""
import json
import math
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exact  # noqa: E402  (vendored cscore/exact.py)

r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 2000
K_EXP = 20
d = 1 << K_EXP
N = (r - 1) // d                  # |H|
assert (r - 1) % d == 0 and N % 6 == 0
M = N // 6                        # classes under mu6
assert r.bit_length() == 255

# Model and stopping rule (proof.md sections 3 and 4).
TAU = Fraction(13, 20)            # model success 1 - exp(-TAU) = 0.4780; fresh toy data (proof.md section 6)
CLAIM_SUCCESS = Fraction(40, 100) # below the Wilson 95% lower bound 0.4227 of the fresh toy data at TAU
DP_BITS = 90                      # distinguished with probability 2^-90
CAP = 1 << 100                    # trail cap, steps
DETECT_SPAN = 1 << 100            # the walk stops at S_model + DETECT_SPAN points in total
S_MODEL = exact.iroot_ceil(-(-(4 * M * TAU.numerator) // TAU.denominator), 2)
BUDGET = S_MODEL + DETECT_SPAN    # total walk points, both families
DP_ENTRIES = -(-BUDGET // 2 ** DP_BITS)

# Walk step (measured, evidence/measure-full-step).
R_STEP = Fraction(713344, 256)    # (Ir(512) - Ir(256)) / 256, full step, 256 walk records
REF_STEP = 925                    # per-step instructions inside ref_add, ref_fp_mul, ref_fr_mul
FR_CALL = REF_STEP // 2           # bound per ref_fr_mul call (two per step)
WINDOWS = 16
ADDS_PER_STEP = WINDOWS           # one comb evaluation per walk point

# Analytic bounds for work that the step measurement does not execute.
SAMPLE = 4000                     # per random element: PRNG output and rejection test
DP_ACTION = 4000                  # per distinguished point: hash, probe, store, compare
SOLVE_BOUND = 2 ** 18             # collision inversion in Fr
PARSE_BOUND = 2 ** 16             # reading T and the 22 powers used
TONELLI_V = 32                    # 2-adicity of r - 1
TONELLI_MULS = 2 * TONELLI_V ** 2 # bound on modular multiplications per square root

# Comb tables (both bases, 16 windows of 2^15 entries, affine points of 192 bytes).
TABLE_ENTRIES = 2 * WINDOWS * 2 ** 15
TABLE_BYTES = TABLE_ENTRIES * 192
TABLE_DOUBLINGS = 2 * (WINDOWS - 1) * 16    # base shifts by 2^16 per window
TABLE_WORDS_MOVED = 2 * 24 * TABLE_ENTRIES  # one load and one store per 8-byte word

# Distinguished-point record: exact high key bits (381-bit canonical x minus the 90 zero bits,
# 291 bits, 37 bytes), exponent (32 bytes), family (1 byte), padded to 72 bytes.
RECORD_BYTES = 72
SLOTS = 2 * DP_ENTRIES                      # load factor 2
ZERO_WORDS = SLOTS * RECORD_BYTES // 8      # one store per word to clear the table

TRAILS = -(-BUDGET // CAP)                  # trail starts over the whole run
MULTIPLIERS = 64
CODE_BYTES = 2 ** 24                        # code, stacks, scratch

# Items: (name, phase, H units, R instructions).
ITEMS = []


def item(name, phase, H=0, R=0):
    ITEMS.append((name, phase, H, Fraction(R)))


# Preprocessing.
item("comb tables: batched additions (2^19 per base)", "preprocessing", H=TABLE_ENTRIES)
item("comb tables: base shifts, doublings", "preprocessing", H=TABLE_DOUBLINGS)
item("comb tables: loads and stores of entries", "preprocessing", R=TABLE_WORDS_MOVED)
item("multipliers: sampling and x^d, 64 elements", "preprocessing",
     R=MULTIPLIERS * (SAMPLE + K_EXP * FR_CALL))
item("DP table zeroing", "preprocessing", R=ZERO_WORDS)
item("input parse (T and powers used)", "preprocessing", R=PARSE_BOUND)

# Walk.
item("walk points: one comb evaluation each (16 additions)", "walk", H=ADDS_PER_STEP * BUDGET)
item("walk points: measured step, instructions", "walk", R=R_STEP * BUDGET)
item("distinguished points: insert and lookup, bound", "walk", R=DP_ENTRIES * DP_ACTION)

# Trails: each start samples c_0 in H (20 squarings), then one comb evaluation.
item("trail starts: comb evaluation (16 additions)", "walk", H=ADDS_PER_STEP * TRAILS)
item("trail starts: sampling and 20 squarings", "walk",
     R=TRAILS * (SAMPLE + K_EXP * FR_CALL))
item("trail starts: step-sized overhead, bound", "walk", R=TRAILS * R_STEP)

# Finish: square-root chain, verification, solve, output.
CHAIN_CHECKS = 2 * K_EXP + 1
item("chain: comb checks (41)", "finish", H=ADDS_PER_STEP * CHAIN_CHECKS)
item("chain: step-sized overhead per check", "finish", R=CHAIN_CHECKS * R_STEP)
item("chain: 20 square roots in Fr (Tonelli-Shanks bound)", "finish",
     R=K_EXP * TONELLI_MULS * FR_CALL)
item("collision solve: inversion bound", "finish", R=SOLVE_BOUND)

H_TOTAL = sum(H for _, _, H, _ in ITEMS)
R_TOTAL = sum(R for _, _, _, R in ITEMS)
W = Fraction(H_TOTAL) + R_TOTAL / C
time_log2 = exact.log2_up(W)

pre_H = sum(H for _, ph, H, _ in ITEMS if ph == "preprocessing")
pre_R = sum(R for _, ph, _, R in ITEMS if ph == "preprocessing")
preprocessing_log2 = exact.log2_up(Fraction(pre_H) + pre_R / C)

peak_bytes = (SLOTS * RECORD_BYTES            # DP table at load factor 2
              + TABLE_BYTES                   # comb tables, both bases
              + TRAILS * 64                   # trail state
              + CODE_BYTES)
memory_log2 = exact.log2_up(Fraction(peak_bytes))

# Stopping rule: a collision at point t is used if its merged tail is detected by the budget.
# The tail costs about 2G points (both trails walk the merged path to its next DP, G geometric
# with mean 2^DP_BITS), so the loss is at most exp(-DETECT_SPAN / (2 * 2^DP_BITS)).
STOP_LOSS = math.exp(-DETECT_SPAN / (2 * 2 ** DP_BITS))
SUCCESS_MODEL_FLOAT = (1 - math.exp(-float(TAU))) * (1 - STOP_LOSS)

out = {
    "r_bits": r.bit_length(),
    "N_log2": float(exact.log2(N)),
    "M_log2": float(exact.log2(M)),
    "tau_model": str(TAU),
    "S_model_log2": float(exact.log2(S_MODEL)),
    "budget_log2": float(exact.log2(BUDGET)),
    "detect_span_log2": float(exact.log2(DETECT_SPAN)),
    "stop_loss_exp": STOP_LOSS,
    "success_model_after_stop": SUCCESS_MODEL_FLOAT,
    "success_claim": str(CLAIM_SUCCESS),
    "dp_entries_log2": float(exact.log2(DP_ENTRIES)),
    "trails": TRAILS,
    "H_total": str(H_TOTAL),
    "R_total_log2": float(exact.log2(R_TOTAL)),
    "R_over_C_log2": float(exact.log2(R_TOTAL / C)),
    "time_log2": str(time_log2),
    "time_log2_float": float(exact.log2(W)),
    "score_rational": {"numerator": str(time_log2.numerator),
                       "denominator": str(time_log2.denominator)},
    "preprocessing_log2": str(preprocessing_log2),
    "memory_log2_bytes": str(memory_log2),
    "memory_bytes": peak_bytes,
    "advice_bytes": 0,
    "items": [{"name": n, "phase": ph, "H_units": H, "R_instructions": str(R)}
              for n, ph, H, R in ITEMS],
}

if __name__ == "__main__":
    print(json.dumps(out, indent=2))
