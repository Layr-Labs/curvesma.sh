"""Exact accounting for the claim in proof.md (stdlib only; prints JSON).

W(n) = H + R/C for n >= 2^20 (the value at n = 2^20 is scored). H counts batched
incomplete affine additions (reference operations); R counts other scalar AArch64
instructions. C = 2000 is the BLS12-381 G1 price from tracks/docs/COST_MODEL.md.
All integers are exact; time_log2 is rounded up at six decimals by exact.log2_up.
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
N = (r - 1) // d                  # |H|, H = d-th powers in Z_r^*
assert (r - 1) % d == 0 and N % 6 == 0
M = N // 6                        # classes under mu6 = {+-1, +-lam, +-lam^2}
assert r.bit_length() == 255

# Fixed-base comb for each of the two bases (P0 and A_d): signed windows.
WINDOW = 16
WINDOWS = -(-r.bit_length() // WINDOW)          # 16
TABLE_POINTS = WINDOWS * 2 ** (WINDOW - 1)       # per base
TABLE_BYTES_PER_POINT = 192                      # affine (x, y), 2 x 48 bytes

# Walk parameters (cheon.py): tame and wild trails, distinguished points, caps.
DP_BITS = 90                       # distinguished with probability 2^-90
CAP_BITS = 10                      # trail cap 2^(DP_BITS + CAP_BITS) = 2^100 steps
CAP = 1 << (DP_BITS + CAP_BITS)
DP_ENTRY_BYTES = 48                # label (8) + exponent (32) + family/bookkeeping
MULTIPLIERS = 64                   # walk multipliers (free in cost)

# Success model: S total points, tame-wild collision probability 1 - exp(-S^2/(4M)).
TAU_MODEL = Fraction(3, 5)         # model success 1 - exp(-0.6) = 0.45119
S_COLL = exact.iroot_ceil(-(-(4 * M * TAU_MODEL.numerator) // TAU_MODEL.denominator), 2)
SUCCESS_MODEL_FLOAT = 1 - math.exp(-float(TAU_MODEL))
SUCCESS_CLAIM = Fraction(40, 100)  # margin 0.02 below the model for H1 finite-size effects

# Trails: 2K trails of length CAP cover S_COLL steps; detection waits 1/theta per trail.
K_TRAILS = -(-S_COLL // (2 * CAP))  # per family
DETECTION_DELAY = 2 * K_TRAILS * 2 ** DP_BITS
S_TOTAL = S_COLL + DETECTION_DELAY

# Per-step cost: one comb evaluation (WINDOWS batched additions) plus other work.
ADDS_PER_STEP = WINDOWS            # reference operations (6M, 1 unit each)
OTHER_PER_STEP = 2800              # measured 2703.56 instructions per step on blst
                                   # (evidence/measure-blst-step/summary.json), plus
                                   # about 100 for walk-loop and trail bookkeeping
                                   # that the kernel does not execute

# Setup and per-trail work, and the square-root chain at the end.
SETUP_ADDS = 2 * TABLE_POINTS      # both comb tables, built by additions
TRAIL_STARTS = 2 * K_TRAILS        # each starts with one comb evaluation
CHAIN_EVALS = 2 * K_EXP + 1        # two candidate checks per level plus the final check
CHAIN_OTHER = 2 ** 20              # Tonelli-Shanks and comparisons, bound
SOLVE_OTHER = 2 ** 16              # collision solve, six candidate checks, inversion

H_total = (S_TOTAL * ADDS_PER_STEP
           + SETUP_ADDS
           + (TRAIL_STARTS + CHAIN_EVALS) * WINDOWS)
R_total = S_TOTAL * OTHER_PER_STEP + CHAIN_OTHER + SOLVE_OTHER
W = Fraction(H_total) + Fraction(R_total, C)
time_log2 = exact.log2_up(W)

pre_H = SETUP_ADDS
pre = Fraction(pre_H) + Fraction(0)
preprocessing_log2 = exact.log2_up(pre)

dp_entries = -(-S_TOTAL // 2 ** DP_BITS)
peak_bytes = (2 * dp_entries * DP_ENTRY_BYTES            # load factor 2
              + 2 * TABLE_POINTS * TABLE_BYTES_PER_POINT  # both comb tables
              + 2 * K_TRAILS * 64                         # trail state
              + 2 ** 24)                                  # code, stacks, scratch
memory_log2 = exact.log2_up(Fraction(peak_bytes))

out = {
    "r_bits": r.bit_length(),
    "N_log2": float(exact.log2(N)),
    "M_log2": float(exact.log2(M)),
    "window": WINDOW,
    "windows": WINDOWS,
    "multipliers": MULTIPLIERS,
    "tau_model": str(TAU_MODEL),
    "success_model": SUCCESS_MODEL_FLOAT,
    "success_claim_lower_bound": str(SUCCESS_CLAIM),
    "S_coll_log2_over_sqrtM": float(exact.log2(S_COLL) - exact.log2(M) / 2),
    "S_coll_log2": float(exact.log2(S_COLL)),
    "S_total_log2": float(exact.log2(S_TOTAL)),
    "trails_per_family_log2": float(exact.log2(K_TRAILS)),
    "H_log2": float(exact.log2(H_total)),
    "R_over_C_log2": float(exact.log2(Fraction(R_total, C))),
    "time_log2": str(time_log2),
    "time_log2_float": float(exact.log2(W)),
    "score_rational": {"numerator": str(time_log2.numerator),
                       "denominator": str(time_log2.denominator)},
    "preprocessing_log2": str(preprocessing_log2),
    "memory_log2_bytes": str(memory_log2),
    "dp_entries_log2": float(exact.log2(dp_entries)),
    "table_bytes": 2 * TABLE_POINTS * TABLE_BYTES_PER_POINT,
    "advice_bytes": 0,
}

if __name__ == "__main__":
    print(json.dumps(out, indent=2))
