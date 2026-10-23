"""ledger.py -- v17 numbers.  DRAFT: the toy does not yet certify (see report).

The v17 step drops v16's per-step canonicalisation and coefficient twist:

    state = (S, (a, c)),  S affine,  S = [a]G + [c]T
    K     = |y(S)|                      # class key, invariant under the six maps
    m     = (K*R) >> (bitlen(P)-1)      # exact top-bit cell, R = 2048
    step: S <- S + J_m,  a <- a + r_m,  c <- c + s_m

There is no beta*x per step (v16 charged 383 for it) and no per-step twist
(v16 charged 138).  The orbit element e is resolved at REPORT time from the two
stored points: two Fp multiplies (~766 instructions) once per record hit.

Batch amortisation, derived (not copied): one Fp inversion is 254 squarings +
11 multiplies = 265 Fp multiplies = 265*383 = 101,495 instructions; over a batch
of B = 2^18 ready walks that is 0.387 instructions per step, and Montgomery's
trick adds 3 multiplies per element, all inside the reference operation's own
12M core, which the cost model forbids charging twice.  The per-step rows below
charge the batch bookkeeping (scan, meter charge, zero-denominator test) as 2,
i.e. 24/2^18 rounded up, so the inversion never appears as a per-step row.
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
N = (r - 1) // 6
C = 2000
R_TAB = 2048
R_EFF = selector.reachable_columns(R_TAB, p)      # 1665 rows are reachable
d = 60
DP = 2 ** d
L = 40 * DP
workers = 2 ** 20
nodes = 2 ** 16
record_bytes = 192
node_cap_records = 2 ** 51
BATCH = 2 ** 18
M_insn = 383
A_insn = 20

SUCCESS_FLOOR = F(39, 100)
MODEL_SUCCESS = F(3903, 10000)

# ---- per-step rows ---------------------------------------------------------
ROWS = [
    ("class_key_abs_y", "6-limb compare y with P-y (7); 6-limb conditional negate (7)", 14),
    ("cell_top_bits", "shift K*R by 64-bit limbs, 2 shifts + 2 moves, no mask (cell < R by construction)", 4),
    ("table_address_and_load", "base+(m*96) address (2); 12 limb loads (12)", 14),
    ("batch_slot_copy_in_out", "12 limb stores into the batch slots (12); 12 loads back (12); 4 slot-address adds (4); 12 batch-index bookkeeping (12)", 40),
    ("coefficient_pair_add", "a+r_m, c+s_m mod r: 4-limb add + conditional subtract (2 x 10); pair store (2)", 22),
    ("dp_mask_and_checkpoint", "DP mask test (3); checkpoint age + 64 compare (4); amortised point confirm (4)", 11),
    ("length_counter_and_cap", "69-bit length increment ldp/add/adc/stp (4); cap test cmp/b.hs (2); subs+b.ne per two unrolled steps (1)", 7),
    ("batch_assembly_inversion_meter_share", "per batch of 2^18: ready-scan (6), zero-denominator OR+cbz (7), meter charge 192-bit add (10) = 23, charged 2 per step (2^-18 * 23 < 0.0001)", 2),
    ("record_store_branch", "DP store taken 2^-60 of steps; 22-instruction write charged 0, branch charged 1", 1),
    ("loop_dispatch_and_spill", "loop counter/branch (5); spill-reload ceiling (10); unitemised control (3)", 18),
]
STEP_INSN = sum(n for _n, _o, n in ROWS)
STEP_PRICE = 1 + F(STEP_INSN, C)

dbl_instructions = 4 * M_insn + 6 * A_insn + 3 * 48 + 56
DOUBLING_PRICE = F(dbl_instructions, C)

# ---- inversion share, derived ---------------------------------------------
INV_MULTS = 254 + 11                    # Fermat: 254 squarings + 11 multiplies
INV_INSN = INV_MULTS * M_insn           # 101,495
INV_PER_STEP = F(INV_INSN, BATCH)       # 0.3873 instructions per step

# ---- success model (unchanged from v16) -----------------------------------
getcontext().prec = 80
_z0 = Decimal(1) - Decimal(MODEL_SUCCESS.numerator) / Decimal(MODEL_SUCCESS.denominator)
_ln0 = -_z0.ln()
L_req = F(int(_ln0 * Decimal(10) ** 70) + 1, 10 ** 70)
Nprime = F(N * R_EFF, R_EFF - 1)


def least_s(success_frac, mu, ln_hi):
    z = Decimal(1) - Decimal(success_frac.numerator) / Decimal(success_frac.denominator)
    ln = -z.ln()
    hi = F(int(ln * Decimal(10) ** 70) + 1, 10 ** 70)
    x2 = 2 * mu * Nprime * hi
    s = math.isqrt(x2.numerator // x2.denominator)
    if F(s * s) < x2:
        s += 1
    return s


def rho_ideal(a):
    v = (-(Decimal(a.numerator) / Decimal(a.denominator)) ** 2 * _ln0).exp()
    return F(int((Decimal(1) - v) * Decimal(10) ** 70), 10 ** 70)


S0 = least_s(MODEL_SUCCESS, F(1), L_req)
TRANSPORT_A_PRIORI = F(2 * DP, S0)
A_MIN = F(1)
while rho_ideal(A_MIN) < SUCCESS_FLOOR + TRANSPORT_A_PRIORI:
    A_MIN += F(1, 1000)
A_H1 = A_MIN + F(5, 1000)
S_ops = math.isqrt((A_H1 * A_H1 * F(S0) * F(S0)).numerator //
                   (A_H1 * A_H1 * F(S0) * F(S0)).denominator)
if F(S_ops * S_ops) < A_H1 * A_H1 * F(S0) * F(S0):
    S_ops += 1

# A CAP, not a mean: the DP count is Binomial(S_ops, 2^-60) with mean S_ops/2^60;
# 4 standard deviations are added explicitly so no printed row rests on an average,
# and the <= 2^-10 probability that the count exceeds this cap is the pool_overflow
# term charged below (same Chernoff object).
mean_records = S_ops // DP
records = mean_records + 4 * math.isqrt(mean_records) + workers

PHASE1 = STEP_PRICE * (S_ops + 2 ** 60)
table_preprocessing = F(6 * R_TAB)
node_setup = F(nodes) * F(2 ** 20)
meter_block_remainders = F(2 ** 61)
routing_report_traffic = F(records) * F(2 ** 10)
verification_output = F(2 ** 40)
fresh_starts = F(workers) * F(16)
drain = F(workers) * F(2 * DP)
report_path = F(records) * F(800)           # e-resolution: 2 Fp muls per hit
# Batch inversion, charged as a TOTAL, not per step and not folded into a row of
# 2 instructions: every batch of BATCH ready walks costs one Fermat inversion
# (254 squarings + 11 multiplies = 265 Fp multiplies = 101,495 instructions).
# ceil((S_ops + 2^60)/BATCH) full batches cover every step of the charged cap, and
# one further whole batch pays the tail batch that closes the last partial group,
# so a run whose final batch holds a single walk is still fully charged.
n_batches = ((S_ops + 2 ** 60) + BATCH - 1) // BATCH
batch_inversion_total = F(INV_INSN * (n_batches + 1), C)   # instructions -> reference units
EXTRAS = (table_preprocessing + node_setup + meter_block_remainders +
          routing_report_traffic + verification_output + fresh_starts + drain +
          report_path)
W = PHASE1 + EXTRAS + batch_inversion_total

transport = F(2 * DP, S_ops)
pool_overflow = F(1, 2 ** 10)
fingerprint_union = F(records * (records - 1), 2 * 2 ** 384)
relation_failure = F(1, r)
detector_abort = F(8 * L * L, N)
LOSSES = [("transport_geometric", transport), ("pool_overflow", pool_overflow),
          ("fingerprint_union", fingerprint_union), ("relation_failure", relation_failure),
          ("detector_abort", detector_abort)]
SUCCESS = rho_ideal(A_H1) - sum(v for _k, v in LOSSES)

MEMORY_BYTES = nodes * node_cap_records * record_bytes + 6 * R_TAB * 8 + 2 ** 40


def log2f(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


CLAIM_LOG2 = math.ceil(log2f(W) * 1000) / 1000.0
CLAIM = F(int(round(CLAIM_LOG2 * 1000)), 1000)

if __name__ == "__main__":
    print("R_eff              ", R_EFF, "of", R_TAB)
    print("log2 S0            ", f"{log2f(S0):.8f}")
    print("A_H1               ", float(A_H1), " log2 S_ops", f"{log2f(S_ops):.8f}")
    print("step insn / price  ", STEP_INSN, f"{float(STEP_PRICE):.6f}")
    print("inversion share    ", f"{float(INV_PER_STEP):.4f}", "instr/step")
    print("doubling price     ", f"{float(DOUBLING_PRICE):.6f}")
    print("log2 phase1        ", f"{log2f(PHASE1):.8f}")
    print("log2 extras        ", f"{log2f(EXTRAS):.8f}")
    print("log2 W             ", f"{log2f(W):.8f}")
    print("claim              ", f"{CLAIM_LOG2:.3f}")
    print("success            ", f"{float(SUCCESS):.10f}")
    print("records cap        ", records)
    print("batches            ", n_batches + 1, " log2 batch inversion", f"{log2f(batch_inversion_total):.6f}")
    print("log2 memory        ", f"{log2f(MEMORY_BYTES):.6f}")
    assert SUCCESS >= SUCCESS_FLOOR
    assert DOUBLING_PRICE < STEP_PRICE
