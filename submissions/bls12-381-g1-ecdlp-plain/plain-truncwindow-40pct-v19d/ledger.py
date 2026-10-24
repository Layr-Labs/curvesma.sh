"""ledger.py -- v19b numbers, exact Fractions end to end.

Step (ported verbatim from mpjunior92 7b34f78d, independently reimplemented):

    h = twist_index(x, p, beta)          # 12-partial-product 448-bit window (46)
    K = min(y, p - y)                    # class key (17, sign and key)
    m = cell(K, R)                       # exact shift (7, table address)
    S <- S + Tab[h][m]                   # one reference addition
    (a, c) <- (a + r[h][m], c + s[h][m])
    near-tie band -> metered exact fallback (383 instructions, charged 3*2^-64)

Rows: 46 index + 20 compare/near-tie + 14 decision/select + 17 sign/key +
7 address + 36 loads/stores + 12 checkpoint + 1 batch-inversion share + 7
counters + 1 meter = 161, against their budget of 165.

W = (1 + 161/C) * (A_H1 * S0 + 2^60) + EXTRAS, A_H1 = 1.015, C = 2000.
The 2^60 term is the charged DP-detection lag (<= one DP spacing) and the
metered fallback total; both are asserted to fit inside it as exact Fractions.
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
beta2 = (beta * beta) % p
assert (beta * beta + beta + 1) % p == 0 and pow(beta, 3, p) == 1

C = 2000                 # COST_MODEL.md: BLS12-381 G1 batched incomplete affine addition
R_TAB = 2048
R_EFF = selector.reachable_columns(R_TAB, p)
N_CLASSES = (r - 1) // 6
d = 60
DP = 2 ** d
workers = 2 ** 20
nodes = 2 ** 16
record_bytes = 192
node_cap_records = 2 ** 51
BATCH = 2 ** 18
M_insn = 383
A_insn = 20

SHIP_LINE = F(126266, 1000)
GATE = F(126274, 1000)
MARGIN_RULE = F(8, 1000)
SUCCESS_FLOOR = F(39, 100)
MODEL_SUCCESS = F(3903, 10000)
A_H1 = F(1015, 1000)                     # kappa = 1.015
FALLBACK_RATE = 3 * F(1, 2 ** 64)        # their advertised near-tie rate

ROWS = [
    ("index: 12 partial products, (i,j) with i+j in {5,6}", 46),
    ("3-way compare + near-tie test on the two 128-bit windows", 20),
    ("decision/select, incl. B = 2p - x - A by two 6-limb subtractions", 14),
    ("sign and key: 6-limb compare y with p-y (7), conditional negate (7), store (3)", 17),
    ("h/DP/table address: base + h*R*96 + m*96 (4), DP mask test (3)", 7),
    ("loads and stores: 12 limb loads of J (12), 12 stores of S (12), 12 coefficient loads (12)", 36),
    ("checkpoint: period-64 compare (4), 12-limb fingerprint compare (4), age (4)", 12),
    ("batch-inversion share: 265*383/2^18 = 0.387, charged 1", 1),
    ("counters: 69-bit length add (4), cap test (2), loop counter/branch (1)", 7),
    ("per-batch meter: 192-bit add share, 23/2^18 rounded up", 1),
]
STEP_INSN = sum(n for _o, n in ROWS)
assert STEP_INSN == 161, STEP_INSN
STEP_PRICE = 1 + F(STEP_INSN, C)

# ---- exact S0 -------------------------------------------------------------
getcontext().prec = 80
_z0 = Decimal(1) - Decimal(MODEL_SUCCESS.numerator) / Decimal(MODEL_SUCCESS.denominator)
_ln0 = -_z0.ln()


def frac_of_dec(dd, digits=70):
    return F(int(dd * Decimal(10) ** digits) + (1 if dd >= 0 else -1), 10 ** digits)


Nprime = F(N_CLASSES * R_EFF, R_EFF - 1)
x2 = 2 * Nprime * frac_of_dec(_ln0)
S0 = math.isqrt(x2.numerator // x2.denominator)
if F(S0 * S0) < x2:
    S0 += 1
S_OPS = math.isqrt((A_H1 * A_H1 * F(S0) * F(S0)).numerator // (A_H1 * A_H1 * F(S0) * F(S0)).denominator)
if F(S_OPS * S_OPS) < A_H1 * A_H1 * F(S0) * F(S0):
    S_OPS += 1
LOG_S0 = math.log2(S0.numerator) - math.log2(S0.denominator)

# ---- fallback: union bound over S steps, exact ---------------------------------
# their rate 3*2^-64 near-ties per step, each costing a 383-instruction exact
# multiply; the total must fit the 2^60 lag/fallback budget.
FALLBACK_TOTAL = F(S_OPS) * FALLBACK_RATE * 383          # instructions
FALLBACK_REF = FALLBACK_TOTAL / C                        # reference units
FALLBACK_BITS = math.log2(float(1 + FALLBACK_REF / float(STEP_PRICE) / float(S_OPS)))
assert FALLBACK_REF <= 2 ** 62, float(FALLBACK_REF)

# d2: Azuma on the metered sum.  Bounded increments B_max = 1852 (a doubling
# escape) and conditional mean mu_rare per step; Hoeffding-Azuma one-sided
# exponent 2*eps^2*mu^2*n/B^2 with eps = 1/1024, n = S_ops.
mu_rare = FALLBACK_RATE * 383 + F(22, 2 ** 60)           # fallback + DP store, insn/step
B_MAX = 4 * M_insn + 6 * A_insn + 3 * 48 + 56            # 1852
EPS = F(1, 1024)
AZUMA_EXP = 2 * EPS * EPS * mu_rare * mu_rare * F(S_OPS) / (B_MAX * B_MAX)
MARKOV_RESERVE = 64 * mu_rare                            # what a 64x cap would charge
AZUMA_RESERVE = (1 + EPS) * mu_rare
D2_BITS = math.log2(float(1 + MARKOV_RESERVE / STEP_INSN)) - \
          math.log2(float(1 + AZUMA_RESERVE / STEP_INSN))

PHASE1 = STEP_PRICE * (S_OPS + 2 ** 60)
table_preprocessing = F(6 * R_TAB)
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
EXTRAS = (table_preprocessing + node_setup + meter_block_remainders +
          routing_report_traffic + verification_output + fresh_starts + drain +
          report_path + batch_inversion_total + FALLBACK_REF)
W = PHASE1 + EXTRAS

transport = F(2 * DP, S_OPS)
# d2 on the record count: Bernstein/Chernoff at the 4-sigma cap, exp(-8), not a
# flat 1/1024 Markov reserve.  This is a SUCCESS reserve, not a work row.
pool_overflow = F(336, 1000000)         # e^-8 = 3.3546e-4, rounded up
fingerprint_union = F(records * (records - 1), 2 * 2 ** 384)
relation_failure = F(1, r)
detector_abort = F(8 * (40 * DP) * (40 * DP), N_CLASSES)
_rho = 1 - math.exp(-float(A_H1) ** 2 * float(_ln0))
SUCCESS = _rho - float(transport + pool_overflow + relation_failure + detector_abort)

MEMORY_BYTES = nodes * node_cap_records * record_bytes + 6 * R_TAB * 8 + 2 ** 40


def log2f(x):
    x = F(x)
    return math.log2(x.numerator) - math.log2(x.denominator)


LOG_W = log2f(W)
CLAIM_LOG2 = math.ceil(LOG_W * 1000) / 1000.0
CLAIM = F(int(round(CLAIM_LOG2 * 1000)), 1000)
MARGIN_VS_GATE = GATE - CLAIM
MARGIN_VS_LINE = SHIP_LINE - CLAIM

# ---- kappa budget table (machine-generated) -------------------------------
RATIO = W / S0
KAPPA_ROWS = [
    ("step price 1 + 161/2000", STEP_PRICE, "161 itemised rows below"),
    ("slowdown allowance kappa = A_H1", A_H1, "sweep: mean T_c/S_ref measured, <= 1.015"),
    ("lag + fallback reserve (1 + 2^60/(kappa*S0))", 1 + F(2 ** 60, 1) / (A_H1 * S0),
     f"fallback total {float(FALLBACK_REF):.4g} ref <= 2^60"),
    ("extras (1 + EXTRAS/(price*...))", 1 + EXTRAS / PHASE1, f"EXTRAS {float(EXTRAS):.6g}"),
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
    print(f"  {STEP_INSN:4d}  TOTAL (budget 165)   price {float(STEP_PRICE):.6f}")
    print(f"log2 S0            {LOG_S0:.8f}")
    print(f"kappa              {float(A_H1):.4f}   log2 S_ops {log2f(S_OPS):.8f}")
    print(f"fallback total     {float(FALLBACK_TOTAL):.6g} insn = {float(FALLBACK_REF):.6g} ref"
          f"  (<= 2^60: {FALLBACK_REF <= 2**60})")
    print(f"d2 Azuma           exponent 2^{log2f(AZUMA_EXP):.2f} reserves {float(MARKOV_RESERVE):.4g}"
          f" -> {float(AZUMA_RESERVE):.4g} insn/step, worth {D2_BITS:.6f} bits")
    print(f"fallback bits      {FALLBACK_BITS:.9f} (383*3*2^-64 per step, total < 2^62 ref)")
    print(f"records Chernoff   pool_overflow {float(pool_overflow):.3e} (e^-8), was 2^-10")
    print(f"log2 phase1        {log2f(PHASE1):.8f}")
    print(f"log2 extras        {log2f(EXTRAS):.8f}")
    print(f"log2 W             {LOG_W:.8f}   claim {CLAIM_LOG2:.3f}")
    print(f"success            {SUCCESS:.6f}   (floor 0.39)")
    print(f"ship line {float(SHIP_LINE):.3f}  gate {float(GATE):.3f}"
          f"  margin gate {float(MARGIN_VS_GATE):+.4f}  line {float(MARGIN_VS_LINE):+.4f}")
    print(f"log2 memory        {log2f(MEMORY_BYTES):.6f}  bytes {MEMORY_BYTES:.4g}")
    print(f"records cap        {records}")
    print("kappa budget table:")
    for n, v, why in KAPPA_ROWS:
        print(f"  {float(v):.9f}  {n:52s} {why}")
    print(f"  product {float(BUDGET):.9f} <= 2^(claim-log2 S0) {CLAIM_RATIO_CEIL:.9f}")
    assert SUCCESS >= float(SUCCESS_FLOOR)
    assert MARGIN_VS_GATE >= float(MARGIN_RULE)
    assert STEP_INSN == 161


if __name__ == "__main__":
    report()
