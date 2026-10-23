"""Exact v13e six-class (negation+GLV) rho ledger; global work meter, no eviction.

Changes from v13c, all forced by the rejection of 18225d19:
  * detector: prime 67-step window with a per-trail phase, full-point-confirmed
    fingerprint compares, and a minimum-fingerprint escape (cycle cap 66);
  * walk inefficiency: pooled over six toy configurations (two curve sizes,
    four table sizes), mu = 1.107 = ceil((1.069618 + 2*0.018344)*1000)/1000;
  * ordinary step: 147 itemised instructions, including the 192-bit key and
    the 20-bit distinguished-point gate built from the truncated products, so
    no exact beta*x is needed outside the metered near-tie fallback;
  * claim convention: the retired +0.10 family floor is replaced by the
    frontier convention claim = W + tiny epsilon, chosen as the smallest
    4-decimal value that also meets the 0.39 success floor.
C and R/C follow tracks/docs/COST_MODEL.md; memory is reported, not scored.
Heuristic H1 and the instruction ceilings are disclosed, not proved.
"""
import math
from decimal import Decimal, getcontext
from fractions import Fraction as F

p = int('fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f', 16)
r = int('fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141', 16)
beta = int('7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee', 16)
lam = int('5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72', 16)

N = (r - 1) // 6
C = 1000
E0 = 8000
E = 9600
R = 2 ** 20
d = 69
DP = 2 ** d
L = 40 * DP
workers = 2 ** 20
nodes = 2 ** 16
record_bytes = 192
node_cap_records = 2 ** 44        # 4x the mean load; Chernoff row below
M = nodes * node_cap_records
S = 2 ** 85
Nprime = F(N * R, R - 1)
WALK_MU = F(1107, 1000)           # pooled toy estimate, see mu_charge below
CHECKPOINT_PERIOD = 67            # prime: gcd(67, L) = 1 for every L <= 66
CHECKPOINT_CYCLE_CAP = 66
checkpoint_detection_bound = CHECKPOINT_PERIOD * CHECKPOINT_CYCLE_CAP  # 4422

# ---- step price: measured aarch64 worst-case Ir / reference addition (6M = 1074) --
# The counted step is the hand-scheduled v15c kernel: 6 fe_mul (4x4 schoolbook +
# Solinas fold, 107 Ir each including the canonical tail) + 6 fe_sub (9 Ir each)
# + 76 Ir of index/DP/parity/store/loop control = 772 Ir, machine-counted by
# count_iter.py over the compiled step (counts-step-v15c.json).  Price is the
# pure ratio to the reference batched affine addition (6M = 1074 Ir).
instruction_primitives = {'insn': 1, 'fp_mul': 179, 'normalize52_pack': 80,
                          'modsub52_clean': 16, 'unpack52': 15, 'cmp256': 5}
itemized_step_rows = [
    # name, routine, primitive, count
    ('asm_fe_mul_x6', 'step', 'insn', 642),
    ('asm_fe_sub_x6', 'step', 'insn', 54),
    ('index_dp_parity_store_loop', 'step', 'insn', 76),
]
step_rows = {q[0]: q[3] * instruction_primitives[q[2]] for q in itemized_step_rows}
step_instructions = sum(step_rows.values())
assert step_instructions == 772, step_instructions
step_price = F(step_instructions, 1074)          # 772/1074 = 386/537
assert step_price == F(386, 537), step_price

# ---- doubling price: dbl-2009-l (7M) + 6A + inversion share (3M) + control ----
M_insn = 200
A_insn = 20
dbl_instructions = 7 * M_insn + 6 * A_insn + 3 * M_insn + 200
doubling_price = F(235, 100)
assert dbl_instructions == 2320

# ---- S_ops: 0.392 success budget at the pooled mu ---------------------------
getcontext().prec = 80
ln_inv_0608 = -(Decimal(608) / Decimal(1000)).ln()
L_hi = F(int(ln_inv_0608 * Decimal(10) ** 70) + 1, 10 ** 70)
L_lo = F(int(ln_inv_0608 * Decimal(10) ** 70), 10 ** 70)
ln_enclosure_upper = L_hi
ln_enclosure_lower = L_lo
X2 = 2 * WALK_MU * Nprime * L_hi
S_ops = math.isqrt(X2.numerator // X2.denominator)
if F(S_ops * S_ops) < X2:
    S_ops += 1
phase1 = step_price * (S_ops + 2 ** 60)

# ---- pooled mu calibration (evidence/mu-sweep.json) -------------------------
# Six census rows, 8,800 trials: q25 R in {256,1024,2048,4096} and q29 R in
# {1024,4096}.  p_i = P(collision within S_ops(mu=1)); mu_i = ln(1/0.58)/(-ln(1-p_i));
# inverse-variance pooled mu_hat = 1.069618, SE = 0.018344, so
# mu_charge = ceil((mu_hat + 2*SE)*1000)/1000 = 1.107, carried above.
MU_HAT = F(1069618, 1000000)
MU_SE_POOLED = F(18344, 1000000)
assert WALK_MU == F(math.ceil(((MU_HAT + 2 * MU_SE_POOLED) * 1000)), 1000)

# ---- hard extras (units), all charged once ----------------------------------
table_preprocessing = F(2 ** 26)
node_setup = F(nodes) * F(2 ** 20)
meter_block_remainders = F(2 ** 61)
routing_report_traffic = F(2 ** 68)
index_fallback_expected = F(2 ** 89)
verification_output = F(2 ** 40)
fresh_starts = F(workers) * F(16)
hard_extras = (table_preprocessing + node_setup + meter_block_remainders +
               routing_report_traffic + index_fallback_expected +
               verification_output + fresh_starts)

# ---- truncated-index fallback: rate derived from the rule's own error bound --
INDEX_GAP = 40
INDEX_ERROR_ULPS = 8
INDEX_WRAP_GUARD = 8
assert 2 * INDEX_ERROR_ULPS + 1 <= INDEX_GAP
index_fallback_rate_bound = F(3 * (2 * INDEX_GAP + 1) + 3 * (2 * INDEX_WRAP_GUARD + 1),
                              2 ** 64)
assert index_fallback_rate_bound < F(1, 2 ** 50)
index_fallback_derived = index_fallback_rate_bound * F(2 ** 12) * S_ops
assert index_fallback_derived < index_fallback_expected

# ---- random extras: expected doubling surcharge -----------------------------
two_cycle_rate = F(1, 6 * R)
other_cycle_rate = F(1, 2 ** 20)
doubling_rate = two_cycle_rate + other_cycle_rate
surcharge = doubling_price - step_price
random_extras = doubling_rate * surcharge * S_ops

raw = phase1 + hard_extras + random_extras

# ---- claim: smallest 4-decimal value meeting work AND success floors --------
# 2^claim is enclosed from below; the Markov headroom is
# 2^claim - phase1, so a smaller claim means a larger early-halt debit.
_POW2 = {}


def pow2_claim_floor(c):
    """Lower Fraction bound on 2^(n/10000) with 60+ correct digits."""
    c = F(c)
    if c not in _POW2:
        getcontext().prec = 100
        v = (Decimal(c.numerator) / Decimal(c.denominator) * Decimal(2).ln()).exp()
        _POW2[c] = F(int(v * Decimal(10) ** 40), 10 ** 40)
    return _POW2[c]


losses = dict(pool_overflow=F(1, 2 ** 40),
              fingerprint_union=F(S * (S - 1), 2 * 2 ** 192),
              relation_failure=F(1, r),
              sampling_rejections=F(4 * (S + R) * (2 ** 256 - r) ** 4, 2 ** 1024),
              unforced_detector_abort=F(8 * L * L, N),
              endpoint_transport=F(64 * DP * DP, N))
success_model = F(392, 1000)


def markov_at(c):
    return random_extras / (pow2_claim_floor(c) - phase1)


def success_at(c):
    return success_model - markov_at(c) - sum(losses.values(), F(0))


getcontext().prec = 100
_log2raw = ((Decimal(raw.numerator) / Decimal(raw.denominator)).ln() /
            Decimal(2).ln())
_start = int(_log2raw * 10000)
cand = F(_start, 10000)
while cand < F(127):
    if (pow2_claim_floor(cand - F(2, 100000)) >= raw
            and success_at(cand) >= F(391, 1000)):
        break
    cand += F(1, 10000)
claim = cand
assert 126 < claim < 127
assert pow2_claim_floor(claim - F(2, 100000)) >= raw
assert success_at(claim) >= F(391, 1000)
assert claim < F(127)
markov_halt = markov_at(claim)
success = success_at(claim)

# ---- inefficiency tolerance: largest mu that still meets 0.39 ----------------
required_main = F(39, 100) + markov_halt + sum(losses.values(), F(0))
_req = Decimal(required_main.numerator) / Decimal(required_main.denominator)
_ln_req = -(Decimal(1) - _req).ln()
L_req_lo = F(int(_ln_req * Decimal(10) ** 70), 10 ** 70)
mu_tolerance_hi = F(S_ops) ** 2 / (2 * Nprime * L_req_lo)
assert mu_tolerance_hi > WALK_MU

# ---- memory (reported, unscored) --------------------------------------------
memory_rows = dict(node_records=M * record_bytes, auxiliary_cap=2 ** 40,
                   jump_table=6 * R * 64)
peak_memory = sum(memory_rows.values())
per_node_bytes = node_cap_records * record_bytes
mean_per_node = F(2 ** 58, nodes)


def frac(x):
    x = F(x)
    return dict(numerator=str(x.numerator), denominator=str(x.denominator))


def log2_display(x):
    x = F(x)
    return f'{math.log2(x.numerator) - math.log2(x.denominator):.10f}'
