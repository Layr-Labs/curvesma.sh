"""Exact v13 six-class (negation+GLV) rho ledger; global work meter, no eviction.

Port of the review-passed plain-track construction 1c600df4 to secp256k1 with a
self-itemized instruction ceiling and an explicit per-node storage partition.
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

N = (r - 1) // 6                 # class count of the six-element unit group
C = 1000                         # secp256k1 reference price (COST_MODEL.md)
E0 = 8000                        # declared per-tick ceiling (repriced below)
E = 9600                         # E0 * 6/5 repricing of the whole ceiling
R = 2 ** 20                      # jump table entries
d = 69                           # distinguished-point spacing 2^d
DP = 2 ** d
L = 40 * DP                      # trail cutoff, 40 expected DP spacings
workers = 2 ** 20
nodes = 2 ** 16                  # storage nodes, routed by bits 128..143 of K
record_bytes = 192
node_cap_records = 2 ** 43        # 2.06x the mean load (Chernoff, see below)
M = nodes * node_cap_records      # total live records
S = 2 ** 85                       # event/start cap
Nprime = F(N * R, R - 1)
WALK_MU = F(529, 500)             # measured first-repeat inefficiency, see S_ops

# ---- step price: one reference addition (6M = 1074, charged 1.0) + control ----
instruction_primitives = {'insn': 1, 'fp_mul': 179, 'normalize52_pack': 80,
                          'modsub52_clean': 16, 'unpack52': 15, 'cmp256': 5}
itemized_step_rows = [
    # name, routine, primitive, count
    ('class_key_192_load', 'step', 'insn', 4),
    ('branch_and_twist_extract', 'step', 'insn', 4),
    ('twisted_table_address', 'step', 'insn', 3),
    ('twisted_jump_load', 'step', 'insn', 8),
    ('twist_alignment_select', 'step', 'insn', 2),
    ('meter_tick', 'meters', 'insn', 3),
    ('index_fraction_one_product', 'canonicalize', 'insn', 17),
    ('index_fraction_two_three_fold', 'canonicalize', 'insn', 6),
    ('near_tie_gap_test', 'canonicalize', 'insn', 4),
    ('argmin_three_way', 'canonicalize', 'insn', 6),
    ('y_parity_table_twist', 'canonicalize', 'insn', 3),
    ('class_key_192_store', 'step', 'insn', 4),
    ('two_cycle_compare', 'escape', 'insn', 5),
    ('checkpoint_compare', 'escape', 'insn', 5),
    ('dp_bit_test', 'step', 'insn', 3),
    ('key_rotate', 'step', 'insn', 4),
    ('escape_flag_pending', 'escape', 'insn', 2),
    ('loop_dispatch', 'step', 'insn', 5),
    ('replay_tag_write', 'step', 'insn', 2),
    ('spill_reload_ceiling', 'step', 'insn', 10),
    ('unitemised_overhead_allowance', 'step', 'insn', 10),
]
step_rows = {q[0]: q[3] * instruction_primitives[q[2]] for q in itemized_step_rows}
step_instructions = sum(step_rows.values())
step_price = F(1) + F(step_instructions, C)          # 1 + itemised/1000
assert step_price == F(111, 100), step_price

# ---- doubling price: dbl-2009-l (7M) + 6A + inversion share (3M) + control ----
M_insn = 200                     # 4x64 Montgomery multiply ceiling, limb-counted
A_insn = 20
dbl_instructions = 7 * M_insn + 6 * A_insn + 3 * M_insn + 200
doubling_price = F(235, 100)     # 2320 instructions rounded up
assert dbl_instructions == 2320

# ---- S_ops: 0.392 success budget in the review-passed shape ----
# The accepted-step count is the least integer S with
#     S^2 >= 2*mu*N'*ln(1000/608),   N' = N*R/(R-1),  mu = 529/500 = 1.058.
# mu is the walk's first-repeat inefficiency MEASURED on the production-shaped
# mechanism at the finest captured spacing (DP 2): that row's all-arrival
# collision probability is 0.4025 at the v13b budget S(0.420), so
#     mu = ln(1/0.58) / (-ln(1-0.4025)) = 1.05796 <= 1.058.
# v13b assumed mu = 1 (the ideal-map value) and reserved 0.42 - 0.40 = 0.02 for
# that model error; the same correction is now carried inside the budget, which
# is what makes the smaller S legal rather than a thinner margin on the same S.
getcontext().prec = 80
ln_inv_0608 = -(Decimal(608) / Decimal(1000)).ln()   # ln(1000/608), success 0.392
L_hi = F(int(ln_inv_0608 * Decimal(10) ** 70) + 1, 10 ** 70)   # upper enclosure
L_lo = F(int(ln_inv_0608 * Decimal(10) ** 70), 10 ** 70)       # lower enclosure
ln_enclosure_upper = L_hi        # ln(1000/608) < L_hi  -> S_ops is large enough
ln_enclosure_lower = L_lo        # ln(1000/608) > L_lo  -> S_ops is minimal
X2 = 2 * WALK_MU * Nprime * L_hi
S_ops = math.isqrt(X2.numerator // X2.denominator)
if F(S_ops * S_ops) < X2:
    S_ops += 1
phase1 = step_price * (S_ops + 2 ** 60)

# ---- the inefficiency constant, derived from the captured DP-2 pipeline row ----
# v13b's budget was S(0.420).  Its DP-2 pipeline row (4,000 trials) measured the
# all-arrival collision in 1,610 of them, so the model 1-exp(-S^2/(2 mu N')) gives
# mu = ln(1/0.58) / (-ln(1-0.4025)) = 1.05772, rounded UP to 1.058 for the budget.
_ln_v13b = -(Decimal(1) - Decimal('0.420')).ln()
_ln_dp2 = -(Decimal(1) - Decimal(1610) / Decimal(4000)).ln()
TOY_MU_DERIVED = F(int((_ln_v13b / _ln_dp2) * Decimal(10) ** 5), 10 ** 5)
assert TOY_MU_DERIVED == F(26443, 25000)          # 1.05772
assert TOY_MU_DERIVED <= WALK_MU

# ---- hard extras (units), all charged once ----
table_preprocessing = F(2 ** 26)         # 6*2^20 twisted entries, batched
node_setup = F(nodes) * F(2 ** 20)       # allocation + index build, per node
meter_block_remainders = F(2 ** 61)
routing_report_traffic = F(2 ** 68)      # 2^57.75 reports * 2^10 transfer
index_fallback_expected = F(2 ** 89)     # rate <= 2^-50, 2^12 units each
verification_output = F(2 ** 40)
fresh_starts = F(workers) * F(16)
hard_extras = (table_preprocessing + node_setup + meter_block_remainders +
               routing_report_traffic + index_fallback_expected +
               verification_output + fresh_starts)

# ---- truncated-index fallback: rate derived from the rule's own error bound ----
# The cheap rule computes limb 5 of x*floor(2^384*beta/p) from eight limb
# products.  The carry from column 4 into column 5 is in [0,7], so F1 is at most
# 7 ulps below the exact 64-bit fraction; F0 is exact (p = 2^256 - 2^32 - 977, so
# the top limb is the fraction up to one ulp) and F2 = -(F0+F1) inherits at most
# 8.  With E <= 8 the argmin is decided correctly whenever the computed gap
# exceeds 2E, so a threshold of 40 is sufficient.  A fraction within 8 of 0 or
# 2^64 carries the sign of the wrap and is added to the fallback condition.
INDEX_GAP = 40
INDEX_ERROR_ULPS = 8
INDEX_WRAP_GUARD = 8
assert 2 * INDEX_ERROR_ULPS + 1 <= INDEX_GAP
index_fallback_rate_bound = F(3 * (2 * INDEX_GAP + 1) + 3 * (2 * INDEX_WRAP_GUARD + 1),
                              2 ** 64)
assert index_fallback_rate_bound < F(1, 2 ** 50)
index_fallback_derived = index_fallback_rate_bound * F(2 ** 12) * S_ops
assert index_fallback_derived < index_fallback_expected      # charged at 2^89

# ---- checkpoint scheduling: a cycle of length L <= 64 is left within 64L steps ----
CHECKPOINT_PERIOD = 64
CHECKPOINT_CYCLE_CAP = 64
checkpoint_detection_bound = CHECKPOINT_PERIOD * CHECKPOINT_CYCLE_CAP   # 2^12

# ---- random extras: expected doubling surcharge (Markov-debited) ----
two_cycle_rate = F(1, 6 * R)
other_cycle_rate = F(1, 2 ** 20)         # length>=3 F-cycles, exact sum ceiling
doubling_rate = two_cycle_rate + other_cycle_rate
surcharge = doubling_price - step_price
random_extras = doubling_rate * surcharge * S_ops

raw = phase1 + hard_extras + random_extras
W = F(27, 25) * raw
claim = F(127)                   # 12700/100: the 0.10-above-raw floor, below gate

# ---- probability bookkeeping ----
success_model = F(392, 1000)             # by construction of S_ops (no pad)
pool_markov = F(1, 2 ** 10)              # Chernoff at 2.06x mean, union over nodes
false_key = F(S * (S - 1), 2 * 2 ** 192)
relation = F(1, r)
sampling = F(4 * (S + R) * (2 ** 256 - r) ** 4, 2 ** 1024)
detector_abort = F(8 * L * L, N)
endpoint_transport = F(64 * DP * DP, N)
losses = dict(pool_overflow=pool_markov, fingerprint_union=false_key,
              relation_failure=relation, sampling_rejections=sampling,
              unforced_detector_abort=detector_abort,
              endpoint_transport=endpoint_transport)
markov_halt = random_extras / (2 ** 127 - phase1)
success = success_model - markov_halt - sum(losses.values(), F(0))

# ---- inefficiency tolerance: largest walk constant that still meets 0.39 ----
# Keep S fixed and ask which mu makes the model success equal
# 0.39 + losses.  Written as mu_max = S^2 / (2*N'*(-ln(1-required))), with the
# logarithm enclosed from below so mu_max_hi is an upper bound on the tolerance.
required_main = F(39, 100) + markov_halt + sum(losses.values(), F(0))
_req = Decimal(required_main.numerator) / Decimal(required_main.denominator)
_ln_req = -(Decimal(1) - _req).ln()
L_req_lo = F(int(_ln_req * Decimal(10) ** 70), 10 ** 70)
mu_tolerance_hi = F(S_ops) ** 2 / (2 * Nprime * L_req_lo)
assert mu_tolerance_hi > WALK_MU

# ---- memory (reported, unscored) ----
memory_rows = dict(node_records=M * record_bytes, auxiliary_cap=2 ** 40,
                   jump_table=6 * R * 64)
peak_memory = sum(memory_rows.values())
per_node_bytes = node_cap_records * record_bytes
mean_per_node = F(2 ** 58, nodes)        # 2^58 expected reports / 2^16 nodes


def frac(x):
    x = F(x)
    return dict(numerator=str(x.numerator), denominator=str(x.denominator))


def log2_display(x):
    x = F(x)
    return f'{math.log2(x.numerator) - math.log2(x.denominator):.10f}'
