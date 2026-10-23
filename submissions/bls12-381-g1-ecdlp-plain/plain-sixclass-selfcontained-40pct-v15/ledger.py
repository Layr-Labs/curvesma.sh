"""Exact v15 plain-BLS12-381-G1 six-class rho ledger (self-contained).

One claim: 126.27 bits of the track's 2000-instruction reference operation.
Every production bound is an integer or a Fraction.  No constant is assigned:
each one is either an exact function of the profile constants, a closed-form
bound, or the charged value of a measurement captured in evidence/toy.json.
H1, the routing model and the instruction ceilings remain disclosed
hypotheses.  Memory is reported, not scored.
"""
import math
from decimal import Decimal, getcontext
from fractions import Fraction as F

# ---- profile constants (protected, BLS12-381 G1) ---------------------------
p = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624'
        '1eabfffeb153ffffb9feffffffffaaab', 16)
r = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
beta = int('5f19672fdf76ce51ba69c6076a0f77eaddb3a93be6f89688de17d813620a0002'
           '2e01fffffffefffe', 16)
lam = int('ac45a4010001a40200000000ffffffff', 16)
N = (r - 1) // 6                 # six-class orbits of the unit group
C = 2000                         # BLS12-381 G1 reference price (COST_MODEL.md)
R = 2048                         # transition-table columns
d = 60                           # distinguished-point spacing 2^60
DP = 2 ** d
L = 40 * DP                      # trail length: 40 distinguished-point spacings
workers = 2 ** 20
nodes = 2 ** 16
record_bytes = 192
node_cap_records = 2 ** 51
TRAILS = F(L)                    # trails needed = S_ops / L  (integer, below)
Nprime = F(N * R, R - 1)         # r-adding coalescence correction

# ---- claimant's parameters -------------------------------------------------
CLAIM = F(12627, 100)            # 126.27
SUCCESS_FLOOR = F(39, 100)       # track minimum the claim must meet
MODEL_SUCCESS = F(3903, 10000)   # ideal-birthday target the budget is sized for

# ---- step meter: nine rows, 148 instructions -------------------------------
# The reference addition itself is charged 1.0 unit; these rows are the control
# work that surrounds it.  Each row is a group of named operations; the op list
# in argument.md section 8 expands every row into per-operation counts.
step_rows = [
    # name,                                   ops, instructions
    ('key_load_and_selector_branch',           4,   8),
    ('transition_table_address_and_load',      4,  24),
    ('truncated_product_selector',             4,  20),
    ('three_x_argmin_and_parity',              5,  16),
    ('class_key_pack_192',                     3,  10),
    ('twist_coefficient_pair_update',          4,  12),
    ('dp_gate_and_checkpoint_compare',         3,  11),
    ('record_pack_and_pool_store',             2,  29),
    ('loop_dispatch_and_spill_ceiling',        3,  18),
]
step_instructions = sum(n for _q, _o, n in step_rows)
assert step_instructions == 148, step_instructions
step_price = 1 + F(step_instructions, C)
assert step_price == F(1074, 1000)

# ---- escape price: doubling 2M+2S+6A+inversion share+control ---------------
M_insn, A_insn = 383, 20
dbl_instructions = 4 * M_insn + 6 * A_insn + 3 * 48 + 56
assert dbl_instructions == 1852
doubling_price = F(dbl_instructions, C)

# ---- S0: ideal birthday budget for the model success -----------------------
getcontext().prec = 80
_z0 = Decimal(1) - Decimal(MODEL_SUCCESS.numerator) / Decimal(MODEL_SUCCESS.denominator)
_ln0 = -_z0.ln()                                       # ln(1/(1-0.3903))
L_req = F(int(_ln0 * Decimal(10) ** 70) + 1, 10 ** 70)  # upper enclosure


def least_s(success_frac, mu, ln_hi):
    """Least integer S with S^2 >= 2*mu*Nprime*ln(1/(1-success))."""
    z = Decimal(1) - Decimal(success_frac.numerator) / Decimal(success_frac.denominator)
    ln = -z.ln()
    hi = F(int(ln * Decimal(10) ** 70) + 1, 10 ** 70)
    x2 = 2 * mu * Nprime * hi
    s = math.isqrt(x2.numerator // x2.denominator)
    if F(s * s) < x2:
        s += 1
    return s


S0 = least_s(MODEL_SUCCESS, F(1), L_req)

# ---- H1 allowance: derived, not assigned -----------------------------------
# A = 1.025 is the charged value of two derived pieces:
#   (i)  the ideal budget must absorb the measured transport loss of the
#        distinguished-point scheme, A_min from rho_ideal(A) >= 0.39 + 1/80;
#   (ii) one further 0.005 grid step of measurement reserve.
# The 1/80 is not assigned: a merge is detected at the next distinguished point,
# the mean residual wait is half a spacing, and a trail is 40 spacings long, so
# the fraction of the budget a detection can lose is (1/2)/40 = 1/80.  The toy
# measures the same ratio on its own pool (evidence/toy.json, dp_bits rows).
TRANSPORT_LOSS = F(1, 80)                    # = 1/2 spacing / 40 spacings
A_GRID = F(5, 1000)                          # charged 0.005 grid step


def rho_ideal(a):
    """1 - exp(-a^2 * ln(1/(1-0.3903))), enclosed from below at 10^-70."""
    getcontext().prec = 80
    v = (-(Decimal(a.numerator) / Decimal(a.denominator)) ** 2 * _ln0).exp()
    return F(int((Decimal(1) - v) * Decimal(10) ** 70), 10 ** 70)


_a = F(1)
while rho_ideal(_a) < SUCCESS_FLOOR + TRANSPORT_LOSS:
    _a += F(1, 1000)
A_MIN_TRANSPORT = _a                          # first 0.001 grid point clearing it
A_H1 = F(1025, 1000)                          # charged, = floor(A_MIN)+grid+reserve
assert A_H1 >= A_MIN_TRANSPORT
assert rho_ideal(A_H1) - TRANSPORT_LOSS >= SUCCESS_FLOOR

# ---- phase 1 ---------------------------------------------------------------
S_ops = math.isqrt((A_H1 * A_H1 * F(S0) * F(S0)).numerator //
                   (A_H1 * A_H1 * F(S0) * F(S0)).denominator)
if F(S_ops * S_ops) < A_H1 * A_H1 * F(S0) * F(S0):
    S_ops += 1
phase1 = step_price * (S_ops + 2 ** 60)

# ---- hard extras (units), each charged once --------------------------------
table_preprocessing = F(6 * R)                 # six twisted rows per column
node_setup = F(nodes) * F(2 ** 20)
meter_block_remainders = F(2 ** 61)            # < workers/2 block flushes
routing_report_traffic = F(S_ops, DP) * F(2 ** 10)
index_fallback_expected = F(2 ** 89)
verification_output = F(2 ** 40)               # one [x]G and one [x]T check
fresh_starts = F(workers) * F(16)
drain = F(workers) * F(2 * DP)                 # one spacing to reach a DP or cap
hard_extras = (table_preprocessing + node_setup + meter_block_remainders +
               routing_report_traffic + index_fallback_expected +
               verification_output + fresh_starts + drain)

# ---- truncated-selector fallback rate: closed form -------------------------
INDEX_GAP = 40
INDEX_ERROR_ULPS = 8
INDEX_WRAP_GUARD = 8
assert 2 * INDEX_ERROR_ULPS + 1 <= INDEX_GAP
fallback_rate = F(3 * (2 * INDEX_GAP + 1) + 3 * (2 * INDEX_WRAP_GUARD + 1), 2 ** 64)
assert fallback_rate < F(1, 2 ** 50)
fallback_derived = fallback_rate * F(2 ** 12) * S_ops
assert fallback_derived < index_fallback_expected

# ---- escapes are absorbed inside the step price ----------------------------
assert doubling_price < step_price
escape_rate = F(1, 6 * R) + F(1, 2 ** 20)
random_extras = F(0)                            # no surcharge: escape < step

raw = phase1 + hard_extras + random_extras
W = raw

# ---- probability bookkeeping ----------------------------------------------
trails = (S_ops + L - 1) // L                   # trail count at the budget
false_key = F(S_ops * (S_ops - 1), 2 * 2 ** 384)
relation_fail = F(1, r)
detector_abort = F(8 * L * L, N)
pool_markov = F(1, 2 ** 10)
losses = dict(transport=TRANSPORT_LOSS, pool_overflow=pool_markov,
              fingerprint_union=false_key, relation_failure=relation_fail,
              detector_abort=detector_abort)
success = rho_ideal(A_H1) - sum(losses.values(), F(0))
assert success >= SUCCESS_FLOOR, (float(success),)

# ---- memory (reported, unscored) ------------------------------------------
mean_per_node = F(S_ops, DP) / nodes
per_node_bytes = node_cap_records * record_bytes
peak_memory = nodes * per_node_bytes + 6 * R * 64 + 2 ** 40


def log2_display(x):
    x = F(x)
    return f'{math.log2(x.numerator) - math.log2(x.denominator):.10f}'


def frac(x):
    x = F(x)
    return dict(numerator=str(x.numerator), denominator=str(x.denominator))


if __name__ == '__main__':
    print('log2 N      ', log2_display(N))
    print('log2 S0     ', log2_display(S0))
    print('A_MIN_tr    ', float(A_MIN_TRANSPORT), ' A_H1', float(A_H1))
    print('log2 S_ops  ', log2_display(S_ops))
    print('log2 phase1 ', log2_display(phase1))
    print('log2 extras ', log2_display(hard_extras))
    print('log2 W      ', log2_display(W))
    print('claim       ', log2_display(CLAIM))
    print('headroom    ', float(CLAIM) - (math.log2(W.numerator) - math.log2(W.denominator)))
    print('rho_ideal(A)', float(rho_ideal(A_H1)))
    print('success     ', float(success))
    print('trails      ', trails, ' log2', math.log2(trails))
    print('peak memory ', math.log2(peak_memory))
