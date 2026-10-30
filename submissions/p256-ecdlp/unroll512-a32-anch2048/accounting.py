"""Bundled experiment: entrypoint accounting_src.py with vendored cscore.exact.
cscore.exact sha256 ad40a116ac9b6b9d430ae9baa9b18c4c218a6fc996c486d18c4815d650d9e638
Derived from the declared predecessor; blinded-corrections revision (metered corrections, blinded prefix chain, trail/round caps, 2^20 batches); small-key y-window adders; checkpoint-only cycle detection; deferred inv/x3/y3' tests; paired reduction rows."""

import sys as _sys
import types as _types

def _install(_modules):
    for _name, _src in _modules:
        _pkg, _, _leaf = _name.rpartition(".")
        if _pkg and _pkg not in _sys.modules:
            _p = _types.ModuleType(_pkg)
            _p.__path__ = []
            _sys.modules[_pkg] = _p
        _m = _types.ModuleType(_name)
        _m.__file__ = "<bundled %s>" % _name
        _sys.modules[_name] = _m
        if _pkg:
            setattr(_sys.modules[_pkg], _leaf, _m)
        exec(compile(_src, _m.__file__, "exec"), _m.__dict__)

_install([
('cscore.exact', r'''"""Exact accounting helpers for curvesma.sh claims (stdlib only).

Work is expressed in target-group-operation equivalents ("units"): one reference
operation is 1 unit and every other scalar AArch64 instruction is 1/C units.
"""
import math
from fractions import Fraction

def iroot_ceil(n, k):
    """Smallest integer x with x^k >= n (integer Newton iteration)."""
    x = 1 << (n.bit_length() // k + 1)
    while True:
        y = ((k - 1) * x + n // x ** (k - 1)) // k
        if y >= x:
            break
        x = y
    while x ** k < n:
        x += 1
    return x

def pow2_ceil(exponent, denominator=100):
    """Smallest integer >= 2^(exponent), exponent given as a decimal (e.g. 96.75).

    Uses 2^(a/d) >= v  <=>  2^a >= v^d, exact on integers.
    """
    e = Fraction(exponent).limit_denominator(denominator)
    return iroot_ceil(2 ** e.numerator, e.denominator) if e.denominator > 1 else 2 ** e.numerator

def log2(v):
    v = Fraction(v)
    return math.log2(v.numerator) - math.log2(v.denominator)

def log2_up(v, places=6):
    """log2 of an exact positive rational, rounded up at the given decimal place.

    math.log2 on Python integers is accurate to ~1e-15 relative; the 1e-4 ulp
    guard keeps the rounding upward.
    """
    approx = log2(v)
    scale = 10 ** places
    q = math.ceil(approx * scale + 1e-4)
    assert q / scale >= approx
    return Fraction(q, scale)

def units(instructions, C):
    """Convert a scalar instruction count to units at reference price C."""
    return Fraction(instructions, C)

class Ledger:
    """Named cost components in units; reports exact totals and rounded-up log2."""

    def __init__(self, C):
        self.C = C
        self.components = {}
        self.preprocessing = set()

    def add(self, name, value_units, preprocessing=False):
        self.components[name] = Fraction(value_units)
        if preprocessing:
            self.preprocessing.add(name)
        return self

    def add_instructions(self, name, instructions, preprocessing=False):
        return self.add(name, units(instructions, self.C), preprocessing)

    def total(self):
        return sum(self.components.values())

    def preprocessing_total(self):
        return sum((self.components[k] for k in self.preprocessing), Fraction(0))

    def report(self):
        total = self.total()
        pre = self.preprocessing_total()
        score = log2_up(total)
        return {
            "components_log2_units": {k: round(log2(v), 6) for k, v in self.components.items()},
            "time_log2": str(score),
            "time_log2_decimal": float(score),
            "time_log2_float": log2(total),
            "score_rational": {"numerator": str(score.numerator), "denominator": str(score.denominator)},
            "preprocessing_log2": str(log2_up(pre)) if pre > 0 else "0",
        }

def birthday_success(T, N, eps=Fraction(0), delta=Fraction(0)):
    """Random-mapping collision bound: 1 - exp(-T^2/(2N) * (1 - eps)) - delta.

    T: total fresh steps, N: number of classes, eps: relative exponent loss
    (walk penalty, fruitless steps), delta: additive failure terms.
    """
    tau = Fraction(T * T, 2 * N)
    p = 1 - math.exp(-float(tau * (1 - Fraction(eps))))
    return {"tau": float(tau), "eps": float(eps), "model": p, "lower_bound": p - float(delta)}

def tau_for_success(p):
    """tau = T^2/(2N) at which 1 - exp(-tau) = p."""
    return -math.log(1 - p)
'''),

])

"""Exact accounting for proof.md (stdlib only; prints JSON). Every instruction is charged at
1/C, C=1000; integers and rationals are exact, logarithms rounded up at 6 decimals."""
import json
import math
from fractions import Fraction

from cscore.exact import Ledger, iroot_ceil, log2, log2_up

r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
C = 1000
ADDERS = 2 ** 32
GAMMA = 2
K = 2 ** 32
THETA_BITS = 60
# P=2048, untested tail U=1994: a cycle of length l<=64 is caught iff P mod l < P-U = 54 (max 53, l=57).
# CAP_PERIODS=2^55 keeps the per-trail cap at 2048*2^55 = the inherited 2^66 steps.
ANCHOR_P = 2048
ANCHOR_U = 1994
ANCHOR_T = ANCHOR_P - ANCHOR_U
CAP_PERIODS = 2 ** 55
CAP = ANCHOR_P * CAP_PERIODS
TRAIL_MAX = CAP + ANCHOR_P
ANCHOR_R_MAX = max(ANCHOR_P % l for l in range(1, 65))
assert ANCHOR_T == 54 and ANCHOR_R_MAX == 53 and ANCHOR_R_MAX < ANCHOR_T
assert CAP == 64 * 2 ** 60 and TRAIL_MAX == 2 ** 66 + ANCHOR_P

N = (r - 1) // 2 + 1

def long_rate(adders):
    return Fraction(3, 4 * adders ** 2) + Fraction(2, adders ** 3)

PERSIST = 2 * ANCHOR_P + 64   # 2P-1 persistent steps per entry (E2) plus 64 look-ahead states, charged 2P+64
eps = Fraction(4, ADDERS) + 2 * PERSIST * (Fraction(1, 2 * ADDERS) + long_rate(ADDERS))
delta = Fraction(1, 2 ** 59)
TAU_TARGET = Fraction(4943024, 10000000)
T_min = iroot_ceil(math.ceil(2 * N * TAU_TARGET), 2)
L = -(-T_min // K)
T = K * L
S_max = K * (L + TRAIL_MAX)

S_FACTORS = {2: 4, 3: 1, 71: 1, 131: 1, 373: 1, 3407: 1, 17449: 1, 38189: 1}
K_ESC = TRAIL_MAX

def order_dividing(g, n, factors):
    """Order of g mod r, given g^n = 1 and the prime factorisation of n as {prime: exponent}."""
    o = n
    for p, e in factors.items():
        for _ in range(e):
            if pow(g, o // p, r) == 1:
                o //= p
            else:
                break
    return o

S_ORDER = math.prod(p ** e for p, e in S_FACTORS.items())
assert all(all(p % q for q in range(2, math.isqrt(p) + 1)) for p in S_FACTORS)
assert (r - 1) % S_ORDER == 0
U2 = pow(2, (r - 1) // S_ORDER, r)
assert pow(U2, S_ORDER, r) == 1
ORD2_LOWER = order_dividing(U2, S_ORDER, S_FACTORS)
assert pow(U2, ORD2_LOWER, r) == 1
assert all(pow(U2, ORD2_LOWER // p, r) != 1 for p in S_FACTORS if ORD2_LOWER % p == 0)
assert ORD2_LOWER > GAMMA * K_ESC

F_MAX = 41190
FINV_SAVING = 11814
F_MAX_R = F_MAX - FINV_SAVING
ADD_MAX = 41185
SAME_MAX = 32
PDP_MAX = 10
ESC2_MAX = 41262
ESC65_MAX = 43089
OVERHEAD = 1024
UNBATCHED = F_MAX
BATCH = 2 ** 20
# 512-copy jump-table blocks: 4*BATCH/8 - 4*BATCH/512 - 3.  The four-instruction end test once per
# 128 lanes and the fixed three-instruction net dispatch are the measured 128-copy terms; the
# 128 -> 512 step quarters the block-end tests, an analytic 4*(BATCH/128 - BATCH/512).
UNROLL512_SAVING = 516093  # 4*BATCH/8 - 4*BATCH/512 - 3, 128-copy measurement + analytic block-count term for 512-copy blocks
PLAIN_FIXED = 653827732 - (21 * (BATCH - 1) + 4 * (BATCH - 2)) + 2 * (BATCH - 1) - 110 * (BATCH - 1) + 90 * (BATCH - 2) + 160 - 3 * (BATCH - 1) - FINV_SAVING - UNROLL512_SAVING
ANCHOR_EXTRA = 3 * (BATCH - 1) + 9
METHOD_ALLOWANCE = 1024
# Three-arm anchor schedule: t % 2048 == 0 is an anchor round, 1 <= t % 2048 <= 1994 is a
# skip round (plain tail without its anchor compare), the rest are plain rounds.
# The plain arm pays the three-way dispatch (3), the skip arm -2*(BATCH-1) from
# dropping `ldr x6` and `ccmp` per copy plus dispatch (2), and the anchor arm
# ANCHOR_EXTRA - 1 over the plain fixed path. With P=2048 the per-period counts are
# 1994 skip, 53 plain and one anchor round; every cycle length l<=64 has a multiple
# in the tested window because 2048 mod l < 54 for all l<=64 (checked below).
PLAIN_ARM = PLAIN_FIXED + 3 + METHOD_ALLOWANCE
SKIP_ARM = PLAIN_FIXED - 2 * (BATCH - 1) + 2 + METHOD_ALLOWANCE
ANCHOR_ARM = PLAIN_FIXED + ANCHOR_EXTRA - 1 + METHOD_ALLOWANCE

def schedule_counts(rounds):
    a = rounds // ANCHOR_P
    rest = rounds - ANCHOR_P * a
    skip = ANCHOR_U * a + min(ANCHOR_U, rest)
    anchor = a
    return (skip, rounds - skip - anchor, anchor)
RHO_DRAWS = 8
RHO_COST = -(-RHO_DRAWS * 2 ** 16 // (K // BATCH)) + 2 ** 11
MAIN_CONTROL = 256 + RHO_COST
FINISH_CONTROL = 1024 + RHO_COST
CONTROL_TOTAL = Fraction(K * (L * MAIN_CONTROL + TRAIL_MAX * FINISH_CONTROL), BATCH)
R_MAX = L + TRAIL_MAX
assert S_max == K * R_MAX
SKIP_ROUNDS, PLAIN_ROUNDS, ANCHOR_ROUNDS = schedule_counts(R_MAX)
COMMON_TOTAL = Fraction(K * (SKIP_ROUNDS * SKIP_ARM + PLAIN_ROUNDS * PLAIN_ARM + ANCHOR_ROUNDS * ANCHOR_ARM), BATCH) + CONTROL_TOTAL
C_COMMON = COMMON_TOTAL / S_max
PLAIN = Fraction(PLAIN_ARM, BATCH)
ANCHOR = Fraction(ANCHOR_ARM, BATCH)
MAX_LOOK = 64

RATE_TWO = Fraction(1, GAMMA * ADDERS)
RATE_LONG = long_rate(ADDERS)
RATE_FALSE = Fraction(1, 2 ** 63)
assert RATE_LONG == (Fraction(3, 2 * ADDERS) + Fraction(4, ADDERS ** 2)) / (GAMMA * ADDERS)
assert RATE_LONG < Fraction(3, ADDERS ** 2)
W_WIN = 64
P_PRIME = Fraction(W_WIN, GAMMA * ADDERS)

W_CHECK = F_MAX_R + SAME_MAX + ESC2_MAX + PDP_MAX + OVERHEAD
EXIT_EXTRA = 13
W_ENTRY = F_MAX_R + ESC65_MAX + MAX_LOOK * (F_MAX_R + PDP_MAX + SAME_MAX + OVERHEAD) + OVERHEAD + EXIT_EXTRA
KM = 970000
RATE_MERGED = Fraction(T * T, 2 * N) * TRAIL_MAX / S_max

P_FIELD = 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1
ADDER_KEY_BITS, ADDER_KEY_FLOOR = 231, 2 ** 192
Q_SMALL = Fraction(2 ** ADDER_KEY_BITS, P_FIELD)
ADDER_TRIALS = 2 ** 40
ADDER_TEST = 2 ** 10
CHAR_SUM = 2 ** 160
Y_LO = (P_FIELD - (P_FIELD >> 7) + 1) // 2
Y_HI = P_FIELD - Y_LO
assert all((abs(2 * y - P_FIELD) <= P_FIELD >> 7) == (Y_LO <= y <= Y_HI) for y in (Y_LO - 1, Y_LO, Y_HI, Y_HI + 1))
U_DEV = max(Fraction(Y_HI, P_FIELD) - Fraction(1, 2), Fraction(1, 2) - Fraction(Y_LO, P_FIELD))
assert U_DEV <= Fraction(1, 2 ** 8)
L_CS = 1 + math.log(P_FIELD)
BOX_ERR = 6 * (math.isqrt(P_FIELD) + 1) * (2 * L_CS + L_CS ** 2)
BOX_MAIN = Fraction((2 ** ADDER_KEY_BITS - ADDER_KEY_FLOOR) * (Y_HI - Y_LO + 1) * (r - 1), P_FIELD ** 2)
M_ACC = 2 ** 224 - 2 ** 185 - CHAR_SUM - 3
assert BOX_ERR < CHAR_SUM // 2 ** 14 and BOX_MAIN - 2 ** 14 * math.ceil(BOX_ERR) >= M_ACC
Q_ACCEPT_LOWER = Fraction(M_ACC, r)
B_LQ8 = (P_FIELD + 2 ** 128) >> 64
assert B_LQ8 < ADDER_KEY_FLOOR and P_FIELD < 2 ** 256 - ADDER_KEY_FLOOR and (2 ** 256 - 1 - ADDER_KEY_FLOOR) * 2 ** 64 + P_FIELD + 2 ** 128 < 2 ** 320
CORR_SITES = [
    ("fd", 5, Q_SMALL), ("fpre", 5, Fraction(1, 4)), ("q", 5, Fraction(1, 4)), ("inv", 5, Fraction(1, 4)),
    ("dy", 5, Fraction(1, 4) + Fraction(1, 2 ** 16)), ("lam", 5, Fraction(1, 4)), ("sq", 5, Fraction(1, 3)), ("x3a", 5, Fraction(1, 2)),
    ("x3b", 5, Q_SMALL), ("t", 5, Fraction(1, 2)), ("y3m", 5, Fraction(1, 4)),
]
CORR_SLOW = [("inv_slow", 11), ("xy_fix", 47), ("y3_tie", 5), ("y3_tie_fix", 5)]
assert Fraction(1, 4) + U_DEV ** 2 <= CORR_SITES[4][2]
Q_MAIN_PAD, Q_SLOW = Fraction(1, 2 ** 10), Fraction(1, 2 ** 24)
Q_DY_PAD = Fraction(0)
NO_PAD = ("fd", "x3b")
C_MAX = sum(w for _, w, _ in CORR_SITES) + sum(w for _, w in CORR_SLOW)
assert C_MAX == 123
CORR_REPR = C_MAX * Fraction(1, 2 ** 31)
CORR = sum(w * (q + (Q_DY_PAD if n in NO_PAD else Q_MAIN_PAD)) for n, w, q in CORR_SITES) + sum(w * Q_SLOW for _, w in CORR_SLOW) + CORR_REPR
TRAIL_RANGE = (C_MAX + math.ceil(CORR)) * TRAIL_MAX
LAMBDA_TR = 2 ** 112
LAMBDA_PARTIAL = 2 * (C_MAX + math.ceil(CORR)) * K * TRAIL_MAX
LAMBDA_RHO = 2 ** 96
LAMBDA_CAP = LAMBDA_TR + LAMBDA_PARTIAL + LAMBDA_RHO
CORR_OVERSHOOT = K * C_MAX
CORR_FIXED = LAMBDA_CAP + CORR_OVERSHOOT

FALSE_DEBIT = 8192   # 66 later exposed steps + 2P-1 <= 4095 persistent steps
def e0_parts(adders):
    return (long_rate(adders) * W_ENTRY, RATE_FALSE * FALSE_DEBIT * (W_CHECK + W_ENTRY), RATE_MERGED * W_ENTRY)

EVENTS0 = [
    ("two_cycle_entries", RATE_TWO, W_ENTRY),
    ("longer_entries", RATE_LONG, W_ENTRY),
    ("false_filter_matches", RATE_FALSE, FALSE_DEBIT * (W_CHECK + W_ENTRY)),
    ("merged_steps", RATE_MERGED, W_ENTRY),
]
ETA_TWO = Fraction(1, 2 ** 16)
MU_TWO = S_max * RATE_TWO
TWO_LIMIT = math.ceil((1 + ETA_TWO) * MU_TWO)
TWO_CHARGE = Fraction(TWO_LIMIT * W_ENTRY, S_max)
CHERNOFF_EXPONENT = MU_TWO * ETA_TWO ** 2 / (2 + ETA_TWO)
assert CHERNOFF_EXPONENT > 100
DELTA_TWO = Fraction(1, 2 ** 100)
E0 = sum(rate * cost for _, rate, cost in EVENTS0[1:])
assert E0 == sum(e0_parts(ADDERS))
LG_DEN = M_ACC - GAMMA * TRAIL_MAX ** 2 // 2
def closure_term(s):
    return (Fraction(ADDERS**2*math.comb(ADDERS,s), ADDERS-s)
            * Fraction(1,GAMMA**(s-1))*Fraction(s,ADDERS)**max(63,2*s-2))
RHO_FAR = Fraction(93,100)
assert Fraction(68,25)**2/8*(1+Fraction(3,ADDERS))**2 < RHO_FAR
Q_FAR_EXACT = sum(closure_term(s) for s in range(1,65))
Q_FAR_BOUND = Q_FAR_EXACT + closure_term(64)*RHO_FAR/(1-RHO_FAR) + Fraction(1,2**600)
assert TRAIL_MAX.bit_length()+ADDERS.bit_length()+601 < ADDERS
assert Q_FAR_BOUND < Fraction(1,2**528) and S_max < 2**128
LOG2_Q_LONG = -528
LOG2_DELTA_LONG = -400
assert GAMMA * TRAIL_MAX ** 2 <= M_ACC // 2
RATE_LOOP_STEP = GAMMA * TRAIL_MAX * (Fraction(2, M_ACC) + Fraction(1, LG_DEN))
RATE_LOOP_TOTAL = RATE_LOOP_STEP + Fraction(2.0 ** LOG2_Q_LONG)
A_LOOP = RATE_LOOP_TOTAL * TRAIL_MAX
E_RARE = E0
E_LOOP = Fraction(0)
RESERVE = TWO_CHARGE + KM * E_RARE
PER_STEP = C_COMMON + CORR + RESERVE
step_units = PER_STEP / C
SCALAR_MUL = 2 * 384 * F_MAX
# Distinguished = the key's top DP_BITS bits zero: key < 2^(256-DP_BITS) = 2^196, so
# theta = 2^196/p >= 2^-60 and the record bound is re-derived at it (not at 2^-60).
P_FIELD_DECL = 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1
assert 2 ** 196 < P_FIELD_DECL
THETA = Fraction(2 ** 196, P_FIELD_DECL)
assert THETA >= Fraction(1, 2 ** THETA_BITS) and THETA < Fraction(1, 2 ** THETA_BITS) * (1 + Fraction(1, 2 ** 31))
DP_MAX = 2 * -(-(S_max * 2 ** 196) // P_FIELD_DECL) + 2 * K
assert DP_MAX >= 2 * S_max * THETA + 2 * K
K_TRAILS = DP_MAX + 2 * K
DY_EXCESS = 5 * (Fraction(2 * K_TRAILS, S_max) + (2 * ANCHOR_P + 2) * (RATE_TWO + RATE_LONG))
assert DY_EXCESS < 5 * Q_MAIN_PAD
LOG2_DELTA_TR = -(2 * LAMBDA_TR ** 2 / (K_TRAILS * TRAIL_RANGE ** 2)) / math.log(2)
assert LOG2_DELTA_TR < -1000
N_ROUNDS = L + TRAIL_MAX
N_BATCH_ROUNDS = N_ROUNDS
LOG2_DELTA_RHO = -(2 * LAMBDA_RHO ** 2 / (N_ROUNDS * (2 * K * C_MAX) ** 2)) / math.log(2)
assert LOG2_DELTA_RHO < -1000
DELTA_TABLE = 4 * math.exp(-2 * ADDERS * 2.0 ** -20)
P_FIELD_RHO = 2 ** 256 - 2 ** 224 + 2 ** 192 + 2 ** 96 - 1
LOG2_DELTA_DRAW = math.log2(N_BATCH_ROUNDS) + RHO_DRAWS * math.log2(float(Fraction(2 ** 256 - P_FIELD_RHO + 1, 2 ** 256)))
DELTA_CORR = 2.0 ** LOG2_DELTA_TR + 2.0 ** LOG2_DELTA_RHO + DELTA_TABLE + 2.0 ** LOG2_DELTA_DRAW
REPLAY_OPS = 2 ** 13
assert (MAX_LOOK + 1) * 70 <= REPLAY_OPS
assert 33 * (MAX_LOOK + 4) < REPLAY_OPS and 33 * (F_MAX + W_ENTRY) <= REPLAY_OPS * (F_MAX + OVERHEAD)
D_CAP = 4
DRAW_COST = 2 ** 16
PROBE_CAP = 2 ** 10
PROBE_COST = 2 ** 11
DP_INDEX_COST = 2 ** 11
DP_COUNTER_COST = 8
DP_OVERFLOW_COST = 8
DP_COST = (PROBE_CAP * PROBE_COST + DP_INDEX_COST + DP_COUNTER_COST + DP_OVERFLOW_COST)
START = D_CAP * (2 * SCALAR_MUL + 2 * D_CAP * DRAW_COST) + 2 ** 14

led = Ledger(C)
led.add("walk_steps", S_max * step_units)
led.add_instructions("correction_meter_slack", CORR_FIXED)
led.add_instructions("adder_table", ADDERS * ADDER_TRIALS * (START + ADDER_TEST), preprocessing=True)
led.add_instructions("initial_walk_starts", K * START, preprocessing=True)
led.add_instructions("trail_restarts", (DP_MAX + K) * START)
led.add_instructions("trail_end_exits", K_TRAILS * EXIT_EXTRA)
led.add_instructions("distinguished_points", DP_MAX * DP_COST)
led.add_instructions("collision_replay", 2 * TRAIL_MAX * REPLAY_OPS * (F_MAX + OVERHEAD))
led.add_instructions("solve_verify_output", 2 * SCALAR_MUL + 2 ** 24)
led.add_instructions("input_validation", 2 ** 24)

B_RARE = S_max * RESERVE
AUX_CAP = B_RARE // 10
N_DP_EXPOSURES = S_max + AUX_CAP + 4*(ADDERS*ADDER_TRIALS + K_TRAILS) + K + 1
MU_DP = THETA * N_DP_EXPOSURES
assert MU_DP < Fraction(101,100)*S_max*THETA
assert DP_MAX-1 > Fraction(3,2)*MU_DP
assert S_max*THETA/10 > 60
DELTA_DP = Fraction(1,2**60)
DELTA_EXCEPTIONAL = Fraction(1,2**95)
assert Fraction(N_DP_EXPOSURES*ADDERS,N) < DELTA_EXCEPTIONAL
assert DELTA_DP + Fraction(1,r) + DELTA_EXCEPTIONAL + Fraction(1,2**80) < delta
J_AB_ALL = 67*TRAIL_MAX + 1
K_AB_ALL = CAP // 89600
P_AB_ENTRY = P_PRIME*(1+Fraction(1,2**31))
D_AB_ALL = CAP - 2*ANCHOR_P*K_AB_ALL - 1
assert Fraction(68,25)*J_AB_ALL*P_AB_ENTRY/K_AB_ALL < Fraction(98,100)
assert K_TRAILS < 2**69 and K_AB_ALL > 8450
INDEX_COUPLING_EXPOSURES = N_DP_EXPOSURES + K_TRAILS*J_AB_ALL
DELTA_INDEX = Fraction(1,2**80)
assert INDEX_COUPLING_EXPOSURES*Fraction(ADDERS,P_FIELD_DECL) < DELTA_INDEX
DELTA_AB_FRESH_UPPER = float(K_TRAILS)*math.exp(-float(Fraction(D_AB_ALL,2**THETA_BITS)))
assert DELTA_AB_FRESH_UPPER < 1.25e-6
DELTA_ANY_ABANDON = Fraction(13,10**7)
assert DELTA_AB_FRESH_UPPER + 2.0**-100 < float(DELTA_ANY_ABANDON)
DELTA_LONG = 2.0 ** LOG2_DELTA_LONG
DELTA_LOOP = float(K_TRAILS * Fraction(TRAIL_MAX ** 2, 2) * GAMMA * (Fraction(2, M_ACC) + Fraction(1, LG_DEN))
                   + Fraction(ADDERS ** 2 * GAMMA, M_ACC))
DELTA_CUT = DELTA_LOOP + DELTA_LONG
TAU_EXACT = Fraction(T * T, 2 * N)
N_SCALARS = 2 * (ADDERS * ADDER_TRIALS + K_TRAILS)
R_REJECT = Fraction(2 ** 256 - r, 2 ** 256)
LOG2_CAP_DRAWS = math.log2(N_SCALARS) + D_CAP * math.log2(float(R_REJECT))
LOG2_CAP_STARTS = math.log2(ADDERS * ADDER_TRIALS + K_TRAILS) - D_CAP * math.log2(r)
LOG2_ADDER_EXHAUST = math.log2(ADDERS) + ADDER_TRIALS * math.log1p(-float(Q_ACCEPT_LOWER)) / math.log(2)
assert LOG2_ADDER_EXHAUST < -300
DELTA_CAP = 2.0 ** LOG2_CAP_DRAWS + 2.0 ** LOG2_CAP_STARTS + 2.0 ** LOG2_ADDER_EXHAUST
C_SLOTS = 2 * DP_MAX
C_TAIL = math.log(2) - 0.5
TAIL_PROBE = math.exp(-C_TAIL * PROBE_CAP) / (1 - math.exp(-C_TAIL))
DELTA_PROBE = float(DP_MAX) * float(C_SLOTS) * TAIL_PROBE
# The DP key is uniform on [0, 2^196) (the test reads the top 60 bits alone), so key >> 32 is uniform on
# exactly M = 2^164 values: the slot input needs no p-granularity refinement.
INDEX_SHIFT = ADDERS.bit_length() - 1
assert INDEX_SHIFT == 32 and 2 ** (256 - 60) >> INDEX_SHIFT == 2 ** 164
M_HI = 2 ** 164
TV_HASH = Fraction(M_HI % C_SLOTS * (C_SLOTS - M_HI % C_SLOTS), C_SLOTS * M_HI)
DELTA_HASH = Fraction(2 * DP_MAX * (C_SLOTS + 1), M_HI)
assert DELTA_HASH >= DP_MAX * TV_HASH
BANK_SLOTS = 2 ** 40
N_BANKS = -(-C_SLOTS // BANK_SLOTS)
INIT_PER_SLOT = 8 + 2
INIT_PER_BANK = 64
TABLE_INIT = INIT_PER_SLOT * C_SLOTS + INIT_PER_BANK * N_BANKS
led.add_instructions("record_table_init", TABLE_INIT, preprocessing=True)

tau = Fraction(T * T, 2 * N)
p_model = 1 - math.exp(-float(tau * (1 - eps)))
DELTA_ADDER_CLASS = float(Fraction(S_max * ADDERS, N))
assert DELTA_ADDER_CLASS < 2 ** -95
p_claim = p_model - float(delta) - DELTA_ADDER_CLASS - float(DELTA_ANY_ABANDON) - DELTA_CUT - DELTA_CAP - DELTA_PROBE - float(DELTA_HASH) - 1 / KM - float(DELTA_TWO) - DELTA_CORR
assert p_claim - 0.39 >= 2.5e-7
assert PLAIN_FIXED == 605065449
assert P_FIELD * 2 ** 64 + P_FIELD + 2 ** 128 < 2 ** 320
X64 = 2 ** 64
assert (2 ** 256 - 1 + (X64 - 1) * P_FIELD) // X64 < P_FIELD + 2 ** 192
assert (P_FIELD - 1 + (X64 - 1) * P_FIELD) // X64 < P_FIELD
for R_BND, B_BND, A3 in ((P_FIELD - 1, P_FIELD - 1, X64 - 1), (2 ** 256 - 1, P_FIELD + 2 ** 192, P_FIELD >> 192)):
    S6 = (X64 - 1) * (2 * B_BND + R_BND + P_FIELD) + X64 * A3 * R_BND
    assert S6 < 2 ** 384 and (S6 + X64 * (X64 - 1) * P_FIELD) // X64 ** 2 < 2 ** 256 + P_FIELD
assert (X64 - 1) ** 2 + X64 * 2 * (X64 - 1) ** 2 < 2 ** 193 + 2 ** 128 and Fraction(2 ** 65 + 2, 2 ** 128) < Fraction(1, 2 ** 62)
assert log2_up(led.total()) == Fraction(126728821, 1000000)
assert log2_up(led.preprocessing_total()) == Fraction(8996119, 100000)
LOG2_HOEFFDING_TR = -(2 * (LAMBDA_TR - 2 ** 109) ** 2 / (K_TRAILS * TRAIL_RANGE ** 2)) / math.log(2)
assert LOG2_HOEFFDING_TR < -1384
assert W_ENTRY == 2021790 and C_MAX == 123 and DELTA_CORR <= 5.07e-14 and log2(N_ROUNDS) < 95.492
assert log2(Fraction(CORR_FIXED, C)) < 102.059 and CORR_REPR * S_max < 2 ** 104 and N_ROUNDS * (K // BATCH) < 2 ** 108
assert Fraction(1, 4) + Fraction(1, 2 ** 10) + Fraction(1, 2 ** 24) + 6 * Fraction(1, 2 ** 31) <= Fraction(250977, 10 ** 6)
assert Fraction(PLAIN_FIXED, BATCH) <= Fraction(578538829, 10 ** 6) and PLAIN <= Fraction(578539806, 10 ** 6)
assert C_COMMON < Fraction(575092878, 10 ** 6) and CORR <= Fraction(14210693, 10 ** 6) and PER_STEP <= Fraction(589305610, 10 ** 6)
# typo guards, not premises (the binding floor is the p_claim assert above)
assert RESERVE < Fraction(3, 1000) and TWO_CHARGE < Fraction(3, 10000) and E0 < Fraction(1, 10 ** 8)
pass
pass
pass
assert DY_EXCESS < Fraction(1938, 10 ** 8) and math.log2(BOX_ERR) < 145.5598 and math.log2(2 ** 224 - BOX_MAIN) < 185.0001
assert Y_HI - Y_LO + 1 < 2 ** 249 and math.log2(2 ** 231) - math.log2(P_FIELD) < -24.9999 and 192 - math.log2(P_FIELD) < -63.99

peak_bytes = DP_MAX * 128 * 2 + K * 512 + 2 * ADDERS * 2 * 64 + 2 ** 24
assert peak_bytes < 2 ** 77

events = [
    {"name": n, "rate": float(rate), "worst_case_instructions": int(cost),
     "expected_instructions_per_walk_step": float(rate * cost)}
    for n, rate, cost in EVENTS0
]

out = led.report()
out.update({
    "L": str(L),
    "log2_L": log2(L),
    "log2_T": log2(T),
    "log2_S_max": log2(S_max),
    "step_units": str(step_units),
    "escape_depth_bound": str(K_ESC),
    "ord_2_lower_bound": {"smooth_part_of_r_minus_1": str(S_ORDER),
                          "ord_r_of_2_to_the_cofactor": str(ORD2_LOWER),
                          "requirement": "ord_r(2) > |Gamma| * escape_depth_bound",
                          "holds": ORD2_LOWER > GAMMA * K_ESC},
    "success": {"tau": float(tau), "eps": float(eps), "delta": float(delta),
                "model": p_model, "lower_bound": p_claim, "delta_any_abandon": float(DELTA_ANY_ABANDON),
                "all_abandon_dp_tail_upper": DELTA_AB_FRESH_UPPER, "entry_bound_p_prime": float(P_PRIME),
                "two_cycle_entry_rate": float(RATE_TWO), "longer_entry_rate": float(RATE_LONG),
                "delta_index": float(DELTA_INDEX), "index_coupling_exposures": str(INDEX_COUPLING_EXPOSURES),
                "delta_dp": float(DELTA_DP), "delta_exceptional": float(DELTA_EXCEPTIONAL),
                "delta_degenerate": float(Fraction(1,r)), "dp_exposures_max": str(N_DP_EXPOSURES),
                "delta_loop": DELTA_LOOP, "delta_cut": DELTA_CUT, "delta_long": DELTA_LONG,
                "delta_cap": DELTA_CAP, "delta_probe": DELTA_PROBE, "delta_hash": float(DELTA_HASH), "hash_tv_per_record": float(TV_HASH), "d_cap": D_CAP, "probe_cap": PROBE_CAP, "start_instructions": START, "expected_merges_tau": float(TAU_EXACT),
                "log2_q_long": LOG2_Q_LONG, "log2_delta_long": LOG2_DELTA_LONG},
    "small_key_adders": {"adder_key_bits": ADDER_KEY_BITS, "adder_key_floor_log2": 192, "q_small_per_site": float(Q_SMALL),
                         "log2_q_small": math.log2(float(Q_SMALL)), "trials_cap_per_adder": ADDER_TRIALS,
                         "trial_instructions": START + ADDER_TEST, "char_sum_bound_log2": math.log2(CHAR_SUM),
                         "log2_accept_lower": math.log2(float(Q_ACCEPT_LOWER)),
                         "log2_all_adders_exhaust": LOG2_ADDER_EXHAUST, "y_window": [str(Y_LO), str(Y_HI)], "y_dev_max": float(U_DEV), "box_main_log2": log2(BOX_MAIN), "box_error_log2": math.log2(BOX_ERR), "dy_unfolded_excess_per_step": float(DY_EXCESS), "m_acc_lower_log2": math.log2(M_ACC), "relation_factor_r_over_m": float(Fraction(r, M_ACC)), "delta_adder_class": DELTA_ADDER_CLASS,
                         "adder_table_instructions": str(ADDERS * ADDER_TRIALS * (START + ADDER_TEST))},
    "memory_bytes_log2": str(log2_up(peak_bytes)),
    "dp_max_log2": log2(DP_MAX),
    "record_table_init": {"slots": str(C_SLOTS), "banks": str(N_BANKS), "per_slot": INIT_PER_SLOT,
                          "per_bank": INIT_PER_BANK, "instructions": str(TABLE_INIT)},
    "table_operation": {"probe_cap": PROBE_CAP, "probe_cost": PROBE_COST,
                        "worst_case_probes": PROBE_CAP * PROBE_COST,
                        "index_reduction": DP_INDEX_COST, "counter_update": DP_COUNTER_COST,
                        "overflow_test": DP_OVERFLOW_COST, "dp_cost": DP_COST,
                        "operations": str(DP_MAX), "charged": str(DP_MAX * DP_COST)},
    "k_trails_log2": log2(K_TRAILS),
    "split_concentration": {"eta": str(ETA_TWO), "two_count_limit": str(TWO_LIMIT),
        "two_charge_per_step": float(TWO_CHARGE), "chernoff_exponent_log2": log2(CHERNOFF_EXPONENT),
        "delta_two": str(DELTA_TWO), "residual_markov_factor": KM,
        "residual_expected_per_step": float(E_RARE), "residual_reserve_per_step": float(KM*E_RARE),
        "e0_at_a16_a20": [float(sum(e0_parts(2 ** b))) for b in (16, 20)]},
    "corrections": {"sites": [{"name": n, "extra": w, "q0": str(q), "charged_q": float(q + (Q_DY_PAD if n in NO_PAD else Q_MAIN_PAD))} for n, w, q in CORR_SITES],
                    "slow_sites": [{"name": n, "extra": w, "charged_q": float(Q_SLOW)} for n, w in CORR_SLOW],
                    "corr_per_step": float(CORR), "corr_repr": float(CORR_REPR), "c_max": C_MAX, "trail_range_log2": log2(TRAIL_RANGE),
                    "lambda_tr_log2": log2(LAMBDA_TR), "lambda_partial_log2": log2(LAMBDA_PARTIAL), "lambda_rho_log2": log2(LAMBDA_RHO),
                    "overshoot": str(CORR_OVERSHOOT), "k_trails_log2": log2(K_TRAILS), "rounds_log2": log2(N_ROUNDS),
                    "log2_delta_tr": LOG2_DELTA_TR, "log2_delta_tr_after_debits": LOG2_HOEFFDING_TR, "log2_delta_rho": LOG2_DELTA_RHO, "log2_delta_draw": LOG2_DELTA_DRAW,
                    "delta_table": DELTA_TABLE, "delta_corr": DELTA_CORR, "rho_cost_per_round": RHO_COST,
                    "batch": BATCH, "plain_fixed_round": PLAIN_FIXED, "anchor_extra_round": ANCHOR_EXTRA,
                    "method_allowance_per_round": METHOD_ALLOWANCE},
    "common_path": {"plain_arm_round": str(PLAIN_ARM), "skip_arm_round": str(SKIP_ARM),
                    "anchor_arm_round": str(ANCHOR_ARM), "skip_rounds_per_slot": str(SKIP_ROUNDS),
                    "plain_rounds_per_slot": str(PLAIN_ROUNDS), "anchor_rounds_per_slot": str(ANCHOR_ROUNDS),
                    "rounds_per_slot": str(R_MAX), "c_common_exact": str(C_COMMON),
                    "anchor_schedule": {"period": ANCHOR_P, "untested_tail": ANCHOR_U,
                        "tested_offsets": [ANCHOR_U + 1, ANCHOR_P], "tested_width": ANCHOR_T,
                        "max_period_mod_l_1_to_64": ANCHOR_R_MAX,
                        "coverage_holds": ANCHOR_R_MAX < ANCHOR_T,
                        "cap_periods": str(CAP_PERIODS), "trail_max_steps": str(TRAIL_MAX),
                        "persistent_states_charged": PERSIST}}, 
    "meter": {"c_common": float(C_COMMON), "e_rare_per_walk_step": float(E_RARE),
              "events": events,
              "forced_maxima": {"f": F_MAX, "add": ADD_MAX, "same": SAME_MAX, "pdp": PDP_MAX, "esc2": ESC2_MAX,
                                "esc65": ESC65_MAX, "overhead": OVERHEAD},
              "reserve_per_walk_step": float(RESERVE), "charge_per_walk_step": float(PER_STEP),
              "abort_probability_bound": "1/%d + 2^-100, on G" % KM, "walk_steps_total": str(S_max)},
})
pass
pass
print(json.dumps(out, indent=2))
