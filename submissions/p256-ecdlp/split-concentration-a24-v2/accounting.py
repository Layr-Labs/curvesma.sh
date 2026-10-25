"""Bundled experiment: entrypoint accounting_src.py with vendored cscore.exact.
cscore.exact sha256 ad40a116ac9b6b9d430ae9baa9b18c4c218a6fc996c486d18c4815d650d9e638
Derived from the declared predecessor; split concentration accounting revision."""

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

# ---- entrypoint: accounting_src.py ----

"""Exact accounting for the claim in proof.md (stdlib only; prints JSON).

The P-256 track has no reference operation for curve arithmetic: every scalar AArch64
instruction of the attack, including all field arithmetic, is charged at 1/C with the
borrowed normalisation C=1000. Integers and rationals are exact; logarithms are rounded
up to 6 decimals.

Revised parameters: A=2^24 adders, tau at the success floor, and the residual Markov
factor re-optimised at that A. Every rate below is an exact function of ADDERS.
"""
import json
import math
from fractions import Fraction

from cscore.exact import Ledger, iroot_ceil, log2, log2_up

r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
C = 1000
ADDERS = 2 ** 24             # r-adding walk; index = low 24 bits of the class key M(x)
GAMMA = 2                    # |Gamma|: the class automorphisms are the signs only
K = 2 ** 32                  # walks (2^16 processors x 2^16-walk batches)
THETA_BITS = 60              # distinguished-point probability 2^-60
CAP_PERIODS = 2 ** THETA_BITS
# A trail ends at its first distinguished point. A trail that crosses more than CAP_PERIODS anchor
# boundaries (64 steps apart) without one is abandoned, after n = 64 (CAP_PERIODS + 1) steps.
CAP = 64 * 2 ** THETA_BITS
TRAIL_MAX = CAP + 64         # every trail has at most TRAIL_MAX steps

# Classes {W, -W}: N nonzero classes, plus one for the identity (never reached by a walk).
N = (r - 1) // 2 + 1

def long_rate(adders):  # E1(c): longer-entry rate E[w']/(|Gamma| A) <= 3/(4A^2) + 2/A^3
    return Fraction(3, 4 * adders ** 2) + Fraction(2, adders ** 3)


# Success model (proof.md 2.6): tau = T^2/(2N), 1 - exp(-tau (1 - eps)), less delta, the
# abandonment and within-trail terms, and the meter's 1/KM. eps is the H1 margin of 2.6, with
# 192 wasted states per 3..64 window charged at the longer-entry rate.
eps = Fraction(4, ADDERS) + 2 * Fraction(3, ADDERS) + 2 * 192 * long_rate(ADDERS)
delta = Fraction(1, 2 ** 20)
# Smallest tau on the 1e-6 grid whose recomputed conditional bound keeps a 1.0e-5 margin above 0.39
# (proof.md, "Complete cost and success ledger"); the ledger asserts p_claim >= 0.390005.
TAU_TARGET = Fraction(494321, 1000000)
T_min = iroot_ceil(math.ceil(2 * N * TAU_TARGET), 2)
L = -(-T_min // K)                 # steps per walk counted towards the collision budget
T = K * L
S_max = K * (L + TRAIL_MAX)        # worst case: after L steps each walk finishes its current trail

# D1 (proof.md 2.1): an escape point's S-coefficient differs from every earlier point of its trail
# when ord_r(2) > |Gamma| K_ESC, where K_ESC bounds the escapes of a trail. Each step escapes at
# most once, so K_ESC = TRAIL_MAX. Let S be the smooth part of r - 1, factored below, and C2 = (r-1)/S.
# ord_r(2^C2) divides ord_r(2), and it is computed exactly from the factorisation of S alone.
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
ORD2_LOWER = order_dividing(U2, S_ORDER, S_FACTORS)     # = ord_r(2^C2), which divides ord_r(2)
assert pow(U2, ORD2_LOWER, r) == 1
assert all(pow(U2, ORD2_LOWER // p, r) != 1 for p in S_FACTORS if ORD2_LOWER % p == 0)
assert ORD2_LOWER > GAMMA * K_ESC

# Forced-path worst cases (Ir per call; measurement/forced/summary.json). Every unbatched operation of a
# rare path is charged its worst case, so no per-call average enters the charge. P-256 has no rotation
# search and no exact-k fallback: measurement/step.c canonicalises by the parity of M(y) alone.
F_MAX = 41190        # f_single: adder selection and the addition
ADD_MAX = 41185      # add_single
SAME_MAX = 32        # same_class
PDP_MAX = 10         # point_is_dp
ESC2_MAX = 41262     # escape of a 2-cycle (least-key member, canon, doubling)
ESC65_MAX = 43089    # escape of a 65-cycle
OVERHEAD = 1024      # control overhead per rare event and per look-ahead step (measured: far below this)
UNBATCHED = F_MAX    # unit for scalar multiplications and replays: an upper bound for an addition or doubling
# Common path charged by the anchor schedule (proof.md). evidence/measurement/seed1.json holds Callgrind totals of
# measurement/step.c over one batch of 2^16 walks: the plain round Ir(2)-Ir(1), equal on both input sets, and the
# anchor round Ir(64)-Ir(63), taken as measured with its two-cycle escapes left in (they are also metered below),
# the larger input set. step.c sets an anchor only in rounds t with t % 64 == 0 of the global batch round counter.
# The 2^24-adder build runs 2 more instructions per round (alignment nops before the first loop), added here.
BATCH = 2 ** 16
PLAIN = Fraction(52337843 + 2, BATCH)          # per walk step of a plain batch round
ANCHOR = Fraction(53305368 + 2, BATCH)         # per walk step of an anchor batch round, rare events included
WIDE_CONTROL = Fraction(1, 64)                 # 1024 per batch round for 128-bit round control, every step
R_MAX = L + TRAIL_MAX                          # batch rounds of a slot: L main-phase steps, then its current trail
assert S_max == K * R_MAX
ANCHOR_ROUNDS = -(-R_MAX // 64)                # at most ceil(R_MAX/64) rounds t <= R_MAX have t % 64 == 0
COMMON_TOTAL = S_max * (PLAIN + WIDE_CONTROL) + K * ANCHOR_ROUNDS * (ANCHOR - PLAIN)
C_COMMON = COMMON_TOTAL / S_max                # per processed step, exact
MAX_LOOK = 64       # look-ahead steps per anchor check (MAX_CYCLE)

# Entry rates (proof.md 2.3, Lemma E1). An entry is a step whose source is on a fresh class and whose successor
# is a revisit. Under M1 and L-ori the source's index and orientation are uniform given the history, so:
#  - the 2-cycle entry (the successor is the predecessor's class) has conditional probability exactly 1/(|Gamma| A);
#  - a longer entry is one of the admissible windows of length 3..64 ending at the source; each window gives one
#    target (index, orientation), so its conditional probability is w'/(|Gamma| A) with w' <= 63 on G. In expectation
#    over the history, only odd windows can be unit vectors, and E[w'] <= 3/(2A) + 4/A^2 (P_3 = 3/(2A) - 3/(4A^2),
#    P_5 <= 15/(4A^2)), so the longer-entry rate is at most 3/(4A^2) + 2/A^3 = 2.664536e-15 at A = 2^24.
#    The ledger charges that derived bound, not the inherited 4x larger allowance c/A^2 with c = 3.
# The false 63-bit filter match is 2^-63 per processed step (rejected by the exact check).
RATE_TWO = Fraction(1, GAMMA * ADDERS)                  # 2-cycle entry, exact (L-ori)
RATE_LONG = long_rate(ADDERS)                           # derived longer-entry bound, charged
RATE_FALSE = Fraction(1, 2 ** 63)                       # false 63-bit filter match, per processed step
assert RATE_LONG == (Fraction(3, 2 * ADDERS) + Fraction(4, ADDERS ** 2)) / (GAMMA * ADDERS)
assert RATE_LONG < Fraction(3, ADDERS ** 2)
W_WIN = 64                                              # admissible windows of length 1..64 (on G, a revisit is within 64 steps)
P_PRIME = Fraction(W_WIN, GAMMA * ADDERS)               # conditional entry bound on G: 64 / (|Gamma| A)

# Metered rare events (proof.md 2.5, Lemma M'). W_CHECK: one 2-cycle check (f_single, same_class, escape,
# point_is_dp, overhead). W_ENTRY: one anchor check with its look-ahead steps (f_single, escape of up to 65
# members, up to 64 look-ahead steps, each f_single, point_is_dp and same_class, with overhead).
W_CHECK = F_MAX + SAME_MAX + ESC2_MAX + PDP_MAX + OVERHEAD
W_ENTRY = F_MAX + ESC65_MAX + MAX_LOOK * (F_MAX + PDP_MAX + SAME_MAX + OVERHEAD) + OVERHEAD
# Markov factor for stopped residual work only: a free analysis parameter. The reserve charges KM
# times the residual expectation and the failure allowance is 1/KM. At this tau, KM=2^18 is the
# cheapest power of two that keeps the declared 1.0e-5 margin above 0.39.
KM = 2 ** 18
# Merged steps (proof.md section 2.8): at most tau TRAIL_MAX of them in expectation, each at the anchor-check cost.
RATE_MERGED = Fraction(T * T, 2 * N) * TRAIL_MAX / S_max


def e0_parts(adders):  # per-step residual expectation parts, proof.md 2.5; only the first depends on A
    return (long_rate(adders) * W_ENTRY, RATE_FALSE * 256 * (W_CHECK + W_ENTRY), RATE_MERGED * W_ENTRY)


EVENTS0 = [
    ("two_cycle_entries", RATE_TWO, W_CHECK),
    ("longer_entries", RATE_LONG, W_ENTRY),
    ("false_filter_matches", RATE_FALSE, 256 * (W_CHECK + W_ENTRY)),   # a false first-filter match and a false anchor match may both run
    ("merged_steps", RATE_MERGED, W_ENTRY),
]
# Separate fresh-source two-cycle entries analytically; runtime retains one total meter.
ETA_TWO = Fraction(1, 2 ** 16)
MU_TWO = S_max * RATE_TWO
TWO_LIMIT = math.ceil((1 + ETA_TWO) * MU_TWO)
TWO_CHARGE = Fraction(TWO_LIMIT * W_CHECK, S_max)
# Chernoff exponent mu*eta^2/(2+eta) is far larger than 100 ln(2).
CHERNOFF_EXPONENT = MU_TWO * ETA_TWO ** 2 / (2 + ETA_TWO)
assert CHERNOFF_EXPONENT > 100  # exp(-100) < 2^-100
DELTA_TWO = Fraction(1, 2 ** 100)
E0 = sum(rate * cost for _, rate, cost in EVENTS0[1:])
assert E0 == sum(e0_parts(ADDERS))   # the quoted component list is the charged expectation
# Within-trail coincidences are excluded by G and charged through DELTA_CUT.
# Compute their probability only; the current reserve does not price a trail
# by a global average work allowance.
LG_DEN = r - GAMMA * TRAIL_MAX ** 2 // 2
# Long closures (proof.md Lemma F): Q_FAR <= sum_{s<A} A C(A,s) |Gamma|^(1-s) (s/A)^(R_s) A/(A-s) + TRAIL_MAX A |Gamma|^(1-A).
def _log2_binom(n, k):
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)) / math.log(2)
def log2_q_far(dmin=64, s_cut=4096):
    """log2 of the Lemma F closure sum; the s < s_cut truncation is an upper bound (proof.md, Complete ledger)."""
    terms = []
    for s in range(1, min(ADDERS, s_cut)):
        R = max(dmin - 1, 2 * s - 2)
        terms.append(math.log2(ADDERS) + _log2_binom(ADDERS, s) + (1 - s) * math.log2(GAMMA)
                     + R * math.log2(s / ADDERS) + math.log2(ADDERS) - math.log2(ADDERS - s))
    rho = (math.e ** 2 / 8) * (1 + 3.0 / ADDERS) ** 2
    assert rho < 0.93 and s_cut >= 64
    terms.append(terms[-1] + math.log2(1.0 / (1.0 - rho)))   # geometric tail bound
    m = max(terms)
    main = m + math.log2(sum(2.0 ** (t - m) for t in terms)) + 2.0 ** -40   # 2^-40 float slack
    full = math.log2(TRAIL_MAX) + math.log2(ADDERS) + (1 - ADDERS) * math.log2(GAMMA)
    hi, lo = max(main, full), min(main, full)
    return hi + math.log2(1 + 2.0 ** (lo - hi))
LOG2_Q_LONG = log2_q_far()                              # per fresh start
assert LOG2_Q_LONG < -200
LOG2_DELTA_LONG = math.log2(S_max) + LOG2_Q_LONG        # union bound over the S_max starts                       # union bound over the S_max fresh starts
assert GAMMA * TRAIL_MAX ** 2 <= r // 2
RATE_LOOP_STEP = GAMMA * TRAIL_MAX * (Fraction(2, r) + Fraction(1, LG_DEN))
RATE_LOOP_TOTAL = RATE_LOOP_STEP + Fraction(2.0 ** LOG2_Q_LONG)   # within-trail coincidences and long closures
A_LOOP = RATE_LOOP_TOTAL * TRAIL_MAX  # retained diagnostic, not used in the reserve
E_RARE = E0  # residual stopped at first G violation; DELTA_CUT pays all bad-G outcomes   # expected metered rare work per walk step, upper bound
E_LOOP = Fraction(0)  # not a work estimate outside G; those outcomes are failed
RESERVE = TWO_CHARGE + KM * E_RARE
PER_STEP = C_COMMON + RESERVE                           # charged instructions per walk step
step_units = PER_STEP / C
SCALAR_MUL = 2 * 384 * F_MAX   # double-and-add over a 256-bit scalar, generous: 768 ops
DP_MAX = 2 * S_max // 2 ** THETA_BITS + 2 * K
K_TRAILS = DP_MAX + 2 * K      # trails started: K initial, at most DP_MAX records, at most K abandonments
REPLAY_OPS = 2 ** 13               # unbatched operations per live step, replay included: 65 replay steps x 70 operations <= 2^13 (proof.md 1.3)
assert (MAX_LOOK + 1) * 70 <= REPLAY_OPS   # one walk step and one look-ahead of MAX_CYCLE = 64 steps per live step
# Sampling caps (proof.md section 1.3, steps 2 and 3): each scalar is drawn by at most D_CAP rejection-sampling draws, and
# each start point (and each adder point) is attempted at most D_CAP times; the next draw or attempt aborts the run.
D_CAP = 4
DRAW_COST = 2 ** 16                          # allowance per draw: OS randomness and the rejection test
# Distinguished-point table (proof.md section 8.2): open addressing, C = 2 DP_MAX slots, probe cap L = PROBE_CAP.
PROBE_CAP = 2 ** 10                          # L: a lookup or insert that examines L slots aborts the run
PROBE_COST = 2 ** 11                         # allowance per probe: slot fetch (memory access) plus key compare and control
# Once-per-operation work, charged inside the same table operation (proof.md 8.2) rather than assumed to sit in the
# probe allowance: the index is reduced from the canonical key by one two-limb long division, the record counter is
# incremented and written back as two limbs, and the incremented counter is tested for overflow against DP_MAX.
DP_INDEX_COST = 2 ** 11                      # index reduction from the canonical key
DP_COUNTER_COST = 8                          # record counter: two-limb load, add, add-with-carry and store
DP_OVERFLOW_COST = 8                         # overflow test: high-limb compare against DP_MAX and branch
DP_COST = (PROBE_CAP * PROBE_COST + DP_INDEX_COST + DP_COUNTER_COST + DP_OVERFLOW_COST)
# one table operation at its forced worst case: L probes plus the itemised once-per-operation work above
# One start point: D_CAP attempts, each with two scalars of at most D_CAP draws and two scalar multiplications.
START = D_CAP * (2 * SCALAR_MUL + 2 * D_CAP * DRAW_COST) + 2 ** 14

led = Ledger(C)
led.add("walk_steps", S_max * step_units)
# the adder table (independent random points), each with two scalar multiplications; the escape is a doubling, no table
led.add_instructions("adder_table", ADDERS * START, preprocessing=True)
led.add_instructions("initial_walk_starts", K * START, preprocessing=True)
# a new trail after each distinguished point (at most DP_MAX) and after each abandonment (at most K)
led.add_instructions("trail_restarts", (DP_MAX + K) * START)
led.add_instructions("distinguished_points", DP_MAX * DP_COST)
# one key match is resolved: both trails are replayed from their coins with coefficients
led.add_instructions("collision_replay", 2 * TRAIL_MAX * REPLAY_OPS * (F_MAX + OVERHEAD))
led.add_instructions("solve_verify_output", 2 * SCALAR_MUL + 2 ** 24)
led.add_instructions("input_validation", 2 ** 24)  # on-curve check; cofactor 1

# Lemma T'' (proof.md 2.4): abandonment, factorial moments over entries. The processed steps before abandonment
# are at most n = TRAIL_MAX. The conditional entry bound P_PRIME holds on G for every history, so the chain rule gives
# E[binom(E_n, k)] <= binom(n, k) P_PRIME^k. With c = 368: e c P_PRIME <= 1/2 and 128/c < 0.35.
C_AB = 368
n_ab = 64 * (CAP_PERIODS + 1)
assert n_ab == TRAIL_MAX
assert math.e * C_AB * float(P_PRIME) <= 0.5
k_ab = n_ab // C_AB
log_entries = k_ab * math.log(math.e * n_ab * float(P_PRIME) / k_ab)   # ln (e n p' / k)^k
n_min_ab = CAP + 1                                                      # abandonment comes at step s + CAP, s in [1, 64]
log_fresh = (n_min_ab - 128 * k_ab) * math.log1p(-2.0 ** -THETA_BITS)     # ln (1 - theta)^(CAP + 1 - 128 k)
P_ABANDON = math.exp(log_entries) + math.exp(log_fresh)
# Abandonment on the good event G (proof.md 2.4): X counts abandoned trails, X <= K_TRAILS, and
# Pr[X > K] <= Pr[not G] + E[X 1_G] / K. The first term is DELTA_CUT below; the second is DELTA_ABANDON.
E_ABANDON = float(K_TRAILS) * P_ABANDON                  # expected abandoned trails on G (at most)
DELTA_ABANDON = E_ABANDON / float(K)                     # Markov on G: P[more than K abandonments, G] <= E / K
DELTA_LONG = 2.0 ** LOG2_DELTA_LONG                     # long closures of length 65 or more (proof.md 2.3)
# Within-trail coincidences, degenerate or not, over all trails started (proof.md 2.2, DELTA_LOOP, Lemma G denominator).
DELTA_LOOP = float(K_TRAILS * Fraction(TRAIL_MAX ** 2, 2) * GAMMA * (Fraction(2, r) + Fraction(1, LG_DEN))
                   + Fraction(ADDERS ** 2 * GAMMA, r))
# Cut-off and within-trail event (proof.md 2.2, DELTA_CUT): not G, which is DELTA_LOOP plus the long-closure term.
DELTA_CUT = DELTA_LOOP + DELTA_LONG
# Merges that do not resolve (proof.md section 2.8, Lemma M): 2 tau (1-theta)^(CAP/2), tau the birthday count.
TAU_EXACT = Fraction(T * T, 2 * N)
DELTA_MERGE = 2 * float(TAU_EXACT) * math.exp(-float(CAP) / 2 ** THETA_BITS / 2)
# Resource aborts (proof.md section 8.2): a sampling cap or the probe cap aborts the run.
N_SCALARS = 2 * (ADDERS + K_TRAILS)          # scalars drawn: two per adder and two per start point
R_REJECT = Fraction(2 ** 256 - r, 2 ** 256)  # a draw of randbelow(r) from [0, 2^256) is rejected with this probability
LOG2_CAP_DRAWS = math.log2(N_SCALARS) + D_CAP * math.log2(float(R_REJECT))
LOG2_CAP_STARTS = math.log2(ADDERS + K_TRAILS) - D_CAP * math.log2(r)
DELTA_CAP = 2.0 ** LOG2_CAP_DRAWS + 2.0 ** LOG2_CAP_STARTS
# Probe cap: C = 2 DP_MAX slots, at most DP_MAX records (load <= 1/2); Chernoff per window, union over operations and windows.
C_SLOTS = 2 * DP_MAX
C_TAIL = math.log(2) - 0.5
TAIL_PROBE = math.exp(-C_TAIL * PROBE_CAP) / (1 - math.exp(-C_TAIL))
DELTA_PROBE = float(DP_MAX) * float(C_SLOTS) * TAIL_PROBE

tau = Fraction(T * T, 2 * N)
p_model = 1 - math.exp(-float(tau * (1 - eps)))
p_claim = p_model - float(delta) - DELTA_ABANDON - DELTA_CUT - DELTA_MERGE - DELTA_CAP - DELTA_PROBE - 1 / KM - float(DELTA_TWO)
assert p_claim >= 0.390005
assert p_claim - 0.39 >= 1.0e-5          # the declared margin, not only the floor
assert log2_up(led.total()) == Fraction(127167780, 1000000)   # the declared score, claims[0].score
assert log2(led.preprocessing_total()) < 50

# Total-work meter (proof section 1.3 step 8): every walk step is charged PER_STEP, every batch
# inversion is inside C_COMMON, and rare work is metered at worst case. The meter aborts when the
# rare-work total exceeds KM times its expectation, which happens with probability at most 1/KM.

peak_bytes = DP_MAX * 128 * 2 + K * 512 + ADDERS * 2 * 64 + 2 ** 24
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
                "model": p_model, "lower_bound": p_claim, "delta_abandon": DELTA_ABANDON,
                "p_abandon_per_trail": P_ABANDON, "entry_bound_p_prime": float(P_PRIME),
                "two_cycle_entry_rate": float(RATE_TWO), "longer_entry_rate": float(RATE_LONG),
                "delta_loop": DELTA_LOOP, "delta_cut": DELTA_CUT, "delta_long": DELTA_LONG,
                "delta_merge": DELTA_MERGE, "delta_cap": DELTA_CAP, "delta_probe": DELTA_PROBE, "d_cap": D_CAP, "probe_cap": PROBE_CAP, "start_instructions": START, "expected_merges_tau": float(TAU_EXACT),
                "log2_q_long": LOG2_Q_LONG, "log2_delta_long": LOG2_DELTA_LONG},
    "memory_bytes_log2": str(log2_up(peak_bytes)),
    "dp_max_log2": log2(DP_MAX),
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
    "common_path": {"plain_per_walk_step": str(PLAIN), "anchor_per_walk_step": str(ANCHOR),
                    "rounds_per_slot": str(R_MAX), "anchor_rounds_per_slot": str(ANCHOR_ROUNDS),
                    "c_common_exact": str(C_COMMON)},
    "meter": {"c_common": float(C_COMMON), "e_rare_per_walk_step": float(E_RARE),
              "events": events,
              "forced_maxima": {"f": F_MAX, "add": ADD_MAX, "same": SAME_MAX, "pdp": PDP_MAX, "esc2": ESC2_MAX,
                                "esc65": ESC65_MAX, "overhead": OVERHEAD},
              "reserve_per_walk_step": float(RESERVE), "charge_per_walk_step": float(PER_STEP),
              "abort_probability_bound": "1/%d + 2^-100, on G" % KM, "walk_steps_total": str(S_max)},
})
print(json.dumps(out, indent=2))
