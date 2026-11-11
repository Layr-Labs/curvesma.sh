"""Planned-phase ledger, finite-width control primitives and G1 output for the BN254 rho.

Part P  exact ledger of the planned search phase (fixed step count S_plan, bounded
        detection tail), the report store, the hard instruction meter and the
        premature-truncation margin. Normal no-hit completion and premature truncation
        are separate events and are reported separately.
Part Q  the control primitives on 64-bit limbs exactly as itemised in machine.md
        section 4a (128-bit planned/tail counters, 192-bit instruction meter with a
        check before every priced block, 64-bit report index), executed limb by limb
        and checked against Python integers at carry, boundary and overflow values,
        including a metered run that must stop before exceeding floor(cap).
Part O  output map x -> ([x]P0, Q0): left-to-right double-and-add on BN254 G1 in
        Jacobian coordinates with counted Fp operations. Gives the analytical worst-case
        bound (all 253 doublings and 253 additions) and, separately, the sample maximum
        over edge and random scalars, each checked against an independent affine
        reference; identity case z = 1 -> (infinity, Q0).
Standard library only. Optional JSON on stdin: {"seed": int, "cases": int}.
This script does not run the rho itself; walk evidence is experiment.py's (see proof.md 7).
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

try:
    cfg = json.loads(sys.stdin.read() or "{}")
except ValueError:
    cfg = {}
SEED = int(cfg.get("seed", 1))
CASES = int(cfg.get("cases", 64))
rng = random.Random(SEED)

base = Path("/snapshot/files")
if not (base / "target-profile.json").is_file():
    base = Path(__file__).resolve().parent
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
P0 = (int(profile["generator_P0"]["x"], 16), int(profile["generator_P0"]["y"], 16))

M64 = (1 << 64) - 1
C_UNIT = 20000
FP_MUL, FP_ADD, FP_SUB, ROUTINE = 162, 26, 27, 32      # leaf listing + call site; non-leaf call
MUL_STEP, SQ_STEP = 18191, 6965                        # machine.md section 4 (COMMON 300 kept)
COMMON_ITEMISED = 153                                  # machine.md section 4a (was 129)
SQ_PROB = Fraction(1, 2)
R_TABLE = 1024
M_WALKS = 2 ** 20
THETA_BITS = 70                                        # DP density 2^-70
TAIL_PER_WALK = 20 * 2 ** THETA_BITS                   # 20 / theta
P_MODEL = Decimal("0.43")
CLAIM = Fraction(2493, 20)                             # 124.65
REPORT_INSTR = 10000                                   # per report, machine.md 5 (incl. meter)
CHAIN_STEP_INSTR = 32                                  # per chain-node visit (incl. meter)
CHAIN_CAP_LOG2 = 70
NODE_BYTES, HEAD_BYTES = 96, 8
START_UNITS = 400
RESTART_CAP_LOG2 = 48
FALSE_CAND_CAP_LOG2 = 16
OUTPUT_UNITS = 101
NONWALK_LOG2 = Decimal("64.5")


def iroot(n, k):
    lo, hi = 0, 1 << (n.bit_length() // k + 1)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid - 1
    return lo


# floor(20000 * 2^(2493/20)) = floor((20000^20 * 2^2493)^(1/20)), exact integer
CAP_INSTR = iroot(C_UNIT ** 20 * 2 ** 2493, 20)


def part_P():
    getcontext().prec = 80
    ln2 = Decimal(2).ln()
    lg = lambda v: Decimal(v).ln() / ln2
    a = Decimal(MUL_STEP) / C_UNIT
    b = Decimal(SQ_STEP) / C_UNIT
    s = Decimal(SQ_PROB.numerator) / SQ_PROB.denominator
    c = (1 - s) * a + s * b
    N = Decimal(r - 1) / 12
    q2 = s * s + (1 - s) ** 2 / R_TABLE
    Np = N / (1 - q2)
    S = (2 * Np * (1 / (1 - P_MODEL)).ln()).sqrt()
    eps = Decimal(5) / (12 * R_TABLE) + Decimal(250) / (R_TABLE * R_TABLE)
    K = int((S * (1 + eps) / M_WALKS).to_integral_value(rounding="ROUND_CEILING"))
    S_plan = K * M_WALKS
    E = M_WALKS * TAIL_PER_WALK
    n = S_plan + E
    theta = Decimal(2) ** -THETA_BITS
    # report store (machine.md 5): 64-bit node indices, one address space below 2^63 bytes
    mu_D = Decimal(n) * theta
    D_cap = int((Decimal("1.1") * S_plan * theta).to_integral_value(rounding="ROUND_CEILING")) + 2 ** 25
    buckets = 1 << (D_cap * 3 // 2).bit_length()          # power of two >= 1.5 D_cap: load <= 2/3
    store_bytes = D_cap * NODE_BYTES + buckets * HEAD_BYTES
    exp_chain_visits = Decimal(D_cap) * D_cap / (2 * buckets)   # expected same-bucket pairs
    chain_cap_markov_log2 = lg(exp_chain_visits) - CHAIN_CAP_LOG2
    delta = Decimal(D_cap) / mu_D - 1
    chernoff_log2 = -(mu_D * delta * delta / (2 + delta)) / ln2
    false_cand_exp = Decimal(D_cap) * D_cap / Decimal(2) ** 129
    exp_restarts = mu_D * Decimal(-20).exp() + M_WALKS   # trails without DP within 20/theta, plus starts
    # non-walk work in units, each row a cap (every row is enforced or deterministic)
    rows = {
        "setup_table_and_starts": Decimal(12288 + M_WALKS) * START_UNITS,
        "bucket_head_zeroing": Decimal(buckets) * Decimal("0.75") / C_UNIT,  # stp 16 B + loop 0.25/head
        "reports": Decimal(D_cap) * REPORT_INSTR / C_UNIT,
        "chain_visits_cap": Decimal(2) ** CHAIN_CAP_LOG2 * CHAIN_STEP_INSTR / C_UNIT,
        "restarts_cap": Decimal(2) ** RESTART_CAP_LOG2 * START_UNITS,
        "false_candidates_cap": Decimal(2) ** FALSE_CAND_CAP_LOG2 * 600,
        "phase_switch": Decimal(M_WALKS) * 50 / C_UNIT,
        "output": Decimal(OUTPUT_UNITS),
    }
    nonwalk = sum(rows.values())
    NONWALK = Decimal(2) ** NONWALK_LOG2
    preprocessing = rows["setup_table_and_starts"] + rows["bucket_head_zeroing"]
    W_nom = c * n + NONWALK
    Wcap_units = Decimal(CAP_INSTR) / C_UNIT
    m_max = (Wcap_units - NONWALK - n * b - a) / (a - b)  # one step of early stop for the pre-check
    tol = m_max / n - (1 - s)
    L_max = 2 * n * tol * tol / (20 * ln2)
    p_model = 1 - (-(S * S) / (2 * Np)).exp()
    losses = {
        "premature_truncation_allowance": Decimal("0.001"),
        "escape_overhead_above_2eps_and_lost_steps": Decimal("0.0005"),
        "report_cap_chernoff_under_H2": Decimal(2) ** chernoff_log2 if chernoff_log2 > -1000 else Decimal(0),
        "chain_visit_cap_markov": Decimal(2) ** chain_cap_markov_log2,
        "restart_cap_markov": exp_restarts / Decimal(2) ** RESTART_CAP_LOG2,
        "false_candidate_cap_markov": false_cand_exp / Decimal(2) ** FALSE_CAND_CAP_LOG2,
        "detection_tail_miss": Decimal(-20).exp(),
        "first_collision_inside_one_walk": Decimal(1) / M_WALKS,
        "beta_tie_and_key_tie": Decimal(2) ** -240,
        "H1_H2_inaccuracy_reserve": Decimal("0.01"),
    }
    success_lb = P_MODEL - sum(losses.values())
    checks = {
        "K_bits_le_128": K.bit_length() <= 128,
        "tail_bits_le_128": TAIL_PER_WALK.bit_length() <= 128,
        "report_cap_bits_le_64": D_cap.bit_length() <= 64,
        "store_below_2^63_bytes": store_bytes < 2 ** 63,
        "cap_instr_bits_le_192": CAP_INSTR.bit_length() <= 192,
        "nonwalk_le_bucket": nonwalk <= NONWALK,
        "itemised_common_le_300": COMMON_ITEMISED <= 300,
        "claim_ok": Wcap_units >= W_nom,
        "success_ok": success_lb >= Decimal("0.39"),
    }
    return {
        "multiply_step_units": float(a), "squaring_step_units": float(b), "mean_step_units": float(c),
        "log2_N_eff": float(lg(Np)), "log2_S_useful_steps": float(lg(S)), "eps": float(eps),
        "K_planned_steps_per_walk": str(K), "K_bits": K.bit_length(),
        "log2_S_plan": float(lg(S_plan)), "tail_per_walk_bits": TAIL_PER_WALK.bit_length(),
        "log2_tail_E": float(lg(E)),
        "cap_instructions_floor": str(CAP_INSTR), "cap_bits": CAP_INSTR.bit_length(),
        "log2_W_nominal_full_run": float(lg(W_nom)),
        "claim_minus_nominal_bits": float(Decimal(CLAIM.numerator) / CLAIM.denominator - lg(W_nom)),
        "tolerated_multiply_fraction_excess_points": float(tol * 100),
        "log2_largest_dependence_block_for_2^-20": float(lg(L_max)),
        "report_cap": D_cap, "log2_report_cap": float(lg(D_cap)), "log2_expected_reports": float(lg(mu_D)),
        "buckets_log2": buckets.bit_length() - 1, "log2_store_bytes": float(lg(store_bytes)),
        "log2_expected_chain_visits": float(lg(exp_chain_visits)),
        "nonwalk_rows_log2_units": {k: float(lg(v)) for k, v in rows.items()},
        "log2_nonwalk_total": float(lg(nonwalk)), "log2_preprocessing_units": float(lg(preprocessing)),
        "model_success_at_S": float(p_model),
        "losses": {k: float(v) for k, v in losses.items()},
        "success_lower_bound": float(success_lb),
        "checks": checks,
    }


# ------------------------------------------------------------------ Part Q: limb primitives
# Each primitive returns (result, executed instruction count). Sequences (machine.md 4a):
#   inc128:   ldp lo,hi ; adds lo,lo,#1 ; adc hi,hi,xzr ; stp            -> 4
#   ge128:    ldp k_lo,k_hi ; cmp lo,k_lo ; sbcs xzr,hi,k_hi ; b.hs      -> 4
#   meter pre-check: ldp a0,a1 ; ldr a2 ; adds t0,a0,#PMAX ; adcs t1,a1,xzr ; adc t2,a2,xzr ;
#             ldp c0,c1 ; ldr c2 ; cmp c0,t0 ; sbcs xzr,c1,t1 ; sbcs xzr,c2,t2 ; b.lo abort -> 11
#   meter charge: mov x,#price (2: movz/movk) ; adds ; adcs ; adc ; stp ; str            -> 7

def inc128(lo, hi):
    lo2 = (lo + 1) & M64
    carry = 1 if lo2 < lo or (lo == M64) else 0
    hi2 = (hi + carry) & M64
    return (lo2, hi2), 4


def ge128(lo, hi, klo, khi):
    # cmp lo,klo sets C = (lo >= klo); sbcs hi - khi - !C sets C = no borrow
    c = 1 if lo >= klo else 0
    c2 = 1 if hi - khi - (1 - c) >= 0 else 0
    return bool(c2), 4


def add192(a, v):
    a0, a1, a2 = a
    t0 = a0 + v
    c = t0 >> 64
    t1 = a1 + c
    c = t1 >> 64
    t2 = (a2 + c) & M64
    return (t0 & M64, t1 & M64, t2)


def le192(x, y):                    # x <= y via cmp/sbcs chain: y - x without borrow
    c = 1 if y[0] >= x[0] else 0
    d1 = y[1] - x[1] - (1 - c)
    c = 1 if d1 >= 0 else 0
    d2 = y[2] - x[2] - (1 - c)
    return d2 >= 0


def limbs192(v):
    return (v & M64, (v >> 64) & M64, v >> 128)


def val192(t):
    return t[0] | (t[1] << 64) | (t[2] << 128)


def meter_block(acc, price, pmax, cap):
    """Check before the block (with its maximal price), then charge the actual price."""
    if not le192(add192(acc, pmax), cap):
        return acc, False, 11
    return add192(acc, price), True, 11 + 7


def part_Q(K, tail):
    fails = []
    # 128-bit counters across the 2^64 boundary, at K - 1, K, and the tail limit
    for v in [0, M64 - 1, M64, 2 ** 64, K - 2, K - 1, K, tail - 1, tail, 2 ** 128 - 2]:
        (lo, hi), n_inc = inc128(v & M64, v >> 64)
        if (lo | (hi << 64)) != (v + 1) % 2 ** 128 or n_inc != 4:
            fails.append(("inc128", v))
        for lim in (K, tail):
            ge, n_ge = ge128(v & M64, v >> 64, lim & M64, lim >> 64)
            if ge != (v >= lim) or n_ge != 4:
                fails.append(("ge128", v, lim))
    # a walk counter run from K - 3 stops after exactly 3 increments at K
    v, steps = K - 3, 0
    while not ge128(v & M64, v >> 64, K & M64, K >> 64)[0]:
        lo, hi = inc128(v & M64, v >> 64)[0]
        v, steps = lo | (hi << 64), steps + 1
    if v != K or steps != 3:
        fails.append(("phase_stop", v, steps))
    # 192-bit meter: carries at 2^64 and 2^128, and stop-before-exceed at floor(cap)
    cap = limbs192(CAP_INSTR)
    for start in [2 ** 64 - 5, 2 ** 128 - 5, CAP_INSTR - 3 * MUL_STEP - 10]:
        acc, n_ok, ok = limbs192(start), 0, True
        while ok:
            price = MUL_STEP if (n_ok % 2 == 0) else SQ_STEP
            acc2, ok, _ = meter_block(acc, price, MUL_STEP, cap)
            if ok:
                if val192(acc2) != val192(acc) + price:
                    fails.append(("add192", start))
                acc, n_ok = acc2, n_ok + 1
            if n_ok > 10 ** 6:
                break
        if val192(acc) > CAP_INSTR:
            fails.append(("meter_exceeded", start))
        if start < CAP_INSTR - 10 ** 6 and n_ok <= 10 ** 6:
            fails.append(("meter_false_abort", start))
        if start == CAP_INSTR - 3 * MUL_STEP - 10 and val192(acc) + MUL_STEP <= CAP_INSTR:
            fails.append(("meter_late_abort", start))
    return {"failures": fails, "cap_limbs_hex": [hex(x) for x in cap],
            "instructions": {"inc128": 4, "ge128": 4, "meter_precheck": 11, "meter_charge": 7}}


# ------------------------------------------------------------------ Part O: G1 output
class Cnt:
    def __init__(self):
        self.mul = self.add = self.sub = self.calls = 0

    def m(self, x, y):
        self.mul += 1
        return x * y % p

    def a(self, x, y):
        self.add += 1
        return (x + y) % p

    def s(self, x, y):
        self.sub += 1
        return (x - y) % p

    def instr(self):
        return self.mul * FP_MUL + self.add * FP_ADD + self.sub * FP_SUB + self.calls * ROUTINE


def jac_dbl(k, X, Y, Z):           # dbl-2009-l, a = 0: 2M + 5S
    k.calls += 1
    A = k.m(X, X); B = k.m(Y, Y); Cc = k.m(B, B)
    t = k.a(X, B); t = k.m(t, t); t = k.s(t, A); t = k.s(t, Cc); D = k.a(t, t)
    E = k.a(A, A); E = k.a(E, A); F = k.m(E, E)
    X3 = k.s(F, k.a(D, D))
    C8 = k.a(Cc, Cc); C8 = k.a(C8, C8); C8 = k.a(C8, C8)
    Y3 = k.s(k.m(E, k.s(D, X3)), C8)
    Z3 = k.m(Y, Z); Z3 = k.a(Z3, Z3)
    return X3, Y3, Z3


def jac_madd(k, X1, Y1, Z1, x2, y2):   # madd-2007-bl: 7M + 4S
    k.calls += 1
    Z1Z1 = k.m(Z1, Z1); U2 = k.m(x2, Z1Z1); S2 = k.m(y2, k.m(Z1, Z1Z1))
    H = k.s(U2, X1); HH = k.m(H, H); I = k.a(HH, HH); I = k.a(I, I); J = k.m(H, I)
    rr = k.s(S2, Y1); rr = k.a(rr, rr); V = k.m(X1, I)
    X3 = k.s(k.s(k.m(rr, rr), J), k.a(V, V))
    t = k.m(Y1, J); t = k.a(t, t)
    Y3 = k.s(k.m(rr, k.s(V, X3)), t)
    Z3 = k.a(Z1, H); Z3 = k.s(k.s(k.m(Z3, Z3), Z1Z1), HH)
    assert H != 0, "exceptional addition"
    return X3, Y3, Z3


def g1_mul_counted(x):
    k = Cnt()
    X, Y, Z = P0[0], P0[1], 1
    for bit in bin(x)[3:]:
        X, Y, Z = jac_dbl(k, X, Y, Z)
        if bit == "1":
            X, Y, Z = jac_madd(k, X, Y, Z, P0[0], P0[1])
    k.calls += 1                    # inversion routine
    zi = Z
    for bit in bin(p - 2)[3:]:
        zi = k.m(zi, zi)
        if bit == "1":
            zi = k.m(zi, Z)
    zi2 = k.m(zi, zi); zi3 = k.m(zi2, zi)
    return (k.m(X, zi2), k.m(Y, zi3)), k


def aff_add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    if P[0] == Q[0] and (P[1] + Q[1]) % p == 0:
        return None
    if P == Q:
        l = 3 * P[0] * P[0] * pow(2 * P[1], -1, p) % p
    else:
        l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, p) % p
    x3 = (l * l - P[0] - Q[0]) % p
    return (x3, (l * (P[0] - x3) - P[1]) % p)


def aff_mul(x, P):
    R_ = None
    while x:
        if x & 1:
            R_ = aff_add(R_, P)
        P = aff_add(P, P)
        x >>= 1
    return R_


def part_O(cases):
    # analytical worst case: per-call counts are input independent; at most 253 doublings
    # and 253 additions since x < r < 2^254
    kd, km = Cnt(), Cnt()
    jac_dbl(kd, P0[0], P0[1], 1)
    X, Y, Z = jac_dbl(Cnt(), P0[0], P0[1], 1)
    jac_madd(km, X, Y, Z, P0[0], P0[1])
    _, k1 = g1_mul_counted(1)                    # inversion + affine only
    analytic = 253 * (kd.instr() + km.instr()) + k1.instr()
    loop = 254 * 12                              # bit extraction, test, branch per bit
    serialise = 2 * (FP_MUL + 64 * 8)            # leave Montgomery form, 64 hex digits per coordinate
    charged = OUTPUT_UNITS * C_UNIT
    edge = [1, 2, 3, r - 1, r - 2, (r - 1) // 2, 2 ** 253 - 1, 3 * 2 ** 252 - 1, (1 << 253) | ((1 << 252) - 1)]
    xs = [x for x in edge if 1 <= x < r] + [rng.randrange(1, r) for _ in range(cases)]
    fails, sample_max, arg = 0, 0, None
    for x in xs:
        P, k = g1_mul_counted(x)
        if P != aff_mul(x, P0) or (P[1] * P[1] - P[0] ** 3 - 3) % p:
            fails += 1
        if k.instr() > analytic:
            fails += 1
        if k.instr() > sample_max:
            sample_max, arg = k.instr(), x
    order_ok = aff_mul(r, P0) is None and aff_mul(r - 1, P0) == (P0[0], (-P0[1]) % p)
    return {"cases": len(xs), "failures": fails, "P0_order_r": order_ok,
            "identity_output": [None, "Q0"],
            "per_doubling_instr": kd.instr(), "per_addition_instr": km.instr(),
            "inversion_and_affine_instr": k1.instr(),
            "analytic_worst_instr": analytic, "loop_instr": loop, "serialise_instr": serialise,
            "analytic_total_instr": analytic + loop + serialise,
            "sample_max_instr": sample_max, "sample_argmax_hex": hex(arg),
            "charged_instr": charged, "charged_units": OUTPUT_UNITS,
            "charge_covers_analytic": charged >= analytic + loop + serialise}


P = part_P()
Q = part_Q(int(P["K_planned_steps_per_walk"]), TAIL_PER_WALK)
O = part_O(CASES)
out = {"seed": SEED, "partP_planned_phase": P, "partQ_control_primitives": Q, "partO_g1_output": O}
out["all_ok"] = (all(P["checks"].values()) and not Q["failures"] and O["failures"] == 0
                 and O["P0_order_r"] and O["charge_covers_analytic"])
print(json.dumps(out, indent=1, default=str))
sys.exit(0 if out["all_ok"] else 1)
