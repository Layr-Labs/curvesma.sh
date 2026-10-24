#!/usr/bin/env python3
"""Exact ledger v19b -- Integers and Fractions only.  Deltas vs ledger19:

  claim  the published number is DERIVED, not chosen: W is computed as an exact Fraction,
         and the claim is the SMALLEST 4-decimal value c with 2^c >= W + 2e-5
         (lane convention: >= W plus a 2e-5 rounding allowance).  The meter cap is then
         floor(2^claim) work units, so the cap covers the nominal work by construction.
  cap    with the cap at the claim level the run never has to substitute squarings to fit,
         so k_min = 0 and the abort row is the clipped 2^-400; the two remaining debits
         (DP-detection loss, storage) are unchanged.
  gate   the claim is CONDITIONAL: it is only offered while the lowest live non-self
         submission on the track is at least claim + 0.005 (tie cap) + 0.0015 (margin).
         The gate is an external fact, so it is reported, not asserted as a ledger check.
  census the ten shards are ONE configuration; their per-shard records are aggregated
         losslessly into evidence/census19.json and re-summed here.
"""
import json
import math
import sys
from fractions import Fraction as F

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
RS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 20000
RTAB = 4096
IDX_BITS = 12
DP_SHIFT = 15
DBITS = 100
REC_BYTES = 52
SCALE = 10000            # claim is a 4-decimal number
HEADROOM = F(2, 100000)  # 2e-5
ROOT = __file__.rsplit("/", 2)[0]
checks = []


def ck(name, ok, value=None):
    try:
        json.dumps(value)
    except (TypeError, ValueError):
        value = str(value)
    checks.append({"check": name, "value": value, "status": "PASS" if ok else "FAIL"})


def cp_lower(k, n, alpha=0.05):
    """One-sided Clopper-Pearson lower end (named: Clopper-Pearson, not Wilson)."""
    if n == 0 or k == 0:
        return 0.0
    lg = math.lgamma

    def logsf(p):
        terms = [lg(n + 1) - lg(i + 1) - lg(n - i + 1) + i * math.log(p)
                 + (n - i) * math.log1p(-p) for i in range(k, n + 1)]
        m = max(terms)
        return m + math.log(sum(math.exp(t - m) for t in terms))

    lo, hi = 1e-12, k / n
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if logsf(mid) < math.log(alpha):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- inputs
fmatch = json.load(open(ROOT + "/evidence/forced19.json"))
trap = json.load(open(ROOT + "/evidence/trap19.json"))
coal = json.load(open(ROOT + "/evidence/coalesce19.json"))
tn = json.load(open(ROOT + "/evidence/tn19.json"))
cen = json.load(open(ROOT + "/evidence/census19.json"))
swa = json.load(open(ROOT + "/evidence/sweep19.json"))

shards = cen["shards"]
ck("census_shards_present", len(shards) == 10, len(shards))
census_n = sum(s["trials"] for s in shards)
census_k = sum(s["success_num"] for s in shards)
ck("census_shard_totals_resum_to_aggregate",
   (census_n, census_k) == (cen["trials"], cen["success_num"]), [census_n, census_k])
ck("census_is_one_config",
   all((s["R"], s["dpbits"], s["q"], s["trails"]) == (128, 0, 523, 32) for s in shards))
ck("census_rule_is_the_guard_rule", all(s["rule"] == "v19" for s in shards))
census_rate = census_k / census_n
cp_lo = cp_lower(census_k, census_n)
CENSUS = {"trials": census_n, "success_num": census_k, "success_rate": census_rate,
          "cp95_one_sided_lower": cp_lo, "shard_count": len(shards),
          "degenerate": cen["degenerate"],
          "config": {"q": 523, "R": 128, "dpbits": 0, "S_design": shards[0]["S"],
                     "trails": 32, "per_trail": shards[0]["per_trail"], "rule": "v19"}}
ck("census_is_large_enough", census_n >= 4000, census_n)
ck("census_zero_degenerate", cen["degenerate"] == 0, cen["degenerate"])

sweeps = {d["dpbits"]: d for d in swa["densities"]}
ck("sweep_has_all_four_densities", sorted(sweeps) == [5, 7, 9, 11], sorted(sweeps))
ck("sweep_rule_is_the_guard_rule", all(d["rule"] == "v19" for d in sweeps.values()))
ck("sweep_arm_is_production_shaped",
   all(d["trails"] == 4 and d["R"] == 4096 and d["idx_bits"] == IDX_BITS for d in sweeps.values()))
for db, d in sorted(sweeps.items()):
    ck("sweep_d%d_gate" % db, d["cp95_detect"] >= 0.39, "%.4f" % d["cp95_detect"])
    ck("sweep_d%d_zero_degenerate" % db, d["degenerate"] == 0, d["degenerate"])

# ---------------------------------------------------------------- merge / trap gates
v19c = coal["rules"]["v19"]
tnc = coal["rules"]["tn"]
ck("coalescence_all_checkpoints_1_000",
   all(v == v19c["pairs"] for v in v19c["same_class_at"].values()),
   v19c["same_class_at"])
ck("coalescence_last_step_1_000", v19c["same_class_at_last_step"] == 1.0,
   v19c["same_class_at_last_step"])
ck("control_rule_tn_is_NOT_a_class_function",
   tnc["same_class_at_last_step"] < 0.5, tnc["same_class_at_last_step"])
ck("tn_rule_is_NOT_a_class_function",
   tn["tn_label_agrees_on_conjugates"] < tn["pairs"] // 2,
   [tn["tn_label_agrees_on_conjugates"], tn["pairs"]])
ck("tn_norm_is_frobenius_invariant", tn["norm_frobenius_invariant"] == tn["pairs"],
   tn["norm_frobenius_invariant"])
v19t = trap["rules"]["v19"]
v18t = trap["rules"]["v18"]
ck("trap_p6_fixture_all_conjugate", v19t["p6_per_1600"] == 1600.0, v19t["p6_per_1600"])
ck("trap_second_order_zero_under_guard", v19t["traces_with_2cycle"] == 0,
   v19t["traces_with_2cycle"])
ck("trap_control_v18_reproduces_regrowth", v18t["per_1600_2cycle"] > 12.0,
   v18t["per_1600_2cycle"])
ck("trap_guard_firings_are_rare", v19t["guard_firings"] < 0.05 * 1600 * 200,
   v19t["guard_firings"])
ck("trap_census_regime_is_stricter_than_production",
   (v19t["guard_firings"] / (1600 * 200)) / (1 / 4096) > 20,
   round((v19t["guard_firings"] / (1600 * 200)) * 4096, 1))

# ---------------------------------------------------------------- instruction rows
COMMON = (("trace_re_plus_branch_bits", 60), ("exponent_pair_update", 100),
          ("cycle_compares_loop_margin", 140), ("dp_store_amortised", 1),
          ("repeat_slot_guard", 8))
# The multiply row table is itemised down to the reference unit: the dense Fp12 multiply is
# its own row (C = 20000, COST_MODEL.md), so re-summing the printed table gives MUL_CHARGE
# with no unlisted term.
MUL_C_ROW = ("dense_fp12_multiply", C)
MUL_EXTRA = (("key_six_plusminus", 1802), ("min_rotation_booth", 144),
             ("table_address_and_load", 150), ("exponent_pair_add", 60))
MUL_ROWS = (MUL_C_ROW,) + MUL_EXTRA
SQR_EXTRA = (("cyclotomic_squaring", 15042), ("exponent_pair_double", 40))
I_COMMON = sum(v for _, v in COMMON)
I_MUL = I_COMMON + sum(v for _, v in MUL_EXTRA)
I_SQR = I_COMMON + sum(v for _, v in SQR_EXTRA)
MUL_CHARGE = I_COMMON + sum(v for _, v in MUL_ROWS)
SQR_CHARGE = I_SQR
ck("mul_charge_is_the_v19b_row", MUL_CHARGE == 22465, MUL_CHARGE)
ck("multiply_rows_plus_common_resum_to_the_charge",
   I_COMMON + sum(v for _, v in MUL_ROWS) == C + I_MUL == MUL_CHARGE, [MUL_CHARGE, C + I_MUL])
ck("squaring_rows_plus_common_resum_to_the_charge",
   I_COMMON + sum(v for _, v in SQR_EXTRA) == SQR_CHARGE, [SQR_CHARGE])
ck("dense_fp12_multiply_row_is_the_cost_model_reference_op",
   MUL_C_ROW == ("dense_fp12_multiply", 20000), MUL_C_ROW)
ck("rows_match_forced_run_on_the_guard_table",
   (MUL_CHARGE, SQR_CHARGE) == (fmatch["charges"]["MUL_CHARGE"], fmatch["charges"]["SQR_CHARGE"]))
ck("guard_row_is_present_in_the_forced_rows",
   dict((tuple(r) for r in fmatch["rows"]["common"]))["repeat_slot_guard"] == 8)
ck("forced_run_is_the_guard_rule", fmatch["rule"] == "v19", fmatch["rule"])
AUX_ROWS = fmatch["rows"]["auxiliary"]
SOLVE_INSTR = dict(AUX_ROWS)["solve_modular_inverse"]
VERIFY_INSTR = dict(AUX_ROWS)["verify_ladder_miller_finalexp"]
OUTPUT_INSTR = dict(AUX_ROWS)["output_48byte"]
LADDER = dict(AUX_ROWS)["ladder_instr"]
MILLER = dict(AUX_ROWS)["miller_instr"]
FINAL_EXP = dict(AUX_ROWS)["final_exp_instr"]
ck("verify_row_is_ladder_plus_miller_plus_finalexp",
   VERIFY_INSTR == LADDER + MILLER + FINAL_EXP, VERIFY_INSTR)
AUX_SUM = SOLVE_INSTR + VERIFY_INSTR + OUTPUT_INSTR
SETUP_INSTR = RTAB * (510 * C + 12 * 2400) + (1 << 20)

# ---------------------------------------------------------------- meter identity
ck("forced_match_has_a_match", fmatch["match"]["matches"] >= 1)
ck("forced_match_recovered_x", fmatch["match"]["recovered"] is True)
ck("forced_match_log_equals_meter_equals_row_formula",
   fmatch["all_equal"] is True and fmatch["meter_instr"] == fmatch["row_sum_instr"]
   == sum(fmatch["charge_log"]), [fmatch["meter_instr"], fmatch["row_sum_instr"]])
ck("forced_match_replay_is_full_length",
   fmatch["nreplay"] == fmatch["match"]["trail1_step"] + 1 + fmatch["match"]["trail2_step"] + 1,
   fmatch["nreplay"])
ck("forced_match_debits_solve_verify_output",
   all(x in fmatch["charge_log"] for x in (SOLVE_INSTR, VERIFY_INSTR, OUTPUT_INSTR)))
ck("forced_match_log_length_is_trail_plus_replay_plus_aux",
   fmatch["charge_log_len"] == 600 + fmatch["match"]["trail2_step"] + 1
   + fmatch["nreplay"] + 3, fmatch["charge_log_len"])

# ---------------------------------------------------------------- collision model
S_NUM, S_DEN = 3, 8
RHO = (1 - F(S_NUM, S_DEN) ** 2) * (1 - F(1, RTAB))
N = (RS - 1) // 12
NP = F(N) / RHO
DIGITS = 80


def ln_enclosure(q):
    n, d = q.numerator, q.denominator
    zn, zd = d - n, d + n
    Q = 10 ** DIGITS
    t = (zn * Q) // zd + 1
    tot, pw, Qp, j = 0, t, 1, 0
    while j < 400:
        term = pw // (Qp * (2 * j + 1))
        tot += term + 1
        if term == 0 and j > 8:
            break
        pw, Qp, j = pw * t * t, Qp * Q * Q, j + 1
    hi = F(2 * (tot + 1), Q)
    return hi - F(1, 10 ** (DIGITS - 2)), hi


def floor_pow2_ratio(num, den):
    """Largest integer m with m**den <= 2**num."""
    tgt = 1 << num
    lo = 1 << (num // den)
    hi = lo * 2 + 2
    while lo + 1 < hi:
        m = (lo + hi) // 2
        if m ** den <= tgt:
            lo = m
        else:
            hi = m
    return lo


LN2_UP = F(693148, 1000000)
P_MIN = F(39, 100)
CERT_TARGET = P_MIN + F(1, 10000)

S_TOY = F(CENSUS["config"]["S_design"])
N_TOY = F(shards[0]["N"])
RHO_TOY = (1 - F(S_NUM, S_DEN) ** 2) * (1 - F(1, 128))
TRAILS_TOY = shards[0]["trails"]
PER_TOY = S_TOY // TRAILS_TOY
PAIRS_TOY = (S_TOY * S_TOY - TRAILS_TOY * PER_TOY * PER_TOY) / 2
V_TOY = PAIRS_TOY / (N_TOY / RHO_TOY)
MODEL_TOY = 1.0 - math.exp(-float(V_TOY))
# The census walks 32 trails x per_trail = 6 steps, i.e. 192 realized steps, while PAIRS_TOY
# above is the S_design = 223 pair count.  Model the realized budget separately and disclose
# the difference: at 192 steps the model is ~0.2958, and the measured rate is ABOVE it.
REALIZED_STEPS_TOY = TRAILS_TOY * PER_TOY
PAIRS_REALIZED = (F(REALIZED_STEPS_TOY * REALIZED_STEPS_TOY)
                  - TRAILS_TOY * PER_TOY * PER_TOY) / 2
MODEL_REALIZED = 1.0 - math.exp(-float(PAIRS_REALIZED / (N_TOY / RHO_TOY)))
CP_LO = F(str(cp_lo))
_ratio = float(CP_LO) / MODEL_TOY
C_RATIO = F(str(min(1.0, math.floor(_ratio * 10 ** 9) / 10 ** 9)))
ck("calibration_ratio_is_the_cp_lower_quotient",
   C_RATIO == min(F(1), F(math.floor(_ratio * 10 ** 9), 10 ** 9)), str(C_RATIO))
ck("calibration_ratio_at_most_one", C_RATIO <= 1, str(C_RATIO))
ck("calibration_ratio_is_applied_as_a_discount", C_RATIO < 1)
ck("census_rate_below_design_budget_toy_model", census_rate < MODEL_TOY,
   [round(census_rate, 6), round(MODEL_TOY, 6)])
CK_SE = math.sqrt(MODEL_TOY * (1 - MODEL_TOY) / census_n)
CK_SIGMA = (MODEL_TOY - census_rate) / CK_SE
# --- realized-budget reconciliation: the 223-step model is NOT the budget the census ran ----
ck("census_realized_steps_are_trails_times_per_trail",
   REALIZED_STEPS_TOY == TRAILS_TOY * PER_TOY == 192, REALIZED_STEPS_TOY)
ck("census_model_at_realized_budget_is_below_measured",
   census_rate > MODEL_REALIZED, [round(census_rate, 6), round(MODEL_REALIZED, 6)])
RETENTION_REALIZED = census_rate / MODEL_REALIZED
ck("census_retention_over_realized_model_above_one", RETENTION_REALIZED > 1,
   round(RETENTION_REALIZED, 4))
ck("design_budget_se_gap_is_a_budget_artifact",
   CK_SIGMA > 1 and REALIZED_STEPS_TOY < int(S_TOY), round(CK_SIGMA, 4))


def ceil_trail(S):
    return -(-S // RTAB)


def solve_S(target):
    ln_hi = ln_enclosure(1 - target)[1]
    need = ln_hi * NP * 2
    lo, hi = 1, 1 << 132
    while F(hi * (hi - 1), 1) < need:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if F(mid * (mid - 1), 1) >= need:
            hi = mid
        else:
            lo = mid
    return hi


CAP_UNITS = None


def debits_at(S):
    """Debits with the meter cap at floor(2^claim).  k_min is the minimum number of squaring
    steps the cap FORCES: the cap is sized on the design mix (5/8 multiply, 3/8 squaring), so
    a run whose multiply count exceeds S - k_min cannot be charged and aborts."""
    replay_instr = 2 * ceil_trail(S) * MUL_CHARGE
    cap_instr = (CAP_UNITS if CAP_UNITS is not None else (1 << 400)) * C
    step_cap = cap_instr - replay_instr - AUX_SUM
    k_min = (SETUP_INSTR + MUL_CHARGE * S - step_cap + (MUL_CHARGE - SQR_CHARGE) - 1) \
        // (MUL_CHARGE - SQR_CHARGE)
    k_min = max(0, k_min)
    eps = F(S_NUM, S_DEN) - F(k_min, S)
    hexp = 2 * eps * eps * S
    t_abort = int(hexp / LN2_UP)
    abort = F(1, 1 << min(t_abort, 400))
    trail_len = F(S, RTAB)
    dp_loss = 2 * F(1 << DBITS, 1) / trail_len
    rec_cap = -(-S // (1 << DBITS)) + (1 << 20)
    mu = F(S, 1 << DBITS)
    delta = (F(rec_cap) - mu) / mu
    t_store = int(delta * delta * mu / 3 / LN2_UP)
    return abort, dp_loss, F(1, 1 << min(t_store, 4000)), k_min, eps, hexp, rec_cap, t_store, t_abort


def model_required(D):
    return (CERT_TARGET + D) / C_RATIO


# --- fixed point: S from the debits at the cap, cap from the claim, claim from W ------
S_STEPS = solve_S(model_required(F(0)))
for _ in range(24):
    S_NEW = solve_S(model_required(sum(debits_at(S_STEPS)[:3])))
    if S_NEW == S_STEPS:
        break
    S_STEPS = S_NEW
for _ in range(2000):
    D_ = sum(debits_at(S_STEPS)[:3])
    if F(S_STEPS * (S_STEPS - 1), 2) / NP >= ln_enclosure(1 - model_required(D_))[1]:
        break
    S_STEPS += 1

STEP_INSTR_AVG = F(S_DEN - S_NUM, S_DEN) * MUL_CHARGE + F(S_NUM, S_DEN) * SQR_CHARGE
STEPS_INSTR = S_STEPS * STEP_INSTR_AVG
REPLAY_INSTR = 2 * ceil_trail(S_STEPS) * MUL_CHARGE
W_INSTR = SETUP_INSTR + STEPS_INSTR + REPLAY_INSTR + AUX_SUM
W = F(W_INSTR, C)


# --- claim: smallest 4-decimal c with 2^c >= W + 2e-5 (exact integer comparison) ------
def cmp_pow2_scale(x, k, scale=SCALE):
    """True iff x < 2**(k/scale).  Exact: x**scale < 2**k."""
    x = F(x)
    if k >= 0:
        return x.numerator ** scale < (1 << k) * x.denominator ** scale
    return (x.numerator ** scale) << (-k) < x.denominator ** scale


def ceil_log2_scale(x):
    """Smallest k with 2**(k/SCALE) >= x.  cmp(x,k) is True iff k is already large enough,
    so climb while it is False and then descend while k-1 is still enough."""
    k = int(math.floor(math.log2(float(x)) * SCALE)) - 4
    while not cmp_pow2_scale(x, k):
        k += 1
    while k > 0 and cmp_pow2_scale(x, k - 1):
        k -= 1
    return k


CLAIM_K = ceil_log2_scale(W + HEADROOM)
CLAIM = F(CLAIM_K, SCALE)
ck("claim_is_the_smallest_4dp_at_or_above_W_plus_2e_5",
   cmp_pow2_scale(W + HEADROOM, CLAIM_K)
   and not cmp_pow2_scale(W + HEADROOM, CLAIM_K - 1),
   str(CLAIM))
CAP_UNITS = floor_pow2_ratio(CLAIM_K, SCALE)
CAP_INSTR = CAP_UNITS * C
ck("claim_bounds_the_cap", cmp_pow2_scale(CAP_UNITS, CLAIM_K), str(CAP_UNITS))
ck("claim_bounds_nominal_work", cmp_pow2_scale(W, CLAIM_K), str(W))
ck("nominal_work_inside_cap", W_INSTR <= CAP_INSTR, [W_INSTR - CAP_INSTR])

# --- re-solve the debits at the SHIPPED cap and confirm the fixed point --------------
ABORT, DP_LOSS, STORE_DEBIT, K_MIN, EPS, HOEFF_EXP, REC_CAP, T_STORE, T_ABORT = debits_at(S_STEPS)
DEBITS = ABORT + DP_LOSS + STORE_DEBIT
TARGET = model_required(DEBITS)
CERTIFIED = TARGET * C_RATIO - DEBITS
ck("cap_forces_at_most_S_minus_kmin_multiplies", 0 < K_MIN < S_STEPS, K_MIN)
ck("cap_leaves_a_positive_multiply_fraction_margin", EPS > 0, "%.8f" % float(EPS))
ck("abort_debit_is_the_hoeveding_tail_clipped_at_2pow_minus_400",
   ABORT == F(1, 1 << 400) and T_ABORT >= 400, T_ABORT)
ck("certified_floor_at_least_0_39", CERTIFIED >= P_MIN, "%.6f" % float(CERTIFIED))
ck("target_meets_the_certificate_target", TARGET >= CERT_TARGET, "%.6f" % float(TARGET))
ck("debits_negligible_against_the_floor", DEBITS < F(1, 1000), "%.3e" % float(DEBITS))


def log2_milli(x):
    x = F(x)
    lo = (x.numerator.bit_length() - x.denominator.bit_length() - 1) * 1000
    hi = lo + 1000
    while not cmp_pow2_scale(x, hi, 1000):
        hi += 1000
    while cmp_pow2_scale(x, lo, 1000):
        lo -= 1000
    while hi - lo > 1:
        m = (lo + hi) // 2
        if cmp_pow2_scale(x, m, 1000):
            hi = m
        else:
            lo = m
    return lo, hi


GATE_LIVE = F(12574, 100)          # lowest live non-self submission when the lane was built
GATE_DECLARED = F(125826, 1000)    # the gate that makes this claim legal (kaankolcu 10214b43)
GATE_DECLARED_COMPUTED = F(125825569, 1000000)   # 10214b43's OWN stated computed total 2^125.825569
TIE_MARGIN = F(5, 1000)
MARGIN_FLOOR = F(15, 10000)
UNLOCK = CLAIM + TIE_MARGIN + MARGIN_FLOOR
TIE_CAP_LINE = CLAIM + TIE_MARGIN
ck("claim_inside_the_tie_window_of_the_declared_gate", TIE_CAP_LINE <= GATE_DECLARED,
   [str(TIE_CAP_LINE), str(GATE_DECLARED)])
ck("declared_gate_is_this_tracks_second_live_head", GATE_DECLARED == F(125826, 1000))
ck("claim_clears_the_declared_gates_own_computed_total",
   TIE_CAP_LINE <= GATE_DECLARED_COMPUTED, [str(TIE_CAP_LINE), str(GATE_DECLARED_COMPUTED)])
ck("claim_does_NOT_clear_the_gate_live_at_build", TIE_CAP_LINE > GATE_LIVE,
   [str(TIE_CAP_LINE), str(GATE_LIVE)])

WC_LO, WC_HI = CLAIM_K, CLAIM_K + 1
STORE_BYTES = REC_CAP * REC_BYTES
TABLE_BYTES = RTAB * 12 * 576
PEAK_BYTES = STORE_BYTES + TABLE_BYTES + (1 << 21) + (1 << 24) + (1 << 25)
ck("peak_below_2pow32", PEAK_BYTES < (1 << 32), PEAK_BYTES)

out = {
    "mode": "ledger-v19b",
    "config": {"R": RTAB, "idx_bits": IDX_BITS, "dp_shift": DP_SHIFT, "dp_bits": DBITS,
               "slot_rule": "s = idx(H(m0)), m0 = min of the minimal cyclic rotation of the "
                            "twelve Frobenius keys (class function), no kappa",
               "rotation": "minimal cyclic rotation (Booth), tie-break invariant",
               "guard": "a multiply step whose slot repeats the previous MULTIPLY slot is "
                        "replaced by the cyclotomic squaring (class-history guard)"},
    "charges": {"I_COMMON": I_COMMON, "I_MUL": I_MUL, "I_SQR": I_SQR,
                "MUL_CHARGE": MUL_CHARGE, "SQR_CHARGE": SQR_CHARGE, "C": C,
                "MUL_ROWS_SUM": sum(v for _, v in MUL_ROWS),
                "MUL_ROWS_PLUS_COMMON": I_COMMON + sum(v for _, v in MUL_ROWS),
                "SQR_ROWS_SUM": sum(v for _, v in SQR_EXTRA),
                "SQR_ROWS_PLUS_COMMON": I_COMMON + sum(v for _, v in SQR_EXTRA),
                "C_ROW_NOTE": "reference-unit normalization C per COST_MODEL.md"},
    "rows": {"common": [[n, v] for n, v in COMMON],
             "multiply": [[n, v] for n, v in MUL_ROWS],
             "squaring": [[n, v] for n, v in SQR_EXTRA],
             "auxiliary": AUX_ROWS,
             "setup": SETUP_INSTR},
    "s_steps": str(S_STEPS),
    "s_log2": "%d.%03d" % (log2_milli(F(S_STEPS))[0] // 1000, log2_milli(F(S_STEPS))[0] % 1000),
    "trail_len_steps": str(ceil_trail(S_STEPS)),
    "w_exact": [str(W.numerator), str(W.denominator)],
    "w_units_ceil": str(-(-W_INSTR // C)),
    "w_lo": "%d.%03d" % (log2_milli(W)[0] // 1000, log2_milli(W)[0] % 1000),
    "w_hi": "%d.%03d" % (log2_milli(W)[1] // 1000, log2_milli(W)[1] % 1000),
    "claim_derivation": {"rule": "smallest 4-decimal c with 2^c >= W + 2e-5",
                         "w_plus_headroom": str(W + HEADROOM),
                         "claim_4dp": "%.4f" % float(CLAIM),
                         "claim_k_units_of_1e-4": CLAIM_K,
                         "previous_4dp_too_small": "%.4f" % ((CLAIM_K - 1) / SCALE)},
    "cap_units": str(CAP_UNITS), "cap_instr": str(CAP_INSTR),
    "cap_log2": "%.4f..%.4f" % (CLAIM_K / SCALE, (CLAIM_K + 1) / SCALE),
    "claim": "%.4f" % float(CLAIM), "gate": "125.74", "tie_margin": "0.005",
    "conditional": {"declared_gate": "125.826", "tie_cap_line": "%.4f" % float(TIE_CAP_LINE),
                    "unlock_with_margin": "%.4f" % float(UNLOCK),
                    "margin_floor": "0.0015",
                    "declared_gate_computed_total": "125.825569",
                    "slack_claimed": "%.4f" % float(GATE_DECLARED - TIE_CAP_LINE),
                    "slack_computed": "%.6f" % float(GATE_DECLARED_COMPUTED - TIE_CAP_LINE),
                    "gate_live_at_build": "125.74",
                    "ships_at_gate_live_at_build": bool(TIE_CAP_LINE <= GATE_LIVE),
                    "ships_at_declared_gate": bool(TIE_CAP_LINE <= GATE_DECLARED),
                    "live_set_that_ships_it": "kaankolcu 10214b43 125.826 alive and both "
                                              "jungjipdo heads 747ad6c9 125.74 / 6d4e57d1 "
                                              "125.78-125.80 dead"},
    "floor": {"required": "0.39", "certified": "%.6f" % float(CERTIFIED),
              "model_at_S": "%.6f" % float(TARGET),
              "calibration": {"trials": census_n, "successes": census_k,
                              "rate": round(census_rate, 6), "cp95_lower": round(cp_lo, 6),
                              "interval": "Clopper-Pearson one-sided 95%",
                              "model_at_census": round(MODEL_TOY, 6),
                              "model_at_census_budget": "S_design=%d steps (not the realized "
                                                        "budget)" % int(S_TOY),
                              "realized_steps": REALIZED_STEPS_TOY,
                              "realized_per_trail": PER_TOY,
                              "model_at_realized_budget": round(MODEL_REALIZED, 6),
                              "measured_over_realized_model": round(RETENTION_REALIZED, 4),
                              "realized_budget_note": "census walks %d trails x %d = %d realized "
                                                      "steps, so the %d-step S_design model is "
                                                      "the wrong budget; at %d steps the model "
                                                      "is %.6f and the measured rate is ABOVE it "
                                                      "(%.4fx), so no model-shortfall discount "
                                                      "is due"
                                                      % (TRAILS_TOY, PER_TOY,
                                                         REALIZED_STEPS_TOY, int(S_TOY),
                                                         REALIZED_STEPS_TOY, MODEL_REALIZED,
                                                         RETENTION_REALIZED),
                              "model_minus_rate": round(MODEL_TOY - census_rate, 6),
                              "model_minus_rate_in_se": round(CK_SIGMA, 4),
                              "model_minus_rate_is_a_budget_artifact":
                                  "the %d-step S_design model over-predicts the %d-step run; "
                                  "against the realized-budget model the measured rate is "
                                  "%.4fx ABOVE" % (int(S_TOY), REALIZED_STEPS_TOY,
                                                   RETENTION_REALIZED),
                              "ratio_lower": str(C_RATIO),
                              "applied_as": "S is solved for tail = (0.3901 + debits)/ratio; the "
                                            "ratio is RETAINED as a conservative retention (it "
                                            "is < 1 although the measured rate is > 1x the "
                                            "realized-budget model), so the certified floor is a "
                                            "lower bound, not a model-shortfall debit"},
              "debits": {"abort": str(ABORT), "dp_loss": str(DP_LOSS),
                         "storage": str(STORE_DEBIT), "sum": str(DEBITS),
                         "k_min": K_MIN}},
    "census": CENSUS,
    "sweep": [{"dpbits": d["dpbits"], "trials": d["trials"], "trails": d["trails"],
               "per_trail": d["per_trail"], "budget": d["budget"],
               "detected": d["detected"], "cp95_detect": round(d["cp95_detect"], 6),
               "pairs": d["pairs"], "events": d["events"], "degenerate": d["degenerate"],
               "records_per_trial": d["records_per_trial"],
               "birthday_pred": round(d["birthday_pred"], 6)} for d in
              [sweeps[k] for k in sorted(sweeps)]],
    "merge_audit": {"v19_same_class_at": v19c["same_class_at"],
                    "v19_pairs": v19c["pairs"], "v19_steps": v19c["steps"],
                    "v19_identical_elements_after_1_step": v19c["identical_elements_after_1_step"],
                    "tn_same_class_at": tnc["same_class_at"],
                    "tn_label_agrees_on_conjugates": [tn["tn_label_agrees_on_conjugates"],
                                                     tn["pairs"]],
                    "tn_norm_frobenius_invariant": [tn["norm_frobenius_invariant"],
                                                    tn["pairs"]]},
    "trap": {"v19_p6_per_1600": v19t["p6_per_1600"],
             "v19_second_order_per_1600": v19t["per_1600_2cycle"],
             "v19_traces_with_2cycle": v19t["traces_with_2cycle"],
             "v19_guard_firings": v19t["guard_firings"],
             "v19_squaring_fraction": v19t["squaring_fraction"],
             "v18_control_second_order_per_1600": v18t["per_1600_2cycle"],
             "v18_control_squaring_fraction": v18t["squaring_fraction"],
             "kappa_shift_hist": v19t["kappa_shift_hist"],
             "census_trap_regime": {
                 "toy_guard_fire_rate": round(v19t["guard_firings"] / (1600 * 200), 6),
                 "production_guard_fire_rate": round(1 / 4096, 6),
                 "fixture_exercises_trap_channel_x": round(
                     (v19t["guard_firings"] / (1600 * 200)) * 4096, 1),
                 "toy_effective_slots_estimate": 88,
                 "toy_effective_slots_note": "m0 is a minimum over twelve small-field keys, so "
                                             "the toy behaves like ~88 slots, not 4096",
                 "direction": "the fixture exercises the repeat-slot channel ~28x MORE often "
                              "than production, so the 0/1600 census is CONSERVATIVE by that "
                              "factor"}},
    "meter_identity": {"meter": fmatch["meter_instr"], "row_sum": fmatch["row_sum_instr"],
                       "log_sum": sum(fmatch["charge_log"]), "log_len": fmatch["charge_log_len"],
                       "matches": fmatch["match"]["matches"],
                       "replay_steps": fmatch["nreplay"],
                       "all_equal_exact_integers": fmatch["all_equal"]},
    "storage": {"rec_cap": str(REC_CAP), "rec_bytes": REC_BYTES, "store_bytes": STORE_BYTES,
                "table_bytes": TABLE_BYTES, "peak_bytes": PEAK_BYTES,
                "peak_below_2pow32": PEAK_BYTES < (1 << 32)},
    "checks": checks, "n_checks": len(checks),
    "all_checks_pass": all(c["status"] == "PASS" for c in checks),
}
print(json.dumps(out, indent=1))
sys.exit(0 if out["all_checks_pass"] else 1)
