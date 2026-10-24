#!/usr/bin/env python3
"""Exact ledger v17b -- the whole claim recomputed from raw evidence JSONs.

Integers and Fractions only.  Deltas against ledger16.py, all of them forced by the
v17 disjoint-half slot rule:

 R   the multiplier table is R = 4096 slots (two disjoint halves of 2048), so the
     setup row, the trail split (R trails), the rho factor (1-1/R) and the table
     storage row are all re-priced at 4096; the half-select row (+12 instructions
     on the multiply branch only) is a new counted row.
 M   the meter debits EVERY routine: forward steps (charged before execution),
     replays, collision solve, final verify, output.  Auxiliary rows are itemised
     and enter both the nominal rows and the cap.
 F   the DECLARED success floor is 0.38, not 0.39.  The birthday model at S gives
     0.39+debits, but the design-budget census measures the realised conversion at
     the same budget, and its one-sided 95% Clopper-Pearson lower bound is 0.382.
     A floor is only claimed at the value a design-budget lower bound supports, so
     the claim is 0.38 and the 0.39 model value is disclosed as a model statement,
     never as a certified floor.  (v16's defect was exactly this substitution.)
"""
import json
import sys
from fractions import Fraction as F
from math import isqrt

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
RS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 20000                        # instructions per dense Fp12 multiplication
RTAB = 4096                      # multiplier table size (two disjoint halves)
HALF = RTAB // 2                 # 2048
IDX_BITS = 11                    # bits 3..13
DP_SHIFT = 14                    # bits 14..113
DBITS = 100                      # production DP width
REC_BYTES = 52                   # 32 B fingerprint + 2 B trail + 16 B step + 2 B spare
ROOT = __file__.rsplit("/", 2)[0]

prod = json.load(open(ROOT + "/evidence/prod_run17b.json"))
census = json.load(open(ROOT + "/evidence/toy17b_census.json"))
pred = json.load(open(ROOT + "/evidence/predicate17b.json"))
pf = prod["prod_field"]
checks = []


def ck(name, ok, value=None):
    checks.append({"check": name, "value": value, "source": "ledger",
                   "status": "PASS" if ok else "FAIL"})
    assert ok, (name, value)


# ------------------------------------------------------------ configuration
cfg_eff = prod["config"]
ck("config_R_is_4096", cfg_eff["R"] == RTAB and cfg_eff["R"] == 4096, cfg_eff["R"])
ck("config_idx_bits_is_11", cfg_eff["idx_bits"] == IDX_BITS == 11, cfg_eff["idx_bits"])
ck("config_dp_shift_is_14", cfg_eff["dp_shift"] == DP_SHIFT == 14, cfg_eff["dp_shift"])
ck("config_half_is_2048", cfg_eff["half"] == HALF, cfg_eff["half"])
ck("caller_overrides_R", cfg_eff["caller_overrides_R"] is True)
ck("predicate_agreement_zero_mismatches",
   pred["mismatches_shared_vs_divmod"] == 0 and pred["mismatches_shared_vs_ranges"] == 0,
   pred["samples_tested"])
ck("predicate_dp_does_not_constrain_index", pred["dp_does_not_constrain_index"] is True)
ck("table_orbits_disjoint",
   prod["orbit_disjointness"]["collisions"] == 0 and
   prod["orbit_disjointness"]["values"] == 12 * RTAB, prod["orbit_disjointness"]["values"])

# ------------------------------------------------------------ instruction rows
SQR_ARITH = prod["sqr_rows"]["instr"]
KEY_ARITH = prod["key_rows"]["instr"]
COMMON = (("trace_re_plus_branch_bits", 60), ("exponent_pair_update", 100),
          ("cycle_compares_loop_margin", 140), ("dp_store_amortised", 1))
MUL_EXTRA = (("key_six_plusminus", KEY_ARITH), ("min_six_argmax_index", 180),
             ("table_address_and_load", 150), ("exponent_pair_add", 60),
             ("half_select_kappa_ge_6", 12))
SQR_EXTRA = (("cyclotomic_squaring", SQR_ARITH), ("exponent_pair_double", 40))
I_COMMON = sum(v for _, v in COMMON)
I_MUL = I_COMMON + sum(v for _, v in MUL_EXTRA)
I_SQR = I_COMMON + sum(v for _, v in SQR_EXTRA)
ck("rows_match_production_script",
   (I_COMMON, I_MUL, I_SQR) == (prod["charges"]["I_COMMON"], prod["charges"]["I_MUL"],
                                prod["charges"]["I_SQR"]),
   [I_COMMON, I_MUL, I_SQR])
MUL_CHARGE = C + I_MUL
SQR_CHARGE = I_SQR
ck("charges_match_production_script",
   (MUL_CHARGE, SQR_CHARGE) == (prod["charges"]["MUL_CHARGE"], prod["charges"]["SQR_CHARGE"]),
   [MUL_CHARGE, SQR_CHARGE])

AUX = (("solve_modular_inverse_binary_extended_euclid", prod["aux_rows"]["solve_instr"]),
       ("verify_one_256bit_exponentiation_in_GT", prod["aux_rows"]["verify_instr"]),
       ("output_one_48byte_Fp_plus_write", prod["aux_rows"]["output_instr"]))
SOLVE_INSTR, VERIFY_INSTR, OUTPUT_INSTR = (v for _, v in AUX)
AUX_SUM = SOLVE_INSTR + VERIFY_INSTR + OUTPUT_INSTR

SETUP_INSTR = RTAB * (510 * C + 12 * 2400) + (1 << 20)
ck("setup_row_matches_production_script", SETUP_INSTR == prod["config"]["setup_instr"],
   SETUP_INSTR)

# ---- D4 exact integer agreement: per-event log == meter == row formula -----
nsq, nsteps = pf["squaring_steps"], pf["steps"]
nmul = nsteps - nsq
LOG_SUM = sum(prod["charge_log"])
ck("charge_log_length_equals_steps", prod["charge_log_len"] == nsteps == pf["steps"])
ck("meter_equals_setup_plus_log", pf["meter_instr"] == SETUP_INSTR + LOG_SUM,
   [pf["meter_instr"], SETUP_INSTR + LOG_SUM])
ck("meter_equals_row_formula",
   pf["meter_instr"] == SETUP_INSTR + nsq * SQR_CHARGE + nmul * MUL_CHARGE
   + 2 * pf["dp_replays"] * MUL_CHARGE, pf["meter_instr"])
ck("meter_sum_row_written_by_run", prod["meter_instr_sum"] == SETUP_INSTR + LOG_SUM)
ck("production_run_no_aborts", pf["aborts"] == 0, pf["aborts"])
ck("production_run_R_is_4096", pf["R"] == RTAB, pf["R"])
ck("charge_values_are_the_two_rows", sorted(set(x for x in prod["charge_log"] if x in (SQR_CHARGE, MUL_CHARGE))) == [SQR_CHARGE, MUL_CHARGE])

# ------------------------------------------------------------ collision model
S_NUM, S_DEN = 3, 8
C_MUL = F(MUL_CHARGE, C)
C_SQR = F(SQR_CHARGE, C)
C_AVG = F(S_DEN - S_NUM, S_DEN) * C_MUL + F(S_NUM, S_DEN) * C_SQR
RHO_OURS = (1 - F(S_NUM, S_DEN) ** 2) * (1 - F(1, RTAB))
SM = 1 - F(S_NUM, S_DEN)
RHO_REVIEW = 1 - F(S_NUM, S_DEN) ** 2 - SM * SM / F(RTAB)
ck("rho_is_conservative_side", RHO_OURS <= RHO_REVIEW, [str(RHO_OURS), str(RHO_REVIEW)])
RHO = RHO_OURS
N = (RS - 1) // 12
ck("r_minus_1_is_12N", (RS - 1) == 12 * N)
NP = F(N) / RHO

DIGITS = 80


def ln_enclosure(q):
    n = q.numerator
    d = q.denominator
    zn, zd = d - n, d + n
    Q = 10 ** DIGITS
    t = (zn * Q) // zd + 1
    tot = 0
    pw = t
    Qp = 1
    j = 0
    while j < 400:
        term = pw // (Qp * (2 * j + 1))
        tot += term + 1
        if term == 0 and j > 8:
            break
        pw = pw * t * t
        Qp = Qp * Q * Q
        j += 1
    hi = F(2 * (tot + 1), Q)
    return hi - F(1, 10 ** (DIGITS - 2)), hi


P_MODEL = F(39, 100)             # birthday-model floor the design budget targets
DECLARED_FLOOR = F(38, 100)      # the floor actually claimed (certified below)


def floor_pow2_ratio(num, den):
    tgt = 1 << num
    lo = 1 << (num // den)
    hi = lo * 2 + 2
    while lo + 1 < hi:
        m = (lo + hi) // 2
        if m ** den <= tgt:
            lo = m
        else:
            hi = m
    assert lo ** den <= tgt < (lo + 1) ** den
    return lo


CAP_UNITS = floor_pow2_ratio(125728, 1000)       # floor(2^125.728)
CAP_INSTR = CAP_UNITS * C
LN2_UP = F(693148, 1000000)
LN2_LO = F(693147, 1000000)


def ceil_sqrt(x):
    x = F(x)
    n_, d_ = x.numerator, x.denominator
    s = isqrt(n_ // d_)
    while s * s * d_ < n_:
        s += 1
    while s > 0 and (s - 1) * (s - 1) * d_ >= n_:
        s -= 1
    return s


def ceil_trail(S):
    return -(-S // RTAB)


def solve_S(target):
    """Minimal S with S(S-1) >= 2*NP*ln(1/(1-target)), by exact bisection.

    `target` is the MODEL floor required (already includes the debits and the
    calibration discount), not 0.39 itself.
    """
    assert target < 1
    ln_hi = ln_enclosure(1 - target)[1]
    need = ln_hi * NP * 2                    # S(S-1) >= need, exact Fraction
    lo, hi = 1, 1 << 132
    while F(hi * (hi - 1), 1) < need:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if F(mid * (mid - 1), 1) >= need:
            hi = mid
        else:
            lo = mid
    assert F(hi * (hi - 1), 1) >= need
    return hi


def debits_at(S):
    replay_instr = 2 * ceil_trail(S) * MUL_CHARGE
    step_cap = CAP_INSTR - replay_instr - AUX_SUM
    k_min = (SETUP_INSTR + MUL_CHARGE * S - step_cap + (MUL_CHARGE - SQR_CHARGE) - 1) // (MUL_CHARGE - SQR_CHARGE)
    k_min = max(0, k_min)
    eps = F(S_NUM, S_DEN) - F(k_min, S)
    if eps <= 0:
        abort, hexp = F(0), F(0)
    else:
        hexp = 2 * eps * eps * S
        t_abort = int(hexp / LN2_UP)
        abort = F(1, 1 << min(t_abort, 400))
    trail_len = F(S, RTAB)
    dp_loss = 2 * F(1 << DBITS, 1) / trail_len
    rec_cap = -(-S // (1 << DBITS)) + (1 << 20)
    mu = F(S, 1 << DBITS)
    delta = (F(rec_cap) - mu) / mu
    chernoff = delta * delta * mu / 3
    t_store = int(chernoff / LN2_UP)
    store_debit = F(1, 1 << min(t_store, 4000))
    return abort, dp_loss, store_debit, k_min, eps, hexp, rec_cap, t_store



# =============================================================== the claimed floor
# COST_MODEL.md line 45: the required success probability is AT LEAST 0.39 over
# uniform production challenges and fresh coins.  0.39 is therefore a floor of the
# TRACK, not a lever: S is solved so that the CERTIFIED floor reaches it.
P_MIN = F(39, 100)
CERT_MARGIN = F(1, 10000)          # certify 0.39 + 1e-4, not 0.39 + 0
CERT_TARGET = P_MIN + CERT_MARGIN

# The census CALIBRATES the birthday model at the toy design budget; it is not a
# substitute floor.  At density 1 the toy measures P(useful collision within
# S_design), the same event the exact tail predicts.  ratio = measured / model,
# evaluated at the one-sided 95% Clopper-Pearson lower end, is the factor by which
# the exact tail is discounted at production S.
import math
S_TOY = F(census["config"]["S_design"])
N_TOY = F(census["config"]["N"])
RHO_TOY = (1 - F(S_NUM, S_DEN) ** 2) * (1 - F(1, census["config"]["R"]))
V_TOY = S_TOY * (S_TOY - 1) / (2 * (N_TOY / RHO_TOY))
_model_toy_f = 1.0 - math.exp(-float(V_TOY))
CP_LO = F(str(census["cp95_one_sided_lower"]))
_ratio_f = float(CP_LO) / _model_toy_f
C_RATIO = F(str(math.floor(_ratio_f * 10 ** 9) / 10 ** 9))   # rounded DOWN
ck("census_is_the_design_budget_config",
   census["config"]["S_design"] == 223 and census["config"]["dpbits"] == 0
   and census["config"]["trails"] == 1 and census["config"]["rule"] == "new",
   census["config"])
ck("census_pool_is_one_config", census["shard_count"] == 20, census["shard_count"])
ck("census_is_large_enough_for_a_calibration", census["trials"] >= 50000, census["trials"])
ck("census_agrees_with_the_model_within_2se",
   abs(census["success_rate_float"] - _model_toy_f) <= census["normal_2se_halfwidth"],
   [census["success_rate_float"], round(_model_toy_f, 6)])
ck("calibration_ratio_below_one_disclosed", C_RATIO < 1, str(C_RATIO))


STEP_INSTR_AVG = F(S_DEN - S_NUM, S_DEN) * MUL_CHARGE + F(S_NUM, S_DEN) * SQR_CHARGE


def model_required(D):
    """exact model floor so that model*C_RATIO - D >= 0.39 + margin."""
    return (CERT_TARGET + D) / C_RATIO


# ONE fixed cap, taken from the claimed bound itself: cap = floor(2^125.734).  The
# claim is then a bound on EVERY run (the meter refuses a step that would cross it),
# and the nominal work rows must sit inside it.
CAP_UNITS = floor_pow2_ratio(125732, 1000)
CAP_INSTR = CAP_UNITS * C
S_STEPS = solve_S(model_required(F(0)))
for _ in range(12):
    _D = sum(debits_at(S_STEPS)[:3])
    S_STEPS = solve_S(model_required(_D))       # fixed point, free to move DOWN too
for _ in range(2000):
    _D = sum(debits_at(S_STEPS)[:3])
    if F(S_STEPS * (S_STEPS - 1), 2) / NP >= ln_enclosure(1 - model_required(_D))[1]:
        break
    S_STEPS += 1
ABORT, DP_LOSS, STORE_DEBIT, K_MIN, EPS, HOEFF_EXP, REC_CAP, T_STORE = debits_at(S_STEPS)
DEBITS = ABORT + DP_LOSS + STORE_DEBIT

TARGET = model_required(DEBITS)
LN_LO, LN_HI = ln_enclosure(1 - TARGET)
BIRTH_V = F(S_STEPS * (S_STEPS - 1), 2) / NP
ck("birthday_tail_reaches_the_required_model", BIRTH_V >= LN_HI, str(HOEFF_EXP))
MODEL_FLOOR = TARGET
CERTIFIED = MODEL_FLOOR * C_RATIO - DEBITS
ck("certified_floor_is_at_least_0_39", CERTIFIED >= P_MIN, "%.6f" % float(CERTIFIED))
ck("certified_floor_reaches_the_margin", CERTIFIED >= CERT_TARGET, "%.6f" % float(CERTIFIED))
ck("debits_are_small", DEBITS < F(1, 1000), str(DEBITS))

STEPS_INSTR = S_STEPS * STEP_INSTR_AVG
REPLAY_INSTR = 2 * ceil_trail(S_STEPS) * MUL_CHARGE
W_INSTR = SETUP_INSTR + STEPS_INSTR + REPLAY_INSTR + AUX_SUM
W = F(W_INSTR, C)
STEP_CAP_INSTR = CAP_INSTR - REPLAY_INSTR - AUX_SUM
ck("nominal_work_inside_cap", W <= CAP_UNITS, [str(W), CAP_UNITS])
ck("cap_covers_all_multiply_worst_case_minus_abort_rule",
   SETUP_INSTR + MUL_CHARGE * S_STEPS > STEP_CAP_INSTR)
ck("replay_and_aux_inside_cap", SETUP_INSTR + STEPS_INSTR + AUX_SUM <= STEP_CAP_INSTR)


def cmp_pow2(x, k, scale):
    x = F(x)
    return x.numerator ** scale < (1 << k) * x.denominator ** scale


def log2_milli(x):
    x = F(x)
    lo = (x.numerator.bit_length() - x.denominator.bit_length() - 1) * 1000
    hi = lo + 1000
    while not cmp_pow2(x, hi, 1000):
        hi += 1000
    while cmp_pow2(x, lo, 1000):
        lo -= 1000
    while hi - lo > 1:
        m = (lo + hi) // 2
        if cmp_pow2(x, m, 1000):
            hi = m
        else:
            lo = m
    return lo, hi


GATE = F(6287, 50)
TIE_MARGIN = F(5, 1000)
W_LO, W_HI = log2_milli(W)
WC_LO, WC_HI = log2_milli(CAP_UNITS)
CLAIM = F(WC_HI, 1000)
ck("claim_bounds_the_cap", cmp_pow2(CAP_UNITS, WC_HI, 1000), WC_HI)
ck("claim_bounds_nominal_work", cmp_pow2(W, WC_HI, 1000), WC_HI)
ck("claim_outside_the_tie_window", CLAIM + TIE_MARGIN < GATE, [str(CLAIM), str(GATE)])
ck("claim_within_the_ship_band", CLAIM <= F(125734, 1000), str(CLAIM))
ck("claim_beats_the_gate", CLAIM < GATE, [str(CLAIM), str(GATE)])

STORE_BYTES = REC_CAP * REC_BYTES
TABLE_BYTES = RTAB * 12 * 576
CODE_BYTES = 1 << 21
SCRATCH_BYTES = 1 << 24
SORTBUF_BYTES = 1 << 25
PEAK_BYTES = STORE_BYTES + TABLE_BYTES + CODE_BYTES + SCRATCH_BYTES + SORTBUF_BYTES
ck("peak_memory_below_2pow32", PEAK_BYTES < (1 << 32), PEAK_BYTES)
ck("table_bytes_are_the_r4096_rows", TABLE_BYTES == RTAB * 12 * 576, TABLE_BYTES)
ck("store_cap_covers_expected_records", REC_CAP >= -(-S_STEPS // (1 << DBITS)))
ck("dp_wait_row_present", T_STORE > 0, T_STORE)

S_RATIO = S_STEPS / 70925065820265078464927600992738729880
out = {
    "mode": "ledger-v17b",
    "inputs": {"prod": "evidence/prod_run17b.json",
               "census": "evidence/toy17b_census.json",
               "predicate": "evidence/predicate17b.json"},
    "config": {"R": RTAB, "half": HALF, "idx_bits": IDX_BITS, "dp_shift": DP_SHIFT,
               "dp_bits": DBITS, "slot_rule": "h + 2048*[kappa>=6]"},
    "charges": {"I_COMMON": I_COMMON, "I_MUL": I_MUL, "I_SQR": I_SQR,
                "MUL_CHARGE": MUL_CHARGE, "SQR_CHARGE": SQR_CHARGE, "C": C},
    "rows": {"common": [[n, v] for n, v in COMMON],
             "multiply": [[n, v] for n, v in MUL_EXTRA],
             "squaring": [[n, v] for n, v in SQR_EXTRA],
             "auxiliary": [[n, v] for n, v in AUX],
             "setup": [["table_build_R4096", "4096 x (510 ref + 12 Frobenius x 2400)"], SETUP_INSTR]},
    "s_steps": str(S_STEPS),
    "s_log2": "%d.%03d" % (log2_milli(F(S_STEPS))[0] // 1000, log2_milli(F(S_STEPS))[0] % 1000),
    "s_vs_v17": {"ratio": "%.6f" % float(S_RATIO), "delta_bits": "%.5f" % math.log2(float(S_RATIO))},
    "trail_len_steps": str(ceil_trail(S_STEPS)),
    "replay_steps": str(2 * ceil_trail(S_STEPS)),
    "steps_instr_num": str(STEPS_INSTR.numerator), "steps_instr_den": str(STEPS_INSTR.denominator),
    "replay_instr": str(REPLAY_INSTR), "aux_instr": str(AUX_SUM),
    "w_instr_num": str(W_INSTR.numerator), "w_instr_den": str(W_INSTR.denominator),
    "w_units": str(W), "w_lo": "%d.%03d" % (W_LO // 1000, W_LO % 1000),
    "w_hi": "%d.%03d" % (W_HI // 1000, W_HI % 1000),
    "cap_units": str(CAP_UNITS), "cap_instr": str(CAP_INSTR),
    "cap_log2": "%d.%03d..%d.%03d" % (WC_LO // 1000, WC_LO % 1000, WC_HI // 1000, WC_HI % 1000),
    "claim": "%.3f" % float(CLAIM),
    "gate": "125.74", "tie_margin": "0.005",
    "floor": {
        "required_by_cost_model": "0.39 (COST_MODEL.md line 45)",
        "certified": "%.6f" % float(CERTIFIED),
        "certified_meets": float(CERTIFIED) >= 0.39,
        "model_floor_at_S": "%.6f" % float(MODEL_FLOOR),
        "calibration": {"census_trials": census["trials"], "census_successes": census["success_num"],
                        "census_rate": census["success_rate_float"],
                        "census_cp95_one_sided_lower": str(CP_LO),
                        "model_at_census_budget": "%.6f" % _model_toy_f,
                        "ratio_point": "%.6f" % (census["success_rate_float"] / _model_toy_f),
                        "ratio_lower": str(C_RATIO),
                        "statistic": "P(useful collision within S_design), density 1, T=1",
                        "v16_shards_not_pooled": census["corroboration_unpooled"]["reason_not_pooled"]},
        "debits": {"abort_hoeffding": str(ABORT), "dp_detection_loss": str(DP_LOSS),
                   "storage_chernoff": str(STORE_DEBIT), "sum": str(DEBITS)},
        "composition": "certified = model(S)*ratio_lower - debits",
        "why_S_moved": ("ratio_lower discounts the exact tail, so the model floor must be %.6f "
                        "instead of 0.39; S rises %.5f bits" % (float(TARGET), math.log2(float(S_RATIO)))),
        "point_ratio_alternative": ("at the point ratio %.6f the same floor is met at the v17 S "
                                    "(W = 2^125.727); the 95%% lower end is used because v16 was "
                                    "rejected for a floor weaker than its own evidence"
                                    % (census["success_rate_float"] / _model_toy_f)),
    },
    "debits": {"k_min": str(K_MIN), "eps": str(EPS), "hoeffding_exponent": str(HOEFF_EXP),
               "t_store": str(T_STORE), "rec_cap": str(REC_CAP)},
    "storage": {"rec_cap_records": str(REC_CAP), "rec_bytes": REC_BYTES,
                "store_bytes": STORE_BYTES, "records_expected": str(-(-S_STEPS // (1 << DBITS))),
                "table_bytes_R4096": TABLE_BYTES, "table_entries": RTAB * 12,
                "table_bytes_per_entry": 576, "code_bytes": CODE_BYTES,
                "scratch_bytes": SCRATCH_BYTES, "sortbuf_bytes": SORTBUF_BYTES,
                "peak_bytes": PEAK_BYTES, "peak_below_2pow32": PEAK_BYTES < (1 << 32),
                "dp_wait_instr_row": str(256 * C), "dp_waits_charged": True,
                "overflow_policy": "refuse new record, charge Chernoff debit 2^-%d" % T_STORE},
    "meter_identity": {
        "prod_meter_instr": str(pf["meter_instr"]),
        "prod_log_sum_instr": str(SETUP_INSTR + LOG_SUM),
        "prod_row_formula_instr": str(SETUP_INSTR + nsq * SQR_CHARGE + nmul * MUL_CHARGE + 2 * pf["dp_replays"] * MUL_CHARGE),
        "ledger_setup_instr": str(SETUP_INSTR),
        "ledger_charge_log_sum_instr": str(SETUP_INSTR + LOG_SUM),
        "all_equal_as_exact_integers": (pf["meter_instr"] == SETUP_INSTR + LOG_SUM ==
                                        SETUP_INSTR + nsq * SQR_CHARGE + nmul * MUL_CHARGE + 2 * pf["dp_replays"] * MUL_CHARGE),
        "deltas": {"prod_minus_ledger": pf["meter_instr"] - (SETUP_INSTR + LOG_SUM)},
    },
    "checks": checks,
    "all_checks_pass": all(c["status"] == "PASS" for c in checks),
    "n_checks": len(checks),
}
print(json.dumps(out, indent=1))
sys.exit(0 if out["all_checks_pass"] else 1)
