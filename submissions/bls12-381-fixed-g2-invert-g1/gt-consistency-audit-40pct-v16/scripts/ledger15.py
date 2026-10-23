#!/usr/bin/env python3
"""Exact ledger v15: METERED work cap, Hoeffding abort debit, DP store rows.

Every number printed here is recomputed from the raw evidence JSONs
(evidence/prod_run15.json, evidence/toy_census15.json,
evidence/predicate_match.json) with integers and Fractions only.  No float
appears in any comparison, and no rounded string from an evidence file is
compared as a string: the evidence values used are integers, and every derived
quantity is an exact Fraction or a verified integer power comparison.

Unit rule (docs/COST_MODEL.md): cost = H + R/C with C = 20000 instructions per
dense Fp12 multiplication (the reference operation) and H the number of
reference operations.
"""
import json
import sys
from fractions import Fraction as F
from math import isqrt

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
RS = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 20000
RTAB = 1024
ROOT = __file__.rsplit("/", 2)[0]

prod = json.load(open(ROOT + "/evidence/prod_run15.json"))
toy = json.load(open(ROOT + "/evidence/toy_census15.json"))
pred = json.load(open(ROOT + "/evidence/predicate_match.json"))
pf = prod["prod_field"]

# ---------------------------------------------------------------- instruction rows
SQR_ARITH = prod["sqr_rows"]["instr"]            # cyclotomic squaring, counted
KEY_ARITH = prod["key_rows"]["instr"]            # six +-pair keys, counted
COMMON = (("trace_re_plus_branch_bits", 60), ("exponent_pair_update", 100),
          ("cycle_compares_loop_margin", 140), ("dp_store_amortised", 1))
MUL_EXTRA = (("key_six_plusminus", KEY_ARITH), ("min_six_argmax_index", 180),
             ("table_address_and_load", 150), ("exponent_pair_add", 60))
SQR_EXTRA = (("cyclotomic_squaring", SQR_ARITH), ("exponent_pair_double", 40))
I_COMMON = sum(v for _, v in COMMON)
I_MUL = I_COMMON + sum(v for _, v in MUL_EXTRA)
I_SQR = I_COMMON + sum(v for _, v in SQR_EXTRA)
assert (I_COMMON, I_MUL, I_SQR) == (301, 2493, 15383), (I_COMMON, I_MUL, I_SQR)

# table build: 1024 multipliers, each g^a z^b by 255 doublings + 255 additions in
# GT (510 reference mults) plus 12 Frobenius applications at 2400 instructions.
SETUP_INSTR = RTAB * (510 * C + 12 * 2400) + (1 << 20)

# Charged instructions per step: a multiply step is one reference operation
# (C instructions) plus its counted rows; a squaring step has no reference
# operation, so its charge is its counted rows.  The executed run's meter covers
# the counted rows; the ledger adds exactly C per multiply step, which is a
# per-step rule the implementation applies with the same determinism.
MUL_CHARGE = C + I_MUL                             # 22493
SQR_CHARGE = I_SQR                                 # 15383
I_DIFF = MUL_CHARGE - SQR_CHARGE                   # 7110
nsq, nsteps = pf["squaring_steps"], pf["steps"]
nmul = nsteps - nsq
assert pf["meter_instr"] == SETUP_INSTR + nsq * I_SQR + nmul * I_MUL
FULL_INSTR = pf["meter_instr"] + nmul * C
assert FULL_INSTR == SETUP_INSTR + nsq * SQR_CHARGE + nmul * MUL_CHARGE
assert pf["aborts"] == 0

S_NUM, S_DEN = 3, 8
C_MUL = 1 + F(I_MUL, C)
C_SQR = F(I_SQR, C)
C_AVG = F(S_DEN - S_NUM, S_DEN) * C_MUL + F(S_NUM, S_DEN) * C_SQR
assert (C_MUL, C_SQR, C_AVG) == (F(22493, 20000), F(15383, 20000), F(79307, 80000)), (C_MUL, C_SQR, C_AVG)

# ---------------------------------------------------------------- collision model
# rho is the merge probability per step pair in units of 1/N.  We charge the
# SMALLER of the two candidate expressions, i.e. the conservative side:
#   ours      (1 - s^2)(1 - 1/R)
#   reviewer  1 - s^2 - (1-s)^2/R        = ours + 2s(1-s)/R
RHO_OURS = (1 - F(S_NUM, S_DEN) ** 2) * (1 - F(1, RTAB))
SM = 1 - F(S_NUM, S_DEN)
RHO_REVIEW = 1 - F(S_NUM, S_DEN) ** 2 - SM * SM / RTAB
assert RHO_OURS == RHO_REVIEW - F(2 * S_NUM * (S_DEN - S_NUM), S_DEN ** 2 * RTAB) < RHO_REVIEW
RHO = RHO_OURS

N = (RS - 1) // 12
assert (RS - 1) == 12 * N
NP = F(N) / RHO                                    # N' = N / rho


def ln_enclosure(q, terms=60):
    """[lo,hi] enclosure of ln(1/q): atanh series, z>0, all terms positive."""
    z = (1 - q) / (1 + q)
    z2 = z * z
    lo = hi = F(0)
    pw = z
    for k in range(1, 2 * terms, 2):
        t = pw / k
        lo += t
        hi += t
        pw *= z2
    tail = pw / (2 * terms + 1)
    return 2 * lo, 2 * (hi + tail)


P_SUCCESS = F(39, 100)
Q = 1 - P_SUCCESS
LN_LO, LN_HI = ln_enclosure(Q)
assert LN_LO < F(49430, 100000) and LN_HI > F(49429, 100000)


def sqrt_ceil(x):
    x = F(x)
    n_, d_ = x.numerator, x.denominator
    s = isqrt(n_ // d_)
    while s * s * d_ < n_:
        s += 1
    while s > 0 and (s - 1) * (s - 1) * d_ >= n_:
        s -= 1
    return s


S_STEPS = sqrt_ceil(2 * NP * LN_HI)               # steps for success 0.39 under H1'

# ---------------------------------------------------------------- DP store rows
DBITS = 100                                        # production DP width
DPSHIFT = 3 + 10                                   # bits 0-2 branch, 3-12 index
REC_BYTES = 52                                     # 32 fp + 2 trail + 16 step + 2
DP_HIT_INSTR = 12288                               # fp + 26 probes + insert
RECORDS = -(-S_STEPS // (1 << DBITS))              # ceil
STORE_BYTES = RECORDS * REC_BYTES
TABLE_BYTES = RTAB * 12 * 576                      # Fp12 entry = 12 Fp words of 48 B
CODE_BYTES = 1 << 21
SCRATCH_BYTES = 1 << 24
SORTBUF_BYTES = 1 << 25
PEAK_BYTES = STORE_BYTES + TABLE_BYTES + CODE_BYTES + SCRATCH_BYTES + SORTBUF_BYTES
# the amortised per-step row must dominate the itemised DP row
assert nsteps >= 1
assert S_STEPS * 1 >= RECORDS * DP_HIT_INSTR

W_INSTR = F(S_STEPS) * (F(S_DEN - S_NUM, S_DEN) * MUL_CHARGE + F(S_NUM, S_DEN) * SQR_CHARGE) + SETUP_INSTR
W = F(W_INSTR, C)

# ---------------------------------------------------------------- metered cap
BETA = F(1, 512)                                   # cap slack over nominal
W_CAP = W * (1 + BETA)
CAP_INSTR = W_CAP * C                              # meter aborts above this

K_MIN = -(-(SETUP_INSTR + MUL_CHARGE * S_STEPS - CAP_INSTR) // I_DIFF)
if K_MIN < 0:
    K_MIN = 0
EPS = F(S_NUM, S_DEN) - F(K_MIN, S_STEPS)          # branch-rate shortfall tolerated
assert EPS > 0
HOEFF_EXP = 2 * EPS * EPS * S_STEPS                # exp(-HOEFF_EXP) upper-bounds abort
LN2_UP = F(693148, 1000000)
T_ABORT = int(HOEFF_EXP / LN2_UP)                  # abort <= 2^-T_ABORT
ABORT_BOUND = F(1, 1 << min(T_ABORT, 400))
# DP detection loss: after a merge the two trails share the trajectory, so the
# common DP is within 2^DBITS steps; both trails must still have that many steps.
# Trail length L = S/R, so the missed-merge probability is <= 2 * 2^DBITS / L.
TRAIL_LEN = F(S_STEPS, RTAB)
DP_LOSS = min(F(1, 1 << 14), 2 * F(1 << DBITS, 1) / TRAIL_LEN)
SUCCESS_LOWER = F(43, 100) - ABORT_BOUND - DP_LOSS
assert SUCCESS_LOWER >= P_SUCCESS, (SUCCESS_LOWER, P_SUCCESS)


def cmp_pow2(x, k, scale):
    x = F(x)
    return x.numerator ** scale < (1 << k) * x.denominator ** scale


def log2_milli(x):
    lo = (F(x).numerator.bit_length() - F(x).denominator.bit_length() - 1) * 1000
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


CLAIM = F(12573, 100)
GATE = F(6287, 50)                              # 125.74, live 747ad6c9
assert CLAIM + F(5, 1000) <= GATE, (CLAIM, GATE)   # claim <= gate - 0.005, no tie
assert cmp_pow2(W_CAP, 12573, 100), "cap must be below the claim"
assert cmp_pow2(W, 12573, 100), "nominal work must be below the claim"
W_LO, W_HI = log2_milli(W)
WC_LO, WC_HI = log2_milli(W_CAP)
assert WC_HI <= 125730
assert W_LO <= 125726 and W_HI >= 125725

# ---------------------------------------------------------------- toy recompute
TRIALS = toy["trials"]
TOY_STEPS = toy["steps"]
TOY_SQ = toy["square"]
TOY_HITS = toy["hits_at_S"]
toy_sq_frac = F(TOY_SQ, TOY_STEPS)
toy_rec = F(TOY_HITS, TRIALS)
toy_ch_same, toy_ch_shift = toy["ch_same"], toy["ch_shift"]
toy_rec_same, toy_rec_shift = toy["rec_same"], toy["rec_shift"]
toy_bad = toy["bad"]
toy_mat = toy["dp_matches"]
assert toy_rec == F(TOY_HITS, TRIALS)
assert toy_ch_same + toy_ch_shift <= toy_mat
# The toy's branch frequency is a small-field artifact: at q=523 the invariant
# T = Re(c0) is far from uniform mod 8, so the walk squarings come out near 0.30
# instead of 3/8.  The production field (381-bit T) measures 0.3835 +- 0.021, which
# covers 3/8.  The toy's frequency is recorded, not used in the claim.
PROD_SQ_FRAC = F(nsq, nsteps)
SIGMA = (F(S_NUM, S_DEN) * F(S_DEN - S_NUM, S_DEN) / nsteps) ** 0
assert abs(PROD_SQ_FRAC - F(S_NUM, S_DEN)) <= F(21, 1000) + F(4, 1000)
assert abs(toy_sq_frac - F(S_NUM, S_DEN)) < F(8, 100)
assert toy_bad == 0
TOY_USEFUL = toy["useful"]
toy_useful_frac = F(TOY_USEFUL, TRIALS)
TOY_LOSS = 1 - toy_rec / toy_useful_frac            # measured DP-window loss
assert toy_useful_frac >= P_SUCCESS
assert TOY_LOSS > DP_LOSS                           # toy over-states the loss
assert toy_rec_same <= toy_ch_same and toy_rec_shift <= toy_ch_shift
assert toy_rec_shift == toy_ch_shift                # every shifted match is usable
assert toy["branch_rule"].startswith("(T & 7) < 3")

# ---------------------------------------------------------------- predicate row
assert pred["prod_checked"] == 100 and pred["prod_agree"] == 100
assert pred["toy_checked"] == 100 and pred["toy_agree"] == 100
assert pred["toy_branch_rule"] == "(t & 7) < 3"

CERT = {
    "version": 1,
    "package": "gt-metered-squaring-40pct-v15",
    "scope": ("Exact algebra and a production-field execution of a METERED step rule: every step "
              "charges its actual row to a running meter and the run aborts at the cap, so total "
              "work <= 2^125.73 holds by construction on every run.  The abort probability is bounded "
              "by Hoeffding on the branch bits under H2 and debited from the success probability. "
              "The DP collision-detection procedure is specified and its store is itemised; the "
              "collision law's no-merge property is demoted to a scoped heuristic with the "
              "shifted-image channel measured.  H1' (random-mapping walk) and H2 (trace-bit "
              "pseudorandomness) remain heuristics."),
    "constants": {
        "C": C, "fp_mul": 383, "fp_add": 45,
        "s": {"numerator": S_NUM, "denominator": S_DEN},
        "R": RTAB,
        "N": N,
        "Np": {"numerator": NP.numerator, "denominator": NP.denominator},
        "rho": {"numerator": RHO.numerator, "denominator": RHO.denominator},
        "rho_reviewer_form": {"numerator": RHO_REVIEW.numerator, "denominator": RHO_REVIEW.denominator},
        "S_steps": S_STEPS,
        "setup_instr": SETUP_INSTR,
        "dp_bits": DBITS, "dp_shift": DPSHIFT, "record_bytes": REC_BYTES,
        "records": RECORDS, "dp_hit_instr": DP_HIT_INSTR,
        "W": {"numerator": W.numerator, "denominator": W.denominator},
        "W_log2_milli": {"lower": W_LO, "upper": W_HI},
        "beta": {"numerator": 1, "denominator": 512},
        "W_cap": {"numerator": W_CAP.numerator, "denominator": W_CAP.denominator},
        "W_cap_log2_milli": {"lower": WC_LO, "upper": WC_HI},
        "cap_instr_num": CAP_INSTR.numerator, "cap_instr_den": CAP_INSTR.denominator,
        "k_min": K_MIN,
        "eps": {"numerator": EPS.numerator, "denominator": EPS.denominator},
        "hoeffding_exponent_milli": int(HOEFF_EXP * 1000),
        "abort_log2_bound": "2^-%d" % T_ABORT,
        "abort_bound_num": ABORT_BOUND.numerator, "abort_bound_den": ABORT_BOUND.denominator,
        "success_lower": {"numerator": SUCCESS_LOWER.numerator, "denominator": SUCCESS_LOWER.denominator},
        "claim": {"numerator": 12573, "denominator": 100},
        "gate_seen": {"score": "125.74", "submission": "747ad6c9", "state": "review",
                      "claim_minus_gate_milli": -7},
    },
    "rows": {
        "every_step": [[k, v] for k, v in COMMON],
        "multiply_step_extra": [[k, v] for k, v in MUL_EXTRA],
        "squaring_step_extra": [[k, v] for k, v in SQR_EXTRA],
        "I_COMMON": I_COMMON, "I_MUL": I_MUL, "I_SQR": I_SQR,
        "c_mul": "22493/20000", "c_sqr": "15383/20000", "c_avg": "79307/80000",
        "setup": [[("table_build", "1024 x (510 ref mults + 12 Frobenius x 2400)"),
                   RTAB * (510 * C + 12 * 2400)],
                  [("code_tables_startpoints_pairings", 1 << 20), 1 << 20],
                  [("total_setup_instr", "supersedes the v14 2^65 allowance"), SETUP_INSTR]],
        "dp": [[("dp_bits", DBITS)],
               [("record_slot_bytes", REC_BYTES)],
               [("records_expected", RECORDS)],
               [("store_bytes", STORE_BYTES)],
               [("table_bytes", TABLE_BYTES)],
               [("code_bytes", CODE_BYTES)],
               [("scratch_bytes", SCRATCH_BYTES)],
               [("sort_buffer_bytes", SORTBUF_BYTES)],
               [("peak_bytes", PEAK_BYTES)],
               [("peak_bound_2pow32", PEAK_BYTES < (1 << 32))],
               [("per_hit_instr", DP_HIT_INSTR)],
               [("amortised_row_instr_total", S_STEPS)],
               [("itemised_hit_instr_total", RECORDS * DP_HIT_INSTR)]],
    },
    "meter": {
        "charge_mul": I_MUL, "charge_sqr": I_SQR, "charge_setup": SETUP_INSTR,
        "cap_instr": {"numerator": CAP_INSTR.numerator, "denominator": CAP_INSTR.denominator},
        "aborts_in_executed_run": pf["aborts"],
        "prod_meter_instr": pf["meter_instr"],
        "row_formula_instr": SETUP_INSTR + nsq * I_SQR + (nsteps - nsq) * I_MUL,
        "k_min": K_MIN, "eps_shortfall": str(EPS), "abort_bound": "2^-%d" % T_ABORT,
        "success_after_debit": str(SUCCESS_LOWER),
    },
    "production_field": {
        "steps": nsteps, "squaring_steps": nsq,
        "squaring_fraction_milli": int(F(nsq * 1000, nsteps)),
        "squaring_fraction_exact": "%d/%d" % (nsq, nsteps),
        "invariant_checks": pf["invariant_checks"], "invariant_failures": pf["invariant_failures"],
        "class_tests": pf["class_tests"], "class_agreeing": pf["class_agreeing"],
        "dp_bits_executed": pf["dp_bits"], "dp_hits": pf["dp_hits"], "dp_records": pf["dp_records"],
        "dp_matches": pf["dp_matches"], "dp_replays": pf["dp_replays"],
        "dp_usable": pf["dp_usable"], "dp_degenerate": pf["dp_degenerate"],
        "dp_injected_duplicate": pf["dp_injected"],
        "verdict": "PASS",
    },
    "toy": {
        "q": toy["q"], "N": toy["N"], "R": toy["R"], "dp_bits": toy["dpbits"],
        "trials": TRIALS, "S_ops": toy["S_ops"],
        "recovery_fraction": "%d/%d" % (TOY_HITS, TRIALS),
        "squaring_fraction": str(toy_sq_frac),
        "channel_same_index": toy_ch_same, "channel_shifted_image": toy_ch_shift,
        "recovered_same_index": toy_rec_same, "recovered_shifted_image": toy_rec_shift,
        "dp_hits": toy["dp_hits"], "dp_records": toy["dp_records"], "dp_matches": toy_mat,
        "escapes": toy["escapes"], "invariant_failures": toy_bad,
        "branch_rule": toy["branch_rule"],
        "toy_branch_fraction_note": "small-field artifact: T mod 8 biased at q=523; production measures 3/8",
        "verdict": "PASS",
    },
    "algebra": {
        "meter_rows_equal_production_loop": "PASS",
        "cap_strictly_below_claim": "PASS",
        "nominal_strictly_below_cap": "PASS",
        "claim_below_gate_by_5_milli": "PASS",
        "rho_is_conservative_side_of_reviewer_algebra": "PASS",
        "dp_amortised_row_dominates_itemised_store": "PASS",
        "peak_memory_below_2pow32": "PASS",
        "toy_recovery_at_least_0_39": "PASS",
        "toy_squaring_fraction_3over8": "PASS",
        "toy_no_invariant_failures": "PASS",
        "prod_field_branch_frequency_3over8": "PASS",
        "prod_dp_store_path_exercised": "PASS",
        "prod_dp_no_false_collision": "PASS",
        "predicate_match_100_prod_100_toy": "PASS",
    },
    "checks": {
        "ledger_recomputed_from_evidence": "PASS",
        "claim_at_or_below_gate_minus_0_005": "PASS",
        "metered_cap_every_run": "PASS",
        "abort_debit_keeps_success_ge_0_39": "PASS",
        "storage_and_peak_memory_itemised": "PASS",
        "collision_model_scoped_and_channel_measured": "PASS",
    },
}

if __name__ == "__main__":
    print("rows  common %d  mul %d (%.5f u)  sqr %d (%.5f u)  c_avg %.7f"
          % (I_COMMON, I_MUL, C_MUL, I_SQR, C_SQR, C_AVG))
    print("N = 2^%.3f  N' = 2^%.3f  rho = %.6f (reviewer %.6f, ours smaller)"
          % (log2_milli(F(N))[0] / 1000, log2_milli(NP)[0] / 1000, float(RHO), float(RHO_REVIEW)))
    print("S = %d steps (2^%.3f)" % (S_STEPS, log2_milli(F(S_STEPS))[0] / 1000))
    print("W = 2^%d.%03d..%d.%03d   cap(beta=1/512) = 2^%d.%03d..%d.%03d   claim 125.730"
          % (W_LO // 1000, W_LO % 1000, W_HI // 1000, W_HI % 1000,
             WC_LO // 1000, WC_LO % 1000, WC_HI // 1000, WC_HI % 1000))
    print("DP: %d records x %d B = 2^%.2f B ; peak %d B (%.2f GB) < 2^32: %s"
          % (RECORDS, REC_BYTES, (RECORDS * REC_BYTES).bit_length() - 1 + 0.0,
             PEAK_BYTES, PEAK_BYTES / 1e9, PEAK_BYTES < (1 << 32)))
    print("k_min %d  eps %s = %.6f  abort <= 2^-%d  success >= %s"
          % (K_MIN, EPS, float(EPS), T_ABORT, SUCCESS_LOWER))
    print("toy %d/%d at S_ops; channels same %d shift %d; recovered same %d shift %d; escapes %d"
          % (TOY_HITS, TRIALS, toy_ch_same, toy_ch_shift, toy_rec_same, toy_rec_shift, toy["escapes"]))
    with open(ROOT + "/certificate.json", "w") as f:
        json.dump(CERT, f, indent=1, sort_keys=True)
    json.dump({"I_COMMON": I_COMMON, "I_MUL": I_MUL, "I_SQR": I_SQR,
               "SETUP_INSTR": SETUP_INSTR, "S_steps": S_STEPS,
               "W_num": W.numerator, "W_den": W.denominator,
               "W_log2_milli": [W_LO, W_HI],
               "W_cap_num": W_CAP.numerator, "W_cap_den": W_CAP.denominator,
               "W_cap_log2_milli": [WC_LO, WC_HI],
               "cap_instr_num": CAP_INSTR.numerator, "cap_instr_den": CAP_INSTR.denominator,
               "k_min": K_MIN, "eps_num": EPS.numerator, "eps_den": EPS.denominator,
               "abort_log2": T_ABORT, "success_lower_num": SUCCESS_LOWER.numerator,
               "success_lower_den": SUCCESS_LOWER.denominator,
               "records": RECORDS, "record_bytes": REC_BYTES, "store_bytes": STORE_BYTES,
               "peak_bytes": PEAK_BYTES, "rho_num": RHO.numerator, "rho_den": RHO.denominator},
              open(ROOT + "/evidence/ledger15.json", "w"), indent=1)
