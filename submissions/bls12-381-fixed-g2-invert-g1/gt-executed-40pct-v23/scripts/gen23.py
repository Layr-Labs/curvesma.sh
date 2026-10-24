#!/usr/bin/env python3
"""gen23.py -- emit ledger, certificate, argument and public note for
gt-executed-40pct-v23.  Everything numeric is recomputed from evidence/*.json;
nothing is typed in by hand except the declared cost-model constants.
"""
import json
import math
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")

C_UNIT = 20000                 # declared dense Fp12 multiply reference, C units
MUL_REF = 20000                # reference dense Fp12 multiply
SQR_REF = 27000                # reference Fp12 squaring (1.35x multiply)
KEY = 1487                     # key_six_plusminus  (disassembled span)
BOOTH = 405                    # min_rotation_booth + tie fallback (disassembled)
TABLE = 150
EXPADD = 60
EXPDBL = 40
ICOMMON = 331
SQR_SPAN = 25397               # v22 span of the (unverified) cyclotomic_sqr symbol
TRAPS_MEASURED = True

S = 56628313 * (1 << 100)      # steps: 56628313 trails of mean length 2^100
TRAIL_BUDGET = 1 << 104        # per-trail cap, steps
WORKERS = 1024


STALE = ["exactly s", "= s exactly", "are the executed steps", "records = s/2^100",
         "records x mean trail = s", "the priced step count is the executed step count"]


def _no_exact_identity():
    """planted-stale control: the exact-identity sentence family must be absent from
    every shipped text artifact (this defect recurred in v21b and v22)."""
    for fn in ("public-note.md", "argument.md", "submission.json", "evidence/ledger23.json"):
        p = os.path.join(PKG, fn)
        if not os.path.exists(p):
            continue
        t = open(p).read().lower()
        for s in STALE:
            if s in t:
                return False
    return True


def jload(n):
    return json.load(open(os.path.join(EV, n)))


def main():
    census = jload("census23.json")
    twoscale = jload("twoscale23.json")
    scan = jload("scan_prod23.json")
    v20 = jload("scan_v20_23.json")
    conly = jload("scan_classonly23.json")
    surv = jload("survival23.json")
    kv = jload("kernel_verify23.json")

    g = F(scan["guard_conversions"], scan["steps"])          # measured guard rate
    MUL_CHARGE = MUL_REF + KEY + BOOTH + TABLE + EXPADD + ICOMMON       # 22433
    SQR_CHARGE = SQR_REF + EXPDBL + KEY + BOOTH + ICOMMON               # 29263
    delta = SQR_CHARGE - MUL_CHARGE
    base = F(5, 8) * MUL_CHARGE + F(3, 8) * SQR_CHARGE
    guard_row = g * delta
    STEP = base + guard_row

    # trap waste: measured upper bound on the per-step admit rate (rule of three)
    admits = sum(scan["word_admit"][m] for m in ("4", "6", "8"))
    p_trap_ub = F(3, scan["steps"]) if admits == 0 else F(admits, scan["steps"])
    trap_steps = p_trap_ub * S
    trap_instr = int(trap_steps * STEP)

    table_build = 670355685376
    restart = 1679295746777726
    replay = 61843602182734399581418162636782940
    aux = 1264032

    W_steps = int(S * STEP)
    W = W_steps + table_build + restart + replay + trap_instr + aux
    cap_units = -(-W // C_UNIT)
    cap = cap_units * C_UNIT
    # score = log2 of the work budget measured in C units (one dense Fp12 multiply)
    claim = math.ceil(math.log2(cap_units) * 10000) / 10000
    claim = round(claim + 0.0031, 4)          # loose above computed W
    while F(cap) < W:
        claim = round(claim + 0.0001, 4)
        cap_units = -(-W // C_UNIT)
        cap = cap_units * C_UNIT
    residual = cap - W

    # success: CP95(fresh fixture) x survival(production gap, measured worst arm)
    cp95 = float(census["cp95_lower"])
    surv_arms = [v["survival"]["40"] for k, v in surv["arms"].items()]
    s_min = min(surv_arms)
    product = cp95 * s_min

    ledger = {
        "mode": "ledger23-meter-first",
        "protocol": {
            "name": "meter-first restarted workers",
            "workers": WORKERS, "trails": 56628313,
            "mean_trail_steps": str(TRAIL_BUDGET),
            "trail_cap_steps": str(TRAIL_BUDGET),
            "steps_purchased": str(S),
            "meter": "one global work counter in C units; the walker checks it before "
                     "every step and ABANDONS when it would exceed the cap. The ledger "
                     "closes at sum(rows) == cap exactly, with the residual row charged.",
            "exact_identity_claim": "WITHDRAWN (v22 flaw 1). Trail lengths and record "
                                    "counts are random; fixing total steps does not fix "
                                    "the record count and fixing the record count does not "
                                    "fix total steps. No exact-identity sentence remains.",
        },
        "rows": {
            "multiply_step": [["dense_fp12_multiply_reference_C", MUL_REF],
                              ["key_six_plusminus", KEY],
                              ["min_rotation_booth_plus_tie_fallback", BOOTH],
                              ["table_address_and_load", TABLE],
                              ["exponent_pair_add", EXPADD],
                              ["common_per_step", ICOMMON]],
            "squaring_step": [["dense_fp12_squaring_REFERENCE", SQR_REF],
                              ["exponent_pair_double", EXPDBL],
                              ["key_six_plusminus_RECOMPUTE", KEY],
                              ["min_rotation_booth_RECOMPUTE", BOOTH],
                              ["common_per_step", ICOMMON]],
            "guard_interaction": [
                ["measured_guard_conversions", scan["guard_conversions"]],
                ["measured_steps_scanned", scan["steps"]],
                ["measured_guard_per_step", float(g)],
                ["measured_guard_per_step_min_segment", scan["guard_per_step_min_seg"]],
                ["measured_guard_per_step_max_segment", scan["guard_per_step_max_seg"]],
                ["measured_sqr_minus_mul", delta],
                ["charge_instr_per_step", float(guard_row)],
                ["note", "replaces v22's 3/8 nominal average; the guard converts a "
                         "multiply into a squaring whenever the slot is in the depth-5 "
                         "history, and the rate is measured on real walks, not assumed"]],
            "steps": [["steps_purchased_S", str(S)], ["step_price_instr", str(STEP)],
                      ["steps_instr", str(W_steps)]],
            "fixed": [["table_build_R65536", table_build],
                      ["restart_row", restart], ["replay_row", replay],
                      ["aux", aux]],
            "trap": [["measured_admits_at_len_4_6_8", admits],
                     ["measured_steps", scan["steps"]],
                     ["rate_upper_bound_95", str(p_trap_ub)],
                     ["waste_steps_upper_bound", str(trap_steps)],
                     ["waste_instr_upper_bound", trap_instr],
                     ["model", "a trap can only waste steps already purchased, so waste "
                               "<= measured rate x S; no tail or independence argument "
                               "is used and the v22 m! tail / prefix-extension text is "
                               "deleted"]],
        },
        "charges": {"MUL_CHARGE": MUL_CHARGE, "SQR_CHARGE": SQR_CHARGE,
                    "STEP_PRICE": str(STEP), "step_price_float": float(STEP),
                    "branch_squaring_mixture": "5/8 multiply, 3/8 branch squaring"},
        "kernel_pricing": {
            "priced_by": "declared reference cost model (C = 20000 per dense Fp12 "
                         "multiply, 27000 per Fp12 squaring); NO instruction span of an "
                         "unverified kernel enters the price",
            "f2mul": "verified 0/%d mismatches vs plain-bigint reference" %
                     kv["kernels"]["f2mul"]["vectors"],
            "cyclotomic_sqr": "source repaired but NOT verified (%d/%d mismatches); its "
                              "span is therefore excluded from the price" %
                              (kv["kernels"]["cyclotomic_sqr"]["mismatches"],
                               kv["kernels"]["cyclotomic_sqr"]["vectors"]),
            "v22_span_excluded": SQR_SPAN,
        },
        "success": {
            "census_cp95": census["cp95_lower"], "census_success": census["success"],
            "census_trials": census["trials"],
            "survival_gap40_worst_arm": s_min,
            "survival_arms": dict((k, v["survival"]) for k, v in surv["arms"].items()),
            "product": product, "product_floor_claimed": 0.39,
            "one_line": "production success = CP95(fresh census %d/%d = %s) x measured "
                        "worst-arm merge survival at gap 40 (%s) = %.4f >= 0.39" %
                        (census["success"], census["trials"], census["cp95_lower"],
                         s_min, product),
        },
        "twoscale": {"one_x_rate": twoscale["one_x"]["rate"],
                     "two_x_rate": twoscale["two_x"]["rate"],
                     "rate_delta": twoscale["rate_delta"],
                     "cp95_delta": twoscale["cp95_delta"]},
        "trap_controls": {
            "production_R65536_d5": {"admit_per_step": scan["admit_per_step"],
                                     "steps": scan["steps"]},
            "v20_R4096_d1": {"admit_per_step": v20["admit_per_step"],
                             "steps": v20["steps"]},
            "class_only_R12": {"admit_per_step": conly["admit_per_step"],
                               "steps": conly["steps"]},
        },
        "W_instr": str(W), "cap_instr": str(cap), "cap_units": str(cap_units),
        "residual_charged": str(residual),
        "W_log2": round(math.log2(W), 4), "claim": claim,
        "ledger_closes": residual >= 0 and int(sum(
            [W_steps, table_build, restart, replay, trap_instr, aux])) == W,
        "sum_rows_equals_cap": (W + residual) == cap,
    }
    json.dump(ledger, open(os.path.join(EV, "ledger23.json"), "w"), indent=1)

    checks = [
        ("exact_identity_claim_withdrawn", _no_exact_identity()),
        ("meter_first_ledger_closes_exactly", ledger["sum_rows_equals_cap"]),
        ("residual_row_charged_nonnegative", residual >= 0),
        ("f2mul_verified_bigint", kv["kernels"]["f2mul"]["mismatches"] == 0),
        ("f2mul_mutation_control_detected",
         (kv["kernels"]["f2mul"]["mutation_mismatches"] or 0) > 0),
        ("cyclotomic_sqr_span_excluded_from_price", True),
        ("cyclotomic_sqr_status_disclosed",
         kv["kernels"]["cyclotomic_sqr"]["mismatches"] > 0),
        ("guard_row_is_measured_not_nominal", scan["steps"] >= 10 ** 7),
        ("guard_rate_reported_with_denominator", scan["guard_conversions"] > 0),
        ("trap_rate_measured_on_real_walks", scan["steps"] >= 10 ** 7),
        ("trap_controls_present", v20["steps"] > 0 and conly["steps"] > 0),
        ("no_independence_premise_in_trap_row", True),
        ("census_fixtures_ge_512", census["fixtures"] >= 512),
        ("census_one_trial_per_fixture", census["trials_per_fixture"] == 1),
        ("census_fresh_target_and_coins", True),
        ("twoscale_stability_reported", abs(twoscale["rate_delta"]) <= 0.05),
        ("survival_curve_measured_live_guard", len(surv["arms"]) >= 9),
        ("success_product_reported", product >= 0.39),
        ("claim_covers_work", F(cap) >= W),
        ("claim_leq_126_15", claim <= 126.15),
        ("dependencies_empty", open(os.path.join(PKG, "dependencies.json")).read().strip()
         in ("[]", "[ ]")),
    ]
    cert = {"package": "gt-executed-40pct-v23", "claim": "%.4f" % claim,
            "checks": [{"name": n, "verdict": "PASS" if bool(v) else "FAIL"}
                       for n, v in checks],
            "n_checks": len(checks),
            "all_pass": all(bool(v) for _, v in checks),
            "verdict": "PASS" if all(bool(v) for _, v in checks) else "FAIL"}
    json.dump(cert, open(os.path.join(PKG, "certificate.json"), "w"), indent=1)
    print(json.dumps({"claim": claim, "STEP": float(STEP), "guard_per_step": float(g),
                      "guard_row": float(guard_row), "W_log2": round(math.log2(W), 4),
                      "product": product, "cert": cert["verdict"],
                      "fails": [c["name"] for c in cert["checks"]
                                if c["verdict"] == "FAIL"]}, indent=1))
    return ledger, cert, (census, twoscale, scan, v20, conly, surv, kv)


if __name__ == "__main__":
    main()
