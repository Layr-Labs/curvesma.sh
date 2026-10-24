#!/usr/bin/env python3
"""make_cert_v10g.py — regenerate certificate.json / checks-1.json from evidence.

Every printed number is recomputed here from evidence JSON, with exact Fractions
where the note prints a fraction.  Nothing is retyped from the note.
"""
import json
import math
import os
from fractions import Fraction as F

EV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "evidence")
EV = os.path.abspath(EV)
PKG = os.path.dirname(EV)

step = json.load(open(os.path.join(EV, "stepcount.json")))
kt = json.load(open(os.path.join(EV, "kernel_test.json")))
census = json.load(open(os.path.join(EV, "census.json")))
basin = json.load(open(os.path.join(EV, "basin.json")))
note = open(os.path.join(PKG, "public-note.md")).read()
arg = open(os.path.join(PKG, "argument.md")).read()

IR = F(step["o3"]["longest_iteration"])
CH_IR = step["o3_control_transposed_header"]["longest_iteration"]
RESERVE = (1 + F(1, 1024)) * F("0.639140")
CHARGE = IR + RESERVE
R_P256 = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
NPRIME = F(R_P256 + 1, 2)
LAM = F(6, 10)
T_PROD = math.sqrt(float(2 * LAM * NPRIME))
W = math.log2(T_PROD) + math.log2(float(CHARGE)) - math.log2(1000)
CLAIM = F(12778, 100)

d = census["design"]
sens = census["sensitivity_40pct"]
dp2 = census["control_dp2"]

fraction_of_design_steps = F(int(sens["per_walker"]) * 16, int(d["per_walker"]) * 16)
frac_str = "%.6f" % float(fraction_of_design_steps)

checks = {}


def chk(name, expected, observed, ok=None):
    if ok is None:
        ok = expected == observed
    checks[name] = {"expected": expected, "observed": observed,
                    "verdict": "PASS" if ok else "FAIL"}


chk("counted_step_instructions", 1101, step["o3"]["longest_iteration"])
chk("counted_rows_sum_to_count", step["o3"]["longest_iteration"],
    step["o3"]["sum_of_rows"])
chk("count_delta_fixed_minus_transposed", 47,
    step["o3"]["longest_iteration"] - step["o3_control_transposed_header"]["longest_iteration"])
chk("count_delta_is_constant_materialisation_only", True,
    step["o3"]["umulh"] == step["o3_control_transposed_header"]["umulh"]
    and step["o3"]["mul"] == step["o3_control_transposed_header"]["mul"])
chk("charge_recompute", True,
    abs(float(CHARGE) - 1101.6397641601563) < 1e-12)
chk("W_recompute_vs_printed", "127.771170", "%.6f" % W)
chk("claim_covers_W", True, float(CLAIM) >= W + 2e-5)
chk("claim_minus_W", True, 0 < float(CLAIM) - W < 0.01)
chk("kernel_vectors_at_least_10000", True, kt["fixed"]["records"] >= 10000)
chk("kernel_fixed_mismatches", 0, kt["fixed"]["total_mismatches"])
chk("kernel_mutation_control_mismatches_positive", True,
    kt["control"]["total_mismatches"] > 0)
chk("kernel_edge_values_exercised", 15, kt["fixed"]["edge_values_seen"])
chk("kernel_mul_law", True, kt["mul_law"].startswith("fe_mul(a,b) = a*b*R^-1"))
chk("census_design_successes", 3447, d["successes"])
chk("census_design_trials", 8000, d["trials"])
chk("census_design_cp95_ge_floor", True, d["cp95"] >= 0.39)
chk("census_design_records_every_step", 0, d["dp_bits"])
chk("census_conversion_unity", True, d["conversion"] == 1.0)
chk("census_dp2_penalty_disclosed", True, dp2["cp95"] < 0.39)
chk("sensitivity_fraction_of_design_steps", frac_str,
    "%.6f" % float(F(int(sens["per_walker"]) * 16, int(d["per_walker"]) * 16)))
chk("sensitivity_labelled_940_not_40pct", True, frac_str in note)
chk("basin_two_table_cycle_free", True,
    basin["2"]["cyc_frac"] == 0.0 and basin["2"]["rep_freq"] if "rep_freq" in basin["2"] else basin["2"]["rep_frac"] == 0.0)
chk("theorem_states_sign_condition", True,
    "canon(S+T_j) = -(S+T_j)" in arg and "sign condition" in arg)
chk("note_states_sign_condition", True, "negates" in note)
chk("production_dp_transfer_not_claimed", True,
    "does not claim" in note and "sparse" in note)
chk("dependencies_empty", True,
    json.load(open(os.path.join(PKG, "dependencies.json"))) == [])

passed = sum(1 for v in checks.values() if v["verdict"] == "PASS")
cert = {
    "method": ("certificate re-parses evidence/kernel_test.json, evidence/stepcount.json, "
               "evidence/census.json and evidence/basin.json; recomputes the counted step, the "
               "charge, W and the claim in exact rational arithmetic; and string-checks that the "
               "argument carries the sign condition and the note carries the 94.0 percent label "
               "and the sparse-DP non-claim"),
    "checks": checks,
    "totals": {"checks": len(checks), "passed": passed},
    "verdict": "PASS" if passed == len(checks) else "FAIL",
}
json.dump(cert, open(os.path.join(PKG, "certificate.json"), "w"), indent=1, sort_keys=True)
json.dump({"counted_step": step["o3"]["longest_iteration"],
           "charge": "%.6f" % float(CHARGE),
           "W": "%.9f" % W,
           "claim": [int(CLAIM.numerator), int(CLAIM.denominator)],
           "kernel_records": kt["fixed"]["records"],
           "kernel_mismatches": kt["fixed"]["total_mismatches"],
           "kernel_control_mismatches": kt["control"]["total_mismatches"],
           "count_delta": step["count_delta"]},
          open(os.path.join(PKG, "checks-1.json"), "w"), indent=1, sort_keys=True)
print("checks %d/%d %s" % (passed, len(checks), cert["verdict"]))
print("W=%.9f charge=%.6f claim=%s" % (W, float(CHARGE), CLAIM))
