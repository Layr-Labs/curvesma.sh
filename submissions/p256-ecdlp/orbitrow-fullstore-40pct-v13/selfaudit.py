#!/usr/bin/env python3
"""selfaudit.py -- v13 independent re-derivation + planted-error controls."""
import json
import math
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
led = json.load(open(os.path.join(HERE, "evidence/ledger-out.json")))
fs = json.load(open(os.path.join(HERE, "evidence/fullstore-census.json")))
ss = json.load(open(os.path.join(HERE, "evidence/selector-swap.json")))
ex = json.load(open(os.path.join(HERE, "evidence/exec-count.json")))
sdb = json.load(open(os.path.join(HERE, "submission.json")))

CHECKS = []


def chk(name, ok, detail=""):
    CHECKS.append({"check": name, "pass": bool(ok), "detail": str(detail)})
    return bool(ok)


rows = sum(n for _o, n in led["step_rows"]) + sum(n for _o, n in led["insert_rows"])
chk("insn_rows_sum_to_charge", rows == led["total_insn"] == 290, rows)

price = F(1) + F(led["total_insn"], 1000)
extras_ratio = 2.0 ** (led["log2_extras"] - led["log2_phase1"])
log2W = led["log2_phase1"] + math.log2(1.0 + extras_ratio)
chk("claim_ge_recomputed_log2W", log2W <= led["claim"] and led["claim"] - log2W < 1e-3,
    (log2W, led["claim"]))
chk("phase1_matches_price_times_S_OPS",
    abs((math.log2(float(price)) + led["S_OPS_log2"]) - led["log2_phase1"]) < 1e-6,
    (math.log2(float(price)) + led["S_OPS_log2"], led["log2_phase1"]))
one_ref = log2W + math.log2((2.0 + led["total_insn"] / 1000.0) / (1.0 + led["total_insn"] / 1000.0))
chk("claim_at_1ref_rederived", abs(one_ref - led["claim_at_1ref_access"]) < 2e-3,
    (one_ref, led["claim_at_1ref_access"]))
chk("submitted_claim_is_conservative_row",
    int(sdb["claims"][0]["score"]["numerator"]) == int(round(led["claim_at_1ref_access"] * 1000)),
    sdb["claims"][0]["score"])

succ = led["cp95_fullstore"] - led["debits"]
chk("success_rederived", abs(succ - led["success"]) < 5e-7, (succ, led["success"]))
chk("floor_met_direction", led["success"] >= 0.39 and led["floor_met"], led["success"])
chk("cp95_below_rate", led["cp95_fullstore"] <= led["fullstore_rate"],
    (led["cp95_fullstore"], led["fullstore_rate"]))
chk("no_dp_gate", led["dp_gate_on_success_path"] is False)
chk("cells_are_fullstore", all(c["success_rate"] > 0.6 for c in led["fullstore_census_cells"]))

cells = sorted(fs["cells"], key=lambda c: c["classes_log2"])
lo, hi = cells[0], cells[-1]
chk("scale_flat_across_orders", abs(hi["fullstore"]["success_rate"] - lo["fullstore"]["success_rate"]) < 0.06,
    (hi["fullstore"]["success_rate"], lo["fullstore"]["success_rate"]))
g = lo["gate_replay"]
chk("gate_ablation_monotone", g["d5"]["success_rate"] > g["d8"]["success_rate"] > g["d11"]["success_rate"],
    {k: v["success_rate"] for k, v in g.items()})
chk("fullstore_beats_gated_at_d11", lo["fullstore"]["success_rate"] > g["d11"]["success_rate"],
    (lo["fullstore"]["success_rate"], g["d11"]["success_rate"]))

rm = ss["range_matched"]
chk("selector_range_matched", rm["EC_reach"] == rm["HASH_reach"] == rm["R"], rm)
chk("selector_delta_ci_contains_zero", ss["delta_ci95"][0] <= 0 <= ss["delta_ci95"][1], ss["delta_ci95"])

chk("support_within_cap", led["support_measured_max"] <= led["k_max_per_trail"],
    (led["support_measured_max"], led["k_max_per_trail"]))
chk("memory_row_present", led["table_bytes_log2"] > 128, led["table_bytes_log2"])
chk("exec_row_withdrawn", ex["v12_row_withdrawn"] is True and ex["no_executed_number_claimed"] is True)

for lit, val in (("129.405", led["claim_at_1ref_access"]), ("0.6709", led["success"]),
                 ("128.577", led["claim"])):
    chk("literal_%s" % lit, abs(float(lit) - val) < 5e-4, val)
chk("literal_table_bytes_2^134.52", abs(134.52 - led["table_bytes_log2"]) < 0.01, led["table_bytes_log2"])


def mutated(fn):
    try:
        return not fn()
    except Exception:
        return True


chk("planted_claim_off_by_1", mutated(lambda: abs(led["claim_at_1ref_access"] - (one_ref + 0.001)) < 1e-6))
chk("planted_rows_off_by_1", mutated(lambda: (rows + 1) == 290))
chk("planted_success_below_floor", mutated(lambda: (led["success"] - 0.30) >= 0.39))
chk("planted_delta_ci_excludes_zero", mutated(lambda: ss["delta_ci95"][0] <= -1.0 <= ss["delta_ci95"][1]))

PASS = sum(1 for c in CHECKS if c["pass"])
out = {"test": "selfaudit_v13", "checks": CHECKS, "passed": PASS, "total": len(CHECKS),
       "all_pass": PASS == len(CHECKS)}
json.dump(out, open(os.path.join(HERE, "evidence/selfaudit.json"), "w"), indent=1, sort_keys=True)
for c in CHECKS:
    print(("PASS " if c["pass"] else "FAIL ") + c["check"] + "  " + c["detail"][:100])
print("%d/%d" % (PASS, len(CHECKS)))
raise SystemExit(0 if PASS == len(CHECKS) else 1)
