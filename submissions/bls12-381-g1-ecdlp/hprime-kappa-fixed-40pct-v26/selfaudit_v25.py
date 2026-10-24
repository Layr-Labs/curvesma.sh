#!/usr/bin/env python3
"""v25 self-audit: certificate <-> documents consistency, CP direction assertions,
and planted-edit ± planted-flip controls (the audit must FAIL when either is
perturbed, otherwise it is vacuous)."""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(HERE + "/certificate.json"))
ARG = open(HERE + "/argument.md").read()
NOTE = open(HERE + "/public-note.md").read()
GATE = json.load(open(HERE + "/gate-v25.json"))
DIRECT = json.load(open(HERE + "/direct-v25.json"))
AUD = DIRECT["summary"]["cp_direction_audit"]


def check(arg, note, planted_claim=None):
    out = {}
    claim = planted_claim or (int(C["claim"]["numerator"]) / int(C["claim"]["denominator"]))
    out["claim_in_argument"] = ("%.4f" % claim) in arg
    out["claim_in_note"] = ("%.4f" % claim) in note
    out["kappa_in_argument"] = ("%.3f" % C["debit"]["kappa"]) in arg
    out["success_in_argument"] = ("%.6f" % C["success"]["after"]) in arg
    out["success_in_note"] = ("%.6f" % C["success"]["after"]) in note
    out["floor_flag_in_note"] = ("meets_floor = %s" % str(C["success"]["meets_floor"]).lower()) in note
    out["ledger_recomputes"] = abs(sum(C["rows"].values()) / 2000.0 - 2 ** C["raw_log2"]) / \
        (sum(C["rows"].values()) / 2000.0) < 1e-6
    out["claim_above_W"] = claim > math.log2(sum(C["rows"].values()) / 2000.0)
    out["claim_below_gate"] = claim < 116.8271
    out["rows_11"] = len(C["rows"]) == 11
    out["gate_passed"] = bool(GATE["passed"] and GATE["ok"] == GATE["n"])
    out["cp_direction_rows_ok"] = all(r["ok"] and r["side_ok"] for r in AUD["rows"])
    out["cp_endpoints_bracket_q"] = all(r["endpoint"] <= r["point"] if r["direction"] == "LOWER"
                                        else r["endpoint"] >= r["point"] for r in AUD["rows"])
    out["cp_planted_flip_caught"] = AUD["planted_flip_all_caught"]
    out["pair_correction_in_every_cell"] = all(
        abs(g["pair_correction"] - (1 - (g["nw"] - 1) / g["nw"] ** 2)) < 1e-12
        for g in DIRECT["groups"])
    out["fixed_W_series"] = sorted(g["nw"] for g in DIRECT["groups"]) == [1, 1, 2, 2, 4, 4]
    out["same_records_per_side"] = len({g["records_per_side"] for g in DIRECT["groups"]}) == 1
    out["success_meets_floor"] = bool(C["success"]["meets_floor"] ==
                                      (C["success"]["after"] >= 0.39))
    out["kappa_slope_row"] = bool("slope_consistency" in C["kappa_168"])
    out["kappa_pinned_fit_rejected"] = bool(C["kappa_168"]["pinned_fit_rejected"]
                                            ["admissibility"].startswith("not admitted"))
    out["kappa_charged_above_band"] = bool(C["kappa_168"]["charged_exceeds_band"])
    out["all_certificate_checks"] = all(C["checks"].values())
    out["PASS"] = all(v for k, v in out.items() if k != "PASS" and isinstance(v, bool))
    return out


base = check(ARG, NOTE)
planted = check(ARG, NOTE, planted_claim=int(C["claim"]["numerator"]) /
                int(C["claim"]["denominator"]) + 1e-3)
for k in ("claim_in_argument", "claim_in_note"):
    base.pop(k)
base["planted_claim_edit_detected"] = not (planted["claim_in_argument"] and planted["claim_in_note"])
base["planted_flip_value"] = "endpoint side swapped in every CP row"
base["evidence_files"] = len([f for f in os.listdir(HERE) if not f.startswith(".")
                              and os.path.isfile(os.path.join(HERE, f))])
base["PASS"] = base["PASS"] and base["planted_claim_edit_detected"] and \
    base["cp_planted_flip_caught"]
json.dump(base, open(HERE + "/selfaudit-v25.json", "w"), indent=1)
print(json.dumps(base, indent=1))
sys.exit(0 if base["PASS"] else 1)
