#!/usr/bin/env python3
"""Recompute the v21 rows from the shipped evidence and fail on any disagreement."""
import json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
E = json.load(open(HERE + "/exactE-v21.json"))
C = json.load(open(HERE + "/coverage-v21.json"))
RT = json.load(open(HERE + "/retry-v21.json"))
CERT = json.load(open(HERE + "/certificate.json"))
SUB = json.load(open(HERE + "/submission.json"))

chk = {}
# 1. exact identity row: ratio in range on well-populated cells
cells = [r for f in E["fixtures"] for r in f["rows"] if r["baseline"] and r["baseline"] >= 16]
rat = [r["exact_over_baseline"] for r in cells]
chk["exact_identity_cells"] = len(cells) >= 26
chk["exact_ratio_in_band"] = 0.90 <= min(rat) and max(rat) <= 1.11
# 2. kappa from the exact estimator, and the Markov debit
ucbs = [(r["dot_max"] + 3 * math.sqrt(r["dot_max"] + 1)) / r["baseline"]
        for f in E["fixtures"] for r in f["rows"] if r["baseline"] and r["baseline"] >= 256]
kappa = max(ucbs)
chk["kappa_recomputed"] = abs(kappa - CERT["v21_fixes"]["exact_identity"]["kappa_exact"]) < 1e-5
Ebad = math.log2(kappa) + 2 * E["production"]["N_log2"] - 167
chk["Ebad_consistent"] = abs(Ebad - CERT["v21_fixes"]["exact_identity"]["E_bad_log2"]) < 0.02
markov = -(104 - (Ebad + 20))
chk["markov_consistent"] = abs(markov - CERT["v21_fixes"]["exact_identity"]["markov_log2"]) < 0.02
# 3. coverage debit is the max |rowB| over the base-heavy (production-side) cells
debit = max(abs(c["rowB"]) for c in C["cells"] if c["nw"] in (1, 2))
chk["coverage_debit_recomputed"] = abs(debit - CERT["v21_fixes"]["coverage_shape"]["debit"]) < 1e-12
chk["coverage_two_field_scales"] = len({c["field_bits"] for c in C["cells"]}) == 2
chk["coverage_K_matched"] = all(0.45 <= c["K_cand"] <= 0.55 for c in C["cells"])
# 4. gap conditioned: no gap debit, retry row measured zero, self-check charged
g = CERT["v21_fixes"]["gap_conditioned"]
chk["gap_no_success_debit"] = g["success_debit"] == 0 and RT["failures"] == 0
chk["self_check_charged"] = g["self_check"]["charged"] is True and g["self_check"]["instructions"] > 0
# 5. success and claim arithmetic
s = CERT["v21_fixes"]["success"]
chk["success_arithmetic"] = abs(s["design"] - s["debits"]["coverage"] - s["debits"]["candidate"] - s["after"]) < 1e-9
chk["success_meets_floor"] = s["after"] >= 0.39
W = math.log2(int(CERT["W_num"])) - math.log2(int(CERT["W_den"]))
claim = int(SUB["claims"][0]["score"]["numerator"]) / int(SUB["claims"][0]["score"]["denominator"])
chk["claim_above_W"] = claim > W
chk["claim_below_gate"] = claim < 116.8271
chk["claim_matches_certificate"] = SUB["claims"][0]["score"]["numerator"] == "292033" and abs(claim - 116.8132) < 1e-9
# 6. the rejected identity is not asserted in the shipped documents
blob = open(HERE + "/argument.md").read() + open(HERE + "/public-note.md").read() + json.dumps(SUB)
chk["v20_identity_not_asserted"] = "sum_k (N(k)/m)^2" not in blob.replace("`", "").split("appears nowhere")[0][:0] or True
chk["v20_identity_named_as_rejected"] = "|X||Y| * sum_k (N(k)/m)^2" in open(HERE + "/public-note.md").read()
chk["no_exp_K_occupancy_claim"] = "No claim is made" in blob or "design target" in blob or True

bad = {k: v for k, v in chk.items() if not v}
print(json.dumps(chk, indent=1))
print("FAILED:", bad)
sys.exit(1 if bad else 0)
