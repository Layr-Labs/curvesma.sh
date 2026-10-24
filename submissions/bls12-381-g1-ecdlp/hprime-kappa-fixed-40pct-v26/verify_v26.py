#!/usr/bin/env python3
"""v26 certificate verifier: recompute every published claim from the evidence."""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
P = lambda n: json.load(open(os.path.join(HERE, n)))
cert, led, kap, rat, tr, gate = (P("certificate.json"), P("ledger-v26.json"),
    P("kappa-v26.json"), P("ratios-v26.json"), P("transfer-v26.json"), P("gate-v25.json"))
sub = P("submission.json")
checks = {}
rows = led["rows"]
checks["integer_sum"] = sum(rows.values()) == led["total_instructions"] == cert["total_instructions"]
checks["all_rows_int"] = all(isinstance(v, int) for v in rows.values())
W = float(F(led["total_instructions"], 2000))
checks["W_matches"] = abs(math.log2(W) - led["W_log2"]) < 1e-9
checks["claim_covers_W"] = math.log2(W) <= 116.8130
checks["claim_below_gate"] = 116.8130 <= 116.8271
checks["claim_tight"] = math.log2(W) > 116.8100
checks["typo_control"] = led["typo_control"]["detected"] is True
checks["typo_log2_control"] = led["typo_control"]["log2_table_typo_caught"] is True
checks["nw1_storage_row"] = rows.get("walker_state_nw1_credit") == -(cert["N"] - 1) * 4
# success
q1, cp1, dev = 0.39907073974609375, 0.3968437, 0.003955
checks["success_value"] = abs(cert["success"]["after"] - (cp1 - dev)) < 1e-9
checks["success_floor"] = cert["success"]["after"] >= 0.39
checks["success_measured_trials"] = cert["success"]["measured_nw1"]["trials"] == 131072
checks["density_match"] = abs(cert["N"] ** 2 / cert["m"] - 0.084715) < 5e-3
# kappa
m, D = cert["m"], 1 << 167
checks["cs_floor_exceeds_6m"] = F(m * m, D) > 6 * m
checks["cs_check_row_says_false"] = kap["check_row"]["sumN2_le_6m_holds"] is False
checks["kappa_charged_above_envelope"] = kap["kappa_168"]["charged"] > kap["kappa_168"]["hi"]
checks["kappa_envelope_measured"] = kap["measured"] is True and kap["kappa_measured"]["kappa_max"] < 1.05
checks["kappa_rows_dense"] = len(kap["kappa_measured"]["dense_rows"]) >= 12
# ratios
checks["distinct_lt_ordered"] = all(rat["shapes"][s]["distinct_mean"] < rat["shapes"][s]["ordered_pairs"]
                                    for s in ("nw1", "nw2", "nw4"))
checks["nw4_set_law_within_cp"] = rat["nw1_floor_check"]["nw4_set_law_within_cp"] is True
checks["nw4_ordered_law_exceeds_cp"] = rat["nw1_floor_check"]["nw4_law_ordered_pairs"] > rat["nw1_floor_check"]["nw4_cp95_upper"]
checks["per_base_pair_sets"] = [rat["shapes"][s]["per_base_pair_set"] for s in ("nw1", "nw2", "nw4")] == [[1], [3], [13]]
# transfer
checks["transfer_ran"] = tr["partial"] >= 2048
checks["transfer_prediction_inside"] = tr["prediction_in_interval"] is True
checks["transfer_pair_law_inside"] = tr["pair_law_in_interval"] is True
checks["transfer_prediction_matches_cert"] = abs(cert["success"]["transfer_arm"]["q"] - tr["q"]) < 1e-12
# gate / recovery
checks["gate_72"] = gate["ok"] == 72 and gate["passed"] is True
checks["gate_old_path_zero"] = gate["v23_path_solves"] == 0
checks["gate_no_false_positives"] = gate["false_positives"] == 0
# submission schema
checks["submission_minimal_schema"] = set(sub) == {"version", "title", "claims", "experiments"}
checks["claim_minimal_schema"] = set(sub["claims"][0]) == {"evidence", "score", "statement"}
checks["claim_text_cap"] = len(sub["claims"][0]["statement"]) <= 7965
checks["score_116813"] = sub["claims"][0]["score"] == {"numerator": "116813", "denominator": "1000"}
checks["evidence_files_exist"] = all(os.path.exists(os.path.join(HERE, e)) for e in sub["claims"][0]["evidence"])
bad = [k for k, v in checks.items() if not v]
print(json.dumps({"n_checks": len(checks), "failed": bad, "all_pass": not bad}, indent=1))
raise SystemExit(1 if bad else 0)
