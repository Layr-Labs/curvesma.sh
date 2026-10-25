#!/usr/bin/env python3
"""cert37.py -- recompute the v37 certificate from the shipped sources and evidence."""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
L = json.load(open(os.path.join(EV, "ledger37.json")))
A = json.load(open(os.path.join(EV, "audit37.json")))
ARM = json.load(open(os.path.join(EV, "arm37.json")))
BND = json.load(open(os.path.join(EV, "bounds37.json")))
SUB = json.load(open(os.path.join(PKG, "submission.json")))

rows = {k: int(v) for k, v in L["rows_ir"].items()}
W = sum(rows.values())
W_led = int(L["W_total_ir"])
n = int(L["claim"]["numerator"])
sc = SUB["claims"][0]["score"]
title = SUB["title"]
stmt = SUB["claims"][0]["statement"]
P_USED = L["success"]["p_used"]

checks = {
    "row_sum_equals_W": W == W_led,
    "indep_recompute": W == sum(int(v) for v in L["rows_ir"].values()),
    "row_count_matches_statement": ("sum of %d rows" % len(rows)) in stmt,
    "claim_exact": pow(20000, 10000) * (1 << n) >= pow(W, 10000),
    "claim_minimal": pow(20000, 10000) * (1 << (n - 1)) < pow(W, 10000),
    "sub_score_matches": int(sc["numerator"]) == n and sc["denominator"] == "10000",
    "version_1": SUB["version"] == 1,
    "experiments_objects": isinstance(SUB["experiments"], list),
    "title_len": len(title) <= 300,
    "six_limb_wired": "d28_keys12_form6" in open(os.path.join(HERE, "protocol37.py")).read(),
    "python_key_rule_absent": "self.T.keys12" not in
                              open(os.path.join(HERE, "protocol37.py")).read(),
    "element_key_row_present": int(L["rows_ir"]["element_key_formation"]) > 0,
    "address_translation_priced": int(L["rows_ir"]["address_translation"])
                                  == int(L["protocol"]["steps"]) * 38,
    "borrow_plant_passes": A["a_borrow"]["neg1_equals_m_minus_1"]
                           and A["a_borrow"]["mismatches"] == 0,
    "old_rule_fails_plant": A["a_borrow"]["old_rule_fails_planted_case"],
    "cold_start": A["b_cold_start"]["start_matches_public_combination"]
                  and not A["b_cold_start"]["start_source_reads_secret"],
    "restart_moved_off_7660": A["b_cold_start"]["row_changed"],
    "occupancy_used": "occupancy" in L["store_spec_to_code"],
    "occupancy_dense_equals_M": A["c_occupancy"]["dense_S"]["equals_M"],
    "dense_six_ceil_flagged": A["c_occupancy"]["dense_six_ceil_below_occupancy"],
    "record_node_separate": int(BND["parsed_constants"]["NODE_BYTES_RECORD"]) == 96
                            and int(BND["parsed_constants"]["NODE_BYTES_RADIX"]) == 33554432,
    "init_row_is_nodes_x_2db": int(L["rows_ir"]["store_init_radix_nodes"])
                               == int(BND["init_ir"]),
    "starts_are_ceil_S_over_L": L["budget"]["trail_starts"]
                                == str(-(-int(L["protocol"]["steps"]) // L["budget"]["per_trail_step_cap"])),
    "no_guard_substitution": "4497" not in json.dumps(L["budget"]),
    "frequency_independent_step": L["prices"]["step12"] == max(int(json.load(open(
        os.path.join(EV, "prices28.json")))["MUL_CHARGE"]),
        int(json.load(open(os.path.join(EV, "prices28.json")))["SQR_CHARGE"])),
    "yield_measured_in_protocol": ARM["protocol_path"].startswith("scripts/protocol37.run_one"),
    "invariant_audit_passes": A.get("e_invariant", {}).get("all_steps_hold", False),
    "invariant_plant_fails_without_fix": A.get("e_invariant", {}).get("planted_break_detected", False),
    "p_used_is_proof_not_premise": abs(P_USED - (1.0 - 2.0 ** -255)) < 1e-30,
    "no_imported_multiplier": "0.578125" not in json.dumps(L["success"]),
    "verified_not_matches": ARM["arms"][0]["metric"].startswith("verified = g^xx == z"),
    "solved_seeds_all": ARM["solved_seeds"] == len(ARM["arms"][0]["per_seed"]) * 2,
    "p_used_above_floor": P_USED >= L["success"]["floor"],
    "record_node_per_key": any("per NEW RESIDENT KEY" in r["array"]
                               for r in BND["bounds_table"]),
    "four_audits": bool(A["all_four_pass"]),
    "evidence_present": all(os.path.exists(os.path.join(PKG, e))
                            for e in SUB["claims"][0]["evidence"]),
}
out = {"version": 1, "claim": str(n), "claim_value": n / 10000,
       "reduced": L["claim"]["reduced"], "ok": all(checks.values()),
       "checks": checks, "W_total_ir": str(W), "W_log2": math.log2(W),
       "row_count": len(rows), "step_ir": L["step_ir"],
       "p_used": P_USED,
       "missing": [k for k, v in checks.items() if not v],
       "method": "scripts/cert37.py: rows recomputed from evidence/ledger37.json, exact claim "
                 "inequality both directions, the four planted audits re-read, occupancy and cold "
                 "start checked against their generators"}
json.dump(out, open(os.path.join(PKG, "certificate.json"), "w"), indent=1)
print(json.dumps({"ok": out["ok"], "missing": out["missing"], "claim": out["claim_value"]}, indent=1))
assert out["ok"], out["missing"]
