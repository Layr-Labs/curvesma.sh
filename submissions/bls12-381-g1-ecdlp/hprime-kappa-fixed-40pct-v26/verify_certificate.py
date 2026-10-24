#!/usr/bin/env python3
"""v25 package verifier: re-derive every headline row from the raw evidence files.

v25 checks (on top of the v24 ledger/retry/population rows):
  * CP DIRECTION: every census cell has cp95_lower <= q <= cp95_upper, the audit's
    rows use the correct side for their declared direction, planted flip caught;
  * PAIR CORRECTION: K_distinct == K * (1 - (nw-1)/nw^2) in every cell, the fixed-W
    series has nw=1,2,4, and the design row is the nw=1 CP lower endpoint;
  * SUCCESS: reported success == min(nw=1 CP lower, production corrected-law lower)
    and meets_floor == (success >= 0.39);
  * KAPPA: slope-consistency row present with the free slope ~0.5, the pinned fit
    evaluated at 168 reported AND rejected, charged kappa above the band top, and
    the structural band matching 6*2^167/m.
"""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, "certificate.json")))
pop = json.load(open(os.path.join(HERE, "population-v22.json")))
ret = json.load(open(os.path.join(HERE, "retry-v25.json")))
direct = json.load(open(os.path.join(HERE, "direct-v25.json")))
gate = json.load(open(os.path.join(HERE, "gate-v25.json")))
kap = json.load(open(os.path.join(HERE, "kappa-v25.json")))
s = direct["summary"]
ok = []

enum_ok = all(x["sum_Nk"] == f["m"] and f["h_enumerated"] == f["m"]
              for f in pop["fixtures"] for x in f["rows"])
ok.append(("population_sum_Nk_equals_m", enum_ok))
t = pop["summary"]["trend"]
ok.append(("weights_trend_monotone_to_1",
           all(t[i]["weight_ratio_mean"] >= t[i + 1]["weight_ratio_mean"] - 0.02
               for i in range(len(t) - 1)) and t[-1]["weight_ratio_mean"] > 0.999))
ok.append(("coverage_debit_deleted", C["coverage"]["debit"] == 0
           and C["success"]["debits"]["coverage"] == 0))

regs = ret["retry"]["per_regime"]
worst = max(regs, key=lambda c: c["q_gap_cp95_upper"])
qg = F(str(worst["q_gap_cp95_upper"]))
cost = C["retry"]["build_draw_cost"]
ok.append(("retry_two_regimes", len(regs) == 2))
ok.append(("retry_draws_per_regime_ge_192", all(c["draws"] >= 192 for c in regs)))
ok.append(("retry_gap_row_matches_worst_regime",
           abs(C["rows"]["gap_retry_worst"] / (cost * (qg / (1 - qg))) - 1) < 1e-9))
ok.append(("retry_worst_is_max", abs(C["retry"]["worst_regime_used"]["q_gap_cp95_upper"]
                                     - worst["q_gap_cp95_upper"]) < 1e-12))
ok.append(("route_row_deleted", "route_retry_worst" not in C["rows"]))

ok.append(("ledger_sums", sum(C["rows"].values()) == C["total_instructions"]))
W = F(C["total_instructions"], 2000)
claim = F(int(C["claim"]["numerator"]), int(C["claim"]["denominator"]))
ok.append(("claim_above_W", claim > math.log2(float(W))))
ok.append(("claim_below_gate", claim < F(1168271, 10000)))
ok.append(("ledger_rows_11", len(C["rows"]) == 11))

ok.append(("gate_known_answer_passed", bool(gate["passed"] and gate["ok"] == gate["n"])))
ok.append(("gate_old_path_zero", gate["v23_path_solves"] == 0))
ok.append(("gate_no_false_positives", gate["false_positives"] == 0
           and gate["false_positive_n"] >= 64))

aud = s["cp_direction_audit"]
ok.append(("cp_audit_passed", bool(aud["passed"])))
ok.append(("cp_audit_planted_flip_caught", bool(aud["planted_flip_all_caught"])))
ok.append(("cp_audit_rows_direction_ok", all(r["ok"] and r["side_ok"] for r in aud["rows"])))
ok.append(("cp_brackets_every_cell",
           all(g["cp95_lower"] <= g["q"] <= g["cp95_upper"] for g in direct["groups"])))

ok.append(("pair_correction_identity",
           all(abs(g["K_distinct"] - g["K"] * g["pair_correction"]) < 1e-9
               for g in direct["groups"])))
ok.append(("pair_correction_values",
           all(abs(g["pair_correction"] - (1 - (g["nw"] - 1) / g["nw"] ** 2)) < 1e-12
               for g in direct["groups"])))
ok.append(("fixed_W_series_nw_1_2_4",
           sorted(v["nw"] for v in s["by_nw"].values()) == [1, 2, 4]))
ok.append(("fixed_W_same_records",
           len({g["records_per_side"] for g in direct["groups"]}) == 1))
ok.append(("design_row_is_nw1_cp_lower",
           abs(s["design_row"]["cp95_lower"] - s["by_nw"]["nw1"]["cp95_lower"]) < 1e-12))
ok.append(("success_row_is_min_of_law_and_nw1",
           abs(s["success_final"] - min(s["by_nw"]["nw1"]["cp95_lower"],
                                        s["production"]["p_lower"])) < 1e-7))
ok.append(("success_flag_matches_floor",
           bool(C["success"]["meets_floor"] == (C["success"]["after"] >= 0.39))))
ok.append(("success_meets_floor_claimed", bool(s["meets_floor"])))
ok.append(("success_after_matches_summary",
           abs(C["success"]["after"] - s["success_final"]) < 1e-12))
ok.append(("trials_ge_4k_per_cell", all(g["trials"] >= 4096 for g in direct["groups"])))
ok.append(("all_points_on_curve",
           all(g["points_on_curve"] == g["records_per_trial"] * g["trials"]
               for g in direct["groups"])))
ok.append(("budget_K_not_above_claim", all(g["K"] <= 0.51 for g in direct["groups"])))

ok.append(("kappa_slope_row_present",
           abs(kap["slope_consistency"]["big"]["measured_slope_free_fit"] - 0.5) < 0.1))
ok.append(("kappa_formula_deviation_row",
           abs(kap["slope_consistency"]["big"]["deviation_from_formula"] - 0.5) < 0.1))
ok.append(("kappa_pinned_fit_at_168_reported_and_rejected",
           kap["kappa_168_pinned_fit_log2"] > 100
           and kap["pinned_fit_rejected"]["admissibility"].startswith("not admitted")))
ok.append(("kappa_structural_band_matches_law",
           abs(kap["structural_excess_168_log2"]
               - (math.log2(6) + 167 - math.log2(C["m"]))) < 1e-6))
ok.append(("kappa_charged_above_band", bool(kap["charged_exceeds_band"])))
ok.append(("kappa_measured_false", kap["measured"] is False))
ok.append(("all_certificate_checks", all(C["checks"].values())))
ok.append(("d_times_m_is_r_minus_1",
           C["d"] * C["m"] == 52435875175126190479447740508185965837690552500527637822603658699938581184512))
ok.append(("m_v2_12", C["m"] % 4096 == 0 and C["m"] % 8192 != 0))

bad = [k for k, v in ok if not v]
for k, v in ok:
    print("  %-44s %s" % (k, "PASS" if v else "FAIL"))
print("verified %d checks" % len(ok))
json.dump(dict(checks={k: v for k, v in ok}, PASS=not bad, count=len(ok), failed=bad,
               W_log2=round(C["raw_log2"], 6), claim=float(claim),
               success=C["success"]["after"], meets_floor=C["success"]["meets_floor"]),
          open(HERE + "/verify-v25.json", "w"), indent=1)
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL PASS  W=2^%.6f  claim=%.4f  success=%.6f (floor_met=%s)"
      % (C["raw_log2"], float(claim), C["success"]["after"], C["success"]["meets_floor"]))
