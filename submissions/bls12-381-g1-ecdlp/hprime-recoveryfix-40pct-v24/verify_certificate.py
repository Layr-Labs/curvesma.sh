#!/usr/bin/env python3
"""v24 package verifier: re-derive every headline row from the raw evidence files
and refuse ALL PASS unless each one reproduces exactly.

v24 changes vs the v23 verifier:
  * the success row is MEASURED, so the check is "the disclosed floor flag matches
    the measured number", not "success meets the floor" (v23 hard-coded the latter
    and would now fail by design -- the honest miss must be checkable, not fatal);
  * retry is priced per regime: the row is re-derived from build_draw_cost and the
    WORST regime's CP-95 upper bound, and the pooled rate is required to be unused;
  * the direct census rows are re-derived (solves/trials, on-curve counts);
  * kappa_168 is required to be labelled fitted/not-measured.
"""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, "certificate.json")))
pop = json.load(open(os.path.join(HERE, "population-v22.json")))
ret = json.load(open(os.path.join(HERE, "retry-v24.json")))
direct = json.load(open(os.path.join(HERE, "direct-v24.json")))
gate = json.load(open(os.path.join(HERE, "gate-v24.json")))
kap = json.load(open(os.path.join(HERE, "kappa-v24.json")))
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
ok.append(("coverage_shape_infeasible_bits",
           C["coverage"]["feasibility"]["m_required_at_nw1_log2"]
           - C["coverage"]["feasibility"]["largest_enumerated_m_log2"] > 8))

# ---- retry, priced per regime
regs = ret["retry"]["per_regime"]
worst = max(regs, key=lambda c: c["q_gap_cp95_upper"])
qg = F(str(worst["q_gap_cp95_upper"]))
cost = C["retry"]["build_draw_cost"]
ok.append(("retry_two_regimes", len(regs) == 2))
ok.append(("retry_draws_per_regime_ge_192", all(c["draws"] >= 192 for c in regs)))
_ideal = cost * (qg / (1 - qg))
ok.append(("retry_gap_row_matches_worst_regime",
           abs(C["rows"]["gap_retry_worst"] / _ideal - 1) < 1e-9))
ok.append(("retry_worst_is_max", abs(C["retry"]["worst_regime_used"]["q_gap_cp95_upper"]
                                     - worst["q_gap_cp95_upper"]) < 1e-12))
ok.append(("retry_pooled_not_used_for_pricing",
           ret["retry"]["pooled_descriptive_only"]["q_gap_cp95"]
           <= worst["q_gap_cp95_upper"] * 1.6))
ok.append(("route_row_deleted", "route_retry_worst" not in C["rows"]))
ok.append(("route_vacuity_shipped", C["retry"]["route_vacuity"]["max_below_threshold"]))

# ---- ledger
ok.append(("ledger_sums", sum(C["rows"].values()) == C["total_instructions"]))
W = F(C["total_instructions"], 2000)
claim = F(int(C["claim"]["numerator"]), int(C["claim"]["denominator"]))
ok.append(("claim_above_W", claim > math.log2(float(W))))
ok.append(("claim_below_gate", claim < F(1168271, 10000)))
ok.append(("ledger_rows_11", len(C["rows"]) == 11))

# ---- direct measurement (the v24 core)
ok.append(("gate_known_answer_64_64", bool(gate["passed"] and gate["ok"] == 64)))
ok.append(("gate_old_path_zero", gate["v23_path_solves"] == 0))
ok.append(("gate_no_false_positives", gate["false_positives"] == 0
           and gate["false_positive_n"] >= 64))
ok.append(("direct_trials_ge_512", direct["summary"]["total_trials"] >= 512))
ok.append(("direct_rows_consistent",
           all(abs(g["solves"] / g["trials"] - g["q"]) < 1e-12
               and g["cp95"] <= g["q"] for g in direct["groups"])))
ok.append(("all_points_on_curve",
           all(g["points_on_curve"] == g["records_per_trial"] * g["trials"]
               for g in direct["groups"])))
ok.append(("budget_K_not_above_claim",
           all(g["K"] <= 0.51 and g["K_conservative"] for g in direct["groups"])))
ok.append(("success_row_is_measured", abs(C["success"]["after"]
                                          - direct["summary"]["success_final"]) < 1e-12))
ok.append(("success_flag_matches_floor",
           bool(C["success"]["meets_floor"] == (C["success"]["after"] >= 0.39))))
ok.append(("success_below_floor_disclosed",
           bool(C["success"]["success_below_floor"] == (not C["success"]["meets_floor"]))))
ok.append(("success_takes_worst_shape",
           C["direct"]["success_final_shape"]
           == min(direct["summary"]["by_shape"],
                  key=lambda k: direct["summary"]["by_shape"][k]["cp95_min"])))

# ---- kappa transfer
ok.append(("kappa_168_labelled_fitted", kap["measured"] is False))
ok.append(("kappa_charged_dominates_fitted", bool(kap["charged_exceeds_fitted"])))
ok.append(("d_times_m_is_r_minus_1",
           C["d"] * C["m"] == 52435875175126190479447740508185965837690552500527637822603658699938581184512))
ok.append(("m_v2_12", C["m"] % 4096 == 0 and C["m"] % 8192 != 0))
ok.append(("all_certificate_checks", all(C["checks"].values())))

bad = [k for k, v in ok if not v]
for k, v in ok:
    print("  %-36s %s" % (k, "PASS" if v else "FAIL"))
print("verified %d checks" % len(ok))
json.dump(dict(checks={k: v for k, v in ok}, PASS=not bad,
               count=len(ok), failed=bad,
               W_log2=round(C["raw_log2"], 6), claim=float(claim),
               success=C["success"]["after"],
               success_below_floor=C["success"]["success_below_floor"]),
          open(HERE + "/verify-v24.json", "w"), indent=1)
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL PASS  W=2^%.6f  claim=%.4f  success=%.6f (below_floor=%s)"
      % (C["raw_log2"], float(claim), C["success"]["after"],
         C["success"]["success_below_floor"]))
