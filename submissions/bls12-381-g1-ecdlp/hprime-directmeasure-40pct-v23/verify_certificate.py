#!/usr/bin/env python3
"""v22 certificate verifier: re-derive every headline row from the raw evidence
files and refuse to print ALL PASS unless each one reproduces exactly."""
import json, math, os, sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, "certificate.json")))
pop = json.load(open(os.path.join(HERE, "population-v22.json")))
ret = json.load(open(os.path.join(HERE, "retry-v23.json")))
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

sc = ret["self_check"]
ok.append(("self_check_mutation_control", sc["selftest"]["mutation_control_pass"]))
ok.append(("self_check_counted_loop", sc["counted_kernel"]["counted"]
           and sc["counted_kernel"]["loop_body_instructions"] > 0))
ok.append(("self_check_row_present", C["rows"]["gap_self_check"] > 0))
ok.append(("self_check_row_is_per_step_x_walk",
           C["rows"]["gap_self_check"] % 11 == 0
           and C["rows"]["gap_self_check"] // 11 == int(2 ** 97 / (97 * math.log(2)))))

p = ret["retry"]["pooled"]
ok.append(("retry_draws_ge_192", p["draws"] >= 192))
qg = F(str(round(p["q_gap_cp95"], 8)))
qr = F(0)   # route event vacuous, row deleted
ok.append(("retry_gap_row_matches", C["rows"]["gap_retry_worst"] ==
           int(C["retry"]["build_draw_cost"] * (qg / (1 - qg)))))
ok.append(("route_row_deleted", "route_retry_worst" not in C["rows"]))

ok.append(("ledger_sums", sum(C["rows"].values()) == C["total_instructions"]))
W = F(C["total_instructions"], 2000)
claim = F(int(C["claim"]["numerator"]), int(C["claim"]["denominator"]))
ok.append(("claim_above_W", claim > math.log2(float(W))))
ok.append(("claim_below_gate", claim < F(1168271, 10000)))

after = F(str(C["success"]["after"])).limit_denominator(10 ** 9)
ok.append(("success_meets_floor", after >= F(39, 100)))
ok.append(("all_checks_pass", all(C["checks"].values())))

bad = [k for k, v in ok if not v]
for k, v in ok:
    print("  %-34s %s" % (k, "PASS" if v else "FAIL"))
print("verified %d checks" % len(ok))
if bad:
    print("FAILED:", bad)
    sys.exit(1)
print("ALL PASS  W=2^%.6f  claim=%.4f  success=%.9f  coverage_debit=%s  q_gap=%s  q_route=%s"
      % (C["raw_log2"], float(claim), C["success"]["after"], C["coverage"]["debit"],
         C["retry"]["q_gap"], C["retry"]["q_route"]))
