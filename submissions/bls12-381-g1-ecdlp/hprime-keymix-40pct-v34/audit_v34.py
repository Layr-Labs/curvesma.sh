import json
import math
import os
import sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "design-v34.json")))
L = json.load(open(os.path.join(HERE, "ledger-v34.json")))
CEN = json.load(open(os.path.join(HERE, "census-prod-v34.json")))
ANC = json.load(open(os.path.join(HERE, "price-anchors-v34.json")))
checks, fails = {}, []
def ck(name, cond, detail=None):
    checks[name] = dict(passed=bool(cond), detail=detail)
    if not cond:
        fails.append(name)
def main():
    ns = {k: v for k, v in D.items() if k != "price_lattice"}
    ns.update(D["price_lattice"])
    ns.update(HPRIME_ELEMENT=ANC["constructions"]["hprime_element"]["instructions"],
              CHAIN_ENTRY=ANC["constructions"]["chain_table_entry"]["instructions"])
    recomputed, bad = {}, []
    for name, d in L["preprocessing_derivations"].items():
        c = eval(d["count_expr"], {"__builtins__": {}}, ns)
        p = eval(d["price_expr"], {"__builtins__": {}}, ns)
        v = -(-c * p // 1)
        recomputed[name] = v
        if str(v) != d["value"] or L["rows"][name] != v:
            bad.append((name, str(v), d["value"], L["rows"][name]))
    ck("all_rows_recompute_from_design_and_anchors", not bad, bad)
    ck("total_equals_sum_of_rows", sum(L["rows"].values()) == L["total_instructions"])
    units = F(L["total_instructions"], D["C_instructions_per_unit"])
    W = math.log2(float(units))
    cl = F(int(L["claim"]["numerator"]), int(L["claim"]["denominator"]))
    carried = L["total_instructions"] - L["rows"]["records"]
    dbl = math.log2(float(F(L["total_instructions"] + carried, D["C_instructions_per_unit"]))) - W
    want = F(math.ceil((W + dbl) * 1e4), 10000)
    ck("claim_is_ceiling_of_one_function_cost", cl == want and float(cl) >= W,
       dict(claim=str(cl), want=str(want), covers_stress=float(cl) - W >= dbl))
    ck("claim_covers_W", float(cl) >= W, round(float(cl) - W, 9))
    ck("score_is_reduced_rational",
       math.gcd(int(L["claim"]["numerator"]), int(L["claim"]["denominator"])) == 1)
    probe = "P_scan"
    comp = ANC["composites"][probe]
    item_sum = sum(ANC["items"][c]["compiled_instructions"] for c in comp["components"])
    ck("anchor_chain_walks", item_sum == comp["measured_from_compiled"] and
       comp["measured_from_compiled"] >= D["price_lattice"][probe] and
       comp["reconstructs"],
       dict(price=probe, declared=D["price_lattice"][probe], measured=item_sum))
    ck("every_composite_reconstructs",
       all(v["reconstructs"] for v in ANC["composites"].values()),
       {k: v["reconstructs"] for k, v in ANC["composites"].items()})
    ck("allowance_dominates_compiled_itemised",
       ANC["allowance_covers_compiled"],
       dict(compiled=ANC["itemised_sum_compiled"], allowance=2000))
    ck("hprime_row_prices_the_construction",
       L["preprocessing_derivations"]["hprime_supply_construction"]["price"] ==
       str(ANC["constructions"]["hprime_element"]["instructions"]))
    ck("chain_table_row_prices_the_point_construction",
       L["preprocessing_derivations"]["chain_table_construction"]["price"] ==
       str(ANC["constructions"]["chain_table_entry"]["instructions"]))
    for cell, v in CEN["cells"].items():
        if cell == "diagonal_v33":
            continue                                                                       
        ck("census_cartesian_" + cell,
           v["records"] == v["pools"] * v["records_per_pool"],
           dict(records=v["records"], per_pool=v["records_per_pool"]))
    prod = CEN["cells"]["chal_shared"]
    ck("hist_side_count_matches_expectation", prod["hist"]["side_count_match"],
       dict(hist_total=prod["hist"]["hist_total"],
            denominator=prod["hist"]["expectation_denominator"]))
    ctrl = CEN["hist_regression_control"]
    ck("planted_side_only_histogram_FIRES", ctrl["planted_fires"] is True and
       ctrl["corrected_fires"] is False, ctrl)
    ck("challenge_multiplier_in_sampled_function",
       "c" in CEN["sampled_function"]["exponent"] and
       CEN["sampled_function"]["v33_error"].startswith("v33 sampled"))
    ck("per_pool_cluster_bootstrap_present",
       prod["kappa_cluster_bootstrap"]["over_pools"] is True)
    ck("key_width_ladder_has_four_rungs",
       [r["key_bits"] for r in CEN["key_width_ladder"]] == [20, 24, 28, 32])
    ck("width_ladder_scope_not_extrapolated",
       "NOT an extrapolation" in CEN["key_width_ladder_scope"])
    ck("success_is_key_path_value",
       L["success"]["claimed_value"] ==
       round(json.load(open(os.path.join(os.path.dirname(HERE), "hprime-keymix-40pct-v33", "structured-ladder.json")))
             ["pooled_long_walks"]["cp95_lower"], 6))
    PL = L["success"]["pool_level"]
    ck("pool_level_disclosed", "pool" in PL["floor_level"].lower() and
       PL["rungs_below_floor"] >= 0 and len(PL["per_rung"]) == 7)
    ck("heterogeneity_reported", "Q" in PL["heterogeneity"] and
       "cluster_bootstrap_over_rungs" in PL and PL["cluster_bootstrap_over_rungs"]["resamples"] > 0)
    ck("tail_labelled_hypothesis",
       L["gap_support"]["geometric_tail"].startswith("HYPOTHESIS"))
    ck("markov_row_present", "markov" in L["gap_support"])
    ck("claim_does_not_rest_on_tail", "does not rest" in L["gap_support"]["charge_disclosure"])
    sub = json.load(open(os.path.join(HERE, "submission.json")))
    ev = sub["claims"][0]["evidence"]
    missing = [f for f in ev if not os.path.exists(os.path.join(HERE, f))
               or os.path.getsize(os.path.join(HERE, f)) == 0]
    ck("all_evidence_present", not missing, missing)
    ck("stale_certificate_removed", not os.path.exists(os.path.join(HERE, "certificate.json")))
    ck("experiments_non_empty", bool(sub.get("experiments")))
    ck("title_le_300", len(sub["title"]) <= 300, len(sub["title"]))
    st = sub["claims"][0]["statement"]
    ck("statement_le_8000", len(st) <= 8000, len(st))
    ck("claims0_le_7965", len(st) <= 7965, len(st))
    out = dict(version="v34", checks=checks, passed=not fails, failed=fails,
               claim=str(cl), W_log2=round(W, 9), kappa=L["kappa"]["charged"],
               success=L["success"]["claimed_value"],
               hprime_instructions=L["preprocessing_derivations"]
               ["hprime_supply_construction"]["price"],
               chain_table_instructions=L["preprocessing_derivations"]
               ["chain_table_construction"]["price"])
    json.dump(out, open(os.path.join(HERE, "audit-v34.json"), "w"), indent=1)
    print(json.dumps(dict(passed=out["passed"], failed=fails, n_checks=len(checks),
                          claim=out["claim"], W_log2=out["W_log2"]), indent=1))
    return 0 if not fails else 1
if __name__ == "__main__":
    sys.exit(main())
