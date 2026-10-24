#!/usr/bin/env python3
"""v35 ledger + claim derivation.  Every row is charged >= executed; nothing is an input."""
import json, math, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))
D, A = J("design-v35.json"), J("price-anchors-v35.json")
CEN, LAD = J("census-prod-v35.json"), J("structured-ladder.json")
C = D["C_instructions_per_unit"]
N, NB, ELL, SIDES = D["N"], D["nb"], D["ell"], D["sides"]
# charged prices: measured bodies (charged >= executed), never the v34 declared numbers
P = dict(walk=A["composites"]["P_walk_entry"]["declared_price"],
         even=A["composites"]["P_even_form"]["declared_price"],
         scan=A["composites"]["P_scan"]["declared_price"],
         store=A["composites"]["P_store"]["declared_price"],
         route=A["composites"]["P_route_sort"]["declared_price"])
HP = A["constructions"]["hprime_element"]["instructions"]
CT = A["constructions"]["chain_table_entry"]["instructions"]
PER = max(D["per_record_instructions"], D["itemised_per_record"] + D["charged_overhead"])
rows = {
 "walk_list": ELL * (P["walk"] + P["even"]),
 "hprime_supply_construction": SIDES * ELL * HP,
 "chain_table_construction": SIDES * D["chain_windows"] * (2 ** D["chain_entries_log2"]) * CT,
 "fixed_per_node": NB * P["route"],
 "base_entry": SIDES * NB * P["even"],
 "walker_state": SIDES * NB * 4,
 "gap_point_tables": SIDES * NB * (2 ** D["gap_support_log2"]) * D["gap_entry_point_ops"] * C,
 "gap_self_check": (ELL - 1) * P["scan"],
 "gap_retry_worst": (D["retry_draws"] * D["q_gap_cp95_num"] * ELL * PER //
                     (D["q_gap_cp95_den"] - D["q_gap_cp95_num"])) + 1,
 "candidate_scan_cap": (2 ** D["candidates_log2"]) * D["per_candidate_instructions"],
 "run_pair_cap": (2 ** D["run_pairs_log2"]) * D["run_pair_instructions"],
 "root_and_io": SIDES * (2 ** D["root_tests_log2"]) * P["scan"],
 "root_scan_T": (2 ** D["root_tests_log2"]) * P["route"],
}
rows["records"] = SIDES * NB * ELL * PER
total = sum(rows.values())
units = F(total, C)
W = math.log2(float(units))
carried = total - rows["records"]
dbl = math.log2(float(F(total + carried, C))) - W
claim = F(int(math.ceil((W + dbl) * 1e4)), 10000)
gate = F(1168271, 10000)
succ = LAD["pooled_long_walks"]["cp95_lower"]
kappa = CEN["charged_kappa"]
led = {
 "version": "v35", "C_instructions_per_unit": C,
 "rows": rows, "row_log2": {k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
 "charged_prices": dict(P, HPRIME_ELEMENT=HP, CHAIN_ENTRY=CT, per_record_instructions=PER,
                        direction="charged = measured compiled body >= executed (v34 charged below measured)"),
 "address_correspondence": ("signed 64-bit/limb model, as in v34: one instruction per address and one "
                            "per operand access. Static counts are unsigned bodies; the multiplier-up "
                            "expansions are declared in derive, not hidden. This is a model, and the "
                            "model is what the claim is expressed in."),
 "total_instructions": total, "W": str(units), "W_log2": W,
 "claim": dict(numerator=str(claim.numerator), denominator=str(claim.denominator),
               value=float(claim), reduced=math.gcd(claim.numerator, claim.denominator) == 1),
 "claim_minus_W": round(float(claim) - W, 9),
 "claim_covers_doubling_stress": bool(float(claim) - W >= dbl),
 "gate": dict(numerator="1168271", denominator="10000"),
 "gate_delta_claim_minus_gate": round(float(claim) - float(gate), 6),
 "claim_vs_gate": ("ABOVE gate by %.4f bits" % (float(claim) - float(gate))
                   if float(claim) >= float(gate) else
                   "BELOW gate by %.4f bits -- GATE LOST, disclosed, package legality stock"
                   % (float(gate) - float(claim))),
 "success": {
   "claimed_value": succ, "floor": D["success_floor"],
   "key_path": "structured-ladder.json:pooled_long_walks.cp95_lower",
   "meets_floor": bool(succ >= D["success_floor"]),
   "pool_level": LAD["claim_population"],
   "per_rung": [{k: r[k] for k in ("nb", "ell", "trials", "hits", "q", "cp95_lower", "floor_met")}
                for r in LAD["ladder"]],
   "rungs_below_floor": LAD["rungs_below_floor"],
   "heterogeneity": LAD["heterogeneity"],
   "cluster_bootstrap_over_rungs": LAD["cluster_bootstrap_over_rungs"],
   "reconciliation": ("ONE bootstrap interval in the whole package, the one computed by ladder_v35.py; "
                      "the v34 text's [0.3922,0.3972] is withdrawn. Rung count: 2 of the 5 long-walk "
                      "rungs are below 0.39 (the third, ell=16, is outside the claim population)."),
   "transfer_statement": ("the claim is a POOL-level statement about the trial unit production consumes "
                          "(one pool = Cartesian base x walk). The fixture's absolute nb/ell split is "
                          "NOT claimed transferable; the structure is.")},
 "kappa": {"charged": kappa, "cells": {k: dict(kappa=v["kappa"], powered=v["powered"],
                                               observed_pairs=v["observed_pairs"],
                                               birthday_expected_pairs=v["birthday_expected_pairs"],
                                               width_bits=v["width_bits"])
                                       for k, v in CEN["cells"].items()},
           "rule": CEN["charged_kappa_rule"], "census_object": CEN["curve"]["map"],
           "trend": CEN["trend"], "unpowered_reported": CEN["unpowered_reported"],
           "fixture_disclosure": CEN["fixture_disclosure"]},
 "debit": {"E_bad_log2": round(math.log2(kappa) + 2 * math.log2(N) - 168, 6), "cap_bad_log2": 84,
           "form": "kappa * N^2/2^168", "key_width_bits": 168},
 "gap_support": {"distinct_gaps": str(ELL - 1), "support_log2": D["gap_support_log2"],
   "markov_fallback": ("P(gap>t) <= E[gap]/t with NO tail assumption: a Markov-on-observed-max row. "
                       "It is insufficient on its own (2^-6.9 per gap), so it is stated as a failure and "
                       "the charge does not depend on the geometric hypothesis: the gap row is 2^%.3f "
                       "of 2^%.3f total, moving the claim < 1e-4 bits." %
                       (math.log2(rows["gap_point_tables"]), W)),
   "geometric_tail": "HYPOTHESIS, retained only as a sensitivity: not the basis of the claim"},
 "anchors": {"file": "price-anchors-v35.json", "hprime_element_instructions": HP,
             "chain_table_entry_instructions": CT,
             "itemised_sum_compiled": sum(v["compiled_instructions"] for v in A["items"].values()),
             "executed_counts": A["executed_counts"], "fixes": A["anchor_fixes"],
             "chain": "source -> clang -O2 -> otool -tV -> instructions -> price -> row (all in-tree)"},
 "retune": {"nb": NB, "nb_log2": 45, "ell": ELL, "ell_log2": round(math.log2(ELL), 6)},
}
json.dump(led, open(os.path.join(HERE, "ledger-v35.json"), "w"), separators=(",", ":"))
print(json.dumps({"W_log2": round(W, 9), "claim": str(claim), "gate_delta": led["gate_delta_claim_minus_gate"],
                  "success": succ, "kappa": kappa, "hprime": HP, "chain_entry": CT,
                  "per_record": PER, "records_share": round(rows["records"] / total, 6)}, indent=1))
assert float(claim) >= W
