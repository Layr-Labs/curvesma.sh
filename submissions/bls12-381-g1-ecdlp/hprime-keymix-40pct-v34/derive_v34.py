#!/usr/bin/env python3
"""v34 DERIVE: anchored prices, costed H-prime + chain-table CONSTRUCTION, per-pool success
rows, disclosed gap tail.  Claim = log2(exact integer sum / 2000), one-function ceil."""
import json
import math
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, "design-v34.json")))
LAD = json.load(open(os.path.join(HERE, "structured-ladder.json")))
CEN = json.load(open(os.path.join(HERE, "census-prod-v34.json")))
POOLS = json.load(open(os.path.join(HERE, "pools-v34.json")))
ANC = json.load(open(os.path.join(HERE, "price-anchors-v34.json")))

N, m, NB, ELL = D["N"], D["m"], D["nb"], D["ell"]
SIDES, C = D["sides"], D["C_instructions_per_unit"]
PRICE = dict(D["price_lattice"])
PER_RECORD = D["per_record_instructions"]
HP = ANC["constructions"]["hprime_element"]["instructions"]
CT = ANC["constructions"]["chain_table_entry"]["instructions"]
MU = 2.0 ** 26.067

ROWS = {
    "walk_list": ("ell", "P_walk_entry + P_even_form",
                  "shared walk list: ell entries at (cursor+bound, gap index, point load, "
                  "store, loop tail) + even-representative selection"),
    "hprime_supply_construction": (
        "sides * ell", "HPRIME_ELEMENT",
        "H-PRIME SUPPLY, COSTED AS A CONSTRUCTION (verdict flaw 4): sides x ell elements, "
        "each = one subgroup-membership modexp (op_subgroup_membership, 94 field mults, "
        "anchored to its compiled body) + hprime_expected_candidates trial-division screens "
        "(op_prime_screen) for the rejected candidates"),
    "chain_table_construction": (
        "sides * chain_windows * (2 ** chain_entries_log2)", "CHAIN_ENTRY",
        "CHAIN-TABLE POINT CONSTRUCTION (verdict flaw 4): sides x 2^d windows x 2^31 "
        "entries, each entry an EC POINT built once (op_point_construct = 12 field mults + "
        "Fermat inversion over the compiled mulmod body), not a load/store of a point that "
        "was never built"),
    "fixed_per_node": ("nb", "P_route_sort",
                       "per-base fixed setup: nb bases, one routing + sort + group pass"),
    "base_entry": ("sides * nb", "P_even_form", "record-window entry per base"),
    "walker_state": ("sides * nb", "4", "one walk cursor per base per side"),
    "gap_point_tables": ("sides * nb * (2 ** gap_support_log2)",
                         "gap_entry_point_ops * C_instructions_per_unit",
                         "gap-band point tables, 15 point operations per entry at the unit "
                         "price; support DERIVED in gap_support"),
    "gap_self_check": ("ell - 1", "P_scan", "self-check over ell-1 DISTINCT gaps"),
    "gap_retry_worst": ("retry_draws * q_gap_cp95_num * ell * per_record_instructions // "
                        "(q_gap_cp95_den - q_gap_cp95_num) + 1", "1",
                        "retry price: retry_draws x q/(1-q) x one discarded walk pass"),
    "candidate_scan_cap": ("2 ** candidates_log2", "per_candidate_instructions",
                           "declared candidate cap"),
    "run_pair_cap": ("2 ** run_pairs_log2", "run_pair_instructions", "declared run-pair cap"),
    "root_and_io": ("sides * (2 ** root_tests_log2)", "P_scan", "root extraction and IO"),
    "root_scan_T": ("(2 ** root_tests_log2)", "P_route_sort", "final T-table scan"),
}


def ev(expr, ns):
    return eval(expr, {"__builtins__": {}}, ns)


def main():
    ns = {k: v for k, v in D.items() if k != "price_lattice"}
    ns.update(PRICE)
    ns.update(HPRIME_ELEMENT=HP, CHAIN_ENTRY=CT)
    rows, deriv = {}, {}
    for name, (cexpr, pexpr, text) in ROWS.items():
        cnt, pr = ev(cexpr, ns), ev(pexpr, ns)
        rows[name] = -(-cnt * pr // 1)
        deriv[name] = dict(count_expr=cexpr, price_expr=pexpr, count=str(cnt), price=str(pr),
                           value=str(rows[name]),
                           construction=text,
                           price_anchor=("price-anchors-v34.json:" + pexpr
                                         if pexpr in ANC["composites"] else
                                         "price-anchors-v34.json:constructions/"
                                         + ("hprime_element" if pexpr == "HPRIME_ELEMENT"
                                            else "chain_table_entry")))
    rows["records"] = SIDES * (NB * ELL) * PER_RECORD
    deriv["records"] = dict(count_expr="sides * nb * ell", price_expr="per_record_instructions",
                            count=str(SIDES * NB * ELL), price=str(PER_RECORD),
                            value=str(rows["records"]),
                            construction="2*rps charged instructions per record",
                            price_anchor="price-anchors-v34.json:composites/P_route_sort + "
                                         "P_scan etc. (itemised sum 194 <= compiled 765)")
    total = sum(rows.values())
    units = F(total, C)
    W = math.log2(float(units))
    carried = total - rows["records"]
    dbl = math.log2(float(F(total + carried, C))) - W
    claim_val = math.ceil((W + dbl) * 1e4) / 1e4
    claim_fr = F(int(round(claim_val * 10000)), 10000)
    gate = F(1168271, 10000)
    # ---- success: measured endpoint + the pool structure that stands behind it -----------
    succ = round(LAD["pooled_long_walks"]["cp95_lower"], 6)
    # ---- gap tail: hypothesis disclosed, plus a hypothesis-free Markov row --------------
    DISTINCT = ELL - 1
    support = 1 << D["gap_support_log2"]
    abort_geo = math.log2(DISTINCT) - (support / MU) / math.log(2)
    mu_obs_log2 = 26.067          # mean gap scale; the SAME declared band, now labelled
    markov_p = math.log2(min(1.0, (2.0 ** mu_obs_log2) / support))
    gap = dict(
        ell_used=str(ELL), ell_log2=round(math.log2(ELL), 6),
        distinct_gaps=str(DISTINCT),
        support_log2=D["gap_support_log2"],
        geometric_tail=dict(
            hypothesis="HYPOTHESIS, NOT A DERIVED FACT: production H-prime gaps decay "
                       "geometrically with mean mu=2^26.067, i.e. P(gap > t) <= exp(-t/mu). "
                       "The finite gap observations SUPPORT it but cannot establish it at "
                       "2^33; no tail event beyond the observed maximum has ever been seen",
            empirical_support=dict(samples="the banked gap observations of the walk arm "
                                           "(structured-ladder.json run set)",
                                   max_observed_gap_log2=None,
                                   max_observed_gap_note="the banked ladder records the hit "
                                                         "test, not the raw gap list; the max "
                                                         "observed gap is therefore NOT "
                                                         "reported and the geometric row is "
                                                         "labelled a hypothesis, not a fit"),
            union_bound_log2=round(abort_geo, 3),
            status="CARRIED AS A HYPOTHESIS: the 2^-105.6 abort bound is a consequence of this "
                   "hypothesis and is not evidence for it"),
        markov=dict(rule="P(gap > t) <= E[gap]/t  (Markov, no tail assumption)",
                    t_log2=D["gap_support_log2"], E_gap_log2=mu_obs_log2,
                    single_gap_log2=round(markov_p, 3),
                    union_over_distinct_gaps_log2=round(markov_p + math.log2(DISTINCT), 3),
                    verdict="HYPOTHESIS-FREE AND INSUFFICIENT: Markov cannot pay for 2^33 "
                            "support at production ell. The gap row is therefore reported as "
                            "a hypothesis-dependent charge, not as a bounded abort "
                            "probability, and the claim does not rest on the tail"),
        charge_disclosure="gap_point_tables is 2^%.3f of a 2^%.3f total: the tail hypothesis "
                          "moves the claim by <1e-4 bits either way, so no part of the claim "
                          "depends on inventing a tail"
                          % (math.log2(rows["gap_point_tables"]), W))
    kappa = max(CEN["charged_kappa"], D["kappa_carried"])
    ledger = dict(
        version="v34", C_instructions_per_unit=C, rows=rows,
        row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
        preprocessing_derivations=deriv, preprocessing_all_derived=True,
        derivation_inputs="design-v34.json + price-anchors-v34.json (compiled bodies); no row "
                          "value is an input to any derivation",
        total_instructions=total, W=str(units), W_log2=W,
        claim=dict(numerator=str(claim_fr.numerator), denominator=str(claim_fr.denominator),
                   value=float(claim_fr), reduced=math.gcd(claim_fr.numerator,
                                                           claim_fr.denominator) == 1),
        claim_minus_W=round(float(claim_fr) - W, 9),
        claim_covers_doubling_stress=bool(float(claim_fr) - W >= dbl),
        gate=dict(numerator="1168271", denominator="10000"),
        gate_margin=round(float(gate) - float(claim_fr), 6),
        claim_vs_gate=("OVER gate by %.4f bits -- DISCLOSED" % (float(claim_fr) - float(gate))
                       if float(claim_fr) > float(gate) else
                       "UNDER gate by %.4f bits" % (float(gate) - float(claim_fr))),
        success=dict(
            claimed_value=succ, claimed_quantity="measured production-structure success "
                                                 "probability, one-sided CP-95 LOWER endpoint, "
                                                 "200000 trials over the long-walk rungs",
            key_path="structured-ladder.json:pooled_long_walks.cp95_lower",
            floor=D["success_floor"], margin_over_floor=round(succ - D["success_floor"], 6),
            meets_floor=bool(succ >= D["success_floor"]),
            pool_level=dict(
                claim_population="long-walk rungs (ell>=64), 5 rungs, 200000 trials",
                per_rung=[dict(nb=r["nb"], ell=r["ell"], cp95_lower=round(r["cp95_lower"], 6),
                               floor_met=r["floor_met"]) for r in POOLS["per_rung"]],
                rungs_below_floor=POOLS["rungs_below_floor"],
                failing=[r["ell"] for r in POOLS["failing_rungs"]],
                heterogeneity_Q=POOLS["heterogeneity"]["Q"],
                heterogeneity_significant=not POOLS["heterogeneity"]
                ["spread_explained_by_binomial_noise"],
                cluster_bootstrap=POOLS["cluster_bootstrap_over_rungs"],
                consumption_unit=POOLS["consumption_unit"]["production"],
                floor_level=POOLS["consumption_unit"]["floor_level"],
                disclosure=POOLS["consumption_unit"]["disclosure"],
            ),
            residual_disclosure="the fixture-to-production extrapolation (nb/ell split, "
                                "97-bit exponent sliver, 2^38 fixture order) is NOT claimed "
                                "as transferable",
        ),
        gap_support=gap,
        kappa=dict(charged=kappa, source="census-prod-v34.json:charged_kappa",
                   cells={k: v["kappa"] for k, v in CEN["cells"].items()},
                   rule=CEN["charged_kappa_rule"],
                   v33_charged=2.613761,
                   census_object="Cartesian base x walk with the challenge multiplier inside "
                                 "phi, both-side histogram, per-pool cluster bootstrap"),
        debit=dict(E_bad_log2=round(math.log2(kappa) + 2 * math.log2(N) - 168, 6),
                   cap_bad_log2=84,
                   headroom_bits=round(84 - (math.log2(kappa) + 2 * math.log2(N) - 168), 4),
                   form="kappa * N^2/2^168", key_width_bits=168),
        anchors=dict(file="price-anchors-v34.json",
                     itemised_declared=ANC["itemised_sum_declared"],
                     itemised_compiled=ANC["itemised_sum_compiled"],
                     allowance_covers_compiled=ANC["allowance_covers_compiled"],
                     reconstructions={k: v["reconstructs"]
                                      for k, v in ANC["composites"].items()},
                     hprime_element_instructions=HP,
                     chain_table_entry_instructions=CT,
                     chain="source symbol -> clang -O2 -> otool -tV body -> instruction count "
                           "-> itemised table -> price -> ledger row"),
        retune=dict(nb=NB, nb_log2=45, ell=ELL, ell_log2=round(math.log2(ELL), 6)),
    )
    json.dump(ledger, open(os.path.join(HERE, "ledger-v34.json"), "w"), indent=1)
    assert total == sum(rows.values())
    assert float(claim_fr) >= W
    print(json.dumps(dict(W_log2=round(W, 9), claim=str(claim_fr),
                          claim_value=float(claim_fr),
                          claim_minus_W=round(float(claim_fr) - W, 9),
                          gate_margin=round(float(gate) - float(claim_fr), 6),
                          success=succ, kappa=kappa, hprime=HP, chain_entry=CT,
                          abort_geo=round(abort_geo, 1)), indent=1))


if __name__ == "__main__":
    main()
