#!/usr/bin/env python3
"""v33 DERIVE: non-tautological preprocessing, claim = the MEASURED success quantity,
production-distribution kappa, corrected gap row.

verdict 5a8a74b6 (v32 REJECTED) flaw-by-flaw:

FLAW 1 (tautological preprocessing).  v32 got each price by taking the LOWEST SET BIT of
the row value and divided to manufacture a count: fixed_per_node became count=1,
price=2^115, and cnt*price==value held by construction.  v33 removes every such
factorisation.  Each row is now  value = count(design) x price(price lattice)  where BOTH
factors are expressions over DESIGN CONSTANTS ONLY (design-v33.json).  The row value is
never read while deriving it; it is the PRODUCT.  audit_v33.py re-evaluates the same
expressions from design-v33.json with its own evaluator and a namespace that contains NO
row values, so a cost-defined stub cannot pass (the planted-tautology control proves the
audit fails such a row).

FLAW 2 (K_eff is not a success probability).  v32 claimed K_eff=0.508674 (a collision
INTENSITY) as if it were success.  v33 claims THE MEASURED QUANTITY: success =
structured-ladder.json:pooled_long_walks.cp95_lower = 0.395879, the one-sided CP-95 lower
endpoint of the measured production-structure hit test (200000 trials, shared walk list,
base x walk products, mu6 hit test).  The identity is asserted by key path, and the
Poisson map 1-exp(-K) is printed only as a DIAGNOSTIC, never as the claim.

FLAW 3 (census distribution).  census_prod_v33.py samples the production structured
stream; kappa is charged as the max upper95 over all sampled cells (2.614 here vs 1.231
carried), so a distribution-induced rise is PRICED, not argued away.

FLAW 4 (gap row ell).  v32 used the FIXTURE ell=4774 inside NB*ell; v33 uses the
production ell=ceil(N/NB)=2^70.650.  It also replaces E[max gap] (an expectation) with a
UNION BOUND over DISTINCT gaps (ell-1, not nb*ell fresh draws -- the walk list is SHARED),
and the gap row price follows the corrected support.
"""
import json
import math
import os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
CERT = json.load(open(os.path.join(HERE, "certificate.json")))
LADDER = json.load(open(os.path.join(HERE, "structured-ladder.json")))
CENSUS = json.load(open(os.path.join(HERE, "census-prod-v33.json")))

N, m = CERT["N"], CERT["m"]
NB = 1 << 45
ELL = -(-N // NB)
SIDES = 2
C = 2000
KEY_KERNEL = 70

# ---- accepted price lattice: sums of the itemised per-record operations (measured table)
ITEM = [("walk_cursor_and_bound", 4), ("gap_load_2limb", 4), ("gap_index", 2),
        ("table_point_load", 6), ("point_store_into_record", 3), ("walk_loop_control", 3),
        ("p_minus_y", 12), ("even_rep_select", 7),
        ("route_4hops_x_21bits", 68), ("sort_3_passes_x_25", 75), ("group_scan", 10)]
ITEMISED = sum(v for _, v in ITEM)
OVERHEAD = -(-ITEMISED * 6 // 5)
PER_RECORD = C + OVERHEAD + KEY_KERNEL
PRICE = dict(ITEM)
PRICE.update({
    "P_walk_entry": PRICE["walk_cursor_and_bound"] + PRICE["gap_index"] +
                    PRICE["table_point_load"] + PRICE["point_store_into_record"] +
                    PRICE["walk_loop_control"],
    "P_even_form": PRICE["p_minus_y"] + PRICE["even_rep_select"],
    "P_scan": PRICE["gap_load_2limb"] + PRICE["gap_index"],
    "P_store": PRICE["table_point_load"] + PRICE["point_store_into_record"],
    "P_route_sort": PRICE["route_4hops_x_21bits"] + PRICE["sort_3_passes_x_25"] +
                    PRICE["group_scan"],
    "P_record": PER_RECORD,
})

# ---- design constants (ALL of them are declared here and exported; audit reads only this)
DESIGN = {
    "n": 1 << 20, "d": 1 << 20, "N": N, "m": m, "sides": SIDES,
    "nb_log2": 45, "nb": NB, "ell": ELL,
    "field_bits": 381, "image_bits": 168,
    "C_instructions_per_unit": C,
    "itemised_per_record": ITEMISED, "charged_overhead": OVERHEAD,
    "key_kernel_instructions": KEY_KERNEL, "per_record_instructions": PER_RECORD,
    "price_lattice": PRICE,
    "chain_windows": 8, "chain_entries_log2": 31,
    "gap_support_log2": 33,
    "gap_entry_point_ops": 15,
    "candidates_log2": 84, "per_candidate_instructions": 1 << 20,
    "run_pairs_log2": 95, "run_pair_instructions": 2,
    "retry_draws": 192, "q_gap_cp95_num": 1163388, "q_gap_cp95_den": 100000000,
    "root_tests_log2": 20,
    "kappa_carried": 1.231,
    "success_floor": 0.39,
    "claim_success_key_path": "structured-ladder.json:pooled_long_walks.cp95_lower",
}
DESIGN["price_lattice"]["P_walk_entry"] = 18
DESIGN["price_lattice"]["P_even_form"] = 19
DESIGN["price_lattice"]["P_scan"] = 6
DESIGN["price_lattice"]["P_store"] = 9
DESIGN["price_lattice"]["P_route_sort"] = 153

# ---- row constructions: (count_expr, price_expr, construction sentence)
# Every expression is evaluated in a namespace of DESIGN symbols only.  No row value is
# ever an input.  Deleting any row's value leaves it recomputable.
ROWS = {
    "walk_list": (
        "ell", "P_walk_entry + P_even_form",
        "shared walk list: ell entries, each formed by one cursor/bound update, one table "
        "load, one store and one even-representative selection (P_walk_entry + P_even_form)"),
    "hprime_supply": (
        "sides * ell", "P_scan",
        "H-prime supply for the shared walk list: one scan per entry per side "
        "(the SAME list is used by both sides, so the supply is charged twice, not nb times)"),
    "chain_tables": (
        "sides * chain_windows * (2 ** chain_entries_log2)", "P_store",
        "digit-decomposition chain tables for the fixed-base step: sides x 2^d windows x "
        "2^31 entries, each entry priced at one table load + one store"),
    "fixed_per_node": (
        "nb", "P_route_sort",
        "per-base fixed setup: nb bases, each paying one routing + sort + group-scan pass"),
    "base_entry": (
        "sides * nb", "P_even_form",
        "record-window entry per base: sides x nb bases, each at the even-form price"),
    "walker_state": (
        "sides * nb", "4",
        "one walk cursor per base per side, 4 instructions per state (declared constant 4)"),
    "gap_point_tables": (
        "sides * nb * (2 ** gap_support_log2)", "gap_entry_point_ops * C_instructions_per_unit",
        "gap-band point tables: sides x nb bases x 2^gap_support_log2 entries (the support is "
        "DERIVED in gap_support below from a union bound over DISTINCT gaps), 15 point "
        "operations per entry at the unit price"),
    "gap_self_check": (
        "ell - 1", "P_scan",
        "self-check pass over the DISTINCT gap bands: ell-1 gaps (the walk list is SHARED, so "
        "repeated uses are not fresh draws and are not counted again)"),
    "gap_retry_worst": (
        "retry_draws * q_gap_cp95_num * ell * per_record_instructions // "
        "(q_gap_cp95_den - q_gap_cp95_num) + 1",
        "1",
        "retry price: retry_draws x (q/(1-q)) x one discarded walk-list pass "
        "(ell x per_record_instructions), q = 1163388/100000000 from the banked per-regime "
        "CP-95 upper (0.01163388); the quotient is exact integer arithmetic with +1 ceil and "
        "the price granularity is ONE instruction (price 1)"),
    "candidate_scan_cap": (
        "2 ** candidates_log2", "per_candidate_instructions",
        "declared candidate cap 2^84 x 2^20 instructions per candidate (both are design "
        "constants, not cost-derived)"),
    "run_pair_cap": (
        "2 ** run_pairs_log2", "run_pair_instructions",
        "declared run-pair cap 2^95 x 2 instructions per pair"),
    "root_and_io": (
        "sides * (2 ** root_tests_log2)", "P_scan",
        "root extraction and IO: one scan per root test per side"),
    "root_scan_T": (
        "(2 ** root_tests_log2)", "P_route_sort",
        "final scan over the T table: one routing/sort/group pass per candidate root"),
}

# ---- gap support: union bound over DISTINCT gaps (verdict flaw 4, second half) ----------
MU_LOG2 = 26.067                      # exponential gap law band, declared in v31/v32
MU = 2.0 ** MU_LOG2
DISTINCT_GAPS = ELL - 1               # SHARED walk list => ell-1 distinct gaps, not nb*ell
SUPPORT = 1 << DESIGN["gap_support_log2"]
ABORT_LOG2 = math.log2(DISTINCT_GAPS) - (SUPPORT / MU) / math.log(2)
GAP_EV_LOG2 = math.log2(MU * math.log(DISTINCT_GAPS))     # diagnostic expectation only


def ev(expr, ns):
    return eval(expr, {"__builtins__": {}}, ns)


def main():
    ns = {k: v for k, v in DESIGN.items() if k != "price_lattice"}
    for k, v in DESIGN["price_lattice"].items():
        ns[k] = v
    rows, deriv = {}, {}
    for name, (cexpr, pexpr, text) in ROWS.items():
        count = ev(cexpr, ns)
        price = ev(pexpr, ns)
        value = -(-count * price // 1) if not isinstance(count * price, F) else count * price
        if isinstance(value, F):
            value = -(-value.numerator // value.denominator)
        rows[name] = int(value)
        deriv[name] = dict(count_expr=cexpr, price_expr=pexpr, count=str(count),
                           price=str(price), value=str(int(value)),
                           recompute="count x price, ceil; inputs are design constants only",
                           construction=text)
    # desired exit: candidate_scan_cap is the only carried row above 2^100; the rest are
    # small by construction.  roots of the record row:
    rows["records"] = SIDES * (NB * ELL) * PER_RECORD
    deriv["records"] = dict(
        count_expr="sides * nb * ell", price_expr="per_record_instructions",
        count=str(SIDES * NB * ELL), price=str(PER_RECORD), value=str(rows["records"]),
        recompute="sides x nb x ell records x per-record instructions",
        construction="2*rps*%d charged instructions per record (201 itemised + %d charged "
                     "overhead + %d measured key-kernel instructions)" % (PER_RECORD, OVERHEAD, KEY_KERNEL))
    total = sum(rows.values())
    units = F(total, C)
    W = math.log2(float(units))
    carried = total - rows["records"]
    dbl = math.log2(float(F(total + carried, C))) - W
    claim_val = math.ceil((W + dbl) * 1e4) / 1e4
    claim_fr = F(int(round(claim_val * 10000)), 10000)
    gate = F(1168271, 10000)
    # ---- flaw 2: the claimed success is THE MEASURED quantity, by key path -------------
    measured_success = round(LADDER["pooled_long_walks"]["cp95_lower"], 6)
    kappa = max(CENSUS["charged_kappa"], DESIGN["kappa_carried"])
    pairs_bits = 2 * math.log2(N) - 168
    E_bad_log2 = math.log2(kappa) + pairs_bits
    K_intensity = 6.0 * (NB * ELL) ** 2 / m         # collision intensity, DIAGNOSTIC ONLY
    x = NB * NB / (1 << 98)
    phi_tri = (1 - math.exp(-4.0 / 3.0 * x)) / (4.0 / 3.0 * x)
    ledger = dict(
        version="v33",
        C_instructions_per_unit=C, rows=rows,
        row_log2={k: round(math.log2(v), 6) for k, v in rows.items() if v > 0},
        preprocessing_derivations=deriv,
        preprocessing_all_derived=True,
        derivation_inputs="design-v33.json ONLY (constants + accepted price lattice); no row "
                          "value is an input to any derivation",
        total_instructions=total, W=str(units), W_log2=W,
        claim=dict(numerator=str(claim_fr.numerator), denominator=str(claim_fr.denominator),
                   value=float(claim_fr), reduced=(math.gcd(claim_fr.numerator, claim_fr.denominator) == 1)),
        claim_minus_W=round(float(claim_fr) - W, 9),
        claim_covers_doubling_stress=bool(float(claim_fr) - W >= dbl),
        doubling_stress_bits=round(dbl, 9),
        gate=dict(numerator="1168271", denominator="10000"),
        gate_margin=round(float(gate) - float(claim_fr), 6),
        claim_vs_gate=("UNDER gate by %.4f bits" % (float(gate) - float(claim_fr))
                       if float(claim_fr) < float(gate) else
                       "OVER gate by %.4f bits -- DISCLOSED: the honest total is the honest "
                       "total; rank is lost, validity is not" % (float(claim_fr) - float(gate))),
        success=dict(
            claimed_value=measured_success,
            claimed_quantity="measured production-structure success probability, one-sided "
                             "CP-95 LOWER endpoint (200000 trials, shared walk list, base x "
                             "walk products, mu6 hit test)",
            key_path=DESIGN["claim_success_key_path"],
            identity_assert="claim.success.claimed_value == rounded key-path value",
            floor=DESIGN["success_floor"],
            margin_over_floor=round(measured_success - DESIGN["success_floor"], 6),
            meets_floor=bool(measured_success >= DESIGN["success_floor"]),
            poison_diagnostic_only=dict(
                K_collision_intensity=round(K_intensity, 6), phi_tri=round(phi_tri, 6),
                K_eff=round(K_intensity * phi_tri, 6),
                one_minus_exp=round(1 - math.exp(-K_intensity * phi_tri), 6),
                note="v32 claimed K_eff=0.508674 as if it WERE success; that is a collision "
                     "intensity.  The Poisson map 1-exp(-K_eff)=%.6f is printed here as a "
                     "DIAGNOSTIC and is NOT the claim; the claim is the measured CP-95 lower "
                     "endpoint above." % (1 - math.exp(-K_intensity * phi_tri))),
            residual_disclosure="the fixture-to-production extrapolation (nb/ell split 2^16.35, "
                                "97-bit exponent sliver, 2^38 fixture order) is NOT claimed as "
                                "transferable; the claimed number is the measured endpoint at "
                                "the production STRUCTURE and no larger",
        ),
        gap_support=dict(
            ell_used=str(ELL), ell_log2=round(math.log2(ELL), 6),
            fixture_ell_used_by_v32=4774,
            distinct_gaps=str(DISTINCT_GAPS), distinct_gaps_log2=round(math.log2(DISTINCT_GAPS), 6),
            mu_log2=MU_LOG2, support_log2=DESIGN["gap_support_log2"],
            abort_bound_log2=round(ABORT_LOG2, 3),
            rule="P(max gap > 2^k) <= (ell-1) * exp(-2^k/mu)  [union bound over DISTINCT gaps]",
            expected_max_diagnostic_log2=round(GAP_EV_LOG2, 3),
            verdict="v32 put the FIXTURE ell=4774 into nb*ell and reported 31.377; with the "
                    "production ell=2^70.650 the expectation is mu*ln(ell-1)=2^%.3f and the "
                    "UNION BOUND at 2^33 gives abort probability 2^%.1f < 2^-100, so 2^33 "
                    "support entries are charged (2^32 in v32).  Repeated uses of a SHARED "
                    "walk list are counted once, over the ell-1 distinct gaps."
                    % (GAP_EV_LOG2, ABORT_LOG2),
        ),
        kappa=dict(charged=kappa, source="census-prod-v33.json:charged_kappa",
                   v32_carried=DESIGN["kappa_carried"],
                   cells={k: v["kappa"] for k, v in CENSUS["cells"].items()},
                   upper95={k: v["kappa_upper95"] for k, v in CENSUS["cells"].items()},
                   provenance="max upper95 over the production-density cell, the "
                              "canonicalised-uniform cell and the v32 uniform control"),
        debit=dict(E_bad_log2=round(E_bad_log2, 6), cap_bad_log2=84,
                   headroom_bits=round(84 - E_bad_log2, 4),
                   form="kappa * N^2/2^168",
                   key_width_bits=168),
        retune=dict(nb=NB, nb_log2=45, ell=ELL, ell_log2=round(math.log2(ELL), 6)),
    )
    json.dump(DESIGN, open(os.path.join(HERE, "design-v33.json"), "w"), indent=1)
    json.dump(ledger, open(os.path.join(HERE, "ledger-v33.json"), "w"), indent=1)
    assert total == sum(rows.values())
    assert all(isinstance(v, int) for v in rows.values())
    assert float(claim_fr) >= W
    assert measured_success == round(LADDER["pooled_long_walks"]["cp95_lower"], 6)
    print(json.dumps(dict(W_log2=round(W, 9), claim=str(claim_fr), claim_value=float(claim_fr),
                          claim_minus_W=round(float(claim_fr) - W, 9), stress=round(dbl, 9),
                          gate_margin=round(float(gate) - float(claim_fr), 6),
                          success=measured_success, kappa=kappa,
                          gap_support_log2=DESIGN["gap_support_log2"],
                          abort_log2=round(ABORT_LOG2, 2),
                          E_bad_log2=round(E_bad_log2, 3),
                          rows_priced=len(deriv)), indent=1))


if __name__ == "__main__":
    main()
