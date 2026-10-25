"""v23 ledger: claim computed once here; the submission cites this file.

WHAT MOVES FROM v22 (verdict v22 fire 2's four numbered flaws, verbatim there):

1. ROW != KEY.  v22 measured 20-bit selector-row occupancy.  v23 measures the
   SHIPPED COLLISION PREDICATE -- bucket probe then canonical-key compare -- on
   a real prime-order curve small enough for COMPLETE x-coordinate keys to
   collide (p=65519, #E=65447, r=32759 sign classes; experiments/key_equal_v23.py).
   Per-pair key hazard measured 2.95e-5 vs 1/r = 3.05e-5; bucket hazard
   measured 9.9e-4 vs 1/R = 9.77e-4; the conditional P(key | same bucket) is
   measured and compared with R/r.  q (recovery per collision) is re-measured
   on the key predicate through the probe path: 3373/3373 nontrivial = 1.0,
   including 1729 the OPPOSITE-point (T = -T') class that a difference-only
   replay would miscount.

2. CARDINALITY != UNION BOUND.  No |distinct rows| is used as an event
   probability anywhere.  The success charge is:
     * Paley-Zygmund P(X>=1) >= E[X]^2/E[X^2] with E[X] = N(N-1)/(2r) and
       E[X^2] = E[X(X-1)] + E[X]; for the Poisson moments this is 0.387584;
     * the Poisson value 1-exp(-E[X]) = 0.468939;
     * the measured factorial moments E[X(X-1)]/E[X]^2 = 0.9766 (carrier),
       1.0236, 0.9877 across the three cells, i.e. Poisson-tight;
     * CHARGED: the minimum Clopper-Pearson 95% lower bound over the three
       seeded key-collision cells -- a lower bound measured on the actual
       predicate, which already carries whatever deficit the fixture's key law
       shows against the uniform model (measured mean 0.590 vs lambda 0.633).

3. LOCATOR ACCESS PRICED.  store_v23.py executes the 128-bit anchor chain and
   charges 4 probes x 6 levels x 10 Ir = 240 Ir/record as its own row
   (`locator_access`), on top of the unchanged v22 sort price.  The page-width
   control moves the locator row 246.518 -> 243.259 Ir/record.

4. WORK.  Rows: the v22 row set plus the locator row.  Claim is the ceiling of
   the exact log2 work at 1e-4, reduced rational, one-ulp check in audit_v23.

Everything else (affine addition, batch inversion, table build, returns
first-moment bound, replay charge, 128-bit split/join, peak) is retained.
"""
import json
import math
import os
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "evidence")

LOG2_V = 127.67
WALKERS = 1 << 20
STEPS = 2 ** (LOG2_V - 20)


def main():
    price = json.load(open(os.path.join(OUT, "insn_price.json")))
    fe = price["prices_ir"]["fe_mul"]
    fa = price["prices_ir"]["fe_add"]
    fs = price["prices_ir"]["fe_sub"]
    store = json.load(open(os.path.join(OUT, "store_v23.json")))
    ret = json.load(open(os.path.join(OUT, "returns_v23.json")))
    key = json.load(open(os.path.join(OUT, "key_equal_v23.json")))

    sort_ir = store["sort_total_ir_per_record"]
    locator_ir = store["locator_access"]["locator_access_ir_per_record"]
    # the locator-access row is its own charge; the sort price above INCLUDES it,
    # so keep the two visible without double-charging:
    sort_only = store["sort_total_ir_per_record_v22"]
    rows = {
        "affine_add": 3 * fe + 4 * fa + 2 * fs,
        "batch_invert_share": 3 * fe,
        "batch_chain_share": (385 * fe) / (1 << 20),
        "table_read_index": 10,
        "canon_sign_rule": 20,
        "record_append_48B": 12,
        "global_sort_and_witness": sort_only,
        "locator_access_128bit_x4probes": locator_ir,
    }
    per_visit = sum(rows.values())
    assert abs(per_visit - (sort_only + locator_ir + rows["affine_add"]
                            + rows["batch_invert_share"] + rows["batch_chain_share"]
                            + 42)) < 1e-9
    units = per_visit / 1000.0
    log2_rows = LOG2_V + math.log2(units)

    # ---------------- success: measured key predicate -------------------------
    cells = key["carrier_cells"]
    carrier = cells[0]
    p_cp = min(c["cp95_lower"] for c in cells)
    lam_meas = carrier["mean"]
    lam_model = carrier["lambda_expected"]
    prod = key["production"]
    lam_prod = prod["lambda"]
    # measured-hazard-adjusted production rate (fixture deficit carried through)
    lam_prod_meas = lam_prod * (lam_meas / lam_model)
    success_poisson_adj = 1.0 - math.exp(-lam_prod_meas)
    pz = prod["paley_zygmund_bound"]
    pz_meas = (lam_meas ** 2) / (lam_meas ** 2 + carrier["moments"]["E_XX1"] + lam_meas)
    charged = min(p_cp, success_poisson_adj)
    attempts = 1.0 if charged >= 0.39 else math.log(1 - 0.39) / math.log(1 - charged)
    assert charged >= 0.39, "charged success %.6f misses the 0.39 floor" % charged
    log2_attempts = math.log2(attempts) if attempts > 1 else 0.0

    lam_coll = carrier["nontrivial"] / (carrier["trials"] * 1.0)
    full_step = 388 * fe + 4 * fa + 2 * fs + 30
    bulk_step = rows["affine_add"] + rows["batch_invert_share"] + 30
    infl = full_step / bulk_step
    trails_per_challenge = 2.0 * (1.0 + lam_coll)
    replay_rel = trails_per_challenge * infl / WALKERS

    rel_restart = ret["production_charge"]["charged_relative_loss"]
    log2_work = log2_rows + math.log2(1.0 + rel_restart) + math.log2(1.0 + replay_rel)
    table_ir = 2 * (1 << 20) * 256 * 3 * fe
    log2_work += math.log2(1.0 + table_ir / 2 ** log2_work)
    log2_work += log2_attempts

    claim = math.ceil(log2_work * 1e4) / 1e4
    assert claim >= log2_work
    num = int(round(claim * 10000))
    g = gcd(num, 10000)
    rational = {"numerator": num // g, "denominator": 10000 // g}

    ev = {
        "rows_ir_per_visit": {k: round(v, 6) for k, v in rows.items()},
        "instructions_per_visit": per_visit,
        "units_per_visit": units,
        "locator": {
            "row_ir_per_record": locator_ir,
            "row_including_page_walk_ir_per_record": store["locator_row_ir_per_record"],
            "per_level_ir": store["locator_access"]["per_level_ir"],
            "levels": store["locator_access"]["levels"],
            "probes_per_record": store["locator_access"]["probes_per_record"],
            "executed_counters": store["locator_access"]["executed_counters"],
            "page_width_control": store["page_width_control"],
            "source": store["locator_access"]["probe_path_source"],
            "note": "this row is charged IN ADDITION to the v22 sort price (which contains "
                    "the per-pass page walk); the row's page-width sensitivity is shown by "
                    "the control, and the ledger charges locator_ir exactly once",
        },
        "sort": {
            "per_pass_ir_per_record": store["per_pass_total_ir"],
            "passes": store["passes"],
            "witness_ir": store["witness_ir_per_record"],
            "directory_construction_ir": store["directory_construction_ir_per_record"],
            "total_sort_only_ir_per_record": sort_only,
            "total_including_locator_ir_per_record": sort_ir,
            "peaks": store["peaks"],
        },
        "success": {
            "floor": 0.39,
            "predicate": "KEY EQUALITY: bucket = key & (R-1); probe reaches entry; "
                         "canonical x compares equal; the replay resolves T=T' vs T=-T' "
                         "and recovers x by the matching linear relation",
            "carrier_cells": cells,
            "key_hazard_measured_per_pair": carrier["key_pairs_per_pair"],
            "key_hazard_model_per_pair": 1.0 / key["r_sign_classes"],
            "bucket_hazard_measured_per_pair": carrier["bucket_hazard_per_pair"],
            "bucket_hazard_model_per_pair": 1.0 / carrier["R"],
            "conditional_key_given_bucket_measured": carrier["conditional_key_given_bucket"],
            "conditional_key_given_bucket_model_R_over_r": carrier["conditional_model_R_over_r"],
            "yield_recovery_given_collision": carrier["yield_nontrivial"],
            "yield_basis": "%d/%d nontrivial key-equal events recovered; %d same-point, "
                           "%d opposite-point, %d trivial duplicates"
                           % (carrier["recovered"], carrier["nontrivial"],
                              carrier["same_point"], carrier["opposite_point"],
                              carrier["trivial"]),
            "production_model": {
                "records_log2": LOG2_V, "sign_classes_log2": 255,
                "E_X": lam_prod, "E_XX1": prod["E_XX1"], "E_X2": prod["E_X2"],
                "factorial_ratio_model": prod["factorial_ratio_model"],
                "paley_zygmund_bound": pz,
                "poisson_success": prod["poisson_success"],
                "union_bound_upper": prod["union_bound_upper"],
                "measured_to_model_rate_ratio": lam_meas / lam_model,
                "lambda_at_measured_rate": lam_prod_meas,
                "poisson_success_at_measured_rate": success_poisson_adj,
            },
            "charged_success": charged,
            "charged_basis": "minimum CP-95 lower bound over the three seeded key-collision "
                             "cells (%d/%d, %d/%d, %d/%d), i.e. a lower bound measured on "
                             "the actual predicate, not a cardinality count"
                             % tuple(sum(([c["wins"], c["trials"]] for c in cells), [])),
            "attempts_charged": attempts,
            "log2_attempts_charged": log2_attempts,
            "moments_note": "E[X(X-1)]/E[X]^2 measured 0.9766 / 1.0236 / 0.9877 vs 1.0 for "
                            "Poisson; Paley-Zygmund at production moments is 0.387584 and is "
                            "reported as the distribution-free floor, below the measured CP "
                            "lower bound that is charged",
            "residual_conjecture": key["residual_conjecture"],
        },
        "returns": {
            "relative_loss": rel_restart,
            "unit": ret["unit"],
            "basis": ret["production_charge"]["basis"],
            "state_space": ret["production_charge"]["state_space"],
            "high_load_cell": {k: ret["cells"]["high_load_R64_S8192"][k] for k in
                               ("seeds", "cycle_entries", "cycle_entry_rate",
                                "first_moment_bound_per_trail", "entry_rate_consistent_with_bound",
                                "repeats_total_descriptive", "fatal_repeats", "identity_breaks")},
            "v21_binomial_retracted": ret["production_charge"]["v21_binomial_retracted"],
        },
        "policy": {
            "record_fields": "32B canonical key || 20B identity pack (20-bit trail || 108-bit step)",
            "rule": "every hit replays both trails from their stored identities; the replay "
                    "reconstructs the exact point (resolving the sign class) and recovers x "
                    "from the matching linear relation",
            "measured_collisions_per_challenge": lam_coll,
            "measured_useful_fraction": carrier["yield_nontrivial"],
            "trail_pairs_charged_per_challenge": trails_per_challenge,
            "replay_inflation_vs_bulk": infl,
            "replay_relative_loss": replay_rel,
        },
        "work": {"log2_visits": LOG2_V, "walkers": WALKERS, "steps_per_walker": STEPS,
                 "log2_rows": log2_rows, "log2_work_final": log2_work},
        "bounds": {"work_log2": claim, "claim": claim,
                   "storage_log2_peak": store["peaks"]["log2_peak_bytes"]},
        "claim_time_log2": claim,
        "claim_reduced_rational": rational,
        "delta_vs_v22": {
            "v22_claim": 128.9035, "v23_claim": claim,
            "locator_access": "+%d Ir/record (4 probes x 6 levels x 10 Ir)" % locator_ir,
            "success": "row-carrier CP 0.424555 -> key-predicate CP %.6f" % charged,
            "yield": "q re-measured on the key predicate: %.1f (%d/%d nontrivial)"
                     % (carrier["yield_nontrivial"], carrier["recovered"], carrier["nontrivial"]),
        },
    }
    json.dump(ev, open(os.path.join(OUT, "ledger_v23.json"), "w"), indent=1, sort_keys=True)
    json.dump({"claim": claim, "rational": rational, "log2_work": log2_work},
              open(os.path.join(OUT, "claim_v23.json"), "w"), indent=1, sort_keys=True)
    print(json.dumps({"per_visit": round(per_visit, 3), "locator": locator_ir,
                      "log2_work": log2_work, "claim": claim, "rational": rational,
                      "charged_success": charged, "attempts": attempts,
                      "key_hazard": carrier["key_pairs_per_pair"],
                      "model_key": 1.0 / key["r_sign_classes"],
                      "cond": carrier["conditional_key_given_bucket"],
                      "cond_model": carrier["conditional_model_R_over_r"],
                      "yield": carrier["yield_nontrivial"]}, indent=1))


if __name__ == "__main__":
    main()
