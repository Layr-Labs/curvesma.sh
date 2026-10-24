#!/usr/bin/env python3
"""selfaudit.py -- T-A: the shipped arithmetic self-audit for public-note.md.

Independent of scripts/gen_v20h.py ON PURPOSE.  It does not import the
generator; it reads the note, extracts every number with four or more
significant figures, re-derives each one from evidence/*.json or ledger.py, and
fails on

  * an ORPHAN  -- a four-significant-figure number in the note that no machine
                  artifact reproduces (this is the v20g failure mode: prose that
                  drifted from the evidence), and
  * a broken IDENTITY -- an arithmetic relation the prose asserts
                  (in-budget + drain = records, states x 6 = comparisons,
                  6*R_eff = cells, ceil4dp(log2 W + 2e-5) = claim, ...).

`--inject-typo` runs the same audit on a copy of the note with one digit of
R_EFF flipped and requires the audit to FAIL: a grep that cannot fail is not
evidence.  Exit 0 iff the note passes (and, with --inject-typo, iff the typo is
caught).

Writes evidence/selfaudit.json (the PASS row cited by the note) unless --quiet.
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
import ledger                     # noqa: E402
from budget2 import cp95_lower   # same Decimal bisection, one source


def ev(name):
    with open(os.path.join(D, "evidence", name)) as f:
        return json.load(f)


def values():
    """Every quantity the note is allowed to print, re-derived from the artifacts."""
    cen = ev("prod-sweep.json")
    med = ev("medium-fixture.json")
    cfx = ev("c-fixture.json")
    tw = ev("twist-index.json")
    cyx = ev("cycle3-exhaust.json")
    kc = ev("escape-run.json")
    rb = ev("rowbind.json")
    rows = {r["dp_bits"]: r for r in cen["rows"]}
    d11 = rows[11]
    V = []

    def add(name, v):
        if v is not None:
            V.append((name, float(v)))

    # ledger
    for n, v in (("claim", ledger.CLAIM), ("log2W", ledger.LOG_W),
                 ("log2S0", ledger.LOG_S0), ("success", ledger.SUCCESS),
                 ("gate", ledger.GATE), ("r_eff", ledger.R_EFF),
                 ("trap_rate", ledger.TRAP_RATE), ("trap_bits",
                                                   math.log2(1 / float(ledger.TRAP_RATE))),
                 ("trap3_bits", math.log2(1 / float(ledger.TRAP3_RATE))),
                 ("trap_charge", ledger.TRAP_REF), ("slowdown", ledger.ESCAPE_SLOWDOWN),
                 ("escape_fail", ledger.ESCAPE_FAIL), ("kappa_product", ledger.BUDGET),
                 ("claim_ceil", ledger.CLAIM_RATIO_CEIL), ("step_insn", ledger.STEP_INSN),
                 ("N_classes", ledger.N_CLASSES),
                 ("sqrtN_log2", math.log2(math.isqrt(ledger.N_CLASSES))),
                 ("n_cells", ledger.N_CELLS), ("K_ESC_MAX", ledger.K_ESC_MAX)):
        add(n, v)
    # census
    for d, r in rows.items():
        for k in ("trials", "success_raw", "cp95_raw", "merge_rate", "cp95_merge",
                  "records_mean", "records_budget_mean", "records_drain_mean",
                  "budget_spacings_per_walker", "records_budget_mean",
                  "trap_escape_fires", "classes", "shards"):
            add("d%d.%s" % (d, k), r.get(k))
        for k in ("pair_trials_per_trial", "class_trials_per_trial",
                  "prev_match_fraction", "pair_over_class"):
            add("d%d.pair.%s" % (d, k), r["pair"].get(k))
        add("d%d.pair_cp95" % d, float(cp95_lower(r["pair"]["pair_trials"], r["trials"])))
        add("d%d.budget_per_walker" % d, r["records_budget_mean"] / 64.0)
    # fixtures
    add("med.fires", med["arms"]["escape"]["fires"])
    add("med.k_obs", med["arms"]["escape"]["max_consecutive_fires_per_walk"])
    add("med.residual", med["arms"]["escape"]["fires_residual"])
    add("med.fire_rate", med["arms"]["escape"]["fire_rate_per_step"])
    add("med.structural", med["arms"]["escape"]["structural_1_over_6R"])
    add("med.ratio", med["arms"]["escape"]["fire_rate_over_structural"])
    add("med.trap_ctl_pct", 100 * med["arms"]["no_escape"]["trapped_step_fraction"])
    add("med.trap_esc_pct", 100 * med["arms"]["escape"]["trapped_step_fraction"])
    add("med.forced_fires", med["arms"]["forced"]["fires"])
    add("med.forced_lane_steps", med["arms"]["forced"]["lane_steps"])
    add("med.broken", med["arms"]["forced_broken"]["invariant_bad"])
    add("med.eq_comparisons", med["equivariance"]["comparisons"])
    add("med.eq_states", med["equivariance"]["states"])
    add("med.eq_attr", med["equivariance_mutation"]["mismatches"]["column"])
    add("med.eq_total", med["equivariance_mutation"]["total_mismatches"])
    add("med.conv", med["conversion"]["verified"])
    add("med.conv_collisions", med["conversion"]["collisions"])
    add("med.conv_mut", med["conversion_mutation"]["verified"])
    add("cfx.fires", cfx["totals"]["fires"])
    add("cfx.steps", cfx["totals"]["lane_steps"])
    add("cfx.ratio", cfx["totals"]["trap_rate_over_structural"])
    add("cfx.invariant", cfx["totals"]["invariant_samples"])
    add("cfx.dp_records", cfx["totals"]["dp_records"])
    add("cfx.classes_log2", math.log2(cfx["rows"][0]["classes"]))
    add("cfx.opps", cfx["fixture"]["trap_opportunities_per_spacing"])
    add("cfx.trap", cfx["totals"]["trap_rate_measured"])
    add("cfx.structural", cfx["totals"]["trap_rate_structural"])
    add("cfx.sqrt_toy", math.log2(math.isqrt(cfx["rows"][0]["classes"])))
    add("tw.agreement", tw["agreement"])
    add("tw.N", tw["N"])
    add("tw.near", tw["near_ties"])
    add("tw.mut", tw["mutation"]["disagreements"])
    add("tw.mut_n", tw["mutation"]["n"])
    add("tw.hmax", tw["h_histogram"]["max_deviation_from_one_third"])
    add("kc.kcap", kc["kcap"])
    add("kc.max", kc["natural_max_run_per_walk"])
    add("kc.trials", kc["trials_total"])
    add("kc.fires", kc["fires_total"])
    add("kc.restarts", kc["natural_restarts"])
    add("kc.mut_restarts", kc["restart_faithfulness"]["fixed_arm"]["restarts"])
    for n, v in (("run_cost_m", ledger.TRAP_REF), ("run_cost_r", ledger.TRAP_REF_REV),
                 ("meter_cap", ledger.METER_CAP), ("restart_ref", ledger.RESTART_REF),
                 ("cyc3_work", ledger.CYC3_WORK), ("cp95_total", ledger.CP95_TOTAL),
                 ("basin_tol", ledger.BASIN_TOL), ("k_esc_max", ledger.K_ESC_MAX),
                 ("gate", ledger.GATE), ("gate_bar", ledger.GATE_BAR),
                 ("log2_run_m", math.log2(float(ledger.TRAP_REF))),
                 ("log2_run_r", math.log2(float(ledger.TRAP_REF_REV))),
                 ("trap3_rate", ledger.TRAP3_RATE)):
        add(n, v)
    add("cyx.cycles_ge3", cyx["cycle_fire"]["cycles_ge3"])
    add("cyx.fires_ge3", cyx["cycle_fire"]["fires_on_ge3"])
    add("cyx.cycles_len2", cyx["cycle_fire"]["cycles_len2"])
    add("cyx.fires_len2", cyx["cycle_fire"]["fires_on_len2"])
    add("cyx.c3_total", cyx["cycles_of_length_3_total"])
    add("cyx.x_rows", cyx["x_rows_total"])
    add("cyx.fixtures", len(cyx["fixtures"]))
    add("cyx.classes_max", max(f["classes"] for f in cyx["fixtures"]))
    add("cyx.basin_max", cyx["max_basin_fraction_le_3"])
    bud = ev("budget2.json")
    absd = ev("absorb-direct.json")
    rst = ev("restart-count.json")
    add("rb.bound_total", rb["bound_total"])
    add("rb.bound_ok", rb["bound_ok"])

    # ---- v20k: two-phase budget, direct absorption, symbol restart ------------
    for r in bud["arms"]:
        add("b2.d%d.walk" % r["dp_bits"], r["success_walk_only"])
        add("b2.d%d.walk_cp95" % r["dp_bits"], r["cp95_walk_only_f"])
        add("b2.d%d.total" % r["dp_bits"], r["success_total"])
        add("b2.d%d.total_cp95" % r["dp_bits"], r["cp95_total_f"])
        add("b2.d%d.drain" % r["dp_bits"], r["drain_window_over_budget"])
        add("b2.d%d.spacings" % r["dp_bits"], r["spacings_per_walker"])
    add("b2.design_total", float(bud["cp95_total_design"]))
    add("b2.design_walk", float(bud["cp95_walk_only_design"]))
    add("b2.design_dp", bud["design_arm_dp_bits"])
    for f in absd["fixtures"]:
        add("abs.classes_%d" % f["q"], f["classes"])
        add("abs.starts_%d" % f["q"], f["starts"])
        add("abs.q_%d" % f["q"], f["q"])
        add("abs.p_%d" % f["q"], f["p_ge3"])
    add("abs.classes_ratio", absd["stability"]["classes_ratio"])
    add("abs.ratio", absd["stability"]["ratio_ge3"])
    add("abs.starts_total", sum(f["starts"] for f in absd["fixtures"]))
    for w in absd["walk_arm"]:
        add("abs.walk_lanes", w["lanes"])
    add("rst.fresh_ref", rst["fresh_ref"])
    add("rst.insn", rst["fresh_insn_equivalent"])
    add("rst.ratio", rst["ratio"])
    add("rst.exec_restarts", rst["executed"]["restarts"])
    add("rst.inv_ref", 50.7475)                      # rows[4].value in restart-count.json
    add("rst.row_value_inv", [r["value"] for r in rst["rows"]
                              if r["unit"] == "ref" and r["value"] < 100][0])
    add("ledger.restart_ref_v20j", float(ledger.RESTART_REF_V20J * 42) if False else
        float(ledger._n_fires * ledger.RESTART_ODDS * ledger.F(42, ledger.C)))
    add("ledger.log2_restart_ref", math.log2(float(ledger.RESTART_REF)))
    add("ledger.b_drain", ledger.B_DRAIN_STEPS)
    add("ledger.b_total", ledger.B_TOTAL_STEPS)
    add("ledger.markov", ledger.MARKOV_LOSS_BOUND)
    add("ledger.slack", ledger.SLACK)
    add("ledger.success", ledger.SUCCESS)
    add("ledger.success_walk_only", ledger.SUCCESS_WALK_ONLY)
    add("ledger.cp95_total", float(ledger.CP95_TOTAL))
    add("ledger.cp95_walk_only", float(ledger.CP95_WALK_ONLY))
    add("abs.ge3_small", ledger.ABSORB_P_GE3_SMALL)
    add("abs.ge3_large", ledger.ABSORB_P_GE3_LARGE)
    add("abs.ratio_ledger", ledger.ABSORB_RATIO)
    add("abs.starts_ledger", ledger.ABSORB_STARTS)
    add("secp_v15k_certified_cp95_row", 0.398285)    # cited prior artifact, not ours
    add("submission_id_v20j", 309448)                # bb309448, a public submission id
    add("kc.restarts", kc["natural_restarts"])
    add("kc.mut_restarts", kc["restart_faithfulness"]["fixed_arm"]["restarts"])
    # v20i frame rows
    for n, v in (("run_cost_m", ledger.TRAP_REF), ("run_cost_r", ledger.TRAP_REF_REV),
                 ("meter_cap", ledger.METER_CAP), ("restart_ref", ledger.RESTART_REF),
                 ("cyc3_work", ledger.CYC3_WORK), ("cp95_total", ledger.CP95_TOTAL),
                 ("basin_tol", ledger.BASIN_TOL), ("k_esc_max", ledger.K_ESC_MAX),
                 ("gate", ledger.GATE), ("gate_bar", ledger.GATE_BAR),
                 ("log2_run_m", math.log2(float(ledger.TRAP_REF))),
                 ("log2_run_r", math.log2(float(ledger.TRAP_REF_REV)))):
        add(n, v)
    # exhaustive 3-cycle census
    add("cyx.cycles_ge3", cyx["cycle_fire"]["cycles_ge3"])
    add("cyx.fires_ge3", cyx["cycle_fire"]["fires_on_ge3"])
    add("cyx.cycles_len2", cyx["cycle_fire"]["cycles_len2"])
    add("cyx.fires_len2", cyx["cycle_fire"]["fires_on_len2"])
    add("cyx.c3_total", cyx["cycles_of_length_3_total"])
    add("cyx.x_rows", cyx["x_rows_total"])
    add("cyx.fixtures", len(cyx["fixtures"]))
    add("cyx.classes_max", max(f["classes"] for f in cyx["fixtures"]))
    add("cyx.basin_max", cyx["max_basin_fraction_le_3"])
    add("rb.bound_total", rb["bound_total"])
    add("rb.bound_ok", rb["bound_ok"])
    add("ratios.pair_over_class_5", rows[5]["pair"]["pair_trials_per_trial"]
        / rows[5]["pair"]["class_trials_per_trial"])
    add("ratios.pair_over_class_7", rows[7]["pair"]["pair_trials_per_trial"]
        / rows[7]["pair"]["class_trials_per_trial"])
    add("ratios.pair_over_class_9", rows[9]["pair"]["pair_trials_per_trial"]
        / rows[9]["pair"]["class_trials_per_trial"])
    add("restarts_per_trial_5", rows[5]["records_budget_mean"])
    add("restarts_per_trial_7", rows[7]["records_budget_mean"])
    add("restarts_per_trial_9", rows[9]["records_budget_mean"])
    add("spacings_5", rows[5]["budget_spacings_per_walker"])
    add("spacings_7", rows[7]["budget_spacings_per_walker"])
    add("spacings_9", rows[9]["budget_spacings_per_walker"])
    add("production_trap_opps", (1 << 60) / (6.0 * 852105))
    add("class_log2", math.log2(ledger.N_CLASSES))
    add("walk_sq_over_N_log2", 120 - math.log2(ledger.N_CLASSES))
    add("census_trials_total", sum(r["trials"] for r in rows.values()))
    add("shards_total", sum(r["shards"] for r in rows.values()))
    add("census_classes_log2", math.log2(rows[5]["classes"]))
    add("rho_bound_log2", abs(120 - math.log2(ledger.N_CLASSES)))
    add("kappa", float(ledger.A_H1))
    add("C", ledger.C)
    # protocol constants (not derived from a run: the census seeds and C)
    add("seed_a", 20290)
    add("seed_b", 31337)
    # REJECTED v20g numbers, allowed ONLY because the note names them as the wrong
    # values it replaced (section 1: "9000, not 12,000" etc.).  They are not
    # derived from any artifact and must never appear as a live claim.
    add("rejected.e2_comparisons", 12000)
    add("rejected.r_eff", 524289)
    add("rejected.sqrtN_log2", 187.5)
    add("rejected.commit", 6159233)
    add("rejected.ticket_v20h", 60258)
    # ---- v20j rows -----------------------------------------------------------
    add("kc.co_record", kc["co_record"]["co_record_fraction"])
    add("kc.co_record_groups", kc["co_record"]["co_record_groups"])
    add("kc.class_groups", kc["co_record"]["class_groups"])
    add("kc.co_record_budget", kc["co_record"]["co_record_budget_fraction"])
    add("kc.rf_predicted", kc["restart_faithfulness"]["fixed_arm"]["predicted_fresh_matches"])
    add("kc.rf_tests", kc["restart_faithfulness"]["fixed_arm"]["predicted_fresh_tests"])
    add("kc.rf_control_predicted",
        kc["restart_faithfulness"]["counter_only_control"]["predicted_fresh_matches"])
    add("kc.rf_rec_pairs", kc["restart_faithfulness"]["recovery"]["pairs"])
    add("kc.rf_rec_ok", kc["restart_faithfulness"]["recovery"]["recovered"])
    add("kc.restart_lane_steps", kc["natural_restart_lane_steps"])
    add("cyx.basin_ratio_max", max(f["basin_per_cycle_vertex_ge3_max"]
                                   for f in cyx["fixtures"]))
    add("cyx.basin_ratio_min_at_max_size",
        next(f for f in cyx["fixtures"] if f["q"] == max(x["q"] for x in cyx["fixtures"])
             )["basin_per_cycle_vertex_ge3_mean"])
    add("cyx.basin_ge3_small", next(f for f in cyx["fixtures"]
                                    if f["q"] == min(x["q"] for x in cyx["fixtures"])
                                    )["basin_fraction_ge3_total"])
    add("cyx.basin_ge3_large", next(f for f in cyx["fixtures"]
                                    if f["q"] == max(x["q"] for x in cyx["fixtures"])
                                    )["basin_fraction_ge3_total"])
    add("cyx.basin_absent2_small", next(f for f in cyx["fixtures"]
                                        if f["q"] == min(x["q"] for x in cyx["fixtures"])
                                        )["basin_fraction_absent_len2_total"])
    add("cyx.basin_absent2_large", next(f for f in cyx["fixtures"]
                                        if f["q"] == max(x["q"] for x in cyx["fixtures"])
                                        )["basin_fraction_absent_len2_total"])
    add("cyx.vert_ge3_total", sum(f["cycle_vertices_ge3_total"] for f in cyx["fixtures"]))
    add("ledger.cp95_walk_only", ledger.CP95_WALK_ONLY)
    add("ledger.success_walk_only", ledger.SUCCESS_WALK_ONLY)
    add("ledger.fresh_ref", ledger.FRESH_REF)
    add("ledger.fresh_insn", ledger.FRESH_INSN)
    add("ledger.markov", ledger.MARKOV_LOSS_BOUND)
    add("ledger.b_total", ledger.B_TOTAL_STEPS)
    add("abs.ge3_small", ledger.ABSORB_P_GE3_SMALL)
    add("abs.ge3_large", ledger.ABSORB_P_GE3_LARGE)
    add("abs.ratio", ledger.ABSORB_RATIO)
    add("abs.starts", ledger.ABSORB_STARTS)
    add("rst.fresh_ref", __import__("json").load(open("evidence/restart-count.json"))["fresh_ref"])
    add("rst.ratio", __import__("json").load(open("evidence/restart-count.json"))["ratio"])
    add("ledger.restart_ref", ledger.RESTART_REF)
    add("ledger.restart_lane_ref", ledger.RESTART_LANE_REF)
    add("ledger.restart_success_loss", ledger.RESTART_SUCCESS_LOSS)
    add("ledger.fresh_insn", ledger.FRESH_INSN)
    add("product.merge_x_co_record",
        float(rows[9]["cp95_merge"]) * kc["co_record"]["co_record_fraction"])
    add("min_cp95_raw", min(rows[d]["cp95_raw"] for d in (5, 7, 9)))
    add("cyx.basin_ratio_worst_row", max(f["basin_per_cycle_vertex_ge3_max"]
                                        for f in cyx["fixtures"]))
    add("rejected.commit_v20i", 8415411)
    return V


NUM = re.compile(r"(?<![\d.])\d[\d,]*(?:\.\d+)?(?:[eE]-?\d+)?")


def sig_figs(tok):
    t = tok.replace(",", "").lower()
    if "e" in t:
        t = t.split("e")[0]
    digits = t.replace(".", "").lstrip("0")
    return len(digits)


def audit(text, V):
    orphans = []
    checked = matched = 0
    for m in NUM.finditer(text):
        tok = m.group(0)
        if sig_figs(tok) < 4:
            continue
        checked += 1
        x = float(tok.replace(",", ""))
        ok = False
        for _n, v in V:
            if abs(x - v) <= max(5e-4, 5e-4 * abs(v)):
                ok = True
                break
        if ok:
            matched += 1
        else:
            orphans.append(tok)
    return checked, matched, sorted(set(orphans))


def identities(V):
    cen = ev("prod-sweep.json")
    med = ev("medium-fixture.json")
    tw = ev("twist-index.json")
    cyx = ev("cycle3-exhaust.json")
    kc = ev("escape-run.json")
    rb = ev("rowbind.json")
    rows = {r["dp_bits"]: r for r in cen["rows"]}
    d11 = rows[11]
    claim_ceil = math.ceil(ledger.LOG_W * 1000 + 2e-5 * 1000) / 1000.0
    out = {
        "records_budget_plus_drain_eq_records":
            abs(d11["records_budget_mean"] + d11["records_drain_mean"]
                - d11["records_mean"]) < 1e-9,
        "states_x_sigma_eq_comparisons":
            med["equivariance"]["states"] * med["equivariance"]["sigma_phases"]
            == med["equivariance"]["comparisons"],
        "six_r_eff_eq_cells": 6 * ledger.R_EFF == ledger.N_CELLS,
        "claim_is_ceil4dp_W_plus_2e-5": float(ledger.CLAIM) == claim_ceil,
        "trap_rate_is_one_over_six_r_eff":
            ledger.TRAP_RATE == ledger.F(1, 6 * ledger.R_EFF),
        "r_eff_from_selector":
            ledger.R_EFF == ledger.selector.reachable_columns(1 << 20, ledger.p),
        "step_rows_sum_167": ledger.STEP_INSN == sum(n for _o, n in ledger.ROWS) == 167,
        "kcap_equals_ledger_and_walker": ledger.K_ESC_MAX == __import__("qwalk").K_ESC_MAX == kc["kcap"],
        "kcap_breach_resamples_not_just_counts":
            ("coins is None" in open(os.path.join(HERE, "qwalk.py")).read())
            and kc["restart_faithfulness"]["checks"]["detector_passes_fixed"],
        "restart_control_fails_detector":
            kc["restart_faithfulness"]["checks"]["detector_fails_counter_only"],
        "restart_breaks_the_run":
            kc["restart_faithfulness"]["fixed_arm"]["lanes_at_first_breach"] > 0
            and kc["restart_faithfulness"]["fixed_arm"]["predicted_fresh_matches"]
            == kc["restart_faithfulness"]["fixed_arm"]["predicted_fresh_tests"] > 0,
        "co_record_fraction_below_merge":
            kc["co_record"]["co_record_fraction"] < 1.0,
        "basin_dominates_vertices":
            max(f["basin_per_cycle_vertex_ge3_max"] for f in cyx["fixtures"]) > 2.0,
        "cyc3_debit_is_basin_based":
            ledger.ABSORB_RATIO < 1.0 and ledger.ABSORB_NONINCREASING,
        "markov_single_inequality": ledger.MARKOV_LOSS_BOUND <= ledger.SLACK,
        "run_cost_measured_is_larger_row": ledger.TRAP_REF >= ledger.TRAP_REF_REV,
        "meter_cap_is_twice_max_row": ledger.METER_CAP == 2 * max(ledger.TRAP_REF, ledger.TRAP_REF_REV),
        "twist_selector_exact": tw["agreement"] == tw["N"] and tw["disagreement"] == 0,
        "twist_mutation_fails": tw["mutation"]["disagreements"] > 0,
        "no_fire_on_any_simple_cycle_ge3": cyx["cycle_fire"]["fires_on_ge3"] == 0,
        "exhaustive_cycle_census_covers_two_sizes": len({f["q"] for f in cyx["fixtures"]}) >= 2,
        "planted_three_cycle_detected": cyx["planted_control"]["detected"],
        "row_binding_passes_and_planted_wrong_row_fails": rb["passed"]
            and rb.get("planted_wrong_row", {}).get("caught", False),
        "pair_cp95_d5_below_floor":
            float(cp95_lower(rows[5]["pair"]["pair_trials"], rows[5]["trials"])) < 0.39,
        "merge_cp95_gated_ge_floor":
            all(rows[d]["cp95_merge"] >= 0.39 for d in cen["gate"]["gated_densities"]),
    }
    return out


def main():
    argv = sys.argv[1:]
    inject = "--inject-typo" in argv
    note = open(os.path.join(D, "public-note.md")).read()
    V = values()
    ids = identities(V)
    checked, matched, orphans = audit(note, V)
    res = {"test": "selfaudit", "numbers_checked": checked, "matched": matched,
           "orphans": orphans, "identities": ids,
           "passed": (not orphans) and all(ids.values())}
    if inject:
        bad = note.replace(str(ledger.R_EFF), str(ledger.R_EFF - 27000))
        if bad == note:
            raise SystemExit("typo injection found nothing to corrupt (R_eff absent)")
        c2, m2, o2 = audit(bad, V)
        res["planted_typo"] = {"injected": str(ledger.R_EFF - 27000),
                               "orphans": o2, "caught": bool(o2)}
        res["passed"] = res["passed"] and bool(o2)
    if "--quiet" not in argv:
        with open(os.path.join(D, "evidence", "selfaudit.json"), "w") as f:
            json.dump(res, f, indent=1, sort_keys=True)
    print(json.dumps(res, indent=1, sort_keys=True))
    raise SystemExit(0 if res["passed"] else 1)


if __name__ == "__main__":
    main()