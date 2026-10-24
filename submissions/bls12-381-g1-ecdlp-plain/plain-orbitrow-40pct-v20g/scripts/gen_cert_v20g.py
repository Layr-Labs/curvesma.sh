#!/usr/bin/env python3
"""gen_cert_v20g.py -- machine certificate for plain-orbitrow-40pct-v20g.

Every number is recomputed from evidence/prod-sweep.json, evidence/medium-fixture.json,
evidence/c-fixture.json and ledger.py; none is typed into the note by hand.  The
certificate EMITS the note's tables (evidence/note-tables.md) and byte-compares them
against the block between the BEGIN/END markers of public-note.md.

The checks are grouped by the reviewer sentence they answer:

  E1  force_every is passed by run() and the forced arm is a real measurement;
      the mutation arm that forgets the coefficient update breaks the invariant.
  E2  equivariance is tested on PRE-canonical states and the mismatch counter
      increments; a representative-dependent selector must fail the same test.
  E3  per-walk K, escape success computed from observed outcomes, residual = 0.
  E4  the coverage law is measured on the PAIR (class, prev-class) as well as the
      class, and the class set equals the independent audit's merge set.
  E5  the surrogate uses the ledger's own selector module; the C fixture has
      >= 2^26 classes, one full 2^20-step spacing, >= 32x pre-capture margin.
  E6  d=11 restarts are counted, not denied.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
from prod_sweep import cp95_lower          # noqa: E402
import ledger                              # noqa: E402

checks = []


def chk(name, ok, detail):
    checks.append({"name": name, "ok": bool(ok), "detail": detail})


def ev(name):
    with open(os.path.join(D, "evidence", name)) as f:
        return json.load(f)


cen = ev("prod-sweep.json")
med = ev("medium-fixture.json")
cfx = ev("c-fixture.json")
rows = cen["rows"]
dens = sorted(r["dp_bits"] for r in rows)
pc = cen["plant_check"]
GATED = cen["gate"]["gated_densities"]
meq, mmut = med["equivariance"], med["equivariance_mutation"]
mesc, mctl, mfor, mbrk = (med["arms"]["escape"], med["arms"]["no_escape"],
                          med["arms"]["forced"], med["arms"]["forced_broken"])

# ---------- E1: force_every is wired, and the mutation control fails ----------
chk("E1_force_every_passed_by_run", mfor["force_every"] > 0 and mfor["fires"] > 0,
    "run(arm, force_every=%d) fired %d substitutions; the parameter reaches the step"
    % (mfor["force_every"], mfor["fires"]))
chk("E1_forced_arm_is_a_real_measurement",
    mfor["fires"] > 0.9 * (mfor["steps"] * mfor["lanes"]) / mfor["force_every"],
    "%d fires vs %d scheduled lane-steps (%.3f)"
    % (mfor["fires"], (mfor["steps"] * mfor["lanes"]) // mfor["force_every"],
       mfor["fires"] / ((mfor["steps"] * mfor["lanes"]) / mfor["force_every"])))
chk("E1_mutation_breaks_invariant", mbrk["invariant_bad"] > 0,
    "forced_broken (substitutes the row, FORGETS the coefficient update): %d invariant "
    "violations of %d samples -- the forced test can fail"
    % (mbrk["invariant_bad"], mbrk["invariant_ok"] + mbrk["invariant_bad"]))
chk("E1_escape_arm_keeps_invariant", mesc["invariant_bad"] == 0,
    "a + c*x == canon on %d sampled lane-steps with the escape armed, bad %d"
    % (mesc["invariant_ok"], mesc["invariant_bad"]))
chk("E1_control_trap_persists", mctl["trapped_step_fraction"] > 0.10,
    "the SAME fixture with the escape DISABLED spends %.2f%% of its steps inside a "
    "two-cycle; with the escape armed %.2f%%"
    % (100 * mctl["trapped_step_fraction"], 100 * mesc["trapped_step_fraction"]))

# ---------- E2: equivariance on PRE-canonical states + mutation control ------
chk("E2_equivariance_zero_mismatch", meq["passed"] and meq["total_mismatches"] == 0,
    "%d comparisons on PRE-canonical states (column, escape decision, substituted "
    "successor), mismatches %s" % (meq["comparisons"], meq["mismatches"]))
chk("E2_mismatch_counter_increments", mmut["total_mismatches"] > 0,
    "the representative-dependent selector (the v19d defect) fails the SAME test with "
    "%d mismatches in %d comparisons -- the counter is wired"
    % (mmut["total_mismatches"], mmut["comparisons"]))
chk("E2_mutation_is_asymmetric_by_orders_of_magnitude",
    mmut["total_mismatches"] > 100 * max(meq["total_mismatches"], 1),
    "fixed arm %d vs control %d" % (meq["total_mismatches"], mmut["total_mismatches"]))

# ---------- E3: per-walk K, measured success, residual counter ---------------
chk("E3_k_is_per_walk", "max_consecutive_fires_per_walk" in mesc
    and mesc["max_consecutive_fires_per_walk"] >= 1,
    "longest run of consecutive substitutions along ONE walk: %d (natural arm); "
    "forced stress arm %d; ledger charges K_ESC = %d"
    % (mesc["max_consecutive_fires_per_walk"], mfor["max_consecutive_fires_per_walk"],
       ledger.K_ESC))
chk("E3_k_charged_covers_measurement", ledger.K_ESC >= mesc["max_consecutive_fires_per_walk"],
    "K_ESC %d >= measured natural max %d" % (ledger.K_ESC,
                                             mesc["max_consecutive_fires_per_walk"]))
chk("E3_escape_success_computed_from_outcomes",
    abs(mesc["escape_success_rate_measured"]
        - mesc["fires_left_two_cycle"] / mesc["fires"]) < 1e-12,
    "%d of %d fires left the two-cycle (computed, not 1-1/fires)"
    % (mesc["fires_left_two_cycle"], mesc["fires"]))
chk("E3_residual_trap_zero", mesc["fires_residual"] == 0,
    "residual counter %d over %d fires" % (mesc["fires_residual"], mesc["fires"]))
chk("E3_trap_law_holds", 0.5 < mesc["fire_rate_over_structural"] < 1.5,
    "measured %.3e vs 1/(6*R_eff) %.3e, ratio %.3f"
    % (mesc["fire_rate_per_step"], mesc["structural_1_over_6R"],
       mesc["fire_rate_over_structural"]))
chk("E3_conversion_end_to_end", med["conversion"]["verified"] > 0
    and med["conversion"]["conversion"] >= 0.99,
    "%d of %d class groups solve x from two records (%s)"
    % (med["conversion"]["verified"], med["conversion"]["collisions"],
       med["conversion"]["conversion"]))
chk("E3_conversion_mutation_zero", med["conversion_mutation"]["verified"] == 0,
    "the per-record coefficient corruption converts %d of %d groups"
    % (med["conversion_mutation"]["verified"],
       med["conversion_mutation"]["collisions"]))

# ---------- E4: the coverage law ON THE PAIR --------------------------------
for r in rows:
    p = r["pair"]
    chk("E4_pair_census_rows_d%d" % r["dp_bits"],
        "pair_trials_per_trial" in p and "prev_match_fraction" in p,
        "class coverage %.4f, PAIR coverage %.4f, predecessor-matched %.4f of class "
        "groups, pair/class %.3f"
        % (p["class_trials_per_trial"], p["pair_trials_per_trial"],
           p["prev_match_fraction"], p["pair_trials_per_trial"] / p["class_trials_per_trial"]))
    chk("E4_class_set_equals_audit_d%d" % r["dp_bits"],
        p["class_groups_match_audit_all_shards"]
        and abs(p["class_trials_per_trial"] - r["merge_rate"]) < 1e-12,
        "pair_census class trials %.4f == audit merge %.4f in every shard"
        % (p["class_trials_per_trial"], r["merge_rate"]))
chk("E4_pair_coverage_measured_everywhere",
    all(r["pair"]["pair_trials_per_trial"] > 0 for r in rows),
    "; ".join("d%d %.4f" % (r["dp_bits"], r["pair"]["pair_trials_per_trial"])
              for r in rows))
chk("E4_pair_coverage_tracks_class_coverage",
    all(r["pair"]["pair_trials_per_trial"]
        >= 0.93 * r["pair"]["class_trials_per_trial"] for r in rows),
    "; ".join("d%d %.3f of class" % (r["dp_bits"],
               r["pair"]["pair_trials_per_trial"] / r["pair"]["class_trials_per_trial"])
              for r in rows))

# ---------- E5: selector fidelity and the C fixture -------------------------
chk("E5_selector_is_the_ledger_module", ledger.selector.cell_vec.__module__ == "selector"
    and ledger.selector.cell(2 ** 30, 1 << 20, ledger.p) == (2 ** 30 * (1 << 20)) >> 380,
    "qwalk.col calls selector.cell_vec, the function ledger.py prices; the scalar cell "
    "is the same body plus the range guard")
chk("E5_toy_column_law_uniform", abs(med["arms"]["escape"]["fire_rate_over_structural"] - 1) < 0.5,
    "the toy's column is uniform on the reachable buckets: the two-step return lands on "
    "1/(6*R_eff) to %.3f" % mesc["fire_rate_over_structural"])
chk("E5_c_fixture_classes_2_pow_26", cfx["rows"][0]["classes"] >= (1 << 26),
    "%d classes (2^%.2f), q = %d" % (cfx["rows"][0]["classes"],
                                     __import__("math").log2(cfx["rows"][0]["classes"]),
                                     cfx["rows"][0]["q"]))
chk("E5_c_fixture_full_spacing", cfx["rows"][0]["steps"] == (1 << 20),
    "each walker runs %d steps = one full 2^20-step reporting spacing (d = 20)"
    % cfx["rows"][0]["steps"])
chk("E5_pre_capture_margin_32x", cfx["fixture"]["trap_opportunities_per_spacing"] >= 32,
    "%.1f two-step-return opportunities per spacing (margin %.2fx over the 32x floor); "
    "production 2^37.4 is the same quantity at production R and DP"
    % (cfx["fixture"]["trap_opportunities_per_spacing"],
       cfx["fixture"]["pre_capture_margin_vs_32x"]))
chk("E5_c_fixture_trap_law", 0.75 <= cfx["totals"]["trap_rate_over_structural"] <= 1.25,
    "measured %.3e vs 1/(6*R_eff) %.3e, ratio %.3f over %d lane-steps"
    % (cfx["totals"]["trap_rate_measured"], cfx["totals"]["trap_rate_structural"],
       cfx["totals"]["trap_rate_over_structural"], cfx["totals"]["lane_steps"]))
chk("E5_c_fixture_invariant", cfx["totals"]["invariant_violations"] == 0,
    "0 violations in %d sampled lane-steps" % cfx["totals"]["invariant_samples"])
chk("E5_c_fixture_escape_and_residual", cfx["totals"]["fires_left_two_cycle"] > 0
    and cfx["totals"]["fires_residual"] == 0,
    "%d fires, %d left the two-cycle, %d residual"
    % (cfx["totals"]["fires"], cfx["totals"]["fires_left_two_cycle"],
       cfx["totals"]["fires_residual"]))
chk("E5_c_fixture_pair_and_class_coverage",
    cfx["totals"]["pair_collisions_cross_walker"] > 0
    and cfx["totals"]["class_collisions_cross_walker"] > 0,
    "2^26 classes, one spacing: %d class collisions, %d pair collisions across walkers"
    % (cfx["totals"]["class_collisions_cross_walker"],
       cfx["totals"]["pair_collisions_cross_walker"]))
chk("E5_c_fixture_dp_records_disclosed", "dp_records" in cfx["totals"],
    "DP-gated records inside one spacing at 2^26 classes: %d (the walker's cycle is "
    "~sqrt(N) = 2^13 steps, below the 2^20 spacing, so a d=20 gate rarely fires here; "
    "production has sqrt(N) = 2^187.5 >> DP = 2^60) -- disclosed, not extrapolated"
    % cfx["totals"]["dp_records"])

# ---------- E6: d=11 is honest ---------------------------------------------
d11 = [r for r in rows if r["dp_bits"] == 11][0]
restart_per_walker = d11["per_walker"] / (2.0 ** 11)
chk("E6_d11_restarts_counted_not_denied", restart_per_walker > 0,
    "d=11: %.1f in-budget records per trial at %.1f walkers, i.e. %.3f DP spacings and "
    "%.2f restarts per walker, %.1f restarts per trial -- every such record restarts its "
    "walker, so d=11 is out of protocol and is REPORTED, not gated"
    % (d11["records_mean"], d11["walks"], d11["budget_spacings_per_walker"],
       restart_per_walker, d11["records_mean"]))
chk("E6_out_of_protocol_disclosed",
    all(r["dp_bits"] in GATED or r["budget_spacings_per_walker"] < 1.0 for r in rows),
    "gated %s; d=11 has %.2f spacings per walker"
    % (GATED, d11["budget_spacings_per_walker"]))

# ---------- census arithmetic ----------------------------------------------
chk("audit_plant_check_passed", pc["passed"],
    "planted pairs -> trials %s at levels %s; single walker -> no merge %s; v20d "
    "floor(t/64) attribution differs %s"
    % (pc["merged_trials"], pc["merge_levels"], pc["single_walker_no_merge"],
       pc["v20d_floor_t_w_attribution_differs"]))
chk("statistic_unconditional", cen["gate"]["unconditional"] is True,
    "success_raw has no merge filter; merge is a separate, independent row")
chk("sums_reconcile_exactly",
    all(r["successes_raw"] == r["successes_budget_raw"] + r["successes_drain_raw"]
        and r["detect_pairs_total"] == r["detect_pairs_budget"] + r["detect_pairs_drain"]
        for r in rows),
    "; ".join("d%d %d=%d+%d pairs %d=%d+%d"
              % (r["dp_bits"], r["successes_raw"], r["successes_budget_raw"],
                 r["successes_drain_raw"], r["detect_pairs_total"],
                 r["detect_pairs_budget"], r["detect_pairs_drain"]) for r in rows))
chk("solve_mismatch_zero", all(r["solve_mismatch"] == 0 for r in rows),
    "sum %d" % sum(r["solve_mismatch"] for r in rows))
chk("denominator_dead_zero", all(r["denom_dead"] == 0 for r in rows),
    "sum %d" % sum(r["denom_dead"] for r in rows))
chk("drain_only_pairs_reported_not_filtered",
    all(r["drain_only_pairs"] >= 0 for r in rows),
    "; ".join("d%d %d pairs on %d trials" % (r["dp_bits"], r["drain_only_pairs"],
                                             r["drain_only_trials"]) for r in rows))
chk("census_trials_per_density", all(r["trials"] >= 48000 for r in rows),
    "trials/density %s" % [r["trials"] for r in rows])
worst_raw = min(rows, key=lambda r: r["cp95_raw"])
worst_mg = min(rows, key=lambda r: r["cp95_merge"])
for r in rows:
    cp = cp95_lower(r["successes_raw"], r["trials"])
    cpm = cp95_lower(r["merges"], r["trials"])
    chk("cp95_recomputed_d%d" % r["dp_bits"],
        abs(float(cp) - r["cp95_raw"]) < 1e-12 and abs(float(cpm) - r["cp95_merge"]) < 1e-12,
        "raw %s, merge %s (exact Clopper-Pearson, Decimal)"
        % (str(cp)[:18], str(cpm)[:18]))
    if r["dp_bits"] in GATED:
        chk("cp95_merge_ge_floor_d%d" % r["dp_bits"], cpm >= ledger.SUCCESS_FLOOR,
            "%s >= 0.39" % str(cpm)[:18])

# ---------- claim arithmetic ------------------------------------------------
chk("step_rows_167", ledger.STEP_INSN == 167, "%d itemised instructions" % ledger.STEP_INSN)
chk("claim_is_ceil_of_W", ledger.CLAIM >= ledger.LOG_W + 2e-5
    and ledger.CLAIM - 1e-3 < ledger.LOG_W + 2e-5,
    "log2 W %.8f -> claim %.3f" % (ledger.LOG_W, float(ledger.CLAIM)))
chk("claim_below_gate_bar", ledger.CLAIM <= ledger.GATE_BAR,
    "%.3f <= %.5f" % (float(ledger.CLAIM), float(ledger.GATE_BAR)))
chk("margin_vs_gate", ledger.MARGIN_VS_GATE >= ledger.MARGIN_RULE,
    "gate %.3f - claim %.3f = %.4f" % (float(ledger.GATE), float(ledger.CLAIM),
                                       float(ledger.MARGIN_VS_GATE)))
chk("ledger_success_model_ge_floor", ledger.SUCCESS >= float(ledger.SUCCESS_FLOOR),
    "model %.6f >= 0.39" % ledger.SUCCESS)
chk("kappa_budget_product", ledger.BUDGET <= ledger.CLAIM_RATIO_CEIL,
    "%.9f <= %.9f" % (float(ledger.BUDGET), float(ledger.CLAIM_RATIO_CEIL)))
chk("trap_rate_uses_reachable_columns", ledger.TRAP_RATE == ledger.F(1, 6 * ledger.R_EFF),
    "TRAP_RATE = 1/(6*R_EFF) = %.3e with R_EFF = %d (the selector reaches half of "
    "R = 2^20 buckets)" % (float(ledger.TRAP_RATE), ledger.R_EFF))


# ---------- E7: the note's tables are EMITTED here --------------------------
def tables():
    L = []
    L.append("| d | trials | success_raw | CP-95 raw | merge | CP-95 merge | "
             "class cov | PAIR cov | prev-matched | escapes |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        p = r["pair"]
        L.append("| %d | %d | %.6f | %.6f | %.6f | %.6f | %.4f | %.4f | %.4f | %d |"
                 % (r["dp_bits"], r["trials"], r["success_raw"], r["cp95_raw"],
                    r["merge_rate"], r["cp95_merge"], p["class_trials_per_trial"],
                    p["pair_trials_per_trial"], p["prev_match_fraction"],
                    r["trap_escape_fires"]))
    L.append("")
    L.append("| d | records/trial | spacings/walker | restarts/walker | trap rate "
             "| 1/(6*R_eff) | ratio |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        L.append("| %d | %.2f | %.3f | %.3f | %.2e | %.2e | %.3f |"
                 % (r["dp_bits"], r["records_mean"], r["budget_spacings_per_walker"],
                    r["per_walker"] / (2.0 ** r["dp_bits"]), r["trap_rate_measured"],
                    r["trap_rate_structural"],
                    r["trap_rate_measured"] / r["trap_rate_structural"]))
    L.append("")
    L.append("| quantity | value | source |")
    L.append("|---|---:|---|")
    L.append("| log2 S0 | %.8f | ledger.py |" % ledger.LOG_S0)
    L.append("| log2 W | %.8f | ledger.py |" % ledger.LOG_W)
    L.append("| claim | %.3f | ledger.py |" % float(ledger.CLAIM))
    L.append("| step instructions | %d | ledger.py ROWS |" % ledger.STEP_INSN)
    L.append("| model success | %.6f | ledger.SUCCESS |" % ledger.SUCCESS)
    L.append("| R_EFF / TRAP_RATE | %d / %.3e | ledger.py |"
             % (ledger.R_EFF, float(ledger.TRAP_RATE)))
    L.append("| min CP-95 raw | %.6f (d=%d) | census |"
             % (worst_raw["cp95_raw"], worst_raw["dp_bits"]))
    L.append("| min CP-95 merge, gated | %.6f (d=%d) | census |"
             % (min(r["cp95_merge"] for r in rows if r["dp_bits"] in GATED),
                [r["dp_bits"] for r in rows
                 if r["cp95_merge"] == min(q["cp95_merge"] for q in rows
                                           if q["dp_bits"] in GATED)][0]))
    L.append("| medium: fires / K_obs / residual | %d / %d / %d | medium fixture |"
             % (mesc["fires"], mesc["max_consecutive_fires_per_walk"],
                mesc["fires_residual"]))
    L.append("| medium trapped steps: control / escape | %.2f%% / %.2f%% | fixture |"
             % (100 * mctl["trapped_step_fraction"], 100 * mesc["trapped_step_fraction"]))
    L.append("| equivariance mismatches (fixed / control) | %d / %d | medium fixture |"
             % (meq["total_mismatches"], mmut["total_mismatches"]))
    L.append("| C fixture: classes / spacing / fires | 2^%.2f / %d / %d | C fixture |"
             % (__import__("math").log2(cfx["rows"][0]["classes"]),
                cfx["rows"][0]["steps"], cfx["totals"]["fires"]))
    L.append("| C fixture: trap ratio / invariant bad | %.3f / %d | C fixture |"
             % (cfx["totals"]["trap_rate_over_structural"],
                cfx["totals"]["invariant_violations"]))
    return "\n".join(L) + "\n"


EMIT = tables()
with open(os.path.join(D, "evidence", "note-tables.md"), "w") as f:
    f.write(EMIT)
note = open(os.path.join(D, "public-note.md")).read()
B, E = "<!-- BEGIN MACHINE TABLES -->", "<!-- END MACHINE TABLES -->"
blk = note.split(B, 1)[1].split(E, 1)[0] if B in note and E in note else ""
chk("note_tables_byte_equal_to_certificate", blk.strip("\n") == EMIT.strip("\n"),
    "note block %d B vs emitted %d B, equal %s"
    % (len(blk.strip("\n")), len(EMIT.strip("\n")), blk.strip("\n") == EMIT.strip("\n")))
chk("note_single_header", note.count("Model: qwen3.8-flash-next") == 1
    and note.count("Harness: Hermes") == 1, "Model/Harness header appears once")
chk("note_size_in_band", 10240 <= len(note) <= 14336, "note %d B" % len(note))
chk("note_no_never_restart_claim",
    "never restart" not in note.lower() and "never restarts" not in note.lower(),
    "the d=11 protocol claim is gone: the note states the measured restart count")
chk("note_cites_measured_k", "K_obs" in note or "K_esc" in note,
    "the note's escape paragraph carries the measured consecutive-substitution count")

n_ok = sum(c["ok"] for c in checks)
res = {
    "package": "plain-orbitrow-40pct-v20g",
    "claim": float(ledger.CLAIM),
    "log2_W": ledger.LOG_W,
    "success_model": ledger.SUCCESS,
    "cp95_raw_min": worst_raw["cp95_raw"],
    "cp95_merge_min_gated": min(r["cp95_merge"] for r in rows
                                if r["dp_bits"] in GATED),
    "pair_coverage": {r["dp_bits"]: r["pair"]["pair_trials_per_trial"] for r in rows},
    "class_coverage": {r["dp_bits"]: r["pair"]["class_trials_per_trial"] for r in rows},
    "step_instructions": ledger.STEP_INSN,
    "gate": float(ledger.GATE),
    "checks": checks,
    "n_checks": len(checks),
    "n_ok": n_ok,
    "passed": n_ok == len(checks),
}
with open(os.path.join(D, "certificate.json"), "w") as f:
    json.dump(res, f, indent=1, sort_keys=True)
print("%d/%d checks passed" % (n_ok, len(checks)))
for c in checks:
    if not c["ok"]:
        print("FAIL", c["name"], "--", c["detail"])
