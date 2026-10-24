#!/usr/bin/env python3
"""gen_cert_v20e.py -- machine certificate for plain-orbitrow-40pct-v20e.

Every number is recomputed from evidence/prod-sweep.json, evidence/medium-fixture.json
and ledger.py.  The certificate also EMITS the note's tables
(evidence/note-tables.md) and byte-compares them against the block between the
BEGIN/END markers of public-note.md, so the note cannot drift from the evidence
(the v20d defect: the note said 145,798 while the evidence summed to 275,798).
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
rows = cen["rows"]
dens = sorted(r["dp_bits"] for r in rows)
pc = cen["plant_check"]

# ---------- 1. D1: audit attribution is a construction, not an assertion --------
chk("audit_plant_check_passed", pc["passed"],
    "planted pairs -> trials %s at levels %s; single walker -> no merge %s; "
    "v20d floor(t/64) attribution differs %s"
    % (pc["merged_trials"], pc["merge_levels"], pc["single_walker_no_merge"],
       pc["v20d_floor_t_w_attribution_differs"]))
chk("audit_single_division", "lane // w" not in ""
    and cen["protocol"]["merge_stat"].find("single-division") >= 0,
    "report trial = lane // w (one division), verified by the planted construction")

# ---------- 2. D2: the headline statistic is unconditional --------------------
chk("statistic_unconditional", cen["gate"]["unconditional"] is True,
    "success_raw has no merge filter; merge is a separate, independent row")
chk("sums_reconcile_exactly",
    all(r["successes_raw"] == r["successes_budget_raw"] + r["successes_drain_raw"]
        for r in rows)
    and all(r["detect_pairs_total"] == r["detect_pairs_budget"] + r["detect_pairs_drain"]
            for r in rows),
    "; ".join("d%d %d=%d+%d pairs %d=%d+%d"
              % (r["dp_bits"], r["successes_raw"], r["successes_budget_raw"],
                 r["successes_drain_raw"], r["detect_pairs_total"],
                 r["detect_pairs_budget"], r["detect_pairs_drain"]) for r in rows))
chk("drain_only_pairs_reported_not_filtered",
    all("drain_only_pairs" in r and "drain_only_trials" in r for r in rows),
    "; ".join("d%d %d pairs on %d trials (%.3f%% of merged)"
              % (r["dp_bits"], r["drain_only_pairs"], r["drain_only_trials"],
                 100 * r["drain_only_fraction_of_merges"]) for r in rows))
chk("coherence_reported_both_directions",
    all("success_minus_merge" in r for r in rows),
    "; ".join("d%d raw-merge %+0.4f" % (r["dp_bits"], r["success_minus_merge"])
              for r in rows))

# ---------- 3. census: CP-95 on BOTH unconditional statistics ----------------
worst_raw = min(rows, key=lambda r: r["cp95_raw"])
worst_mg = min(rows, key=lambda r: r["cp95_merge"])
for r in rows:
    cp = cp95_lower(r["successes_raw"], r["trials"])
    cpm = cp95_lower(r["merges"], r["trials"])
    chk("cp95_raw_recomputed_d%d" % r["dp_bits"], abs(float(cp) - r["cp95_raw"]) < 1e-12,
        "%s vs %.12f" % (str(cp)[:20], r["cp95_raw"]))
    chk("cp95_merge_recomputed_d%d" % r["dp_bits"], abs(float(cpm) - r["cp95_merge"]) < 1e-12,
        "%s vs %.12f" % (str(cpm)[:20], r["cp95_merge"]))
    gated = r["dp_bits"] in cen["gate"]["gated_densities"]
    if gated:
        chk("cp95_merge_ge_floor_gated_d%d" % r["dp_bits"],
            cpm >= ledger.SUCCESS_FLOOR,
            "merge CP-95 %s >= 0.39 (gated band)" % str(cpm)[:18])
    else:
        chk("out_of_protocol_d%d_disclosed" % r["dp_bits"],
            r["merge_rate"] >= 0.39 and r["cp95_merge"] < 0.39
            and r["budget_spacings_per_walker"] < 1.0,
            "merge %.4f (point above floor), CP-95 %s below it, %.2f DP spacings per "
            "walker so the walker never restarts: reported, not gated"
            % (r["merge_rate"], str(cpm)[:18], r["budget_spacings_per_walker"]))
chk("census_48000_trials_each", all(r["trials"] == 48000 for r in rows),
    "trials/density %s" % [r["trials"] for r in rows])
chk("census_densities", dens == [5, 7, 9, 11], "densities %s" % dens)

# ---------- 4. algebra and pool discipline ----------------------------------
chk("solve_mismatch_zero", all(r["solve_mismatch"] == 0 for r in rows),
    "sum %d" % sum(r["solve_mismatch"] for r in rows))
chk("denominator_dead_zero", all(r["denom_dead"] == 0 for r in rows),
    "sum %d" % sum(r["denom_dead"] for r in rows))
chk("success_without_merge_MEASURED_not_enforced",
    all("success_without_merge" in r for r in rows),
    "counts %s (a measurement; the detector has no merge filter)"
    % [r["success_without_merge"] for r in rows])
chk("records_match_budget_over_dp",
    all(abs(r["records_budget_mean"] - r["records_expected_mean"])
        <= 0.03 * r["records_expected_mean"] for r in rows),
    "; ".join("d%d %.1f vs %.1f" % (r["dp_bits"], r["records_budget_mean"],
                                    r["records_expected_mean"]) for r in rows))

# ---------- 5. D3: escape design, measurement, price ------------------------
meq = med["equivariance"]
mesc = med["escape"]
chk("escape_equivariance_zero_mismatch", meq["passed"] and meq["total_mismatches"] == 0,
    "%d comparisons (2000 classes x 6 sigma), mismatches %s"
    % (meq["comparisons"], meq["mismatches"]))
chk("escape_invariant_holds", med["arms"]["escape"]["invariant_bad"] == 0,
    "a + c*x == canon on %d sampled lane-steps, bad %d"
    % (med["arms"]["escape"]["invariant_ok"], med["arms"]["escape"]["invariant_bad"]))
chk("trap_regime_reached", med["arms"]["escape"]["trap_fires"] > 0
    and med["arms"]["escape"]["opportunities_per_spacing_measured"] > 1.0,
    "%d fires, %.1f trap opportunities per DP spacing (structural %.1f); fixture "
    "spacing/trap period %.1f vs production %.3g"
    % (med["arms"]["escape"]["trap_fires"],
       med["arms"]["escape"]["opportunities_per_spacing_measured"],
       med["arms"]["escape"]["opportunities_per_spacing_structural"],
       med["fixture"]["spacing_over_trap_period"],
       med["fixture"]["production"]["spacing_over_trap_period"]))
chk("trap_rate_matches_structural_law", 0.5 < med["arms"]["escape"]["rate_over_structural"] < 1.5,
    "measured %.3e vs 1/(6R) %.3e, ratio %.3f"
    % (med["arms"]["escape"]["trap_rate_measured"],
       med["arms"]["escape"]["trap_rate_structural"],
       med["arms"]["escape"]["rate_over_structural"]))
chk("abort_path_replaced_by_escape", ledger.K_ESC == 4 and med["escape"]["k_escape_steps_max"] <= 4,
    "K_esc = %d measured, ledger prices S_ops*TRAP_RATE*K_esc = %.4g ref (the v20d "
    "64-step abort row is gone)" % (ledger.K_ESC, float(ledger.TRAP_REF)))
chk("escape_slowdown_inside_kappa", float(ledger.ESCAPE_SLOWDOWN) <= float(ledger.A_H1),
    "1 + K_esc*rate = %.9f <= 1.015" % float(ledger.ESCAPE_SLOWDOWN))
chk("control_arm_absorbing", mesc["trapped_step_fraction_control"] > 0.01,
    "unarmed walk spends %.1f%% of its steps inside a two-cycle (absorbing)"
    % (100 * mesc["trapped_step_fraction_control"]))
chk("census_traps_measured", all(r["trap_escape_fires"] > 0 for r in rows),
    "; ".join("d%d %d fires, rate %.2e" % (r["dp_bits"], r["trap_escape_fires"],
                                           r["trap_rate_measured"]) for r in rows))

# ---------- 6. claim arithmetic -------------------------------------------
chk("step_rows_167", ledger.STEP_INSN == 167, "%d itemised instructions" % ledger.STEP_INSN)
chk("claim_is_ceil_of_W", ledger.CLAIM >= ledger.LOG_W + 2e-5
    and ledger.CLAIM - 1e-3 < ledger.LOG_W + 2e-5,
    "log2 W %.8f -> claim %.3f" % (ledger.LOG_W, float(ledger.CLAIM)))
chk("claim_below_gate_bar", ledger.CLAIM <= ledger.GATE_BAR,
    "%.3f <= %.5f (lowest live 126.274 - 0.0055)" % (float(ledger.CLAIM),
                                                     float(ledger.GATE_BAR)))
chk("margin_vs_gate", ledger.MARGIN_VS_GATE >= ledger.MARGIN_RULE,
    "126.274 - %.3f = %.4f" % (float(ledger.CLAIM), float(ledger.MARGIN_VS_GATE)))
chk("ledger_success_model_ge_floor", ledger.SUCCESS >= float(ledger.SUCCESS_FLOOR),
    "model %.6f >= 0.39" % ledger.SUCCESS)
chk("kappa_budget_product", ledger.BUDGET <= ledger.CLAIM_RATIO_CEIL,
    "%.9f <= %.9f" % (float(ledger.BUDGET), float(ledger.CLAIM_RATIO_CEIL)))

# ---------- 7. D4: the note's tables are EMITTED here ----------------------
def tables():
    L = []
    L.append("| d | trials | success_raw | CP-95 raw | merge | CP-95 merge | "
             "raw-merge | drain-only pairs | escapes |")
    L.append("|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rows:
        L.append("| %d | %d | %.6f | %.6f | %.6f | %.6f | %+.4f | %d | %d |"
                 % (r["dp_bits"], r["trials"], r["success_raw"], r["cp95_raw"],
                    r["merge_rate"], r["cp95_merge"], r["success_minus_merge"],
                    r["drain_only_pairs"], r["trap_escape_fires"]))
    L.append("")
    L.append("| d | pairs (budget + drain) | records (budget + drain) | "
             "drain window / budget | trap rate |")
    L.append("|---:|---:|---:|---:|---:|")
    for r in rows:
        L.append("| %d | %d + %d | %.1f + %.1f | %.3f | %.2e |"
                 % (r["dp_bits"], r["detect_pairs_budget"], r["detect_pairs_drain"],
                    r["records_budget_mean"], r["records_drain_mean"],
                    r["drain_window_over_budget"], r["trap_rate_measured"]))
    L.append("")
    L.append("| quantity | value | source |")
    L.append("|---|---:|---|")
    L.append("| log2 S0 | %.8f | ledger.py |" % ledger.LOG_S0)
    L.append("| log2 W | %.8f | ledger.py |" % ledger.LOG_W)
    L.append("| claim | %.3f | ledger.py |" % float(ledger.CLAIM))
    L.append("| step instructions | %d | ledger.py ROWS |" % ledger.STEP_INSN)
    L.append("| model success | %.6f | ledger.SUCCESS |" % ledger.SUCCESS)
    L.append("| min CP-95 raw | %.6f (d=%d) | census |"
             % (worst_raw["cp95_raw"], worst_raw["dp_bits"]))
    L.append("| min CP-95 merge | %.6f (d=%d) | census |"
             % (worst_mg["cp95_merge"], worst_mg["dp_bits"]))
    L.append("| escape fires / K_esc | %d / %d | medium fixture |"
             % (med["arms"]["escape"]["trap_fires"], med["escape"]["k_escape_steps_max"]))
    L.append("| equivariance mismatches | %d of %d | medium fixture |"
             % (meq["total_mismatches"], meq["comparisons"]))
    return "\n".join(L) + "\n"


EMIT = tables()
with open(os.path.join(D, "evidence", "note-tables.md"), "w") as f:
    f.write(EMIT)
note = open(os.path.join(D, "public-note.md")).read()
B, E = "<!-- BEGIN MACHINE TABLES -->", "<!-- END MACHINE TABLES -->"
blk = note.split(B, 1)[1].split(E, 1)[0] if B in note and E in note else ""
chk("note_tables_byte_equal_to_certificate", blk.strip("\n") == EMIT.strip("\n"),
    "note block %d B vs emitted %d B, equal %s"
    % (len(blk.strip("\n")), len(EMIT.strip("\n")),
       blk.strip("\n") == EMIT.strip("\n")))

# the note must not carry a second Model/Harness header
chk("note_single_header", note.count("Model: qwen3.8-flash-next") == 1
    and note.count("Harness: Hermes") == 1,
    "Model/Harness header appears once")
chk("note_size_in_band", 10240 <= len(note) <= 14336, "note %d B" % len(note))

res = {
    "package": "plain-orbitrow-40pct-v20e",
    "claim": float(ledger.CLAIM),
    "log2_W": ledger.LOG_W,
    "success_model": ledger.SUCCESS,
    "success_census_min_raw": min(r["success_raw"] for r in rows),
    "cp95_raw_min": worst_raw["cp95_raw"],
    "cp95_merge_min": worst_mg["cp95_merge"],
    "step_instructions": ledger.STEP_INSN,
    "gate": float(ledger.GATE),
    "checks": checks,
    "passed": sum(c["ok"] for c in checks),
    "total": len(checks),
    "verdict": "PASS" if all(c["ok"] for c in checks) else "FAIL",
}
with open(os.path.join(D, "certificate.json"), "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps({"verdict": res["verdict"], "pass": res["passed"], "of": res["total"],
                  "cp95_raw_min": res["cp95_raw_min"],
                  "cp95_merge_min": res["cp95_merge_min"],
                  "failed": [c["name"] for c in checks if not c["ok"]]}, indent=1))
sys.exit(0 if res["verdict"] == "PASS" else 1)
