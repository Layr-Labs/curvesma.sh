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
from prod_sweep import cp95_lower  # noqa: E402


def ev(name):
    with open(os.path.join(D, "evidence", name)) as f:
        return json.load(f)


def values():
    """Every quantity the note is allowed to print, re-derived from the artifacts."""
    cen = ev("prod-sweep.json")
    med = ev("medium-fixture.json")
    cfx = ev("c-fixture.json")
    tw = ev("twist-index.json")
    cy = ev("cycle3.json")
    kc = ev("escape-kcap.json")
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
                 ("n_cells", ledger.N_CELLS), ("K_ESC", ledger.K_ESC)):
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
    add("cy.structural", cy["escape_armed"]["three_cycle_structural"])
    add("cy.steps", cy["escape_armed"]["lane_steps"])
    add("cy.visits", cy["escape_armed"]["three_cycle_visits"])
    add("cy.esc2", cy["escape_armed"]["two_cycle_visits"])
    add("cy.ctl2", cy["control_no_escape"]["two_cycle_visits"])
    add("cy.L2", cy["control_no_escape"]["settled_cycle_lengths"].get("2", 0))
    add("cy.L3", cy["control_no_escape"]["settled_cycle_lengths"].get("3", 0))
    add("kc.kcap", kc["kcap"])
    add("kc.max", kc["kcap_max_seen"])
    add("kc.trials", kc["trials_total"])
    add("kc.fires", kc["fires_total"])
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
    cy = ev("cycle3.json")
    kc = ev("escape-kcap.json")
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
        "kcap_charged_equals_enforced": ledger.K_ESC == ledger.qwalk.K_ESC == kc["kcap"],
        "kcap_never_tripped": kc["kcap_max_seen"] <= kc["kcap"] and not kc["trips"],
        "twist_selector_exact": tw["agreement"] == tw["N"] and tw["disagreement"] == 0,
        "twist_mutation_fails": tw["mutation"]["disagreements"] > 0,
        "three_cycle_structural_zero": cy["escape_armed"]["three_cycle_structural"] == 0,
        "control_shows_long_cycles": any(
            int(k) >= 3 for k in cy["control_no_escape"]["settled_cycle_lengths"]
            if int(k) <= 64),
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
