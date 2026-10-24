#!/usr/bin/env python3
"""gen_cert_v20d.py -- machine certificate for plain-orbitrow-40pct-v20d.

Every number is recomputed here from evidence/prod-sweep.json, the low-n
disclosure evidence/prod-sweep-8000.json and the ledger machinery.  Nothing is
hand-typed except the package name.  Writes certificate.json; exits nonzero if
any check is false.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
from prod_sweep import cp95_lower          # noqa: E402  (same function the census used)
import ledger                               # noqa: E402

checks = []


def chk(name, ok, detail):
    checks.append({"name": name, "ok": bool(ok), "detail": detail})


def ev(name):
    with open(os.path.join(D, "evidence", name)) as f:
        return json.load(f)


cen = ev("prod-sweep.json")
low = ev("prod-sweep-8000.json")
rows = cen["rows"]
dens = sorted(r["dp_bits"] for r in rows)

# ---------- 1. census shape ----------
chk("census_densities_5_7_9_11", dens == [5, 7, 9, 11], "densities %s" % dens)
chk("census_48000_trials_each",
    all(r["trials"] == 48000 for r in rows),
    "trials/density %s, total %d" % ([r["trials"] for r in rows],
                                     sum(r["trials"] for r in rows)))

# ---------- 2. CP-95 recomputed exactly, floor at every density ----------
worst = None
for r in rows:
    cp = cp95_lower(r["successes"], r["trials"])
    r["_cp"] = cp
    if worst is None or cp < worst[1]:
        worst = (r["dp_bits"], cp)
    chk("cp95_recomputed_d%d" % r["dp_bits"],
        abs(float(cp) - r["cp95_lower"]) < 1e-12,
        "Decimal %s vs stored %.12f" % (str(cp)[:20], r["cp95_lower"]))
    chk("cp95_ge_0p39_d%d" % r["dp_bits"], cp >= ledger.SUCCESS_FLOOR,
        "CP-95 lower %s >= 0.39 (k=%d n=%d)" % (str(cp)[:18], r["successes"], r["trials"]))
chk("cp95_floor_all_densities", worst is not None and worst[1] >= ledger.SUCCESS_FLOOR,
    "min CP-95 %s at d=%d" % (str(worst[1])[:18], worst[0]))

# ---------- 3. success <= merge, recovery verified, no algebra residue ----------
chk("success_le_merge_all", all(r["successes"] <= r["merges"] for r in rows),
    "; ".join("d%d %d<=%d" % (r["dp_bits"], r["successes"], r["merges"]) for r in rows))
chk("success_without_merge_zero", all(r["success_without_merge"] == 0 for r in rows),
    "sum %d" % sum(r["success_without_merge"] for r in rows))
chk("solve_mismatch_zero", all(r["solve_mismatch"] == 0 for r in rows),
    "sum %d" % sum(r["solve_mismatch"] for r in rows))
chk("denominator_dead_zero", all(r["denom_dead"] == 0 for r in rows),
    "sum %d" % sum(r["denom_dead"] for r in rows))
chk("recoveries_verified_cover_successes",
    all(r["detect_verified"] >= r["successes"] for r in rows),
    "; ".join("d%d %d>=%d" % (r["dp_bits"], r["detect_verified"], r["successes"])
              for r in rows))

# ---------- 4. records match the priced pool ----------
chk("records_match_budget_over_dp",
    all(abs(r["records_budget_mean"] - r["records_expected_mean"])
        <= 0.03 * r["records_expected_mean"] for r in rows),
    "; ".join("d%d %.1f vs %.1f" % (r["dp_bits"], r["records_budget_mean"],
                                    r["records_expected_mean"]) for r in rows))
chk("drain_reports_le_one_per_walker",
    all(r["records_drain_mean"] <= r["walks"] for r in rows),
    "; ".join("d%d %.1f <= %d" % (r["dp_bits"], r["records_drain_mean"], r["walks"])
              for r in rows))

# ---------- 5. drain is a flush, not a search ----------
tot_ex = sum(r["drain_only_excluded"] for r in rows)
tot_mg = sum(r["merges"] for r in rows)
chk("drain_only_collisions_excluded_and_small", tot_ex <= 0.05 * tot_mg,
    "%d excluded duplicate pairs vs %d merged trials (%.3f%%)"
    % (tot_ex, tot_mg, 100.0 * tot_ex / tot_mg))
chk("drain_delivers_reports",
    all(r["records_drain_mean"] > 0 for r in rows),
    "; ".join("d%d %.1f" % (r["dp_bits"], r["records_drain_mean"]) for r in rows))
chk("in_budget_success_reported_separately",
    all(0 <= r["success_budget"] <= r["success"] for r in rows),
    "; ".join("d%d %.4f b + %.4f d" % (r["dp_bits"], r["success_budget"],
                                       r["success_drain"]) for r in rows))

# ---------- 6. low-n disclosure is present and is NOT the gate ----------
lowrows = low["rows"]
lowmin = min(r["cp95_lower"] for r in lowrows)
chk("low_n_8000_disclosed", lowmin < 0.39,
    "8000/density census min CP-95 %.4f, disclosed as the precision shortfall" % lowmin)
chk("low_n_8000_same_arm", all(r["trials"] == 8000 for r in lowrows)
    and sorted(r["dp_bits"] for r in lowrows) == dens,
    "same four densities, 8000 trials each")

# ---------- 7. claim arithmetic, exact ----------
chk("step_rows_161", ledger.STEP_INSN == 161, "%d itemised instructions" % ledger.STEP_INSN)
chk("claim_is_ceil_of_W_plus_2e5",
    ledger.CLAIM >= ledger.LOG_W + 2e-5 and ledger.CLAIM - 1e-3 < ledger.LOG_W + 2e-5,
    "log2 W %.8f -> claim %.3f" % (ledger.LOG_W, float(ledger.CLAIM)))
chk("claim_126p262", abs(float(ledger.CLAIM) - 126.262) < 1e-12, "claim %.3f" % float(ledger.CLAIM))
chk("claim_below_ship_line", ledger.CLAIM <= ledger.SHIP_LINE,
    "%.3f <= %.3f" % (float(ledger.CLAIM), float(ledger.SHIP_LINE)))
chk("margin_vs_gate_at_least_0p008", ledger.MARGIN_VS_GATE >= ledger.MARGIN_RULE,
    "gate %.3f - claim %.3f = %.4f >= 0.008" % (float(ledger.GATE),
                                                float(ledger.CLAIM),
                                                float(ledger.MARGIN_VS_GATE)))
chk("ledger_success_model_ge_floor", ledger.SUCCESS >= float(ledger.SUCCESS_FLOOR),
    "model %.6f >= 0.39" % ledger.SUCCESS)
chk("kappa_budget_product", ledger.BUDGET <= ledger.CLAIM_RATIO_CEIL,
    "%.9f <= %.9f" % (float(ledger.BUDGET), float(ledger.CLAIM_RATIO_CEIL)))

res = {
    "package": "plain-orbitrow-40pct-v20d",
    "claim": float(ledger.CLAIM),
    "log2_W": ledger.LOG_W,
    "success_model": ledger.SUCCESS,
    "success_census_min": min(r["success"] for r in rows),
    "cp95_lower": float(worst[1]),
    "checks": checks,
    "passed": sum(c["ok"] for c in checks),
    "total": len(checks),
    "verdict": "PASS" if all(c["ok"] for c in checks) else "FAIL",
}
with open(os.path.join(D, "certificate.json"), "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps({"verdict": res["verdict"], "pass": res["passed"], "of": res["total"],
                  "cp95_lower": res["cp95_lower"],
                  "failed": [c["name"] for c in checks if not c["ok"]]}, indent=1))
sys.exit(0 if res["verdict"] == "PASS" else 1)
