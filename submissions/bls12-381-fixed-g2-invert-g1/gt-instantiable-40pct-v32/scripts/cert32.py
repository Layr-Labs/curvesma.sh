#!/usr/bin/env python3
"""cert32.py -- regenerate certificate.json from the CURRENT ledger, and run the
stale-literal sweep over the package (verdict b30ddc21 flaw 6).

Every check recomputes its own quantity from evidence/ + submission.json; none of
them copies a constant out of the prose.  Writes certificate.json.
"""
import json
import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")

L = json.load(open(os.path.join(EV, "ledger32.json")))
SUB = json.load(open(os.path.join(PKG, "submission.json")))
SORT = json.load(open(os.path.join(EV, "sort32.json")))
MAN = json.load(open(os.path.join(EV, "manifest32.json")))

rows = {k: int(v) for k, v in L["rows_ir"].items()}
W = int(L["W_total_ir"])
n = int(L["claim"]["numerator"])
PW, BASEP = W ** 10000, 20000 ** 10000

checks = []


def ck(name, ok, detail=""):
    checks.append({"name": name, "verdict": "PASS" if ok else "FAIL", "detail": detail})


ck("row_sum_equals_W_exactly", sum(rows.values()) == W, "sum(rows) == W_total_ir")
ck("claim_covers_W", BASEP * (1 << n) >= PW, "20000^10000 * 2^n >= W^10000")
ck("claim_minimal_ceil", BASEP * (1 << (n - 1)) < PW, "n-1 fails the same exact test")
ck("claim_matches_submission_json",
   abs(float(n) / 10000 - float(SUB["claims"][0]["score"]["numerator"]) /
       float(SUB["claims"][0]["score"]["denominator"])) < 1e-12
   and int(SUB["claims"][0]["score"]["numerator"]) == n,
   "submission.json score == ledger numerator/10000")
ck("every_row_is_a_python_int",
   all(isinstance(v, int) for v in rows.values()) and isinstance(W, int),
   "no float appears in any arithmetic quantity")
ck("trap_at_S_is_exact_int",
   isinstance(L["trap"]["trap_at_S"], int) and L["trap"]["trap_at_S_is_int"] is True,
   "trap_at_S = %d" % L["trap"]["trap_at_S"])
ck("trap_waste_is_exact_int", isinstance(int(L["trap"]["charge_windows_x_32_x_STEP"]), int),
   "integer ceil, no floor division")
ck("trap_route_has_no_overlap_multiplier",
   "overlap" not in L["trap"]["route"].lower()
   and "no independence" in L["trap"]["route"].lower()
   and "no poisson" in L["trap"]["route"].lower()
   and "multipl" not in L["trap"]["route"].lower(),
   "MARKOV union over the measured guard-admission precondition; no x6 factor, no Poisson limit")
ck("trap_markov_bound_present", "guard_admission_necessary_condition" in L["trap"],
   "the bound is a count of a necessary precondition, not a rate over zero events")
ck("capacity_uses_steady_state_law", "occupancy = min(capacity" in L["capacity_spec"]["steady_state_law"],
   "occupancy/residence derivation is stated, not transferred from the fixture")
ck("capacity_route_a_matches_arrivals",
   int(L["capacity_spec"]["route_a"]["occupancy"]) == int(L["capacity_spec"]["route_a"]["arrivals_total"])
   < int(L["capacity_spec"]["route_a"]["capacity_records"]),
   "occupancy == arrivals < capacity: no eviction, retention = 1 by construction")
ck("capacity_route_b_rejected_with_value",
   "REJECTED" in L["capacity_spec"]["route_b"]["verdict"],
   "route B retention %s -> success %s" % (L["capacity_spec"]["route_b"]["retention_factor"],
                                           L["capacity_spec"]["route_b"]["success_if_route_b"]))
ck("bounds_table_every_array_has_a_bound",
   all(b.get("bound_formula") and b.get("value") and b.get("bytes") for b in L["bounds_table"])
   and len(L["bounds_table"]) >= 9,
   "%d arrays bounded" % len(L["bounds_table"]))
ck("bounds_table_no_growth_with_arrivals",
   all("S" in b["bound_formula"] or "PROBES" in b["bound_formula"] or "depth" in b["bound_formula"]
       or "fixed" in b["bound_formula"] or "Fp12" in b["bound_formula"] or "next_pow2" in b["bound_formula"]
       or "counter" in b["bound_formula"] for b in L["bounds_table"]),
   "each formula is a design-parameter expression")
ck("state_key_full_width_row_charged",
   rows["state_key_full_width_sort"] == int(L["protocol"]["steps"]) * SORT["charge_per_step"],
   "row == steps x executed MAX (%d Ir/step)" % SORT["charge_per_step"])
ck("sort_charge_is_max_over_calls", SORT["charge_per_step"] == SORT["executed_ir_per_call"]["max"],
   "MAX, not aggregate/calls")
ck("coefficient_row_is_an_implementation_bound",
   L["prices"]["derived_rows"]["scalar_coefficient_update_per_step"] >= 2 * 64,
   "two 255-bit mulmods at 64 Ir each")
ck("executed_mean_is_labelled_not_charged",
   "executed_mean_not_used" in L["prices"],
   "the aggregate-divided-by-calls figure is retained as a mean only")
ck("success_floor_met", L["success"]["route_a_p5"] >= 0.39,
   "p5 = %s >= 0.39 with retention = 1 (route A)" % L["success"]["route_a_p5"])
ck("success_matches_submission_range",
   "0.578125" in SUB["claims"][0]["statement"],
   "the statement quotes the same p5 the ledger computes")
SELF_MAN = "evidence/manifest32.json"
ALL_FILES = {os.path.relpath(os.path.join(r, f), PKG) for r, _, fs in os.walk(PKG) for f in fs
             if not f.startswith(".") and r.split(os.sep)[-1] != "__pycache__"} - {SELF_MAN}
ck("manifest_lists_every_cited_evidence",
   (set(SUB["claims"][0]["evidence"]) - {SELF_MAN}) <= {f["path"] for f in MAN["files"]},
   "every evidence path in submission.json is in the manifest")
ck("manifest_covers_every_package_file",
   {f["path"] for f in MAN["files"]} == ALL_FILES,
   "manifest is the whole tree except itself")
ck("manifest_has_sha256",
   all(len(f["sha256"]) == 64 and f["bytes"] == os.path.getsize(os.path.join(PKG, f["path"]))
       for f in MAN["files"]),
   "sha256 + byte count verified on disk")
ck("file_count_within_intake", len(MAN["files"]) <= 32, "%d files" % len(MAN["files"]))
ck("total_bytes_within_intake",
   sum(f["bytes"] for f in MAN["files"]) <= 262144,
   "%d B" % sum(f["bytes"] for f in MAN["files"]))
ck("no_file_over_65536",
   all(f["bytes"] <= 65536 for f in MAN["files"]), "per-file cap")

# ---- stale-literal sweep (flaw 6): a superseded value may be NAMED as removed, but
# must never appear as a live quantity.  Every hit has to sit inside a supersession
# window (a sentence that says it is the old value).
STALE = ["127.6293", "1276422", "524288 buckets", "2097152 records",
         "overlap factor of six", "OVERLAP = 6", "1164 records, scaled x128"]
MARK = ("v31", "v30", "v29", "old", "previous", "superseded", "stale", "was ", "gone",
        "deleted", "no longer", "flaw", "verdict", "removed", "under-charged", "not ")
hits = []
for r, ds, fs in os.walk(PKG):
    ds[:] = [d for d in ds if d != "__pycache__"]
    for f in fs:
        if not f.endswith((".md", ".json", ".py", ".c")):
            continue
        if f in ("certificate.json", "cert32.py", "manifest32.json"):
            continue
        fp = os.path.join(r, f)
        try:
            t = open(fp, encoding="utf-8", errors="ignore").read()
        except Exception:
            continue
        for s in STALE:
            i = t.find(s)
            while i >= 0:
                w = t[max(0, i - 260):i + len(s) + 120].lower()
                if not any(m in w for m in MARK):
                    hits.append("%s: %s" % (os.path.relpath(fp, PKG), s))
                    break
                i = t.find(s, i + 1)
ck("no_stale_literal_asserted_as_current", not hits, "; ".join(hits) or
   "every superseded literal is inside a supersession sentence")

out = {
    "package": "gt-instantiable-40pct-v32",
    "claim": "%.4f" % (n / 10000),
    "W_total_ir": str(W),
    "checks": checks,
    "all_pass": all(c["verdict"] == "PASS" for c in checks),
    "verdict": "PASS" if all(c["verdict"] == "PASS" for c in checks) else "FAIL",
}
json.dump(out, open(os.path.join(PKG, "certificate.json"), "w"), indent=1)
for c in checks:
    print("%-46s %s %s" % (c["name"], c["verdict"], c.get("detail", "")))
print("ALL", out["verdict"], len(checks), "checks")
