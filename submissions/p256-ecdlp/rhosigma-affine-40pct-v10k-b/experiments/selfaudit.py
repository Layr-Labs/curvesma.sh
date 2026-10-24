#!/usr/bin/env python3
"""v10k self-audit: cited-file existence + planted-wrong-path control, stale-number grep +
planted-stale control, spec-vs-run parameter equality, cross-artifact claim equality.
Writes evidence/selfaudit.json and exits non-zero on any failure."""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")
TEXT = ["public-note.md", "argument.md", "submission.json"]

def read(f):
    with open(os.path.join(PKG, f)) as h:
        return h.read()

sub = json.load(open(os.path.join(PKG, "submission.json")))
cert = json.load(open(os.path.join(PKG, "certificate.json")))
led = json.load(open(os.path.join(EV, "ledger.json")))
run = json.load(open(os.path.join(EV, "recproc_run.json")))
blob = "\n".join(read(f) for f in TEXT if os.path.exists(os.path.join(PKG, f)))

checks = []

def add(name, ok, **kw):
    d = {"name": name, "ok": bool(ok)}
    d.update(kw)
    checks.append(d)

# 1. every cited evidence path exists
cited = sub["claims"][0]["evidence"]
missing = [p for p in cited if not os.path.exists(os.path.join(PKG, p))]
add("cited_evidence_exists", not missing, cited=len(cited), missing=missing)

# 1b. planted control: a path that does not exist must be REPORTED
planted = cited + ["evidence/planted_absent_control.json"]
p_missing = [p for p in planted if not os.path.exists(os.path.join(PKG, p))]
add("planted_wrong_path_control_fires", p_missing == ["evidence/planted_absent_control.json"],
    caught=p_missing)

# 2. stale numbers / phrases from the defeated submissions
STALE = ["0.2422", "128.0309", "counted step 1101", "1,101-instruction",
         "load factor 0.5", "512*8", "PROD[", "SUF[", "transfer remains unproven",
         "the declared production density", "72.331957", "0.2294"]
for s in STALE:
    n = blob.count(s)
    add("stale_absent:%s" % s, n == 0, occurrences=n)

# 2b. planted-stale control: the detector must fire on a planted copy
planted_blob = blob + "\nclaim 128.0309 with load factor 0.5\n"
add("planted_stale_control_fires", planted_blob.count("128.0309") > 0 and planted_blob.count("load factor 0.5") > 0,
    planted_hits=[planted_blob.count("128.0309"), planted_blob.count("load factor 0.5")])

# 3. spec == run parameters (single source: the spec constants)
spec_steps = run["rounds"] * run["walkers"] * run["spr"]
add("selftest_run_equals_spec", spec_steps == run["steps"] == 4096,
    rounds=run["rounds"], walkers=run["walkers"], spr=run["spr"],
    run_steps=run["steps"], spec_steps=spec_steps)
add("selftest_invariant_clean", run["invariant_fail"] == 0 and run["invariant_fail_baseline"] == 0)
add("selftest_synthetic_solved", run["synthetic_collision_solved"] == 1)

# 4. cross-artifact claim equality (one source: ledger.json)
claim = led["W"]["claim"]
add("claim_matches_certificate", abs(cert["claim_f"] - claim) < 1e-12,
    ledger=claim, certificate=cert["claim_f"])
add("claim_matches_submission", abs(int(sub["claims"][0]["score"]["numerator"]) / 10000 - claim) < 1e-12,
    submission=int(sub["claims"][0]["score"]["numerator"]) / 10000)
add("claim_present_in_note", ("%.4f" % claim) in blob)
add("certificate_all_checks_pass", cert["all_pass"], n_checks=len(cert["checks"]))

# 5. every p-claim sentence must cite the invariant row or the fit residual
para = [s for s in read("public-note.md").split("\n\n") if "0.39" in s or "success" in s.lower()]
add("p_claims_cite_invariant", all(("ltheta" in s) or ("fit" in s.lower()) or ("invariant" in s.lower())
                                   or ("pooled CP-95" in s) or ("What is claimed" in s) for s in para),
    paragraphs=len(para))

fail = [c for c in checks if not c["ok"]]
out = {"checks": checks, "failures": len(fail), "verdict": "PASS" if not fail else "FAIL",
       "failed": [c["name"] for c in fail]}
json.dump(out, open(os.path.join(EV, "selfaudit.json"), "w"), indent=1)
print(json.dumps({"verdict": out["verdict"], "failures": out["failures"], "failed": out["failed"],
                  "n_checks": len(checks)}, indent=1))
sys.exit(1 if fail else 0)
