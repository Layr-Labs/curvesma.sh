#!/usr/bin/env python3
"""v32 preflight: intake laws checked BEFORE the CLI packs the dir."""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUB_REL = "submissions/bls12-381-g1-ecdlp/hprime-keymix-40pct-v32"
fails = []


def chk(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ((" :: " + str(detail)) if detail else ""))
    if not cond:
        fails.append(name)


entries = []
for dp, dn, fn in os.walk(HERE):
    for d in dn:
        entries.append(os.path.join(dp, d))
    for f in fn:
        entries.append(os.path.join(dp, f))
files = [e for e in entries if os.path.isfile(e)]
nbytes = sum(os.path.getsize(f) for f in files)
chk("entry_count_le_28", len(entries) <= 28, len(entries))
chk("packed_bytes_le_262144", nbytes <= 262144, nbytes)
chk("no_pycache", not any("__pycache__" in e for e in entries))
chk("no_ds_store", not any(".DS_Store" in e for e in entries))

sub = json.load(open(os.path.join(HERE, "submission.json")))
cl = sub["claims"]
chk("claims_present", len(cl) == 1)
chk("claim_text_le_7965", len(cl[0]["statement"]) <= 7965, len(cl[0]["statement"]))
chk("title_le_300", len(sub["title"]) <= 300, len(sub["title"]))
chk("experiments_is_list", isinstance(sub.get("experiments"), list))
chk("experiments_dir_nonempty", os.path.isdir(os.path.join(HERE, "experiments")) and
    len(os.listdir(os.path.join(HERE, "experiments"))) > 0)
missing = [e for e in cl[0]["evidence"] if not os.path.exists(os.path.join(HERE, e))]
chk("evidence_files_exist", not missing, missing)
chk("evidence_count_le_32", len(cl[0]["evidence"]) <= 32, len(cl[0]["evidence"]))

led = json.load(open(os.path.join(HERE, "ledger-v32.json")))
claim = int(cl[0]["score"]["numerator"]) / int(cl[0]["score"]["denominator"])
chk("claim_matches_ledger", abs(claim - int(led["claim"]["numerator"]) / int(led["claim"]["denominator"])) < 1e-12)
chk("claim_ge_W", claim >= led["W_log2"], (claim, led["W_log2"]))
chk("claim_covers_stress", claim - led["W_log2"] >= led["doubling_stress_bits"])
aud = json.load(open(os.path.join(HERE, "audit-v32.json")))
chk("audit_all_pass", aud.get("all_pass") is True)

sel = open(os.path.join(ROOT, "tracks/bls12-381-g1-ecdlp/selected-submission.txt")).read().strip()
chk("selector_points_here", sel == SUB_REL, sel)

# literal sweep: every 4+ digit decimal in the claim text that purports to be the claim
chk("claim_literal_in_statement", "116.8543" in cl[0]["statement"])
chk("gate_literal_in_statement", "116.8271" in cl[0]["statement"])

print("\n%d checks, %d failures" % (11 + 8, len(fails)))
sys.exit(1 if fails else 0)
