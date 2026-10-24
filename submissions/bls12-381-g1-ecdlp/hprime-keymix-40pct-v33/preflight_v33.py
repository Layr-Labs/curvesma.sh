#!/usr/bin/env python3
"""v33 preflight: intake laws checked BEFORE the CLI packs the directory."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUB_REL = "submissions/bls12-381-g1-ecdlp/hprime-keymix-40pct-v33"
fails = []


def chk(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ((" :: " + str(detail)) if detail else ""))
    if not cond:
        fails.append(name)


entries, files = [], []
for dp, dn, fn in os.walk(HERE):
    for d in dn:
        entries.append(os.path.join(dp, d))
    for f in fn:
        p = os.path.join(dp, f)
        entries.append(p)
        files.append(p)
nbytes = sum(os.path.getsize(f) for f in files)
biggest = max((os.path.getsize(f), f) for f in files)
chk("entry_count_le_28", len(entries) <= 28, len(entries))
chk("file_count_le_28", len(files) <= 28, len(files))
chk("packed_bytes_le_262144", nbytes <= 262144, nbytes)
chk("per_file_le_60kb", biggest[0] <= 60000, biggest)
chk("no_pycache", not any("__pycache__" in e for e in entries))
chk("no_ds_store", not any(".DS_Store" in e for e in entries))
chk("no_symlinks", not any(os.path.islink(e) for e in entries))
bom = [f for f in files if open(f, "rb").read(3) == b"\xef\xbb\xbf"]
chk("no_bom", not bom, bom)

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

led = json.load(open(os.path.join(HERE, "ledger-v33.json")))
num, den = int(cl[0]["score"]["numerator"]), int(cl[0]["score"]["denominator"])
claim = num / den
import math
chk("claim_reduced_rational", math.gcd(num, den) == 1, (num, den))
chk("claim_matches_ledger", claim == int(led["claim"]["numerator"]) / int(led["claim"]["denominator"]))
chk("claim_ge_W", claim >= led["W_log2"], (claim, led["W_log2"]))
chk("claim_covers_stress", claim - led["W_log2"] >= led["doubling_stress_bits"])
aud = json.load(open(os.path.join(HERE, "audit-v33.json")))
chk("audit_all_pass", aud.get("all_pass") is True)

sel_path = os.path.join(ROOT, "tracks/bls12-381-g1-ecdlp/selected-submission.txt")
sel = open(sel_path).read().strip()
chk("selector_points_here", sel == SUB_REL, sel)

stmt = cl[0]["statement"]
chk("claim_literal_in_statement", "116.8539" in stmt)
chk("gate_literal_in_statement", "116.8271" in stmt)
chk("success_literal_in_statement", "0.395879" in stmt)
chk("gap_ell_literal_in_statement", "2^70.650" in stmt)
print("\n%d failures" % len(fails))
sys.exit(1 if fails else 0)
