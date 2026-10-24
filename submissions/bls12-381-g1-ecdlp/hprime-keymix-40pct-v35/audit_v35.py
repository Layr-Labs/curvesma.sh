#!/usr/bin/env python3
"""audit_v35.py -- the four controls the repair asked for, each able to FAIL.

1. reduction control   : recompile op_fp381_mulmod with the reduction deleted and assert the
                         disassembled instruction count CHANGES.  If it does not change, the
                         reduction is not in the body (v34's exact defect) and this fails.
2. evidence control    : every claims[0].evidence path must resolve INSIDE this directory --
                         no '..', no absolute path, file must exist.  This is the intake killer
                         from ticket 3056f3d5, asserted instead of trusted.
3. rung control        : the number of rungs below the floor inside the claim population must
                         equal the number the claim text states.
4. literal control     : the bootstrap interval in submission.json must equal the interval in
                         structured-ladder.json, to the last quoted digit.
"""
import json, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))
fails = []


def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name + ((" :: " + str(detail)) if detail else ""))
    if not ok:
        fails.append(name)


# 1. reduction removal control ------------------------------------------------
src = open(os.path.join(HERE, "price_body_v35.c")).read()
broken = re.sub(r"for \(int bit = 767; bit >= 0; bit--\) \{.*?\n    \}\n", "",
                src, count=1, flags=re.S)
assert "bit = 767" in src and "bit = 767" not in broken, "could not construct the control"
tmp = tempfile.mkdtemp()
open(os.path.join(tmp, "b.c"), "w").write(broken)
subprocess.run(["clang", "-O2", "-c", os.path.join(tmp, "b.c"), "-o", os.path.join(tmp, "b.o")],
               check=True, capture_output=True)
dis = subprocess.run(["otool", "-tV", os.path.join(tmp, "b.o")], capture_output=True, text=True).stdout


def count_symbol(text, sym):
    cur, n = None, 0
    for ln in text.splitlines():
        s = ln.strip()
        m = re.match(r"^(_[A-Za-z0-9_]+):$", s)
        if m:
            cur = m.group(1)
            continue
        if cur == sym and re.match(r"^[0-9a-f]{16}\s", s):
            n += 1
    return n


real = count_symbol(open(os.path.join(HERE, "disasm_v35.txt")).read(), "_op_fp381_mulmod")
ctl = count_symbol(dis, "_op_fp381_mulmod")
check("reduction_removed_control", real != ctl,
      "with reduction %d insns vs without %d (delta %d)" % (real, ctl, real - ctl))

# 2. evidence paths resolve inside the packed directory -----------------------
sub = J("submission.json")
bad = []
for ev in sub["claims"][0].get("evidence", []):
    if ".." in ev or os.path.isabs(ev) or not os.path.exists(os.path.join(HERE, ev)):
        bad.append(ev)
check("evidence_inside_tree", not bad, bad or "%d/%d resolve inside" % (
    len(sub["claims"][0]["evidence"]), len(sub["claims"][0]["evidence"])))

# 3. rung count control -------------------------------------------------------
lad = J("structured-ladder.json")
stated = sub["claims"][0]["statement"]
m = re.search(r"(two|three|four|2|3|4) of the five long-walk", stated)
words = {"two": 2, "three": 3, "four": 4}
claim_n = words.get(m.group(1), int(m.group(1)) if m and m.group(1).isdigit() else None) if m else None
actual = lad["rungs_below_floor"]["in_claim_population"]
check("rung_count_control", claim_n == actual, "stated %s vs measured %s" % (claim_n, actual))

# 4. bootstrap literal control -----------------------------------------------
led = J("ledger-v35.json")
bl = lad["cluster_bootstrap_over_rungs"]
quoted = re.findall(r"\[(0\.\d+), (0\.\d+)\]", stated)
ok_lit = bool(quoted) and all(abs(float(a) - bl["lo"]) < 5e-7 and abs(float(b) - bl["hi"]) < 5e-7
                             for a, b in quoted)
check("bootstrap_literal_control", ok_lit,
      "quoted %s vs ladder [%s, %s]" % (quoted, bl["lo"], bl["hi"]))
check("single_interval_control", len(quoted) == 1, "%d intervals quoted" % len(quoted))
check("claim_matches_ledger",
      str(sub["claims"][0]["score"]["numerator"]) == led["claim"]["numerator"] and
      str(sub["claims"][0]["score"]["denominator"]) == led["claim"]["denominator"],
      "%s/%s" % (led["claim"]["numerator"], led["claim"]["denominator"]))
check("success_matches_ladder", abs(led["success"]["claimed_value"] - bl["lo"]) < 1e-9 or
      abs(led["success"]["claimed_value"] - lad["pooled_long_walks"]["cp95_lower"]) < 1e-9,
      str(led["success"]["claimed_value"]))
print("\naudit: %d checks, %d failed" % (7, len(fails)))
sys.exit(1 if fails else 0)
