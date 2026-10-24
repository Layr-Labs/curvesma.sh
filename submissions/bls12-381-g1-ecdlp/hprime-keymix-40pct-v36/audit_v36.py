#!/usr/bin/env python3
"""audit_v36.py -- the controls that can be run against THIS package, plus the packing
checks.  Controls (b) executed-vs-static and (c) borrow propagation from the lane brief
are NOT applicable: v36 makes no compiled-body or field-arithmetic claim (see
ledger-v36.anchors), so there is no static price to unroll and no limb subtractor to
test.  Saying that is the audit, not skipping it.

Runs:
  (a) lifting control       corrupt one bit of s0 -> zero passing candidates
  (d) hard-coded-ladder     perturb a raw counter -> the ladder output must change
  (e) independent recompute of the ladder counts from the written cell data
  (f) evidence resolution + packed-byte sum + reduced rational
"""
import json, math, os, random, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ec, keymod, ladder_v36 as LAD

J = lambda f: json.load(open(os.path.join(HERE, f)))
out = {"kind": "v36 audit"}

# (a) lifting control
cv, G, r = ec.Curve(ec.P37), ec.G37, ec.R37
rng = random.Random(7)
D = 10
m = (r - 1) >> D
bad_pass, good_pass = 0, 0
for _ in range(20):
    s = rng.randrange(1, r)
    T = cv.mul(s, G)
    good_pass += bool(keymod.lift(cv, G, T, s % m, m, D)["found"])
    bad_pass += bool(keymod.lift(cv, G, T, (s % m) ^ (1 << rng.randrange(0, m.bit_length())), m, D)["found"])
out["lifting_control"] = dict(trials=20, correct_seed_found=good_pass, corrupted_seed_found=bad_pass,
                             verdict="PASS" if (good_pass == 20 and bad_pass == 0) else "FAIL")

# (d) hard-coded-ladder control: perturb a raw pool count and require the output to move
cells = J("structured-ladder.json")["cells"]
base = cells["nb40_ell50_b12"]
pert = dict(base, successes=base["successes"] - 1, pools=base["pools"])
q0 = base["successes"] / base["pools"]
q1 = pert["successes"] / pert["pools"]
cp0 = LAD.cp_lower(base["successes"], base["pools"])
cp1 = LAD.cp_lower(pert["successes"], pert["pools"])
out["hardcoded_ladder_control"] = dict(
    cell="nb40_ell50_b12", original=round(q0, 6), perturbed=round(q1, 6),
    cp95_moves=bool(abs(cp0 - cp1) > 0), delta=round(cp0 - cp1, 9),
    verdict=("PASS: the endpoint is computed from the counted successes at run time, so moving a "
             "raw count moves the output" if abs(cp0 - cp1) > 0 else "FAIL: endpoint is not a function of the counts"))

# (e) independent recompute of every ladder cell from its own recorded counts
bad = []
for k, v in cells.items():
    q = v["successes"] / v["pools"]
    if abs(q - v["q"]) > 1e-9 or abs(LAD.cp_lower(v["successes"], v["pools"]) - v["cp95_lower"]) > 1e-6:
        bad.append(k)
out["independent_recompute"] = dict(cells=len(cells), mismatches=bad,
                                    verdict="PASS" if not bad else "FAIL: %r" % bad)

# (f) packing
sub = J("submission.json")
missing = [e for e in sub["claims"][0]["evidence"] if not os.path.exists(os.path.join(HERE, e))]
files = sorted(f for f in os.listdir(HERE) if os.path.isfile(os.path.join(HERE, f)))
total_bytes = sum(os.path.getsize(os.path.join(HERE, f)) for f in files)
led = J("ledger-v36.json")
from fractions import Fraction as F
num, den = int(led["claim"]["numerator"]), int(led["claim"]["denominator"])
out["packing"] = dict(evidence_missing=missing, files=len(files), packed_bytes=total_bytes,
                      under_250KB=bool(total_bytes <= 250 * 1024), under_28_files=bool(len(files) <= 28),
                      claim_reduced=bool(math.gcd(num, den) == 1),
                      claim_value=round(num / den, 4),
                      claim_matches_statement=bool(("(%d/%d)" % (num, den)) in sub["claims"][0]["statement"]))
out["not_applicable"] = dict(
    executed_vs_static="no compiled-body price is claimed in v36 (anchors carried=False)",
    borrow_propagation="no limb-subtractor is claimed in v36 (anchors carried=False)")
ok = (out["lifting_control"]["verdict"] == "PASS" and out["hardcoded_ladder_control"]["verdict"].startswith("PASS")
      and not bad and not missing and out["packing"]["under_250KB"] and out["packing"]["under_28_files"]
      and out["packing"]["claim_reduced"] and out["packing"]["claim_matches_statement"])
out["verdict"] = "PASS" if ok else "FAIL"
json.dump(out, open(os.path.join(HERE, "audit-v36.json"), "w"), indent=1, sort_keys=True)
print(json.dumps(out, indent=1))
sys.exit(0 if ok else 1)
