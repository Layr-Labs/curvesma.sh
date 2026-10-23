#!/usr/bin/env python3
"""Render the prose artifacts from render.json + toy-summary.json.

Single source of truth: every number in public-note.md, argument.md and submission.json is
substituted from build_certificate.py's render.json or from the raw toy summary.  Any token
left unfilled is a hard error.
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
r = json.load(open(HERE + "/render.json"))
ts = json.load(open(HERE + "/toy-summary.json"))
fit = ts["gapFit"]["sizes"]

val = {k.upper(): str(v) for k, v in r.items()}
val["SCALEMAX"] = str(r["scale_max"])
for a, b in (("CELLROWS","cell_rows"),("CLAIMNUM","claim_num"),("CLAIMDEN","claim_den"),("DEBITSUM","debit_sum"),("FITROWS","fit_rows"),("RESIDROWS","resid_rows"),("SENSROWS","sens_rows"),("GAPMEAN","gapmean"),("GAPMAX","gapmax"),("GAPLAW","gaplaw"),("ALLOWA_FRAC","allowa_frac"),("ALLOWB_FRAC","allowb_frac"),("ALLOW_SUM","allow_sum")):
    val[a] = str(r[b])
val["FITMEANS"] = "/".join("%.3f" % s["mean_over_law"] for s in fit)
val["FITKS"] = "/".join("%.4f" % s["ks"] for s in fit)

for src, dst in (("submission.json.tpl", "submission.json"),
                 ("public-note.md.tpl", "public-note.md"),
                 ("argument.md.tpl", "argument.md")):
    txt = open(HERE + "/" + src).read()
    for k, v in val.items():
        txt = txt.replace("__%s__" % k, v)
    left = sorted(set(re.findall(r"__[A-Z0-9_]+__", txt)))
    if left:
        print("UNFILLED in %s: %s" % (src, left))
        sys.exit(1)
    open(HERE + "/" + dst, "w").write(txt)
    print("%-18s -> %-18s %d bytes" % (src, dst, len(txt)))

if len(open(HERE + "/public-note.md").read()) < 6000:
    print("NOTE TOO SMALL")
    sys.exit(1)
print("render ok")
