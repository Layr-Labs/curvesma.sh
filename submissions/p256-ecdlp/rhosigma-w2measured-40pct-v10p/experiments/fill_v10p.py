#!/usr/bin/env python3
"""fill_v10p.py -- substitute the MEASURED census numbers into every token in the package.

Single source of truth: evidence/pair_census_v10p.json (written by pair_census_pool.py).
Any token left unfilled is a hard error, so a document can never ship with a placeholder.
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
EV = os.path.join(PKG, "evidence")

pc = json.load(open(os.path.join(EV, "pair_census_v10p.json")))
S = pc["shipped_guarded|a=0"]
C = pc["control_pre_parity|a=0"]
S3 = pc["shipped_guarded|a=-3"]
C3 = pc["control_pre_parity|a=-3"]
PER = 50764
W = 16


def runs(a):
    return a["trials"] + (a["trials"] - a["run1_solved"])


def steps(a):
    return W * PER * runs(a)


cal = pc["cp_calibration"]
tokens = {
    "SHIP_N": "%d" % S["trials"],
    "CTRL_N": "%d" % C["trials"],
    "SHIP_P1": "%.5f" % S["p1_single_run"],
    "SHIP_P2": "%.5f" % S["p2_set"],
    "SHIP_CP": "%.5f" % S["p2_cp95_lower"],
    "SHIP_AGREE": "%.5f" % S["formula_1_minus_1mp1_sq"],
    "SHIP_R2C": "%.5f" % (S["run2_given_run1_failed"] or 0.0),
    "CTRL_P1": "%.5f" % C["p1_single_run"],
    "CTRL_P2": "%.5f" % C["p2_set"],
    "CTRL_CP": "%.5f" % C["p2_cp95_lower"],
    "CTRL_AGREE": "%.5f" % C["formula_1_minus_1mp1_sq"],
    "A3_SHIP_P1": "%.5f" % S3["p1_single_run"],
    "A3_SHIP_P2": "%.5f" % S3["p2_set"],
    "A3_SHIP_CP": "%.5f" % S3["p2_cp95_lower"],
    "A3_CTRL_P1": "%.5f" % C3["p1_single_run"],
    "A3_CTRL_P2": "%.5f" % C3["p2_set"],
    "A3_CTRL_CP": "%.5f" % C3["p2_cp95_lower"],
    "A3_N": "%d" % S3["trials"],
    "SHIP_STEPS": "%d" % steps(S),
    "CTRL_STEPS": "%d" % steps(C),
    "DBL_STEPS_TOTAL": "%d" % (steps(S) + steps(C)),
    "A3_GATE": "%.4f (%s of %s)" % (C3["verified"] / C3["meets"], format(C3["verified"], ","),
                                    format(C3["meets"], ",")),
    "A3_STEPS": "%d" % steps(S3),
    "SHIP_DBL": "%d" % S["doubles"],
    "A3_DBL": "%d" % S3["doubles"],
    "DBL_RATE": "%.3e" % (S["doubles"] / steps(S)),
    "DBL_UB": "%.3e" % (3.0 / steps(S)),
    "DBL3_RATE": "%.3e" % (S3["doubles"] / steps(S3)),
    "CAL_COV": "%.3f" % cal["shipped_arm"]["coverage_one_sided"],
    "CAL_SIMS": "%d" % cal["shipped_arm"]["sims"],
    "CAL_ERR": "%.1e" % max(cal["shipped_arm"]["k_equals_n_abs_err"],
                            cal["control_arm"]["k_equals_n_abs_err"]),
    "SHIP_MEETS": "%d" % S["meets"],
    "SHIP_VERIFIED": "%d" % S["verified"],
    "PAIR_SHARDS": "%d" % S["shards"],
}
# audit tokens are filled after the audit runs; keep them out of the doc pass
AUDIT_TOKENS = ("AUDIT_N", "AUDIT_FAIL", "CUR_N", "AUDIT_VERDICT")

docs = ["public-note.md", "argument.md", "submission.json"]
if len(sys.argv) > 1 and sys.argv[1] == "--audit":
    au = json.load(open(os.path.join(EV, "arith_audit_v10p.json")))
    tokens.update({
        "AUDIT_N": "%d" % au["regex_pass"]["n_expressions"],
        "AUDIT_FAIL": "%d" % au["regex_pass"]["n_fail"],
        "CUR_N": "%d" % au["curated_pass"]["n_assertions"],
        "AUDIT_VERDICT": au["verdict"],
    })

second = len(sys.argv) > 1 and sys.argv[1] == "--audit"
for fn in docs:
    p = os.path.join(PKG, fn)
    t = open(p).read()
    for k, v in tokens.items():
        t = t.replace(k, v)
    open(p, "w").write(t)
    bad = [k for k in sorted(set(tokens) | set(AUDIT_TOKENS)) if k in t]
    if not second:
        bad = [k for k in bad if k not in AUDIT_TOKENS]
    print(fn, "remaining tokens:", bad if bad else "none")
    if bad:
        sys.exit(1)

stmt = json.load(open(os.path.join(PKG, "submission.json")))["claims"][0]["statement"]
print("claims[0] chars:", len(stmt), "cap 7965", "OK" if len(stmt) <= 7965 else "TOO LONG")
print("shipped p2 %.5f cp %.5f | control p2 %.5f cp %.5f" %
      (S["p2_set"], S["p2_cp95_lower"], C["p2_set"], C["p2_cp95_lower"]))