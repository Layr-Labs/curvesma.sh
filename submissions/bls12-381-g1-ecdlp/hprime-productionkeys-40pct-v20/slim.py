#!/usr/bin/env python3
"""Trim the v20 package to the cited fields and compact JSON (size budget).

Removes the per-challenge duplicate rows from the key census (the averaged rows
and the first challenge's multiplicity histogram are what the argument cites),
trims the gap rows to their cited fields, drops deep-evidence.json (strict subset
of certificate.json) and re-serialises every JSON compactly.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))


def dump(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, separators=(",", ":"), sort_keys=True)


ck = json.load(open(HERE + "/census-key-v20.json"))
for sz in ck["faithful"].values():
    sz.pop("challenges", None)
    for r0 in sz["exact_rows"]:
        r0.pop("within_side_b", None)
dump(HERE + "/census-key-v20.json", ck)

cr = json.load(open(HERE + "/census-rest-v20.json"))
KEEP = ("scale_bits", "d", "seed", "n", "mean_gap", "mean_gap_bits", "max_gap",
        "max_gap_bits", "max_over_mean", "ks", "ks_cut", "ks_pass")
cr["gap"]["rows"] = [{k: r[k] for k in KEEP} for r in cr["gap"]["rows"]]
dump(HERE + "/census-rest-v20.json", cr)

cert = json.load(open(HERE + "/certificate.json"))
dump(HERE + "/certificate.json", cert)
for p in ("deep-evidence.json",):
    if os.path.exists(HERE + "/" + p):
        os.remove(HERE + "/" + p)
print("trimmed")
