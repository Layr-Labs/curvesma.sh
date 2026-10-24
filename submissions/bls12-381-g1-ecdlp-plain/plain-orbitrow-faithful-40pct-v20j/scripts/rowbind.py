#!/usr/bin/env python3
"""rowbind.py -- ROW BINDING for the note's evidence assertions (v20i).

WHY THIS EXISTS (reviewer point 4 of 2ba60258).  v20h's self-audit matched note
numbers against a FLAT BAG of every value in the package: any token that hit any
artifact passed.  The note therefore said "the equivariance mutation control
mismatches on 0 comparisons" while evidence/medium-fixture.json held 7494 column
+ 7494 successor = 14988 mismatches in the CONTROL row, and 0 in the FIXED row.
The arithmetic was right and the ASSERTION was bound to the wrong row: a bag
cannot tell "fixed = 0" from "control = 0".

WHAT v20i DOES.  Every evidence row carries (test, arm) as explicit FIELDS, and
the note cites assertions by FIELD PATH -- [E2.fixed.mismatches] -- which this
module resolves as a query, not as a number match:

    rows[(test, arm)][field] == token

A token bound to a path is checked ONLY against that path (never against the
bag), so a right number in the wrong row FAILS.  `--plant-wrong-row` swaps the
fixed and control rows' mismatches and requires the binding to FAIL; a binding
that cannot fail is not evidence.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
OUT = os.path.join(D, "evidence", "rowbind.json")

# note field paths -> (evidence file, test, arm, field)
PATHS = {
    "E2.fixed.mismatches": ("medium-fixture.json", "equivariance", "fixed", "mismatches"),
    "E2.fixed.comparisons": ("medium-fixture.json", "equivariance", "fixed", "comparisons"),
    "E2.control.mismatches": ("medium-fixture.json", "equivariance", "control", "mismatches"),
    "E2.control.comparisons": ("medium-fixture.json", "equivariance", "control", "comparisons"),
    "E2.fixed.column": ("medium-fixture.json", "equivariance", "fixed", "mismatch_column"),
    "E2.control.column": ("medium-fixture.json", "equivariance", "control", "mismatch_column"),
    "E2.fixed.alt": ("medium-fixture.json", "equivariance", "fixed", "mismatch_alt_successor"),
    "E2.control.alt": ("medium-fixture.json", "equivariance", "control", "mismatch_alt_successor"),
    "E3.fixed.mismatches": ("equivariance-rows.json", "equivariance", "fixed", "mismatches"),
    "E3.control.mismatches": ("equivariance-rows.json", "equivariance", "control", "mismatches"),
    "E3.fixed.comparisons": ("equivariance-rows.json", "equivariance", "fixed", "comparisons"),
    "E3.control.comparisons": ("equivariance-rows.json", "equivariance", "control", "comparisons"),
    "E4.fixed.verified": ("equivariance-rows.json", "conversion", "fixed", "verified"),
    "E4.control.verified": ("equivariance-rows.json", "conversion", "control", "verified"),
}


def load(name):
    with open(os.path.join(D, "evidence", name)) as f:
        return json.load(f)


def rows_medium(med=None):
    """The two equivariance arms of medium-fixture.json, as (test, arm) ROWS."""
    med = med or load("medium-fixture.json")
    out = {}
    for arm, key in (("fixed", "equivariance"), ("control", "equivariance_mutation")):
        r = med[key]
        out[("equivariance", arm)] = {
            "test": "equivariance",
            "arm": arm,
            "selector": r["selector"],
            "states": r["states"],
            "sigma_phases": r["sigma_phases"],
            "comparisons": r["comparisons"],
            "mismatches": r["total_mismatches"],
            "mismatch_column": r["mismatches"]["column"],
            "mismatch_alt_successor": r["mismatches"]["alt_successor"],
            "mismatch_decision": r["mismatches"]["decision"],
            "passed": r["passed"],
        }
    for arm, key, want in (("fixed", "conversion", True), ("control", "conversion_mutation", False)):
        r = med[key]
        out[("conversion", arm)] = {
            "test": "conversion",
            "arm": arm,
            "collisions": r["collisions"],
            "verified": r["verified"],
            "corrupt": r["corrupt"],
            "passed": (r["verified"] == r["collisions"]) == want,
        }
    return out


def resolver(rows):
    def get(path):
        _f, test, arm, field = PATHS[path]
        # ALWAYS resolve from the rows dict handed in, so an injected row swap
        # propagates (v20h's bug was resolving against a fresh read of the file).
        return rows[(test, arm)][field]
    return get


BIND = re.compile(r"\[([A-Z]\d+\.[a-z_]+(?:\.[a-z_]+)?)\]\s*([0-9][0-9,]*)")


def bound_assertions(note, get):
    """Check every [path] token against THAT row field only."""
    out = []
    for m in BIND.finditer(note):
        path, tok = m.group(1), m.group(2).replace(",", "")
        try:
            want = get(path)
        except KeyError:
            out.append({"path": path, "token": tok, "value": None,
                        "ok": False, "why": "unknown path"})
            continue
        out.append({"path": path, "token": tok, "value": want,
                    "ok": int(tok) == int(want), "why": ""})
    return out


def main():
    argv = sys.argv[1:]
    note = open(os.path.join(D, "public-note.md")).read()
    rows = rows_medium()
    with open(os.path.join(D, "evidence", "equivariance-rows.json"), "w") as f:
        json.dump({"test": "equivariance+conversion", "arm": "field-labelled",
                   "rows": [v for _k, v in sorted(rows.items())]}, f,
                  indent=1, sort_keys=True)
    res = {"test": "rowbind", "bound": bound_assertions(note, resolver(rows))}
    res["bound_total"] = len(res["bound"])
    res["bound_ok"] = sum(1 for b in res["bound"] if b["ok"])
    res["passed"] = res["bound_total"] > 0 and res["bound_ok"] == res["bound_total"]
    if "--plant-wrong-row" in argv:
        bad = dict(rows)
        f, c = rows[("equivariance", "fixed")], rows[("equivariance", "control")]
        bad[("equivariance", "fixed")] = dict(f, mismatches=c["mismatches"],
                                              mismatch_column=c["mismatch_column"])
        bad[("equivariance", "control")] = dict(c, mismatches=f["mismatches"],
                                                mismatch_column=f["mismatch_column"])
        b2 = bound_assertions(note, resolver(bad))
        res["planted_wrong_row"] = {
            "injected": "fixed<->control mismatches swapped",
            "failed_paths": [b["path"] for b in b2 if not b["ok"]],
            "caught": any(not b["ok"] for b in b2),
        }
        res["passed"] = res["passed"] and res["planted_wrong_row"]["caught"]
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    print(json.dumps(res, indent=1, sort_keys=True))
    raise SystemExit(0 if res["passed"] else 1)


if __name__ == "__main__":
    main()
