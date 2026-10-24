#!/usr/bin/env python3
"""cap_recount.py -- v20: recount the cap rows from the RAW census JSONs.

The success debit that depends on the cap was carried by hand: 1033 in one place
and 1070 in another (would_solve_meet in the escaped_lag3 / cap_1687 cell).  No
number in this file is typed from a note; every count is summed over the cells of
evidence/census-v20n.json and evidence/census-seed11.json:

    successes(cell) = would_solve_meet(cell) - lost(cell)
    pooled          = sum over cells of the SAME arm (never pooled across arms)

It also re-derives the two ledger quantities that consume the row:
    kills_per_1e6_meets = 1e6 * sum(lost) / sum(would_solve_meet)
    p_loss_upper        = Clopper-Pearson one-sided upper on sum(lost)/sum(meets)

Output: evidence/cap-recount.json
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import census_v20n as C   # noqa: E402

FILES = ["evidence/census-v20n.json", "evidence/census-seed11.json"]


def cp_upper(k, n):
    return float(C.cp95_upper(k, n)) if k < n else 1.0


def main():
    raw = {}
    for rel in FILES:
        with open(os.path.join(D, rel)) as f:
            raw[rel] = json.load(f)

    cells, per_arm = [], {}
    for rel, doc in raw.items():
        for arm, caps in sorted(doc["cap_loss"].items()):
            for capname, c in sorted(caps.items()):
                meets = int(c["would_solve_meet"])
                lost = int(c["lost"])
                succ = meets - lost
                cells.append({"file": rel, "arm": arm, "cap": int(c["cap"]),
                              "trials": int(c["trials"]), "would_solve_meet": meets,
                              "lost": lost, "successes": succ,
                              "p_loss": (lost / meets) if meets else None})
                a = per_arm.setdefault(arm, {"meets": 0, "lost": 0, "successes": 0,
                                             "cells": 0, "trials": 0})
                a["meets"] += meets
                a["lost"] += lost
                a["successes"] += succ
                a["cells"] += 1
                a["trials"] += int(c["trials"])

    for arm, a in per_arm.items():
        a["kills_per_1e6_meets"] = 1e6 * a["lost"] / a["meets"] if a["meets"] else None
        a["p_loss_pooled"] = a["lost"] / a["meets"] if a["meets"] else None
        a["p_loss_upper_cp95"] = cp_upper(a["lost"], a["meets"]) if a["meets"] else None

    # the shipped cap is 2*K_ESC_MAX*per ... the PRODUCTION cap row is the cap_1687
    # column (the value the ledger's RESTART_WORST_COUNT arithmetic uses); it is the
    # ONLY column with zero kills in every armed cell, and the only one quotable.
    shipped = {}
    for rel, doc in raw.items():
        for arm, caps in doc["cap_loss"].items():
            for capname, c in caps.items():
                if int(c["cap"]) == 1687:
                    shipped.setdefault(arm, []).append(
                        {"file": rel, "would_solve_meet": int(c["would_solve_meet"]),
                         "lost": int(c["lost"]),
                         "successes": int(c["would_solve_meet"]) - int(c["lost"])})

    # every literal that appears in a note or in ledger.py must be one of these
    literals = sorted({int(c["successes"]) for c in cells}
                      | {int(c["would_solve_meet"]) for c in cells})

    out = {
        "test": "cap_rows_recounted_from_raw_census",
        "sources": FILES,
        "per_arm_pooled": per_arm,
        "cap_1687_cells": shipped,
        "success_literals_present_in_raw": literals,
        "cells": cells,
        "mismatch_note": (
            "1033 is NOT a cap_1687 success count in either census file; the raw "
            "escaped_lag3 cap_1687 cells give 1070 (seed 20250) and 1098 (seed 11). "
            "1033 is the escaped_lag1 successes-at-cap row of census-v20n, i.e. a "
            "DIFFERENT arm. The success row is rebuilt as sum(would_solve_meet - lost) "
            "over the cells of the arm it belongs to."),
    }
    with open(os.path.join(D, "evidence", "cap-recount.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    for arm, a in sorted(per_arm.items()):
        print("%-16s cells=%d meets=%-6d lost=%-5d succ=%-6d p_loss=%.5f upper=%.5f"
              % (arm, a["cells"], a["meets"], a["lost"], a["successes"],
                 a["p_loss_pooled"] or 0.0, a["p_loss_upper_cp95"] or 0.0))
    print("literals:", literals)


if __name__ == "__main__":
    main()
