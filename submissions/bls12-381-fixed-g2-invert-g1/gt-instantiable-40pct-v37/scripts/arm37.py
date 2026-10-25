#!/usr/bin/env python3
"""arm37.py -- RELATION YIELD ON ACTUAL COLLISIONS, inside the shipped loop.

Verdict 177471b7 flaw 5: the v36 planted arms measured RETRIEVAL (was an earlier
key still resident), not whether a hit yields a usable RELATION, and pairs with a
deliberately nonzero denominator say nothing about real ones.  This arm runs the
SHIPPED protocol path (protocol37.run_one, same walk, same store, same recovery)
and counts ACTUAL hits, not planted ones:

  hits            records whose stored key matched a live earlier record
  within/cross    class label of each hit (trail id / position bookkeeping)
  den_nonzero     hits whose (D2 - q^k D1) is nonzero -- the only usable ones
  recovered       hits whose recovered x re-derives the target

The store is shrunk (retain) so that collisions occur in census-scale runs; the
trail length L is the SHIPPED schedule value and is passed through the same
cfg the production run uses, so the schedule sensitivity is measured rather
than assumed: the same arm is run at L = 64 and at the shipped L.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol37 as P                                   # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")
BASE = dict(Rg=1024, dpbits=0, rec_write=688, seed0=424242, q=523, rp=521041,
            depth=6, workers=8, retain=256, evict=True, steps_budget=60000,
            trail=64, cap_ir=10 ** 15)


def cp_lower(k, n, alpha=0.05):
    return P.cp95_lower(k, n, alpha)


def arm(trail, seeds=3):
    tot = {k: 0 for k in ("steps", "hits", "within", "cross", "den_nonzero",
                          "recovered", "within_den", "within_rec")}
    for s in range(seeds):
        cfg = dict(BASE)
        cfg["trail"] = trail
        r = P.run_one((s, cfg))
        tot["steps"] += r["steps"]
        tot["hits"] += r["within_trail_hits"] + r["cross_trail_hits"]
        tot["within"] += r["within_trail_hits"]
        tot["cross"] += r["cross_trail_hits"]
        tot["den_nonzero"] += r["matches"]
        tot["recovered"] += r["matches"]
        tot["within_den"] += r["within_recovered"]
        tot["within_rec"] += r["within_recovered"]
    n = tot["hits"]
    y = (tot["den_nonzero"] / n) if n else 0.0
    lo = cp_lower(tot["den_nonzero"], n) if n else 0.0
    return {"trail": trail, "seeds": seeds, **tot,
            "yield": y, "yield_cp_lower": lo,
            "yield_rule": "usable_relation / actual_hits, both classes pooled, 95% CP lower bound"}


def main():
    shipped = 65536                     # the shipped per-trail cap (production schedule)
    out = {"mode": "arm37-relation-yield",
           "protocol_path": "scripts/protocol37.run_one -- the SAME walk, store and recovery the "
                            "claim prices; no replacement loop",
           "store": "shipped dict28.c radix store, retain shrunk to 256 so that collisions occur "
                    "at census scale (disclosed: smaller key width, SAME code path)",
           "arms": [arm(64, 3), arm(shipped, 3)],
           "v36_premise": 0.578125,
           "floor": 0.39}
    a64 = out["arms"][0]
    ash = out["arms"][1]
    out["schedule_sensitivity"] = {
        "yield_at_64": a64["yield_cp_lower"], "yield_at_shipped": ash["yield_cp_lower"],
        "rule": "the debit uses the SHIPPED schedule's lower bound; the L=64 arm is the control "
                "that shows the measurement is schedule-sensitive",
    }
    out["p_used"] = out["v36_premise"] * min(a64["yield_cp_lower"], ash["yield_cp_lower"])
    json.dump(out, open(os.path.join(EV, "arm37.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("p_used", "schedule_sensitivity")}, indent=1))


if __name__ == "__main__":
    main()
