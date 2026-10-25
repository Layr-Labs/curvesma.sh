#!/usr/bin/env python3
"""arm37.py -- v38: VERIFIED recoveries only, independent seeds, no fitted constant.

Verdict 4558f762 flaw 2: v37 reported `matches` (denominator nonzero AND
recover() non-None) as "recovered", while protocol37 incremented matches BEFORE
testing g^xx == z, and run_one stopped at the first success -- so at most three
verified recoveries existed in total.  Here `verified` is incremented ONLY when
g^xx == z, the loop runs the full budget (success does not terminate it), and the
tally is per independent seed with the seed value quoted.

The run is EVIDENCE THAT THE PATH WORKS (does the shipped loop solve instances?),
not an estimate of a transfer constant: success given a genuine collision is
1 - 1/r by the relation itself, and that is what the ledger prices.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol37 as P                                   # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")
BASE = dict(Rg=1024, dpbits=0, rec_write=688, seed0=424242, q=523, rp=521041,
            depth=6, workers=8, retain=256, evict=True, steps_budget=40000,
            trail=64, cap_ir=10 ** 17)


def arm(trail, seeds=4):
    rows = []
    for s in range(seeds):
        cfg = dict(BASE)
        cfg["trail"] = trail
        r = P.run_one((s, cfg))
        rows.append({"seed": BASE["seed0"] + s, "steps": r["steps"],
                     "hits": r["within_trail_hits"] + r["cross_trail_hits"],
                     "den_nonzero": r["matches"], "verified": r["verified"],
                     "solved": r["solved"], "dead_hits": r["dead_hits"],
                     "aborts": r["aborts"],
                     "addr_walks": r["addr_walks"], "addr_writes": r["addr_writes"]})
    tot = {k: sum(x[k] for x in rows) for k in ("steps", "hits", "den_nonzero", "verified", "aborts")}
    return {"trail": trail, "seeds": [x["seed"] for x in rows], "per_seed": rows, **tot,
            "solved_seeds": sum(1 for x in rows if x["solved"]),
            "metric": "verified = g^xx == z, counted per hit, loop NOT terminated by success"}


def main():
    out = {"mode": "arm38-verified-recovery",
           "protocol_path": "scripts/protocol37.run_one -- the SHIPPED walk, store, recovery and "
                            "the SAME loop the claim prices (no replacement loop)",
           "store": "shipped dict28.c radix store; retain shrunk to 256 so collisions occur at "
                    "census scale (disclosed; same code path)",
           "arms": [arm(64, 4), arm(65536, 4)],
           "metric_note": "matches (denominator nonzero) is reported SEPARATELY from verified",
           "floor": 0.39}
    out["success_rule"] = ("success given a genuine collision = 1 - 1/r (exact relation); this arm "
                           "is executed evidence, not a fitted constant")
    out["verified_total"] = sum(a["verified"] for a in out["arms"])
    out["solved_seeds"] = sum(a["solved_seeds"] for a in out["arms"])
    json.dump(out, open(os.path.join(EV, "arm37.json"), "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("verified_total", "solved_seeds")}, indent=1))


if __name__ == "__main__":
    main()
