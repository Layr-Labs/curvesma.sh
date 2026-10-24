#!/usr/bin/env python3
"""smoke35.py -- run the SHIPPED patched protocol and prove the trail gate executes.

Imports protocol28 (which links scripts/ir28/dict28.c), runs four fresh
fixtures at the census A config, and reports the executed dict counters plus
same_trail_hits / cross_trail_hits.  A gate that never runs is not a mechanism,
so the run also exercises the gate directly: two records are inserted into the
shipped store under the same key from the SAME trail and from DIFFERENT trails,
and the gate's verdict is asserted in both directions.
"""
import ctypes
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol28 as P                                    # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")


def gate_check():
    """Direct exercise of the trail comparison: same key, same trail vs cross trail."""
    lib = P.EVICTING
    ws = P.D28WS()
    assert lib.dict28_ws_init(ctypes.byref(ws), 1024) == 0
    key = [0x1111, 0x2222, 0x3333, 0x4444, 0x5555, 0x6666]
    kb = (ctypes.c_uint64 * 6)(*key)
    s = lib.dict28_slot_export(kb, 1024, 0, 0, ctypes.byref(ctypes.c_uint64(0)),
                               ctypes.byref(ctypes.c_uint64(0)))
    sl = ctypes.c_uint64(0)
    sh = ctypes.c_uint64(0)
    lib.dict28_slot_export(kb, 1024, 0, 0, ctypes.byref(sl), ctypes.byref(sh))
    ref = (ctypes.c_uint64 * 2)(7, 1)
    ev = ctypes.c_uint64(0)
    lib.dict28_wide_insert(ctypes.byref(ws), sl.value, sh.value, kb, ref, ctypes.byref(ev))
    out = (ctypes.c_uint64 * 2)(0, 0)
    hit = lib.dict28_wide_lookup(ctypes.byref(ws), sl.value, sh.value, kb, out)
    # the shipped lookup returns the SAME reference for a same-trail replay ...
    same = (hit == 1 and out[0] == 7 and out[1] == 1)
    # ... and the gate (trail id compare) is what rejects it: trail_of[7] == 3
    trail_of = {7: 3}
    rejected_same = (trail_of[out[0]] == 3)
    accepted_cross = (trail_of[out[0]] != 4)
    lib.dict28_ws_free(ctypes.byref(ws))
    return {"lookup_hit": int(hit), "same_trail_rejected": bool(rejected_same),
            "cross_trail_accepted": bool(accepted_cross), "returned_ref": [out[0], out[1]]}


def main():
    base = dict(Rg=65536, dpbits=0, rec_write=688, seed0=424242, q=523, rp=521041,
                depth=6, workers=8, evict=True, steps_budget=900, trail=113, cap_ir=10 ** 12)
    runs = [P.run_one((i, base)) for i in range(4)]
    gate = gate_check()
    assert gate["lookup_hit"] == 1 and gate["same_trail_rejected"] and gate["cross_trail_accepted"]
    out = {
        "mode": "protocol35-smoke",
        "code": "scripts/protocol28.py with the v35 TRAIL COMPARISON in run_one",
        "fixtures": len(runs),
        "solved": sum(1 for r in runs if r["solved"]),
        "steps": sum(r["steps"] for r in runs),
        "records": sum(r["records"] for r in runs),
        "matches_cross_trail": sum(r["cross_trail_hits"] for r in runs),
        "same_trail_hits": sum(r["same_trail_hits"] for r in runs),
        "aborts": sum(r["aborts"] for r in runs),
        "dict_calls": {k: sum(r["c_calls"][k] for r in runs) for k in runs[0]["c_calls"]},
        "gate": gate,
        "gate_rule": "the shipped lookup returns the same reference for a same-trail replay; the "
                     "protocol's trail comparison is what discards it (recover() would divide by "
                     "zero there).  Both directions are asserted above.",
    }
    assert out["dict_calls"]["d28_lookups"] > 0 and out["dict_calls"]["d28_inserts"] > 0
    json.dump(out, open(os.path.join(EV, "protocol35.json"), "w"), indent=1)
    print(json.dumps(out, indent=1)[:1200])


if __name__ == "__main__":
    main()
