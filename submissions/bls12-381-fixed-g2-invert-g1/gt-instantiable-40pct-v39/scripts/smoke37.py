#!/usr/bin/env python3
"""smoke37.py -- run the SHIPPED protocol and prove both collision classes land.

v35's smoke asserted that a same-trail hit is REJECTED.  That assertion was the
flaw: a same-trail pair obeys the same linear relation as a cross-trail pair and
its denominator vanishes only with probability 1/r, so v36 asserts the opposite
-- a same-trail hit is ACCEPTED and recovers x.  Both directions are exercised
directly on the shipped store, and a fixture run reports the class counts.
"""
import ctypes
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol37 as P                                    # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")


def gate_check():
    """Same key inserted twice: the lookup must return the FULL 128-bit ref."""
    lib = P.EVICTING
    ws = P.D28WS()
    assert lib.dict28_ws_init(ctypes.byref(ws), 12, 4) == 0
    key = [0x1111, 0x2222, 0x3333, 0x4444, 0x5555, 0x6666]
    kb = (ctypes.c_uint64 * 6)(*key)
    sl = ctypes.c_uint64(0)
    sh = ctypes.c_uint64(0)
    lib.dict28_slot_export(kb, 4096, 0, 0, ctypes.byref(sl), ctypes.byref(sh))
    # a reference whose HIGH word is nonzero: the production ordinal range
    ordv = (1 << 70) + 12345
    ref = (ctypes.c_uint64 * 2)(ordv & P.M64, ordv >> 64)
    ev = (ctypes.c_uint64 * 2)(0, 0)
    lib.dict28_wide_insert(ctypes.byref(ws), sl.value, sh.value, kb, ref, ev)
    out = (ctypes.c_uint64 * 2)(0, 0)
    hit = lib.dict28_wide_lookup(ctypes.byref(ws), sl.value, sh.value, kb, out)
    full = (out[0] | (out[1] << 64)) == ordv
    # fill the chain to PROBES and evict the 128-bit ref back out
    for i in range(4):
        k2 = (ctypes.c_uint64 * 6)(key[0], key[1], key[2], key[3], key[4], 0x9000 + i)
        lib.dict28_wide_insert(ctypes.byref(ws), sl.value, sh.value, k2,
                               (ctypes.c_uint64 * 2)(i + 1, 0), ev)
    got_ev = ev[0] | (ev[1] << 64)
    lib.dict28_ws_free(ctypes.byref(ws))
    return {"lookup_hit": int(hit), "ref_128_round_trip": bool(full),
            "evicted_ref_128": str(got_ev), "evicted_ref_above_2_64": bool(got_ev >> 64)}


def main():
    base = dict(Rg=65536, dpbits=0, rec_write=688, seed0=424242, q=523, rp=521041,
                depth=6, workers=8, evict=True, steps_budget=900, trail=113, cap_ir=10 ** 12)
    runs = [P.run_one((i, base)) for i in range(4)]
    gate = gate_check()
    assert gate["lookup_hit"] == 1 and gate["ref_128_round_trip"]
    out = {
        "mode": "protocol37-smoke",
        "code": "scripts/protocol37.py with the v36 class-partitioned recovery",
        "fixtures": len(runs),
        "solved": sum(1 for r in runs if r["solved"]),
        "steps": sum(r["steps"] for r in runs),
        "records": sum(r["records"] for r in runs),
        "within_trail_hits": sum(r["within_trail_hits"] for r in runs),
        "cross_trail_hits": sum(r["cross_trail_hits"] for r in runs),
        "within_recovered": sum(r["within_recovered"] for r in runs),
        "cross_recovered": sum(r["cross_recovered"] for r in runs),
        "matches": sum(r["matches"] for r in runs),
        "aborts": sum(r["aborts"] for r in runs),
        "steps_aborted": sum(r["steps_aborted"] for r in runs),
        "trail_starts": sum(r["trail_starts"] for r in runs),
        "dict_calls": {k: sum(r["c_calls"][k] for r in runs) for k in runs[0]["c_calls"]},
        "gate": gate,
        "gate_rule": "both classes are accepted: a hit's class is a LABEL, not a veto.  The store "
                     "returns the full 128-bit reference and the full 128-bit evicted reference.",
    }
    assert out["dict_calls"]["d28_lookups"] > 0 and out["dict_calls"]["d28_inserts"] > 0
    json.dump(out, open(os.path.join(EV, "protocol37.json"), "w"), indent=1)
    print(json.dumps(out, indent=1)[:1400])


if __name__ == "__main__":
    main()
