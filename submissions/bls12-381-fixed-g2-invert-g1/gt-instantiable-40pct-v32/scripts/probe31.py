#!/usr/bin/env python3
"""probe31.py -- MEASURED live occupancy of the shipped wide store.

The production bucket count is sized from the steady-state occupancy the store
actually holds, not from the number of arrivals: with the shipped chain-full
eviction the store's live set is bounded by 4 x buckets, so the sizing question
is "how many records must be live to reproduce the census retention", and that
number is measured here by d28_ws_live_max on the same fixtures.
"""
import ctypes
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol28 as P                      # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")
BASE = dict(Rg=65536, dpbits=0, rec_write=688, seed0=424242)
SPECS = [("A", dict(q=523, rp=521041, depth=6, workers=8, ratio=1.4493, evict=True)),
         ("B", dict(q=547, rp=10566001, depth=6, workers=8, ratio=1.4493, evict=True))]


def cfg_of(spec):
    c = dict(BASE, **spec)
    N = (c["rp"] - 1) // 12
    c["steps_budget"] = int(c["ratio"] * math.isqrt(N))
    c["trail"] = max(1, -(-c["steps_budget"] // c["workers"]))
    c["cap_ir"] = int(c["steps_budget"] * (P.MUL_CHARGE + P.LOOKUP + c["rec_write"]) +
                      4 * (2 * c["trail"] * P.MUL_CHARGE) + P.VERIFY_RESERVE)
    return c


out = {"mode": "probe31-occupancy",
       "what": "peak live records (d28_ws_live_max) and end occupancy of the SHIPPED "
               "wide store on the census fixtures, per scale",
       "scales": []}
for tag, spec in SPECS:
    c = cfg_of(spec)
    lib = P.EVICTING if c["evict"] else P.CONTROL
    lives, alive, slots, ok = [], [], [], 0
    for gid in range(8):
        r = P.run_one((gid, c))
        lm = ctypes.c_uint64.in_dll(lib, "d28_ws_live_max").value
        live = ctypes.c_uint64.in_dll(lib, "d28_ws_live").value
        lives.append(lm); alive.append(live); slots.append(r["tab_slots"])
        ok += 1 if r["solved"] else 0
    row = {"scale": tag, "fixtures": 8, "solved": ok,
           "tab_slots": slots[0], "records_per_fixture": r["records"] // 8,
           "live_max": max(lives), "live_max_all": lives,
           "live_end": max(alive),
           "deep": ctypes.c_uint64.in_dll(lib, "d28_ws_deep").value,
           "evictions_last": ctypes.c_uint64.in_dll(lib, "d28_evictions").value}
    out["scales"].append(row)
    print(row)
json.dump(out, open(os.path.join(EV, "probe31.json"), "w"), indent=1)
print("wrote evidence/probe31.json")
