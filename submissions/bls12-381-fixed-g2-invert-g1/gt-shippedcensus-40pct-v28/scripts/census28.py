#!/usr/bin/env python3
"""census28.py -- THE CENSUS AT THE SHIPPED CONFIG (v28).

Repairs verdict 082e0906 flaw 1 and flaw 5: v27's census ran depth 5, eight
workers and a PYTHON dictionary, while the claim declared depth 6, 1024 workers
and an evicting dictionary, and v27 left the eviction loss unmodelled.

What this driver measures, all through the shipped objects (scripts/ir28/*.dylib
linked from dict28.o, ctypes, call counters asserted nonzero):

  A  rp=521041  depth 6  W=8    evicting      ratio 1.4493  320 fixtures
  B  rp=10566001 depth 6 W=8    evicting      ratio 1.4493  320 fixtures
  C  rp=521041  depth 5  W=8    evicting      ratio 1.4493  128   (depth contrast)
  D  rp=521041  depth 6  W=64   evicting      ratio 4.0     128   (partition)
  E  rp=521041  depth 6  W=128  evicting      ratio 4.0     128   (partition)
  F  rp=521041  depth 6  W=8    evicting      ratio 4.0     128   (partition baseline)
  G  rp=521041  depth 6  W=8    NO-EVICT ctl  ratio 1.4493  320   (eviction loss)

Worker-count fidelity is stated, not asserted: at fixture scale the trail is
ceil(ratio*sqrt(N)/W), so W is bounded by the requirement trail >= ring (6) for
the guard statistics to be at stationarity; the largest feasible W is 128 at
ratio 4.0 (trail 7).  Production's trail is 9.36e34.  The partition series
D/E/F measures whether p_hat depends on W at FIXED total step budget S, which is
the only W-dependence the census can carry; the transfer of the per-step
statistics to production's trail rests on stationarity, measured separately by
the guard/trap rates in trap28.json.
"""
import json
import math
import multiprocessing as mp
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import protocol28 as P                                    # noqa: E402

EV = os.path.join(os.path.dirname(HERE), "evidence")
RING = 6


def _work(job):
    cfg, i, name = job
    r = P.run_one((i, cfg))
    r["cfg"] = name
    return r


def main():
    t0 = time.time()
    base = dict(Rg=65536, dpbits=0, rec_write=688, seed0=424242)
    specs = [
        ("A_rp521041_d6_w8_evict", dict(q=523, rp=521041, depth=6, workers=8,
                                        ratio=1.4493, evict=True), 320),
        ("B_rp10566001_d6_w8_evict", dict(q=547, rp=10566001, depth=6, workers=8,
                                          ratio=1.4493, evict=True), 320),
        ("C_rp521041_d5_w8_evict", dict(q=523, rp=521041, depth=5, workers=8,
                                        ratio=1.4493, evict=True), 128),
        ("D_rp521041_d6_w64_evict_r4", dict(q=523, rp=521041, depth=6, workers=64,
                                            ratio=4.0, evict=True), 128),
        ("E_rp521041_d6_w128_evict_r4", dict(q=523, rp=521041, depth=6, workers=128,
                                             ratio=4.0, evict=True), 128),
        ("F_rp521041_d6_w8_evict_r4", dict(q=523, rp=521041, depth=6, workers=8,
                                           ratio=4.0, evict=True), 128),
        ("G_rp521041_d6_w8_NOEVICT", dict(q=523, rp=521041, depth=6, workers=8,
                                          ratio=1.4493, evict=False), 320),
    ]
    cfgs, all_jobs, meta = {}, [], []
    for name, spec, nfix in specs:
        c = dict(base, **spec)
        c["name"] = name
        c["fixtures"] = nfix
        N = (c["rp"] - 1) // 12
        c["steps_budget"] = int(c["ratio"] * math.isqrt(N))
        c["trail"] = max(1, -(-c["steps_budget"] // c["workers"]))
        c["cap_ir"] = int(c["steps_budget"] * (P.MUL_CHARGE + P.LOOKUP + c["rec_write"]) +
                          4 * (2 * c["trail"] * P.MUL_CHARGE) + P.VERIFY_RESERVE)
        cfgs[name] = c
        all_jobs += [(c, i, name) for i in range(nfix)]
        meta.append({"name": name, "q": c["q"], "rp": c["rp"], "N": N,
                     "depth": c["depth"], "workers": c["workers"],
                     "trail": c["trail"], "steps_budget": c["steps_budget"],
                     "ratio": c["ratio"], "evict": c["evict"], "fixtures": nfix,
                     "cap_ir": c["cap_ir"]})

    nproc = min(20, mp.cpu_count())
    ctx = mp.get_context("fork")
    res = []
    with ctx.Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(_work, all_jobs, chunksize=2)):
            res.append(r)
            if (i + 1) % 64 == 0:
                print("  %d/%d fixtures" % (i + 1, len(all_jobs)), flush=True)

    out = {"mode": "census28-shipped-config",
           "protocol": "protocol28.run_one: meter-first, depth and worker count from the "
                       "config, every record lookup/insert a call into the shipped "
                       "scripts/ir28/dict28.o (evicting) or dict28_nolru.o (control) via "
                       "ctypes; key built by the shipped state_key48w on twelve 4-word "
                       "(255-bit) field lanes; 32-squaring abort implemented",
           "worker_fidelity": {
               "production_trail": "S/1024 = 9.36e34 steps",
               "fixture_trail_formula": "ceil(ratio*sqrt(N)/W)",
               "max_feasible_W_here": 128,
               "why": "trail >= guard ring (6) is required for the guard statistics to be "
                      "at stationarity; at ratio 1.4493 the fixture admits W <= 8 before "
                      "trail < ring.  W enters only through trail; the partition series "
                      "D/E/F measures p_hat against W at FIXED S, which is the only "
                      "W-dependence the fixture can carry.",
               "production_1024_not_reachable": True},
           "eviction_loss": "measured as the A vs G contrast at the same config; G loads "
                            "dict28_nolru.o (drop on four full probes) and therefore loses "
                            "records differently, so the difference is the measured "
                            "eviction cost, not an assumption",
           "procs": nproc, "scales": []}
    for m in meta:
        rs = [r for r in res if r["cfg"] == m["name"]]
        ok = sum(1 for r in rs if r["solved"])
        steps = sum(r["steps"] for r in rs)
        ev = sum(r["c_calls"]["d28_evictions"] for r in rs)
        dr = sum(r["c_calls"]["d28_drops"] for r in rs)
        row = dict(m)
        row.update({
            "success": ok, "rate": ok / m["fixtures"],
            "cp95_lower": P.cp95_lower(ok, m["fixtures"]),
            "steps_executed": steps,
            "records": sum(r["records"] for r in rs),
            "matches": sum(r["matches"] for r in rs),
            "verifies": sum(r["verifies"] for r in rs),
            "aborts": sum(r["aborts"] for r in rs),
            "abandoned": sum(1 for r in rs if r["abandoned"]),
            "exhausted": sum(1 for r in rs if r["exhausted"]),
            "evicted_unmatched": sum(r["evicted_unmatched"] for r in rs),
            "dead_hits": sum(r["dead_hits"] for r in rs),
            "mean_ir_used": sum(r["ir_used"] for r in rs) // m["fixtures"],
            "dict_calls": {k: sum(r["c_calls"][k] for r in rs) for k in rs[0]["c_calls"]},
            "evictions": ev, "drops": dr,
            "tab_slots": rs[0]["tab_slots"],
            "rows_ir": {k: sum(r["rows"][k] for r in rs) for k in rs[0]["rows"]},
        })
        out["scales"].append(row)
        print("scale %-28s ok %3d/%3d cp95 %.4f ev %6d steps %7d" %
              (m["name"], ok, m["fixtures"], row["cp95_lower"], ev, steps), flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    # call-counter assertion: the shipped object really ran
    for s in out["scales"]:
        assert s["dict_calls"]["d28_lookups"] > 0 and s["dict_calls"]["d28_inserts"] > 0, s["name"]
    json.dump(out, open(os.path.join(EV, "census28.json"), "w"), indent=1)
    print("wrote census28.json in %ss" % out["seconds"])


if __name__ == "__main__":
    main()
