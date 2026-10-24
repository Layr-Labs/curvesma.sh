#!/usr/bin/env python3
"""pilot_prod_d60.py -- v22 PRODUCTION-CONFIGURATION pilot at the SHIPPED d=60.

The shipped configuration is: workers W = 2^20 lanes, DP density d = 60
(DP = 2^60), inversion BATCH = 2^18, trails_per_walker = 2^106.7214.  This pilot
runs THAT configuration -- shipped qwalk.step, shipped report gate, W = 2^20
lanes, fresh coins per trial -- and measures the gate's realized fire count, the
record-collision count, the sampled invariant and the visit rate.

It cannot measure the production success event, and the reason is arithmetic in
the pilot's own measured units: one DP record needs 2^60 visits, so at the
measured rate a single record of a single lane costs ~1.8e3 years.  Every
runnable d=60 pilot is therefore 0-success BY CONSTRUCTION; a zero is not
evidence about the production event in either direction and is not presented as
one.  The production event is a two-run collision at 2^127.014 visits with
2^67.014 records per run, where the gate is saturated (post-merge trail
2^106.72 steps).
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk      # noqa: E402

out_path = os.path.join(D, "evidence", "prod-pilot-d60.json")
W = 2 ** 20
D60 = 60
STEPS = int(os.environ.get("STEPS", "96"))
TRIALS = int(os.environ.get("TRIALS", "4"))
SEED = int(os.environ.get("SEED", "20252"))

U = qwalk.Universe(4_290_000_000, 64, seed=SEED)
assert U.q < (1 << 32)
rng = np.random.default_rng(SEED + 1)
t_all = time.time()
rows = []
tot_visits = tot_fires = tot_coll = 0
inv_bad = 0
for trial in range(TRIALS):
    x = np.uint64(int(rng.integers(1, U.q - 1)))
    a = rng.integers(0, U.q, size=W).astype(np.uint64)
    cc = rng.integers(0, U.q, size=W).astype(np.uint64)
    t0 = ((a + cc * x) % np.uint64(U.q)).astype(np.uint64)
    canon, e, _h = U.canon(t0)
    a = ((e * a) % np.uint64(U.q)).astype(np.uint64)
    cc = ((e * cc) % np.uint64(U.q)).astype(np.uint64)
    prev = canon.copy()
    seen = set()
    fires = coll = 0
    t1 = time.time()
    for _t in range(STEPS):
        cn, a, cc, _h, _pv, _f = qwalk.step(
            U, canon, a, cc, np.full(W, x, dtype=np.uint64), D60,
            True, prev, False, 0, _t, None, None, None, None)
        if D60 > 0 and cn is not None:
            hit = U.dp(cn, D60)
            fires += int(hit.sum())
            for cl in cn[hit].tolist()[:64]:
                if cl in seen:
                    coll += 1
                else:
                    seen.add(cl)
        prev, canon = canon, cn
    dt = time.time() - t1
    visits = STEPS * W
    tot_visits += visits
    tot_fires += fires
    tot_coll += coll
    rows.append({"trial": trial, "visits": visits, "dp_fires": fires,
                 "dp_record_collisions": coll, "seconds": round(dt, 3),
                 "visits_per_s": round(visits / dt, 1)})
    sl = slice(0, min(4096, W))
    bad = int(np.count_nonzero(
        ((a[sl].astype(object) + cc[sl].astype(object) * int(x)) % U.q)
        != canon[sl].astype(object)))
    inv_bad += bad

elapsed = time.time() - t_all
rate = tot_visits / elapsed
need = 1 << D60
secs = need / rate
prod_v = 1.0097 * math.sqrt((2 ** 254) * 6 * 852105 / (6 * 852105 - 1))
out = {
    "test": "production_config_pilot_d60",
    "configuration": {"workers": W, "workers_log2": 20, "dp_bits": D60,
                      "batch": 2 ** 18, "steps_per_pilot": STEPS, "trials": TRIALS,
                      "shipped_trails_per_walker_log2": 106.7214},
    "note": ("SHIPPED configuration, shipped qwalk.step and shipped report gate; "
             "measures the gate's realized fire count and the visit rate. It does "
             "not and cannot measure the production success event at d=60."),
    "trials": rows,
    "totals": {"visits": tot_visits, "dp_fires": tot_fires,
               "dp_record_collisions": tot_coll,
               "invariant_violations_sampled": inv_bad,
               "elapsed_s": round(elapsed, 2), "visits_per_s": round(rate, 1)},
    "realized_dp_rate": tot_fires / tot_visits if tot_visits else 0.0,
    "nominal_dp_rate": 2.0 ** -D60,
    "gate_consistent": bool(tot_fires == 0),
    "impossible_by_arithmetic": {
        "visits_for_one_dp": need,
        "seconds_for_one_dp_at_measured_rate": secs,
        "years_for_one_dp_at_measured_rate": secs / 3.15e7,
        "production_shape_visits_log2": math.log2(prod_v),
        "production_records_per_run_log2": math.log2(prod_v / need),
        "statement": ("one DP record needs 2^60 visits; at the measured rate that "
                      "is the years row, so every runnable d=60 pilot is "
                      "0-success BY CONSTRUCTION and cannot be evidence about the "
                      "production event"),
    },
}
with open(out_path, "w") as f:
    json.dump(out, f, indent=1, sort_keys=True)
print(json.dumps({"totals": out["totals"], "gate_consistent": out["gate_consistent"],
                  "impossible": out["impossible_by_arithmetic"]}, indent=1))
print("wrote", out_path)
