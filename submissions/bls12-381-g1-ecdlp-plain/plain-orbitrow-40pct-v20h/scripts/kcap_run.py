#!/usr/bin/env python3
"""kcap_run.py -- T-C: run the escape census with K_ESC ENFORCED in the walker loop.

prod_sweep.shard now passes the per-lane run-length array to qwalk.step with
kcap = K_ESC, so a walk that substituted more than K_ESC times in a row RAISES
instead of finishing.  This is the assertion run: the same shard code as the
336,000-trial census, on a fast path (a few thousand trials), proving the cap
never trips on the natural arm and that a walk which exceeds it is caught.

The mutation control is a forced arm at force_every=1 (a substitution on EVERY
step, which necessarily exceeds the cap): it must RAISE.  A cap that cannot fire
is not an assertion.

Writes evidence/escape-kcap.json.
"""
import json
import os
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, D)
import numpy as np            # noqa: E402
import prod_sweep             # noqa: E402
import qwalk                  # noqa: E402


def job(args):
    return prod_sweep.shard(args)


def forced_trips():
    """MUTATION: force a substitution on every lane-step; the cap must raise."""
    U = qwalk.Universe((1 << 28) + 1, 1 << 10, seed=20250)
    q = U.q
    lanes = 8
    rng = np.random.default_rng(1)
    x = rng.integers(1, q - 1, size=lanes).astype(np.uint64)
    a = rng.integers(0, q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, q, size=lanes).astype(np.uint64)
    c = U.canon((a + cc * x) % q)[0]
    prev = c.copy()
    run = np.zeros(lanes, dtype=np.int64)
    try:
        for t in range(20):
            c, a, cc, _h, pv, _f = qwalk.step(U, c, a, cc, x, 0, want_dp=False,
                                              prev=prev, esc=True, force_every=1,
                                              tick=t, esc_run=run, kcap=qwalk.K_ESC)
            prev = c.copy()
    except qwalk.K_ESCExceeded as e:
        return {"raised": True, "message": str(e), "max_run": int(run.max())}
    return {"raised": False, "message": "", "max_run": int(run.max())}


def main():
    trials_per_shard = 150
    dens = [5, 7, 9, 11]
    jobs = [(d, trials_per_shard, 20290 * 100000 + d * 1000 + i, prod_sweep.KAPPA,
             prod_sweep.W, prod_sweep.DRAIN_SPACINGS)
            for d in dens for i in range(4)]
    trips = []
    fires = 0
    trials = 0
    kcap_max = 0
    with Pool(6) as pool:
        for r in pool.imap_unordered(job, jobs):
            trials += r["trials"]
            fires += r["trap_escape_fires"]
            kcap_max = max(kcap_max, r["kcap_max_seen"])
            if r.get("kcap_tripped"):
                trips.append((r["dp_bits"], r["kcap_max_seen"]))
    mut = forced_trips()
    res = {
        "test": "K_ESC enforced in the walker loop",
        "kcap": int(qwalk.K_ESC),
        "trials_total": int(trials),
        "densities": dens,
        "trials_per_density": int(trials_per_shard * 4),
        "fires_total": int(fires),
        "kcap_max_seen": int(kcap_max),
        "trips": trips,
        "fast_path_of": "scripts/prod_sweep.py shard (same code as the 336,000-trial census)",
        "mutation": {"name": "force_every=1 (a substitution on every step)",
                     **mut},
        "checks": {
            "cap_never_tripped_on_natural_arm": not trips,
            "observed_run_within_cap": kcap_max <= qwalk.K_ESC,
            "mutation_raises": bool(mut["raised"]),
        },
    }
    res["passed"] = all(res["checks"].values())
    with open(os.path.join(D, "evidence", "escape-kcap.json"), "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
