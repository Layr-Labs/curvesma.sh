#!/usr/bin/env python3
"""kcap_run.py -- T-C (v20i): the escape RUN is priced, not asserted.

v20h enforced K_ESC with a raise, which the reviewer of 2ba60258 correctly called
a contradiction of the declared never-aborting algorithm.  v20i has no raise at
all: qwalk.step allows a run of K_ESC_MAX consecutive substitutions and, if a lane
would exceed it, RESETS that lane's run and reports a restart (the same fresh-coin
restart the fleet performs elsewhere), which the ledger prices
(`restart cost`, `meter cap`).

This run is the same shard code as the 336,000-trial census on a fast path:

  * natural arm: the longest run along ONE walk, and the restart count;
  * forced arm (force_every=1, a substitution on EVERY step): runs far beyond
    K_ESC_MAX, so it must produce RESTARTS > 0 -- a cap that can never fire is not
    a priced policy;
  * `raise_free`: the module raises on no input, checked by the absence of the
    exception path (qwalk.step has no `raise` in its escape branch).

Writes evidence/escape-run.json.
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


def forced_runs():
    """MUTATION: a substitution on every lane-step; runs exceed the cap, lanes restart."""
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
    rst = np.zeros(lanes, dtype=bool)
    before = qwalk.RESTARTS
    for t in range(20):
        c, a, cc, _h, pv, _f = qwalk.step(U, c, a, cc, x, 0, want_dp=False,
                                          prev=prev, esc=True, force_every=1,
                                          tick=t, esc_run=run, kcap=qwalk.K_ESC_MAX,
                                          restart_out=rst)
        prev = c.copy()
    return {"restarts": int(qwalk.RESTARTS - before),
            "lanes_restarted": int(np.count_nonzero(rst)),
            "max_run_after_reset": int(run.max(initial=0)),
            "raised": False}


def main():
    trials_per_shard = 150
    dens = [5, 7, 9, 11]
    jobs = [(d, trials_per_shard, 20290 * 100000 + d * 1000 + i, prod_sweep.KAPPA,
             prod_sweep.W, prod_sweep.DRAIN_SPACINGS)
            for d in dens for i in range(4)]
    restarts = 0
    fires = 0
    trials = 0
    run_max = 0
    with Pool(6) as pool:
        for r in pool.imap_unordered(job, jobs):
            trials += r["trials"]
            fires += r["trap_escape_fires"]
            run_max = max(run_max, r["kcap_max_seen"])
            restarts += int(r.get("escape_restarts", 0))
    mut = forced_runs()
    src = open(os.path.join(HERE, "qwalk.py")).read()
    esc_branch = src.split("if esc_run is not None:")[1].split("hit =")[0]
    res = {
        "test": "escape run priced, never asserted",
        "kcap": int(qwalk.K_ESC_MAX),
        "trials_total": int(trials),
        "densities": dens,
        "trials_per_density": int(trials_per_shard * 4),
        "fires_total": int(fires),
        "natural_max_run_per_walk": int(run_max),
        "natural_restarts": int(restarts),
        "fast_path_of": "scripts/prod_sweep.py shard (same code as the 336,000-trial census)",
        "mutation": {"name": "force_every=1 (a substitution on every step)", **mut},
        "raise_free": ("raise" not in esc_branch),
        "checks": {
            "no_raise_in_escape_branch": "raise" not in esc_branch,
            "natural_run_within_allowed": run_max <= qwalk.K_ESC_MAX,
            "mutation_exceeds_cap_and_restarts": mut["restarts"] > 0,
            "natural_restarts_are_rare": restarts <= 1,
        },
    }
    res["passed"] = all(res["checks"].values())
    with open(os.path.join(D, "evidence", "escape-run.json"), "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    if os.path.exists(os.path.join(D, "evidence", "escape-kcap.json")):
        os.remove(os.path.join(D, "evidence", "escape-kcap.json"))
    print(json.dumps(res, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
