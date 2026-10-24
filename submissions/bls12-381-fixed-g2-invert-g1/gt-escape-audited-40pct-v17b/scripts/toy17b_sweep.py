#!/usr/bin/env python3
"""v17b extended design-budget census (D5 discipline: machine-generated only).

Runs the PRODUCTION toy step/slot rule at the DESIGN budget S = s_ops(N,R) only,
T=1 trail, DP density 1 (every state recorded), rule="new" -- the exact config of
v17's main_a..main_d shards -- with fresh disjoint seed ranges, sharded over
worker processes.  Nothing is run at 4x the budget.

Config equality is enforced by toy17.py's own recorded fields (q,R,dpbits,trails,
rule,use_sqr,S_design,dp_shift,idx_bits); toy17b_merge.py refuses to pool shards
whose config dict differs from the main_a..main_d shards' config, so the old-rule
control and the S=225 side-run can never leak into this pool.

usage: python3 toy17b_sweep.py            (writes evidence/toy17b_shards/*.json)
"""
import json
import os
import subprocess
import sys

ROOT = __file__.rsplit("/", 2)[0]
SHARD_DIR = ROOT + "/evidence/toy17b_shards"
PY = sys.executable

WORKERS = 20
TRIALS_PER_WORKER = 4000          # 20 x 4000 = 80000 fresh design-budget trials
SEED_BASE = 200001                # disjoint from v17 shards (seeds 1..1600, 5001.., 7001..)
SEED_STRIDE = 10007               # > trials-per-worker, so seed ranges never overlap

CFG_MAIN = dict(q=523, R=128, dpbits=0, trails=1, rule="new", use_sqr=True)


def main():
    os.makedirs(SHARD_DIR, exist_ok=True)
    only = sys.argv[1:]
    procs = []
    for w in range(WORKERS):
        name = "ext_%02d" % w
        if only and name not in only:
            continue
        out = SHARD_DIR + "/" + name + ".json"
        if os.path.exists(out):
            print("exists", name)
            continue
        cfg = dict(CFG_MAIN, trials=TRIALS_PER_WORKER, seed=SEED_BASE + w * SEED_STRIDE)
        p = subprocess.Popen(["nice", "-n", "10", PY, ROOT + "/scripts/toy17.py"],
                             stdin=subprocess.PIPE, stdout=open(out + ".tmp", "w"))
        p.stdin.write(json.dumps(cfg).encode())
        p.stdin.close()
        procs.append((name, out, p))
    for name, out, p in procs:
        p.wait()
        os.replace(out + ".tmp", out)
        d = json.load(open(out))
        print(name, "succ", d["success_num"], "/", d["success_den"],
              "degen", d["degenerate"], "max_fwd", d["max_fwd"], flush=True)
    print("ALL DONE", len(procs), flush=True)


if __name__ == "__main__":
    main()
