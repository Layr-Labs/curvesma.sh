#!/usr/bin/env python3
"""v17b design-budget census merge -- machine-generated, exact integers only.

Pools the v17 main_a..main_d shards (design budget, T=1, DP density 1, rule=new)
with the v17b ext_* shards (same config, fresh disjoint seeds) and refuses to pool
anything whose recorded configuration differs.  Emits, from the parsed JSON only:

  * success_num / success_den  (exact integer counts)
  * the one-sided 95% Clopper-Pearson LOWER bound (exact bisection on the exact
    binomial tail computed with log-gamma; no scipy)
  * the two-sided 95% CP interval and the normal-approx +-2SE, for disclosure
  * per-shard sha256 digests

usage: python3 toy17b_merge.py > ../evidence/toy17b_census.json
"""
import glob
import hashlib
import json
import math
import os
import sys

ROOT = __file__.rsplit("/", 2)[0]
CONFIG_KEYS = ("mode", "rule", "use_sqr", "q", "rp", "N", "R", "dpbits", "trails",
               "S_design", "dp_shift", "idx_bits", "slot_rule")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def log_binom_row(n):
    lg = math.lgamma
    lnf = [lg(i + 1) for i in range(n + 1)]
    return [lnf[n] - lnf[i] - lnf[n - i] for i in range(n + 1)]


def logsf(k, n, p, logc):
    """log P(X >= k) for X ~ Bin(n,p), by log-sum-exp with the row cached."""
    if k <= 0:
        return 0.0
    lp, lq = math.log(p), math.log1p(-p)
    terms = []
    for i in range(k, n + 1):
        terms.append(logc[i] + i * lp + (n - i) * lq)
    m = max(terms)
    return m + math.log(sum(math.exp(t - m) for t in terms))


def cp_lower(k, n, alpha=0.05, logc=None):
    """p such that P(X >= k; n, p) = alpha, i.e. the 1-alpha CP lower bound."""
    if k == 0:
        return 0.0
    logc = logc or log_binom_row(n)
    lo, hi = 0.0, k / n
    for _ in range(80):
        mid = (lo + hi) / 2
        if logsf(k, n, mid, logc) > math.log(alpha):
            hi = mid            # tail too big for alpha -> p too large
        else:
            lo = mid
    return (lo + hi) / 2


def cp_upper(k, n, alpha=0.05, logc=None):
    """p such that P(X <= k; n, p) = alpha."""
    logc = logc or log_binom_row(n)
    if k >= n:
        return 1.0
    lo, hi = k / n, 1.0
    for _ in range(80):
        mid = (lo + hi) / 2
        s = logsf(k + 1, n, mid, logc)          # P(X >= k+1) = 1 - P(X <= k)
        if s > math.log(alpha):
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def main():
    # POOL = the v17b ext_* shards only: one command, one config, fresh disjoint
    # seeds.  The v17 main_* shards are reported as UN-POOLED corroboration: they
    # were written by the pre-v17b toy17.py and do not record `use_sqr`, so their
    # config tuple cannot be proven equal and they are not allowed into the pool.
    paths = sorted(glob.glob(ROOT + "/evidence/toy17b_shards/ext_*.json"))
    if not paths:
        print(json.dumps({"error": "no shards"}))
        return 1
    shards, cfgs = [], {}
    num = den = deg = 0
    for pth in paths:
        d = json.load(open(pth))
        cfg = tuple(d.get(k) for k in CONFIG_KEYS)
        cfgs.setdefault(cfg, []).append(os.path.basename(pth))
        shards.append({"file": os.path.basename(pth), "sha256": sha(pth),
                       "trials": d["trials"], "success_num": d["success_num"],
                       "success_den": d["success_den"],
                       "degenerate": d["degenerate"], "seed": d.get("seed"),
                       "max_fwd": d.get("max_fwd"), "S_design": d["S_design"],
                       "over_budget": d.get("over_budget")})
        num += d["success_num"]
        den += d["success_den"]
        deg += d["degenerate"]
    if len(cfgs) != 1:
        print(json.dumps({"error": "CONFIG MISMATCH -- refusing to pool",
                          "configs": {str(k): v for k, v in cfgs.items()}}, indent=1))
        return 2
    cfgs_keys = list(cfgs)[0]
    cfg = dict(zip(CONFIG_KEYS, cfgs_keys))
    assert num <= den
    logc = log_binom_row(den)
    lo1 = cp_lower(num, den, 0.05, logc)        # one-sided 95% lower (alpha=0.05)
    lo2 = cp_lower(num, den, 0.025, logc)       # two-sided 95% lower
    hi2 = cp_upper(num, den, 0.025, logc)
    p = num / den
    se = math.sqrt(p * (1 - p) / den)
    corr = {"files": [], "trials": 0, "success_num": 0, "pooled": False}
    for pth in sorted(glob.glob(ROOT + "/evidence/toy17_shards/main_*.json")):
        d = json.load(open(pth))
        k = tuple(d.get(x) for x in CONFIG_KEYS if x != "use_sqr")
        ref = tuple(cfg2 for cfg2, kk in zip(cfgs_keys, CONFIG_KEYS) if kk != "use_sqr")
        if k != ref:
            continue
        corr["files"].append(os.path.basename(pth))
        corr["trials"] += d["trials"]
        corr["success_num"] += d["success_num"]
    corr["reason_not_pooled"] = ("shards predate the v17b toy17.py and record no use_sqr "
                                 "field, so config equality with the pool cannot be proven")
    out = {
        "mode": "census-merge-v17b",
        "corroboration_unpooled": corr,
        "config": cfg,
        "config_shards": {str(k): v for k, v in cfgs.items()},
        "shard_count": len(shards),
        "trials": den, "success_num": num, "success_den": den,
        "success_rate_float": round(p, 6),
        "degenerate_total": deg,
        "cp95_one_sided_lower": round(lo1, 6),
        "cp95_two_sided": [round(lo2, 6), round(hi2, 6)],
        "normal_2se_halfwidth": round(2 * se, 6),
        "normal_95_lower": round(p - 1.96 * se, 6),
        "alpha_one_sided": 0.05,
        "bound_method": "Clopper-Pearson by bisection on the exact binomial tail "
                        "(log-sum-exp with lgamma row); no scipy",
        "shards": shards,
        "command": "python3 scripts/toy17b_merge.py > evidence/toy17b_census.json",
    }
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
