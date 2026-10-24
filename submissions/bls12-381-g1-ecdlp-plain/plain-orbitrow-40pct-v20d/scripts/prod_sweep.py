#!/usr/bin/env python3
"""prod_sweep.py -- v20d: the v20c pipeline plus the LEDGER-CHARGED DRAIN ARM.

v20c replayed the production pipeline on the algebra and measured success 0.360
against merge 0.418 at d = 5.  The deficit is the REPORT LAG: a merge inside the
last DP spacing of the budget has no report before the budget expires.  The
ledger prices that lag as work (`drain = workers*2*DP`, i.e. 2*DP per walker,
plus the one-spacing `+2^60` in phase 1) and as the tiny `transport = 2*DP/S_OPS`
success debit.

WHAT THE DRAIN ARM IS, AND WHAT IT REFUSES TO BE
------------------------------------------------
The queue is real: after the budget, every walker with a PENDING report (its
trial has an audited in-budget merge but no in-budget duplicate report) is
stepped to its next DP and the report is delivered to the same pool; recovery
still has to solve x = (a1-a2)/(c2-c1) and verify against the trial's x.

Two things are deliberately NOT done, because at toy scale they would fabricate
success rather than measure it:

  * the drain is not a free extension of the search.  A duplicate pair is only
    counted on a trial that has an INDEPENDENT in-budget merge (the audit over
    visited classes).  Drain-only collisions -- two lanes that first meet inside
    the drain -- are counted and reported as `drain_only_excluded`, never as
    success.  At d = 11 the toy's drain window is ~11x its budget, so those
    would otherwise drive success to ~1.0; in production the drain is
    2*DP/S_OPS = 2^-45 of the search and they are negligible, so excluding them
    is the conservative direction.
  * only pending reports are flushed.  Walking every lane to its next report
    would add fresh search in proportion to the toy's drain/budget ratio; the
    production drain flushes 2^20 workers whose budget is 2^46 spacings, so it
    only ever delivers reports already earned.

The flush cap is `--flush-cap-spacings` (default 8) spacings.  The *charged*
drain is 2 spacings per lane; the cap is larger only so that the toy's coarse DP
grid (budget = 0.18..11.7 spacings per walker, against 2^46 in production) does
not truncate the queue.  `success_budget` (v20c's number, both records inside
the budget) is reported alongside, so the reader can see exactly how much of the
success is drain-delivered.  The production lag debit is 2*DP/S_OPS = 2^-65 in
the ledger's own success model, i.e. below the fourth decimal of 0.399.

CP-95 lower is the Decimal-exact Clopper-Pearson copied from gen_cert_v15k.py.
Protocol is otherwise v20c (qwalk.Universe, invariant a + c*x == canon).
"""
import argparse
import glob
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np          # noqa: E402
import qwalk                # noqa: E402

LAMF = math.log(1 / 0.6097)
SQRTS = 23700.0             # class universe N = SQRTS^2 (rival scale)
R_TAB = 1 << 16
KAPPA = 1.015
W = 64                      # walkers per trial
DRAIN_SPACINGS = 2.0        # ledger charge: drain = workers * 2 * DP
FLUSH_CAP_SPACINGS = 8.0    # queue cap for the toy's coarse DP grid
getcontext().prec = 60


def cp95_lower(k, n, alpha=Decimal('0.05')):
    """Exact one-sided Clopper-Pearson lower bound (copied from gen_cert_v15k.py).

    Smallest p with P[X >= k | p, n] >= alpha, by bisection on the exact Decimal
    binomial tail.  Bracket [0, 1]; gen_cert's [0.30, 0.45] is only valid for its
    own census, the arithmetic is identical.
    """
    k = int(k)
    n = int(n)
    if k == 0:
        return Decimal(0)
    if k >= n:
        return alpha ** (Decimal(1) / Decimal(n))
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p):
        p = Decimal(p)
        q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal('1e-40') * s:
                break
        return s

    lo, hi = Decimal(0), Decimal(1)
    for _ in range(90):
        mid = (lo + hi) / 2
        if mid <= 0 or mid >= 1:
            break
        if tail(mid) > alpha:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal('1e-14'):
            break
    return lo


def shard(job):
    d, bt, seed, kappa, w, flush_cap_sp = job
    U = qwalk.Universe(int(round(6 * SQRTS * SQRTS)), R_TAB, seed=20250)
    q = U.q
    N = U.classes
    s_ref = math.sqrt(2 * N * LAMF)
    budget = int(math.ceil(kappa * s_ref))
    per = budget // w
    lanes = bt * w
    rng = np.random.default_rng(seed & 0x7FFFFFFF)
    xt = rng.integers(1, q - 1, size=bt).astype(np.uint64)
    xv = np.repeat(xt, w).astype(np.uint64)
    a = rng.integers(0, q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, q, size=lanes).astype(np.uint64)
    t0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                   for ai, ci, xi in zip(a.tolist(), cc.tolist(), xv.tolist())],
                  dtype=np.uint64)
    canon, e, _h = U.canon(t0)
    a = (e * a) % q
    cc = (e * cc) % q

    trace = np.empty((per, lanes), dtype=np.uint64)
    rl, rcls, ra, rcc, rw = [], [], [], [], []
    for r in range(per):
        trace[r] = canon
        canon, a, cc, hit = qwalk.step(U, canon, a, cc, xv, d, want_dp=True)
        if hit.any():
            idx = np.flatnonzero(hit)
            rl.append(np.full(idx.size, r, dtype=np.int64))
            rcls.append(canon[idx])
            ra.append(a[idx])
            rcc.append(cc[idx])
            rw.append(idx)
            a0 = rng.integers(0, q, size=idx.size).astype(np.uint64)
            c0 = rng.integers(0, q, size=idx.size).astype(np.uint64)
            s0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                           for ai, ci, xi in zip(a0.tolist(), c0.tolist(),
                                                 xv[idx].tolist())], dtype=np.uint64)
            cn2, e2, _h2 = U.canon(s0)
            canon[idx] = cn2
            a[idx] = (e2 * a0) % q
            cc[idx] = (e2 * c0) % q

    ctr = np.repeat(np.arange(bt, dtype=np.int64), w)
    wid = np.tile(np.arange(lanes, dtype=np.int64), per)

    # ---- independent merge audit: same class, two different walkers, in budget
    key = (ctr[None, :] * np.int64(q + 1) + trace.astype(np.int64)).reshape(-1)
    lvl = np.repeat(np.arange(per, dtype=np.int64)[:, None], lanes,
                    axis=1).reshape(-1)
    if rl:
        rtr = np.repeat(np.arange(bt, dtype=np.int64), w)[np.concatenate(rw)] // w
        key = np.concatenate([key, rtr * np.int64(q + 1)
                              + np.concatenate(rcls).astype(np.int64)])
        lvl = np.concatenate([lvl, np.concatenate(rl)])
        wid = np.concatenate([wid, np.concatenate(rw)])
    o = np.argsort(key, kind="stable")
    sk = key[o]
    dup = np.flatnonzero(sk[1:] == sk[:-1])
    merge = np.zeros(bt, dtype=bool)
    merge_lvl = np.full(bt, 10 ** 9, dtype=np.int64)
    if dup.size:
        sw = wid[o]
        sl = lvl[o]
        for kk in np.unique(sk[dup]).tolist():
            g0 = int(np.searchsorted(sk, kk, "left"))
            g1 = int(np.searchsorted(sk, kk, "right"))
            first = {}
            for ww, ll in zip(sw[g0:g1].tolist(), sl[g0:g1].tolist()):
                if ww not in first or ll < first[ww]:
                    first[ww] = ll
            if len(first) > 1:
                tr = int(kk // (q + 1))
                merge[tr] = True
                sv = sorted(first.values())
                if sv[1] < merge_lvl[tr]:
                    merge_lvl[tr] = sv[1]

    # ---- pass 1: in-budget detection (v20c statistic) ----------------------
    def detect(rlv, rc, rav, rcv, rwv, only_merged=False):
        succ = np.zeros(bt, dtype=bool)
        leven = np.full(bt, 10 ** 9, dtype=np.int64)
        n_ok = n_b = n_d = 0
        bad = dead = excluded = 0
        if rlv.size == 0:
            return succ, leven, n_ok, n_b, n_d, bad, dead, excluded
        tri = rwv // w
        o2 = np.lexsort((rwv, tri, rc))
        cls_s, tri_s, lvl_s = rc[o2], tri[o2], rlv[o2]
        a_s, c_s, w_s = rav[o2], rcv[o2], rwv[o2]
        m = rlv.size
        i = 0
        while i < m:
            j = i + 1
            while j < m and tri_s[j] == tri_s[i] and cls_s[j] == cls_s[i]:
                j += 1
            if j - i >= 2:
                tr = int(tri_s[i])
                k1 = i
                k2 = i + 1
                while k2 < j and w_s[k2] == w_s[k1]:
                    k2 += 1
                if k2 < j:
                    if only_merged and not merge[tr]:
                        excluded += 1
                        i = j
                        continue
                    a1, c1 = int(a_s[k1]), int(c_s[k1])
                    a2, c2 = int(a_s[k2]), int(c_s[k2])
                    L = max(int(lvl_s[k1]), int(lvl_s[k2]))
                    den = (c2 - c1) % q
                    if den:
                        xr = (a1 - a2) * pow(den, -1, q) % q
                        if xr == int(xt[tr]):
                            n_ok += 1
                            if L < per:
                                n_b += 1
                            else:
                                n_d += 1
                            succ[tr] = True
                            if L < leven[tr]:
                                leven[tr] = L
                        else:
                            bad += 1
                    else:
                        dead += 1
            i = j
        return succ, leven, n_ok, n_b, n_d, bad, dead, excluded

    if rl:
        rlv_b = np.concatenate(rl)
        rc_b = np.concatenate(rcls)
        ra_b = np.concatenate(ra)
        rcc_b = np.concatenate(rcc)
        rw_b = np.concatenate(rw)
    else:
        rlv_b = np.zeros(0, dtype=np.int64)
        rc_b = np.zeros(0, dtype=np.uint64)
        ra_b = np.zeros(0, dtype=np.uint64)
        rcc_b = np.zeros(0, dtype=np.uint64)
        rw_b = np.zeros(0, dtype=np.int64)
    succ_b, lev_b, ok_b, nb_b, _nd_b, bad_b, dead_b, _ex_b = detect(
        rlv_b, rc_b, ra_b, rcc_b, rw_b)

    # ---- pass 2: the DRAIN.  Flush only trials with a pending report --------
    need = merge & ~succ_b
    drain_cap = int(round(flush_cap_sp * (1 << d)))
    dl, dcls, da, dcc, dw = [], [], [], [], []
    drain_steps = 0
    all_lanes = np.arange(lanes, dtype=np.int64)
    if need.any():
        sel = np.flatnonzero(need[ctr])
        c0 = canon[sel].copy()
        aa0 = a[sel].copy()
        cc0 = cc[sel].copy()
        x0 = xv[sel]
        live = np.arange(sel.size, dtype=np.int64)
        for r in range(drain_cap):
            if live.size == 0:
                break
            drain_steps += live.size
            cn, ax, cx, hit = qwalk.step(U, c0[live], aa0[live], cc0[live], x0[live],
                                         d, want_dp=True)
            c0[live] = cn
            aa0[live] = ax
            cc0[live] = cx
            if hit.any():
                jj = live[hit]
                dl.append(np.full(jj.size, per + r, dtype=np.int64))
                dcls.append(c0[jj])
                da.append(aa0[jj])
                dcc.append(cc0[jj])
                dw.append(sel[jj])
                live = live[~hit]
        drain_reports = int(sum(x.size for x in dl))
    else:
        drain_reports = 0

    # ---- unified detection over budget + drain, counting only audited merges
    rlv = np.concatenate([rlv_b] + dl) if dl else rlv_b
    rc = np.concatenate([rc_b] + dcls)
    rav = np.concatenate([ra_b] + da)
    rcv = np.concatenate([rcc_b] + dcc)
    rwv = np.concatenate([rw_b] + dw)
    succ, lev, ok, nb, nd, bad, dead, excluded = detect(
        rlv, rc, rav, rcv, rwv, only_merged=True)
    success_d = succ & ~succ_b

    swm = int(np.count_nonzero(succ & ~merge))
    both = merge & succ
    lag = (lev[both] - merge_lvl[both]).astype(float) if both.any() else np.zeros(0)
    return {
        "dp_bits": d, "trials": bt, "seed": seed, "walks": w, "per_walker": per,
        "budget_steps": budget, "s_ref": s_ref, "classes": int(N), "q": int(q),
        "successes": int(succ.sum()),
        "successes_budget": int(succ_b.sum()),
        "successes_drain": int(success_d.sum()),
        "merges": int(merge.sum()),
        "records": int(rlv.size), "records_budget": int(rlv_b.size),
        "records_drain": int(sum(x.size for x in dl)),
        "records_expected": lanes * per / (2.0 ** d),
        "detect_verified": ok, "detect_budget": nb, "detect_drain": nd,
        "denom_dead": dead, "solve_mismatch": bad,
        "drain_only_excluded": excluded,
        "success_without_merge": swm,
        "flush_trials": int(need.sum()), "flush_cap": drain_cap,
        "drain_lane_steps": drain_steps,
        "drain_lane_steps_charged": int(need.sum()) * w * drain_cap,
        "lag_mean": float(lag.mean()) if lag.size else None,
        "lag_max": float(lag.max()) if lag.size else None,
        "merge_mean_over_sref": float(merge_lvl[merge].mean() / s_ref) if merge.any() else None,
    }


def reduce_all(adir, out_path, tag="prod-sweep"):
    rows = {}
    for p in sorted(glob.glob(os.path.join(adir, tag + "-*.json"))):
        with open(p) as f:
            r = json.load(f)
        d = r["dp_bits"]
        a = rows.setdefault(d, dict(trials=0, successes=0, succ_b=0, succ_d=0, merges=0,
                                    records=0, rec_b=0, rec_d=0, detect_verified=0,
                                    det_b=0, det_d=0, denom_dead=0, mismatch=0, swm=0,
                                    excluded=0, flush_trials=0, drain_steps=0,
                                    drain_charged=0, lag_sum=0.0, lag_n=0, lag_max=0.0,
                                    ml_sum=0.0, ml_n=0, s_ref=r["s_ref"],
                                    classes=r["classes"], q=r["q"], walks=r["walks"],
                                    per_walker=r["per_walker"],
                                    budget_steps=r["budget_steps"],
                                    flush_cap=r.get("flush_cap"), shards=0))
        a["trials"] += r["trials"]
        a["successes"] += r["successes"]
        a["succ_b"] += r.get("successes_budget", 0)
        a["succ_d"] += r.get("successes_drain", 0)
        a["merges"] += r["merges"]
        a["records"] += r["records"]
        a["rec_b"] += r.get("records_budget", 0)
        a["rec_d"] += r.get("records_drain", 0)
        a["detect_verified"] += r["detect_verified"]
        a["det_b"] += r.get("detect_budget", 0)
        a["det_d"] += r.get("detect_drain", 0)
        a["denom_dead"] += r["denom_dead"]
        a["mismatch"] += r.get("solve_mismatch", 0)
        a["swm"] += r["success_without_merge"]
        a["excluded"] += r.get("drain_only_excluded", 0)
        a["flush_trials"] += r.get("flush_trials", 0)
        a["drain_steps"] += r.get("drain_lane_steps", 0)
        a["drain_charged"] += r.get("drain_lane_steps_charged", 0)
        if r["lag_mean"] is not None:
            a["lag_sum"] += r["lag_mean"] * r["successes"]
            a["lag_n"] += r["successes"]
        a["lag_max"] = max(a["lag_max"], r["lag_max"] or 0.0)
        if r["merge_mean_over_sref"] is not None:
            a["ml_sum"] += r["merge_mean_over_sref"] * r["merges"]
            a["ml_n"] += r["merges"]
        a["shards"] += 1
    out = []
    for d in sorted(rows):
        a = rows[d]
        n = a["trials"]
        k = a["successes"]
        km = a["merges"]
        cp = cp95_lower(k, n)
        out.append({
            "dp_bits": d, "trials": n, "successes": k, "success": k / n,
            "cp95_lower": float(cp), "cp95_lower_exact": str(cp),
            "successes_budget": a["succ_b"], "successes_drain": a["succ_d"],
            "success_budget": a["succ_b"] / n, "success_drain": a["succ_d"] / n,
            "merges": km, "merge_rate": km / n,
            "merge_cp95_lower": float(cp95_lower(km, n)),
            "records_mean": a["records"] / n,
            "records_budget_mean": a["rec_b"] / n,
            "records_drain_mean": a["rec_d"] / n,
            "records_expected_mean": a["per_walker"] * a["walks"] / (2.0 ** d),
            "detect_verified": a["detect_verified"],
            "detect_budget": a["det_b"], "detect_drain": a["det_d"],
            "denom_dead": a["denom_dead"], "solve_mismatch": a["mismatch"],
            "drain_only_excluded": a["excluded"],
            "success_without_merge": a["swm"], "coherence_ok": a["swm"] == 0,
            "flush_trials_mean": a["flush_trials"] / n,
            "lag_mean_steps": (a["lag_sum"] / a["lag_n"]) if a["lag_n"] else None,
            "lag_max_steps": a["lag_max"] if a["lag_n"] else None,
            "lag_over_spacing": ((a["lag_sum"] / a["lag_n"]) / (2 ** d)) if a["lag_n"] else None,
            "merge_mean_over_sref": (a["ml_sum"] / a["ml_n"]) if a["ml_n"] else None,
            "budget_spacings_per_walker": a["per_walker"] / (2.0 ** d),
            "flush_cap_steps": a["flush_cap"],
            "flush_cap_spacings": a["flush_cap"] / (2.0 ** d) if a["flush_cap"] else None,
            "drain_lane_steps_mean": a["drain_steps"] / n,
            "drain_multiplicity": (a["drain_steps"] / a["drain_charged"]
                                   if a["drain_charged"] else None),
            "bayes_model": 1 - math.exp(-(KAPPA ** 2) * LAMF),
            "s_ref": a["s_ref"], "classes": a["classes"], "q": a["q"],
            "walks": a["walks"], "per_walker": a["per_walker"],
            "budget_steps": a["budget_steps"], "shards": a["shards"]})
    gate = {"rule": "CP-95 lower >= 0.39 at every density",
            "min_cp95_lower": min(r["cp95_lower"] for r in out) if out else None,
            "min_at_dp_bits": min(out, key=lambda r: r["cp95_lower"])["dp_bits"] if out else None,
            "coherence_all": all(r["coherence_ok"] for r in out),
            "success_le_merge_all": all(r["successes"] <= r["merges"] for r in out),
            "passed": bool(out) and all(r["cp95_lower"] >= 0.39 for r in out)}
    res = {"package": "plain-orbitrow-40pct-v20d", "arm": "prod-sweep+drain",
           "protocol": {
               "walk": "qwalk quotient walk, invariant a + c*x == canon",
               "dp": "theta=2^-d gates REPORTING; walker restarts at a DP",
               "pool": "records only; in-budget records = budget/2^d",
               "drain": "after the budget each pending walker is flushed to its next "
                        "report (cap 8 spacings); only audited in-budget merges may "
                        "be recovered, drain-only collisions are excluded",
               "detect": "duplicate class among pool records, two walkers",
               "recovery": "(a1-a2)/(c2-c1) mod q from the two records, verified",
               "merge_stat": "independent audit over in-budget visited classes",
               "budget": "kappa*S_ref total steps; lag charged separately, as ledger"},
           "rows": out, "gate": gate}
    with open(out_path, "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trials", type=int, default=8000)
    ap.add_argument("--chunk", type=int, default=100)
    ap.add_argument("--dens", default="5,7,9,11")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--kappa", type=float, default=KAPPA)
    ap.add_argument("--walks", type=int, default=W)
    ap.add_argument("--flush-cap-spacings", type=float, default=FLUSH_CAP_SPACINGS)
    ap.add_argument("--seed", type=int, default=20270)
    ap.add_argument("--raw", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--reduce", action="store_true")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    adir = os.path.abspath(a.raw or os.path.join(here, "..", "evidence", "raw"))
    os.makedirs(adir, exist_ok=True)
    out_path = a.out or os.path.join(here, "..", "evidence", "prod-sweep.json")
    if a.reduce:
        res = reduce_all(adir, out_path)
        print(json.dumps(res["gate"], indent=1))
        for r in res["rows"]:
            print("d=%2d n=%6d succ %.4f CP95 %.4f (b %.4f + d %.4f) merge %.4f "
                  "rec %.0f ex %d swm %d"
                  % (r["dp_bits"], r["trials"], r["success"], r["cp95_lower"],
                     r["success_budget"], r["success_drain"], r["merge_rate"],
                     r["records_mean"], r["drain_only_excluded"],
                     r["success_without_merge"]))
        return
    dens = [int(t) for t in a.dens.split(",") if t]
    jobs = []
    for d in dens:
        done = 0
        i = 0
        while done < a.trials:
            take = min(a.chunk, a.trials - done)
            jobs.append((d, take, a.seed * 100000 + d * 1000 + i, a.kappa, a.walks,
                         a.flush_cap_spacings))
            done += take
            i += 1
    t0 = time.time()
    print("jobs", len(jobs), "dens", dens, "trials/density", a.trials,
          "flush_cap_spacings", a.flush_cap_spacings, flush=True)
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(shard, jobs):
            r["succ_rate"] = r["successes"] / r["trials"]
            name = "prod-sweep-d%02d-s%05d.json" % (r["dp_bits"], r["seed"] % 100000)
            with open(os.path.join(adir, name), "w") as f:
                json.dump(r, f, indent=1, sort_keys=True)
            print("  d=%2d n=%4d succ %.4f (b %.4f d %.4f) merge %.4f rec %6d/%5d "
                  "fl %4d ex %3d swm %d %.1fs"
                  % (r["dp_bits"], r["trials"], r["succ_rate"],
                     r["successes_budget"] / r["trials"],
                     r["successes_drain"] / r["trials"], r["merges"] / r["trials"],
                     r["records_budget"], r["records_drain"], r["flush_trials"],
                     r["drain_only_excluded"], r["success_without_merge"],
                     time.time() - t0), flush=True)
    res = reduce_all(adir, out_path)
    print(json.dumps(res["gate"], indent=1))
    for r in res["rows"]:
        print("d=%2d n=%6d succ %.4f CP95 %.4f (b %.4f + d %.4f) merge %.4f rec %.0f "
              "ex %d swm %d"
              % (r["dp_bits"], r["trials"], r["success"], r["cp95_lower"],
                 r["success_budget"], r["success_drain"], r["merge_rate"],
                 r["records_mean"], r["drain_only_excluded"],
                 r["success_without_merge"]))
    print("wrote", out_path, "%.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
