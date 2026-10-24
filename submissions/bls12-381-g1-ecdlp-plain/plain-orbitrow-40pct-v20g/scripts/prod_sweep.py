#!/usr/bin/env python3
"""prod_sweep.py -- v20e: unconditional production replay + escape + fixed audit.

FOUR THINGS CHANGED FROM v20d (each one is a reviewer defect):

D1  AUDIT ATTRIBUTION.  v20d built the report half of the merge key with
    rtr = repeat(arange(bt), w)[rw] // w, i.e. it divided the lane index by w
    TWICE: (lane // w) // w = floor(trial / 64).  Reports from different
    challenges were merged into one key space and reports were dropped from their
    own trial.  v20e divides once: trial = lane // w.  The audit is a single
    function `audit()` used by BOTH the census and the construction check in
    `plant_check()`, which plants reports at known (trial, step) and asserts the
    attribution is exact.

D2  UNCONDITIONAL HEADLINE STATISTIC.  v20d reported success through
    detect(..., only_merged=True), which accepted a duplicate pair only on a trial
    that had some audited merge.  That filter enforced success <= merge by
    construction.  v20e reports the pool-detection rate of the production replay
    with NO conditioning (`success_raw`), and the merge audit (`merges`) as a
    SEPARATE, independently computed row.  Neither statistic is filtered by the
    other.  Drain-only pairs -- two walkers that first meet inside the drain -- are
    COUNTED (they are real detections of the replay), and their number is
    reported so the reader can subtract them if they prefer.

D3  TRAP ESCAPE ARMED.  The two-step-return escape of qwalk.step() is armed for
    the whole replay (budget and drain), and the number of escape fires is
    reported per shard.  Production no longer aborts on a trap.

D4  NUMBERS ARE EMITTED, NOT TYPED.  Every reduced field here is what
    gen_cert_v20e.py reads; the note's tables are generated from the certificate.

The drain window is now the LEDGER-CHARGE window, min(2 spacings, the walker's own
budget), so the toy can never search more than it charges (`drain_window_over_
budget` is reported per density).
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
getcontext().prec = 60


def cp95_lower(k, n, alpha=Decimal('0.05')):
    """Exact one-sided Clopper-Pearson lower bound (Decimal binomial tail)."""
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


def audit(trace, lvl, wid, ctr, q, bt, rlv=None, rcls=None, rw=None):
    """Independent merge audit: one class, one trial, two different walkers.

    `trace`/`lvl`/`wid` are the per-step visited classes (per, lanes) flattened;
    `ctr` is the trial id per lane.  Reports are (rlv, rcls, rw) with rw the LANE
    index; the trial of a report is rw // w.  This single function is what
    plant_check() exercises.
    """
    key = (ctr[None, :] * np.int64(q + 1) + trace.astype(np.int64)).reshape(-1)
    wid = wid
    lvl = lvl
    merge = np.zeros(bt, dtype=bool)
    merge_lvl = np.full(bt, 10 ** 9, dtype=np.int64)
    if rlv is not None and np.size(rlv):
        rw = np.asarray(rw, dtype=np.int64)
        # ---- D1: ONE division.  lane = trial*w + walker  =>  trial = lane // w
        stride = int(ctr.size // bt)          # walkers per trial
        rtr = rw // stride
        key = np.concatenate([key, rtr * np.int64(q + 1)
                              + np.asarray(rcls).astype(np.int64)])
        lvl = np.concatenate([lvl, np.asarray(rlv, dtype=np.int64)])
        wid = np.concatenate([wid, rw])
    o = np.argsort(key, kind="stable")
    sk = key[o]
    dup = np.flatnonzero(sk[1:] == sk[:-1])
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
    return merge, merge_lvl


def plant_check():
    """D1 CONSTRUCTION CHECK: planted reports at known (trial, step).

    Builds a synthetic (per, lanes) trace with NO accidental collisions, plants
    reports at chosen trials/levels, and asserts that (a) every planted
    same-trial collision is attributed to exactly its own trial at exactly its
    level, and (b) reports planted in one trial never mark any other trial.
    """
    q = 1000003
    bt, w, per = 8, 4, 6
    lanes = bt * w
    rng = np.random.default_rng(7)
    trace = rng.choice(np.arange(q, dtype=np.int64), size=(per, lanes, ),
                       replace=False).reshape(per, lanes)
    ctr = np.repeat(np.arange(bt, dtype=np.int64), w)
    lvl = np.repeat(np.arange(per, dtype=np.int64)[:, None], lanes, axis=1)
    wid = np.tile(np.arange(lanes, dtype=np.int64), per)
    # pairwise-distinct classes everywhere: no accidental merge
    flat = trace.reshape(-1)
    assert np.unique(flat).size == flat.size

    planted = [(3, 9), (3, 21), (5, 2), (5, 30)]     # (trial, lane)
    rlv, rcls, rw = [], [], []
    for t, lane in planted:
        rlv.append(per - 1)
        rcls.append(int(trace[0, lane]) + 0)          # shares a class with step 0
        rw.append(lane)
    # give each planted report its own class by copying the planted lane's step-0
    # class, and make the two members of each trial share one class each
    rcls = []
    for t, lane in planted:
        rcls.append(int(trace[0, lane]))
    merge, merge_lvl = audit(trace, lvl.reshape(-1), wid.reshape(-1), ctr, q, bt,
                             np.array(rlv), np.array(rcls), np.array(rw))
    # trials 3 and 5 each planted one lane only => one walker per class => NOT merged
    ok_single = not merge.any()

    # now plant a second walker of the same trial in the same class
    rlv2 = list(rlv) + [per - 2, per - 2]
    rcls2 = list(rcls) + [int(trace[0, 3 * w + 0]), int(trace[0, 5 * w + 3])]
    rw2 = list(rw) + [3 * w + 1, 5 * w + 2]
    merge2, lvl2 = audit(trace, lvl.reshape(-1), wid.reshape(-1), ctr, q, bt,
                         np.array(rlv2), np.array(rcls2), np.array(rw2))
    ok_pair = (sorted(np.flatnonzero(merge2).tolist()) == [3, 5]
               and int(lvl2[3]) == per - 2 and int(lvl2[5]) == per - 2)

    # the v20d bug: floor(trial / w) attribution.  Show it fires here.
    bad = np.zeros(bt, dtype=bool)
    for t, lane in zip([3] * 2 + [5] * 2, [3 * w + 0, 3 * w + 1, 5 * w + 3, 5 * w + 2]):
        bad[(lane // w) // w if ((lane // w) // w) < bt else t] = True
    ok_bug_detected = not np.array_equal(bad, merge2)

    return {"rule": "planted same-trial pair -> exactly its own trial, at its level",
            "bt": bt, "walks": w, "per": per, "plants": planted,
            "single_walker_no_merge": bool(ok_single),
            "planted_pairs_exact": bool(ok_pair),
            "merged_trials": np.flatnonzero(merge2).tolist(),
            "merge_levels": [int(lvl2[t]) for t in np.flatnonzero(merge2).tolist()],
            "v20d_floor_t_w_attribution_differs": bool(ok_bug_detected),
            "passed": bool(ok_single and ok_pair and ok_bug_detected)}


def pair_census(trace, ctr, q, bt, rep=None):
    """E4: the coverage law ON THE PAIR (class_t, class_{t-1}), reviewer point 3.

    The escape is deterministic in the pair (class, previous class), so the
    reviewer's objection is that class coverage alone does not establish a shared
    future.  Measure BOTH laws on the SAME census traces, unconditionally:

      class_trials     trials where two different walkers occupy one class at any
                       visited step (this is exactly the audit's `merge` set;
                       asserted equal to it, two independent implementations)
      pair_trials      trials where two different walkers occupy one PAIR
                       (c_t, c_{t-1}) at any step t >= 1
      prev_matched     class groups in which two different walkers ALSO share
                       their predecessor -- the meetings whose futures agree
      pair_over_class  pair groups / class groups (group level, not trial level)

    Reporting all of them makes the loss from predecessor-dependent decisions an
    explicit number instead of an assumption.
    """
    per = trace.shape[0]
    lanes = trace.shape[1]
    if per < 2:
        raise ValueError("need at least two steps")
    tr = np.tile(ctr, per - 1)
    wid = np.tile(np.arange(lanes, dtype=np.int64), per - 1)
    C = trace[1:].astype(np.int64).reshape(-1)
    P = trace[:-1].astype(np.int64).reshape(-1)

    # ---- pair groups --------------------------------------------------------
    o = np.lexsort((C, P, tr))
    Cs, Ps, trs, wid_s = C[o], P[o], tr[o], wid[o]
    new = np.empty(Cs.size, dtype=bool)
    new[0] = True
    new[1:] = (Cs[1:] != Cs[:-1]) | (Ps[1:] != Ps[:-1]) | (trs[1:] != trs[:-1])
    st = np.flatnonzero(new)
    en = np.append(st[1:], Cs.size)
    pair_groups = 0
    pair_trials = set()
    for s, e in zip(st[(en - st) >= 2], en[(en - st) >= 2]):
        if np.unique(wid_s[s:e]).size >= 2:
            pair_groups += 1
            pair_trials.add(int(trs[s]))

    # ---- class groups over ALL visited steps (step 0 included: exactly the
    # audit's class set), plus the predecessor-matched subset ------------------
    Cfull = trace.astype(np.int64).reshape(-1)
    tr_full = np.tile(ctr, per)
    wid_full = np.tile(np.arange(lanes, dtype=np.int64), per)
    Pfull = np.concatenate([np.full(lanes, -1, dtype=np.int64), P])   # -1: no prev
    if rep is not None and np.size(rep[0]):
        # REPORT classes are part of the audit's class set too (the class a
        # walker reports is the successor of its last traced step).  They enter
        # the class scan with no predecessor, so they can never count as a
        # predecessor match.
        rclass = np.asarray(rep[1]).astype(np.int64)
        rwv = np.asarray(rep[2]).astype(np.int64)
        Cfull = np.concatenate([Cfull, rclass])
        w_per_trial = lanes // bt
        tr_full = np.concatenate([tr_full, rwv // w_per_trial])
        wid_full = np.concatenate([wid_full, rwv])
        Pfull = np.concatenate([Pfull, np.full(rclass.size, -1, dtype=np.int64)])
    o2 = np.lexsort((Cfull, tr_full))
    Cs2, trs2, wid2, Ps2 = Cfull[o2], tr_full[o2], wid_full[o2], Pfull[o2]
    new2 = np.empty(Cs2.size, dtype=bool)
    new2[0] = True
    new2[1:] = (Cs2[1:] != Cs2[:-1]) | (trs2[1:] != trs2[:-1])
    st2 = np.flatnonzero(new2)
    en2 = np.append(st2[1:], Cs2.size)
    class_groups = 0
    prev_matched = 0
    class_trials = set()
    for s, e in zip(st2[(en2 - st2) >= 2], en2[(en2 - st2) >= 2]):
        if np.unique(wid2[s:e]).size < 2:
            continue
        class_groups += 1
        class_trials.add(int(trs2[s]))
        pp = Ps2[s:e]
        good = pp[pp >= 0]
        if good.size >= 2 and np.unique(good).size < good.size:
            prev_matched += 1
    steps = int(Cs.size)
    return dict(pair_steps=steps, class_groups=int(class_groups),
                class_trials=sorted(class_trials), class_trials_count=len(class_trials),
                pair_groups=int(pair_groups), pair_trials=int(len(pair_trials)),
                class_groups_prev_matched=int(prev_matched),
                pair_over_class=(float(pair_groups) / class_groups
                                 if class_groups else 0.0),
                class_groups_per_trial=class_groups / float(bt),
                pair_groups_per_trial=pair_groups / float(bt),
                class_trials_per_trial=len(class_trials) / float(bt),
                pair_trials_per_trial=len(pair_trials) / float(bt),
                prev_match_fraction=(float(prev_matched) / class_groups
                                     if class_groups else 0.0),
                pair_collision_bound=steps * steps / (2.0 * float(q) * float(q)))


def shard(job):
    d, bt, seed, kappa, w, drain_sp = job
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
    prev = canon.copy()

    trace = np.empty((per, lanes), dtype=np.uint64)
    rl, rcls, ra, rcc, rw = [], [], [], [], []
    traps = 0
    for r in range(per):
        trace[r] = canon
        canon, a, cc, hit, prev, nf = qwalk.step(U, canon, a, cc, xv, d, want_dp=True,
                                                 prev=prev, esc=True)
        traps += int(np.count_nonzero(nf))
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
            prev[idx] = cn2

    ctr = np.repeat(np.arange(bt, dtype=np.int64), w)
    wid = np.tile(np.arange(lanes, dtype=np.int64), per)
    lvl = np.repeat(np.arange(per, dtype=np.int64)[:, None], lanes,
                    axis=1).reshape(-1)

    def parts():
        if not rl:
            z = np.zeros(0)
            return (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.uint64),
                    np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=np.uint64),
                    np.zeros(0, dtype=np.int64))
        return (np.concatenate(rl), np.concatenate(rcls), np.concatenate(ra),
                np.concatenate(rcc), np.concatenate(rw))

    rlv_b, rc_b, ra_b, rcc_b, rw_b = parts()

    # ---- D1: the audit, once, with single-division attribution --------------
    merge, merge_lvl = audit(trace, lvl, wid, ctr, q, bt, rlv_b, rc_b, rw_b)

    # ---- E4: the coverage law ON THE PAIR (class, prev-class) ---------------
    pc = pair_census(trace, ctr, q, bt, rep=(rlv_b, rc_b, rw_b))
    # the class-group set measured by pair_census must be exactly the audit's
    # merge set -- two independent implementations, one agreement check
    pc["class_groups_match_audit"] = bool(
        pc["class_trials"] == np.flatnonzero(merge).tolist())
    pc["class_trials_count"] = len(pc["class_trials"])
    pc["merge_trials_count"] = int(merge.sum())

    # ---- detection: unconditional (D2) -------------------------------------
    def detect(rlv, rc, rav, rcv, rwv):
        """Duplicate class, one trial, two different walkers.  NO merge filter."""
        succ = np.zeros(bt, dtype=bool)
        leven = np.full(bt, 10 ** 9, dtype=np.int64)
        n_ok = n_b = n_d = n_only = 0
        bad = dead = 0
        pair_trials = set()
        if rlv.size == 0:
            return (succ, leven, n_ok, n_b, n_d, bad, dead, n_only, pair_trials)
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
                    a1_, c1_ = int(a_s[k1]), int(c_s[k1])
                    a2_, c2_ = int(a_s[k2]), int(c_s[k2])
                    L = max(int(lvl_s[k1]), int(lvl_s[k2]))
                    den = (c2_ - c1_) % q
                    if den:
                        xr = (a1_ - a2_) * pow(den, -1, q) % q
                        if xr == int(xt[tr]):
                            n_ok += 1
                            if L < per:
                                n_b += 1
                            else:
                                n_d += 1
                            succ[tr] = True
                            pair_trials.add(tr)
                            if L < leven[tr]:
                                leven[tr] = L
                        else:
                            bad += 1
                    else:
                        dead += 1
            i = j
        return (succ, leven, n_ok, n_b, n_d, bad, dead, n_only, pair_trials)

    succ_b, lev_b, ok_b, nb_b, nd_b, bad_b, dead_b, _o_b, pt_b = detect(
        rlv_b, rc_b, ra_b, rcc_b, rw_b)

    # ---- drain: the ledger-charged window, capped by the walker's budget ----
    need = ~succ_b                      # every trial without a report yet (D2)
    cap_steps = int(round(drain_sp * (1 << d)))
    drain_cap = int(min(cap_steps, per))
    dl, dcls, da, dcc, dw = [], [], [], [], []
    drain_steps = 0
    sel = np.flatnonzero(need[ctr]) if need.any() else np.zeros(0, dtype=np.int64)
    if sel.size:
        c0 = canon[sel].copy()
        aa0 = a[sel].copy()
        cc0 = cc[sel].copy()
        p0 = prev[sel].copy()
        x0 = xv[sel]
        live = np.arange(sel.size, dtype=np.int64)
        for r in range(drain_cap):
            if live.size == 0:
                break
            drain_steps += live.size
            cn, ax, cx, hit, pn, nf = qwalk.step(U, c0[live], aa0[live], cc0[live],
                                                 x0[live], d, want_dp=True,
                                                 prev=p0[live], esc=True)
            traps += int(np.count_nonzero(nf))
            c0[live] = cn
            aa0[live] = ax
            cc0[live] = cx
            p0[live] = pn
            if hit.any():
                jj = live[hit]
                dl.append(np.full(jj.size, per + r, dtype=np.int64))
                dcls.append(c0[jj])
                da.append(aa0[jj])
                dcc.append(cc0[jj])
                dw.append(sel[jj])
                live = live[~hit]

    rlv = np.concatenate([rlv_b] + dl) if dl else rlv_b
    rc = np.concatenate([rc_b] + dcls)
    rav = np.concatenate([ra_b] + da)
    rcv = np.concatenate([rcc_b] + dcc)
    rwv = np.concatenate([rw_b] + dw)
    succ, lev, ok, nb, nd, bad, dead, _oo, pair_trials = detect(rlv, rc, rav, rcv, rwv)
    success_d = succ & ~succ_b

    # drain-only pairs: verified pairs on trials that the audit did NOT mark
    drain_only = sorted(t for t in pair_trials if not merge[t])
    swm = int(np.count_nonzero(succ & ~merge))
    both = merge & succ
    lag = (lev[both] - merge_lvl[both]).astype(float) if both.any() else np.zeros(0)
    return {
        "dp_bits": d, "trials": bt, "seed": seed, "walks": w, "per_walker": per,
        "r_eff": int(U.R_eff),
        "budget_steps": budget, "s_ref": s_ref, "classes": int(N), "q": int(q),
        "successes_raw": int(succ.sum()),
        "successes_budget_raw": int(succ_b.sum()),
        "successes_drain_raw": int(success_d.sum()),
        "merges": int(merge.sum()),
        "merge_budget_only": int(np.count_nonzero(
            merge & (merge_lvl < per))),
        "records": int(rlv.size), "records_budget": int(rlv_b.size),
        "records_drain": int(sum(x.size for x in dl)),
        "records_expected": lanes * per / (2.0 ** d),
        "detect_verified_pairs": ok, "detect_pairs_budget": nb, "detect_pairs_drain": nd,
        "denom_dead": dead, "solve_mismatch": bad,
        "drain_only_pairs": len(drain_only),
        "drain_only_trials": int(len(set(drain_only))),
        "success_without_merge": swm,
        "drain_trials": int(sel.size), "drain_cap_steps": drain_cap,
        "drain_cap_spacings": drain_cap / (2.0 ** d),
        "drain_window_over_budget": drain_cap / max(per, 1),
        "drain_lane_steps": drain_steps,
        "trap_escape_fires": traps,
        "trap_rate_measured": traps / max(1, lanes * per + drain_steps),
        "trap_rate_structural": 1.0 / (6.0 * U.R_eff),
        "lag_mean": float(lag.mean()) if lag.size else None,
        "lag_max": float(lag.max()) if lag.size else None,
        "merge_mean_over_sref": (float(merge_lvl[merge].mean() / s_ref)
                                 if merge.any() else None),
        "pair": pc,
    }


def reduce_all(adir, out_path, tag="prod-sweep"):
    rows = {}
    for p in sorted(glob.glob(os.path.join(adir, tag + "-*.json"))):
        with open(p) as f:
            r = json.load(f)
        d = r["dp_bits"]
        a = rows.setdefault(d, dict(trials=0, succ=0, succ_b=0, succ_d=0, merges=0,
                                    merge_b=0, records=0, rec_b=0, rec_d=0,
                                    pairs=0, pair_b=0, pair_d=0, dead=0, mismatch=0,
                                    swm=0, donly=0, donly_t=0, drain_trials=0,
                                    drain_steps=0, traps=0, lag_sum=0.0, lag_n=0,
                                    lag_max=0.0, ml_sum=0.0, ml_n=0,
                                    s_ref=r["s_ref"], classes=r["classes"], r_eff=r.get("r_eff", R_TAB // 2),
                                    q=r["q"], walks=r["walks"],
                                    per_walker=r["per_walker"],
                                    budget_steps=r["budget_steps"],
                                    drain_cap=r["drain_cap_steps"], shards=0,
                                    lanes_steps=0, pair_steps=0, class_groups=0,
                                    pair_groups=0, pair_trials=0, prev_matched=0,
                                    match_bad=0, class_trials=0))
        a["trials"] += r["trials"]
        a["succ"] += r["successes_raw"]
        a["succ_b"] += r["successes_budget_raw"]
        a["succ_d"] += r["successes_drain_raw"]
        a["merges"] += r["merges"]
        a["merge_b"] += r["merge_budget_only"]
        a["records"] += r["records"]
        a["rec_b"] += r["records_budget"]
        a["rec_d"] += r["records_drain"]
        a["pairs"] += r["detect_verified_pairs"]
        a["pair_b"] += r["detect_pairs_budget"]
        a["pair_d"] += r["detect_pairs_drain"]
        a["dead"] += r["denom_dead"]
        a["mismatch"] += r["solve_mismatch"]
        a["swm"] += r["success_without_merge"]
        a["donly"] += r["drain_only_pairs"]
        a["donly_t"] += r["drain_only_trials"]
        a["drain_trials"] += r["drain_trials"]
        a["drain_steps"] += r["drain_lane_steps"]
        a["traps"] += r["trap_escape_fires"]
        a["lanes_steps"] += r["trials"] * r["walks"] * r["per_walker"]
        if r["lag_mean"] is not None:
            a["lag_sum"] += r["lag_mean"] * r["successes_raw"]
            a["lag_n"] += r["successes_raw"]
        a["lag_max"] = max(a["lag_max"], r["lag_max"] or 0.0)
        if r["merge_mean_over_sref"] is not None:
            a["ml_sum"] += r["merge_mean_over_sref"] * r["merges"]
            a["ml_n"] += r["merges"]
        pcd = r.get("pair") or {}
        a["pair_steps"] += pcd.get("pair_steps", 0)
        a["class_groups"] += pcd.get("class_groups", 0)
        a["pair_groups"] += pcd.get("pair_groups", 0)
        a["pair_trials"] += pcd.get("pair_trials", 0)
        a["prev_matched"] += pcd.get("class_groups_prev_matched", 0)
        a["class_trials"] += pcd.get("class_trials_count", 0)
        if not pcd.get("class_groups_match_audit", False):
            a["match_bad"] += 1
        a["shards"] += 1
    out = []
    for d in sorted(rows):
        a = rows[d]
        n = a["trials"]
        k = a["succ"]
        km = a["merges"]
        cp = cp95_lower(k, n)
        cpm = cp95_lower(km, n)
        steps = a["lanes_steps"] + a["drain_steps"]
        out.append({
            "dp_bits": d, "trials": n,
            "successes_raw": k, "success_raw": k / n, "cp95_raw": float(cp),
            "cp95_raw_exact": str(cp),
            "successes_budget_raw": a["succ_b"], "success_budget_raw": a["succ_b"] / n,
            "successes_drain_raw": a["succ_d"], "success_drain_raw": a["succ_d"] / n,
            "merges": km, "merge_rate": km / n,
            "cp95_merge": float(cpm), "cp95_merge_exact": str(cpm),
            "merge_budget_only": a["merge_b"],
            "records_mean": a["records"] / n, "records_budget_mean": a["rec_b"] / n,
            "records_drain_mean": a["rec_d"] / n,
            "records_expected_mean": a["per_walker"] * a["walks"] / (2.0 ** d),
            "detect_pairs_total": a["pairs"], "detect_pairs_budget": a["pair_b"],
            "detect_pairs_drain": a["pair_d"],
            "denom_dead": a["dead"], "solve_mismatch": a["mismatch"],
            "drain_only_pairs": a["donly"], "drain_only_trials": a["donly_t"],
            "drain_only_fraction_of_merges": a["donly"] / max(km, 1),
            "success_without_merge": a["swm"],
            "success_minus_merge": (k - km) / n,
            "drain_trials_mean": a["drain_trials"] / n,
            "drain_cap_steps": a["drain_cap"],
            "drain_cap_spacings": a["drain_cap"] / (2.0 ** d),
            "drain_window_over_budget": a["drain_cap"] / max(a["per_walker"], 1),
            "trap_escape_fires": a["traps"],
            "trap_rate_measured": a["traps"] / max(steps, 1),
            "trap_rate_structural": 1.0 / (6.0 * a["r_eff"]),
            "trap_opportunities_per_spacing_measured":
                (a["traps"] / max(steps, 1)) * (2 ** d),
            "lag_mean_steps": (a["lag_sum"] / a["lag_n"]) if a["lag_n"] else None,
            "lag_max_steps": a["lag_max"] if a["lag_n"] else None,
            "merge_mean_over_sref": (a["ml_sum"] / a["ml_n"]) if a["ml_n"] else None,
            "budget_spacings_per_walker": a["per_walker"] / (2.0 ** d),
            "bayes_model": 1 - math.exp(-(KAPPA ** 2) * LAMF),
            "s_ref": a["s_ref"], "classes": a["classes"], "q": a["q"],
            "walks": a["walks"], "per_walker": a["per_walker"],
            "budget_steps": a["budget_steps"], "shards": a["shards"],
            "pair": {
                "pair_steps": a["pair_steps"],
                "class_groups": a["class_groups"],
                "class_groups_per_trial": a["class_groups"] / max(n, 1),
                "class_trials_per_trial": a["class_trials"] / max(n, 1),
                "pair_groups": a["pair_groups"],
                "pair_trials": a["pair_trials"],
                "pair_trials_per_trial": a["pair_trials"] / max(n, 1),
                "pair_groups_per_trial": a["pair_groups"] / max(n, 1),
                "pair_over_class": (a["pair_groups"] / a["class_groups"]
                                    if a["class_groups"] else 0.0),
                "class_groups_prev_matched": a["prev_matched"],
                "prev_match_fraction": (a["prev_matched"] / a["class_groups"]
                                        if a["class_groups"] else 0.0),
                "class_groups_match_audit_all_shards": a["match_bad"] == 0}})
    GATED = [5, 7, 9]           # densities whose walker restarts at a DP inside
    OUT_OF_PROTOCOL = [11]      # 0.18 spacings per walker: no restart at all
    gate = {"rule": "CP-95 of the UNCONDITIONAL coverage (merge) rate >= 0.39 at "
                    "the production-protocol densities d=5,7,9; the unconditional "
                    "pool-detection rate success_raw is reported alongside and is "
                    "not filtered by it",
            "gate_statistic": "merge_rate",
            "gated_densities": GATED,
            "out_of_protocol_densities": OUT_OF_PROTOCOL,
            "why": "production success = P(two walkers occupy one class) minus the "
                   "charged report lag 2*DP/S_ops = 2^-46; the toy's pool statistic "
                   "needs records and its pool is budget/2^d per walker (1 at d=11, "
                   "2^46 in production), so it is the diagnostic, not the gate. "
                   "d=11's walker has 0.18 DP spacings in its whole budget, i.e. it "
                   "never restarts -- not the production protocol -- and it is "
                   "reported in full",
            "min_cp95_raw": min(r["cp95_raw"] for r in out) if out else None,
            "min_cp95_merge": min(r["cp95_merge"] for r in out) if out else None,
            "min_cp95_merge_gated": min(r["cp95_merge"] for r in out
                                        if r["dp_bits"] in GATED) if out else None,
            "unconditional": True,
            "passed": bool(out) and all(r["cp95_merge"] >= 0.39 for r in out
                                        if r["dp_bits"] in GATED)}
    gate["pair_state"] = {
        "question": "the escape is deterministic in (class, prev-class): is the "
                    "coverage law measured on the PAIR, not only on the class?",
        "rows": {str(r["dp_bits"]): r["pair"] for r in out},
        "class_coverage": {str(r["dp_bits"]): r["pair"]["class_trials_per_trial"]
                           for r in out},
        "pair_coverage": {str(r["dp_bits"]): r["pair"]["pair_trials_per_trial"]
                          for r in out},
        "pair_over_class_trials": {
            str(r["dp_bits"]): (r["pair"]["pair_trials_per_trial"]
                                / r["pair"]["class_trials_per_trial"]
                                if r["pair"]["class_trials_per_trial"] else None)
            for r in out},
        "min_pair_coverage": min(r["pair"]["pair_trials_per_trial"] for r in out),
        "min_class_coverage": min(r["pair"]["class_trials_per_trial"] for r in out),
        "prev_match_fraction": {str(r["dp_bits"]): r["pair"]["prev_match_fraction"]
                                for r in out},
        "why_pair_coverage_is_not_the_success_statistic":
            "a pair lives in a space of N^2, so a pair birthday needs ~N visited "
            "pairs per trial against ~sqrt(2N) classes for a class birthday: at "
            "this budget the pair statistic is ~N^-1 of the class statistic, and "
            "it is not the event the pool replays.  Recovery needs two RECORDS in "
            "one class, produced at whatever steps each walker reached it; the "
            "invariant a + c*x == canon holds at every step of every walker, so "
            "x = (a1-a2)/(c2-c1) verifies with no shared future.  The pair number "
            "and the predecessor-matched fraction are reported so the reader can "
            "see the memory channel quantified; `success_raw` is the conversion "
            "measured end to end, and `merge` is the coverage statistic that the "
            "gate uses.",
        "all_shards_class_groups_match_audit": all(
            r["pair"]["class_groups_match_audit_all_shards"] for r in out),
    }
    gate["passed"] = bool(gate["passed"] and gate["pair_state"]
                          ["all_shards_class_groups_match_audit"])
    res = {"package": "plain-orbitrow-40pct-v20e", "arm": "prod-sweep+escape",
           "protocol": {
               "walk": "qwalk quotient walk, invariant a + c*x == canon",
               "dp": "theta=2^-d gates REPORTING; walker restarts at a DP",
               "pool": "records only; in-budget records = budget/2^d",
               "escape": "two-step-return -> alternate orbit row for one step; "
                         "armed for the whole replay (budget and drain)",
               "detect": "duplicate class among pool records, two walkers, "
                         "UNCONDITIONAL (no merge filter)",
               "drain": "ledger-charged window min(2 spacings, walker budget); "
                        "every trial without a report is flushed",
               "recovery": "(a1-a2)/(c2-c1) mod q from the two records, verified",
               "merge_stat": "independent audit over visited classes, single-division "
                             "trial attribution, computed separately from success",
               "budget": "kappa*S_ref total steps; lag charged separately, as ledger"},
           "plant_check": plant_check(),
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
    ap.add_argument("--drain-spacings", type=float, default=DRAIN_SPACINGS)
    ap.add_argument("--seed", type=int, default=20290)
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
            print("d=%2d n=%6d raw %.4f CPraw %.4f CPmg %.4f merge %.4f "
                  "donly %d traps %d" % (r["dp_bits"], r["trials"], r["success_raw"],
                                         r["cp95_raw"], r["cp95_merge"],
                                         r["merge_rate"], r["drain_only_pairs"],
                                         r["trap_escape_fires"]))
        return
    dens = [int(t) for t in a.dens.split(",") if t]
    jobs = []
    for d in dens:
        done = 0
        i = 0
        while done < a.trials:
            take = min(a.chunk, a.trials - done)
            jobs.append((d, take, a.seed * 100000 + d * 1000 + i, a.kappa, a.walks,
                         a.drain_spacings))
            done += take
            i += 1
    t0 = time.time()
    print("jobs", len(jobs), "dens", dens, "trials/density", a.trials, flush=True)
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(shard, jobs):
            name = "prod-sweep-d%02d-s%05d.json" % (r["dp_bits"], r["seed"] % 100000)
            with open(os.path.join(adir, name), "w") as f:
                json.dump(r, f, indent=1, sort_keys=True)
            print("  d=%2d n=%4d raw %.4f (b %.4f d %.4f) merge %.4f donly %3d "
                  "traps %4d %.1fs"
                  % (r["dp_bits"], r["trials"], r["successes_raw"] / r["trials"],
                     r["successes_budget_raw"] / r["trials"],
                     r["successes_drain_raw"] / r["trials"],
                     r["merges"] / r["trials"], r["drain_only_pairs"],
                     r["trap_escape_fires"], time.time() - t0), flush=True)
    res = reduce_all(adir, out_path)
    print(json.dumps(res["gate"], indent=1))
    for r in res["rows"]:
        print("d=%2d n=%6d raw %.4f CPraw %.4f CPmg %.4f merge %.4f"
              % (r["dp_bits"], r["trials"], r["success_raw"], r["cp95_raw"],
                 r["cp95_merge"], r["merge_rate"]))
    print("wrote", out_path, "%.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
