#!/usr/bin/env python3
"""trails_v25.py -- RAW TRAIL-LENGTH HISTOGRAMS and the EXACT-SUM success.

Verdict da40eb9c point 2 killed the exp-of-mean substitution:
    miss = E[(1-2^-60)^T]  is NOT  exp(-2^-60 * E[T]).
f(T) = (1-p)^T is CONVEX, so Jensen gives
    E[(1-p)^T] >= (1-p)^{E[T]} = exp(E[T]*log(1-p)),
i.e. the exp-of-mean number is a LOWER bound on the miss, hence an UPPER bound
on the success.  v24 shipped the exp-of-mean value as the success.  This module
records T for every trial (visits actually walked before the verified solve) and
sums the exact finite sum

    success_exact(p) = 1 - (1/n) * sum_j (1-p)^{T_j}

over the EMPIRICAL distribution -- no moment shortcut, no fitted tail.

SCALE / p.  p is fixed by the ledger's own design point, not fitted here: the
design row is "expected DP opportunities in a trail is 1", i.e. p = 1/E[T] with
E[T] the measured mean trail length in the SAME cell.  With that p the
exp-of-mean value is the constant 1-e^-1 = 0.632121, so the two routes are
compared at identical mean by construction and the only thing that can move the
number is the SHAPE of the measured T distribution.  The absolute p of a given
cell is reported alongside; the correction factor E[(1-p)^T]/(1-p)^{E[T]} is
scale-free and is what transfers.

WITHIN-HISTORY CELLS (verdict point 1).  Every cell is estimated on its own
data; nothing is pooled across cells for the bound.  Cells vary (L, R, d) -- the
three coordinates the verdict said the transfer held fixed -- and each cell also
splits its own trials by the escape history (esc_run fired / clean), so the
history-dependent part of the walk is a separate stratum with its own count,
its own CP lower bound, and its own exact sum.

Output: evidence/trail-hist-v25.json (histogram buckets + exact sums + CP bounds).
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dladder as D  # noqa: E402


def cp_lower(k, n, alpha=0.05):
    return D.cp_lower(k, n, alpha)


def run_cell(tc, trials, seed, L=16, R=32, dgate=0, pairs=False, esc_depth=3,
             V_mult=1.0):
    """One cell, recording the RAW trail length T of every run (visits walked)."""
    import numpy as np
    C = D.find_fixture(tc)
    G = D.make_G(C)
    D.FIXED_G[C.p] = G
    C.resolve_omega(G)
    S_OPS = 1.009696 * math.sqrt(C.N_CLASSES)
    V = max(L, int(round(V_mult * S_OPS)))
    rng = np.random.default_rng(seed)
    rows = []          # per challenge: (T1, ok1, fires1, T2, ok2, fires2)
    t0 = time.time()
    for _tr in range(trials):
        x = int(rng.integers(1, C.N))
        Q = C.mul(G, x)
        tab = D.Table(C, R, Q, seed * 100003 + _tr)
        r1 = D.run_trial(C, tab, G, Q, V, L, dgate=dgate, esc=True,
                         esc_depth=esc_depth, rng=rng, R=R)
        if pairs:
            tab = D.Table(C, R, Q, seed * 100003 + _tr)
            r2 = D.run_trial(C, tab, G, Q, V, L, dgate=dgate, esc=True,
                             esc_depth=esc_depth, rng=rng, R=R)
            f1 = sum(r1["fires"].values()) if isinstance(r1["fires"], dict) else r1["fires"]
            f2 = sum(r2["fires"].values()) if isinstance(r2["fires"], dict) else r2["fires"]
            rows.append((r1["visits"], int(r1["ok"]), f1, r1["restarts"],
                         r2["visits"], int(r2["ok"]), f2, r2["restarts"]))
        else:
            f1 = sum(r1["fires"].values()) if isinstance(r1["fires"], dict) else r1["fires"]
            rows.append((r1["visits"], int(r1["ok"]), f1, r1["restarts"]))
    el = time.time() - t0
    return dict(cell=dict(target_classes=tc, trials=trials, seed=seed, L=L, R=R,
                          dgate=dgate, esc_depth=esc_depth, V=V, N=C.N,
                          N_CLASSES=C.N_CLASSES, pairs=pairs, elapsed_s=round(el, 3)),
                rows=rows)


def hist_of(vals, nb=48):
    """Bucket histogram, ASCII/JSON-safe, counts only."""
    if not vals:
        return dict(edges=[], counts=[], n=0)
    lo, hi = min(vals), max(vals)
    if hi == lo:
        return dict(edges=[lo, hi + 1], counts=[len(vals)], n=len(vals))
    w = (hi - lo) / nb
    counts = [0] * nb
    for v in vals:
        i = min(nb - 1, int((v - lo) / w))
        counts[i] += 1
    return dict(edges=[round(lo + i * w, 6) for i in range(nb)] + [round(hi, 6)],
                counts=counts, n=len(vals), lo=lo, hi=hi, width=round(w, 6))


def moments(vals):
    n = len(vals)
    if not n:
        return dict(n=0)
    m = sum(vals) / n
    var = sum((v - m) ** 2 for v in vals) / n
    sd = math.sqrt(var)
    return dict(n=n, mean=m, var=var, sd=sd, cv=(sd / m if m else None),
                min=min(vals), max=max(vals))


def exact_success(vals, p):
    """The exact finite sum 1 - (1/n) sum (1-p)^{T_j}.  No moment shortcut."""
    n = len(vals)
    miss = sum((1.0 - p) ** t for t in vals) / n
    return 1.0 - miss, miss


def exact_sum_block(vals, censored=None, tag=""):
    """Exact sums at the cell's own design point plus the exp-of-mean twin."""
    mo = moments(vals)
    if not mo["n"]:
        return dict(tag=tag, n=0)
    E = mo["mean"]
    p = 1.0 / E                      # design point: E[pT] = 1
    s_exact, miss_exact = exact_success(vals, p)
    miss_eom = (1.0 - p) ** E        # the v24 substitution
    return dict(tag=tag, n=mo["n"], mean_T=E, sd_T=mo["sd"], cv_T=mo["cv"],
                p_design=p, log10_p_design=math.log10(p),
                success_exact=s_exact, miss_exact=miss_exact,
                success_exp_of_mean=1.0 - miss_eom,
                miss_exp_of_mean=miss_eom,
                miss_ratio_exact_over_eom=(miss_exact / miss_eom),
                jensen_gap_success=(1.0 - miss_eom) - s_exact,
                censored=censored)


def pair_block(rows):
    """Per-challenge two-run composition, exact vs exp-of-mean, + CP lower bounds."""
    T1 = [r[0] for r in rows]
    T2 = [r[4] for r in rows]
    ok = [1 if (r[1] or r[5]) else 0 for r in rows]
    n = len(rows)
    k = sum(ok)
    both = [t for r in rows for t in (r[0], r[4])]
    E = sum(both) / len(both)
    p = 1.0 / E
    s_j = [1.0 - (1.0 - p) ** r[0] * (1.0 - p) ** r[4] for r in rows]
    a_j = [1.0 - (1.0 - p) ** (r[0] + r[4]) for r in rows]   # exp-of-mean twin
    frac_k = sum(s_j)
    # fractional-count Clopper-Pearson: k_eff = sum s_j, n = challenges.  This is
    # the conservative reading of a mean of non-0/1 per-challenge successes.
    cp_meas = cp_lower(k, n)
    cp_frac_exact = cp_lower(max(0.0, math.floor(frac_k)), n) if frac_k >= 1 else 0.0
    cp_frac_exact_hi = cp_lower(min(n, math.ceil(frac_k)), n) if frac_k >= 1 else 0.0
    return dict(n=n, measured_set_success=k / n, cp95_measured=cp_meas,
                mean_T_run=E, p_design=p,
                mean_challenge_success_exact=sum(s_j) / n,
                mean_challenge_success_eom=sum(a_j) / n,
                frac_k_exact=frac_k,
                cp95_frac_floor=cp_frac_exact, cp95_frac_ceil=cp_frac_exact_hi,
                censored_run1=sum(1 for r in rows if not r[1]),
                censored_run2=sum(1 for r in rows if not r[5]),
                hist_T_run1=hist_of(T1), hist_T_run2=hist_of(T2))


def cell_block(res, do_pairs=False):
    rows = res["rows"]
    c = res["cell"]
    if do_pairs:
        Tsucc = [tt for r in rows for tt, okk in ((r[0], r[1]), (r[4], r[5])) if okk]
        Tall = [tt for r in rows for tt in (r[0], r[4])]
        esc = [tt for r in rows for tt, ff in ((r[0], r[2]), (r[4], r[6])) if ff > 0]
        clean = [tt for r in rows for tt, ff in ((r[0], r[2]), (r[4], r[6])) if ff == 0]
    else:
        Tsucc = [r[0] for r in rows if r[1]]
        Tall = [r[0] for r in rows]
        esc = [r[0] for r in rows if r[2] > 0]
        clean = [r[0] for r in rows if r[2] == 0]
    out = dict(cell=c, hist_T_success=hist_of(Tsucc), hist_T_all=hist_of(Tall),
               all=(exact_sum_block(Tall, censored=len(Tall) - len(Tsucc), tag="all_runs")),
               success_only=exact_sum_block(Tsucc, tag="successful_runs_only"),
               history_escape_fired=exact_sum_block(esc, tag="esc_run_fired"),
               history_clean=exact_sum_block(clean, tag="esc_run_clean"))
    if do_pairs:
        out["pairs"] = pair_block(rows)
    return out


def main():
    tc17, tc20 = 2 ** 17, 2 ** 20
    t = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
    out = dict(test="trail_hist_exact_sum_v25",
               law="success_exact = 1 - (1/n) sum_j (1-p)^{T_j}, p = 1/mean(T)",
               cells=[])
    # (A) the shipped w2 design cell, per-challenge trails
    out["cells"].append(cell_block(run_cell(tc20, t, 13, L=16, R=32, dgate=0,
                                            pairs=True), do_pairs=True))
    # (B) within-history re-ladder: same class count, L/R/d moved one at a time
    for (L, R, d) in ((16, 32, 4), (32, 64, 0), (64, 32, 0)):
        out["cells"].append(cell_block(run_cell(tc17, t, 12, L=L, R=R, dgate=d)))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
