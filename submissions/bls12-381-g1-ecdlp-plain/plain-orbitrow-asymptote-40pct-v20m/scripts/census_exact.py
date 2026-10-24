#!/usr/bin/env python3
"""census_exact.py -- v20m: escaped-walk census at the EXACT ledger-normalized budget.

Reviewer 2d98d2a8 point 2: the v20l census ran 1.1342 visits/sqrt(N) while the
ledger budget is S_OPS/sqrt(N_CLASSES) = 1.009696311282301.  This script solves
the fixture shape FROM the ledger's own S_OPS (same Decimal enclosure +
Fractions recipe as ledger.py) and asserts, per arm, that the realized
visits/sqrt(N) equals it to 1e-3 (hard assert, no arm ships otherwise).

Mechanism decomposition: success(d) = P(two walkers' class sets intersect) x
P(a DP class lies in the shared tail | intersection).  The second factor is
EXACTLY computable from the same trajectories at d=0 (every class reports),
where success = intersection itself.  Per size we run (a) the production arm at
d=5 and (b) a d=0 arm on the SAME seed and per-trial coins, store its class
trace (uint32: q < 2^32), and compute the exact conditional detection factor.

Restart counter is accumulated from the ACTUAL restart_out mask (the v20l
`restarts = 0` never-updated bug); per-lane peak escape run is histogrammed.
"""
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = HERE
sys.path.insert(0, HERE)
import qwalk  # noqa: E402
import selector  # noqa: E402

# ---- EXACT ledger budget (same recipe as ledger.py) ------------------------
p = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624'
        '1eabfffeb153ffffb9feffffffffaaab', 16)
r = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
R_TAB = 2 ** 20
R_EFF = selector.reachable_columns(R_TAB, p)
N_CELLS = 6 * R_EFF
N_CLASSES = (r - 1) // 6
A_H1 = F(1015, 1000)
MODEL_SUCCESS = F(3903, 10000)
getcontext().prec = 80
_z0 = Decimal(1) - Decimal(MODEL_SUCCESS.numerator) / Decimal(MODEL_SUCCESS.denominator)
_ln0 = -_z0.ln()


def frac_of_dec(dd, digits=70):
    return F(int(dd * Decimal(10) ** digits) + (1 if dd >= 0 else -1), 10 ** digits)


Nprime = F(N_CLASSES * N_CELLS, N_CELLS - 1)
_x2 = 2 * Nprime * frac_of_dec(_ln0)
S0 = math.isqrt(_x2.numerator // _x2.denominator)
if F(S0 * S0) < _x2:
    S0 += 1
_S = A_H1 * A_H1 * F(S0) * F(S0)
S_OPS = math.isqrt(_S.numerator // _S.denominator)
if F(S_OPS * S_OPS) < _S:
    S_OPS += 1
SHAPE = S_OPS / math.sqrt(N_CLASSES)          # = 1.009696311282301
assert abs(SHAPE - 1.009696311282301) < 1e-12, SHAPE

ALPHA = Decimal("0.05")


def cp95_lower(k, n):
    """One-sided Clopper-Pearson lower bound (bisection; v15k helper)."""
    k, n = int(k), int(n)
    if k == 0:
        return Decimal(0)
    if k >= n:
        return ALPHA ** (Decimal(1) / Decimal(n))
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
            if t < Decimal("1e-40") * s:
                break
        return s

    lo, hi = Decimal("0.001"), Decimal("0.999")
    for _ in range(300):
        mid = (lo + hi) / 2
        if tail(mid) > ALPHA:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-14"):
            break
    assert Decimal("0.001") < lo < Decimal("0.999"), ("cp95 clipped", lo, k, n)
    return lo


def make_start(U, bt, w, seed):
    q = U.q
    rng = np.random.default_rng(seed & 0x7FFFFFFF)
    xt = rng.integers(1, q - 1, size=bt).astype(np.uint64)
    xv = np.repeat(xt, w).astype(np.uint64)
    lanes = bt * w
    a = rng.integers(0, q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, q, size=lanes).astype(np.uint64)
    t0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                   for ai, ci, xi in zip(a.tolist(), cc.tolist(), xv.tolist())],
                  dtype=np.uint64)
    canon, e, _h = U.canon(t0)
    a = (e * a) % q
    cc = (e * cc) % q
    return xt, xv, a, cc, canon


def run_prod(U, bt, d, w, seed, per, want_trace=False, label=""):
    """The shipped qwalk.step, escape armed, kcap live, restart counter real."""
    q = U.q
    xt, xv, a, cc, canon = make_start(U, bt, w, seed)
    lanes = bt * w
    prev = canon.copy()
    esc_run = np.zeros(lanes, dtype=np.int64)
    max_run = np.zeros(lanes, dtype=np.int64)
    coins = np.random.default_rng((seed ^ 0x5EED5EED) & 0x7FFFFFFF)
    restarts = 0
    fires = 0
    first = [dict() for _ in range(bt)]
    succ_step = np.full(bt, -1, dtype=np.int64)
    trace = np.empty((per, lanes), dtype=np.uint32) if want_trace else None
    t_start = time.time()
    for t in range(per):
        if want_trace:
            trace[t] = canon.astype(np.uint32)
        rst = np.zeros(lanes, dtype=bool)
        canon, a, cc, hit, prev, fired = qwalk.step(
            U, canon, a, cc, xv, d, want_dp=True, prev=prev, esc=True, tick=t,
            esc_run=esc_run, kcap=qwalk.K_ESC, restart_out=rst, coins=coins)
        restarts += int(rst.sum())
        fires += int(fired.sum())
        np.maximum(max_run, esc_run, out=max_run)
        if d > 0 and hit is not None and hit.any():
            idx = np.flatnonzero(hit)
            for i in idx.tolist():
                b = i // w
                if succ_step[b] >= 0:
                    continue
                cl = int(canon[i])
                fd = first[b]
                rec = fd.get(cl)
                if rec is None:
                    fd[cl] = (int(a[i]), int(cc[i]), i)
                elif rec[2] != i:
                    a1, c1, _l1 = rec
                    den = (int(cc[i]) - c1) % q
                    if den:
                        cand = ((a1 - int(a[i])) * pow(den, -1, q)) % q
                        if cand == int(xt[b]):
                            succ_step[b] = t
    elapsed = time.time() - t_start
    succ = int(np.count_nonzero(succ_step >= 0))
    cp = cp95_lower(succ, bt) if d > 0 else None
    visit_shape = per * w / math.sqrt(U.classes)
    dev = abs(visit_shape - float(SHAPE))
    assert dev <= 1e-3, ("BUDGET MISMATCH", label, visit_shape, float(SHAPE), dev)
    step_hist = {}
    for t in succ_step[succ_step >= 0].tolist():
        bucket = int(10 * t / per)
        step_hist[bucket] = step_hist.get(bucket, 0) + 1
    out = {
        "label": label, "q": int(q), "classes": int(U.classes),
        "walks": int(w), "per_walker_steps": int(per), "trials": int(bt),
        "dp_bits": int(d), "visits_per_trial": int(per * w),
        "visits_over_sqrtN": visit_shape,
        "shape_assert_dev": dev, "ledger_shape": float(SHAPE),
        "successes": succ, "success_rate": succ / float(bt),
        "cp95": str(cp) if cp is not None else None,
        "cp95_f": float(cp) if cp is not None else None,
        "restarts": restarts, "escape_fires": fires,
        "max_escape_run": int(max_run.max(initial=0)),
        "peak_run_hist": {int(k): int(v) for k, v in
                          enumerate(np.bincount(max_run, minlength=10).tolist()) if v},
        "success_step_decile_hist": step_hist,
        "elapsed_s": round(elapsed, 2),
    }
    return out, trace, succ_step


def coalescence_from_trace(U, trace, w, dp_bits):
    """Exact offline decomposition on the d=0 trajectories."""
    steps, lanes = trace.shape
    bt = lanes // w
    coalesced = np.zeros(bt, dtype=bool)
    detect = np.zeros(bt, dtype=bool)
    meet = np.full(bt, -1, dtype=np.int64)
    tail_len = np.zeros(bt, dtype=np.int64)
    for b in range(bt):
        cols = trace[:, b * w:(b + 1) * w]
        flat = cols.ravel()
        order = np.argsort(flat, kind="stable")
        sv = flat[order]
        brk = np.flatnonzero(np.r_[True, sv[1:] != sv[:-1]])
        ends = np.r_[brk[1:], sv.size]
        g = np.flatnonzero(ends - brk >= 2)
        if g.size == 0:
            continue
        best = None
        for gi in g.tolist():
            idxs = order[brk[gi]:ends[gi]]
            lanes_g = idxs % w
            if np.unique(lanes_g).size < 2:
                continue          # a walker's own revisit: not a meeting
            l0 = int(lanes_g[0])
            i2 = int(idxs[lanes_g != l0].min())   # earliest second-walker arrival
            if best is None or i2 < best:
                best = i2
        if best is None:
            continue
        coalesced[b] = True
        t2, l2 = divmod(best, w)
        meet[b] = t2
        tail_len[b] = steps - t2
        detect[b] = bool(U.dp(cols[t2:, l2].astype(np.uint64), dp_bits).any())
    return coalesced, detect, meet, tail_len


def main():
    t_all = time.time()
    # class-space sizes: 6.4e7 (2^25.9), 3.24e8 (2^28.3), 5.62e8 (2^29.1).
    # q < 2^32 is required by the uint64 canon (lam*c must not wrap).
    arms = [
        (int(round(6 * 8000.0 ** 2)), 8192, 2048),
        (int(round(6 * 18000.0 ** 2)), 5120, 2048),
        (int(round(6 * 23700.0 ** 2)), 4096, 2048),
    ]
    W = 16
    d = 5
    seed = 20250
    out = {"test": "walk_census_exact_budget", "ledger_shape": float(SHAPE),
           "protocol": ("escape armed (kcap=8, shipped step, restart counter "
                        "accumulated from restart_out), DP density 2^-5, success "
                        "= two DP records of one class from two walkers of one "
                        "challenge solving x; d=0 arm on the same seed measures "
                        "coalescence and the exact detection factor"),
           "walks": W, "dp_bits": d, "arms": [], "controls": []}
    for q_min, bt, bt0 in arms:
        U = qwalk.Universe(q_min, 1 << 16, seed=seed)
        per = int(round(float(SHAPE) * math.sqrt(U.classes) / W))
        print(f"[arm] q={U.q} classes={U.classes} per={per} bt={bt} "
              f"shape={per*W/math.sqrt(U.classes):.9f}", flush=True)
        r, _tr, ss = run_prod(U, bt, d, W, seed, per, label=f"d5_q{U.q}")
        r["successes_first_bt0"] = int(np.count_nonzero(ss[:bt0] >= 0))
        out["arms"].append(r)
        print(json.dumps(r), flush=True)
        # d=0 control on the same seed: same trajectories, every class reports
        r0, trace, _ss0 = run_prod(U, bt0, 0, W, seed, per, want_trace=True,
                                   label=f"d0_q{U.q}")
        coal, det, meet, tail = coalescence_from_trace(U, trace, W, d)
        del trace
        n = bt0
        coal_n = int(np.count_nonzero(coal))
        det_n = int(np.count_nonzero(coal & det))
        r0.update({
            "coalesced": coal_n, "coalescence_rate": coal_n / n,
            "cp95_coalescence": str(cp95_lower(coal_n, n)),
            "detected": det_n,
            "detection_given_coalescence": det_n / coal_n if coal_n else None,
            "predicted_success": det_n / n,
            "predicted_cp95": str(cp95_lower(det_n, n)),
            "measured_success_first_bt0": r["successes_first_bt0"],
            "measured_rate_first_bt0": r["successes_first_bt0"] / n,
            "mean_meet_step": float(meet[coal].mean()) if coal_n else None,
            "mean_shared_tail": float(tail[coal].mean()) if coal_n else None,
            "mean_tail_over_per": float(tail[coal].mean() / per) if coal_n else None,
            "birthday_pred": 1 - math.exp(-(per * W) * (per * W - 1) / (2.0 * U.classes)),
        })
        out["controls"].append(r0)
        print(json.dumps(r0), flush=True)
    out["total_elapsed_s"] = round(time.time() - t_all, 1)
    with open(os.path.join(D, "evidence", "walk-census-exact.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("[done]", out["total_elapsed_s"], "s", flush=True)


if __name__ == "__main__":
    main()
