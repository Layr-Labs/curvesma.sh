#!/usr/bin/env python3
"""census_v20n.py -- v20n matched-pair census at the EXACT ledger budget.

Fixes verdict e9ee1dbc point 2 exactly:
  * make_start draws a FIXED BT_MAX challenge/lane block, so every arm of trial
    b starts from the SAME (x, a, c) and walks the SAME trajectory (v20m drew
    `bt` challenges first: a 5120-trial arm never shared a start with its
    2048-trial control).
  * d=0 is not a second walk: it is the SAME trace read with the report gate
    off, so coalescence and the detection factor are on identical trajectories.
  * escape vs no-escape is a MATCHED PAIR: same start, same coins, identical
    until the first escape fire; the treatment substitutes the alternate row at
    the trigger, the control does not.  Prefix equality is ASSERTED and the
    per-pair trajectory hashes are shipped.

Also emitted: a lag-2/lag-3 escape variant (3-cycle channel), the cycle-length
census {2,3,>=4} of residual absorption, and the cap-loss row
P(would-solve step > cap | would-solve) read straight off the trail histogram.
"""
import hashlib
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction as F

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402
import selector  # noqa: E402

p = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624'
        '1eabfffeb153ffffb9feffffffffaaab', 16)
r = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
N_CLASSES = (r - 1) // 6
R_TAB = 2 ** 20
R_EFF = selector.reachable_columns(R_TAB, p)
N_CELLS = 6 * R_EFF
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
SHAPE = S_OPS / math.sqrt(N_CLASSES)
assert abs(SHAPE - 1.009696311282301) < 1e-12, SHAPE
ALPHA = Decimal("0.05")


def cp95_lower(k, n):
    """One-sided Clopper-Pearson lower bound (bisection; v15k helper)."""
    k, n = int(k), int(n)
    if n == 0 or k == 0:
        return Decimal(0)
    if k >= n:
        return ALPHA ** (Decimal(1) / Decimal(n))
    lnC = (sum(Decimal(i).ln() for i in range(n - k + 1, n + 1))
           - sum(Decimal(i).ln() for i in range(1, k + 1)))

    def tail(p_):
        p_ = Decimal(p_)
        q_ = Decimal(1) - p_
        t = (lnC + k * p_.ln() + (n - k) * q_.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p_ / ((j + 1) * q_)
            s += t
            if t < Decimal("1e-40") * s:
                break
        return s

    lo, hi = Decimal("0.0005"), Decimal("0.999")
    for _ in range(300):
        mid = (lo + hi) / 2
        if tail(mid) > ALPHA:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal("1e-14"):
            break
    return lo


def cp95_upper(k, n):
    """One-sided 95% upper bound on a rate, from the same bisection."""
    if n == 0:
        return Decimal(1)
    if k == 0:
        return Decimal(1) - ALPHA ** (Decimal(1) / Decimal(n))
    return Decimal(1) - cp95_lower(n - k, n)


def make_start(U, bt_max, w, seed):
    """FIXED-size draw: every arm takes a prefix of these lanes (verdict pt 2)."""
    q = U.q
    rng = np.random.default_rng(seed & 0x7FFFFFFF)
    xt = rng.integers(1, q - 1, size=bt_max).astype(np.uint64)
    xv = np.repeat(xt, w).astype(np.uint64)
    lanes = bt_max * w
    a = rng.integers(0, q, size=lanes).astype(np.uint64)
    cc = rng.integers(0, q, size=lanes).astype(np.uint64)
    t0 = np.array([(int(ai) + int(ci) * int(xi)) % q
                   for ai, ci, xi in zip(a.tolist(), cc.tolist(), xv.tolist())],
                  dtype=np.uint64)
    canon, e, _h = U.canon(t0)
    return xt, xv, (e * a) % q, (e * cc) % q, canon


def walk_arm(U, canon0, a0, cc0, xv, d, per, esc, lag, w, bt):
    """One arm; returns trace (per, lanes) uint32, solve steps, first-fire steps."""
    lanes = canon0.size
    canon = canon0.copy(); a = a0.copy(); cc = cc0.copy()
    prev = canon.copy(); prev2 = canon.copy(); prev3 = canon.copy()
    esc_run = np.zeros(lanes, dtype=np.int64)
    coins = np.random.default_rng(0xC0FFEE)
    trace = np.empty((per, lanes), dtype=np.uint32)
    first_fire = np.full(bt, -1, dtype=np.int64)
    first = [dict() for _ in range(bt)]
    succ_step = np.full(bt, -1, dtype=np.int64)
    fires = restarts = 0
    for t in range(per):
        trace[t] = canon.astype(np.uint32)
        rst = np.zeros(lanes, dtype=bool)
        cn, a, cc, hit, pv, fired = qwalk.step(
            U, canon, a, cc, xv, d, True, prev, esc, 0, t, esc_run,
            (qwalk.K_ESC if esc else None), rst, coins)
        if lag >= 2:
            for other in (prev2, prev3) if lag >= 3 else (prev2,):
                extra = (cn == other) & ~fired
                if extra.any():
                    idx = np.flatnonzero(extra)
                    m = (U.col(canon[idx]) + 1) % U.R
                    cn2, a3, cc3, _e = qwalk._apply(U, canon[idx], a[idx], cc[idx], xv[idx], m)
                    cn = cn.copy(); a = a.copy(); cc = cc.copy()
                    cn[idx] = cn2; a[idx] = a3; cc[idx] = cc3
                    fired = fired | extra
        canon = cn                      # ADVANCE the class (v20n shipped without this)
        fires += int(fired.sum()); restarts += int(rst.sum())
        if fired.any():
            fb = np.flatnonzero(np.bincount(np.flatnonzero(fired) // w, minlength=bt) > 0)
            fresh = fb[first_fire[fb] < 0]
            first_fire[fresh] = t
        prev3 = prev2; prev2 = prev; prev = pv
        if d > 0 and hit is not None and hit.any():
            for i in np.flatnonzero(hit).tolist():
                b = i // w
                if succ_step[b] >= 0:
                    continue
                cl = int(cn[i])
                rec = first[b].get(cl)
                if rec is None:
                    first[b][cl] = (int(a[i]), int(cc[i]), i)
                elif rec[2] != i:
                    a1, c1, _l1 = rec
                    den = (int(cc[i]) - c1) % U.q
                    if den:
                        cand = ((a1 - int(a[i])) * pow(den, -1, U.q)) % U.q
                        if cand == int(xv[i]):
                            succ_step[b] = t
    return {"trace": trace, "succ_step": succ_step, "fires": fires,
            "restarts": restarts, "first_fire": first_fire}


def decomposition(U, trace, w, bt, d):
    """d=0 coalescence, detection factor, meeting step, terminal cycle length."""
    steps = trace.shape[0]
    coal = np.zeros(bt, dtype=bool)
    det = np.zeros(bt, dtype=bool)
    cyc = np.zeros(bt, dtype=np.int64)
    meet = np.full(bt, -1, dtype=np.int64)
    for b in range(bt):
        cols = trace[:, b * w:(b + 1) * w].astype(np.uint64)
        flat = cols.ravel()
        order = np.argsort(flat, kind="stable")
        sv = flat[order]
        brk = np.flatnonzero(np.r_[True, sv[1:] != sv[:-1]])
        ends = np.r_[brk[1:], sv.size]
        g = np.flatnonzero(ends - brk >= 2)
        if g.size == 0:
            continue
        best_meet = None
        term_len = 0
        for gi in g.tolist():
            idxs = order[brk[gi]:ends[gi]]
            lg = idxs % w
            two = np.sort(idxs)
            # cycle period = STEPS between two revisits by the SAME lane.
            # Flattened index gaps are w x the step gap (a same-lane two-step
            # return has index gap 2w), and the trial's basin is the SHORTEST
            # lane period (the first lane to be absorbed), not the longest.
            for lane in np.unique(two % w).tolist():
                sel = two[two % w == lane]
                if sel.size >= 2:
                    gperiod = int((sel[-1] - sel[-2]) // w)
                    term_len = gperiod if term_len == 0 else min(term_len, gperiod)
            if np.unique(lg).size < 2:
                continue
            l0 = int(lg[0])
            i2 = int(idxs[lg != l0].min())
            if best_meet is None or i2 < best_meet:
                best_meet = i2
        cyc[b] = term_len
        if best_meet is not None:
            coal[b] = True
            t2, l2 = divmod(best_meet, w)
            meet[b] = t2
            det[b] = bool(U.dp(cols[t2:, l2], d).any())
    return coal, det, meet, cyc


def main():
    t_all = time.time()
    W, d, BT_MAX, CHUNK, seed = 16, 5, 3200, 400, 20250
    U = qwalk.Universe(4_290_000_000, 64, seed=seed)   # R_eff=32: trap channel exercised at 1/(6*32)/step, a DISCLOSED stress rate (production 1/(6*852105))
    assert U.q < (1 << 32), ("uint64 canon wrap", U.q)
    per = int(round(float(SHAPE) * math.sqrt(U.classes) / W))
    shape_real = per * W / math.sqrt(U.classes)
    assert abs(shape_real - float(SHAPE)) <= 1e-3, (shape_real, float(SHAPE))
    xt, xv, a_all, cc_all, canon_all = make_start(U, BT_MAX, W, seed)
    arms = {"escaped_lag1": (True, 1), "escaped_lag3": (True, 3),
            "control_noescape": (False, 0)}
    acc = {k: {"succ": 0, "coal": 0, "det": 0, "len2": 0, "len3": 0, "ge4": 0,
               "norep": 0, "fires": 0, "restarts": 0, "solved": [],
               "meets": [], "fire_trials": 0} for k in arms}
    pairs = {"pairs_with_trigger": 0, "pairs_prefix_equal": 0, "rows": []}
    for c0 in range(0, BT_MAX, CHUNK):
        bt = min(CHUNK, BT_MAX - c0)
        sl = slice(c0 * W, (c0 + bt) * W)
        res = {}
        for k, (esc, lag) in arms.items():
            res[k] = walk_arm(U, canon_all[sl], a_all[sl], cc_all[sl], xv[sl],
                              d, per, esc, lag, W, bt)
            coal, det, meet, cyc = decomposition(U, res[k]["trace"], W, bt, d)
            A = acc[k]
            A["succ"] += int(np.count_nonzero(res[k]["succ_step"] >= 0))
            A["coal"] += int(np.count_nonzero(coal))
            A["det"] += int(np.count_nonzero(coal & det))
            A["len2"] += int(np.count_nonzero(cyc == 2))
            A["len3"] += int(np.count_nonzero(cyc == 3))
            A["ge4"] += int(np.count_nonzero(cyc >= 4))
            A["norep"] += int(np.count_nonzero(cyc == 0))
            A["fires"] += res[k]["fires"]
            A["restarts"] += res[k]["restarts"]
            A["solved"] += res[k]["succ_step"][res[k]["succ_step"] >= 0].tolist()
            A["meets"] += meet[meet >= 0].tolist()
            A["fire_trials"] += int(np.count_nonzero(res[k]["first_fire"] >= 0))
        # matched-pair prefix equality (treatment vs control, same start)
        tt, tc = res["escaped_lag1"]["trace"], res["control_noescape"]["trace"]
        ff = res["escaped_lag1"]["first_fire"]
        for b in np.flatnonzero(ff >= 0).tolist():
            t = int(ff[b])
            pt = tt[:t, b * W:(b + 1) * W]
            pc = tc[:t, b * W:(b + 1) * W]
            eq = bool(pt.shape == pc.shape and np.array_equal(pt, pc))
            pairs["pairs_with_trigger"] += 1
            pairs["pairs_prefix_equal"] += int(eq)
            if len(pairs["rows"]) < 6:
                pairs["rows"].append({
                    "trial": c0 + b, "branch_step": t,
                    "prefix_sha256_treat": hashlib.sha256(pt.tobytes()).hexdigest()[:16],
                    "prefix_sha256_control": hashlib.sha256(pc.tobytes()).hexdigest()[:16],
                    "equal": eq, "prefix_lanes": int(pt.size)})
    assert pairs["pairs_with_trigger"] > 0 and \
        pairs["pairs_prefix_equal"] == pairs["pairs_with_trigger"], pairs
    out = {"test": "census_v20n_matched_pairs",
           "protocol": ("shipped qwalk.step; escape armed kcap=8; DP 2^-5; one "
                        "class universe; FIXED seed block so every arm shares "
                        "starts; treatment substitutes the alt row at the "
                        "trigger, control does not; d=0 is the same trace with "
                        "the gate off"),
           "ledger_shape": float(SHAPE), "q": int(U.q), "classes": int(U.classes),
           "walks_per_trial": W, "trials": BT_MAX, "dp_bits": d,
           "per_walker_steps": per, "visits_per_trial": per * W,
           "visits_over_sqrtN": shape_real,
           "shape_assert_dev": abs(shape_real - float(SHAPE)),
           "arms": {}, "cap_loss": {}, "cycle_length": {}, "matched_pairs": {}}
    n = BT_MAX
    for k in arms:
        A = acc[k]
        cap = {}
        meets = sorted(A["meets"])
        for frac in (0.25, 0.5, 1.0):
            cstep = int(per * frac)
            over = sum(1 for s_ in meets if s_ > cstep)
            cap["cap_%d" % cstep] = {
                "cap": cstep, "would_solve": n, "would_solve_meet": len(meets),
                "lost": over, "p_loss_given_would_solve": over / n if n else None,
                "cp95_upper_loss": str(cp95_upper(over, n)),
                "cp95_upper_loss_f": float(cp95_upper(over, n))}
        out["arms"][k] = {
            "successes": A["succ"], "trials": n, "success_rate": A["succ"] / n,
            "cp95_success": str(cp95_lower(A["succ"], n)),
            "cp95_success_f": float(cp95_lower(A["succ"], n)),
            "coalesced_d0": A["coal"], "coalescence_rate": A["coal"] / n,
            "cp95_coalescence": str(cp95_lower(A["coal"], n)),
            "detected_given_coalescence": (A["det"] / A["coal"]) if A["coal"] else None,
            "d0_implied_success": A["det"] / n,
            "escape_fires": A["fires"], "restarts": A["restarts"],
            "trials_with_escape_fire": A["fire_trials"],
            "solve_step_median": (sorted(A["solved"])[len(A["solved"]) // 2]
                                  if A["solved"] else None),
            "solve_step_max": (max(A["solved"]) if A["solved"] else None),
        }
        out["cap_loss"][k] = cap
        out["cycle_length"][k] = {
            "trials": n, "len2": A["len2"], "len3": A["len3"], "len_ge4": A["ge4"],
            "no_repeat_within_budget": A["norep"],
            "p_len2": A["len2"] / n, "p_len3": A["len3"] / n, "p_ge4": A["ge4"] / n,
            "p_absorbed": (A["len2"] + A["len3"] + A["ge4"]) / n,
            "three_cycle_share_of_residual": (
                A["len3"] / (A["len2"] + A["len3"] + A["ge4"])
                if (A["len2"] + A["len3"] + A["ge4"]) else None)}
    E, C, L3 = out["arms"]["escaped_lag1"], out["arms"]["control_noescape"], out["arms"]["escaped_lag3"]
    out["matched_pairs"] = {
        "note": ("treatment/control share the fixed start block; prefix trace "
                 "equality asserted per pair up to the first escape fire"),
        "pairs_with_trigger": pairs["pairs_with_trigger"],
        "pairs_prefix_equal": pairs["pairs_prefix_equal"],
        "rows": pairs["rows"],
        "success_treat": E["success_rate"], "success_control": C["success_rate"],
        "delta_escape_matched": E["success_rate"] - C["success_rate"],
        "d0_treat": E["d0_implied_success"], "d5_treat": E["success_rate"],
        "delta_d0_d5": E["success_rate"] - E["d0_implied_success"],
        "lag3_success": L3["success_rate"],
        "lag3_minus_lag1": L3["success_rate"] - E["success_rate"],
        "lag3_cp95": L3["cp95_success_f"],
    }
    out["total_elapsed_s"] = round(time.time() - t_all, 1)
    os.makedirs(os.path.join(D, "evidence"), exist_ok=True)
    with open(os.path.join(D, "evidence", "census-v20n.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("[done]", out["total_elapsed_s"], "s", flush=True)
    for k in arms:
        print(k, json.dumps(out["arms"][k]), flush=True)


if __name__ == "__main__":
    main()
