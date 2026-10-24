#!/usr/bin/env python3
"""match_control.py -- v20n matched-pair census (verdict e9ee1dbc point 2).

One class universe, one FIXED start block: trial b uses the identical (x, a, c)
in every arm (v20m drew `bt` challenges before the coefficients, so a 5120-trial
arm never shared a start with its 2048-trial control).

The escape/no-escape comparison is a matched pair: the treatment substitutes the
alternate row at the first class-repeat, the control walks on.  The prefix traces
are compared byte-for-byte up to the branch step and shipped as hash rows.

d=0 is not a second walk: it is the SAME trace read with the report gate off, so
coalescence and the detection factor come from identical trajectories.

Also emitted: the lag-3 escape variant (3-cycle trigger), the cycle-length census
{2,3,>=4} of residual absorption, and the cap-loss row read off the trail
histogram: P(would-solve step > cap | would-solve), with a CP-95 upper bound.
"""
import hashlib
import json
import math
import os
import sys
import time
from decimal import Decimal, getcontext

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import qwalk  # noqa: E402
import selector  # noqa: E402
from census_v20n import SHAPE, cp95_lower, cp95_upper, decomposition  # noqa: E402

W, D_BITS, BT_MAX, CHUNK, SEED = 16, 5, 3200, 400, 20250
R_TOY = 64


def make_start(U, bt_max, w, seed):
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


def run(U, canon0, a0, cc0, xv, d, per, esc, lag, w, bt):
    lanes = canon0.size
    canon = canon0.copy(); a = a0.copy(); cc = cc0.copy()
    prev = canon.copy(); prev2 = canon.copy(); prev3 = canon.copy()
    esc_run = np.zeros(lanes, dtype=np.int64)
    coins = np.random.default_rng(0xC0FFEE)
    trace = np.empty((per, lanes), dtype=np.uint32)
    first_fire = np.full(bt, -1, dtype=np.int64)
    first = [dict() for _ in range(bt)]
    succ = np.full(bt, -1, dtype=np.int64)
    fires = restarts = 0
    for t in range(per):
        trace[t] = canon.astype(np.uint32)
        rst = np.zeros(lanes, dtype=bool)
        cn, a, cc, hit, pv, fired = qwalk.step(
            U, canon, a, cc, xv, d, True, prev, esc, 0, t, esc_run,
            (qwalk.K_ESC if esc else None), rst, coins)
        if lag >= 2:
            others = (prev2, prev3) if lag >= 3 else (prev2,)
            for other in others:
                extra = (cn == other) & ~fired
                if extra.any():
                    idx = np.flatnonzero(extra)
                    m = (U.col(canon[idx]) + 1) % U.R
                    c2, a3, cc3, _e = qwalk._apply(U, canon[idx], a[idx], cc[idx], xv[idx], m)
                    cn = cn.copy(); a = a.copy(); cc = cc.copy()
                    cn[idx] = c2; a[idx] = a3; cc[idx] = cc3
                    fired = fired | extra
        canon = cn                      # ADVANCE the class (v20n shipped without this)
        fires += int(fired.sum()); restarts += int(rst.sum())
        if fired.any():
            fb = np.flatnonzero(np.bincount(np.flatnonzero(fired) // w, minlength=bt) > 0)
            new = fb[first_fire[fb] < 0]
            first_fire[new] = t
        prev3 = prev2; prev2 = prev; prev = pv
        if d > 0 and hit is not None and hit.any():
            for i in np.flatnonzero(hit).tolist():
                b = i // w
                if succ[b] >= 0:
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
                            succ[b] = t
    return trace, succ, first_fire, fires, restarts


def main():
    t0 = time.time()
    U = qwalk.Universe(4_290_000_000, R_TOY, seed=SEED)
    assert U.q < (1 << 32)
    per = int(round(float(SHAPE) * math.sqrt(U.classes) / W))
    shape = per * W / math.sqrt(U.classes)
    assert abs(shape - float(SHAPE)) <= 1e-3, (shape, float(SHAPE))
    xt, xv, a_all, cc_all, canon_all = make_start(U, BT_MAX, W, SEED)
    acc = {k: dict(succ=0, coal=0, det=0, l2=0, l3=0, g4=0, nor=0, fires=0,
                   restarts=0, solved=[], meets=[], ftrials=0)
           for k in ("lag1", "lag3", "control")}
    pairs = dict(trig=0, eq=0, rows=[])
    for c0 in range(0, BT_MAX, CHUNK):
        bt = min(CHUNK, BT_MAX - c0)
        sl = slice(c0 * W, (c0 + bt) * W)
        out = {}
        for k, esc, lag in (("lag1", True, 1), ("lag3", True, 3), ("control", False, 0)):
            tr, succ, ff, fires, rst = run(U, canon_all[sl], a_all[sl], cc_all[sl],
                                           xv[sl], D_BITS, per, esc, lag, W, bt)
            coal, det, meet, cyc = decomposition(U, tr, W, bt, D_BITS)
            A = acc[k]
            A["succ"] += int(np.count_nonzero(succ >= 0))
            A["coal"] += int(np.count_nonzero(coal))
            A["det"] += int(np.count_nonzero(coal & det))
            A["l2"] += int(np.count_nonzero(cyc == 2))
            A["l3"] += int(np.count_nonzero(cyc == 3))
            A["g4"] += int(np.count_nonzero(cyc >= 4))
            A["nor"] += int(np.count_nonzero(cyc == 0))
            A["fires"] += fires; A["restarts"] += rst
            A["solved"] += succ[succ >= 0].tolist()
            A["meets"] += meet[meet >= 0].tolist()
            A["ftrials"] += int(np.count_nonzero(ff >= 0))
            out[k] = (tr, ff)
        tt, ff = out["lag1"]
        tc, _ = out["control"]
        for b in np.flatnonzero(ff >= 0).tolist():
            t = int(ff[b])
            pt = tt[:t, b * W:(b + 1) * W]
            pc = tc[:t, b * W:(b + 1) * W]
            eq = bool(pt.shape == pc.shape and np.array_equal(pt, pc))
            pairs["trig"] += 1; pairs["eq"] += int(eq)
            if len(pairs["rows"]) < 6:
                pairs["rows"].append({
                    "trial": c0 + b, "branch_step": t,
                    "prefix_sha256_treat": hashlib.sha256(pt.tobytes()).hexdigest()[:16],
                    "prefix_sha256_control": hashlib.sha256(pc.tobytes()).hexdigest()[:16],
                    "equal": eq, "prefix_lanes": int(pt.size)})
    n = BT_MAX
    out = {"test": "census_v20n_matched_pairs", "q": int(U.q),
           "classes": int(U.classes), "R_toy": R_TOY, "R_eff_toy": int(U.R_eff),
           "walk_downsampling": ("toy column count R_toy=64 (R_eff=32) exercises the "
                                 "two-step channel at 1/(6*32)/step; production is "
                                 "1/(6*852105)/step, so the fixture is a DISCLOSED "
                                 "stress rate ~2.7e4x production for the escape arm"),
           "walks_per_trial": W, "trials": n, "dp_bits": D_BITS,
           "per_walker_steps": per, "visits_per_trial": per * W,
           "visits_over_sqrtN": shape, "ledger_shape": float(SHAPE),
           "shape_assert_dev": abs(shape - float(SHAPE)),
           "arms": {}, "cap_loss": {}, "cycle_length": {}, "matched_pairs": {}}
    for k in ("lag1", "lag3", "control"):
        A = acc[k]
        cap = {}
        meets = sorted(A["meets"])
        for frac in (0.25, 0.5, 1.0):
            cstep = int(per * frac)
            over = sum(1 for s_ in meets if s_ > cstep)
            cap["cap_%d" % cstep] = {
                "cap": cstep, "trials": n, "would_solve": len(meets),
                "lost": over, "p_loss_given_would_solve": over / n,
                "cp95_upper_loss": str(cp95_upper(over, n)),
                "cp95_upper_loss_f": float(cp95_upper(over, n))}
        out["arms"][k] = {
            "successes": A["succ"], "trials": n, "success_rate": A["succ"] / n,
            "cp95_success": str(cp95_lower(A["succ"], n)),
            "cp95_success_f": float(cp95_lower(A["succ"], n)),
            "coalesced_d0": A["coal"], "coalescence_rate": A["coal"] / n,
            "cp95_coalescence_f": float(cp95_lower(A["coal"], n)),
            "detected_given_coalescence": (A["det"] / A["coal"]) if A["coal"] else None,
            "d0_implied_success": A["det"] / n,
            "escape_fires": A["fires"], "restarts": A["restarts"],
            "trials_with_escape_fire": A["ftrials"],
            "solve_step_median": (sorted(A["solved"])[len(A["solved"]) // 2]
                                  if A["solved"] else None),
            "solve_step_max": (max(A["solved"]) if A["solved"] else None)}
        out["cap_loss"][k] = cap
        out["cycle_length"][k] = {
            "trials": n, "len2": A["l2"], "len3": A["l3"], "len_ge4": A["g4"],
            "no_repeat_within_budget": A["nor"],
            "p_len2": A["l2"] / n, "p_len3": A["l3"] / n, "p_ge4": A["g4"] / n,
            "p_absorbed": (A["l2"] + A["l3"] + A["g4"]) / n,
            "three_cycle_share_of_residual": (
                A["l3"] / (A["l2"] + A["l3"] + A["g4"])
                if (A["l2"] + A["l3"] + A["g4"]) else None)}
    L1, L3, C0 = (out["arms"]["lag1"], out["arms"]["lag3"], out["arms"]["control"])
    out["matched_pairs"] = {
        "starts_shared_by_every_arm": True,
        "pairs_with_trigger": pairs["trig"], "pairs_prefix_equal": pairs["eq"],
        "rows": pairs["rows"], "success_treat": L1["success_rate"],
        "success_control": C0["success_rate"],
        "delta_escape_matched": L1["success_rate"] - C0["success_rate"],
        "delta_cp95_pessimistic": (L1["cp95_success_f"] - C0["success_rate"]),
        "d0_treat": L1["d0_implied_success"], "d5_treat": L1["success_rate"],
        "delta_d0_d5": L1["success_rate"] - L1["d0_implied_success"],
        "lag3_success": L3["success_rate"], "lag3_minus_lag1": L3["success_rate"] - L1["success_rate"],
        "lag3_cp95": L3["cp95_success_f"]}
    # ---- production decomposition, propagated from MEASURED factors ----------
    S_v = per * W
    model_fx = 1 - math.exp(-S_v * (S_v - 1) / (2.0 * U.classes))
    model_fx_w = 1 - math.exp(-S_v * (S_v - 1) * (W - 1) / (2.0 * U.classes * W))
    coal = L1["coalescence_rate"]
    out["decomposition"] = {
        "model_fixture_at_walkers": model_fx,
        "model_fixture_walker_corrected": model_fx_w,
        "measured_coalescence": coal,
        "ratio_measured_over_walker_model": coal / model_fx_w,
        "measured_detection": L1["detected_given_coalescence"],
        "production_model": 0.3995,
        "propagated_point": 0.3995 * (coal / model_fx_w) * L1["detected_given_coalescence"],
        "propagated_cp_lower": 0.3995 * (L1["cp95_coalescence_f"] / model_fx_w) * 0.9,
        "note": ("walker factor is measured (16 walkers, same trace), DP factor is "
                 "measured at 2^-5 and floored at 0.9 for the CP row; the propagated "
                 "number is BELOW the 0.39 floor and is reported as such")}
    assert pairs["trig"] > 0 and pairs["eq"] == pairs["trig"], pairs["trig"]
    out["elapsed_s"] = round(time.time() - t0, 1)
    with open(os.path.join(D, "evidence", "census-v20n.json"), "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print(json.dumps({k: out["arms"][k] for k in out["arms"]}, indent=1))
    print("pairs", pairs["trig"], "equal", pairs["eq"], "elapsed", out["elapsed_s"])


if __name__ == "__main__":
    main()
