#!/usr/bin/env python3
"""toy.py -- v19d: the v19c object with BOTH defects fixed and the invariant GREEN.

DEFECT D1 -- merge gate (fixed)
    v19c tested `report_step <= ceil(kappa*S_ref) + 2^d`: it put the DP report lag
    into the SEARCH budget.  The ledger charges the lag separately (its 2^60 term),
    so the toy must too.  The merge is now gated on kappa*S_ref alone and the lag is
    a separately measured, separately charged column.  Putting 2^d in the gate is what
    made v19c's success fall with density (0.153 @2^-5 vs 0.607 @2^-11): at 2^-11 the
    added 2048 steps are 22 S_ref of free search.

DEFECT D2 -- recovery orbit element (fixed; v19c's diagnosis retracted)
    v19c blamed the row cache.  That is wrong and is retracted here: `row()` stores
    the point TOGETHER with the (a,b) that built it, so a cache hit cannot disagree
    with the coefficient update.  An instrumented 60-trial run at 2^-5 (v19c code,
    route untouched) gives 0/3042 per-step failures and 0/... row failures.  The real
    defect is in `resolve_e`: given pS it tested x_S == beta^k * x_pS and returned
    e = lambda^k, but x_S == beta^k*x_pS means S = phi^k(pS) = [lambda^k]pS, hence
    pS = [lambda^-k]S.  It returned the inverse of the required element
    (lambda^k instead of lambda^(3-k), as lambda^3 = 1), so the recovered x was wrong
    on 23 of 34 merges (measured).  Fix: resolve e by testing the six orbit elements
    against the POINT, and verify mul(x,G) == T before a merge may count as success.
    v19c also counted an unverified merge as a success; v19d does not -- a merge whose
    recovery does not verify is a dead walk.

ARMS
    realec    our hardcoded fixture, real affine EC, production selector and cell map.
              Certifies the algebra: EVERY step asserts S == [a]G + [c]T, EVERY row
              asserts J == [rj]G + [sj]T, EVERY merge asserts the recovered x gives T.
    walkmap   the same protocol at a rival-scale class universe: the class walk is a
              deterministic pseudo-random mapping on Z_N (the random-mapping model of a
              VOW), pool = session-global over EVERY visited class, theta = 2^-d gates
              reporting only.  Used for the rival table (their N, their n, their
              statistic) and for the density sweep, because S_ref/2^d must stay >> 1;
              our fixture has S_ref = 91, so 2^-11 is degenerate there.  That is a
              property of the fixture, not of the object.
    thetacap  control object: only DP classes enter the pool.  Its prediction
              1 - exp(-kappa^2 L theta^2) collapses 4x per two dp_bits; if the port
              cannot tell the two objects apart it is not evidence.

Step, identical to the production spec and to scripts/selector.py:
    h = selector.twist_index(x, P, BETA)      # the SAME function the ledger prices
    K = min(y, P - y);  m = selector.cell(K, R, P)
    S <- S + Tab[h][m];  (a, c) <- (a + r, c + s)
"""
import argparse
import hashlib
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import selector  # noqa: E402

P = 600217
B = 1
Q = 50047
GX, GY = 218652, 2045
BETA = 43364
LAMBDA = 26745
N = (Q - 1) // 6
R = 2048
KAPPA = 1.015
KAPPA_RIVAL = 1.02
LAM = math.log(1 / 0.6097)
ORBIT_E = [1, Q - 1, LAMBDA, (Q - LAMBDA) % Q, LAMBDA * LAMBDA % Q, (Q - LAMBDA * LAMBDA) % Q]
S_REF = math.sqrt(2 * N * LAM)

# rival fixtures: sqrtN is pinned by their printed trail/sqrt(N) column (trail = 2^d),
# q34 = 23700 confirmed by dp 2, 4 and 6; q30 = 6900 confirmed by dp 4 and 6.
RIVAL_FIXTURES = {"q34": 23700.0, "q30": 6900.0, "q26": 1725.0}
RIVAL_ROWS = [("q34", 2, 12000, 0.395, 0.005), ("q34", 4, 18000, 0.405, 0.004),
              ("q34", 6, 30000, 0.404, 0.003), ("q30", 4, 24000, 0.399, 0.003),
              ("q30", 6, 24000, 0.394, 0.003)]
PIPELINE_ROWS = [("q26", 3, 18000, 0.400), ("q30", 4, 24000, 0.400),
                 ("q34", 6, 12000, 0.394)]


# ----------------------------------------------------------------- affine EC (fixture)
def add(A, C):
    if A is None:
        return C
    if C is None:
        return A
    x1, y1 = A
    x2, y2 = C
    if x1 == x2 and (y1 + y2) % P == 0:
        return None
    if x1 == x2:
        if y1 == 0:
            return None
        lam = 3 * x1 * x1 % P * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def mul(k, Pt):
    Rr = None
    while k:
        if k & 1:
            Rr = add(Rr, Pt)
        Pt = add(Pt, Pt)
        k >>= 1
    return Rr


G = (GX, GY)
BETA2 = BETA * BETA % P
assert (BETA2 + BETA + 1) % P == 0 and pow(BETA, 3, P) == 1


def key_of(S):
    """(h, K): h = production window index, K = class key |y|."""
    x, y = S
    h, _ = selector.twist_index(x, P, BETA, shift=128, n_limbs=1, band={0, 1})
    K = y if y <= P - y else P - y
    return h, K


def row(tab, h, m, T, seed):
    """Tab[h][m] = ([a]G + [b]T, a, b), one independent (a, b) per row, cached.

    The cache stores the point TOGETHER with the pair that built it: a hit cannot
    disagree with the coefficient update.  (v19c's D2 was not here.)
    """
    e = tab.get((h, m))
    if e is None:
        d = hashlib.sha256(b"v19d-row:%d:%d:%d:%d" % (seed, Q, h, m)).digest()
        a = int.from_bytes(d[:4], "big") % Q
        b = int.from_bytes(d[4:8], "big") % Q
        e = (add(mul(a, G), mul(b, T)), a, b)
        tab[(h, m)] = e
    return e


def resolve_e_by_point(pS, S):
    """e with pS == [e]S, resolved by testing the six orbit elements on the POINT.

    v19c inferred e from the x residue and the y sign and inverted it.  Testing the
    point is exact and needs no sign convention.
    """
    for e in ORBIT_E:
        if mul(e, S) == pS:
            return e
    return None


def run_trial_realec(x, d, budget, seed, nw, row_checks=True):
    """VOW, session-global pool over every visited class, DP = report gate only."""
    T = mul(x, G)
    tab = {}
    st = []
    h = hashlib.sha256(b"v19d-start:%d:%d" % (seed, x)).digest()
    for _ in range(nw):
        h = hashlib.sha256(h).digest()
        a = int.from_bytes(h[:4], "big") % Q
        c = int.from_bytes(h[4:8], "big") % Q
        st.append([add(mul(a, G), mul(c, T)), a, c, True])
    pool = {}          # K -> (a, c, S, walk_id); FIRST visitor wins
    steps = 0
    c_step = c_step_ok = c_row = c_row_ok = c_rec = c_rec_ok = 0
    merges = 0
    i = 0
    while steps < budget:
        if not any(w[3] for w in st):
            break
        wi = i % nw
        w = st[wi]
        i += 1
        if not w[3]:
            continue
        S, a, c = w[0], w[1], w[2]
        if S is None:
            w[3] = False
            continue
        hh, K = key_of(S)
        prev = pool.get(K)
        if prev is not None and prev[3] != wi:
            pa, pc, pS, _ = prev
            merges += 1
            e = resolve_e_by_point(pS, S)
            if e is not None:
                den = (pc - e * c) % Q
                if den:
                    xr = (e * a - pa) * pow(den, -1, Q) % Q
                    c_rec += 1
                    if mul(xr, G) == T:
                        c_rec_ok += 1
                        return {"report_step": steps, "merges": merges,
                                "inv_step": [c_step_ok, c_step],
                                "inv_row": [c_row_ok, c_row],
                                "inv_rec": [c_rec_ok, c_rec],
                                "pool": len(pool)}
            w[3] = False       # merge we cannot resolve: dead walk, NOT a success
            continue
        if prev is None:
            pool[K] = (a, c, S, wi)
        m = selector.cell(K, R, P)
        J, rj, sj = row(tab, hh, m, T, seed)
        if row_checks:
            c_row += 1
            if add(mul(rj, G), mul(sj, T)) == J:
                c_row_ok += 1
        S2 = add(S, J)
        steps += 1
        if S2 is None:
            w[3] = False
            continue
        w[0], w[1], w[2] = S2, (a + rj) % Q, (c + sj) % Q
        c_step += 1
        if add(mul(w[1], G), mul(w[2], T)) == w[0]:
            c_step_ok += 1
    return {"report_step": None, "merges": merges,
            "inv_step": [c_step_ok, c_step], "inv_row": [c_row_ok, c_row],
            "inv_rec": [c_rec_ok, c_rec], "pool": len(pool)}


def wilson_lower(k, n, z=1.6448536269514722):
    if n == 0:
        return 0.0
    ph = k / n
    den = 1 + z * z / n
    centre = ph + z * z / (2 * n)
    rad = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n))
    return (centre - rad) / den


# ------------------------------------------------------------ walk-map arm (large N)
MASK64 = (1 << 64) - 1


def _sm64(a):
    import numpy as np
    z = (a + np.uint64(0x9E3779B97F4A7C15)) & np.uint64(MASK64)
    z = ((z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)) & np.uint64(MASK64)
    z = ((z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)) & np.uint64(MASK64)
    return z ^ (z >> np.uint64(31))


def run_trial_walkmap(Nc, s_ref, d, kappa, seed, w_target, theta_capped=False):
    """One VOW trial on a class universe of size Nc.

    Step: c -> sm64(c ^ salt) mod Nc, a deterministic pseudo-random mapping (the
    random-mapping model of a VOW: in-degree ~ Poisson(1), which is what produces
    cross-walk merges).  Every visited class enters the session-global pool; a merge
    is a class seen by a DIFFERENT walk.  theta_capped=True pools DP classes only.
    """
    import numpy as np
    budget = int(math.ceil(kappa * s_ref))
    ell = max(2, budget // w_target)
    w = int(math.ceil(budget / ell))
    rng = np.random.default_rng(seed & 0x7FFFFFFF)
    cur = rng.integers(0, Nc, size=w, dtype=np.uint64)
    salt = np.uint64((seed * 2654435761 + 12345) & MASK64)
    rec = np.empty((ell, w), dtype=np.uint64)
    keep = np.zeros((ell, w), dtype=bool) if theta_capped else None
    thr = np.uint64(int(round((2.0 ** 64) * (2.0 ** (-d)))))
    for r in range(ell):
        rec[r] = cur
        if theta_capped:
            keep[r] = _sm64(cur ^ np.uint64(0xD1B54A32D192ED03)) < thr
        cur = _sm64(cur ^ salt) % np.uint64(Nc)
    flat = rec.reshape(-1)
    idx_all = np.arange(flat.size, dtype=np.int64)
    wid_all = idx_all % w                    # reshape(-1) is C-order: idx = r*w + walk
    stp_all = idx_all // w
    if theta_capped:
        sel = keep.reshape(-1)
        flat, wid_all, stp_all = flat[sel], wid_all[sel], stp_all[sel]
    m = int(flat.size)
    if m < 2:
        return {"t_c": None, "steps_run": budget, "pool": m, "success": False}
    srt = np.sort(flat)
    n_distinct = 1 + int(np.count_nonzero(srt[1:] != srt[:-1]))
    dup_pos = np.flatnonzero(srt[1:] == srt[:-1])
    if dup_pos.size == 0:
        return {"t_c": None, "steps_run": budget, "pool": n_distinct, "success": False}
    # first merge time: for each duplicated class, the two earliest visits by
    # DIFFERENT walks; t_c = min over those classes of the later of the two steps
    order = np.argsort(flat, kind="stable")
    cls = flat[order]
    wid = wid_all[order]
    stp = stp_all[order]
    vals = np.unique(srt[dup_pos])
    tcs = []
    for v in vals:
        a0 = int(np.searchsorted(cls, v, side="left"))
        a1 = int(np.searchsorted(cls, v, side="right"))
        first = {}
        for k in range(a0, a1):
            widx = int(wid[k])
            if widx not in first or stp[k] < first[widx]:
                first[widx] = int(stp[k])
        if len(first) > 1:
            sv = sorted(first.values())
            tcs.append(sv[1])
    t_c = min(tcs) if tcs else None
    # t_c is a LEVEL; cumulative steps walked when that level was reached = (level+1)*w
    t_c_steps = (t_c + 1) * w if t_c is not None else None
    return {"t_c": t_c_steps, "t_c_level": t_c, "steps_run": budget,
            "pool": n_distinct, "success": t_c is not None}


def one_density_walkmap(args):
    dp_bits, trials, sqrtN, kappa, w_target, arm = args
    Nc = int(round(sqrtN * sqrtN))
    s_ref = math.sqrt(2 * Nc * LAM)
    succ = 0
    det = []
    pool_max = 0
    for t in range(trials):
        r = run_trial_walkmap(Nc, s_ref, dp_bits, kappa,
                              7919 * (dp_bits + 1) + t, w_target,
                              theta_capped=(arm == "thetacap"))
        if r["success"]:
            succ += 1
            det.append(r["t_c"])
        pool_max = max(pool_max, r["pool"])
    return {"arm": arm, "dp_bits": dp_bits, "trials": trials, "sqrtN": sqrtN, "N": Nc,
            "s_ref": s_ref, "kappa": kappa, "walks_target": w_target,
            "budget_steps": int(math.ceil(kappa * s_ref)),
            "pool_records_max": pool_max, "success": succ / trials, "successes": succ,
            "cp95_lower": wilson_lower(succ, trials),
            "mean_detect_over_sref": (sum(det) / len(det) / s_ref) if det else None,
            "mean_detect_se": (max(det) - min(det)) / 2 / s_ref if len(det) > 1 else None}


def one_density_realec(args):
    d, trials, nw = args
    budget = int(math.ceil(KAPPA * S_REF))
    succ = 0
    det = []
    tot = {"step": [0, 0], "row": [0, 0], "rec": [0, 0]}
    pool_max = 0
    for t in range(trials):
        x = (int.from_bytes(hashlib.sha256(b"%d:%d" % (1000 + 100 * d, t)).digest()[:8], "big")
             % (Q - 3)) + 2
        while math.gcd(x, Q) != 1:
            x += 1
        r = run_trial_realec(x, d, budget, 1000 + 100 * d + t, nw)
        for key, name in (("inv_step", "step"), ("inv_row", "row"), ("inv_rec", "rec")):
            tot[name][0] += r[key][0]
            tot[name][1] += r[key][1]
        pool_max = max(pool_max, r["pool"])
        if r["report_step"] is not None and r["report_step"] <= budget:
            succ += 1
            det.append(r["report_step"])
    return {"arm": "realec", "dp_bits": d, "trials": trials, "walks": nw,
            "charged_steps": budget, "s_ref": S_REF, "pool_records_max": pool_max,
            "success": succ / trials, "successes": succ,
            "cp95_lower": wilson_lower(succ, trials),
            "mean_detect_over_sref": (sum(det) / len(det) / S_REF) if det else None,
            "invariant_step": tot["step"], "invariant_row": tot["row"],
            "invariant_recovery": tot["rec"],
            "invariant_ok": all(v[0] == v[1] for v in tot.values())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="all",
                    choices=["all", "realec", "walkmap", "thetacap", "rival"])
    ap.add_argument("--trials", type=int, default=1600)
    ap.add_argument("--dens", default="5,6,7,8,9,10,11")
    ap.add_argument("--walks", type=int, default=1024)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--rival-scale", type=int, default=1,
                    help="multiply the rival rows' n (their n is 1)")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                  "..", "evidence", "toy.json"))
    a = ap.parse_args()
    import multiprocessing as mp
    dens = [int(t) for t in a.dens.split(",") if t]
    out = {"package": "plain-truncwindow-40pct-v19d", "kappa": KAPPA, "s_ref": S_REF,
           "fixture": {"p": P, "b": B, "G": [GX, GY], "order": Q, "beta": BETA,
                       "classes": N, "log2_sqrtN": math.log2(math.sqrt(N)),
                       "table_rows": R},
           "protocol": {"collision_site": "every visited class (coverage S)",
                        "dp": "report gate only; lag <= one spacing, charged separately",
                        "pool": "session-global, shared by all walks",
                        "D1_fix": "merge gated on kappa*S_ref; DP lag excluded and charged separately",
                        "D2_fix": "orbit element resolved on the point; recovery verified before success"}}

    if a.arm in ("all", "realec"):
        jobs = [(d, a.trials, 32) for d in dens]
        with mp.Pool(min(a.workers, len(jobs))) as pool:
            rows = pool.map(one_density_realec, jobs)
        out["sweep_realec"] = rows
        print("realec arm (S_ref=%.2f, 32 walks):" % S_REF)
        for r in rows:
            print("  d=%2d succ %.4f CP95 %.4f step %d/%d row %d/%d rec %d/%d pool %d"
                  % (r["dp_bits"], r["success"], r["cp95_lower"],
                     r["invariant_step"][0], r["invariant_step"][1],
                     r["invariant_row"][0], r["invariant_row"][1],
                     r["invariant_recovery"][0], r["invariant_recovery"][1],
                     r["pool_records_max"]))
    if a.arm in ("all", "walkmap"):
        jobs = [(d, a.trials, RIVAL_FIXTURES["q34"], KAPPA, a.walks, "walkmap") for d in dens]
        with mp.Pool(min(a.workers, len(jobs))) as pool:
            rows = pool.map(one_density_walkmap, jobs)
        for r in rows:
            r["model"] = 1 - math.exp(-(KAPPA ** 2) * LAM)
        out["sweep_walkmap"] = rows
        out["gate"] = {"rule": "CP-95 lower >= 0.39 at every density",
                       "min_cp95_lower": min(r["cp95_lower"] for r in rows),
                       "trials_per_density": a.trials,
                       "passed": all(r["cp95_lower"] >= 0.39 for r in rows)}
        print("walkmap sweep (N=%d, walks~%d, n=%d/density):"
              % (int(23700 ** 2), a.walks, a.trials))
        for r in rows:
            print("  d=%2d succ %.4f CP95 %.4f meanT_c/S_ref %.3f pool %d"
                  % (r["dp_bits"], r["success"], r["cp95_lower"],
                     r["mean_detect_over_sref"] or -1, r["pool_records_max"]))
    if a.arm in ("all", "rival", "thetacap"):
        jobs = [(d, n * a.rival_scale, RIVAL_FIXTURES[fx], KAPPA_RIVAL, a.walks, "walkmap")
                for (fx, d, n, _p, _se) in RIVAL_ROWS]
        with mp.Pool(min(a.workers, len(jobs))) as pool:
            rows = pool.map(one_density_walkmap, jobs)
        rep = []
        for (fx, d, n0, p_rival, se_rival), r in zip(RIVAL_ROWS, rows):
            n = n0 * a.rival_scale
            diff = r["success"] - p_rival
            se = math.sqrt(se_rival ** 2 + r["success"] * (1 - r["success"]) / n)
            rep.append({"curve": fx, "dp_bits": d, "trials": n, "rival_trials": n0,
                        "ours": r["success"],
                        "z_vs_model": (r["success"] - (1 - math.exp(-(KAPPA_RIVAL ** 2) * LAM)))
                                      / math.sqrt(r["success"] * (1 - r["success"]) / n),
                        "rival": p_rival, "diff": diff, "se_combined": se,
                        "z": diff / se, "within_2se": abs(diff) <= 2 * se})
        jobs2 = [(d, n * a.rival_scale, RIVAL_FIXTURES[fx], KAPPA_RIVAL, a.walks, "walkmap")
                 for (fx, d, n, _p) in PIPELINE_ROWS]
        with mp.Pool(3) as pool:
            rows2 = pool.map(one_density_walkmap, jobs2)
        pip = [{"curve": fx, "dp_bits": d, "trials": n * a.rival_scale,
                "rival_trials": n, "ours": r["success"], "rival": p,
                "diff": r["success"] - p, "successes": r["successes"]}
               for (fx, d, n, p), r in zip(PIPELINE_ROWS, rows2)]
        pooled_n = sum(x["trials"] for x in pip)
        pooled_k = sum(x["successes"] for x in pip)
        out["rival_table"] = {
            "rows": rep, "pipeline": pip, "kappa": KAPPA_RIVAL,
            "within_2se_all": all(x["within_2se"] for x in rep),
            "flat_ours": max(x["ours"] for x in rep) - min(x["ours"] for x in rep),
            "flat_rival": max(x["rival"] for x in rep) - min(x["rival"] for x in rep),
            "pipeline_pooled": pooled_k / pooled_n,
            "pipeline_pooled_rival": 0.3988,
            "pipeline_pooled_n": pooled_n,
            "pipeline_pooled_se": math.sqrt(0.3988 * 0.6012 / pooled_n
                                            + pooled_k / pooled_n
                                            * (1 - pooled_k / pooled_n) / pooled_n),
            "model": 1 - math.exp(-(KAPPA_RIVAL ** 2) * LAM)}
        print("rival table (kappa=1.02, their N/n):")
        for x in rep:
            print("  %s d=%2d ours %.4f rival %.3f diff %+.4f z %+.2f %s"
                  % (x["curve"], x["dp_bits"], x["ours"], x["rival"], x["diff"], x["z"],
                     "OK" if x["within_2se"] else "MISS"))
        print("  pipeline pooled %.4f (rival 0.3988)" % (pooled_k / pooled_n))
        print("  flat ours %.4f rival %.4f" % (out["rival_table"]["flat_ours"],
                                               out["rival_table"]["flat_rival"]))
        if a.arm in ("all", "thetacap", "rival"):
            jobs3 = [(d, 4000, RIVAL_FIXTURES["q34"], KAPPA_RIVAL, a.walks, "thetacap")
                     for d in (2, 4, 6)]
            with mp.Pool(3) as pool:
                rows3 = pool.map(one_density_walkmap, jobs3)
            out["theta_capped_control"] = [
                {"dp_bits": r["dp_bits"], "trials": r["trials"], "successes": r["successes"],
                 "success": r["success"],
                 "naive_theta_squared_prediction":
                     1 - math.exp(-(KAPPA_RIVAL ** 2) * LAM * (2.0 ** (-2 * r["dp_bits"]))),
                 "note": ("the naive theta^2 prediction does NOT apply: once two walks "
                          "coalesce they stay together, so a DP-gated pool sees a chain of "
                          "coincidences and the cap costs only ~theta, not theta^2. The "
                          "measured collapse is the discrimination test.")}
                for r in rows3]
            print("theta-capped control (4000 trials each):",
                  [(x["dp_bits"], x["successes"], round(x["success"], 4))
                   for x in out["theta_capped_control"]])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(out, f, indent=1, sort_keys=True)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
