#!/usr/bin/env python3
"""
filtsim_z12.py — Pure-Python-stdlib Z[zeta_12] filtering & product-form simulator
for BLS12-381 CM-STNFS (B71.5_d04_e109: Lf=71.29, Lg=71.69, norms=548.9/677.8, excess=1.09, W=200).

Generates relations from the Dickman/Kummer prime-occurrence profile over Z[zeta_12]
(chi^4 - chi^2 + 1 = 0, 12 units), runs singleton purge, Markowitz merge, and Boudot-style
order-labeled recorded-pivot selection, returning the actual live filtered product-form matrix.

Usage (stdin JSON, default {"log2U":12,"seed":1}):
  echo '{"log2U":12,"seed":1}' | python3 -I -B filtsim_z12.py
"""
import heapq
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE
for _p in (str(ROOT), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from cm_models import LN2, primes_upto, rho_exact

# The 12 units of Z[zeta_12] = Z[chi]/(chi^4 - chi^2 + 1)
U12 = (
    (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1),
    (-1, 0, 1, 0), (0, -1, 0, 1),
    (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1),
    (1, 0, -1, 0), (0, 1, 0, -1),
)
U12_SET = set(U12)
U12_INV = {
    U12[k]: U12[(-k) % 12] for k in range(12)
}


def z12_mul(a, b):
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    c0 = a0 * b0
    c1 = a0 * b1 + a1 * b0
    c2 = a0 * b2 + a1 * b1 + a2 * b0
    c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0
    c4 = a1 * b3 + a2 * b2 + a3 * b1
    c5 = a2 * b3 + a3 * b2
    c6 = a3 * b3
    # Reduce chi^6 = -1, chi^5 = chi^3 - chi, chi^4 = chi^2 - 1
    return (c0 - c4 - c6, c1 - c5, c2 + c4, c3 + c5)


def z12_add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def z12_sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3])


def z12_neg(a):
    return (-a[0], -a[1], -a[2], -a[3])


def z12_is_zero(a):
    return a == (0, 0, 0, 0)


def z12_maxabs(a):
    return max(abs(a[0]), abs(a[1]), abs(a[2]), abs(a[3]))


def z12_eval_mod(v, chi, q):
    c2 = (chi * chi) % q
    c3 = (c2 * chi) % q
    return (v[0] + v[1] * chi + v[2] * c2 + v[3] * c3) % q


def _poisson(rng, lam):
    if lam <= 0.0:
        return 0
    if lam < 25.0:
        L = math.exp(-lam)
        p = 1.0
        k = 0
        while True:
            p *= rng.random()
            if p <= L:
                return k
            k += 1
    return max(0, round(rng.gauss(lam, math.sqrt(lam))))


def run_filtsim_z12(log2U=13, seed=1, lf=71.29, lg=71.69, nf=548.9, ng=677.8, excess=1.09, W_target=200):
    import bisect
    rng = random.Random(seed)
    U_samp = 1 << log2U
    uf, ug = nf / lf, ng / lg
    LiBf = (2.0 ** lf) / (lf * LN2)
    LiBg = (2.0 ** lg) / (lg * LN2)
    LiB_sum = LiBf + LiBg
    Up_f = int(round(U_samp * LiBf / LiB_sum))
    Up_g = U_samp - Up_f
    n_rels = int(round(excess * U_samp))
    t_dense = 0.80
    lsmall = min(t_dense * min(lf, lg), max(8.0, math.log2(max(16.0, n_rels / 25.0)) + math.log2(6.0)))
    ps_all = primes_upto(int(math.floor(2.0 ** lsmall)))
    ps_side = ps_all[::12]  # 12-orbit explicit small band per side
    n_small_tot = 2 * len(ps_side)

    rf0, rg0 = rho_exact(uf), rho_exact(ug)
    kap_f = lambda t: rho_exact(uf - t) / rf0 if rf0 > 0 else 0.0
    kap_g = lambda t: rho_exact(ug - t) / rg0 if rg0 > 0 else 0.0

    # Sampled large-band columns [0.80*L, L] matching filtsim_z12.c exactly!
    nt = 1000
    ts = [t_dense + (1.0 - t_dense) * (i + 0.5) / nt for i in range(nt)]
    cum_f, cum_g = [], []
    tf = tg = 0.0
    for t in ts:
        tf += (2.0 ** (lf * (t - 1.0))) / t
        tg += (2.0 ** (lg * (t - 1.0))) / t
        cum_f.append(tf)
        cum_g.append(tg)
    cum_f = [x / tf for x in cum_f]
    cum_g = [x / tg for x in cum_g]

    # Columns 0..n_small_tot-1 are explicit small (unpivoted); n_small_tot..n_cols-1 are pivotable large columns!
    col_lams = []
    for p in ps_side:
        col_lams.append(n_rels * min(0.99, 12.0 * kap_f(math.log2(p) / lf) / p))
        col_lams.append(n_rels * min(0.99, 12.0 * kap_g(math.log2(p) / lg) / p))
    dt_step = (1.0 - t_dense) / nt
    for _ in range(Up_f):
        idx = min(nt - 1, bisect.bisect_left(cum_f, rng.random()))
        t = ts[idx] + (rng.random() - 0.5) * dt_step
        lam = excess * LiB_sum * kap_f(t) / (2.0 ** (lf * t))
        col_lams.append(min(float(n_rels), lam))
    for _ in range(Up_g):
        idx = min(nt - 1, bisect.bisect_left(cum_g, rng.random()))
        t = ts[idx] + (rng.random() - 0.5) * dt_step
        lam = excess * LiB_sum * kap_g(t) / (2.0 ** (lg * t))
        col_lams.append(min(float(n_rels), lam))
    n_small = n_small_tot
    n_cols = len(col_lams)

    # Implicit mid-band weight on [lsmall/L, 0.80]
    w_mid_0 = 0.0
    for side_L, side_k in ((lf, kap_f), (lg, kap_g)):
        a0, b0 = lsmall / side_L, t_dense
        if b0 > a0:
            for i in range(64):
                t0 = a0 + (b0 - a0) * i / 64.0
                t1 = a0 + (b0 - a0) * (i + 1) / 64.0
                w_mid_0 += side_k(0.5 * (t0 + t1)) * math.log(t1 / t0)

    # Generate sparse relation rows with Z[zeta_12] unit coefficients
    rows = [dict() for _ in range(n_rels)]
    row_wmid = [w_mid_0] * n_rels
    for c, lam in enumerate(col_lams):
        k = _poisson(rng, lam)
        if k == 0:
            continue
        if k > n_rels:
            k = n_rels
        chosen = rng.sample(range(n_rels), k)
        for r_id in chosen:
            u = U12[rng.randrange(12)]
            prev = rows[r_id].get(c)
            if prev is None:
                rows[r_id][c] = u
            else:
                nv = z12_add(prev, u)
                if z12_is_zero(nv):
                    del rows[r_id][c]
                else:
                    rows[r_id][c] = nv

    # 1. Iterative singleton purge + clique/excess trim to PEXC = 0.0001 * U (matches filtsim_z12.c!)
    alive_r = [bool(r) for r in rows]
    col_rows = [set() for _ in range(n_cols)]
    for r_id, r in enumerate(rows):
        if alive_r[r_id]:
            for c in r:
                col_rows[c].add(r_id)

    def purge_singletons():
        q_sing = [c for c in range(n_small, n_cols) if len(col_rows[c]) == 1]
        while q_sing:
            c = q_sing.pop()
            if len(col_rows[c]) != 1:
                continue
            r_id = next(iter(col_rows[c]))
            if not alive_r[r_id]:
                continue
            alive_r[r_id] = False
            for c2 in list(rows[r_id].keys()):
                col_rows[c2].discard(r_id)
                if c2 >= n_small and len(col_rows[c2]) == 1:
                    q_sing.append(c2)

    purge_singletons()
    pexc = max(1, int(round(0.0001 * U_samp)))
    for _ in range(12):
        n_act_r = sum(alive_r)
        n_act_c = sum(1 for c in range(n_cols) if col_rows[c])
        drop = n_act_r - (n_act_c + pexc)
        if drop <= 0:
            break
        # Drop heaviest rows first (clique-weight approximation)
        cand_drop = sorted((i for i in range(n_rels) if alive_r[i]), key=lambda i: len(rows[i]), reverse=True)
        for r_id in cand_drop[: max(1, drop // 2 if drop > 4 else drop)]:
            alive_r[r_id] = False
            for c2 in list(rows[r_id].keys()):
                col_rows[c2].discard(r_id)
        purge_singletons()

    # 2. Markowitz merge on pivotable columns (c >= n_small) until mean total weight (explicit + w_mid) >= W_target
    n_alive_r = sum(alive_r)
    tot_w = sum(len(rows[i]) + row_wmid[i] for i in range(n_rels) if alive_r[i])
    pq = []
    for c in range(n_small, n_cols):
        deg = len(col_rows[c])
        if 2 <= deg <= 24:
            heapq.heappush(pq, (deg, c))

    while pq and n_alive_r > 24 and (tot_w / n_alive_r) < W_target:
        deg0, c = heapq.heappop(pq)
        deg = len(col_rows[c])
        if deg < 2 or deg != deg0 or deg > 24:
            continue
        piv = None
        piv_w = 1e18
        for r_id in col_rows[c]:
            if rows[r_id][c] in U12_SET:
                rw = len(rows[r_id]) + row_wmid[r_id]
                if rw < piv_w:
                    piv = r_id
                    piv_w = rw
        if piv is None:
            continue
        piv_row = rows[piv]
        piv_wm = row_wmid[piv]
        inv_p = U12_INV[piv_row[c]]
        others = [r_id for r_id in col_rows[c] if r_id != piv]
        ok = True
        new_rows = []
        for r_id in others:
            fac = z12_mul(rows[r_id][c], inv_p)
            nr = dict(rows[r_id])
            for cc, vv in piv_row.items():
                nv = z12_sub(nr.get(cc, (0, 0, 0, 0)), z12_mul(fac, vv))
                if z12_is_zero(nv):
                    nr.pop(cc, None)
                else:
                    if z12_maxabs(nv) > 255:
                        ok = False
                        break
                    nr[cc] = nv
            if not ok:
                break
            new_rows.append((r_id, nr))
        if not ok:
            continue
        alive_r[piv] = False
        n_alive_r -= 1
        tot_w -= (len(piv_row) + piv_wm)
        for cc in list(piv_row.keys()):
            col_rows[cc].discard(piv)
        touched = set()
        for r_id, nr in new_rows:
            old_r = rows[r_id]
            tot_w += (len(nr) - len(old_r)) + piv_wm
            row_wmid[r_id] += piv_wm
            for cc in old_r:
                if cc not in nr:
                    col_rows[cc].discard(r_id)
                    touched.add(cc)
            for cc in nr:
                if cc not in old_r:
                    col_rows[cc].add(r_id)
                    touched.add(cc)
            rows[r_id] = nr
            if not nr:
                alive_r[r_id] = False
                n_alive_r -= 1
                tot_w -= row_wmid[r_id]
        for cc in touched:
            if cc >= n_small:
                d = len(col_rows[cc])
                if 2 <= d <= 24:
                    heapq.heappush(pq, (d, cc))

    # 3. Boudot-style order-labeled product-form recorded pivots
    active_r = [i for i in range(n_rels) if alive_r[i]]
    N_sw = len(active_r)
    rec_label = {}
    row_max_lbl = {r_id: 0.0 for r_id in active_r}
    row_min_lbl = {r_id: 1.0 for r_id in active_r}
    recorded_rows = set()
    cand_cols = sorted(
        (c for c in range(n_small, n_cols) if 2 <= len(col_rows[c]) <= 64),
        key=lambda c: (len(col_rows[c]), c),
    )
    for c in cand_cols:
        c_rows = [r_id for r_id in col_rows[c] if alive_r[r_id]]
        if len(c_rows) < 2:
            continue
        hi_rec = min((row_min_lbl[r_id] for r_id in c_rows if r_id in recorded_rows), default=1.0)
        best_r = None
        best_w = 1e18
        best_mid = None
        for r_id in c_rows:
            if r_id in recorded_rows or rows[r_id][c] not in U12_SET:
                continue
            lo = row_max_lbl[r_id]
            if lo < hi_rec and (len(rows[r_id]) + row_wmid[r_id]) < best_w:
                mid = 0.5 * (lo + hi_rec)
                if lo < mid < hi_rec:
                    best_r, best_w, best_mid = r_id, len(rows[r_id]) + row_wmid[r_id], mid
        if best_r is not None:
            rec_label[c] = best_mid
            recorded_rows.add(best_r)
            for r2 in c_rows:
                if best_mid > row_max_lbl[r2]:
                    row_max_lbl[r2] = best_mid
                if best_mid < row_min_lbl[r2]:
                    row_min_lbl[r2] = best_mid

    live_rows = [r_id for r_id in active_r if r_id not in recorded_rows]
    tot_ent = sum(len(rows[r_id]) for r_id in active_r)
    unit_ent = sum(1 for r_id in active_r for v in rows[r_id].values() if v in U12_SET)
    exp_unit = unit_ent / max(tot_ent, 1)
    k5_meas = exp_unit * 7.80 + (1.0 - exp_unit) * 33.0
    mean_w = sum(len(rows[r_id]) + row_wmid[r_id] for r_id in active_r) / max(N_sw, 1)

    live_col_set = sorted({c for r_id in live_rows for c in rows[r_id] if c not in rec_label})
    col_map = {c: idx for idx, c in enumerate(live_col_set)}
    live_matrix = [
        [(col_map[c], v) for c, v in rows[r_id].items() if c in col_map]
        for r_id in live_rows
    ]
    live_matrix = [r for r in live_matrix if r]
    return dict(
        log2U=log2U, U=U_samp, seed=seed, n_rels=n_rels,
        N_sw=N_sw, N_l=len(live_rows), N_rec=len(recorded_rows),
        N_sw_U=round(N_sw / U_samp, 6), N_l_U=round(len(live_rows) / U_samp, 6),
        mean_w=round(mean_w, 2), exp_unit=round(exp_unit, 5), k5_meas=round(k5_meas, 4),
        n_live_cols=len(live_col_set),
        live_matrix=live_matrix,
    )


def main():
    raw = sys.stdin.read().strip()
    req = json.loads(raw) if raw else {}
    log2U = int(req.get("log2U", 12))
    seed = int(req.get("seed", 1))
    if not (10 <= log2U <= 15):
        raise ValueError("log2U must be in [10, 15]")
    res = run_filtsim_z12(log2U=log2U, seed=seed)
    res.pop("live_matrix")
    res["checks_ok"] = (res["N_l"] > 0 and res["N_rec"] > 0 and res["exp_unit"] > 0.50)
    print(json.dumps(res, separators=(",", ":")))
    if not res["checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
