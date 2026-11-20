#!/usr/bin/env python3
"""
hx_z12.py — On-demand H-R & H-X verification on REAL filtered product-form Z[zeta_12]
matrices produced LIVE by filtsim_z12.py (mod Q = 1,000,033, chi^4 - chi^2 + 1 == 0,
4 dense Schirokauer columns, block-diagonal X with 4 pad coords in distinct blocks vs
1-block pad negative control).

Usage (stdin JSON, default {"log2U":13,"seed":1,"m":12,"n":8}):
  echo '{"log2U":13,"seed":1}' | python3 -I -B hx_z12.py
"""
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

from filtsim_z12 import run_filtsim_z12, z12_eval_mod


def _rank_cols(vecs, dim, q):
    basis = {}
    for v in vecs:
        cur = [x % q for x in v]
        for i in range(dim):
            if cur[i]:
                if i in basis:
                    b = basis[i]
                    f = cur[i]
                    for j in range(i, dim):
                        cur[j] = (cur[j] - f * b[j]) % q
                else:
                    inv = pow(cur[i], -1, q)
                    basis[i] = [(x * inv) % q for x in cur]
                    break
    return len(basis)


def run_hx_on_filtered_z12(log2U=13, seed=1, m=12, n=8, s_sm=4):
    Q = 1000033  # == 1 mod 12
    chi = 385966
    assert (pow(chi, 4, Q) - pow(chi, 2, Q) + 1) % Q == 0
    fres = run_filtsim_z12(log2U=log2U, seed=seed)
    raw_mat = fres["live_matrix"]
    ncols_raw = fres["n_live_cols"]

    # Evaluate Z[zeta_12] entries mod (Q, chi) and extract a full-column-rank square-plus-s_sm live system!
    rng = random.Random(seed * 10007 + log2U)
    # Build column-evaluated rows
    ev_rows = []
    for r in raw_mat:
        row_q = []
        for c, v in r:
            vq = z12_eval_mod(v, chi, Q)
            if vq:
                row_q.append((c, vq))
        if row_q:
            ev_rows.append(row_q)

    # Extract maximal linearly independent set of rows (up to cap 96 for fast Krylov rank check)
    cap_D = min(96, len(ev_rows))
    basis = {}
    indep_rows = []
    used_cols = set()
    for r in ev_rows:
        if len(indep_rows) >= cap_D:
            break
        cur = dict(r)
        while cur:
            p = min(cur)
            if p in basis:
                b = basis[p]
                f = cur[p]
                for cc, vv in b.items():
                    nv = (cur.get(cc, 0) - f * vv) % Q
                    if nv:
                        cur[cc] = nv
                    else:
                        cur.pop(cc, None)
            else:
                inv = pow(cur[p], -1, Q)
                basis[p] = {cc: (vv * inv) % Q for cc, vv in cur.items()}
                indep_rows.append(r)
                used_cols.update(c for c, _ in r)
                break

    # Pick N = len(indep_rows) - s_sm independent columns from basis pivots to form exact D x N live M!
    D = len(indep_rows)
    assert D >= s_sm + 12, f"filtered live system too small (D={D})"
    N = D - s_sm
    piv_cols = sorted(basis.keys())[:N]
    cmap = {c: i for i, c in enumerate(piv_cols)}
    M_rows = [[(cmap[c], v) for c, v in r if c in cmap] for r in indep_rows]
    # Verify M has full rank N (else top up from remaining cols if needed)
    S_cols = [[rng.randrange(1, Q) for _ in range(D)] for _ in range(s_sm)]

    def matvec_B(v):
        return [sum(val * v[c] for c, val in row) % Q for row in M_rows]

    Y = [list(col) for col in S_cols] + [matvec_B([rng.randrange(1, Q) for _ in range(D)]) for _ in range(n - s_sm)]
    d = math.ceil(D / n) + 1
    e = math.ceil(D / m) + 1

    cols_C = []
    cur_Y = [matvec_B(y) for y in Y]
    for _ in range(d):
        cols_C.extend(cur_Y)
        cur_Y = [matvec_B(y) for y in cur_Y]
    rk_Cd = _rank_cols(cols_C, D, Q)

    A_cols = []
    for c in range(N):
        col = [0] * D
        for r_idx in range(D):
            for cc, vv in M_rows[r_idx]:
                if cc == c:
                    col[r_idx] = vv
        A_cols.append(col)
    A_cols.extend(S_cols)
    rk_A = _rank_cols(A_cols, D, Q)
    rk_M = _rank_cols(A_cols[:N], D, Q)

    def obs_rank(one_block_pad):
        blocks = [[] for _ in range(m)]
        perm = list(range(N))
        rng.shuffle(perm)
        for idx, c in enumerate(perm):
            blocks[idx % m].append((c, rng.randrange(1, Q)))
        for k in range(s_sm):
            blk = 0 if one_block_pad else k
            blocks[blk].append((N + k, rng.randrange(1, Q)))
        obs_vecs = []
        for col in A_cols:
            out = []
            cur = col
            for _ in range(e):
                for blk in blocks:
                    out.append(sum(w * cur[idx] for idx, w in blk) % Q)
                cur = matvec_B(cur)
            obs_vecs.append(out)
        return _rank_cols(obs_vecs, m * e, Q)

    rk_obs = obs_rank(False)
    rk_neg = obs_rank(True)
    ok = (rk_M == N and rk_A == D and rk_Cd == N and rk_obs == D and rk_neg == D - (s_sm - 1))
    return dict(
        checks_ok=ok,
        source="live_filtsim_z12",
        log2U=log2U, seed=seed,
        filt_N_sw=fres["N_sw"], filt_N_l=fres["N_l"], filt_N_rec=fres["N_rec"],
        filt_exp_unit=fres["exp_unit"],
        D=D, N=N, s_sm=s_sm, m=m, n=n, d=d, e=e,
        rank_M=rk_M, rank_A=rk_A, rank_Cd=rk_Cd, rank_OeA=rk_obs, neg_1blk_rank=rk_neg,
    )


def main():
    raw = sys.stdin.read().strip()
    req = json.loads(raw) if raw else {}
    log2U = int(req.get("log2U", 13))
    seed = int(req.get("seed", 1))
    m = int(req.get("m", 12))
    n = int(req.get("n", 8))
    if not (12 <= log2U <= 15):
        raise ValueError("log2U must be in [12, 15]")
    res = run_hx_on_filtered_z12(log2U=log2U, seed=seed, m=m, n=n)
    print(json.dumps(res, separators=(",", ":")))
    if not res["checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
