#!/usr/bin/env python3
"""
hx_z12.py — On-demand H-R & H-X verification on the ACTUAL recorded-pivot product-form
operator A_prod = [M_live - E_live L_rec^(-1) B_rec | S_live - E_live L_rec^(-1) S_rec]
produced LIVE by filtsim_z12.py over Z[zeta_12] mod Q = 1,000,033 (chi^4 - chi^2 + 1 == 0),
with 4 merged Schirokauer columns (modeled generically mod Q under H-SM) and block-diagonal X
placing the 4 pad coords in distinct blocks vs 1-block pad negative control.

Every Krylov step evaluates the actual 2-step recorded-pivot product operator:
  1) forward substitution through ALL N_rec strictly lower-triangular recorded pivot rows,
  2) evaluation of the D live rows using BOTH live and back/forward-substituted recorded coords.

Usage (stdin JSON, default {"log2U":12,"seed":1,"m":12,"n":8}):
  echo '{"log2U":12,"seed":1}' | python3 -I -B hx_z12.py
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

from filtsim_z12 import CHI_HX, Q_HX, run_filtsim_z12, z12_eval_mod


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


def run_hx_on_filtered_z12(log2U=12, seed=1, m=12, n=8, s_sm=4):
    Q = Q_HX
    chi = CHI_HX
    assert (pow(chi, 4, Q) - pow(chi, 2, Q) + 1) % Q == 0
    fres = run_filtsim_z12(log2U=log2U, seed=seed, s_sm=s_sm)

    rec_order = fres["rec_order"]
    rec_col_idx = {c: k for k, c in enumerate(rec_order)}
    live_col_set = fres["live_col_set"]
    live_col_idx = {c: j for j, c in enumerate(live_col_set)}
    n_rec = len(rec_order)
    n_lcols = len(live_col_set)
    assert n_rec > 0, "Expected non-empty recorded pivot block"

    # Compile recorded pivot rows mod (Q, chi) and verify strict lower-triangularity!
    rec_compiled = []
    for k, (piv_c, piv_u, entries, sm_v) in enumerate(fres["recorded_pivots"]):
        assert rec_col_idx[piv_c] == k
        piv_q = z12_eval_mod(piv_u, chi, Q)
        assert piv_q != 0
        neg_inv_u = (-pow(piv_q, -1, Q)) % Q
        rec_terms = []
        live_terms = []
        for cc, vv in entries:
            vq = z12_eval_mod(vv, chi, Q)
            if not vq:
                continue
            if cc in rec_col_idx:
                rk = rec_col_idx[cc]
                assert rk < k, f"Upper-triangular recorded entry {rk} >= {k}!"
                rec_terms.append((rk, vq))
            else:
                live_terms.append((live_col_idx[cc], vq))
        rec_compiled.append((neg_inv_u, rec_terms, live_terms, list(sm_v)))

    # Compile ALL live rows mod (Q, chi), retaining BOTH recorded AND live column references!
    live_compiled = []
    for entries, sm_v in fres["live_rows_full"]:
        rec_terms = []
        live_terms = []
        for cc, vv in entries:
            vq = z12_eval_mod(vv, chi, Q)
            if not vq:
                continue
            if cc in rec_col_idx:
                rec_terms.append((rec_col_idx[cc], vq))
            else:
                live_terms.append((live_col_idx[cc], vq))
        if rec_terms or live_terms:
            live_compiled.append((rec_terms, live_terms, list(sm_v)))
    n_lrows = len(live_compiled)
    assert n_lrows >= s_sm + 12, f"Live system too small ({n_lrows})"

    def apply_full_product(x_live, x_sm):
        """Exact 2-step recorded-pivot product operator: forward-sub L_rec, then eval live rows."""
        x_rec = [0] * n_rec
        for k, (neg_inv_u, r_terms, l_terms, sm_v) in enumerate(rec_compiled):
            acc = 0
            for rk, vq in r_terms:
                acc += vq * x_rec[rk]
            for lc, vq in l_terms:
                acc += vq * x_live[lc]
            if x_sm is not None:
                for j in range(s_sm):
                    acc += sm_v[j] * x_sm[j]
            x_rec[k] = (neg_inv_u * (acc % Q)) % Q
        out = [0] * n_lrows
        for i, (r_terms, l_terms, sm_v) in enumerate(live_compiled):
            acc = 0
            for rk, vq in r_terms:
                acc += vq * x_rec[rk]
            for lc, vq in l_terms:
                acc += vq * x_live[lc]
            if x_sm is not None:
                for j in range(s_sm):
                    acc += sm_v[j] * x_sm[j]
            out[i] = acc % Q
        return out

    # Compute Schur-complemented S_prod = S_live - E_live L_rec^(-1) S_rec via actual product operator!
    z_live = [0] * n_lcols
    S_full = []
    for j in range(s_sm):
        e_sm = [0] * s_sm
        e_sm[j] = 1
        S_full.append(apply_full_product(z_live, e_sm))

    # Extract maximal full-rank (D x N) live sub-operator (D <= 84) WITHOUT dropping any recorded row/col!
    cap_N = min(80, n_lcols, n_lrows - s_sm)
    col_basis = {}
    sel_cols = []
    sel_col_vecs = []
    for lc in range(n_lcols):
        if len(sel_cols) >= cap_N:
            break
        e_l = [0] * n_lcols
        e_l[lc] = 1
        col_vec = apply_full_product(e_l, None)
        cur = list(col_vec)
        for i in range(n_lrows):
            if cur[i]:
                if i in col_basis:
                    b = col_basis[i]
                    f = cur[i]
                    for j in range(i, n_lrows):
                        cur[j] = (cur[j] - f * b[j]) % Q
                else:
                    inv = pow(cur[i], -1, Q)
                    col_basis[i] = [(x * inv) % Q for x in cur]
                    sel_cols.append(lc)
                    sel_col_vecs.append(col_vec)
                    break

    # Pick D = len(sel_cols) + s_sm linearly independent rows of [M_prod | S_prod]
    all_vecs = sel_col_vecs + S_full
    N = len(sel_cols)
    D = N + s_sm
    assert N >= 12 and D <= n_lrows, f"Insufficient independent columns (N={N}, D={D}, n_lrows={n_lrows})"
    row_basis = {}
    sel_rows = []
    for r_idx in range(n_lrows):
        if len(sel_rows) >= D:
            break
        cur = [vec[r_idx] for vec in all_vecs]
        for j in range(D):
            if cur[j]:
                if j in row_basis:
                    b = row_basis[j]
                    f = cur[j]
                    for t in range(j, D):
                        cur[t] = (cur[t] - f * b[t]) % Q
                else:
                    inv = pow(cur[j], -1, Q)
                    row_basis[j] = [(x * inv) % Q for x in cur]
                    sel_rows.append(r_idx)
                    break
    assert len(sel_rows) == D, f"Could not select D={D} independent rows (got {len(sel_rows)})"

    # Restrict live_compiled to sel_rows and build fast column scatter for sel_cols (ALL N_rec rows active!)
    sel_live_compiled = [live_compiled[r_idx] for r_idx in sel_rows]
    sel_col_set = {lc: idx for idx, lc in enumerate(sel_cols)}
    rec_fast = [
        (neg_inv_u, r_terms, [(sel_col_set[lc], vq) for lc, vq in l_terms if lc in sel_col_set])
        for neg_inv_u, r_terms, l_terms, _ in rec_compiled
    ]
    live_fast = [
        (r_terms, [(sel_col_set[lc], vq) for lc, vq in l_terms if lc in sel_col_set])
        for r_terms, l_terms, _ in sel_live_compiled
    ]
    S_cols = [[col[r_idx] for r_idx in sel_rows] for col in S_full]
    A_cols = [[col[r_idx] for r_idx in sel_rows] for col in sel_col_vecs] + S_cols

    def matvec_B_prod(v):
        """ACTUAL 2-step recorded-pivot product operator B(v) = (M_live - E_live L_rec^(-1) B_rec) v[:N]."""
        x_rec = [0] * n_rec
        for k, (neg_inv_u, r_terms, l_terms) in enumerate(rec_fast):
            acc = 0
            for rk, vq in r_terms:
                acc += vq * x_rec[rk]
            for j, vq in l_terms:
                acc += vq * v[j]
            x_rec[k] = (neg_inv_u * (acc % Q)) % Q
        out = [0] * D
        for i, (r_terms, l_terms) in enumerate(live_fast):
            acc = 0
            for rk, vq in r_terms:
                acc += vq * x_rec[rk]
            for j, vq in l_terms:
                acc += vq * v[j]
            out[i] = acc % Q
        return out

    # Cross-check matvec_B_prod against dense Schur complement A_cols[:N] on a random vector!
    rng = random.Random(seed * 10007 + log2U)
    v_test = [rng.randrange(1, Q) for _ in range(D)]
    y_prod = matvec_B_prod(v_test)
    y_dense = [sum(A_cols[j][i] * v_test[j] for j in range(N)) % Q for i in range(D)]
    assert y_prod == y_dense, "Product-form operator != explicit Schur complement!"

    rk_M = _rank_cols(A_cols[:N], D, Q)
    rk_A = _rank_cols(A_cols, D, Q)

    Y = [list(col) for col in S_cols] + [matvec_B_prod([rng.randrange(1, Q) for _ in range(D)]) for _ in range(n - s_sm)]
    d = math.ceil(D / n) + 1
    e = math.ceil(D / m) + 1

    cols_C = []
    cur_Y = [matvec_B_prod(y) for y in Y]
    for _ in range(d):
        cols_C.extend(cur_Y)
        cur_Y = [matvec_B_prod(y) for y in cur_Y]
    rk_Cd = _rank_cols(cols_C, D, Q)

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
                cur = matvec_B_prod(cur)
            obs_vecs.append(out)
        return _rank_cols(obs_vecs, m * e, Q)

    rk_obs = obs_rank(False)
    rk_neg = obs_rank(True)
    ok = (rk_M == N and rk_A == D and rk_Cd == N and rk_obs == D and rk_neg == D - (s_sm - 1))
    return dict(
        checks_ok=ok,
        source="live_filtsim_z12_recorded_product",
        log2U=log2U, seed=seed,
        filt_N_sw=fres["N_sw"], filt_N_l=fres["N_l"], filt_N_rec=n_rec,
        l_rec_strict_lower_tri=True, prod_eq_schur=True,
        filt_exp_unit=fres["exp_unit"],
        D=D, N=N, s_sm=s_sm, m=m, n=n, d=d, e=e,
        rank_M=rk_M, rank_A=rk_A, rank_Cd=rk_Cd, rank_OeA=rk_obs, neg_1blk_rank=rk_neg,
    )


def main():
    raw = sys.stdin.read().strip()
    req = json.loads(raw) if raw else {}
    log2U = int(req.get("log2U", 12))
    seed = int(req.get("seed", 1))
    m = int(req.get("m", 12))
    n = int(req.get("n", 8))
    if not (11 <= log2U <= 14):
        raise ValueError("log2U must be in [11, 14]")
    res = run_hx_on_filtered_z12(log2U=log2U, seed=seed, m=m, n=n)
    print(json.dumps(res, separators=(",", ":")))
    if not res["checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
