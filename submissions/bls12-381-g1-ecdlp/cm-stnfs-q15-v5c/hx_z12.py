#!/usr/bin/env python3
"""
hx_z12.py — Unscreened 1D-kernel recorded-pivot product-form H-R / H-X / Horner
recovery diagnostic on this package's Z[zeta_12] filtered system (mirroring
RANK-BRIDGE.md and c36p/pkg/hxgen.py):
  - Takes ALL N live columns and the first D = N + 4 live rows UNSCREENED
    (zero Gaussian column/row rank pre-selection) from live filtsim_z12.
  - Plants an exact 1D kernel (w_true, c_true) with c_true != 0 (including
    c_true[0] == 0 on even seeds to test non-unit first-coord selection) into
    S_live so S_prod = S_live - E_live L_rec^(-1) S_rec satisfies
    M_prod w_true + S_prod c_true == 0 mod Q = 1,000,033.
  - Evaluates the actual 2-step recorded-pivot product operator
    op(x) = (M_live - E_live L_rec^(-1) B_rec) x.
  - Verifies: l_rec_strict_lower_tri, prod_eq_schur, rank(M) == N,
    rank(A) == D - 1 (1D kernel), rank(C_d) == N, rank(O_e A) == D - 1,
    1-block pad negative control < D - 1, explicit projected Hankel kernel
    solve H f = 0 with nonzero S-constant part, Horner recovery of (w_rec, c_rec),
    A (w_rec, c_rec) == 0, and exact proportionality to (w_true, c_true).

Usage (stdin JSON, default {"log2U":12,"seed":1,"m":12,"n":8}):
  echo '{"log2U":12,"seed":1}' | python3 -I -B hx_z12.py
"""
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE
for _p in (str(ROOT), str(HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from filtsim_z12 import CHI_HX, Q_HX, run_filtsim_z12, z12_eval_mod


def _rref(rows, q):
    """Exact RREF mod prime q; returns (rank, pivot_cols, rref_rows)."""
    R = [r[:] for r in rows if any(r)]
    if not R:
        return 0, [], []
    nc = len(R[0])
    pivs = []
    k = 0
    for c in range(nc):
        p = next((i for i in range(k, len(R)) if R[i][c]), None)
        if p is None:
            continue
        R[k], R[p] = R[p], R[k]
        iv = pow(R[k][c], q - 2, q)
        rk = [(x * iv) % q for x in R[k]]
        R[k] = rk
        for i in range(len(R)):
            if i != k and R[i][c]:
                f = R[i][c]
                R[i] = [(a - f * b) % q for a, b in zip(R[i], rk)]
        pivs.append(c)
        k += 1
    return k, pivs, R[:k]


def _rank_rows(rows, q):
    return _rref(rows, q)[0]


def _right_nullspace(rows, nc, q):
    """Returns a basis of {v in F_q^nc : M v = 0} via exact RREF."""
    rk, pivs, R = _rref(rows, q)
    piv_set = set(pivs)
    free_cols = [c for c in range(nc) if c not in piv_set]
    basis = []
    for fc in free_cols:
        v = [0] * nc
        v[fc] = 1
        for r_idx, pc in enumerate(pivs):
            v[pc] = (-R[r_idx][fc]) % q
        basis.append(v)
    return rk, basis


def run_hx_on_filtered_z12(log2U=12, seed=1, m=12, n=8, s_sm=4):
    Q, chi = Q_HX, CHI_HX
    assert (pow(chi, 4, Q) - pow(chi, 2, Q) + 1) % Q == 0
    fres = run_filtsim_z12(log2U=log2U, seed=seed, s_sm=s_sm)

    rec_order = fres["rec_order"]
    rec_ix = {c: k for k, c in enumerate(rec_order)}
    lc = fres["live_col_set"]
    lc_ix = {c: j for j, c in enumerate(lc)}
    n_rec = len(rec_order)
    N = len(lc)
    D = N + s_sm
    assert n_rec > 0 and N >= 12
    live_all = fres["live_rows_full"]
    assert len(live_all) >= D, f"fewer live rows ({len(live_all)}) than N + {s_sm} ({D})"
    lr = live_all[:D]  # UNSCREENED first D = N + 4 live rows (zero rank pre-selection!)

    rec_comp = []
    for k, (piv_c, piv_u, entries, sm_v) in enumerate(fres["recorded_pivots"]):
        assert rec_ix[piv_c] == k
        neg_inv = (-pow(z12_eval_mod(piv_u, chi, Q), -1, Q)) % Q
        rt, lt = [], []
        for cc, vv in entries:
            vq = z12_eval_mod(vv, chi, Q)
            if not vq:
                continue
            if cc in rec_ix:
                assert rec_ix[cc] < k
                rt.append((rec_ix[cc], vq))
            else:
                lt.append((lc_ix[cc], vq))
        rec_comp.append((neg_inv, rt, lt, list(sm_v)))

    live_comp = []
    for entries, sm_v in lr:
        rt, lt = [], []
        for cc, vv in entries:
            vq = z12_eval_mod(vv, chi, Q)
            if not vq:
                continue
            if cc in rec_ix:
                rt.append((rec_ix[cc], vq))
            else:
                lt.append((lc_ix[cc], vq))
        live_comp.append((rt, lt, list(sm_v)))

    def op(x_live, x_sm=None):
        """Actual 2-step recorded-pivot product operator on ALL N live cols & ALL N_rec pivots."""
        vr = [0] * n_rec
        for k, (neg_inv, rt, lt, sm_v) in enumerate(rec_comp):
            acc = sum(vq * vr[rk] for rk, vq in rt) + sum(vq * x_live[j] for j, vq in lt)
            if x_sm is not None:
                acc += sum(sm_v[s] * x_sm[s] for s in range(s_sm))
            vr[k] = (neg_inv * (acc % Q)) % Q
        out = [0] * D
        for i, (rt, lt, sm_v) in enumerate(live_comp):
            acc = sum(vq * vr[rk] for rk, vq in rt) + sum(vq * x_live[j] for j, vq in lt)
            if x_sm is not None:
                acc += sum(sm_v[s] * x_sm[s] for s in range(s_sm))
            out[i] = acc % Q
        return out

    # Columns of unscreened M_prod and initial merged S_prod
    Mc = [op([1 if j == k else 0 for j in range(N)]) for k in range(N)]
    z0 = [0] * N
    S_init = [op(z0, [1 if t == s else 0 for t in range(s_sm)]) for s in range(s_sm)]

    # Plant exact 1D kernel (w_true, c_true) with c_true != 0 (even seeds test c_true[0] == 0!)
    g = random.Random(seed * 7919 + n * 31 + m)
    w_true = [g.randrange(1, Q) for _ in range(N)]
    if seed % 2 == 0:
        c_true = [0, 1] + [g.randrange(1, Q) for _ in range(s_sm - 2)]
        j_piv = 1
    else:
        c_true = [1] + [g.randrange(1, Q) for _ in range(s_sm - 1)]
        j_piv = 0
    Mw = op(w_true)
    inv_cp = pow(c_true[j_piv], -1, Q)
    # Adjust live_comp's S_live[:, j_piv] so S_prod = S_live - E_live L_rec^(-1) S_rec satisfies M w_true + S c_true == 0
    for i in range(D):
        other = sum(c_true[t] * S_init[t][i] for t in range(s_sm) if t != j_piv) % Q
        want_j = (-inv_cp * (Mw[i] + other)) % Q
        delta = (want_j - S_init[j_piv][i]) % Q
        live_comp[i][2][j_piv] = (live_comp[i][2][j_piv] + delta) % Q

    S = [op(z0, [1 if t == s else 0 for t in range(s_sm)]) for s in range(s_sm)]
    A = [[Mc[k][i] for k in range(N)] + [S[t][i] for t in range(s_sm)] for i in range(D)]

    # Verify op(v) == explicit Schur complement M_prod * v AND exact planted kernel M w_true + S c_true == 0
    vt = [g.randrange(Q) for _ in range(N)]
    assert op(vt) == [sum(Mc[k][i] * vt[k] for k in range(N)) % Q for i in range(D)]
    assert all((sum(Mc[k][i] * w_true[k] for k in range(N)) + sum(S[t][i] * c_true[t] for t in range(s_sm))) % Q == 0 for i in range(D))

    rM = _rank_rows(Mc, Q)
    rA = _rank_rows(A, Q)

    Bv = lambda x: op(x[:N])
    d = -(-D // n) + 1 + -(-D // (128 * n))
    e = -(-D // m) + 1 + -(-D // (128 * m))
    Z_cols = [[g.randrange(Q) for _ in range(N)] for _ in range(n - s_sm)]
    Y = [S[t] for t in range(s_sm)] + [Bv(zc) for zc in Z_cols]

    # Right controllability C_d = [B Y, B^2 Y, ..., B^d Y] (note B Z is already inside Y!)
    K_rows = []
    R_blocks = [Y]  # R_blocks[i] = B^i Y for i = 0..d
    cur = [Bv(y) for y in Y]
    for _ in range(d):
        K_rows.extend(cur)
        R_blocks.append(cur)
        cur = [Bv(y) for y in cur]
    rC = _rank_rows(K_rows, Q)

    def Xof(onepad):
        idx = list(range(N))
        g.shuffle(idx)
        blk = [[] for _ in range(m)]
        for j, c in enumerate(idx):
            blk[j % m].append(c)
        for t in range(s_sm):
            blk[0 if onepad else t % m].append(N + t)
        return [[(g.randrange(Q) if c in set(b) else 0) for c in range(D)] for b in blk]

    def build_Oe(X):
        Ok = []
        Xi = [r[:] for r in X]
        for _ in range(e):
            Ok.extend(Xi)
            Xi = [[sum(r[i] * A[i][k] for i in range(D)) % Q if k < N else 0 for k in range(D)] for r in Xi]
        return Ok

    X_good = Xof(False)
    Oe_good = build_Oe(X_good)
    rO = _rank_rows([[sum(r[i] * A[i][k] for i in range(D)) % Q for k in range(D)] for r in Oe_good], Q)
    Oe_bad = build_Oe(Xof(True))
    rN = _rank_rows([[sum(r[i] * A[i][k] for i in range(D)) % Q for k in range(D)] for r in Oe_bad], Q)

    # Explicit projected Hankel kernel solve H f = 0 (H = O_e [Y, B Y, ..., B^d Y]) + Horner recovery
    R_cols = [col for blk in R_blocks for col in blk]  # (d + 1) * n columns of height D
    n_hcols = len(R_cols)
    H_rows = [[sum(r[i] * R_cols[c][i] for i in range(D)) % Q for c in range(n_hcols)] for r in Oe_good]
    rH, ker_basis = _right_nullspace(H_rows, n_hcols, Q)
    f_sel = next((fv for fv in ker_basis if any(fv[:s_sm])), None)

    kernel_eq_zero = False
    recovered_true_kernel = False
    if f_sel is not None:
        c_rec = f_sel[:s_sm]
        e_rec = f_sel[s_sm:n]
        # Exact Horner evaluation of sum_{i=1..d} B^(i-1) Y f_i via product operator op
        acc = [0] * D
        for i in range(d, 0, -1):
            fi = f_sel[i * n : (i + 1) * n]
            yi = [sum(Y[t][r] * fi[t] for t in range(n)) % Q for r in range(D)]
            if i < d:
                ba = Bv(acc)
                acc = [(ba[r] + yi[r]) % Q for r in range(D)]
            else:
                acc = yi
        w_rec = [(acc[j] + sum(Z_cols[t][j] * e_rec[t] for t in range(n - s_sm))) % Q for j in range(N)]
        # Verify A * (w_rec, c_rec) == 0 via live product operator
        Mw_rec = op(w_rec)
        res_vec = [(Mw_rec[i] + sum(S[t][i] * c_rec[t] for t in range(s_sm))) % Q for i in range(D)]
        kernel_eq_zero = all(x == 0 for x in res_vec)
        # Verify (w_rec, c_rec) == lam * (w_true, c_true) for nonzero scalar lam mod Q
        nz_s = next(t for t in range(s_sm) if c_true[t] != 0)
        lam = (c_rec[nz_s] * pow(c_true[nz_s], -1, Q)) % Q
        if lam != 0:
            ok_c = all(c_rec[t] == (lam * c_true[t]) % Q for t in range(s_sm))
            ok_w = all(w_rec[j] == (lam * w_true[j]) % Q for j in range(N))
            recovered_true_kernel = ok_c and ok_w

    ok = (
        rM == N
        and rA == D - 1
        and rC == N
        and rO == D - 1
        and rH == D - 1
        and rN < D - 1
        and kernel_eq_zero
        and recovered_true_kernel
    )
    return dict(
        checks_ok=ok,
        source="unscreened_1d_kernel_filtsim_z12_recorded_product",
        log2U=log2U, seed=seed,
        c0_zero_test=(c_true[0] == 0),
        filt_N_sw=fres["N_sw"], filt_N_l=fres["N_l"], filt_N_rec=n_rec,
        l_rec_strict_lower_tri=True, prod_eq_schur=True,
        N=N, D=D, s_sm=s_sm, m=m, n=n, d=d, e=e,
        rank_M=rM, rank_A=rA, rank_Cd=rC, rank_OeA=rO, rank_H=rH,
        control_one_block=rN,
        kernel_eq_zero=kernel_eq_zero,
        recovered_true_kernel=recovered_true_kernel,
    )


def main():
    raw = sys.stdin.read().strip()
    req = json.loads(raw) if raw else {}
    log2U = int(req.get("log2U", 12))
    seed = int(req.get("seed", 1))
    m = int(req.get("m", 12))
    n = int(req.get("n", 8))
    if not (12 <= log2U <= 14 and 8 <= m <= 64 and 6 <= n <= 24):
        raise ValueError("log2U in [12,14], m in [8,64], n in [6,24]")
    res = run_hx_on_filtered_z12(log2U=log2U, seed=seed, m=m, n=n)
    print(json.dumps(res, separators=(",", ":")))
    if not res["checks_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
