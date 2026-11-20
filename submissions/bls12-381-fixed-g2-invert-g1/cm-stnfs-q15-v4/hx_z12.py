#!/usr/bin/env python3
"""
hx_z12.py — Unscreened H-R & H-X check on this package's Z[zeta_12] product-form live system
(mirroring c36p/pkg/hxgen.py): takes ALL N live columns and the first D = N + 4 live rows
UNSCREENED (zero Gaussian column/row rank pre-selection), evaluates the actual 2-step
recorded-pivot product operator op(x) = (M_live - E_live L_rec^(-1) B_rec) x with 4 merged
Schirokauer columns S = S_live - E_live L_rec^(-1) S_rec mod Q = 1,000,033, and checks at the
charged d, e: rank(C_d) == N, rank(O_e A) == rank(A), and one-block pad control < rank(A).

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


def _rank_rows(rows, q):
    R = [r[:] for r in rows if any(r)]
    k = 0
    nc = len(R[0]) if R else 0
    for c in range(nc):
        p = next((i for i in range(k, len(R)) if R[i][c]), None)
        if p is None:
            continue
        R[k], R[p] = R[p], R[k]
        iv = pow(R[k][c], q - 2, q)
        R[k] = [(x * iv) % q for x in R[k]]
        rk = R[k]
        for i in range(len(R)):
            if i != k and R[i][c]:
                f = R[i][c]
                R[i] = [(a - f * b) % q for a, b in zip(R[i], rk)]
        k += 1
    return k


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

    # Columns of unscreened M_prod and Schur-complemented merged S_prod
    Mc = [op([1 if j == k else 0 for j in range(N)]) for k in range(N)]
    z0 = [0] * N
    S = [op(z0, [1 if t == s else 0 for t in range(s_sm)]) for s in range(s_sm)]
    A = [[Mc[k][i] for k in range(N)] + [S[t][i] for t in range(s_sm)] for i in range(D)]

    # Verify op(v) == explicit Schur complement M_prod * v on a random vector
    g = random.Random(seed * 7919 + n * 31 + m)
    vt = [g.randrange(Q) for _ in range(N)]
    assert op(vt) == [sum(Mc[k][i] * vt[k] for k in range(N)) % Q for i in range(D)]

    Bv = lambda x: op(x[:N])
    d = -(-D // n) + 1 + -(-D // (128 * n))
    e = -(-D // m) + 1 + -(-D // (128 * m))
    Y = [S[t] for t in range(s_sm)] + [Bv([g.randrange(Q) for _ in range(N)]) for _ in range(n - s_sm)]
    K_rows = []
    cur = [Bv(y) for y in Y]
    for _ in range(d):
        K_rows.extend(cur)
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

    def obs(X):
        Ok = []
        Xi = [r[:] for r in X]
        for _ in range(e):
            Ok.extend(Xi)
            Xi = [[sum(r[i] * A[i][k] for i in range(D)) % Q if k < N else 0 for k in range(D)] for r in Xi]
        return _rank_rows([[sum(r[i] * A[i][k] for i in range(D)) % Q for k in range(D)] for r in Ok], Q)

    rA = _rank_rows(A, Q)
    rO = obs(Xof(False))
    rN = obs(Xof(True))
    ok = (rC == N and rO == rA and rN < rA)
    return dict(
        checks_ok=ok,
        source="unscreened_live_filtsim_z12_recorded_product",
        log2U=log2U, seed=seed,
        filt_N_sw=fres["N_sw"], filt_N_l=fres["N_l"], filt_N_rec=n_rec,
        l_rec_strict_lower_tri=True, prod_eq_schur=True,
        N=N, D=D, s_sm=s_sm, m=m, n=n, d=d, e=e,
        rank_A=rA, rank_Cd=rC, rank_OeA=rO, control_one_block=rN,
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
