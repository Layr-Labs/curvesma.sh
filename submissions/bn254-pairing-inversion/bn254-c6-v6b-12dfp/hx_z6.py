#!/usr/bin/env python3
"""
hx_z6.py — Exact 1D-Kernel H-R / H-X / Projected-Hankel Recovery Check on Z[omega] Recorded-Pivot System:
 1. Runs the live pure-Python Z[omega] filtering simulator at (Lf=61, Lg=60, exc=0.82, W=160) to obtain
  the real recorded-pivot lower-triangular L_rec, B_rec, E_live, M_live (UNSCREENED: zero rank pre-selection).
 2. Constructs the Schur-complement product operator op(x) = (M_live - E_live L_rec^(-1) B_rec) x mod Q = 1,000,033
  (with chi^2 - chi + 1 == 0 mod Q) on N live columns and D = N + s_sm rows (s_sm = 5, matching BN254 C6).
 3. Plants an EXACT 1D right kernel (w_true, c_true) (c_true != 0) into S_prod so A = [M_prod | S_prod] satisfies
  A * [w_true; c_true] == 0 mod Q and rank(A) == D - 1 (exact H-R 1D kernel!).
 4. Checks H-X at the charged Block-Wiedemann depths d, e on sparse left projection X:
   - Controllability: rank(C_d) == N
   - Observability:   rank(O_e A) == rank(A) == D - 1
   - Negative control (one-block Schirokauer pad): rank(O_e^{one} A) < rank(A)
 5. Forms the projected Block-Wiedemann Hankel matrix H = [X B^{i+j} Y] (0 <= i < e, 0 <= j <= d),
  solves H f = 0 mod Q, maps f -> (w(f), c(f)), and verifies EXACT recovery (w(f), c(f)) == alpha * (w_true, c_true)
  with alpha != 0 mod Q!
Usage (stdin JSON, default {"log2U":12,"seed":1,"m":12,"n":8}):
 echo '{"log2U":12,"seed":1}' | python3 -I -B hx_z6.py
"""
import json, random, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE
for _p in (str(ROOT), str(HERE)):
  if _p not in sys.path:
    sys.path.insert(0, _p)
import verify_bn254_v6 as V
Q_HX = 1000033  # Prime with 6 | (Q_HX - 1)
CHI_HX = next(x for g in range(2, 50) for x in [pow(g, (Q_HX - 1) // 6, Q_HX)] if (x * x - x + 1) % Q_HX == 0)
def _z6_eval(v, chi=CHI_HX, q=Q_HX):
  return (v[0] + v[1] * chi) % q
def _rref_and_nullspace(rows, ncols, q):
  """Exact Gauss-Jordan RREF over F_q returning (rank, basis of right nullspace)."""
  R = [r[:] for r in rows]
  nrows = len(R)
  piv_cols = []
  r = 0
  for c in range(ncols):
    piv = next((i for i in range(r, nrows) if R[i][c]), None)
    if piv is None:
      continue
    R[r], R[piv] = R[piv], R[r]
    iv = pow(R[r][c], q - 2, q)
    rr = [(x * iv) % q for x in R[r]]
    R[r] = rr
    for i in range(nrows):
      if i != r and R[i][c]:
        f = R[i][c]
        ri = R[i]
        for j in range(c, ncols):
          if rr[j]:
            ri[j] = (ri[j] - f * rr[j]) % q
    piv_cols.append(c)
    r += 1
    if r == nrows:
      break
  piv_set = set(piv_cols)
  free_cols = [c for c in range(ncols) if c not in piv_set]
  basis = []
  for fc in free_cols:
    vec = [0] * ncols
    vec[fc] = 1
    for row_i, pc in enumerate(piv_cols):
      if R[row_i][fc]:
        vec[pc] = (-R[row_i][fc]) % q
    basis.append(vec)
  return r, basis
def _rank_rows(rows, q):
  if not rows:
    return 0
  return _rref_and_nullspace(rows, len(rows[0]), q)[0]
def run_hx_1d_kernel_z6(log2U=12, seed=1, m=12, n=8, s_sm=5):
  Q, chi = Q_HX, CHI_HX
  assert (chi * chi - chi + 1) % Q == 0 and (pow(chi, 3, Q) + 1) % Q == 0
  fres = V.verify_live_filtsim_z6(log2U=log2U, seed=seed, W_target=160, check_kernel=False, export_hx=True)
  hx = fres["hx_export"]
  rec_order = hx["rec_order"]
  rec_ix = {c: k for k, c in enumerate(rec_order)}
  lc = hx["live_col_set"]
  lc_ix = {c: j for j, c in enumerate(lc)}
  n_rec = len(rec_order)
  N = len(lc)
  D = N + s_sm
  assert n_rec > 0 and N >= 8, f"too small: N={N}, n_rec={n_rec}"
  live_all = hx["live_rows_full"]
  assert len(live_all) >= D, f"fewer live rows ({len(live_all)}) than N+{s_sm} ({D})"
  lr = live_all[:D]  # UNSCREENED first D = N + s_sm live rows!
  rec_comp = []
  for k, (piv_c, piv_u, entries) in enumerate(hx["recorded_pivots"]):
    assert rec_ix[piv_c] == k
    neg_inv = (-pow(_z6_eval(piv_u), -1, Q)) % Q
    rt, lt = [], []
    for cc, vv in entries:
      vq = _z6_eval(vv)
      if not vq: continue
      if cc in rec_ix:
        assert rec_ix[cc] < k  # strictly lower triangular L_rec!
        rt.append((rec_ix[cc], vq))
      else:
        lt.append((lc_ix[cc], vq))
    rec_comp.append((neg_inv, rt, lt))
  live_comp = []
  for entries in lr:
    rt, lt = [], []
    for cc, vv in entries:
      vq = _z6_eval(vv)
      if not vq: continue
      if cc in rec_ix: rt.append((rec_ix[cc], vq))
      else: lt.append((lc_ix[cc], vq))
    live_comp.append((rt, lt))
  def op(x_live):
    vr = [0] * n_rec
    for k, (neg_inv, rt, lt) in enumerate(rec_comp):
      acc = sum(vq * vr[rk] for rk, vq in rt) + sum(vq * x_live[j] for j, vq in lt)
      vr[k] = (neg_inv * (acc % Q)) % Q
    return [(sum(vq * vr[rk] for rk, vq in rt) + sum(vq * x_live[j] for j, vq in lt)) % Q
        for rt, lt in live_comp]
  # Explicit columns of unscreened Schur-complement operator M_prod (D x N)
  Mc = [op([1 if j == k else 0 for j in range(N)]) for k in range(N)]
  g = random.Random(seed * 7919 + n * 31 + m)
  vt = [g.randrange(Q) for _ in range(N)]
  prod_eq_schur = (op(vt) == [sum(Mc[k][i] * vt[k] for k in range(N)) % Q for i in range(D)])
  assert prod_eq_schur
  # Construct S_prod (D x s_sm) with an EXACT 1D right kernel (w_true, c_true), c_true[0] = 1 != 0!
  w_true = [g.randrange(1, Q) for _ in range(N)]
  c_true = [1] + [g.randrange(1, Q) for _ in range(s_sm - 1)]
  Mw = op(w_true)
  S_cols = [None] + [[g.randrange(Q) for _ in range(D)] for _ in range(s_sm - 1)]
  S_cols[0] = [(-(Mw[i] + sum(c_true[t] * S_cols[t][i] for t in range(1, s_sm)))) % Q for i in range(D)]
  A = [[Mc[k][i] for k in range(N)] + [S_cols[t][i] for t in range(s_sm)] for i in range(D)]
  # Verify H-R (exact rank(M_prod) == N and exact 1D kernel of A!)
  rM = _rank_rows([[Mc[k][i] for k in range(N)] for i in range(D)], Q)
  rA, kerA = _rref_and_nullspace(A, D, Q)
  v_true = w_true + c_true
  hr_1d_ok = (rM == N and rA == D - 1 and len(kerA) == 1)
  if hr_1d_ok:
    iv_k = pow(kerA[0][N], Q - 2, Q)
    hr_1d_ok = ([(x * iv_k) % Q for x in kerA[0]] == v_true)
  # Verify H-X controllability & observability at charged depths d, e on sparse X
  Bv = lambda x: op(x[:N])
  d = -(-D // n) + 1 + -(-D // (128 * n))
  e = -(-D // m) + 1 + -(-D // (128 * m))
  U_rand = [[g.randrange(Q) for _ in range(N)] for _ in range(n - s_sm)]
  Y = [S_cols[t] for t in range(s_sm)] + [op(u) for u in U_rand]  # D x n
  # Krylov columns B^j Y for j = 0..d+e-1
  BY = [Y]
  for _ in range(d + e - 1):
    BY.append([Bv(col) for col in BY[-1]])
  # Controllability C_d = [B Y, B^2 Y, ..., B^d Y]
  K_rows = [col for j in range(1, d + 1) for col in BY[j]]
  rC = _rank_rows(K_rows, Q)
  def make_X(onepad):
    idx = list(range(N)); g.shuffle(idx)
    blk = [[] for _ in range(m)]
    for j, c in enumerate(idx): blk[j % m].append(c)
    for t in range(s_sm): blk[0 if onepad else t % m].append(N + t)
    return [[(g.randrange(1, Q) if c in set(b) else 0) for c in range(D)] for b in blk]
  X_good = make_X(False)
  X_bad = make_X(True)
  def obs_rank(X):
    Ok = []; Xi = [r[:] for r in X]
    for _ in range(e):
      Ok.extend(Xi)
      Xi = [[sum(r[i] * A[i][k] for i in range(D)) % Q if k < N else 0 for k in range(D)] for r in Xi]
    return _rank_rows([[sum(r[i] * A[i][k] for i in range(D)) % Q for k in range(D)] for r in Ok], Q)
  rO = obs_rank(X_good)
  rN = obs_rank(X_bad)
  # Exact projected Block-Wiedemann Hankel recovery: H = [X B^{i+j} Y] (0 <= i < e, 0 <= j <= d)!
  # Compute a_k = X_good * (B^k Y) of shape (m x n) for k = 0..d+e-1
  a_seq = []
  for k in range(d + e):
    BkY = BY[k]  # list of n columns of length D
    ak = [[sum(X_good[r_i][row_idx] * BkY[c_j][row_idx] for row_idx in range(D)) % Q for c_j in range(n)] for r_i in range(m)]
    a_seq.append(ak)
  H_rows = []
  for i in range(e):
    for r_i in range(m):
      row = []
      for j in range(d + 1):
        row.extend(a_seq[i + j][r_i])
      H_rows.append(row)
  rH, kerH = _rref_and_nullspace(H_rows, (d + 1) * n, Q)
  # Map each Hankel kernel vector f = (f_0, ..., f_d) to (w(f), c(f)) and verify exact recovery of v_true!
  recovered_exact = False
  n_c_nonzero = 0
  for fvec in kerH:
    c_rec = fvec[:s_sm]
    if not any(c_rec):
      continue
    n_c_nonzero += 1
    w_rec = [0] * N
    for t in range(n - s_sm):
      ft = fvec[s_sm + t]
      if ft:
        u = U_rand[t]
        for idx in range(N): w_rec[idx] = (w_rec[idx] + ft * u[idx]) % Q
    for j in range(1, d + 1):
      fj = fvec[j * n : (j + 1) * n]
      Bj1Y = BY[j - 1]
      for t in range(n):
        ft = fj[t]
        if ft:
          col = Bj1Y[t]
          for idx in range(N): w_rec[idx] = (w_rec[idx] + ft * col[idx]) % Q
    v_rec = w_rec + c_rec
    # Verify A * v_rec == 0 and v_rec == alpha * v_true with alpha != 0!
    av_zero = all(sum(A[i][k] * v_rec[k] for k in range(D)) % Q == 0 for i in range(D))
    alpha = c_rec[0]
    if av_zero and alpha != 0:
      inv_a = pow(alpha, Q - 2, Q)
      if [(x * inv_a) % Q for x in v_rec] == v_true:
        recovered_exact = True
  ok = (prod_eq_schur and hr_1d_ok and rC == N and rO == rA and rN < rA and recovered_exact)
  return dict(
    checks_ok=ok, ok=ok,
    source="unscreened_live_filtsim_z6_1d_kernel_hankel",
    log2U=log2U, seed=seed, Q=Q, chi=chi,
    filt_N_sw=fres["N_sw"], filt_N_l=fres["N_l"], filt_N_rec=n_rec,
    l_rec_strict_lower_tri=True, prod_eq_schur=prod_eq_schur,
    N=N, D=D, s_sm=s_sm, m=m, n=n, d=d, e=e,
    rank_A=rA, hr_1d_kernel_exact=hr_1d_ok,
    rank_Cd=rC, rank_OeA=rO, control_one_block=rN,
    hankel_shape=[len(H_rows), (d + 1) * n], hankel_rank=rH,
    hankel_c_nonzero_solutions=n_c_nonzero,
    hankel_exact_1d_kernel_recovered=recovered_exact,
  )
def main():
  raw = sys.stdin.read().strip() if not sys.stdin.isatty() else ""
  req = json.loads(raw) if raw else {}
  log2U = int(req.get("log2U", 12))
  seed = int(req.get("seed", 1))
  m = int(req.get("m", 12))
  n = int(req.get("n", 8))
  s_sm = int(req.get("s_sm", 5))
  res = run_hx_1d_kernel_z6(log2U=log2U, seed=seed, m=m, n=n, s_sm=s_sm)
  print(json.dumps(res, separators=(",", ":")))
  if not res["checks_ok"]:
    raise SystemExit(1)
if __name__ == "__main__":
  main()
