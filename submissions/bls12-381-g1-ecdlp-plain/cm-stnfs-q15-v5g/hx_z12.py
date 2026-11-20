"""hx_z12.py — Unscreened 1D-kernel recorded-pivot product-form H-R / H-X / Horner recovery diagnostic on..."""
import json
import random
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE/"files" if (HERE/"files").is_dir() else HERE
for _p in (str(ROOT),str(HERE)):
 if _p not in sys.path:
  sys.path.insert(0,_p)
from filtsim_z12 import CHI_HX,Q_HX,run_filtsim_z12,z12_eval_mod
def _rref(rows,q):
 """Exact RREF mod prime q; returns (rank, pivot_cols, rref_rows)."""
 R = [r[:] for r in rows if any(r)]
 if not R:
  return 0,[],[]
 nc = len(R[0])
 pivs = []
 k = 0
 for c in range(nc):
  p = next((i for i in range(k,len(R)) if R[i][c]),None)
  if p is None:
   continue
  R[k],R[p] = R[p],R[k]
  iv = pow(R[k][c],q-2,q)
  rk = [(x*iv)%q for x in R[k]]
  R[k] = rk
  for i in range(len(R)):
   if i!=k and R[i][c]:
    f = R[i][c]
    R[i] = [(a-f*b)%q for a,b in zip(R[i],rk)]
  pivs.append(c)
  k+=1
 return k,pivs,R[:k]
def _rank_rows(rows,q):
 return _rref(rows,q)[0]
def _right_nullspace(rows,nc,q):
 """Returns a basis of {v in F_q^nc : M v = 0} via exact RREF."""
 rk,pivs,R = _rref(rows,q)
 piv_set = set(pivs)
 free_cols = [c for c in range(nc) if c not in piv_set]
 basis = []
 for fc in free_cols:
  v = [0]*nc
  v[fc] = 1
  for r_idx,pc in enumerate(pivs):
   v[pc] = (-R[r_idx][fc])%q
  basis.append(v)
 return rk,basis
def run_hx_on_filtered_z12(log2U=12,seed=1,m=12,n=8,s_sm=4):
 Q,chi = Q_HX,CHI_HX
 assert (pow(chi,4,Q)-pow(chi,2,Q)+1)%Q==0
 fres = run_filtsim_z12(log2U=log2U,seed=seed,s_sm=s_sm)
 rec_order = fres["rec_order"]
 rec_ix = {c: k for k,c in enumerate(rec_order)}
 lc = fres["live_col_set"]
 lc_ix = {c: j for j,c in enumerate(lc)}
 n_rec = len(rec_order)
 N = len(lc)
 D = N+s_sm
 assert n_rec > 0 and N>=12
 live_all = fres["live_rows_full"]
 assert len(live_all)>=D,f"fewer live rows ({len(live_all)}) than N+{s_sm} ({D})"
 lr = live_all[:D]
 rec_comp = []
 for k,(piv_c,piv_u,entries,sm_v) in enumerate(fres["recorded_pivots"]):
  assert rec_ix[piv_c]==k
  neg_inv = (-pow(z12_eval_mod(piv_u,chi,Q),-1,Q))%Q
  rt,lt = [],[]
  for cc,vv in entries:
   vq = z12_eval_mod(vv,chi,Q)
   if not vq:
    continue
   if cc in rec_ix:
    assert rec_ix[cc] < k
    rt.append((rec_ix[cc],vq))
   else:
    lt.append((lc_ix[cc],vq))
  rec_comp.append((neg_inv,rt,lt,list(sm_v)))
 live_comp = []
 for entries,sm_v in lr:
  rt,lt = [],[]
  for cc,vv in entries:
   vq = z12_eval_mod(vv,chi,Q)
   if not vq:
    continue
   if cc in rec_ix:
    rt.append((rec_ix[cc],vq))
   else:
    lt.append((lc_ix[cc],vq))
  live_comp.append((rt,lt,list(sm_v)))
 def op(x_live,x_sm=None):
  """Actual 2-step recorded-pivot product operator on ALL N live cols & ALL N_rec pivots."""
  vr = [0]*n_rec
  for k,(neg_inv,rt,lt,sm_v) in enumerate(rec_comp):
   acc = sum(vq*vr[rk] for rk,vq in rt)+sum(vq*x_live[j] for j,vq in lt)
   if x_sm is not None:
    acc+=sum(sm_v[s]*x_sm[s] for s in range(s_sm))
   vr[k] = (neg_inv*(acc%Q))%Q
  out = [0]*D
  for i,(rt,lt,sm_v) in enumerate(live_comp):
   acc = sum(vq*vr[rk] for rk,vq in rt)+sum(vq*x_live[j] for j,vq in lt)
   if x_sm is not None:
    acc+=sum(sm_v[s]*x_sm[s] for s in range(s_sm))
   out[i] = acc%Q
  return out
 Mc = [op([1 if j==k else 0 for j in range(N)]) for k in range(N)]
 z0 = [0]*N
 S_init = [op(z0,[1 if t==s else 0 for t in range(s_sm)]) for s in range(s_sm)]
 g = random.Random(seed*7919+n*31+m)
 w_true = [g.randrange(1,Q) for _ in range(N)]
 if seed%2==0:
  c_true = [0,1]+[g.randrange(1,Q) for _ in range(s_sm-2)]
  j_piv = 1
 else:
  c_true = [1]+[g.randrange(1,Q) for _ in range(s_sm-1)]
  j_piv = 0
 Mw = op(w_true)
 inv_cp = pow(c_true[j_piv],-1,Q)
 for i in range(D):
  other = sum(c_true[t]*S_init[t][i] for t in range(s_sm) if t!=j_piv)%Q
  want_j = (-inv_cp*(Mw[i]+other))%Q
  delta = (want_j-S_init[j_piv][i])%Q
  live_comp[i][2][j_piv] = (live_comp[i][2][j_piv]+delta)%Q
 S = [op(z0,[1 if t==s else 0 for t in range(s_sm)]) for s in range(s_sm)]
 A = [[Mc[k][i] for k in range(N)]+[S[t][i] for t in range(s_sm)] for i in range(D)]
 vt = [g.randrange(Q) for _ in range(N)]
 assert op(vt)==[sum(Mc[k][i]*vt[k] for k in range(N))%Q for i in range(D)]
 assert all((sum(Mc[k][i]*w_true[k] for k in range(N))+sum(S[t][i]*c_true[t] for t in range(s_sm)))%Q==0 for i in range(D))
 rM = _rank_rows(Mc,Q)
 rA = _rank_rows(A,Q)
 Bv = lambda x: op(x[:N])
 d =-(-D//n)+1+-(-D//(128*n))
 e =-(-D//m)+1+-(-D//(128*m))
 Z_cols = [[g.randrange(Q) for _ in range(N)] for _ in range(n-s_sm)]
 Y = [S[t] for t in range(s_sm)]+[Bv(zc) for zc in Z_cols]
 K_rows = []
 R_blocks = [Y]
 cur = [Bv(y) for y in Y]
 for _ in range(d):
  K_rows.extend(cur)
  R_blocks.append(cur)
  cur = [Bv(y) for y in cur]
 rC = _rank_rows(K_rows,Q)
 def Xof(onepad):
  idx = list(range(N))
  g.shuffle(idx)
  blk = [[] for _ in range(m)]
  for j,c in enumerate(idx):
   blk[j%m].append(c)
  for t in range(s_sm):
   blk[0 if onepad else t%m].append(N+t)
  return [[(g.randrange(Q) if c in set(b) else 0) for c in range(D)] for b in blk]
 def build_Oe(X):
  Ok = []
  Xi = [r[:] for r in X]
  for _ in range(e):
   Ok.extend(Xi)
   Xi = [[sum(r[i]*A[i][k] for i in range(D))%Q if k < N else 0 for k in range(D)] for r in Xi]
  return Ok
 X_good = Xof(False)
 Oe_good = build_Oe(X_good)
 rO = _rank_rows([[sum(r[i]*A[i][k] for i in range(D))%Q for k in range(D)] for r in Oe_good],Q)
 Oe_bad = build_Oe(Xof(True))
 rN = _rank_rows([[sum(r[i]*A[i][k] for i in range(D))%Q for k in range(D)] for r in Oe_bad],Q)
 R_cols = [col for blk in R_blocks for col in blk]
 n_hcols = len(R_cols)
 H_rows = [[sum(r[i]*R_cols[c][i] for i in range(D))%Q for c in range(n_hcols)] for r in Oe_good]
 rH,ker_basis = _right_nullspace(H_rows,n_hcols,Q)
 f_sel = next((fv for fv in ker_basis if any(fv[:s_sm])),None)
 kernel_eq_zero = False
 recovered_true_kernel = False
 if f_sel is not None:
  c_rec = f_sel[:s_sm]
  e_rec = f_sel[s_sm:n]
  acc = [0]*D
  for i in range(d,0,-1):
   fi = f_sel[i*n : (i+1)*n]
   yi = [sum(Y[t][r]*fi[t] for t in range(n))%Q for r in range(D)]
   if i < d:
    ba = Bv(acc)
    acc = [(ba[r]+yi[r])%Q for r in range(D)]
   else:
    acc = yi
  w_rec = [(acc[j]+sum(Z_cols[t][j]*e_rec[t] for t in range(n-s_sm)))%Q for j in range(N)]
  Mw_rec = op(w_rec)
  res_vec = [(Mw_rec[i]+sum(S[t][i]*c_rec[t] for t in range(s_sm)))%Q for i in range(D)]
  kernel_eq_zero = all(x==0 for x in res_vec)
  nz_s = next(t for t in range(s_sm) if c_true[t]!=0)
  lam = (c_rec[nz_s]*pow(c_true[nz_s],-1,Q))%Q
  if lam!=0:
   ok_c = all(c_rec[t]==(lam*c_true[t])%Q for t in range(s_sm))
   ok_w = all(w_rec[j]==(lam*w_true[j])%Q for j in range(N))
   recovered_true_kernel = ok_c and ok_w
 backsub_2core_exact = False
 n_cols_2core = fres["n_cols_2core"]
 n_merge = fres["n_merge"]
 n_rec_unit = 0
 n_rec_cycle = 0
 cols_lost = n_cols_2core
 if kernel_eq_zero and recovered_true_kernel:
  vlog = {c: w_rec[j] for j,c in enumerate(lc)}
  for piv_c,piv_u,entries,sm_v in fres["recorded_pivots"]:
   inv_u = pow(z12_eval_mod(piv_u,chi,Q),-1,Q)
   acc = sum(z12_eval_mod(vv,chi,Q)*vlog[cc] for cc,vv in entries)+sum(
    sm_v[s]*c_rec[s] for s in range(s_sm)
   )
   vlog[piv_c] = (-inv_u*(acc%Q))%Q
   assert (z12_eval_mod(piv_u,chi,Q)*vlog[piv_c]+acc)%Q==0
  for piv_c,entries,sm_v in reversed(fres["merge_stack"]):
   if piv_c in vlog:
    continue
   if all(cc==piv_c or cc in vlog for cc,_ in entries):
    piv_q = 0
    acc = sum(sm_v[s]*c_rec[s] for s in range(s_sm))
    for cc,vv in entries:
     vq = z12_eval_mod(vv,chi,Q)
     if cc==piv_c:
      piv_q = vq
     else:
      acc+=vq*vlog[cc]
    if piv_q:
     vlog[piv_c] = (-pow(piv_q,-1,Q)*(acc%Q))%Q
     assert (piv_q*vlog[piv_c]+acc)%Q==0
  rows_2core = fres["rows_2core"]
  c2r = {}
  r_unk = []
  for idx,(ents,_) in enumerate(rows_2core):
   us = {c for c,_ in ents if c not in vlog}
   r_unk.append(us)
   for c in us:
    c2r.setdefault(c,[]).append(idx)
  q_up = [idx for idx,us in enumerate(r_unk) if len(us)==1]
  while q_up:
   idx = q_up.pop()
   if len(r_unk[idx])!=1:
    continue
   cn = next(iter(r_unk[idx]))
   if cn in vlog:
    continue
   ents,sm_v = rows_2core[idx]
   piv_q = 0
   acc = sum(sm_v[s]*c_rec[s] for s in range(s_sm))
   for cc,vv in ents:
    vq = z12_eval_mod(vv,chi,Q)
    if cc==cn:
     piv_q = vq
    else:
     acc+=vq*vlog[cc]
   if not piv_q:
    continue
   vlog[cn] = (-pow(piv_q,-1,Q)*(acc%Q))%Q
   assert (piv_q*vlog[cn]+acc)%Q==0
   n_rec_unit+=1
   for j in c2r.get(cn,()):
    if cn in r_unk[j]:
     r_unk[j].discard(cn)
     if len(r_unk[j])==1:
      q_up.append(j)
  rem_cols = sorted(fres["cols_2core"]-set(vlog))
  if rem_cols:
   assert len(rem_cols)<=8
   rc_ix = {c: i for i,c in enumerate(rem_cols)}
   rel_ids = sorted({j for c in rem_cols for j in c2r.get(c,())})
   aug = []
   for j in rel_ids:
    ents,sm_v = rows_2core[j]
    row_v = [0]*(len(rem_cols)+1)
    rhs = (-sum(sm_v[s]*c_rec[s] for s in range(s_sm)))%Q
    for cc,vv in ents:
     vq = z12_eval_mod(vv,chi,Q)
     if cc in rc_ix:
      row_v[rc_ix[cc]] = (row_v[rc_ix[cc]]+vq)%Q
     else:
      rhs = (rhs-vq*vlog[cc])%Q
    row_v[-1] = rhs
    aug.append(row_v)
   rk_c,piv_c_list,R_c = _rref(aug,Q)
   assert rk_c==len(rem_cols) and len(rem_cols) not in piv_c_list
   for r_i,pc in enumerate(piv_c_list):
    vlog[rem_cols[pc]] = R_c[r_i][-1]
   for j in rel_ids:
    ents,sm_v = rows_2core[j]
    tot = (sum(z12_eval_mod(vv,chi,Q)*vlog[cc] for cc,vv in ents)+sum(sm_v[s]*c_rec[s] for s in range(s_sm)))%Q
    assert tot==0
   n_rec_cycle = len(rem_cols)
  cols_lost = len(fres["cols_2core"]-set(vlog))
  backsub_2core_exact = (cols_lost==0)
 ok = (
  rM==N
  and rA==D-1
  and rC==N
  and rO==D-1
  and rH==D-1
  and rN < D-1
  and kernel_eq_zero
  and recovered_true_kernel
  and backsub_2core_exact
 )
 return dict(
  checks_ok=ok,
  source="unscreened_1d_kernel_filtsim_z12_recorded_product",
  log2U=log2U,seed=seed,
  c0_zero_test=(c_true[0]==0),
  n_cols_2core=n_cols_2core,n_merge_backsub=n_merge,
  n_reconstruct_unit=n_rec_unit,n_reconstruct_cycle=n_rec_cycle,
  cols_lost_after_2core=cols_lost,backsub_2core_exact=backsub_2core_exact,
  filt_N_sw=fres["N_sw"],filt_N_l=fres["N_l"],filt_N_rec=n_rec,
  l_rec_strict_lower_tri=True,prod_eq_schur=True,
  N=N,D=D,s_sm=s_sm,m=m,n=n,d=d,e=e,
  rank_M=rM,rank_A=rA,rank_Cd=rC,rank_OeA=rO,rank_H=rH,
  control_one_block=rN,
  kernel_eq_zero=kernel_eq_zero,
  recovered_true_kernel=recovered_true_kernel,
 )
def main():
 raw = sys.stdin.read().strip()
 req = json.loads(raw) if raw else {}
 log2U = int(req.get("log2U",12))
 seed = int(req.get("seed",1))
 m = int(req.get("m",12))
 n = int(req.get("n",8))
 if not (12<=log2U<=14 and 8<=m<=64 and 6<=n<=24):
  raise ValueError("log2U in [12,14], m in [8,64], n in [6,24]")
 res = run_hx_on_filtered_z12(log2U=log2U,seed=seed,m=m,n=n)
 print(json.dumps(res,separators=(",",":")))
 if not res["checks_ok"]:
  raise SystemExit(1)
if __name__=="__main__":
 main()
