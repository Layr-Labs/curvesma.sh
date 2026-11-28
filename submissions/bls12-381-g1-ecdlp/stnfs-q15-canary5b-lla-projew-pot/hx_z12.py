"""hx_z12.py — Unscreened 1D-kernel recorded-pivot product-form H-R / H-X / structured PM-Basis / Horner recovery diagnostic on Z[zeta_12]."""
import json
import random
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
for _p in (str(ROOT),str(HERE)):
 if _p not in sys.path:sys.path.insert(0,_p)
from filtsim_z12 import CHI_HX,Q_HX,run_filtsim_z12,z12_eval_mod
def _rref(rows,q):
 R=[r[:]for r in rows if any(r)]
 if not R:return 0,[],[]
 nc=len(R[0])
 pivs=[]
 k=0
 for c in range(nc):
  p=next((i for i in range(k,len(R))if R[i][c]),None)
  if p is None:continue
  R[k],R[p]=R[p],R[k]
  iv=pow(R[k][c],q-2,q)
  rk=[(x*iv)%q for x in R[k]]
  R[k]=rk
  for i in range(len(R)):
   if i!=k and R[i][c]:
    f=R[i][c]
    R[i]=[(a-f*b)%q for a,b in zip(R[i],rk)]
  pivs.append(c)
  k+=1
 return k,pivs,R[:k]
def _rank_rows(rows,q):
 return _rref(rows,q)[0]
def _right_nullspace(rows,nc,q):
 rk,pivs,R=_rref(rows,q)
 piv_set=set(pivs)
 free_cols=[c for c in range(nc)if c not in piv_set]
 basis=[]
 for fc in free_cols:
  v=[0]*nc
  v[fc]=1
  for r_idx,pc in enumerate(pivs):
   v[pc]=(-R[r_idx][fc])%q
  basis.append(v)
 return rk,basis
def run_hx_on_filtered_z12(log2U=12,seed=1,m=12,n=8,s_sm=4,lf=71.50,lg=71.86,nf=543.4,ng=664.4,excess=0.84):
 Q,chi=Q_HX,CHI_HX
 assert(pow(chi,4,Q)-pow(chi,2,Q)+1)%Q==0
 g=random.Random(seed*7919+n*31+m)
 if seed%2==0:
  c_true=[0,1]+[g.randrange(1,Q)for _ in range(s_sm-2)]
  j_piv=1
 else:
  c_true=[1]+[g.randrange(1,Q)for _ in range(s_sm-1)]
  j_piv=0
 fres=run_filtsim_z12(log2U=log2U,seed=seed,lf=lf,lg=lg,nf=nf,ng=ng,excess=excess,s_sm=s_sm,planted_kernel=(c_true,j_piv,seed*7919+n*31+m+101))
 w_all=fres["w_all"]
 rec_order=fres["rec_order"]
 rec_ix={c:k for k,c in enumerate(rec_order)}
 lc=fres["live_col_set"]
 lc_ix={c:j for j,c in enumerate(lc)}
 n_rec=len(rec_order)
 N=len(lc)
 D=N+s_sm
 assert n_rec>0 and N>=12
 live_all=fres["live_rows_full"]
 assert len(live_all)>=D,f"fewer live rows ({len(live_all)}) than N+{s_sm} ({D})"
 _ldeg={}
 for _ents,_ in live_all:
  for _cc,_ in _ents:
   if _cc in lc_ix:_ldeg[_cc]=_ldeg.get(_cc,0)+1
 _cov0=set(_ldeg)
 _keep=[True]*len(live_all);_rem=len(live_all)
 for _idx in sorted(range(len(live_all)),key=lambda i:(-len(live_all[i][0]),-i)):
  if _rem<=D:break
  _cs=[_cc for _cc,_ in live_all[_idx][0]if _cc in lc_ix]
  if all(_ldeg[_cc]>=2 for _cc in _cs):
   _keep[_idx]=False;_rem-=1
   for _cc in _cs:_ldeg[_cc]-=1
 assert _rem==D and {c for c,dv in _ldeg.items()if dv>=1}==_cov0,"Degree-preserving D-row trim lost a live column!"
 lr=[live_all[i]for i in range(len(live_all))if _keep[i]]
 rec_comp=[]
 for k,(piv_c,piv_u,entries,sm_v) in enumerate(fres["recorded_pivots"]):
  assert rec_ix[piv_c]==k
  neg_inv=(-pow(z12_eval_mod(piv_u,chi,Q),-1,Q))%Q
  rt,lt=[],[]
  for cc,vv in entries:
   vq=z12_eval_mod(vv,chi,Q)
   if not vq:continue
   if cc in rec_ix:
    assert rec_ix[cc]<k
    rt.append((rec_ix[cc],vq))
   else:lt.append((lc_ix[cc],vq))
  rec_comp.append((neg_inv,rt,lt,list(sm_v)))
 live_comp=[]
 for entries,sm_v in lr:
  rt,lt=[],[]
  for cc,vv in entries:
   vq=z12_eval_mod(vv,chi,Q)
   if not vq:continue
   if cc in rec_ix:rt.append((rec_ix[cc],vq))
   else:lt.append((lc_ix[cc],vq))
  live_comp.append((rt,lt,list(sm_v)))
 def op(x_live,x_sm=None):
  vr=[0]*n_rec
  for k,(neg_inv,rt,lt,sm_v) in enumerate(rec_comp):
   acc=sum(vq*vr[rk]for rk,vq in rt)+sum(vq*x_live[j]for j,vq in lt)
   if x_sm is not None:acc+=sum(sm_v[s]*x_sm[s]for s in range(s_sm))
   vr[k]=(neg_inv*(acc%Q))%Q
  out=[0]*D
  for i,(rt,lt,sm_v) in enumerate(live_comp):
   acc=sum(vq*vr[rk]for rk,vq in rt)+sum(vq*x_live[j]for j,vq in lt)
   if x_sm is not None:acc+=sum(sm_v[s]*x_sm[s]for s in range(s_sm))
   out[i]=acc%Q
  return out
 Mc=[op([1 if j==k else 0 for j in range(N)])for k in range(N)]
 z0=[0]*N
 w_true=[w_all[c]for c in lc]
 S=[op(z0,[1 if t==s else 0 for t in range(s_sm)])for s in range(s_sm)]
 A=[[Mc[k][i]for k in range(N)]+[S[t][i]for t in range(s_sm)]for i in range(D)]
 # Independent block-Gaussian Schur complement [R_ll - R_lr P_rr^{-1} P_rl | R_ls - R_lr P_rr^{-1} P_rs] directly from raw recorded_pivots & lr (never calling op()):
 P_rr=[[0]*n_rec for _ in range(n_rec)]
 P_rl=[[0]*N for _ in range(n_rec)]
 P_rs=[[0]*s_sm for _ in range(n_rec)]
 for k,(piv_c,piv_u,entries,sm_v) in enumerate(fres["recorded_pivots"]):
  P_rr[k][k]=z12_eval_mod(piv_u,chi,Q)
  for cc,vv in entries:
   vq=z12_eval_mod(vv,chi,Q)
   if cc in rec_ix:P_rr[k][rec_ix[cc]]=(P_rr[k][rec_ix[cc]]+vq)%Q
   else:P_rl[k][lc_ix[cc]]=(P_rl[k][lc_ix[cc]]+vq)%Q
  for s in range(s_sm):P_rs[k][s]=sm_v[s]%Q
 R_lr=[[0]*n_rec for _ in range(D)]
 R_ll=[[0]*N for _ in range(D)]
 R_ls=[[0]*s_sm for _ in range(D)]
 for i,(entries,sm_v) in enumerate(lr):
  for cc,vv in entries:
   vq=z12_eval_mod(vv,chi,Q)
   if cc in rec_ix:R_lr[i][rec_ix[cc]]=(R_lr[i][rec_ix[cc]]+vq)%Q
   else:R_ll[i][lc_ix[cc]]=(R_ll[i][lc_ix[cc]]+vq)%Q
  for s in range(s_sm):R_ls[i][s]=sm_v[s]%Q
 X_rl=[[0]*N for _ in range(n_rec)]
 X_rs=[[0]*s_sm for _ in range(n_rec)]
 for k in range(n_rec):
  inv_d=pow(P_rr[k][k],-1,Q)
  for j in range(N):
   X_rl[k][j]=((P_rl[k][j]-sum(P_rr[k][t]*X_rl[t][j]for t in range(k)))*inv_d)%Q
  for s in range(s_sm):
   X_rs[k][s]=((P_rs[k][s]-sum(P_rr[k][t]*X_rs[t][s]for t in range(k)))*inv_d)%Q
 M_schur=[[(R_ll[i][j]-sum(R_lr[i][k]*X_rl[k][j]for k in range(n_rec)))%Q for j in range(N)]for i in range(D)]
 S_schur=[[(R_ls[i][s]-sum(R_lr[i][k]*X_rs[k][s]for k in range(n_rec)))%Q for s in range(s_sm)]for i in range(D)]
 assert all(M_schur[i][k]==Mc[k][i]for i in range(D)for k in range(N)),"Product-form op M mismatch vs block-Gaussian Schur complement!"
 assert all(S_schur[i][t]==S[t][i]for i in range(D)for t in range(s_sm)),"Product-form op S mismatch vs block-Gaussian Schur complement!"
 vt=[g.randrange(Q)for _ in range(N)]
 assert op(vt)==[sum(M_schur[i][k]*vt[k]for k in range(N))%Q for i in range(D)]
 assert all((sum(Mc[k][i]*w_true[k]for k in range(N))+sum(S[t][i]*c_true[t]for t in range(s_sm)))%Q==0 for i in range(D))
 rM=_rank_rows(Mc,Q)
 rA=_rank_rows(A,Q)
 Bv=lambda x:op(x[:N])
 d=-(-D//n)+1+-(-D//(128*n))
 e=-(-D//m)+1+-(-D//(128*m))
 Z_cols=[[g.randrange(Q)for _ in range(N)]for _ in range(n-s_sm)]
 Y=[S[t]for t in range(s_sm)]+[Bv(zc)for zc in Z_cols]
 K_rows=[]
 R_blocks=[Y]
 cur=[Bv(y)for y in Y]
 for _ in range(d):
  K_rows.extend(cur)
  R_blocks.append(cur)
  cur=[Bv(y)for y in cur]
 rC=_rank_rows(K_rows,Q)
 def Xof(onepad):
  idx=list(range(N))
  g.shuffle(idx)
  blk=[[]for _ in range(m)]
  for j,c in enumerate(idx):
   blk[j%m].append(c)
  for t in range(s_sm):
   blk[0 if onepad else t%m].append(N+t)
  return[[(g.randrange(Q)if c in set(b)else 0)for c in range(D)]for b in blk]
 def build_Oe(X):
  Ok=[]
  Xi=[r[:]for r in X]
  for _ in range(e):
   Ok.extend(Xi)
   Xi=[[sum(r[i]*A[i][k]for i in range(D))%Q if k<N else 0 for k in range(D)]for r in Xi]
  return Ok
 X_good=Xof(False)
 Oe_good=build_Oe(X_good)
 rO=_rank_rows([[sum(r[i]*A[i][k]for i in range(D))%Q for k in range(D)]for r in Oe_good],Q)
 Oe_bad=build_Oe(Xof(True))
 rN=_rank_rows([[sum(r[i]*A[i][k]for i in range(D))%Q for k in range(D)]for r in Oe_bad],Q)
 R_cols=[col for blk in R_blocks for col in blk]
 n_hcols=len(R_cols)
 H_rows=[[sum(r[i]*R_cols[c][i]for i in range(D))%Q for c in range(n_hcols)]for r in Oe_good]
 rH,ker_basis=_right_nullspace(H_rows,n_hcols,Q)
 f_sel=next((fv for fv in ker_basis if any(fv[:s_sm])),None)
 kernel_eq_zero=False
 recovered_true_kernel=False
 if f_sel is not None:
  c_rec=f_sel[:s_sm]
  e_rec=f_sel[s_sm:n]
  acc=[0]*D
  for i in range(d,0,-1):
   fi=f_sel[i*n:(i+1)*n]
   yi=[sum(Y[t][r]*fi[t]for t in range(n))%Q for r in range(D)]
   if i<d:
    ba=Bv(acc)
    acc=[(ba[r]+yi[r])%Q for r in range(D)]
   else:acc=yi
  w_rec=[(acc[j]+sum(Z_cols[t][j]*e_rec[t]for t in range(n-s_sm)))%Q for j in range(N)]
  Mw_rec=op(w_rec)
  res_vec=[(Mw_rec[i]+sum(S[t][i]*c_rec[t]for t in range(s_sm)))%Q for i in range(D)]
  kernel_eq_zero=all(x==0 for x in res_vec)
  nz_s=next(t for t in range(s_sm)if c_true[t]!=0)
  lam=(c_rec[nz_s]*pow(c_true[nz_s],-1,Q))%Q
  if lam!=0:
   ok_c=all(c_rec[t]==(lam*c_true[t])%Q for t in range(s_sm))
   ok_w=all(w_rec[j]==(lam*w_true[j])%Q for j in range(N))
   recovered_true_kernel=ok_c and ok_w
 # Run production-aligned recursive PM-Basis directly on the recorded-pivot product-form [M|0] operator Bv, k>=0, Y=[S, B Z], n=8 > s_sm=4!
 L_pm=e+d
 seq0=[]
 cur_b=Y
 for _ in range(L_pm):
  seq0.append([[sum(X_good[ri][k]*cur_b[bi][k]for k in range(D))%Q for bi in range(n)]for ri in range(m)])
  cur_b=[Bv(col)for col in cur_b]
 Y_rows=[[Y[bi][ri]for bi in range(n)]for ri in range(D)]
 Z_rows=[[Z_cols[t][j]for t in range(n-s_sm)]for j in range(N)]
 def _Av_op(v):
  mw=op(v[:N])
  return[(mw[i]+sum(S[t][i]*v[N+t]for t in range(s_sm)))%Q for i in range(D)]
 pm_prod=pm_basis_verify(Bv,_Av_op,seq0,Y_rows,Z_rows,Q,e,d+1,s_sm,w_true+c_true)
 backsub_2core_exact=False
 n_cols_2core=fres["n_cols_2core"]
 n_merge=fres["n_merge"]
 n_rec_unit=0
 n_repurge_backsub=0
 n_rec_cycle=0
 cols_lost=n_cols_2core
 unrec_frac=1.0
 occ_unrec_frac=1.0
 n_2core_rels_checked=0
 all_2core_rels_zero=False
 recovered_match_true=False
 if kernel_eq_zero and recovered_true_kernel:
  vlog={c:w_rec[j]for j,c in enumerate(lc)}
  for piv_c,piv_u,entries,sm_v in fres["recorded_pivots"]:
   inv_u=pow(z12_eval_mod(piv_u,chi,Q),-1,Q)
   acc=sum(z12_eval_mod(vv,chi,Q)*vlog[cc]for cc,vv in entries)+sum(sm_v[s]*c_rec[s]for s in range(s_sm))
   vlog[piv_c]=(-inv_u*(acc%Q))%Q
   assert(z12_eval_mod(piv_u,chi,Q)*vlog[piv_c]+acc)%Q==0
  for piv_c,entries,sm_v in reversed(fres["merge_stack"]):
   if piv_c in vlog:continue
   if all(cc==piv_c or cc in vlog for cc,_ in entries):
    piv_q=0
    acc=sum(sm_v[s]*c_rec[s]for s in range(s_sm))
    for cc,vv in entries:
     vq=z12_eval_mod(vv,chi,Q)
     if cc==piv_c:piv_q=vq
     else:acc+=vq*vlog[cc]
    if piv_q:
     vlog[piv_c]=(-pow(piv_q,-1,Q)*(acc%Q))%Q
     assert(piv_q*vlog[piv_c]+acc)%Q==0
  n_repurge_backsub=0
  for piv_c,entries,sm_v in reversed(fres["repurge_stack"]):
   if piv_c in vlog:continue
   if all(cc==piv_c or cc in vlog for cc,_ in entries):
    piv_q=0
    acc=sum(sm_v[s]*c_rec[s]for s in range(s_sm))
    for cc,vv in entries:
     vq=z12_eval_mod(vv,chi,Q)
     if cc==piv_c:piv_q=vq
     else:acc+=vq*vlog[cc]
    if piv_q:
     vlog[piv_c]=(-pow(piv_q,-1,Q)*(acc%Q))%Q
     assert(piv_q*vlog[piv_c]+acc)%Q==0
     n_repurge_backsub+=1
  rows_2core=fres["rows_2core"]
  while True:
   _prog_u=False
   for ents,sm_v in rows_2core:
    unk=[cc for cc,_ in ents if cc not in vlog]
    if len(unk)==1:
     uc=unk[0];pq=0;acc=sum(sm_v[s]*c_rec[s]for s in range(s_sm))
     for cc,vv in ents:
      vq=z12_eval_mod(vv,chi,Q)
      if cc==uc:pq=vq
      else:acc+=vq*vlog[cc]
     if pq:
      vlog[uc]=(-pow(pq,-1,Q)*(acc%Q))%Q;n_rec_unit+=1;_prog_u=True
   if not _prog_u:break
  rem_cols=sorted(fres["cols_2core"]-set(vlog))
  rem_set=set(rem_cols)
  tot_occ=sum(len(ents)for ents,_ in rows_2core)
  unrec_occ=sum(1 for ents,_ in rows_2core for cc,_ in ents if cc in rem_set)
  occ_unrec_frac=unrec_occ/max(1,tot_occ)
  n_2core_rels_checked=0
  all_2core_rels_zero=True
  for ents,sm_v in rows_2core:
   if all(cc in vlog for cc,_ in ents):
    tot=(sum(z12_eval_mod(vv,chi,Q)*vlog[cc]for cc,vv in ents)+sum(sm_v[s]*c_rec[s]for s in range(s_sm)))%Q
    if tot!=0:
     all_2core_rels_zero=False
     break
    n_2core_rels_checked+=1
  recovered_match_true=all(vlog[c]==(lam*w_all[c])%Q for c in vlog)
  cols_lost=len(rem_cols)
  unrec_frac=cols_lost/max(1,n_cols_2core)
  n_purge_backsub=0;purge_backsub_errors=0
  for piv_c,entries,sm_v in reversed(fres["purge_stack"]):
   if piv_c in vlog:continue
   if all(cc==piv_c or cc in vlog for cc,_ in entries):
    piv_q=0;acc=sum(sm_v[s]*c_rec[s]for s in range(s_sm))
    for cc,vv in entries:
     vq=z12_eval_mod(vv,chi,Q)
     if cc==piv_c:piv_q=vq
     else:acc+=vq*vlog[cc]
    if piv_q:
     vc=(-pow(piv_q,-1,Q)*(acc%Q))%Q
     vlog[piv_c]=vc;n_purge_backsub+=1
     if vc!=(lam*w_all[piv_c])%Q or (piv_q*vc+acc)%Q!=0:purge_backsub_errors+=1
  backsub_2core_exact=(all_2core_rels_zero and recovered_match_true and fres["post_peel_cols_lost"]==0 and fres["post_merge_cols_lost"]==0 and fres["repurge_backsub_errors"]==0 and n_repurge_backsub==fres["repurge_stack_len"]and n_rec_unit==0 and cols_lost==0 and n_2core_rels_checked==len(rows_2core)and purge_backsub_errors==0 and n_purge_backsub>0)
 ok=(rM==N and rA==D-1 and rC==N and rO==D-1 and rH==D-1 and rN<D-1 and kernel_eq_zero and recovered_true_kernel and backsub_2core_exact and pm_prod["pm_basis_matches_true_ker"])
 return dict(checks_ok=ok,source="unscreened_1d_kernel_filtsim_z12_recorded_product",log2U=log2U,seed=seed,lf=lf,lg=lg,excess=excess,final_surplus=len(live_all)-N,c0_zero_test=(c_true[0]==0),n_cols_2core=n_cols_2core,n_merge_backsub=n_merge,repurge_stack_len=fres["repurge_stack_len"],n_repurge_backsub=n_repurge_backsub,purge_stack_len=fres["purge_stack_len"],n_purge_backsub=n_purge_backsub,purge_backsub_errors=purge_backsub_errors,post_peel_cols_lost=fres["post_peel_cols_lost"],post_merge_cols_lost=fres["post_merge_cols_lost"],repurge_backsub_errors=fres["repurge_backsub_errors"],n_reconstruct_unit=n_rec_unit,cols_unrecovered_tri=cols_lost,unrecovered_frac=round(unrec_frac,6),occ_unrecovered_frac=round(occ_unrec_frac,6),n_2core_rels_verified_zero=n_2core_rels_checked,all_2core_rels_zero=all_2core_rels_zero,recovered_cols_match_true=recovered_match_true,backsub_2core_exact=backsub_2core_exact,filt_N_sw=fres["N_sw"],filt_N_l=fres["N_l"],filt_N_rec=n_rec,l_rec_strict_lower_tri=True,prod_eq_schur=True,N=N,D=D,s_sm=s_sm,m=m,n=n,d=d,e=e,rank_M=rM,rank_A=rA,rank_Cd=rC,rank_OeA=rO,rank_H=rH,control_one_block=rN,kernel_eq_zero=kernel_eq_zero,recovered_true_kernel=recovered_true_kernel,**pm_prod)
def main():
 raw=sys.stdin.read().strip()
 req=json.loads(raw)if raw else{}
 log2U=int(req.get("log2U",13))
 seed=int(req.get("seed",2))
 m=int(req.get("m",12))
 n=int(req.get("n",8))
 if not(12<=log2U<=14 and 8<=m<=64 and 6<=n<=24):raise ValueError("log2U in [12,14], m in [8,64], n in [6,24]")
 res=run_hx_on_filtered_z12(log2U=log2U,seed=seed,m=m,n=n)
 print(json.dumps(res,separators=(",",":")))
 if not res["checks_ok"]:raise SystemExit(1)

def _padd(a,b,p):
 n=max(len(a),len(b));o=[0]*n
 for i,v in enumerate(a):o[i]=v
 for i,v in enumerate(b):o[i]=(o[i]+v)%p
 while o and o[-1]==0:o.pop()
 return o
def _pmul(a,b,p,cap=None):
 if not a or not b:return[]
 la,lb=len(a),len(b);n=la+lb-1 if cap is None else min(la+lb-1,cap);o=[0]*n
 for i,ai in enumerate(a):
  if ai:
   for j in range(min(lb,n-i)):o[i+j]=(o[i+j]+ai*b[j])%p
 while o and o[-1]==0:o.pop()
 return o
def _pmat_mul(A,B,p,cap=None):
 r,k,c=len(A),len(A[0]),len(B[0])
 C=[[[]for _ in range(c)]for _ in range(r)]
 for i in range(r):
  for t in range(k):
   ait=A[i][t]
   if not ait:continue
   bt=B[t]
   for j in range(c):
    if bt[j]:C[i][j]=_padd(C[i][j],_pmul(ait,bt[j],p,cap),p)
 return C
def _mbasis1(F0,s,p):
 m,t=len(F0),len(F0[0]);order=sorted(range(t),key=lambda j:(s[j],j))
 pivs={};Pcols=[[0]*t for _ in range(t)];dinc=[0]*t
 for cj in order:
  vec=[F0[r][cj]%p for r in range(m)];pc=[0]*t;pc[cj]=1
  for r in range(m):
   v=vec[r]
   if not v:continue
   if r in pivs:
    pv,pp=pivs[r]
    for rr in range(r,m):
     if pv[rr]:vec[rr]=(vec[rr]-v*pv[rr])%p
    for cc in range(t):
     if pp[cc]:pc[cc]=(pc[cc]-v*pp[cc])%p
   else:
    inv=pow(v,-1,p);pivs[r]=([(x*inv)%p for x in vec],[(x*inv)%p for x in pc]);dinc[cj]=1;break
  Pcols[cj]=pc
 P=[[[]for _ in range(t)]for _ in range(t)]
 for j in range(t):
  for i in range(t):
   v=Pcols[j][i]
   if v:P[i][j]=[0,v]if dinc[j]else[v]
 return P,[s[j]+dinc[j]for j in range(t)]
def pm_basis(F,sigma,s,p):
 """Giorgi-Jeannerod-Villard (2003) / Thome (2002) divide-and-conquer PM-Basis over F_p[z]."""
 m,t=len(F),len(F[0])
 if sigma==1:return _mbasis1([[F[i][j][0]if F[i][j]else 0 for j in range(t)]for i in range(m)],s,p)
 s1=sigma//2;s2=sigma-s1
 P1,sm=pm_basis([[F[i][j][:s1]for j in range(t)]for i in range(m)],s1,s,p)
 FP1=_pmat_mul(F,P1,p,cap=sigma)
 G=[[FP1[i][j][s1:sigma]for j in range(t)]for i in range(m)]
 P2,so=pm_basis(G,s2,sm,p)
 return _pmat_mul(P1,P2,p),so
def pm_basis_verify(Bv_fn,Av_fn,seq0,Y_rows,Z_rows,r,hr,hc,s_sm,v_ref):
 """Production-aligned divide-and-conquer PM-Basis on S_k = X B^k Y (k >= 0) with structured Y = [S, B Z], n > s_sm, exact leading-coeff [z^{s_out[j]}] q_j(z)[:s_sm] != 0 selector, zero-SM-lead negative control, and 1-pass Horner recovery."""
 m=len(seq0[0]);n=len(seq0[0][0]);D=len(Y_rows);N=D-s_sm;L=hr+hc-1
 assert n>s_sm and len(Z_rows)==N and len(Z_rows[0])==n-s_sm and len(seq0)>=L
 S=seq0[:L]
 M=[[[]for _ in range(n+m)]for _ in range(m)]
 for i in range(m):
  for j in range(n):
   poly=[S[k][i][j]for k in range(L)]
   while poly and poly[-1]==0:poly.pop()
   M[i][j]=poly
  M[i][n+i]=[r-1]
 P_b,s_out=pm_basis(M,L,[0]*n+[1]*m,r)
 def _eval_horner(q_poly,h_deg):
  lead_v=[q_poly[bi][h_deg]if h_deg<len(q_poly[bi])else 0 for bi in range(n)]
  c_v=lead_v[:s_sm];e_v=lead_v[s_sm:n]
  acc=[0]*D
  for t in range(h_deg):
   qt=[q_poly[bi][t]if t<len(q_poly[bi])else 0 for bi in range(n)]
   yc=[sum(Y_rows[ri][bi]*qt[bi]for bi in range(n))%r for ri in range(D)]
   acc=[(a+y)%r for a,y in zip(Bv_fn(acc),yc)]
  w_v=[(acc[j]+sum(Z_rows[j][t]*e_v[t]for t in range(n-s_sm)))%r for j in range(N)]
  return w_v+c_v,lead_v
 v_pm=None;deg_used=None;horner_passes=0;cands_valid=[];ctrl_zero_ok=False
 for col in sorted(range(n+m),key=lambda j:s_out[j]):
  h=s_out[col]
  if h>hc-1:continue
  q=[P_b[i][col]for i in range(n)]
  lead=[q[bi][h]if h<len(q[bi])else 0 for bi in range(n)]
  if not any(lead):continue
  cands_valid.append((h,q,lead))
  if not any(lead[:s_sm]):
   # Naturally occurring zero-SM-lead column (any(lead) True, any(lead[:s_sm]) False): verify Horner yields [0]*D!
   z_cand,_=_eval_horner(q,h)
   assert all(x==0 for x in z_cand),"Zero-SM-lead column unexpectedly produced nonzero vector!"
   ctrl_zero_ok=True
   continue
  if v_pm is None:
   dg=max((len(poly)-1 for poly in q),default=-1)
   assert dg==h and any(lead[:s_sm]),"Leading [z^{s_out[j]}] q_j(z)[:s_sm] selector mismatch!"
   horner_passes+=1
   cand,_=_eval_horner(q,h)
   if any(cand)and all(x==0 for x in Av_fn(cand)):
    v_pm=cand;deg_used=dg
 # If no column had lead[:s_sm]==0 naturally, construct an exact degree-aligned PM-Basis module combination q_ctrl with lead[:s_sm]==0 and lead[s_sm:]!=0 (guaranteed since n > s_sm) and verify Horner gives [0]*D!
 if not ctrl_zero_ok:
  d_max=max(h for h,_,_ in cands_valid)
  sm_mat=[[lead[s]for _,_,lead in cands_valid]for s in range(s_sm)]
  _,ker_sm=_right_nullspace(sm_mat,len(cands_valid),r)
  assert ker_sm,"Must exist nonzero combination annihilating lead[:s_sm]!"
  alpha=ker_sm[0]
  q_ctrl=[[]for _ in range(n)]
  for coef,(h_u,q_u,_) in zip(alpha,cands_valid):
   if not coef:continue
   sh=d_max-h_u
   for bi in range(n):
    if q_u[bi]:
     shifted=[0]*sh+[(coef*x)%r for x in q_u[bi]]
     q_ctrl[bi]=_padd(q_ctrl[bi],shifted,r)
  z_cand,z_lead=_eval_horner(q_ctrl,d_max)
  assert(not any(z_lead[:s_sm]))and any(z_lead[s_sm:n])and all(x==0 for x in z_cand),"Negative control: zero-SM-lead nonzero-tail PM-Basis vector must map to [0]*D!"
  ctrl_zero_ok=True
 assert v_pm is not None and horner_passes==1 and ctrl_zero_ok,f"PM-Basis failed 1-pass [:s_sm] recovery (passes={horner_passes})!"
 idx0=next(i for i,x in enumerate(v_ref)if x)
 ratio=(v_pm[idx0]*pow(v_ref[idx0],-1,r))%r
 assert ratio!=0 and all((v_ref[i]*ratio)%r==v_pm[i]for i in range(D))
 return dict(pm_basis_order=L,pm_basis_cols=n+m,pm_basis_n=n,pm_basis_s_sm=s_sm,pm_basis_deg=deg_used,pm_basis_lead_nonzero=True,pm_basis_lead_sm_nonzero=True,pm_basis_zero_sm_ctrl_zero_vec=ctrl_zero_ok,pm_basis_horner_passes=horner_passes,pm_basis_matches_true_ker=True)

def unplanted_two_sided_bw_check():
 import math, random
 p,r,m=467,233,381
 f_c,g_c=[1,1,-5,0,1],[-10,-1,1]
 assert sum(ci*pow(m,i,p)for i,ci in enumerate(f_c))%p==0
 assert sum(ci*pow(m,i,p)for i,ci in enumerate(g_c))%p==0
 eps_f,eps_g,r2=r-1,r*r-1,r*r
 def peval(c,x,mod):
  s=0
  for ci in reversed(c):
   s=(s*x+ci)%mod
  return s
 def pmul(a,b,mod_p,q):
  d=len(mod_p)-1
  out=[0]*(len(a)+len(b)-1)
  for i,ai in enumerate(a):
   if ai:
    for j,bj in enumerate(b):
     out[i+j]=(out[i+j]+ai*bj)%q
  for k in range(len(out)-1,d-1,-1):
   ck=out[k]
   if ck:
    for j in range(d):
     out[k-d+j]=(out[k-d+j]-ck*mod_p[j])%q
  return out[:d]
 def ppow(base,exp,mod_p,q):
  d=len(mod_p)-1
  res=[1]+[0]*(d-1)
  b=[x%q for x in base]+[0]*(d-len(base))
  e=exp
  while e:
   if e&1:res=pmul(res,b,mod_p,q)
   e>>=1
   if e:b=pmul(b,b,mod_p,q)
  return res
 primes=[q for q in range(2,115)if all(q%d for d in range(2,int(q**0.5)+1))]
 def unram(c,q):
  for rt in range(q):
   if peval(c,rt,q)==0:
    if sum(i*c[i]*pow(rt,i-1,q)for i in range(1,len(c)))%q==0:return False
  return True
 fb_f,fb_g={},{}
 for q in primes:
  if unram(f_c,q):
   for rt in range(q):
    if peval(f_c,rt,q)==0:fb_f[(q,rt)]=len(fb_f)
  if unram(g_c,q):
   for rt in range(q):
    if peval(g_c,rt,q)==0:fb_g[(q,rt)]=len(fb_g)
 pf_act=sorted({q for(q,_) in fb_f})
 pg_act=sorted({q for(q,_) in fb_g})
 def sfac(val,plist):
  n=abs(val)
  if n==0:return None
  fac={}
  for q in plist:
   if n%q==0:
    e=0
    while n%q==0:
     n//=q
     e+=1
    fac[q]=e
    if n==1:return fac
  return fac if n==1 else None
 rels=[]
 for b in range(1,55):
  for a in range(-450,451):
   if a==0 or math.gcd(abs(a),b)!=1 or(a-b*m)%p==0:continue
   fg=sfac(a*a-a*b-10*b*b,pg_act)
   if fg is None:continue
   ff=sfac(a**4-5*a*a*b*b+a*b**3+b**4,pf_act)
   if ff is None:continue
   vf={fb_f[(q,(a*pow(b,-1,q))%q)]:e for q,e in ff.items()}
   vg={fb_g[(q,(a*pow(b,-1,q))%q)]:e for q,e in fg.items()}
   rels.append((a,b,vf,vg))
 nf_cnt=len(fb_f)
 cur=list(range(len(rels)))
 while True:
  ccnt={}
  for ri in cur:
   for k in rels[ri][2]:
    ccnt[k]=ccnt.get(k,0)+1
   for k in rels[ri][3]:
    ccnt[nf_cnt+k]=ccnt.get(nf_cnt+k,0)+1
  nxt=[ri for ri in cur if all(ccnt[k]>=2 for k in rels[ri][2])and all(ccnt[nf_cnt+k]>=2 for k in rels[ri][3])]
  if len(nxt)==len(cur):break
  cur=nxt
 used=sorted({k for ri in cur for k in rels[ri][2]}|{nf_cnt+k for ri in cur for k in rels[ri][3]})
 cmap={c:i for i,c in enumerate(used)}
 N_2core=len(used)
 def smap(a,b,mod_p,eps):
  el=ppow([a%r2,(-b)%r2],eps,mod_p,r2)
  assert el[0]%r==1 and all(x%r==0 for x in el[1:])
  return[((el[0]-1)//r)%r]+[(x//r)%r for x in el[1:]]
 rows_M,rows_S=[],[]
 for ri in cur:
  a,b,vf,vg=rels[ri]
  row=[0]*N_2core
  for k,e in vf.items():
   row[cmap[k]]=e%r
  for k,e in vg.items():
   row[cmap[nf_cnt+k]]=(-e)%r
  rows_M.append(row)
  rows_S.append(smap(a,b,f_c,eps_f)[1:4]+smap(a,b,g_c,eps_g)[1:2])
 m_rows=[row[:]for row in rows_M];s_rows=[row[:]for row in rows_S]
 alive_r=[True]*len(cur);alive_c=[True]*N_2core
 merge_stack=[]
 while True:
  piv_c=next((c for c in range(N_2core)if alive_c[c]and sum(1 for i in range(len(cur))if alive_r[i]and m_rows[i][c]!=0)==2),-1)
  if piv_c<0:break
  r_pair=[i for i in range(len(cur))if alive_r[i]and m_rows[i][piv_c]!=0]
  rp,ro=r_pair[0],r_pair[1]
  merge_stack.append((piv_c,m_rows[rp][:],s_rows[rp][:]))
  fac=(m_rows[ro][piv_c]*pow(m_rows[rp][piv_c],-1,r))%r
  m_rows[ro]=[(m_rows[ro][j]-fac*m_rows[rp][j])%r for j in range(N_2core)]
  s_rows[ro]=[(s_rows[ro][j]-fac*s_rows[rp][j])%r for j in range(4)]
  alive_r[rp]=False;alive_c[piv_c]=False
 core_cols=[c for c in range(N_2core)if alive_c[c]]
 N,D=len(core_cols),len(core_cols)+4
 active=[i for i in range(len(cur))if alive_r[i]]
 cc={c:sum(1 for i in active if m_rows[i][c]!=0)for c in core_cols}
 while len(active)>D:
  cands=[i for i in active if all(cc[c]>=3 for c in core_cols if m_rows[i][c]!=0)]or[i for i in active if all(cc[c]>=2 for c in core_cols if m_rows[i][c]!=0)]
  drop=max(cands,key=lambda i:(sum(1 for c in core_cols if m_rows[i][c]!=0),i))
  active.remove(drop)
  for c in core_cols:
   if m_rows[drop][c]!=0:cc[c]-=1
 def rref_null(mat):
  mr=[row[:]for row in mat]
  R,C=len(mr),len(mr[0])
  pivs,rk=[],0
  for c in range(C):
   piv=next((i for i in range(rk,R)if mr[i][c]%r!=0),-1)
   if piv<0:continue
   mr[rk],mr[piv]=mr[piv],mr[rk]
   inv=pow(mr[rk][c],-1,r)
   rr=mr[rk]
   for j in range(c,C):
    rr[j]=(rr[j]*inv)%r
   for i in range(R):
    if i!=rk and mr[i][c]:
     fv=mr[i][c]
     ri_r=mr[i]
     for j in range(c,C):
      ri_r[j]=(ri_r[j]-fv*rr[j])%r
   pivs.append(c)
   rk+=1
  pset=set(pivs)
  vecs=[]
  for fc in [c for c in range(C)if c not in pset]:
   v=[0]*C
   v[fc]=1
   for ri_i,pc in enumerate(pivs):
    v[pc]=(-mr[ri_i][fc])%r
   vecs.append(v)
  return rk,vecs
 M_sq=[[m_rows[i][c]for c in core_cols]for i in active]
 S_sq=[s_rows[i][:]for i in active]
 A_sq=[M_sq[i]+S_sq[i]for i in range(D)]
 rk_M,_=rref_null(M_sq)
 rk_A,ker_A=rref_null(A_sq)
 rk_3sm,ker_3sm=rref_null([M_sq[i]+S_sq[i][:3]for i in range(D)])
 assert rk_M==N and rk_A==D-1 and len(ker_A)==1 and len(ker_3sm)==0
 def mmul(P,Q):
  R,K,C=len(P),len(P[0]),len(Q[0])
  out=[[0]*C for _ in range(R)]
  for i in range(R):
   Pi,oi=P[i],out[i]
   for k in range(K):
    v=Pi[k]
    if v:
     Qk=Q[k]
     for j in range(C):
      oi[j]=(oi[j]+v*Qk[j])%r
  return out
 def mvec(P,v):
  return[sum(P[i][j]*v[j]for j in range(len(v)))%r for i in range(len(P))]
 rng=random.Random(20261119)
 s_sm=4;m_bw=4;n_bw=8
 B_sq=[M_sq[i]+[0]*s_sm for i in range(D)]
 idx_c=list(range(N));rng.shuffle(idx_c)
 blk_c=[[]for _ in range(m_bw)]
 for j,c in enumerate(idx_c):blk_c[j%m_bw].append(c)
 for t in range(s_sm):blk_c[t%m_bw].append(N+t)
 X_good=[[rng.randrange(1,r)if c in set(b)else 0 for c in range(D)]for b in blk_c]
 Z_mat=[[rng.randrange(1,r)for _ in range(n_bw-s_sm)]for _ in range(N)]
 Y_st=[S_sq[i]+mmul(M_sq,Z_mat)[i]for i in range(D)]
 hr=-(-D//m_bw)+1+-(-D//(128*m_bw))
 hc=-(-D//n_bw)+1+-(-D//(128*n_bw))
 BY,seq0=[row[:]for row in Y_st],[]
 for _ in range(hr+hc+2):
  seq0.append(mmul(X_good,BY))
  BY=mmul(B_sq,BY)
 H_rows=[[v for j in range(hc)for v in seq0[i+j][ri_b]]for i in range(hr)for ri_b in range(m_bw)]
 rk_H,ker_H=rref_null(H_rows)
 v_rec=None
 for f_vec in ker_H:
  if not any(f_vec[:s_sm]):continue
  c_rec=f_vec[:s_sm];e_rec=f_vec[s_sm:n_bw]
  acc=[0]*D
  for j in range(hc-1,0,-1):
   fj=f_vec[j*n_bw:(j+1)*n_bw]
   yc=[sum(Y_st[ri_i][bi]*fj[bi]for bi in range(n_bw))%r for ri_i in range(D)]
   acc=[(a_i+y_i)%r for a_i,y_i in zip(mvec(B_sq,acc),yc)]
  w_r=[(acc[j]+sum(Z_mat[j][t]*e_rec[t]for t in range(n_bw-s_sm)))%r for j in range(N)]
  cand=w_r+c_rec
  if any(cand)and all(x==0 for x in mvec(A_sq,cand)):
   v_rec=cand;break
 assert v_rec is not None
 O_rows,cX=[],[row[:]for row in X_good]
 for _ in range(hr):
  O_rows.extend(cX)
  cX=mmul(cX,B_sq)
 rk_OA,_=rref_null(mmul(O_rows,A_sq))
 C_cols,cY=[],mmul(B_sq,Y_st)
 for _ in range(hc-1):
  for bi in range(n_bw):
   C_cols.append([cY[ri_i][bi]for ri_i in range(D)])
  cY=mmul(B_sq,cY)
 rk_Cd,_=rref_null(C_cols)
 X_bad=[[((k+1)*X_good[0][j])%r for j in range(D)]for k in range(m_bw)]
 O_bad,cXb=[],[row[:]for row in X_bad]
 for _ in range(hr):
  O_bad.extend(cXb)
  cXb=mmul(cXb,B_sq)
 rk_OA_bad,_=rref_null(mmul(O_bad,A_sq))
 vlog_full={c:v_rec[j]for j,c in enumerate(core_cols)}
 sm_rec=v_rec[N:]
 for piv_c,pm,ps in reversed(merge_stack):
  acc=(sum(pm[j]*vlog_full[j]for j in range(N_2core)if j!=piv_c and pm[j]!=0)+sum(ps[s]*sm_rec[s]for s in range(4)))%r
  vlog_full[piv_c]=(-acc*pow(pm[piv_c],-1,r))%r
 while True:
  prog=False
  for idx in range(len(cur)):
   unk=[j for j in range(N_2core)if rows_M[idx][j]!=0 and j not in vlog_full]
   if len(unk)==1:
    cn=unk[0]
    acc=(sum(rows_M[idx][j]*vlog_full[j]for j in range(N_2core)if j!=cn and rows_M[idx][j]!=0)+sum(rows_S[idx][s]*sm_rec[s]for s in range(4)))%r
    vlog_full[cn]=(-acc*pow(rows_M[idx][cn],-1,r))%r
    prog=True
  if not prog:break
 assert len(vlog_full)==N_2core
 gen=next(g for g in range(2,p)if pow(g,(p-1)//r,p)!=1 and pow(g,2,p)!=1)
 g_r,dlog,ag=pow(gen,(p-1)//r,p),{},1
 for e in range(r):
  dlog[ag]=e
  ag=(ag*g_r)%p
 ratios=set()
 for idx,ri in enumerate(cur):
  a,b,vf,vg=rels[ri]
  lhs=(sum(e*vlog_full[cmap[k]]for k,e in vf.items())+sum(rows_S[idx][j]*sm_rec[j]for j in range(3)))%r
  rhs=(sum(e*vlog_full[cmap[nf_cnt+k]]for k,e in vg.items())-rows_S[idx][3]*sm_rec[3])%r
  assert lhs==rhs
  dl=dlog[pow((a-b*m)%p,(p-1)//r,p)]
  if dl!=0:ratios.add((lhs*pow(dl,-1,r))%r)
 pm_res=pm_basis_verify(lambda v:mvec(B_sq,v),lambda v:mvec(A_sq,v),seq0,Y_st,Z_mat,r,hr,hc,s_sm,v_rec)
 assert rk_H==D-1 and rk_OA==D-1 and rk_Cd==N and rk_OA_bad<D-1 and len(ratios)==1
 return dict(p=p,r=r,raw_rels=len(rels),purged_rels=len(cur),N_2core=N_2core,n_merged_cols=len(merge_stack),N=N,D=D,rank_M=rk_M,rank_A=rk_A,ker_dim_A=len(ker_A),ker_dim_3sm_ctrl=len(ker_3sm),rank_Cd=rk_Cd,rank_OeA=rk_OA,rank_H=rk_H,rank_OeA_1block_ctrl=rk_OA_bad,**pm_res,unplanted_2core_backsub_exact=True,horner_all_purged_rels_ok=True,distinct_Fp_dlp_ratios=len(ratios),checks_ok=True)
if __name__=="__main__":main()
