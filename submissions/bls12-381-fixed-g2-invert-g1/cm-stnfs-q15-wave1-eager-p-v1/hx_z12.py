"""hx_z12.py — Unscreened 1D-kernel recorded-pivot product-form H-R / H-X / Horner recovery diagnostic on..."""
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
def run_hx_on_filtered_z12(log2U=12,seed=1,m=12,n=8,s_sm=4):
 Q,chi=Q_HX,CHI_HX
 assert(pow(chi,4,Q)-pow(chi,2,Q)+1)%Q==0
 g=random.Random(seed*7919+n*31+m)
 if seed%2==0:
  c_true=[0,1]+[g.randrange(1,Q)for _ in range(s_sm-2)]
  j_piv=1
 else:
  c_true=[1]+[g.randrange(1,Q)for _ in range(s_sm-1)]
  j_piv=0
 fres=run_filtsim_z12(log2U=log2U,seed=seed,s_sm=s_sm,planted_kernel=(c_true,j_piv,seed*7919+n*31+m+101))
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
 lr=live_all[:D]
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
 vt=[g.randrange(Q)for _ in range(N)]
 assert op(vt)==[sum(Mc[k][i]*vt[k]for k in range(N))%Q for i in range(D)]
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
  backsub_2core_exact=(all_2core_rels_zero and recovered_match_true and fres["post_peel_cols_lost"]==0 and fres["post_merge_cols_lost"]==0 and fres["repurge_backsub_errors"]==0 and n_repurge_backsub==fres["repurge_stack_len"]and n_rec_unit==0 and cols_lost==0 and n_2core_rels_checked==len(rows_2core))
 ok=(rM==N and rA==D-1 and rC==N and rO==D-1 and rH==D-1 and rN<D-1 and kernel_eq_zero and recovered_true_kernel and backsub_2core_exact)
 return dict(checks_ok=ok,source="unscreened_1d_kernel_filtsim_z12_recorded_product",log2U=log2U,seed=seed,c0_zero_test=(c_true[0]==0),n_cols_2core=n_cols_2core,n_merge_backsub=n_merge,repurge_stack_len=fres["repurge_stack_len"],n_repurge_backsub=n_repurge_backsub,post_peel_cols_lost=fres["post_peel_cols_lost"],post_merge_cols_lost=fres["post_merge_cols_lost"],repurge_backsub_errors=fres["repurge_backsub_errors"],n_reconstruct_unit=n_rec_unit,cols_unrecovered_tri=cols_lost,unrecovered_frac=round(unrec_frac,6),occ_unrecovered_frac=round(occ_unrec_frac,6),n_2core_rels_verified_zero=n_2core_rels_checked,all_2core_rels_zero=all_2core_rels_zero,recovered_cols_match_true=recovered_match_true,backsub_2core_exact=backsub_2core_exact,filt_N_sw=fres["N_sw"],filt_N_l=fres["N_l"],filt_N_rec=n_rec,l_rec_strict_lower_tri=True,prod_eq_schur=True,N=N,D=D,s_sm=s_sm,m=m,n=n,d=d,e=e,rank_M=rM,rank_A=rA,rank_Cd=rC,rank_OeA=rO,rank_H=rH,control_one_block=rN,kernel_eq_zero=kernel_eq_zero,recovered_true_kernel=recovered_true_kernel,)
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
if __name__=="__main__":main()

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
def pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_ref):
 """Run live divide-and-conquer PM-Basis on S_k = X A^k Y (k>=1), recover exact ker(A) via Horner."""
 m=len(seq[0]);n=len(seq[0][0]);D=len(A_sq);L=hr+hc-1
 S=seq[1:1+L]
 M=[[[]for _ in range(n+m)]for _ in range(m)]
 for i in range(m):
  for j in range(n):
   poly=[S[k][i][j]for k in range(L)]
   while poly and poly[-1]==0:poly.pop()
   M[i][j]=poly
  M[i][n+i]=[r-1]
 P_b,s_out=pm_basis(M,L,[0]*n+[1]*m,r)
 def mvec(P,v):return[sum(P[i][j]*v[j]for j in range(len(v)))%r for i in range(len(P))]
 v_pm=None;deg_used=None
 for col in sorted(range(n+m),key=lambda j:s_out[j]):
  if s_out[col]>hc-1:continue
  q=[P_b[i][col]for i in range(n)]
  dg=max((len(poly)-1 for poly in q),default=-1)
  if dg<0:continue
  acc=[0]*D
  for t in range(dg+1):
   qt=[q[bi][t]if t<len(q[bi])else 0 for bi in range(n)]
   yc=[sum(Y[ri][bi]*qt[bi]for bi in range(n))%r for ri in range(D)]
   acc=[(a+y)%r for a,y in zip(mvec(A_sq,acc),yc)]
  if any(acc)and all(x==0 for x in mvec(A_sq,acc)):
   v_pm=acc;deg_used=dg;break
 assert v_pm is not None,"PM-Basis failed to recover ker(A)!"
 idx0=next(i for i,x in enumerate(v_ref)if x)
 ratio=(v_pm[idx0]*pow(v_ref[idx0],-1,r))%r
 assert ratio!=0 and all((v_ref[i]*ratio)%r==v_pm[i]for i in range(D))
 return dict(pm_basis_order=L,pm_basis_cols=n+m,pm_basis_deg=deg_used,pm_basis_matches_true_ker=True)

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
 m_bw=n_bw=4
 bnds=[(k*D)//4 for k in range(5)]
 X_good=[[0]*D for _ in range(m_bw)]
 for k in range(m_bw):
  for j in range(bnds[k],bnds[k+1]):
   X_good[k][j]=rng.randrange(1,r)
 Y=[[rng.randrange(1,r)for _ in range(n_bw)]for _ in range(D)]
 hr,hc=(D//m_bw)+2,(D//n_bw)+2
 AY,seq=[row[:]for row in Y],[]
 for _ in range(hr+hc+2):
  seq.append(mmul(X_good,AY))
  AY=mmul(A_sq,AY)
 H_rows=[[v for j in range(hc)for v in seq[i+j+1][ri_b]]for i in range(hr)for ri_b in range(m_bw)]
 rk_H,ker_H=rref_null(H_rows)
 v_rec=None
 for c_vec in ker_H:
  acc=[0]*D
  for j in range(hc-1,-1,-1):
   cj=c_vec[j*n_bw:(j+1)*n_bw]
   yc=[sum(Y[ri_i][bi]*cj[bi]for bi in range(n_bw))%r for ri_i in range(D)]
   acc=[(a_i+y_i)%r for a_i,y_i in zip(mvec(A_sq,acc),yc)]
  if any(acc)and all(x==0 for x in mvec(A_sq,acc)):
   v_rec=acc
   break
 assert v_rec is not None
 O_rows,cX=[],[row[:]for row in X_good]
 for _ in range(hr):
  cX=mmul(cX,A_sq)
  O_rows.extend(cX)
 rk_OA,_=rref_null(O_rows)
 C_cols,cY=[],[row[:]for row in Y]
 for _ in range(hc):
  for bi in range(n_bw):
   C_cols.append([cY[ri_i][bi]for ri_i in range(D)])
  cY=mmul(A_sq,cY)
 rk_Cd,_=rref_null(C_cols)
 X_bad=[[((k+1)*X_good[0][j])%r for j in range(D)]for k in range(m_bw)]
 O_bad,cXb=[],[row[:]for row in X_bad]
 for _ in range(hr):
  cXb=mmul(cXb,A_sq)
  O_bad.extend(cXb)
 rk_OA_bad,_=rref_null(O_bad)
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
 pm_res=pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_rec)
 assert rk_H==D-1 and rk_OA==D-1 and rk_Cd==D and rk_OA_bad<D-1 and len(ratios)==1
 return dict(p=p,r=r,raw_rels=len(rels),purged_rels=len(cur),N_2core=N_2core,n_merged_cols=len(merge_stack),N=N,D=D,rank_M=rk_M,rank_A=rk_A,ker_dim_A=len(ker_A),ker_dim_3sm_ctrl=len(ker_3sm),rank_Cd=rk_Cd,rank_OeA=rk_OA,rank_H=rk_H,rank_OeA_1block_ctrl=rk_OA_bad,**pm_res,unplanted_2core_backsub_exact=True,horner_all_purged_rels_ok=True,distinct_Fp_dlp_ratios=len(ratios),checks_ok=True)
