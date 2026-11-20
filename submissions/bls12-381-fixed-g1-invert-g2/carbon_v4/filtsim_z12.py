"""filtsim_z12.py — Pure-Python-stdlib Z[zeta_12] filtering & product-form simulator for BLS12-381 CM-STNFS..."""
import bisect
import collections
import heapq
import json
import math
import random
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
for _p in (str(ROOT),str(HERE)):
 if _p not in sys.path:sys.path.insert(0,_p)
from cm_models import LN2,primes_upto,rho_exact
Q_HX=1000033
CHI_HX=385966
U12=((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,0,1,0),(0,-1,0,1),(-1,0,0,0),(0,-1,0,0),(0,0,-1,0),(0,0,0,-1),(1,0,-1,0),(0,1,0,-1),)
U12_SET=set(U12)
U12_INV={U12[k]:U12[(-k)%12]for k in range(12)}
def z12_mul(a,b):
 a0,a1,a2,a3=a
 b0,b1,b2,b3=b
 c0=a0*b0
 c1=a0*b1+a1*b0
 c2=a0*b2+a1*b1+a2*b0
 c3=a0*b3+a1*b2+a2*b1+a3*b0
 c4=a1*b3+a2*b2+a3*b1
 c5=a2*b3+a3*b2
 c6=a3*b3
 return(c0-c4-c6,c1-c5,c2+c4,c3+c5)
def z12_add(a,b):
 return(a[0]+b[0],a[1]+b[1],a[2]+b[2],a[3]+b[3])
def z12_sub(a,b):
 return(a[0]-b[0],a[1]-b[1],a[2]-b[2],a[3]-b[3])
def z12_is_zero(a):
 return a==(0,0,0,0)
def z12_maxabs(a):
 return max(abs(a[0]),abs(a[1]),abs(a[2]),abs(a[3]))
def z12_eval_mod(v,chi=CHI_HX,q=Q_HX):
 c2=(chi*chi)%q
 c3=(c2*chi)%q
 return(v[0]+v[1]*chi+v[2]*c2+v[3]*c3)%q
def _poisson(rng,lam):
 if lam<=0.0:return 0
 if lam<25.0:
  L=math.exp(-lam)
  p=1.0
  k=0
  while True:
   p*=rng.random()
   if p<=L:return k
   k+=1
 return max(0,round(rng.gauss(lam,math.sqrt(lam))))
def run_filtsim_z12(log2U=12,seed=1,lf=71.29,lg=71.69,nf=548.9,ng=677.8,excess=1.09,W_target=200,s_sm=4,w_mid_scale=1.0,planted_kernel=None):
 rng=random.Random(seed)
 U_samp=1<<log2U
 uf,ug=nf/lf,ng/lg
 LiBf=(2.0**lf)/(lf*LN2)
 LiBg=(2.0**lg)/(lg*LN2)
 LiB_sum=LiBf+LiBg
 Up_f=int(round(U_samp*LiBf/LiB_sum))
 Up_g=U_samp-Up_f
 n_rels=int(round(excess*U_samp))
 t_dense=0.80
 lsmall=min(t_dense*min(lf,lg),max(8.0,math.log2(max(16.0,n_rels/25.0))+math.log2(6.0)))
 ps_all=primes_upto(int(math.floor(2.0**lsmall)))
 ps_side=ps_all[: :12]
 n_small_tot=2*len(ps_side)
 rf0,rg0=rho_exact(uf),rho_exact(ug)
 kap_f=lambda t:rho_exact(uf-t)/rf0 if rf0>0 else 0.0
 kap_g=lambda t:rho_exact(ug-t)/rg0 if rg0>0 else 0.0
 nt=1000
 ts=[t_dense+(1.0-t_dense)*(i+0.5)/nt for i in range(nt)]
 cum_f,cum_g=[],[]
 tf=tg=0.0
 for t in ts:
  tf+=(2.0**(lf*(t-1.0)))/t
  tg+=(2.0**(lg*(t-1.0)))/t
  cum_f.append(tf)
  cum_g.append(tg)
 cum_f=[x/tf for x in cum_f]
 cum_g=[x/tg for x in cum_g]
 col_lams=[]
 for p in ps_side:
  col_lams.append(n_rels*min(0.99,12.0*kap_f(math.log2(p)/lf)/p))
  col_lams.append(n_rels*min(0.99,12.0*kap_g(math.log2(p)/lg)/p))
 dt_step=(1.0-t_dense)/nt
 for _ in range(Up_f):
  idx=min(nt-1,bisect.bisect_left(cum_f,rng.random()))
  t=ts[idx]+(rng.random()-0.5)*dt_step
  lam=excess*LiB_sum*kap_f(t)/(2.0**(lf*t))
  col_lams.append(min(float(n_rels),lam))
 for _ in range(Up_g):
  idx=min(nt-1,bisect.bisect_left(cum_g,rng.random()))
  t=ts[idx]+(rng.random()-0.5)*dt_step
  lam=excess*LiB_sum*kap_g(t)/(2.0**(lg*t))
  col_lams.append(min(float(n_rels),lam))
 n_small=n_small_tot
 n_cols=len(col_lams)
 w_mid_0=0.0
 for side_L,side_k in ((lf,kap_f),(lg,kap_g)):
  a0,b0=lsmall/side_L,t_dense
  if b0>a0:
   for i in range(64):
    t0=a0+(b0-a0)*i/64.0
    t1=a0+(b0-a0)*(i+1)/64.0
    w_mid_0+=side_k(0.5*(t0+t1))*math.log(t1/t0)
 rows=[dict()for _ in range(n_rels)]
 row_wmid=[w_mid_0*w_mid_scale]*n_rels
 row_sm=[[rng.randrange(1,Q_HX)for _ in range(s_sm)]for _ in range(n_rels)]
 for c,lam in enumerate(col_lams):
  k=_poisson(rng,lam)
  if k==0:continue
  if k>n_rels:k=n_rels
  chosen=rng.sample(range(n_rels),k)
  for r_id in chosen:
   u=U12[rng.randrange(12)]
   prev=rows[r_id].get(c)
   if prev is None:rows[r_id][c]=u
   else:
    nv=z12_add(prev,u)
    if z12_is_zero(nv):del rows[r_id][c]
    else:rows[r_id][c]=nv
 alive_r=[bool(r)for r in rows]
 col_rows=[set()for _ in range(n_cols)]
 for r_id,r in enumerate(rows):
  if alive_r[r_id]:
   for c in r:
    col_rows[c].add(r_id)
 def purge_singletons():
  q_sing=[c for c in range(n_small,n_cols)if len(col_rows[c])==1]
  while q_sing:
   c=q_sing.pop()
   if len(col_rows[c])!=1:continue
   r_id=next(iter(col_rows[c]))
   if not alive_r[r_id]:continue
   alive_r[r_id]=False
   for c2 in list(rows[r_id].keys()):
    col_rows[c2].discard(r_id)
    if c2>=n_small and len(col_rows[c2])==1:q_sing.append(c2)
 purge_singletons()
 cols_2core={c for c in range(n_cols)if col_rows[c]}
 w_all=None
 if planted_kernel is not None:
  c_true,j_piv,w_seed=planted_kernel
  gw=random.Random(w_seed)
  w_all={c:gw.randrange(1,Q_HX)for c in sorted(cols_2core)}
  inv_cp=pow(c_true[j_piv],-1,Q_HX)
  for r_id in range(n_rels):
   if alive_r[r_id]:
    acc=sum(z12_eval_mod(vv,CHI_HX,Q_HX)*w_all[cc]for cc,vv in rows[r_id].items())
    acc+=sum(row_sm[r_id][s]*c_true[s]for s in range(s_sm)if s!=j_piv)
    row_sm[r_id][j_piv]=(-inv_cp*(acc%Q_HX))%Q_HX
 rows_2core=[(list(rows[i].items()),list(row_sm[i]))for i in range(n_rels)if alive_r[i]]
 n_cols_2core=len(cols_2core)
 pexc=max(6,int(round(0.0001*U_samp)))
 repurge_stack=[];active_cols=set(cols_2core);n_act_r=sum(alive_r)
 def try_drop_row(r_id,na):
  if not alive_r[r_id]or any(len(col_rows[cc])<2 for cc in rows[r_id]):return False,na
  alive_r[r_id]=False;rm_rows=[r_id];q1=[]
  for cc in rows[r_id]:
   col_rows[cc].discard(r_id)
   if len(col_rows[cc])==1:q1.append(cc)
  delta=[];bad=False
  while q1:
   c1=q1.pop()
   if len(col_rows[c1])==0:bad=True;break
   if len(col_rows[c1])!=1:continue
   r0=next(iter(col_rows[c1]))
   if not alive_r[r0]:bad=True;break
   d0=rows[r0]
   if any(len(col_rows[cc])<2 for cc in d0 if cc!=c1)or z12_eval_mod(d0[c1],CHI_HX,Q_HX)==0:bad=True;break
   alive_r[r0]=False;rm_rows.append(r0);delta.append((c1,list(d0.items()),list(row_sm[r0])))
   for cc in d0:
    col_rows[cc].discard(r0)
    if cc!=c1 and len(col_rows[cc])==1:q1.append(cc)
  if bad:
   for rr in reversed(rm_rows):
    alive_r[rr]=True
    for cc in rows[rr]:col_rows[cc].add(rr)
   return False,na
  repurge_stack.extend(delta)
  for c1,_,_ in delta:active_cols.discard(c1)
  return True,na-len(rm_rows)
 while n_act_r-len(active_cols)>pexc:
  par=list(range(n_rels))
  def _find(x):
   while par[x]!=x:par[x]=par[par[x]];x=par[x]
   return x
  for c in active_cols:
   cr=col_rows[c]
   if len(cr)==2:
    r1,r2=cr;a,b=_find(r1),_find(r2)
    if a!=b:par[a]=b
  comp_w={}
  for i in range(n_rels):
   if alive_r[i]:rt=_find(i);comp_w[rt]=comp_w.get(rt,0)+len(rows[i])+row_wmid[i]
  comps=sorted(comp_w.items(),key=lambda kv:-kv[1])
  need=n_act_r-len(active_cols)-pexc
  kill={rt for rt,_ in comps[:max(1,min(need//2+1,max(1,len(comps)//100)))]}
  dropped_any=False
  for i in range(n_rels):
   if alive_r[i]and _find(i) in kill:
    ok_d,n_act_r=try_drop_row(i,n_act_r)
    if ok_d:
     dropped_any=True
     if n_act_r-len(active_cols)<=pexc:break
  if not dropped_any:
   safe=sorted([(len(rows[i])+row_wmid[i],i)for i in range(n_rels)if alive_r[i]and all(len(col_rows[cc])>=2 for cc in rows[i])],reverse=True)
   if not safe:break
   for _,i in safe[:max(1,need)]:
    ok_d,n_act_r=try_drop_row(i,n_act_r)
    if ok_d:
     dropped_any=True
     if n_act_r-len(active_cols)<=pexc:break
   if not dropped_any:break
 post_peel_cols_lost=len(cols_2core-(active_cols|{c for c,_,_ in repurge_stack}))
 assert post_peel_cols_lost==0
 repurge_backsub_errors=0
 if planted_kernel is not None:
  c_true,_,_=planted_kernel;v_rp={c:w_all[c]for c in active_cols}
  for c1,ents,sm_v in reversed(repurge_stack):
   piv_q=0;acc=sum(sm_v[s]*c_true[s]for s in range(s_sm))
   for cc,vv in ents:
    vq=z12_eval_mod(vv,CHI_HX,Q_HX)
    if cc==c1:piv_q=vq
    else:acc+=vq*v_rp[cc]
   vc=(-pow(piv_q,-1,Q_HX)*(acc%Q_HX))%Q_HX
   v_rp[c1]=vc
   if vc!=w_all[c1]:repurge_backsub_errors+=1
  assert repurge_backsub_errors==0
 merge_stack=[]
 n_alive_r=sum(alive_r)
 tot_w=sum(len(rows[i])+row_wmid[i]for i in range(n_rels)if alive_r[i])
 pq=[]
 for c in range(n_small,n_cols):
  deg=len(col_rows[c])
  if 2<=deg<=24:heapq.heappush(pq,(deg,c))
 while pq and n_alive_r>24 and(tot_w/n_alive_r)<W_target:
  deg0,c=heapq.heappop(pq)
  deg=len(col_rows[c])
  if deg<2 or deg!=deg0 or deg>24:continue
  piv=None
  piv_w=1e18
  for r_id in col_rows[c]:
   if rows[r_id][c] in U12_SET:
    rw=len(rows[r_id])+row_wmid[r_id]
    if rw<piv_w:
     piv=r_id
     piv_w=rw
  if piv is None:continue
  piv_row=rows[piv]
  piv_wm=row_wmid[piv]
  piv_sm=row_sm[piv]
  inv_p=U12_INV[piv_row[c]]
  others=[r_id for r_id in col_rows[c]if r_id!=piv]
  ok=True
  new_rows=[]
  for r_id in others:
   fac=z12_mul(rows[r_id][c],inv_p)
   nr=dict(rows[r_id])
   for cc,vv in piv_row.items():
    nv=z12_sub(nr.get(cc,(0,0,0,0)),z12_mul(fac,vv))
    if z12_is_zero(nv):nr.pop(cc,None)
    else:
     if z12_maxabs(nv)>255:
      ok=False
      break
     nr[cc]=nv
   if not ok:break
   new_rows.append((r_id,nr,fac))
  if not ok:continue
  _dc={}
  for cc in piv_row:
   if cc!=c:_dc[cc]=_dc.get(cc,0)+1
  for r_id,nr,_ in new_rows:
   for cc in rows[r_id]:
    if cc!=c and cc not in nr:_dc[cc]=_dc.get(cc,0)+1
   for cc in nr:
    if cc not in rows[r_id]:_dc[cc]=_dc.get(cc,0)-1
  if any(len(col_rows[cc])<=d0 for cc,d0 in _dc.items()if d0>0):continue
  merge_stack.append((c,list(piv_row.items()),list(piv_sm)))
  alive_r[piv]=False
  n_alive_r-=1
  tot_w-=(len(piv_row)+piv_wm)
  for cc in list(piv_row.keys()):
   col_rows[cc].discard(piv)
  touched=set()
  for r_id,nr,fac in new_rows:
   old_r=rows[r_id]
   tot_w+=(len(nr)-len(old_r))+piv_wm
   row_wmid[r_id]+=piv_wm
   fac_q=z12_eval_mod(fac)
   rsm=row_sm[r_id]
   for j in range(s_sm):
    rsm[j]=(rsm[j]-fac_q*piv_sm[j])%Q_HX
   for cc in old_r:
    if cc not in nr:
     col_rows[cc].discard(r_id)
     touched.add(cc)
   for cc in nr:
    if cc not in old_r:
     col_rows[cc].add(r_id)
     touched.add(cc)
   rows[r_id]=nr
   if not nr:
    alive_r[r_id]=False
    n_alive_r-=1
    tot_w-=row_wmid[r_id]
  q_ms=[cc for cc in (touched|set(piv_row.keys()))if len(col_rows[cc])==1]
  while q_ms:
   c1=q_ms.pop()
   if len(col_rows[c1])!=1:continue
   rc=next(iter(col_rows[c1]))
   if not alive_r[rc]or any(len(col_rows[cc])<2 for cc in rows[rc]if cc!=c1):continue
   if z12_eval_mod(rows[rc][c1],CHI_HX,Q_HX)==0:continue
   merge_stack.append((c1,list(rows[rc].items()),list(row_sm[rc])))
   alive_r[rc]=False
   n_alive_r-=1
   tot_w-=(len(rows[rc])+row_wmid[rc])
   for cc in list(rows[rc].keys()):
    col_rows[cc].discard(rc)
    touched.add(cc)
    if cc!=c1 and len(col_rows[cc])==1:q_ms.append(cc)
  for cc in touched:
   if cc>=n_small:
    d=len(col_rows[cc])
    if 2<=d<=24:heapq.heappush(pq,(d,cc))
 active_r=[i for i in range(n_rels)if alive_r[i]]
 N_sw=len(active_r)
 rec_label={}
 rec_piv_row={}
 rec_piv_lbl={}
 row_max_lbl={r_id:0.0 for r_id in active_r}
 cand_cols=sorted((c for c in range(n_small,n_cols)if 2<=len(col_rows[c])<=64),key=lambda c:(len(col_rows[c]),c),)
 for c in cand_cols:
  c_rows=[r_id for r_id in col_rows[c]if alive_r[r_id]]
  if len(c_rows)<2:continue
  hi_rec=min((rec_piv_lbl[r_id]for r_id in c_rows if r_id in rec_piv_lbl),default=1.0)
  best_r=None
  best_w=1e18
  best_mid=None
  for r_id in c_rows:
   if r_id in rec_piv_lbl or rows[r_id][c]not in U12_SET:continue
   lo=row_max_lbl[r_id]
   if lo<hi_rec and(len(rows[r_id])+row_wmid[r_id])<best_w:
    mid=0.5*(lo+hi_rec)
    if lo<mid<hi_rec:best_r,best_w,best_mid=r_id,len(rows[r_id])+row_wmid[r_id],mid
  if best_r is not None:
   rec_label[c]=best_mid
   rec_piv_row[c]=best_r
   rec_piv_lbl[best_r]=best_mid
   for r2 in c_rows:
    if best_mid>row_max_lbl[r2]:row_max_lbl[r2]=best_mid
 rec_order=sorted(rec_label.keys(),key=lambda c:rec_label[c])
 rec_col_idx={c:k for k,c in enumerate(rec_order)}
 recorded_pivots=[]
 for k,c in enumerate(rec_order):
  r_id=rec_piv_row[c]
  r_dict=rows[r_id]
  for cc in r_dict:
   if cc in rec_col_idx:assert rec_col_idx[cc]<=k,"Non-triangular recorded pivot!"
  recorded_pivots.append((c,r_dict[c],[(cc,vv)for cc,vv in r_dict.items()if cc!=c],list(row_sm[r_id]),))
 live_rows=[r_id for r_id in active_r if r_id not in rec_piv_lbl]
 tot_ent=sum(len(rows[r_id])for r_id in active_r)
 unit_ent=sum(1 for r_id in active_r for v in rows[r_id].values()if v in U12_SET)
 exp_unit=unit_ent/max(tot_ent,1)
 k5_meas=exp_unit*7.80+(1.0-exp_unit)*33.0
 mean_w=sum(len(rows[r_id])+row_wmid[r_id]for r_id in active_r)/max(N_sw,1)
 live_col_set=sorted({c for r_id in active_r for c in rows[r_id]if c not in rec_col_idx})
 live_rows_full=[(list(rows[r_id].items()),list(row_sm[r_id]))for r_id in live_rows]
 known_sim=set(live_col_set)|set(rec_order)|{c for c,_,_ in merge_stack}|{c for c,_,_ in repurge_stack}
 post_merge_cols_lost=len(cols_2core-known_sim)
 return dict(log2U=log2U,U=U_samp,seed=seed,n_rels=n_rels,n_cols_2core=n_cols_2core,n_merge=len(merge_stack),repurge_stack_len=len(repurge_stack),post_peel_cols_lost=post_peel_cols_lost,post_merge_cols_lost=post_merge_cols_lost,repurge_backsub_errors=repurge_backsub_errors,n_reconstruct_unit=0,cols_unrecovered_tri=post_merge_cols_lost,unrecovered_frac=round(post_merge_cols_lost/max(1,n_cols_2core),6),N_sw=N_sw,N_l=len(live_rows),N_rec=len(recorded_pivots),N_sw_U=round(N_sw/U_samp,6),N_l_U=round(len(live_rows)/U_samp,6),mean_w=round(mean_w,2),exp_unit=round(exp_unit,5),k5_meas=round(k5_meas,4),n_live_cols=len(live_col_set),live_col_set=live_col_set,rec_order=rec_order,recorded_pivots=recorded_pivots,live_rows_full=live_rows_full,merge_stack=merge_stack,repurge_stack=repurge_stack,cols_2core=cols_2core,rows_2core=rows_2core,w_all=w_all,)
def main():
 raw=sys.stdin.read().strip()
 req=json.loads(raw)if raw else{}
 log2U=int(req.get("log2U",13))
 seed=int(req.get("seed",2))
 if not(10<=log2U<=15):raise ValueError("log2U must be in [10, 15]")
 res=run_filtsim_z12(log2U=log2U,seed=seed)
 for k in ("live_col_set","rec_order","recorded_pivots","live_rows_full","merge_stack","repurge_stack","cols_2core","rows_2core","w_all"):
  res.pop(k,None)
 res["checks_ok"]=(res["N_l"]>0 and res["N_rec"]>0 and res["exp_unit"]>0.50 and res["post_peel_cols_lost"]==0 and res["post_merge_cols_lost"]==0 and res["cols_unrecovered_tri"]==0)
 print(json.dumps(res,separators=(",",":")))
 if not res["checks_ok"]:raise SystemExit(1)
if __name__=="__main__":main()
