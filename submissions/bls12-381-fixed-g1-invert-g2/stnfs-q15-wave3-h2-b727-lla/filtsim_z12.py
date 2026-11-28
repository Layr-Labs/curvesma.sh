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
def run_filtsim_z12(log2U=12,seed=1,lf=70.79,lg=71.19,nf=547.0,ng=671.0,excess=1.09,W_target=200,s_sm=4,w_mid_scale=1.0,planted_kernel=None,w_max=4096):
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
 purge_stack_raw=[]
 def purge_singletons():
  q_sing=[c for c in range(n_small,n_cols)if len(col_rows[c])==1]
  while q_sing:
   c=q_sing.pop()
   if len(col_rows[c])!=1:continue
   r_id=next(iter(col_rows[c]))
   if not alive_r[r_id]:continue
   alive_r[r_id]=False;purge_stack_raw.append((c,r_id))
   for c2 in list(rows[r_id].keys()):
    col_rows[c2].discard(r_id)
    if c2>=n_small and len(col_rows[c2])==1:q_sing.append(c2)
 purge_singletons()
 cols_2core={c for c in range(n_cols)if col_rows[c]}
 cols_purge_all=set(cols_2core)
 for _,r_id in purge_stack_raw:cols_purge_all.update(rows[r_id].keys())
 w_all=None
 if planted_kernel is not None:
  c_true,j_piv,w_seed=planted_kernel
  gw=random.Random(w_seed)
  w_all={c:gw.randrange(1,Q_HX)for c in sorted(cols_2core)}
  for c in sorted(cols_purge_all-cols_2core):w_all[c]=gw.randrange(1,Q_HX)
  inv_cp=pow(c_true[j_piv],-1,Q_HX)
  for r_id in range(n_rels):
   if alive_r[r_id]:
    acc=sum(z12_eval_mod(vv,CHI_HX,Q_HX)*w_all[cc]for cc,vv in rows[r_id].items())
    acc+=sum(row_sm[r_id][s]*c_true[s]for s in range(s_sm)if s!=j_piv)
    row_sm[r_id][j_piv]=(-inv_cp*(acc%Q_HX))%Q_HX
  for _,r_id in purge_stack_raw:
   acc=sum(z12_eval_mod(vv,CHI_HX,Q_HX)*w_all[cc]for cc,vv in rows[r_id].items())
   acc+=sum(row_sm[r_id][s]*c_true[s]for s in range(s_sm)if s!=j_piv)
   row_sm[r_id][j_piv]=(-inv_cp*(acc%Q_HX))%Q_HX
 purge_stack=[(c,list(rows[r_id].items()),list(row_sm[r_id]))for c,r_id in purge_stack_raw]
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
 tot_w=sum(sum(_vinfo(vv)[0]for vv in rows[i].values())+row_wmid[i]for i in range(n_rels)if alive_r[i])
 pq=[]
 for c in range(n_small,n_cols):
  deg=len(col_rows[c])
  if 2<=deg<=24:heapq.heappush(pq,(deg,c))
 assert all((len(rows[i])+row_wmid[i])<=w_max for i in range(n_rels)if alive_r[i]),f"Initial row exceeded W_max={w_max}!"
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
   mnpt=_vinfo(fac)[0]
   nr=dict(rows[r_id])
   for cc,vv in piv_row.items():
    nv=z12_sub(nr.get(cc,(0,0,0,0)),z12_mul(fac,vv))
    if z12_is_zero(nv):nr.pop(cc,None)
    else:
     if z12_maxabs(nv)>255:
      ok=False
      break
     nr[cc]=nv
   exp_terms_nr=sum(_vinfo(vv)[0]for vv in nr.values())
   dm_mid=piv_wm if mnpt==1 else mnpt*piv_wm
   if not ok or (exp_terms_nr+row_wmid[r_id]+dm_mid)>w_max:ok=False;break
   new_rows.append((r_id,nr,fac,dm_mid))
  if not ok:continue
  _dc={}
  for cc in piv_row:
   if cc!=c:_dc[cc]=_dc.get(cc,0)+1
  for r_id,nr,_,_ in new_rows:
   for cc in rows[r_id]:
    if cc!=c and cc not in nr:_dc[cc]=_dc.get(cc,0)+1
   for cc in nr:
    if cc not in rows[r_id]:_dc[cc]=_dc.get(cc,0)-1
  if any(len(col_rows[cc])<=d0 for cc,d0 in _dc.items()if d0>0):continue
  merge_stack.append((c,list(piv_row.items()),list(piv_sm),piv_wm))
  alive_r[piv]=False
  n_alive_r-=1
  tot_w-=(sum(_vinfo(vv)[0]for vv in piv_row.values())+piv_wm)
  for cc in list(piv_row.keys()):
   col_rows[cc].discard(piv)
  touched=set()
  for r_id,nr,fac,dm_mid in new_rows:
   old_r=rows[r_id]
   tot_w+=(sum(_vinfo(vv)[0]for vv in nr.values())-sum(_vinfo(vv)[0]for vv in old_r.values()))+dm_mid
   row_wmid[r_id]+=dm_mid
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
   merge_stack.append((c1,list(rows[rc].items()),list(row_sm[rc]),row_wmid[rc]))
   alive_r[rc]=False
   n_alive_r-=1
   tot_w-=(sum(_vinfo(vv)[0]for vv in rows[rc].values())+row_wmid[rc])
   for cc in list(rows[rc].keys()):
    col_rows[cc].discard(rc)
    touched.add(cc)
    if cc!=c1 and len(col_rows[cc])==1:q_ms.append(cc)
  for cc in touched:
   if cc>=n_small:
    d=len(col_rows[cc])
    if 2<=d<=24:heapq.heappush(pq,(d,cc))
 active_r=[i for i in range(n_rels)if alive_r[i]]
 assert all((sum(_vinfo(vv)[0]for vv in rows[i].values())+row_wmid[i])<=w_max for i in active_r),f"Final row exceeded W_max={w_max}!"
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
 tot_ent=sum(_vinfo(v)[0]for r_id in active_r for v in rows[r_id].values())
 unit_ent=sum(1 for r_id in active_r for v in rows[r_id].values()for tv in _vinfo(v)[3]if abs(tv)==1)
 exp_unit=unit_ent/max(tot_ent,1)
 k5_meas=exp_unit*7.65+(1.0-exp_unit)*22.0
 mean_w=sum(sum(_vinfo(vv)[0]for vv in rows[r_id].values())+row_wmid[r_id]for r_id in active_r)/max(N_sw,1)
 live_w_sum=sum(sum(_vinfo(vv)[0]for vv in rows[r_id].values())+row_wmid[r_id]for r_id in active_r)
 piv_w_sum=sum(sum(_vinfo(vv)[0]for _,vv in ents)+pw for _,ents,_,pw in merge_stack)
 merge_stack_mean_w=(live_w_sum+piv_w_sum)/max(N_sw+len(merge_stack),1)
 live_col_set=sorted({c for r_id in active_r for c in rows[r_id]if c not in rec_col_idx})
 live_rows_full=[(list(rows[r_id].items()),list(row_sm[r_id]))for r_id in live_rows]
 known_sim=set(live_col_set)|set(rec_order)|{c for c,_,_,_ in merge_stack}|{c for c,_,_ in repurge_stack}
 post_merge_cols_lost=len(cols_2core-known_sim)
 merge_backsub_errors=0
 if planted_kernel is not None:
  c_true,_,_=planted_kernel;v_m={c:w_all[c]for c in live_col_set}
  for c1,pv,oth,sm_v in recorded_pivots:
   piv_q=z12_eval_mod(pv,CHI_HX,Q_HX);acc=sum(sm_v[s]*c_true[s]for s in range(s_sm))
   for cc,vv in oth:acc+=z12_eval_mod(vv,CHI_HX,Q_HX)*v_m[cc]
   vc=(-pow(piv_q,-1,Q_HX)*(acc%Q_HX))%Q_HX;v_m[c1]=vc
   if vc!=w_all[c1]:merge_backsub_errors+=1
  for c1,ents,sm_v,_ in reversed(merge_stack):
   piv_q=0;acc=sum(sm_v[s]*c_true[s]for s in range(s_sm))
   for cc,vv in ents:
    vq=z12_eval_mod(vv,CHI_HX,Q_HX)
    if cc==c1:piv_q=vq
    else:acc+=vq*v_m[cc]
   vc=(-pow(piv_q,-1,Q_HX)*(acc%Q_HX))%Q_HX;v_m[c1]=vc
   if vc!=w_all[c1]:merge_backsub_errors+=1
  assert merge_backsub_errors==0,f"merge_backsub_errors={merge_backsub_errors}"
 return dict(log2U=log2U,U=U_samp,seed=seed,n_rels=n_rels,n_cols_2core=n_cols_2core,n_merge=len(merge_stack),repurge_stack_len=len(repurge_stack),purge_stack_len=len(purge_stack),post_peel_cols_lost=post_peel_cols_lost,post_merge_cols_lost=post_merge_cols_lost,repurge_backsub_errors=repurge_backsub_errors,n_reconstruct_unit=0,cols_unrecovered_tri=post_merge_cols_lost,unrecovered_frac=round(post_merge_cols_lost/max(1,n_cols_2core),6),N_sw=N_sw,N_l=len(live_rows),N_rec=len(recorded_pivots),N_sw_U=round(N_sw/U_samp,6),N_l_U=round(len(live_rows)/U_samp,6),mean_w=round(mean_w,2),merge_stack_mean_w=round(merge_stack_mean_w,4),exp_unit=round(exp_unit,5),k5_meas=round(k5_meas,4),n_live_cols=len(live_col_set),live_col_set=live_col_set,rec_order=rec_order,recorded_pivots=recorded_pivots,live_rows_full=live_rows_full,merge_stack=[(c,e,s)for c,e,s,_ in merge_stack],repurge_stack=repurge_stack,purge_stack=purge_stack,cols_2core=cols_2core,rows_2core=rows_2core,w_all=w_all,)
def main():
 raw=sys.stdin.read().strip()
 req=json.loads(raw)if raw else{}
 log2U=int(req.get("log2U",13))
 seed=int(req.get("seed",2))
 if not(10<=log2U<=15):raise ValueError("log2U must be in [10, 15]")
 res=run_filtsim_z12(log2U=log2U,seed=seed)
 for k in ("live_col_set","rec_order","recorded_pivots","live_rows_full","merge_stack","repurge_stack","purge_stack","cols_2core","rows_2core","w_all"):
  res.pop(k,None)
 res["checks_ok"]=(res["N_l"]>0 and res["N_rec"]>0 and res["exp_unit"]>0.50 and res["post_peel_cols_lost"]==0 and res["post_merge_cols_lost"]==0 and res["cols_unrecovered_tri"]==0)
 print(json.dumps(res,separators=(",",":")))
 if not res["checks_ok"]:raise SystemExit(1)

import base64 as _b64, cmath as _cm, lzma as _lz
_M64=(1<<64)-1
_INV2P53=1.0/9007199254740992.0
_VCACHE={}
def _vinfo(v):
 r=_VCACHE.get(v)
 if r is not None:return r
 a0,a1,a2,a3=v
 if(a0>=0 and a2>=0)or(a0<=0 and a2<=0):t0,t2,t4=a0,a2,0
 elif abs(a0)<=abs(a2):t0,t2,t4=0,a0+a2,-a0
 else:t0,t2,t4=a0+a2,0,a2
 if(a1>=0 and a3>=0)or(a1<=0 and a3<=0):t1,t3,t5=a1,a3,0
 elif abs(a1)<=abs(a3):t1,t3,t5=0,a1+a3,-a1
 else:t1,t3,t5=a1+a3,0,a3
 t6=(t0,t1,t2,t3,t4,t5);npt=n1=mx=0
 for x in t6:
  if x:
   ax=-x if x<0 else x;npt+=1;n1+=ax
   if ax>mx:mx=ax
 r=(npt,n1,mx,t6);_VCACHE[v]=r;return r
for _u0 in U12:_vinfo(_u0)
def _z12_canon(v):
 best=v
 for u in U12:
  w=z12_mul(v,u)
  if w<best:best=w
 return w if False else best
class _Xo256SS:
 __slots__=("s0","s1","s2","s3")
 def __init__(self,seed):
  x=((seed*0x9E3779B97F4A7C15)+0x123456789ABCDEF0)&_M64;st=[]
  for _ in range(4):
   x=(x+0x9E3779B97F4A7C15)&_M64;z=((x^(x>>30))*0xBF58476D1CE4E5B9)&_M64
   z=((z^(z>>27))*0x94D049BB133111EB)&_M64;st.append((z^(z>>31))&_M64)
  self.s0,self.s1,self.s2,self.s3=st
 def u64(self):
  s0,s1,s2,s3=self.s0,self.s1,self.s2,self.s3
  t5=(s1*5)&_M64;res=((((t5<<7)|(t5>>57))&_M64)*9)&_M64;t=(s1<<17)&_M64
  s2^=s0;s3^=s1;s1^=s2;s0^=s3;s2^=t;s3=((s3<<45)|(s3>>19))&_M64
  self.s0,self.s1,self.s2,self.s3=s0,s1,s2,s3;return res
 def dbl(self):return(self.u64()>>11)*_INV2P53
 def poisson(self,lam):
  if lam<=0.0:return 0
  if lam<30.0:
   L=math.exp(-lam);p=1.0;k=0
   while True:
    p*=self.dbl()
    if p<=L:return k
    k+=1
  u1=max(1e-15,self.dbl());u2=self.dbl()
  xv=lam+math.sqrt(lam)*math.sqrt(-2.0*math.log(u1))*math.cos(6.28318530717958647692*u2)
  k=int(math.floor(xv+0.5))if xv>=0.0 else-int(math.floor(-xv+0.5))
  return 0 if k<0 else k
def run_filtsim_z12_ladder(log2U=12,seed=0,Lf=71.28,Lg=71.56,Nf=547.0,Ng=668.0,excess=1.03,pexc=0.0001,t_dense=0.80,targets=(200,250),kmax=64,cmax=255,w_max=4096,pf_s=60.5,pg_s=63.5,ret_f=1.0,ret_g=1.0,sbf=345.0,sbg=290.0):
 rng=_Xo256SS(seed);uf,ug=Nf/Lf,Ng/Lg;rf0,rg0=rho_exact(uf),rho_exact(ug)
 LiBf=math.exp(Lf*LN2)/(Lf*LN2);LiBg=math.exp(Lg*LN2)/(Lg*LN2);LiB_sum=LiBf+LiBg
 Up=1<<log2U;Up_f=int(math.floor(Up*(LiBf/LiB_sum)+0.5));Up_g=Up-Up_f;Rp=int(math.floor(excess*Up+0.5))
 lsmall=min(t_dense*min(Lf,Lg),max(8.0,math.log(Rp/25.0)/LN2+math.log(6.0)/LN2))
 ps=primes_upto(int(math.floor(math.exp(lsmall*LN2))));nprimes=len(ps);nsmall_side=(nprimes+11)//12;U=Up+2*nsmall_side
 nt=1000;ts=[0.0]*nt;cum_f=[0.0]*nt;cum_g=[0.0]*nt;tot_f=tot_g=0.0
 for i in range(nt):
  t=t_dense+(1.0-t_dense)*(i+0.5)/nt;ts[i]=t
  tot_f+=math.exp(Lf*(t-1.0)*LN2)/t;tot_g+=math.exp(Lg*(t-1.0)*LN2)/t;cum_f[i]=tot_f;cum_g[i]=tot_g
 cum_f=[x/tot_f for x in cum_f];cum_g=[x/tot_g for x in cum_g]
 col_lams=[0.0]*U;col_log2p=[0.0]*U;col_q=[0]*U;col_side=[0]*U;col_lb_uns=[0.0]*U;dt_step=(1.0-t_dense)/nt;cid=0
 for side_id,n_s,cm_s,u_s,L_s,r0_s,p_s,r_m in((0,Up_f,cum_f,uf,Lf,rf0,pf_s,ret_f),(1,Up_g,cum_g,ug,Lg,rg0,pg_s,ret_g)):
  for _ in range(n_s):
   idx=min(nt-1,bisect.bisect_left(cm_s,rng.dbl()));t=ts[idx]+(rng.dbl()-0.5)*dt_step
   lb=L_s*t;rm=r_m if lb>p_s else 1.0
   kap=(rho_exact(u_s-t)/r0_s)if r0_s>0.0 else 0.0
   lam=min(0.90*Rp,excess*LiB_sum*kap*rm/math.exp(lb*LN2))
   col_lams[cid]=lam;col_log2p[cid]=lb;col_side[cid]=side_id;col_lb_uns[cid]=lb if lb>p_s else 0.0;cid+=1
 for side_id,(u_s,L_s,r0_s) in enumerate(((uf,Lf,rf0),(ug,Lg,rg0))):
  for i in range(0,nprimes,12):
   q=ps[i];lp=math.log(float(q))/LN2;kap=(rho_exact(u_s-lp/L_s)/r0_s)if r0_s>0.0 else 0.0
   col_lams[cid]=Rp*min(0.90,12.0*kap/q);col_log2p[cid]=lp;col_q[cid]=q;col_side[cid]=side_id;col_lb_uns[cid]=0.0;cid+=1
 w_mid=w_mid_28=0.0
 for u_s,L_s,r0_s in((uf,Lf,rf0),(ug,Lg,rg0)):
  tmin=lsmall/L_s;dt=(t_dense-tmin)/256.0
  for i in range(256):
   t=tmin+dt*(i+0.5);dw=((rho_exact(u_s-t)/r0_s)if r0_s>0.0 else 0.0)*dt/t;w_mid+=dw
   if L_s*t<=28.0:w_mid_28+=dw
 Rp_raw=Rp;n_surv=Rp;max_row_bf=0.0;max_row_bg=0.0;stageb_row_keep_frac=1.0;ret_f_realized=1.0;ret_g_realized=1.0
 row_e=[[]for _ in range(Rp)]
 if sbf>0.0 or sbg>0.0:
  p_keep_est=0.55;rng_p1=(rng.s0,rng.s1,rng.s2,rng.s3)
  while True:
   rng.s0,rng.s1,rng.s2,rng.s3=rng_p1
   Rp_raw=max(Rp,int(math.ceil(Rp/p_keep_est)));scale_raw=Rp_raw/Rp
   row_bf=[0.0]*Rp_raw;row_bg=[0.0]*Rp_raw;stamp=[0]*Rp_raw;uns_buf=[];raw_f_inc=raw_g_inc=0
   for c in range(U):
    lb_c=col_lb_uns[c]
    if lb_c<=0.0:continue
    lam_r=min(0.90*Rp_raw,col_lams[c]*scale_raw)
    k=rng.poisson(lam_r)
    if k<=0:continue
    if k>Rp_raw:k=Rp_raw
    tag=c+1;rb_arr=row_bf if col_side[c]==0 else row_bg
    if col_side[c]==0:raw_f_inc+=k
    else:raw_g_inc+=k
    for _ in range(k):
     s0,s1,s2,s3=rng.s0,rng.s1,rng.s2,rng.s3
     while True:
      t5=(s1*5)&_M64;res=((((t5<<7)|(t5>>57))&_M64)*9)&_M64;t=(s1<<17)&_M64
      s2^=s0;s3^=s1;s1^=s2;s0^=s3;s2^=t;s3=((s3<<45)|(s3>>19))&_M64
      r=res%Rp_raw
      if stamp[r]!=tag:break
     rng.s0,rng.s1,rng.s2,rng.s3=s0,s1,s2,s3
     stamp[r]=tag;rb_arr[r]+=lb_c;uns_buf.append((r,c))
   n_surv=0
   for r in range(Rp_raw):
    if(sbf<=0.0 or row_bf[r]<=sbf)and(sbg<=0.0 or row_bg[r]<=sbg):n_surv+=1
   if n_surv>=Rp:
    raw_to_surv=[-1]*Rp_raw;kept=0;max_row_bf=0.0;max_row_bg=0.0
    for r in range(Rp_raw):
     if(sbf<=0.0 or row_bf[r]<=sbf)and(sbg<=0.0 or row_bg[r]<=sbg)and kept<Rp:
      raw_to_surv[r]=kept;kept+=1
      if row_bf[r]>max_row_bf:max_row_bf=row_bf[r]
      if row_bg[r]>max_row_bg:max_row_bg=row_bg[r]
    assert kept==Rp
    break
   p_keep_est*=0.5;assert p_keep_est>1e-6
  if sbf>0.0:assert max_row_bf<=sbf+1e-9,f"Stage-B f-cap violated: {max_row_bf} > {sbf}"
  if sbg>0.0:assert max_row_bg<=sbg+1e-9,f"Stage-B g-cap violated: {max_row_bg} > {sbg}"
  stageb_row_keep_frac=n_surv/Rp_raw
  stamp2=[0]*Rp;u_pos=0;u_len=len(uns_buf);surv_f_inc=surv_g_inc=0
  for c in range(U):
   q=col_q[c]
   if col_lb_uns[c]>0.0:
    while u_pos<u_len and uns_buf[u_pos][1]==c:
     sr=raw_to_surv[uns_buf[u_pos][0]];u_pos+=1
     if sr<0:continue
     if col_side[c]==0:surv_f_inc+=1
     else:surv_g_inc+=1
     u=U12[rng.u64()%12]
     if q>0:
      v=1;invq=1.0/q
      while rng.dbl()<invq:v+=1
      if v>1:u=(v*u[0],v*u[1],v*u[2],v*u[3])
     row_e[sr].append((c,u))
   else:
    k=rng.poisson(col_lams[c])
    if k<=0:continue
    if k>Rp:k=Rp
    tag=c+1;invq=(1.0/q)if q>0 else 0.0
    for _ in range(k):
     s0,s1,s2,s3=rng.s0,rng.s1,rng.s2,rng.s3
     while True:
      t5=(s1*5)&_M64;res=((((t5<<7)|(t5>>57))&_M64)*9)&_M64;t=(s1<<17)&_M64
      s2^=s0;s3^=s1;s1^=s2;s0^=s3;s2^=t;s3=((s3<<45)|(s3>>19))&_M64
      sr=res%Rp
      if stamp2[sr]!=tag:break
     t5=(s1*5)&_M64;res=((((t5<<7)|(t5>>57))&_M64)*9)&_M64;t=(s1<<17)&_M64
     s2^=s0;s3^=s1;s1^=s2;s0^=s3;s2^=t;s3=((s3<<45)|(s3>>19))&_M64
     rng.s0,rng.s1,rng.s2,rng.s3=s0,s1,s2,s3
     stamp2[sr]=tag;u=U12[res%12]
     if q>0:
      v=1
      while rng.dbl()<invq:v+=1
      if v>1:u=(v*u[0],v*u[1],v*u[2],v*u[3])
     row_e[sr].append((c,u))
  assert u_pos==u_len
  raw_f_mean=raw_f_inc/Rp_raw;surv_f_mean=surv_f_inc/Rp
  ret_f_realized=(surv_f_mean/raw_f_mean)if raw_f_mean>0.0 else 1.0
  raw_g_mean=raw_g_inc/Rp_raw;surv_g_mean=surv_g_inc/Rp
  ret_g_realized=(surv_g_mean/raw_g_mean)if raw_g_mean>0.0 else 1.0
 else:
  stamp=[0]*Rp
  for c in range(U):
   k=rng.poisson(col_lams[c])
   if k<=0:continue
   if k>Rp:k=Rp
   q=col_q[c];tag=c+1;invq=(1.0/q)if q>0 else 0.0
   for _ in range(k):
    while True:
     r=rng.u64()%Rp
     if stamp[r]!=tag:break
    stamp[r]=tag;u=U12[rng.u64()%12]
    if q>0:
     v=1
     while rng.dbl()<invq:v+=1
     if v>1:u=(v*u[0],v*u[1],v*u[2],v*u[3])
    row_e[r].append((c,u))
 row_alive=[bool(e)for e in row_e];row_rw=[0.0]*Rp;colw=[0]*U;crows=[[]for _ in range(U)];nalive=0
 for i in range(Rp):
  if not row_alive[i]:continue
  nalive+=1;w=w_mid
  for c,v in row_e[i]:w+=_vinfo(v)[0];colw[c]+=1;crows[c].append(i)
  assert w<=w_max,f"Initial row {i} weight {w} exceeded W_max={w_max}!"
  row_rw[i]=w
 sq=[c for c in range(U)if colw[c]==1];ncols=sum(1 for w in colw if w>0)
 def peel_bfs():
  nonlocal nalive,ncols
  qh=0
  while qh<len(sq):
   c=sq[qh];qh+=1
   if colw[c]!=1:continue
   for ri in crows[c]:
    if row_alive[ri]:
     row_alive[ri]=False;nalive-=1
     for cc,_ in row_e[ri]:
      colw[cc]-=1
      if colw[cc]==0:ncols-=1
      elif colw[cc]==1:sq.append(cc)
     break
  sq.clear()
 peel_bfs();ncols_2core_init=ncols;n_repurge_pivs=0;texc=6 if pexc<=0.0 else max(2,int(pexc*Up));act=[i for i in range(Rp)if row_alive[i]];par=list(range(Rp));col_first=[-1]*U
 def uf_find(x):
  while par[x]!=x:par[x]=par[par[x]];x=par[x]
  return x
 def safe_peel_bfs():
  nonlocal nalive,ncols,n_repurge_pivs
  qh=0
  while qh<len(sq):
   c=sq[qh];qh+=1
   if colw[c]!=1:continue
   for ri in crows[c]:
    if row_alive[ri]:
     ok_r=True
     for cc,cv in row_e[ri]:
      if cc==c:
       if cv not in U12_SET:ok_r=False;break
      elif colw[cc]<2:ok_r=False;break
     if not ok_r:break
     row_alive[ri]=False;nalive-=1;n_repurge_pivs+=1
     for cc,_ in row_e[ri]:
      colw[cc]-=1
      if colw[cc]==0:ncols-=1
      elif colw[cc]==1:sq.append(cc)
     break
  sq.clear()
 while nalive-ncols>texc:
  act=[i for i in act if row_alive[i]]
  for i in act:par[i]=i
  touched=[]
  for i in act:
   for c,_ in row_e[i]:
    if colw[c]==2:
     cf=col_first[c]
     if cf<0:col_first[c]=i;touched.append(c)
     else:
      u0,v0=uf_find(cf),uf_find(i)
      if u0!=v0:par[u0]=v0
  for c in touched:col_first[c]=-1
  comp_w={}
  for i in act:
   ok2=True
   for c,_ in row_e[i]:
    if colw[c]<2:ok2=False;break
   if not ok2:continue
   rt=uf_find(i);comp_w[rt]=comp_w.get(rt,0.0)+row_rw[i]
  if not comp_w:break
  comps=sorted(comp_w.items(),key=lambda kv:(-kv[1],kv[0]))
  need=nalive-ncols-texc;nb=max(1,min(need//2+1,max(1,len(comps)//100)));kill={rt for rt,_ in comps[:nb]}
  dropped_any=False
  for i in act:
   if not row_alive[i]or uf_find(i)not in kill:continue
   ok2=True
   for c,_ in row_e[i]:
    if colw[c]<2:ok2=False;break
   if not ok2:continue
   row_alive[i]=False;nalive-=1;dropped_any=True
   for cc,_ in row_e[i]:
    colw[cc]-=1
    if colw[cc]==1:sq.append(cc)
   safe_peel_bfs()
   if nalive-ncols<=texc:break
  if not dropped_any:
   for i in act:
    if not row_alive[i]:continue
    ok2=True
    for c,_ in row_e[i]:
     if colw[c]<2:ok2=False;break
    if not ok2:continue
    row_alive[i]=False;nalive-=1;dropped_any=True
    for cc,_ in row_e[i]:
     colw[cc]-=1
     if colw[cc]==1:sq.append(cc)
    safe_peel_bfs()
    if nalive-ncols<=texc:break
   if not dropped_any:break
 post_peel_cols_lost=ncols_2core_init-(ncols+n_repurge_pivs);n_merge_pivs=0;merge_piv_w_sum=0.0
 r_e=[row_e[i]for i in range(Rp)if row_alive[i]];nrows_tot=len(r_e);r_alive=[True]*nrows_tot
 r_rw=[row_rw[i]for i in range(Rp)if row_alive[i]];r_wmid_t=[w_mid]*nrows_tot;r_wmid_c=[w_mid]*nrows_tot
 r_mid=[[]for _ in range(nrows_tot)];r_imul=[1]*nrows_tot;cols=[[]for _ in range(U)]
 for i in range(nrows_tot):
  for c,v in r_e[i]:cols[c].append([i,v in U12_SET])
 def col_del(c,row):
  cl=cols[c]
  for idx in range(len(cl)):
   if cl[idx][0]==row:cl[idx]=cl[-1];cl.pop();return
 def col_upd(c,row,is_u):
  cl=cols[c]
  for idx in range(len(cl)):
   if cl[idx][0]==row:cl[idx][1]=is_u;return
 def col_key(c):
  cl=cols[c];wgt=len(cl)
  if wgt<=2:return float(wgt-3)
  best=1e18
  for ri,is_u in cl:
   if is_u and r_rw[ri]<best:best=r_rw[ri]
  return float(wgt-2)*best
 def rfind(ents,col):
  lo,hi=0,len(ents)-1
  while lo<=hi:
   mid=(lo+hi)>>1;cc=ents[mid][0]
   if cc==col:return ents[mid][1]
   if cc<col:lo=mid+1
   else:hi=mid-1
  return None
 def mid_add(ml,ck,x):
  if x<=0.0:return
  for it in ml:
   if it[0]==ck:it[1]+=x;return
  ml.append([ck,x])
 def mid_combine(di,pi,m):
  mu=m in U12_SET;pm=r_mid[pi]
  if mu and not pm:return
  dm=r_mid[di];pnu=sum(it[1]for it in pm)
  if not mu:
   pu=r_wmid_c[pi]-pnu
   if pu>0.0:mid_add(dm,_z12_canon(m),pu)
  for pk,pc in pm:
   nk=pk if mu else _z12_canon(z12_mul(m,pk))
   if nk not in U12_SET:mid_add(dm,nk,pc)
 def eval_snap(W_tgt,maxc_cur,sseed):
  U5,G1=7.65,22.0;alive_idx=[i for i in range(nrows_tot)if r_alive[i]];N_sw=len(alive_idx)
  _vi=_vinfo;_U12S=U12_SET
  tot_rw=sum(r_rw[i]for i in alive_idx);crs=[[]for _ in range(U)]
  for a,ri in enumerate(alive_idx):
   for c,_ in r_e[ri]:crs[c].append(a)
  order=sorted((len(crs[c]),c)for c in range(U)if crs[c]);dc=[-1]*N_sw;lab=[0.0]*U;mn=[1.0]*U;hl=[False]*U;nrec=0
  for _,c in order:
   if c>=Up:continue
   hi=mn[c];cands=[]
   for a in crs[c]:
    if dc[a]>=0:continue
    ri=alive_idx[a];vv=rfind(r_e[ri],c)
    if vv is not None and vv in _U12S:cands.append((r_rw[ri],a))
   cands.sort()
   for _,a in cands:
    ri=alive_idx[a];lo=0.0
    for x,_ in r_e[ri]:
     if x!=c and hl[x]and lab[x]>lo:lo=lab[x]
    v=lo+(hi-lo)*0.5
    if lo<v<hi:
     dc[a]=c;lab[c]=v;hl[c]=True;nrec+=1
     for x,_ in r_e[ri]:
      if x!=c and mn[x]>v:mn[x]=v
     break
  recs=sorted((lab[dc[a]],a)for a in range(N_sw)if dc[a]>=0);col_dp=[0]*U;max_dp=0
  for _,a in recs:
   c=dc[a];ri=alive_idx[a];mx=0
   for x,_ in r_e[ri]:
    if x!=c and col_dp[x]>mx:mx=col_dp[x]
   col_dp[c]=mx+1
   if mx+1>max_dp:max_dp=mx+1
  et=eo=ec=es4=it=io=ic=is4=tp=tg=0.0;srng=_Xo256SS((sseed^0xABCDEF1234567890)&_M64);m28f=(w_mid_28/w_mid)if w_mid>0.0 else 0.0
  for a,ri in enumerate(alive_idx):
   pc=dc[a];nb4=nb3=nb2=0;rp=0.0
   for c,val in r_e[ri]:
    t6=_vi(val)[3];lp=col_log2p[c];ip=(c==pc)
    for tv in t6:
     if not tv:continue
     v=-tv if tv<0 else tv;et+=1.0
     if v==1:eo+=1.0
     elif v<=4:es4+=1.0
     tc=(U5*v)if v<=4 else G1;ec+=tc
     if not ip:
      rp+=tc
      if v<=4 and lp<=28.0:
       if lp<=14.0:nb4+=v
       elif lp<=18.25:nb3+=v
       else:nb2+=v
   nuc=sum(x[1]for x in r_mid[ri]);uc=max(0.0,r_wmid_c[ri]-nuc);io+=uc;it+=uc;ic+=uc*U5;rp+=uc*U5
   if m28f>0.0 and uc>0.0:nb2+=srng.poisson(uc*m28f)
   for mk,xc in r_mid[ri]:
    t6=_vi(mk)[3];npm=0;vc=0.0
    for tv in t6:
     if not tv:continue
     v=-tv if tv<0 else tv;npm+=1;vc+=(U5*v)if v<=4 else G1
     if 1<v<=4:is4+=xc
     if v<=4 and m28f>0.0:nb2+=v*srng.poisson(xc*m28f)
    it+=xc*npm;ic+=xc*vc;rp+=xc*vc
   g4,r4=divmod(nb4,4);g3,r3=divmod(nb3+r4,3);g2=(nb2+r3)//2
   tp+=rp;tg+=rp-(3*g4+2*g3+g2)*U5
  tt=et+it;to=eo+io;tc_all=ec+ic;max_rw_s=max((r_rw[i]for i in alive_idx),default=0.0)
  assert tt<=tot_rw+1e-6,f"Expanded term count {tt} exceeded upper bound tot_rw={tot_rw}!"
  n_ret=N_sw+n_merge_pivs
  ms_mw=(tot_rw+merge_piv_w_sum)/max(n_ret,1)
  mp_mw=(merge_piv_w_sum/n_merge_pivs)if n_merge_pivs>0 else 0.0
  return dict(W_target=W_tgt,N_sw=N_sw,N_l=N_sw-nrec,depth=max_dp,avg_w=tt/max(N_sw,1),max_row_w=max_rw_s,merge_stack_mean_w=ms_mw,merge_piv_mean_w=mp_mw,unit_frac=(to/tt)if tt>0 else 1.0,exp_unit=(eo/et)if et>0 else 1.0,small4_nu_frac=((es4+is4)/max(tt-to,1e-9)),k5=(tc_all/tt)if tt>0 else U5,tab_ratio=(tg/tp)if tp>0 else 1.0)
 nrows=nrows_tot;totw=sum(r_rw);heap=[]
 _heappush=heapq.heappush;_heappop=heapq.heappop;_vi=_vinfo;_U12S=U12_SET;_U12I=U12_INV;_zmul=z12_mul;_zsub=z12_sub
 for c in range(U):
  if 0<len(cols[c])<=kmax:_heappush(heap,(col_key(c),c))
 snaps=[];ti=0;ntg=len(targets);maxc=1
 while heap and ti<ntg:
  ik,c=_heappop(heap);cl=cols[c];wgt=len(cl)
  if wgt==0 or wgt>kmax:continue
  ck=col_key(c)
  if ck>ik+1e-9:_heappush(heap,(ck,c));continue
  if wgt==1:
   if not cl[0][1]:continue
   ri=cl[0][0]
   if any(len(cols[cc])<2 for cc,_ in r_e[ri]if cc!=c):continue
   n_merge_pivs+=1;merge_piv_w_sum+=r_rw[ri];r_alive[ri]=False;nrows-=1;totw-=r_rw[ri]
   for cc,_ in r_e[ri]:
    col_del(cc,ri)
    if 0<len(cols[cc])<=kmax:_heappush(heap,(col_key(cc),cc))
   r_e[ri]=[]
   # Snapshot threshold check after singleton peel (aligned with C ladder generator; zero observed col loss):
   avg=totw/max(nrows,1)
   while ti<ntg and avg>=targets[ti]:
    assert max((r_rw[j]for j in range(nrows_tot)if r_alive[j]),default=0.0)<=w_max
    snaps.append(eval_snap(targets[ti],maxc,seed+targets[ti]));ti+=1
   continue
  piv=-1;brw=1e18
  for ri,iu in cl:
   if iu and(r_rw[ri]<brw or(r_rw[ri]==brw and ri<piv)):brw=r_rw[ri];piv=ri
  if piv<0:continue
  pre=r_e[piv];pc_inv=_U12I[rfind(pre,c)];prl=len(pre);pri=r_imul[piv];ok=True
  _dc={cb:1 for cb,_ in pre if cb!=c}
  for ri,_ in cl:
   if ri==piv:continue
   rei=r_e[ri];m=_zmul(rfind(rei,c),pc_inv);vim=_vi(m)
   if r_imul[ri]+vim[1]*pri>cmax or r_rw[ri]+vim[0]*r_rw[piv]>w_max:ok=False;break
   a=0;rel=len(rei)
   for cb,pv in pre:
    while a<rel and rei[a][0]<cb:a+=1
    has_a=(a<rel and rei[a][0]==cb)
    ov=rei[a][1]if has_a else(0,0,0,0)
    nv=_zsub(ov,_zmul(m,pv))
    if _vi(nv)[2]>cmax:ok=False;break
    if cb!=c:
     if has_a and nv==(0,0,0,0):_dc[cb]=_dc.get(cb,0)+1
     elif not has_a:_dc[cb]=_dc.get(cb,0)-1
   if not ok:break
  if ok and any(len(cols[cb])<=d0 for cb,d0 in _dc.items()if d0>0):ok=False
  if not ok:continue
  n_merge_pivs+=1;merge_piv_w_sum+=r_rw[piv]
  csnap=[x[0]for x in cl]
  for ri in csnap:
   if ri==piv:continue
   rei=r_e[ri];orw=r_rw[ri];m=_zmul(rfind(rei,c),pc_inv);mnpt,mn1,_,_= _vi(m)
   r_imul[ri]+=mn1*pri;maxc=max(maxc,r_imul[ri]);mid_combine(ri,piv,m)
   dm=r_wmid_t[piv]if mnpt==1 else mnpt*r_wmid_t[piv];r_wmid_t[ri]+=dm;r_wmid_c[ri]+=r_wmid_c[piv];totw+=dm
   oe=[];a=b=0;rel=len(rei);ntm=0
   while a<rel or b<prl:
    if b==prl or(a<rel and rei[a][0]<pre[b][0]):
     ent=rei[a];oe.append(ent);ntm+=_vi(ent[1])[0];a+=1
    elif a==rel or pre[b][0]<rei[a][0]:
     cc,pv=pre[b];mv=_zmul(m,pv);nv=(-mv[0],-mv[1],-mv[2],-mv[3]);vi=_vi(nv)
     oe.append((cc,nv));cols[cc].append([ri,nv in _U12S]);totw+=vi[0];ntm+=vi[0];maxc=max(maxc,vi[2]);b+=1
    else:
     cc,ov=rei[a];_,pv=pre[b];nv=_zsub(ov,_zmul(m,pv));no=_vi(ov)[0]
     if nv==(0,0,0,0):
      totw-=no;col_del(cc,ri)
      if 0<len(cols[cc])<=kmax:_heappush(heap,(col_key(cc),cc))
     else:
      vi=_vi(nv);totw+=vi[0]-no;ntm+=vi[0];oe.append((cc,nv));col_upd(cc,ri,nv in _U12S);maxc=max(maxc,vi[2])
     a+=1;b+=1
   r_e[ri]=oe;nrw=r_wmid_t[ri]+ntm;r_rw[ri]=nrw
   if nrw<orw:
    for cc,_ in oe:
     if len(cols[cc])<=kmax:_heappush(heap,(col_key(cc),cc))
  r_alive[piv]=False;nrows-=1;totw-=r_rw[piv]
  for cc,_ in pre:
   if cc==c:continue
   col_del(cc,piv)
   if 0<len(cols[cc])<=kmax:_heappush(heap,(col_key(cc),cc))
  cols[c].clear();r_e[piv]=[]
  avg=totw/max(nrows,1)
  while ti<ntg and avg>=targets[ti]:
   assert max((r_rw[j]for j in range(nrows_tot)if r_alive[j]),default=0.0)<=w_max
   snaps.append(eval_snap(targets[ti],maxc,seed+targets[ti]));ti+=1
 while ti<ntg:
  assert max((r_rw[j]for j in range(nrows_tot)if r_alive[j]),default=0.0)<=w_max
  snaps.append(eval_snap(targets[ti],maxc,seed+targets[ti]));ti+=1
 s200=snaps[0];s250=snaps[1]if len(snaps)>1 else s200
 max_rw=s200["max_row_w"]
 assert max_rw<=w_max,f"Row weight {max_rw} exceeded W_max={w_max}!"
 # End-of-merge column-loss accounting (identical zero-loss invariant to C ladder generator):
 ncols_final=sum(1 for cl0 in cols if cl0)
 post_merge_cols_lost=ncols_2core_init-(ncols_final+n_repurge_pivs+n_merge_pivs)
 return dict(log2U=log2U,seed=seed,U=Up,U_tot=U,w_mid=round(w_mid,4),N_sw=s200["N_sw"],N_l=s200["N_l"],N_sw_250=s250["N_sw"],depth=s200["depth"],avg_w=round(s200["avg_w"],4),max_row_w=round(max_rw,2),w_max_guard=w_max,merge_stack_mean_w=round(s200["merge_stack_mean_w"],4),merge_piv_mean_w=round(s200["merge_piv_mean_w"],4),unit_frac=round(s200["unit_frac"],6),exp_unit=round(s200["exp_unit"],6),small4_nu_frac=round(s200["small4_nu_frac"],6),k5=round(s200["k5"],6),tab_ratio=round(s200["tab_ratio"],6),max_row_bf=round(max_row_bf,2),max_row_bg=round(max_row_bg,2),Rp_raw=Rp_raw,n_surv=n_surv,stageb_row_keep_frac=round(stageb_row_keep_frac,6),ret_f_realized=round(ret_f_realized,6),ret_g_realized=round(ret_g_realized,6),post_peel_cols_lost=post_peel_cols_lost,post_merge_cols_lost=post_merge_cols_lost)
def verify_h3_ladder_and_fit(la,log2_Up,run_u12=False,op=None,nf_cond=None,ng_cond=None):
 lad=la["ladder"]
 assert[x["log2U"]for x in lad]==list(range(16,27)),"Ladder must span 2^16..2^26!"
 for i in range(1,len(lad)):
  a,b=lad[i-1],lad[i]
  assert b["nl"]<a["nl"]and b["sw"]<a["sw"]and b["exp_unit"]>a["exp_unit"],"Active LA ladder metrics (nl, sw, exp_unit) must be strictly monotonic across 2^16..2^26!"
 assert all(x["tab_ratio"]<=0.83 for x in lad),"tab_ratio must stay <= 0.83 << 1.0 = tab_ratio_charged across 2^16..2^26!"
 r26=lad[-1]
 assert all(x["post_peel_cols_lost_max"]==0 and x["post_merge_cols_lost_max"]==0 for x in lad),"Zero column loss required across 2^16..2^26!"
 # Live binding of archived ladder parameters & row-width invariants to operating point (H1/H2) + Stage-B caps & retained merge_stack_mean_w:
 assert all(x["max_row_w"]<=4096.0 for x in lad),"All ladder runs must satisfy max_row_w <= W_max = 4,096!"
 assert all(x["avg_w"]<201.0 for x in lad if x["log2U"]>=19),"All U >= 2^19 ladder runs must satisfy avg_w < 201.0!"
 assert all(x.get("merge_stack_mean_w_max",0.0)<=201.0 for x in lad),"Aggregate retained merge_stack mean weight must be <= 201.0!"
 assert all(abs(x.get("sbf",355.0)-355.0)<1e-9 and abs(x.get("sbg",294.0)-294.0)<1e-9 and abs(x.get("ret_f",1.0)-1.0)<1e-9 and abs(x.get("ret_g",1.0)-1.0)<1e-9 and x.get("stageb_mode")=="row_rejection" and 0.50<x.get("stageb_row_keep_frac",0.0)<0.85 for x in lad),"Stage-B row-rejection caps (355, 294) & unconditioned raw rows (ret=1.0) must be bound across ladder!"
 assert all(x.get("max_row_bf_max",0.0)<=355.0+1e-9 and x.get("max_row_bg_max",0.0)<=294.0+1e-9 for x in lad),"Stage-B row residual caps violated!"
 if op is not None:
  assert all(abs(x["excess"]-float(op["exc"]))<1e-9 for x in lad),"Ladder excess must match op['exc']!"
 if nf_cond is not None and ng_cond is not None:
  assert all(abs(x["nf_sim"]-nf_cond)<=2.0 and abs(x["ng_sim"]-ng_cond)<=2.0 for x in lad),"Ladder conditioned norms must match operating point (two-sided +-2.0b)!"
 assert abs(la["eu_26"]-r26["exp_unit"])<1e-9 and abs(la["se_eu_26"]-r26["se_eu"])<1e-9,"cfg_la eu_26/se_eu_26 must match 2^26 ladder!"
 assert la["nl_cap"]>1.45*r26["nl"]and la["sw_cap"]>1.45*r26["sw"],"H3 caps must have >1.45x headroom over 2^26 means!"
 assert la["nl_cap"]>(r26["nl"]+2.0*r26["se_l"])*1.45 and la["sw_cap"]>(r26["sw"]+2.0*r26["se_sw"])*1.45
 assert all(la["nl_cap"]>x["nl"]and la["sw_cap"]>x["sw"]for x in lad[7:]),"H3 caps must exceed EVERY 2^23..2^26 mean!"
 assert r26["unit_frac"]>=0.970 and r26["exp_unit"]>=0.900 and la["tab_ratio_charged"]==1.0
 eu_lo=la["eu_26"]-2.0*la["se_eu_26"]
 _nonunit_ir=float(la.get("gen_ir",22.0))if la.get("k5_mode")=="hex6" else 116.0
 _unit_ir=float(la.get("unit_ir",7.65))
 k5_live=eu_lo*_unit_ir+(1.0-eu_lo)*_nonunit_ir
 assert la["k5_charged"]>=k5_live+0.60 and all(la["k5_charged"]>(x["exp_unit"]*_unit_ir+(1.0-x["exp_unit"])*_nonunit_ir)for x in lad[2:])
 xs=[math.log(x["log2U"])for x in lad];mx=sum(xs)/len(xs);sxx=sum((x-mx)**2 for x in xs)
 fits={}
 for k,cap in(("nl",la["nl_cap"]),("sw",la["sw_cap"])):
  ys=[math.log(x[k])for x in lad];my=sum(ys)/len(ys);sl=sum((x-mx)*(y-my)for x,y in zip(xs,ys))/sxx
  al=-sl;c=math.exp(my-sl*mx);merr=max(abs((c/(x["log2U"]**al))-x[k])/x[k]for x in lad);pr=c/(log2_Up**al)
  assert al>1.50 and merr<0.25 and cap/pr>5.0
  fits[k]=dict(c=round(c,6),alpha=round(al,6),max_rel_err=round(merr,4),proj=round(pr,7),cap_over_proj=round(cap/pr,2))
 u12=None;u12_lla=None
 paired=la.get("paired_base_vs_lla")
 if paired is not None:
  assert len(paired)>=4 and all(p["lla_nl"]<p["base_nl"]and p["lla_sw"]<p["base_sw"]and p["lla_lost_max"]==0 and p["base_lost_max"]==0 for p in paired),"Paired L-LA trimming must strictly reduce nl and sw with zero column loss!"
 if run_u12:
  u12=run_filtsim_z12_ladder(12,0)
  # Check true Stage-B row-rejection caps (345, 290), retained merge_stack_mean_w <= 201, BOTH avg_w bounds, and bit-exact C/Python parity!
  assert u12["N_sw"]==16 and u12["N_l"]==13 and u12["N_sw_250"]==15 and u12["depth"]==2 and u12["avg_w"]<=200.0+64.0*4096.0/u12["N_sw"] and u12["avg_w"]<=200.0+2.0*u12["max_row_w"]/u12["N_sw"] and u12["max_row_w"]<=4096 and u12["merge_stack_mean_w"]<=201.0 and abs(u12["exp_unit"]-0.588101)<1e-6 and u12["Rp_raw"]==7671 and u12["n_surv"]==5530 and abs(u12["stageb_row_keep_frac"]-0.720897)<1e-6 and u12["max_row_bf"]<=345.0+1e-9 and u12["max_row_bg"]<=290.0+1e-9 and u12["post_peel_cols_lost"]==0 and u12["post_merge_cols_lost"]==0 and u12["small4_nu_frac"]>0.25
  u12_lla=run_filtsim_z12_ladder(12,0,Lf=72.70,Lg=73.06,Nf=539.2,Ng=663.2,excess=0.62,pf_s=60.4,pg_s=63.6,sbf=355.0,sbg=294.0,pexc=0.0)
  assert u12_lla["N_sw"]==26 and u12_lla["N_l"]==19 and u12_lla["N_sw_250"]==24 and u12_lla["depth"]==5 and u12_lla["max_row_w"]<=4096 and u12_lla["merge_stack_mean_w"]<=201.0 and abs(u12_lla["exp_unit"]-0.700542)<1e-6 and u12_lla["max_row_bf"]<=355.0+1e-9 and u12_lla["max_row_bg"]<=294.0+1e-9 and u12_lla["post_peel_cols_lost"]==0 and u12_lla["post_merge_cols_lost"]==0
 return dict(scales="16..26",n_seeds_26=r26["n_seeds"],stage_b_capped=True,stage_b_mode="row_rejection",stageb_row_keep_frac_26=r26["stageb_row_keep_frac"],ret_f_realized_26=r26["ret_f_realized"],ret_g_realized_26=r26["ret_g_realized"],all_active_monotonic=True,max_row_w_all_le_4096=True,avg_w_ge19_lt_201=True,merge_stack_mean_w_max=round(max(x.get("merge_stack_mean_w_max",0.0)for x in lad),4),nl_26=r26["nl"],se_l_26=r26["se_l"],sw_26=r26["sw"],se_sw_26=r26["se_sw"],eu_26=r26["exp_unit"],se_eu_26=r26["se_eu"],nl_headroom_26=round(la["nl_cap"]/r26["nl"],4),sw_headroom_26=round(la["sw_cap"]/r26["sw"],4),tab_ratio_charged=1.0,unit_ir=_unit_ir,gen_ir=_nonunit_ir,k5_live=round(k5_live,4),k5_charged=la["k5_charged"],inv_log_fit=fits,paired_base_vs_lla=paired,live_u12_seed0=u12,live_u12_lla_seed0=u12_lla)
if __name__=="__main__":main()
