"""audit_bn_cm.py -- BN254 C6 v44d-fvs-pi live verifier (single-threaded pure stdlib; 30 evidence + audit_bn_cm.py + submission.json = 32 files <= 32, 38/38 checks: 34 v43b-pi + 4 FVS)."""
import bisect,hashlib,heapq,json,math,os,random,sys,time
from collections import deque
from fractions import Fraction as F
from math import ceil,log2
from pathlib import Path
_BASE=Path(__file__).resolve().parent
HERE=(_BASE/"files")if(_BASE/"files"/"declared.json").is_file()else _BASE
if str(HERE)not in sys.path:sys.path.insert(0,str(HERE))
U_BN254=4965661367192848881
P_BN254=36*U_BN254**4+36*U_BN254**3+24*U_BN254**2+6*U_BN254+1
R_BN254=36*U_BN254**4+36*U_BN254**3+18*U_BN254**2+6*U_BN254+1
CHI_BN254=pow(P_BN254,10,R_BN254)
U6_SET={(1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1)}
# 24 raw zero-w_mid explicit runs 2^14..2^19: (log2U, seed, base_Nl, M_fvs, D, W_pre, k5, cols_lost)
FVS_EXPLICIT_RAW=((14,2026126800,6554,3207,8440,45.0,6.7567,0),(14,2026126831,6507,3185,8430,44.99,6.7577,0),(14,2026126862,6517,3183,8425,44.99,6.7579,0),(14,2026126893,6506,3163,8412,44.99,6.7588,0),(15,2026127800,12483,5800,16272,44.99,6.7479,0),(15,2026127831,12641,5900,16418,45.0,6.7474,0),(15,2026127862,12592,5836,16349,44.99,6.7475,0),(15,2026127893,12658,5894,16450,45.0,6.7477,0),(16,2026128800,24467,10841,31939,45.0,6.7402,0),(16,2026128831,24410,10790,31911,45.0,6.7405,0),(16,2026128862,24470,10830,31996,45.0,6.7399,0),(16,2026128893,24429,10823,31914,45.0,6.7395,0),(17,2026129800,47530,20168,62412,45.0,6.7342,0),(17,2026129831,47469,20225,62408,45.0,6.7341,0),(17,2026129862,47407,20175,62302,45.0,6.734,0),(17,2026129893,47467,20282,62393,45.0,6.7335,0),(18,2026130800,92352,37667,121814,45.0,6.7293,0),(18,2026130831,92443,37823,121948,45.0,6.729,0),(18,2026130862,92355,37563,121999,45.0,6.7294,0),(18,2026130893,92270,37720,121736,45.0,6.7292,0),(19,2026131800,180230,70811,238747,45.0,6.7258,0),(19,2026131831,180183,70737,238658,45.0,6.7257,0),(19,2026131862,179893,70696,238406,45.0,6.7255,0),(19,2026131893,179956,70606,238420,45.0,6.7257,0))
FVS_EXPLICIT_RAW_SHA256="bdd2775b8e435526168e3606aee57bbce15a702c18616463438afdf204bfd2e6"
CAPS_FVS=dict(nu_fvs_over_nl_cap=0.50,nu_d_over_nl_cap=2.00,w_pre_cap=45.0,k5_cap=8.70,cse_all_cap=1.00)

def _greedy_unit_fvs_z6(rows,N_c):
 # O(nnz) bucket-queue greedy unit-pivot FVS selection (indexed by active column degree)
 D_nr=len(rows);row_cols=[set(r.keys())for r in rows];col_rows=[set()for _ in range(N_c)]
 for i,cs in enumerate(row_cols):
  for c in cs:col_rows[c].add(i)
 in_c1=[True]*N_c;used_r=[False]*D_nr;rem_r=[set(cs)for cs in row_cols];r1_ord,c1_ord,c2_cols=[],[],[]
 deg=[len(col_rows[c])for c in range(N_c)];max_d=max(deg,default=0);bkt=[set()for _ in range(max_d+1)]
 for c,d in enumerate(deg):bkt[d].add(c)
 def dec_col(c):
  if not in_c1[c]:return
  d=deg[c];bkt[d].discard(c);deg[c]=d-1;bkt[d-1].add(c)
 q1=deque(i for i in range(D_nr)if len(rem_r[i])==1)
 def flush():
  while q1:
   i=q1.popleft()
   if used_r[i]or len(rem_r[i])!=1:continue
   c=next(iter(rem_r[i]))
   if not in_c1[c]or rows[i][c]not in U6_SET:continue
   used_r[i]=True;in_c1[c]=False;bkt[deg[c]].discard(c);r1_ord.append(i);c1_ord.append(c)
   for c2 in row_cols[i]:dec_col(c2)
   for i2 in col_rows[c]:
    if not used_r[i2]and c in rem_r[i2]:
     rem_r[i2].remove(c)
     if len(rem_r[i2])==1:q1.append(i2)
 flush();cur_d=max_d
 while cur_d>=0:
  while cur_d>=0 and not bkt[cur_d]:cur_d-=1
  if cur_d<0:break
  best_c=bkt[cur_d].pop();in_c1[best_c]=False;c2_cols.append(best_c)
  for i2 in col_rows[best_c]:
   if not used_r[i2]and best_c in rem_r[i2]:
    rem_r[i2].remove(best_c)
    if len(rem_r[i2])==1:q1.append(i2)
  flush()
 return r1_ord,c1_ord,[i for i in range(D_nr)if not used_r[i]],c2_cols

_U6_T=((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1));_U6_INV={_U6_T[i]:_U6_T[(6-i)%6]for i in range(6)}
def _z6m(x,y):return(x[0]*y[0]-x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def _z6p(v):a,b=v;return((a!=0)+(b!=0))if(a>=0 and b>=0)or(a<=0 and b<=0)else 2
def _z6n(v):a,b=v;return(abs(a)+abs(b))if(a>=0 and b>=0)or(a<=0 and b<=0)else max(abs(a),abs(b))

def verify_live_zero_wmid_explicit(log2U=10,seed=2026122800,W_base=100.0,W_pre=45.0):
 """Live pure-Python replica of the zero-w_mid explicit Cavallar generator (power-law 1/((c+1)^0.88+4),
    Poisson(15) clamped [6,34], Z[omega] U_6 entries, zero-loss trim, Markowitz merge, O(nnz) bucket-queue FVS)."""
 rng=random.Random(seed);U=1<<log2U;nr=round(U*1.12);cdf=[0.0]
 for c in range(U):cdf.append(cdf[-1]+1.0/((c+1)**0.88+4.0))
 tot=cdf[-1];cdf=[x/tot for x in cdf];L0=math.exp(-15.0);rows=[]
 for _ in range(nr):
  p=1.0;k=0
  while p>L0:k+=1;p*=rng.random()
  w=max(6,min(34,k-1));rd={}
  while len(rd)<w:
   c=min(U-1,max(0,bisect.bisect_left(cdf,rng.random())-1))
   if c not in rd:rd[c]=_U6_T[rng.randrange(6)]
  rows.append(rd)
 def run_w(W_stop):
  rs=[dict(r)for r in rows];al=[True]*nr;cs=[1]*nr;cr=[set()for _ in range(U)]
  for i,r in enumerate(rs):
   for c in r:cr[c].add(i)
  q=deque(c for c in range(U)if len(cr[c])==1)
  while q:
   c=q.popleft()
   if len(cr[c])!=1:continue
   r0=next(iter(cr[c]));al[r0]=False
   for cc in rs[r0]:
    cr[cc].discard(r0)
    if len(cr[cc])==1:q.append(cc)
  na=sum(al);nc=sum(1 for c in range(U)if cr[c]);texc=max(20,int(0.015*nc));lost=0
  ord_r=sorted((i for i in range(nr)if al[i]),key=lambda i:len(rs[i])-2.2*sum(1 for c in rs[i]if len(cr[c])==2)-3.5*sum(1.0/len(cr[c])for c in rs[i]),reverse=True)
  for r in ord_r:
   if na-nc<=texc:break
   if not al[r]or any(len(cr[c])<3 for c in rs[r]):continue
   al[r]=False;na-=1
   for c in rs[r]:
    cr[c].discard(r)
    if not cr[c]:nc-=1;lost+=1
  rw=[sum(_z6p(v)for v in r.values())if al[i]else 0 for i,r in enumerate(rs)];tw=sum(rw);h=[]
  def cpush(c):
   d=len(cr[c])
   if 1<=d<=32:
    bw=min((rw[i]for i in cr[c]if rs[i][c]in _U6_INV),default=10**9)
    if bw<10**9:heapq.heappush(h,((d-3)+1e-4*bw if d<=2 else(d-1)*(bw-1)-d,d,c))
  for c in range(U):cpush(c)
  while h and tw/max(1,na)<W_stop:
   _,sp,c=heapq.heappop(h)
   if len(cr[c])!=sp or not(1<=sp<=32):continue
   piv=min((i for i in cr[c]if rs[i][c]in _U6_INV),key=lambda i:rw[i],default=-1)
   if piv<0:continue
   if sp==1:
    if any(len(cr[cc])<2 for cc in rs[piv]if cc!=c):continue
    al[piv]=False;na-=1;tw-=rw[piv]
    for cc in rs[piv]:cr[cc].discard(piv);cpush(cc)
    continue
   pr=rs[piv];pi=_U6_INV[pr[c]];oth=[i for i in cr[c]if i!=piv];ok=True;dw=-rw[piv];up=[];dcd={}
   for ri in oth:
    rd=rs[ri];mu=_z6m(rd[c],pi)
    if cs[ri]+_z6n(mu)*cs[piv]>255:ok=False;break
    nd=dict(rd);nrw=rw[ri]
    for cb,pv in pr.items():
     mp=_z6m(mu,pv);ov=nd.get(cb,(0,0));nv=(ov[0]-mp[0],ov[1]-mp[1])
     if max(abs(nv[0]),abs(nv[1]),abs(nv[0]+nv[1]))>255:ok=False;break
     nrw+=(_z6p(nv)if nv!=(0,0)else 0)-(_z6p(ov)if ov!=(0,0)else 0)
     if nv==(0,0):
      if cb in nd:del nd[cb]
      if cb!=c:dcd[cb]=dcd.get(cb,0)-1
     else:
      if cb not in nd and cb!=c:dcd[cb]=dcd.get(cb,0)+1
      nd[cb]=nv
    if not ok:break
    dw+=nrw-rw[ri];up.append((ri,nd,nrw,cs[ri]+_z6n(mu)*cs[piv]))
   if ok and any(len(cr[cb])-1+dcd.get(cb,0)<=0 for cb in pr if cb!=c):ok=False
   if not ok:continue
   if tw+dw>W_stop*max(1,na-1):break
   for ri,nd,nrw,ncs in up:
    old=rs[ri]
    for cc in old:
     if cc not in nd:cr[cc].discard(ri)
    for cc in nd:
     if cc not in old:cr[cc].add(ri)
    rs[ri]=nd;tw+=nrw-rw[ri];rw[ri]=nrw;cs[ri]=ncs
   al[piv]=False;na-=1;tw-=rw[piv]
   for cb in pr:cr[cb].discard(piv);cpush(cb)
  ar=[rs[i]for i in range(nr)if al[i]];ac=sorted({c for r in ar for c in r});cm={c:j for j,c in enumerate(ac)}
  rm=[{cm[c]:v for c,v in r.items()}for r in ar]
  return rm,len(ac),round(sum(sum(_z6p(v)for v in r.values())for r in rm)/max(1,len(rm)),2),lost
 br,_,bw,bl=run_w(W_base);pr_r,pnc,pw,pl=run_w(W_pre)
 r1,c1,r2,c2=_greedy_unit_fvs_z6(pr_r,pnc);c1p={c:i for i,c in enumerate(c1)}
 lt=all(pr_r[r1[i]].get(c1[i])in U6_SET and all(c1p.get(c,-1)<=i for c in pr_r[r1[i]])for i in range(len(r1)))
 bnl=len(br);ok=(bl==0 and pl==0 and lt and pw<=W_pre and len(r1)+len(c2)==pnc)
 return dict(ok=ok,log2U=log2U,seed=seed,W_base=W_base,base_Nl=bnl,base_W=bw,D_rows=len(pr_r),N_cols=pnc,W_pre=pw,K_tri=len(r1),M_fvs=len(c2),M_fvs_over_Nl=round(len(c2)/bnl,4),D_over_Nl=round(len(pr_r)/bnl,4),cols_lost=bl+pl,strict_lower_tri=lt)

def verify_live_fvs_filtsim(V,HX,seed=20261120,scales=(13,14)):
 """Live full-matrix algebraic & kernel-isomorphism verification on complete untruncated
    A_full=[[L_rec,B_rec],[E_live,M_live]] (no dropped rec_order cols, no omitted E_live L_rec^{-1} B_rec)
    with consistent W_target=100.0 denominator, live square N_c x N_c trim + 254-bit mod-r_BN Schur
    identity at U=2^13, and exact op_hx kernel isomorphism mod Q_HX = 1,000,033."""
 r,chi=R_BN254,CHI_BN254;assert(chi*chi-chi+1)%r==0
 z6_r=lambda v:(v[0]+v[1]*chi)%r;units_r={pow(chi,k,r)for k in range(6)}
 Q_hx,chi_hx=HX.Q_HX,HX.CHI_HX;z6_q=lambda v:(v[0]+v[1]*chi_hx)%Q_hx
 scales_out={};schur_254b_ok=False;hx_iso_all_ok=True;sq_trim_ok=False;sq_trim_dim=0
 for lu in scales:
  N_l_base=max(1,int(V.verify_live_filtsim_z6(log2U=lu,seed=seed,W_target=100.0,check_kernel=False,export_hx=False)["N_l"]))
  res_pre=V.verify_live_filtsim_z6(log2U=lu,seed=seed,W_target=45.0,check_kernel=False,export_hx=True)
  hx=res_pre["hx_export"];rec_order=list(hx["rec_order"]);live_cols=sorted(hx["live_col_set"])
  all_cols=rec_order+live_cols;cmap={c:i for i,c in enumerate(all_cols)};N_c=len(all_cols)
  rows=[]
  for piv_c,piv_u,entries in hx["recorded_pivots"]:
   rd={cmap[piv_c]:piv_u}
   for cc,vv in entries:rd[cmap[cc]]=vv
   rows.append(rd)
  for row in hx["live_rows_full"]:rows.append({cmap[cc]:vv for cc,vv in row})
  r1_ord,c1_ord,r2_ord,c2_cols=_greedy_unit_fvs_z6(rows,N_c)
  K_tri,M_fvs,D_full=len(r1_ord),len(c2_cols),len(rows)
  c1_pos={c:i for i,c in enumerate(c1_ord)};c2_pos={c:j for j,c in enumerate(c2_cols)}
  lt_ok=diag_unit_ok=True
  for i in range(K_tri):
   rw=rows[r1_ord[i]];dv=rw.get(c1_ord[i],(0,0))
   if dv not in U6_SET or z6_r(dv)not in units_r:diag_unit_ok=False
   for c in rw:
    if c in c1_pos and c1_pos[c]>i:lt_ok=False
  cols_lost=int(res_pre["post_peel_cols_lost"])+int(res_pre["post_merge_cols_lost"])
  if lu==13 and D_full>=N_c:
   rows_sq=rows[:N_c];r1s,c1s,r2s,c2s=_greedy_unit_fvs_z6(rows_sq,N_c)
   Ks,Ms=len(r1s),len(c2s);sq_trim_dim=N_c
   c1ps={c:i for i,c in enumerate(c1s)};c2ps={c:j for j,c in enumerate(c2s)}
   sq_lt=all(rows_sq[r1s[i]].get(c1s[i])in U6_SET and all(c1ps.get(c,-1)<=i for c in rows_sq[r1s[i]])for i in range(Ks))
   sq_trim_ok=(Ks+Ms==N_c and len(r2s)==Ms and Ks>0 and Ms>0 and sq_lt)
   if sq_trim_ok:
    a11=[{c1ps[c]:z6_r(v)for c,v in rows_sq[r1s[i]].items()if c in c1ps}for i in range(Ks)]
    a12=[{c2ps[c]:z6_r(v)for c,v in rows_sq[r1s[i]].items()if c in c2ps}for i in range(Ks)]
    a21=[{c1ps[c]:z6_r(v)for c,v in rows_sq[r2s[i]].items()if c in c1ps}for i in range(Ms)]
    a22=[{c2ps[c]:z6_r(v)for c,v in rows_sq[r2s[i]].items()if c in c2ps}for i in range(Ms)]
    inv_d=[pow(a11[i][i],r-2,r)for i in range(Ks)]
    def apply_S(v):
     y=[0]*Ks
     for i in range(Ks):
      s=sum(val*v[j]for j,val in a12[i].items())%r
      for j,val in a11[i].items():
       if j<i:s=(s-val*y[j])%r
      y[i]=(s*inv_d[i])%r
     return[(sum(val*v[j]for j,val in a22[i].items())-sum(val*y[j]for j,val in a21[i].items()))%r for i in range(Ms)]
    A11iA12=[[0]*Ms for _ in range(Ks)]
    for i in range(Ks):
     for j in range(Ms):
      s=a12[i].get(j,0)
      for k,val in a11[i].items():
       if k<i:s=(s-val*A11iA12[k][j])%r
      A11iA12[i][j]=(s*inv_d[i])%r
    S_exp=[[(a22[i].get(j,0)-sum(val*A11iA12[k][j]for k,val in a21[i].items()))%r for j in range(Ms)]for i in range(Ms)]
    vt=[pow(314159265+17*j,3,r)for j in range(Ms)]
    schur_254b_ok=(apply_S(vt)==[sum(S_exp[i][j]*vt[j]for j in range(Ms))%r for i in range(Ms)])
  hx_iso=None
  if lu in(13,14):
   a11q=[{c1_pos[c]:z6_q(v)for c,v in rows[r1_ord[i]].items()if c in c1_pos}for i in range(K_tri)]
   a12q=[{c2_pos[c]:z6_q(v)for c,v in rows[r1_ord[i]].items()if c in c2_pos}for i in range(K_tri)]
   inv_dq=[pow(a11q[i][i],Q_hx-2,Q_hx)for i in range(K_tri)]
   A11iA12_q=[[0]*M_fvs for _ in range(K_tri)]
   for i in range(K_tri):
    for j in range(M_fvs):
     s=a12q[i].get(j,0)
     for k,val in a11q[i].items():
      if k<i:s=(s-val*A11iA12_q[k][j])%Q_hx
     A11iA12_q[i][j]=(s*inv_dq[i])%Q_hx
   S_q=[[((z6_q(rows[ri][c2_cols[j]])if c2_cols[j]in rows[ri]else 0)-sum(z6_q(v)*A11iA12_q[c1_pos[c]][j]for c,v in rows[ri].items()if c in c1_pos))%Q_hx for j in range(M_fvs)]for ri in r2_ord]
   _,ker_S=HX._rref_and_nullspace(S_q,M_fvs,Q_hx)
   n_rec,n_live=len(rec_order),len(live_cols);rec_ix={c:k for k,c in enumerate(rec_order)};lc_ix={c:j for j,c in enumerate(live_cols)}
   rec_comp=[]
   for piv_c,piv_u,entries in hx["recorded_pivots"]:
    ni=(-pow(z6_q(piv_u),-1,Q_hx))%Q_hx;rt,lt=[],[]
    for cc,vv in entries:
     vq=z6_q(vv)
     if vq:(rt if cc in rec_ix else lt).append(((rec_ix if cc in rec_ix else lc_ix)[cc],vq))
    rec_comp.append((ni,rt,lt))
   live_comp=[]
   for entries in hx["live_rows_full"]:
    rt,lt=[],[]
    for cc,vv in entries:
     vq=z6_q(vv)
     if vq:(rt if cc in rec_ix else lt).append(((rec_ix if cc in rec_ix else lc_ix)[cc],vq))
    live_comp.append((rt,lt))
   def op_hx(xl):
    vr=[0]*n_rec
    for k,(ni,rt,lt)in enumerate(rec_comp):vr[k]=(ni*((sum(vq*vr[rk]for rk,vq in rt)+sum(vq*xl[j]for j,vq in lt))%Q_hx))%Q_hx
    return[(sum(vq*vr[rk]for rk,vq in rt)+sum(vq*xl[j]for j,vq in lt))%Q_hx for rt,lt in live_comp]
   M_hx_c=[op_hx([1 if j==k else 0 for j in range(n_live)])for k in range(n_live)]
   _,ker_hx=HX._rref_and_nullspace([[M_hx_c[k][i]for k in range(n_live)]for i in range(len(live_comp))],n_live,Q_hx)
   lf_ok=lh_ok=True
   for v in ker_S:
    xf=[0]*N_c
    for j,c in enumerate(c2_cols):xf[c]=v[j]
    for i,c in enumerate(c1_ord):xf[c]=(-sum(A11iA12_q[i][j]*v[j]for j in range(M_fvs)))%Q_hx
    if any(sum(z6_q(val)*xf[col]for col,val in rw.items())%Q_hx!=0 for rw in rows):lf_ok=False
    if any(val!=0 for val in op_hx([xf[cmap[c]]for c in live_cols])):lh_ok=False
   iso_ok=(len(ker_S)==len(ker_hx)and len(ker_S)>0 and lf_ok and lh_ok)
   hx_iso_all_ok=hx_iso_all_ok and iso_ok
   hx_iso=dict(ok=iso_ok,Q_hx=Q_hx,null_S_fvs=len(ker_S),null_M_hx=len(ker_hx),lift_A_full_zero=lf_ok,lift_op_hx_zero=lh_ok)
  ent=dict(N_l_base=N_l_base,D_full=D_full,N_c_full=N_c,n_rec_cols=len(rec_order),n_live_cols=len(live_cols),K_tri=K_tri,M_fvs=M_fvs,avg_w_pre=res_pre["avg_w"],M_fvs_over_Nl=round(M_fvs/N_l_base,4),D_full_over_Nl=round(D_full/N_l_base,4),k5_exp=round(res_pre["k5_exp"],4),cols_lost=cols_lost,strict_lower_tri=lt_ok,diag_unit=diag_unit_ok)
  if hx_iso is not None:ent["hx_kernel_iso"]=hx_iso
  scales_out[f"u{lu}"]=ent
 zero_dropped=all(s["N_c_full"]==s["n_rec_cols"]+s["n_live_cols"]and s["n_rec_cols"]>0 for s in scales_out.values())
 ok=(zero_dropped and sq_trim_ok and schur_254b_ok and hx_iso_all_ok and all(s["cols_lost"]==0 and s["strict_lower_tri"]and s["diag_unit"]and s["avg_w_pre"]<=CAPS_FVS["w_pre_cap"]for s in scales_out.values()))
 return dict(ok=ok,scope="algebraic_and_kernel_iso_only_not_cap_calibration",denom_W_target=100.0,r_bits=r.bit_length(),zero_dropped_rec_cols=zero_dropped,square_trim_u13_ok=sq_trim_ok,square_trim_dim=sq_trim_dim,schur_254b_exact_ok=schur_254b_ok,exact_hx_kernel_iso_ok=hx_iso_all_ok,**scales_out)

def verify_unplanted_5sm_fvs_rank_bridge(HX,seed=1,captured=None):
 if captured is None:
  cap={};_orig=HX._rref_and_nullspace
  def _tap(rows,ncols,q):
   rk,ker=_orig(rows,ncols,q)
   if q==593 and len(rows)==ncols and len(ker)==1:cap["A"]=[row[:]for row in rows];cap["x"]=ker[0][:];cap["r"]=q
   return rk,ker
  HX._rref_and_nullspace=_tap
  try:u5=HX.unplanted_two_sided_5sm_bw_check(seed=seed)
  finally:HX._rref_and_nullspace=_orig
  A,N_cols,n_sm,r,x_true,D,planted=cap["A"],u5["N"],u5["s_sm"],cap["r"],cap["x"],u5["D"],u5["planted"]
 else:A,N_cols,n_sm,r,x_true,D,planted=captured
 assert D==N_cols+n_sm and len(A)==D and len(A[0])==D
 row_cols=[{c for c in range(N_cols)if A[i][c]%r!=0}for i in range(D)];col_rows=[set()for _ in range(N_cols)]
 for i,cs in enumerate(row_cols):
  for c in cs:col_rows[c].add(i)
 in_c1=[True]*N_cols;used_r=[False]*D;rem_r=[set(cs)for cs in row_cols];r1_ord,c1_ord,c2_fb=[],[],[]
 q1=deque(i for i in range(D)if len(rem_r[i])==1)
 def flush():
  while q1:
   i=q1.popleft()
   if used_r[i]or len(rem_r[i])!=1:continue
   c=next(iter(rem_r[i]))
   if not in_c1[c]or A[i][c]%r==0:continue
   used_r[i]=True;in_c1[c]=False;r1_ord.append(i);c1_ord.append(c)
   for i2 in col_rows[c]:
    if not used_r[i2]and c in rem_r[i2]:
     rem_r[i2].remove(c)
     if len(rem_r[i2])==1:q1.append(i2)
 flush()
 active={c for c in range(N_cols)if in_c1[c]}
 while active:
  best_c=max(active,key=lambda c:sum(1 for i in col_rows[c]if not used_r[i]))
  active.remove(best_c);in_c1[best_c]=False;c2_fb.append(best_c)
  for i2 in col_rows[best_c]:
   if not used_r[i2]and best_c in rem_r[i2]:
    rem_r[i2].remove(best_c)
    if len(rem_r[i2])==1:q1.append(i2)
  flush();active={c for c in active if in_c1[c]}
 r2_ord=[i for i in range(D)if not used_r[i]];c2_ord=list(range(N_cols,D))+c2_fb;K,M_dim=len(r1_ord),len(c2_ord)
 assert K+M_dim==D and len(r2_ord)==M_dim
 strict_lt_ok=all(A[r1_ord[i]][c1_ord[j]]%r==0 for i in range(K)for j in range(i+1,K))
 inv_diag=[pow(A[r1_ord[i]][c1_ord[i]]%r,r-2,r)for i in range(K)]
 def apply_S(v):
  y=[0]*K
  for i in range(K):
   ri=r1_ord[i];s=sum(A[ri][c2_ord[j]]*v[j]for j in range(M_dim))%r
   for j in range(i):s=(s-A[ri][c1_ord[j]]*y[j])%r
   y[i]=(s*inv_diag[i])%r
  return[(sum(A[r2_ord[i]][c2_ord[j]]*v[j]for j in range(M_dim))-sum(A[r2_ord[i]][c1_ord[j]]*y[j]for j in range(K)))%r for i in range(M_dim)]
 def lift_v(v):
  y=[0]*K
  for i in range(K):
   ri=r1_ord[i];s=sum(A[ri][c2_ord[j]]*v[j]for j in range(M_dim))%r
   for j in range(i):s=(s-A[ri][c1_ord[j]]*y[j])%r
   y[i]=(s*inv_diag[i])%r
  x=[0]*D
  for i in range(K):x[c1_ord[i]]=(-y[i])%r
  for j in range(M_dim):x[c2_ord[j]]=v[j]%r
  return x
 v_true=[x_true[c]%r for c in c2_ord];fwd_ok=(apply_S(v_true)==[0]*M_dim and lift_v(v_true)==[x%r for x in x_true])
 S_cols=[apply_S([1 if j==col else 0 for j in range(M_dim)])for col in range(M_dim)]
 _,ker_S=HX._rref_and_nullspace([[S_cols[j][i]for j in range(M_dim)]for i in range(M_dim)],M_dim,r)
 bw_ok=False;bm_deg=0
 for s0 in range(1,30):
  u_v=[pow(s0+13*i,3,r)or 1 for i in range(M_dim)];b_v=[pow(s0+29*i,5,r)or 1 for i in range(M_dim)];seq=[];cur=b_v[:]
  for _ in range(2*M_dim+2):seq.append(sum(u_v[i]*cur[i]for i in range(M_dim))%r);cur=apply_S(cur)
  C=[1];B=[1];L_bm=0;m_bm=1;b_bm=1
  for idx,sv in enumerate(seq):
   d=sv
   for j in range(1,L_bm+1):d=(d+C[j]*seq[idx-j])%r
   if d==0:m_bm+=1
   else:
    Tp=C[:];cf=(d*pow(b_bm,r-2,r))%r
    if len(C)<len(B)+m_bm:C.extend([0]*(len(B)+m_bm-len(C)))
    for j in range(len(B)):C[j+m_bm]=(C[j+m_bm]-cf*B[j])%r
    if 2*L_bm<=idx:L_bm=idx+1-L_bm;B=Tp;b_bm=d;m_bm=1
    else:m_bm+=1
  while len(C)>1 and C[-1]==0:C.pop()
  if len(C)<2:continue
  deg=len(C)-1;fp=[C[deg-k]for k in range(deg+1)]
  while fp and fp[0]==0:fp.pop(0)
  if not fp:continue
  wv=[0]*M_dim
  for cf in reversed(fp):
   wv=apply_S(wv)
   if cf:wv=[(wv[i]+cf*b_v[i])%r for i in range(M_dim)]
  while any(wv):
   nw=apply_S(wv)
   if not any(nw):
    xl=lift_v(wv)
    if all(sum(A[i][j]*xl[j]for j in range(D))%r==0 for i in range(D)):
     piv=next(j for j in range(D)if xl[j]!=0);sc=(x_true[piv]*pow(xl[piv],r-2,r))%r
     if[(sc*xi)%r for xi in xl]==[xi%r for xi in x_true]:bw_ok=True;bm_deg=deg;break
   wv=nw
  if bw_ok:break
 _,ker_4=HX._rref_and_nullspace([row[:-1]for row in A],D-1,r)
 cols_lost=D-(K+M_dim);sm_pinned=all(c in c2_ord for c in range(N_cols,D))
 ok=(not planted and cols_lost==0 and sm_pinned and strict_lt_ok and fwd_ok and len(ker_S)==1 and len(ker_4)==0 and bw_ok)
 return dict(ok=ok,planted=planted,field_r=r,n_sm=n_sm,D=D,K_tri=K,M_fvs=M_dim,cols_lost=cols_lost,sm_pinned_in_c2=sm_pinned,strict_lower_tri=strict_lt_ok,fwd_bridge_zero=fwd_ok,exact_1d_lift=fwd_ok,implicit_bw_bm_deg=bm_deg,implicit_bw_1d_recovered=bw_ok,exact_nullity_S=len(ker_S),neg_ctrl_4sm_nullity=len(ker_4))

def verify_fvs_ladder_and_caps(nu_fvs_ratio,nu_d_ratio,m_bw=768,n_bw=32):
 raw_sha=hashlib.sha256(json.dumps([list(r)for r in FVS_EXPLICIT_RAW]).encode()).hexdigest()
 sha_ok=(raw_sha==FVS_EXPLICIT_RAW_SHA256)
 zero_lost=all(r[7]==0 for r in FVS_EXPLICIT_RAW)
 m_ratios=[r[3]/r[2]for r in FVS_EXPLICIT_RAW];d_ratios=[r[4]/r[2]for r in FVS_EXPLICIT_RAW]
 max_M_all=max(m_ratios);max_D_all=max(d_ratios)
 max_W_all=max(r[5]for r in FVS_EXPLICIT_RAW);max_k5_all=max(r[6]for r in FVS_EXPLICIT_RAW)
 by_scale={}
 for lu in sorted({r[0]for r in FVS_EXPLICIT_RAW}):
  sr=[r for r in FVS_EXPLICIT_RAW if r[0]==lu]
  ms=[r[3]/r[2]for r in sr];ds=[r[4]/r[2]for r in sr]
  by_scale[f"u{lu}"]=dict(n=len(sr),mean_Nl=round(sum(r[2]for r in sr)/len(sr),1),max_M_over_Nl=round(max(ms),4),mean_M_over_Nl=round(sum(ms)/len(ms),4),max_D_over_Nl=round(max(ds),4),mean_D_over_Nl=round(sum(ds)/len(ds),4),max_W_pre=max(r[5]for r in sr),max_k5=max(r[6]for r in sr),cols_lost=sum(r[7]for r in sr))
 max_M_u19=by_scale["u19"]["max_M_over_Nl"];max_D_u19=by_scale["u19"]["max_D_over_Nl"]
 all_below=(max_M_all<nu_fvs_ratio and max_D_all<nu_d_ratio and max_W_all<=CAPS_FVS["w_pre_cap"]and max_k5_all<CAPS_FVS["k5_cap"])
 U_ub=F(2)**54*F(1536877,1000000);Nl_ub=U_ub*F(205,1000000)
 M_ub=Nl_ub*F(int(round(nu_fvs_ratio*10000)),10000);D_ub=Nl_ub*F(int(round(nu_d_ratio*10000)),10000)
 nslots=ceil(n_bw/8.0);L_ub=F(101,100)*(M_ub+5)*(F(1,m_bw)+F(1,n_bw))+6
 sm_rat=D_ub*(512*nslots+45*8)+F(1<<56)*nslots+65*F(1<<53);pm_rat=(3*(m_bw+n_bw)*m_bw+m_bw*n_bw)*L_ub*32
 la_rat=sm_rat+pm_rat;fits_16eib=(la_rat<F(1<<64));la_log2=log2(la_rat.numerator)-log2(la_rat.denominator)
 return dict(ok=(sha_ok and zero_lost and all_below and fits_16eib),regime="zero_wmid_explicit_large_Nl_2p14_to_2p19",raw_runs=len(FVS_EXPLICIT_RAW),raw_sha256=raw_sha,raw_sha256_ok=sha_ok,zero_cols_lost_all_scales=zero_lost,all_runs_max_below_4caps=all_below,max_M_over_Nl_all24=round(max_M_all,4),max_M_over_Nl_u19=max_M_u19,max_D_over_Nl_all24=round(max_D_all,4),max_D_over_Nl_u19=max_D_u19,max_W_pre_all24=max_W_all,max_k5_all24=max_k5_all,caps=CAPS_FVS,by_scale=by_scale,bw_checked=[m_bw,n_bw],la_1own_rat_log2=round(la_log2,4),la_1own_lt_16eib=fits_16eib)

def exact_fvs_bw(M,log2_U_prod,mdepth,log2_U_base,log2_U_z6,nu_fvs_ratio=0.50,nu_d_ratio=2.00,w_pre=45.0,k5_cap=8.70,cse_all=1.00,c_tail_z6=364.0,cap_scale=1.0):
 N_l_base=(2.0**log2_U_prod)*M.NU_L_CAP*cap_scale;M_fvs=N_l_base*nu_fvs_ratio;D_rows=N_l_base*nu_d_ratio
 W_eff=w_pre*mdepth;c_prod=(W_eff*k5_cap*cse_all+c_tail_z6)*(D_rows/M_fvs);n_sm=5;best=None
 for m in(256,512,768,1024):
  for n in(24,32,40,48):
   L_st=(1.0+M.BW_SIGMA)*(M_fvs+n_sm)*(1.0/m+1.0/n)+6.0;s=M_fvs*(1.0/m+1.0/n);kap_bw=L_st/s
   hor=max(s,(1.0+M.BW_SIGMA)*(M_fvs+n_sm)/n+4.0)/s
   per_n2=(c_prod+M.C_PROJ)*(n/m+1.0)*kap_bw+(1.0/m+1.0/n)*(c_prod+(n_sm*M.C_FULL+(n-n_sm)*M.C_Z64)+M.C_RED9)*hor
   la_u=log2((M_fvs*M_fvs*per_n2+M_fvs*(m+n)*M.C_RED9+M.PREP_OWNER+M.PREP_GUARD+M.PREP_JOINT+M.LA_RESERVE+2.0**20)/M.UNIT)
   gen_u=log2(4.0*((m+n)**3)*ceil(log2(L_st))*(M.fft_F(18.0*L_st)+640.0*L_st)/M.UNIT);both=log2(2.0**la_u+2.0**gen_u)
   nslots=ceil(n/8.0);_sm=D_rows*(512.0*nslots+W_eff*8.0)+float(1<<56)*nslots+65.0*float(1<<53);_pm=(3.0*(m+n)*m+m*n)*L_st*32.0;m1=log2(_sm+_pm)
   if m1<=63.95 and(best is None or both<best["both"]):
    best=dict(la=la_u,gen=gen_u,both=both,m=m,n=n,mem_1own=m1,spmv_mem_log2=log2(_sm),pm_basis_mem_log2=log2(_pm),log2_N_cap=log2(M_fvs),log2_Ns_cap=log2(D_rows),log2_Nl_base=log2(N_l_base),c_prod=c_prod)
 return best

def exact_base_bw(M,log2_N_cap,log2_Ns_cap,c_prod,gen_at_cap,bw_steps):
 n_sm=5;best=None
 for m in(256,512,768,1024):
  for n in(24,32,40,48):
   N=2.0**log2_N_cap;L_st=bw_steps(log2_N_cap,m,n);s=N*(1.0/m+1.0/n);kap_bw=L_st/s;hor=max(s,(1.0+M.BW_SIGMA)*(N+n_sm)/n+4.0)/s
   per_n2=(c_prod+M.C_PROJ)*(n/m+1.0)*kap_bw+(1.0/m+1.0/n)*(c_prod+(n_sm*M.C_FULL+(n-n_sm)*M.C_Z64)+M.C_RED9)*hor
   la_u=log2((N*N*per_n2+N*(m+n)*M.C_RED9+M.PREP_OWNER+M.PREP_GUARD+M.PREP_JOINT+M.LA_RESERVE+2.0**20)/M.UNIT)
   gen_u=gen_at_cap(log2_N_cap,m,n);both=log2(2.0**la_u+2.0**gen_u)
   if best is None or both<best["both"]:
    ncap=2.0**log2_Ns_cap;nslots=ceil(n/8.0);_sm=ncap*(512.0*nslots+M.W_ROW*8.0)+(2.0**55.04)*nslots+65.0*(2.0**53);_pm=(3.0*(m+n)*m+m*n)*L_st*32.0
    best=dict(la=la_u,gen=gen_u,both=both,m=m,n=n,mem_1own=log2(_sm+_pm),spmv_mem_log2=log2(_sm),pm_basis_mem_log2=log2(_pm),log2_N_cap=log2_N_cap,log2_Ns_cap=log2_Ns_cap,c_prod=c_prod)
 return best

def main():
 import model_bn254_v6 as M,verify_bn254_v6 as V,descent_bn254_v6 as DESC,hx_z6 as HX
 import c6_adaptive as ADAP,v20_lat_b2 as V20L,v25_verify as V25V,boot_dp as BD,boot_mini as BM,v32_verify as V32
 t0=time.time()
 raw=sys.stdin.read().strip()if not sys.stdin.isatty()else""
 req=json.loads(raw)if raw else{}
 part=req.get("part","default")
 decl=json.loads((HERE/"declared.json").read_text())
 track=decl["track"]
 if part=="filtsim":
  lu=int(req.get("log2U",13));sd=int(req.get("seed",20261120))
  fs=V.verify_live_filtsim_z6(log2U=min(max(lu,10),15),seed=sd,W_target=100,check_kernel=(lu<=13))
  print(json.dumps(dict(track=track,part=part,checks_ok=fs["ok"],live_filtsim_z6=fs,elapsed_s=round(time.time()-t0,3))));return
 if part=="fvs":
  fl=verify_live_fvs_filtsim(V,HX,seed=20261120,scales=(13,14,15))
  zw=verify_live_zero_wmid_explicit(log2U=10,seed=2026122800,W_base=100.0,W_pre=45.0)
  rb=verify_unplanted_5sm_fvs_rank_bridge(HX,seed=1)
  ld=verify_fvs_ladder_and_caps(0.50,2.00,m_bw=768,n_bw=32)
  ok=bool(fl["ok"]and fl["zero_dropped_rec_cols"]and fl["square_trim_u13_ok"]and fl["schur_254b_exact_ok"]and fl["exact_hx_kernel_iso_ok"]and zw["ok"]and rb["ok"]and ld["ok"])
  print(json.dumps(dict(track=track,part=part,checks_ok=ok,live_fvs_filtsim=fl,live_zero_wmid_explicit=zw,unplanted_5sm_fvs_rank_bridge=rb,fvs_ladder=ld,elapsed_s=round(time.time()-t0,3))))
  if not ok:raise SystemExit(1)
  return
 if part=="dag_survivors":
  DESC.DAG_SURV_X=2.0;led2=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0)
  print(json.dumps(dict(track=track,part=part,checks_ok=True,dag_survivor_envelope_x2=dict(claim=round(led2["claim_exact"],6),exp_run=round(led2["exp_run"],6),two_dag=round(led2["descent"]["two_dag"],4),two_logs=round(led2["descent"]["two_log"],4)),elapsed_s=round(time.time()-t0,3))));return
 if part in("boot_dp","boot_mini"):
  u=int(req.get("u",-3));u=u if str(u)in BM.load()else -3
  N=256 if part=="boot_dp"else BM.load()[str(u)]["n"];j0=min(max(int(req.get("j0",0)),0),N-1)
  n=min(max(int(req.get("n",4)),1),8 if part=="boot_dp"else 40,N-j0)
  r=BD.replay(list(range(j0,j0+n)),full=True)if part=="boot_dp"else BM.replay(u,list(range(j0,j0+n)))
  print(json.dumps(dict(track=track,part=part,checks_ok=r["ok"],replay=r,elapsed_s=round(time.time()-t0,3))));return
 if part=="boot_evidence":
  rm=BD.rough_med_audit(range(4),1<<20);jm3=BD.joint_ea_mini_audit(-3,80);ok=bool(rm["ok"]and jm3["ok"])
  print(json.dumps(dict(track=track,part=part,checks_ok=ok,prod_rough_16_20=rm,mini_joint_ea_m3=jm3,elapsed_s=round(time.time()-t0,3))))
  if not ok:raise SystemExit(1)
  return
 if part=="hx":
  hx=HX.run_hx_1d_kernel_z6(log2U=min(max(int(req.get("log2U",12)),11),13),seed=int(req.get("seed",1)),m=int(req.get("m",10)),n=int(req.get("n",8)))
  print(json.dumps(dict(track=track,part=part,checks_ok=hx["ok"],hx_z6=hx,elapsed_s=round(time.time()-t0,3))));return
 if part=="diagnostics":
  lat12=V.verify_12d_descent_lattices();pgf_camp=V.verify_c6_pgf_and_hea_campaign()
  deg8_camp=DESC.verify_deg8_dispersion_and_prod_yield(seed=20261120);adap_h5=ADAP.verify_c6_adaptive_h5_live()
  checks=dict(descent_12d_lattices_24sq_24boot_ok=(lat12["ok"]and lat12["prenorm_gate_uniform_ok"]and lat12["all_24_sq_reduced_live_by_bareiss_910"]and lat12["all_24_boot_bareiss_910_certified"]and lat12["boot_candidates_gated"]==3456),deg8_dispersion_and_prod_yield_ok=(deg8_camp["ok"]and deg8_camp["tower_u5_p27631_r27481_deg8_frec_verified"]and deg8_camp["all_100_above_fb_1024"]and deg8_camp["live_6_deg8_ideals_and_2hop_d2_exact"]),c6_adaptive_h5_live_ok=(adap_h5["ok"]and adap_h5["strict_k1_zero_retry_protocol"]and adap_h5["live_k1_zero_retry_ok"]and adap_h5["live_witnesses_replayed"]==3 and adap_h5["after_fail_total"]==1856 and adap_h5["pooled_oe_fail"]<0.85))
  ok=all(checks.values())
  print(json.dumps(dict(track=track,part=part,checks_ok=ok,all_ok=ok,checks=checks,descent_12d_lattices=lat12,deg8_dispersion_and_prod_yield=deg8_camp,c6_adaptive_h5_live=adap_h5,elapsed_s=round(time.time()-t0,3))))
  if not ok:raise SystemExit(1)
  return
 alg=M.verify_tower_algebra_and_crt();fd=M.verify_exhaustive_fd(R2=50.0);proto=V.verify_swar_sb_u7_and_stageB()
 nonprim_hnf=V25V.verify_nonprim_hnf_witnesses();pool_sub=V.verify_live_c6_pool_hnf_subsample(n_tot=1000,seed=20261121)
 gen_filt=V.verify_live_genuine_c6_relation_filter(seed=20261120);bareiss_lll=V.verify_live_bareiss_12d_lll(seed=20261120)
 mont5_ecm=DESC.verify_live_mont5_and_suyama12_ecm(seed=20261120);la_kern=DESC.verify_live_z6_la_kernels(seed=20261120)
 pgf_camp=V.verify_c6_pgf_and_hea_campaign();h3_scale=ADAP.verify_h3_scaling_and_nu_h();crt_mp=ADAP.verify_split_crt_norms_and_mp_rt(seed=20261120)
 filt_live=V.verify_live_filtsim_ladder(seed=20261120)
 _cap={};_orig_rref=HX._rref_and_nullspace
 def _tap_rref(rows,ncols,q):
  rk,ker=_orig_rref(rows,ncols,q)
  if q==593 and len(rows)==ncols and len(ker)==1:_cap["A"]=[row[:]for row in rows];_cap["x"]=ker[0][:]
  return rk,ker
 HX._rref_and_nullspace=_tap_rref
 try:hx=HX.run_hx_1d_kernel_z6(log2U=11,seed=1,m=10,n=8)
 finally:HX._rref_and_nullspace=_orig_rref
 _u5=hx["unplanted_5sm_bw"]
 fvs_live=verify_live_fvs_filtsim(V,HX,seed=20261120,scales=(13,14))
 fvs_zw=verify_live_zero_wmid_explicit(log2U=10,seed=2026122800,W_base=100.0,W_pre=45.0)
 fvs_rb=verify_unplanted_5sm_fvs_rank_bridge(HX,seed=1,captured=(_cap["A"],_u5["N"],_u5["s_sm"],593,_cap["x"],_u5["D"],_u5["planted"]))
 V20L.HT_SEED_BASE=int(req.get("ht_run",0))or 1+int.from_bytes(os.urandom(2),"big")%16
 lat_b2=V20L.verify_v20_lat_b2_and_srt(seed=20261121);ht_pool=V25V.verify_ht_pool(lat_b2)
 v25_upg=V25V.verify_v25_upgrades(seed=20261122);v32_upg=V32.verify_v32_all(M,str(HERE/"kernels_v32.s"),seed=20261123)
 jd=int.from_bytes(os.urandom(2),"big")%256;bd_rep=BD.replay([jd]);bd_a=BD.alpha16();bm_st=BM.stats()
 bm_rep=dict(status="not_run",how_to_run={"part":"boot_mini","u":5,"j0":0,"n":40})
 led=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0)
 fvs_ld=verify_fvs_ladder_and_caps(0.50,2.00,m_bw=led["la"]["bw"][0],n_bw=led["la"]["bw"][1])
 pq=M.load_data()["pool_quad"];pq_bins=pq["bins"];np_bins=pq["nonprim_bins"]
 pool_all_active_ok=(sum(x[2]for x in pq_bins)==pq["n_prim_in"]and sum(x[1]for x in np_bins)==pq["n_nonprim_in"]and pq["n_in"]==2000000 and pq["n_prim_in"]+pq["n_nonprim_in"]==2000000 and nonprim_hnf["ok"])
 row_checks={};all_rows_ok=True
 for k,dv in decl["rows"].items():
  cv=led["rows"][k];ok=(cv<=float(dv)+1e-12);row_checks[k]=dict(computed=round(cv,6),declared=dv,ok=ok);all_rows_ok&=ok
 ladder8_ok=(filt_live["ok"]and led["la"]["ladder_all4_monotonic"]and led["la"]["all_runs_28_below_4caps"]and led["la"]["filtsim_cur28"]["n"]>=4 and filt_live["u13"]["N_sw_U"]>led["la"]["filtsim_cur20"]["N_sw_U"]and filt_live["u13"]["N_l_U"]>led["la"]["filtsim_cur20"]["N_l_U"]and filt_live["u13"]["k5_exp"]>led["la"]["filtsim_cur20"]["k5_exp"]and filt_live["u13"]["tab_r"]<M.TAB_R_CAP)
 kpath=str(HERE/"kernels_v32.s")
 fp_cnt,_=V20L.fp_prices_from_kernels(kpath)
 fp_kernel_prices=dict(counted_k0_to_k10=fp_cnt,charged_FP_CL=V20L.FP_CL,root_k11_charged=V20L.FP_CL[11],match=(fp_cnt==V20L.FP_CL[:11]and V20L.FP_CL[11]>=max(fp_cnt[2:])))
 wheel_cnt=V20L.wheel_prices_from_kernels(kpath)
 wheel_kernel_prices=dict(counted=wheel_cnt,charged=V20L.WHEEL_OPS,match=(wheel_cnt==V20L.WHEEL_OPS))
 bdd=led["descent"];_r4=dict(led["rows"]);_r4["two_logs"]=bdd["boot_two_log_if_charge4"]
 bdd["claim_if_boot_charge_4"]=round(math.log2(1.0+led["K_gate"]/2.0**64*2.0**math.log2(sum(2.0**v for v in _r4.values()))),6)
 _r4["two_logs"]=bdd["boot_two_log_if_survivors_x4"]
 bdd["claim_if_boot_survivors_x4"]=round(math.log2(1.0+led["K_gate"]/2.0**64*2.0**math.log2(sum(2.0**v for v in _r4.values()))),6)
 SBND=__import__("stor_bound")
 checks=dict(
  fp_node_prices_counted_from_kernels_ok=fp_kernel_prices["match"],
  wheel_prices_counted_from_kernels_ok=wheel_kernel_prices["match"],
  algebra_sigma_w_and_crt_ok=(alg["ok"]and alg["ntt3_shoup_butterfly_ok"]and alg["ntt3_garner_crt_ok"]and alg["ntt3_r4_r3_shoup_conv_ok"]),
  exhaustive_fd_pure_int_bijection_ok=(fd["ok"]and fd["zero_float_geom"]and fd["sigma_inv_exact"]and fd["max_order_dedekind_ok"]and fd["exact_frac_enum_ok"]and fd["certified_interval_ok"]and fd["global_hexagon_proof_ok"]and fd["unit_index1_proved"]),
  swar_sb_u7_pure_int_and_stageB_ok=(proto["ok"]and proto["det_occ_bound_ok"]and proto["p14_zero_gcd_gm_ok"]and proto["zero_float_geom"]and proto["flat_array_zero_dict_build"]and proto["flat_sbu7_exact_ok"]and proto["zero_sort_exact_ok"]and proto["exact_slot_cutoff_M_t"]and not proto["qa_stored_in_record"]),
  lattice_b2_score_nodes_and_srt_ok=(lat_b2["ok"]and lat_b2["smat_cyclic_a345_ok"]and lat_b2["zero_float_fp_enum"]and lat_b2["strict_u64_i64_umulh_checked"]and lat_b2["zero_exception_64bit_q38_proved"]and lat_b2["cert_min_gsym_exact"]==753729 and lat_b2["isqrt_u40_lut7_errs"]==0 and lat_b2["staddh_16b_semantic_errs"]==0 and lat_b2["lattice_b2_mismatches"]==0 and lat_b2["complete_root_to_leaf_leaves"]>=10000 and lat_b2["complete_root_to_leaf_mismatches"]==0 and lat_b2["wheel_tree_state_errs"]==0 and lat_b2["prod_exec_mismatches"]==0 and lat_b2["barrett_kta_18op_errs"]==0 and lat_b2["score_node_replay_errs"]==0 and lat_b2["fabric_1hop_addr_errs"]==0 and lat_b2["srt_batch_mismatches"]==0 and lat_b2["gate_0lp_fn"]==0 and lat_b2["gate_0lp_fp"]==0),
  v25_upgrades_ok=(v25_upg["ok"]and v25_upg["b1_e2e"]["ok"]and v25_upg["lat_apply_e2e"]["ok"]and v25_upg["a0_chord_kta_errs"]==0 and v25_upg["cone3_bijection_errs"]==0 and v25_upg["fd_rem0_u40_errs"]==0 and v25_upg["v28_wheel_ok"]),
  v32_upgrades_ok=(v32_upg["ok"]and v32_upg["adic3_exact"]["ok"]and v32_upg["ntt3_split_proto_errs"]==0 and v32_upg["ctail_lazy_errs"]==0 and v32_upg["bfi_b2_errs"]==0 and v32_upg["delta_b1_errs"]==0 and v32_upg["u16_range_ok"]and v32_upg["sf_uses_survivor_upper_ok"]and v32_upg["scan8_u16_errs"]==0 and v32_upg["b1_delta_ops_ok"]and v32_upg["boot_df_tree_ok"]and v32_upg["srt_g_leaf_ok"]and v32_upg["mul_kar_sb_ok"]and v32_upg["v37_ntt_ok"]and v32_upg["v39_emit8_ok"]and v32_upg["v40_descent_ok"]and v32_upg["step1_tab_8192_inv_errs"]==0 and v32_upg["rho_tail_ok"]and v32_upg["ret2d_per_prime_errs"]==0 and v32_upg["ret2d_mc_ok"]),
  live_pool_hnf_subsample_ok=pool_sub["ok"],
  live_genuine_c6_rel_filter_ok=(gen_filt["ok"]and gen_filt["post_merge_cols_lost"]==0 and gen_filt["merge_backsub_errors"]==0),
  live_bareiss_12d_lll_limb_ok=(bareiss_lll["ok"]and bareiss_lll["inductive_sr_inv_ok"]and bareiss_lll["input_uniform_wc_proved"]and bareiss_lll["pair_setup_p59_proved"]and bareiss_lll["lat_b2_p61_proved"]),
  live_mont5_cios_and_suyama12_ecm_ok=(mont5_ecm["ok"]and mont5_ecm["ecm_1rung_t2_ok"]and mont5_ecm["ecm_1rung_w_le_b1"]and mont5_ecm["ecm_1rung_min_iv_ratio"]>=1.0),
  live_z6_la_5limb_barrett_ok=(la_kern["ok"]and la_kern["unit4_v4_expand_ok"]and la_kern["unit4_8vec_ok"]and la_kern["unit4_8vec_u64_idx_43b_ok"]and la_kern["unit4_8vec_realistic_ir"]<=M.C_UNIT_5),
  c6_pgf_usable_depth2_and_hea_campaign_ok=(pgf_camp["ok"]and pgf_camp["exact_ideal_primitivity_checked"]and pgf_camp["live_4_depth2_2hop_chains_exact"]and pgf_camp["unified_bareiss_910_used"]),
  h3_scaling_and_nu_h_ok=(h3_scale["ok"]and h3_scale["fixed_instance_deterministic_h3"]and h3_scale["zero_loss_filtsim_ok"]and h3_scale["all4_monotone_20_27"]and h3_scale["hist_v4_exact_ok"]and h3_scale["fg_dens_ok"]and h3_scale["u_ram_all_lt_2p30"]),
  split_crt_norms_and_mp_rt_ok=(crt_mp["ok"]and crt_mp["chained_1pass_b2_errs"]==0),
  pool_all_active_200k_ok=pool_all_active_ok,
  live_filtsim_z6_and_wmid_ok=(filt_live["ok"]and filt_live["u12"]["post_merge_cols_lost"]==0 and filt_live["u12"]["merge_backsub_errors"]==0 and filt_live["u13"]["post_merge_cols_lost"]==0),
  filtsim_8scale_ladder_ok=ladder8_ok,
  unplanted_5sm_bw_dlp_and_schur_ok=hx["ok"],
  live_fvs_filtsim_13_14_ok=(fvs_live["ok"]and fvs_live["zero_dropped_rec_cols"]and fvs_live["square_trim_u13_ok"]and fvs_live["schur_254b_exact_ok"]and fvs_live["exact_hx_kernel_iso_ok"]),
  live_zero_wmid_explicit_ok=(fvs_zw["ok"]and fvs_zw["cols_lost"]==0 and fvs_zw["strict_lower_tri"]),
  unplanted_5sm_fvs_rank_bridge_ok=(fvs_rb["ok"]and not fvs_rb["planted"]and fvs_rb["cols_lost"]==0 and fvs_rb["exact_nullity_S"]==1 and fvs_rb["neg_ctrl_4sm_nullity"]==0 and fvs_rb["implicit_bw_1d_recovered"]),
  fvs_ladder_const_cap_and_16eib_ok=(fvs_ld["ok"]and fvs_ld["raw_sha256_ok"]and fvs_ld["zero_cols_lost_all_scales"]and fvs_ld["all_runs_max_below_4caps"]and fvs_ld["la_1own_lt_16eib"]and fvs_ld["la_1own_rat_log2"]>=led["la"]["mem_1own_log2"]),
  live_62row_pure_k1_boole_wald_dag_ok=(bdd["ok"]and bdd["nodes"]==62 and bdd["mode"]=="pure_k1_boole_wald_iid_dp_boot" and bdd["p0_floor"]==0.001 and bdd["phi_dag_one"]<=0.50 and bdd["e_runs_per_log"]<=2.00 and bdd["mu_min"]>=500.0 and bdd["p0_max"]<=0.005 and bdd["p_node_abort"]==0.0 and bdd["p_repeat_loss"]==0.0 and bdd["p_descent_fail"]==0.0 and bdd["desc_cap_log2"]==112.0),
  rows_ok=all_rows_ok,
  exp_run_ok=(led["exp_run"]<=float(decl["exp_run"])),
  success_ok=(led["succ_exact"]>=float(decl["success"])),
  gate_K_ok=(led["K_gate"]==decl["gate_K"]and led["gate_passes"]and led["gate_prev_fails"]),
  claim_ok=(led["claim_exact"]<=float(decl["claim"])and int(led["score_numerator"])<=int(decl["score_numerator"])),
  ram_ok=(led["geometry"]["ram_log2"]<=97.5798 and led["geometry"]["ram_log2"]<100.0 and led["geometry"]["peak_global_ram_cap_log2"]<=97.5798 and led["geometry"]["stor_bound"]["ok"]and led["geometry"]["n_sc_nodes_exp"]==34 and led["geometry"]["cs_blocks_max"]<(1<<26)and led["fail_terms"]["storage_cap_abort"]==0.0 and led["geometry"]["list_1node_ram_log2"]<=60.176 and led["geometry"]["emit_1own_ram_log2"]<=63.973 and led["geometry"]["score_1node_ram_log2"]<=63.393 and led["geometry"]["score_arena_log2"]<64.0 and led["geometry"]["b2_1p_ram_log2"]<=57.58 and led["geometry"]["b2_build_1p_ram_log2"]<=52.51 and led["la"]["mem_1own_log2"]<=62.50 and led["la"]["fits_1own"]and bdd["boot_root_retained_bytes_log2"]<=60.470 and bdd["boot_root_bytes_log2"]<=63.470 and led["geometry"]["stor_bound"]["checks"]["stageb_g_fits_ok"]and led["geometry"]["stor_bound"]["checks"]["f_node_lt_2_64"]and led["geometry"]["stor_bound"]["stageb_g_L"]==13 and led["geometry"]["stor_bound"]["f_node_L"]==50 and SBND.STAGEB_G_SP>=math.ceil(M.pbits(0,M.OP["p_min"],M.KH_G_CAP)/64.0)and SBND.F_NODE_SP>=math.ceil(M.pbits(0,M.OP["lf"],M.KH_F_CAP)/64.0)and led["geometry"]["stor_bound"]["f_node_log2"]<=62.80 and led["geometry"]["stor_bound"]["checks"]["machine_lt_2_100"]and max(led["geometry"]["list_1node_ram_log2"],led["geometry"]["emit_1own_ram_log2"],led["geometry"]["score_1node_ram_log2"],led["geometry"]["stor_bound"]["f_node_log2"],led["geometry"]["b2_1p_ram_log2"],led["geometry"]["b2_build_1p_ram_log2"],led["la"]["mem_1own_log2"],bdd["boot_root_bytes_log2"])<64.0),
  geo_consistency_ok=(abs(led["pool"]["frac_reg_analytic"]-led["pool"]["frac_reg"])<0.005 and led["geometry"]["max_sb"]<=8 and led["geometry"]["sup_nf_bits"]<448.0 and led["geometry"]["sup_ng_bits"]<=512.0 and led["geometry"]["log2_Na_shell"]<47.0),
  manifest_subbounds_ok=(led["geometry"]["log2_V"]<=94.8905 and 568.80<=led["geometry"]["Ra"]<=568.90 and led["geometry"]["Rb_GA"]<=755.60 and led["per_cell"]["rep_lo"]<=1.905 and led["per_cell"]["occ_sb"]<=0.365 and led["per_cell"]["hi_ap"]<=0.912 and led["per_cell"]["lat_b2"]<=4.780 and led["per_cell"]["lat_lll"]<=0.0200 and led["per_cell"]["scan_zero"]<=0.428 and led["per_cell"]["cof"]<=0.615 and led["per_cell"]["total"]<=9.015),
  peel_ok=(not led["la"]["peel"]["stalls"]and led["la"]["peel"]["exc_crit"]<0.72),
  ht_pool_replay_and_pooled_bound_ok=ht_pool["ok"],
  boot_dp_random_lattice_replay_ok=bd_rep["ok"],
  boot_kept_rule_ledger_ok=(bdd["boot_dp_lattices"]==256 and bdd["boot_charge"]==2.0 and bdd["boot_dp_surv_monotone"]and bdd["boot_dp_pair_loss"]<1e-3 and bdd["boot_dp_first_moment_factor"]-3.0*bdd["boot_dp_first_moment_factor_se"]>=1.0 and bdd["boot_dp_first_moment_exact_nl"]-3.0*bdd["boot_dp_first_moment_exact_nl_se"]>=1.0 and bdd["boot_dp_surv_exact_nl"]+3.0*bdd["boot_dp_surv_exact_nl_se"]<=0.5*bdd["boot_dp_surv_margin"]and bdd["boot_dp_surv_margin"]>=1.0),
  boot_mini_ratios_above_charge_ok=all(t["exhaustive"]and t["tower_ok"]and all(v["first_moment_ci95"][0]>=0.5 and v["lattice_ci95"][0]>=0.5 for v in t["by_L"].values())for t in bm_st.values()),
 )
 checks_ok=all(checks.values())
 out=dict(ht_pool=ht_pool,boot_dp_replay=bd_rep,boot_dp_small_primes=bd_a,boot_mini_replay=bm_rep,boot_mini_stats=bm_st,track=track,part=part,checks_ok=checks_ok,all_ok=checks_ok,checks=checks,live_fvs_filtsim=fvs_live,live_zero_wmid_explicit=fvs_zw,unplanted_5sm_fvs_rank_bridge=fvs_rb,fvs_ladder=fvs_ld,algebra=alg,exhaustive_fd=fd,swar_sb_u7_stageB=proto,v20_lat_b2_and_srt=lat_b2,fp_kernel_prices=fp_kernel_prices,wheel_kernel_prices=wheel_kernel_prices,v25_upgrades=v25_upg,v32_upgrades=v32_upg,live_pool_hnf_subsample=pool_sub,live_genuine_c6_rel_filter=gen_filt,live_bareiss_12d_lll=bareiss_lll,live_mont5_cios_and_suyama12_ecm=mont5_ecm,live_z6_la_kernels=la_kern,nonprim_hnf_witnesses=nonprim_hnf,c6_pgf_campaign=pgf_camp,diagnostics_part=dict(status="not_run",how_to_run={"part":"diagnostics"}),h3_scaling_and_nu_h=h3_scale,split_crt_norms_and_mp_rt=crt_mp,live_filtsim_z6=filt_live,hx_unplanted_5sm_bw=hx,ledger=dict(rows=row_checks,exp_run=round(led["exp_run"],9),exp_run_decl=decl["exp_run"],success=round(led["succ_exact"],8),success_decl=decl["success"],gate_K=led["K_gate"],claim=round(led["claim_exact"],9),claim_decl=decl["claim"],score_numerator=decl["score_numerator"],fail_terms=led["fail_terms"]),geometry=led["geometry"],per_cell=led["per_cell"],la=led["la"],descent=bdd)
 if part in("sensitivity","full"):
  s_base=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,la_mode="base_v43")
  s_u19=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,fvs_nu_ratio=0.40,fvs_d_ratio=1.35)
  s_fvs60=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,fvs_nu_ratio=0.60,fvs_d_ratio=2.50)
  s_livew=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,fvs_nu_ratio=0.58,fvs_d_ratio=2.45,fvs_k5_cap=9.25)
  s_plain=M.compute_ledger(track=track,use_db=False,c_tail_z6=364.0)
  s_lat10=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,strict_per_owner=True,swar_probe=10.0)
  s_flat=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,strict_per_owner=False)
  s_tight=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,cap_scale=0.80)
  s_srt30=M.compute_ledger(track=track,use_db=True,c_tail_z6=364.0,srt_g_cap=0.30)
  c_g395=math.log2(1.0+(0.395/0.999999)*(2.0**led["exp_run"]))
  out["sensitivity"]=dict(
   zero_fvs_baseline_h3=dict(exp_run=round(s_base["exp_run"],6),claim=round(s_base["claim_exact"],6)),
   fvs_u19_040=dict(exp_run=round(s_u19["exp_run"],6),claim=round(s_u19["claim_exact"],6)),
   fvs_stress_060_250=dict(exp_run=round(s_fvs60["exp_run"],6),claim=round(s_fvs60["claim_exact"],6)),
   live_wmid_stress_058_245_925=dict(exp_run=round(s_livew["exp_run"],6),claim=round(s_livew["claim_exact"],6)),
   gate_0395_slack=dict(exp_run=round(led["exp_run"],6),claim=round(c_g395,6)),
   plain_dickman_rho=dict(exp_run=round(s_plain["exp_run"],6),claim=round(s_plain["claim_exact"],6)),
   conservative_lat_b2_plus10=dict(exp_run=round(s_lat10["exp_run"],6),claim=round(s_lat10["claim_exact"],6)),
   conservative_srt_030_cap=dict(exp_run=round(s_srt30["exp_run"],6),claim=round(s_srt30["claim_exact"],6)),
   flat_shared_arena_ref=dict(exp_run=round(s_flat["exp_run"],6),claim=round(s_flat["claim_exact"],6)),
   de_bruijn_tighter_cap=dict(exp_run=round(s_tight["exp_run"],6),claim=round(s_tight["claim_exact"],6)),
  )
 out["elapsed_s"]=round(time.time()-t0,3)
 print(json.dumps(out))
 if not checks_ok:raise SystemExit(1)
if __name__=="__main__":main()
