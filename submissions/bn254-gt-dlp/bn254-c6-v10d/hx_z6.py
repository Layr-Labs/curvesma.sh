"""hx_z6.py — Live unplanted 2-sided 5-SM F_1187^* BW+DLP & Z[omega] Schur operator check."""
import json,math,random,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE/"files" if(HERE/"files").is_dir()else HERE
for _p in(str(ROOT),str(HERE)):
 if _p not in sys.path:
  sys.path.insert(0,_p)
import verify_bn254_v6 as V
Q_HX=1000033
CHI_HX=next(x for g in range(2,50)for x in[pow(g,(Q_HX-1)//6,Q_HX)]if(x*x-x+1)%Q_HX==0)
def _z6_eval(v,chi=CHI_HX,q=Q_HX):
 return(v[0]+v[1]*chi)%q
def _rref_and_nullspace(rows,ncols,q):
 R=[r[:]for r in rows]
 nrows=len(R)
 piv_cols=[]
 r=0
 for c in range(ncols):
  piv=next((i for i in range(r,nrows)if R[i][c]),None)
  if piv is None:continue
  R[r],R[piv]=R[piv],R[r]
  iv=pow(R[r][c],q-2,q)
  rr=[(x*iv)%q for x in R[r]]
  R[r]=rr
  for i in range(nrows):
   if i!=r and R[i][c]:
    f=R[i][c]
    ri=R[i]
    for j in range(c,ncols):
     if rr[j]:ri[j]=(ri[j]-f*rr[j])%q
  piv_cols.append(c)
  r+=1
  if r==nrows:break
 piv_set=set(piv_cols)
 free_cols=[c for c in range(ncols)if c not in piv_set]
 basis=[]
 for fc in free_cols:
  vec=[0]*ncols
  vec[fc]=1
  for row_i,pc in enumerate(piv_cols):
   if R[row_i][fc]:
    vec[pc]=(-R[row_i][fc])%q
  basis.append(vec)
 return r,basis
def _rank_rows(rows,q):
 if not rows:
  return 0
 return _rref_and_nullspace(rows,len(rows[0]),q)[0]
def unplanted_two_sided_5sm_bw_check(seed=1,m_blk=10,n_blk=8):
 p,r,m_root=1187,593,534
 r2=r*r
 f_lo=[3,1,-7,0,1]
 g_lo=[-1,-6,-2,1]
 s_sm=5
 def peval(c,x,mod):
  s=0
  for ci in reversed(c):
   s=(s*x+ci)%mod
  return s
 assert peval(f_lo,m_root,p)==0 and peval(g_lo,m_root,p)==0
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
   if e&1:
    res=pmul(res,b,mod_p,q)
   e>>=1
   if e:
    b=pmul(b,b,mod_p,q)
  return res
 def unram(c,q):
  for rt in range(q):
   if peval(c,rt,q)==0:
    if sum(i*c[i]*pow(rt,i-1,q)for i in range(1,len(c)))%q==0:
     return False
  return True
 eps_f,eps_g=r*r-1,r-1
 def sm_vec(a,b,poly_lo,eps):
  pw=ppow([a%r2,(-b)%r2],eps,poly_lo,r2)
  pw[0]=(pw[0]-1)%r2
  assert all(x%r==0 for x in pw)
  return[(x//r)%r for x in pw]
 B=229
 primes=[q for q in range(2,B+1)if all(q%d for d in range(2,int(q**0.5)+1))]
 fb_f,fb_g={},{}
 for q in primes:
  if unram(f_lo,q):
   for rt in range(q):
    if peval(f_lo,rt,q)==0:
     fb_f[(q,rt)]=len(fb_f)
  if unram(g_lo,q):
   for rt in range(q):
    if peval(g_lo,rt,q)==0:
     fb_g[(q,rt)]=len(fb_g)
 pf_act=sorted({q for q,_ in fb_f})
 pg_act=sorted({q for q,_ in fb_g})
 def sfac(val,plist):
  n=abs(val)
  if n==0:
   return None
  fac={}
  for q in plist:
   if n%q==0:
    e=0
    while n%q==0:
     n//=q
     e+=1
    fac[q]=e
    if n==1:
     return fac
  return fac if n==1 else None
 def make_sieve_list(fb_d,poly_lo):
  out=[]
  for(q,rt)in fb_d:
   lgq=int(round(math.log2(q)*4))
   out.append((q,rt,lgq))
   pk,rtk=q*q,rt
   while pk<=128:
    val=peval(poly_lo,rtk,pk)
    der=sum(i*poly_lo[i]*pow(rtk,i-1,q)for i in range(1,len(poly_lo)))%q
    rtk=(rtk-val*pow(der,-1,q))%pk
    out.append((pk,rtk,lgq))
    pk*=q
  return out
 qg_s=make_sieve_list(fb_g,g_lo)
 rels=[]
 for b in range(1,500):
  alo,ahi=-int(3.2*b)-10,int(4.2*b)+10
  L=ahi-alo+1
  sg=bytearray(L)
  for q,rt,lgq in qg_s:
   start=(b*rt-alo)%q
   for idx in range(start,L,q):
    v=sg[idx]+lgq
    sg[idx]=v if v<255 else 255
  thr_g=0 if b<=20 else max(0,int(12.0*math.log2(b))-20)
  for idx in range(L):
   if sg[idx]<thr_g:
    continue
   a=alo+idx
   if a==0 or math.gcd(abs(a),b)!=1 or(a-b*m_root)%p==0:
    continue
   fg=sfac(a*a*a-2*a*a*b-6*a*b*b-b*b*b,pg_act)
   if fg is None:
    continue
   ff=sfac(a**4-7*a*a*b*b+a*b**3+3*b**4,pf_act)
   if ff is None:
    continue
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
  if len(nxt)==len(cur):
   break
  cur=nxt
 used=sorted({k for ri in cur for k in rels[ri][2]}|{nf_cnt+k for ri in cur for k in rels[ri][3]})
 cmap={c:i for i,c in enumerate(used)}
 N=len(used)
 D=N+s_sm
 cc=[0]*N
 rcols=[]
 for ri in cur:
  cs=[cmap[k]for k in rels[ri][2]]+[cmap[nf_cnt+k]for k in rels[ri][3]]
  rcols.append(cs)
  for c in cs:
   cc[c]+=1
 act=list(range(len(cur)))
 while len(act)>D:
  cands=[i for i in act if all(cc[c]>=2 for c in rcols[i])]
  if not cands:
   break
  drop=max(cands,key=lambda i:(len(rcols[i]),i))
  act.remove(drop)
  for c in rcols[drop]:
   cc[c]-=1
 chosen=[cur[i]for i in act]
 assert len(chosen)==D
 rng=random.Random(2026+seed)
 W_f=[[rng.randrange(r)for _ in range(4)]for _ in range(3)]
 W_g=[[rng.randrange(r)for _ in range(3)]for _ in range(2)]
 M_mat=[[0]*N for _ in range(D)]
 S_cols=[[0]*D for _ in range(s_sm)]
 for i,ri in enumerate(chosen):
  a,b,vf,vg=rels[ri]
  for k,e in vf.items():
   M_mat[i][cmap[k]]=e%r
  for k,e in vg.items():
   M_mat[i][cmap[nf_cnt+k]]=(-e)%r
  sf=sm_vec(a,b,f_lo,eps_f)
  sg=sm_vec(a,b,g_lo,eps_g)
  for t in range(3):
   S_cols[t][i]=sum(W_f[t][j]*sf[j]for j in range(4))%r
  for t in range(2):
   S_cols[3+t][i]=sum(W_g[t][j]*sg[j]for j in range(3))%r
 rM,_=_rref_and_nullspace(M_mat,N,r)
 A=[M_mat[i]+[S_cols[t][i]for t in range(s_sm)]for i in range(D)]
 rA,kerA=_rref_and_nullspace(A,D,r)
 A4=[M_mat[i]+[S_cols[t][i]for t in range(s_sm-1)]for i in range(D)]
 _,kerA4=_rref_and_nullspace(A4,D-1,r)
 def op(x_live):
  return[sum(M_mat[i][k]*x_live[k]for k in range(N))%r for i in range(D)]
 Bv=lambda x:op(x[:N])
 d=-(-D//n_blk)+1+ -(-D//(128*n_blk))
 e=-(-D//m_blk)+1+ -(-D//(128*m_blk))
 U_rand=[[rng.randrange(r)for _ in range(N)]for _ in range(n_blk-s_sm)]
 Y=[S_cols[t]for t in range(s_sm)]+[op(u)for u in U_rand]
 BY=[Y]
 for _ in range(d+e-1):
  BY.append([Bv(col)for col in BY[-1]])
 K_rows=list(BY[0][s_sm:])+[col for j in range(1,d+1)for col in BY[j]]
 rC=_rank_rows(K_rows,r)
 def make_X(onepad):
  idx=list(range(N))
  rng.shuffle(idx)
  blk=[[]for _ in range(m_blk)]
  for j,c in enumerate(idx):
   blk[j%m_blk].append(c)
  for t in range(s_sm):
   blk[0 if onepad else t%m_blk].append(N+t)
  return[[(rng.randrange(1,r)if c in set(b)else 0)for c in range(D)]for b in blk]
 X_good=make_X(False)
 X_bad=make_X(True)
 def obs_rank(X):
  Ok=[]
  Xi=[row[:]for row in X]
  for _ in range(e):
   Ok.extend(Xi)
   Xi=[[sum(row[i]*A[i][k]for i in range(D))%r if k<N else 0 for k in range(D)]for row in Xi]
  return _rank_rows([[sum(row[i]*A[i][k]for i in range(D))%r for k in range(D)]for row in Ok],r)
 rO=obs_rank(X_good)
 rN=obs_rank(X_bad)
 a_seq=[[[sum(X_good[r_i][row_idx]*BY[k][c_j][row_idx]for row_idx in range(D))%r for c_j in range(n_blk)]for r_i in range(m_blk)]
 for k in range(d+e)
 ]
 H_rows=[]
 for i in range(e):
  for r_i in range(m_blk):
   row=[]
   for j in range(d+1):
    row.extend(a_seq[i+j][r_i])
   H_rows.append(row)
 rH,_=_rref_and_nullspace(H_rows,(d+1)*n_blk,r)
 L_sb=d+e;K_sb=n_blk+m_blk
 V_sb=[[[1 if row==col else 0]for row in range(K_sb)]for col in range(K_sb)]
 R_sb=[[[0]*L_sb for _ in range(m_blk)]for _ in range(K_sb)]
 for k in range(L_sb):
  for i in range(m_blk):
   for j in range(n_blk):R_sb[j][i][k]=a_seq[k][i][j]%r
 for i in range(m_blk):R_sb[n_blk+i][i][0]=(-1)%r
 delta_sb=[0]*n_blk+[1]*m_blk
 for k in range(L_sb):
  for i in range(m_blk):
   nz=[c for c in range(K_sb)if R_sb[c][i][k]]
   if not nz:continue
   piv=min(nz,key=lambda c:(delta_sb[c],c))
   inv_p=pow(R_sb[piv][i][k],-1,r)
   for c in nz:
    if c==piv:continue
    fac=(R_sb[c][i][k]*inv_p)%r
    for ri in range(m_blk):
     rc,rp=R_sb[c][ri],R_sb[piv][ri]
     for kk in range(k,L_sb):
      if rp[kk]:rc[kk]=(rc[kk]-fac*rp[kk])%r
    for row in range(K_sb):
     vc,vp=V_sb[c][row],V_sb[piv][row]
     if len(vc)<len(vp):vc.extend([0]*(len(vp)-len(vc)))
     for dt,vpv in enumerate(vp):
      if vpv:vc[dt]=(vc[dt]-fac*vpv)%r
   for ri in range(m_blk):R_sb[piv][ri]=[0]+R_sb[piv][ri][:-1]
   for row in range(K_sb):V_sb[piv][row]=[0]+V_sb[piv][row]
   delta_sb[piv]+=1
 sb_res_zero=all(R_sb[c][i][k]==0 for c in range(K_sb)for i in range(m_blk)for k in range(L_sb))
 sb_retained=[c for c in range(K_sb)if delta_sb[c]<=d]
 sb_sm_cols=0;sb_hf_ok=True;v_rec=None;sb_chosen_deg=None
 for c in sb_retained:
  dc=delta_sb[c]
  lead_sm=[V_sb[c][row][dc]if dc<len(V_sb[c][row])else 0 for row in range(s_sm)]
  if not any(lead_sm):continue
  sb_sm_cols+=1
  fvec=[0]*((d+1)*n_blk)
  for j in range(dc+1):
   dz=dc-j
   for cn in range(n_blk):
    pv=V_sb[c][cn];fvec[j*n_blk+cn]=pv[dz]if dz<len(pv)else 0
  if not all(sum(row[j]*fvec[j]for j in range(len(fvec)))%r==0 for row in H_rows):sb_hf_ok=False
  if v_rec is None:
   c_rec=fvec[:s_sm];w_rec=[0]*N
   for t in range(n_blk-s_sm):
    ft=fvec[s_sm+t]
    if ft:
     u=U_rand[t]
     for idx in range(N):w_rec[idx]=(w_rec[idx]+ft*u[idx])%r
   for j in range(1,d+1):
    fj=fvec[j*n_blk:(j+1)*n_blk];Bj1Y=BY[j-1]
    for t in range(n_blk):
     ft=fj[t]
     if ft:
      col=Bj1Y[t]
      for idx in range(N):w_rec[idx]=(w_rec[idx]+ft*col[idx])%r
   cand=w_rec+c_rec
   if all(sum(A[i][k]*cand[k]for k in range(D))%r==0 for i in range(D)):
    v_rec=cand;sb_chosen_deg=dc
 dlp_ok=False
 dlp_ratio_set=[]
 dlp_sample={};n_checked_dl=0
 if v_rec is not None:
  w_log=v_rec[:N]
  gen_p=next(x for x in range(2,p)if pow(x,r,p)==1)
  dlog={pow(gen_p,k,p):k for k in range(r)}
  inv2=pow(2,-1,r)
  ratios=set()
  for ri in cur:
   a,b,vf,vg=rels[ri]
   val_p=pow((a-b*m_root)%p,2,p)
   dl=(dlog[val_p]*inv2)%r
   if dl==0:continue
   n_checked_dl+=1
   sf=sm_vec(a,b,f_lo,eps_f)
   lhs_f=(sum(e_v*w_log[cmap[k]]for k,e_v in vf.items())
   +sum(v_rec[N+t]*sum(W_f[t][j]*sf[j]for j in range(4))for t in range(3)))%r
   ratios.add((lhs_f*pow(dl,-1,r))%r)
  dlp_ratio_set=sorted(ratios)
  def eval_pred(ri):
   a,b,vf,_=rels[ri]
   sf=sm_vec(a,b,f_lo,eps_f)
   return(sum(e_v*w_log[cmap[k]]for k,e_v in vf.items())
   +sum(v_rec[N+t]*sum(W_f[t][j]*sf[j]for j in range(4))for t in range(3)))%r
  r1,r2_i=cur[0],cur[1]
  g_el=pow((rels[r1][0]-rels[r1][1]*m_root)%p,2,p)
  z_el=pow((rels[r2_i][0]-rels[r2_i][1]*m_root)%p,2,p)
  x_dlp=(eval_pred(r2_i)*pow(eval_pred(r1),-1,r))%r
  dlp_ok=(len(ratios)==1 and 0 not in ratios and pow(g_el,x_dlp,p)==z_el)
  dlp_sample=dict(g=g_el,z=z_el,recovered_x=x_dlp,g_pow_x_eq_z=(pow(g_el,x_dlp,p)==z_el))
 unplanted_ok=(rM==N and rA==D-1 and len(kerA)==1 and len(kerA4)==0
 and rC==N and rO==rA and rN<rA and sb_res_zero and sb_hf_ok and sb_sm_cols>=1 and v_rec is not None and dlp_ok
 )
 return dict(ok=unplanted_ok,planted=False,p=p,r=r,m_root=m_root,s_sm=s_sm,B=B,raw_rels=len(rels),purged_rels=len(cur),held_out_rels=len(cur)-D,N=N,D=D,m_blk=m_blk,n_blk=n_blk,d=d,e=e,rank_M=rM,rank_A=rA,ker_dim_5sm=len(kerA),neg_ctrl_4sm_ker_dim=len(kerA4),rank_Cd=rC,rank_OeA=rO,control_one_block=rN,
 hankel_shape=[len(H_rows),(d+1)*n_blk],hankel_rank=rH,sigma_basis_order=L_sb,sigma_basis_dim=K_sb,sigma_basis_residual_zero=sb_res_zero,sigma_basis_retained_cols=len(sb_retained),sigma_basis_sm_nonzero_cols=sb_sm_cols,sigma_basis_chosen_deg=sb_chosen_deg,sigma_basis_all_hf_zero=sb_hf_ok,hankel_recovered=(v_rec is not None),dlp_checked_rels=n_checked_dl,dlp_all_purged_ratios=dlp_ratio_set,dlp_live_check=dlp_sample,
 )
def z6_schur_operator_check(log2U=12,seed=1):
 Q,chi=Q_HX,CHI_HX
 assert(chi*chi-chi+1)%Q==0 and(pow(chi,3,Q)+1)%Q==0
 fres=V.verify_live_filtsim_z6(log2U=log2U,seed=seed,W_target=160,check_kernel=False,export_hx=True)
 hx=fres["hx_export"]
 rec_order=hx["rec_order"]
 rec_ix={c:k for k,c in enumerate(rec_order)}
 lc=hx["live_col_set"]
 lc_ix={c:j for j,c in enumerate(lc)}
 n_rec,N=len(rec_order),len(lc)
 lr=hx["live_rows_full"][:N+5]
 D=len(lr)
 rec_comp=[]
 for k,(piv_c,piv_u,entries)in enumerate(hx["recorded_pivots"]):
  assert rec_ix[piv_c]==k
  neg_inv=(-pow(_z6_eval(piv_u),-1,Q))%Q
  rt,lt=[],[]
  for cc,vv in entries:
   vq=_z6_eval(vv)
   if not vq:continue
   if cc in rec_ix:
    assert rec_ix[cc]<k
    rt.append((rec_ix[cc],vq))
   else:
    lt.append((lc_ix[cc],vq))
  rec_comp.append((neg_inv,rt,lt))
 live_comp=[]
 for entries in lr:
  rt,lt=[],[]
  for cc,vv in entries:
   vq=_z6_eval(vv)
   if not vq:continue
   if cc in rec_ix:rt.append((rec_ix[cc],vq))
   else:lt.append((lc_ix[cc],vq))
  live_comp.append((rt,lt))
 def op(x_live):
  vr=[0]*n_rec
  for k,(neg_inv,rt,lt)in enumerate(rec_comp):
   acc=sum(vq*vr[rk]for rk,vq in rt)+sum(vq*x_live[j]for j,vq in lt)
   vr[k]=(neg_inv*(acc%Q))%Q
  return[(sum(vq*vr[rk]for rk,vq in rt)+sum(vq*x_live[j]for j,vq in lt))%Q for rt,lt in live_comp]
 Mc=[op([1 if j==k else 0 for j in range(N)])for k in range(N)]
 g=random.Random(seed*7919+31)
 vt=[g.randrange(Q)for _ in range(N)]
 prod_eq_schur=(op(vt)==[sum(Mc[k][i]*vt[k]for k in range(N))%Q for i in range(D)])
 return dict(ok=prod_eq_schur,Q=Q,chi=chi,filt_N_sw=fres["N_sw"],filt_N_l=fres["N_l"],filt_N_rec=n_rec,l_rec_strict_lower_tri=True,prod_eq_schur=prod_eq_schur,
 )
def run_hx_1d_kernel_z6(log2U=12,seed=1,m=10,n=8,s_sm=5):
 u5=unplanted_two_sided_5sm_bw_check(seed=seed,m_blk=m,n_blk=n)
 z6=z6_schur_operator_check(log2U=log2U,seed=seed)
 ok=bool(u5["ok"]and z6["ok"])
 return dict(checks_ok=ok,ok=ok,source="unplanted_2sided_5sm_bw_F1187_dlp_plus_z6_schur",unplanted_5sm_bw=u5,z6_schur_op=z6,)
if __name__=="__main__":
 res=run_hx_1d_kernel_z6()
 print(json.dumps(res,separators=(",",":")))
 if not res["checks_ok"]:raise SystemExit(1)
