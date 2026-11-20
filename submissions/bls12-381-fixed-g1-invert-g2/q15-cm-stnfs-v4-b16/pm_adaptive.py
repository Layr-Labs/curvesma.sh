"""Self-contained classic PM-Basis live verifier & Q(zeta_15) real-2-core-FB adaptive H5 verifier."""
import base64,cmath,json,lzma,math,random

def _padd(a,b,p):
 n=max(len(a),len(b));o=[0]*n
 for i,v in enumerate(a):o[i]=v
 for i,v in enumerate(b):o[i]=(o[i]+v)%p
 while o and o[-1]==0:o.pop()
 return o

def _pmul(a,b,p,cap=None):
 if not a or not b:return []
 la,lb=len(a),len(b);n=la+lb-1 if cap is None else min(la+lb-1,cap);o=[0]*n
 for i,ai in enumerate(a):
  if ai:
   for j in range(min(lb,n-i)):o[i+j]=(o[i+j]+ai*b[j])%p
 while o and o[-1]==0:o.pop()
 return o

def _pmat_mul(A,B,p,cap=None):
 r,k,c=len(A),len(A[0]),len(B[0])
 C=[[[] for _ in range(c)] for _ in range(r)]
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
 P=[[[] for _ in range(t)] for _ in range(t)]
 for j in range(t):
  for i in range(t):
   v=Pcols[j][i]
   if v:P[i][j]=[0,v] if dinc[j] else [v]
 return P,[s[j]+dinc[j] for j in range(t)]

def pm_basis(F,sigma,s,p):
 """Giorgi-Jeannerod-Villard (2003) / Thome (2002) divide-and-conquer PM-Basis over F_p[z]."""
 m,t=len(F),len(F[0])
 if sigma==1:return _mbasis1([[F[i][j][0] if F[i][j] else 0 for j in range(t)] for i in range(m)],s,p)
 s1=sigma//2;s2=sigma-s1
 P1,sm=pm_basis([[F[i][j][:s1] for j in range(t)] for i in range(m)],s1,s,p)
 FP1=_pmat_mul(F,P1,p,cap=sigma)
 G=[[FP1[i][j][s1:sigma] for j in range(t)] for i in range(m)]
 P2,so=pm_basis(G,s2,sm,p)
 return _pmat_mul(P1,P2,p),so

def pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_ref):
 """Run live divide-and-conquer PM-Basis on S_k = X A^k Y (k>=1), recover exact ker(A) via Horner."""
 m=len(seq[0]);n=len(seq[0][0]);D=len(A_sq);L=hr+hc-1
 S=seq[1:1+L]
 M=[[[] for _ in range(n+m)] for _ in range(m)]
 for i in range(m):
  for j in range(n):
   poly=[S[k][i][j] for k in range(L)]
   while poly and poly[-1]==0:poly.pop()
   M[i][j]=poly
  M[i][n+i]=[r-1]
 P_b,s_out=pm_basis(M,L,[0]*n+[1]*m,r)
 def mvec(P,v):return [sum(P[i][j]*v[j] for j in range(len(v)))%r for i in range(len(P))]
 v_pm=None;deg_used=None
 for col in sorted(range(n+m),key=lambda j:s_out[j]):
  if s_out[col]>hc-1:continue
  q=[P_b[i][col] for i in range(n)]
  dg=max((len(poly)-1 for poly in q),default=-1)
  if dg<0:continue
  acc=[0]*D
  for t in range(dg+1):
   qt=[q[bi][t] if t<len(q[bi]) else 0 for bi in range(n)]
   yc=[sum(Y[ri][bi]*qt[bi] for bi in range(n))%r for ri in range(D)]
   acc=[(a+y)%r for a,y in zip(mvec(A_sq,acc),yc)]
  if any(acc) and all(x==0 for x in mvec(A_sq,acc)):
   v_pm=acc;deg_used=dg;break
 assert v_pm is not None,"PM-Basis failed to recover ker(A)!"
 idx0=next(i for i,x in enumerate(v_ref) if x)
 ratio=(v_pm[idx0]*pow(v_ref[idx0],-1,r))%r
 assert ratio!=0 and all((v_ref[i]*ratio)%r==v_pm[i] for i in range(D))
 return dict(pm_basis_order=L,pm_basis_cols=n+m,pm_basis_deg=deg_used,pm_basis_matches_true_ker=True)

def verify_threshold_counterexample():
 """Verify THRESHOLD-REPAIR counterexample showing PGF(psi) fails without independence while T_k+psi^k holds."""
 mu_raw,r_v,psi=60.0,1.0,0.10
 beta=0.5*mu_raw/r_v
 G=lambda t:(1.0+beta*(1.0-t))**(-r_v)
 for i in range(81):
  t=i/100.0
  assert 0.10*(t**1)+0.90*(t**100)<=G(t)+1e-12
 assert 0.10>G(psi)+0.05
 best=G(0.0)+psi
 for k in range(2,65):
  tk=min(0.8,(k-1)*(1.0+1.0/beta)/(k-1+r_v))
  best=min(best,min(1.0,G(tk)/(tk**(k-1)))+(psi**k))
 assert best>=0.10-1e-12
 return dict(pgf_psi_invalid=round(G(psi),5),true_fail=0.10,tr_bound=round(best,5),tr_valid=True)

def verify_adaptive_h5_live(b85_path,camp,primes_upto_fn,y_pow_fn,add_fn,sub_fn,mul_fn,scale_fn,norm_fn,lll_fn):
 """Replay all 1,178 Q(zeta_15) real-2-core-FB after-failure rows AND live-regenerate all 3 witness roots."""
 ce=verify_threshold_counterexample()
 raw=lzma.decompress(base64.b85decode(open(b85_path,"rb").read()))
 dat=json.loads(raw.decode())
 names=dat["names"];cal=dat["cal_tbl"];after_by=dat["after_by_set"]
 fstats=dat["first_stats"];foe=dat["first_oe"]
 def _pw(si,cs,cb):
  k=f"{si}:{cs}:{cb}"
  if k in cal:
   n,sf,sw=cal[k];return sf/n,sw/n
  n=sf=sw=0
  for kk,(nn,sff,sww) in cal.items():
   _,ccs,ccb=map(int,kk.split(":"))
   if ccs==cs and ccb==cb:n+=nn;sf+=sff;sw+=sww
  return sf/n,sw/n
 per_set={};tot_f=tot_ef=tot_w=tot_ew=tot_an=tot_af=tot_shp=tot_shr=0
 tot_f1_o=tot_f1_e=0.0
 st_acc={k:[0,0,0.0] for k in ("d2","d3p","c1","c2","top","sub","p1","p2")}
 for si,nm in enumerate(names):
  ar=after_by[si];an=len(ar);af=sum(x[1] for x in ar)
  shp=sum(x[3] for x in ar);shr=sum(x[4] for x in ar)
  nw_s=dat["set_cfg"][si][2]
  of=ef=ow=ew=0.0
  for flg,fl,wk,sp,sr,c1,c2 in ar:
   chs=[c1]+([c2] if c2 else [])
   ph=sum(_pw(si,c>>5,c&31)[0] for c in chs)
   wh=sum(_pw(si,c>>5,c&31)[1] for c in chs)
   of+=fl;ef+=ph;ow+=wk;ew+=wh
   hf=flg>>1
   for sk in ("d3p" if (hf&1) else "d2","c2" if c2 else "c1","top" if any((c&31)>=nw_s for c in chs) else "sub","p2" if (hf&2) else "p1"):
    st_acc[sk][0]+=1;st_acc[sk][1]+=fl;st_acc[sk][2]+=ph
  oe_f=of/ef;oe_w=ow/ew;rat1=oe_f/(foe[si][0]/foe[si][1])
  assert oe_f<1.0 and rat1<1.0 and shp==0 and shr==0
  per_set[nm]=dict(first_n=fstats[si][0],after_n=an,after_fail=af,
                   oe_fail=round(oe_f,3),rat1=round(rat1,3),oe_work=round(oe_w,3),sh_prev=shp,sh_run=shr)
  tot_an+=an;tot_af+=af;tot_f+=of;tot_ef+=ef;tot_w+=ow;tot_ew+=ew;tot_shp+=shp;tot_shr+=shr
  tot_f1_o+=foe[si][0];tot_f1_e+=foe[si][1]
 assert tot_an==1178 and tot_af==194 and tot_shp==0 and tot_shr==0
 pooled_oe_f=tot_f/tot_ef;pooled_rat1=pooled_oe_f/(tot_f1_o/tot_f1_e);pooled_oe_w=tot_w/tot_ew
 assert pooled_oe_f<0.85 and pooled_rat1<0.85 and pooled_oe_w<1.05
 assert st_acc["d3p"][1]/st_acc["d3p"][2]<1.0 and st_acc["top"][1]/st_acc["top"][2]<1.0 and st_acc["c2"][1]/st_acc["c2"][2]<1.0 and max(v[1]/v[2] for v in st_acc.values())<1.15
 conc=dat["conc_rows"]
 assert len(conc)==12 and max(x[6] for x in conc)<20.0 and max(x[7] for x in conc)<3.0
 lx_b=int(camp["lx_b"]);B_FB=1<<lx_b;small_ns_max=int(camp["small_ns_max"])
 w_el=y_pow_fn(5);w2_el=y_pow_fn(10);c_f=add_fn([2]+[0]*7,w_el);c_g=scale_fn(2,w2_el)
 sig_m=[y_pow_fn((7*j)%15) for j in range(8)]
 def _sig(v):
  o=[0]*8
  for j,c in enumerate(v):
   if c:
    sj=sig_m[j]
    for i in range(8):o[i]+=c*sj[i]
  return o
 def _nzw(el):
  s1=_sig(el);s2=_sig(s1);s3=_sig(s2);res=mul_fn(mul_fn(el,s1),mul_fn(s2,s3))
  return res[0],res[5]
 def _wr(p):
  for g in range(2,p):
   w=pow(g,(p-1)//3,p)
   if w!=1:w2=(w*w)%p;return (w,w2) if w<w2 else (w2,w)
 def _isp(n):
  if n<2:return False
  for p in (2,3,5,7,11,13,17,19,23,29,31,37):
   if n==p:return True
   if n%p==0:return False
  d,s=n-1,0
  while d%2==0:d//=2;s+=1
  for a in (2,3,5,7,11,13,17,19,23,29,31,37):
   x=pow(a,d,n)
   if x in (1,n-1):continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:break
   else:return False
  return True
 all_ps=primes_upto_fn(B_FB);sp_tbl={};ns_tbl={}
 for p in all_ps:
  if p%15==1:
   w0,w1=_wr(p)
   sp_tbl[p]=[(oi,w,pow((2+w)%p,(p-1)//3,p)==1,pow((-2*w*w)%p,(p-1)//3,p)==1) for oi,w in enumerate((w0,w1))]
  elif p<=small_ns_max:ns_tbl[p]=_wr(p) if p%3==1 else None
 ns_ps=sorted(ns_tbl.keys());sp_ps=sorted(sp_tbl.keys())
 def _pr(xs):
  xs=list(xs)
  while len(xs)>1:xs=[xs[i]*xs[i+1] if i+1<len(xs) else xs[i] for i in range(0,len(xs),2)]
  return xs[0]
 P1=_pr(ns_ps+[p for p in sp_ps if p<=256]);P2=_pr([p for p in sp_ps if p>256])
 ret_f={(int(p),int(oi)) for p,oi in camp["ret_f_split"]+camp["ret_f_ns"]}
 ret_g={(int(p),int(oi)) for p,oi in camp["ret_g_split"]+camp["ret_g_ns"]}
 bts=[];cur=[0]*8
 def _gb(idx,rem):
  if idx==8:
   if any(cur):bts.append(tuple(cur))
   return
  r=math.isqrt(rem)
  for v in range(-r,r+1):cur[idx]=v;_gb(idx+1,rem-v*v)
 _gb(0,int(dat["ball_r2"]))
 assert len(bts)==int(dat["ball_pts"])==576
 def _strip(n):
  for P in (P1,P2):
   while n>1:
    g=math.gcd(n,P)
    if g==1:break
    n//=g
  return n
 def _chk_fb(sm,c0,c1,sd):
  rem=sm;rf=ret_f if sd==0 else ret_g
  for p in ns_ps:
   if rem%p==0:
    while rem%p==0:rem//=p
    ws=ns_tbl[p]
    if ws is None:
     if (p,0) not in rf:return False
    else:
     if not any((c0+c1*w)%p==0 and (p,oi) in rf for oi,w in enumerate(ws)):return False
  for p in sp_ps:
   if rem%p==0:
    while rem%p==0:rem//=p
    ok_h=False
    for oi,w,rfok,rgok in sp_tbl[p]:
     if (c0+c1*w)%p==0:
      if not (rfok if sd==0 else rgok) or (p,oi) not in rf:return False
      ok_h=True
    if not ok_h:return False
  return rem==1
 def _fac(cof,nw,km):
  if cof==1:return []
  if cof.bit_length()>nw*km:return None
  lim=1<<nw;out=[];st=[cof]
  while st:
   v=st.pop()
   if v==1:continue
   if _isp(v):
    if v>lim or v<=B_FB or v%15!=1:return None
    if v not in out:
     if len(out)>=km:return None
     out.append(v)
   else:
    d=None
    for c in range(1,32):
     x=y=2;g=1
     while g==1:x=(x*x+c)%v;y=(y*y+c)%v;y=(y*y+c)%v;g=math.gcd(abs(x-y),v)
     if 1<g<v:d=g;break
    if d is None:return None
    st.extend((d,v//d))
  return out
 def _cry(p,oi):
  w0,w1=_wr(p);wt=w0 if oi==0 else w1
  for g in range(2,p):
   y0=pow(g,(p-1)//15,p)
   if pow(y0,5,p)!=1 and pow(y0,3,p)!=1:break
  return min(pow(y0,k,p) for k in (1,2,4,7,8,11,13,14) if pow(pow(y0,k,p),5,p)==wt)
 def _scan(sd,q,ry,rt,nw,km):
  B=[[0]*16 for _ in range(16)];B[0][0]=q;ryp=[pow(ry,i,q) for i in range(8)]
  for i in range(1,8):B[i][0]=(-ryp[i])%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*ryp[i])%q;B[8+i][8+i]=1
  R,dsq,_=lll_fn(B,99,100);assert dsq[16]==q*q;a1,b1=R[0][:8],R[0][8:]
  a2,b2=next((v[:8],v[8:]) for v in R[1:] if any(sub_fn(mul_fn(a1,v[8:]),mul_fn(b1,v[:8]))))
  a2s=[mul_fn(y_pow_fn(j),a2) for j in range(8)];b2s=[mul_fn(y_pow_fn(j),b2) for j in range(8)]
  rels=[]
  for bt in bts:
   a=list(a1);b=list(b1)
   for j,c in enumerate(bt):
    if c:
     aj,bj=a2s[j],b2s[j]
     for i in range(8):a[i]+=c*aj[i];b[i]+=c*bj[i]
   a3=mul_fn(mul_fn(a,a),a);b3=mul_fn(mul_fn(b,b),b)
   ef=sub_fn(a3,mul_fn(c_f,b3));eg=add_fn(a3,mul_fn(c_g,b3))
   c0f,c1f=_nzw(ef);c0g,c1g=_nzw(eg)
   nf=c0f*c0f-c0f*c1f+c1f*c1f;ng=c0g*c0g-c0g*c1g+c1g*c1g
   if sd==0:assert nf%q==0;nf//=q
   else:assert ng%q==0;ng//=q
   if nf==0 or ng==0:continue
   cf=_strip(nf)
   if cf>1 and cf.bit_length()>nw*km:continue
   cg=_strip(ng)
   if cg>1 and cg.bit_length()>nw*km:continue
   pf=_fac(cf,nw,km)
   if pf is None:continue
   pg=_fac(cg,nw,km-len(pf))
   if pg is None:continue
   if not _chk_fb(nf//cf,c0f,c1f,0) or not _chk_fb(ng//cg,c0g,c1g,1):continue
   ch=[];ok=True
   for csd,ps,c0,c1 in ((0,pf,c0f,c1f),(1,pg,c0g,c1g)):
    for p in ps:
     w0,w1=_wr(p);h=0
     for oi,w in enumerate((w0,w1)):
      if (c0+c1*w)%p==0:
       rhs=(2+w)%p if csd==0 else (-2*w*w)%p
       if pow(rhs,(p-1)//3,p)!=1:ok=False;break
       ryc=_cry(p,oi);rp=1;av=bv=0
       for j in range(8):av=(av+a[j]*rp)%p;bv=(bv+b[j]*rp)%p;rp=(rp*ryc)%p
       if bv==0:ok=False;break
       rtc=(av*pow(bv,-1,p))%p
       if pow(rtc,3,p)!=rhs:ok=False;break
       ch.append((csd,p,oi,ryc,rtc));h+=1
     if not ok or h==0:ok=False;break
    if not ok:break
   if ok and len(ch)<=km:rels.append(ch)
  return rels
 def _desc(sd,q,ry,rt,nw,km,dsh,rec=False):
  rels=_scan(sd,q,ry,rt,nw,km);tw=1;prev=set();seq=[]
  for ri,ch in enumerate(rels):
   if not ch:return True,tw,seq
   aok=True;cw=0;sh=0;cids=[]
   for csd,cp,coi,cry,crt in ch:
    cid=(csd,cp,coi);cids.append(cid)
    if ri>0 and cid in prev:sh=1
    cb=cp.bit_length();cnw=max(lx_b,cb-dsh);ckm=0 if cnw==lx_b else 1
    cok,cwk,_=_desc(csd,cp,cry,crt,cnw,ckm,dsh,False)
    cw+=cwk
    if not cok:aok=False
   prev.update(cids);tw+=cw
   if rec:seq.append([ri,int(not aok),cw,len(ch),sh])
   if aok:return True,tw,seq
  return False,tw,seq
 live_w_ok=0
 for si,sd,rq,roi,rry,rrt,nw,kr,dsh,exp_seq in dat["live_witnesses"]:
  _,_,got_seq=_desc(sd,rq,rry,rrt,nw,kr,dsh,True)
  assert got_seq==exp_seq,f"Adaptive witness mismatch q={rq}: {got_seq} != {exp_seq}"
  live_w_ok+=1
 return dict(
  threshold_counterexample=ce,
  first_attempts_total=sum(x[0] for x in fstats),
  after_fail_total=tot_an,after_fail_failures=tot_af,
  pooled_oe_fail=round(pooled_oe_f,4),pooled_rat1=round(pooled_rat1,4),pooled_oe_work=round(pooled_oe_w,4),
  shared_prev_total=tot_shp,shared_run_total=tot_shr,
  conc_bins_checked=len(conc),max_C_emp=max(x[6] for x in conc),max_C_bday=max(x[7] for x in conc),
  history_strata={k:dict(n=v[0],fail=v[1],oe_fail=round(v[1]/v[2],3)) for k,v in st_acc.items()},
  live_witnesses_replayed=live_w_ok,per_setting=per_set,checks_ok=True,
 )

def verify_b16_dsieve_live(data,desc,dlive,m_node,U_bls,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 _Rm=160;_M=_Rm*_Rm
 def _ksq_b(poly,bpc):
  _x=int.from_bytes(b"".join(v.to_bytes(bpc,"little") for v in poly),"little");_p=_x*_x
  _nb=max((_M+1)*bpc,(_p.bit_length()+7)//8);_bs=_p.to_bytes(_nb,"little")
  return [int.from_bytes(_bs[i*bpc:(i+1)*bpc],"little") for i in range(_M+1)]
 _t1=[0]*(_M+1);_t1[0]=1
 for _k in range(1,_Rm+1):_t1[_k*_k]=2
 _t16=_ksq_b(_ksq_b(_ksq_b(_ksq_b(_t1,4),8),8),16)
 for _rk,_cv in desc["cells_rank2"].items():
  _R=int(_rk);assert int(_cv)==(sum(_t16[:_R*_R+1])-1)//2,f"B_16(R)/+-1 mismatch R={_R}"
 assert 4*max(int(desc["cells_rank2"][str(r["R"])]) for r in dlive["rows"])<(1<<111)
 _pl=[cmath.exp(2j*math.pi*k/15.0) for k in (1,7,4,2)];_yp=[y_pow_z(k) for k in range(15)]
 _ug=[[1,0,0,1,0,0,0,0],[-1,0,1,0,0,0,0,0],[-1,1,0,0,0,0,0,0]]
 _ui=[[1,0,-1,0,0,0,1,-1],[0,-1,0,0,-1,1,-1,0],[0,-1,-1,0,-1,0,0,-1]]
 for _a,_b in zip(_ug,_ui):assert mul_phi15_z(_a,_b)==[1]+[0]*7
 _Lg=[[math.log(abs(sum(u[j]*(z**j) for j in range(8)))) for z in _pl] for u in _ug]
 _Mm=[[_Lg[i][k] for i in range(3)] for k in range(3)]
 _a,_b,_c=_Mm[0];_d,_e,_f=_Mm[1];_g,_h,_ii=_Mm[2]
 _dt=_a*(_e*_ii-_f*_h)-_b*(_d*_ii-_f*_g)+_c*(_d*_h-_e*_g)
 assert abs(8.0*abs(_dt)-4.6618207772966835)<1e-9
 _Mi=[[(_e*_ii-_f*_h)/_dt,(_c*_h-_b*_ii)/_dt,(_b*_f-_c*_e)/_dt],
      [(_f*_g-_d*_ii)/_dt,(_a*_ii-_c*_g)/_dt,(_c*_d-_a*_f)/_dt],
      [(_d*_h-_e*_g)/_dt,(_b*_g-_a*_h)/_dt,(_a*_e-_b*_d)/_dt]]
 _ym=[tuple(tuple(mul_phi15_z(y,[1 if i==j else 0 for i in range(8)])[r] for j in range(8)) for r in range(8)) for y in _yp]
 _pc=[{0:[1]+[0]*7} for _ in range(3)]
 def _upow(_idx,_n):
  _dd=_pc[_idx]
  if _n in _dd:return _dd[_n]
  _rv=mul_phi15_z(_upow(_idx,_n-1),_ug[_idx]) if _n>0 else mul_phi15_z(_upow(_idx,_n+1),_ui[_idx])
  _dd[_n]=_rv;return _rv
 def _cnt_orb(_va,_vb,_R2):
  _eks=[]
  for _z in _pl:
   _sa=sum(_va[j]*(_z**j) for j in range(8));_sb=sum(_vb[j]*(_z**j) for j in range(8))
   _eks.append(2.0*(_sa.real*_sa.real+_sa.imag*_sa.imag+_sb.real*_sb.real+_sb.imag*_sb.imag))
  _bd=15.001*_R2;_Bk=[0.5*math.log(_bd/_ek) for _ek in _eks];_SB=sum(_Bk)
  _lo=[-(_SB-_Bk[k]) for k in range(3)]
  _nlo=[math.floor(sum(_Mi[i][k]*(_lo[k] if _Mi[i][k]>0 else _Bk[k]) for k in range(3))-1e-6) for i in range(3)]
  _nhi=[math.ceil(sum(_Mi[i][k]*(_Bk[k] if _Mi[i][k]>0 else _lo[k]) for k in range(3))+1e-6) for i in range(3)]
  _mo=0
  for _n0 in range(_nlo[0],_nhi[0]+1):
   for _n1 in range(_nlo[1],_nhi[1]+1):
    for _n2 in range(_nlo[2],_nhi[2]+1):
     _ok=True;_se=0.0
     for _k in range(4):
      _xk=_n0*_Lg[0][_k]+_n1*_Lg[1][_k]+_n2*_Lg[2][_k]
      if _xk>_Bk[_k]+1e-6:_ok=False;break
      _se+=_eks[_k]*math.exp(2.0*_xk)
     if not _ok or _se>_bd*1.0001:continue
     _u=mul_phi15_z(mul_phi15_z(_upow(0,_n0),_upow(1,_n1)),_upow(2,_n2))
     a0,a1,a2,a3,a4,a5,a6,a7=mul_phi15_z(_u,_va);b0,b1,b2,b3,b4,b5,b6,b7=mul_phi15_z(_u,_vb)
     for _mt in _ym:
      _s=0
      for r0,r1,r2,r3,r4,r5,r6,r7 in _mt:
       _xa=r0*a0+r1*a1+r2*a2+r3*a3+r4*a4+r5*a5+r6*a6+r7*a7;_s+=_xa*_xa
       if _s>_R2:break
       _xb=r0*b0+r1*b1+r2*b2+r3*b3+r4*b4+r5*b5+r6*b6+r7*b7;_s+=_xb*_xb
       if _s>_R2:break
      else:_mo+=1
  assert _mo>=1;return _mo
 def _is_cp(_va,_vb):
  _g0=math.gcd(abs(norm_q15_z(_va)),abs(norm_q15_z(_vb)))
  if _g0<=1:return True
  _t=_g0;_d=2;_ps=[]
  while _d*_d<=_t and _d<=10000:
   if _t%_d==0:
    _ps.append(_d)
    while _t%_d==0:_t//=_d
   _d+=1 if _d==2 else 2
  if _t>1:_ps.append(_t)
  for _pp in _ps:
   _g1=pgcd_mod(_va,PHI15,_pp)
   if len(_g1)>1 and len(pgcd_mod(_g1,_vb,_pp))>1:return False
  return True
 _rng=random.Random(20261120);_im_s=0.0;_r2_tot=0
 for _Rr in (32,64,112):
  _R2=_Rr*_Rr
  for _ in range(30):
   while True:
    _g=[_rng.gauss(0,1) for _ in range(16)];_rn=math.sqrt(sum(x*x for x in _g))
    _rd=_Rr*(_rng.random()**(1.0/16.0))/_rn;_vv=[int(round(x*_rd)) for x in _g]
    if 0<sum(x*x for x in _vv)<=_R2:break
   _r2_tot+=1
   if _is_cp(_vv[:8],_vv[8:]):_im_s+=1.0/_cnt_orb(_vv[:8],_vv[8:],_R2)
 _udup_live=_im_s/_r2_tot;assert _r2_tot==90 and _udup_live>float(desc["params"]["udup_node"])==0.40
 _p_ds=31;_c_ds=[1]+[_rng.randrange(_p_ds) for _ in range(15)]
 _B_ds=[[0]*16 for _ in range(16)];_B_ds[0][0]=_p_ds
 for _j in range(1,16):_B_ds[_j][0]=(-_c_ds[_j])%_p_ds;_B_ds[_j][_j]=1
 _Br_ds,_dg_ds,_lm_ds=lll_16d_int(_B_ds,3,4);assert _dg_ds[16]==_p_ds*_p_ds
 _bs2=[_dg_ds[i+1]/_dg_ds[i] for i in range(16)]
 _mu_ds=[[1.0 if i==j else (_lm_ds[i][j]/_dg_ds[j+1] if j<i else 0.0) for j in range(16)] for i in range(16)]
 _fp_h=[];_z_ds=[0]*16
 def _fp_rec(_k,_rm,_nz):
  _cc=-sum(_mu_ds[j][_k]*_z_ds[j] for j in range(_k+1,16));_rd=math.sqrt(max(0.0,_rm/_bs2[_k]))
  _lo=math.ceil(_cc-_rd);_hi=math.floor(_cc+_rd)
  if not _nz:_lo=max(0,_lo)
  for _zk in range(_lo,_hi+1):
   _z_ds[_k]=_zk;_nz2=_nz or (_zk>0)
   if _k==0:
    if not _nz2:continue
    _v=[sum(_z_ds[i]*_Br_ds[i][t] for i in range(16)) for t in range(16)]
    if 0<sum(x*x for x in _v)<=4:
     for _x in _v:
      if _x!=0:
       if _x<0:_v=[-_t for _t in _v]
       break
     _fp_h.append(tuple(_v))
   else:
    _nr=_rm-_bs2[_k]*(_zk-_cc)*(_zk-_cc)
    if _nr>=-1e-9:_fp_rec(_k-1,max(0.0,_nr),_nz2)
 _fp_rec(15,4.0+1e-9,False)
 _bf_h=[];_vb_ds=[0]*16
 def _bf_rec(_idx,_rm,_nz):
  if _idx==16:
   if _nz and sum(_c_ds[j]*_vb_ds[j] for j in range(16))%_p_ds==0:_bf_h.append(tuple(_vb_ds))
   return
  _r=math.isqrt(_rm)
  for _x in range(-_r,_r+1):
   if not _nz and _x<0:continue
   _vb_ds[_idx]=_x;_bf_rec(_idx+1,_rm-_x*_x,_nz or (_x>0))
 _bf_rec(0,4,False)
 assert len(_fp_h)==len(set(_fp_h))==len(_bf_h)>400 and set(_fp_h)==set(_bf_h)
 _bc=[[0]*5 for _ in range(17)];_bc[0]=[1]*5
 for _d in range(1,17):
  for _s in range(5):_bc[_d][_s]=sum(_bc[_d-1][_s-_x*_x] for _x in range(-math.isqrt(_s),math.isqrt(_s)+1))
 _s16_2=(_bc[16][4]-1)//2;_off2=(_bc[16][4]+1)//2;_rks=set()
 for _hv in _fp_h:
  _rk=0;_rm=4
  for _d in range(16):
   for _x in range(-math.isqrt(_rm),_hv[_d]):_rk+=_bc[15-_d][_rm-_x*_x]
   _rm-=_hv[_d]*_hv[_d]
  _idx=_rk-_off2;assert 0<=_idx<_s16_2;_rks.add(_idx)
 assert len(_rks)==len(_fp_h)
 _ds_nph={};_vb16=lambda _d,_R:(math.pi**(_d/2.0)/math.gamma(_d/2.0+1.0))*(_R**_d)
 for _bp in (16,48,70):
  _pb=(1<<(_bp-1))|(_rng.getrandbits(_bp-2)<<1)|1;_cb=[1]+[_rng.randrange(_pb) for _ in range(15)]
  _Bb=[[0]*16 for _ in range(16)];_Bb[0][0]=_pb
  for _j in range(1,16):_Bb[_j][0]=(-_cb[_j])%_pb;_Bb[_j][_j]=1
  _,_dgb,_=lll_16d_int(_Bb,3,4);assert _dgb[16]==_pb*_pb
  _ls=[math.sqrt(_dgb[i+1]/_dgb[i]) for i in range(16)]
  _rt=sum((_vb16(16-k,64)/_vb16(16,64))*math.prod(_ls[:k]) for k in range(1,16))
  assert (_rt*60.0+205.0)<550.0;_ds_nph[f"p2_{_bp}"]=round(_rt,3)
 _w_el=y_pow_z(5);_wp1=add_phi15_z(_w_el,[1]+[0]*7);_w2=y_pow_z(10)
 for _ in range(20):
  _aa=[_rng.randint(-50,50) for _ in range(8)];_bb=[_rng.randint(-50,50) for _ in range(8)]
  if not any(_aa) or not any(_bb):continue
  _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
  _ng_ex=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3)))))
  _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3);_a9,_b9=mul_phi15_z(_a6,_a3),mul_phi15_z(_b6,_b3)
  _nf_ex=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(_a9,mul_phi15_z(_wp1,mul_phi15_z(_a6,_b3))),mul_phi15_z(_a3,_b6)),mul_phi15_z(_w_el,_b9))))
  _bg_e=_bf_e=1.0
  for _z in _pl:
   _sa=sum(_aa[j]*(_z**j) for j in range(8));_sb=sum(_bb[j]*(_z**j) for j in range(8))
   _sa3=_sa**3;_sb3=_sb**3;_wz=_z**5;_w2z=_wz*_wz
   _bg_e*=abs(_sa3-float(U_bls)*_w2z*_sb3)**2
   _bf_e*=abs(_sa3**3+(_wz+1.0)*(_sa3**2)*_sb3+_sa3*(_sb3**2)-_wz*(_sb3**3))**2
  assert math.log2(_bg_e)+1.0>=math.log2(_ng_ex) and math.log2(_bf_e)+1.0>=math.log2(_nf_ex)
 _yr_min=1e99;_ss_yr=dlive["_ss_eval"]
 for _yi in (0,5,12,17):
  _sd,_qh,_ryh,_rth=data["prod_sq24"][_yi];_sn="g" if _sd==0 else "f"
  _qq,_ry,_rt=int(_qh,16),int(_ryh,16),int(_rth,16);_nb=_qq.bit_length()-1
  _prw=next(r for r in dlive["rows"] if r["side"]==_sn and r["lo"]<=_nb<r["hi"])
  _Rp=_prw["R"];_Rp2=_Rp*_Rp;_Bq=[[0]*16 for _ in range(16)];_Bq[0][0]=_qq;_rp=[pow(_ry,i,_qq) for i in range(8)]
  for _i in range(1,8):_Bq[_i][0]=(-_rp[_i])%_qq;_Bq[_i][_i]=1
  for _i in range(8):_Bq[8+_i][0]=(_rt*_rp[_i])%_qq;_Bq[8+_i][8+_i]=1
  _Rr,_,_=lll_16d_int(_Bq,99,100);_a1,_b1=_Rr[0][:8],_Rr[0][8:]
  _a2,_b2=next((v[:8],v[8:]) for v in _Rr[1:] if any(sub_phi15_z(mul_phi15_z(_a1,v[8:]),mul_phi15_z(_b1,v[:8]))))
  _nfs,_ngs,_wts=[],[],[]
  for _ in range(24):
   while True:
    _g=[_rng.gauss(0,1) for _ in range(16)];_rn=math.sqrt(sum(x*x for x in _g))
    _rd=_Rp*(_rng.random()**(1.0/16.0))/_rn;_vv=[int(round(x*_rd)) for x in _g]
    if 0<sum(x*x for x in _vv)<=_Rp2:break
   _al,_bt=_vv[:8],_vv[8:]
   if not _is_cp(_al,_bt):_wts.append(0.0);_nfs.append(9999.0);_ngs.append(9999.0);continue
   _wts.append(1.0/_cnt_orb(_al,_bt,_Rp2))
   _aa=add_phi15_z(mul_phi15_z(_al,_a1),mul_phi15_z(_bt,_a2));_bb=add_phi15_z(mul_phi15_z(_al,_b1),mul_phi15_z(_bt,_b2))
   _a3,_b3=pow_phi15_z(_aa,3),pow_phi15_z(_bb,3)
   _ng=abs(norm_q15_z(sub_phi15_z(_a3,scale_phi15_z(U_bls,mul_phi15_z(_w2,_b3)))))
   _a6,_b6=mul_phi15_z(_a3,_a3),mul_phi15_z(_b3,_b3);_a9,_b9=mul_phi15_z(_a6,_a3),mul_phi15_z(_b6,_b3)
   _nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(_a9,mul_phi15_z(_wp1,mul_phi15_z(_a6,_b3))),mul_phi15_z(_a3,_b6)),mul_phi15_z(_w_el,_b9))))
   _nfs.append(math.log2(_nf)-(math.log2(_qq) if _sn=="f" else 0.0));_ngs.append(math.log2(_ng)-(math.log2(_qq) if _sn=="g" else 0.0))
  _pf=_ss_yr(_nfs,float(_prw["n"]),_prw["iA"],m_node,float(desc["params"]["alpha_f"]),_prw["kA"],_prw["kB"])
  _pg=_ss_yr(_ngs,float(_prw["n"]),_prw["iA"],m_node,float(desc["params"]["alpha_g"]),_prw["kA"],_prw["kB"])
  _yr=(sum(w*a*b for w,a,b in zip(_wts,_pf,_pg))/len(_wts))/_prw["pr"]
  if _yr<_yr_min:_yr_min=_yr
 assert _yr_min>=10.0*float(desc["params"]["udup_node"])
 return round(_udup_live,4),round(_yr_min,2),len(_fp_h),_ds_nph

def verify_prod_cond_q15(enc12,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn,lll_16d_int):
 assert len(enc12)==12
 w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10)
 log2_u8=float(desc["params"]["log2_u8"]);S8=math.log2(desc["cells_rank2"]["8"])
 qb_f=[x[1] for x in desc["qbins"]];qb_g=[x[2] for x in desc["qbins"]]
 rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<16:
  al=[rng.randint(-4,4) for _ in range(8)];bt=[rng.randint(-4,4) for _ in range(8)]
  if any(al) and sum(x*x for x in al)+sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
 cfs,cgs,pr_rats=[],[],[]
 _ss=dlive["_ss_eval"];_mn=dlive["mesh_node"];_af=float(desc["params"]["alpha_f"]);_ag=float(desc["params"]["alpha_g"])
 for sid,q_hex,pa,pb,ry_hex,rt_hex in enc12:
  q,ry,rt=int(q_hex,16),int(ry_hex,16),int(rt_hex,16)
  pa3,pb3=pow_fn(pa,3),pow_fn(pb,3)
  if sid==0:pn=abs(norm_fn(sub_fn(pa3,scale_fn(U_bls,mul_fn(w2_el,pb3)))))
  else:
   pa6,pb6=mul_fn(pa3,pa3),mul_fn(pb3,pb3);pa9,pb9=mul_fn(pa6,pa3),mul_fn(pb6,pb3)
   pn=abs(norm_fn(sub_fn(add_fn(add_fn(pa9,mul_fn(wp1_el,mul_fn(pa6,pb3))),mul_fn(pa3,pb6)),mul_fn(w_el,pb9))))
  assert pn>0 and pn%q==0
  ry_p=[pow(ry,i,q) for i in range(8)]
  assert (sum(pa[i]*ry_p[i] for i in range(8))-rt*sum(pb[i]*ry_p[i] for i in range(8)))%q==0
  B=[[0]*16 for _ in range(16)];B[0][0]=q
  for i in range(1,8):B[i][0]=(-ry_p[i])%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*ry_p[i])%q;B[8+i][8+i]=1
  R,dsq,lsq=lll_16d_int(B,99,100);assert dsq[16]==q*q
  for kk in range(1,16):assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
  a1,b1=R[0][:8],R[0][8:]
  a2,b2=next((v[:8],v[8:]) for v in R[1:] if any(sub_fn(mul_fn(a1,v[8:]),mul_fn(b1,v[:8]))))
  sname="g" if sid==0 else "f";nb=q.bit_length()-1;lq=math.log2(q)
  prow=next(r for r in dlive["rows"] if r["side"]==sname and r["lo"]<=nb<r["hi"])
  S_bin=math.log2(desc["cells_rank2"][str(prow["R"])]);c_act=(lq+S_bin)/16.0
  c_val=(lq+S8)/16.0;nfas,ngas=[],[]
  for al,bt in ab_pts:
   a=add_fn(mul_fn(al,a1),mul_fn(bt,a2));b=add_fn(mul_fn(al,b1),mul_fn(bt,b2))
   a3,b3=pow_fn(a,3),pow_fn(b,3)
   ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
   a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
   nf=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))))
   assert (ng if sid==0 else nf)%q==0
   cf=math.log2(nf)-72.0*c_val;cg=math.log2(ng)-24.0*c_val-log2_u8
   cfs.append(cf);cgs.append(cg)
   nfas.append(72.0*c_act+cf-(lq if sname=="f" else 0.0));ngas.append(24.0*c_act+log2_u8+cg-(lq if sname=="g" else 0.0))
  _pfe=_ss(nfas,float(prow["n"]),prow["iA"],_mn,_af,prow["kA"],prow["kB"])
  _pge=_ss(ngas,float(prow["n"]),prow["iA"],_mn,_ag,prow["kA"],prow["kB"])
  pr_rats.append((sum(x*y for x,y in zip(_pfe,_pge))/len(ab_pts))/prow["pr"])
 cfs.sort();cgs.sort();n=len(cfs)
 sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
 sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
 assert all(sf6[i]<qb_f[i] for i in range(6)) and all(sg6[i]<qb_g[i] for i in range(6)) and min(pr_rats)>=10.0
 _prs=sorted(pr_rats)
 return dict(n_enc_live=12,parent_divisibility_live="12/12",live_bareiss_99_100_all12=True,
  presupplied_vectors=0,enc_sext_f=sf6,enc_sext_g=sg6,
  min_headroom_f=round(min(qb_f[i]-sf6[i] for i in range(6)),2),
  min_headroom_g=round(min(qb_g[i]-sg6[i] for i in range(6)),2),
  joint_pr_over_charged_min=round(_prs[0],2),joint_pr_over_charged_med=round(_prs[len(_prs)//2],2),
  all_sextiles_below=True,checks_ok=True)
