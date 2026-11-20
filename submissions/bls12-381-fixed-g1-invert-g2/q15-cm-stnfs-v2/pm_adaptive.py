"""Self-contained classic PM-Basis live verifier & Q(zeta_15) real-2-core-FB adaptive H5 verifier."""
import base64,json,lzma,math

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
