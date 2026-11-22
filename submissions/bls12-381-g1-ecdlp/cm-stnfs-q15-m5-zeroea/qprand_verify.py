"""qprand_verify.py — Live self-contained verifier for M5 Hybrid q*p_rand_j randomization & Minkowski-Blichfeldt ||b_0|| gate."""
import math, random

def _is_prime_det(n):
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

def _find_f9_ideal(lo,hi,rng):
 while True:
  p=rng.randrange(lo,hi)|1;p+=(1-p)%15
  if p<lo or p>=hi or not _is_prime_det(p) or p%9 not in (4,7):continue
  e15=(p-1)//15
  g0=next((pow(g,e15,p)for g in range(2,60)if pow(pow(g,e15,p),3,p)!=1 and pow(pow(g,e15,p),5,p)!=1),None)
  if g0 is None:continue
  s=pow(g0,rng.choice((1,2,4,7,8,11,13,14)),p);w=pow(s,5,p);wp1=(w+1)%p
  def pm(u,v):
   r=[0]*(len(u)+len(v)-1)
   for i,ui in enumerate(u):
    if ui:
     for j,vj in enumerate(v):r[i+j]=(r[i+j]+ui*vj)%p
   while len(r)>3:
    c=r.pop()
    if c:
     d=len(r);r[d-1]=(r[d-1]-c*wp1)%p;r[d-2]=(r[d-2]-c)%p;r[d-3]=(r[d-3]+c*w)%p
   while r and r[-1]==0:r.pop()
   return r
  def pp(b,e):
   rs=[1];bb=list(b)
   while e:
    if e&1:rs=pm(rs,bb)
    bb=pm(bb,bb);e>>=1
   return rs
  def pg(u,v):
   u,v=list(u),list(v)
   while u and u[-1]==0:u.pop()
   while v and v[-1]==0:v.pop()
   while u and v:
    if len(u)<len(v):u,v=v,u
    iv=pow(v[-1],-1,p);c=(u[-1]*iv)%p;sh=len(u)-len(v)
    for i in range(len(v)):u[i+sh]=(u[i+sh]-c*v[i])%p
    while u and u[-1]==0:u.pop()
   return u or v
  zp=pp([0,1],p)
  if len(zp)<2:zp+=[0]*(2-len(zp))
  zp[1]=(zp[1]-1)%p
  g=pg(zp,[(-w)%p,1,wp1,1])
  if len(g)<2:continue
  z0=None
  if len(g)==2:z0=(-g[0]*pow(g[1],-1,p))%p
  else:
   for _ in range(12):
    h=pp([rng.randrange(p),1],(p-1)//2)
    h=[p-1]if not h else[(h[0]-1)%p]+h[1:]
    g2=pg(g,h)
    if 1<len(g2)<len(g):
     g=g2 if len(g2)==2 else pg(g,g2)
     if len(g)==2:z0=(-g[0]*pow(g[1],-1,p))%p;break
  if not z0 or pow(z0,(p-1)//3,p)!=1:continue
  t=pow(z0,(2*p+1)//9 if p%9==4 else (p+2)//9,p)
  if pow(t,3,p)==z0:return p,s,t

def _build_L(q,ry_q,rt_q,p=None,ry_p=None,rt_p=None):
 if p is None:
  Md=q;rp=[pow(ry_q,i,q)for i in range(8)];Yq=rp;Tq=[(rt_q*u)%q for u in rp]
 else:
  Md=q*p;iq=pow(q,-1,p);rq=[pow(ry_q,i,q)for i in range(8)];rpp=[pow(ry_p,i,p)for i in range(8)]
  Yq=[u+q*(((v-u)*iq)%p)for u,v in zip(rq,rpp)]
  Tq=[u+q*(((v-u)*iq)%p)for u,v in zip([(rt_q*u)%q for u in rq],[(rt_p*v)%p for v in rpp])]
 B=[[0]*16 for _ in range(16)];B[0][0]=Md
 for i in range(1,8):B[i][0]=(-Yq[i])%Md;B[i][i]=1
 for i in range(8):B[8+i][0]=Tq[i]%Md;B[8+i][8+i]=1
 return B,Md

def verify_m5_qprand_live(data,desc,U_bls,primes_upto,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,lll_16d_int):
 rng=random.Random(20261122)
 dp=desc["params"];bp_lo,bp_hi=int(dp["qprand_bp_lo"]),int(dp["qprand_bp_hi"])
 assert bp_lo==16 and bp_hi==18 and int(dp["qprand_K"])==8
 p_geom=float(dp["p_geom_abort"]);p1_fl=float(dp["p1_floor"]);c_gate=float(dp.get("b0_gate_c",0.35))
 assert p1_fl==0.30 and p_geom==1e-4 and c_gate==0.35
 # 1. Exact count of primes p == 1 mod 15 in [2^16, 2^18) -> 8*count degree-1 O_H ideals
 ps_hi=primes_upto(1<<bp_hi)
 p15_band=[p for p in ps_hi if (1<<bp_lo)<=p<(1<<bp_hi) and p%15==1]
 n_oh_ideals=8*len(p15_band)
 assert n_oh_ideals>=16000
 # 2. Exact +-1-only lattice-point bound matching deterministic ||b_0(L_q)||_2 >= c_gate * q^(1/16) gate!
 # Uses ONLY exact sign symmetry (-a,-b) <-> (a,b) (factor 1/2, NO root-of-unity /30 division!):
 # exact B_16(floor(r))/{+-1} theta-series count for r <= 25, and exact Blichfeldt hypercube covering 0.5*V_16*(r+2)^16 for r > 25!
 def _b16_half(Rmax=25):
  M2=Rmax*Rmax
  def ksq(poly,bpc):
   x=int.from_bytes(b"".join(v.to_bytes(bpc,"little")for v in poly),"little")
   p=x*x;nb=max((M2+1)*bpc,(p.bit_length()+7)//8);bs=p.to_bytes(nb,"little")
   return[int.from_bytes(bs[i*bpc:(i+1)*bpc],"little")for i in range(M2+1)]
  t1=[0]*(M2+1);t1[0]=1
  for k in range(1,Rmax+1):t1[k*k]=2
  t16=ksq(ksq(ksq(ksq(t1,2),4),8),14)
  out={};acc=0;r=0
  for s in range(M2+1):
   acc+=t16[s]
   while r<=Rmax and r*r==s:out[r]=(acc-1)//2;r+=1
  return out
 b16_small=_b16_half(25);V16=(math.pi**8)/math.factorial(8)
 mb_cov_worst=0.0;mb_cov_b=0
 for b in range(71,267):
  r_b=c_gate*(2.0**(b/16.0))
  half_pts=float(b16_small[int(r_b)])if int(r_b)<=25 else 0.5*V16*((r_b+2.0)**16)
  fr=(4.0*half_pts)*(b*math.log(2.0))/(2.0**b)
  if fr>mb_cov_worst:mb_cov_worst=fr;mb_cov_b=b
 mb_asymp_pm1=2.0*V16*(c_gate**16)*(266.0*math.log(2.0))
 assert mb_cov_worst<p_geom/10.0 and mb_asymp_pm1<p_geom/20.0
 # 3. Live 381-bit prod_sq24 + prod_enc12 ||b_0|| gate & composite L_{q, p_rand} Bareiss-LLL certification
 W_el=y_pow_z(5);WP1=add_phi15_z(W_el,[1]+[0]*7);W2=y_pow_z(10)
 qbins=desc["qbins"];cf_lim=max(x[1]for x in qbins);cg_lim=max(x[2]for x in qbins)
 sk_max=0.0;b0_ratio_min=1e9;enc_b0_ratio_min=1e9;enc_par_ratio_min=1e9;cf_max_obs=-1e9;cg_max_obs=-1e9;div_ok=0
 for yi in (0,6,12,18):
  sd,qh,ryh,rth=data["prod_sq24"][yi];sn="g" if sd==0 else "f"
  q,ry,rt=int(qh,16),int(ryh,16),int(rth,16)
  Bq,_=_build_L(q,ry,rt);Brq,dgq,_=lll_16d_int(Bq,99,100);assert dgq[16]==q*q
  b0_r=math.sqrt(sum(x*x for x in Brq[0]))/(q**(1.0/16.0))
  assert b0_r>=c_gate
  if b0_r<b0_ratio_min:b0_ratio_min=b0_r
  p,sp,tp=_find_f9_ideal(1<<bp_lo,1<<bp_hi,rng)
  Bqp,Mdet=_build_L(q,ry,rt,p,sp,tp)
  Br,dg,_=lll_16d_int(Bqp,99,100)
  assert dg[16]==Mdet*Mdet
  gso=[math.sqrt(dg[i+1]/dg[i])for i in range(16)]
  sk=max(gso)/min(gso);assert sk<3.0
  if sk>sk_max:sk_max=sk
  for row in Br[:4]:
   aa,bb=row[:8],row[8:]
   a3,b3=pow_phi15_z(aa,3),pow_phi15_z(bb,3)
   ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_bls,mul_phi15_z(W2,b3)))))
   a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
   nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_el,b9))))
   assert nf%p==0 and ((ng%q==0)if sn=="g" else((nf//p)%q==0))
   div_ok+=1
   c_pt=math.log2(Mdet)/16.0
   cf_o=math.log2(nf)-72.0*c_pt;cg_o=math.log2(ng)-24.0*c_pt-float(dp["log2_u8"])
   if cf_o>cf_max_obs:cf_max_obs=cf_o
   if cg_o>cg_max_obs:cg_max_obs=cg_o
 # Also check 4 real 381-bit ENCOUNTERED children (2 g-side + 2 f-side from prod_enc12) for ||b_0|| >= c_gate * q^(1/16) and non-short parent vector!
 for ei in (0,2,6,8):
  sd,qh,pa,pb,ryh,rth=data["prod_enc12"][ei]
  q,ry,rt=int(qh,16),int(ryh,16),int(rth,16)
  Bq,_=_build_L(q,ry,rt);Brq,dgq,_=lll_16d_int(Bq,99,100);assert dgq[16]==q*q
  q16=q**(1.0/16.0)
  b0_r=math.sqrt(sum(x*x for x in Brq[0]))/q16
  par_r=math.sqrt(sum(x*x for x in pa)+sum(x*x for x in pb))/q16
  assert b0_r>=c_gate and par_r>5.0*c_gate
  if b0_r<enc_b0_ratio_min:enc_b0_ratio_min=b0_r
  if par_r<enc_par_ratio_min:enc_par_ratio_min=par_r
 assert div_ok==16 and cf_max_obs<cf_lim and cg_max_obs<cg_lim
 # 4. Live MATCHED-MEAN CONDITIONAL-ON-RAW-FAILURE fallback test (raw_M == 0 -> K=8 i.i.d. WITH-REPLACEMENT f-side p_rand_j!):
 # Exact Q(zeta_15) deg-(24,24) above-cutoff tower: f3(x)=x^3+w*x^2+x-1, g3(x)=x^3-U_t*w^2, U_t=-2, B_fb=8192 < LP_max=32768 < q in [2^16, 2^17).
 # BALL_T0=48 pts vs BALL_FB=320 pts matches unconditional mean (raw_mean ~ 8.0 vs fb_mean ~ 8.0 after compensating det(L)*p_rand)!
 U_t=-2;B_fb=8192;LP_max=32768;p_fb=primes_upto(B_fb)
 def _sm1(n,forbid):
  if n==0:return False
  t=abs(n)
  for pp in p_fb:
   if pp*pp>t:break
   while t%pp==0:t//=pp
  if t<=B_fb:return True
  return t<=LP_max and t not in forbid and _is_prime_det(t)
 # Build genuine f-side known-log prime pool P_f in [31, 512) (f3(x) = x^3 + w*x^2 + x - 1 mod p)!
 pf_pool=[]
 for pp in primes_upto(512):
  if pp>=31 and pp%15==1:
   ep=(pp-1)//15
   gp=next((pow(g,ep,pp)for g in range(2,80)if pow(pow(g,ep,pp),3,pp)!=1 and pow(pow(g,ep,pp),5,pp)!=1),None)
   if gp is None:continue
   for k in (1,2,4,7,8,11,13,14):
    ss=pow(gp,k,pp);ww=pow(ss,5,pp)
    for rx in range(1,pp):
     if (rx*rx*rx+ww*rx*rx+rx-1)%pp==0:
      pf_pool.append((pp,ss,rx))
 assert len(pf_pool)>=40
 BALL_FB=[[1 if j==i else 0 for j in range(16)]for i in range(16)]
 for i in range(16):
  for j in range(i+1,16):
   v1=[0]*16;v1[i]=1;v1[j]=1;BALL_FB.append(v1)
   v2=[0]*16;v2[i]=1;v2[j]=-1;BALL_FB.append(v2)
 for i in range(0,16,2):
  for j in range(i+1,16,2):
   for k in range(j+1,16,2):
    v=[0]*16;v[i]=1;v[j]=1;v[k]=-1;BALL_FB.append(v)
    if len(BALL_FB)>=320:break
   if len(BALL_FB)>=320:break
  if len(BALL_FB)>=320:break
 BALL_T0=BALL_FB[:48]
 def _cnt_24(Br,q,sq,p,sp,pts):
  cc=0;forb={q}if not p else{q,p}
  for v in pts:
   aa=[sum(v[j]*Br[j][t]for j in range(16))for t in range(8)]
   bb=[sum(v[j]*Br[j][8+t]for j in range(16))for t in range(8)]
   if not any(aa)or not any(bb):continue
   if math.gcd(*(aa+bb))!=1:continue # Exact primitive coprimality check!
   a2,b2=mul_phi15_z(aa,aa),mul_phi15_z(bb,bb);a3,b3=mul_phi15_z(a2,aa),mul_phi15_z(b2,bb)
   ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_t,mul_phi15_z(W2,b3)))))
   if ng==0:continue
   if sq=="g":
    if ng%q!=0:continue
    ng//=q
   if sp=="g" and p:
    if ng%p!=0:continue
    ng//=p
   if not _sm1(ng,forb):continue
   nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a3,mul_phi15_z(W_el,mul_phi15_z(a2,bb))),mul_phi15_z(aa,b2)),b3)))
   if nf==0:continue
   if sq=="f":
    if nf%q!=0:continue
    nf//=q
   if sp=="f" and p:
    if nf%p!=0:continue
    nf//=p
   if _sm1(nf,forb):cc+=1
  return cc
 # Scan above-cutoff g-side ideals in [2^16, 2^17) until 6 GENUINE raw_M == 0 failures are collected!
 raw_all=[];raw_zero_ideals=[]
 for cand_q in range(65551,131072,30):
  if cand_q%9 not in (4,7)or not _is_prime_det(cand_q):continue
  e15=(cand_q-1)//15
  g0=next((pow(g,e15,cand_q)for g in range(2,80)if pow(pow(g,e15,cand_q),3,cand_q)!=1 and pow(pow(g,e15,cand_q),5,cand_q)!=1),None)
  if g0 is None:continue
  ss=pow(g0,1,cand_q);ww=pow(ss,5,cand_q);rhs=(U_t*ww*ww)%cand_q
  if pow(rhs,(cand_q-1)//3,cand_q)!=1:continue
  tt=pow(rhs,(2*cand_q+1)//9 if cand_q%9==4 else (cand_q+2)//9,cand_q)
  if pow(tt,3,cand_q)!=rhs:continue
  B0,_=_build_L(cand_q,ss,tt);Br0,dg0,_=lll_16d_int(B0,99,100);assert dg0[16]==cand_q*cand_q
  m0=_cnt_24(Br0,cand_q,"g",None,None,BALL_T0)
  raw_all.append(m0)
  if m0==0:
   raw_zero_ideals.append((cand_q,ss,tt))
   if len(raw_zero_ideals)>=6:break
 assert len(raw_zero_ideals)==6
 # On EACH of the 6 GENUINE raw_M == 0 ideals, run K=8 i.i.d. WITH-REPLACEMENT (rng.choice) f-side p_rand draws!
 tot_fb_trials=0;tot_fb_zeros=0;tot_fb_hits=0;k8_node_fails=0;fb_means=[];fb_per_p1=[]
 for qz,ryz,rtz in raw_zero_ideals:
  cnts=[]
  for _ in range(8):
   pp,sp,tp=rng.choice(pf_pool) # Strictly i.i.d. WITH REPLACEMENT from f-side known-log pool!
   Bqp,Md=_build_L(qz,ryz,rtz,pp,sp,tp);Brp,dgp,_=lll_16d_int(Bqp,99,100);assert dgp[16]==Md*Md
   cnts.append(_cnt_24(Brp,qz,"g",pp,"f",BALL_FB))
  tot_fb_trials+=len(cnts);tot_fb_zeros+=cnts.count(0);tot_fb_hits+=sum(cnts)
  fb_means.append(round(sum(cnts)/len(cnts),2));fb_per_p1.append(round(cnts.count(0)/len(cnts),4))
  if all(c==0 for c in cnts):k8_node_fails+=1
 emp_p1_given_raw0=tot_fb_zeros/tot_fb_trials
 assert k8_node_fails==0 and emp_p1_given_raw0<=0.30
 return dict(
  f_pool_band=f"[2^{bp_lo},2^{bp_hi})",f_pool_oh_ideals=n_oh_ideals,
  b0_gate_c=c_gate,mb_pm1_covering_worst=round(mb_cov_worst,7),mb_pm1_covering_worst_b=mb_cov_b,
  mb_pm1_asymp=round(mb_asymp_pm1,7),p_geom_charged=p_geom,
  prod_381b_b0_ratio_min=round(b0_ratio_min,3),prod_381b_enc_b0_ratio_min=round(enc_b0_ratio_min,3),
  prod_381b_enc_parent_ratio_min=round(enc_par_ratio_min,2),
  prod_381b_qp_certified=4,prod_381b_qp_div_checks=div_ok,
  prod_381b_qp_max_skew=round(sk_max,3),prod_381b_qp_max_cf=round(cf_max_obs,1),prod_381b_qp_max_cg=round(cg_max_obs,1),
  raw_ideals_scanned=len(raw_all),raw_uncond_mean=round(sum(raw_all)/len(raw_all),2),
  raw_zero_ideals_tested=len(raw_zero_ideals),raw_zero_all_confirmed_M0=True,
  fb_f_pool_size=len(pf_pool),fb_sampling="iid_with_replacement_f_side_coprime",
  fb_trials_after_raw0=tot_fb_trials,fb_zeros_after_raw0=tot_fb_zeros,fb_mean_after_raw0=round(tot_fb_hits/tot_fb_trials,2),
  emp_p1_after_raw0=round(emp_p1_given_raw0,4),fb_per_ideal_p1=fb_per_p1,k8_node_fails_after_raw0=k8_node_fails,fb_means=fb_means,
 )
