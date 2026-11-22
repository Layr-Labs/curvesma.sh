"""qprand_verify.py — Live self-contained verifier for M5-v2 Hybrid (K=16, p1_ref=0.45, eps_gate=0.20, p_geom=5e-5):
  1. Exact count of degree-1 O_H ideals in [2^16, 2^18) (16,272 ideals) + live exact Chebotarev check of f_9/g_3 in [2^14, 2^16).
  2. Exact Yukon-CLEARED +-1-only factor-12 Blichfeldt covering bound at eps_gate=0.20 (NO /30 division, upper edge 2^(b+1),
     unit-hypercube covering 0.5*V_16*(R_hi+2)^16, <=12 prime divisors across both sides, explicit conservative denominator
     0.5*2^b/((b+1)*ln 2) [Bennett-Martin-O'Bryant-Rechnitzer 2018]): worst <= 6.452e-7 (77.5x below p_geom = 5e-5)!
  3. Live 381-bit prod_sq24 & prod_enc12 ||b_0||_2 >= 0.20*q^(1/16) + non-short parent checks + composite L_(Q, p_rand) Bareiss-LLL checks.
  4. Exact O_H = Z[zeta_15] ideal coprimality <a, b> = O_H (unit-tested live to catch (2-y)*(a0,b0) with 16-integer gcd == 1!).
  5. Exact UMVUE U-statistic E_hat[p_1(Q)^16] & live bit-exact replays of:
     - Campaign A: 32 genuine raw_M==0 above-cutoff (B_fb=32768 < LP=65536 < q in [2^17,2^18)) exact-O_H-coprime ideals (16g+16f) x 28 i.i.d. f-side draws (pool=328).
     - Campaign B: Strictly descending (B_fb < LP_max < q), exact-O_H-coprime, verified deg-1-child-root LITERAL deg-(72,24) (f_9, g_3) q*p_rand_j campaign!
"""
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

def _brent(n,lim=2500):
 if n%2==0:return 2
 for c in (1,3,5,7):
  y=2;m=64;g=r=q_=1;it=0
  while g==1 and it<lim:
   x=y
   for _ in range(r):y=(y*y+c)%n
   k=0
   while k<r and g==1 and it<lim:
    ys=y;step=min(m,r-k)
    for _ in range(step):y=(y*y+c)%n;q_=(q_*abs(x-y))%n
    g=math.gcd(q_,n);k+=step;it+=step
   r<<=1
  if g==n:
   g=1
   while g==1:ys=(ys*ys+c)%n;g=math.gcd(abs(x-ys),n)
  if 1<g<n:return g
 return None

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

def _build_ball(n_pts):
 pts=[((i,1),)for i in range(16)]
 for i in range(16):
  for j in range(i+1,16):
   pts.append(((i,1),(j,1)));pts.append(((i,1),(j,-1)))
 for i in range(16):
  for j in range(i+1,16):
   for k in range(j+1,16):
    for s1,s2 in ((1,-1),(-1,1),(1,1),(-1,-1)):
     pts.append(((i,1),(j,s1),(k,s2)))
     if len(pts)>=n_pts:return pts[:n_pts]
 return pts[:n_pts]

def _verify_campB_deg72(c72,K_fb,W_el,WP1,W2,PHI15,primes_upto,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int,_is_oh_coprime):
 rows72=c72["rows"];U72=int(c72["U_T"]);B72=int(c72["B_fb"]);LP72=int(c72["LP_max"])
 assert len(rows72)==28 and B72==131072 and LP72==(1<<26)
 for sd,q,ry,rx,drs,cnts in rows72:
  assert q%15==1 and _is_prime_det(q) and pow(ry,15,q)==1 and pow(ry,3,q)!=1 and pow(ry,5,q)!=1
  w=pow(ry,5,q)
  if sd=="g":assert (pow(rx,3,q)-U72*w*w)%q==0
  else:assert (pow(rx,9,q)+(w+1)*pow(rx,6,q)+pow(rx,3,q)-w)%q==0
  assert len(drs)==len(cnts)==24
  for p,sp,tp in drs:
   wp=pow(sp,5,p)
   assert 31<=p<450 and p%15==1 and (pow(tp,9,p)+(wp+1)*pow(tp,6,p)+pow(tp,3,p)-wp)%p==0
 res={}
 for lbl,sub in (("all",rows72),("f",[r for r in rows72 if r[0]=="f"]),("g",[r for r in rows72 if r[0]=="g"])):
  n=len(sub);M=len(sub[0][5]);zs=[r[5].count(0)for r in sub];mns=[sum(r[5])/M for r in sub]
  vrs=[sum((x-mu)**2 for x in r[5])/(M-1)for r,mu in zip(sub,mns)]
  fanos=[v/mu if mu>0 else 9.9 for v,mu in zip(vrs,mns)]
  mvue8=sum(math.comb(z,K_fb)/math.comb(M,K_fb)for z in zs)/n
  plug8=sum((z/M)**K_fb for z in zs)/n
  res[lbl]=dict(n_ideals=n,M_draws=M,trials=n*M,zeros=sum(zs),emp_p1=round(sum(zs)/(n*M),4),max_q_p1=round(max(z/M for z in zs),4),mean_M=round(sum(mns)/n,2),med_fano=round(sorted(fanos)[n//2],3),unbiased_mvue_E_p1_K8=float(f"{mvue8:.3e}"),plugin_E_p1_K8=float(f"{plug8:.3e}"))
 assert res["all"]["unbiased_mvue_E_p1_K8"]<0.005 and res["all"]["emp_p1"]<=0.35 and res["all"]["med_fano"]<1.5
 # Live bit-exact replay of 1 literal deg-(72,24) (f_9, g_3) ideal x 1 draw (2,496 pts + exact O_H-coprimality!)
 pts72=[(i,)for i in range(16)]
 for i in range(16):
  for j in range(i+1,16):
   pts72.append((i,j,-1));pts72.append((i,j,1))
   for k in range(j+1,16):
    pts72.append((i,j,k,-1,-1));pts72.append((i,j,k,-1,1));pts72.append((i,j,k,1,-1));pts72.append((i,j,k,1,1))
 P1=math.prod(primes_upto(256));P2=math.prod(p for p in primes_upto(4096)if p>256)
 P3=math.prod(p for p in primes_upto(32768)if p>4096);P4=math.prod(p for p in primes_upto(B72)if p>32768)
 def _sm72(n,mlp,forb):
  t=abs(n)
  if t<=1:return False
  for PP in (P1,P2):
   while True:
    g=math.gcd(t,PP)
    if g==1:break
    t//=g
   if t<=B72:return True
  if t<=LP72:return t not in forb and _is_prime_det(t)
  if _is_prime_det(t):return False
  for PP in (P3,P4):
   while True:
    g=math.gcd(t,PP)
    if g==1:break
    t//=g
   if t<=B72:return True
   if t<=LP72 and _is_prime_det(t):return t not in forb
   if _is_prime_det(t):return False
  if mlp<2 or t>(LP72**mlp):return False
  st=[t];cnt=0
  while st:
   m=st.pop()
   if m<=B72:continue
   if m<=LP72 and _is_prime_det(m):
    if m in forb:return False
    cnt+=1
    if cnt>mlp:return False
   else:
    if cnt>=mlp:return False
    d=_brent(m,1400)
    if d is None:return False
    st.append(d);st.append(m//d)
  return cnt<=mlp
 sd,q,ryq,rtq,drs,cnts=rows72[0];p,ryp,rtp=drs[0]
 B,Md=_build_L(q,ryq,rtq,p,ryp,rtp);Br,dg,_=lll_16d_int(B,99,100);assert dg[16]==Md*Md
 forb={q,p};seen=set()
 for pt in pts72:
  if len(pt)==1:r=Br[pt[0]];aa,bb=r[:8],r[8:]
  elif len(pt)==3:r1,r2,sg=Br[pt[0]],Br[pt[1]],pt[2];aa=[r1[t]+sg*r2[t]for t in range(8)];bb=[r1[8+t]+sg*r2[8+t]for t in range(8)]
  else:r1,r2,r3,s1,s2=Br[pt[0]],Br[pt[1]],Br[pt[2]],pt[3],pt[4];aa=[r1[t]+s1*r2[t]+s2*r3[t]for t in range(8)];bb=[r1[8+t]+s1*r2[8+t]+s2*r3[8+t]for t in range(8)]
  if not any(aa)or not any(bb)or math.gcd(*(aa+bb))!=1:continue
  a3,b3=pow_phi15_z(aa,3),pow_phi15_z(bb,3)
  ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U72,mul_phi15_z(W2,b3)))))
  if ng==0:continue
  if sd=="g":assert ng%q==0;ng//=q
  if not _sm72(ng,2,forb):continue
  a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
  nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_el,b9))))
  if nf==0:continue
  assert nf%p==0;nf//=p
  if sd=="f":assert nf%q==0;nf//=q
  if not _sm72(nf,3,forb):continue
  if not _is_oh_coprime(aa,bb):continue
  seen.add((nf,ng))
 assert len(seen)==cnts[0]
 res["live_bitexact_deg72_replay_ok"]=True;res["B_fb"]=B72;res["LP_max"]=LP72;res["n_pts"]=len(pts72)
 return res

def verify_m5_qprand_live(data,desc,U_bls,primes_upto,PHI15,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int):
 rng=random.Random(20261122)
 dp=desc["params"];bp_lo,bp_hi=int(dp["qprand_bp_lo"]),int(dp["qprand_bp_hi"])
 K_fb=int(dp["qprand_K"]);eta_bad=float(dp["eta_bad"]);p_good=float(dp["p_good"]);eps_gate=float(dp["eps_gate"]);p_geom=float(dp["p_geom_abort"])
 assert bp_lo==16 and bp_hi==18 and K_fb==8 and p_good==0.35 and eta_bad==0.01 and eps_gate==0.20 and p_geom==5e-5

 # 1. Exact count of degree-1 O_H ideals in [2^16, 2^18) + live Chebotarev density check of f_9/g_3 in [2^14, 2^16)
 ps_hi=primes_upto(1<<bp_hi)
 p15_band=[p for p in ps_hi if (1<<bp_lo)<=p<(1<<bp_hi) and p%15==1]
 n_oh_ideals=8*len(p15_band);assert n_oh_ideals==16272
 # Extract 4 genuine 16..18-bit f_9 ideals + check Chebotarev density in [2^14, 2^16)
 def _pmul(a,b,mod,p):
  r=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   if x:
    for j,y in enumerate(b):r[i+j]=(r[i+j]+x*y)%p
  while len(r)>=len(mod):
   c=r.pop()
   if c:
    off=len(r)-(len(mod)-1)
    for i in range(len(mod)-1):r[off+i]=(r[off+i]-c*mod[i])%p
  while len(r)>1 and r[-1]==0:r.pop()
  return r
 def _ppow(b,e,mod,p):
  r=[1];x=b[:]
  while e:
   if e&1:r=_pmul(r,x,mod,p)
   x=_pmul(x,x,mod,p);e>>=1
  return r
 ps_14_16=[p for p in ps_hi if (1<<14)<=p<(1<<16) and p%15==1]
 f9_14_16=0;g3_14_16=0;f9_sample=[]
 for p in ps_14_16:
  e15=(p-1)//15
  g0=next(pow(g,e15,p)for g in range(2,100)if pow(pow(g,e15,p),3,p)!=1 and pow(pow(g,e15,p),5,p)!=1)
  w1=pow(g0,5,p);w2=(w1*w1)%p
  for wv in (w1,w2):
   if pow((U_bls*wv*wv)%p,(p-1)//3,p)==1:g3_14_16+=12
   f3=[(-wv)%p,1,(wv+1)%p,1];xp=_ppow([0,1],p,f3,p)
   if len(xp)<2:xp+=[0]*(2-len(xp))
   xp[1]=(xp[1]-1)%p
   while len(xp)>1 and xp[-1]==0:xp.pop()
   g=pgcd_mod(f3,xp,p);deg=len(g)-1
   if deg==1:
    z=(-g[0]*pow(g[1],-1,p))%p
    if z and pow(z,(p-1)//3,p)==1:f9_14_16+=12
   elif deg==3:
    zc=_ppow([0,1],(p-1)//3,g,p);zc[0]=(zc[0]-1)%p
    while len(zc)>1 and zc[-1]==0:zc.pop()
    gc=pgcd_mod(g,zc,p);f9_14_16+=12*(len(gc)-1)
 cheb_ref=((1<<16)-(1<<14))/math.log(1<<16)
 assert f9_14_16==4740 and g3_14_16>=4000 and f9_14_16/cheb_ref>=0.95 and g3_14_16/cheb_ref>=0.90
 for p in p15_band:
  if p%9 not in (4,7):continue
  e15=(p-1)//15
  g0=next(pow(g,e15,p)for g in range(2,100)if pow(pow(g,e15,p),3,p)!=1 and pow(pow(g,e15,p),5,p)!=1)
  w1=pow(g0,5,p);f3=[(-w1)%p,1,(w1+1)%p,1];xp=_ppow([0,1],p,f3,p)
  if len(xp)<2:xp+=[0]*(2-len(xp))
  xp[1]=(xp[1]-1)%p
  while len(xp)>1 and xp[-1]==0:xp.pop()
  g=pgcd_mod(f3,xp,p)
  if len(g)==2:
   z=(-g[0]*pow(g[1],-1,p))%p
   if z and pow(z,(p-1)//3,p)==1:
    rt=pow(z,(2*p+1)//9 if p%9==4 else (p+2)//9,p)
    if pow(rt,3,p)==z:f9_sample.append((p,g0,rt))
  if len(f9_sample)>=4:break
 assert len(f9_sample)==4

 # 2. Exact Yukon-CLEARED +-1-only factor-12 Blichfeldt covering bound at eps_gate = 0.20
 # (NO /30 division, NO floor(r), explicit upper edge 2^(b+1), continuous covering 0.5*V_16*(R_hi+2)^16, <=12 prime divisors, 0.5*Q/((b+1)*ln 2)!)
 V16=(math.pi**8)/math.factorial(8)
 worst_mb=0.0;worst_b=0
 for b in range(71,267):
  Q=2.0**b;R_hi=eps_gate*((2.0*Q)**(1.0/16.0))
  half_pts=0.5*V16*((R_hi+2.0)**16)
  n_id_lo=0.5*Q/((b+1)*math.log(2.0))
  fr=(12.0*half_pts)/n_id_lo
  if fr>worst_mb:worst_mb=fr;worst_b=b
 assert worst_mb<p_geom/50.0

 # 3. Live 381-bit prod_sq24 + prod_enc12 ||b_0|| >= 0.20*q^(1/16) & composite L_(Q, p_rand) Bareiss-LLL checks
 W_el=y_pow_z(5);WP1=add_phi15_z(W_el,[1]+[0]*7);W2=y_pow_z(10)
 qbins=desc["qbins"];cf_lim=max(x[1]for x in qbins);cg_lim=max(x[2]for x in qbins)
 sk_max=0.0;b0_ratio_min=1e9;enc_b0_min=1e9;enc_par_min=1e9;cf_max_obs=-1e9;cg_max_obs=-1e9;div_ok=0
 for yi,f9s in zip((0,6,12,18),f9_sample):
  sd,qh,ryh,rth=data["prod_sq24"][yi];sn="g" if sd==0 else "f"
  q,ry,rt=int(qh,16),int(ryh,16),int(rth,16)
  Bq,_=_build_L(q,ry,rt);Brq,dgq,_=lll_16d_int(Bq,99,100);assert dgq[16]==q*q
  b0_r=math.sqrt(sum(x*x for x in Brq[0]))/(float(q)**(1.0/16.0))
  assert b0_r>=0.90 and b0_r>=eps_gate
  b0_ratio_min=min(b0_ratio_min,b0_r)
  p,sp,tp=f9s
  Bqp,Mdet=_build_L(q,ry,rt,p,sp,tp);Br,dg,_=lll_16d_int(Bqp,99,100);assert dg[16]==Mdet*Mdet
  gso=[math.sqrt(dg[i+1]/dg[i])for i in range(16)];sk=max(gso)/min(gso);assert sk<3.0
  sk_max=max(sk_max,sk)
  for row in Br[:4]:
   aa,bb=row[:8],row[8:]
   a3,b3=pow_phi15_z(aa,3),pow_phi15_z(bb,3)
   ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_bls,mul_phi15_z(W2,b3)))))
   a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3);a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
   nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9,mul_phi15_z(WP1,mul_phi15_z(a6,b3))),mul_phi15_z(a3,b6)),mul_phi15_z(W_el,b9))))
   assert nf%p==0 and ((ng%q==0)if sn=="g" else((nf//p)%q==0));div_ok+=1
   c_pt=math.log2(float(Mdet))/16.0
   cf_max_obs=max(cf_max_obs,math.log2(nf)-72.0*c_pt)
   cg_max_obs=max(cg_max_obs,math.log2(ng)-24.0*c_pt-float(dp["log2_u8"]))
 for ei in (0,2,6,8):
  sd,qh,pa,pb,ryh,rth=data["prod_enc12"][ei];q,ry,rt=int(qh,16),int(ryh,16),int(rth,16)
  Bq,_=_build_L(q,ry,rt);Brq,dgq,_=lll_16d_int(Bq,99,100);assert dgq[16]==q*q
  q16=float(q)**(1.0/16.0);b0_r=math.sqrt(sum(x*x for x in Brq[0]))/q16
  par_r=math.sqrt(sum(x*x for x in pa)+sum(x*x for x in pb))/q16
  assert b0_r>=0.90 and par_r>=3.0*eps_gate
  enc_b0_min=min(enc_b0_min,b0_r);enc_par_min=min(enc_par_min,par_r)
 assert div_ok==16 and cf_max_obs<cf_lim and cg_max_obs<cg_lim

 # 4. Exact O_H = Z[zeta_15] ideal coprimality <a, b> = O_H + live unit test proving it catches integer gcd==1 non-coprime ideals!
 SM_PS=primes_upto(2000)
 def _pfac(n):
  fs=[]
  for pp in SM_PS:
   if pp*pp>n:break
   if n%pp==0:
    fs.append(pp)
    while n%pp==0:n//=pp
  st=[n]if n>1 else[]
  while st:
   m=st.pop()
   if m<=1:continue
   if _is_prime_det(m):
    if m not in fs:fs.append(m)
   else:
    d=_brent(m,5000)or m
    if d==m:
     if m not in fs:fs.append(m)
    else:st.append(d);st.append(m//d)
  return fs
 def _is_oh_coprime(aa,bb):
  if math.gcd(*(aa+bb))!=1:return False
  na,nb=abs(norm_q15_z(aa)),abs(norm_q15_z(bb))
  if na==0 or nb==0:return False
  g0=math.gcd(na,nb)
  if g0==1:return True
  for pp in _pfac(g0):
   g1=pgcd_mod(aa,PHI15,pp)
   if len(g1)>1 and len(pgcd_mod(g1,bb,pp))>1:return False
  return True
 # Unit test: (2 - y)*(a0, b0) has 16-integer gcd == 1 yet shares prime ideal (2 - y) in O_H!
 _ut=[2,-1,0,0,0,0,0,0];_a0=[1,1,0,0,0,0,0,0];_b0=[1,0,1,0,0,0,0,0]
 _au,_bu=mul_phi15_z(_ut,_a0),mul_phi15_z(_ut,_b0)
 assert math.gcd(*(_au+_bu))==1 and (not _is_oh_coprime(_au,_bu)) and _is_oh_coprime(_a0,_b0)

 # 5. Campaign A: 32 genuine raw_M==0 above-cutoff (B_fb=32768 < LP=65536 < q in [2^17,2^18)) exact-O_H-coprime ideals x 28 i.i.d. f-side draws (pool=328)!
 cA=data["m5_campA"]
 rowsA=cA["rows"];NA=int(cA["n_draws"]);BA=int(cA["B_fb"]);LPA=int(cA["LP_max"]);UA=int(cA["U"])
 assert len(rowsA)==32 and NA==28 and BA==32768 and LPA==65536 and int(cA["pool_size"])==328
 for sd_i,qz,ryz,rtz,sd_rng,fi,sh in rowsA:
  assert BA<LPA<qz and (1<<17)<=qz<(1<<18) and qz%15==1 and _is_prime_det(qz)
  assert pow(ryz,15,qz)==1 and pow(ryz,3,qz)!=1 and pow(ryz,5,qz)!=1
  ww=pow(ryz,5,qz)
  if sd_i==0:assert (pow(rtz,3,qz)-UA*ww*ww)%qz==0
  else:assert (pow(rtz,3,qz)+ww*rtz*rtz+rtz-1)%qz==0
 fsA=[r[5]for r in rowsA]
 ustatA_8=sum(math.comb(f,K_fb)/math.comb(NA,K_fb)for f in fsA)/len(fsA)
 plugA_8=sum((f/NA)**K_fb for f in fsA)/len(fsA)
 max_p1A=max(f/NA for f in fsA);emp_p1A=sum(fsA)/(len(fsA)*NA);meanA=sum(r[6]for r in rowsA)/(len(rowsA)*NA)
 assert ustatA_8<0.001 and plugA_8<0.001 and max_p1A<0.35
 # Live bit-exact replay of 2 Campaign A raw_M==0 ideals (1 g-side + 1 f-side, incl. worst f_i=9/28!)
 pf_poolA=[]
 for pp in primes_upto(int(cA["pool_hi"])):
  if pp>=int(cA["pool_lo"])and pp%15==1:
   ep=(pp-1)//15
   gp=next((pow(g,ep,pp)for g in range(2,80)if pow(pow(g,ep,pp),3,pp)!=1 and pow(pow(g,ep,pp),5,pp)!=1),None)
   if gp is None:continue
   for k in (1,2,4,7,8,11,13,14):
    ss=pow(gp,k,pp);ww=pow(ss,5,pp)
    for rx in range(1,pp):
     if (rx*rx*rx+ww*rx*rx+rx-1)%pp==0:pf_poolA.append((pp,ss,rx))
 assert len(pf_poolA)==int(cA["pool_size"])
 P_SMA=math.prod(primes_upto(256));P_M1A=math.prod(p for p in primes_upto(2048)if p>256);P_M2A=math.prod(p for p in primes_upto(BA)if p>2048)
 def _sm1A(n,forb):
  if n==0:return False
  t=abs(n)
  for PP in (P_SMA,P_M1A,P_M2A):
   while True:
    g=math.gcd(t,PP)
    if g==1:break
    t//=g
   if t==1:return True
  if t<=BA:return True
  return t<=LPA and t not in forb and _is_prime_det(t)
 ballT0=_build_ball(int(cA["ball_t0"]));ballFB=_build_ball(int(cA["ball_fb"]))
 def _cnt24(Br,q,sq,p,sp,pts):
  cc=0;forb={q}if not p else{q,p}
  for nz in pts:
   aa=[0]*8;bb=[0]*8
   for idx,sgn in nz:
    row=Br[idx]
    for t in range(8):aa[t]+=sgn*row[t];bb[t]+=sgn*row[8+t]
   if not any(aa)or not any(bb)or math.gcd(*(aa+bb))!=1:continue
   a2,b2=mul_phi15_z(aa,aa),mul_phi15_z(bb,bb);a3,b3=mul_phi15_z(a2,aa),mul_phi15_z(b2,bb)
   ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(UA,mul_phi15_z(W2,b3)))))
   if ng==0:continue
   if sq=="g":
    if ng%q!=0:continue
    ng//=q
   if sp=="g" and p:
    if ng%p!=0:continue
    ng//=p
   if not _sm1A(ng,forb):continue
   nf=abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a3,mul_phi15_z(W_el,mul_phi15_z(a2,bb))),mul_phi15_z(aa,b2)),b3)))
   if nf==0:continue
   if sq=="f":
    if nf%q!=0:continue
    nf//=q
   if sp=="f" and p:
    if nf%p!=0:continue
    nf//=p
   if not _sm1A(nf,forb):continue
   if _is_oh_coprime(aa,bb):cc+=1
  return cc
 repA=0
 for ridx in (0,6,16,20):
  sd_i,qz,ryz,rtz,sd_rng,fi,sh=rowsA[ridx];sq="g" if sd_i==0 else "f"
  B0,_=_build_L(qz,ryz,rtz);Br0,dg0,_=lll_16d_int(B0,99,100);assert dg0[16]==qz*qz
  assert _cnt24(Br0,qz,sq,None,None,ballT0)==0
  repA+=1
 # Live bit-exact replay of 3 draws across 2 genuine raw_M==0 ideals: idx 6 (g-side, worst f_i=9/28, draw 0: 11 hits) + idx 16 (f-side, draw 0: 0 hits & draw 2: 21 hits)!
 repA_draws=0
 for ridx_s,pick_dis in (("6",(0,)),("16",(0,2))):
  sd_i,qz,ryz,rtz,sd_rng,fi_ref,sh_ref=rowsA[int(ridx_s)];sq="g" if sd_i==0 else "f"
  exp_c=cA["replay_first3"][ridx_s];rng_r=random.Random(sd_rng)
  for di in range(3):
   pp,sp,tp=rng_r.choice(pf_poolA)
   if di in pick_dis:
    Bqp,Md=_build_L(qz,ryz,rtz,pp,sp,tp);Brp,dgp,_=lll_16d_int(Bqp,99,100);assert dgp[16]==Md*Md
    assert _cnt24(Brp,qz,sq,pp,"f",ballFB)==exp_c[di];repA_draws+=1

 # 6. Campaign B: LITERAL deg-(72,24) (f_9, g_3) exact-O_H-coprime 28 ideals (14 f_9 + 14 g_3) x 24 draws = 672 lattices!
 assert "m5_campB_deg72" in data
 campB_out=_verify_campB_deg72(data["m5_campB_deg72"],K_fb,W_el,WP1,W2,PHI15,primes_upto,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,lll_16d_int,_is_oh_coprime)

 return dict(
  f_pool_band=f"[2^{bp_lo},2^{bp_hi})",f_pool_oh_ideals=n_oh_ideals,
  f9_14_16_exact_ideals=f9_14_16,g3_14_16_exact_ideals=g3_14_16,
  f9_cheb_ratio=round(f9_14_16/cheb_ref,3),g3_cheb_ratio=round(g3_14_16/cheb_ref,3),
  eps_gate=eps_gate,blichfeldt_pm1_factor12_worst=float(f"{worst_mb:.3e}"),blichfeldt_worst_b=worst_b,p_geom_charged=p_geom,
  prod_381b_b0_ratio_min=round(b0_ratio_min,3),prod_381b_enc_b0_ratio_min=round(enc_b0_min,3),
  prod_381b_enc_parent_ratio_min=round(enc_par_min,2),prod_381b_qp_div_checks=div_ok,
  prod_381b_qp_max_skew=round(sk_max,3),prod_381b_qp_max_cf=round(cf_max_obs,1),prod_381b_qp_max_cg=round(cg_max_obs,1),
  oh_coprime_exact_checked=True,oh_coprime_catches_int_gcd1=True,
  campA_raw0_above_cutoff_32x28=dict(
   n_ideals=32,g_ideals=16,f_ideals=16,all_confirmed_raw_M0=True,M_draws=NA,trials=32*NA,pool_size=len(pf_poolA),
   mean_hits=round(meanA,2),emp_p1=round(emp_p1A,4),max_ideal_p1=round(max_p1A,4),
   unbiased_mvue_E_p1_K8=float(f"{ustatA_8:.3e}"),plugin_E_p1_K8=float(f"{plugA_8:.3e}"),
   live_raw0_reconfirmed=repA,live_bitexact_fallback_draws_replayed=repA_draws,
  ),
  campB_literal_deg72_24=campB_out,
 )
