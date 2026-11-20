"""c6_adaptive.py — Live C6 O_h real-2-core-FB adaptive H5 & H3 scaling verifier."""
import base64,json,lzma,math
from pathlib import Path
from fractions import Fraction as F
def verify_threshold_counterexample():
 mu,rv,psi=60.0,1.0,0.10;bt=0.5*mu/rv;G=lambda t:(1.0+bt*(1.0-t))**(-rv)
 for i in range(81):t=i/100.0;assert 0.10*t+0.90*(t**100)<=G(t)+1e-12
 assert 0.10>G(psi)+0.05
 best=G(0.0)+psi
 for k in range(2,65):tk=min(0.8,(k-1)*(1.0+1.0/bt)/(k-1+rv));best=min(best,min(1.0,G(tk)/(tk**(k-1)))+(psi**k))
 assert best>=0.10-1e-12
 return dict(pgf_psi_invalid=round(G(psi),5),true_fail=0.10,tr_bound=round(best,5),tr_valid=True)
def verify_c6_adaptive_h5_live():
 import model_bn254_v6 as M,verify_bn254_v6 as V,descent_bn254_v6 as DESC
 ce=verify_threshold_counterexample()
 dat=json.loads(lzma.decompress(base64.b85decode((Path(__file__).resolve().parent/"c6_adaptive_h5.b85").read_bytes())).decode())
 names,cal,after_by,fstats=dat["names"],dat["cal_tbl"],dat["after_by_set"],dat["first_stats"]
 def _pw(si,cs,cb):
  k=f"{si}:{cs}:{cb}"
  if k in cal:n,sf,sw=cal[k];return sf/n,sw/n
  n=sf=sw=0
  for kk,(nn,sff,sww)in cal.items():
   _,ccs,ccb=map(int,kk.split(":"))
   if ccs==cs and ccb==cb:n+=nn;sf+=sff;sw+=sww
  return sf/n,sw/n
 per_set={};tot_f=tot_ef=tot_w=tot_ew=tot_an=tot_af=tot_shp=tot_shr=0
 for si,nm in enumerate(names):
  ar=after_by[si];an=len(ar);af=sum(x[1]for x in ar);shp=sum(x[3]for x in ar);shr=sum(x[4]for x in ar);of=ef=ow=ew=0.0
  for _,fl,wk,_,_,c1,c2 in ar:
   chs=[c1]+([c2]if c2 else[])
   of+=fl;ef+=sum(_pw(si,c>>5,c&31)[0]for c in chs);ow+=wk;ew+=sum(_pw(si,c>>5,c&31)[1]for c in chs)
  oe_f,oe_w=of/ef,ow/ew;assert oe_f<0.90 and 0.90<oe_w<1.10
  per_set[nm]=dict(roots=fstats[si][2],first_n=fstats[si][0],first_fail=fstats[si][1],after_n=an,after_fail=af,oe_fail=round(oe_f,3),oe_work=round(oe_w,3),sh_prev=shp,sh_run=shr)
  tot_an+=an;tot_af+=af;tot_f+=of;tot_ef+=ef;tot_w+=ow;tot_ew+=ew;tot_shp+=shp;tot_shr+=shr
 pooled_oe_f,pooled_oe_w=tot_f/tot_ef,tot_w/tot_ew
 assert tot_an==1856 and pooled_oe_f<0.85 and 0.95<pooled_oe_w<1.05
 conc=dat["conc_rows"];assert len(conc)==10 and max(x[5]for x in conc)<15.0 and max(x[6]for x in conc)<3.0
 U_T,P_T,LX_B,NS_MAX=int(dat["u_t"]),int(dat["p_t"]),int(dat["lx_b"]),int(dat["ns_max"]);B_FB=1<<LX_B
 ZKN,W_P=M.ZK_NUM,M.W_COEFFS;WBAR_P=[-4-c if i==0 else -c for i,c in enumerate(W_P)]
 y1=M._ppow_h([0,1,0,0,0,0],P_T,P_T);y2=M._ppow_h(y1,P_T,P_T);y3=M._ppow_h(y2,P_T,P_T);y6=M._ppow_h(y3,P_T**3,P_T)
 uw=[(U_T*c)%P_T for c in W_P];cur=uw;prd=[1,0,0,0,0,0]
 for _ in range(6):prd=M._pmul_h(prd,cur,P_T);cur=M._ppow_h(cur,P_T,P_T)
 tower_ok=(P_T==1+6*U_T and y6==[0,1,0,0,0,0]and y2!=[0,1,0,0,0,0]and y3!=[0,1,0,0,0,0]and pow(prd[0],(P_T-1)//2,P_T)==P_T-1);assert tower_ok
 sig_q=[F(c,3)for c in M.SIGMA_NUM];pw_s=[[F(1)]+[F(0)]*5,sig_q]
 for _ in range(4):pw_s.append(M._pmul_h_q(pw_s[-1],sig_q))
 SMY=tuple(tuple(int(3*pw_s[k][j])for k in range(6))for j in range(6));I6=[[1 if i==j else 0 for j in range(6)]for i in range(6)]
 _isp=M._isp;_sp=math.sumprod
 def _sqrtm(n,p):
  n%=p
  if n==0:return 0
  if pow(n,(p-1)//2,p)!=1:return -1
  if p%4==3:r=pow(n,(p+1)//4,p);return min(r,p-r)
  q,s=p-1,0
  while q%2==0:q//=2;s+=1
  z=2
  while pow(z,(p-1)//2,p)!=p-1:z+=1
  m,c,t,r=s,pow(z,q,p),pow(n,q,p),pow(n,(q+1)//2,p)
  while t!=1:
   i,tt=1,(t*t)%p
   while tt!=1:tt=(tt*tt)%p;i+=1
   b=pow(c,1<<(m-i-1),p);m,c,t,r=i,(b*b)%p,(t*b*b)%p,(r*b)%p
  return min(r,p-r)
 _mh=M._mh_int
 def _nc0c1(p):
  s1=tuple(_sp(r,p)for r in SMY);s2=tuple(_sp(r,s1)for r in SMY);n3=_mh(_mh(p,s1),s2);c1=-(n3[2]//729);c0=(n3[0]//729)-2*c1
  return c0,c1,c0*c0+2*c1*c1
 _hsp_c={}
 def _pinfo(p):
  v=_hsp_c.get(p)
  if v is not None:return v
  if p<=3 or(p&7)not in(1,3)or M._ppow_h([0,1,0,0,0,0],p,p)!=[0,1,0,0,0,0]:_hsp_c[p]=None;return None
  z0=_sqrtm(-2,p);z1=p-z0;rts=M._roots_mod(M.H_COEFFS,p)
  s0=next(r for r in rts if(2+V._eval_poly_mod(W_P,r,p))%p==z0);s1=next(r for r in rts if(2+V._eval_poly_mod(W_P,r,p))%p==z1)
  v=(z0,z1,(s0,s1),(_sqrtm(2+z0,p),_sqrtm(2+z1,p)),(_sqrtm(pow((U_T*(z0-2))%p,-1,p),p),_sqrtm(pow((U_T*(z1-2))%p,-1,p),p)));_hsp_c[p]=v;return v
 ret_f_sp={(int(p),int(o))for p,o in dat["ret_f_split"]};ret_g_sp={(int(p),int(o))for p,o in dat["ret_g_split"]}
 ret_f_ns={(int(p),int(o))for p,o in dat["ret_f_ns"]};ret_g_ns={(int(p),int(o))for p,o in dat["ret_g_ns"]}
 sp_ps=sorted({p for p,_ in ret_f_sp}|{p for p,_ in ret_g_sp});all_ps=M._primes_upto(B_FB);sp_set=set(sp_ps);ns_ps=[];ns_z={}
 for p in all_ps:
  if p<=NS_MAX and p not in sp_set:
   if _pinfo(p)is None:ns_ps.append(p);ns_z[p]=_sqrtm(-2,p)if(p>3 and(p&7)in(1,3))else -1
   else:sp_ps.append(p);sp_set.add(p)
 for p in sp_ps:_pinfo(p)
 def _pr(xs):
  xs=list(xs)
  while len(xs)>1:xs=[xs[i]*xs[i+1]if i+1<len(xs)else xs[i]for i in range(0,len(xs),2)]
  return xs[0]
 P_NS,P_SP,P_FB_ALL=_pr(ns_ps),_pr(sp_ps),_pr([p for p in all_ps if p>NS_MAX and p not in sp_set])
 bts=[];cur=[0]*6
 def _gb(idx,rem):
  if idx==6:
   if any(cur):bts.append(tuple(cur))
   return
  r=math.isqrt(rem)
  for v in range(-r,r+1):cur[idx]=v;_gb(idx+1,rem-v*v)
 _gb(0,int(dat["ball_r2"]));assert len(bts)==int(dat["ball_pts"])==484
 def _chk_fb(n,c0,c1,sd):
  if math.gcd(n,P_FB_ALL)>1:return None
  rem=n;r_ns=ret_f_ns if sd==0 else ret_g_ns;r_sp=ret_f_sp if sd==0 else ret_g_sp
  g_ns=math.gcd(rem,P_NS)
  if g_ns>1:
   for p in ns_ps:
    if g_ns%p==0:
     while rem%p==0:rem//=p
     z0=ns_z[p]
     if z0>=0:
      h0=((c0+c1*z0)%p==0);h1=((c0+c1*(p-z0))%p==0)
      if(not h0 and not h1)or(h0 and(p,0)not in r_ns)or(h1 and(p,1)not in r_ns):return None
     elif(p,0)not in r_ns:return None
  g_sp=math.gcd(rem,P_SP)
  if g_sp>1:
   for p in sp_ps:
    if g_sp%p==0:
     while rem%p==0:rem//=p
     z0,z1=_pinfo(p)[0],_pinfo(p)[1];h0=((c0+c1*z0)%p==0);h1=((c0+c1*z1)%p==0)
     if(not h0 and not h1)or(h0 and(p,0)not in r_sp)or(h1 and(p,1)not in r_sp):return None
  return rem
 def _fac(cof,nw,km):
  if cof==1:return[]
  if km<=0 or cof.bit_length()>nw*km:return None
  lim=1<<nw;out=[];st=[cof]
  while st:
   v=st.pop()
   if v==1:continue
   if _isp(v):
    if v>lim or v<=B_FB or _pinfo(v)is None:return None
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
 def _scan(sd,q,oi,nw,km):
  z0,z1,ss,rfs,rgs=_pinfo(q);s=ss[oi];rho=rfs[oi]if sd==0 else rgs[oi]
  zks=[V._eval_zk_mod([1 if j==k else 0 for j in range(6)],s,q)for k in range(6)];iz0=pow(zks[0],-1,q)
  B=[[0]*12 for _ in range(12)];B[0][0]=q
  for k in range(1,6):B[k][k]=1;B[k][0]=(-zks[k]*iz0)%q
  for k in range(6):B[6+k][6+k]=1;B[6+k][0]=(rho*zks[k]*iz0)%q
  DESC.bareiss_lll_12d(B,I6,1,1,False)
  a1=tuple(sum(B[0][c]*ZKN[c][j]for c in range(6))for j in range(6));b1=tuple(sum(B[0][6+c]*ZKN[c][j]for c in range(6))for j in range(6));a2=b2=None
  for r in range(1,12):
   ta=tuple(sum(B[r][c]*ZKN[c][j]for c in range(6))for j in range(6));tb=tuple(sum(B[r][6+c]*ZKN[c][j]for c in range(6))for j in range(6))
   if _mh(a1,tb)!=_mh(b1,ta):a2,b2=ta,tb;break
  if a2 is None:return[]
  a2s=[tuple(x//3 for x in _mh(ZKN[k],a2))for k in range(6)];b2s=[tuple(x//3 for x in _mh(ZKN[k],b2))for k in range(6)];rels=[]
  for bt in bts:
   pa=list(a1);pb=list(b1)
   for k,c in enumerate(bt):
    if c:
     ak,bk=a2s[k],b2s[k]
     for j in range(6):pa[j]+=c*ak[j];pb[j]+=c*bk[j]
   if not any(pa)or not any(pb):continue
   pa2=tuple(x//3 for x in _mh(pa,pa));pb2=tuple(x//3 for x in _mh(pb,pb))
   def _ev(s):
    if s==0:wb2=_mh(WBAR_P,pb2);return _nc0c1(tuple(pa2[j]+wb2[j]for j in range(6)))
    wa2=_mh(W_P,pa2);return _nc0c1(tuple(U_T*wa2[j]-pb2[j]for j in range(6)))
   c0a,c1a,na=_ev(sd)
   if na%q!=0 or na==0:continue
   ca=_chk_fb(na//q,c0a,c1a,sd);pa_f=_fac(ca,nw,km)if ca is not None else None
   if pa_f is None:continue
   c0b,c1b,nb=_ev(1-sd);cb=_chk_fb(nb,c0b,c1b,1-sd)if nb else None
   pb_f=_fac(cb,nw,km-len(pa_f))if cb is not None else None
   if pb_f is None:continue
   (pf,c0f,c1f),(pg,c0g,c1g)=((pa_f,c0a,c1a),(pb_f,c0b,c1b))if sd==0 else((pb_f,c0b,c1b),(pa_f,c0a,c1a))
   ch=[];ok=True
   for csd,ps,c0,c1 in((0,pf,c0f,c1f),(1,pg,c0g,c1g)):
    for p in ps:
     z0p,z1p,_,rfp,rgp=_pinfo(p);h=0
     for coi,zp in((0,z0p),(1,z1p)):
      if(c0+c1*zp)%p==0:
       if(rfp[coi]if csd==0 else rgp[coi])<=0 or len(ch)>=km:ok=False;break
       ch.append((csd,p,coi));h+=1
     if not ok or h==0:ok=False;break
    if not ok:break
   if ok:rels.append(ch)
  return rels
 def _desc(sd,q,oi,nw,km,dsh,rec=False):
  rels=_scan(sd,q,oi,nw,km);tw=1;prev=set();seq=[]
  for ri,ch in enumerate(rels):
   if not ch:return True,tw,seq
   aok=True;cw=0;sh=0
   for csd,cp,coi in ch:
    if ri>0 and(csd,cp,coi)in prev:sh=1
    cb=cp.bit_length();cnw=max(LX_B,cb-dsh);ckm=0 if cnw==LX_B else 1
    cok,cwk,_=_desc(csd,cp,coi,cnw,ckm,dsh,False);cw+=cwk
    if not cok:aok=False
   prev.update(ch);tw+=cw
   if rec:seq.append([ri,int(not aok),cw,len(ch),sh])
   if aok:return True,tw,seq
  return False,tw,seq
 live_w=0
 for si,sd,rq,roi,nw,km,dsh,exp_seq in dat["live_witnesses"]:
  _,_,got=_desc(sd,rq,roi,nw,km,dsh,True);assert got==exp_seq;live_w+=1
 return dict(ok=(live_w==3 and tower_ok and ce["tr_valid"]),threshold_counterexample=ce,tower_u2_p13_verified=tower_ok,roots_total=sum(x[2]for x in fstats),first_attempts_total=sum(x[0]for x in fstats),after_fail_total=tot_an,after_fail_failures=tot_af,pooled_oe_fail=round(pooled_oe_f,4),pooled_oe_work=round(pooled_oe_w,4),shared_prev_total=tot_shp,shared_run_total=tot_shr,conc_bins_checked=len(conc),max_C_emp=max(x[5]for x in conc),max_C_bday=max(x[6]for x in conc),live_witnesses_replayed=live_w,per_setting=per_set)
def verify_h3_scaling_and_nu_h():
 import model_bn254_v6 as M
 fr=M.load_data()["filtsim_raw"];scales=[20,21,22,23,24]
 m_sw=[];m_nl=[];m_k5=[];m_tr=[]
 for lu in scales:
  rs=fr[f"m184_z6_61.0_60.0_e0.82_W160_{lu}"];n=len(rs)
  m_sw.append(sum(r[1]for r in rs)/n);m_nl.append(sum(r[2]for r in rs)/n)
  m_k5.append(sum(r[6]*M.C_UNIT_5+(1.0-r[6])*M.C_GEN_1 for r in rs)/n);m_tr.append(sum(r[4]for r in rs)/n)
 all4_mono=all(arr[i]>arr[i+1]for arr in(m_sw,m_nl,m_k5,m_tr)for i in range(4))
 lx=[math.log(lu)for lu in scales];ly=[math.log(m)for m in m_sw];mx,my=sum(lx)/len(lx),sum(ly)/len(ly)
 beta=-sum((x-mx)*(y-my)for x,y in zip(lx,ly))/sum((x-mx)**2 for x in lx);A=math.exp(my+beta*mx)
 max_rel_err=max(abs(A/(lu**beta)-m)/m for lu,m in zip(scales,m_sw));pred_prod=A/(53.641526**beta)
 ps=M._primes_upto(16384);d3=d6=0
 for p in ps:
  if p in(2,3,19):continue
  xp=M._ppow_h([0,1,0,0,0,0],p,p);xp[1]=(xp[1]-1)%p
  if not any(xp):d6+=1
  else:
   d=len(M._pgcd(xp,M.H_COEFFS,p))-1
   if d==3:d3+=1
   elif d==6:d6+=1
   elif d!=0:raise AssertionError(f"d={d} p={p}")
 tot_orb=d3+2*d6;nu_h=tot_orb/(len(ps)/3.0)
 ok=(all4_mono and beta>3.0 and max_rel_err<0.03 and pred_prod<M.NU_SW_CAP/20.0 and M.NU_SW_CAP/m_sw[-1]>1.75 and M.NU_L_CAP/m_nl[-1]>1.75 and m_k5[-1]<M.K5_CAP and m_tr[-1]<M.TAB_R_CAP and 0.97<=nu_h<=1.03)
 return dict(ok=ok,fixed_instance_deterministic_h3=True,all4_monotone_20_24=all4_mono,inv_log_A=round(A,3),inv_log_beta=round(beta,3),max_rel_fit_err=round(max_rel_err,4),pred_nsw_prod=round(pred_prod,6),nu_sw_cap=M.NU_SW_CAP,cap_over_mean24=round(M.NU_SW_CAP/m_sw[-1],3),nu_l_cap_over_mean24=round(M.NU_L_CAP/m_nl[-1],3),nu_h_upto_2p14=round(nu_h,4),c6_d3_primes_2p14=d3,c6_d6_primes_2p14=d6,c6_deg1_orbits_2p14=tot_orb)
