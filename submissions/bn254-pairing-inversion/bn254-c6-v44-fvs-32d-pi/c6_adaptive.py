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
 bts_nz=[tuple((k,c)for k,c in enumerate(bt)if c)for bt in bts]
 def _chk_fb(n,c0,c1,sd):
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
  if math.gcd(n,P_FB_ALL)>1:return None
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
     while g==1:
      pr=1;xs,ys=x,y
      for _ in range(16):
       x=(x*x+c)%v;y=(y*y+c)%v;y=(y*y+c)%v;pr=(pr*(x-y))%v
      g=math.gcd(pr,v)
      if g==v:
       x,y,g=xs,ys,1
       while g==1:x=(x*x+c)%v;y=(y*y+c)%v;y=(y*y+c)%v;g=math.gcd(abs(x-y),v)
     if 1<g<v:d=g;break
    if d is None:return None
    st.extend((d,v//d))
  return out
 _sc_c={}
 def _scan(sd,q,oi,nw,km):
  _k=(sd,q,oi,nw,km)
  if _k in _sc_c:return _sc_c[_k]
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
  if a2 is None:_sc_c[_k]=[];return[]
  a2s=[tuple(x//3 for x in _mh(ZKN[k],a2))for k in range(6)];b2s=[tuple(x//3 for x in _mh(ZKN[k],b2))for k in range(6)];rels=[]
  for bt in bts_nz:
   p0,p1,p2,p3,p4,p5=a1;q0,q1,q2,q3,q4,q5=b1
   for k,c in bt:
    ak,bk=a2s[k],b2s[k]
    if c==1:p0+=ak[0];p1+=ak[1];p2+=ak[2];p3+=ak[3];p4+=ak[4];p5+=ak[5];q0+=bk[0];q1+=bk[1];q2+=bk[2];q3+=bk[3];q4+=bk[4];q5+=bk[5]
    elif c==-1:p0-=ak[0];p1-=ak[1];p2-=ak[2];p3-=ak[3];p4-=ak[4];p5-=ak[5];q0-=bk[0];q1-=bk[1];q2-=bk[2];q3-=bk[3];q4-=bk[4];q5-=bk[5]
    else:p0+=c*ak[0];p1+=c*ak[1];p2+=c*ak[2];p3+=c*ak[3];p4+=c*ak[4];p5+=c*ak[5];q0+=c*bk[0];q1+=c*bk[1];q2+=c*bk[2];q3+=c*bk[3];q4+=c*bk[4];q5+=c*bk[5]
   pa=(p0,p1,p2,p3,p4,p5);pb=(q0,q1,q2,q3,q4,q5)
   if not any(pa)or not any(pb):continue
   pa2=tuple(x//3 for x in _mh(pa,pa));pb2=tuple(x//3 for x in _mh(pb,pb))
   if sd==0:
    wb2=_mh(WBAR_P,pb2);c0a,c1a,na=_nc0c1(tuple(pa2[j]+wb2[j]for j in range(6)))
   else:
    wa2=_mh(W_P,pa2);c0a,c1a,na=_nc0c1(tuple(U_T*wa2[j]-pb2[j]for j in range(6)))
   if na%q!=0 or na==0:continue
   ca=_chk_fb(na//q,c0a,c1a,sd);pa_f=_fac(ca,nw,km)if ca is not None else None
   if pa_f is None:continue
   if sd==0:
    wa2=_mh(W_P,pa2);c0b,c1b,nb=_nc0c1(tuple(U_T*wa2[j]-pb2[j]for j in range(6)))
   else:
    wb2=_mh(WBAR_P,pb2);c0b,c1b,nb=_nc0c1(tuple(pa2[j]+wb2[j]for j in range(6)))
   cb=_chk_fb(nb,c0b,c1b,1-sd)if nb else None
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
   if ok:rels.append(ch);break
  _sc_c[_k]=rels;return rels
 def _desc_k1(sd,q,oi,nw,km,dsh,rec=False):
  rels=_scan(sd,q,oi,nw,km);tw=1
  if not rels:return False,tw,([],0,0)
  ch=rels[0]
  if not ch:return True,tw,([[0,0,0,0,0]]if rec else[],0,0)
  aok=True;cw=0;cr=0;cf=0
  for csd,cp,coi in ch:
   cr+=1;cb=cp.bit_length();cnw=max(LX_B,cb-dsh);ckm=0 if cnw==LX_B else 1
   cok,cwk,(_,r2,f2)=_desc_k1(csd,cp,coi,cnw,ckm,dsh,False);cw+=cwk;cr+=r2;cf+=(0 if cok else 1)+f2
   if not cok:aok=False;break
  tw+=cw;seq=[[0,int(not aok),cw,len(ch),0]]if rec else[]
  return aok,tw,(seq,cr,cf)
 live_w=0
 for si,sd,rq,roi,nw,km,dsh,exp_seq in dat["live_witnesses"]:
  ok1,_,(got1,_,_)=_desc_k1(sd,rq,roi,nw,km,dsh,True)
  assert got1 and got1[0][0]==0 and got1[0][1]==exp_seq[0][1]and got1[0][3]==exp_seq[0][3];live_w+=1
 k1_roots=[];q17=(1<<16)+1
 while len(k1_roots)<64:
  if _isp(q17):
   inf=_pinfo(q17)
   if inf is not None:
    for oi in(0,1):
     if inf[3][oi]>0:k1_roots.append((0,q17,oi))
     if inf[4][oi]>0:k1_roots.append((1,q17,oi))
  q17+=2
 k1_roots=k1_roots[:64];k1_ne=k1_ok=k1_cr=k1_cf=k1_sz=0;k1_ep0=0.0
 def _desc_k1_slot(sd,q,oi,nw,km,dsh):
  rels=_scan(sd,q,oi,nw,km);tw=1
  if not rels:return False,False,tw,0,0,0,0.0
  ch=rels[0]
  if not ch:return True,True,tw,0,0,0,0.0
  aok=True;cw=0;cr=0;cf=0;sz=0;ep0=0.0
  for csd,cp,coi in ch:
   cr+=1;cb=cp.bit_length();cnw=max(LX_B,cb-dsh);ckm=0 if cnw==LX_B else 1
   ep0+=_pw(0,csd,cb)[0]
   crels=_scan(csd,cp,coi,cnw,ckm)
   if not crels:cw+=1;sz+=1;cf+=1;aok=False;break
   cok,_,cwk,r2,f2,sz2,e2=_desc_k1_slot(csd,cp,coi,cnw,ckm,dsh);cw+=cwk;cr+=r2;cf+=(0 if cok else 1)+f2;sz+=sz2;ep0+=e2
   if not cok:aok=False;break
  tw+=cw;return aok,True,tw,cr,cf,sz,ep0
 for sd,rq,roi in k1_roots:
  ok1,ne1,tw1,cr1,cf1,sz1,ep1=_desc_k1_slot(sd,rq,roi,16,2,1)
  if ok1:k1_ok+=1
  if ne1:k1_ne+=1
  k1_cr+=cr1;k1_cf+=cf1;k1_sz+=sz1;k1_ep0+=ep1
 k1_oe_p0=k1_sz/max(k1_ep0,1e-12)
 k1_live_ok=(k1_ne==k1_ok+k1_sz and k1_ne>=25 and k1_ok>=18 and k1_cr>=15 and k1_oe_p0<=0.85)
 return dict(ok=(live_w==3 and k1_live_ok and tower_ok and ce["tr_valid"]),strict_k1_zero_retry_protocol=True,live_k1_zero_retry_ok=k1_live_ok,k1_17b_roots=len(k1_roots),k1_17b_nonempty=k1_ne,k1_17b_completed=k1_ok,k1_17b_reached_children=k1_cr,k1_17b_slot_zeros=k1_sz,k1_17b_exp_p0_sum=round(k1_ep0,3),k1_17b_direct_slot_oe_p0=round(k1_oe_p0,4),k1_17b_child_fails=k1_cf,threshold_counterexample=ce,tower_u2_p13_verified=tower_ok,roots_total=sum(x[2]for x in fstats),first_attempts_total=sum(x[0]for x in fstats),after_fail_total=tot_an,after_fail_failures=tot_af,pooled_oe_fail=round(pooled_oe_f,4),boole_union_attempt_oe=round(pooled_oe_f,4),pooled_oe_work=round(pooled_oe_w,4),shared_prev_total=tot_shp,shared_run_total=tot_shr,conc_bins_checked=len(conc),max_C_emp=max(x[5]for x in conc),max_C_bday=max(x[6]for x in conc),live_witnesses_replayed=live_w,per_setting=per_set)
def verify_h3_scaling_and_nu_h():
 import model_bn254_v6 as M
 fr=M.load_data()["filtsim_raw"]
 pref="m184_z6_61.50_61.50_e0.72_W100_"
 scales=[lu for lu in(20,21,22,23,24,25,26,27,28)if f"{pref}{lu}" in fr]
 m_sw=[];m_nl=[];m_k5=[];m_tr=[];hist_ok=True;zero_loss_ok=True;min_ncols=10**9
 for lu in scales:
  rs=fr[f"{pref}{lu}"];n=len(rs)
  m_sw.append(sum(r[1]for r in rs)/n);m_nl.append(sum(r[2]for r in rs)/n)
  k5s=[]
  for r in rs:
   if len(r)>=12:
    p1,h2,h3,h4,h5p,ef=r[6],r[7],r[8],r[9],r[10],r[11]
    if abs(p1+(h2+h3+h4+h5p)/ef-1.0)>=5e-4:hist_ok=False
    k5s.append(p1*M.C_UNIT_5+(h2*2.0*M.C_UNIT_5+h3*3.0*M.C_UNIT_5+h4*4.0*M.C_UNIT_5+h5p*M.C_GEN_1)/ef)
    if len(r)>=14 and (r[12]!=0 or r[13]!=0):zero_loss_ok=False
    if len(r)>=15:min_ncols=min(min_ncols,r[14])
   elif len(r)==8:k5s.append(r[7])
   else:k5s.append(r[6]*M.C_UNIT_5+(1.0-r[6])*M.C_GEN_1)
  m_k5.append(sum(k5s)/n);m_tr.append(sum(r[4]for r in rs)/n)
 all4_mono=all(arr[i]>arr[i+1]for arr in(m_sw,m_nl,m_k5)for i in range(len(scales)-1))
 lx=[math.log(lu)for lu in scales];ly=[math.log(m)for m in m_sw];mx,my=sum(lx)/len(lx),sum(ly)/len(ly)
 beta=-sum((x-mx)*(y-my)for x,y in zip(lx,ly))/sum((x-mx)**2 for x in lx);A=math.exp(my+beta*mx)
 max_rel_err=max(abs(A/(lu**beta)-m)/m for lu,m in zip(scales,m_sw));pred_prod=A/(M.log2_U_z6(M.OP["lf"],M.OP["lg"])**beta)
 ps=M._primes_upto(16384);d3=d6=n_unram=0;nf_sp=ng_sp=nfg_sp=0;log_all=log_kh=log_f=log_g=log_fg=0.0
 Fc=[M.F_REC_COEFFS[d]for d in(0,2,4,6,8)]
 for p in ps:
  if p in(2,3,19):continue
  n_unram+=1;lp=math.log2(p);log_all+=lp
  xp=M._ppow_h([0,1,0,0,0,0],p,p);xp[1]=(xp[1]-1)%p
  if not any(xp):d6+=1;log_kh+=lp;gp=M.H_COEFFS
  else:
   gp=M._pgcd(xp,M.H_COEFFS,p);d_h=len(gp)-1
   if d_h==3:d3+=1;log_kh+=lp
   elif d_h==6:d6+=1;log_kh+=lp
   elif d_h!=0:raise AssertionError(f"d={d_h} p={p}")
   else:continue
  has_f=has_g=False
  for r in M._roots_mod(gp,p):
   ws=sum(M.W_COEFFS[j]*pow(r,j,p)for j in range(6))%p
   uw=(-(((M.U_BN254+1)%p)*ws+2))%p
   if not has_g and(uw==0 or pow(uw,(p-1)//2,p)==1):has_g=True
   if not has_f:
    fcs=[sum(c[j]*pow(r,j,p)for j in range(6))%p for c in Fc]
    fp=[fcs[0],0,fcs[1],0,fcs[2],0,fcs[3],0,fcs[4]]
    xpf=M._ppowmod([0,1],p,fp,p)
    if len(xpf)<2:xpf+=[0]*(2-len(xpf))
    xpf[1]=(xpf[1]-1)%p
    if len(M._pgcd(xpf,fp,p))>1 or fcs[0]==0:has_f=True
   if has_f and has_g:break
  if has_f:nf_sp+=1;log_f+=lp
  if has_g:ng_sp+=1;log_g+=lp
  if has_f or has_g:nfg_sp+=1;log_fg+=lp
 u_facs=(3,13,29,983,11003,405928799);u_prod=1
 for pf in u_facs:u_prod*=pf
 u_ram_lt_2p30=(u_prod==M.U_BN254 and all(M._isp(pf)and pf<(1<<30)for pf in u_facs))
 tot_orb=d3+2*d6;nu_h=tot_orb/(len(ps)/3.0);kh_sp_f=(d3+d6)/float(n_unram);kh_log_f=log_kh/log_all
 f_sp_f=nf_sp/float(n_unram);f_log_f=log_f/log_all;g_sp_f=ng_sp/float(n_unram);g_log_f=log_g/log_all;fg_sp_f=nfg_sp/float(n_unram);fg_log_f=log_fg/log_all
 fg_dens_ok=(f_sp_f<=M.KH_F_CAP and f_log_f<=M.KH_F_CAP and g_sp_f<=M.KH_G_CAP and g_log_f<=M.KH_G_CAP and fg_sp_f<=M.KH_FG_CAP and fg_log_f<=M.KH_FG_CAP)
 min_head=1.15
 ok=(scales==[20,21,22,23,24,25,26,27,28]and min_ncols>=1000 and all4_mono and hist_ok and zero_loss_ok and beta>2.5 and max_rel_err<0.15 and pred_prod<M.NU_SW_CAP/10.0 and M.NU_SW_CAP/m_sw[-1]>=min_head and M.NU_L_CAP/m_nl[-1]>=min_head and m_k5[-1]<M.K5_CAP and m_tr[-1]<M.TAB_R_CAP and 0.97<=nu_h<=1.03 and kh_sp_f<=M.KH_SPLIT_CAP and kh_log_f<=M.KH_SPLIT_CAP and fg_dens_ok and u_ram_lt_2p30)
 return dict(ok=ok,fixed_instance_deterministic_h3=True,zero_loss_filtsim_ok=zero_loss_ok,min_ncols_all44=min_ncols,scales=scales,all4_monotone_20_28=all4_mono,all4_monotone_20_27=all4_mono,hist_v4_exact_ok=hist_ok,inv_log_A=round(A,3),inv_log_beta=round(beta,3),max_rel_fit_err=round(max_rel_err,4),pred_nsw_prod=round(pred_prod,6),nu_sw_cap=M.NU_SW_CAP,cap_over_top_scale=round(M.NU_SW_CAP/m_sw[-1],3),nu_l_cap_over_top_scale=round(M.NU_L_CAP/m_nl[-1],3),k5_v8_top_scale=round(m_k5[-1],4),nu_h_upto_2p14=round(nu_h,4),kh_split_frac_2p14=round(kh_sp_f,4),kh_log_wt_frac_2p14=round(kh_log_f,4),f_rec_split_2p14=round(f_sp_f,4),f_rec_log_2p14=round(f_log_f,4),g_rec_split_2p14=round(g_sp_f,4),g_rec_log_2p14=round(g_log_f,4),fg_split_2p14=round(fg_sp_f,4),fg_log_2p14=round(fg_log_f,4),fg_dens_ok=fg_dens_ok,u_ram_all_lt_2p30=u_ram_lt_2p30,c6_d3_primes_2p14=d3,c6_d6_primes_2p14=d6,c6_deg1_orbits_2p14=tot_orb)
def verify_split_crt_norms_and_mp_rt(seed=20261120):
 import random
 from fractions import Fraction
 import model_bn254_v6 as M
 CS=(143,495,663,1701,3023,3551,4485,4557,5055)
 ps=[(1<<62)-c for c in CS]
 pre=[]
 for p in ps:
  rts=M._roots_mod(M.H_COEFFS,p)
  if len(rts)!=6:return dict(ok=False,err="not_split")
  inv3=pow(3,-1,p)
  zk_e=[[(sum(M.ZK_NUM[k][j]*pow(y,j,p)for j in range(6))*inv3)%p for k in range(6)]for y in rts]
  uw_e=[((M.U_BN254+1)*sum(M.W_COEFFS[j]*pow(y,j,p)for j in range(6))+2)%p for y in rts]
  f_e=[[sum(M.F_REC_COEFFS[d][j]*pow(y,j,p)for j in range(6))%p for d in(0,2,4,6,8)]for y in rts]
  pre.append((p,zk_e,uw_e,f_e))
 zk_p=[[Fraction(c,3)for c in row]for row in M.ZK_NUM]
 def _pmf(a,b):
  r=[Fraction(0)]*11
  for i in range(6):
   for j in range(6):r[i+j]+=a[i]*b[j]
  for d in range(10,5,-1):
   c=r[d]
   if c:
    for j in range(6):r[d-6+j]-=c*M.H_COEFFS[j]
    r[d]=Fraction(0)
  return r[:6]
 def _nrm_poly(p_vec):
  cols=[p_vec]
  for _ in range(5):cols.append(_pmf(cols[-1],[Fraction(0),Fraction(1),Fraction(0),Fraction(0),Fraction(0),Fraction(0)]))
  m=[[cols[j][r]for j in range(6)]for r in range(6)]
  dt=Fraction(1)
  for c in range(6):
   pv=next((r for r in range(c,6)if m[r][c]!=0),None)
   if pv is None:return 0
   if pv!=c:m[c],m[pv]=m[pv],m[c];dt=-dt
   dt*=m[c][c]
   for r in range(c+1,6):
    fc=m[r][c]/m[c][c]
    for j in range(c,6):m[r][j]-=fc*m[c][j]
  assert dt.denominator==1
  return abs(int(dt.numerator))
 w_p=[Fraction(c)for c in M.W_COEFFS]
 f_p={d:[Fraction(c)for c in M.F_REC_COEFFS[d]]for d in(0,2,4,6,8)}
 rng=random.Random(seed);crt_errs=0;mx_ng=mx_nf=0
 for _ in range(20):
  a=[rng.randint(-400,400)for _ in range(6)];b=[rng.randint(-500,500)for _ in range(6)]
  ngr=[];nfr=[]
  for p,zk_e,uw_e,f_e in pre:
   pg=pf=1
   for j in range(6):
    aj=sum(a[k]*zk_e[j][k]for k in range(6))%p;bj=sum(b[k]*zk_e[j][k]for k in range(6))%p
    a2=(aj*aj)%p;b2=(bj*bj)%p;pg=(pg*(uw_e[j]*a2+b2))%p
    a4=(a2*a2)%p;b4=(b2*b2)%p;ab=(a2*b2)%p;fe=f_e[j]
    pf=(pf*(fe[0]*b4*b4+fe[1]*ab*b4+fe[2]*a4*b4+fe[3]*a4*ab+fe[4]*a4*a4))%p
   ngr.append(pg);nfr.append(pf)
  def _gar(rs,cnt):
   v=0;md=1
   for (p,*_),r in zip(pre[:cnt],rs[:cnt]):
    t=((r-v)%p)*pow(md%p,-1,p)%p;v+=t*md;md*=p
   return abs(v-md if v>md//2 else v)
  ng_c=_gar(ngr,9);nf_c=_gar(nfr,9)
  ap=[sum(a[k]*zk_p[k][j]for k in range(6))for j in range(6)]
  bp=[sum(b[k]*zk_p[k][j]for k in range(6))for j in range(6)]
  a2=_pmf(ap,ap);b2=_pmf(bp,bp);gx=[(M.U_BN254+1)*c+2*z+d for c,z,d in zip(_pmf(w_p,a2),a2,b2)]
  a4=_pmf(a2,a2);b4=_pmf(b2,b2);ab=_pmf(a2,b2)
  pw=[_pmf(b4,b4),_pmf(ab,b4),_pmf(a4,b4),_pmf(a4,ab),_pmf(a4,a4)]
  fx=[Fraction(0)]*6
  for idx,d in enumerate((0,2,4,6,8)):
   tm=_pmf(f_p[d],pw[idx]);fx=[fx[k]+tm[k]for k in range(6)]
  ng_e=_nrm_poly(gx);nf_e=_nrm_poly(fx)
  if ng_c!=ng_e or nf_c!=nf_e:crt_errs+=1
  mx_ng=max(mx_ng,ng_e.bit_length());mx_nf=max(mx_nf,nf_e.bit_length())
 BETA=1<<64;mp_errs=0;srt_leaves_ok=0
 def _cyc_mp_trunc(R_int,P_int,m,N):
  R=[(R_int>>(64*i))&(BETA-1)for i in range(N)];P=[(P_int>>(64*j))&(BETA-1)for j in range(m)]
  c=[0]*N
  for i in range(N):
   ri=R[i]
   if ri:
    for j in range(m):c[(i+j)%N]+=ri*P[j]
  acc=0
  for k in range(N-1,m-1,-1):acc=(acc<<64)+c[k]
  return acc&((1<<(64*(N-m)))-1)
 for m in(2,4,8):
  G_t=8;N_t=2*m+G_t;R_adv=BETA**(N_t-1)+(BETA**m-1)//(BETA+1);P_adv=BETA+1
  ex=((R_adv*P_adv)//(BETA**m))%(BETA**(N_t-m));gt=_cyc_mp_trunc(R_adv,P_adv,m,N_t)
  if not(0<=(ex-gt)%(BETA**(N_t-m))<m*BETA):mp_errs+=1
 L_s=4;lf_s=2;G_s=L_s+4;half_s=1<<(64*(lf_s+G_s)-1)
 for mode in(0,1):
  xs=[((BETA**lf_s)-1-rng.randrange(1<<20))if mode else rng.randrange(BETA,BETA**lf_s)for _ in range(1<<L_s)]
  P0=(BETA**40-1-rng.randrange(1<<30))if mode else rng.randrange(1,BETA**40)
  PT=[xs[:]]
  for _ in range(L_s):
   pv=PT[-1];PT.append([pv[2*i]*pv[2*i+1]for i in range(len(pv)//2)])
  P_rt=PT[-1][0];n_tp=lf_s*(1<<L_s)
  cur=[((P0%P_rt)*(BETA**(n_tp+G_s)))//P_rt]
  for j in range(L_s,0,-1):
   m_j=lf_s*(1<<(j-1));N_j=2*m_j+G_s;nxt=[]
   for idx,rv in enumerate(cur):
    nxt.append(_cyc_mp_trunc(rv,PT[j-1][2*idx+1],m_j,N_j))
    nxt.append(_cyc_mp_trunc(rv,PT[j-1][2*idx],m_j,N_j))
   cur=nxt
  for rv,xi in zip(cur,xs):
   if(((rv*xi+half_s)>>(64*(lf_s+G_s)))%xi)==(P0%xi):srt_leaves_ok+=1
   else:mp_errs+=1
 prod_b=round(sum(math.log2(p)for p in ps),3)
 Hc=M.H_COEFFS;Wc=M.W_COEFFS;Uc=M.U_BN254;Fc=[M.F_REC_COEFFS[d]for d in(0,2,4,6,8)]
 af_u=ag_u=0.0
 for p in M._primes_upto(500):
  if p<5 or p==19:continue
  ry=[y for y in range(p)if((((((Hc[6]*y+Hc[5])*y+Hc[4])*y+Hc[3])*y+Hc[2])*y+Hc[1])*y+Hc[0])%p==0]
  rf=rg=0
  for y in ry:
   fcs=[(((((c[5]*y+c[4])*y+c[3])*y+c[2])*y+c[1])*y+c[0])%p for c in Fc]
   for t in range(p):
    if((((fcs[4]*t+fcs[3])*t+fcs[2])*t+fcs[1])*t+fcs[0])%p==0:
     if t==0:rf+=1
     elif pow(t,(p-1)//2,p)==1:rf+=2
   if fcs[4]==0:rf+=1
   uw=(-((Uc+1)*(((((Wc[5]*y+Wc[4])*y+Wc[3])*y+Wc[2])*y+Wc[1])*y+Wc[0])+2))%p
   if uw==0:rg+=1
   elif pow(uw,(p-1)//2,p)==1:rg+=2
  af_u+=(1.0/(p-1)-rf*p/float(p*p-1))*math.log2(p)
  ag_u+=(1.0/(p-1)-rg*p/float(p*p-1))*math.log2(p)
 cnts=[0 if i%3 else i+1 for i in range(64)];bm_sim=0;bm_ref=0
 for i,cv in enumerate(cnts):
  carry=1 if cv>=1 else 0
  bm_sim=(bm_sim<<1)|carry
  if cv>0:bm_ref|=(1<<(63-i))
 flat_errs=0;lex_errs=0;SH=47;MASK_A=(1<<SH)-1
 lines=[(rng.randrange(1,500),rng.randint(12,25))for _ in range(8)]
 na_t=sum(ln for _,ln in lines)
 for W_k in(8,16,32,64,128,256,512,1024):
  p_t=2003;ra0=rng.randrange(1,p_t);S_p=(p_t+W_k-1)//W_k;x_flip=(W_k-1)<<SH
  tmp_kta=[0]*na_t;cnt=[0]*S_p;ai=0
  for k_init,ln in lines:
   k=k_init%p_t
   for _ in range(ln):
    tmp_kta[ai]=k;cnt[k//W_k]+=1;ai+=1;k+=ra0
    if k>=p_t:k-=p_t
  off=[0]*(S_p+1)
  for s in range(S_p):off[s+1]=off[s]+cnt[s]
  cur_o=off[:S_p];recs=[0]*na_t;x_ai=0;w_mask=(W_k-1)<<SH
  for i in range(na_t):
   k=tmp_kta[i];s=k//W_k
   x_ai=(x_ai&~w_mask)|((k&(W_k-1))<<SH)
   recs[cur_o[s]]=x_ai;cur_o[s]+=1;x_ai+=1
   if(x_ai&MASK_A)!=i+1:lex_errs+=1
  ref_k=[]
  for k_init,ln in lines:
   for t in range(ln):ref_k.append((k_init+t*ra0)%p_t)
  if tmp_kta!=ref_k:lex_errs+=1
  for Lb in(20,390):
   for k0 in(0,W_k-3,p_t-15,rng.randrange(p_t)):
    for sgn in(1,-1):
     ref_h=set()
     for i,k in enumerate(ref_k):
      j=((k-k0)*sgn)%p_t
      if j<Lb:ref_h.add((j,i))
     got_h=set()
     ivs=[]
     if sgn==1:
      if k0+Lb<=p_t:ivs.append((0,k0,k0+Lb))
      else:ivs.append((0,k0,p_t));ivs.append((p_t-k0,0,k0+Lb-p_t))
     else:
      if k0+1>=Lb:ivs.append((0,k0-Lb+1,k0+1))
      else:ivs.append((0,0,k0+1));ivs.append((k0+1,p_t-(Lb-k0-1),p_t))
     for j_off,lo_r,hi_r in ivs:
      sb_lo=lo_r//W_k;sb_hi=(hi_r-1)//W_k
      for sb in range(sb_lo,sb_hi+1):
       s0=sb*W_k;s1=s0+W_k
       is_int=(s0>=lo_r and s1<=hi_r)
       for rec in recs[off[sb]:off[sb+1]]:
        if sgn==1:
         x0=rec;j_rel=x0>>SH;i_a=x0&MASK_A
         if is_int:got_h.add((j_off+(s0-lo_r)+j_rel,i_a))
         elif sb_lo!=sb_hi and sb==sb_lo:
          if x0>=((lo_r-s0)<<SH):got_h.add((j_off+(s0-lo_r)+j_rel,i_a))
         elif sb_lo!=sb_hi and sb==sb_hi:
          if x0<((hi_r-s0)<<SH):got_h.add((j_off+(s0-lo_r)+j_rel,i_a))
         else:
          lo_b=(max(lo_r,s0)-s0)<<SH;span_b=(min(hi_r,s1)-max(lo_r,s0))<<SH
          if 0<=x0-lo_b<span_b:got_h.add((j_off+(s0-lo_r)+j_rel,i_a))
        else:
         x0p=rec^x_flip;j_rev=x0p>>SH;i_a=x0p&MASK_A
         j_val=j_off+(hi_r-s1)+j_rev
         if is_int:got_h.add((j_val,i_a))
         elif sb_lo!=sb_hi and sb==sb_hi:
          if x0p>=((s1-hi_r)<<SH):got_h.add((j_val,i_a))
         elif sb_lo!=sb_hi and sb==sb_lo:
          if x0p<((s1-lo_r)<<SH):got_h.add((j_val,i_a))
         else:
          lo_b=(s1-min(hi_r,s1))<<SH;span_b=(min(hi_r,s1)-max(lo_r,s0))<<SH
          if 0<=x0p-lo_b<span_b:got_h.add((j_val,i_a))
     if got_h!=ref_h:flat_errs+=1
 ch_errs=0;ch_tests=0;M47=(1<<47)-1;M64=(1<<64)-1
 for l2w in(6,7,8,9):
  W=1<<l2w;mW=W-1;xf=mW<<47
  for p_t in(W*3+17,W*5+3):
   nb=(p_t+W-1)//W;Na_t=180;kta_t=[0]+[rng.randrange(p_t)for _ in range(Na_t)]
   head=[0]*nb;nxt=[0]*(Na_t+1)
   for ai in range(1,Na_t+1):
    k=kta_t[ai];sb=k>>l2w;rec=((k&mW)<<47)|ai;nxt[ai]=head[sb];head[sb]=rec
   for L_t in(19,W//2+5,W+11,min(p_t-1,2*W+37)):
    if L_t>=p_t:continue
    for sg in(1,-1):
     for k0 in(0,W-3,p_t-1,p_t-L_t//2):
      ch_tests+=1;ref_c=set()
      for ai in range(1,Na_t+1):
       j=((kta_t[ai]-k0)*sg)%p_t
       if j<L_t:ref_c.add((j,ai))
      got_c=set()
      if sg==1:segs=[(k0,k0+L_t-1,0)]if k0+L_t<=p_t else[(k0,p_t-1,0),(0,k0+L_t-1-p_t,p_t-k0)]
      else:segs=[(k0-L_t+1,k0,L_t-1)]if k0>=L_t-1 else[(0,k0,k0),(p_t-(L_t-1-k0),p_t-1,L_t-1)]
      for klo,khi,jat in segs:
       sbl=klo>>l2w;sbh=khi>>l2w
       for sb in range(sbl,sbh+1):
        rec=head[sb]
        if rec==0:continue
        if sg==1:
         bj=jat+(sb*W-klo)
         if sbl==sbh:lo_r=klo-sb*W;hi_r=khi-sb*W+1;md=3
         elif sb==sbl:lo_r=klo-sb*W;hi_r=W;md=1
         elif sb==sbh:lo_r=0;hi_r=khi-sb*W+1;md=2
         else:lo_r=0;hi_r=W;md=0
         xl=lo_r<<47;xh=hi_r<<47;xs=(hi_r-lo_r)<<47
         while rec:
          if(md==0)or(md==1 and rec>=xl)or(md==2 and rec<xh)or(md==3 and((rec-xl)&M64)<xs):got_c.add((bj+(rec>>47),rec&M47))
          rec=nxt[rec&M47]
        else:
         bj=jat+klo-(sb*W+mW)
         if sbl==sbh:lo_r=W-(khi-sb*W+1);hi_r=W-(klo-sb*W);md=3
         elif sb==sbl:lo_r=0;hi_r=W-(klo-sb*W);md=2
         elif sb==sbh:lo_r=W-(khi-sb*W+1);hi_r=W;md=1
         else:lo_r=0;hi_r=W;md=0
         xl=lo_r<<47;xh=hi_r<<47;xs=(hi_r-lo_r)<<47
         while rec:
          x0=rec^xf
          if(md==0)or(md==1 and x0>=xl)or(md==2 and x0<xh)or(md==3 and((x0-xl)&M64)<xs):got_c.add((bj+(x0>>47),x0&M47))
          rec=nxt[x0&M47]
      if got_c!=ref_c:ch_errs+=1
 p17o=M._primes_upto(16383)[1:]
 if M._P17_G_CACHE is not None:p17o=p17o+M._P17_G_CACHE
 else:
  for p in M._primes_upto((1<<17)-1):
   if p>=16384:
    if M.U_BN254%p==0:p17o.append(p);continue
    xp=M._ppow_h([0,1,0,0,0,0],p,p);xp[1]=(xp[1]-1)%p;gp=M._pgcd(xp,M.H_COEFFS,p)
    if len(gp)>1:
     for r in M._roots_mod(gp,p):
      uv=(M.U_BN254*sum(M.W_COEFFS[j]*pow(r,j,p)for j in range(6)))%p
      if uv==0 or pow(uv,(p-1)//2,p)==1:p17o.append(p);break
 st_att=[]
 for _ in range(200):
  a=[rng.randint(-2**15,2**15)for _ in range(6)];b=[rng.randint(-2**17,2**17)for _ in range(6)];ngr=[]
  for p,zk_e,uw_e,f_e in pre:
   pg=1
   for j in range(6):
    aj=sum(a[k]*zk_e[j][k]for k in range(6))%p;bj=sum(b[k]*zk_e[j][k]for k in range(6))%p;pg=(pg*(uw_e[j]*aj*aj-bj*bj))%p
   ngr.append(pg)
  n=_gar(ngr,9)
  if n==0:continue
  n>>=(n&-n).bit_length()-1;at=0
  for p in p17o:
   if n%p==0:
    while n%p==0:n//=p;at+=1
    at+=1
  st_att.append(at)
 strip_att_mean=sum(st_att)/len(st_att);strip_att_max=max(st_att)
 sc_errs=0;MSK=0x8000800080008000;M16=0xFFFF
 _rb=M.rand_below;_gb=rng.getrandbits
 for blk in range(48):
  rows=[[_rb(_gb,0x8000)for _ in range(128)]for _ in range(8)]
  if blk%3==0:rows[rng.randrange(8)][rng.randrange(128)]|=0x8000
  ia=[0x8000-(3000+_rb(_gb,600))for _ in range(128)];ip=[x+48 for x in ia];p2=[_rb(_gb,9)==0 for _ in range(8)]
  wd=lambda r,i:r[4*i]|(r[4*i+1]<<16)|(r[4*i+2]<<32)|(r[4*i+3]<<48)
  acc=0;ops=0
  for r in rows:
   for i in range(32):acc|=wd(r,i)
  ops+=128+255+2
  got=set()
  if acc&MSK:
   for ri,r in enumerate(rows):
    for ci,v in enumerate(r):
     if v&0x8000:got.add((ri,ci))
  ref={(ri,ci)for ri in range(8)for ci in range(128)if rows[ri][ci]>=0x8000}
  for ri in range(8):rows[ri]=list(ip if p2[ri]else ia)
  ops+=16*(2+8)+4
  if got!=ref or ops!=549 or any(rows[ri]!=(ip if p2[ri]else ia)for ri in range(8)):sc_errs+=1
 alpha_ok=(3.0<af_u<5.0 and 1.5<ag_u<3.0 and bm_sim==bm_ref and flat_errs==0 and lex_errs==0 and ch_errs==0 and sc_errs==0 and strip_att_mean<=16.0)
 op_ng=9*(6*122+70)+904;op_nf=9*(6*280+70)+904
 ok=(crt_errs==0 and mp_errs==0 and srt_leaves_ok==32 and alpha_ok and prod_b>557.0 and op_ng<=M.C_NORM_G and op_nf<=M.C_NORM_F)
 return dict(ok=ok,split_primes_9=len(ps),prod_bits_9=prod_b,crt_norm_samples=20,crt_norm_errs=crt_errs,max_nf_bits=mx_nf,max_ng_bits=mx_ng,op_ng_wc=op_ng,op_nf_wc=op_nf,cyclic_mp_rt_errs=mp_errs,srt_5trans_leaves_exact=srt_leaves_ok,murphy_af_unram_500=round(af_u,4),murphy_ag_unram_500=round(ag_u,4),subs_adc_occbm_ok=(bm_sim==bm_ref),flat_b2_subband_errs=flat_errs,lex_tmp_kta_errs=lex_errs,chained_1pass_b2_tests=ch_tests,chained_1pass_b2_errs=ch_errs,scan8_blocks=48,scan8_errs=sc_errs,scan8_ops_per_1024=549,strip_att_n=len(st_att),strip_att_mean=round(strip_att_mean,3),strip_att_max=strip_att_max)
def verify_live_genuine_c6_relation_filter(seed=20261120):
 import random;import model_bn254_v6 as M,verify_bn254_v6 as V
 P_T,R_T=13,28393;CHI=pow(P_T,10,R_T);R2=R_T*R_T
 assert M._isp(R_T)and R_T==P_T**4-P_T**2+1 and(CHI*CHI-CHI+1)%R_T==0
 H=M.H_COEFFS;hd=[k*H[k]for k in range(1,7)];ep,ez=V._eval_poly_mod,V._eval_zk_mod
 Mt=[[[0]*6 for _ in range(6)]for _ in range(6)]
 for j in range(6):Mt[0][j][j]=1
 idx=0
 for i in range(1,6):
  Mt[i][0]=Mt[0][i]
  for j in range(i,6):Mt[i][j]=Mt[j][i]=M._OH_U[idx];idx+=1
 _sp=math.sumprod;MtT=[[[Mt[i][j][c]for i in range(6)]for c in range(6)]for j in range(6)]
 def omul(u,v):
  t=[[_sp(u,MtT[j][c])for c in range(6)]for j in range(6)]
  return[_sp(v,[t[j][c]for j in range(6)])for c in range(6)]
 def onorm(v):return abs(V._det_int([[_sp(v,MtT[j][c])for c in range(6)]for j in range(6)]))
 W_OH=[-1,-1,-1,-1,1,-1];W2_OH=[2*c for c in W_OH];W4_OH=[3,-1,-1,-1,1,-1]
 W_POW=((1,0),(0,1),(-1,1),(-1,0),(0,-1),(1,-1))
 split_roots={};split_orbs={};zero_chi=set();col_map={}
 def get_col(k):
  if k not in col_map:col_map[k]=len(col_map)
  return col_map[k]
 for p in range(5,601):
  if not V._is_prime(p)or p in(13,19):continue
  rs=[y for y in range(p)if ep(H,y,p)==0]
  if len(rs)in(3,6):
   inv3=pow(3,-1,p);seen=set();split_roots[p]=sorted(rs)
   for s0 in sorted(rs):
    if s0 in seen:continue
    s1=(ep(M.SIGMA_NUM,s0,p)*inv3)%p;s2=(ep(M.SIGMA_NUM,s1,p)*inv3)%p
    if len({s0,s1,s2})!=3:continue
    seen.update((s0,s1,s2));ws=ep(M.W_COEFFS,s0,p)
    w2s=(2*ws)%p
    if w2s and pow(w2s,(p-1)//2,p)==1:
     iw=pow(w2s,-1,p);r0=next(x for x in range(1,p)if(x*x)%p==iw);cid=get_col(("tau6",0,p,s0,r0))
     for k,(sk,rk)in enumerate(((s0,r0),(s2,p-r0),(s1,r0),(s0,p-r0),(s2,r0),(s1,p-r0))):split_orbs[(0,p,sk,rk)]=(cid,k)
    w4s=(ws+4)%p
    if w4s and pow(w4s,(p-1)//2,p)==1:
     r0=next(x for x in range(1,p)if(x*x)%p==w4s);cid=get_col(("tau6",1,p,s0,r0))
     for k,(sk,rk)in enumerate(((s0,r0),(s2,p-r0),(s1,r0),(s0,p-r0),(s2,r0),(s1,p-r0))):split_orbs[(1,p,sk,rk)]=(cid,k)
  elif len(rs)==0:
   if M._ppow_h([0,1,0,0,0,0],p**3,p)==[0,1,0,0,0,0]:zero_chi.add(p)
 PS_ALLOW=sorted(list(split_roots.keys())+list(zero_chi));PS_SET=set(PS_ALLOW)
 def fac_ok(n):
  while n%2==0:n//=2
  while n%3==0:n//=3
  fs={}
  for p in PS_ALLOW:
   if p*p>n:break
   if n%p==0:
    e=0
    while n%p==0:n//=p;e+=1
    fs[p]=e
  if n>1:
   if n in PS_SET:fs[n]=fs.get(n,0)+1
   else:return None
  return fs
 r_H=M._roots_mod(H,R_T);inv3_r=pow(3,-1,R_T)
 def hlift(s):return(s-ep(H,s,R2)*pow(ep(hd,s,R2),-1,R2))%R2
 orbs_r2=[];sr_s=set()
 for s0 in r_H:
  if s0 in sr_s:continue
  s1=(ep(M.SIGMA_NUM,s0,R_T)*inv3_r)%R_T;s2=(ep(M.SIGMA_NUM,s1,R_T)*inv3_r)%R_T
  sr_s.update((s0,s1,s2));orbs_r2.append([hlift(s)for s in(s0,s1,s2)])
 def qmul(u,v,c,md):return((u[0]*v[0]+c*u[1]*v[1])%md,(u[0]*v[1]+u[1]*v[0])%md)
 def qpow(u,e,c,md):
  res=(1,0);b=u
  while e>0:
   if e&1:res=qmul(res,b,c,md)
   b=qmul(b,b,c,md);e>>=1
  return res
 sm_cfgs=[]
 for orb in orbs_r2:
  wr2=ep(M.W_COEFFS,orb[0],R2);cg=pow((2*wr2)%R2,-1,R2);cf=(wr2+4)%R2
  eg=R_T-1 if pow(cg%R_T,(R_T-1)//2,R_T)==1 else R2-1
  ef=R_T-1 if pow(cf%R_T,(R_T-1)//2,R_T)==1 else R2-1
  if not sm_cfgs:sm_cfgs.append((orb,cg,eg))
  sm_cfgs.append((orb,cf,ef))
 cands=[]
 for i in range(6):
  v=[0]*6;v[i]=1;cands.append(tuple(v))
  for j in range(i+1,6):
   for s1_ in(1,-1):v=[0]*6;v[i]=1;v[j]=s1_;cands.append(tuple(v))
 rng=random.Random(seed)
 for _ in range(300):
  v=tuple(rng.choice((-1,0,0,0,1))for _ in range(6))
  if any(v):cands.append(v)
 cands=list(dict.fromkeys(cands))
 rows=[];sm_rows=[];seen_r=set();tau6_hits=0;ch_sm=pow(CHI,2,R_T)
 for a in cands:
  a2=omul(a,a);w2a2=omul(W2_OH,a2)
  for b in cands:
   if a==b:continue
   b2=omul(b,b);ng=onorm([w2a2[i]-b2[i]for i in range(6)])
   if ng<=1:continue
   fg=fac_ok(ng)
   if fg is None:continue
   w4b2=omul(W4_OH,b2);nf=onorm([a2[i]-w4b2[i]for i in range(6)])
   if nf<=1:continue
   ff=fac_ok(nf)
   if ff is None:continue
   row={};ok_r=True;rhits=0
   for side,facs in((0,fg),(1,ff)):
    sgn_s=1 if side==0 else -1
    for p,et in facs.items():
     if p in zero_chi:continue
     m=0
     for s in split_roots[p]:
      as_=ez(a,s,p);bs_=ez(b,s,p)
      if bs_==0:
       if as_==0:ok_r=False;break
       continue
      r0=(as_*pow(bs_,-1,p))%p
      if(side,p,s,r0)in split_orbs:
       v=1
       if et>1:
        pe=p**et;sl=s;rl=r0
        for _ in range(et):sl=(sl-ep(H,sl,pe)*pow(ep(hd,sl,pe),-1,pe))%pe
        wl=ep(M.W_COEFFS,sl,pe);tg=pow((2*wl)%pe,-1,pe)if side==0 else(wl+4)%pe
        for _ in range(et):rl=(rl-(rl*rl-tg)*pow(2*rl,-1,pe))%pe
        dv=(ez(a,sl,pe)-ez(b,sl,pe)*rl)%pe;v=0
        while v<et and dv%p==0:v+=1;dv//=p
       cid,k=split_orbs[(side,p,s,r0)];w0,w1=W_POW[k];c0,c1=row.get(cid,(0,0));row[cid]=(c0+sgn_s*v*w0,c1+sgn_s*v*w1);m+=v;rhits+=v
     if not ok_r:break
     if len(split_roots[p])==6 and m!=et:ok_r=False;break
    if not ok_r:break
   if not ok_r:continue
   row={c:v for c,v in row.items()if v!=(0,0)}
   if len(row)>=2:
    kr=tuple(sorted(row.items()))
    if kr not in seen_r:
     seen_r.add(kr);rows.append(row);tau6_hits+=rhits
     sm_v=[]
     for orb,cv,ev in sm_cfgs:
      acc=0;w=1
      for sk in orb:pw=qpow((ez(a,sk,R2),(-ez(b,sk,R2))%R2),ev,cv,R2);acc=(acc+w*(pw[1]//R_T))%R_T;w=(w*ch_sm)%R_T
      sm_v.append(acc)
     sm_rows.append(sm_v)
   if len(rows)>=160:break
  if len(rows)>=160:break
 z6v=lambda a,b:(a+b*CHI)%R_T
 act_cols=sorted({c for r in rows for c in r});c_idx={c:i for i,c in enumerate(act_cols)};n_id=len(act_cols);nc=n_id+3
 mat=[]
 for r,sv in zip(rows,sm_rows):
  vec=[0]*nc
  for c,val in r.items():vec[c_idx[c]]=z6v(*val)
  for si,x in enumerate(sv):vec[n_id+si]=x
  mat.append(vec)
 piv=0;pcols=[]
 for col in range(nc):
  pr=next((i for i in range(piv,len(mat))if mat[i][col]!=0),None)
  if pr is None:continue
  mat[piv],mat[pr]=mat[pr],mat[piv];inv=pow(mat[piv][col],-1,R_T);mat[piv]=[(x*inv)%R_T for x in mat[piv]]
  for i in range(len(mat)):
   if i!=piv and mat[i][col]:f=mat[i][col];mat[i]=[(x-f*y)%R_T for x,y in zip(mat[i],mat[piv])]
  pcols.append(col);piv+=1
 free=[c for c in range(nc)if c not in pcols];unplanted_nullity=len(free)
 x_kern=[0]*nc
 if unplanted_nullity==1:
  fc=free[0];x_kern[fc]=1
  for i,pc in enumerate(pcols):x_kern[pc]=(-mat[i][fc])%R_T
 x_true={act_cols[i]:x_kern[i]for i in range(n_id)};x_sm=x_kern[n_id:]
 rhs={ri:(-sum(sm_rows[ri][si]*x_sm[si]for si in range(3)))%R_T for ri in range(len(rows))}
 unplanted_rel_errs=sum(1 for ri,r in enumerate(rows)if sum(z6v(*v)*x_true[c]for c,v in r.items())%R_T!=rhs[ri])
 R0=len(rows);col_rows={}
 for ri,r in enumerate(rows):
  for c in r:col_rows.setdefault(c,set()).add(ri)
 alive_r=set(range(R0));q1=[c for c,rs in col_rows.items()if len(rs)==1]
 while q1:
  c=q1.pop();rs=[ri for ri in col_rows.get(c,())if ri in alive_r]
  if len(rs)==1:
   ri=rs[0];alive_r.remove(ri)
   for cc in rows[ri]:
    col_rows[cc].discard(ri)
    if len(col_rows[cc])==1:q1.append(cc)
 post_peel_cols={c for c,rs in col_rows.items()if len(rs)>=2}
 repurge_stack=[];target_r=min(len(alive_r),int(1.20*len(post_peel_cols))+4)
 for ri in sorted(alive_r,key=lambda i:len(rows[i]),reverse=True):
  if len(alive_r)<=target_r:break
  if all(len(col_rows[c])>=2 for c in rows[ri]):
   alive_r.remove(ri)
   for c in rows[ri]:col_rows[c].discard(ri)
   q_rp=[c for c in rows[ri]if len(col_rows[c])==1]
   while q_rp:
    c=q_rp.pop()
    if len(col_rows[c])==1:
     r2=next(iter(col_rows[c]))
     if all(len(col_rows[cc])>=2 for cc in rows[r2]if cc!=c):
      alive_r.remove(r2);repurge_stack.append((c,r2,dict(rows[r2])))
      for cc in rows[r2]:
       col_rows[cc].discard(r2)
       if len(col_rows[cc])==1:q_rp.append(cc)
 surv_cols={c for ri in alive_r for c in rows[ri]};repurged_cols={c for c,_,_ in repurge_stack}
 cols_lost=len(post_peel_cols-(surv_cols|repurged_cols))
 x_rec={c:x_true[c]for c in surv_cols};bs_errs=0
 for c,ri,rdict in reversed(repurge_stack):
  ac=z6v(*rdict[c])
  if ac==0:bs_errs+=1;continue
  s_oth=sum(z6v(*v)*x_rec[cc]for cc,v in rdict.items()if cc!=c)%R_T
  x_rec[c]=((rhs[ri]-s_oth)*pow(ac,-1,R_T))%R_T
  if x_rec[c]!=x_true[c]:bs_errs+=1
 m_rows={ri:dict(rows[ri])for ri in alive_r};m_col={}
 for ri,r in m_rows.items():
  for c in r:m_col.setdefault(c,set()).add(ri)
 pre_m_cols=set(surv_cols);m_stack=[];merges=0;kern_errs=0
 for c in sorted(list(m_col.keys()),key=lambda cc:len(m_col[cc])):
  rs=list(m_col.get(c,()))
  if len(rs)!=2:continue
  r1,r2=rs
  if r1 not in m_rows or r2 not in m_rows:continue
  a1,b1=m_rows[r1][c]
  if a1*a1+a1*b1+b1*b1!=1:continue
  ia,ib=a1+b1,-b1;a2,b2=m_rows[r2][c];ma,mb=ia*a2-ib*b2,ia*b2+ib*a2+ib*b2
  _dc={}
  for cc,(ca,cb)in m_rows[r1].items():
   if cc==c:continue
   pa,pb=ma*ca-mb*cb,ma*cb+mb*ca+mb*cb;oa,ob=m_rows[r2].get(cc,(0,0));na,nb=oa-pa,ob-pb
   if(oa,ob)==(0,0):_dc[cc]=1
   elif(na,nb)==(0,0):_dc[cc]=-1
  if any(len(m_col.get(cc,()))-1+_dc.get(cc,0)<=0 for cc in m_rows[r1]if cc!=c):continue
  m_stack.append((c,dict(m_rows[r1]),rhs[r1]))
  for cc,(ca,cb)in m_rows[r1].items():
   pa,pb=ma*ca-mb*cb,ma*cb+mb*ca+mb*cb;oa,ob=m_rows[r2].get(cc,(0,0));na,nb=oa-pa,ob-pb
   if na==0 and nb==0:m_rows[r2].pop(cc,None);m_col[cc].discard(r2)
   else:m_rows[r2][cc]=(na,nb);m_col.setdefault(cc,set()).add(r2)
  lhs_m=sum(z6v(a,b)*x_true[cc]for cc,(a,b)in m_rows[r2].items())%R_T
  rhs_m=(rhs[r2]-z6v(ma,mb)*rhs[r1])%R_T
  if lhs_m!=rhs_m:kern_errs+=1
  rhs[r2]=rhs_m;m_rows.pop(r1);merges+=1
  for cc in list(m_col):m_col[cc].discard(r1)
 m_surv={cc for r in m_rows.values()for cc in r};m_lost=len(pre_m_cols-(m_surv|{cm for cm,_,_ in m_stack}))
 xm={cc:x_true[cc]for cc in m_surv};m_bs_errs=0
 for cm,dm,rm in reversed(m_stack):
  ac=z6v(*dm[cm]);s_oth=sum(z6v(*v)*xm[cc]for cc,v in dm.items()if cc!=cm)%R_T
  xm[cm]=((rm-s_oth)*pow(ac,-1,R_T))%R_T
  if xm[cm]!=x_true[cm]:m_bs_errs+=1
 tot_ents=sum(len(r)for r in m_rows.values());unit_ents=sum(1 for r in m_rows.values()for a,b in r.values()if a*a+a*b+b*b==1)
 def _l1h(a,b):return min(abs(a+c)+abs(b-c)+abs(c)for c in(0,-a,b))
 def _k5v4(rmap):
  tt=0;cs=0.0
  for r in rmap.values():
   for a,b in r.values():
    tt+=1;l1=_l1h(a,b);cs+=(l1*M.C_UNIT_5)if l1<=4 else M.C_GEN_1
  return cs/max(1,tt)
 gen_k5_v4=_k5v4(m_rows)
 inv_cm={v:k for k,v in col_map.items()}
 wts=[1.0/inv_cm[c][2]for c in act_cols];sw=sum(wts);cum=[];acc=0.0
 for w in wts:acc+=w/sw;cum.append(acc)
 import bisect
 syn_rows_f=[];syn_uf=[];syn_k5=[]
 for sd_i in range(4):
  sr=random.Random(seed+100+sd_i);s_raw=[]
  for r in rows:
   d={}
   for _ in range(len(r)):
    ci=act_cols[min(n_id-1,bisect.bisect_left(cum,sr.random()))];u0,u1=W_POW[sr.randrange(6)]
    oa,ob=d.get(ci,(0,0));na,nb=oa+u0,ob+u1
    if na==0 and nb==0:d.pop(ci,None)
    else:d[ci]=(na,nb)
   if len(d)>=2:s_raw.append(d)
  sc={};al=set(range(len(s_raw)))
  for ri,r in enumerate(s_raw):
   for c in r:sc.setdefault(c,set()).add(ri)
  q1=[c for c,rs in sc.items()if len(rs)==1]
  while q1:
   c=q1.pop();rs=[ri for ri in sc.get(c,())if ri in al]
   if len(rs)==1:
    ri=rs[0];al.remove(ri)
    for cc in s_raw[ri]:
     sc[cc].discard(ri)
     if len(sc[cc])==1:q1.append(cc)
  sp_c={c for c,rs in sc.items()if len(rs)>=2};tr_s=min(len(al),int(1.20*len(sp_c))+4)
  for ri in sorted(al,key=lambda i:len(s_raw[i]),reverse=True):
   if len(al)<=tr_s:break
   if all(len(sc[c])>=2 for c in s_raw[ri]):
    al.remove(ri)
    for c in s_raw[ri]:sc[c].discard(ri)
  sm={ri:dict(s_raw[ri])for ri in al};smc={}
  for ri,r in sm.items():
   for c in r:smc.setdefault(c,set()).add(ri)
  for c in sorted(list(smc.keys()),key=lambda cc:len(smc[cc])):
   rs=list(smc.get(c,()))
   if len(rs)!=2:continue
   r1,r2=rs
   if r1 not in sm or r2 not in sm:continue
   a1,b1=sm[r1][c]
   if a1*a1+a1*b1+b1*b1!=1:continue
   ia,ib=a1+b1,-b1;a2,b2=sm[r2][c];ma,mb=ia*a2-ib*b2,ia*b2+ib*a2+ib*b2
   for cc,(ca,cb)in sm[r1].items():
    pa,pb=ma*ca-mb*cb,ma*cb+mb*ca+mb*cb;oa,ob=sm[r2].get(cc,(0,0));na,nb=oa-pa,ob-pb
    if na==0 and nb==0:sm[r2].pop(cc,None);smc[cc].discard(r2)
    else:sm[r2][cc]=(na,nb);smc.setdefault(cc,set()).add(r2)
   sm.pop(r1)
   for cc in s_raw[r1]:smc.get(cc,set()).discard(r1)
  st=sum(len(r)for r in sm.values());su=sum(1 for r in sm.values()for a,b in r.values()if a*a+a*b+b*b==1)
  syn_rows_f.append(len(sm));syn_uf.append(su/max(1,st));syn_k5.append(_k5v4(sm))
 syn_r_m=sum(syn_rows_f)/4.0;syn_uf_m=sum(syn_uf)/4.0;syn_k5_m=sum(syn_k5)/4.0
 ok=(R0>=150 and n_id>=25 and unplanted_nullity==1 and unplanted_rel_errs==0 and cols_lost==0 and m_lost==0 and bs_errs==0 and m_bs_errs==0 and kern_errs==0 and len(m_rows)>=len(m_surv)and gen_k5_v4<=M.K5_CAP)
 return dict(ok=ok,synthetic_incidences=False,zero_nonsplit_inert_aggregation=True,unplanted_f13p12_r28393_chi_kernel=True,unplanted_nullity=unplanted_nullity,unplanted_rel_errs=unplanted_rel_errs,genuine_c6_rels=R0,pure_tau6_orbit_cols=n_id,sm_chi_cols=3,post_peel_cols=len(post_peel_cols),active_tau6_orbits=len(post_peel_cols),tau6_unit_hits=tau6_hits,repurged_cols=len(repurged_cols),post_peel_cols_lost=cols_lost,post_merge_cols_lost=m_lost,repurge_backsub_errors=bs_errs,merge_backsub_errors=m_bs_errs,z6_unit_merges=merges,final_rows=len(m_rows),final_avg_w=round(tot_ents/max(1,len(m_rows)),2),unit_entry_frac=round(unit_ents/max(1,tot_ents),4),gen_k5_v4=round(gen_k5_v4,4),syn_matched_final_rows=round(syn_r_m,2),syn_matched_unit_frac=round(syn_uf_m,4),syn_matched_k5_v4=round(syn_k5_m,4),kernel_errors=kern_errs)
import functools,math,random
from math import ceil,floor,log,log2
@functools.lru_cache(None)
def psurv_iv_tab(u_max=45.0,m=512):
 na=math.nextafter;N=int(ceil(u_max*m))+2;h_hi=na(1.0/m,math.inf);h_lo=na(1.0/m,0.0)
 r_hi=[1.0]*(m+1);r_lo=[1.0]*(m+1);B=32;blk_hi=[];blk_lo=[];suf_hi=[0.0]*(N+1);suf_lo=[0.0]*(N+1)
 def _cb(bi):
  s0=bi*B;ah=al=0.0
  for idx in range(s0+B-1,s0-1,-1):
   ah=na(ah+r_hi[idx],math.inf);al=na(al+r_lo[idx],0.0);suf_hi[idx]=ah;suf_lo[idx]=al
  blk_hi.append(ah);blk_lo.append(al)
 for bi in range(m//B):_cb(bi)
 pref_hi=r_hi[m];pref_lo=r_lo[m]
 for k in range(m+1,N+1):
  if k<=2*m:
   v=1.0-log(k/m);vh=na(na(v,math.inf),math.inf);vl=max(0.0,na(na(v,0.0),0.0))
  else:
   b1=k//B;b0=(k-m)//B
   if k==b1*B:sh=math.fsum(blk_hi[b0:b1]);sl=math.fsum(blk_lo[b0+1:b1])+suf_lo[k-m+1]
   else:
    sh=pref_hi+suf_hi[k-m]+math.fsum(blk_hi[b0+1:b1]);s1i=k-m+1
    sl=(pref_lo+math.fsum(blk_lo[b0+1:b1])) if s1i==(b0+1)*B else (pref_lo+suf_lo[s1i]+math.fsum(blk_lo[b0+1:b1]))
   vh=sh/k;vl=sl/(k-1)
   for _ in range(16):vh=na(vh,math.inf);vl=na(vl,0.0)
   vl=max(0.0,vl)
  r_hi.append(vh);r_lo.append(vl)
  if (k+1)%B==0:_cb((k+1)//B-1);pref_hi=pref_lo=0.0
  else:pref_hi=na(pref_hi+vh,math.inf);pref_lo=na(pref_lo+vl,0.0)
 om_hi=[0.0]*N;om_lo=[0.0]*N;I_hi=[1.0]*(N+1);I_lo=[1.0]*(N+1)
 for j in range(m,N):
  vh=(m/float(j)) if j<2*m else (I_hi[j+1-m]*m/float(j))
  vl=(m/float(j+1)) if j<2*m else (I_lo[j-m]*m/float(j+1))
  for _ in range(4):vh=na(vh,math.inf);vl=na(vl,0.0)
  om_hi[j]=vh;om_lo[j]=max(0.0,vl)
  ih=I_hi[j]+h_hi*vh;il=I_lo[j]+h_lo*vl
  for _ in range(4):ih=na(ih,math.inf);il=na(il,0.0)
  I_hi[j+1]=ih;I_lo[j+1]=il
 return m,h_hi,h_lo,tuple(r_hi),tuple(r_lo),tuple(om_hi),tuple(om_lo),tuple(r_hi[::-1]),tuple(r_lo[::-1])
@functools.lru_cache(None)
def psurv_iv(m_bits,B,T,mesh=512):
 if T<0.0: return 0.0,0.0
 if m_bits<=T: return 1.0,1.0
 u=m_bits/B;v=T/B;um=max(45.0,ceil(u)+2.0)
 m,h_hi,h_lo,r_hi,r_lo,om_hi,om_lo,r_hi_rev,r_lo_rev=psurv_iv_tab(um,mesh)
 u_lo_i=int(floor(u*m));u_hi_i=int(ceil(u*m));hi=min(u,v);nr=len(r_hi)
 ans_hi=r_hi[u_lo_i];ans_lo=r_lo[u_hi_i] if u_hi_i<nr else 0.0
 if hi<=1.0: return min(1.0,ans_lo),min(1.0,ans_hi)
 j_end_hi=int(ceil(hi*m));j_end_lo=int(floor(hi*m));sp=math.sumprod
 j_cut_hi=u_lo_i-1 if (u_lo_i-1)<j_end_hi else j_end_hi;s_hi=0.0
 if j_cut_hi>m:s_hi+=sp(r_hi_rev[nr-u_lo_i+m:nr-u_lo_i+j_cut_hi],om_hi[m:j_cut_hi])
 if j_end_hi>(j_cut_hi if j_cut_hi>m else m):s_hi+=math.fsum(om_hi[max(m,j_cut_hi):j_end_hi])
 j_start_lo=m if (u_hi_i-m)<nr else (u_hi_i-nr+1);j_cut_lo=u_hi_i if u_hi_i<j_end_lo else j_end_lo;s_lo=0.0
 if j_cut_lo>j_start_lo:s_lo+=sp(r_lo_rev[nr-1-u_hi_i+j_start_lo:nr-1-u_hi_i+j_cut_lo],om_lo[j_start_lo:j_cut_lo])
 if j_end_lo>(j_cut_lo if j_cut_lo>j_start_lo else j_start_lo):s_lo+=math.fsum(om_lo[max(j_start_lo,j_cut_lo):j_end_lo])
 na=math.nextafter
 for _ in range(8):s_hi=na(s_hi,math.inf);s_lo=na(s_lo,0.0)
 return min(1.0,max(0.0,na(ans_lo+na(h_lo*max(0.0,s_lo),0.0),0.0))),min(1.0,na(ans_hi+na(h_hi*s_hi,math.inf),math.inf))
def psurv512(m_bits,B,T,mesh=512): return psurv_iv(m_bits,B,T,mesh)[1]
def verify_v43_descent(D,M,segs):
 import re
 z_ops=segs.get("desc_zero16",-1);sc_ops=segs.get("desc_scan8_u8",-1);ht_ops=segs.get("desc_sieve_hit",-1)
 m5p=segs.get("mont5_cios_prologue",-1);m5r=segs.get("mont5_cios_round",-1);m5e=segs.get("mont5_cios_epilogue",-1)
 m2=segs.get("mont2_cios",-1);bg_ops=segs.get("bingcd_gen_loop",-1);mi_ops=segs.get("modinv_gen_loop",-1);cios5_ir=m5p+5*m5r+m5e
 cg_mul=segs.get("cios_gen_mul_inner",-1);cg_red=segs.get("cios_gen_red_inner",-1)
 m1_fp=D.MONT5;a1_fp=25.0;m2_fp=3*m1_fp+5*a1_fp;s2_fp=2*m1_fp+3*a1_fp;a2_fp=50.0
 pm2=M.P_BN254-2;i2_fp=(4+(pm2.bit_length()-1)+(pm2.bit_count()-1))*m1_fp+2*a1_fp
 g2_mul_wc=254*(3*m2_fp+4*s2_fp+10*a2_fp)+254*(12*m2_fp+4*s2_fp+10*a2_fp)+i2_fp+3*m2_fp+s2_fp
 m6_fp=8*m2_fp+15*a2_fp;m12_fp=3*m6_fp+m2_fp+30*a2_fp;i6_fp=12*m2_fp+3*s2_fp+15*a2_fp+i2_fp;i12_fp=4*m6_fp+m2_fp+i6_fp+10*a2_fp
 H_num=M.P_BN254**4-M.P_BN254**2+1;H_hard=H_num//M.R_BN254;h_m12=(H_hard.bit_length()-1)+(H_hard.bit_count()-1)
 ate_ver_wc=65*(2*m12_fp+6*m2_fp+6*s2_fp+20*a2_fp)+67*(m12_fp+12*m2_fp+4*s2_fp+20*a2_fp)+((h_m12+4)*m12_fp+i12_fp+200000.0)
 g2_ver_ok=(pm2.bit_length()==254 and pm2.bit_count()==110 and H_num%M.R_BN254==0 and H_hard.bit_length()==761 and H_hard.bit_count()==389 and g2_mul_wc<=5422256.0 and g2_mul_wc<=M.G2_OUTPUT_BOUND and ate_ver_wc<=40297005.0 and ate_ver_wc<=20000.0*(2.0**11))
 segs_ok=(z_ops==18 and sc_ops==22 and ht_ops==9 and m5p==18 and m5r==48 and m5e==21 and cios5_ir==279 and cios5_ir+7<=D.MONT5
and m2==47 and bg_ops==49 and mi_ops==91 and cg_mul==11 and cg_red==11 and z_ops/256.0<=D.ZERO and sc_ops/64.0<=D.SCAN and g2_ver_ok)
 kpath=Path(__file__).resolve().parent/"kernels_v32.s"
 klines=kpath.read_text(encoding="ascii").splitlines()
 m5_macro=[];in_mac=False;fn5=[];fn2=[];fn_bg=[];fn_cg=[];fn_mi=[];lab_bg={};lab_cg={};lab_mi={};cur_fn=None;rep_n=0;rep_buf=[]
 for raw in klines:
  ln=raw.split("//",1)[0].strip()
  if not ln: continue
  if ln.startswith(".macro MONT5_ROUND"):in_mac=True; continue
  if in_mac:
   if ln.startswith(".endm"):in_mac=False
   else:m5_macro.append(ln)
   continue
  if ln=="seg_mont5_cios:":cur_fn="m5"; continue
  if ln=="seg_mont2_cios:":cur_fn="m2"; continue
  if ln=="seg_bingcd_gen:":cur_fn="bg"; continue
  if ln=="seg_cios_gen:":cur_fn="cg"; continue
  if ln=="seg_modinv_gen:":cur_fn="mi"; continue
  if cur_fn=="m5":
   if ln.startswith(".rept "):rep_n=int(ln.split()[1]);rep_buf=[]; continue
   if ln==".endr":
    for _ in range(rep_n):
     for rl in rep_buf:fn5.extend(m5_macro if rl=="MONT5_ROUND" else [rl])
    rep_n=0; continue
   if rep_n:rep_buf.append(ln); continue
   if ln.endswith(":") or ln.startswith("."): continue
   fn5.extend(m5_macro if ln=="MONT5_ROUND" else [ln])
   if ln=="ret":cur_fn=None
  elif cur_fn=="m2":
   if ln.endswith(":") or ln.startswith("."): continue
   fn2.append(ln)
   if ln=="ret":cur_fn=None
  elif cur_fn=="bg":
   if ln.endswith(":"):lab_bg[ln[:-1]]=len(fn_bg); continue
   if ln.startswith("."): continue
   fn_bg.append(ln)
   if ln=="ret":cur_fn=None
  elif cur_fn=="cg":
   if ln.endswith(":"):lab_cg[ln[:-1]]=len(fn_cg); continue
   if ln.startswith("."): continue
   fn_cg.append(ln)
   if ln=="ret":cur_fn=None
  elif cur_fn=="mi":
   if ln.endswith(":"):lab_mi[ln[:-1]]=len(fn_mi); continue
   if ln.startswith("."): continue
   fn_mi.append(ln)
   if ln=="ret":cur_fn=None
 M64=(1<<64)-1
 def _reg_r(s):
  s=s.strip().rstrip(",")
  return 32 if s=="xzr" else(31 if s=="sp" else int(s[1:]))
 def _reg_w(s):
  s=s.strip().rstrip(",")
  return 33 if s=="xzr" else(31 if s=="sp" else int(s[1:]))
 _pcache={};_OP3={"adcs":1,"adds":2,"mul":3,"umulh":4,"sbcs":5,"subs":6,"add":7,"sub":8,"adc":9,"orr":10,"csel":11};_OPB={"b":20,"b.ne":21,"b.lo":22,"b.hs":23,"b.pl":24};_OPI={"add":13,"sub":14,"lsr":15,"lsl":16,"subs":17}
 def _compile_arm64(prog,labels):
  key=id(prog)
  if key in _pcache:return _pcache[key]
  code=[]
  for ins in prog:
   op,_,rest=ins.partition(" ");rest=rest.strip()
   if op=="ret":code.append((0,0,0,0,0))
   elif op in _OPB:code.append((_OPB[op],labels[rest],0,0,0))
   elif op in("cbz","cbnz"):r,lb=rest.split(",");code.append((25 if op=="cbz" else 26,_reg_r(r),labels[lb.strip()],0,0))
   elif op in("tbz","tbnz"):r,bt,lb=rest.split(",");code.append((27 if op=="tbz" else 28,_reg_r(r),int(bt.strip().lstrip("#")),labels[lb.strip()],0))
   elif op=="cmp":a,b=rest.split(",");code.append((29,_reg_r(a),_reg_r(b),0,0))
   elif op=="mov":a,b=rest.split(",");bs=b.strip();code.append((30,_reg_w(a),int(bs.lstrip("#")),0,0) if bs.startswith("#") else (31,_reg_w(a),_reg_r(bs),0,0))
   elif op in _OPI and "#" in rest:a,b,imm=rest.split(",");code.append((_OPI[op],_reg_w(a),_reg_r(b),int(imm.strip().lstrip("#")),0))
   elif op in("stp","ldp"):m=re.match(r"(\w+),\s*(\w+),\s*\[(\w+),\s*#(\d+)\]$",rest);code.append((32,_reg_r(m.group(1)),_reg_r(m.group(2)),_reg_r(m.group(3)),int(m.group(4))) if op=="stp" else (33,_reg_w(m.group(1)),_reg_w(m.group(2)),_reg_r(m.group(3)),int(m.group(4))))
   elif op in("str","ldr"):
    is_st=(op=="str");m=re.match(r"(\w+),\s*\[(\w+),\s*(\w+),\s*lsl\s*#(\d+)\]$",rest)
    if m:code.append((34 if is_st else 12,_reg_r(m.group(1)) if is_st else _reg_w(m.group(1)),_reg_r(m.group(2)),_reg_r(m.group(3)),int(m.group(4))))
    else:
     m=re.match(r"(\w+),\s*\[(\w+)\],\s*#(\d+)$",rest)
     if m:code.append((35 if is_st else 18,_reg_r(m.group(1)) if is_st else _reg_w(m.group(1)),_reg_r(m.group(2)),int(m.group(3)),0))
     else:m=re.match(r"(\w+),\s*\[(\w+)(?:,\s*#(\d+))?\]$",rest);code.append((36 if is_st else 19,_reg_r(m.group(1)) if is_st else _reg_w(m.group(1)),_reg_r(m.group(2)),int(m.group(3) or 0),0))
   else:
    a,b,c=[t.strip() for t in rest.split(",")[:3]];oc=_OP3.get(op)
    if oc is None:raise ValueError(f"unsupported op: {ins}")
    code.append((oc,_reg_w(a),_reg_r(b),_reg_r(c),0))
  _pcache[key]=code;return code
 def _run_arm64(prog,x_init,mem,labels=None,max_ops=2000000):
  code=_compile_arm64(prog,labels or {});x=list(x_init)+[0,0];C=Z=N_fl=pc=ops=0;m64=M64
  # Fast flat array mirror for contiguous 8-byte-aligned addresses in [0, 0x9800] (100% exact sync back to mem!)
  use_flat=bool(mem) and all((a&7)==0 and 0<=a<0x9800 for a in mem) and x[31]<0x9800
  if use_flat:
   fm=[0]*4864
   for k_,v_ in mem.items():fm[k_>>3]=v_
   while ops<max_ops:
    oc,a,b,c,d=code[pc];pc+=1;ops+=1
    if oc<=12:
     if oc==1:s=x[b]+x[c]+C;C=s>>64;v=s&m64;Z=v==0;N_fl=v>>63;x[a]=v
     elif oc==2:s=x[b]+x[c];C=s>>64;v=s&m64;Z=v==0;N_fl=v>>63;x[a]=v
     elif oc==3:x[a]=(x[b]*x[c])&m64
     elif oc==4:x[a]=((x[b]*x[c])>>64)&m64
     elif oc==5:s=x[b]+(x[c]^m64)+C;C=s>>64;v=s&m64;Z=v==0;N_fl=v>>63;x[a]=v
     elif oc==6:s=x[b]+(x[c]^m64)+1;C=s>>64;v=s&m64;Z=v==0;N_fl=v>>63;x[a]=v
     elif oc==12:x[a]=fm[(x[b]+(x[c]<<d))>>3]
     elif oc==7:x[a]=(x[b]+x[c])&m64
     elif oc==8:x[a]=(x[b]-x[c])&m64
     elif oc==9:x[a]=(x[b]+x[c]+C)&m64
     elif oc==10:x[a]=x[b]|x[c]
     elif oc==11:x[a]=x[b] if C else x[c]
     else:
      for k_ in range(4864):
       v_=fm[k_]
       if v_ or (k_<<3) in mem:mem[k_<<3]=v_
      return x[:32],ops
    elif oc<=24:
     if oc==15:x[a]=(x[b]>>c)&m64
     elif oc==16:x[a]=(x[b]<<c)&m64
     elif oc==13:x[a]=(x[b]+c)&m64
     elif oc==14:x[a]=(x[b]-c)&m64
     elif oc==17:s=x[b]+(c^m64)+1;C=s>>64;v=s&m64;Z=v==0;N_fl=v>>63;x[a]=v
     elif oc==18:ad=x[b];x[a]=fm[ad>>3];x[b]=(ad+c)&m64
     elif oc==19:x[a]=fm[(x[b]+c)>>3]
     elif oc==20:pc=a
     elif oc==21:
      if not Z:pc=a
     elif oc==22:
      if not C:pc=a
     elif oc==23:
      if C:pc=a
     else:
      if not N_fl:pc=a
    else:
     if oc==34:fm[(x[b]+(x[c]<<d))>>3]=x[a]
     elif oc==35:ad=x[b];fm[ad>>3]=x[a];x[b]=(ad+c)&m64
     elif oc==36:fm[(x[b]+c)>>3]=x[a]
     elif oc==29:s=x[a]+(x[b]^m64)+1;C=s>>64;v=s&m64;Z=v==0;N_fl=v>>63
     elif oc==31:x[a]=x[b]
     elif oc==30:x[a]=b&m64
     elif oc==25:
      if x[a]==0:pc=b
     elif oc==26:
      if x[a]!=0:pc=b
     elif oc==27:
      if not((x[a]>>b)&1):pc=c
     elif oc==28:
      if(x[a]>>b)&1:pc=c
     elif oc==32:ad=(x[c]+d)>>3;fm[ad]=x[a];fm[ad+1]=x[b]
     else:ad=(x[c]+d)>>3;x[a]=fm[ad];x[b]=fm[ad+1]
   raise RuntimeError("max_ops exceeded")
  while ops<max_ops:
   oc,a,b,c,d=code[pc];pc+=1;ops+=1
   if oc<=12:
    if oc==1:s=x[b]+x[c]+C;C=s>>64;v=s&m64;Z=v==0;N_fl=(v>>63);x[a]=v
    elif oc==2:s=x[b]+x[c];C=s>>64;v=s&m64;Z=v==0;N_fl=(v>>63);x[a]=v
    elif oc==3:x[a]=(x[b]*x[c])&m64
    elif oc==4:x[a]=((x[b]*x[c])>>64)&m64
    elif oc==5:s=x[b]+(x[c]^m64)+C;C=s>>64;v=s&m64;Z=v==0;N_fl=(v>>63);x[a]=v
    elif oc==6:s=x[b]+(x[c]^m64)+1;C=s>>64;v=s&m64;Z=v==0;N_fl=(v>>63);x[a]=v
    elif oc==12:x[a]=mem[x[b]+(x[c]<<d)]
    elif oc==7:x[a]=(x[b]+x[c])&m64
    elif oc==8:x[a]=(x[b]-x[c])&m64
    elif oc==9:x[a]=(x[b]+x[c]+C)&m64
    elif oc==10:x[a]=x[b]|x[c]
    elif oc==11:x[a]=x[b] if C else x[c]
    else:return x[:32],ops
   elif oc<=24:
    if oc==15:x[a]=(x[b]>>c)&m64
    elif oc==16:x[a]=(x[b]<<c)&m64
    elif oc==13:x[a]=(x[b]+c)&m64
    elif oc==14:x[a]=(x[b]-c)&m64
    elif oc==17:s=x[b]+(c^m64)+1;C=s>>64;v=s&m64;Z=v==0;N_fl=(v>>63);x[a]=v
    elif oc==18:ad=x[b];x[a]=mem[ad];x[b]=(ad+c)&m64
    elif oc==19:x[a]=mem[x[b]+c]
    elif oc==20:pc=a
    elif oc==21:
     if not Z:pc=a
    elif oc==22:
     if not C:pc=a
    elif oc==23:
     if C:pc=a
    else:
     if not N_fl:pc=a
   else:
    if oc==34:mem[x[b]+(x[c]<<d)]=x[a]
    elif oc==35:ad=x[b];mem[ad]=x[a];x[b]=(ad+c)&m64
    elif oc==36:mem[x[b]+c]=x[a]
    elif oc==29:s=x[a]+(x[b]^m64)+1;C=s>>64;v=s&m64;Z=v==0;N_fl=(v>>63)
    elif oc==31:x[a]=x[b]
    elif oc==30:x[a]=b&m64
    elif oc==25:
     if x[a]==0:pc=b
    elif oc==26:
     if x[a]!=0:pc=b
    elif oc==27:
     if not((x[a]>>b)&1):pc=c
    elif oc==28:
     if(x[a]>>b)&1:pc=c
    elif oc==32:ad=x[c]+d;mem[ad]=x[a];mem[ad+8]=x[b]
    else:ad=x[c]+d;x[a]=mem[ad];x[b]=mem[ad+8]
  raise RuntimeError("max_ops exceeded")
 rng=random.Random(20261143)
 asm_m5_ok=(len(fn5)==279);asm_m2_ok=(len(fn2)==47);asm_bg_ok=(len(fn_bg)==49);asm_mi_ok=(len(fn_mi)==mi_ops);asm_cg_ok=True
 for ti in range(32):
  N=(rng.getrandbits(318)|(1<<318)|1);n0p=(-pow(N,-1,1<<64))&M64
  a=(N-1 if ti==0 else rng.randrange(N));b=(N-1 if ti==0 else rng.randrange(N))
  mem={};PA,PB,PO,PP,SP0=0x1000,0x2000,0x3000,0x4000,0x8000
  for j in range(5):mem[PA+8*j]=(a>>(64*j))&M64;mem[PB+8*j]=(b>>(64*j))&M64;mem[PP+8*j]=(N>>(64*j))&M64
  x0=[0]*32;x0[5],x0[6],x0[7],x0[8],x0[9],x0[31]=PA,PB,n0p,PO,PP,SP0
  for r_s in range(19,28):x0[r_s]=0xDEADBEEF0000+r_s
  x1,op5=_run_arm64(fn5,x0,mem)
  got5=sum(mem[PO+8*j]<<(64*j) for j in range(5))
  if op5!=279 or got5!=(a*b*pow(1<<320,-1,N))%N or x1[31]!=SP0 or any(x1[r_s]!=x0[r_s] for r_s in range(19,28)):asm_m5_ok=False
  N2=(1<<128)-1 if ti==0 else((1<<127)+1 if ti==1 else(rng.getrandbits(126)|(1<<126)|1))
  n0p2=(-pow(N2,-1,1<<64))&M64;a2=(N2-1 if ti==0 else rng.randrange(N2));b2=(N2-1 if ti==0 else rng.randrange(N2))
  x2=[0]*32;x2[0],x2[1],x2[2],x2[3],x2[4],x2[5],x2[6]=a2&M64,(a2>>64)&M64,b2&M64,(b2>>64)&M64,N2&M64,(N2>>64)&M64,n0p2
  x3,op2=_run_arm64(fn2,x2,{})
  if op2!=47 or (x3[0]|(x3[1]<<64))!=(a2*b2*pow(1<<128,-1,N2))%N2:asm_m2_ok=False
 for Lg in (2,5,6,8,12,20):
  nb_g=64*Lg;max_bg_ops=int(D.lehmer_gcd_ir(nb_g))
  for ti in range(3):
   ug=(rng.getrandbits(nb_g)|1);vg=(0 if ti==0 else (ug if ti==1 else rng.getrandbits(nb_g)))
   if ti==2:
    cf=(rng.getrandbits(min(60,nb_g//2))|1);ug=(ug//cf)*cf;vg=(vg//cf)*cf
    if (ug&1)==0:ug|=1
   mb={};PU,PV=0x1000,0x3000
   for j in range(Lg):mb[PU+8*j]=(ug>>(64*j))&M64;mb[PV+8*j]=(vg>>(64*j))&M64
   xb=[0]*32;xb[0],xb[1],xb[2]=PU,PV,Lg
   xb1,opb=_run_arm64(fn_bg,xb,mb,lab_bg)
   gotb=sum(mb[xb1[0]+8*j]<<(64*j) for j in range(Lg))
   if gotb!=math.gcd(ug,vg) or opb>max_bg_ops:asm_bg_ok=False
 for Lg in (6,8,12,20):
  R_inv_bit=64*Lg;max_cg_ops=int(D.mul_L(Lg))
  for ti in range(4):
   Ng=(rng.getrandbits(R_inv_bit-2)|(1<<(R_inv_bit-2))|1);n0g=(-pow(Ng,-1,1<<64))&M64
   ag=(Ng-1 if ti==0 else rng.randrange(Ng));bg=(Ng-1 if ti==0 else rng.randrange(Ng))
   mg={};PA,PB,PN,PT,PO=0x1000,0x3000,0x5000,0x7000,0x9000
   for j in range(Lg):mg[PA+8*j]=(ag>>(64*j))&M64;mg[PB+8*j]=(bg>>(64*j))&M64;mg[PN+8*j]=(Ng>>(64*j))&M64
   xg=[0]*32;xg[0],xg[1],xg[2],xg[3],xg[4],xg[5],xg[6]=PA,PB,PN,PT,Lg,n0g,PO
   _,opg=_run_arm64(fn_cg,xg,mg,lab_cg)
   gotg=sum(mg[PO+8*j]<<(64*j) for j in range(Lg))
   if gotg!=(ag*bg*pow(1<<R_inv_bit,-1,Ng))%Ng or opg>max_cg_ops:asm_cg_ok=False
 for Lg in (2,5,6,8,12,20):
  nb_g=64*Lg;max_mi_ops=int(D.modinv_ir(nb_g))
  for ti in ((0,2) if Lg==20 else (0,1,2)):
   Nm=(rng.getrandbits(nb_g-1)|(1<<(nb_g-1))|1)
   am=(rng.randrange(1,Nm) if ti<2 else ((rng.randrange(1,Nm//3)*3)%Nm or 3))
   if ti==2 and Nm%3!=0:Nm=(Nm//3)*3;Nm=Nm if (Nm&1) else Nm+3
   if am>=Nm:am%=Nm
   if am==0:am=1
   mm={};PU,PV,PX1,PX2,PN=0x1000,0x3000,0x5000,0x7000,0x9000
   for j in range(Lg):
    mm[PU+8*j]=(Nm>>(64*j))&M64;mm[PV+8*j]=(am>>(64*j))&M64
    mm[PX1+8*j]=0;mm[PX2+8*j]=(1 if j==0 else 0);mm[PN+8*j]=(Nm>>(64*j))&M64
   xm=[0]*32;xm[0],xm[1],xm[2],xm[3],xm[4],xm[5]=PU,PV,Lg,PX1,PX2,PN
   xm1,opm=_run_arm64(fn_mi,xm,mm,lab_mi)
   gotg=sum(mm[xm1[0]+8*j]<<(64*j) for j in range(Lg))
   goti=sum(mm[xm1[3]+8*j]<<(64*j) for j in range(Lg))
   eg=math.gcd(am,Nm)
   if gotg!=eg or opm>max_mi_ops or (eg==1 and (goti*am)%Nm!=1):asm_mi_ok=False
 segs_ok=(segs_ok and asm_m5_ok and asm_m2_ok and asm_bg_ok and asm_mi_ok and asm_cg_ok)
 lh=bg_ops
 data=M.load_data()["descent"];PLAN=data["plan59"];CELLS={int(k):v for k,v in data["cells"].items()}
 cfmax=max(b[1] for b in D.QBINS);cgmax=max(b[2] for b in D.QBINS);mx_nf=mx_ng=mx_q=mx_S=0.0
 for blo,bhi,sd,_,_,_,_,R,_ in PLAN:
  Sc=log2(CELLS[R]);cc=(bhi+Sc)/12.0
  mx_nf=max(mx_nf,D.DF*cc+cfmax-(blo if sd=="f" else 0.0))
  mx_ng=max(mx_ng,12.0*cc+D.LOG2U6+cgmax-(blo if sd=="g" else 0.0))
  mx_q=max(mx_q,float(bhi));mx_S=max(mx_S,Sc)
 mnb=max(ceil(mx_nf),ceil(mx_ng),ceil(data["boot_dp"]["ng_hi"]));mnb_gate=max(mnb,2560);ncrt=ceil(mnb_gate/60.0)+1
 w_norms=(ncrt/8.0)*(M.C_NORM_G+M.C_NORM_F)+60.0*ncrt*ncrt
 w_tdiv=309*(ceil(mnb_gate/64.0)*5+10);mq=D.mul_L(max(4,ceil(mx_q/64.0)))
 wc_ok=(ncrt>=44 and w_norms<D.NORMS and w_tdiv<D.TDIV and 64.0*mx_q*mq<D.PRIME and 6.0*mx_q*45.0*mq<D.IDEAL and mx_S<=82.0)
 ps_ok=True;prev_hi=1.0
 for mb in range(400,1900,100):
  lo,hi=psurv_iv(mb,59.0,500.0,512)
  if not(0.0<=lo<=hi<=1.0 and hi<=prev_hi+1e-15):ps_ok=False
  prev_hi=hi
 N_p1=(rng.randrange(1<<592,1<<593)|1);N_p2=(rng.randrange(1<<591,1<<592)|1);N=N_p1*N_p2
 def _pm(A,B):
  out=[0]*(len(A)+len(B)-1)
  for i,a in enumerate(A):
   if a:
    for j,b in enumerate(B):out[i+j]+=a*b
  return [x%N for x in out]
 def _ntt_p(arr,w,p,inv=False):
  n=len(arr);wb=pow(w,-1,p) if inv else w;kb=n.bit_length()-1
  a=[arr[int(f"{i:0{kb}b}"[::-1],2)]%p for i in range(n)];m=1
  while m<n:
   st=m;m<<=1;ws=pow(wb,n//m,p);w1=1;tw=[]
   for _ in range(st):tw.append((w1,(w1<<64)//p));w1=(w1*ws)%p
   for k0 in range(0,n,m):
    for j in range(st):
     w1,wp=tw[j];u=a[k0+j];v0=a[k0+j+st]
     v=((v0*w1-((v0*wp)>>64)*p)&M64)%p;a[k0+j]=(u+v)%p;a[k0+j+st]=(u-v)%p
  return [(x*pow(n,-1,p))%p for x in a] if inv else a
 A16=[rng.randrange(N) for _ in range(16)];B16=[rng.randrange(N) for _ in range(16)]
 S_k=ceil((2*N.bit_length()+5)/64.0);out_l=32*S_k;n_ntt=1<<ceil(log2(out_l));p1,p2,p3=M.NTT3_PRIMES
 def _pk(P):
  ar=[0]*n_ntt
  for i,c in enumerate(P):
   v=c%N;off=i*S_k
   for l in range(S_k):ar[off+l]=v&M64;v>>=64
  return ar
 ap,bp=_pk(A16),_pk(B16);inv12=pow(p1,-1,p2);inv123=pow((p1*p2)%p3,-1,p3);rs=[]
 for pi in M.NTT3_PRIMES:
  gi=2
  while pow(gi,(pi-1)//2,pi)==1 or pow(gi,(pi-1)//3,pi)==1:gi+=1
  wi=pow(gi,(pi-1)//n_ntt,pi);fa,fb=_ntt_p(ap,wi,pi),_ntt_p(bp,wi,pi)
  rs.append(_ntt_p([(fa[i]*fb[i])%pi for i in range(n_ntt)],wi,pi,True))
 cv=[0]*out_l;cry=0
 for i in range(out_l):
  r1,r2,r3=rs[0][i],rs[1][i],rs[2][i];v12=r1+(((r2-r1)%p2*inv12)%p2)*p1
  tot=v12+(((r3-v12%p3)%p3*inv123)%p3)*p1*p2+cry;cv[i]=tot&M64;cry=tot>>64
 got16=[]
 for k in range(31):
  v=0
  for l in range(S_k-1,-1,-1):v=(v<<64)|cv[k*S_k+l]
  got16.append(v%N)
 kron_ok=(cry==0 and got16==_pm(A16,B16))
 def _u_z():
  while True:
   v=rng.randrange(1,N)
   if math.gcd(v,N)==1: return v
 pb_X=[rng.randrange(1,N) for _ in range(32)];pb_Z=[_u_z() for _ in range(32)]
 pg_X=[rng.randrange(1,N) for _ in range(32)];pg_Z=[_u_z() for _ in range(32)]
 all_Z=pb_Z+pg_Z;pref_Z=[all_Z[0]];bn_muls=0
 for z_ in all_Z[1:]:pref_Z.append((pref_Z[-1]*z_)%N);bn_muls+=1
 inv_all=pow(pref_Z[-1],-1,N);inv_Z=[0]*64
 for i_ in range(63,0,-1):inv_Z[i_]=(inv_all*pref_Z[i_-1])%N;inv_all=(inv_all*all_Z[i_])%N;bn_muls+=2
 inv_Z[0]=inv_all
 bab=[(x_*iz_)%N for x_,iz_ in zip(pb_X,inv_Z[:32])];gia=[(x_*iz_)%N for x_,iz_ in zip(pg_X,inv_Z[32:])];bn_muls+=64;Fm=[1]
 bad_Z=list(all_Z);bad_Z[3]=(bad_Z[3]*N_p1)%N;bad_Z[47]=(bad_Z[47]*N_p2)%N;bpr=[bad_Z[0]]
 for z_ in bad_Z[1:]:bpr.append((bpr[-1]*z_)%N)
 lo_b,hi_b,b_gcds=0,63,0
 while lo_b<hi_b:
  md_b=(lo_b+hi_b)//2;b_gcds+=1
  if math.gcd(bpr[md_b],N)>1:hi_b=md_b
  else:lo_b=md_b+1
 bsplit_ok=(math.gcd(bpr[-1],N)==N and b_gcds<=6 and 1<math.gcd(bpr[lo_b],N)<N)
 for b in bab:Fm=_pm(Fm,[-b%N,1])
 def _bt(pts):
  if len(pts)==1: return ([-pts[0]%N,1],None,None,1)
  md=len(pts)//2;L=_bt(pts[:md]);R=_bt(pts[md:]); return (_pm(L[0],R[0]),L,R,len(pts))
 rt_g=_bt(gia);Gm=rt_g[0];Fr=[(Fm[i]-Gm[i])%N for i in range(32)];rg=Gm[::-1];invg=[1];mk=1
 while mk<32:
  mk2=min(2*mk,32);hi_=_pm(rg[:mk2],invg)[:mk2];t2=[(-x)%N for x in hi_];t2[0]=(2-hi_[0])%N
  invg=_pm(invg,t2)[:mk2];mk=mk2
 M0=_pm(Fr[::-1],invg)[:32]
 def _desc(Mc,nd):
  _,L,R,sz=nd
  if L is None: return [Mc[0]%N]
  return _desc(_pm(Mc,R[0][::-1])[R[3]:sz],L)+_desc(_pm(Mc,L[0][::-1])[L[3]:sz],R)
 bls_ok=(bsplit_ok and bn_muls==4*64-3 and bn_muls<=4*64 and all(v==math.prod((g-b)%N for b in bab)%N for g,v in zip(gia,_desc(M0,rt_g))))
 max_J=max(D._stage_ops(w,w,w)[3] for w in D.WHEELS);L_mx=ceil(max(mnb,1225)/64.0)
 S_mx=ceil((2*max(mnb,1225)+20)/64.0);n_ntt_mx=1<<ceil(log2(2*(max_J+1)*S_mx))
 ram_b=6*n_ntt_mx*8+4*max_J*L_mx*8*4+4*max_J*L_mx*8
 j_ok=(max_J==46080 and all(D._stage_ops(b1,b2,wh)[3]<=max_J for _,b1,b2,wh,_ in D._ECM_BEST_C.values() if wh>0))
 ok=(segs_ok and wc_ok and ps_ok and kron_ok and bls_ok and j_ok and n_ntt_mx<(1<<54) and ram_b<(1<<30))
 return dict(ok=ok,cios5=cios5_ir,mont2=m2,bingcd_gen=lh,psurv_iv_ok=ps_ok,kron_ntt3_ok=kron_ok,bls_mpe_ok=bls_ok,ecm_ram_mib=round(ram_b/(1<<20),2))
