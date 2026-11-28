"""v32_verify.py -- BN254 C6 verifier."""
import math,random,re
from operator import add as _oadd, mul as _omul
from itertools import repeat as _orep
from fractions import Fraction as Fr
from math import log,log2,exp,sqrt,ceil,floor
from pathlib import Path
LN2=math.log(2.0)
L3=log2(3.0)
M64=(1<<64)-1
S_U8=0.35
B1_DELTA_HI=29.5
C_HIT_DELTA_BASE=2.82
_RET2D_CACHE={}
def get_ret2d_u8(M,pmin=24,B=61.50,odd=2,scale=S_U8,bstep=1.0,dstep=0.50,sl_step=0.25,max_bits=488):
 """Discretized 2D generating-series scoring model (H1/H2) for uint8_t STADDB (scale=0.35); see KERNELS.md."""
 key=(pmin,B,odd,scale,bstep,dstep,sl_step,max_bits)
 if key in _RET2D_CACHE:
  return _RET2D_CACHE[key]
 D_LO,D_HI=-40.0,40.0
 L=ceil(max_bits/bstep)
 D0=int(round(D_LO/dstep));D1=int(round(D_HI/dstep))
 ND=D1-D0+1
 tilt=0.04
 ret_1d_rg=M.Retention(B,B,prime_min=1<<pmin,explicit_bits=14,step=bstep,max_bits=max_bits,tilt=tilt,odd=odd)
 ret_1d_sf=M.Retention(float(max_bits),B,prime_min=1<<pmin,explicit_bits=14,step=bstep,max_bits=max_bits+1,tilt=0.005,odd=odd)
 def _build_sieve_2d(tilt_v,is_sf):
  w2=[dict() for _ in range(L+1)]
  tot=0.0
  k_lo=int(floor(scale*pmin+0.5));k_hi=int(ceil(scale*B+0.5))+1
  pts=sorted(set([float(pmin),float(B)]+[(k-0.5)/scale for k in range(k_lo,k_hi+1) if pmin<(k-0.5)/scale<B]))
  for a,c in zip(pts[:-1],pts[1:]):
   wp=int(floor(scale*0.5*(a+c)+0.5))
   n_sub=max(4,int(ceil((c-a)/0.05)))
   dw=(c-a)/n_sub
   for i in range(n_sub):
    b=a+(i+0.5)*dw
    bi=round(b/bstep)
    if bi>L: continue
    delta_bin=bi*bstep-wp/scale
    di=ceil(delta_bin/dstep) if is_sf else floor(delta_bin/dstep)
    amt1=exp(tilt_v*bi*bstep)*dw/b
    tot+=amt1
    dm=w2[bi];dm[di]=dm.get(di,0.0)+amt1
    for kp in range(2,5):
     idxk=kp*bi
     if idxk>L: break
     amtk=exp(kp*tilt_v*bi*bstep-(kp-1)*LN2*b)*dw/(kp*b)
     if amtk<=1e-25: break
     dmk=w2[idxk]
     if is_sf:
      sgn=-1.0 if(kp%2==0) else 1.0
      tot+=sgn*amtk
      dik=max(D0,min(D1,kp*di))
      dmk[dik]=dmk.get(dik,0.0)+sgn*amtk
     else:
      tot+=amtk
      for jp in range(1,kp+1):
       coef=math.comb(kp,jp)*(1.0 if(jp%2==1) else-1.0)
       djk=max(D0,min(D1,floor((kp*bi*bstep-jp*wp/scale)/dstep)))
       dmk[djk]=dmk.get(djk,0.0)+coef*amtk
  nz=[(j,tuple((dd,j*val) for dd,val in w2[j].items())) for j in range(1,L+1) if w2[j]]
  out=[[0.0]*ND for _ in range(L+1)]
  out[0][-D0]=exp(-tot)
  rng_d=[(-D0,-D0)]+[(ND,-1)]*L
  for n in range(nz[0][0],L+1):
   row=out[n];inv_n=1.0/n
   d_min=ND;d_max=-1
   for j,items in nz:
    if j>n: break
    prev=out[n-j];p_lo,p_hi=rng_d[n-j]
    if p_lo>p_hi: continue
    for dd,jw in items:
     lo_in=p_lo if p_lo+dd>=0 else-dd
     hi_in=p_hi if p_hi+dd<ND else ND-1-dd
     if lo_in<=hi_in:
      a0=lo_in+dd;b0=hi_in+1+dd
      row[a0:b0]=map(_oadd,row[a0:b0],map(_omul,_orep(jw),prev[lo_in:hi_in+1]))
      if a0<d_min:d_min=a0
      if b0-1>d_max:d_max=b0-1
     if not is_sf and p_lo+dd<0:
      s_lo=sum(prev[p_lo:min(p_hi+1,-dd)])
      if s_lo:
       row[0]+=jw*s_lo;d_min=0
       if 0>d_max:d_max=0
     if p_hi+dd>=ND:
      s_hi=sum(prev[max(p_lo,ND-dd):p_hi+1])
      if s_hi:
       row[ND-1]+=jw*s_hi
       if ND-1<d_min:d_min=ND-1
       d_max=ND-1
   if d_min<=d_max:
    for idx in range(d_min,d_max+1):row[idx]*=inv_n
    while d_min<=d_max and abs(row[d_min])<1e-26:d_min+=1
    while d_max>=d_min and abs(row[d_max])<1e-26:d_max-=1
    rng_d[n]=(d_min,d_max)
  return out,tot
 sv_rg,tot_rg=_build_sieve_2d(tilt,True)
 sv_sf,tot_sf=_build_sieve_2d(0.005,False)
 sel_norm=exp(tot_rg-ret_1d_rg.normalizers["sieve_geo"])
 sm_rg=ret_1d_rg.co("small");sm_sf=ret_1d_sf.co("small");lg_sf=ret_1d_sf.co("large_geo")
 m_sf=[0.0]*(L+1);tsm=L
 while tsm>0 and sm_sf[tsm]<1e-22:tsm-=1
 for l,v in enumerate(lg_sf):
  if v>1e-22:
   for t in range(min(L-l,tsm)+1):m_sf[l+t]+=v*sm_sf[t]
 p_rg=[[0.0]*ND for _ in range(L+1)];p_sf=[[0.0]*ND for _ in range(L+1)]
 for y in range(L+1):
  acc=0.0
  for i,v in enumerate(sv_rg[y]):acc+=v;p_rg[y][i]=acc
  acc=0.0
  for i,v in enumerate(sv_sf[y]):acc+=v;p_sf[y][i]=acc
 m_margin=0.5/scale+0.016+0.5*bstep
 den_rg=ret_1d_rg.denominator;den_sf=ret_1d_sf.denominator
 rg_c={};sf_c={};s0_c={};ub_c={}
 def get_rg(n_int,sl_c):
  r_idx=int(floor(sl_c/sl_step))
  if r_idx<0: return 0.0
  k=(n_int,r_idx);v=rg_c.get(k)
  if v is not None: return v
  ni=round(n_int/bstep)
  lim=r_idx*sl_step-m_margin
  x_max=min(ni,ceil((lim-D_LO)/bstep)+1)
  acc=0.0
  for x in range(max(0,x_max+1)):
   wx=sm_rg[x]
   if wx>1e-22:
    dl=floor((lim-x*bstep)/dstep)-D0
    if dl>=0:acc+=wx*p_rg[ni-x][dl if dl<ND else ND-1]
  v=max(0.0,min(1.0,acc*sel_norm/den_rg[ni]))
  rg_c[k]=v; return v
 def get_sf_s0(n_int,sl_c):
  r_idx=int(ceil(sl_c/sl_step))
  k=(n_int,r_idx);v=sf_c.get(k)
  if v is not None: return v,s0_c[k]
  ni=round(n_int/bstep)
  if ni not in ub_c:
   ub_c[ni]=ret_1d_sf.two(ni,ni,0,None,False)/den_sf[ni]
  lim=r_idx*sl_step+m_margin
  x_max=min(ni,ceil((lim-D_LO)/bstep)+1)
  acc=0.0
  for x in range(max(0,x_max+1)):
   wx=m_sf[x]
   if wx>1e-22:
    dl=floor((lim-x*bstep)/dstep)-D0
    if dl>=0:acc+=wx*p_sf[ni-x][dl if dl<ND else ND-1]
  sv=acc/den_sf[ni]
  s0=ub_c[ni] if ub_c[ni]<sv else sv
  sf_c[k]=sv;s0_c[k]=s0; return sv,s0
 out=dict(get_rg=get_rg,get_sf_s0=get_sf_s0,ret_1d_rg=ret_1d_rg,ret_1d_sf=ret_1d_sf,
scale=scale,bstep=bstep,dstep=dstep,sl_step=sl_step,m_margin=m_margin,
sv_rg_sum_350=sum(sv_rg[round(350/bstep)])if 350<=max_bits else 0.0,sv_rg_ref_350=ret_1d_rg.co("sieve_sf")[round(350/bstep)]if 350<=max_bits else 0.0,
sv_sf_sum_350=sum(sv_sf[round(350/bstep)])if 350<=max_bits else 0.0,sv_sf_ref_350=ret_1d_sf.co("sieve_geo")[round(350/bstep)]if 350<=max_bits else 0.0)
 _RET2D_CACHE[key]=out
 return out
_C23_CACHE=None
def get_exact_23_classes():
 global _C23_CACHE
 if _C23_CACHE is not None:
  return _C23_CACHE
 K_TAIL=6
 p_d={(0,0):Fr(1,2)}
 rem=Fr(1,2)
 for k in range(1,K_TAIL):
  pr=Fr(1,3**k)
  p_d[(2,k)]=pr
  rem-=pr
 p_d[(2,K_TAIL)]=rem
 assert sum(p_d.values())==1
 def _conv(d1,d2):
  out={}
  for(f1,g1),w1 in d1.items():
   for(f2,g2),w2 in d2.items():
    k=(f1+f2,g1+g2)
    out[k]=out.get(k,Fr(0))+w1*w2
  return out
 p3_d=_conv(_conv(p_d,p_d),p_d)
 q3_d={(0,0):Fr(27,28),(6,0):Fr(1,28)}
 full_3=_conv(p3_d,q3_d)
 assert sum(full_3.values())==1
 c3_cap={}
 for(f3,g3),w in full_3.items():
  k=(min(f3,12),min(g3,12))
  c3_cap[k]=c3_cap.get(k,Fr(0))+w
 CLASSES_2=((Fr(7,9),0,0),(Fr(1,9),6,0),(Fr(7,72),0,6),(Fr(1,72),0,9))
 exact_ev2_f,exact_ev2_g=Fr(2,3),Fr(17,24)
 exact_ev3_f,exact_ev3_g=Fr(45,14),Fr(9,4)
 ex_f_bits=float(exact_ev2_f)+float(exact_ev3_f)*L3
 ex_g_bits=float(exact_ev2_g)+float(exact_ev3_g)*L3
 classes=[];w_zero=0.0
 tot_w=Fr(0);mean_f=0.0;mean_g=0.0
 for w2,v2f,v2g in CLASSES_2:
  for(f3,g3),w3 in sorted(c3_cap.items()):
   w=w2*w3
   tot_w+=w
   wf=float(w)
   if(v2f==0 and v2g==0 and f3==0 and g3==0) or wf<1e-5:
    w_zero+=wf
   else:
    vf_b=v2f+f3*L3;vg_b=v2g+g3*L3
    mean_f+=wf*vf_b;mean_g+=wf*vg_b
    classes.append((wf,vf_b,vg_b))
 classes.insert(0,(w_zero,0.0,0.0))
 u_vf=sorted({c[1] for c in classes})
 u_vg=sorted({c[2] for c in classes})
 vf_i={v:i for i,v in enumerate(u_vf)}
 vg_i={v:i for i,v in enumerate(u_vg)}
 sparse_fg=tuple((w,vf_i[vf],vg_i[vg]) for w,vf,vg in classes)
 assert tot_w==1
 assert mean_f<=ex_f_bits+1e-12 and mean_g<=ex_g_bits+1e-12
 _C23_CACHE=dict(
classes=tuple(classes),
u_vf=tuple(u_vf),u_vg=tuple(u_vg),sparse_fg=sparse_fg,
ex_f_bits=ex_f_bits,
ex_g_bits=ex_g_bits,
trunc_f_bits=mean_f,
trunc_g_bits=mean_g,
n_classes=len(classes),
)
 return _C23_CACHE
def verify_3adic_exact(M,seed=20261123):
 c23=get_exact_23_classes()
 P_MAPS=(
(1,473,344,58,71,434),
(1,642,473,17,681,58),
(1,344,642,434,199,17),
)
 if M._OHM_CACHE is None:
  M.is_ideal_primitive_oh([1,0,0,0,0,0],[0,1,0,0,0,0])
 ohm=M._OHM_CACHE
 hom_errs=0
 for mp in P_MAPS:
  for i in range(6):
   for j in range(6):
    lhs=(mp[i]*mp[j])%729
    rhs=sum(ohm[i][j][k]*mp[k] for k in range(6))%729
    if lhs!=rhs:
     hom_errs+=1
 distinct_ok=(len({tuple(x%3 for x in mp) for mp in P_MAPS})==3)
 from fractions import Fraction as F
 ZK_f=[[F(c,3) for c in row] for row in M.ZK_NUM]
 aug=[ZK_f[i][:]+[F(1 if i==j else 0) for j in range(6)] for i in range(6)]
 for col in range(6):
  pr=next(r for r in range(col,6) if aug[r][col]!=0)
  aug[col],aug[pr]=aug[pr],aug[col]
  pv=aug[col][col];aug[col]=[x/pv for x in aug[col]]
  for r in range(6):
   if r!=col and aug[r][col]:
    fac=aug[r][col];aug[r]=[aug[r][c]-fac*aug[col][c] for c in range(12)]
 ZK_inv=[row[6:] for row in aug]
 def _to_zk(poly):
  co=[sum(F(poly[m])*ZK_inv[m][k] for m in range(6)) for k in range(6)]
  assert all(x.denominator==1 for x in co)
  return[x.numerator for x in co]
 f_zk={d:_to_zk(M.F_REC_COEFFS[d]) for d in(8,6,4,2,0)}
 gl_zk=_to_zk(M.G_LEAD_COEFFS)
 gc_zk=[1,0,0,0,0,0]
 coeff_errs=0
 for mp in P_MAPS:
  fc=[sum(f_zk[d][k]*mp[k] for k in range(6))%729 for d in(8,6,4,2,0)]
  gc=[sum(gl_zk[k]*mp[k] for k in range(6))%729,sum(gc_zk[k]*mp[k] for k in range(6))%729]
  if fc!=[1,374,615,650,286] or gc!=[68,1]:
   coeff_errs+=1
 rng=random.Random(seed)
 v3=lambda x:6 if x%729==0 else(next(k for k in range(6) if(x//(3**k))%3!=0))
 law_errs=0;tested=0
 while tested<1500:
  a=[rng.randrange(729) for _ in range(6)]
  b=[rng.randrange(729) for _ in range(6)]
  ais=[sum(a[k]*mp[k] for k in range(6))%729 for mp in P_MAPS]
  bis=[sum(b[k]*mp[k] for k in range(6))%729 for mp in P_MAPS]
  if any(ais[i]%3==0 and bis[i]%3==0 for i in range(3)):
   continue
  tested+=1
  for i in range(3):
   ai,bi=ais[i],bis[i]
   a2,b2=(ai*ai)%729,(bi*bi)%729
   fv=(pow(a2,4,729)+374*pow(a2,3,729)*b2+615*a2*a2*b2*b2+650*a2*pow(b2,3,729)+286*pow(b2,4,729))%729
   gv=(68*a2+b2)%729
   vf_obs,vg_obs=v3(fv),v3(gv)
   if ai%3==0 or bi%3==0:
    if vf_obs!=0 or vg_obs!=0:
     law_errs+=1
   else:
    if vf_obs!=2 or vg_obs<1:
     law_errs+=1
 ok=(hom_errs==0 and distinct_ok and coeff_errs==0 and law_errs==0
and c23["trunc_f_bits"]<=c23["ex_f_bits"]+1e-12
and c23["trunc_g_bits"]<=c23["ex_g_bits"]+1e-12)
 return dict(ok=ok,hom_errs=hom_errs,distinct_3_primes_ok=distinct_ok,coeff_mod729_errs=coeff_errs,
law_errs=law_errs,tested=tested,n_classes=c23["n_classes"],
ex_f_bits=round(c23["ex_f_bits"],5),ex_g_bits=round(c23["ex_g_bits"],5),
trunc_f_bits=round(c23["trunc_f_bits"],5),trunc_g_bits=round(c23["trunc_g_bits"],5))
def split_ntt3_parts(M,R,lf,lg):
 Q=float(1<<46)
 W_b=ceil(R/Q)
 out={}
 tot=0.0
 max_ram=0.0
 for side,cap,lb,nb,lleaf in(("f",M.KH_F_CAP,lf,448,7),("g",M.KH_G_CAP,lg,512,8)):
  pb=M.pbits(0,lb,cap)
  SP=ceil(pb/64.0)
  np_=1.25506*(2.0**(0.5*lb))/(0.5*lb*LN2)+cap*1.25506*(2.0**lb)/(lb*LN2)
  h=ceil(log2(np_))
  K_max=nb//24
  M_tot=(lleaf+K_max)*Q
  L_bt=ceil(log2(M_tot))
  c_root=3.0*(M.ntt3_F(SP)+M.ntt3_F(M_tot))
  c_ptrt=L_bt*(M.ntt3_F(M_tot)+M.mul2_F(M_tot)+16.0*M_tot)
  c_gcd=(K_max*Q)*(2.0*(2.0*nb*(12.0*lleaf+16.0))+20000.0)
  c_route=(K_max*Q)*100.0
  c_ptree=h*(M.ntt3_F(SP)+16.0*SP)
  c_small=R*600000.0
  c_ideal=R*35.0*(2.0**19)
  side_ir=W_b*h*(c_root+c_ptrt+c_gcd+c_route)+c_ptree+c_small+c_ideal
  out[side]=side_ir
  tot+=side_ir
  ks_=max(2,ceil(log2(SP)));n2_=float(1<<ks_);n3_=3.0*float(1<<(ks_-2))
  N_sp=n3_ if(ks_>=3 and n3_>=SP) else n2_
  max_ram=max(max_ram,28.0*N_sp+24.0*SP+32.0*M_tot+8.0*Q)
 out["total_ir"]=tot
 out["ram_log2"]=log2(max_ram)
 return out
def b1_delta_save_per_cell(M,raw_mult,p_min=24.0,p_hi=B1_DELTA_HI,c_base=C_HIT_DELTA_BASE,c_hit_lo=None):
 c_ref=M.C_HIT_LO if c_hit_lo is None else c_hit_lo
 w_k=[math.comb(3,k)*((2.0/9.0)**k)*((2.0/3.0)**(3-k))/(512.0/729.0) for k in range(4)]
 tot_save_log=0.0
 steps=int(round((p_hi-p_min)*4))
 dp=(p_hi-p_min)/steps
 for i in range(steps):
  pb=p_min+(i+0.5)*dp
  p_hi_cell=2.0**(pb+0.5*dp)
  p_ov=sum(w_k[k]*exp(-(2.0**32)*((2.0/3.0)**k)/p_hi_cell) for k in range(4))
  eff_hit=c_base+24.0*p_ov
  if eff_hit<c_ref:
   tot_save_log+=(c_ref-eff_hit)*log((pb+0.5*dp)/(pb-0.5*dp))
 return M._WKF*M.B1_PREF_COV*raw_mult*tot_save_log
def verify_v32_all(M,kernels_path,seed=20261123):
 rng=random.Random(seed)
 t3=verify_3adic_exact(M,seed=seed)
 ps=M._primes_upto(200)
 def _bpt(arr):
  if len(arr)==1:
   return(arr[0],None,None)
  mid=len(arr)//2
  L,Rt=_bpt(arr[:mid]),_bpt(arr[mid:])
  return(L[0]*Rt[0],L,Rt)
 pt=_bpt(ps)
 nums=[]
 for _ in range(64):
  v=1
  for _ in range(rng.randint(3,9)):
   v*=rng.choice(ps)
  nums.append(v)
 def _rec(node,frags):
  if not frags:
   return[]
  prod,L,Rt=node
  if L is None:
   out=[]
   for idx,val in frags:
    e=0;x=val
    while x%prod==0:
     x//=prod;e+=1
    assert x==1 and e>=1
    out.append((idx,prod,e))
   return out
  P_L=L[0];fl,fr=[],[]
  for idx,val in frags:
   g=math.gcd(val,P_L%val)
   if g==1:
    fr.append((idx,val))
   elif g==val:
    fl.append((idx,val))
   else:
    left=g;right=val//g
    while True:
     g2=math.gcd(right,g)
     if g2==1:
      break
     left*=g2;right//=g2
    fl.append((idx,left))
    if right>1:
     fr.append((idx,right))
  return _rec(L,fl)+_rec(Rt,fr)
 got=_rec(pt,list(enumerate(nums)))
 by_i={}
 for idx,p_v,e_v in got:
  by_i.setdefault(idx,[]).append((p_v,e_v))
 ntt3_split_errs=0
 for idx,v in enumerate(nums):
  pr=1
  for p_v,e_v in by_i[idx]:
   pr*=p_v**e_v
  if pr!=v:
   ntt3_split_errs+=1
 r=M.R_BN254;chi=pow(M.P_BN254,10,r)
 mu_r=(1<<317)//r
 chi_R=(chi*(1<<256))%r
 r_inv=(-pow(r,-1,1<<64))&M64
 def _mont4(a,b):
  t=0
  for i in range(4):
   ai=(a>>(64*i))&M64
   t+=ai*b
   u0=(t*r_inv)&M64
   t=(t+u0*r)>>64
  if t>=r:
   t-=r
  return t
 def _bar5(A):
  neg=(A<0);U=-A if neg else A
  top=(((U>>256)<<3)|((U>>192)>>61))&M64
  qh=(top*mu_r)>>64
  rem=U-qh*r
  if rem>=r:
   rem-=r
  if rem>=r:
   rem-=r
  assert 0<=rem<r
  return(r-rem) if(neg and rem) else rem
 ctail_errs=0
 for _ in range(200):
  B0=rng.randrange(-(1<<313),1<<313)
  B1=rng.randrange(-(1<<313),1<<313)
  B2=rng.randrange(-(1<<313),1<<313)
  got_c=(_bar5(B0-B2)+_mont4(_bar5(B1+B2),chi_R))%r
  ref_c=((B0-B2)+chi*(B1+B2))%r
  if got_c!=ref_c:
   ctail_errs+=1
 ctail_ops=43+112+154+14
 import stor_bound as _SB_v38
 _sb38=_SB_v38.main();_SB1_38=_SB_v38.SB(_SB_v38.QB_MAX)
 _na_p38=int(_SB_v38.cnt_a_up(_SB_v38.QA_MAX+_SB_v38.DELTA))+1
 _act_ks=list(range(_SB_v38.K_MIN_STRIDE,_sb38["kappa_max"]+1))
 _F38=_SB_v38.F
 _om6_dn=_SB_v38.dn(_F38(314159265,100000000)**3/6)
 _isd6_dn=_SB_v38.dn(1/_SB_v38.sqrt_up(_SB_v38.det(_SB_v38.G)))
 _r_lo=_SB_v38.sqrt_dn(_SB_v38.QA_MAX*_SB_v38.G_SC*(_F38(56880,100)**2))-_SB_v38.E6.h
 _na_lo38=int((_om6_dn*(_r_lo**6)*_isd6_dn)/6)-1
 _b2_pref_mult=sum(min(1.0,(1<<_k)/float(_na_lo38)) for _k in _act_ks)
 _stx={};_cs=None;_nbs=("x16_nobfi_pos","x16_nobfi_bnd1_hi","x16_nobfi_bnd1_lo","x16_nobfi_bnd2")
 _nbs_all={s+sf for s in _nbs for sf in("_pre","","_tail")}
 for _rl in open(kernels_path,encoding="ascii"):
  _st=_rl.strip()
  if _st.startswith("//@seg "):_cs=_st.split()[1];_stx[_cs]=[]; continue
  if _st.startswith("//@end"):_cs=None; continue
  if _cs in _nbs_all:
   _ln=_rl.split("//",1)[0].strip()
   if _ln:_stx[_cs].append(_ln)
 _nb_c={}
 def _ex_nb(sn,ents,bp,wa,r24=0,r25=0,r26=0):
  cc=_nb_c.get(sn)
  if cc is None:
   rl=_stx.get(sn+"_pre",[])+_stx.get(sn,[])+_stx.get(sn+"_tail",[])
   el,i=[],0
   while i<len(rl):
    if rl[i].startswith(".rept "):
     rc=int(rl[i].split()[1]);i+=1;bk=[]
     while i<len(rl) and rl[i]!=".endr":bk.append(rl[i]);i+=1
     el.extend(bk*rc)
    elif rl[i]!=".endr":el.append(rl[i])
    i+=1
   cc=[]
   for idx_,ins in enumerate(el):
    if ins.endswith(":"):cc.append((None,));continue
    tk=[t for t in ins.replace("["," ").replace("]"," ").replace(","," ").split() if t];op=tk[0]
    if op in("tbz","b.lo","b.hi","b.hs"):
     tg=tk[3] if op=="tbz" else tk[1]
     jt=next(j for j in(range(idx_+1,len(el)) if tg.endswith("f") else range(len(el))) if el[j]==(tg[:-1] if tg.endswith("f") else tg)+":")
     cc.append((op,tk[1],int(tk[2][1:]),jt) if op=="tbz" else (op,jt))
    else:cc.append((op,tk))
   _nb_c[sn]=cc
  R={f"x{k}":0 for k in range(32)}
  R.update(x10=bp,x19=wa&255,x22=0,x23=len(ents)*8,x24=r24,x25=r25,x26=r26)
  fc=fz=False;mm={};pc=ops=0;nc=len(cc)
  rv=lambda s:R["x"+s[1:]]&(0xFFFFFFFF if s[0]=="w" else M64)
  def wv(s,v):R["x"+s[1:]]=v&(0xFFFFFFFF if s[0]=="w" else M64)
  while pc<nc:
   it=cc[pc];op=it[0]
   if op is None:pc+=1;continue
   ops+=1
   if op=="tbz":
    if((rv(it[1])>>it[2])&1)==0:pc=it[3];continue
   elif op in("b.lo","b.hi","b.hs"):
    if((not fc) if op=="b.lo" else(fc if op=="b.hs" else(fc and not fz))):pc=it[1];continue
   elif op=="ret":break
   else:
    tk=it[1]
    if op=="ldp" and tk[1:5]==["x0","x1","x22","#16"]:
     ix=R["x22"]>>3;R["x0"],R["x1"]=ents[ix],ents[ix+1];R["x22"]+=16
    elif op=="ldr" and tk[1:4]==["x0","x22","#8"]:
     R["x0"]=ents[R["x22"]>>3];R["x22"]+=8
    elif op=="cmp" and len(tk)==3:
     a,b=rv(tk[1]),rv(tk[2]);fc=(a>=b);fz=(a==b)
    elif op=="lsr" and len(tk)==4:
     if tk[3]=="#48":raise ValueError(tk)
     wv(tk[1],rv(tk[2])>>(int(tk[3][1:]) if tk[3][0]=="#" else(rv(tk[3])&63)))
    elif op=="bic" and len(tk)==4:wv(tk[1],rv(tk[2])&~int(tk[3][1:]))
    elif op=="sub":wv(tk[1],rv(tk[2])-rv(tk[3]))
    elif op=="add":wv(tk[1],rv(tk[2])+(rv(tk[3])<<int(tk[5][1:]) if len(tk)==6 and tk[4]=="lsl" else rv(tk[3])))
    elif op=="staddb":ad=rv(tk[2]);mm[ad]=(mm.get(ad,0)+(rv(tk[1])&255))&255
    else:raise ValueError(tk)
   pc+=1
  return mm,ops
 _no48=all("#48" not in "\n".join(_stx.get(k,[])) for k in _nbs_all)
 bfi_b2_errs=0 if(_no48 and _act_ks==list(range(41,49)) and _b2_pref_mult<=M.B2_PREFIX_MULT and M.C_LIST_B2>=29.70+2.0*(_b2_pref_mult-1.0)*3.0) else 1
 bfi_b2_cases=0
 for k_s in _act_ks:
  stride=1<<k_s
  for W_t in(8,16,32):
   for sgn_v in(1,-1):
    bp=rng.randrange(1<<40,(1<<63)-W_t*stride-16);recs=[]
    for _ in range(16):
     jo=rng.randrange(W_t);ai=rng.randrange(1,stride)
     js=jo if sgn_v>0 else(W_t-1-jo);recs.append((js,ai,(js<<k_s)|ai))
    ed=(k_s==41 and W_t==8 and sgn_v==1)
    if ed:recs[0],recs[1],recs[2]=(1,1,(1<<41)|1),(0,(1<<41)-1,(1<<41)-1),(1,0,1<<41)
    j_max,j_min=(0,1) if ed else(rng.randrange(W_t),rng.randrange(W_t))
    j_lo=rng.randrange(W_t);j_hi=rng.randrange(j_lo,W_t);ents=[r[2] for r in recs]
    (mp,op),(mh,oh),(ml,ol),(mb,ob)=(
_ex_nb("x16_nobfi_pos",ents,bp,13),
_ex_nb("x16_nobfi_bnd1_hi",ents,bp,13,r25=((j_max+1)<<k_s)-1),
_ex_nb("x16_nobfi_bnd1_lo",ents,bp,13,r25=(j_min<<k_s)),
_ex_nb("x16_nobfi_bnd2",ents,bp,13,r24=k_s,r25=j_lo,r26=j_hi-j_lo))
    rp,rh,rl,rb={},{},{},{}
    for js,ai,ent in recs:
     bfi_b2_cases+=1;ad=bp+js*stride+ai
     if ent>=(1<<53) or bp+ent!=ad or ad>=(1<<64):bfi_b2_errs+=1
     rp[ad]=(rp.get(ad,0)+13)&255
     if js<=j_max:rh[ad]=(rh.get(ad,0)+13)&255
     if js>=j_min:rl[ad]=(rl.get(ad,0)+13)&255
     if j_lo<=js<=j_hi:rb[ad]=(rb.get(ad,0)+13)&255
    if mp!=rp or mh!=rh or ml!=rl or mb!=rb:bfi_b2_errs+=1
 for k_s in range(41,49):
  stride=1<<k_s;bp=1<<50;j_max,j_min,j_lo,j_hi=3,2,1,5
  for n_t in(0,1,7,15,16,17,19,31,33):
   tr=[(rng.randrange(8),rng.randrange(1,stride)) for _ in range(n_t)]
   te=[(js<<k_s)|ai for js,ai in tr]
   for sn,kw,cond,rate in(
("x16_nobfi_pos",{},lambda j: True,42.0/16.0),
("x16_nobfi_bnd1_hi",dict(r25=((j_max+1)<<k_s)-1),lambda j:j<=j_max,74.0/16.0),
("x16_nobfi_bnd1_lo",dict(r25=j_min<<k_s),lambda j:j>=j_min,74.0/16.0),
("x16_nobfi_bnd2",dict(r24=k_s,r25=j_lo,r26=j_hi-j_lo),lambda j:j_lo<=j<=j_hi,106.0/16.0)):
    mg,og=_ex_nb(sn,te,bp,13,**kw);rf={}
    for js,ai in tr:
     if cond(js):ad=bp+js*stride+ai;rf[ad]=(rf.get(ad,0)+13)&255
    if mg!=rf or og>rate*n_t+11.55:bfi_b2_errs+=1
 _cp_ins=[];_cp_lbl={};_cp_on=False;_cp_buf=None
 for _ln in open(kernels_path,encoding="ascii"):
  _s=_ln.strip()
  if _s in("//@seg x16_bfi_copy_page","//@seg x16_bfi_copy"):_cp_on=True; continue
  if _s=="//@end":_cp_on=False; continue
  if not _cp_on or not _s or _s.startswith("//"): continue
  if _s.startswith(".rept"):_cp_rep=int(_s.split()[1]);_cp_buf=[]; continue
  if _s.startswith(".endr"):
   for _ in range(_cp_rep):_cp_ins.extend(_cp_buf)
   _cp_buf=None; continue
  if _s.endswith(":"):_cp_lbl[_s[:-1]]=len(_cp_ins); continue
  if "_cp_buf" in locals() and _cp_buf is not None:_cp_buf.append(_s)
  else:_cp_ins.append(_s)
 for n_c in(0,1,2,3,15,16,127,128,129,255,256,2045,2046):
  src_b=1<<40;dst_b=16;src_m={src_b+8*i:(0xCAFE0000+i) for i in range(n_c)};dst_m={}
  R={11:n_c,15:src_b,25:dst_b};pc=0;ops=0;c_lt=False
  while pc<len(_cp_ins):
   st=_cp_ins[pc];ops+=1;pc+=1
   if st.startswith("cbz x11,"):
    if R[11]==0:pc=_cp_lbl[st.split(",")[1].strip()]
   elif st=="add x16, x15, x11, lsl #3":R[16]=R[15]+(R[11]<<3)
   elif st.startswith("tbz x11, #0,"):
    if(R[11]&1)==0:pc=_cp_lbl[st.split(",")[2].strip()]
   elif st=="ldr x0, [x15], #8":R[0]=src_m[R[15]];R[15]+=8
   elif st=="str x0, [x25], #8":dst_m[R[25]]=R[0];R[25]+=8
   elif st=="ubfx x12, x11, #1, #6":R[12]=(R[11]>>1)&63
   elif st=="adr x13, copy_end":R[13]=_cp_lbl["copy_end"]*4
   elif st=="sub x13, x13, x12, lsl #3":R[13]-=(R[12]<<3)
   elif st=="br x13":pc=R[13]//4
   elif st=="ldp x0, x1, [x15], #16":R[0]=src_m[R[15]];R[1]=src_m[R[15]+8];R[15]+=16
   elif st=="stp x0, x1, [x25], #16":dst_m[R[25]]=R[0];dst_m[R[25]+8]=R[1];R[25]+=16
   elif st=="cmp x15, x16":c_lt=(R[15]<R[16])
   elif st.startswith("b.lo "):
    if c_lt:pc=_cp_lbl[st.split()[1]]
   elif st=="ret": break
   else: raise ValueError(st)
  exp_d={dst_b+8*i:src_m[src_b+8*i] for i in range(n_c)}
  if dst_m!=exp_d or(dst_m and max(dst_m)+8>16384) or ops>n_c*(130.0/128.0)+11.0+1e-9:
   bfi_b2_errs+=1
  if n_c==2046 and(ops+59.0)/2046.0>1.05:bfi_b2_errs+=1
 _sp_tr=[(rng.randrange(8),rng.randrange(1,1<<41)) for _ in range(33)]
 _sp_te=[(js<<41)|ai for js,ai in _sp_tr]
 _m1,_o1=_ex_nb("x16_nobfi_pos",_sp_te[:19],1<<50,13)
 _m2,_o2=_ex_nb("x16_nobfi_pos",_sp_te[19:],1<<50,13)
 _mc=dict(_m1)
 for _ad,_v in _m2.items():_mc[_ad]=(_mc.get(_ad,0)+_v)&255
 _mr={}
 for js,ai in _sp_tr:_ad=(1<<50)+(js<<41)+ai;_mr[_ad]=(_mr.get(_ad,0)+13)&255
 if _mc!=_mr or(_o1+_o2)>33*(42.0/16.0)+2.0*11.55:bfi_b2_errs+=1
 _db1_ins=[];_db1_tins=[];_db1_cur=None;_db1_buf=None
 for _ln in open(kernels_path,encoding="ascii"):
  _s=_ln.strip()
  if _s=="//@seg delta_b1_apply":_db1_cur=0;continue
  if _s=="//@seg delta_b1_apply_tail":_db1_cur=1;continue
  if _s=="//@end":_db1_cur=None;continue
  if _db1_cur is None or not _s or _s.startswith("//"):continue
  if _s.startswith(".rept"):_db1_rep=int(_s.split()[1]);_db1_buf=[];continue
  if _s.startswith(".endr"):
   for _ in range(_db1_rep):(_db1_ins if _db1_cur==0 else _db1_tins).extend(_db1_buf)
   _db1_buf=None;continue
  (_db1_buf if _db1_buf is not None else(_db1_ins if _db1_cur==0 else _db1_tins)).append(_s)
 _db1_bod=[s for s in _db1_ins if not s.endswith(":")];_db1_tbod=[s for s in _db1_tins if not s.endswith(":")]
 def _ex_db1(ds,rbp,lbp,w8):
  nd=len(ds);nq=nd//4;lm={}
  for ii in range(0,nd,2):lm[lbp+4*ii]=(ds[ii]&0xFFFFFFFF)|(((ds[ii+1]&0xFFFFFFFF)<<32)if ii+1<nd else 0)
  for ii in range(nq*4,nd):lm[lbp+4*ii]=ds[ii]&0xFFFFFFFF
  R={0:0,1:0,4:lbp,5:rbp,6:w8&255,8:lbp+16*nq,9:lbp+4*nd};cm={};ops=4
  if nq>0:
   qr=nq%32;pc=(32-qr)*9 if qr else 0
   while pc<len(_db1_bod):
    st=_db1_bod[pc];ops+=1;pc+=1
    if st=="ldp x0, x1, [x4], #16":R[0]=lm[R[4]];R[1]=lm[R[4]+8];R[4]+=16
    elif st=="add x5, x5, w0, uxtw":R[5]=(R[5]+(R[0]&0xFFFFFFFF))&M64
    elif st=="add x5, x5, x0, lsr #32":R[5]=(R[5]+(R[0]>>32))&M64
    elif st=="add x5, x5, w1, uxtw":R[5]=(R[5]+(R[1]&0xFFFFFFFF))&M64
    elif st=="add x5, x5, x1, lsr #32":R[5]=(R[5]+(R[1]>>32))&M64
    elif st=="staddb w6, [x5]":cm[R[5]]=(cm.get(R[5],0)+R[6])&255
    elif st=="cmp x4, x8":clo=(R[4]<R[8])
    elif st=="b.lo delta_b1_apply":
     if clo:pc=0
    elif st=="ret":break
    else:raise ValueError(st)
  pc=0
  while pc<len(_db1_tbod):
   st=_db1_tbod[pc];ops+=1;pc+=1
   if st=="cmp x4, x9":clo=(R[4]<R[9])
   elif st=="b.hs .Ldb1_done":
    if not clo:break
   elif st=="ldr w0, [x4], #4":R[0]=lm[R[4]]&0xFFFFFFFF;R[4]+=4
   elif st=="add x5, x5, w0, uxtw":R[5]=(R[5]+(R[0]&0xFFFFFFFF))&M64
   elif st=="staddb w6, [x5]":cm[R[5]]=(cm.get(R[5],0)+R[6])&255
   elif st=="b.lo .Ldb1_tloop":
    if clo:pc=2
   elif st=="ret":break
   else:raise ValueError(st)
  return cm,ops
 win_e2e_errs=0;win_e2e_hits=0;C0=0;L0=14<<60
 if C0+int(2.0**_sb38["score_node_log2"])+1>=L0 or L0+int(2.0**_sb38["list_node_log2"])+1>=(1<<64):win_e2e_errs+=1
 for nd in(0,1,3,4,7,15,16,127,128,129,131,259):
  ds=[rng.randrange(1,1<<20)for _ in range(nd)];rbp=C0+(3<<38)+17;mg,og=_ex_db1(ds,rbp,L0+(1<<40),19);rf={};ax=rbp
  for dv in ds:ax=(ax+dv)&M64;rf[ax]=(rf.get(ax,0)+19)&255;win_e2e_hits+=1
  if mg!=rf or og>2.29*nd+24.0 or(mg and max(mg)>=L0):win_e2e_errs+=1
 b1a=sorted(rng.sample(range(1,1<<36),43))
 for rb in(0,7,19,31,41278121):
  rbp=C0+(rb<<38);mg,og=_ex_nb("x16_nobfi_pos",b1a,rbp,11);win_e2e_hits+=len(b1a)
  if mg!={rbp+ai:11 for ai in b1a}or og>(42.0/16.0)*len(b1a)+11.55 or max(mg)>=L0:win_e2e_errs+=1
 kw=41;sw=1<<kw;cr={};cg={}
 for si,ms,sn,kwa,cnd,rt in((0,"lo","x16_nobfi_bnd1_lo",dict(r25=2<<kw),lambda j:j>=2,74.0/16.0),(1,"int","x16_nobfi_pos",{},lambda j:True,42.0/16.0),(2,"hi","x16_nobfi_bnd1_hi",dict(r25=(6<<kw)-1),lambda j:j<=5,74.0/16.0),(3,"b2","x16_nobfi_bnd2",dict(r24=kw,r25=1,r26=4),lambda j:1<=j<=5,106.0/16.0)):
  sbp=C0+(1<<52)+si*(8*sw);sp=[(rng.randrange(8),rng.randrange(1,1<<36))for _ in range(25)]
  mg,og=_ex_nb(sn,[(j<<kw)|a for j,a in sp],sbp,15,**kwa)
  for ad,v in mg.items():cg[ad]=(cg.get(ad,0)+v)&255
  for j,a in sp:
   win_e2e_hits+=1
   if cnd(j):ad=sbp+(j<<kw)+a;cr[ad]=(cr.get(ad,0)+15)&255
  if og>rt*len(sp)+11.55:win_e2e_errs+=1
 if cg!=cr or(cg and max(cg)>=L0):win_e2e_errs+=1
 win_map_e2e_ok=(win_e2e_errs==0 and win_e2e_hits>=1000)
 delta_b1_errs=0
 ais_big=[10,250,(1<<32)+500,(1<<32)+900,(1<<33)+1200,(1<<33)+1205,(1<<33)+1299,(1<<33)+1400,(1<<33)+1555]
 segs_d=[]
 cur_base=ais_big[0];cur_ds=[0]
 prev=ais_big[0]
 for ai in ais_big[1:]:
  d=ai-prev
  if d>=(1<<32):
   segs_d.append((cur_base,cur_ds))
   cur_base=ai;cur_ds=[0];prev=ai
  else:
   cur_ds.append(d);prev=ai
 segs_d.append((cur_base,cur_ds))
 row0=1<<50;w_u8=17
 got_d={};ref_d={}
 for ai in ais_big:
  ref_d[row0+ai]=(ref_d.get(row0+ai,0)+w_u8)&0xFF
 for sb,ds in segs_d:
  x5=row0+sb
  for i in range(0,len(ds),2):
   w64=ds[i]|((ds[i+1]<<32) if i+1<len(ds) else 0)
   x5=(x5+(w64&0xFFFFFFFF))&M64
   got_d[x5]=(got_d.get(x5,0)+w_u8)&0xFF
   if i+1<len(ds):
    x5=(x5+(w64>>32))&M64
    got_d[x5]=(got_d.get(x5,0)+w_u8)&0xFF
 if got_d!=ref_d or len(segs_d)!=3:
  delta_b1_errs+=1
 T4=M.load_data()["t4_geom"]
 cA_emb=T4["GA"][0][0]/3.0
 gl_min_abs=min(abs(float(M.U_BN254+1)*complex(r_e,im_e)+2.0) for r_e,im_e in T4["EMB_W"])
 a_emb2_max=(572.0**2)*M.QA_MAX/cA_emb
 b_emb2_max=(782.76**2)*M.QB_MAX/cA_emb
 a_emb2_min=1.0/(a_emb2_max*a_emb2_max)
 ratio_max=b_emb2_max/(gl_min_abs*a_emb2_min)
 ng_diff_ub=max(6.0*log2(1.0+ratio_max),-6.0*log2(1.0-ratio_max))
 sl_u8=M.OP["slack"]
 sup_ng_cert=478.0
 thr_max=int(floor(S_U8*(sup_ng_cert-sl_u8)+0.5))
 init_min=128-thr_max
 max_w_ratio=max(int(floor(S_U8*(24.0+i*0.01)+0.5))/(24.0+i*0.01) for i in range(3701))
 max_u8=int(ceil(128.0+S_U8*sl_u8+0.5+(max_w_ratio-S_U8)*sup_ng_cert+0.02))
 r2d_chk=get_ret2d_u8(M,pmin=M.OP["p_min"],B=M.OP["lg"],odd=2)
 rg_marg_err=abs(r2d_chk["sv_rg_sum_350"]-r2d_chk["sv_rg_ref_350"])/max(r2d_chk["sv_rg_ref_350"],1e-30)
 sf_marg_err=abs(r2d_chk["sv_sf_sum_350"]-r2d_chk["sv_sf_ref_350"])/max(r2d_chk["sv_sf_ref_350"],1e-30)
 sf_468,s0_468=r2d_chk["get_sf_s0"](468,sl_u8)
 rg_468=r2d_chk["get_rg"](468,sl_u8)
 sf_468_1d=r2d_chk["ret_1d_sf"].survivor_upper(468.0,sl_u8,ex=0.0)
 bstep_c,dstep_c=r2d_chk["bstep"],r2d_chk["dstep"]
 ret2d_per_prime_errs=0
 for b_v in[24.51]+[24.0+i*0.01 for i in range(3701)]:
  wp_v=int(floor(S_U8*b_v+0.5))
  bi_v=round(b_v/bstep_c)
  db_v=bi_v*bstep_c-wp_v/S_U8
  if bi_v*bstep_c-ceil(db_v/dstep_c)*dstep_c>wp_v/S_U8+1e-12:
   ret2d_per_prime_errs+=1
  if bi_v*bstep_c-floor(db_v/dstep_c)*dstep_c<wp_v/S_U8-1e-12:
   ret2d_per_prime_errs+=1
 rho_335=M.rho(33.5);rho_400=M.rho(40.0)
 rho_tail_ok=(1e-60<rho_335<1e-55 and 1e-75<rho_400<1e-70 and len(M.rho_tab_1024())>=41*1024)
 r2d_mc=get_ret2d_u8(M,pmin=14,B=30.0,odd=2,max_bits=60)
 ps_mc=M._primes_upto(1<<14)
 _mc_small=ps_mc[:54];_mc_sp=1
 for _q in _mc_small:_mc_sp*=_q
 _mc_chunks=[]
 for _ci in range(54,len(ps_mc),64):
  _csub=ps_mc[_ci:_ci+64];_cp=1
  for _q in _csub:_cp*=_q
  _mc_chunks.append((_cp,_csub))
 _gcd_m=math.gcd;_isp_m=M._isp;_br_m=M._brent
 def _fac48(n_v):
  facs={};g=_gcd_m(n_v,_mc_sp)
  if g>1:
   for p_v in _mc_small:
    if g%p_v==0:
     c_v=0
     while n_v%p_v==0:n_v//=p_v;c_v+=1
     facs[p_v]=c_v;g//=p_v
     if g==1:break
  if n_v==1:return facs
  if _isp_m(n_v):facs[n_v]=1;return facs
  for _cp,_csub in _mc_chunks:
   g=_gcd_m(n_v,_cp)
   if g>1:
    for p_v in _csub:
     if g%p_v==0:
      c_v=0
      while n_v%p_v==0:n_v//=p_v;c_v+=1
      facs[p_v]=c_v;g//=p_v
      if g==1:break
    if n_v==1:return facs
    if _isp_m(n_v):facs[n_v]=1;return facs
  stk=[n_v]
  while stk:
   m_v=stk.pop()
   if m_v==1:continue
   if _isp_m(m_v):facs[m_v]=facs.get(m_v,0)+1
   else:d_v=_br_m(m_v);stk.append(d_v);stk.append(m_v//d_v)
  return facs
 mc_n=4000;mc_sl=22.0;mc_sm_tot=mc_sm_hit=mc_all_hit=0
 mc_lo=int(2.0**49.5);mc_hi=1<<50
 rng_mc=random.Random(seed^0x5A3C)
 for _ in range(mc_n):
  while True:
   v_mc=rng_mc.randrange(mc_lo,mc_hi)
   if v_mc%6 in(1,5): break
  facs=_fac48(v_mc)
  sc_mc=sum(int(floor(S_U8*log2(p_v)+0.5)) for p_v in facs if(1<<14)<=p_v<(1<<30))
  t_mc=int(floor(S_U8*(log2(v_mc)-mc_sl)+0.5))
  hit_mc=(sc_mc>=t_mc)
  if hit_mc:mc_all_hit+=1
  if all(p_v<(1<<30) for p_v in facs):
   mc_sm_tot+=1
   if hit_mc:mc_sm_hit+=1
 mc_emp_rg=mc_sm_hit/max(1,mc_sm_tot);mc_se_rg=sqrt(mc_emp_rg*(1.0-mc_emp_rg)/max(1,mc_sm_tot))
 mc_emp_sf=mc_all_hit/float(mc_n);mc_se_sf=sqrt(mc_emp_sf*(1.0-mc_emp_sf)/float(mc_n))
 mc_mod_rg=r2d_mc["get_rg"](50,mc_sl);mc_mod_sf,_=r2d_mc["get_sf_s0"](49,mc_sl)
 ret2d_mc_ok=(mc_sm_tot>=1400 and mc_emp_rg>=mc_mod_rg and mc_emp_sf<=mc_mod_sf)
 rot_tab_errs=0
 step1_inv_errs=0
 step1_tab_sim=[(((st*7+13)%729)<<3)|((st*3+5)&7) for st in range(729) for _ in range(8)]
 for st_t in range(729):
  ref_ent=step1_tab_sim[st_t<<3]
  for any3 in range(8):
   if step1_tab_sim[(st_t<<3)|any3]!=ref_ent:
    step1_inv_errs+=1
 seen_rot_idx={}
 for bmask_t in range(8):
  for l0_t in(-4194304,-4194303,-2100000,-1000001,-3,-2,-1,0,1,2,3,1234567,2100000,4194303):
   w6=l0_t&0xFFFFFFFF
   st_t=(bmask_t*89+abs(l0_t))%729
   w14=(st_t<<3)|bmask_t
   w6_after=(w6&~(0x7<<23))|((w14&0x7)<<23)
   sx=w6_after-(1<<32) if(w6_after&(1<<31)) else w6_after
   idx=(1<<26)+sx
   r_ref=((bmask_t|(bmask_t<<3))>>(l0_t%3))&7
   if not(0<=idx<(1<<27)):
    rot_tab_errs+=1
    continue
   sx_rec=idx-(1<<26)
   bm_rec=(sx_rec>>23)&7
   l0_lo23=sx_rec&0x7FFFFF
   l0_rec=(l0_lo23-(1<<23)) if(sx_rec<0) else l0_lo23
   r_tab=((bm_rec|(bm_rec<<3))>>(l0_rec%3))&7
   if bm_rec!=bmask_t or l0_rec!=l0_t or r_tab!=r_ref or idx in seen_rot_idx:
    rot_tab_errs+=1
   seen_rot_idx[idx]=r_tab
 sf_uses_survivor_upper_ok=(rg_marg_err<1e-3 and sf_marg_err<1e-3 and sf_468>=sf_468_1d and 0.90<rg_468<1.0 and s0_468<=sf_468 and ret2d_per_prime_errs==0 and ret2d_mc_ok and rho_tail_ok)
 u8_range_ok=(ng_diff_ub<=0.02 and thr_max<=128 and init_min>=0 and max_u8<=255
and sf_uses_survivor_upper_ok and rot_tab_errs==0 and step1_inv_errs==0)
 scan8_u8_errs=0
 n_c=128
 init_u8=[128-rng.randint(90,127) for _ in range(n_c)]
 rows_8=[init_u8[:] for _ in range(8)]
 exp_s=set()
 for rr in range(8):
  for cc in range(n_c):
   if rng.random()<0.03:
    rows_8[rr][cc]=(rows_8[rr][cc]+rng.randint(100,130))&0xFF
   if rows_8[rr][cc]&0x80:
    exp_s.add((rr,cc))
 fnd_s=set()
 for c0 in range(0,n_c,8):
  orr=0
  for rr in range(8):
   w=0
   for dc in range(8):
    w|=rows_8[rr][c0+dc]<<(8*dc)
   orr|=w
  if orr&0x8080808080808080:
   for rr in range(8):
    for dc in range(8):
     if rows_8[rr][c0+dc]&0x80:
      fnd_s.add((rr,c0+dc))
 if fnd_s!=exp_s:
  scan8_u8_errs+=1
 ins_rx=re.compile(r"^\s+[a-z][a-z0-9.]*\b")
 segs={};cur=None;rep_mult=1
 for raw in open(kernels_path,encoding="ascii"):
  st=raw.strip()
  if st.startswith("//@seg"):
   cur=st.split()[1];segs.setdefault(cur,0);rep_mult=1; continue
  if st.startswith("//@end"):
   cur=None;rep_mult=1; continue
  ln=raw.split("//",1)[0].rstrip()
  if cur is None or not ln.strip() or ln.endswith(":"):
   continue
  if ln.strip().startswith(".rept"):
   rep_mult=int(ln.strip().split()[1]); continue
  if ln.strip().startswith(".endr"):
   rep_mult=1; continue
  if ln.strip().startswith("."):
   continue
  if ins_rx.match(ln):
   segs[cur]+=rep_mult
 _b2_pg_ov=130.0/128.0+(59.0+11.0+22.0)/2046.0
 b1_delta_ops_ok=(segs.get("delta_b1_copy")==35 and segs.get("delta_b1_apply")==291 and segs.get("delta_b1_apply_tail")==8
and segs.get("x16_bfi_copy_page")==9 and segs.get("x16_bfi_copy")==131
and(130.0/128.0+(59.0+11.0)/2046.0)<=1.05 and M.C_SB_ENTER_X4>=42.60
and segs.get("x16_nobfi_pos")==42 and segs.get("x16_nobfi_bnd1_hi")==74
and segs.get("x16_nobfi_bnd1_lo")==74 and segs.get("x16_nobfi_bnd2")==106
and(42.0/16.0+_b2_pg_ov)<=M.C_B2_INT_X4
and(74.0/16.0+_b2_pg_ov)<=M.C_B2_BND1_HIT and(42.0/16.0+_b2_pg_ov)<=M.C_B2_BND1_MISS
and(106.0/16.0+_b2_pg_ov)<=M.C_B2_BND_HIT and(74.0/16.0+_b2_pg_ov)<=M.C_B2_BND_MISS
and(291.0/128.0)<=2.29 and(42.0/16.0)<=M.C_HIT_LO_WIN and M.C_SB_ENTER_WIN>=42.60
and(42.0/16.0)<=M.C_B2_INT_WIN and(74.0/16.0)<=M.C_B2_BND1_HIT_WIN and(42.0/16.0)<=M.C_B2_BND1_MISS_WIN
and(106.0/16.0)<=M.C_B2_BND_HIT_WIN and(74.0/16.0)<=M.C_B2_BND_MISS_WIN and win_map_e2e_ok
and segs.get("init8_u8")==75 and segs.get("scan8_u8")==199
and(74.0+197.0)/1024.0<=M.C_SCAN_G8 and 0.375<=M.C_SCAN_ZERO)
 pb_prod_df=M.pbits(0,62.9375,1.0)
 m_root_prod=ceil(pb_prod_df/64.0)
 root_bytes_prod=8*m_root_prod
 def _ntt3_len(m_v):
  kv=max(2,ceil(log2(m_v)))
  n2v=1<<kv
  n3v=3*(1<<(kv-2)) if kv>=3 else n2v
  return n3v if(kv>=3 and n3v>=m_v) else n2v
 def _kar_ntt3_scratch(m_v,m_base=(1<<55)):
  if m_v<=m_base:
   return 2*(8*_ntt3_len(m_v))+2*(8*m_v)
  hv=(m_v+1)//2
  return max(_kar_ntt3_scratch(hv,m_base),2*8*(hv+2)+_kar_ntt3_scratch(hv+2,m_base))
 prod_peak_bytes=0
 for d_lv in range(61):
  st_l=sum(ceil(m_root_prod/(1<<i)) for i in range(1,d_lv+1))
  m_dv=max(1,ceil(m_root_prod/(1<<d_lv)))
  sb_dv=8*st_l+8*(2*m_dv)+_kar_ntt3_scratch(m_dv)
  if sb_dv>prod_peak_bytes:
   prod_peak_bytes=sb_dv
 p1_n,p2_n,p3_n=M.NTT3_PRIMES
 ntt3_crt_no_ovf=(((1<<55)//2)*(M64*M64)<p1_n*p2_n*p3_n)
 inv12_n=pow(p1_n,-1,p2_n);inv123_n=pow((p1_n*p2_n)%p3_n,-1,p3_n)
 _ntt_roots={}
 for _pn in M.NTT3_PRIMES:
  _g=3
  while pow(_g,(_pn-1)//2,_pn)==1:_g+=1
  _ntt_roots[_pn]=_g
 def _ntt_r2(A_v,N_v,p_v,inv_f=False):
  wb=pow(_ntt_roots[p_v],(p_v-1)//N_v,p_v)
  if inv_f:wb=pow(wb,-1,p_v)
  a_v=[x%p_v for x in A_v]+[0]*(N_v-len(A_v))
  jv=0
  for iv in range(1,N_v):
   bit=N_v>>1
   while jv&bit:jv^=bit;bit>>=1
   jv^=bit
   if iv<jv:a_v[iv],a_v[jv]=a_v[jv],a_v[iv]
  ln=2
  while ln<=N_v:
   wl=pow(wb,N_v//ln,p_v);hf=ln>>1
   for iv in range(0,N_v,ln):
    wv=1
    for kv in range(hf):
     uv=a_v[iv+kv];vv=(a_v[iv+kv+hf]*wv)%p_v
     a_v[iv+kv]=(uv+vv)%p_v;a_v[iv+kv+hf]=(uv-vv)%p_v
     wv=(wv*wl)%p_v
   ln<<=1
  if inv_f:
   invN=pow(N_v,-1,p_v)
   a_v=[(x*invN)%p_v for x in a_v]
  return a_v
 _live_b=[0];_peak_b=[0]
 def _ab(nb):
  _live_b[0]+=nb
  if _live_b[0]>_peak_b[0]:_peak_b[0]=_live_b[0]
 def _fb(nb):_live_b[0]-=nb
 def _mul_ntt3_seq(A_v,B_v,m_o):
  Nv=1<<max(2,ceil(log2(max(1,m_o))))
  _ab(2*8*Nv)
  ta=_ntt_r2(A_v,Nv,p1_n);tb=_ntt_r2(B_v,Nv,p1_n)
  r0=_ntt_r2([(ta[i]*tb[i])%p1_n for i in range(Nv)],Nv,p1_n,True)[:m_o]
  _ab(8*m_o)
  ta=_ntt_r2(A_v,Nv,p2_n);tb=_ntt_r2(B_v,Nv,p2_n)
  r1=_ntt_r2([(ta[i]*tb[i])%p2_n for i in range(Nv)],Nv,p2_n,True)[:m_o]
  _ab(8*m_o)
  ta=_ntt_r2(A_v,Nv,p3_n);tb=_ntt_r2(B_v,Nv,p3_n)
  r2=_ntt_r2([(ta[i]*tb[i])%p3_n for i in range(Nv)],Nv,p3_n,True)[:m_o]
  out_v=[0]*m_o;cy=0
  for i in range(m_o):
   x1=r0[i];t2=((r1[i]-x1)*inv12_n)%p2_n;x12=x1+p1_n*t2
   t3=((r2[i]-x12)*inv123_n)%p3_n;val=x12+(p1_n*p2_n)*t3+cy
   out_v[i]=val&M64;cy=val>>64
  _fb(2*8*m_o+2*8*Nv)
  return out_v
 def _mul_kar_ntt3(A_v,B_v,m_o,m_base=8):
  if m_o<=m_base:
   return _mul_ntt3_seq(A_v,B_v,m_o)
  hv=(m_o+1)//2;hf=(hv+1)//2
  A0,A1=A_v[:hf],A_v[hf:];B0,B1=B_v[:hf],B_v[hf:]
  z0=_mul_kar_ntt3(A0,B0,2*hf,m_base)
  z2=_mul_kar_ntt3(A1,B1,max(1,m_o-2*hf),m_base)
  sal=max(len(A0),len(A1))+1;sbl=max(len(B0),len(B1))+1;z1l=sal+sbl
  _ab(8*(sal+sbl+z1l))
  SA=[0]*sal;cy=0
  for i in range(sal):
   s=(A0[i] if i<len(A0) else 0)+(A1[i] if i<len(A1) else 0)+cy;SA[i]=s&M64;cy=s>>64
  SB=[0]*sbl;cy=0
  for i in range(sbl):
   s=(B0[i] if i<len(B0) else 0)+(B1[i] if i<len(B1) else 0)+cy;SB[i]=s&M64;cy=s>>64
  Z1=_mul_kar_ntt3(SA,SB,z1l,m_base)
  for src in(z0,z2):
   bw=0
   for i in range(z1l):
    d=Z1[i]-(src[i] if i<len(src) else 0)-bw
    if d<0:Z1[i]=d+(1<<64);bw=1
    else:Z1[i]=d;bw=0
  out_v=[0]*m_o
  for i,v in enumerate(z0):out_v[i]=v
  for i,v in enumerate(z2):
   if 2*hf+i<m_o:out_v[2*hf+i]=v
  cy=0
  for i,v in enumerate(Z1):
   if hf+i<m_o:
    s=out_v[hf+i]+v+cy;out_v[hf+i]=s&M64;cy=s>>64
  kv=hf+z1l
  while cy and kv<m_o:
   s=out_v[kv]+cy;out_v[kv]=s&M64;cy=s>>64;kv+=1
  _fb(8*(sal+sbl+z1l))
  return out_v
 _kA=[rng.randrange(1<<64) for _ in range(32)]
 _kB=[rng.randrange(1<<64) for _ in range(32)]
 _ab(8*(32+32+64))
 _kP=_mul_kar_ntt3(_kA,_kB,64,m_base=8)
 _fb(8*(32+32+64))
 def _l2u(L_v):
  acc=0
  for i in range(len(L_v)-1,-1,-1):acc=(acc<<64)|L_v[i]
  return acc
 kar_ntt3_live_ok=(_l2u(_kP)==_l2u(_kA)*_l2u(_kB) and _live_b[0]==0 and _peak_b[0]<8*(8*64))
 df_leaves=M._primes_upto(8200)[:1024]
 assert len(df_leaves)==1024
 lvl=list(df_leaves)
 full_tree_bits=sum(x.bit_length() for x in lvl)
 while len(lvl)>1:
  lvl=[lvl[i]*lvl[i+1] for i in range(0,len(lvl),2)]
  full_tree_bits+=sum(x.bit_length() for x in lvl)
 ref_root=lvl[0]
 root_bits=ref_root.bit_length()
 live_stack_bits=0;max_stack_bits=0
 def _df_prod(arr):
  nonlocal live_stack_bits,max_stack_bits
  if len(arr)==1:
   return arr[0]
  mid=len(arr)//2
  pl=_df_prod(arr[:mid]);bl=pl.bit_length()
  live_stack_bits+=bl
  if live_stack_bits>max_stack_bits:max_stack_bits=live_stack_bits
  pr=_df_prod(arr[mid:])
  live_stack_bits-=bl
  return pl*pr
 df_root=_df_prod(df_leaves)
 boot_df_prod_peak_ratio=prod_peak_bytes/root_bytes_prod
 boot_df_tree_ok=(df_root==ref_root and live_stack_bits==0
and max_stack_bits<1.02*root_bits
and ntt3_crt_no_ovf and kar_ntt3_live_ok
and prod_peak_bytes<5.0*root_bytes_prod
and 8*root_bytes_prod<(1<<64)
and full_tree_bits>9.0*root_bits)
 def _u64(x): return x&M64
 def _mul64(a,b):
  p=(a&M64)*(b&M64)
  return p&M64,(p>>64)&M64
 def _int_to_u64(x,n):
  x&=(1<<(64*n))-1
  return[(x>>(64*i))&M64 for i in range(n)]
 def _u64_to_int(L):
  acc=0
  for i in range(len(L)-1,-1,-1):
   acc=(acc<<64)|(L[i]&M64)
  return acc
 def _srt_odd_norm8(A):
  k=0
  while k<8 and A[k]==0:k+=1
  v=A[k];s=(v&-v).bit_length()-1
  out=[0]*8
  for i in range(8):
   if i+k<8:
    lo=A[i+k]>>s
    hi=(A[i+k+1]<<(64-s))&M64 if(s>0 and i+k+1<8) else 0
    out[i]=lo|hi
  return out
 def _inv64_odd(a):
  a&=M64;x=1
  for _ in range(6):
   x=_u64(x*_u64(2-_u64(a*x)))
  return x
 _ops_cios8_one=(segs.get("cios8_setup",0)
+8*(segs.get("cios8_outer_pre",0)+8*segs.get("cios8_inner_mul",0)
+segs.get("cios8_outer_mid",0)+7*segs.get("cios8_inner_red",0)
+segs.get("cios8_outer_post",0))
+segs.get("cios8_csub_init",0)+8*segs.get("cios8_csub_loop",0)
+segs.get("cios8_csel_init",0)+8*segs.get("cios8_csel_loop",0)+20)
 _cios8_ins=[];_cios8_lbl={};_in_c8=False
 for _ln in Path(__file__).with_name("kernels_v32.s").read_text().splitlines():
  _s=_ln.strip()
  if _s=="//@seg cios8_setup":_in_c8=True; continue
  if _in_c8:
   if _s.startswith("//@seg ") or _s.startswith("//@end") or not _s or _s.startswith("//"):
    if _s=="//@end" and _cios8_lbl.get("cios8_csel_loop") is not None: break
    continue
   if _s.endswith(":"):_cios8_lbl[_s[:-1]]=len(_cios8_ins)
   else:
    _op,_rest=_s.split(None,1);_cios8_ins.append((_op,_rest))
 def _c8r(k): return 32 if k=="xzr" else int(k[1:])
 _c8_code=[]
 for _op,_rest in _cios8_ins:
  if _op=="stp":
   a,b,m=_rest.split(",",2);mp=[x.strip() for x in m.strip()[1:-1].split(",")]
   _c8_code.append((0,_c8r(a.strip()),_c8r(b.strip()),_c8r(mp[0]),int(mp[1][1:])//8 if len(mp)>1 else 0))
  elif _op in("ldr","str"):
   d,m=_rest.split(",",1);mp=[x.strip() for x in m.strip()[1:-1].split(",")]
   is_s=(_op=="str")
   if len(mp)==1:_c8_code.append((1 if is_s else 2,_c8r(d.strip()),_c8r(mp[0]),0,0))
   elif mp[1].startswith("#"):_c8_code.append((1 if is_s else 2,_c8r(d.strip()),_c8r(mp[0]),int(mp[1][1:])//8,0))
   else:_c8_code.append((3 if is_s else 4,_c8r(d.strip()),_c8r(mp[0]),_c8r(mp[1]),0))
  elif _op=="mov":d,im=[x.strip() for x in _rest.split(",")];_c8_code.append((5,_c8r(d),int(im[1:]),0,0))
  elif _op=="mul":d,a,b=[x.strip() for x in _rest.split(",")];_c8_code.append((6,_c8r(d),_c8r(a),_c8r(b),0))
  elif _op=="umulh":d,a,b=[x.strip() for x in _rest.split(",")];_c8_code.append((7,_c8r(d),_c8r(a),_c8r(b),0))
  elif _op=="adds":d,a,b=[x.strip() for x in _rest.split(",")];_c8_code.append((8,_c8r(d),_c8r(a),_c8r(b),0))
  elif _op=="adc":d,a,b=[x.strip() for x in _rest.split(",")];_c8_code.append((9,_c8r(d),_c8r(a),_c8r(b),0))
  elif _op=="add":d,a,im=[x.strip() for x in _rest.split(",")];_c8_code.append((10,_c8r(d),_c8r(a),int(im[1:]),0))
  elif _op=="sub":d,a,im=[x.strip() for x in _rest.split(",")];_c8_code.append((11,_c8r(d),_c8r(a),int(im[1:]),0))
  elif _op=="subs":d,a,b=[x.strip() for x in _rest.split(",")];_c8_code.append((12,_c8r(d),_c8r(a),_c8r(b),0))
  elif _op=="sbcs":d,a,b=[x.strip() for x in _rest.split(",")];_c8_code.append((13,_c8r(d),_c8r(a),_c8r(b),0))
  elif _op=="cmp":a,im=[x.strip() for x in _rest.split(",")];_c8_code.append((14,_c8r(a),int(im[1:]),0,0))
  elif _op=="b.lo":_c8_code.append((15,_cios8_lbl[_rest.strip()],0,0,0))
  elif _op=="cbnz":a,lb=[x.strip() for x in _rest.split(",")];_c8_code.append((16,_c8r(a),_cios8_lbl[lb],0,0))
  elif _op=="csel":d,a,b,cc=[x.strip() for x in _rest.split(",")];assert cc=="hs";_c8_code.append((17,_c8r(d),_c8r(a),_c8r(b),0))
  else: raise ValueError(_op)
 def _cios8_mul(A,B,N,n0inv):
  mem=[None]*18;mem[0]=list(A);mem[1]=list(B);mem[2]=list(N);mem[4]=[0]*10;mem[17]=[0]*8
  R=[0]*33;R[0]=0;R[1]=1;R[2]=2;R[3]=n0inv&M64;R[4]=4;R[7]=0xdeadbeef;R[8]=99;R[17]=17
  C=0;pc=0;ops=0;ni=len(_c8_code);m64=M64;cd=_c8_code
  while pc<ni:
   oc,a,b,c,d=cd[pc];pc+=1;ops+=1
   if oc==6:R[a]=(R[b]*R[c])&m64
   elif oc==7:R[a]=((R[b]*R[c])>>64)&m64
   elif oc==8:s=R[b]+R[c];C=1 if s>m64 else 0;R[a]=s&m64 if a<32 else 0
   elif oc==9:R[a]=(R[b]+R[c]+C)&m64 if a<32 else 0
   elif oc==4:m=mem[b];i=R[c];R[a]=m[i] if a<32 else 0
   elif oc==3:mem[b][R[c]]=R[a]
   elif oc==10:R[a]=(R[b]+c)&m64
   elif oc==14:C=1 if R[a]>=b else 0
   elif oc==15:
    if C==0:pc=a
   elif oc==2:m=mem[b];R[a]=m[c] if a<32 else 0
   elif oc==1:mem[b][c]=R[a]
   elif oc==0:m=mem[c];m[d]=R[a];m[d+1]=R[b]
   elif oc==5:R[a]=b
   elif oc==11:R[a]=(R[b]-c)&m64
   elif oc==12:va,vb=R[b],R[c];C=1 if va>=vb else 0;R[a]=(va-vb)&m64 if a<32 else 0
   elif oc==13:s=R[b]+(R[c]^m64)+C;C=1 if s>m64 else 0;R[a]=s&m64 if a<32 else 0
   elif oc==16:
    if R[a]!=0:pc=b
   elif oc==17:R[a]=R[b] if C else R[c]
  assert ops<=_ops_cios8_one
  return mem[17]
 def _sub8_limbs(A,B):
  D=[0]*8;bw=0
  for j in range(8):
   d=A[j]-B[j]-bw
   if d<0:D[j]=d+(1<<64);bw=1
   else:D[j]=d;bw=0
  return D,bw
 def _shr1_8_limbs(X):
  return[((X[j]>>1)|((X[j+1]&1)<<63)) for j in range(7)]+[X[7]>>1]
 def _bingcd8(U_in,V_in):
  U=U_in[:];V=V_in[:]
  if not any(x!=0 for x in V):
   return U,0
  iters=0;state=0
  while True:
   iters+=1
   if iters>1024:
    return None,iters
   if state==0:
    if(V[0]&1)==0:
     V=_shr1_8_limbs(V)
     if not any(x!=0 for x in V): return None,iters
    else:
     D,bw=_sub8_limbs(V,U)
     if not bw:
      V=_shr1_8_limbs(D)
      if not any(x!=0 for x in V): return U,iters
     else:
      D,bw2=_sub8_limbs(U,V)
      if bw2: return None,iters
      U=_shr1_8_limbs(D)
      if not any(x!=0 for x in U): return None,iters
      state=1
   else:
    if(U[0]&1)==0:
     U=_shr1_8_limbs(U)
     if not any(x!=0 for x in U): return None,iters
    else:
     D,bw=_sub8_limbs(U,V)
     if not bw:
      U=_shr1_8_limbs(D)
      if not any(x!=0 for x in U): return V,iters
     else:
      D,bw2=_sub8_limbs(V,U)
      if bw2: return None,iters
      V=_shr1_8_limbs(D)
      if not any(x!=0 for x in V): return None,iters
      state=0
 def _hensel_div8(N_odd_L,g_L):
  ginv=_inv64_odd(g_L[0]);rem=N_odd_L[:];Q=[0]*8
  for i in range(8):
   qi=_u64(rem[i]*ginv);Q[i]=qi
   cm=0;bs=0
   for j in range(8-i):
    lo,hi=_mul64(qi,g_L[j])
    p=lo+cm;cm=hi+(p>>64);p&=M64
    d=rem[i+j]-p-bs
    if d<0:rem[i+j]=d+(1<<64);bs=1
    else:rem[i+j]=d;bs=0
   if rem[i]!=0: return None
  return Q
 srt_9sq_errs=0;bingcd8_max_iters=0
 ps_small_s=M._primes_upto(64)
 P_small_s=1
 for p_v in ps_small_s:
  P_small_s*=p_v
 fa,fb=1,1
 while fa+fb<(1<<512):fa,fb=fb,fa+fb
 fib_odd,fib_prev=(fb,fa) if(fb&1) else(fa,fb)
 srt_cases=[
1<<480,
1<<511,
3**302,
(2**199)*(3**120)*(5**40)*1000000007,
(1<<479)-1,
(1<<512)-2,
fib_odd,
(fib_odd*2) if(fib_odd*2<(1<<512)) else(fib_odd-1),
]+[rng.randrange(1<<470,1<<512) for _ in range(56)]
 for xv in srt_cases:
  N_odd_L=_srt_odd_norm8(_int_to_u64(xv,8))
  N_odd=_u64_to_int(N_odd_L)
  if(N_odd&1)!=1:
   srt_9sq_errs+=1; continue
  r0=P_small_s%N_odd
  if N_odd==1:
   cof_got=1
  else:
   n0inv=_u64(-_inv64_odd(N_odd_L[0]))
   X_L=_int_to_u64(r0,8)
   for _ in range(9):
    X_L=_cios8_mul(X_L,X_L,N_odd_L,n0inv)
   g_L,it_v=_bingcd8(N_odd_L,X_L)
   if it_v>bingcd8_max_iters:bingcd8_max_iters=it_v
   if g_L is None or it_v>1024:
    srt_9sq_errs+=1; continue
   cof_L=_hensel_div8(N_odd_L,g_L)
   if cof_L is None:
    srt_9sq_errs+=1; continue
   cof_got=_u64_to_int(cof_L)
  cur=xv
  for p_v in ps_small_s:
   while cur%p_v==0:
    cur//=p_v
  if cof_got!=cur:
   srt_9sq_errs+=1
 for ua,va in(
(fib_odd,0),
(fib_odd,fib_odd),
(fib_odd,fib_prev),
((1<<512)-1,(1<<512)-3),
((1<<511)+1,1),
(((1<<255)-19)*((1<<255)-31),((1<<255)-19)*((1<<255)-101)),
):
  g_tL,it_t=_bingcd8(_int_to_u64(ua,8),_int_to_u64(va,8))
  if it_t>bingcd8_max_iters:bingcd8_max_iters=it_t
  if g_tL is None or it_t>1024 or _u64_to_int(g_tL)!=math.gcd(ua,va):
   srt_9sq_errs+=1
 _req_srt=(
"srt_odd_norm8_find","srt_odd_norm8_mid","srt_odd_norm8_shift","inv64_odd",
"cios8_setup","cios8_outer_pre","cios8_inner_mul","cios8_outer_mid",
"cios8_inner_red","cios8_outer_post","cios8_csub_init","cios8_csub_loop","cios8_csel_init","cios8_csel_loop",
"reg_bingcd8_io","reg_bingcd8_u_odd_chk","reg_bingcd8_u_odd_ge","reg_bingcd8_u_odd_lt","reg_bingcd8_v_even",
"reg_bingcd8_v_odd_chk","reg_bingcd8_v_odd_ge","reg_bingcd8_v_odd_lt","reg_bingcd8_u_even",
"reg_bingcd8_done_u","reg_bingcd8_done_v",
"hensel_div8_outer","hensel_div8_inner",
)
 _srt_segs_ok=all(segs.get(k,0)>0 for k in _req_srt) and(
segs.get("reg_bingcd8_u_odd_chk")==segs.get("reg_bingcd8_v_odd_chk")==1
and segs.get("reg_bingcd8_u_odd_ge")==segs.get("reg_bingcd8_v_odd_ge")==26
and segs.get("reg_bingcd8_u_odd_lt")==segs.get("reg_bingcd8_v_odd_lt")==17
and segs.get("reg_bingcd8_v_even")==segs.get("reg_bingcd8_u_even")==9
and segs.get("reg_bingcd8_done_u")==segs.get("reg_bingcd8_done_v")==4
)
 _ops_odd_norm=8*segs.get("srt_odd_norm8_find",0)+segs.get("srt_odd_norm8_mid",0)+8*segs.get("srt_odd_norm8_shift",0)+10
 _ops_cios8_9=segs.get("inv64_odd",0)+2+9*_ops_cios8_one
 for _Na,_Aa,_Ba in(((1<<512)-1,(1<<512)-2,(1<<512)-2),((1<<512)-3,(1<<512)-4,(1<<512)-4),((1<<511)+1,(1<<511),(1<<511))):
  _NL=_int_to_u64(_Na,8);_n0=_u64(-_inv64_odd(_NL[0]));_got=_u64_to_int(_cios8_mul(_int_to_u64(_Aa,8),_int_to_u64(_Ba,8),_NL,_n0))
  if _got!=(_Aa*_Ba*pow(1<<512,-1,_Na))%_Na:srt_9sq_errs+=1
 _bg_iter_wc=max(
segs.get("reg_bingcd8_u_odd_chk",0)+segs.get("reg_bingcd8_u_odd_ge",0),
segs.get("reg_bingcd8_u_odd_chk",0)+9+segs.get("reg_bingcd8_u_odd_lt",0),
segs.get("reg_bingcd8_u_odd_chk",0)+segs.get("reg_bingcd8_v_even",0),
)
 _ops_bingcd8=segs.get("reg_bingcd8_io",0)+1024*_bg_iter_wc+max(segs.get("reg_bingcd8_done_u",0),segs.get("reg_bingcd8_done_v",0))+1
 _ops_hensel=segs.get("inv64_odd",0)+8*segs.get("hensel_div8_outer",0)+36*segs.get("hensel_div8_inner",0)+20
 srt_leaf_9sq_wc_ops=_ops_odd_norm+_ops_cios8_9+_ops_bingcd8+_ops_hensel
 srt_g_leaf_ok=(_srt_segs_ok and srt_9sq_errs==0 and bingcd8_max_iters<=1024 and srt_leaf_9sq_wc_ops<=M.C_SRT_G_LEAF)
 def _limbs_to_int(L):
  acc=0
  for i in range(len(L)-1,-1,-1):
   acc=(acc<<64)|L[i]
  return acc
 def _sb_limbs(A,B):
  na,nb=len(A),len(B)
  C=[0]*(na+nb)
  for i in range(na):
   ai=A[i];carry=0
   for j in range(nb):
    p=ai*B[j]+C[i+j]+carry
    C[i+j]=p&M64
    carry=p>>64
   C[i+nb]=carry
  return C
 def _add_limbs(A,B):
  n=max(len(A),len(B))
  out=[0]*(n+1);c=0
  for i in range(n):
   s=(A[i] if i<len(A) else 0)+(B[i] if i<len(B) else 0)+c
   out[i]=s&M64;c=s>>64
  out[n]=c
  while len(out)>1 and out[-1]==0:out.pop()
  return out
 def _sub_limbs(A,B):
  out=[0]*len(A);b=0
  for i in range(len(A)):
   d=A[i]-(B[i] if i<len(B) else 0)-b
   if d<0:out[i]=d+(1<<64);b=1
   else:out[i]=d;b=0
  assert b==0
  while len(out)>1 and out[-1]==0:out.pop()
  return out
 def _cmp_limbs(A,B):
  la,lb=len(A),len(B)
  while la>0 and A[la-1]==0:la-=1
  while lb>0 and B[lb-1]==0:lb-=1
  if la!=lb: return 1 if la>lb else -1
  for i in range(la-1,-1,-1):
   if A[i]!=B[i]: return 1 if A[i]>B[i] else -1
  return 0
 def _kar_limbs(A,B):
  n=max(len(A),len(B))
  if n<=16:
   return _sb_limbs(A,B)
  m=n//2
  A0,A1=A[:m],A[m:]
  B0,B1=B[:m],B[m:]
  z0=_kar_limbs(A0,B0)
  z2=_kar_limbs(A1,B1)
  ca=_cmp_limbs(A1,A0);cb=_cmp_limbs(B1,B0)
  dA=_sub_limbs(A1,A0) if ca>=0 else _sub_limbs(A0,A1)
  dB=_sub_limbs(B1,B0) if cb>=0 else _sub_limbs(B0,B1)
  assert len(dA)<=max(len(A0),len(A1)) and len(dB)<=max(len(B0),len(B1))
  zd=_kar_limbs(dA,dB)
  z02=_add_limbs(z0,z2)
  z1=_sub_limbs(z02,zd) if ca*cb>=0 else _add_limbs(z02,zd)
  out=[0]*(len(A)+len(B))
  for off,arr in((0,z0),(m,z1),(2*m,z2)):
   c=0
   for i,v in enumerate(arr):
    if off+i<len(out):
     s=out[off+i]+v+c
     out[off+i]=s&M64;c=s>>64
   k=off+len(arr)
   while c and k<len(out):
    s=out[k]+c;out[k]=s&M64;c=s>>64;k+=1
  return out
 mul_kar_sb_errs=0
 for idx_k in range(200):
  n_l=(4,8,12,16,20,24,32,48,64)[idx_k%9]
  A_l=[M64]*n_l if idx_k<9 else [rng.randrange(1<<64) for _ in range(n_l)]
  B_l=[M64]*n_l if idx_k<9 else [rng.randrange(1<<64) for _ in range(n_l)]
  got_p=_limbs_to_int(_kar_limbs(A_l,B_l))
  ref_p=_limbs_to_int(A_l)*_limbs_to_int(B_l)
  if got_p!=ref_p:
   mul_kar_sb_errs+=1
 mul4_sb_ops=segs.get("mul4_sb",-1)
 mul_kar_sb_ok=(mul_kar_sb_errs==0 and mul4_sb_ops==105 and mul4_sb_ops<=int(M._kar_sb(4)))
 _v37=_verify_v37_ntt(M,segs,rng)
 v37_ntt_ok=_v37["ok"]
 _v39=_verify_v39_emit8(M,kernels_path)
 v39_emit8_ok=_v39["ok"]
 import descent_bn254_v6 as _D40
 _v40=_D40.verify_v40_descent(M,segs)
 v40_descent_ok=_v40["ok"]
 ok=(t3["ok"] and ntt3_split_errs==0 and ctail_errs==0 and ctail_ops<=M.OP["c_tail_z6"]
and bfi_b2_errs==0 and bfi_b2_cases==768 and win_map_e2e_ok and delta_b1_errs==0 and u8_range_ok and scan8_u8_errs==0
and b1_delta_ops_ok and boot_df_tree_ok and srt_g_leaf_ok and mul_kar_sb_ok and v37_ntt_ok and v39_emit8_ok and v40_descent_ok)
 return dict(
ok=ok,
win_map_e2e_ok=win_map_e2e_ok,
win_map_e2e_hits=win_e2e_hits,
win_map_e2e_errs=win_e2e_errs,
adic3_exact=t3,
ntt3_split_proto_errs=ntt3_split_errs,
ctail_lazy_errs=ctail_errs,
ctail_lazy_ops=ctail_ops,
ctail_budget=M.OP["c_tail_z6"],
bfi_b2_errs=bfi_b2_errs,
bfi_b2_cases=bfi_b2_cases,
b2_prefix_mult=round(_b2_pref_mult,4),
b2_active_ks=_act_ks,
delta_b1_errs=delta_b1_errs,
ng_diff_ub_bits=round(ng_diff_ub,6),
u8_scale=S_U8,u8_thr_max=thr_max,u8_init_min=init_min,u8_max_val=max_u8,
ret2d_m_margin_bits=round(r2d_chk["m_margin"],6),
ret2d_rg_marg_err=round(rg_marg_err,8),
ret2d_sf_marg_err=round(sf_marg_err,8),
ret2d_per_prime_errs=ret2d_per_prime_errs,
ret2d_mc_ok=ret2d_mc_ok,
ret2d_mc_emp_rg=round(mc_emp_rg,4),
ret2d_mc_mod_rg=round(mc_mod_rg,4),
ret2d_mc_emp_sf=round(mc_emp_sf,4),
ret2d_mc_mod_sf=round(mc_mod_sf,4),
rho_tail_ok=rho_tail_ok,
rho_335=f"{rho_335:.6e}",
rho_400=f"{rho_400:.6e}",
rot_tab_2s_errs=rot_tab_errs,
step1_tab_8192_inv_errs=step1_inv_errs,
sf_uses_survivor_upper_ok=sf_uses_survivor_upper_ok,
u16_range_ok=u8_range_ok,
scan8_u16_errs=scan8_u8_errs,
b1_delta_ops_ok=b1_delta_ops_ok,
boot_df_tree_ok=boot_df_tree_ok,
boot_df_stack_ratio=round(max_stack_bits/root_bits,4),
boot_df_peak_ratio=round(boot_df_prod_peak_ratio,4),
boot_df_prod_peak_log2=round(log2(prod_peak_bytes),4),
boot_kar_ntt3_live_peak_ratio=round(_peak_b[0]/(8*64),4),
ntt3_crt_no_ovf=ntt3_crt_no_ovf,
boot_full_tree_ratio=round(full_tree_bits/root_bits,2),
srt_g_leaf_ok=srt_g_leaf_ok,
srt_g_leaf_9sq_wc_ops=int(srt_leaf_9sq_wc_ops),
srt_g_leaf_budget=M.C_SRT_G_LEAF,
bingcd8_max_iters=bingcd8_max_iters,
bingcd8_iter_wc_ops=_bg_iter_wc,
mul_kar_sb_ok=mul_kar_sb_ok,
mul_kar_sb_errs=mul_kar_sb_errs,
mul4_sb_ops=mul4_sb_ops,
v37_ntt_ok=v37_ntt_ok,
v37_r4_inner=_v37["r4_in"],v37_r4_hoist_x2=_v37["r4_ho"],v37_r3=_v37["r3"],v37_r5=_v37["r5"],
v37_r4_blend=_v37["c_bl"],v37_r3_lin=_v37["c_r3"],v37_r5_lin=_v37["c_r5"],
v37_bf_errs=_v37["bf_errs"],v37_conv_errs=_v37["conv_errs"],
v39_emit8_ok=v39_emit8_ok,
v39_emit8_flush_amort=_v39["flush_amort"],
v39_emit8_b1cpu_errs=_v39["cpu_errs"],
v39_max_rows_per_Cs=_v39["max_rows"],
v39_kmax_525=_v39["k525"],v39_kmax_53=_v39["k53"],v39_kmax_555=_v39["k555"],
v39_sink_ub=_v39["sink_ub"],v39_sink_reset_ok=_v39["sink_reset_ok"],
v40_descent_ok=v40_descent_ok,
v43_cios5_ir=_v40["cios5"],v43_mont2_ir=_v40["mont2"],v43_bingcd_gen_ir=_v40["bingcd_gen"],
v43_psurv_iv_ok=_v40["psurv_iv_ok"],v43_kron_ntt3_ok=_v40["kron_ntt3_ok"],
v43_bls_mpe_ok=_v40["bls_mpe_ok"],v43_ecm_ram_mib=_v40["ecm_ram_mib"],
)
def _verify_v37_ntt(M,segs,rng):
 M64=(1<<64)-1
 r4_in,r4_ho=segs.get("r4_4clamp_inner_tw",-1),segs.get("r4_4clamp_hoist_tw_x2",-1)
 r3_bf,r5_bf=segs.get("r3_shoup_bfly",-1),segs.get("r5_shoup_bfly",-1)
 c_bl=0.5*(9.0*r4_in/8.0+9.0*r4_ho/16.0)
 c_r3,c_r5=9.0*r3_bf/3.0+220.0,9.0*r5_bf/5.0+220.0
 def is_p(n):
  if n<2: return False
  for q in(2,3,5,7,11,13,17,19,23,29,31,37):
   if n==q: return True
   if n%q==0: return False
  d,s=n-1,0
  while(d&1)==0:d>>=1;s+=1
  for a in(2,3,5,7,11,13,17,19,23,29,31,37):
   x=pow(a,d,n)
   if x in(1,n-1): continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1: break
   else: return False
  return True
 def proot(p):
  phi,fs,t=p-1,set(),p-1
  for q in(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61):
   while t%q==0:fs.add(q);t//=q
  d=67
  while d*d<=t:
   while t%d==0:fs.add(d);t//=d
   d+=2
  if t>1:fs.add(t)
  g=2
  while any(pow(g,phi//q,p)==1 for q in fs):g+=1
  return g
 def sm(v,w,wp,p): return(v*w-((v*wp)>>64)*p)&M64
 def c2(x,p2): return x-p2 if x>=p2 else x
 def bf4(a,ws,J,Jp,p):
  p2=2*p
  t0=c2(a[0],p2);t1=sm(a[1],ws[0],(ws[0]<<64)//p,p)
  t2=sm(a[2],ws[1],(ws[1]<<64)//p,p);t3=sm(a[3],ws[2],(ws[2]<<64)//p,p)
  u0,u1=c2(t0+t2,p2),c2(t0+p2-t2,p2)
  v0,v1=c2(t1+t3,p2),sm(t1+p2-t3,J,Jp,p)
  return(u0+v0,u1+v1,u0+p2-v0,u1+p2-v1)
 def bf3(a,ws,c1,c1p,c2v,c2p,p):
  p2=2*p
  t0=c2(a[0],p2);t1=sm(a[1],ws[0],(ws[0]<<64)//p,p);t2=sm(a[2],ws[1],(ws[1]<<64)//p,p)
  u=c2(t1+t2,p2);m1=sm(u,c1,c1p,p);m2=sm(t1+p2-t2,c2v,c2p,p);r=c2(t0+m1,p2)
  return(t0+u,r+m2,r+p2-m2)
 def r3_c(w3,p):
  i2=(p+1)//2;c1=(p-i2)%p;c2v=((w3-(w3*w3)%p)*i2)%p
  return c1,(c1<<64)//p,c2v,(c2v<<64)//p
 def r5_c(w5,p):
  w2=(w5*w5)%p;w3=(w2*w5)%p;w4=(w3*w5)%p;i2=(p+1)//2;i4=(i2*i2)%p
  c0=(-i4)%p;c1=((w5+w4-w2-w3)*i4)%p;d1=((w5-w4)*i2)%p;d2=((w2-w3)*i2)%p
  k0,k1,k2=d2,(d1-d2)%p,(-(d1+d2))%p
  return tuple((x,(x<<64)//p) for x in(c0,c1,k0,k1,k2))
 def bf5(a,ws,cs,p):
  p2=2*p;(c0,c0p),(c1,c1p),(k0,k0p),(k1,k1p),(k2,k2p)=cs
  t0=c2(a[0],p2);ts=[sm(a[i+1],ws[i],(ws[i]<<64)//p,p) for i in range(4)]
  s1,s2=c2(ts[0]+ts[3],p2),c2(ts[1]+ts[2],p2)
  d1,d2=c2(ts[0]+p2-ts[3],p2),c2(ts[1]+p2-ts[2],p2)
  s12=c2(s1+s2,p2);u1=c2(t0+sm(s12,c0,c0p,p),p2);m1=sm(s1+p2-s2,c1,c1p,p)
  r1,r2=c2(u1+m1,p2),c2(u1+p2-m1,p2)
  q0=sm(d1+d2,k0,k0p,p);z1,z2=c2(q0+sm(d1,k1,k1p,p),p2),c2(q0+sm(d2,k2,k2p,p),p2)
  return(t0+s12,r1+z1,r2+z2,r2+p2-z2,r1+p2-z1)
 p1,p2=M.NTT2_52_PRIMES
 primes_ok=(all(is_p(p) and(p-1)%(3*(1<<55))==0 and 4*p<(1<<64) for p in M.NTT3_PRIMES)
and all(is_p(p) and(p-1)%(15*(1<<40))==0 and 4*p<(1<<64) for p in M.NTT2_52_PRIMES)
and(p1*p2)>(1<<18)*((1<<52)-1)**2*3)
 proots={p:proot(p) for p in M.NTT3_PRIMES+M.NTT2_52_PRIMES}
 bf_errs=0
 for p in M.NTT3_PRIMES+M.NTT2_52_PRIMES:
  p4,g=4*p,proots[p];J=pow(g,(p-1)//4,p);Jp=(J<<64)//p
  for idx in range(516):
   xs=[p4-1 if(idx>>b)&1 else 0 for b in range(4)] if idx<16 else[rng.randrange(p4) for _ in range(4)]
   ws=[rng.randrange(p) for _ in range(3)];ys=bf4(xs,ws,J,Jp,p)
   aa=[xs[0]%p]+[(xs[i+1]*ws[i])%p for i in range(3)]
   if max(ys)>=p4 or[y%p for y in ys]!=[(aa[0]+aa[2]+aa[1]+aa[3])%p,(aa[0]-aa[2]+J*(aa[1]-aa[3]))%p,(aa[0]+aa[2]-aa[1]-aa[3])%p,(aa[0]-aa[2]-J*(aa[1]-aa[3]))%p]:
    bf_errs+=1
  w3=pow(g,(p-1)//3,p);c3=r3_c(w3,p)
  for idx in range(258):
   xs=[p4-1 if(idx>>b)&1 else 0 for b in range(3)] if idx<8 else[rng.randrange(p4) for _ in range(3)]
   ws=[rng.randrange(p) for _ in range(2)];ys=bf3(xs,ws,*c3,p)
   aa=[xs[0]%p]+[(xs[i+1]*ws[i])%p for i in range(2)]
   if max(ys)>=p4 or[y%p for y in ys]!=[(aa[0]+aa[1]+aa[2])%p,(aa[0]+w3*aa[1]+w3*w3*aa[2])%p,(aa[0]+w3*w3*aa[1]+w3*aa[2])%p]:
    bf_errs+=1
  if(p-1)%5==0:
   w5=pow(g,(p-1)//5,p);c5=r5_c(w5,p)
   for idx in range(282):
    xs=[p4-1 if(idx>>b)&1 else 0 for b in range(5)] if idx<32 else[rng.randrange(p4) for _ in range(5)]
    ws=[rng.randrange(p) for _ in range(4)];ys=bf5(xs,ws,c5,p)
    aa=[xs[0]%p]+[(xs[i+1]*ws[i])%p for i in range(4)]
    if max(ys)>=p4 or[y%p for y in ys]!=[sum(aa[j]*pow(w5,i*j,p) for j in range(5))%p for i in range(5)]:
     bf_errs+=1
 def mr_ntt(a,p,inv=False):
  n=len(a);wn=pow(proots[p],(p-1)//n,p)
  if inv:wn=pow(wn,p-2,p)
  def rec(v,ws):
   m=len(v)
   if m==1: return[v[0]]
   r=5 if m%5==0 else(3 if m%3==0 else(4 if m%4==0 else 2))
   sb=m//r;pts=[rec(v[i::r],pow(ws,r,p)) for i in range(r)];out=[0]*m;wk=1
   w_r=pow(ws,sb,p)
   c3=r3_c(w_r,p) if r==3 else None;c5=r5_c(w_r,p) if r==5 else None;Jp=(w_r<<64)//p if r==4 else 0
   for j in range(sb):
    pw=[wk]
    for _ in range(r-2):pw.append((pw[-1]*wk)%p)
    inp=[pts[i][j] for i in range(r)]
    if r==5:ys=bf5(inp,pw,c5,p)
    elif r==3:ys=bf3(inp,pw,*c3,p)
    elif r==4:ys=bf4(inp,pw,w_r,Jp,p)
    else:
     t0,t1=c2(inp[0],2*p),sm(inp[1],wk,(wk<<64)//p,p)
     ys=(t0+t1,t0+2*p-t1)
    for i in range(r):out[j+i*sb]=ys[i]
    wk=(wk*ws)%p
   return out
  res=rec([x%p for x in a],wn)
  if inv:
   iv=pow(n,p-2,p)
   return[(x*iv)%p for x in res]
  return[x%p for x in res]
 conv_errs=0;inv12=pow(p1,p2-2,p2);B52=1<<52
 for nv in(12,20,24,30,48,60):
  hf=nv//2
  for md in(0,1):
   A=([B52-1]*hf if md==0 else[rng.randrange(B52) for _ in range(hf)])+[0]*(nv-hf)
   B=([B52-1]*hf if md==0 else[rng.randrange(B52) for _ in range(hf)])+[0]*(nv-hf)
   c1=mr_ntt([(x*y)%p1 for x,y in zip(mr_ntt(A,p1),mr_ntt(B,p1))],p1,True)
   c2v=mr_ntt([(x*y)%p2 for x,y in zip(mr_ntt(A,p2),mr_ntt(B,p2))],p2,True)
   pi=cy=0
   for i in range(nv):
    s=c1[i]+p1*(((c2v[i]-c1[i])*inv12)%p2)+cy
    pi|=(s&(B52-1))<<(52*i);cy=s>>52
   pi|=cy<<(52*nv)
   if pi!=sum(A[i]<<(52*i) for i in range(hf))*sum(B[i]<<(52*i) for i in range(hf)):
    conv_errs+=1
 for p in M.NTT3_PRIMES:
  for nv in(12,24):
   A,B=[rng.randrange(p) for _ in range(nv)],[rng.randrange(p) for _ in range(nv)]
   cc=mr_ntt([(x*y)%p for x,y in zip(mr_ntt(A,p),mr_ntt(B,p))],p,True)
   if cc!=[sum(A[j]*B[(i-j)%nv] for j in range(nv))%p for i in range(nv)]:
    conv_errs+=1
 ok=(r4_in==45 and r4_ho==74 and r3_bf==36 and r5_bf==93
and c_bl<=46.5 and c_r3<=364.0 and c_r5<=420.0
and primes_ok and bf_errs==0 and conv_errs==0)
 return dict(ok=ok,r4_in=r4_in,r4_ho=r4_ho,r3=r3_bf,r5=r5_bf,
c_bl=round(c_bl,4),c_r3=round(c_r3,2),c_r5=round(c_r5,2),
bf_errs=bf_errs,conv_errs=conv_errs)
def _verify_v39_emit8(M,kernels_path):
 import stor_bound as S,v20_lat_b2 as V20,v25_verify as V25
 from fractions import Fraction as F
 from math import isqrt
 gchk=M.verify_exhaustive_fd(R2=50.0)
 eps_cap=F(3,100000)
 ga00_lo=F(S.G[0][0],S.G_SC)-eps_cap*F(88,100)
 c_amgm=(F(1)-eps_cap)*ga00_lo*S.G_SC
 r2_max=(F(1)+S.DELTA)*S.G_SC*(S.R_A_MAX**2)
 kmap={}
 for num,den in((105,2),(53,1),(107,2),(54,1),(109,2),(111,2),(123,2)):
  pw=3*den;sh=num+pw*64;lo,hi=0,1<<(sh//pw+2)
  while hi-lo>1:
   mid=(lo+hi)//2
   if mid**pw<=(1<<sh):lo=mid
   else:hi=mid
  d2=c_amgm*F(lo,1<<64)
  if d2<=2*r2_max: return dict(ok=False,flush_amort=99.0,cpu_errs=1,max_rows=0,k525=0,k53=0,k555=0,sink_ub=0,sink_reset_ok=False)
  km=int(d2//(d2-2*r2_max))
  if d2>4*r2_max and km!=1: return dict(ok=False,flush_amort=99.0,cpu_errs=1,max_rows=0,k525=0,k53=0,k555=0,sink_ub=0,sink_reset_ok=False)
  kmap[round(num/den,4)]=km
 sb=S.main();max_rows=sb["max_rows_per_Cs"];ep_recs=sb["epoch_buf_bytes"]//8
 bands=((52.5,53.0,13),(53.0,53.5,5),(53.5,54.0,3),(54.0,55.5,2),(55.5,61.50,1))
 tot_ep=0.0;bands_ok=(V20.EMIT8_P_MIN==6369051672525773 and V20.EMIT8_P_MIN*V20.EMIT8_P_MIN>(1<<105))
 for blo,bhi,km in bands:
  b_ep=ep_recs//(6*km*max_rows)
  if kmap[blo]>km or b_ep<1 or b_ep*6*km*max_rows>ep_recs:bands_ok=False
  tot_ep+=ceil((2.0*M.KH_G_CAP)*(2.0**M.log2_li(bhi)-2.0**M.log2_li(blo))/b_ep)
 cells_lb=2.0**94.889
 flush_amort=(tot_ep*float(S.N_SC)*50.0)/(cells_lb*1.07*log(61.50/52.5))
 sgs8={};b2_lines=[];cur=None;after_e8=[];saw_e8_end=False
 for raw in open(kernels_path,encoding="ascii"):
  st=raw.strip()
  if st.startswith("//@seg"):cur=st.split()[1];saw_e8_end=False; continue
  if st.startswith("//@end"):
   if cur=="emit8_last":saw_e8_end=True
   cur=None; continue
  ln=raw.split("//",1)[0].strip()
  if not ln or ln.endswith(":") or ln.startswith("."): continue
  if saw_e8_end and len(after_e8)<2:after_e8.append(ln)
  if cur in("wheel_disp4","wheel_skip2_stub","emit8_step1","emit8_step2","emit8_step3","emit8_last"):sgs8.setdefault(cur,[]).append(ln)
  if cur=="body_2":b2_lines.append(ln)
 t12_eff=(F(1)+S.DELTA)*S.W_A*S.G_SC*(S.R_A_MAX**2)
 b0_min_sq=F(1000*753729);b1_min_sq=F(74,100)*b0_min_sq
 n0_ub=2*(isqrt(int(t12_eff//b0_min_sq))+1)+3
 n1_ub=2*(isqrt(int(t12_eff//b1_min_sq))+1)+3
 sink_reset_ok=(len(b2_lines)>=2 and b2_lines[-2]=="ldr x8, [sp, #416]" and b2_lines[-1]=="str x8, [x17]" and n0_ub*n1_ub<=ep_recs)
 e8_tail_ok=(after_e8==["b.pl .Lk0_top","b .Ltop_2"])
 wd8=sgs8.get("wheel_disp4",[]);sk2=sgs8.get("wheel_skip2_stub",[]);em_lines=sgs8.get("emit8_step1",[])
 c8=wd8+sk2+em_lines+sgs8.get("emit8_step2",[])+sgs8.get("emit8_step3",[])+sgs8.get("emit8_last",[])+after_e8
 Pe=V25._b1_decode(c8,{".Lk0_top":len(c8),".Ltop_2":len(c8)+1})
 m8={};AT=1<<48;BT=1<<50;PT=1<<52;RT=1<<45;RR=RT+(1<<20);sm=(1<<34)-1
 sl=(0,1,sm);pb={s_:(1<<54)+(1<<27)*(i+1) for i,s_ in enumerate(sl)}
 for s_ in sl:m8[PT+8*s_]=pb[s_]
 kk=(((-5)&0xFFFFFFFF)&~(7<<23))|(1<<23);ks=kk-(1<<32) if(kk>>31) else kk
 m8[(RT+(ks<<3))&((1<<64)-1)]=RR;m8[RR+6*8]=4*(len(wd8)+len(sk2));m8[RR+8*8]=4*len(wd8)
 cpu=V25._B1CPU(Pe,m8,{});x=cpu.x
 x[1],x[10],x[11],x[12],x[26],x[27],x[17]=8,RT,48,3,16,1,PT
 x[4],x[13]=0xDEADBEEFCAFE1234,0x123456789ABCDEF0
 ais=(1,(1<<37)-1,1<<38,(1<<47)+12345,(1<<48)-1);exp={s_:[] for s_ in sl}
 offs=(0,1,3,6)
 for ch in range(128):
  for st_i,off_i in enumerate(offs):
   i=ch*4+st_i;pos=ch*7+off_i
   s_i=sl[i%3];rb_i=0 if s_i==0 else(43963475 if s_i==sm else 1+(i%7));ai_i=ais[i%5]
   bk=BT+(i&3)*(1<<40);m8[AT+16*pos]=ai_i;m8[AT+16*pos+8]=bk;m8[bk+8*pos]=(s_i<<30)|rb_i
   exp[s_i].append((rb_i<<38)+ai_i)
 tally={};cpu_errs=0;x[23]=127
 for ch in range(128):
  x[2],x[5],x[6],x[14],x[0]=AT+16*(ch*7-(2 if ch==0 else 0)),ch*7-(2 if ch==0 else 0),(-5)&0xFFFFFFFF,1,(8 if ch==0 else 6)
  x[23]=(x[23]-1)&((1<<64)-1);cpu.N=bool(x[23]>>63)
  cpu.run(0,{len(c8),len(c8)+1},{},tally)
 if sum(tally.values())!=128*35+2+1 or not e8_tail_ok or len(sk2)!=2:cpu_errs+=1
 if x[4]!=0xDEADBEEFCAFE1234 or x[13]!=0x123456789ABCDEF0 or x[1]!=8:cpu_errs+=1
 for s_ in sl:
  if(m8[PT+8*s_]-pb[s_])//8!=len(exp[s_]):cpu_errs+=1
  for j,ev in enumerate(exp[s_]):
   if m8.get(pb[s_]+8*j)!=ev:cpu_errs+=1
 fold_ok=all((n1+1)**2<=937*n1+1 for n1 in range(n1_ub+1)) and(((n1_ub+2)*n1_ub+1)<<32)<(1<<64) and n1_ub==909
 ok=(gchk["ok"] and F(str(gchk["eps_G_hi"]))<eps_cap and sb["ok"] and bands_ok
and flush_amort<=V20.EMIT8_FLUSH_AMORT and len(em_lines)==8 and cpu_errs==0 and fold_ok and sink_reset_ok)
 return dict(ok=ok,flush_amort=round(flush_amort,5),cpu_errs=cpu_errs,max_rows=max_rows,k525=kmap[52.5],k53=kmap[53.0],k555=kmap[55.5],sink_ub=n0_ub*n1_ub,sink_reset_ok=sink_reset_ok)
