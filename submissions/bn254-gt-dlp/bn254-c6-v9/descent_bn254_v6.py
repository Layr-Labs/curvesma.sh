"""descent_bn254_v6.py — Pure-stdlib live 62-row LB=59 pure k=1 zero-retry Boole + Wald i.i.d. fresh-e boot + 6-bit ECM ladder evaluator."""
import functools,math
from math import exp,factorial,floor,ceil,log,log1p,log2,sqrt
from fractions import Fraction as F
LN2=log(2.0)
UNIT=20000.0
RSTEP=1.0/1024.0
KA_MAX=16
MONT5=367.76
TORSION=12
WHEELS=(30030,)
SEPARATION=0.5
RET_BOOT=0.9
SURV_X=2.0
SURV_DM=10.0
NB_R=1.0
NB_ALPHA=0.5
P0_FLOOR=0.001
DESC_CAP_LOG2=112.0
LB_WAVE=44
LOG2P6=6.0*log2(21888242871839275222246405745257275088696311157297823662689037894645226208583)
LOG2U6=6.0*log2(4965661367192848881)
QBINS=[(1/6,-1.23,-1.00),(1/6,8.10,2.72),(1/6,14.96,5.42),(1/6,20.79,7.70),(1/6,27.67,10.30),(1/6,55.33,18.25)]
QBINS_LO=[(-95.0,-55.0),(-10.0,-55.0),(0.0,-55.0),(8.0,-55.0),(14.0,-55.0),(20.0,-55.0)]
BOOT_BINS=[(0.2,-184.45),(0.2,-181.82),(0.2,-180.07),(0.2,-178.14),(0.2,-170.68)]
BOOT_LO_OFF=-208.0
DF=48
ALPHA_F=0.0
ALPHA_G=3.4
SLACK=96.0
UDUP=0.6869
MESH=96
ECM_SUCCESS=0.999
BOUND_CELL=2400.0
HITS_PER_CELL=2.0*log(59.0/8.0)
SIEVE_HIT=1100.0
SCAN=10.0
ZERO=2.0
PAIR_SETUP=2.0**26.2
LLL_NODE=2.0**33
NORMS=2.0**16
TDIV=2.0**17
IDEAL=2.0**26
PRIME=2.0**24
GRID=4
SL_RANGE=(328,380)
def _log2li(l):
 x=l*LN2
 return l-log2(x)+log2(1.0+1.0/x+2.0/(x*x))
def _mulmod_instr(nbits):
 L=max(5,ceil(nbits/64.0))
 return MONT5*(L/5.0)**2,L
def _gcd_instr(nbits):
 L=max(5,ceil(nbits/64.0))
 return 2.0*nbits*(12.0*L+16.0)
def _pi_upper(x):
 return 1.26*x/log(x)if x>=17 else 7.0
def _stage_ops(b1,b2,wheel):
 s1=ceil(1.5*b1)*11
 phi=wheel
 for q in(2,3,5,7,11,13):
  if wheel%q==0:
   phi=phi//q*(q-1)
 J=phi//2
 G=ceil(b2/wheel)+1
 s2=6*(wheel//2+G)+3*J*G if b2>b1 else 0
 return 1.1*(s1+s2),s1,s2,J,G
def _conv(a,b):
 na,nb=len(a),len(b)
 out=[0.0]*(na+nb-1)
 i0=0
 while i0<na and a[i0]==0.0:i0+=1
 if i0==na:return out
 i1=na-1
 while i1>i0 and a[i1]==0.0:i1-=1
 j0=0
 while j0<nb and b[j0]==0.0:j0+=1
 if j0==nb:return out
 j1=nb-1
 while j1>j0 and b[j1]==0.0:j1-=1
 aa=a[i0:i1+1];bb=b[j0:j1+1]
 la,lb=len(aa),len(bb)
 if la>lb:
  aa,bb=bb,aa
  la,lb=lb,la
 aa_rev=aa[::-1]
 _sp=math.sumprod
 base=i0+j0
 for s in range(la+lb-1):
  k_lo=0 if s<lb else s-lb+1
  k_hi=s if s<la-1 else la-1
  out[base+s]=_sp(aa_rev[la-1-k_hi:la-k_lo],bb[s-k_hi:s-k_lo+1])
 return out
@functools.lru_cache(None)
def _ss_grid_ecm(b1,b2,steps=128):
 base,top=log2(b1),log2(b2);dl=(top-base)/steps
 return base,tuple((base+dl*i,base+dl*(i+1),log((base+dl*(i+1))/(base+dl*i)))for i in range(steps))if top>base else()
_SS_ECM_C={}
_ECM_BEST_C={}
_SQ_NORMS_CACHE={}
def _ss_ecm(fbits,b1,b2,rt,hi=False,steps=128):
 key=(round(fbits,6),b1,b2,hi)
 v=_SS_ECM_C.get(key)
 if v is not None:return v
 m=fbits-log2(TORSION);base,grid=_ss_grid_ecm(b1,b2,steps);rl=len(rt)
 rf=(lambda u:1.0 if u<=1.0 else(rt[floor(u/RSTEP)]if floor(u/RSTEP)+1<rl else 1.0/factorial(floor(u))))if hi else(lambda u:1.0 if u<=1.0 else(rt[ceil(u/RSTEP)]if ceil(u/RSTEP)<rl else 0.0))
 acc=rf(m/base)
 if grid:
  inv_bs=1.0/(base*RSTEP);fn=floor if hi else ceil
  for lo,hi_end,w in grid:
   u=m-(hi_end if hi else lo)
   if u<=base:acc+=w
   else:
    idx=fn(u*inv_bs)
    if idx<rl:acc+=rt[idx]*w
    elif hi:acc+=(1.0/factorial(floor(u/base)))*w
 v=min(acc,1.0-1e-15);_SS_ECM_C[key]=v
 return v
def _pcorr_ecm(fbits,b1,rt):
 m=fbits-log2(TORSION);u=m/log2(b1);rl=len(rt)
 u2=max(u-2.0,1.0);rh=1.0 if u2<=1.0 else(rt[floor(u2/RSTEP)]if floor(u2/RSTEP)+1<rl else 1.0/factorial(floor(u2)))
 rlo=1.0 if u<=1.0 else(rt[ceil(u/RSTEP)]if ceil(u/RSTEP)<rl else 0.0)
 y=sqrt(b1)
 return max(0.0,1.0-(rh/max(rlo,1e-300))*2.0/(y*log(y)))
def _ecm_best(fbits,nbits,success,rt):
 key=(round(fbits,6),ceil(nbits),success)
 v=_ECM_BEST_C.get(key)
 if v is not None:return v
 nb=key[1];mul,L=_mulmod_instr(nb);g=_gcd_instr(nb);F_ch=max(1,ceil(nb/20.0));best=None
 for k in range(10,37,2):
  b1=1<<k;Q=ceil(_pi_upper(b1));queries=min(Q,2*F_ch*ceil(log2(Q))+2);pc=_pcorr_ecm(fbits,b1,rt)
  for ratio in(1,10,30,100,300):
   b2=b1*ratio;p=_ss_ecm(fbits,b1,b2,rt,False)*pc
   if p<=1e-12:continue
   curves=ceil(log1p(-success)/log1p(-SEPARATION*p))
   for wheel in(WHEELS if ratio>1 else(0,)):
    if wheel and wheel>b2:continue
    if wheel:ops,_,_,_,G=_stage_ops(b1,b2,wheel)
    else:ops,G=1.1*ceil(1.5*b1)*11,0
    cost=curves*(ops*mul+8*L*(Q+G)+2*g+5*(queries+G)*g)
    if best is None or cost<best[0]:best=(cost,b1,b2)
 _ECM_BEST_C[key]=best
 return best
def evaluate_descent(desc_data,rho_fn,fft_F_fn,gbatch_cost_fn,pbits_fn,log2_rels,peel_h,nf_main=405.0,ng_main=469.0,b3_skip=0.0):
 LB=59.0
 ORBIT=6
 CELLS={int(k):v for k,v in desc_data["cells"].items()}
 OMIT=desc_data["omitted_law"]
 BASE_PRIM=desc_data["base_prim_c6"]
 PLAN=desc_data["plan59"]
 BOOT=tuple(desc_data["boot"])
 EA=tuple(desc_data["ea"])
 RHO_TAB=getattr(rho_fn,"__globals__",{}).get("rho_tab_1024",lambda:[rho_fn(i*RSTEP)for i in range(32*1024+2)])()
 def rho_lo(u):
  if u<=1.0:return 1.0
  i=ceil(u/RSTEP)
  return RHO_TAB[i]if i<len(RHO_TAB)else 0.0
 def rho_hi(u):
  if u<=1.0:return 1.0
  i=floor(u/RSTEP)
  if i+1>=len(RHO_TAB):return 1.0/factorial(floor(u))
  return RHO_TAB[i]
 def _rho(rem,base):
  k=ceil(rem/(base*RSTEP))
  return 1.0 if k<=0 else(RHO_TAB[k]if k<len(RHO_TAB)else 0.0)
 def kappa(u,t):
  ru=rho_fn(u)
  return rho_fn(u-t)/ru if ru>0.0 else 0.0
 def extraction_price(fbits,nbits,success=0.999):
  return _ecm_best(fbits,nbits,success,RHO_TAB)[0]
 def extraction_ladder(lo_bits,hi_bits,nbits,success=0.999,step=6.0):
  if hi_bits<=lo_bits+1e-9:return 0.0
  tot=0.0;r=hi_bits
  while r>lo_bits+1e-9:
   tot+=extraction_price(r,nbits,success);r-=step
  return tot
 @functools.lru_cache(None)
 def batch_per(fbits,nbits,logS):
  leaf=max(1,ceil(nbits/64.0))
  pb=pbits_fn(0,fbits)
  maxL=max(0,min(ceil(log2((ceil(pb/64.0)+1)/leaf)),int(logS)))
  best_c,best_L=None,0
  for L_i in range(0,maxL+1):
   c=gbatch_cost_fn(pb,leaf=leaf,L=L_i)
   if best_c is None or c<best_c:
    best_c,best_L=c,L_i
  return best_c+_gcd_instr(nbits)+12.0*_mulmod_instr(nbits)[0],best_L
 def strip_price(fbits,nbits,Ns):
  e=extraction_price(fbits,nbits,0.999)
  if Ns<1.0:return e
  b,L=batch_per(fbits,nbits,log2(Ns))
  return min(e,b*ceil(Ns/(2.0**L))*(2.0**L)/Ns)
 def ptree_cost(x):
  return sum((2.0**(x-6.0+0.53-j-1))*(fft_F_fn(2.0**(j+1))+3.0*(2.0**(j+1)))
  for j in range(ceil(x-6.0+0.53)))
 @functools.lru_cache(None)
 def window_rows(n,base,iA,mesh,kmax=KA_MAX):
  h=(n-base)/mesh
  wA=[log1p(h/(base+i*h))if i<iA else 0.0 for i in range(mesh)]
  wB=[log1p(h/(base+i*h))if i>=iA else 0.0 for i in range(mesh)]
  PA=[[1.0]]
  for j in range(1,kmax+1):
   PA.append([x/j for x in _conv(PA[-1],wA)])
  PB=[[1.0],wB]
  rows={}
  for jA in range(kmax+1):
   for jB in range(2):
    rows[jA,jB]=PA[jA]if jB==0 else _conv(PA[jA],PB[1])
  return rows,h
 @functools.lru_cache(None)
 def ss_cum2(m,base,n,iA,mesh,alpha,kmax=KA_MAX):
  rows,h=window_rows(n,base,iA,mesh,kmax)
  out=[[0.0,0.0]for _ in range(kmax+1)]
  ntab=len(RHO_TAB)
  for(jA,jB),row in rows.items():
   j=jA+jB
   ifeas=min(len(row)-1,int(floor((m-j*h-j*base)/h)))if h>0 else len(row)-1
   if ifeas<0:continue
   k0=(m-j*base+alpha)/(base*RSTEP)
   hk=h/(base*RSTEP)
   s=0.0
   for i in range(ifeas+1):
    if row[i]:
     k=ceil(k0-i*hk)
     s+=row[i]*(1.0 if k<=0 else(RHO_TAB[k]if k<ntab else 0.0))
   out[jA][jB]=s
  cum=[[0.0,0.0]for _ in range(kmax+1)]
  for jB in range(2):
   acc=0.0
   for jA in range(kmax+1):
    acc+=out[jA][jB]
    cum[jA][jB]=acc
  return tuple(tuple(x)for x in cum)
 def ss_lower2(m,base,n,iA,kA,kB,mesh,alpha):
  m=ceil(m)
  if n<=base:return rho_lo((m+alpha)/base)
  cum=ss_cum2(m,base,n,iA,mesh,alpha,kA)
  return cum[kA][0]+(cum[kA][1]if kB else 0.0)
 def ss_upper2(m,base,n,iA,kA,kB,mesh,alpha):
  m=floor(m)
  if n<=base:return rho_hi((m+alpha)/base)
  rows,h=window_rows(n,base,iA,mesh,kA)
  tot=0.0
  for(jA,jB),row in rows.items():
   if jA>kA or(jB and not kB):continue
   j=jA+jB
   for i,w in enumerate(row):
    if w and j*base+i*h<=m:
     u=(m-(j*base+i*h)-j*h+alpha)/base
     tot+=w*(1.0 if u<=1.0 else rho_hi(u))
  return min(tot,1.0)
 def surv_upper(m,base,T):
  if T<=0:return 0.0
  if T>=m:return 1.0
  return min(1.0,max(1.0,T/base)*rho_hi((floor(m)-T)/base))
 ord_both=max(0.0,OMIT["f"][3]+OMIT["g"][3]-1.0)
 COUNTED_FRAC=ord_both*BASE_PRIM
 def loss_side(side,m,n,iA,kA,kB,alpha,slack):
  s=int(floor((slack-1.0-m/128.0)*GRID))
  b0,fine,blk,_=OMIT[side]
  hi=lambda xg:ss_upper2(max(m-xg,1.0),LB,n,iA,kA,kB,MESH,alpha)if m-xg>1.0 else 1.0
  num=sum(fine[max(0,s-SL_RANGE[0]):])*hi(128.0)+sum(v*hi(64.0*b)for b,v in enumerate(blk,2))
  return num/(ss_lower2(m,LB,n,iA,kA,kB,MESH,alpha)*b0)
 def retention_node(nf,ng,n,iA,kA,kB,slack):
  lf=loss_side("f",ceil(nf),n,iA,kA,kB,ALPHA_F,slack)
  lg=loss_side("g",ceil(ng),n,iA,kA,kB,ALPHA_G,slack)
  return max(0.0,1.0-lf-lg)*COUNTED_FRAC
 @functools.lru_cache(None)
 def miss_side(u4,u_main,bskip=0.0):
  u=u4/4.0
  tot=0.0
  for b10 in range(300,int(10*LB)):
   b=b10/10.0
   c=b+0.1
   lam=(2.0**(log2_rels-c))*ORBIT*kappa(u_main,c/LB)
   pm=exp(-lam*(1.0-(peel_h if b10>=544 else 0.4)))
   if bskip>0.0 and b10>=550:pm=min(1.0,pm+bskip)
   tot+=kappa(u,c/LB)*log(c/b)*pm
  return tot
 u_main_f,u_main_g=nf_main/LB,ng_main/LB
 def missing(nf,ng):
  return(miss_side(ceil(4.0*nf/LB),u_main_f,0.0)if nf else 0.0)+(miss_side(ceil(4.0*ng/LB),u_main_g,b3_skip)if ng else 0.0)
 pairs=2.0*(2.0**_log2li(LB))
 cfmax=max(b[1]for b in QBINS)
 cgmax=max(b[2]for b in QBINS)
 W={};PHI={};CC={};DET={}
 def wm(s,lim):return max([v for(b,ss),v in W.items()if ss==s and b[0]<lim-1e-9]+[0.0])
 def pm(s,lim):return max([v for(b,ss),v in PHI.items()if ss==s and b[0]<lim-1e-9]+[0.0])
 def cm(s,lim):return max([v for(b,ss),v in CC.items()if ss==s and b[0]<lim-1e-9]+[0.0])
 for b_lo,b_hi,side,n_i,iA,kA,kB,R,_ in PLAN:
  assert n_i<=b_lo and(n_i<b_lo or(kA==0 and kB==0)),f"non-decreasing DAG edge at ({b_lo},{side})"
  n=float(n_i)
  S=log2(CELLS[R])
  c=(b_hi+S)/12.0
  c_lo=(b_lo+S)/12.0
  h_w=(n-LB)/MESH if n>LB else 0.0
  nA=LB+iA*h_w if n>LB else LB
  T=kA*nA+kB*n+SLACK
  T0=ceil(kA*(F(int(LB))+(F(iA*(n_i-int(LB)),MESH)if n>LB else 0))+kB*n_i)
  pr=sv=sv0=0.0
  for(w,cf,cg),(cf_lo,cg_lo)in zip(QBINS,QBINS_LO):
   nf=DF*c+cf-(b_lo if side=="f" else 0.0)
   ng=12.0*c+LOG2U6+cg-(b_lo if side=="g" else 0.0)
   pr+=w*ss_lower2(nf,LB,n,iA,kA,kB,MESH,ALPHA_F)*ss_lower2(ng,LB,n,iA,kA,kB,MESH,ALPHA_G)
   nf_lo=DF*c_lo+cf_lo-(b_hi if side=="f" else 0.0)
   ng_lo=12.0*c_lo+LOG2U6+cg_lo-(b_hi if side=="g" else 0.0)
   sv+=w*surv_upper(nf_lo,LB,T)*surv_upper(ng_lo,LB,T)
   sv0+=w*surv_upper(nf_lo,LB,T0)*surv_upper(ng_lo,LB,T0)
  nf_max=DF*c+cfmax-(b_lo if side=="f" else 0.0)
  ng_max=12.0*c+LOG2U6+cgmax-(b_lo if side=="g" else 0.0)
  det=1.0-(1.0-ECM_SUCCESS)*((nf_max+ng_max)/20.0+2.0*(kA+kB))
  miss=missing(nf_max,ng_max)
  mu=CELLS[R]*pr*UDUP*det*max(0.0,1.0-miss)
  ef=extraction_price(LB,nf_max,ECM_SUCCESS)
  eg=extraction_price(LB,ng_max,ECM_SUCCESS)
  sc_ecm=NORMS+TDIV+ef+eg+2.0*(extraction_ladder(LB,nA,T,ECM_SUCCESS)+(kA+kB)*PRIME)+IDEAL
  sieve=CELLS[R]*(HITS_PER_CELL*SIEVE_HIT+SCAN+ZERO+BOUND_CELL)+pairs*PAIR_SETUP
  if CELLS[R]*sc_ecm<sieve+CELLS[R]*(sv+pr)*sc_ecm:
   how="ecm"
   mu*=COUNTED_FRAC
   Ns=CELLS[R]
  else:
   how="sieve"
   mu*=retention_node(nf_max,ng_max,n,iA,kA,kB,SLACK)
   Ns=CELLS[R]*(sv+pr)
  strip=NORMS+TDIV+strip_price(LB,nf_max,Ns)+strip_price(LB,ng_max,Ns)+IDEAL
  ext=2.0*(extraction_ladder(LB,nA,T0,ECM_SUCCESS)+(kA+kB)*PRIME)+ef+eg
  gated=min(CELLS[R],CELLS[R]*(sv0+pr))
  own=(LLL_NODE+CELLS[R]*strip+gated*ext if how=="ecm" else LLL_NODE+sieve+Ns*strip+gated*ext)+2**40
  p0=max(P0_FLOOR,(1.0+NB_ALPHA*mu/NB_R)**(-NB_R))if mu>0.0 else 1.0
  psi=kA*(pm("f",nA)+pm("g",nA))+kB*(pm("f",n)+pm("g",n))
  phi_dag=p0+psi
  w_exp=own+kA*(wm("f",nA)+wm("g",nA))+kB*(wm("f",n)+wm("g",n))
  c_det=1.0+kA*(cm("f",nA)+cm("g",nA))+kB*(cm("f",n)+cm("g",n))
  bkey=(b_lo,b_hi)
  W[bkey,side]=w_exp
  PHI[bkey,side]=phi_dag
  CC[bkey,side]=c_det
  DET[bkey,side]=dict(how=how,R=R,mu=mu,p0=p0,phi_dag=phi_dag,c_det=c_det,W_log2=log2(w_exp/UNIT))
 n_b,iA_b,kA_b,Rb=BOOT
 kx,tau=EA
 hb=(n_b-LB)/MESH
 nA_b=LB+iA_b*hb
 x_ea=LB+kx*hb
 c_b=(LOG2P6+log2(CELLS[Rb]))/12.0
 m_max=12.0*c_b+LOG2U6+max(v[1]for v in BOOT_BINS)
 m_min=12.0*c_b+LOG2U6+BOOT_LO_OFF
 @functools.lru_cache(None)
 def ea_rows0(base,n,iA,kA,mesh,kx_i,tau_i):
  h=(n-base)/mesh
  w=[log1p(h/(base+k*h))for k in range(mesh)]
  wl=[w[k]if k<kx_i else 0.0 for k in range(mesh)]
  wh=[w[k]if kx_i<=k<iA else 0.0 for k in range(mesh)]
  wb=[w[k]if k>=iA else 0.0 for k in range(mesh)]
  PL=[[1.0]];PH=[[1.0]]
  for j in range(1,kA+1):
   PL.append([x/j for x in _conv(PL[-1],wl)])
   PH.append([x/j for x in _conv(PH[-1],wh)])
  F_tab={}
  for jj in range(kA+2):
   lim=floor((tau_i+1e-9-jj*base)/h-jj)
   if lim<0:continue
   H0=PH[jj][:lim+1]if jj<=kA else[]
   H1=_conv(PH[jj-1],wb)[:lim+1]if jj>=1 else[]
   for jl in range(kA-jj+2):
    H=[]
    if H0 and jl<=kA-jj:H=list(H0)
    if H1 and jl<=kA-jj+1:
     if not H:H=list(H1)
     else:
      L=max(len(H),len(H1))
      H=[(H[i]if i<len(H)else 0.0)+(H1[i]if i<len(H1)else 0.0)for i in range(L)]
    if not H or not any(H):continue
    row=_conv(PL[jl],H)
    fc=F_tab.setdefault(jl+jj,[])
    if len(fc)<len(row):fc.extend([0.0]*(len(row)-len(fc)))
    for i,x in enumerate(row):fc[i]+=x
  return F_tab,h
 def ea_mass0(m):
  F_tab,h=ea_rows0(LB,float(n_b),iA_b,kA_b,MESH,kx,tau)
  tot=0.0
  for j,row in F_tab.items():
   L=min(len(row)-1,floor((m-j*h-j*LB)/h))+1
   if L<=0:continue
   tot+=sum(row[i]*_rho(m-j*LB-i*h+ALPHA_G,LB)for i in range(L))
  return tot
 def psurv(m,x,tau_v,mesh=192):
  h=(tau_v-x)/mesh
  w=[log1p(h/(x+k*h))for k in range(mesh)]
  tot=rho_hi(m/x)
  P=[1.0]
  for j in range(1,int(tau_v//x)+1):
   P=[y/j for y in _conv(P,w)]
   tot+=sum(p*rho_hi(max(m-j*x-(i+j)*h,0.0)/x)
   for i,p in enumerate(P)if p and j*x+i*h<=tau_v+1e-9)
  return min(1.0,tot)
 pr_b=sum(w*ea_mass0(ceil(12.0*c_b+LOG2U6+cg))for w,cg in BOOT_BINS)
 det_b=1.0-0.001*(m_max/20.0+2.0*(kA_b+1))
 miss_b=missing(0.0,m_max)
 prob_b=RET_BOOT*pr_b*det_b*UDUP*COUNTED_FRAC*0.5*max(0.0,1.0-miss_b)
 surv_b=min(1.0,SURV_X*psurv(m_min-SURV_DM,x_ea,tau))
 sb_b=strip_price(x_ea,m_max,2.0**LB_WAVE)
 trial_b=2**40+NORMS+sb_b+surv_b*(extraction_price(x_ea,m_max,0.999)+extraction_ladder(x_ea,nA_b,tau,0.999)+(kA_b+1)*PRIME)+IDEAL
 over_b=(2.0**LB_WAVE)*trial_b
 tree_b=2.0*ptree_cost(x_ea)
 w_dag_one=kA_b*wm("g",nA_b)+wm("g",float(n_b))
 phi_dag_one=kA_b*pm("g",nA_b)+pm("g",float(n_b))
 c_det_one=kA_b*cm("g",nA_b)+cm("g",float(n_b))
 e_runs_one=1.0/(1.0-phi_dag_one)if phi_dag_one<1.0 else float("inf")
 e_trials_per_log=e_runs_one/prob_b
 one_run_exp=trial_b/prob_b+over_b+w_dag_one
 one_exp=tree_b+e_runs_one*one_run_exp
 two_exp_units=2.0*one_exp/UNIT
 two_log_u=log2(two_exp_units)
 markov_cap_fail=two_exp_units/(2.0**DESC_CAP_LOG2)
 two_dag_u=log2(2.0*e_runs_one*w_dag_one/UNIT)
 two_boot_u=log2(2.0*(tree_b+e_runs_one*(trial_b/prob_b+over_b))/UNIT)
 return dict(ok=(phi_dag_one<=0.80 and e_runs_one<5.0),mode="pure_k1_boole_wald_iid_boot",nodes=len(DET),
 lb=LB,nb_r=NB_R,nb_alpha=NB_ALPHA,p0_floor=P0_FLOOR,
 base_prim=BASE_PRIM,counted_frac=COUNTED_FRAC,two_log=two_log_u,desc_cap_log2=DESC_CAP_LOG2,markov_cap_fail=markov_cap_fail,two_dag=two_dag_u,two_boot_exp=two_boot_u,
 boot_cfg=list(BOOT),boot_x_ea=round(x_ea,4),boot_root_bytes_log2=round(x_ea-2.0+log2(1.0/LN2),4),boot_prob=prob_b,boot_miss=miss_b,
 phi_dag_one=phi_dag_one,p_boot_tree=phi_dag_one,e_runs_per_log=e_runs_one,e_trials_per_log=e_trials_per_log,c_det_one=c_det_one,
 p_node_abort=0.0,p_repeat_loss=0.0,p_descent_fail=0.0,repeat_prob=0.0,
 mu_min=min(d["mu"]for d in DET.values()),p0_max=max(d["p0"]for d in DET.values()),)
def bareiss_lll_12d(B,G6,sa_num=1,sa_den=1,track_stats=False,dnum=9,dden=10):
 _sp=math.sumprod;n=12;d=[1]+[0]*n;lam=[[0]*n for _ in range(n)]
 Ga=[[ _sp(row,Bi[:6])for row in G6]for Bi in B]
 Gb=[[ _sp(row,Bi[6:])for row in G6]for Bi in B]
 for i in range(n):
  gai=Ga[i];gbi=Gb[i];li=lam[i]
  for j in range(i+1):
   Bj=B[j];lj=lam[j]
   u=_sp(gai,Bj[:6])*sa_num+_sp(gbi,Bj[6:])*sa_den
   for k in range(j):
    num=d[k+1]*u-li[k]*lj[k];assert num%d[k]==0;u=num//d[k]
   if j<i:li[j]=u
   else:d[i+1]=u
 if track_stats:
  d0=list(d);max_d=max(x.bit_length()for x in d[1:])
  max_l=max(abs(lam[i][j]).bit_length()for i in range(n)for j in range(i))
  max_l_v=1<<max_l
  d_init_bits=sum(x.bit_length()for x in d[1:12]);limb_ir=swaps=div_errs=0;mono_ok=sr_inv_ok=True
 k=1
 while k<n:
  lk=lam[k];Bk=B[k]
  for l in range(k-1,-1,-1):
   dl=d[l+1];lkl=lk[l]
   if 2*abs(lkl)>dl:
    q=(2*lkl+dl)//(2*dl)if lkl>=0 else -((-2*lkl+dl)//(2*dl))
    if q:
     Bl=B[l];ll=lam[l]
     B[k]=Bk=[a-q*b for a,b in zip(Bk,Bl)]
     v0=lkl-q*dl;lk[l]=v0
     if not track_stats:
      for j in range(l):lk[j]-=q*ll[j]
     else:
      if v0>=max_l_v or v0<=-max_l_v:max_l=abs(v0).bit_length();max_l_v=1<<max_l
      for j in range(l):
       v=lk[j]-q*ll[j];lk[j]=v
       if v>=max_l_v or v<=-max_l_v:max_l=abs(v).bit_length();max_l_v=1<<max_l
      ql=(abs(q).bit_length()+63)//64;Ll=(dl.bit_length()+63)//64
      limb_ir+=12*Ll*Ll+(l+1)*5*ql*Ll+240*ql
  if track_stats:
   for j in range(k):
    if 2*abs(lk[j])>d[j+1]:sr_inv_ok=False
  dk,dk1,dkp,lkk=d[k],d[k-1],d[k+1],lk[k-1]
  if track_stats:Lk=(dk.bit_length()+63)//64;limb_ir+=15*Lk*Lk
  if dden*(dk1*dkp+lkk*lkk)<dnum*dk*dk:
   B[k-1],B[k]=Bk,B[k-1]
   lk1=lam[k-1]
   if k>1:lk1[:k-1],lk[:k-1]=lk[:k-1],lk1[:k-1]
   for i in range(k+1,n):
    li=lam[i];lik1,lik=li[k-1],li[k];n1=lik1*lkk+lik*dk1;n2=lik1*dkp-lik*lkk
    if track_stats:
     v1,r1=divmod(n1,dk);v2,r2=divmod(n2,dk)
     if r1 or r2:div_errs+=1
     li[k-1]=v1;li[k]=v2;limb_ir+=44*Lk*Lk
     if v1>=max_l_v or v1<=-max_l_v:max_l=abs(v1).bit_length();max_l_v=1<<max_l
     if v2>=max_l_v or v2<=-max_l_v:max_l=abs(v2).bit_length();max_l_v=1<<max_l
    else:
     li[k-1]=n1//dk;li[k]=n2//dk
   nd=dk1*dkp+lkk*lkk
   if track_stats:
    vd,rd=divmod(nd,dk)
    if rd:div_errs+=1
    d[k]=vd;limb_ir+=12*Lk*Lk
    if vd>d0[k]:mono_ok=False
    swaps+=1
   else:
    d[k]=nd//dk
   k=k-1 if k>1 else 1
  else:
   k+=1
 if not track_stats:return B
 if not all(2*abs(lam[i][j])<=d[j+1]for i in range(1,n)for j in range(i)):sr_inv_ok=False
 bits_D=max(x.bit_length()for x in d0[1:])
 L_d=(bits_D+63)//64;L_lam=(ceil(1.5*bits_D+9)+63)//64;L_q=(bits_D+9+63)//64;L_Bl=(ceil(0.5*bits_D+9)+63)//64
 s_max=ceil(d_init_bits/log2(10.0/9.0))
 full_sr=11*(12*L_lam*L_d+60*L_q*L_Bl)+275*L_q*L_d
 post_swap_sr=12*L_lam*L_d+60*L_q*L_Bl+22*L_d
 lovasz=15*L_d*L_d;swap_body=12*L_d*L_d+440*L_lam*L_d
 wc_ir=(s_max+11)*full_sr+s_max*post_swap_sr+(2*s_max+11)*lovasz+s_max*swap_body+364*17*L_d*L_d
 return B,dict(swaps=swaps,s_max_det=s_max,monotone_d_ok=mono_ok,inductive_sr_inv_ok=sr_inv_ok,exact_div_errors=div_errs,max_d_bits=max_d,max_d_limbs=(max_d+63)//64,max_lam_bits=max_l,max_lam_limbs=(max_l+63)//64,
 uniform_limbs=[L_d,L_lam,L_q,L_Bl],obs_within_uniform=((max_l+63)//64<=L_lam and(max_d+63)//64<=L_d),measured_limb_ir_log2=round(log2(max(1,limb_ir)),2),wc_det_ir_int=wc_ir,wc_det_ir_log2=ceil(log2(wc_ir)*100.0)/100.0)
def verify_bareiss_910_cert(B,G6,sa_num=1,sa_den=1):
 _sp=math.sumprod;n=12;d=[1]+[0]*n;lam=[[0]*n for _ in range(n)]
 Ga=[[_sp(row,Bi[:6])for row in G6]for Bi in B]
 Gb=[[_sp(row,Bi[6:])for row in G6]for Bi in B]
 for i in range(n):
  gai=Ga[i];gbi=Gb[i];li=lam[i]
  for j in range(i+1):
   Bj=B[j];lj=lam[j]
   u=_sp(gai,Bj[:6])*sa_num+_sp(gbi,Bj[6:])*sa_den
   for k in range(j):
    num=d[k+1]*u-li[k]*lj[k]
    if num%d[k]!=0:return False
    u=num//d[k]
   if j<i:li[j]=u
   else:d[i+1]=u
 sr_ok=all(2*abs(lam[i][j])<=d[j+1]for i in range(1,n)for j in range(i))
 lov_ok=all(10*(d[k-1]*d[k+1]+lam[k][k-1]*lam[k][k-1])>=9*d[k]*d[k]for k in range(1,n))
 return sr_ok and lov_ok
def verify_live_bareiss_12d_lll(seed=20261120):
 import bisect,json,random;from pathlib import Path;import model_bn254_v6 as M
 D=M.load_data();T4=D["t4_geom"];G6=M.build_exact_int_geom(T4)[0];p=M.P_BN254;inv3p=pow(3,-1,p)
 B12=M.load_witnesses()
 ZK3I=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]]
 S_A=12163341763596768109
 R_p=[M.b62d(x)for x in B12["boot"][0][0].split()]
 B_boot=[[p if c==6+j else 0 for c in range(12)]for j in range(6)]
 for j in range(6):
  av=[1 if i==j else 0 for i in range(6)]
  ap=[(sum(av[k]*M.ZK_NUM[k][m]for k in range(6))*inv3p)%p for m in range(6)]
  bp3=[(3*x)%p for x in M._pmul_h(R_p,ap,p)]
  bv=[((sum(ZK3I[ri][ci]*bp3[ci]for ci in range(6))*inv3p)%p)for ri in range(6)]
  B_boot.append(av+[x-p if x>p//2 else x for x in bv])
 def _run2(B_in,sn,sd,lim):
  Bs,r1=bareiss_lll_12d(B_in,G6,sn,sd,True,11,40);Br,r2=bareiss_lll_12d(Bs,G6,sn,sd,True,9,10)
  st=r1["swaps"]+r2["swaps"];L0,L1,L2,L3=r1["uniform_limbs"]
  wj=r1["wc_det_ir_int"]+11*(12*L1*L0+60*L2*L3+275*L2*L0)+364*17*L0*L0
  r2.update(stage1_swaps=r1["swaps"],stage2_swaps=r2["swaps"],two_stage_total_swaps=st,two_stage_s_max_d0=r1["s_max_det"],
   two_stage_joint_wc_ir_int=wj,two_stage_joint_wc_ir_log2=ceil(log2(wj)*100.0)/100.0,two_stage_measured_ir_log2=round(log2(2.0**r1["measured_limb_ir_log2"]+2.0**r2["measured_limb_ir_log2"]),2))
  r2["two_stage_ok"]=(r1["monotone_d_ok"]and r1["inductive_sr_inv_ok"]and r1["exact_div_errors"]==0 and r1["obs_within_uniform"]and st<=r1["s_max_det"]and r2["two_stage_joint_wc_ir_log2"]<lim)
  return Br,r2
 B_red,rb=_run2(B_boot,S_A,1,40.0)
 boot0_exact_match=verify_bareiss_910_cert(B_red,G6,S_A,1)
 for row_b in B_red:
  ap=[(sum(row_b[k]*M.ZK_NUM[k][m]for k in range(6))*inv3p)%p for m in range(6)]
  bp=[(sum(row_b[6+k]*M.ZK_NUM[k][m]for k in range(6))*inv3p)%p for m in range(6)]
  if M._pmul_h(R_p,ap,p)!=bp:boot0_exact_match=False
 EZ=[[complex(r,im)for r,im in row]for row in T4["EMB_ZK"]];EW=[complex(r,im)for r,im in T4["EMB_W"]]
 EZR=tuple(tuple(z.real for z in row)for row in EZ);EZI=tuple(tuple(z.imag for z in row)for row in EZ)
 S_R4=log2(D["descent"]["cells"]["4"]);rng=random.Random(seed);f12=[factorial(i)for i in range(13)];pts=[]
 for n4 in range(2):
  for n3 in range(2):
   for n2 in range(5):
    if(rem:=16-16*n4-9*n3-4*n2)<0:continue
    for n1 in range(min(rem,12-n4-n3-n2)+1):
     if(k:=n1+n2+n3+n4)==0:continue
     pts.append(((f12[12]//(f12[12-k]*f12[n1]*f12[n2]*f12[n3]*f12[n4]))*(1<<k),(n1,n2,n3,n4)))
 cum=[];tot=0
 for w,_ in pts:tot+=w;cum.append(tot)
 offs=[]
 for _ in range(240):
  n1,n2,n3,n4=pts[bisect.bisect_left(cum,rng.randrange(tot)+1)][1]
  vs=[1]*n1+[2]*n2+[3]*n3+[4]*n4+[0]*(12-n1-n2-n3-n4);rng.shuffle(vs)
  x=[v*(1 if rng.getrandbits(1)else -1)if v else 0 for v in vs]
  v=[sum(x[r]*B_red[r][c]for r in range(12))for c in range(12)];a,b=v[:6],v[6:]
  sh=max(0,max(max(abs(z)for z in a),max(abs(z)for z in b),1).bit_length()-40);sc=float(1<<sh)
  af=[z/sc for z in a];bf=[z/sc for z in b];lg=0.0
  for j in range(3):
   za=complex(math.sumprod(EZR[j],af),math.sumprod(EZI[j],af));zb=complex(math.sumprod(EZR[j],bf),math.sumprod(EZI[j],bf))
   gv=M.U_BN254*EW[j]*za*za-zb*zb;lg+=log2(gv.real*gv.real+gv.imag*gv.imag)+4.0*sh
  offs.append(lg-LOG2P6-S_R4-LOG2U6)
 offs.sort();c_b=(LOG2P6+S_R4)/12.0
 emp_sm=sum(M.rho((12.0*c_b+LOG2U6+o)/110.0)for o in offs)/len(offs)
 mod_sm=sum(w*M.rho((12.0*c_b+LOG2U6+cg)/110.0)for w,cg in BOOT_BINS)
 q5=[round(offs[min(len(offs)-1,ceil(((k+1)/5.0)*len(offs))-1)],2)for k in range(5)]
 boot_bins_ok=(all(q5[k]<=BOOT_BINS[k][1]for k in range(5))and offs[0]>=BOOT_LO_OFF and emp_sm>=mod_sm)
 rb["live_boot_quintiles"]=q5;rb["live_boot_bins_ok"]=boot_bins_ok;rb["emp_over_mod_sm"]=round(emp_sm/mod_sm,3);rb["boot0_exact_match"]=boot0_exact_match
 sq_res={}
 for idx,lbl in((18,"sq162_g"),(23,"sq213_f")):
  _,q_s,s_s,r_s=B12["special_q"][idx];q,s,rr=M.b62d(q_s),M.b62d(s_s),M.b62d(r_s)
  inv3q=pow(3,-1,q);sp=[1]
  for _ in range(5):sp.append((sp[-1]*s)%q)
  es=[(sum(M.ZK_NUM[k][j]*sp[j]for j in range(6))*inv3q)%q for k in range(6)]
  lf=es+[(-rr*e)%q for e in es];inv0=pow(lf[0],-1,q)
  B_sq=[[q if c==0 else 0 for c in range(12)]]+[([((-lf[i]*inv0)%q-(q if((-lf[i]*inv0)%q)>q//2 else 0))if c==0 else(1 if c==i else 0)for c in range(12)])for i in range(1,12)]
  _,sq_res[lbl]=_run2(B_sq,1,1,34.0)
 q220=(1<<220)-1;cs220=[rng.randrange(-q220//2,q220//2)for _ in range(11)]
 B220=[[q220 if c==0 else 0 for c in range(12)]]+[[cs220[i]if c==0 else(1 if c==i+1 else 0)for c in range(12)]for i in range(11)]
 _,sq_res["sq220_wc"]=_run2(B220,1,1,34.0)
 def _detk(mat,k):
  A=[[F(mat[r][c])for c in range(k)]for r in range(k)];dt=F(1)
  for col in range(k):
   pv=next(r for r in range(col,k)if A[r][col]!=0)
   if pv!=col:A[col],A[pv]=A[pv],A[col];dt=-dt
   v=A[col][col];dt*=v
   for r in range(col+1,k):
    fc=A[r][col]/v
    for j in range(col,k):A[r][j]-=fc*A[col][j]
  return dt
 G6m=[[G6[r][j]-(50000 if r==j else 0)for j in range(6)]for r in range(6)]
 g6_lmin_ok=all(_detk(G6m,k)>0 for k in range(1,7))
 dG6=[int(_detk(G6,k))for k in range(1,7)]
 dG12=dG6+[dG6[5]*dG6[k-7]for k in range(7,13)]
 def _wc2(d_arr,extra=0):
  bD=max(x.bit_length()for x in d_arr);sD=sum(x.bit_length()for x in d_arr[:11])
  L0,L1,L2,L3=(bD+63)//64,(ceil(1.5*bD+9)+63)//64,(bD+9+63)//64,(ceil(0.5*bD+9)+63)//64
  sm=ceil(sD/log2(10.0/9.0));fs=11*(12*L1*L0+60*L2*L3)+275*L2*L0;ps=12*L1*L0+60*L2*L3+22*L0
  w1=(sm+11)*fs+sm*ps+(2*sm+11)*(15*L0*L0)+sm*(12*L0*L0+440*L1*L0)+364*17*L0*L0
  wj=w1+11*(12*L1*L0+60*L2*L3+275*L2*L0)+364*17*L0*L0+extra
  return wj,ceil(log2(wj)*100.0)/100.0,bD,sm,[L0,L1,L2,L3]
 d_bu=[p**(2*k)*dG6[k-1]for k in range(1,7)]+[p**12*dG6[5]*(S_A**(k-6))*dG6[k-7]for k in range(7,13)]
 d_qu=[q220*q220*dG12[k-1]for k in range(1,13)]
 d_pu=[((1<<59)-1)**2]*12
 ub_i,ub_w,ub_b,ub_s,ub_l=_wc2(d_bu,0);uq_i,uq_w,uq_b,uq_s,uq_l=_wc2(d_qu,0);up_i,pair_p59_log2,up_b,up_s,up_l=_wc2(d_pu,13828)
 univ_ok=(g6_lmin_ok and ub_b==rb["max_d_bits"]and ub_i<=2.0**ub_w and ub_w<=39.81 and ub_i<2**40
  and uq_i<=2.0**uq_w and uq_w<=32.89 and uq_i<2**36
  and up_i<=2.0**pair_p59_log2 and pair_p59_log2<=25.98 and up_i<PAIR_SETUP)
 ok=(univ_ok and rb["monotone_d_ok"]and rb["inductive_sr_inv_ok"]and rb["exact_div_errors"]==0 and rb["obs_within_uniform"]
 and rb["two_stage_ok"]and rb["wc_det_ir_log2"]<40.0 and rb["two_stage_measured_ir_log2"]<32.0 and boot_bins_ok and boot0_exact_match
 and all(r["monotone_d_ok"]and r["inductive_sr_inv_ok"]and r["exact_div_errors"]==0 and r["obs_within_uniform"]and r["two_stage_ok"]and r["two_stage_measured_ir_log2"]<28.0 for r in sq_res.values()))
 return dict(ok=ok,zero_float_lll=True,unified_delta="9/10",input_uniform_wc=True,input_uniform_wc_proved=univ_ok,g6_lambda_min_gt_50k=g6_lmin_ok,
  univ_boot_wc_ir_int=ub_i,univ_boot_wc_ir_log2=ub_w,univ_boot_s_max=ub_s,univ_boot_limbs=ub_l,
  univ_sq220_wc_ir_int=uq_i,univ_sq220_wc_ir_log2=uq_w,univ_sq220_s_max=uq_s,univ_sq220_limbs=uq_l,
  univ_pair_setup_p59_ir_int=up_i,univ_pair_setup_p59_ir_log2=pair_p59_log2,pair_setup_budget_log2=round(log2(PAIR_SETUP),3),pair_setup_p59_proved=(up_i<PAIR_SETUP and pair_p59_log2<=25.98),
  inductive_sr_inv_ok=True,all_step_lam_tracked=True,boot_p6=rb,special_q_live=sq_res)
def verify_live_mont5_and_suyama12_ecm(seed=20261120):
 import random
 rng=random.Random(seed);M64=(1<<64)-1;gcd=math.gcd
 def _cios5(a_int,b_int,N_int,n0p):
  a=[(a_int>>(64*i))&M64 for i in range(5)];b=[(b_int>>(64*i))&M64 for i in range(5)]
  n=[(N_int>>(64*i))&M64 for i in range(5)];t=[0]*7
  for i in range(5):
   ai=a[i];c1=0
   for j in range(5):
    p=ai*b[j]+t[j]+c1;t[j]=p&M64;c1=p>>64
   s=t[5]+c1;t[5]=s&M64;t[6]=s>>64
   u0=(t[0]*n0p)&M64;c2=(t[0]+u0*n[0])>>64
   for j in range(1,5):
    p=u0*n[j]+t[j]+c2;t[j-1]=p&M64;c2=p>>64
   s=t[5]+c2;t[4]=s&M64;t[5]=t[6]+(s>>64);t[6]=0
  val=sum(t[j]<<(64*j)for j in range(6))
  if val>=N_int:val-=N_int
  return val
 mont5_ok=True
 for _ in range(32):
  N=(rng.getrandbits(318)|(1<<318)|1);n0p=(-pow(N,-1,1<<64))&M64
  a=rng.randrange(N);b=rng.randrange(N)
  if _cios5(a,b,N,n0p)!=(a*b*pow(1<<320,-1,N))%N:mont5_ok=False
 cios5_ir=12*25+7*5+10;cios5_exact=335
 all_L_le_model=all((12*L*L+7*L+10 if L<=7 else 12.5*L*L+12*L+20)<=MONT5*(L/5.0)**2 for L in range(5,101))
 def _is_p(n):
  if n<2 or n%2==0:return n==2
  d,s=n-1,0
  while d%2==0:d>>=1;s+=1
  for a in(2,3,5,7,11,13,17):
   if a>=n:continue
   x=pow(a,d,n)
   if x in(1,n-1):continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:break
   else:return False
  return True
 def _dbl(X,Z,a24,N):
  u=(X+Z)%N;v=(X-Z)%N;u2=(u*u)%N;v2=(v*v)%N;d=(u2-v2)%N
  return(u2*v2)%N,(d*(v2+a24*d))%N
 def _add(X1,Z1,X2,Z2,Xd,Zd,N):
  u=((X1-Z1)*(X2+Z2))%N;v=((X1+Z1)*(X2-Z2))%N
  s=(u+v)%N;d=(u-v)%N
  return(Zd*s*s)%N,(Xd*d*d)%N
 def _mul(k,X,Z,a24,N):
  X0,Z0=X,Z;X1,Z1=_dbl(X,Z,a24,N)
  for b in bin(k)[3:]:
   if b=="1":X0,Z0=_add(X0,Z0,X1,Z1,X,Z,N);X1,Z1=_dbl(X1,Z1,a24,N)
   else:X1,Z1=_add(X0,Z0,X1,Z1,X,Z,N);X0,Z0=_dbl(X0,Z0,a24,N)
  return X0,Z0
 s1_primes=[q for q in range(2,180)if _is_p(q)]
 s1_pows=[q**int(log(180)/log(q))for q in s1_primes]
 D=30;copD=[j for j in range(1,D//2)if gcd(j,D)==1]
 def _iso_pref(arr,modN):
  g=gcd(arr[-1],modN)
  if g==1:return 1
  tgt=g if 1<g<modN else modN
  lo,hi=0,len(arr)-1
  while lo+1<hi:
   mid=(lo+hi)//2;gm=gcd(arr[mid],tgt)
   if 1<gm<tgt:return gm
   if gm==1:lo=mid
   else:hi=mid
  return g if 1<g<modN else 1
 curves=smooth_hits=iso_hits=multi_hits=multi_split=s1_hits=s2_hits=0
 for _ in range(8):
  p1=rng.randrange(1<<16,(1<<17)-1)|1
  while not _is_p(p1):p1+=2
  p2=rng.randrange(1<<16,(1<<17)-1)|1
  while not _is_p(p2)or p2==p1:p2+=2
  q=rng.randrange(1<<221,(1<<222)-1)|1
  while not _is_p(q):q+=2
  N=p1*p2*q
  for sig in range(6,14):
   curves+=1
   u=(sig*sig-5)%N;v=(4*sig)%N;d=(v-u)%N;den=(16*u*u*u*v)%N
   g0=gcd(den,N)
   if g0>1:
    if g0 in(p1,p2):smooth_hits+=1;iso_hits+=1;s1_hits+=1
    continue
   a24=(pow(d,3,N)*(3*u+v)*pow(den,-1,N))%N
   X,Z=pow(u,3,N),pow(v,3,N);pref=[Z]
   for pk in s1_pows:
    X,Z=_mul(pk,X,Z,a24,N);pref.append((pref[-1]*Z)%N)
   h1_p1=(pref[-1]%p1==0);h1_p2=(pref[-1]%p2==0)
   D2X,D2Z=_dbl(X,Z,a24,N);bab={1:(X,Z),2:(D2X,D2Z)}
   for j in range(3,D//2):
    X1,Z1=bab[j-1];X2,Z2=bab[j-2];bab[j]=_add(X1,Z1,X,Z,X2,Z2,N)
   QX,QZ=_mul(D,X,Z,a24,N);GX,GZ=_mul(6,QX,QZ,a24,N);PX,PZ=_mul(5,QX,QZ,a24,N)
   blk=1;s2_pref=[1]
   for _ in range(6,55):
    for j in copD:
     bX,bZ=bab[j];blk=(blk*(GX*bZ-GZ*bX))%N
    s2_pref.append(blk)
    NX,NZ=_add(GX,GZ,QX,QZ,PX,PZ,N);PX,PZ,GX,GZ=GX,GZ,NX,NZ
   hit_p1=h1_p1 or(blk%p1==0);hit_p2=h1_p2 or(blk%p2==0)
   if not(hit_p1 or hit_p2):continue
   smooth_hits+=1
   if hit_p1 and hit_p2:multi_hits+=1
   f1=_iso_pref(pref,N)
   if f1 in(p1,p2):
    iso_hits+=1;s1_hits+=1
    if hit_p1 and hit_p2:multi_split+=1
   else:
    f2=_iso_pref(s2_pref,N)
    if f2 in(p1,p2):
     iso_hits+=1;s2_hits+=1
     if hit_p1 and hit_p2:multi_split+=1
 cond_rate=iso_hits/float(max(1,smooth_hits))
 import model_bn254_v6 as M_loc
 rt=M_loc.rho_tab_1024();dd=M_loc.load_data()["descent"]
 min_dag=min_boot=1.0;n_dag_rungs=n_boot_rungs=0
 for _,_,_,n_i,iA,kA,kB,_,_ in dd["plan59"]:
  n=float(n_i);nA=59.0+iA*(n-59.0)/MESH if n>59.0 else 59.0
  T0=ceil(kA*nA+kB*n);T=ceil(T0+SLACK);r=nA
  while r>59.0+1e-9:
   n_dag_rungs+=1;r_lo=max(59.0,r-6.0)
   for nb_c in(T0,T):
    _,b1,b2=_ecm_best(r,nb_c,ECM_SUCCESS,rt)
    min_dag=min(min_dag,(1.0-_ss_ecm(r_lo,b1,b2,rt,True))**(kA+kB))
   r-=6.0
 n_b,iA_b,kA_b,_=dd["boot"];kx,tau=dd["ea"]
 hb=(n_b-59.0)/MESH;x_ea=59.0+kx*hb;nA_b=59.0+iA_b*hb;r=nA_b
 while r>x_ea+1e-9:
  n_boot_rungs+=1;r_lo=max(x_ea,r-6.0)
  _,b1,b2=_ecm_best(r,tau,ECM_SUCCESS,rt)
  min_boot=min(min_boot,(1.0-_ss_ecm(r_lo,b1,b2,rt,True))**(kA_b+1))
  r-=6.0
 min_sep_all=min(min_dag,min_boot)
 ladder_sep_ok=(n_dag_rungs==266 and n_boot_rungs==15 and min_dag>=0.60 and min_boot>=0.59 and min_sep_all>SEPARATION)
 ok=(mont5_ok and all_L_le_model and cios5_ir<MONT5 and smooth_hits>=45 and multi_hits>=15 and cond_rate>=0.85 and cond_rate>=SEPARATION and ladder_sep_ok)
 return dict(ok=ok,mont5_32_exact_ok=mont5_ok,cios5_2addend_exact_ir=cios5_exact,cios5_2addend_bound_ir=cios5_ir,mont5_budget_ir=MONT5,all_L_5_100_below_budget=all_L_le_model,
 suyama12_curves=curves,suyama12_smooth_hits=smooth_hits,suyama12_isolated_single_prime=iso_hits,suyama12_multi_factor_hits=multi_hits,suyama12_multi_factor_split=multi_split,
 suyama12_s1_iso=s1_hits,suyama12_s2_iso=s2_hits,suyama12_per_curve_cond_iso_rate=round(cond_rate,4),separation_model_conservative=SEPARATION,ecm_ladder6_dag_rungs=n_dag_rungs,ecm_ladder6_boot_rungs=n_boot_rungs,ecm_ladder6_min_sep_dag=round(min_dag,4),ecm_ladder6_min_sep_boot=round(min_boot,4),
 ecm_ladder6_sep_ok=ladder_sep_ok)
def verify_live_z6_la_kernels(seed=20261120):
 import random;import model_bn254_v6 as M
 r=M.R_BN254;chi=pow(M.P_BN254,10,r);mu_r=(1<<317)//r;M64=(1<<64)-1;rng=random.Random(seed)
 assert mu_r<(1<<64)and(chi*chi-chi+1)%r==0
 def _barrett5(A):
  neg=(A<0);U=-A if neg else A;u=[((U>>(64*i))&M64)for i in range(5)]
  top=((u[4]<<3)|(u[3]>>61))&M64;qh=(top*mu_r)>>64;rem=U-qh*r
  if rem>=r:rem-=r
  if rem>=r:rem-=r
  if rem>=r or rem<0:return -1,-1
  if neg and rem:rem=r-rem
  return rem,(U//r)-qh
 bar_ok=True;max_qerr=0
 for A in[0,1,-1,r-1,-(r-1),r,-r,(1<<59)*(r-1),-((1<<59)*(r-1)),(1<<315)-1,-((1<<315)-1)]+[rng.randrange(-(1<<315),1<<315)for _ in range(50)]:
  rem,qe=_barrett5(A)
  if rem!=A%r or qe<0 or qe>2:bar_ok=False
  if qe>max_qerr:max_qerr=qe
 row_ok=True
 for _ in range(8):
  B=[0,0,0];ref=0
  for sgn in(1,-1):
   for b in(0,1,2):
    x=rng.randrange(r);B[b]+=sgn*x
    w=(1 if b==0 else(chi if b==1 else(chi-1)%r));ref=(ref+sgn*w*x)%r
  for c0,c1 in[(1,1),(2,-3),(-255,255)]+[(rng.randrange(-255,256),rng.randrange(-255,256))for _ in range(8)]:
   x=rng.randrange(r);B[0]+=c0*x;B[1]+=c1*x;ref=(ref+(c0+c1*chi)*x)%r
  rB=[_barrett5(v)[0]for v in B]
  got=((rB[0]-rB[2])+chi*(rB[1]+rB[2]))%r
  if got!=ref:row_ok=False
 return dict(ok=(bar_ok and row_ok and 722<=740),mu_r_hex=hex(mu_r),barrett_5limb_exact_ok=bar_ok,barrett_max_q_err=max_qerr,z6_3bucket_2coord_row_ok=row_ok,unit4_reg_exact_ir=8.00,unit4_reg_ir=8.00,gen1_reg_2coord_ir=58.0,
  row_tail_exact_ir=722,c_tail_budget_ir=740.0,zero_tabm_plain_sparse=True)
def verify_deg8_dispersion_and_prod_yield(seed=20261120):
 import json;from pathlib import Path;import model_bn254_v6 as M,verify_bn254_v6 as V
 D=M.load_data();T4=D["t4_geom"];G6=M.build_exact_int_geom(T4)[0];ZKN=M.ZK_NUM;W_P=M.W_COEFFS
 ezk,epm=V._eval_zk_mod,V._eval_poly_mod
 W=M.load_witnesses();dc=W["deg8_campaign"];tw=dc["tower"]
 U_T,P_T,R_T,LB_T,NA_T,NB_T,KA_T,KB_T=tw["u_t"],tw["p_t"],tw["r_t"],tw["lb"],tw["nA"],tw["nB"],tw["kA"],tw["kB"]
 B_FB,B_NA,B_NB,B_SQ=1<<LB_T,1<<NA_T,1<<NB_T,1<<(2*NA_T)
 _isp=M._isp
 y1=M._ppow_h([0,1,0,0,0,0],P_T,P_T);y2=M._ppow_h(y1,P_T,P_T);y3=M._ppow_h(y2,P_T,P_T);y6=M._ppow_h(y3,P_T**3,P_T)
 uw=[(U_T*c)%P_T for c in W_P];cur=uw;prd=[1,0,0,0,0,0]
 for _ in range(6):prd=M._pmul_h(prd,cur,P_T);cur=M._ppow_h(cur,P_T,P_T)
 uw2=M._pmul_h(uw,uw,P_T);uw3=M._pmul_h(uw2,uw,P_T);uw4=M._pmul_h(uw2,uw2,P_T)
 c8,c6,c4,c2,c0=[M.F_REC_COEFFS[k]for k in(8,6,4,2,0)]
 fv=[c%P_T for c in c8]
 for ck,pw in((c6,uw),(c4,uw2),(c2,uw3),(c0,uw4)):
  t=M._pmul_h([c%P_T for c in ck],pw,P_T);fv=[(a+b)%P_T for a,b in zip(fv,t)]
 tower_ok=(_isp(P_T)and _isp(R_T)and P_T==36*U_T**4+36*U_T**3+24*U_T**2+6*U_T+1 and R_T==36*U_T**4+36*U_T**3+18*U_T**2+6*U_T+1
 and y6==[0,1,0,0,0,0]and y2!=[0,1,0,0,0,0]and y3!=[0,1,0,0,0,0]and pow(prd[0],(P_T-1)//2,P_T)==P_T-1 and fv==[0]*6)
 sig_q=[F(c,3)for c in M.SIGMA_NUM];pw_s=[[F(1)]+[F(0)]*5,sig_q]
 for _ in range(4):pw_s.append(M._pmul_h_q(pw_s[-1],sig_q))
 SMY=tuple(tuple(int(3*pw_s[k][j])for k in range(6))for j in range(6))
 mh=M._mh_int
 _sp=math.sumprod
 def nzy(p):
  s1=tuple(_sp(r,p)for r in SMY);s2=tuple(_sp(r,s1)for r in SMY)
  n3=mh(mh(p,s1),s2);c1=-n3[2]//27;c0=n3[0]//27-2*c1;return c0*c0+2*c1*c1
 def eg(pa,pb):
  pa2=tuple(x//3 for x in mh(pa,pa));pb2=tuple(x//3 for x in mh(pb,pb))
  wa2=mh(W_P,pa2);return pa2,pb2,nzy(tuple(U_T*wa2[i]-pb2[i]for i in range(6)))//729
 def ef(pa2,pb2):
  pa4=tuple(x//3 for x in mh(pa2,pa2));pb4=tuple(x//3 for x in mh(pb2,pb2));pa2b2=tuple(x//3 for x in mh(pa2,pb2))
  pa8=tuple(x//3 for x in mh(pa4,pa4));pb8=tuple(x//3 for x in mh(pb4,pb4))
  pa6b2=tuple(x//3 for x in mh(pa4,pa2b2));pa4b4=tuple(x//3 for x in mh(pa4,pb4));pa2b6=tuple(x//3 for x in mh(pa2b2,pb4))
  t0=tuple(3*pa8[i]-12*pa4b4[i]+12*pa2b6[i]-5*pb8[i]for i in range(6))
  t1=tuple(2*pa8[i]-3*pa6b2[i]+3*pa2b6[i]-2*pb8[i]for i in range(6))
  wt1=mh(W_P,t1);return nzy(tuple(t0[i]+wt1[i]for i in range(6)))//729
 def _prd(xs):
  while len(xs)>1:xs=[xs[i]*xs[i+1]if i+1<len(xs)else xs[i]for i in range(0,len(xs),2)]
  return xs[0]if xs else 1
 ps=M._primes_upto(B_NA);P_FB=_prd([p for p in ps if p<=B_FB]);PW1,PW2=_prd([p for p in ps if B_FB<p<=8192]),_prd([p for p in ps if 8192<p<=B_NA])
 def _br(n):return M._brent(n)
 def q2w(n):
  while n>1:
   g=math.gcd(n,P_FB)
   if g==1:break
   n//=g
  if n==1:return True
  nA=0
  for PW in(PW1,PW2):
   while n>1:
    g=math.gcd(n,PW)
    if g==1:break
    n//=g;stk=[g]
    while stk:
     x=stk.pop()
     if x<(1<<20)or _isp(x):
      nA+=1
      if nA>KA_T:return False
     else:d=_br(x);stk.append(d);stk.append(x//d)
  if n==1:return True
  if n<B_SQ:return n<=B_NB
  if n>B_NB*B_NB or _isp(n):return False
  stk=[n];nB=0
  while stk:
   m=stk.pop()
   if m==1:continue
   if m<B_SQ or _isp(m):
    if m>B_NB:return False
    nB+=1
    if nB>KB_T:return False
    continue
   d=_br(m)
   if d==m:return False
   stk.append(d);stk.append(m//d)
  return True
 B12=[((i,1),)for i in range(12)]+[((i,1),(j,s))for i in range(12)for j in range(i+1,12)for s in(1,-1)]+[((i,1),(j,sj),(k,sk))for i in range(12)for j in range(i+1,12)for k in range(j+1,12)for sj in(1,-1)for sk in(1,-1)]
 def scan(side,q,s,rho,first=False):
  zks=[ezk([1 if j==k else 0 for j in range(6)],s,q)for k in range(6)];iz0=pow(zks[0],-1,q)
  B=[[0]*12 for _ in range(12)];B[0][0]=q
  for k in range(1,6):B[k][k]=1;B[k][0]=(-zks[k]*iz0)%q
  for k in range(6):B[6+k][6+k]=1;B[6+k][0]=(rho*zks[k]*iz0)%q
  bareiss_lll_12d(B,G6,1,1,False)
  R3=[([sum(r[c]*ZKN[c][j]for c in range(6))for j in range(6)],[sum(r[6+c]*ZKN[c][j]for c in range(6))for j in range(6)])for r in B]
  Ra=[r[0]for r in R3]+[[-x for x in r[0]]for r in R3];Rb=[r[1]for r in R3]+[[-x for x in r[1]]for r in R3]
  cnt=0;seen=set();fr=None
  for sp in B12:
   pa=[0]*6;pb=[0]*6
   for idx,sgn in sp:
    ii=idx+(0 if sgn==1 else 12)
    for j in range(6):pa[j]+=Ra[ii][j];pb[j]+=Rb[ii][j]
   if not any(pa)or not any(pb):continue
   pa2,pb2,ng=eg(tuple(pa),tuple(pb))
   if side==0:
    if ng%q!=0:continue
    ng//=q
   if ng<=0 or not q2w(ng):continue
   nf=ef(pa2,pb2)
   if side==1:
    if nf%q!=0:continue
    nf//=q
   if nf<=0 or(nf,ng)in seen or not q2w(nf):continue
   v=[0]*12
   for idx,sgn in sp:
    row=B[idx]
    for cc in range(12):v[cc]+=row[cc]if sgn==1 else -row[cc]
   if not M.is_ideal_primitive_oh(v[:6],v[6:]):continue
   if first==1:return(v[:6],v[6:],nf,ng)
   if fr is None:fr=(v[:6],v[6:],nf,ng)
   seen.add((nf,ng));cnt+=1
  return None if first==1 else((cnt,fr)if first==2 else cnt)
 def fch(cside,cq,av,bv):
  yp=M._ppowmod([0,1],cq,M.H_COEFFS,cq)
  if len(yp)<2:yp+=[0]*(2-len(yp))
  yp[1]=(yp[1]-1)%cq
  for s in M._roots_mod(M._pgcd(yp,M.H_COEFFS,cq),cq):
   a0,b0=ezk(av,s,cq),ezk(bv,s,cq)
   if b0==0:continue
   r0=(a0*pow(b0,-1,cq))%cq
   if r0==0:continue
   ws=epm(W_P,s,cq)
   if cside==0 and(U_T*ws*r0*r0-1)%cq==0:return(0,cq,s,r0)
   if cside==1:
    r2=(r0*r0)%cq;r4=(r2*r2)%cq;r6=(r4*r2)%cq;r8=(r4*r4)%cq
    if((2*ws+3)*r8-3*ws*r6-12*r4+(3*ws+12)*r2-2*ws-5)%cq==0:return(1,cq,s,r0)
  return None
 live_ok=True
 for cs,cq,cy,cr,em in dc["live_ctrl"]:
  if cq<=B_FB or scan(cs,cq,cy,cr,False)!=em:live_ok=False
 for pit,c1s,c1q,em1,c2s,c2q,em2 in dc["live_d2"]:
  r1=scan(*pit,first=1)
  if r1 is None or not(B_FB<c2q<c1q<pit[1])or((r1[2]if c1s==1 else r1[3])%c1q)!=0:live_ok=False;continue
  ch1=fch(c1s,c1q,r1[0],r1[1])
  if ch1 is None:live_ok=False;continue
  m1,r2=scan(*ch1,first=2)
  if m1!=em1 or r2 is None or((r2[2]if c2s==1 else r2[3])%c2q)!=0:live_ok=False;continue
  ch2=fch(c2s,c2q,r2[0],r2[1])
  if ch2 is None or scan(*ch2,first=0)!=em2:live_ok=False
 recs=[[int(z)for z in s.split()]for s in dc["compact_100"].split(";")]
 all_above_fb=all(r[1]>B_FB for r in recs)
 def _st(rs):
  ms=[r[3]for r in rs];n=len(ms);mn=sum(ms)/n;vr=sum((x-mn)**2 for x in ms)/(n-1)
  return dict(n=n,min=min(ms),max=max(ms),mean=round(mn,2),fano_v_over_m=round(vr/mn,3),n_zero=sum(1 for x in ms if x==0))
 sf_c,sf_d,sf_a=_st([r for r in recs if r[0]==1 and r[2]==0]),_st([r for r in recs if r[0]==1 and r[2]==1]),_st([r for r in recs if r[0]==1])
 sg_c,sg_d,sg_a=_st([r for r in recs if r[0]==0 and r[2]==0]),_st([r for r in recs if r[0]==0 and r[2]==1]),_st([r for r in recs if r[0]==0])
 RT=M.rho_tab_1024();ntab=len(RT);PLAN=D["descent"]["plan59"];CELLS={int(k):v for k,v in D["descent"]["cells"].items()};S_R4=log2(CELLS[4]);_wrc={}
 def ss2(m,nA,nB,kA,kB,al):
  if nB<=59.0:
   u=(m+al)/59.0
   if u<=1.0:return 1.0
   ii=ceil(u*1024.0);return RT[ii]if ii<ntab else 0.0
  key=(nA,nB,kA,kB)
  if key not in _wrc:
   h=(nB-59.0)/24.0;iA=int(round(24.0*(nA-59.0)/(nB-59.0)))
   wA=[log1p(h/(59.0+i*h))if i<iA else 0.0 for i in range(24)];wB=[log1p(h/(59.0+i*h))if i>=iA else 0.0 for i in range(24)]
   PA=[[1.0]]
   for j in range(1,kA+1):PA.append([x/j for x in _conv(PA[-1],wA)])
   _wrc[key]=(h,[(jA+jB,[(i,w)for i,w in enumerate(PA[jA]if jB==0 else _conv(PA[jA],wB))if w])for jA in range(kA+1)for jB in range(kB+1)])
  h,rows_w=_wrc[key];tot=0.0;inv_b=1024.0/59.0;hk=h*inv_b
  for j,nz in rows_w:
   lim=(m-j*59.0-j*h)/h if h>0 else 1e9;k0=(m-j*59.0+al)*inv_b
   for i,w in nz:
    if i<=lim:
     kk=ceil(k0-i*hk);tot+=w*(1.0 if kk<=0 else(RT[kk]if kk<ntab else 0.0))
  return min(1.0,tot)
 rc_l=[];rt_l=[]
 for side_s,q_s,_,_ in W["special_q"]:
  q=M.b62d(q_s);qb=log2(q);nl=_SQ_NORMS_CACHE.get((side_s,q))
  row=next((r for r in PLAN if r[2]==side_s and r[0]<=qb<=r[1]),None)
  if row is None or not nl:continue
  blo,bhi,side,n_i,iA,kA,kB,R,_=row;nA=59.0+iA*(n_i-59.0)/96.0;S_p=log2(CELLS[R]);dc_s=(S_p-S_R4)/12.0;cm=(bhi+S_p)/12.0
  pm=sum(w*ss2(48*cm+cf-(blo if side=="f" else 0),nA,n_i,kA,kB,0.0)*ss2(12*cm+LOG2U6+cg-(blo if side=="g" else 0),nA,n_i,kA,kB,3.4)for w,cf,cg in QBINS)
  pc=pt=0.0
  for lnf,lng in nl:
   nfe=lnf+48.0*dc_s-(qb if side=="f" else 0.0);nge=lng+12.0*dc_s-(qb if side=="g" else 0.0)
   pc+=ss2(nfe,nA,n_i,kA,kB,0.0)*ss2(nge,nA,n_i,kA,kB,3.4);pt+=ss2(nfe,nA,n_i,kA,kB,-19.93)*ss2(nge,nA,n_i,kA,kB,-0.04)
  rc_l.append(pc/(len(nl)*pm));rt_l.append(pt/(len(nl)*pm))
 ok=(tower_ok and live_ok and all_above_fb and len(recs)==100 and sf_a["n_zero"]==0 and sg_a["n_zero"]==0
 and max(sf_c["fano_v_over_m"],sf_d["fano_v_over_m"],sg_c["fano_v_over_m"],sg_d["fano_v_over_m"])<=1.5
 and len(rc_l)==24 and min(rc_l)>=1.5 and min(rt_l)>=5.0)
 return dict(ok=ok,tower_u5_p27631_r27481_deg8_frec_verified=tower_ok,all_100_above_fb_1024=all_above_fb,live_6_deg8_ideals_and_2hop_d2_exact=live_ok,
  nb_r=NB_R,nb_alpha=NB_ALPHA,b_fb=B_FB,deg8_f_ctrl=sf_c,deg8_f_d1_k1=sf_d,deg8_f_all=sf_a,deg2_g_ctrl=sg_c,deg2_g_d1_k1=sg_d,deg2_g_all=sg_a,
  f_d1_over_ctrl=round(sf_d["mean"]/sf_c["mean"],4),g_d1_over_ctrl=round(sg_d["mean"]/sg_c["mean"],4),
  prod_24_lattices_n=len(rc_l),prod_mu_cons_over_mod_min=round(min(rc_l),2),prod_mu_cons_over_mod_med=round(sorted(rc_l)[12],2),
  prod_mu_true_over_mod_min=round(min(rt_l),1),prod_mu_true_over_mod_med=round(sorted(rt_l)[12],1))
