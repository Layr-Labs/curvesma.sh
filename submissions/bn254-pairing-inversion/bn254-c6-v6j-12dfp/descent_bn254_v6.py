"""descent_bn254_v6.py — Pure-stdlib live 62-row LB=59 zero-retry-DAG + B_max=10 root-only i.i.d. boot-retry metered evaluator."""

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

NB_ALPHA=1.0

BOOT_RETRIES_MAX=10

BOOT_CAP_FAIL_ONE=1.0e-4

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

PAIR_SETUP=2.0 ** 25.1

LLL_NODE=2.0 ** 36

NORMS=2.0 ** 26

TDIV=2.0 ** 32

IDEAL=2.0 ** 34

PRIME=2.0 ** 24

GRID=4

SL_RANGE=(328,380)

def _log2li(l):

 x=l*LN2

 return l-log2(x)+log2(1.0+1.0/x+2.0/(x*x))

def _mulmod_instr(nbits):

 L=max(5,ceil(nbits/64.0))

 return MONT5*(L/5.0)** 2,L

def _gcd_instr(nbits):

 L=max(5,ceil(nbits/64.0))

 return 2.0*nbits*(12.0*L+16.0)

def _pi_upper(x):

 return 1.26*x/log(x)if x >= 17 else 7.0

def _stage_ops(b1,b2,wheel):

 s1=ceil(1.5*b1)*11

 phi=wheel

 for q in(2,3,5,7,11,13):

  if wheel%q == 0:

   phi=phi // q*(q-1)

 J=phi // 2

 G=ceil(b2/wheel)+1

 s2=6*(wheel // 2+G)+3*J*G if b2>b1 else 0

 return 1.1*(s1+s2),s1,s2,J,G

def _conv(a,b):

 na,nb=len(a),len(b)

 out=[0.0]*(na+nb-1)

 i0=0

 while i0<na and a[i0]== 0.0:i0 += 1

 if i0 == na:return out

 i1=na-1

 while i1>i0 and a[i1]== 0.0:i1 -= 1

 j0=0

 while j0<nb and b[j0]== 0.0:j0 += 1

 if j0 == nb:return out

 j1=nb-1

 while j1>j0 and b[j1]== 0.0:j1 -= 1

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

def evaluate_descent(desc_data,rho_fn,fft_F_fn,gbatch_cost_fn,pbits_fn,

                     log2_rels,peel_h,nf_main=405.0,ng_main=469.0):

 """Execute the COMPLETE 62-row LB=59 ZERO-RETRY DETERMINISTIC-DAG descent & boot LIVE!"""

 LB=59.0

 ORBIT=6

 CELLS={int(k):v for k,v in desc_data["cells"].items()}

 OMIT=desc_data["omitted_law"]

 BASE_PRIM=desc_data["base_prim_c6"]

 PLAN=desc_data["plan59"]

 BOOT=tuple(desc_data["boot"])

 EA=tuple(desc_data["ea"])

 RHO_TAB=[rho_fn(i*RSTEP)for i in range(32*1024+2)]

 def rho_lo(u):

  if u <= 1.0:return 1.0

  i=ceil(u/RSTEP)

  return RHO_TAB[i]if i<len(RHO_TAB)else 0.0

 def rho_hi(u):

  if u <= 1.0:return 1.0

  i=floor(u/RSTEP)

  if i+1 >= len(RHO_TAB):return 1.0/factorial(floor(u))

  return RHO_TAB[i]

 def _rho(rem,base):

  k=ceil(rem/(base*RSTEP))

  return 1.0 if k <= 0 else(RHO_TAB[k]if k<len(RHO_TAB)else 0.0)

 def kappa(u,t):

  ru=rho_fn(u)

  return rho_fn(u-t)/ru if ru>0.0 else 0.0

 @functools.lru_cache(None)

 def _ss_grid(b1,b2,steps=256):

  base,top=log2(b1),log2(b2)

  dl=(top-base)/steps

  return base,tuple((base+dl*i,log((base+dl*(i+1))/(base+dl*i)))for i in range(steps))if top>base else()

 @functools.lru_cache(None)

 def semismooth(fbits,b1,b2,steps=256):

  m=fbits-log2(TORSION)

  base,grid=_ss_grid(b1,b2,steps)

  acc=rho_lo(m/base)

  if grid:

   inv_bs=1.0/(base*RSTEP)

   rt=RHO_TAB;rl=len(rt);_ceil=ceil

   for lo,w in grid:

    u=m-lo

    if u <= base:acc += w

    else:

     idx=_ceil(u*inv_bs)

     if idx<rl:acc += rt[idx]*w

  return min(acc,1.0-1e-15)

 @functools.lru_cache(None)

 def power_corr(fbits,b1,b2):

  m=fbits-log2(TORSION)

  u=m/log2(b1)

  tilt=rho_hi(max(u-2.0,1.0))/max(rho_lo(u),1e-300)

  y=sqrt(b1)

  return max(0.0,1.0-tilt*2.0/(y*log(y)))

 @functools.lru_cache(None)

 def extraction_price(fbits,nbits,success=0.999):

  nb=ceil(nbits)

  mul,L=_mulmod_instr(nb)

  g=_gcd_instr(nb)

  F_ch=max(1,ceil(nb/20.0))

  best=None

  for k in range(10,37,2):

   b1=1 << k

   Q=ceil(_pi_upper(b1))

   queries=min(Q,2*F_ch*ceil(log2(Q))+2)

   for ratio in(1,10,30,100,300):

    b2=b1*ratio

    p=semismooth(fbits,b1,b2)*power_corr(fbits,b1,b2)

    if p <= 1e-12:continue

    curves=ceil(log1p(-success)/log1p(-SEPARATION*p))

    for wheel in(WHEELS if ratio>1 else(0,)):

     if wheel and wheel>b2:continue

     if wheel:ops,_,_,_,G=_stage_ops(b1,b2,wheel)

     else:s1=ceil(1.5*b1)*11;G=0;ops=1.1*s1

     cost=curves*(ops*mul+8*L*(Q+G)+2*g+5*(queries+G)*g)

     if best is None or cost<best:

      best=cost

  return best

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

  wB=[log1p(h/(base+i*h))if i >= iA else 0.0 for i in range(mesh)]

  PA=[[1.0]]

  for j in range(1,kmax+1):

   PA.append([x/j for x in _conv(PA[-1],wA)])

  PB=[[1.0],wB]

  rows={}

  for jA in range(kmax+1):

   for jB in range(2):

    rows[jA,jB]=PA[jA]if jB == 0 else _conv(PA[jA],PB[1])

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

     s += row[i]*(1.0 if k <= 0 else(RHO_TAB[k]if k<ntab else 0.0))

   out[jA][jB]=s

  cum=[[0.0,0.0]for _ in range(kmax+1)]

  for jB in range(2):

   acc=0.0

   for jA in range(kmax+1):

    acc += out[jA][jB]

    cum[jA][jB]=acc

  return tuple(tuple(x)for x in cum)

 def ss_lower2(m,base,n,iA,kA,kB,mesh,alpha):

  m=ceil(m)

  if n <= base:return rho_lo((m+alpha)/base)

  cum=ss_cum2(m,base,n,iA,mesh,alpha,kA)

  return cum[kA][0]+(cum[kA][1]if kB else 0.0)

 def ss_upper2(m,base,n,iA,kA,kB,mesh,alpha):

  m=floor(m)

  if n <= base:return rho_hi((m+alpha)/base)

  rows,h=window_rows(n,base,iA,mesh,kA)

  tot=0.0

  for(jA,jB),row in rows.items():

   if jA>kA or(jB and not kB):continue

   j=jA+jB

   for i,w in enumerate(row):

    if w and j*base+i*h <= m:

     u=(m-(j*base+i*h)-j*h+alpha)/base

     tot += w*(1.0 if u <= 1.0 else rho_hi(u))

  return min(tot,1.0)

 def surv_upper(m,base,T):

  if T <= 0:return 0.0

  if T >= m:return 1.0

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

 def miss_side(u4,u_main):

  u=u4/4.0

  tot=0.0

  for b10 in range(300,int(10*LB)):

   b=b10/10.0

   c=b+0.1

   lam=(2.0 **(log2_rels-c))*ORBIT*kappa(u_main,c/LB)

   tot += kappa(u,c/LB)*log(c/b)*exp(-lam*(1.0-(peel_h if b10 >= 544 else 0.4)))

  return tot

 u_main_f,u_main_g=nf_main/LB,ng_main/LB

 def missing(nf,ng):

  return(miss_side(ceil(4.0*nf/LB),u_main_f)if nf else 0.0)+(miss_side(ceil(4.0*ng/LB),u_main_g)if ng else 0.0)

 pairs=2.0*(2.0 ** _log2li(LB))

 cfmax=max(b[1]for b in QBINS)

 cgmax=max(b[2]for b in QBINS)

 W={};PHI={};CC={};DET={}

 def wm(s,lim):return max([v for(b,ss),v in W.items()if ss == s and b[0]<lim-1e-9]+[0.0])

 def pm(s,lim):return max([v for(b,ss),v in PHI.items()if ss == s and b[0]<lim-1e-9]+[0.0])

 def cm(s,lim):return max([v for(b,ss),v in CC.items()if ss == s and b[0]<lim-1e-9]+[0.0])

 for b_lo,b_hi,side,n_i,iA,kA,kB,R,_ in PLAN:

  assert n_i <= b_lo and(n_i<b_lo or(kA == 0 and kB == 0)),f"non-decreasing DAG edge at ({b_lo},{side})"

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

   nf=DF*c+cf-(b_lo if side == "f" else 0.0)

   ng=12.0*c+LOG2U6+cg-(b_lo if side == "g" else 0.0)

   pr += w*ss_lower2(nf,LB,n,iA,kA,kB,MESH,ALPHA_F)*ss_lower2(ng,LB,n,iA,kA,kB,MESH,ALPHA_G)

   nf_lo=DF*c_lo+cf_lo-(b_hi if side == "f" else 0.0)

   ng_lo=12.0*c_lo+LOG2U6+cg_lo-(b_hi if side == "g" else 0.0)

   sv += w*surv_upper(nf_lo,LB,T)*surv_upper(ng_lo,LB,T)

   sv0 += w*surv_upper(nf_lo,LB,T0)*surv_upper(ng_lo,LB,T0)

  nf_max=DF*c+cfmax-(b_lo if side == "f" else 0.0)

  ng_max=12.0*c+LOG2U6+cgmax-(b_lo if side == "g" else 0.0)

  det=1.0-(1.0-ECM_SUCCESS)*((nf_max+ng_max)/20.0+2.0*(kA+kB))

  miss=missing(nf_max,ng_max)

  mu=CELLS[R]*pr*UDUP*det*max(0.0,1.0-miss)

  ef=extraction_price(LB,nf_max,ECM_SUCCESS)

  eg=extraction_price(LB,ng_max,ECM_SUCCESS)

  sc_ecm=NORMS+TDIV+ef+eg+2.0*(extraction_price(nA,T,ECM_SUCCESS)+PRIME)+IDEAL

  sieve=CELLS[R]*(HITS_PER_CELL*SIEVE_HIT+SCAN+ZERO+BOUND_CELL)+pairs*PAIR_SETUP

  if CELLS[R]*sc_ecm<sieve+CELLS[R]*(sv+pr)*sc_ecm:

   how="ecm"

   mu *= COUNTED_FRAC

   Ns=CELLS[R]

  else:

   how="sieve"

   mu *= retention_node(nf_max,ng_max,n,iA,kA,kB,SLACK)

   Ns=CELLS[R]*(sv+pr)

  strip=NORMS+TDIV+strip_price(LB,nf_max,Ns)+strip_price(LB,ng_max,Ns)+IDEAL

  ext=2.0*(extraction_price(nA,T0,ECM_SUCCESS)+PRIME)+ef+eg

  gated=min(CELLS[R],CELLS[R]*(sv0+pr))

  own=(LLL_NODE+CELLS[R]*strip+gated*ext if how == "ecm" else LLL_NODE+sieve+Ns*strip+gated*ext)+2**40

  p0=(1.0+NB_ALPHA*mu/NB_R)**(-NB_R)if mu>0.0 else 1.0

  w_exp=own+kA*(wm("f",nA)+wm("g",nA))+kB*(wm("f",n)+wm("g",n))

  phi_dag=p0+kA*(pm("f",nA)+pm("g",nA))+kB*(pm("f",n)+pm("g",n))

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

  wh=[w[k]if kx_i <= k<iA else 0.0 for k in range(mesh)]

  wb=[w[k]if k >= iA else 0.0 for k in range(mesh)]

  PL=[[1.0]];PH=[[1.0]]

  for j in range(1,kA+1):

   PL.append([x/j for x in _conv(PL[-1],wl)])

   PH.append([x/j for x in _conv(PH[-1],wh)])

  F_tab={}

  for jj in range(kA+2):

   lim=floor((tau_i+1e-9-jj*base)/h-jj)

   if lim<0:continue

   H0=PH[jj][:lim+1]if jj <= kA else[]

   H1=_conv(PH[jj-1],wb)[:lim+1]if jj >= 1 else[]

   for jl in range(kA-jj+2):

    H=[]

    if H0 and jl <= kA-jj:H=list(H0)

    if H1 and jl <= kA-jj+1:

     if not H:H=list(H1)

     else:

      L=max(len(H),len(H1))

      H=[(H[i]if i<len(H)else 0.0)+(H1[i]if i<len(H1)else 0.0)for i in range(L)]

    if not H or not any(H):continue

    row=_conv(PL[jl],H)

    fc=F_tab.setdefault(jl+jj,[])

    if len(fc)<len(row):fc.extend([0.0]*(len(row)-len(fc)))

    for i,x in enumerate(row):fc[i]+= x

  return F_tab,h

 def ea_mass0(m):

  F_tab,h=ea_rows0(LB,float(n_b),iA_b,kA_b,MESH,kx,tau)

  tot=0.0

  for j,row in F_tab.items():

   L=min(len(row)-1,floor((m-j*h-j*LB)/h))+1

   if L <= 0:continue

   tot += sum(row[i]*_rho(m-j*LB-i*h+ALPHA_G,LB)for i in range(L))

  return tot

 def psurv(m,x,tau_v,mesh=192):

  h=(tau_v-x)/mesh

  w=[log1p(h/(x+k*h))for k in range(mesh)]

  tot=rho_hi(m/x)

  P=[1.0]

  for j in range(1,int(tau_v // x)+1):

   P=[y/j for y in _conv(P,w)]

   tot += sum(p*rho_hi(max(m-j*x-(i+j)*h,0.0)/x)

              for i,p in enumerate(P)if p and j*x+i*h <= tau_v+1e-9)

  return min(1.0,tot)

 pr_b=sum(w*ea_mass0(ceil(12.0*c_b+LOG2U6+cg))for w,cg in BOOT_BINS)

 det_b=1.0-0.001*(m_max/20.0+2.0*(kA_b+1))

 miss_b=missing(0.0,m_max)

 prob_b=RET_BOOT*pr_b*det_b*UDUP*COUNTED_FRAC*0.5*max(0.0,1.0-miss_b)

 surv_b=min(1.0,SURV_X*psurv(m_min-SURV_DM,x_ea,tau))

 sb_b=strip_price(x_ea,m_max,2.0 ** LB_WAVE)

 trial_b=2**40+NORMS+sb_b+surv_b*(extraction_price(x_ea,m_max,0.999)+extraction_price(nA_b,tau,0.999))+PRIME+IDEAL

 over_b=(2.0 ** LB_WAVE)*trial_b

 tree_b=2.0*ptree_cost(x_ea)

 T_boot_max=ceil(log(1.0/BOOT_CAP_FAIL_ONE)/prob_b)

 w_dag_one=kA_b*wm("g",nA_b)+wm("g",float(n_b))

 phi_dag_one=kA_b*pm("g",nA_b)+pm("g",float(n_b))

 c_det_one=kA_b*cm("g",nA_b)+cm("g",float(n_b))

 phi_run_one=min(0.999999,phi_dag_one+BOOT_CAP_FAIL_ONE)

 phi_log_one=phi_run_one ** BOOT_RETRIES_MAX

 e_runs_one=(1.0-phi_log_one)/(1.0-phi_run_one)

 one_run_exp=trial_b/prob_b+over_b+w_dag_one

 one_exp=tree_b+e_runs_one*one_run_exp

 two_exp_units=2.0*one_exp/UNIT

 two_log_u=log2(two_exp_units)

 markov_cap_fail=two_exp_units/(2.0 ** DESC_CAP_LOG2)

 two_dag_u=log2(2.0*e_runs_one*w_dag_one/UNIT)

 two_boot_u=log2(2.0*(tree_b+e_runs_one*(trial_b/prob_b+over_b))/UNIT)

 return dict(

  ok=(phi_run_one <= 0.35 and 2.0*phi_log_one<1e-4),

  mode="zero_retry_dag_k1_root_iid_b10_metered",

  nodes=len(DET),

  lb=LB,

  nb_r=NB_R,

  nb_alpha=NB_ALPHA,

  boot_retries_max=BOOT_RETRIES_MAX,

  base_prim=BASE_PRIM,

  counted_frac=COUNTED_FRAC,

  two_log=two_log_u,

  desc_cap_log2=DESC_CAP_LOG2,

  markov_cap_fail=markov_cap_fail,

  two_dag=two_dag_u,

  two_boot_exp=two_boot_u,

  boot_cfg=list(BOOT),

  boot_x_ea=round(x_ea,4),

  boot_root_bytes_log2=round(x_ea-2.0+log2(1.0/LN2),4),

  boot_prob=prob_b,

  boot_miss=miss_b,

  boot_T_max=T_boot_max,

  boot_cap_fail_one=BOOT_CAP_FAIL_ONE,

  phi_dag_one=phi_dag_one,

  phi_run_one=phi_run_one,

  phi_log_one=phi_log_one,

  phi_log_two=2.0*phi_log_one,

  e_runs_per_log=e_runs_one,

  c_det_one=c_det_one,

  c_det_two_max=2.0*BOOT_RETRIES_MAX*c_det_one,

  mu_min=min(d["mu"]for d in DET.values()),

  p0_max=max(d["p0"]for d in DET.values()),

  repeat_prob=0.0,

 )
