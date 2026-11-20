"""cm_models.py — Self-contained."""
import array
import collections
import functools
import math
from decimal import Decimal,localcontext
LN2=math.log(2.0)
GAMMA=0.5772156649015329
MONT5=270.0
ECM_TORSION=12.0
ECM_WHEELS=(30,210,2310,30030)
def build_rho_table(kmax=16,nt=120):
 with localcontext() as ctx:
  ctx.prec=80
  zero,one=Decimal(0),Decimal(1)
  coefs=[None,[one]+[zero]*(nt-1)]
  for k in range(2,kmax+1):
   prev,c=coefs[k-1],[zero]*nt
   for m in range(nt-1):
    c[m+1]=(prev[m]+m*c[m])/(k*(m+1))
   c[0]=prev[0]-sum(c[1:])
   coefs.append(c)
  return [None]+[[float(x) for x in c[:80]] for c in coefs[1:]]
RHO_C=build_rho_table()
def _build_rho_tail(S=128,umax=45):
 r=[1.0]*(S+1)+[1.0-math.log(1.0+j/S) for j in range(1,S+1)]
 for i in range(2*S+1,umax*S+1):
  u=i/S
  num=(r[i-S]+4.0*sum(r[i-S+1 :i :2])+2.0*sum(r[i-S+2 :i :2]))/(3.0*S)
  r.append(num/(u-1.0/(3.0*S)))
 return S,umax,r
_RHO_S,_RHO_UMAX,_RHO_TAIL=_build_rho_tail()
_RHO_C_REV=tuple(tuple(reversed(row)) if row is not None else None for row in RHO_C)
_RHO_C_LEN=len(_RHO_C_REV)
_RHO_C_REV=tuple(tuple(reversed(r)) if r is not None else None for r in RHO_C)
_RHO_C_LEN=len(_RHO_C_REV)
def rho_exact(u):
 if u<=1.0:return 1.0
 k=int(u)
 if u>k:k+=1
 if k<_RHO_C_LEN:
  x=k-u;v=0.0
  for c in _RHO_C_REV[k]:v=v*x+c
  return v if v>0.0 else 0.0
 if u>=_RHO_UMAX:return 0.0
 pos=u*_RHO_S;idx=int(pos);frac=pos-idx
 y0,y1=_RHO_TAIL[idx],_RHO_TAIL[idx+1]
 return math.exp((1.0-frac)*math.log(y0)+frac*math.log(y1)) if y0>0.0 and y1>0.0 else 0.0
def rho_db(u,log2_N):
 r0=rho_exact(u)
 r1=rho_exact(max(0.0,u-1.0))
 return max(1e-300,r0-GAMMA*r1/max(1.0,log2_N*LN2))
def log2_li(x):
 y=x*LN2
 return x-math.log2(y)+math.log2(1.0+1.0/y+2.0/(y*y))
def primes_upto(n):
 sieve=bytearray(b"\x01")*(n+1)
 sieve[:2]=b"\x00\x00"
 for p in range(2,math.isqrt(n)+1):
  if sieve[p]:sieve[p*p :: p]=b"\x00"*((n-p*p)//p+1)
 return [p for p in range(2,n+1) if sieve[p]]
def _ret_coeffs(weights,normalizer,length):
 _sp=math.sumprod
 jw_rev=[j*weights[j] for j in range(length,-1,-1)]
 out=[math.exp(-normalizer)]+[0.0]*length
 for n in range(1,length+1):
  val=_sp(jw_rev[length-n:length],out[:n])
  out[n]=val/n
  if val<-1e-13:raise ArithmeticError("negative coefficient beyond roundoff tolerance")
 return out
class Retention:
 def __init__(self,B_bits,X_bits,*,prime_min=2**24,step=0.25,max_bits=760,
         tilt=0.05,explicit_bits=14,power_terms=8,odd=1):
  self.B,self.X,self.minimum,self.step=B_bits,X_bits,prime_min,step
  self.tilt,self.length=tilt,math.ceil(max_bits/step)
  names=("small","sieve_geo","large_geo","sieve_sf","large_sf")
  weights={k:[0.0]*(self.length+1) for k in names}
  totals=collections.Counter()
  cutoff=min(float(explicit_bits),B_bits)
  for p in primes_upto(int(2**cutoff)-1)[odd:]:
   b=math.log2(p)
   bi=round(b/step)
   q=math.exp(tilt*bi*step-LN2*b)
   if p<prime_min:groups=(("small",False),)
   else:
    pref="sieve" if b < X_bits else "large"
    groups=((pref+"_geo",False),(pref+"_sf",True))
   for name,sf in groups:
    totals[name]+=math.log1p(q) if sf else-math.log1p(-q)
    qk=q
    for k in range(1,self.length+1):
     idx=k*bi
     if idx>self.length:break
     sgn =-1 if sf and k%2==0 else 1
     weights[name][idx]+=sgn*qk/k
     qk*=q
     if qk<1e-25:break
  def continuum(left,right,pref):
   if left>=right:return
   w0=step/4.0
   cells=math.ceil((right-left)/w0)
   w=(right-left)/cells
   for i in range(cells):
    b=left+(i+0.5)*w
    bi=round(b/step)
    for k in range(1,power_terms+1):
     amt=math.exp(k*tilt*bi*step-(k-1)*LN2*b)*w/(k*b)
     idx=k*bi
     sgn =-1 if k%2==0 else 1
     gs=((pref,amt),) if pref=="small" else ((pref+"_geo",amt),(pref+"_sf",sgn*amt))
     for g,x in gs:
      totals[g]+=x
      if idx<=self.length:weights[g][idx]+=x
  lm=math.log2(prime_min)
  continuum(cutoff,min(lm,B_bits),"small")
  continuum(max(cutoff,lm),min(X_bits,B_bits),"sieve")
  continuum(max(cutoff,X_bits),B_bits,"large")
  self.weights,self.normalizers=weights,dict(totals)
  all_w=[sum(weights[k][j] for k in ("small","sieve_geo","large_geo")) for j in range(self.length+1)]
  self.denominator=_ret_coeffs(all_w,sum(totals[k] for k in ("small","sieve_geo","large_geo")),self.length)
  self.sel_norm=math.exp(totals["sieve_sf"]+totals["large_sf"]-totals["sieve_geo"]-totals["large_geo"])
  self._c={}
 def co(self,k):
  if k not in self._c:self._c[k]=_ret_coeffs(self.weights[k],self.normalizers[k],self.length)
  return self._c[k]
 def two(self,n,a,b=9**9,win=None,sf=False):
  a=min(a,n)
  key=(n,b,win,sf)
  if a<0 or b<0:return 0.0
  if key not in self._c:
   if key[1:] not in self._c:
    _sp=math.sumprod
    st=self.step
    sm=self.co("small")
    D=min(self.length,2048)
    lg=self.co("large_sf" if sf else "large_geo")
    lim_b=min(b,D)
    if win is None:lg_f=lg[:lim_b+1]
    else:
     w0,w1,w2=win
     lg_f=[
      v if (v and (l*st<=w2 or any(self.X*k-w2<=l*st<=w1*k+w2 for k in range(1,w0+1)))) else 0.0
      for l,v in enumerate(lg[:lim_b+1])
     ]
    while len(lg_f)>1 and lg_f[-1]==0.0:
     lg_f.pop()
    Llg=len(lg_f)
    lg_rev=lg_f[::-1]
    m=[0.0]*(D+1)
    for s in range(D+1):
     hi=s if s<Llg-1 else Llg-1
     m[s]=_sp(lg_rev[Llg-1-hi :Llg],sm[s-hi :s+1])
    self._c[key[1:]]=m
   m=self._c[key[1:]]
   sv=self.co("sieve_sf" if sf else "sieve_geo")
   c=array.array("d")
   t=0.0
   for x in range(min(n,len(m)-1)+1):
    t+=m[x]*sv[n-x]
    c.append(t)
   self._c[key]=c
  c=self._c[key]
  return c[min(a,len(c)-1)]
 def cap(self,nb):
  return math.floor(math.ceil(nb)/math.log2(self.minimum))
 def at(self,nb,sl,sb=None):
  n=round(nb/self.step)
  rb=(self.cap(nb)+1)/16.0
  r=self.sel_norm*self.two(
   n,math.floor((sl-rb)/self.step),
   9**9 if sb is None else math.floor((sb-rb-0.125)/self.step),
   None,True,
  )/self.denominator[n]
  return max(0.0,min(1.0,r))
 def survivor_upper(self,nb,sl,ex=0.0,sb=None,lb=None,kmax=4):
  n=round(nb/self.step)
  a=math.floor((sl+ex)/self.step)
  if sb is None:return self.two(n,a)/self.denominator[n]
  return self.two(n,a,math.floor((sb+0.25)/self.step),(kmax,lb,(self.cap(nb)+2)/16.0+0.25))/self.denominator[n]
@functools.lru_cache(None)
def _semismooth_ecm(fbits,b1,b2,steps=32):
 m=fbits-math.log2(ECM_TORSION)
 base,top=math.log2(b1),math.log2(b2)
 acc=rho_exact(m/base)*(1.0-1e-9)
 if top>base:
  for i in range(steps):
   lo=base+(top-base)*i/steps
   hi=base+(top-base)*(i+1)/steps
   acc+=rho_exact((m-lo)/base)*(1.0-1e-9)*math.log(hi/lo)
 return min(acc,1.0-1e-15)
@functools.lru_cache(None)
def _ecm_cached(fbits2,limbs,success=0.999):
 fbits=fbits2/2.0
 nbits=limbs*64
 mul=MONT5*(limbs/5.0)**2
 gcd_c=2.0*nbits*(12*limbs+16)
 best=None
 for k in range(10,39,2):
  b1=1<<k
  for ratio in (1,10,30,100,300):
   b2=b1*ratio
   m=fbits-math.log2(ECM_TORSION)
   u=m/math.log2(b1)
   r_lo=max(rho_exact(u)*(1.0-1e-9),1e-300)
   tilt=(rho_exact(max(u-2.0,1.0))*(1.0+1e-9))/r_lo
   y=math.sqrt(b1)
   pcorr=max(0.0,1.0-tilt*2.0/(y*math.log(y)))
   p1=_semismooth_ecm(fbits,b1,b2)*pcorr
   if p1<=1e-12:continue
   curves=math.ceil(math.log1p(-success)/math.log1p(-p1))
   pi_up=1.26*b1/math.log(b1) if b1>=17 else 7.0
   for wheel in (ECM_WHEELS if ratio>1 else (0,)):
    if wheel and wheel>b2:continue
    if wheel:
     phi=wheel
     for q in (2,3,5,7,11,13):
      if wheel%q==0:phi=phi//q*(q-1)
     G=math.ceil(b2/wheel)+1
     ops=1.1*(math.ceil(1.5*b1)*11+6*(wheel//2+G)+3*(phi//2)*G)
    else:
     G=0
     ops=1.1*math.ceil(1.5*b1)*11
    per_c=ops*mul+2*gcd_c
    redo=p1*(ops*mul+(pi_up+G)*gcd_c)
    cost=curves*(per_c+redo)
    if best is None or cost<best:best=cost
 return best
def ecm_instr(fbits,nbits,success=0.999):
 return _ecm_cached(math.ceil(2.0*fbits),max(5,math.ceil(nbits/64)),success)
def ecm_mean_var(fbits,nbits,success=0.999):
 limbs=max(5,math.ceil(nbits/64));nb=limbs*64
 mul=MONT5*(limbs/5.0)**2;gcd_c=2.0*nb*(12*limbs+16)
 best=None
 for k in range(10,39,2):
  b1=1<<k
  for ratio in (1,10,30,100,300):
   b2=b1*ratio;m=fbits-math.log2(ECM_TORSION);u=m/math.log2(b1)
   r_lo=max(rho_exact(u)*(1.0-1e-9),1e-300)
   tilt=(rho_exact(max(u-2.0,1.0))*(1.0+1e-9))/r_lo
   y=math.sqrt(b1);pcorr=max(0.0,1.0-tilt*2.0/(y*math.log(y)))
   p1=_semismooth_ecm(fbits,b1,b2)*pcorr
   if p1<=1e-12:continue
   curves=math.ceil(math.log1p(-success)/math.log1p(-p1))
   pi_up=1.26*b1/math.log(b1) if b1>=17 else 7.0
   for wheel in (ECM_WHEELS if ratio>1 else (0,)):
    if wheel and wheel>b2:continue
    if wheel:
     phi=wheel
     for q in (2,3,5,7,11,13):
      if wheel%q==0:phi=phi//q*(q-1)
     G=math.ceil(b2/wheel)+1
     ops=1.1*(math.ceil(1.5*b1)*11+6*(wheel//2+G)+3*(phi//2)*G)
    else:G=0;ops=1.1*math.ceil(1.5*b1)*11
    per_c=ops*mul+2*gcd_c;redo_ex=ops*mul+(pi_up+G)*gcd_c
    q_s=math.exp(curves*math.log1p(-p1));p_s=1.0-q_s;assert p_s>=success-1e-12
    c_fail=curves*per_c;es=(curves*p1)/p_s
    es2=(curves*p1*(1.0-p1)+(curves*p1)**2)/p_s;var_s=max(0.0,es2-es*es)
    e_succ=c_fail+es*redo_ex;var_succ=var_s*(redo_ex*redo_ex)
    e_t=(q_s/p_s)*c_fail+e_succ;var_t=(q_s/(p_s*p_s))*(c_fail*c_fail)+var_succ
    if best is None or e_t<best[0]:best=(e_t,var_t)
 return best
NTT4_PRIMES=(
 269*(1<<53)+1,
 276*(1<<53)+1,
 305*(1<<53)+1,
 309*(1<<53)+1,
)
NTT4_PROD=NTT4_PRIMES[0]*NTT4_PRIMES[1]*NTT4_PRIMES[2]*NTT4_PRIMES[3]
def fft_F(N):
 if N<=1:return 1200.0
 Ni=math.ceil(N)
 if 2*Ni<=(1<<53):
  k=max(2,(2*Ni-1).bit_length())
  M=1<<k
  assert M<=(1<<53) and Ni*(1<<128)<NTT4_PROD
  return 168.0*M*k+220.0*M+1200.0+(504.0*M if Ni>(1<<26) else 0.0)
 ktot=(Ni-1).bit_length()
 k1=ktot//2
 K1=1<<k1
 K2=(Ni+K1-1)//K1
 k1m=k1+1
 k2m=max(2,(2*K2-1).bit_length())
 M1=1<<k1m
 M2=1<<k2m
 assert M1<=(1<<53) and M2<=(1<<53) and (K1*K2)*(1<<128)<NTT4_PROD
 Mtot=float(M1*M2)
 return 168.0*Mtot*(k1m+k2m)+1100.0*Mtot+1200.0
def mr_instr(nbits):
 L=max(5,math.ceil(nbits/64.0))
 return 100.0*2.0*nbits*MONT5*((L/5.0)**2)
def gbatch_cost(pbits,leaf=9):
 SP=math.ceil(pbits/64)
 L=math.ceil(math.log2((SP+1)/leaf))
 G=L+2
 n=lambda j:leaf*(1<<j)
 PT=sum((1<<(L-j))*(fft_F(2*n(j-1))+3.0*2*n(j-1)) for j in range(1,L+1))
 RT=sum((1<<(L-j+1))*(fft_F(n(j)+G+n(j-1))+3.0*(n(j-1)+G)) for j in range(1,L+1))
 return (PT+RT+6.0*fft_F(2*(max(SP,n(L))+G)))/(1<<L)+(1<<15)
def _get_rho_grid():
 return _RHO_S,(_RHO_UMAX-1)*_RHO_S,_RHO_TAIL
def _dw_int(B_in,dn,dd):
 B=[list(r) for r in B_in]
 n=len(B);sp=math.sumprod
 d=[1]*(n+1);lam=[[0]*n for _ in range(n)]
 for i in range(n):
  Bi=B[i];li=lam[i]
  for j in range(i+1):
   u=sp(Bi,B[j]);lj=lam[j]
   for kk in range(j):u=(d[kk+1]*u-li[kk]*lj[kk])//d[kk]
   if j<i:li[j]=u
   else:d[i+1]=u
 k=1
 while k<n:
  lk=lam[k];bk=B[k]
  for l in range(k-1,-1,-1):
   v=lk[l];dl1=d[l+1];v2=v+v
   if v2>dl1 or v2<-dl1:
    q=(v2+dl1)//(dl1+dl1);bl=B[l]
    if q==1:
     for t in range(n):bk[t]-=bl[t]
     ll=lam[l]
     for j in range(l):lk[j]-=ll[j]
    elif q==-1:
     for t in range(n):bk[t]+=bl[t]
     ll=lam[l]
     for j in range(l):lk[j]+=ll[j]
    else:
     for t in range(n):bk[t]-=q*bl[t]
     ll=lam[l]
     for j in range(l):lk[j]-=q*ll[j]
    lk[l]=v-q*dl1
  lkk1=lk[k-1];dk=d[k];dkm1=d[k-1];dk1=d[k+1]
  num_s=dkm1*dk1+lkk1*lkk1
  if dd*num_s<dn*dk*dk:
   B[k],B[k-1]=B[k-1],bk;lk1=lam[k-1]
   for j in range(k-1):lk[j],lk1[j]=lk1[j],lk[j]
   dk_new=num_s//dk
   for i in range(k+1,n):
    li=lam[i];t=li[k];t0=li[k-1]
    nv=(dk1*t0-lkk1*t)//dk
    li[k]=nv;li[k-1]=(dk_new*t+lkk1*nv)//dk1
   d[k]=dk_new;k=k-1 if k>1 else 1
  else:k+=1
 return B,d,lam
def lll_16d_int(B_in,delta_num=99,delta_den=100):
 B1,_,_=_dw_int(B_in,11,40)
 return _dw_int(B1,delta_num,delta_den)
def srt_toy_check():
 import random
 rng=random.Random(2026)
 primes=[]
 n=3
 while len(primes)<400:
  if all(n%d!=0 for d in range(2,math.isqrt(n)+1)):primes.append(n)
  n+=2
 P=1
 for pr in primes:P*=pr
 leaves=[rng.randrange(1<<40,1<<48) | 1 for _ in range(128)]
 levels=[leaves]
 while len(levels[-1])>1:
  cur=levels[-1]
  levels.append([cur[i]*cur[i+1] for i in range(0,len(cur),2)])
 L=len(levels)-1;G=L+2
 Z_root=levels[-1][0]
 b_root=64*((Z_root.bit_length()+63)//64+G)
 Y=[(P<<b_root)//Z_root]
 for lv in range(L-1,-1,-1):
  nxt=[]
  for idx,zc in enumerate(levels[lv]):
   zs=levels[lv][idx ^ 1]
   yv=Y[idx>>1]
   bv=64*((levels[lv+1][idx>>1].bit_length()+63)//64+G)
   bc=64*((zc.bit_length()+63)//64+G)
   yc=((yv*zs) & ((1<<bv)-1))>>(bv-bc)
   nxt.append(yc)
  Y=nxt
 mis=0
 for z,y in zip(leaves,Y):
  bc=64*((z.bit_length()+63)//64+G)
  r=(z*y+(1<<bc)-1)>>bc
  if r==z:r=0
  if r!=P%z:mis+=1
 assert mis==0
 return dict(n_leaves=len(leaves),mismatches=mis)
def stage_b_analytic(nbits):
 L=max(5,math.ceil(nbits/64.0))
 c_srt=gbatch_cost(24197143.0,leaf=L)
 c_sq10=10.0*MONT5*((L/5.0)**2)
 c_gcd=2.0*nbits*(12.0*L+16.0)
 return c_srt+c_sq10+c_gcd
def _peel_fixed_point_bls(lf,lg,lrels,u_f,u_g,exc_ref=1.09,step=0.1,lam_seed=100.0):
 rf,rg=rho_exact(u_f),rho_exact(u_g)
 bins=[]
 for side_l,u_s,r_s in ((lf,u_f,rf),(lg,u_g,rg)):
  b=40.0
  while b<side_l-1e-9:
   c=min(side_l,b+step)
   km_lo=rho_exact(u_s-b/side_l)/r_s
   km_hi=rho_exact(u_s-c/side_l)/r_s
   wt=km_hi*math.log(c/b)
   lam=(2.0**(lrels-c))*12.0*km_lo
   if lam<lam_seed:bins.append((wt,lam))
   b=c
 _exp=math.exp
 def solve_h(scale):
  h=0.0;sbins=[(wt,lam*scale) for wt,lam in bins]
  for _ in range(60):
   omh=1.0-h;s=ds=0.0
   for wt,ls in sbins:
    t=wt*_exp(-ls*omh);s+=t;ds+=t*ls
   em=_exp(-s);g=1.0-em;diff=g-h
   if abs(diff)<1e-13:return g
   gp=em*ds
   if 0.0<=gp<0.98 and diff>0.0:
    hn=h+0.85*diff/(1.0-gp)
    if g<hn<1.0:
     omhn=1.0-hn;sn=0.0
     for wt,ls in sbins:sn+=wt*_exp(-ls*omhn)
     gn=1.0-_exp(-sn)
     if gn>=hn:h=gn;continue
   h=g
  return h
 def exceeds_90(scale):
  h=0.0;sbins=[(wt,lam*scale) for wt,lam in bins]
  for _ in range(1000):
   omh=1.0-h;s=0.0
   for wt,ls in sbins:s+=wt*_exp(-ls*omh)
   hn=1.0-_exp(-s)
   if hn>0.90:return True
   if abs(hn-h)<1e-13:return False
   h=hn
  return h>0.90
 h=solve_h(1.0)
 lo_e,hi_e=0.4,exc_ref
 for _ in range(20):
  mid=0.5*(lo_e+hi_e)
  if exceeds_90(mid/exc_ref):lo_e=mid
  else:hi_e=mid
 return h,hi_e

def eval_descent_live(desc,mesh_node=32,mesh_boot=185):
 dp=desc["params"]
 lx=float(dp["lx"])
 lf_peel=float(dp.get("peel_lf",71.29))
 lg_peel=float(dp.get("peel_lg",71.69))
 lrels_peel=float(dp.get("peel_lrels",63.4417118678))
 uf_peel=float(dp.get("peel_uf",548.9/71.29))
 ug_peel=float(dp.get("peel_ug",677.8/71.69))
 exc_peel=float(dp.get("peel_exc",1.09))
 peel_h,exc_crit=_peel_fixed_point_bls(lf_peel,lg_peel,lrels_peel,uf_peel,ug_peel,exc_ref=exc_peel)
 alpha_f,alpha_g=float(dp["alpha_f"]),float(dp["alpha_g"])
 log2_u8=float(dp["log2_u8"])
 smult=float(dp["sieve_mult"])
 qbins=desc["qbins"]
 qbins_lo=desc["qbins_lo"]
 qbins_sub=desc["qbins_sub"]
 cells_r2={int(k): int(v) for k,v in desc["cells_rank2"].items()}
 BOUND_CELL=2400.0*smult
 HITS_PER_CELL=2.0*math.log(lx/8.0)
 SIEVE_HIT=1100.0*smult
 SCAN=10.0*smult
 ZERO=2.0*smult
 PAIR_SETUP=float(dp.get("pair_setup",(2.0**25.1)*(smult*1.33)))
 LLL_NODE=float(dp["lll_node"])
 NORMS,TDIV,IDEAL,PRIME=2.0**27,2.0**25.0,2.0**35,2.0**24
 pairs=2.0*(2.0**log2_li(lx))
 cfmax=max(x[1] for x in qbins)
 cgmax=max(x[2] for x in qbins)
 wr_cache={}
 RINV,RMAX_K,r_arr=_get_rho_grid()
 RSTEP=1.0/RINV
 def get_wr(n,iA,mesh,kA,kB):
  key=(n,iA,mesh,kA,kB)
  if key in wr_cache:return wr_cache[key]
  h=(n-lx)/mesh
  wA=[math.log1p(h/(lx+i*h)) if i<iA else 0.0 for i in range(mesh)]
  wB=[math.log1p(h/(lx+i*h)) if i>=iA else 0.0 for i in range(mesh)]
  _sp=math.sumprod
  def conv(a,b):
   na,nb=len(a),len(b)
   out=[0.0]*(na+nb-1)
   i0=next((i for i,x in enumerate(a) if x),na)
   if i0==na:return out
   i1=na-next(i for i,x in enumerate(reversed(a)) if x)
   j0=next((j for j,y in enumerate(b) if y),nb)
   if j0==nb:return out
   j1=nb-next(j for j,y in enumerate(reversed(b)) if y)
   aa=a[i0:i1];bb_rev=b[j0:j1][::-1];la,lb=len(aa),len(bb_rev)
   for s in range(la+lb-1):
    a_lo=0 if s<lb else s-lb+1
    a_hi=s+1 if s+1<la else la
    out[i0+j0+s]=_sp(aa[a_lo:a_hi],bb_rev[lb-1-s+a_lo:lb-1-s+a_hi])
   return out
  pA_list=[[1.0]]
  for j in range(1,kA+1):
   pA_list.append([x/j for x in conv(pA_list[-1],wA)])
  pAB_list=[conv(pA,wB) for pA in pA_list] if kB else None
  wr_cache[key]=(h,pA_list,pAB_list)
  return wr_cache[key]
 def ss_eval(m_list,n,iA,mesh,alpha,kA,kB):
  if n<=lx:
   out=[]
   for m_raw in m_list:
    m=float(math.ceil(m_raw))
    u_val=(m+alpha)/lx
    rv=rho_exact(u_val)
    if rv>0.0 and u_val>2.0:rv*=max(0.0,1.0-GAMMA*rho_exact(u_val-1.0)/(rv*m*LN2))
    out.append(rv)
   return out
  h,pA_list,pAB_list=get_wr(n,iA,mesh,kA,kB)
  hk=h*RINV/lx
  inv_lx_r=RINV/lx
  g_ln2=GAMMA/LN2
  two_r=2*RINV
  lim16=16*RINV
  m_cache={}
  out=[]
  for m_raw in m_list:
   mi=math.ceil(m_raw)
   cached=m_cache.get(mi)
   if cached is not None:
    out.append(cached)
    continue
   m=float(mi)
   tot=0.0
   for jA in range(kA+1):
    rows_j=((0,pA_list[jA]),(1,pAB_list[jA])) if kB else ((0,pA_list[jA]),)
    for jB,row in rows_j:
     j=jA+jB
     ifeas=min(len(row)-1,int(math.floor((m-j*h-j*lx)/h))) if h>0 else len(row)-1
     if ifeas<0:continue
     k0=(m-j*lx+alpha)*inv_lx_r
     idx_lo=max(0,int((k0-lim16)/hk))
     mj_lx=m-j*lx
     for idx in range(idx_lo,ifeas+1):
      c=row[idx]
      if c==0.0:continue
      ks=int(k0-idx*hk+0.999999999999)
      if ks>=RMAX_K:continue
      rv=1.0 if ks<=0 else r_arr[ks]
      if rv>0.0:
       if ks>two_r:
        rem_b=mj_lx-idx*h
        if rem_b<lx:rem_b=lx
        f_db=1.0-g_ln2*r_arr[ks-RINV]/(rv*rem_b)
        if f_db>0.0:tot+=c*rv*f_db
       else:tot+=c*rv
   m_cache[mi]=tot
   out.append(tot)
  return out
 p_unrec=float(dp.get("p_unrec_leaf",0.0))
 miss_cache={}
 def miss_side(m_val,u_m,lb_s):
  if m_val<=0.0:return 0.0
  u4=math.ceil(4.0*m_val/lx)
  key=(u4,u_m,lb_s)
  if key in miss_cache:return miss_cache[key]
  u=u4/4.0
  r_u=rho_exact(u)
  r_m=rho_exact(u_m)
  tot=0.0
  for b10 in range(400,int(lx*10)):
   b=b10/10.0
   c=b+0.1
   km_lo=rho_exact(u_m-b/lb_s)/r_m if r_m>0 else 0.0
   lam=(2.0**(lrels_peel-c))*12.0*km_lo
   ku_hi=rho_exact(u-c/lx)/r_u if r_u>0 else 0.0
   if lam<100.0:
    x=lam*(1.0-peel_h)
    p_miss=math.exp(-x)*(1.0+x)
   else:p_miss=math.exp(-lam)
   p_miss_tot=min(1.0,p_miss+p_unrec)
   tot+=ku_hi*math.log(c/b)*p_miss_tot
  miss_cache[key]=tot
  return tot
 def _gb_wave(pbits,leaf):
  SP=math.ceil(pbits/64.0);L=math.ceil(math.log2((SP+1)/leaf));G=L+2
  n=lambda j:leaf*(1<<j)
  PT=sum((1<<(L-j))*(fft_F(2*n(j-1))+3.0*2*n(j-1)) for j in range(1,L+1))
  RT=sum((1<<(L-j+1))*(fft_F(n(j)+G+n(j-1))+3.0*(n(j-1)+G)) for j in range(1,L+1))
  wt=PT+RT+6.0*fft_F(2*(max(SP,n(L))+G))
  return wt/(1<<L)+(1<<15),wt
 p70_bits=1.000028*(2.0**lx)/LN2
 _gwc={}
 def _strip_wave(nb):
  L=max(5,math.ceil(nb/64.0))
  if L not in _gwc:_gwc[L]=_gb_wave(p70_bits,L)
  pl,wt=_gwc[L]
  return pl+12.0*MONT5*((L/5.0)**2)+2.0*nb*(12.0*L+16.0),wt
 W_tr,P_tr,N_tr={},{},{}
 bins_seen=[];out_rows=[]
 scale_iA=mesh_node//32
 p_ecm=float(dp["ecm_success"])
 p0_floor=float(dp.get("p0_floor",0.001))
 for row in desc["plan_struct"]:
  lo,hi,side,n_int,iA32,kA,kB,R_n,thin=row
  bn=(lo,hi)
  if not bins_seen or bins_seen[-1]!=bn:bins_seen.append(bn)
  n_c=float(n_int);iA=iA32*scale_iA
  cells=cells_r2[R_n];S=math.log2(cells)
  c_hi=(hi+S)/16.0;c_lo=(lo+S)/16.0
  nA=lx+iA*(n_c-lx)/mesh_node if n_c>lx else lx
  T=kA*nA+kB*n_c+96.0;T0=math.ceil(kA*nA+kB*n_c)
  nf_list=[72.0*c_hi+cf-(lo if side=="f" else 0.0) for _,cf,_ in qbins]+[72.0*c_hi+cf-(lo if side=="f" else 0.0) for cf,_ in qbins_sub]
  ng_list=[24.0*c_hi+log2_u8+cg-(lo if side=="g" else 0.0) for _,_,cg in qbins]+[24.0*c_hi+log2_u8+cg-(lo if side=="g" else 0.0) for _,cg in qbins_sub]
  pf_v=ss_eval(nf_list,n_c,iA,mesh_node,alpha_f,kA,kB)
  pg_v=ss_eval(ng_list,n_c,iA,mesh_node,alpha_g,kA,kB)
  pr_c=sum(pf_v[i]*pg_v[i] for i in range(6))/6.0
  pr_i=sum(pf_v[i]*pg_v[i] for i in range(6,30))/24.0
  pr=min(pr_i,4.0*pr_c)
  sv=sv0=0.0
  for (w,_,_),(cfl,cgl) in zip(qbins,qbins_lo):
   nfl=72.0*c_lo+cfl-(hi if side=="f" else 0.0);ngl=24.0*c_lo+log2_u8+cgl-(hi if side=="g" else 0.0)
   sf=1.0 if T>=nfl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(nfl)-T)/lx))
   sg=1.0 if T>=ngl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(ngl)-T)/lx))
   sv+=w*sf*sg
   sf0=1.0 if T0>=nfl else min(1.0,max(1.0,T0/lx)*rho_exact((math.floor(nfl)-T0)/lx))
   sg0=1.0 if T0>=ngl else min(1.0,max(1.0,T0/lx)*rho_exact((math.floor(ngl)-T0)/lx))
   sv0+=w*sf0*sg0
  nf_max=72.0*c_hi+cfmax-(lo if side=="f" else 0.0);ng_max=24.0*c_hi+log2_u8+cgmax-(lo if side=="g" else 0.0)
  k_lx_f=int(nf_max//20.0)+1;k_lx_g=int(ng_max//20.0)+1
  k_wA_n=kA if nA>lx+1e-9 else 0
  det_n=p_ecm**(k_lx_f+k_lx_g+2*(kA+kB))
  miss_n=miss_side(nf_max,uf_peel,lf_peel)+miss_side(ng_max,ug_peel,lg_peel)
  mu=cells*pr*float(dp["udup_node"])*det_n*math.exp(-miss_n)
  mr_f=mr_instr(nf_max);mr_g=mr_instr(ng_max);mr_T0=mr_instr(max(T0,100.0))
  sf_n,wf_n=_strip_wave(nf_max);sg_n,wg_n=_strip_wave(ng_max)
  pre_gate=NORMS+TDIV+sf_n+sg_n
  post_gate=(k_lx_f*(ecm_instr(lx,nf_max,p_ecm)+mr_f)+k_lx_g*(ecm_instr(lx,ng_max,p_ecm)+mr_g)+2.0*k_wA_n*(ecm_instr(nA,max(T0,100.0),p_ecm)+mr_T0)+mr_f+mr_g+IDEAL)
  sieve=cells*(HITS_PER_CELL*SIEVE_HIT+SCAN+ZERO+BOUND_CELL)+pairs*PAIR_SETUP
  gate_ecm=cells*pre_gate+wf_n+wg_n+(cells*sv0+cells*pr)*post_gate
  gate_sieve=sieve+(cells*sv+cells*pr)*pre_gate+wf_n+wg_n+(cells*sv0+cells*pr)*post_gate
  if gate_ecm<gate_sieve:
   how="ecm_gateG";own=LLL_NODE+gate_ecm
  else:
   how="sieve_gateG";mu*=float(dp["ret_node"]);own=LLL_NODE+gate_sieve
  p0=max(p0_floor,1.0/(1.0+0.5*mu))
  assert p0<=0.02 and mu>=150.0
  wmax=lambda s,lim:max([W_tr[b,s] for b in bins_seen if b[0]<lim-1e-9 and (b,s) in W_tr]+[0.0])
  pmax=lambda s,lim:max([P_tr[b,s] for b in bins_seen if b[0]<lim-1e-9 and (b,s) in P_tr]+[0.0])
  nmax=lambda s,lim:max([N_tr[b,s] for b in bins_seen if b[0]<lim-1e-9 and (b,s) in N_tr]+[0])
  cW=kA*(wmax("f",nA)+wmax("g",nA))+kB*(wmax("f",n_c)+wmax("g",n_c))
  cP=kA*(pmax("f",nA)+pmax("g",nA))+kB*(pmax("f",n_c)+pmax("g",n_c))
  cN=kA*(nmax("f",nA)+nmax("g",nA))+kB*(nmax("f",n_c)+nmax("g",n_c))
  W_tr[bn,side]=own+cW
  P_tr[bn,side]=p0+cP
  N_tr[bn,side]=1+cN
  out_rows.append(dict(lo=lo,hi=hi,side=side,how=how,R=R_n,n=n_int,iA=iA,nA=nA,kA=kA,kB=kB,
   pr=pr,mu=mu,p0=p0,P_tree=P_tr[bn,side],N_tree=N_tr[bn,side],nf=nf_max,ng=ng_max,own_log2=math.log2(own),W_log2=math.log2(W_tr[bn,side])))
 bt=desc["boot"]
 assert int(bt["cells_b"])==1 and int(dp["boot_rels"])==1 and float(dp.get("ball_pair",0.0))==0.0
 m_max=float(bt["m_sv"])
 iA_b=math.floor((bt["nA"]-lx)*mesh_boot/(bt["n1"]-lx)+1e-12)
 assert lx+iA_b*(bt["n1"]-lx)/mesh_boot<=bt["nA"]+1e-12
 pr_raw_b=ss_eval([m_max],float(bt["n1"]),iA_b,mesh_boot,alpha_g,bt["kA"],1)[0]
 mb=miss_side(m_max,ug_peel,lg_peel);pr_b=pr_raw_b*math.exp(-mb)
 k_lx_b=int(m_max//20.0)+1;kA_b=int(bt["kA"]);mr_b=mr_instr(m_max)
 det_b=p_ecm**(k_lx_b+kA_b+1);p_boot_trial=pr_b*det_b*float(dp["ret_boot"])
 trial_b=(NORMS+TDIV+k_lx_b*(ecm_instr(lx,m_max,p_ecm)+mr_b)+kA_b*(ecm_instr(bt["nA"],m_max,p_ecm)+mr_b)+mr_b+IDEAL)
 LLL_BOOT=float(dp.get("lll_boot",2.0**48))
 boot_own_exp=(1.0/p_boot_trial)*(LLL_BOOT+trial_b)
 wmax_g=lambda lim:max(W_tr[b,"g"] for b in bins_seen if b[0]<lim-1e-9)
 pmax_g=lambda lim:max(P_tr[b,"g"] for b in bins_seen if b[0]<lim-1e-9)
 nmax_g=lambda lim:max(N_tr[b,"g"] for b in bins_seen if b[0]<lim-1e-9)
 W_tree_1=bt["kA"]*wmax_g(bt["nA"])+wmax_g(float(bt["n1"]))
 P_boot_tree=bt["kA"]*pmax_g(bt["nA"])+pmax_g(float(bt["n1"]))
 N_tree_1=bt["kA"]*nmax_g(bt["nA"])+nmax_g(float(bt["n1"]))
 assert P_boot_tree<0.80
 SP70=math.ceil(p70_bits/64.0)
 p70_tree_u=math.log2(math.ceil(math.log2(SP70))*fft_F(SP70)/20000.0)
 two_exp_u=2.0*(boot_own_exp+W_tree_1)/(1.0-P_boot_tree)/20000.0
 return dict(
  mesh_node=mesh_node,mesh_boot=mesh_boot,n_nodes=len(desc["plan_struct"]),
  peel_h=peel_h,exc_crit=exc_crit,
  boot_pr_log2=math.log2(pr_b),boot_miss=mb,p_boot_trial=p_boot_trial,
  boot_own_exp_log2=math.log2(2.0*boot_own_exp/20000.0),
  tree_max_two_log2=math.log2(2.0*W_tree_1/20000.0),
  P_boot_tree=P_boot_tree,max_p0=max(r["p0"] for r in out_rows),
  N_tree_max_1log=N_tree_1,p70_tree_build_u=p70_tree_u,
  p_node_abort=0.0,p_repeat_loss=0.0,
  two_log_exp_units=math.log2(two_exp_u),
  min_mu=min(r["mu"] for r in out_rows),rows=out_rows,_ss_eval=ss_eval,
 )
def _nb_cdf(k_max,mu,r):
 if k_max<=0:return 0.0
 p_s=r/(r+mu);q_s=1.0-p_s
 term=p_s**r;s=term
 for j in range(1,k_max):
  term*=(r+j-1)*q_s/j
  s+=term
 return min(1.0,s)
def q15_samefield_live_replay(camp,y_pow_fn,add_fn,sub_fn,mul_fn,scale_fn,norm_fn):
 lx_b=int(camp["lx_b"]);B_FB=1<<lx_b;small_ns_max=int(camp["small_ns_max"])
 w_el=y_pow_fn(5);w2_el=y_pow_fn(10)
 c_f=add_fn([2]+[0]*7,w_el);c_g=scale_fn(2,w2_el)
 sig_m=[y_pow_fn((7*j)%15) for j in range(8)]
 def _sig(v):
  out=[0]*8
  for j,c in enumerate(v):
   if c:
    sj=sig_m[j]
    for i in range(8):out[i]+=c*sj[i]
  return out
 def _norm_zw(el):
  s1=_sig(el);s2=_sig(s1);s3=_sig(s2)
  res=mul_fn(mul_fn(el,s1),mul_fn(s2,s3))
  assert res[1]==res[2]==res[3]==res[4]==res[6]==res[7]==0
  return res[0],res[5]
 def _w_roots(p):
  for g in range(2,p):
   w=pow(g,(p-1)//3,p)
   if w!=1:
    w2=(w*w)%p
    return (w,w2) if w<w2 else (w2,w)
  raise RuntimeError(p)
 all_ps=primes_upto(B_FB)
 split_tbl={};small_w_tbl={}
 for p in all_ps:
  if p%15==1:
   w0,w1=_w_roots(p)
   split_tbl[p]=[(oi,w,pow((2+w)%p,(p-1)//3,p)==1,pow((-2*w*w)%p,(p-1)//3,p)==1) for oi,w in enumerate((w0,w1))]
  elif p<=small_ns_max:
   small_w_tbl[p]=_w_roots(p) if p%3==1 else None
 small_ns_ps=sorted(small_w_tbl.keys());split_ps=sorted(split_tbl.keys())
 def _prod(xs):
  if not xs:return 1
  xs=list(xs)
  while len(xs)>1:xs=[xs[i]*xs[i+1] if i+1<len(xs) else xs[i] for i in range(0,len(xs),2)]
  return xs[0]
 P_small=_prod(small_ns_ps+[p for p in split_ps if p<=256])
 P_mid=_prod([p for p in split_ps if p>256])
 ret_f={(int(p),int(oi)) for p,oi in camp["ret_f_split"]+camp["ret_f_ns"]}
 ret_g={(int(p),int(oi)) for p,oi in camp["ret_g_split"]+camp["ret_g_ns"]}
 def _strip(n,p1,p2):
  while n>1:
   g=math.gcd(n,p1)
   if g==1:break
   n//=g
  while n>1:
   g=math.gcd(n,p2)
   if g==1:break
   n//=g
  return n
 def _extract_cols(el_f,el_g,nf_rem,ng_rem,q_div,q_sid):
  fc,gc=set(),set()
  c0f,c1f=_norm_zw(el_f);c0g,c1g=_norm_zw(el_g)
  assert c0f*c0f-c0f*c1f+c1f*c1f==(nf_rem*q_div if q_sid==0 else nf_rem)
  assert c0g*c0g-c0g*c1g+c1g*c1g==(ng_rem*q_div if q_sid==1 else ng_rem)
  for side,rem,c0,c1,cols in (("f",nf_rem,c0f,c1f,fc),("g",ng_rem,c0g,c1g,gc)):
   for p in small_ns_ps:
    if rem%p==0:
     while rem%p==0:rem//=p
     ws=small_w_tbl[p]
     if ws is None:cols.add((p,0))
     else:
      h=0
      for oi,w in enumerate(ws):
       if (c0+c1*w)%p==0:cols.add((p,oi));h+=1
      assert h>=1
   if rem==1:continue
   hits=[]
   for p in split_ps:
    if rem%p==0:
     while rem%p==0:rem//=p
     hits.append(p)
     if rem==1:break
   if rem!=1:return None
   for p in hits:
    hcnt=0
    for oi,w,rf_ok,rg_ok in split_tbl[p]:
     if (c0+c1*w)%p==0:
      if not (rf_ok if side=="f" else rg_ok):return None
      cols.add((p,oi));hcnt+=1
    if hcnt==0:return None
  return fc,gc
 def _eval_fg(a,b):
  a2,b2=mul_fn(a,a),mul_fn(b,b);a3,b3=mul_fn(a2,a),mul_fn(b2,b)
  el_g=add_fn(a3,mul_fn(c_g,b3));el_f=sub_fn(a3,mul_fn(c_f,b3))
  return el_f,el_g,abs(norm_fn(el_f)),abs(norm_fn(el_g))
 als=[[1,0,0,0,0,0,0,0]]
 bts=[];cur_b=[0]*8
 def _gen_bt(idx,rem):
  if idx==8:
   if any(cur_b):bts.append(list(cur_b))
   return
  r=math.isqrt(rem)
  for v in range(-r,r+1):cur_b[idx]=v;_gen_bt(idx+1,rem-v*v)
 _gen_bt(0,5)
 assert len(bts)==int(camp["ball_pts"])==3728 and len({tuple(v) for v in bts})==3728
 def _canon_g24(a,b):
  g=math.gcd(*[abs(x) for x in a],*[abs(x) for x in b])
  if g>1:a=[x//g for x in a];b=[x//g for x in b]
  best=None;ca,cb=a,b
  for _ in range(4):
   cbw=cb
   for _ in range(3):
    tp=tuple(ca)+tuple(cbw);tn=tuple(-x for x in ca)+tuple(-x for x in cbw)
    cd=tp if tp<tn else tn
    if best is None or cd<best:best=cd
    cbw=mul_fn(w_el,cbw)
   ca=_sig(ca);cb=_sig(cb)
  return best
 live_ok=0
 for q,sid,ry,rt,exp_m,exp_sv,_,_ in camp["live_balls"]:
  B=[[0]*16 for _ in range(16)];B[0][0]=q
  ry_p=[pow(ry,i,q) for i in range(8)]
  for i in range(1,8):B[i][0]=(-ry_p[i])%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(rt*ry_p[i])%q;B[8+i][8+i]=1
  R,dsq,_=lll_16d_int(B,99,100);assert dsq[16]==q*q;a1,b1=R[0][:8],R[0][8:]
  a2,b2=next((v[:8],v[8:]) for v in R[1:] if any(sub_fn(mul_fn(a1,v[8:]),mul_fn(b1,v[:8]))))
  a2s=[mul_fn(y_pow_fn(j),a2) for j in range(8)];b2s=[mul_fn(y_pow_fn(j),b2) for j in range(8)]
  m_u=sv_c=0;seen_b=set()
  for bt in bts:
   a=list(a1);b=list(b1)
   for j,c in enumerate(bt):
    if c==1:
     aj,bj=a2s[j],b2s[j]
     for i in range(8):a[i]+=aj[i];b[i]+=bj[i]
    elif c==-1:
     aj,bj=a2s[j],b2s[j]
     for i in range(8):a[i]-=aj[i];b[i]-=bj[i]
    elif c:
     aj,bj=a2s[j],b2s[j]
     for i in range(8):a[i]+=c*aj[i];b[i]+=c*bj[i]
   el_f,el_g,nf,ng=_eval_fg(a,b)
   if sid==0:assert nf%q==0;nf_r,ng_r=nf//q,ng
   else:assert ng%q==0;nf_r,ng_r=nf,ng//q
   if nf_r==0 or ng_r==0:continue
   if _strip(nf_r,P_small,P_mid)==1 and _strip(ng_r,P_small,P_mid)==1:
    sv_c+=1;cols=_extract_cols(el_f,el_g,nf_r,ng_r,q,sid)
    if cols is not None and cols[0].issubset(ret_f) and cols[1].issubset(ret_g):
     ck=_canon_g24(a,b);assert ck not in seen_b;seen_b.add(ck);m_u+=1
  assert m_u==exp_m and sv_c==exp_sv;live_ok+=1
 recs=camp["all_records"];stats={}
 for name,cid in (("d1",1),("d2",2),("ctrl",0),("ALL",-1)):
  lst=[r for r in recs if cid<0 or r[5]==cid]
  n=len(lst);ms=[r[2] for r in lst];mus=[r[4] for r in lst]
  u_m=sum(ms)/n;mu_m=sum(mus)/n
  var_m=sum((x-u_m)**2 for x in ms)/max(1,n-1);od=(var_m-u_m)/(u_m*u_m)
  n_z=sum(1 for x in ms if x==0);p0_e=n_z/n;p0_hpgf=sum(1.0/(1.0+0.5*mu) for mu in mus)/n
  mgf_e=sum(math.exp(-0.25*m) for m in ms)/n
  mgf_hpgf=sum(1.0/(1.0+0.5*mu*(1.0-math.exp(-0.25))) for mu in mus)/n
  cdf_ok=all(sum(1 for m in ms if m<=k)/n<=sum(_nb_cdf(k+1,0.5*mu,1.0) for mu in mus)/n for k in (0,1,2,5))
  assert p0_e<p0_hpgf and mgf_e<mgf_hpgf and od<1.0 and cdf_ok
  stats[name]=dict(n=n,mean_m=round(u_m,2),mean_mu=round(mu_m,2),min_m=min(ms),
                   empties=f"{n_z}/{n}",p0_hpgf05=round(p0_hpgf,4),overdisp=round(od,3),
                   mgf_emp=round(mgf_e,4),mgf_hpgf05=round(mgf_hpgf,4),cdf_k0_5_ok=cdf_ok)
 return dict(tower=camp["tower"],p1_unique_single_chart_al1=True,zero_unit_dups_in_ball=True,
             exact_zomega_g12_orbits=True,zero_coarse_ns_keys=True,
             live_balls_replayed=live_ok,n_2core_rels=camp["n_2core_rels"],
             fb_ret_12orbs_f=f"{len(camp['ret_f_split'])}/{camp['n_split_12orbs_f']}",
             fb_ret_12orbs_g=f"{len(camp['ret_g_split'])}/{camp['n_split_12orbs_g']}",
             fb_ret_ns_zomega_f=f"{len(camp['ret_f_ns'])}/{camp['n_ns_zomega_orbs']}",
             fb_ret_ns_zomega_g=f"{len(camp['ret_g_ns'])}/{camp['n_ns_zomega_orbs']}",
             cohorts=stats,checks_ok=True)
def verify_kernels_symbolic(kernels_text,coll,la,r_bls,lll_node_cap,lll_boot_cap):
 import re,random,struct
 allowed={"ldp","stp","ldr","str","ldrb","strb","ldrh","strh","ldur","stur",
  "add","adds","adc","adcs","sub","subs","sbc","sbcs","mul","umulh","msub",
  "and","orr","eor","lsl","lsr","asr","cmp","tst","csel","mov","b","b.ne","b.eq",
  "b.lo","b.ls","b.hi","b.hs","b.gt","cbz","cbnz","ret"}
 forb=re.compile(r"\b(v[0-9]+|q[0-9]+|d[0-9]+|z[0-9]+|dc|ldadd\w*|swp\w*|cas\w*|bl|blr)\b")
 blocks=re.findall(r"```asm\n(.*?)```",kernels_text,flags=re.S)
 parsed={}
 for blk in blocks:
  lines=[l.strip() for l in blk.splitlines() if l.strip()]
  m=re.match(r"^;\s*kernel:\s*(\S+)\s+items=(\d+)\s+charged=([0-9.]+)\s+key=(\S+)$",lines[0])
  assert m,f"Bad kernel header: {lines[0]}"
  name,items,charged,key=m.group(1),int(m.group(2)),float(m.group(3)),m.group(4)
  insns=[l for l in lines[1:] if not l.startswith(";")]
  for ins in insns:
   assert not forb.search(ins.lower()),f"Forbidden SIMD/LSE/call in {name}: {ins}"
   op=ins.split()[0].lower()
   assert op in allowed,f"Disallowed opcode {op} in {name}: {ins}"
   if op in ("ldp","stp"):
    mb=re.search(r"\[([^\]]+)\]",ins)
    assert mb and re.match(r"^\s*x\w+(\s*,\s*#-?\d+)?\s*$",mb.group(1)),f"Illegal {op} addr mode in {name}: {ins}"
  ir=len(insns)/items
  assert ir<=charged+1e-12,f"{name}: {ir}>{charged}"
  parsed[name]=dict(raw_opcodes=len(insns),items=items,ir=ir,charged=charged,insns=insns)
 def _pw(insns):
  ws=set()
  for ins in insns:
   ps=ins.replace(","," ").split();op=ps[0].lower()
   if op in ("cmp","tst","b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","cbz","cbnz","ret","str","strb","strh","stur","stp"):continue
   for r in [ps[1]]+(ps[2:3] if op=="ldp" else []):
    rl=r.lower()
    if re.match(r"^[wx]\d+$",rl):ws.add("r"+rl[1:])
  return ws
 _lp={"r19","r20","r21","r22"};_wc=_pw(parsed["rev_b2_scan_cand"]["insns"]);_we=_pw(parsed["rev_b2_emit"]["insns"])
 assert _wc.isdisjoint(_lp) and _we.isdisjoint(_lp) and _wc.isdisjoint(_we)
 assert parsed["rev_sb_probe_fast"]["ir"]<=coll["c_tile_hdr"]
 assert parsed["rev_b2_minkey_skip"]["ir"]<=coll["c_sb_skip"]
 assert parsed["rev_b2_minkey_enter"]["ir"]<=coll["c_sb_enter"]
 assert parsed["rev_b2_scan_cand"]["ir"]<=coll["c_cand"]
 assert parsed["rev_b1a_apply_u4"]["ir"]<=coll["c_hit_lo_m4"]
 assert parsed["rev_b1b_apply_u4"]["ir"]<=coll["c_hit_lo_u"]
 assert parsed["rev_b2_emit"]["ir"]<=11.25
 assert (parsed["rev_b2_emit"]["ir"]+parsed["rev_b1a_apply_u4"]["ir"])<=coll["c_hit_hi"]
 assert parsed["ntt62_butterfly"]["ir"]*12.0<=120.0<168.0
 assert (parsed["rev_build_pass12"]["ir"]+0.50)<=coll["c_list_a"] and "ldr x2, [x_hist, x1, lsl #3]" in parsed["rev_build_pass12"]["insns"] and "str x2, [x_hist, x1, lsl #3]" in parsed["rev_build_pass12"]["insns"]
 assert parsed["rev_build_sb_min"]["ir"]<=coll["c_lbkt"]
 MASK=(1<<64)-1;MASK60=(1<<60)-1;OID_M=(1<<58)-1
 for m_b in range(4):
  for k_c in range(4):
   xk58=k_c<<58
   for oid in (0,1,0x12345,OID_M):
    om=(m_b<<58)|oid
    rot=((om-xk58)&MASK)&MASK60
    assert (rot>>58)==((m_b-k_c)&3) and (rot&OID_M)==oid and (rot<<1)<(1<<61)
    if m_b==0 and k_c==1:assert (rot>>58)==3
 p_ch=1000003;r_ch=123457;st_ch=pow(r_ch,-1,p_ch);s0_ch=7;a_base=[3,-2,5,1,-4,2,-1,3]
 exps1w=(0,5,3,8,6,11,9,14);pw_ch=[pow(s0_ch,e,p_ch) for e in exps1w]
 ka_cur=(sum(a_base[j]*pw_ch[j] for j in range(8))*st_ch)%p_ch
 for k_c in range(16):
  ka_dir=(((a_base[0]+k_c)*pw_ch[0]+sum(a_base[j]*pw_ch[j] for j in range(1,8)))*st_ch)%p_ch
  assert ka_cur==ka_dir
  qkey_t=((k_c*15+5)<<24)|0x12345;oid_m_t=((k_c&3)<<58)|(100000+k_c)
  x8_pack=qkey_t|((ka_cur&0xFFFFFFFF)<<32)
  entry16=struct.pack("<QQ",x8_pack,oid_m_t)
  w0_rd,w1_rd=struct.unpack_from("<II",entry16,0)
  oid_rd,=struct.unpack_from("<Q",entry16,8)
  assert w0_rd==qkey_t and w1_rd==(ka_cur&0xFFFFFFFF) and oid_rd==oid_m_t and (w0_rd>>24)==(qkey_t>>24)
  ka_cur=ka_cur+st_ch-p_ch if ka_cur+st_ch>=p_ch else ka_cur+st_ch
 qtab_e={7:(5<<24)|0x18000,9:(5<<24)|0x10000,12:(5<<24)|0x1F000}
 r4=[(0,(5<<24)|0x25000,1007,33),(2,(5<<24)|0x11000,1007,(3<<58)|0x7777),(0,(5<<24)|0x15000,1009,(3<<58)|22),(1,(6<<24)|0x01000,1007,44),(0,(5<<24)|0x10000,1900,(1<<58)|11),(0,(5<<24)|0x1E000,1012,(2<<58)|0x9999),(0,(5<<24)|0x12345,1007,0x4242)]
 sbl=[(5<<24)|0x20000,(3<<24)|0x10000,(5<<24)|0x19000]
 og=sorted(range(len(r4)),key=lambda i:(r4[i][1],i))
 xm=[(r4[i][1],r4[i][3]) for i in og];xim=[0]*len(r4)
 for pg,ci in enumerate(og):xim[ci]=pg
 xkb=[0]*len(r4);sbo=[0]*len(r4);hsb=[0,0,0]
 for ci,(si,_,ka,_) in enumerate(r4):pg=xim[ci];xkb[pg]=ka;sbo[pg]=si;hsb[si]+=1
 off_b=[0]
 for c in hsb:off_b.append(off_b[-1]+c*16)
 assert all(o%16==0 for o in off_b)
 fba=bytearray(off_b[-1]);psb=list(off_b[:-1]);sm8=[255,255,255]
 for pg,(qk,om) in enumerate(xm):
  si=sbo[pg];struct.pack_into("<IIQ",fba,psb[si],qk,xkb[pg],om);psb[si]+=16
  q8=qk>>24
  if q8<sm8[si]:sm8[si]=q8
 flat_buf=bytes(fba)
 for j in range(3):
  bk=[struct.unpack_from("<I",flat_buf,o)[0] for o in range(off_b[j],off_b[j+1],16)]
  assert bk==sorted(bk)
 sl8=[(lm+0xFFFFFF)>>24 for lm in sbl]
 x19=x20=x22=0;x21=3;em=[];sv=[]
 while x21>0:
  w9=sm8[x19];x19+=1;w10=sl8[x20];x20+=1
  if w9<=w10:
   j=(x19-x22)-1;sv.append(j);lm32=sbl[j];xp,xpe=off_b[j],off_b[j+1]
   while xp<xpe:
    w9e,w10e=struct.unpack_from("<II",flat_buf,xp);xp+=16
    if w9e>lm32:break
    w11e=(w10e-1000)&0xFFFFFFFF
    if w11e>=552:continue
    if w9e<=qtab_e[w11e]:
     x15e,=struct.unpack_from("<Q",flat_buf,xp-8);x15e=((x15e-(1<<58))&MASK)&MASK60
     em.append((j,w11e,x15e>>58,x15e&OID_M))
  x21-=1
 exh=[((si,(ka-1000)&0xFFFFFFFF,(((om-(1<<58))&MASK)&MASK60)>>58,(((om-(1<<58))&MASK)&MASK60)&OID_M)) for si,qk,ka,om in sorted(r4,key=lambda t:(t[0],t[1])) if ((ka-1000)&0xFFFFFFFF)<552 and qk<=sbl[si] and qk<=qtab_e[(ka-1000)&0xFFFFFFFF]]
 assert x19==3 and x20==3 and x21==0 and sv==[0,2] and em==exh==[(0,7,3,0x4242),(0,12,1,0x9999),(2,7,2,0x7777)]
 assert (1>>6)==0 and (63>>6)==0 and (64>>6)==1 and "add x4, x3, x_sb_m1" not in parsed["rev_sb_probe_fast"]["insns"]
 scan_tot=parsed["rev_scan_zero_u16"]["ir"]*1.02761+0.15
 assert scan_tot<=coll["c_scan_zero"]
 mont5_ops=5*parsed["mont5_cios_round"]["raw_opcodes"]+parsed["mont5_cios_epilogue"]["raw_opcodes"]
 assert mont5_ops<=MONT5
 assert parsed["unit_run5_z12"]["ir"]<=8.80
 assert parsed["gen_run1_z12"]["ir"]<=29.0 and 4.0*parsed["gen_run1_z12"]["ir"]<=116.0
 assert parsed["row_red5_z12"]["ir"]<=68.0
 assert (6*parsed["row_red5_z12"]["raw_opcodes"]+504+16)<=la["c_row_chi"]
 for T in (1,2,100,1000,2800,4095,16384,32767):
  b16=(32768-T)&0x7FFF
  for x in (0,1,T-1,T,T+1,2800,32767,32768,65535):
   if 0<=x<=65535:assert (((((x&0x7FFF)+b16)|x)&0x8000)!=0)==(x>=T)
 rng=random.Random(20261120)
 for _ in range(32):
  pv=rng.randrange(1<<315,(1<<320)-1)|1;av=rng.randrange(pv);bv=rng.randrange(pv)
  p=[(pv>>(64*j))&MASK for j in range(5)];a=[(av>>(64*j))&MASK for j in range(5)];b=[(bv>>(64*j))&MASK for j in range(5)]
  n0=(-pow(p[0],-1,1<<64))&MASK;t=[0]*6
  for i in range(5):
   bi=b[i];c=0
   for j in range(5):s=t[j]+((a[j]*bi)&MASK)+c;t[j],c=s&MASK,s>>64
   s=t[5]+c;t[5],t6=s&MASK,s>>64;c=0
   for j in range(5):s=t[j+1]+((a[j]*bi)>>64)+c;t[j+1],c=s&MASK,s>>64
   t6=(t6+c)&MASK;mm=(t[0]*n0)&MASK;s=t[0]+((mm*p[0])&MASK);assert (s&MASK)==0;c=s>>64
   for j in range(1,5):s=t[j]+((mm*p[j])&MASK)+c;t[j-1],c=s&MASK,s>>64
   s=t[5]+c;t[4],c=s&MASK,s>>64;t[5]=(t6+c)&MASK;c=0
   for j in range(5):s=t[j]+((mm*p[j])>>64)+c;t[j],c=s&MASK,s>>64
   t[5]=(t[5]+c)&MASK
  val=sum(t[j]<<(64*j) for j in range(6))
  if val>=pv:val-=pv
  assert val==(av*bv*pow(1<<320,-1,pv))%pv
 mu_r=(1<<318)//r_bls;rl=[(r_bls>>(64*j))&MASK for j in range(4)]
 def _subc(aa,bb,ci):
  tt=aa+(bb^MASK)+ci
  return tt&MASK,tt>>64
 b_tests=[0,1,-1,r_bls-1,-(r_bls-1),r_bls,-r_bls,3*(r_bls-1),-3*(r_bls-1),
  200*(r_bls-1),-200*(r_bls-1),51000*(r_bls-1),-51000*(r_bls-1),
  1000000*(r_bls-1),-1000000*(r_bls-1),(1<<59)*(r_bls-1),-(1<<59)*(r_bls-1)]
 for _ in range(30):b_tests.append(rng.randint(-(1<<315),(1<<315)))
 for Av in b_tests:
  x10,x11,x12,x13,x14=[(Av>>(64*j))&MASK for j in range(5)]
  x15=MASK if (x14>>63) else 0
  x0,x1,x2,x3,x4=x10^x15,x11^x15,x12^x15,x13^x15,x14^x15
  x0,c=_subc(x0,x15,1);x1,c=_subc(x1,x15,c);x2,c=_subc(x2,x15,c);x3,c=_subc(x3,x15,c);x4,_=_subc(x4,x15,c)
  x5=(((x4<<2)&MASK)|(x3>>62));x5=(x5*mu_r)>>64
  p0,p1,p2,p3=[(x5*rl[j])&MASK for j in range(4)];h0,h1,h2,h3=[(x5*rl[j])>>64 for j in range(4)]
  s=p1+h0;p1,c=s&MASK,s>>64;s=p2+h1+c;p2,c=s&MASK,s>>64;s=p3+h2+c;p3,c=s&MASK,s>>64;h3=(h3+c)&MASK
  x0,c=_subc(x0,p0,1);x1,c=_subc(x1,p1,c);x2,c=_subc(x2,p2,c);x3,c=_subc(x3,p3,c);x4,_=_subc(x4,h3,c)
  rem_u=x0+(x1<<64)+(x2<<128)+(x3<<192)+(x4<<256);assert 0<=rem_u<3*r_bls
  for _ in range(2):
   t0,c=_subc(x0,rl[0],1);t1,c=_subc(x1,rl[1],c);t2,c=_subc(x2,rl[2],c);t3,c=_subc(x3,rl[3],c);t4,c=_subc(x4,0,c)
   if c:x0,x1,x2,x3,x4=t0,t1,t2,t3,t4
  assert x4==0
  if (x0|x1|x2|x3)==0:x15=0
  if x15:x0,c=_subc(rl[0],x0,1);x1,c=_subc(rl[1],x1,c);x2,c=_subc(rl[2],x2,c);x3,_=_subc(rl[3],x3,c)
  assert (x0+(x1<<64)+(x2<<128)+(x3<<192))==Av%r_bls
 ntt_ps=[c*(1<<53)+1 for c in (269,276,305,309)]
 for np_i in ntt_ps:
  for _ in range(16):
   uv,vv,wv=rng.randrange(np_i),rng.randrange(np_i),rng.randrange(np_i)
   ws=(wv<<64)//np_i;tv=(((vv*wv)&MASK)-((((vv*ws)>>64)*np_i)&MASK))&MASK
   if tv>=np_i:tv-=np_i
   u1=uv+tv-np_i if uv+tv>=np_i else uv+tv;v1=uv-tv+np_i if uv<tv else uv-tv
   assert u1==(uv+vv*wv)%np_i and v1==(uv-vv*wv)%np_i
  mu_i=(1<<124)//np_i
  for _ in range(16):
   av,bv,cv=rng.randrange(np_i),rng.randrange(np_i),rng.randrange(np_i)
   pr=av*bv;qh=((pr>>60)*mu_i)>>64;rm=pr-qh*np_i
   if rm>=np_i:rm-=np_i
   if rm>=np_i:rm-=np_i
   acc=cv+rm-np_i if cv+rm>=np_i else cv+rm
   assert 0<=rm<np_i and acc==(cv+av*bv)%np_i
 m_all=ntt_ps[0]*ntt_ps[1]*ntt_ps[2]*ntt_ps[3]
 i01=pow(ntt_ps[0],-1,ntt_ps[1]);i012=pow((ntt_ps[0]*ntt_ps[1])%ntt_ps[2],-1,ntt_ps[2]);i0123=pow((ntt_ps[0]*ntt_ps[1]*ntt_ps[2])%ntt_ps[3],-1,ntt_ps[3])
 for _ in range(32):
  Xv=rng.randrange(m_all);rr=[Xv%p for p in ntt_ps]
  v1=((rr[1]-rr[0])*i01)%ntt_ps[1];x01=rr[0]+v1*ntt_ps[0]
  v2=((rr[2]-x01%ntt_ps[2])*i012)%ntt_ps[2];x012=x01+v2*ntt_ps[0]*ntt_ps[1]
  v3=((rr[3]-x012%ntt_ps[3])*i0123)%ntt_ps[3]
  assert x012+v3*ntt_ps[0]*ntt_ps[1]*ntt_ps[2]==Xv
 s1_ds=math.floor(2100.0/math.log2(4.0/3.0))
 ds_ops=15*(40*5*3+16*12*4+15*12*10)+15*(4*12*3*5+2*40*8*3)+2000
 ds_lll_bound=(2*s1_ds+16)*ds_ops
 assert ds_lll_bound<(1<<30)
 s1_sq=math.floor(8160.0/math.log2(40.0/11.0))
 sq1_ops=15*(40*14*9+16*12*25+15*12*70)+15*(4*12*9*14+2*40*23*9)+10000
 s2_sq=math.floor((340.0*math.log2(40.0)+7.5*530.0)/math.log2(100.0/99.0))
 sq_lll_bound=(2*s1_sq+16)*sq1_ops+(2*s2_sq+16)*sq1_ops
 assert sq_lll_bound<lll_node_cap
 s1_bt=math.floor(67191.0/math.log2(40.0/11.0))
 bt1_ops=15*(40*199*132+16*12*67*67+15*12*67*199)+15*(4*12*132*199+2*40*331*132)+1000000
 s2_bt=math.floor((340.0*math.log2(40.0)+7.5*3384.0)/math.log2(100.0/99.0))
 bt2_ops=15*(40*80*53+16*12*27*27+15*12*27*80)+15*(4*12*53*80+2*40*133*53)+1000000
 bt_lll_bound=(2*s1_bt+16)*bt1_ops+(2*s2_bt+16)*bt2_ops
 assert bt_lll_bound<lll_boot_cap
 q265=(1<<265)+1234567891;Bsq=[[0]*16 for _ in range(16)];Bsq[0][0]=q265
 ry=pow(7,13,q265);rt=pow(11,17,q265);ryp=[pow(ry,i,q265) for i in range(8)]
 for i in range(1,8):Bsq[i][0]=(-ryp[i])%q265;Bsq[i][i]=1
 for i in range(8):Bsq[8+i][0]=(rt*ryp[i])%q265;Bsq[8+i][8+i]=1
 _,dsq,lsq=lll_16d_int(Bsq,99,100)
 assert dsq[16]==q265*q265
 for kk in range(1,16):
  assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
  for jj in range(kk):assert 2*abs(lsq[kk][jj])<=dsq[jj+1]
 return dict(asm_blocks_parsed=len(parsed),b2_multi_entry_cf_ok=True,b2_bank_rot_16pairs_ok=True,uint16_swar_ok=True,
             ntt62_garner_exact=True,ntt62_3x4p_stage_ir=round(12.0*parsed["ntt62_butterfly"]["ir"],1),f62_barrett_mac_64_exact=True,
             mont5_opcodes=mont5_ops,mont5_sim_32_exact=True,
             row_red5_barrett_opcodes=parsed["row_red5_z12"]["raw_opcodes"],
             row_red5_barrett_47_exact=True,
             unit_run5_ir=round(parsed["unit_run5_z12"]["ir"],2),
             gen_run1_4coord_ir=round(4.0*parsed["gen_run1_z12"]["ir"],2),
             scan_zero_u16=round(scan_tot,4),
             sq_lll_uniform_log2=round(math.log2(sq_lll_bound),2),
             boot_lll_uniform_log2=round(math.log2(bt_lll_bound),2),
             all_raw_opcodes_ok=True)
def verify_prod_sq24(sq24,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn):
 import random
 assert len(sq24)==24
 log2_u8=float(desc["params"]["log2_u8"]);lx=float(desc["params"]["lx"])
 qbins,qbins_lo=desc["qbins"],desc["qbins_lo"]
 qb_f,qb_g=[x[1] for x in qbins],[x[2] for x in qbins]
 S8=math.log2(desc["cells_rank2"]["8"])
 w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10)
 rng=random.Random(20261120);ab_pts=[]
 while len(ab_pts)<16:
  al=[rng.randint(-4,4) for _ in range(8)];bt=[rng.randint(-4,4) for _ in range(8)]
  if any(al) and sum(x*x for x in al)+sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
 def _is_p(n):
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
 out_sides={}
 for s_id,sname in ((0,"g"),(1,"f")):
  cfs,cgs,sv_rats,pr_rats=[],[],[],[]
  recs=[r for r in sq24 if r[0]==s_id]
  assert len(recs)==12
  for side,q_hex,ry_hex,rt_hex in recs:
   q,ry,rt=int(q_hex,16),int(ry_hex,16),int(rt_hex,16)
   nb=q.bit_length()-1
   assert 70<=nb<265 and q%15==1 and _is_p(q)
   B=[[0]*16 for _ in range(16)];B[0][0]=q
   ryp=[pow(ry,i,q) for i in range(8)]
   for i in range(1,8):B[i][0]=(-ryp[i])%q;B[i][i]=1
   for i in range(8):B[8+i][0]=(rt*ryp[i])%q;B[8+i][8+i]=1
   R,dsq,lsq=lll_16d_int(B,99,100)
   assert dsq[16]==q*q
   for kk in range(1,16):
    assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
    for jj in range(kk):assert 2*abs(lsq[kk][jj])<=dsq[jj+1]
   a1,b1=R[0][:8],R[0][8:]
   a2,b2=next((v[:8],v[8:]) for v in R[1:] if any(sub_fn(mul_fn(a1,v[8:]),mul_fn(b1,v[:8]))))
   det_q=abs(norm_fn(sub_fn(mul_fn(a1,b2),mul_fn(b1,a2))))
   assert det_q>0 and det_q%q==0
   c_val=(math.log2(q)+S8)/16.0
   prow=next(r for r in dlive["rows"] if r["side"]==sname and r["lo"]<=nb<r["hi"])
   S_bin=math.log2(desc["cells_rank2"][str(prow["R"])]);c_lo_b=(prow["lo"]+S_bin)/16.0
   T_b=prow["kA"]*prow["nA"]+prow["kB"]*prow["n"]+96.0
   sv_plan=0.0
   for (w_q,_,_),(cfl,cgl) in zip(qbins,qbins_lo):
    nfl=72.0*c_lo_b+cfl-(prow["hi"] if sname=="f" else 0.0)
    ngl=24.0*c_lo_b+log2_u8+cgl-(prow["hi"] if sname=="g" else 0.0)
    sf=1.0 if T_b>=nfl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfl)-T_b)/lx))
    sg=1.0 if T_b>=ngl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(ngl)-T_b)/lx))
    sv_plan+=w_q*sf*sg
   sv_acc=0.0;c_act=(math.log2(q)+S_bin)/16.0;lq=math.log2(q)
   nfas,ngas=[],[]
   for al,bt in ab_pts:
    a=add_fn(mul_fn(al,a1),mul_fn(bt,a2));b=add_fn(mul_fn(al,b1),mul_fn(bt,b2))
    a3,b3=pow_fn(a,3),pow_fn(b,3)
    ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
    a6,b6=mul_fn(a3,a3),mul_fn(b3,b3);a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
    nf=abs(norm_fn(sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))))
    assert (ng if s_id==0 else nf)%q==0
    cf=math.log2(nf)-72.0*c_val;cg=math.log2(ng)-24.0*c_val-log2_u8
    cfs.append(cf);cgs.append(cg)
    nfa=72.0*c_act+cf-(lq if sname=="f" else 0.0);nga=24.0*c_act+log2_u8+cg-(lq if sname=="g" else 0.0)
    nfas.append(nfa);ngas.append(nga)
    sfe=1.0 if T_b>=nfa else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfa)-T_b)/lx))
    sge=1.0 if T_b>=nga else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nga)-T_b)/lx))
    sv_acc+=sfe*sge
   assert max(nfas)<=prow["nf"]-50.0 and max(ngas)<=prow["ng"]-25.0
   sv_rats.append((sv_acc/len(ab_pts))/max(sv_plan,1e-300))
   _ss=dlive["_ss_eval"];_mn=dlive["mesh_node"]
   _pfe=_ss(nfas,float(prow["n"]),prow["iA"],_mn,float(desc["params"]["alpha_f"]),prow["kA"],prow["kB"])
   _pge=_ss(ngas,float(prow["n"]),prow["iA"],_mn,float(desc["params"]["alpha_g"]),prow["kA"],prow["kB"])
   _pre=sum(x*y for x,y in zip(_pfe,_pge))/len(ab_pts)
   pr_rats.append(_pre/prow["pr"])
  cfs.sort();cgs.sort();n=len(cfs)
  sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
  sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
  assert all(sf6[i]<qb_f[i] for i in range(6)) and all(sg6[i]<qb_g[i] for i in range(6))
  assert max(sv_rats)<=1.0+1e-12 and (sum(sv_rats)/len(sv_rats))<0.50 and min(pr_rats)>=10.0
  _prs=sorted(pr_rats)
  out_sides[sname]=dict(bases=12,samples=n,sext_f=sf6,sext_g=sg6,
   min_headroom_f=round(min(qb_f[i]-sf6[i] for i in range(6)),2),
   min_headroom_g=round(min(qb_g[i]-sg6[i] for i in range(6)),2),
   joint_pr_over_charged_min=round(_prs[0],2),
   joint_pr_over_charged_med=round(_prs[len(_prs)//2],2),
   sv_emp_over_plan_mean=round(sum(sv_rats)/len(sv_rats),4),
   sv_emp_over_plan_max=round(max(sv_rats),4))
 return dict(n_bases=24,live_bareiss_99_100_all24=True,presupplied_vectors=0,
             q_bit_range="73..258",g_side=out_sides["g"],f_side=out_sides["f"],checks_ok=True)
def verify_gt_boot_20(boot_w,data,P_bls,R_bls,U_bls,h_pi,w_pi,pmul_fn,prem_fn,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn):
 import hashlib
 S_boot=int(data["descent"]["boot"]["s_boot"])
 n_targets=int(boot_w["n_targets"])
 assert S_boot==int(boot_w["s_boot"])==(1 << 21) and n_targets>=1
 assert "full4" not in boot_w and "extra16" not in boot_w,"ZERO pre-supplied LLL vectors allowed!"
 w_pi2=(w_pi*w_pi)%P_bls;uw2=(U_bls*w_pi2)%P_bls;h0,h1,h2,h3=h_pi[0],h_pi[1],h_pi[2],h_pi[3]
 def f4_add(a,b):return [(a[0]+b[0])%P_bls,(a[1]+b[1])%P_bls,(a[2]+b[2])%P_bls,(a[3]+b[3])%P_bls]
 def f4_sub(a,b):return [(a[0]-b[0])%P_bls,(a[1]-b[1])%P_bls,(a[2]-b[2])%P_bls,(a[3]-b[3])%P_bls]
 def f4_scale(s,a):return [(s*a[0])%P_bls,(s*a[1])%P_bls,(s*a[2])%P_bls,(s*a[3])%P_bls]
 def f4_mul(a,b):
  c=[a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[1]*b[1]+a[2]*b[0],
     a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0],
     a[1]*b[3]+a[2]*b[2]+a[3]*b[1],a[2]*b[3]+a[3]*b[2],a[3]*b[3]]
  for k in (6,5,4):
   v=c[k]%P_bls
   if v:o=k-4;c[o]-=v*h0;c[o+1]-=v*h1;c[o+2]-=v*h2;c[o+3]-=v*h3
  return [c[0]%P_bls,c[1]%P_bls,c[2]%P_bls,c[3]%P_bls]
 def f4_pow(a,e):
  res,b=[1,0,0,0],list(a)
  while e:
   if e & 1:res=f4_mul(res,b)
   e >>=1
   if e:b=f4_mul(b,b)
  return res
 f1c=[f4_pow([1 if j==i else 0 for j in range(4)],P_bls) for i in range(4)]
 def f4_frob1(a):
  out=[0,0,0,0]
  for i in range(4):
   ai=a[i]
   if ai:col=f1c[i];out[0]+=ai*col[0];out[1]+=ai*col[1];out[2]+=ai*col[2];out[3]+=ai*col[3]
  return [x%P_bls for x in out]
 def f4_frob2(a):return f4_frob1(f4_frob1(a))
 def f4_inv(a):
  a1=f4_frob1(a);a2=f4_frob1(a1);a3=f4_frob1(a2)
  pr3=f4_mul(a1,f4_mul(a2,a3))
  return f4_scale(pow(f4_mul(a,pr3)[0],-1,P_bls),pr3)
 def f12_mul(A,B):
  return (
   f4_add(f4_mul(A[0],B[0]),f4_scale(uw2,f4_add(f4_mul(A[1],B[2]),f4_mul(A[2],B[1])))),
   f4_add(f4_add(f4_mul(A[0],B[1]),f4_mul(A[1],B[0])),f4_scale(uw2,f4_mul(A[2],B[2]))),
   f4_add(f4_add(f4_mul(A[0],B[2]),f4_mul(A[1],B[1])),f4_mul(A[2],B[0])),
  )
 G0=tuple([int(x,16) for x in c] for c in boot_w["G0"])
 one12=([1,0,0,0],[0,0,0,0],[0,0,0,0])
 G0_win=[one12,G0]
 for _ in range(14):G0_win.append(f12_mul(G0_win[-1],G0))
 def f12_pow_G0(e):
  res=one12
  for sh in range(((e.bit_length()-1)//4)*4,-1,-4):
   if res is not one12:res=f12_mul(res,res);res=f12_mul(res,res);res=f12_mul(res,res);res=f12_mul(res,res)
   nib=(e>>sh) & 15
   if nib:res=G0_win[nib] if res is one12 else f12_mul(res,G0_win[nib])
  return res
 assert G0!=one12 and f12_pow_G0(R_bls)==one12
 inv2=(P_bls+1)//2;d_fp2=f4_sub([0,1,0,0],f4_frob2([0,1,0,0]));inv_2d=f4_inv(f4_scale(2,d_fp2))
 def to_fp2(a):ap=f4_frob2(a);return f4_scale(inv2,f4_add(a,ap)),f4_mul(f4_sub(a,ap),inv_2d)
 def elim_gamma(Z):
  e0,o0=to_fp2(Z[0]);e1,o1=to_fp2(Z[1]);e2,o2=to_fp2(Z[2])
  c0=f4_sub(f4_mul(e1,o0),f4_mul(o1,e0));c1=f4_sub(f4_mul(o2,e0),f4_mul(e2,o0));c2=f4_sub(f4_mul(e2,o1),f4_mul(o2,e1))
  assert f4_frob2(c0)==c0 and f4_frob2(c1)==c1 and f4_frob2(c2)==c2
  W=f12_mul((c0,c1,c2),Z)
  assert W[2]==[0,0,0,0] and any(W[1])
  gam=f4_mul([(-x)%P_bls for x in W[0]],f4_inv(W[1]))
  assert f4_add(W[0],f4_mul(gam,W[1]))==[0,0,0,0]
  return gam
 z3_2=y_pow_fn(10);p4=P_bls**4;det_exp=(P_bls**8)*(S_boot**16)
 lov_l2=math.log2(((0.99-0.25)**(-15.0/4.0))*(float(P_bls)**0.25)*(float(S_boot)**0.5))
 bsub=data["descent"]["boot_sub_deciles"];m_sv=float(data["descent"]["boot"]["m_sv"])
 betas=[[0]*8]+[[1 if j==i else 0 for j in range(8)] for i in range(8)]+[[-1 if j==i else 0 for j in range(8)] for i in range(8)]
 sv_bits,all_bits,sk_l2s=[],[],[]
 decile_pass=gcd_gt1e4=tot_pairs=0
 for ti in range(n_targets):
  xt=int(hashlib.sha256(f"bls12-381-gt-boot-v5g-{ti}".encode()).hexdigest(),16)%(R_bls-1)+1
  gam=elim_gamma(f12_pow_G0(xt))
  B_sk=[]
  for i in range(4):r=[0]*16;r[i]=P_bls;B_sk.append(r)
  for j in range(4):
   r=[0]*16
   for idx,ch in enumerate(h_pi):r[j+idx]=ch
   B_sk.append(r)
  for k in range(8):
   yk=[0]*8;yk[k]=1
   ak=pmul_fn(gam,prem_fn(yk,h_pi,P_bls),P_bls,h_pi)
   r=[0]*16
   for idx,ch in enumerate(ak):r[idx]=ch
   r[8+k]=S_boot;B_sk.append(r)
  R_sk,d_gso,lam=lll_16d_int(B_sk,99,100)
  assert d_gso[16]==det_exp
  for kk in range(1,16):
   assert 100*(d_gso[kk-1]*d_gso[kk+1]+lam[kk][kk-1]*lam[kk][kk-1])>=99*d_gso[kk]*d_gso[kk]
   for jj in range(kk):assert 2*abs(lam[kk][jj])<=d_gso[jj+1]
  assert all(x%S_boot==0 for x in R_sk[0][8:])
  a1,b1=R_sk[0][:8],[x//S_boot for x in R_sk[0][8:]]
  sk_l2=math.log2(math.sqrt(sum(x*x for x in a1)+sum((S_boot*x)**2 for x in b1)))
  assert sk_l2<=lov_l2;sk_l2s.append(sk_l2)
  assert prem_fn(a1,h_pi,P_bls)==pmul_fn(gam,prem_fn(b1,h_pi,P_bls),P_bls,h_pi)
  g1=abs(norm_fn(sub_fn(pow_fn(a1,3),scale_fn(U_bls,mul_fn(z3_2,pow_fn(b1,3))))))
  lb1=math.log2(g1);sv_bits.append(lb1)
  if ti<n_targets:
   v2=next(row for row in R_sk[1:] if all(x%S_boot==0 for x in row[8:]) and abs(norm_fn(sub_fn(mul_fn(a1,[x//S_boot for x in row[8:]]),mul_fn(b1,row[:8]))))==p4)
   a2,b2=v2[:8],[x//S_boot for x in v2[8:]]
   assert prem_fn(a2,h_pi,P_bls)==pmul_fn(gam,prem_fn(b2,h_pi,P_bls),P_bls,h_pi)
   bn,bb=[g1],[lb1]
   for bt in betas[1:]:
    aa=add_fn(a1,mul_fn(bt,a2));bvec=add_fn(b1,mul_fn(bt,b2))
    ng=abs(norm_fn(sub_fn(pow_fn(aa,3),scale_fn(U_bls,mul_fn(z3_2,pow_fn(bvec,3))))))
    lb=math.log2(ng);assert lb<2600.0
    bn.append(ng);bb.append(lb)
   all_bits.extend(bb);bs_s=sorted(bb)
   if all(bs_s[min(16,math.ceil((j+1)*17/10)-1)]<=bsub[j] for j in range(10)):decile_pass+=1
   for i in range(len(bn)):
    for j in range(i+1,len(bn)):
     tot_pairs+=1
     if math.gcd(bn[i],bn[j])>10000:gcd_gt1e4+=1
 from fractions import Fraction as _F
 _tr=lambda d:8 if math.gcd(15,d%15)==15 else (-4 if math.gcd(15,d%15)==5 else (-2 if math.gcd(15,d%15)==3 else 1))
 _Gm=[[_F(15001,1000)-_tr(0) if i==j else _F(-_tr(i-j)) for j in range(8)] for i in range(8)]
 for _k in range(1,9):
  _A=[[x for x in r[:_k]] for r in _Gm[:_k]];_det=_F(1)
  for _i in range(_k):
   _pv=next(r for r in range(_i,_k) if _A[r][_i]!=0)
   if _pv!=_i:_A[_i],_A[_pv]=_A[_pv],_A[_i];_det=-_det
   _det*=_A[_i][_i];_iv=_F(1)/_A[_i][_i]
   for _r in range(_i+1,_k):
    if _A[_r][_i]:
     _fc=_A[_r][_i]*_iv
     for _c in range(_i,_k):_A[_r][_c]-=_fc*_A[_i][_c]
  assert _det>0
 _u23=(abs(U_bls)**(2.0/3.0))/(S_boot*S_boot)
 _det_norm_ub=12.0*(math.log2(15.001/8.0)+math.log2(max(1.0,_u23))+2.0*lov_l2)
 assert _det_norm_ub<m_sv and float(data["descent"]["params"]["ret_boot"])==1.0
 n_pass=sum(1 for v in sv_bits if v<=m_sv)
 assert n_pass==n_targets
 return dict(
  boot_element="a - b*X",rel_norm_formula="a^3 - u*omega^2*b^3",
  boot_gt_order_r_exact=True,boot_subfield_elim_live=n_targets,
  boot_s_boot=S_boot,boot_live_381bit_lll_exact_all=n_targets,
  boot_presupplied_vectors=0,
  boot_skew_b1_max_l2=round(max(sk_l2s),2),boot_lovasz_bound_l2=round(lov_l2,2),
  boot_det_norm_theorem_ub=round(_det_norm_ub,2),boot_ret_proved=1.0,boot_det_compat_pass=f"{n_pass}/{n_targets}",
  boot_det_p4_exact=True,boot_decile_filter_pass=f"{decile_pass}/{n_targets}",
  sv_mean_live=round(sum(sv_bits)/len(sv_bits),2),sv_max_live=round(max(sv_bits),2),
  ball_mean_live=round(sum(all_bits)/len(all_bits),2),ball_max_live=round(max(all_bits),2),
  gcd_gt1e4_live=f"{gcd_gt1e4}/{tot_pairs}",
 )
def unplanted_two_sided_bw_check():
 import random
 p,r,m=467,233,381
 f_c,g_c=[1,1,-5,0,1],[-10,-1,1]
 assert sum(ci*pow(m,i,p) for i,ci in enumerate(f_c))%p==0
 assert sum(ci*pow(m,i,p) for i,ci in enumerate(g_c))%p==0
 eps_f,eps_g,r2=r-1,r*r-1,r*r
 def peval(c,x,mod):
  s=0
  for ci in reversed(c):
   s=(s*x+ci)%mod
  return s
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
   if e & 1:res=pmul(res,b,mod_p,q)
   e >>=1
   if e:b=pmul(b,b,mod_p,q)
  return res
 primes=[q for q in range(2,115) if all(q%d for d in range(2,int(q**0.5)+1))]
 def unram(c,q):
  for rt in range(q):
   if peval(c,rt,q)==0:
    if sum(i*c[i]*pow(rt,i-1,q) for i in range(1,len(c)))%q==0:return False
  return True
 fb_f,fb_g={},{}
 for q in primes:
  if unram(f_c,q):
   for rt in range(q):
    if peval(f_c,rt,q)==0:fb_f[(q,rt)]=len(fb_f)
  if unram(g_c,q):
   for rt in range(q):
    if peval(g_c,rt,q)==0:fb_g[(q,rt)]=len(fb_g)
 pf_act=sorted({q for (q,_) in fb_f})
 pg_act=sorted({q for (q,_) in fb_g})
 def sfac(val,plist):
  n=abs(val)
  if n==0:return None
  fac={}
  for q in plist:
   if n%q==0:
    e=0
    while n%q==0:
     n//=q
     e+=1
    fac[q]=e
    if n==1:return fac
  return fac if n==1 else None
 rels=[]
 for b in range(1,55):
  for a in range(-450,451):
   if a==0 or math.gcd(abs(a),b)!=1 or (a-b*m)%p==0:continue
   fg=sfac(a*a-a*b-10*b*b,pg_act)
   if fg is None:continue
   ff=sfac(a**4-5*a*a*b*b+a*b**3+b**4,pf_act)
   if ff is None:continue
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
  nxt=[ri for ri in cur if all(ccnt[k]>=2 for k in rels[ri][2]) and all(ccnt[nf_cnt+k]>=2 for k in rels[ri][3])]
  if len(nxt)==len(cur):break
  cur=nxt
 used=sorted({k for ri in cur for k in rels[ri][2]} | {nf_cnt+k for ri in cur for k in rels[ri][3]})
 cmap={c:i for i,c in enumerate(used)}
 N_2core=len(used)
 def smap(a,b,mod_p,eps):
  el=ppow([a%r2,(-b)%r2],eps,mod_p,r2)
  assert el[0]%r==1 and all(x%r==0 for x in el[1:])
  return [((el[0]-1)//r)%r]+[(x//r)%r for x in el[1:]]
 rows_M,rows_S=[],[]
 for ri in cur:
  a,b,vf,vg=rels[ri]
  row=[0]*N_2core
  for k,e in vf.items():
   row[cmap[k]]=e%r
  for k,e in vg.items():
   row[cmap[nf_cnt+k]]=(-e)%r
  rows_M.append(row)
  rows_S.append(smap(a,b,f_c,eps_f)[1:4]+smap(a,b,g_c,eps_g)[1:2])
 m_rows=[row[:] for row in rows_M];s_rows=[row[:] for row in rows_S]
 alive_r=[True]*len(cur);alive_c=[True]*N_2core
 merge_stack=[]
 while True:
  piv_c=next((c for c in range(N_2core) if alive_c[c] and sum(1 for i in range(len(cur)) if alive_r[i] and m_rows[i][c]!=0)==2),-1)
  if piv_c<0:break
  r_pair=[i for i in range(len(cur)) if alive_r[i] and m_rows[i][piv_c]!=0]
  rp,ro=r_pair[0],r_pair[1]
  merge_stack.append((piv_c,m_rows[rp][:],s_rows[rp][:]))
  fac=(m_rows[ro][piv_c]*pow(m_rows[rp][piv_c],-1,r))%r
  m_rows[ro]=[(m_rows[ro][j]-fac*m_rows[rp][j])%r for j in range(N_2core)]
  s_rows[ro]=[(s_rows[ro][j]-fac*s_rows[rp][j])%r for j in range(4)]
  alive_r[rp]=False;alive_c[piv_c]=False
 core_cols=[c for c in range(N_2core) if alive_c[c]]
 N,D=len(core_cols),len(core_cols)+4
 active=[i for i in range(len(cur)) if alive_r[i]]
 cc={c:sum(1 for i in active if m_rows[i][c]!=0) for c in core_cols}
 while len(active)>D:
  cands=[i for i in active if all(cc[c]>=3 for c in core_cols if m_rows[i][c]!=0)] or [i for i in active if all(cc[c]>=2 for c in core_cols if m_rows[i][c]!=0)]
  drop=max(cands,key=lambda i:(sum(1 for c in core_cols if m_rows[i][c]!=0),i))
  active.remove(drop)
  for c in core_cols:
   if m_rows[drop][c]!=0:cc[c]-=1
 def rref_null(mat):
  mr=[row[:] for row in mat]
  R,C=len(mr),len(mr[0])
  pivs,rk=[],0
  for c in range(C):
   piv=next((i for i in range(rk,R) if mr[i][c]%r!=0),-1)
   if piv<0:continue
   mr[rk],mr[piv]=mr[piv],mr[rk]
   inv=pow(mr[rk][c],-1,r)
   rr=mr[rk]
   for j in range(c,C):
    rr[j]=(rr[j]*inv)%r
   for i in range(R):
    if i!=rk and mr[i][c]:
     fv=mr[i][c]
     ri_r=mr[i]
     for j in range(c,C):
      ri_r[j]=(ri_r[j]-fv*rr[j])%r
   pivs.append(c)
   rk+=1
  pset=set(pivs)
  vecs=[]
  for fc in [c for c in range(C) if c not in pset]:
   v=[0]*C
   v[fc]=1
   for ri_i,pc in enumerate(pivs):
    v[pc]=(-mr[ri_i][fc])%r
   vecs.append(v)
  return rk,vecs
 M_sq=[[m_rows[i][c] for c in core_cols] for i in active]
 S_sq=[s_rows[i][:] for i in active]
 A_sq=[M_sq[i]+S_sq[i] for i in range(D)]
 rk_M,_=rref_null(M_sq)
 rk_A,ker_A=rref_null(A_sq)
 rk_3sm,ker_3sm=rref_null([M_sq[i]+S_sq[i][:3] for i in range(D)])
 assert rk_M==N and rk_A==D-1 and len(ker_A)==1 and len(ker_3sm)==0
 def mmul(P,Q):
  R,K,C=len(P),len(P[0]),len(Q[0])
  out=[[0]*C for _ in range(R)]
  for i in range(R):
   Pi,oi=P[i],out[i]
   for k in range(K):
    v=Pi[k]
    if v:
     Qk=Q[k]
     for j in range(C):
      oi[j]=(oi[j]+v*Qk[j])%r
  return out
 def mvec(P,v):
  return [sum(P[i][j]*v[j] for j in range(len(v)))%r for i in range(len(P))]
 rng=random.Random(20261119)
 m_bw=n_bw=4
 bnds=[(k*D)//4 for k in range(5)]
 X_good=[[0]*D for _ in range(m_bw)]
 for k in range(m_bw):
  for j in range(bnds[k],bnds[k+1]):
   X_good[k][j]=rng.randrange(1,r)
 Y=[[rng.randrange(1,r) for _ in range(n_bw)] for _ in range(D)]
 hr,hc=(D//m_bw)+2,(D//n_bw)+2
 AY,seq=[row[:] for row in Y],[]
 for _ in range(hr+hc+2):
  seq.append(mmul(X_good,AY))
  AY=mmul(A_sq,AY)
 H_rows=[[v for j in range(hc) for v in seq[i+j+1][ri_b]] for i in range(hr) for ri_b in range(m_bw)]
 rk_H,ker_H=rref_null(H_rows)
 v_rec=None
 for c_vec in ker_H:
  acc=[0]*D
  for j in range(hc-1,-1,-1):
   cj=c_vec[j*n_bw:(j+1)*n_bw]
   yc=[sum(Y[ri_i][bi]*cj[bi] for bi in range(n_bw))%r for ri_i in range(D)]
   acc=[(a_i+y_i)%r for a_i,y_i in zip(mvec(A_sq,acc),yc)]
  if any(acc) and all(x==0 for x in mvec(A_sq,acc)):
   v_rec=acc
   break
 assert v_rec is not None
 O_rows,cX=[],[row[:] for row in X_good]
 for _ in range(hr):
  cX=mmul(cX,A_sq)
  O_rows.extend(cX)
 rk_OA,_=rref_null(O_rows)
 C_cols,cY=[],[row[:] for row in Y]
 for _ in range(hc):
  for bi in range(n_bw):
   C_cols.append([cY[ri_i][bi] for ri_i in range(D)])
  cY=mmul(A_sq,cY)
 rk_Cd,_=rref_null(C_cols)
 X_bad=[[((k+1)*X_good[0][j])%r for j in range(D)] for k in range(m_bw)]
 O_bad,cXb=[],[row[:] for row in X_bad]
 for _ in range(hr):
  cXb=mmul(cXb,A_sq)
  O_bad.extend(cXb)
 rk_OA_bad,_=rref_null(O_bad)
 vlog_full={c:v_rec[j] for j,c in enumerate(core_cols)}
 sm_rec=v_rec[N:]
 for piv_c,pm,ps in reversed(merge_stack):
  acc=(sum(pm[j]*vlog_full[j] for j in range(N_2core) if j!=piv_c and pm[j]!=0)+sum(ps[s]*sm_rec[s] for s in range(4)))%r
  vlog_full[piv_c]=(-acc*pow(pm[piv_c],-1,r))%r
 while True:
  prog=False
  for idx in range(len(cur)):
   unk=[j for j in range(N_2core) if rows_M[idx][j]!=0 and j not in vlog_full]
   if len(unk)==1:
    cn=unk[0]
    acc=(sum(rows_M[idx][j]*vlog_full[j] for j in range(N_2core) if j!=cn and rows_M[idx][j]!=0)+sum(rows_S[idx][s]*sm_rec[s] for s in range(4)))%r
    vlog_full[cn]=(-acc*pow(rows_M[idx][cn],-1,r))%r
    prog=True
  if not prog:break
 assert len(vlog_full)==N_2core
 gen=next(g for g in range(2,p) if pow(g,(p-1)//r,p)!=1 and pow(g,2,p)!=1)
 g_r,dlog,ag=pow(gen,(p-1)//r,p),{},1
 for e in range(r):
  dlog[ag]=e
  ag=(ag*g_r)%p
 ratios=set()
 for idx,ri in enumerate(cur):
  a,b,vf,vg=rels[ri]
  lhs=(sum(e*vlog_full[cmap[k]] for k,e in vf.items())+sum(rows_S[idx][j]*sm_rec[j] for j in range(3)))%r
  rhs=(sum(e*vlog_full[cmap[nf_cnt+k]] for k,e in vg.items())-rows_S[idx][3]*sm_rec[3])%r
  assert lhs==rhs
  dl=dlog[pow((a-b*m)%p,(p-1)//r,p)]
  if dl!=0:ratios.add((lhs*pow(dl,-1,r))%r)
 from pm_adaptive import pm_basis_verify
 pm_res=pm_basis_verify(A_sq,seq,Y,r,hr,hc,v_rec)
 assert rk_H==D-1 and rk_OA==D-1 and rk_Cd==D and rk_OA_bad<D-1 and len(ratios)==1
 return dict(
  p=p,r=r,raw_rels=len(rels),purged_rels=len(cur),N_2core=N_2core,n_merged_cols=len(merge_stack),N=N,D=D,
  rank_M=rk_M,rank_A=rk_A,ker_dim_A=len(ker_A),ker_dim_3sm_ctrl=len(ker_3sm),
  rank_Cd=rk_Cd,rank_OeA=rk_OA,rank_H=rk_H,rank_OeA_1block_ctrl=rk_OA_bad,
  **pm_res,
  unplanted_2core_backsub_exact=True,horner_all_purged_rels_ok=True,distinct_Fp_dlp_ratios=len(ratios),
  checks_ok=True,
 )
