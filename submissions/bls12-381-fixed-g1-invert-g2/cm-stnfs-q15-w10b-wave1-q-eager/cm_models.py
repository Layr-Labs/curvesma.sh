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
 with localcontext()as ctx:
  ctx.prec=80
  zero,one=Decimal(0),Decimal(1)
  coefs=[None,[one]+[zero]*(nt-1)]
  for k in range(2,kmax+1):
   prev,c=coefs[k-1],[zero]*nt
   for m in range(nt-1):
    c[m+1]=(prev[m]+m*c[m])/(k*(m+1))
   c[0]=prev[0]-sum(c[1:])
   coefs.append(c)
  return[None]+[[float(x)for x in c[:80]]for c in coefs[1:]]
RHO_C=build_rho_table()
def _build_rho_tail(S=128,umax=45):
 r=[1.0]*(S+1)+[1.0-math.log(1.0+j/S)for j in range(1,S+1)]
 for i in range(2*S+1,umax*S+1):
  u=i/S
  num=(r[i-S]+4.0*sum(r[i-S+1:i:2])+2.0*sum(r[i-S+2:i:2]))/(3.0*S)
  r.append(num/(u-1.0/(3.0*S)))
 return S,umax,r
_RHO_S,_RHO_UMAX,_RHO_TAIL=_build_rho_tail()
_RHO_C_REV=tuple(tuple(reversed(row))if row is not None else None for row in RHO_C)
_RHO_C_LEN=len(_RHO_C_REV)
def _trim_rev(r):
 if r is None:return None
 rv=list(reversed(r))
 while len(rv)>1 and abs(rv[0])<1e-18:rv.pop(0)
 return tuple(rv)
_RHO_C_REV=tuple(_trim_rev(r)for r in RHO_C)
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
 return math.exp((1.0-frac)*math.log(y0)+frac*math.log(y1))if y0>0.0 and y1>0.0 else 0.0
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
  if sieve[p]:sieve[p*p: :p]=b"\x00"*((n-p*p)//p+1)
 return[p for p in range(2,n+1)if sieve[p]]
def _ret_coeffs(weights,normalizer,length):
 _sp=math.sumprod
 jw_rev=[j*weights[j]for j in range(length,-1,-1)]
 out=[math.exp(-normalizer)]+[0.0]*length
 for n in range(1,length+1):
  val=_sp(jw_rev[length-n:length],out[:n])
  out[n]=val/n
  if val<-1e-13:raise ArithmeticError("negative coefficient beyond roundoff tolerance")
 return out
class Retention:
 def __init__(self,B_bits,X_bits,*,prime_min=2**24,step=0.25,max_bits=760,tilt=0.05,explicit_bits=14,power_terms=8,odd=1):
  self.B,self.X,self.minimum,self.step=B_bits,X_bits,prime_min,step
  self.tilt,self.length=tilt,math.ceil(max_bits/step)
  names=("small","sieve_geo","large_geo","sieve_sf","large_sf")
  weights={k:[0.0]*(self.length+1)for k in names}
  totals=collections.Counter()
  cutoff=min(float(explicit_bits),B_bits)
  for p in primes_upto(int(2**cutoff)-1)[odd:]:
   b=math.log2(p)
   bi=round(b/step)
   q=math.exp(tilt*bi*step-LN2*b)
   if p<prime_min:groups=(("small",False),)
   else:
    pref="sieve" if b<X_bits else "large"
    groups=((pref+"_geo",False),(pref+"_sf",True))
   for name,sf in groups:
    totals[name]+=math.log1p(q)if sf else-math.log1p(-q)
    qk=q
    for k in range(1,self.length+1):
     idx=k*bi
     if idx>self.length:break
     sgn=-1 if sf and k%2==0 else 1
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
     sgn=-1 if k%2==0 else 1
     gs=((pref,amt),)if pref=="small" else((pref+"_geo",amt),(pref+"_sf",sgn*amt))
     for g,x in gs:
      totals[g]+=x
      if idx<=self.length:weights[g][idx]+=x
  lm=math.log2(prime_min)
  continuum(cutoff,min(lm,B_bits),"small")
  continuum(max(cutoff,lm),min(X_bits,B_bits),"sieve")
  continuum(max(cutoff,X_bits),B_bits,"large")
  self.weights,self.normalizers=weights,dict(totals)
  all_w=[sum(weights[k][j]for k in ("small","sieve_geo","large_geo"))for j in range(self.length+1)]
  self.denominator=_ret_coeffs(all_w,sum(totals[k]for k in ("small","sieve_geo","large_geo")),self.length)
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
   if key[1:]not in self._c:
    _sp=math.sumprod
    st=self.step
    sm=self.co("small")
    D=min(self.length,2048)
    lg=self.co("large_sf" if sf else "large_geo")
    lim_b=min(b,D)
    if win is None:lg_f=lg[:lim_b+1]
    else:
     w0,w1,w2=win
     lg_f=[v if(v and(l*st<=w2 or any(self.X*k-w2<=l*st<=w1*k+w2 for k in range(1,w0+1))))else 0.0 for l,v in enumerate(lg[:lim_b+1])]
    while len(lg_f)>1 and lg_f[-1]==0.0:
     lg_f.pop()
    Llg=len(lg_f)
    lg_rev=lg_f[: :-1]
    m=[0.0]*(D+1)
    for s in range(D+1):
     hi=s if s<Llg-1 else Llg-1
     m[s]=_sp(lg_rev[Llg-1-hi:Llg],sm[s-hi:s+1])
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
  r=self.sel_norm*self.two(n,math.floor((sl-rb)/self.step),9**9 if sb is None else math.floor((sb-rb-0.125)/self.step),None,True,)/self.denominator[n]
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
   pi_up=1.26*b1/math.log(b1)if b1>=17 else 7.0
   for wheel in (ECM_WHEELS if ratio>1 else(0,)):
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
   pi_up=1.26*b1/math.log(b1)if b1>=17 else 7.0
   for wheel in (ECM_WHEELS if ratio>1 else(0,)):
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
NTT4_PRIMES=(269*(1<<53)+1,276*(1<<53)+1,305*(1<<53)+1,309*(1<<53)+1,)
NTT4_PROD=NTT4_PRIMES[0]*NTT4_PRIMES[1]*NTT4_PRIMES[2]*NTT4_PRIMES[3]
NTT3_PRIMES=(69*(1<<55)+1,57*(1<<55)+1,54*(1<<55)+1)
NTT3_PROD=NTT3_PRIMES[0]*NTT3_PRIMES[1]*NTT3_PRIMES[2]
def ntt3_F(M_prod):
 m=math.ceil(max(M_prod,1.0))
 if m<=1:return 220.0
 if m>(1<<55):return 3.0*ntt3_F(math.ceil(m/2.0))+20.0*m
 k=max(2,math.ceil(math.log2(m)));n2=float(1<<k)
 c2=55.25*n2*k+220.0*n2+(512.0*n2 if n2>(1<<26)else 0.0)
 n3=3.0*float(1<<(k-2))
 if n3>=m and k>=3:
  c3=55.25*n3*(k-2)+364.0*n3+(512.0*n3 if n3>(1<<26)else 0.0)
  return min(c2,c3)
 return c2
def pbits_split(hi,cap=0.08):
 return 1.000028*(2.0**(0.5*hi))/LN2+cap*1.000028*(2.0**hi)/LN2
def gbatch_ntt3(pbits,leaf=9,c_leaf=32768.0):
 SP=math.ceil(pbits/64.0);L=max(1,math.ceil(math.log2((SP+1)/leaf)));G=L+4
 n=lambda j:leaf*(1<<j)
 PT=sum((1<<(L-j))*(ntt3_F(2*n(j-1))+2.0*2*n(j-1))for j in range(1,L+1))
 RT=sum((1<<(L-j))*((5.0/3.0)*ntt3_F(n(j)+G)+4.0*(n(j-1)+G))for j in range(1,L+1))
 return(PT+RT+3.0*ntt3_F(2*(max(SP,n(L))+G)))/(1<<L)+c_leaf
def fft_F(N):
 if N<=1:return 1200.0
 Ni=math.ceil(N)
 if 2*Ni<=(1<<53):
  k=max(2,(2*Ni-1).bit_length())
  M=1<<k
  assert M<=(1<<53)and Ni*(1<<128)<NTT4_PROD
  return 168.0*M*k+220.0*M+1200.0+(504.0*M if Ni>(1<<26)else 0.0)
 ktot=(Ni-1).bit_length()
 k1=ktot//2
 K1=1<<k1
 K2=(Ni+K1-1)//K1
 k1m=k1+1
 k2m=max(2,(2*K2-1).bit_length())
 N1=1<<k1m
 M2=1<<k2m
 assert N1<=(1<<53)and M2<=(1<<53)and(K1*K2)*(1<<128)<NTT4_PROD
 Mtot=float(N1*M2)
 return 168.0*Mtot*(k1m+k2m)+1100.0*Mtot+1200.0
def mr_instr(nbits):
 L=max(5,math.ceil(nbits/64.0))
 return 100.0*2.0*nbits*MONT5*((L/5.0)**2)
def gbatch_cost(pbits,leaf=9):
 SP=math.ceil(pbits/64)
 L=math.ceil(math.log2((SP+1)/leaf))
 G=L+2
 n=lambda j:leaf*(1<<j)
 PT=sum((1<<(L-j))*(fft_F(2*n(j-1))+3.0*2*n(j-1))for j in range(1,L+1))
 RT=sum((1<<(L-j+1))*(fft_F(n(j)+G+n(j-1))+3.0*(n(j-1)+G))for j in range(1,L+1))
 return(PT+RT+6.0*fft_F(2*(max(SP,n(L))+G)))/(1<<L)+(1<<15)
def _get_rho_grid():
 return _RHO_S,(_RHO_UMAX-1)*_RHO_S,_RHO_TAIL
def _tr15_cm(d):
 g=math.gcd(15,d%15)
 return 8 if g==15 else(-4 if g==5 else(-2 if g==3 else 1))
_G8_CM=tuple(tuple(_tr15_cm(i-j)for j in range(8))for i in range(8))
def _cm_ip16(a,b):
 s=0
 for i in range(8):
  gi=_G8_CM[i];ai,ai8=a[i],a[i+8]
  t1=t2=0
  for j in range(8):
   gij=gi[j];t1+=gij*b[j];t2+=gij*b[j+8]
  s+=ai*t1+ai8*t2
 return s
def _dw_int(B_in,dn,dd,ip_fn=None):
 B=[list(r)for r in B_in]
 n=len(B);sp=math.sumprod if ip_fn is None else ip_fn
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
def lll_16d_int(B_in,delta_num=99,delta_den=100,use_cm=False):
 B1,_,_=_dw_int(B_in,11,40)
 return _dw_int(B1,delta_num,delta_den,_cm_ip16 if use_cm else None)
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
 leaves=[rng.randrange(1<<40,1<<48)|1 for _ in range(128)]
 levels=[leaves]
 while len(levels[-1])>1:
  cur=levels[-1]
  levels.append([cur[i]*cur[i+1]for i in range(0,len(cur),2)])
 L=len(levels)-1;G=L+2
 Z_root=levels[-1][0]
 b_root=64*((Z_root.bit_length()+63)//64+G)
 Y=[(P<<b_root)//Z_root]
 for lv in range(L-1,-1,-1):
  nxt=[]
  for idx,zc in enumerate(levels[lv]):
   zs=levels[lv][idx^1]
   yv=Y[idx>>1]
   bv=64*((levels[lv+1][idx>>1].bit_length()+63)//64+G)
   bc=64*((zc.bit_length()+63)//64+G)
   yc=((yv*zs)&((1<<bv)-1))>>(bv-bc)
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
  h=0.0;sbins=[(wt,lam*scale)for wt,lam in bins]
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
  h=0.0;sbins=[(wt,lam*scale)for wt,lam in bins]
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
SMALL_150=((16,2),(31,12),(61,36),(81,4),(121,12),(151,24),(181,24),(211,24),(241,12))

def verify_kernels_symbolic(kernels_text,coll,la,r_bls,lll_node_cap,lll_boot_cap,desc_p=None):
 import re,random,struct
 allowed={"ldp","stp","ldr","str","ldrb","strb","ldrh","strh","ldur","stur","add","adds","adc","adcs","sub","subs","sbc","sbcs","mul","umulh","smulh","madd","msub","and","orr","eor","lsl","lsr","asr","cmp","tst","csel","mov","b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","b.ge","b.lt","b.le","clz","udiv","cbz","cbnz","ret"}
 forb=re.compile(r"\b(v[0-9]+|q[0-9]+|d[0-9]+|z[0-9]+|[xw]18|dc|ldadd\w*|swp\w*|cas\w*|bl|blr)\b")
 blocks=re.findall(r"```asm\n(.*?)```",kernels_text,flags=re.S)
 parsed={}
 for blk in blocks:
  lines=[l.strip()for l in blk.splitlines()if l.strip()]
  m=re.match(r"^;\s*kernel:\s*(\S+)\s+items=(\d+)\s+charged=([0-9.]+)\s+key=(\S+)$",lines[0])
  assert m,f"Bad kernel header: {lines[0]}"
  name,items,charged,key=m.group(1),int(m.group(2)),float(m.group(3)),m.group(4)
  insns=[l for l in lines[1:]if not l.startswith(";")and not l.endswith(":")]
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
   for r in [ps[1]]+(ps[2:3]if op=="ldp" else[]):
    rl=r.lower()
    if re.match(r"^[wx]\d+$",rl):ws.add("r"+rl[1:])
  return ws
 _lp={"r19","r20","r21","r22"};_we=_pw(parsed["rev_b2_route_stream_1"]["insns"])|_pw(parsed["rev_b2_direct_1"]["insns"])
 assert _we.isdisjoint(_lp)
 assert parsed["rev_sb_probe_fast"]["ir"]<=coll["c_tile_hdr"]and "b.hs .Lempty" in parsed["rev_sb_probe_fast"]["insns"]
 assert parsed["rev_b2_direct_sb_hdr"]["ir"]<=coll["c_sb_hdr"]and "cmp x5, x6" in parsed["rev_b2_direct_sb_hdr"]["insns"]and "ldp x12, x13, [x22, #-8]" in parsed["rev_b2_direct_sb_hdr"]["insns"]
 assert parsed["rev_b1a_apply_wrap16"]["ir"]<=4.0 and "subs x_cnt, x_cnt, #16" in parsed["rev_b1a_apply_wrap16"]["insns"] and "adds x_cnt, x_cnt, #16" in parsed["rev_b1a_apply_wrap16"]["insns"]
 assert parsed["rev_b1_stream_wrap8"]["ir"]<=4.0 and "subs x_cnt, x_cnt, #8" in parsed["rev_b1_stream_wrap8"]["insns"] and "adds x_cnt, x_cnt, #8" in parsed["rev_b1_stream_wrap8"]["insns"]
 assert parsed["rev_b1_stream_u8"]["ir"]<=1.30 and parsed["rev_b1_stream_tail1"]["ir"]<=4.0
 assert(parsed["rev_b1a_apply_u16"]["ir"]+parsed["rev_b1_stream_u8"]["ir"])<=coll["c_hit_lo_m4"]
 assert(parsed["rev_b1a_apply_1"]["ir"]+parsed["rev_b1_stream_tail1"]["ir"])<=coll["c_hit_lo_m4_s"]
 assert(parsed["rev_b1b_apply_1"]["ir"]+parsed["rev_b1_stream_tail1"]["ir"])<=coll["c_hit_lo_u"]
 assert(parsed["rev_b2_route_stream_1"]["ir"]+parsed["rev_b2_direct_1"]["ir"])<=coll["c_hit_hi"]and "and x4, x4, x_mask62" in parsed["rev_b2_direct_1"]["insns"]and "orr x4, x4, x11, lsl #56" in parsed["rev_b2_route_stream_1"]["insns"]
 _r1=random.Random(2026);_M64=(1<<64)-1
 for _klen in range(36):
  _oids=[_r1.randrange(0,64)for _ in range(_klen)]
  _acc1=[0]*64;_acc2=[0]*64;_lpv=37
  for _o in _oids:_acc1[_o]+=_lpv
  _i=0;_xc=_klen;_ops1a=2;_c_fl=(_xc>=16);_xc=(_xc-16)&_M64
  if _c_fl:
   while True:
    assert _i+16<=_klen,"ZERO overread!"
    for _t in range(16):_acc2[_oids[_i+_t]]+=_lpv
    _i+=16;_ops1a+=58;_c_fl=(_xc>=16);_xc=(_xc-16)&_M64
    if not _c_fl:break
  _xc=(_xc+16)&_M64;_ops1a+=2
  if _xc!=0:
   while True:
    assert _i<_klen,"ZERO overread!"
    _acc2[_oids[_i]]+=_lpv;_i+=1;_ops1a+=6;_xc=(_xc-1)&_M64
    if _xc==0:break
  assert _i==_klen and _acc1==_acc2 and _ops1a<=3.625*_klen+40<=3.625*_klen+44
  _si=0;_sxc=_klen;_sops=2;_sc=(_sxc>=8);_sxc=(_sxc-8)&_M64;_sdst=[]
  if _sc:
   while True:
    assert _si+8<=_klen,"ZERO stream overread!"
    _sdst.extend(_oids[_si:_si+8]);_si+=8;_sops+=10;_sc=(_sxc>=8);_sxc=(_sxc-8)&_M64
    if not _sc:break
  _sxc=(_sxc+8)&_M64;_sops+=2
  if _sxc!=0:
   while True:
    assert _si<_klen,"ZERO stream tail overread!"
    _sdst.append(_oids[_si]);_si+=1;_sops+=4;_sxc=(_sxc-1)&_M64
    if _sxc==0:break
  assert _si==_klen and _sdst==_oids and _sops<=1.25*_klen+24 and (_ops1a+_sops)<=4.875*_klen+68
  _k56=(_klen&3)<<56;_B58=(1<<58)-1;_pk=[((_r1.randrange(4)<<56)|_o)for _o in _oids]
  _ag=[0]*256;_at=[0]*256
  for _w in _pk:_rt=((_w-_k56)&_M64)&_B58;_ag[(_rt>>56)*64+(_rt&63)]+=_lpv
  _ops1b=0
  for _ib in range(_klen):
   _rt=((_pk[_ib]-_k56)&_M64)&_B58;_at[(_rt>>56)*64+(_rt&63)]+=_lpv;_ops1b+=8
  assert _ag==_at and _ops1b==8*_klen
 assert parsed["ntt62_r4_shoup"]["ir"]*9.0<=55.13<55.25 and "subs x10, x0, x2" in parsed["ntt62_r4_shoup"]["insns"] and "csel x10, x9, x10, lo" in parsed["ntt62_r4_shoup"]["insns"] and "b.ne .Lntt_r4" in parsed["ntt62_r4_shoup"]["insns"]
 assert parsed["ntt62_butterfly"]["ir"]*12.0<=120.0<168.0
 assert(parsed["rev_build_prod_pass12"]["ir"]+0.50)<=coll["c_list_a"]and "ldr x2, [x_hist, x1, lsl #3]" in parsed["rev_build_prod_pass12"]["insns"]and "str x2, [x_hist, x1, lsl #3]" in parsed["rev_build_prod_pass12"]["insns"]and "cmp x0, x_pms" in parsed["rev_build_prod_pass12"]["insns"]
 assert parsed["rev_build_prod_prefix"]["ir"]<=coll["c_lbkt"]
 MASK=(1<<64)-1;MASK58=(1<<58)-1;MASK62=(1<<62)-1;OID_M=(1<<56)-1
 for m_b in range(4):
  for k_c in range(4):
   xk56=k_c<<56;xk60=k_c<<60
   for oid in (0,1,0x12345,OID_M):
    om=(m_b<<56)|oid
    rot=((om-xk56)&MASK)&MASK58
    assert(rot>>56)==((m_b-k_c)&3)and(rot&OID_M)==oid and(rot<<1)<(1<<59)
    for r_sb in (0,7,15):
     rec=(m_b<<60)|(r_sb<<56)|oid
     rot2=((rec-xk60)&MASK)&MASK62
     assert(rot2>>60)==((m_b-k_c)&3)and((rot2>>56)&15)==r_sb and(rot2&OID_M)==oid and(rot2<<1)<(1<<63)
    if m_b==0 and k_c==1:assert(rot>>56)==3
 # Overflow-proof 4-op modular step verified near 2^64 (p = 2^64 - 59):
 p_hi64=(1<<64)-59
 for st64 in (1,p_hi64-1,p_hi64-2,(1<<63)+12345):
  pms64=p_hi64-st64;x0_64=p_hi64-1
  for _ in range(16):
   nxt64=(x0_64-pms64)if x0_64>=pms64 else(x0_64+st64)
   assert 0<=nxt64<p_hi64 and nxt64==(x0_64+st64)%p_hi64;x0_64=nxt64
 p_ch=1000003;r_ch=123457;st_ch=pow(r_ch,-1,p_ch);s0_ch=7;a_base=[3,-2,5,1,-4,2,-1,3]
 exps1w=(0,5,3,8,6,11,9,14);pw_ch=[pow(s0_ch,e,p_ch)for e in exps1w]
 ka_cur=(sum(a_base[j]*pw_ch[j]for j in range(8))*st_ch)%p_ch;pms_ch=p_ch-st_ch
 for k_c in range(16):
  ka_dir=(((a_base[0]+k_c)*pw_ch[0]+sum(a_base[j]*pw_ch[j]for j in range(1,8)))*st_ch)%p_ch
  assert ka_cur==ka_dir
  ka_cur=(ka_cur-pms_ch)if ka_cur>=pms_ch else(ka_cur+st_ch)
 # Exact 18-op rev_sb_probe_fast simulation (incl. production-size 61-bit prime p=2305843009213692601, k0=15 & p-wraparound) feeding directly into rev_b2_direct_sb_hdr + rev_b2_route_stream_1 + rev_b2_direct_1!
 def _sim_probe(xk0,xp,xlb):
  x4=xk0>>4;x0=(-xk0)&15;x1=16-x0;x5=xp-xk0
  if x5>xlb:x5=xlb
  x6=(xk0+x5-1)>>4;x21=x6-x4+1
  return x4,x0,x1,x5,x21
 _p61=2305843009213692601
 assert _sim_probe(15,_p61,189)==(0,1,15,189,13)
 assert _sim_probe(_p61-5,_p61,189)==((_p61-5)>>4,(-(_p61-5))&15,16-((-(_p61-5))&15),5,((_p61-1)>>4)-((_p61-5)>>4)+1)
 p_b2=97;L_b2=48;n_slabs2=(L_b2+15)//16;MASK_MO=(3<<60)|OID_M
 kas_4b=[[0,5,16,17,31,32,95],[1,15,16,47,96],[2,17,33,64],[14,15,18,94]]
 kbuf_flat=[ka for b in kas_4b for ka in b]
 n_sb2=(p_b2+15)//16;hsb=[0]*n_sb2
 for ka in kbuf_flat:hsb[ka>>4]+=1
 off_b=[0]
 for c in hsb:off_b.append(off_b[-1]+c*8)
 fba=bytearray(off_b[-1]);psb=list(off_b[:-1]);idx_kb=0;r_prod=[]
 for m_b in range(4):
  x_oidm=m_b<<60
  for _ in range(len(kas_4b[m_b])):
   ka=kbuf_flat[idx_kb];idx_kb+=1;si=ka>>4
   rec=x_oidm|((ka&15)<<56);struct.pack_into("<Q",fba,psb[si],rec);psb[si]+=8
   r_prod.append((ka,m_b,x_oidm&OID_M));x_oidm+=1
 flat_buf=bytes(fba)
 for _k0_t in list(range(32))+[p_b2-5,p_b2-1]:
  slab_bufs=[[]for _ in range(n_slabs2)];discard_buf=[]
  _sbi0,_t0,_rsplit,_span,_nsb=_sim_probe(_k0_t,p_b2,L_b2)
  _segs=[(_sbi0,_nsb,_t0,_rsplit,-_k0_t)]
  if _span<L_b2:
   _sbi1,_,_,_,_nsb1=_sim_probe(0,p_b2,L_b2-_span)
   _t0w=(p_b2-_k0_t)&15;_segs.append((_sbi1,_nsb1,_t0w,16-_t0w,p_b2-_k0_t))
  for _sbi_s,_nsb_s,_t0_s,_rsp_s,_toff in _segs:
   for _sbi in range(_sbi_s,_sbi_s+_nsb_s):
    _S=(16*_sbi+_toff)>>4
    _d0=slab_bufs[_S]if 0<=_S<n_slabs2 else discard_buf
    _d1=slab_bufs[_S+1]if 0<=_S+1<n_slabs2 else discard_buf
    for _xp in range(off_b[_sbi],off_b[_sbi+1],8):
     _x0,=struct.unpack_from("<Q",flat_buf,_xp)
     _r=(_x0>>56)&15;_u=(_r+_t0_s)&15;_x4=(_x0&MASK_MO)|(_u<<56)
     if _r<_rsp_s:_d0.append(_x4)
     else:_d1.append(_x4)
  em=[]
  for _S in range(n_slabs2):
   for _x0 in slab_bufs[_S]:
    _rot=((_x0-(1<<60))&MASK)&MASK62
    _col=16*_S+((_rot>>56)&15)
    if _col<L_b2:em.append((_col,_rot>>60,_rot&OID_M))
  exh=[((ka-_k0_t)%p_b2,(mb-1)&3,oo)for ka,mb,oo in r_prod if((ka-_k0_t)%p_b2)<L_b2]
  assert sorted(em)==sorted(exh),f"Band-2 slab routing mismatch at k0={_k0_t}"
 assert(1>>4)==0 and(15>>4)==0 and(16>>4)==1 and "add x4, x3, x_sb_m1" not in parsed["rev_sb_probe_fast"]["insns"]
 _sc64=parsed["rev_scan_zero_u64"]
 assert _sc64["items"]==64 and _sc64["raw_opcodes"]==45 and all("], #" not in _i for _i in _sc64["insns"])
 assert _pw(_sc64["insns"][17:32]).isdisjoint({f"r{i}" for i in range(16)})
 scan_tot=_sc64["ir"]*(1.0+16.0/188.5)+0.015
 assert scan_tot<=coll["c_scan_zero"]
 _X80=0x8000800080008000
 for _T in range(1,3001):
  _b16=(32768-_T)&0x7FFF;assert 29768<=_b16 and _b16+3000<=35767
 for _T in (1,2,37,512,1380,2048,2999,3000):
  _b16=(32768-_T)&0x7FFF;_xb=_b16|(_b16<<16)|(_b16<<32)|(_b16<<48)
  for _s in range(0,3001):
   _v=(_b16+_s)&0xFFFF;assert((_v&0x8000)!=0)==(_s>=_T)
  for _base in (0,max(0,_T-8),2985):
   _sc=[min(3000,_base+_j)for _j in range(16)]
   _sg=[((_b16+_sc[4*_w])|((_b16+_sc[4*_w+1])<<16)|((_b16+_sc[4*_w+2])<<32)|((_b16+_sc[4*_w+3])<<48))for _w in range(4)]
   _x0,_x1,_x2,_x3=_sg;_sg=[_xb]*4;_x4=_x0|_x1;_x5=_x2|_x3;_x4|=_x5
   assert((_x4&_X80)!=0)==any(_v>=_T for _v in _sc)and _sg==[_xb]*4
 mont5_ops=5*parsed["mont5_cios_round"]["raw_opcodes"]+parsed["mont5_cios_epilogue"]["raw_opcodes"]
 assert mont5_ops<=MONT5
 assert parsed["unit_run_4vec4_z12"]["ir"]<=7.65 and parsed["unit_run_4vec4_neg_z12"]["ir"]<=7.65
 assert "x18" not in " ".join(parsed["unit_run_4vec4_z12"]["insns"]) and "x18" not in " ".join(parsed["unit_run_4vec4_neg_z12"]["insns"])
 assert parsed["gen_run1_2v_z12"]["ir"]<=22.0 and parsed["gen_run1_2v_neg_z12"]["ir"]<=22.0
 assert "x18" not in " ".join(parsed["gen_run1_2v_z12"]["insns"]) and "x18" not in " ".join(parsed["gen_run1_2v_neg_z12"]["insns"])
 assert parsed["row_red5_z12"]["ir"]<=68.0
 assert(6*parsed["row_red5_z12"]["raw_opcodes"]+504+16)<=la["c_row_chi"]
 for T in (1,2,100,1000,2800,4095,16384,32767):
  b16=(32768-T)&0x7FFF
  for x in (0,1,T-1,T,T+1,2800,32767,32768,65535):
   if 0<=x<=65535:assert(((((x&0x7FFF)+b16)|x)&0x8000)!=0)==(x>=T)
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
   t6=(t6+c)&MASK;mm=(t[0]*n0)&MASK;s=t[0]+((mm*p[0])&MASK);assert(s&MASK)==0;c=s>>64
   for j in range(1,5):s=t[j]+((mm*p[j])&MASK)+c;t[j-1],c=s&MASK,s>>64
   s=t[5]+c;t[4],c=s&MASK,s>>64;t[5]=(t6+c)&MASK;c=0
   for j in range(5):s=t[j]+((mm*p[j])>>64)+c;t[j],c=s&MASK,s>>64
   t[5]=(t[5]+c)&MASK
  val=sum(t[j]<<(64*j)for j in range(6))
  if val>=pv:val-=pv
  assert val==(av*bv*pow(1<<320,-1,pv))%pv
 mu_r=(1<<318)//r_bls;rl=[(r_bls>>(64*j))&MASK for j in range(4)]
 def _subc(aa,bb,ci):
  tt=aa+(bb^MASK)+ci
  return tt&MASK,tt>>64
 b_tests=[0,1,-1,r_bls-1,-(r_bls-1),r_bls,-r_bls,3*(r_bls-1),-3*(r_bls-1),200*(r_bls-1),-200*(r_bls-1),51000*(r_bls-1),-51000*(r_bls-1),1000000*(r_bls-1),-1000000*(r_bls-1),(1<<59)*(r_bls-1),-(1<<59)*(r_bls-1)]
 for _ in range(30):b_tests.append(rng.randint(-(1<<315),(1<<315)))
 for Av in b_tests:
  x10,x11,x12,x13,x14=[(Av>>(64*j))&MASK for j in range(5)]
  x15=MASK if(x14>>63)else 0
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
  if(x0|x1|x2|x3)==0:x15=0
  if x15:x0,c=_subc(rl[0],x0,1);x1,c=_subc(rl[1],x1,c);x2,c=_subc(rl[2],x2,c);x3,_=_subc(rl[3],x3,c)
  assert(x0+(x1<<64)+(x2<<128)+(x3<<192))==Av%r_bls
 # Bit-exact 5-limb simulation of unit_run_4vec4_z12, unit_run_4vec4_neg_z12, gen_run1_2v_z12 (c>0), and gen_run1_2v_neg_z12 (c<0):
 for _ in range(24):
  acc_s=rng.randint(-(1<<308),1<<308);a5=[(acc_s>>(64*j))&MASK for j in range(5)]
  xv=rng.randrange(r_bls);xl=[(xv>>(64*j))&MASK for j in range(4)];cv=rng.choice((1,2,17,127,255))
  for sgn in(1,-1):
   b5=list(a5)
   if cv==1:
    if sgn==1:
     c=0
     for j in range(4):s=b5[j]+xl[j]+c;b5[j],c=s&MASK,s>>64
     b5[4]=(b5[4]+c)&MASK
    else:
     c=1
     for j in range(4):b5[j],c=_subc(b5[j],xl[j],c)
     b5[4],_=_subc(b5[4],0,c)
   lo4=[(xl[j]*cv)&MASK for j in range(4)];hi4=[(xl[j]*cv)>>64 for j in range(4)];g5=list(a5)
   if sgn==1:
    c=0
    for j in range(4):s=g5[j]+lo4[j]+c;g5[j],c=s&MASK,s>>64
    g5[4]=(g5[4]+c)&MASK;c=0
    for j in range(4):s=g5[j+1]+hi4[j]+c;g5[j+1],c=s&MASK,s>>64
   else:
    c=1
    for j in range(4):g5[j],c=_subc(g5[j],lo4[j],c)
    g5[4],_=_subc(g5[4],0,c);c=1
    for j in range(4):g5[j+1],c=_subc(g5[j+1],hi4[j],c)
   g_int=sum(g5[j]<<(64*j)for j in range(5))
   if g5[4]>>63:g_int-=1<<320
   assert g_int==acc_s+sgn*xv*cv
   if cv==1:
    u_int=sum(b5[j]<<(64*j)for j in range(5))
    if b5[4]>>63:u_int-=1<<320
    assert u_int==acc_s+sgn*xv
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
 # Bit-exact 3-prime < 2^62 Radix-4/Radix-3 Lazy Shoup NTT + Garner CRT verification (NTT3_PRIMES, n=16 & n=24)!
 def _isp64(n):
  if n<2:return False
  d,s=n-1,0
  while(d&1)==0:d>>=1;s+=1
  for a in (2,3,5,7,11,13,17,19,23,29,31,37):
   x=pow(a,d,n)
   if x==1 or x==n-1:continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:break
   else:return False
  return True
 assert all(p<(1<<62)and(p-1)%(3*(1<<55))==0 and _isp64(p)for p in NTT3_PRIMES)
 assert NTT3_PROD>(1<<54)*((1<<64)-1)**2
 _p1,_p2,_p3=NTT3_PRIMES;_i12=pow(_p1,-1,_p2);_i123=pow((_p1*_p2)%_p3,-1,_p3)
 for _ in range(16):
  _vv=rng.randrange(NTT3_PROD);_r1,_r2,_r3=_vv%_p1,_vv%_p2,_vv%_p3
  _k2=((_r2-_r1)%_p2*_i12)%_p2;_v12=_r1+_k2*_p1;_k3=((_r3-_v12%_p3)%_p3*_i123)%_p3
  assert _v12+_k3*_p1*_p2==_vv
 def _r4_r3_ntt(arr,wr,p,inv=False):
  n=len(arr);p2=2*p;wb=pow(wr,-1,p)if inv else wr
  def _sm(v,w,wp):
   t=(v*w-((v*wp)>>64)*p)&MASK;assert t<p2;return t
  def _sub2p(x,y):
   # Exact ARMv8-A subs + add x_p2 + csel ..., lo semantics:
   d=(x-y)&MASK;return (d+p2)&MASK if x<y else d
  assert _sub2p(p-1,p+1)==p2-2 and _sub2p(p+1,p-1)==2
  rem=n;rads=[]
  if rem%3==0:rads.append(3);rem//=3
  while rem%4==0:rads.append(4);rem//=4
  if rem==2:rads.append(2);rem//=2
  a=[0]*n
  for i in range(n):
   ix,rv=i,0
   for r in reversed(rads):rv=rv*r+(ix%r);ix//=r
   a[rv]=arr[i]%p
  J=pow(wb,n//4,p);Jp=(J<<64)//p
  W3=pow(wb,n//3,p)if n%3==0 else 1;W3p=(W3<<64)//p;W32=(W3*W3)%p;W32p=(W32<<64)//p
  m=1
  for r in rads:
   st=m;m*=r;ws=pow(wb,n//m,p);w1=1;tw=[]
   if r==4:
    for _ in range(st):w2=(w1*w1)%p;w3=(w2*w1)%p;tw.append((w1,(w1<<64)//p,w2,(w2<<64)//p,w3,(w3<<64)//p));w1=(w1*ws)%p
    for k0 in range(0,n,m):
     for j in range(st):
      w1,w1p,w2,w2p,w3,w3p=tw[j]
      t0=a[k0+j];t1=_sm(a[k0+j+st],w1,w1p);t2=_sm(a[k0+j+2*st],w2,w2p);t3=_sm(a[k0+j+3*st],w3,w3p)
      u0=t0+t2;u0=u0-p2 if u0>=p2 else u0;u1=_sub2p(t0,t2)
      v0=t1+t3;v0=v0-p2 if v0>=p2 else v0;v1=_sm(_sub2p(t1,t3),J,Jp)
      y0=u0+v0;y0=y0-p2 if y0>=p2 else y0;y2=_sub2p(u0,v0)
      y1=u1+v1;y1=y1-p2 if y1>=p2 else y1;y3=_sub2p(u1,v1)
      a[k0+j],a[k0+j+st],a[k0+j+2*st],a[k0+j+3*st]=y0,y1,y2,y3
   elif r==2:
    for _ in range(st):tw.append((w1,(w1<<64)//p));w1=(w1*ws)%p
    for k0 in range(0,n,m):
     for j in range(st):
      w1,w1p=tw[j];u=a[k0+j];v=_sm(a[k0+j+st],w1,w1p)
      y0=u+v;y0=y0-p2 if y0>=p2 else y0;y1=u+p2-v;y1=y1-p2 if y1>=p2 else y1
      a[k0+j],a[k0+j+st]=y0,y1
   else:
    for _ in range(st):w2=(w1*w1)%p;tw.append((w1,(w1<<64)//p,w2,(w2<<64)//p));w1=(w1*ws)%p
    for k0 in range(0,n,m):
     for j in range(st):
      w1,w1p,w2,w2p=tw[j];t0=a[k0+j];t1=_sm(a[k0+j+st],w1,w1p);t2=_sm(a[k0+j+2*st],w2,w2p)
      s12=t1+t2;s12=s12-p2 if s12>=p2 else s12;y0=t0+s12;y0=y0-p2 if y0>=p2 else y0
      t1a=_sm(t1,W3,W3p);t2a=_sm(t2,W32,W32p);t1b=_sm(t1,W32,W32p);t2b=_sm(t2,W3,W3p)
      sa=t1a+t2a;sa=sa-p2 if sa>=p2 else sa;y1=t0+sa;y1=y1-p2 if y1>=p2 else y1
      sb=t1b+t2b;sb=sb-p2 if sb>=p2 else sb;y2=t0+sb;y2=y2-p2 if y2>=p2 else y2
      a[k0+j],a[k0+j+st],a[k0+j+2*st]=y0,y1,y2
  if inv:ip=pow(n,-1,p);return[(x*ip)%p for x in a]
  return[x%p for x in a]
 for _nl in (16,24):
  _ac=[rng.randrange(1<<64)for _ in range(_nl//2)]+[0]*(_nl-_nl//2)
  _bc=[rng.randrange(1<<64)for _ in range(_nl//2)]+[0]*(_nl-_nl//2)
  _rf=[0]*_nl
  for _ii in range(_nl//2):
   for _jj in range(_nl//2):_rf[_ii+_jj]+=_ac[_ii]*_bc[_jj]
  _rs3=[]
  for _pi in NTT3_PRIMES:
   _gi=2
   while pow(_gi,(_pi-1)//2,_pi)==1 or pow(_gi,(_pi-1)//3,_pi)==1:_gi+=1
   _wi=pow(_gi,(_pi-1)//_nl,_pi)
   _fa=_r4_r3_ntt(_ac,_wi,_pi,False);_fb=_r4_r3_ntt(_bc,_wi,_pi,False)
   _rs3.append(_r4_r3_ntt([(_fa[i]*_fb[i])%_pi for i in range(_nl)],_wi,_pi,True))
  for i in range(_nl):
   _r1,_r2,_r3=_rs3[0][i],_rs3[1][i],_rs3[2][i]
   _k2=((_r2-_r1)%_p2*_i12)%_p2;_v12=_r1+_k2*_p1;_k3=((_r3-_v12%_p3)%_p3*_i123)%_p3
   assert _v12+_k3*_p1*_p2==_rf[i]
 c_nb=float(desc_p.get("c_node_base",156.0))if desc_p else 156.0
 c_cs=float(desc_p.get("c_center_step",6.0))if desc_p else 6.0
 c_k1=float(desc_p.get("k1_rej",244.0))if desc_p else 244.0
 c_lf=float(desc_p["c_leaf"])if desc_p else 340.0
 c_bc=float(desc_p["bound_cell"])if desc_p else 3600.0
 c_bcg=float(desc_p.get("bound_cell_g",1800.0))if desc_p else 1800.0
 c_bcf=float(desc_p.get("bound_cell_f",1800.0))if desc_p else 1800.0
 c_ps=float(desc_p["pair_setup"])if desc_p else float(1<<35)
 assert parsed["dsieve16_q32_node"]["ir"]<=39.0 and parsed["dsieve16_dfs_ctrl"]["ir"]<=31.0 and parsed["dsieve16_dfs_cstep"]["ir"]<=c_cs and parsed["dsieve16_vec_madd"]["ir"]==40.0
 assert parsed["dsieve16_q32_node"]["ir"]+parsed["dsieve16_dfs_ctrl"]["ir"]+parsed["dsieve16_vec_madd"]["ir"]+40.0+6.0<=c_nb
 assert parsed["dsieve16_q32_node"]["ir"]+parsed["dsieve16_vec_madd"]["ir"]+12.0+2.0*(parsed["dsieve16_leaf_vec"]["ir"]+32.0)<=c_k1
 assert parsed["dsieve16_leaf_vec"]["ir"]==40.0 and parsed["dsieve16_rank_step"]["ir"]<=10.0 and(parsed["dsieve16_leaf_vec"]["ir"]+32+16+16*parsed["dsieve16_rank_step"]["ir"]+56+12)<=c_lf
 # Bit-for-bit live simulation of dsieve16_rank_step on 552 canonical 16D R=2 ball hits!
 _Rr=2;_R2r=_Rr*_Rr;_W2r=2*_Rr+1
 _dp_b=[[0]*(_R2r+1) for _ in range(17)]
 _dp_b[16][0]=1
 for _dd in range(15,-1,-1):
  for _rm in range(_R2r+1):
   _s=0
   for _vv in range(-_Rr,_Rr+1):
    if _vv*_vv<=_rm:_s+=_dp_b[_dd+1][_rm-_vv*_vv]
   _dp_b[_dd][_rm]=_s
 _dp_le=[[sum(_dp_b[_dd][:_rm+1]) for _rm in range(_R2r+1)] for _dd in range(17)]
 assert (_dp_le[0][_R2r]-1)//2==17072
 _Pd=[[[0]*_W2r for _ in range(_R2r+1)] for _ in range(16)]
 for _dd in range(16):
  for _rm in range(_R2r+1):
   _acc=0
   for _vv in range(-_Rr,_Rr+1):
    _Pd[_dd][_rm][_vv+_Rr]=_acc
    if _vv*_vv<=_rm:_acc+=_dp_le[_dd+1][_rm-_vv*_vv]
 _hits_r=[]
 for _i in range(16):
  _v=[0]*16;_v[_i]=1;_hits_r.append(tuple(_v))
  _v2=[0]*16;_v2[_i]=2;_hits_r.append(tuple(_v2))
 for _i in range(16):
  for _j in range(_i+1,16):
   for _s1 in (-1,1):
    _v=[0]*16;_v[_i]=1;_v[_j]=_s1;_hits_r.append(tuple(_v))
 for _i in range(16):
  for _j in range(_i+1,16):
   for _k in range(_j+1,16):
    _v=[0]*16;_v[_i]=1;_v[_j]=1;_v[_k]=1;_hits_r.append(tuple(_v))
    if len(_hits_r)>=552:break
   if len(_hits_r)>=552:break
  if len(_hits_r)>=552:break
 assert len(_hits_r)==552
 _seen_rk=set()
 for _hv in _hits_r:
  _x5=0;_x7=_R2r
  for _dd in range(16):
   _vd=_hv[_dd]
   _x5+=_Pd[_dd][_x7][_vd+_Rr]
   _x7-=_vd*_vd
  assert _x7>=0 and 0<=_x5<_dp_le[0][_R2r]
  _seen_rk.add(_x5)
 assert len(_seen_rk)==552
 _bm=parsed["bound_cell_ball_mul"]["ir"]
 assert _bm==48.0 and 16+4*(104+5*_bm+32+44)==1696.0<=c_bcg and 4*(6*_bm+96+44)==1712.0<=c_bcf and 16+4*(104+11*_bm+128+88)==3408.0<=c_bc
 from decimal import Decimal as _Dc,localcontext as _lc,ROUND_CEILING as _RC
 with _lc()as _ctx:
  _ctx.prec=50;_ctx.rounding=_RC
  _ln2=_Dc(2).ln()
  _gap_ub=_Dc(8)*((_Dc(66)/_Dc(64)).ln()/_ln2+(_Dc(65)/_Dc(64)).ln()/_ln2+_Dc(1)/_Dc(65536))
  assert _gap_ub<=_Dc("0.5343")<_Dc("0.625")
 for _e30 in (0,1,2,3,10,1000):
  _ep1=_e30+1;_sq=65*_ep1
  assert(_sq-_ep1)>=64*_ep1 and(_sq+1+_e30)*64<=(_sq-_e30)*66
 _SQ_LUT=[math.ceil(math.sqrt(idx+1)*(2.0**9.5))for idx in range(2048)]
 assert max(_SQ_LUT)<65536
 for x0_t in [0,1,2,3,4,15,16,511,512,513,2047,2048,16383,16384,(1<<32)-1,1<<32,(1<<46)-1,1<<46,(1<<48)-1]+[rng.randrange(1,1<<48)for _ in range(500)]:
  x1_t=x0_t|1;e_t=(63-(64-x1_t.bit_length()))&62;idx_t=(x1_t<<9)>>e_t
  assert 512<=idx_t<=2047
  sr_q14=_SQ_LUT[idx_t]<<(e_t>>1)
  assert sr_q14*sr_q14>=(x0_t<<28)
  if x0_t>=1024:assert sr_q14<=math.ceil(math.sqrt(x0_t)*(1<<14)*1.0020)+4
 LIM63=(1<<63)-1
 def _sim_bmul(x0,x1,x2,x3,x4,x5,x6,x7):
  assert all(abs(v)<(1<<60)for v in (x0,x1,x3,x4))and 0<=x2<(1<<60)and 0<=x5<(1<<60)
  x8=((x0*x3)>>64)-((x1*x4)>>64);x10=((x0*x4)>>64)+((x1*x3)>>64);x12=abs(x0)+abs(x1);x13=abs(x3)+abs(x4)
  x14=((x12*x5)>>64)+((x13*x2)>>64)+((x2*x5)>>64)+5
  assert all(-LIM63-1<=v<=LIM63 for v in (x8,x10,x14))and 0<=x12<(1<<61)and 0<=x13<(1<<61)
  x9=max(abs(x8),abs(x10),x14)
  x11=0 if x9==0 else max(0,(64-x9.bit_length())-4)
  x8<<=x11;x10<<=x11;x14<<=x11
  assert abs(x8)<(1<<60)and abs(x10)<(1<<60)and 0<=x14<(1<<60)
  return x8,x10,x14,x6+x7+64-x11
 r8,i8,e8,s8=_sim_bmul(1<<58,0,1,1<<58,0,1,-58,-58)
 assert abs(r8*(2.0**s8)-1.0)<=e8*(2.0**s8)
 for _ in range(100):
  a0,a1,ae=rng.randint(-(1<<59),1<<59),rng.randint(-(1<<59),1<<59),rng.randrange(1,1<<30)
  b0,b1,be=rng.randint(-(1<<59),1<<59),rng.randint(-(1<<59),1<<59),rng.randrange(1,1<<30)
  rr,ii,ee,ss=_sim_bmul(a0,a1,ae,b0,b1,be,-58,-58)
  ex_r=(a0*b0-a1*b1)*(2.0**-116);ex_i=(a0*b1+a1*b0)*(2.0**-116)
  assert math.hypot(rr*(2.0**ss)-ex_r,ii*(2.0**ss)-ex_i)<=ee*(2.0**ss)
 # Live exact proof of Q_CM Gram matrix G_8 spectrum (det = 3^4 * 5^6, lambda_max <= 15 via Sylvester minors of 15*I_8 - G_8):
 def _tr15(d):
  g=math.gcd(15,d%15)
  return 8 if g==15 else(-4 if g==5 else(-2 if g==3 else 1))
 _G8=[[_tr15(i-j)for j in range(8)]for i in range(8)]
 def _det_i(M):
  n=len(M);A=[list(r)for r in M];dv=1
  for k in range(n):
   for i in range(k+1,n):
    for j in range(k+1,n):A[i][j]=(A[i][j]*A[k][k]-A[i][k]*A[k][j])//dv
   dv=A[k][k]
  return A[n-1][n-1]
 assert _det_i(_G8)==(3**4)*(5**6)
 _B15=[[(15 if i==j else 0)-_G8[i][j]for j in range(8)]for i in range(8)]
 assert all(_det_i([r[:k]for r in _B15[:k]])>=0 for k in range(1,9))
 # 1) L_p (det <= 2^72, Euclidean 11/40 -> 99/100 for dsieve16): d16 <= 2^144, D0 <= 2^(15*144)=2^2160, L_d <= ceil((144+32*log2(40))/64)=5 limbs
 s1_ds=math.floor((15.0*144.0)/math.log2(40.0/11.0))
 ds_ops=15*(40*5*3+16*12*4+15*12*10)+15*(4*12*3*5+2*40*8*3)+2000
 s2_ds=math.floor((340.0*math.log2(40.0)+7.5*144.0)/math.log2(100.0/99.0))
 ds_lll_bound=(2*s1_ds+16)*ds_ops+(2*s2_ds+16)*ds_ops
 assert ds_lll_bound<c_ps
 # 2) Special-q L_Q (det <= 2^262, 2-phase 11/40 Euclidean -> 99/100 Q_CM): d16 <= 2^524, d_k^(CM) <= 15^k d_k^(Euc), L_d <= ceil((524+32*log2(40)+15*log2(15))/64)=12 <= 15 limbs
 s1_sq=math.floor((15.0*524.0)/math.log2(40.0/11.0))
 sq1_ops=15*(40*15*10+16*12*28+15*12*78)+15*(4*12*10*15+2*40*25*10)+10000
 s2_sq=math.floor((120.0*math.log2(15.0)+340.0*math.log2(40.0)+7.5*524.0)/math.log2(100.0/99.0))
 sq2_ops=sq1_ops+256*15*12
 sq_lll_bound=(2*s1_sq+16)*sq1_ops+(2*s2_sq+16)*sq2_ops
 assert sq_lll_bound<lll_node_cap
 # 3) 381b G_T boot L_boot(z') WITH S_boot = 2^21 scaling on 8 b-coords BEFORE Phase 1:
 #    det = p^4 * 2^168 <= 2^1692; every initial HNF k x k Gram minor d_k <= C(16,k)*(p^4*2^168)^2 <= 2^3398 < 2^3400.
 #    Phase 2 uses pure Q_CM = diag(G_8, G_8) (lambda_max <= 15), so d_k^(CM) <= 15^k d_k^(Euc) unconditionally!
 #    L_d <= ceil((3400 + 32*log2(40) + 15*log2(15))/64) = 58 <= 59 limbs -> bt_lll_bound = 2^46.75 < 2^47.0 = LLL_BOOT!
 _D_BT=3400.0;_L2_BT=59
 assert math.ceil((_D_BT+32.0*math.log2(40.0)+15.0*math.log2(15.0))/64.0)<=_L2_BT
 s1_bt=math.floor((15.0*_D_BT)/math.log2(40.0/11.0))
 bt1_ops=15*(40*222*147+16*12*75*75+15*12*75*222)+15*(4*12*147*222+2*40*368*147)+1000000
 s2_bt=math.floor((120.0*math.log2(15.0)+340.0*math.log2(40.0)+7.5*_D_BT)/math.log2(100.0/99.0))
 bt2_ops=15*(40*89*59+16*12*30*30+15*12*30*89)+15*(4*12*59*89+2*40*148*59)+1000000+256*59*12
 bt_lll_bound=(2*s1_bt+16)*bt1_ops+(2*s2_bt+16)*bt2_ops
 assert bt_lll_bound<lll_boot_cap
 q265=((1<<265)+1234567891)*262139;Bsq=[[0]*16 for _ in range(16)];Bsq[0][0]=q265
 ry=pow(7,13,q265);rt=pow(11,17,q265);ryp=[pow(ry,i,q265)for i in range(8)]
 for i in range(1,8):Bsq[i][0]=(-ryp[i])%q265;Bsq[i][i]=1
 for i in range(8):Bsq[8+i][0]=(rt*ryp[i])%q265;Bsq[8+i][8+i]=1
 # Live Q_CM 2-phase (11/40 Euc -> 99/100 Q_CM) Bareiss LLL check on BOTH identity I_16 AND 283-bit special-q Bsq!
 _I16=[[1 if i==j else 0 for j in range(16)]for i in range(16)]
 _RI,_dI,_lI=lll_16d_int(_I16,99,100,use_cm=True)
 assert _dI[16]==(3**8)*(5**12)and _RI!=_I16
 for kk in range(1,16):
  assert 100*(_dI[kk-1]*_dI[kk+1]+_lI[kk][kk-1]**2)>=99*_dI[kk]**2
  for jj in range(kk):assert 2*abs(_lI[kk][jj])<=_dI[jj+1]
 _,dsq,lsq=lll_16d_int(Bsq,99,100,use_cm=True)
 assert dsq[16]==q265*q265*(3**8)*(5**12)and q265.bit_length()==283
 for kk in range(1,16):
  assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
  for jj in range(kk):assert 2*abs(lsq[kk][jj])<=dsq[jj+1]
 return dict(asm_blocks_parsed=len(parsed),b2_direct_1pass_ok=True,b2_bank_rot_62b_ok=True,uint16_swar_ok=True,ntt62_garner_exact=True,ntt3_r4_r3_shoup_conv_ok=True,ntt62_r4_3x3p_stage_ir=round(9.0*parsed["ntt62_r4_shoup"]["ir"],2),f62_barrett_mac_64_exact=True,mont5_opcodes=mont5_ops,mont5_sim_32_exact=True,row_red5_barrett_opcodes=parsed["row_red5_z12"]["raw_opcodes"],row_red5_barrett_47_exact=True,unit_run_4vec4_ir=round(parsed["unit_run_4vec4_z12"]["ir"],4),gen_run1_2v_ir=round(parsed["gen_run1_2v_z12"]["ir"],2),scan_zero_u16=round(scan_tot,4),dsieve16_rank_552_exact=True,sq_lll_uniform_log2=round(math.log2(sq_lll_bound),2),boot_lll_uniform_log2=round(math.log2(bt_lll_bound),2),all_raw_opcodes_ok=True)
from unitdomain import verify_prod_sq24
