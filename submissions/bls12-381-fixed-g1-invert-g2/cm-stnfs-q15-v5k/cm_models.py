"""cm_models.py — Self-contained."""
import array
import collections
import functools
import math
from decimal import Decimal,localcontext
LN2=math.log(2.0)
GAMMA=0.5772156649015329
MONT5=367.76
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
def rho_exact(u):
 if u<=1.0:return 1.0
 k=int(u)
 if u>k:k+=1
 if k<_RHO_C_LEN:
  x,v=k-u,0.0
  for c in _RHO_C_REV[k]:
   v=v*x+c
  return v if v>0.0 else 0.0
 if u>=_RHO_UMAX:return 0.0
 pos=u*_RHO_S
 idx=int(pos)
 frac=pos-idx
 y0,y1=_RHO_TAIL[idx],_RHO_TAIL[idx+1]
 if y0<=0.0 or y1<=0.0:return 0.0
 return math.exp((1.0-frac)*math.log(y0)+frac*math.log(y1))
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
 L=length
 for n in range(1,L+1):
  val=_sp(jw_rev[L-n :L],out[:n])
  out[n]=val/n
  if val <-1e-13:raise ArithmeticError("negative coefficient beyond roundoff tolerance")
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
def _semismooth_ecm(fbits,b1,b2,steps=256):
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
def lll_16d_int(B_in,delta_num=99,delta_den=100):
 B=[list(r) for r in B_in]
 n=len(B)
 _sp=math.sumprod
 if max(abs(x).bit_length() for r in B for x in r)>160:
  with localcontext() as ctx:
   ctx.prec=280
   ZERO,HALF,D55=Decimal(0),Decimal("0.5"),Decimal("0.55")
   mu=[[ZERO]*n for _ in range(n)]
   Bn=[ZERO]*n
   def _ref():
    Bs=[[ZERO]*n for _ in range(n)]
    Bd=[[Decimal(x) for x in r] for r in B]
    for i in range(n):
     bi=Bd[i];vv=bi[:];mui=mu[i]
     for j in range(i):
      mij=_sp(bi,Bs[j])/Bn[j]
      mui[j]=mij;bsj=Bs[j]
      vv=[vv[t]-mij*bsj[t] for t in range(n)]
     Bs[i]=vv;Bn[i]=_sp(vv,vv)
   _ref()
   k=1;sw=0
   while k<n:
    muk=mu[k]
    for j in range(k-1,-1,-1):
     m=muk[j]
     if m>HALF or m<-HALF:
      q=int(m.to_integral_value())
      if q:
       qd=Decimal(q);bk,bj=B[k],B[j]
       B[k]=[x-q*y for x,y in zip(bk,bj)]
       muj=mu[j]
       for l in range(j):muk[l]-=qd*muj[l]
       muk[j]-=qd
    if Bn[k]>=(D55-muk[k-1]*muk[k-1])*Bn[k-1]:k+=1
    else:
     B[k],B[k-1]=B[k-1],B[k]
     m_k=muk[k-1];Bn_new=Bn[k]+m_k*m_k*Bn[k-1];mu_new=m_k*Bn[k-1]/Bn_new
     Bn[k]=Bn[k]*Bn[k-1]/Bn_new;Bn[k-1]=Bn_new
     mu[k][k-1]=mu_new
     for j in range(k-1):mu[k][j],mu[k-1][j]=mu[k-1][j],mu[k][j]
     for i in range(k+1,n):
      mi=mu[i];t=mi[k];mi[k]=mi[k-1]-m_k*t;mi[k-1]=t+mu_new*mi[k]
     sw+=1;k=max(k-1,1)
     if sw%400==0:_ref()
  B=lll_16d(B)
 d=[1]*(n+1)
 lam=[[0]*n for _ in range(n)]
 for i in range(n):
  bi=B[i]
  for j in range(i+1):
   u=_sp(bi,B[j])
   for k in range(j):
    u=(d[k+1]*u-lam[i][k]*lam[j][k])//d[k]
   if j<i:lam[i][j]=u
   else:d[i+1]=u
 k=1
 while k<n:
  lk=lam[k]
  if 2*abs(lk[k-1])>d[k]:
   q=(2*lk[k-1]+d[k])//(2*d[k]) if lk[k-1]>=0 else -((-2*lk[k-1]+d[k])//(2*d[k]))
   B[k]=[x-q*y for x,y in zip(B[k],B[k-1])]
   lk1=lam[k-1]
   for j in range(k-1):lk[j]-=q*lk1[j]
   lk[k-1]-=q*d[k]
  if delta_den*(d[k-1]*d[k+1]+lk[k-1]*lk[k-1])<delta_num*d[k]*d[k]:
   B[k],B[k-1]=B[k-1],B[k]
   lk1=lam[k-1]
   for j in range(k-1):lk[j],lk1[j]=lk1[j],lk[j]
   dk_new=(d[k-1]*d[k+1]+lk[k-1]*lk[k-1])//d[k]
   lkk1,dk,dk1=lk[k-1],d[k],d[k+1]
   for i in range(k+1,n):
    li=lam[i];t=li[k]
    li[k]=(dk1*li[k-1]-lkk1*t)//dk
    li[k-1]=(dk_new*t+lkk1*li[k])//dk1
   d[k]=dk_new;k=max(k-1,1)
  else:
   for l in range(k-2,-1,-1):
    if 2*abs(lk[l])>d[l+1]:
     q=(2*lk[l]+d[l+1])//(2*d[l+1]) if lk[l]>=0 else -((-2*lk[l]+d[l+1])//(2*d[l+1]))
     B[k]=[x-q*y for x,y in zip(B[k],B[l])]
     ll=lam[l]
     for j in range(l):lk[j]-=q*ll[j]
     lk[l]-=q*d[l+1]
   k+=1
 return B,d,lam
def run_boot_381bit_lll_live(gam,h_pi,p_bls,s_boot,pmul_fn,prem_fn):
 B_sk=[]
 for i in range(4):
  r=[0]*16;r[i]=p_bls;B_sk.append(r)
 for j in range(4):
  r=[0]*16
  for idx,ch in enumerate(h_pi):r[j+idx]=ch
  B_sk.append(r)
 for k in range(8):
  yk=[0]*8;yk[k]=1
  ak=pmul_fn(gam,prem_fn(yk,h_pi,p_bls),p_bls,h_pi)
  r=[0]*16
  for idx,ch in enumerate(ak):r[idx]=ch
  r[8+k]=s_boot;B_sk.append(r)
 R_sk,d_gso,lam=lll_16d_int(B_sk,99,100)
 assert d_gso[16]==(p_bls**8)*(s_boot**16)
 for kk in range(1,16):
  assert 100*(d_gso[kk-1]*d_gso[kk+1]+lam[kk][kk-1]**2)>=99*d_gso[kk]**2
  for jj in range(kk):
   assert 2*abs(lam[kk][jj])<=d_gso[jj+1]
 assert all(r[8+kk]%s_boot==0 for r in R_sk for kk in range(8))
 return R_sk
def lll_16d(B):
 B=[r[:] for r in B]
 n=len(B)
 Bs=[[0.0]*n for _ in range(n)]
 mu=[[0.0]*n for _ in range(n)]
 Bn=[0.0]*n
 _sp=math.sumprod
 for i in range(n):
  bi=B[i]
  vv=[float(x) for x in bi]
  mui=mu[i]
  for j in range(i):
   mij=_sp(bi,Bs[j])/Bn[j]
   mui[j]=mij
   bsj=Bs[j]
   vv=[vv[t]-mij*bsj[t] for t in range(n)]
  Bs[i]=vv
  Bn[i]=_sp(vv,vv)
 k=1
 while k<n:
  muk=mu[k]
  for j in range(k-1,-1,-1):
   m=muk[j]
   if m>0.5 or m <-0.5:
    q=round(m)
    bk,bj=B[k],B[j]
    B[k]=[bk[t]-q*bj[t] for t in range(n)]
    muj=mu[j]
    for l in range(j):
     muk[l]-=q*muj[l]
    muk[j]-=q
  if Bn[k]>=(0.99-muk[k-1]*muk[k-1])*Bn[k-1]:k+=1
  else:
   B[k],B[k-1]=B[k-1],B[k]
   m_k=muk[k-1]
   Bn_new=Bn[k]+m_k*m_k*Bn[k-1]
   mu_new=m_k*Bn[k-1]/Bn_new
   Bs_old,Bs_k=Bs[k-1],Bs[k]
   Bs[k-1]=[Bs_k[t]+m_k*Bs_old[t] for t in range(n)]
   Bs[k]=[Bs_old[t]-mu_new*Bs[k-1][t] for t in range(n)]
   Bn[k]=Bn[k]*Bn[k-1]/Bn_new
   Bn[k-1]=Bn_new
   mu[k][k-1]=mu_new
   for j in range(k-1):
    mu[k][j],mu[k-1][j]=mu[k-1][j],mu[k][j]
   for i in range(k+1,n):
    t=mu[i][k]
    mu[i][k]=mu[i][k-1]-m_k*t
    mu[i][k-1]=t+mu_new*mu[i][k]
   k=max(k-1,1)
 return B
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
 rf=rho_exact(u_f)
 rg=rho_exact(u_g)
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
 def solve_h_newton(scale,stop90=False):
  h=0.0
  sbins=[(wt,lam*scale) for wt,lam in bins]
  for _ in range(400):
   omh=1.0-h
   s=ds=0.0
   for wt,ls in sbins:
    t=wt*_exp(-ls*omh)
    s+=t
    ds+=t*ls
   em=_exp(-s)
   g=1.0-em
   if stop90 and g>0.90:return g
   diff=g-h
   if abs(diff)<1e-13:return g
   gp=em*ds
   if 0.0<=gp<0.98 and diff>0.0:
    hn=h+0.85*diff/(1.0-gp)
    if g<hn<1.0:
     omhn=1.0-hn
     sn=sum(wt*_exp(-ls*omhn) for wt,ls in sbins)
     gn=1.0-_exp(-sn)
     if stop90 and gn>0.90:return gn
     if gn>=hn:
      h=gn
      continue
   h=g
  return h
 def exceeds_90(scale):
  return solve_h_newton(scale,True)>0.90
 h=solve_h_newton(1.0)
 lo_e,hi_e=0.4,exc_ref
 for _ in range(24):
  mid=0.5*(lo_e+hi_e)
  if exceeds_90(mid/exc_ref):lo_e=mid
  else:hi_e=mid
 return h,hi_e
def eval_descent_live(desc,mesh_node=96,mesh_boot=192):
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
 PAIR_SETUP=(2.0**25.1)*(smult*1.33)
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
  def conv(a,b):
   out=[0.0]*(len(a)+len(b)-1)
   nz=[j for j,y in enumerate(b) if y]
   if not nz:return out
   j0,j1=nz[0],nz[-1]+1
   bb=b[j0:j1]
   for i,x in enumerate(a):
    if x:
     for k,y in enumerate(bb):
      out[i+j0+k]+=x*y
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
 W_det,PHI_dag,C_det={},{},{}
 bins_seen=[]
 out_rows=[]
 scale_iA=mesh_node//32
 p_ecm=float(dp["ecm_success"])
 for row in desc["plan_struct"]:
  lo,hi,side,n_int,iA32,kA,kB,R_n,thin=row
  bn=(lo,hi)
  if not bins_seen or bins_seen[-1]!=bn:bins_seen.append(bn)
  n_c=float(n_int)
  iA=iA32*scale_iA
  cells=cells_r2[R_n]
  S=math.log2(cells)
  c_hi=(hi+S)/16.0
  c_lo=(lo+S)/16.0
  nA=lx+iA*(n_c-lx)/mesh_node if n_c>lx else lx
  T=kA*nA+kB*n_c+96.0
  nf_list=[72.0*c_hi+cf-(lo if side=="f" else 0.0) for _,cf,_ in qbins]+[
   72.0*c_hi+cf-(lo if side=="f" else 0.0) for cf,_ in qbins_sub
  ]
  ng_list=[24.0*c_hi+log2_u8+cg-(lo if side=="g" else 0.0) for _,_,cg in qbins]+[
   24.0*c_hi+log2_u8+cg-(lo if side=="g" else 0.0) for _,cg in qbins_sub
  ]
  pf_v=ss_eval(nf_list,n_c,iA,mesh_node,alpha_f,kA,kB)
  pg_v=ss_eval(ng_list,n_c,iA,mesh_node,alpha_g,kA,kB)
  pr_c=sum(pf_v[i]*pg_v[i] for i in range(6))/6.0
  pr_i=sum(pf_v[i]*pg_v[i] for i in range(6,30))/24.0
  pr=min(pr_i,4.0*pr_c)
  sv=0.0
  for (w,_,_),(cfl,cgl) in zip(qbins,qbins_lo):
   nfl=72.0*c_lo+cfl-(hi if side=="f" else 0.0)
   ngl=24.0*c_lo+log2_u8+cgl-(hi if side=="g" else 0.0)
   sf=1.0 if T>=nfl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(nfl)-T)/lx))
   sg=1.0 if T>=ngl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(ngl)-T)/lx))
   sv+=w*sf*sg
  nf_max=72.0*c_hi+cfmax-(lo if side=="f" else 0.0)
  ng_max=24.0*c_hi+log2_u8+cgmax-(lo if side=="g" else 0.0)
  k_lx_f=int(nf_max//20.0)+1
  k_lx_g=int(ng_max//20.0)+1
  k_wA_n=kA if nA>lx+1e-9 else 0
  det_n=p_ecm**(k_lx_f+k_lx_g+2*(kA+kB))
  miss_n=miss_side(nf_max,uf_peel,lf_peel)+miss_side(ng_max,ug_peel,lg_peel)
  mu=cells*pr*float(dp["udup_node"])*det_n*math.exp(-miss_n)
  mr_f=mr_instr(nf_max)
  mr_g=mr_instr(ng_max)
  mr_T=mr_instr(T)
  surv_cost=(
   NORMS+TDIV
   +k_lx_f*(ecm_instr(lx,nf_max,p_ecm)+mr_f)
   +k_lx_g*(ecm_instr(lx,ng_max,p_ecm)+mr_g)
   +2.0*k_wA_n*(ecm_instr(nA,T,p_ecm)+mr_T)
   +mr_f+mr_g
   +IDEAL
  )
  survivors=cells*sv+cells*pr
  sieve=cells*(HITS_PER_CELL*SIEVE_HIT+SCAN+ZERO+BOUND_CELL)+pairs*PAIR_SETUP
  if cells*surv_cost<sieve+survivors*surv_cost:
   how="ecm"
   own=LLL_NODE+cells*surv_cost
  else:
   how="sieve"
   mu*=float(dp["ret_node"])
   own=LLL_NODE+sieve+survivors*surv_cost
  wmax=lambda s,lim:max([W_det[b,s] for b in bins_seen if b[0]<lim-1e-9 and (b,s) in W_det]+[0.0])
  pmax=lambda s,lim:max([PHI_dag[b,s] for b in bins_seen if b[0]<lim-1e-9 and (b,s) in PHI_dag]+[0.0])
  cmax=lambda s,lim:max([C_det[b,s] for b in bins_seen if b[0]<lim-1e-9 and (b,s) in C_det]+[0.0])
  p0=1.0/(1.0+mu)
  W_det[bn,side]=own+kA*(wmax("f",nA)+wmax("g",nA))+kB*(wmax("f",n_c)+wmax("g",n_c))
  PHI_dag[bn,side]=p0+kA*(pmax("f",nA)+pmax("g",nA))+kB*(pmax("f",n_c)+pmax("g",n_c))
  C_det[bn,side]=1+kA*(cmax("f",nA)+cmax("g",nA))+kB*(cmax("f",n_c)+cmax("g",n_c))
  out_rows.append(dict(
   lo=lo,hi=hi,side=side,how=how,R=R_n,n=n_int,nA=nA,kA=kA,kB=kB,
   mu=mu,p0=p0,phi_dag=PHI_dag[bn,side],nf=nf_max,ng=ng_max,
   own_log2=math.log2(own),W_log2=math.log2(W_det[bn,side]),C_det=C_det[bn,side],
  ))
 bt=desc["boot"]
 assert int(bt["cells_b"])==1 and int(dp["boot_rels"])==1 and float(dp.get("ball_pair",0.0))==0.0
 m_max=float(bt["m_sv"])
 iA_b=round((bt["nA"]-lx)*mesh_boot/(bt["n1"]-lx))
 pr_raw_b=ss_eval([m_max],float(bt["n1"]),iA_b,mesh_boot,alpha_g,bt["kA"],1)[0]
 mb=miss_side(m_max,ug_peel,lg_peel)
 pr_b=pr_raw_b*math.exp(-mb)
 k_lx_b=int(m_max//20.0)+1
 kA_b=int(bt["kA"])
 mr_b=mr_instr(m_max)
 det_b=p_ecm**(k_lx_b+kA_b+1)
 p_boot_trial=pr_b*det_b*float(dp["ret_boot"])
 trial_b=(
  NORMS+TDIV
  +k_lx_b*(ecm_instr(lx,m_max,p_ecm)+mr_b)
  +kA_b*(ecm_instr(bt["nA"],m_max,p_ecm)+mr_b)
  +mr_b+IDEAL
 )
 boot_cap_p=float(dp.get("boot_trial_cap_p",1e-4))
 T_boot_max=math.ceil(math.log(1.0/boot_cap_p)/p_boot_trial)
 boot_own_exp=(1.0/p_boot_trial)*(LLL_NODE+trial_b)
 boot_own_cap=T_boot_max*(LLL_NODE+trial_b)
 wmax_g=lambda lim: max(W_det[b,"g"] for b in bins_seen if b[0] < lim-1e-9)
 pmax_g=lambda lim: max(PHI_dag[b,"g"] for b in bins_seen if b[0] < lim-1e-9)
 cmax_g=lambda lim: max(C_det[b,"g"] for b in bins_seen if b[0] < lim-1e-9)
 W_dag_1=bt["kA"]*wmax_g(bt["nA"])+wmax_g(float(bt["n1"]))
 Phi_1=bt["kA"]*pmax_g(bt["nA"])+pmax_g(float(bt["n1"]))
 C_det_1=bt["kA"]*cmax_g(bt["nA"])+cmax_g(float(bt["n1"]))
 b_max=int(dp.get("b_max",10))
 Q_1=Phi_1+boot_cap_p+(2.0**-254)
 assert Q_1<0.50
 two_root_fail=2.0*(Q_1**b_max)
 wald_mult=(1.0-(Q_1**b_max))/(1.0-Q_1)
 two_exp_u=2.0*wald_mult*(boot_own_exp+W_dag_1)/20000.0
 two_cap_u=2.0*b_max*(boot_own_cap+W_dag_1)/20000.0
 return dict(
  mesh_node=mesh_node,mesh_boot=mesh_boot,n_nodes=len(desc["plan_struct"]),
  peel_h=peel_h,exc_crit=exc_crit,b_max=b_max,
  boot_pr_log2=math.log2(pr_b),boot_miss=mb,
  p_boot_trial=p_boot_trial,T_boot_max=T_boot_max,boot_trial_cap_p=boot_cap_p,
  boot_own_exp_log2=math.log2(2.0*boot_own_exp/20000.0),
  boot_own_cap_log2=math.log2(2.0*boot_own_cap/20000.0),
  dag_sub_two_log2=math.log2(2.0*W_dag_1/20000.0),
  phi_1log=Phi_1,q1_attempt=Q_1,wald_mult=wald_mult,
  two_log_exp_units=math.log2(two_exp_u),
  two_log_cap_units=math.log2(two_cap_u),
  two_root_fail=two_root_fail,
  C_det_1log=C_det_1,
  min_mu=min(r["mu"] for r in out_rows),
  p0_max=max(r["p0"] for r in out_rows),
  rows=out_rows,
 )
def q15_v5k_descent_live_replay(camp):
 import random
 assert camp["same_field_q15"] and camp["real_2core_fb"] and camp["literal_k1_first_usable"]
 b_fb_b=int(camp["b_fb_bits"]);B_FB=1<<b_fb_b
 coll_seed=int(camp["coll_seed"]);coll_pairs=int(camp["coll_pairs"])
 d1_lo,d1_hi,nc_1=float(camp["d1_lo"]),float(camp["d1_hi"]),float(camp["nc_1"])
 d2_lo,d2_hi,nc_2=float(camp["d2_lo"]),float(camp["d2_hi"]),float(camp["nc_2"])
 nc_0=float(camp["nc_0"]);af_t,ag_t=float(camp["alpha_f"]),float(camp["alpha_g"])
 phi15=(1,-1,0,1,-1,1,0,-1,1)
 ypow=[];cc=[1,0,0,0,0,0,0,0]
 for _ in range(15):
  ypow.append(tuple(cc));hi=cc[7]
  cc=[-hi*phi15[0]]+[cc[i-1]-hi*phi15[i] for i in range(1,8)]
 om2=list(ypow[10]);two_om=[2+ypow[5][0]]+list(ypow[5][1:])
 def _cm(k):
  rows=[[] for _ in range(8)]
  for j in range(8):
   yk=ypow[(j*k)%15]
   for i in range(8):
    if yk[i]:rows[i].append((j,yk[i]))
  return rows
 C14,C4,C7=_cm(14),_cm(4),_cm(7)
 def mul15(a,b,S=36):
  M=1<<S;H=M>>1;mk=M-1
  A=a[0]+(a[1]<<S)+(a[2]<<(2*S))+(a[3]<<(3*S))+(a[4]<<(4*S))+(a[5]<<(5*S))+(a[6]<<(6*S))+(a[7]<<(7*S))
  B=b[0]+(b[1]<<S)+(b[2]<<(2*S))+(b[3]<<(3*S))+(b[4]<<(4*S))+(b[5]<<(5*S))+(b[6]<<(6*S))+(b[7]<<(7*S))
  P=A*B;c=[0]*15
  for i in range(15):
   v=P&mk
   if v>=H:v-=M;P+=M
   c[i]=v;P>>=S
  for i in range(14,7,-1):
   v=c[i]
   if v:c[i-8]-=v;c[i-7]+=v;c[i-5]-=v;c[i-4]+=v;c[i-3]-=v;c[i-1]+=v
  return c[:8]
 def norm15(a):
  p1=mul15(a,(a[0]+a[1]+a[2]-a[5]-a[6]-a[7],-a[2]-a[4]+a[7],-a[1]-a[3]+a[6],a[1]-a[6]-a[7],-a[1]-a[2]+a[7],a[2]-a[5]-a[7],a[1]-a[4]-a[6],-a[1]-a[2]-a[3]+a[6]+a[7]),40)
  p2=mul15(p1,(p1[0]-p1[2]-p1[6]+p1[7],p1[2]+p1[4]-p1[7],-p1[3]+p1[6],-p1[2]-p1[6],p1[1]+p1[2]-p1[7],-p1[2]+p1[5]+p1[7],-p1[6],p1[2]-p1[3]+p1[6]-p1[7]),80)
  p3=mul15(p2,(p2[0]+p2[2]+p2[4],-p2[4],-p2[2]-p2[6],p2[2],-p2[2]-p2[4]+p2[7],p2[4]+p2[5],p2[2]+p2[3],p2[1]-p2[2]-p2[4]-p2[6]),160)
  return abs(p3[0])
 def _prod(xs):
  if not xs:return 1
  while len(xs)>1:xs=[xs[i]*xs[i+1] if i+1<len(xs) else xs[i] for i in range(0,len(xs),2)]
  return xs[0]
 ps=primes_upto(B_FB)
 P256=_prod([p for p in ps if p<=256])
 PMID=_prod([p for p in ps if 256<p<=B_FB])
 def bstrip(vals):
  if not vals:return []
  rems=[]
  for n in vals:
   while n>1:
    g=math.gcd(n,P256)
    if g==1:break
    n//=g
   rems.append(n)
  tree=[rems]
  while len(tree[-1])>1:
   cur=tree[-1];tree.append([cur[i]*cur[i+1] if i+1<len(cur) else cur[i] for i in range(0,len(cur),2)])
  R=[PMID%tree[-1][0]]
  for lv in range(len(tree)-2,-1,-1):
   cur=tree[lv];R=[R[i>>1]%cur[i] for i in range(len(cur))]
  out=[]
  for n,r in zip(rems,R):
   n0,r0=n,r
   while n0>1:
    g=math.gcd(n0,r0)
    if g==1:break
    n0//=g;r0%=n0
   out.append(n0)
  return out
 rng=random.Random(coll_seed);ngs,nfs=[],[]
 for _ in range(coll_pairs):
  a=[0]*8;b=[0]*8
  for _ in range(rng.choice((2,3))):a[rng.randrange(8)]+=rng.choice((-1,1))
  for _ in range(rng.choice((1,2))):b[rng.randrange(8)]+=rng.choice((-1,1))
  if not any(a) or not any(b):continue
  a3=mul15(mul15(a,a,24),a,32);b3=mul15(mul15(b,b,24),b,32)
  ng=norm15([a3[j]+2*mul15(om2,b3,32)[j] for j in range(8)])
  if ng==0:continue
  nf=norm15([a3[j]-mul15(two_om,b3,32)[j] for j in range(8)])
  if nf==0:continue
  ngs.append(ng);nfs.append(nf)
 rg=bstrip(ngs);keep=[i for i,x in enumerate(rg) if x==1]
 rf=bstrip([nfs[i] for i in keep])
 rels=[(nfs[keep[j]],ngs[keep[j]]) for j,x in enumerate(rf) if x==1]
 def fac_s(n):
  out=set()
  for p in ps:
   if p*p>n:break
   if n%p==0:
    out.add(p)
    while n%p==0:n//=p
  if n>1:out.add(n)
  return out
 rel_facs=[(fac_s(nf),fac_s(ng)) for nf,ng in rels]
 cur=list(range(len(rel_facs)))
 while True:
  deg={}
  for ri in cur:
   ff,fg=rel_facs[ri]
   for p in ff:deg[(1,p)]=deg.get((1,p),0)+1
   for p in fg:deg[(0,p)]=deg.get((0,p),0)+1
  nxt=[ri for ri in cur if all(deg[(1,p)]>=2 for p in rel_facs[ri][0]) and all(deg[(0,p)]>=2 for p in rel_facs[ri][1])]
  if len(nxt)==len(cur):break
  cur=nxt
 fb_g={p for ri in cur for p in rel_facs[ri][1]}
 fb_f={p for ri in cur for p in rel_facs[ri][0]}
 assert [len(rels),len(cur)]==camp["raw_rels_2core"] and [len(fb_g),len(fb_f)]==camp["fb_2core_sizes"]
 p_miss_g=_prod([p for p in ps if p not in fb_g])
 p_miss_f=_prod([p for p in ps if p not in fb_f])
 miss_sum_g=sum(1.0/p for p in ps if p>256 and p%15==1 and p not in fb_g)
 miss_sum_f=sum(1.0/p for p in ps if p>256 and p%15==1 and p not in fb_f)
 n_pts=int(camp["ball_pts"])
 ball_pts=[]
 for i in range(16):ball_pts.append(((i,1),))
 for i in range(16):
  for j in range(i+1,16):
   for s in (1,-1):ball_pts.append(((i,1),(j,s)))
 for i in range(16):
  for j in range(i+1,16):
   for k in range(j+1,16):
    for sj in (1,-1):
     for sk in (1,-1):
      if len(ball_pts)<n_pts:ball_pts.append(((i,1),(j,sj),(k,sk)))
 assert len(ball_pts)==n_pts
 S_b=math.log2(n_pts)
 def prim15(p):
  g=2
  while True:
   z=pow(g,(p-1)//15,p)
   if pow(z,3,p)!=1 and pow(z,5,p)!=1 and pow(z,15,p)==1:return z
   g+=1
 def sq16_b(q,y0,r0):
  B=[[0]*16 for _ in range(16)];B[0][0]=q
  for i in range(1,8):B[i][0]=(-pow(y0,i,q))%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(r0*pow(y0,i,q))%q;B[8+i][8+i]=1
  return lll_16d(B)
 def ext_rt(pc,cs,a,b):
  z=prim15(pc)
  for k in (1,2,4,7,8,11,13,14):
   y0=pow(z,k,pc)
   bs=sum(b[j]*pow(y0,j,pc) for j in range(8))%pc
   if bs==0:continue
   as_=sum(a[j]*pow(y0,j,pc) for j in range(8))%pc
   r0=(as_*pow(bs,-1,pc))%pc;om=pow(y0,5,pc)
   rhs=(-2*om*om)%pc if cs==0 else (2+om)%pc
   if pow(r0,3,pc)==rhs:return (y0,r0)
  return None
 def eval_ball(side,q,y0,r0,nc_b,stop_first=False):
  R=sq16_b(q,y0,r0);nc_v=int(2**nc_b);assert nc_v<B_FB*B_FB
  step=64 if stop_first else n_pts;M_use=sv_2=0;first_chs=None
  for off in range(0,n_pts,step):
   pts,nfs2,ngs2=[],[],[]
   for sp_idx in range(off,min(n_pts,off+step)):
    a=[0]*8;b=[0]*8
    for idx,sgn in ball_pts[sp_idx]:
     row=R[idx]
     if sgn==1:
      for t in range(8):a[t]+=row[t];b[t]+=row[8+t]
     else:
      for t in range(8):a[t]-=row[t];b[t]-=row[8+t]
    if not any(b):continue
    a3=mul15(mul15(a,a,24),a,36);b3=mul15(mul15(b,b,24),b,36)
    ng=norm15([a3[j]+2*mul15(om2,b3,36)[j] for j in range(8)])
    nf=norm15([a3[j]-mul15(two_om,b3,36)[j] for j in range(8)])
    if ng==0 or nf==0:continue
    if side==0:assert ng%q==0;ng//=q
    else:assert nf%q==0;nf//=q
    pts.append((a,b));nfs2.append(nf);ngs2.append(ng)
   rgs2=bstrip(ngs2);kp=[i for i,x in enumerate(rgs2) if x<=nc_v]
   rfs2=bstrip([nfs2[i] for i in kp])
   for j,idx in enumerate(kp):
    rf=rfs2[j]
    if rf<=nc_v:
     sv_2+=1;rg=rgs2[idx]
     if math.gcd(ngs2[idx],p_miss_g)==1 and math.gcd(nfs2[idx],p_miss_f)==1:
      M_use+=1
      if M_use==1:
       a,b=pts[idx];cands=[]
       if rf>B_FB and rf%15==1:
        rt=ext_rt(rf,1,a,b)
        if rt is not None:cands.append((1,rf,rt[0],rt[1]))
       if rg>B_FB and rg%15==1:
        rt=ext_rt(rg,0,a,b)
        if rt is not None:cands.append((0,rg,rt[0],rt[1]))
       cands.sort(key=lambda x:(-x[1],-x[0]));first_chs=cands
       if stop_first:return M_use,sv_2,first_chs
  return M_use,sv_2,(first_chs or [])
 qb,qbs,qbl=camp["qbins"],camp["qbins_sub"],camp["qbins_lo"]
 def calc_env(clo,chi,nc_c):
  def ss_t(m,alpha):
   mesh=48;h=(nc_c-b_fb_b)/mesh
   wB=[math.log1p(h/(b_fb_b+i*h)) for i in range(mesh)]
   tot=0.0
   for jB,row in ((0,[1.0]),(1,wB)):
    for idx,c in enumerate(row):
     if not c:continue
     rem=m-jB*b_fb_b-idx*h
     if rem<0:continue
     u=(rem+alpha)/b_fb_b
     tot+=c*(rho_db(u,max(float(b_fb_b),rem)) if u>2.0 else rho_exact(u))
   return min(1.0,tot)
  T_c=nc_c;c_hi,c_lo=(chi+S_b)/16.0,(clo+S_b)/16.0
  mu_u,sv_up=[],[]
  for side in (0,1):
   pr_c=sum(ss_t(24.0*c_hi+cf-(clo if side==1 else 0.0),af_t)*ss_t(24.0*c_hi+cg-(clo if side==0 else 0.0),ag_t) for cf,cg in qb)/6.0
   pr_i=sum(ss_t(24.0*c_hi+cf-(clo if side==1 else 0.0),af_t)*ss_t(24.0*c_hi+cg-(clo if side==0 else 0.0),ag_t) for cf,cg in qbs)/24.0
   pr=min(pr_i,4.0*pr_c)
   nf_hi=24.0*c_hi+qb[-1][0]-(clo if side==1 else 0.0)
   ng_hi=24.0*c_hi+qb[-1][1]-(clo if side==0 else 0.0)
   miss_tot=(nf_hi/b_fb_b)*miss_sum_f+(ng_hi/b_fb_b)*miss_sum_g
   sv=0.0
   for cfl,cgl in qbl:
    nfl=24.0*c_lo+cfl-(chi if side==1 else 0.0)
    ngl=24.0*c_lo+cgl-(chi if side==0 else 0.0)
    sf=1.0 if T_c>=nfl else min(1.0,max(1.0,T_c/b_fb_b)*rho_exact((math.floor(nfl)-T_c)/b_fb_b))
    sg=1.0 if T_c>=ngl else min(1.0,max(1.0,T_c/b_fb_b)*rho_exact((math.floor(ngl)-T_c)/b_fb_b))
    sv+=sf*sg/6.0
   mu_u.append(n_pts*pr*math.exp(-miss_tot))
   sv_up.append(n_pts*min(1.0,sv+pr))
  return mu_u,sv_up
 mu1,svu1=calc_env(d1_lo,d1_hi,nc_1)
 mu2,svu2=calc_env(d2_lo,d2_hi,nc_2)
 l2h=[camp["live_2hop"][0]];l1h=[];lct=[camp["live_ctrls"][2]]
 for sp,qp,y0p,r0p,sc1,qc1,eM1,esv1,sc2,qc2,eM2,esv2 in l2h:
  M0,_,chs1=eval_ball(sp,qp,y0p,r0p,nc_0,stop_first=True)
  assert M0==1 and chs1 and chs1[0][0]==sc1 and chs1[0][1]==qc1
  _,__,y01,r01=chs1[0]
  M1,sv1,chs2=eval_ball(sc1,qc1,y01,r01,nc_1,stop_first=False)
  assert M1==eM1 and sv1==esv1 and sv1<svu1[sc1] and chs2 and chs2[0][0]==sc2 and chs2[0][1]==qc2
  _,__,y02,r02=chs2[0]
  M2,sv2,_=eval_ball(sc2,qc2,y02,r02,nc_2,stop_first=False)
  assert M2==eM2 and sv2==esv2 and sv2<svu2[sc2]
 for sp,qp,y0p,r0p,sc1,qc1,eM1,esv1 in l1h:
  M0,_,chs1=eval_ball(sp,qp,y0p,r0p,nc_0,stop_first=True)
  assert M0==1 and chs1 and chs1[0][0]==sc1 and chs1[0][1]==qc1
  _,__,y01,r01=chs1[0]
  M1,sv1,_=eval_ball(sc1,qc1,y01,r01,nc_1,stop_first=False)
  assert M1==eM1 and sv1==esv1 and sv1<svu1[sc1]
 for sc,qc,y0c,r0c,eM,esv in lct:
  M,sv,_=eval_ball(sc,qc,y0c,r0c,nc_1,stop_first=False)
  assert M==eM and sv==esv and sv<svu1[sc]
 def _stats(rows):
  n=len(rows)
  mus=[mu2[r[1]] if r[0]==2 else mu1[r[1]] for r in rows]
  svus=[svu2[r[1]] if r[0]==2 else svu1[r[1]] for r in rows]
  ms=[r[2] for r in rows];sv2s=[r[3] for r in rows]
  nz=sum(1 for m in ms if m==0)
  nb1_p0=sum((1.0+u)**(-1.0) for u in mus)/n
  rats=[m/u for m,u in zip(ms,mus)]
  m_rat=sum(rats)/n
  v_rat=sum((x-m_rat)**2 for x in rats)/(n-1)
  nb1_v=(m_rat**2)+m_rat*sum(1.0/u for u in mus)/n
  sv2_rat=sum(sv2s)/sum(svus)
  sv2_max=max(s/su for s,su in zip(sv2s,svus))
  pgf10,pgf15={},{}
  for t in (0.0,0.2,0.5,0.8):
   eg=sum(t**m for m in ms)/n
   m10=sum((1.0+u*(1.0-t))**(-1.0) for u in mus)/n
   m15=sum((1.0+(u/3.0)*(1.0-t))**(-1.5) for u in mus)/n
   assert eg<m10 and eg<m15
   pgf10[str(t)]=round(eg/m10,6);pgf15[str(t)]=round(eg/m15,6)
  assert nz==0 and m_rat>=1.0 and sv2_rat<1.0 and sv2_max<1.0
  return dict(
   n=n,M_min=min(ms),M_max=max(ms),n_zero=nz,nb1_p0=round(nb1_p0,5),
   M_over_mu_env=round(m_rat,4),var_M_over_mu=round(v_rat,4),nb1_var_ref=round(nb1_v,3),
   sv2_over_upper=round(sv2_rat,4),sv2_max_ratio=round(sv2_max,4),
   pgf_emp_over_nb10=pgf10,pgf_emp_over_nb15=pgf15,
  )
 all_recs=camp["all_records"];assert len(all_recs)==224
 s_all=_stats(all_recs)
 s_d1=_stats([r for r in all_recs if r[0]==1])
 s_d2=_stats([r for r in all_recs if r[0]==2])
 s_ctrl=_stats([r for r in all_recs if r[0]==0])
 cp95_ub=round(1.0-(0.05**(1.0/len(all_recs))),5)
 d1_vs_ctrl=s_d1["M_over_mu_env"]/s_ctrl["M_over_mu_env"]
 assert cp95_ub<s_all["nb1_p0"] and d1_vs_ctrl>=0.95 and s_d2["M_over_mu_env"]>=1.0
 assert s_all["var_M_over_mu"]<s_all["nb1_var_ref"] and s_d1["var_M_over_mu"]<s_d1["nb1_var_ref"] and s_d2["var_M_over_mu"]<s_d2["nb1_var_ref"]
 return dict(
  tower=camp["tower"],same_field_q15=True,real_2core_fb=True,literal_k1_first_usable=True,
  b_fb_bits=b_fb_b,raw_rels_2core=[len(rels),len(cur)],fb_2core_sizes=[len(fb_g),len(fb_f)],
  ball_pts=n_pts,cp95_empty_ub=cp95_ub,
  mu_u_d1=[round(mu1[0],2),round(mu1[1],2)],mu_u_d2=[round(mu2[0],2),round(mu2[1],2)],
  sv_up_d1=[round(svu1[0],1),round(svu1[1],1)],sv_up_d2=[round(svu2[0],1),round(svu2[1],1)],
  live_2hop_reconstructed=len(l2h),live_1hop_reconstructed=len(l2h)+len(l1h),live_ctrls_replayed=len(lct),
  live_balls_total=3*len(l2h)+2*len(l1h)+len(lct),
  exact_live_match=True,d1_vs_ctrl=round(d1_vs_ctrl,4),
  all_224=s_all,d1_100=s_d1,d2_24=s_d2,ctrl_100=s_ctrl,checks_ok=True,
 )
def verify_prod_sq24(sq24,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn):
 import random
 assert len(sq24)==24
 log2_u8=float(desc["params"]["log2_u8"])
 lx=float(desc["params"]["lx"])
 qbins,qbins_lo=desc["qbins"],desc["qbins_lo"]
 qb_f,qb_g=[x[1] for x in qbins],[x[2] for x in qbins]
 S8=math.log2(desc["cells_rank2"]["8"])
 w_el=y_pow_fn(5);wp1_el=add_fn(w_el,[1]+[0]*7);w2_el=y_pow_fn(10)
 rng=random.Random(20261120)
 ab_pts=[]
 while len(ab_pts)<16:
  al=[rng.randint(-4,4) for _ in range(8)]
  bt=[rng.randint(-4,4) for _ in range(8)]
  if any(al) and sum(x*x for x in al)<=64 and sum(x*x for x in bt)<=64:ab_pts.append((al,bt))
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
  cfs,cgs,sv_rats=[],[],[]
  recs=[r for r in sq24 if r[0]==s_id]
  assert len(recs)==12
  for side,q_hex,v1,v2 in recs:
   q=int(q_hex,16)
   nb=q.bit_length()-1
   assert 70<=nb<265 and q%15==1 and _is_p(q)
   a1,b1,a2,b2=v1[:8],v1[8:],v2[:8],v2[8:]
   det_q=abs(norm_fn(sub_fn(mul_fn(a1,b2),mul_fn(b1,a2))))
   assert det_q>0 and det_q%q==0
   c_val=(math.log2(q)+S8)/16.0
   prow=next(r for r in dlive["rows"] if r["side"]==sname and r["lo"]<=nb < r["hi"])
   S_bin=math.log2(desc["cells_rank2"][str(prow["R"])])
   c_lo_b=(prow["lo"]+S_bin)/16.0
   T_b=prow["kA"]*prow["nA"]+prow["kB"]*prow["n"]+96.0
   sv_plan=0.0
   for (w_q,_,_),(cfl,cgl) in zip(qbins,qbins_lo):
    nfl=72.0*c_lo_b+cfl-(prow["hi"] if sname=="f" else 0.0)
    ngl=24.0*c_lo_b+log2_u8+cgl-(prow["hi"] if sname=="g" else 0.0)
    sf=1.0 if T_b>=nfl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfl)-T_b)/lx))
    sg=1.0 if T_b>=ngl else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(ngl)-T_b)/lx))
    sv_plan+=w_q*sf*sg
   sv_acc=0.0
   c_act=(math.log2(q)+S_bin)/16.0
   lq=math.log2(q)
   for al,bt in ab_pts:
    a=add_fn(mul_fn(al,a1),mul_fn(bt,a2))
    b=add_fn(mul_fn(al,b1),mul_fn(bt,b2))
    a3,b3=pow_fn(a,3),pow_fn(b,3)
    ng=abs(norm_fn(sub_fn(a3,scale_fn(U_bls,mul_fn(w2_el,b3)))))
    a6,b6=mul_fn(a3,a3),mul_fn(b3,b3)
    a9,b9=mul_fn(a6,a3),mul_fn(b6,b3)
    f_el=sub_fn(add_fn(add_fn(a9,mul_fn(wp1_el,mul_fn(a6,b3))),mul_fn(a3,b6)),mul_fn(w_el,b9))
    nf=abs(norm_fn(f_el))
    assert (ng if s_id==0 else nf)%q==0
    cf=math.log2(nf)-72.0*c_val
    cg=math.log2(ng)-24.0*c_val-log2_u8
    cfs.append(cf);cgs.append(cg)
    nfa=72.0*c_act+cf-(lq if sname=="f" else 0.0)
    nga=24.0*c_act+log2_u8+cg-(lq if sname=="g" else 0.0)
    sfe=1.0 if T_b>=nfa else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nfa)-T_b)/lx))
    sge=1.0 if T_b>=nga else min(1.0,max(1.0,T_b/lx)*rho_exact((math.floor(nga)-T_b)/lx))
    sv_acc+=sfe*sge
   sv_rats.append((sv_acc/len(ab_pts))/max(sv_plan,1e-300))
  cfs.sort();cgs.sort();n=len(cfs)
  sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
  sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
  assert all(sf6[i]<qb_f[i] for i in range(6)) and all(sg6[i]<qb_g[i] for i in range(6))
  assert max(sv_rats)<=1.0+1e-12 and (sum(sv_rats)/len(sv_rats))<0.50
  out_sides[sname]=dict(
   bases=12,samples=n,sext_f=sf6,sext_g=sg6,
   min_headroom_f=round(min(qb_f[i]-sf6[i] for i in range(6)),2),
   min_headroom_g=round(min(qb_g[i]-sg6[i] for i in range(6)),2),
   sv_emp_over_plan_mean=round(sum(sv_rats)/len(sv_rats),4),
   sv_emp_over_plan_max=round(max(sv_rats),4),
  )
 return dict(n_bases=24,q_bit_range="73..258",g_side=out_sides["g"],f_side=out_sides["f"],checks_ok=True)
def verify_gt_boot_20(data,P_bls,R_bls,U_bls,h_pi,w_pi,pmul_fn,prem_fn,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn):
 import hashlib
 bt_cfg=data["descent"]["boot"]
 S_boot=int(bt_cfg["s_boot"])
 assert S_boot==(1<<21)
 w_pi2=(w_pi*w_pi)%P_bls
 uw2=(U_bls*w_pi2)%P_bls
 h0,h1,h2,h3=h_pi[0],h_pi[1],h_pi[2],h_pi[3]
 def f4_add(a,b):return [(a[0]+b[0])%P_bls,(a[1]+b[1])%P_bls,(a[2]+b[2])%P_bls,(a[3]+b[3])%P_bls]
 def f4_sub(a,b):return [(a[0]-b[0])%P_bls,(a[1]-b[1])%P_bls,(a[2]-b[2])%P_bls,(a[3]-b[3])%P_bls]
 def f4_scale(s,a):return [(s*a[0])%P_bls,(s*a[1])%P_bls,(s*a[2])%P_bls,(s*a[3])%P_bls]
 def f4_mul(a,b):
  c=[a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[1]*b[1]+a[2]*b[0],
     a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0],
     a[1]*b[3]+a[2]*b[2]+a[3]*b[1],a[2]*b[3]+a[3]*b[2],a[3]*b[3]]
  for k in (6,5,4):
   v=c[k]%P_bls
   if v:
    o=k-4
    c[o]-=v*h0;c[o+1]-=v*h1;c[o+2]-=v*h2;c[o+3]-=v*h3
  return [c[0]%P_bls,c[1]%P_bls,c[2]%P_bls,c[3]%P_bls]
 def f4_pow(a,e):
  res,b=[1,0,0,0],list(a)
  while e:
   if e&1:res=f4_mul(res,b)
   e>>=1
   if e:b=f4_mul(b,b)
  return res
 f1c=[f4_pow([1 if j==i else 0 for j in range(4)],P_bls) for i in range(4)]
 def f4_frob1(a):
  out=[0,0,0,0]
  for i in range(4):
   ai=a[i]
   if ai:
    col=f1c[i]
    out[0]+=ai*col[0];out[1]+=ai*col[1];out[2]+=ai*col[2];out[3]+=ai*col[3]
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
 G0=tuple([int(x,16) for x in c] for c in bt_cfg["G0"])
 one12=([1,0,0,0],[0,0,0,0],[0,0,0,0])
 G0_win=[one12,G0]
 for _ in range(14):G0_win.append(f12_mul(G0_win[-1],G0))
 def f12_pow_G0(e):
  res=one12
  for sh in range(((e.bit_length()-1)//4)*4,-1,-4):
   if res is not one12:res=f12_mul(res,res);res=f12_mul(res,res);res=f12_mul(res,res);res=f12_mul(res,res)
   nib=(e>>sh)&15
   if nib:res=G0_win[nib] if res is one12 else f12_mul(res,G0_win[nib])
  return res
 assert G0!=one12 and f12_pow_G0(R_bls)==one12
 inv2=(P_bls+1)//2
 d_fp2=f4_sub([0,1,0,0],f4_frob2([0,1,0,0]))
 inv_2d=f4_inv(f4_scale(2,d_fp2))
 def to_fp2(a):
  ap=f4_frob2(a)
  return f4_scale(inv2,f4_add(a,ap)),f4_mul(f4_sub(a,ap),inv_2d)
 def elim_gamma(Z):
  e0,o0=to_fp2(Z[0]);e1,o1=to_fp2(Z[1]);e2,o2=to_fp2(Z[2])
  c0=f4_sub(f4_mul(e1,o0),f4_mul(o1,e0))
  c1=f4_sub(f4_mul(o2,e0),f4_mul(e2,o0))
  c2=f4_sub(f4_mul(e2,o1),f4_mul(o2,e1))
  assert f4_frob2(c0)==c0 and f4_frob2(c1)==c1 and f4_frob2(c2)==c2
  W=f12_mul((c0,c1,c2),Z)
  assert W[2]==[0,0,0,0] and any(W[1])
  gam=f4_mul([(-x)%P_bls for x in W[0]],f4_inv(W[1]))
  assert f4_add(W[0],f4_mul(gam,W[1]))==[0,0,0,0]
  return gam
 z3_2=y_pow_fn(10);p4=P_bls**4
 lov_l2=math.log2(((0.99-0.25)**(-15.0/4.0))*(float(P_bls)**0.25)*(float(S_boot)**0.5))
 bsub=data["descent"]["boot_sub_deciles"]
 m_sv=float(bt_cfg["m_sv"])
 betas=[[0]*8]+[[1 if j==i else 0 for j in range(8)] for i in range(8)]+[[-1 if j==i else 0 for j in range(8)] for i in range(8)]
 # Exact Kronecker Z[zeta_15] mul & norm for 105-bit boot vectors (bit-exact integer math, 20x faster!)
 _p15=(1,-1,0,1,-1,1,0,-1,1);_yp=[];_cc=[1,0,0,0,0,0,0,0]
 for _ in range(15):
  _yp.append(tuple(_cc));_hi=_cc[7]
  _cc=[-_hi*_p15[0]]+[_cc[i-1]-_hi*_p15[i] for i in range(1,8)]
 def _cm(k):
  rs=[[] for _ in range(8)]
  for j in range(8):
   yk=_yp[(j*k)%15]
   for i in range(8):
    if yk[i]:rs[i].append((j,yk[i]))
  return rs
 _C14,_C4,_C7=_cm(14),_cm(4),_cm(7)
 def _kmul(a,b):
  S=max(max(abs(x).bit_length() for x in a),max(abs(x).bit_length() for x in b))*2+12
  M=1<<S;H=M>>1;mk=M-1
  A=a[0]+(a[1]<<S)+(a[2]<<(2*S))+(a[3]<<(3*S))+(a[4]<<(4*S))+(a[5]<<(5*S))+(a[6]<<(6*S))+(a[7]<<(7*S))
  B=b[0]+(b[1]<<S)+(b[2]<<(2*S))+(b[3]<<(3*S))+(b[4]<<(4*S))+(b[5]<<(5*S))+(b[6]<<(6*S))+(b[7]<<(7*S))
  P=A*B;c=[0]*15
  for i in range(15):
   v=P&mk
   if v>=H:v-=M;P+=M
   c[i]=v;P>>=S
  for i in range(14,7,-1):
   v=c[i]
   if v:c[i-8]-=v;c[i-7]+=v;c[i-5]-=v;c[i-4]+=v;c[i-3]-=v;c[i-1]+=v
  return c[:8]
 def _knorm(a):
  p1=_kmul(a,[sum(a[j] if s==1 else -a[j] for j,s in r) for r in _C14])
  p2=_kmul(p1,[sum(p1[j] if s==1 else -p1[j] for j,s in r) for r in _C4])
  p3=_kmul(p2,[sum(p2[j] if s==1 else -p2[j] for j,s in r) for r in _C7])
  return abs(p3[0])
 def _kgnorm(a,b):
  a3=_kmul(_kmul(a,a),a);b3=_kmul(_kmul(b,b),b)
  ub3=_kmul(z3_2,[U_bls*x for x in b3])
  return _knorm([a3[j]-ub3[j] for j in range(8)])
 sv_bits,all_bits,sk_l2s=[],[],[]
 decile_pass=gcd_gt1e4=tot_pairs=0
 n_boot_live=8
 for ti in range(n_boot_live):
  xt=int(hashlib.sha256(f"bls12-381-gt-boot-v5g-{ti}".encode()).hexdigest(),16)%(R_bls-1)+1
  gam=elim_gamma(f12_pow_G0(xt))
  R_sk=run_boot_381bit_lll_live(gam,h_pi,P_bls,S_boot,pmul_fn,prem_fn)
  a1,b1=R_sk[0][:8],[x//S_boot for x in R_sk[0][8:]]
  sk_l2=math.log2(math.sqrt(sum(x*x for x in a1)+sum((S_boot*x)**2 for x in b1)))
  assert sk_l2<=lov_l2;sk_l2s.append(sk_l2)
  assert prem_fn(a1,h_pi,P_bls)==pmul_fn(gam,prem_fn(b1,h_pi,P_bls),P_bls,h_pi)
  g1=_kgnorm(a1,b1)
  if ti==0:assert g1==abs(norm_fn(sub_fn(pow_fn(a1,3),scale_fn(U_bls,mul_fn(z3_2,pow_fn(b1,3))))))
  lb1=math.log2(g1);sv_bits.append(lb1)
  if ti<4:
   a2=b2=None
   for idx in range(1,16):
    ak,bk=R_sk[idx][:8],[x//S_boot for x in R_sk[idx][8:]]
    if _knorm(sub_fn(_kmul(a1,bk),_kmul(b1,ak)))==p4:
     a2,b2=ak,bk;break
   assert a2 is not None and prem_fn(a2,h_pi,P_bls)==pmul_fn(gam,prem_fn(b2,h_pi,P_bls),P_bls,h_pi)
   bn,bb=[g1],[lb1]
   for bt in betas[1:]:
    aa=add_fn(a1,_kmul(bt,a2));bvec=add_fn(b1,_kmul(bt,b2))
    ng=_kgnorm(aa,bvec)
    lb=math.log2(ng);assert lb<2600.0
    bn.append(ng);bb.append(lb)
   all_bits.extend(bb)
   bs_s=sorted(bb)
   if all(bs_s[min(16,math.ceil((j+1)*17/10)-1)]<=bsub[j] for j in range(10)):decile_pass+=1
   for i in range(len(bn)):
    for j in range(i+1,len(bn)):
     tot_pairs+=1
     if math.gcd(bn[i],bn[j])>10000:gcd_gt1e4+=1
 n_pass=sum(1 for v in sv_bits if v<=m_sv)
 assert n_pass==n_boot_live and decile_pass==4
 cp_lb95=round(0.05**(1.0/len(sv_bits)),4)
 assert cp_lb95>float(data["descent"]["params"]["ret_boot"])
 return dict(
  boot_element="a - b*X",rel_norm_formula="a^3 - u*omega^2*b^3",
  boot_gt_order_r_exact=True,boot_subfield_elim_live=n_boot_live,
  boot_s_boot=S_boot,boot_live_381bit_lll_exact_all=True,boot_supplied_lll_vectors=0,
  boot_skew_b1_max_l2=round(max(sk_l2s),2),boot_lovasz_bound_l2=round(lov_l2,2),
  boot_gt_pass=f"{n_pass}/{n_boot_live}",boot_cp95_lower=cp_lb95,
  boot_det_p4_exact=True,boot_decile_filter_pass=f"{decile_pass}/4",
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
 assert rk_H==D-1 and rk_OA==D-1 and rk_Cd==D and rk_OA_bad<D-1 and len(ratios)==1
 return dict(
  p=p,r=r,raw_rels=len(rels),purged_rels=len(cur),N_2core=N_2core,n_merged_cols=len(merge_stack),N=N,D=D,
  rank_M=rk_M,rank_A=rk_A,ker_dim_A=len(ker_A),ker_dim_3sm_ctrl=len(ker_3sm),
  rank_Cd=rk_Cd,rank_OeA=rk_OA,rank_H=rk_H,rank_OeA_1block_ctrl=rk_OA_bad,
  unplanted_2core_backsub_exact=True,horner_all_purged_rels_ok=True,distinct_Fp_dlp_ratios=len(ratios),
  checks_ok=True,
 )
