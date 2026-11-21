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
_RHO_C_REV=tuple(tuple(reversed(r))if r is not None else None for r in RHO_C)
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
 M1=1<<k1m
 M2=1<<k2m
 assert M1<=(1<<53)and M2<=(1<<53)and(K1*K2)*(1<<128)<NTT4_PROD
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
 PT=sum((1<<(L-j))*(fft_F(2*n(j-1))+3.0*2*n(j-1))for j in range(1,L+1))
 RT=sum((1<<(L-j+1))*(fft_F(n(j)+G+n(j-1))+3.0*(n(j-1)+G))for j in range(1,L+1))
 return(PT+RT+6.0*fft_F(2*(max(SP,n(L))+G)))/(1<<L)+(1<<15)
def _get_rho_grid():
 return _RHO_S,(_RHO_UMAX-1)*_RHO_S,_RHO_TAIL
def _dw_int(B_in,dn,dd):
 B=[list(r)for r in B_in]
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
def dsieve16_wc_pc(R,cells,lx=70.0,c_leaf=340.0,c_node=196.0,bnd=1800.0,scan_zero=6.0,gate_frac=1.0,gate_frac_hi=None,incl_small150=True,incl_pp=True,c_node_base=172.0,c_center_step=6.0,k1_rej=260.0):
 alpha=100.0/74.0;Re=R+0.05;steps=60;db=(lx-8.0)/steps;tot_ir=0.0
 gfh=gate_frac if gate_frac_hi is None else gate_frac_hi
 vb=[1.0]+[(math.pi**(d/2.0)/math.gamma(d/2.0+1.0))*(Re**d)for d in range(1,17)]
 def _wc_det_ir(det):
  c=(det**(1.0/16.0))*(alpha**(60.0/16.0))
  inv=[1.0/(c*(alpha**(-j/2.0)))for j in range(16)];dp=[1.0];s_ir=140.0
  for k in range(15,0,-1):
   u=inv[k];ndp=[0.0]*(len(dp)+1)
   for s,val in enumerate(dp):ndp[s]+=val;ndp[s+1]+=val*u
   dp=ndp;nk=0.5*(sum(vb[s]*dp[s]for s in range(len(dp)))+1.0)
   s_ir+=(c_node_base+c_center_step*k+(k1_rej if k==1 else 0.0))*nk
  return s_ir
 for i in range(steps):
  b=8.0+(i+0.5)*db;p=2.0**b;wt=2.0*1.10*(p/b)*db;gf=gate_frac if b<48.0 else gfh
  tot_ir+=wt*gf*_wc_det_ir(p)
 hits_pc=2.0*1.10*math.log(lx/8.0)
 if incl_small150:
  for Np,cnt in SMALL_150:
   hits_pc+=cnt/float(Np);tot_ir+=cnt*_wc_det_ir(float(Np))
 if incl_pp:
  lim_pp=2.0**lx
  for Np,cnt in SMALL_150:
   pe=float(Np)*float(Np)
   while pe<=lim_pp:
    hits_pc+=cnt/pe;tot_ir+=cnt*_wc_det_ir(pe);pe*=float(Np)
  for e in range(2,9):
   b_hi=lx/e
   if b_hi<=8.0:break
   st_pp=20;db_pp=(b_hi-8.0)/st_pp
   for i in range(st_pp):
    b=8.0+(i+0.5)*db_pp;p=2.0**b;det=p**e;wt=2.0*1.10*(p/b)*db_pp
    hits_pc+=wt/det;tot_ir+=wt*_wc_det_ir(det)
 return c_leaf*hits_pc+(tot_ir/cells)+bnd+scan_zero
def eval_descent_live(desc,mesh_node=32,mesh_boot=185):
 dp=desc["params"];lx=float(dp["lx"]);assert lx==71.0
 lf_peel=float(dp.get("peel_lf",71.29));lg_peel=float(dp.get("peel_lg",71.69));assert lx<lf_peel<lg_peel
 lrels_peel=float(dp.get("peel_lrels",63.441700678096815))
 uf_peel=float(dp.get("peel_uf",548.9/71.29));ug_peel=float(dp.get("peel_ug",677.8/71.69))
 exc_peel=float(dp.get("peel_exc",1.09))
 peel_h,exc_crit=_peel_fixed_point_bls(lf_peel,lg_peel,lrels_peel,uf_peel,ug_peel,exc_ref=exc_peel)
 alpha_f,alpha_g=float(dp["alpha_f"]),float(dp["alpha_g"])
 log2_u8=float(dp["log2_u8"]);smult=float(dp["sieve_mult"])
 qbins,qbins_lo,qbins_sub=desc["qbins"],desc["qbins_lo"],desc["qbins_sub"]
 cells_r2={int(k):int(v)for k,v in desc["cells_rank2"].items()}
 BOUND_CELL=2400.0*smult;HITS_PER_CELL=2.0*math.log(lx/8.0)
 SIEVE_HIT,SCAN,ZERO=1100.0*smult,10.0*smult,2.0*smult
 PAIR_SETUP=float(dp.get("pair_setup",(2.0**25.1)*(smult*1.33)))
 LLL_NODE=float(dp["lll_node"]);NORMS,TDIV,IDEAL=2.0**27,2.0**25.0,2.0**35
 pairs=2.0*(2.0**log2_li(lx))+sum(2.0*(2.0**log2_li(lx/e))for e in range(2,9))+150.0*18.0
 cfmax=max(x[1]for x in qbins);cgmax=max(x[2]for x in qbins)
 wr_cache={};RINV,RMAX_K,r_arr=_get_rho_grid()
 def get_wr(n,iA,mesh,kA,kB):
  key=(n,iA,mesh,kA,kB)
  if key in wr_cache:return wr_cache[key]
  h=(n-lx)/mesh
  wA=[math.log1p(h/(lx+i*h))if i<iA else 0.0 for i in range(mesh)]
  wB=[math.log1p(h/(lx+i*h))if i>=iA else 0.0 for i in range(mesh)]
  _sp=math.sumprod
  def conv(a,b):
   na,nb=len(a),len(b);out=[0.0]*(na+nb-1)
   i0=next((i for i,x in enumerate(a)if x),na)
   if i0==na:return out
   i1=na-next(i for i,x in enumerate(reversed(a))if x)
   j0=next((j for j,y in enumerate(b)if y),nb)
   if j0==nb:return out
   j1=nb-next(j for j,y in enumerate(reversed(b))if y)
   aa=a[i0:i1];bb_rev=b[j0:j1][::-1];la,lb=len(aa),len(bb_rev)
   for s in range(la+lb-1):
    a_lo=0 if s<lb else s-lb+1;a_hi=s+1 if s+1<la else la
    out[i0+j0+s]=_sp(aa[a_lo:a_hi],bb_rev[lb-1-s+a_lo:lb-1-s+a_hi])
   return out
  pA_list=[[1.0]]
  for j in range(1,kA+1):pA_list.append([x/j for x in conv(pA_list[-1],wA)])
  pAB_list=[conv(pA,wB)for pA in pA_list]if kB else None
  wr_cache[key]=(h,pA_list,pAB_list);return wr_cache[key]
 ss_1_cache={}
 def ss_eval(m_list,n,iA,mesh,alpha,kA,kB):
  if n<=lx:
   out=[]
   for m_raw in m_list:
    m=float(math.ceil(m_raw));u_val=(m+alpha)/lx;rv=rho_exact(u_val)
    if rv>0.0 and u_val>2.0:rv*=max(0.0,1.0-GAMMA*rho_exact(u_val-1.0)/(rv*m*LN2))
    out.append(rv)
   return out
  h,pA_list,pAB_list=get_wr(n,iA,mesh,kA,kB)
  hk=h*RINV/lx;inv_lx_r=RINV/lx;g_ln2=GAMMA/LN2;two_r=2*RINV;lim16=16*RINV
  m_cache=ss_1_cache.setdefault((n,iA,mesh,alpha,kA,kB),{})
  out=[]
  for m_raw in m_list:
   mi=math.ceil(m_raw);cached=m_cache.get(mi)
   if cached is not None:out.append(cached);continue
   m=float(mi);tot=0.0
   for jA in range(kA+1):
    rows_j=((0,pA_list[jA]),(1,pAB_list[jA]))if kB else((0,pA_list[jA]),)
    for jB,row in rows_j:
     j=jA+jB;ifeas=min(len(row)-1,int(math.floor((m-j*h-j*lx)/h)))if h>0 else len(row)-1
     if ifeas<0:continue
     k0=(m-j*lx+alpha)*inv_lx_r;idx_lo=max(0,int((k0-lim16)/hk));mj_lx=m-j*lx
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
   m_cache[mi]=tot;out.append(tot)
  return out
 p_unrec=float(dp.get("p_unrec_leaf",0.0));miss_cache={}
 def miss_side(m_val,u_m,lb_s):
  if m_val<=0.0:return 0.0
  u4=math.ceil(4.0*m_val/lx);key=(u4,u_m,lb_s)
  if key in miss_cache:return miss_cache[key]
  u=u4/4.0;r_u=rho_exact(u);r_m=rho_exact(u_m);tot=0.0
  for b10 in range(400,int(lx*10)):
   b=b10/10.0;c=b+0.1
   km_lo=rho_exact(u_m-b/lb_s)/r_m if r_m>0 else 0.0
   lam=(2.0**(lrels_peel-c))*12.0*km_lo
   ku_hi=rho_exact(u-c/lx)/r_u if r_u>0 else 0.0
   if lam<100.0:x=lam*(1.0-peel_h);p_miss=math.exp(-x)*(1.0+x)
   else:p_miss=math.exp(-lam)
   tot+=ku_hi*math.log(c/b)*min(1.0,p_miss+p_unrec)
  miss_cache[key]=tot;return tot
 def _gb_wave(pbits,leaf):
  SP=math.ceil(pbits/64.0);L=math.ceil(math.log2((SP+1)/leaf));G=L+2;n=lambda j:leaf*(1<<j)
  PT=sum((1<<(L-j))*(fft_F(2*n(j-1))+3.0*2*n(j-1))for j in range(1,L+1))
  RT=sum((1<<(L-j+1))*(fft_F(n(j)+G+n(j-1))+3.0*(n(j-1)+G))for j in range(1,L+1))
  wt=PT+RT+6.0*fft_F(2*(max(SP,n(L))+G));return wt/(1<<L)+(1<<15),wt
 p71_bits=1.000028*(2.0**lx)/LN2;_gwc={};_wpc={}
 def _strip_wave(nb):
  L=max(5,math.ceil(nb/64.0))
  if L not in _gwc:_gwc[L]=_gb_wave(p71_bits,L)
  pl,wt=_gwc[L];return pl+12.0*MONT5*((L/5.0)**2)+2.0*nb*(12.0*L+16.0),wt
 assert desc.get("mode")=="oneshot_boole_gateG"
 r_sh=float(dp["nb_r_shape"]);assert r_sh in (1.5,2.0) and int(dp.get("K_boot",1))==1
 W_tr,P_tr,N_tr={},{},{};bins_seen=[];out_rows=[];scale_iA=mesh_node//32;p_ecm=float(dp["ecm_success"])
 def _knap(sd_k,m_cap,nA_k,nB_k,kA_k,kB_k,vd):
  if kA_k==0 and kB_k==0:return 0.0
  MC=int(math.floor(m_cap))
  def _par(lim):
   its=sorted((int(b[0]),vd[b,sd_k])for b in bins_seen if b[0]<lim-1e-9 and(b,sd_k)in vd)
   o=[];mx=-1.0
   for l0,v0 in its:
    if v0>mx:o.append((l0,v0));mx=v0
   return o
  Ab,Bb=_par(nA_k),_par(nB_k);dpk=[0.0]*(MC+1)
  for _ in range(kA_k):
   nd=list(dpk)
   for l0,v0 in Ab:
    for b in range(l0,MC+1):
     vv=dpk[b-l0]+v0
     if vv>nd[b]:nd[b]=vv
   dpk=nd
  for _ in range(kB_k):
   nd=list(dpk)
   for l0,v0 in Bb:
    for b in range(l0,MC+1):
     vv=dpk[b-l0]+v0
     if vv>nd[b]:nd[b]=vv
   dpk=nd
  return dpk[MC]
 for row in desc["plan_struct"]:
  lo,hi,side,n_int,iA32,kA,kB,R_n,thin=row;bn=(lo,hi)
  if not bins_seen or bins_seen[-1]!=bn:bins_seen.append(bn)
  n_c=float(n_int);iA=iA32*scale_iA;cells=cells_r2[R_n];S=math.log2(cells)
  c_hi=(hi+S)/16.0;c_lo=(lo+S)/16.0
  s_f_hi=lo if side=="f" else 0.0;s_g_hi=lo if side=="g" else 0.0
  s_f_lo=hi if side=="f" else 0.0;s_g_lo=hi if side=="g" else 0.0
  nA=lx+iA*(n_c-lx)/mesh_node if n_c>lx else lx
  T=kA*nA+kB*n_c+96.0;T0=math.ceil(kA*nA+kB*n_c)
  nf_list=[72.0*c_hi+cf-s_f_hi for _,cf,_ in qbins]+[72.0*c_hi+cf-s_f_hi for cf,_ in qbins_sub]
  ng_list=[24.0*c_hi+log2_u8+cg-s_g_hi for _,_,cg in qbins]+[24.0*c_hi+log2_u8+cg-s_g_hi for _,cg in qbins_sub]
  pf_v=ss_eval(nf_list,n_c,iA,mesh_node,alpha_f,kA,kB);pg_v=ss_eval(ng_list,n_c,iA,mesh_node,alpha_g,kA,kB)
  pr_c=sum(pf_v[i]*pg_v[i]for i in range(6))/6.0;pr_i=sum(pf_v[i]*pg_v[i]for i in range(6,30))/24.0
  pr=min(pr_i,4.0*pr_c);sv=sv0=sg_sv=0.0
  for(w,_,_),(cfl,cgl)in zip(qbins,qbins_lo):
   nfl=72.0*c_lo+cfl-s_f_lo;ngl=24.0*c_lo+log2_u8+cgl-s_g_lo
   sf=1.0 if T>=nfl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(nfl)-T)/lx))
   sg=1.0 if T>=ngl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(ngl)-T)/lx))
   sv+=w*sf*sg;sg_sv+=w*sg
   sf0=1.0 if T0>=nfl else min(1.0,max(1.0,T0/lx)*rho_exact((math.floor(nfl)-T0)/lx))
   sg0=1.0 if T0>=ngl else min(1.0,max(1.0,T0/lx)*rho_exact((math.floor(ngl)-T0)/lx))
   sv0+=w*sf0*sg0
  nf_max=72.0*c_hi+cfmax-s_f_hi;ng_max=24.0*c_hi+log2_u8+cgmax-s_g_hi
  k_lx_f=int(nf_max//20.0)+1;k_lx_g=int(ng_max//20.0)+1;k_wA_n=kA if nA>lx+1e-9 else 0
  det_n=p_ecm**(k_lx_f+k_lx_g+2*(kA+kB))
  miss_n=miss_side(nf_max,uf_peel,lf_peel)+miss_side(ng_max,ug_peel,lg_peel)
  mu=cells*pr*float(dp["udup_node"])*det_n*math.exp(-miss_n)
  mr_f,mr_g,mr_T0=mr_instr(nf_max),mr_instr(ng_max),mr_instr(max(T0,100.0))
  sf_n,wf_n=_strip_wave(nf_max);sg_n,wg_n=_strip_wave(ng_max)
  pre_gate=NORMS+TDIV+sf_n+sg_n
  post_gate=(k_lx_f*(ecm_instr(lx,nf_max,p_ecm)+mr_f)+k_lx_g*(ecm_instr(lx,ng_max,p_ecm)+mr_g)+2.0*k_wA_n*(ecm_instr(nA,max(T0,100.0),p_ecm)+mr_T0)+mr_f+mr_g+IDEAL)
  bnd_g=float(dp.get("bound_cell_g",dp.get("bound_cell",3600.0)));bnd_f=float(dp.get("bound_cell_f",0.0))
  if dp.get("dsieve_mode")=="wc_lll99_det":
   if R_n not in _wpc:_wpc[R_n]=dsieve16_wc_pc(R_n,cells,lx,float(dp.get("c_leaf",340.0)),float(dp.get("c_node",196.0)),bnd_g,6.0,float(dp.get("gate_frac",1.0)),dp.get("gate_frac_hi"),True,True,float(dp.get("c_node_base",172.0)),float(dp.get("c_center_step",6.0)),float(dp.get("k1_rej",260.0)))
   sieve=cells*_wpc[R_n]+pairs*PAIR_SETUP
  else:sieve=cells*(HITS_PER_CELL*SIEVE_HIT+SCAN+ZERO+BOUND_CELL)+pairs*PAIR_SETUP
  gate_ecm=cells*pre_gate+wf_n+wg_n+(cells*sv0+cells*pr)*post_gate;p_nc=float(dp.get("p_near_cancel",0.0))
  gate_sieve=(sieve+cells*(sg_sv+p_nc)*bnd_f+(cells*(sv+p_nc)+cells*pr)*pre_gate+wf_n+wg_n+(cells*sv0+cells*pr)*post_gate)
  if gate_ecm<gate_sieve:how="ecm_gateG";own_1=LLL_NODE+gate_ecm
  else:how="sieve_gateG";mu*=float(dp["ret_node"]);own_1=LLL_NODE+gate_sieve
  # 1-shot zero-sibling-retry pure Boole union bound: NO retries, NO independence!
  cP=_knap("f",nf_max,nA,n_c,kA,kB,P_tr)+_knap("g",ng_max,nA,n_c,kA,kB,P_tr)
  cW=_knap("f",nf_max,nA,n_c,kA,kB,W_tr)+_knap("g",ng_max,nA,n_c,kA,kB,W_tr)
  cN=_knap("f",nf_max,nA,n_c,kA,kB,N_tr)+_knap("g",ng_max,nA,n_c,kA,kB,N_tr)
  p0=(1.0+mu/r_sh)**(-r_sh)
  assert mu>=30.0 and p0<0.01,f"Node ({lo},{hi},{side}) mu={mu} p0={p0}"
  W_tr[bn,side]=own_1+cW;P_tr[bn,side]=p0+cP;N_tr[bn,side]=1.0+cN
  out_rows.append(dict(lo=lo,hi=hi,side=side,how=how,R=R_n,n=n_int,iA=iA,nA=nA,kA=kA,kB=kB,r_shape=r_sh,pr=pr,mu=mu,p0=p0,P_sub=P_tr[bn,side],N_det=int(round(N_tr[bn,side])),nf=nf_max,ng=ng_max,own1_log2=math.log2(own_1),W_log2=math.log2(W_tr[bn,side])))
 assert len(bins_seen)==39 and bins_seen[0]==(71,76)and bins_seen[-1]==(261,266)
 bt=desc["boot"];assert int(bt["cells_b"])==1 and int(dp["boot_rels"])==1 and float(dp.get("ball_pair",0.0))==0.0
 m_max=float(bt["m_sv"]);iA_b=math.floor((bt["nA"]-lx)*mesh_boot/(bt["n1"]-lx)+1e-12)
 assert lx+iA_b*(bt["n1"]-lx)/mesh_boot<=bt["nA"]+1e-12
 pr_raw_b=ss_eval([m_max],float(bt["n1"]),iA_b,mesh_boot,alpha_g,bt["kA"],1)[0]
 mb=miss_side(m_max,ug_peel,lg_peel);pr_b=pr_raw_b*math.exp(-mb)
 k_lx_b=int(m_max//20.0)+1;kA_b=int(bt["kA"]);mr_b=mr_instr(m_max)
 p_ecm_b=float(dp.get("ecm_success_boot",p_ecm));det_b=p_ecm_b**(k_lx_b+kA_b+1)
 p_boot_trial=pr_b*det_b*float(dp["ret_boot"])
 trial_b=NORMS+TDIV+k_lx_b*(ecm_instr(lx,m_max,p_ecm_b)+mr_b)+kA_b*(ecm_instr(bt["nA"],m_max,p_ecm_b)+mr_b)+mr_b+IDEAL
 LLL_BOOT=float(dp.get("lll_boot",2.0**48));boot_own_1=(1.0/p_boot_trial)*(LLL_BOOT+trial_b)
 P_b_ch=_knap("g",m_max,bt["nA"],float(bt["n1"]),bt["kA"],1,P_tr)
 W_b_ch=_knap("g",m_max,bt["nA"],float(bt["n1"]),bt["kA"],1,W_tr)
 N_b_ch=int(round(_knap("g",m_max,bt["nA"],float(bt["n1"]),bt["kA"],1,N_tr)))
 # Boot draws fresh i.i.d. uniform s in Z_r^* up to T_boot = ceil(20 / p_boot_trial) -> fail <= e^-20 per log, then runs 1 tree (K_boot=1)!
 p_boot_find_fail_2=2.0*math.exp(-20.0)
 P_2log_boole=2.0*P_b_ch
 p_desc_total_loss=P_2log_boole+p_boot_find_fail_2
 two_exp_u=2.0*(boot_own_1+W_b_ch)/20000.0
 SP71=math.ceil(p71_bits/64.0);p71_tree_u=math.log2(math.ceil(math.log2(SP71))*fft_F(SP71)/20000.0)
 p71_owner_bytes=(p71_bits/8.0)/128.0;assert p71_owner_bytes<(2.0**62.0)and p71_tree_u<=74.5
 return dict(mesh_node=mesh_node,mesh_boot=mesh_boot,lx=lx,r_shape=r_sh,n_nodes=len(desc["plan_struct"]),peel_h=peel_h,exc_crit=exc_crit,boot_pr_log2=math.log2(pr_b),boot_miss=mb,p_boot_trial=p_boot_trial,boot_own_exp_log2=math.log2(2.0*boot_own_1/20000.0),tree_det_two_log2=math.log2(2.0*W_b_ch/20000.0),N_tree_1log_det=N_b_ch,N_tree_2log_det=2*N_b_ch,max_p0=max(r["p0"]for r in out_rows),min_mu=min(r["mu"]for r in out_rows),P_1log_boole=P_b_ch,P_2log_boole=P_2log_boole,p_boot_find_fail_2=p_boot_find_fail_2,p_desc_total_loss=p_desc_total_loss,p70_tree_build_u=p71_tree_u,p71_tree_build_u=p71_tree_u,p71_owner_log2_B=math.log2(p71_owner_bytes),two_log_exp_units=math.log2(two_exp_u),rows=out_rows,_ss_eval=ss_eval,)

def verify_kernels_symbolic(kernels_text,coll,la,r_bls,lll_node_cap,lll_boot_cap,desc_p=None):
 import re,random,struct
 allowed={"ldp","stp","ldr","str","ldrb","strb","ldrh","strh","ldur","stur","add","adds","adc","adcs","sub","subs","sbc","sbcs","mul","umulh","smulh","madd","msub","and","orr","eor","lsl","lsr","asr","cmp","tst","csel","mov","b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","b.ge","b.lt","b.le","clz","udiv","cbz","cbnz","ret"}
 forb=re.compile(r"\b(v[0-9]+|q[0-9]+|d[0-9]+|z[0-9]+|dc|ldadd\w*|swp\w*|cas\w*|bl|blr)\b")
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
 _lp={"r19","r20","r21","r22","r18","r3","r4"};_we=_pw(parsed["rev_b2_prod_emit_1"]["insns"])
 assert _we.isdisjoint(_lp)
 assert parsed["rev_sb_probe_fast"]["ir"]<=coll["c_tile_hdr"]and "b.hs .Lempty" in parsed["rev_sb_probe_fast"]["insns"]
 assert parsed["rev_b2_prod_sb_hdr"]["ir"]<=coll["c_sb_hdr"]and "add x22, x22, #256" in parsed["rev_b2_prod_sb_hdr"]["insns"]
 assert parsed["rev_b1a_apply_wrap8"]["ir"]<=4.0 and "subs x_cnt, x_cnt, #8" in parsed["rev_b1a_apply_wrap8"]["insns"] and "adds x_cnt, x_cnt, #8" in parsed["rev_b1a_apply_wrap8"]["insns"]
 assert parsed["rev_b1a_apply_u8"]["ir"]<=coll["c_hit_lo_m4"]and parsed["rev_b1a_apply_tail1"]["ir"]<=6.0
 assert parsed["rev_b1b_apply_u8"]["ir"]<=coll["c_hit_lo_u"]and parsed["rev_b1b_apply_tail1"]["ir"]<=8.0
 assert parsed["rev_b2_prod_emit_1"]["ir"]<=coll["c_b2_emit"]
 assert(parsed["rev_b2_prod_emit_1"]["ir"]+parsed["rev_b1a_apply_u8"]["ir"]+0.10)<=coll["c_hit_hi"]
 _r1=random.Random(2026);_M64=(1<<64)-1
 for _klen in range(26):
  _oids=[_r1.randrange(0,64)for _ in range(_klen)]
  _acc1=[0]*64;_acc2=[0]*64;_lpv=37
  for _o in _oids:_acc1[_o]+=_lpv
  _i=0;_xc=_klen;_ops1a=2;_c_fl=(_xc>=8);_xc=(_xc-8)&_M64
  if _c_fl:
   while True:
    assert _i+8<=_klen,"ZERO overread!"
    for _t in range(8):_acc2[_oids[_i+_t]]+=_lpv
    _i+=8;_ops1a+=30;_c_fl=(_xc>=8);_xc=(_xc-8)&_M64
    if not _c_fl:break
  _xc=(_xc+8)&_M64;_ops1a+=2
  if _xc!=0:
   while True:
    assert _i<_klen,"ZERO overread!"
    _acc2[_oids[_i]]+=_lpv;_i+=1;_ops1a+=6;_xc=(_xc-1)&_M64
    if _xc==0:break
  assert _i==_klen and _acc1==_acc2 and _ops1a<=3.75*_klen+20
  _k57=(_klen&3)<<57;_M59=(1<<59)-1;_pk=[((_r1.randrange(4)<<57)|_o)for _o in _oids]
  _ag=[0]*256;_at=[0]*256
  for _w in _pk:_rt=((_w-_k57)&_M64)&_M59;_ag[(_rt>>57)*64+(_rt&63)]+=_lpv
  _ib=0;_xb=_klen;_ops1b=2;_cb=(_xb>=8);_xb=(_xb-8)&_M64
  if _cb:
   while True:
    assert _ib+8<=_klen,"ZERO Band-1b overread!"
    for _t in range(8):_rt=((_pk[_ib+_t]-_k57)&_M64)&_M59;_at[(_rt>>57)*64+(_rt&63)]+=_lpv
    _ib+=8;_ops1b+=46;_cb=(_xb>=8);_xb=(_xb-8)&_M64
    if not _cb:break
  _xb=(_xb+8)&_M64;_ops1b+=2
  if _xb!=0:
   while True:
    assert _ib<_klen,"ZERO Band-1b tail overread!"
    _rt=((_pk[_ib]-_k57)&_M64)&_M59;_at[(_rt>>57)*64+(_rt&63)]+=_lpv;_ib+=1;_ops1b+=8;_xb=(_xb-1)&_M64
    if _xb==0:break
  assert _ib==_klen and _ag==_at and _ops1b<=5.75*_klen+20<=5.75*_klen+22
 assert parsed["ntt62_butterfly"]["ir"]*12.0<=120.0<168.0
 assert(parsed["rev_build_prod_pass12"]["ir"]+0.50)<=coll["c_list_a"]and "ldr x2, [x_hist, x1, lsl #3]" in parsed["rev_build_prod_pass12"]["insns"]and "str x2, [x_hist, x1, lsl #3]" in parsed["rev_build_prod_pass12"]["insns"]and "cmp x0, x_pms" in parsed["rev_build_prod_pass12"]["insns"]
 assert parsed["rev_build_prod_prefix"]["ir"]<=coll["c_lbkt"]
 MASK=(1<<64)-1;MASK59=(1<<59)-1;OID_M=(1<<57)-1
 for m_b in range(4):
  for k_c in range(4):
   xk57=k_c<<57
   for oid in (0,1,0x12345,OID_M):
    om=(m_b<<57)|oid
    rot=((om-xk57)&MASK)&MASK59
    assert(rot>>57)==((m_b-k_c)&3)and(rot&OID_M)==oid and(rot<<1)<(1<<60)
    if m_b==0 and k_c==1:assert(rot>>57)==3
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
 # Flat-array 2-pass 4-bank build (explicit x_oidm = m<<57 bank reset at V_a=2 << 2^57!) + 1-rec-at-a-time Band-2 multi-SB replay (incl. empty SB=1, odd SB=2 len=1, even SB=0 len=4):
 kas_4b=[[5,7],[31,0],[64+19,12],[9,28]] # 4 banks x V_a=2 entries each (V_a << 2^57!)
 kbuf_flat=[ka for b in kas_4b for ka in b]
 hsb=[0,0,0]
 for ka in kbuf_flat:hsb[ka>>5]+=1
 off_b=[0]
 for c in hsb:off_b.append(off_b[-1]+c*8)
 assert off_b==[0,56,56,64]
 fba=bytearray(off_b[-1]);psb=list(off_b[:-1]);idx_kb=0;r_prod=[]
 for m_b in range(4):
  x_oidm=m_b<<57 # explicit bank-start reset!
  for _ in range(len(kas_4b[m_b])):
   ka=kbuf_flat[idx_kb];idx_kb+=1;si=ka>>5
   rec=((ka&31)<<59)|x_oidm;struct.pack_into("<Q",fba,psb[si],rec);psb[si]+=8
   r_prod.append((ka,x_oidm));x_oidm+=1
 assert r_prod[1][1]==(0<<57)|1 and r_prod[2][1]==(1<<57)|0 and r_prod[4][1]==(2<<57)|0 and r_prod[6][1]==(3<<57)|0
 flat_buf=bytes(fba)
 x19=0;x21=3;x22=0;em=[]
 while x21>0:
  x5,x6=off_b[x19],off_b[x19+1];x19+=1
  if x5!=x6:
   xp,xpe=x5,x6
   while xp<xpe:
    x15,=struct.unpack_from("<Q",flat_buf,xp);xp+=8
    rr=x15>>59;rot=((x15-(1<<57))&MASK)&MASK59
    em.append((x22//256,rr,rot>>57,rot&OID_M))
  x22+=256;x21-=1
 exh=[(ka>>5,ka&31,(((om-(1<<57))&MASK)&MASK59)>>57,om&OID_M)for ka,om in r_prod]
 assert x19==3 and x21==0 and x22==768 and sorted(em)==sorted(exh)and len(em)==8
 assert(1>>5)==0 and(31>>5)==0 and(32>>5)==1 and "add x4, x3, x_sb_m1" not in parsed["rev_sb_probe_fast"]["insns"]
 _sc32=parsed["rev_scan_zero_u32"]
 assert _sc32["items"]==32 and _sc32["raw_opcodes"]==25 and all("], #" not in _i for _i in _sc32["insns"])
 assert _pw(_sc32["insns"][4:7]+_sc32["insns"][12:18]).isdisjoint({"r0","r1","r2","r3","r6","r7","r8","r9"})
 scan_tot=_sc32["ir"]*1.02761+0.12
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
 assert parsed["unit_run5_z12"]["ir"]<=8.80 and parsed["unit_run5_neg_z12"]["ir"]<=8.80
 assert parsed["gen_run1_z12"]["ir"]<=29.0 and parsed["gen_run1_neg_z12"]["ir"]<=29.0 and 4.0*max(parsed["gen_run1_z12"]["ir"],parsed["gen_run1_neg_z12"]["ir"])<=116.0
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
 # Bit-exact 5-limb simulation of unit_run5_z12, unit_run5_neg_z12, gen_run1_z12 (c>0), and gen_run1_neg_z12 (c<0):
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
 c_nb=float(desc_p.get("c_node_base",172.0))if desc_p else 172.0
 c_cs=float(desc_p.get("c_center_step",6.0))if desc_p else 6.0
 c_k1=float(desc_p.get("k1_rej",260.0))if desc_p else 260.0
 c_lf=float(desc_p["c_leaf"])if desc_p else 340.0
 c_bc=float(desc_p["bound_cell"])if desc_p else 3600.0
 c_bcg=float(desc_p.get("bound_cell_g",1800.0))if desc_p else 1800.0
 c_bcf=float(desc_p.get("bound_cell_f",1800.0))if desc_p else 1800.0
 c_ps=float(desc_p["pair_setup"])if desc_p else float(1<<35)
 assert parsed["dsieve16_q32_node"]["ir"]<=55.0 and parsed["dsieve16_dfs_ctrl"]["ir"]<=31.0 and parsed["dsieve16_dfs_cstep"]["ir"]<=c_cs and parsed["dsieve16_vec_madd"]["ir"]==40.0
 assert parsed["dsieve16_q32_node"]["ir"]+parsed["dsieve16_dfs_ctrl"]["ir"]+parsed["dsieve16_vec_madd"]["ir"]+40.0+6.0<=c_nb
 assert parsed["dsieve16_q32_node"]["ir"]+parsed["dsieve16_vec_madd"]["ir"]+12.0+2.0*(parsed["dsieve16_leaf_vec"]["ir"]+32.0)<=c_k1
 assert parsed["dsieve16_leaf_vec"]["ir"]==40.0 and parsed["dsieve16_rank_step"]["ir"]<=10.0 and(parsed["dsieve16_leaf_vec"]["ir"]+32+16+16*parsed["dsieve16_rank_step"]["ir"]+56+12)<=c_lf
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
 for x0_t in [0,1,2,3,4,15,16,16383,16384,16385,(1<<32)-1,1<<32,(1<<46)-1,1<<46]+[rng.randrange(1,1<<46)for _ in range(100)]:
  if x0_t==0:x2_t=1
  else:
   clz_t=64-x0_t.bit_length();x2_t=1<<((65-clz_t)>>1)
   for _ in range(5):x2_t=(x2_t+(x0_t//x2_t))>>1
   x2_t+=1
  assert x2_t>=math.isqrt(x0_t)+1 and x2_t<=math.isqrt(x0_t)+2
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
 s1_ds=math.floor(2100.0/math.log2(10.0/9.0))
 ds_ops=15*(40*5*3+16*12*4+15*12*10)+15*(4*12*3*5+2*40*8*3)+2000
 ds_lll_bound=(2*s1_ds+16)*ds_ops
 assert ds_lll_bound<c_ps
 s1_sq=math.floor(8400.0/math.log2(40.0/11.0))
 sq1_ops=15*(40*15*10+16*12*28+15*12*78)+15*(4*12*10*15+2*40*25*10)+10000
 s2_sq=math.floor((340.0*math.log2(40.0)+7.5*560.0)/math.log2(100.0/99.0))
 sq_lll_bound=(2*s1_sq+16)*sq1_ops+(2*s2_sq+16)*sq1_ops
 assert sq_lll_bound<lll_node_cap
 s1_bt=math.floor(67191.0/math.log2(40.0/11.0))
 bt1_ops=15*(40*199*132+16*12*67*67+15*12*67*199)+15*(4*12*132*199+2*40*331*132)+1000000
 s2_bt=math.floor((340.0*math.log2(40.0)+7.5*3384.0)/math.log2(100.0/99.0))
 bt2_ops=15*(40*80*53+16*12*27*27+15*12*27*80)+15*(4*12*53*80+2*40*133*53)+1000000
 bt_lll_bound=(2*s1_bt+16)*bt1_ops+(2*s2_bt+16)*bt2_ops
 assert bt_lll_bound<lll_boot_cap
 q265=((1<<265)+1234567891)*8191;Bsq=[[0]*16 for _ in range(16)];Bsq[0][0]=q265
 ry=pow(7,13,q265);rt=pow(11,17,q265);ryp=[pow(ry,i,q265)for i in range(8)]
 for i in range(1,8):Bsq[i][0]=(-ryp[i])%q265;Bsq[i][i]=1
 for i in range(8):Bsq[8+i][0]=(rt*ryp[i])%q265;Bsq[8+i][8+i]=1
 _,dsq,lsq=lll_16d_int(Bsq,99,100)
 assert dsq[16]==q265*q265 and q265.bit_length()==278
 for kk in range(1,16):
  assert 100*(dsq[kk-1]*dsq[kk+1]+lsq[kk][kk-1]**2)>=99*dsq[kk]**2
  for jj in range(kk):assert 2*abs(lsq[kk][jj])<=dsq[jj+1]
 return dict(asm_blocks_parsed=len(parsed),b2_multi_entry_cf_ok=True,b2_bank_rot_16pairs_ok=True,uint16_swar_ok=True,ntt62_garner_exact=True,ntt62_3x4p_stage_ir=round(12.0*parsed["ntt62_butterfly"]["ir"],1),f62_barrett_mac_64_exact=True,mont5_opcodes=mont5_ops,mont5_sim_32_exact=True,row_red5_barrett_opcodes=parsed["row_red5_z12"]["raw_opcodes"],row_red5_barrett_47_exact=True,unit_run5_ir=round(parsed["unit_run5_z12"]["ir"],2),gen_run1_4coord_ir=round(4.0*parsed["gen_run1_z12"]["ir"],2),scan_zero_u16=round(scan_tot,4),sq_lll_uniform_log2=round(math.log2(sq_lll_bound),2),boot_lll_uniform_log2=round(math.log2(bt_lll_bound),2),all_raw_opcodes_ok=True)
from unitdomain import verify_prod_sq24
def verify_gt_boot_20(boot_w,data,P_bls,R_bls,U_bls,h_pi,w_pi,pmul_fn,prem_fn,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn):
 import hashlib
 S_boot=int(data["descent"]["boot"]["s_boot"])
 n_targets=int(boot_w["n_targets"])
 assert S_boot==int(boot_w["s_boot"])==(1<<21)and n_targets>=1
 assert "full4" not in boot_w and "extra16" not in boot_w,"ZERO pre-supplied LLL vectors allowed!"
 w_pi2=(w_pi*w_pi)%P_bls;uw2=(U_bls*w_pi2)%P_bls;h0,h1,h2,h3=h_pi[0],h_pi[1],h_pi[2],h_pi[3]
 def f4_add(a,b):return[(a[0]+b[0])%P_bls,(a[1]+b[1])%P_bls,(a[2]+b[2])%P_bls,(a[3]+b[3])%P_bls]
 def f4_sub(a,b):return[(a[0]-b[0])%P_bls,(a[1]-b[1])%P_bls,(a[2]-b[2])%P_bls,(a[3]-b[3])%P_bls]
 def f4_scale(s,a):return[(s*a[0])%P_bls,(s*a[1])%P_bls,(s*a[2])%P_bls,(s*a[3])%P_bls]
 def f4_mul(a,b):
  c=[a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[1]*b[1]+a[2]*b[0],a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0],a[1]*b[3]+a[2]*b[2]+a[3]*b[1],a[2]*b[3]+a[3]*b[2],a[3]*b[3]]
  for k in (6,5,4):
   v=c[k]%P_bls
   if v:o=k-4;c[o]-=v*h0;c[o+1]-=v*h1;c[o+2]-=v*h2;c[o+3]-=v*h3
  return[c[0]%P_bls,c[1]%P_bls,c[2]%P_bls,c[3]%P_bls]
 def f4_pow(a,e):
  res,b=[1,0,0,0],list(a)
  while e:
   if e&1:res=f4_mul(res,b)
   e>>=1
   if e:b=f4_mul(b,b)
  return res
 f1c=[f4_pow([1 if j==i else 0 for j in range(4)],P_bls)for i in range(4)]
 def f4_frob1(a):
  out=[0,0,0,0]
  for i in range(4):
   ai=a[i]
   if ai:col=f1c[i];out[0]+=ai*col[0];out[1]+=ai*col[1];out[2]+=ai*col[2];out[3]+=ai*col[3]
  return[x%P_bls for x in out]
 def f4_frob2(a):return f4_frob1(f4_frob1(a))
 def f4_inv(a):
  a1=f4_frob1(a);a2=f4_frob1(a1);a3=f4_frob1(a2)
  pr3=f4_mul(a1,f4_mul(a2,a3))
  return f4_scale(pow(f4_mul(a,pr3)[0],-1,P_bls),pr3)
 def f12_mul(A,B):
  return(f4_add(f4_mul(A[0],B[0]),f4_scale(uw2,f4_add(f4_mul(A[1],B[2]),f4_mul(A[2],B[1])))),f4_add(f4_add(f4_mul(A[0],B[1]),f4_mul(A[1],B[0])),f4_scale(uw2,f4_mul(A[2],B[2]))),f4_add(f4_add(f4_mul(A[0],B[2]),f4_mul(A[1],B[1])),f4_mul(A[2],B[0])),)
 G0=tuple([int(x,16)for x in c]for c in boot_w["G0"])
 one12=([1,0,0,0],[0,0,0,0],[0,0,0,0])
 G0_win=[one12,G0]
 for _ in range(14):G0_win.append(f12_mul(G0_win[-1],G0))
 def f12_pow_G0(e):
  res=one12
  for sh in range(((e.bit_length()-1)//4)*4,-1,-4):
   if res is not one12:res=f12_mul(res,res);res=f12_mul(res,res);res=f12_mul(res,res);res=f12_mul(res,res)
   nib=(e>>sh)&15
   if nib:res=G0_win[nib]if res is one12 else f12_mul(res,G0_win[nib])
  return res
 assert G0!=one12 and f12_pow_G0(R_bls)==one12
 inv2=(P_bls+1)//2;d_fp2=f4_sub([0,1,0,0],f4_frob2([0,1,0,0]));inv_2d=f4_inv(f4_scale(2,d_fp2))
 def to_fp2(a):ap=f4_frob2(a);return f4_scale(inv2,f4_add(a,ap)),f4_mul(f4_sub(a,ap),inv_2d)
 def elim_gamma(Z):
  e0,o0=to_fp2(Z[0]);e1,o1=to_fp2(Z[1]);e2,o2=to_fp2(Z[2])
  c0=f4_sub(f4_mul(e1,o0),f4_mul(o1,e0));c1=f4_sub(f4_mul(o2,e0),f4_mul(e2,o0));c2=f4_sub(f4_mul(e2,o1),f4_mul(o2,e1))
  assert f4_frob2(c0)==c0 and f4_frob2(c1)==c1 and f4_frob2(c2)==c2
  W=f12_mul((c0,c1,c2),Z)
  assert W[2]==[0,0,0,0]and any(W[1])
  gam=f4_mul([(-x)%P_bls for x in W[0]],f4_inv(W[1]))
  assert f4_add(W[0],f4_mul(gam,W[1]))==[0,0,0,0]
  return gam
 z3_2=y_pow_fn(10);p4=P_bls**4;det_exp=(P_bls**8)*(S_boot**16)
 lov_l2=math.log2(((0.99-0.25)**(-15.0/4.0))*(float(P_bls)**0.25)*(float(S_boot)**0.5))
 bsub=data["descent"]["boot_sub_deciles"];m_sv=float(data["descent"]["boot"]["m_sv"])
 betas=[[0]*8]+[[1 if j==i else 0 for j in range(8)]for i in range(8)]+[[-1 if j==i else 0 for j in range(8)]for i in range(8)]
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
   v2=next(row for row in R_sk[1:]if all(x%S_boot==0 for x in row[8:])and abs(norm_fn(sub_fn(mul_fn(a1,[x//S_boot for x in row[8:]]),mul_fn(b1,row[:8]))))==p4)
   a2,b2=v2[:8],[x//S_boot for x in v2[8:]]
   assert prem_fn(a2,h_pi,P_bls)==pmul_fn(gam,prem_fn(b2,h_pi,P_bls),P_bls,h_pi)
   bn,bb=[g1],[lb1]
   for bt in betas[1:]:
    aa=add_fn(a1,mul_fn(bt,a2));bvec=add_fn(b1,mul_fn(bt,b2))
    ng=abs(norm_fn(sub_fn(pow_fn(aa,3),scale_fn(U_bls,mul_fn(z3_2,pow_fn(bvec,3))))))
    lb=math.log2(ng);assert lb<2600.0
    bn.append(ng);bb.append(lb)
   all_bits.extend(bb);bs_s=sorted(bb)
   if all(bs_s[min(16,math.ceil((j+1)*17/10)-1)]<=bsub[j]for j in range(10)):decile_pass+=1
   for i in range(len(bn)):
    for j in range(i+1,len(bn)):
     tot_pairs+=1
     if math.gcd(bn[i],bn[j])>10000:gcd_gt1e4+=1
 from fractions import Fraction as _F
 _tr=lambda d:8 if math.gcd(15,d%15)==15 else(-4 if math.gcd(15,d%15)==5 else(-2 if math.gcd(15,d%15)==3 else 1))
 _Gm=[[_F(15001,1000)-_tr(0)if i==j else _F(-_tr(i-j))for j in range(8)]for i in range(8)]
 for _k in range(1,9):
  _A=[[x for x in r[:_k]]for r in _Gm[:_k]];_det=_F(1)
  for _i in range(_k):
   _pv=next(r for r in range(_i,_k)if _A[r][_i]!=0)
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
 return dict(boot_element="a - b*X",rel_norm_formula="a^3 - u*omega^2*b^3",boot_gt_order_r_exact=True,boot_subfield_elim_live=n_targets,boot_s_boot=S_boot,boot_live_381bit_lll_exact_all=n_targets,boot_presupplied_vectors=0,boot_skew_b1_max_l2=round(max(sk_l2s),2),boot_lovasz_bound_l2=round(lov_l2,2),boot_det_norm_theorem_ub=round(_det_norm_ub,2),boot_ret_proved=1.0,boot_det_compat_pass=f"{n_pass}/{n_targets}",boot_det_p4_exact=True,boot_decile_filter_pass=f"{decile_pass}/{n_targets}",sv_mean_live=round(sum(sv_bits)/len(sv_bits),2),sv_max_live=round(max(sv_bits),2),ball_mean_live=round(sum(all_bits)/len(all_bits),2),ball_max_live=round(max(all_bits),2),gcd_gt1e4_live=f"{gcd_gt1e4}/{tot_pairs}",)
