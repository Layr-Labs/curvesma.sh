"""cm_models.py — Self-contained."""
import array
import collections
import functools
import math
from decimal import Decimal,localcontext
LN2=math.log(2.0);GAMMA=0.5772156649015329;MONT5=367.76;ECM_TORSION=12.0;ECM_WHEELS=(30,210,2310,30030)
def build_rho_table(kmax=16,nt=120):
 with localcontext() as ctx:
  ctx.prec=80;zero,one=Decimal(0),Decimal(1);coefs=[None,[one]+[zero]*(nt-1)]
  for k in range(2,kmax+1):
   prev,c=coefs[k-1],[zero]*nt
   for m in range(nt-1):
    c[m+1]=(prev[m]+m*c[m])/(k*(m+1))
   c[0]=prev[0]-sum(c[1:]);coefs.append(c)
  return [None]+[[float(x) for x in c[:80]] for c in coefs[1:]]
RHO_C=build_rho_table()
def _build_rho_tail(S=128,umax=45):
 r=[1.0]*(S+1)+[1.0-math.log(1.0+j/S) for j in range(1,S+1)]
 for i in range(2*S+1,umax*S+1):
  u=i/S;num=(r[i-S]+4.0*sum(r[i-S+1:i:2])+2.0*sum(r[i-S+2:i:2]))/(3.0*S);r.append(num/(u-1.0/(3.0*S)))
 return S,umax,r
_RHO_S,_RHO_UMAX,_RHO_TAIL=_build_rho_tail()
def _trim_rc(row):
 if row is None:return None
 thr=1e-18*abs(row[0]);hi=len(row)-1
 while hi>0 and abs(row[hi])<thr:hi-=1
 return tuple(reversed(row[:hi+1]))
_RHO_C_REV=tuple(_trim_rc(row) for row in RHO_C);_RHO_C_LEN=len(_RHO_C_REV)
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
 pos=u*_RHO_S;idx=int(pos);frac=pos-idx;y0,y1=_RHO_TAIL[idx],_RHO_TAIL[idx+1]
 if y0<=0.0 or y1<=0.0:return 0.0
 return math.exp((1.0-frac)*math.log(y0)+frac*math.log(y1))
def rho_db(u,log2_N):
 r0=rho_exact(u);r1=rho_exact(max(0.0,u-1.0))
 return max(1e-300,r0-GAMMA*r1/max(1.0,log2_N*LN2))
def log2_li(x):
 y=x*LN2
 return x-math.log2(y)+math.log2(1.0+1.0/y+2.0/(y*y))
@functools.lru_cache(None)
def primes_upto(n):
 sieve=bytearray(b"\x01")*(n+1);sieve[:2]=b"\x00\x00"
 for p in range(2,math.isqrt(n)+1):
  if sieve[p]:sieve[p*p::p]=b"\x00"*((n-p*p)//p+1)
 return [p for p in range(2,n+1) if sieve[p]]
def _ret_coeffs(weights,normalizer,length):
 _sp=math.sumprod;jw_rev=[j*weights[j] for j in range(length,-1,-1)];out=[math.exp(-normalizer)]+[0.0]*length;L=length
 jm=next((j for j in range(1,L+1) if weights[j]),L+1);e_rev=L-jm+1
 for n in range(jm,L+1):
  val=_sp(jw_rev[L-n:e_rev],out[:n-jm+1]);out[n]=val/n
  if val<-1e-13:raise ArithmeticError("negative coefficient beyond roundoff tolerance")
 return out
_RET_SM_W={};_RET_SM_C={}
class Retention:
 def __init__(self,B_bits,X_bits,*,prime_min=2**24,step=0.25,max_bits=760,tilt=0.05,explicit_bits=14,power_terms=8,odd=1):
  self.B,self.X,self.minimum,self.step=B_bits,X_bits,prime_min,step;self.tilt,self.length=tilt,math.ceil(max_bits/step);names=("small","sieve_geo","large_geo","sieve_sf","large_sf")
  weights={k:[0.0]*(self.length+1) for k in names};af=collections.Counter();bo=min(float(explicit_bits),B_bits);lm=math.log2(prime_min)
  self._skey=(prime_min,step,tilt,explicit_bits,power_terms,odd) if (bo==float(explicit_bits) and B_bits>=lm) else None
  sm_hit=self._skey is not None and self._skey in _RET_SM_W and len(_RET_SM_W[self._skey][0])>self.length
  if sm_hit:
   sw0,st0=_RET_SM_W[self._skey];weights["small"]=sw0[:self.length+1];af["small"]=st0
  else:
   for p in primes_upto(int(2**bo)-1)[odd:]:
    b=math.log2(p);bi=round(b/step);q=math.exp(tilt*bi*step-LN2*b)
    if p<prime_min:groups=(("small",False),)
    else:
     pref="sieve" if b<X_bits else "large";groups=((pref+"_geo",False),(pref+"_sf",True))
    for name,sf in groups:
     af[name]+=math.log1p(q) if sf else-math.log1p(-q);qk=q
     for k in range(1,self.length+1):
      idx=k*bi
      if idx>self.length:break
      sgn=-1 if sf and k%2==0 else 1;weights[name][idx]+=sgn*qk/k;qk*=q
      if qk<1e-25:break
  def cnt(left,right,pref):
   if left>=right:return
   w0=step/4.0;cells=math.ceil((right-left)/w0);w=(right-left)/cells;_exp=math.exp
   for i in range(cells):
    b=left+(i+0.5)*w;bi=round(b/step);tb=tilt*bi*step;lb=LN2*b
    for k in range(1,power_terms+1):
     amt=_exp(k*tb-(k-1)*lb)*w/(k*b);idx=k*bi
     if k>=2 and amt<1e-20:break
     sgn=-1 if k%2==0 else 1;gs=((pref,amt),) if pref=="small" else((pref+"_geo",amt),(pref+"_sf",sgn*amt))
     for g,x in gs:
      af[g]+=x
      if idx<=self.length:weights[g][idx]+=x
  if not sm_hit:
   cnt(bo,min(lm,B_bits),"small")
   if self._skey is not None:_RET_SM_W[self._skey]=(weights["small"],af["small"])
  cnt(max(bo,lm),min(X_bits,B_bits),"sieve");cnt(max(bo,X_bits),B_bits,"large")
  self.weights,self.normalizers=weights,dict(af);all_w=[sum(weights[k][j] for k in ("small","sieve_geo","large_geo")) for j in range(self.length+1)]
  self.denominator=_ret_coeffs(all_w,sum(af[k] for k in ("small","sieve_geo","large_geo")),self.length)
  self.sel_norm=math.exp(af["sieve_sf"]+af["large_sf"]-af["sieve_geo"]-af["large_geo"]);self._c={}
 def co(self,k,lim=None):
  L=self.length if lim is None else min(self.length,lim)
  if k=="small" and self._skey is not None:
   ck=(self._skey,L)
   if ck not in _RET_SM_C:_RET_SM_C[ck]=_ret_coeffs(self.weights[k],self.normalizers[k],L)
   return _RET_SM_C[ck]
  ck=(k,L)
  if ck not in self._c:self._c[ck]=_ret_coeffs(self.weights[k],self.normalizers[k],L)
  return self._c[ck]
 def two(self,n,a,b=9**9,win=None,sf=False):
  a=min(a,n)
  if a<0 or b<0:return 0.0
  D=min(self.length,1080 if a<=1080 else 2048);mkey=(D,b,win,sf)
  if mkey not in self._c:
   _sp=math.sumprod;st=self.step;sm=self.co("small",D);lim_b=min(b,D);lg=self.co("large_sf" if sf else "large_geo",D)
   if win is None:lg_f=lg[:lim_b+1]
   else:
    w0,w1,w2=win;lg_f=[v if (v and (l*st<=w2 or any(self.X*k-w2<=l*st<=w1*k+w2 for k in range(1,w0+1)))) else 0.0 for l,v in enumerate(lg[:lim_b+1])]
   while len(lg_f)>1 and lg_f[-1]==0.0:
    lg_f.pop()
   chunks=[];i=0;Llg=len(lg_f)
   while i<Llg:
    while i<Llg and lg_f[i]==0.0:i+=1
    if i>=Llg:break
    i0=i
    while i<Llg and lg_f[i]!=0.0:i+=1
    chunks.append((i0,i,lg_f[i0:i][::-1]))
   m=[0.0]*(D+1)
   for s in range(D+1):
    acc=0.0
    for i0,i1,lgr in chunks:
     if i0>s:break
     hi=s if s<i1-1 else i1-1
     if hi==i0:acc+=lgr[-1]*sm[s-i0]
     else:acc+=_sp(lgr[i1-1-hi:i1-i0],sm[s-hi:s-i0+1])
    m[s]=acc
   self._c[mkey]=m
  m=self._c[mkey];sv=self.co("sieve_sf" if sf else "sieve_geo");L=min(n,a,len(m)-1)+1
  return math.sumprod(m[:L],sv[n:n-L:-1] if n>=L else sv[n::-1])
 def cap(self,nb):
  return math.floor(math.ceil(nb)/math.log2(self.minimum))
 def at(self,nb,sl,sb=None):
  n=round(nb/self.step);rb=(self.cap(nb)+1)/16.0;r=self.sel_norm*self.two(n,math.floor((sl-rb)/self.step),9**9 if sb is None else math.floor((sb-rb-0.125)/self.step),None,True,)/self.denominator[n]
  return max(0.0,min(1.0,r))
 def survivor_upper(self,nb,sl,ex=0.0,sb=None,lb=None,kmax=4):
  n=round(nb/self.step);a=math.floor((sl+ex)/self.step)
  if sb is None:return self.two(n,a)/self.denominator[n]
  return self.two(n,a,math.floor((sb+0.25)/self.step),(kmax,lb,(self.cap(nb)+2)/16.0+0.25))/self.denominator[n]
@functools.lru_cache(None)
def _semismooth_ecm(fbits,b1,b2,steps=256):
 m=fbits-math.log2(ECM_TORSION);base,top=math.log2(b1),math.log2(b2);RINV,cz,ra=_get_rho_grid()
 def _rf(u):
  if u<=1.0:return 1.0
  p=u*RINV;i=int(p)
  if i+1<cz:return ra[i]+(p-i)*(ra[i+1]-ra[i])
  return rho_exact(u)
 acc=_rf(m/base)*(1.0-1e-9)
 if top>base:
  d=(top-base)/steps;ib=1.0/base;l0=math.log(base);c9=1.0-1e-9
  for i in range(steps):
   lo=base+d*i;l1=math.log(lo+d);acc+=_rf((m-lo)*ib)*c9*(l1-l0);l0=l1
 return min(acc,1.0-1e-15)
@functools.lru_cache(None)
def _ecm_cached_moments(fbits2,limbs,success=0.999):
 fbits=fbits2/2.0;nbits=limbs*64;mul=MONT5*(limbs/5.0)**2;gcd_c=2.0*nbits*(12*limbs+16);best=None;k_opt=2*round(4.5+fbits/19.0)
 for k in sorted(range(10,39,2),key=lambda x:abs(x-k_opt)):
  b1=1<<k
  for ratio in (300,100,30,10,1):
   b2=b1*ratio;m=fbits-math.log2(ECM_TORSION);u=m/math.log2(b1);r_lo=max(rho_exact(u)*(1.0-1e-9),1e-300);tilt=(rho_exact(max(u-2.0,1.0))*(1.0+1e-9))/r_lo;y=math.sqrt(b1)
   pcorr=max(0.0,1.0-tilt*2.0/(y*math.log(y)))
   if pcorr<=0.0:continue
   b1l,b2l=math.log2(b1),math.log2(b2)
   p_ub=min(0.999999,rho_exact(max(0.0,(m-b2l)/b1l))*(1.0+math.log(b2l/b1l))*pcorr)
   if p_ub<=1e-12:continue
   if best is not None and (math.log1p(-success)/math.log1p(-p_ub))*(1.1*1.5*b1*11*mul+2*gcd_c)>=best[0]:continue
   p1=_semismooth_ecm(fbits,b1,b2)*pcorr
   if p1<=1e-12:continue
   bp=math.ceil(math.log1p(-success)/math.log1p(-p1));cp=(1.0-p1)**bp;pi_up=1.26*b1/math.log(b1) if b1>=17 else 7.0
   for wheel in (ECM_WHEELS if ratio>1 else(0,)):
    if wheel and wheel>b2:continue
    if wheel:
     phi=wheel
     for q in (2,3,5,7,11,13):
      if wheel%q==0:phi=phi//q*(q-1)
     G=math.ceil(b2/wheel)+1;ops=1.1*(math.ceil(1.5*b1)*11+6*(wheel//2+G)+3*(phi//2)*G)
    else:
     G=0;ops=1.1*math.ceil(1.5*b1)*11
    per_c=ops*mul+2*gcd_c;b5=ops*mul+(pi_up+G)*gcd_c;redo=p1*b5;cost=bp*(per_c+redo);bu=bp*(per_c+b5);et=cost/(1.0-cp)
    et2=(bu*bu+2.0*cp*(bp*per_c)*et)/(1.0-cp)
    if best is None or cost<best[0]:
     best=(cost,et,et2,cp,bp,per_c,b5)
 return best
def ecm_instr(fbits,nbits,success=0.999):
 return _ecm_cached_moments(math.ceil(2.0*fbits),max(5,math.ceil(nbits/64)),success)[0]
def ecm_instr_moments(fbits,nbits,success=0.999):
 return _ecm_cached_moments(math.ceil(2.0*fbits),max(5,math.ceil(nbits/64)),success)
NTT4_PRIMES=(269*(1<<53)+1,276*(1<<53)+1,305*(1<<53)+1,309*(1<<53)+1,);NTT4_PROD=NTT4_PRIMES[0]*NTT4_PRIMES[1]*NTT4_PRIMES[2]*NTT4_PRIMES[3]
def fft_F(N):
 if N<=1:return 1200.0
 Ni=math.ceil(N)
 if 2*Ni<=(1<<53):
  k=max(2,(2*Ni-1).bit_length());M=1<<k;assert M<=(1<<53) and Ni*(1<<128)<NTT4_PROD
  return 168.0*M*k+220.0*M+1200.0+(504.0*M if Ni>(1<<26) else 0.0)
 ktot=(Ni-1).bit_length();k1=ktot//2;K1=1<<k1;K2=(Ni+K1-1)//K1;k1m=k1+1;k2m=max(2,(2*K2-1).bit_length());M1=1<<k1m;M2=1<<k2m;assert M1<=(1<<53) and M2<=(1<<53) and (K1*K2)*(1<<128)<NTT4_PROD
 Mtot=float(M1*M2)
 return 168.0*Mtot*(k1m+k2m)+1100.0*Mtot+1200.0
def mr_instr(nbits):
 L=max(5,math.ceil(nbits/64.0))
 return 100.0*2.0*nbits*MONT5*((L/5.0)**2)
def gbatch_cost(pbits,leaf=9):
 SP=math.ceil(pbits/64);L=math.ceil(math.log2((SP+1)/leaf));G=L+2;n=lambda j:leaf*(1<<j);PT=sum((1<<(L-j))*(fft_F(2*n(j-1))+3.0*2*n(j-1)) for j in range(1,L+1))
 RT=sum((1<<(L-j+1))*(fft_F(n(j)+G+n(j-1))+3.0*(n(j-1)+G)) for j in range(1,L+1))
 return (PT+RT+6.0*fft_F(2*(max(SP,n(L))+G)))/(1<<L)+(1<<15)
def _get_rho_grid():
 return _RHO_S,(_RHO_UMAX-1)*_RHO_S,_RHO_TAIL
def _dw_int(B,dn,dd):
 n=len(B);_sp=math.sumprod;d=[1]*(n+1);lam=[[0]*n for _ in range(n)]
 for i in range(n):
  bi=B[i]
  for j in range(i+1):
   u=_sp(bi,B[j])
   for k in range(j):u=(d[k+1]*u-lam[i][k]*lam[j][k])//d[k]
   if j<i:lam[i][j]=u
   else:d[i+1]=u
 k=1
 while k<n:
  lk=lam[k];dk=d[k];lkm1=lk[k-1]
  if 2*abs(lkm1)>dk:
   dk2=dk<<1;q=(2*lkm1+dk)//dk2 if lkm1>=0 else-((-2*lkm1+dk)//dk2);B[k]=[x-q*y for x,y in zip(B[k],B[k-1])];lk1=lam[k-1]
   if k>1:lk[:k-1]=[x-q*y for x,y in zip(lk[:k-1],lk1[:k-1])]
   lkm1-=q*dk;lk[k-1]=lkm1
  if dd*(d[k-1]*d[k+1]+lkm1*lkm1)<dn*dk*dk:
   B[k],B[k-1]=B[k-1],B[k];lk1=lam[k-1]
   if k>1:lk[:k-1],lk1[:k-1]=lk1[:k-1],lk[:k-1]
   dk_new=(d[k-1]*d[k+1]+lkm1*lkm1)//dk;dk1=d[k+1]
   for i in range(k+1,n):
    li=lam[i];t=li[k];v1=(dk1*li[k-1]-lkm1*t)//dk;li[k]=v1;li[k-1]=(dk_new*t+lkm1*v1)//dk1
   d[k]=dk_new;k=max(k-1,1)
  else:
   for l in range(k-2,-1,-1):
    dl1=d[l+1];lkl=lk[l]
    if 2*abs(lkl)>dl1:
     dl2=dl1<<1;q=(2*lkl+dl1)//dl2 if lkl>=0 else-((-2*lkl+dl1)//dl2);B[k]=[x-q*y for x,y in zip(B[k],B[l])];ll=lam[l]
     if l>0:lk[:l]=[x-q*y for x,y in zip(lk[:l],ll[:l])]
     lk[l]-=q*dl1
   k+=1
 return B,d,lam
def _int_fp_pre(B,F=256):
 n=len(B);_sp=math.sumprod;F2=2*F;S=1<<F;H=S>>1;S2=1<<F2
 Bs=[[0]*n for _ in range(n)];mu=[[0]*n for _ in range(n)];Bn=[0]*n
 for i in range(n):
  bi=B[i];vv=[x<<F for x in bi];mui=mu[i]
  for j in range(i):
   bnj=Bn[j]
   if bnj==0:mui[j]=0;continue
   num=_sp(bi,Bs[j])<<F2;m=(num+(bnj>>1))//bnj if num>=0 else -((-num+(bnj>>1))//bnj);mui[j]=m;bsj=Bs[j]
   if m:vv=[vv[t]-((m*bsj[t])>>F) for t in range(n)]
  Bs[i]=vv;Bn[i]=_sp(vv,vv)
 k=1;sw=0
 while k<n and sw<6000:
  muk=mu[k];m=muk[k-1]
  if abs(m)>H:
   q=(m+H)>>F if m>=0 else -((-m+H)>>F)
   if q:
    B[k]=[x-q*y for x,y in zip(B[k],B[k-1])];m1=mu[k-1]
    for l in range(k-1):muk[l]-=q*m1[l]
    m-=q*S;muk[k-1]=m
  bnk,bnk1=Bn[k],Bn[k-1]
  if bnk1==0 or 4*(bnk<<F2)>=(3*S2-4*m*m)*bnk1:
   for j in range(k-2,-1,-1):
    mj=muk[j]
    if abs(mj)>H:
     q=(mj+H)>>F if mj>=0 else -((-mj+H)>>F)
     if q:
      B[k]=[x-q*y for x,y in zip(B[k],B[j])];mj0=mu[j]
      for l in range(j):muk[l]-=q*mj0[l]
      muk[j]-=q*S
   k+=1
  else:
   sw+=1;B[k],B[k-1]=B[k-1],B[k];bn_n=bnk+((m*m*bnk1)>>F2)
   if bn_n==0:break
   mn=(m*bnk1)//bn_n;Bn[k]=(bnk*bnk1)//bn_n;Bn[k-1]=bn_n;muk[k-1]=mn;m1=mu[k-1]
   for j in range(k-1):muk[j],m1[j]=m1[j],muk[j]
   for i in range(k+1,n):
    mi=mu[i];tv=mi[k];v1=mi[k-1]-((m*tv)>>F);mi[k]=v1;mi[k-1]=tv+((mn*v1)>>F)
   k=max(k-1,1)
 return B,sw
def lll_16d_int(B_in,delta_num=99,delta_den=100):
 B=[list(r) for r in B_in]
 if delta_num*40>delta_den*11 and max(abs(x).bit_length() for r in B for x in r)>160:
  for _ in range(14):
   B,sw=_int_fp_pre(B,256)
   if sw==0:break
 return _dw_int(B,delta_num,delta_den)
def run_boot_381bit_lll_live(gam,h_pi,p_bls,s_boot,pmul_fn,prem_fn):
 B_sk=[]
 for i in range(4):
  r=[0]*16;r[i]=p_bls;B_sk.append(r)
 for j in range(4):
  r=[0]*16
  for idx,ch in enumerate(h_pi):r[j+idx]=ch
  B_sk.append(r)
 for k in range(8):
  yk=[0]*8;yk[k]=1;ak=pmul_fn(gam,prem_fn(yk,h_pi,p_bls),p_bls,h_pi);r=[0]*16
  for idx,ch in enumerate(ak):r[idx]=ch
  r[8+k]=s_boot;B_sk.append(r)
 R_sk,d_gso,lam=lll_16d_int(B_sk,99,100);assert d_gso[16]==(p_bls**8)*(s_boot**16)
 for kk in range(1,16):
  assert 100*(d_gso[kk-1]*d_gso[kk+1]+lam[kk][kk-1]**2)>=99*d_gso[kk]**2
  for jj in range(kk):
   assert 2*abs(lam[kk][jj])<=d_gso[jj+1]
 assert all(r[8+kk]%s_boot==0 for r in R_sk for kk in range(8))
 return R_sk
def srt_toy_check():
 import random
 rng=random.Random(2026);bf=[];n=3
 while len(bf)<400:
  if all(n%d!=0 for d in range(2,math.isqrt(n)+1)):bf.append(n)
  n+=2
 P=1
 for pr in bf:P*=pr
 da=[rng.randrange(1<<40,1<<48)|1 for _ in range(128)];am=[da]
 while len(am[-1])>1:
  cur=am[-1];am.append([cur[i]*cur[i+1] for i in range(0,len(cur),2)])
 L=len(am)-1;G=L+2;Z_root=am[-1][0];b_root=64*((Z_root.bit_length()+63)//64+G);Y=[(P<<b_root)//Z_root]
 for lv in range(L-1,-1,-1):
  nxt=[]
  for idx,zc in enumerate(am[lv]):
   zs=am[lv][idx^1];yv=Y[idx>>1];bv=64*((am[lv+1][idx>>1].bit_length()+63)//64+G);bc=64*((zc.bit_length()+63)//64+G);yc=((yv*zs)&((1<<bv)-1))>>(bv-bc);nxt.append(yc)
  Y=nxt
 mis=0
 for z,y in zip(da,Y):
  bc=64*((z.bit_length()+63)//64+G);r=(z*y+(1<<bc)-1)>>bc
  if r==z:r=0
  if r!=P%z:mis+=1
 assert mis==0
 return dict(n_leaves=len(da),mismatches=mis)
def stage_b_analytic(nbits):
 L=max(5,math.ceil(nbits/64.0));c_srt=gbatch_cost(24197143.0,leaf=L);c_sq10=10.0*MONT5*((L/5.0)**2);c_gcd=2.0*nbits*(12.0*L+16.0)
 return c_srt+c_sq10+c_gcd
def _peel_fixed_point_bls(lf,lg,lrels,u_f,u_g,exc_ref=1.09,step=0.1,lam_seed=100.0):
 rf=rho_exact(u_f);rg=rho_exact(u_g);bins=[]
 for cm,u_s,r_s in ((lf,u_f,rf),(lg,u_g,rg)):
  b=40.0
  while b<cm-1e-9:
   c=min(cm,b+step);km_lo=rho_exact(u_s-b/cm)/r_s;km_hi=rho_exact(u_s-c/cm)/r_s;wt=km_hi*math.log(c/b);lam=(2.0**(lrels-c))*12.0*km_lo
   if lam<lam_seed:bins.append((wt,lam))
   b=c
 _exp=math.exp
 def aj(scale,stop90=False):
  h=0.0;sbins=[(wt,lam*scale) for wt,lam in bins]
  for _ in range(400):
   omh=1.0-h;s=ds=0.0
   for wt,ls in sbins:
    t=wt*_exp(-ls*omh);s+=t;ds+=t*ls
   em=_exp(-s);g=1.0-em
   if stop90 and g>0.90:return g
   diff=g-h
   if abs(diff)<1e-13:return g
   gp=em*ds
   if 0.0<=gp<0.98 and diff>0.0:
    hn=h+0.85*diff/(1.0-gp)
    if g<hn<1.0:
     omhn=1.0-hn;sn=sum(wt*_exp(-ls*omhn) for wt,ls in sbins);gn=1.0-_exp(-sn)
     if stop90 and gn>0.90:return gn
     if gn>=hn:
      h=gn
      continue
   h=g
  return h
 def exceeds_90(scale):
  return aj(scale,True)>0.90
 h=aj(1.0);lo_e,hi_e=0.4,exc_ref
 for _ in range(16):
  mid=0.5*(lo_e+hi_e)
  if exceeds_90(mid/exc_ref):lo_e=mid
  else:hi_e=mid
 return h,hi_e
def eval_descent_live(desc,mesh_node=96,mesh_boot=192):
 dp=desc["params"];lx=float(dp["lx"]);di=float(dp.get("peel_lf",71.29));cq=float(dp.get("peel_lg",71.69));bl=float(dp.get("peel_lrels",63.4417118678))
 dc=float(dp.get("peel_uf",548.9/71.29));cb=float(dp.get("peel_ug",677.8/71.69));exc_peel=float(dp.get("peel_exc",1.09))
 peel_h,exc_crit=_peel_fixed_point_bls(di,cq,bl,dc,cb,exc_ref=exc_peel);alpha_f,alpha_g=float(dp["alpha_f"]),float(dp["alpha_g"]);ba=float(dp["log2_u8"])
 smult=float(dp["sieve_mult"]);qbins=desc["qbins"];qbins_lo=desc["qbins_lo"];bx=desc["qbins_sub"];cells_r2={int(k):int(v) for k,v in desc["cells_rank2"].items()};BOUND_CELL=2400.0*smult
 HITS_PER_CELL=2.0*math.log(lx/8.0);SIEVE_HIT=1100.0*smult;SCAN=10.0*smult;ZERO=2.0*smult;PAIR_SETUP=(2.0**25.1)*(smult*1.33);d3=float(dp["lll_node"])
 NORMS,TDIV,IDEAL,PRIME=2.0**27,2.0**25.0,2.0**35,2.0**24;pairs=2.0*(2.0**log2_li(lx));cfmax=max(x[1] for x in qbins);cgmax=max(x[2] for x in qbins);wrc={};RINV,cz,r_arr=_get_rho_grid()
 RSTEP=1.0/RINV
 def get_wr(n,iA,mesh,kA,kB):
  key=(n,iA,mesh,kA,kB)
  if key in wrc:return wrc[key]
  h=(n-lx)/mesh;wA=[math.log1p(h/(lx+i*h)) if i<iA else 0.0 for i in range(mesh)];wB=[math.log1p(h/(lx+i*h)) if i>=iA else 0.0 for i in range(mesh)]
  def conv(a,b):
   _sp=math.sumprod;la,lb=len(a),len(b);out=[0.0]*(la+lb-1);nz=[j for j,y in enumerate(b) if y]
   if not nz:return out
   j0,j1=nz[0],nz[-1]+1;br=b[j0:j1][::-1];lbr=len(br)
   for s in range(la+lbr-1):
    lo=max(0,s-lbr+1);hi=min(la,s+1);out[s+j0]=_sp(a[lo:hi],br[lbr-1-(s-lo):lbr-(s-hi+1)])
   return out
  pAl=[[1.0]]
  for j in range(1,kA+1):
   pAl.append([x/j for x in conv(pAl[-1],wA)])
  pABl=[conv(pA,wB) for pA in pAl] if kB else None
  rows_meta=[]
  for jA in range(kA+1):
   for jB,row in (((0,pAl[jA]),(1,pABl[jA])) if kB else((0,pAl[jA]),)):
    nz=[i for i,c in enumerate(row) if c]
    if nz:rows_meta.append((jA+jB,nz[0],nz[-1],row))
  wrc[key]=(h,rows_meta)
  return wrc[key]
 def cd(m_list,n,iA,mesh,alpha,kA,kB):
  if n<=lx:
   out=[]
   for m_raw in m_list:
    m=float(math.ceil(m_raw));u_val=(m+alpha)/lx;rv=rho_exact(u_val)
    if rv>0.0 and u_val>2.0:rv*=max(0.0,1.0-GAMMA*rho_exact(u_val-1.0)/(rv*m*LN2))
    out.append(rv)
   return out
  h,rows_meta=get_wr(n,iA,mesh,kA,kB);hk=h*RINV/lx;ilr=RINV/lx;g_ln2=GAMMA/LN2;two_r=2*RINV;lim16=16*RINV;dh={};out=[];ra=r_arr
  for m_raw in m_list:
   mi=math.ceil(m_raw);cached=dh.get(mi)
   if cached is not None:
    out.append(cached)
    continue
   m=float(mi);tot=0.0
   for j,nz_lo,nz_hi,row in rows_meta:
    ifeas=min(nz_hi,int(math.floor((m-j*h-j*lx)/h))) if h>0 else nz_hi
    if ifeas<nz_lo:continue
    k0=(m-j*lx+alpha)*ilr;aw=max(nz_lo,int((k0-lim16)/hk));mj_lx=m-j*lx;idx_z=int((k0+0.999999999999)/hk)+1
    if idx_z<=ifeas:
     if idx_z<=aw:tot+=sum(row[aw:ifeas+1]);continue
     tot+=sum(row[idx_z:ifeas+1]);ifeas=idx_z-1
    cr=max(aw,int((k0-two_r+0.999999999999)/hk)+1)
    if cr<=ifeas:
     k_cur=k0-cr*hk+0.999999999999
     for idx in range(cr,ifeas+1):
      ks=int(k_cur);k_cur-=hk
      if ks<=0:tot+=row[idx]
      elif ks<cz:tot+=row[idx]*ra[ks]
     ifeas=cr-1
    if aw<=ifeas:
     k_cur=k0-aw*hk+0.999999999999;rem_b=mj_lx-aw*h
     for idx in range(aw,ifeas+1):
      ks=int(k_cur);k_cur-=hk
      if ks<cz:
       c=row[idx]
       if c:
        rv=ra[ks]
        if rv>0.0:
         rb=rem_b if rem_b>=lx else lx;v=rv-g_ln2*ra[ks-RINV]/rb
         if v>0.0:tot+=c*v
      rem_b-=h
   dh[mi]=tot;out.append(tot)
  return out
 p_unrec=float(dp.get("p_unrec_leaf",0.0));mc_={}
 def b8(m_val,u_m,lb_s):
  if m_val<=0.0:return 0.0
  u4=math.ceil(4.0*m_val/lx);key=(u4,u_m,lb_s)
  if key in mc_:return mc_[key]
  u=u4/4.0;r_u=rho_exact(u);r_m=rho_exact(u_m);tot=0.0
  for b10 in range(400,int(lx*10)):
   b=b10/10.0;c=b+0.1;km_lo=rho_exact(u_m-b/lb_s)/r_m if r_m>0 else 0.0;lam=(2.0**(bl-c))*12.0*km_lo;ku_hi=rho_exact(u-c/lx)/r_u if r_u>0 else 0.0
   if lam<100.0:
    x=lam*(1.0-peel_h);p_miss=math.exp(-x)*(1.0+x)
   else:p_miss=math.exp(-lam)
   p_miss_tot=min(1.0,p_miss+p_unrec);tot+=ku_hi*math.log(c/b)*p_miss_tot
  mc_[key]=tot
  return tot
 _TG=tuple(ti*0.01 for ti in range(1,81))
 def _tphi(mu_v,psi_v,r_v):
  hm=0.5*mu_v/r_v;bst=min(1.0,(1.0+hm)**(-r_v)+psi_v)
  if bst<=0.01:return 0.01
  gs=[(1.0+hm*(1.0-t))**(-r_v) for t in _TG];tp=[1.0]*80;pk=psi_v
  for k in range(2,65):
   pk*=psi_v
   if pk>=bst:break
   lim=bst-pk
   for i in range(80):
    tp[i]*=_TG[i];gv=gs[i]/tp[i]
    if gv<lim:bst=gv+pk;lim=gv
   if bst<=0.01:return 0.01
  return max(0.01,bst)
 W_det,be,EN_dag={},{},{};bsn=[];orw=[];siA=mesh_node//32;p_ecm=float(dp["ecm_success"])
 for row in desc["plan_struct"]:
  lo,hi,side,n_int,iA32,kA,kB,R_n,thin=row;bn=(lo,hi)
  if not bsn or bsn[-1]!=bn:bsn.append(bn)
  n_c=float(n_int);iA=iA32*siA;cells=cells_r2[R_n];S=math.log2(cells);c_hi=(hi+S)/16.0;c_lo=(lo+S)/16.0;nA=lx+iA*(n_c-lx)/mesh_node if n_c>lx else lx;T=kA*nA+kB*n_c+96.0
  nf_list=[72.0*c_hi+cf-(lo if side=="f" else 0.0) for _,cf,_ in qbins]+[72.0*c_hi+cf-(lo if side=="f" else 0.0) for cf,_ in bx]
  ng_list=[24.0*c_hi+ba+cg-(lo if side=="g" else 0.0) for _,_,cg in qbins]+[24.0*c_hi+ba+cg-(lo if side=="g" else 0.0) for _,cg in bx]
  pf_v=cd(nf_list,n_c,iA,mesh_node,alpha_f,kA,kB);pg_v=cd(ng_list,n_c,iA,mesh_node,alpha_g,kA,kB);pr_c=sum(pf_v[i]*pg_v[i] for i in range(6))/6.0
  pr_i=sum(pf_v[i]*pg_v[i] for i in range(6,30))/24.0;pr=min(pr_i,4.0*pr_c);sv=0.0
  for (w,_,_),(cfl,cgl) in zip(qbins,qbins_lo):
   nfl=72.0*c_lo+cfl-(hi if side=="f" else 0.0);ngl=24.0*c_lo+ba+cgl-(hi if side=="g" else 0.0);sf=1.0 if T>=nfl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(nfl)-T)/lx))
   sg=1.0 if T>=ngl else min(1.0,max(1.0,T/lx)*rho_exact((math.floor(ngl)-T)/lx));sv+=w*sf*sg
  bh=72.0*c_hi+cfmax-(lo if side=="f" else 0.0);bg=24.0*c_hi+ba+cgmax-(lo if side=="g" else 0.0);k_lx_f=int(bh//20.0)+1;k_lx_g=int(bg//20.0)+1;k_wA_n=kA if nA>lx+1e-9 else 0
  det_n=p_ecm**(k_lx_f+k_lx_g+2*(kA+kB));miss_n=b8(bh,dc,di)+b8(bg,cb,cq);mu=cells*pr*float(dp["udup_node"])*det_n*math.exp(-miss_n);mr_f=mr_instr(bh)
  mr_g=mr_instr(bg);mr_T=mr_instr(T)
  svc=(NORMS+TDIV+k_lx_f*(ecm_instr(lx,bh,p_ecm)+mr_f)+k_lx_g*(ecm_instr(lx,bg,p_ecm)+mr_g)+2.0*k_wA_n*(ecm_instr(nA,T,p_ecm)+mr_T)+mr_f+mr_g+IDEAL);svs=cells*sv+cells*pr
  sieve=cells*(HITS_PER_CELL*SIEVE_HIT+SCAN+ZERO+BOUND_CELL)+pairs*PAIR_SETUP
  if cells*svc<sieve+svs*svc:
   how="ecm";own=d3+cells*svc
  else:
   how="sieve";mu*=float(dp["ret_node"]);own=d3+sieve+svs*svc
  wmax=lambda s,lim:max([W_det[b,s] for b in bsn if b[0]<lim-1e-9 and (b,s) in W_det]+[0.0]);pmax=lambda s,lim:max([be[b,s] for b in bsn if b[0]<lim-1e-9 and (b,s) in be]+[0.0])
  nmax=lambda s,lim:max([EN_dag[b,s] for b in bsn if b[0]<lim-1e-9 and (b,s) in EN_dag]+[0.0])
  psi=kA*(pmax("f",nA)+pmax("g",nA))+kB*(pmax("f",n_c)+pmax("g",n_c));assert psi<0.80
  phi=_tphi(mu,psi,0.5 if (kA+kB<=1) else 1.0);assert 0.01<=phi<=0.04;retry=1.0/(1.0-psi)
  W_det[bn,side]=own+retry*(kA*(wmax("f",nA)+wmax("g",nA))+kB*(wmax("f",n_c)+wmax("g",n_c)))
  be[bn,side]=phi;EN_dag[bn,side]=1.0+retry*(kA*(nmax("f",nA)+nmax("g",nA))+kB*(nmax("f",n_c)+nmax("g",n_c)))
  orw.append(dict(lo=lo,hi=hi,side=side,how=how,R=R_n,n=n_int,nA=nA,kA=kA,kB=kB,mu=mu,psi=psi,phi=phi,nf=bh,ng=bg,own_log2=math.log2(own),W_log2=math.log2(W_det[bn,side]),EN=EN_dag[bn,side],))
 bt=desc["boot"];assert int(bt["cells_b"])==1 and int(dp["boot_rels"])==1 and float(dp.get("ball_pair",0.0))==0.0;m_max=float(bt["m_sv"]);iA_b=round((bt["nA"]-lx)*mesh_boot/(bt["n1"]-lx))
 pr_raw_b=cd([m_max],float(bt["n1"]),iA_b,mesh_boot,alpha_g,bt["kA"],1)[0];mb=b8(m_max,cb,cq);pr_b=pr_raw_b*math.exp(-mb);k_lx_b=int(m_max//20.0)+1;kA_b=int(bt["kA"])
 mr_b=mr_instr(m_max);det_b=p_ecm**(k_lx_b+kA_b+1);p_boot_trial=pr_b*det_b*float(dp["ret_boot"])
 trial_b=(NORMS+TDIV+k_lx_b*(ecm_instr(lx,m_max,p_ecm)+mr_b)+kA_b*(ecm_instr(bt["nA"],m_max,p_ecm)+mr_b)+mr_b+IDEAL);LLL_BOOT=2.0**48;boe=(1.0/p_boot_trial)*(LLL_BOOT+trial_b)
 wmax_g=lambda lim:max(W_det[b,"g"] for b in bsn if b[0]<lim-1e-9);pmax_g=lambda lim:max(be[b,"g"] for b in bsn if b[0]<lim-1e-9)
 nmax_g=lambda lim:max(EN_dag[b,"g"] for b in bsn if b[0]<lim-1e-9);do=bt["kA"]*wmax_g(bt["nA"])+wmax_g(float(bt["n1"]));Phi_1=bt["kA"]*pmax_g(bt["nA"])+pmax_g(float(bt["n1"]))
 EN_1=bt["kA"]*nmax_g(bt["nA"])+nmax_g(float(bt["n1"]));assert Phi_1<0.50;wald_mult=1.0/(1.0-Phi_1)
 two_exp_u=2.0*wald_mult*(boe+do)/20000.0;EN_2log=2.0*wald_mult*EN_1
 return dict(mesh_node=mesh_node,mesh_boot=mesh_boot,n_nodes=len(desc["plan_struct"]),peel_h=peel_h,exc_crit=exc_crit,boot_pr_log2=math.log2(pr_b),boot_miss=mb,p_boot_trial=p_boot_trial,boot_own_exp_log2=math.log2(2.0*boe/20000.0),dag_sub_two_log2=math.log2(2.0*do/20000.0),phi_boot=Phi_1,wald_mult=wald_mult,two_log_exp_units=math.log2(two_exp_u),log2_EN_2log=math.log2(EN_2log),min_mu=min(r["mu"] for r in orw),max_phi=max(r["phi"] for r in orw),max_psi=max(r["psi"] for r in orw),rows=orw,)
_KMB={}
def _kmb(S):
 e=_KMB.get(S)
 if e is None:
  M=1<<S;H=M>>1;mk=M-1;e=(H,mk,(((1<<(15*S))-1)//mk)*H,2*S,3*S,4*S,5*S,6*S,7*S,8*S,9*S,10*S,11*S,12*S,13*S,14*S);_KMB[S]=e
 return e
def _kmul15(a,b,S=0):
 if not S:S=max(max(a),-min(a),max(b),-min(b)).bit_length()*2+10
 H,mk,bi,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14=_kmb(S)
 P=a[0]+(a[1]<<S)+(a[2]<<s2)+(a[3]<<s3)+(a[4]<<s4)+(a[5]<<s5)+(a[6]<<s6)+(a[7]<<s7)
 Q=(P*P if a is b else P*(b[0]+(b[1]<<S)+(b[2]<<s2)+(b[3]<<s3)+(b[4]<<s4)+(b[5]<<s5)+(b[6]<<s6)+(b[7]<<s7)))+bi
 c0=(Q&mk)-H;c1=((Q>>S)&mk)-H;c2=((Q>>s2)&mk)-H;c3=((Q>>s3)&mk)-H;c4=((Q>>s4)&mk)-H;c5=((Q>>s5)&mk)-H;c6=((Q>>s6)&mk)-H;c7=((Q>>s7)&mk)-H
 c8=((Q>>s8)&mk)-H;c9=((Q>>s9)&mk)-H;c10=((Q>>s10)&mk)-H;c11=((Q>>s11)&mk)-H;c12=((Q>>s12)&mk)-H;c13=((Q>>s13)&mk)-H;c14=(Q>>s14)-H
 return [c0-c8-c9-c10+c13+c14,c1+c8-c11-c13,c2+c9-c12-c14,c3-c8-c9+c14,c4+c8-c13-c14,c5-c8-c10+c13,c6-c9-c11+c14,c7+c8+c9-c12-c13-c14]
def _knorm_zomega(a,S1=0,S2=0):
 p1=_kmul15(a,(a[0]-a[2]-a[6]+a[7],a[2]+a[4]-a[7],-a[3]+a[6],-a[2]-a[6],a[1]+a[2]-a[7],-a[2]+a[5]+a[7],-a[6],a[2]-a[3]+a[6]-a[7]),S1)
 q=(p1[0]+p1[2]+p1[4],-p1[4],-p1[2]-p1[6],p1[2],-p1[2]-p1[4]+p1[7],p1[4]+p1[5],p1[2]+p1[3],p1[1]-p1[2]-p1[4]-p1[6])
 S=S2 or max(max(p1),-min(p1),max(q),-min(q)).bit_length()*2+10;H,mk,bi,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14=_kmb(S)
 Q=(p1[0]+(p1[1]<<S)+(p1[2]<<s2)+(p1[3]<<s3)+(p1[4]<<s4)+(p1[5]<<s5)+(p1[6]<<s6)+(p1[7]<<s7))*(q[0]+(q[1]<<S)+(q[2]<<s2)+(q[3]<<s3)+(q[4]<<s4)+(q[5]<<s5)+(q[6]<<s6)+(q[7]<<s7))+bi
 c0=(Q&mk)-H;c5=((Q>>s5)&mk)-H;c8=((Q>>s8)&mk)-H;c9=((Q>>s9)&mk)-H;c10=((Q>>s10)&mk)-H;c13=((Q>>s13)&mk)-H;c14=(Q>>s14)-H
 u=c0-c8-c9-c10+c13+c14;v=c5-c8-c10+c13;return u,v,u*u-u*v+v*v
def q15_v5l_descent_live_replay(camp):
 import random
 assert camp["same_field_q15"] and camp["exact_tau_orbit_keyed"] and camp["real_2core_fb"] and camp["literal_k1_first_usable"];a8=int(camp["b_fb_bits"]);B_FB=1<<a8
 coll_seed=int(camp["coll_seed"]);coll_pairs=int(camp["coll_pairs"]);d1_lo,d1_hi,nc_1=float(camp["d1_lo"]),float(camp["d1_hi"]),float(camp["nc_1"])
 d2_lo,d2_hi,nc_2=float(camp["d2_lo"]),float(camp["d2_hi"]),float(camp["nc_2"]);nc_0=float(camp["nc_0"]);af_t,ag_t=float(camp["alpha_f"]),float(camp["alpha_g"]);mul15=_kmul15
 _cb3={}
 def _gonly(a,b):
  a3=mul15(mul15(a,a,20),a,32);kb=tuple(b);b3=_cb3.get(kb)
  if b3 is None:
   b3=mul15(mul15(b,b,20),b,32)
   if len(_cb3)<8192:_cb3[kb]=b3
  return a3,b3,_knorm_zomega((a3[0]+2*(-b3[0]+b3[3]+b3[4]+b3[5]),a3[1]+2*(-b3[1]-b3[3]+b3[6]),a3[2]+2*(-b3[2]-b3[4]+b3[7]),a3[3]+2*b3[4],a3[4]+2*(-b3[3]-b3[4]),a3[5]+2*(-b3[0]+b3[3]),a3[6]+2*(-b3[1]+b3[4]),a3[7]+2*(-b3[2]-b3[3]-b3[4])),40,80)
 def _fonly(a3,b3):
  return _knorm_zomega((a3[0]-(2*b3[0]-b3[3]-b3[4]-b3[5]),a3[1]-(2*b3[1]+b3[3]-b3[6]),a3[2]-(2*b3[2]+b3[4]-b3[7]),a3[3]-(b3[3]-b3[4]),a3[4]-(b3[3]+2*b3[4]),a3[5]-(b3[0]-b3[3]+b3[5]),a3[6]-(b3[1]-b3[4]+b3[6]),a3[7]-(b3[2]+b3[3]+b3[4]+b3[7])),40,80)
 def _prod(xs):
  if not xs:return 1
  while len(xs)>1:xs=[xs[i]*xs[i+1] if i+1<len(xs) else xs[i] for i in range(0,len(xs),2)]
  return xs[0]
 ps=primes_upto(B_FB);P256=_prod([p for p in ps if p<=256]);PMID=_prod([p for p in ps if 256<p<=B_FB]);omr={}
 for p in ps:
  if p%3==1:
   g=2
   while pow(g,(p-1)//3,p)==1:g+=1
   r1=pow(g,(p-1)//3,p);r2=(p-1-r1)%p;omr[p]=(min(r1,r2),max(r1,r2))
 def d9(p):
  g=2
  while True:
   z=pow(g,(p-1)//15,p)
   if pow(z,3,p)!=1 and pow(z,5,p)!=1 and pow(z,15,p)==1:return z
   g+=1
 for p in ps:
  if p>200 or p%15!=1:continue
  z=d9(p);ys=[pow(z,k,p) for k in (1,2,4,7,8,11,13,14)]
  for sd in (0,1):
   ids=set()
   for y0 in ys:
    om=pow(y0,5,p);rhs=(-2*om*om)%p if sd==0 else(2+om)%p
    if pow(rhs,(p-1)//3,p)==1:
     r0=next(r for r in range(1,p) if pow(r,3,p)==rhs)
     for m in (0,1,2):ids.add((y0,(r0*pow(om,m,p))%p))
   rem=set(ids);oms=[]
   while rem:
    st=next(iter(rem));orb=set();cu=st
    for _ in range(12):orb.add(cu);cu=(pow(cu[0],13,p),(cu[1]*pow(cu[0],5,p))%p)
    assert len(orb)==12;s_om={pow(y0,5,p) for y0,_ in orb};assert len(s_om)==1;oms.append(next(iter(s_om)));rem-=orb
   assert len(oms)==len(set(oms))
 def cw(vals):
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
 Q_BIG=1152921504606847081;_z15b=d9(Q_BIG);_omb=pow(_z15b,5,Q_BIG);_ominv=(1,pow(_omb,-1,Q_BIG),pow(_omb,-2,Q_BIG))
 _ypb=tuple(tuple(pow(pow(_z15b,pow(7,s,15),Q_BIG),j,Q_BIG) for j in range(8)) for s in range(4))
 _yp0,_yp1,_yp2,_yp3=_ypb;_oi1,_oi2=_ominv[1],_ominv[2]
 def rtk(a,b):
  _sp=math.sumprod
  b0=_sp(b,_yp0)%Q_BIG;b1=_sp(b,_yp1)%Q_BIG;b2=_sp(b,_yp2)%Q_BIG;b3=_sp(b,_yp3)%Q_BIG
  p1=(b0*b1)%Q_BIG;p2=(p1*b2)%Q_BIG;p3=(p2*b3)%Q_BIG
  if p3==0:return None
  iv=pow(p3,-1,Q_BIG);i3=(iv*p2)%Q_BIG;iv=(iv*b3)%Q_BIG;i2=(iv*p1)%Q_BIG;iv=(iv*b2)%Q_BIG;i1=(iv*b0)%Q_BIG;i0=(iv*b1)%Q_BIG
  r0=(_sp(a,_yp0)%Q_BIG*i0)%Q_BIG;r1=(_sp(a,_yp1)%Q_BIG*i1)%Q_BIG;r2=(_sp(a,_yp2)%Q_BIG*i2)%Q_BIG;r3=(_sp(a,_yp3)%Q_BIG*i3)%Q_BIG
  return min(r0,(r0*_oi1)%Q_BIG,(r0*_oi2)%Q_BIG,r1,(r1*_oi1)%Q_BIG,(r1*_oi2)%Q_BIG,r2,(r2*_oi1)%Q_BIG,(r2*_oi2)%Q_BIG,r3,(r3*_oi1)%Q_BIG,(r3*_oi2)%Q_BIG)
 _ta=[1,-1,2,0,0,1,0,0];_tb=[0,1,0,-1,0,0,1,0];_k0=rtk(_ta,_tb);assert _k0==rtk([-x for x in _ta],[-x for x in _tb])==rtk([2*x for x in _ta],[2*x for x in _tb])
 _u1=[1,1,0,0,0,0,0,0];_om1=[0,0,0,0,0,1,0,0];assert _k0==rtk(mul15(_ta,_u1),mul15(_tb,_u1))==rtk(mul15(_ta,_om1),mul15(_tb,_om1));_sig_a=[0]*8;_sig_b=[0]*8
 for _j in range(8):
  _e=(7*_j)%15;_v=[1]+[0]*7
  for _ in range(_e):_v=[-_v[7],_v[0]+_v[7],_v[1],_v[2]-_v[7],_v[3]+_v[7],_v[4]-_v[7],_v[5],_v[6]+_v[7]]
  for _i in range(8):_sig_a[_i]+=_ta[_j]*_v[_i];_sig_b[_i]+=_tb[_j]*_v[_i]
 assert _k0==rtk(_sig_a,mul15(_om1,_sig_b));rng=random.Random(coll_seed);ngs,uvg,ab3=[],[],[];d5=set();_rb=rng._randbelow
 for _ in range(coll_pairs):
  a=[0]*8;b=[0]*8
  for _ in range(2+_rb(3)):a[_rb(8)]+=1 if _rb(2) else -1
  for _ in range(1+_rb(3)):b[_rb(8)]+=1 if _rb(2) else -1
  if not any(a) or not any(b):continue
  rk=rtk(a,b)
  if rk is None or rk in d5:continue
  d5.add(rk);a3,b3,(ug,vg,ng)=_gonly(a,b)
  if ng==0:continue
  ngs.append(ng);uvg.append((ug,vg));ab3.append((a3,b3))
 assert len(ngs)<=len(d5);rg=cw(ngs);keep=[i for i,x in enumerate(rg) if x==1];kf=[];kidx=[]
 for i in keep:
  uf,vf,nf=_fonly(*ab3[i])
  if nf:kf.append((nf,uf,vf));kidx.append(i)
 rf=cw([x[0] for x in kf])
 rels=[(kf[j][0],ngs[kidx[j]],(kf[j][1],kf[j][2]),uvg[kidx[j]]) for j,x in enumerate(rf) if x==1]
 def d1(n,u,v):
  out=set();rem=n;ps_h=[]
  for p in ps:
   if p*p>rem:break
   if rem%p==0:
    while rem%p==0:rem//=p
    ps_h.append(p)
  if rem>1:ps_h.append(rem)
  for p in ps_h:
   if p%3==1:
    r1,r2=omr[p]
    if (u+v*r1)%p==0:out.add((p,r1))
    if (u+v*r2)%p==0:out.add((p,r2))
   else:out.add((p,1 if p==3 else 0))
  return out
 ror=[(d1(nf,uf,vf),d1(ng,ug,vg)) for nf,ng,(uf,vf),(ug,vg) in rels];cur=list(range(len(ror)))
 while True:
  deg={}
  for ri in cur:
   ff,fg=ror[ri]
   for o in ff:deg[(1,o)]=deg.get((1,o),0)+1
   for o in fg:deg[(0,o)]=deg.get((0,o),0)+1
  nxt=[ri for ri in cur if all(deg[(1,o)]>=2 for o in ror[ri][0]) and all(deg[(0,o)]>=2 for o in ror[ri][1])]
  if len(nxt)==len(cur):break
  cur=nxt
 fb_g={o for ri in cur for o in ror[ri][1]};fb_f={o for ri in cur for o in ror[ri][0]}
 def bz(fb_s):
  ao,pl=[],[]
  for p in ps:
   if p%3==1:
    r1,r2=omr[p];i1,i2=(p,r1) in fb_s,(p,r2) in fb_s
    if not i1 and not i2:ao.append(p)
    elif not i1:pl.append((p,r1))
    elif not i2:pl.append((p,r2))
   elif p==3:
    if (3,1) not in fb_s:ao.append(3)
   else:
    if (p,0) not in fb_s:ao.append(p)
  return _prod(ao),_prod([p for p,_ in pl]),pl
 flt_g=bz(fb_g);flt_f=bz(fb_f);assert [len(rels),len(cur)]==camp["raw_rels_2core"] and [len(fb_g),len(fb_f)]==camp["fb_2core_orbits"]
 assert [len(flt_g[2]),len(flt_f[2])]==camp["partial_split_primes"];msg=msf=0.0
 for p in ps:
  if p<=256 or p%15!=1:continue
  for om in omr[p]:
   if pow((-2*om*om)%p,(p-1)//3,p)==1 and (p,om) not in fb_g:msg+=1.0/p
   if pow((2+om)%p,(p-1)//3,p)==1 and (p,om) not in fb_f:msf+=1.0/p
 def dg(nv,u,v,flt):
  pao,pp,pl=flt
  if math.gcd(nv,pao)!=1:return False
  g=math.gcd(nv,pp)
  if g==1:return True
  for p,rm in pl:
   if g%p==0:
    if (u+v*rm)%p==0:return False
    while g%p==0:g//=p
    if g==1:break
  return True
 n_pts=int(camp["ball_pts"]);_yp0=_ypb[0]
 def _rkh(a,b):
  av=(a[0]*_yp0[0]+a[1]*_yp0[1]+a[2]*_yp0[2]+a[3]*_yp0[3]+a[4]*_yp0[4]+a[5]*_yp0[5]+a[6]*_yp0[6]+a[7]*_yp0[7])%Q_BIG
  bv=(b[0]*_yp0[0]+b[1]*_yp0[1]+b[2]*_yp0[2]+b[3]*_yp0[3]+b[4]*_yp0[4]+b[5]*_yp0[5]+b[6]*_yp0[6]+b[7]*_yp0[7])%Q_BIG
  return (-1 if av else-2) if bv==0 else(av*pow(bv,-1,Q_BIG))%Q_BIG
 _e1,_e2,_e3=[],[],[]
 for i in range(8):
  for s in (1,-1):v=[0]*8;v[i]=s;_e1.append(v)
 for i in range(8):
  for j in range(i+1,8):
   for si in (1,-1):
    for sj in (1,-1):v=[0]*8;v[i]=si;v[j]=sj;_e2.append(v)
 for i in range(8):
  for j in range(i+1,8):
   for k in range(j+1,8):
    for si in (1,-1):
     for sj in (1,-1):
      for sk in (1,-1):v=[0]*8;v[i]=si;v[j]=sj;v[k]=sk;_e3.append(v)
 _one=[1]+[0]*7;_zer=[0]*8;ball_pts=[];_sr=set()
 def _addp(a,b):
  if len(ball_pts)>=n_pts:return
  r=_rkh(a,b)
  if r!=-2 and r not in _sr:_sr.add(r);ball_pts.append((a,b))
 _addp(_one,_zer);_addp(_zer,_one)
 for b in _e1+_e2+_e3:_addp(_one,b)
 for a in _e2+_e3:_addp(a,_one)
 for a in _e2:
  for b in _e2:
   if len(ball_pts)>=n_pts:break
   _addp(a,b)
  if len(ball_pts)>=n_pts:break
 assert len(ball_pts)==n_pts and len(_sr)==n_pts;S_b=math.log2(n_pts)
 bsp=[([(i,v) for i,v in enumerate(al) if v],[(i,v) for i,v in enumerate(bt) if v]) for al,bt in ball_pts];_E8=[[1 if j==i else 0 for j in range(8)] for i in range(8)]
 def sq16_b(q,y0,r0):
  B=[[0]*16 for _ in range(16)];B[0][0]=q
  for i in range(1,8):B[i][0]=(-pow(y0,i,q))%q;B[i][i]=1
  for i in range(8):B[8+i][0]=(r0*pow(y0,i,q))%q;B[8+i][8+i]=1
  return lll_16d_int(B,99,100)[0]
 def ext_rt(pc,cs,a,b):
  z=d9(pc)
  for k in (1,2,4,7,8,11,13,14):
   y0=pow(z,k,pc);bs=sum(b[j]*pow(y0,j,pc) for j in range(8))%pc
   if bs==0:continue
   as_=sum(a[j]*pow(y0,j,pc) for j in range(8))%pc;r0=(as_*pow(bs,-1,pc))%pc;om=pow(y0,5,pc);rhs=(-2*om*om)%pc if cs==0 else(2+om)%pc
   if pow(r0,3,pc)==rhs:return (y0,r0)
  return None
 def au(side,q,y0,r0,nc_b,stop_first=False):
  R=sq16_b(q,y0,r0);nc_v=int(2**nc_b);assert nc_v<B_FB*B_FB;a1,b1=R[0][:8],R[0][8:];a2=b2=None
  for idx in range(1,16):
   ak,bk=R[idx][:8],R[idx][8:]
   if mul15(a1,bk)!=mul15(b1,ak):a2,b2=ak,bk;break
  assert a2 is not None;a1s=[mul15(e,a1) for e in _E8];a2s=[mul15(e,a2) for e in _E8];b1s=[mul15(e,b1) for e in _E8];b2s=[mul15(e,b2) for e in _E8]
  step=64 if stop_first else n_pts;M_use=sv_2=0;fch=None;seen_b=set()
  for off in range(0,n_pts,step):
   pts,ngs2,uvg2,ab3_2=[],[],[],[]
   for sp_idx in range(off,min(n_pts,off+step)):
    al_nz,bt_nz=bsp[sp_idx];a=[0]*8;b=[0]*8
    for i,s in al_nz:
     v1,v2=a1s[i],b1s[i]
     if s==1:
      for t in range(8):a[t]+=v1[t];b[t]+=v2[t]
     else:
      for t in range(8):a[t]-=v1[t];b[t]-=v2[t]
    for i,s in bt_nz:
     v1,v2=a2s[i],b2s[i]
     if s==1:
      for t in range(8):a[t]+=v1[t];b[t]+=v2[t]
     else:
      for t in range(8):a[t]-=v1[t];b[t]-=v2[t]
    if not any(b):continue
    a3,b3,(ug,vg,ng)=_gonly(a,b)
    if ng==0:continue
    if side==0:assert ng%q==0;ng//=q
    pts.append((a,b));ngs2.append(ng);uvg2.append((ug,vg));ab3_2.append((a3,b3))
   rgs2=cw(ngs2);kp=[i for i,x in enumerate(rgs2) if x<=nc_v];kpf=[];kidx=[]
   for i in kp:
    uf,vf,nf=_fonly(*ab3_2[i])
    if nf==0:continue
    if side==1:assert nf%q==0;nf//=q
    kpf.append((nf,uf,vf));kidx.append(i)
   rfs2=cw([x[0] for x in kpf])
   for j,idx in enumerate(kidx):
    rf=rfs2[j]
    if rf<=nc_v:
     sv_2+=1;rg=rgs2[idx];ug,vg=uvg2[idx];nf_v,uf,vf=kpf[j]
     if dg(ngs2[idx],ug,vg,flt_g) and dg(nf_v,uf,vf,flt_f):
      a,b=pts[idx];rk=rtk(a,b)
      if rk is None or rk in seen_b:continue
      seen_b.add(rk);M_use+=1
      if M_use==1:
       cands=[]
       if rf>B_FB and rf%15==1:
        rt=ext_rt(rf,1,a,b)
        if rt is not None:cands.append((1,rf,rt[0],rt[1]))
       if rg>B_FB and rg%15==1:
        rt=ext_rt(rg,0,a,b)
        if rt is not None:cands.append((0,rg,rt[0],rt[1]))
       cands.sort(key=lambda x:(-x[1],-x[0]));fch=cands
       if stop_first:return M_use,sv_2,fch
  return M_use,sv_2,(fch or [])
 qb,qbs,qbl=camp["qbins"],camp["qbins_sub"],camp["qbins_lo"]
 def d2(clo,chi,nc_c):
  def ss_t(m,alpha):
   mesh=48;h=(nc_c-a8)/mesh;wB=[math.log1p(h/(a8+i*h)) for i in range(mesh)];tot=0.0
   for jB,row in ((0,[1.0]),(1,wB)):
    for idx,c in enumerate(row):
     if not c:continue
     rem=m-jB*a8-idx*h
     if rem<0:continue
     u=(rem+alpha)/a8;tot+=c*(rho_db(u,max(float(a8),rem)) if u>2.0 else rho_exact(u))
   return min(1.0,tot)
  T_c=nc_c;c_hi,c_lo=(chi+S_b)/16.0,(clo+S_b)/16.0;mu_u,sv_up=[],[]
  for side in (0,1):
   pr_c=sum(ss_t(24.0*c_hi+cf-(clo if side==1 else 0.0),af_t)*ss_t(24.0*c_hi+cg-(clo if side==0 else 0.0),ag_t) for cf,cg in qb)/6.0
   pr_i=sum(ss_t(24.0*c_hi+cf-(clo if side==1 else 0.0),af_t)*ss_t(24.0*c_hi+cg-(clo if side==0 else 0.0),ag_t) for cf,cg in qbs)/24.0;pr=min(pr_i,4.0*pr_c)
   nf_hi=24.0*c_hi+qb[-1][0]-(clo if side==1 else 0.0);ng_hi=24.0*c_hi+qb[-1][1]-(clo if side==0 else 0.0);miss_tot=(nf_hi/a8)*msf+(ng_hi/a8)*msg;sv=0.0
   for cfl,cgl in qbl:
    nfl=24.0*c_lo+cfl-(chi if side==1 else 0.0);ngl=24.0*c_lo+cgl-(chi if side==0 else 0.0);sf=1.0 if T_c>=nfl else min(1.0,max(1.0,T_c/a8)*rho_exact((math.floor(nfl)-T_c)/a8))
    sg=1.0 if T_c>=ngl else min(1.0,max(1.0,T_c/a8)*rho_exact((math.floor(ngl)-T_c)/a8));sv+=sf*sg/6.0
   mu_u.append(n_pts*pr*math.exp(-miss_tot));sv_up.append(n_pts*min(1.0,sv+pr))
  return mu_u,sv_up
 mu1,svu1=d2(d1_lo,d1_hi,nc_1);mu2,svu2=d2(d2_lo,d2_hi,nc_2);l2h=camp["live_2hop"][:1];lct=camp["live_ctrls"][:1]
 for sp,qp,y0p,r0p,sc1,qc1,eM1,esv1,sc2,qc2,eM2,esv2 in l2h:
  M0,_,chs1=au(sp,qp,y0p,r0p,nc_0,stop_first=True);assert M0==1 and chs1 and chs1[0][0]==sc1 and chs1[0][1]==qc1;_,__,y01,r01=chs1[0]
  M1,sv1,chs2=au(sc1,qc1,y01,r01,nc_1,stop_first=False);assert M1==eM1 and sv1==esv1 and sv1<svu1[sc1] and chs2 and chs2[0][0]==sc2 and chs2[0][1]==qc2;_,__,y02,r02=chs2[0]
  M2,sv2,_=au(sc2,qc2,y02,r02,nc_2,stop_first=False);assert M2==eM2 and sv2==esv2 and sv2<svu2[sc2]
 for sc,qc,y0c,r0c,eM,esv in lct:
  M,sv,_=au(sc,qc,y0c,r0c,nc_1,stop_first=False);assert M==eM and sv==esv and sv<svu1[sc]
 def cx(rows):
  n=len(rows);mus=[mu2[r[1]] if r[0]==2 else mu1[r[1]] for r in rows];svus=[svu2[r[1]] if r[0]==2 else svu1[r[1]] for r in rows];ms=[r[2] for r in rows];sv2s=[r[3] for r in rows]
  nz=sum(1 for m in ms if m==0);hpgf_p0=sum((1.0+0.5*u)**(-1.0) for u in mus)/n;rats=[m/u for m,u in zip(ms,mus)];m_rat=sum(rats)/n;v_rat=sum((x-m_rat)**2 for x in rats)/(n-1)
  de=sum(sv2s)/sum(svus);df=max(s/su for s,su in zip(sv2s,svus));pgf_r1,pgf_r05={},{}
  for t in (0.0,0.2,0.5,0.8):
   eg=sum(t**m for m in ms)/n;mr1=sum((1.0+0.5*u*(1.0-t))**(-1.0) for u in mus)/n;mr05=sum((1.0+u*(1.0-t))**(-0.5) for u in mus)/n;assert eg<mr1 and eg<mr05
   pgf_r1[str(t)]=round(eg/mr1,6);pgf_r05[str(t)]=round(eg/mr05,6)
  assert nz==0 and min(ms)>=1 and m_rat>=1.0 and de<1.0 and df<1.0
  return dict(n=n,M_min=min(ms),M_max=max(ms),n_zero=nz,hpgf_r1_p0=round(hpgf_p0,5),M_over_mu_env=round(m_rat,4),var_M_over_mu=round(v_rat,4),sv2_over_upper=round(de,4),sv2_max_ratio=round(df,4),pgf_emp_over_hpgf_r1=pgf_r1,pgf_emp_over_hpgf_r05=pgf_r05,)
 ah=camp["all_records"];assert len(ah)==224;s_all=cx(ah);s_d1=cx([r for r in ah if r[0]==1]);s_d2=cx([r for r in ah if r[0]==2])
 d6=cx([r for r in ah if r[0]==0]);dj=round(1.0-(0.05**(1.0/len(ah))),5);d1_vs_ctrl=s_d1["M_over_mu_env"]/d6["M_over_mu_env"]
 d2_vs_ctrl=s_d2["M_over_mu_env"]/d6["M_over_mu_env"];assert dj<s_all["hpgf_r1_p0"] and d1_vs_ctrl>=0.90 and d2_vs_ctrl>=0.90
 return dict(tower=camp["tower"],same_field_q15=True,exact_tau_orbit_keyed=True,tau_orbit_bijection_verified_live=True,exact_kh_tau_dedup=True,rank2_ring_multiple_enum=True,stage1_unique_ab=len(ngs),real_2core_fb=True,literal_k1_first_usable=True,b_fb_bits=a8,raw_rels_2core=[len(rels),len(cur)],fb_2core_orbits=[len(fb_g),len(fb_f)],partial_split_primes=[len(flt_g[2]),len(flt_f[2])],ball_pts=n_pts,cp95_empty_ub=dj,mu_u_d1=[round(mu1[0],2),round(mu1[1],2)],mu_u_d2=[round(mu2[0],2),round(mu2[1],2)],sv_up_d1=[round(svu1[0],1),round(svu1[1],1)],sv_up_d2=[round(svu2[0],1),round(svu2[1],1)],live_2hop_reconstructed=len(l2h),live_1hop_reconstructed=len(l2h),live_ctrls_replayed=len(lct),live_balls_total=3*len(l2h)+len(lct),exact_live_match=True,d1_vs_ctrl=round(d1_vs_ctrl,4),d2_vs_ctrl=round(d2_vs_ctrl,4),all_224=s_all,d1_100=s_d1,d2_24=s_d2,ctrl_100=d6,checks_ok=True,)
def verify_prod_sq36(sq24,enc12,desc,dlive,U_bls,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn):
 import random
 assert len(sq24)==24 and len(enc12)==12;l2u=float(desc["params"]["log2_u8"]);lx=float(desc["params"]["lx"]);qb,qbl=desc["qbins"],desc["qbins_lo"];qbf,qbg=[x[1] for x in qb],[x[2] for x in qb]
 S8=math.log2(desc["cells_rank2"]["8"]);we=y_pow_fn(5);wp1=add_fn(we,[1]+[0]*7);w2=y_pow_fn(10);rng=random.Random(20261120);ab=[]
 while len(ab)<16:
  al=[rng.randint(-4,4) for _ in range(8)];bt=[rng.randint(-4,4) for _ in range(8)]
  if any(al) and sum(x*x for x in al)<=64 and sum(x*x for x in bt)<=64:ab.append((al,bt))
 def _ip(n):
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
 def _fg(a,b):
  a3=_kmul15(_kmul15(a,a),a);b3=_kmul15(_kmul15(b,b),b);ub=_kmul15(w2,[U_bls*x for x in b3]);ng=_knorm_zomega([a3[i]-ub[i] for i in range(8)])[2]
  a6,b6=_kmul15(a3,a3),_kmul15(b3,b3);a9,b9=_kmul15(a6,a3),_kmul15(b6,b3);t1,t2,t3=_kmul15(wp1,_kmul15(a6,b3)),_kmul15(a3,b6),_kmul15(we,b9)
  return _knorm_zomega([a9[i]+t1[i]+t2[i]-t3[i] for i in range(8)])[2],ng
 def _ec(recs,ie):
  os={}
  for si,sn in ((0,"g"),(1,"f")):
   cfs,cgs,sr=[],[],[];rr=[r for r in recs if r[0]==si];assert len(rr)==(6 if ie else 12)
   for rc in rr:
    if ie:
     sd,qh,pa,pb,v1,v2=rc;q=int(qh,16);assert sum(x*x for x in pa)+sum(x*x for x in pb)>(q**(1.0/8.0));pnf,png=_fg(pa,pb);assert (png if sd==0 else pnf)%q==0
    else:sd,qh,v1,v2=rc;q=int(qh,16)
    nb=q.bit_length()-1;assert 70<=nb<265 and q%15==1 and _ip(q);a1,b1,a2,b2=v1[:8],v1[8:],v2[:8],v2[8:];dq=abs(norm_fn(sub_fn(mul_fn(a1,b2),mul_fn(b1,a2))));assert dq>0 and dq%q==0
    cv=(math.log2(q)+S8)/16.0;pr=next(r for r in dlive["rows"] if r["side"]==sn and r["lo"]<=nb<r["hi"])
    Sb=math.log2(desc["cells_rank2"][str(pr["R"])]);clb=(pr["lo"]+Sb)/16.0;Tb=pr["kA"]*pr["nA"]+pr["kB"]*pr["n"]+96.0;sp=0.0
    for (wq,_,_),(cfl,cgl) in zip(qb,qbl):
     nfl=72.0*clb+cfl-(pr["hi"] if sn=="f" else 0.0);ngl=24.0*clb+l2u+cgl-(pr["hi"] if sn=="g" else 0.0);sf=1.0 if Tb>=nfl else min(1.0,max(1.0,Tb/lx)*rho_exact((math.floor(nfl)-Tb)/lx))
     sg=1.0 if Tb>=ngl else min(1.0,max(1.0,Tb/lx)*rho_exact((math.floor(ngl)-Tb)/lx));sp+=wq*sf*sg
    sa=0.0;ca=(math.log2(q)+Sb)/16.0;lq=math.log2(q)
    for al,bt in ab:
     aa1,aa2=_kmul15(al,a1),_kmul15(bt,a2);bb1,bb2=_kmul15(al,b1),_kmul15(bt,b2)
     a=[aa1[i]+aa2[i] for i in range(8)];b=[bb1[i]+bb2[i] for i in range(8)];nf,ng=_fg(a,b);assert (ng if si==0 else nf)%q==0
     cf=math.log2(nf)-72.0*cv;cg=math.log2(ng)-24.0*cv-l2u;cfs.append(cf);cgs.append(cg);nfa=72.0*ca+cf-(lq if sn=="f" else 0.0);nga=24.0*ca+l2u+cg-(lq if sn=="g" else 0.0)
     sfe=1.0 if Tb>=nfa else min(1.0,max(1.0,Tb/lx)*rho_exact((math.floor(nfa)-Tb)/lx));sge=1.0 if Tb>=nga else min(1.0,max(1.0,Tb/lx)*rho_exact((math.floor(nga)-Tb)/lx));sa+=sfe*sge
    sr.append((sa/len(ab))/max(sp,1e-300))
   n=len(cfs);mcf,mcg=sum(cfs)/n,sum(cgs)/n;cfs.sort();cgs.sort();sf6=[round(cfs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)];sg6=[round(cgs[min(n-1,int((i+1)*n/6)-1)],2) for i in range(6)]
   assert all(sf6[i]<qbf[i] for i in range(6)) and all(sg6[i]<qbg[i] for i in range(6)) and max(sr)<=1.0+1e-12 and (sum(sr)/len(sr))<0.50
   os[sn]=dict(bases=len(rr),samples=n,mean_cf=round(mcf,2),mean_cg=round(mcg,2),sext_f=sf6,sext_g=sg6,min_headroom_f=round(min(qbf[i]-sf6[i] for i in range(6)),2),min_headroom_g=round(min(qbg[i]-sg6[i] for i in range(6)),2),sv_emp_over_plan_mean=round(sum(sr)/len(sr),4),sv_emp_over_plan_max=round(max(sr),4))
  return os
 cs=_ec(sq24,False);es=_ec(enc12,True)
 d7=dict(g_delta_cf=round(es["g"]["mean_cf"]-cs["g"]["mean_cf"],2),g_delta_cg=round(es["g"]["mean_cg"]-cs["g"]["mean_cg"],2),f_delta_cf=round(es["f"]["mean_cf"]-cs["f"]["mean_cf"],2),f_delta_cg=round(es["f"]["mean_cg"]-cs["f"]["mean_cg"],2),parent_divisibility_verified="12/12")
 assert max(d7["g_delta_cf"],d7["g_delta_cg"],d7["f_delta_cf"],d7["f_delta_cg"])<=5.0;pw=1
 for _p in primes_upto(4096):pw*=_p
 bts=[[]]+[[(i,1)] for i in range(8)]+[[(i,-1)] for i in range(8)]
 for i in range(8):
  for j in range(i+1,8):
   for si in (1,-1):
    for sj in (1,-1):
     if len(bts)<64:bts.append([(i,si),(j,sj)])
 _E8=[[1 if j==i else 0 for j in range(8)] for i in range(8)]
 def _rw(recs,ie):
  rs=[]
  for rc in recs:
   sd=rc[0];q=int(rc[1],16);v1,v2=(rc[4],rc[5]) if ie else(rc[2],rc[3]);a1,b1,a2,b2=v1[:8],v1[8:],v2[:8],v2[8:];sfs,sgs,rfs,rgs=[],[],[],[]
   a2s=[_kmul15(e,a2) for e in _E8];b2s=[_kmul15(e,b2) for e in _E8]
   for bt in bts:
    a,b=list(a1),list(b1)
    for idx,sgn in bt:
     va,vb=a2s[idx],b2s[idx]
     if sgn==1:
      for t in range(8):a[t]+=va[t];b[t]+=vb[t]
     else:
      for t in range(8):a[t]-=va[t];b[t]-=vb[t]
    nf,ng=_fg(a,b)
    if sd==0:ng//=q
    else:nf//=q
    for val,dst,rdst in ((nf,sfs,rfs),(ng,sgs,rgs)):
     n0=val;g=math.gcd(val,pw)
     while g>1:val//=g;g=math.gcd(val,g)
     dst.append(math.log2(n0//val) if n0>val else 0.0);rdst.append(val.bit_length())
   rs.append((sfs,sgs,sum(rfs)/len(rfs),sum(rgs)/len(rgs)))
  return rs
 sc=_rw(sq24,False);se=_rw(enc12,True);nw={}
 for thr in (12,13,14,15):
  cc=[sum(1 for sf,sg in zip(x,y) if sf>=thr and sg>=thr) for x,y,_,_ in sc];ce=[sum(1 for sf,sg in zip(x,y) if sf>=thr and sg>=thr) for x,y,_,_ in se]
  ca=cc+ce;mu=sum(ca)/len(ca);vr=sum((x-mu)**2 for x in ca)/(len(ca)-1);assert min(ca)>=1
  pgr1={str(t):round((sum(t**m for m in ca)/len(ca))/((1.0+0.5*mu*(1.0-t))**(-1.0)),4) for t in (0.0,0.2,0.5,0.8)}
  pgr05={str(t):round((sum(t**m for m in ca)/len(ca))/((1.0+mu*(1.0-t))**(-0.5)),4) for t in (0.0,0.2,0.5,0.8)}
  assert all(v<1.0 for v in pgr1.values()) and all(v<1.0 for v in pgr05.values())
  nw[f"thr{thr}b"]=dict(n_zero=0,M_min=min(ca),mean=round(mu,2),enc_over_ctrl=round((sum(ce)/len(ce))/(sum(cc)/len(cc)),4),var=round(vr,2),pgf_over_hpgf_r1=pgr1,pgf_over_hpgf_r05=pgr05)
 rem_st=dict(ctrl_rem_f=round(sum(z[2] for z in sc)/len(sc),1),enc_rem_f=round(sum(z[2] for z in se)/len(se),1),ctrl_rem_g=round(sum(z[3] for z in sc)/len(sc),1),enc_rem_g=round(sum(z[3] for z in se)/len(se),1))
 return dict(n_bases=36,n_ctrl=24,n_enc=12,q_bit_range="71..258",ctrl_g=cs["g"],ctrl_f=cs["f"],enc_g=es["g"],enc_f=es["f"],live_7e_enc_vs_ctrl=d7,prod_f9_g3_hpgf_window=nw,prod_cofactor_rem_bits=rem_st,checks_ok=True)
def verify_gt_boot_20(data,P_bls,R_bls,U_bls,h_pi,w_pi,pmul_fn,prem_fn,y_pow_fn,add_fn,sub_fn,mul_fn,pow_fn,scale_fn,norm_fn):
 import hashlib
 db=data["descent"]["boot"];ay=int(db["s_boot"]);assert ay==(1<<21);w_pi2=(w_pi*w_pi)%P_bls;uw2=(U_bls*w_pi2)%P_bls;h0,h1,h2,h3=h_pi[0],h_pi[1],h_pi[2],h_pi[3]
 def a5(a,b):return [(a[0]+b[0])%P_bls,(a[1]+b[1])%P_bls,(a[2]+b[2])%P_bls,(a[3]+b[3])%P_bls]
 def ac(a,b):return [(a[0]-b[0])%P_bls,(a[1]-b[1])%P_bls,(a[2]-b[2])%P_bls,(a[3]-b[3])%P_bls]
 def ad(s,a):return [(s*a[0])%P_bls,(s*a[1])%P_bls,(s*a[2])%P_bls,(s*a[3])%P_bls]
 def a4(a,b):
  c=[a[0]*b[0],a[0]*b[1]+a[1]*b[0],a[0]*b[2]+a[1]*b[1]+a[2]*b[0],a[0]*b[3]+a[1]*b[2]+a[2]*b[1]+a[3]*b[0],a[1]*b[3]+a[2]*b[2]+a[3]*b[1],a[2]*b[3]+a[3]*b[2],a[3]*b[3]]
  for k in (6,5,4):
   v=c[k]%P_bls
   if v:
    o=k-4;c[o]-=v*h0;c[o+1]-=v*h1;c[o+2]-=v*h2;c[o+3]-=v*h3
  return [c[0]%P_bls,c[1]%P_bls,c[2]%P_bls,c[3]%P_bls]
 def f4_pow(a,e):
  res,b=[1,0,0,0],list(a)
  while e:
   if e&1:res=a4(res,b)
   e>>=1
   if e:b=a4(b,b)
  return res
 f1c=[f4_pow([1 if j==i else 0 for j in range(4)],P_bls) for i in range(4)]
 def aq(a):
  out=[0,0,0,0]
  for i in range(4):
   ai=a[i]
   if ai:
    col=f1c[i];out[0]+=ai*col[0];out[1]+=ai*col[1];out[2]+=ai*col[2];out[3]+=ai*col[3]
  return [x%P_bls for x in out]
 def an(a):return aq(aq(a))
 def f4_inv(a):
  a1=aq(a);a2=aq(a1);a3=aq(a2);pr3=a4(a1,a4(a2,a3))
  return ad(pow(a4(a,pr3)[0],-1,P_bls),pr3)
 def ct(A,B):
  a0,a1,a2=A;b0,b1,b2=B;v0=a4(a0,b0);v1=a4(a1,b1);v2=a4(a2,b2)
  t01=ac(ac(a4(a5(a0,a1),a5(b0,b1)),v0),v1);t02=ac(ac(a4(a5(a0,a2),a5(b0,b2)),v0),v2);t12=ac(ac(a4(a5(a1,a2),a5(b1,b2)),v1),v2)
  return (a5(v0,ad(uw2,t12)),a5(t01,ad(uw2,v2)),a5(t02,v1))
 def bd(A):
  a0,a1,a2=A;v0=a4(a0,a0);v1=a4(a1,a1);v2=a4(a2,a2);m01=a4(a0,a1);m02=a4(a0,a2);m12=a4(a1,a2)
  return (a5(v0,ad(uw2,a5(m12,m12))),a5(a5(m01,m01),ad(uw2,v2)),a5(a5(m02,m02),v1))
 G0=tuple([int(x,16) for x in c] for c in db["G0"]);one12=([1,0,0,0],[0,0,0,0],[0,0,0,0]);d0=[one12,G0]
 for _ in range(14):d0.append(ct(d0[-1],G0))
 def bm(e):
  res=one12
  for sh in range(((e.bit_length()-1)//4)*4,-1,-4):
   if res is not one12:res=bd(res);res=bd(res);res=bd(res);res=bd(res)
   nib=(e>>sh)&15
   if nib:res=d0[nib] if res is one12 else ct(res,d0[nib])
  return res
 assert G0!=one12 and bm(R_bls)==one12;inv2=(P_bls+1)//2;d_fp2=ac([0,1,0,0],an([0,1,0,0]));inv_2d=f4_inv(ad(2,d_fp2))
 def d4(a):
  ap=an(a)
  return ad(inv2,a5(a,ap)),a4(ac(a,ap),inv_2d)
 def elim_gamma(Z):
  e0,o0=d4(Z[0]);e1,o1=d4(Z[1]);e2,o2=d4(Z[2]);c0=ac(a4(e1,o0),a4(o1,e0));c1=ac(a4(o2,e0),a4(e2,o0));c2=ac(a4(e2,o1),a4(o2,e1))
  assert an(c0)==c0 and an(c1)==c1 and an(c2)==c2;W=ct((c0,c1,c2),Z);assert W[2]==[0,0,0,0] and any(W[1]);gam=a4([(-x)%P_bls for x in W[0]],f4_inv(W[1]))
  assert a5(W[0],a4(gam,W[1]))==[0,0,0,0]
  return gam
 z3_2=y_pow_fn(10);p4=P_bls**4;lov_l2=math.log2(((0.99-0.25)**(-15.0/4.0))*(float(P_bls)**0.25)*(float(ay)**0.5));bsub=data["descent"]["boot_sub_deciles"];m_sv=float(db["m_sv"])
 betas=[[0]*8]+[[1 if j==i else 0 for j in range(8)] for i in range(8)]+[[-1 if j==i else 0 for j in range(8)] for i in range(8)];_kmul=_kmul15
 def _knorm(a):return _knorm_zomega(a)[2]
 def dm(a,b):
  a3=_kmul(_kmul(a,a),a);b3=_kmul(_kmul(b,b),b);ub3=_kmul(z3_2,[U_bls*x for x in b3])
  return _knorm([a3[j]-ub3[j] for j in range(8)])
 az,b4,sk_l2s=[],[],[];ar=by=bq=0;ae=4
 for ti in range(ae):
  xt=int(hashlib.sha256(f"bls12-381-gt-boot-v5g-{ti}".encode()).hexdigest(),16)%(R_bls-1)+1;gam=elim_gamma(bm(xt));R_sk=run_boot_381bit_lll_live(gam,h_pi,P_bls,ay,pmul_fn,prem_fn)
  a1,b1=R_sk[0][:8],[x//ay for x in R_sk[0][8:]];sk_l2=math.log2(math.sqrt(sum(x*x for x in a1)+sum((ay*x)**2 for x in b1)));assert sk_l2<=lov_l2;sk_l2s.append(sk_l2)
  assert prem_fn(a1,h_pi,P_bls)==pmul_fn(gam,prem_fn(b1,h_pi,P_bls),P_bls,h_pi);g1=dm(a1,b1)
  if ti==0:assert g1==abs(norm_fn(sub_fn(pow_fn(a1,3),scale_fn(U_bls,mul_fn(z3_2,pow_fn(b1,3))))))
  lb1=math.log2(g1);az.append(lb1)
  if ti<4:
   a2=b2=None
   for idx in range(1,16):
    ak,bk=R_sk[idx][:8],[x//ay for x in R_sk[idx][8:]]
    if _knorm(sub_fn(_kmul(a1,bk),_kmul(b1,ak)))==p4:
     a2,b2=ak,bk;break
   assert a2 is not None and prem_fn(a2,h_pi,P_bls)==pmul_fn(gam,prem_fn(b2,h_pi,P_bls),P_bls,h_pi);bn,bb=[g1],[lb1]
   for bt in betas[1:]:
    aa=add_fn(a1,_kmul(bt,a2));bvec=add_fn(b1,_kmul(bt,b2));ng=dm(aa,bvec);lb=math.log2(ng);assert lb<2600.0;bn.append(ng);bb.append(lb)
   b4.extend(bb);bs_s=sorted(bb)
   if all(bs_s[min(16,math.ceil((j+1)*17/10)-1)]<=bsub[j] for j in range(10)):ar+=1
   for i in range(len(bn)):
    for j in range(i+1,len(bn)):
     bq+=1
     if math.gcd(bn[i],bn[j])>10000:by+=1
 n_pass=sum(1 for v in az if v<=m_sv);assert n_pass==ae and ar==4
 return dict(boot_element="a - b*X",rel_norm_formula="a^3 - u*omega^2*b^3",boot_gt_order_r_exact=True,boot_subfield_elim_live=ae,boot_s_boot=ay,boot_live_381bit_lll_exact_all=True,boot_supplied_lll_vectors=0,boot_skew_b1_max_l2=round(max(sk_l2s),2),boot_lovasz_bound_l2=round(lov_l2,2),boot_gt_pass=f"{n_pass}/{ae}",boot_det_compat_check=True,boot_det_p4_exact=True,boot_decile_filter_pass=f"{ar}/4",sv_mean_live=round(sum(az)/len(az),2),sv_max_live=round(max(az),2),ball_mean_live=round(sum(b4)/len(b4),2),ball_max_live=round(max(b4),2),gcd_gt1e4_live=f"{by}/{bq}",)
def unplanted_two_sided_bw_check():
 import random
 p,r,m=467,233,381;f_c,g_c=[1,1,-5,0,1],[-10,-1,1];assert sum(ci*pow(m,i,p) for i,ci in enumerate(f_c))%p==0;assert sum(ci*pow(m,i,p) for i,ci in enumerate(g_c))%p==0;eps_f,eps_g,r2=r-1,r*r-1,r*r
 def peval(c,x,mod):
  s=0
  for ci in reversed(c):
   s=(s*x+ci)%mod
  return s
 def pmul(a,b,mod_p,q):
  d=len(mod_p)-1;out=[0]*(len(a)+len(b)-1)
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
  d=len(mod_p)-1;res=[1]+[0]*(d-1);b=[x%q for x in base]+[0]*(d-len(base));e=exp
  while e:
   if e&1:res=pmul(res,b,mod_p,q)
   e>>=1
   if e:b=pmul(b,b,mod_p,q)
  return res
 bf=[q for q in range(2,115) if all(q%d for d in range(2,int(q**0.5)+1))]
 def unram(c,q):
  for rt in range(q):
   if peval(c,rt,q)==0:
    if sum(i*c[i]*pow(rt,i-1,q) for i in range(1,len(c)))%q==0:return False
  return True
 fb_f,fb_g={},{}
 for q in bf:
  if unram(f_c,q):
   for rt in range(q):
    if peval(f_c,rt,q)==0:fb_f[(q,rt)]=len(fb_f)
  if unram(g_c,q):
   for rt in range(q):
    if peval(g_c,rt,q)==0:fb_g[(q,rt)]=len(fb_g)
 pf_act=sorted({q for (q,_) in fb_f});pg_act=sorted({q for (q,_) in fb_g})
 def sfac(val,plist):
  n=abs(val)
  if n==0:return None
  fac={}
  for q in plist:
   if n%q==0:
    e=0
    while n%q==0:
     n//=q;e+=1
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
   vf={fb_f[(q,(a*pow(b,-1,q))%q)]:e for q,e in ff.items()};vg={fb_g[(q,(a*pow(b,-1,q))%q)]:e for q,e in fg.items()};rels.append((a,b,vf,vg))
 b7=len(fb_f);cur=list(range(len(rels)))
 while True:
  ccnt={}
  for ri in cur:
   for k in rels[ri][2]:
    ccnt[k]=ccnt.get(k,0)+1
   for k in rels[ri][3]:
    ccnt[b7+k]=ccnt.get(b7+k,0)+1
  nxt=[ri for ri in cur if all(ccnt[k]>=2 for k in rels[ri][2]) and all(ccnt[b7+k]>=2 for k in rels[ri][3])]
  if len(nxt)==len(cur):break
  cur=nxt
 used=sorted({k for ri in cur for k in rels[ri][2]}|{b7+k for ri in cur for k in rels[ri][3]});cmap={c:i for i,c in enumerate(used)};N_2core=len(used)
 def smap(a,b,mod_p,eps):
  el=ppow([a%r2,(-b)%r2],eps,mod_p,r2);assert el[0]%r==1 and all(x%r==0 for x in el[1:])
  return [((el[0]-1)//r)%r]+[(x//r)%r for x in el[1:]]
 rM,rows_S=[],[]
 for ri in cur:
  a,b,vf,vg=rels[ri];row=[0]*N_2core
  for k,e in vf.items():
   row[cmap[k]]=e%r
  for k,e in vg.items():
   row[cmap[b7+k]]=(-e)%r
  rM.append(row);rows_S.append(smap(a,b,f_c,eps_f)[1:4]+smap(a,b,g_c,eps_g)[1:2])
 mrw=[row[:] for row in rM];srw=[row[:] for row in rows_S];alr=[True]*len(cur);alc=[True]*N_2core;mst=[]
 while True:
  piv_c=next((c for c in range(N_2core) if alc[c] and sum(1 for i in range(len(cur)) if alr[i] and mrw[i][c]!=0)==2),-1)
  if piv_c<0:break
  r_pair=[i for i in range(len(cur)) if alr[i] and mrw[i][piv_c]!=0];rp,ro=r_pair[0],r_pair[1];mst.append((piv_c,mrw[rp][:],srw[rp][:]))
  fac=(mrw[ro][piv_c]*pow(mrw[rp][piv_c],-1,r))%r;mrw[ro]=[(mrw[ro][j]-fac*mrw[rp][j])%r for j in range(N_2core)];srw[ro]=[(srw[ro][j]-fac*srw[rp][j])%r for j in range(4)]
  alr[rp]=False;alc[piv_c]=False
 ccl=[c for c in range(N_2core) if alc[c]];N,D=len(ccl),len(ccl)+4;ax=[i for i in range(len(cur)) if alr[i]]
 cc={c:sum(1 for i in ax if mrw[i][c]!=0) for c in ccl}
 while len(ax)>D:
  cands=[i for i in ax if all(cc[c]>=3 for c in ccl if mrw[i][c]!=0)] or [i for i in ax if all(cc[c]>=2 for c in ccl if mrw[i][c]!=0)]
  drop=max(cands,key=lambda i:(sum(1 for c in ccl if mrw[i][c]!=0),i));ax.remove(drop)
  for c in ccl:
   if mrw[drop][c]!=0:cc[c]-=1
 def a7(mat):
  mr=[row[:] for row in mat];R,C=len(mr),len(mr[0]);pivs,rk=[],0
  for c in range(C):
   piv=next((i for i in range(rk,R) if mr[i][c]%r!=0),-1)
   if piv<0:continue
   mr[rk],mr[piv]=mr[piv],mr[rk];inv=pow(mr[rk][c],-1,r);rr=mr[rk]
   for j in range(c,C):
    rr[j]=(rr[j]*inv)%r
   for i in range(R):
    if i!=rk and mr[i][c]:
     fv=mr[i][c];ri_r=mr[i]
     for j in range(c,C):
      ri_r[j]=(ri_r[j]-fv*rr[j])%r
   pivs.append(c);rk+=1
  pset=set(pivs);vecs=[]
  for fc in [c for c in range(C) if c not in pset]:
   v=[0]*C;v[fc]=1
   for ri_i,pc in enumerate(pivs):
    v[pc]=(-mr[ri_i][fc])%r
   vecs.append(v)
  return rk,vecs
 M_sq=[[mrw[i][c] for c in ccl] for i in ax];S_sq=[srw[i][:] for i in ax];A_sq=[M_sq[i]+S_sq[i] for i in range(D)];rk_M,_=a7(M_sq);rk_A,ker_A=a7(A_sq)
 rk_3sm,ker_3sm=a7([M_sq[i]+S_sq[i][:3] for i in range(D)]);assert rk_M==N and rk_A==D-1 and len(ker_A)==1 and len(ker_3sm)==0
 def mmul(P,Q):
  R,K,C=len(P),len(P[0]),len(Q[0]);out=[[0]*C for _ in range(R)]
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
 def _sfr(rg,m,kmax=64):
  b=m.bit_length()
  for _ in range(kmax):
   x=rg.getrandbits(b)
   if 1<=x<m:return x
  raise AssertionError("F_r rejection cap")
 rng=random.Random(20261119);m_bw=n_bw=4;bnds=[(k*D)//4 for k in range(5)];cy=[[0]*D for _ in range(m_bw)]
 for k in range(m_bw):
  for j in range(bnds[k],bnds[k+1]):
   cy[k][j]=_sfr(rng,r)
 Y=[[_sfr(rng,r) for _ in range(n_bw)] for _ in range(D)];hr,hc=(D//m_bw)+2,(D//n_bw)+2;AY,seq=[row[:] for row in Y],[]
 for _ in range(hr+hc+2):
  seq.append(mmul(cy,AY));AY=mmul(A_sq,AY)
 H_rows=[[v for j in range(hc) for v in seq[i+j+1][ri_b]] for i in range(hr) for ri_b in range(m_bw)];rk_H,ker_H=a7(H_rows);v_rec=None
 for c_vec in ker_H:
  acc=[0]*D
  for j in range(hc-1,-1,-1):
   cj=c_vec[j*n_bw:(j+1)*n_bw];yc=[sum(Y[ri_i][bi]*cj[bi] for bi in range(n_bw))%r for ri_i in range(D)];acc=[(a_i+y_i)%r for a_i,y_i in zip(mvec(A_sq,acc),yc)]
  if any(acc) and all(x==0 for x in mvec(A_sq,acc)):
   v_rec=acc
   break
 assert v_rec is not None;O_rows,cX=[],[row[:] for row in cy]
 for _ in range(hr):
  cX=mmul(cX,A_sq);O_rows.extend(cX)
 rk_OA,_=a7(O_rows);C_cols,cY=[],[row[:] for row in Y]
 for _ in range(hc):
  for bi in range(n_bw):
   C_cols.append([cY[ri_i][bi] for ri_i in range(D)])
  cY=mmul(A_sq,cY)
 rk_Cd,_=a7(C_cols);X_bad=[[((k+1)*cy[0][j])%r for j in range(D)] for k in range(m_bw)];O_bad,cXb=[],[row[:] for row in X_bad]
 for _ in range(hr):
  cXb=mmul(cXb,A_sq);O_bad.extend(cXb)
 bw,_=a7(O_bad);vlf={c:v_rec[j] for j,c in enumerate(ccl)};cl=v_rec[N:]
 for piv_c,pm,ps in reversed(mst):
  acc=(sum(pm[j]*vlf[j] for j in range(N_2core) if j!=piv_c and pm[j]!=0)+sum(ps[s]*cl[s] for s in range(4)))%r;vlf[piv_c]=(-acc*pow(pm[piv_c],-1,r))%r
 while True:
  prog=False
  for idx in range(len(cur)):
   unk=[j for j in range(N_2core) if rM[idx][j]!=0 and j not in vlf]
   if len(unk)==1:
    cn=unk[0];acc=(sum(rM[idx][j]*vlf[j] for j in range(N_2core) if j!=cn and rM[idx][j]!=0)+sum(rows_S[idx][s]*cl[s] for s in range(4)))%r
    vlf[cn]=(-acc*pow(rM[idx][cn],-1,r))%r;prog=True
  if not prog:break
 assert len(vlf)==N_2core;gen=next(g for g in range(2,p) if pow(g,(p-1)//r,p)!=1 and pow(g,2,p)!=1);g_r,dlog,ag=pow(gen,(p-1)//r,p),{},1
 for e in range(r):
  dlog[ag]=e;ag=(ag*g_r)%p
 d8=set()
 for idx,ri in enumerate(cur):
  a,b,vf,vg=rels[ri];lhs=(sum(e*vlf[cmap[k]] for k,e in vf.items())+sum(rows_S[idx][j]*cl[j] for j in range(3)))%r
  rhs=(sum(e*vlf[cmap[b7+k]] for k,e in vg.items())-rows_S[idx][3]*cl[3])%r;assert lhs==rhs;dl=dlog[pow((a-b*m)%p,(p-1)//r,p)]
  if dl!=0:d8.add((lhs*pow(dl,-1,r))%r)
 assert rk_H==D-1 and rk_OA==D-1 and rk_Cd==D and bw<D-1 and len(d8)==1
 return dict(p=p,r=r,raw_rels=len(rels),purged_rels=len(cur),N_2core=N_2core,n_merged_cols=len(mst),N=N,D=D,rank_M=rk_M,rank_A=rk_A,ker_dim_A=len(ker_A),ker_dim_3sm_ctrl=len(ker_3sm),rank_Cd=rk_Cd,rank_OeA=rk_OA,rank_H=rk_H,rank_OeA_1block_ctrl=bw,unplanted_2core_backsub_exact=True,horner_all_purged_rels_ok=True,distinct_Fp_dlp_ratios=len(d8),checks_ok=True,)
