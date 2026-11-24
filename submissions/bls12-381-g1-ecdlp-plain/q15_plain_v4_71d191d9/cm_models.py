"""cm_models.py — Self-contained."""
import array
import collections
import functools
import math
from decimal import Decimal,localcontext
LN2=math.log(2.0)
_SUMPROD=getattr(math,"sumprod",lambda a,b:sum(x*y for x,y in zip(a,b)))
GAMMA=0.5772156649015329
MONT5=280.0
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
 while len(rv)>1 and abs(rv[0])<1e-18*abs(r[0]):rv.pop(0)
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
 _sp=_SUMPROD
 jm=length
 while jm>0 and weights[jm]==0.0:jm-=1
 if jm==0:return[math.exp(-normalizer)]+[0.0]*length
 j_lo=1
 while j_lo<=jm and weights[j_lo]==0.0:j_lo+=1
 jw_rev=[j*weights[j]for j in range(jm,j_lo-1,-1)]
 span=jm-j_lo+1
 out=[math.exp(-normalizer)]+[0.0]*length
 for n in range(j_lo,length+1):
  k=n if n<jm else jm
  cnt=k-j_lo+1
  val=_sp(jw_rev[span-cnt:span],out[n-k:n-j_lo+1])
  out[n]=val/n
  if val<-1e-13:raise ArithmeticError("negative coefficient beyond roundoff tolerance")
 return out
class Retention:
 def __init__(self,B_bits,X_bits,*,prime_min=2**24,step=0.25,max_bits=760,tilt=0.05,explicit_bits=14,power_terms=8,odd=1):
  self.B,self.X,self.minimum,self.step=B_bits,X_bits,prime_min,step
  self._lmin=math.log2(prime_min)
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
   _exp=math.exp;_len=self.length;_pt=power_terms
   if pref=="small":
    ws=weights["small"];ts=totals["small"]
    for i in range(cells):
     b=left+(i+0.5)*w;bi=round(b/step);idx=bi
     for k in range(1,_pt+1):
      ex=k*tilt*bi*step-(k-1)*LN2*b
      if ex<-75.0 and idx>_len:break
      amt=_exp(ex)*w/(k*b);ts+=amt
      if idx<=_len:ws[idx]+=amt
      idx+=bi
    totals["small"]=ts
   else:
    wg=weights[pref+"_geo"];wf=weights[pref+"_sf"];tg=totals[pref+"_geo"];tf=totals[pref+"_sf"]
    for i in range(cells):
     b=left+(i+0.5)*w;bi=round(b/step);idx=bi
     for k in range(1,_pt+1):
      ex=k*tilt*bi*step-(k-1)*LN2*b
      if ex<-75.0 and idx>_len:break
      amt=_exp(ex)*w/(k*b);tg+=amt;sa=-amt if(k&1)==0 else amt;tf+=sa
      if idx<=_len:wg[idx]+=amt;wf[idx]+=sa
      idx+=bi
    totals[pref+"_geo"]=tg;totals[pref+"_sf"]=tf
  lm=self._lmin
  continuum(cutoff,min(lm,B_bits),"small")
  continuum(max(cutoff,lm),min(X_bits,B_bits),"sieve")
  continuum(max(cutoff,X_bits),B_bits,"large")
  self.weights,self.normalizers=weights,dict(totals)
  all_w=[sum(weights[k][j]for k in ("small","sieve_geo","large_geo"))for j in range(self.length+1)]
  self.denominator=_ret_coeffs(all_w,sum(totals[k]for k in ("small","sieve_geo","large_geo")),self.length)
  self.sel_norm=math.exp(totals["sieve_sf"]+totals["large_sf"]-totals["sieve_geo"]-totals["large_geo"])
  self._c={}
  self._bkeys={False:[],True:[]}
 def co(self,k):
  v=self._c.get(k)
  if v is None:
   ln=self.length if k.startswith("sieve") else(self.length if self.length<2048 else 2048)
   v=_ret_coeffs(self.weights[k],self.normalizers[k],ln)
   self._c[k]=v
  return v
 def two(self,n,a,b=9**9,win=None,sf=False):
  a=n if n<a else a
  if a<0 or b<0:return 0.0
  c=self._c
  D=self.length if self.length<2048 else 2048
  lim=a if a<D else D
  lim_b=b if b<D else D
  if win is None and lim_b>=lim:lim_b=D
  mk=(lim_b if win is None else b,win,sf)
  m=c.get(mk)
  if m is None:
   _sp=_SUMPROD
   sm=self.co("small")
   lg=self.co("large_sf" if sf else "large_geo")
   if win is None:
    import bisect as _bs
    bks=self._bkeys[sf]
    ip=_bs.bisect_left(bks,lim_b)
    if ip>0:
     pb=bks[ip-1]
     m=list(c[(pb,None,sf)])
     for j in range(pb+1,lim_b+1):
      cb=lg[j]
      if cb!=0.0:
       for s in range(j,D+1):m[s]+=cb*sm[s-j]
    else:
     lg_f=lg[:lim_b+1]
     while len(lg_f)>1 and lg_f[-1]==0.0:lg_f.pop()
     Llg=len(lg_f);l0=lg_f[0]
     l_min=next((i for i in range(1,Llg)if lg_f[i]),Llg)
     m=[l0*sm[s]for s in range(D+1)]
     if l_min<Llg:
      lg_tail_rev=lg_f[l_min:][::-1];Lt=len(lg_tail_rev)
      for s in range(l_min,D+1):
       hi=s if s<Llg-1 else Llg-1
       m[s]+=_sp(lg_tail_rev[Lt-1-(hi-l_min):Lt],sm[s-hi:s-l_min+1])
    c[mk]=m
    _bs.insort(bks,lim_b)
   else:
    st=self.step;w0,w1,w2=win
    lg_f=[v if(v and(l*st<=w2 or any(self.X*k-w2<=l*st<=w1*k+w2 for k in range(1,w0+1))))else 0.0 for l,v in enumerate(lg[:lim_b+1])]
    while len(lg_f)>1 and lg_f[-1]==0.0:lg_f.pop()
    Llg=len(lg_f);l0=lg_f[0]
    l_min=next((i for i in range(1,Llg)if lg_f[i]),Llg)
    m=[l0*sm[s]for s in range(D+1)]
    if l_min<Llg:
     lg_tail_rev=lg_f[l_min:][::-1];Lt=len(lg_tail_rev)
     for s in range(l_min,D+1):
      hi=s if s<Llg-1 else Llg-1
      m[s]+=_sp(lg_tail_rev[Lt-1-(hi-l_min):Lt],sm[s-hi:s-l_min+1])
    c[mk]=m
  sk="sieve_sf_rev" if sf else "sieve_geo_rev"
  sv_r=c.get(sk)
  if sv_r is None:
   sv_r=self.co("sieve_sf" if sf else "sieve_geo")[::-1]
   c[sk]=sv_r
  off=self.length-n
  if win is None:
   tk=("two",n,lim,lim_b,sf)
   tv=c.get(tk)
   if tv is not None:return tv
   pv=c.get(("two",n,lim-1,lim_b,sf))
   if pv is not None:
    tv=pv+m[lim]*sv_r[off+lim]
    c[tk]=tv
    return tv
   tv=_SUMPROD(m[:lim+1],sv_r[off:off+lim+1])
   c[tk]=tv
   return tv
  return _SUMPROD(m[:lim+1],sv_r[off:off+lim+1])
 def cap(self,nb):
  return math.floor(math.ceil(nb)/self._lmin)
 def at(self,nb,sl,sb=None):
  st=self.step
  n=round(nb/st)
  rb=(math.floor(math.ceil(nb)/self._lmin)+1)*0.0625
  a=math.floor((sl-rb)/st)
  if a<0 or n<0:return 0.0
  lim=a if a<n else n
  b=9**9 if sb is None else math.floor((sb-rb-0.125)/st)
  if b<0:return 0.0
  if b>=lim:b=9**9
  ck=("at",n,lim,b)
  v=self._c.get(ck)
  if v is None:
   r=self.sel_norm*self.two(n,lim,b,None,True)/self.denominator[n]
   v=0.0 if r<0.0 else(1.0 if r>1.0 else r)
   self._c[ck]=v
  return v
 def survivor_upper(self,nb,sl,ex=0.0,sb=None,lb=None,kmax=4):
  n=round(nb/self.step)
  a=math.floor((sl+ex)/self.step)
  if sb is None:return self.two(n,a)/self.denominator[n]
  return self.two(n,a,math.floor((sb+0.25)/self.step),(kmax,lb,(self.cap(nb)+2)/16.0+0.25))/self.denominator[n]
@functools.lru_cache(None)
def ecm_instr(fbits,nbits,success=0.999999,lo_r=48.0):
 import q15_ecm16 as _EC
 return _EC.ladder_instr(lo_r,float(fbits),int(math.ceil(nbits)),success=success,step=12.0,b_min=96,l_min=2)
def ecm_mean_var(fbits,nbits,success=0.999999,lo_r=48.0):
 """Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree post-2^48 collection ECM schedule (H-ECM).
 Fixed deterministic curve count per cofactor with Stage-1 stepped gcd and Stage-2 J-leaf binary-product-tree
 giant-step bisection on every curve (ZERO curve replay, ZERO p1-as-upper-bound, ZERO variance => var_t = 0.0);
 per-relation survival discount 1 - 29*(1 - success) is charged directly inside collection yield."""
 import q15_ecm16 as _EC
 cost=_EC.ladder_instr(lo_r,float(fbits),int(math.ceil(nbits)),success=success,step=12.0,b_min=96,l_min=2)
 return (cost,0.0)
NTT4_PRIMES=(269*(1<<53)+1,276*(1<<53)+1,305*(1<<53)+1,309*(1<<53)+1,)
NTT4_PROD=NTT4_PRIMES[0]*NTT4_PRIMES[1]*NTT4_PRIMES[2]*NTT4_PRIMES[3]
NTT3_PRIMES=(69*(1<<55)+1,57*(1<<55)+1,54*(1<<55)+1)
NTT3_PROD=NTT3_PRIMES[0]*NTT3_PRIMES[1]*NTT3_PRIMES[2]
def ntt3_F(M_prod):
 m=math.ceil(max(M_prod,1.0))
 if m<=1:return 220.0
 if m>(1<<55):return 3.0*ntt3_F(math.ceil(m/2.0))+20.0*m
 k=max(2,(m-1).bit_length());n2=float(1<<k)
 c2=55.25*n2*k+(26.0*n2 if(k&1)else 0.0)+220.0*n2+(512.0*n2 if n2>(1<<26)else 0.0)
 n3=3.0*float(1<<(k-2))
 if n3>=m and k>=3:
  c3=55.25*n3*(k-2)+(26.0*n3 if((k-2)&1)else 0.0)+364.0*n3+(512.0*n3 if n3>(1<<26)else 0.0)
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
 n=len(B);sp=_SUMPROD if ip_fn is None else ip_fn
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
SMALL_150=((16,2),(31,12),(61,36),(81,4),(121,12),(151,24),(181,24),(211,24),(241,12))

def verify_kernels_symbolic(kernels_text,coll,la,r_bls,lll_node_cap,lll_boot_cap,desc_p=None,sga_mean=None):
 import re,random,struct
 allowed={"ldp","stp","ldr","str","ldrb","strb","ldrh","strh","ldur","stur","add","adds","adc","adcs","sub","subs","sbc","sbcs","mul","umulh","smulh","madd","msub","and","orr","eor","lsl","lsr","asr","bfi","cmp","tst","csel","mov","b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","b.ge","b.lt","b.le","clz","udiv","cbz","cbnz","ret"}
 forb=re.compile(r"\b(v[0-9]+|q[0-9]+|d[0-9]+|z[0-9]+|[xw]18|dc|ldadd\w*|swp\w*|cas\w*|bl|blr)\b")
 blocks=re.findall(r"```asm\n(.*?)```",kernels_text,flags=re.S)
 parsed={};all_defs=set();all_refs=set()
 for blk in blocks:
  lines=[l.strip()for l in blk.splitlines()if l.strip()]
  m=re.match(r"^;\s*kernel:\s*(\S+)\s+items=(\d+)\s+charged=([0-9.]+)\s+key=(\S+)$",lines[0])
  assert m,f"Bad kernel header: {lines[0]}"
  name,items,charged,key=m.group(1),int(m.group(2)),float(m.group(3)),m.group(4)
  for _l in lines[1:]:
   if _l.endswith(":"):all_defs.add(_l[:-1])
  insns=[l for l in lines[1:]if not l.startswith(";")and not l.endswith(":")]
  for ins in insns:
   assert not forb.search(ins.lower()),f"Forbidden SIMD/LSE/call in {name}: {ins}"
   ps=ins.replace(","," ").split();op=ps[0].lower()
   assert op in allowed,f"Disallowed opcode {op} in {name}: {ins}"
   if op in ("b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","b.ge","b.lt","b.le"):all_refs.add(ps[1])
   elif op in ("cbz","cbnz"):all_refs.add(ps[2])
   if op in ("ldp","stp"):
    mb=re.search(r"\[([^\]]+)\]",ins)
    assert mb and re.match(r"^\s*(?:x\w+|sp)(\s*,\s*#-?\d+)?\s*$",mb.group(1)),f"Illegal {op} addr mode in {name}: {ins}"
  ir=len(insns)/items
  assert ir<=charged+1e-12,f"{name}: {ir}>{charged}"
  parsed[name]=dict(raw_opcodes=len(insns),items=items,ir=ir,charged=charged,insns=insns)
 undef_labels=all_refs-all_defs
 assert not undef_labels,f"Undefined ASM labels in KERNELS.md: {sorted(undef_labels)}"
 _KERN_REG_MAPS={
  "rev_b2_row_k0_step4":{"x_st":"r19","x_p":"r20","x_d1":"r21","x_pms":"r22","x_b0":"r23"},
  "rev_b1_kt_step_shared":{"x_pms":"r21","x_st":"r22"},
  "rev_b1_bkt_hdr_u1":{"x_off":"r23","x_p":"r19","x_end":"r24"},
  "rev_b1_bkt_hdr_cnt":{"x_off":"r23","x_p":"r19","x_end":"r24","x_cnt":"r25"},
  "rev_b1_stream_wrap32":{"x_cnt":"r25"},
  "rev_b1_stream_wrap4":{"x_cnt":"r25"},
  "rev_b1_stream_u32":{"x_src":"r19","x_dst":"r20","x_cnt":"r25"},
  "rev_b1_stream_u4":{"x_src":"r19","x_dst":"r20","x_cnt":"r25"},
  "rev_b1_stream_tail1":{"x_src":"r19","x_dst":"r20","x_cnt":"r25"},
  "rev_b1a_apply_wrap32":{"x_cnt":"r25"},
  "rev_b1a_apply_wrap4":{"x_cnt":"r25"},
  "rev_b1a_apply_u32":{"x_p":"r19","x_acc":"r20","x_cnt":"r25","w_lp":"r26"},
  "rev_b1a_apply_u4":{"x_p":"r19","x_acc":"r20","x_cnt":"r25","w_lp":"r26"},
  "rev_b1a_apply_1":{"x_p":"r19","x_acc":"r20","x_cnt":"r25","w_lp":"r26"},
  "rev_b1b_eor_u4":{"x_p":"r19","x_acc":"r20","x_cnt":"r25","w_lp":"r26","x_bit57":"r27"},
  "rev_b1b_sub_u4":{"x_p":"r19","x_acc":"r20","x_k56":"r22","x_mask58":"r23","x_cnt":"r25","w_lp":"r26"},
  "rev_b1b_eor_1":{"x_p":"r19","x_acc":"r20","x_cnt":"r25","w_lp":"r26","x_bit57":"r27"},
  "rev_b1b_apply_1":{"x_p":"r19","x_acc":"r20","x_k56":"r22","x_mask58":"r23","x_cnt":"r25","w_lp":"r26"},
  "rev_sb_probe_fast":{"x_k0":"r7","x_p":"r8","x_lb":"r9","x_col0":"r15","x_krot56":"r16","x_off":"r12","x_rot":"r23","x_tx_tbl":"r25"},
  "rev_b2_direct_sb_hdr":{},
  "rev_b2_emit_c0":{"x_rot":"r23"},
  "rev_b2_emit_c1":{},
  "rev_b2_flush_s0":{"x_rx":"r12","x_rx_end":"r13","x_tx_tbl":"r25"},
  "rev_b2_flush_s1":{"x_rx":"r12","x_rx_end":"r13","x_tx_tbl":"r25"},
  "rev_scan_zero_u64":{"x_g":"r21","x_f":"r22","x_b16":"r23","x_8000":"r24","x_cnt":"r25"},
  "mont5_cios_round":{"x_b":"r6","x_n0":"r7"},
  "mont5_cios_epilogue":{"x_a":"r5","x_out":"r8","x_p":"r9"},
  "ntt62_butterfly":{"x_u":"r19","x_v":"r20","x_w":"r21","x_p":"r22","x_cnt":"r23"},
  "ntt62_r4_shoup":{"x_a0":"r19","x_a1":"r20","x_a2":"r21","x_a3":"r22","x_tw":"r23","x_p":"r24","x_p2":"r25","x_j":"r26","x_jp":"r27","x_cnt":"r28"},
  "row_red5_z12":{"x_out":"r21","x_mu":"r26"},
 }
 for _kn,_kd in parsed.items():
  _kmap=_KERN_REG_MAPS.get(_kn,{})
  _nums_k={f"r{m0}" for _ins in _kd["insns"] for m0 in re.findall(r"\b[xw](\d+)\b",_ins.lower())}
  _syms_k={m0 for _ins in _kd["insns"] for m0 in re.findall(r"\b([xw]_[a-z0-9_]+)\b",_ins.lower())}
  for _sk in _syms_k:
   assert _sk in _kmap,f"Unmapped symbolic register {_sk} in {_kn}!"
  _mapped_vals=[_kmap[_sk] for _sk in _syms_k]
  assert len(_mapped_vals)==len(set(_mapped_vals)),f"Symbolic-to-symbolic register collision in {_kn}: {_kmap}"
  assert set(_mapped_vals).isdisjoint(_nums_k),f"Symbolic-to-numbered register collision in {_kn}: {set(_mapped_vals)&_nums_k}"
  assert len(_nums_k|set(_mapped_vals))<=28 and "r18" not in (_nums_k|set(_mapped_vals)) and "r29" not in (_nums_k|set(_mapped_vals)),f"GPR limit exceeded in {_kn}"
 def _pw(kname,insns=None):
  _kmap=_KERN_REG_MAPS.get(kname,{})
  if insns is None:insns=parsed[kname]["insns"]
  ws=set()
  def _norm_r(tok):
   t=tok.lower()
   if t in ("xzr","wzr"):return None
   if re.match(r"^[wx]\d+$",t):return "r"+t[1:]
   if t in _kmap:return _kmap[t]
   if re.match(r"^[wx]_\w+$",t)or t=="sp":return t
   return None
  for ins in insns:
   wb=re.search(r"\[([xw]\w+|sp)(?:,\s*#-?\d+)?\](?:\s*,\s*#-?\d+|!)",ins.lower())
   if wb:
    nr=_norm_r(wb.group(1))
    if nr:ws.add(nr)
   ps=ins.replace(","," ").split();op=ps[0].lower()
   if op in ("cmp","tst","b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","b.ge","b.lt","b.le","cbz","cbnz","ret","str","strb","strh","stur","stp"):continue
   for r in [ps[1]]+(ps[2:3]if op=="ldp" else[]):
    nr=_norm_r(r)
    if nr:ws.add(nr)
  return ws
 assert _pw("rev_b2_row_k0_step4").isdisjoint({"r19","r20","r21","r22","r23"})
 assert parsed["rev_b2_row_k0_step4"]["raw_opcodes"]==33 and parsed["rev_b2_row_k0_step4"]["items"]==4 and (parsed["rev_b2_row_k0_step4"]["ir"]+256.0/189.0)<=float(coll.get("c_b2_row_k0",9.65))
 assert _pw("rev_b2_emit_c0").isdisjoint({"r3","r14","r19","r21","r23","r25"}) and _pw("rev_b2_emit_c1").isdisjoint({"r3","r14","r19","r21","r23","r25"})
 assert _pw("rev_b2_flush_s0").isdisjoint({"r2","r3","r11","r14","r19","r21","r23","r25"}) and _pw("rev_b2_flush_s1").isdisjoint({"r2","r3","r10","r14","r19","r21","r23","r25"})
 assert parsed["rev_sb_probe_fast"]["raw_opcodes"]<=23 and parsed["rev_sb_probe_fast"]["ir"]<=coll["c_tile_hdr"]and "cbz x_lb, .Lempty_chord" in parsed["rev_sb_probe_fast"]["insns"]and "asr x22, x22, #4" in parsed["rev_sb_probe_fast"]["insns"]and "ldp x10, x11, [x14]" in parsed["rev_sb_probe_fast"]["insns"]
 assert parsed["rev_b2_direct_sb_hdr"]["raw_opcodes"]==10 and (8+1)<=coll["c_sb_hdr"]and parsed["rev_b2_direct_sb_hdr"]["ir"]<=coll["c_sb_hdr"]and "cbz x21, .Lb2_sb_done" in parsed["rev_b2_direct_sb_hdr"]["insns"]and "str x10, [x14], #8" in parsed["rev_b2_direct_sb_hdr"]["insns"]and "ldr x11, [x14, #8]" in parsed["rev_b2_direct_sb_hdr"]["insns"]
 assert parsed["rev_b1_kt_step_shared"]["raw_opcodes"]==3 and parsed["rev_b1_bkt_hdr_cnt"]["raw_opcodes"]==4 and parsed["rev_b1_bkt_hdr_u1"]["raw_opcodes"]==4
 assert parsed["rev_b1a_apply_wrap32"]["ir"]<=4.0 and parsed["rev_b1a_apply_wrap4"]["ir"]<=4.0 and parsed["rev_b1_stream_wrap32"]["ir"]<=4.0 and parsed["rev_b1_stream_wrap4"]["ir"]<=4.0
 assert parsed["rev_b1_stream_u32"]["ir"]<=1.065 and parsed["rev_b1_stream_u4"]["ir"]<=1.50 and parsed["rev_b1_stream_tail1"]["ir"]<=4.0
 assert(parsed["rev_b1a_apply_u32"]["ir"]+parsed["rev_b1_stream_u32"]["ir"])<=coll["c_hit_lo_m4"]
 assert(parsed["rev_b1a_apply_u4"]["ir"]+parsed["rev_b1_stream_u4"]["ir"])<=coll["c_hit_lo_m4_s"]
 _orb_1bl=(parsed["rev_b1a_apply_u4"]["ir"]+2.0*parsed["rev_b1b_sub_u4"]["ir"]+parsed["rev_b1b_eor_u4"]["ir"])/4.0+parsed["rev_b1_stream_u4"]["ir"]
 _orb_1bs=(parsed["rev_b1a_apply_1"]["ir"]+2.0*parsed["rev_b1b_apply_1"]["ir"]+parsed["rev_b1b_eor_1"]["ir"])/4.0+parsed["rev_b1_stream_tail1"]["ir"]
 assert _orb_1bl<=float(coll.get("c_hit_1bl",6.75))and _orb_1bs<=float(coll.get("c_hit_1bs",11.25))
 assert parsed["rev_b2_emit_c0"]["raw_opcodes"]==10 and parsed["rev_b2_emit_c1"]["raw_opcodes"]==6 and parsed["rev_b2_flush_s0"]["raw_opcodes"]<=10 and parsed["rev_b2_flush_s1"]["raw_opcodes"]<=10
 assert(parsed["rev_b2_emit_c0"]["ir"]+parsed["rev_b1a_apply_u32"]["ir"]+16.0/8192.0)<=coll["c_hit_hi"]
 assert "bfi x4, x9, #58, #6" in parsed["rev_b2_emit_c0"]["insns"]and "b.hs .Lb2_c1" in parsed["rev_b2_emit_c0"]["insns"]and "str x4, [x10], #8" in parsed["rev_b2_emit_c0"]["insns"]and "str x4, [x11], #8" in parsed["rev_b2_emit_c1"]["insns"]
 _r1=random.Random(2026);_M64=(1<<64)-1;_BIT57=1<<57;_B58=(1<<58)-1;_max_1al=_max_1as=_max_1bl=-1e9
 for _klen in range(128):
  _oids=[_r1.randrange(0,64)for _ in range(_klen)];_acc1=[0]*64;_lpv=37
  for _o in _oids:_acc1[_o]+=_lpv
  def _sim_b1(u32_en,u4_op,u1_op):
   _si=0;_sd=[];_so=0;_xc=8*_klen
   if _xc:
    if u32_en:
     _so+=2;_cf=(_xc>=256);_xc=(_xc-256)&_M64
     if _cf:
      while True:
       assert _si+32<=_klen;_sd.extend(_oids[_si:_si+32]);_si+=32;_so+=34;_cf=(_xc>=256);_xc=(_xc-256)&_M64
       if not _cf:break
     _xc=(_xc+256)&_M64;_so+=2
    if _xc:
     _so+=2;_cf=(_xc>=32);_xc=(_xc-32)&_M64
     if _cf:
      while True:
       assert _si+4<=_klen;_sd.extend(_oids[_si:_si+4]);_si+=4;_so+=6;_cf=(_xc>=32);_xc=(_xc-32)&_M64
       if not _cf:break
     _xc=(_xc+32)&_M64;_so+=2
     if _xc:
      while True:
       assert _si<_klen;_sd.append(_oids[_si]);_si+=1;_so+=4;_xc=(_xc-8)&_M64
       if not _xc:break
   assert _si==_klen and _sd==_oids
   _ai=0;_ac=[0]*64;_ao=0;_xc=8*_klen
   if _xc:
    if u32_en:
     _ao+=2;_cf=(_xc>=256);_xc=(_xc-256)&_M64
     if _cf:
      while True:
       assert _ai+32<=_klen
       for _t in range(32):_ac[_sd[_ai+_t]]+=_lpv
       _ai+=32;_ao+=114;_cf=(_xc>=256);_xc=(_xc-256)&_M64
       if not _cf:break
     _xc=(_xc+256)&_M64;_ao+=2
    if _xc:
     _ao+=2;_cf=(_xc>=32);_xc=(_xc-32)&_M64
     if _cf:
      while True:
       assert _ai+4<=_klen
       for _t in range(4):_ac[_sd[_ai+_t]]+=_lpv
       _ai+=4;_ao+=u4_op;_cf=(_xc>=32);_xc=(_xc-32)&_M64
       if not _cf:break
     _xc=(_xc+32)&_M64;_ao+=2
     if _xc:
      while True:
       assert _ai<_klen;_ac[_sd[_ai]]+=_lpv;_ai+=1;_ao+=u1_op;_xc=(_xc-8)&_M64
       if not _xc:break
   assert _ai==_klen and _ac==_acc1
   return _so+_ao
  _max_1al=max(_max_1al,_sim_b1(True,16,6)-4.625*_klen)
  _max_1as=max(_max_1as,_sim_b1(False,16,6)-5.50*_klen)
  if _klen<64:
   _pk=[((_r1.randrange(4)<<56)|_o)for _o in _oids]
   for _kc,_r4,_u1 in ((0,16,6),(1,24,8),(2,20,7),(3,24,8)):
    _ag=[0]*256;_at=[0]*256
    for _w in _pk:_ag[(((_w>>56)-_kc)&3)*64+(_w&63)]+=_lpv
    _ai=0;_ao=0;_xc=8*_klen
    if _xc:
     _ao+=2;_cf=(_xc>=32);_xc=(_xc-32)&_M64
     if _cf:
      while True:
       for _t in range(4):
        _w=_pk[_ai+_t];_x4=_w if _kc==0 else(_w^_BIT57 if _kc==2 else((_w-(_kc<<56))&_M64)&_B58);_at[(_x4>>56)*64+(_x4&63)]+=_lpv
       _ai+=4;_ao+=_r4;_cf=(_xc>=32);_xc=(_xc-32)&_M64
       if not _cf:break
     _xc=(_xc+32)&_M64;_ao+=2
     if _xc:
      while True:
       _w=_pk[_ai];_x4=_w if _kc==0 else(_w^_BIT57 if _kc==2 else((_w-(_kc<<56))&_M64)&_B58);_at[(_x4>>56)*64+(_x4&63)]+=_lpv;_ai+=1;_ao+=_u1;_xc=(_xc-8)&_M64
       if not _xc:break
    assert _ag==_at
    _max_1bl=max(_max_1bl,(_ao-(_r4/4.0)*_klen)+11.5)
 assert 3.0+4.0*(4.0+_max_1al)<=float(coll.get("b1a_long_lk",248.0))and 3.0+4.0*(4.0+_max_1as)<=float(coll.get("b1a_short_lk",106.0))and 3.0+4.0+_max_1bl<=float(coll.get("b1bl_lk",29.0))
 assert parsed["ntt62_r4_shoup"]["ir"]*9.0<=55.13<55.25 and "subs x10, x0, x2" in parsed["ntt62_r4_shoup"]["insns"] and "csel x10, x9, x10, lo" in parsed["ntt62_r4_shoup"]["insns"] and "b.ne .Lntt_r4" in parsed["ntt62_r4_shoup"]["insns"]
 assert parsed["ntt62_butterfly"]["ir"]*12.0<=120.0<168.0
 assert coll["c_list_a"]>=32.0 and coll["c_lbkt"]>=10.0
 _p_b1=17;_va_b1=8;_st_b1=5;_pms_b1=_p_b1-_st_b1;_kb1=[];_h1=[0]*_p_b1;_x0=3
 for _ in range(4*_va_b1):
  _x0=(_x0-_pms_b1)if _x0>=_pms_b1 else(_x0+_st_b1);_kb1.append(_x0);_h1[_x0]+=1
 _o1=[0]
 for _c in _h1:_o1.append(_o1[-1]+_c)
 _pos1=list(_o1[:-1]);_out1=[None]*_o1[-1];_ix1=0
 for _mb in range(4):
  _xom=_mb<<56
  for _ in range(_va_b1):
   _ka=_kb1[_ix1];_ix1+=1;_out1[_pos1[_ka]]=_xom;_pos1[_ka]+=1;_xom+=1
 assert sorted((_w>>56,_w&((1<<56)-1))for _w in _out1)==[(m,o)for m in range(4)for o in range(_va_b1)]
 MASK=(1<<64)-1;MASK58=(1<<58)-1;OID_M=(1<<56)-1
 for m_b in range(4):
  for k_c in range(4):
   xk56=k_c<<56;krot56=((4-k_c)&3)<<56
   for oid in (0,1,0x12345,OID_M):
    om=(m_b<<56)|oid
    rot=((om-xk56)&MASK)&MASK58
    assert(rot>>56)==((m_b-k_c)&3)and(rot&OID_M)==oid and(rot<<1)<(1<<59)
    if m_b==0 and k_c==1:assert(rot>>56)==3
    for r_sb in (0,7,15):
     rec=(r_sb<<60)|(m_b<<56)|oid
     for t0_c in (0,1,7,15):
      xrot=(t0_c<<60)|krot56;s65=rec+xrot;c_fl=(s65>>64);x4=s65&MASK
      u_ex=x4>>60;x4_bfi=(x4&MASK58)|((u_ex&0x3F)<<58);x4=(x4&MASK58)|(u_ex<<58)
      assert x4_bfi==x4 and c_fl==(1 if r_sb+t0_c>=16 else 0)and(x4>>58)==((r_sb+t0_c)&15)and((x4>>56)&3)==((m_b-k_c)&3)and(x4&OID_M)==oid and(x4<<1)<(1<<63)
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
 def _sim_b2_k0_step4(kt4,xp,xd1,xb0):
  assert -xp<xb0<xp;xpms=(xp-xd1)&MASK;xb0u=xb0&MASK
  x0=((xb0u+xp)&MASK)if(xb0u>>63)else xb0u;x1=(xp-x0)&MASK;nkt=[];nk0=[]
  for kt in kt4:
   kn=((kt-xpms)&MASK)if kt>=xpms else((kt+xd1)&MASK);nkt.append(kn)
   nk0.append(((kn-x1)&MASK)if kn>=x1 else((kn+x0)&MASK))
  return nkt,nk0
 for _d1_64 in (0,1,p_hi64-1,(1<<63)+1234567):
  _kt64=[0,1,p_hi64-1,(1<<63)+999]
  for _b0_64 in (-448,-220,-1,0,1,220,448):
   _kt_e=[(_v+_d1_64)%p_hi64 for _v in _kt64];_k0_e=[(_v+_b0_64)%p_hi64 for _v in _kt_e]
   _kt64,_k0_g=_sim_b2_k0_step4(_kt64,p_hi64,_d1_64,_b0_64)
   assert _kt64==_kt_e and _k0_g==_k0_e and all(0<=_v<p_hi64 for _v in _kt64+_k0_g)
 for _pq in (31,61,151,211,241,1000081):
  _e15=(_pq-1)//15;_s0=next(pow(g,_e15,_pq)for g in range(2,80)if pow(pow(g,_e15,_pq),3,_pq)!=1 and pow(pow(g,_e15,_pq),5,_pq)!=1)
  _sk=[_s0,pow(_s0,7,_pq),pow(_s0,4,_pq),pow(_s0,13,_pq)]
  _pws=[(1,pow(s,5,_pq),pow(s,3,_pq),pow(s,8,_pq),pow(s,6,_pq),pow(s,11,_pq),pow(s,9,_pq),pow(s,14,_pq))for s in _sk]
  _w0=_pws[0][1];assert all(pw[1]==_w0 and pw[0]==1 for pw in _pws)
  _t27=[2,-1,3,0,-2,1];_b1l=-4;_ktq=[sum(([0,_b1l-1]+_t27)[j]*pw[j]for j in range(8))%_pq for pw in _pws]
  for _b1 in range(_b1l,_b1l+9):
   _b0s=((_b1*7)%min(25,_pq-1))-min(12,(_pq-1)//2)
   _ktq,_k0q=_sim_b2_k0_step4(_ktq,_pq,_w0,_b0s)
   assert _k0q==[sum(([_b0s,_b1]+_t27)[j]*pw[j]for j in range(8))%_pq for pw in _pws]
 def _sim_probe(xk0,xp,xlb,xcol0,krot56,xoff=0):
  if xlb==0:return 0,0,0,0,0,0
  x4=xk0>>4;x1=xcol0-xk0;x0=x1&15;xrot=krot56|(x0<<60);x22=((x4<<4)+x1)>>4
  x5=xp-xk0
  if x5>xlb:x5=xlb
  x6=(xk0+x5-1)>>4;x21=x6-x4+1;x19=xoff+(x4<<3)
  return x4,xrot,x5,x21,x19,x22
 _p61=2305843009213692601
 assert _sim_probe(15,_p61,0,0,3<<56,0x1000)[3]==0,"Empty chord L_b==0 must produce n_sb==0!"
 assert _sim_probe(15,_p61,189,0,3<<56,0x1000)==(0,(1<<60)|(3<<56),189,13,0x1000,-1)
 assert _sim_probe(0,_p61,184,21,3<<56,0x1000)==(0,(5<<60)|(3<<56),184,12,0x1000,1),"Wraparound x_col0=21 must yield slab 1, t0=5!"
 p_b2=509;L_b2=440;n_slabs2=28;krot_t=((4-1)&3)<<56
 kas_4b=[[0,1,15,16,31,32,439,447,508],[1,15,16,47,432,508],[2,17,33,64,446,505],[14,15,18,94,435,508]]
 kbuf_flat=[ka for b in kas_4b for ka in b]
 n_sb2=(p_b2+15)//16;hsb=[0]*n_sb2
 for ka in kbuf_flat:hsb[ka>>4]+=1
 off_b=[0]
 for c in hsb:off_b.append(off_b[-1]+c*8)
 fba=bytearray(off_b[-1]);psb=list(off_b[:-1]);idx_kb=0;r_prod=[]
 for m_b in range(4):
  x_oidm=m_b<<56
  for _ in range(len(kas_4b[m_b])):
   ka=kbuf_flat[idx_kb];idx_kb+=1;si=ka>>4
   rec=x_oidm|((ka&15)<<60);struct.pack_into("<Q",fba,psb[si],rec);psb[si]+=8
   r_prod.append((ka,m_b,x_oidm&OID_M));x_oidm+=1
 flat_buf=bytes(fba)
 _B2_SZ=65536;_B2_MASK=_B2_SZ-1
 for _k0_t in (0,1,7,15,16,31,68,69,440,p_b2-16,p_b2-5,p_b2-1):
  _bufs=[bytearray(_B2_SZ)for _ in range(30)];_tx=[0]*31;_slab_sg=[{}for _ in range(n_slabs2)]
  def _drain_u32(_S,_end,_wlp):
   assert 0<=_S<n_slabs2 and(_end&7)==0;_rx=0;_xc=_end
   while _xc>=256:
    for _ in range(32):_x4,=struct.unpack_from("<Q",_bufs[_S+1],_rx);_rx+=8;_slab_sg[_S][_x4]=_slab_sg[_S].get(_x4,0)+_wlp
    _xc-=256
   while _xc>=32:
    for _ in range(4):_x4,=struct.unpack_from("<Q",_bufs[_S+1],_rx);_rx+=8;_slab_sg[_S][_x4]=_slab_sg[_S].get(_x4,0)+_wlp
    _xc-=32
   while _xc>0:_x4,=struct.unpack_from("<Q",_bufs[_S+1],_rx);_rx+=8;_slab_sg[_S][_x4]=_slab_sg[_S].get(_x4,0)+_wlp;_xc-=8
   _tx[_S+1]=0
  _sbi0,_xrot0,_span,_nsb,_x19_0,_x22_0=_sim_probe(_k0_t,p_b2,L_b2,0,krot_t,0)
  assert _x19_0==_sbi0*8
  _segs=[(_sbi0,_nsb,_xrot0,_x22_0)]
  if _span<L_b2:
   _sbi1,_xrot1,_,_nsb1,_,_x22_1=_sim_probe(0,p_b2,L_b2-_span,p_b2-_k0_t,krot_t,0)
   _segs.append((_sbi1,_nsb1,_xrot1,_x22_1))
  for _sbi_s,_nsb_s,_xrot_s,_x22_s in _segs:
   assert -1<=_x22_s<=27;_x2=off_b[_sbi_s];_x14=_x22_s;_x10=_tx[_x14+1];_x11=_tx[_x14+2]
   for _st in range(_nsb_s):
    _x3=off_b[_sbi_s+_st+1]
    while _x2<_x3:
     _x0,=struct.unpack_from("<Q",flat_buf,_x2);_x2+=8
     _s65=_x0+_xrot_s;_cfl=_s65>>64;_x4=_s65&MASK;_x9=_x4>>60;_x4=(_x4&MASK58)|((_x9&0x3F)<<58)
     if _cfl==0:
      struct.pack_into("<Q",_bufs[_x14+1],_x10,_x4);_x10+=8
      if(_x10&_B2_MASK)==0:
       if _x14<0:_x10=0
       else:_drain_u32(_x14,_B2_SZ,232);_x10=0
     else:
      struct.pack_into("<Q",_bufs[_x14+2],_x11,_x4);_x11+=8
      if(_x11&_B2_MASK)==0:
       if _x14+1>=n_slabs2:_x11=0
       else:_drain_u32(_x14+1,_B2_SZ,232);_x11=0
    _tx[_x14+1]=_x10;_x14+=1;_x10=_x11;_x11=_tx[_x14+2]
   _tx[_x14+1]=_x10
  for _S in range(n_slabs2):
   if _tx[_S+1]>0:_drain_u32(_S,_tx[_S+1],232)
  em=[]
  for _S in range(n_slabs2):
   for _x4,_vl in _slab_sg[_S].items():
    assert _vl==232;_col=16*_S+(_x4>>58)
    if _col<L_b2:em.append((_col,(_x4>>56)&3,_x4&OID_M))
  exh=[((ka-_k0_t)%p_b2,(mb-1)&3,oo)for ka,mb,oo in r_prod if((ka-_k0_t)%p_b2)<L_b2]
  assert sorted(em)==sorted(exh),f"Band-2 rolling-cursor mismatch at k0={_k0_t}"
 _buf4=[bytearray(_B2_SZ)for _ in range(30)];_tx4=[0]*31;_sg4=[{}for _ in range(n_slabs2)];_exp_sg4=[{}for _ in range(n_slabs2)]
 _n_full_drains=0;_n_sink_drains=0;_n_cls_drains=0
 def _drain4(_S,_end,_wlp):
  for _off in range(0,_end,8):_x4,=struct.unpack_from("<Q",_buf4[_S+1],_off);_sg4[_S][_x4]=_sg4[_S].get(_x4,0)+_wlp
  _tx4[_S+1]=0
 _xrot4=(9<<60)|(3<<56)
 for _wlp in (232,236):
  for _x22c in (-1,0,27):
   _x14=_x22c;_x10=_tx4[_x14+1];_x11=_tx4[_x14+2]
   for _j in range(18000):
    _rec=((_j&15)<<60)|(((_j>>4)&3)<<56)|(_wlp*100000+(_x22c+2)*20000+_j)
    _s65=_rec+_xrot4;_cfl=_s65>>64;_x4=_s65&MASK;_x9=_x4>>60;_x4=(_x4&MASK58)|((_x9&0x3F)<<58)
    if _cfl==0:
     struct.pack_into("<Q",_buf4[_x14+1],_x10,_x4);_x10+=8
     if 0<=_x14<n_slabs2:_exp_sg4[_x14][_x4]=_exp_sg4[_x14].get(_x4,0)+_wlp
     if(_x10&_B2_MASK)==0:
      if _x14<0:_n_sink_drains+=1;_x10=0
      else:_drain4(_x14,_B2_SZ,_wlp);_n_full_drains+=1;_x10=0
    else:
     struct.pack_into("<Q",_buf4[_x14+2],_x11,_x4);_x11+=8
     if 0<=_x14+1<n_slabs2:_exp_sg4[_x14+1][_x4]=_exp_sg4[_x14+1].get(_x4,0)+_wlp
     if(_x11&_B2_MASK)==0:
      if _x14+1>=n_slabs2:_n_sink_drains+=1;_x11=0
      else:_drain4(_x14+1,_B2_SZ,_wlp);_n_full_drains+=1;_x11=0
   _tx4[_x14+1]=_x10;_tx4[_x14+2]=_x11
  for _S in range(n_slabs2):
   if _tx4[_S+1]>0:_drain4(_S,_tx4[_S+1],_wlp);_n_cls_drains+=1
  _tx4[0]=0;_tx4[29]=0
 assert _n_full_drains>=4 and _n_sink_drains>=2 and _n_cls_drains>=6 and _sg4==_exp_sg4,"Multi-w_lp-class 64 KiB rolling-cursor & sink Band-2 exactness!"
 _sc64=parsed["rev_scan_zero_u64"]
 assert _sc64["items"]==64 and _sc64["raw_opcodes"]==45 and _sc64["insns"][0]=="cbz x_cnt, .Lscan_done" and _sc64["insns"][-1]=="ret" and "stp x_b16, x_b16, [x_g], #128" in _sc64["insns"]and "stp xzr, xzr, [x_f], #128" in _sc64["insns"]
 assert _pw("rev_scan_zero_u64",_sc64["insns"][17:32]).isdisjoint({f"r{i}" for i in range(16)})
 _sga_eff=float(sga_mean)if sga_mean is not None else 5e-5
 _p_hit64=min(1.0,64.0*_sga_eff)
 scan_cold_ir=_p_hit64*(1000.0/64.0)
 scan_tot=_sc64["ir"]*(1.0+16.0/189.0)+scan_cold_ir
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
 _sg192=[_b16]*192;_sf192=[100]*192
 _sg192[13]=(_b16+_T)&0xFFFF;_sf192[13]=500       # block 0 lane 13 -> global oid 13 (ACCEPT)
 _sg192[64+13]=(_b16+_T)&0xFFFF;_sf192[64+13]=600 # block 1 lane 13 -> global oid 77 (ACCEPT)
 _sg192[64+47]=(_b16+_T)&0xFFFF;_sf192[64+47]=400 # block 1 lane 47 -> S_f < 472 (REJECT)
 _sg192[128+47]=(_b16+_T)&0xFFFF;_sf192[128+47]=472 # block 2 lane 47 -> global oid 175 (ACCEPT)
 _q_out=[]
 for _blk in range(3):
  _xf_off=_blk*128;_x2=_xf_off>>1
  _w16=[(_sg192[64*_blk+4*w]|(_sg192[64*_blk+4*w+1]<<16)|(_sg192[64*_blk+4*w+2]<<32)|(_sg192[64*_blk+4*w+3]<<48))for w in range(16)]
  for _wi in range(16):
   _x17=_w16[_wi]
   if _x17&_X80:
    for _ln in range(4):
     if _x17&0x8000:
      _x19=4*_wi+_ln
      if _sf192[64*_blk+_x19]>=472:_q_out.append(_x2+_x19)
     _x17>>=16
  for _jj in range(64):_sf192[64*_blk+_jj]=0
 assert _q_out==[13,77,175]and all(v==0 for v in _sf192)
 mont5_ops=5*parsed["mont5_cios_round"]["raw_opcodes"]+parsed["mont5_cios_epilogue"]["raw_opcodes"]
 assert mont5_ops<=MONT5 and "sub sp, sp, #80" in parsed["mont5_cios_epilogue"]["insns"]and "add sp, sp, #80" in parsed["mont5_cios_epilogue"]["insns"]and "stp x19, x20, [sp, #0]" in parsed["mont5_cios_epilogue"]["insns"]and "ldp x19, x20, [sp, #0]" in parsed["mont5_cios_epilogue"]["insns"]
 mont2_ops=parsed["mont2_cios"]["raw_opcodes"]
 assert mont2_ops<=48 and mont2_ops<=367.76*(2.0/5.0)**2 and "csel x0, x10, x7, hs" in parsed["mont2_cios"]["insns"]and "csel x1, x11, x8, hs" in parsed["mont2_cios"]["insns"]and parsed["mont2_cios"]["insns"][-1]=="ret"
 assert parsed["unit_run_4vec4_z12"]["ir"]<=7.65 and parsed["unit_run_4vec4_neg_z12"]["ir"]<=7.65
 assert parsed["unit_run_1vec4_z12"]["ir"]<=9.60 and parsed["unit_run_1vec4_neg_z12"]["ir"]<=9.60 and "add x7, x6, x7, lsr #2" in parsed["unit_run_1vec4_z12"]["insns"]
 assert parsed["gen_run1_1v_z12"]["ir"]<=24.0 and parsed["gen_run1_1v_neg_z12"]["ir"]<=24.0 and "add x8, x9, x8, lsr #2" in parsed["gen_run1_1v_z12"]["insns"]
 for _cj in (0,1,7,1023,(1<<55)-1):
  _off128=_cj*128;assert(_off128>>2)==_cj*32
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
 for _pv2 in [(1<<128)-1,(1<<128)-159,(1<<127)+1,(1<<64)+1,3]+[rng.randrange(1<<65,1<<128)|1 for _ in range(64)]:
  _n0,_n1=_pv2&MASK,_pv2>>64;_ni=(-pow(_n0,-1,1<<64))&MASK
  for _av2,_bv2 in (((_pv2-1,_pv2-1),(0,0),(1,_pv2-1))if _pv2 in((1<<128)-1,3)else((rng.randrange(_pv2),rng.randrange(_pv2)),)):
   _x0,_x1,_x2,_x3=_av2&MASK,_av2>>64,_bv2&MASK,_bv2>>64
   _x7,_x8=(_x0*_x2)&MASK,(_x0*_x2)>>64;_x9,_x10=(_x1*_x2)&MASK,(_x1*_x2)>>64
   _s=_x8+_x9;_x8,_c=_s&MASK,_s>>64;_x9=(_x10+_c)&MASK
   _m=(_x7*_ni)&MASK;_s=_x7+((_m*_n0)&MASK);assert(_s&MASK)==0;_c=_s>>64
   _s=_x8+((_m*_n1)&MASK)+_c;_x7,_c=_s&MASK,_s>>64;_s=_x9+_c;_x8,_x9=_s&MASK,_s>>64
   _s=_x7+((_m*_n0)>>64);_x7,_c=_s&MASK,_s>>64;_s=_x8+((_m*_n1)>>64)+_c;_x8,_c=_s&MASK,_s>>64;_x9=(_x9+_c)&MASK
   _s=_x7+((_x0*_x3)&MASK);_x7,_c=_s&MASK,_s>>64;_s=_x8+((_x1*_x3)&MASK)+_c;_x8,_c=_s&MASK,_s>>64;_s=_x9+_c;_x9,_x10=_s&MASK,_s>>64
   _s=_x8+((_x0*_x3)>>64);_x8,_c=_s&MASK,_s>>64;_s=_x9+((_x1*_x3)>>64)+_c;_x9,_c=_s&MASK,_s>>64;_x10=(_x10+_c)&MASK
   _m=(_x7*_ni)&MASK;_s=_x7+((_m*_n0)&MASK);assert(_s&MASK)==0;_c=_s>>64
   _s=_x8+((_m*_n1)&MASK)+_c;_x7,_c=_s&MASK,_s>>64;_s=_x9+_c;_x8,_c=_s&MASK,_s>>64;_x9=(_x10+_c)&MASK
   _s=_x7+((_m*_n0)>>64);_x7,_c=_s&MASK,_s>>64;_s=_x8+((_m*_n1)>>64)+_c;_x8,_c=_s&MASK,_s>>64;_x9=(_x9+_c)&MASK
   _t0,_c=(_x7+(_n0^MASK)+1)&MASK,(_x7+(_n0^MASK)+1)>>64;_t1,_c=(_x8+(_n1^MASK)+_c)&MASK,(_x8+(_n1^MASK)+_c)>>64;_c=(_x9+MASK+_c)>>64
   _r0,_r1=(_t0,_t1)if _c else(_x7,_x8)
   assert(_r0+(_r1<<64))==(_av2*_bv2*pow(1<<128,-1,_pv2))%_pv2
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
 assert parsed["unit_run_4vec4_wrap"]["ir"]<=4.0 and "subs x5, x5, #4" in parsed["unit_run_4vec4_wrap"]["insns"] and "adds x5, x5, #4" in parsed["unit_run_4vec4_wrap"]["insns"]
 for _uk in range(16):
  _ui=0;_ux=_uk;_uc=(_ux>=4);_ux=(_ux-4)&MASK
  if _uc:
   while True:
    assert _ui+4<=_uk,"ZERO unit_run_4vec4 overread at k=4/8/12!"
    _ui+=4;_uc=(_ux>=4);_ux=(_ux-4)&MASK
    if not _uc:break
  _ux=(_ux+4)&MASK
  assert _ui+_ux==_uk and 0<=_ux<4
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
 c_bc=float(desc_p["bound_cell"])if desc_p else 4200.0
 c_ps=float(desc_p["pair_setup"])if desc_p else float(1<<35)
 assert c_bc>=4200.0 and 16+4*(104+13*48+128+88)==3792.0<=c_bc-400.0
 _SQ_LUT=[math.ceil(math.sqrt(idx+1)*(2.0**9.5))for idx in range(2048)]
 assert max(_SQ_LUT)<65536
 _B4=[[4,1,0,0],[1,4,1,0],[0,1,4,1],[1,-1,1,4]]
 _G4=[[sum(_B4[i][t]*_B4[j][t]for t in range(4))for j in range(4)]for i in range(4)]
 _mu4=[[0.0]*4 for _ in range(4)];_d4=[0.0]*4
 for _i in range(4):
  for _j in range(_i):
   _mu4[_i][_j]=(_G4[_i][_j]-sum(_mu4[_i][_t]*_mu4[_j][_t]*_d4[_t]for _t in range(_j)))/_d4[_j]
  _d4[_i]=_G4[_i][_i]-sum(_mu4[_i][_t]**2*_d4[_t]for _t in range(_i))
 _Q32=1<<32
 _mu_q=[[int(round(_mu4[i][j]*_Q32))for j in range(4)]for i in range(4)]
 _d_q=[int(math.floor(_d4[i]*_Q32))for i in range(4)]
 _isd_q=[int(math.ceil((1.0/math.sqrt(_d4[i]))*(1<<32)))for i in range(4)]
 _R4=35;_R4sq=_R4*_R4
 _bf4=set()
 for _z3 in range(0,13):
  _v30,_v31,_v32,_v33=_z3*_B4[3][0],_z3*_B4[3][1],_z3*_B4[3][2],_z3*_B4[3][3]
  for _z2 in range(-12,13):
   if _z3==0 and _z2<0:continue
   _v20,_v21,_v22,_v23=_v30,_v31+_z2,_v32+4*_z2,_v33+_z2
   for _z1 in range(-12,13):
    if _z3==0 and _z2==0 and _z1<0:continue
    _v10,_v11,_v12,_v13=_v20+_z1,_v21+4*_z1,_v22+_z1,_v23
    _n3=_v12*_v12+_v13*_v13
    if _n3>_R4sq:continue
    for _z0 in range(-12,13):
     if _z3==0 and _z2==0 and _z1==0 and _z0<=0:continue
     _u0=_v10+4*_z0;_u1=_v11+_z0
     if _u0*_u0+_u1*_u1+_n3<=_R4sq:_bf4.add((_z0,_z1,_z2,_z3))
 _m4=4;_C_stk=[[0]*4 for _ in range(4)];_V_stk=[[0]*4 for _ in range(4)];_stk=[];_zcur=[0]*4;_enum4=set();_pruned_adv=0;_pruned_k0=0
 _x20=_m4-1;_x0=_R4sq*_Q32;_x8=0;_x13=0;_x10=0;_x7=_isd_q[_x20];_x17=_d_q[_x20]
 _state="ENTER"
 while _state!="DONE":
  if _state=="ENTER":
   _x1=_x0|1;_x2=(64-_x1.bit_length())^63;_x3=_x2&62;_x4=(_x1<<9)>>_x3
   _sr=_SQ_LUT[_x4]<<(_x3>>1)
   _rad=((_sr*_x7)>>30)+1+_x8
   _x9=(_x10-_rad+0xFFFFFFFF)>>32;_x12=(_x10+_rad)>>32
   if _x13==0 and _x9<0:_x9=0
   _state="SIB"
  elif _state=="SIB":
   if _x9>_x12:
    _state="POP";continue
   _diff=abs((_x9<<32)-_x10)
   _dx=0 if _diff<=_x8 else(_diff-_x8)
   _dx2=(_dx*_dx)>>32
   _x16=(_dx2*_x17)>>32
   if _x0<_x16:
    _pruned_adv+=1
    if _x20==0:_pruned_k0+=1
    _x9+=1;_state="SIB";continue
   _x19=_x0-_x16
   _zcur[_x20]=_x9
   if _x20==0:
    _v=[_V_stk[0][t]+_x9*_B4[0][t]for t in range(4)]
    _nz=sum(x*x for x in _v)
    if _nz!=0 and _nz<=_R4sq:_enum4.add(tuple(_zcur))
    _x9+=1;_state="SIB"
   else:
    _stk.append((_x9,_x12,_x0,_x13,_x8))
    _x8+=abs(_x9)
    if _x9>0:_x13=1
    _x22=_x20;_k_parent=_x20;_x20-=1
    for _t in range(4):_V_stk[_x20][_t]=_V_stk[_k_parent][_t]+_x9*_B4[_k_parent][_t]
    while _x22>0:
     _x22-=1;_C_stk[_x20][_x22]=_C_stk[_k_parent][_x22]-_x9*_mu_q[_k_parent][_x22]
    _x10=_C_stk[_x20][_x20];_x7=_isd_q[_x20];_x17=_d_q[_x20];_x0=_x19
    _state="ENTER"
  elif _state=="POP":
   _x20+=1
   if _x20==_m4:_state="DONE"
   else:
    _x9,_x12,_x0,_x13,_x8=_stk.pop()
    _x10=_C_stk[_x20][_x20];_x7=_isd_q[_x20];_x17=_d_q[_x20]
    _x9+=1;_state="SIB"
 assert _enum4==_bf4 and len(_bf4)>10000 and _pruned_adv>=2 and _pruned_k0>=1,f"DFS mismatch: enum={len(_enum4)}, bf={len(_bf4)}, pruned_adv={_pruned_adv}, pruned_k0={_pruned_k0}"
 def _canon_ambient_sign(vec):
  for _c in vec:
   if _c!=0:return tuple(-_x for _x in vec)if _c<0 else tuple(vec)
  return tuple(vec)
 _zw_neg=(-1,1,0,0)
 assert _zw_neg in _enum4
 _vw_neg=tuple(sum(_zw_neg[i]*_B4[i][t]for i in range(4))for t in range(4))
 assert _vw_neg==(-3,3,1,0) and _canon_ambient_sign(_vw_neg)==(3,-3,-1,0)
 _sub4=list(_enum4)[:1000]
 _canon4={_canon_ambient_sign((z[0]*4+z[1]+z[3],z[0]+z[1]*4+z[2]-z[3],z[1]+z[2]*4+z[3],z[2]+z[3]*4))for z in _sub4}
 assert len(_canon4)==len(_sub4) and all(next(c for c in v if c!=0)>0 for v in _canon4)
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
 _zero_v=(0,)*16;_x_nz_0=sum(c*c for c in _zero_v);assert _x_nz_0==0
 for _hv_raw in [_zero_v]+_hits_r+[tuple(-c for c in _h)for _h in _hits_r[:32]]:
  if sum(c*c for c in _hv_raw)==0:continue
  _hv=_canon_ambient_sign(_hv_raw)
  _x5=0;_x7=_R2r
  for _dd in range(16):
   _vd=_hv[_dd]
   _x5+=_Pd[_dd][_x7][_vd+_Rr]
   _x7-=_vd*_vd
  _half_off=(_dp_le[0][_R2r]+1)//2;_hrk=_x5-_half_off
  assert _x7>=0 and 0<=_hrk<17072,"Half-ball rank must lie in [0, C_R-1]!"
  _seen_rk.add(_hrk)
 assert len(_seen_rk)==552
 from decimal import Decimal as _Dc,localcontext as _lc,ROUND_CEILING as _RC
 with _lc()as _ctx:
  _ctx.prec=50;_ctx.rounding=_RC
  _ln2=_Dc(2).ln()
  _gap_ub=_Dc(8)*((_Dc(66)/_Dc(64)).ln()/_ln2+(_Dc(65)/_Dc(64)).ln()/_ln2+_Dc(1)/_Dc(65536))
  assert _gap_ub<=_Dc("0.5343")<_Dc("0.625")
 for _e30 in (0,1,2,3,10,1000):
  _ep1=_e30+1;_sq=65*_ep1
  assert(_sq-_ep1)>=64*_ep1 and(_sq+1+_e30)*64<=(_sq-_e30)*66
 for x0_t in [0,1,2,3,4,15,16,511,512,513,2047,2048,16383,16384,(1<<32)-1,1<<32,(1<<48)-1,(1<<53)-1,1<<54,(1<<55)-1]+[rng.randrange(1,1<<55)for _ in range(500)]:
  x1_t=x0_t|1;e_t=(63-(64-x1_t.bit_length()))&62;idx_t=(x1_t<<9)>>e_t
  assert 512<=idx_t<=2047 and (x1_t<<9)<(1<<64)
  sr_q14=_SQ_LUT[idx_t]<<(e_t>>1)
  assert sr_q14*sr_q14>=(x0_t<<28)
  if x0_t>=1024:assert sr_q14<=math.ceil(math.sqrt(x0_t)*(1<<14)*1.0020)+4
 LIM63=(1<<63)-1
 def _sim_bmul(x0,x1,x2,x3,x4,x5,x6,x7):
  assert all(abs(v)<(1<<60)for v in (x0,x1,x3,x4))and 0<=x2<(1<<60)and 0<=x5<(1<<60)
  x8=((x0*x3)>>64)-((x1*x4)>>64);x10=((x0*x4)>>64)+((x1*x3)>>64);x12=abs(x0)+abs(x1);x13=abs(x3)+abs(x4)
  x14=((x12*x5)>>64)+((x13*x2)>>64)+((x2*x5)>>64)+6
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
 _tw=(1<<32)-1;_wr,_wi,_we,_ws=_sim_bmul(_tw,-1,_tw,1,_tw,_tw,0,0)
 with _lc()as _ctx2:
  _ctx2.prec=50;_D=_Dc
  _p1r=_D(_tw)+_D(_tw)*_D(12)/_D(13);_p1i=_D(-1)+_D(_tw)*_D(5)/_D(13)
  _p2r=_D(1)-_D(_tw)*_D(5)/_D(13);_p2i=_D(_tw)+_D(_tw)*_D(12)/_D(13)
  _pr=(_p1r*_p2r-_p1i*_p2i);_pi=(_p1r*_p2i+_p1i*_p2r);_sc2=_D(1<<_ws)
  _dr=_pr-_D(_wr)*_sc2;_di=_pi-_D(_wi)*_sc2;_rr2=_D(_we)*_sc2
  assert _dr*_dr+_di*_di<=_rr2*_rr2 and (_D(5).sqrt()+_D(3))<=_D(6)
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
 s1_ds=math.floor((15.0*144.0)/math.log2(40.0/11.0))
 ds_ops=15*(40*5*3+16*12*4+15*12*10)+15*(4*12*3*5+2*40*8*3)+2000
 s2_ds=math.floor((340.0*math.log2(40.0)+7.5*144.0)/math.log2(100.0/99.0))
 ds_lll_bound=(2*s1_ds+16)*ds_ops+(2*s2_ds+16)*ds_ops
 assert ds_lll_bound<c_ps
 assert max(k*(16-k)/2.0*math.log2(40.0/11.0)+(2.0*k/16.0)*276.0+k*math.log2(15.0)for k in range(1,17))<=640.0
 s1_sq=math.floor((15.0*552.0)/math.log2(40.0/11.0))
 sq1_ops=15*(40*15*10+16*12*28+15*12*78)+15*(4*12*10*15+2*40*25*10)+10000
 s2_sq=math.floor((120.0*math.log2(15.0)+340.0*math.log2(40.0)+7.5*552.0)/math.log2(100.0/99.0))
 sq2_ops=sq1_ops+256*15*12
 sq_lll_bound=(2*s1_sq+16)*sq1_ops+(2*s2_sq+16)*sq2_ops
 assert sq_lll_bound<lll_node_cap
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
 return dict(asm_blocks_parsed=len(parsed),undef_asm_labels=len(undef_labels),per_kernel_reg_disjoint_ok=True,b2_buffered_64k_direct_ok=True,b2_full_and_rowend_drains_ok=True,b2_bank_rot_62b_ok=True,b1_pass12_orbit0_ok=True,dsieve16_dfs_4d_ok=True,dsieve16_ambient_sign_canon_ok=True,uint16_swar_ok=True,ntt62_garner_exact=True,ntt3_r4_r3_shoup_conv_ok=True,ntt62_r4_3x3p_stage_ir=round(9.0*parsed["ntt62_r4_shoup"]["ir"],2),f62_barrett_mac_64_exact=True,mont5_opcodes=mont5_ops,mont5_aapcs64_ok=True,mont5_sim_32_exact=True,mont2_opcodes=mont2_ops,mont2_sim_exact=True,row_red5_barrett_opcodes=parsed["row_red5_z12"]["raw_opcodes"],row_red5_barrett_47_exact=True,unit_run_4vec4_ir=round(parsed["unit_run_4vec4_z12"]["ir"],4),unit_run_1vec4_ir=round(parsed["unit_run_1vec4_z12"]["ir"],4),gen_run1_2v_ir=round(parsed["gen_run1_2v_z12"]["ir"],2),gen_run1_1v_ir=round(parsed["gen_run1_1v_z12"]["ir"],2),scan_zero_u16=round(scan_tot,4),dsieve16_rank_552_exact=True,sq_lll_uniform_log2=round(math.log2(sq_lll_bound),2),boot_lll_uniform_log2=round(math.log2(bt_lll_bound),2),all_raw_opcodes_ok=True)
def missing_recon_z12(log2U=14,seed=1,lf=70.79,lg=71.19,nf=546.0,ng=670.0,excess=1.09,p_f=58.0,p_g=63.5,ret_f=0.90,ret_g=0.82,sbf=375.0,sbg=315.0):
 import bisect as _bs,random as _rd,filtsim_z12 as _FS
 rg=_rd.Random(seed);Us=1<<log2U;uf,ug=nf/lf,ng/lg
 LiBf=(2.0**lf)/(lf*LN2);LiBg=(2.0**lg)/(lg*LN2);LiBs=LiBf+LiBg
 Upf=int(round(Us*LiBf/LiBs));Upg=Us-Upf;nr=int(round(excess*Us))
 rf0,rg0=_FS.rho_exact(uf),_FS.rho_exact(ug)
 kf=lambda t:_FS.rho_exact(uf-t)/rf0 if rf0>0 else 0.0
 kg=lambda t:_FS.rho_exact(ug-t)/rg0 if rg0>0 else 0.0
 nt=1000;ts=[0.80+0.20*(i+0.5)/nt for i in range(nt)];cf,cg=[],[];tf=tg=0.0
 for t in ts:tf+=(2.0**(lf*(t-1.0)))/t;tg+=(2.0**(lg*(t-1.0)))/t;cf.append(tf);cg.append(tg)
 cf=[x/tf for x in cf];cg=[x/tg for x in cg];cl_raw=[];cl=[];c_meta=[];ds=0.20/nt
 for side_id,cnt,cm,L,kp,ps,rm in ((0,Upf,cf,lf,kf,p_f,ret_f),(1,Upg,cg,lg,kg,p_g,ret_g)):
  for _ in range(cnt):
   ix=min(nt-1,_bs.bisect_left(cm,rg.random()));t=ts[ix]+(rg.random()-0.5)*ds
   lb=L*t;base_l=excess*LiBs*kp(t)/(2.0**lb)
   cl_raw.append(min(float(nr),base_l))
   cl.append(min(float(nr),base_l*(rm if lb>ps else 1.0)))
   c_meta.append((side_id,lb if lb>ps else 0.0))
 nc=len(cl);p_keep_est=0.618
 while True:
  nr_raw=max(nr,int(math.ceil(nr/p_keep_est)));sc_raw=nr_raw/nr
  row_bf=[0.0]*nr_raw;row_bg=[0.0]*nr_raw;uns_buf=[]
  for c in range(nc):
   sd_c,lb_c=c_meta[c]
   if lb_c<=0.0:continue
   k=min(nr_raw,_FS._poisson(rg,min(0.90*nr_raw,cl_raw[c]*sc_raw)))
   if not k:continue
   rb_arr=row_bf if sd_c==0 else row_bg
   if k==1:
    r=rg._randbelow(nr_raw);rb_arr[r]+=lb_c;uns_buf.append((r,c))
   else:
    for r in rg.sample(range(nr_raw),k):rb_arr[r]+=lb_c;uns_buf.append((r,c))
  surv=[r for r in range(nr_raw)if row_bf[r]<=sbf and row_bg[r]<=sbg]
  if len(surv)>=nr:
   raw_to_surv=[-1]*nr_raw
   for idx_s,r in enumerate(surv[:nr]):raw_to_surv[r]=idx_s
   assert all(row_bf[r]<=sbf+1e-9 and row_bg[r]<=sbg+1e-9 for r in surv[:nr]),"Stage-B residual cap violated in missing_recon_z12!"
   break
  p_keep_est*=0.5;assert p_keep_est>1e-6
 rc=[[]for _ in range(nr)];cr=[[]for _ in range(nc)];rd=[0]*nc
 for r,c in uns_buf:
  sr=raw_to_surv[r]
  if sr>=0:rc[sr].append(c);cr[c].append(sr);rd[c]+=1
 for c in range(nc):
  if c_meta[c][1]>0.0:continue
  k=min(nr,_FS._poisson(rg,cl_raw[c]))
  if not k:continue
  if k==1:
   sr=rg._randbelow(nr);rc[sr].append(c);cr[c].append(sr);rd[c]+=1
  else:
   for sr in rg.sample(range(nr),k):rc[sr].append(c);cr[c].append(sr);rd[c]+=1
 dg=list(rd);al=[True]*nr;st=[];qq=[c for c in range(nc)if dg[c]==1]
 while qq:
  c=qq.pop()
  if dg[c]!=1:continue
  r=next((rr for rr in cr[c]if al[rr]),-1)
  if r<0:continue
  al[r]=False;st.append((r,c))
  for c2 in rc[r]:
   dg[c2]-=1
   if dg[c2]==1:qq.append(c2)
 kn=[dg[c]>=2 for c in range(nc)]
 for r,c in reversed(st):
  if all(kn[c2]or c2==c for c2 in rc[r]):kn[c]=True
 out_b=[]
 for lo,hi in ((1.5,3.0),(3.0,5.0),(5.0,8.0),(8.0,14.0),(14.0,30.0)):
  cs=[c for c,lm in enumerate(cl)if lo<=lm<hi]
  if not cs:continue
  n_=len(cs);mf=sum(1 for c in cs if not kn[c])/n_
  md=sum(max(math.exp(-cl[c])*(1.0+cl[c]+0.5*cl[c]**2+(cl[c]**3)/6.0),math.exp(-0.75*cl[c])*(1.0+0.75*cl[c]))for c in cs)/n_
  l3=sum(1 for c in cs if(not kn[c])and rd[c]>=3)
  assert mf<=md,f"missing_recon band [{lo},{hi}): meas {mf} > mod {md}"
  out_b.append(dict(band=[lo,hi],n=n_,meas_miss=round(mf,5),mod_miss=round(md,5),model_over_meas=round(md/mf,2)if mf>0 else 999.0,lost_ge3=l3))
 return out_b
