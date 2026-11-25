"""cm_models.py — Self-contained."""
import array
import collections
import functools
import math
import operator as _OP
from decimal import Decimal,localcontext
LN2=math.log(2.0)
_SUMPROD=getattr(math,"sumprod",lambda a,b,_m=_OP.mul:sum(map(_m,a,b)))
GAMMA=0.5772156649015329
MONT5=280.0
ECM_TORSION=16.0
ECM_WHEELS=(30,210,2310,30030,510510)
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
class QuantizedRetention(Retention):
 def __init__(self,B_bits,X_bits,*,q_step=3.75,prime_min=2**33,step=0.25,max_bits=800,tilt=0.05,explicit_bits=14,power_terms=8,odd=1):
  super().__init__(B_bits,X_bits,prime_min=prime_min,step=step,max_bits=max_bits,tilt=tilt,explicit_bits=explicit_bits,power_terms=power_terms,odd=odd)
  self.q_step=float(q_step);qi=round(q_step/step)
  assert abs(qi*step-q_step)<1e-12,"q_step must be an exact multiple of step!"
  self.q_i=qi;self._dmax=72;self._qc={}
  hi=min(float(X_bits),float(B_bits));lo=max(float(explicit_bits),self._lmin)
  self._wmax=int(math.ceil(max_bits/q_step*1.12))+8
  ag={};asf={};D=self._dmax;W=self._wmax
  if hi>lo:
   w0=step/4.0;cells=math.ceil((hi-lo)/w0);w=(hi-lo)/cells;_exp=math.exp
   for i in range(cells):
    b=lo+(i+0.5)*w;bi=round(b/step);wq=round(b/q_step);di=bi-wq*qi
    for k in range(1,3):
     wk=k*wq;dk=k*di
     if wk>W or dk<-D or dk>D:break
     ex=k*tilt*bi*step-(k-1)*LN2*b
     if ex<-75.0:break
     amt=_exp(ex)*w/(k*b);key=(wk,dk)
     ag[key]=ag.get(key,0.0)+amt
     asf[key]=asf.get(key,0.0)+(-amt if(k&1)==0 else amt)
  self._ag=ag;self._asf=asf;self._C_geo=None;self._C_sf=None
 def _get_C(self,sf):
  if sf and self._C_sf is not None:return self._C_sf
  if (not sf) and self._C_geo is not None:return self._C_geo
  adict=self._asf if sf else self._ag;norm=self.normalizers["sieve_sf" if sf else "sieve_geo"]
  D=self._dmax;qi=self.q_i;ln=self.length;W=min(self._wmax,(ln+D)//qi);sp=_SUMPROD
  by_w=collections.defaultdict(list)
  for (wk,dk),av in adict.items():
   if av!=0.0:by_w[wk].append((dk,wk*av))
  items=[];pad=0
  for wk,lst in sorted(by_w.items()):
   dlo=min(dk for dk,_ in lst);dhi=max(dk for dk,_ in lst)
   if -dlo>pad:pad=-dlo
   if dhi>pad:pad=dhi
   ker=[0.0]*(dhi-dlo+1)
   for dk,wa in lst:ker[dk-dlo]+=wa
   items.append((wk,dlo,dhi,tuple(ker[::-1])))
  span=2*D+1;zpad=[0.0]*pad
  C=[[0.0]*span for _ in range(W+1)];C[0][D]=math.exp(-norm)
  Cp=[None]*(W+1);Cp[0]=zpad+C[0]+zpad
  lo_a=[span]*(W+1);hi_a=[-1]*(W+1);lo_a[0]=D;hi_a[0]=D
  for wv in range(1,W+1):
   row=C[wv];inv_w=1.0/wv;rlo=span;rhi=-1
   for wk,dlo,dhi,krev in items:
    if wk>wv:break
    pw=wv-wk;plo=lo_a[pw];phi=hi_a[pw]
    if plo>phi:continue
    prev=Cp[pw];d0=plo+dlo;d1=phi+dhi
    if d0<0:d0=0
    if d1>=span:d1=span-1
    if d0<=d1:
     if d0<rlo:rlo=d0
     if d1>rhi:rhi=d1
     o0=pad-dhi;o1=pad-dlo+1
     for d in range(d0,d1+1):
      row[d]+=sp(krev,prev[d+o0:d+o1])
   if rlo<=rhi:
    for d in range(rlo,rhi+1):row[d]*=inv_w
    lo_a[wv]=rlo;hi_a[wv]=rhi;Cp[wv]=zpad+row+zpad
  sv=self.co("sieve_sf" if sf else "sieve_geo");mg=[0.0]*(ln+1)
  for wv in range(W+1):
   rw=C[wv];base=wv*qi-D;rlo=lo_a[wv];rhi=hi_a[wv]
   if rlo<=rhi:
    for j in range(rlo,rhi+1):
     n=base+j
     if 0<=n<=ln:mg[n]+=rw[j]
  assert max(abs(mg[n]-sv[n])for n in range(ln+1))<1e-9,"QuantizedRetention 2D marginal mismatch!"
  if sf:self._C_sf=C
  else:self._C_geo=C
  return C
 def _spass_rev(self,T_u8,sf):
  key=("sp",T_u8,sf);v=self._qc.get(key)
  if v is not None:return v
  C=self._get_C(sf);ln=self.length;qi=self.q_i;D=self._dmax
  arr=[0.0]*(ln+1)
  for wv in range(max(0,int(T_u8)),len(C)):
   row=C[wv];base=wv*qi-D
   j0=0 if base>=0 else -base;j1=2*D if base+2*D<=ln else ln-base
   for j in range(j0,j1+1):arr[base+j]+=row[j]
  s_lo=next((i for i,x in enumerate(arr)if x!=0.0),ln+1)
  s_hi=ln-next((i for i,x in enumerate(reversed(arr))if x!=0.0),ln+1)
  v=(arr[::-1],s_lo,s_hi);self._qc[key]=v;return v
 def _m_conv(self,lim_b,win,sf):
  mk=(lim_b,win,sf);c=self._c;m=c.get(mk)
  if m is None:
   D=self.length if self.length<2048 else 2048
   self.two(D,D,lim_b,win,sf);m=c[mk]
  bk=("mb",lim_b,win,sf);b=self._qc.get(bk)
  if b is None:
   ml=len(m)-1
   m_lo=next((i for i,x in enumerate(m)if x!=0.0),ml+1)
   m_hi=ml-next((i for i,x in enumerate(reversed(m))if x!=0.0),ml+1)
   b=(m,m_lo,m_hi);self._qc[bk]=b
  return b
 def at_u8(self,nb,T_u8,sb=None):
  st=self.step;n=round(nb/st)
  if n<0 or n>self.length:return 0.0
  jmin=int(T_u8)*self.q_i-self._dmax
  if jmin<0:jmin=0
  u_lim=n-jmin
  if u_lim<0:return 0.0
  rb=(math.floor(math.ceil(nb)/self._lmin)+1)*0.0625
  D=self.length if self.length<2048 else 2048
  b=9**9 if sb is None else math.floor((sb-rb-0.125)/st)
  if b<0:return 0.0
  lim=n if n<D else D
  if u_lim<lim:lim=u_lim
  lim_b=D if b>=lim else b
  ck=("atu8",n,T_u8,lim_b);v=self._qc.get(ck)
  if v is not None:return v
  m,m_lo,m_hi=self._m_conv(lim_b,None,True);sv_r,s_lo,s_hi=self._spass_rev(T_u8,True);off=self.length-n
  lo_k=m_lo if m_lo>n-s_hi else n-s_hi
  hi_k=lim if lim<m_hi else m_hi
  if n-s_lo<hi_k:hi_k=n-s_lo
  if lo_k>hi_k:v=0.0
  else:
   r=self.sel_norm*_SUMPROD(m[lo_k:hi_k+1],sv_r[off+lo_k:off+hi_k+1])/self.denominator[n]
   v=0.0 if r<0.0 else(1.0 if r>1.0 else r)
  self._qc[ck]=v;return v
 def survivor_upper_u8(self,nb,T_u8,sb=None,lb=None,kmax=4):
  st=self.step;n=round(nb/st)
  if n<0 or n>self.length:return 0.0
  jmin=int(T_u8)*self.q_i-self._dmax
  if jmin<0:jmin=0
  u_lim=n-jmin
  if u_lim<0:return 0.0
  D=self.length if self.length<2048 else 2048;lim=n if n<D else D
  if u_lim<lim:lim=u_lim
  if sb is None:
   lim_b=D;win=None
  else:
   lim_b=math.floor((sb+0.25)/st)
   if lim_b<0:return 0.0
   lim_b=lim_b if lim_b<D else D
   win=(kmax,lb,(self.cap(nb)+2)/16.0+0.25)
  ck=("su8",n,T_u8,lim_b,win);v=self._qc.get(ck)
  if v is not None:return v
  m,m_lo,m_hi=self._m_conv(lim_b,win,False)
  sv_r,s_lo,s_hi=self._spass_rev(T_u8,False);off=self.length-n
  lo_k=m_lo if m_lo>n-s_hi else n-s_hi
  hi_k=lim if lim<m_hi else m_hi
  if n-s_lo<hi_k:hi_k=n-s_lo
  if lo_k>hi_k:v=0.0
  else:
   r=_SUMPROD(m[lo_k:hi_k+1],sv_r[off+lo_k:off+hi_k+1])/self.denominator[n]
   v=0.0 if r<0.0 else(1.0 if r>1.0 else r)
  self._qc[ck]=v;return v
@functools.lru_cache(None)
def ecm_instr(fbits,nbits,success=0.999999,lo_r=48.0):
 import q15_ecm16 as _EC
 return _EC.ladder_instr(lo_r,float(fbits),int(math.ceil(nbits)),success=success,step=12.0,b_min=96,l_min=2)
def ecm_mean_var(fbits,nbits,success=0.999999,lo_r=48.0):
 import q15_ecm16 as _EC
 cost=_EC.ladder_instr(lo_r,float(fbits),int(math.ceil(nbits)),success=success,step=12.0,b_min=96,l_min=2)
 return (cost,0.0)
NTT4_PRIMES=(269*(1<<53)+1,276*(1<<53)+1,305*(1<<53)+1,309*(1<<53)+1,)
NTT4_PROD=NTT4_PRIMES[0]*NTT4_PRIMES[1]*NTT4_PRIMES[2]*NTT4_PRIMES[3]
NTT3_PRIMES=(456*(1<<52)+1,447*(1<<52)+1,432*(1<<52)+1)
NTT3_PROD=NTT3_PRIMES[0]*NTT3_PRIMES[1]*NTT3_PRIMES[2]
def ntt3_F(M_prod):
 m=math.ceil(max(M_prod,1.0))
 if m<=1:return 12.0
 sb=(2.0*m*m+4.0*m)if m<=256 else 1e30
 if m>(1<<52):return 3.0*ntt3_F(math.ceil(m/2.0))+20.0*m
 k=max(2,(m-1).bit_length());n2=float(1<<k)
 c2=37.50*n2*k+(26.0*n2 if(k&1)else 0.0)+100.0*n2+(512.0*n2 if n2>(1<<26)else 0.0)
 n3=3.0*float(1<<(k-2))
 if n3>=m and k>=3:
  c3=37.50*n3*(k-2)+(26.0*n3 if((k-2)&1)else 0.0)+244.0*n3+(512.0*n3 if n3>(1<<26)else 0.0)
  c2=min(c2,c3)
 return min(sb,c2)
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
def verify_kernels_symbolic(kernels_text,coll,la,r_bls,lll_node_cap,lll_boot_cap,desc_p=None,sga_mean=None,lb_live=156.13):
 import re,random,struct
 allowed={"ldp","stp","ldr","str","ldrb","strb","ldrh","strh","staddb","stadd","ldur","stur","add","adds","adc","adcs","sub","subs","sbc","sbcs","mul","umulh","smulh","madd","msub","and","orr","eor","bic","lsl","lsr","asr","bfi","ubfiz","ubfx","cmp","cmn","tst","csel","cset","csinc","mov","b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","b.ge","b.lt","b.le","clz","udiv","cbz","cbnz","tbz","tbnz","ret"}
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
   elif op in ("tbz","tbnz"):all_refs.add(ps[3])
   if op in ("ldp","stp"):
    mb=re.search(r"\[([^\]]+)\]",ins)
    assert mb and re.match(r"^\s*(?:x\w+|sp)(\s*,\s*#-?\d+)?\s*$",mb.group(1)),f"Illegal {op} addr mode in {name}: {ins}"
   if op=="staddb":
    assert re.match(r"^staddb\s+w\w+\s*,\s*\[x\w+\]$",ins.lower()),f"Illegal {op} form in {name}: {ins}"
   if op=="stadd":
    assert re.match(r"^stadd\s+x\w+\s*,\s*\[x\w+\]$",ins.lower()),f"Illegal {op} form in {name}: {ins}"
  ir=len(insns)/items
  assert ir<=charged+1e-12,f"{name}: {ir}>{charged}"
  parsed[name]=dict(raw_opcodes=len(insns),items=items,ir=ir,charged=charged,insns=insns)
 undef_labels=all_refs-all_defs
 assert not undef_labels,f"Undefined ASM labels in KERNELS.md: {sorted(undef_labels)}"
 _KERN_REG_MAPS={
 "rev_b2_row_k0_step4":{"x_st":"r19","x_p":"r20","x_d1":"r21","x_pms":"r22","x_b0":"r23"},
 "rev_b1_kt_step_shared":{"x_pms":"r21","x_st":"r22"},
 "rev_b1_bkt_hdr_cnt":{"x_off":"r23","x_p":"r19","x_end":"r24","x_cnt":"r25"},
 "rev_b1_stream_pre64":{"x_cnt":"r25"},
 "rev_b1_stream_u64":{"x_src":"r19","x_dst":"r20","x_cnt":"r25"},
 "rev_b1_stream_bin32":{"x_src":"r19","x_dst":"r20","x_cnt":"r25"},
 "rev_b1a_apply_pre64":{"x_cnt":"r25"},
 "rev_b1a_apply_f_u64":{"x_p":"r19","x_cnt":"r25","w_lp":"r26"},
 "rev_b1a_apply_f_bin32":{"x_p":"r19","x_cnt":"r25","w_lp":"r26"},
 "rev_w1_probe_fast":{"x_lb":"r12","x_p":"r9","x_k0":"r8","x_col0":"r14","x_t0":"r23","x_crot":"r27","x_tx":"r20","x_tx_base":"r22","x_oc":"r24","x_off":"r19","x_oe":"r25"},
 "rev_b2_w1_wrap":{"x_col0":"r14","x_oc":"r24","x_oe":"r25","x_on":"r26"},
 "rev_b2_w1_slab_hdr":{"x_oc":"r24","x_cnt":"r13","x_tx":"r20"},
 "rev_b2_w1_tail_slab":{"x_oc":"r24","x_oe":"r25","x_cnt":"r13","x_tx":"r20"},
 "rev_b2_w1_prefix_emit":{"x_oc":"r24","x_on":"r26","x_cnt":"r13","x_tx":"r20"},
 "rev_b2_flush_x10":{"x_rx":"r15","x_fcnt":"r14"},
 "rev_b2_flush_w8_x10":{"x_rx":"r15","x_fcnt":"r14"},
 "rev_b2_flush_x11":{"x_rx":"r15","x_fcnt":"r14"},
 "rev_b2_rx_drain_pre64":{"x_fcnt":"r14"},
 "rev_b2_rx_drain_u64":{"x_rx":"r15","x_fcnt":"r14","w_lp":"r28"},
 "rev_b2_rx_drain_bin32":{"x_rx":"r15","x_fcnt":"r14","w_lp":"r28"},
 "rev_b2_emit_id_pre16":{"x_cnt":"r13"},
 "rev_b2_emit_id_u16":{"x_t0":"r23","x_cnt":"r13"},
 "rev_b2_emit_id_bin8421":{"x_t0":"r23","x_cnt":"r13"},
 "rev_b2_emit_c123_pre16":{"x_cnt":"r13"},
 "rev_b2_emit_c123_u16":{"x_t0":"r23","x_cnt":"r13"},
 "rev_b2_emit_c123_bin8421":{"x_t0":"r23","x_cnt":"r13"},
 "rev_w8_probe_fast":{"x_lb":"r12","x_p":"r9","x_k0":"r26","x_col0":"r11","x_t0":"r23","x_crot":"r27","x_off_base":"r24","x_off":"r19","x_oe":"r25"},
 "rev_b2_w8_wrap":{"x_col0":"r11","x_tx":"r20","x_tx_base":"r22","x_off":"r19","x_oe":"r25"},
 "rev_b2_w8_slab_hdr":{"x_off":"r19","x_tx":"r20","x_cnt":"r13"},
 "rev_b2_w8_tail_slab":{"x_off":"r19","x_oe":"r25","x_tx":"r20","x_cnt":"r13"},
 "rev_b2_w8_prefix_pure":{"x_col0":"r11","x_k0":"r26","x_off":"r19","x_oe":"r25","x_on":"r12","x_tx":"r20","x_tx_base":"r22","x_t0":"r23","x_cnt":"r13"},
 "rev_b2_emit_split_id_pre4":{"x_cnt":"r13"},
 "rev_b2_emit_split_id_u4":{"x_t0":"r23","x_cnt":"r13"},
 "rev_b2_emit_split_id_bin21":{"x_t0":"r23","x_cnt":"r13"},
 "rev_b2_emit_split_c123_pre4":{"x_cnt":"r13"},
 "rev_b2_emit_split_c123_u4":{"x_t0":"r23","x_cnt":"r13"},
 "rev_b2_emit_split_c123_bin21":{"x_t0":"r23","x_cnt":"r13"},
 "rev_sb_probe_fast":{"x_lb":"r12","x_k0":"r8","x_col0":"r14","x_tx":"r20","x_tx_base":"r22","x_t0":"r23","x_crot":"r27","x_p":"r9","x_slot":"r19","x_slot_base":"r25"},
 "rev_b2_sb_wrap":{},
 "rev_b2_sieve_direct1_u2":{"x_slot":"r19","x_tx":"r20"},
 "rev_b2_sieve_direct1_hit2":{"x_t0":"r23","x_tx":"r20","x_one":"r22","x_ovf_adj":"r24","x_rx":"r15","x_fcnt":"r14"},
 "rev_b2_sieve_direct1_tail1":{"x_slot":"r19","x_t0":"r23","x_tx":"r20","x_ovf_adj":"r24","x_rx":"r15","x_fcnt":"r14"},
 "rev_scan_zero_wrap":{"x_cnt":"r25"},
 "rev_scan_g_u8_fused_128":{"x_b1":"r21","x_b2":"r22","x_mg":"r23","x_b8":"r24","x_8080":"r25","x_cnt":"r26"},
 "rev_scan_f_u8_128":{"x_b1":"r21","x_b2":"r22","x_mg":"r23","x_init_f":"r24","x_cnt":"r26"},
 "rev_pad_zero_u8_128":{"x_b2":"r22","x_cnt":"r26"},
 "rev_b2_build_subchord_hdr":{"x_sc":"r19","x_tab0":"r20","x_k":"r21","x_kbase":"r22","x_p":"r23","x_cnt":"r24","x_d":"r25","x_d4":"r26","x_d2":"r27","x_pms":"r28","x_pms4":"r16","x_pms2":"r17"},
 "rev_b2_build_fiber_hdr":{"x_tab1":"r20","x_kbase":"r22","x_klvl2":"r15","x_p":"r23"},
 "rev_b2_build_hist_w1_stadd":{"x_k":"r19","x_d":"r20","x_pms":"r21","x_kbuf":"r22","x_hist":"r23","x_eight":"r24","x_cnt":"r25"},
 "rev_b2_build_scatter_w1_u8":{"x_kbuf":"r19","x_pos":"r20","x_om":"r21","x_cnt":"r22"},
 "rev_b2_build_hist_u8_stadd":{"x_k":"r19","x_d":"r20","x_pms":"r21","x_kbuf":"r22","x_hist":"r23","x_eight":"r24","x_cnt":"r25"},
 "rev_build_prefix_u64_8x":{"x_pos":"r20","x_cur":"r21","x_cnt":"r22"},
 "rev_b2_build_scatter_r1_u8":{"x_kbuf":"r19","x_pos":"r20","x_om":"r21","x_cnt":"r22"},
 "rev_b2_memset_ff_u16":{"x_slot":"r19","x_cnt":"r20","x_ff":"r21"},
 "rev_b2_build_direct1_u8":{"x_k":"r19","x_d":"r20","x_pms":"r21","x_slot":"r22","x_om":"r23","x_cnt":"r24"},
 "rev_b2_build_direct1_ovf":{"x_slot":"r22","x_ovf":"r25","x_ovf_tag":"r26"},
 "mont5_cios_round":{"x_b":"r6","x_n0":"r7"},
 "mont5_cios_epilogue":{"x_a":"r5","x_out":"r8","x_p":"r9"},
 "ntt62_butterfly":{"x_u":"r19","x_v":"r20","x_w":"r21","x_p":"r22","x_cnt":"r23"},
 "ntt61_r4_u2":{"x_a0":"r19","x_a1":"r20","x_a2":"r21","x_a3":"r22","x_w1":"r23","x_s1":"r24","x_w2":"r25","x_s2":"r26","x_w3":"r27","x_s3":"r28","x_p":"r4","x_j":"r5","x_jp":"r9","x_cnt":"r13"},
 "lehmer32_step":{},
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
   if op in ("cmp","cmn","tst","b","b.ne","b.eq","b.lo","b.ls","b.hi","b.hs","b.gt","b.ge","b.lt","b.le","cbz","cbnz","tbz","tbnz","ret","str","strb","strh","staddb","stadd","stur","stp"):continue
   for r in [ps[1]]+(ps[2:3]if op=="ldp" else[]):
    nr=_norm_r(r)
    if nr:ws.add(nr)
  return ws
 assert _pw("rev_b2_row_k0_step4").isdisjoint({"r19","r20","r21","r22","r23"})
 assert _pw("rev_b2_flush_x10").isdisjoint({"r13"}) and _pw("rev_b2_flush_w8_x10").isdisjoint({"r13"}) and _pw("rev_b2_flush_x11").isdisjoint({"r13"})
 _drain_w=_pw("rev_b2_rx_drain_pre64")|_pw("rev_b2_rx_drain_u64")|_pw("rev_b2_rx_drain_bin32")
 _b2_live={"r1","r2","r3","r8","r10","r11","r12","r13","r19","r20","r21","r22","r23","r24","r25","r26","r27","r28"}
 assert _drain_w.isdisjoint(_b2_live),f"rev_b2_rx_drain clobbers live Band-2 regs: {_drain_w&_b2_live}"
 _emit_w=_pw("rev_b2_emit_id_u16")|_pw("rev_b2_emit_id_bin8421")|_pw("rev_b2_emit_c123_u16")|_pw("rev_b2_emit_c123_bin8421")
 _w8_pre_live={"r11","r19","r21","r22","r23","r25","r26","r28"}
 assert (_drain_w|_emit_w).isdisjoint(_w8_pre_live),f"W8 prefix emit/drain clobbers live suffix regs: {(_drain_w|_emit_w)&_w8_pre_live}"
 assert "lsr x1, x_t0, #59" in parsed["rev_b2_w8_prefix_pure"]["insns"]
 assert parsed["rev_b2_rx_drain_pre64"]["raw_opcodes"]==2 and parsed["rev_b2_rx_drain_u64"]["raw_opcodes"]==98 and parsed["rev_b2_rx_drain_bin32"]["raw_opcodes"]==103
 assert parsed["rev_b2_rx_drain_u64"]["ir"]+(4.0+2.0+103.0)/16384.0<=1.54
 assert parsed["rev_b2_row_k0_step4"]["raw_opcodes"]==33 and parsed["rev_b2_row_k0_step4"]["items"]==4 and (parsed["rev_b2_row_k0_step4"]["ir"]+64.0/lb_live)<=float(coll["c_b2_row_k0"])
 assert parsed["rev_w1_probe_fast"]["raw_opcodes"]==12 and parsed["rev_b2_w1_wrap"]["raw_opcodes"]==11 and (12.0+5.0+2.0)<=coll["c_tile_hdr_w1"]
 assert "lsr x21, x7, #8" in parsed["rev_b2_w1_wrap"]["insns"] and "cbz x21, .Lw1_tail_slab" in parsed["rev_b2_w1_wrap"]["insns"]
 assert "cmp x_oc, x_oe" in parsed["rev_b2_w1_tail_slab"]["insns"] and "b.hs .Lw1_done" in parsed["rev_b2_w1_tail_slab"]["insns"]
 assert parsed["rev_sb_probe_fast"]["raw_opcodes"]==13 and parsed["rev_b2_sb_wrap"]["raw_opcodes"]==5 and (parsed["rev_sb_probe_fast"]["ir"]+parsed["rev_b2_sb_wrap"]["ir"])<=coll["c_tile_hdr"]
 assert "subs x21, x21, #1" in parsed["rev_b2_sb_wrap"]["insns"] and "b.lo .Lsd1_tail1" in parsed["rev_b2_sb_wrap"]["insns"] and ".Lsd1_tail1" in all_defs
 assert "asr x1, x0, #5" in parsed["rev_sb_probe_fast"]["insns"] and "add x_tx, x_tx_base, x1, lsl #3" in parsed["rev_sb_probe_fast"]["insns"]
 assert parsed["rev_b2_w1_slab_hdr"]["raw_opcodes"]==11 and (parsed["rev_b2_w1_slab_hdr"]["ir"]+8.875)<=coll["c_sb_hdr_u4"]
 assert parsed["rev_w8_probe_fast"]["raw_opcodes"]==13 and parsed["rev_b2_w8_wrap"]["raw_opcodes"]==6 and (13.0+6.0)<=coll["c_tile_hdr_w8"]
 assert "lsr x21, x7, #5" in parsed["rev_b2_w8_wrap"]["insns"] and "cbz x21, .Lw8_tail_slab" in parsed["rev_b2_w8_wrap"]["insns"]
 assert "cmp x_off, x_oe" in parsed["rev_b2_w8_tail_slab"]["insns"] and "ldr x8, [x_oe]" in parsed["rev_b2_w8_tail_slab"]["insns"]
 assert parsed["rev_b2_w8_slab_hdr"]["raw_opcodes"]<=16 and (parsed["rev_b2_w8_slab_hdr"]["ir"]+8.875+6.0)<=coll["c_sb_hdr_u2"]
 assert "ldp x11, x10, [x_tx]" in parsed["rev_b2_w8_slab_hdr"]["insns"] and "stp x11, x10, [x_tx]" in parsed["rev_b2_w8_slab_hdr"]["insns"] and _emit_w.isdisjoint({"r11"})
 assert "ldp x11, x10, [x_tx]" in parsed["rev_b2_w8_tail_slab"]["insns"] and "stp x11, x10, [x_tx]" in parsed["rev_b2_w8_tail_slab"]["insns"]
 assert "tbnz x10, #17, .Lw8_flush_x10" in parsed["rev_b2_w8_slab_hdr"]["insns"] and "b .Lw8_chk_x11" in parsed["rev_b2_flush_w8_x10"]["insns"]
 assert parsed["rev_b1_kt_step_shared"]["raw_opcodes"]==3 and parsed["rev_b1_bkt_hdr_cnt"]["raw_opcodes"]==4
 assert parsed["rev_b1_stream_pre64"]["raw_opcodes"]==2 and parsed["rev_b1_stream_u64"]["ir"]<=1.035 and parsed["rev_b1_stream_bin32"]["raw_opcodes"]==72
 assert parsed["rev_b1a_apply_pre64"]["raw_opcodes"]==2 and parsed["rev_b1a_apply_f_u64"]["ir"]<=1.535 and parsed["rev_b1a_apply_f_bin32"]["raw_opcodes"]==103
 assert (parsed["rev_b1a_apply_f_u64"]["ir"]+parsed["rev_b1_stream_u64"]["ir"])<=coll["c_hit_lo_m4_f"] and (parsed["rev_b1a_apply_f_u64"]["ir"]+parsed["rev_b1_stream_u64"]["ir"])<=coll["c_hit_lo_m4_g"]
 _w_u16=(parsed["rev_b2_emit_id_u16"]["ir"]+3.0*parsed["rev_b2_emit_c123_u16"]["ir"])/4.0
 _c_w1=_w_u16+parsed["rev_b2_rx_drain_u64"]["ir"]+0.015
 _c_spl=((parsed["rev_b2_emit_split_id_u4"]["ir"]-1.0)+3.0*(parsed["rev_b2_emit_split_c123_u4"]["ir"]-1.0))/4.0+parsed["rev_b2_rx_drain_u64"]["ir"]+0.015
 assert _c_w1<=coll["c_hit_hi_u4"] and (0.75*_c_w1+0.25*_c_spl)<=coll["c_hit_hi_u2"]
 assert parsed["rev_b2_flush_x10"]["raw_opcodes"]==20 and "b .Lb2i_fl16_ret" in parsed["rev_b2_flush_x10"]["insns"] and "b .Lb2c_fl16_ret" in parsed["rev_b2_flush_x10"]["insns"]
 assert parsed["rev_b2_flush_w8_x10"]["raw_opcodes"]==28 and parsed["rev_b2_flush_x11"]["raw_opcodes"]==16 and "b .Lspi_fl10_r" in parsed["rev_b2_flush_w8_x10"]["insns"] and "b .Lspi_fl11_r" in parsed["rev_b2_flush_x11"]["insns"]
 assert parsed["rev_b2_sieve_direct1_u2"]["raw_opcodes"]==7 and parsed["rev_b2_sieve_direct1_u2"]["items"]==2 and "add x_tx, x_tx, #16" in parsed["rev_b2_sieve_direct1_u2"]["insns"]
 assert "csinc x6, x_one, x_one, lo" in parsed["rev_b2_sieve_direct1_hit2"]["insns"] and "cset x6, hs" in parsed["rev_b2_sieve_direct1_hit2"]["insns"]
 assert parsed["rev_b2_build_subchord_hdr"]["raw_opcodes"]==12 and parsed["rev_b2_build_fiber_hdr"]["raw_opcodes"]==8
 assert "sub x3, x_p, x2" in parsed["rev_b2_build_subchord_hdr"]["insns"] and "sub x3, x_p, x2" in parsed["rev_b2_build_fiber_hdr"]["insns"]
 assert (parsed["rev_b2_build_hist_w1_stadd"]["ir"]+parsed["rev_b2_build_scatter_w1_u8"]["ir"])<=coll["c_list_a_w1"]
 assert (parsed["rev_b2_build_hist_u8_stadd"]["ir"]+parsed["rev_b2_build_scatter_r1_u8"]["ir"])<=coll["c_list_a_b2"]
 assert "ubfiz x3, x0, #59, #5" in parsed["rev_b2_build_scatter_w1_u8"]["insns"] and "ubfiz x3, x0, #59, #5" in parsed["rev_b2_build_scatter_r1_u8"]["insns"] and "ubfiz x11, x_k, #59, #5" in parsed["rev_b2_build_direct1_u8"]["insns"]
 assert (parsed["rev_build_prefix_u64_8x"]["ir"]+parsed["rev_b2_memset_ff_u16"]["ir"])<=coll["c_lbkt"] and parsed["rev_b2_memset_ff_u16"]["ir"]<=coll["c_memset_dir"]
 assert parsed["rev_b2_build_direct1_u8"]["ir"]<=11.25 and parsed["rev_b2_build_direct1_ovf"]["ir"]<=4.25 and "b .Lbdir_done" in parsed["rev_b2_build_direct1_ovf"]["insns"] and "ret" in parsed["rev_b2_build_direct1_ovf"]["insns"] and all(f".Lbdir_ovf{j}" in all_defs and f".Lbdir_ret{j}" in all_defs for j in range(8))
 assert "add x4, x_ovf_adj, x1, lsl #4" in parsed["rev_b2_sieve_direct1_hit2"]["insns"] and "add x4, x_ovf_adj, x2, lsl #4" in parsed["rev_b2_sieve_direct1_hit2"]["insns"] and "add x4, x_ovf_adj, x1, lsl #4" in parsed["rev_b2_sieve_direct1_tail1"]["insns"]
 assert parsed["ntt61_r4_u2"]["ir"]*9.0<=37.20<37.50 and "b.ne .Lntt_r4_u2" in parsed["ntt61_r4_u2"]["insns"]
 assert parsed["lehmer32_step"]["ir"]<=16.0 and parsed["ntt62_butterfly"]["ir"]*12.0<=120.0<168.0
 import w8_verifiers as _W8V
 _v06_k=_W8V.verify_v06_kernels_extra(parsed,coll,la,r_bls,NTT3_PRIMES)
 MASK=(1<<64)-1;OID_M=(1<<56)-1;BIT58=1<<58;CLR58=MASK^BIT58
 _oid=0x123456789abc
 for _m in range(4):
  for _c in range(4):
   for _delta in range(32):
    _xt0=((_delta<<59)|(_c<<56))&MASK
    for _r in range(32):
     _rec=((_r<<59)|(_m<<56)|_oid)&MASK
     _s=_rec+_xt0
     _carry=1 if _s>MASK else 0
     _out=(_s&CLR58)&MASK
     assert _carry==(1 if (_r+_delta)>=32 else 0)
     assert (_out>>59)==((_r+_delta)&31) and ((_out>>58)&1)==0 and ((_out>>56)&3)==((_m+_c)&3) and (_out&OID_M)==_oid
 _scg=parsed["rev_scan_g_u8_fused_128"];_scf=parsed["rev_scan_f_u8_128"];_spz=parsed["rev_pad_zero_u8_128"]
 assert parsed["rev_scan_zero_wrap"]["raw_opcodes"]==2 and _scg["items"]==128 and _scg["raw_opcodes"]==68
 assert _scf["items"]==128 and _scf["raw_opcodes"]==21 and _spz["items"]==128 and _spz["raw_opcodes"]==10
 assert _pw("rev_scan_g_u8_fused_128",_scg["insns"][48:63]).isdisjoint({f"r{i}" for i in range(16)})
 _sga_eff=float(sga_mean)if sga_mean is not None else 5e-5
 _p_hit128=min(1.0,128.0*_sga_eff)
 scan_cold_ir=_p_hit128*(1100.0/128.0)+_p_hit128*(1200.0/128.0)
 scan_tot=(68.0/128.0)+(21.0/128.0)+(31.0/lb_live)*(10.0/128.0)+2.0/lb_live+scan_cold_ir
 assert scan_tot<=coll["c_scan_zero"],f"scan_tot={scan_tot:.4f} > c_scan_zero={coll['c_scan_zero']}"
 _qg=float(coll["det"]["b1_q_g"]);_qf=float(coll["det"]["b1_q_f"])
 _rmax_g=max(round(b/_qg)/b for b in [33.0+i/16.0 for i in range(625)])
 _rmax_f=max(round(b/_qf)/b for b in [33.0+i/16.0 for i in range(625)])
 _tg_u8=tuple(coll["det"]["tg_u8"]);_tf_u8=int(coll["det"]["tf_u8"])
 _bg_u8=tuple(128-t for t in _tg_u8)
 assert max(_bg_u8[c]+math.floor((690-12*c)*_rmax_g)for c in range(3))<=252<=255,"g-side uint8_t overflow!"
 assert math.floor(630*_rmax_f)<=240<=255,"f-side uint8_t overflow!"
 _X8080=0x8080808080808080
 for _c in range(3):
  _T=_tg_u8[_c];_b8=_bg_u8[_c]
  for _s1 in (0,20,60,100):
   for _s2 in (0,max(0,_T-_s1-1),max(0,_T-_s1),min(90,195-_s1)):
    _v=_b8+_s1+_s2;assert _v<=255 and ((_v&0x80)!=0)==(_s1+_s2>=_T)
 _NC=3*128;_S_b1=[_bg_u8[i//128]for i in range(_NC)];_S_b2=[[0]*_NC for _ in range(32)];_rng_2p=random.Random(20261125)
 for _n_act in (16,11,16):
  _r_last=_n_act-1;_mg_sb=[[0]*(_NC//64)for _ in range(_n_act)]
  _g1=[[_rng_2p.randrange(0,55)for _ in range(_NC)]for _ in range(_n_act)]
  _g2=[[_rng_2p.randrange(0,35)for _ in range(_NC)]for _ in range(16)]
  _f1=[[_rng_2p.randrange(0,35)for _ in range(_NC)]for _ in range(_n_act)]
  _f2=[[_rng_2p.randrange(0,20)for _ in range(_NC)]for _ in range(16)]
  for _r in range(16):
   for _i in range(_NC):_S_b2[_r][_i]+=_g2[_r][_i]
  for _r in range(_n_act):
   for _i in range(_NC):_S_b1[_i]+=_g1[_r][_i]
   for _c in range(3):
    _xb_nxt=0 if _r==_r_last else _bg_u8[_c]
    for _sub in range(2):
     _blk=2*_c+_sub;_mb_w=0
     for _jj in range(64):
      _idx=64*_blk+_jj;_v=_S_b1[_idx]+_S_b2[_r][_idx];assert _v<=255
      _S_b1[_idx]=_xb_nxt;_S_b2[_r][_idx]=0
      if _v&0x80:_mb_w|=1<<_jj
     _mg_sb[_r][_blk]=_mb_w
  for _r in range(_n_act,16):_S_b2[_r]=[0]*_NC
  assert all(v==0 for v in _S_b1)and all(v==0 for row in _S_b2[:16] for v in row)
  for _r in range(16):
   for _i in range(_NC):_S_b2[_r][_i]+=_f2[_r][_i]
  for _r in range(_n_act):
   for _i in range(_NC):_S_b1[_i]+=_f1[_r][_i]
   _got_r=[]
   for _c in range(3):
    _xf_nxt=_bg_u8[_c]if _r==_r_last else 0;_w0,_w1=_mg_sb[_r][2*_c],_mg_sb[_r][2*_c+1]
    if (_w0|_w1)!=0:
     for _sub,_ww in ((0,_w0),(64,_w1)):
      for _jj in range(64):
       if (_ww>>_jj)&1:
        _idx=128*_c+_sub+_jj
        if _S_b1[_idx]+_S_b2[_r][_idx]>=_tf_u8:_got_r.append(_idx)
    for _jj in range(128):_S_b1[128*_c+_jj]=_xf_nxt;_S_b2[_r][128*_c+_jj]=0
   _exp_r=[_i for _i in range(_NC)if(_g1[_r][_i]+_g2[_r][_i]>=_tg_u8[_i//128])and(_f1[_r][_i]+_f2[_r][_i]>=_tf_u8)]
   assert _got_r==_exp_r
  for _r in range(_n_act,16):_S_b2[_r]=[0]*_NC
  assert all(_S_b1[_i]==_bg_u8[_i//128]for _i in range(_NC))and all(v==0 for row in _S_b2[:16] for v in row)
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
 assert parsed["row_red5_z12"]["ir"]<=58.0
 assert(6*parsed["row_red5_z12"]["raw_opcodes"]+240+16)<=la["c_row_chi"]
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
 assert all(p<(1<<61)and 8*p<(1<<64)and(p-1)%(3*(1<<52))==0 and _isp64(p)for p in NTT3_PRIMES)
 assert NTT3_PROD>(1<<54)*((1<<64)-1)**2
 _p1,_p2,_p3=NTT3_PRIMES;_i12=pow(_p1,-1,_p2);_i123=pow((_p1*_p2)%_p3,-1,_p3)
 for _ in range(16):
  _vv=rng.randrange(NTT3_PROD);_r1,_r2,_r3=_vv%_p1,_vv%_p2,_vv%_p3
  _k2=((_r2-_r1)%_p2*_i12)%_p2;_v12=_r1+_k2*_p1;_k3=((_r3-_v12%_p3)%_p3*_i123)%_p3
  assert _v12+_k3*_p1*_p2==_vv
 def _r4_r3_ntt(arr,g,p,inv=False):
  n=len(arr);p2=2*p;p4=4*p;p8=8*p;a=[x%p for x in arr]
  def _sm(v,w,wp):
   assert 0<=v<p8;t=(v*w-((v*wp)>>64)*p)&MASK;assert t<p2;return t
  def _c2(x):
   x=x-p4 if x>=p4 else x;return x-p2 if x>=p2 else x
  rads=[];rem=n
  while rem%4==0:rads.append(4);rem//=4
  if rem in (2,3):rads.append(rem);rem=1
  assert rem==1
  revs=[[0]];nbs=[1]
  for r in rads:
   nb=nbs[-1];revs.append([idx+d*nb for idx in revs[-1] for d in range(r)]);nbs.append(nb*r)
  if not inv:
   jv=pow(g,(p-1)//4,p);jp=(jv<<64)//p;m=n
   for st,r in enumerate(rads):
    m//=r;nb=nbs[st];rv=revs[st];ws=pow(g,(p-1)//(r*nb),p)
    if r==4:
     for k in range(nb):
      w1=pow(ws,rv[k],p);s1=(w1<<64)//p;w2=(w1*w1)%p;s2=(w2<<64)//p;w3=(w2*w1)%p;s3=(w3<<64)//p;b0=4*k*m
      for j0 in range(0,m,2):
       for jj in range(j0,min(j0+2,m)):
        a0,a1,a2,a3=a[b0+jj],a[b0+m+jj],a[b0+2*m+jj],a[b0+3*m+jj]
        x0=_c2(a0);x1=_sm(a1,w1,s1);x2=_sm(a2,w2,s2);x3=_sm(a3,w3,s3)
        x8=x0+x2;x10=x0+p2-x2;x11=x1+x3;x12=_sm(x1+p2-x3,jv,jp)
        a[b0+jj],a[b0+m+jj],a[b0+2*m+jj],a[b0+3*m+jj]=x8+x11,x10+x12,x8+p4-x11,x10+p2-x12
    elif r==3:
     w3=pow(g,(p-1)//3,p);w32=(w3*w3)%p
     for k in range(nb):
      w1=pow(ws,rv[k],p);w2=(w1*w1)%p;b0=3*k*m
      for j in range(m):
       x0,x1,x2=a[b0+j]%p,(a[b0+m+j]*w1)%p,(a[b0+2*m+j]*w2)%p
       a[b0+j],a[b0+m+j],a[b0+2*m+j]=(x0+x1+x2)%p,(x0+x1*w3+x2*w32)%p,(x0+x1*w32+x2*w3)%p
    else:
     for k in range(nb):
      w1=pow(ws,rv[k],p);b0=2*k*m
      for j in range(m):x0,x1=a[b0+j]%p,(a[b0+m+j]*w1)%p;a[b0+j],a[b0+m+j]=(x0+x1)%p,(x0-x1)%p
   return[x%p for x in a]
  jv=pow(g,3*(p-1)//4,p);jp=(jv<<64)//p;m=1
  for st in range(len(rads)-1,-1,-1):
   r=rads[st];nb=nbs[st];rv=revs[st];ws=pow(g,(p-1)-(p-1)//(r*nb),p)
   if r==4:
    for k in range(nb):
     w1=pow(ws,rv[k],p);s1=(w1<<64)//p;w2=(w1*w1)%p;s2=(w2<<64)//p;w3=(w2*w1)%p;s3=(w3<<64)//p;b0=4*k*m
     for j in range(m):
      y0,y1,y2,y3=a[b0+j],a[b0+m+j],a[b0+2*m+j],a[b0+3*m+j];assert max(y0,y1,y2,y3)<p2
      t0=y0+y2;t2=y0+p2-y2;t1=y1+y3;t3=_sm(y1+p2-y3,jv,jp)
      a[b0+j],a[b0+m+j],a[b0+2*m+j],a[b0+3*m+j]=_c2(t0+t1),_sm(t2+t3,w1,s1),_sm(t0+p4-t1,w2,s2),_sm(t2+p2-t3,w3,s3)
   elif r==3:
    w3=pow(g,2*(p-1)//3,p);w32=(w3*w3)%p
    for k in range(nb):
     w1=pow(ws,rv[k],p);w2=(w1*w1)%p;b0=3*k*m
     for j in range(m):
      y0,y1,y2=a[b0+j]%p,a[b0+m+j]%p,a[b0+2*m+j]%p
      a[b0+j],a[b0+m+j],a[b0+2*m+j]=(y0+y1+y2)%p,((y0+y1*w3+y2*w32)*w1)%p,((y0+y1*w32+y2*w3)*w2)%p
   else:
    for k in range(nb):
     w1=pow(ws,rv[k],p);b0=2*k*m
     for j in range(m):y0,y1=a[b0+j]%p,a[b0+m+j]%p;a[b0+j],a[b0+m+j]=(y0+y1)%p,((y0-y1)*w1)%p
   m*=r
  ip=pow(n,-1,p);return[(x*ip)%p for x in a]
 for _nl in (16,48,64):
  _ac=[rng.randrange(1<<64)for _ in range(_nl//2)]+[0]*(_nl-_nl//2)
  _bc=[rng.randrange(1<<64)for _ in range(_nl//2)]+[0]*(_nl-_nl//2)
  _rf=[0]*_nl
  for _ii in range(_nl//2):
   for _jj in range(_nl//2):_rf[_ii+_jj]+=_ac[_ii]*_bc[_jj]
  _rs3=[]
  for _pi in NTT3_PRIMES:
   _gi=2
   while pow(_gi,(_pi-1)//2,_pi)==1 or pow(_gi,(_pi-1)//3,_pi)==1:_gi+=1
   _fa=_r4_r3_ntt(_ac,_gi,_pi,False);_fb=_r4_r3_ntt(_bc,_gi,_pi,False)
   _rs3.append(_r4_r3_ntt([(_fa[i]*_fb[i])%_pi for i in range(_nl)],_gi,_pi,True))
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
 return dict(asm_blocks_parsed=len(parsed),undef_asm_labels=len(undef_labels),per_kernel_reg_disjoint_ok=True,b2_buffered_128k_defer8_ok=True,b2_full_and_rowend_drains_ok=True,b2_bank_rot_62b_ok=True,b1_pass12_orbit0_ok=True,v06_extra=_v06_k,dsieve16_dfs_4d_ok=True,dsieve16_ambient_sign_canon_ok=True,uint16_swar_ok=True,ntt62_garner_exact=True,ntt3_r4_r3_shoup_conv_ok=True,ntt61_r4_3x3p_stage_ir=round(9.0*parsed["ntt61_r4_u2"]["ir"],2),f62_barrett_mac_64_exact=True,mont5_opcodes=mont5_ops,mont5_aapcs64_ok=True,mont5_sim_32_exact=True,mont2_opcodes=mont2_ops,mont2_sim_exact=True,row_red5_barrett_opcodes=parsed["row_red5_z12"]["raw_opcodes"],row_red5_barrett_47_exact=True,unit_run_4vec4_ir=round(parsed["unit_run_4vec4_z12"]["ir"],4),unit_run_1vec4_ir=round(parsed["unit_run_1vec4_z12"]["ir"],4),gen_run1_2v_ir=round(parsed["gen_run1_2v_z12"]["ir"],2),gen_run1_1v_ir=round(parsed["gen_run1_1v_z12"]["ir"],2),scan_zero_u8_fused=round(scan_tot,4),u8_no_overflow_det_ok=True,dsieve16_rank_552_exact=True,sq_lll_uniform_log2=round(math.log2(sq_lll_bound),2),boot_lll_uniform_log2=round(math.log2(bt_lll_bound),2),all_raw_opcodes_ok=True)
from w8_verifiers import missing_recon_z12
