import base64,json,lzma,math,random
from math import log,log2,exp,sqrt,ceil,floor,pi,log1p,gamma
from fractions import Fraction as F
from decimal import Decimal,localcontext
from collections import Counter
from array import array
from pathlib import Path
import descent_bn254_v6 as DESC
LN2=math.log(2.0)
GAMMA=0.5772156649015329
UNIT=20000.0
U_BN254=4965661367192848881
P_BN254=36*U_BN254**4+36*U_BN254**3+24*U_BN254**2+6*U_BN254+1
R_BN254=36*U_BN254**4+36*U_BN254**3+18*U_BN254**2+6*U_BN254+1
H_COEFFS=[3,-4,0,2,1,-2,1]
W_COEFFS=[0,-2,-1,0,1,-1]
SIGMA_NUM=[9,-4,-5,1,3,-2]
ZK_NUM=[[3,0,0,0,0,0],[0,3,0,0,0,0],[9,-4,-5,1,3,-2],[0,2,1,-2,0,1],[6,1,-1,-1,3,-1],[-6,6,6,0,-3,3],]
F_REC_COEFFS={8:[3,-4,-2,0,2,-2],6:[0,6,3,0,-3,3],4:[-12,0,0,0,0,0],2:[12,-6,-3,0,3,-3],0:[-5,4,2,0,-2,2],}
CONTENT_BITS=12.0+6.0*log2(3.0)
EXACT_PRIM_T4=0.683542
ALPHA_F_RAW= -20.010+0.08
ALPHA_G_RAW= -0.119+0.08
_DATA=None
_WIT=None
_P17_G_CACHE=None
rand_below=DESC.rand_below
_B62={c:i for i,c in enumerate("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")}
def b62d(s):
 neg=s.startswith("-");v=0
 for c in(s[1:]if neg else s):v=v*62+_B62[c]
 return-v if neg else v
def load_data():
 global _DATA
 if _DATA is None:
  p=Path(__file__).resolve().parent/"c6_data.b85"
  _DATA=json.loads(lzma.decompress(base64.b85decode(p.read_bytes())).decode())
  import boot_dp;_bs=boot_dp.load();_bs["off_hi"]=max(o for r in _bs["rec"]for o in r[3])/100.0+1.0;_bs["primes16"]=list(_primes_upto(1<<16));_bs["ng_hi"]=boot_dp.gate_norm_bits();_DATA["descent"]["boot_dp"]=_bs
  pq=_DATA["pool_quad"]
  if "bins_b62" in pq and "bins" not in pq:
   acc=0;bs=[]
   for s in pq["bins_b62"].split(";"):
    a,b,c,d=s.split();acc+=b62d(a);cp=b62d(c)
    bs.append([acc/100.0,(b62d(b)+42000)/100.0,cp,cp+b62d(d)])
   pq["bins"]=bs
 return _DATA
def load_witnesses():
 global _WIT
 if _WIT is None:
  p=Path(__file__).resolve().parent/"c6_witnesses.b85"
  _WIT=json.loads(lzma.decompress(base64.b85decode(p.read_bytes())).decode())
 return _WIT
C_HIT_LO=3.5699
C_B1_ROW=26.2
C_SB_ENTER=8.05
C_SB_ENTER_DN=8.05
C_SB_ENTER_X4=20.05
C_SB_ENTER_LL_SLOPE=3.50
C_SB_ENTER_LL_CONST=7.05
C_B2_INT=5.55
C_B2_INT_X4=6.55
C_B2_INT_LL=6.55
C_B2_BND_HIT=9.55
C_B2_BND_MISS=7.55
C_B2_BND1_HIT=8.55
C_B2_BND1_MISS=6.55
C_CAND=6.5
C_HIT_HI=3.0
C_SCAN_ZERO=0.65
SCAN_G=8
C_SCAN_G8=549.0/1024.0
C_LIST_B1=24.70
C_LIST_B2=27.70
C_LIST_LL=21.90
C_LBUCKET=4.75
C_LBUCKET_DN=2.75
C_LBUCKET_LL=0.60
B2_BANDS=((41.5,43.5,8,"x4"),(43.5,45.5,16,"x4"))
C_B2_PREF_EX=2.0
C_LAT_EMIT=13.15
C_LAT_APPLY=1.55
C_LAT_LEAF=C_LAT_EMIT+C_LAT_APPLY
C_LAT_FALLBACK_RES=0.0
SHARD_PAGES_EXP=46
SHARD_BYTES=float(1<<(SHARD_PAGES_EXP+14))
N_SC_EXP=35
N_SC_NODES=float(1<<N_SC_EXP)
SC_BLOCK_EXP=38
SC_BLOCK_BYTES=float(1<<SC_BLOCK_EXP)
SC_W_SUBBAND=64
SC_GUARD_DELTA=0.005
N_LIST_STREAM_KINDS=8
SC_DUMMY_SINK_BYTES=float(1<<50)
STORAGE_CAP_ABORT_PROB=0.0
C_LAT_LLL=2.0**34.5
C_LAT_EXC=2.0**31.2
A_TAB_SLOTS_A3=284
B0_ARENA_MAX=560
QA_MAX=0.85
QB_MAX=0.85
B_A_SHELLS=2
K_B1_PREF=256
B1_PREF_COV=1.0075
_WKF=512.0/729.0;_WBR=125.0/27.0;_WBL=8.0
C_UNIT_5=7.20
C_GEN_1=42.0
C_NORM_G=10200.0
C_NORM_F=18000.0
C_SRT_G_LEAF=115000.0
C_SB=68000.0
C_TDIV=2.0**23.45
C_BATCH=2.0**22.12
LI_PAD=1.001
Q_F_SAFETY=2.0
KMAX=4
FFT_A,FFT_B,FFT_XFER,FFT_LOCAL,FFT_COPY,FFT_ROOT,FFT_LEAF=8.5,1200.0,504.0,67108864,2.0,3.0,32768.0
NTT3=(69.75,119.0,90.0,151.5,371.0,139.25,271.5)
NTT3_XFER=512.0
SPLIT_Q,SPLIT_REC,SPLIT_C,SPLIT_IDEAL,NORM_BITS=70368744177664.0,56,33024.0,67108864.0,576
SPLIT_CAP,SPLIT_RES=112.0,84.0
W_ROW=160.5
NU_SW_CAP=0.000312
NU_L_CAP=0.000261
K5_CAP=9.211739
TAB_R_CAP=1.00
P_CAP_H3=0.0
T_RUN_LOG2=112.0
C_PROJ=28.0
C_FULL=97.0
C_Z64=40.0
C_RED9=270.0
BW_SIGMA=0.0078125
PREP_OWNER=2.0**80
PREP_GUARD=2.0**80
PREP_JOINT=2.0**80
LA_RESERVE=2.0**78
G2_OUTPUT_BOUND=69456896.0
OP=dict(lf=61.0,lg=61.0,exc=0.76,W=160,p_min=24,p_cut=41.5,p_mid=45.5,p_lat=45.5,lx2=61.0,slack=114.0,slack_b=None,w_sb=64,c_tail_z6=740.0,margin=0.005,work_allow=1.0,)
_RHO_COEFS=None
def _build_rho_coefs(kmax=33,nt=140):
 with localcontext()as ctx:
  ctx.prec=70
  zero,one=Decimal(0),Decimal(1)
  coefs=[None,[one]+[zero]*(nt-1)]
  for k in range(2,kmax+1):
   prev,c=coefs[k-1],[zero]*nt
   for m in range(nt-1):
    c[m+1]=(prev[m]+m*c[m])/(k*(m+1))
   c[0]=prev[0]-sum(c[1:])
   coefs.append(c)
  out=[None]
  for c in coefs[1:]:
   big=max(abs(x)for x in c)
   while len(c)>1 and abs(c[-1])<big*Decimal("1e-18"):
    c=c[:-1]
   out.append(tuple(float(x)for x in reversed(c)))
  return tuple(out)
def rho(u):
 global _RHO_COEFS
 if u<=1.0:
  return 1.0
 if _RHO_COEFS is None:
  _RHO_COEFS=_build_rho_coefs()
 k=int(u)
 if u>k:k+=1
 if k>=len(_RHO_COEFS):
  return 0.0
 x=k-u
 v=0.0
 for c in _RHO_COEFS[k]:
  v=v*x+c
 return v if v>0.0 else 0.0
_RHO_TAB_1024=None
def rho_tab_1024():
 global _RHO_COEFS,_RHO_TAB_1024
 if _RHO_TAB_1024 is not None:return _RHO_TAB_1024
 if _RHO_COEFS is None:_RHO_COEFS=_build_rho_coefs()
 rc=_RHO_COEFS;out=[1.0]*1025
 for k in range(2,33):
  cs=rc[k]
  for j in range(1,1025):
   x=(1024-j)/1024.0;v=0.0
   for c in cs:v=v*x+c
   out.append(v if v>0.0 else 0.0)
 out.append(rho(32.0+1.0/1024.0))
 _RHO_TAB_1024=out
 return out
def _rho_fast(u):
 if u<=1.0:return 1.0
 rt=rho_tab_1024()
 f=u*1024.0;i=int(f)
 if i+1>=len(rt):return rho(u)
 return rt[i]+(f-i)*(rt[i+1]-rt[i])
def kappa_rho(u,t):
 ru=rho(u)
 return rho(u-t)/ru if ru>0.0 else 0.0
def _primes_upto(n):
 sieve=bytearray(b"\x01")*(n+1)
 sieve[:2]=b"\x00\x00"
 for p in range(2,int(sqrt(n))+1):
  if sieve[p]:
   sieve[p*p::p]=b"\x00"*((n-p*p)//p+1)
 return[p for p in range(2,n+1)if sieve[p]]
def _series_coeffs(weights,normalizer,length):
 _sp=math.sumprod;L=length
 jw_rev=[j*weights[j]for j in range(L,-1,-1)]
 out=[exp(-normalizer)]+[0.0]*L
 j_lo=1
 while j_lo<=L and weights[j_lo]==0.0:j_lo+=1
 j_hi=L
 while j_hi>=j_lo and weights[j_hi]==0.0:j_hi-=1
 if j_lo>L:return out
 for n in range(j_lo,L+1):
  k0=0 if n<=j_hi else n-j_hi
  out[n]=_sp(jw_rev[L-(n-k0):L-j_lo+1],out[k0:n-j_lo+1])/n
 return out
class Retention:
 def __init__(self,B_bits,X_bits,prime_min=2**14,step=0.25,max_bits=600,tilt=0.05,explicit_bits=14,power_terms=8,odd=1):
  self.B,self.X,self.minimum,self.step=B_bits,X_bits,prime_min,step
  self.tilt,self.length=tilt,ceil(max_bits/step)
  self.names=("small","sieve_geo","large_geo","sieve_sf","large_sf")
  weights={k:[0.0]*(self.length+1)for k in self.names}
  totals=Counter()
  cutoff=min(float(explicit_bits),B_bits)
  for p in _primes_upto(int(2**cutoff)-1)[odd:]:
   b=log2(p)
   bi=round(b/step)
   q=exp(tilt*bi*step-LN2*b)
   groups=[("small",False)]if p<prime_min else[(("sieve" if b<X_bits else "large")+"_geo",False),(("sieve" if b<X_bits else "large")+"_sf",True)]
   for name,sf in groups:
    totals[name]+=log1p(q)if sf else-log1p(-q)
    qk=q
    for k in range(1,self.length+1):
     idx=k*bi
     if idx>self.length:break
     weights[name][idx]+=(-1 if sf and k%2==0 else 1)*qk/k
     qk*=q
     if qk<1e-25:break
  def continuum(left,right,pref):
   if left>=right:return
   w=step/4.0
   cells=ceil((right-left)/w)
   w=(right-left)/cells
   for i in range(cells):
    b=left+(i+0.5)*w
    bi=round(b/step)
    for k in range(1,power_terms+1):
     idx=k*bi
     if idx>self.length:break
     amt=exp(k*tilt*bi*step-(k-1)*LN2*b)*w/(k*b)
     if amt<1e-25:break
     sign= -1 if k%2==0 else 1
     grps=((pref,amt),)if pref=="small" else((pref+"_geo",amt),(pref+"_sf",sign*amt))
     for g,x in grps:
      totals[g]+=x
      weights[g][idx]+=x
  lm=log2(prime_min)
  continuum(cutoff,min(lm,B_bits),"small")
  continuum(max(cutoff,lm),min(X_bits,B_bits),"sieve")
  continuum(max(cutoff,X_bits),B_bits,"large")
  self.weights,self.normalizers=weights,{k:totals[k]for k in self.names}
  all_w=[sum(weights[k][j]for k in("small","sieve_geo","large_geo"))for j in range(self.length+1)]
  self.denominator=_series_coeffs(all_w,sum(totals[k]for k in("small","sieve_geo","large_geo")),self.length)
  self.selection_normalizer=exp(totals["sieve_sf"]+totals["large_sf"]-totals["sieve_geo"]-totals["large_geo"])
  self._c={};self._ac={}
 def co(self,k):
  if k not in self._c:
   self._c[k]=_series_coeffs(self.weights[k],self.normalizers[k],self.length)
  return self._c[k]
 def two(self,n,a,b=9**9,win=None,sf=False):
  a=min(a,n)
  if a<0 or b<0:return 0.0
  key=(n,b,win,sf)
  if key not in self._c:
   if key[1:]not in self._c:
    st=self.step;sm=self.co("small");D=min(self.length,2048);m=[0.0]*(D+1);tsm=D
    while tsm>0 and sm[tsm]<1e-22:tsm-=1
    for l,v in enumerate(self.co("large_sf" if sf else "large_geo")[:min(b,D)+1]):
     if v>1e-22 and(win is None or l*st<=win[2]or any(self.X*k-win[2]<=l*st<=win[1]*k+win[2]for k in range(1,win[0]+1))):
      for t in range(min(D-l,tsm)+1):m[l+t]+=v*sm[t]
    self._c[key[1:]]=m
   m=self._c[key[1:]];sv=self.co("sieve_sf" if sf else "sieve_geo");c=array("d");t=0.0
   for x in range(min(n,len(m)-1)+1):
    t+=m[x]*sv[n-x];c.append(t)
   self._c[key]=c
  c=self._c[key]
  return c[min(a,len(c)-1)]
 def cap(self,nb):return floor(ceil(nb)/log2(self.minimum))
 def at(self,nb,sl,sb=None):
  ak=(nb,sl,sb);v=self._ac.get(ak)
  if v is not None:return v
  n=round(nb/self.step);rb=(self.cap(nb)+1)/16.0
  r=self.selection_normalizer*self.two(n,floor((sl-rb)/self.step),9**9 if sb is None else floor((sb-rb-0.125)/self.step),None,True)/self.denominator[n]
  v=max(0.0,min(1.0,r));self._ac[ak]=v;return v
 def survivor_upper(self,nb,sl,ex=0.0,sb=None,lb=None,kmax=4):
  n=round(nb/self.step);a=floor((sl+ex)/self.step)
  return(self.two(n,a)if sb is None else self.two(n,a,floor((sb+0.25)/self.step),(kmax,lb,(self.cap(nb)+2)/16.0+0.25)))/self.denominator[n]
def _ptrim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def _pmul(a,b,mod):
 la=len(a);lb=len(b);c=[0]*(la+lb-1)
 for i in range(la):
  x=a[i]
  if x:
   for j in range(lb):
    y=b[j]
    if y:c[i+j]+=x*y
 c=[x%mod for x in c]
 while len(c)>1 and c[-1]==0:c.pop()
 return c
def _pmod(a,b,mod):
 a=[x%mod for x in a]
 while len(a)>1 and a[-1]==0:a.pop()
 lb=len(b)
 while lb>1 and b[lb-1]%mod==0:lb-=1
 bl=b[lb-1]%mod;inv=1 if bl==1 else pow(bl,-1,mod)
 bm=[b[i]%mod for i in range(lb-1)]
 while len(a)>=lb and not(len(a)==1 and a[0]==0):
  d=len(a)-lb;fac=a.pop()
  if inv!=1:fac=(fac*inv)%mod
  if fac:
   for i in range(lb-1):a[d+i]=(a[d+i]-fac*bm[i])%mod
  while len(a)>1 and a[-1]==0:a.pop()
 return a if a else[0]
def _pdiv(a,b,mod):
 a=_ptrim([x%mod for x in a]);b=_ptrim([x%mod for x in b])
 q=[0]*max(1,len(a)-len(b)+1);inv=pow(b[-1],-1,mod)
 while len(a)>=len(b)and not(len(a)==1 and a[0]==0):
  d=len(a)-len(b);fac=(a[-1]*inv)%mod;q[d]=fac
  for i in range(len(b)):a[d+i]=(a[d+i]-fac*b[i])%mod
  _ptrim(a)
 return _ptrim(q),a
def _pgcd(a,b,mod):
 a=_ptrim([x%mod for x in a]);b=_ptrim([x%mod for x in b])
 while not(len(b)==1 and b[0]==0):
  a,b=b,_pmod(a,b,mod)
 inv=pow(a[-1],-1,mod)
 return[(x*inv)%mod for x in a]
def _ppowmod(base,e,m,mod):
 if len(m)==7 and len(base)<=6 and m[6]%mod==1 and m[5]%mod==mod-2 and m[4]%mod==1 and m[3]%mod==2 and m[2]%mod==0 and m[1]%mod==mod-4 and m[0]%mod==3:return _ptrim(_ppow_h(base,e,mod))
 if len(m)==3 and m[2]%mod==1:
  m0,m1=(-m[0])%mod,(-m[1])%mod;r0,r1=1,0;b0,b1=(base[0]%mod if base else 0),(base[1]%mod if len(base)>1 else 0)
  while e>0:
   if e&1:t=r1*b1;r0,r1=(r0*b0+t*m0)%mod,(r0*b1+r1*b0+t*m1)%mod
   t=b1*b1;b0,b1=(b0*b0+t*m0)%mod,(2*b0*b1+t*m1)%mod;e>>=1
  return[r0,r1]if r1 else[r0]
 if len(m)==4 and m[3]%mod==1:
  m0,m1,m2=(-m[0])%mod,(-m[1])%mod,(-m[2])%mod;r0,r1,r2=1,0,0
  b0,b1,b2=(base[0]%mod if base else 0),(base[1]%mod if len(base)>1 else 0),(base[2]%mod if len(base)>2 else 0)
  while e>0:
   if e&1:c4=r2*b2;c3=r1*b2+r2*b1+c4*m2;r0,r1,r2=(r0*b0+c3*m0)%mod,(r0*b1+r1*b0+c4*m0+c3*m1)%mod,(r0*b2+r1*b1+r2*b0+c4*m1+c3*m2)%mod
   c4=b2*b2;c3=2*b1*b2+c4*m2;b0,b1,b2=(b0*b0+c3*m0)%mod,(2*b0*b1+c4*m0+c3*m1)%mod,(2*b0*b2+b1*b1+c4*m1+c3*m2)%mod;e>>=1
  return _ptrim([r0,r1,r2])
 res=[1];b=_pmod(base,m,mod)
 while e>0:
  if e&1:res=_pmod(_pmul(res,b,mod),m,mod)
  b=_pmod(_pmul(b,b,mod),m,mod);e>>=1
 return res
def _ddf_degrees(f,mod):
 f=_ptrim([x%mod for x in f]);inv=pow(f[-1],-1,mod);f=[(x*inv)%mod for x in f]
 df=[(i*f[i])%mod for i in range(1,len(f))]
 if len(_pgcd(f,df,mod))!=1:return None
 degs=[];h=[0,1];rem=f[:];d=1
 while len(rem)-1>=2*d:
  h=_ppowmod(h,mod,rem,mod)
  hx=h[:]+[0]*max(0,2-len(h));hx[1]=(hx[1]-1)%mod
  g=_pgcd(hx,rem,mod);gd=len(g)-1
  if gd>0:
   degs.extend([d]*(gd//d))
   rem,_=_pdiv(rem,g,mod)
   h=_pmod(h,rem,mod)
  d+=1
 if len(rem)-1>0:degs.append(len(rem)-1)
 return sorted(degs)
def _roots_mod(f,mod):
 rng=random.Random(12345)
 f=_ptrim([x%mod for x in f]);inv=pow(f[-1],-1,mod);f=[(x*inv)%mod for x in f]
 via_h=len(f)<=7 and _pmod(H_COEFFS,f,mod)==[0]
 out,stack=[],[f]
 while stack:
  cur=stack.pop();deg=len(cur)-1
  if deg==0:continue
  if deg==1:out.append((-cur[0])%mod);continue
  while True:
   a=rand_below(rng.getrandbits,mod)
   if via_h and 4<=deg<=5:pw=_pmod(_ppow_h([a,1],(mod-1)//2,mod),cur,mod)
   else:pw=_ppowmod([a,1],(mod-1)//2,cur,mod)
   pw[0]=(pw[0]-1)%mod
   g=_pgcd(pw,cur,mod)
   if 0<len(g)-1<deg:
    q,_=_pdiv(cur,g,mod)
    stack.append(g);stack.append(q);break
 return sorted(out)
def _brent(n):
 if n%2==0:return 2
 for c in(1,2,3,5,7,11,13):
  y,m,g,r,q=2,64,1,1,1
  while g==1:
   x=y
   for _ in range(r):y=(y*y+c)%n
   k=0
   while k<r and g==1:
    ys=y
    for _ in range(min(m,r-k)):y=(y*y+c)%n;q=(q*abs(x-y))%n
    g=math.gcd(q,n);k+=m
   r<<=1
  if g==n:
   while True:
    ys=(ys*ys+c)%n;g=math.gcd(abs(x-ys),n)
    if g>1:break
  if 1<g<n:return g
 return n
def _isp(n):
 if n<2:return False
 for p in(2,3,5,7,11,13,17,19,23,29,31,37):
  if n==p:return True
  if n%p==0:return False
 d,s=n-1,0
 while d%2==0:d>>=1;s+=1
 for a in((2,7,61)if n<4759123141 else(2,325,9375,28178,450775,9780504,1795265022)):
  if a%n==0:continue
  x=pow(a,d,n)
  if x in(1,n-1):continue
  for _ in range(s-1):
   x=(x*x)%n
   if x==n-1:break
  else:return False
 return True
_MHP={S:(1<<S,1<<(S-1),(1<<S)-1,S,2*S,3*S,4*S,5*S)for S in(32,44,88)}
def _mh_k(a,b,S=44):
 M,H,k,s1,s2,s3,s4,s5=_MHP[S]
 P=(a[0]+(a[1]<<s1)+(a[2]<<s2)+(a[3]<<s3)+(a[4]<<s4)+(a[5]<<s5))*(b[0]+(b[1]<<s1)+(b[2]<<s2)+(b[3]<<s3)+(b[4]<<s4)+(b[5]<<s5))
 c0=P&k
 if c0>=H:c0-=M;P+=M
 P>>=s1;c1=P&k
 if c1>=H:c1-=M;P+=M
 P>>=s1;c2=P&k
 if c2>=H:c2-=M;P+=M
 P>>=s1;c3=P&k
 if c3>=H:c3-=M;P+=M
 P>>=s1;c4=P&k
 if c4>=H:c4-=M;P+=M
 P>>=s1;c5=P&k
 if c5>=H:c5-=M;P+=M
 P>>=s1;c6=P&k
 if c6>=H:c6-=M;P+=M
 P>>=s1;c7=P&k
 if c7>=H:c7-=M;P+=M
 P>>=s1;c8=P&k
 if c8>=H:c8-=M;P+=M
 P>>=s1;c9=P&k
 if c9>=H:c9-=M;P+=M
 c10=P>>s1
 if c10:c9+=2*c10;c8-=c10;c7-=2*c10;c5+=4*c10;c4-=3*c10
 if c9:c8+=2*c9;c7-=c9;c6-=2*c9;c4+=4*c9;c3-=3*c9
 if c8:c7+=2*c8;c6-=c8;c5-=2*c8;c3+=4*c8;c2-=3*c8
 if c7:c6+=2*c7;c5-=c7;c4-=2*c7;c2+=4*c7;c1-=3*c7
 if c6:c5+=2*c6;c4-=c6;c3-=2*c6;c1+=4*c6;c0-=3*c6
 return(c0,c1,c2,c3,c4,c5)
def _mh_int(a,b):
 a0,a1,a2,a3,a4,a5=a;b0,b1,b2,b3,b4,b5=b
 c0=a0*b0;c1=a0*b1+a1*b0;c2=a0*b2+a1*b1+a2*b0;c3=a0*b3+a1*b2+a2*b1+a3*b0;c4=a0*b4+a1*b3+a2*b2+a3*b1+a4*b0
 c5=a0*b5+a1*b4+a2*b3+a3*b2+a4*b1+a5*b0;c6=a1*b5+a2*b4+a3*b3+a4*b2+a5*b1;c7=a2*b5+a3*b4+a4*b3+a5*b2;c8=a3*b5+a4*b4+a5*b3;c9=a4*b5+a5*b4;c10=a5*b5
 if c10:c9+=2*c10;c8-=c10;c7-=2*c10;c5+=4*c10;c4-=3*c10
 if c9:c8+=2*c9;c7-=c9;c6-=2*c9;c4+=4*c9;c3-=3*c9
 if c8:c7+=2*c8;c6-=c8;c5-=2*c8;c3+=4*c8;c2-=3*c8
 if c7:c6+=2*c7;c5-=c7;c4-=2*c7;c2+=4*c7;c1-=3*c7
 if c6:c5+=2*c6;c4-=c6;c3-=2*c6;c1+=4*c6;c0-=3*c6
 return(c0,c1,c2,c3,c4,c5)
def _pmul_h(a,b,mod):return[x%mod for x in _mh_int(a,b)]
def _ppow_h(a,e,mod):
 if e>0 and len(a)>=2 and a[0]==0 and a[1]==1 and not any(a[2:]):
  r0,r1,r2,r3,r4,r5=1,0,0,0,0,0
  for bit in bin(e)[2:]:
   c0=r0*r0;c1=2*r0*r1;c2=2*r0*r2+r1*r1;c3=2*(r0*r3+r1*r2);c4=2*(r0*r4+r1*r3)+r2*r2
   c5=2*(r0*r5+r1*r4+r2*r3);c6=2*(r1*r5+r2*r4)+r3*r3;c7=2*(r2*r5+r3*r4);c8=2*r3*r5+r4*r4;c9=2*r4*r5;c10=r5*r5
   if c10:c9+=2*c10;c8-=c10;c7-=2*c10;c5+=4*c10;c4-=3*c10
   if c9:c8+=2*c9;c7-=c9;c6-=2*c9;c4+=4*c9;c3-=3*c9
   if c8:c7+=2*c8;c6-=c8;c5-=2*c8;c3+=4*c8;c2-=3*c8
   if c7:c6+=2*c7;c5-=c7;c4-=2*c7;c2+=4*c7;c1-=3*c7
   if c6:c5+=2*c6;c4-=c6;c3-=2*c6;c1+=4*c6;c0-=3*c6
   r0,r1,r2,r3,r4,r5=c0%mod,c1%mod,c2%mod,c3%mod,c4%mod,c5%mod
   if bit=="1":
    c=r5;r0,r1,r2,r3,r4,r5=(-3*c)%mod,(r0+4*c)%mod,r1,(r2-2*c)%mod,(r3-c)%mod,(r4+2*c)%mod
  return[r0,r1,r2,r3,r4,r5]
 r=[1,0,0,0,0,0];b=[(a[i]%mod if i<len(a)else 0)for i in range(6)]
 while e>0:
  if e&1:r=_pmul_h(r,b,mod)
  b=_pmul_h(b,b,mod);e>>=1
 return r
def _pinv_h(a,mod):
 r0,r1=[x%mod for x in H_COEFFS],_ptrim([x%mod for x in a])
 t0,t1=[0],[1]
 while not(len(r1)==1 and r1[0]==0):
  q,rem=_pdiv(r0,r1,mod)
  qt=_pmul(q,t1,mod)
  L=max(len(t0),len(qt))
  t2=_ptrim([((t0[i]if i<len(t0)else 0)-(qt[i]if i<len(qt)else 0))%mod for i in range(L)])
  r0,r1,t0,t1=r1,rem,t1,t2
 inv=pow(r0[0],-1,mod)
 return([(x*inv)%mod for x in t0]+[0]*6)[:6]
def _norm_Kh_exact(coeffs_int):
 M=[[F(0)]*6 for _ in range(6)];cur=[F(c)for c in coeffs_int]
 for col in range(6):
  for r in range(6):M[r][col]=cur[r]
  top=cur[5];nxt=[F(0)]+cur[:5]
  for j in range(6):nxt[j]-=top*H_COEFFS[j]
  cur=nxt
 det=F(1)
 for col in range(6):
  piv=next(r for r in range(col,6)if M[r][col]!=0)
  if piv!=col:M[col],M[piv]=M[piv],M[col];det= -det
  pv=M[col][col];det*=pv
  for r in range(col+1,6):
   if M[r][col]!=0:
    fac=M[r][col]/pv
    for j in range(col,6):M[r][j]-=fac*M[col][j]
 return int(det)
def _pmul_h_q(a,b):
 res=[F(0)]*11
 for i in range(6):
  for j in range(6):res[i+j]+=a[i]*b[j]
 for k in range(10,5,-1):
  c=res[k]
  if c:
   res[k]=F(0)
   for j in range(6):res[k-6+j]-=c*H_COEFFS[j]
 return res[:6]
def verify_tower_algebra_and_crt():
 ord_ok=(pow(P_BN254,12,R_BN254)==1)and(pow(P_BN254,6,R_BN254)!=1)and(pow(P_BN254,4,R_BN254)!=1)
 y=[0,1,0,0,0,0]
 inv3=pow(3,-1,P_BN254)
 sig_p=[(c*inv3)%P_BN254 for c in SIGMA_NUM]
 y_p1=_ppow_h(y,P_BN254,P_BN254)
 y_p2=_ppow_h(y_p1,P_BN254,P_BN254)
 y_p3=_ppow_h(y_p2,P_BN254,P_BN254)
 y_p4=_ppow_h(y_p3,P_BN254,P_BN254)
 y_p6=_ppow_h(_ppow_h(y_p4,P_BN254,P_BN254),P_BN254,P_BN254)
 h_irred=(y_p6==y)and(y_p2!=y)and(y_p3!=y)
 sig_frob_ok=(y_p4==sig_p)
 sig_q=[F(c,3)for c in SIGMA_NUM]
 pows=[[F(1)]+[F(0)]*5,sig_q]
 for _ in range(4):pows.append(_pmul_h_q(pows[-1],sig_q))
 w_sig=[sum(W_COEFFS[k]*pows[k][j]for k in range(6))for j in range(6)]
 sigma_fixes_w=(w_sig==[F(c)for c in W_COEFFS])
 uw=[(U_BN254*c)%P_BN254 for c in W_COEFFS]
 uw_half=_ppow_h(uw,(P_BN254-1)//2,P_BN254)
 cur,acc=uw_half,[1,0,0,0,0,0]
 val_p6=val_p10=None
 for j in range(10):
  acc=_pmul_h(acc,cur,P_BN254)
  if j==5:val_p6=acc[:]
  if j==9:val_p10=acc[:]
  cur=_ppow_h(cur,P_BN254,P_BN254)
 minus_one=[P_BN254-1,0,0,0,0,0]
 g_irred=(val_p6==minus_one)
 tau_frob10=sig_frob_ok and(val_p10==minus_one)
 chi=pow(P_BN254,10,R_BN254)
 chi_expr=36*U_BN254**3+18*U_BN254**2+6*U_BN254+2
 chi_ok=(chi==chi_expr)and((chi*chi-chi+1)%R_BN254==0)and((pow(chi,3,R_BN254)+1)%R_BN254==0)
 inv_uw=_pinv_h(uw,P_BN254)
 f_eval=[0]*6;x2_pow=[1,0,0,0,0,0]
 for d in(0,2,4,6,8):
  term=_pmul_h([c%P_BN254 for c in F_REC_COEFFS[d]],x2_pow,P_BN254)
  f_eval=[(a+b)%P_BN254 for a,b in zip(f_eval,term)]
  x2_pow=_pmul_h(x2_pow,inv_uw,P_BN254)
 f_vanishes=all(c==0 for c in f_eval)
 nw,n_lead,n_const=abs(_norm_Kh_exact(W_COEFFS)),abs(_norm_Kh_exact(F_REC_COEFFS[8])),abs(_norm_Kh_exact(F_REC_COEFFS[0]))
 h_mod_r_degs=_ddf_degrees(H_COEFFS,R_BN254)
 rts_r=_roots_mod(H_COEFFS,R_BN254)
 g_mod_r_degs,f_mod_r_degs=[],[]
 for s in rts_r:
  ws=sum(W_COEFFS[j]*pow(s,j,R_BN254)for j in range(6))%R_BN254
  g_mod_r_degs.append(_ddf_degrees([(-1)%R_BN254,0,(U_BN254*ws)%R_BN254],R_BN254))
  fp=[0]*9
  for d,cs in F_REC_COEFFS.items():
   fp[d]=sum(cs[j]*pow(s,j,R_BN254)for j in range(6))%R_BN254
  f_mod_r_degs.append(_ddf_degrees(fp,R_BN254))
 crt_ok=(h_mod_r_degs==[1,1,1,1,1,1]and len(rts_r)==6
 and all(gd==[2]for gd in g_mod_r_degs)
 and all(fd==[2,2,2,2]for fd in f_mod_r_degs))
 T4=load_data()["t4_geom"]
 S=T4["SMAT"]
 def iso_err(G):
  return max(abs(sum(S[a][i]*G[a][b]*S[b][j]for a in range(6)for b in range(6))-G[i][j])for i in range(6)for j in range(6))
 iso_a,iso_b=iso_err(T4["GA"]),iso_err(T4["GB"])
 gcd_r_p6_ok=(math.gcd(R_BN254,pow(P_BN254,6,R_BN254)-1)==1)
 c0_t=[(17*k+3)%P_BN254 for k in range(6)]
 c1_t=[(29*k+11)%P_BN254 for k in range(6)]
 a_t=[(7*k+5)%P_BN254 for k in range(6)]
 inv_c0=_pinv_h(c0_t,P_BN254)
 R_boot=_pmul_h(c1_t,inv_c0,P_BN254)
 b_t=_pmul_h(R_boot,a_t,P_BN254)
 cof_t=_pmul_h(a_t,inv_c0,P_BN254)
 boot_id_ok=(_pmul_h(cof_t,c0_t,P_BN254)==a_t and _pmul_h(cof_t,c1_t,P_BN254)==b_t)
 def _is_prime64(n):
  d,s=n-1,0
  while d%2==0:d//=2;s+=1
  for a in(2,325,9375,28178,450775,9780504,1795265022):
   if a%n==0:continue
   x=pow(a,d,n)
   if x in(1,n-1):continue
   for _ in range(s-1):
    x=(x*x)%n
    if x==n-1:break
   else:return False
  return True
 ntt3_prod=1
 for p_ntt in NTT3_PRIMES:ntt3_prod*=p_ntt
 def _mont_u64_ok(p):
  M64=(1<<64)-1;R2=pow(1<<64,2,p);pinv=(-pow(p,-1,1<<64))&M64
  rng_m=random.Random(p&0xFFFFFFFF)
  pairs=[(p-1,p-1),(p-2,p-2),(2*p-1,2*p-1),(p-1,2*p-1),(0,p-1),(1,p-2)]+[(rng_m.randrange(2*p),rng_m.randrange(2*p))for _ in range(32)]
  for a,b in pairs:
   ab=a*b;m=((ab&M64)*pinv)&M64;s=ab+m*p
   if s>=(1<<128):return False
   t=s>>64
   if t>=2*p:return False
   if t>=p:t-=p
   if t!=(a*b*pow(1<<64,-1,p))%p:return False
  return True
 ntt3_mont_ok=all(_mont_u64_ok(p)for p in NTT3_PRIMES)
 def _shoup_bf_ok(p):
  M64=(1<<64)-1;rng_s=random.Random((p^0x55AA)&0xFFFFFFFF)
  for _ in range(64):
   u=rng_s.randrange(p);v=rng_s.randrange(p);w=rng_s.randrange(p)
   ws=(w<<64)//p;q=(v*ws)>>64;t=(v*w-q*p)&M64
   if t>=2*p:return False
   if t>=p:t-=p
   a=u+t
   if a>=p:a-=p
   b=u-t
   if b<0:b+=p
   if a!=(u+v*w)%p or b!=(u-v*w)%p:return False
  return True
 ntt3_shoup_ok=all(_shoup_bf_ok(p)for p in NTT3_PRIMES)
 p1,p2,p3=NTT3_PRIMES;inv12=pow(p1,-1,p2);inv123=pow((p1*p2)%p3,-1,p3);rng_g=random.Random(20261120)
 ntt3_garner_ok=True
 for _ in range(32):
  val=rng_g.randrange(ntt3_prod)
  r1,r2,r3=val%p1,val%p2,val%p3
  k2=((r2-r1)%p2*inv12)%p2;v12=r1+k2*p1
  k3=((r3-v12%p3)%p3*inv123)%p3;rec=v12+k3*p1*p2
  if rec!=val:ntt3_garner_ok=False
 def _r4_r3_shoup_ntt(arr,w_root,p,inv=False):
  n=len(arr);p2=2*p;M64=(1<<64)-1;wb=pow(w_root,-1,p)if inv else w_root
  def _sm(v,w,wp):
   q=(v*wp)>>64;t=(v*w-q*p)&M64
   if t>=p2:return-1
   return t
  rem=n;rads=[]
  if rem%3==0:rads.append(3);rem//=3
  while rem%4==0:rads.append(4);rem//=4
  if rem==2:rads.append(2);rem//=2
  a=[0]*n
  for i in range(n):
   ix=i;rv=0
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
      if t1<0 or t2<0 or t3<0:return None
      u0=t0+t2;u0=u0-p2 if u0>=p2 else u0;u1=t0+p2-t2;u1=u1-p2 if u1>=p2 else u1
      v0=t1+t3;v0=v0-p2 if v0>=p2 else v0;v1=_sm(t1+p2-t3,J,Jp)
      if v1<0:return None
      y0=u0+v0;y0=y0-p2 if y0>=p2 else y0;y2=u0+p2-v0;y2=y2-p2 if y2>=p2 else y2
      y1=u1+v1;y1=y1-p2 if y1>=p2 else y1;y3=u1+p2-v1;y3=y3-p2 if y3>=p2 else y3
      a[k0+j],a[k0+j+st],a[k0+j+2*st],a[k0+j+3*st]=y0,y1,y2,y3
   elif r==2:
    for _ in range(st):tw.append((w1,(w1<<64)//p));w1=(w1*ws)%p
    for k0 in range(0,n,m):
     for j in range(st):
      w1,w1p=tw[j];u=a[k0+j];v=_sm(a[k0+j+st],w1,w1p)
      if v<0:return None
      y0=u+v;y0=y0-p2 if y0>=p2 else y0;y1=u+p2-v;y1=y1-p2 if y1>=p2 else y1
      a[k0+j],a[k0+j+st]=y0,y1
   else:
    for _ in range(st):w2=(w1*w1)%p;tw.append((w1,(w1<<64)//p,w2,(w2<<64)//p));w1=(w1*ws)%p
    for k0 in range(0,n,m):
     for j in range(st):
      w1,w1p,w2,w2p=tw[j];t0=a[k0+j];t1=_sm(a[k0+j+st],w1,w1p);t2=_sm(a[k0+j+2*st],w2,w2p)
      if t1<0 or t2<0:return None
      s12=t1+t2;s12=s12-p2 if s12>=p2 else s12;y0=t0+s12;y0=y0-p2 if y0>=p2 else y0
      t1a=_sm(t1,W3,W3p);t2a=_sm(t2,W32,W32p);t1b=_sm(t1,W32,W32p);t2b=_sm(t2,W3,W3p)
      if min(t1a,t2a,t1b,t2b)<0:return None
      sa=t1a+t2a;sa=sa-p2 if sa>=p2 else sa;y1=t0+sa;y1=y1-p2 if y1>=p2 else y1
      sb=t1b+t2b;sb=sb-p2 if sb>=p2 else sb;y2=t0+sb;y2=y2-p2 if y2>=p2 else y2
      a[k0+j],a[k0+j+st],a[k0+j+2*st]=y0,y1,y2
  if inv:in_p=pow(n,-1,p);return[(x*in_p)%p for x in a]
  return[x%p for x in a]
 def _ntt_r4_r3_conv_exact(n_len):
  a_c=[rng_g.randrange(1<<64)for _ in range(n_len//2)]+[0]*(n_len-n_len//2)
  b_c=[rng_g.randrange(1<<64)for _ in range(n_len//2)]+[0]*(n_len-n_len//2)
  ref=[0]*n_len
  for i in range(n_len//2):
   for j in range(n_len//2):ref[i+j]+=a_c[i]*b_c[j]
  rs=[]
  for p_i in NTT3_PRIMES:
   g_i=2
   while pow(g_i,(p_i-1)//2,p_i)==1 or pow(g_i,(p_i-1)//3,p_i)==1:g_i+=1
   w_i=pow(g_i,(p_i-1)//n_len,p_i)
   fa=_r4_r3_shoup_ntt(a_c,w_i,p_i,False);fb=_r4_r3_shoup_ntt(b_c,w_i,p_i,False)
   if fa is None or fb is None:return False
   fc=[(fa[i]*fb[i])%p_i for i in range(n_len)]
   ic=_r4_r3_shoup_ntt(fc,w_i,p_i,True)
   if ic is None:return False
   rs.append(ic)
  for i in range(n_len):
   r1,r2,r3=rs[0][i],rs[1][i],rs[2][i]
   k2=((r2-r1)%p2*inv12)%p2;v12=r1+k2*p1;k3=((r3-v12%p3)%p3*inv123)%p3
   if v12+k3*p1*p2!=ref[i]:return False
  return True
 ntt3_r4_r3_shoup_conv_ok=(_ntt_r4_r3_conv_exact(16)and _ntt_r4_r3_conv_exact(24))
 ntt3_ok=(len(NTT3_PRIMES)==3 and all(p<(1<<62)and(p-1)%(3*(1<<55))==0 and _is_prime64(p)for p in NTT3_PRIMES)
 and ntt3_prod>(1<<54)*((1<<64)-1)**2 and ntt3_mont_ok and ntt3_shoup_ok and ntt3_garner_ok and ntt3_r4_r3_shoup_conv_ok)
 ok=(ord_ok and h_irred and g_irred and tau_frob10 and sigma_fixes_w and chi_ok and f_vanishes
 and crt_ok and gcd_r_p6_ok and boot_id_ok and ntt3_ok and(nw==216)and(n_lead==729)and(n_const==729)
 and(iso_a<1e-10)and(iso_b<1e-10))
 return dict(ok=ok,ord_r_p_12=ord_ok,h_irred_mod_p=h_irred,g_rec_irred_mod_p6=g_irred,tau_is_frob10=tau_frob10,sigma_fixes_w_exact=sigma_fixes_w,chi_192bit_ok=chi_ok,f_rec_vanishes_at_inv_uw=f_vanishes,gcd_r_p6_minus_1_is_1=gcd_r_p6_ok,boot_identity_ok=boot_id_ok,
 ntt3_primes_verified=ntt3_ok,ntt3_montgomery_no_carry_ok=ntt3_mont_ok,ntt3_shoup_butterfly_ok=ntt3_shoup_ok,ntt3_garner_crt_ok=ntt3_garner_ok,ntt3_r4_r3_shoup_conv_ok=ntt3_r4_r3_shoup_conv_ok,
 ntt3_prod_bits=round(log2(ntt3_prod),3),norm_w=nw,norm_f_lead=n_lead,norm_f_const=n_const,mod_r_h_degs=h_mod_r_degs,mod_r_g_degs=g_mod_r_degs,mod_r_f_degs=f_mod_r_degs,mod_r_crt_all_deg2_ok=crt_ok,sigma_iso_err_Ga=iso_a,sigma_iso_err_Gb=iso_b,
 )
def build_exact_int_geom(T4=None):
 if T4 is None:T4=load_data()["t4_geom"]
 GA,UMATS,S=T4["GA"],T4["UMATS"],T4["SMAT"]
 def mm(A,B):return[[sum(A[i][k]*B[k][j]for k in range(6))for j in range(6)]for i in range(6)]
 def tr(A):return[[A[j][i]for j in range(6)]for i in range(6)]
 SCALE=1<<16
 G0=[[round(GA[i][j]*SCALE)for j in range(6)]for i in range(6)]
 for i in range(6):
  for j in range(i+1,6):G0[j][i]=G0[i][j]
 S2=mm(S,S)
 G1=mm(tr(S),mm(G0,S));G2=mm(tr(S2),mm(G0,S2))
 G_sym=[[G0[i][j]+G1[i][j]+G2[i][j]for j in range(6)]for i in range(6)]
 assert mm(tr(S),mm(G_sym,S))==G_sym
 D_mats=[[[mm(tr(U),mm(G_sym,U))[i][j]-G_sym[i][j]for j in range(6)]for i in range(6)]for U in UMATS]
 return G_sym,D_mats,UMATS,S,S2,3*SCALE
def verify_exhaustive_fd(R2=50.0):
 T4=load_data()["t4_geom"]
 G_sym,D_mats,UMATS,S,S2,G_SC=build_exact_int_geom(T4)
 def _sp(M):return tuple(tuple((c,M[r][c])for c in range(6)if M[r][c])for r in range(6))
 S_sp,S2_sp=_sp(S),_sp(S2);U_sps=tuple(_sp(U)for U in UMATS)
 def mv_sp(sp,v):return tuple(sum(mv*v[c]for c,mv in row)for row in sp)
 def qf_i(Mat,v):
  s=0
  for i in range(6):
   vi=v[i]
   if vi:
    s+=Mat[i][i]*vi*vi
    for j in range(i+1,6):
     if v[j]:s+=2*Mat[i][j]*vi*v[j]
  return s
 T_lim=int(round(R2*G_SC))
 from fractions import Fraction as Frac
 from math import isqrt
 D_f=[Frac(0)]*6;L_f=[[Frac(1 if i==j else 0)for j in range(6)]for i in range(6)]
 for c in range(6):
  d=Frac(G_sym[c][c])-sum(L_f[c][k]*L_f[c][k]*D_f[k]for k in range(c))
  assert d>0;D_f[c]=d
  for r in range(c+1,6):L_f[r][c]=(Frac(G_sym[r][c])-sum(L_f[r][k]*L_f[c][k]*D_f[k]for k in range(c)))/d
 pts,v=[],[0]*6
 def enum_q(idx,rem):
  sh=sum(L_f[j][idx]*v[j]for j in range(idx+1,6));rad2=rem/D_f[idx]
  r_hi=Frac(isqrt(rad2.numerator*rad2.denominator)+1,rad2.denominator);di=D_f[idx]
  for x in range((-sh-r_hi).__floor__(),(-sh+r_hi).__ceil__()+1):
   d=x+sh;nrem=rem-di*d*d
   if nrem>=0:
    v[idx]=x
    if idx==0:
     if any(v):pts.append(tuple(v))
    else:enum_q(idx-1,nrem)
 enum_q(5,Frac(T_lim))
 def is_can_int(p):
  if(-p[0],-p[1],-p[2],-p[3],-p[4],-p[5])<p:return False
  for fi in range(6):
   q=qf_i(D_mats[fi],p)
   if q<0:return False
   if q==0:
    uv=mv_sp(U_sps[fi],p)
    for w in(uv,mv_sp(S_sp,uv),mv_sp(S2_sp,uv)):
     nw=(-w[0],-w[1],-w[2],-w[3],-w[4],-w[5])
     if(nw if nw<w else w)<p:return False
  for Sm in(S_sp,S2_sp):
   vk=mv_sp(Sm,p);nvk=(-vk[0],-vk[1],-vk[2],-vk[3],-vk[4],-vk[5])
   if(nvk if nvk<vk else vk)<p:return False
  return True
 pt_set=set(pts);seen=set();orbits=[]
 for p in pts:
  if p in seen:continue
  q=[p];orb={p};seen.add(p)
  while q:
   u=q.pop()
   for Msp in U_sps+(S_sp,):
    w=mv_sp(Msp,u)
    for cnd in(w,(-w[0],-w[1],-w[2],-w[3],-w[4],-w[5])):
     if cnd in pt_set and cnd not in orb:orb.add(cnd);seen.add(cnd);q.append(cnd)
  orbits.append(orb)
 exact1=missed=dups=0
 for members in orbits:
  c=sum(1 for m in members if is_can_int(m))
  if c==1:exact1+=1
  elif c==0:missed+=1
  else:dups+=c-1
 rel_err=max(abs(G_sym[i][j]/float(G_SC)-T4["GA"][i][j])for i in range(6)for j in range(6))
 def _pmul_h_z(a,b):
  t=[0]*11
  for i in range(6):
   for j in range(6):t[i+j]+=a[i]*b[j]
  for k in range(10,5,-1):
   c=t[k]
   for j in range(6):t[k-6+j]-=c*H_COEFFS[j]
  return t[:6]
 def _tr_h_z(p):
  cols,cur=[],list(p)
  for j in range(6):
   cols.append(cur)
   if j<5:
    nxt=[0]+cur[:5];c=cur[5]
    if c:
     for k in range(6):nxt[k]-=c*H_COEFFS[k]
    cur=nxt
  return sum(cols[i][i]for i in range(6))
 if _OHM_CACHE is None:is_ideal_primitive_oh([1,0,0,0,0,0],[0,1,0,0,0,0])
 ohm=_OHM_CACHE
 ring_ok=all(_pmul_h_z(ZK_NUM[i],ZK_NUM[j])==[3*sum(ohm[i][j][k]*ZK_NUM[k][m]for k in range(6))for m in range(6)]for i in range(6)for j in range(6))
 zk3i=[[1,0,1,-1,-3,-3],[0,1,-1,0,-1,-1],[0,0,-1,1,1,3],[0,0,-1,-1,0,2],[0,0,1,0,1,-1],[0,0,0,1,1,2]]
 zy_ok=([[sum(zk3i[k][i]*ZK_NUM[k][j]for k in range(6))for j in range(6)]for i in range(6)]==[[3 if i==j else 0 for j in range(6)]for i in range(6)])
 disc_O=_d6i([[_tr_h_z(_pmul_h_z(ZK_NUM[i],ZK_NUM[j]))//9 for j in range(6)]for i in range(6)])
 def _pmp(a,p):
  r=[c%p for c in a]
  while len(r)>1 and r[-1]==0:r.pop()
  return r
 def _pmz(a,b):
  t=[0]*(len(a)+len(b)-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):t[i+j]+=x*y
  return t
 def _pg(a,b,p):
  a,b=_pmp(a,p),_pmp(b,p)
  while b!=[0]:
   r=a[:];inv=pow(b[-1],-1,p)
   while len(r)>=len(b)and r!=[0]:
    f=(r[-1]*inv)%p;sh=len(r)-len(b)
    for i in range(len(b)):r[sh+i]=(r[sh+i]-f*b[i])%p
    r=_pmp(r,p)
   a,b=b,r
  inv=pow(a[-1],-1,p);return[(c*inv)%p for c in a]
 g2=[1,0,1,1];p2z=_pmz(g2,g2)
 ded2=(_pmp(p2z,2)==_pmp(H_COEFFS,2)and _pg(g2,[(H_COEFFS[k]-p2z[k])//2 for k in range(7)],2)==[1])
 c19=[14,6,18,1];g19=[6,1];p19z=_pmz(_pmz(g19,g19),_pmz(g19,c19))
 ded19=(_pmp(p19z,19)==_pmp(H_COEFFS,19)and _pg(c19,[k*c19[k]for k in range(1,4)],19)==[1]and sum(c19[k]*pow(13,k,19)for k in range(4))%19!=0 and _pg(g19,[(H_COEFFS[k]-p19z[k])//19 for k in range(7)],19)==[1])
 max_ord_ok=(ring_ok and zy_ok and disc_O==-(2**9)*(19**2)and ded2 and ded19)
 from decimal import Decimal as D,Context,ROUND_FLOOR,ROUND_CEILING
 cd,cu=Context(prec=55,rounding=ROUND_FLOOR),Context(prec=55,rounding=ROUND_CEILING)
 da,ua,ds,us,dm,um,dd,ud=cd.add,cu.add,cd.subtract,cu.subtract,cd.multiply,cu.multiply,cd.divide,cu.divide
 Ic=lambda x:(D(x),D(x));Ia=lambda a,b:(da(a[0],b[0]),ua(a[1],b[1]));Is=lambda a,b:(ds(a[0],b[1]),us(a[1],b[0]))
 Im=lambda a,b:(min(dm(x,y)for x in a for y in b),max(um(x,y)for x in a for y in b))
 Iq=lambda a:(D(0),max(um(a[0],a[0]),um(a[1],a[1])))if a[0]<=0<=a[1]else((dm(a[0],a[0]),um(a[1],a[1]))if a[0]>0 else(dm(a[1].copy_negate(),a[1].copy_negate()),um(a[0].copy_negate(),a[0].copy_negate())))
 def Id(a,b):assert b[0]>0 or b[1]<0;return(min(dd(x,y)for x in a for y in b),max(ud(x,y)for x in a for y in b))
 def Ir(a):
  lo,hi=cd.sqrt(a[0]),cu.sqrt(a[1])
  while um(lo,lo)>a[0]:lo=cd.next_minus(lo)
  while dm(hi,hi)<a[1]:hi=cu.next_plus(hi)
  return(lo,hi)
 PT=D("1e-40");EL,EH=ds(D(1),PT),ua(D(1),PT)
 Il=lambda a:(ds(cd.ln(a[0]),PT),ua(cu.ln(a[1]),PT));Ie=lambda a:(dm(cd.exp(a[0]),EL),um(cu.exp(a[1]),EH))
 Ca=lambda a,b:(Ia(a[0],b[0]),Ia(a[1],b[1]));Cm=lambda a,b:(Is(Im(a[0],b[0]),Im(a[1],b[1])),Ia(Im(a[0],b[1]),Im(a[1],b[0])))
 HC=[D(c)for c in H_COEFFS];DHC=[D(k*H_COEFFS[k])for k in range(1,7)];DEL=D("1e-40");roots_I=[];rouche_ok=True
 for j in range(3):
  zr,zi=D(str(T4["EMB_ZK"][j][1][0])),D(str(T4["EMB_ZK"][j][1][1]))
  for _ in range(7):
   fr,fi=HC[6],D(0);dfr,dfi=DHC[5],D(0)
   for k in range(5,-1,-1):fr,fi=da(ds(dm(fr,zr),dm(fi,zi)),HC[k]),da(dm(fr,zi),dm(fi,zr))
   for k in range(4,-1,-1):dfr,dfi=da(ds(dm(dfr,zr),dm(dfi,zi)),DHC[k]),da(dm(dfr,zi),dm(dfi,zr))
   den=da(dm(dfr,dfr),dm(dfi,dfi))
   zr=ds(zr,dd(da(dm(fr,dfr),dm(fi,dfi)),den))
   zi=ds(zi,dd(ds(dm(fi,dfr),dm(fr,dfi)),den))
  zI=((zr,zr),(zi,zi));ck=[(Ic(0),Ic(0))]*7;ck[0]=(Ic(H_COEFFS[6]),Ic(0))
  for m in range(5,-1,-1):
   for k in range(6,0,-1):ck[k]=Ca(Cm(ck[k],zI),ck[k-1])
   ck[0]=Ca(Cm(ck[0],zI),(Ic(H_COEFFS[m]),Ic(0)))
  lhs=Ir(Ia(Iq(ck[0][0]),Iq(ck[0][1])))[1];dp=um(DEL,DEL)
  for k in range(2,7):lhs=ua(lhs,um(Ir(Ia(Iq(ck[k][0]),Iq(ck[k][1])))[1],dp));dp=um(dp,DEL)
  if not(lhs<dm(Ir(Ia(Iq(ck[1][0]),Iq(ck[1][1])))[0],DEL)):rouche_ok=False
  roots_I.append(((ds(zr,DEL),ua(zr,DEL)),(ds(zi,DEL),ua(zi,DEL))))
 sep_ok=all(r[1][0]>D("0.05")for r in roots_I)and all(max(ds(roots_I[i][0][0],roots_I[j][0][0]).copy_abs(),ds(roots_I[i][1][0],roots_I[j][1][0]).copy_abs())>D("0.05")for i in range(3)for j in range(i+1,3))
 EZ_I=[]
 for zI in roots_I:
  pw=[(Ic(1),Ic(0))]
  for _ in range(5):pw.append(Cm(pw[-1],zI))
  row=[]
  for k in range(6):
   sr,si=Ic(0),Ic(0)
   for m in range(6):
    if ZK_NUM[k][m]:cI=Ic(ZK_NUM[k][m]);sr=Ia(sr,Im(cI,pw[m][0]));si=Ia(si,Im(cI,pw[m][1]))
   row.append((Id(sr,Ic(3)),Id(si,Ic(3))))
  EZ_I.append(row)
 c_ga_I=Id(Ic(str(T4["GA"][0][0])),Ic(3));GA_I=[]
 for i in range(6):
  row=[]
  for j in range(6):
   acc=Ic(0)
   for k in range(3):acc=Ia(acc,Ia(Im(EZ_I[k][i][0],EZ_I[k][j][0]),Im(EZ_I[k][i][1],EZ_I[k][j][1])))
   row.append(Im(c_ga_I,acc))
  GA_I.append(row)
 lb_I=Ic("0.88");LI=[[Ic(1 if i==j else 0)for j in range(6)]for i in range(6)];DI=[Ic(0)]*6;chol_I=True
 for c in range(6):
  d=Is(GA_I[c][c],lb_I)
  for k in range(c):d=Is(d,Im(Iq(LI[c][k]),DI[k]))
  if d[0]<=0:chol_I=False;break
  DI[c]=d
  for r in range(c+1,6):
   s=GA_I[r][c]
   for k in range(c):s=Is(s,Im(Im(LI[r][k],LI[c][k]),DI[k]))
   LI[r][c]=Id(s,d)
 sq_a=Ic(0)
 for i in range(6):
  for j in range(6):sq_a=Ia(sq_a,Iq(Is(Id(Ic(G_sym[i][j]),Ic(G_SC)),GA_I[i][j])))
 eps_G_hi=Id(Ir(sq_a),lb_I)[1];ei_I=Ic("1e-4")
 un_I=[]
 for fi in range(6):
  row=[]
  for j in range(3):
   sr,si=Ic(0),Ic(0)
   for k in range(6):
    if UMATS[fi][k][0]:cI=Ic(UMATS[fi][k][0]);sr=Ia(sr,Im(cI,EZ_I[j][k][0]));si=Ia(si,Im(cI,EZ_I[j][k][1]))
   row.append(Ia(Iq(sr),Iq(si)))
  un_I.append(row)
 def _d3I(m):return Ia(Is(Im(m[0][0],Is(Im(m[1][1],m[2][2]),Im(m[1][2],m[2][1]))),Im(m[0][1],Is(Im(m[1][0],m[2][2]),Im(m[1][2],m[2][0])))),Im(m[0][2],Is(Im(m[1][0],m[2][1]),Im(m[1][1],m[2][0]))))
 verts_I=[]
 for i in range(6):
  for j in range(i+1,6):
   M3=[[Is(un_I[i][k],Ic(1))for k in range(3)],[Is(un_I[j][k],Ic(1))for k in range(3)],[Ic(1)]*3];D0=_d3I(M3)
   if D0[0]<=D("1e-20")and D0[1]>=D("-1e-20"):continue
   w=[]
   for col in range(3):
    Mc=[r[:]for r in M3];Mc[0][col]=Mc[1][col]=(ei_I[1].copy_negate(),ei_I[0].copy_negate());Mc[2][col]=Ic(1);w.append(Id(_d3I(Mc),D0))
   if all(x[0]>0 for x in w):
    feas=True
    for f in range(6):
     val=Ic(0)
     for k in range(3):val=Ia(val,Im(w[k],Is(un_I[f][k],Ic(1))))
     if val[1]<ds(ei_I[0].copy_negate(),D("1e-20")):feas=False;break
    if feas:verts_I.append(w)
 V_min_lo=min(min(x[0]for x in w)for w in verts_I);B_max_hi=ud(ua(D(1),ei_I[1]),V_min_lo);log_B_hi=ua(cu.ln(B_max_hi),PT)
 l1_I=[Il(x)for x in un_I[0]];l2_I=[Il(x)for x in un_I[2]];detL_I=Is(Im(l1_I[0],l2_I[1]),Im(l1_I[1],l2_I[0]))
 c1_hi=max(x.copy_abs()for x in Id(l2_I[1],detL_I));c2_hi=max(x.copy_abs()for x in Id(l2_I[0],detL_I))
 M_bnd=int(ceil(float(um(um(ua(c1_hi,c2_hi),D(2)),log_B_hi))))+2
 facets={(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)};inner_cnt=0;min_gap_lo=D("1e9")
 for m1 in range(-M_bnd,M_bnd+1):
  for m2 in range(-M_bnd,M_bnd+1):
   if(m1,m2)in facets:continue
   ms_I=[Ie(Ia(Im(Ic(m1),l1_I[j]),Im(Ic(m2),l2_I[j])))for j in range(3)]
   if max(x[0]for x in ms_I)>=B_max_hi:continue
   inner_cnt+=1
   for w in verts_I:
    g_I=Ic(0)
    for j in range(3):g_I=Ia(g_I,Im(w[j],Is(ms_I[j],Ic(1))))
    if g_I[0]<min_gap_lo:min_gap_lo=g_I[0]
 max_t2_hi=D(0)
 for s1 in(D("-0.5"),D("0.5")):
  for s2 in(D("-0.5"),D("0.5")):
   s_I=Ic(0)
   for j in range(3):s_I=Ia(s_I,Ie(Ia(Im((s1,s1),l1_I[j]),Im((s2,s2),l2_I[j]))))
   if s_I[1]>max_t2_hi:max_t2_hi=s_I[1]
 R2_u_box_hi=um(um(c_ga_I[1],max_t2_hi),D("1.01"))
 cert_I_ok=(rouche_ok and sep_ok and chol_I and len(verts_I)==6 and eps_G_hi<D("3.0e-5")
  and dm(dd(ds(D(1),eps_G_hi),ua(D(1),eps_G_hi)),da(D(1),ei_I[0]))>D("1.00004")
  and V_min_lo>D("0.09")and inner_cnt==3 and min_gap_lo>D("0.40"))
 T_u_box=int(um(R2_u_box_hi,D(G_SC)).to_integral_value(rounding=ROUND_CEILING))
 def _norm_oh(vec):
  poly=[sum(vec[k]*ZK_NUM[k][j]for k in range(6))for j in range(6)];cols,cur=[],poly[:]
  for j in range(6):
   cols.append(cur)
   if j<5:
    nxt=[0]+cur[:5];c=cur[5]
    if c:
     for k in range(6):nxt[k]-=c*H_COEFFS[k]
    cur=nxt
  d=abs(_d6i([[cols[c][r]for c in range(6)]for r in range(6)]));assert d%729==0;return d//729
 ubox_pts=[p for p in pts if qf_i(G_sym,p)<=T_u_box]
 ubox_units=[p for p in ubox_pts if _norm_oh(p)==1]
 unit_idx1_ok=(max_ord_ok and T_u_box<T_lim and set(ubox_units)=={(1,0,0,0,0,0),(-1,0,0,0,0,0)})
 return dict(ok=(exact1==len(orbits)and missed==0 and dups==0 and len(orbits)==382 and rel_err<1e-5 and max_ord_ok and cert_I_ok and unit_idx1_ok),
 zero_float_geom=True,sigma_inv_exact=True,max_order_dedekind_ok=max_ord_ok,disc_Oh=disc_O,exact_frac_enum_ok=True,certified_interval_ok=cert_I_ok,rouche_roots_ok=rouche_ok,
 G_sym_scale=G_SC,G_sym_rel_err=float(f"{rel_err:.2e}"),unit_index1_proved=unit_idx1_ok,unit_box_R2_hi=round(float(R2_u_box_hi),4),unit_box_pts=len(ubox_pts),
 unit_box_units=len(ubox_units),roots_of_unity_exact=len(ubox_units),global_hexagon_proof_ok=cert_I_ok,eps_G_hi=round(float(eps_G_hi),8),eps_inf=1.0e-4,
 V_min_lo=round(float(V_min_lo),6),B_max_hi=round(float(B_max_hi),4),M_box_bound=M_bnd,inner_nonfacet_units=inner_cnt,min_nonfacet_gap_lo=round(float(min_gap_lo),6),
 R2=R2,T_lim_int=T_lim,nonzero_pts=len(pts),total_orbits=len(orbits),exact1=exact1,missed_orbits=missed,duplicate_orbits=dups,)
def log2_li(log2x):
 x=log2x*LN2
 return log2x-log2(x)+log2(1.0+1.0/x+2.0/(x*x))
def log2_U_z6(lf,lg):
 uf=(2.0**log2_li(lf))/6.0
 ug=(2.0**log2_li(lg))/6.0
 deg2=(2.0**(0.5*lf)/(0.5*lf*LN2)+2.0**(0.5*lg)/(0.5*lg*LN2))/6.0
 return log2(uf+ug+deg2)
NTT3_PRIMES=(69*(1<<55)+1,57*(1<<55)+1,54*(1<<55)+1)
def ntt3_F(M_prod):
 m=ceil(max(M_prod,1.0))
 if m<=1:return 220.0
 if m>(1<<55):
  return 3.0*ntt3_F(ceil(m/2.0))+20.0*m
 k=max(2,ceil(log2(m)))
 n2=float(1<<k)
 c2=55.0*n2*k+220.0*n2+(NTT3_XFER*n2 if n2>FFT_LOCAL else 0.0)
 n3=3.0*float(1<<(k-2))
 if n3>=m and k>=3:
  c3=55.0*n3*(k-2)+364.0*n3+(NTT3_XFER*n3 if n3>FFT_LOCAL else 0.0)
  return min(c2,c3)
 return c2
def fft_F(M_prod):return ntt3_F(M_prod)
def mul_F(M_prod):return ntt3_F(M_prod)
KH_SPLIT_CAP=0.30
KH_F_CAP=0.08
KH_G_CAP=0.17
KH_FG_CAP=0.20
def pbits_all(lo,hi):
 return 1.000028*(2.0**hi-(2.0**lo if lo else 0))/LN2 if lo==0 else 1.000028*(2.0**hi)/LN2
def pbits(lo,hi,cap=None):
 c_use=KH_SPLIT_CAP if cap is None else cap
 cut=0.5*hi
 base=pbits_all(lo,cut)if(not lo or lo<cut)else 0.0
 return base+c_use*pbits_all(max(lo,cut)if lo else cut,hi)
_GB_CACHE={}
def gbatch_cost(pb,leaf=9,L=None,c_leaf=FFT_LEAF):
 SP=ceil(pb/64.0)
 L=L if L is not None else ceil(log2((SP+1)/leaf))
 key=(SP,leaf,L,c_leaf)
 val=_GB_CACHE.get(key)
 if val is not None:return val
 G=L+4;n=lambda j:leaf*(2**j)
 PT=sum((2**(L-j))*(mul_F(2*n(j-1))+FFT_COPY*2*n(j-1))for j in range(1,L+1))
 RT=sum((2**(L-j))*((5.0/3.0)*mul_F(n(j)+G)+2.0*FFT_COPY*(n(j-1)+G))for j in range(1,L+1))
 val=(PT+RT+FFT_ROOT*mul_F(2*(max(SP,n(L))+G)))/(2**L)+c_leaf
 _GB_CACHE[key]=val
 return val
def split_parts(R,lb):
 np_lo=1.25506*(2.0**(0.5*lb))/(0.5*lb*LN2)
 np_hi=1.25506*(2.0**lb)/(lb*LN2)
 np_=np_lo+KH_FG_CAP*max(0.0,np_hi-np_lo)
 h=ceil(log2(np_));D=h+1
 SP=2.0**ceil(log2(np_*(lb/64.0+1.0)))
 L=ceil(log2(SPLIT_REC*SPLIT_Q));G=L+2;C=SPLIT_C
 W=R/SPLIT_Q+1.0;B=SPLIT_REC*R*D;occ=2*NORM_BITS
 return dict(trees=B*(2*(C*(36*L+2*G)+3*FFT_B)+2**15),roots=6*C*(2*SP*W*D+136*B)+6*FFT_B*B,gcd=B*2.0**24,small=R*2.0**26,output=occ*R*2.0**20,
 ideal=occ*R*SPLIT_IDEAL,control=B*2.0**16,ptree=D*SP*C+2*FFT_B*(2.0**h),debit=128*(W+2),)
def _det_mat(M):
 n=len(M);A=[row[:]for row in M];det=1.0
 for col in range(n):
  piv=max(range(col,n),key=lambda r:abs(A[r][col]))
  if piv!=col:A[col],A[piv]=A[piv],A[col];det= -det
  pv=A[col][col];det*=pv
  for r in range(col+1,n):
   fac=A[r][col]/pv
   for j in range(col,n):A[r][j]-=fac*A[col][j]
 return det
def _peel_fixed_point(lf,lg,exc,nf=403.5,ng=468.2):
 lU=log2_U_z6(lf,lg)
 terms=[]
 for m,lb in((nf,lf),(ng,lg)):
  u=m/lb;steps=int(round((lb-30.0)*10))
  for i in range(steps):
   b=30.0+i*0.1;c=b+0.1
   k=kappa_rho(u,c/lb)
   terms.append((2.0**(lU-c)*6.0*k,k*log(c/b)))
 def run(e_val):
  def step(h):
   tot=0.0;omh=1.0-h
   for lam0,wt in terms:
    lam=lam0*e_val
    if lam<100.0:tot+=wt*exp(-lam*omh)
   return 1.0-exp(-tot)
  h=1.0
  for _ in range(100):
   hn=step(h)
   if abs(hn-h)<1e-10:return hn
   h=hn
  return h
 h_val=run(exc)
 lo,hi=0.50,1.00
 for _ in range(24):
  mid=0.5*(lo+hi)
  if run(mid)>0.75:lo=mid
  else:hi=mid
 return dict(h=h_val,exc_crit=hi,stalls=(h_val>0.75))
_RET_CACHE={}
_PEEL_CACHE={}
_DESC_CACHE={}
_PROF_CACHE={}
def _exact_b2_cov(W_k,Rb,g00):
 dw=W_k*sqrt(g00)/Rb;num=den=0.0;nr=120
 for ir in range(nr):
  r5=(ir+0.5)/nr;wt=5.0*(r5**4)/nr;h=sqrt(max(0.0,1.0-r5*r5))
  den+=wt*(32.0/35.0)*(h**7);acc=0.0
  for ph in range(4):
   z= -h-(ph/4.0)*dw
   while z<h:
    z0=max(-h,z);z1=min(h,z+dw)
    if z1>z0:
     zmin=0.0 if(z0<=0.0<=z1)else min(abs(z0),abs(z1))
     acc+=(z1-z0)*((max(0.0,1.0-r5*r5-zmin*zmin))**3)
    z+=dw
  num+=wt*(acc/4.0)
 return num/den
def compute_ledger(track="bn254-gt-dlp",use_db=True,c_tail_z6=None,cap_scale=1.0,strict_per_owner=True,swar_probe=None,srt_g_cap=None):
 D=load_data()
 PQ=D["pool_quad"];T4=D["t4_geom"];FILT=D["filtsim_raw"]
 lf,lg,exc=OP["lf"],OP["lg"],OP["exc"]
 p_min,p_cut,p_lat,lx2=OP["p_min"],OP["p_cut"],OP["p_lat"],OP["lx2"]
 assert p_lat==OP["p_mid"]and lx2==lg
 sl,sb,w_sb=OP["slack"],OP["slack_b"],OP["w_sb"]
 if c_tail_z6 is None:c_tail_z6=OP["c_tail_z6"]
 if swar_probe is None:swar_probe=0.0
 rkey=(lg,lx2,p_min)
 if rkey not in _RET_CACHE:
  _RET_CACHE[rkey]=(Retention(lg,lx2,prime_min=2**p_min,explicit_bits=min(p_min,14),step=0.25,max_bits=530,tilt=0.05,odd=1),Retention(530.0,lx2,prime_min=2**p_min,explicit_bits=min(p_min,14),step=0.25,max_bits=531,tilt=0.005,odd=1),)
 ret_obj,cost_obj=_RET_CACHE[rkey]
 prim_tail_corr=min(1.0,EXACT_PRIM_T4/(PQ["n_prim_in"]/float(PQ["n_in"])))
 lU=log2_U_z6(lf,lg)
 target_lr=lU+log2(exc)+OP["margin"]
 frac_reg=PQ["frac_reg"];frac_reg_tr=PQ["frac_reg_trunc"];raw_mult_lat=1.0/frac_reg;raw_mult=1.0/frac_reg_tr
 log2_V_sub0=PQ["log2_V0"]+log2(frac_reg)
 n_in=PQ["n_in"];bins=PQ["bins"]
 assert sum(x[2]for x in bins)==PQ["n_prim_in"]and sum(x[1]for x in PQ["nonprim_bins"])==PQ["n_nonprim_in"]and PQ["n_prim_in"]+PQ["n_nonprim_in"]==n_in
 ALPHA_F_ODD=ALPHA_F_RAW+2.0/3.0
 ALPHA_G_ODD=ALPHA_G_RAW+1.0/3.0
 CLASSES_2=((7.0/9.0,0.0,0.0),(1.0/9.0,0.0,3.0),(1.0/9.0,6.0,0.0))
 _rt=rho_tab_1024();_rtl=len(_rt)-1;_ret_at=ret_obj.at;_inv_lf=1.0/lf;_inv_lg=1.0/lg;_g_ln2=GAMMA/LN2
 _pbins=[(lnf0,lng0,cnt_prim)for lnf0,lng0,cnt_prim,_ in bins if cnt_prim]
 def _rdb(u,ln):
  if u<=1.0:return 1.0
  ff=u*1024.0;i0=int(ff);rv0=_rt[i0]+(ff-i0)*(_rt[i0+1]-_rt[i0])if i0<_rtl else rho(u)
  u1=u-1.0
  if u1<=1.0:rv1=1.0
  else:f1=u1*1024.0;i1=int(f1);rv1=_rt[i1]+(f1-i1)*(_rt[i1+1]-_rt[i1])if i1<_rtl else rho(u1)
  r=rv0-_g_ln2*rv1/ln;return r if r>1e-300 else 1e-300
 def eval_d(d):
  ys=ys2=nf_w=ng_w=rg_w=0.0;d48=48.0*d;d12=12.0*d;scale=(1.0-0.0014)*prim_tail_corr
  for lnf0,lng0,cnt_prim in _pbins:
   nf_raw=lnf0+d48;ng_raw=lng0+d12;ng_q=round(ng_raw*2.0)*0.5
   r0=_ret_at(ng_q,sl,sb);r3=_ret_at(ng_q-3.0,sl,sb);y=rg_num=0.0
   for w2,v2f,v2g in CLASSES_2:
    nf_c=nf_raw-v2f+ALPHA_F_ODD;ng_c=ng_raw-v2g+ALPHA_G_ODD;r_c=r3 if v2g else r0
    yc=w2*(_rdb(nf_c*_inv_lf,nf_c)*_rdb(ng_c*_inv_lg,ng_c)if use_db else rho(nf_c*_inv_lf)*rho(ng_c*_inv_lg))*r_c
    y+=yc;rg_num+=yc*r_c
   r_g_bin=(rg_num/y)if y>0.0 else r0;y*=scale
   wy=cnt_prim*y;ys+=wy;ys2+=cnt_prim*(y*y);nf_w+=wy*nf_raw;ng_w+=wy*ng_raw;rg_w+=wy*r_g_bin
  my=ys/n_in
  return log2_V_sub0+d12+log2(max(my,1e-300)),my,ys2/n_in,nf_w/ys-CONTENT_BITS,ng_w/ys,rg_w/ys
 lo,hi=-0.1,0.1;fl,fh=eval_d(lo)[0]-target_lr,eval_d(hi)[0]-target_lr
 for _ in range(9):
  xm=hi-fh*(hi-lo)/(fh-fl)if fh!=fl else 0.5*(lo+hi)
  fm=eval_d(xm)[0]-target_lr
  if fm<0.0:lo,fl=xm,fm
  else:hi,fh=xm,fm
  if abs(fm)<1e-10:break
 d_opt=xm
 lr,my,my2,nf_rel,ng_rel,r_g=eval_d(d_opt)
 se_y=sqrt(max(my2-my*my,0.0)/n_in)/max(my,1e-300)/LN2
 se_f=sqrt((1.0-frac_reg)/(PQ["n_total"]*frac_reg))/LN2
 se_bits=sqrt(se_y*se_y+se_f*se_f)
 log2_V=log2_V_sub0+12.0*d_opt
 cells=2.0**log2_V
 GA=T4["GA"];GB=T4["GB"];g00=GB[0][0];det_GA=_det_mat(GA);det_GB=_det_mat(GB)
 ratio_B_A=GB[0][0]/GA[0][0]
 scale_r=2.0**d_opt
 Ra=PQ["ra0"]*scale_r;Rb_GA=PQ["rb0_ga"]*scale_r;Rb=Rb_GA*sqrt(ratio_B_A)
 Ra_tr=Ra*sqrt(QA_MAX);Rb_tr=Rb*sqrt(QB_MAX)
 Na=(1.0/6.0)*(pi**3/6.0)*(Ra_tr**6)/sqrt(det_GA)
 Vb=0.5*(pi**3/6.0)*(Rb_tr**6)/sqrt(det_GB)
 log2_Na=log2(Na);log2_Vb=log2(Vb)
 S5=[[GB[r][c]-GB[r][0]*GB[0][c]/g00 for c in range(1,6)]for r in range(1,6)]
 det_S5=_det_mat(S5)
 L_max=2.0*Rb_tr/sqrt(g00)
 n5=((4.0*pi**2)/15.0)*(Rb_tr**5)/sqrt(det_S5)
 L_b_avg=Vb/n5
 B_a=B_A_SHELLS
 us=[QA_MAX*((j/B_a)**(1.0/3.0))for j in range(B_a+1)]
 Na_sh=Na/B_a
 assert Na_sh+1.0<(1<<47)
 Na_sh_pad=(1<<ceil(log2(Na_sh+1.0)))
 score_arena_bytes=(L_max+2.0)*(2.0*Na_sh_pad)+Na_sh*4.0
 score_arena_log2=log2(score_arena_bytes)
 vb_score_log2=log2(Vb*2.0)
 max_sp_sb_ll=max(((2.0**b_hi)/w for _,b_hi,w,m in B2_BANDS if m=="ll"),default=0.0)
 max_sp_sb_x4=max(((2.0**b_hi)/w for _,b_hi,w,m in B2_BANDS if m=="x4"),default=0.0)
 b2_1p_ram_log2=log2(max(8.0*(Na_sh+1.0)+32.0*Na_sh+8.0*max_sp_sb_ll,48.0*Na_sh+8.0*(max_sp_sb_x4+1.0))+8.25*Na+score_arena_bytes)
 b2_build_1p_ram_log2=log2(max(8.0*(Na_sh+1.0)+32.0*Na_sh+8.0*max_sp_sb_ll,48.0*Na_sh+8.0*(max_sp_sb_x4+1.0))+8.25*Na)
 cA_emb=GA[0][0]/3.0;w_max_abs=max(abs(complex(r,im))for r,im in T4["EMB_W"]);eps_g_pad=1.00003
 sup_nf_bits=6.0*log2(36.0)+24.0*log2(eps_g_pad*max(Rb_GA*Rb_GA*QB_MAX/cA_emb,1.5*Ra*Ra*QA_MAX/cA_emb)/3.0)-CONTENT_BITS
 sup_ng_bits=6.0*log2(eps_g_pad*max(U_BN254*w_max_abs*(Ra*Ra*QA_MAX/cA_emb),Rb_GA*Rb_GA*QB_MAX/cA_emb)/3.0)
 sb_row=L_b_avg/w_sb+1.0
 max_sb=0
 cr2_chg=ceil(((sb_row*w_sb)/L_b_avg)*1.005*1000.0)/1000.0
 prof=dict(cr2=(sb_row*w_sb)/L_b_avg,sb_row=sb_row,max_sb=max_sb,frac_reg_analytic=0.6586)
 c_gb=gbatch_cost(pbits(0,lf,KH_F_CAP),leaf=7)
 c_srt_g=gbatch_cost(pbits(0,p_min,KH_G_CAP if srt_g_cap is None else srt_g_cap),leaf=8,c_leaf=C_SRT_G_LEAF)
 cc=(c_gb,C_TDIV+C_BATCH)
 np_bins=PQ["nonprim_bins"]
 cof_sum=sf_sum=s0_sum=0.0;tot_p=0
 for lnf0,lng0,cnt_prim,cnt_all in bins:
  tot_p+=cnt_prim
  if not cnt_prim:continue
  nf_raw=lnf0+48.0*d_opt;ng_raw=lng0+12.0*d_opt
  ng_q=round(ng_raw*2.0)/2.0
  sf0=cost_obj.survivor_upper(ng_q,sl,0.25);sf3=cost_obj.survivor_upper(ng_q-3.0,sl,0.25)
  if sb is None:
   n0=round(ng_q/cost_obj.step);n3=round((ng_q-3.0)/cost_obj.step);a0=floor((sl+0.25)/cost_obj.step)
   s00=cost_obj.two(n0,a0,0)/cost_obj.denominator[n0]
   s03=cost_obj.two(n3,a0,0)/cost_obj.denominator[n3]
  else:
   s00=cost_obj.survivor_upper(ng_q,sl,0.25,sb,lg,KMAX);s03=cost_obj.survivor_upper(ng_q-3.0,sl,0.25,sb,lg,KMAX)
  cof_p=sf_p=s0_p=0.0
  for w2,v2f,v2g in CLASSES_2:
   sf_c=sf3 if v2g else sf0;s0_c=s03 if v2g else s00
   q_f_c=min(1.0,Q_F_SAFETY*rho((nf_raw-v2f+ALPHA_F_ODD)/lf))
   cof_p+=w2*(sf_c*(C_NORM_G+c_srt_g)+s0_c*(C_NORM_F+cc[0]+q_f_c*cc[1]))
   sf_p+=w2*sf_c;s0_p+=w2*s0_c
  cof_sum+=cnt_prim*cof_p;sf_sum+=cnt_prim*sf_p;s0_sum+=cnt_prim*s0_p
 np_sum_cnt=0;cof_np_sum=sf_np_sum=0.0
 for ng_red_q0,cnt_np in np_bins:
  np_sum_cnt+=cnt_np
  ng_eff=floor((ng_red_q0+12.0*d_opt)*2.0)/2.0
  sf_np=cost_obj.survivor_upper(max(20.0,ng_eff),sl,0.25)
  cof_np_sum+=cnt_np*(sf_np*C_NORM_G);sf_np_sum+=cnt_np*sf_np
 assert tot_p==PQ["n_prim_in"]and np_sum_cnt==PQ["n_nonprim_in"]and tot_p+np_sum_cnt==n_in
 cof_sum+=cof_np_sum;sf_sum+=sf_np_sum
 cof_b=cof_sum/n_in;cof_np_cell=cof_np_sum/n_in
 sf=sf_sum/n_in;s0=s0_sum/n_in
 c_hit_lo=C_HIT_LO if strict_per_owner else 2.65
 n1_tot=LI_PAD*(2.0**log2_li(p_cut)-2.0**log2_li(p_min))
 vb_sum=sum(Vb*((min(QB_MAX,1.0-us[j])/QB_MAX)**3)for j in range(B_a))
 n5_sum=sum(n5*((min(QB_MAX,1.0-us[j])/QB_MAX)**2.5)for j in range(B_a))
 rep_lo=_WKF*B1_PREF_COV*raw_mult*log(p_cut/p_min)*(c_hit_lo*1.004+0.06)+(n1_tot*Vb*C_B1_ROW+n1_tot*n5_sum*64.0)/cells
 rep_b0=0.0
 assert B2_BANDS[0][0]==p_cut and B2_BANDS[-1][1]==p_lat
 n2_tot=n2_occ=probe=occ_sb=hi_ap=ins_hi=byt_hi=0.0
 for b_lo,b_hi,w_k,mode_k in B2_BANDS:
  bt=LI_PAD*(2.0**log2_li(b_hi)-2.0**log2_li(b_lo))
  cov_w=_exact_b2_cov(w_k,Rb_tr,g00)*raw_mult
  for j in range(B_a):
   sbf=sqrt(min(QB_MAX,1.0-us[j])/QB_MAX);L_b_j=L_b_avg*sbf;n5_j=n5*(sbf**5);sb_k=L_b_j/float(w_k)+1.0
   occ_sb+=(bt*n5_j*(sb_k*C_SB_ENTER_X4 if mode_k=="x4" else(sb_k*C_SB_ENTER_LL_SLOPE+C_SB_ENTER_LL_CONST)))/cells
  b_in=min(L_b_avg,2.0*w_k);i_in=max(0.0,L_b_avg-2.0*w_k);oob=1.005*w_k
  c_int=C_B2_INT_X4 if mode_k=="x4" else(C_B2_INT_LL+C_B2_PREF_EX)
  rat=min(1.0,(2.0*w_k)/L_max);f_sh=1.0-(1.0-rat*rat)**2.5
  c_bh=(1.0-f_sh)*C_B2_BND1_HIT+f_sh*C_B2_BND_HIT+C_B2_PREF_EX
  c_bm=(1.0-f_sh)*C_B2_BND1_MISS+f_sh*C_B2_BND_MISS+C_B2_PREF_EX
  eff_h=(i_in*c_int+b_in*c_bh+oob*c_bm)/L_b_avg
  hi_ap+=cov_w*log(b_hi/b_lo)*eff_h
  n2_tot+=bt;n2_occ+=bt
  sp_k=LI_PAD*0.5*((2.0**(2.0*b_hi))/(2.0*b_hi*LN2-1.0)-(2.0**(2.0*b_lo))/(2.0*b_lo*LN2-1.0))
  c_lst=C_LIST_B2 if mode_k=="x4" else C_LIST_LL
  c_lbk=C_LBUCKET_DN if mode_k=="x4" else C_LBUCKET_LL
  ins_hi+=0.5*bt*Na*c_lst+B_a*c_lbk*(sp_k/w_k)
  byt_hi+=8.0*(Na+B_a)*(0.5*bt)+B_a*8.0*(sp_k/w_k)
 n2_ent=n2_occ
 lat_hits=log(lx2/p_lat)
 c_lat_use=C_LAT_LEAF+(swar_probe if swar_probe else 0.0)
 lat_b2=lat_hits*c_lat_use*raw_mult_lat+C_LAT_FALLBACK_RES
 n_orb_lat=KH_SPLIT_CAP*(2.0**log2_li(lx2)-2.0**log2_li(p_lat))/3.0
 lat_lll=(n_orb_lat*(C_LAT_LLL+C_LAT_EXC))/cells
 scan_c=1.05*raw_mult*(C_SCAN_G8+(SCAN_G-1)*(C_SCAN_ZERO-C_SCAN_G8)/L_b_avg)+1024.0*sf*1.05*raw_mult
 per_cell=rep_lo+probe+occ_sb+hi_ap+lat_b2+lat_lll+scan_c+cof_b
 coll_u=log2(cells*per_cell*OP["work_allow"]/UNIT)
 n_lo=LI_PAD*0.5*(2.0**log2_li(p_cut)-2.0**log2_li(p_min))
 sp_lo=LI_PAD*0.5*(2.0**(2.0*p_cut))/(2.0*p_cut*LN2-1.0)
 ins_lo=(n_lo+1.0)*(C_LIST_B1*Na*_WBR)+C_LBUCKET*sp_lo*(B_a*K_B1_PREF)*_WBL
 byt_lo=8.0*Na*n_lo*_WBR+8.0*sp_lo*(B_a*K_B1_PREF)*_WBL
 avg_b2_list_ir=C_LIST_LL
 atab_entries=float(79151674016786875)
 btab_entries=float((3<<59)+435563320543855875)
 Na_plus=Na*(1.0+SC_GUARD_DELTA)**3
 ins_pre=Na_plus*(1000.0+12.0*log2(Na_plus))+2.0*atab_entries+2.0*btab_entries+50.0*(6.0*Na_plus+8.0*Vb)+20.0*Na+6.4*float(1<<36)*Na+128.0*Vb+2**36*(4.4*Vb+8.4*n5_sum)
 lists_u=log2((ins_lo+ins_hi+ins_pre)/UNIT)
 import stor_bound as SBND
 _sb=SBND.main()
 assert _sb["ok"]
 live_score_total_bytes=2.0**_sb["score_rows_total_log2"]
 cs_score_max_bytes=2.0**_sb["score_node_log2"]
 cs_blocks_max=int(ceil(cs_score_max_bytes/SC_BLOCK_BYTES))
 assert cs_blocks_max<(1<<26)
 score_1node_ram_log2=_sb["score_node_log2"]
 phase1_active_bytes=2.0**_sb["emitter_log2"]
 emit_1own_ram_log2=_sb["emitter_log2"]
 _b1r=18.0*Na+256.0*(2.0**p_cut)
 SB1=2.0**_sb["SB1_log2"];assert 7.0*SB1<=0.2*Vb
 list_1node_ram_log2=log2(SHARD_BYTES+max(2.0**b2_build_1p_ram_log2,float(1<<56))+N_SC_NODES*16400.0+1064.0*SB1)
 assert _b1r<float(1<<56)
 _top=64.0*Na+65536.0*(2.0**p_cut)+8.0*(2.0**p_lat)
 assert (byt_lo*8/_WBR+byt_hi)/2**36+_top<SHARD_BYTES
 assert list_1node_ram_log2<=_sb["list_node_log2"]
 n_logp_bands=int(round(16.0*(lx2-p_lat)))
 ram_log2=_sb["machine_total_log2"]
 peak_global_ram_cap_log2=ram_log2
 stor_bound_out=_sb
 max_a_box_emb=sum(510.0*max(abs(complex(T4["EMB_ZK"][j][k][0],T4["EMB_ZK"][j][k][1]))for j in range(3))for k in range(6))
 max_b_box_emb=sum((1022.0 if k==1 else 510.0)*max(abs(complex(T4["EMB_ZK"][j][k][0],T4["EMB_ZK"][j][k][1]))for j in range(3))for k in range(6))
 sup_ng_box_bits=6.0*log2(U_BN254*w_max_abs*(max_a_box_emb**2)+max_b_box_emb**2)
 max_pf_lat=int(floor(sup_ng_box_bits/p_lat))
 delta5=sqrt(max(sum((0.5 if(m>>i)&1 else-0.5)*S5[i][j]*(0.5 if(m>>j)&1 else-0.5)for i in range(5)for j in range(5))for m in range(32)))
 n5_disc_ub=2.0*n5*((1.0015+delta5/Rb)**5)
 def _k5_row(x):
  if len(x)>=12:
   p1,h2,h3,h4,h5p,ef=x[6],x[7],x[8],x[9],x[10],x[11]
   assert abs(p1+(h2+h3+h4+h5p)/ef-1.0)<5e-4
   return p1*C_UNIT_5+(h2*2.0*C_UNIT_5+h3*3.0*C_UNIT_5+h4*4.0*C_UNIT_5+h5p*C_GEN_1)/ef
  if len(x)==8:return x[7]
  return x[6]*C_UNIT_5+(1.0-x[6])*C_GEN_1
 def stats_of(key):
  rs=FILT[key];n=len(rs)
  def ms(idx):
   vs=[x[idx]for x in rs];m=sum(vs)/n
   s2=sum((v-m)**2 for v in vs)/(n-1)if n>1 else 0.0
   return m,sqrt(s2/n),max(vs)
  sw,sw_se,sw_mx=ms(1);nl,nl_se,nl_mx=ms(2);k5,_,k5_mx=ms(3);tr,_,tr_mx=ms(4);uf,_,_=ms(5);eu,_,_=ms(6)
  k5_exp_vs=[_k5_row(x)for x in rs]
  k5_exp_mean=sum(k5_exp_vs)/n
  min_nc=min(x[14]for x in rs)if len(rs[0])>=15 else 0
  return dict(n=n,N_sw_U=sw,N_sw_se=sw_se,N_sw_max=sw_mx,N_l_U=nl,N_l_se=nl_se,N_l_max=nl_mx,k5_ev=k5,k5_exp=k5_exp_mean,k5_exp_max=max(k5_exp_vs),tab_r=tr,tab_r_max=tr_mx,unit_frac=uf,exp_unit=eu,min_ncols=min_nc)
 fpref="m184_z6_61.0_61.0_e0.76_W160_"
 cur20=stats_of(f"{fpref}20");cur21=stats_of(f"{fpref}21");cur22=stats_of(f"{fpref}22")
 cur23=stats_of(f"{fpref}23");cur24=stats_of(f"{fpref}24");cur25=stats_of(f"{fpref}25");cur26=stats_of(f"{fpref}26")
 assert all(c["n"]>=4 for c in(cur20,cur21,cur22,cur23,cur24,cur25,cur26))
 nu_sw_cap=NU_SW_CAP*cap_scale
 nu_l_cap=NU_L_CAP*cap_scale
 k5_cap=K5_CAP;tab_r_cap=TAB_R_CAP
 pf_rl_cap=min(1.0,nu_l_cap/nu_sw_cap)
 log2_Ns_cap=lU+log2(nu_sw_cap);log2_N_cap=lU+log2(nu_l_cap)
 n_sm=5
 c_prod=(W_ROW*k5_cap*tab_r_cap+c_tail_z6)/pf_rl_cap
 def bw_steps(lN,m,n):return(1.0+BW_SIGMA)*(2.0**lN+n_sm)*(1.0/m+1.0/n)+6.0
 def gen_at_cap(lN_cap,m,n):
  L_cap=bw_steps(lN_cap,m,n);k=m+n
  return log2(4.0*(k**3)*log2(L_cap)*(fft_F(18.0*L_cap)+640.0*L_cap)/UNIT)
 best_la=None
 for m in(256,512,768,1024):
  for n in(24,32,40,48):
   N=2.0**log2_N_cap;L_st=bw_steps(log2_N_cap,m,n);s=N*(1.0/m+1.0/n)
   kap_bw=L_st/s;hor=max(s,(1.0+BW_SIGMA)*(N+n_sm)/n+4.0)/s
   per_n2=(c_prod+C_PROJ)*(n/m+1.0)*kap_bw+(1.0/m+1.0/n)*(c_prod+(n_sm*C_FULL+(n-n_sm)*C_Z64)+C_RED9)*hor
   per_n1=(m+n)*C_RED9
   la_u=log2((N*N*per_n2+N*per_n1+PREP_OWNER+PREP_GUARD+PREP_JOINT+LA_RESERVE+2.0**20)/UNIT)
   gen_u=gen_at_cap(log2_N_cap,m,n)
   both=log2(2.0**la_u+2.0**gen_u)
   if best_la is None or both<best_la["both"]:
    ncap=2.0**log2_Ns_cap;nslots=ceil(n/8.0)
    _sm=ncap*(256.0*nslots+256.0*nslots+W_ROW*8.0)+(2.0**55.04)*nslots+65.0*(2.0**53)
    _pm=(3.0*(m+n)*m+m*n)*L_st*32.0
    mem_1own=log2(_sm+_pm)
    best_la=dict(la=la_u,gen=gen_u,both=both,m=m,n=n,mem_1own=mem_1own,spmv_mem_log2=log2(_sm),pm_basis_mem_log2=log2(_pm),log2_N_cap=log2_N_cap,log2_Ns_cap=log2_Ns_cap,c_prod=c_prod)
 pkey_peel=(lf,lg,exc)
 if pkey_peel not in _PEEL_CACHE:
  _PEEL_CACHE[pkey_peel]=_peel_fixed_point(lf,lg,exc,nf=nf_rel,ng=ng_rel)
 peel=_PEEL_CACHE[pkey_peel]
 rels_low=lr-2.25*se_bits
 dkey=(round(rels_low,3),round(peel["h"],4))
 if dkey not in _DESC_CACHE:
  _DESC_CACHE[dkey]=DESC.evaluate_descent(D["descent"],rho,fft_F,gbatch_cost,pbits,log2_rels=rels_low,peel_h=peel["h"],nf_main=nf_rel,ng_main=ng_rel,)
 desc=_DESC_CACHE[dkey]
 rels_up=lr+2.25*se_bits;u_log=log2(UNIT);lb_max=max(lf,lg)
 canon_u=log2_Na+42.8-u_log
 r_bits=R_BN254.bit_length()
 bit_f=48*48*512+7*4*36*512+2**17+3*48*1024+2*48*6*8*40+48*1024
 bit_g=12*12*512+1*6*6*512+2**15+3*12*1024+2*12*6*8*40+12*1024
 sm_u=rels_up+log2(48*r_bits*bit_f+12*r_bits*bit_g)-u_log
 misc_u=log2(2.0**53+3.0*N_SC_NODES*Na_plus+2.0**(rels_up+30.4)+2.0**(rels_up+14.0)+n_logp_bands*N_SC_NODES*4130.0+float(N_LIST_STREAM_KINDS)*float(1<<36)*N_SC_NODES*4130.0+n1_tot*N_SC_NODES*4130.0)-u_log
 pi_b=log2(1.25506*(2.0**lb_max)/(lb_max*LN2))
 fb_u=pi_b+23.4-u_log
 split_u=log2(sum(split_parts(2.0**rels_up,lb=lb_max).values()))-u_log+1.0
 L_gt=ceil(lb_max-6.0+0.53)
 pg=sum((2.0**(lb_max-6.0+0.53-j-1))*(fft_F(2.0**(j+1))+3.0*(2.0**(j+1)))for j in range(L_gt))
 gtree_u=log2(2.0*pg)-u_log
 rows=dict(collection=coll_u,chi_setup=65.712288,la=best_la["la"],gen=best_la["gen"],lists=lists_u,two_logs=desc["two_log"],canon=canon_u,sm=sm_u,misc=misc_u,fb=fb_u,split=split_u,gtree=gtree_u,
 cert=50.000001,verify=11.000001,wrapper=53.712288,)
 if track=="bn254-pairing-inversion":
  rows["output"]=log2(G2_OUTPUT_BOUND/UNIT)
 exp_run=log2(sum(2.0**v for v in rows.values()))
 ncap=ceil(2.0**best_la["log2_Ns_cap"])
 unified_meter_fail=(2.0**exp_run)/(2.0**T_RUN_LOG2)
 fail_terms=dict(main_matrix=4.0*ncap/(2.0**64),filt_cap_h3_given_H3=P_CAP_H3,split_abort=(2.0**split_u)/(2.0**SPLIT_CAP-(2.0**SPLIT_RES)/UNIT),descent_wald_retry=0.0,descent_node_abort=0.0,descent_repeat_ideal=0.0,unified_run_meter=unified_meter_fail,storage_cap_abort=STORAGE_CAP_ABORT_PROB,
 miller_rabin_64r=(2.0**69)*(2.0**(-128.0)),group_order_boot_and_sm=(8.0+2.0*desc["e_trials_per_log"])/R_BN254,numerical_pure_int=0.0,)
 succ_exact=1.0-sum(fail_terms.values())
 succ_decl=f"{floor(succ_exact*1e6)/1e6:.6f}"
 s_frac=F(succ_decl)
 K_gate=int(-(-(F(39,100)*(2**64))//s_frac))
 gate_passes=(F(K_gate,2**64)*s_frac>=F(39,100))
 gate_prev_fails=(F(K_gate-1,2**64)*s_frac<F(39,100))
 exp_run_decl=ceil(exp_run*1e9)/1e9
 claim_exact=log2(1.0+float(F(K_gate,2**64))*(2.0**exp_run_decl))
 claim_num=ceil(claim_exact*1e9)
 claim_decl=f"{claim_num/1e9:.9f}"
 all4_mono=(cur20["N_sw_U"]>cur21["N_sw_U"]>cur22["N_sw_U"]>cur23["N_sw_U"]>cur24["N_sw_U"]>cur25["N_sw_U"]>cur26["N_sw_U"]and
 cur20["N_l_U"]>cur21["N_l_U"]>cur22["N_l_U"]>cur23["N_l_U"]>cur24["N_l_U"]>cur25["N_l_U"]>cur26["N_l_U"]and
 cur20["k5_exp"]>cur21["k5_exp"]>cur22["k5_exp"]>cur23["k5_exp"]>cur24["k5_exp"]>cur25["k5_exp"]>cur26["k5_exp"]
 )
 all_runs_26_below_4caps=(cur26["N_sw_max"]<nu_sw_cap and
 cur26["N_l_max"]<nu_l_cap and
 max(cur20["k5_exp_max"],cur21["k5_exp_max"],cur22["k5_exp_max"],cur23["k5_exp_max"],cur24["k5_exp_max"],cur25["k5_exp_max"],cur26["k5_exp_max"])<k5_cap and
 max(cur20["tab_r_max"],cur21["tab_r_max"],cur22["tab_r_max"],cur23["tab_r_max"],cur24["tab_r_max"],cur25["tab_r_max"],cur26["tab_r_max"])<tab_r_cap
 )
 return dict(track=track,use_db=use_db,c_tail_z6=c_tail_z6,cap_scale=cap_scale,strict_per_owner=strict_per_owner,swar_probe=swar_probe,op=OP,
 pool=dict(n_total=PQ["n_total"],n_in=n_in,n_prim_in=PQ["n_prim_in"],frac_reg=frac_reg,frac_reg_analytic=prof["frac_reg_analytic"]),geometry=dict(log2_U=lU,log2_V=log2_V,se_bits=se_bits,log2_Vb=log2_Vb,log2_Na=log2_Na,log2_Na_shell=log2(Na_sh),log2_N5=log2(n5),L_b_avg=L_b_avg,L_max=L_max,Ra=round(Ra,2),Rb_GA=round(Rb_GA,2),sb_row=prof["sb_row"],
 max_sb=prof["max_sb"],cr2_tight=cr2_chg,p_occ_b2=n2_occ/max(n2_tot,1.0),p_ent_b2=n2_ent/max(n2_tot,1.0),nf_rel=nf_rel,ng_rel=ng_rel,r_g=r_g,sf=sf,s0=s0,score_arena_log2=score_arena_log2,vb_score_log2=vb_score_log2,
 b0_blocks=B_a,n_sc_nodes_exp=N_SC_EXP,sc_block_exp=SC_BLOCK_EXP,sc_w_subband=SC_W_SUBBAND,sc_guard_delta=SC_GUARD_DELTA,stor_bound=stor_bound_out,live_score_total_log2=round(log2(live_score_total_bytes),4),cs_blocks_max=cs_blocks_max,lat_8b_direct_ptr=True,delta5_blichfeldt=round(delta5,6),n5_disc_ub_log2=round(log2(n5_disc_ub),4),list_1node_ram_log2=list_1node_ram_log2,b1_lists_log2=round(log2(byt_lo/2**36),4),lists_8x_log2=round(log2((byt_lo*8/_WBR+byt_hi)/2**36),4),lists_node_max_log2=round(log2((byt_lo*8/_WBR+byt_hi)/2**36+_top),4),b1_build_scratch_log2=round(log2(_b1r),4),emit_1own_ram_log2=emit_1own_ram_log2,score_1node_ram_log2=score_1node_ram_log2,det_max_pf_lat=max_pf_lat,b2_1p_ram_log2=b2_1p_ram_log2,b2_build_1p_ram_log2=b2_build_1p_ram_log2,avg_b2_list_ir=round(avg_b2_list_ir,3),ram_log2=ram_log2,peak_global_ram_cap_log2=peak_global_ram_cap_log2,sup_nf_bits=round(sup_nf_bits,2),sup_ng_bits=round(sup_ng_bits,2)),
 per_cell=dict(total=per_cell,rep_b0=rep_b0,rep_lo=rep_lo,probe=probe,occ_sb=occ_sb,hi_ap=hi_ap,lat_b2=lat_b2,lat_lll=lat_lll,scan_zero=scan_c,cof=cof_b,cof_nonprim=round(cof_np_cell,6),
 n_prim_exact=tot_p,n_nonprim_exact=np_sum_cnt,nonprim_bins_n=len(np_bins)),la=dict(bw=[best_la["m"],best_la["n"]],log2_N_cap=best_la["log2_N_cap"],log2_Ns_cap=best_la["log2_Ns_cap"],
 W_row=W_ROW,nu_sw_cap=nu_sw_cap,nu_l_cap=nu_l_cap,k5_cap=k5_cap,tab_r_cap=tab_r_cap,pf_rl_cap=pf_rl_cap,c_prod=c_prod,spmv_mem_log2=round(best_la["spmv_mem_log2"],4),pm_basis_mem_log2=round(best_la["pm_basis_mem_log2"],4),mem_1own_log2=best_la["mem_1own"],fits_1own=(best_la["mem_1own"]<64.0),h3_mode="fixed_instance_deterministic_4cap",
 ladder_all4_monotonic=all4_mono,all_runs_26_below_4caps=all_runs_26_below_4caps,filtsim_cur20=cur20,filtsim_cur21=cur21,filtsim_cur22=cur22,filtsim_cur23=cur23,filtsim_cur24=cur24,filtsim_cur25=cur25,filtsim_cur26=cur26,peel=peel),
 descent=desc,rows=rows,exp_run=exp_run,exp_run_decl=f"{exp_run_decl:.9f}",fail_terms=fail_terms,succ_exact=succ_exact,succ_decl=succ_decl,K_gate=K_gate,gate_passes=gate_passes,gate_prev_fails=gate_prev_fails,claim_exact=claim_exact,claim_decl=claim_decl,score_numerator=str(claim_num),)
_OH_U=((1,-1,-1,-1,1,0),(-1,1,1,1,0,0),(1,1,0,1,-1,0),(0,0,1,0,0,1),(-1,0,0,0,0,1),(1,0,0,-1,0,-1),(0,-1,0,-1,1,0),(1,0,0,0,1,-1),(0,1,1,1,0,1),(-3,0,2,1,0,1),(0,1,0,0,0,0),(1,0,-1,0,0,0),(-1,2,2,1,0,1),(-1,0,0,1,0,2),(-1,-1,-1,0,-1,1))
_OHM_CACHE=None;_OHM_T=None
def _d6i(m):
 A=[r[:]for r in m];pv=1;sg=1
 for k in range(6):
  if A[k][k]==0:
   for r in range(k+1,6):
    if A[r][k]:A[k],A[r]=A[r],A[k];sg=-sg;break
   else:return 0
  p1=A[k][k];ak=A[k]
  for i in range(k+1,6):
   ai=A[i];aik=ai[k]
   for j in range(k+1,6):ai[j]=(ai[j]*p1-aik*ak[j])//pv
  pv=p1
 return sg*A[5][5]
def is_ideal_primitive_oh(a,b):
 if math.gcd(*a,*b)>1:return False
 global _OHM_CACHE,_OHM_T
 if _OHM_T is None:
  if _OHM_CACHE is None:
   m=[[[1 if c==j else 0 for c in range(6)]for j in range(6)]for _ in range(6)];ix=0
   for i in range(1,6):
    m[i][0]=m[0][i]
    for j in range(i,6):m[i][j]=m[j][i]=_OH_U[ix];ix+=1
   _OHM_CACHE=m
  m=_OHM_CACHE;_OHM_T=[[[m[i][j][c]for i in range(6)]for c in range(6)]for j in range(6)]
 M=_OHM_T;_sp=math.sumprod
 Ra=[[_sp(a,M[j][c])for c in range(6)]for j in range(6)]
 Rb=[[_sp(b,M[j][c])for c in range(6)]for j in range(6)]
 g=math.gcd(_d6i(Ra),_d6i(Rb))
 if g==1:return True
 R=([[x%g for x in r]for r in Ra+Rb]+[[g if c==i else 0 for c in range(6)]for i in range(6)])if g else Ra+Rb
 nr=len(R);r0=0
 for col in range(6):
  while True:
   nz=[i for i in range(r0,nr)if R[i][col]]
   if not nz:return False
   if len(nz)==1:
    p=nz[0]
    if p!=r0:R[r0],R[p]=R[p],R[r0]
    if abs(R[r0][col])!=1:return False
    r0+=1;break
   p=min(nz,key=lambda i:abs(R[i][col]))
   if R[p][col]<0:R[p]=[-x for x in R[p]]
   pv=R[p][col];pr=R[p]
   for i in nz:
    if i==p:continue
    ri=R[i];q=ri[col]//pv
    for c in range(col,6):ri[c]-=q*pr[c]
 return True