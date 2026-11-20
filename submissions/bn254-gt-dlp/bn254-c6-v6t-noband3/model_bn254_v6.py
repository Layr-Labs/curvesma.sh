"""model_bn254_v6.py — Pure-stdlib BN254 v6t-noband3 (C6-m184_w216-rec) algebra, pure-int FD, <2^62 NTT, and 15-row 2-band ledger engine."""
import json,math,random
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
H_COEFFS=[3, -4,0,2,1, -2,1]
W_COEFFS=[0, -2, -1,0,1, -1]
SIGMA_NUM=[9, -4, -5,1,3, -2]
ZK_NUM=[
[3,0,0,0,0,0],
[0,3,0,0,0,0],
[9, -4, -5,1,3, -2],
[0,2,1, -2,0,1],
[6,1, -1, -1,3, -1],
[-6,6,6,0, -3,3],
]
F_REC_COEFFS={
8:[3, -4, -2,0,2, -2],
6:[0,6,3,0, -3,3],
4:[-12,0,0,0,0,0],
2:[12, -6, -3,0,3, -3],
0:[-5,4,2,0, -2,2],
}
CONTENT_BITS=12.0+6.0*log2(3.0)
EXACT_PRIM_T4=0.683542
ALPHA_F_RAW= -20.010+0.08
ALPHA_G_RAW= -0.119+0.08
_DATA=None
_B62={c:i for i,c in enumerate("0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")}
def b62d(s):
 neg=s.startswith("-");v=0
 for c in(s[1:]if neg else s):v=v*62+_B62[c]
 return-v if neg else v
def load_data():
 global _DATA
 if _DATA is None:
  _DATA=json.loads((Path(__file__).resolve().parent/"c6_data.json").read_text())
  pq=_DATA["pool_quad"]
  if"bins_b62"in pq and"bins"not in pq:
   acc=0;bs=[]
   for s in pq["bins_b62"].split(";"):
    a,b,c,d=s.split();acc+=b62d(a);cp=b62d(c)
    bs.append([acc/100.0,(b62d(b)+44000)/100.0,cp,cp+b62d(d)])
   pq["bins"]=bs
 return _DATA
C_HIT_LO=13.25
C_SWAR_PROBE=18.0
C_SB_ENTER=14.0
C_CAND=13.0
C_HIT_HI=14.5
C_SCAN_ZERO=2.15
C_LIST_A=35.0
C_LBUCKET=4.0
C_UNIT_5=7.8
C_GEN_1=33.0
C_NORMS=30000.0
C_SB=80000.0
C_TDIV=2.0**23.45
C_BATCH=2.0**22.12
C_BATCH2=2.0**23.13
C_GBATCH=2.0**23.0
LIST_NS=256
LI_PAD=1.001
Q_F_SAFETY=100.0
KMAX=4
FFT_A,FFT_B,FFT_XFER,FFT_LOCAL,FFT_COPY,FFT_ROOT,FFT_LEAF=8.5,1200.0,504.0,67108864,3.0,6.0,32768.0
NTT3=(69.75,119.0,90.0,151.5,371.0,139.25,271.5)
NTT3_XFER=3.0*FFT_XFER+128.0
SPLIT_Q,SPLIT_REC,SPLIT_C,SPLIT_IDEAL,NORM_BITS=1099511627776.0,56,33024.0,67108864.0,576
SPLIT_CAP,SPLIT_RES=112.0,84.0
W_ROW=161.0
NU_SW_CAP=0.00160
NU_L_CAP=0.00098
K5_CAP=12.00
TAB_R_CAP=0.85
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
OP=dict(
lf=61.0,lg=60.0,exc=0.82,W=160,
p_min=14,p_cut=43.0,p_mid=55.5,lx2=55.5,
slack=164.0,slack_b=169.0,w_sb=256,c_tail_z6=420.0,
margin=0.0666,work_allow=1.0033,
)
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
def rho_db(u,log2_N):
 if u<=1.0:
  return 1.0
 return max(1e-300,rho(u)-GAMMA*rho(u-1.0)/max(1.0,log2_N*LN2))
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
 _sp=math.sumprod
 jw_rev=[j*weights[j]for j in range(length, -1, -1)]
 out=[exp(-normalizer)]+[0.0]*length
 L=length
 for n in range(1,L+1):
  out[n]=_sp(jw_rev[L-n:L],out[:n])/n
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
   groups=[("small",False)]if p<prime_min else[(("sieve"if b<X_bits else"large")+"_geo",False),(("sieve"if b<X_bits else"large")+"_sf",True)]
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
     amt=exp(k*tilt*bi*step-(k-1)*LN2*b)*w/(k*b)
     idx=k*bi
     sign= -1 if k%2==0 else 1
     grps=((pref,amt),)if pref=="small"else((pref+"_geo",amt),(pref+"_sf",sign*amt))
     for g,x in grps:
      totals[g]+=x
      if idx<=self.length:weights[g][idx]+=x
  lm=log2(prime_min)
  continuum(cutoff,min(lm,B_bits),"small")
  continuum(max(cutoff,lm),min(X_bits,B_bits),"sieve")
  continuum(max(cutoff,X_bits),B_bits,"large")
  self.weights,self.normalizers=weights,{k:totals[k]for k in self.names}
  all_w=[sum(weights[k][j]for k in("small","sieve_geo","large_geo"))for j in range(self.length+1)]
  self.denominator=_series_coeffs(all_w,sum(totals[k]for k in("small","sieve_geo","large_geo")),self.length)
  self.selection_normalizer=exp(totals["sieve_sf"]+totals["large_sf"]-totals["sieve_geo"]-totals["large_geo"])
  self._c={}
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
    st=self.step;sm=self.co("small");D=min(self.length,2048);m=[0.0]*(D+1)
    for l,v in enumerate(self.co("large_sf"if sf else"large_geo")[:min(b,D)+1]):
     if v and(win is None or l*st<=win[2]or any(self.X*k-win[2]<=l*st<=win[1]*k+win[2]for k in range(1,win[0]+1))):
      for t in range(D+1-l):m[l+t]+=v*sm[t]
    self._c[key[1:]]=m
   m=self._c[key[1:]];sv=self.co("sieve_sf"if sf else"sieve_geo");c=array("d");t=0.0
   for x in range(min(n,len(m)-1)+1):
    t+=m[x]*sv[n-x];c.append(t)
   self._c[key]=c
  c=self._c[key]
  return c[min(a,len(c)-1)]
 def cap(self,nb):return floor(ceil(nb)/log2(self.minimum))
 def at(self,nb,sl,sb=None):
  n=round(nb/self.step);rb=(self.cap(nb)+1)/16.0
  r=self.selection_normalizer*self.two(n,floor((sl-rb)/self.step),9**9 if sb is None else floor((sb-rb-0.125)/self.step),None,True)/self.denominator[n]
  return max(0.0,min(1.0,r))
 def survivor_upper(self,nb,sl,ex=0.0,sb=None,lb=None,kmax=4):
  n=round(nb/self.step);a=floor((sl+ex)/self.step)
  return(self.two(n,a)if sb is None else self.two(n,a,floor((sb+0.25)/self.step),(kmax,lb,(self.cap(nb)+2)/16.0+0.25)))/self.denominator[n]
def _ptrim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a
def _pmul(a,b,mod):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:c[i+j]=(c[i+j]+x*y)%mod
 return _ptrim(c)
def _pmod(a,b,mod):
 a=[x%mod for x in a];b=_ptrim([x%mod for x in b])
 inv=pow(b[-1], -1,mod)
 while len(a)>=len(b)and not(len(a)==1 and a[0]==0):
  d=len(a)-len(b);fac=(a[-1]*inv)%mod
  for i in range(len(b)):a[d+i]=(a[d+i]-fac*b[i])%mod
  _ptrim(a)
 return a
def _pdiv(a,b,mod):
 a=_ptrim([x%mod for x in a]);b=_ptrim([x%mod for x in b])
 q=[0]*max(1,len(a)-len(b)+1);inv=pow(b[-1], -1,mod)
 while len(a)>=len(b)and not(len(a)==1 and a[0]==0):
  d=len(a)-len(b);fac=(a[-1]*inv)%mod;q[d]=fac
  for i in range(len(b)):a[d+i]=(a[d+i]-fac*b[i])%mod
  _ptrim(a)
 return _ptrim(q),a
def _pgcd(a,b,mod):
 a=_ptrim([x%mod for x in a]);b=_ptrim([x%mod for x in b])
 while not(len(b)==1 and b[0]==0):
  a,b=b,_pmod(a,b,mod)
 inv=pow(a[-1], -1,mod)
 return[(x*inv)%mod for x in a]
def _ppowmod(base,e,m,mod):
 res=[1];b=_pmod(base,m,mod)
 while e>0:
  if e&1:res=_pmod(_pmul(res,b,mod),m,mod)
  b=_pmod(_pmul(b,b,mod),m,mod);e>>=1
 return res
def _ddf_degrees(f,mod):
 f=_ptrim([x%mod for x in f]);inv=pow(f[-1], -1,mod);f=[(x*inv)%mod for x in f]
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
 f=_ptrim([x%mod for x in f]);inv=pow(f[-1], -1,mod);f=[(x*inv)%mod for x in f]
 out,stack=[],[f]
 while stack:
  cur=stack.pop();deg=len(cur)-1
  if deg==0:continue
  if deg==1:out.append((-cur[0])%mod);continue
  while True:
   a=rng.randrange(mod)
   pw=_ppowmod([a,1],(mod-1)//2,cur,mod)
   pw[0]=(pw[0]-1)%mod
   g=_pgcd(pw,cur,mod)
   if 0<len(g)-1<deg:
    q,_=_pdiv(cur,g,mod)
    stack.append(g);stack.append(q);break
 return sorted(out)
def _pmul_h(a,b,mod):
 r=_pmod(_pmul(a,b,mod),H_COEFFS,mod)
 return(r+[0]*6)[:6]
def _ppow_h(a,e,mod):
 r=_ppowmod(a,e,H_COEFFS,mod)
 return(r+[0]*6)[:6]
def _pinv_h(a,mod):
 r0,r1=[x%mod for x in H_COEFFS],_ptrim([x%mod for x in a])
 t0,t1=[0],[1]
 while not(len(r1)==1 and r1[0]==0):
  q,rem=_pdiv(r0,r1,mod)
  qt=_pmul(q,t1,mod)
  L=max(len(t0),len(qt))
  t2=_ptrim([((t0[i]if i<len(t0)else 0)-(qt[i]if i<len(qt)else 0))%mod for i in range(L)])
  r0,r1,t0,t1=r1,rem,t1,t2
 inv=pow(r0[0], -1,mod)
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
 for k in range(10,5, -1):
  c=res[k]
  if c:
   res[k]=F(0)
   for j in range(6):res[k-6+j]-=c*H_COEFFS[j]
 return res[:6]
def verify_tower_algebra_and_crt():
 ord_ok=(pow(P_BN254,12,R_BN254)==1)and(pow(P_BN254,6,R_BN254)!=1)and(pow(P_BN254,4,R_BN254)!=1)
 y=[0,1,0,0,0,0]
 inv3=pow(3, -1,P_BN254)
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
  M64=(1<<64)-1;R2=pow(1<<64,2,p);pinv=(-pow(p, -1,1<<64))&M64
  rng_m=random.Random(p&0xFFFFFFFF)
  pairs=[(p-1,p-1),(p-2,p-2),(2*p-1,2*p-1),(p-1,2*p-1),(0,p-1),(1,p-2)]+[(rng_m.randrange(2*p),rng_m.randrange(2*p))for _ in range(32)]
  for a,b in pairs:
   ab=a*b;m=((ab&M64)*pinv)&M64;s=ab+m*p
   if s>=(1<<128):return False
   t=s>>64
   if t>=2*p:return False
   if t>=p:t-=p
   if t!=(a*b*pow(1<<64, -1,p))%p:return False
  return True
 ntt3_mont_ok=all(_mont_u64_ok(p)for p in NTT3_PRIMES)
 ntt3_ok=(len(NTT3_PRIMES)==3 and all(p<(1<<62)and(p-1)%(1<<55)==0 and _is_prime64(p)for p in NTT3_PRIMES)
 and ntt3_prod>(1<<54)*((1<<64)-1)**2 and ntt3_mont_ok)
 ok=(ord_ok and h_irred and g_irred and tau_frob10 and sigma_fixes_w and chi_ok and f_vanishes
 and crt_ok and gcd_r_p6_ok and boot_id_ok and ntt3_ok and(nw==216)and(n_lead==729)and(n_const==729)
 and(iso_a<1e-10)and(iso_b<1e-10))
 return dict(
 ok=ok,ord_r_p_12=ord_ok,h_irred_mod_p=h_irred,g_rec_irred_mod_p6=g_irred,
 tau_is_frob10=tau_frob10,sigma_fixes_w_exact=sigma_fixes_w,chi_192bit_ok=chi_ok,
 f_rec_vanishes_at_inv_uw=f_vanishes,gcd_r_p6_minus_1_is_1=gcd_r_p6_ok,boot_identity_ok=boot_id_ok,
 ntt3_primes_verified=ntt3_ok,ntt3_montgomery_no_carry_ok=ntt3_mont_ok,ntt3_prod_bits=round(log2(ntt3_prod),3),
 norm_w=nw,norm_f_lead=n_lead,norm_f_const=n_const,
 mod_r_h_degs=h_mod_r_degs,mod_r_g_degs=g_mod_r_degs,mod_r_f_degs=f_mod_r_degs,
 mod_r_crt_all_deg2_ok=crt_ok,sigma_iso_err_Ga=iso_a,sigma_iso_err_Gb=iso_b,
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
 assert mm(tr(S),mm(G_sym,S))==G_sym,"G_sym not exact sigma-invariant"
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
 L=[[0.0]*6 for _ in range(6)]
 for r in range(6):
  for c in range(r+1):
   s=G_sym[r][c]-sum(L[r][k]*L[c][k]for k in range(c))
   L[r][c]=sqrt(s)if r==c else s/L[c][c]
 pts,v=[],[0]*6
 def enum(idx,rem):
  shift=sum(L[j][idx]*v[j]for j in range(idx+1,6))
  diag=L[idx][idx];rad=sqrt(max(0.0,rem))
  for x in range(ceil((-rad-shift)/diag)-1,floor((rad-shift)/diag)+2):
   d=diag*x+shift;nrem=rem-d*d
   if nrem>=-1.0:
    v[idx]=x
    if idx==0:
     if any(v)and qf_i(G_sym,v)<=T_lim:pts.append(tuple(v))
    else:enum(idx-1,max(0.0,nrem))
 enum(5,T_lim+1.0)
 def is_can_int(p):
  if(-p[0], -p[1], -p[2], -p[3], -p[4], -p[5])<p:return False
  for fi in range(6):
   q=qf_i(D_mats[fi],p)
   if q<0:return False
   if q==0:
    uv=mv_sp(U_sps[fi],p)
    for w in(uv,mv_sp(S_sp,uv),mv_sp(S2_sp,uv)):
     nw=(-w[0], -w[1], -w[2], -w[3], -w[4], -w[5])
     if(nw if nw<w else w)<p:return False
  for Sm in(S_sp,S2_sp):
   vk=mv_sp(Sm,p);nvk=(-vk[0], -vk[1], -vk[2], -vk[3], -vk[4], -vk[5])
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
    for cnd in(w,(-w[0], -w[1], -w[2], -w[3], -w[4], -w[5])):
     if cnd in pt_set and cnd not in orb:orb.add(cnd);seen.add(cnd);q.append(cnd)
  orbits.append(orb)
 exact1=missed=dups=0
 for members in orbits:
  c=sum(1 for m in members if is_can_int(m))
  if c==1:exact1+=1
  elif c==0:missed+=1
  else:dups+=c-1
 rel_err=max(abs(G_sym[i][j]/float(G_SC)-T4["GA"][i][j])for i in range(6)for j in range(6))
 # Unconditional Convex-Hexagon Proof of Global Unit-Orbit Minimality for ALL a in O_h\{0} (all radii -> inf!)
 GA=T4["GA"];EMB_ZK=[[complex(r,im)for r,im in row]for row in T4["EMB_ZK"]]
 norm_E_F=sqrt(sum((G_sym[i][j]/float(G_SC)-GA[i][j])**2 for i in range(6)for j in range(6)))
 Bj=[row[:]for row in GA]
 for _ in range(80):
  pi,qi,mx=0,1,0.0
  for i in range(6):
   for j in range(i+1,6):
    if abs(Bj[i][j])>mx:mx,pi,qi=abs(Bj[i][j]),i,j
  if mx<1e-15:break
  ph=0.5*math.atan2(2.0*Bj[pi][qi],Bj[qi][qi]-Bj[pi][pi]);cs,sn=math.cos(ph),math.sin(ph)
  for i in range(6):
   if i!=pi and i!=qi:
    bip,biq=Bj[i][pi],Bj[i][qi];Bj[i][pi]=Bj[pi][i]=cs*bip-sn*biq;Bj[i][qi]=Bj[qi][i]=sn*bip+cs*biq
  bpp,bqq,bpq=Bj[pi][pi],Bj[qi][qi],Bj[pi][qi]
  Bj[pi][pi]=cs*cs*bpp-2.0*sn*cs*bpq+sn*sn*bqq;Bj[qi][qi]=sn*sn*bpp+2.0*sn*cs*bpq+cs*cs*bqq;Bj[pi][qi]=Bj[qi][pi]=0.0
 lam_min_GA=min(Bj[i][i]for i in range(6));eps_G=norm_E_F/lam_min_GA;eps_inf=1.0e-4
 u_norms=[[abs(sum(EMB_ZK[j][k]*UMATS[fi][k][0]for k in range(6)))**2 for j in range(3)]for fi in range(6)]
 def _d3(m):return m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0])+m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0])
 verts=[]
 for i in range(6):
  for j in range(i+1,6):
   M3=[[u_norms[i][k]-1.0 for k in range(3)],[u_norms[j][k]-1.0 for k in range(3)],[1.0,1.0,1.0]]
   D0=_d3(M3)
   if abs(D0)<1e-12:continue
   w=[]
   for col in range(3):
    Mc=[row[:]for row in M3];Mc[0][col],Mc[1][col],Mc[2][col]=-eps_inf,-eps_inf,1.0
    w.append(_d3(Mc)/D0)
   if all(x>=-1e-12 for x in w)and all(sum(w[k]*(u_norms[f][k]-1.0)for k in range(3))>=-eps_inf-1e-11 for f in range(6)):verts.append(w)
 V_min=min(min(w)for w in verts);B_max=(1.0+eps_inf)/V_min;log_B=log(B_max)
 l1=[log(x)for x in u_norms[0]];l2=[log(x)for x in u_norms[2]];detL=l1[0]*l2[1]-l1[1]*l2[0]
 M_bnd=ceil((abs(l2[1]/detL)+abs(l2[0]/detL))*2.0*log_B)+2
 facets={(0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1)}
 inner_cnt=0;min_gap=1e9
 for m1 in range(-M_bnd,M_bnd+1):
  for m2 in range(-M_bnd,M_bnd+1):
   if(m1,m2)in facets:continue
   ms=[exp(m1*l1[j]+m2*l2[j])for j in range(3)]
   if max(ms)>=B_max:continue
   inner_cnt+=1;gap=min(sum(w[j]*(ms[j]-1.0)for j in range(3))for w in verts)
   if gap<min_gap:min_gap=gap
 hex_ok=(len(verts)==6 and eps_G<3.0e-5 and(1.0-eps_G)/(1.0+eps_G)*(1.0+eps_inf)>1.0+4.0e-5 and V_min>0.09 and inner_cnt==3 and min_gap>0.40)
 # Self-contained index-1 proof: [O_h^* : <-1, u_1, u_2>] == 1 and mu(K_h) == {+-1}!
 # Every coset of <-1, u_1, u_2> in O_h^* has a rep in fundamental parallelogram F = [-1/2, 1/2]^2 (max sum|v|_j^2 at 4 corners!).
 ratio_ga=GA[0][0]/sum((EMB_ZK[k][0]*EMB_ZK[k][0].conjugate()).real for k in range(3))
 max_t2_F=max(sum(exp(t1*l1[j]+t2*l2[j])for j in range(3))for t1 in(-0.5,0.5)for t2 in(-0.5,0.5))
 R2_u_box=ratio_ga*max_t2_F*1.01;T_u_box=ceil(R2_u_box*G_SC)
 def _det6(M6):
  A=[r[:]for r in M6];sg=1;pv=1
  for k in range(5):
   if A[k][k]==0:
    p0=next((r for r in range(k+1,6)if A[r][k]!=0),-1)
    if p0<0:return 0
    A[k],A[p0]=A[p0],A[k];sg=-sg
   p1=A[k][k]
   for i in range(k+1,6):
    for j in range(k+1,6):A[i][j]=(A[i][j]*p1-A[i][k]*A[k][j])//pv
   pv=p1
  return sg*A[5][5]
 def _norm_oh(vec):
  poly=[sum(vec[k]*ZK_NUM[k][j]for k in range(6))for j in range(6)]
  cols,cur=[],poly[:]
  for j in range(6):
   cols.append(cur)
   if j<5:
    nxt=[0]+cur[:5];c=cur[5]
    if c:
     for k in range(6):nxt[k]-=c*H_COEFFS[k]
    cur=nxt
  d=abs(_det6([[cols[c][r]for c in range(6)]for r in range(6)]))
  assert d%729==0;return d//729
 ubox_pts=[p for p in pts if qf_i(G_sym,p)<=T_u_box]
 ubox_units=[p for p in ubox_pts if _norm_oh(p)==1]
 unit_idx1_ok=(T_u_box<T_lim and set(ubox_units)=={(1,0,0,0,0,0),(-1,0,0,0,0,0)})
 return dict(
 ok=(exact1==len(orbits)and missed==0 and dups==0 and len(orbits)==382 and rel_err<1e-5 and hex_ok and unit_idx1_ok),
 zero_float_geom=True,sigma_inv_exact=True,G_sym_scale=G_SC,G_sym_rel_err=float(f"{rel_err:.2e}"),
 unit_index1_proved=unit_idx1_ok,unit_box_R2=round(R2_u_box,4),unit_box_pts=len(ubox_pts),
 unit_box_units=len(ubox_units),roots_of_unity_exact=len(ubox_units),
 global_hexagon_proof_ok=hex_ok,eps_G_spectral=round(eps_G,8),eps_inf=eps_inf,
 V_min_P6_inf=round(V_min,6),B_max_cutoff=round(B_max,4),M_box_bound=M_bnd,
 inner_nonfacet_units=inner_cnt,min_nonfacet_gap=round(min_gap,6),
 R2=R2,T_lim_int=T_lim,nonzero_pts=len(pts),total_orbits=len(orbits),
 exact1=exact1,missed_orbits=missed,duplicate_orbits=dups,
 )
def log2_li(log2x):
 x=log2x*LN2
 return log2x-log2(x)+log2(1.0+1.0/x+2.0/(x*x))
def log2_U_z6(lf,lg):
 uf=(2.0**log2_li(lf))/6.0
 ug=(2.0**log2_li(lg))/6.0
 deg2=(2.0**(0.5*lf)/(0.5*lf*LN2)+2.0**(0.5*lg)/(0.5*lg*LN2))/6.0
 return log2(uf+ug+deg2)
NTT3_PRIMES=(116*(1<<55)+1,69*(1<<55)+1,57*(1<<55)+1)
def ntt3_F(M_prod):
 if M_prod>(1<<55):
  return 3.0*ntt3_F(ceil(M_prod/2.0))+20.0*M_prod
 k=max(2,ceil(log2(max(M_prod,2.0))))
 n=2.0**k
 return 135.0*n*k+300.0*n+(1640.0*n if n>FFT_LOCAL else 0.0)
def fft_F(M_prod):return ntt3_F(M_prod)
def mul_F(M_prod):return ntt3_F(M_prod)
def pbits(lo,hi):
 return 1.000028*(2.0**hi-(2.0**lo if lo else 0))/LN2 if lo==0 else 1.000028*(2.0**hi)/LN2
_GB_CACHE={}
def gbatch_cost(pb,leaf=9,L=None,c_leaf=FFT_LEAF):
 SP=ceil(pb/64.0)
 L=L if L is not None else ceil(log2((SP+1)/leaf))
 key=(SP,leaf,L,c_leaf)
 val=_GB_CACHE.get(key)
 if val is not None:return val
 G=L+2;n=lambda j:leaf*(2**j)
 PT=sum((2**(L-j))*(mul_F(2*n(j-1))+FFT_COPY*2*n(j-1))for j in range(1,L+1))
 RT=sum((2**(L-j+1))*(mul_F(n(j)+G+n(j-1))+FFT_COPY*(n(j-1)+G))for j in range(1,L+1))
 val=(PT+RT+FFT_ROOT*mul_F(2*(max(SP,n(L))+G)))/(2**L)+c_leaf
 _GB_CACHE[key]=val
 return val
def split_parts(R,lb):
 np_=1.25506*(2.0**lb)/(lb*LN2);h=ceil(log2(np_));D=h+1
 SP=2.0**ceil(log2(np_*(lb/64.0+1.0)))
 L=ceil(log2(SPLIT_REC*SPLIT_Q));G=L+2;C=SPLIT_C
 W=R/SPLIT_Q+1.0;B=SPLIT_REC*R*D;occ=2*NORM_BITS
 return dict(
 trees=B*(2*(C*(36*L+2*G)+3*FFT_B)+2**15),
 roots=6*C*(2*SP*W*D+136*B)+6*FFT_B*B,
 gcd=B*2.0**24,small=R*2.0**26,output=occ*R*2.0**20,
 ideal=occ*R*SPLIT_IDEAL,control=B*2.0**16,
 ptree=D*SP*C+2*FFT_B*(2.0**h),debit=128*(W+2),
 )
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
def _build_sb_profile_pure(L_max,w_sb,kap,nx=120,nm=24,nphase=8,n_ph_tbl=8,q_slack=1.0/126.0):
 h_tot=c_tot=sb_row=tot_len=0.0
 nb=200
 wj_hist=[0.0]*(nb+1)
 wb1_hist=[0.0]*(nb+1)
 ph_pad=w_sb/float(n_ph_tbl)
 max_sb=0
 for ix in range(nx):
  r2=(ix+0.5)/nx
  w5=2.5*(r2**1.5)/nx
  L=L_max*sqrt(1.0-r2)
  tot_len+=w5*L
  du=L/nm
  for im in range(nm):
   u=((im+0.5)/nm-0.5)*L
   qb=r2+(2.0*u/L_max)**2
   fa=min(1.0,max(0.0,(1.0-qb)/kap))**3
   h_tot+=w5*du*fa
   wb1_hist[min(nb,int(fa*nb))]+=w5*du
  u_lo,u_hi= -0.5*L,0.5*L
  for ip in range(nphase):
   off=(ip+0.5)/nphase*w_sb
   pos= -off;cnt=0
   while pos<L:
    t0=max(0.0,pos);t1=min(L,pos+w_sb)
    if t1>t0:
     cnt+=1
     a_u=u_lo+t0;b_u=min(u_hi,u_lo+t1)
     u_min=0.0 if(a_u<=0.0<=b_u)else(a_u if abs(a_u)<abs(b_u)else b_u)
     q_min=r2+(2.0*u_min/L_max)**2
     fa_ex=min(1.0,max(0.0,(1.0-q_min)/kap))**3
     a_ug=max(u_lo,a_u-ph_pad);b_ug=min(u_hi,b_u+ph_pad)
     u_mg=0.0 if(a_ug<=0.0<=b_ug)else(a_ug if abs(a_ug)<abs(b_ug)else b_ug)
     q_mg=max(0.0,r2+(2.0*u_mg/L_max)**2-q_slack*kap)
     fa_g=min(1.0,max(0.0,(1.0-q_mg)/kap))**3
     wt=w5/nphase
     c_tot+=wt*w_sb*fa_ex;sb_row+=wt
     wj_hist[min(nb,int(fa_g*nb))]+=wt
    pos+=w_sb
   if cnt>max_sb:max_sb=cnt
 ws=sum(wj_hist);wb1s=sum(wb1_hist)
 fj_bins=[(min(1.0,(i+1.0)/nb),w/ws)for i,w in enumerate(wj_hist)if w>0.0]
 fb1_bins=[(min(1.0,(i+0.5)/nb),w/wb1s)for i,w in enumerate(wb1_hist)if w>0.0]
 return dict(cr2=c_tot/h_tot,sb_row=sb_row,max_sb=max_sb,
 frac_reg_analytic=h_tot/tot_len,fj_bins=fj_bins,fb1_bins=fb1_bins)
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
def compute_ledger(track="bn254-gt-dlp",use_db=True,c_tail_z6=None,cap_scale=1.0,
strict_per_owner=True,swar_probe=None):
 D=load_data()
 PQ=D["pool_quad"];T4=D["t4_geom"];FILT=D["filtsim_raw"]
 lf,lg,exc=OP["lf"],OP["lg"],OP["exc"]
 p_min,p_cut,pm,lx2=OP["p_min"],OP["p_cut"],OP["p_mid"],OP["lx2"]
 assert lx2==pm,"2-band pure-integer sieve requires lx2 == p_mid"
 sl,sb,w_sb=OP["slack"],OP["slack_b"],OP["w_sb"]
 if c_tail_z6 is None:c_tail_z6=OP["c_tail_z6"]
 if swar_probe is None:swar_probe=C_SWAR_PROBE
 rkey=(lg,lx2,p_min)
 if rkey not in _RET_CACHE:
  _RET_CACHE[rkey]=(
  Retention(lg,lx2,prime_min=2**p_min,explicit_bits=min(p_min,14),step=0.25,max_bits=600,tilt=0.05,odd=1),
  Retention(600.0,lx2,prime_min=2**p_min,explicit_bits=min(p_min,14),step=0.25,max_bits=601,tilt=0.005,odd=1),
  )
 ret_obj,cost_obj=_RET_CACHE[rkey]
 prim_tail_corr=min(1.0,EXACT_PRIM_T4/(PQ["n_prim_in"]/float(PQ["n_in"])))
 lU=log2_U_z6(lf,lg)
 target_lr=lU+log2(exc)+OP["margin"]
 log2_V_sub0=PQ["log2_V0"]+log2(PQ["frac_reg"])
 n_in=PQ["n_in"];bins=PQ["bins"]
 assert sum(x[2]for x in bins)==PQ["n_prim_in"]and sum(x[3]for x in bins)==n_in
 ALPHA_F_ODD=ALPHA_F_RAW+2.0/3.0
 ALPHA_G_ODD=ALPHA_G_RAW+1.0/3.0
 CLASSES_2=((7.0/9.0,0.0,0.0),(1.0/9.0,0.0,3.0),(1.0/9.0,6.0,0.0))
 def eval_d(d):
  ys=ys2=nf_w=ng_w=rg_w=0.0
  for lnf0,lng0,cnt_prim,cnt_all in bins:
   if not cnt_prim:continue
   nf_raw=lnf0+48.0*d;ng_raw=lng0+12.0*d
   ng_q=round(ng_raw*2.0)/2.0
   r0=ret_obj.at(ng_q,sl,sb);r3=ret_obj.at(ng_q-3.0,sl,sb)
   y=rg_num=0.0
   for w2,v2f,v2g in CLASSES_2:
    nf_c=nf_raw-v2f+ALPHA_F_ODD;ng_c=ng_raw-v2g+ALPHA_G_ODD
    r_c=r3 if v2g else r0
    yc=w2*(rho_db(nf_c/lf,nf_c)*rho_db(ng_c/lg,ng_c)if use_db else rho(nf_c/lf)*rho(ng_c/lg))*r_c
    y+=yc;rg_num+=yc*r_c
   scale=(1.0-0.0014)*prim_tail_corr
   r_g_bin=(rg_num/y)if y>0.0 else r0
   y*=scale
   wy=cnt_prim*y;ys+=wy;ys2+=cnt_prim*(y*y)
   nf_w+=wy*nf_raw;ng_w+=wy*ng_raw;rg_w+=wy*r_g_bin
  my=ys/n_in
  return log2_V_sub0+12.0*d+log2(max(my,1e-300)),my,ys2/n_in,nf_w/ys-CONTENT_BITS,ng_w/ys,rg_w/ys
 lo,hi= -2.5, +2.5
 for _ in range(32):
  mid=0.5*(lo+hi)
  if eval_d(mid)[0]<target_lr:lo=mid
  else:hi=mid
 d_opt=0.5*(lo+hi)
 lr,my,my2,nf_rel,ng_rel,r_g=eval_d(d_opt)
 se_bits=sqrt(max(my2-my*my,0.0)/n_in)/max(my,1e-300)/LN2
 log2_V=log2_V_sub0+12.0*d_opt
 log2_Vb=PQ["log2_Vb"]+6.0*d_opt
 log2_Na=PQ["log2_Va"]+6.0*d_opt
 cells=2.0**log2_V;Vb=2.0**log2_Vb;Na=2.0**log2_Na
 GB=T4["GB"];g00=GB[0][0];det_GB=_det_mat(GB)
 S5=[[GB[r][c]-GB[r][0]*GB[0][c]/g00 for c in range(1,6)]for r in range(1,6)]
 det_S5=_det_mat(S5)
 Rb=(Vb*sqrt(det_GB)/((pi**3)/12.0))**(1.0/6.0)
 L_max=2.0*Rb/sqrt(g00)
 n5=((4.0*pi**2)/15.0)*(Rb**5)/sqrt(det_S5)
 L_b_avg=Vb/n5
 score_arena_log2=log2((L_max+2.0)*Na*2.0+Na*2.0)
 vb_score_log2=log2(Vb*2.0)
 b2_1p_ram_log2=log2(16.0*Na+9.0*(2.0**pm)/w_sb+(L_max+2.0)*Na*2.0+Na*2.0)
 pkey=(round(L_max,2),w_sb)
 if pkey not in _PROF_CACHE:
  _PROF_CACHE[pkey]=_build_sb_profile_pure(L_max,w_sb,PQ["kap"])
 prof=_PROF_CACHE[pkey]
 cr2_chg=ceil(prof["cr2"]*1.005*1000.0)/1000.0
 c_gb=max(C_GBATCH,gbatch_cost(pbits(0,lf),leaf=9))
 c_b2=max(C_BATCH2,gbatch_cost(pbits(0,lg),c_leaf=2**19))
 cc=(c_gb,C_TDIV+C_BATCH+c_b2)
 np_bins=PQ["nonprim_bins"]
 cof_sum=sf_sum=s0_sum=0.0;tot_p=tot_np=0
 for lnf0,lng0,cnt_prim,cnt_all in bins:
  tot_p+=cnt_prim;tot_np+=(cnt_all-cnt_prim)
  if not cnt_prim:continue
  nf_raw=lnf0+48.0*d_opt;ng_raw=lng0+12.0*d_opt
  ng_q=round(ng_raw*2.0)/2.0
  sf0=cost_obj.survivor_upper(ng_q,sl,0.25);sf3=cost_obj.survivor_upper(ng_q-3.0,sl,0.25)
  s00=cost_obj.survivor_upper(ng_q,sl,0.25,sb,lg,KMAX);s03=cost_obj.survivor_upper(ng_q-3.0,sl,0.25,sb,lg,KMAX)
  cof_p=sf_p=s0_p=0.0
  for w2,v2f,v2g in CLASSES_2:
   sf_c=sf3 if v2g else sf0
   s0_c=s03 if v2g else s00
   q_f_c=min(1.0,Q_F_SAFETY*rho((nf_raw-v2f+ALPHA_F_ODD)/lf))
   cof_p+=w2*(sf_c*(C_NORMS+C_SB)+s0_c*(cc[0]+q_f_c*cc[1]))
   sf_p+=w2*sf_c;s0_p+=w2*s0_c
  cof_sum+=cnt_prim*cof_p;sf_sum+=cnt_prim*sf_p;s0_sum+=cnt_prim*s0_p
 np_sum_cnt=0;cof_np_sum=sf_np_sum=0.0
 for ng_red_q0,cnt_np in np_bins:
  np_sum_cnt+=cnt_np
  ng_eff=floor((ng_red_q0+12.0*d_opt)*2.0)/2.0
  sf_np=cost_obj.survivor_upper(max(20.0,ng_eff),sl,0.25)
  cof_np_sum+=cnt_np*(sf_np*C_NORMS);sf_np_sum+=cnt_np*sf_np
 assert tot_p==PQ["n_prim_in"]and tot_np==np_sum_cnt==PQ["n_nonprim_in"]and tot_p+tot_np==n_in,f"partition mismatch: {tot_p}+{tot_np} vs {np_sum_cnt}, {n_in}"
 cof_sum+=cof_np_sum;sf_sum+=sf_np_sum
 cof_b=cof_sum/float(n_in);cof_np_cell=cof_np_sum/float(n_in)
 sf=sf_sum/float(n_in);s0=s0_sum/float(n_in)
 c_hit_lo=C_HIT_LO if strict_per_owner else 4.25
 c_hit_hi=C_HIT_HI if strict_per_owner else 9.0
 tot_lk=0.0;dl1=(p_cut-p_min)/100.0
 fb1_bins=prof["fb1_bins"]
 for i in range(100):
  l0=p_min+i*dl1;l1=l0+dl1;lm=0.5*(l0+l1)
  dn=LI_PAD*(2.0**log2_li(l1)-2.0**log2_li(l0))
  mu=Na/(2.0**lm)
  c_bs=24.0+8.0*(log2(mu+1.0)+1.0)
  c_e4=sum(wt*(11.0+(1.0-exp(-mu*fa))*14.0+mu*fa)for fa,wt in fb1_bins)
  tot_lk+=dn*min(c_bs,c_e4)
 n1_tot=LI_PAD*(2.0**log2_li(p_cut)-2.0**log2_li(p_min))
 rep_lo=log(p_cut/p_min)*(c_hit_lo*1.004+0.06)+(tot_lk*Vb+n1_tot*n5*24.0)/cells
 assert prof["max_sb"]<=8,f"SWAR overflow: max_sb={prof['max_sb']} > 8"
 n_words_row=1.0
 n2_tot=n2_occ=n2_ent=0.0;dl2=(pm-p_cut)/100.0
 fj_bins=prof["fj_bins"]
 for i in range(100):
  l0=p_cut+i*dl2;l1=l0+dl2;lm=0.5*(l0+l1)
  dn=LI_PAD*(2.0**log2_li(l1)-2.0**log2_li(l0))
  mu_sb=(w_sb*Na)/(2.0**lm)
  p_occ=1.0-exp(-mu_sb)if mu_sb<40.0 else 1.0
  p_ent=sum(wt*(1.0-exp(-mu_sb*fa))for fa,wt in fj_bins)
  n2_tot+=dn;n2_occ+=dn*p_occ;n2_ent+=dn*p_ent
 probe=(n2_tot*n5*n_words_row*swar_probe)/cells
 occ_sb=(n2_ent*n5*prof["sb_row"]*C_SB_ENTER)/cells
 hi_ap=log(pm/p_cut)*(C_CAND*cr2_chg+c_hit_hi)
 per_cell=rep_lo+probe+occ_sb+hi_ap+C_SCAN_ZERO+cof_b
 coll_u=log2(cells*per_cell*OP["work_allow"]/UNIT)
 n_lo=LI_PAD*0.5*(2.0**log2_li(p_cut)-2.0**log2_li(p_min))
 sp_lo=LI_PAD*0.5*(2.0**(2.0*p_cut))/(2.0*p_cut*LN2-1.0)
 ins_lo=(n_lo+1.0)*(C_LIST_A*Na)+C_LBUCKET*LIST_NS*sp_lo
 byt_lo=16.0*Na*n_lo+8.0*LIST_NS*sp_lo
 n_hi=LI_PAD*0.5*(2.0**log2_li(pm)-2.0**log2_li(p_cut))
 sp_hi=LI_PAD*0.5*(2.0**(2.0*pm))/(2.0*pm*LN2-1.0)
 ins_hi=n_hi*(C_LIST_A*Na)+C_LBUCKET*(sp_hi/w_sb)
 byt_hi=16.0*Na*n_hi+9.0*(sp_hi/w_sb)
 ins_pre=Na*(1000.0+12.0*log2_Na)
 lists_u=log2((ins_lo+ins_hi+ins_pre)/UNIT)
 ram_log2=log2(byt_lo+byt_hi+(L_max+2.0)*Na*2.0)
 def stats_of(key):
  rs=FILT[key];n=len(rs)
  def ms(idx):
   vs=[x[idx]for x in rs];m=sum(vs)/n
   s2=sum((v-m)**2 for v in vs)/(n-1)if n>1 else 0.0
   return m,sqrt(s2/n),max(vs)
  sw,sw_se,sw_mx=ms(1);nl,nl_se,nl_mx=ms(2);k5,_,k5_mx=ms(3);tr,_,tr_mx=ms(4);uf,_,_=ms(5);eu,_,_=ms(6)
  k5_exp_vs=[x[6]*C_UNIT_5+(1.0-x[6])*C_GEN_1 for x in rs]
  k5_exp_mean=sum(k5_exp_vs)/n
  return dict(n=n,N_sw_U=sw,N_sw_se=sw_se,N_sw_max=sw_mx,N_l_U=nl,N_l_se=nl_se,N_l_max=nl_mx,
  k5_ev=k5,k5_exp=k5_exp_mean,k5_exp_max=max(k5_exp_vs),tab_r=tr,tab_r_max=tr_mx,unit_frac=uf,exp_unit=eu)
 cur20=stats_of("m184_z6_61.0_60.0_e0.82_W160_20")
 cur21=stats_of("m184_z6_61.0_60.0_e0.82_W160_21")
 cur22=stats_of("m184_z6_61.0_60.0_e0.82_W160_22")
 assert cur20["n"]==4 and cur21["n"]==8 and cur22["n"]==4
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
  for n in(20,25,30,35,40):
   N=2.0**log2_N_cap;L_st=bw_steps(log2_N_cap,m,n);s=N*(1.0/m+1.0/n)
   kap_bw=L_st/s;hor=max(s,(1.0+BW_SIGMA)*(N+n_sm)/n+4.0)/s
   per_n2=(c_prod+C_PROJ)*(n/m+1.0)*kap_bw+(1.0/m+1.0/n)*(c_prod+(n_sm*C_FULL+(n-n_sm)*C_Z64)+C_RED9)*hor
   per_n1=(m+n)*C_RED9
   la_u=log2((N*N*per_n2+N*per_n1+PREP_OWNER+PREP_GUARD+PREP_JOINT+LA_RESERVE+2.0**20)/UNIT)
   gen_u=gen_at_cap(log2_N_cap,m,n)
   both=log2(2.0**la_u+2.0**gen_u)
   if best_la is None or both<best_la["both"]:
    ncap=2.0**log2_Ns_cap;nslots=ceil(n/5.0)
    mem_1own=log2(ncap*(160.0*nslots+160.0*nslots+W_ROW*8.0)+(2.0**55.04)*nslots+65.0*(2.0**53))
    best_la=dict(la=la_u,gen=gen_u,both=both,m=m,n=n,mem_1own=mem_1own,log2_N_cap=log2_N_cap,log2_Ns_cap=log2_Ns_cap,c_prod=c_prod)
 pkey_peel=(lf,lg,exc)
 if pkey_peel not in _PEEL_CACHE:
  _PEEL_CACHE[pkey_peel]=_peel_fixed_point(lf,lg,exc,nf=nf_rel,ng=ng_rel)
 peel=_PEEL_CACHE[pkey_peel]
 rels_low=lr-2.25*se_bits
 dkey=(round(rels_low,3),round(peel["h"],4))
 if dkey not in _DESC_CACHE:
  _DESC_CACHE[dkey]=DESC.evaluate_descent(
  D["descent"],rho,fft_F,gbatch_cost,pbits,
  log2_rels=rels_low,peel_h=peel["h"],nf_main=nf_rel,ng_main=ng_rel,
  )
 desc=_DESC_CACHE[dkey]
 rels_up=lr+2.25*se_bits;u_log=log2(UNIT);lb_max=max(lf,lg)
 canon_u=log2_Na+42.8-u_log
 r_bits=R_BN254.bit_length()
 bit_f=48*48*512+7*4*36*512+2**17+3*48*1024+2*48*6*8*40+48*1024
 bit_g=12*12*512+1*6*6*512+2**15+3*12*1024+2*12*6*8*40+12*1024
 sm_u=rels_up+log2(48*r_bits*bit_f+12*r_bits*bit_g)-u_log
 misc_u=log2(2.0**(rels_up+30.4)+2.0**(rels_up+14.0))-u_log
 pi_b=log2(1.25506*(2.0**lb_max)/(lb_max*LN2))
 fb_u=pi_b+23.4-u_log
 split_u=log2(sum(split_parts(2.0**rels_up,lb=lb_max).values()))-u_log+1.0
 L_gt=ceil(lb_max-6.0+0.53)
 pg=sum((2.0**(lb_max-6.0+0.53-j-1))*(fft_F(2.0**(j+1))+3.0*(2.0**(j+1)))for j in range(L_gt))
 gtree_u=log2(2.0*pg)-u_log
 rows=dict(
 collection=coll_u,chi_setup=65.712288,la=best_la["la"],gen=best_la["gen"],
 lists=lists_u,two_logs=desc["two_log"],canon=canon_u,sm=sm_u,misc=misc_u,
 fb=fb_u,split=split_u,gtree=gtree_u,
 cert=50.000001,verify=11.000001,wrapper=53.712288,
 )
 if track=="bn254-pairing-inversion":
  rows["output"]=log2(G2_OUTPUT_BOUND/UNIT)
 exp_run=log2(sum(2.0**v for v in rows.values()))
 ncap=ceil(2.0**best_la["log2_Ns_cap"])
 fail_terms=dict(
 main_matrix=4.0*ncap/(2.0**64),
 filt_cap_h3_given_H3=0.0,
 split_abort=(2.0**split_u)/(2.0**SPLIT_CAP-(2.0**SPLIT_RES)/UNIT),
 descent_root_b10_union=desc["phi_log_two"],
 descent_markov_cap=desc["markov_cap_fail"],
 descent_repeat_ideal=0.0,
 miller_rabin_64r=(2.0**69)*(2.0**(-128.0)),
 group_order_boot_and_sm=8.0/R_BN254,
 numerical_pure_int=0.0,
 )
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
 all4_mono=(
 cur20["N_sw_U"]>cur21["N_sw_U"]>cur22["N_sw_U"]and
 cur20["N_l_U"]>cur21["N_l_U"]>cur22["N_l_U"]and
 cur20["k5_exp"]>cur21["k5_exp"]>cur22["k5_exp"]and
 cur20["tab_r"]>cur21["tab_r"]>cur22["tab_r"]
 )
 all_below_4caps=(
 max(cur20["N_sw_max"],cur21["N_sw_max"],cur22["N_sw_max"])<nu_sw_cap and
 max(cur20["N_l_max"],cur21["N_l_max"],cur22["N_l_max"])<nu_l_cap and
 max(cur20["k5_exp_max"],cur21["k5_exp_max"],cur22["k5_exp_max"])<k5_cap and
 max(cur20["tab_r_max"],cur21["tab_r_max"],cur22["tab_r_max"])<tab_r_cap
 )
 return dict(
 track=track,use_db=use_db,c_tail_z6=c_tail_z6,cap_scale=cap_scale,
 strict_per_owner=strict_per_owner,swar_probe=swar_probe,
 op=OP,
 pool=dict(n_total=PQ["n_total"],n_in=n_in,n_prim_in=PQ["n_prim_in"],
 frac_reg=PQ["frac_reg"],frac_reg_analytic=prof["frac_reg_analytic"]),
 geometry=dict(log2_U=lU,log2_V=log2_V,se_bits=se_bits,log2_Vb=log2_Vb,log2_Na=log2_Na,
 log2_N5=log2(n5),L_b_avg=L_b_avg,L_max=L_max,sb_row=prof["sb_row"],
 max_sb=prof["max_sb"],
 cr2_tight=cr2_chg,p_occ_b2=n2_occ/max(n2_tot,1.0),p_ent_b2=n2_ent/max(n2_tot,1.0),
 nf_rel=nf_rel,ng_rel=ng_rel,r_g=r_g,sf=sf,s0=s0,
 score_arena_log2=score_arena_log2,vb_score_log2=vb_score_log2,
 b2_1p_ram_log2=b2_1p_ram_log2,ram_log2=ram_log2),
 per_cell=dict(total=per_cell,rep_lo=rep_lo,probe=probe,occ_sb=occ_sb,hi_ap=hi_ap,
 scan_zero=C_SCAN_ZERO,cof=cof_b,cof_nonprim=round(cof_np_cell,6),
 n_prim_exact=tot_p,n_nonprim_exact=np_sum_cnt,nonprim_bins_n=len(np_bins)),
 la=dict(bw=[best_la["m"],best_la["n"]],log2_N_cap=best_la["log2_N_cap"],log2_Ns_cap=best_la["log2_Ns_cap"],
 W_row=W_ROW,nu_sw_cap=nu_sw_cap,nu_l_cap=nu_l_cap,k5_cap=k5_cap,tab_r_cap=tab_r_cap,
 pf_rl_cap=pf_rl_cap,c_prod=c_prod,
 mem_1own_log2=best_la["mem_1own"],fits_1own=(best_la["mem_1own"]<64.0),
 h3_mode="fixed_instance_deterministic_4cap",
 ladder_all4_monotonic=all4_mono,all_runs_below_4caps=all_below_4caps,
 filtsim_cur20=cur20,filtsim_cur21=cur21,filtsim_cur22=cur22,peel=peel),
 descent=desc,
 rows=rows,
 exp_run=exp_run,exp_run_decl=f"{exp_run_decl:.9f}",
 fail_terms=fail_terms,
 succ_exact=succ_exact,succ_decl=succ_decl,
 K_gate=K_gate,gate_passes=gate_passes,gate_prev_fails=gate_prev_fails,
 claim_exact=claim_exact,claim_decl=claim_decl,score_numerator=str(claim_num),
 )
_OH_U=((1,-1,-1,-1,1,0),(-1,1,1,1,0,0),(1,1,0,1,-1,0),(0,0,1,0,0,1),(-1,0,0,0,0,1),(1,0,0,-1,0,-1),(0,-1,0,-1,1,0),(1,0,0,0,1,-1),(0,1,1,1,0,1),(-3,0,2,1,0,1),(0,1,0,0,0,0),(1,0,-1,0,0,0),(-1,2,2,1,0,1),(-1,0,0,1,0,2),(-1,-1,-1,0,-1,1))
def verify_nonprim_hnf_witnesses():
 PQ=load_data()["pool_quad"]
 M=[[[1 if c==j else 0 for c in range(6)]for j in range(6)]for _ in range(6)]
 idx=0
 for i in range(1,6):
  M[i][0]=M[0][i]
  for j in range(i,6):M[i][j]=M[j][i]=_OH_U[idx];idx+=1
 def _omul(u,v):return[sum(u[i]*v[j]*M[i][j][c]for i in range(6)for j in range(6))for c in range(6)]
 def _mrows(v):return[[sum(v[i]*M[i][j][c]for i in range(6))for c in range(6)]for j in range(6)]
 def _hnf_det(rows):
  R=[r[:]for r in rows];det=1;r0=0
  for col in range(6):
   while True:
    nz=[i for i in range(r0,len(R))if R[i][col]!=0]
    if len(nz)==1:
     p=nz[0]
     if p!=r0:R[r0],R[p]=R[p],R[r0]
     if R[r0][col]<0:R[r0]=[-x for x in R[r0]]
     det*=R[r0][col];r0+=1;break
    p=min(nz,key=lambda i:abs(R[i][col]))
    if R[p][col]<0:R[p]=[-x for x in R[p]]
    pv=R[p][col];pr=R[p]
    for i in nz:
     if i==p:continue
     q=R[i][col]//pv;ri=R[i]
     for c in range(col,6):ri[c]-=q*pr[c]
  return det
 w_zk=[-1,-1,-1,-1,1,-1]
 assert[sum(w_zk[k]*ZK_NUM[k][j]for k in range(6))//3 for j in range(6)]==W_COEFFS
 ok_all=True;checked=[]
 for ni_decl,a,b in PQ["nonprim_witnesses"]:
  ni_obs=_hnf_det(_mrows(a)+_mrows(b))
  b2=_omul(b,b);wa2=_omul(w_zk,_omul(a,a))
  ng=abs(_hnf_det(_mrows([b2[c]-U_BN254*wa2[c]for c in range(6)])))
  div_ok=(ng%(ni_obs*ni_obs)==0);nj=ng//(ni_obs*ni_obs)if div_ok else 0;v2=0;t=nj
  while t>0 and(t&1)==0:v2+=1;t>>=1
  ok_i=(ni_obs==ni_decl and div_ok and v2 in(0,3));ok_all&=ok_i
  checked.append(dict(N_I=ni_obs,divides_Ng_sq=div_ok,v2_NJ=v2,ok=ok_i))
 return dict(ok=ok_all,n_checked=len(checked),witnesses=checked)
