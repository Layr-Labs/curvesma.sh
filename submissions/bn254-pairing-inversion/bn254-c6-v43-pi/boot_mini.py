"""boot_mini.py: exhaustive miniature MONIC_M1 boot census at u=-3 (p=2143,r=2089) & u=5 (p=27631,r=27481) for H-EB0..H-EB2.
e = j+1 for j = 0..(r-3)/2 covers one representative per sign class (e and r-e give equivalent lattices under b -> -b); exact factoring decides 2^L-smoothness vs p = rho((log2|N_g|+ALPHA_G)/L).
boot_mini*.b85 store records; stats() recomputes all figures; replay(u,js) recomputes from scratch. Scope: plain 2^L-smoothness of 60..120-bit norms; it
tests the per-vector random-integer model and the dependence between kept vectors of one lattice, not production
sizes."""
import base64,itertools,json,lzma,math
from math import log2
from pathlib import Path
import model_bn254_v6 as M,descent_bn254_v6 as DESC,boot_dp as BD
TOWERS=(-3,5)
LS=(18,20,22)
_PR=[None]
def _chunks():
 """Primes <= 2^max(LS) in ascending chunks of 256 with their products (gcd screen, same hits as testing each q)."""
 if _PR[0] is None:
  n=1<<max(LS);sv=bytearray(b"\x01")*(n+1);sv[:2]=b"\x00\x00"
  for q in range(2,int(n**0.5)+1):
   if sv[q]:sv[q*q::q]=b"\x00"*((n-q*q)//q+1)
  pr=list(itertools.compress(range(n+1),sv));ch=[]
  for i in range(0,len(pr),256):
   c=1
   for q in pr[i:i+256]:c*=q
   ch.append((c,pr[i:i+256]))
  _PR[0]=ch
 return _PR[0]
def tower_ok(u):
 """p prime, h irreducible mod p (no factor of degree 1, 2 or 3), MONIC_M1 radicand -((u+1)w+2) a non-square in F_p[y]/h."""
 p=36*u**4+36*u**3+24*u**2+6*u+1;pm=M._pmul_h
 if p<5 or any(p%q==0 for q in range(2,int(p**0.5)+1)):return False
 def pw(a,e):
  r=[1,0,0,0,0,0]
  for bt in bin(e)[2:]:
   r=pm(r,r,p)
   if bt=="1":r=pm(r,a,p)
  return r
 yk=[0,1,0,0,0,0]
 for _ in range(3):
  yk=pw(yk,p)
  d=[(yk[i]-(1 if i==1 else 0))%p for i in range(6)]
  h=[x%p for x in M.H_COEFFS]
  a,b=h,d
  while any(b):
   while b and b[-1]==0:b=b[:-1]
   if not b:break
   while len(a)>=len(b)and any(a):
    t=a[-1]*pow(b[-1],-1,p)%p;s=len(a)-len(b)
    a=[(a[i]-t*b[i-s])%p if i>=s else a[i]for i in range(len(a))]
    while a and a[-1]==0:a=a[:-1]
   a,b=b,a
  if len(a)>1:return False
 return pw([(-((u+1)*x+(2 if i==0 else 0)))%p for i,x in enumerate(M.W_COEFFS)],(p**6-1)//2)==[p-1,0,0,0,0,0]
def record(j,u,sizes=None):
 c=BD.ctx(u);B=BD.lattice_e(c,j+1)
 seen=set();keep=[]
 for v in BD.cands(B):
  if max(abs(z)for z in v[:6]).bit_length()>=128 or max(abs(z)for z in v[6:]).bit_length()>=160:continue
  N=BD.ng(c,v)
  if N==0:continue
  new=N not in seen;seen.add(N)
  if new and M.is_ideal_primitive_oh(v[:6],v[6:]):keep.append(N)
 if sizes is not None:sizes+=[log2(N)for N in keep]
 P=1
 for N in keep:P*=N
 hit=[]
 for cp,qs in _chunks():
  g=math.gcd(P,cp)
  if g>1:hit+=[q for q in qs if g%q==0]
 lp=[];s16=0.0
 for N in keep:
  x=N;big=1
  for q in hit:
   if x%q==0:
    big=q
    while x%q==0:
     x//=q
     if q<(1<<16):s16+=log2(q)
  lp.append(big if x==1 else None)
 out=[]
 for L in LS:
  ps=[M.rho((log2(N)+DESC.ALPHA_G)/L)for N in keep]
  m=sum(1 for b in lp if b is not None and b<(1<<L));E=sum(ps);Q=(E*E-sum(x*x for x in ps))/2.0
  out.append([m,m*(m-1)//2,round(100*log2(E)),round(100*log2(Q))if Q>0 else -10**9])
 return out+[[len(keep),round(10*s16)]]
_S=[None]
def _dec(T):
 """Column store -> rows: M/C sparse [[j, v]], e2 and K, s plain, e_i - e2 (i < 2), q_i - 2 e_i."""
 n=T["n"];nL=len(LS);rows=[[[0,0,0,0]for _ in range(nL)]+[[0,0]]for _ in range(n)]
 for i in range(nL):
  for j,v in T["M%d"%i]:rows[j][i][0]=v
  for j,v in T["C%d"%i]:rows[j][i][1]=v
 for j in range(n):
  e2=T["e2"][j];rows[j][nL-1][2]=e2
  for i in range(nL-1):rows[j][i][2]=e2+T["d%d"%i][j]
  for i in range(nL):rows[j][i][3]=T["g%d"%i][j]+2*rows[j][i][2]
  rows[j][nL]=[T["K"][j],T["s"][j]]
 return dict(p=T["p"],r=T["r"],n=n,rec=rows)
def load():
 if _S[0] is None:
  parts={}
  for f in sorted((Path(__file__).resolve().parent).glob("boot_mini*.b85")):
   for u,T in json.loads(lzma.decompress(base64.b85decode(f.read_bytes())).decode()).items():
    if u in parts:
     for k,v in T.items():
      if isinstance(v,list):parts[u][k]+=v
    else:parts[u]=T
  _S[0]={u:_dec(T)for u,T in parts.items()}
 return _S[0]
def replay(u,js):
 S=load()[str(u)];sz=[];bad=[j for j in js if record(j,u,sizes=sz)!=S["rec"][j]]
 return dict(u=u,replayed=list(js),mismatches=bad,kept_norm_bits_min_mean_max=[round(min(sz),2),round(sum(sz)/len(sz),2),round(max(sz),2)],ok=not bad)
def _pois(k):
 """Two-sided 95% interval for a Poisson mean given k (exact for k <= 300, normal beyond)."""
 if k>300:return k-1.96*math.sqrt(k),k+1.96*math.sqrt(k)
 def cdf(n,m):
  t=s=math.exp(-m)
  for i in range(1,n+1):t*=m/i;s+=t
  return s
 def bis(f,lo,hi):
  for _ in range(80):
   mid=(lo+hi)/2
   if f(mid):hi=mid
   else:lo=mid
  return (lo+hi)/2
 lo=0.0 if k==0 else bis(lambda m:1-cdf(k-1,m)>=0.025,0.0,3.0*k+50)
 return lo,bis(lambda m:cdf(k,m)<=0.025,0.0,3.0*k+50)
def _wide(k,num,den,nl):
 """95% interval for sum(num)/sum(den): the wider of the Poisson interval of k = sum(num) and the lattice-clustered
 normal interval of the ratio estimator (lattices are the independent units)."""
 D=sum(den);r=sum(num)/D;lo,hi=_pois(k)
 se=math.sqrt(sum((a-r*b)**2 for a,b in zip(num,den))/(nl-1)/nl)/(D/nl)
 return [round(min(lo/D,r-1.96*se),4),round(max(hi/D,r+1.96*se),4)]
def stats():
 out={};er=sum(log2(q)/(q-1)for q in M._primes_upto(1<<16))
 for u,S in sorted(load().items(),key=lambda t:int(t[0])):
  R=S["rec"];nl=len(R);by={}
  for i,L in enumerate(LS):
   Ms=[r[i][0]for r in R];Cs=[r[i][1]for r in R];Hs=[1 if r[i][0]>0 else 0 for r in R]
   Es=[2.0**(r[i][2]/100.0)for r in R];Qs=[2.0**(r[i][3]/100.0)if r[i][3]>-10**8 else 0.0 for r in R]
   Bs=[e-q for e,q in zip(Es,Qs)];E,Q,B=sum(Es),sum(Qs),sum(Bs)
   by[str(L)]=dict(kept_smooth=sum(Ms),model_sum=round(E,2),first_moment_ratio=round(sum(Ms)/E,4),first_moment_ci95=_wide(sum(Ms),Ms,Es,nl),
lattices_hit=sum(Hs),bonferroni_indep=round(B,2),lattice_ratio=round(sum(Hs)/B,4),lattice_ci95=_wide(sum(Hs),Hs,Bs,nl),
pairs=sum(Cs),pairs_model=round(Q,2),kappa_hat=round(sum(Cs)/Q,4)if Q>0 else None,kappa_ci95=_wide(sum(Cs),Cs,Qs,nl)if Q>0 else None)
  K=sum(r[3][0]for r in R);s16=sum(r[3][1]for r in R)/10.0/K
  out[u]=dict(lattices=nl,classes_of_e=(S["r"]-1)//2,exhaustive=(nl==(S["r"]-1)//2),tower_ok=tower_ok(int(u)),p=S["p"],kept=K,mean_s16=round(s16,4),alpha_emp=round(er-s16,4),by_L=by)
 return out
