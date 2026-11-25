"""exact_levels.py -- test of the exact-integer level body of KERNELS.md section 1.3 (stdlib, ~0.7 s): a transcription
of v20_lat_b2.py full_tree in which every level with B_k >= T runs the exact body, on synthetic 12-D Gram matrices with
such levels at k >= 3 (B_0..B_2 < T), compared with exact-rational enumeration of {x != 0 : Q(x) <= T} (one of +-x):
missing points and the largest Q/T - 1 emitted.  Negative control: fixed-point body at every level (b2 clipped)."""
import math,random,json,sys,time
from fractions import Fraction as F
Q32,Q38,M64=1<<32,1<<38,(1<<64) -1
def lut4(r2):
 r=r2-1
 if r==0:
  return math.isqrt(1<<52) +1
 c=64-r.bit_length();top= (r<<c) >>53
 hi= ((top+1) << (53-c)) -1 if c<=53 else top>> (c-53)
 return math.isqrt((hi+1) <<52) +1
def gs(G):
 n=len(G);D= [1] + [0] *n;lam= [[0] *n for _ in range(n)]
 for i in range(n):
  for j in range(i+1):
   u=G[i][j]
   for k in range(j):
    u= (D[k+1] *u-lam[i][k] *lam[j][k]) //D[k]
   if j<i:
    lam[i][j] =u
   else:
    D[i+1] =u
 return D,lam
def tree(G,T,exact):
 n=12;D,lam=gs(G)
 mu= [[(lam[j][k] <<32) //D[k+1] if j>k else 0 for k in range(n)] for j in range(n)]
 b2= [min(M64,(D[k+1] <<64) // (D[k] *T)) for k in range(n)]
 inv= [((D[k] *T<<38) +D[k+1] -1) //D[k+1] for k in range(n)]
 Xm= [0] *n
 for k in range(n-1,-1,-1):
  Xm[k] =math.isqrt((D[k] *T+D[k+1] -1) //D[k+1]) +1+sum((Xm[j] +1) //2 for j in range(k+1,n))
 ec= [sum(Xm[j] for j in range(k+1,n)) +1 for k in range(n)]
 A1= (D[2] <<38) // (D[1] *T)
 exc= [exact and D[k+1] >=T*D[k] for k in range(n)]
 xs= [0] *n;out= [];st=dict(exact_nodes=0)
 def node(k,rem,az):
  ck= -sum(mu[j][k] *xs[j] for j in range(k+1,n))
  if exc[k]:
   st["exact_nodes"] +=1
   S=sum(lam[j][k] *xs[j] for j in range(k+1,n));W=rem*T*D[k] *D[k+1];den=T*D[k] *D[k+1]
   f= (-S) //D[k+1]
   for x in range(max(f-1,0) if az else f-1,f+3):
    N=D[k+1] *x+S;t=Q38*N*N
    if t<=W:
     xs[k] =x;node(k-1,-((t-W) //den),az and x==0)
   xs[k] =0; return
  r2= ((rem*inv[k]) >>64) +1;rad=lut4(r2);e=ec[k]
  lo= (1 if k==0 else 0) if az else (ck-e-rad+Q32-1) >>32
  hi= (ck+e+rad) >>32
  if k==0:
   for x in range(lo,hi+1):
    xs[0] =x;out.append(tuple(xs))
  elif k==1:
   R1=hi-lo+1;d19=abs((lo<<32) -ck) >>13
   rem0=rem+ ((2*R1*e*A1) >>32) +2+937* (hi-lo) +1- ((d19*d19*b2[1]) >>64)
   dr1= (((2*lo+1) *Q32-2*ck) *A1) >>32
   for x in range(lo,hi+1):
    if rem0>=0:
     xs[1] =x;node(0,rem0,az and x==0)
    rem0-=dr1;dr1+=2*A1
  else:
   for x in range(lo,hi+1):
    ad=abs((x<<32) -ck);d19= (ad-e if ad>e else 0) >>13
    sub= (d19*d19*b2[k]) >>64
    if sub<=rem:
     xs[k] =x;node(k-1,rem-sub,az and x==0)
  xs[k] =0
 node(n-1,Q38,True)
 return out,st,D
def brute(G,T):
 n=12;D,lam=gs(G);B= [F(D[k+1],D[k]) for k in range(n)]
 mu= [[F(lam[j][k],D[k+1]) for k in range(n)] for j in range(n)];out= [];xs= [0] *n
 def rec(k,rem):
  c= -sum(mu[j][k] *xs[j] for j in range(k+1,n));s=math.isqrt(math.floor(rem/B[k])) +1
  for x in range(math.floor(c) -s-1,math.ceil(c) +s+2):
   d= (x-c) **2*B[k]
   if d<=rem:
    xs[k] =x
    if k:rec(k-1,rem-d)
    else:out.append(tuple(xs))
  xs[k] =0
 rec(n-1,F(T)); return out
def canon(x):
 for v in x:
  if v: return x if v>0 else tuple(-t for t in x)
 return x
PROFS= [[1/30,1/25,1/28,1.3,1.0,.8,2.5,1.9,1.5,1.15,.9,3.0],[1/20,1/22,1/18,1.0,.78,2.2,1.7,1.3,1.0,4.0,3.0,2.3],
[1/40,1/35,1/30,.9,1.6,1.25,.95,.75,1.8,1.4,1.05,2.6],[1/12,1/15,1/10,2.9,2.2,1.7,1.3,1.0,.76,.6,.5,.45],
[1/60,1/50,1/45,1.01,.76,.6,.5,1.2,.9,3.5,2.6,2.0]]
def main(runs=16,seed=7):
 t0=time.time();rng=random.Random(seed)
 r=dict(runs=0,ref_points=0,exact_levels=0,exact_nodes=0,missing=0,extra=0,zero_or_dup=0,max_over=-1.0,
ctrl_runs_over_0024=0,ctrl_max_over=-1.0)
 for i in range(runs):
  prof=PROFS[i%len(PROFS)];T=rng.randrange(1<<23,1<<25)
  a= [max(1,math.isqrt(int(f*T))) for f in prof]
  L= [[(a[i2] if j==i2 else rng.randint(-(a[j] //2),a[j] //2) if j<i2 else 0) for j in range(12)] for i2 in range(12)]
  G= [[sum(L[u][t] *L[v][t] for t in range(12)) for v in range(12)] for u in range(12)]
  ref=set(canon(x) for x in brute(G,T) if any(x))
  lv,st,D=tree(G,T,True);got=set(lv)
  Qf= lambda x:sum(G[u][v] *x[u] *x[v] for u in range(12) for v in range(12))
  r["runs"] +=1;r["ref_points"] +=len(ref);r["exact_nodes"] +=st["exact_nodes"]
  r["exact_levels"] +=sum(D[k+1] >=T*D[k] for k in range(3,12))
  r["missing"] +=len(ref-set(canon(x) for x in got));r["extra"] +=len(set(canon(x) for x in got) -ref)
  r["zero_or_dup"] +=sum(not any(x) for x in lv) +len(lv) -len(set(canon(x) for x in lv))
  r["max_over"] =max(r["max_over"],max((Qf(x) /T-1 for x in lv),default=-1.0))
  lc,_,_=tree(G,T,False);oc=max((Qf(x) /T-1 for x in lc),default=-1.0)
  r["ctrl_runs_over_0024"] +=oc>0.0024;r["ctrl_max_over"] =max(r["ctrl_max_over"],oc)
 r["ok"] = (r["missing"] ==0 and r["zero_or_dup"] ==0 and r["max_over"] <=0.0024 and r["exact_levels"] >0
               and r["ctrl_runs_over_0024"] >0)
 r["seconds"] =round(time.time() -t0,3)
 return r
if __name__=="__main__":
 o=main();print(json.dumps(o));sys.exit(0 if o["ok"] else 1)
