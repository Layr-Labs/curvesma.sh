"""Toy GLV class-rho quality measurement (own code, stdlib python3).
Curve y^2=x^3+13 over F_524347: 3|P-1 (beta), prime Q=523969, 3|Q-1
(lambda); asserts 3|P-1, 3|Q-1, point order Q (Schoof-free: add G Q
times == O on a 300-step ladder + order divides Q and Q prime), and
(beta*x,y)=[lambda](x,y) via full dlog table. Walk = canonical 6-orbit
reps, f(c)=canon(c+T_h(c)), h = bits above DP bits, ring<=8 + least-
canonical doubling escape. rho = sqrt(pi*DK/2)/E[t]. 8 tables x 2000
starts; ~10 s. Supports mechanism + collision model only, NOT
production partition/costs. Input: {"tables":8,"trials":2000}."""
import json, math, random, sys
P,B,Q = 524347,13,523969; DPB=4
cfg=json.loads(sys.argv[1]) if len(sys.argv)>1 else {}
TB,TR = cfg.get("tables",8), cfg.get("trials",2000)
def pr(n):
 if n<2:return 0
 for d in(2,3,5,7,11,13,17,19,23,29,31,37):
  if n%d==0:return n==d
 d,s=n-1,0
 while d%2==0:d//=2;s+=1
 for a in(2,3,5,7,11,13,17,19,23,29,31,37):
  x=pow(a,d,n)
  if x in(1,n-1):continue
  for _ in range(s-1):
   x=x*x%n
   if x==n-1:break
  else:return 0
 return 1
def sd(v,p):
 if v==0:return 0
 if pow(v,(p-1)//2,p)!=1:return
 y=pow(v,(p+1)//4,p)
 return min(y,p-y)
assert pr(P) and pr(Q) and P%3==1 and Q%3==1 and pow(B,(P-1)//2,P)==P-1
def r3(m):
 g=2
 while pow(g,(m-1)//3,m)==1:g+=1
 return pow(g,(m-1)//3,m)
BE,LA=r3(P),r3(Q);B2=BE*BE%P;DK=(Q-1)//6
def ad(X,Y):
 if X is None:return Y
 if Y is None:return X
 x1,y1=X;x2,y2=Y
 if x1==x2 and(y1+y2)%P==0:return
 s=(3*x1*x1*pow(2*y1,-1,P)if x1==x2 else(y2-y1)*pow((x2-x1)%P,-1,P))%P
 x3=(s*s-x1-x2)%P
 return(x3,(s*(x1-x3)-y1)%P)
ys={x:y for x in range(P)if(y:=sd((x*x%P*x+B)%P,P))is not None}
assert len(ys)*2==Q-1  # point count: every QR x gives 2 points, B nres
pts=list(ys.items());G=pts[0]
ob=lambda x:frozenset((x,x*BE%P,x*B2%P))
cls=sorted({ob(x)for x,_ in pts});assert len(cls)==DK
rep={m:(m,ys[m]if ys[m]%2==0 else(P-ys[m])%P)for m in(map(min,cls))}
XS=sorted(rep)
p2k={};k2p={};cu=None
for k in range(1,Q):cu=ad(cu,G);p2k[cu]=k;k2p[k]=cu
assert ad(cu,G) is None  # order Q
X7=k2p[7];phi0=(X7[0]*BE%P,X7[1])
if p2k.get(phi0)==LA*7%Q: pass
elif p2k.get(phi0)==LA*LA%Q*7%Q: LA=LA*LA%Q
else: raise AssertionError("endo mismatch")
assert all(p2k[(x*BE%P,y)]==LA*k%Q
 for k in random.Random(1).sample(range(1,Q),500)
 for x,y in[k2p[k]])
out=[]
for ts in range(TB):
 rr=random.Random(8007+ts)
 def mu(k,X):
  Rr,Aa=None,X
  while k:
   if k&1:Rr=ad(Rr,Aa)
   Aa=ad(Aa,Aa);k>>=1
  return Rr
 Tt=[ad(mu(rr.randrange(1,Q),G),mu(rr.randrange(1,Q),G))for _ in range(256)]
 f={xm:min(ob((ad(rep[xm],Tt[(xm>>DPB)&255])or G)[0]))for xm in XS}
 es=lambda xm:min(ob((ad(rep[xm],rep[xm])or G)[0]))
 tt=[]
 for i in range(TR):
  x=XS[random.Random(1000+ts*100003+i).randrange(DK)]
  seen,ring,t=set(),[],0
  while x not in seen:
   seen.add(x);x=f[x];t+=1
   if x in ring:x=es(min(ring+[x]));ring=[]
   ring.append(x)
   if len(ring)>8:ring.pop(0)
  tt.append(t)
 idl=math.sqrt(math.pi*DK/2);m=sum(tt)/TR
 out.append((round(idl/m,4), m))
pm=sum(mm for _,mm in out)/len(out)
oo=[r for r,_ in out]
print(json.dumps(dict(curve=dict(P=P,B=B,Q=Q,DK=DK),dp_bits=DPB,
 tables=TB,trials=TR,total_trials=TB*TR,rho_tables=oo,
 rho_pooled=round(math.sqrt(math.pi*DK/2)/pm,4),
 rho_mean=round(sum(oo)/len(oo),4),rho_worst=min(oo),
 scope="supports mechanism and collision model, not exact production partition/costs")))
