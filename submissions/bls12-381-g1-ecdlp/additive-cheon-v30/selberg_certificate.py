#!/usr/bin/env python3
"""Exact finite main/exception certificate only. NO analytic remainder bound."""
from math import isqrt, gcd, comb
from fractions import Fraction as F
from pathlib import Path
import itertools, json

def factors(n):
 p=2; out=[]
 while p*p<=n:
  if n%p==0:
   out.append(p)
   while n%p==0:n//=p
  p+=1
 if n>1:out.append(n)
 return out
def sf(n):return all(n%(p*p) for p in factors(n))
def mob(n):return 0 if not sf(n) else (-1 if len(factors(n))%2 else 1)
def g(n):
 ans=F(1)
 for p in factors(n):ans*=F(2*p-1,p*p)
 return ans
def h(n):
 if not sf(n):return F(0)
 ans=F(1)
 for p in factors(n):ans*=F(2*p-1,(p-1)**2)
 return ans

weight_checks=[]
for Rtest in (2,4,8,16,32):
 G=sum((h(n) for n in range(1,Rtest+1)),F())
 lam={}
 for v in range(1,Rtest+1):
  if not sf(v):lam[v]=F();continue
  fac=F(1)
  for p in factors(v):fac*=F(p*p,(p-1)**2)
  tail=sum((h(k) for k in range(1,Rtest//v+1) if gcd(v,k)==1),F())
  lam[v]=mob(v)*fac*tail/G
 assert lam[1]==1
 Q=sum((lam[v]*lam[w]*g(v*w//gcd(v,w))
        for v in lam for w in lam),F())
 assert Q==1/G
 weight_checks.append(dict(R=Rtest,G=str(G),quadratic_form=str(Q)))

# Exhaustive verification of the two graph inequalities used in the proof.
# A full nonreal-U graph on <=3 representatives is contained in K3.
# At most two representatives are odd. All even vertices are marked.
graph_checks=0
V=range(3)
for oddmask in range(8):
 if oddmask.bit_count()>2:continue
 evenmask=7^oddmask
 opts=[]
 for supp in range(8):
  for marked in range(8):
   if marked&~supp:continue
   if (supp&evenmask)&~marked:continue
   opts.append((supp,marked))
 for supp,marked in opts:
  ev=(supp&evenmask).bit_count()
  edges=sum(bool(supp>>i&1 and supp>>j&1 and
                 ((marked>>i&1) or (marked>>j&1)))
            for i in V for j in V if i!=j)
  assert edges<=2*marked.bit_count()+2*ev
  for supp2,marked2 in opts:
   ev2=(supp2&evenmask).bit_count()
   cross=sum(bool(supp>>i&1 and supp2>>j&1 and
                  ((marked>>i&1) or (marked2>>j&1)))
             for i in V for j in V if i!=j)
   assert cross<=marked.bit_count()+marked2.bit_count()+ev+ev2
   graph_checks+=1

# Toy exact prefix-star audit over every pair of two-row prefix choices.
lists=((2,3,5),(7,11))
options=[(c,tuple(L[:t])) for c,L in enumerate(lists) for t in range(len(L)+1)]
star_checks=0
for rowsB in itertools.product(options,repeat=2):
 WB=[sum(len(P) for c,P in rowsB if c==j) for j in range(2)]
 RB=[sum(c==j for c,P in rowsB) for j in range(2)]
 b={p:sum(p in P for c,P in rowsB) for L in lists for p in L}
 nB=sum(b.values()); Bsquare=sum(v*v for v in b.values())
 assert Bsquare<=sum(WB[j]*RB[j] for j in range(2))
 for rowsG in itertools.product(options,repeat=2):
  WG=[sum(len(P) for c,P in rowsG if c==j) for j in range(2)]
  RG=[sum(c==j for c,P in rowsG) for j in range(2)]
  gg={p:sum(p in P for c,P in rowsG) for L in lists for p in L}
  Gsquare=sum(v*v for v in gg.values())
  assert Gsquare<=sum(WG[j]*RG[j] for j in range(2))
  cross=sum(b[p]*gg[p] for p in b)
  assert cross*cross<=Bsquare*Gsquare
  star_checks+=1

# Small exact audit of the dyadic-binomial divisibility argument.
def primes_to(v):
 sieve=bytearray(b"\x01")*(v+1);sieve[:2]=b"\x00\x00"
 for p in range(2,isqrt(v)+1):
  if sieve[p]:sieve[p*p:v+1:p]=b"\x00"*(((v-p*p)//p)+1)
 return [p for p in range(2,v+1) if sieve[p]]
primes=primes_to(4096)
for j in range(2,13):
 mid=1<<(j-1);end=1<<j
 ps=[p for p in primes if mid<p<=end]
 C=comb(end,mid)
 for p in ps:assert C%p==0
 assert C < 1<<end
 assert len(ps)<=(end-1)//(j-1)

# Production inputs.
r=52435875175126190479447740508185965837690552500527637822603658699938581184513
d=1<<20;m=(r-1)//d
X1=1<<82;X2=1<<70;Z=X1*X2
n=67704887904593045852627984777072135
M=5893133993817406318491
L=228988810152649578064853576960394133503
assert L*L+L+1==r
assert (r-1)%d==0 and m%6==0
assert 120*(n-1)**2<11*m<=120*n*n
NA=(196*n*d)//X2+2
assert (NA+1)//2==M
assert X2<X1//2 and 3*X1*X1*X2<r
baseline=F(36*n**4,m)

# Elementary finite G lower certificate; all assertions use integers/Fractions.
R=1<<35;J=isqrt(R)
bits=128;B=1<<bits;terms=24
z=B//3;zz=z*z//B;power=z;LN2=0
for j in range(terms):
 LN2+=(2*power)//(2*j+1)
 power=power*zz//B
def log_ratio_lower_scaled(num,den):
 assert num>=den
 k=0
 while num>=2*den:den*=2;k+=1
 z=((num-den)*B)//(num+den)
 zz=z*z//B;power=z;value=0
 for j in range(terms):
  value+=(2*power)//(2*j+1)
  power=power*zz//B
 return k*LN2+value
Hlo=sum(B//j for j in range(1,J+1))
Hhi=Hlo+J
minus=(B+2*J-1)//(2*J)
total=0
for u in range(1,J+1):
 v=R//u
 Hnlo=Hlo+log_ratio_lower_scaled(v,J)-minus+B//(2*v)
 total+=Hnlo//u
squarehi=(Hhi*Hhi+B-1)//B
Dlo=F(2*total-squarehi,B)
assert J==185363
assert Dlo==F(109830543788030032342585375633314908186445,
              340282366920938463463374607431768211456)

# Rational P1 max-class upper. The endpoint sum is 0.57+0.64=1.21.
C=1<<128
q=isqrt(500*(d-1)*C*C)+1
assert q*q>500*(d-1)*C*C
wup=F(C+q,d*C)
area=F(24*n*M*121**2*Z*Z,10000*r)*wup
A=area/baseline
main=A/Dlo
S35=1+sum(((1<<j)-1)//(j-1) for j in range(2,36))
assert S35==2086667226
exception=F(48*n*M*(S35+1),1)/baseline
assert main<F(3897665922609721,10**17)
assert exception<F(264184166416664,10**17)
