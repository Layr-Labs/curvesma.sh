import json, math, sys
from fractions import Fraction as Q
from pathlib import Path
from itertools import product

H=[-1,0,-2,-2,2,0,1]; F=[1,1,0,2,0,-2,1]
R=52435875175126190479447740508185965837690552500527637822603658699938581184513
U=15132376222941642752
CERTIFICATE=json.loads('{"basis":[[2,0,0,0,0,0],[-1,3,-2,1,-1,0],[1,-5,2,-1,1,0],[0,-4,-4,4,0,2],[3,0,3,-3,0,-1],[0,1,-5,2,-1,1]],"factors":{"11":[[[[6,3,1],1],[[10,3,1],1],[[2,5,1],1]],[[[8,10,1],1],[[7,1,2,10,1],1]]],"13":[[[[12,0,11,11,2,0,1],1]],[[[4,3,12,1],1],[[10,9,12,1],1]]],"17":[[[[7,8,6,1],1],[[12,13,11,1],1]],[[[1,1,0,2,0,15,1],1]]],"23":[[[[6,7,5,1],1],[[19,20,18,1],1]],[[[1,1,0,2,0,21,1],1]]],"29":[[[[28,0,27,27,2,0,1],1]],[[[19,1,1],1],[[26,17,13,26,1],1]]],"3":[[[[2,0,1,1,2,0,1],1]],[[[2,1],2],[[1,0,1],2]]],"31":[[[[17,1],1],[[20,1],1],[[16,6,1],1],[[2,19,1],1]],[[[20,1],1],[[22,1],1],[[4,8,1],1],[[22,10,1],1]]],"37":[[[[36,19,1],1],[[36,23,1],1],[[36,32,1],1]],[[[2,1],1],[[32,34,1],1],[[11,10,36,1],1]]],"41":[[[[11,1],1],[[26,1],1],[[6,13,1],1],[[7,32,1],1]],[[[1,1,0,2,0,39,1],1]]],"43":[[[[42,4,1],1],[[38,11,1],1],[[17,28,1],1]],[[[40,1],1],[[12,2,1],1],[[37,36,42,1],1]]],"47":[[[[9,1],1],[[15,1],1],[[25,1],1],[[26,1],1],[[30,1],1],[[36,1],1]],[[[1,1,0,2,0,45,1],1]]],"5":[[[[2,1,1],1],[[3,2,1],1],[[4,2,1],1]],[[[3,2,1],1],[[2,4,0,1,1],1]]],"53":[[[[52,0,51,51,2,0,1],1]],[[[1,1,0,2,0,51,1],1]]],"59":[[[[58,0,57,57,2,0,1],1]],[[[27,57,1],1],[[35,7,32,0,1],1]]],"61":[[[[20,32,1],1],[[58,35,1],1],[[60,55,1],1]],[[[3,1],1],[[18,1],1],[[49,1],1],[[54,1],1],[[25,57,1],1]]],"67":[[[[66,0,65,65,2,0,1],1]],[[[12,1],1],[[23,1],1],[[27,1],2],[[45,43,1],1]]],"7":[[[[4,5,3,1],1],[[5,6,4,1],1]],[[[4,1],1],[[3,2,1],1],[[3,2,6,1],1]]],"71":[[[[36,1],1],[[69,1],1],[[26,47,1],1],[[41,61,1],1]],[[[2,43,1],1],[[36,43,16,26,1],1]]],"73":[[[[33,34,32,1],1],[[42,43,41,1],1]],[[[17,1],1],[[42,1],1],[[45,30,1],1],[[22,55,1],1]]],"79":[[[[66,1],1],[[73,1],1],[[60,22,1],1],[[54,76,1],1]],[[[24,23,78,1],1],[[56,55,78,1],1]]],"83":[[[[48,3,1],1],[[82,23,1],1],[[64,57,1],1]],[[[1,1,0,2,0,81,1],1]]],"89":[[[[57,1],1],[[64,1],1],[[23,7,1],1],[[31,50,1],1]],[[[7,34,1],1],[[51,83,60,53,1],1]]],"97":[[[[27,1],1],[[79,1],1],[[60,4,1],1],[[76,84,1],1]],[[[41,1],1],[[72,1],1],[[80,1],1],[[36,35,96,1],1]]]},"special_charges":{"19":[72,1],"2":[5,6]},"witness41":[[[[11,1],1],[[26,1],1],[[6,13,1],1],[[7,32,1],1]],[[[1,1,0,2,0,39,1],1]]]}')

def trim(a):
 while a and a[-1]==0:a.pop()
 return a
def mul(a,b,p=0):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return trim([v%p for v in c] if p else c)
def sub(a,b,p):return trim([((a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0))%p for i in range(max(len(a),len(b)))])
def rem(a,b,p):
 a=[x%p for x in a];b=trim([x%p for x in b]);inv=pow(b[-1],-1,p)
 while len(a)>=len(b):
  c=a[-1]*inv%p;k=len(a)-len(b)
  for i,x in enumerate(b):a[i+k]=(a[i+k]-c*x)%p
  trim(a)
 return a
def gcd(a,b,p):
 while b:a,b=b,rem(a,b,p)
 return [x*pow(a[-1],-1,p)%p for x in a] if a else []
def power(a,n,f,p):
 r=[1]
 while n:
  if n&1:r=rem(mul(r,a,p),f,p)
  a=rem(mul(a,a,p),f,p);n//=2
 return r
def irreducible(f,p):
 n=len(f)-1
 if n<1 or f[-1]!=1:return False
 x=rem([0,1],f,p);t=x
 for k in range(1,n+1):
  t=power(t,p,f,p)
  if k<=n//2 and n%k==0 and len(gcd(sub(t,x,p),f,p))>1:return False
 return t==x
def factors(a,fs,p):
 v=[1];seen=[]
 for f,e in fs:
  assert isinstance(e,int) and e>0 and all(isinstance(x,int) and 0<=x<p for x in f)
  assert irreducible(f,p) and f not in seen
  seen.append(f)
  for _ in range(e):v=mul(v,f,p)
 assert v==[x%p for x in a]

def det(a):
 a=[[Q(x) for x in row] for row in a];d=Q(1);n=len(a)
 for j in range(n):
  k=next((i for i in range(j,n) if a[i][j]),None)
  if k is None:return Q(0)
  if j!=k:a[j],a[k]=a[k],a[j];d=-d
  v=a[j][j];d*=v
  for i in range(j+1,n):
   z=a[i][j]/v
   for k in range(j+1,n):a[i][k]-=z*a[j][k]
 return d
def solve(a,b):
 a=[[Q(x) for x in row]+[Q(y)] for row,y in zip(a,b)];n=len(a)
 for j in range(n):
  k=next(i for i in range(j,n) if a[i][j]);a[j],a[k]=a[k],a[j];v=a[j][j];a[j]=[x/v for x in a[j]]
  for i in range(n):
   if i!=j:
    v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[j])]
 return [row[-1] for row in a]
def kmul(a,b):
 c=mul(a,b)
 for i in range(len(c)-1,5,-1):
  v=c[i]
  for j in range(6):c[i-6+j]-=v*H[j]
 return (c[:6]+[0]*6)[:6]
def matmul(a,b):return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def resultant(a,b):
 m=len(a)-1;n=len(b)-1
 return det([[0]*i+a+[0]*(n-1-i) for i in range(n)]+[[0]*i+b+[0]*(m-1-i) for i in range(m)])
def exp_lo(x,n):return sum((x**k/Q(math.factorial(k)) for k in range(n+1)),Q(0))
def exp_hi(x,n):return exp_lo(x,n)+x**(n+1)/math.factorial(n+1)/(1-x/Q(n+2))

def main():
 assert json.load(sys.stdin)=={}
 data=CERTIFICATE
 basis=data['basis'];assert len(basis)==6 and all(len(row)==6 for row in basis)
 B=[[Q(x,2) for x in row] for row in basis];bt=[list(v) for v in zip(*B)]
 assert abs(det(B))==Q(1,4)
 def coords(v):return solve(bt,v)
 def integr(v):return all(x.denominator==1 for x in coords(v))
 def matrix(v):return [list(row) for row in zip(*(kmul(v,b) for b in B))]
 def norm(v):return det(matrix(v))/det(bt)
 assert all(integr(kmul(a,b)) for a in B for b in B)
 assert all(integr([int(i==j) for i in range(6)]) for j in range(6))
 assert integr([1,1,0,0,0,0])
 traces=[]
 for b in B:
  M=[coords(kmul(b,c)) for c in B];traces.append(sum(M[i][i] for i in range(6)))
 gram=[[sum(coords(kmul(a,b))[k]*traces[k] for k in range(6)) for b in B] for a in B]
 assert abs(det(gram))==184832
 assert abs(resultant(H,[i*H[i] for i in range(1,7)]))==2957312
 assert resultant(F,[i*F[i] for i in range(1,7)])==34371
 w=[1,1,0,0,0,0];pi=[Q(x,2) for x in [-1,-5,2,-1,1,0]]
 assert norm(w)==2 and abs(norm(pi))==4
 def divide(a,b):return solve(matrix(b),a)
 def member(a,b):return all(x.denominator==1 for x in divide(a,b))
 wp=kmul(kmul(w,w),kmul(pi,pi));u=divide([2,0,0,0,0,0],wp)
 # divide returns coordinates in B, convert back to the power basis.
 unit=[sum(u[i]*B[i][j] for i in range(6)) for j in range(6)]
 assert all(x.denominator==1 for x in u) and kmul(unit,wp)==[2,0,0,0,0,0] and abs(norm(unit))==1
 assert not member(pi,w) and not member(w,pi)
 one=[1,0,0,0,0,0]
 assert member([2,0,0,0,0,0],w) and member([2,0,0,0,0,0],pi)
 root=None
 for bits in product(range(2),repeat=6):
  t=[sum(bits[i]*B[i][j] for i in range(6)) for j in range(6)]
  v=kmul(t,t);v=[v[i]+t[i]+one[i] for i in range(6)]
  if member(v,pi):root=t;break
 assert root is not None and not member(root,pi) and not member([root[i]+one[i] for i in range(6)],pi)
 assert irreducible([x%2 for x in F],2)
 primes=[p for p in range(3,100) if all(p%d for d in range(2,math.isqrt(p)+1))]
 fs=data['factors'];assert set(fs)=={str(p) for p in primes if p!=19}
 for p in primes:
  if p!=19:
   h,f=fs[str(p)];factors(H,h,p);factors(F,f,p)
 assert len(fs['3'][0])==1 and len(fs['3'][0][0][0])==7 and fs['3'][0][0][1]==1
 # Minkowski <15; only the already-principal primes over2 have norm below15.
 assert 184832<431**2 and Q(40,31)**3*Q(math.factorial(6),6**6)*431<15
 assert all(p**(len(h)-1)>=15 for p in (3,5,7,11,13) for h,e in fs[str(p)][0])
 h41,f41=data['witness41'];factors(H,h41,41);factors(F,f41,41)
 assert [[11,1],1] in h41 and len(f41)==1 and len(f41[0][0])==7 and f41[0][1]==1
 # Dedekind criterion at 19, using monic lifts of the exact factorization.
 # Generate its factors by trial monic division, avoiding reliance on a saved index flag.
 def dividepoly(a,b,p):
  a=[x%p for x in a];q=[0]*max(0,len(a)-len(b)+1)
  while len(a)>=len(b):
   z=a[-1];k=len(a)-len(b);q[k]=z
   for i,x in enumerate(b):a[i+k]=(a[i+k]-z*x)%p
   trim(a)
  return trim(q),a
 a=[x%19 for x in H];hs=[]
 for degree in range(1,4):
  for co in product(range(19),repeat=degree):
   g=list(co)+[1]
   if not irreducible(g,19):continue
   e=0
   while len(a)>=len(g):
    q,rm=dividepoly(a,g,19)
    if rm:break
    a=q;e+=1
   if e:hs.append((g,e))
  if a==[1]:break
 assert a==[1]
 lift=[1]
 for g,e in hs:
  for _ in range(e):lift=mul(lift,g)
 dif=[(H[i]-lift[i])//19 for i in range(7)]
 assert all((H[i]-lift[i])%19==0 for i in range(7))
 assert all(len(gcd(dif,g,19))==1 for g,e in hs if e>1)
 # F mod 3 and the two translated Eisenstein clusters.
 assert mul(mul([2,1],[2,1],3),mul([1,0,1],[1,0,1],3),3)==[x%3 for x in F]
 der=[i*F[i] for i in range(1,7)]
 assert sum(F)==3 and sum(der)%3==0
 assert rem([x%3 for x in der],[1,0,1],3)==[]
 assert rem([x%9 for x in F],[1,0,1],3)==[]
 # F(i)=-3i by exact reduction over Z[i].
 fi=[0,0]
 for i,c in enumerate(F):fi[i%2]+=c*(-1)**(i//2)
 assert fi==[0,-3]
 assert 34371==3**3*19*67
 charges=data['special_charges']
 assert charges=={'2':[5,6],'19':[72,1]}
 partial=sum((Q(n,p**d-1) for p,(n,d) in ((2,charges['2']),(19,charges['19']))),Q(0))
 for p in primes:
  if p==19:continue
  hf,ff=fs[str(p)]
  local_degree=0
  for h,e in hf:
   a=len(h)-1;nw=(-1)**a*sum(c*(-1)**i for i,c in enumerate(h))%p
   for f,ef in ff:
    b=len(f)-1;L=math.lcm(a,b);g=math.gcd(a,b)
    v=pow(nw,L//a,p)*pow((-1)**b*f[0]%p,L//b,p)%p
    char=pow(v,(p-1)//2,p);assert char in (1,p-1)
    local_degree+=2*g*L
    partial+=Q(2*g,p**L-1) if char==1 else Q(g,p**(2*L)-1)
  assert local_degree<=72
 assert partial<6
 assert exp_hi(Q(23,5),50)<100 and exp_lo(Q(69,50),30)>Q(250,63) and exp_hi(Q(69,100),30)<2
 D=(184832*34371)**12*2**78
 assert Q(40,31)**36*Q(math.factorial(72),72**72)*(math.isqrt(D)+1)<2**149
 assert Q(12)**2*(2-Q(97,56))>Q(31,5)**2 and Q(97,56)**2>3
 tail=72*(Q(125506,100000)*Q(13,10)*Q(21,115)+Q(1,99))
 bound=Q(149*13,10)+(6+tail)/Q(69,100)
 assert bound<Q(23454,100)<240<254 and 2**254<R
 Dg=184832**2*4**6*(U//2**16)**6*2
 gb=(Dg.bit_length()+1)//2
 assert U%2**16==0 and Q(40,31)**6*Q(math.factorial(12),12**12)<1 and 2**gb*(1+gb)**11<R
 assert Q(40,31)**3*Q(math.factorial(6),6**6)*431<512 and 512*513**5<R
 print(json.dumps({'status':'ok','partial_lt_6':True,'ef_log2_upper':str(bound),'ef_log2_lt':240,'k_disc':184832,'k_class_number':1,'field_degree':72,'eg_minkowski_bits':gb,'primes_below_100':len(primes)+1},sort_keys=True,separators=(',',':')))
main()
