import random
P=101

def trim(a):
 while len(a)>1 and a[-1]==0:a.pop()
 return a or [0]
def add(a,b):return trim([((a[i]if i<len(a)else 0)+(b[i]if i<len(b)else 0))%P for i in range(max(len(a),len(b)))])
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]=(c[i+j]+x*y)%P
 return trim(c)
def mm(A,B):
 return [[sum_poly([mul(A[i][l],B[l][j])for l in range(len(B))])for j in range(len(B[0]))]for i in range(len(A))]
def sum_poly(a):
 v=[0]
 for x in a:v=add(v,x)
 return v
def deg(a):return len(trim(a[:]))-1 if any(a)else -10**6
def degrees(A,s):return [max(deg(v)+s[j]for j,v in enumerate(row))for row in A]
def coeff(a,j):return a[j]if 0<=j<len(a)else 0
def basis(F,s):
 k=len(F);nc=len(F[0]);order=sorted(range(k),key=lambda i:(s[i],i));U=[[int(j==i)for j in range(k)]for i in order];R=[[v[0]for v in F[i]]for i in order];piv=[];flags=[]
 for i in range(k):
  for row,c in piv:
   v=R[i][c]
   R[i]=[(x-v*y)%P for x,y in zip(R[i],R[row])];U[i]=[(x-v*y)%P for x,y in zip(U[i],U[row])]
  c=next((j for j in range(nc)if R[i][j]),None)
  if c is None:flags.append(False);continue
  inv=pow(R[i][c],-1,P);R[i]=[x*inv%P for x in R[i]];U[i]=[x*inv%P for x in U[i]];piv.append((i,c));flags.append(True)
 B=[[[0,v]if flags[i]else[v]for v in row]for i,row in enumerate(U)]
 return B,degrees(B,s)
def pm(F,sigma,s):
 if sigma==1:return basis(F,s)
 h=sigma//2;B,t=pm([[v[:h]or[0]for v in row]for row in F],h,s);G=mm(B,F);G=[[v[h:sigma]or[0]for v in row]for row in G];C,u=pm(G,sigma-h,t);return mm(C,B),u

def rank(A):
 if not A:return 0
 A=[r[:]for r in A];i=0
 for j in range(len(A[0])):
  z=next((z for z in range(i,len(A))if A[z][j]%P),None)
  if z is None:continue
  A[i],A[z]=A[z],A[i];iv=pow(A[i][j]%P,-1,P);A[i]=[x*iv%P for x in A[i]]
  for z in range(i+1,len(A)):
   f=A[z][j]%P;A[z]=[(x-f*y)%P for x,y in zip(A[z],A[i])]
  i+=1
  if i==len(A):break
 return i

def check(seed=1):
 r=random.Random(seed);runs=0;max_order=12
 for k in range(2,7):
  for nc in range(1,k):
   for sig in (1,2,3,4,7,12):
    s=[r.randrange(4)for _ in range(k)];F=[[[r.randrange(P)for _ in range(sig)]for _ in range(nc)]for _ in range(k)];B,u=pm(F,sig,s);BF=mm(B,F)
    assert all(not any(v[:sig])for row in BF for v in row)
    assert u==degrees(B,s)
    lead=[[coeff(B[i][j],u[i]-s[j])for j in range(k)]for i in range(k)]
    assert rank(lead)==k
    constraints=[[coeff(F[i][j],t-v)for i in range(k)for v in range(sig)]for j in range(nc)for t in range(sig)]
    assert sum(u)-sum(s)==rank(constraints)
    runs+=1
 return dict(ok=True,cases=runs,field=P,largest_order=max_order,largest_rows=6,checks='annihilation, shifted degrees, leading rank, module index')

def nullvector(A):
 A=[row[:]for row in A];nr=len(A);nc=len(A[0]);piv=[];r=0
 for c in range(nc):
  j=next((j for j in range(r,nr)if A[j][c]%P),None)
  if j is None:continue
  A[r],A[j]=A[j],A[r];iv=pow(A[r][c]%P,-1,P);A[r]=[v*iv%P for v in A[r]]
  for j in range(nr):
   if j!=r:
    z=A[j][c]%P;A[j]=[(x-z*y)%P for x,y in zip(A[j],A[r])]
  piv.append(c);r+=1
  if r==nr:break
 free=next((c for c in range(nc)if c not in piv),None)
 if free is None:return None
 v=[0]*nc;v[free]=1
 for i,c in reversed(list(enumerate(piv))):v[c]=-sum(A[i][j]*v[j]for j in range(c+1,nc))%P
 return v

def bw_check(seed=1):
 global P
 old=P;P=65537;r=random.Random(seed);success=0;tried=40
 for case in range(tried):
  d=8;m=5;n=3;L=(d+m-1)//m+(d+n-1)//n+2
  X=[[r.randrange(P)for _ in range(m)]for _ in range(d)]
  if case>=20:
   ix=list(range(d));r.shuffle(ix);X=[[0]*m for _ in range(d)]
   for z,h in enumerate(ix):X[h][z%m]=r.getrandbits(128)%P
  Y=[[r.randrange(P)for _ in range(n)]for _ in range(d)]
  S=[[[sum(X[h][i]*pow(h,t,P)*Y[h][j]for h in range(d))%P for j in range(n)]for i in range(m)]for t in range(L)]
  F=[[[S[t][j][i]for t in range(L)]for j in range(m)]for i in range(n)]+[[[int(i==j)]for j in range(m)]for i in range(m)]
  B,u=pm(F,L,[0]*n+[1]*m);ix=sorted(range(m+n),key=lambda i:(u[i],i))[:n];C=[[coeff(B[i][j],u[i])for j in range(n)]for i in ix];c=nullvector([list(row)for row in zip(*C)])
  if c is None:continue
  v=[0]*d
  for z,i in zip(c,ix):
   D=u[i]
   for h in range(d):
    rel=sum(pow(h,D-a,P)*Y[h][j]*coeff(B[i][j],a)for j in range(n)for a in range(D+1))%P
    assert rel==0
    v[h]=(v[h]+z*sum(pow(h,D-a-1,P)*Y[h][j]*coeff(B[i][j],a)for j in range(n)for a in range(D)))%P
  assert all(h*v[h]%P==0 for h in range(d))
  success+=any(v)
 P=old
 assert success>=36
 return dict(cases=tried,success=success,field=65537,dimension=8,block=[5,3],dense_cases=20,block_diagonal_cases=20)

def packing_check(seed=1):
 r=random.Random(seed);mod=int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001',16);base=1<<576
 for n in (1,2,3,8,32,128):
  a=[r.randrange(mod)for _ in range(n)];b=[r.randrange(mod)for _ in range(n)];u=sum(x<<(576*i)for i,x in enumerate(a));v=sum(x<<(576*i)for i,x in enumerate(b));prod=u*v
  got=[((prod>>(576*i))&(base-1))%mod for i in range(2*n-1)];ref=[sum(a[j]*b[i-j]for j in range(n)if 0<=i-j<n)%mod for i in range(2*n-1)]
  assert got==ref
 return dict(ok=True,cases=6,slot_bits=576,largest_degree=127,scope='algebraic packing, not instructions')

if __name__=='__main__':
 import json,sys
 x=json.load(sys.stdin);seed=x.get('seed',1);print(json.dumps(dict(pm=check(seed),bw=bw_check(seed),packing=packing_check(seed))))
