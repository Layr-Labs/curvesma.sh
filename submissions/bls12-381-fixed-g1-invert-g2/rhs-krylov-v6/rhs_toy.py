import json,random,sys
from pathlib import Path
_root=Path(__file__).parent
sys.path.insert(0,str(_root/"files" if(_root/"files").is_dir()else _root))
import pm_toy as PM

def mv(A,v):return [sum(x*y for x,y in zip(row,v))%PM.P for row in A]
def transpose(A):return[list(r)for r in zip(*A)]
def imul(A,B):return[[sum(x*y for x,y in zip(row,col))%PM.P for col in zip(*B)]for row in A]
def nullbasis(A):
 A=[row[:]for row in A];nr=len(A);nc=len(A[0]);piv=[];r=0
 for c in range(nc):
  j=next((j for j in range(r,nr)if A[j][c]%PM.P),None)
  if j is None:continue
  A[r],A[j]=A[j],A[r];iv=pow(A[r][c]%PM.P,-1,PM.P);A[r]=[x*iv%PM.P for x in A[r]]
  for j in range(nr):
   if j!=r:
    z=A[j][c]%PM.P;A[j]=[(x-z*y)%PM.P for x,y in zip(A[j],A[r])]
  piv.append(c);r+=1
  if r==nr:break
 out=[]
 for free in range(nc):
  if free in piv:continue
  v=[0]*nc;v[free]=1
  for i,c in reversed(list(enumerate(piv))):v[c]=-sum(A[i][j]*v[j]for j in range(c+1,nc))%PM.P
  out.append(v)
 return out

def trial(r,zero_s=False,zero_ideal=False):
 N=6;D=N+10;m=n=12;L=(D+m-1)//m+(D+n-2)//(n-1)+2;rand=lambda a,b:[[r.randrange(PM.P)for _ in range(b)]for _ in range(a)]
 if zero_s:
  Q=rand(D,N-1);v=[r.randrange(PM.P)for _ in range(N-1)]+[1];M=[row+[-sum(x*y for x,y in zip(row,v))%PM.P]for row in Q];S=rand(D,10);want=v+[0]*10
 else:
  M=rand(D,N);Q=rand(D,9);v=([0]*N if zero_ideal else[r.randrange(PM.P)for _ in range(N)]);s=[r.randrange(PM.P)for _ in range(9)]+[1];Mv=mv(M,v);S=[row+[(-x-sum(a*b for a,b in zip(row,s)))%PM.P]for row,x in zip(Q,Mv)];want=v+s
 A=[u+v for u,v in zip(M,S)];assert not any(mv(A,want));assert PM.rank(A)==D-1
 alpha=pow(2,-320,PM.P);assert alpha!=1
 B=[[alpha*x%PM.P for x in row]+[0]*10 for row in M];Y=[[alpha*x%PM.P for x in s]+z for s,z in zip(S,rand(D,n-10))];X=[[0]*m for _ in range(D)]
 for h in range(10):X[N+h][h]=r.getrandbits(128)%PM.P
 ix=list(range(N))
 for j,h in enumerate(ix):X[h][10+j%(m-10)]=r.getrandbits(128)%PM.P
 seq=[];powers=[];cur=Y
 for t in range(L):powers.append(cur);seq.append(imul(transpose(X),cur));cur=imul(B,cur)
 pp=L//2;qq=L-pp;K=[sum((powers[t][i]for t in range(pp)),[])for i in range(D)];O=[];curO=transpose(X)
 for t in range(qq):O+=curO;curO=imul(curO,B)
 assert PM.rank(K)==PM.rank(O)==D
 F=[[[seq[t][j][i]for t in range(L)]for j in range(m)]for i in range(n)]+[[[int(i==j)]for j in range(m)]for i in range(m)]
 G,u=PM.pm(F,L,[0]*n+[1]*m);ix=sorted(range(m+n),key=lambda i:(u[i],i))[:n];C=[[PM.coeff(G[i][j],u[i])for j in range(n)]for i in ix]
 # Independent explicit matrix powers for each Horner reconstruction.
 vecs=[]
 for i in ix:
  d=u[i];acc=[0]*D
  for a in range(d):
   q=[PM.coeff(G[i][j],a)for j in range(n)];z=mv(Y,q)
   for t in range(d-a-1):z=mv(B,z)
   acc=[(x+y)%PM.P for x,y in zip(acc,z)]
  assert not any((x+y)%PM.P for x,y in zip(mv(B,acc),mv(Y,[PM.coeff(G[i][j],d)for j in range(n)])))
  vecs.append(acc)
 recovered=[];W=nullbasis(transpose([row[10:]for row in C]));assert len(W)==11
 assert PM.rank(C)==n-(D-PM.rank(B))
 assert sum(u[i]for i in ix)==D
 for c in W:
  z=[sum(c[i]*vecs[i][h]for i in range(n))%PM.P for h in range(D)];s=[sum(c[i]*C[i][j]for i in range(n))%PM.P for j in range(10)];w=z[:N]+s
  assert not any(mv(A,w))
  if any(w):recovered.append(w)
 assert recovered
 j=next(j for j,x in enumerate(want)if x);v=recovered[0];scale=v[j]*pow(want[j],-1,PM.P)%PM.P;assert scale and all(x==scale*y%PM.P for x,y in zip(v,want))
 return dict(field=PM.P,alpha=alpha,rank_K=PM.rank(K),rank_O=PM.rank(O),p=pp,q=qq,L=L,nonzero_ideal=any(v[:N]),nonzero_s=any(v[N:]),rank_M=PM.rank(M),rank_A=PM.rank(A),kernel_B=D-PM.rank(B),dimension=D,blocks=[m,n],selected=len(recovered))

def run(seed=1,cases=4):
 old=PM.P;PM.P=65521;r=random.Random(seed);rows=[]
 try:
  for zero,ideal_zero in ((False,False),(True,False),(False,True)):
   for j in range(cases):rows.append(trial(r,zero,ideal_zero))
  return dict(ok=True,cases=len(rows),nonzero_s=sum(x['nonzero_s']for x in rows),zero_s=sum(not x['nonzero_s']for x in rows),ideal_zero=sum(not x['nonzero_ideal']for x in rows),details=rows)
 finally:PM.P=old
if __name__=='__main__':
 x=json.load(sys.stdin);print(json.dumps(run(x.get('seed',1),x.get('cases',4))))
