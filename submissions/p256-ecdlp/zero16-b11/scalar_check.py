import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
MR=(1<<256)%P;MRI=pow(MR,-1,P)
G=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
W=(1<<64)-1
def limbs(x,n=4):return [(x>>(64*i))&W for i in range(n)]
def integer(a):return sum(v<<(64*i) for i,v in enumerate(a))
def add4(a,b):
 out=[];c=0
 for x,y in zip(a,b):z=x+y+c;out.append(z&W);c=z>>64
 return out+[c]
def reduce5(z):
 out=[];borrow=0
 for a,b in zip(z,limbs(R,5)):
  t=a-b-borrow;out.append(t&W);borrow=int(t<0)
 out=z if borrow else out
 assert out[4]==0 and integer(out)<R
 return out[:4]
def modmul(a,b):
 assert 0<=a<R and 0<=b<R
 aa=limbs(a);bb=limbs(b);x=[0]*4
 for i in range(255,-1,-1):
  x=reduce5(add4(x,x))
  if (bb[i//64]>>(i%64))&1:x=reduce5(add4(x,aa))
 return integer(x)
def halfcoef(x):
 assert 0<=x<R
 z=add4(limbs(x),limbs(R)) if x&1 else limbs(x,5)
 out=[(z[i]>>1)|((z[i+1]&1)<<63) for i in range(4)]
 y=integer(out);assert 0<=y<R and (2*y-x)%R==0
 return y
def inv(a,stats=None):
 assert 0<a<R
 u,v=a,R;x,y=1,0;halves=subs=0
 while u!=1 and v!=1:
  while u%2==0:u//=2;x=halfcoef(x);halves+=1
  while v%2==0:v//=2;y=halfcoef(y);halves+=1
  if u==1 or v==1:break
  if u>=v:u-=v;x=(x-y)%R
  else:v-=u;y=(y-x)%R
  subs+=1
  assert u and v and u%R==a*x%R and v%R==a*y%R
 out=x if u==1 else y
 assert halves<=512 and subs<=512 and modmul(a,out)==1
 if stats is not None:stats.append((halves,subs))
 return out
def add(A,B):
 if A is None:return B
 if B is None:return A
 x,y=A;u,v=B
 if x==u:
  if (y+v)%P==0:return None
  m=(3*x*x-3)*pow(2*y,-1,P)%P
 else:m=(v-y)*pow(u-x,-1,P)%P
 z=(m*m-x-u)%P;return z,(m*(x-z)-y)%P
def mul(k,A=G):
 z=None
 while k:
  if k&1:z=add(z,A)
  A=add(A,A);k>>=1
 return z
def neg(A):return None if A is None else (A[0],-A[1]%P)
def leaf(vals,typ):
 xx,xy,yx,yy,a,b,c=vals;sign=1 if typ==0 else -1
 X=(xx*MRI%P,xy*MRI%P);Y=(yx*MRI%P,yy*MRI%P)
 return add(X,Y if sign==1 else neg(Y)),(a+sign*c)%R,b
def recover(left,right,typ1,typ2):
 A,m,n=leaf(left,typ1);B,mm,nn=leaf(right,typ2)
 if A is None or B is None:return None
 for sign in [1,-1]:
  if (A if sign==1 else neg(A))==B:
   den=(nn-sign*n)%R
   if den:return modmul((sign*m-mm)%R,inv(den))
 return None
def check():
 rng=random.Random(931);samples=[0,1,2,R-2,R-1]+[rng.randrange(R) for _ in range(12)];stats=[]
 for a in samples:
  for b in [0,1,R-1,rng.randrange(R)]:assert modmul(a,b)==a*b%R
  assert halfcoef(a)==((a+R)//2 if a&1 else a//2)
  if a:assert inv(a,stats)==pow(a,-1,R)
 assert R-2+(R)>2**256
 return dict(modular_products=len(samples)*4,inversions=len(stats),max_halvings=max(h for h,s in stats),max_subtractions=max(s for h,s in stats),fifth_word_halving_checked=True,all_checks=True,scope='Word-level modular product and carry-aware coefficient halving plus binary-GCD invariants; conservative instruction envelope is proved separately.')
