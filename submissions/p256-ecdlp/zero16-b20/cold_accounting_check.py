import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
G=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
C=[0,0,0]
def fm(a,b):C[0]+=1;return a*b%P
def fs(a):C[1]+=1;return a*a%P
def fa(a,b):C[2]+=1;return (a+b)%P
def fd(a,b):C[2]+=1;return (a-b)%P
def jdbl(X1,Y1,Z1):
 delta=fs(Z1);gamma=fs(Y1);beta=fm(X1,gamma);t2=fm(fd(X1,delta),fa(X1,delta));alpha=fa(fa(t2,t2),t2)
 b2=fa(beta,beta);b4=fa(b2,b2);b8=fa(b4,b4);X3=fd(fs(alpha),b8)
 yz=fa(Y1,Z1);Z3=fd(fd(fs(yz),gamma),delta)
 g2=fs(gamma);g4=fa(g2,g2);g8=fa(g4,g4);g16=fa(g8,g8);Y3=fd(fm(alpha,fd(b4,X3)),g16)
 return X3,Y3,Z3
def jmadd(X1,Y1,Z1,X2,Y2):
 Z1Z1=fs(Z1);U2=fm(X2,Z1Z1);S2=fm(Y2,fm(Z1,Z1Z1));H=fd(U2,X1);assert H
 HH=fs(H);I2=fa(HH,HH);I=fa(I2,I2);J=fm(H,I);rr=fd(S2,Y1);r=fa(rr,rr);V=fm(X1,I)
 V2=fa(V,V);X3=fd(fd(fs(r),J),V2);zh=fa(Z1,H);Z3=fd(fd(fs(zh),Z1Z1),HH)
 YJ=fm(Y1,J);YJ2=fa(YJ,YJ);Y3=fd(fm(r,fd(V,X3)),YJ2)
 return X3,Y3,Z3
def ladder(k,A):
 acc=None;d=m=0
 for i in range(255,-1,-1):
  if acc is not None:acc=jdbl(*acc);d+=1
  if (k>>i)&1:
   if acc is None:acc=(A[0],A[1],1)
   else:acc=jmadd(*acc,*A);m+=1
 if acc is None:return None,d,m
 X,Y,Z=acc;zi=pow(Z,-1,P);zi2=fs(zi);zi3=fm(zi2,zi);return (fm(X,zi2),fm(Y,zi3)),d,m
def aff_add(A,B):
 if A is None:return B
 if B is None:return A
 x,y=A;u,v=B
 if x==u:
  if (y+v)%P==0:return None
  m=(3*x*x-3)*pow(2*y,-1,P)%P
 else:m=(v-y)*pow(u-x,-1,P)%P
 z=(m*m-x-u)%P;return z,(m*(x-z)-y)%P
def aff_mul(k,A):
 z=None
 while k:
  if k&1:z=aff_add(z,A)
  A=aff_add(A,A);k>>=1
 return z
def check():
 affine=39479+15*116+512
 jd=3*116+5*92+16*23;jm=7*116+4*92+14*23
 scalar=256*(jd+jm+2*1024)+(39479+92+3*116+512)+65536
 generator=scalar+45056+262144
 product=64+256*(2*64+16)
 inversion=512*(product+32)+512
 resolver=inversion+scalar+2*45056+(1<<18)
 assert affine==41731<45056
 assert 24*16+256<=1024 and 25*16+256<=1024 and 5*16+256<=512
 assert jd==1176 and jm==1502 and scalar==1315823 and generator==1623023<(1<<21)
 assert product==36928 and inversion==18924032<(1<<25)
 assert resolver==20592111<(1<<25)
 rng=random.Random(20261101);v=[0,1,2,R-1,R-2,R//2,(1<<128)-1,1<<128,(1<<255)-1,1<<255]
 pairs=[(a,b) for a in v for b in v]+[(rng.randrange(R),rng.randrange(R)) for _ in range(128)]
 counts=[]
 for a,b in pairs:
  acc=0;count=64
  for i in range(255,-1,-1):
   acc+=acc
   if acc>=R:acc-=R
   count+=64+16
   if (b>>i)&1:
    acc+=a
    if acc>=R:acc-=R
    count+=64
   assert 0<=acc<R
  assert acc==(a*b)%R and count<=product;counts.append(count)
 jc=0;bases=[G,aff_mul(7,G),aff_mul(rng.randrange(1,R),G),aff_mul(rng.randrange(1,R),G)]
 for A in bases:
  for k in [0,1,2,3,R-1,R-2,(R-1)//2,(R+1)//2,1<<255,(1<<255)+1]+[rng.randrange(R) for _ in range(24)]:
   C[:]=[0,0,0];got,d,m=ladder(k,A);ref=aff_mul(k,A)
   assert got==ref and d<=255 and m<=255 and C[0]==3*d+7*m+(0 if got is None else 3) and C[1]==5*d+4*m+(0 if got is None else 1) and C[2]==16*d+14*m,(k,C,d,m)
   jc+=1
 C[:]=[0,0,0];jdbl(*(G[0],G[1],1));assert C==[3,5,16];C[:]=[0,0,0];jmadd(G[0],G[1],1,*bases[1]);assert C==[7,4,14]
 return {'checks_ok':True,'mod_product_cases':len(pairs),'jacobian_ladders':jc,'jacobian_ops':{'dbl':[3,5,16],'madd':[7,4,14]},'bounds':{'affine':affine,'jdbl':jd,'jmadd':jm,'scalar_mul':scalar,'gen_header':generator,'scalar_mod_product':product,'scalar_inversion':inversion,'resolver':resolver},'scope':'Exact inequalities, finite modular checks and exact Jacobian-ladder op counts for explicit analytic construction bounds, not executed assembly counts or production attack.'}
