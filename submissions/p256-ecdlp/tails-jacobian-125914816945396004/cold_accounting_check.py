import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
A=P-3
G=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
   0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
def inv(v):return pow(v,P-2,P)
def dbl(Q):
 X,Y,Z=Q
 if Z==0 or Y==0:return(0,1,0)
 XX=X*X%P;YY=Y*Y%P;C=YY*YY%P;ZZ=Z*Z%P
 D=2*(((X+YY)**2-XX-C)%P)%P;E=3*((XX-ZZ*ZZ)%P)%P;F=E*E%P;X3=(F-2*D)%P
 return (X3,(E*((D-X3)%P)-8*C)%P,2*Y*Z%P)
def madd(Q,T):
 X1,Y1,Z1=Q;X2,Y2=T
 if Z1==0:return(X2,Y2,1)
 Z2=Z1*Z1%P;U=Y2*Z1%P*Z2%P;H=(X2*Z2-X1)%P;I=4*H*H%P;J=H*I%P
 r=2*((U-Y1)%P)%P;V=X1*I%P
 if H==0:return dbl(Q) if r==0 else(0,1,0)
 X3=(r*r-J-2*V)%P
 return (X3,(r*((V-X3)%P)-2*Y1*J)%P,((Z1+H)**2-Z2-H*H)%P)
def aff(Q):
 X,Y,Z=Q
 if Z==0:return None
 i=inv(Z);i2=i*i%P
 return (X*i2%P,Y*i2*i%P)
def win(k,T):
 tab=[None,T]
 for d in range(2,16):tab.append(aff(madd((tab[-1][0],tab[-1][1],1),T)))
 Q=(0,1,0)
 for i in range(63,-1,-1):
  for _ in range(4):Q=dbl(Q)
  d=(k>>(4*i))&15
  if d:Q=madd(Q,tab[d])
 return aff(Q)
def aa(X,Y):
 if X is None:return Y
 if Y is None:return X
 if X[0]==Y[0]:
  if (X[1]+Y[1])%P==0:return None
  l=(3*X[0]*X[0]+A)*inv(2*X[1])%P
 else:l=(Y[1]-X[1])*inv(Y[0]-X[0])%P
 x=(l*l-X[0]-Y[0])%P;return(x,(l*(X[0]-x)-X[1])%P)
def mul(k,T):
 Q=None;S=T
 while k:
  if k&1:Q=aa(Q,S)
  S=aa(S,S);k>>=1
 return Q
def check():
 affine=39479+15*116+512;dbl_items=6*92+3*116+16*23;madd_items=4*92+7*116+12*23
 jop=4096;scalar=(15+256+64)*jop+65536+45056
 generator=scalar+45056+262144;yitem=scalar+262144
 product=64+256*(2*64+16);inversion=512*(product+32)+512
 resolver=inversion+scalar+2*45056+(1<<18)
 assert affine==41731<45056 and dbl_items==1268<=jop and madd_items==1456<=jop
 assert (scalar,generator,yitem)==(1482752,1789952,1744896) and generator<(1<<21)
 assert product==36928 and inversion==18924032<(1<<25) and resolver==20759040<(1<<25)
 rg=random.Random(20261105)
 for k in (1,15,R-1,rg.randrange(1,R)):
  Q=win(k,G);assert Q==mul(k,G)
 return {'checks_ok':True,'bounds':{'affine':affine,'doubling_items':dbl_items,'mixed_add_items':madd_items,'per_op_bound':jop,'scalar_mul':scalar,'gen_header':generator,'y_header':yitem,'scalar_mod_product':product,'scalar_inversion':inversion,'resolver':resolver},'scope':'Exact inequalities and finite Jacobian point checks for the cold-generation construction; not executed assembly counts.'}
