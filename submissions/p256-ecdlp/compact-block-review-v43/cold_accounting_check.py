import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
def check():
 doubling=39479+2*116+2*92+(4*23+4*19)+(13*16+256)
 addition=39479+2*116+92+6*19+(10*16+256)
 # Projective fixed-window-4 scalar multiplication (Jacobian, a=-3):
 # doubling 2M+7S+16A, mixed addition 7M+4S+13A, 16/wrapper-call + 256 fixed.
 pdoubling=2*116+7*92+16*23+(9*16+256)
 paddition=7*116+4*92+13*23+(11*16+256)
 conversion=39479+92+4*116+16+256
 scalar=256*pdoubling+64*paddition+64*5+384+8192+conversion
 table=2*15*40333
 bookkeeping=64+2*116+512
 generator_x=scalar+addition+bookkeeping
 generator_y=scalar+bookkeeping
 product=64+256*(2*64+16)
 inversion=512*(product+32)+512
 resolver=inversion+scalar+2*doubling+(1<<18)
 assert doubling==40527 and addition==40333
 assert pdoubling==1644 and paddition==1911 and conversion==40307
 assert scalar==592371 and generator_x==633512 and generator_y==593179<(1<<26)
 assert table==1209990
 assert product==36928 and inversion==18924032<(1<<25)
 assert resolver==19859601<(1<<26)
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
 p256=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
 Gx=0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
 Gy=0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
 def padd(P,Q):
  if P is None:return Q
  if Q is None:return P
  x1,y1=P;x2,y2=Q
  if x1==x2:return None if (y1+y2)%p256==0 else pdoub(P)
  l=(y2-y1)*pow(x2-x1,-1,p256)%p256
  x3=(l*l-x1-x2)%p256;return (x3,(l*(x1-x3)-y1)%p256)
 def pdoub(P):
  if P is None:return None
  x1,y1=P;l=(3*x1*x1-3)*pow(2*y1,-1,p256)%p256
  x3=(l*l-2*x1)%p256;return (x3,(l*(x1-x3)-y1)%p256)
 def jdbl(X1,Y1,Z1):
  A=X1*X1%p256;B=Y1*Y1%p256;C=B*B%p256;D=2*((X1+B)**2-A-C)%p256
  Z1Z1=Z1*Z1%p256;Z1Q=Z1Z1*Z1Z1%p256;E=3*(A-Z1Q)%p256;F=E*E%p256
  return ((F-2*D)%p256,(E*(D-(F-2*D)%p256)-8*C)%p256,2*Y1*Z1%p256)
 def jmadd(X1,Y1,Z1,x2,y2):
  Z1Z1=Z1*Z1%p256;U2=x2*Z1Z1%p256;S2=y2*(Z1*Z1Z1%p256)%p256
  H=(U2-X1)%p256;HH=H*H%p256;I=4*HH%p256;J=H*I%p256
  r=2*(S2-Y1)%p256;V=X1*I%p256
  return ((r*r-J-2*V)%p256,(r*(V-(r*r-J-2*V)%p256)-2*Y1*J)%p256,((Z1+H)**2-Z1Z1-HH)%p256)
 def jaff(X,Y,Z):
  if Z%p256==0:return None
  zi=pow(Z,-1,p256);z2=zi*zi%p256;return (X*z2%p256,Y*z2*zi%p256)
 def jmul(k):
  if k==0:return None
  T=[None]+[(Gx,Gy)]
  for i in range(2,16):
   T.append(padd(T[-1],(Gx,Gy)))
  acc=None
  for d in [(k>>(4*i))&15 for i in range((k.bit_length()+3)//4-1,-1,-1)]:
   if acc is not None:
    X,Y,Z=acc;acc=jdbl(X,Y,Z)
    X,Y,Z=acc;acc=jdbl(X,Y,Z)
    X,Y,Z=acc;acc=jdbl(X,Y,Z)
    X,Y,Z=acc;acc=jdbl(X,Y,Z)
   if d:
    if acc is None:acc=(T[d][0],T[d][1],1)
    else:X,Y,Z=acc;acc=jmadd(X,Y,Z,T[d][0],T[d][1])
  return jaff(*acc)
 def amul(k):
  R=None;Q=(Gx,Gy)
  while k:
   if k&1:R=padd(R,Q)
   Q=pdoub(Q);k>>=1
  return R
 agree=0
 for k in [1,2,3,15,16,17,255,256,257]+[rng.randrange(1,1<<48) for _ in range(24)]:
  assert jmul(k)==amul(k);agree+=1
 return {'checks_ok':True,'modular_product_cases':len(pairs),'projective_agreement_cases':agree,'construction_bounds':{'doubling':doubling,'addition':addition,'projective_doubling':pdoubling,'projective_addition':paddition,'conversion':conversion,'scalar_multiplication':scalar,'window_table':table,'generator_x':generator_x,'generator_y':generator_y,'scalar_modular_product':product,'scalar_inversion':inversion,'resolver':resolver},'scope':'Exact inequalities and finite modular checks for explicit analytic construction bounds, not executed assembly counts or production attack.'}
