import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
def check():
 affine=39479+15*116+512
 dbl_i=6*92+3*116+16*23;madd_i=4*92+7*116+12*23;add_i=5*92+11*116+14*23;jop=4096
 scalar=(15+256+64)*jop+65536+45056
 generator=scalar+45056+262144;yitem=scalar+262144
 product=64+256*(2*64+16)
 inversion=512*(product+32)+512
 resolver=inversion+scalar+2*45056+(1<<18)
 assert affine==41731<45056 and(dbl_i,madd_i,add_i)==(1268,1456,2058)and max(dbl_i,madd_i,add_i)<=jop
 assert(scalar,generator,yitem)==(1482752,1789952,1744896)and generator<(1<<21)
 assert product==36928 and inversion==18924032<(1<<25)
 assert resolver==20759040<(1<<25)
 PF=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff;BC=0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
 GP=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
 O=(1,1,0)
 def dbl(Q):
  X,Y,Z=Q
  if Z==0 or Y==0:return O
  A_=X*X%PF;B_=Y*Y%PF;C=B_*B_%PF;Z2=Z*Z%PF;D=2*((X+B_)**2-A_-C)%PF;E=3*(A_-Z2*Z2)%PF;X3=(E*E-2*D)%PF
  return(X3,(E*(D-X3)-8*C)%PF,2*Y*Z%PF)
 def add(Q,R_):
  if Q[2]==0:return R_
  if R_[2]==0:return Q
  X1,Y1,Z1=Q;X2,Y2,Z2=R_;Z11=Z1*Z1%PF;Z22=Z2*Z2%PF;U1=X1*Z22%PF;U2=X2*Z11%PF;S1=Y1*Z2*Z22%PF;S2=Y2*Z1*Z11%PF
  H=(U2-U1)%PF;r_=2*(S2-S1)%PF
  if H==0:return dbl(Q)if r_==0 else O
  I=4*H*H%PF;J=H*I%PF;V=U1*I%PF;X3=(r_*r_-J-2*V)%PF
  return(X3,(r_*(V-X3)-2*S1*J)%PF,((Z1+Z2)**2-Z11-Z22)*H%PF)
 def aff(Q):
  if Q[2]==0:return None
  i=pow(Q[2],-1,PF);return(Q[0]*i*i%PF,Q[1]*i*i*i%PF)
 def ref(k,P_):
  acc=None
  for b in bin(k)[2:]if k else'':
   if acc:
    l=3*(acc[0]**2-1)*pow(2*acc[1],-1,PF)%PF;x=(l*l-2*acc[0])%PF;acc=(x,(l*(acc[0]-x)-acc[1])%PF)
   if b=='1':
    if acc is None:acc=P_
    elif acc[0]==P_[0]:acc=None if(acc[1]+P_[1])%PF==0 else 0
    else:
     l=(P_[1]-acc[1])*pow(P_[0]-acc[0],-1,PF)%PF;x=(l*l-acc[0]-P_[0])%PF;acc=(x,(l*(acc[0]-x)-acc[1])%PF)
  return acc
 def window(k,P_):
  ops=0;T=[O,(P_[0],P_[1],1)]
  for d in range(2,16):T.append(add(T[-1],T[1]));ops+=1
  ops+=1;acc=O
  for i in range(63,-1,-1):
   for _ in range(4):acc=dbl(acc);ops+=1
   n=(k>>(4*i))&15
   if n:acc=add(acc,T[n]);ops+=1
  assert ops<=335;return aff(acc)
 rng=random.Random(20261101);v=[0,1,2,R-1,R-2,R//2,(1<<128)-1,1<<128,(1<<255)-1,1<<255]
 assert(GP[1]**2-GP[0]**3+3*GP[0]-BC)%PF==0 and ref(R,GP)is None and ref(R+1,GP)==GP and R<PF
 import json,pathlib;tp=[q/'tracks/p256-ecdlp/target-profile.json'for q in pathlib.Path(__file__).resolve().parents];tp=[q for q in tp if q.exists()]
 if tp:
  d=json.load(open(tp[0]));h=lambda k:int(d[k],16)
  assert(h('field_modulus'),h('subgroup_order'),h('coefficient_a'),h('coefficient_b'),int(d['generator_P']['x'],16),int(d['generator_P']['y'],16))==(PF,R,PF-3,BC,GP[0],GP[1])
 Q5=ref(5,GP);wins=0
 for P_ in(GP,Q5):
  for k in[0,1,2,3,15,16,17,255,256,R-1,R-2,R//2]+[rng.randrange(R)for _ in range(6)]:
   assert window(k,P_)==ref(k,P_);wins+=1
 pairs=[(a,b)for a in v for b in v]+[(rng.randrange(R),rng.randrange(R))for _ in range(128)]
 counts=[]
 for a,b in pairs:
  acc=0;count=64
  for i in range(255,-1,-1):
   acc+=acc
   if acc>=R:acc-=R
   count+=64+16
   if(b>>i)&1:
    acc+=a
    if acc>=R:acc-=R
    count+=64
   assert 0<=acc<R
  assert acc==(a*b)%R and count<=product;counts.append(count)
 PF=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
 def frac_div(v,s):
  rem=v;q=0;count=64
  for _ in range(s):
   rem+=rem;q+=q;count+=32
   if rem>=PF:rem-=PF;q+=1
  assert q==(v<<s)//PF and 0<=rem<PF;return count
 fcases=[0,1,PF-1,PF//2,(1<<255)%PF,(1<<192)%PF]+[rng.randrange(PF)for _ in range(64)]
 dmax=max(max(frac_div(v,256),frac_div(v,192))for v in fcases)
 per_y=5*8256+8*116+1024;per_x=9*8256+16*116+1024
 assert dmax==8256<=8256 and per_y<(1<<16) and per_x<(1<<17)
 return{'checks_ok':True,'mod_product_cases':len(pairs),'protected_constants':'intrinsic+profile'if tp else'intrinsic','jacobian_window_cases':wins,'fraction_division_cases':2*len(fcases),'fraction_division':dmax,'fraction_entry_y':per_y,'fraction_area_x':per_x,'bounds':{'affine':affine,'jacobian_doubling':dbl_i,'jacobian_mixed_add':madd_i,'jacobian_add':add_i,'per_op':jop,'scalar_mul':scalar,'gen_header':generator,'y_header':yitem,'scalar_mod_product':product,'scalar_inversion':inversion,'resolver':resolver},'scope':'Exact inequalities and finite modular checks for explicit analytic construction bounds, not executed assembly counts or production attack.'}
