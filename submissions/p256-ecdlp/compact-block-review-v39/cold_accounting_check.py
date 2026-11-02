import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
def check():
 affine=39479+15*116+512
 scalar=511*41731+8192
 generator_x=scalar+41731+16384
 generator_y=scalar+16384
 product=64+256*(2*64+16)
 inversion=512*(product+32)+512
 resolver=inversion+scalar+2*41731+(1<<18)
 assert affine==41731
 assert scalar==21332733 and generator_x==21390848 and generator_y==21349117<(1<<26)
 assert product==36928 and inversion==18924032<(1<<25)
 assert resolver==40602371<(1<<26)
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
 return {'checks_ok':True,'modular_product_cases':len(pairs),'construction_bounds':{'affine':affine,'scalar_multiplication':scalar,'generator_x':generator_x,'generator_y':generator_y,'scalar_modular_product':product,'scalar_inversion':inversion,'resolver':resolver},'scope':'Exact inequalities and finite modular checks for explicit analytic construction bounds, not executed assembly counts or production attack.'}
