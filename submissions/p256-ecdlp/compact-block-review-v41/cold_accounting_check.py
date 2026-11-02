import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
def check():
 doubling=39479+2*116+2*92+(4*23+4*19)+(13*16+256)
 addition=39479+2*116+92+6*19+(10*16+256)
 scalar=255*doubling+256*addition+8192
 bookkeeping=64+2*116+512
 generator_x=scalar+addition+bookkeeping
 generator_y=scalar+bookkeeping
 product=64+256*(2*64+16)
 inversion=512*(product+32)+512
 resolver=inversion+scalar+2*doubling+(1<<18)
 assert doubling==40527 and addition==40333
 assert scalar==20667825 and generator_x==20708966 and generator_y==20668633<(1<<26)
 assert product==36928 and inversion==18924032<(1<<25)
 assert resolver==39935055<(1<<26)
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
 return {'checks_ok':True,'modular_product_cases':len(pairs),'construction_bounds':{'doubling':doubling,'addition':addition,'scalar_multiplication':scalar,'generator_x':generator_x,'generator_y':generator_y,'scalar_modular_product':product,'scalar_inversion':inversion,'resolver':resolver},'scope':'Exact inequalities and finite modular checks for explicit analytic construction bounds, not executed assembly counts or production attack.'}
