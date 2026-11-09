import random
R=0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
def check():
 affine=39479+15*116+512
 scalar=512*45056+65536
 generator=scalar+45056+262144
 product=64+256*(2*64+16)
 inversion=512*(product+32)+512
 resolver=inversion+scalar+2*45056+(1<<18)
 assert affine==41731<45056
 assert scalar==23134208 and generator==23441408<(1<<26)
 assert product==36928 and inversion==18924032<(1<<25)
 assert resolver==42410496<(1<<26)
 rng=random.Random(20261101);v=[0,1,2,R-1,R-2,R//2,(1<<128)-1,1<<128,(1<<255)-1,1<<255]
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
 return{'checks_ok':True,'mod_product_cases':len(pairs),'fraction_division_cases':2*len(fcases),'fraction_division':dmax,'fraction_entry_y':per_y,'fraction_area_x':per_x,'bounds':{'affine':affine,'scalar_mul':scalar,'gen_header':generator,'scalar_mod_product':product,'scalar_inversion':inversion,'resolver':resolver},'scope':'Exact inequalities and finite modular checks for explicit analytic construction bounds, not executed assembly counts or production attack.'}
