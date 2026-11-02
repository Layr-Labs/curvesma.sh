"""Finite comparison to frozen canonical loops, preserving raw labels."""
import random
import a64,range_witness
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
RR=1<<256;R=RR%P;RI=pow(R,-1,P);MASK=(1<<64)-1;D=RR-P;L=RR-(1<<224)
def check(root,build,descriptors):
 old,ol=a64.a64_load((root/'pair_inline_reference.s').read_text())
 new,nl=a64.a64_load((root/'p256_fields.s').read_text()+range_witness.generated_assembly(root))
 rng=random.Random(61202);cases=outputs=0;variants=set();noncanonical=0
 def put(m,p,v):
  for i in range(4):m[p+8*i]=(v>>(64*i))&MASK
 def get(m,p):return sum(m[p+8*i]<<(64*i) for i in range(4))
 for n in [1,2,7,16]:
  for trial in range(16):
   xx=[0,1,D,D+1,RR//2,P//2,P-1,rng.randrange(P)][trial%8];xy=rng.randrange(1,P//2+1)
   pool=[0,1,P-1,xx-1,xx+1,RR-xx-1,RR-xx,RR-xx+1]
   pool=[z for z in pool if 0<=z<P and z!=xx]
   ys=[(pool[j%len(pool)] if trial<7 else rng.randrange(P),1+rng.randrange(P//2),j) for j in range(n)]
   ys=[((x+1)%P if x==xx else x,y,j) for x,y,j in ys]
   if trial>=8:
    # Algebraic fixtures, not necessarily curve points. Force one of the
    # two affine formulas to a chosen Montgomery output near rare labels.
    ax=rng.randrange(P);by=rng.randrange(1,P);xx=ax*R%P;xy=by*R%P
    if xy>P//2:xy=P-xy
    ys=[]
    targets=[0,1,D-1,D,D+1,P-D-1,P-D+1,P-1]
    for j in range(n):
     target=targets[(trial-8+j)%len(targets)];m=rng.randrange(1,P)
     cx=(m*m-ax-target*RI)%P
     while cx==ax:
      m=(m+1)%P;cx=(m*m-ax-target*RI)%P
     dy=(by+m*(cx-ax))%P;ym=dy*R%P
     if ym>P//2:ym=P-ym
     ys.append((cx*R%P,ym,j))
   tr,seq,nodes=build(ys);ds=descriptors(tr,seq,xx,xy)
   xp=0x30000;bp=0x40000;scratch=0x60000;yp=0xA0000;sp=0x80000
   m={};put(m,xp,xx);put(m,xp+32,xy);put(m,bp-32,R)
   for nd in nodes:
    for j,(x,y,idx) in enumerate(nd['items']):
     put(m,nd['base']+72*j,x);put(m,nd['base']+72*j+32,y);m[nd['base']+72*j+64]=idx<<36
   order=[];segments=[];off=0
   def regs(**kw):
    z={f'x{i}':rng.getrandbits(64) for i in range(31)};z.update(sp=sp,x30=12345);z.update(kw);return z
   for nd,a,b,dx,yl,carry in ds:
    cnt=b-a;ptr=nd['base']+72*a;out={}
    a64.a64_run(new,nl,'pair_fwd_lo' if dx else 'pair_fwd_hi',regs(x0=cnt,x1=ptr,x2=bp+64*off,x3=xp,x4=scratch),m,limit=300000,out_regs=out)
    assert out['x0']==1 and out['sp']==sp and out['x30']==12345
    segments.append((off,cnt,ptr,yl,carry));order.extend(nd['items'][a:b]);off+=cnt
   mo={};put(mo,xp,xx);put(mo,xp+32,xy);put(mo,bp-32,R)
   for j,(x,y,idx) in enumerate(order):put(mo,yp+64*j,x);put(mo,yp+64*j+32,y)
   out={};a64.a64_run(old,ol,'pair_fwd',regs(x0=n,x1=yp,x2=bp,x3=xp,x4=scratch),mo,limit=300000,out_regs=out);assert out['x0']==1
   for mm in [m,mo]:put(mm,scratch,pow(get(mm,bp+64*(n-1)+32)*RI%P,-1,P)*R%P)
   for off,cnt,ptr,yl,carry in reversed(segments):
    variants.add((yl,carry));out={}
    a64.a64_run(new,nl,f'pair_bwd_y{yl}_c{carry}',regs(x0=cnt,x1=ptr+72*(cnt-1),x2=bp+64*(off+cnt-1),x3=xp,x4=scratch),m,limit=300000,out_regs=out)
    assert out['sp']==sp and out['x30']==12345
   a64.a64_run(old,ol,'pair_bwd',regs(x0=n,x1=yp+64*(n-1),x2=bp+64*(n-1),x3=xp,x4=scratch),mo,limit=300000)
   for j in range(2*n):
    a=get(m,bp+32*j);b=get(mo,bp+32*j);allowed={b}
    if b<D:allowed.add(b+P)
    if b>P-D:allowed.add(b+D)
    assert a in allowed,(trial,n,j,a,b)
    if a<L:assert a==b
    else:assert b<D or b>=P-D
    noncanonical+=a!=b;outputs+=1
   assert get(m,scratch)==R and get(mo,scratch)==R;cases+=1
 assert variants=={(0,0),(0,1),(1,0),(1,1)}
 assert noncanonical>0, 'planted fixtures must exercise noncanonical labels'
 return dict(cases=cases,raw_outputs=outputs,noncanonical_observed=noncanonical,all_four_variants=True,raw_label_membership=True,scope='finite boundary/random raw differential, not universal proof or success claim')
