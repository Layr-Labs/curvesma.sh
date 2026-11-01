"""Compare retained results to original canonical inline loops; prefixes may differ in sign."""
import random
import a64
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
RR=1<<256;R=RR%P;RI=pow(R,-1,P);MASK=(1<<64)-1;L=RR-(1<<224);D=RR-P

def check(root):
 old,ol=a64.a64_load((root/'pair_inline_reference.s').read_text())
 new,nl=a64.a64_load((root/'pair_batch.s').read_text())
 rng=random.Random(59017);cases=0
 def put(m,p,v):
  for i in range(4):m[p+8*i]=(v>>(64*i))&MASK
 def get(m,p):return sum(m[p+8*i]<<(64*i) for i in range(4))
 for n in [0,1,2,3,7,15,16,31]:
  for trial in range(4):
   yp=0x10000+64*rng.randrange(128);bp=0x30000+64*rng.randrange(128);xp=0x50000+64*rng.randrange(128);sptr=0x70000+64*rng.randrange(128)
   m={};xx=rng.randrange(P);xy=rng.randrange(P);put(m,xp,xx);put(m,xp+32,xy);put(m,bp-32,R)
   for j in range(n):
    yx=rng.randrange(P)
    if yx==xx:yx=(yx+1)%P
    put(m,yp+64*j,yx);put(m,yp+64*j+32,rng.randrange(P))
   regs={f'x{i}':rng.getrandbits(64) for i in range(31)};regs.update(x0=n,x1=yp,x2=bp,x3=xp,x4=sptr,sp=0x90000,x30=0x12345)
   mo=dict(m);mn=dict(m)
   for pr,la,mm in [(old,ol,mo),(new,nl,mn)]:
    rr={};a64.a64_run(pr,la,'pair_fwd',regs,mm,limit=100000,out_regs=rr);assert rr['x0']==1
    total=get(mm,bp+64*(n-1)+32) if n else R
    put(mm,sptr,pow(total*RI%P,-1,P)*R%P)
   dead={sptr+32+8*j:rng.getrandbits(64) for j in range(4)};mo.update(dead);mn.update(dead)
   regs.update(x0=n,x1=yp+64*(n-1),x2=bp+64*(n-1),x3=xp,x4=sptr)
   for pr,la,mm in [(old,ol,mo),(new,nl,mn)]:
    rr={};a64.a64_run(pr,la,'pair_bwd',regs,mm,limit=100000,out_regs=rr);assert rr['sp']==0x90000 and rr['x30']==0x12345
   for j in range(2*n):
    a=get(mn,bp+32*j);b=get(mo,bp+32*j)
    if a<L:assert a==b
    else:assert b<D or b>=P-D
   for start,size in [(yp,64*n),(xp,64)]:assert all(mn.get(k)==m.get(k) for k in range(start,start+size,8))
   assert get(mn,sptr)==R and all(mn[k]==v for k,v in dead.items())
   cases+=1
 return dict(cases=cases,retained_outputs_match=True,filtered_true_keys_in_proved_set=True,random_registers_and_addresses=True,SP_LR_preserved=True,inputs_and_dead_scratch_preserved=True)
