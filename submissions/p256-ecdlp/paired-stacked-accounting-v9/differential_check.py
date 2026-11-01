"""Randomized-memory differential checks against the published inline loops."""
import random
import a64
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
R=pow(2,256,P);RI=pow(R,-1,P);MASK=(1<<64)-1

def check(root):
    old,ol=a64.a64_load((root/'pair_inline_reference.s').read_text())
    new,nl=a64.a64_load((root/'pair_batch.s').read_text())
    rng=random.Random(59017);cases=0
    def put(m,p,v):
        for i in range(4):m[p+8*i]=(v>>(64*i))&MASK
    def get(m,p):return sum(m[p+8*i]<<(64*i) for i in range(4))
    for n in [0,1,2,3,7,15,16,31]:
      for trial in range(4):
        yp=0x10000+64*rng.randrange(128);bp=0x30000+64*rng.randrange(128)
        xp=0x50000+64*rng.randrange(128);sptr=0x70000+64*rng.randrange(128)
        m={};xx=rng.randrange(P);xy=rng.randrange(P);put(m,xp,xx);put(m,xp+32,xy);put(m,bp-32,R)
        for j in range(n):
            yx=rng.randrange(P)
            if yx==xx:yx=(yx+1)%P
            put(m,yp+64*j,yx);put(m,yp+64*j+32,rng.randrange(P))
        regs={f'x{i}':rng.getrandbits(64) for i in range(31)}
        regs.update(x0=n,x1=yp,x2=bp,x3=xp,x4=sptr,sp=0x90000,x30=0x12345)
        mo=dict(m);mn=dict(m);ro={};rn={}
        co=a64.a64_run(old,ol,'pair_fwd',regs,mo,limit=100000,out_regs=ro)
        cn=a64.a64_run(new,nl,'pair_fwd',regs,mn,limit=100000,out_regs=rn)
        assert co-cn==5*n and ro['x0']==rn['x0']==1
        assert all(mo.get(p)==mn.get(p) for p in range(bp-32,bp+64*n,8))
        total=get(mo,bp+64*(n-1)+32) if n else R
        inv=pow(total*RI%P,-1,P)*R%P
        put(mo,sptr,inv);put(mn,sptr,inv)
        dead_words={sptr+32+8*j:rng.getrandbits(64) for j in range(4)}
        mo.update(dead_words);mn.update(dead_words)
        regs.update(x0=n,x1=yp+64*(n-1),x2=bp+64*(n-1),x3=xp,x4=sptr)
        co=a64.a64_run(old,ol,'pair_bwd',regs,mo,limit=100000)
        cn=a64.a64_run(new,nl,'pair_bwd',regs,mn,limit=100000)
        assert co-cn==34*n
        for start,size in [(bp,64*n),(yp,64*n),(xp,64),(sptr,32)]:
            assert all(mo.get(p)==mn.get(p) for p in range(start,start+size,8))
        assert get(mn,sptr)==R
        assert all(mn[k]==v for k,v in dead_words.items())
        cases+=1
    return {'cases':cases,'random_initial_registers':True,'random_disjoint_addresses':True,'retained_outputs_match':True,'eliminated_scratch_stores_leave_canaries_unchanged':True,'forward_saving_per_item':5,'backward_saving_per_item':34}
