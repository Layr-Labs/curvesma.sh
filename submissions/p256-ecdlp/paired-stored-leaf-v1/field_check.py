import json,random,sys
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import a64
P=0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
R=pow(2,256,P);RI=pow(R,-1,P);MASK=(1<<64)-1
G=(0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296,
   0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5)
prog,labels=a64.a64_load((root/'p256_fields.s').read_text())
counts={k:set() for k in ['k_mul','k_sqr','k_add','k_sub']}
def call(k,a,b=0):
    mem={}
    for ptr,v in [(4096,a),(8192,b)]:
        for i in range(4):mem[ptr+8*i]=(v>>(64*i))&MASK
    n=a64.a64_run(prog,labels,k,dict(x0=12288,x1=4096,x2=8192,x30=0x123456),mem,limit=1000)
    counts[k].add(n)
    return sum(mem[12288+8*i]<<(64*i) for i in range(4))
cfg=json.load(sys.stdin);rng=random.Random(cfg['seed'])
values=[0,1,2,P-1,P-2,P//2,(1<<192)-1,1<<192,(1<<224)-1,1<<224]
samples=[(a,b) for a in values for b in values]
samples += [(rng.randrange(P),rng.randrange(P)) for _ in range(cfg['samples'])]
for a,b in samples:
    for k,w in [('k_mul',a*b*RI%P),('k_sqr',a*a*RI%P),('k_add',(a+b)%P),('k_sub',(a-b)%P)]:
        got=call(k,a,b)
        assert got==w,(k,a,b,got,w)
        assert 0<=got<P
assert counts=={'k_mul':{116},'k_sqr':{92},'k_add':{23},'k_sub':{19}},counts
def add(A,B):
    if A is None:return B
    if B is None:return A
    x,y=A;u,v=B
    if x==u:
        if (y+v)%P==0:return None
        slope=(3*x*x-3)*pow(2*y,-1,P)%P
    else:slope=(v-y)*pow(u-x,-1,P)%P
    xx=(slope*slope-x-u)%P
    return xx,(slope*(x-xx)-y)%P
A=G;B=add(G,G);pairs=0;skipped=0
for j in range(80):
    A=add(A,G);B=add(B,add(G,G))
    if A is None or B is None or A[0]==B[0]:skipped+=1;continue
    x,y=[v*R%P for v in A];u,v=[v*R%P for v in B]
    dx=call('k_sub',u,x)
    # Reference inverse only. No assembly-inversion cost is claimed here.
    inv=pow(dx*RI%P,-1,P)*R%P
    dp=call('k_sub',v,y);dm=call('k_add',v,y)
    lp=call('k_mul',dp,inv);lm=call('k_mul',dm,inv)
    sm=call('k_add',x,u)
    xp=call('k_sub',call('k_sqr',lp),sm)
    xm=call('k_sub',call('k_sqr',lm),sm)
    plus=add(A,B);minus=add(A,(B[0],(-B[1])%P))
    assert plus is not None and minus is not None
    assert xp*RI%P==plus[0] and xm*RI%P==minus[0]
    pairs+=1
prog2,labels2=a64.a64_load((root/'p256_fields.s').read_text()+'\n'+(root/'pair_batch.s').read_text()+'\n'+(root/'p256_inv.s').read_text())
batch_rows=[]
def put(mem,ptr,v):
    for i in range(4):mem[ptr+8*i]=(v>>(64*i))&MASK
def get(mem,ptr):return sum(mem[ptr+8*i]<<(64*i) for i in range(4))
for n in [1,2,3,4,8,16,32,64]:
    mem={};ys=[];X=add(G,G);Y=add(X,G)
    for j in range(n):
        Y=add(Y,G);ys.append(Y)
        put(mem,0x10000+64*j,Y[0]*R%P);put(mem,0x10020+64*j,Y[1]*R%P)
    put(mem,0x30000,X[0]*R%P);put(mem,0x30020,X[1]*R%P)
    put(mem,0x20000-32,R)
    regs={}
    fc=a64.a64_run(prog2,labels2,'pair_fwd',dict(x0=n,x1=0x10000,x2=0x20000,x3=0x30000,x4=0x40000,sp=0x50000,x30=12345),mem,limit=100000,out_regs=regs)
    assert regs['x0']==1
    total=get(mem,0x20000+64*(n-1)+32)
    ic=a64.a64_run(prog2,labels2,'k_inv',dict(x0=0x40000,x1=0x20000+64*(n-1)+32,x2=0x60000,sp=0x50000,x30=12345),mem,limit=200000)
    assert get(mem,0x40000)==pow(total*RI%P,-1,P)*R%P
    bc=a64.a64_run(prog2,labels2,'pair_bwd',dict(x0=n,x1=0x10000+64*(n-1),x2=0x20000+64*(n-1),x3=0x30000,x4=0x40000,sp=0x50000,x30=12345),mem,limit=200000)
    for j,Y in enumerate(ys):
        assert get(mem,0x20000+64*j)*RI%P==add(X,Y)[0]
        assert get(mem,0x20020+64*j)*RI%P==add(X,(Y[0],-Y[1]%P))[0]
    assert get(mem,0x40000)==R
    batch_rows.append({'n':n,'forward':fc,'backward':bc,'inverse':ic,'all_points_match':True})
# Deliberate zero denominator: the forward pass must abort before inversion.
mem={}
for ptr in [0x10000,0x30000]:
    put(mem,ptr,G[0]*R%P);put(mem,ptr+32,G[1]*R%P)
put(mem,0x20000-32,R);regs={}
zc=a64.a64_run(prog2,labels2,'pair_fwd',dict(x0=1,x1=0x10000,x2=0x20000,x3=0x30000,x4=0x40000,sp=0x50000,x30=12345),mem,limit=1000,out_regs=regs)
assert regs['x0']==0 and zc<=165
kp,kl=a64.a64_load((root/'pair_keys.s').read_text())
key_rows=[]
for n in [0,1,2,3,8,16]:
    for cap in [0,1,2,100]:
        mem={};expected=[];seen={};port=0x70000;dstbase=0x90000
        w2=17|(123<<36)
        for j in range(2*n):
            key=rng.randrange(P-(1<<36))
            if cap<100:key=(key & ~((1<<36)-1)) | (j%3)
            put(mem,0x10000+32*j,key)
            d1=key&((1<<36)-1);d2=(key>>64)&((1<<36)-1)
            cursor=port+64+8*d1
            mem[cursor]=dstbase
            idx=seen.get(d1,0);seen[d1]=idx+1
            if idx<cap:expected.append((d1,dstbase+32*idx,(key>>128)&MASK,key>>192,w2,d2|((j%2)<<36)|((j//2)<<37)))
        packet=[0]*6;got=[]
        def device(offset,value):
            packet[offset//8]=value
            if offset==40:got.append(tuple(packet))
        kc=a64.a64_run(kp,kl,'pair_keys',dict(x0=n,x1=0x10000,x2=port,x3=w2,x4=dstbase+32*cap,x30=12345),mem,limit=100000,dev=(port,48,device))
        assert got==expected,(n,cap)
        assert kc <= 29*n+11
        key_rows.append({'n':n,'cap':cap,'instructions':kc,'packets':len(got)})
import contextlib,io,runpy
buf=io.StringIO()
with contextlib.redirect_stdout(buf):runpy.run_path(str(root/'success_probe.py'),run_name='__main__')
probe=json.loads(buf.getvalue())
import routing_check,toy_recovery
toy=toy_recovery.check()
routing=routing_check.check()
assert all(routing[k] for k in ['messages_equal_emit_reference','slots_equal_reference','source_tags_ok','buckets_pointer_order_and_scan_ok','candidates_equal_exact_order_reference','formulas_ok'])
print(json.dumps({'probe':probe,'toy_recovery':toy,'routing':routing,'checks_ok':True,'field_operand_pairs':len(samples),
 'instruction_counts_including_RET':{k:sorted(v) for k,v in counts.items()},
 'zero_denominator_abort_verified':True,'paired_point_cases':pairs,'excluded_zero_denominators':skipped,'batch_rows':batch_rows,'key_rows':key_rows,
 'scope':'Field, complete batch inversion and key-emission prototype; no full attack, two-hop routing or success bound.'}))
