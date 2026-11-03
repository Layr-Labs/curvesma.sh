"""Full source/new AArch64 word interpretation; no new native Callgrind claim."""
import sys,json,re,random,time
from pathlib import Path
root=Path(__file__).resolve().parent
if (root/'files').is_dir():root=root/'files'
sys.path.insert(0,str(root))
import a64,lazy_witness,kern2_witness
params=json.loads((root/'params.json').read_text());p=int(params['p'],16);R=pow(2,384,p);RI=pow(R,-1,p);MASK=(1<<64)-1
assert (1<<384)>8*p
BASE=0x80000

def prepare(text):
    code,data=text.split('.section .rodata',1)
    vals=[]
    for line in data.splitlines():
        if '.quad' in line:vals += [int(x.strip(),0) for x in line.split('.quad',1)[1].split(',')]
    assert sum(vals[i]<<(64*i) for i in range(6))==p
    assert sum(vals[i+6]<<(64*i) for i in range(6))==2*p
    assert vals[12]==(-pow(p,-1,1<<64))&MASK
    code=code.replace('__SIZEOF_POINTER__','8')
    code=re.sub(r'(?<![A-Za-z_])(-?\d+)\*(\d+)',lambda m:str(int(m[1])*int(m[2])),code)
    code,n=re.subn(r'\badrp\s+x4,\.Lsl_consts','mov x4,#'+str(BASE),code);assert n==1
    code,n=re.subn(r'\badd\s+x4,x4,:lo12:\.Lsl_consts','add x4,x4,#0',code);assert n==1
    code=re.sub(r'\b(lsl|lsr|asr)#',r'\1 #',code)
    code='\n'.join(l for l in code.splitlines() if not l.lstrip().startswith('#'))
    pr,labels=a64.a64_load(code)
    return pr,labels,{BASE+8*i:v for i,v in enumerate(vals)}

def put(m,a,v):
    for i in range(6):m[a+8*i]=(v>>(64*i))&MASK

def get(m,a):return sum(m[a+8*i]<<(64*i) for i in range(6))

def add(A,B):
    if A is None:return B
    if B is None:return A
    x,y=A;u,v=B
    if x==u:
        if (y+v)%p==0:return None
        t=3*x*x*pow(2*y,-1,p)%p
    else:t=(v-y)*pow(u-x,-1,p)%p
    xx=(t*t-x-u)%p
    return xx,(t*(x-xx)-y)%p

def run(cfg):
    start=time.time();meta=kern2_witness.check(root);rng=random.Random(cfg['seed'])
    old,ol,cm=prepare((root/'measure/side_lanes.S').read_text());new,nl,cm2=prepare((root/'side_k2.S').read_text());assert cm==cm2
    vals=[0,1,p-1,p-2,R];cases=[]
    for i,a in enumerate(vals):
      for j,b in enumerate(vals):
       for iv in vals:cases.append((vals[(i+j)%5],a,vals[(2*i+j)%5],b,iv))
    cases += [tuple(rng.randrange(p) for _ in range(5)) for _ in range(cfg['samples'])]
    G=(int(params['gx'],16),int(params['gy'],16));assert G[1]**2%p==(G[0]**3+4)%p
    A=G;B=add(G,G);curve_expected=[]
    for _ in range(32):
        A=add(A,G);B=add(B,add(G,G))
        if A[0]==B[0]:continue
        out=add(A,(B[0],-B[1]%p));assert out is not None
        cases.append((A[0]*R%p,A[1]*R%p,B[0]*R%p,B[1]*R%p,pow((B[0]-A[0])%p,-1,p)*R%p))
        curve_expected.append((len(cases)-1,(out[0]*R%p,out[1]*R%p)))
    cases.append((0,1,1,p-1,R))
    rows=[];fixed=set();evens=0;outputs=[];loose_seen=0;negative_total=0
    for begin in range(0,len(cases),8):
        group=cases[begin:begin+8];m=dict(cm);n=len(group);want=[];negatives=0
        for i,(x1,y1,x2,y2,iv) in enumerate(group):
            out=0x20000+256*i;pp=0x30000+256*i;tt=0x40000+256*i;ii=0x50000+256*i
            for j,v in enumerate([out,pp,tt,ii]):m[0x10000+32*i+8*j]=v
            put(m,pp,x1);put(m,pp+48,y1);put(m,tt,x2);put(m,tt+48,y2);put(m,ii,iv)
            mu=(y1+y2)*iv*RI%p;xx=(mu*mu*RI-x1-x2)%p;yy=(mu*((xx-x1)%p)*RI-y1)%p
            key=yy if yy%2==0 else p-yy
            radix=1<<384
            def raw_mont(a,b):return (a*b+((-a*b*pow(p,-1,radix))%radix)*p)//radix
            raw_mu=raw_mont(y1+y2,iv);raw_y=raw_mont(xx+p-x1,raw_mu)
            assert raw_mu<2*p and raw_y<2*p
            negatives+=raw_y<y1
            want.append((xx,yy,[(key>>(64*j))&MASK for j in [1,2,3]]))
        regs={f'x{i}':rng.getrandbits(64) for i in range(31)};regs.update(x0=0x10000,x1=n,x30=0x12345,sp=0x7000000)
        mo=dict(m);mn=dict(m);ro={};rn={}
        co=a64.a64_run(old,ol,'side_lanes',regs,mo,limit=20000,out_regs=ro)
        cn=a64.a64_run(new,nl,'side_lanes',regs,mn,limit=20000,out_regs=rn)
        saving=co-cn-5*negatives+2  # +2: once-per-call prologue pair (call-row funded)
        assert 55*n<=saving<=61*n,(co,cn,n,negatives)
        assert (saving-55*n)%6==0,(co,cn,n,negatives)
        negative_total+=negatives
        for addr,value in m.items():assert mo[addr]==mn[addr]==value
        mu_old=get(mo,regs['sp']-288);mu_new=get(mn,regs['sp']-288)
        assert 0<=mu_old<p and 0<=mu_new<2*p and mu_old==mu_new%p
        loose_seen+=mu_new>=p
        e=sum(yy%2==0 for _,yy,_ in want);evens+=e
        fixed.add(co-(1051*n-4*e));assert cn<=996*n+512
        for i,(xx,yy,key) in enumerate(want):
            o=0x20000+256*i
            assert get(mo,o)==get(mn,o)==xx
            assert get(mo,o+48)==get(mn,o+48)==yy
            assert [mo[o+96+8*j] for j in range(3)]==[mn[o+96+8*j] for j in range(3)]==key
            outputs.append((xx,yy))
        for reg in ['sp']+[f'x{i}' for i in range(18,30)]:assert ro[reg]==rn[reg]==regs[reg],reg
        rows.append({'n':n,'old':co,'new':cn,'saved':co-cn})
    assert len(fixed)==1 and next(iter(fixed))<=512,fixed
    assert loose_seen>0
    assert all(outputs[i]==v for i,v in curve_expected)
    S,C,ST,nn=128,0xA5A5A5A55A5A5A5A,0x60000,8
    for stride,slots,pad in ((128,25,16),(96,9,12)):
        m=dict(cm2)
        for i in range(slots):
            for w in range(pad):m[ST+stride*i+8*w]=C
        for i in range(nn):
            x1,y1,x2,y2,iv=cases[i]
            o=ST+stride*i;pp=ST+S*(16+i) if stride==128 else 0x30000+256*i;tt=0x40000+256*i;ii=0x50000+256*i
            for j,a in enumerate((o,pp,tt,ii)):m[0x10000+32*i+8*j]=a
            put(m,pp,x1);put(m,pp+48,y1);put(m,tt,x2);put(m,tt+48,y2);put(m,ii,iv)
        rg={f'x{i}':7 for i in range(31)};rg.update(x0=0x10000,x1=nn,x30=9,sp=0x7000000)
        pre=set(m)
        a64.a64_run(new,nl,'side_lanes',rg,m,limit=20000,out_regs={})
        assert all(0x6FFFE80<=a<0x7000000 for a in m if a not in pre)
        if stride==128:
            for i in range(nn):
                x1,y1,x2,y2,iv=cases[i];o,q=ST+S*i,ST+S*(16+i)
                mu=(y1+y2)*iv*RI%p;xx=(mu*mu*RI-x1-x2)%p;yy=(mu*((xx-x1)%p)*RI-y1)%p;ky=yy if yy%2==0 else p-yy
                assert get(m,o)==xx and get(m,o+48)==yy and [m[o+96+8*j] for j in range(3)]==[(ky>>(64*j))&MASK for j in [1,2,3]]
                assert m[o+120]==C and m[q+120]==C and get(m,q)==x1 and get(m,q+48)==y1
            assert all(m[ST+S*24+8*w]==C for w in range(16))
        else:
            assert any(m[ST+96*8+8*w]!=C for w in range(3))
    return {'status':'passed','witness':meta,'lanes':len(cases),'genuine_curve_cases':len(curve_expected),'lazy_representatives_observed':loose_seen,'input_memory_unchanged':True,'negative_y_difference_cases':negative_total,'even_keys':evens,'odd_keys':len(cases)-evens,'fixed_function_instructions':list(fixed),'new_worst_per_lane':996,'production_charge_per_lane_with_fixed_overhead':997,'rows':rows,'seconds':time.time()-start,'layout_audit':{'stride':128,'lanes':nn,'pad_trailing_inputs_intact':True,'footprint_exact':True,'neg96_control_detected':True},'scope':'Word interpreter + integer reference on side_k2.S; ADRP/lo12 relocated. Fixed calls in ledger row.'}
if __name__=='__main__':print(json.dumps(run(json.load(sys.stdin)),separators=(',',':')))
