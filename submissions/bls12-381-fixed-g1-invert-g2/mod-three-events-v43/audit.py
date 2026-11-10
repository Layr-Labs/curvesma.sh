"""Only shared llm-verifier credential-free Docker may execute this candidate.
Scoped scalar-instruction interpretation; virtual table values are a test oracle,
not a claim that production tables can be replaced by arithmetic for free.
"""
import ast,json,random,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent/'files'
P=int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab',16)
MASK=2**64-1
R=pow(2,384,P);RI=pow(R,-1,P)

def integer(s):
    node=ast.parse(s.strip().removeprefix('#'),mode='eval').body
    def go(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return n.value
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -go(n.operand)
        if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub,ast.Mult)):
            a,b=go(n.left),go(n.right)
            return a+b if isinstance(n.op,ast.Add) else a-b if isinstance(n.op,ast.Sub) else a*b
        raise AssertionError('noninteger immediate')
    return go(node)

def check_exact_ledger():
    # Reconstruct full-size asymmetric event/price formulas,not group-theorem evidence.
    from fractions import Fraction as F
    text=(ROOT/'proof.md').read_text()
    m,n,L=map(int,re.search(r'm=(\d+), n=(\d+), L=m\+n=(\d+)\.',text).groups())
    assert m==34021169017091701373279488505416777728
    assert n==52656969102635639157598878336710344704 and L==m+n
    assert m%(3*2**96)==n%(3*2**96)==0 and m<n<2**126
    u=0xd201000000010000;r=u**4-u*u+1;q=2**61-1
    ss=12*(m-m//u);M=n*ss;G0=5*n*n//3-n
    J=(n-1)//2;W=(n+u-1)//u+1
    windows=[(W,u-W)]
    assert 2*W<=u and u*u%3==1 and 2*n*n<r and 3*n+u<u*u
    assert 3*n+u*u*n+3*u**3*J<u*r
    Q,T=divmod(J,u);K=S=0
    for lo,hi in windows:
        assert W<=lo<=hi<=u-W
        width=hi-lo+1;last=max(0,min(hi,T)-lo+1)
        K+=Q*width+last
        S+=width*u*Q*(Q-1)//2+Q*(lo+hi)*width//2+last*Q*u+(lo+min(hi,T))*last//2
    Eeven=2*(n*K-2*S)
    Wodd=(2*n+u-1)//u+2
    assert 2*Wodd<=u and 8*n*n<r and 6*n+u<2*u*u
    assert 6*n+2*u*u*n+3*u**3*(n-1)<u*r
    rawwindows=[(0,(u-Wodd)//3),((u+Wodd+2)//3,u-Wodd),(u+Wodd,(5*u-Wodd)//3),((5*u+Wodd+2)//3,2*u-1)]
    Qo,To=divmod(n-1,2*u);Ko=So=0;oddwindows=[]
    for lo,hi in rawwindows:
        lo+=1-lo%2;hi-=1-hi%2;assert 0<=lo<=hi<2*u
        assert hi<=u-Wodd or lo>=u+Wodd
        # Entire3x interval lies in one component excluding odd-multiple neighborhoods.
        assert any(a<=3*lo<=3*hi<=b for a,b in [(0,u-Wodd),(u+Wodd,3*u-Wodd),(3*u+Wodd,5*u-Wodd),(5*u+Wodd,6*u)])
        oddwindows.append([lo,hi]);width=(hi-lo)//2+1;last=max(0,(min(hi,To)-lo)//2+1)
        Ko+=Qo*width+last
        So+=width*2*u*Qo*(Qo-1)//2+Qo*(lo+hi)*width//2+last*Qo*2*u+last*(lo+last-1)
    assert str(oddwindows) in text
    Eodd=2*(n*Ko-So);E=Eeven+Eodd;G=G0+E
    for name,value in [('Eeven',Eeven),('Eodd',Eodd),('Wodd',Wodd),('Ko',Ko),('So',So)]:
        assert int(re.search(r'^'+name+r'=(\d+)\.$',text,re.M)[1])==value,('odd-gap count',name)
    for name,value in [('W',W),('K',K),('S',S),('G0',G0),('E',E)]:
        assert int(re.search(r'^'+name+r'=(\d+)\.$',text,re.M)[1])==value,('order-four count',name)
    for name,value in [('s',ss),('M',M),('G',G)]:
        assert int(re.search(r'^'+name+r'=(\d+)[,.]$',text,re.M)[1])==value,('printed event count',name)
    assert m<u*u-1 and 2*m*m<r and m*u*u<r-m and u*m+m<r
    assert 3*n<u*u and 4*n*(u*u+1)<r
    mu=F(M,r)+F(G*(r-1),r*r);U2=F(M*(M-1)//2+M*G,r*r)
    delta=(F(1,2**63)+F(1,q))**3
    loss=(L*(L-1)//2*(F(12,r)+delta)+m//u)/(L//2**29+1)
    success=mu-U2-loss-F(n,r)-F(1,2**970)
    assert success>F(3900031,10000000)>F(39,100)
    sn,sd=map(int,re.search(r'Proposed time_log2=(\d+)/(\d+)=',text).groups())
    assert F(sn,sd)==F(1243913,10000)
    formula='N=6414L+12m+28n+2^124'
    assert formula in text
    N=6414*L+12*m+28*n+2**124
    assert N==557857494710942619131912428293860739776512
    assert N**sd<20000**sd*2**sn
    for doc in ['architecture.md','radix80_certificate.md','fp2_certificate.md']:
        stated=(ROOT/doc).read_text();assert formula in stated
        assert f'N^{sd}<20000^{sd}*2^{sn}' in stated,('current ledger transcription',doc)
    transport=(ROOT/'metadata40_certificate.md').read_text()
    got=tuple(map(int,re.search(r'payload4 costs(\d+),metadata6 costs(\d+),total(\d+)/512',transport).groups()))
    payload=512*4+264+64*3;metadata=512*6+264+64*5
    assert got==(payload,metadata,payload+metadata)==(2504,3656,6160)
    return {'printed_constants_reconstructed':True,'asymmetric_counts':True,'order_four_counts':True,'odd_gap_counts':True,'odd_quotient_windows':True,'even_mod_three_separation':True,'order_four_separation_inequalities':True,
            'exact_success_above_03900031':True,'weighted_scalar_cost':True,
            'transport_formula_checked':True,'score_numerator':sn,'score_denominator':sd,
            'exact_score_inequality':True,
            'scope':'exact full-size integer/count/price transcription only;not group-theorem or physical attack certification'}


def parse(name):
    rows=[s.split('//')[0].strip() for s in (ROOT/name).read_text().splitlines()]
    macros={};ordinary=[];i=0
    while i<len(rows):
        s=rows[i];i+=1
        if s.startswith('.macro '):
            head=s.split(None,2);formals=[] if len(head)==2 else head[2].split(',');body=[]
            while rows[i]!='.endm':body.append(rows[i]);i+=1
            i+=1;macros[head[1]]=([x.strip() for x in formals],body)
        else:ordinary.append(s)
    out=[]
    def expand(rows,depth=0):
        assert depth<=4
        expanded=[];at=0
        while at<len(rows):
            if rows[at].startswith('.rept '):
                times=integer(rows[at][6:]);assert 0<=times<=512
                end=at+1;body=[]
                while end<len(rows) and rows[end]!='.endr':
                    assert not rows[end].startswith('.rept ')
                    body.append(rows[end]);end+=1
                assert end<len(rows)
                expanded.extend(body*times);at=end+1
            else:expanded.append(rows[at]);at+=1
        rows=expanded
        active=[]
        for s in rows:
            if s.startswith('.if '):active.append(bool(integer(s[4:])));continue
            if s=='.endif':assert active;active.pop();continue
            if not all(active) or not s:continue
            if s.endswith(':'):out.append(('@label',[s[:-1]]));continue
            if s.startswith('.'):continue
            op,*tail=s.split(None,1)
            if op in macros:
                names,body=macros[op];vals=[] if not tail else [v.strip() for v in tail[0].split(',')];assert len(names)==len(vals)
                mapped=[]
                for line in body:
                    for k,v in sorted(zip(names,vals),key=lambda t:-len(t[0])):line=line.replace('\\'+k,v)
                    mapped.append(line)
                expand(mapped,depth+1)
            else:out.append((op,[] if not tail else [v.strip() for v in tail[0].split(',')]))
        assert not active
    expand(ordinary);return out


RODATA={};RODATA_VALUES={}
ro_at=0x260000
ro_text=(ROOT/'canonical_b1.S').read_text().split('.section .rodata',1)[1]
for ro_line in ro_text.splitlines():
    ro_line=ro_line.split('//')[0].strip()
    if ro_line.endswith(':'):RODATA[ro_line[:-1]]=ro_at
    elif ro_line.startswith('.quad '):
        for literal in ro_line[6:].split(','):
            value=integer(literal);assert 0<=value<=MASK
            RODATA_VALUES[ro_at]=value;ro_at+=8
assert sum(RODATA_VALUES[RODATA['conjugate_p']+8*i]<<(64*i) for i in range(6))==P
assert sum(RODATA_VALUES[RODATA['conjugate_halfceil']+8*i]<<(64*i) for i in range(6))==(P+1)//2

class Machine:
    def __init__(self,coeffs):
        self.reg=[0]*32;self.mem=dict(RODATA_VALUES);self.coeffs=coeffs;self.count=0;self.calls={};self.sp=0x100000
        self.out=0x200000;self.prev1=0x210000;self.prev2=0x220000;self.const=0x230000;self.mod=0x240000
        self.table=0x408*2**32;self.block=96*(11*2**32+2**30);self.wide=57*2**39;self.vector_bytes=2*self.wide+4*self.block;self.fixed_base=0x4de8*2**32;self.fixed_bytes=48*(11*2**32+2**29)
        self.reg[31]=self.sp
        self.reg[18]=0xfeedfacedeadbeef
        for i in range(19,31):self.reg[i]=0xabc000+i
        self.saved=self.reg[18:]
    def readmem(self,at):
        assert at%8==0
        if hasattr(self,'fixed_s') and self.fixed_base<=at<self.fixed_base+self.fixed_bytes:
            off=at-self.fixed_base;start=0
            for width in [32]*11+[29]:
                size=48*(1<<width)
                if off<size:
                    digit,part=divmod(off,48);value=self.fixed_s*RI*pow(2,start,P)*digit%P
                    return (value>>(8*part))&MASK
                off-=size;start+=width
            raise AssertionError('fixed lookup range')
        if self.table<=at<self.table+self.vector_bytes:
            off=at-self.table
            if off<2*self.wide:
                col,off=divmod(off,self.wide);start=0
                for width in [35]*7+[34]*3+[35]:
                    size=96*(1<<width)
                    if off<size:
                        index,part=divmod(off,96);break
                    off-=size;start+=width
                else:raise AssertionError('wide column range')
                assert self.coeffs[col]==self.coeffs[col+4]
            else:
                slot,off=divmod(off-2*self.wide,self.block);col=[2,3,6,7][slot]
                digit,off=divmod(off,96*2**32);index,part=divmod(off,96)
                assert digit<12 and index<(2**30 if digit==11 else 2**32)
                start=32*digit
            row,limb=divmod(part,48)
            value=self.coeffs[col][row]*RI*pow(2,start,P)*(index-(2**35 if col in (0,1) and start==347 and index>=2**34 else 0))%P
            return (value>>(8*limb))&MASK
        assert at in self.mem,('uninitialized',hex(at))
        return self.mem[at]
    def writemem(self,at,value):
        assert at%8==0 and ((self.sp-816<=at<self.sp) or (self.out<=at<self.out+getattr(self,'out_bytes',192)) or any(lo<=at<hi for lo,hi in getattr(self,'extra_write_ranges',[]))),('write',hex(at))
        self.mem[at]=value&MASK
    def field(self,at,value):
        for i in range(6):self.mem[at+8*i]=(value>>(64*i))&MASK
    def loadfield(self,at):return sum(self.readmem(at+8*i)<<(64*i) for i in range(6))
    def execute(self,name):
        ops=OPS[name];carry=0;zero=0;before=self.count
        def ix(s):return 31 if s=='sp' else int(s[1:])
        def rd(s):
            if s.startswith(':lo12:'):return RODATA[s.removeprefix(':lo12:')]&4095
            value=0 if s in ('xzr','wzr') else integer(s) if s.startswith('#') else self.reg[ix(s)]
            return value&0xffffffff if s.startswith('w') else value
        def wr(s,v):
            if s not in ('xzr','wzr'):self.reg[ix(s)]=v&MASK
        def shifted(a):
            value=rd(a[2])
            if len(a)==4:
                kind,n=a[3].split();n=integer(n);assert kind in ('lsl','lsr') and 0<=n<64
                value=(value<<n)&MASK if kind=='lsl' else value>>n
            return value
        labels={a[0]:i for i,(op,a) in enumerate(ops) if op=='@label'}
        pc=labels.get(name,0)
        while pc<len(ops):
            op,a=ops[pc];pc+=1
            if op=='@label':continue
            self.count+=1;assert self.count<10000
            assert not any(re.search(r'\b[wx]18\b',v) for v in a)
            if op=='ret':self.calls.setdefault(name,[]).append(self.count-before);return
            if op=='bl':self.reg[30]=pc*4;self.execute(a[0]);carry=0;zero=0;continue
            if op=='b':pc=labels[a[0]];continue
            if op in ('cbz','cbnz'):
                take=(rd(a[0])==0)
                if op=='cbnz':take=not take
                if take:pc=labels[a[1]]
                continue
            if op.startswith('b.'):
                cond=op[2:];assert cond in ('lo','hi','eq','ne','hs')
                take={'lo':not carry,'hi':carry and not zero,'eq':bool(zero),'ne':not zero,'hs':bool(carry)}[cond]
                if take:pc=labels[a[0]]
                continue
            if op=='cmp':
                delta=rd(a[0])-rd(a[1]);carry=int(delta>=0);zero=int((delta&MASK)==0);continue
            if op in ('ldr','ldp','str','stp','ldrb','strb'):
                n=2 if op in ('ldp','stp') else 1
                addr=','.join(a[n:]).replace(' ','');m=re.fullmatch(r'\[(x\d+|sp)(?:,#(.+?))?\](!)?(?:,#(.+))?',addr);assert m,addr
                base,off,pre,post=m.groups();assert not(pre and post)
                displacement=integer(post if post else off if off else '0')
                if op in ('ldp','stp'):
                    assert displacement%8==0 and -512<=displacement<=504,('pair immediate',op,displacement)
                elif pre or post:
                    assert -256<=displacement<=255,('unscaled immediate',op,displacement)
                else:
                    if op in ('ldrb','strb'):assert 0<=displacement<=4095,('byte immediate',op,displacement)
                    else:assert displacement%8==0 and 0<=displacement<=32760,('scaled immediate',op,displacement)
                old=rd(base);at=old+(integer(off) if off else 0)
                if pre:wr(base,at)
                if op in ('ldrb','strb'):
                    word_at=at&~7;shift=(at&7)*8;word=self.readmem(word_at)
                    if op=='ldrb':wr(a[0],(word>>shift)&255)
                    else:self.writemem(word_at,(word&~(255<<shift))|((rd(a[0])&255)<<shift))
                    if post:wr(base,old+integer(post))
                    continue
                if op in ('ldr','ldp'):
                    vals=[self.readmem(at+8*j) for j in range(n)]
                    for j,v in enumerate(vals):wr(a[j],v)
                else:
                    vals=[rd(a[j]) for j in range(n)]
                    for j,v in enumerate(vals):self.writemem(at+8*j,v)
                if post:wr(base,old+integer(post))
                continue
            if op=='adrp':value=RODATA[a[1]]&~4095
            elif op=='mov':value=rd(a[1])
            elif op=='movz':
                value=rd(a[1]);assert len(a) in (2,3)
                if len(a)==3:
                    kind,n=a[2].split();assert kind=='lsl';value<<=integer(n)
            elif op=='ubfx':value=(rd(a[1])>>rd(a[2]))&((1<<rd(a[3]))-1)
            elif op=='bfi':
                lsb,width=rd(a[2]),rd(a[3]);assert 0<=lsb<64 and 1<=width<=64-lsb
                mask=((1<<width)-1)<<lsb
                value=(rd(a[0])&~mask)|((rd(a[1])<<lsb)&mask)
            elif op=='and':value=rd(a[1])&shifted(a)
            elif op=='orr':value=rd(a[1])|shifted(a)
            elif op=='lsr':value=rd(a[1])>>rd(a[2])
            elif op=='madd':value=rd(a[1])*rd(a[2])+rd(a[3])
            elif op in ('mul','umulh'):
                value=rd(a[1])*rd(a[2]);value=value>>64 if op=='umulh' else value
            elif op=='extr':value=((rd(a[1])<<64)|rd(a[2]))>>rd(a[3])
            elif op=='udiv':assert rd(a[2]);value=rd(a[1])//rd(a[2])
            elif op in ('add','adds','adc','adcs'):
                value=rd(a[1])+shifted(a)+(carry if op in ('adc','adcs') else 0)
                if op in ('adds','adcs'):carry=int(value>MASK);zero=int((value&MASK)==0)
            elif op in ('sub','subs','sbcs'):
                value=rd(a[1])-shifted(a)-(1-carry if op=='sbcs' else 0)
                if op!='sub':carry=int(value>=0);zero=int((value&MASK)==0)
            elif op=='csel':
                assert a[3] in ('lo','hs');yes=bool(carry) if a[3]=='hs' else not carry;value=rd(a[1] if yes else a[2])
            else:raise AssertionError(('opcode',op))
            wr(a[0],value)
        raise AssertionError('missing ret')

OPS={n:parse(n+'.S') for n in ['recurrence_fp2_step','recurrence_fp2_map','fp381_add','fp381_sub','nh_fingerprint','fixed_s_lookup','canonical_twelve_bfirst','odd_trace_compact','fr255_add']}
OPS['fill_pair_offset']=parse('fill_pair.S')
OPS['fill_pair']=parse('fill_pair.S')

def own_acyclic_bound(name):
    # All branch outcomes,including infeasible ones. BL counts one own instruction.
    ops=OPS[name];labels={a[0]:i for i,(op,a) in enumerate(ops) if op=='@label'}
    memo={};visiting=set()
    def visit(i):
        assert 0<=i<len(ops)
        if i in memo:return memo[i]
        assert i not in visiting,'cyclic parent control flow'
        visiting.add(i);op,a=ops[i]
        if op=='@label':v=visit(i+1)
        elif op=='ret':v=1
        elif op=='b':v=1+visit(labels[a[0]])
        elif op.startswith('b.') or op in ('cbz','cbnz'):
            v=1+max(visit(i+1),visit(labels[a[-1]]))
        else:v=1+visit(i+1)
        visiting.remove(i);memo[i]=v;return v
    return visit(labels[name])

PARENT_OWN_BOUND=own_acyclic_bound('canonical_twelve_bfirst')
assert PARENT_OWN_BOUND<=94,('parent own bound',PARENT_OWN_BOUND)

OPS['canonical_conjugate_b1']=parse('canonical_b1.S')

for io_name in ['prefix_copy4','prefix_copy6','prefix_copy384','prefix_copy386','prefix_copy512','metadata_pending_get','metadata_pending_set','metadata_pending_clear','metadata_ordinal','metadata_color','baby_scalar_fill']:OPS[io_name]=parse('metadata_io.S')

def cadd(a,b):return tuple((x+y)%P for x,y in zip(a,b))
def cneg(a):return tuple(-x%P for x in a)
def cmul(a,b):return ((a[0]*b[0]-a[1]*b[1])%P,(a[0]*b[1]+a[1]*b[0])%P)
def oracle(c,x,y,z):
    A,B=c[:2],c[2:];a,b=x[:2],x[2:];d,e=y[:2],y[2:]
    # Independent dense Fp4 products, not the factored native transform.
    def mul(v,w):
        aa,bb=v[:2],v[2:];cc,dd=w[:2],w[2:]
        return cadd(cmul(aa,cc),cmul((1,1),cmul(bb,dd)))+cadd(cmul(aa,dd),cmul(bb,cc))
    left=mul(c,x);right=mul(A+cneg(B),y)
    return tuple((left[i]-right[i]+z[i])%P for i in range(4))

def run(inputs):
    count=inputs['count'];assert 1<=count<=16;rng=random.Random(inputs['seed']);results=[]
    for case in range(count):
        # Fixed full-p boundaries plus full-p random cases; no reduced field.
        vals=[rng.randrange(P) for _ in range(16)]
        if case==0:vals=[0,P-1,1,(P-1)//2]*4
        if case==1:vals=[P-1]*16
        if case==2:vals=[0]*16
        if case==3:vals=[(P-1)*RI%P]*16
        if case==4:vals[4:8]=[(P-1)*RI%P]*4;vals[8:12]=[0]*4
        if case==5:vals[4:8]=[0]*4;vals[8:12]=[(P-1)*RI%P]*4
        c,x,y,z=[tuple(vals[i:i+4]) for i in range(0,16,4)]
        A,B=c[:2],c[2:];Bxi=cmul(B,(1,1));coeffs=[]
        for C in (A,Bxi,A,B):coeffs.extend([(C[0]*R%P,C[1]*R%P),(-C[1]*R%P,C[0]*R%P)])
        machine=Machine(coeffs)
        for base,vec in [(machine.out,z),(machine.prev1,x),(machine.prev2,y)]:
            for j,v in enumerate(vec):machine.field(base+48*j,v*R%P)
        for at in range(machine.const,machine.const+3728,8):machine.mem[at]=0
        machine.field(machine.const+448,P);machine.mem[machine.const+512]=P//2**325;machine.field(machine.mod,P)
        for quotient in range(50):
            for limb in range(7):machine.mem[machine.const+528+64*quotient+8*limb]=(quotient*P>>(64*limb))&MASK
        machine.reg[:6]=[machine.out,machine.prev1,machine.prev2,machine.table,machine.const,machine.mod]
        machine.execute('recurrence_fp2_step')
        rawx=tuple(v*R%P for v in x);rawy=tuple(v*R%P for v in y)
        expected_temp=tuple((rawx[j]-rawy[j])%2**384 for j in (0,1))+tuple(rawx[j]+rawy[j] for j in (2,3))+tuple((rawx[j]-rawy[j])%2**384 for j in (2,3))+tuple(rawx[j]+rawy[j] for j in (0,1))
        actual_temp=tuple(machine.loadfield(machine.sp-384+48*j) for j in range(8))
        assert actual_temp==expected_temp,('raw transform integers',case)
        for at,j in zip((0,1,4,5),(0,1,2,3)):
            code=actual_temp[at]%2**382;signed=code-2**382 if code>=2**381 else code
            assert signed==rawx[j]-rawy[j]
            assert actual_temp[at]>>382==(3 if signed<0 else 0)
        assert all(actual_temp[j]<2*P<2**382 for j in (2,3,6,7))
        if case==3:assert all(actual_temp[j]==2*P-2 and actual_temp[j]>=2**381 for j in (2,3,6,7))
        if case==4:assert all(actual_temp[j]==P-1 and actual_temp[j]<2**381 for j in (0,1,4,5))
        if case==5:assert all(actual_temp[j]==2**384-P+1 and actual_temp[j]%2**382>=2**381 for j in (0,1,4,5))
        got=tuple(machine.loadfield(machine.out+48*j)*RI%P for j in range(4));assert got==oracle(c,x,y,z),(case,got)
        assert all(machine.loadfield(machine.out+48*j)<P for j in range(4))
        assert machine.reg[18:]==machine.saved
        assert machine.count==2485 and not machine.calls.get('recurrence_fp2_map'),('native_counts',machine.count,machine.calls.get('recurrence_fp2_map'))
        assert all(v==31 for n in ('fp381_add','fp381_sub') for v in machine.calls.get(n,[]))
        results.append({'case':case,'step_instructions':machine.count,'fused_leaf':True,'exact_third_alias':True,'abi_and_initialized_reads':True})
    return {'scope':'full-p native fused recurrence map/step instruction interpretation against independent dense Fp4 oracle; virtual read-only table test fixture only; no production fill, GT coverage, hardware or score qualification','cases':results}


def run_offset_fill(inputs):
    rng=random.Random(inputs['seed']+381);cases=[]
    for case in range(inputs['count']):
        n=2+case%4;seed=(rng.randrange(P),rng.randrange(P))
        if case==0:seed=(0,0)
        if case==1:seed=(P-1,P-1)
        for name,offset in [('fill_pair',(0,0)),('fill_pair_offset',tuple(-v*2**34%P for v in seed))]:
            machine=Machine([]);machine.out_bytes=96*n
            for j in range(2):
                machine.field(machine.prev1+48*j,seed[j]);machine.field(machine.prev2+48*j,offset[j])
            machine.field(machine.mod,P)
            machine.reg[:5]=[machine.out,machine.prev1,machine.mod,n,machine.prev2]
            machine.execute(name)
            for digit in range(n):
                got=tuple(machine.loadfield(machine.out+96*digit+48*j) for j in range(2))
                assert got==tuple((offset[j]+digit*seed[j])%P for j in range(2))
                assert all(v<P for v in got)
            assert machine.reg[18:]==machine.saved
            assert machine.count<256*n+128
            cases.append({'case':case,'helper':name,'entries':n,'instructions':machine.count,'canonical_progression':True,'abi_and_initialized_reads':True})
    return {'scope':'actual scalar zero/offset-fill sources on bounded2..5 full-p entries vs independent progression;production2^34 half-fills/construction/cost are analytical/paid,not virtual/free tables','cases':cases}


def run_nh(inputs):
    rng=random.Random(inputs['seed']+72);q=2**61-1;cases=[]
    for case in range(inputs['count']):
        words=[rng.getrandbits(64) for _ in range(24)]
        if case==0:words=[0]*24
        if case==1:words=[MASK]*24
        keys=[[rng.getrandbits(64) for _ in range(24)] for _ in range(3)]
        affine=[[rng.randrange(q) for _ in range(4)] for _ in range(3)]
        if case==0:affine=[[q-1]*4 for _ in range(3)]
        if case==1:affine=[[0,q-1,q-1,q-1] for _ in range(3)]
        machine=Machine([]);base=machine.out if case%2 else machine.prev1
        for j,v in enumerate(words):machine.mem[base+8*j]=v
        for pair in range(12):
            for row in range(3):
                for j in range(2):machine.mem[machine.const+48*pair+16*row+8*j]=keys[row][2*pair+j]
        for row in range(3):
            for j,v in enumerate(affine[row]):machine.mem[machine.const+576+32*row+8*j]=v
        machine.reg[:3]=[machine.out,base,machine.const];before=dict(machine.mem)
        machine.execute('nh_fingerprint')
        ys=[]
        for row in range(3):
            nh=sum(((words[2*i]+keys[row][2*i])&MASK)*((words[2*i+1]+keys[row][2*i+1])&MASK) for i in range(12))%(1<<128)
            chunks=[nh&((1<<60)-1),(nh>>60)&((1<<60)-1),nh>>120]
            y=(affine[row][0]+sum(affine[row][j+1]*chunks[j] for j in range(3)))%q
            ys.append(y)
        want=ys[0]|(ys[1]<<61)|(ys[2]<<122)
        got=sum(machine.readmem(machine.out+8*j)<<(64*j) for j in range(4));assert got==want
        assert machine.count==341 and machine.reg[18:]==machine.saved
        assert all(machine.out<=at<machine.out+32 or before.get(at)==v for at,v in machine.mem.items())
        cases.append({'case':case,'instructions':341,'exact_input_alias':base==machine.out,'abi_and_initialized_reads':True})
    return {'scope':'actual encodable scalar NH plus affine-compression instructions against independent integer oracle; no empirical universal-hash probability or hardware/score qualification','cases':cases}


def run_fixed(inputs):
    rng=random.Random(inputs['seed']+77);power=(1,1);e=(P-1)//2;C=(1,0)
    while e:
        if e&1:C=cmul(C,power)
        power=cmul(power,power);e>>=1
    assert C[0]==C[1] and 2*C[0]*C[0]%P==P-1
    scalar=C[0]*R%P;cases=[]
    for case in range(inputs['count']):
        x=[0,P-1,(P+1)//2][case] if case<3 else rng.randrange(P)
        machine=Machine([]);machine.fixed_s=scalar;base=machine.out if case%2 else machine.prev1
        machine.field(base,x)
        for at in range(machine.const,machine.const+3728,8):machine.mem[at]=0
        machine.field(machine.const+448,P);machine.mem[machine.const+512]=P//2**325
        for q in range(50):
            for i in range(7):machine.mem[machine.const+528+64*q+8*i]=(q*P>>(64*i))&MASK
        machine.reg[:4]=[machine.out,base,machine.fixed_base,machine.const];before=dict(machine.mem)
        machine.execute('fixed_s_lookup');got=machine.loadfield(machine.out)
        assert got==x*scalar*RI%P and got<P
        assert machine.count==200 and machine.reg[18:]==machine.saved
        assert all(machine.out<=at<machine.out+48 or machine.sp-96<=at<machine.sp or before.get(at)==v for at,v in machine.mem.items())
        cases.append({'case':case,'instructions':200,'exact_input_alias':base==machine.out,'abi_and_initialized_reads':True})
    return {'scope':'actual wider-digit constant lookup plus small-quotient native correction vs independent full-p product for actual fixed s; virtual initialized table fixture only, not production fill or GT/key probability qualification','cases':cases}


def run_key(inputs):
    rng=random.Random(inputs['seed']+83);power=(1,1);e=(P-1)//2;C=(1,0)
    while e:
        if e&1:C=cmul(C,power)
        power=cmul(power,power);e>>=1
    assert C[0]==C[1] and 2*C[0]*C[0]%P==P-1
    scalar=C[0]*R%P;cases=[]
    def sign(B):
        first=B[0] if B[0] else B[1]
        return tuple(-x%P for x in B) if first and (first&MASK)>((P-first)&MASK) else B
    def words(B):return tuple((v>>(64*j))&MASK for v in B for j in range(6))
    for case in range(inputs['count']):
        phi=[rng.randrange(P) for _ in range(4)]
        if case==0:phi=[R,0,0,0]
        if case==1:phi=[0,P-1,0,P-1]
        if case==2:phi=[P-1,(P+1)//2,(P+1)//2,0]
        if case and phi[2:]==[0,0]:phi[2]=1
        machine=Machine([]);machine.fixed_s=scalar
        for j,v in enumerate(phi):machine.field(machine.prev1+48*j,v)
        for at in range(machine.const,machine.const+3728,8):machine.mem[at]=0
        machine.field(machine.const+448,P);machine.mem[machine.const+512]=P//2**325;machine.field(machine.mod,P)
        for q in range(50):
            for j in range(7):machine.mem[machine.const+528+64*q+8*j]=(q*P>>(64*j))&MASK
        machine.reg[:5]=[machine.out,machine.prev1,machine.fixed_base,machine.mod,machine.const]
        machine.execute('canonical_twelve_bfirst')
        assert machine.reg[18:]==machine.saved and machine.count<=710
        assert [machine.loadfield(machine.prev1+48*j) for j in range(4)]==phi
        assert all(v<=46 for v in machine.calls['canonical_conjugate_b1'])
        assert all(v<=523 for v in machine.calls['odd_trace_compact'])
        assert machine.calls['fixed_s_lookup']==[200,200]
        if case==0:
            assert machine.reg[0]==1 and not any(machine.out<=at<machine.out+192 for at in machine.mem)
        else:
            even=sign(tuple(phi[2:]));Br,Bi=phi[2:];odd=sign((C[0]*(Br+Bi)%P,C[0]*(Br-Bi)%P))
            assert even!=odd
            expected=(tuple(phi[:2])+even) if words(even)<words(odd) else ((phi[0],-phi[1]%P)+odd)
            got=tuple(machine.loadfield(machine.out+48*j) for j in range(4))
            assert machine.reg[0]==0 and got==expected
        cases.append({'case':case,'instructions':machine.count,'key_bound':710,'immutable_raw_phi':True,'abi_and_initialized_reads':True,'identity_case':case==0})
    return {'scope':'actual complete key/sign/odd/lookup scalar instructions against independent full-p field/canonical-order formula; nonzero-b1 field cases plus true GT identity only; not empirical GT twelve-fiber theorem, production inversion, probability, hardware or qualification','cases':cases}

# CANDIDATE TEST TEXT: only the credential-free shared verifier may execute.
# Sparse dictionaries are bounded initialized-zero test fixtures, not production
# replacements for the explicitly constructed/scanned 2^80 bucket arrays.

def add128(a,b):
    assert 0<=a<1<<128 and 0<=b<1<<128
    lo=(a&MASK)+(b&MASK);hi=(a>>64)+(b>>64)+(lo>>64)
    assert hi<1<<64
    return (hi<<64)|(lo&MASK)

def digit_radix(raw,stage):
    w0,w1,w2=[int.from_bytes(raw[8*j:8*j+8],'little') for j in range(3)]
    U=int.from_bytes(raw[:24],'little');assert U<1<<183
    if stage==0:
        high=(w1>>1)&((1<<31)-1);packed=w0|((w1&1)<<64)
    else:
        high=(w2>>33)&((1<<31)-1);packed=(((w2<<64)|w1)>>32)&((1<<65)-1)
    middle=packed>>33;low=packed&((1<<33)-1)
    assert (high<<65)|(middle<<33)|low==(U>>(96*stage))&((1<<96)-1)
    assert high<1<<31 and middle<1<<32 and low<1<<33 and 2+high<1<<32
    return high,middle,low,packed

def meta_fixture(U,ordinal,color):
    assert 0<=U<1<<183 and 0<=ordinal<1<<128 and color in (0,1)
    return U.to_bytes(24,'little')+bytes([color])+bytes(7)+ordinal.to_bytes(16,'little')

def packet_fixture(at,data):
    assert len(data) in (48,3072)
    result=at.to_bytes(16,'little')+data
    assert len(result)+32<=6144
    if len(data)==3072:
        assert at%64==0 and at//(1<<34)==(at+63)//(1<<34)
    return result

class LeafRadix:
    def __init__(self,stage,high,middle):
        self.stage=stage;self.high=high;self.middle=middle
        self.counts={};self.next={};self.pending={}
    def count(self,words):
        assert 0<len(words)<=512
        for low in words:
            assert 0<=low<1<<33
            self.counts[low]=add128(self.counts.get(low,0),1)
    def prefix(self,start):
        for low in sorted(self.counts):
            self.next[low]=start;start=add128(start,self.counts[low])
        return start
    def receive(self,records):
        assert 0<len(records)<=64
        queue=[]
        for raw in records:
            high,middle,low,_=digit_radix(raw,self.stage)
            assert (high,middle)==(self.high,self.middle) and raw[25:32]==bytes(7)
            at=self.next[low];self.next[low]=add128(at,1)
            pending=self.pending.get(low,[])
            if pending:assert pending[0][25]==len(pending)
            if not pending and at%64:
                queue.append(packet_fixture(at,raw));continue
            pending.append(bytearray(raw))
            if len(pending)<64:
                pending[0][25:26]=len(pending).to_bytes(1,'little')
                self.pending[low]=pending
            else:
                assert len(pending)==64
                pending[0][25:26]=bytes(1)
                queue.append(packet_fixture(at-63,b''.join(pending)))
                self.pending.pop(low,None)
        assert len(queue)<=len(records)<=64
        return queue
    def finish(self):
        for low in sorted(self.pending):
            records=self.pending[low];count=len(records)
            assert 1<=count<=63 and records[0][25]==count
            records[0][25:26]=bytes(1)
            start=self.next[low]-count
            for j,raw in enumerate(records):yield packet_fixture(add128(start,j),bytes(raw))
        self.pending.clear()

class Relay96:
    def __init__(self,stage,high):
        self.stage=stage;self.high=high
        self.leaves={};self.totals={};self.count_cache={};self.scatter_cache={};self.queue=[];self.maximum=0
    def leaf(self,middle):
        if middle not in self.leaves:self.leaves[middle]=LeafRadix(self.stage,self.high,middle)
        return self.leaves[middle]
    def count(self,words):
        assert 0<len(words)<=512
        for word in words:
            middle,low=word>>33,word&((1<<33)-1)
            self.totals[middle]=add128(self.totals.get(middle,0),1)
            cache=self.count_cache.setdefault(middle,[]);cache.append(low)
            if len(cache)==512:self.leaf(middle).count(cache);cache.clear()
    def count_finish(self):
        for middle,cache in self.count_cache.items():
            if cache:self.leaf(middle).count(cache);cache.clear()
    def prefix(self,start):
        for middle in sorted(self.totals):
            end=self.leaf(middle).prefix(start)
            assert end==add128(start,self.totals[middle]);start=end
        return start
    def receive(self,records):
        assert 0<len(records)<=64 and not self.queue
        for raw in records:
            high,middle,_,_=digit_radix(raw,self.stage);assert high==self.high
            cache=self.scatter_cache.setdefault(middle,[]);cache.append(raw)
            if len(cache)==64:self.queue.extend(self.leaf(middle).receive(cache));cache.clear()
        self.maximum=max(self.maximum,len(self.queue));assert len(self.queue)<=4096

def sort_radix_fixture(records,stage):
    nodes={};totals={};count_cache={};scatter_cache={};out={};full_count_frames=0
    entry_bytes=16
    count_capacity=4096//entry_bytes
    def node(high):
        if high not in nodes:nodes[high]=Relay96(stage,high)
        return nodes[high]
    def consume_queue(r):
        for packet in r.queue:
            at=int.from_bytes(packet[:16],'little');data=packet[16:]
            assert len(data) in (48,3072)
            for j in range(len(data)//48):
                dest=add128(at,j);raw=data[48*j:48*j+48]
                assert 0<=dest<len(records) and dest not in out and raw[25:32]==bytes(7)
                out[dest]=raw
        r.queue.clear()
    for raw in records:
        high,_,_,word=digit_radix(raw,stage);totals[high]=add128(totals.get(high,0),1)
        cache=count_cache.setdefault(high,[]);cache.append(word)
        if len(cache)==count_capacity:
            frame=b''.join(x.to_bytes(entry_bytes,'little') for x in cache);assert len(frame)==4096
            node(high).count([int.from_bytes(frame[entry_bytes*j:entry_bytes*j+entry_bytes],'little') for j in range(count_capacity)])
            full_count_frames+=1;cache.clear()
    for high,cache in count_cache.items():
        if cache:node(high).count(cache);cache.clear()
    for r in nodes.values():r.count_finish()
    start=0
    for high in sorted(totals):
        end=node(high).prefix(start);assert end==add128(start,totals[high]);start=end
    assert start==len(records)
    for raw in records:
        high,_,_,_=digit_radix(raw,stage);cache=scatter_cache.setdefault(high,[]);cache.append(raw)
        if len(cache)==64:node(high).receive(cache);consume_queue(node(high));cache.clear()
    for high,cache in scatter_cache.items():
        if cache:node(high).receive(cache);consume_queue(node(high));cache.clear()
    for r in nodes.values():
        for middle,cache in r.scatter_cache.items():
            if cache:
                r.queue.extend(r.leaf(middle).receive(cache));r.maximum=max(r.maximum,len(r.queue))
                assert len(r.queue)<=4096;consume_queue(r);cache.clear()
        for leaf in r.leaves.values():
            for packet in leaf.finish():
                r.queue.append(packet)
                if len(r.queue)==64:consume_queue(r)
            consume_queue(r)
    assert len(out)==len(records)
    return [out[i] for i in range(len(records))],max(r.maximum for r in nodes.values()),full_count_frames

def run_protocol(inputs):
    rng=random.Random(inputs['seed']+85);results=[]
    for case in range(inputs['count']):
        records=[];current=0;baby_records=1536
        for i in range(4096):
            if case%4==0:U=0
            elif case%4==1:
                if i%4==0:current=rng.getrandbits(183)
                U=current
            elif case%4==2:U=((i%2)<<182)|((i%2)<<160)|((i%16)<<129)|((i%7)<<96)|((i%2)<<95)|((i%16)<<65)|((i%7)<<33)|((i%2)<<32)|(i%7)
            else:U=((i%3)<<144)|((i%64)<<112)|((i//64%7)<<80)|((i%2)<<64)|((i%64)<<32)|(i//64%7)
            records.append(meta_fixture(U,(1<<100)+i,0 if i<baby_records else 1))
        reconstructed=[];batches=0
        for start in range(0,len(records),512):
            batch=records[start:start+512];assert len(batch)==512
            payload=b''.join(raw[32:48]*2 for raw in batch);metadata=b''.join(batch);pending=True
            assert len(payload)==16384 and len(metadata)==24576
            payload_frames=[payload[j:j+4096] for j in range(0,len(payload),4096)]
            metadata_frames=[metadata[j:j+4096] for j in range(0,len(metadata),4096)]
            assert len(payload_frames)==4 and len(metadata_frames)==6
            assert b''.join(payload_frames)==payload and pending
            reply=b''.join(metadata_frames);pending=False
            assert reply==metadata and not pending
            reconstructed.extend(reply[j:j+48] for j in range(0,len(reply),48));batches+=1
        assert reconstructed==records and batches==8
        assert sum(raw[24]==0 for raw in reconstructed)==baby_records
        assert sum(raw[24]==1 for raw in reconstructed)==4096-baby_records
        for group in range(0,len(records),2048):
            bulk=b''.join(reconstructed[group:group+2048])
            frames=[bulk[j:j+4096] for j in range(0,len(bulk),4096)]
            assert len(bulk)==98304 and len(frames)==24 and b''.join(frames)==bulk
        stage0,q0,f0=sort_radix_fixture(reconstructed,0)
        stage1,q1,f1=sort_radix_fixture(stage0,1)
        expected=sorted(records,key=lambda raw:int.from_bytes(raw[:24],'little'))
        assert stage1==expected
        # Standalone full-width destination carry and bank-boundary case.
        boundary=LeafRadix(0,0,0);boundary.count([0]*200)
        base=(1<<100)+(1<<34)-17;assert boundary.prefix(base)==add128(base,200)
        test=[meta_fixture(0,(1<<110)+i,i&1) for i in range(200)];packets=[]
        for i in range(0,200,64):packets.extend(boundary.receive(test[i:i+64]))
        packets.extend(boundary.finish());written={}
        for packet in packets:
            at=int.from_bytes(packet[:16],'little');data=packet[16:]
            for j in range(len(data)//48):written[add128(at,j)]=data[48*j:48*j+48]
        assert [written[add128(base,j)] for j in range(200)]==test
        assert add128((1<<64)-1,1)==1<<64
        results.append({'case':case,'records':4096,'baby_source_records':baby_records,'giant_source_records':4096-baby_records,'worker_batches':8,'payload_bytes':32,'payload_frames_per_batch':4,'metadata_frames_per_batch':6,'stable_full183_sort':True,'encoded192_padding_zero':True,'digit_widths':[96,96],'root_count_entry_bytes':[16,16],'root_count_frame_records':[256,256],'full128_ordinals_preserved':True,'padding_restored':True,'maximum_relay_queue':max(q0,q1),'full_count_frames':f0+f1,'bank_and_counter_boundary':True})
    return {'scope':'bounded sparse initialized-zero packet/state-model fixture at full183-bit fingerprints in192-bit encoding/full128-bit ordinals and real48-byte metadata/32-byte payloads; two96-bit two-edge passes and6 metadata frames; not production bucket initialization, native dispatch timing, GT success, probability or qualification','cases':results}

def run_prefix_copy(inputs):
    rng=random.Random(inputs['seed']+86);cases=[]
    for case in range(inputs['count']):
        for size in (4,6,384,386,512):
            machine=Machine([]);machine.out_bytes=size*8
            base=machine.out if case%2 else machine.prev1;words=[rng.getrandbits(64) for _ in range(size)]
            for j,v in enumerate(words):machine.mem[base+8*j]=v
            machine.reg[:2]=[machine.out,base];before=dict(machine.mem)
            machine.execute('prefix_copy'+str(size))
            assert [machine.readmem(machine.out+8*j) for j in range(size)]==words
            assert [machine.readmem(base+8*j) for j in range(size)]==words
            assert machine.count==size+(size%2)+1 and machine.reg[18:]==machine.saved
            assert machine.reg[0]==machine.out+size*8 and machine.reg[1]==base+size*8
            assert all(machine.out<=at<machine.out+size*8 or before.get(at)==v for at,v in machine.mem.items())
            cases.append({'case':case,'words':size,'instructions_including_ret':size+(size%2)+1,'exact_alias':base==machine.out,'abi_and_initialized_reads':True})
    return {'scope':'actual post-indexed scalar4/6/384/386/512-word frame-copy source, including encodable immediates, full bytes, alias/ABI/read/write/count constraints; not a production routing/controller timing certificate','cases':cases}


def run_pending_byte(inputs):
    rng=random.Random(inputs['seed']+87);cases=[]
    for case in range(inputs['count']):
        machine=Machine([]);machine.out_bytes=48
        raw=bytearray(rng.getrandbits(8) for _ in range(48));raw[24]=case&1;raw[25:32]=bytes(7)
        for j in range(6):machine.mem[machine.out+8*j]=int.from_bytes(raw[8*j:8*j+8],'little')
        q=[1,63,32,0][case%4];machine.reg[:2]=[machine.out,q]
        machine.execute('metadata_pending_set');assert machine.count==2 and machine.reg[18:]==machine.saved
        actual=b''.join(machine.readmem(machine.out+8*j).to_bytes(8,'little') for j in range(6));expected=bytearray(raw);expected[25]=q
        assert actual==expected and actual[24]==raw[24] and actual[32:48]==raw[32:48]
        machine.reg[0]=machine.out;machine.execute('metadata_color');assert machine.reg[0]==raw[24]
        machine.reg[0]=machine.out;machine.execute('metadata_ordinal')
        assert machine.reg[0]|(machine.reg[1]<<64)==int.from_bytes(raw[32:48],'little')
        machine.reg[0]=machine.out;machine.execute('metadata_pending_get');assert machine.reg[0]==q
        machine.reg[0]=machine.out;machine.execute('metadata_pending_clear')
        cleared=b''.join(machine.readmem(machine.out+8*j).to_bytes(8,'little') for j in range(6));assert cleared==raw
        assert machine.count==10 and machine.reg[18:]==machine.saved
        cases.append({'case':case,'occupancy':q,'only_byte25_changes':True,'color_and_full128_ordinal_preserved':True,'byte_and_ordinal_wrappers_instructions':2})
    return {'scope':'actual unsigned-offset LDRB/STRB occupancy primitives against complete48-byte row oracle; production inlines one byte instruction,wrapper RET is separately control; not a production controller timing or probability certificate','cases':cases}


def run_scalar_streams(inputs):
    order=int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001',16)
    baby_count=34021169017091701373279488505416777728
    assert baby_count<1<<126 and baby_count<order and 2*order<1<<256
    rng=random.Random(inputs['seed']+91);cases=[]
    for case in range(inputs['count']):
        j=[1,baby_count,(1<<64)-1,1<<64][case%4]
        machine=Machine([]);machine.out_bytes=32;machine.extra_write_ranges=[(machine.prev1,machine.prev1+16)]
        machine.mem[machine.prev1]=j&MASK;machine.mem[machine.prev1+8]=j>>64
        machine.reg[:2]=[machine.out,machine.prev1];before=dict(machine.mem)
        machine.execute('baby_scalar_fill')
        got=sum(machine.readmem(machine.out+8*k)<<(64*k) for k in range(4))
        state=machine.readmem(machine.prev1)|(machine.readmem(machine.prev1+8)<<64)
        assert got==j and got<order and state==j+1 and machine.count==7 and machine.reg[18:]==machine.saved
        assert all(machine.out<=at<machine.out+32 or machine.prev1<=at<machine.prev1+16 or before.get(at)==v for at,v in machine.mem.items())
        a,b=[(0,0),(order-1,1),(order-1,order-1),(rng.randrange(order),rng.randrange(order))][case%4]
        giant=Machine([]);giant.out_bytes=32;ap=giant.out if case%2 else giant.prev1
        for base,value in [(ap,a),(giant.prev2,b),(giant.mod,order)]:
            for k in range(4):giant.mem[base+8*k]=(value>>(64*k))&MASK
        giant.reg[:4]=[giant.out,ap,giant.prev2,giant.mod]
        giant.execute('fr255_add')
        got=sum(giant.readmem(giant.out+8*k)<<(64*k) for k in range(4))
        assert got==(a+b)%order and got<order and giant.count==21 and giant.reg[18:]==giant.saved
        cases.append({'case':case,'baby_j':str(j),'baby_instructions':7,'full32_encoding_and_exact128_next':True,'giant_instructions':21,'giant_full_r_oracle':True,'giant_exact_alias':ap==giant.out})
    return {'scope':'actual small public baby counter/full32 Fr encoding and full-r giant modular-add instructions with boundaries,carry,alias,read/write/ABI/count checks;no reduced production order,GT success,probability or qualification assertion','cases':cases}

def run_scalar_payload32(inputs):
    rng=random.Random(inputs['seed']+183);cases=[]
    order=int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001',16)
    for case in range(inputs['count']):
        ordinal=(1<<100)+(1<<32)-1+case;scalar=rng.randrange(order)
        raw=meta_fixture(rng.getrandbits(183),ordinal,case&1);payload=scalar.to_bytes(32,'little')
        data_id=(1<<94)+(ordinal>>32);offset=32*(ordinal&((1<<32)-1))
        assert data_id<1<<96 and offset+32<=32*(1<<32)<1<<40
        assert int.from_bytes(payload,'little')==scalar
        assert raw[24]==case&1 and int.from_bytes(raw[32:48],'little')==ordinal
        assert raw[25:32]==bytes(7) and int.from_bytes(raw[:24],'little')<1<<183
        start=ordinal&~511;assert start>>32==(start+511)>>32
        assert (32*(start&((1<<32)-1)))%4096==0
        stream=b''.join(j.to_bytes(32,'little') for j in range(512))
        frames=[stream[j:j+4096] for j in range(0,len(stream),4096)]
        assert len(stream)==16384 and len(frames)==4 and b''.join(frames)==stream
        cases.append({'case':case,'full_scalar32':True,'no_stored_full_key':True,'color_from_metadata48':True,'full128_ordinal':True,'payload_frames':4,'bank_boundary':True})
    return {'scope':'bounded scalar-only32-byte payload/full48-byte metadata183-bit fingerprint fixture;not production GT verification,physical addressing,probability or qualification','cases':cases}

def run_verified_candidates(inputs):
    rng=random.Random(inputs['seed']+6144);cases=[]
    order=int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001',16)
    phases=[pow(P,i,order) for i in range(12)];assert phases[0]==1 and all(x!=1 for x in phases[1:]) and len(set(phases))==12 and pow(P,12,order)==1
    def candidates(pair):
        color,x,y=pair
        if color=='BB' or (color=='GG' and x==y):return []
        if color=='BG':return [(y+v*x)%order for v in phases]
        return [((y-v*x)*pow(1-v,-1,order))%order for v in phases[1:]]
    def verify_stream(pairs,target,cap):
        checks=0;processed=0
        for pair in pairs:
            if processed==cap:return None,checks,True
            processed+=1
            for value in candidates(pair):
                checks+=1
                # Test-only exact exponent oracle; production performs h^value==z.
                if value==target:return value,checks,False
        return None,checks,False
    for case in range(inputs['count']):
        b=rng.randrange(1,order);j=case+1;bi=b*(case+17)%order
        target=(bi+phases[case%12]*j)%order
        bad=j+1
        while target in candidates(('BG',bad,bi)):bad+=1
        pairs=[('BB',j,bad),('BG',bad,bi),('BG',j,bi)]
        result,checks,aborted=verify_stream(pairs,target,3)
        assert result==target and checks<=24 and not aborted
        assert verify_stream(pairs,target,2)[2]
        bt=b*(case+19)%order;v=phases[1+case%11]
        target2=(bt-v*bi)*pow(1-v,-1,order)%order
        wrong=bt+1
        while target2 in candidates(('GG',bi,wrong)):wrong+=1
        pairs2=[('GG',bi,bi),('GG',bi,wrong),('GG',bi,bt)]
        result2,checks2,aborted2=verify_stream(pairs2,target2,3)
        assert result2==target2 and checks2<=22 and not aborted2
        assert verify_stream(pairs2,target2,2)[2]
        assert 20000*6144+(1<<25)<1<<28<1<<29
        cases.append({'case':case,'full_r_mixed_and_giant_formulas':True,'continues_after_false_fingerprint':True,'BB_equalGG_skipped':True,'cap_applies_to_every_pair':True,'at_most12_phases_per_pair':True})
    return {'scope':'full-production-r algebra/control fixture using a TEST-ONLY exact exponent oracle;not production discrete logs or a GT oracle,not empirical success or instruction/cost certification','cases':cases}

for line in sys.stdin:
    if line.strip():
        try:
            inputs=json.loads(line)
            report={'exact_ledger':check_exact_ledger(),'offset_fill':run_offset_fill(inputs),'recurrence':run(inputs),'fingerprint':run_nh(inputs),'constant_lookup':run_fixed(inputs),'orbit_key':run_key(inputs),'prefix_copy':run_prefix_copy(inputs),'radix192_protocol':run_protocol(inputs),'metadata_pending_byte':run_pending_byte(inputs),'scalar_streams':run_scalar_streams(inputs),'scalar_payload32':run_scalar_payload32(inputs),'verified_candidates':run_verified_candidates(inputs)}
            for name,entry in report.items():
                cases=entry.get('cases',[]);entry['executed_cases']=len(cases)
                fields=set().union(*(set(c) for c in cases)) if cases else set()
                entry['maxima_over_all_cases']={k:max(c[k] for c in cases if k in c) for k in fields if all(k not in c or isinstance(c[k],int) and not isinstance(c[k],bool) for c in cases)}
                entry['cases']=cases[:10 if name=='prefix_copy' else 4]
                entry['additional_cases_passed']=len(cases)-len(entry['cases'])
            print(json.dumps(report),flush=True)
        except Exception as exc:print(json.dumps({'audit_error':type(exc).__name__,'message':str(exc)}),flush=True);raise
