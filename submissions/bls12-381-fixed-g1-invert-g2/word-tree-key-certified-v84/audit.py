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
        self.table=2**43;self.block=96*(11*2**32+2**29)
        self.reg[31]=self.sp
        self.reg[18]=0xfeedfacedeadbeef
        for i in range(19,31):self.reg[i]=0xabc000+i
        self.saved=self.reg[18:]
    def readmem(self,at):
        assert at%8==0
        if hasattr(self,'fixed_s') and 56*2**40<=at<69.5*2**40:
            off=at-56*2**40;start=0
            for width in [35]*7+[34]*4:
                size=48*(1<<width)
                if off<size:
                    digit,part=divmod(off,48);value=self.fixed_s*RI*pow(2,start,P)*digit%P
                    return (value>>(8*part))&MASK
                off-=size;start+=width
            raise AssertionError('fixed lookup range')
        if self.table<=at<self.table+8*self.block:
            col,off=divmod(at-self.table,self.block);digit,off=divmod(off,96*2**32);index,part=divmod(off,96);row,limb=divmod(part,48)
            assert digit<12 and index<(2**29 if digit==11 else 2**32)
            value=self.coeffs[col][row]*RI*pow(2,32*digit,P)*index%P
            return (value>>(8*limb))&MASK
        assert at in self.mem,('uninitialized',hex(at))
        return self.mem[at]
    def writemem(self,at,value):
        assert at%8==0 and ((self.sp-816<=at<self.sp) or (self.out<=at<self.out+192)),('write',hex(at))
        self.mem[at]=value&MASK
    def field(self,at,value):
        for i in range(6):self.mem[at+8*i]=(value>>(64*i))&MASK
    def loadfield(self,at):return sum(self.readmem(at+8*i)<<(64*i) for i in range(6))
    def execute(self,name):
        ops=OPS[name];carry=0;zero=0;before=self.count
        def ix(s):return 31 if s=='sp' else int(s[1:])
        def rd(s):
            if s.startswith(':lo12:'):return RODATA[s.removeprefix(':lo12:')]&4095
            value=0 if s=='xzr' else integer(s) if s.startswith('#') else self.reg[ix(s)]
            return value&0xffffffff if s.startswith('w') else value
        def wr(s,v):
            if s!='xzr':self.reg[ix(s)]=v&MASK
        def shifted(a):
            value=rd(a[2])
            if len(a)==4:
                kind,n=a[3].split();n=integer(n);assert kind in ('lsl','lsr') and 0<=n<64
                value=(value<<n)&MASK if kind=='lsl' else value>>n
            return value
        labels={a[0]:i for i,(op,a) in enumerate(ops) if op=='@label'}
        pc=0
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
            if op in ('ldr','ldp','str','stp'):
                n=2 if op in ('ldp','stp') else 1
                addr=','.join(a[n:]).replace(' ','');m=re.fullmatch(r'\[(x\d+|sp)(?:,#(.+?))?\](!)?(?:,#(.+))?',addr);assert m,addr
                base,off,pre,post=m.groups();assert not(pre and post)
                displacement=integer(post if post else off if off else '0')
                if op in ('ldp','stp'):
                    assert displacement%8==0 and -512<=displacement<=504,('pair immediate',op,displacement)
                elif pre or post:
                    assert -256<=displacement<=255,('unscaled immediate',op,displacement)
                else:
                    assert displacement%8==0 and 0<=displacement<=32760,('scaled immediate',op,displacement)
                old=rd(base);at=old+(integer(off) if off else 0)
                if pre:wr(base,at)
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

OPS={n:parse(n+'.S') for n in ['recurrence_fp2_step','recurrence_fp2_map','fp381_add','fp381_sub','nh_fingerprint','fixed_s_lookup','canonical_twelve_bfirst','odd_trace_compact']}
OPS['canonical_conjugate_b1']=parse('canonical_b1.S')

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
        got=tuple(machine.loadfield(machine.out+48*j)*RI%P for j in range(4));assert got==oracle(c,x,y,z),(case,got)
        assert all(machine.loadfield(machine.out+48*j)<P for j in range(4))
        assert machine.reg[18:]==machine.saved
        assert machine.count==2711 and machine.calls['recurrence_fp2_map']==[2461]
        assert all(v==31 for n in ('fp381_add','fp381_sub') for v in machine.calls.get(n,[]))
        results.append({'case':case,'step_instructions':machine.count,'map_instructions':2461,'exact_third_alias':True,'abi_and_initialized_reads':True})
    return {'scope':'full-p native fused recurrence map/step instruction interpretation against independent dense Fp4 oracle; virtual read-only table test fixture only; no production fill, GT coverage, hardware or score qualification','cases':results}


def run_nh(inputs):
    rng=random.Random(inputs['seed']+72);q=2**61-1;cases=[]
    for case in range(inputs['count']):
        words=[rng.getrandbits(64) for _ in range(24)]
        if case==0:words=[0]*24
        if case==1:words=[MASK]*24
        keys=[[rng.getrandbits(64) for _ in range(24)] for _ in range(3)]
        affine=[[rng.randrange(q) for _ in range(4)] for _ in range(3)]
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
            ys.append(y&((1<<(52 if row==2 else 54))-1))
        want=ys[0]|(ys[1]<<54)|(ys[2]<<108)
        got=sum(machine.readmem(machine.out+8*j)<<(64*j) for j in range(4));assert got==want
        assert machine.count==375 and machine.reg[18:]==machine.saved
        assert all(machine.out<=at<machine.out+32 or before.get(at)==v for at,v in machine.mem.items())
        cases.append({'case':case,'instructions':375,'exact_input_alias':base==machine.out,'abi_and_initialized_reads':True})
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
        machine.reg[:4]=[machine.out,base,56*2**40,machine.const];before=dict(machine.mem)
        machine.execute('fixed_s_lookup');got=machine.loadfield(machine.out)
        assert got==x*scalar*RI%P and got<P
        assert machine.count==212 and machine.reg[18:]==machine.saved
        assert all(machine.out<=at<machine.out+48 or machine.sp-96<=at<machine.sp or before.get(at)==v for at,v in machine.mem.items())
        cases.append({'case':case,'instructions':212,'exact_input_alias':base==machine.out,'abi_and_initialized_reads':True})
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
        return tuple(-x%P for x in B) if first>=(P+1)//2 else B
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
        machine.reg[:5]=[machine.out,machine.prev1,56*2**40,machine.mod,machine.const]
        machine.execute('canonical_twelve_bfirst')
        assert machine.reg[18:]==machine.saved and machine.count<=813
        assert [machine.loadfield(machine.prev1+48*j) for j in range(4)]==phi
        assert all(v<=65 for v in machine.calls['canonical_conjugate_b1'])
        assert machine.calls['odd_trace_compact']==[572]
        assert machine.calls['fixed_s_lookup']==[212,212]
        if case==0:
            assert machine.reg[0]==1 and not any(machine.out<=at<machine.out+192 for at in machine.mem)
        else:
            even=sign(tuple(phi[2:]));Br,Bi=phi[2:];odd=sign((C[0]*(Br+Bi)%P,C[0]*(Br-Bi)%P))
            assert even!=odd
            expected=(tuple(phi[:2])+even) if words(even)<words(odd) else ((phi[0],-phi[1]%P)+odd)
            got=tuple(machine.loadfield(machine.out+48*j) for j in range(4))
            assert machine.reg[0]==0 and got==expected
        cases.append({'case':case,'instructions':machine.count,'key_bound':813,'immutable_raw_phi':True,'abi_and_initialized_reads':True,'identity_case':case==0})
    return {'scope':'actual complete key/sign/odd/lookup scalar instructions against independent full-p field/canonical-order formula; nonzero-b1 field cases plus true GT identity only; not empirical GT twelve-fiber theorem, production inversion, probability, hardware or qualification','cases':cases}

for line in sys.stdin:
    if line.strip():
        try:
            inputs=json.loads(line)
            print(json.dumps({'recurrence':run(inputs),'fingerprint':run_nh(inputs),'constant_lookup':run_fixed(inputs),'orbit_key':run_key(inputs)}),flush=True)
        except Exception as exc:print(json.dumps({'audit_error':type(exc).__name__,'message':str(exc)}),flush=True);raise
