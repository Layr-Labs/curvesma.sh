"""ONLY the shared verifier credential-free Docker runner may execute this file.
Bounded interpretation of actual full-width scalar-core/interpolation instructions.
No native hardware run, production inversion, hash/storage certification is claimed.
"""
import json
from pathlib import Path
import random
import re
import sys
ROOT = Path(__file__).resolve().parent / 'files'
MASK = (1 << 64)-1
P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab',16)
N0 = 0x89f3fffcfffcfffd

def integer(s):
    s = s.strip().removeprefix('#').strip()
    if s.startswith('(') and s.endswith(')'):
        s = s[1:-1]
    if '+' in s:
        parts = s.split('+')
        assert all(re.fullmatch(r'[0-9]+', t.strip()) for t in parts)
        return sum(int(t) for t in parts)
    assert re.fullmatch(r'-?(?:0x[0-9a-fA-F]+|[0-9]+)', s), s
    return int(s, 0)


def parse(name):
    lines = [s.split('//', 1)[0].strip().replace('12*__SIZEOF_POINTER__','96') for s in (ROOT / name).read_text().splitlines() if not s.lstrip().startswith('#')]
    macros, ordinary = {}, []
    pos = 0
    while pos < len(lines):
        line = lines[pos]
        pos += 1
        if line.startswith('.macro '):
            head = line.split(None, 2)
            args = [s.strip() for s in head[2].split(',')]
            body = []
            while lines[pos] != '.endm':
                body.append(lines[pos])
                pos += 1
            pos += 1
            macros[head[1]] = (args, body)
        else:
            ordinary.append(line)
    out = []

    def emit(rows, depth=0):
        assert depth <= 2
        conditions = []
        for line in rows:
            if line.startswith('.if '):
                v = integer(line[4:])
                assert v in (0, 1)
                conditions.append(bool(v))
                continue
            if line == '.else':
                assert conditions
                conditions[-1] = not conditions[-1]
                continue
            if line == '.endif':
                assert conditions
                conditions.pop()
                continue
            if not all(conditions):
                continue
            if not line or line.startswith('.') or line.endswith(':'):
                continue
            parts = line.split(None, 1)
            if parts[0] in macros:
                formals, template = macros[parts[0]]
                actual = [s.strip() for s in parts[1].split(',')]
                assert len(actual) == len(formals)
                expanded = []
                for row in template:
                    for k, v in sorted(zip(formals, actual), key=lambda kv: -len(kv[0])):
                        row = row.replace('\\' + k, v)
                    expanded.append(row)
                emit(expanded, depth + 1)
            else:
                out.append((parts[0], [] if len(parts) == 1 else [s.strip() for s in parts[1].split(',')]))
        assert not conditions

    emit(ordinary)
    return out


def execute(ops, regs, memory=None, readable=None, callbacks=None):
    carry = 0

    def index(s):
        return 31 if s == 'sp' else int(s[1:])

    def read(s):
        if s == 'xzr':
            return 0
        if s.startswith('#'):
            return integer(s)
        return regs[index(s)]

    def shifted(args):
        v = read(args[2])
        if len(args) == 4:
            kind, amount = args[3].split()
            n = integer(amount)
            assert 0 <= n < 64
            assert kind in ('lsl', 'lsr')
            v = ((v << n) & MASK) if kind == 'lsl' else v >> n
        else:
            assert len(args) == 3
        return v

    for pc, (op, a) in enumerate(ops):
        if op == 'ret':
            return
        if op in ('ldr', 'ldp', 'str', 'stp'):
            assert memory is not None
            n = 2 if op in ('ldp', 'stp') else 1
            address = ','.join(a[n:]).replace(' ', '')
            match = re.fullmatch(r'\[(x[0-9]+|sp)(?:,#(.+?))?\](!)?(?:,#(.+))?', address)
            assert match, address
            base, offset, pre, post = match.groups()
            assert not (pre and post)
            before_base = read(base)
            at = before_base + (integer(offset) if offset is not None else 0)
            if pre:
                regs[index(base)] = at & MASK
            if op in ('ldr', 'ldp'):
                for j in range(n):
                    assert readable is None or at + 8*j in readable
                    regs[index(a[j])] = memory[at + 8*j]
            else:
                for j in range(n):
                    memory[at + 8*j] = read(a[j])
                    if readable is not None:
                        readable.add(at + 8*j)
            if post:
                regs[index(base)] = (before_base + integer(post)) & MASK
            continue
        if op == 'bl':
            assert callbacks is not None and a[0] in callbacks
            regs[30] = (pc+1)*4
            callbacks[a[0]](regs, memory, readable)
            carry = 0
            continue
        if op == 'mov':
            value = read(a[1])
        elif op == 'movn':
            value = ~read(a[1])
        elif op == 'and':
            value = read(a[1]) & read(a[2])
        elif op == 'orr':
            value = read(a[1]) | shifted(a)
        elif op == 'lsr':
            value = read(a[1]) >> read(a[2])
        elif op == 'extr':
            value = ((read(a[1]) << 64) | read(a[2])) >> read(a[3])
        elif op in ('mul', 'umulh'):
            value = read(a[1]) * read(a[2])
            if op == 'umulh':
                value >>= 64
        elif op in ('add', 'adds', 'adc', 'adcs'):
            value = read(a[1]) + shifted(a) + (carry if op in ('adc', 'adcs') else 0)
            if op in ('adds', 'adcs'):
                carry = int(value > MASK)
        elif op in ('sub', 'subs', 'sbcs'):
            value = read(a[1]) - shifted(a) - (1-carry if op == 'sbcs' else 0)
            if op != 'sub':
                carry = int(value >= 0)
        elif op == 'csel':
            assert a[3] in ('hs', 'lo')
            choose = bool(carry) if a[3] == 'hs' else not carry
            value = read(a[1] if choose else a[2])
        else:
            raise AssertionError('unsupported opcode ' + op)
        if a[0] != 'xzr':
            regs[index(a[0])] = value & MASK
    raise AssertionError('missing RET')



def store_field(memory, at, value):
    for i in range(6):
        memory[at+8*i] = (value >> (64*i)) & MASK

def load_field(memory, at):
    return sum(memory[at+8*i] << (64*i) for i in range(6))


POINTS = [0,1,-1,2,-2,3,None]
DENOMINATORS = [-12,12,24,-24,-120,120,1]
SHIFT_SCALES = [1,2,2,1,1,2,2]

def evaluations(raw, weighted):
    # Independent high-level cold representation from native GT coefficients.
    pairs = [(raw[0],raw[4]),(raw[3],raw[2]),(raw[1],raw[5])]
    values = []
    for group,((ar,ai),(br,bi)) in enumerate(pairs):
        coefficients = [(ar-ai)%P,(br-bi)%P,ai,bi]
        for i,x in enumerate(POINTS):
            if x is None:
                value = coefficients[3]
            else:
                value = 0
                for coefficient in reversed(coefficients):
                    value = (value*x+coefficient)%P
            if weighted:
                scale = 1 if group == 0 else SHIFT_SCALES[i]
                value = value*scale*pow(DENOMINATORS[i],-1,P)%P
            values.append(value)
    assert len(values) == 21
    return values

def dense(x,y):
    # Full independent schoolbook multiplication, w^6=1+u.
    order = [0,3,1,4,2,5]
    x,y = [x[i] for i in order],[y[i] for i in order]
    z = [(0,0) for _ in range(11)]
    for i,(ar,ai) in enumerate(x):
        for j,(br,bi) in enumerate(y):
            cr,ci = z[i+j]
            z[i+j] = ((cr+ar*br-ai*bi)%P,(ci+ar*bi+ai*br)%P)
    for i in range(10,5,-1):
        ar,ai = z[i]
        br,bi = z[i-6]
        z[i-6] = ((br+ar-ai)%P,(bi+ar+ai)%P)
    rinv = pow(1 << 384,-1,P)
    return tuple(value*rinv%P for pair in [z[0],z[3]] for value in pair)

def cases(rng,count):
    ops = {n:parse(n+'.S') for n in ['project_toom','pinned_mont_core','fp381_add','fp381_sub']}
    assert [len(ops[n]) for n in ops] == [905,341,31,31]
    runs = 0
    for case in range(count+2):
        x = [(rng.randrange(P),rng.randrange(P)) for _ in range(6)]
        y = [(rng.randrange(P),rng.randrange(P)) for _ in range(6)]
        if case == 0:
            x,y = [(0,0)]*6,[(0,0)]*6
        if case == 1:
            x,y = [(P-1,P-1)]*6,[(P-1,0)]*6
        expected = dense(x,y)
        for out in [0x700000,0x100000,0x200000]:
            memory = {}
            for base,data in [(0x100000,evaluations(x,False)),(0x200000,evaluations(y,True))]:
                for i,value in enumerate(data):
                    store_field(memory,base+48*i,value)
            store_field(memory,0x500000,P)
            readable = set(memory)
            before = dict(memory)
            regs = [rng.getrandbits(64) for _ in range(32)]
            regs[:5] = [out,0x100000,0x200000,0x500000,N0]
            regs[31] = 0x900000
            initial = list(regs)
            steps = [0]
            callbacks = {}
            def run(name,r,m,rd):
                steps[0] += len(ops[name])
                execute(ops[name],r,m,rd,callbacks)
            callbacks.update({
                '__mul_mont_384':lambda r,m,rd:run('pinned_mont_core',r,m,rd),
                'fp381_add':lambda r,m,rd:run('fp381_add',r,m,rd),
                'fp381_sub':lambda r,m,rd:run('fp381_sub',r,m,rd)})
            run('project_toom',regs,memory,readable)
            observed = tuple(load_field(memory,out+d) for d in [0,48,384,432])
            assert observed == expected
            assert steps[0] == 11631
            assert regs[18] == initial[18] and regs[19:32] == initial[19:32]
            allowed = set(range(out,out+96,8)) | set(range(out+384,out+480,8))
            allowed |= set(range(initial[31]-4128,initial[31],8))
            assert all(at in allowed or memory[at] == value for at,value in before.items())
            assert all(at in before or at in allowed for at in memory)
            runs += 1
    return {'bounded_interpreted_cases':runs,'complete_leaf_instructions':11631,
            'inclusive_bound':11776,'scope':'actual native instruction semantics vs independent dense convolution;full-width inputs,ABI,aliases,stack bounds'}

request = json.loads(sys.stdin.read())
count,seed = request['count'],request['seed']
assert isinstance(count,int) and 1 <= count <= 32
assert isinstance(seed,int) and 0 <= seed <= 4294967295
print(json.dumps({'status':'passed','projection':cases(random.Random(seed),count),
                  'scope':'no native hardware,production inversion,hash or bank certification'},sort_keys=True))
