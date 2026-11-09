"""ONLY the shared verifier credential-free Docker runner may execute this file.
Bounded interpretation of actual full-width scalar-core/fused-Fp2 instructions.
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
            value = read(a[1]) - read(a[2]) - (1-carry if op == 'sbcs' else 0)
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

def cases(rng,count):
    wrapper,core,subtract = parse('native_fp2_fused.S'),parse('pinned_mont_core.S'),parse('fp381_sub.S')
    assert (len(wrapper),len(core),len(subtract)) == (81,341,31)
    assert (N0*(P & MASK)+1) & MASK == 0
    rinv = pow(1 << 384,-1,P)
    values = [(0,0,0,0),(1,1,1,1),(P-1,P-1,P-1,P-1),
              (P-1,1,1,P-1),(0,P-1,P-1,0),(P//2,P//2,1,P-1)]
    values += [tuple(rng.randrange(P) for _ in range(4)) for _ in range(count)]
    runs = 0
    for ar,ai,br,bi in values:
        expected = ((ar*br-ai*bi)*rinv % P,(ar*bi+ai*br)*rinv % P)
        for out in (0x700000,0x100000,0x200000):
            memory = {}
            for at,value in [(0x100000,ar),(0x100030,ai),(0x200000,br),(0x200030,bi),
                             (0x300000,(ar+ai)%P),(0x400000,(br+bi)%P),(0x500000,P)]:
                store_field(memory,at,value)
            readable = set(memory)
            before = dict(memory)
            regs = [rng.getrandbits(64) for _ in range(32)]
            regs[:7] = [out,0x100000,0x200000,0x300000,0x400000,0x500000,N0]
            regs[31] = 0x900000
            initial = list(regs)
            steps = [0]
            def call_core(r,m,readable):
                steps[0] += len(core)
                execute(core,r,m,readable)
            def call_sub(r,m,readable):
                steps[0] += len(subtract)
                execute(subtract,r,m,readable)
            execute(wrapper,regs,memory,readable,{'__mul_mont_384':call_core,'fp381_sub':call_sub})
            steps[0] += len(wrapper)
            assert steps[0] == 1197
            assert (load_field(memory,out),load_field(memory,out+48)) == expected
            assert regs[18] == initial[18]
            assert regs[19:31] == initial[19:31]
            assert regs[31] == initial[31]
            allowed = set(range(out,out+96,8)) | set(range(initial[31]-336,initial[31],8))
            assert all(at in allowed or memory[at] == value for at,value in before.items())
            assert all(at in before or at in allowed for at in memory)
            runs += 1
    return {'bounded_interpreted_cases':runs,'complete_leaf_instructions':1197,
            'inclusive_outer_bound':1216,'scope':'actual instruction semantics and full-width arithmetic/ABI/aliases only'}

def projection_cases(rng,count):
    ops = {n:parse(n+'.S') for n in ['project_quartic','native_fp2_fused','pinned_mont_core','fp381_add','fp381_sub']}
    assert [len(ops[n]) for n in ops] == [353,81,341,31,31]
    rinv = pow(1 << 384,-1,P)
    runs = 0
    for case in range(count):
        x = [(rng.randrange(P),rng.randrange(P)) for _ in range(6)]
        y = [(rng.randrange(P),rng.randrange(P)) for _ in range(6)]
        if case == 0:
            x,y = [(0,0)]*6,[(0,0)]*6
        if case == 1:
            x,y = [(P-1,P-1)]*6,[(P-1,0)]*6
        # Independent dense convolution in increasing w degree, w^6=1+u.
        order = [0,3,1,4,2,5]
        xc,yc = [x[i] for i in order],[y[i] for i in order]
        z = [(0,0) for _ in range(11)]
        for i,(ar,ai) in enumerate(xc):
            for j,(br,bi) in enumerate(yc):
                cr,ci = z[i+j]
                z[i+j] = ((cr+ar*br-ai*bi)%P,(ci+ar*bi+ai*br)%P)
        for i in range(10,5,-1):
            ar,ai = z[i]
            br,bi = z[i-6]
            z[i-6] = ((br+ar-ai)%P,(bi+ar+ai)%P)
        expected = tuple(v*rinv % P for pair in [z[0],z[3]] for v in pair)
        memory = {}
        for base,data in [(0x100000,x),(0x200000,y)]:
            for i,(ar,ai) in enumerate(data):
                store_field(memory,base+96*i,ar)
                store_field(memory,base+96*i+48,ai)
                store_field(memory,base+640+48*i,(ar+ai)%P)
        store_field(memory,0x500000,P)
        readable = set(memory)
        before = dict(memory)
        regs = [rng.getrandbits(64) for _ in range(32)]
        regs[:5] = [0x700000,0x100000,0x200000,0x500000,N0]
        regs[31] = 0x900000
        initial = list(regs)
        steps = [0]
        callbacks = {}
        def run(name,r,m,rd):
            steps[0] += len(ops[name])
            execute(ops[name],r,m,rd,callbacks)
        callbacks.update({
            'native_fp2_fused':lambda r,m,rd:run('native_fp2_fused',r,m,rd),
            '__mul_mont_384':lambda r,m,rd:run('pinned_mont_core',r,m,rd),
            'fp381_add':lambda r,m,rd:run('fp381_add',r,m,rd),
            'fp381_sub':lambda r,m,rd:run('fp381_sub',r,m,rd)})
        run('project_quartic',regs,memory,readable)
        observed = tuple(load_field(memory,0x700000+d) for d in [0,48,384,432])
        assert observed == expected
        assert steps[0] == 12738
        assert regs[18] == initial[18] and regs[19:32] == initial[19:32]
        allowed = set(range(0x700000,0x700060,8)) | set(range(0x700180,0x7001e0,8))
        allowed |= set(range(initial[31]-2144,initial[31],8))
        assert all(at in allowed or memory[at] == value for at,value in before.items())
        assert all(at in before or at in allowed for at in memory)
        runs += 1
    return {'bounded_interpreted_cases':runs,'complete_leaf_instructions':12738,
            'inclusive_bound':13376,'scope':'actual full native instruction interpretation vs independent dense polynomial convolution'}

request = json.loads(sys.stdin.read())
count,seed = request['count'],request['seed']
assert isinstance(count,int) and 1 <= count <= 64
assert isinstance(seed,int) and 0 <= seed <= 4294967295
print(json.dumps({'status':'passed','fused_fp2':cases(random.Random(seed),count),
                  'projection':projection_cases(random.Random(seed+1),min(count,16)),
                  'scope':'no physical/native run,production inversion,hash or storage certification'},sort_keys=True))
