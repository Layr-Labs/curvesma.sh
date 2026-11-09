"""Run only in the shared verifier's credential-free isolated runner."""
import itertools
import json
import math
from pathlib import Path
import random
import re
import sys
from fractions import Fraction

MASK = (1 << 64) - 1
ROOT = Path(__file__).resolve().parent / 'files'

def execute(ops, regs, memory=None):
    carry = 0
    def index(x):
        return 31 if x == 'sp' else int(x[1:])
    def read(x):
        if x == 'xzr':
            return 0
        return int(x[1:], 0) if x.startswith('#') else regs[index(x)]
    for op, z in ops:
        if op == 'ret':
            return regs
        if op in ('ldr','ldp','str','stp'):
            assert memory is not None
            n = 2 if op in ('ldp','stp') else 1
            address_text = ','.join(z[n:]).replace(' ','')
            match = re.fullmatch(r'\[(sp|x[0-9]+)(?:,#(-?[0-9]+))?\](!)?(?:,#(-?[0-9]+))?',address_text)
            assert match, address_text
            base, offset, pre, post = match.groups()
            assert not (pre and post)
            address = read(base) + int(offset or 0)
            if pre:
                regs[index(base)] = address & MASK
            if op in ('ldr','ldp'):
                for j in range(n):
                    regs[index(z[j])] = memory[address + 8*j]
            else:
                for j in range(n):
                    memory[address + 8*j] = read(z[j])
            if post:
                regs[index(base)] = (read(base) + int(post)) & MASK
            continue
        if op == 'mov':
            value = read(z[1])
        elif op == 'movn':
            value = ~read(z[1])
        elif op == 'and':
            value = read(z[1]) & read(z[2])
        elif op == 'orr':
            value = read(z[1]) | read(z[2])
        elif op == 'lsr':
            value = read(z[1]) >> read(z[2])
        elif op == 'extr':
            value = ((read(z[1]) << 64) | read(z[2])) >> read(z[3])
        elif op in ('mul','umulh'):
            value = read(z[1]) * read(z[2])
            if op == 'umulh':
                value >>= 64
        elif op in ('add','adds','adcs','adc'):
            value = read(z[1]) + read(z[2]) + (carry if op in ('adcs','adc') else 0)
            if op in ('adds','adcs'):
                carry = int(value > MASK)
        elif op in ('subs','sbcs'):
            value = read(z[1]) - read(z[2]) - (1-carry if op == 'sbcs' else 0)
            carry = int(value >= 0)
        elif op == 'csel':
            assert z[3] == 'hs'
            value = read(z[1]) if carry else read(z[2])
        else:
            raise AssertionError('Unsupported opcode: ' + op)
        assert z[0] == 'sp' or (z[0].startswith('x') and z[0] != 'xzr')
        regs[index(z[0])] = value & MASK
    raise AssertionError('Missing RET')


def cios_program():
    lines = [x.split('//',1)[0].strip() for x in (ROOT/'cios381_mul.S').read_text().splitlines()]
    template = None
    output = []
    def emit(line):
        if not line or line.startswith('.') or line.endswith(':'):
            return
        parts = line.split(None,1)
        if parts[0] == 'CIOS_ACCUM':
            assert template is not None
            args = [x.strip() for x in parts[1].split(',')]
            assert len(args) == len(formals) == 3
            for row in template:
                for formal, arg in zip(formals,args):
                    row = row.replace('\\'+formal,arg)
                emit(row)
        else:
            output.append((parts[0],[] if len(parts)==1 else [x.strip() for x in parts[1].split(',')]))
    pos = 0
    while pos < len(lines):
        line = lines[pos]
        if line.startswith('.macro '):
            parts = line.split(None,2)
            assert parts[1] == 'CIOS_ACCUM'
            formals = [x.strip() for x in parts[2].split(',')]
            template = []
            pos += 1
            while lines[pos] != '.endm':
                template.append(lines[pos]); pos += 1
        elif line.startswith('.irp '):
            data = [x.strip() for x in line.split(None,1)[1].split(',')]
            variable = data[0]
            values = [int(x) for x in data[1:]]
            assert values == list(range(6))
            body = []
            pos += 1
            while lines[pos] != '.endr':
                body.append(lines[pos]); pos += 1
            for value in values:
                for row in body:
                    row = row.replace('\\'+variable,str(value))
                    row = row.replace('#(8*'+str(value)+')','#'+str(8*value))
                    emit(row)
        else:
            emit(line)
        pos += 1
    assert len(output) == 551
    assert output[-1][0] == 'ret'
    return output


def cios_cases(rng,count):
    p = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab',16)
    n0 = 0x89f3fffcfffcfffd
    inverse_R = pow(1<<384,-1,p)
    ops = cios_program()
    edges = [0,1,p-1,p-2,(1<<64)-1,1<<64,(1<<380)-1,1<<380]
    tests = list(itertools.product(edges,edges))
    tests.extend((rng.randrange(p),rng.randrange(p)) for _ in range(min(count,64)))
    for a,b in tests:
        initial = [rng.getrandbits(64) for _ in range(31)] + [0x800000]
        regs = initial.copy()
        regs[0:5] = [0x400000,0x100000,0x200000,0x300000,n0]
        memory = {}
        for base,value in [(0x100000,a),(0x200000,b),(0x300000,p)]:
            for j in range(6):
                memory[base+8*j] = (value>>(64*j)) & MASK
        execute(ops,regs,memory)
        actual = sum(memory[0x400000+8*j]<<(64*j) for j in range(6))
        assert actual == a*b*inverse_R % p
        assert actual < p
        assert regs[18:31] == initial[18:31]
        assert regs[31] == initial[31]
    return {'cases':len(tests),'expanded_instructions':len(ops),'native_scale':'R=2^384'}

def moment_certificate():
    r = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
    W = 1 << 92
    R = math.isqrt((r + 10) // 11)
    while 11 * R * R < r:
        R += 1
    L = (1031 * R + 999) // 1000
    b2 = math.comb(L, 2)
    b3 = math.comb(L, 3)
    b4 = math.comb(L, 4)
    mu = Fraction(b2 * (W * W - 1) * 11 * (r - 1), W * W * r * r)
    U = 3 * 121 * (b3 * (Fraction(1, r*r) + Fraction(1, r*W*W))
                   + b4 * (Fraction(1, r*r) + Fraction(4, r*W**4)))
    q = (1 << 127) - 1
    qh = (1 << 255) - 19
    amin = Fraction(q - W + 1, q)
    amax = Fraction(q + 1, q)
    delta = Fraction(1 << 95, qh)
    eq = b2 * (amax**4 * (Fraction(1, W*W) + Fraction(12, r)) + delta)
    probability = amin**4 * mu - amax**8 * U - eq / (L//16 + 1) - Fraction(L, r) - Fraction(1, 1 << 970)
    assert probability > Fraction(39024, 100000)
    N = 28640 * L + (1 << 120)
    assert N**40 < 20000**40 * (1 << 5051)
    return {'success_lower_diagnostic': float(probability), 'score_cap': '5051/40',
            'exact_success_test': True, 'exact_cost_test': True, 'sample_count': str(L)}

request = json.loads(sys.stdin.read())
count = request['count']
rng = random.Random(request['seed'])
q127 = (1 << 127) - 1
cios_result = cios_cases(rng,count)
ll = 4
for _ in range(125):
    ll = (ll * ll - 2) % q127
assert ll == 0
print(json.dumps({'status':'passed','scope':'bounded full-width native-kernel interpretation and exact public-parameter certificates; no production inversion',
                  'cios':cios_result,'lucas_lehmer_127':ll,
                  'certificate':moment_certificate()},sort_keys=True))
