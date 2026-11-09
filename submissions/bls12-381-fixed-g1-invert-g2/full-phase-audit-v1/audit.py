"""Run ONLY in the shared verifier's credential-free isolated runner."""
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import random
import re
import sys

ROOT = Path(__file__).resolve().parent / 'files'
MASK = (1 << 64) - 1
R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
M = (1 << 61) - 1
Q = (1 << 127) - 1


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
    lines = [s.split('//', 1)[0].strip() for s in (ROOT / name).read_text().splitlines()]
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


def execute(ops, regs, memory=None, readable=None):
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

    for op, a in ops:
        if op == 'ret':
            return
        if op in ('ldr', 'ldp', 'str', 'stp'):
            assert memory is not None
            n = 2 if op in ('ldp', 'stp') else 1
            address = ','.join(a[n:]).replace(' ', '')
            match = re.fullmatch(r'\[(x[0-9]+|sp)(?:,#(.+))?\]', address)
            assert match, address
            base, offset = match.groups()
            at = read(base) + (integer(offset) if offset is not None else 0)
            if op in ('ldr', 'ldp'):
                for j in range(n):
                    assert readable is None or at + 8*j in readable
                    regs[index(a[j])] = memory[at + 8*j]
            else:
                for j in range(n):
                    memory[at + 8*j] = read(a[j])
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
        elif op in ('subs', 'sbcs'):
            value = read(a[1]) - read(a[2]) - (1-carry if op == 'sbcs' else 0)
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


def fingerprint_cases(rng, count):
    ops = parse('mersenne_fingerprint.S')
    assert len(ops) == 515
    tests = [(bytes(192), [0]*81), (bytes([255])*192, [M-1]*81)]
    for bit in (0, 59, 60, 63, 64, 121, 127, 1535):
        tests.append(((1 << bit).to_bytes(192, 'little'), [M-1 if i & 1 else 1 for i in range(81)]))
    for _ in range(count):
        tests.append((rng.getrandbits(1536).to_bytes(192, 'little'), [rng.randrange(M) for _ in range(81)]))
    runs = 0
    for key, coefficients in tests:
        raw = int.from_bytes(key, 'little')
        chunks = [(raw >> (60*i)) & ((1 << (36 if i == 25 else 60))-1) for i in range(26)]
        lanes = [(coefficients[27*j] + sum(coefficients[27*j+i+1]*chunks[i] for i in range(26))) % M for j in range(3)]
        expected = ((lanes[0] & ((1 << 54)-1)) | ((lanes[1] & ((1 << 54)-1)) << 54)
                    | ((lanes[2] & ((1 << 52)-1)) << 108)).to_bytes(32, 'little')
        for out in (0x300000, 0x100000, 0x200000):
            memory = {0x100000 + 8*i: int.from_bytes(key[8*i:8*i+8], 'little') for i in range(24)}
            memory.update({0x200000 + 8*i: v for i, v in enumerate(coefficients)})
            readable = set(memory)
            before = dict(memory)
            initial = [rng.getrandbits(64) for _ in range(32)]
            initial[:3] = [out, 0x100000, 0x200000]
            regs = list(initial)
            execute(ops, regs, memory, readable)
            actual = b''.join(memory[out+8*i].to_bytes(8, 'little') for i in range(4))
            assert actual == expected
            assert regs[18:] == initial[18:]
            changed = {out + 8*i for i in range(4)}
            assert set(memory) <= set(before) | changed
            assert all(memory[k] == v for k, v in before.items() if k not in changed)
            runs += 1
    return {'cases_with_aliases': runs, 'expanded_instructions': len(ops), 'inclusive_cap': 640}


def auxiliary_cases(rng, count):
    answers = {}
    values = [0, 1, 2, Q-1, Q-2, (1 << 64)-1, 1 << 64, 1 << 126]
    tests = list(product(values, values)) + [(rng.randrange(Q), rng.randrange(Q)) for _ in range(count)]
    for filename, opcode_count, kind, preserved in [('aux127_mul.S', 40, 'mul', 12), ('aux127_add.S', 13, 'add', 9)]:
        ops = parse(filename)
        assert len(ops) == opcode_count
        for a, b in tests:
            initial = [rng.getrandbits(64) for _ in range(32)]
            initial[:4] = [a & MASK, a >> 64, b & MASK, b >> 64]
            regs = list(initial)
            execute(ops, regs)
            actual = regs[0] + (regs[1] << 64)
            assert actual == ((a*b) if kind == 'mul' else (a+b)) % Q
            assert actual < Q
            assert regs[preserved:] == initial[preserved:]
        answers[kind] = {'cases': len(tests), 'expanded_instructions': len(ops)}
    return answers


def support_patterns(phases):
    def support(labels, edge, phase):
        s, t = edge
        row = {labels[s]: (-phase) % R}
        row[labels[t]] = (row.get(labels[t], 0) + 1) % R
        return {k for k, v in row.items() if v}
    four = three = 0
    for labels in product(range(4), repeat=4):
        allowed = ((labels[0] == labels[2] and labels[1] == labels[3])
                   or (labels[0] == labels[3] and labels[1] == labels[2])
                   or (labels[0] == labels[1] and labels[2] == labels[3]))
        for x, y in product(phases, repeat=2):
            if support(labels, (0, 1), x) == support(labels, (2, 3), y):
                assert allowed
            four += 1
    shapes = [((0, 1), (0, 2), (1, 2)), ((0, 1), (1, 2), (0, 2)), ((0, 2), (1, 2), (0, 1))]
    for labels in product(range(3), repeat=3):
        for x, y in product(phases, repeat=2):
            for e, f, outer in shapes:
                if support(labels, e, x) == support(labels, f, y):
                    assert labels[outer[0]] == labels[outer[1]]
                three += 1
    return {'four_time_phase_patterns': four, 'three_time_phase_patterns': three,
            'scope': 'necessary one-family support conditions, including empty supports'}


def relation_cases(rng, count, phases):
    tested = 0
    for phase in phases:
        for _ in range(count):
            a = rng.randrange(1, R)
            As, Bs, Bt = (rng.randrange(R) for _ in range(3))
            xs = (As + a*Bs) % R
            At = (phase*As + a*(phase*Bs-Bt)) % R
            det = (At*Bs-As*Bt) % R
            assert det == xs*(phase*Bs-Bt) % R
            if xs and det:
                denominator = (Bt-phase*Bs) % R
                assert denominator
                recovered = (phase*As-At)*pow(denominator, -1, R) % R
                assert recovered == a
            tested += 1
    return {'planted_scalar_relations': tested, 'scope': 'full Fr algebra, not GT inversion or pairing execution'}


def certificate():
    for exponent, modulus in [(61, M), (127, Q)]:
        residue = 4
        for _ in range(exponent-2):
            residue = (residue*residue-2) % modulus
        assert residue == 0
    assert 26*(M-1)*((1 << 60)-1)+(M-1) < 1 << 126
    W = 1 << 92
    root = math.isqrt(R // 12)
    while 12*root*root < R:
        root += 1
    assert 12*(root-1)**2 < R <= 12*root*root
    L = (1031*root+999)//1000
    assert L == 68152580234859358892647710588466745853
    b2, b3, b4 = (math.comb(L, k) for k in (2, 3, 4))
    amin, amax = F(Q-W+1, Q), F(Q+1, Q)
    delta = F(1 << 23, M**3)
    mu = F(b2*(W*W-1)*12*(R-1)**2, W*W*R**3)
    u2 = 3*144*(b3*(F(1, R*R)+F(1, R*W*W)) + b4*(F(1, R*R)+F(9, R*W**4)))
    caploss = F(b2, L//(1 << 17)+1)*(amax**4*(F(1, W*W)+F(12, R))+delta)
    success = amin**4*mu-amax**8*u2-caploss-F(L, R)-F(1, 1 << 970)
    assert success > F(39024, 100000)
    N = 25728*L+(1 << 120)
    assert N**20 < 20000**20*(1 << 2521)
    assert 6*(1 << 20)+1024*20000+65536 < 1 << 25
    assert 1024*(16*255+8*127+81*61) < 1 << 24
    return {'exact_success_test': True, 'success_lower_diagnostic': float(success),
            'exact_work_test': True, 'score_cap': '2521/20', 'sample_count': str(L),
            'scope': 'conditional on bounded fair-coin and finite-bank admissibility'}


request = json.loads(sys.stdin.read())
count, seed = request['count'], request['seed']
assert isinstance(count, int) and 1 <= count <= 64
assert isinstance(seed, int) and 0 <= seed <= 4294967295
rng = random.Random(seed)
phases = [pow(P, i, R) for i in range(12)]
assert len(set(phases)) == 12 and pow(P, 12, R) == 1
print(json.dumps({'status': 'passed', 'native_fingerprint': fingerprint_cases(rng, count),
                  'native_q127': auxiliary_cases(rng, count), 'support': support_patterns(phases),
                  'relations': relation_cases(rng, count, phases), 'certificate': certificate(),
                  'scope': 'bounded native interpretation and finite combinatorial/rational checks; no production attack, coin or storage certification'}, sort_keys=True))
