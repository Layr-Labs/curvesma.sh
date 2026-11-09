"""Run only in the shared verifier's credential-free isolated runner."""
import itertools
import json
import math
from pathlib import Path
import random
import sys
from fractions import Fraction

MASK = (1 << 64) - 1
ROOT = Path(__file__).resolve().parent / 'files'
EXPECTED = {'aux127_mul.S': 40, 'aux127_add.S': 13,
            'aux255_mul128.S': 87, 'aux255_add.S': 24}

def program(name):
    out = []
    for raw in (ROOT / name).read_text().splitlines():
        line = raw.split('//', 1)[0].strip()
        if not line or line.startswith('.') or line.endswith(':'):
            continue
        fields = line.split(None, 1)
        out.append((fields[0], [] if len(fields) == 1 else [x.strip() for x in fields[1].split(',')]))
    assert len(out) == EXPECTED[name]
    assert out[-1][0] == 'ret'
    return out

def execute(ops, regs):
    carry = 0
    def read(x):
        if x == 'xzr':
            return 0
        return int(x[1:], 0) if x.startswith('#') else regs[int(x[1:])]
    for op, z in ops:
        if op == 'ret':
            return regs
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
        elif op in ('mul', 'umulh'):
            value = read(z[1]) * read(z[2])
            if op == 'umulh':
                value >>= 64
        elif op in ('adds', 'adcs', 'adc'):
            value = read(z[1]) + read(z[2]) + (carry if op != 'adds' else 0)
            if op != 'adc':
                carry = int(value > MASK)
        elif op in ('subs', 'sbcs'):
            value = read(z[1]) - read(z[2]) - (1 - carry if op == 'sbcs' else 0)
            carry = int(value >= 0)
        elif op == 'csel':
            assert z[3] == 'hs'
            value = read(z[1]) if carry else read(z[2])
        else:
            raise AssertionError('Unsupported opcode: ' + op)
        assert z[0].startswith('x') and z[0] != 'xzr'
        regs[int(z[0][1:])] = value & MASK
    raise AssertionError('Missing RET')

def kernel_cases(name, bits_a, bits_b, modulus, words_a, words_b, words_out, preserved, rng, count):
    ops = program(name)
    vals_a = [0, 1, modulus - 1, modulus - 2, (1 << 63) - 1,
              1 << 63, 1 << 64, (1 << (bits_a - 1)) - 1]
    bound_b = (1 << bits_b) if name == 'aux255_mul128.S' else modulus
    vals_b = [0, 1, bound_b - 1, bound_b - 2, (1 << 63) - 1,
              1 << 63, 1 << 64, (1 << (bits_b - 1)) - 1]
    tests = list(itertools.product(vals_a, vals_b))
    tests.extend((rng.randrange(modulus), rng.randrange(bound_b)) for _ in range(count))
    for a, b in tests:
        initial = [rng.getrandbits(64) for _ in range(31)]
        regs = initial.copy()
        for j in range(words_a):
            regs[j] = (a >> (64 * j)) & MASK
        for j in range(words_b):
            regs[words_a + j] = (b >> (64 * j)) & MASK
        result = execute(ops, regs)
        actual = sum(result[j] << (64 * j) for j in range(words_out))
        expected = ((a * b) if 'mul' in name else (a + b)) % modulus
        assert actual == expected, (name, hex(a), hex(b), hex(actual), hex(expected))
        assert actual < modulus
        assert result[preserved:] == initial[preserved:]
    return len(tests)

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
    N = 31200 * L + (1 << 120)
    assert N**5 < 20000**5 * (1 << 632)
    return {'success_lower_diagnostic': float(probability), 'score_cap': '632/5',
            'exact_success_test': True, 'exact_cost_test': True, 'sample_count': str(L)}

def toy_independence():
    outputs = set()
    for coeff in itertools.product(range(7), repeat=4):
        outputs.add(tuple(sum(coeff[j] * pow(x, j, 7) for j in range(4)) % 7 for x in range(4)))
    assert len(outputs) == 7**4
    counts = {}
    for a, b, c in itertools.product(range(7), repeat=3):
        pair = (c, (a + 2*b + c) % 7)
        counts[pair] = counts.get(pair, 0) + 1
    assert len(counts) == 49 and set(counts.values()) == {7}
    return {'four_wise_toy_outputs': len(outputs), 'pairwise_toy_outputs': len(counts)}

request = json.loads(sys.stdin.read())
count = request['count']
rng = random.Random(request['seed'])
q127 = (1 << 127) - 1
q255 = (1 << 255) - 19
cases = {}
cases['aux127_mul.S'] = kernel_cases('aux127_mul.S',127,127,q127,2,2,2,12,rng,count)
cases['aux127_add.S'] = kernel_cases('aux127_add.S',127,127,q127,2,2,2,9,rng,count)
cases['aux255_mul128.S'] = kernel_cases('aux255_mul128.S',255,128,q255,4,2,4,18,rng,count)
cases['aux255_add.S'] = kernel_cases('aux255_add.S',255,255,q255,4,4,4,15,rng,count)
ll = 4
for _ in range(125):
    ll = (ll * ll - 2) % q127
assert ll == 0
print(json.dumps({'status':'passed','scope':'bounded interpreted arithmetic and exact public-parameter certificates; no production inversion',
                  'kernel_cases':cases,'instruction_counts':EXPECTED,'lucas_lehmer_127':ll,
                  'certificate':moment_certificate(),'toy':toy_independence()},sort_keys=True))
