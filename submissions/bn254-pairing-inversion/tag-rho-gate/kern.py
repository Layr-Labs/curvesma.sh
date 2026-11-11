"""Instruction interpreter, generated field kernels and reference tower arithmetic
(imported by experiment.py and toy.py). Python standard library only."""
import math
import random

M64 = (1 << 64) - 1
P = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
R = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
U = 0x44e992b44a6909f1
G_PINNED = [int(c, 16) for c in (
    "28c6e04df059260df7d2d2a1f9b5f77676d1939847852c4ed50d2318744c1d5f",
    "17bb74adab1705c26133af1dac87044a3833ac011018e8158da48382bbd2dcd6",
    "0d3bd72f54d742f78ea9e6015c8ea2f2e7fbb728c9c905ec531dcf7de5b246f0",
    "090cb8ee97e091a667af03882b06c3ecb4e437993cbd1b05b98c7f9dfcfe9c40",
    "016b6d855b5cbf76f9829a309db52f5c442f65ae29f996af59d65f85f4afe78a",
    "0a0272204db51dadc0342bd318b9302a44faec12ff500bdd4d4b012ffe45f36f",
    "084f330485b09e866bc2f2ea2b897394deaf3f12aa31f28cb0552990967d4704",
    "27ed208e7a0b55ae6e710bbfbd2fd922669c026360e37cc5b2ab862411536104",
    "2067586885c3318eeffa1938c754fe3c60224ee5ae15e66af6b5104c47c8c5d8",
    "279db296f9d479292532c7c493d8e0722b6efae42158387564889c79fc038ee3",
    "2b03614464f04dd772d86df88674c270ffc8747ea13e72da95e3594468f222c4",
    "108c19d15f9446f744d0f110405d3856d6cc3bda6c4d537663729f5257628417")]

# ---------------------------------------------------------------- interpreter

FORBIDDEN = {'x18'}


class Program:
    def __init__(self):
        self.code, self.labels = [], {}

    def label(self, name):
        self.labels[name] = len(self.code)

    def emit(self, op, *a):
        for x in a:
            assert x not in FORBIDDEN, x
        if op == 'add':
            assert 0 <= a[2] <= 4095, a                      # ADD (immediate), 12 bits
        if op in ('ldp', 'stp', 'stp_pre', 'ldp_post'):
            assert -512 <= a[3] <= 504 and a[3] % 8 == 0, a  # LDP/STP signed imm7 * 8
        if op in ('ldr', 'str'):
            assert 0 <= a[2] <= 32760 and a[2] % 8 == 0, a   # LDR unsigned imm12 * 8
        if op in ('sub_sp', 'add_sp'):
            assert 0 <= a[0] <= 4095 << 12 and a[0] % 4096 == 0, a  # SUB/ADD sp, #imm12, LSL #12
        if op in ('lsr', 'lsl', 'add_lsl', 'extr'):
            assert 0 <= a[-1] <= 63, a                       # shift amount, EXTR lsb
        if op == 'subi':
            assert 0 <= a[2] <= 4095, a                      # SUBS (immediate), 12 bits
        if op in ('b', 'bc', 'cbz', 'bl'):
            assert isinstance(a[-1], str), a                 # branch and call targets are labels
        if op == 'movi':
            assert 0 <= a[1] < 65536, a                      # MOVZ
        self.code.append((op, a))


class Machine:
    """Executes the programs on 64-bit words; every instruction executed counts one."""

    def __init__(self, prog):
        self.p = prog
        self.r = {'x%d' % i: 0 for i in range(31)}
        self.r['sp'] = 1 << 40
        self.C = self.Z = 0
        self.mem = {}
        self.count = 0
        self.depth = 0

    def g(self, a):
        return 0 if a == 'xzr' else self.r[a]

    def s(self, a, v):
        if a != 'xzr':
            self.r[a] = v & M64

    def ld(self, ad):
        assert ad % 8 == 0
        return self.mem.get(ad, (ad * 0x9E3779B97F4A7C15 + 0x5851F42D4C957F2D) & M64)   # unwritten: a scramble

    def st(self, ad, v):
        assert ad % 8 == 0
        self.mem[ad] = v & M64

    def cond(self, c):
        return {'cs': self.C == 1, 'cc': self.C == 0, 'eq': self.Z == 1, 'ne': self.Z == 0}[c]

    def call(self, label):
        start, pc, stack = self.count, self.p.labels[label], []
        g, s = self.g, self.s
        while True:
            op, a = self.p.code[pc]
            self.count += 1
            pc += 1
            if op == 'ldr':
                s(a[0], self.ld(g(a[1]) + a[2]))
            elif op == 'ldp':
                b = g(a[2]) + a[3]; s(a[0], self.ld(b)); s(a[1], self.ld(b + 8))
            elif op == 'stp':
                b = g(a[2]) + a[3]; self.st(b, g(a[0])); self.st(b + 8, g(a[1]))
            elif op == 'stp_pre':
                self.r['sp'] += a[3]; b = self.r['sp']
                self.st(b, g(a[0])); self.st(b + 8, g(a[1]))
                self.depth = max(self.depth, (1 << 40) - self.r['sp'])
            elif op == 'ldp_post':
                b = self.r['sp']; s(a[0], self.ld(b)); s(a[1], self.ld(b + 8)); self.r['sp'] += a[3]
            elif op == 'sub_sp':
                self.r['sp'] -= a[0]; self.depth = max(self.depth, (1 << 40) - self.r['sp'])
            elif op == 'add_sp':
                self.r['sp'] += a[0]
            elif op == 'mul':
                s(a[0], g(a[1]) * g(a[2]))
            elif op == 'umulh':
                s(a[0], (g(a[1]) * g(a[2])) >> 64)
            elif op in ('adds', 'adcs', 'adc'):
                t = g(a[1]) + g(a[2]) + (self.C if op != 'adds' else 0)
                s(a[0], t)
                if op != 'adc':
                    self.C, self.Z = t >> 64, int(t & M64 == 0)
            elif op in ('subs', 'sbcs'):
                t = g(a[1]) - g(a[2]) - ((1 - self.C) if op == 'sbcs' else 0)
                s(a[0], t)
                self.C, self.Z = int(t >= 0), int(t & M64 == 0)
            elif op == 'csel':
                s(a[0], g(a[1]) if self.cond(a[3]) else g(a[2]))
            elif op == 'cset':
                s(a[0], 1 if self.cond(a[1]) else 0)
            elif op == 'orr':
                s(a[0], g(a[1]) | g(a[2]))
            elif op == 'mov':
                s(a[0], g(a[1]))
            elif op == 'movi':
                s(a[0], a[1])
            elif op == 'add':
                s(a[0], g(a[1]) + a[2])
            elif op == 'str':
                self.st(g(a[1]) + a[2], g(a[0]))
            elif op == 'and_r':
                s(a[0], g(a[1]) & g(a[2]))
            elif op == 'lsr':
                s(a[0], g(a[1]) >> a[2])
            elif op == 'lsl':
                s(a[0], g(a[1]) << a[2])
            elif op == 'add_r':
                s(a[0], g(a[1]) + g(a[2]))
            elif op == 'add_lsl':                # ADD (shifted register)
                s(a[0], g(a[1]) + (g(a[2]) << a[3]))
            elif op == 'extr':                   # EXTR: bits lsb.. of the pair (a[1] high, a[2] low)
                s(a[0], ((g(a[1]) << 64 | g(a[2])) >> a[3]))
            elif op == 'madd':
                s(a[0], g(a[1]) * g(a[2]) + g(a[3]))
            elif op == 'subi':                   # SUBS (immediate)
                t = g(a[1]) - a[2]; s(a[0], t); self.C, self.Z = int(t >= 0), int(t & M64 == 0)
            elif op == 'cmp':
                t = g(a[0]) - g(a[1]); self.C, self.Z = int(t >= 0), int(t & M64 == 0)
            elif op == 'b':                      # real branches: code addresses are instruction indices
                pc = self.p.labels[a[0]]
            elif op == 'bc':
                if self.cond(a[0]): pc = self.p.labels[a[1]]
            elif op == 'cbz':
                if g(a[0]) == 0: pc = self.p.labels[a[1]]
            elif op == 'bl':
                stack.append(pc); self.r['x30'] = pc; pc = self.p.labels[a[0]]   # BL writes the link register
            elif op == 'ret':                    # RET branches to x30; check it against the call stack
                if not stack:
                    return self.count - start
                assert self.r['x30'] == stack[-1], ('RET target is not the return address', pc)
                pc = stack.pop()
            else:
                raise ValueError(op)


# ---------------------------------------------------------------- priced kernels

CONST, SLOT = 1 << 32, 32


def limbs(x):
    return [(x >> (64 * i)) & M64 for i in range(4)]


def gen_fp_mul(pr, name, co):
    """4-limb CIOS Montgomery product. x0=dst, x1=a, x2=b, modulus at [x28+co], n0 at [x28+co+32]."""
    e = pr.emit
    pr.label(name)
    e('stp_pre', 'x19', 'x20', 'sp', -16)
    e('ldp', 'x3', 'x4', 'x1', 0); e('ldp', 'x5', 'x6', 'x1', 16); e('ldp', 'x7', 'x8', 'x2', 0)
    e('ldp', 'x9', 'x10', 'x28', co); e('ldp', 'x11', 'x12', 'x28', co + 16); e('ldr', 'x13', 'x28', co + 32)
    A, Pm, t, m, L = ['x3', 'x4', 'x5', 'x6'], ['x9', 'x10', 'x11', 'x12'], ['x14', 'x15', 'x16', 'x17', 'x1'], 'x19', 'x20'
    for i in range(4):
        if i == 2:
            e('ldp', 'x7', 'x8', 'x2', 16)
        b = 'x7' if i % 2 == 0 else 'x8'
        if i == 0:
            for j in range(4):
                e('mul', t[j], A[j], b)
            e('umulh', L, A[0], b); e('adds', t[1], t[1], L)
            e('umulh', L, A[1], b); e('adcs', t[2], t[2], L)
            e('umulh', L, A[2], b); e('adcs', t[3], t[3], L)
            e('umulh', L, A[3], b); e('adc', t[4], L, 'xzr')
        else:
            e('mul', L, A[0], b); e('adds', t[0], t[0], L)
            e('mul', L, A[1], b); e('adcs', t[1], t[1], L)
            e('mul', L, A[2], b); e('adcs', t[2], t[2], L)
            e('mul', L, A[3], b); e('adcs', t[3], t[3], L)
            e('adc', t[4], 'xzr', 'xzr')
            e('umulh', L, A[0], b); e('adds', t[1], t[1], L)
            e('umulh', L, A[1], b); e('adcs', t[2], t[2], L)
            e('umulh', L, A[2], b); e('adcs', t[3], t[3], L)
            e('umulh', L, A[3], b); e('adc', t[4], t[4], L)
        e('mul', m, t[0], 'x13')
        e('mul', L, m, Pm[0]); e('adds', t[0], t[0], L)
        e('mul', L, m, Pm[1]); e('adcs', t[1], t[1], L)
        e('mul', L, m, Pm[2]); e('adcs', t[2], t[2], L)
        e('mul', L, m, Pm[3]); e('adcs', t[3], t[3], L)
        e('adc', t[4], t[4], 'xzr')
        e('umulh', L, m, Pm[0]); e('adds', t[1], t[1], L)
        e('umulh', L, m, Pm[1]); e('adcs', t[2], t[2], L)
        e('umulh', L, m, Pm[2]); e('adcs', t[3], t[3], L)
        e('umulh', L, m, Pm[3]); e('adc', t[4], t[4], L)
        t = t[1:] + t[:1]
    e('subs', A[0], t[0], Pm[0]); e('sbcs', A[1], t[1], Pm[1])
    e('sbcs', A[2], t[2], Pm[2]); e('sbcs', A[3], t[3], Pm[3])
    for j in range(4):
        e('csel', t[j], A[j], t[j], 'cs')
    e('stp', t[0], t[1], 'x0', 0); e('stp', t[2], t[3], 'x0', 16)
    e('ldp_post', 'x19', 'x20', 'sp', 16)
    e('ret')


def gen_fp_add(pr, name, co):
    e = pr.emit
    pr.label(name)
    e('ldp', 'x3', 'x4', 'x1', 0); e('ldp', 'x5', 'x6', 'x1', 16)
    e('ldp', 'x7', 'x8', 'x2', 0); e('ldp', 'x9', 'x10', 'x2', 16)
    e('adds', 'x3', 'x3', 'x7'); e('adcs', 'x4', 'x4', 'x8'); e('adcs', 'x5', 'x5', 'x9'); e('adc', 'x6', 'x6', 'x10')
    e('ldp', 'x11', 'x12', 'x28', co); e('ldp', 'x13', 'x14', 'x28', co + 16)
    e('subs', 'x7', 'x3', 'x11'); e('sbcs', 'x8', 'x4', 'x12'); e('sbcs', 'x9', 'x5', 'x13'); e('sbcs', 'x10', 'x6', 'x14')
    for a, b in (('x3', 'x7'), ('x4', 'x8'), ('x5', 'x9'), ('x6', 'x10')):
        e('csel', a, b, a, 'cs')
    e('stp', 'x3', 'x4', 'x0', 0); e('stp', 'x5', 'x6', 'x0', 16)
    e('ret')


def gen_fp_sub(pr, name, co):
    e = pr.emit
    pr.label(name)
    e('ldp', 'x3', 'x4', 'x1', 0); e('ldp', 'x5', 'x6', 'x1', 16)
    e('ldp', 'x7', 'x8', 'x2', 0); e('ldp', 'x9', 'x10', 'x2', 16)
    e('subs', 'x3', 'x3', 'x7'); e('sbcs', 'x4', 'x4', 'x8'); e('sbcs', 'x5', 'x5', 'x9'); e('sbcs', 'x6', 'x6', 'x10')
    e('ldp', 'x11', 'x12', 'x28', co); e('ldp', 'x13', 'x14', 'x28', co + 16)
    for a in ('x11', 'x12', 'x13', 'x14'):
        e('csel', a, 'xzr', a, 'cs')
    e('adds', 'x3', 'x3', 'x11'); e('adcs', 'x4', 'x4', 'x12'); e('adcs', 'x5', 'x5', 'x13'); e('adc', 'x6', 'x6', 'x14')
    e('stp', 'x3', 'x4', 'x0', 0); e('stp', 'x5', 'x6', 'x0', 16)
    e('ret')


class Tracer:
    """Emits Fp operations on 32-byte memory slots (base register, byte offset <= 480, so every LDP/STP
    immediate is encodable). Multiplications call fp_mul; additions and subtractions are inlined.
    Every slot carries a static bound B (value < B*p), checked before each operation."""

    def __init__(self, pr, bases, skip=0):
        self.pr, self.used, self.peak = pr, 0, 0
        self.free = [(b, k * SLOT) for b in bases for k in range(16)][::-1]
        self.free = [s for s in self.free if not (s[0] == bases[0] and s[1] < skip)]
        self.scratch = set(bases)
        self.calls = {'M': 0, 'A': 0, 'U': 0}
        self.bound = {}

    def B(self, s):
        return self.bound.get(s, 1)          # inputs, outputs and constants are canonical

    def alloc(self):
        self.used += 1; self.peak = max(self.peak, self.used)
        return self.free.pop()

    def rel(self, s):
        if s[0] in self.scratch:
            self.used -= 1; self.free.append(s); self.bound.pop(s, None)

    def addr(self, reg, sl):
        if sl[1] == 0:
            self.pr.emit('mov', reg, sl[0])
        else:
            self.pr.emit('add', reg, sl[0], sl[1])

    def mul(self, a, b, d=None):
        ba, bb = self.B(a), self.B(b)
        assert (ba * P) * (bb * P) < P << 256, (ba, bb)       # CIOS output < 2p, one subtraction
        assert (max(ba, bb) + 1) * P < 1 << 256                # accumulator <= 3p - 1 fits the 5-word bound
        d = d if d is not None else self.alloc()
        self.addr('x0', d); self.addr('x1', a); self.addr('x2', b)
        self.pr.emit('bl', 'fp_mul')
        self.calls['M'] += 1
        self.bound[d] = 1
        return d

    def _load(self, a, b):
        e = self.pr.emit
        e('ldp', 'x3', 'x4', a[0], a[1]); e('ldp', 'x5', 'x6', a[0], a[1] + 16)
        e('ldp', 'x7', 'x8', b[0], b[1]); e('ldp', 'x9', 'x10', b[0], b[1] + 16)

    def _store(self, d):
        self.pr.emit('stp', 'x3', 'x4', d[0], d[1]); self.pr.emit('stp', 'x5', 'x6', d[0], d[1] + 16)

    def add(self, a, b, d=None):
        """(a + b) mod p for canonical a, b: 20 instructions."""
        assert self.B(a) == 1 and self.B(b) == 1
        e = self.pr.emit
        d = d if d is not None else self.alloc()
        self._load(a, b)
        e('adds', 'x3', 'x3', 'x7'); e('adcs', 'x4', 'x4', 'x8'); e('adcs', 'x5', 'x5', 'x9'); e('adc', 'x6', 'x6', 'x10')
        e('ldp', 'x11', 'x12', 'x28', 0); e('ldp', 'x13', 'x14', 'x28', 16)
        e('subs', 'x7', 'x3', 'x11'); e('sbcs', 'x8', 'x4', 'x12'); e('sbcs', 'x9', 'x5', 'x13'); e('sbcs', 'x10', 'x6', 'x14')
        for r_, t_ in (('x3', 'x7'), ('x4', 'x8'), ('x5', 'x9'), ('x6', 'x10')):
            e('csel', r_, t_, r_, 'cs')
        self._store(d)
        self.calls['A'] += 1; self.bound[d] = 1
        return d

    def sub(self, a, b, d=None):
        """(a - b) mod p for canonical a, b: 20 instructions."""
        assert self.B(a) == 1 and self.B(b) == 1
        e = self.pr.emit
        d = d if d is not None else self.alloc()
        self._load(a, b)
        e('subs', 'x3', 'x3', 'x7'); e('sbcs', 'x4', 'x4', 'x8'); e('sbcs', 'x5', 'x5', 'x9'); e('sbcs', 'x6', 'x6', 'x10')
        e('ldp', 'x11', 'x12', 'x28', 0); e('ldp', 'x13', 'x14', 'x28', 16)
        for r_ in ('x11', 'x12', 'x13', 'x14'):
            e('csel', r_, 'xzr', r_, 'cs')
        e('adds', 'x3', 'x3', 'x11'); e('adcs', 'x4', 'x4', 'x12'); e('adcs', 'x5', 'x5', 'x13'); e('adc', 'x6', 'x6', 'x14')
        self._store(d)
        self.calls['A'] += 1; self.bound[d] = 1
        return d

    def addu(self, a, b, d=None):
        """Unreduced a + b (no modular correction): 10 instructions; result bound B(a)+B(b)."""
        nb = self.B(a) + self.B(b)
        assert nb * P < 1 << 256
        e = self.pr.emit
        d = d if d is not None else self.alloc()
        self._load(a, b)
        e('adds', 'x3', 'x3', 'x7'); e('adcs', 'x4', 'x4', 'x8'); e('adcs', 'x5', 'x5', 'x9'); e('adc', 'x6', 'x6', 'x10')
        self._store(d)
        self.calls['U'] += 1; self.bound[d] = nb
        return d


def tower_slots(base):
    s = [(base, k * SLOT) for k in range(12)]
    return (((s[0], s[1]), (s[2], s[3]), (s[4], s[5])), ((s[6], s[7]), (s[8], s[9]), (s[10], s[11])))


def gen_f12_mul(pr):
    """x0=out, x1=a, x2=b: tower Karatsuba (Fp2, Fp6, Fp12), xi=9+u by additions; additions inlined,
    the two Karatsuba sums inside each Fp2 product left unreduced (< 2p)."""
    e = pr.emit
    pr.label('f12_mul')
    e('stp_pre', 'x19', 'x20', 'sp', -16); e('stp_pre', 'x21', 'x22', 'sp', -16)
    e('stp_pre', 'x23', 'x24', 'sp', -16); e('stp_pre', 'x25', 'x30', 'sp', -16)
    e('sub_sp', 4096)
    e('mov', 'x19', 'x1'); e('mov', 'x20', 'x2'); e('mov', 'x21', 'x0'); e('mov', 'x22', 'sp')
    e('add', 'x23', 'sp', 512); e('add', 'x24', 'sp', 1024)
    T = Tracer(pr, ['x22', 'x23', 'x24'])
    mul, add, sub, addu = T.mul, T.add, T.sub, T.addu

    def f2add(a, b, d=(None, None)):
        return (add(a[0], b[0], d[0]), add(a[1], b[1], d[1]))

    def f2sub(a, b, d=(None, None)):
        return (sub(a[0], b[0], d[0]), sub(a[1], b[1], d[1]))

    def free2(x):
        T.rel(x[0]); T.rel(x[1])

    def f2mul(a, b):
        v0 = mul(a[0], b[0]); v1 = mul(a[1], b[1])
        sa = addu(a[0], a[1]); sb = addu(b[0], b[1]); t = mul(sa, sb); T.rel(sa); T.rel(sb)
        r0 = sub(v0, v1); t = sub(t, v0, t); r1 = sub(t, v1, t)
        T.rel(v0); T.rel(v1)
        return (r0, r1)

    def f2xi(a, d=(None, None)):
        def nine(x):
            y = add(x, x); y = add(y, y, y); y = add(y, y, y)
            return add(y, x, y)
        n0 = nine(a[0]); n1 = nine(a[1])
        r0 = sub(n0, a[1], d[0]); r1 = add(n1, a[0], d[1])
        T.rel(n0); T.rel(n1)
        return (r0, r1)

    def f6mul(a, b):
        v0, v1, v2 = f2mul(a[0], b[0]), f2mul(a[1], b[1]), f2mul(a[2], b[2])
        out = []
        for (i, j) in ((1, 2), (0, 1), (0, 2)):
            s1 = f2add(a[i], a[j]); s2 = f2add(b[i], b[j]); out.append(f2mul(s1, s2)); free2(s1); free2(s2)
        t0, t1, t2 = out
        t0 = f2sub(t0, v1, t0); t0 = f2sub(t0, v2, t0); x = f2xi(t0); free2(t0)
        c0 = f2add(v0, x, x)
        t1 = f2sub(t1, v0, t1); t1 = f2sub(t1, v1, t1); y = f2xi(v2); c1 = f2add(t1, y, t1); free2(y)
        t2 = f2sub(t2, v0, t2); t2 = f2sub(t2, v2, t2); c2 = f2add(t2, v1, t2)
        free2(v0); free2(v1); free2(v2)
        return (c0, c1, c2)

    A, B, O = tower_slots('x19'), tower_slots('x20'), tower_slots('x21')
    v0 = f6mul(A[0], B[0]); v1 = f6mul(A[1], B[1])
    sa = tuple(f2add(A[0][k], A[1][k]) for k in range(3))
    sb = tuple(f2add(B[0][k], B[1][k]) for k in range(3))
    t = f6mul(sa, sb)
    for k in range(3):
        f2sub(t[k], v0[k], t[k]); f2sub(t[k], v1[k], O[1][k])
    xv = f2xi(v1[2])
    f2add(v0[0], xv, O[0][0]); f2add(v0[1], v1[0], O[0][1]); f2add(v0[2], v1[1], O[0][2])
    e('add_sp', 4096)
    e('ldp_post', 'x25', 'x30', 'sp', 16); e('ldp_post', 'x23', 'x24', 'sp', 16)
    e('ldp_post', 'x21', 'x22', 'sp', 16); e('ldp_post', 'x19', 'x20', 'sp', 16)
    e('ret')
    return T


# ---------------------------------------------------------------- reference tower arithmetic

class Tower:
    """Fq2=Fq[u]/(u^2+1), xi=a0+u, Fq6=Fq2[v]/(v^3-xi), Fq12=Fq6[w]/(w^2-v); y=((e0,e2,e4),(e1,e3,e5))."""

    def __init__(self, q, a0):
        self.q, self.a0, self.xi = q, a0, (a0 % q, 1)
        self.g1 = self.f2pow(self.xi, (q - 1) // 6)
        self.gk = [self.f2pow(self.g1, k) for k in range(6)]
        z = self.f2pow(self.xi, (q * q - 1) // 6)
        assert z[1] == 0
        self.zeta = z[0]
        self.one = (((1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)))

    def f2add(self, a, b): q = self.q; return ((a[0] + b[0]) % q, (a[1] + b[1]) % q)
    def f2sub(self, a, b): q = self.q; return ((a[0] - b[0]) % q, (a[1] - b[1]) % q)
    def f2mul(self, a, b): q = self.q; return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

    def f2pow(self, a, n):
        r, b = (1, 0), a
        while n:
            if n & 1: r = self.f2mul(r, b)
            b = self.f2mul(b, b); n >>= 1
        return r

    def f6mul(self, a, b):
        m, ad, sb, xi = self.f2mul, self.f2add, self.f2sub, lambda x: self.f2mul(x, self.xi)
        v0, v1, v2 = m(a[0], b[0]), m(a[1], b[1]), m(a[2], b[2])
        t0 = m(ad(a[1], a[2]), ad(b[1], b[2])); t1 = m(ad(a[0], a[1]), ad(b[0], b[1])); t2 = m(ad(a[0], a[2]), ad(b[0], b[2]))
        return (ad(v0, xi(sb(sb(t0, v1), v2))), ad(sb(sb(t1, v0), v1), xi(v2)), ad(sb(sb(t2, v0), v2), v1))

    def mul(self, a, b):
        v0, v1 = self.f6mul(a[0], b[0]), self.f6mul(a[1], b[1])
        t = self.f6mul(tuple(map(self.f2add, a[0], a[1])), tuple(map(self.f2add, b[0], b[1])))
        vv1 = (self.f2mul(v1[2], self.xi), v1[0], v1[1])
        return (tuple(map(self.f2add, v0, vv1)), tuple(map(self.f2sub, tuple(map(self.f2sub, t, v0)), v1)))

    def pow(self, a, n):
        r, b = self.one, a
        while n:
            if n & 1: r = self.mul(r, b)
            b = self.mul(b, b); n >>= 1
        return r

    def frob(self, y):
        e = [y[0][0], y[1][0], y[0][1], y[1][1], y[0][2], y[1][2]]
        f = [self.f2mul((x[0], (-x[1]) % self.q), self.gk[k]) for k, x in enumerate(e)]
        return ((f[0], f[2], f[4]), (f[1], f[3], f[5]))

    def frob_k(self, y, k):
        for _ in range(k % 12):
            y = self.frob(y)
        return y

    def to_tower(self, c):
        q, a0 = self.q, self.a0
        e = [((c[k] + a0 * c[k + 6]) % q, c[k + 6] % q) for k in range(6)]
        return ((e[0], e[2], e[4]), (e[1], e[3], e[5]))

    def from_tower(self, y):
        e = [y[0][0], y[1][0], y[0][1], y[1][1], y[0][2], y[1][2]]
        return [(x[0] - self.a0 * x[1]) % self.q for x in e] + [x[1] for x in e]

    def schoolbook(self, a, b):
        q, a0, c = self.q, self.a0, [0] * 23
        for i in range(12):
            for j in range(12):
                c[i + j] += a[i] * b[j]
        for k in range(22, 11, -1):
            c[k - 6] += 2 * a0 * c[k]; c[k - 12] -= (a0 * a0 + 1) * c[k]
        return [x % q for x in c[:12]]


def mem_order(y):
    return [y[i][j][k] for i in range(2) for j in range(3) for k in range(2)]


def from_mem(v):
    it = iter(v)
    return tuple(tuple((next(it), next(it)) for _ in range(3)) for _ in range(2))


# ---------------------------------------------------------------- checks

def is_prime(n):
    if n < 2: return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0: return n == sp
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True
