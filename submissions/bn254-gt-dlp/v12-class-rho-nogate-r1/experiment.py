"""Supporting checks for the BN254 GT class-rho claim. Python standard library only.

Modes (input JSON on stdin):
  kernels  generate the priced AArch64 programs (Fp multiply/add/sub, Fp12 product,
           class selector), execute them on an instruction interpreter, check them
           against big-integer arithmetic, and print the executed instruction counts.
  algebra  check the algebraic facts of proof.md on the pinned BN254 constants.
  ledger   recompute every integer of the cost ledger, the score, the trap sums.
  toy      run the complete algorithm of proof.md in a BN-shaped toy field.
  all      everything above.
Scope: toy runs are evidence at toy scale only. Nothing here runs a production attack.
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext, ROUND_CEILING

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
        if op in ('lsr', 'eor_lsr'):
            assert 0 <= a[-1] <= 63, a                       # shift amount
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
        return self.mem.get(ad, 0)

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
            elif op == 'eor':
                s(a[0], g(a[1]) ^ g(a[2]))
            elif op == 'add_r':
                s(a[0], g(a[1]) + g(a[2]))
            elif op == 'eor_lsr':
                s(a[0], g(a[1]) ^ (g(a[2]) >> a[3]))
            elif op == 'and_r':
                s(a[0], g(a[1]) & g(a[2]))
            elif op == 'lsr':
                s(a[0], g(a[1]) >> a[2])
            elif op == 'madd':
                s(a[0], g(a[1]) * g(a[2]) + g(a[3]))
            elif op == 'bcond':                  # conditional branch on the measured path: counted, must fall through
                assert not self.cond(a[0]), ('branch taken on the accept path', a)
            elif op == 'cbnz':                   # counted; the exit itself is charged in other ledger rows
                pass
            elif op == 'cmp':
                t = g(a[0]) - g(a[1]); self.C, self.Z = int(t >= 0), int(t & M64 == 0)
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


def gen_selector(pr):
    """x0=record (K: 4 words, k*, tie), x1=y (tower). Keys V_k = Re(c1 coefficient of sigma^k(y))
    from 4 multiplications and 9 additions/subtractions; then the minimum, its index and a tie flag."""
    e = pr.emit
    pr.label('selector')
    e('stp_pre', 'x19', 'x20', 'sp', -16); e('stp_pre', 'x21', 'x22', 'sp', -16); e('stp_pre', 'x26', 'x30', 'sp', -16)
    e('sub_sp', 4096)
    e('mov', 'x19', 'x1'); e('mov', 'x21', 'x0'); e('mov', 'x22', 'sp'); e('add', 'x26', 'x28', 256)
    T = Tracer(pr, ['x22'], skip=12 * SLOT)
    mul, add, sub = T.mul, T.add, T.sub
    a, b = ('x19', 6 * SLOT), ('x19', 7 * SLOT)
    ZETA, Dc, Ec, ZERO = ('x26', 0), ('x26', 32), ('x26', 64), ('x26', 96)
    V = [('x22', k * SLOT) for k in range(12)]
    t1 = mul(a, Dc); t2 = mul(b, Ec); add(t1, t2, V[1]); T.rel(t1); T.rel(t2)
    mul(a, ZETA, V[2]); mul(V[1], ZETA, V[3])
    sub(V[2], a, V[4]); sub(V[3], V[1], V[5])
    for k, src in ((6, a), (7, V[1]), (8, V[2]), (9, V[3]), (10, V[4]), (11, V[5])):
        sub(ZERO, src, V[k])
    e('ldp', 'x3', 'x4', 'x19', 6 * SLOT); e('ldp', 'x5', 'x6', 'x19', 6 * SLOT + 16)
    e('movi', 'x7', 0); e('mov', 'x8', 'xzr')
    for k in range(1, 12):
        e('ldp', 'x9', 'x10', 'x22', k * SLOT); e('ldp', 'x11', 'x12', 'x22', k * SLOT + 16)
        e('subs', 'x13', 'x9', 'x3'); e('sbcs', 'x14', 'x10', 'x4'); e('sbcs', 'x15', 'x11', 'x5'); e('sbcs', 'x16', 'x12', 'x6')
        e('orr', 'x13', 'x13', 'x14'); e('orr', 'x13', 'x13', 'x15'); e('orr', 'x13', 'x13', 'x16')
        e('cset', 'x17', 'cc')
        e('cmp', 'x13', 'xzr'); e('cset', 'x14', 'eq'); e('orr', 'x8', 'x8', 'x14')
        e('cmp', 'x17', 'xzr'); e('csel', 'x8', 'x8', 'xzr', 'eq')
        e('csel', 'x3', 'x3', 'x9', 'eq'); e('csel', 'x4', 'x4', 'x10', 'eq')
        e('csel', 'x5', 'x5', 'x11', 'eq'); e('csel', 'x6', 'x6', 'x12', 'eq')
        e('movi', 'x15', k); e('csel', 'x7', 'x7', 'x15', 'eq')
    e('stp', 'x3', 'x4', 'x21', 0); e('stp', 'x5', 'x6', 'x21', 16); e('stp', 'x7', 'x8', 'x21', 32)
    e('add_sp', 4096)
    e('ldp_post', 'x26', 'x30', 'sp', 16); e('ldp_post', 'x21', 'x22', 'sp', 16); e('ldp_post', 'x19', 'x20', 'sp', 16)
    e('ret')
    return T


GLUE_POOL = 512           # offset of the walk constants in the pool
GLUE_ST = 2048            # offset of the walk state block in the pool region (counters, window, buffers)


def gen_glue(pr):
    """One walk slot outside f12_mul and selector, branch-free: the same instructions run whether
    the probe is accepted, hits the window, ties or reaches a distinguished key; results are
    committed with CSEL. Registers: x19 = state record (K, k*, tie) with the candidate record at
    +64, x20 = window (24 hash words, keys at +256), x21 = probe index, x22 = y, x23 = next buffer,
    x24 = alpha, beta (temporaries at +64), x25 = counters/flags, x26 = pool (x28 + 512)."""
    e = pr.emit
    pr.label('glue')
    e('stp_pre', 'x29', 'x30', 'sp', -16); e('stp_pre', 'x26', 'x27', 'sp', -16)
    e('add', 'x26', 'x28', GLUE_POOL)
    e('ldp', 'x13', 'x14', 'x26', 32); e('ldp', 'x15', 'x16', 'x26', 48); e('ldp', 'x6', 'x7', 'x26', 64)
    e('ldr', 'x0', 'x25', 48); e('ldr', 'x1', 'x19', 32)
    e('add_r', 'x2', 'x0', 'x21'); e('and_r', 'x2', 'x2', 'x13'); e('madd', 'x2', 'x2', 'x14', 'x1')
    e('madd', 'x3', 'x2', 'x15', 'x6'); e('movi', 'x4', 64); e('madd', 'x27', 'x2', 'x4', 'x7')
    e('mov', 'x0', 'x23'); e('mov', 'x1', 'x22'); e('mov', 'x2', 'x3'); e('bl', 'f12_mul')
    e('add', 'x0', 'x19', 64); e('mov', 'x1', 'x23'); e('bl', 'selector')
    e('ldp', 'x9', 'x10', 'x26', 0); e('ldp', 'x11', 'x12', 'x26', 16)
    e('ldp', 'x2', 'x3', 'x19', 64); e('ldp', 'x4', 'x5', 'x19', 80)
    for hreg, salt in (('x0', 'x9'), ('x1', 'x10')):
        e('mov', hreg, salt)
        for limb in ('x2', 'x3', 'x4', 'x5'):
            e('eor', hreg, hreg, limb); e('eor_lsr', hreg, hreg, hreg, 30); e('mul', hreg, hreg, 'x11')
            e('eor_lsr', hreg, hreg, hreg, 27); e('mul', hreg, hreg, 'x12'); e('eor_lsr', hreg, hreg, hreg, 31)
    e('str', 'x0', 'x25', 72); e('str', 'x1', 'x25', 64)
    e('mov', 'x8', 'xzr')                              # window hit flag over all 24 stored hash words
    for j in range(24):
        e('ldr', 'x6', 'x20', 8 * j); e('cmp', 'x6', 'x0'); e('cset', 'x7', 'eq'); e('orr', 'x8', 'x8', 'x7')
    e('ldr', 'x9', 'x19', 104); e('orr', 'x9', 'x9', 'x8'); e('cmp', 'x9', 'xzr'); e('cset', 'x10', 'eq')
    e('str', 'x10', 'x25', 96)                         # accept flag
    e('add', 'x0', 'x24', 64); e('mov', 'x1', 'x24'); e('mov', 'x2', 'x27'); e('bl', 'fr_add')
    e('add', 'x0', 'x24', 96); e('add', 'x1', 'x24', 32); e('add', 'x2', 'x27', 32); e('bl', 'fr_add')
    e('ldr', 'x10', 'x25', 96); e('cmp', 'x10', 'xzr')  # eq: reject, keep the old value
    for off in (0, 16, 32, 48):
        e('ldp', 'x2', 'x3', 'x24', 64 + off); e('ldp', 'x4', 'x5', 'x24', off)
        e('csel', 'x2', 'x4', 'x2', 'eq'); e('csel', 'x3', 'x5', 'x3', 'eq'); e('stp', 'x2', 'x3', 'x24', off)
    for off in (0, 16, 32):                            # commit the candidate record (K, k*, tie)
        e('ldp', 'x2', 'x3', 'x19', 64 + off); e('ldp', 'x4', 'x5', 'x19', off)
        e('csel', 'x2', 'x4', 'x2', 'eq'); e('csel', 'x3', 'x5', 'x3', 'eq'); e('stp', 'x2', 'x3', 'x19', off)
    e('ldr', 'x0', 'x25', 56); e('movi', 'x7', 8); e('madd', 'x1', 'x0', 'x7', 'x20')
    e('movi', 'x7', 32); e('madd', 'x2', 'x0', 'x7', 'x20'); e('add', 'x2', 'x2', 256)
    e('ldr', 'x3', 'x1', 0); e('ldr', 'x4', 'x25', 72); e('csel', 'x3', 'x3', 'x4', 'eq'); e('str', 'x3', 'x1', 0)
    for off in (0, 16):
        e('ldp', 'x3', 'x4', 'x2', off); e('ldp', 'x5', 'x6', 'x19', off)
        e('csel', 'x3', 'x3', 'x5', 'eq'); e('csel', 'x4', 'x4', 'x6', 'eq'); e('stp', 'x3', 'x4', 'x2', off)
    e('add', 'x3', 'x0', 1); e('movi', 'x7', 24); e('cmp', 'x3', 'x7'); e('csel', 'x3', 'xzr', 'x3', 'cs')
    e('cmp', 'x10', 'xzr'); e('csel', 'x0', 'x0', 'x3', 'eq'); e('str', 'x0', 'x25', 56)
    e('ldr', 'x13', 'x26', 32)                         # reload R-1 (caller-saved registers were clobbered)
    e('ldr', 'x4', 'x25', 72); e('and_r', 'x5', 'x4', 'x13'); e('ldr', 'x6', 'x25', 48)
    e('csel', 'x5', 'x6', 'x5', 'eq'); e('str', 'x5', 'x25', 48)
    e('mov', 'x0', 'x22'); e('csel', 'x22', 'x22', 'x23', 'eq'); e('csel', 'x23', 'x23', 'x0', 'eq')
    e('add', 'x0', 'x21', 1)                           # next probe index: a+1, or 0 on accept
    e('csel', 'x21', 'x0', 'xzr', 'eq')
    e('ldr', 'x8', 'x26', 80)                          # distinguished flag of the new key, only if accepted
    e('lsr', 'x2', 'x4', 14); e('cmp', 'x2', 'xzr'); e('cset', 'x3', 'eq')
    e('ldr', 'x5', 'x25', 64); e('and_r', 'x5', 'x5', 'x8'); e('cmp', 'x5', 'xzr'); e('cset', 'x5', 'eq')
    e('and_r', 'x3', 'x3', 'x5'); e('and_r', 'x3', 'x3', 'x10'); e('str', 'x3', 'x25', 104)
    e('mov', 'x9', 'xzr')                              # cap flags
    for off, cap, inc in ((0, 96, None), (16, 112, None), (32, 128, 'x10')):
        e('ldp', 'x0', 'x1', 'x25', off)
        if inc is None:
            e('movi', 'x2', 1); inc_r = 'x2'
        else:
            inc_r = inc
        e('adds', 'x0', 'x0', inc_r); e('adc', 'x1', 'x1', 'xzr'); e('stp', 'x0', 'x1', 'x25', off)
        e('ldp', 'x5', 'x6', 'x26', cap); e('subs', 'xzr', 'x0', 'x5'); e('sbcs', 'xzr', 'x1', 'x6')
        e('cset', 'x7', 'cs'); e('orr', 'x9', 'x9', 'x7')
    e('movi', 'x2', 8); e('cmp', 'x21', 'x2'); e('cset', 'x7', 'cs'); e('orr', 'x9', 'x9', 'x7')
    e('ldr', 'x7', 'x19', 104); e('orr', 'x9', 'x9', 'x7'); e('orr', 'x9', 'x9', 'x3')
    e('str', 'x9', 'x25', 112); e('cbnz', 'x9')        # end of trail / report / cap: one counted branch
    e('ldp_post', 'x26', 'x27', 'sp', 16); e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')


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

    def keys(self, y):
        q, z = self.q, self.zeta
        a, b = y[1][0]
        B = (a * self.g1[0] + b * self.g1[1]) % q
        az, Bz = a * z % q, B * z % q
        s1, s2 = (az - a) % q, (Bz - B) % q
        return [a, B, az, Bz, s1, s2, -a % q, -B % q, -az % q, -Bz % q, -s1 % q, -s2 % q]

    def select(self, y):
        V = self.keys(y)
        K = min(V)
        return None if V.count(K) != 1 else (K, V.index(K))

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


def kernels(seed):
    rng = random.Random(seed)
    pr = Program()
    gen_fp_mul(pr, 'fp_mul', 0); gen_fp_add(pr, 'fp_add', 0); gen_fp_sub(pr, 'fp_sub', 0)
    gen_fp_mul(pr, 'fr_mul', 64); gen_fp_add(pr, 'fr_add', 64)
    t12 = gen_f12_mul(pr); tsel = gen_selector(pr); gen_glue(pr)
    mc = Machine(pr)
    TW = Tower(P, 9)
    RM = 1 << 256
    mont, unmont = (lambda x: x * RM % P), (lambda x: x * pow(RM, -1, P) % P)

    def setc(addr, x):
        for k, w in enumerate(limbs(x)): mc.st(addr + 8 * k, w)

    def getv(addr):
        return sum(mc.ld(addr + 8 * k) << (64 * k) for k in range(4))
    for k, w in enumerate(limbs(P)): mc.st(CONST + 8 * k, w)
    mc.st(CONST + 32, (-pow(P, -1, 1 << 64)) & M64)
    for k, w in enumerate(limbs(R)): mc.st(CONST + 64 + 8 * k, w)
    mc.st(CONST + 96, (-pow(R, -1, 1 << 64)) & M64)
    for j, v in enumerate((mont(TW.zeta), mont(TW.g1[0]), mont(TW.g1[1]), 0)):
        setc(CONST + 256 + 32 * j, v)
    mc.r['x28'] = CONST
    A, B, O = 1 << 20, (1 << 20) + 64, (1 << 20) + 128
    edge = [0, 1, 2, P - 1, P - 2, (1 << 255) % P, RM % P, RM * RM % P]
    bad, cnt = 0, {}
    for i in range(400):
        x = edge[i % 8] if i < 64 else rng.randrange(P)
        y = edge[(i // 8) % 8] if i < 64 else rng.randrange(P)
        setc(A, x); setc(B, y)
        for fn, ref in (('fp_mul', x * y * pow(RM, -1, P) % P), ('fp_add', (x + y) % P), ('fp_sub', (x - y) % P)):
            mc.r['x0'], mc.r['x1'], mc.r['x2'] = O, A, B
            cnt[fn] = mc.call(fn)
            bad += getv(O) != ref
        xr, yr = x % R, y % R
        setc(A, xr); setc(B, yr)
        for fn, ref in (('fr_mul', xr * yr * pow(RM, -1, R) % R), ('fr_add', (xr + yr) % R)):
            mc.r['x0'], mc.r['x1'], mc.r['x2'] = O, A, B
            cnt[fn] = mc.call(fn)
            bad += getv(O) != ref
        xu = (2 * P - 1 - (i % 3)) if i < 24 else rng.randrange(2 * P)      # unreduced inputs < 2p
        yu = (2 * P - 1 - (i % 5)) if i < 24 else rng.randrange(2 * P)
        setc(A, xu); setc(B, yu)
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = O, A, B
        mc.call('fp_mul'); bad += getv(O) != xu * yu * pow(RM, -1, P) % P
        setc(A, x); setc(B, y)
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = A, A, B      # destination aliases an input
        mc.call('fp_mul'); bad += getv(A) != x * y * pow(RM, -1, P) % P
    g = TW.to_tower(G_PINNED)
    YA, YB, YO, KO = 1 << 24, (1 << 24) + 4096, (1 << 24) + 8192, (1 << 24) + 12288
    bad12 = badsel = 0
    for i in range(12):
        xa, xb = [rng.randrange(P) for _ in range(12)], [rng.randrange(P) for _ in range(12)]
        a, b = TW.to_tower(xa), TW.to_tower(xb)
        for k, v in enumerate(mem_order(a)): setc(YA + 32 * k, mont(v))
        for k, v in enumerate(mem_order(b)): setc(YB + 32 * k, mont(v))
        mc.r['x0'], mc.r['x1'], mc.r['x2'] = YO, YA, YB
        n12 = mc.call('f12_mul')
        out = from_mem([unmont(getv(YO + 32 * k)) for k in range(12)])
        bad12 += out != TW.mul(a, b) or TW.from_tower(out) != TW.schoolbook(xa, xb)
        y = TW.pow(g, rng.randrange(1, R))
        for k, v in enumerate(mem_order(y)): setc(YA + 32 * k, mont(v))
        mc.r['x0'], mc.r['x1'] = KO, YA
        nsel = mc.call('selector')
        Km, idx, tie = getv(KO), mc.ld(KO + 32), mc.ld(KO + 40)
        Vx = [mont(TW.frob_k(y, k)[1][0][0]) for k in range(12)]      # from explicit conjugates
        K = min(Vx)
        badsel += tie != (Vx.count(K) > 1) or (not tie and (Km != K or idx != Vx.index(K)))
    # glue: one full slot on its longest path, checked against Python
    G0 = 1 << 28
    Y, NX, REC, WIN, AB, CNT, TWB, EXB = G0, G0 + 4096, G0 + 8192, G0 + 12288, G0 + 16384, G0 + 20480, 1 << 30, 1 << 31
    salt1, salt2 = rng.getrandbits(64), rng.getrandbits(64)
    C1, C2 = 0xbf58476d1ce4e5b9, 0x94d049bb133111eb
    pool = [salt1, salt2, C1, C2, (1 << 14) - 1, 12, 384, 40, TWB, EXB, (1 << 37) - 1, 0,
            0, 1, 0, 1, 0, 1]                               # caps 2^64 for slots, trail slots, steps (lo, hi words)
    for k_, w_ in enumerate(pool):
        mc.st(CONST + GLUE_POOL + 8 * k_, w_)
    glue_bad, counts, glue_err = 0, set(), None
    for trial in range(4):
        y = TW.pow(g, rng.randrange(1, R)); T_ = TW.pow(g, rng.randrange(1, R))
        V0 = [mont(v) for v in TW.keys(y)]; Kc = min(V0); kc = V0.index(Kc)
        h, a_ = rng.randrange(1 << 14), rng.randrange(7)
        t_ = ((h + a_) & ((1 << 14) - 1)) * 12 + kc
        for k, v in enumerate(mem_order(y)): setc(Y + 32 * k, mont(v))
        for k, v in enumerate(mem_order(T_)): setc(TWB + 384 * t_ + 32 * k, mont(v))
        u_, v_, al, be = (rng.randrange(R) for _ in range(4))
        setc(EXB + 64 * t_, u_); setc(EXB + 64 * t_ + 32, v_); setc(AB, al); setc(AB + 32, be)
        setc(REC, Kc); mc.st(REC + 32, kc); mc.st(REC + 40, 0)
        nxt = TW.mul(y, T_)
        Vn = [mont(v) for v in TW.keys(nxt)]; Kn = min(Vn); kn = Vn.index(Kn)
        w1, w2 = khash(Kn, salt1), khash(Kn, salt2)
        hit = trial % 2 == 1                              # odd trials plant a window hit
        win = [rng.getrandbits(64) for _ in range(24)]
        if hit:
            win[rng.randrange(24)] = w1
        for j in range(24):
            mc.st(WIN + 8 * j, win[j]); setc(WIN + 256 + 32 * j, rng.getrandbits(250))
        c0 = [rng.getrandbits(40) for _ in range(3)]
        for j, c_ in enumerate(c0):
            mc.st(CNT + 16 * j, c_); mc.st(CNT + 16 * j + 8, 0)
        ring = rng.randrange(24)
        mc.st(CNT + 48, h); mc.st(CNT + 56, ring)
        for reg, val in (('x19', REC), ('x20', WIN), ('x21', a_), ('x22', Y), ('x23', NX), ('x24', AB), ('x25', CNT)):
            mc.r[reg] = val
        mc.r['x30'] = 0x6c756500                          # glue entry address: RET must return here
        counts.add(mc.call('glue') - n12 - nsel)
        glue_bad += mc.r['x30'] != 0x6c756500 or mc.r['sp'] != 1 << 40   # next slot starts at glue's first instruction
        acc = not hit
        dp = (w1 >> 14) == 0 and (w2 & ((1 << 37) - 1)) == 0 and acc
        ok = (mc.ld(CNT + 96) == int(acc) and mc.ld(CNT + 104) == int(dp)
              and [mc.ld(CNT + 16 * j) for j in range(3)] == [c0[0] + 1, c0[1] + 1, c0[2] + int(acc)]
              and mc.r['x21'] == (0 if acc else a_ + 1))
        if acc:
            ok &= (from_mem([unmont(getv(NX + 32 * k)) for k in range(12)]) == nxt and getv(REC) == Kn
                   and mc.ld(REC + 32) == kn and mc.ld(CNT + 48) == w1 & ((1 << 14) - 1)
                   and mc.ld(WIN + 8 * ring) == w1 and getv(WIN + 256 + 32 * ring) == Kn
                   and mc.ld(CNT + 56) == (ring + 1) % 24 and getv(AB) == (al + u_) % R
                   and getv(AB + 32) == (be + v_) % R and mc.r['x22'] == NX and mc.r['x23'] == Y)
        else:
            ok &= (getv(REC) == Kc and mc.ld(REC + 32) == kc and mc.ld(CNT + 48) == h
                   and mc.ld(WIN + 8 * ring) == win[ring] and mc.ld(CNT + 56) == ring
                   and getv(AB) == al and getv(AB + 32) == be and mc.r['x22'] == Y and mc.r['x23'] == NX)
        glue_bad += not ok
    nglue = max(counts)
    glue_bad += len(counts) != 1                          # the count must not depend on the path
    regs = set()
    for op, a in pr.code:
        regs |= {x for x in a if isinstance(x, str) and x.startswith('x') and x[1:].isdigit()}
    return {"fp_mul": cnt['fp_mul'], "fp_add": cnt['fp_add'], "fp_sub": cnt['fp_sub'],
            "fr_mul": cnt['fr_mul'], "fr_add": cnt['fr_add'], "fp_mismatches": bad,
            "f12_mul_executed": n12, "f12_calls": t12.calls, "f12_scratch_slots": t12.peak, "f12_mismatches": bad12,
            "selector_executed": nsel, "selector_calls": tsel.calls, "selector_mismatches": badsel,
            "stack_bytes": mc.depth, "registers_used": sorted(regs, key=lambda s: int(s[1:])),
            "glue_executed": nglue, "glue_mismatches": glue_bad, "glue_paths_tested": "2 accepted, 2 window hits",
            "program_instructions": len(pr.code)}


def algebra(seed):
    rng = random.Random(seed)
    TW = Tower(P, 9)
    g = TW.to_tower(G_PINNED)
    lam = P % R
    out = {"p_r_prime": is_prime(P) and is_prime(R),
           "bn_polynomials": P == 36 * U ** 4 + 36 * U ** 3 + 24 * U ** 2 + 6 * U + 1 and R == 36 * U ** 4 + 36 * U ** 3 + 18 * U ** 2 + 6 * U + 1,
           "p_mod_12": P % 12,
           "ord_r_p_is_12": all(pow(P, k, R) != 1 for k in range(1, 12)) and pow(P, 12, R) == 1,
           "lambda6_is_minus1": pow(lam, 6, R) == R - 1,
           "zeta_order_6": (TW.zeta ** 2 - TW.zeta + 1) % P == 0,
           "r_quarter_power_log2": math.log2(R) / 4,
           "g_order_r": TW.pow(g, R) == TW.one and g != TW.one,
           "frobenius_is_pow_p": TW.frob(g) == TW.pow(g, P),
           "frobenius_is_pow_lambda_on_GT": TW.frob(g) == TW.pow(g, lam)}
    keys_ok = shift_ok = step_ok = ties = 0
    trials = 6
    for _ in range(trials):
        y = TW.pow(g, rng.randrange(1, R))
        V = TW.keys(y)
        keys_ok += V == [TW.frob_k(y, k)[1][0][0] for k in range(12)]
        shift_ok += TW.keys(TW.frob(y)) == V[1:] + V[:1]
        # class-invariance of the twisted step: y and sigma^j(y) give the same next class
        T = TW.pow(g, rng.randrange(1, R))
        tw = [TW.frob_k(T, (12 - k) % 12) for k in range(12)]
        sel = TW.select(y)
        if sel is None:
            ties += 1
            continue
        j = rng.randrange(1, 12)
        y2 = TW.frob_k(y, j)
        s2 = TW.select(y2)
        n1, n2 = TW.mul(y, tw[sel[1]]), TW.mul(y2, tw[s2[1]])
        c1, c2 = TW.frob_k(n1, TW.select(n1)[1]), TW.frob_k(n2, TW.select(n2)[1])
        step_ok += c1 == c2 and s2[1] == (sel[1] - j) % 12 and TW.frob_k(y, sel[1]) == TW.frob_k(y2, s2[1])
    T = TW.pow(g, 12345)
    out.update({"keys_match_explicit_conjugates": f"{keys_ok}/{trials}", "keys_shift_under_frobenius": f"{shift_ok}/{trials}",
                "step_class_invariant": f"{step_ok}/{trials - ties}",
                "twist_exponent_identity": TW.frob_k(T, 12 - 5) == TW.pow(T, pow(lam, -5, R)),
                "inverse_is_sigma6": TW.mul(TW.frob_k(T, 6), T) == TW.one})
    return out


COIN = 16                  # instructions per 64-bit word of fresh coins (one RNDR and its status branch, store, loop)
SCALAR = 128 * (4 * COIN + 64)   # one scalar mod r: at most 128 draws of 4 coin words, mask and compare (64)


def nu_bound(Rr, Ncl, w=24):
    """Per-slot window-rejection bound: sum_{l=2..w} E(l) + w/N (+ w 2^-64 for hash-word coincidences)."""
    def E(L):
        return sum(math.comb(Rr, d) * (d / Rr) ** L * 12.0 ** (-d) for d in range(1, min(L // 2, Rr) + 1))
    return sum(E(L) for L in range(2, w + 1)) + w / Ncl + w * 2.0 ** -64


# ---------------------------------------------------------------- ledger

PARAM = {"R_index": 1 << 14, "window": 24, "probes": 8, "dp_bits": 87, "L_cap_log2": 94,
         "ST_max_log2": 40, "D_max_log2": 40, "REC_max": 16,
         "budget": "least S with S^2 >= 2 (6/5) N ln(5/3)", "mu": "6/5", "loss": "3/10000", "gate": "none"}


def ledger(k):
    getcontext().prec = 120
    C = 20000
    F12, SEL, FRM = k["f12_mul_executed"], k["selector_executed"], k["fr_mul"] + 4
    SLOTC = F12 + SEL + k["glue_executed"]
    Ri = PARAM["R_index"]
    N = (R - 1) // 12
    X_ = 2 * Decimal(6) / 5 * Decimal(N) * (Decimal(5) / 3).ln()
    S_max = int(X_.sqrt().to_integral_value(rounding=ROUND_CEILING))    # least S with S^2 >= 2 (6/5) N ln(5/3)
    assert Decimal(S_max - 1) ** 2 < X_ <= Decimal(S_max) ** 2
    M_max = S_max + S_max // (1 << 13) + 1
    L_cap, ST, DM, REC = 1 << 94, 1 << 40, 1 << 40, 16
    rows = {"walk_slots": SLOTC * (M_max + L_cap),
            "constants": 1 << 22,
            "input_check": 508 * F12 + (1 << 16),
            "table": Ri * (1028 * F12 + (1 << 13) + 2 * SCALAR),
            "trail_starts": ST * (1017 * F12 + SEL + 2 * SCALAR + (1 << 10)),
            "reports": DM * (1 << 15),
            "recoveries": REC * (508 * F12 + 512 * FRM + (1 << 10)),
            "output": 1 << 16, "fixed": 1 << 32, "coins": 2 * COIN}
    run_total = sum(rows.values())
    mu, loss = Decimal(6) / 5, Decimal(3) / 10000
    p_main = 1 - (-(Decimal(S_max) ** 2 / (2 * mu * Decimal(N)))).exp()   # H1 at the budget
    p_run = p_main - loss
    assert p_run >= Decimal(39) / 100
    total = Decimal(run_total)                             # every-run cap; no gate
    pre = rows["constants"] + rows["input_check"] + rows["table"]
    mem = DM * 128 + 12 * Ri * 384 + 12 * Ri * 64 + (1 << 20)
    ln2 = Decimal(2).ln()
    lg = lambda v: Decimal(v).ln() / ln2
    up = lambda v: (lg(v) * 10 ** 6).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** 6)
    a = Decimal(S_max) ** 2 / (2 * Decimal(N))
    score = up(total / C)

    def lE(L):
        ts = []
        for d in range(1, min(L // 2, Ri) + 1):
            t = (math.lgamma(Ri + 1) - math.lgamma(d + 1) - math.lgamma(Ri - d + 1)) / math.log(2)
            ts.append(t + L * math.log2(d / Ri) - d * math.log2(12))
        mx = max(ts)
        return mx + math.log2(sum(2 ** (t - mx) for t in ts))

    def lsum(v):
        mx = max(v)
        return mx + math.log2(sum(2 ** (t - mx) for t in v))
    nu_rest = lsum([lE(L) for L in range(3, 25)])
    tau = lsum([lE(L) for L in range(25, 513)])
    tail = lsum([math.log2(L / 2) - L / 2 for L in range(513, 4096)])
    return {"slot_instructions": SLOTC, "classes_N": str(N), "S_max": str(S_max), "S_max_log2": float(lg(S_max)),
            "M_max": str(M_max), "rows": {k2: str(v) for k2, v in rows.items()},
            "rows_log2_units_rounded_up": {k2: math.ceil(float(lg(Decimal(v) / C)) * 10 ** 4) / 10 ** 4 for k2, v in rows.items()},
            "rows_log2_exact": {k2: float(lg(Decimal(v) / C)) for k2, v in rows.items()},
            "run_total_instructions": str(run_total),
            "p_main_H1_mu_6_5": str(p_main)[:22], "success_lower": str(p_run)[:22],
            "log2_W_unrounded": str(lg(total / C))[:24],
            "time_log2": str(score), "score_rational": [str(int(score * 10 ** 6)), "1000000"],
            "preprocessing_log2": str(up(Decimal(pre) / C)), "peak_bytes": str(mem), "memory_log2_bytes": str(up(mem)),
            "memory_below_2^48": mem < (1 << 48),
            "S2_over_2N": float(a), "ideal_mapping_P_main": float(1 - (-a).exp()),
            "E2_log2": lE(2), "inv_12R_log2": -math.log2(12 * Ri), "sum_E3_to_E24_over_E2": 2 ** (nu_rest - lE(2)),
            "trap_sum_25_512_log2": tau, "trap_tail_513_4095_log2_bound": tail}


def stated_figures(L, k):
    """Compare every figure stated in submission.json and proof.md with the computed ledger."""
    from pathlib import Path
    here = Path(__file__).resolve().parent
    if (here / "files" / "submission.json").exists():      # verifier sandbox: snapshot under files/
        here = here / "files"
    man = json.loads((here / "submission.json").read_text())
    st = man["claims"][0]["statement"]
    pr_ = (here / "proof.md").read_text()
    sc = man["claims"][0]["score"]
    getcontext().prec = 120
    claimed = Decimal(sc["numerator"]) / Decimal(sc["denominator"])
    fmt = lambda n: f"{int(n):,}"
    need_claim = ["2^" + L["time_log2"], fmt(L["slot_instructions"]), fmt(k["f12_mul_executed"]),
                  fmt(k["selector_executed"]), str(k["glue_executed"]), "2^" + L["memory_log2_bytes"],
                  "2^" + L["preprocessing_log2"], "at least 0.39"]
    need_proof = ["**time_log2 = " + L["time_log2"] + "**", L["S_max"], "Preprocessing 2^" + L["preprocessing_log2"],
                  L["run_total_instructions"], fmt(L["slot_instructions"]),
                  str(k["glue_executed"]), L["rows"]["walk_slots"], "2^" + L["memory_log2_bytes"], "2^" + L["preprocessing_log2"]]
    for key in ("input_check", "table", "trail_starts", "recoveries"):
        need_proof.append(fmt(L["rows"][key]))
    mach = (here / "machine.md").read_text()
    need_mach = [fmt(k["f12_mul_executed"]), fmt(k["selector_executed"]), str(k["glue_executed"]),
                 fmt(L["slot_instructions"]), fmt(k["stack_bytes"]) + " bytes", str(k["fp_mul"])]
    need_proof.append(fmt(SCALAR))
    wt = math.log2((1 << 94) * (int(L["S_max"]) + (1 << 94)) / (2 * int(L["classes_N"])))
    need_proof.append("2^%.2f" % wt)
    toy = json.loads((here / "toy_offline.json").read_text())
    for fam in ("bn", "generic"):
        F_ = toy[fam]
        for R_, v in F_["runs"].items():
            sci = lambda x: ("%.2e" % x).replace("e-0", "e-")
            need_proof.append("| %s (%s) |" % (sci(v["rejections"] / v["slots"]), sci(nu_bound(int(R_), F_["classes"]))))
    gv = toy["generic"]["runs"]["64"]
    need_proof.append("%s rejections in %s slots" % (fmt(gv["rejections"]), fmt(gv["slots"])))
    for fam in ("bn", "generic"):                         # toy table recomputed at the production budget rule
        F_ = toy[fam]; sb = math.isqrt(int(2 * 6 * F_["classes"] * math.log(5 / 3) / 5)) + 1
        for R_, v in F_["runs"].items():
            n_ = v["instances"]
            fr = sum(1 for x in v["v"] if x is not None and x <= sb) / n_
            su = sum(1 for t in v["t"] if t is not None and t < sb) / n_
            need_proof.append(f"| {R_} | {n_} | {fr:.3f} | {su:.3f} | {v['wrong_outputs']} |")
    miss = ([x for x in need_claim if x not in st] + [x for x in need_proof if x not in pr_]
            + [x for x in need_mach if x not in mach])
    # every stated log2 in the ledger table must be rounded up at 4 decimals
    rows_ok = True
    for key, val in L["rows_log2_exact"].items():
        up4 = math.ceil(val * 10 ** 4) / 10 ** 4
        rows_ok &= f"| {up4:.4f} |" in pr_ or key in ("output", "fixed")
    of = math.ceil(math.log2((int(L["rows"]["output"]) + int(L["rows"]["fixed"])) / 20000) * 10 ** 4) / 10 ** 4
    rows_ok &= f"| {of:.4f} |" in pr_
    return {"claim_score_equals_ledger": claimed == Decimal(L["time_log2"]), "missing_figures": miss,
            "ledger_log2_rounded_up": rows_ok,
            "stated_figures_match": claimed == Decimal(L["time_log2"]) and not miss and rows_ok}


# ---------------------------------------------------------------- toy field and the full algorithm

def find_toy(rbits, qbits):
    rp = (1 << rbits) - ((1 << rbits) % 12) + 1
    while True:
        rp += 12
        if is_prime(rp): break
    for h in range(2, rp):
        lam = pow(h, (rp - 1) // 12, rp)
        if pow(lam, 6, rp) != 1 and pow(lam, 4, rp) != 1: break
    m = 12 * rp
    q0 = next(lam + k * rp for k in range(12) if (lam + k * rp) % 12 == 7)
    q = q0 + (((1 << qbits) - q0) // m + 1) * m
    while not is_prime(q):
        q += m
    for a0 in range(1, 1000):
        T = Tower.__new__(Tower); T.q = q
        xi = (a0, 1)
        if T.f2pow(xi, (q * q - 1) // 2) != (1, 0) and T.f2pow(xi, (q * q - 1) // 3) != (1, 0): break
    return q, rp, a0


def bn_toy(x):
    """BN-family toy: q = p(x), r' = r(x) from the BN polynomials (x odd and negative here), so q ~ r'
    as in production; xi = a0 + u with a0 the least value making xi neither a square nor a cube."""
    q = 36 * x ** 4 + 36 * x ** 3 + 24 * x ** 2 + 6 * x + 1
    rp = 36 * x ** 4 + 36 * x ** 3 + 18 * x ** 2 + 6 * x + 1
    assert is_prime(q) and is_prime(rp) and q % 12 == 7 and all(pow(q, k, rp) != 1 for k in range(1, 12))
    for a0 in range(1, 1000):
        T = Tower.__new__(Tower); T.q = q
        xi = (a0, 1)
        if T.f2pow(xi, (q * q - 1) // 2) != (1, 0) and T.f2pow(xi, (q * q - 1) // 3) != (1, 0): break
    return q, rp, a0


def mix(x):
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9 & M64
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb & M64
    return x ^ (x >> 31)


def khash(K, salt):
    h = salt
    for k in range(4):
        h = mix(h ^ ((K >> (64 * k)) & M64))
    return h


class ClassRho:
    """proof.md's algorithm with toy sizes: R indices, window w, A probes, dbits distinguished bits,
    trail cap 2^(dbits+7) slots. Starts fresh trails while S < S_max; the last trail runs to its end."""

    def __init__(self, T, rp, g, R_idx, dbits, rng, w=24, A=8):
        self.T, self.rp, self.g, self.R, self.dbits, self.rng, self.w, self.A = T, rp, g, R_idx, dbits, rng, w, A
        self.rb = R_idx.bit_length() - 1
        self.lam = T.q % rp

    def setup(self, z):
        T, rp = self.T, self.rp
        self.z, self.s1 = z, self.rng.getrandbits(64)
        li = pow(self.lam, -1, rp)
        self.TW, self.EX = [], []
        for h in range(self.R):
            u, v = self.rng.randrange(rp), self.rng.randrange(rp)
            Th = T.mul(T.pow(self.g, u), T.pow(z, v))
            self.TW.append([T.frob_k(Th, (12 - k) % 12) for k in range(12)])
            self.EX.append([(u * pow(li, k, rp) % rp, v * pow(li, k, rp) % rp) for k in range(12)])

    def run(self, budget, track):
        T, rp, lam = self.T, self.rp, self.lam
        cap = 1 << (self.dbits + 7)
        store, seen = {}, {}
        S = slots = starts = rejects = kills = ties = trail = capped = 0
        first, tstart = None, 0
        while S < budget:
            trail += 1; starts += 1; tstart = S
            al, be = self.rng.randrange(rp), self.rng.randrange(rp)
            y = T.mul(T.pow(self.g, al), T.pow(self.z, be))
            sel = T.select(y)
            if sel is None:
                ties += 1; continue
            K, ks = sel
            win, n = [K], 0
            while True:
                w1 = khash(K, self.s1)
                h, dp = w1 & (self.R - 1), ((w1 >> self.rb) & ((1 << self.dbits) - 1)) == 0
                if track:
                    c = tuple(sorted(T.keys(y)))         # the 12-key multiset: a class invariant
                    prev = seen.get(c)
                    if prev is None:
                        seen[c] = trail
                    elif prev != trail and first is None:
                        first = S + starts
                if dp:
                    f = pow(lam, ks, rp)
                    ac, bc = al * f % rp, be * f % rp
                    if K in store:
                        a2, b2 = store[K]
                        if b2 != bc:
                            x = (ac - a2) * pow((b2 - bc) % rp, -1, rp) % rp
                            if T.pow(self.g, x) == self.z:
                                return dict(x=x, start=tstart, first=first, slots=slots, rejects=rejects,
                                            kills=kills, ties=ties, starts=starts, capped=capped)
                    else:
                        store[K] = (ac, bc)
                    break
                if n >= cap:
                    capped += 1
                    break
                ok = False
                for a in range(self.A):
                    hh = (h + a) % self.R
                    y2 = T.mul(y, self.TW[hh][ks]); slots += 1; n += 1
                    s2 = T.select(y2)
                    if s2 is None:
                        ties += 1; break
                    if s2[0] in win:
                        rejects += 1; continue
                    ua, va = self.EX[hh][ks]
                    al, be, y, (K, ks) = (al + ua) % rp, (be + va) % rp, y2, s2
                    win.append(K)
                    if len(win) > self.w: win.pop(0)
                    S += 1; ok = True
                    break
                if not ok:
                    kills += 1; break
        return dict(x=None, start=tstart, first=first, slots=slots, rejects=rejects, kills=kills, ties=ties,
                    starts=starts, capped=capped)


def wilson(k, n, z=1.96):
    """Wilson 95% score interval for a binomial proportion k/n."""
    p_ = k / n; d = 1 + z * z / n
    c = (p_ + z * z / (2 * n)) / d; hw = z * math.sqrt(p_ * (1 - p_) / n + z * z / (4 * n * n)) / d
    return [round(c - hw, 4), round(c + hw, 4)]


def toy(seed, instances, R_idx, family="bn", rbits=22, qbits=40, bn_x=-15, dbits=4):
    q, rp, a0 = bn_toy(bn_x) if family == "bn" else find_toy(rbits, qbits)
    T = Tower(q, a0)
    rng = random.Random(10 ** 6)
    while True:
        y = tuple(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(3)) for _ in range(2))
        g = T.pow(y, (q ** 12 - 1) // rp)
        if g != T.one: break
    N = (rp - 1) // 12
    S_max = math.isqrt(int(2 * 6 * N * math.log(5 / 3) / 5)) + 1     # same budget rule as production
    rows = []
    for i in range(instances):
        rr = random.Random((seed << 20) + (R_idx << 8) + i)
        x = rr.randrange(1, rp)
        alg = ClassRho(T, rp, g, R_idx, dbits, rr)
        alg.setup(T.pow(g, x))
        res = alg.run(4 * S_max, track=True)
        res["correct"] = res["x"] == x if res["x"] is not None else None
        rows.append(res)
    found = [o for o in rows if o["x"] is not None]
    fv = [o["first"] for o in rows if o["first"] is not None]
    return {"family": family, "q": q, "r_prime": rp, "xi": f"{a0}+u", "classes": N, "R": R_idx, "window": 24, "probes": 8,
            "dp_bits": dbits, "trail_cap": 1 << (dbits + 7), "S_budget": S_max, "instances": instances, "seed": seed,
            "H1_requirement_P_main": 0.3903,
            "success_within_budget": sum(1 for o in found if o["start"] < S_max and o["correct"]),
            "wrong_outputs": sum(1 for o in found if not o["correct"]),
            "first_repeat_within_budget": sum(1 for v in fv if v <= S_max),
            "mu_hat": round(sum(v * v for v in fv) / max(1, len(fv)) / (2 * N), 4),
            "mu_hat_se": round(math.sqrt(sum((v * v / (2 * N)) ** 2 for v in fv) / max(1, len(fv))
                                         - (sum(v * v for v in fv) / max(1, len(fv)) / (2 * N)) ** 2) / math.sqrt(max(1, len(fv))), 4),
            "success_wilson95": wilson(sum(1 for o in found if o["start"] < S_max and o["correct"]), instances),
            "P_from_mu_hat": round(1 - math.exp(-(S_max ** 2 / (2 * N)) / max(1e-9, sum(v * v for v in fv) / max(1, len(fv)) / (2 * N))), 4),
            "reject_rate": round(sum(o["rejects"] for o in rows) / max(1, sum(o["slots"] for o in rows)), 6),
            "inv_12R": round(1 / (12 * R_idx), 6), "nu_bound": round(nu_bound(R_idx, N), 6),
            "killed_trails": sum(o["kills"] for o in rows), "ties": sum(o["ties"] for o in rows),
            "trails_at_cap": sum(o["capped"] for o in rows),
            "first_repeat_visited": fv}


def main():
    args = json.load(sys.stdin)
    mode, seed = args.get("mode", "all"), args.get("seed", 1)
    out = {}
    k = kernels(seed)
    L = ledger(k)
    out["stated_figures"] = stated_figures(L, k)          # every mode: declared figures against the computed ledger
    if mode in ("all", "kernels"):
        out["kernels"] = k
    if mode in ("all", "algebra"):
        out["algebra"] = algebra(seed)
    if mode in ("all", "ledger"):
        out["parameters"] = PARAM
        out["ledger"] = L
    if mode in ("all", "toy"):
        out["toy"] = toy(seed, args.get("toy_instances", 48), args.get("toy_R", 16), "bn")
    out["scope"] = ("kernel counts are executed instructions of the generated programs; algebra on the pinned "
                    "constants; toy runs are evidence at toy scale only; no production attack is run")
    print(json.dumps(out, indent=1))
    bad = [n for n in ("fp_mismatches", "f12_mismatches", "selector_mismatches", "glue_mismatches") if k[n]]
    if bad or not out["stated_figures"]["stated_figures_match"]:
        print("FAIL: " + ", ".join(bad + ["stated_figures"] * (not out["stated_figures"]["stated_figures_match"])), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
