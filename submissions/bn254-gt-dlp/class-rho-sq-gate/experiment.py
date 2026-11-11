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
from fractions import Fraction

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
            elif op == 'adr':                    # ADR: PC-relative address of a label (code addresses are indices)
                s(a[0], self.p.labels[a[1]])
            elif op == 'blr':                    # BLR: call through a register, writes the link register
                stack.append(pc); self.r['x30'] = pc; pc = g(a[0])
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


def gen_cyc_sqr(pr):
    """x0=out, x1=a (a in GT, out distinct from a): Granger-Scott cyclotomic squaring. With s = w^3 (s^2 = xi),
    A = (e0, e3), B = (e1, e4), C = (e2, e5) are squared in Fp4 by 3 Fp2 squarings each, then
    e0' = 3A0 - 2e0, e3' = 3A1 + 2e3, e1' = 3 xi C1 + 2e1, e4' = 3C0 - 2e4, e2' = 3B0 - 2e2, e5' = 3B1 + 2e5."""
    e = pr.emit
    pr.label('cyc_sqr')
    e('stp_pre', 'x19', 'x20', 'sp', -16); e('stp_pre', 'x21', 'x22', 'sp', -16)
    e('stp_pre', 'x23', 'x24', 'sp', -16); e('stp_pre', 'x25', 'x30', 'sp', -16)
    e('sub_sp', 4096)
    e('mov', 'x19', 'x1'); e('mov', 'x21', 'x0'); e('mov', 'x22', 'sp')
    e('add', 'x23', 'sp', 512); e('add', 'x24', 'sp', 1024)
    T = Tracer(pr, ['x22', 'x23', 'x24'])
    mul, add, sub, addu = T.mul, T.add, T.sub, T.addu

    def f2add(a, b, d=(None, None)):
        return (add(a[0], b[0], d[0]), add(a[1], b[1], d[1]))

    def f2sub(a, b, d=(None, None)):
        return (sub(a[0], b[0], d[0]), sub(a[1], b[1], d[1]))

    def free2(x):
        T.rel(x[0]); T.rel(x[1])

    def f2sq(a):                                   # (x + yu)^2 = (x + y)(x - y) + 2xy u: 2 products
        sa = addu(a[0], a[1]); da = sub(a[0], a[1]); r0 = mul(sa, da); T.rel(sa); T.rel(da)
        tw = addu(a[0], a[0]); r1 = mul(tw, a[1]); T.rel(tw)
        return (r0, r1)

    def f2xi(a, d=(None, None)):
        def nine(x):
            y = add(x, x); y = add(y, y, y); y = add(y, y, y)
            return add(y, x, y)
        n0 = nine(a[0]); n1 = nine(a[1])
        r0 = sub(n0, a[1], d[0]); r1 = add(n1, a[0], d[1])
        T.rel(n0); T.rel(n1)
        return (r0, r1)

    def sq4(a, b):                                 # (a + b s)^2 = (a^2 + xi b^2) + ((a + b)^2 - a^2 - b^2) s
        t0 = f2sq(a); t1 = f2sq(b); ab = f2add(a, b); t2 = f2sq(ab); free2(ab)
        t2 = f2sub(t2, t0, t2); r1 = f2sub(t2, t1, t2)
        x = f2xi(t1); free2(t1); r0 = f2add(t0, x, t0); free2(x)
        return r0, r1

    Y, O = tower_slots('x19'), tower_slots('x21')
    E_ = [Y[0][0], Y[1][0], Y[0][1], Y[1][1], Y[0][2], Y[1][2]]   # e0..e5
    O_ = [O[0][0], O[1][0], O[0][1], O[1][1], O[0][2], O[1][2]]

    def out(x, c, d, minus):                       # d = 3x - 2c (minus) or 3x + 2c
        t = f2sub(x, c) if minus else f2add(x, c)
        t = f2add(t, t, t); f2add(t, x, d); free2(t)

    A0, A1 = sq4(E_[0], E_[3])
    out(A0, E_[0], O_[0], True); out(A1, E_[3], O_[3], False); free2(A0); free2(A1)
    B0, B1 = sq4(E_[1], E_[4])
    out(B0, E_[2], O_[2], True); out(B1, E_[5], O_[5], False); free2(B0); free2(B1)
    C0, C1 = sq4(E_[2], E_[5])
    X = f2xi(C1); free2(C1)
    out(X, E_[1], O_[1], False); out(C0, E_[4], O_[4], True); free2(X); free2(C0)
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
    """One walk slot outside f12_mul, cyc_sqr and selector, branch-free: the same instructions run whether
    the probe squares or multiplies, is accepted, hits the window, ties or reaches a distinguished key;
    the product routine is chosen by CSEL of two ADR addresses and called by BLR; results are committed
    with CSEL. Registers: x19 = state record (K, k*, tie) with the candidate record at
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
    e('ldr', 'x8', 'x25', 80)                          # branch bit of the current class (1: squaring class)
    e('cmp', 'x21', 'xzr'); e('cset', 'x9', 'eq'); e('and_r', 'x8', 'x8', 'x9')   # squaring only on probe 0
    e('adr', 'x16', 'f12_mul'); e('adr', 'x17', 'cyc_sqr')
    e('cmp', 'x8', 'xzr'); e('csel', 'x16', 'x17', 'x16', 'ne')
    e('csel', 'x27', 'x24', 'x27', 'ne')               # exponent addend: (alpha, beta) itself, which doubles them
    e('mov', 'x0', 'x23'); e('mov', 'x1', 'x22'); e('mov', 'x2', 'x3'); e('blr', 'x16')
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
    e('ldr', 'x2', 'x25', 64); e('lsr', 'x2', 'x2', 63); e('ldr', 'x3', 'x25', 80)
    e('csel', 'x2', 'x3', 'x2', 'eq'); e('str', 'x2', 'x25', 80)   # branch bit of the new class: bit 63 of w2
    e('mov', 'x0', 'x22'); e('csel', 'x22', 'x22', 'x23', 'eq'); e('csel', 'x23', 'x23', 'x0', 'eq')
    e('add', 'x0', 'x21', 1)                           # next probe index: a+1, or 0 on accept
    e('csel', 'x21', 'x0', 'xzr', 'eq')
    e('ldr', 'x8', 'x26', 80)                          # distinguished flag of the new key, only if accepted
    e('lsr', 'x2', 'x4', 14); e('cmp', 'x2', 'xzr'); e('cset', 'x3', 'eq')
    e('ldr', 'x5', 'x25', 64); e('and_r', 'x5', 'x5', 'x8'); e('cmp', 'x5', 'xzr'); e('cset', 'x5', 'eq')
    e('and_r', 'x3', 'x3', 'x5'); e('and_r', 'x3', 'x3', 'x10'); e('str', 'x3', 'x25', 104)
    e('mov', 'x9', 'xzr')                              # cap flags
    e('ldp', 'x0', 'x1', 'x25', 0); e('ldr', 'x4', 'x25', 32); e('movi', 'x2', 1)   # run slot counter: 3 words
    e('adds', 'x0', 'x0', 'x2'); e('adcs', 'x1', 'x1', 'xzr'); e('adc', 'x4', 'x4', 'xzr')
    e('stp', 'x0', 'x1', 'x25', 0); e('str', 'x4', 'x25', 32)
    e('ldp', 'x5', 'x6', 'x26', 96); e('ldr', 'x7', 'x26', 128)                    # its 3-word cap
    e('subs', 'xzr', 'x0', 'x5'); e('sbcs', 'xzr', 'x1', 'x6'); e('sbcs', 'xzr', 'x4', 'x7')
    e('cset', 'x7', 'cs'); e('orr', 'x9', 'x9', 'x7')
    e('ldp', 'x0', 'x1', 'x25', 16)                    # trail slot counter: 2 words, cap 2^94
    e('adds', 'x0', 'x0', 'x2'); e('adc', 'x1', 'x1', 'xzr'); e('stp', 'x0', 'x1', 'x25', 16)
    e('ldp', 'x5', 'x6', 'x26', 112); e('subs', 'xzr', 'x0', 'x5'); e('sbcs', 'xzr', 'x1', 'x6')
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
    t12 = gen_f12_mul(pr); tsq = gen_cyc_sqr(pr); tsel = gen_selector(pr); gen_glue(pr)
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
    bad12 = badsel = badsq = 0
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
        mc.r['x0'], mc.r['x1'] = YO, YA
        nsq = mc.call('cyc_sqr')                          # GT input: equals the dense square
        badsq += from_mem([unmont(getv(YO + 32 * k)) for k in range(12)]) != TW.mul(y, y)
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
    M_cap, L_cap = slot_bounds()[4], 1 << 94               # production caps: run slots (3 words), trail slots (2 words)
    pool = [salt1, salt2, C1, C2, (1 << 14) - 1, 12, 384, 40, TWB, EXB, (1 << 37) - 1, 0,
            M_cap & M64, (M_cap >> 64) & M64, L_cap & M64, L_cap >> 64, M_cap >> 128]
    assert M_cap >> 192 == 0 and M_cap >> 128 > 0
    # counter values before the slot (run slots, trail slots): carries into each word and both sides of each cap
    presets = [(((M_cap >> 128) + 1) << 128, rng.getrandbits(40)), (M64, rng.getrandbits(40)), ((1 << 128) - 1, M64),
               (M_cap - 2, L_cap - 2), (M_cap - 1, rng.getrandbits(40)), (rng.getrandbits(130), L_cap - 1),
               ((M_cap & ((1 << 128) - 1)) + 5, (1 << 64) + 5), (M_cap + 7, L_cap + 7)]
    for k_, w_ in enumerate(pool):
        mc.st(CONST + GLUE_POOL + 8 * k_, w_)
    glue_bad, counts, cap_cases = 0, set(), 0
    for trial in range(8):
        y = TW.pow(g, rng.randrange(1, R)); T_ = TW.pow(g, rng.randrange(1, R))
        V0 = [mont(v) for v in TW.keys(y)]; Kc = min(V0); kc = V0.index(Kc)
        sqp = (trial // 2) % 2 == 1                       # trials 2, 3, 6, 7 take the squaring path
        sqbit = 1 if sqp else rng.randrange(2)
        h, a_ = rng.randrange(1 << 14), (0 if sqp else (rng.randrange(1, 7) if sqbit else rng.randrange(7)))
        t_ = ((h + a_) & ((1 << 14) - 1)) * 12 + kc
        for k, v in enumerate(mem_order(y)): setc(Y + 32 * k, mont(v))
        for k, v in enumerate(mem_order(T_)): setc(TWB + 384 * t_ + 32 * k, mont(v))
        u_, v_, al, be = (rng.randrange(R) for _ in range(4))
        setc(EXB + 64 * t_, u_); setc(EXB + 64 * t_ + 32, v_); setc(AB, al); setc(AB + 32, be)
        setc(REC, Kc); mc.st(REC + 32, kc); mc.st(REC + 40, 0)
        nxt = TW.mul(y, y) if sqp else TW.mul(y, T_)
        Vn = [mont(v) for v in TW.keys(nxt)]; Kn = min(Vn); kn = Vn.index(Kn)
        w1, w2 = khash(Kn, salt1), khash(Kn, salt2)
        hit = trial % 2 == 1                              # odd trials plant a window hit
        win = [rng.getrandbits(64) for _ in range(24)]
        if hit:
            win[rng.randrange(24)] = w1
        for j in range(24):
            mc.st(WIN + 8 * j, win[j]); setc(WIN + 256 + 32 * j, rng.getrandbits(250))
        cs_, ct_ = presets[trial]
        mc.st(CNT, cs_ & M64); mc.st(CNT + 8, (cs_ >> 64) & M64); mc.st(CNT + 32, cs_ >> 128)
        mc.st(CNT + 16, ct_ & M64); mc.st(CNT + 24, ct_ >> 64)
        ring = rng.randrange(24)
        mc.st(CNT + 48, h); mc.st(CNT + 56, ring); mc.st(CNT + 80, sqbit)
        for reg, val in (('x19', REC), ('x20', WIN), ('x21', a_), ('x22', Y), ('x23', NX), ('x24', AB), ('x25', CNT)):
            mc.r[reg] = val
        mc.r['x30'] = 0x6c756500                          # glue entry address: RET must return here
        counts.add(mc.call('glue') - (nsq if sqp else n12) - nsel)
        glue_bad += mc.r['x30'] != 0x6c756500 or mc.r['sp'] != 1 << 40   # next slot starts at glue's first instruction
        acc = not hit
        dp = (w1 >> 14) == 0 and (w2 & ((1 << 37) - 1)) == 0 and acc
        slots_after = mc.ld(CNT) + (mc.ld(CNT + 8) << 64) + (mc.ld(CNT + 32) << 128)
        trail_after = mc.ld(CNT + 16) + (mc.ld(CNT + 24) << 64)
        exit_flag = int(cs_ + 1 >= M_cap or ct_ + 1 >= L_cap or dp)
        ok = (mc.ld(CNT + 96) == int(acc) and mc.ld(CNT + 104) == int(dp)
              and slots_after == cs_ + 1 and trail_after == ct_ + 1 and mc.ld(CNT + 112) == exit_flag
              and mc.r['x21'] == (0 if acc else a_ + 1))
        cap_cases += exit_flag
        da, db = (al, be) if sqp else (u_, v_)            # squaring doubles the exponents
        if acc:
            ok &= (from_mem([unmont(getv(NX + 32 * k)) for k in range(12)]) == nxt and getv(REC) == Kn
                   and mc.ld(REC + 32) == kn and mc.ld(CNT + 48) == w1 & ((1 << 14) - 1)
                   and mc.ld(CNT + 80) == w2 >> 63
                   and mc.ld(WIN + 8 * ring) == w1 and getv(WIN + 256 + 32 * ring) == Kn
                   and mc.ld(CNT + 56) == (ring + 1) % 24 and getv(AB) == (al + da) % R
                   and getv(AB + 32) == (be + db) % R and mc.r['x22'] == NX and mc.r['x23'] == Y)
        else:
            ok &= (getv(REC) == Kc and mc.ld(REC + 32) == kc and mc.ld(CNT + 48) == h and mc.ld(CNT + 80) == sqbit
                   and mc.ld(WIN + 8 * ring) == win[ring] and mc.ld(CNT + 56) == ring
                   and getv(AB) == al and getv(AB + 32) == be and mc.r['x22'] == Y and mc.r['x23'] == NX)
        glue_bad += not ok
    nglue = max(counts)
    glue_bad += len(counts) != 1                          # the count must not depend on the path
    glue_bad += cap_cases != 4                            # trials 0, 4, 5 and 7 reach a cap and must set the exit flag
    regs = set()
    for op, a in pr.code:
        regs |= {x for x in a if isinstance(x, str) and x.startswith('x') and x[1:].isdigit()}
    return {"fp_mul": cnt['fp_mul'], "fp_add": cnt['fp_add'], "fp_sub": cnt['fp_sub'],
            "fr_mul": cnt['fr_mul'], "fr_add": cnt['fr_add'], "fp_mismatches": bad,
            "f12_mul_executed": n12, "f12_calls": t12.calls, "f12_scratch_slots": t12.peak, "f12_mismatches": bad12,
            "cyc_sqr_executed": nsq, "cyc_sqr_calls": tsq.calls, "cyc_sqr_scratch_slots": tsq.peak, "cyc_sqr_mismatches": badsq,
            "selector_executed": nsel, "selector_calls": tsel.calls, "selector_mismatches": badsel,
            "stack_bytes": mc.depth, "registers_used": sorted(regs, key=lambda s: int(s[1:])),
            "glue_executed": nglue, "glue_mismatches": glue_bad, "glue_paths_tested": "multiply and squaring paths, each 2 accepted and 2 window hits; counter carries and both sides of both caps",
            "slot_cap_words": 3,
            "program_instructions": len(pr.code)}


R1_FACTORS = ((2, 28), (3, 2), (13, 1), (29, 1), (983, 1), (11003, 1), (237073, 1), (405928799, 1),
              (1670836401704629, 1), (13818364434197438864469338081, 1))


def alg_t_min_log2():
    """log2 of the least t >= 1 with 2^t in the order-12 group generated by lambda mod r (4 decimals, rounded down)."""
    o2 = R - 1
    for q_, e_ in R1_FACTORS:
        for _ in range(e_):
            if pow(2, o2 // q_, R) == 1:
                o2 //= q_
            else:
                break
    return math.floor(math.log2(o2 // math.gcd(o2, 12)) * 10 ** 4) / 10 ** 4


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
    # squaring branch: the class of y^2 depends only on the class of y (squaring commutes with sigma)
    sq_ok = 0
    for _ in range(4):
        y = TW.pow(g, rng.randrange(1, R)); j = rng.randrange(1, 12)
        y2 = TW.frob_k(y, j)
        s1, s2 = TW.select(TW.mul(y, y)), TW.select(TW.mul(y2, y2))
        sq_ok += (TW.frob(TW.mul(y, y)) == TW.mul(TW.frob(y), TW.frob(y)) and s1 is not None and s2 is not None
                  and s1[0] == s2[0] and TW.frob_k(TW.mul(y, y), s1[1]) == TW.frob_k(TW.mul(y2, y2), s2[1]))
    # order of 2 mod r: 2^t lies in <lambda> (order 12) only if ord(2) divides 12 t
    prod = 1
    for q_, e_ in R1_FACTORS:
        prod *= q_ ** e_
    o2 = R - 1
    for q_, e_ in R1_FACTORS:
        for _ in range(e_):
            if pow(2, o2 // q_, R) == 1:
                o2 //= q_
            else:
                break
    t_min = o2 // math.gcd(o2, 12)
    T = TW.pow(g, 12345)
    out.update({"squaring_class_function": f"{sq_ok}/4",
                "r_minus_1_factorisation": prod == R - 1 and all(is_prime(q_) for q_, _ in R1_FACTORS),
                "order_of_2_mod_r_log2": round(math.log2(o2), 4), "order_2_is_order": pow(2, o2, R) == 1,
                "least_t_with_2^t_in_lambda_group_log2": round(math.log2(t_min), 4),
                "no_squaring_return_below_2^94_steps": t_min > 1 << 94,
                "keys_match_explicit_conjugates": f"{keys_ok}/{trials}", "keys_shift_under_frobenius": f"{shift_ok}/{trials}",
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
         "ST_max_log2": 50, "D_max_log2": 50, "REC_max": 16, "squaring": "1/2 (bit 63 of w2; probe 0 only)",
         "mean_allowance": "11/10", "slot_cap": "2^8 x ceil((11/10) sqrt(pi N'/2)), N' = 4RN/(3R-1)",
         "gate": "one fresh 64-bit coin before any work; the run opens when it is below K"}
C_MEAN = Fraction(11, 10)   # H1: an opened run uses at most C_MEAN sqrt(pi N'/2) slots in expectation
CAP_LOG2 = 8                # run slot cap M_cap = 2^8 times that bound (Markov: reached with probability <= 2^-8)
GATE = 64                   # the gate: one coin word (16), compare with K, branch, return
REJ_LOG2 = 13               # allowance for multiplying slots after a rejected probe: 2^-13 of the slots
REVISIT = 1 << 96           # allowance for slots at classes visited before, all priced as multiplying slots
LO_RUN = 1 - Fraction(1, 1 << 8) - Fraction(1, 1 << 10) - Fraction(1, 1 << 20)   # success of an opened run


def dec_pi():
    """pi to the current Decimal precision (Machin: 16 atan(1/5) - 4 atan(1/239))."""
    def atan_inv(n):
        x = Decimal(1) / n; t, s_, k = x, x, 1
        while abs(t) > Decimal(10) ** -(getcontext().prec + 5):
            t = -t / (n * n); k += 2; s_ += t / k
        return s_
    return 16 * atan_inv(5) - 4 * atan_inv(239)


def slot_bounds():
    """N, N' (exact and decimal), E_slots = ceil(C_MEAN sqrt(pi N'/2)) and the run slot cap M_cap = 2^8 E_slots."""
    getcontext().prec = 120
    N = (R - 1) // 12
    Ne = n_eff(N, PARAM["R_index"], True)                  # N' = 4RN/(3R-1)
    NeD = Decimal(Ne.numerator) / Decimal(Ne.denominator)
    X_ = Decimal(C_MEAN.numerator) / C_MEAN.denominator * (dec_pi() * NeD / 2).sqrt()
    E_slots = int(X_.to_integral_value(rounding=ROUND_CEILING))
    assert Decimal(E_slots - 1) < X_ <= Decimal(E_slots)
    return N, Ne, NeD, E_slots, E_slots << CAP_LOG2


ROWS = ("gate", "walk_slots", "multiply_surcharge", "input_check", "constants", "table", "trail_starts", "reports",
        "recoveries", "output_fixed", "coins")


def ledger(k):
    getcontext().prec = 120
    C = 20000
    F12, SQR, SEL, FRM = k["f12_mul_executed"], k["cyc_sqr_executed"], k["selector_executed"], k["fr_mul"] + 4
    SLOTM = F12 + SEL + k["glue_executed"]                 # multiplying slot
    SLOTS = SQR + SEL + k["glue_executed"]                 # squaring slot
    Ri = PARAM["R_index"]
    N, Ne, NeD, E_slots, M_cap = slot_bounds()
    # multiplying slots of an opened run, in expectation: half of the slots at new classes, the slots after a
    # rejected probe (2^-13 of the slots), and every slot at a class visited before (REVISIT)
    E_mult = (E_slots + 1) // 2 + (E_slots >> REJ_LOG2) + 1 + REVISIT
    L_cap, ST, DM, REC = 1 << 94, 1 << 50, 1 << 50, 16
    other = {"input_check": 508 * F12 + (1 << 16),
             "constants": 1 << 22,
             "table": Ri * (1028 * F12 + (1 << 13) + 2 * SCALAR),
             "trail_starts": ST * (1017 * F12 + SEL + 2 * SCALAR + (1 << 10)),
             "reports": DM * (1 << 15),
             "recoveries": REC * (508 * F12 + 512 * FRM + (1 << 10)),
             "output_fixed": (1 << 16) + (1 << 32), "coins": 2 * COIN}
    exp_rows = dict(other, walk_slots=SLOTS * E_slots, multiply_surcharge=(F12 - SQR) * E_mult)
    cap_rows = dict(other, walk_slots=SLOTS * M_cap, multiply_surcharge=(F12 - SQR) * M_cap)
    E_run, cap_run = sum(exp_rows.values()), sum(cap_rows.values())
    lim = Fraction(39, 100) * (1 << 64) / LO_RUN
    K = lim.numerator // lim.denominator + (lim.numerator % lim.denominator != 0)   # least K with (K/2^64) LO_RUN >= 0.39
    assert Fraction(K, 1 << 64) * LO_RUN >= Fraction(39, 100) > Fraction(K - 1, 1 << 64) * LO_RUN and K < 1 << 64
    exp_total = GATE + Fraction(K * E_run, 1 << 64)        # expected total work: gate + P(open) x E_run
    every_run = GATE + cap_run                             # every-run bound (an opened run that reaches its caps)
    # failure budget of an opened run: slot cap 2^-8 (Markov, H1); trail-start and report caps (Markov, H2);
    # everything else below 2^-20
    rep_mean = Fraction(E_slots, 1 << 87) + 8              # expected reports, and expected trail starts
    assert rep_mean / ST + rep_mean / DM <= Fraction(1, 1 << 10)
    assert LO_RUN == 1 - Fraction(1, 1 << CAP_LOG2) - Fraction(1, 1 << 10) - Fraction(1, 1 << 20)
    pre = other["constants"] + other["input_check"] + other["table"]
    mem = DM * 128 + 12 * Ri * 384 + 12 * Ri * 64 + (1 << 20)
    ln2 = Decimal(2).ln()
    lg = lambda v: (Decimal(v.numerator) / Decimal(v.denominator) if isinstance(v, Fraction) else Decimal(v)).ln() / ln2
    up = lambda v: (lg(v) * 10 ** 6).to_integral_value(rounding=ROUND_CEILING) / Decimal(10 ** 6)
    up4 = lambda v: "%.4f" % (math.ceil(float(lg(v) - lg(C)) * 10 ** 4) / 10 ** 4)
    score = up(exp_total / C)
    sc = Fraction(int(score * 10 ** 6), 10 ** 6)
    exp_rows["gate"] = cap_rows["gate"] = GATE
    fmt = lambda n: f"{int(n):,}"
    table = [f"| {key} | {fmt(exp_rows[key])} | {up4(exp_rows[key])} | {fmt(cap_rows[key])} | {up4(cap_rows[key])} |" for key in ROWS]
    table.append(f"| **opened run** | **{fmt(GATE + E_run)}** | {up4(GATE + E_run)} | **{fmt(every_run)}** | {up4(every_run)} |")

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
    return {"slot_instructions_multiply": SLOTM, "slot_instructions_square": SLOTS,
            "classes_N": str(N), "effective_classes_N_prime": "%d/%d" % (Ne.numerator, Ne.denominator),
            "E_slots": str(E_slots), "E_slots_log2": float(lg(E_slots)), "M_cap": str(M_cap), "M_cap_log2": float(lg(M_cap)),
            "E_mult": str(E_mult),
            "expectation_rows": {k2: str(v) for k2, v in exp_rows.items()},
            "every_run_rows": {k2: str(v) for k2, v in cap_rows.items()},
            "table": table,
            "opened_run_expected_instructions": str(GATE + E_run), "every_run_instructions": str(every_run),
            "every_run_log2": str(up(Fraction(every_run, C))),
            "opened_run_expected_log2": str(up(Fraction(GATE + E_run, C))),
            "gate_K": str(K), "gate_open_probability": float(Fraction(K, 1 << 64)),
            "opened_run_success_lower": "%d/%d" % (LO_RUN.numerator, LO_RUN.denominator),
            "opened_run_success_lower_decimal": float(LO_RUN),
            "success_lower": str(Decimal(K) / Decimal(1 << 64) * Decimal(LO_RUN.numerator) / Decimal(LO_RUN.denominator))[:22],
            "expected_reports_log2": float(lg(rep_mean)),
            "cap_failure_starts_reports_log2": float(lg(rep_mean / ST + rep_mean / DM)),
            "expected_total_instructions_rounded_up": str(exp_total.numerator // exp_total.denominator + 1),
            "log2_W_unrounded": str(lg(exp_total / C))[:24],
            "time_log2": str(score), "score_rational": [str(sc.numerator), str(sc.denominator)],
            "preprocessing_log2": str(up(Fraction(pre, C))), "peak_bytes": str(mem), "memory_log2_bytes": str(up(mem)),
            "memory_below_2^64": mem < (1 << 64),
            "E2_log2": lE(2), "inv_12R_log2": -math.log2(12 * Ri), "sum_E3_to_E24_over_E2": 2 ** (nu_rest - lE(2)),
            "trap_sum_25_512_log2": tau, "trap_tail_513_4095_log2_bound": tail}


def toy_stats(toy):
    """Per-row and pooled statistics of toy_offline.json: slots at the stop over the model mean sqrt(pi N'/2)."""
    rows, pool = [], []
    for fam in ("bn", "generic"):
        F_ = toy[fam]
        for R_, v in F_["runs"].items():
            m0 = math.sqrt(math.pi * float(n_eff(F_["classes"], int(R_), True)) / 2)
            x = [s_ / m0 for s_ in v["slots"]]
            n_ = len(x); mean = sum(x) / n_
            se = math.sqrt(sum((t - mean) ** 2 for t in x) / (n_ - 1) / n_)
            rows.append((fam, int(R_), n_, m0, mean, se, max(x), v))
            pool += x
    n_ = len(pool); mean = sum(pool) / n_
    se = math.sqrt(sum((t - mean) ** 2 for t in pool) / (n_ - 1) / n_)
    return rows, (n_, mean, se, max(pool))


def toy_rows(toy):
    """The proof.md toy table, recomputed from toy_offline.json."""
    sci = lambda x: ("%.2e" % x).replace("e-0", "e-")
    rows, (n_, mean, se, mx) = toy_stats(toy)
    out = []
    for fam, Ri, n1, m0, m1, s1, x1, v in rows:
        tot = sum(v["slots"])
        out.append(f"| {fam} | {Ri} | {n1} | {m0:,.0f} | {m1:.4f} | {s1:.4f} | {x1:.2f} | {v['correct']} | {v['wrong']} | "
                   f"{v['at_slot_cap']} | {v['squaring_steps'] / tot:.4f} | "
                   f"{sci(v['rejections'] / tot)} ({sci(nu_bound(Ri, toy[fam]['classes']))}) | {v['trail_starts'] / tot:.4f} |")
    out.append(f"Pooled over {n_:,} runs: mean {mean:.4f}, standard error {se:.4f}, largest {mx:.2f}.")
    return out


def stated_figures(L, k):
    """Compare every figure stated in submission.json, proof.md and machine.md with the computed ledger."""
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
    need_claim = ["2^" + L["time_log2"], "2^" + L["every_run_log2"], "K = " + L["gate_K"],
                  fmt(L["slot_instructions_multiply"]), fmt(L["slot_instructions_square"]),
                  fmt(k["f12_mul_executed"]), fmt(k["cyc_sqr_executed"]), fmt(k["selector_executed"]),
                  str(k["glue_executed"]), "2^" + L["memory_log2_bytes"], "2^" + L["preprocessing_log2"], "at least 0.39"]
    need_proof = ["**time_log2 = " + L["time_log2"] + "**", "(exact rational %s/%s)" % tuple(L["score_rational"]),
                  "K = " + L["gate_K"], "E_slots = " + L["E_slots"], "M_cap = " + L["M_cap"],
                  "(2^%.4f)" % L["E_slots_log2"], "(2^%.4f)" % L["M_cap_log2"],
                  "2^" + L["every_run_log2"], "2^" + L["opened_run_expected_log2"],
                  "Preprocessing 2^" + L["preprocessing_log2"], L["opened_run_success_lower"],
                  fmt(L["slot_instructions_multiply"]), fmt(L["slot_instructions_square"]),
                  fmt(k["cyc_sqr_executed"]), str(k["glue_executed"]),
                  "2^" + L["memory_log2_bytes"], "2^" + L["preprocessing_log2"],
                  "2^%.4f" % alg_t_min_log2(), fmt(SCALAR), "%.6f" % (math.ceil(L["gate_open_probability"] * 10 ** 6) / 10 ** 6),
                  "2^%.2f" % L["expected_reports_log2"], "2^%.2f" % L["cap_failure_starts_reports_log2"]]
    need_proof += L["table"]
    mach = (here / "machine.md").read_text()
    need_mach = [fmt(k["f12_mul_executed"]), fmt(k["cyc_sqr_executed"]), fmt(k["selector_executed"]), str(k["glue_executed"]),
                 fmt(L["slot_instructions_multiply"]), fmt(L["slot_instructions_square"]), fmt(k["stack_bytes"]) + " bytes",
                 str(k["fp_mul"]), "%d fp_mul calls, %d inline reduced additions or subtractions and %d inline unreduced"
                 % (k["cyc_sqr_calls"]["M"], k["cyc_sqr_calls"]["A"], k["cyc_sqr_calls"]["U"]),
                 fmt(L["peak_bytes"]) + " = 2^" + L["memory_log2_bytes"]]
    toy = json.loads((here / "toy_offline.json").read_text())
    need_proof += toy_rows(toy)
    n_, mean, se, mx = toy_stats(toy)[1]
    h1_ok = mean + 5 * se <= float(C_MEAN)                 # H1 support: pooled toy mean plus 5 standard errors
    g_ = k["glue_executed"]
    need_mach += ["glue: %d executed" % g_, "the same %d\ninstructions" % g_, "+ %d = **%s**" % (g_, fmt(L["slot_instructions_multiply"])),
                  "+ %d = **%s**" % (g_, fmt(L["slot_instructions_square"]))]
    miss = ([x for x in need_claim if x not in st] + [x for x in need_proof if x not in pr_]
            + [x for x in need_mach if x not in mach])
    # every 6-decimal power of two and every fraction over 2^20 in the texts must be a computed figure
    import re
    allowed = {L["time_log2"], L["every_run_log2"], L["opened_run_expected_log2"], L["memory_log2_bytes"], L["preprocessing_log2"]}
    for name, txt in (("claim", st), ("proof.md", pr_), ("machine.md", mach), ("README.md", (here / "README.md").read_text())):
        miss += ["%s: 2^%s" % (name, t) for t in re.findall(r"2\^(\d+\.\d{6})", txt) if t not in allowed]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"\d+/1048576", txt) if t != L["opened_run_success_lower"]]
        miss += ["%s: %s" % (name, t) for t in re.findall(r"time_log2 = (\d+\.\d+)", txt) if t != L["time_log2"]]
    return {"claim_score_equals_ledger": claimed == Decimal(L["time_log2"]), "missing_figures": miss,
            "toy_mean_plus_5se_within_allowance": h1_ok,
            "stated_figures_match": claimed == Decimal(L["time_log2"]) and not miss and h1_ok}


# ---------------------------------------------------------------- hashing and budget rule (toy code in toyrho.py)

def mix(x):
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9 & M64
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb & M64
    return x ^ (x >> 31)


def khash(K, salt):
    h = salt
    for k in range(4):
        h = mix(h ^ ((K >> (64 * k)) & M64))
    return h


def n_eff(N, R_idx, sq=True):
    """Effective class count of the walk: N / (1 - s^2 - (1 - s)^2 / R) with s = 1/2 (squaring branch);
    the plain walk (sq False) is budgeted on N itself, as before."""
    return Fraction(4 * R_idx * N, 3 * R_idx - 1) if sq else Fraction(N)


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
        from pathlib import Path                          # toyrho.py sits beside this file (sandbox: under files/)
        here = Path(__file__).resolve().parent
        sys.path.insert(0, str(here / "files" if (here / "files" / "toyrho.py").exists() else here))
        from toyrho import toy
        out["toy"] = toy(seed, args.get("toy_instances", 48), args.get("toy_R", 16), "bn")
    out["scope"] = ("kernel counts are executed instructions of the generated programs; algebra on the pinned "
                    "constants; toy runs are evidence at toy scale only; no production attack is run")
    print(json.dumps(out, indent=1))
    bad = [n for n in ("fp_mismatches", "f12_mismatches", "cyc_sqr_mismatches", "selector_mismatches", "glue_mismatches") if k[n]]
    if bad or not out["stated_figures"]["stated_figures_match"]:
        print("FAIL: " + ", ".join(bad + ["stated_figures"] * (not out["stated_figures"]["stated_figures_match"])), file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
