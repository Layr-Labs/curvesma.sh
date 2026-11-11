"""Priced AArch64 programs and the instruction interpreter (imported by experiment.py)."""

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
RET_SENTINEL = 0x6c696e6b            # the caller's return address, as seen by a called routine


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
        if op == 'movi':
            assert 0 <= a[1] < 65536, a                      # MOVZ
        if op in ('lsl_i', 'lsr'):
            assert 0 <= a[2] <= 63, a                        # LSL/LSR (immediate) shift 0..63
        if op in ('eor_lsr', 'add_lsl'):
            assert 0 <= a[3] <= 63, a                        # shifted-register operand, shift 0..63
        if op == 'adds_i':
            assert 0 <= a[2] <= 4095, a                      # ADDS (immediate), 12 bits
        if op == 'extr':
            assert 0 <= a[3] <= 63, a                        # EXTR lsb 0..63
        if op == 'and_mask':
            assert 1 <= a[2] <= 63, a                        # AND with a bitmask immediate 2^w - 1, w in 1..63
        if op == 'subs_i':
            assert 0 <= a[2] <= 4095, a                      # SUBS (immediate), 12 bits
        if op == 'b':
            assert isinstance(a[0], str), a                  # B label
        if op == 'bc':
            assert a[0] in ('eq', 'ne', 'cs', 'cc') and isinstance(a[1], str), a   # B.cond label
        if op in ('cbz', 'cbnz_l'):
            assert isinstance(a[1], str), a                  # CBZ / CBNZ register, label
        if op in ('csel', 'cset'):
            assert a[-1] in ('eq', 'ne', 'cs', 'cc'), a      # the interpreter's conditions
        if op == 'bl':
            assert isinstance(a[0], str), a                  # BL label
        if op == 'movz_lsl':
            assert 0 <= a[1] < 65536 and a[2] in (0, 16, 32, 48), a   # MOVZ imm16, LSL #0/16/32/48
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
        self.on_stlr = None                  # transport hook: the completing store runs the receiver's handler
        self.sp0 = 1 << 40

    def g(self, a):
        return 0 if a == 'xzr' else self.r[a]

    def s(self, a, v):
        if a != 'xzr':
            self.r[a] = v & M64

    def ld(self, ad):
        assert ad % 8 == 0
        return self.mem.get(ad, 0)

    def ld_acq(self, ad):
        return self.ld(ad)

    def st(self, ad, v):
        assert ad % 8 == 0
        self.mem[ad] = v & M64

    def cond(self, c):
        return {'cs': self.C == 1, 'cc': self.C == 0, 'eq': self.Z == 1, 'ne': self.Z == 0}[c]

    def call(self, label):
        """Runs from `label` to its final RET. The caller's BL is modelled by a sentinel link value: the final
        RET must find it in x30 and sp unchanged; BL writes x30 and every inner RET must return to it."""
        start, pc, stack = self.count, self.p.labels[label], []
        self.r['x30'] = RET_SENTINEL
        sp_in = self.r['sp']
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
                self.depth = max(self.depth, self.sp0 - self.r['sp'])
            elif op == 'ldp_post':
                b = self.r['sp']; s(a[0], self.ld(b)); s(a[1], self.ld(b + 8)); self.r['sp'] += a[3]
            elif op == 'sub_sp':
                self.r['sp'] -= a[0]; self.depth = max(self.depth, self.sp0 - self.r['sp'])
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
            elif op == 'movz_lsl':
                s(a[0], a[1] << a[2])
            elif op == 'add':
                s(a[0], g(a[1]) + a[2])
            elif op == 'extr':
                s(a[0], ((g(a[1]) << 64 | g(a[2])) >> a[3]))
            elif op == 'and_mask':
                s(a[0], g(a[1]) & ((1 << a[2]) - 1))
            elif op == 'stlr':
                self.st(g(a[1]), g(a[0]))
                if self.on_stlr is not None:     # the store retires only after the receiver's handler returns
                    self.on_stlr(self, g(a[1]))
            elif op in ('ldar', 'ldaxr'):
                s(a[0], self.ld_acq(g(a[1])))
            elif op == 'wfe':
                pass
            elif op == 'b':
                pc = self.p.labels[a[0]]
            elif op == 'bc':
                if self.cond(a[0]):
                    pc = self.p.labels[a[1]]
            elif op == 'cbz':
                if g(a[0]) == 0:
                    pc = self.p.labels[a[1]]
            elif op == 'cbnz_l':
                if g(a[0]) != 0:
                    pc = self.p.labels[a[1]]
            elif op == 'subs_i':
                t = g(a[1]) - a[2]; s(a[0], t); self.C, self.Z = int(t >= 0), int(t & M64 == 0)
            elif op == 'stadd':
                self.st(g(a[1]), self.ld(g(a[1])) + g(a[0]))
            elif op == 'lsl_i':
                s(a[0], g(a[1]) << a[2])
            elif op == 'add_lsl':
                s(a[0], g(a[1]) + (g(a[2]) << a[3]))
            elif op == 'mov32':
                s(a[0], g(a[1]) & 0xffffffff)
            elif op == 'adds_i':
                t = g(a[1]) + a[2]; s(a[0], t); self.C, self.Z = t >> 64, int(t & M64 == 0)
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
            elif op == 'cmp':
                t = g(a[0]) - g(a[1]); self.C, self.Z = int(t >= 0), int(t & M64 == 0)
            elif op == 'bl':
                stack.append(pc); self.r['x30'] = pc; pc = self.p.labels[a[0]]     # BL writes the link register
            elif op == 'ret':
                if not stack:
                    assert self.r['x30'] == RET_SENTINEL and self.r['sp'] == sp_in, ('bad return', label)
                    return self.count - start
                assert self.r['x30'] == stack[-1], ('RET target is not the return address', label, pc)
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
        assert (max(ba, bb) + 1) * P < 1 << 256                # accumulator <= 3p fits the 5-word bound
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



def gen_emit(pr, S1=28, A=36, B=28, C=56):
    """Walker: trace t in x0..x3. State block at x29 (words): [0] s1 s2, [16] C1 C2, [32] j_lo j_hi,
    [48] BLK stride, [64] CUR listbit, [80] 40 K. Writes record (t0, t1, d2|list, j_lo, j_hi) into relay d1's
    block at +HDR + 40c; on the K-th record of a block it tail-calls callsend_w."""
    e = pr.emit
    pr.label('emit')
    e('ldp', 'x9', 'x10', 'x29', 0); e('ldp', 'x11', 'x12', 'x29', 16)
    for hreg, salt in (('x4', 'x9'), ('x5', 'x10')):
        e('mov', hreg, salt)
        for limb in ('x0', 'x1', 'x2', 'x3'):
            e('eor', hreg, hreg, limb); e('eor_lsr', hreg, hreg, hreg, 30); e('mul', hreg, hreg, 'x11')
            e('eor_lsr', hreg, hreg, hreg, 27); e('mul', hreg, hreg, 'x12'); e('eor_lsr', hreg, hreg, hreg, 31)
    e('lsr', 'x6', 'x4', S1)                                   # d1 (production: top 36 bits of w1)
    e('lsl_i', 'x7', 'x4', A); e('lsr', 'x7', 'x7', B); e('lsr', 'x8', 'x5', C); e('orr', 'x7', 'x7', 'x8')  # d2
    e('ldp', 'x13', 'x14', 'x29', 64); e('orr', 'x7', 'x7', 'x14')            # d2 | list bit
    e('ldp', 'x9', 'x10', 'x29', 32); e('adds_i', 'x11', 'x9', 1); e('adc', 'x12', 'x10', 'xzr')
    e('stp', 'x11', 'x12', 'x29', 32)                          # index j used, j+1 stored
    e('ldp', 'x15', 'x16', 'x29', 48); e('madd', 'x15', 'x6', 'x16', 'x15')     # block = BLK + stride d1
    e('add_lsl', 'x13', 'x13', 'x6', 3); e('ldr', 'x17', 'x13', 0)            # cursor word of d1
    e('ldp', 'x11', 'x12', 'x29', 80); e('madd', 'x15', 'x17', 'x11', 'x15')    # slot = block + 40 c
    e('stp', 'x0', 'x1', 'x15', HDR); e('stp', 'x7', 'x9', 'x15', HDR + 16); e('str', 'x10', 'x15', HDR + 32)
    e('add', 'x17', 'x17', 1); e('str', 'x17', 'x13', 0)
    e('cmp', 'x17', 'x12'); e('bc', 'eq', 'callsend_w')        # full block: tail call, returns to the caller
    e('ret')


def gen_radix(pr, word, shift):
    """One LSD pass on a 32-bit digit of key word `word` (0 or 1) at bit `shift` (0 or 32), over the records
    [x1, x24). count_<w><s>: x19 = counters (2^32 x 8 B). scatter_<w><s>: x19 = cursors, x20 = destination,
    x22 = 40."""
    e = pr.emit
    tag = '%d%d' % (word, shift)
    pr.label('count' + tag)
    e('ldr', 'x2', 'x1', 8 * word)
    e('lsr', 'x2', 'x2', 32) if shift else e('mov32', 'x2', 'x2')
    e('add_lsl', 'x3', 'x19', 'x2', 3); e('ldr', 'x4', 'x3', 0); e('add', 'x4', 'x4', 1); e('str', 'x4', 'x3', 0)
    e('add', 'x1', 'x1', 40); e('cmp', 'x1', 'x24'); e('bc', 'ne', 'count' + tag)
    e('ret')
    pr.label('scatter' + tag)
    e('ldp', 'x5', 'x6', 'x1', 0); e('ldp', 'x7', 'x8', 'x1', 16); e('ldr', 'x9', 'x1', 32)
    e('mov', 'x2', 'x5' if word == 0 else 'x6')
    e('lsr', 'x2', 'x2', 32) if shift else e('mov32', 'x2', 'x2')
    e('add_lsl', 'x3', 'x19', 'x2', 3); e('ldr', 'x4', 'x3', 0); e('add', 'x10', 'x4', 1); e('str', 'x10', 'x3', 0)
    e('madd', 'x11', 'x4', 'x22', 'x20')
    e('stp', 'x5', 'x6', 'x11', 0); e('stp', 'x7', 'x8', 'x11', 16); e('str', 'x9', 'x11', 32)
    e('add', 'x1', 'x1', 40); e('cmp', 'x1', 'x24'); e('bc', 'ne', 'scatter' + tag)
    e('ret')


def gen_table_ops(pr):
    """zero: clear [x19, x24) (a multiple of 32 bytes). prefix: exclusive prefix sum over the counter words
    [x19, x24), running total in x3 (set to 0 by the caller)."""
    e = pr.emit
    pr.label('zero')
    e('stp', 'xzr', 'xzr', 'x19', 0); e('stp', 'xzr', 'xzr', 'x19', 16)
    e('add', 'x19', 'x19', 32); e('cmp', 'x19', 'x24'); e('bc', 'ne', 'zero')
    e('ret')
    pr.label('prefix')
    e('ldr', 'x2', 'x19', 0); e('str', 'x3', 'x19', 0); e('add_r', 'x3', 'x3', 'x2')
    e('add', 'x19', 'x19', 8); e('cmp', 'x19', 'x24'); e('bc', 'ne', 'prefix')
    e('ret')


def gen_merge(pr):
    """Merge over sorted baby (x1) and giant (x2) lists by the 128-bit key, ends x3, x4. Advances the smaller
    head until the keys are equal (returns x9 = 1: the pair routine expands the run, charged per pair) or a
    list ends (returns x9 = 0)."""
    e = pr.emit
    pr.label('merge')
    e('ldp', 'x5', 'x6', 'x1', 0); e('ldp', 'x7', 'x8', 'x2', 0)
    e('subs', 'xzr', 'x5', 'x7'); e('sbcs', 'xzr', 'x6', 'x8'); e('cset', 'x10', 'cc')   # baby < giant
    e('eor', 'x11', 'x5', 'x7'); e('eor', 'x12', 'x6', 'x8'); e('orr', 'x11', 'x11', 'x12')
    e('cbz', 'x11', 'merge_eq')
    e('movi', 'x13', 40); e('cmp', 'x10', 'xzr')
    e('csel', 'x15', 'x13', 'xzr', 'ne'); e('csel', 'x16', 'xzr', 'x13', 'ne')
    e('add_r', 'x1', 'x1', 'x15'); e('add_r', 'x2', 'x2', 'x16')
    e('cmp', 'x1', 'x3'); e('bc', 'cs', 'merge_end'); e('cmp', 'x2', 'x4'); e('bc', 'cc', 'merge')
    pr.label('merge_end')
    e('movi', 'x9', 0); e('ret')
    pr.label('merge_eq')
    e('movi', 'x9', 1); e('ret')


# ---------------------------------------------------------------- table-lookup trace step

KB, DIG = 51, 5                    # digit bits; digits per pair sum (< 2^255)
TB = 1 << 58                       # walker table region: entry (v, t) at TB + 960 v + 32 t, t = 5(i-1) + d
FOLD = 1 << 57                     # fold table: v 2^254 mod p, v < 32


def gen_tabstep(pr, phase):
    """Normal-form traces in a ring at x19; p at x28; x26 = TB, x24 = TB + 480, x23 = 960, x25 = FOLD."""
    e = pr.emit
    sl = lambda k: 32 * ((phase + k) % 12)
    A, S, Tt = ['x0', 'x1', 'x2', 'x3', 'x4'], ['x5', 'x6', 'x7', 'x8'], ['x9', 'x10', 'x11', 'x12']
    e('ldp', 'x13', 'x14', 'x28', 0); e('ldp', 'x15', 'x16', 'x28', 16)
    e('ldp', S[0], S[1], 'x19', sl(0)); e('ldp', S[2], S[3], 'x19', sl(0) + 16)
    e('subs', A[0], 'x13', S[0]); e('sbcs', A[1], 'x14', S[1]); e('sbcs', A[2], 'x15', S[2]); e('sbcs', A[3], 'x16', S[3])
    e('mov', A[4], 'xzr')
    for i in range(1, 7):
        if i < 6:
            e('ldp', S[0], S[1], 'x19', sl(i)); e('ldp', S[2], S[3], 'x19', sl(i) + 16)
            e('ldp', Tt[0], Tt[1], 'x19', sl(12 - i)); e('ldp', Tt[2], Tt[3], 'x19', sl(12 - i) + 16)
            e('adds', S[0], S[0], Tt[0]); e('adcs', S[1], S[1], Tt[1]); e('adcs', S[2], S[2], Tt[2]); e('adc', S[3], S[3], Tt[3])
        else:
            e('ldp', S[0], S[1], 'x19', sl(6)); e('ldp', S[2], S[3], 'x19', sl(6) + 16)
        for d in range(DIG):
            w, b = divmod(KB * d, 64)
            if b + KB <= 64:
                e('lsr', 'x17', S[w], b) if b else e('mov', 'x17', S[w])
            else:
                e('extr', 'x17', S[w + 1], S[w], b)
            e('and_mask', 'x17', 'x17', KB)
            t = (i - 1) * DIG + d
            e('madd', 'x20', 'x17', 'x23', 'x26' if t < 15 else 'x24')
            e('ldp', Tt[0], Tt[1], 'x20', 32 * (t % 15)); e('ldp', Tt[2], Tt[3], 'x20', 32 * (t % 15) + 16)
            e('adds', A[0], A[0], Tt[0]); e('adcs', A[1], A[1], Tt[1]); e('adcs', A[2], A[2], Tt[2])
            e('adcs', A[3], A[3], Tt[3]); e('adc', A[4], A[4], 'xzr')
    e('extr', 'x17', A[4], A[3], 62); e('and_mask', A[3], A[3], 62)
    e('add_lsl', 'x20', 'x25', 'x17', 5)
    e('ldp', Tt[0], Tt[1], 'x20', 0); e('ldp', Tt[2], Tt[3], 'x20', 16)
    e('adds', A[0], A[0], Tt[0]); e('adcs', A[1], A[1], Tt[1]); e('adcs', A[2], A[2], Tt[2]); e('adc', A[3], A[3], Tt[3])
    for _ in range(2):
        e('subs', S[0], A[0], 'x13'); e('sbcs', S[1], A[1], 'x14'); e('sbcs', S[2], A[2], 'x15'); e('sbcs', S[3], A[3], 'x16')
        for k in range(4):
            e('csel', A[k], S[k], A[k], 'cs')
    e('stp', A[0], A[1], 'x19', sl(0)); e('stp', A[2], A[3], 'x19', sl(0) + 16)


def gen_tabbuild(pr):
    """Fills table t from entry (1, t) to the end x21: x20 = address of entry (v, t); x5..x8 = entry (v-1, t);
    x9..x12 = c = (c_i 2^(51 d)) mod p. Each pass writes entry = previous - c mod p and advances x20 by 960."""
    e = pr.emit
    pr.label('tabbuild')
    e('subs', 'x5', 'x5', 'x9'); e('sbcs', 'x6', 'x6', 'x10'); e('sbcs', 'x7', 'x7', 'x11'); e('sbcs', 'x8', 'x8', 'x12')
    e('ldp', 'x13', 'x14', 'x28', 0); e('ldp', 'x15', 'x16', 'x28', 16)
    for r_ in ('x13', 'x14', 'x15', 'x16'):
        e('csel', r_, 'xzr', r_, 'cs')
    e('adds', 'x5', 'x5', 'x13'); e('adcs', 'x6', 'x6', 'x14'); e('adcs', 'x7', 'x7', 'x15'); e('adc', 'x8', 'x8', 'x16')
    e('stp', 'x5', 'x6', 'x20', 0); e('stp', 'x7', 'x8', 'x20', 16)
    e('add', 'x20', 'x20', 960); e('cmp', 'x20', 'x21'); e('bc', 'ne', 'tabbuild')
    e('ret')


# ---------------------------------------------------------------- executed transport

HDR = 64                            # block header words: [0] record count, [8] sequence, [16] END flag
PARW = 52                           # parameter words carried by END messages: z (48 words), delta (4 words)
PARB = 8 * PARW
WIN = 1 << 46                       # receive window of a delivery processor (maps the block being delivered)
assert WIN == (WIN >> 32) << 32 and WIN >> 32 < 1 << 16
PAR = 512                           # parameter copy inside every node's state block
REP = 1024                          # report fields: +0 x (4 words), +32 found, +40 done, +48 expected,
XS, FND, DONE, NEXP, READY, RBLK, MPORT, ENDC, NH, SBLK = 0, 32, 40, 48, 3072, 64, 72, 80, 88, 96   # (+REP)
# READY: the wake word at state + 4096, the only word written in its aligned 2^11-byte block, so a store
# to another word never clears the exclusive monitor of a wait on it (no reservation granule exceeds 2^11 bytes)
DLT = PAR + 48 * 8                  # delta inside the parameter copy


def stride(k):
    return HDR + 40 * k + PARB


def _entry(e):
    """Handler entry: state at x27 (fixed for the processor); the window base WIN is loaded by MOVZ, so nothing
    from the previous message is reused. The acquire load of the sequence word waits for a complete message."""
    for a_, b_ in (('x29', 'x30'), ('x19', 'x20'), ('x21', 'x22'), ('x23', 'x24'), ('x25', 'x26')):
        e('stp_pre', a_, b_, 'sp', -16)
    e('movz_lsl', 'x16', WIN >> 32, 32)                    # the window base is a constant of this processor
    e('add', 'x2', 'x16', 8); e('ldar', 'x2', 'x2')
    e('ldr', 'x3', 'x16', 0)
    e('movi', 'x4', 40); e('madd', 'x24', 'x3', 'x4', 'x16'); e('add', 'x24', 'x24', HDR)
    e('add', 'x1', 'x16', HDR)
    e('ldp', 'x19', 'x20', 'x27', 0); e('ldp', 'x21', 'x22', 'x27', 16)
    e('ldp', 'x23', 'x25', 'x27', 32); e('ldr', 'x26', 'x27', 48)


def _exit(e):
    """Handler exit. The final RET stands for the handler loop's branch back to its entry (one instruction)."""
    e('str', 'x25', 'x27', 40)
    e('movi', 'x11', 1); e('add', 'x12', 'x27', 56); e('stadd', 'x11', 'x12')
    for a_, b_ in (('x25', 'x26'), ('x23', 'x24'), ('x21', 'x22'), ('x19', 'x20'), ('x29', 'x30')):
        e('ldp_post', a_, b_, 'sp', 16)
    e('ret')


def _copy(e, src, dst):
    for q in range(PARW // 2):
        e('ldp', 'x4', 'x5', src, 16 * q); e('stp', 'x4', 'x5', dst, 16 * q)


def gen_transport(pr, K, F=28):
    """Send routine, call sites, the relay and home handlers, flushes, reports and waits. K = records per
    block; F = relay routing shift (production 28: d2 = low 36 bits of record word 2)."""
    e = pr.emit
    paroff = HDR + 40 * K
    assert paroff < 1 << 16
    # SEND: x13 = block, x17 = count, x15 = END flag, x14 = port words (sequence, messages sent)
    pr.label('send')
    e('str', 'x17', 'x13', 0); e('str', 'x15', 'x13', 16)
    e('ldr', 'x10', 'x14', 0); e('add', 'x10', 'x10', 1); e('str', 'x10', 'x14', 0)
    e('add', 'x12', 'x13', 8); e('stlr', 'x10', 'x12')      # completing store: retires after the handler returns
    e('movi', 'x11', 1); e('add', 'x12', 'x14', 8); e('stadd', 'x11', 'x12')
    e('ret')
    # walker call site (tail-called by emit on a full block: x6 = d1, x16 = stride)
    pr.label('callsend_w')
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    e('ldr', 'x13', 'x29', 48); e('madd', 'x13', 'x6', 'x16', 'x13'); e('mov', 'x15', 'xzr'); e('ldr', 'x14', 'x29', 96)
    e('bl', 'send')
    e('ldr', 'x13', 'x29', 64); e('add_lsl', 'x13', 'x13', 'x6', 3); e('str', 'xzr', 'x13', 0)
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')
    # relay call site (x7 = d2, x10 = count, x9 = cursor address)
    pr.label('callsend_r')
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    e('madd', 'x13', 'x7', 'x20', 'x19'); e('mov', 'x17', 'x10'); e('mov', 'x15', 'xzr'); e('ldr', 'x14', 'x27', 64)
    e('bl', 'send')
    e('str', 'xzr', 'x9', 0)
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')
    # RELAY HANDLER (delivery processor 1: links from the two walkers)
    pr.label('h_relay')
    _entry(e)
    e('cmp', 'x1', 'x24'); e('bc', 'eq', 'r_end')
    pr.label('r_loop')
    e('ldp', 'x2', 'x3', 'x1', 0); e('ldp', 'x4', 'x5', 'x1', 16); e('ldr', 'x6', 'x1', 32)
    e('lsl_i', 'x7', 'x4', F); e('lsr', 'x7', 'x7', F)                      # d2
    e('madd', 'x8', 'x7', 'x20', 'x19')
    e('add_lsl', 'x9', 'x21', 'x7', 3); e('ldr', 'x10', 'x9', 0); e('madd', 'x8', 'x10', 'x22', 'x8')
    e('stp', 'x2', 'x3', 'x8', HDR); e('stp', 'x4', 'x5', 'x8', HDR + 16); e('str', 'x6', 'x8', HDR + 32)
    e('add', 'x10', 'x10', 1); e('str', 'x10', 'x9', 0); e('cmp', 'x10', 'x23'); e('bc', 'ne', 'r_next')
    e('bl', 'callsend_r')
    pr.label('r_next')
    e('add', 'x1', 'x1', 40); e('cmp', 'x1', 'x24'); e('bc', 'ne', 'r_loop')
    pr.label('r_end')
    e('ldr', 'x2', 'x16', 16); e('cbz', 'x2', 'r_exit')
    e('movi', 'x2', paroff); e('add_r', 'x2', 'x16', 'x2'); e('add', 'x3', 'x27', PAR)
    _copy(e, 'x2', 'x3')                                                     # parameters from the END message
    e('ldr', 'x4', 'x27', REP + ENDC); e('add', 'x4', 'x4', 1); e('str', 'x4', 'x27', REP + ENDC)
    e('movi', 'x5', 2); e('cmp', 'x4', 'x5'); e('bc', 'ne', 'r_exit')
    e('mov', 'x7', 'xzr'); e('ldr', 'x8', 'x27', REP + NH)
    pr.label('r_flush')                                                      # both walkers ended: flush every home
    e('madd', 'x13', 'x7', 'x20', 'x19'); e('add_lsl', 'x9', 'x21', 'x7', 3); e('ldr', 'x17', 'x9', 0)
    e('movi', 'x2', paroff); e('add_r', 'x2', 'x13', 'x2')
    _copy(e, 'x3', 'x2')
    e('movi', 'x15', 1); e('ldr', 'x14', 'x27', 64); e('bl', 'send'); e('str', 'xzr', 'x9', 0)
    e('add', 'x7', 'x7', 1); e('cmp', 'x7', 'x8'); e('bc', 'ne', 'r_flush')
    pr.label('r_exit')
    _exit(e)
    # HOME HANDLER (delivery processor: the link from its relay)
    pr.label('h_home')
    _entry(e)
    e('cmp', 'x1', 'x24'); e('bc', 'eq', 'h_end')
    pr.label('h_loop')
    e('ldp', 'x2', 'x3', 'x1', 0); e('ldp', 'x4', 'x5', 'x1', 16); e('ldr', 'x6', 'x1', 32)
    e('lsr', 'x7', 'x4', 63); e('madd', 'x8', 'x7', 'x20', 'x19')           # list base
    e('add_lsl', 'x9', 'x21', 'x7', 3); e('ldr', 'x10', 'x9', 0)
    e('cmp', 'x10', 'x23'); e('cset', 'x11', 'cs'); e('orr', 'x25', 'x25', 'x11')   # cap reached -> flag
    e('madd', 'x8', 'x10', 'x22', 'x8'); e('csel', 'x8', 'x26', 'x8', 'cs')         # never write past the cap
    e('stp', 'x2', 'x3', 'x8', 0); e('stp', 'x4', 'x5', 'x8', 16); e('str', 'x6', 'x8', 32)
    e('add', 'x10', 'x10', 1); e('str', 'x10', 'x9', 0)
    e('add', 'x1', 'x1', 40); e('cmp', 'x1', 'x24'); e('bc', 'ne', 'h_loop')
    pr.label('h_end')
    e('ldr', 'x2', 'x16', 16); e('cbz', 'x2', 'h_exit')
    e('movi', 'x2', paroff); e('add_r', 'x2', 'x16', 'x2'); e('add', 'x3', 'x27', PAR)
    _copy(e, 'x2', 'x3')
    e('movi', 'x2', 1); e('str', 'x2', 'x27', REP + READY)                   # wakes the home's main processor
    pr.label('h_exit')
    _exit(e)
    # REPORT HANDLER (relay delivery processor 2: links from its homes; walker 0: links from the relays)
    pr.label('h_rep')
    e('add', 'x2', 'x1', 8); e('ldar', 'x2', 'x2')
    e('ldp', 'x3', 'x4', 'x1', HDR); e('ldp', 'x5', 'x6', 'x1', HDR + 16); e('ldr', 'x7', 'x1', HDR + 32)
    e('add', 'x14', 'x27', REP)
    e('ldp', 'x8', 'x9', 'x14', XS); e('ldp', 'x10', 'x11', 'x14', XS + 16)
    e('cmp', 'x7', 'xzr')
    for a_, b_ in (('x8', 'x3'), ('x9', 'x4'), ('x10', 'x5'), ('x11', 'x6')):
        e('csel', a_, b_, a_, 'ne')                                          # keep a reported x
    e('stp', 'x8', 'x9', 'x14', XS); e('stp', 'x10', 'x11', 'x14', XS + 16)
    e('ldr', 'x12', 'x14', FND); e('orr', 'x12', 'x12', 'x7'); e('str', 'x12', 'x14', FND)
    e('ldr', 'x12', 'x14', DONE); e('add', 'x12', 'x12', 1); e('str', 'x12', 'x14', DONE)
    e('ldr', 'x13', 'x14', NEXP); e('cmp', 'x12', 'x13'); e('bc', 'ne', 'rep_x')
    e('movi', 'x12', 1); e('str', 'x12', 'x14', READY)                        # last report: wake main
    pr.label('rep_x')
    e('ret')
    # REPORT SEND (home or relay main processor): x0..x3 = x, x4 = found
    pr.label('report')
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    e('ldr', 'x13', 'x27', REP + RBLK)
    e('stp', 'x0', 'x1', 'x13', HDR); e('stp', 'x2', 'x3', 'x13', HDR + 16); e('str', 'x4', 'x13', HDR + 32)
    e('movi', 'x17', 1); e('mov', 'x15', 'xzr'); e('ldr', 'x14', 'x27', REP + MPORT); e('bl', 'send')
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')
    pr.label('relay_up')                                                     # relay main after its wait
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    e('add', 'x14', 'x27', REP); e('ldp', 'x0', 'x1', 'x14', XS); e('ldp', 'x2', 'x3', 'x14', XS + 16)
    e('ldr', 'x4', 'x14', FND); e('bl', 'report')
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')
    pr.label('output')                                                       # walker 0 main after its wait
    e('add', 'x14', 'x27', REP); e('ldp', 'x0', 'x1', 'x14', XS); e('ldp', 'x2', 'x3', 'x14', XS + 16)
    e('ldr', 'x4', 'x14', FND)
    e('ret')
    # WAIT on the word at x1 (exclusive load; WFE sleeps until a store to that word or a stale event)
    pr.label('wait')
    e('ldaxr', 'x2', 'x1'); e('cbnz_l', 'x2', 'wait_d'); e('wfe'); e('b', 'wait')
    pr.label('wait_d')
    e('ret')
    # WALKER FLUSH: every relay gets one END message with the remaining records and the parameters
    pr.label('flush_w')
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    e('mov', 'x7', 'xzr'); e('ldr', 'x8', 'x29', 104); e('add', 'x3', 'x29', PAR)
    pr.label('fw_loop')
    e('ldp', 'x13', 'x16', 'x29', 48); e('madd', 'x13', 'x7', 'x16', 'x13')
    e('ldr', 'x9', 'x29', 64); e('add_lsl', 'x9', 'x9', 'x7', 3); e('ldr', 'x17', 'x9', 0)
    e('movi', 'x2', paroff); e('add_r', 'x2', 'x13', 'x2')
    _copy(e, 'x3', 'x2')
    e('movi', 'x15', 1); e('ldr', 'x14', 'x29', 96); e('bl', 'send'); e('str', 'xzr', 'x9', 0)
    e('add', 'x7', 'x7', 1); e('cmp', 'x7', 'x8'); e('bc', 'ne', 'fw_loop')
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')
    # SETUP MESSAGE walker 0 -> walker 1 (delta, s1, s2) and its handler
    pr.label('setup_send')
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    e('ldr', 'x13', 'x29', REP + SBLK); e('add', 'x2', 'x29', DLT)
    e('ldp', 'x3', 'x4', 'x2', 0); e('ldp', 'x5', 'x6', 'x2', 16); e('ldp', 'x7', 'x8', 'x29', 0)
    e('stp', 'x3', 'x4', 'x13', HDR); e('stp', 'x5', 'x6', 'x13', HDR + 16); e('stp', 'x7', 'x8', 'x13', HDR + 32)
    e('movi', 'x17', 1); e('mov', 'x15', 'xzr'); e('ldr', 'x14', 'x29', 96); e('bl', 'send')
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')
    pr.label('h_setup')
    e('add', 'x2', 'x1', 8); e('ldar', 'x2', 'x2')
    e('ldp', 'x3', 'x4', 'x1', HDR); e('ldp', 'x5', 'x6', 'x1', HDR + 16); e('ldp', 'x7', 'x8', 'x1', HDR + 32)
    e('add', 'x9', 'x27', DLT); e('stp', 'x3', 'x4', 'x9', 0); e('stp', 'x5', 'x6', 'x9', 16); e('stp', 'x7', 'x8', 'x27', 0)
    e('movi', 'x2', 1); e('str', 'x2', 'x27', REP + READY)
    e('ret')


def gen_walk(pr):
    """Walker main: walk_init emits the 12 ring traces; walk runs x22 >= 1 trace steps, each followed by emit."""
    e = pr.emit
    pr.label('walk_init')
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    for k in range(12):
        e('ldp', 'x0', 'x1', 'x19', 32 * k); e('ldp', 'x2', 'x3', 'x19', 32 * k + 16); e('bl', 'emit')
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')
    pr.label('walk')
    e('stp_pre', 'x29', 'x30', 'sp', -16)
    pr.label('walk_l')
    for ph in range(12):
        e('bl', 'tstep%d' % ph); e('bl', 'emit'); e('subs_i', 'x22', 'x22', 1); e('bc', 'eq', 'walk_d')
    e('b', 'walk_l')
    pr.label('walk_d')
    e('ldp_post', 'x29', 'x30', 'sp', 16)
    e('ret')


class Proc(Machine):
    """One processor of a node: its own registers, the node's memory. A delivery processor's receive window
    maps, read-only, the block of the message being delivered (win = (sender memory, base, bytes)). With
    coefficients, reads in the walker's table region return the entries the table-build listing writes."""
    def __init__(self, prog, mem=None, coeffs=None, sp0=1 << 40):
        super().__init__(prog)
        self.mem = {} if mem is None else mem
        self.c, self.win = coeffs, None
        self.sp0 = sp0; self.r['sp'] = sp0

    def ld(self, ad):
        assert ad % 8 == 0
        if self.win is not None and WIN <= ad < WIN + (1 << 20):
            smem, base, size = self.win
            assert ad - WIN < size, 'read outside the delivered block'
            return smem.get(base + ad - WIN, 0)
        if self.c is not None:
            if TB <= ad < TB + (1 << KB) * 960:
                v, rem = divmod(ad - TB, 960)
                t, w = divmod(rem, 32)
                i, d = divmod(t, DIG)
                return ((-self.c[i + 1] * v * (1 << (KB * d))) % P >> (64 * (w // 8))) & M64
            if FOLD <= ad < FOLD + 32 * 32:
                v, w = divmod(ad - FOLD, 32)
                return (((v << 254) % P) >> (64 * (w // 8))) & M64
        return self.mem.get(ad, 0)

    def st(self, ad, v):
        assert ad % 8 == 0
        assert not (self.win is not None and WIN <= ad < WIN + (1 << 20)), 'store into a receive window'
        self.mem[ad] = v & M64


def TabMachine(prog, coeffs):
    return Proc(prog, None, coeffs)


# ---------------------------------------------------------------- helpers shared with experiment.py

TARGET = Fraction(3901, 10000)                       # coverage target: 0.39 plus every loss term
DUP_C = 18                                           # expected repeated trace values in all lists <= DUP_C (m+1)^2/r


class WakeMem(dict):
    """A node's memory; records every store into its wake block [lo, lo + 2^11) (set lo before use)."""
    lo = None

    def __setitem__(self, a, v):
        if self.lo <= a < self.lo + 2048:
            self.wake = getattr(self, 'wake', []) + [a]
        dict.__setitem__(self, a, v)


def covers(m, r):
    """Exact rational test mu - mu^2/2 >= 0.3901 for mu = m(12m+1)/r (h = m)."""
    mu = Fraction(m * (12 * m + 1), r)
    return mu - mu * mu / 2 >= TARGET


def stride_params(r):
    """Least m = h with covers(m, r), by exact bisection (the test is monotone while mu < 1)."""
    lo, hi = 1, math.isqrt(r // 12)
    assert not covers(lo, r) and covers(hi, r)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if covers(mid, r) else (mid, hi)
    assert covers(hi, r) and not covers(hi - 1, r)
    return hi, hi


def charpoly(TW, y, q):
    P_ = [None] + [12 * TW.pow(y, k)[0][0][0] % q for k in range(1, 13)]
    e = [1]
    for k in range(1, 13):
        e.append(sum((-1) ** (i - 1) * e[k - i] * P_[i] for i in range(1, k + 1)) * pow(k, -1, q) % q)
    return [((-1) ** (12 - i) * e[12 - i]) % q for i in range(13)]


def mix(x):
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9 & M64
    x = (x ^ (x >> 27)) * 0x94d049bb133111eb & M64
    return x ^ (x >> 31)


def khash4(t, salt):
    h = salt
    for k in range(4):
        h = mix(h ^ ((t >> (64 * k)) & M64))
    return h


def lagrange_min(b1, b2):
    """Lagrange-Gauss reduction of a 2-D integer lattice basis; returns a shortest vector."""
    n = lambda v: v[0] * v[0] + v[1] * v[1]
    if n(b1) > n(b2):
        b1, b2 = b2, b1
    while True:
        mu = (b1[0] * b2[0] + b1[1] * b2[1] + n(b1) // 2) // n(b1)
        b2 = (b2[0] - mu * b1[0], b2[1] - mu * b1[1])
        if n(b2) >= n(b1):
            return b1
        b1, b2 = b2, b1


def walk(c, init, n, q):
    t = list(init[:12])
    while len(t) < n:
        k_ = len(t) - 12
        t.append(-(t[k_] + sum(c[i] * (t[k_ + i] + t[k_ + 12 - i]) for i in range(1, 6)) + c[6] * t[k_ + 6]) % q)
    return t[:n]


def flat(y):
    return [c_ for half in y for f2 in half for c_ in f2]


def unflat(v):
    return tuple(tuple((v[6 * i + 2 * j], v[6 * i + 2 * j + 1]) for j in range(3)) for i in range(2))


def regions_ok(regs):
    regs = sorted(regs)
    return all(a + s <= regs[i + 1][0] for i, (a, s) in enumerate(regs[:-1])) and regs[-1][0] + regs[-1][1] <= 1 << 63
