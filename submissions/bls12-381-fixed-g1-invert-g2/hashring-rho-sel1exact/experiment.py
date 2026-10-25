"""Class-folded distinguished-point rho in GT (BLS12-381, fixed-G1 pairing inversion).

Standard library only; one JSON object in, one out; modes ledger, algebra, toy and sim (`all` runs the four).
Nothing here is a production-scale run; see proof.md section 9 for each mode's scope.
"""
import itertools
import json
from collections import Counter
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

M64 = (1 << 64) - 1
P_FP = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORD = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
MU = Fraction(11, 10)
P_MAIN = Fraction(3903, 10000)

def ln_bounds(p, d, terms=64):
    """Exact rational bounds on ln(d/(d-p)): the partial sum below, the geometric tail above; no floats."""
    x = Fraction(p, d)
    s = Fraction(0)
    xk = Fraction(1)
    for k in range(1, terms + 1):
        xk *= x
        s += xk / k
    xk *= x
    return s, s + xk / ((terms + 1) * (1 - x))

LN_MAIN_LO, LN_MAIN_UP = ln_bounds(P_MAIN.numerator, P_MAIN.denominator)
T_CLAIM = float(MU) * math.log(1 / (1 - float(P_MAIN)))


SLOT_ITEMS = (
    ('selector: 3 blst_fp_mul at 383; kernel price', 3 * 383),
    ('selector: 3 modular subtractions at 22, negation core 13, three 3-LDP loads of p; executed', 3 * 22 + 13 + 9),
    ('selector: extrema 8 + 5*30, eleventh comparison 14; executed', 8 + 5 * 30 + 14),
    ('khash: 2*6*9 + 4 loads; written out', 2 * 6 * 9 + 4),
    ('distinguished-point test, 80 bits; executed', 3),
    ('look-back membership, 4 a miss and 5 a hit; executed', 5),
    ('look-back counted ring write; executed', 10),
    ('bookkeeping: step counter 7, slot update 3 (executed), two branches; written out', 7 + 3 + 2),
    ('table entry address and the dense kernel arguments; written out', 10),
    ('accept tail; written out', 4),
    ('probe ladder, reject path only; written out', 3),
    ('selector operand base, twice; written out', 2),
    ('three blst_fp_mul argument setups; written out', 9),
    ('trail slot counter; written out', 5),
    ('global slot counter; written out', 5),
    ('next-slot branch; written out', 1),
)
C_ITEM = sum(n for _, n in SLOT_ITEMS)
C_MAX_SLOT = C_ITEM - 3
C_SLOT = 1596
assert C_ITEM == 1590, C_ITEM
assert C_MAX_SLOT == 1587, C_MAX_SLOT
SEL_ITEM = sum(n for _, n in SLOT_ITEMS[:3])
assert SEL_ITEM == 1149 + 88 + 172 == 1409, SEL_ITEM
assert C_ITEM - SEL_ITEM == sum(n for _, n in SLOT_ITEMS[3:])
assert C_SLOT - C_MAX_SLOT == 9, "explicit margin (9/1587 = 0.5671%), as proof.md section 5 states"
assert Fraction(C_SLOT - C_MAX_SLOT, C_MAX_SLOT) >= Fraction(5, 1000), "at least 0.5% of it"
W_WINDOW = 24


def mix(x):
    x &= M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)

def khash(kappa, salt, nlimbs):
    """64-bit fixed mixing of the nlimbs 64-bit limbs of a class key."""
    acc = mix(salt)
    for i in range(nlimbs):
        acc = mix(acc ^ ((kappa >> (64 * i)) & M64))
    return acc


class A64:
    """A minimal AArch64 subset machine with architectural flags: a SUBS/SBCS chain's carry is the whole borrow, Z
    the last limb's only, and ADC does not carry."""
    M64 = (1 << 64) - 1

    def __init__(self, sp=0x100000):
        self.x = [0] * 31
        self.sp = sp
        self.mem = {}
        self.N = self.Z = self.C = self.V = 0
        self.taken = []
        self.labels = {}
        self.executed = 0

    def get(self, r):
        if r == 'xzr': return 0
        if r == 'sp': return self.sp
        v = self.x[int(r[1:])]
        return v & 0xFFFFFFFF if r[0] == 'w' else v

    def put(self, r, v):
        if r == 'xzr': return
        if r == 'sp':
            self.sp = v & self.M64
        elif r[0] == 'w':
            self.x[int(r[1:])] = v & 0xFFFFFFFF
        else:
            self.x[int(r[1:])] = v & self.M64

    def word(self, addr): return self.mem.get(addr & self.M64, 0)

    def setword(self, addr, v): self.mem[addr & self.M64] = v & self.M64

    def byte(self, addr):
        """The byte lane of an unaligned address, from the aligned 64-bit word that holds it."""
        addr &= self.M64
        return (self.mem.get(addr & ~7, 0) >> (8 * (addr & 7))) & 0xFF

    def setbyte(self, addr, v):
        addr &= self.M64
        wa = addr & ~7
        sh = 8 * (addr & 7)
        w = self.mem.get(wa, 0)
        self.mem[wa] = (w & ~(0xFF << sh)) | ((v & 0xFF) << sh)

    def lane32(self, addr):
        assert addr % 4 == 0, 'unaligned 32-bit access'
        return (self.word(addr & ~7) >> (8 * (addr & 4))) & 0xFFFFFFFF

    def setlane32(self, addr, v):
        assert addr % 4 == 0, 'unaligned 32-bit access'
        sh = 8 * (addr & 4)
        w = self.word(addr & ~7)
        self.setword(addr & ~7, (w & ~(0xFFFFFFFF << sh)) | ((v & 0xFFFFFFFF) << sh))

    def add(self, a, b, cin=0, flags=False):
        full = a + b + cin
        res = full & self.M64
        if flags:
            self.N = (res >> 63) & 1
            self.Z = 1 if res == 0 else 0
            self.C = 1 if full > self.M64 else 0
        return res

    def sub(self, a, b, cin=1, flags=False):
        full = a - b - (1 - cin)
        res = full & self.M64
        if flags:
            self.N = (res >> 63) & 1
            self.Z = 1 if res == 0 else 0
            self.C = 1 if full >= 0 else 0
        return res

    def cond(self, c):
        return {'EQ': self.Z == 1, 'NE': self.Z == 0, 'CS': self.C == 1, 'HS': self.C == 1,
                'CC': self.C == 0, 'LO': self.C == 0, 'HI': self.C == 1 and self.Z == 0,
                'LS': self.C == 0 or self.Z == 1, 'MI': self.N == 1, 'PL': self.N == 0,
                'GE': self.N == self.V, 'LT': self.N != self.V}[c]

    def val(self, a):
        """An operand group: `#imm`, `xm`, or `xm, LSL #s`."""
        if a[0].startswith('#'): return int(a[0][1:], 0)
        v = self.get(a[0])
        if len(a) > 2 and a[1] == 'LSL': v <<= int(a[2][1:], 0)
        elif len(a) > 2 and a[1] == 'LSR': v >>= int(a[2][1:], 0)
        return v

    def addr(self, a):
        return self.get(a[0]) + self.val(a[1:]) if len(a) > 1 else self.get(a[0])

    def jump(self, label):
        self.taken.append(label)
        return self.labels[label]

    def _exec(self, op, a):
        if op == 'MOVZ': self.put(a[0], int(a[1][1:], 0))
        elif op == 'MOV': self.put(a[0], self.get(a[1]))
        elif op == 'ADD': self.put(a[0], self.get(a[1]) + self.val(a[2:]))
        elif op == 'ADDS': self.put(a[0], self.add(self.get(a[1]), self.val(a[2:]), 0, True))
        elif op == 'ADC': self.put(a[0], self.add(self.get(a[1]), self.val(a[2:]), self.C))
        elif op == 'ADCS': self.put(a[0], self.add(self.get(a[1]), self.val(a[2:]), self.C, True))
        elif op == 'SUB': self.put(a[0], self.get(a[1]) - self.val(a[2:]))
        elif op == 'SUBS': self.put(a[0], self.sub(self.get(a[1]), self.val(a[2:]), 1, True))
        elif op == 'SBC': self.put(a[0], self.sub(self.get(a[1]), self.val(a[2:]), self.C))
        elif op == 'SBCS': self.put(a[0], self.sub(self.get(a[1]), self.val(a[2:]), self.C, True))
        elif op == 'CMP': self.sub(self.get(a[0]), self.val(a[1:]), 1, True)
        elif op == 'AND': self.put(a[0], self.get(a[1]) & self.val(a[2:]))
        elif op == 'ORR': self.put(a[0], self.get(a[1]) | self.val(a[2:]))
        elif op == 'EOR': self.put(a[0], self.get(a[1]) ^ self.val(a[2:]))
        elif op == 'LSL': self.put(a[0], self.get(a[1]) << self.val(a[2:]))
        elif op == 'LSR': self.put(a[0], self.get(a[1]) >> self.val(a[2:]))
        elif op == 'CSEL': self.put(a[0], self.get(a[1]) if self.cond(a[3]) else self.get(a[2]))
        elif op == 'CSETM': self.put(a[0], self.M64 if self.cond(a[1]) else 0)
        elif op == 'LDRB': self.put(a[0], self.byte(self.addr(a[1:])))
        elif op == 'LDR' and a[0][0] == 'w': self.put(a[0], self.lane32(self.addr(a[1:])))
        elif op == 'LDR': self.put(a[0], self.word(self.addr(a[1:])))
        elif op == 'STRB': self.setbyte(self.addr(a[1:]), self.get(a[0]))
        elif op == 'STR' and a[0][0] == 'w': self.setlane32(self.addr(a[1:]), self.get(a[0]))
        elif op == 'STR': self.setword(self.addr(a[1:]), self.get(a[0]))
        elif op == 'LDP':
            ad = self.addr(a[2:])
            self.put(a[0], self.word(ad))
            self.put(a[1], self.word(ad + 8))
        elif op == 'STP':
            ad = self.addr(a[2:])
            self.setword(ad, self.get(a[0]))
            self.setword(ad + 8, self.get(a[1]))
        elif op == 'B': return self.jump(a[0])
        elif op.startswith('B.'): return self.jump(a[0]) if self.cond(op[2:]) else None
        elif op == 'CBZ': return self.jump(a[1]) if self.get(a[0]) == 0 else None
        elif op == 'CBNZ': return self.jump(a[1]) if self.get(a[0]) != 0 else None
        else: raise AssertionError('unsupported instruction %r' % op)
        return None

    def run(self, prog, labels, label=None, limit=100000):
        """Execute from `label` and return the number of instructions executed."""
        self.labels = labels
        self.taken = []
        i = labels[label] if label else 0
        n = 0
        while i < len(prog):
            op, a = prog[i]
            n += 1
            if n > limit: raise AssertionError('the priced list did not terminate')
            tgt = self._exec(op, a)
            i = tgt if tgt is not None else i + 1
        self.executed += n
        return n


def asm(text):
    """Assemble the priced subset: one instruction per line, `label:` alone, `;` starts a comment."""
    prog, labels = [], {}
    for line in text.splitlines():
        line = line.split(';')[0].strip()
        if not line: continue
        if line.endswith(':'):
            labels[line[:-1].strip()] = len(prog)
            continue
        parts = line.replace(',', ' ').replace('[', ' ').replace(']', ' ').split()
        encodable(parts[0], parts[1:])
        prog.append((parts[0], parts[1:]))
    return prog, labels

def bitmask_imm(v, n):
    """An AArch64 logical immediate: a rotated run of ones in a repeated 2..n-bit element, not 0 or all ones."""
    if v <= 0 or v >= (1 << n) - 1: return False
    e = n
    while e > 2 and (v >> (e // 2)) & ((1 << (e // 2)) - 1) == v & ((1 << (e // 2)) - 1): e //= 2
    m = (1 << e) - 1; q = v & m
    rots = [((q >> s) | (q << (e - s))) & m for s in range(e)]
    return any(r & (r + 1) == 0 for r in rots)

def encodable(op, a):
    """Assert an encoding exists: LDP/STP imm7 scaled (X -512..504), LDR/STR imm12 scaled (X 8, W 4, B 1),
    register offsets shifted by 0 or log2(size), ADD/SUB/CMP imm12, MOVZ imm16, bitmask logical immediates."""
    imm = [int(t[1:], 0) for t in a if t.startswith('#')]
    w = 32 if a and a[0][0] == 'w' else 64
    if op in ('LDP', 'STP'):
        sz = w // 8
        assert len(imm) <= 1 and all(-64 * sz <= i <= 63 * sz and i % sz == 0 for i in imm), (op, a)
    elif op in ('LDR', 'STR', 'LDRB', 'STRB'):
        sz = 1 if op.endswith('B') else w // 8
        if imm and 'LSL' not in a:
            assert len(imm) == 1 and 0 <= imm[0] <= 4095 * sz and imm[0] % sz == 0, (op, a)
        elif len(a) > 2:
            assert a[2][0] == 'x' and (len(a) == 3 or (a[3] == 'LSL' and imm[0] in {0, sz.bit_length() - 1})), (op, a)
    elif op in ('ADD', 'SUB', 'ADDS', 'SUBS', 'CMP'):
        assert all(0 <= i <= 4095 for i in imm) and 'LSL' not in a, (op, a)
    elif op in ('AND', 'ORR', 'EOR'):
        assert all(bitmask_imm(i, w) for i in imm), (op, a)
    elif op in ('LSL', 'LSR'):
        assert len(imm) == 1 and 0 <= imm[0] < w, (op, a)
    elif op == 'MOVZ':
        assert len(imm) == 1 and 0 <= imm[0] <= 0xFFFF, (op, a)
    elif op in ('MOV', 'ADC', 'ADCS', 'SBC', 'SBCS', 'CSEL', 'CSETM', 'MUL', 'MADD'):
        assert not imm, (op, a)
    else:
        assert op in ('B', 'CBZ', 'CBNZ') or op.startswith('B.'), ('no encoding rule for', op)


T_OFF, V_OFF, P_OFF = 0, 96, 384
AB_OFF, H_OFF, KS_OFF = 432, 480, 488
Z_OFF, G_OFF = 1680, 1728
S_OFF, KM_OFF, CAP_OFF, CB_OFF = 1808, 1824, 1840, 1944
YA, YB, YX = 2048, 3072, 0x400
C1_A, C1_B = 288, 336

def limbs(x, k=6):
    """The six 64-bit limbs of a non-negative integer (canonical representative, low limb first)."""
    return [(x >> (64 * i)) & M64 for i in range(k)]

def dp_test_prog():
    """The distinguished-point test (3): AND/ORR/CBNZ on 64-bit registers, i.e. 80 bits."""
    return asm("""
    AND x3, x2, #0xFFFF ; the low 16 bits of w_2
    ORR x3, x3, x1 ; the whole of w_1
    CBNZ x3, not_dp
    not_dp:
    """)

def member_prog():
    """The tag-count membership test: 5 executed, 4 on a miss."""
    return asm("""
    member:
    LSR x3, x1, #32
    LDR x4, [sp, #1944]
    LDRB w5, [x4, x3]
    CBZ w5, miss
    B reject
    miss:
    reject:
    """)


def write_prog():
    """The counted ring write (10): old tag load, new tag, base, old count - 1, new count + 1, tag store."""
    return asm("""
    LDR w3, [sp, x27, LSL #2]
    LSR x4, x1, #32
    LDR x5, [sp, #1944]
    LDRB w12, [x5, x3]
    SUB w12, w12, #1
    STRB w12, [x5, x3]
    LDRB w13, [x5, x4]
    ADD w13, w13, #1
    STRB w13, [x5, x4]
    STR w4, [sp, x27, LSL #2]
    """)


def seed_init_prog(w):
    """The first seed (28): the array is zeroed, so set count[tag] = w and T[s] = tag for every slot."""
    lines = ['LSR x4, x1, #32', 'LDR x5, [sp, #1944]', 'MOVZ w6, #%d' % w, 'STRB w6, [x5, x4]']
    lines += ['STR w4, [sp, #%d]' % (4 * s) for s in range(w)]
    return asm('\n'.join(lines))


def reseed_prog(w):
    """The per-trail clear+set (124): decrement the w old tags, then set the new tag's count to w."""
    lines = ['LDR x5, [sp, #1944]']
    for s in range(w):
        lines += ['LDR w3, [sp, #%d]' % (4 * s), 'LDRB w4, [x5, x3]', 'SUB w4, w4, #1', 'STRB w4, [x5, x3]']
    lines += ['LSR x4, x1, #32', 'MOVZ w6, #%d' % w, 'STRB w6, [x5, x4]']
    lines += ['STR w4, [sp, #%d]' % (4 * s) for s in range(w)]
    return asm('\n'.join(lines))

def slot_update_prog(w):
    """The bounded ring-slot update (3): `ADD p, p, #1 ; CMP p, #w ; CSEL p, p, xzr, NE`."""
    return asm("""
    ADD x27, x27, #1
    CMP x27, #%d
    CSEL x27, x27, xzr, NE
    """ % w)

def selector_tail_prog():
    """The selector tail: 188 (182 when M = 0). Min x6-x11 (index x29), max x12-x17 (x30), candidate x0-x5
    (index x18), p in x0-x5 once the candidate is dead; only carries are read, never Z."""
    def ldp(dst, off):
        return ['LDP x%d, x%d, [sp, #%d]' % (dst, dst + 1, off),
                'LDP x%d, x%d, [sp, #%d]' % (dst + 2, dst + 3, off + 16),
                'LDP x%d, x%d, [sp, #%d]' % (dst + 4, dst + 5, off + 32)]
    lines = []
    lines += ldp(6, V_OFF)
    lines += ldp(12, V_OFF)
    lines += ['MOVZ x29, #0', 'MOVZ x30, #0']
    for k in range(1, 6):
        lines.append('MOVZ x18, #%d' % k)
        lines += ldp(0, V_OFF + 48 * k)
        lines.append('SUBS xzr, x0, x6')
        for i in range(1, 6):
            lines.append('SBCS xzr, x%d, x%d' % (i, 6 + i))
        for i in range(6):
            lines.append('CSEL x%d, x%d, x%d, LO' % (6 + i, i, 6 + i))
        lines.append('CSEL x29, x18, x29, LO')
        lines.append('SUBS xzr, x12, x0')
        for i in range(1, 6):
            lines.append('SBCS xzr, x%d, x%d' % (12 + i, i))
        for i in range(6):
            lines.append('CSEL x%d, x%d, x%d, LO' % (12 + i, i, 12 + i))
        lines.append('CSEL x30, x18, x30, LO')
    lines += ldp(0, P_OFF)
    lines.append('SUBS x12, xzr, x12')
    for i in range(1, 6):
        lines.append('SBCS x%d, xzr, x%d' % (12 + i, 12 + i))
    lines.append('B.CS neg_done')
    for i in range(6):
        lines.append('%s x%d, x%d, x%d' % ('ADDS' if i == 0 else 'ADCS', 12 + i, 12 + i, i))
    lines.append('neg_done:')
    lines.append('SUBS xzr, x12, x6')
    for i in range(1, 6):
        lines.append('SBCS xzr, x%d, x%d' % (12 + i, 6 + i))
    lines.append('ADD x30, x30, #6')
    for i in range(6):
        lines.append('CSEL x%d, x%d, x%d, CS' % (6 + i, 6 + i, 12 + i))
    lines.append('CSEL x29, x29, x30, CS')
    return asm('\n'.join(lines))

def p_load_prog():
    """The block's six-limb load of the modulus into x12-x17 by three LDP; the subtraction leaves them intact."""
    return asm('\n'.join(['LDP x12, x13, [sp, #%d]' % P_OFF,
                          'LDP x14, x15, [sp, #%d]' % (P_OFF + 16),
                          'LDP x16, x17, [sp, #%d]' % (P_OFF + 32)]))

X_ROW, Y_ROW = 0, 6
P_LOAD_ROW = 3
MODSUB_ROW = 2 * P_LOAD_ROW + 6 + 1 + 6 + 3
MODSUB_SKIP_ROW = 2 * P_LOAD_ROW + 6 + 1 + 3
NEG_ROW = 6 + 1 + 6
NEG_ROW_WITH_LOAD = NEG_ROW + P_LOAD_ROW
assert (MODSUB_ROW, MODSUB_SKIP_ROW, NEG_ROW, NEG_ROW_WITH_LOAD) == (22, 16, 13, 16)

def modsub_prog(x_off=V_OFF, y_off=V_OFF + 48, out_off=V_OFF + 96):
    """x - y mod p inside x0-x17: 22, or 16 when B.CS skips the correction; p in x12-x17 is never written."""
    lines = []
    for i in range(0, 6, 2):
        lines.append('LDP x%d, x%d, [sp, #%d]' % (X_ROW + i, X_ROW + i + 1, x_off + 8 * i))
        lines.append('LDP x%d, x%d, [sp, #%d]' % (Y_ROW + i, Y_ROW + i + 1, y_off + 8 * i))
    lines.append('SUBS x0, x0, x6')
    for i in range(1, 6):
        lines.append('SBCS x%d, x%d, x%d' % (i, i, Y_ROW + i))
    lines.append('B.CS sub_done')
    for i in range(6):
        lines.append('%s x%d, x%d, x%d' % ('ADDS' if i == 0 else 'ADCS', i, i, 12 + i))
    lines.append('sub_done:')
    for i in range(0, 6, 2):
        lines.append('STP x%d, x%d, [sp, #%d]' % (i, i + 1, out_off + 8 * i))
    return asm('\n'.join(lines))

def legacy_modsub_prog():
    """The rejected shape (p into x20-x25, a mask in x19): the negative control."""
    lines = ['LDR x%d, [sp, #%d]' % (20 + i, P_OFF + 8 * i) for i in range(6)]
    lines.append('CSETM x19, CC')
    return asm('\n'.join(lines))


def plumbing_texts():
    """The rows priced from written-out sequences and not executed (SLOT_ITEMS index: list)."""
    kh = 'LDR x4,[sp,#%d]|LDR x5,[sp,#%d]' % (KM_OFF, KM_OFF + 8)
    for A, o in (('x1', S_OFF), ('x2', S_OFF + 8)):
        kh += '|LDR %s,[sp,#%d]' % (A, o) + ''.join(('|EOR A,A,x%d|LSR x3,A,#30|EOR A,A,x3|MUL A,A,x4|LSR x3,A,#27'
              '|EOR A,A,x3|MUL A,A,x5|LSR x3,A,#31|EOR A,A,x3' % (6 + i)).replace('A', A) for i in range(6))
    return {3: kh,
            7: 'LDR x12,[sp,#%d]|LDR x13,[sp,#%d]|ADDS x25,x25,#1|ADCS x26,x26,xzr|SUBS xzr,x25,x12|SBCS xzr,x26,x13'
               '|B.CS halt|ADD x27,x27,#1|CMP x27,#24|CSEL x27,x27,xzr,NE|B.LO tail|B next' % (CAP_OFF + 32, CAP_OFF + 40),
            8: 'LDR x0,[sp,#%d]|ADD x0,x0,x28|AND x0,x0,#0x3FFF|MOVZ x1,#6912|MADD x2,x0,x1,x19|LDR x3,[sp,#%d]'
               '|MOVZ x4,#576|MADD x2,x3,x4,x2|MOV x1,x20|EOR x0,x20,#0x400' % (H_OFF, KS_OFF),
            9: 'AND x0,x1,#0x3FFF|STP x0,x29,[sp,#%d]|MOVZ x28,#0|EOR x20,x20,#0x400' % H_OFF,
            10: 'ADD x28,x28,#1|CMP x28,#8|B.LT probe', 11: 'EOR x29,x20,#0x400|EOR x29,x20,#0x400',
            12: 'ADD x0,sp,#%d|ADD x1,sp,#%d|ADD x2,sp,#%d|ADD x0,sp,#%d|ADD x1,x29,#%d|ADD x2,sp,#%d|ADD x0,sp,#%d'
                '|ADD x1,sp,#%d|ADD x2,sp,#%d' % (V_OFF + 48, AB_OFF, G_OFF, V_OFF + 96, C1_A, Z_OFF, V_OFF + 144,
                                                  V_OFF + 48, Z_OFF),
            13: 'LDR x12,[sp,#%d]|ADDS x23,x23,#1|ADCS x24,x24,xzr|SUBS xzr,x24,x12|B.CS end' % (CAP_OFF + 8),
            14: 'LDR x12,[sp,#%d]|ADDS x21,x21,#1|ADCS x22,x22,xzr|SUBS xzr,x22,x12|B.CS halt' % (CAP_OFF + 24),
            15: 'B slot'}

def encoding_check():
    """Every priced list through encodable(), with the c1 loads from the y' base, and three refused controls."""
    progs = [dp_test_prog(), member_prog(), write_prog(), seed_init_prog(24), reseed_prog(24), slot_update_prog(24),
             selector_tail_prog(), p_load_prog(), modsub_prog(V_OFF, V_OFF + 48, AB_OFF),
             modsub_prog(V_OFF + 96, V_OFF, V_OFF + 192), modsub_prog(V_OFF + 144, V_OFF + 48, V_OFF + 240),
             asm('\n'.join('LDP x0, x1, [x29, #%d]' % (c + i) for c in (C1_A, C1_B) for i in (0, 16, 32)))]
    for i, t in plumbing_texts().items():
        progs.append(asm(t.replace('|', '\n')))
        assert len(progs[-1][0]) == SLOT_ITEMS[i][1], ('a written-out row is not its priced length', i)
    ins = [x for prog, _ in progs for x in prog]
    top = max(int(a[-1][1:]) for op, a in ins if op in ('LDP', 'STP') and a[-1][0] == '#')
    refused = []
    for t in ('LDP x12 x13 sp #1632', 'STR x4 sp x27 LSL #2', 'STR w4 sp x27 LSL #3'):
        try: encodable(t.split()[0], t.split()[1:])
        except AssertionError: refused.append(t)
    assert len(refused) == 3 and top <= 504, ('the encoding check is vacuous', refused, top)
    return {'instructions_encoding_checked': len(ins), 'largest_ldp_stp_offset': top, 'controls_refused': refused}


RING_SLOTS = W_WINDOW

def ring_slot_update(p, w=W_WINDOW):
    """The reference for the executed slot update: (p+1) mod w."""
    p = p + 1
    if p == w: p = 0
    return p


MEMBER_HIT_ROW = 5
MEMBER_MISS_ROW = 4
WRITE_ROW = 10
SEED_INIT_ROW = 4 + RING_SLOTS
RESEED_ROW = 1 + 4 * RING_SLOTS + 3 + RING_SLOTS
MEMB_ROW = MEMBER_HIT_ROW
assert (MEMBER_HIT_ROW, MEMBER_MISS_ROW, WRITE_ROW, SEED_INIT_ROW, RESEED_ROW, MEMB_ROW) == (5, 4, 10, 28, 124, 5)
assert SLOT_ITEMS[5][1] == MEMB_ROW, "the membership row is the sequence member_prog executes"
assert SLOT_ITEMS[6][1] == WRITE_ROW, "the ring-write row is the sequence write_prog executes"


class HashRing:
    """The production window on toy keys: 24 32-bit tags and a 2^32-byte count array; count > 0 rejects."""
    COUNT_BASE = 0x100000000

    def __init__(self, w, word):
        self.w = w
        self.word = word
        self.m = A64()
        self.m.setword(self.m.sp + CB_OFF, self.COUNT_BASE)
        self.p = 0
        self.seeded = False

    def tag(self, key):
        return (self.word(key) >> 32) & 0xFFFFFFFF

    def _tags(self):
        return [self.m.lane32(self.m.sp + 4 * s) for s in range(self.w)]

    def _count(self, t):
        return self.m.byte(self.COUNT_BASE + t)

    def keys(self):
        """The stored tags in slot order (the multiset the count array mirrors)."""
        return self._tags()

    def seed(self, key):
        """The first seed (28): T all tag(key), count[tag] = w; the array is assumed zeroed."""
        prog, labels = seed_init_prog(self.w)
        self.m.put('x1', self.word(key))
        n = self.m.run(prog, labels)
        assert n == 4 + self.w, n
        self.p = 0
        self.seeded = True
        return n

    def reseed(self, key):
        """The per-trail clear+set (124): decrement the w old tags, then set the new tag's count to w."""
        assert self.seeded, 'a reseed needs a previously seeded ring'
        prog, labels = reseed_prog(self.w)
        self.m.put('x1', self.word(key))
        n = self.m.run(prog, labels)
        assert n == 1 + 4 * self.w + 3 + self.w, n
        self.p = 0
        return n

    def member(self, key):
        """Execute the tag-count test: (rejected, instructions); 5 on a hit, 4 on a miss."""
        prog, labels = member_prog()
        self.m.put('x1', self.word(key))
        n = self.m.run(prog, labels, 'member')
        return ('reject' in self.m.taken), n

    def insert(self, key):
        """Execute the counted ring write (10) into slot p = (p+1) mod w."""
        self.p = ring_slot_update(self.p, self.w)
        prog, labels = write_prog()
        self.m.put('x27', self.p)
        self.m.put('x1', self.word(key))
        n = self.m.run(prog, labels)
        assert n == WRITE_ROW, n
        return n


def dp_test_check(seed=20261025):
    """Execute the 80-bit test, with the w_1 = 2^32, w_2 = 0 witness against the 48-bit form."""
    prog, labels = dp_test_prog()
    rng = random.Random(seed)
    cases = [(0, 0), (0, 0xFFFF), (0, 0x10000), (1, 0), (0, 1), (0xFFFFFFFF, 0), (0, 0xFFFFFFFF),
             (1 << 32, 0), (0, 1 << 32), (1 << 63, 0), ((1 << 64) - 1, 0), (0, (1 << 64) - 1),
             (0x1234567890ABCDEF, 0x1234), (0x1234567890ABCDEF, 0x1235)]
    for _ in range(100):
        cases.append((rng.getrandbits(64), rng.getrandbits(64)))
        cases.append((0, rng.getrandbits(16)))
        cases.append((0, rng.getrandbits(20)))
    distinguished = 0
    for w1, w2 in cases:
        m = A64()
        m.put('x1', w1); m.put('x2', w2)
        n = m.run(prog, labels)
        assert n == 3, n
        got = 'not_dp' not in m.taken
        assert got == (w1 == 0 and (w2 & 0xFFFF) == 0), (w1, w2)
        distinguished += got
    m = A64(); m.put('x1', 1 << 32); m.put('x2', 0)
    m.run(prog, labels)
    witness_64 = 'not_dp' not in m.taken
    witness_48 = ((1 << 32) % (1 << 32) == 0) and (0 & 0xFFFF) == 0
    assert witness_64 is False and witness_48 is True, 'the width witness does not separate the two forms'
    return {'distinguished_point_test_instructions': 3, 'distinguished_point_test_bits': 80,
            'distinguished_point_test_cases_checked': len(cases), 'distinguished_point_test_hits': distinguished,
            'distinguished_point_test_64_bit_form_says_dp_for_w1_2^32_w2_0': witness_64,
            'distinguished_point_test_48_bit_form_says_dp_for_w1_2^32_w2_0': witness_48}

def ring_check(accepts=200, seed=20261024):
    """Execute every window list (slot update, membership, ring write, seed, clear+set) on every fill 0..24,
    with a forced duplicate tag and a forced tag collision."""
    rng = random.Random(seed)
    w = W_WINDOW
    m = A64()
    m.put('x27', 0)
    prog, labels = slot_update_prog(w)
    K = [None] * w
    K[0] = 0
    executed = set()
    for k in range(1, accepts + 1):
        executed.add(m.run(prog, labels))
        p = m.get('x27')
        assert executed == {3}, 'the slot update is not 3 instructions'
        assert p == ring_slot_update((p - 1) % w, w) == k % w, ('ring slot', k, p)
        assert K[p] == (None if k < w else k - w), ('the evicted slot does not hold the oldest key', k, K[p])
        K[p] = k
        live = sorted(x for x in K if x is not None)
        assert live == list(range(k + 1 - min(k + 1, w), k + 1)), ('live set', k, live)
    m.put('x27', 0)
    seq = []
    for _ in range(2 * w):
        m.run(prog, labels)
        seq.append(m.get('x27'))
    assert seq == list(range(1, w)) + [0] + list(range(1, w)) + [0], seq

    hw = lambda kk: khash(kk, 0x5EED0F24, 6)
    tag = lambda kk: hw(kk) >> 32
    scanned = hits = 0
    max_instr = 0
    miss_instr = None
    count_assertions = 0
    collisions = 0
    seed_row = reseed_row = None
    for fill in range(0, w + 1):
        keys = [rng.randrange(1 << 40) for _ in range(fill)]
        ring = HashRing(w, hw)
        seed_row = ring.seed(0xAAAA)
        assert seed_row == SEED_INIT_ROW == 28
        t0 = tag(0xAAAA)
        assert ring._tags() == [t0] * w, ('the first seed does not fill T', fill)
        assert ring._count(t0) == w, ('the first seed does not count w', fill)
        for kk in keys:
            assert ring.insert(kk) == WRITE_ROW == 10, 'the ring write is not 10 instructions'
        touched = {t0} | {tag(kk) for kk in keys}
        counts = Counter(ring._tags())
        for t in touched:
            assert ring._count(t) == counts.get(t, 0), ('the counts do not mirror T', fill, t)
            count_assertions += 1
        window = ([0xAAAA] + keys)[max(0, len(keys) + 1 - w):]
        window_tags = {tag(kk) for kk in window}
        assert set(ring._tags()) == window_tags, ('the tag ring does not hold the window', fill)
        cands = list(window) + [0xBBBB, (1 << 40) + 7, rng.randrange(1 << 40)]
        for cand in cands:
            rejected, instr = ring.member(cand)
            assert rejected == (tag(cand) in window_tags), ('the count ring and the tag set disagree', cand, fill)
            assert instr == (MEMBER_HIT_ROW if rejected else MEMBER_MISS_ROW), ('membership row', instr, rejected)
            if rejected:
                hits += 1
                max_instr = max(max_instr, instr)
            else:
                miss_instr = instr if miss_instr is None else min(miss_instr, instr)
            scanned += 1
    assert max_instr == MEMB_ROW == 5, ('the maximum is a hit at 5', max_instr)
    assert miss_instr == MEMBER_MISS_ROW == 4, ('a miss is 4 instructions', miss_instr)
    assert hits > 0, 'no hit was exercised'

    t0 = tag(0xAAAA)
    def word_syn(kk):
        """A synthetic tag collision: three keys share the seed's tag; all three must reject."""
        return (t0 << 32) if kk in (0xC0FFEE, 0xDEADBEEF, 0xF00D) else hw(kk)
    ring = HashRing(w, word_syn)
    ring.seed(0xAAAA)
    assert ring.insert(0xC0FFEE) == WRITE_ROW
    assert ring.insert(0xDEADBEEF) == WRITE_ROW
    assert ring._count(t0) == w, ('the duplicate tag lost a count', ring._count(t0))
    rejected, instr = ring.member(0xF00D)
    assert (rejected, instr) == (True, MEMBER_HIT_ROW), ('a tag collision must reject conservatively', rejected)
    assert 0xF00D not in (0xAAAA, 0xC0FFEE, 0xDEADBEEF), 'the collision candidate must not be stored'
    collisions += 1
    before = set(ring._tags())
    reseed_row = ring.reseed(0xDDDD)
    assert reseed_row == RESEED_ROW == 124, ('the clear+set is 124', reseed_row)
    tn = tag(0xDDDD)
    assert ring._tags() == [tn] * w, 'the clear+set did not replace every tag'
    assert ring._count(tn) == w, ('the clear+set did not set the new count', ring._count(tn))
    for t_old in before - {tn}:
        assert ring._count(t_old) == 0, ('the clear left an old tag', t_old, ring._count(t_old))

    return {'ring_slot_accepts_checked': accepts, 'ring_slot_assertions': 3 * accepts,
            'ring_slot_update_instructions': 3, 'ring_slot_cycle_length': w,
            'ring_seed_instructions': seed_row, 'ring_reseed_instructions': reseed_row,
            'ring_write_instructions': WRITE_ROW,
            'membership_scan_candidates_checked': scanned,
            'membership_scan_hits_checked': hits,
            'membership_scan_collisions_exercised': collisions,
            'membership_scan_miss_instructions': miss_instr,
            'membership_scan_max_instructions': max_instr,
            'membership_instructions_executed_bound': MEMB_ROW,
            'membership_instructions_charged': SLOT_ITEMS[5][1],
            'tag_count_assertions': count_assertions,
            'distinguished_point_test': dp_test_check(seed)}


def make_field(p):
    """Fp2 = Fp[u]/(u^2+1); Fp12 = Fp2[w]/(w^6 - xi), xi = 1+u; element = 6 Fp2 coefficients of w^i."""
    xi = (1, 1)
    def f2mul(a, b):
        t0 = a[0] * b[0]; t1 = a[1] * b[1]
        return ((t0 - t1) % p, ((a[0] + a[1]) * (b[0] + b[1]) - t0 - t1) % p)
    def f2conj(a): return (a[0], (-a[1]) % p)
    def f2pow(a, e):
        r = (1, 0)
        while e:
            if e & 1: r = f2mul(r, a)
            a = f2mul(a, a); e >>= 1
        return r
    def f12mul(x, y):
        t = [[0, 0] for _ in range(11)]
        for i in range(6):
            xa, xb = x[i]
            for j in range(6):
                ya, yb = y[j]
                t0 = xa * ya; t1 = xb * yb
                c = t[i + j]
                c[0] += t0 - t1
                c[1] += (xa + xb) * (ya + yb) - t0 - t1
        out = []
        for k in range(6):
            a = t[k]
            if k + 6 < 11:
                b = t[k + 6]
                out.append(((a[0] + b[0] - b[1]) % p, (a[1] + b[0] + b[1]) % p))
            else:
                out.append((a[0] % p, a[1] % p))
        return tuple(out)
    one12 = ((1, 0),) + ((0, 0),) * 5
    def f12pow(x, e):
        r = one12
        for bit in bin(e)[2:]:
            r = f12mul(r, r)
            if bit == '1': r = f12mul(r, x)
        return r
    gamma = f2pow(xi, (p - 1) // 6)
    gam = [(1, 0)]
    for _ in range(5): gam.append(f2mul(gam[-1], gamma))
    def frob(x):
        return tuple(f2mul(f2conj(x[i]), gam[i]) for i in range(6))
    zeta = f2mul(gamma, f2conj(gamma))
    assert zeta[1] == 0 and (zeta[0] * zeta[0] - zeta[0] + 1) % p == 0

    assert (gamma[0] + gamma[1]) % p == 0, "gamma is not g0*(1-u) in this tower"
    assert (2 * gamma[0] * gamma[0] - zeta[0]) % p == 0, "zeta = 2*g0^2"
    return dict(p=p, mul=f12mul, pow=f12pow, one=one12, frob=frob, gamma=gamma, zeta=zeta[0],
                f2mul=f2mul, f2pow=f2pow, nlimbs=(p.bit_length() + 63) // 64)

def frobk(F, x, k):
    for _ in range(k % 12): x = F['frob'](x)
    return x

def sel12(F, y):
    """The twelve first coordinates of c1(sigma^k(y)) in 3 Fp multiplications."""
    p = F['p']; zeta = F['zeta']; g0 = F['gamma'][0]
    a, b = y[1]
    t = (a - b) % p
    d = g0 * t % p
    za = zeta * a % p
    zd = zeta * d % p
    h2 = (za - a) % p; hd = (zd - d) % p
    return (a, d, za, zd, h2, hd, (-a) % p, (-d) % p, (-za) % p, (-zd) % p, (-h2) % p, (-hd) % p)

def sel(F, y):
    """Class selector: 3 Fp multiplications; (kappa, kstar) over the twelve first coordinates, a tie to the stored half."""
    p = F['p']; zeta = F['zeta']; g0 = F['gamma'][0]
    a, b = y[1]
    t = (a - b) % p
    d = g0 * t % p
    za = zeta * a % p
    zd = zeta * d % p
    h2 = (za - a) % p; hd = (zd - d) % p
    V = (a, d, za, zd, h2, hd)
    vi = vx = 0
    for k in range(1, 6):
        if V[k] < V[vi]: vi = k
        if V[k] > V[vx]: vx = k
    neg = (p - V[vx]) % p

    if V[vi] <= neg:
        return V[vi], vi
    return neg, vx + 6

def sel_ref(F, y):
    """Direct minimum over the twelve stored representatives: the definition, used to check sel()."""
    V = sel12(F, y)
    best, kb = V[0], 0
    for k in range(1, 12):
        if V[k] < best: best, kb = V[k], k
    return best, kb

def sel12_general(F, y):
    """The same twelve values from the general form d = g0*a + g1*b, as the cross-check; not charged."""
    p = F['p']; zeta = F['zeta']; g0, g1 = F['gamma']
    a, b = y[1]; za = zeta * a % p; d = (g0 * a + g1 * b) % p; zd = zeta * d % p
    h2 = (za - a) % p; hd = (zd - d) % p
    return (a, d, za, zd, h2, hd, (-a) % p, (-d) % p, (-za) % p, (-zd) % p, (-h2) % p, (-hd) % p)


def isqrt_ceil(x):
    """ceil(sqrt(x)) for a non-negative int or Fraction, by exact integer comparisons."""
    if isinstance(x, int):
        r = math.isqrt(x)
        return r if r * r == x else r + 1
    num, den = x.numerator, x.denominator
    r = math.isqrt(num // den)
    while r * r * den < num:
        r += 1
    while r > 0 and (r - 1) * (r - 1) * den >= num:
        r -= 1
    return r

def trap_bound(R, ell_min, ell_max):
    """E(l) = sum_{d<=min(l//2,R)} C(R,d) (d/R)^l 12^-d, evaluated in log2; returns (max, sum) over l."""
    lg = math.lgamma
    def lc(n, k): return (lg(n + 1) - lg(k + 1) - lg(n - k + 1)) / math.log(2)
    mx = -1e9; tot = 0.0
    for ell in range(ell_min, ell_max + 1):
        terms = []
        for d in range(1, min(ell // 2, R) + 1):
            terms.append(lc(R, d) + ell * math.log2(d / R) - d * math.log2(12))
        m = max(terms)
        s = m + math.log2(sum(2.0 ** (t - m) for t in terms))
        mx = max(mx, s); tot += 2.0 ** s
    return mx, tot

def sel_extrema(V):
    """The extrema of the six stored values with the smallest attaining index (strict < and >)."""
    vi = vx = 0
    for k in range(1, 6):
        if V[k] < V[vi]: vi = k
        if V[k] > V[vx]: vx = k
    return vi, vx

def sel_two_half(V, p):
    """The reference for the executed eleventh comparison: min against p - max, a tie to the stored half."""
    vi, vx = sel_extrema(V)
    nb = (p - V[vx]) % p
    return (V[vi], vi) if V[vi] <= nb else (nb, vx + 6)

def sel_direct12(V, p):
    """The definition: the first index of the minimum of the twelve values V_k and (p - V_k) mod p."""
    twelve = tuple(V) + tuple((p - v) % p for v in V)
    best, kb = twelve[0], 0
    for k in range(1, 12):
        if twelve[k] < best: best, kb = twelve[k], k
    return best, kb

def a64_selector(V, p):
    """Execute the priced tail on one six-tuple: (executed, kappa, k*, taken labels)."""
    m = A64()
    for k, v in enumerate(V):
        for i, l in enumerate(limbs(v)):
            m.setword(m.sp + V_OFF + 48 * k + 8 * i, l)
    for i, l in enumerate(limbs(p)):
        m.setword(m.sp + P_OFF + 8 * i, l)
    prog, labels = selector_tail_prog()
    n = m.run(prog, labels)
    kappa = sum(m.get('x%d' % (6 + i)) << (64 * i) for i in range(6))
    return n, kappa, m.get('x29'), m.taken

def live_snapshot(m):
    """x18-x30 as the machine holds them: everything a subtraction block must leave alone."""
    return {r: m.get('x%d' % r) for r in range(18, 31)}

def live_changed(snap, m):
    return [r for r in sorted(snap) if snap[r] != m.get('x%d' % r)]

def modsub_block_check(pairs=40, seed=20261026):
    """Execute blocks (i) and (ii) (two subtractions under ONE load of p) with x19-x29 live and asserted unchanged."""
    F = make_field(P_FP); p0 = F['p']; zeta = F['zeta']; g0 = F['gamma'][0]
    rng = random.Random(seed)
    live = ((19, 0x0000000F00B40000), (20, YA), (21, 0x0123456789ABCDEF), (22, 0xABCDEF01),
            (23, 0xDEADBEEFCAFEBABE), (24, 7), (25, 0x1111222233334444), (26, 0x55555555),
            (27, 17), (28, 5), (29, YA ^ YX))

    def new_machine(pp, V):
        m = A64()
        for k, v in enumerate(V):
            for i, l in enumerate(limbs(v)):
                m.setword(m.sp + V_OFF + 48 * k + 8 * i, l)
        for i, l in enumerate(limbs(pp)):
            m.setword(m.sp + P_OFF + 8 * i, l)
        for r, v in live:
            m.put('x%d' % r, v)
        return m

    def untouched(m, snap, pp):
        bad = live_changed(snap, m)
        assert not bad, ('a subtraction wrote live registers x%s' % bad)
        assert [m.get('x%d' % (12 + i)) for i in range(6)] == limbs(pp), 'the subtraction destroyed p'

    block_i, block_ii = set(), set()
    corrected = skipped = 0
    printed = []
    checked = 0
    load, load_labels = p_load_prog()
    sub1, lab1 = modsub_prog(V_OFF, V_OFF + 48, AB_OFF)
    sub2, lab2 = modsub_prog(V_OFF + 96, V_OFF, V_OFF + 192)
    sub3, lab3 = modsub_prog(V_OFF + 144, V_OFF + 48, V_OFF + 240)
    for pp in (p0, 101):
        cases = [(0, 0), (1, 0), (0, 1), (pp - 1, pp - 1), (pp - 1, 0), (0, pp - 1), (pp - 1, pp - 2),
                 (pp - 2, pp - 1)]
        for _ in range(max(8, pairs // 2)):
            a = rng.randrange(pp); b = rng.randrange(pp)
            cases += [(a, b), (a, (a + 1) % pp), ((a + 1) % pp, a), (a, a)]
        for a, b in cases:
            d = g0 * ((a - b) % pp) % pp; za = zeta * a % pp; zd = zeta * d % pp
            m = new_machine(pp, (a, b, 0, 0, 0, 0)); snap = live_snapshot(m)
            assert m.run(load, load_labels) == P_LOAD_ROW, 'the load of p is three LDP'
            untouched(m, snap, pp)
            mem_before = dict(m.mem)
            n1 = m.run(sub1, lab1)
            block_i.add(n1)
            corrected += n1 == MODSUB_ROW
            skipped += n1 == MODSUB_SKIP_ROW
            assert n1 == (MODSUB_ROW if a < b else MODSUB_SKIP_ROW), ('the branch path', a, b, n1)
            assert sum(m.word(m.sp + AB_OFF + 8 * i) << (64 * i) for i in range(6)) == (a - b) % pp
            untouched(m, snap, pp)
            touched = sorted(x for x in set(mem_before) | set(m.mem) if mem_before.get(x) != m.mem.get(x))
            assert touched == [m.sp + AB_OFF + 8 * i for i in range(6)], ('the frame writes', touched)
            m = new_machine(pp, (a, d, za, zd, 0, 0)); snap = live_snapshot(m)
            assert m.run(load, load_labels) == P_LOAD_ROW, 'the load of p is three LDP'
            na = m.run(sub2, lab2)
            untouched(m, snap, pp)
            nb = m.run(sub3, lab3)
            untouched(m, snap, pp)
            block_ii.add(P_LOAD_ROW + na + nb)
            r1 = sum(m.word(m.sp + V_OFF + 192 + 8 * i) << (64 * i) for i in range(6))
            r2 = sum(m.word(m.sp + V_OFF + 240 + 8 * i) << (64 * i) for i in range(6))
            assert (r1, r2) == ((za - a) % pp, (zd - d) % pp), ('block (ii)', a, b, pp)
            checked += 1
            if pp == p0 and a != b:
                printed = (printed + [dict(x1=hex(za), y1=hex(a), r1=hex(r1), x2=hex(zd), y2=hex(d), r2=hex(r2),
                                          executed=6 + na + nb)])[-2:]
    assert block_i == {MODSUB_SKIP_ROW, MODSUB_ROW} and corrected > 0 and skipped > 0, \
        ('both branch paths must be exercised', sorted(block_i), corrected, skipped)
    assert block_ii == {P_LOAD_ROW + 2 * MODSUB_SKIP_ROW, P_LOAD_ROW + MODSUB_SKIP_ROW + MODSUB_ROW, P_LOAD_ROW + 2 * MODSUB_ROW}, sorted(block_ii)
    m = new_machine(p0, (1, 1, 0, 0, 0, 0))
    snap = live_snapshot(m)
    control_prog, control_labels = legacy_modsub_prog()
    m.run(control_prog, control_labels)
    control = live_changed(snap, m)
    assert set(range(19, 26)) <= set(control), ('the live-register check is vacuous', control)
    return {'block_i_counts': sorted(block_i), 'block_ii_counts': sorted(block_ii),
            'block_i_row_with_load': max(block_i) + P_LOAD_ROW, 'block_ii_row_with_load': max(block_ii),
            'subtraction_row': MODSUB_ROW, 'subtraction_skip_row': MODSUB_SKIP_ROW,
            'negation_row': NEG_ROW_WITH_LOAD, 'modulus_load_instructions': P_LOAD_ROW,
            'subtractions_correcting': corrected, 'subtractions_skipping': skipped,
            'live_registers_preserved_after_each_subtraction': True, 'modulus_registers_unchanged_after_each': True,
            'frame_words_written_per_subtraction': 6, 'destructive_control_detected': control,
            'tuples_executed': checked, 'tuples_printed': printed}

def selector_min_check(tuples_=240, pairs=40, seed=20261025):
    """Check the two-half minimum against the direct scan and execute the tail and the subtraction blocks."""
    ps = 101
    vs = (0, 1, 2, 3, 50, 99, 100)
    n = 0
    for base in itertools.product(vs, repeat=6):

        assert sel_extrema(base) == (min(range(6), key=lambda k: base[k]),
                                     max(range(6), key=lambda k: base[k])), base
        assert sel_two_half(base, ps) == sel_direct12(base, ps), base
        n += 1
    rng = random.Random(seed)

    sample = [(0, 0, 0, 0, 0, 0), (ps - 1, 1, 1, 1, 1, 1), (ps - 1, 0, 0, 0, 0, 0),
              (1, ps - 1, 2, 3, 4, 5), (50, 51, 50, 51, 50, 51)]
    for _ in range(tuples_):
        sample.append(tuple(rng.choice(vs) for _ in range(6)))
    tail_row = 8 + 5 * 30 + 14
    ties = zeros = 0
    printed = []
    tail_rows = set()
    for V in sample:
        n_ex, kappa, kstar, taken = a64_selector(V, ps)
        zero_max = 'neg_done' in taken
        assert zero_max == all(v == 0 for v in V), ('the negation path', V, taken)
        assert n_ex == 188 - (6 if zero_max else 0), ('the tail is 188 instructions, 182 when M = 0', n_ex)
        tail_rows.add(n_ex)
        assert (kappa, kstar) == sel_two_half(V, ps) == sel_direct12(V, ps), ('the executed tail', V)
        vi, vx = sel_extrema(V)
        if V[vi] == (ps - V[vx]) % ps:
            ties += 1
            if len(printed) < 2: printed.append(('tie', V, kappa, kstar))
        if all(v == 0 for v in V):
            zeros += 1
            printed.append(('all-zero', V, kappa, kstar))
    assert ties > 0, 'the executed eleventh comparison never met a tie between the halves'
    assert zeros > 0, "the executed negation never met M = 0"
    assert tail_rows == {182, 188}, ('both negation paths must be executed', sorted(tail_rows))

    F = make_field(P_FP)
    p = P_FP
    zeta = F['zeta']; g0 = F['gamma'][0]
    for _ in range(pairs):
        a = rng.randrange(p); b = rng.randrange(p)
        d = g0 * ((a - b) % p) % p
        za = zeta * a % p; zd = zeta * d % p
        V = (a, d, za, zd, (za - a) % p, (zd - d) % p)
        n_ex, kappa, kstar, taken = a64_selector(V, p)
        assert n_ex == 188, n_ex
        assert 'neg_done' not in taken, 'a nonzero maximum takes the correction'
        assert (kappa, kstar) == sel_two_half(V, p) == sel_direct12(V, p), (a, b)
        if len(printed) < 4: printed.append(('production', V, kappa, kstar))


    sub_res = modsub_block_check(pairs=pairs, seed=seed + 1)
    assert 3 * MODSUB_ROW + NEG_ROW + 3 * P_LOAD_ROW == SLOT_ITEMS[1][1], 'the alloc row is the executed sequences'
    assert sub_res['block_i_row_with_load'] == 25 and sub_res['block_ii_row_with_load'] == 47
    return {'exhaustive_tuples': n,
            'small_field_tuples_executed': len(sample),
            'production_tuples_executed': pairs,
            'tied_halves_exercised': ties,
            'all_zero_tuples_exercised': zeros,
            'tail_instructions_executed': 188,
            'tail_counts_executed': sorted(tail_rows),
            'extrema_row_instructions': 8 + 5 * 30,
            'negation_instructions': NEG_ROW,
            'eleventh_comparison_instructions': 14,
            'tail_instructions_priced': tail_row,
            'modular_subtraction_blocks': sub_res,
            'modular_subtraction_row': MODSUB_ROW,
            'modulus_load_instructions': P_LOAD_ROW,
            'tuples_printed': [{'kind': k, 'V': [hex(v) for v in V], 'kappa': hex(kappa), 'kstar': ks}
                               for k, V, kappa, ks in printed]}

def ledger():
    getcontext().prec = 80
    N = (R_ORD - 1) // 12
    mu = MU


    s2 = 2 * mu * N * LN_MAIN_UP
    S_max = isqrt_ceil(s2)
    assert Fraction(S_max) ** 2 >= s2 > Fraction(S_max - 1) ** 2

    xt = Fraction(P_MAIN.numerator, P_MAIN.denominator)
    part, xtk = Fraction(0), Fraction(1)
    for k in range(1, 71):
        xtk *= xt
        part += xtk / k
    assert 0 < LN_MAIN_UP - part < Fraction(1, 1 << 99)
    up70 = part + xtk * xt / (71 * (1 - xt))
    assert LN_MAIN_UP > up70 > part, "U is the wider upper bound"
    assert S_max - isqrt_ceil(2 * mu * N * up70) > (1 << 25)

    assert S_max - isqrt_ceil(2 * mu * N * LN_MAIN_LO) <= (S_max >> 90)
    M_max = S_max + (S_max >> 13) + 1
    c_over = C_SLOT
    unit = Fraction(20000 + c_over, 20000)
    assert unit == Fraction(5399, 5000)
    sel_res = selector_min_check()
    enc_res = encoding_check()
    ring_res = ring_check()
    R = 1 << 14
    theta_bits = 80
    L_cap = 1 << 87
    D_max = 1 << 47
    ST_max = 1 << 47
    REC_max = 16
    add = {}
    add['pairings (g and one spare), 2*2^30 instr'] = Fraction(2 * (1 << 30), 20000)
    add['table: R*(1100 + 12) units'] = R * 1112
    add['starts: ST_max*1100 units'] = ST_max * 1100
    add['reports: D_max*2^16 instr'] = Fraction(D_max * (1 << 16), 20000)
    add['replays: REC_max*2*L_cap*(5/4)'] = Fraction(REC_max * 2 * L_cap * 5, 4)
    add['finish last trail: L_cap*unit'] = L_cap * unit
    add['G2 scalar mult + GT check (<= 2^12 units)'] = 1 << 12
    add['tag count array: 2^32 byte stores'] = Fraction(1 << 32, 20000)
    A = sum(add.values())
    W = unit * M_max + A

    S_lb = isqrt_ceil(2 * mu * N * part)
    assert S_lb <= S_max and W / (unit * (S_lb + (S_lb >> 13) + 1) + A) - 1 < Fraction(1, 1 << 99)

    log2W = (Decimal(W.numerator) / Decimal(W.denominator)).ln() / Decimal(2).ln()

    lo, hi = 0, 100 * W.numerator.bit_length()
    while lo < hi:
        mid = (lo + hi) // 2
        if W.numerator ** 100 <= (1 << mid) * W.denominator ** 100:
            hi = mid
        else:
            lo = mid + 1
    c100 = lo
    exact_le = W.numerator ** 100 <= (1 << c100) * W.denominator ** 100
    exact_gt_prev = W.numerator ** 100 > (1 << (c100 - 1)) * W.denominator ** 100
    claim = Fraction(c100, 100)
    pre = (add['table: R*(1100 + 12) units'] + add['pairings (g and one spare), 2*2^30 instr']
           + add['tag count array: 2^32 byte stores'])
    mem = {
        'reports_bytes_max': D_max * 96,
        'twist_table_bytes': 12 * R * 576,
        'uv_scalars_bytes': 1 << 20,
        'window_and_state_bytes': 1 << 20,
        'tag_count_array_bytes': 1 << 32,
    }
    peak = sum(mem.values())
    trap_max, trap_sum = trap_bound(R, 25, 512)
    expected_dp = Decimal(S_max) * Decimal(2) ** (-theta_bits)

    log2_S = Decimal(S_max).ln() / Decimal(2).ln()
    assert S_max < (1 << 127) and S_max > (1 << 95), 'the two report rates either side of the 2^47 cap'
    expected_reports_48 = log2_S - 48

    def hundredth(c):
        """The least integer m with W(c)^100 <= 2^m, by exact integer binary search, at charge c."""
        Wc = Fraction(20000 + c, 20000) * M_max + A
        lo2, hi2 = 0, 100 * Wc.numerator.bit_length()
        while lo2 < hi2:
            mid2 = (lo2 + hi2) // 2
            if Wc.numerator ** 100 <= (1 << mid2) * Wc.denominator ** 100: hi2 = mid2
            else: lo2 = mid2 + 1
        return lo2
    def lg2(v):
        v = Fraction(v)
        return float((Decimal(v.numerator) / Decimal(v.denominator)).ln() / Decimal(2).ln())
    ceiling = c_over
    while hundredth(ceiling + 1) == hundredth(c_over): ceiling += 1
    return {
        'log2_r': float(Decimal(R_ORD).ln() / Decimal(2).ln()), 'log2_N_classes': float(Decimal(N).ln() / Decimal(2).ln()),
        'mu': str(mu), 'S_max': S_max, 'log2_S_max': float(Decimal(S_max).ln() / Decimal(2).ln()),
        'M_max': M_max, 'unit_per_mult': '%d/%d (= %s)' % (20000 + c_over, 20000, unit),
        'additive_items_log2': {k: lg2(v) for k, v in add.items()},
        'log2_additive_total': lg2(A),
        'log2_W_max': str(log2W)[:12], 'claim_time_log2': '%.2f' % float(claim),
        'claim_from_exact_integer_search_hundredths': c100,
        'claim_check_W_le_2^claim': bool(Decimal(c100) >= log2W * 100),
        'exact_integer_check_W^100_le_2^%d' % c100: exact_le,
        'exact_integer_check_W^100_gt_2^%d' % (c100 - 1): exact_gt_prev,
        'slot_items_instructions': {k: n for k, n in SLOT_ITEMS},
        'selector_instructions_total': SEL_ITEM,
        'slot_instructions_itemised': C_ITEM, 'slot_instructions_max_slot': C_MAX_SLOT,
        'slot_instructions_charged': c_over, 'slot_margin_instructions': c_over - C_ITEM,
        'slot_margin_instructions_on_max_slot': c_over - C_MAX_SLOT,

        'sensitivity_bits_per_instr_charged': float(Fraction(M_max, 20000) / W) / math.log(2),
        'sensitivity_bits_per_instr_at_zero_charge': 1.0 / (20000 * math.log(2)),
        'selector_two_half_minimum_tuples_checked': sel_res['exhaustive_tuples'],
        'selector_sequences_executed': sel_res,
        'ring_and_equality_sequences': ring_res,
        'encodability': enc_res,
        'log2_preprocessing': lg2(pre),
        'expected_reports_log2': float(expected_dp.ln() / Decimal(2).ln()), 'report_cap_log2': 47,
        'expected_reports_log2_if_48_bit_test': float(expected_reports_48),
        'charge_ceiling_leaving_claim_unchanged': ceiling, 'charge_room_above_charged': ceiling - c_over,
        'peak_memory_log2_bytes': math.log2(peak),
        'trap_lemma_log2_max_E_l_for_l_25_to_512': trap_max, 'trap_lemma_log2_sum': math.log2(trap_sum),
        'trap_cost_exponent_log2(L_cap*sum)': 87 + math.log2(trap_sum),
        'success_main_event': '3903/10000 under mu', 'success_claim': '390029/1000000 (>= 0.39)',
    }


def algebra(elements, seed):
    F = make_field(P_FP); p = P_FP
    rng = random.Random(seed)
    lam = p % R_ORD
    out = {'lambda_order_12': pow(lam, 12, R_ORD) == 1 and pow(lam, 6, R_ORD) == R_ORD - 1 and pow(lam, 4, R_ORD) != 1,
           'zeta_primitive_6th_root_in_Fp': (F['zeta'] ** 6) % p == 1 and (F['zeta'] ** 3) % p == p - 1,
           'gamma_hex': [format(F['gamma'][0], 'x'), format(F['gamma'][1], 'x')], 'zeta_hex': format(F['zeta'], 'x'),
           'gamma_equals_g0_times_1_minus_u(g0_plus_g1_is_zero)': (F['gamma'][0] + F['gamma'][1]) % p == 0,
           'zeta_equals_2_g0_squared': (2 * F['gamma'][0] ** 2 - F['zeta']) % p == 0,
           'g0_is_a_6th_root_of_minus_one_eighth': (8 * pow(F['gamma'][0], 6, p) + 1) % p == 0,
           'selector_fp_multiplications_charged': 3,
           'selector_fast_form_equals_general_form_on_GT': True,
           'checked_elements': 0}
    e = (p ** 12 - 1) // R_ORD
    first = True
    for _ in range(elements):
        x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
        y = F['pow'](x, e)
        assert F['pow'](y, R_ORD) == F['one'] and y != F['one']
        if first:
            assert F['frob'](y) == F['pow'](y, p) == F['pow'](y, lam)
            first = False
        kap0, k0 = sel(F, y)
        assert sel(F, y) == sel_ref(F, y)
        Vfast = sel12(F, y)
        out['selector_fast_form_equals_general_form_on_GT'] &= (Vfast == sel12_general(F, y))
        assert Vfast == sel12_general(F, y)
        for k in range(12):
            assert frobk(F, y, k)[1][0] == Vfast[k]
        for m in range(12):
            ym = frobk(F, y, m)
            kap, km = sel(F, ym)
            assert (kap, km) == sel_ref(F, ym)
            assert kap == kap0 and km == (k0 - m) % 12
        t = F['pow'](y, rng.randrange(2, R_ORD))
        for m in range(12):
            ym = frobk(F, y, m); _, km = sel(F, ym)
            y1 = F['mul'](ym, frobk(F, t, (-km) % 12))
            y0 = F['mul'](y, frobk(F, t, (-k0) % 12))
            assert any(frobk(F, y0, j) == y1 for j in range(12))
        out['checked_elements'] += 1
    out['all_assertions_passed'] = True
    return out


def is_prime(n):
    if n < 2: return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0: return n == q
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True

def find_order(bits, rng):
    while True:
        c = rng.randrange(1 << (bits - 1), 1 << bits) // 12 * 12 + 1
        if is_prime(c): return c

def prim12(r, rng):
    while True:
        lam = pow(rng.randrange(2, r - 1), (r - 1) // 12, r)
        if pow(lam, 6, r) == r - 1 and pow(lam, 4, r) != 1: return lam

def toy_field(rp, lam, qbits=60):
    """prime q = lam (mod rp), q = 7 (mod 12), xi = 1+u neither a square nor a cube in Fq2."""
    s0 = ((7 - lam) * pow(rp, -1, 12)) % 12
    j = (1 << qbits) // (12 * rp)
    while True:
        q = lam + rp * (s0 + 12 * j)
        j += 1
        if q % 12 == 7 and is_prime(q) and pow(2, (q - 1) // 2, q) == q - 1:
            F = make_field_checked(q)
            if F: return F

def make_field_checked(q):
    try:
        F = make_field(q)
    except AssertionError:
        return None
    one2 = (1, 0)
    if F['f2pow']((1, 1), (q * q - 1) // 3) == one2: return None
    return F

class Walk:
    """The emulated walk's index and distinguished-point bits, from independent salts."""
    def __init__(self, nlimbs, R, w, db, A=8, seed=0):
        rng = random.Random(seed)
        self.R, self.w, self.A, self.db, self.nl = R, w, A, db, nlimbs
        self.s_h = rng.getrandbits(64); self.s_dp = [rng.getrandbits(64) for _ in range((db + 63) // 64)]
    def word1(self, kappa): return khash(kappa, self.s_h, self.nl)
    def hidx(self, kappa): return self.word1(kappa) % self.R
    def is_dp(self, kappa):
        left = self.db
        for s in self.s_dp:
            take = min(64, left); left -= take
            if khash(kappa, s, self.nl) & ((1 << take) - 1): return False
        return True


def toy(bits, trials, seed):
    rng = random.Random(seed)
    res = []
    for t in range(trials):
        rp = find_order(bits, rng); lam = prim12(rp, rng)
        F = toy_field(rp, lam)
        q = F['p']
        assert (q % rp) == lam and pow(lam, 4, rp) != 1
        e = (q ** 12 - 1) // rp
        while True:
            x = tuple((rng.randrange(q), rng.randrange(q)) for _ in range(6))
            g = F['pow'](x, e)
            if g != F['one']: break
        a = rng.randrange(1, rp)
        z = F['pow'](g, a)
        assert F['pow'](z, rp) == F['one']
        N = (rp - 1) // 12
        res.append(toy_instance(F, g, z, a, rp, lam, N, rng.getrandbits(32), rng))
    ok = sum(r['recovered'] for r in res)
    return {'trials': trials, 'bits': bits, 'recovered_correct': ok,
            'selector_form': 'g0*(a-b): 3 Fp multiplications, exercised by every toy instance (its tower forces gamma = g0*(1-u))',
            'wrong_outputs': sum(r['wrong'] for r in res),
            'within_budget': sum(1 for r in res if r['recovered'] and r['S'] <= r['S_max']),
            'mean_steps_over_sqrt_N': round(sum(r['S'] for r in res) / trials / math.sqrt(sum(r['N'] for r in res) / trials), 3),
            'replay_matched_original_path': all(r['replay_ok'] for r in res),
            'probes_total': sum(r['probes'] for r in res), 'escapes_total': sum(r['escapes'] for r in res),
            'dead_trails': sum(r['dead'] for r in res), 'budget_S_max_mean': round(sum(r['S_max'] for r in res) / trials, 1),
            'scope': 'toy tower of the production shape (q ~ 2^60, r ~ 2^%d); no production cost or success claim' % bits}

def toy_instance(F, g, z, a, rp, lam, N, wseed, rng, R=64, w=12, db=3):
    """One toy instance: the complete algorithm with the production window."""
    mul = F['mul']; nl = F['nlimbs']
    W = Walk(nl, R, w, db, seed=wseed)
    lamp = [pow(lam, k, rp) for k in range(12)]
    lami = [pow(lam, (-k) % 12, rp) for k in range(12)]
    U = [rng.randrange(rp) for _ in range(R)]; V = [rng.randrange(rp) for _ in range(R)]
    TW = []
    for h in range(R):
        Th = mul(F['pow'](g, U[h]), F['pow'](z, V[h]))
        fr = [Th]
        for _ in range(11): fr.append(F['frob'](fr[-1]))
        TW.append([fr[(12 - k) % 12] for k in range(12)])
    S_max = isqrt_ceil(2 * MU * N * LN_MAIN_UP)
    stats = dict(probes=0, escapes=0, dead=0)

    def trail(tid, track):
        r2 = random.Random(wseed * 7919 + tid)
        al = r2.randrange(rp); be = r2.randrange(rp)
        y = mul(F['pow'](g, al), F['pow'](z, be))
        kap, k = sel(F, y)
        if tid == 1:

            Vf = sel12(F, y)
            assert all(frobk(F, y, j)[1][0] == Vf[j] for j in range(12))
            assert Vf == sel12_general(F, y)
        win = [kap]
        ring = HashRing(W.w, W.word1)
        tag = lambda kk: W.word1(kk) >> 32
        ring.seed(kap)
        n = 0; steps = 0; slots = 0
        while True:
            steps += 1
            if W.is_dp(kap): return ('dp', kap, y, al, be, steps)
            if slots >= 64 * (1 << db): return ('cap', kap, y, al, be, steps)
            h0 = W.hidx(kap); ok = False
            for ai in range(W.A):
                hh = (h0 + ai) % R
                slots += 1
                stats['probes'] += 1 if track is False else 0
                y2 = mul(y, TW[hh][k])
                kap2, k2 = sel(F, y2)
                hit, _ = ring.member(kap2)
                assert hit == (tag(kap2) in {tag(k) for k in win}), 'the tag ring and the tag list disagree'
                if hit:
                    stats['escapes'] += 1 if track is False else 0
                    continue
                ok = True; break
            if not ok:
                stats['dead'] += 1 if track is False else 0
                return ('dead', kap, y, al, be, steps)
            if track:
                al = (al + lami[k] * U[hh]) % rp; be = (be + lami[k] * V[hh]) % rp
            y, kap, k = y2, kap2, k2
            win.append(kap)
            if len(win) > W.w: win.pop(0)
            ring.insert(kap)
            assert sorted(set(ring.keys())) == sorted({tag(k) for k in win}), 'the tag ring does not hold the window tags'
            n += 1

    store = {}; S = 0; tid = 0; recovered = False; wrong = 0; replay_ok = True
    while True:
        tid += 1
        status, kap, y, al, be, steps = trail(tid, False)
        S += steps
        if status == 'dp':
            if kap in store:
                t1 = store[kap]
                A_ = trail_track(trail, t1); B_ = trail_track(trail, tid)

                replay_ok = replay_ok and A_[1] == kap and B_[1] == kap
                yA, alA, beA = A_[2], A_[3], A_[4]; yB, alB, beB = B_[2], B_[3], B_[4]
                for m in range(12):
                    if frobk(F, yA, m) == yB:
                        lm = lamp[m]; den = (lm * beA - beB) % rp
                        if den:
                            cand = ((alB - lm * alA) * pow(den, -1, rp)) % rp
                            if F['pow'](g, cand) == z:
                                recovered = True
                                if cand != a: wrong += 1
                        break
                if recovered: break
            else:
                store[kap] = tid
        if S >= 6 * S_max: break
    return dict(recovered=recovered, wrong=wrong, S=S, S_max=S_max, N=N, replay_ok=replay_ok, **stats)

def trail_track(trail, tid): return trail(tid, True)


def sim(bits, R, trials, seed, w=24, db=5):
    rng = random.Random(seed)
    tot = 0.0; n_first = 0; esc = 0; att = 0; dead = 0; merged = 0; coalesced = 0; within = 0
    for _ in range(trials):
        rp = find_order(bits, rng); lam = prim12(rp, rng)
        N = (rp - 1) // 12
        W = Walk(1, R, w, db, seed=rng.getrandbits(32))
        salt = rng.getrandbits(64)
        ip = [pow(lam, (-k) % 12, rp) for k in range(12)]
        a = rng.randrange(1, rp)
        U = [rng.randrange(rp) for _ in range(R)]; V = [rng.randrange(rp) for _ in range(R)]
        Mh = [(U[h] + a * V[h]) % rp for h in range(R)]
        def ssel(x):
            best = None; kb = 0; xk = x
            for k in range(12):
                hv = mix(xk ^ salt)
                if best is None or hv < best: best = hv; kb = k
                xk = xk * lam % rp
            return best, kb
        visited = {}; end_key = {}; steps = 0; first = None; tid = 0; limit = int(12 * math.sqrt(N))
        follow = None
        while (first is None or follow is not None) and steps < limit + 4 * (1 << db) * 16:
            tid += 1
            x = rng.randrange(1, rp)
            kap, k = ssel(x); win = [kap]; sl = 0
            while True:
                if first is None:
                    steps += 1
                    pv = visited.get(kap)
                    if pv is not None and pv != tid:
                        first = steps; follow = pv
                    else:
                        visited.setdefault(kap, tid)
                if W.is_dp(kap) or sl >= 16 * (1 << db) or (first is None and steps >= limit): break
                h0 = W.hidx(kap); good = False
                for ai in range(W.A):
                    hh = (h0 + ai) % R
                    att += 1; sl += 1
                    x2 = (x + ip[k] * Mh[hh]) % rp
                    kap2, k2 = ssel(x2)
                    if kap2 in win: esc += 1; continue
                    good = True; break
                if not good: dead += 1; break
                x, kap, k = x2, kap2, k2
                win.append(kap)
                if len(win) > w: win.pop(0)
            end_key[tid] = kap if W.is_dp(kap) else None
            if follow is not None:
                if tid != follow:
                    merged += 1
                    if end_key[tid] is not None and end_key[tid] == end_key.get(follow): coalesced += 1
                    follow = None
        if first is not None:
            n_first += 1; tot += first * first / (2 * N)
            if first * first / (2 * N) <= T_CLAIM: within += 1

    ctl = 0.0; cn = min(trials, 300)
    for _ in range(cn):
        rp = find_order(bits, rng); N = (rp - 1) // 12; seen = set(); s = 0
        while True:
            v = rng.randrange(N); s += 1
            if v in seen: break
            seen.add(v)
        ctl += s * s / (2 * N)
    return {'bits': bits, 'R': R, 'w': w, 'dp_bits': db, 'trials': trials, 'repeats_found': n_first,
            'mu_hat_first_inter_trail_repeat': round(tot / max(1, n_first), 3),
            'std_error_approx': round(1 / math.sqrt(max(1, n_first)), 3),
            'groups_first_repeat_within_S_max_budget': within,
            'fraction_first_repeat_within_S_max_budget': round(within / trials, 4),
            'budget_threshold_S^2/(2N)': round(T_CLAIM, 5),
            'ideal_map_fraction_1-(6097/10000)^1.1': round(1 - float(1 - P_MAIN) ** 1.1, 4), 'required_fraction': 0.3903,
            'control_plain_birthday_mu_hat': round(ctl / cn, 3), 'control_trials': cn,
            'merged_trails_followed': merged, 'merged_trails_reaching_the_same_DP': coalesced,
            'escape_rate_per_attempt': esc / max(1, att), 'predicted_two_cycle_rate_1_over_12R': 1 / (12 * R),
            'dead_trails': dead, 'ideal_birthday_value': 1.0,
            'scope': 'exponent-domain emulation of the identical walk; keys are pseudorandom hashes of the 12 conjugates'}


def main():
    cfg = json.load(sys.stdin)
    mode = cfg['mode']
    if mode == 'ledger': out = ledger()
    elif mode == 'algebra': out = algebra(cfg['elements'], cfg['seed'])
    elif mode == 'toy': out = toy(cfg['bits'], cfg['trials'], cfg['seed'])
    elif mode == 'sim': out = sim(cfg['bits'], cfg['R'], cfg['trials'], cfg['seed'], db=cfg.get('db', 5))
    elif mode == 'all':
        out = {'ledger': ledger(),
               'algebra': algebra(cfg['algebra_elements'], cfg['seed']),
               'toy': toy(cfg['toy_bits'], cfg['toy_trials'], cfg['seed']),
               'sim': sim(cfg['sim_bits'], cfg['sim_R'], cfg['sim_trials'], cfg['seed'], db=cfg['sim_db'])}
    else: raise SystemExit('unknown mode')
    print(json.dumps(out, indent=1))

if __name__ == '__main__':
    main()
