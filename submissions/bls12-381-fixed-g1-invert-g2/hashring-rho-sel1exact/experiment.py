"""Class-folded distinguished-point rho in GT (BLS12-381, fixed-G1 pairing inversion).

Standard library only; one JSON object in, one out; modes ledger, algebra, toy and sim (`all` runs the four).
Nothing here is a production-scale run; see proof.md section 9 for each mode's scope.
"""
import itertools
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction

M64 = (1 << 64) - 1
P_FP = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R_ORD = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
MU = Fraction(11, 10)                                    # allowed walk inefficiency (H1), unchanged
P_MAIN = Fraction(3903, 10000)                           # main-event probability at S_max under H1: 0.39 plus the loss bound

def ln_bounds(p, d, terms=64):
    """Exact rational bounds on ln(d/(d-p)), d > p > 0: the partial sum of sum x^k/k, x = p/d, is a lower bound and the
    geometric tail an upper one. No floats."""
    x = Fraction(p, d)
    s = Fraction(0)
    xk = Fraction(1)
    for k in range(1, terms + 1):
        xk *= x
        s += xk / k
    xk *= x
    return s, s + xk / ((terms + 1) * (1 - x))

LN_MAIN_LO, LN_MAIN_UP = ln_bounds(P_MAIN.numerator, P_MAIN.denominator)  # ln(10000/6097) = 0.494788, bracketed;
# the ledger's budget uses the upper bound U = LN_MAIN_UP, with 0 < U - ln(10000/6097) < 2^-99 (checked in ledger())
T_CLAIM = float(MU) * math.log(1 / (1 - float(P_MAIN)))  # mu * ln(10000/6097) = 0.544267 (U is within 2^-99 of it)
# Per-slot work other than the reference multiplication (proof.md section 5): every row is the length of a written-out
# AArch64 sequence. Rows marked EXECUTED are assembled and run on class A64 below; the `alloc` rows are priced from the
# sequences proof.md writes out, since they are pointer arithmetic and single loads, not kernels.
SLOT_ITEMS = (
    ('selector: 3 blst_fp_mul at 383 (a*zeta, zeta*d, g0*(a-b)); the organizer kernel price, priced not executed', 3 * 383),
    ('alloc: selector: 3 modular subtractions at 37 (12 LDR, SUBS/SBCS, CSETM, 6 AND, ADDS + 5 ADCS, 6 STR) + 1 negation at '
     '19 (6 SUBS/SBCS, CSETM, 6 AND, ADDS + 5 ADCS, no store) + 18: three blocks, one six-limb load of p each; EXECUTED per '
     'block', 3 * 37 + 19 + 18),
    ('selector: one-pass min and max, 20 + 5*33, each value paying its MOVZ index supply, and the eleventh comparison '
     '(6 CSEL + 1 CSEL + 1 ADD = 14); EXECUTED, 224 instructions with the negation', 20 + 5 * 33 + 14),
    ('alloc: khash: two words over six limbs held in registers, 2*6*9, plus the two mix(salt) seeds and the two '
     'splitmix64 multiplicands (no callee-saved register holds them across a slot); priced, not executed', 2 * 6 * 9 + 4),
    ('distinguished-point test: AND x3, x2, #0xFFFF ; ORR x3, x3, x1 ; CBNZ x3, 80 bits; EXECUTED', 3),
    ('look-back membership: 24 unrolled LDR/CMP/B.EQ slots = 72 plus at most one 20-instruction hit path (key address, '
     'six limbs, 6 EOR + 5 ORR + CBZ, one branch); EXECUTED, 72 a miss and 92 a hit', 24 * 3 + 20),
    ('look-back ring write: three shifted adds, three STP of the key limbs, ADD x11, sp, p, LSL #3 ; STR x1, [x11]; EXECUTED', 3 + 3 + 2),
    ('alloc: bookkeeping: accepted-step counter with its 126-bit bound loaded (2 LDR + 5), the ring-slot update '
     'ADD/CMP/CSEL (3, EXECUTED) and two branches', 7 + 3 + 2),
    ('alloc: table entry address and the dense kernel\'s three arguments: LDR h, ADD h+a, AND (h+a) mod R, MOVZ 6912, '
     'MADD, LDR k_n, MOVZ 576, MADD, MOV x1,y, EOR x0,y,#0x400', 10),
    ('alloc: accept tail: AND h\',w1,#0x3FFF ; STR h\' ; STR k_n ; MOVZ a,#0 ; EOR y,y,#0x400', 5),
    ('alloc: probe ladder on the reject path: ADD a,#1 ; CMP a,#8 ; B.LT', 3),
    ('alloc: selector operand base, twice: EOR x29, x20, #0x400, after the dense kernel and after the fp_mul kernels', 2),
    ('alloc: three blst_fp_mul argument setups, three ADD each (the kernel clobbers x0-x18)', 9),
    ('alloc: trail slot counter with its 126-bit bound loaded (2 LDR + 5)', 7),
    ('alloc: global slot counter with its 126-bit bound loaded (2 LDR + 5)', 7),
    ('alloc: next-slot branch', 1),
)
C_ITEM = sum(n for _, n in SLOT_ITEMS)                   # itemised per-slot instructions
C_MAX_SLOT = C_ITEM - 3                                  # over slot types (no ladder on accept)
C_SLOT = 1773                                            # charged: the ceiling, C_MAX_SLOT + 9
assert C_ITEM == 1767, C_ITEM
assert C_MAX_SLOT == 1764, C_MAX_SLOT
SEL_ITEM = sum(n for _, n in SLOT_ITEMS[:3])             # the selector's three rows, quoted in proof.md section 3
assert SEL_ITEM == 1149 + 148 + 199 == 1496, SEL_ITEM
assert C_ITEM - SEL_ITEM == sum(n for _, n in SLOT_ITEMS[3:])
assert C_SLOT - C_MAX_SLOT == 9, "explicit margin (9/1764 = 0.5102%), as proof.md section 5 states"
assert Fraction(C_SLOT - C_MAX_SLOT, C_MAX_SLOT) >= Fraction(5, 1000), "at least 0.5% of it"
W_WINDOW = 24                                            # look-back window slots (section 4.1)

# ---------------------------------------------------------------- hashing of keys
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

# --------------------------------------------- the priced AArch64 lists, executed on a register file
class A64:
    """A minimal AArch64 subset machine: each priced list is run here, so a row is a count of instructions, not a
    constant. Flags are architectural: a SUBS/SBCS chain's carry is the whole subtraction's borrow, Z is the last
    limb's only, ADC does not carry."""
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
        return 0 if r == 'xzr' else (self.sp if r == 'sp' else self.x[int(r[1:])])

    def put(self, r, v):
        v &= self.M64
        if r == 'xzr': return
        if r == 'sp': self.sp = v
        else: self.x[int(r[1:])] = v

    def word(self, addr): return self.mem.get(addr & self.M64, 0)

    def setword(self, addr, v): self.mem[addr & self.M64] = v & self.M64

    def add(self, a, b, cin=0, flags=False):
        full = a + b + cin
        res = full & self.M64
        if flags:
            self.N = (res >> 63) & 1
            self.Z = 1 if res == 0 else 0
            self.C = 1 if full > self.M64 else 0
        return res

    def sub(self, a, b, cin=1, flags=False):
        full = a - b - (1 - cin)                    # the carry is the not-borrow
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
        elif op == 'CSEL': self.put(a[0], self.get(a[1]) if self.cond(a[3]) else self.get(a[2]))
        elif op == 'CSETM': self.put(a[0], self.M64 if self.cond(a[1]) else 0)
        elif op == 'LDR': self.put(a[0], self.word(self.addr(a[1:])))
        elif op == 'STR': self.setword(self.addr(a[1:]), self.get(a[0]))
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
        prog.append((parts[0], parts[1:]))
    return prog, labels


K_OFF, V_OFF, P_OFF = 192, 1344, 1632       # frame: H[24] at 0, K[24] at 192, the six values at 1344, p at 1632
Z_OFF, G_OFF = 1680, 1728                   # the two tower constants zeta and g0 (48 bytes each)
H_OFF, KS_OFF, ST_OFF = 1776, 1784, 1792    # the state's table index h, its twist index k_n, the two strides (6912, 576)
S_OFF, KM_OFF, W_OFF = 1808, 1824, 1888     # the two accumulator seeds, the two splitmix64 multiplicands, the a-b scratch
YA, YB, YX = 2048, 3072, 0x400              # the two y buffers, 0x400 apart with bit 10 clear, so `EOR xd, x20, #0x400` toggles

def limbs(x, k=6):
    """The six 64-bit limbs of a non-negative integer (canonical representative, low limb first)."""
    return [(x >> (64 * i)) & M64 for i in range(k)]

def dp_test_prog():
    """The distinguished-point test (3): `AND x3, x2, #0xFFFF ; ORR x3, x3, x1 ; CBNZ x3`, 64-bit, so 80 bits; the
    32-bit form is a 48-bit predicate (`dp_test_check` witnesses it with `w_1 = 2^32`)."""
    return asm("""
        AND  x3, x2, #0xFFFF        ; the low 16 bits of w_2
        ORR  x3, x3, x1             ; the whole of w_1
        CBNZ x3, not_dp
    not_dp:
    """)

def seed_prog(w):
    """Seeding both rings (1 + 5w): `MOVZ p, #0` and, per slot, the key-slot `ADD`, three `STP` and `STR x1, [sp, #8s]`
    (a 64-bit store whose scaled immediate `8s <= 184` is encodable)."""
    lines = ['MOVZ x27, #0']
    for s in range(w):
        lines += ['ADD x10, sp, #%d' % (K_OFF + 48 * s),
                  'STP x4, x5, [x10]', 'STP x6, x7, [x10, #16]', 'STP x8, x9, [x10, #32]',
                  'STR x1, [sp, #%d]' % (8 * s)]
    return asm('\n'.join(lines))

def scan_prog(w):
    """The membership scan (3w) and its hit stub (20): `LDR x10, [sp, #8s]`, `CMP x10, x1`, `B.EQ hit_s`; the stub is an
    address, 6 `LDR`, 6 `EOR` folded by 5 `ORR`, `CBZ` and one branch. A miss falls off the last `B.EQ` at 3*24 = 72."""
    lines = []
    for s in range(w):
        lines += ['hit_%d:' % s, 'ADD x11, sp, #%d' % (K_OFF + 48 * s),
                  'LDR x12, [x11]', 'EOR x12, x4, x12']
        for i in range(1, 6):
            lines += ['LDR x13, [x11, #%d]' % (8 * i), 'EOR x13, x%d, x13' % (4 + i), 'ORR x12, x12, x13']
        lines += ['CBZ x12, equal_%d' % s, 'B end']
    for s in range(w):
        lines += ['equal_%d:' % s, 'B end']
    lines.append('scan:')
    for s in range(w):
        lines += ['LDR x10, [sp, #%d]' % (8 * s), 'CMP x10, x1', 'B.EQ hit_%d' % s]
    lines.append('end:')
    return asm('\n'.join(lines))

def write_prog():
    """The ring write (8): three shifted adds for the key slot, three `STP`, then `ADD x11, sp, p, LSL #3` and
    `STR x1, [x11]` — `STR w1, [sp, p, LSL #3]` is not encodable and would store half a word."""
    return asm("""
        ADD x10, sp, #192
        ADD x10, x10, x27, LSL #4
        ADD x10, x10, x27, LSL #5
        STP x4, x5, [x10]
        STP x6, x7, [x10, #16]
        STP x8, x9, [x10, #32]
        ADD x11, sp, x27, LSL #3
        STR x1, [x11]
    """)

def slot_update_prog(w):
    """The bounded ring-slot update (3): `ADD p, p, #1 ; CMP p, #w ; CSEL p, p, xzr, NE`."""
    return asm("""
        ADD  x27, x27, #1
        CMP  x27, #%d
        CSEL x27, x27, xzr, NE
    """ % w)

def selector_tail_prog():
    """The priced selector tail: 20 + 5*33 + 25 + 14 = 224. Minimum in x4-x9, maximum in x10-x15 (overwritten in place by
    its negation), candidate in x16, x17, x0-x3, indices in x18 and x29, index supply in x30: the six value `CSEL`s of
    the eleventh comparison write the key into x4-x9, the set the khash, the scan and the ring write read, with no move
    between. Every chain writes xzr and only its carry is read (`cand - min` borrows iff `cand < min`, `nb - m` leaves
    C = 1 iff `m <= nb`); Z is never read."""
    cand = (16, 17, 0, 1, 2, 3)                  # one candidate's six limbs
    lines = []
    for i in range(6):
        lines.append('LDR x%d, [sp, #%d]' % (cand[i], V_OFF + 8 * i))
    for i in range(6):
        lines.append('MOV x%d, x%d' % (4 + i, cand[i]))
    for i in range(6):
        lines.append('MOV x%d, x%d' % (10 + i, cand[i]))
    lines += ['MOVZ x18, #0', 'MOVZ x29, #0']
    for k in range(1, 6):
        lines.append('MOVZ x30, #%d' % k)
        for i in range(6):
            lines.append('LDR x%d, [sp, #%d]' % (cand[i], V_OFF + 48 * k + 8 * i))
        lines.append('SUBS xzr, x%d, x4' % cand[0])
        for i in range(1, 6):
            lines.append('SBCS xzr, x%d, x%d' % (cand[i], 4 + i))
        for i in range(6):
            lines.append('CSEL x%d, x%d, x%d, LO' % (4 + i, cand[i], 4 + i))
        lines.append('CSEL x18, x30, x18, LO')
        lines.append('SUBS xzr, x10, x%d' % cand[0])
        for i in range(1, 6):
            lines.append('SBCS xzr, x%d, x%d' % (10 + i, cand[i]))
        for i in range(6):
            lines.append('CSEL x%d, x%d, x%d, LO' % (10 + i, cand[i], 10 + i))
        lines.append('CSEL x29, x30, x29, LO')
    lines.append('SUBS x10, xzr, x10')
    for i in range(1, 6):
        lines.append('SBCS x%d, xzr, x%d' % (10 + i, 10 + i))
    lines.append('CSETM x16, CC')
    for i in range(6):
        lines += ['LDR x17, [sp, #%d]' % (P_OFF + 8 * i), 'AND x17, x17, x16',
                  '%s x%d, x%d, x17' % ('ADDS' if i == 0 else 'ADCS', 10 + i, 10 + i)]
    lines.append('SUBS xzr, x10, x4')
    for i in range(1, 6):
        lines.append('SBCS xzr, x%d, x%d' % (10 + i, 4 + i))
    lines.append('ADD x29, x29, #6')
    for i in range(6):
        lines.append('CSEL x%d, x%d, x%d, CS' % (4 + i, 4 + i, 10 + i))
    lines.append('CSEL x18, x18, x29, CS')
    return asm('\n'.join(lines))

def modsub_prog():
    """One modular subtraction x - y mod p (37) after the block's six loads of p (the alloc row's 6): 12 `LDR`, the
    chain, `CSETM`, 6 `AND` with one `ADDS` and five `ADCS`, 6 `STR`."""
    lines = []
    for i in range(6):
        lines.append('LDR x%d, [sp, #%d]' % (20 + i, P_OFF + 8 * i))
    for i in range(6):
        lines.append('LDR x%d, [sp, #%d]' % (i, V_OFF + 8 * i))
    for i in range(6):
        lines.append('LDR x%d, [sp, #%d]' % (6 + i, V_OFF + 48 + 8 * i))
    lines.append('SUBS x12, x0, x6')
    for i in range(1, 6):
        lines.append('SBCS x%d, x%d, x%d' % (12 + i, i, 6 + i))
    lines.append('CSETM x19, CC')
    for i in range(6):
        lines += ['AND x%d, x%d, x19' % (20 + i, 20 + i),
                  '%s x%d, x%d, x%d' % ('ADDS' if i == 0 else 'ADCS', 12 + i, 12 + i, 20 + i)]
    for i in range(6):
        lines.append('STR x%d, [sp, #%d]' % (12 + i, V_OFF + 96 + 8 * i))
    return asm('\n'.join(lines))

# --------------------------------------- the section 5 window sequences and their executed checks
RING_SLOTS = W_WINDOW                                    # 24 stored 64-bit hashes and 24 stored 48-byte keys

def ring_slot_update(p, w=W_WINDOW):
    """The reference for the executed slot update: `(p+1) mod w`, three register operands only (a mask is not a
    reduction modulo 24 and is not used)."""
    p = p + 1
    if p == w: p = 0
    return p

# The priced look-back membership of proof.md section 5: 24 unrolled hash slots plus at most one hit path.
SCAN_SLOT = 3                    # LDR the 8-byte stored hash ; CMP that hash with w1 in x1 ; B.EQ this slot's stub
SCAN_ROW = RING_SLOTS * SCAN_SLOT                        # 72: all 24 hashes compared and no hash matched
HIT_ROW = 1 + 6 + 6 + 5 + 1 + 1                          # 20: key-slot address, six LDR, 6 EOR, 5 ORR, CBZ, one B
MEMB_ROW = SCAN_ROW + HIT_ROW                            # 92: a miss is 72, a hit at slot 23 is 92 (the maximum)
assert (SCAN_ROW, HIT_ROW, MEMB_ROW) == (72, 20, 92)
assert SLOT_ITEMS[5][1] == MEMB_ROW, "the membership row is the sequence scan_prog executes"

def key_equal(a, b, k=6):
    """The exact key test, as the reference: six `EOR` folded by five `ORR`; the keys are canonical and below p, so the
    residues vanish exactly when they are equal."""
    acc = 0
    for i in range(k):
        acc |= (a[i] ^ b[i])
    return acc == 0

def ring_frame(H, K, w):
    """A machine whose frame holds the two rings: H[s] at 8s, K[s] at K_OFF + 48s."""
    m = A64()
    for s in range(w):
        m.setword(m.sp + 8 * s, H[s])
        kk = K[s] if K[s] is not None else (0,) * 6
        for i in range(6):
            m.setword(m.sp + K_OFF + 48 * s + 8 * i, kk[i])
    return m

def hash_scan(H, K, w1, kk, w=RING_SLOTS):
    """The membership scan, executed as `scan_prog` builds it: `w_1` in x1, the candidate's limbs in x4..x9, at most one
    stub, at most 3w + 20 = 92. Returns (rejected, executed, hit slot, equality)."""
    m = ring_frame(H, K, w)
    m.put('x1', w1)
    for i, l in enumerate(kk[:6]): m.put('x%d' % (4 + i), l)
    prog, labels = scan_prog(w)
    instr = m.run(prog, labels, 'scan')
    hit = [t for t in m.taken if t.startswith('hit_') or t.startswith('equal_')]
    if not hit: return False, instr, -1, False
    return True, instr, int(hit[0].split('_')[1]), hit[-1].startswith('equal_')

def ring_write(H, K, p, kk, w1):
    """The ring write, executed as `write_prog` builds it: 8 instructions, the whole 64-bit word stored."""
    m = ring_frame(H, K, len(H))
    m.put('x27', p)
    m.put('x1', w1)
    for i, l in enumerate(kk[:6]): m.put('x%d' % (4 + i), l)
    prog, labels = write_prog()
    n = m.run(prog, labels)
    assert n == 8, n
    K[p] = tuple(m.word(m.sp + K_OFF + 48 * p + 8 * i) for i in range(6))
    H[p] = m.word(m.sp + 8 * p)
    assert H[p] == w1 and K[p] == tuple(kk[:6]), 'the ring write did not store the key and the word'
    return n

RING_SEED_ROW = RING_SLOTS * (1 + 3 + 1) + 1             # 121: per slot address + three STP + hash store, and p = 0

class HashRing:
    """The production window on toy keys: a ring of keys (state k at slot k mod 24) and a ring of their words, seeded
    with `kappa_0` in all slots (121), so the window is the last `min(k+1, 24)` keys."""
    def __init__(self, w, hashword):
        self.w = w
        self.H = [0] * w
        self.K = [None] * w
        self.p = 0
        self.hashword = hashword

    def seed(self, key):
        """Execute the priced seeding list: `kappa_0` in all w key slots, its word in all w hash slots, p = 0."""
        h = self.hashword(key)
        kk = tuple(limbs(key))
        m = A64()
        m.put('x1', h)
        for i, l in enumerate(kk): m.put('x%d' % (4 + i), l)
        prog, labels = seed_prog(self.w)
        n = m.run(prog, labels)
        for s in range(self.w):
            self.K[s] = tuple(m.word(m.sp + K_OFF + 48 * s + 8 * i) for i in range(6))
            self.H[s] = m.word(m.sp + 8 * s)
        self.p = 0
        assert n == 1 + 5 * self.w, n
        return n

    def keys(self):
        """The stored keys in slot order, as integers (the set they form is the window)."""
        return [sum(l << (64 * i) for i, l in enumerate(kk)) for kk in self.K]

    def member(self, key):
        rejected, _instr, _slot, _eq = hash_scan(self.H, self.K, self.hashword(key), limbs(key), self.w)
        return rejected

    def insert(self, key):
        self.p = ring_slot_update(self.p, self.w)
        return ring_write(self.H, self.K, self.p, limbs(key), self.hashword(key))

def dp_test_check(seed=20261025):
    """Execute the distinguished-point test and check it is the 80-bit predicate. `w_1 = 2^32`, `w_2 = 0` is the witness
    against the 32-bit form, whose `S_max*2^-48 = 2^77.7` reports are far above the `2^47` cap."""
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
        got = 'not_dp' not in m.taken                       # CBNZ taken means the state is not distinguished
        assert got == (w1 == 0 and (w2 & 0xFFFF) == 0), (w1, w2)
        distinguished += got
    m = A64(); m.put('x1', 1 << 32); m.put('x2', 0)
    m.run(prog, labels)
    witness_64 = 'not_dp' not in m.taken                    # False: not distinguished
    witness_48 = ((1 << 32) % (1 << 32) == 0) and (0 & 0xFFFF) == 0
    assert witness_64 is False and witness_48 is True, 'the width witness does not separate the two forms'
    return {'distinguished_point_test_instructions': 3, 'distinguished_point_test_bits': 80,
            'distinguished_point_test_cases_checked': len(cases), 'distinguished_point_test_hits': distinguished,
            'distinguished_point_test_64_bit_form_says_dp_for_w1_2^32_w2_0': witness_64,
            'distinguished_point_test_48_bit_form_says_dp_for_w1_2^32_w2_0': witness_48}

def ring_check(accepts=200, pairs=600, seed=20261024):
    """Identity checks for the window sequences, each executed as its priced list: the slot update, the key test (3,412
    pairs), the seed (121), the scan over every fill 0..24 (72 a miss, 92 the maximum), the write (8), `dp_test_check`."""
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
    P = P_FP
    cases = [(0, 0), (0, 1), (1, 0), (P - 1, P - 1), (P - 1, P - 2), (P - 2, P - 1), (P - 1, 0), (0, P - 1),
             (1 << 380, (1 << 380) + 1), ((1 << 320) - 1, 1 << 320), ((1 << 64) - 1, 1 << 64),
             (0x1234567890abcdef, 0x1234567890abcdef)]
    for _ in range(pairs):
        x = rng.randrange(P)
        y = rng.randrange(P)
        cases += [(x, x), (x, y), (x, (x + 1) % P), (x, (x - 1) % P), (x & ~M64, y & ~M64)]
    checked = 0
    for x, y in cases:
        assert key_equal(limbs(x), limbs(y)) == (x == y), ('equality', x, y)
        checked += 1
    toy_cases = 0
    for _ in range(200):
        q = rng.randrange(1 << 60)
        x = rng.randrange(q)
        y = rng.randrange(q)
        for a, b in ((x, x), (x, y)):
            assert key_equal(limbs(a), limbs(b)) == (a == b), ('toy equality', a, b)
            toy_cases += 1
    hw = lambda kk: khash(kk, 0x5EED0F24, 6)          # the window hash word of this toy ring
    scanned = hits = 0
    max_instr = 0
    max_slot = -1
    miss_instr = None
    collisions = 0
    seed_row = None
    for fill in range(0, w + 1):                         # 0 keys (the seeded ring only) .. 24 accepts
        keys = [rng.randrange(1 << 40) for _ in range(fill)]
        ring = HashRing(w, hw)
        seed_row = ring.seed(0xAAAA)
        assert seed_row == RING_SEED_ROW == 121
        assert all(k == 0xAAAA for k in ring.keys()) and all(ring.H[s] == hw(0xAAAA) for s in range(w))
        for kk in keys:
            assert ring.insert(kk) == 8, 'the ring write is not 8 instructions'
        window = ([0xAAAA] + keys)[max(0, len(keys) + 1 - w):]      # the last min(k+1, 24) states
        assert sorted(set(ring.keys())) == sorted(set(window)), ('the ring does not hold the window', fill)
        candidates = list(ring.keys()) + [0xBBBB, (1 << 40) + 7, rng.randrange(1 << 40)]
        present = set(window)
        for cand in candidates:
            rejected, instr, slot, eq = hash_scan(ring.H, ring.K, hw(cand), limbs(cand), w)
            assert rejected == (cand in present), ('the hash ring and the plain set disagree', cand)
            assert rejected == eq or not eq, 'a rejection without a hash hit is impossible'
            if rejected:
                assert eq == (cand in present), ('the hit path is not exact', cand, fill)
                hits += 1
                if slot > max_slot or instr > max_instr:
                    max_slot, max_instr = slot, instr
            else:
                miss_instr = instr if miss_instr is None else min(miss_instr, instr)
            scanned += 1
    assert max_instr == MEMB_ROW == 92 and max_slot == w - 1, ('the maximum is a hit at the last slot', max_instr, max_slot)
    assert miss_instr == SCAN_ROW == 72, ('a miss is 72 instructions', miss_instr)
    ring = HashRing(w, hw)
    ring.seed(0xAAAA)
    keys = [rng.randrange(1 << 40) for _ in range(w)]
    for kk in keys:
        ring.insert(kk)
    assert ring.p == 0, 'after 24 accepts the slot counter is back at 0'
    top = keys[w - 2]                                   # the key written into slot 23
    assert sum(l << (64 * i) for i, l in enumerate(ring.K[w - 1])) == top
    rejected, instr, slot, eq = hash_scan(ring.H, ring.K, hw(top), limbs(top), w)
    assert (rejected, eq, instr, slot) == (True, True, MEMB_ROW, w - 1), ('the 92-instruction maximum', instr, slot)
    ring.H[w - 1] = hw(0xCCCC)                          # a stored word matching, a stored key that does not
    rejected, instr, slot, eq = hash_scan(ring.H, ring.K, hw(0xCCCC), limbs(0xCCCC), w)
    assert (rejected, eq, instr, slot) == (True, False, MEMB_ROW, w - 1), ('the collision path', rejected, eq)
    assert 0xCCCC not in set(ring.keys()), 'the collision case must not be a stored key'
    collisions += 1
    return {'ring_slot_accepts_checked': accepts, 'ring_slot_assertions': 3 * accepts,
            'ring_slot_update_instructions': 3, 'ring_slot_cycle_length': w,
            'key_equality_pairs_checked': checked + toy_cases,
            'ring_seed_instructions': seed_row, 'ring_write_instructions': 8,
            'membership_scan_candidates_checked': scanned + 1, 'membership_scan_hits_checked': hits,
            'membership_scan_collisions_exercised': collisions,
            'membership_scan_miss_instructions': miss_instr, 'membership_scan_max_instructions': max_instr,
            'membership_scan_max_hit_slot': max_slot,
            'membership_instructions_executed_bound': MEMB_ROW,
            'membership_instructions_charged': SLOT_ITEMS[5][1],
            'distinguished_point_test': dp_test_check(seed)}

# ---------------------------------------------------------------- tower Fp2 / Fp12
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
                b = t[k + 6]                       # multiply by xi = 1 + u
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
    def frob(x):                                    # x -> x^p
        return tuple(f2mul(f2conj(x[i]), gam[i]) for i in range(6))
    zeta = f2mul(gamma, f2conj(gamma))
    assert zeta[1] == 0 and (zeta[0] * zeta[0] - zeta[0] + 1) % p == 0
    # gamma = g0*(1-u) exactly (proof.md section 3), which makes the selector's d one multiplication, g0*(a-b).
    assert (gamma[0] + gamma[1]) % p == 0, "gamma is not g0*(1-u) in this tower"
    assert (2 * gamma[0] * gamma[0] - zeta[0]) % p == 0, "zeta = 2*g0^2"
    return dict(p=p, mul=f12mul, pow=f12pow, one=one12, frob=frob, gamma=gamma, zeta=zeta[0],
                f2mul=f2mul, f2pow=f2pow, nlimbs=(p.bit_length() + 63) // 64)

def frobk(F, x, k):
    for _ in range(k % 12): x = F['frob'](x)
    return x

def sel12(F, y):
    """The twelve first-coordinates of c1(sigma^k(y)) in 3 Fp multiplications: V[2j] = zeta^j*a, V[2j+1] = zeta^j*d,
    V[k+6] = -V[k], d = g0*(a-b) since gamma = g0*(1-u) exactly."""
    p = F['p']; zeta = F['zeta']; g0 = F['gamma'][0]
    a, b = y[1]
    t = (a - b) % p                        # add/sub 1
    d = g0 * t % p                         # mult 1
    za = zeta * a % p                      # mult 2
    zd = zeta * d % p                      # mult 3
    h2 = (za - a) % p; hd = (zd - d) % p   # zeta^2 = zeta - 1, add/sub 2 and 3
    return (a, d, za, zd, h2, hd, (-a) % p, (-d) % p, (-za) % p, (-zd) % p, (-h2) % p, (-hd) % p)

def sel(F, y):
    """Class selector: 3 Fp multiplications; returns (kappa, kstar), the minimum over the twelve first coordinates of
    c1(sigma^k(y)) and the smallest k attaining it (proof.md section 5's priced sequence). A minimiser among the six
    stored values has index <= 5 and the negated half carries i_max + 6 in 6..11, so a tie goes to the stored half."""
    p = F['p']; zeta = F['zeta']; g0 = F['gamma'][0]
    a, b = y[1]
    t = (a - b) % p                        # add/sub 1
    d = g0 * t % p                         # mult 1
    za = zeta * a % p                      # mult 2
    zd = zeta * d % p                      # mult 3
    h2 = (za - a) % p; hd = (zd - d) % p   # zeta^2 = zeta - 1, add/sub 2 and 3
    V = (a, d, za, zd, h2, hd)
    vi = vx = 0                            # index of the minimum / of the maximum over the six
    for k in range(1, 6):                  # 10 comparisons with conditional move of value and index
        if V[k] < V[vi]: vi = k
        if V[k] > V[vx]: vx = k
    neg = (p - V[vx]) % p                  # one modular negation, canonical (M = 0 gives 0, which wins only
                                           # when all six are 0, i.e. never on a well-formed input)
    if V[vi] <= neg:                       # the eleventh comparison; on a tie the smaller index wins
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
    """The same twelve values from the general form d = g0*a + g1*b: the cross-check that the 3-multiplication form is
    the same quantity. Not charged."""
    p = F['p']; zeta = F['zeta']; g0, g1 = F['gamma']
    a, b = y[1]; za = zeta * a % p; d = (g0 * a + g1 * b) % p; zd = zeta * d % p
    h2 = (za - a) % p; hd = (zd - d) % p
    return (a, d, za, zd, h2, hd, (-a) % p, (-d) % p, (-za) % p, (-zd) % p, (-h2) % p, (-hd) % p)

# ---------------------------------------------------------------- exact ledger
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
    """The reference for the executed eleventh comparison: `(min, i_min)` when `min <= nb` (`nb = p - max`), else
    `(nb, i_max + 6)`; a tie goes to the stored half."""
    vi, vx = sel_extrema(V)
    nb = (p - V[vx]) % p
    return (V[vi], vi) if V[vi] <= nb else (nb, vx + 6)

def sel_direct12(V, p):
    """The definition: the smallest index attaining the minimum of the twelve values V_k and (p - V_k) mod p."""
    twelve = tuple(V) + tuple((p - v) % p for v in V)
    best, kb = twelve[0], 0
    for k in range(1, 12):
        if twelve[k] < best: best, kb = twelve[k], k
    return best, kb

def a64_selector(V, p):
    """Execute the priced tail on one six-tuple: (executed, kappa, k*), the key in x4-x9 and k* in x18, the set the scan
    and the ring write read."""
    m = A64()
    for k, v in enumerate(V):
        for i, l in enumerate(limbs(v)):
            m.setword(m.sp + V_OFF + 48 * k + 8 * i, l)
    for i, l in enumerate(limbs(p)):
        m.setword(m.sp + P_OFF + 8 * i, l)
    prog, labels = selector_tail_prog()
    n = m.run(prog, labels)
    kappa = sum(m.get('x%d' % (4 + i)) << (64 * i) for i in range(6))
    return n, kappa, m.get('x18')

def selector_min_check(tuples_=240, pairs=40, seed=20261025):
    """The two-half minimum equals the direct 12-value scan, and the lists are executed: the 7^6 identity; the tail (224)
    on 245 small-field tuples (51 ties) and 40 production tuples; the 37-instruction subtraction with its six loads."""
    ps = 101
    vs = (0, 1, 2, 3, 50, 99, 100)
    n = 0
    for base in itertools.product(vs, repeat=6):
        # the one-pass network the ledger prices: same extrema and smallest indices as two separate scans
        assert sel_extrema(base) == (min(range(6), key=lambda k: base[k]),
                                     max(range(6), key=lambda k: base[k])), base
        assert sel_two_half(base, ps) == sel_direct12(base, ps), base
        n += 1
    rng = random.Random(seed)
    # (b) the priced tail, executed on sampled tuples of the small field (ties and zeros included)
    sample = [(0, 0, 0, 0, 0, 0), (ps - 1, 1, 1, 1, 1, 1), (ps - 1, 0, 0, 0, 0, 0),
              (1, ps - 1, 2, 3, 4, 5), (50, 51, 50, 51, 50, 51)]
    for _ in range(tuples_):
        sample.append(tuple(rng.choice(vs) for _ in range(6)))
    tail_row = 20 + 5 * 33 + 14
    ties = zeros = 0
    printed = []
    for V in sample:
        n_ex, kappa, kstar = a64_selector(V, ps)
        assert n_ex == 224, ('the tail is not 224 instructions', n_ex)
        assert (kappa, kstar) == sel_two_half(V, ps) == sel_direct12(V, ps), ('the executed tail', V)
        vi, vx = sel_extrema(V)
        if V[vi] == (ps - V[vx]) % ps:
            ties += 1
            if len(printed) < 2: printed.append(('tie', V, kappa, kstar))
        if all(v == 0 for v in V):
            zeros += 1
            printed.append(('all-zero', V, kappa, kstar))
    assert ties > 0, 'the executed eleventh comparison never met a tie between the halves'
    assert zeros > 0, 'the executed negation never met M = 0'
    # (c) production tuples: the six values the three fp_mul kernels produce, on the real field and constants
    F = make_field(P_FP)
    p = P_FP
    zeta = F['zeta']; g0 = F['gamma'][0]
    for _ in range(pairs):
        a = rng.randrange(p); b = rng.randrange(p)
        d = g0 * ((a - b) % p) % p                       # the kernel g0*(a-b)
        za = zeta * a % p; zd = zeta * d % p             # the kernels zeta*a, zeta*d
        V = (a, d, za, zd, (za - a) % p, (zd - d) % p)   # and zeta^2 = zeta - 1 on each
        n_ex, kappa, kstar = a64_selector(V, p)
        assert n_ex == 224, n_ex
        assert (kappa, kstar) == sel_two_half(V, p) == sel_direct12(V, p), (a, b)
        if len(printed) < 4: printed.append(('production', V, kappa, kstar))
    # (d) one modular subtraction with its block's modulus load (the 148 row is three blocks), executed
    sub_exec = set()
    for pp in (p, ps):
        for _ in range(max(1, pairs // 2)):
            x = rng.randrange(pp); y = rng.randrange(pp)
            m = A64()
            for i, l in enumerate(limbs(x)): m.setword(m.sp + V_OFF + 8 * i, l)
            for i, l in enumerate(limbs(y)): m.setword(m.sp + V_OFF + 48 + 8 * i, l)
            for i, l in enumerate(limbs(pp)): m.setword(m.sp + P_OFF + 8 * i, l)
            prog, labels = modsub_prog()
            sub_exec.add(m.run(prog, labels))
            r = sum(m.word(m.sp + V_OFF + 96 + 8 * i) << (64 * i) for i in range(6))
            assert r == (x - y) % pp, ('the executed subtraction', x, y, pp)
    assert sub_exec == {43}, ('the subtraction is 37 + 6 loads of p', sub_exec)
    return {'exhaustive_tuples': n,
            'small_field_tuples_executed': len(sample),
            'production_tuples_executed': pairs,
            'tied_halves_exercised': ties,
            'all_zero_tuples_exercised': zeros,
            'tail_instructions_executed': 224,
            'extrema_row_instructions': 20 + 5 * 33,
            'negation_instructions': 25,
            'eleventh_comparison_instructions': 14,
            'tail_instructions_priced': tail_row,
            'modular_subtraction_executed': 43,
            'modular_subtraction_row': 37,
            'modulus_load_instructions': 6,
            'tuples_printed': [{'kind': k, 'V': [hex(v) for v in V], 'kappa': hex(kappa), 'kstar': ks}
                               for k, V, kappa, ks in printed]}

def ledger():
    getcontext().prec = 80
    N = (R_ORD - 1) // 12
    mu = MU
    # S^2 >= 2 mu N U at the main-event probability 3903/10000, U = LN_MAIN_UP being the exact rational upper bound on
    # L = ln(10000/6097); this U-based S_max is the one the ledger, the claim and proof.md section 5 use
    s2 = 2 * mu * N * LN_MAIN_UP
    S_max = isqrt_ceil(s2)
    assert Fraction(S_max) ** 2 >= s2 > Fraction(S_max - 1) ** 2
    # 0 < U - L < 2^-99 (proof.md section 5) in exact rationals: the 70-term partial sum is below L.
    xt = Fraction(P_MAIN.numerator, P_MAIN.denominator)
    part, xtk = Fraction(0), Fraction(1)
    for k in range(1, 71):
        xtk *= xt
        part += xtk / k
    assert 0 < LN_MAIN_UP - part < Fraction(1, 1 << 99)
    up70 = part + xtk * xt / (71 * (1 - xt))
    assert LN_MAIN_UP > up70 > part, "U is the wider upper bound"
    assert S_max - isqrt_ceil(2 * mu * N * up70) > (1 << 25)
    # the lower Taylor bracket is loose by far less than the budget: the check at 2^-90
    assert S_max - isqrt_ceil(2 * mu * N * LN_MAIN_LO) <= (S_max >> 90)
    M_max = S_max + (S_max >> 13) + 1                    # probes allowed: 2^-13 of the accepted steps
    c_over = C_SLOT                                      # instructions per Fp12 multiplication slot (section 5)
    unit = Fraction(20000 + c_over, 20000)               # 21773/20000 per multiplication slot
    assert unit == Fraction(21773, 20000)
    sel_res = selector_min_check()                       # the priced two-half minimum = the direct scan, executed
    ring_res = ring_check()                              # the window's ring-slot, scan, write, seed and DP sequences
    R = 1 << 14
    theta_bits = 80
    L_cap = 1 << 87
    D_max = 1 << 47                                      # reports cap
    ST_max = 1 << 47                                     # starts cap
    REC_max = 16                                         # recoveries cap
    add = {}
    add['pairings (g and one spare), 2*2^30 instr'] = Fraction(2 * (1 << 30), 20000)
    add['table: R*(1100 + 12) units'] = R * 1112
    add['starts: ST_max*1100 units'] = ST_max * 1100
    add['reports: D_max*2^16 instr'] = Fraction(D_max * (1 << 16), 20000)
    add['replays: REC_max*2*L_cap*(5/4)'] = Fraction(REC_max * 2 * L_cap * 5, 4)
    add['finish last trail: L_cap*unit'] = L_cap * unit
    add['G2 scalar mult + GT check (<= 2^12 units)'] = 1 << 12
    A = sum(add.values())
    W = unit * M_max + A
    # U's total is within a relative 2^-99 of L's (proof.md section 5).
    S_lb = isqrt_ceil(2 * mu * N * part)
    assert S_lb <= S_max and W / (unit * (S_lb + (S_lb >> 13) + 1) + A) - 1 < Fraction(1, 1 << 99)
    # exact log2 with Decimal for display; the claim is found by exact integer comparisons below
    log2W = (Decimal(W.numerator) / Decimal(W.denominator)).ln() / Decimal(2).ln()
    # claim = least hundredth at or above log2 W: the least integer m with W^100 <= 2^m, by binary search on m
    lo, hi = 0, 100 * W.numerator.bit_length()
    while lo < hi:
        mid = (lo + hi) // 2
        if W.numerator ** 100 <= (1 << mid) * W.denominator ** 100:
            hi = mid
        else:
            lo = mid + 1
    c100 = lo
    exact_le = W.numerator ** 100 <= (1 << c100) * W.denominator ** 100           # W^100 <= 2^(100*claim)
    exact_gt_prev = W.numerator ** 100 > (1 << (c100 - 1)) * W.denominator ** 100  # W^100 > 2^(100*claim - 1)
    claim = Fraction(c100, 100)
    pre = add['table: R*(1100 + 12) units'] + add['pairings (g and one spare), 2*2^30 instr']
    mem = {
        'reports_bytes_max': D_max * 96,                # key 48 B + seed 16 B + 3 pointers + colour
        'twist_table_bytes': 12 * R * 576,
        'uv_scalars_bytes': 1 << 20,
        'window_and_state_bytes': 1 << 20,
    }
    peak = sum(mem.values())
    trap_max, trap_sum = trap_bound(R, 25, 512)
    expected_dp = Decimal(S_max) * Decimal(2) ** (-theta_bits)
    # 80 bits against the report cap: S_max*2^-80 = 2^45.697 fits D_max = 2^47; 48 bits give 2^77.697.
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
    ceiling = c_over                                     # the largest charge whose hundredth is still this claim's
    while hundredth(ceiling + 1) == hundredth(c_over): ceiling += 1
    return {
        'log2_r': float(Decimal(R_ORD).ln() / Decimal(2).ln()), 'log2_N_classes': float(Decimal(N).ln() / Decimal(2).ln()),
        'mu': str(mu), 'S_max': S_max, 'log2_S_max': float(Decimal(S_max).ln() / Decimal(2).ln()),
        'M_max': M_max, 'unit_per_mult': '%d/%d (= %s)' % (20000 + c_over, 20000, unit),
        'additive_items_log2': {k: float((Decimal(v.numerator if isinstance(v, Fraction) else v) / Decimal(v.denominator if isinstance(v, Fraction) else 1)).ln() / Decimal(2).ln()) for k, v in add.items()},
        'log2_additive_total': float((Decimal(A.numerator if isinstance(A, Fraction) else A) / Decimal(A.denominator if isinstance(A, Fraction) else 1)).ln() / Decimal(2).ln()),
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
        # log2 W bits per charged instruction: M_max/(20000*W*ln 2) at c = 1773; 1/(20000 ln 2) bounds it
        'sensitivity_bits_per_instr_charged': float(Fraction(M_max, 20000) / W) / math.log(2),
        'sensitivity_bits_per_instr_at_zero_charge': 1.0 / (20000 * math.log(2)),
        'selector_two_half_minimum_tuples_checked': sel_res['exhaustive_tuples'],
        'selector_sequences_executed': sel_res,
        'ring_and_equality_sequences': ring_res,
        'log2_preprocessing': float((Decimal(pre.numerator if isinstance(pre, Fraction) else pre) / Decimal(pre.denominator if isinstance(pre, Fraction) else 1)).ln() / Decimal(2).ln()),
        'expected_reports_log2': float(expected_dp.ln() / Decimal(2).ln()), 'report_cap_log2': 47,
        'expected_reports_log2_if_48_bit_test': float(expected_reports_48),
        'charge_ceiling_leaving_claim_unchanged': ceiling, 'charge_room_above_charged': ceiling - c_over,
        'peak_memory_log2_bytes': math.log2(peak),
        'trap_lemma_log2_max_E_l_for_l_25_to_512': trap_max, 'trap_lemma_log2_sum': math.log2(trap_sum),
        'trap_cost_exponent_log2(L_cap*sum)': 87 + math.log2(trap_sum),
        'success_main_event': '3903/10000 under mu', 'success_claim': '39/100',
    }

# ---------------------------------------------------------------- algebra (production field)
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
            assert F['frob'](y) == F['pow'](y, p) == F['pow'](y, lam)      # Frobenius is x -> x^lambda on GT
            first = False
        kap0, k0 = sel(F, y)
        assert sel(F, y) == sel_ref(F, y)                                    # priced two-half minimum = direct scan
        Vfast = sel12(F, y)
        out['selector_fast_form_equals_general_form_on_GT'] &= (Vfast == sel12_general(F, y))
        assert Vfast == sel12_general(F, y)           # 3-multiplication form = 4-multiplication form
        for k in range(12):                                                  # candidates are the true coordinates
            assert frobk(F, y, k)[1][0] == Vfast[k]
        for m in range(12):
            ym = frobk(F, y, m)
            kap, km = sel(F, ym)
            assert (kap, km) == sel_ref(F, ym)                               # and on every conjugate's key
            assert kap == kap0 and km == (k0 - m) % 12
        t = F['pow'](y, rng.randrange(2, R_ORD))
        for m in range(12):                                                  # twisted step is class-equivariant
            ym = frobk(F, y, m); _, km = sel(F, ym)
            y1 = F['mul'](ym, frobk(F, t, (-km) % 12))
            y0 = F['mul'](y, frobk(F, t, (-k0) % 12))
            assert any(frobk(F, y0, j) == y1 for j in range(12))
        out['checked_elements'] += 1
    out['all_assertions_passed'] = True
    return out

# ---------------------------------------------------------------- number-theory helpers (toy)
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
        if q % 12 == 7 and is_prime(q) and pow(2, (q - 1) // 2, q) == q - 1:       # N(xi)=2 non-residue
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
    """Shared parameters of the emulated walk: the index and the distinguished-point bits come from independent salts, the
    idealised production rule; a distinguished state is reported and never stepped from."""
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

# ---------------------------------------------------------------- toy: real tower, end to end
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
    """One toy instance: the complete algorithm with the production window, on a tower of the production shape."""
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
        TW.append([fr[(12 - k) % 12] for k in range(12)])               # TW[h][k] = sigma^-k(T_h)
    S_max = isqrt_ceil(2 * MU * N * LN_MAIN_UP)          # the same exact budget as the ledger
    stats = dict(probes=0, escapes=0, dead=0)

    def trail(tid, track):
        r2 = random.Random(wseed * 7919 + tid)
        al = r2.randrange(rp); be = r2.randrange(rp)
        y = mul(F['pow'](g, al), F['pow'](z, be))
        kap, k = sel(F, y)
        if tid == 1:
            # the production selector returns the true first coordinates of c1(sigma^k(y)) in this tower
            Vf = sel12(F, y)
            assert all(frobk(F, y, j)[1][0] == Vf[j] for j in range(12))
            assert Vf == sel12_general(F, y)
        win = [kap]
        ring = HashRing(W.w, W.word1)               # the window hash is the word the index already uses (w_1)
        ring.seed(kap)                              # the trail start seeds all 24 slots (proof.md section 5)
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
                hit = ring.member(kap2)             # the priced hash scan and exact six-limb test
                assert hit == (kap2 in win), 'the hash ring and the key list disagree'
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
            assert sorted(set(ring.keys())) == sorted(win), 'the ring does not hold the last min(n+1, w) keys'
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
                # replay reproduces the untracked path (same endpoint and key)
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

# ---------------------------------------------------------------- sim: exponent-domain walk
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
        follow = None                                   # (earlier trail id) while the merged trail is followed to its end
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
    # control: plain birthday on N classes with the same estimator
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

# ---------------------------------------------------------------- main
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
