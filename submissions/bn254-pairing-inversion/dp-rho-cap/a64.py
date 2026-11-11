"""Scalar AArch64 subset: program generator and emulator for the walk step.

Every instruction emitted here is one ordinary A64 instruction (64-bit general
registers only).  The emulator executes the exact straight-line sequence, checks
immediate encodability, and the instruction count of a sequence is its length.
Used by experiment.py; standard library only.  Tower Fp12 = Fp6[w]/(w^2-v),
Fp6 = Fp2[v]/(v^3-(9+u)), Fp2 = Fp[u]/(u^2+1), with w equal to the pinned
basis variable (so v = w^2 and u = w^6 - 9).
"""
MASK = (1 << 64) - 1
ZR = 31                      # xzr
P_REGS = (25, 26, 27, 28)    # modulus p limbs, loaded once per run
N0 = 29                      # -p^-1 mod 2^64
CUR, TAB, NXT, S0, S1, CTL = 19, 20, 21, 22, 23, 24   # base registers

PAIR_OPS = ('ldp', 'stp')


class Prog:
    def __init__(self):
        self.ins = []

    def e(self, *t):
        op = t[0]
        if op in PAIR_OPS:
            assert t[1] != ZR or op == 'stp'
            if op == 'ldp':
                assert t[1] != t[2], 'ldp with equal destinations'
            assert t[4] % 8 == 0 and -512 <= t[4] <= 504, ('pair imm', t)
        elif op in ('ldr', 'str'):
            assert t[3] % 8 == 0 and 0 <= t[3] <= 32760, ('ldr imm', t)
        elif op in ('adds_imm', 'add_imm'):
            assert 0 <= t[3] <= 4095
        elif op == 'and_imm':
            v = t[3]
            assert isinstance(v, int) and 0 < v < MASK and (v & (v + 1)) == 0, 'contiguous low mask of 1..63 bits'
        elif op in ('add_lsl', 'sub_lsl', 'orr_lsl'):
            assert 0 <= t[4] <= 63
        self.ins.append(t)


def run(ins, X, mem):
    """Execute; returns list of (index, taken) for branches seen."""
    C = 0
    br = []

    def R(i):
        return 0 if i == ZR else X[i]

    def W(i, v):
        if i != ZR:
            X[i] = v & MASK
    for k, t in enumerate(ins):
        op = t[0]
        if op == 'ldp':
            a = R(t[3]) + t[4]
            W(t[1], mem[a]); W(t[2], mem[a + 8])
        elif op == 'stp':
            a = R(t[3]) + t[4]
            mem[a] = R(t[1]); mem[a + 8] = R(t[2])
        elif op == 'ldr':
            W(t[1], mem[R(t[2]) + t[3]])
        elif op == 'str':
            mem[R(t[2]) + t[3]] = R(t[1])
        elif op == 'mul':
            W(t[1], R(t[2]) * R(t[3]))
        elif op == 'umulh':
            W(t[1], (R(t[2]) * R(t[3])) >> 64)
        elif op == 'adds':
            s = R(t[2]) + R(t[3]); C = s >> 64; W(t[1], s)
        elif op == 'adds_imm':
            s = R(t[2]) + t[3]; C = s >> 64; W(t[1], s)
        elif op == 'add_imm':
            W(t[1], R(t[2]) + t[3])
        elif op == 'adcs':
            s = R(t[2]) + R(t[3]) + C; C = s >> 64; W(t[1], s)
        elif op == 'adc':
            W(t[1], R(t[2]) + R(t[3]) + C)
        elif op == 'subs':
            s = R(t[2]) - R(t[3]); C = int(s >= 0); W(t[1], s)
        elif op == 'sbcs':
            s = R(t[2]) - R(t[3]) - (1 - C); C = int(s >= 0); W(t[1], s)
        elif op == 'csel':
            take = (C == 1) if t[4] == 'cs' else (C == 0)
            W(t[1], R(t[2]) if take else R(t[3]))
        elif op == 'csetm':
            W(t[1], MASK if ((C == 0) if t[2] == 'cc' else (C == 1)) else 0)
        elif op == 'and':
            W(t[1], R(t[2]) & R(t[3]))
        elif op == 'and_imm':
            W(t[1], R(t[2]) & t[3])
        elif op == 'orr_lsl':
            W(t[1], R(t[2]) | ((R(t[3]) << t[4]) & MASK))
        elif op == 'add_lsl':
            W(t[1], R(t[2]) + ((R(t[3]) << t[4]) & MASK))
        elif op == 'sub_lsl':
            W(t[1], R(t[2]) - ((R(t[3]) << t[4]) & MASK))
        elif op == 'cbz':
            br.append((k, R(t[1]) == 0))
        elif op == 'bcs':
            br.append((k, C == 1))
        elif op == 'b':
            br.append((k, True))
        else:
            raise ValueError(op)
    return br


# ---------------------------------------------------------------- Fp macros
# An Fp location is (base_register, byte_offset); 4 little-endian words.

def fp_mul(P, d, a, b, mod=P_REGS, n0=N0):
    """CIOS Montgomery product d = a*b/2^256 mod m, operands canonical (< m < 2^254)."""
    P.e('ldp', 0, 1, a[0], a[1]); P.e('ldp', 2, 3, a[0], a[1] + 16)
    P.e('ldp', 4, 5, b[0], b[1]); P.e('ldp', 6, 7, b[0], b[1] + 16)
    A, B, h, m = (0, 1, 2, 3), (4, 5, 6, 7), (13, 14, 15, 16), 17
    t = [8, 9, 10, 11, 12]
    for j in range(4):
        P.e('mul', t[j], A[j], B[0])
    for j in range(4):
        P.e('umulh', h[j], A[j], B[0])
    P.e('adds', t[1], t[1], h[0]); P.e('adcs', t[2], t[2], h[1])
    P.e('adcs', t[3], t[3], h[2]); P.e('adc', t[4], h[3], ZR)
    for i in range(4):
        if i:
            for j in range(4):
                P.e('mul', h[j], A[j], B[i])
            P.e('adds', t[0], t[0], h[0]); P.e('adcs', t[1], t[1], h[1])
            P.e('adcs', t[2], t[2], h[2]); P.e('adcs', t[3], t[3], h[3])
            P.e('adc', t[4], ZR, ZR)
            for j in range(4):
                P.e('umulh', h[j], A[j], B[i])
            P.e('adds', t[1], t[1], h[0]); P.e('adcs', t[2], t[2], h[1])
            P.e('adcs', t[3], t[3], h[2]); P.e('adc', t[4], t[4], h[3])
        P.e('mul', m, t[0], n0)
        for j in range(4):
            P.e('mul', h[j], m, mod[j])
        P.e('adds', ZR, t[0], h[0]); P.e('adcs', t[1], t[1], h[1])
        P.e('adcs', t[2], t[2], h[2]); P.e('adcs', t[3], t[3], h[3])
        P.e('adc', t[4], t[4], ZR)
        for j in range(4):
            P.e('umulh', h[j], m, mod[j])
        P.e('adds', t[1], t[1], h[0]); P.e('adcs', t[2], t[2], h[1])
        P.e('adcs', t[3], t[3], h[2]); P.e('adc', t[4], t[4], h[3])
        t = t[1:] + t[:1]
    P.e('subs', h[0], t[0], mod[0])
    for j in (1, 2, 3):
        P.e('sbcs', h[j], t[j], mod[j])
    for j in range(4):
        P.e('csel', h[j], h[j], t[j], 'cs')
    P.e('stp', h[0], h[1], d[0], d[1]); P.e('stp', h[2], h[3], d[0], d[1] + 16)


def fp_add(P, d, a, b, mod=P_REGS):
    """d = a + b mod m for canonical inputs (a + b < 2m < 2^255)."""
    if a == b:
        P.e('ldp', 0, 1, a[0], a[1]); P.e('ldp', 2, 3, a[0], a[1] + 16)
        P.e('adds', 0, 0, 0); P.e('adcs', 1, 1, 1); P.e('adcs', 2, 2, 2); P.e('adc', 3, 3, 3)
    else:
        P.e('ldp', 0, 1, a[0], a[1]); P.e('ldp', 2, 3, a[0], a[1] + 16)
        P.e('ldp', 4, 5, b[0], b[1]); P.e('ldp', 6, 7, b[0], b[1] + 16)
        P.e('adds', 0, 0, 4); P.e('adcs', 1, 1, 5); P.e('adcs', 2, 2, 6); P.e('adc', 3, 3, 7)
    P.e('subs', 4, 0, mod[0]); P.e('sbcs', 5, 1, mod[1])
    P.e('sbcs', 6, 2, mod[2]); P.e('sbcs', 7, 3, mod[3])
    for j in range(4):
        P.e('csel', 4 + j, 4 + j, j, 'cs')
    P.e('stp', 4, 5, d[0], d[1]); P.e('stp', 6, 7, d[0], d[1] + 16)


def fp_sub(P, d, a, b, mod=P_REGS):
    """d = a - b mod m for canonical inputs."""
    P.e('ldp', 0, 1, a[0], a[1]); P.e('ldp', 2, 3, a[0], a[1] + 16)
    P.e('ldp', 4, 5, b[0], b[1]); P.e('ldp', 6, 7, b[0], b[1] + 16)
    P.e('subs', 0, 0, 4); P.e('sbcs', 1, 1, 5); P.e('sbcs', 2, 2, 6); P.e('sbcs', 3, 3, 7)
    P.e('csetm', 8, 'cc')
    for j in range(4):
        P.e('and', 4 + j, mod[j], 8)
    P.e('adds', 0, 0, 4); P.e('adcs', 1, 1, 5); P.e('adcs', 2, 2, 6); P.e('adc', 3, 3, 7)
    P.e('stp', 0, 1, d[0], d[1]); P.e('stp', 2, 3, d[0], d[1] + 16)


# ---------------------------------------------------------------- scratch
class Scratch:
    """Hands out 32-byte Fp slots in the two scratch windows (offsets -512..480)."""

    def __init__(self):
        self.free = [(reg, off) for reg in (S0, S1) for off in range(-512, 481, 32)]
        self.peak = 0
        self.total = len(self.free)

    def fp(self):
        loc = self.free.pop(0)
        self.peak = max(self.peak, self.total - len(self.free))
        return loc

    def fp2(self):
        return (self.fp(), self.fp())

    def fp6(self):
        return (self.fp2(), self.fp2(), self.fp2())

    def put(self, x):
        if isinstance(x[0], int):
            self.free.insert(0, x)
        else:
            for y in x:
                self.put(y)


# ---------------------------------------------------------------- Fp2 = Fp[u]/(u^2+1)
def f2_add(P, d, a, b):
    fp_add(P, d[0], a[0], b[0]); fp_add(P, d[1], a[1], b[1])


def f2_sub(P, d, a, b):
    fp_sub(P, d[0], a[0], b[0]); fp_sub(P, d[1], a[1], b[1])


def f2_mul(P, S, d, a, b):
    """Karatsuba: 3 products, 2 additions, 3 subtractions.  d may alias a or b."""
    v0, v1, s, t = S.fp(), S.fp(), S.fp(), S.fp()
    fp_mul(P, v0, a[0], b[0]); fp_mul(P, v1, a[1], b[1])
    fp_add(P, s, a[0], a[1]); fp_add(P, t, b[0], b[1])
    fp_mul(P, s, s, t)
    fp_sub(P, d[0], v0, v1)
    fp_sub(P, s, s, v0); fp_sub(P, d[1], s, v1)
    S.put([v0, v1, s, t])


def f2_mulxi(P, S, d, a):
    """d = a*(9+u) = (9a0 - a1) + (a0 + 9a1)u; d must not alias a."""
    t = S.fp()
    fp_add(P, t, a[0], a[0]); fp_add(P, t, t, t); fp_add(P, t, t, t); fp_add(P, t, t, a[0])
    fp_sub(P, d[0], t, a[1])
    fp_add(P, t, a[1], a[1]); fp_add(P, t, t, t); fp_add(P, t, t, t); fp_add(P, t, t, a[1])
    fp_add(P, d[1], a[0], t)
    S.put(t)


# ---------------------------------------------------------------- Fp6 = Fp2[v]/(v^3 - xi)
def f6_add(P, d, a, b):
    for k in range(3):
        f2_add(P, d[k], a[k], b[k])


def f6_sub(P, d, a, b):
    for k in range(3):
        f2_sub(P, d[k], a[k], b[k])


def f6_mul(P, S, d, a, b):
    """Karatsuba-style 6 Fp2 products.  All input reads precede output writes."""
    v0, v1, v2 = S.fp2(), S.fp2(), S.fp2()
    f2_mul(P, S, v0, a[0], b[0]); f2_mul(P, S, v1, a[1], b[1]); f2_mul(P, S, v2, a[2], b[2])
    q = []
    for i, j in ((1, 2), (0, 1), (0, 2)):
        s, t = S.fp2(), S.fp2()
        f2_add(P, s, a[i], a[j]); f2_add(P, t, b[i], b[j]); f2_mul(P, S, s, s, t)
        S.put(t)
        q.append(s)
    s12, s01, s02 = q
    x = S.fp2()
    f2_sub(P, s12, s12, v1); f2_sub(P, s12, s12, v2); f2_mulxi(P, S, x, s12)
    f2_add(P, d[0], x, v0)
    f2_sub(P, s01, s01, v0); f2_sub(P, s01, s01, v1); f2_mulxi(P, S, x, v2)
    f2_add(P, d[1], s01, x)
    f2_sub(P, s02, s02, v0); f2_add(P, s02, s02, v1); f2_sub(P, d[2], s02, v2)
    S.put([v0, v1, v2, s12, s01, s02, x])


# ---------------------------------------------------------------- Fp12 = Fp6[w]/(w^2 - v)
def f12_mul(P, S, d, a, b):
    t0, t1, sa, sb = S.fp6(), S.fp6(), S.fp6(), S.fp6()
    f6_mul(P, S, t0, a[0], b[0]); f6_mul(P, S, t1, a[1], b[1])
    f6_add(P, sa, a[0], a[1]); f6_add(P, sb, b[0], b[1]); f6_mul(P, S, sa, sa, sb)
    f6_sub(P, sa, sa, t0); f6_sub(P, d[1], sa, t1)
    x = S.fp2()
    f2_mulxi(P, S, x, t1[2])
    f2_add(P, d[0][0], t0[0], x); f2_add(P, d[0][1], t0[1], t1[0]); f2_add(P, d[0][2], t0[2], t1[1])
    S.put([t0, t1, sa, sb, x])


def f12_loc(reg):
    """Tower Fp12 at offsets 0..383 of reg: part k, Fp2 index m, component c."""
    return tuple(tuple(((reg, 192 * k + 64 * m), (reg, 192 * k + 64 * m + 32)) for m in range(3))
                 for k in range(2))


# ---------------------------------------------------------------- walk step
# CTL window (non-negative offsets): labels, r, counter, cap, table pointers.
CTL_A, CTL_B, CTL_R, CTL_CNT, CTL_CAP, CTL_TAB, CTL_NDP, CTL_DPB, CTL_DPMAX = 0, 32, 64, 96, 112, 128, 136, 144, 152
REC_BYTES = 448              # 384-byte element, c_j at 384, d_j at 416
REC_C, REC_D = 384, 416
PART_OFF = 32                # word 0 of A0.im  (partition: low 10 bits)
DP_SHIFT = 34                # DP: word0 of A0.re == 0 and low 30 bits of word1 == 0 (94 bits)


def dp_test(P, nxt):
    """Branch taken iff word0 of A0.re is 0 and the low 30 bits of word1 are 0."""
    P.e('ldp', 0, 1, nxt, 0)
    P.e('orr_lsl', 0, 0, 1, DP_SHIFT)
    P.e('cbz', 0)


def step_program(cur, nxt):
    P, S = Prog(), Scratch()
    P.e('ldr', 0, cur, PART_OFF)
    P.e('and_imm', 0, 0, 0x3ff)
    P.e('ldr', 1, CTL, CTL_TAB)
    P.e('add_lsl', TAB, 1, 0, 9)
    P.e('sub_lsl', TAB, TAB, 0, 6)
    f12_mul(P, S, f12_loc(nxt), f12_loc(cur), f12_loc(TAB))
    rr = (9, 10, 11, 12)
    P.e('ldp', 9, 10, CTL, CTL_R); P.e('ldp', 11, 12, CTL, CTL_R + 16)
    fp_add(P, (CTL, CTL_A), (CTL, CTL_A), (TAB, REC_C), mod=rr)
    fp_add(P, (CTL, CTL_B), (CTL, CTL_B), (TAB, REC_D), mod=rr)
    dp_test(P, nxt)
    P.e('ldp', 2, 3, CTL, CTL_CNT)
    P.e('adds_imm', 2, 2, 1); P.e('adc', 3, 3, ZR)
    P.e('stp', 2, 3, CTL, CTL_CNT)
    P.e('ldp', 4, 5, CTL, CTL_CAP)
    P.e('subs', ZR, 2, 4); P.e('sbcs', ZR, 3, 5)
    P.e('bcs',)
    P.e('b',)              # loop back edge (charged to every step)
    return P, S


def dp_program(nxt):
    """Rare path: append (Y', a, b) as one 448-byte record; stop when the table is full."""
    P = Prog()
    P.e('ldr', 2, CTL, CTL_NDP); P.e('ldr', 3, CTL, CTL_DPB)
    P.e('add_lsl', 4, 3, 2, 9); P.e('sub_lsl', 4, 4, 2, 6)
    for k in range(0, 384, 16):
        P.e('ldp', 5, 6, nxt, k); P.e('stp', 5, 6, 4, k)
    for k in range(0, 64, 16):
        P.e('ldp', 5, 6, CTL, CTL_A + k); P.e('stp', 5, 6, 4, 384 + k)
    P.e('add_imm', 2, 2, 1); P.e('str', 2, CTL, CTL_NDP)
    P.e('ldr', 3, CTL, CTL_DPMAX); P.e('subs', ZR, 2, 3); P.e('bcs',)
    P.e('b',)
    return P


def counts():
    out = {}
    for name, fn in (('fp_mul', lambda P, S: fp_mul(P, (S0, 0), (S0, 32), (S0, 64))),
                     ('fp_add', lambda P, S: fp_add(P, (S0, 0), (S0, 32), (S0, 64))),
                     ('fp_dbl', lambda P, S: fp_add(P, (S0, 0), (S0, 32), (S0, 32))),
                     ('fp_sub', lambda P, S: fp_sub(P, (S0, 0), (S0, 32), (S0, 64)))):
        P = Prog(); fn(P, None); out[name] = len(P.ins)
    for name, fn in (('fp2_mul', lambda P, S: f2_mul(P, S, S.fp2(), S.fp2(), S.fp2())),
                     ('fp2_mulxi', lambda P, S: f2_mulxi(P, S, S.fp2(), S.fp2())),
                     ('fp6_mul', lambda P, S: f6_mul(P, S, S.fp6(), S.fp6(), S.fp6())),
                     ('fp12_mul', lambda P, S: f12_mul(P, S, f12_loc(NXT), f12_loc(CUR), f12_loc(TAB)))):
        P, S = Prog(), Scratch(); fn(P, S); out[name] = len(P.ins)
    sp, ss = step_program(CUR, NXT)
    out['step'] = len(sp.ins)
    out['step_scratch_slots_peak'] = ss.peak
    out['dp_path'] = len(dp_program(NXT).ins)
    hist = {}
    for t in sp.ins:
        hist[t[0]] = hist.get(t[0], 0) + 1
    out['step_opcode_histogram'] = dict(sorted(hist.items()))
    return out
