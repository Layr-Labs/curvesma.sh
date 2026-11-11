"""Lazy-reduction BN254 Fp12 routines and the gated expected-work ledger (standard library only).

Part L  the routines priced by this package -- dense Fp12 product with a precomputed table operand,
        Granger-Scott cyclotomic squaring, the twelve Frobenius keys with their argmin -- are generated as
        straight-line scalar AArch64 programs (64-bit general registers only), executed instruction by
        instruction on a simulator, and compared with big-integer arithmetic on the protected BN254 GT.
        The executed instruction counts are the prices. While generating, the generator asserts a static
        worst-case bound for every stored value (operands below 2^256, products below 2^512, signed 9-word
        intermediates inside (-2^575, 2^575), REDC inputs nonnegative, canonicalisation bound below 2p).
Part G  the gate and expected-work ledger with these prices, in exact decimal arithmetic.
Base    the unchanged experiment.py of submission 2ece9f04 is then executed in-process with its own input
        (identities on the protected GT, production step rule, toy-field runs of the complete algorithm).
Inputs (stdin JSON): seed, lazy_cases, lazy_gt_cases, mean_allowance_percent, cap_multiple, listings,
base_run, base_seed, base_seconds, base_production_steps.
"""
import json
import random
import sys
import time
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction
from pathlib import Path

T_START = time.time()
cfg = json.load(sys.stdin)
SEED = cfg.get("seed", 1)
SHOW_LISTINGS = cfg.get("listings", False)
rng = random.Random(SEED)
base = Path("/snapshot/files")
if not (base / "target-profile.json").is_file():
    base = Path(__file__).resolve().parent
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
G_PROFILE = [int(c, 16) for c in profile["generator_g"]["coefficients"]]
RMONT = 1 << 256


class Tower:
    """Reference arithmetic (big integers): Fp12 = Fp2[w]/(w^6 - xi), xi = 9 + u, Fp2 = Fp[u]/(u^2+1)."""

    def __init__(self, q):
        self.q = q
        self.one = ((1, 0),) + ((0, 0),) * 5

    def f2m(self, a, b):
        q = self.q
        return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

    def f2pow(self, a, e):
        res = (1, 0)
        while e:
            if e & 1:
                res = self.f2m(res, a)
            a = self.f2m(a, a)
            e >>= 1
        return res

    def mul(self, x, y):
        q = self.q
        acc = [[0, 0] for _ in range(11)]
        for i in range(6):
            for j in range(6):
                t = acc[i + j]
                t[0] += x[i][0] * y[j][0] - x[i][1] * y[j][1]
                t[1] += x[i][0] * y[j][1] + x[i][1] * y[j][0]
        out = []
        for m in range(6):
            c0, c1 = acc[m]
            if m + 6 <= 10:
                h0, h1 = acc[m + 6]
                c0 += 9 * h0 - h1
                c1 += 9 * h1 + h0
            out.append((c0 % q, c1 % q))
        return tuple(out)

    def pow(self, x, e):
        res = self.one
        while e:
            if e & 1:
                res = self.mul(res, x)
            x = self.mul(x, x)
            e >>= 1
        return res

    def frob(self, y, k):
        """y^(q^k) by repeated Frobenius: (sum c_m w^m)^q = sum conj(c_m) xi^(m(q-1)/6) w^m"""
        q = self.q
        d = [self.f2pow((9, 1), m * (q - 1) // 6) for m in range(6)]
        for _ in range(k):
            y = tuple(self.f2m((c[0], (-c[1]) % q), d[m]) for m, c in enumerate(y))
        return y

    def keys(self, y):
        """Re(w^1 coefficient of y^(q^k)), k = 0..11, computed from the twelve explicit conjugates"""
        out, z = [], y
        for k in range(12):
            out.append(z[1][0])
            z = self.frob(z, 1)
        return out


def profile_to_tower(c, q):
    return tuple(((c[m] + 9 * c[m + 6]) % q, c[m + 6] % q) for m in range(6))


# =============================================================== Part L: generator and simulator
# ---------------------------------------------------------------------------------------
# Lazy-reduction AArch64 program generator and instruction-level simulator for BN254 Fp12.
# Every routine is emitted as one straight-line program of scalar AArch64 instructions on
# 64-bit general registers; the simulator executes it and counts every instruction.
# Static worst-case bounds are tracked for every stored value and asserted before use.
# ---------------------------------------------------------------------------------------

M64 = (1 << 64) - 1
ZR = 31          # xzr
SP = 32          # stack pointer (only used as a base register)
RNAME = {ZR: "xzr", SP: "sp"}


def rn(r):
    return RNAME.get(r, "x%d" % r)


# register roles
OUT, Y, TA, TB, CST = 0, 1, 2, 3, 4          # output, input y, table windows, constants
SCR = [5, 6, 7, 8]                           # scratch windows
POOL = [9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
A = POOL[0:9]       # 9-word accumulator
V = POOL[9:18]      # term / operand temporaries
X1, X2, X3 = POOL[18:21]
WLO, WHI = -512, 504


class Val:
    """A stored value: kind 'S' (4 words, integer in [0, mx]) or 'D' (9-word two's complement,
    or 8 nonnegative words for a raw product), true integer in [mn, mx]."""
    __slots__ = ("base", "off", "kind", "width", "mn", "mx", "name")

    def __init__(self, base, off, kind, width, mn, mx, name=""):
        self.base, self.off, self.kind, self.width, self.mn, self.mx, self.name = base, off, kind, width, mn, mx, name


class Gen:
    def __init__(self, p, consts_layout):
        self.p = p
        self.code = []
        self.cst = consts_layout          # dict name -> offset from CST
        self.free = {b: list(range(WLO, 440 - 79 + 80, 80))[::-1] for b in SCR}   # 80-byte slots
        self.free = {b: [o for o in self.free[b] if o + 72 <= WHI] for b in SCR}
        self.peak = 0
        self.live = 0
        self.acc_mn = self.acc_mx = 0
        self.acc8_zero = False
        self.kconsts = {}

    # ---------------------------------------------------------------- emission helpers
    def e(self, *ins):
        self.code.append(ins)

    def alloc(self, kind, width, mn, mx, name=""):
        for b in SCR:
            if self.free[b]:
                off = self.free[b].pop()
                self.live += 1
                self.peak = max(self.peak, self.live)
                return Val(b, off, kind, width, mn, mx, name)
        raise RuntimeError("out of scratch slots")

    def release(self, v):
        if v.base in SCR:
            self.free[v.base].append(v.off)
            self.live -= 1

    def ld2(self, r0, r1, base, off):
        assert WLO <= off <= WHI and off % 8 == 0
        self.e("ldp", r0, r1, base, off)

    def st2(self, r0, r1, base, off):
        assert WLO <= off <= WHI and off % 8 == 0
        self.e("stp", r0, r1, base, off)

    def ldc(self, r0, r1, name):
        self.ld2(r0, r1, CST, self.cst.get(name, 0))

    # ---------------------------------------------------------------- single-width values
    def load_s(self, v, regs):
        self.ld2(regs[0], regs[1], v.base, v.off)
        self.ld2(regs[2], regs[3], v.base, v.off + 16)

    def store_s(self, regs, v):
        self.st2(regs[0], regs[1], v.base, v.off)
        self.st2(regs[2], regs[3], v.base, v.off + 16)

    def sadd(self, a, b, reduce):
        """a+b (unreduced: < 2^256 asserted) or (a+b) mod p for canonical inputs."""
        r = V[0:4]
        t = V[4:8]
        self.load_s(a, r)
        self.load_s(b, t)
        self.e("adds", r[0], r[0], t[0])
        for j in (1, 2, 3):
            self.e("adcs" if j < 3 else "adc", r[j], r[j], t[j])
        if reduce:
            assert a.mx < self.p and b.mx < self.p
            self.ldc(t[0], t[1], "p0")
            self.ldc(t[2], t[3], "p2")
            self.e("subs", t[0], r[0], t[0])
            for j in (1, 2, 3):
                self.e("sbcs", t[j], r[j], t[j])
            for j in range(4):
                self.e("csel", r[j], t[j], r[j], "cs")
            mx = self.p - 1
        else:
            mx = a.mx + b.mx
            assert mx < 1 << 256
        out = self.alloc("S", 4, 0, mx)
        self.store_s(r, out)
        return out

    def ssub_kp(self, a, b, k):
        """a - b + k*p (unreduced), requires k*p >= b.mx; result < 2^256 asserted."""
        assert k * self.p >= b.mx
        r = V[0:4]
        t = V[4:8]
        self.load_s(a, r)
        self.load_s(b, t)
        self.e("subs", r[0], r[0], t[0])
        for j in (1, 2, 3):
            self.e("sbcs", r[j], r[j], t[j])
        self.ldc(t[0], t[1], "kp%d_0" % k)
        self.ldc(t[2], t[3], "kp%d_2" % k)
        self.e("adds", r[0], r[0], t[0])
        for j in (1, 2, 3):
            self.e("adcs" if j < 3 else "adc", r[j], r[j], t[j])
        mx = a.mx + k * self.p
        assert mx < 1 << 256
        out = self.alloc("S", 4, 0, mx)
        self.store_s(r, out)
        return out

    # ---------------------------------------------------------------- 256x256 product
    def mul_to_acc(self, a, b):
        """acc[0..7] <- a*b (8 words, nonnegative); acc[8] implicitly zero."""
        assert a.mx < 1 << 256 and b.mx < 1 << 256
        B = V[0:4]
        ai = V[4]
        L = V[5:9]
        self.load_s(b, B)
        self.e("ldr", ai, a.base, a.off)
        for j in range(4):
            self.e("mul", A[j], ai, B[j])
        for j in range(4):
            self.e("umulh", L[j], ai, B[j])
        self.e("adds", A[1], A[1], L[0])
        self.e("adcs", A[2], A[2], L[1])
        self.e("adcs", A[3], A[3], L[2])
        self.e("adc", A[4], L[3], ZR)
        for i in (1, 2, 3):
            self.e("ldr", ai, a.base, a.off + 8 * i)
            for j in range(4):
                self.e("mul", L[j], ai, B[j])
            self.e("adds", A[i], A[i], L[0])
            for j in (1, 2, 3):
                self.e("adcs", A[i + j], A[i + j], L[j])
            self.e("adc", A[i + 4], ZR, ZR)
            for j in range(4):
                self.e("umulh", L[j], ai, B[j])
            self.e("adds", A[i + 1], A[i + 1], L[0])
            self.e("adcs", A[i + 2], A[i + 2], L[1])
            self.e("adcs", A[i + 3], A[i + 3], L[2])
            self.e("adc", A[i + 4], A[i + 4], L[3])
        self.acc_mn, self.acc_mx = 0, a.mx * b.mx
        assert self.acc_mx < 1 << 512
        self.acc8_zero = True

    def mul(self, a, b):
        self.mul_to_acc(a, b)
        return self.store_acc(width=8)

    # ---------------------------------------------------------------- accumulator
    def store_acc(self, width=9):
        if width == 8:
            assert self.acc_mn >= 0 and self.acc_mx < 1 << 512
        else:
            assert -(1 << 575) < self.acc_mn and self.acc_mx < 1 << 575
        v = self.alloc("D", width, self.acc_mn, self.acc_mx)
        for k in range(4):
            self.st2(A[2 * k], A[2 * k + 1], v.base, v.off + 16 * k)
        if width == 9:
            self.st2(ZR if self.acc8_zero else A[8], ZR, v.base, v.off + 64)
        return v

    def acc_set(self, v):
        for k in range(4):
            self.ld2(A[2 * k], A[2 * k + 1], v.base, v.off + 16 * k)
        if v.width == 9:
            self.ld2(A[8], X3, v.base, v.off + 64)
            self.acc8_zero = False
        else:
            self.acc8_zero = True
        self.acc_mn, self.acc_mx = v.mn, v.mx

    @staticmethod
    def naf(c):
        out, k = [], 0
        while c:
            if c & 1:
                d = 2 - (c & 3)
                out.append((d, k))
                c -= d
            c >>= 1
            k += 1
        return out          # list of (+-1, shift)

    def _chain(self, sign, words):
        """acc += / -= the 9 words (registers; ZR allowed)"""
        a8 = ZR if self.acc8_zero else A[8]
        if sign > 0:
            self.e("adds", A[0], A[0], words[0])
            for j in range(1, 8):
                self.e("adcs", A[j], A[j], words[j])
            self.e("adc", A[8], a8, words[8])
        else:
            self.e("subs", A[0], A[0], words[0])
            for j in range(1, 8):
                self.e("sbcs", A[j], A[j], words[j])
            self.e("sbc", A[8], a8, words[8])
        self.acc8_zero = False

    def acc_add(self, c, v):
        """acc += c * v for a small nonzero integer c (signed binary digits; one pass per digit)."""
        for k in range(4):
            self.ld2(V[2 * k], V[2 * k + 1], v.base, v.off + 16 * k)
        if v.width == 9:
            self.ld2(V[8], X3, v.base, v.off + 64)
            w = V[:9]
        else:
            w = V[:8] + [ZR]
        for d, k in self.naf(abs(c)):
            sign = d * (1 if c > 0 else -1)
            if k == 0:
                self._chain(sign, w)
            else:
                # shifted words are formed one at a time into X1/X2 (EXTR/LSL do not touch flags)
                sh = []
                ins = [("lsl", X1, w[0], k)] + [("extr", X1 if j % 2 == 0 else X2, w[j], w[j - 1], 64 - k) for j in range(1, 9)]
                a8 = ZR if self.acc8_zero else A[8]
                for j in range(9):
                    self.e(*ins[j])
                    src = ins[j][1]
                    if j == 0:
                        self.e("adds" if sign > 0 else "subs", A[0], A[0], src)
                    elif j < 8:
                        self.e("adcs" if sign > 0 else "sbcs", A[j], A[j], src)
                    else:
                        self.e("adc" if sign > 0 else "sbc", A[8], a8, src)
                self.acc8_zero = False
        if c > 0:
            self.acc_mn += c * v.mn
            self.acc_mx += c * v.mx
        else:
            self.acc_mn += c * v.mx
            self.acc_mx += c * v.mn
        assert -(1 << 575) < self.acc_mn and self.acc_mx < 1 << 575

    def acc_add_high(self, c, v):
        """acc += c * v * 2^256 for a single-width value v (Montgomery-domain additive term)."""
        for k in range(2):
            self.ld2(V[2 * k], V[2 * k + 1], v.base, v.off + 16 * k)
        w = V[:4] + [ZR]
        for d, k in self.naf(abs(c)):
            sign = d * (1 if c > 0 else -1)
            a8 = ZR if self.acc8_zero else A[8]
            if k == 0:
                srcs = w
                pre = [None] * 5
            else:
                pre = [("lsl", X1, w[0], k)] + [("extr", X1 if j % 2 == 0 else X2, w[j], w[j - 1], 64 - k) for j in range(1, 4)] + [("lsr", X1, w[3], 64 - k)]
                srcs = [x[1] for x in pre]
            for j in range(5):
                if pre[j]:
                    self.e(*pre[j])
                if j == 0:
                    self.e("adds" if sign > 0 else "subs", A[4], A[4], srcs[0])
                elif j < 4:
                    self.e("adcs" if sign > 0 else "sbcs", A[4 + j], A[4 + j], srcs[j])
                else:
                    self.e("adc" if sign > 0 else "sbc", A[8], a8, srcs[4])
            self.acc8_zero = False
        sh = 1 << 256
        if c > 0:
            self.acc_mn += c * v.mn * sh
            self.acc_mx += c * v.mx * sh
        else:
            self.acc_mn += c * v.mx * sh
            self.acc_mx += c * v.mn * sh

    def acc_add_k(self):
        """add the routine's single constant K (a multiple of p) making every accumulator nonnegative.
        First build pass: records the largest multiple needed; second pass: uses that one K everywhere."""
        if self.acc_mn >= 0:
            return 0
        k = (-self.acc_mn + self.p - 1) // self.p
        if self.kconsts:
            K = max(self.kconsts)
            assert K >= -self.acc_mn
        else:
            self.kneed = max(getattr(self, "kneed", 0), k * self.p)
            K = k * self.p
        for j in range(4):
            self.ldc(V[2 * j], V[2 * j + 1], "K0_%d" % (2 * j))
        self.ldc(V[8], X3, "K0_8")
        self._chain(1, V[:9])
        self.acc_mn += K
        self.acc_mx += K
        return K

    # ---------------------------------------------------------------- reduction
    def redc_canon(self, dest):
        """dest <- REDC(acc) reduced to [0, p): word-level Montgomery reduction, quotient by UDIV of the
        top 64 bits of t >> 200, one subtraction of q*p and a final conditional subtraction."""
        p = self.p
        assert self.acc_mn >= 0
        Tmx = self.acc_mx
        assert Tmx < 1 << 575
        PP = V[0:4]
        NN, DD = V[4], V[5]
        m = V[6]
        Lt = [V[7], V[8], X1]
        H = [X2, X3, Lt[0], Lt[1]]
        self.ldc(PP[0], PP[1], "p0")
        self.ldc(PP[2], PP[3], "p2")
        self.ldc(NN, DD, "n0")
        for i in range(4):
            self.e("subs", ZR, A[i], ("#", 1))
            self.e("mul", m, A[i], NN)
            for j in (1, 2, 3):
                self.e("mul", Lt[j - 1], m, PP[j])
            for j in (1, 2, 3):
                self.e("adcs", A[i + j], A[i + j], Lt[j - 1])
            for j in range(4):
                self.e("umulh", H[j], m, PP[j])
            self.e("adc", H[3], H[3], ZR)
            self.e("adds", A[i + 1], A[i + 1], H[0])
            self.e("adcs", A[i + 2], A[i + 2], H[1])
            self.e("adcs", A[i + 3], A[i + 3], H[2])
            if i + 4 < 8:
                self.e("adcs", A[i + 4], A[i + 4], H[3])
                for k in range(i + 5, 8):
                    self.e("adcs", A[k], A[k], ZR)
                self.e("adc", A[8], ZR if self.acc8_zero else A[8], ZR)
            else:
                self.e("adc", A[8], ZR if self.acc8_zero else A[8], H[3])
            self.acc8_zero = False
        # t = A4..A8, t < Tmx/2^256 + p
        tmx = (Tmx + ((1 << 256) - 1) * p) >> 256
        assert tmx < 1 << 264
        dd = (p >> 200) + 1
        # q = floor(floor(t/2^200)/dd) <= t/p, so t - q p >= 0 ; bound t - q p:
        hi_mn = Fraction(tmx, 1 << 200) - 1
        bound = Fraction(tmx) - p * (hi_mn / dd - 1)
        bound = max(bound, Fraction(p))  # small t: q may be 0
        assert bound < 2 * p, "canonicalisation bound fails"
        q = m
        self.e("extr", Lt[0], A[8], A[7], 8)
        self.e("udiv", q, Lt[0], DD)
        G = [X1, X2, X3, NN]
        Lq = [Lt[0], Lt[1], A[8], DD]
        for j in range(4):
            self.e("mul", Lq[j], q, PP[j])
        for j in range(3):
            self.e("umulh", G[j], q, PP[j])
        self.e("adds", Lq[1], Lq[1], G[0])
        self.e("adcs", Lq[2], Lq[2], G[1])
        self.e("adc", Lq[3], Lq[3], G[2])
        self.e("subs", A[4], A[4], Lq[0])
        for j in (1, 2, 3):
            self.e("sbcs", A[4 + j], A[4 + j], Lq[j])
        # A4..A7 now hold t - q p in [0, 2p); conditional subtraction of p
        U = [X1, X2, X3, NN]
        self.e("subs", U[0], A[4], PP[0])
        for j in (1, 2, 3):
            self.e("sbcs", U[j], A[4 + j], PP[j])
        for j in range(4):
            self.e("csel", A[4 + j], U[j], A[4 + j], "cs")
        self.st2(A[4], A[5], dest.base, dest.off)
        self.st2(A[6], A[7], dest.base, dest.off + 16)
        self.acc8_zero = False
        return dest

    # ---------------------------------------------------------------- linear combination
    def lin(self, terms, highs=()):
        """acc <- sum c*v (terms: list of (c, Val)), plus sum c*v*2^256 for single-width highs."""
        i0 = [k for k, t in enumerate(terms) if t[0] == 1]
        if not i0:
            raise RuntimeError("a term with coefficient +1 is required")
        self.acc_set(terms[i0[0]][1])
        rest = terms[:i0[0]] + terms[i0[0] + 1:]
        for c, v in rest:
            self.acc_add(c, v)
        for c, v in highs:
            self.acc_add_high(c, v)


def run(code, mem, regs):
    """Execute a straight-line program; returns the number of executed instructions."""
    C = 0
    n = 0
    for ins in code:
        n += 1
        op = ins[0]
        if op == "mul":
            regs[ins[1]] = (regs[ins[2]] * regs[ins[3]]) & M64
        elif op == "umulh":
            regs[ins[1]] = (regs[ins[2]] * regs[ins[3]]) >> 64
        elif op == "adds" or op == "adcs" or op == "adc" or op == "add":
            s = regs[ins[2]] + regs[ins[3]] + (C if op != "adds" and op != "add" else 0)
            regs[ins[1]] = s & M64
            if op != "adc" and op != "add":
                C = s >> 64
        elif op == "subs" or op == "sbcs" or op == "sbc":
            b = ins[3][1] if isinstance(ins[3], tuple) else regs[ins[3]]
            d = regs[ins[2]] - b - (0 if op == "subs" else 1 - C)
            regs[ins[1]] = d & M64
            if op != "sbc":
                C = 1 if d >= 0 else 0
        elif op == "ldp":
            a = regs[ins[3]] + ins[4]
            regs[ins[1]] = mem[a]
            regs[ins[2]] = mem.get(a + 8, 0)
        elif op == "stp":
            a = regs[ins[3]] + ins[4]
            mem[a] = regs[ins[1]]
            mem[a + 8] = regs[ins[2]]
        elif op == "ldr":
            regs[ins[1]] = mem[regs[ins[2]] + ins[3]]
        elif op == "str":
            mem[regs[ins[2]] + ins[3]] = regs[ins[1]]
        elif op == "csel":
            regs[ins[1]] = regs[ins[2]] if (C == 1) == (ins[4] == "cs") else regs[ins[3]]
        elif op == "csetm":
            regs[ins[1]] = M64 if (C == 1) == (ins[2] == "cs") else 0
        elif op == "and":
            regs[ins[1]] = regs[ins[2]] & regs[ins[3]]
        elif op == "extr":
            regs[ins[1]] = (((regs[ins[2]] << 64) | regs[ins[3]]) >> ins[4]) & M64
        elif op == "lsl":
            regs[ins[1]] = (regs[ins[2]] << ins[3]) & M64
        elif op == "lsr":
            regs[ins[1]] = regs[ins[2]] >> ins[3]
        elif op == "udiv":
            regs[ins[1]] = regs[ins[2]] // regs[ins[3]] if regs[ins[3]] else 0
        elif op == "addi":
            regs[ins[1]] = (regs[ins[2]] + ins[3]) & M64
        elif op == "movi":
            regs[ins[1]] = ins[2]
        elif op == "ret" or op == "bl":
            pass
        else:
            raise RuntimeError("unsupported " + op)
        regs[ZR] = 0
    return n


P_BN = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47
RM = 1 << 256

def consts_layout(p, kconsts):
    """constant pool layout: p, n0/dd, k*p for k=1..4, K constants"""
    lay, vals = {}, {}
    off = WLO
    def put(name, words):
        nonlocal off
        for i in range(0, len(words), 2):
            lay[name + "_%d" % i] = off
            vals[off] = words[i]; vals[off + 8] = words[i + 1] if i + 1 < len(words) else 0
            off += 16
    def limbs(v, n): return [(v >> (64 * i)) & M64 for i in range(n)]
    put("p", limbs(p, 4))
    lay["p0"], lay["p2"] = lay["p_0"], lay["p_2"]
    n0 = (-pow(p, -1, 1 << 64)) % (1 << 64)
    lay["n0"] = off; vals[off] = n0; vals[off + 8] = (p >> 200) + 1; off += 16
    for k in (1, 2, 4):
        put("kp%d" % k, limbs(k * p, 4))
    for K, name in kconsts.items():
        put(name, limbs(K, 10))
    return lay, vals

class Lazy:
    """Builds the three routines; returns Gen objects (code, counts) and the input/output layout."""
    def __init__(self, p):
        self.p = p

    def new_gen(self, kconsts=None):
        # two-pass: first pass collects K constants, second pass lays them out
        lay, vals = consts_layout(self.p, kconsts or {})
        g = Gen(self.p, lay)
        if kconsts: g.kconsts = dict(kconsts)
        g.cvals = vals
        return g

    def prologue(self, g, table=False):
        # frame: save x19-x30 (6 STP), allocate 4 scratch windows of 1 KiB below the saves
        g.e("addi", SP, SP, -4096 - 112)
        for k in range(6):
            g.e("stp", 19 + 2 * k, 20 + 2 * k, SP, 4096 + 16 * k)
        for i, b in enumerate(SCR):
            g.e("addi", b, SP, 512 + 1024 * i)
        if table:
            g.e("addi", TB, TA, 1024)
        g.e("addi", CST, CST, 0)     # constants base (ADR in the real routine)

    def epilogue(self, g):
        for k in range(6):
            g.e("ldp", 19 + 2 * k, 20 + 2 * k, SP, 4096 + 16 * k)
        g.e("addi", SP, SP, 4096 + 112)
        g.e("ret")

    # ------------------------------------------------------------------ Fp12 product
    def build_mul(self, kconsts=None):
        g = self.new_gen(kconsts)
        p = self.p
        self.prologue(g, table=True)
        ys = [Val(Y, 32 * i - 0, "S", 4, 0, p - 1, "y%d" % i) for i in range(12)]   # y slot i: coef w^(i//2), re/im
        # y[m] (Fp2 coefficient of w^m) = (ys[2m], ys[2m+1]); Fp4 Y_j = (y[j], y[j+3])
        tslots = [Val(TA if i < 32 else TB, -512 + 32 * (i % 32), "S", 4, 0, p - 1, "t%d" % i) for i in range(54)]
        tcount = [0]
        def tnext():
            v = tslots[tcount[0]]; tcount[0] += 1; return v
        def yf2(m): return (ys[2 * m], ys[2 * m + 1])
        def yf4(j): return (yf2(j), yf2(j + 3))
        def f2prod(a):
            """Fp2 product of y-side a=(a0,a1) with the next three table operands -> (re, im) D values"""
            b0, b1, bs = tnext(), tnext(), tnext()
            s = g.sadd(a[0], a[1], reduce=False)
            P0 = g.mul(a[0], b0)
            P1 = g.mul(a[1], b1)
            g.mul_to_acc(s, bs)
            g.release(s)
            g.acc_add(-1, P0); g.acc_add(-1, P1)
            im = g.store_acc()
            g.lin([(1, P0), (-1, P1)])
            re = g.store_acc()
            g.release(P0); g.release(P1)
            return (re, im)
        def f4prod(X):
            """Fp4 product: X=(X0,X1) y-side Fp2 pairs -> (R0, R1) each an Fp2 pair of D values"""
            u0 = f2prod(X[0]); u1 = f2prod(X[1])
            s = (g.sadd(X[0][0], X[1][0], False), g.sadd(X[0][1], X[1][1], False))
            u2 = f2prod(s)
            for v in s: g.release(v)
            # R0 = u0 + xi u1 = (u0re + 9u1re - u1im, u0im + 9u1im + u1re); R1 = u2 - u0 - u1
            out = []
            g.lin([(1, u0[0]), (9, u1[0]), (-1, u1[1])]); r0re = g.store_acc()
            g.lin([(1, u0[1]), (9, u1[1]), (1, u1[0])]); r0im = g.store_acc()
            g.lin([(1, u2[0]), (-1, u0[0]), (-1, u1[0])]); r1re = g.store_acc()
            g.lin([(1, u2[1]), (-1, u0[1]), (-1, u1[1])]); r1im = g.store_acc()
            for v in u0 + u1 + u2: g.release(v)
            return ((r0re, r0im), (r1re, r1im))
        def ysum4(i, j):
            Yi, Yj = yf4(i), yf4(j)
            return tuple(tuple(g.sadd(Yi[a][b], Yj[a][b], True) for b in range(2)) for a in range(2))
        V4 = {}
        for j in range(3):
            V4["v%d" % j] = f4prod(yf4(j))
        for (i, j) in ((1, 2), (0, 1), (0, 2)):
            X = ysum4(i, j)
            V4["t%d%d" % (i, j)] = f4prod(X)
            for a in X:
                for b in a: g.release(b)
        # symbolic Fp4 combination: each Fp4 value is dict over (name, part, comp) -> coef
        def F(name): return {(name, a, b): 1 for a in range(2) for b in range(2)}
        def f4(name): return [[{(name, a, b): 1} for b in range(2)] for a in range(2)]
        def add(x, y, c=1):
            return [[{k: x[a][b].get(k, 0) + c * y[a][b].get(k, 0) for k in set(x[a][b]) | set(y[a][b])} for b in range(2)] for a in range(2)]
        def xi2(e):   # xi on an Fp2 symbolic pair (re, im)
            re, im = e
            nre = {k: 9 * re.get(k, 0) - im.get(k, 0) for k in set(re) | set(im)}
            nim = {k: 9 * im.get(k, 0) + re.get(k, 0) for k in set(re) | set(im)}
            return [nre, nim]
        def mul_s(x):  # s*(A0 + A1 s) = xi A1 + A0 s
            return [xi2(x[1]), x[0]]
        v0, v1, v2, t12, t01, t02 = (f4(n) for n in ("v0", "v1", "v2", "t12", "t01", "t02"))
        Z0 = add(v0, mul_s(add(add(t12, v1, -1), v2, -1)))
        Z1 = add(add(add(t01, v0, -1), v1, -1), mul_s(v2))
        Z2 = add(add(add(t02, v0, -1), v1, 1), v2, -1)
        Z = [Z0, Z1, Z2]
        outs = []
        for j in range(3):
            for a in range(2):          # Z_j = (coef w^j, coef w^(j+3))
                m = j + 3 * a
                for b in range(2):
                    terms = [(c, V4[n][pa][pb]) for (n, pa, pb), c in Z[j][a][b].items() if c]
                    g.lin(terms)
                    g.acc_add_k()
                    dest = Val(OUT, 32 * (2 * m + b), "S", 4, 0, p - 1)
                    g.redc_canon(dest)
        for n in V4:
            for a in V4[n]:
                for v in a: g.release(v)
        self.epilogue(g)
        return g

    # ------------------------------------------------------------------ cyclotomic squaring
    def build_sqr(self, kconsts=None):
        g = self.new_gen(kconsts)
        p = self.p
        self.prologue(g)
        ys = [Val(Y, 32 * i, "S", 4, 0, p - 1) for i in range(12)]
        def yf2(m): return (ys[2 * m], ys[2 * m + 1])
        def fp4sq(a, b):
            """products for (a + b s)^2: returns X0=(re,im), X1=(re,im) as D values"""
            def f2sq_prods(x, kp):
                s = g.sadd(x[0], x[1], False)
                d = g.ssub_kp(x[0], x[1], kp)
                Pm = g.mul(s, d); g.release(s); g.release(d)
                Qm = g.mul(x[0], x[1])
                return Pm, Qm          # x^2 = (Pm, 2 Qm)
            Pa, Qa = f2sq_prods(a, 1)
            Pb, Qb = f2sq_prods(b, 1)
            sab = (g.sadd(a[0], b[0], False), g.sadd(a[1], b[1], False))
            Ps, Qs = f2sq_prods(sab, 2)
            for v in sab: g.release(v)
            g.lin([(1, Pa), (9, Pb), (-2, Qb)]); x0re = g.store_acc()
            g.lin([(1, Pb), (2, Qa), (18, Qb)]); x0im = g.store_acc()
            g.lin([(1, Ps), (-1, Pa), (-1, Pb)]); x1re = g.store_acc()
            g.lin([(1, Qs), (-1, Qa), (-1, Qb)]); h1im = g.store_acc()   # X1im = 2*h1im
            for v in (Pa, Qa, Pb, Qb, Ps, Qs): g.release(v)
            return (x0re, x0im), (x1re, h1im)
        A0, A1 = fp4sq(yf2(0), yf2(3))
        B0, B1 = fp4sq(yf2(1), yf2(4))
        C0, C1 = fp4sq(yf2(2), yf2(5))
        def emit(m, b, terms, cterm):
            g.lin(terms, highs=[cterm])
            g.acc_add_k()
            g.redc_canon(Val(OUT, 32 * (2 * m + b), "S", 4, 0, p - 1))
        # c0' = 3A0 - 2c0 ; c3' = 3A1 + 2c3 ; c1' = 3 xi C1 + 2c1 ; c4' = 3C0 - 2c4 ; c2' = 3B0 - 2c2 ; c5' = 3B1 + 2c5
        # X1 = (x1re, 2*h1im)
        for (m, X, sgn) in ((0, A0, -1), (4, C0, -1), (2, B0, -1)):
            emit(m, 0, [(1, X[0]), (2, X[0])], (2 * sgn, ys[2 * m]))
            emit(m, 1, [(1, X[1]), (2, X[1])], (2 * sgn, ys[2 * m + 1]))
        for (m, X) in ((3, A1), (5, B1)):
            emit(m, 0, [(1, X[0]), (2, X[0])], (2, ys[2 * m]))
            emit(m, 1, [(1, X[1]), (5, X[1])], (2, ys[2 * m + 1]))      # 6*h1im = 3*X1im
        # c1' = 3 xi C1 + 2 c1: xi(x, y) = (9x - y, 9y + x), with x = x1re, y = 2 h1im
        emit(1, 0, [(1, C1[0]), (26, C1[0]), (-6, C1[1])], (2, ys[2]))  # 27 x - 6 h
        emit(1, 1, [(1, C1[0]), (2, C1[0]), (54, C1[1])], (2, ys[3]))   # 3x + 54 h
        for X in (A0, A1, B0, B1, C0, C1):
            for v in X: g.release(v)
        self.epilogue(g)
        return g

def f2pow(a, e, p):
    r = (1, 0)
    while e:
        if e & 1: r = ((r[0]*a[0]-r[1]*a[1]) % p, (r[0]*a[1]+r[1]*a[0]) % p)
        a = ((a[0]*a[0]-a[1]*a[1]) % p, (2*a[0]*a[1]) % p); e >>= 1
    return r

def key_consts(p):
    d1 = f2pow((9, 1), (p - 1) // 6, p)
    z = f2pow((9, 1), (p * p - 1) // 6, p)
    assert z[1] == 0
    return {"d1r": d1[0] * RM % p, "d1i": d1[1] * RM % p, "zeta": z[0] * RM % p}

def build_keys(self, kconsts=None):
    """kappa = argmin over the 12 stored keys; writes kappa to [OUT]."""
    g = self.new_gen(kconsts)
    p = self.p
    kc = key_consts(p)
    base = max(g.cvals) + 16
    for name, v in kc.items():
        g.cst[name] = base
        for i in range(4): g.cvals[base + 8 * i] = (v >> (64 * i)) & M64
        base += 32
    assert base <= WHI + 8
    self.prologue(g)
    a = Val(Y, 64, "S", 4, 0, p - 1)
    b = Val(Y, 96, "S", 4, 0, p - 1)
    cv = {n: Val(CST, g.cst[n], "S", 4, 0, p - 1) for n in kc}
    P = g.mul(a, cv["d1r"])
    g.mul_to_acc(b, cv["d1i"])
    g.acc_add(1, P); g.release(P)
    Bv = g.alloc("S", 4, 0, p - 1); g.redc_canon(Bv)
    g.mul_to_acc(a, cv["zeta"]); za = g.alloc("S", 4, 0, p - 1); g.redc_canon(za)
    g.mul_to_acc(Bv, cv["zeta"]); zb = g.alloc("S", 4, 0, p - 1); g.redc_canon(zb)
    MN = A[0:4]; IDX = A[4]; CND = A[5:9]; T = V[0:4]; PP = V[4:8]; MSK = V[8]; TMP = X1
    def load_to(v, regs):
        g.ld2(regs[0], regs[1], v.base, v.off); g.ld2(regs[2], regs[3], v.base, v.off + 16)
    def modsub(x, y, regs):        # regs <- (x - y) mod p, x may be None (= 0)
        if x is not None: load_to(x, regs)
        load_to(y, T)
        src = regs if x is not None else [ZR] * 4
        g.e("subs", regs[0], src[0], T[0])
        for j in (1, 2, 3): g.e("sbcs", regs[j], src[j], T[j])
        g.e("csetm", MSK, "cc")
        g.ldc(PP[0], PP[1], "p0"); g.ldc(PP[2], PP[3], "p2")
        for j in range(4): g.e("and", PP[j], PP[j], MSK)
        g.e("adds", regs[0], regs[0], PP[0])
        for j in (1, 2): g.e("adcs", regs[j], regs[j], PP[j])
        g.e("adc", regs[3], regs[3], PP[3])
    def six(x, zx):
        return [("v", x, None), ("v", zx, None), ("s", zx, x), ("s", None, x), ("s", None, zx), ("s", x, zx)]
    A6, B6 = six(a, za), six(Bv, zb)
    keys = []
    for j in range(6): keys += [A6[j], B6[j]]
    load_to(a, MN); g.e("movi", IDX, 0)
    for k in range(1, 12):
        kind, x, y = keys[k]
        if kind == "v": load_to(x, CND)
        else: modsub(x, y, CND)
        g.e("subs", ZR, CND[0], MN[0])
        for j in (1, 2, 3): g.e("sbcs", ZR, CND[j], MN[j])
        for j in range(4): g.e("csel", MN[j], MN[j], CND[j], "cs")
        g.e("movi", TMP, k)
        g.e("csel", IDX, IDX, TMP, "cs")
    g.e("stp", IDX, ZR, OUT, 0)
    for v in (Bv, za, zb): g.release(v)
    self.epilogue(g)
    return g
Lazy.build_keys = build_keys


def build_two_pass(fn):
    g = fn()
    kn = getattr(g, "kneed", None)
    return fn({kn: "K0"}) if kn else g


def table_operands(M, q):
    """the 54 canonical table-side operands of the Fp12 product, in the order build_mul consumes them
    (precomputed once per table entry at setup; M is in Montgomery form)"""
    ops = []

    def f4(j):
        return (M[j], M[j + 3])

    def f2ops(b):
        ops.extend([b[0], b[1], (b[0] + b[1]) % q])

    def f4ops(W):
        f2ops(W[0])
        f2ops(W[1])
        f2ops(((W[0][0] + W[1][0]) % q, (W[0][1] + W[1][1]) % q))
    for j in range(3):
        f4ops(f4(j))
    for (i, j) in ((1, 2), (0, 1), (0, 2)):
        Wi, Wj = f4(i), f4(j)
        f4ops(tuple(((Wi[a][0] + Wj[a][0]) % q, (Wi[a][1] + Wj[a][1]) % q) for a in range(2)))
    return ops


def lazy_env(g, y, M=None):
    mem = {}
    regs = [0] * 33
    for k, v in g.cvals.items():
        mem[0x100000 + 512 + k] = v
    regs[CST] = 0x100000 + 512
    regs[SP], regs[Y], regs[OUT], regs[TA] = 0x200000, 0x300000, 0x400000, 0x500000 + 512
    for i in range(33):
        if i not in (CST, SP, Y, OUT, TA, ZR):
            regs[i] = rng.getrandbits(64)       # garbage in every other register

    def put(addr, v):
        for i in range(4):
            mem[addr + 8 * i] = (v >> (64 * i)) & M64
    for m in range(6):
        for b in range(2):
            put(regs[Y] + 32 * (2 * m + b), y[m][b] * RMONT % p)
    if M is not None:
        Mm = [tuple(c * RMONT % p for c in M[m]) for m in range(6)]
        for i, v in enumerate(table_operands(Mm, p)):
            put(regs[TA] - 512 + 32 * i, v)
    return mem, regs


def lazy_out(mem, regs):
    raw = []
    for m in range(6):
        c = []
        for b in range(2):
            a = regs[OUT] + 32 * (2 * m + b)
            c.append(sum(mem[a + 8 * i] << (64 * i) for i in range(4)))
        raw.append(tuple(c))
    return raw


def part_L(cases, gt_cases):
    Rinv = pow(RMONT, -1, p)
    Lz = Lazy(p)
    gm, gs, gk = build_two_pass(Lz.build_mul), build_two_pass(Lz.build_sqr), build_two_pass(Lz.build_keys)
    P = Tower(p)
    g = profile_to_tower(G_PROFILE, p)
    res = {}
    edge = [0, 1, 2, p - 1, p - 2, (p >> 1), RMONT % p, (1 << 192) % p]
    fails = noncanon = 0
    counts = set()
    for t in range(cases):
        if t % 3 == 0:
            y = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
        elif t % 3 == 1:
            y = tuple((rng.choice(edge), rng.choice(edge)) for _ in range(6))
        else:
            y = tuple((p - 1 - rng.getrandbits(40), p - 1 - rng.getrandbits(40)) for _ in range(6))
        M = tuple((rng.choice([p - 1, rng.randrange(p)]), rng.choice([0, p - 1, rng.randrange(p)])) for _ in range(6))
        mem, regs = lazy_env(gm, y, M)
        counts.add(("fp12_mul", run(gm.code, mem, regs)))
        raw = lazy_out(mem, regs)
        noncanon += any(c >= p for cc in raw for c in cc)
        fails += tuple(tuple(c * Rinv % p for c in cc) for cc in raw) != P.mul(y, M)
    res["fp12_mul"] = {"cases": cases, "failures": fails, "noncanonical_outputs": noncanon}
    fs = fk = nc = 0
    for t in range(gt_cases):
        y = P.pow(g, rng.randrange(1, r))
        mem, regs = lazy_env(gs, y)
        counts.add(("cyclotomic_square", run(gs.code, mem, regs)))
        raw = lazy_out(mem, regs)
        nc += any(c >= p for cc in raw for c in cc)
        fs += tuple(tuple(c * Rinv % p for c in cc) for cc in raw) != P.mul(y, y)
        for yy in (y, (y[0], (rng.choice(edge), rng.choice(edge))) + tuple(y[2:])):
            kst = [v * RMONT % p for v in P.keys(yy)]
            mem, regs = lazy_env(gk, yy)
            counts.add(("keys_argmin", run(gk.code, mem, regs)))
            fk += mem[regs[OUT]] != kst.index(min(kst))
    res["cyclotomic_square"] = {"GT_cases": gt_cases, "failures": fs, "noncanonical_outputs": nc}
    res["keys_argmin"] = {"cases": 2 * gt_cases, "failures": fk}
    ex = {}
    for k, n in counts:
        ex.setdefault(k, set()).add(n)
    res["executed_instructions"] = {k: sorted(v) for k, v in ex.items()}
    res["straight_line_lengths"] = {"fp12_mul": len(gm.code), "cyclotomic_square": len(gs.code), "keys_argmin": len(gk.code)}
    res["peak_scratch_slots"] = {"fp12_mul": gm.peak, "cyclotomic_square": gs.peak}
    res["K_constant_multiple_of_p"] = {"fp12_mul": max(gm.kconsts) // p, "cyclotomic_square": max(gs.kconsts) // p}
    hist = {}
    for i in gm.code:
        hist[i[0]] = hist.get(i[0], 0) + 1
    res["opcode_histogram_fp12_mul"] = hist
    if SHOW_LISTINGS:
        res["listing_fp12_mul"] = [" ".join(str(a) for a in i) for i in gm.code]
    return res, len(gm.code), len(gs.code), len(gk.code)


# =============================================================== Part G: gate and expected-work ledger
CALLS_MUL = 8          # two routine call sites per multiply step (BL + pointer set-ups)
CALLS_SQ = 4           # one routine call site per squaring step
COMMON = 300           # per-step bookkeeping, itemised 137 in machine.md section 4 (charged with margin)
R_TABLE = 1024
S_SQ = Fraction(1, 2)
THETA_BITS, WALKS = 60, 2 ** 20


def ledger(mul_len, sq_len, key_len, allow_percent=5, cap_multiple=64):
    getcontext().prec = 80
    ln2 = Decimal(2).ln()
    lg = lambda v: Decimal(v).ln() / ln2
    dec = lambda fr: Decimal(fr.numerator) / Decimal(fr.denominator)

    def arctan_inv(n):
        x = Decimal(1) / n
        total, term, k, n2 = x, x, 1, n * n
        while abs(term) > Decimal(10) ** -75:
            term = -term / n2
            k += 2
            total += term / k
        return total
    PI = 4 * (4 * arctan_inv(5) - arctan_inv(239))
    mul_step = mul_len + key_len + CALLS_MUL + COMMON
    sq_step = sq_len + CALLS_SQ + COMMON
    C = 20000
    c_avg = dec((1 - S_SQ) * Fraction(mul_step, C) + S_SQ * Fraction(sq_step, C))
    c_mul = dec(Fraction(mul_step, C))
    N = (r - 1) // 12
    sum_q2 = S_SQ ** 2 + (1 - S_SQ) ** 2 / R_TABLE
    N_eff = Decimal(N) / dec(1 - sum_q2)
    eps = dec(Fraction(5, 12 * R_TABLE) + Fraction(250, R_TABLE ** 2))
    E_delay = Decimal(WALKS * 20) * Decimal(2) ** THETA_BITS
    model_mean = (PI * N_eff / 2).sqrt()
    ET = Decimal(100 + allow_percent) / 100 * model_mean
    S_cap = cap_multiple * N_eff.sqrt()
    markov = ET / S_cap
    escape_cap = eps * (ET + E_delay) / (S_cap + E_delay)
    fail = markov + Decimal(2) ** -8 + Decimal("0.001") + escape_cap
    succ = Decimal("0.97")
    K = int((Decimal("0.39") / succ * Decimal(2) ** 64).to_integral_value(rounding=ROUND_CEILING))
    gate = Decimal(K) / Decimal(2) ** 64
    walk = c_avg * ET * (1 + eps)
    delay = c_mul * E_delay
    NONWALK = Decimal(2) ** 66
    reports_expected = ET * Decimal(2) ** -THETA_BITS * (1 + eps) + Decimal(WALKS)
    run_ = walk + delay + NONWALK
    total = 1 + gate * run_
    every_run = c_mul * 2 * (S_cap + E_delay) + NONWALK * 2 ** 8
    claim = CLAIM
    return {
        "multiply_step_instr": mul_step, "squaring_step_instr": sq_step,
        "multiply_step_parts": {"fp12_mul": mul_len, "keys_argmin": key_len, "call_sites": CALLS_MUL, "common": COMMON},
        "squaring_step_parts": {"cyclotomic_square": sq_len, "call_site": CALLS_SQ, "common": COMMON},
        "multiply_step_units": str(c_mul), "squaring_step_units": str(Decimal(sq_step) / C),
        "average_step_units": str(c_avg),
        "log2_N_classes": float(lg(N)), "log2_N_effective": float(lg(N_eff)), "eps": float(eps),
        "log2_model_mean_steps": float(lg(model_mean)), "log2_charged_mean_steps": float(lg(ET)),
        "log2_step_cap": float(lg(S_cap)), "markov_bound_step_cap": float(markov),
        "all_step_counter_cap_bound": float(escape_cap),
        "failure_budget_past_gate": float(fail), "success_past_gate_lower": float(1 - fail),
        "success_past_gate_claimed": str(succ), "budget_ok": 1 - fail >= succ,
        "gate_K": str(K), "gate_probability": float(gate), "total_success": float(gate * succ),
        "success_ok": gate * succ >= Decimal("0.39"),
        "K_is_least": Decimal(K - 1) / Decimal(2) ** 64 * succ < Decimal("0.39"),
        "log2_expected_walk_units": float(lg(walk)), "log2_delay_units": float(lg(delay)),
        "log2_expected_reports": float(lg(reports_expected)),
        "nonwalk_bucket_covers_reports": reports_expected / 2 + Decimal(2) ** 60 <= NONWALK,
        "log2_expected_run_units": float(lg(run_)), "log2_expected_total_units": float(lg(total)),
        "claim": str(claim), "claim_ok": claim >= lg(total), "claim_minus_computed_bits": float(claim - lg(total)),
        "log2_every_run_bound_units": float(lg(every_run)),
    }



# =============================================================== main
CLAIM = Decimal("122.85")
out = {"seed": SEED}
lres, mul_len, sq_len, key_len = part_L(cfg.get("lazy_cases", 30), cfg.get("lazy_gt_cases", 4))
out["partL_lazy_routines"] = lres
out["partG_ledger"] = ledger(mul_len, sq_len, key_len, cfg.get("mean_allowance_percent", 5), cfg.get("cap_multiple", 64))
# the base package's experiment.py (submission 2ece9f04, unchanged) runs in the same process with its own
# input: tower identities, Frobenius/key/trace identities, production step rule in the real GT, toy runs
if cfg.get("base_run", True):
    import contextlib
    import io
    bpath = base / "experiment.py"
    bcfg = {"seed": cfg.get("base_seed", 20261111), "seconds": cfg.get("base_seconds", 9),
            "production_steps": cfg.get("base_production_steps", 1500)}
    buf = io.StringIO()
    saved = sys.stdin
    sys.stdin = io.StringIO(json.dumps(bcfg))
    try:
        with contextlib.redirect_stdout(buf):
            exec(compile(bpath.read_text(), str(bpath), "exec"), {"__name__": "base_experiment", "__file__": str(bpath)})
    finally:
        sys.stdin = saved
    bout = json.loads(buf.getvalue())
    bout.pop("partD_ledger", None)      # the base package's own 124.65 ledger, superseded by part G
    out["base_experiment_input"] = bcfg
    out["base_experiment"] = bout
out["elapsed_seconds"] = round(time.time() - T_START, 2)
out["scope"] = ("Generated routines executed instruction by instruction on a simulator of the scalar AArch64 subset "
                "they use and compared with big-integer arithmetic on the protected BN254 GT; ledger arithmetic; "
                "base identities, production step rule and toy-field runs. No production DLP is run.")
print(json.dumps(out, indent=1, default=str))
