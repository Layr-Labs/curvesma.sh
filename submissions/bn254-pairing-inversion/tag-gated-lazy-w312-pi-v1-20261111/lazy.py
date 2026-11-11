"""Gated tag-tracing rho on BN254 GT: lazy window-end product, tag lemma in its layout, toy mean, ledger.

Part L  the window-end routine priced by this package -- the dense Fp12 product Y0 * T_e with the
        table side stored as its 54 Karatsuba operands -- is generated as a straight-line scalar AArch64
        program (64-bit general registers only), executed instruction by instruction on a simulator and
        compared with big-integer arithmetic; static worst-case bounds are asserted while generating.
Part T  the tag certificate (proof.md section 2) checked in this package's coordinate layout on window
        bases produced by the generated routine itself.
Part H  toy-scale runs of the complete walk (8-adding, selector floor(8t/q), distinguished points,
        sequential walks, run to the first recovery) measuring the mean collision time against
        sqrt(pi kappa n / 2), kappa = 8/7 (premise H-mean).
Part G  the gate and expected-work ledger in exact decimal arithmetic.
Base    the unchanged experiment.py of submission 73463878 (tag-tracing package) is then executed
        in-process: tower checks, tag row identity, certificate on 260 cases, 120-step production walk,
        toy runs at the cap, listing recount.
"""
import json
import random
import sys
import time
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction
from pathlib import Path

T_START = time.time()
DP_BITS_PRODUCTION = 84   # distinguished iff t(Y)*2^84 < p (proof.md section 1); density below 2^-84 + 2^-254
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


class Tower:  # reference big-integer arithmetic
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
        AI = V[4:8]
        L = [V[8], X1, X2, X3]
        self.load_s(b, B)
        self.load_s(a, AI)
        ai = AI[0]
        for j in range(4):
            self.e("mul", A[j], ai, B[j])
        for j in range(4):
            self.e("umulh", L[j], ai, B[j])
        self.e("adds", A[1], A[1], L[0])
        self.e("adcs", A[2], A[2], L[1])
        self.e("adcs", A[3], A[3], L[2])
        self.e("adc", A[4], L[3], ZR)
        for i in (1, 2, 3):
            ai = AI[i]
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
        elif op == "adr":
            regs[ins[1]] = regs[ins[1]]          # the harness preloads the constant-pool address
        elif op == "ret" or op == "bl":
            pass
        else:
            raise RuntimeError("unsupported " + op)
        regs[ZR] = 0
    return n


GPR = set(range(31))


def validate(code):
    """Reject any instruction that is not an encodable scalar AArch64 instruction of the forms used.
    Returns a list of violations (empty when every instruction is legal)."""
    bad = []
    for n, ins in enumerate(code):
        op = ins[0]
        regs = []
        ok = True
        if op in ("ldp", "stp"):
            r0, r1, base, off = ins[1:5]
            regs = [r0, r1]
            ok = off % 8 == 0 and -512 <= off <= 504 and (base in GPR or base == SP) and base != 18
            if op == "ldp":
                ok = ok and r0 != r1 and r0 in GPR and r1 in GPR
            else:
                ok = ok and all(x in GPR or x == ZR for x in (r0, r1))
        elif op in ("ldr", "str"):
            r0, base, off = ins[1:4]
            regs = [r0]
            ok = (base in GPR or base == SP) and ((off % 8 == 0 and 0 <= off <= 32760) or -256 <= off <= 255)
        elif op in ("mul", "umulh", "and", "udiv"):
            ok = ins[1] in GPR and all(x in GPR or x == ZR for x in ins[2:4])
        elif op in ("adds", "adcs", "sbcs"):
            ok = (ins[1] in GPR or ins[1] == ZR) and all(x in GPR or x == ZR for x in ins[2:4])
        elif op in ("adc", "sbc", "add"):
            ok = ins[1] in GPR and all(x in GPR or x == ZR for x in ins[2:4])
        elif op == "subs":
            if isinstance(ins[3], tuple):
                ok = (ins[1] in GPR or ins[1] == ZR) and ins[2] in GPR and 0 <= ins[3][1] <= 4095
            else:
                ok = (ins[1] in GPR or ins[1] == ZR) and all(x in GPR or x == ZR for x in ins[2:4])
        elif op == "csel":
            ok = ins[1] in GPR and all(x in GPR or x == ZR for x in ins[2:4]) and ins[4] in ("cs", "cc")
        elif op == "csetm":
            ok = ins[1] in GPR and ins[2] in ("cs", "cc")
        elif op == "extr":
            ok = ins[1] in GPR and all(x in GPR or x == ZR for x in ins[2:4]) and 0 <= ins[4] <= 63
        elif op in ("lsl", "lsr"):
            ok = ins[1] in GPR and (ins[2] in GPR or ins[2] == ZR) and 1 <= ins[3] <= 63
        elif op == "addi":
            d, src, imm = ins[1:4]
            enc = abs(imm) <= 4095 or (imm % 4096 == 0 and abs(imm) // 4096 <= 4095)
            ok = (d in GPR or d == SP) and (src in GPR or src == SP) and enc
        elif op == "movi":
            ok = ins[1] in GPR and 0 <= ins[2] <= 65535
        elif op in ("adr", "ret", "bl"):
            ok = True
        else:
            ok = False
        if ok and any(x == 18 for x in ins[1:] if isinstance(x, int)):
            ok = False
        if not ok:
            bad.append((n, " ".join(str(x) for x in ins)))
    return bad


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
        # frame of 4,208 bytes: saves of x19-x30 at [sp, 0..95], four 1 KiB scratch windows above them
        g.e("addi", SP, SP, -4096)          # SUB sp, sp, #1, LSL #12
        g.e("addi", SP, SP, -112)           # SUB sp, sp, #112
        for k in range(6):
            g.e("stp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k)
        for i, b in enumerate(SCR):
            g.e("addi", b, SP, 112 + 512 + 1024 * i)     # window i covers [sp+112+1024i, sp+112+1024i+1016]
        if table:
            g.e("addi", TB, TA, 1024)
        g.e("adr", CST)                      # ADR x4, constant pool

    def epilogue(self, g):
        for k in range(6):
            g.e("ldp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k)
        g.e("addi", SP, SP, 112)
        g.e("addi", SP, SP, 4096)
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


def part_L(cases):
    Rinv = pow(RMONT, -1, p)
    gm = build_two_pass(Lazy(p).build_mul)
    P = Tower(p)
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
        counts.add(run(gm.code, mem, regs))
        raw = lazy_out(mem, regs)
        noncanon += any(c >= p for cc in raw for c in cc)
        fails += tuple(tuple(c * Rinv % p for c in cc) for cc in raw) != P.mul(y, M)
    hist = {}
    for i in gm.code:
        hist[i[0]] = hist.get(i[0], 0) + 1
    viol = validate(gm.code)
    res = {"encoding_violations": len(viol), "first_violations": viol[:5],
           "cases": cases, "failures": fails, "noncanonical_outputs": noncanon,
           "executed_instructions": sorted(counts), "straight_line_length": len(gm.code),
           "peak_scratch_slots": gm.peak, "K_multiple_of_p": max(gm.kconsts) // p, "opcode_histogram": hist}
    if SHOW_LISTINGS:
        # bounded so the run stays inside the runner's output allowance: the first and last 1,000 instructions
        # in full, plus the SHA-256 of the complete listing text (regenerated deterministically by part L)
        import hashlib
        lines = [" ".join(str(a) for a in i) for i in gm.code]
        res["listing_sha256"] = hashlib.sha256("\n".join(lines).encode()).hexdigest()
        res["listing_head"] = lines[:1000]
        res["listing_tail"] = lines[-1000:]
        res["listing_omitted_middle_instructions"] = max(0, len(lines) - 2000)
    return res, gm


# =============================================================== Part T: tag certificate in this layout
def t_of(Z):
    """protected-basis coefficient of w^0: tower coefficient 0 is (c0 + 9 c6) + c6 u"""
    return (Z[0][0] - 9 * Z[0][1]) % p


def row(P, T):
    """v_k with t(Y*T) = sum_k y_k v_k, coordinate k = 2m + b of the stored layout"""
    v = []
    for k in range(12):
        e = [[0, 0] for _ in range(6)]
        e[k // 2][k % 2] = 1
        v.append(t_of(P.mul(tuple(tuple(x) for x in e), T)))
    return v


def fwords(v):
    F = []
    for k in range(12):
        for i in range(4):
            u = v[k] * pow(2, 64 * i, p) * pow(RMONT, -1, p) % p
            Fv = (u << 128) // p
            F.append((Fv & M64, Fv >> 64))
    return F


def part_T(gm, cases):
    P = Tower(p)
    g = profile_to_tower(G_PROFILE, p)
    worst, bad, cert, sel_ok = 0, 0, 0, 0
    for c in range(cases):
        Y = P.pow(g, rng.randrange(1, r)) if c % 2 == 0 else tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
        Tp = P.pow(g, rng.randrange(1, r))
        mem, regs = lazy_env(gm, Y, Tp)               # window base Y0 = Y * Tp by the generated routine
        run(gm.code, mem, regs)
        raw = lazy_out(mem, regs)
        Y0 = tuple(tuple(c_ * pow(RMONT, -1, p) % p for c_ in cc) for cc in raw)
        T = P.pow(g, rng.randrange(1, r))             # current table product
        F = fwords(row(P, T))
        A = 0
        j = 0
        for m in range(6):
            for b in range(2):
                for i in range(4):
                    limb = (raw[m][b] >> (64 * i)) & M64
                    lo, hi = F[j]
                    A = (A + ((limb * lo) >> 64) + limb * hi) & M64
                    j += 1
        t = t_of(P.mul(Y0, T))
        off = ((t << 64) // p - A) % (1 << 64)        # floor(tau) - A; tau - A lies in [off, off + 1)
        worst = max(worst, off)
        bad += off >= 96
        if ((A + 96) % (1 << 64)) % (1 << 61) > 96:
            cert += 1
            sel_ok += (A >> 61) == (8 * t) // p
    return {"cases": cases, "offsets_outside_[0,96)": bad, "largest_offset": worst,
            "certified": cert, "certified_selector_equals_floor_8t_over_p": sel_ok}

# =============================================================== Part H: toy mean collision time (8-adding walk)
import math as _m


def _is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def toy_group(bits, coins):
    """a fresh 61-bit prime q = 2 n h + 1 with prime n of the given size; generator of the order-n subgroup"""
    while True:
        n = coins.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not _is_prime(n):
            continue
        for _ in range(200):
            h = coins.getrandbits(60 - bits) | (1 << (59 - bits))
            q = 2 * n * h + 1
            if _is_prime(q):
                while True:
                    g = pow(coins.randrange(2, q - 1), 2 * h, q)
                    if g != 1:
                        return q, n, g


def toy_run(bits, dbits, coins, m=8):
    """one target, the walk of proof.md section 3 at toy scale: s = floor(m t / q), distinguished iff
    t 2^dbits < q, sequential walks to their distinguished point, run to the first recovery.
    Returns (touch steps, total steps to recovery, walks, n)."""
    q, n, g = toy_group(bits, coins)
    x = coins.randrange(1, n)
    z = pow(g, x, q)
    cs = [coins.randrange(n) for _ in range(m)]
    ds = [coins.randrange(n) for _ in range(m)]
    Ms = [pow(g, cs[i], q) * pow(z, ds[i], q) % q for i in range(m)]
    lim = q >> dbits
    store = {}
    seen = set()
    steps, touch, walks = 0, None, 0
    len_cap = 64 << dbits
    while True:
        walks += 1
        a, b = coins.randrange(n), coins.randrange(n)
        y = pow(g, a, q) * pow(z, b, q) % q
        trail = []
        k = 0
        while y >= lim and k < len_cap:
            if touch is None:
                trail.append(y)
                if y in seen:
                    touch = steps
            s = y * m // q
            y = y * Ms[s] % q
            a = (a + cs[s]) % n
            b = (b + ds[s]) % n
            steps += 1
            k += 1
        if touch is None:
            seen.update(trail)
            if y in seen:
                touch = steps
            seen.add(y)
        if k >= len_cap:
            continue
        old = store.get(y)
        if old is not None and old[1] != b:
            xr = (old[0] - a) * pow(b - old[1], -1, n) % n
            assert pow(g, xr, q) == z
            return (touch if touch is not None else steps), steps, walks, n
        store[y] = (a, b)


def part_H(settings, coins, budget):
    """mean of touch/sqrt(pi kappa n/2) and of recovery/sqrt(pi kappa n/2), kappa = 8/7, per setting"""
    t0 = time.time()
    out = []
    kappa = 8 / 7
    for bits, dbits, trials in settings:
        RT, RS, Wk = [], [], []
        while len(RT) < trials and time.time() - t0 < budget:
            a, b, c, n = toy_run(bits, dbits, coins)
            unit = (_m.pi * kappa * n / 2) ** 0.5
            RT.append(a / unit)
            RS.append(b / unit)
            Wk.append(c)
        if not RT:
            continue
        k = len(RT)
        mt, ms = sum(RT) / k, sum(RS) / k
        st = (sum((v - mt) ** 2 for v in RT) / max(1, k - 1)) ** 0.5 / k ** 0.5
        ss = (sum((v - ms) ** 2 for v in RS) / max(1, k - 1)) ** 0.5 / k ** 0.5
        out.append({"n_bits": bits, "dp_bits": dbits, "runs": k,
                    "mean_touch_over_model": round(mt, 4), "se": round(st, 4),
                    "mean_recovery_steps_over_model": round(ms, 4), "se_recovery": round(ss, 4),
                    "mean_walks": round(sum(Wk) / k, 2)})
    return out

# =============================================================== Part G: gate and expected-work ledger
STEP = 273             # one certified step: ttstep.disasm.txt (base package, recounted by experiment.py)
TT_ENTRY = 65          # tt_run entry/exit per call
DRIVER = 160           # driver per call with a three-word global step counter (base: two words, 128)
LABELS = 2 * 34 + 20   # two fr_add (tower.disasm.txt) and rebase glue
WINDOW = 312
FIXED = 2 ** 110       # capped work outside steps and window ends (machine.md), every run


def ledger(mul_len, allow_percent=5, cap_bits=133):
    getcontext().prec = 80
    ln2 = Decimal(2).ln()
    lg = lambda v: Decimal(v).ln() / ln2

    def arctan_inv(n):
        x = Decimal(1) / n
        total, term, k, n2 = x, x, 1, n * n
        while abs(term) > Decimal(10) ** -75:
            term = -term / n2
            k += 2
            total += term / k
        return total
    PI = 4 * (4 * arctan_inv(5) - arctan_inv(239))
    rebase = mul_len + LABELS
    wend_exact = rebase + TT_ENTRY + DRIVER
    WEND = 9216
    assert wend_exact <= WEND
    kappa = Decimal(8) / 7
    model = (PI * kappa * Decimal(r) / 2).sqrt()                     # sqrt(pi kappa r / 2)
    E_touch = Decimal(100 + allow_percent) / 100 * model             # H-mean
    E_steps = E_touch + Decimal(2) ** 85                             # plus the touching walk's way to its DP
    S_cap = Decimal(2) ** cap_bits
    E_walks = E_steps / Decimal(2) ** 84 * 2 + 1                     # generous: twice the DP count, plus one
    E_exact = E_steps * 97 / Decimal(2) ** 61
    # every walk: no DP within 2^90 = 64 x 2^84 steps (e^-64) or a self-cycle first (E[len^2]/(2r), len ~ geometric
    # with mean 2^84, E[len^2] <= 2 x 2^168); union over the 2^52-walk cap
    walk_union = Decimal(2) ** 52 * (Decimal(-64).exp() + 2 * Decimal(2) ** 168 / (2 * Decimal(r)))
    losses = {
        "step_cap_markov": E_steps / S_cap,
        "walk_cap_2^52_markov": E_walks / Decimal(2) ** 52,
        "exact_meter_2^90_markov": E_exact / Decimal(2) ** 90,
        "store_probe_meter_2^70_markov": Decimal(2) ** -20,
        "walk_length_cap_and_self_cycles": Decimal(2) ** -33,   # 2^52 x (e^-64 + 2^-85.6) < 2^-33
        "degeneracy_and_sampling": Decimal(2) ** -200,
    }
    fail = sum(losses.values())
    succ = Decimal("0.97")
    K = int((Decimal("0.39") / succ * Decimal(2) ** 64).to_integral_value(rounding=ROUND_CEILING))
    gate = Decimal(K) / Decimal(2) ** 64
    R_run = STEP * E_steps + WEND * (E_steps / WINDOW + E_walks + 1) + FIXED
    W_run = R_run / 20000
    total = 1 + gate * W_run
    every_run = (STEP * S_cap + WEND * (S_cap / WINDOW + Decimal(2) ** 52 + 1) + FIXED) / 20000
    entries = _m.comb(WINDOW + 8, 8)
    table_bytes = entries * 2624
    store_bytes = 2 ** 54 * 96
    peak_bytes = table_bytes + store_bytes + 2 ** 30
    assert entries == 2496767727736920
    assert peak_bytes < 2 ** 63
    return {
        "rebase_instr": rebase, "window_end_instr": wend_exact, "window_end_charged": WEND,
        "window_steps": WINDOW, "table_entries": entries, "table_bytes": table_bytes,
        "store_bytes": store_bytes, "peak_bytes_bound": peak_bytes, "peak_below_2^63": peak_bytes < 2 ** 63,
        "per_step_instr": f"{STEP * WINDOW + WEND}/{WINDOW}",
        "log2_model_mean_touch": float(lg(model)), "log2_charged_mean_steps": float(lg(E_steps)), "log2_expected_exact_evaluations": float(lg(E_exact)),
        "log2_walk_cap_union_bound": float(lg(walk_union)),
        "log2_step_cap": cap_bits, "losses": {k: float(v) for k, v in losses.items()},
        "failure_budget_past_gate": float(fail), "success_past_gate_lower": float(1 - fail),
        "budget_ok": 1 - fail >= succ, "gate_K": str(K), "gate_probability": float(gate),
        "total_success": float(gate * succ), "success_ok": gate * succ >= Decimal("0.39"),
        "K_is_least": Decimal(K - 1) / Decimal(2) ** 64 * succ < Decimal("0.39"),
        "log2_expected_run_units": float(lg(W_run)), "log2_expected_total_units": float(lg(total)),
        "claim": str(CLAIM), "claim_ok": CLAIM >= lg(total), "claim_minus_computed_bits": float(CLAIM - lg(total)),
        "log2_every_run_bound_units": float(lg(every_run)),
    }


# =============================================================== main
CLAIM = Decimal("119.94")
out = {"seed": SEED}
TIMING = {}
GUARD_SECONDS = 12          # optional parts are skipped (and reported) once this much time has elapsed
SKIPPED = []


def time_left(part):
    if time.time() - T_START > GUARD_SECONDS:
        SKIPPED.append(part)
        return False
    return True
_t = time.time()
lres, gm = part_L(cfg.get("lazy_cases", 30))
out["partL_window_end_product"] = lres
TIMING["L"] = round(time.time() - _t, 2); _t = time.time()
if time_left("T"):
    out["partT_tag_certificate_in_this_layout"] = part_T(gm, cfg.get("tag_cases", 6))
TIMING["T"] = round(time.time() - _t, 2); _t = time.time()
if time_left("H"):
    out["partH_toy_mean_collision_time"] = part_H([(16, 3, cfg.get("toy_runs", 600)), (20, 4, cfg.get("toy_runs", 600) // 2),
                                                   (24, 5, cfg.get("toy_runs", 600) // 8)], random.Random(SEED + 1), cfg.get("toy_seconds", 0.5))
TIMING["H"] = round(time.time() - _t, 2); _t = time.time()
out["partG_ledger"] = ledger(lres["straight_line_length"], cfg.get("mean_allowance_percent", 5), cfg.get("cap_bits", 133))
DECLARED_LEDGER = out["partG_ledger"]
if (cfg.get("mean_allowance_percent", 5), cfg.get("cap_bits", 133)) != (5, 133):
    # sensitivity inputs: partG_ledger shows their effect; the claim is checked against the declared parameters
    DECLARED_LEDGER = ledger(lres["straight_line_length"], 5, 133)
    out["partG_ledger_declared_parameters"] = DECLARED_LEDGER
    out["partG_ledger_note"] = ("partG_ledger uses the requested sensitivity inputs; the claim, the stated figures "
                                "and part C refer to the declared parameters (mean allowance 5%, cap 2^133)")
TIMING["G"] = round(time.time() - _t, 2); _t = time.time()
if cfg.get("base_run", True) and time_left("base"):
    import contextlib
    import io
    bpath = base / "experiment.py"
    bcfg = {"seed": cfg.get("base_seed", 20261111), "trials": cfg.get("base_trials", 16)}
    buf = io.StringIO()
    saved = sys.stdin
    sys.stdin = io.StringIO(json.dumps(bcfg))
    bns = {"__name__": "base_experiment", "__file__": str(bpath)}
    try:
        with contextlib.redirect_stdout(buf):
            exec(compile(bpath.read_text(), str(bpath), "exec"), bns)
    finally:
        sys.stdin = saved
    txt = buf.getvalue()
    try:
        out["base_experiment"] = json.loads(txt)
    except ValueError:
        out["base_experiment_text"] = txt[-20000:]
    out["base_experiment_input"] = bcfg
    # Production predicate (this package): Y is distinguished iff t(Y)*2^84 < p.  The base experiment's
    # full-width walk instantiates the same predicate family with its own parameter 80; here the same
    # reference solve loop is run on the full-width group with the production parameter 84, and every
    # step's distinguished flag is checked against the definition t*2^84 < p.
    def production_predicate_walk(dbits=DP_BITS_PRODUCTION, steps_wanted=cfg.get("production_steps", 16)):
        ref, G, g_, r_, p_ = bns["ref"], bns["G"], bns["g"], bns["r"], bns["p"]
        rngw = random.Random(SEED + 84)
        xw = rngw.randrange(1, r_)
        zw = G.power(g_, xw)
        st = {"steps": 0, "dp_flag_checks": 0, "dp_true": 0, "selector_checks": 0}
        def trace(table, Y0, e, s_, dp, a, b):
            t = G.coord0(G.mul(Y0, e.T))
            assert dp == ((t << dbits) < p_), "distinguished flag differs from t*2^dbits < p"
            assert s_ == ref.exact_selector(G, t)
            st["dp_flag_checks"] += 1; st["selector_checks"] += 1; st["dp_true"] += int(dp)
        got, steps, outcome, stats = ref.solve(G, zw, 5, dbits, 3, steps_wanted, 1, 1 << 30, rngw, trace)
        st.update({"steps": steps, "outcome": outcome, "dbits": dbits, "exact_paths": stats["exact_paths"],
                   "rebases": stats["rebases"]})
        return st
    TIMING["base"] = round(time.time() - _t, 2); _t = time.time()
    try:
        if time_left("P"):
            out["partP_production_predicate_full_width"] = production_predicate_walk()
    except Exception as exc:  # reported, not hidden
        out["partP_production_predicate_full_width"] = {"error": repr(exc)}
TIMING["P"] = round(time.time() - _t, 2); _t = time.time()
# Part C: cross-document parameter consistency (consistency.py), on this run's ledger and part P/L outputs
cns = {"__name__": "consistency", "__file__": str(base / "consistency.py")}
exec(compile((base / "consistency.py").read_text(), str(base / "consistency.py"), "exec"), cns)
try:
    out["partC_consistency"] = cns["run"](base, DECLARED_LEDGER, out["partL_window_end_product"],
                                          out.get("partP_production_predicate_full_width", {}), out.get("base_experiment"))
except Exception as exc:  # reported, not hidden
    out["partC_consistency"] = {"all_ok": False, "error": repr(exc)}
TIMING["C"] = round(time.time() - _t, 2)
out["timing_seconds"] = TIMING
out["skipped_parts_by_time_guard"] = SKIPPED
out["elapsed_seconds"] = round(time.time() - T_START, 2)
out["scope"] = ("Generated window-end routine executed on a simulator of the scalar AArch64 subset it uses and compared "
                "with big integers; tag certificate checked in this layout; toy-scale mean collision times; ledger "
                "arithmetic; base package checks. No production DLP is run.")
print(json.dumps(out, indent=1, default=str))
