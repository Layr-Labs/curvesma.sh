"""Frobenius-class rho on BN254 GT: kernel listings, counted tower, identities, walks and ledger.

Part K  the straight-line AArch64 listings of machine.md (Montgomery multiplication, modular
        addition and subtraction for 4-limb moduli) executed by a small instruction-level
        simulator on edge and random inputs for both p and r; instruction counts are the
        lengths of the executed listings.
Part A  exact identities on the protected BN254 GT: g from the profile has order r, the
        profile basis maps to the priced tower, Frobenius formula, 4-multiplication key
        formula, canonical conjugate, trace invariant, Granger-Scott squaring, operation
        counts of the priced Fp12 multiplication and squaring, orders of p and 2 modulo r.
Part C  the production step rule executed in the real BN254 GT for a few thousand steps.
Part B  the complete algorithm in toy fields with the same tower (u^2=-1, w^6=9+u).
Part E  the pair-collision law of the mixed step in toy fields.
Part D  exact recomputation of the instruction ledger and the score.
Toy runs are evidence at their own scale only. Standard library only.
"""
import json
import math
import random
import sys
import time
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction
from pathlib import Path

T_START = time.time()
cfg = json.load(sys.stdin)
SEED = cfg.get("seed", 1)
KERNEL_CASES = cfg.get("kernel_cases", 300)
PROD_STEPS = cfg.get("production_steps", 1500)
Q_TOY = cfg.get("q", 151)
TRIALS = cfg.get("trials", 60)
R_TOY = cfg.get("R", 64)
SQ8 = cfg.get("sq8", 4)
DP_BITS = cfg.get("dp_bits", 3)
WALKS = cfg.get("walks", 4)
BUDGET = cfg.get("seconds", 6)
PAIR_SAMPLES = cfg.get("pair_samples", 0)
SHOW_LISTINGS = cfg.get("listings", False)
rng = random.Random(SEED)

base = Path("/snapshot/files")
if not (base / "target-profile.json").is_file():
    base = Path(__file__).resolve().parent
profile = json.loads((base / "target-profile.json").read_text())
p = int(profile["field_modulus"], 16)
r = int(profile["subgroup_order"], 16)
G_PROFILE = [int(c, 16) for c in profile["generator_g"]["coefficients"]]
M64 = (1 << 64) - 1
RMONT = 1 << 256
R_FACT = {2: 28, 3: 2, 13: 1, 29: 1, 983: 1, 11003: 1, 237073: 1, 405928799: 1,
          1670836401704629: 1, 13818364434197438864469338081: 1}

# =============================================================== Part K: AArch64 listings
# Register use: x0 = result pointer, x1 = a pointer, x2 = b pointer, x3..x17 scratch.
# Only caller-saved registers x0..x17 are written, so no spills are needed. xzr reads 0.
# Constants block at label K: p0, p1, p2, p3, n0 = -m^-1 mod 2^64 (one block per modulus).


def listing_fp_mul():
    """Montgomery product a*b*2^-256 mod m for a, b < m < 2^254: product scanning by rows,
    then four word-level REDC rounds, then one conditional subtraction."""
    L = []
    b = ["x3", "x4", "x5", "x6"]
    t = ["x7", "x8", "x9", "x10", "x11", "x12", "x13", "x14"]
    T = ["x15", "x16", "x17", "x2"]
    L += ["ldp x3, x4, [x2, #0]", "ldp x5, x6, [x2, #16]"]
    # row 0, a0 held in x14 (t7 is not live before row 3)
    L += ["ldr x14, [x1, #0]"]
    L += ["mul %s, x14, %s" % (t[j], b[j]) for j in range(4)]
    L += ["umulh %s, x14, %s" % (T[j], b[j]) for j in range(4)]
    L += ["adds x8, x8, x15", "adcs x9, x9, x16", "adcs x10, x10, x17", "adc x11, x2, xzr"]
    for i in (1, 2, 3):
        ai = "x14" if i < 3 else "x1"
        L += ["ldr %s, [x1, #%d]" % (ai, 8 * i)]
        L += ["mul %s, %s, %s" % (T[j], ai, b[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i], t[i], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["adc %s, xzr, xzr" % t[i + 4]]
        L += ["umulh %s, %s, %s" % (T[j], ai, b[j]) for j in range(4)]
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2)]
        L += ["adc %s, %s, %s" % (t[i + 4], t[i + 4], T[3])]
    # constants: p0..p3 -> x3..x6, n0 -> x1
    L += ["adr x2, K", "ldp x3, x4, [x2, #0]", "ldp x5, x6, [x2, #16]", "ldr x1, [x2, #32]"]
    pr = ["x3", "x4", "x5", "x6"]
    for i in range(4):
        m = t[i]
        L += ["subs xzr, %s, #1" % t[i]]              # C = (t_i != 0) = carry of t_i + lo(m*p0)
        L += ["mul %s, %s, x1" % (m, t[i])]           # m = t_i * n0 (flags untouched)
        L += ["mul %s, %s, %s" % (T[j], m, pr[j]) for j in (1, 2, 3)]
        L += ["adcs %s, %s, %s" % (t[i + j], t[i + j], T[j]) for j in (1, 2, 3)]
        L += ["umulh %s, %s, %s" % (T[j], m, pr[j]) for j in range(4)]
        L += ["adc x2, x2, xzr"]                      # hi3 + carry <= 2^64-1
        L += ["adds %s, %s, %s" % (t[i + 1], t[i + 1], T[0])]
        L += ["adcs %s, %s, %s" % (t[i + 1 + j], t[i + 1 + j], T[j]) for j in (1, 2, 3)]
        for k in range(i + 5, 8):
            L += ["adcs %s, %s, xzr" % (t[k], t[k])]
    U = t[4:]
    L += ["subs x15, %s, x3" % U[0], "sbcs x16, %s, x4" % U[1], "sbcs x17, %s, x5" % U[2],
          "sbcs x2, %s, x6" % U[3]]
    L += ["csel %s, %s, %s, cs" % (U[j], T[j], U[j]) for j in range(4)]
    L += ["stp %s, %s, [x0, #0]" % (U[0], U[1]), "stp %s, %s, [x0, #16]" % (U[2], U[3]), "ret"]
    return L


def listing_fp_add():
    """(a+b) mod m for a, b < m < 2^254 (a+b < 2^255, no carry out of 256 bits)."""
    return ["ldp x3, x4, [x1, #0]", "ldp x5, x6, [x1, #16]",
            "ldp x7, x8, [x2, #0]", "ldp x9, x10, [x2, #16]",
            "adds x3, x3, x7", "adcs x4, x4, x8", "adcs x5, x5, x9", "adc x6, x6, x10",
            "adr x11, K", "ldp x7, x8, [x11, #0]", "ldp x9, x10, [x11, #16]",
            "subs x7, x3, x7", "sbcs x8, x4, x8", "sbcs x9, x5, x9", "sbcs x10, x6, x10",
            "csel x3, x7, x3, cs", "csel x4, x8, x4, cs", "csel x5, x9, x5, cs", "csel x6, x10, x6, cs",
            "stp x3, x4, [x0, #0]", "stp x5, x6, [x0, #16]", "ret"]


def listing_fp_sub():
    """(a-b) mod m for a, b < m: subtract, then add m masked by the borrow."""
    return ["ldp x3, x4, [x1, #0]", "ldp x5, x6, [x1, #16]",
            "ldp x7, x8, [x2, #0]", "ldp x9, x10, [x2, #16]",
            "subs x3, x3, x7", "sbcs x4, x4, x8", "sbcs x5, x5, x9", "sbcs x6, x6, x10",
            "csetm x11, cc",
            "adr x12, K", "ldp x7, x8, [x12, #0]", "ldp x9, x10, [x12, #16]",
            "and x7, x7, x11", "and x8, x8, x11", "and x9, x9, x11", "and x10, x10, x11",
            "adds x3, x3, x7", "adcs x4, x4, x8", "adcs x5, x5, x9", "adc x6, x6, x10",
            "stp x3, x4, [x0, #0]", "stp x5, x6, [x0, #16]", "ret"]


ALLOWED = {"x%d" % i for i in range(18)} | {"xzr"}


def simulate(listing, mem, regs, kaddr):
    """Execute a straight-line listing. Returns the executed instruction count."""
    C = 0
    count = 0

    def rd(x):
        if x not in ALLOWED:
            raise RuntimeError("register outside x0..x17: " + x)
        return 0 if x == "xzr" else regs[x]

    def wr(x, v):
        if x not in ALLOWED:
            raise RuntimeError("register outside x0..x17: " + x)
        if x != "xzr":
            regs[x] = v & M64

    def addr(s):
        s = s.strip("[]")
        bn, off = [u.strip() for u in s.split(",")]
        return rd(bn) + int(off.lstrip("#"))

    for line in listing:
        count += 1
        op, _, rest = line.partition(" ")
        args = [a.strip() for a in rest.replace(", [", "|[").split("|")]
        ops = [a.strip() for a in args[0].split(",")] if args[0] else []
        if op == "ret":
            break
        if op == "ldr":
            wr(ops[0], mem[addr(args[1])])
        elif op == "ldp":
            a = addr(args[1])
            wr(ops[0], mem[a]); wr(ops[1], mem[a + 8])
        elif op == "stp":
            a = addr(args[1])
            mem[a] = rd(ops[0]); mem[a + 8] = rd(ops[1])
        elif op == "adr":
            wr(ops[0], kaddr)
        elif op == "mul":
            wr(ops[0], rd(ops[1]) * rd(ops[2]))
        elif op == "umulh":
            wr(ops[0], (rd(ops[1]) * rd(ops[2])) >> 64)
        elif op in ("adds", "adcs", "adc", "add"):
            cin = C if op in ("adcs", "adc") else 0
            s = rd(ops[1]) + rd(ops[2]) + cin
            wr(ops[0], s)
            if op in ("adds", "adcs"):
                C = s >> 64
        elif op in ("subs", "sbcs"):
            bv = int(ops[2].lstrip("#")) if ops[2].startswith("#") else rd(ops[2])
            bin_ = (1 - C) if op == "sbcs" else 0
            d = rd(ops[1]) - bv - bin_
            wr(ops[0], d)
            C = 1 if d >= 0 else 0
        elif op == "csel":
            wr(ops[0], rd(ops[1]) if (C == 1) == (ops[3] == "cs") else rd(ops[2]))
        elif op == "csetm":
            wr(ops[0], M64 if (C == 1) == (ops[1] == "cs") else 0)
        elif op == "and":
            wr(ops[0], rd(ops[1]) & rd(ops[2]))
        else:
            raise RuntimeError("unsupported " + line)
    return count


def limbs(v):
    return [(v >> (64 * i)) & M64 for i in range(4)]


def run_kernel(listing, m, a, b):
    n0 = (-pow(m, -1, 1 << 64)) % (1 << 64)
    mem = {}
    A, B, O, K = 0x1000, 0x2000, 0x3000, 0x4000
    for i, w in enumerate(limbs(a)):
        mem[A + 8 * i] = w
    for i, w in enumerate(limbs(b)):
        mem[B + 8 * i] = w
    for i, w in enumerate(limbs(m) + [n0]):
        mem[K + 8 * i] = w
    regs = {"x%d" % i: rng.getrandbits(64) for i in range(18)}
    regs.update({"x0": O, "x1": A, "x2": B})
    n = simulate(listing, mem, regs, K)
    return sum(mem[O + 8 * i] << (64 * i) for i in range(4)), n


def part_K(cases):
    out = {}
    Ls = {"fp_mul": listing_fp_mul(), "fp_add": listing_fp_add(), "fp_sub": listing_fp_sub()}
    for name, mod in (("p", p), ("r", r)):
        Rinv = pow(RMONT, -1, mod)
        edge = [0, 1, 2, mod - 2, mod - 1, (1 << 64) - 1, 1 << 64, mod >> 1, (mod >> 1) + 1,
                (1 << 192) - 1, 1 << 192, mod - (1 << 64)]
        pairs = [(x, y) for x in edge for y in edge]
        pairs += [(rng.randrange(mod), rng.randrange(mod)) for _ in range(cases)]
        # inputs whose limbs are all large exercise every carry chain
        pairs += [(mod - 1 - rng.getrandbits(32), mod - 1 - rng.getrandbits(32)) for _ in range(cases // 4)]
        fails = 0
        counts = set()
        for x, y in pairs:
            want = {"fp_mul": x * y * Rinv % mod, "fp_add": (x + y) % mod, "fp_sub": (x - y) % mod}
            for k, L in Ls.items():
                got, n = run_kernel(L, mod, x, y)
                counts.add((k, n))
                fails += got != want[k]
        out["modulus_" + name] = {"cases": len(pairs), "failures": fails}
    out["executed_instructions"] = {k: n for k, n in sorted(counts)}
    out["listing_lengths"] = {k: len(L) for k, L in Ls.items()}
    if SHOW_LISTINGS:
        out["listings"] = Ls
    return out, out["listing_lengths"]


# =============================================================== counted tower
class Tower:
    """F_q2 = F_q[u]/(u^2+1), F_q12 = F_q2[w]/(w^6 - xi), xi = 9+u.
    Grouping for the priced multiplication: F_q6 = F_q2[v]/(v^3 - xi), v = w^2,
    F_q12 = F_q6[w]/(w^2 - v): y = (c0 + c2 v + c4 v^2) + (c1 + c3 v + c5 v^2) w.
    With count=dict, every base-field operation and every non-leaf routine call is counted."""

    def __init__(self, q):
        self.q = q
        self.count = None
        self.xi = (9 % q, 1)
        self.d = [self.f2pow(self.xi, m * (q - 1) // 6) for m in range(6)]
        z = self.f2pow(self.xi, (q * q - 1) // 6)
        assert z[1] == 0 and (z[0] * z[0] - z[0] + 1) % q == 0
        self.zeta = z[0]
        self.one = ((1, 0),) + ((0, 0),) * 5
        self._fk = {}

    # ---- base field, counted
    def c(self, k, n=1):
        if self.count is not None:
            self.count[k] = self.count.get(k, 0) + n

    def fadd(self, a, b):
        self.c("fp_add")
        return (a + b) % self.q

    def fsub(self, a, b):
        self.c("fp_sub")
        return (a - b) % self.q

    def fmul(self, a, b):
        self.c("fp_mul")
        return a * b % self.q

    # ---- Fp2, counted (add/sub/xi are inline pairs of leaf calls; mul/sqr are routines)
    def f2add(self, a, b):
        return (self.fadd(a[0], b[0]), self.fadd(a[1], b[1]))

    def f2sub(self, a, b):
        return (self.fsub(a[0], b[0]), self.fsub(a[1], b[1]))

    def f2mulxi(self, a):
        # (a0 + a1 u)(9 + u) = (9a0 - a1) + (9a1 + a0) u ; 9t = 8t + t by three doublings and an add
        def nine(t):
            t2 = self.fadd(t, t)
            t4 = self.fadd(t2, t2)
            t8 = self.fadd(t4, t4)
            return self.fadd(t8, t)
        n0, n1 = nine(a[0]), nine(a[1])
        return (self.fsub(n0, a[1]), self.fadd(n1, a[0]))

    def f2mul(self, a, b):
        self.c("routine")
        v0 = self.fmul(a[0], b[0])
        v1 = self.fmul(a[1], b[1])
        s = self.fadd(a[0], a[1])
        t = self.fadd(b[0], b[1])
        v2 = self.fmul(s, t)
        c0 = self.fsub(v0, v1)
        c1 = self.fsub(self.fsub(v2, v0), v1)
        return (c0, c1)

    def f2sqr(self, a):
        self.c("routine")
        s = self.fadd(a[0], a[1])
        d = self.fsub(a[0], a[1])
        c0 = self.fmul(s, d)
        m = self.fmul(a[0], a[1])
        return (c0, self.fadd(m, m))

    def f2conj(self, a):
        return (a[0], self.fsub(0, a[1]))

    # ---- uncounted helpers
    def f2mul_raw(self, a, b):
        q = self.q
        return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

    def f2pow(self, a, e):
        res = (1, 0)
        while e:
            if e & 1:
                res = self.f2mul_raw(res, a)
            a = self.f2mul_raw(a, a)
            e >>= 1
        return res

    # ---- Fp6 / Fp12 multiplication (Karatsuba), counted
    def f6mul(self, a, b):
        self.c("routine")
        a0, a1, a2 = a
        b0, b1, b2 = b
        v0, v1, v2 = self.f2mul(a0, b0), self.f2mul(a1, b1), self.f2mul(a2, b2)
        t = self.f2mul(self.f2add(a1, a2), self.f2add(b1, b2))
        c0 = self.f2add(v0, self.f2mulxi(self.f2sub(self.f2sub(t, v1), v2)))
        t = self.f2mul(self.f2add(a0, a1), self.f2add(b0, b1))
        c1 = self.f2add(self.f2sub(self.f2sub(t, v0), v1), self.f2mulxi(v2))
        t = self.f2mul(self.f2add(a0, a2), self.f2add(b0, b2))
        c2 = self.f2sub(self.f2add(self.f2sub(t, v0), v1), v2)
        return (c0, c1, c2)

    def f6add(self, a, b):
        return tuple(self.f2add(x, y) for x, y in zip(a, b))

    def f6sub(self, a, b):
        return tuple(self.f2sub(x, y) for x, y in zip(a, b))

    def mul_counted(self, x, y):
        self.c("routine")
        X0, X1 = (x[0], x[2], x[4]), (x[1], x[3], x[5])
        Y0, Y1 = (y[0], y[2], y[4]), (y[1], y[3], y[5])
        t0, t1 = self.f6mul(X0, Y0), self.f6mul(X1, Y1)
        vt1 = (self.f2mulxi(t1[2]), t1[0], t1[1])
        C0 = self.f6add(t0, vt1)
        t2 = self.f6mul(self.f6add(X0, X1), self.f6add(Y0, Y1))
        C1 = self.f6sub(self.f6sub(t2, t0), t1)
        return (C0[0], C1[0], C0[1], C1[1], C0[2], C1[2])

    def mul(self, x, y):
        """uncounted schoolbook product, used for speed in setup and toy runs"""
        q = self.q
        acc = [[0, 0] for _ in range(11)]
        for i in range(6):
            a0, a1 = x[i]
            for j in range(6):
                b0, b1 = y[j]
                t = acc[i + j]
                t[0] += a0 * b0 - a1 * b1
                t[1] += a0 * b1 + a1 * b0
        out = []
        for m in range(6):
            c0, c1 = acc[m]
            if m + 6 <= 10:
                h0, h1 = acc[m + 6]            # w^6 = 9 + u
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

    def frobk(self, y, k):
        k %= 12
        if k == 0:
            return y
        D = self._fk.get(k)
        if D is None:
            D = [self.f2pow(self.xi, m * (self.q ** k - 1) // 6) for m in range(6)]
            self._fk[k] = D
        q = self.q
        if k & 1:
            return tuple(self.f2mul_raw((c[0], (-c[1]) % q), D[m]) for m, c in enumerate(y))
        return tuple(self.f2mul_raw(c, D[m]) for m, c in enumerate(y))

    # ---- keys and cyclotomic squaring, counted
    def keys(self, y):
        """Re(w^1 coefficient of y^(q^k)), k = 0..11: 4 multiplications, 1 addition, 8 sub/neg."""
        a, b1 = y[1]
        B = self.fadd(self.fmul(a, self.d[1][0]), self.fmul(b1, self.d[1][1]))
        za, zb = self.fmul(a, self.zeta), self.fmul(B, self.zeta)
        A6 = [a, za, self.fsub(za, a), self.fsub(0, a), self.fsub(0, za), self.fsub(a, za)]
        B6 = [B, zb, self.fsub(zb, B), self.fsub(0, B), self.fsub(0, zb), self.fsub(B, zb)]
        out = [0] * 12
        for j in range(6):
            out[2 * j], out[2 * j + 1] = A6[j], B6[j]
        return out

    def keys_direct(self, y):
        return [self.frobk(y, k)[1][0] for k in range(12)]

    def fp4sqr(self, a, b):
        self.c("routine")
        a2, b2 = self.f2sqr(a), self.f2sqr(b)
        ab2 = self.f2sub(self.f2sub(self.f2sqr(self.f2add(a, b)), a2), b2)
        return self.f2add(a2, self.f2mulxi(b2)), ab2

    def cyc_sqr(self, y):
        """Granger-Scott squaring on the cyclotomic subgroup. A=(c0,c3), B=(c1,c4), C=(c2,c5)
        in F_q4 = F_q2[s]/(s^2 - xi), s = w^3:
        c0'=3A0-2c0, c3'=3A1+2c3, c1'=3 xi C1+2c1, c4'=3C0-2c4, c2'=3B0-2c2, c5'=3B1+2c5."""
        self.c("routine")
        c0, c1, c2, c3, c4, c5 = y
        A0, A1 = self.fp4sqr(c0, c3)
        B0, B1 = self.fp4sqr(c1, c4)
        C0, C1 = self.fp4sqr(c2, c5)
        xC1 = self.f2mulxi(C1)

        def minus(t, u):        # 3t - 2u = t + 2(t - u)
            dd = self.f2sub(t, u)
            return self.f2add(t, self.f2add(dd, dd))

        def plus(t, u):         # 3t + 2u = t + 2(t + u)
            ss = self.f2add(t, u)
            return self.f2add(t, self.f2add(ss, ss))

        return (minus(A0, c0), plus(xC1, c1), minus(B0, c2), plus(A1, c3), minus(C0, c4), plus(B1, c5))

    def rand_gt(self, rng_, rp):
        hc = (self.q ** 12 - 1) // rp
        while True:
            x = tuple((rng_.randrange(self.q), rng_.randrange(self.q)) for _ in range(6))
            y = self.pow(x, hc)
            if y != self.one:
                return y


def profile_to_tower(c, q):
    """sum c_i w^i with w^12 = 18 w^6 - 82, u = w^6 - 9: tower coeff m = (c_m + 9 c_{m+6}) + c_{m+6} u."""
    return tuple(((c[m] + 9 * c[m + 6]) % q, c[m + 6] % q) for m in range(6))


def poly_mul_profile(a, b, q):
    acc = [0] * 23
    for i in range(12):
        for j in range(12):
            acc[i + j] += a[i] * b[j]
    for k in range(22, 11, -1):
        h = acc[k]
        acc[k - 6] += 18 * h
        acc[k - 12] -= 82 * h
        acc[k] = 0
    return [v % q for v in acc[:12]]


def mult_order(a, n, fact):
    o = n - 1
    for pf, e in fact.items():
        for _ in range(e):
            if pow(a, o // pf, n) == 1:
                o //= pf
            else:
                break
    return o


# =============================================================== Part A
def part_A(checks):
    P = Tower(p)
    out = {}
    prod = 1
    for k, e in R_FACT.items():
        prod *= k ** e
    o2 = mult_order(2, r, R_FACT)
    op_ = mult_order(p % r, r, R_FACT)
    out.update({
        "p_mod_12": p % 12,
        "r_minus_1_factorisation_verified": prod == r - 1,
        "order_of_p_mod_r": op_,
        "p6_is_minus1_mod_r": pow(p, 6, r) == r - 1,
        "phi12_p_is_0_mod_r": (pow(p, 4, r) - pow(p, 2, r) + 1) % r == 0,
        "log2_order_of_2_mod_r": round(math.log2(o2), 3),
        "smallest_t_with_2^t_in_<p>_log2": round(math.log2(o2 // math.gcd(o2, 12)), 3),
        "zeta_in_Fp_primitive_6th_root": (P.zeta * P.zeta - P.zeta + 1) % p == 0,
        "xi_not_square_not_cube_in_Fp2": P.f2pow(P.xi, (p * p - 1) // 2) != (1, 0) and P.f2pow(P.xi, (p * p - 1) // 3) != (1, 0),
    })
    # profile basis -> tower: ring homomorphism check on random elements, and g of order r
    a = [rng.randrange(p) for _ in range(12)]
    b = [rng.randrange(p) for _ in range(12)]
    out["profile_to_tower_is_ring_map"] = profile_to_tower(poly_mul_profile(a, b, p), p) == P.mul(profile_to_tower(a, p), profile_to_tower(b, p))
    g = profile_to_tower(G_PROFILE, p)
    out["g_order_r"] = P.pow(g, r) == P.one and g != P.one
    ok = dict.fromkeys(["frob", "keys", "distinct", "canon", "inv", "sq", "mul", "trace_inv", "stored_canon"], True)
    counts_mul, counts_sq, counts_keys = None, None, None
    Rp = RMONT % p
    for _ in range(checks):
        y = P.pow(g, rng.randrange(1, r))
        x2 = P.pow(g, rng.randrange(1, r))
        ok["frob"] &= P.frobk(y, 1) == P.pow(y, p)
        kd, kf = P.keys_direct(y), P.keys(y)
        ok["keys"] &= kd == kf
        ok["distinct"] &= len(set(kf)) == 12
        conj = [P.frobk(y, k) for k in range(12)]
        canon, canon_st = set(), set()
        for cj in conj:
            kk = P.keys(cj)
            canon.add(P.frobk(cj, kk.index(min(kk))))
            ks = [v * Rp % p for v in kk]          # stored Montgomery limbs
            canon_st.add(P.frobk(cj, ks.index(min(ks))))
        ok["canon"] &= len(canon) == 1
        ok["stored_canon"] &= len(canon_st) == 1
        ok["inv"] &= P.mul(conj[6], y) == P.one
        ok["trace_inv"] &= len({cj[0][0] for cj in conj}) == 1
        P.count = {}
        s = P.cyc_sqr(y)
        counts_sq = dict(P.count)
        P.count = {}
        m = P.mul_counted(y, x2)
        counts_mul = dict(P.count)
        P.count = {}
        P.keys(y)
        counts_keys = dict(P.count)
        P.count = None
        ok["sq"] &= s == P.mul(y, y)
        ok["mul"] &= m == P.mul(y, x2)
    out.update({"random_GT_elements": checks, **{k + "_ok": v for k, v in ok.items()}})
    out["op_counts_fp12_mul"] = counts_mul
    out["op_counts_cyclotomic_square"] = counts_sq
    out["op_counts_keys"] = counts_keys
    return out, counts_mul, counts_sq, counts_keys


# =============================================================== Part C: production walk
def part_C(steps, rprod=16, verify=150):
    """Production step rule in the real BN254 GT; bits read from stored T = Re(c_0)*2^256 mod p."""
    P = Tower(p)
    Rp = RMONT % p
    lam = [pow(p, k, r) for k in range(12)]
    g = profile_to_tower(G_PROFILE, p)
    x0 = rng.randrange(1, r)
    z0 = P.pow(g, x0)
    ah = [rng.randrange(r) for _ in range(rprod)]
    bh = [rng.randrange(r) for _ in range(rprod)]
    Mh = [P.mul(P.pow(g, ah[i]), P.pow(z0, bh[i])) for i in range(rprod)]
    Tab = [[P.frobk(Mh[i], k) for k in range(12)] for i in range(rprod)]
    a0, b0 = rng.randrange(r), rng.randrange(r)
    y = P.mul(P.pow(g, a0), P.pow(z0, b0))
    al, be = a0, b0
    n_sq = n_mul = bad = 0
    bins_sq, bins_h = [0] * 8, [0] * 16
    two = 0
    hist = []
    t0 = time.time()
    done = 0
    for t in range(steps):
        if time.time() - t0 > 4:
            break
        done += 1
        u = y[0][0] * Rp % p
        if len(hist) >= 2 and hist[-2][0] == u:
            two += 1
            _, yy, aa, bb = min(hist[-1], (u, y, al, be), key=lambda s_: s_[0])
            y, al, be = P.cyc_sqr(yy), 2 * aa % r, 2 * bb % r
            hist = []
            continue
        hist.append((u, y, al, be))
        if len(hist) > 2:
            hist.pop(0)
        bins_sq[u & 7] += 1
        if (u & 7) < SQ8:
            y2 = P.cyc_sqr(y)
            if t < verify:
                kk = [v * Rp % p for v in P.keys(y)]
                kap = kk.index(min(kk))
                bad += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), P.frobk(y, kap))
            y, al, be = y2, 2 * al % r, 2 * be % r
            n_sq += 1
        else:
            kk = [v * Rp % p for v in P.keys(y)]
            kap = kk.index(min(kk))
            h = (u >> 3) % rprod
            bins_h[(u >> 3) & 15] += 1
            idx = (12 - kap) % 12
            y2 = P.mul_counted(y, Tab[h][idx]) if t < 20 else P.mul(y, Tab[h][idx])
            if t < verify:
                bad += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), Mh[h])
            y = y2
            al = (al + ah[h] * lam[idx]) % r
            be = (be + bh[h] * lam[idx]) % r
            n_mul += 1
    exp_ok = y == P.mul(P.pow(g, al), P.pow(z0, be))
    n_all = sum(bins_sq)
    n_br = sum(bins_h)
    chi_sq = sum((c_ - n_all / 8) ** 2 / (n_all / 8) for c_ in bins_sq) if n_all else None
    chi_h = sum((c_ - n_br / 16) ** 2 / (n_br / 16) for c_ in bins_h) if n_br else None
    return {"steps_run": done, "R_table": rprod, "squaring_steps": n_sq, "multiply_steps": n_mul,
            "squaring_fraction": round(n_sq / max(1, n_all), 4), "target_fraction": SQ8 / 8,
            "explicit_class_check_failures_in_first_%d" % verify: int(bad),
            "exponent_invariant_y_equals_g^al_z^be": exp_ok,
            "chi2_low3_bits_of_stored_T_7dof": round(chi_sq, 2) if chi_sq is not None else None,
            "chi2_next4_bits_of_stored_T_15dof": round(chi_h, 2) if chi_h is not None else None,
            "two_cycle_events": two}


# =============================================================== Part B: toy algorithm
def toy_prime(q):
    n = q ** 4 - q * q + 1
    d, best = 2, 1
    while d * d <= n:
        while n % d == 0:
            best = max(best, d)
            n //= d
        d += 1
    return max(best, n)


def mix64(x):
    x &= M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)


def classinfo(T, y):
    kf = T.keys(y)
    m = min(kf)
    ks = [k for k in range(12) if kf[k] == m]
    if len(ks) == 1:
        return ks[0], T.frobk(y, ks[0]), m
    can, kap = min((T.frobk(y, k), k) for k in ks)
    return kap, can, m


def toy_hash(can):
    h = 0x9E3779B97F4A7C15
    for c_ in can:
        h = mix64(h ^ (c_[0] * 0x100000001B3 + c_[1]))
    return h


def run_toy(T, q, rp, walks, Rm, sq8, dp_bits, cap_factor=40):
    qk = [pow(q, k, rp) for k in range(12)]
    g = T.rand_gt(rng, rp)
    xs = rng.randrange(1, rp)
    z0 = T.pow(g, xs)
    ah = [rng.randrange(rp) for _ in range(Rm)]
    bh = [rng.randrange(rp) for _ in range(Rm)]
    Mh = [T.mul(T.pow(g, ah[i]), T.pow(z0, bh[i])) for i in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    At = [[ah[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    Bt = [[bh[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    mask = (1 << dp_bits) - 1

    def new_walk():
        a, b = rng.randrange(rp), rng.randrange(rp)
        return {"y": T.mul(T.pow(g, a), T.pow(z0, b)), "a": a, "b": b, "hist": [], "ck": None, "n": 0, "since": 0}

    def advance(y, a, b):
        kap, can, K = classinfo(T, y)
        u = toy_hash(can)
        if (u & 7) < sq8:
            return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp
        h = (u >> 3) % Rm
        idx = (12 - kap) % 12
        return T.mul(y, Tab[h][idx]), (a + At[h][idx]) % rp, (b + Bt[h][idx]) % rp

    W = [new_walk() for _ in range(walks)]
    visited, dps = {}, {}
    steps, Tc, found = 0, None, None
    st = {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0, "sq": 0, "ord": 0}
    while found is None:
        for w in W:
            y, a, b = w["y"], w["a"], w["b"]
            kap, can, K = classinfo(T, y)
            u = toy_hash(can)
            ca, cb = a * qk[kap] % rp, b * qk[kap] % rp
            if Tc is None:
                old = visited.get(can)
                if old is None:
                    visited[can] = (ca, cb)
                elif old != (ca, cb):
                    Tc = steps
            if (u >> 13) & mask == 0:
                old = dps.get(can)
                if old is not None and old[1] != cb:
                    xr = (old[0] - ca) * pow(cb - old[1], -1, rp) % rp
                    if T.pow(g, xr) != z0:
                        raise RuntimeError("recovered logarithm failed verification")
                    found = xr
                    break
                if old is not None:
                    st["useless"] += 1
                dps[can] = (ca, cb)
                w["since"] = 0
            hist = w["hist"]
            key = (K, can)
            if len(hist) >= 2 and hist[-2][0][1] == can:
                st["two_cycles"] += 1
                _, yy, aa, bb = min(hist[-1], (key, y, a, b), key=lambda s_: s_[0])
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["ck"] is not None and w["ck"][0] == can and w["ck"][1] < w["n"]:
                st["long_cycles"] += 1
                best = (key, y, a, b)
                yy, aa, bb = y, a, b
                for _ in range(w["n"] - w["ck"][1]):
                    yy, aa, bb = advance(yy, aa, bb)
                    steps += 1
                    k2, c2, K2 = classinfo(T, yy)
                    if (K2, c2) < best[0]:
                        best = ((K2, c2), yy, aa, bb)
                _, yy, aa, bb = best
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["n"] % 64 == 0:
                w["ck"] = (can, w["n"])
            hist.append((key, y, a, b))
            if len(hist) > 2:
                hist.pop(0)
            st["ord"] += 1
            st["sq"] += (u & 7) < sq8
            w["y"], w["a"], w["b"] = advance(y, a, b)
            w["n"] += 1
            w["since"] += 1
            steps += 1
            if w["since"] > cap_factor << dp_bits:
                st["restarts"] += 1
                w.update(new_walk())
    return Tc, steps, st


def part_B(q, trials, Rm, sq8, dp_bits, walks, budget):
    T = Tower(q)
    rp = toy_prime(q)
    if not (pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))):
        raise RuntimeError("toy field does not have embedding degree 12")
    N = (rp - 1) // 12
    s = sq8 / 8
    sumq2 = s * s + (1 - s) ** 2 / Rm
    Np = N / (1 - sumq2)
    Tcs, ev = [], {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0}
    t0 = time.time()
    allsteps = 0
    while len(Tcs) < trials and time.time() - t0 < budget:
        Tc, steps, st = run_toy(T, q, rp, walks, Rm, sq8, dp_bits)
        Tcs.append(Tc)
        allsteps += steps
        for k in ev:
            ev[k] += st[k]
    n = len(Tcs)
    sN = math.sqrt(N)
    S = {x: math.sqrt(2 * Np * math.log(1 / (1 - x))) for x in (0.40, 0.42, 0.43)}
    return {"q": q, "r_prime": rp, "classes_N": N, "R": Rm, "squaring_prob": s,
            "effective_N_over_N": round(1 / (1 - sumq2), 4), "trials_completed": n,
            "logs_recovered_and_verified": n,
            "mean_Tc_over_sqrtN": round(sum(Tcs) / n / sN, 4) if n else None,
            "model_mean_over_sqrtN": round(math.sqrt(math.pi * Np / 2) / sN, 4),
            **{"fraction_Tc_le_S_model_%.2f" % x: round(sum(t <= S[x] for t in Tcs) / max(1, n), 3) for x in S},
            "events": ev, "two_cycles_per_step": round(ev["two_cycles"] / max(1, allsteps), 6),
            "predicted_two_cycles_per_step": round((1 - s) ** 2 / (12 * Rm), 6)}


def part_E(q, m, Rm, sq8):
    T = Tower(q)
    rp = toy_prime(q)
    N = (rp - 1) // 12
    g = T.rand_gt(rng, rp)
    Tab = []
    for _ in range(Rm):
        Mh = T.pow(g, rng.randrange(1, rp))
        Tab.append([T.frobk(Mh, k) for k in range(12)])
    classes = {}
    tries = 0
    while len(classes) < m and tries < 4 * m:
        tries += 1
        y = T.pow(g, rng.randrange(1, rp))
        kap, can, _ = classinfo(T, y)
        if can not in classes:
            classes[can] = (y, kap)
    images = {}
    for can, (y, kap) in classes.items():
        u = toy_hash(can)
        y2 = T.cyc_sqr(y) if (u & 7) < sq8 else T.mul(y, Tab[(u >> 3) % Rm][(12 - kap) % 12])
        _, can2, _ = classinfo(T, y2)
        images[can2] = images.get(can2, 0) + 1
    pairs = sum(n * (n - 1) // 2 for n in images.values())
    mm = len(classes)
    s = sq8 / 8
    return {"q": q, "classes_N": N, "sampled_classes": mm, "R": Rm, "squaring_prob": s,
            "N_times_pair_collision_rate_observed": round(pairs / (mm * (mm - 1) / 2) * N, 4),
            "image_coincidence_pairs": pairs,
            "predicted_exact": round(1 - s * s - (1 - s) ** 2 / (12 * Rm), 4),
            "predicted_conservative_used": round(1 - s * s - (1 - s) ** 2 / Rm, 4)}


# =============================================================== Part D: ledger
CALL_SITE = 4          # BL plus up to three pointer set-ups (ADD/MOV) per leaf call
ROUTINE = 32           # per non-leaf routine call: call site, frame, 4 callee-saved, epilogue, slack
COMMON = 300           # every step (itemised in machine.md section 4: 129)
KEYMIN = 160           # minimum of 12 four-limb keys, kappa, twist index (itemised 138, tie path included)
SQ_PROB = Fraction(1, 2)
R_TABLE = 1024
P_MODEL = Decimal("0.43")
CLAIM = Decimal("124.65")


def ledger(lens, cmul, csq, ckeys):
    getcontext().prec = 60
    ln2 = Decimal(2).ln()
    lg = lambda v: Decimal(v).ln() / ln2
    leaf = {"fp_mul": lens["fp_mul"] + CALL_SITE, "fp_add": lens["fp_add"] + CALL_SITE,
            "fp_sub": lens["fp_sub"] + CALL_SITE, "routine": ROUTINE}

    def price(cnt):
        return sum(leaf[k] * v for k, v in cnt.items())

    MULI, SQI = price(cmul), price(csq)
    KEYI = price(ckeys) + KEYMIN
    mul_step = MULI + KEYI + COMMON
    sq_step = SQI + COMMON
    C = 20000
    s = SQ_PROB
    cbar = (1 - s) * Fraction(mul_step, C) + s * Fraction(sq_step, C)
    cb = Decimal(cbar.numerator) / Decimal(cbar.denominator)
    N = Fraction(r - 1, 12)
    sumq2 = s * s + (1 - s) ** 2 / R_TABLE
    Np = N / (1 - sumq2)
    Npd = Decimal(Np.numerator) / Decimal(Np.denominator)
    S = (2 * Npd * (1 / (1 - P_MODEL)).ln()).sqrt()
    eps = Decimal(5) / Decimal(12 * R_TABLE) + Decimal(250) / Decimal(R_TABLE * R_TABLE)
    E = Decimal(2 ** 20) * 20 * Decimal(2 ** 60)
    NONWALK = Decimal(2) ** Decimal("64.5")
    W = cb * S * (1 + eps) + cb * E + NONWALK
    Wcap = Decimal(2) ** CLAIM
    n = S * (1 + eps)
    slack = (Wcap - cb * E - NONWALK) / n - cb
    gap = Decimal(mul_step - sq_step) / C
    tol_points = slack / gap * 100

    def success_with_drift(points):
        delta = Decimal(points) / 100 * gap / cb
        steps = min(n, (Wcap - cb * E - NONWALK) / (cb * (1 + delta))) / (1 + eps)
        return round(float(1 - (-(steps * steps) / (2 * Npd)).exp()), 4)

    lw = lg(W)
    ceil2 = (lw * 100).to_integral_value(rounding=ROUND_CEILING) / 100
    return {"leaf_prices_incl_call_site": leaf,
            "fp12_mul_instr": MULI, "cyclotomic_square_instr": SQI,
            "keys_instr": KEYI, "common_instr": COMMON,
            "multiply_step_instr": mul_step, "squaring_step_instr": sq_step,
            "multiply_step_units": float(Fraction(mul_step, C)), "squaring_step_units": float(Fraction(sq_step, C)),
            "squaring_probability": str(s), "R_table": R_TABLE,
            "cbar_units_per_step": float(cbar), "log2_N": float(lg(Decimal(N.numerator) / N.denominator)),
            "sum_q2": float(sumq2), "log2_N_eff": float(lg(Npd)), "log2_S": float(lg(S)), "eps": float(eps),
            "log2_E": float(lg(E)), "log2_nonwalk_bucket": 64.5,
            "log2_W_nominal": float(lw), "least_hundredth_at_or_above": str(ceil2),
            "claim": str(CLAIM), "claim_minus_nominal_bits": float(CLAIM - lw),
            "claim_ok": CLAIM >= lw,
            "cap_slack_units_per_step": float(slack),
            "tolerated_drop_of_squaring_fraction_points": float(tol_points),
            "success_if_squaring_fraction_drops_by_points": {k: success_with_drift(k) for k in ("1", "2", "3", "4", "5")}}


# =============================================================== main
out = {"seed": SEED}
kres, lens = part_K(KERNEL_CASES)
out["partK_kernels"] = kres
ares, cmul, csq, ckeys = part_A(max(1, cfg.get("production_checks", 2)))
out["partA_identities"] = ares
out["partD_ledger"] = ledger(lens, cmul, csq, ckeys)
out["partC_production_walk"] = part_C(max(1, PROD_STEPS))
if PAIR_SAMPLES:
    out["partE_pair_law"] = part_E(Q_TOY, PAIR_SAMPLES, R_TOY, SQ8)
# declared runs (seconds <= 14) stay inside the 20 s runner limit; longer budgets are offline only
budget_B = BUDGET if BUDGET > 14 else max(0.5, min(BUDGET, 15 - (time.time() - T_START)))
out["partB_toy"] = part_B(Q_TOY, TRIALS, R_TOY, SQ8, DP_BITS, WALKS, budget_B)
out["elapsed_seconds"] = round(time.time() - T_START, 2)
out["scope"] = ("Kernel listings are executed instruction by instruction on a simulator of the AArch64 "
                "subset they use; tower identities and the step rule run on the protected BN254 GT; "
                "collision statistics are toy-field evidence only. No production DLP is run.")
print(json.dumps(out, indent=1, default=str))
