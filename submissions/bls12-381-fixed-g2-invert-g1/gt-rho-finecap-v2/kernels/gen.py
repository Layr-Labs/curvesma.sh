"""Emit straight-line scalar AArch64 for two BLS12-381 Fp12 kernels (blst's tower, Montgomery form,
6 x 64-bit limbs per Fp element, outputs fully reduced):

  mul.s  out = x * m, m a fixed table entry stored with its Karatsuba pre-sums (exact integers)
  sqr.s  out = x^2 for x in the cyclotomic subgroup (Granger-Scott, as blst's cyclotomic_sqr_fp12)

Calling convention of both: x0 = out (12 x 48 bytes), x1 = x (12 x 48 bytes), x2 = table entry (mul),
x3 = constant pool (layout.json). General registers only, no branches, no calls; ends in ret.

Lazy reduction: every 384x384 product goes to 768 bits; each output coordinate is one integer
T = K + sum(+-products) (+ 2 a R terms for the squaring) taken mod 2^768, where every sum of
single-width operands is exact (no modular reduction), so T equals a fixed integer polynomial in the
inputs whose range model.py bounds monomial by monomial; K is the least multiple of p that makes T >= 0.
One Montgomery reduction per output coordinate, then a branch-free subtraction of an estimated q p and of p.
"""
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import model as M

P, R = M.P, M.R
N0 = (-pow(P, -1, 1 << 64)) % (1 << 64)
Z = [f"x{i}" for i in range(4, 16)]
A = ["x16", "x17", "x19", "x20", "x21", "x22"]
B = ["x23", "x24"]
T = [f"x{i}" for i in range(25, 31)]
PR, N0R, MR = A, "x23", "x24"
SAVE = 112


POOL = {"qt": 0, "n0": 48, "c": 56, "p": 64, "2p": 128, "4p": 256, "K": 512}
QC = (1 << 64) // ((P >> 320) + 1)
assert QC == 9 and 8 * P < (1 << 384)
assert all(k - 1 <= (k * P >> 320) * QC >> 64 and ((k + 1) * P - 1 >> 320) * QC >> 64 <= k for k in range(8))


class Asm:
    def __init__(self):
        self.lines, self.cats, self.pool, self.addr, self.clock = [], [], [], {}, 0

    def ins(self, cat, text):
        op = text.split()[0]
        if op in ("ldp", "stp", "ldr", "str"):
            cat = "mem"
        self.lines.append("    " + text)
        self.cats.append(cat)

    def comment(self, text):
        for t in text.split("\n"):
            self.lines.append("// " + t)
            self.cats.append(None)

    def set_pool(self, regs, shift=None):
        self.pool, self.addr = list(regs), {}
        self.shift = shift if shift is not None else getattr(self, "shift", 48)

    def keep_pool(self, regs, shift, drop=()):
        """Phase transition that keeps derived bases; `drop` invalidates regs restored between phases."""
        for r in drop:
            self.addr.pop(r, None)
        self.pool = list(regs)
        self.shift = shift

    def pin_sp(self, mid):
        """Pin x18 (unused, caller-saved) to sp+mid for the whole phase; one derive, never evicted."""
        assert mid % 8 == 0
        if mid > 4095:
            self.ins("addr", f"add x18, sp, #{mid >> 12}, lsl #12")
            if mid & 4095:
                self.ins("addr", f"add x18, x18, #{mid & 4095}")
        else:
            self.ins("addr", f"add x18, sp, #{mid}")
        self.addr["x18"] = ("sp", mid, 10 ** 9)

    def pin_base(self, base, mid):
        """Pin x18 to base+mid (x18 unused caller-saved); one derive, never evicted."""
        assert mid % 8 == 0
        if mid > 4095:
            self.ins("addr", f"add x18, {base}, #{mid >> 12}, lsl #12")
            if mid & 4095:
                self.ins("addr", f"add x18, x18, #{mid & 4095}")
        else:
            self.ins("addr", f"add x18, {base}, #{mid}")
        self.addr["x18"] = (base, mid, 10 ** 9)

    def ref(self, base, off, pair=True):
        lo, hi = (-512, 504) if pair else (0, 32760)
        if lo <= off <= hi:
            return f"[{base}, #{off}]"
        self.clock += 1
        for r, (b, o, _) in list(self.addr.items()):
            if b == base and lo <= off - o <= hi:
                if r != "x18":
                    self.addr[r] = (b, o, self.clock)
                return f"[{r}, #{off - o}]"
        free = [r for r in self.pool if r not in self.addr]
        same = [k for k in self.addr if k != "x18" and self.addr[k][0] == base]
        r = free[0] if free else (min(same, key=lambda k: self.addr[k][2]) if same else min((k for k in self.addr if k != "x18"), key=lambda k: self.addr[k][2]))
        new = off + self.shift
        if new > 4095:
            self.ins("addr", f"add {r}, {base}, #{new >> 12}, lsl #12")
            if new & 4095:
                self.ins("addr", f"add {r}, {r}, #{new & 4095}")
        else:
            self.ins("addr", f"add {r}, {base}, #{new}")
        self.addr[r] = (base, new, self.clock)
        return f"[{r}, #{off - new}]"

    def load6(self, regs, base, off):
        for i in range(0, 6, 2):
            self.ins("mem", f"ldp {regs[i]}, {regs[i + 1]}, {self.ref(base, off + 8 * i)}")

    def store6(self, regs, base, off):
        for i in range(0, 6, 2):
            self.ins("mem", f"stp {regs[i]}, {regs[i + 1]}, {self.ref(base, off + 8 * i)}")

    def chain(self, cat, sign, dst, x, y, first=True):
        """dst = x +- y over len(dst) limbs; first=False continues the carry already in C."""
        for i in range(len(dst)):
            if sign > 0:
                op = "adds" if (i == 0 and first) else "adcs"
            else:
                op = "subs" if (i == 0 and first) else "sbcs"
            self.ins(cat, f"{op} {dst[i]}, {x[i]}, {y[i]}")

    def acc12(self, sign, base, off, cat="addsub"):
        """Z (12 limbs) +-= the 768-bit value at base+off, through T in two halves (flags survive ldp)."""
        self.load6(T, base, off)
        self.chain(cat, sign, Z[:6], Z[:6], T)
        self.load6(T, base, off + 48)
        self.chain(cat, sign, Z[6:], Z[6:], T, first=False)

    def text(self):
        return "\n".join(self.lines) + "\n"

    def count(self):
        out = {}
        for c in self.cats:
            if c is not None:
                out[c] = out.get(c, 0) + 1
        return out


def prologue(a, frame):
    left = frame
    while left:
        step = min(left, 4080)
        a.ins("glue", f"sub sp, sp, #{step}")
        left -= step
    regs = ["x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30"]
    for i in range(0, 12, 2):
        a.ins("mem", f"stp {regs[i]}, {regs[i + 1]}, [sp, #{8 * i}]")
    a.ins("mem", "stp x0, x3, [sp, #96]")


def epilogue(a, frame):
    regs = ["x19", "x20", "x21", "x22", "x23", "x24", "x25", "x26", "x27", "x28", "x29", "x30"]
    for i in range(0, 12, 2):
        a.ins("mem", f"ldp {regs[i]}, {regs[i + 1]}, [sp, #{8 * i}]")
    left = frame
    while left:
        step = min(left, 4080)
        a.ins("glue", f"add sp, sp, #{step}")
        left -= step
    a.ins("glue", "ret")


def load_lin(a, dst, terms):
    """dst = sum of sign * vector over `terms` [(sign, base, off)], exact (range checked by the caller)."""
    s0, b0, o0 = terms[0]
    assert s0 > 0
    a.load6(dst, b0, o0)
    for s, b, o in terms[1:]:
        a.load6(T, b, o)
        a.chain("addsub", s, dst, dst, T)


def product(a, aterms, bbase, boff, double_a=False):
    """Z = (sum of aterms) * (vector at bbase+boff), schoolbook by rows of b (143 arithmetic instructions)."""
    load_lin(a, A, aterms)
    if double_a:
        a.chain("addsub", 1, A, A, A)
    for jj in range(0, 6, 2):
        a.ins("mem", f"ldp {B[0]}, {B[1]}, {a.ref(bbase, boff + 8 * jj)}")
        for j, b in ((jj, B[0]), (jj + 1, B[1])):
            if j == 0:
                for i in range(6):
                    a.ins("prod", f"mul {Z[i]}, {A[i]}, {b}")
                for i in range(6):
                    a.ins("prod", f"umulh {T[i]}, {A[i]}, {b}")
                a.chain("prod", 1, Z[1:6], Z[1:6], T[:5])
                a.ins("prod", f"adc {Z[6]}, {T[5]}, xzr")
            else:
                for i in range(6):
                    a.ins("prod", f"mul {T[i]}, {A[i]}, {b}")
                a.chain("prod", 1, Z[j:j + 6], Z[j:j + 6], T)
                a.ins("prod", f"adc {Z[j + 6]}, xzr, xzr")
                for i in range(6):
                    a.ins("prod", f"umulh {T[i]}, {A[i]}, {b}")
                a.chain("prod", 1, Z[j + 1:j + 6], Z[j + 1:j + 6], T[:5])
                a.ins("prod", f"adc {Z[j + 6]}, {Z[j + 6]}, {T[5]}")


def redc_store(a, k):
    """Montgomery-reduce Z (T < 2^768) to r = (T + M p)/2^384 < T/2^384 + p, then subtract the
    estimated multiple q p (r < 8p) and p branch-free, and store coordinate k. Window registers
    rotate; the low limb of each round becomes the next top (its sum with lo(m p0) is 0 mod 2^64,
    carry = [limb != 0])."""
    W = Z[:6]
    for _ in range(6):
        a.ins("redc", f"mul {MR}, {W[0]}, {N0R}")
        for j in range(1, 6):
            a.ins("redc", f"mul {T[j]}, {MR}, {PR[j]}")
        a.ins("redc", f"subs xzr, {W[0]}, #1")
        for j in range(1, 6):
            a.ins("redc", f"adcs {W[j]}, {W[j]}, {T[j]}")
        a.ins("redc", f"adc {W[0]}, xzr, xzr")
        for j in range(6):
            a.ins("redc", f"umulh {T[j]}, {MR}, {PR[j]}")
        a.chain("redc", 1, W[1:6], W[1:6], T[:5])
        a.ins("redc", f"adc {W[0]}, {W[0]}, {T[5]}")
        W = W[1:] + W[:1]
    a.chain("redc", 1, W, W, Z[6:])
    assert W == Z[:6]
    q, ad = Z[6], Z[7]
    a.ins("csub", f"ldr {q}, [x3, #{POOL['c']}]")
    a.ins("csub", f"umulh {q}, {W[5]}, {q}")
    a.ins("csub", f"add {ad}, x3, {q}, lsl #6")
    a.load6(T, ad, POOL["qt"])
    a.chain("csub", -1, W, W, T)
    a.chain("csub", -1, T, W, PR)
    for j in range(6):
        a.ins("csub", f"csel {W[j]}, {W[j]}, {T[j]}, lo")
    a.store6(W, "x0", 48 * k)


def out_bounds(polys, vmax):
    rows = []
    for k, poly in enumerate(polys):
        lo, hi = M.poly_range(poly, vmax)
        K = -(-(-lo) // P) * P if lo < 0 else 0
        top = K + hi
        assert top < (1 << 768)
        bres = top // R + 1 + P
        assert bres <= 8 * P
        rows.append((k, lo, hi, K, bres))
    return rows


def bound_table(rows):
    t = ["output  lo/p^2    hi/p^2    K/p^2   T_max/p^2  r bound/p",
         "(T = K + polynomial in [K + lo, K + hi], lo/hi summed monomial by monomial, inputs in [0, p))"]
    for k, lo, hi, K, bres in rows:
        t.append(f"c{k:<6} {lo / P / P:+8.3f} {hi / P / P:+9.3f} {K / P / P:8.3f} {(K + hi) / P / P:10.3f} "
                 f"{bres / P:10.4f}")
    t.append("(each r < 8p: subtract q p, q = umulh(r5, 9) read from the q*p table, then p unless that borrows)")
    return "\n".join(t)



class AtomBuilder(M.Builder):
    """Records each Fp2 multiplication; its two results are atoms combined at the Fp6/Fp12 levels."""
    def __init__(self):
        super().__init__()
        self.fp2 = []

    def fp2_mul(self, Ax, Bm):
        f = len(self.fp2)
        self.fp2.append((Ax, Bm))
        return ({2 * f: 1}, {2 * f + 1: 1})


def mul_kernel():
    a = Asm()
    bd = AtomBuilder()
    x0, x1 = M.fp6_vars("x", 0), M.fp6_vars("x", 6)
    m0, m1 = M.fp6_vars("m", 0), M.fp6_vars("m", 6)
    t0 = M.fp6_mul(bd, x0, m0)
    t1 = M.fp6_mul(bd, x1, m1)
    t2 = M.fp6_mul(bd, M.fp6_lsum(x0, x1), M.fp6_lsum(m0, m1))
    c1 = [M.f2add(M.f2add(t2[j], t0[j], -1), t1[j], -1) for j in range(3)]
    vt1 = [M.xi(t1[2]), t1[0], t1[1]]
    c0 = [M.f2add(t0[j], vt1[j]) for j in range(3)]
    outs = [c for f6 in (c0, c1) for f2 in f6 for c in f2]
    nf = len(bd.fp2)
    assert nf == 18


    prods = []
    atom_comb = []
    for Ax, Bm in bd.fp2:
        i = len(prods)
        prods += [(Ax[0], Bm[0]), (Ax[1], Bm[1]), (M.ladd(Ax[0], Ax[1]), M.ladd(Bm[0], Bm[1]))]
        atom_comb += [{i: 1, i + 1: -1}, {i + 2: 1, i: -1, i + 1: -1}]
    flat = []
    for c in outs:
        d = {}
        for at, s in c.items():
            for pi, v in atom_comb[at].items():
                d[pi] = d.get(pi, 0) + s * v
        flat.append({k: v for k, v in d.items() if v})
    vmax = {f"x{i}": P - 1 for i in range(12)}
    vmax.update({f"m{i}": P - 1 for i in range(12)})
    polys = [M.expand(prods, c) for c in flat]

    fprods, fouts = M.mul_model()
    assert [M.expand(fprods, c) for c in fouts] == polys
    rows = out_bounds(polys, vmax)


    derived = {}
    need = sorted({tuple(sorted(L.items())) for Ax, _ in bd.fp2 for L in Ax if len(L) > 1}, key=lambda k: (len(k), k))
    scr = SAVE
    for key in need:
        derived[key] = scr
        scr += 48

    def decompose(L):
        L = dict(L)
        terms = []
        for key in sorted(derived, key=lambda k: (len(k), k), reverse=True):
            if all(L.get(n, 0) >= c for n, c in key):
                if True:
                    for n, c in key:
                        L[n] -= c
                    L = {n: c for n, c in L.items() if c}
                    terms.append((1, "sp", derived[key]))
        for n, c in sorted(L.items()):
            assert c == 1
            terms.append((1, "x1", 48 * int(n[1:])))
        return terms

    entry = []
    for Ax, Bm in bd.fp2:
        entry += [Bm[0], Bm[1], M.ladd(Bm[0], Bm[1])]
    tmp0, tmpS = scr, scr + 96
    scr += 192
    atoms = scr
    scr += 96 * 2 * nf
    frame = (scr + 15) // 16 * 16


    rng = lambda L: M.lin_range(L, vmax)[1]
    opmax = max(max(rng(L) for L in list(Ax) + [M.ladd(*Ax)]) for Ax, _ in bd.fp2)
    bmax = max(rng(L) for L in entry)
    assert opmax < (1 << 384) and bmax < (1 << 384)

    a.comment("BLS12-381 Fp12 multiplication by a fixed table entry, Montgomery form, outputs in [0, p)\n"
              "generated by gen.py; args x0 out, x1 x, x2 entry, x3 constant pool\n"
              f"x-side operands are exact sums of reduced inputs, at most {opmax / P:.0f}p - 1 < 2^384; "
              f"entry operands at most {bmax / P:.0f}p - 1 < 2^384\n"
              f"entry: {len(entry)} vectors x 48 bytes = {48 * len(entry)} bytes per table entry\n"
              + bound_table(rows))
    prologue(a, frame)
    a.set_pool(["x0", "x3"], shift=400)
    a.pin_sp(952)
    for key in need:
        L = dict(key)
        saved = dict(derived)
        del derived[key]
        terms = decompose(L)
        derived.update(saved)
        load_lin(a, A, terms)
        a.store6(A, "sp", derived[key])
    for f, (Ax, Bm) in enumerate(bd.fp2):
        e = 144 * f
        product(a, decompose(Ax[0]), "x2", e)
        a.store6(Z[:6], "sp", tmp0)
        a.store6(Z[6:], "sp", tmp0 + 48)
        product(a, decompose(Ax[1]), "x2", e + 48)

        re = atoms + 96 * (2 * f)
        for h in (0, 1):
            a.load6(T, "sp", tmp0 + 48 * h)
            a.chain("addsub", -1, A, T, Z[6 * h:6 * h + 6], first=(h == 0))
            a.store6(A, "sp", re + 48 * h)
        a.acc12(1, "sp", tmp0)
        a.store6(Z[:6], "sp", tmpS)
        a.store6(Z[6:], "sp", tmpS + 48)
        product(a, decompose(M.ladd(Ax[0], Ax[1])), "x2", e + 96)
        a.acc12(-1, "sp", tmpS)
        im = atoms + 96 * (2 * f + 1)
        a.store6(Z[:6], "sp", im)
        a.store6(Z[6:], "sp", im + 48)
    a.ins("mem", "ldp x0, x3, [sp, #96]")
    a.load6(PR, "x3", POOL["p"])
    a.ins("mem", f"ldr {N0R}, [x3, #{POOL['n0']}]")
    a.set_pool(["x1", "x2"], shift=288)
    a.pin_sp(1000)

    tc = [[c for f2 in t for c in f2] for t in (t0, t1, t2)]
    sacc = [[SAVE + 96 * (6 * h + k) for k in range(6)] for h in (0, 1)]
    assert sacc[1][5] + 96 <= tmp0 + 192
    for h in (0, 1):
        for k, c in enumerate(tc[h]):
            items = sorted(c.items())
            at, s = items.pop(next(i for i, (at, s) in enumerate(items) if s > 0))
            a.load6(Z[:6], "sp", atoms + 96 * at)
            a.load6(Z[6:], "sp", atoms + 96 * at + 48)
            for at, s in items:
                assert abs(s) == 1
                a.acc12(s, "sp", atoms + 96 * at)
            a.store6(Z[:6], "sp", sacc[h][k])
            a.store6(Z[6:], "sp", sacc[h][k] + 48)
    vt = [(1, 4, -1, 5), (1, 4, 1, 5), (1, 0), (1, 1), (1, 2), (1, 3)]
    plans = []
    for k in range(6):
        v = vt[k]
        plans.append([(1, sacc[0][k])] + [(v[i], sacc[1][v[i + 1]]) for i in range(0, len(v), 2)])
    for k in range(6):
        plans.append([(s, atoms + 96 * at) for at, s in sorted(tc[2][k].items())] + [(-1, sacc[0][k]), (-1, sacc[1][k])])
    sub = [dict(tc[0][k]) for k in range(6)], [dict(tc[1][k]) for k in range(6)]
    for k, c in enumerate(outs):
        e = {}
        for s, off in plans[k]:
            if off >= atoms:
                src = {(off - atoms) // 96: 1}
            else:
                h, j = divmod((off - SAVE) // 96, 6)
                src = sub[h][j]
            for at, v in src.items():
                e[at] = e.get(at, 0) + s * v
        assert {at: v for at, v in e.items() if v} == {at: v for at, v in c.items() if v}
    consts = {}
    for k, c in enumerate(outs):
        K = rows[k][3]
        items = plans[k]
        if K:
            consts[k] = K
            koff = POOL["K"] + 96 * k
            a.load6(Z[:6], "x3", koff)
            a.load6(Z[6:], "x3", koff + 48)
            rest = items
        else:
            s, off = items.pop(next(i for i, (s, off) in enumerate(items) if s > 0))
            a.load6(Z[:6], "sp", off)
            a.load6(Z[6:], "sp", off + 48)
            rest = items
        for s, off in rest:
            assert abs(s) == 1
            a.acc12(s, "sp", off)
        redc_store(a, k)
    epilogue(a, frame)
    layout = {"entry": [{k: v for k, v in L.items()} for L in entry],
              "K": {str(k): v for k, v in consts.items()}, "frame": frame}
    return a, layout, rows



def sqr_kernel():
    a = Asm()
    bd = M.Builder()
    Av = [[(M.lin(f"a{6*i+2*j}"), M.lin(f"a{6*i+2*j+1}")) for j in range(3)] for i in range(2)]
    vmax = {f"a{i}": P - 1 for i in range(12)}
    prods, outs = M.sqr_model()
    polys = [M.expand(prods, c) for c in outs]
    rows = out_bounds(polys, vmax)


    stored = {}
    scr = SAVE
    order = []
    for Aop, Bop in prods:
        for L in (Bop,):
            key = tuple(sorted(L.items()))
            if len(L) > 1 and key not in stored:
                order.append(L)
                stored[key] = None

    sums = []
    for Aop, Bop in prods:
        for L in (Aop, Bop):
            two = {n: c for n, c in L.items() if n != "1"}
            if len(two) == 2 and all(c == 1 for c in two.values()) and "1" not in L:
                k = tuple(sorted(two.items()))
                if k not in [tuple(sorted(s.items())) for s in sums]:
                    sums.append(two)
    allv = {}
    for L in sums:
        allv[tuple(sorted(L.items()))] = scr
        scr += 48
    for L in order:
        key = tuple(sorted(L.items()))
        if key not in allv:
            allv[key] = scr
            scr += 48

    def decompose(L):
        L = dict(L)
        terms = []
        if "1" in L:
            c = L.pop("1")
            mult = c // P
            assert c == mult * P and mult in (1, 2, 4)
            terms.append((1, "x3", POOL["p" if mult == 1 else f"{mult}p"]))
        neg = []
        for key in sorted(allv, key=lambda k: (len(k), k), reverse=True):
            if len(key) == 2 and all(L.get(n, 0) >= c > 0 for n, c in key):
                for n, c in key:
                    L[n] -= c
                terms.append((1, "sp", allv[key]))
            elif len(key) == 2 and all(L.get(n, 0) <= -c < 0 for n, c in key):
                for n, c in key:
                    L[n] += c
                neg.append((-1, "sp", allv[key]))
        L = {n: c for n, c in L.items() if c}
        for n, c in sorted(L.items(), key=lambda kv: -kv[1]):
            for _ in range(abs(c)):
                (terms if c > 0 else neg).append((1 if c > 0 else -1, "x1", 48 * int(n[1:])))
        return terms + neg

    prod_off = scr
    scr += 96 * len(prods)
    frame = (scr + 15) // 16 * 16
    opmax = 0
    for Aop, Bop in prods:
        for L in (Aop, Bop):
            lo, hi = M.lin_range(L, vmax)
            assert lo >= 0 and hi < (1 << 384)
            opmax = max(opmax, hi)
    a.comment("BLS12-381 cyclotomic squaring (Granger-Scott, blst's cyclotomic_sqr_fp12 formula), Montgomery\n"
              "form, outputs in [0, p); generated by gen.py; args x0 out, x1 x, x3 constant pool\n"
              f"every product operand is an exact non-negative sum below {opmax / P:.2f}p < 2^384\n"
              + bound_table(rows))
    prologue(a, frame)
    a.set_pool(["x0", "x2"], shift=96)
    a.pin_sp(536)
    for key, off in list(allv.items()):
        L = dict(key)
        keep = allv.pop(key)
        terms = decompose(L)
        allv[key] = keep
        load_lin(a, A, terms)
        a.store6(A, "sp", off)
    a.keep_pool(["x0", "x2", "x3"], 96)
    a.pin_sp(1000)
    for i, (Aop, Bop) in enumerate(prods):
        bkey = tuple(sorted(Bop.items()))
        if bkey in allv:
            bb, bo = "sp", allv[bkey]
        else:
            (n, c), = Bop.items()
            assert c == 1
            bb, bo = "x1", 48 * int(n[1:])
        double = all(v == 2 for v in Aop.values())
        aterms = decompose({n: 1 for n in Aop}) if double else decompose(Aop)
        product(a, aterms, bb, bo, double_a=double)
        a.store6(Z[:6], "sp", prod_off + 96 * i)
        a.store6(Z[6:], "sp", prod_off + 96 * i + 48)
    a.ins("mem", "ldp x0, x3, [sp, #96]")
    a.load6(PR, "x3", POOL["p"])
    a.ins("mem", f"ldr {N0R}, [x3, #{POOL['n0']}]")
    a.keep_pool(["x2"], 96, drop=("x0", "x3"))
    a.pin_base("x3", 1000)
    consts = {}
    for k in [10, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 11]:
        c = outs[k]
        prodterms = sorted((pi, v // 3) for pi, v in c.items() if not isinstance(pi, tuple))
        assert all(v * 3 == c[pi] for pi, v in prodterms)
        hi_terms = [(key[1], v) for key, v in c.items() if isinstance(key, tuple)]
        first = next(i for i, (pi, s) in enumerate(prodterms) if s > 0)
        pi, s = prodterms.pop(first)
        a.load6(Z[:6], "sp", prod_off + 96 * pi)
        a.load6(Z[6:], "sp", prod_off + 96 * pi + 48)
        for pi, s in prodterms:
            for _ in range(abs(s)):
                a.acc12(1 if s > 0 else -1, "sp", prod_off + 96 * pi)
        a.ins("addsub", f"add {T[0]}, {Z[0]}, {Z[0]}")
        for i in range(12):
            if i < 11:
                a.ins("addsub", f"extr {T[1 - i % 2]}, {Z[i + 1]}, {Z[i]}, #63")
            a.ins("addsub", f"{'adcs' if i else 'adds'} {Z[i]}, {Z[i]}, {T[i % 2]}")
        for name, v in hi_terms:
            a.load6(T, "x1", 48 * int(name[1:]))
            for _ in range(abs(v)):
                a.chain("addsub", 1 if v > 0 else -1, Z[6:], Z[6:], T)
        K = rows[k][3]
        if K:
            consts[k] = K
            a.acc12(1, "x3", POOL["K"] + 96 * k)
        redc_store(a, k)
    epilogue(a, frame)
    return a, {"K": {str(k): v for k, v in consts.items()}, "frame": frame}, rows


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    rep = {}
    for name, fn in (("mul", mul_kernel), ("sqr", sqr_kernel)):
        a, layout, rows = fn()
        with open(os.path.join(here, f"{name}.s"), "w") as f:
            f.write(a.text())
        layout["count"] = a.count()
        layout["executed"] = sum(layout["count"].values())
        layout["pool"] = {"n0": N0, "qc": QC, "offsets": POOL}
        with open(os.path.join(here, f"{name}.layout.json"), "w") as f:
            json.dump(layout, f, indent=1, default=str)
        rep[name] = (layout["executed"], layout["count"])
        print(name, layout["executed"], layout["count"])
    return rep


if __name__ == "__main__":
    main()
