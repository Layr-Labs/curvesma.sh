"""Emit straight-line scalar AArch64 for two BLS12-381 Fp12 kernels (blst's tower, Montgomery form,
6 x 64-bit limbs per Fp element, outputs fully reduced):

  mul.s  out = x * m, m a fixed table entry stored as 54 vectors, each the residue mod p of a linear form in m
  sqr.s  out = x^2 for x in the cyclotomic subgroup (Granger-Scott, as blst's cyclotomic_sqr_fp12)

mul: table-side forms at all three tower levels (proof.md section 5): Fp12 c0 = k1 + k2, c1 = k1 + k3 with
k1 = (X0 + X1) M0, k2 = X1 (v M1 - M0), k3 = X0 (M1 - M0); each Fp6 product y b in six Fp2 products
PA = (y1 + y2) xi b2, PB = (y0 + y2) b0, PC = (y0 + y1) b1, Q0 = y2 (xi b1 - xi b2 - b0), Q1 = y1 (b0 - xi b2 - b1),
Q2 = y0 (b2 - b0 - b1), o0 = Q0 + PA + PB, o1 = Q1 + PA + PC, o2 = Q2 + PB + PC; each Fp2 product a c in three
products p0 = (a0 + a1) c0, p1 = a1 (-(c0 + c1)), p2 = a0 (c1 - c0), re = p0 + p1, im = p0 + p2. Every m-side value
is stored reduced, so each output's T is a sum of 12 non-negative products (K = 0, asserted) congruent to the
reference mod p (asserted). entry_program gives the 113 modular add/sub that build the 54 vectors from m.

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
A_REGS = A
ENTRY_UNITS = 12
ENTRY_OPS = 113
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
        r = free[0] if free else min((k for k in self.addr if k != "x18"), key=lambda k: self.addr[k][2])
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



def l2(a, b, s=1):
    return (M.ladd(a[0], b[0], s), M.ladd(a[1], b[1], s))


def neg(L):
    return {n: -c for n, c in L.items()}


def fp6_form(b):
    """y * b for a fixed Fp6 b in 6 Fp2 products: [(x-side operand as y index pair or one index, table Fp2
    value, outputs)]. o0 = Q0 + PA + PB, o1 = Q1 + PA + PC, o2 = Q2 + PB + PC with PA = (y1 + y2) xi b2,
    PB = (y0 + y2) b0, PC = (y0 + y1) b1, Q0 = y2 (xi b1 - xi b2 - b0), Q1 = y1 (b0 - xi b2 - b1),
    Q2 = y0 (b2 - b0 - b1)."""
    xb1, xb2 = M.xi(b[1]), M.xi(b[2])
    return [((1, 2), xb2, (0, 1)), ((0, 2), b[0], (0, 2)), ((0, 1), b[1], (1, 2)),
            ((2,), l2(l2(xb1, xb2, -1), b[0], -1), (0,)), ((1,), l2(l2(b[0], xb2, -1), b[1], -1), (1,)),
            ((0,), l2(l2(b[2], b[0], -1), b[1], -1), (2,))]


def entry_program(entry):
    """Node construction of the 54 entry vectors from the 12 reduced coordinates: op i forms value 12 + i as
    (value a, or 0 when a = -1) + s * (value b) mod p from earlier values (values 0..11 are m0..m11), following the
    kernel's own recipe and reusing every value already formed. Returns the ops and the value index of each vector;
    every vector's form is checked exactly against the entry."""
    key = lambda d: tuple(sorted((n, c) for n, c in d.items() if c))
    forms, idx, ops = [], {}, []

    def put(k, op):
        if op:
            ops.append(op)
        idx[k] = len(forms)
        forms.append(k)
        return k

    def get(k, neg_ok=False):
        if k in idx:
            return k
        for a in list(idx):
            for s in (1, -1):
                b = key({n: s * c for n, c in M.ladd(dict(k), dict(a), -1).items()})
                if b in idx:
                    return put(k, (idx[a], idx[b], s))
        assert neg_ok, k
        b = get(key(neg(dict(k))))
        return put(k, (-1, idx[b], -1))

    add = lambda a, b, s=1: get(key(M.ladd(dict(a), dict(b), s)))
    f2 = lambda a, b, s=1: (add(a[0], b[0], s), add(a[1], b[1], s))
    xi = lambda u: (add(u[0], u[1], -1), add(u[0], u[1]))
    m = [put(key({f"m{i}": 1}), None) for i in range(12)]
    m0 = [(m[2 * j], m[2 * j + 1]) for j in range(3)]
    m1 = [(m[6 + 2 * j], m[7 + 2 * j]) for j in range(3)]
    vm1 = [xi(m1[2]), m1[0], m1[1]]
    out = []
    for b in (m0, [f2(vm1[j], m0[j], -1) for j in range(3)], [f2(m1[j], m0[j], -1) for j in range(3)]):
        xb1, xb2 = xi(b[1]), xi(b[2])
        for c in (xb2, b[0], b[1], f2(f2(xb1, xb2, -1), b[0], -1), f2(f2(b[0], xb2, -1), b[1], -1),
                  f2(f2(b[2], b[0], -1), b[1], -1)):
            s = key(M.ladd(dict(c[0]), dict(c[1])))
            if key(neg(dict(s))) not in idx:
                s = get(s, True) if s not in idx else s
            out += [c[0], get(key(neg(dict(s))), True), add(c[1], c[0], -1)]
    assert out == [key(L) for L in entry]
    return ops, [idx[k] for k in out]


def mul_kernel():
    a = Asm()
    m0, m1 = M.fp6_vars("m", 0), M.fp6_vars("m", 6)
    vm1 = [M.xi(m1[2]), m1[0], m1[1]]
    sixes = [("Y", m0, (0, 1)), ("x1", [l2(vm1[j], m0[j], -1) for j in range(3)], (0,)),
             ("x0", [l2(m1[j], m0[j], -1) for j in range(3)], (1,))]

    src = {}                                   # name -> (base, offset, Lin over x)
    for n in range(12):
        src[f"x{n}"] = ("x1", 48 * n, {f"x{n}": 1})
    scr = SAVE
    derived = []
    for n in range(6):
        derived.append((f"Y{n}", [f"x{n}", f"x{n + 6}"]))
    for op, base in (("x0", 0), ("x1", 6)):
        for k in range(3):
            derived.append((f"s{op}{k}", [f"x{base + 2 * k}", f"x{base + 2 * k + 1}"]))
    for k in range(3):
        derived.append((f"sY{k}", [f"sx0{k}", f"sx1{k}"]))
    for name, parts in derived:
        L = {}
        for q in parts:
            L = M.ladd(L, src[q][2])
        src[name] = ("sp", scr, L)
        scr += 48
    coord = {"Y": [f"Y{n}" for n in range(6)], "x0": [f"x{n}" for n in range(6)],
             "x1": [f"x{n}" for n in range(6, 12)]}
    ssum = {op: [f"s{op}{k}" for k in range(3)] for op in coord}
    tmp = scr
    scr += 96
    S = [scr + 96 * c for c in range(6)]
    scr += 96 * 6
    C = [scr + 96 * c for c in range(12)]
    scr += 96 * 12
    frame = (scr + 15) // 16 * 16

    prods, entry, plan = [], [], []            # plan: per product, its A source names and (re, im) targets
    outs = [[] for _ in range(12)]
    for f6, (op, b, halves) in enumerate(sixes):
        for idx, c, js in fp6_form(b):
            y = coord[op]
            a0 = [y[2 * i] for i in idx]
            a1 = [y[2 * i + 1] for i in idx]
            s01 = [ssum[op][i] for i in idx]
            tre = [(f6, 2 * j) for j in js]
            tim = [(f6, 2 * j + 1) for j in js]
            i0 = len(prods)
            for names, L in ((s01, c[0]), (a1, neg(M.ladd(c[0], c[1]))), (a0, M.ladd(c[1], c[0], -1))):
                A = {}
                for q in names:
                    A = M.ladd(A, src[q][2])
                prods.append((A, len(entry)))
                entry.append(L)
                plan.append(names)
            for h in halves:
                for j in js:
                    outs[6 * h + 2 * j] += [i0, i0 + 1]
                    outs[6 * h + 2 * j + 1] += [i0, i0 + 2]
            plan[i0] = (plan[i0], None)
            plan[i0 + 1] = (plan[i0 + 1], tre)
            plan[i0 + 2] = (plan[i0 + 2], tim)
    assert len(prods) == 54 and all(len(o) == 12 for o in outs)

    sym = [M.expand([(A, entry[e]) for A, e in prods], {pi: 1 for pi in o}) for o in outs]
    fprods, fouts = M.mul_model()
    for want, got in zip([M.expand(fprods, c) for c in fouts], sym):
        assert all(v % P == 0 for v in M.ladd(got, want, -1).values())     # congruent to x * m mod p
    vmax = {f"x{i}": P - 1 for i in range(12)}
    vmax.update({f"e{i}": P - 1 for i in range(54)})
    polys = [M.expand([(A, {f"e{e}": 1}) for A, e in prods], {pi: 1 for pi in o}) for o in outs]
    rows = out_bounds(polys, vmax)
    assert all(r[1] == 0 and r[3] == 0 and r[2] < 45 * P * P for r in rows)   # every term >= 0: no K offset
    units = max(sum(abs(v) for v in L.values()) for L in entry)
    assert all("1" not in L for L in entry) and units <= ENTRY_UNITS
    opmax = max(M.lin_range(A, vmax)[1] for A, _ in prods)
    assert opmax <= 8 * (P - 1) and 8 * P < (1 << 384)

    a.comment("BLS12-381 Fp12 multiplication by a fixed table entry, Montgomery form, outputs in [0, p)\n"
              "generated by gen.py; args x0 out, x1 x, x2 entry, x3 constant pool\n"
              f"x-side operands are exact sums of reduced inputs, at most {opmax // (P - 1)}(p - 1) < 2^384; "
              "entry values reduced mod p\n"
              f"entry: {len(entry)} vectors x 48 bytes = {48 * len(entry)} bytes per table entry, each the residue "
              f"of an integer form of at most {units} units in m\n"
              + bound_table(rows))
    prologue(a, frame)
    a.set_pool(["x0", "x3"], shift=400)
    a.pin_sp(1000)
    for name, parts in derived:
        load_lin(a, A_REGS, [(1, src[q][0], src[q][1]) for q in parts])
        a.store6(A_REGS, "sp", src[name][1])

    slot = {(0, c): S[c] for c in range(6)}
    slot.update({(1, c): C[c] for c in range(6)})
    slot.update({(2, c): C[6 + c] for c in range(6)})
    held = {}                                  # slot offset -> product indices summed there
    zsum = []

    def deposit(t):
        off = slot[t]
        if off not in held:
            a.store6(Z[:6], "sp", off)
            a.store6(Z[6:], "sp", off + 48)
            held[off] = list(zsum)
            return
        for h in (0, 1):
            a.load6(T, "sp", off + 48 * h)
            a.chain("addsub", 1, A_REGS, Z[6 * h:6 * h + 6], T, first=(h == 0))
            a.store6(A_REGS, "sp", off + 48 * h)
        held[off] = held[off] + zsum

    for pi, ((A, e), (names, tg)) in enumerate(zip(prods, plan)):
        product(a, [(1, src[q][0], src[q][1]) for q in names], "x2", 48 * e)
        zsum = [pi]
        if tg is None:
            a.store6(Z[:6], "sp", tmp)
            a.store6(Z[6:], "sp", tmp + 48)
            held[tmp] = [pi]
            continue
        a.acc12(1, "sp", tmp)
        zsum = zsum + held[tmp]
        for t in tg:
            deposit(t)
    a.ins("mem", "ldp x0, x3, [sp, #96]")
    a.load6(PR, "x3", POOL["p"])
    a.ins("mem", f"ldr {N0R}, [x3, #{POOL['n0']}]")
    a.set_pool(["x1", "x2"], shift=400)
    a.pin_sp(2000)
    for k in range(12):
        h, c = divmod(k, 6)
        off = slot[(h + 1, c)]
        a.load6(Z[:6], "sp", off)
        a.load6(Z[6:], "sp", off + 48)
        a.acc12(1, "sp", slot[(0, c)])
        assert sorted(held[off] + held[slot[(0, c)]]) == sorted(outs[k])
        redc_store(a, k)
    epilogue(a, frame)
    ops, index = entry_program(entry)
    assert len(ops) == ENTRY_OPS and all(-1 <= a < 12 + i and 0 <= b < 12 + i for i, (a, b, s) in enumerate(ops))
    layout = {"entry": [{n: v for n, v in L.items()} for L in entry], "entry_reduced_mod_p": True,
              "entry_units_max": units, "entry_ops": ops, "entry_index": index, "K": {}, "frame": frame}
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
    a.set_pool(["x0", "x2", "x3"], shift=96)
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
    a.set_pool(["x2"], shift=96)
    a.pin_sp(2200)
    consts = {}
    for k, c in enumerate(outs):
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
