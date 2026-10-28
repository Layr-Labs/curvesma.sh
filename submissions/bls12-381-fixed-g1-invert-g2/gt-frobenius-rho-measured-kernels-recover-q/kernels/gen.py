"""Emit straight-line scalar AArch64 for two BLS12-381 Fp12 kernels (blst's tower, Montgomery form,
6 x 64-bit limbs per Fp element, outputs fully reduced):

  mul.s  out = x * m, m a fixed table entry stored as 54 table-side linear forms reduced mod p
  sqr.s  out = x^2 for x in the cyclotomic subgroup (Granger-Scott, as blst's cyclotomic_sqr_fp12)

Calling convention of both: x0 = out (12 x 48 bytes), x1 = x (12 x 48 bytes), x2 = table entry (mul),
x3 = constant pool (layout.json). General registers only, no branches, no calls; ends in ret.

Lazy reduction: each output coordinate is one 768-bit integer T = K + sum(+-products) (+ 2 a R terms
for the squaring), operand sums exact, so T is an integer polynomial whose range model.py bounds monomial
by monomial (K: least multiple of p making T >= 0); one Montgomery reduction, then subtract q p and p.
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
PR, N0R, MR = A, "x23", "x24"          # registers of the reduction phase
SAVE = 112                               # callee-saved x19..x30 (96 bytes) and x0, x3 (16 bytes)

# constant pool: the q*p table at 64q (q = 0..7), n0, c, p, 2p, 4p, then the per-kernel offsets K_k
POOL = {"qt": 0, "n0": 48, "c": 56, "p": 64, "2p": 128, "4p": 256, "K": 512}
QC = (1 << 64) // ((P >> 320) + 1)      # q = r5 * QC >> 64 is monotone in r: check both ends of each [kp, (k+1)p)
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

    def ref(self, base, off, pair=True):
        lo, hi = (-512, 504) if pair else (0, 32760)
        if lo <= off <= hi:
            return f"[{base}, #{off}]"
        self.clock += 1
        for r, (b, o, _) in list(self.addr.items()):
            if b == base and lo <= off - o <= hi:
                self.addr[r] = (b, o, self.clock)
                return f"[{r}, #{off - o}]"
        free = [r for r in self.pool if r not in self.addr]
        r = free[0] if free else min(self.addr, key=lambda k: self.addr[k][2])
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

def product(a, aterms, bbase, boff, double_a=False, bregs=None):
    """Z = (sum of aterms) * (vector at bbase+boff), schoolbook by rows of b (143 arithmetic instructions);
    with bregs, A and b are already in A and bregs (b_j must survive rows 0..j: row j writes Z[j..j+6], T)."""
    if not bregs:
        load_lin(a, A, aterms)
    if double_a:
        a.chain("addsub", 1, A, A, A)
    for jj in range(0, 6, 2):
        if not bregs:
            a.ins("mem", f"ldp {B[0]}, {B[1]}, {a.ref(bbase, boff + 8 * jj)}")
        for j, b in (zip((jj, jj + 1), bregs[jj:jj + 2]) if bregs else ((jj, B[0]), (jj + 1, B[1]))):
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
    assert W == Z[:6]                             # r < 8p; q = floor(r5 * QC / 2^64) in [floor(r/p) - 1, floor(r/p)]
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

def out_bounds(polys, vmax, unit=P, share=None):
    """K = least multiple of unit with K + lo >= 0; outputs in one group of `share` take the group's largest K."""
    lohi = [M.poly_range(poly, vmax) for poly in polys]
    Ks = [-(-(-lo) // unit) * unit if lo < 0 else 0 for lo, hi in lohi]
    for g in share or []:
        for k in g:
            Ks[k] = max(Ks[i] for i in g)
    rows = []
    for k, ((lo, hi), K) in enumerate(zip(lohi, Ks)):
        top = K + hi
        assert K + lo >= 0 and top < (1 << 768)
        bres = top // R + 1 + P            # r < T/R + p
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

# -- multiplication by a table entry
def fp2_fixed(Ax, Bm):
    """A * B, B fixed: products (A0+A1) B1, A0 (B0+B1), A1 (B0-B1); re = P1 - P0, im = P2 + P0."""
    return [M.ladd(*Ax), Ax[0], Ax[1]], [Bm[1], M.ladd(*Bm), M.ladd(Bm[0], Bm[1], -1)]

def m2(*terms):
    """Fp2 linear combination of table-side pairs, (sign, pair) terms."""
    s, (b0, b1) = terms[0]
    for t, (c0, c1) in terms[1:]:
        b0, b1 = M.ladd(b0, c0, t), M.ladd(b1, c1, t)
    return b0, b1

def xim(b):
    return M.ladd(b[0], b[1], -1), M.ladd(b[0], b[1])

FP6_OUT = ((3, 0, 1), (4, 0, 2), (5, 1, 2))     # o_j = Q_j + two of PA, PB, PC (Fp2 groups of an instance)

def mul_plan():
    """Fp12 x * m with every table-side form free: k1 = (X0+X1) M0, k2 = X1 (v M1 - M0), k3 = X0 (M1 - M0),
    c0 = k1 + k2, c1 = k1 + k3; Fp6 y * n by six Fp2 products PA = (y1+y2) xi n2, PB = (y0+y2) n0,
    PC = (y0+y1) n1, Q0 = y2 (xi n1 - xi n2 - n0), Q1 = y1 (n0 - xi n2 - n1), Q2 = y0 (n2 - n0 - n1):
    o0 = Q0+PA+PB, o1 = Q1+PA+PC, o2 = Q2+PB+PC. Returns the 54 products (x-side Lin, table-side Lin)."""
    x0, x1 = M.fp6_vars("x", 0), M.fp6_vars("x", 6)
    n0, n1 = M.fp6_vars("m", 0), M.fp6_vars("m", 6)
    vn1 = [xim(n1[2]), n1[0], n1[1]]
    insts = [(M.fp6_lsum(x0, x1), n0), (x1, [m2((1, vn1[j]), (-1, n0[j])) for j in range(3)]),
             (x0, [m2((1, n1[j]), (-1, n0[j])) for j in range(3)])]
    prods = []
    for y, n in insts:
        xn1, xn2 = xim(n[1]), xim(n[2])
        for Ax, Bm in ((M.l2add(y[1], y[2]), xn2), (M.l2add(y[0], y[2]), n[0]), (M.l2add(y[0], y[1]), n[1]),
                       (y[2], m2((1, xn1), (-1, xn2), (-1, n[0]))), (y[1], m2((1, n[0]), (-1, xn2), (-1, n[1]))),
                       (y[0], m2((1, n[2]), (-1, n[0]), (-1, n[1])))):
            prods += list(zip(*fp2_fixed(Ax, Bm)))
    return prods

def mul_outs():
    """Each output coordinate as {product index: sign}."""
    def kc(i, o, k):
        d = {}
        for g in FP6_OUT[o]:
            p = 3 * (6 * i + g)
            d.update({p + 1: 1, p: -1} if k == 0 else {p + 2: 1, p: 1})
        return d
    return [{**kc(0, o, k), **kc(h, o, k)} for h in (1, 2) for o in range(3) for k in range(2)]

def mul_kernel():
    a = Asm()
    prods, outs = mul_plan(), mul_outs()
    assert len(prods) == 54 and all(len(c) == 12 for c in outs)
    mod = lambda q: {k: v % P for k, v in q.items() if v % P}
    fprods, fouts = M.mul_model()                # the reference Karatsuba model, coefficient by coefficient mod p
    assert all(mod(M.expand(prods, c)) == mod(M.expand(fprods, f)) for c, f in zip(outs, fouts))
    # ranges: each stored table value is reduced, an independent input in [0, p)
    eprods = [(xs, {f"e{i}": 1}) for i, (xs, _) in enumerate(prods)]
    vmax = {f"x{i}": P - 1 for i in range(12)}
    vmax.update({f"e{i}": P - 1 for i in range(54)})
    rows = out_bounds([M.expand(eprods, c) for c in outs], vmax, P * R, [(j, j + 6) for j in range(6)])
    assert all(r[3] == 0 for r in rows[1::2]) and all(r[3] <= 7 * P * R for r in rows)

    need = sorted({tuple(sorted(xs.items())) for p, (xs, _) in enumerate(prods) if p % 3 and len(xs) > 1},
                  key=lambda k: (len(k), k))
    derived = {key: SAVE + 48 * i for i, key in enumerate(need)}

    def decompose(L, skip=None):
        L, terms = dict(L), []
        for key in sorted(derived, key=lambda k: (len(k), k), reverse=True):
            if key != skip and all(L.get(n, 0) >= c for n, c in key):
                for n, c in key:
                    L[n] -= c
                L = {n: c for n, c in L.items() if c}
                terms.append((1, "sp", derived[key]))
        assert all(c == 1 for c in L.values())
        return terms + [(1, "x1", 48 * int(n[1:])) for n in sorted(L)]

    tmp = SAVE + 48 * len(need)
    atoms = tmp + 96                             # PA, PB, PC re/im of the current instance
    kslot = atoms + 6 * 96                       # o_j re/im of k1 (k2, k3 outputs reduce in place)
    frame = (kslot + 6 * 96 + 15) // 16 * 16
    opmax = max(M.lin_range(xs, vmax)[1] for xs, _ in prods)
    assert opmax < (1 << 384)
    a.comment("BLS12-381 Fp12 multiplication by a fixed table entry, Montgomery form, outputs in [0, p)\n"
              "generated by gen.py; args x0 out, x1 x, x2 entry, x3 constant pool\n"
              f"x-side operands are exact sums of reduced inputs, at most {opmax / P:.0f}p - 1 < 2^384; "
              "entry values reduced mod p\n"
              f"entry: {len(prods)} vectors x 48 bytes = {48 * len(prods)} bytes per table entry\n"
              + bound_table(rows))
    prologue(a, frame)
    a.set_pool(["x0", "x3"])
    for key in need:                             # x-side pre-sums
        load_lin(a, A, decompose(dict(key), key))
        a.store6(A, "sp", derived[key])

    def prod(p):
        product(a, decompose(prods[p][0]), "x2", 48 * p)

    def store12(off):
        a.store6(Z[:6], "sp", off)
        a.store6(Z[6:], "sp", off + 48)

    for i in range(3):
        for g in range(3):
            p = 3 * (6 * i + g)
            prod(p)
            store12(tmp)
            for k in (0, 1):
                prod(p + 1 + k)
                a.acc12(1 if k else -1, "sp", tmp)
                store12(atoms + 96 * (2 * g + k))
        for o, (q, x, y) in enumerate(FP6_OUT):
            p = 3 * (6 * i + q)
            prod(p)
            store12(tmp)
            for k in (0, 1):
                prod(p + 1 + k)
                a.acc12(1 if k else -1, "sp", tmp)
                a.acc12(1, "sp", atoms + 96 * (2 * x + k))
                a.acc12(1, "sp", atoms + 96 * (2 * y + k))
                if i == 0 and k == 0:            # K = j p 2^384 into the upper half, shared by c0 and c1
                    a.ins("mem", "ldr x3, [sp, #104]")
                    a.addr.pop("x3", None)
                    a.load6(T, "x3", 64 * (rows[2 * o][3] // (P * R)))
                    a.chain("addsub", 1, Z[6:], Z[6:], T)
                if i == 0:
                    store12(kslot + 96 * (2 * o + k))
                    continue
                a.acc12(1, "sp", kslot + 96 * (2 * o + k))   # c = k1 + k2 (i = 1) or k1 + k3 (i = 2)
                a.ins("mem", "ldp x0, x3, [sp, #96]")
                a.load6(PR, "x3", POOL["p"])
                a.ins("mem", f"ldr {N0R}, [x3, #{POOL['n0']}]")
                redc_store(a, 6 * (i - 1) + 2 * o + k)
                a.set_pool(["x0", "x3"])
    epilogue(a, frame)
    entry = [ms for _, ms in prods]
    assert sum(abs(c) for f in entry for c in f.values()) - len(entry) == 270  # modular adds per entry at setup
    return a, {"entry": entry, "xs": [xs for xs, _ in prods], "outs": outs, "frame": frame, "K": {}}, rows

# -- cyclotomic squaring
def sqr_kernel():
    a = Asm()
    bd = M.Builder()
    Av = [[(M.lin(f"a{6*i+2*j}"), M.lin(f"a{6*i+2*j+1}")) for j in range(3)] for i in range(2)]
    vmax = {f"a{i}": P - 1 for i in range(12)}
    prods, outs = M.sqr_model()
    polys = [M.expand(prods, c) for c in outs]
    rows = out_bounds(polys, vmax)

    # stored vectors: Fp4 sums s = a + b (two Fp coordinates each) and every B operand sum
    stored = {}
    scr = SAVE
    order = []
    for Aop, Bop in prods:
        for L in (Bop,):
            key = tuple(sorted(L.items()))
            if len(L) > 1 and key not in stored:
                order.append(L)
                stored[key] = None
    # the Fp4 sums are inside B operands of the t2 squarings; store them first
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
    a.set_pool(["x0", "x2"])
    for key, off in list(allv.items()):          # stored sums (x3 still the constant pool here)
        L = dict(key)
        keep = allv.pop(key)
        terms = decompose(L)
        allv[key] = keep
        load_lin(a, A, terms)
        a.store6(A, "sp", off)
    a.set_pool(["x0", "x2", "x3"])
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
    a.set_pool(["x2"])
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
        a.ins("addsub", f"add {T[0]}, {Z[0]}, {Z[0]}")   # 3Z = Z + 2Z, 2Z by extr (no flags)
        for i in range(12):
            if i < 11:
                a.ins("addsub", f"extr {T[1 - i % 2]}, {Z[i + 1]}, {Z[i]}, #63")
            a.ins("addsub", f"{'adcs' if i else 'adds'} {Z[i]}, {Z[i]}, {T[i % 2]}")
        for name, v in hi_terms:                   # +-2 a in the upper half
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
