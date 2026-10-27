"""Symbolic models: the squaring kernel's products and the Karatsuba reference the multiply plan of gen.py
is checked against, and how each of the 12 output coordinates combines them, before any reduction. Tower as blst: Fp2 = Fp[u]/(u^2+1),
Fp6 = Fp2[v]/(v^3-(1+u)), Fp12 = Fp6[w]/(w^2-v); coordinate index 6*i + 2*j + k for c[i][j][k].

Lin: dict name -> int, a linear form over single-width input vectors ('x3', 'm7', 'a0', and '1' for
the constant 1, scaled by a multiple of p). Comb: dict key -> int, an output's linear combination of
products (key = product index) and of input vectors placed in the upper half (key = ('hi', name)).
"""
P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 1 << 384

def lin(*names):
    d = {}
    for n in names:
        d[n] = d.get(n, 0) + 1
    return d

def ladd(a, b, s=1):
    d = dict(a)
    for k, v in b.items():
        d[k] = d.get(k, 0) + s * v
        if d[k] == 0:
            del d[k]
    return d

def cadd(a, b, s=1):
    return ladd(a, b, s)

class Builder:
    def __init__(self):
        self.prods = []      # list of (A Lin, B Lin)
        self.index = {}

    def prod(self, A, B):
        key = (tuple(sorted(A.items())), tuple(sorted(B.items())))
        if key not in self.index:
            self.index[key] = len(self.prods)
            self.prods.append((A, B))
        return {self.index[key]: 1}

    # Fp2 combos (pairs of Comb)
    def fp2_mul(self, A, B):
        t0 = self.prod(A[0], B[0])
        t1 = self.prod(A[1], B[1])
        t2 = self.prod(ladd(A[0], A[1]), ladd(B[0], B[1]))
        return (cadd(t0, t1, -1), cadd(cadd(t2, t0, -1), t1, -1))

def f2add(a, b, s=1):
    return (cadd(a[0], b[0], s), cadd(a[1], b[1], s))

def xi(c):
    return (cadd(c[0], c[1], -1), cadd(c[0], c[1]))

def l2add(a, b):
    return (ladd(a[0], b[0]), ladd(a[1], b[1]))

def fp6_mul(bd, a, b):
    v0, v1, v2 = bd.fp2_mul(a[0], b[0]), bd.fp2_mul(a[1], b[1]), bd.fp2_mul(a[2], b[2])
    t12 = bd.fp2_mul(l2add(a[1], a[2]), l2add(b[1], b[2]))
    t01 = bd.fp2_mul(l2add(a[0], a[1]), l2add(b[0], b[1]))
    t02 = bd.fp2_mul(l2add(a[0], a[2]), l2add(b[0], b[2]))
    c0 = f2add(v0, xi(f2add(f2add(t12, v1, -1), v2, -1)))
    c1 = f2add(f2add(f2add(t01, v0, -1), v1, -1), xi(v2))
    c2 = f2add(f2add(f2add(t02, v0, -1), v2, -1), v1)
    return [c0, c1, c2]

def fp6_vars(pfx, base):
    return [(lin(f"{pfx}{base + 2 * j}"), lin(f"{pfx}{base + 2 * j + 1}")) for j in range(3)]

def fp6_lsum(a, b):
    return [l2add(a[j], b[j]) for j in range(3)]

def mul_model():
    """x * m, Karatsuba at the Fp12, Fp6 and Fp2 levels: 3 x 6 x 3 = 54 products."""
    bd = Builder()
    x0, x1 = fp6_vars("x", 0), fp6_vars("x", 6)
    m0, m1 = fp6_vars("m", 0), fp6_vars("m", 6)
    t0 = fp6_mul(bd, x0, m0)
    t1 = fp6_mul(bd, x1, m1)
    t2 = fp6_mul(bd, fp6_lsum(x0, x1), fp6_lsum(m0, m1))
    c1 = [f2add(f2add(t2[j], t0[j], -1), t1[j], -1) for j in range(3)]
    vt1 = [xi(t1[2]), t1[0], t1[1]]
    c0 = [f2add(t0[j], vt1[j]) for j in range(3)]
    outs = [c for f6 in (c0, c1) for f2 in f6 for c in f2]
    return bd.prods, outs

def sqr_model():
    """Granger-Scott cyclotomic squaring as blst's cyclotomic_sqr_fp12: three Fp4 squarings of
    (a00, a11), (a10, a02), (a01, a12), each from three Fp2 squarings; an Fp2 square of (s0, s1) is
    ((s0+s1)(s0-s1), 2 s0 s1), with s0 - s1 taken as s0 + (k p - s1) to stay non-negative."""
    bd = Builder()
    A = [[(lin(f"a{6*i+2*j}"), lin(f"a{6*i+2*j+1}")) for j in range(3)] for i in range(2)]

    def kp_minus(s, bound_mult):
        return ladd({"1": bound_mult * P}, s, -1)

    def fp2_sqr(s, smax):
        # smax: s0, s1 < smax * p as integers
        c0 = bd.prod(ladd(s[0], s[1]), ladd(s[0], kp_minus(s[1], smax)))
        c1 = bd.prod(ladd(s[0], s[0]), s[1])
        return (c0, c1)

    def fp4_sqr(a0, a1):
        t0 = fp2_sqr(a0, 1)
        t1 = fp2_sqr(a1, 1)
        s = l2add(a0, a1)
        t2 = fp2_sqr(s, 2)
        r0 = f2add(xi(t1), t0)
        r1 = f2add(f2add(t2, t0, -1), t1, -1)
        return r0, r1

    t0 = fp4_sqr(A[0][0], A[1][1])
    t1 = fp4_sqr(A[1][0], A[0][2])
    t2 = fp4_sqr(A[0][1], A[1][2])

    def three(c):
        return {k: 3 * v for k, v in c.items()}

    def out(t, a, sign):
        # 3 t + 2 sign a, a placed in the upper half (a * R reduces to a)
        res = []
        for k in range(2):
            c = three(t[k])
            name = next(iter(a[k]))
            c[("hi", name)] = 2 * sign
            res.append(c)
        return res

    o = {}
    o[(0, 0)] = out(t0[0], A[0][0], -1)
    o[(0, 1)] = out(t1[0], A[0][1], -1)
    o[(0, 2)] = out(t2[0], A[0][2], -1)
    o[(1, 0)] = out(xi(t2[1]), A[1][0], +1)
    o[(1, 1)] = out(t0[1], A[1][1], +1)
    o[(1, 2)] = out(t1[1], A[1][2], +1)
    outs = [c for i in range(2) for j in range(3) for c in o[(i, j)]]
    return bd.prods, outs

def lin_range(L, vmax):
    """[lo, hi] of a Lin with each named input in [0, vmax[name]] and '1' the constant 1."""
    lo = hi = 0
    for k, c in L.items():
        if k == "1":
            lo += c
            hi += c
        else:
            (hi := hi)
            if c > 0:
                hi += c * vmax[k]
            else:
                lo += c * vmax[k]
    return lo, hi

def expand(prods, comb):
    """The exact polynomial of one output: dict monomial -> coefficient; monomial = sorted tuple of
    input names ('1' dropped). Exact over the integers because every sum is taken without reduction."""
    poly = {}

    def addm(mono, c):
        mono = tuple(sorted(n for n in mono if n != "1"))
        poly[mono] = poly.get(mono, 0) + c
        if poly[mono] == 0:
            del poly[mono]

    for key, c in comb.items():
        if isinstance(key, tuple) and key[0] == "hi":
            addm((key[1],), c * R)
            continue
        A, B = prods[key]
        for ka, ca in A.items():
            for kb, cb in B.items():
                addm((ka, kb), c * ca * cb)
    return poly

def poly_range(poly, vmax):
    """Rigorous [lo, hi] of a polynomial with non-negative inputs bounded by vmax: every monomial
    lies in [0, its maximum], so positive coefficients bound hi and negative ones bound lo."""
    lo = hi = 0
    for mono, c in poly.items():
        mx = 1
        for n in mono:
            mx *= vmax[n]
        if c > 0:
            hi += c * mx
        else:
            lo += c * mx
    return lo, hi

if __name__ == "__main__":
    for name, (prods, outs) in (("mul", mul_model()), ("sqr", sqr_model())):
        inc = sum(sum(abs(v) for k, v in c.items() if not isinstance(k, tuple)) for c in outs)
        nz = sum(sum(1 for k, v in c.items() if not isinstance(k, tuple)) for c in outs)
        uses = {}
        for c in outs:
            for k in c:
                if not isinstance(k, tuple):
                    uses[k] = uses.get(k, 0) + 1
        print(name, "products", len(prods), "nonzero incidences", nz, "abs-coef sum", inc,
              "max uses", max(uses.values()))
        vmax = {}
        for A, B in prods:
            for n in list(A) + list(B):
                vmax[n] = P - 1
        for c in outs:
            for k in c:
                if isinstance(k, tuple):
                    vmax[k[1]] = P - 1
        for i, c in enumerate(outs):
            lo, hi = poly_range(expand(prods, c), vmax)
            print(f"  out{i}: lo {lo / P / P:+.2f} p^2  hi {hi / P / P:+.2f} p^2  terms {len(c)}")
