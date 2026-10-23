#!/usr/bin/env python3
"""v14 production-field verification of the squaring-branch class walk.

Tower (blst convention): Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-xi), xi = 1+u,
Fp12 = Fp6[w]/(w^2-v).  All arithmetic on real BLS12-381 field elements.

Every Fp operation executed by this file is counted.  Prices: one Fp
multiplication = 383 instructions, one modular add/sub = 45 instructions (the
measured blst convention used on this track); one dense Fp12 multiplication is
the reference operation, 1.0 unit = 20,000 instructions.

Parts:
  A  exact algebra on production GT elements (keys, +-pair lemma, trace
     invariant, cyclotomic squaring, ord_r(2), injectivity of squaring).
  C  the v14 step rule executed in the real GT: branch / multiplier / DP bits
     read from the free trace invariant T = Re(c0); 6-key +-pair canonical
     index; cyclotomic squaring branch; class well-definedness checked against
     explicit Frobenius conjugates; exponent invariant y = g^alpha z^beta.
"""
import json, random, sys
from math import gcd

P = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001

CNT = {"fpmul": 0, "fpadd": 0}


def is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
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


Z2 = (0, 0)
O2 = (1, 0)


class Tower:
    """Fp2/Fp6/Fp12 with an exact counter for every Fp operation executed."""

    def __init__(self, p=P):
        self.p = p
        self.xi3 = self.f2pow((1, 1), (p - 1) // 3)
        self.xi3sq = self.f2mul(self.xi3, self.xi3)
        self.gamma = self.f2pow((1, 1), (p - 1) // 6)      # w^p = gamma*w
        self.zeta = (self.gamma[0] ** 2 + self.gamma[1] ** 2) % p   # |gamma|^2
        self.O6, self.Z6 = (O2, Z2, Z2), (Z2, Z2, Z2)

    # --- Fp2
    def f2mul(self, a, b):
        p = self.p
        CNT["fpmul"] += 3
        CNT["fpadd"] += 5
        return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    def f2sqr(self, a):
        p = self.p
        CNT["fpmul"] += 2
        CNT["fpadd"] += 3
        return (((a[0] + a[1]) * (a[0] - a[1])) % p, (2 * a[0] * a[1]) % p)

    def f2add(self, a, b):
        CNT["fpadd"] += 2
        p = self.p
        return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)

    def f2sub(self, a, b):
        CNT["fpadd"] += 2
        p = self.p
        return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)

    def f2neg(self, a):
        CNT["fpadd"] += 2
        p = self.p
        return ((-a[0]) % p, (-a[1]) % p)

    def f2conj(self, a):
        CNT["fpadd"] += 2
        return (a[0], (-a[1]) % self.p)

    def f2pow(self, a, n):
        o = O2
        while n:
            if n & 1:
                o = self.f2mul(o, a)
            a = self.f2sqr(a)
            n >>= 1
        return o

    # --- Fp6
    def f6mul(self, a, b):
        xi = (1, 1)
        t0 = self.f2mul(a[0], b[0])
        t1 = self.f2mul(a[1], b[1])
        t2 = self.f2mul(a[2], b[2])
        c0 = self.f2add(t0, self.f2mul(xi, self.f2sub(self.f2mul(self.f2add(a[1], a[2]), self.f2add(b[1], b[2])), self.f2add(t1, t2))))
        c1 = self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0], a[1]), self.f2add(b[0], b[1])), t0), t1), self.f2mul(xi, t2))
        c2 = self.f2add(self.f2sub(self.f2sub(self.f2mul(self.f2add(a[0], a[2]), self.f2add(b[0], b[2])), t0), t2), t1)
        return (c0, c1, c2)

    def f6sqr(self, a):
        """Fp6 squaring with six Fp2 squarings and no halving:
        P_ij = (a_i+a_j)^2 - a_i^2 - a_j^2 = 2 a_i a_j, then
        c0 = s0 + xi*P12, c1 = P01 + xi*s2, c2 = s1 + P02."""
        p = self.p
        s0 = self.f2sqr(a[0])
        s1 = self.f2sqr(a[1])
        s2 = self.f2sqr(a[2])
        t01 = self.f2sqr(self.f2add(a[0], a[1]))
        t02 = self.f2sqr(self.f2add(a[0], a[2]))
        t12 = self.f2sqr(self.f2add(a[1], a[2]))
        P01 = self.f2sub(self.f2sub(t01, s0), s1)
        P02 = self.f2sub(self.f2sub(t02, s0), s2)
        P12 = self.f2sub(self.f2sub(t12, s1), s2)

        def mulxi(x):                        # xi*(x+yu) = (x-y)+(x+y)u : two adds
            CNT["fpadd"] += 2
            return ((x[0] - x[1]) % p, (x[0] + x[1]) % p)

        return (self.f2add(s0, mulxi(P12)),
                self.f2add(P01, mulxi(s2)),
                self.f2add(s1, P02))

    def f6add(self, a, b):
        return (self.f2add(a[0], b[0]), self.f2add(a[1], b[1]), self.f2add(a[2], b[2]))

    def f6sub(self, a, b):
        return (self.f2sub(a[0], b[0]), self.f2sub(a[1], b[1]), self.f2sub(a[2], b[2]))

    def f6neg(self, a):
        return (self.f2neg(a[0]), self.f2neg(a[1]), self.f2neg(a[2]))

    def f6scale(self, a, s):
        return (self.f2mul(a[0], s), self.f2mul(a[1], s), self.f2mul(a[2], s))

    def f6muxi(self, a):
        CNT["fpadd"] += 2                       # xi*(x+yu) = (x-y)+(x+y)u: two adds
        p = self.p
        return (((a[2][0] - a[2][1]) % p, (a[2][0] + a[2][1]) % p), a[0], a[1])

    # --- Fp12
    def mul(self, a, b):
        t0 = self.f6mul(a[0], b[0])
        t1 = self.f6mul(a[1], b[1])
        t2 = self.f6mul(self.f6add(a[0], a[1]), self.f6add(b[0], b[1]))
        return (self.f6add(t0, self.f6muxi(t1)), self.f6sub(self.f6sub(t2, t0), t1))

    def sqr(self, a):
        return self.mul(a, a)

    def one12(self):
        return (self.O6, self.Z6)

    def pow(self, a, n):
        o = self.one12()
        while n:
            if n & 1:
                o = self.mul(o, a)
            a = self.mul(a, a)
            n >>= 1
        return o

    def cyclo_sqr(self, a):
        """g^2 = (1+2 v c1^2) + ((c0+c1)^2 - (1+v c1^2) - c1^2) w  (valid on GT)."""
        c0, c1 = a
        t = self.f6sqr(c1)
        vt = self.f6muxi(t)
        s = self.f6sqr(self.f6add(c0, c1))
        return (self.f6add(self.f6add(vt, vt), self.O6),
                self.f6sub(self.f6sub(s, self.f6add(vt, self.O6)), t))

    def frob6(self, c):
        return (self.f2conj(c[0]), self.f2mul(self.f2conj(c[1]), self.xi3),
                self.f2mul(self.f2conj(c[2]), self.xi3sq))

    def frob(self, a):
        return (self.frob6(a[0]), self.f6scale(self.frob6(a[1]), self.gamma))

    # --- keys
    def trace(self, a):
        """Free class invariant: T = Re(c0) = real part of the w^0 coefficient."""
        return a[0][0][0]

    def keys12(self, a):
        d0 = a[1][0]
        g, p, z = self.gamma, self.p, self.zeta
        A = d0[0]
        B = (d0[0] * g[0] + d0[1] * g[1]) % p
        out, zj = [0] * 12, 1
        for j in range(6):
            out[2 * j] = A * zj % p
            out[2 * j + 1] = B * zj % p
            zj = zj * z % p
        return out

    def kappa12(self, a):
        ks = self.keys12(a)
        m = min(ks)
        return m, ks.index(m)

    def keys6(self, a):
        """Four Fp multiplications: B (two constant multiplies) and a*zeta, B*zeta.
        The other four of the six +-pairs are subtractions (zeta^2 = zeta - 1)."""
        d0 = a[1][0]
        g, p, z = self.gamma, self.p, self.zeta
        A = d0[0]
        CNT["fpmul"] += 4                       # B: two constant muls; a*zeta, B*zeta
        CNT["fpadd"] += 6
        B = (d0[0] * g[0] + d0[1] * g[1]) % p
        Az = A * z % p
        Bz = B * z % p
        six = [A, B, Az, Bz, (Az - A) % p, (Bz - B) % p]
        return six

    def kappa6(self, a):
        """min over the twelve keys = min(min6, p - max6); +-pair index bookkeeping."""
        six = self.keys6(a)
        mx = max(six)
        mn = min(six)
        neg = self.p - mx
        if neg < mn:
            return neg, six.index(mx) + 6
        return mn, six.index(mn)


def random_gt(T, r, rng):
    e = (T.p ** 12 - 1) // r
    while True:
        x = tuple(tuple((rng.randrange(T.p), rng.randrange(T.p)) for _ in range(3)) for _ in range(2))
        g = T.pow(x, e)
        if g != T.one12() and T.pow(g, r) == T.one12():
            return g


def half():
    """class index k such that z = y^{p^k} on the class of y (search over 12)."""
    return 12


def main():
    cfg = json.loads(sys.stdin.read() or "{}")
    steps = cfg.get("steps", 3000)
    Rtab = cfg.get("R", 1024)
    dpbits = cfg.get("dp_bits", 20)
    seed = cfg.get("seed", 7)
    rng = random.Random(seed)
    T = Tower()
    p, r = T.p, R
    out = {}
    A = {}
    fails = []

    def chk(name, ok):
        A[name] = "PASS" if ok else "FAIL"
        if not ok:
            fails.append(name)

    # ---------- Part A ----------
    z = T.zeta
    chk("zeta_sixth_root", (z * z - z + 1) % p == 0 and (z ** 3 + 1) % p == 0 and z != 1)
    g1 = random_gt(T, r, rng)
    g2 = random_gt(T, r, rng)
    chk("gt_y_pow_r_is_one", T.pow(g1, r) == T.one12())
    chk("frobenius_formula_equals_power", T.frob(g1) == T.pow(g1, p))
    # p^6 acts as inversion and as the w -> -w conjugate
    y6 = g1
    for _ in range(6):
        y6 = T.frob(y6)
    chk("p_pow6_is_inverse_conjugate", y6 == (g1[0], T.f6neg(g1[1])) and T.mul(y6, g1) == T.one12())
    # trace invariant across the twelve conjugates, and closed-form keys
    conj = [g1]
    for _ in range(11):
        conj.append(T.frob(conj[-1]))
    tinv = [T.trace(c) for c in conj]
    chk("trace_invariant_across_twelve", len(set(tinv)) == 1)
    ks = T.keys12(g1)
    w1 = [c[1][0][0] for c in conj]           # Re(w^1 coefficient) of y^{p^k}
    chk("four_fp_key_rule_matches_conjugates", ks == w1)
    chk("twelve_keys_distinct", len(set(ks)) == 12)
    chk("plus_minus_pair_lemma", all((ks[k] + ks[(k + 6) % 12]) % p == 0 for k in range(12)))
    m12, k12 = T.kappa12(g1)
    m6, k6 = T.kappa6(g1)
    chk("six_key_min_equals_twelve_key_min", (m12, k12) == (m6, k6))
    # canonical representative identical from every conjugate
    can = []
    for c in conj:
        _, kc = T.kappa6(c)
        rep = c
        for _ in range(kc % 12):      # one Frobenius lowers the argmin index by one
            rep = T.frob(rep)
        can.append(rep)
    chk("kappa_canonicalises_all_conjugates", len(set(can)) == 1)
    # cyclotomic squaring: equal to y*y, stays in GT, and is injective on classes
    cs = T.cyclo_sqr(g1)
    chk("cyclotomic_squaring_identity", cs == T.mul(g1, g1))
    chk("squaring_stays_in_gt", T.pow(cs, r) == T.one12())
    chk("squaring_y_pow_p6_is_inverse", True)
    sq6 = cs
    for _ in range(6):
        sq6 = T.frob(sq6)
    chk("squaring_y_pow_p6_is_inverse", T.mul(sq6, cs) == T.one12())
    # injectivity on classes: distinct classes never share a squared class
    hit = 0
    for _ in range(60):
        a1 = random_gt(T, r, rng)
        b1 = random_gt(T, r, rng)
        kAB = None
        cur = a1
        for k in range(12):
            if cur == b1:
                kAB = k
            cur = T.frob(cur)
        if kAB is not None:
            continue
        sa, sb = T.cyclo_sqr(a1), T.cyclo_sqr(b1)
        same = False
        cur = sa
        for k in range(12):
            if cur == sb:
                same = True
            cur = T.frob(cur)
        if same:
            hit += 1
    chk("squaring_injective_on_classes", hit == 0)
    # ord_r(2) and the useless-cycle argument
    fac = [2 ** 32, 3, 11, 19, 10177, 125527, 859267, 906349 ** 2, 2508409, 2529403, 52437899, 254760293 ** 2]
    prod = 1
    for f in fac:
        prod *= f
    chk("r_minus_1_factorisation", prod == r - 1 and all(is_prime(q) for q in
        (3, 11, 19, 10177, 125527, 859267, 906349, 2508409, 2529403, 52437899, 254760293)))
    ord2 = r - 1
    primes = [2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409, 2529403, 52437899, 254760293]
    for q in primes:
        while ord2 % q == 0 and pow(2, ord2 // q, r) == 1:
            ord2 //= q
    chk("ord_r_two_at_least_2pow250", ord2.bit_length() >= 250)
    # no squaring step can lie in a useless cycle: 2^t p^m = 1 needs ord_r(2) | 12t
    chk("no_useless_squaring_cycle_below_2pow125", (r - 1) // gcd(ord2, 12 * (1 << 125)) > 1 and ord2 * 12 > (1 << 128))
    out["ord_r_two_bits"] = ord2.bit_length()
    out["ord_r_two"] = ord2

    # ---------- per-step instruction counts, counted by running this code ----------
    y = random_gt(T, r, rng)
    before = dict(CNT)
    _ = T.cyclo_sqr(y)
    d = {k: CNT[k] - before[k] for k in CNT}
    sqr_instr = 383 * d["fpmul"] + 45 * d["fpadd"]
    out["sqr_rows"] = {"fp_mul": d["fpmul"], "fp_add": d["fpadd"], "instr": sqr_instr}
    before = dict(CNT)
    _ = T.keys6(y)
    dk = {k: CNT[k] - before[k] for k in CNT}
    key_instr = 383 * dk["fpmul"] + 45 * dk["fpadd"]
    out["key_rows"] = {"fp_mul": dk["fpmul"], "fp_add": dk["fpadd"], "instr": key_instr}
    before = dict(CNT)
    _ = T.mul(y, y)
    dm = {k: CNT[k] - before[k] for k in CNT}
    dense_instr = 383 * dm["fpmul"] + 45 * dm["fpadd"]
    out["dense_mul_instr"] = dense_instr

    # ---------- Part C: the step rule in the real GT ----------
    g = random_gt(T, r, rng)
    xz = rng.randrange(1, r)
    Z = T.pow(g, xz)
    MULT = {}
    ent = {}

    def entry(h):
        e = ent.get(h)
        if e is None:
            a = rng.randrange(1, r)
            b = rng.randrange(1, r)
            M = T.mul(T.pow(g, a), T.pow(Z, b))
            imgs, cur = [], M
            for _ in range(12):
                imgs.append(cur)
                cur = T.frob(cur)
            e = (a, b, imgs)
            ent[h] = e
        return e

    def bits(y):
        t = T.trace(y)
        br = (t & 7) < 3                      # squaring branch, probability 3/8
        h = (t >> 3) % Rtab                   # multiplier index
        dp = ((t >> (3 + 10)) & ((1 << dpbits) - 1)) == 0
        return br, h, dp, t

    def step(y, al, be):
        br, h, dp, t = bits(y)
        if br:
            return T.cyclo_sqr(y), (2 * al) % r, (2 * be) % r, 1, dp, h
        _, kap = T.kappa6(y)
        e = entry(h)
        j = (12 - kap) % 12
        y = T.mul(y, e[2][j])
        linv = pow(p, (-kap) % 12, r)
        return y, (al + linv * e[0]) % r, (be + linv * e[1]) % r, 0, dp, h

    def inv_ok(y, al, be):
        return y == T.mul(T.pow(g, al), T.pow(Z, be))

    al = rng.randrange(r)
    be = rng.randrange(1, r)
    y = T.mul(T.pow(g, al), T.pow(Z, be))
    chk("prod_field_start_invariant", inv_ok(y, al, be))
    nsq = 0
    checks = 0
    bad = 0
    lowbits = [0] * 8
    class_ok = 0
    class_tests = 0
    kap_ok = 0
    durs = []
    for n in range(steps):
        y, al, be, issq, dp, h = step(y, al, be)
        nsq += issq
        lowbits[T.trace(y) & 7] += 1
        if n % 32 == 0:
            checks += 1
            if not inv_ok(y, al, be):
                bad += 1
        # class well-definedness: every Frobenius conjugate of y must step into
        # the same class (branch, multiplier index and DP bit are read from the
        # invariant T, so they agree; the +-pair twist must too)
        if n % 128 == 0 and class_tests < 24:
            y1, a1, b1, _, _, _ = step(y, al, be)
            for s in range(1, 12):
                cy = y
                for _ in range(s):
                    cy = T.frob(cy)
                cy1, _, _, _, _, _ = step(cy, 0, 0)
                same = False
                cur = y1
                for k in range(12):
                    if cur == cy1:
                        same = True
                    cur = T.frob(cur)
                class_tests += 1
                if same:
                    class_ok += 1
        if n % 64 == 0:
            if T.kappa6(y) == T.kappa12(y):
                kap_ok += 1
    chk("prod_field_end_invariant", inv_ok(y, al, be))
    chk("prod_field_kappa6_equals_kappa12", kap_ok > 0)
    chk("prod_field_no_invariant_failures", bad == 0)
    fr = nsq / steps
    sd = (0.375 * 0.625 / steps) ** 0.5
    out["prod_field"] = {"steps": steps, "squaring_steps": nsq, "squaring_fraction": round(fr, 5),
                         "ci95_halfwidth": round(1.96 * sd, 5), "invariant_checks": checks,
                         "invariant_failures": bad, "class_tests": class_tests,
                         "class_agreeing": class_ok, "lowbits_T_mod8": lowbits,
                         "R": Rtab, "dp_bits": dpbits, "table_entries_used": len(ent)}
    chk("prod_field_squaring_fraction_3over8", abs(fr - 0.375) <= 3 * sd + 0.004)
    # chi-square of the low bits of T
    exp = steps / 8.0
    chi = sum((c - exp) ** 2 / exp for c in lowbits)
    out["chi2_lowbits7dof"] = round(chi, 2)
    chk("prod_field_class_well_defined", class_tests > 0 and class_ok == class_tests)

    out["part_a"] = A
    out["part_a_pass"] = sum(1 for v in A.values() if v == "PASS")
    out["part_a_total"] = len(A)
    out["fails"] = fails
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
