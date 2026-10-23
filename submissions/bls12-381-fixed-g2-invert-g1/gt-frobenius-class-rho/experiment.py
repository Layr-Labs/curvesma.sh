#!/usr/bin/env python3
"""Self-contained checks for the Frobenius-class rho package.

Part A runs on the real BLS12-381 field and the pinned profile basis.  It
checks the tower, the order-12 Frobenius orbit, the twelve-key invariant and
its three-Fp-multiplication construction, the phase convention, exponent
tracking for one walk step, the collision equation, the squaring escape, and
the isomorphism to Fp[w]/(w^12 - 2*w^6 + 2).

Part B is a toy class-rho with distinguished points on q=67 (by default),
ell=1549741, N=(ell-1)/12=129145.  It exercises the DP-terminated multi-walk
loop, the full-element canonical key, the collision equation, and the
rolling exact-state escape.  It is evidence only at toy scope.

Reads a JSON object from stdin: {"order": toy field prime, "trials": int,
"t": float}.  Writes one JSON object to stdout.
"""
import json
import math
import random
import sys
import time

BLS_P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
BLS_R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001


class Field:
    """Fp2 = Fp[u]/(u^2+1), Fp6 = Fp2[v]/(v^3-xi), Fp12 = Fp6[w]/(w^2-v),
    xi = 1+u, p = 7 mod 12.  The tower variable w is the profile basis
    variable of Fp[w]/(w^12 - 2*w^6 + 2), since u = w^6 - 1."""

    def __init__(self, p):
        self.p = p
        if p % 12 != 7:
            raise AssertionError("p must be 7 mod 12")
        self.xi = (1, 1)
        self.d = self.fp2_pow(self.xi, (p - 1) // 6)
        self.d2 = self.fp2_mul(self.d, self.d)
        self.d4 = self.fp2_mul(self.d2, self.d2)
        z = self.fp2_pow(self.d, 1 + p)
        if z[1] != 0:
            raise AssertionError("zeta not in Fp")
        self.zeta = z[0]
        if pow(self.zeta, 6, p) != 1 or pow(self.zeta, 3, p) != p - 1:
            raise AssertionError("zeta is not a primitive sixth root")
        if self.d[1] != (-self.d[0]) % p:
            raise AssertionError("d is not d0*(1-u)")
        if self.fp2_pow(self.d, 6) != ((0, p - 1)):
            raise AssertionError("d^6 != -u")

    # ---------- Fp2 ----------
    def fp2_mul(self, a, b):
        p = self.p
        a0, a1 = a
        b0, b1 = b
        return ((a0 * b0 - a1 * b1) % p, (a0 * b1 + a1 * b0) % p)

    def fp2_sqr(self, a):
        p = self.p
        a0, a1 = a
        return ((a0 * a0 - a1 * a1) % p, (2 * a0 * a1) % p)

    def fp2_add(self, a, b):
        p = self.p
        return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)

    def fp2_sub(self, a, b):
        p = self.p
        return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)

    def fp2_conj(self, a):
        return (a[0], (-a[1]) % self.p)

    def fp2_neg(self, a):
        return ((-a[0]) % self.p, (-a[1]) % self.p)

    def fp2_pow(self, a, e):
        r = (1, 0)
        b = a
        while e:
            if e & 1:
                r = self.fp2_mul(r, b)
            b = self.fp2_sqr(b)
            e >>= 1
        return r

    # ---------- Fp6 ----------
    def fp6_add(self, a, b):
        return tuple(self.fp2_add(x, y) for x, y in zip(a, b))

    def fp6_sub(self, a, b):
        return tuple(self.fp2_sub(x, y) for x, y in zip(a, b))

    def fp6_mul(self, a, b):
        F = self
        a0, a1, a2 = a
        b0, b1, b2 = b
        c0 = F.fp2_mul(a0, b0)
        c1 = F.fp2_mul(a1, b1)
        c2 = F.fp2_mul(a2, b2)
        t12 = F.fp2_mul(F.fp2_add(a1, a2), F.fp2_add(b1, b2))
        t01 = F.fp2_mul(F.fp2_add(a0, a1), F.fp2_add(b0, b1))
        t02 = F.fp2_mul(F.fp2_add(a0, a2), F.fp2_add(b0, b2))
        xi = F.xi
        d0 = F.fp2_add(c0, F.fp2_mul(xi, F.fp2_sub(F.fp2_sub(t12, c1), c2)))
        d1 = F.fp2_add(F.fp2_sub(F.fp2_sub(t01, c0), c1), F.fp2_mul(xi, c2))
        d2 = F.fp2_add(F.fp2_sub(F.fp2_sub(t02, c0), c2), c1)
        return (d0, d1, d2)

    def fp6_mul_v(self, a):
        c0, c1, c2 = a
        return (self.fp2_mul(self.xi, c2), c0, c1)

    def fp6_conj(self, a):
        c0, c1, c2 = a
        F = self
        return (F.fp2_conj(c0), F.fp2_mul(F.fp2_conj(c1), F.d2),
                F.fp2_mul(F.fp2_conj(c2), F.d4))

    def fp6_scalar(self, a, c):
        F = self
        return tuple(F.fp2_mul(x, c) for x in a)

    def fp6_one(self):
        return ((1, 0), (0, 0), (0, 0))

    # ---------- Fp12 ----------
    def fp12_mul(self, a, b):
        F = self
        a0, a1 = a
        b0, b1 = b
        c0 = F.fp6_mul(a0, b0)
        c1 = F.fp6_mul(a1, b1)
        t = F.fp6_mul(F.fp6_add(a0, a1), F.fp6_add(b0, b1))
        d0 = F.fp6_add(c0, F.fp6_mul_v(c1))
        d1 = F.fp6_sub(F.fp6_sub(t, c0), c1)
        return (d0, d1)

    def fp12_sqr(self, a):
        return self.fp12_mul(a, a)

    def fp12_one(self):
        return (self.fp6_one(), ((0, 0), (0, 0), (0, 0)))

    def fp12_is_one(self, a):
        return a == self.fp12_one()

    def fp12_is_zero(self, a):
        return all(x == (0, 0) for x in a[0]) and all(x == (0, 0) for x in a[1])

    def fp12_frob(self, y):
        F = self
        A0, A1 = y
        return (F.fp6_conj(A0), F.fp6_scalar(F.fp6_conj(A1), F.d))

    def fp12_pow(self, y, e):
        r = self.fp12_one()
        b = y
        while e:
            if e & 1:
                r = self.fp12_mul(r, b)
            b = self.fp12_sqr(b)
            e >>= 1
        return r

    def b0(self, y):
        return y[1][0]

    def keys12_formula(self, y):
        """Twelve first-coefficient keys in Frobenius phase order.
        K_{2j}=Re(b0*zeta^j), K_{2j+1}=Re(conj(b0)*d*zeta^j), j=0..5."""
        F = self
        p = F.p
        b0 = F.b0(y)
        zpow = [1]
        for _ in range(5):
            zpow.append(zpow[-1] * F.zeta % p)
        e = F.fp2_mul(F.fp2_conj(b0), F.d)
        ks = []
        for j in range(6):
            ks.append(b0[0] * zpow[j] % p)
            ks.append(e[0] * zpow[j] % p)
        return tuple(ks)

    def keys12_explicit(self, y):
        ks = []
        cur = y
        for _ in range(12):
            ks.append(self.b0(cur)[0])
            cur = self.fp12_frob(cur)
        return tuple(ks)

    def from_parts(self, parts):
        p = self.p
        c = [x % p for x in parts]
        return (((c[0], c[1]), (c[2], c[3]), (c[4], c[5])),
                ((c[6], c[7]), (c[8], c[9]), (c[10], c[11])))

    def fp12_to_coeffs(self, y):
        F = self
        p = F.p
        A0, A1 = y
        cs = [0] * 12
        for idx, c in enumerate(A0):
            x, yy = c
            cs[2 * idx] = (x - yy) % p
            cs[6 + 2 * idx] = (cs[6 + 2 * idx] + yy) % p
        for idx, c in enumerate(A1):
            x, yy = c
            cs[2 * idx + 1] = (x - yy) % p
            cs[7 + 2 * idx] = (cs[7 + 2 * idx] + yy) % p
        return cs

    def fp12_from_coeffs(self, cs):
        F = self
        p = F.p
        xi = F.xi

        def add_xi(c, coef):
            return ((c[0] + coef * xi[0]) % p, (c[1] + coef * xi[1]) % p)

        A0 = (add_xi((cs[0] % p, 0), cs[6]),
              add_xi((cs[2] % p, 0), cs[8]),
              add_xi((cs[4] % p, 0), cs[10]))
        A1 = (add_xi((cs[1] % p, 0), cs[7]),
              add_xi((cs[3] % p, 0), cs[9]),
              add_xi((cs[5] % p, 0), cs[11]))
        return (A0, A1)


def profile_mul(cs, ds, p):
    prod = [0] * 23
    for i in range(12):
        if cs[i]:
            for j in range(12):
                prod[i + j] = (prod[i + j] + cs[i] * ds[j]) % p
    for deg in range(22, 11, -1):
        c = prod[deg]
        if c:
            prod[deg - 6] = (prod[deg - 6] + 2 * c) % p
            prod[deg - 12] = (prod[deg - 12] - 2 * c) % p
    return prod[:12]


def is_probable_prime(n, bases=None):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    for q in small:
        if n % q == 0:
            return n == q
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (bases or small):
        if a >= n:
            continue
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


def production_checks():
    t0 = time.time()
    p, r = BLS_P, BLS_R
    assert p % 12 == 7 and r % 12 == 1
    assert is_probable_prime(r)
    assert pow(p, 12, r) == 1
    assert pow(p, 6, r) == r - 1
    for m in (1, 2, 3, 4, 6):
        assert pow(p, m, r) != 1
    F = Field(p)
    # d properties and the zeta relations from the construction
    assert F.fp2_pow(F.d, 6) == (0, p - 1)
    assert F.fp2_pow(F.d, 1 + p) == (F.zeta, 0)
    assert F.zeta * F.zeta % p == (F.zeta - 1) % p
    assert pow(F.zeta, 3, p) == p - 1
    # tower multiplication matches the profile polynomial basis
    rng = random.Random(20261022)
    y1 = F.from_parts([rng.randrange(p) for _ in range(12)])
    y2 = F.from_parts([rng.randrange(p) for _ in range(12)])
    c1 = F.fp12_to_coeffs(y1)
    c2 = F.fp12_to_coeffs(y2)
    assert F.fp12_from_coeffs(c1) == y1
    assert F.fp12_from_coeffs(c2) == y2
    assert F.fp12_to_coeffs(F.fp12_mul(y1, y2)) == profile_mul(c1, c2, p)
    # a GT generator
    cof = (p ** 12 - 1) // r
    while True:
        y = F.from_parts([rng.randrange(p) for _ in range(12)])
        g = F.fp12_pow(y, cof)
        if not F.fp12_is_one(g) and not F.fp12_is_zero(g):
            break
    assert F.fp12_is_one(F.fp12_pow(g, r))
    assert not F.fp12_is_one(F.fp12_pow(g, 1))
    # the twelve-key formula, phases, and the B1 formula
    y = F.fp12_pow(g, rng.randrange(1, r))
    kf = F.keys12_formula(y)
    ke = F.keys12_explicit(y)
    assert kf == ke
    assert len(set(kf)) == 12
    b0 = F.b0(y)
    e = F.fp2_mul(F.fp2_conj(b0), F.d)
    assert e[0] == (F.d[0] * (b0[0] - b0[1])) % p
    j0 = min(range(12), key=lambda i: kf[i])
    canon0 = tuple(kf[(j0 + i) % 12] for i in range(12))
    cur = y
    for m in range(12):
        ks = F.keys12_formula(cur)
        jm = min(range(12), key=lambda i: ks[i])
        canonm = tuple(ks[(jm + i) % 12] for i in range(12))
        assert canonm == canon0
        assert jm == (j0 - m) % 12
        cur = F.fp12_frob(cur)
    assert cur == y
    # one walk step: phase-corrected multiplier and exponent tracking
    x = rng.randrange(1, r)
    z = F.fp12_pow(g, x)
    a = rng.randrange(r)
    b = rng.randrange(r)
    ye = F.fp12_mul(F.fp12_pow(g, a), F.fp12_pow(z, b))
    alpha = rng.randrange(r)
    beta = rng.randrange(r)
    M = F.fp12_mul(F.fp12_pow(g, alpha), F.fp12_pow(z, beta))
    ks = F.keys12_formula(ye)
    j = min(range(12), key=lambda i: ks[i])
    kap = (-j) % 12
    Mj = M
    for _ in range(kap):
        Mj = F.fp12_frob(Mj)
    ystep = F.fp12_mul(ye, Mj)
    inv = pow(p, kap, r)
    a2 = (a + alpha * inv) % r
    b2 = (b + beta * inv) % r
    lhs = F.fp12_mul(F.fp12_pow(g, a2), F.fp12_pow(z, b2))
    assert lhs == ystep
    # collision equation: y2 = sigma^k(y1) with the same canonical key
    a1 = rng.randrange(r)
    b1 = rng.randrange(r)
    yA = F.fp12_mul(F.fp12_pow(g, a1), F.fp12_pow(z, b1))
    ktrue = rng.randrange(1, 12)
    yB = yA
    for _ in range(ktrue):
        yB = F.fp12_frob(yB)
    ksA = F.keys12_formula(yA)
    ksB = F.keys12_formula(yB)
    jA = min(range(12), key=lambda i: ksA[i])
    jB = min(range(12), key=lambda i: ksB[i])
    assert tuple(ksA[(jA + i) % 12] for i in range(12)) == \
           tuple(ksB[(jB + i) % 12] for i in range(12))
    k = (jA - jB) % 12
    assert k == ktrue
    aB = (pow(p, k, r) * (a1 + x * b1)) % r
    b2c = rng.randrange(r)
    # choose a current record (a,b) with a + x b = q^k (a1 + x b1)
    bcur = rng.randrange(r)
    acur = (aB - x * bcur) % r
    den = (pow(p, k, r) * b1 - bcur) % r
    assert den != 0
    xrec = (acur - pow(p, k, r) * a1) % r * pow(den, r - 2, r) % r
    assert xrec == x
    # the squaring escape preserves y = g^a z^b
    ysq = F.fp12_sqr(ye)
    asq = (2 * a) % r
    bsq = (2 * b) % r
    assert ysq == F.fp12_mul(F.fp12_pow(g, asq), F.fp12_pow(z, bsq))
    return {
        "p_bits": p.bit_length(),
        "r_bits": r.bit_length(),
        "ord_r_p": 12,
        "p6_eq_minus1": True,
        "tower_is_profile_basis": True,
        "keys_formula_equals_frobenius": True,
        "keys_distinct": True,
        "phase_rotation_consistent": True,
        "b1_equals_d0_times_x_minus_y": True,
        "step_exponent_tracking": True,
        "collision_equation_recovers_x": True,
        "squaring_escape_preserves_exponents": True,
        "zeta_squared_eq_zeta_minus_1": True,
        "seconds": round(time.time() - t0, 2),
    }


# ---------- toy class-rho with distinguished points ----------

def factor_phi12(q):
    n = q ** 4 - q ** 2 + 1
    fac = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            fac.append(d)
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        fac.append(n)
    return fac


def toy_run(q, trials, t):
    t0 = time.time()
    if not is_probable_prime(q) or q % 12 != 7:
        raise SystemExit(json.dumps({"error": "order must be a prime 7 mod 12",
                                     "order": q}))
    fac = sorted(set(factor_phi12(q)))
    primes = [f for f in fac if is_probable_prime(f) and f % 12 == 1]
    if not primes:
        raise SystemExit(json.dumps({"error": "no prime 1 mod 12 factor",
                                     "order": q, "factors": fac}))
    ell = max(primes)
    F = Field(q)
    N = (ell - 1) // 12
    sqrtN = math.isqrt(N)
    T = int(t * sqrtN)
    rng = random.Random(0xC0FFEE ^ q)
    g = None
    for _ in range(100):
        y = F.from_parts([rng.randrange(q) for _ in range(12)])
        cand = F.fp12_pow(y, (q ** 12 - 1) // ell)
        if not F.fp12_is_one(cand) and not F.fp12_is_zero(cand):
            g = cand
            break
    if g is None:
        raise SystemExit(json.dumps({"error": "no toy GT generator", "order": q}))
    assert F.fp12_is_one(F.fp12_pow(g, ell))
    R = 64
    theta_mod = 6
    L = 16
    M64 = (1 << 64) - 1

    def mix64(vals):
        v = 0
        for i, val in enumerate(vals):
            v = (v * 0x9E3779B97F4A7C15 + val * 0xBF58476D1CE4E5B9
                 + i * 0x94D049BB133111EB) & M64
            v = (v ^ (v >> 29)) * 0xBF58476D1CE4E5B9 & M64
        return v ^ (v >> 32)

    def key_state(y):
        ks = F.keys12_formula(y)
        j = min(range(12), key=lambda i: ks[i])
        return j, tuple(ks[(j + i) % 12] for i in range(12))

    def canon_elem(y, j):
        for _ in range(j):
            y = F.fp12_frob(y)
        return y

    def multipliers(z):
        conj = [[None] * R for _ in range(12)]
        da = [[0] * R for _ in range(12)]
        db = [[0] * R for _ in range(12)]
        for h in range(R):
            alpha = rng.randrange(ell)
            beta = rng.randrange(ell)
            M = F.fp12_mul(F.fp12_pow(g, alpha), F.fp12_pow(z, beta))
            cur = M
            for k in range(12):
                conj[k][h] = cur
                qk = pow(q, k, ell)
                da[k][h] = alpha * qk % ell
                db[k][h] = beta * qk % ell
                cur = F.fp12_frob(cur)
        return conj, da, db

    ok = 0
    walks_total = 0
    steps_total = 0
    esc_total = 0
    zeroden_total = 0
    false_total = 0
    hits_total = 0
    for _ in range(trials):
        x = rng.randrange(1, ell)
        z = F.fp12_pow(g, x)
        conj, da, db = multipliers(z)
        table = {}
        steps = 0
        success = False
        while steps < T:
            a = rng.randrange(ell)
            b = rng.randrange(ell)
            y = F.fp12_mul(F.fp12_pow(g, a), F.fp12_pow(z, b))
            walks_total += 1
            hist = [None] * L
            while steps < T:
                j, canon = key_state(y)
                hv = mix64(canon[:4])
                if hv % theta_mod == 0:
                    key = canon_elem(y, j)
                    old = table.get(key)
                    if old is None:
                        table[key] = (a, b, j)
                    else:
                        hits_total += 1
                        a1, b1, j1 = old
                        k = (j1 - j) % 12
                        qk = pow(q, k, ell)
                        den = (qk * b1 - b) % ell
                        if den == 0:
                            zeroden_total += 1
                        else:
                            xa = (a - qk * a1) % ell * pow(den, ell - 2, ell) % ell
                            if xa == x:
                                success = True
                            else:
                                false_total += 1
                    break
                h = hv % R
                kap = (-j) % 12
                y = F.fp12_mul(y, conj[kap][h])
                a = (a + da[kap][h]) % ell
                b = (b + db[kap][h]) % ell
                steps += 1
                j2, canon2 = key_state(y)
                fp = (j2, canon2, a, b)
                lag = None
                for li in range(L):
                    if hist[li] == fp:
                        lag = li + 1
                        break
                if lag is not None:
                    y = F.fp12_sqr(y)
                    a = (2 * a) % ell
                    b = (2 * b) % ell
                    esc_total += 1
                    hist = [None] * L
                    continue
                hist = [fp] + hist[:-1]
            if success:
                break
        ok += success
        steps_total += steps
    return {
        "q": q,
        "ell": ell,
        "N": N,
        "sqrtN": sqrtN,
        "T": T,
        "trials": trials,
        "R": R,
        "theta": "1/%d" % theta_mod,
        "success": ok,
        "rate": ok / trials,
        "mean_steps": steps_total / trials,
        "walks": walks_total,
        "escapes": esc_total,
        "zeroden": zeroden_total,
        "false": false_total,
        "dp_hits": hits_total,
        "scope": "toy only; the small q makes the class key and the walk "
                 "less random than the production instance",
        "seconds": round(time.time() - t0, 2),
    }


def main():
    try:
        req = json.load(sys.stdin)
    except Exception:
        req = {}
    order = int(req.get("order", 67))
    trials = int(req.get("trials", 120))
    t = float(req.get("t", 1.05))
    trials = max(1, min(300, trials))
    t = max(0.8, min(1.2, t))
    out = {"production": production_checks()}
    out["toy"] = toy_run(order, trials, t)
    print(json.dumps(out, sort_keys=True))


if __name__ == "__main__":
    main()
