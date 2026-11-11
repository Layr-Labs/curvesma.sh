#!/usr/bin/env python3
"""Exact lattice / norm toolkit for the BN254 tower-NFS sample sets.

Pure standard library, Python 3.9+.  Everything is exact integer arithmetic
and deterministic: the same call gives the same result on every machine.

Elements of Z[y]/h, h = y^6 + y^4 - 2y - 1, are lists of 6 ints, low degree
first.  A cell is a pair (a, b) of such elements.
"""
import hashlib
import json
import math
import os
import time

U = 0x44E992B44A6909F1
P = 36 * U**4 + 36 * U**3 + 24 * U**2 + 6 * U + 1
H = [-1, -2, 0, 0, 1, 0, 1]                    # h, low degree first
BAD_F = 2 * 3 * 197
BAD_G = 2 * 3 * 197 * 13 * 29 * 983 * 11003 * 405928799
BAD_Q = (2, 3, 197, 13, 29, 983, 11003, 405928799)
E2 = [-1, 0, 1, 1, 0, 1]                       # y^5 + y^3 + y^2 - 1
E3 = [-1, 1, 0, 1, 0, 1]                       # y^5 + y^3 + y - 1
Y = [0, 1, 0, 0, 0, 0]
ONE = [1, 0, 0, 0, 0, 0]


# ---------------------------------------------------------------- Z[y]/h
def hmul(a, b, mod=None):
    """Product in Z[y]/h (optionally with coefficients reduced mod `mod`)."""
    c = [0] * 11
    for i in range(6):
        ai = a[i]
        if ai:
            for j in range(6):
                c[i + j] += ai * b[j]
    for k in (10, 9, 8, 7, 6):                 # y^6 = -y^4 + 2y + 1
        t = c[k]
        if t:
            c[k - 2] -= t
            c[k - 5] += 2 * t
            c[k - 6] += t
    if mod is None:
        return c[:6]
    return [x % mod for x in c[:6]]


def det(m):
    """Exact determinant of a square integer matrix (fraction-free)."""
    m = [list(r) for r in m]
    n = len(m)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if m[k][k] == 0:
            for r in range(k + 1, n):
                if m[r][k]:
                    m[k], m[r] = m[r], m[k]
                    sign = -sign
                    break
            else:
                return 0
        pk = m[k][k]
        rk = m[k]
        for i in range(k + 1, n):
            ri = m[i]
            f = ri[k]
            for j in range(k + 1, n):
                ri[j] = (ri[j] * pk - f * rk[j]) // prev
        prev = pk
    return sign * m[n - 1][n - 1]


def hmat(e):
    """Rows e*y^i, i = 0..5 (matrix of multiplication by e)."""
    rows = [list(e)]
    for _ in range(5):
        r = rows[-1]
        t = r[5]                               # shift, then reduce y^6
        rows.append([t, r[0] + 2 * t, r[1], r[2], r[3] - t, r[4]])
    return rows


def hnorm(e):
    """Signed norm N_{K_h/Q}(e)."""
    return det(hmat(e))


def f_elem(a, b):
    """36a^8 + 36y a^6b^2 + 24y^2 a^4b^4 + 6y^3 a^2b^6 + y^4 b^8 in Z[y]/h."""
    A = hmul(a, a)
    B = hmul(b, b)
    yB = hmul(Y, B)
    A2 = hmul(A, A)
    A3 = hmul(A2, A)
    A4 = hmul(A2, A2)
    yB2 = hmul(yB, yB)
    yB3 = hmul(yB2, yB)
    yB4 = hmul(yB2, yB2)
    t1 = hmul(A3, yB)
    t2 = hmul(A2, yB2)
    t3 = hmul(A, yB3)
    return [36 * A4[i] + 36 * t1[i] + 24 * t2[i] + 6 * t3[i] + yB4[i]
            for i in range(6)]


def g_elem(a, b):
    A = hmul(a, a)
    B = hmul(Y, hmul(b, b))
    return [A[i] - U * B[i] for i in range(6)]


def norm_f(a, b):
    return abs(hnorm(f_elem(a, b)))


def norm_g(a, b):
    return abs(hnorm(g_elem(a, b)))


def hinv_unit(e):
    """Inverse of a unit of Z[y]/h (exact; asserts integrality)."""
    from fractions import Fraction
    m = hmat(e)
    # solve x * m = (1,0,0,0,0,0)  <=>  m^T x^T = e0
    aug = [[Fraction(m[r][c]) for r in range(6)] + [Fraction(1 if c == 0 else 0)]
           for c in range(6)]
    for c in range(6):
        piv = next(r for r in range(c, 6) if aug[r][c] != 0)
        aug[c], aug[piv] = aug[piv], aug[c]
        pv = aug[c][c]
        aug[c] = [x / pv for x in aug[c]]
        for r in range(6):
            if r != c and aug[r][c] != 0:
                fct = aug[r][c]
                aug[r] = [x - fct * z for x, z in zip(aug[r], aug[c])]
    x = [aug[r][6] for r in range(6)]
    assert all(v.denominator == 1 for v in x)
    x = [int(v) for v in x]
    assert hmul(x, e) == ONE
    return x


def hpow(e, k):
    if k < 0:
        e = hinv_unit(e)
        k = -k
    r = ONE
    for _ in range(k):
        r = hmul(r, e)
    return r


def unit_table(bound=4):
    """Distinct units y^i e2^j e3^k, |i|,|j|,|k| <= bound, up to sign.

    Returns a list of column tables: cols[c][i] = coefficient c of eps*y^i.
    The identity is first.
    """
    seen = {}
    py = {i: hpow(Y, i) for i in range(-bound, bound + 1)}
    p2 = {i: hpow(E2, i) for i in range(-bound, bound + 1)}
    p3 = {i: hpow(E3, i) for i in range(-bound, bound + 1)}
    order = sorted(range(-bound, bound + 1), key=lambda v: (abs(v), v))
    for i in order:
        for j in order:
            for k in order:
                e = hmul(hmul(py[i], p2[j]), p3[k])
                key = tuple(e)
                neg = tuple(-x for x in e)
                if key in seen or neg in seen:
                    continue
                seen[key] = (i, j, k)
    out = []
    for key in seen:
        rows = hmat(list(key))
        out.append([[rows[i][c] for i in range(6)] for c in range(6)])
    return out


def unit_multiplicity(a, b, A, table):
    """Number of units eps in `table` with eps*a and eps*b both in [-A,A]^6."""
    m = 0
    a0, a1, a2, a3, a4, a5 = a
    b0, b1, b2, b3, b4, b5 = b
    for cols in table:
        ok = True
        for c in cols:
            v = a0 * c[0] + a1 * c[1] + a2 * c[2] + a3 * c[3] + a4 * c[4] + a5 * c[5]
            if v > A or v < -A:
                ok = False
                break
            v = b0 * c[0] + b1 * c[1] + b2 * c[2] + b3 * c[3] + b4 * c[4] + b5 * c[5]
            if v > A or v < -A:
                ok = False
                break
        if ok:
            m += 1
    return m


# ------------------------------------------------------------ small tools
def _sieve(limit):
    s = bytearray([1]) * limit
    s[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(limit) if s[i]]


PRIMES_1000 = _sieve(1000)


def smooth_part(n, primes=PRIMES_1000):
    """Largest divisor of n composed of the given primes (n > 0)."""
    sp = 1
    for l in primes:
        while n % l == 0:
            n //= l
            sp *= l
    return sp


def log2i(n):
    """log2 of a positive integer of any size."""
    return math.log2(n)


def t10(n):
    """round(10*log2 n) for a positive integer."""
    return int(round(10.0 * math.log2(n)))


def admissible(a, b, nf=None, ng=None):
    if nf is None:
        nf = norm_f(a, b)
    if ng is None:
        ng = norm_g(a, b)
    if math.gcd(nf, BAD_F) != 1 or math.gcd(ng, BAD_G) != 1:
        return False
    return math.gcd(hnorm(a), hnorm(b)) == 1


class Rng:
    """Deterministic PRNG: SHA-256 in counter mode over a text seed.

    Block i is sha256(seed + "|" + str(i)); blocks are concatenated big
    endian into a bit stream.  below(n) takes (n-1).bit_length() bits and
    rejects values >= n.
    """

    def __init__(self, seed):
        self.key = str(seed).encode() + b"|"
        self.ctr = 0
        self.buf = 0
        self.nb = 0

    def bits(self, k):
        while self.nb < k:
            d = hashlib.sha256(self.key + str(self.ctr).encode()).digest()
            self.ctr += 1
            self.buf = (self.buf << 256) | int.from_bytes(d, "big")
            self.nb += 256
        self.nb -= k
        r = self.buf >> self.nb
        self.buf &= (1 << self.nb) - 1
        return r

    def below(self, n):
        k = (n - 1).bit_length()
        while True:
            r = self.bits(k)
            if r < n:
                return r

    def rint(self, lo, hi):
        return lo + self.below(hi - lo + 1)


def uniform_cell(rng, A):
    a = [rng.rint(-A, A) for _ in range(6)]
    b = [rng.rint(-A, A) for _ in range(6)]
    return a, b


_MR_BASES = _sieve(180)[:40]


def is_prime(n):
    """Miller-Rabin with the first 40 primes as fixed bases."""
    if n < 2:
        return False
    for l in PRIMES_1000:
        if n == l:
            return True
        if n % l == 0:
            return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def random_prime(rng, nbits):
    """Uniform prime with exactly nbits bits (rejection on odd candidates)."""
    while True:
        c = (1 << (nbits - 1)) | (rng.bits(nbits - 2) << 1) | 1
        if is_prime(c):
            return c


# ---------------------------------------------- polynomials over F_q
def _ptrim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def _pmod(a, m, q):
    """a mod m over F_q, m monic."""
    a = [x % q for x in a]
    dm = len(m) - 1
    for k in range(len(a) - 1, dm - 1, -1):
        t = a[k]
        if t:
            for j in range(dm):
                a[k - dm + j] = (a[k - dm + j] - t * m[j]) % q
    return _ptrim(a[:dm])


def _pmulmod(a, b, m, q):
    if not a or not b:
        return []
    c = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                c[i + j] += ai * bj
    return _pmod(c, m, q)


def _ppow(base, e, m, q):
    r = [1]
    base = _pmod(base, m, q)
    for bit in bin(e)[2:]:
        r = _pmulmod(r, r, m, q)
        if bit == "1":
            r = _pmulmod(r, base, m, q)
    return r


def _pmonic(a, q):
    a = _ptrim([x % q for x in a])
    if not a:
        return a
    inv = pow(a[-1], -1, q)
    return [x * inv % q for x in a]


def _pgcd(a, b, q):
    a = _pmonic(a, q)
    b = _pmonic(b, q)
    while b:
        a, b = b, _pmonic(_pmod(a, b, q), q)
    return a


def _psub(a, b, q):
    n = max(len(a), len(b))
    a = a + [0] * (n - len(a))
    b = b + [0] * (n - len(b))
    return _ptrim([(x - y) % q for x, y in zip(a, b)])


def poly_eval(poly, x, q):
    r = 0
    for c in reversed(poly):
        r = (r * x + c) % q
    return r


def poly_roots(poly, q):
    """Sorted distinct roots in F_q of poly (low degree first), q odd prime.

    Cantor-Zassenhaus; the splitting shifts come from Rng("cz|<q hex>"), so
    the routine is deterministic (and the output, a sorted list, does not
    depend on them anyway).
    """
    f = _pmonic(poly, q)
    if len(f) <= 1:
        return []
    xq = _ppow([0, 1], q, f, q)
    g = _pgcd(_psub(xq, [0, 1], q), f, q)
    rng = Rng("cz|%x" % q)
    out = []
    stack = [g]
    while stack:
        g = stack.pop()
        d = len(g) - 1
        if d <= 0:
            continue
        if d == 1:
            out.append((-g[0]) % q)
            continue
        while True:
            c = rng.below(q)
            w = _ppow([c, 1], (q - 1) // 2, g, q)
            w = _psub(w, [1], q)
            k = _pgcd(w, g, q) if w else g
            if 0 < len(k) - 1 < d:
                break
        # quotient g / k
        quo = _pdiv(g, k, q)
        stack.append(k)
        stack.append(quo)
    return sorted(out)


def _pdiv(a, b, q):
    a = list(a)
    db = len(b) - 1
    out = [0] * (len(a) - db)
    for k in range(len(a) - 1, db - 1, -1):
        t = a[k]
        out[k - db] = t
        if t:
            for j in range(db + 1):
                a[k - db + j] = (a[k - db + j] - t * b[j]) % q
    return out


def _deriv(poly):
    return [i * poly[i] for i in range(1, len(poly))]


def simple_roots(poly, q):
    d = _deriv(poly)
    return [r for r in poly_roots(poly, q) if poly_eval(d, r, q) != 0]


def f_at(s, q):
    """f(s, x) as a polynomial in x over F_q."""
    return [pow(s, 4, q), 0, 6 * pow(s, 3, q) % q, 0, 24 * s * s % q, 0,
            36 * s % q, 0, 36 % q]


def g_at(s, q):
    return [(-U * s) % q, 0, 1]


def degree_one_ideals(q, side):
    """All (s, t): s simple root of h mod q, t simple root of f(s,x) or g(s,x)."""
    out = []
    for s in simple_roots(H, q):
        poly = f_at(s, q) if side == "f" else g_at(s, q)
        for t in simple_roots(poly, q):
            out.append((s, t))
    return out


# ------------------------------------------------------------------ LLL
def lll(basis):
    """Exact integer LLL, delta = 3/4 (Cohen, Algorithm 2.6.7, de Weger form).

    Input: list of linearly independent integer rows.  Returns the reduced
    rows.  Only integer operations are used, so the result is a function of
    the input alone.
    """
    n = len(basis)
    b = [list(r) for r in basis]
    D = [1] + [0] * n                           # D[i+1] = d_i of vector i
    L = [[0] * n for _ in range(n)]
    D[1] = sum(x * x for x in b[0])
    k = 1
    kmax = 0
    while k < n:
        if k > kmax:
            kmax = k
            bk = b[k]
            Lk = L[k]
            for j in range(k + 1):
                bj = b[j]
                u = 0
                for x, z in zip(bk, bj):
                    u += x * z
                Lj = L[j]
                for i in range(j):
                    u = (D[i + 1] * u - Lk[i] * Lj[i]) // D[i]
                if j < k:
                    Lk[j] = u
                else:
                    if u == 0:
                        raise ValueError("rows are dependent")
                    D[k + 1] = u
        while True:
            # RED(k, k-1)
            Lk = L[k]
            l = k - 1
            lam = Lk[l]
            dl = D[k]
            if 2 * abs(lam) > dl:
                r = (2 * lam + dl) // (2 * dl)
                bl = b[l]
                b[k] = [x - r * z for x, z in zip(b[k], bl)]
                Ll = L[l]
                for i in range(l):
                    Lk[i] -= r * Ll[i]
                lam -= r * dl
                Lk[l] = lam
            if 4 * (D[k + 1] * D[k - 1] + lam * lam) < 3 * dl * dl:
                # SWAP(k)
                b[k], b[k - 1] = b[k - 1], b[k]
                Lp = L[k - 1]
                for j in range(k - 1):
                    Lk[j], Lp[j] = Lp[j], Lk[j]
                dk1 = D[k + 1]
                Bn = (D[k - 1] * dk1 + lam * lam) // dl
                for i in range(k + 1, kmax + 1):
                    Li = L[i]
                    t = Li[k]
                    v = (dk1 * Li[k - 1] - lam * t) // dl
                    Li[k] = v
                    Li[k - 1] = (Bn * t + lam * v) // dk1
                D[k] = Bn
                if k > 1:
                    k -= 1
            else:
                for l in range(k - 2, -1, -1):
                    lam = Lk[l]
                    dl = D[l + 1]
                    if 2 * abs(lam) > dl:
                        r = (2 * lam + dl) // (2 * dl)
                        bl = b[l]
                        b[k] = [x - r * z for x, z in zip(b[k], bl)]
                        Ll = L[l]
                        for i in range(l):
                            Lk[i] -= r * Ll[i]
                        Lk[l] = lam - r * dl
                k += 1
                break
    return b


def is_lll_reduced(basis):
    """Independent exact check (rational Gram-Schmidt) of delta = 3/4 LLL."""
    from fractions import Fraction
    n = len(basis)
    bs = []
    mu = [[Fraction(0)] * n for _ in range(n)]
    nb = []
    for i in range(n):
        v = [Fraction(x) for x in basis[i]]
        for j in range(i):
            mu[i][j] = sum(Fraction(x) * z for x, z in zip(basis[i], bs[j])) / nb[j]
            v = [x - mu[i][j] * z for x, z in zip(v, bs[j])]
        bs.append(v)
        nb.append(sum(x * x for x in v))
    for i in range(n):
        for j in range(i):
            if abs(mu[i][j]) > Fraction(1, 2):
                return False
    for i in range(1, n):
        if nb[i] < (Fraction(3, 4) - mu[i][i - 1] ** 2) * nb[i - 1]:
            return False
    return True


# --------------------------------------------------- lattice constructors
def specialq_basis(q, s, t):
    """Row basis of {(a,b) in Z^12 : a(s) - t*b(s) = 0 mod q}; index q.

    Coordinates are (a0..a5, b0..b5).  Row 0 is q*e0; row j >= 1 is
    e_j - c_j*e0 with c_j the coefficient of coordinate j in the linear form.
    """
    c = [pow(s, i, q) for i in range(6)] + [(-t * pow(s, i, q)) % q for i in range(6)]
    rows = [[q] + [0] * 11]
    for j in range(1, 12):
        r = [0] * 12
        r[0] = (-c[j]) % q
        r[j] = 1
        rows.append(r)
    return rows


def specialq_member(q, s, t, v):
    a = poly_eval(v[:6], s, q)
    b = poly_eval(v[6:], s, q)
    return (a - t * b) % q == 0


def specialq_reduced(q, s, t):
    return lll(specialq_basis(q, s, t))


def boot_basis(tau, weight=1):
    """Row basis of {(a,b) : a + tau*b = 0 in F_p[y]/h}, b columns * weight."""
    rows = []
    for i in range(6):
        r = [0] * 12
        r[i] = P
        rows.append(r)
    ty = [x % P for x in tau]
    for i in range(6):
        r = [(-x) % P for x in ty] + [0] * 6
        r[6 + i] = weight
        rows.append(r)
        ty = hmul(ty, Y, P)
    return rows


def boot_member(tau, v):
    tb = hmul([x % P for x in tau], [x % P for x in v[6:]], P)
    return all((v[i] + tb[i]) % P == 0 for i in range(6))


def boot_reduced(tau, wbits=31):
    """LLL of the boot lattice with b weighted by 2^wbits, then un-weighted."""
    w = 1 << wbits
    red = lll(boot_basis(tau, w))
    out = []
    for r in red:
        assert all(x % w == 0 for x in r[6:])
        out.append(r[:6] + [x // w for x in r[6:]])
    return out


def combine(coeffs, basis):
    n = len(basis[0])
    v = [0] * n
    for c, row in zip(coeffs, basis):
        if c:
            for i in range(n):
                v[i] += c * row[i]
    return v


def draw_coeffs(rng, R, n=12):
    while True:
        c = [rng.rint(-R, R) for _ in range(n)]
        if any(c):
            return c


MAX_IDEALS = {"f": 48, "g": 12}    # at most 6 roots of h, times at most 8 (f) or 2 (g) roots above each


def draw_ideal(rng, side, nbits):
    """Uniformly random degree-one prime ideal (q, s, t) of the given side among those above the nbits-bit primes.

    A uniform nbits-bit prime q is drawn and kept with probability (number of degree-one ideals above q) /
    MAX_IDEALS[side]; then one of its ideals is taken uniformly.  Every ideal therefore has the same probability.
    Returns (q, s, t, number of primes tried)."""
    tries = 0
    while True:
        tries += 1
        q = random_prime(rng, nbits)
        if q in BAD_Q or q == P:
            continue
        ids = degree_one_ideals(q, side)
        k = rng.below(MAX_IDEALS[side])
        if k >= len(ids):
            continue
        s, t = ids[k]
        return q, s, t, tries


# ------------------------------------------------------------- packing
_B36 = "0123456789abcdefghijklmnopqrstuvwxyz"


def pack36(vals, width=3):
    """Fixed-width base-36 packing of non-negative ints (< 36^width)."""
    out = []
    lim = 36 ** width
    for v in vals:
        if not 0 <= v < lim:
            raise ValueError("value %r does not fit" % (v,))
        s = ""
        for _ in range(width):
            v, r = divmod(v, 36)
            s = _B36[r] + s
        out.append(s)
    return "".join(out)


def unpack36(s, width=3):
    return [int(s[i:i + width], 36) for i in range(0, len(s), width)]


def pack_cells(recs):
    """recs: list of tuples of ints; each field 3 chars except 0/1 flags.

    Layout is given by the caller through unpack_cells(fmt)."""
    return "".join("".join(_B36[v] if w == 1 else pack36([v], w) for v, w in rec)
                   for rec in recs)


def unpack_cells(s, fmt):
    """fmt: string of field widths, e.g. "33133" -> list of lists of ints."""
    ws = [int(c) for c in fmt]
    step = sum(ws)
    out = []
    for i in range(0, len(s), step):
        rec = []
        p = i
        for w in ws:
            rec.append(int(s[p:p + w], 36))
            p += w
        out.append(rec)
    return out


def hexrow(v):
    return " ".join(("-%x" % -x) if x < 0 else ("%x" % x) for x in v)


def unhexrow(s):
    return [int(x, 16) for x in s.split()]


# --------------------------------------------------------------- replay
# ------------------------------------------------------------ F_(p^12) = F_p[y,x]/(h, x^2 - u y) and GT
def f12_mul(a, b):
    """(a0 + a1 x)(b0 + b1 x) with x^2 = u y; a0, a1, b0, b1 in F_p[y]/h."""
    a0, a1 = a
    b0, b1 = b
    t = hmul(a1, b1, P)
    uy = [(U * v) % P for v in ([0] + t[:5])]          # y * t, before reduction of the y^6 term
    # y * t = t0 y + ... + t5 y^6 and y^6 = -y^4 + 2y + 1
    c = t[5]
    uy[0] = (uy[0] + U * c) % P
    uy[1] = (uy[1] + 2 * U * c) % P
    uy[4] = (uy[4] - U * c) % P
    m0 = hmul(a0, b0, P)
    return ([(x + z) % P for x, z in zip(m0, uy)],
            [(x + z) % P for x, z in zip(hmul(a0, b1, P), hmul(a1, b0, P))])


def f12_pow(a, e):
    r = ([1, 0, 0, 0, 0, 0], [0] * 6)
    while e:
        if e & 1:
            r = f12_mul(r, a)
        a = f12_mul(a, a)
        e >>= 1
    return r


def hpow_mod(a, e):
    r = [1, 0, 0, 0, 0, 0]
    a = [v % P for v in a]
    while e:
        if e & 1:
            r = hmul(r, a, P)
        a = hmul(a, a, P)
        e >>= 1
    return r


R_ORDER = 36 * U**4 + 36 * U**3 + 18 * U**2 + 6 * U + 1


def gt_element(rng):
    """Uniform element of GT (the order-r subgroup of F_(p^12)^*) other than 1: w^((p^12-1)/r) for uniform w."""
    while True:
        w = ([rng.below(P) for _ in range(6)], [rng.below(P) for _ in range(6)])
        t = f12_pow(w, (P**12 - 1) // R_ORDER)
        if t != ([1, 0, 0, 0, 0, 0], [0] * 6):
            return t


def boot_tau(t):
    """Slope tau = t0/t1 of T = t0 + t1 x: a - b x is an F_(p^6)-multiple of T iff a + tau b = 0."""
    t0, t1 = t
    return hmul(t0, hpow_mod(t1, P**6 - 2), P)
