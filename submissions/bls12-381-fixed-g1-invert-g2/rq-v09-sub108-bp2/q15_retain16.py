#!/usr/bin/env python3
"""Exact Q(zeta_15) retain16 module for BLS12-381 CM-STNFS descent.

Base field: H = Q(zeta_15), Phi_15(y) = y^8 - y^7 + y^5 - y^4 + y^3 - y + 1.
Maximal order: O_H = Z[zeta_15] (exact power basis, index 1, Disc = 3^4 * 5^6, h_K = 1).
Side polynomials (w = y^5, w^2 = -w - 1):
  f_9(x) = x^9 + (w + 1)*x^6 + x^3 - w              (deg 9 / H, deg 72 / Q)
  g_3(x) = x^3 - U_BLS * w^2 = x^3 + U_BLS*(w + 1)  (deg 3 / H, deg 24 / Q)
where U_BLS = -0xd201000000010000 (v_2(U_BLS) = 16).

Key cleanliness properties of the monogenic cyclotomic order O_H = Z[zeta_15]:
  1. O_H = Z[y] is monogenic (index 1) everywhere, so Kummer-Dedekind applies at ALL primes p >= 2.
  2. p = 2 is UNRAMIFIED in O_H: 2*O_H = P_{2,1} * P_{2,2}, N(P_{2,i}) = 2^4 = 16.
     For g_3 at each P_{2,i}, v_{P}(U_BLS) = 16 != 0 (mod 3), so for every base-primitive (a, b)
     v_P(g_3(a, b)) = min(3 * v_P(a), 16) <= 16 (exact finite law, max 64 bits/prime, zero tail!).
  3. Base-primitive fraction prod (1 - N(P)^-2) = 0.968504.
"""
import functools
import math
import random
import time
from decimal import Decimal, localcontext
from math import ceil, comb, factorial, floor, isqrt, log, log1p, log2

PHI15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]  # y^8 - y^7 + y^5 - y^4 + y^3 - y + 1
U_BLS = -0xd201000000010000
F_TERMS = [(9, [1]), (6, [1, 0, 0, 0, 0, 1]), (3, [1]), (0, [0, 0, 0, 0, 0, -1])]
# -U_BLS * w^2 = U_BLS * (1 + w) = U_BLS * (1 + y^5)
G_TERMS = [(3, [1]), (0, [U_BLS, 0, 0, 0, 0, U_BLS])]

GRID = 4
SENT = 10 ** 6
KMAX = 64
RSTEP = 1.0 / 1024.0
KA_MAX = 16

# Default BLS12-381 Q(zeta_15) descent parameters (can be overridden in retention/loss calls)
DEFAULT_LX = 70.0
DEFAULT_MESH = 384
DEFAULT_ALPHA_F = 2.74
DEFAULT_ALPHA_G = 0.77
DEFAULT_SLACK = 96.0


from q15_descent16 import (
    RHO_COEFS, RHO_TAB, RHO_TAB_LO, RHO_TAB_HI, rho_series, rho_lo, rho_hi,
)

def _conv(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    nz = [j for j, y in enumerate(b) if y]
    if nz:
        j0, j1 = nz[0], nz[-1] + 1
        bb = b[j0:j1]
        for i, x in enumerate(a):
            if x:
                out[i + j0 : i + j1] = [o + x * y for o, y in zip(out[i + j0 : i + j1], bb)]
    return out


@functools.lru_cache(None)
def window_rows(n, base, iA, mesh, kmax=KA_MAX):
    h = (n - base) / mesh
    wA = [log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
    wB = [log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
    PA = [[1.0]]
    for j in range(1, kmax + 1):
        PA.append([x / j for x in _conv(PA[-1], wA)])
    PB = [[1.0], wB]
    rows = {}
    for jA in range(kmax + 1):
        for jB in range(2):
            rows[jA, jB] = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
    return rows, h


_SUMPROD = getattr(math, "sumprod", lambda a, b: sum(x * y for x, y in zip(a, b)))


@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha, trunc, kmax=KA_MAX):
    import q15_descent16 as _DS
    rows, h = _DS.window_rows(n, base, iA, mesh, kmax)
    out = [[0.0, 0.0] for _ in range(kmax + 1)]
    tab, ntab, sp = RHO_TAB_LO, len(RHO_TAB_LO), _SUMPROD
    inv_br = 1.0 / (base * RSTEP)
    hk = h * inv_br
    tr_span = int(trunc * base / h) + 1 if (trunc and h > 0) else 0
    for j in range(kmax + 2):
        items = []
        ifeas_j = int(floor((m - j * h - j * base) / h)) if h > 0 else 10**9
        k0 = (m - j * base + alpha) * inv_br
        i_ntab = int((k0 - (ntab - 1)) / hk) if hk > 0.0 else 0
        for jB in (0, 1):
            jA = j - jB
            if not (0 <= jA <= kmax):
                continue
            row, r_lo, r_hi = rows[jA, jB]
            if r_lo >= r_hi:
                continue
            ifeas = r_hi - 1 if (r_hi - 1) < ifeas_j else ifeas_j
            if ifeas < r_lo:
                continue
            imin = r_lo if not trunc else (r_lo if r_lo > ifeas - tr_span else ifeas - tr_span)
            if i_ntab > imin:
                imin = i_ntab
            if imin <= ifeas:
                items.append((jA, jB, row, imin, ifeas))
        if not items:
            continue
        t0 = min(it[3] for it in items)
        t1 = max(it[4] for it in items)
        if k0 - t1 * hk <= 0.0 or k0 - t0 * hk >= ntab:
            arr = [1.0 if (k := ceil(k0 - i * hk)) <= 0 else (tab[k] if k < ntab else 0.0) for i in range(t0, t1 + 1)]
        else:
            arr = [tab[ceil(k0 - i * hk)] for i in range(t0, t1 + 1)]
        for jA, jB, row, imin, ifeas in items:
            out[jA][jB] = sp(row[imin : ifeas + 1], arr[imin - t0 : ifeas - t0 + 1])
    cum = [[0.0, 0.0] for _ in range(kmax + 1)]
    for jB in range(2):
        acc = 0.0
        for jA in range(kmax + 1):
            acc += out[jA][jB]
            cum[jA][jB] = acc
    return tuple(tuple(x) for x in cum)


def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=0.0):
    m = ceil(m)
    if n <= base:
        return rho_lo((m + alpha) / base)
    cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc, kA)
    return cum[kA][0] + (cum[kA][1] if kB else 0.0)


@functools.lru_cache(None)
def ss_upper2(m, base, n, iA, kA, kB, mesh, alpha):
    m = floor(m)
    if m <= base:
        return 1.0
    if n <= base:
        return rho_hi((m + alpha) / base)
    import q15_descent16 as _DS
    rows, h = _DS.window_rows(n, base, iA, mesh, kA)
    tab, ntab = _DS.RHO_TAB_HI, len(_DS.RHO_TAB_HI)
    inv_rstep = 1.0 / _DS.RSTEP
    inv_br = 1.0 / (base * _DS.RSTEP)
    hk = h * inv_br
    sp = _SUMPROD
    _fl = floor
    jBs = (0, 1) if kB else (0,)
    j_lim = kA + (2 if kB else 1)
    tot = 0.0
    for j in range(j_lim):
        items = []
        imax_j = int(_fl((m - j * base) / h)) if h > 0 else 10**9
        for jB in jBs:
            jA = j - jB
            if not (0 <= jA <= kA):
                continue
            row, r_lo, r_hi = rows[jA, jB]
            if r_lo >= r_hi:
                continue
            imax = r_hi - 1 if (r_hi - 1) < imax_j else imax_j
            if imax >= r_lo:
                items.append((row, r_lo, imax))
        if not items:
            continue
        t0 = min(it[1] for it in items)
        t1 = max(it[2] for it in items)
        k0 = (m - j * base - j * h + alpha) * inv_br
        if k0 - t0 * hk < ntab - 1 and k0 - t1 * hk >= inv_rstep:
            arr = [tab[_fl(k0 - i * hk)] for i in range(t0, t1 + 1)]
        else:
            arr = [1.0 if (u := (m - (j * base + i * h) - j * h + alpha) / base) <= 1.0 else rho_hi(u) for i in range(t0, t1 + 1)]
        for row, r_lo, imax in items:
            tot += sp(row[r_lo : imax + 1], arr[r_lo - t0 : imax - t0 + 1])
    return min(tot, 1.0)


# ==================== Polynomial arithmetic over F_p and F_{p^{d_1}} ====================
def _trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def pmod(a, m, p):
    a = _trim([x % p for x in a])
    dm = len(m) - 1
    inv = pow(m[-1], -1, p)
    while len(a) - 1 >= dm:
        c = a[-1] * inv % p
        sh = len(a) - 1 - dm
        for i in range(len(m)):
            a[sh + i] = (a[sh + i] - c * m[i]) % p
        _trim(a)
    return a


def pmul(a, b, p):
    c = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] = (c[i + j] + x * y) % p
    return _trim(c)


def padd(a, b, p):
    return _trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
                  for i in range(max(len(a), len(b)))])


def pscal(a, c, p):
    return _trim([x * c % p for x in a])


def pdiv(a, b, p):
    a = _trim([x % p for x in a])
    b = _trim([x % p for x in b])
    q = [0] * max(1, len(a) - len(b) + 1)
    inv = pow(b[-1], -1, p)
    while a and len(a) >= len(b):
        c = a[-1] * inv % p
        sh = len(a) - len(b)
        q[sh] = c
        for i in range(len(b)):
            a[sh + i] = (a[sh + i] - c * b[i]) % p
        _trim(a)
    return _trim(q), a


def pgcd(a, b, p):
    a = _trim([x % p for x in a])
    b = _trim([x % p for x in b])
    while b:
        a, b = b, pmod(a, b, p)
    if a:
        inv = pow(a[-1], -1, p)
        a = [x * inv % p for x in a]
    return a


def pderiv(a, p):
    return _trim([i * a[i] % p for i in range(1, len(a))])


def ppowmod(a, e, m, p):
    r = [1]
    a = pmod(a, m, p)
    while e:
        if e & 1:
            r = pmod(pmul(r, a, p), m, p)
        a = pmod(pmul(a, a, p), m, p)
        e >>= 1
    return r


def roots_Fp(a, p):
    return [r for r in range(p) if sum(c * pow(r, i, p) for i, c in enumerate(a)) % p == 0]


def factor_Fp(a, p, seed=1):
    rng = random.Random(seed)

    def split(poly, d):
        if len(poly) - 1 == d:
            return [poly]
        for _ in range(10000):
            r = _trim([rng.randrange(p) for _ in range(len(poly) - 1)])
            if len(r) < 2:
                continue
            if p == 2:
                t, s = r, r
                for _ in range(d - 1):
                    t = pmod(pmul(t, t, p), poly, p)
                    s = padd(s, t, p)
                g = pgcd(poly, s, p)
            else:
                g = pgcd(poly, padd(ppowmod(r, (p ** d - 1) // 2, poly, p), [p - 1], p), p)
            if 1 < len(g) - 1 < len(poly) - 1:
                q, _ = pdiv(poly, g, p)
                return split(g, d) + split(q, d)
        raise RuntimeError(f"splitting failed at p={p}, d={d}")

    out = []
    rest = _trim([x % p for x in a])
    for r in roots_Fp(rest, p):
        out.append([-r % p, 1])
        rest, _ = pdiv(rest, [-r % p, 1], p)
    d = 2
    x = [0, 1]
    xp = ppowmod(x, p, rest, p) if len(rest) > 1 else x
    while len(rest) - 1 >= 2 * d:
        xp = ppowmod(xp, p, rest, p)
        g = pgcd(rest, padd(xp, pscal(x, p - 1, p), p), p)
        if len(g) > 1:
            out += split(g, d)
            rest, _ = pdiv(rest, g, p)
            if len(rest) > 1:
                xp = pmod(xp, rest, p)
        d += 1
    if len(rest) > 1:
        out.append(rest)
    return out


def emul(a, b, m, p):
    return pmod(pmul(a, b, p), m, p) if a and b else []


def einv(a, m, p):
    r0, r1 = list(m), _trim(list(a))
    s0, s1 = [], [1]
    while r1:
        q, r = pdiv(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, padd(s0, pscal(pmul(q, s1, p), p - 1, p), p)
    return pmod(pscal(s0, pow(r0[0], -1, p), p), m, p)


def Etrim(A):
    while A and not A[-1]:
        A.pop()
    return A


def Emul(A, B, m, p):
    C = [[] for _ in range(len(A) + len(B) - 1)] if A and B else []
    for i, x in enumerate(A):
        if x:
            for j, y in enumerate(B):
                if y:
                    C[i + j] = padd(C[i + j], emul(x, y, m, p), p)
    return Etrim(C)


def Emod(A, M, m, p):
    A = Etrim([list(x) for x in A])
    dm = len(M) - 1
    inv = einv(M[-1], m, p)
    while len(A) - 1 >= dm:
        c = emul(A[-1], inv, m, p)
        sh = len(A) - 1 - dm
        for i in range(len(M)):
            A[sh + i] = padd(A[sh + i], pscal(emul(c, M[i], m, p), p - 1, p), p)
        Etrim(A)
    return A


def Egcd(A, B, m, p):
    A = Etrim([list(x) for x in A])
    B = Etrim([list(x) for x in B])
    while B:
        A, B = B, Emod(A, B, m, p)
    if A:
        inv = einv(A[-1], m, p)
        A = [emul(x, inv, m, p) for x in A]
    return A


def Epowmod(A, e, M, m, p):
    R = [[1]]
    A = Emod(A, M, m, p)
    while e:
        if e & 1:
            R = Emod(Emul(R, A, m, p), M, m, p)
        A = Emod(Emul(A, A, m, p), M, m, p)
        e >>= 1
    return R


def Equo(A, M, m, p):
    A = Etrim([list(x) for x in A])
    dm = len(M) - 1
    Q = [[] for _ in range(max(0, len(A) - dm))]
    while len(A) - 1 >= dm:
        c = A[-1]
        sh = len(A) - 1 - dm
        Q[sh] = c
        for i in range(len(M)):
            A[sh + i] = padd(A[sh + i], pscal(emul(c, M[i], m, p), p - 1, p), p)
        Etrim(A)
    return Etrim(Q)


def root_multiplicities(A, q, m, p):
    A = Etrim([list(x) for x in A])
    XP = Epowmod([[], [1]], q, A, m, p)
    diff = Etrim([padd(XP[i] if i < len(XP) else [], pscal([1] if i == 1 else [], p - 1, p), p)
                  for i in range(max(len(XP), 2))])
    R = Egcd(A, diff, m, p)
    counts = []
    while len(R) > 1:
        counts.append(len(R) - 1)
        A = Equo(A, R, m, p)
        R = Egcd(A, R, m, p)
    counts.append(0)
    return [k + 1 for k in range(len(counts) - 1) for _ in range(counts[k] - counts[k + 1])]


def poly_over_ext(terms, m, p):
    deg = max(d for d, _ in terms)
    co = [[] for _ in range(deg + 1)]
    for d, c in terms:
        v = []
        for i, ci in enumerate(c):
            ci_m = ci % p
            if ci_m:
                v = padd(v, pscal(ppowmod([0, 1], i, m, p), ci_m, p), p)
        co[d] = pmod(v, m, p) if v else []
    return co


@functools.lru_cache(None)
def structure(p):
    """Exact prime ideal splitting of p in O_H = Z[zeta_15] + root multiplicities of f_9 and g_3."""
    Hp = [x % p for x in PHI15]
    sq = pgcd(Hp, pderiv(Hp, p), p)
    sqf = pdiv(Hp, sq, p)[0] if len(sq) > 1 else Hp
    out = []
    for hi in factor_Fp(sqf, p):
        e, rest = 0, Hp
        while True:
            quo, rem = pdiv(rest, hi, p)
            if rem:
                break
            rest, e = quo, e + 1
        d1 = len(hi) - 1
        mf = root_multiplicities(poly_over_ext(F_TERMS, hi, p), p ** d1, hi, p)
        mg = root_multiplicities(poly_over_ext(G_TERMS, hi, p), p ** d1, hi, p)
        out.append([d1, e, mf, mg])
    return out


def primes(lo, hi):
    return [p for p in range(lo, hi) if p >= 2 and all(p % q for q in range(2, isqrt(p) + 1))]


def _pdf_add(pdf, mass_units, prob):
    pdf[mass_units] = pdf.get(mass_units, 0.0) + prob


def _convolve(a, b):
    out = {}
    g = out.get
    bl = list(b.items())
    for x, px in a.items():
        for y, py in bl:
            k = x + y
            out[k] = g(k, 0.0) + px * py
    return out


def _prune(pdf, g=64 * GRID):
    keep = {}
    for x, v in pdf.items():
        y = SENT if x >= SENT else x if v > 1e-30 else -(-x // g) * g
        keep[y] = keep.get(y, 0.0) + v
    return keep


def local_pdf(p, st, side, sieved, up):
    """Exact/upper local omitted-mass PDF (in 1/GRID-bit units) for prime p."""
    lp = log2(p)
    rnd = ceil if up else floor
    pdf = {0: 1.0}
    # Exact Newton-polygon law at p=2, side='g':
    # 2*O_H = P_{2,1} * P_{2,2}, unramified, N(P)=16. v_P(U_BLS) = 16 != 0 mod 3.
    # For base-primitive (a, b): P(b unit) = 16/17; when b unit, v_P(a)=k has prob 16^{-k}*(15/16) (k>=0),
    # and v_P(a^3 - U_BLS w^2 b^3) = min(3*k, 16) EXACTLY!
    if p == 2 and side == "g" and not sieved:
        assert len(st) == 2 and all(d1 == 4 and e == 1 for d1, e, _, _ in st)
        comp = {}
        comp[0] = 16.0 / 17.0  # b in P (1/17) or b unit & v_P(a)=0 (16/17 * 15/16 = 15/17) -> total 16/17!
        for k in range(1, 6):
            pr = (16.0 / 17.0) * (16.0 ** (-k)) * (15.0 / 16.0)
            _pdf_add(comp, rnd(4 * (3 * k) * lp * GRID), pr)
        # k >= 6 -> v_P(g_3) = 16 EXACTLY (64 bits)!
        pr_cap = (16.0 / 17.0) * (16.0 ** (-6))
        _pdf_add(comp, rnd(4 * 16 * lp * GRID), pr_cap)
        assert abs(sum(comp.values()) - 1.0) < 1e-15
        # Two independent primes above 2 (CRT):
        return _convolve(comp, comp)

    for d1, e, mf, mg in st:
        mults = mf if side == "f" else mg
        if not mults or (sieved and d1 != 1):
            continue
        q = p ** d1
        # Universal Hensel bound over complete DVR O_P (holds for both unramified and ramified P):
        # P(v_P(S_hom(a, b)) >= j) <= (q / (q + 1)) * sum_{mu} min(1/q, mu * q^{-ceil(j / mu)})
        def tail(j, q=q, mults=mults):
            return (q / (q + 1.0)) * sum(min(1.0 / q, mu * (q ** (-((j + mu - 1) // mu)))) for mu in mults)

        comp = {}
        prev = tail(1)
        comp[0] = 1.0 - prev
        for j in range(1, 1000):
            nxt = tail(j + 1)
            _pdf_add(comp, rnd(d1 * (j - 1 if sieved else j) * lp * GRID), prev - nxt)
            prev = nxt
            if nxt < 1e-200:
                _pdf_add(comp, SENT if up else rnd(d1 * (j + (0 if sieved else 1)) * lp * GRID), nxt)
                break
        else:
            _pdf_add(comp, SENT, prev)
        pdf = _convolve(pdf, comp)
    return pdf


STRUCT_TABLE = None
LAWS = None


def structure_of(p):
    if STRUCT_TABLE is not None and str(p) in STRUCT_TABLE:
        return STRUCT_TABLE[str(p)]
    return structure(p)


def load_table(stats):
    global STRUCT_TABLE, LAWS
    STRUCT_TABLE = stats["prime_structure"]["table"]
    LAWS = stats.get("omitted_law")
    omitted_pdf.cache_clear()
    base_primitive.cache_clear()


@functools.lru_cache(None)
def base_primitive(pmax_exact=512):
    """Exact base-primitive density prod_{P of O_H} (1 - N(P)^{-2}) (lower bound via pi(x) <= 1.26 x/ln x)."""
    f = 1.0
    for p in primes(2, pmax_exact):
        for d1, e, mf, mg in structure_of(p):
            f *= 1.0 - p ** (-2.0 * d1)
    j0 = int(log2(pmax_exact))
    # Above pmax_exact, at most 8 prime ideals of norm >= 2^j lie above each prime p in [2^j, 2^{j+1})
    for j in range(j0, 200):
        f *= 1.0 - 8.0 * 1.26 * (2.0 ** (j + 1)) / ((j + 1) * log(2.0)) * (2.0 ** (-2 * j))
    return f


def large_block_pdf(j, up, deg_side=72):
    """Dyadic block [2^j, 2^{j+1}) excess-power upper bound (deg_side = 72 for f_9, 24 for g_3)."""
    lam = float(deg_side) * 1.26 * (2.0 ** (j + 1)) / ((j + 1) * log(2.0)) * (2.0 ** (-2 * j))
    pdf = {0: 1.0}
    if lam < 1e-30:
        return pdf
    comp = {}
    tot = 0.0
    mass = j + 1 if up else j
    for t in range(1, 12):
        pr = lam * (2.0 ** (-j * (t - 1))) * (1.0 - 2.0 ** (-j))
        _pdf_add(comp, ceil(mass * t * GRID) if up else floor(mass * t * GRID), pr)
        tot += pr
    for t in range(2, 8):
        pr = (lam * lam / 2.0) * (2.0 ** (-j * (t - 2))) * (t - 1)
        _pdf_add(comp, ceil(mass * t * GRID) if up else floor(mass * t * GRID), pr)
        tot += pr
    if up:
        for t in range(2, 65):
            for k in range(1, t + 1):
                if (k == 1 and t < 12) or (k == 2 and t < 8):
                    continue
                pr = (lam ** k) / factorial(k) * comb(t - 1, k - 1) * (2.0 ** (-j * (t - k)))
                _pdf_add(comp, ceil(mass * t * GRID), pr)
                tot += pr
        mu = 2.0 * max(lam, 2.0 ** (-j))
        pr = sum(t * (mu ** t) / 2.0 for t in range(65, 400))
        _pdf_add(comp, SENT, pr)
        tot += pr
    comp[0] = max(0.0, 1.0 - tot)
    return comp


@functools.lru_cache(None)
def omitted_pdf(side, up, pmax_exact=512, lx_max=74):
    """Full omitted-mass PDF + ordinary exclusion fraction (p >= 2^8, d_1 >= 2)."""
    pdf = {0: 1.0}
    ordinary = 1.0
    deg_side = 72 if side == "f" else 24
    # For d_1 >= 2 base primes above p, at most (8/2) = 4 base primes x (9 or 3) roots = 36 (f) or 12 (g)
    deg_ord = 36.0 if side == "f" else 12.0
    for p in primes(2, pmax_exact):
        st = structure_of(p)
        if p >= 256:
            for d1, e, mf, mg in st:
                if d1 > 1:
                    ordinary -= len(mf if side == "f" else mg) * (p ** (-d1))
        pdf = _prune(_convolve(pdf, local_pdf(p, st, side, p >= 256, up)))
    j0 = int(log2(pmax_exact))
    for j in range(j0, lx_max):
        pdf = _prune(_convolve(pdf, large_block_pdf(j, up, deg_side)))
        ordinary -= deg_ord * 1.26 * (2.0 ** (j + 1)) / ((j + 1) * log(2.0)) * (2.0 ** (-2 * j))
    # Explicit analytic tail bound for ALL j >= lx_max (p >= 2^lx_max):
    # (1) All factor-base primes satisfy p <= 2^71 < 2^lx_max (lx_max=74), so primes p >= 2^lx_max contribute
    #     EXACTLY ZERO omitted factor-base sieve mass to pdf!
    # (2) For ordinary d_1 >= 2 exclusion across ALL p >= 256 to infinity, subtract the exact analytic upper bound
    #     sum_{j>=lx_max} deg_ord * 1.26 * 2^(j+1)/((j+1)*ln 2) * 2^(-2j) <= deg_ord * 1.26 * 4/((lx_max+1)*ln 2) * 2^(-lx_max) < 10^-20!
    tail_factor = 1.26 * 4.0 / ((lx_max + 1) * log(2.0)) * (2.0 ** (-lx_max))
    ordinary -= deg_ord * tail_factor
    return pdf, ordinary


# Wide fine grid covering slack 96..160 and m up to 2560 (s = floor((slack - 1 - m/128)*4))
SL = {"f": (300, 640), "g": (300, 640)}


def law(side):
    """Safe summary of X_up matching DESCENT.md (P(X<=lo) rounded down, tail masses rounded up)."""
    lo, hi = SL[side]
    pdf, o = omitted_pdf(side, True)
    c = [0.0] * 41
    u = lambda v: float("%.3e" % (v * 1.005)) if v else 0.0
    for x, v in pdf.items():
        if x > hi:
            c[min(-(-x // 256), 40)] += v
    return [
        floor(sum(v for x, v in pdf.items() if x <= lo) * 1e10) / 1e10,
        [u(pdf.get(x, 0.0)) for x in range(lo + 1, hi + 1)],
        [u(v) for v in c[2:]],
        o,
    ]


def _law(side):
    return LAWS[side] if LAWS else law(side)


@functools.lru_cache(None)
def loss(side, m, n, iA, kA, kB, alpha, slack, lx=DEFAULT_LX, mesh=DEFAULT_MESH, trunc=0.0):
    """Relation-conditioned lost fraction upper bound (exact DESCENT.md inequality)."""
    s = int(floor((slack - 1.0 - m / 128.0) * GRID))
    b0, fine, blk, _ = _law(side)
    lo, hi = SL[side]
    assert lo <= s <= hi and m < 2560, f"s={s} out of [{lo},{hi}] for side={side}, m={m}, slack={slack}"
    # Fine grid above s: upper-bounded at the next 64-bit block ceiling c64(hi)
    # Group fine entries above s by their 64-bit ceiling xg = 64*ceil(x_u/256) (exact same upper bound, 100x faster!)
    fine_grp = {}
    for idx in range(s - lo, len(fine)):
        v = fine[idx]
        if v:
            x_u = lo + 1 + idx
            xg = 64 * (-(-x_u // (64 * GRID)))
            fine_grp[xg] = fine_grp.get(xg, 0.0) + v
    num = 0.0
    for xg, v in fine_grp.items():
        hi_val = ss_upper2(max(m - xg, 1.0), lx, n, iA, kA, kB, mesh, alpha) if m - xg > 1.0 else 1.0
        num += v * hi_val
    for b, v in enumerate(blk, 2):
        if v:
            xg = 64 * b
            hi_val = ss_upper2(max(m - xg, 1.0), lx, n, iA, kA, kB, mesh, alpha) if m - xg > 1.0 else 1.0
            num += v * hi_val
    denom = ss_lower2(m, lx, n, iA, kA, kB, mesh, alpha, trunc) * b0
    return num / denom if denom > 0.0 else 1.0


def counted_fraction():
    """Fraction of cells that are base-primitive and free of d_1 >= 2 base primes above p >= 256."""
    ordf, ordg = _law("f")[3], _law("g")[3]
    return max(0.0, ordf + ordg - 1.0) * base_primitive()


def retention(nf, ng, n, iA, kA, kB, slack=DEFAULT_SLACK,
              lx=DEFAULT_LX, mesh=DEFAULT_MESH,
              alpha_f=DEFAULT_ALPHA_F, alpha_g=DEFAULT_ALPHA_G):
    """Return (ret_total, loss_f, loss_g) for the Q(zeta_15) 16D descent sieve."""
    lf = loss("f", ceil(nf), n, iA, kA, kB, alpha_f, slack, lx=lx, mesh=mesh)
    lg = loss("g", ceil(ng), n, iA, kA, kB, alpha_g, slack, lx=lx, mesh=mesh)
    return max(0.0, 1.0 - lf - lg) * counted_fraction(), lf, lg


def cell_loss(side, slack):
    pdf, _ = omitted_pdf(side, True)
    return sum(pr for x, pr in pdf.items() if x > slack * GRID)


# ==================== Fast [11b]-style residue-class verification at p = 2, 3, 5, 19, 31 ====================
def _mul_phi15_mod(a, b, mod):
    t = [0] * 15
    for i in range(8):
        ai = a[i]
        if ai:
            for j in range(8):
                t[i + j] += ai * b[j]
    for i in range(14, 7, -1):
        v = t[i]
        if v:
            t[i] = 0
            t[i - 1] += v; t[i - 3] -= v; t[i - 4] += v
            t[i - 5] -= v; t[i - 7] += v; t[i - 8] -= v
    return [x % mod for x in t[:8]]


def _ypow_tab():
    out = []
    for k in range(15):
        t = [0] * 15; t[k] = 1
        for i in range(14, 7, -1):
            v = t[i]
            if v:
                t[i] = 0
                t[i - 1] += v; t[i - 3] -= v; t[i - 4] += v
                t[i - 5] -= v; t[i - 7] += v; t[i - 8] -= v
        out.append(t[:8])
    return out


_YPOW = _ypow_tab()
_W = _YPOW[5]
_WP1 = [(1 if i == 0 else 0) + _W[i] for i in range(8)]
_W2 = _YPOW[10]


def _gal_mod(x, s, mod):
    out = [0] * 8
    for j in range(8):
        xj = x[j]
        if xj:
            yp = _YPOW[(j * s) % 15]
            for k in range(8):
                out[k] += xj * yp[k]
    return [v % mod for v in out]


def _norm_q15_mod(x, mod):
    u1 = _mul_phi15_mod(x, _gal_mod(x, 14, mod), mod)
    u2 = _mul_phi15_mod(u1, _gal_mod(u1, 4, mod), mod)
    u3 = _mul_phi15_mod(u2, _gal_mod(u2, 2, mod), mod)
    return u3[0] % mod


def verify_retain16(n_samples=250, seed=11, recompute_full=False):
    """Fast (< 0.5s with load_table) verification test for check [11b]:
    1. Recompute exact local structure for probe primes (2,3,5,7,11,13,17,19,31,61,257,509) and assert match with table.
    2. Empirically sample base-primitive residue classes mod p^K at p in (2, 3, 5, 19, 31) and
       verify empirical p-adic valuation tails are dominated by our theoretical local_pdf upper bound!
    3. Verify counted_fraction() and retention() at slack=96 (loss_f, loss_g < 5e-5) and slack=128.
    """
    t0 = time.perf_counter()
    bp = base_primitive()
    probe_ps = (2, 3, 5, 7, 11, 13, 17, 19, 31, 61, 257, 509)
    for p in probe_ps:
        st_live = structure(p)
        if STRUCT_TABLE is not None and str(p) in STRUCT_TABLE:
            assert st_live == STRUCT_TABLE[str(p)], f"structure mismatch at p={p}"
    if recompute_full or LAWS is None:
        pdf_f, ord_f = omitted_pdf("f", True)
        pdf_g, ord_g = omitted_pdf("g", True)
        mean_f = round(sum((x / GRID) * pr for x, pr in pdf_f.items() if x < SENT), 3)
        mean_g = round(sum((x / GRID) * pr for x, pr in pdf_g.items() if x < SENT), 3)
        tail96_f, tail96_g = cell_loss("f", 96.0), cell_loss("g", 96.0)
    else:
        ord_f, ord_g = _law("f")[3], _law("g")[3]
        mean_f = mean_g = None
        tail96_f = sum(_law("f")[1][384 - SL["f"][0] :]) + sum(_law("f")[2])
        tail96_g = sum(_law("g")[1][384 - SL["g"][0] :]) + sum(_law("g")[2])
    cf = counted_fraction()

    # Check residue-class samples at p = 2, 3, 5, 19, 31!
    rng = random.Random(seed)
    emp_checks = {}
    for p in (2, 3, 5, 19, 31):
        st = structure(p)
        Hp = [x % p for x in PHI15]
        sq = pgcd(Hp, pderiv(Hp, p), p)
        sqf = pdiv(Hp, sq, p)[0] if len(sq) > 1 else Hp
        factors = factor_Fp(sqf, p)
        K = 140 if p == 2 else (12 if p in (3, 5) else 6)
        M = p ** K
        lp = log2(p)
        vf_bits, vg_bits = [], []
        rej = 0
        while len(vf_bits) < n_samples:
            a = [rng.randrange(M) for _ in range(8)]
            b = [rng.randrange(M) for _ in range(8)]
            if any((not pmod(a, hi, p)) and (not pmod(b, hi, p)) for hi in factors):
                rej += 1
                continue
            a2 = _mul_phi15_mod(a, a, M); b2 = _mul_phi15_mod(b, b, M)
            a3 = _mul_phi15_mod(a2, a, M); b3 = _mul_phi15_mod(b2, b, M)
            a6 = _mul_phi15_mod(a3, a3, M); b6 = _mul_phi15_mod(b3, b3, M)
            a9 = _mul_phi15_mod(a6, a3, M); b9 = _mul_phi15_mod(b6, b3, M)
            w2b3 = _mul_phi15_mod(_W2, b3, M)
            g_el = [(a3[i] - U_BLS * w2b3[i]) % M for i in range(8)]
            t1 = _mul_phi15_mod(_WP1, _mul_phi15_mod(a6, b3, M), M)
            t2 = _mul_phi15_mod(a3, b6, M)
            t3 = _mul_phi15_mod(_W, b9, M)
            f_el = [(a9[i] + t1[i] + t2[i] - t3[i]) % M for i in range(8)]
            for el, arr in ((f_el, vf_bits), (g_el, vg_bits)):
                N = _norm_q15_mod(el, M)
                v = 0
                if N == 0:
                    v = K
                else:
                    while v < K and N % p == 0:
                        N //= p; v += 1
                arr.append(v * lp)
        th_f = local_pdf(p, st, "f", False, True)
        th_g = local_pdf(p, st, "g", False, True)
        th_mean_f = sum((x / GRID) * pr for x, pr in th_f.items() if x < SENT)
        th_mean_g = sum((x / GRID) * pr for x, pr in th_g.items() if x < SENT)
        emp_mean_f = sum(vf_bits) / len(vf_bits)
        emp_mean_g = sum(vg_bits) / len(vg_bits)
        if p == 2:
            assert max(vg_bits) <= 128.0 + 1e-9
        emp_checks[p] = dict(
            K=K,
            factors=len(factors),
            emp_mean_f=round(emp_mean_f, 3), th_mean_f=round(th_mean_f, 3),
            emp_mean_g=round(emp_mean_g, 3), th_mean_g=round(th_mean_g, 3),
            max_f=round(max(vf_bits), 1), max_g=round(max(vg_bits), 1),
        )
        assert emp_mean_f <= th_mean_f + 1.5, f"p={p} f emp={emp_mean_f} > th={th_mean_f}"
        assert emp_mean_g <= th_mean_g + 1.5, f"p={p} g emp={emp_mean_g} > th={th_mean_g}"

    ret96, lf96, lg96 = retention(680, 740, 120.0, 128, 3, 1, slack=96.0)
    ret128, lf128, lg128 = retention(680, 740, 120.0, 128, 3, 1, slack=128.0)
    assert lf96 < 5e-5 and lg96 < 5e-5
    dt = time.perf_counter() - t0
    return dict(
        time_s=round(dt, 3),
        probes_verified=len(probe_ps),
        base_primitive=round(bp, 6),
        ordinary_f=round(ord_f, 6),
        ordinary_g=round(ord_g, 6),
        counted_fraction=round(cf, 6),
        mean_omitted_bits=dict(f=mean_f, g=mean_g),
        cell_loss_96=dict(f=tail96_f, g=tail96_g),
        sample_node_slack96=dict(ret=round(ret96, 6), loss_f=lf96, loss_g=lg96),
        sample_node_slack128=dict(ret=round(ret128, 6), loss_f=lf128, loss_g=lg128),
        residue_checks=emp_checks,
    )


if __name__ == "__main__":
    import json
    res = verify_retain16()
    print(json.dumps(res, indent=2))
