import math, json, functools
import descent16 as D
from common12 import H, U, TERMS
GRID = 4
G_TERMS = [(2, [1]), (0, [0, U])]
DEG_F, DEG_G = 12, 2
BAD = (2, 3, 19, 23, 67)
def _trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a
def pmod(a, m, p):
    a = _trim([x % p for x in a]); dm = len(m) - 1; inv = pow(m[-1], -1, p)
    while len(a) - 1 >= dm:
        c = a[-1] * inv % p; sh = len(a) - 1 - dm
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
    return _trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(max(len(a), len(b)))])
def pscal(a, c, p):
    return _trim([x * c % p for x in a])
def pdiv(a, b, p):
    a = _trim([x % p for x in a]); b = _trim([x % p for x in b]); q = [0] * max(1, len(a) - len(b) + 1); inv = pow(b[-1], -1, p)
    while a and len(a) >= len(b):
        c = a[-1] * inv % p; sh = len(a) - len(b); q[sh] = c
        for i in range(len(b)):
            a[sh + i] = (a[sh + i] - c * b[i]) % p
        _trim(a)
    return _trim(q), a
def pgcd(a, b, p):
    a = _trim([x % p for x in a]); b = _trim([x % p for x in b])
    while b:
        a, b = b, pmod(a, b, p)
    if a:
        inv = pow(a[-1], -1, p); a = [x * inv % p for x in a]
    return a
def pderiv(a, p):
    return _trim([(i * a[i]) % p for i in range(1, len(a))])
def ppowmod(a, e, m, p):
    r = [1]; a = pmod(a, m, p)
    while e:
        if e & 1:
            r = pmod(pmul(r, a, p), m, p)
        a = pmod(pmul(a, a, p), m, p); e >>= 1
    return r
def roots_Fp(a, p):
    return [r for r in range(p) if sum(c * pow(r, i, p) for i, c in enumerate(a)) % p == 0]
def factor_Fp(a, p, seed=1):
    import random
    rng = random.Random(seed)
    def split(a, d):
        if len(a) - 1 == d:
            return [a]
        for _ in range(10000):
            r = _trim([rng.randrange(p) for _ in range(len(a) - 1)])
            if len(r) < 2:
                continue
            if p == 2:
                t = r; s = r
                for _ in range(d - 1):
                    t = pmod(pmul(t, t, p), a, p); s = padd(s, t, p)
                g = pgcd(a, s, p)
            else:
                g = pgcd(a, padd(ppowmod(r, (p ** d - 1) // 2, a, p), [p - 1], p), p)
            if 1 < len(g) - 1 < len(a) - 1:
                q, _ = pdiv(a, g, p); return split(g, d) + split(q, d)
        raise RuntimeError('splitting failed')
    out = []; rest = _trim([x % p for x in a])
    for r in roots_Fp(rest, p):
        out.append([(-r) % p, 1]); rest, _ = pdiv(rest, [(-r) % p, 1], p)
    d = 2; x = [0, 1]
    xp = ppowmod(x, p, rest, p) if len(rest) > 1 else x
    while len(rest) - 1 >= 2 * d:
        xp = ppowmod(xp, p, rest, p)
        g = pgcd(rest, padd(xp, pscal(x, p - 1, p), p), p)
        if len(g) > 1:
            out += split(g, d); rest, _ = pdiv(rest, g, p)
            if len(rest) > 1:
                xp = pmod(xp, rest, p)
        d += 1
    if len(rest) > 1:
        out.append(rest)
    return out
def emul(a, b, m, p):
    return pmod(pmul(a, b, p), m, p) if a and b else []
def einv(a, m, p):
    r0, r1 = list(m), _trim(list(a)); s0, s1 = [], [1]
    while r1:
        q, r = pdiv(r0, r1, p); r0, r1 = r1, r; s0, s1 = s1, padd(s0, pscal(pmul(q, s1, p), p - 1, p), p)
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
    A = Etrim([list(x) for x in A]); dm = len(M) - 1; inv = einv(M[-1], m, p)
    while len(A) - 1 >= dm:
        c = emul(A[-1], inv, m, p); sh = len(A) - 1 - dm
        for i in range(len(M)):
            A[sh + i] = padd(A[sh + i], pscal(emul(c, M[i], m, p), p - 1, p), p)
        Etrim(A)
    return A
def Egcd(A, B, m, p):
    A = Etrim([list(x) for x in A]); B = Etrim([list(x) for x in B])
    while B:
        A, B = B, Emod(A, B, m, p)
    if A:
        inv = einv(A[-1], m, p); A = [emul(x, inv, m, p) for x in A]
    return A
def Ederiv(A, p):
    return Etrim([pscal(A[i], i % p, p) for i in range(1, len(A))])
def Epowmod(A, e, M, m, p):
    R = [[1]]; A = Emod(A, M, m, p)
    while e:
        if e & 1:
            R = Emod(Emul(R, A, m, p), M, m, p)
        A = Emod(Emul(A, A, m, p), M, m, p); e >>= 1
    return R
def degree_one_count(A, q, m, p):
    A = Etrim([list(x) for x in A])
    XP = Epowmod([[], [1]], q, A, m, p)
    diff = Etrim([padd(XP[i] if i < len(XP) else [], pscal([1] if i == 1 else [], p - 1, p), p) for i in range(max(len(XP), 2))])
    return len(Egcd(A, diff, m, p)) - 1
def poly_over_ext(terms, m, p):
    deg = max(d for d, _ in terms); co = [[] for _ in range(deg + 1)]
    for d, c in terms:
        v = []
        for i, ci in enumerate(c):
            if ci:
                v = padd(v, pscal(ppowmod([0, 1], i, m, p), ci % p, p), p)
        co[d] = pmod(v, m, p) if v else []
    return co
@functools.lru_cache(None)
def structure(p):
    Hp = [x % p for x in H]
    if len(pgcd(Hp, pderiv(Hp, p), p)) != 1:
        return None
    out = []
    for h in factor_Fp(Hp, p):
        d1 = len(h) - 1
        F_ = poly_over_ext(TERMS, h, p); G_ = poly_over_ext(G_TERMS, h, p)
        if len(Egcd(F_, Ederiv(F_, p), h, p)) != 1 or len(Egcd(G_, Ederiv(G_, p), h, p)) != 1:
            return None
        out.append(dict(d1=d1, f=[(1, degree_one_count(F_, p ** d1, h, p))], g=[(1, degree_one_count(G_, p ** d1, h, p))]))
    return out
def primes(lo, hi):
    return [p for p in range(lo, hi) if all(p % q for q in range(2, math.isqrt(p) + 1))]
def _pdf_add(pdf, mass_units, prob):
    pdf[mass_units] = pdf.get(mass_units, 0.0) + prob
def good_prime_pdf(p, st, side, sieved, up):
    lp = math.log2(p)
    rnd = math.ceil if up else math.floor
    pdf = {0: 1.0}
    for b in st:
        d1 = b['d1']; n1 = sum(n for dd, n in b[side] if dd == 1)
        if sieved and d1 != 1:
            continue
        qv = p ** d1
        comp = {}; tot = 0.0
        for v in range(1, 120):
            pv = n1 * qv ** (-v) / (1 + 1.0 / qv) * (1 - 1.0 / qv)
            if pv < 1e-40:
                break
            mass = d1 * (v - 1 if sieved else v) * lp
            _pdf_add(comp, rnd(mass * GRID), pv); tot += pv
        _pdf_add(comp, 0, 1 - tot)
        pdf = _convolve(pdf, comp)
    return pdf
def base_factor_degrees(p):
    Hp = [x % p for x in H]
    g = pgcd(Hp, pderiv(Hp, p), p)
    sqf = pdiv(Hp, g, p)[0] if len(g) > 1 else Hp
    return [len(h) - 1 for h in factor_Fp(sqf, p)]
@functools.lru_cache(None)
def base_primitive(pmax_exact=512):
    f = 1.0
    for p in primes(2, pmax_exact):
        for d in base_factor_degrees(p):
            f *= 1 - p ** (-2.0 * d)
    j0 = int(math.log2(pmax_exact))
    for j in range(j0, 200):
        f *= 1 - 6.0 * 1.26 * 2.0 ** (j + 1) / ((j + 1) * math.log(2)) * 2.0 ** (-2 * j)
    return f
def _convolve(a, b):
    out = {}
    for x, px in a.items():
        for y, py in b.items():
            out[x + y] = out.get(x + y, 0.0) + px * py
    return out
def sampled_pdf(p, hist, K, n, side, up):
    lp = math.log2(p); rnd = math.ceil if up else math.floor
    pdf = {}
    for v, c in hist.items():
        v = int(v)
        _pdf_add(pdf, rnd(v * lp * GRID), c / n)
    return pdf
def large_block_pdf(j, up):
    lam = 72.0 * 1.26 * 2.0 ** (j + 1) / ((j + 1) * math.log(2)) * 2.0 ** (-2 * j)
    pdf = {0: 1.0}
    if lam < 1e-30:
        return pdf
    comp = {}; tot = 0.0
    mass = (j + 1) if up else j
    for t in range(1, 12):
        pr = lam * 2.0 ** (-j * (t - 1)) * (1 - 2.0 ** -j)
        _pdf_add(comp, math.ceil(mass * t * GRID) if up else math.floor(mass * t * GRID), pr); tot += pr
    for t in range(2, 8):
        pr = lam * lam / 2 * 2.0 ** (-j * (t - 2)) * (t - 1)
        _pdf_add(comp, math.ceil(mass * t * GRID) if up else math.floor(mass * t * GRID), pr); tot += pr
    comp[0] = max(0.0, 1 - tot)
    return comp
@functools.lru_cache(None)
def omitted_pdf(side, up, pmax_exact=512, bad=None):
    bad = bad or BAD_DATA
    pdf = {0: 1.0}; ordinary = 1.0
    for p in primes(2, pmax_exact):
        if p in BAD:
            d = bad[str(p)]
            comp = sampled_pdf(p, d[side], d['K'], d['n'], side, up)
        else:
            st = structure_of(p)
            if st is None:
                raise ValueError('unexpected non-squarefree prime %d' % p)
            if p >= 256:
                for b in st:
                    d1 = b['d1']; n1 = sum(n for dd, n in b[side] if dd == 1)
                    ordinary -= n1 * p ** (-d1) if d1 > 1 else 0.0
            comp = good_prime_pdf(p, st, side, p >= 256, up)
        pdf = _convolve(pdf, comp)
        pdf = {x: v for x, v in pdf.items() if v > 1e-30}
    j0 = int(math.log2(pmax_exact))
    for j in range(j0, 69):
        pdf = _convolve(pdf, large_block_pdf(j, up))
        pdf = {x: v for x, v in pdf.items() if v > 1e-30}
        ordinary -= (36.0 if side == 'f' else 6.0) * 1.26 * 2.0 ** (j + 1) / ((j + 1) * math.log(2)) * 2.0 ** (-2 * j)
    return pdf, ordinary
BAD_DATA = None
STRUCT_TABLE = None
def load_bad(stats):
    global BAD_DATA, STRUCT_TABLE
    BAD_DATA = stats['bad_prime_valuations']
    STRUCT_TABLE = stats.get('prime_structure', {}).get('table')
    omitted_pdf.cache_clear()


def structure_of(p):
    if STRUCT_TABLE is not None and str(p) in STRUCT_TABLE:
        t = STRUCT_TABLE[str(p)]
        return None if t is None else [dict(d1=d1, f=[(1, nf)], g=[(1, ng)]) for d1, nf, ng in t]
    return structure(p)
@functools.lru_cache(None)
def ss_upper2(m, base, n, iA, kA, kB, mesh, alpha):
    m = math.floor(m)
    if n <= base:
        return D.rho_hi((m + alpha) / base)
    rows, h = D.window_rows(n, base, iA, mesh)
    tot = 0.0
    for (jA, jB), row in rows.items():
        if jA > kA or (jB and not kB):
            continue
        j = jA + jB
        for i, weight in enumerate(row):
            if not weight:
                continue
            left = j * base + i * h
            if left > m:
                continue
            u = (m - left - j * h + alpha) / base
            tot += weight * (1.0 if u <= 1 else D.rho_hi(u))
    return min(tot, 1.0)
@functools.lru_cache(None)
def loss(side, m, n, iA, kA, kB, alpha, slack, step=64):
    Pm = D.P
    s_units = int(math.floor((slack - 1.0 - m / 128.0) * GRID))
    pdf_up, _ = omitted_pdf(side, True)
    base = Pm.LX
    def hi(x_bits):
        xg = math.ceil(x_bits / step) * step
        return ss_upper2(max(m - xg, 1.0), base, n, iA, kA, kB, Pm.MESH, alpha) if m - xg > 1 else 1.0
    def lo(x_bits):
        xg = math.floor(x_bits / step) * step
        return D.ss_lower2(max(m - xg, 1.0), base, n, iA, kA, kB, Pm.MESH, alpha, Pm.TRUNC) if m - xg > 1 else 1.0
    num = sum(pr * hi(x / GRID) for x, pr in pdf_up.items() if x > s_units)
    den = sum(pr * lo(x / GRID) for x, pr in pdf_up.items() if x <= s_units)
    return num / den if den > 0 else 1.0
def counted_fraction():
    _, ordf = omitted_pdf('f', True); _, ordg = omitted_pdf('g', True)
    return max(0.0, ordf + ordg - 1) * base_primitive()


def retention(nf, ng, n, iA, kA, kB, slack):
    lf = loss('f', math.ceil(nf), n, iA, kA, kB, D.P.ALPHA_F, slack)
    lg = loss('g', math.ceil(ng), n, iA, kA, kB, D.P.ALPHA_G, slack)
    _, ordf = omitted_pdf('f', True); _, ordg = omitted_pdf('g', True)
    return max(0.0, (1 - lf - lg)) * max(0.0, ordf + ordg - 1) * base_primitive(), lf, lg
def cell_loss(side, slack):
    pdf, _ = omitted_pdf(side, True)
    return sum(pr for x, pr in pdf.items() if x > slack * GRID)
