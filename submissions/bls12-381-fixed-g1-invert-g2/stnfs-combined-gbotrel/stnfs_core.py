"""Core algebraic, Murphy-alpha, Rouché interval, canonical-a, and main-phase
relation-collection routines for the combined STNFS submission."""

import bisect
import functools
import itertools
import math
import random
from array import array
from collections import Counter
from decimal import Decimal as D, localcontext
from fractions import Fraction as F

LN2 = math.log(2.0)
UNIT = 20000.0

P_MOD = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R_ORD = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
U_ABS = 0xD201000000010000
LOG2U6 = 6.0 * math.log2(U_ABS)

H_POLY = [-1, -1, 1, 1, 0, 0, 1]
H_NZ = [(j, H_POLY[j]) for j in range(6) if H_POLY[j]]
F_TERMS = [(12, [1]), (10, [0, -2]), (6, [0, 0, 0, 2]), (2, [0, 0, 0, 0, 0, 1]), (0, [1, 1, -1, -1])]
G_TERMS = [(2, [1]), (0, [0, U_ABS])]

MATRIX_A = [
    ["0.5528759885", "-0.1260799858", "0.1447295864", "-0.002614363166", "-0.02018451304", "0.0956632886"],
    ["-0.1260799858", "0.4730025995", "-0.06826093202", "0.08498958762", "-0.001395387797", "-0.06294781792"],
    ["0.1447295864", "-0.06826093202", "0.4552066288", "-0.05662522189", "0.04120452268", "-0.02031101183"],
    ["-0.002614363166", "0.08498958762", "-0.05662522189", "0.4367157553", "0.0001512689972", "-0.02187465135"],
    ["-0.02018451304", "-0.001395387797", "0.04120452268", "0.0001512689972", "0.4322681594", "0.02026484702"],
    ["0.0956632886", "-0.06294781792", "-0.02031101183", "-0.02187465135", "0.02026484702", "0.4746279595"],
]
MATRIX_B = [
    ["0.5907186542", "-0.1192505666", "0.1631112482", "-0.003050921392", "-0.04000276543", "0.1352719521"],
    ["-0.1192505666", "0.4974396776", "-0.08398310872", "0.09084451399", "-0.01109118934", "-0.08908672312"],
    ["0.1631112482", "-0.08398310872", "0.4933090946", "-0.0460977154", "0.03924332861", "-0.03828634712"],
    ["-0.003050921392", "0.09084451399", "-0.0460977154", "0.4744115241", "0.02143547529", "-0.04387611416"],
    ["-0.04000276543", "-0.01109118934", "0.03924332861", "0.02143547529", "0.4924220451", "0.05071856806"],
    ["0.1352719521", "-0.08908672312", "-0.03828634712", "-0.04387611416", "0.05071856806", "0.5676976811"],
]


def zy_mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] += x * y
    for i in range(len(c) - 1, 5, -1):
        t = c[i]
        if t:
            for j, hj in H_NZ:
                c[i - 6 + j] -= t * hj
    return (c + [0] * 6)[:6]


def bareiss_det(m_in):
    m = [row[:] for row in m_in]
    n = len(m)
    prev, sign = 1, 1
    for k in range(n - 1):
        if m[k][k] == 0:
            row = next((j for j in range(k + 1, n) if m[j][k]), None)
            if row is None:
                return 0
            m[k], m[row] = m[row], m[k]
            sign = -sign
        piv = m[k][k]
        for i in range(k + 1, n):
            mik = m[i][k]
            for j in range(k + 1, n):
                val = m[i][j] * piv - mik * m[k][j]
                if val % prev:
                    raise AssertionError("inexact Bareiss division")
                m[i][j] = val // prev
        prev = piv
    return sign * m[-1][-1] if n else 1


def zy_norm(a):
    cols = [zy_mul(a, [0] * i + [1]) for i in range(6)]
    m = [[cols[j][i] for j in range(6)] for i in range(6)]
    return abs(bareiss_det(m))


def exact_norms(a, b):
    ap, bp = [[1]], [[1]]
    for _ in range(12):
        ap.append(zy_mul(ap[-1], a))
        bp.append(zy_mul(bp[-1], b))
    f_el = [0] * 6
    for d, c in F_TERMS:
        term = zy_mul(c, zy_mul(ap[d], bp[12 - d]))
        f_el = [x + y for x, y in zip(f_el, term)]
    g_el = [x + y for x, y in zip(ap[2], zy_mul([0, U_ABS], bp[2]))]
    return zy_norm(f_el), zy_norm(g_el)


def inverse_q(matrix):
    n = len(matrix)
    a = [[F(x) for x in row] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(matrix)]
    for k in range(n):
        j = next(j for j in range(k, n) if a[j][k])
        a[k], a[j] = a[j], a[k]
        p = a[k][k]
        a[k] = [x / p for x in a[k]]
        for j in range(n):
            if j != k and a[j][k]:
                p = a[j][k]
                a[j] = [x - p * y for x, y in zip(a[j], a[k])]
    return [row[n:] for row in a]


def integer_metric(matrix):
    den = math.lcm(*(F(x).denominator for row in matrix for x in row))
    q = [[int(F(x) * den) for x in row] for row in matrix]
    return q, den, inverse_q(q)


def quadratic_int(q, a):
    s = 0
    for i in range(6):
        ai = a[i]
        qi = q[i]
        s += ai * (qi[0] * a[0] + qi[1] * a[1] + qi[2] * a[2] + qi[3] * a[3] + qi[4] * a[4] + qi[5] * a[5])
    return s


def integer_interval(a, t, c):
    a, t, c = F(a), F(t), F(c)
    if a <= 0:
        raise ValueError("non-positive pivot")
    centre = -t / a
    radius2 = (t * t - a * c) / (a * a)
    if radius2 < 0:
        return range(0)
    n, d = centre.numerator, centre.denominator
    z = math.isqrt((radius2.numerator * d * d) // radius2.denominator)
    lo = -((-(n - z)) // d)
    hi = (n + z) // d
    return range(lo, hi + 1)


KMAX_RHO = 64
RSTEP = 1.0 / 128.0


@functools.lru_cache(None)
def _build_rho_coefs(kmax=KMAX_RHO, nt=600):
    with localcontext() as ctx:
        ctx.prec = 200
        zero, one = D(0), D(1)
        coefs = [None, [one] + [zero] * (nt - 1)]
        for k in range(2, kmax + 1):
            prev, c = coefs[k - 1], [zero] * nt
            for m in range(nt - 1):
                c[m + 1] = (prev[m] + m * c[m]) / (k * (m + 1))
            c[0] = prev[0] - sum(c[1:])
            coefs.append(c)
        out = [None]
        for c in coefs[1:]:
            big = max(abs(x) for x in c)
            while len(c) > 1 and abs(c[-1]) < big * D("1e-25"):
                c = c[:-1]
            out.append([float(x) for x in c])
        return out


RHO_COEFS = _build_rho_coefs()


def rho_series(u):
    if u <= 1.0:
        return 1.0
    k = math.ceil(u)
    if k > KMAX_RHO:
        return 0.0
    x, val = k - u, 0.0
    for c in reversed(RHO_COEFS[k]):
        val = val * x + c
    return max(val, 0.0)


RHO_TAB = [rho_series(i * RSTEP) for i in range(int(KMAX_RHO / RSTEP) + 2)]


def rho_lo(u):
    if u <= 1.0:
        return 1.0
    i = math.ceil(u / RSTEP)
    return RHO_TAB[i] if i < len(RHO_TAB) else 0.0


def rho_hi(u):
    if u <= 1.0:
        return 1.0
    i = math.floor(u / RSTEP)
    if i + 1 >= len(RHO_TAB):
        return 1.0 / math.factorial(math.floor(u))
    return RHO_TAB[i]


def kappa_rho(u, t):
    r = rho_series(u)
    return rho_series(u - t) / r if r > 0.0 else 0.0


def p_trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def p_mod(a, m, p):
    a = p_trim([x % p for x in a])
    m = p_trim([x % p for x in m])
    dm = len(m) - 1
    inv = pow(m[-1], -1, p)
    while len(a) - 1 >= dm:
        c = a[-1] * inv % p
        sh = len(a) - 1 - dm
        for i, mi in enumerate(m):
            a[sh + i] = (a[sh + i] - c * mi) % p
        p_trim(a)
    return a


def p_mul(a, b, p):
    if not a or not b:
        return []
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] = (c[i + j] + x * y) % p
    return p_trim(c)


def p_add(a, b, p):
    return p_trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(max(len(a), len(b)))])


def p_sub(a, b, p):
    return p_trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(max(len(a), len(b)))])


def p_scal(a, c, p):
    return p_trim([(x * c) % p for x in a])


def p_div(a, b, p):
    a = p_trim([x % p for x in a])
    b = p_trim([x % p for x in b])
    q = [0] * max(1, len(a) - len(b) + 1)
    inv = pow(b[-1], -1, p)
    while a and len(a) >= len(b):
        c = a[-1] * inv % p
        sh = len(a) - len(b)
        q[sh] = c
        for i, bi in enumerate(b):
            a[sh + i] = (a[sh + i] - c * bi) % p
        p_trim(a)
    return p_trim(q), a


def p_gcd(a, b, p):
    a = p_trim([x % p for x in a])
    b = p_trim([x % p for x in b])
    while b:
        a, b = b, p_mod(a, b, p)
    if a:
        inv = pow(a[-1], -1, p)
        a = [(x * inv) % p for x in a]
    return a


def p_deriv(a, p):
    return p_trim([(i * a[i]) % p for i in range(1, len(a))])


def p_powmod(a, e, m, p):
    r = [1]
    a = p_mod(a, m, p)
    while e:
        if e & 1:
            r = p_mod(p_mul(r, a, p), m, p)
        a = p_mod(p_mul(a, a, p), m, p)
        e >>= 1
    return r


def verify_tower_isomorphism_and_fields():
    f0 = [0] * 6
    for d, c in F_TERMS:
        term = (c + [0] * 6)[:6]
        for _ in range(d // 2):
            term = zy_mul(term, [0, -U_ABS])
        f0 = [x + y for x, y in zip(f0, term)]
    res_val = zy_norm(f0)
    assert res_val != 0 and res_val % P_MOD == 0

    hp = [x % P_MOD for x in H_POLY]
    yp = p_powmod([0, 1], P_MOD, hp, P_MOD)
    cur = [0, 1]
    for k in range(1, 7):
        cur = p_powmod(cur, P_MOD, hp, P_MOD) if k > 1 else yp
        diff = p_sub(cur, [0, 1], P_MOD)
        if k in (1, 2, 3):
            assert p_gcd(hp, diff, P_MOD) == [1]
        elif k == 6:
            assert diff == []
    uy_chi = p_powmod([0, (-U_ABS) % P_MOD], (P_MOD**6 - 1) // 2, hp, P_MOD)
    assert uy_chi == [P_MOD - 1]

    deriv = [i * H_POLY[i] for i in range(1, 7)]
    sylv = [[0] * i + list(reversed(H_POLY)) + [0] * (4 - i) for i in range(5)] + [
        [0] * i + list(reversed(deriv)) + [0] * (5 - i) for i in range(6)
    ]
    disc_h = (-1) ** 15 * bareiss_det(sylv)
    assert disc_h == 58673 and 23 * 2551 == 58673

    hr = [x % R_ORD for x in H_POLY]

    def res_mod_r(a, b):
        a = p_trim([x % R_ORD for x in a])
        b = p_trim([x % R_ORD for x in b])
        s = 1
        while b and len(b) > 1:
            da, db = len(a) - 1, len(b) - 1
            rem = p_mod(a, b, R_ORD)
            if not rem:
                return 0
            s = s * ((-1) ** (da * db)) * pow(b[-1], da - (len(rem) - 1), R_ORD) % R_ORD
            a, b = b, rem
        return s * pow(b[0], len(a) - 1, R_ORD) % R_ORD if b else 0

    def norm_poly_mod_r(terms, deg):
        n = 6 * deg + 1
        xs = list(range(n))
        ys = [res_mod_r(hr, [sum(c[i] * pow(k, d, R_ORD) for d, c in terms if i < len(c)) % R_ORD for i in range(6)]) for k in xs]
        for j in range(1, n):
            for i in range(n - 1, j - 1, -1):
                ys[i] = (ys[i] - ys[i - 1]) * pow(xs[i] - xs[i - j], -1, R_ORD) % R_ORD
        poly = [ys[-1]]
        for i in range(n - 2, -1, -1):
            poly = [(-xs[i] * poly[0] + ys[i]) % R_ORD] + [(poly[k - 1] - xs[i] * poly[k]) % R_ORD for k in range(1, len(poly))] + [poly[-1]]
        return poly

    pf_r = norm_poly_mod_r(F_TERMS, 12)
    pg_r = norm_poly_mod_r(G_TERMS, 2)
    sqf_f = len(pf_r) == 73 and len(p_gcd(pf_r, p_deriv(pf_r, R_ORD), R_ORD)) == 1
    sqf_g = len(pg_r) == 13 and len(p_gcd(pg_r, p_deriv(pg_r, R_ORD), R_ORD)) == 1
    assert sqf_f and sqf_g

    m_odd = U_ABS >> 16
    d_kg = (4**6) * (m_odd**6) * (disc_h**2)
    log2_d = math.log2(d_kg)
    log2_kappa = 11.0 * math.log2(math.e * math.log(d_kg) / 22.0)
    w_tors = max(w for w in range(1, 1000) if sum(1 for k in range(1, w + 1) if math.gcd(k, w) == 1) <= 12)
    log2_h_kg = math.log2(w_tors) + 0.5 * log2_d + log2_kappa - math.log2(0.2052)
    assert log2_h_kg < math.log2(R_ORD)

    return {
        "resultant_divisible_by_p": True,
        "h_irreducible_mod_p": True,
        "uy_non_square_in_Fp6": True,
        "disc_h": disc_h,
        "disc_h_factors": [23, 2551],
        "maximal_order_power_basis": True,
        "unramified_mod_r_f_deg72": sqf_f,
        "unramified_mod_r_g_deg12": sqf_g,
        "log2_h_Kg_upper": round(log2_h_kg, 2),
        "log2_r": round(math.log2(R_ORD), 2),
    }


def primes_range(lo, hi):
    return [p for p in range(lo, hi) if all(p % q for q in range(2, math.isqrt(p) + 1))]


def roots_fp(a, p):
    return [r for r in range(p) if sum(c * pow(r, i, p) for i, c in enumerate(a)) % p == 0]


def factor_fp(a, p, seed=1):
    rng = random.Random(seed)

    def split(poly, d):
        if len(poly) - 1 == d:
            return [poly]
        for _ in range(10000):
            r = p_trim([rng.randrange(p) for _ in range(len(poly) - 1)])
            if len(r) < 2:
                continue
            if p == 2:
                t, s = r, r
                for _ in range(d - 1):
                    t = p_mod(p_mul(t, t, p), poly, p)
                    s = p_add(s, t, p)
                g = p_gcd(poly, s, p)
            else:
                g = p_gcd(poly, p_add(p_powmod(r, (p**d - 1) // 2, poly, p), [p - 1], p), p)
            if 1 < len(g) - 1 < len(poly) - 1:
                q, _ = p_div(poly, g, p)
                return split(g, d) + split(q, d)
        raise RuntimeError("factor_fp split failed")

    out = []
    rest = p_trim([x % p for x in a])
    for r in roots_fp(rest, p):
        out.append([(-r) % p, 1])
        rest, _ = p_div(rest, [(-r) % p, 1], p)
    d = 2
    x = [0, 1]
    xp = p_powmod(x, p, rest, p) if len(rest) > 1 else x
    while len(rest) - 1 >= 2 * d:
        xp = p_powmod(xp, p, rest, p)
        g = p_gcd(rest, p_add(xp, p_scal(x, p - 1, p), p), p)
        if len(g) > 1:
            out += split(g, d)
            rest, _ = p_div(rest, g, p)
            if len(rest) > 1:
                xp = p_mod(xp, rest, p)
        d += 1
    if len(rest) > 1:
        out.append(rest)
    return out


def e_mul(a, b, m, p):
    return p_mod(p_mul(a, b, p), m, p) if a and b else []


def e_inv(a, m, p):
    r0, r1 = list(m), p_trim(list(a))
    s0, s1 = [], [1]
    while r1:
        q, r = p_div(r0, r1, p)
        r0, r1 = r1, r
        s0, s1 = s1, p_add(s0, p_scal(p_mul(q, s1, p), p - 1, p), p)
    return p_mod(p_scal(s0, pow(r0[0], -1, p), p), m, p)


def E_trim(A):
    while A and not A[-1]:
        A.pop()
    return A


def E_mul(A, B, m, p):
    C = [[] for _ in range(len(A) + len(B) - 1)] if A and B else []
    for i, x in enumerate(A):
        if x:
            for j, y in enumerate(B):
                if y:
                    C[i + j] = p_add(C[i + j], e_mul(x, y, m, p), p)
    return E_trim(C)


def E_mod(A, M, m, p):
    A = E_trim([list(x) for x in A])
    dm = len(M) - 1
    inv = e_inv(M[-1], m, p)
    while len(A) - 1 >= dm:
        c = e_mul(A[-1], inv, m, p)
        sh = len(A) - 1 - dm
        for i in range(len(M)):
            A[sh + i] = p_add(A[sh + i], p_scal(e_mul(c, M[i], m, p), p - 1, p), p)
        E_trim(A)
    return A


def E_gcd(A, B, m, p):
    A = E_trim([list(x) for x in A])
    B = E_trim([list(x) for x in B])
    while B:
        A, B = B, E_mod(A, B, m, p)
    if A:
        inv = e_inv(A[-1], m, p)
        A = [e_mul(x, inv, m, p) for x in A]
    return A


def E_powmod(A, e, M, m, p):
    R = [[1]]
    A = E_mod(A, M, m, p)
    while e:
        if e & 1:
            R = E_mod(E_mul(R, A, m, p), M, m, p)
        A = E_mod(E_mul(A, A, m, p), M, m, p)
        e >>= 1
    return R


def E_quo(A, M, m, p):
    A = E_trim([list(x) for x in A])
    dm = len(M) - 1
    Q = [[] for _ in range(max(0, len(A) - dm))]
    while len(A) - 1 >= dm:
        c = A[-1]
        sh = len(A) - 1 - dm
        Q[sh] = c
        for i in range(len(M)):
            A[sh + i] = p_add(A[sh + i], p_scal(e_mul(c, M[i], m, p), p - 1, p), p)
        E_trim(A)
    return E_trim(Q)


def poly_over_ext(terms, m, p):
    deg = max(d for d, _ in terms)
    co = [[] for _ in range(deg + 1)]
    for d, c in terms:
        v = []
        for i, ci in enumerate(c):
            if ci:
                v = p_add(v, p_scal(p_powmod([0, 1], i, m, p), ci % p, p), p)
        co[d] = p_mod(v, m, p) if v else []
    return co


def root_multiplicities_ext(A, q, m, p):
    A = E_trim([list(x) for x in A])
    XP = E_powmod([[], [1]], q, A, m, p)
    diff = E_trim([p_add(XP[i] if i < len(XP) else [], p_scal([1] if i == 1 else [], p - 1, p), p) for i in range(max(len(XP), 2))])
    R = E_gcd(A, diff, m, p)
    if p not in (2, 19, 23):
        return [1] * (len(R) - 1)
    counts = []
    while len(R) > 1:
        counts.append(len(R) - 1)
        A = E_quo(A, R, m, p)
        R = E_gcd(A, R, m, p)
    counts.append(0)
    return [k + 1 for k in range(len(counts) - 1) for _ in range(counts[k] - counts[k + 1])]


@functools.lru_cache(None)
def structure_of_prime(p):
    Hp = [x % p for x in H_POLY]
    sq = p_gcd(Hp, p_deriv(Hp, p), p)
    out = []
    for hi in factor_fp(p_div(Hp, sq, p)[0] if len(sq) > 1 else Hp, p):
        e, rest = 0, Hp
        while True:
            quo, rem = p_div(rest, hi, p)
            if rem:
                break
            rest, e = quo, e + 1
        d1 = len(hi) - 1
        q = p**d1
        mf = root_multiplicities_ext(poly_over_ext(F_TERMS, hi, p), q, hi, p)
        mg = root_multiplicities_ext(poly_over_ext(G_TERMS, hi, p), q, hi, p)
        out.append([d1, e, mf, mg])
    return out


def compute_prime_structure_and_alphas(pmax=512):
    inv_zeta_2 = 1.0
    alpha_f_e = 0.0
    alpha_g_e = 0.0
    special_primes = {}
    for p in primes_range(2, pmax):
        st = structure_of_prime(p)
        if any(e != 1 or any(m != 1 for m in mf + mg) for d1, e, mf, mg in st):
            special_primes[p] = st
        vf = vg = 0.0
        for d1, e, mf, mg in st:
            q = float(p**d1)
            inv_zeta_2 *= 1.0 - 1.0 / (q * q)
            for side_idx, mults in enumerate((mf, mg)):
                v_id = sum(q / (q * q - 1.0) if mu == 1 else 1.0 / (q + 1.0) for mu in mults)
                if side_idx == 0:
                    vf += d1 * v_id
                else:
                    vg += d1 * v_id
        alpha_f_e += math.log(p) * (1.0 / (p - 1.0) - vf)
        alpha_g_e += math.log(p) * (1.0 / (p - 1.0) - vg)

    return {
        "pmax_exact": pmax,
        "inv_zeta_Kh_2": round(inv_zeta_2, 6),
        "alpha_f_bits_computed": round(alpha_f_e / LN2, 3),
        "alpha_g_bits_computed": round(alpha_g_e / LN2, 3),
        "alpha_f_bits_charged": 5.0,
        "alpha_g_bits_charged": 1.0,
        "special_primes_below_512": special_primes,
    }


def c_add(a, b):
    return a[0] + b[0], a[1] + b[1]


def c_mul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def c_div(a, b):
    d = b[0] * b[0] + b[1] * b[1]
    return (a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d


def c_eval(coeffs, z):
    out = (F(0), F(0))
    for c in reversed(coeffs):
        out = c_add(c_mul(out, z), (F(c), F(0)))
    return out


@functools.lru_cache(None)
def root_disks():
    roots = [1.3 * complex(math.cos(2 * math.pi * (i + 0.2) / 6), math.sin(2 * math.pi * (i + 0.2) / 6)) for i in range(6)]
    for _ in range(150):
        new = []
        for i, z in enumerate(roots):
            num = sum(c * (z**j) for j, c in enumerate(H_POLY))
            den = math.prod(z - w for j, w in enumerate(roots) if i != j)
            new.append(z - num / den)
        if max(abs(x - y) for x, y in zip(new, roots)) < 1e-15:
            roots = new
            break
        roots = new
    quantum, radius = 1 << 128, F(1, 1 << 100)
    deriv = [j * H_POLY[j] for j in range(1, 7)]
    disks = []
    for z in roots:
        centre = (F(z.real), F(0) if abs(z.imag) < 1e-12 else F(z.imag))
        for _ in range(8):
            corr = c_div(c_eval(H_POLY, centre), c_eval(deriv, centre))
            upd = (centre[0] - corr[0], centre[1] - corr[1])
            centre = tuple(F(round(x * quantum), quantum) for x in upd)
        coeffs = [c_eval([math.comb(j, k) * H_POLY[j] for j in range(k, 7)], centre) for k in range(7)]
        lower_lin = max(abs(x) for x in coeffs[1]) * radius
        upper_rem = sum(abs(x) for x in coeffs[0]) + sum(sum(abs(x) for x in coeffs[k]) * (radius**k) for k in range(2, 7))
        assert lower_lin > upper_rem
        disks.append((centre, radius))
    return sorted(disks, key=lambda x: (x[0][1] != 0, x[0][0], x[0][1]))


def atan_interval(q, terms=80):
    val = sum((F((-1) ** k, (2 * k + 1) * (q ** (2 * k + 1))) for k in range(terms)), F(0))
    rem = F((-1) ** terms, (2 * terms + 1) * (q ** (2 * terms + 1)))
    return min(val, val + rem), max(val, val + rem)


def certify_ellipsoid_counts():
    a5, b239 = atan_interval(5), atan_interval(239)
    pi_lo, pi_hi = 16 * a5[0] - 4 * b239[1], 16 * a5[1] - 4 * b239[0]
    v6_lo = (pi_lo**3) / 6
    v5_hi = F(8, 15) * (pi_hi**2)
    out = {}
    for name, mat in (("a", MATRIX_A), ("b", MATRIX_B)):
        q, den, _ = integer_metric(mat)
        d6 = F(bareiss_det(q), den**6)
        d5 = F(bareiss_det([r[:5] for r in q[:5]]), den**5)
        scale = 1 << 64
        sd6_hi = F(math.isqrt((d6.numerator * scale * scale) // d6.denominator) + 1, scale)
        sd5_lo = F(math.isqrt((d5.numerator * scale * scale) // d5.denominator), scale)
        sec_hi = (v5_hi / v6_lo) * (sd6_hi / sd5_lo)
        r, t = F(1024), F(2)
        lo = (1 - t / r) ** 6 - sec_hi * ((1 + t / r) ** 5) / r
        hi = (1 + t / r) ** 6
        assert lo >= F(98, 100) and hi <= F(102, 100)
        out[name] = (float(lo), float(hi))
    return out


class Interval:
    def __init__(self, lo, hi=None):
        self.lo = float(lo)
        self.hi = float(lo if hi is None else hi)

    @classmethod
    def rational(cls, lo, hi=None):
        lo, hi = F(lo), F(lo if hi is None else hi)
        return cls(math.nextafter(float(lo), -math.inf), math.nextafter(float(hi), math.inf))

    def __add__(self, other):
        if not isinstance(other, Interval):
            other = Interval.rational(other)
        return Interval(math.nextafter(self.lo + other.lo, -math.inf), math.nextafter(self.hi + other.hi, math.inf))

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not isinstance(other, Interval):
            other = Interval.rational(other)
        cands = (self.lo * other.lo, self.lo * other.hi, self.hi * other.lo, self.hi * other.hi)
        return Interval(math.nextafter(min(cands), -math.inf), math.nextafter(max(cands), math.inf))

    __rmul__ = __mul__

    def absolute(self):
        lo = 0.0 if self.lo <= 0.0 <= self.hi else min(abs(self.lo), abs(self.hi))
        return Interval(lo, max(abs(self.lo), abs(self.hi)))

    def square(self):
        a = self.absolute()
        return Interval(max(0.0, math.nextafter(a.lo * a.lo, -math.inf)), math.nextafter(a.hi * a.hi, math.inf))


class ComplexInterval:
    def __init__(self, real, imag):
        self.real, self.imag = real, imag

    @classmethod
    def constant(cls, val):
        return cls(Interval.rational(val), Interval(0.0))

    def __add__(self, other):
        if not isinstance(other, ComplexInterval):
            other = ComplexInterval.constant(other)
        return ComplexInterval(self.real + other.real, self.imag + other.imag)

    __radd__ = __add__

    def __mul__(self, other):
        if not isinstance(other, ComplexInterval):
            other = ComplexInterval.constant(other)
        return ComplexInterval(self.real * other.real - self.imag * other.imag, self.real * other.imag + self.imag * other.real)

    __rmul__ = __mul__

    def absolute_square(self):
        v = self.real.square() + self.imag.square()
        return Interval(max(0.0, v.lo), v.hi)


def _build_log_table():
    lower, upper = [], []
    for j in range(17):
        val = 1 << (j + 16 * 52)
        k = val
        for _ in range(4):
            k = math.isqrt(k)
        lower.append(k / (1 << 52))
        upper.append((k + int(k**16 != val)) / (1 << 52))
    return lower, upper


LOG_LOWER, LOG_UPPER = _build_log_table()


def scaled_log_float_bounds(lo, hi):
    if lo <= 0.0 or not math.isfinite(lo) or not math.isfinite(hi) or hi < lo:
        return None
    m, e = math.frexp(lo)
    lower = 16 * (e - 1) + bisect.bisect_right(LOG_UPPER, 2.0 * m) - 1
    m, e = math.frexp(hi)
    upper = 16 * (e - 1) + bisect.bisect_left(LOG_LOWER, 2.0 * m)
    return lower, upper


class NormBounds:
    def __init__(self, a, disks=None):
        self.a = tuple(a)
        self.evals = 0
        self.places = []
        for (r, i), radius in (root_disks() if disks is None else disks):
            if i < -radius:
                continue
            real = i == 0
            root = Interval.rational(r - radius, r + radius) if real else ComplexInterval(Interval.rational(r - radius, r + radius), Interval.rational(i - radius, i + radius))
            one = Interval(1.0) if real else ComplexInterval.constant(1)
            av = one * 0
            for c in reversed(a):
                av = av * root + c
            ap, yp = [one], [one]
            for _ in range(12):
                ap.append(ap[-1] * av)
            for _ in range(6):
                yp.append(yp[-1] * root)
            coeffs = [ap[12], -2 * root * ap[10], one * 0, 2 * yp[3] * ap[6], one * 0, yp[5] * ap[2], yp[6]]
            self.places.append((real, root, one, ap[2], U_ABS * root, coeffs))

    def adaptive_row(self, suffix, lo, hi, width=128):
        contexts = []
        for real, root, one, a2, uy, coeffs in self.places:
            bv = one * 0
            for c in reversed(suffix):
                bv = bv * root + c
            contexts.append((real, bv * root, a2, uy, coeffs))
        pending, leaves = [(lo, hi)], []
        while pending:
            l_cur, h_cur = pending.pop()
            self.evals += 1
            logs = [[0, 0], [0, 0]]
            for real, base, a2, uy, coeffs in contexts:
                b0 = Interval(l_cur, h_cur)
                bv = base + b0 if real else base + ComplexInterval(b0, Interval(0.0))
                w = bv * bv
                fv = coeffs[6]
                for k in range(5, -1, -1):
                    fv = fv * w
                    if k not in (2, 4):
                        fv = fv + coeffs[k]
                gv = a2 + uy * w
                for side, val in enumerate((fv, gv)):
                    mag = val.absolute() if real else val.absolute_square()
                    bnd = scaled_log_float_bounds(mag.lo, mag.hi)
                    if bnd is None:
                        logs[side] = None
                    elif logs[side] is not None:
                        logs[side][0] += bnd[0]
                        logs[side][1] += bnd[1]
            if all(v is not None and v[1] - v[0] <= width for v in logs):
                leaves.append((l_cur, h_cur, logs))
            elif l_cur == h_cur:
                leaves.append((l_cur, h_cur, None))
            else:
                mid = (l_cur + h_cur) // 2
                pending.extend(((mid + 1, h_cur), (l_cur, mid)))
        return leaves


@functools.lru_cache(None)
def get_small_units(bound=2):
    out = []
    for coeffs in itertools.product(range(-bound, bound + 1), repeat=6):
        if any(coeffs) and coeffs not in ((1, 0, 0, 0, 0, 0), (-1, 0, 0, 0, 0, 0)):
            if zy_norm(list(coeffs)) == 1:
                out.append(list(coeffs))
    return out


def row_hnf_6(rows, d):
    pending = [[x % d for x in r] for r in rows]
    hnf = []
    for col in range(6):
        pending.append([d if j == col else 0 for j in range(6)])
        k = next(j for j, r in enumerate(pending) if r[col])
        piv = pending.pop(k)
        for j, r in enumerate(pending):
            if r[col]:
                a, b = piv[col], r[col]
                x, y, u, v = 1, 0, 0, 1
                aa, bb = a, b
                while bb:
                    q, rem = divmod(aa, bb)
                    aa, bb = bb, rem
                    x, u, y, v = u, x - q * u, v, y - q * v
                if aa < 0:
                    aa, x, y = -aa, -x, -y
                pnew = [x * px + y * rx for px, rx in zip(piv, r)]
                rnew = [(b // aa) * px - (a // aa) * rx for px, rx in zip(piv, r)]
                for m in range(col + 1, 6):
                    pnew[m] %= d
                    rnew[m] %= d
                piv, pending[j] = pnew, rnew
        if piv[col] < 0:
            piv = [-x for x in piv]
        hnf.append(piv)
        pending = [r for r in pending if any(r)]
    return hnf


def coprime_pair(a, b):
    da, db = zy_norm(a), zy_norm(b)
    if da == 0 or db == 0:
        return False
    d = math.gcd(da, db)
    if d == 1:
        return True
    rows = [zy_mul(v, [0] * i + [1]) for v in (a, b) for i in range(6)]
    key = row_hnf_6(rows, d)
    return all(key[i][i] == 1 for i in range(6))


def canonical_a_fast_and_exact(a, qa):
    a_tup = tuple(a)
    qa_val = quadratic_int(qa, a)
    for eps in get_small_units(2):
        v = zy_mul(eps, a)
        if v[5] < 0:
            v = [-x for x in v]
        if v[5] > 0 and tuple(v) != a_tup:
            qv = quadratic_int(qa, v)
            if (qv, tuple(v)) < (qa_val, a_tup):
                return False

    na = zy_norm(a)
    if na == 0:
        return False
    b = [zy_mul(a, [0] * i + [1]) for i in range(6)]

    def inner_qa(u, v):
        s = 0
        for i in range(6):
            ui, qi = u[i], qa[i]
            s += ui * (qi[0] * v[0] + qi[1] * v[1] + qi[2] * v[2] + qi[3] * v[3] + qi[4] * v[4] + qi[5] * v[5])
        return s

    def gs():
        gram = [[inner_qa(u, v) for v in b] for u in b]
        mu, lengths = [], []
        for i in range(6):
            row = []
            for j in range(i):
                row.append((F(gram[i][j]) - sum(row[k] * mu[j][k] * lengths[k] for k in range(j))) / lengths[j])
            length = F(gram[i][i]) - sum((row[k] ** 2) * lengths[k] for k in range(i))
            mu.append(row)
            lengths.append(length)
        return mu, lengths

    k = 1
    while k < 6:
        mu, lengths = gs()
        for j in range(k - 1, -1, -1):
            z = (mu[k][j] + F(1, 2)).numerator // (mu[k][j] + F(1, 2)).denominator
            if z:
                b[k] = [x - z * y for x, y in zip(b[k], b[j])]
                mu, lengths = gs()
        if lengths[k] >= (F(3, 4) - mu[k][k - 1] ** 2) * lengths[k - 1]:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            k = max(1, k - 1)

    q_trans = [[F(inner_qa(u, v)) for v in b] for u in b]
    projs, q_cur = [q_trans], q_trans
    while len(q_cur) > 1:
        q_cur = [[q_cur[i][j] - q_cur[i][0] * q_cur[0][j] / q_cur[0][0] for j in range(1, len(q_cur))] for i in range(1, len(q_cur))]
        projs.append(q_cur)

    x_vec = [0] * 6

    def walk(idx):
        q_p = projs[idx]
        suf = x_vec[idx + 1 :]
        t_val = sum(q_p[0][j + 1] * v for j, v in enumerate(suf))
        c_val = sum(q_p[j + 1][m + 1] * v * w for j, v in enumerate(suf) for m, w in enumerate(suf)) - F(qa_val)
        for val in integer_interval(q_p[0][0], t_val, c_val):
            x_vec[idx] = val
            if idx:
                if not walk(idx - 1):
                    return False
            else:
                v = tuple(sum(x_vec[i] * b[i][j] for i in range(6)) for j in range(6))
                if v[5] > 0 and v != a_tup and zy_norm(list(v)) == na:
                    if (quadratic_int(qa, v), v) < (qa_val, a_tup):
                        return False
        return True

    return walk(5)


def _series_coeffs(weights, normalizer, length):
    terms = [(j, j * v) for j, v in enumerate(weights) if j and abs(v) >= 1e-24]
    out = [math.exp(-normalizer)] + [0.0] * length
    for n in range(1, length + 1):
        val = 0.0
        for j, w in terms:
            if j > n:
                break
            val += w * out[n - j]
        out[n] = val / n
    return out


class RetentionModel:
    def __init__(self, B_bits, X_bits, step=0.25, max_bits=850, tilt=0.05):
        self.B, self.X, self.minimum, self.step = B_bits, X_bits, 256, step
        self.tilt, self.length = tilt, math.ceil(max_bits / step)
        names = ("small", "sieve_geo", "large_geo", "sieve_sf", "large_sf")
        weights = {k: [0.0] * (self.length + 1) for k in names}
        totals = Counter()
        cutoff = min(12.0, float(B_bits))
        for p in primes_range(2, int(2**cutoff)):
            b = math.log2(p)
            b_idx = round(b / step)
            q = math.exp(tilt * b_idx * step - LN2 * b)
            groups = [("small", False)] if p < 256 else ([("sieve_geo", False), ("sieve_sf", True)] if b < X_bits else [("large_geo", False), ("large_sf", True)])
            for name, sqf in groups:
                totals[name] += math.log1p(q) if sqf else -math.log1p(-q)
                qk = q
                for k in range(1, self.length + 1):
                    idx = k * b_idx
                    if idx > self.length:
                        break
                    sgn = -1 if sqf and k % 2 == 0 else 1
                    weights[name][idx] += sgn * qk / k
                    qk *= q
                    if qk < 1e-25:
                        break

        def continuum(left, right, pref):
            if left >= right:
                return
            w_step = step / 4.0
            cells = math.ceil((right - left) / w_step)
            w_step = (right - left) / cells
            for i in range(cells):
                b = left + (i + 0.5) * w_step
                b_idx = round(b / step)
                for k in range(1, 9):
                    amt = math.exp(k * tilt * b_idx * step - (k - 1) * LN2 * b) * w_step / (k * b)
                    idx = k * b_idx
                    sgn = -1 if k % 2 == 0 else 1
                    totals[pref + "_geo"] += amt
                    totals[pref + "_sf"] += sgn * amt
                    if idx <= self.length:
                        weights[pref + "_geo"][idx] += amt
                        weights[pref + "_sf"][idx] += sgn * amt

        continuum(cutoff, min(X_bits, B_bits), "sieve")
        continuum(max(cutoff, X_bits), B_bits, "large")
        self.pdf = {name: _series_coeffs(weights[name], totals[name], self.length) for name in ("sieve_sf", "sieve_geo")}
        self.missing_sf = _series_coeffs([x + y for x, y in zip(weights["small"], weights["large_sf"])], totals["small"] + totals["large_sf"], self.length)
        self.missing_geo = _series_coeffs([x + y for x, y in zip(weights["small"], weights["large_geo"])], totals["small"] + totals["large_geo"], self.length)
        all_w = [sum(weights[k][j] for k in ("small", "sieve_geo", "large_geo")) for j in range(self.length + 1)]
        self.denom = _series_coeffs(all_w, sum(totals[k] for k in ("small", "sieve_geo", "large_geo")), self.length)
        self.sel_norm = math.exp(totals["sieve_sf"] + totals["large_sf"] - totals["sieve_geo"] - totals["large_geo"])
        self._cdfs = {}

    def _mass(self, n, stop, sqf):
        if stop < 0:
            return 0.0
        key = (n, sqf)
        if key not in self._cdfs:
            miss = self.missing_sf if sqf else self.missing_geo
            sieve = self.pdf["sieve_sf" if sqf else "sieve_geo"]
            vals, tot = array("d"), 0.0
            for d in range(n + 1):
                tot += miss[d] * sieve[n - d]
                vals.append(tot)
            self._cdfs[key] = vals
        return self._cdfs[key][min(n, stop)]

    def retention_at(self, norm_bits, slack_bits, scale=16):
        n = round(norm_bits / self.step)
        cnt_cap = math.floor(math.ceil(norm_bits) / 8.0)
        allowed = slack_bits - (cnt_cap + 1) / scale
        stop = min(n, math.floor(allowed / self.step))
        return max(0.0, min(1.0, self.sel_norm * self._mass(n, stop, True) / self.denom[n]))

    def survivor_upper(self, norm_bits, slack_bits, extra_bits=8):
        n = round(norm_bits / self.step)
        stop = min(n, math.floor((slack_bits + extra_bits) / self.step))
        return self._mass(n, stop, False) / self.denom[n]


LINE_TAB = [(13.82, 24.10), (9.82, 31.11), (5.82, 88.31), (1.82, 310.70)]
REC_TAB = [(18.91, 44.90), (10.91, 78.21), (3.91, 266.97), (0.91, 938.08), (-1.09, 2231.66)]
ECM_TAB = {
    20: (11e3, 86), 25: (5e4, 214), 30: (25e4, 430), 35: (1e6, 904),
    40: (3e6, 2350), 45: (11e6, 4480), 50: (43e6, 7553), 55: (11e7, 17769),
    60: (26e7, 42017), 65: (85e7, 69408), 70: (29e8, 102212), 75: (76e8, 188056), 80: (25e9, 265557),
}
MONT5 = 367.76


def log2_li(log2x):
    x = log2x * LN2
    return log2x - math.log2(x) + math.log2(1.0 + 1.0 / x + 2.0 / (x * x))


def interp_tab(t, log2h):
    if log2h >= t[0][0]:
        return t[0][1]
    for (h1, c1), (h2, c2) in zip(t, t[1:]):
        if log2h >= h2:
            return c2 + (c1 - c2) * (log2h - h2) / (h1 - h2)
    return t[-1][1]


def main_ecm_instr(fbits, nbits, success=0.99):
    d = fbits / 3.3219
    ks = sorted(ECM_TAB)
    k = min(x for x in ks if x >= d) if d <= ks[-1] else ks[-1]
    b1, cur = ECM_TAB[k]
    curves = cur * math.log(1.0 / (1.0 - success))
    L = max(5, math.ceil(nbits / 64))
    return curves * 2.0 * 1.44 * b1 * 11.0 * MONT5 * ((L / 5.0) ** 2)


def lam_top(u, t0, m=40):
    r = rho_series(u)
    if r <= 0.0 or t0 >= 1.0:
        return 0.0
    return sum(rho_series(u - t) / r / t * ((1.0 - t0) / m) for t in (t0 + (1.0 - t0) * (i + 0.5) / m for i in range(m)))


def uniform_half_ellipsoid(q, inv_q, rad2, rng):
    cap = [math.isqrt((F(rad2) * inv_q[j][j]).numerator // (F(rad2) * inv_q[j][j]).denominator) for j in range(6)]
    while True:
        a = [rng.randrange(-c, c + 1) for c in cap[:-1]] + [rng.randrange(1, cap[-1] + 1)]
        q_val = quadratic_int(q, a)
        if q_val <= rad2:
            return a, q_val


def evaluate_main_phase_live(n_draws=2500, n_bound_checks=24):
    A2, theta = 2097152, F(1, 2)
    lb, lx2, excess = 69.0, 63.0, 1.3
    Sf, Sg, W = 292, 220, 8
    alpha_f, alpha_g = 5.0, 1.0
    sel_cells, sel_yield = 0.8078, 0.9617

    qa, da, ia = integer_metric(MATRIX_A)
    qb, db, ib = integer_metric(MATRIX_B)
    det_a, det_b = F(bareiss_det(qa), da**6), F(bareiss_det(qb), db**6)
    v6 = (math.pi**3) / 6.0
    va = v6 * ((float(theta) * A2) ** 3) / (2.0 * math.sqrt(float(det_a)))
    vb = v6 * (A2**3) / (2.0 * math.sqrt(float(det_b)))
    inv55 = float(ib[5][5]) * db
    lower_cnt, upper_cnt = float(F(2401, 2500)), float(F(2601, 2500))

    disks = root_disks()
    rows = []
    n_bound_verified = 0
    rad_a = int(A2 * theta * da)

    for idx in range(n_draws):
        rng = random.Random(f"v11-production-ideal-key:{1 + (idx // 1000)}:{idx % 1000}")
        a, aq = uniform_half_ellipsoid(qa, ia, rad_a, rng)
        if not canonical_a_fast_and_exact(a, qa):
            continue
        b, _ = uniform_half_ellipsoid(qb, ib, (A2 - F(aq, da)) * db, rng)
        qa_frac = float(F(aq, da)) / A2
        nf, ng = exact_norms(a, b)
        prim = coprime_pair(a, b)

        rad_b = (A2 - F(aq, da)) * db
        t_val = sum(qb[0][j] * b[j] for j in range(1, 6))
        c_val = sum(qb[j][k] * b[j] * b[k] for j in range(1, 6) for k in range(1, 6)) - rad_b
        iv = integer_interval(qb[0][0], t_val, c_val)
        blo = (b[0] // 64) * 64
        lo, hi = max(iv.start, blo), min(iv.stop - 1, blo + 63)
        rlen = 16 * max(1, (iv.stop - iv.start) // 16)
        blen = max(1, hi - lo + 1)

        if n_bound_verified < n_bound_checks:
            est = NormBounds(a, disks)
            leaves = est.adaptive_row(b[1:], lo, hi, width=128)
            at = [l for l in leaves if l[0] <= b[0] <= l[1]][0]
            assert at[2] is not None
            (fl, fu), (gl, gu) = at[2]
            assert fl <= 16.0 * math.log2(nf) <= fu and gl <= 16.0 * math.log2(ng) <= gu
            evals_cnt, leaves_cnt = est.evals, len(leaves)
            n_bound_verified += 1
        else:
            evals_cnt, leaves_cnt = max(1, round(0.04232 * blen)), max(1, round(0.02932 * blen))

        rows.append(dict(qa=qa_frac, nf=round(math.log2(nf) * 16.0) / 16.0, ng=round(math.log2(ng) * 16.0) / 16.0, primitive=prim, evals=evals_cnt, leaves=leaves_cnt, blen=blen, rlen=rlen))

    retain = RetentionModel(int(lb), int(lx2), step=0.25, max_bits=850, tilt=0.05)
    cost_ret = RetentionModel(850, int(lx2), step=0.25, max_bits=851, tilt=0.005)
    c_surv = (2.0**22) + (2.0**24.5) + main_ecm_instr(lb, Sf + 8) + main_ecm_instr(lb, Sg + 8)
    pairs_all = 2.0**log2_li(lx2)

    wsum = ysum = y2 = nf_m = ng_m = nf_r = ng_r = work_sum = 0.0
    for r in rows:
        w = (1.0 - r["qa"]) ** 3
        wsum += w
        if r["primitive"]:
            rf = retain.retention_at(r["nf"], Sf)
            rg = retain.retention_at(r["ng"], Sg)
            uf_r, ug_r = (r["nf"] + alpha_f) / lb, (r["ng"] + alpha_g) / lb
            miss = 0.99 ** (lam_top(uf_r, lx2 / lb) + lam_top(ug_r, lx2 / lb))
            y = rho_series(uf_r) * rho_series(ug_r) * rf * rg * miss
            ysum += w * y
            y2 += (w * y) ** 2
            nf_r += w * y * r["nf"]
            ng_r += w * y * r["ng"]
        nf_m += w * r["nf"]
        ng_m += w * r["ng"]

        S = vb * ((1.0 - r["qa"]) ** 3) * upper_cnt
        K = max(1.0, math.floor(math.sqrt((1.0 - r["qa"]) * A2 * inv55)))
        p_all = 12.325 + pairs_all * 500.0 / S + 5900.0 * r["evals"] / r["blen"] + 4200.0 / r["rlen"] + 15300.0 / S

        lp, line = 8.0, 0.0
        while lp < 20.0:
            lp2 = min(lp + 0.25, 20.0)
            line += 2.0 * math.log(lp2 / lp) * interp_tab(LINE_TAB, math.log2(S / K) - 0.5 * (lp + lp2))
            lp = lp2
        lp, rec, empty = 20.0, 0.0, 0.0
        while lp < lx2:
            lp2 = min(lp + 0.25, lx2)
            log2h = math.log2(S) - 0.5 * (lp + lp2)
            if log2h >= REC_TAB[-1][0]:
                rec += 2.0 * math.log(lp2 / lp) * interp_tab(REC_TAB, log2h)
            else:
                empty += (2.0**log2_li(lp2) - 2.0**log2_li(lp)) * 1600.0 / S
            lp = lp2

        hits = 2.0 * math.log(lx2 / 20.0)
        sf = cost_ret.survivor_upper(r["nf"], Sf, W)
        sg = cost_ret.survivor_upper(r["ng"], Sg, W)
        p_keep = line + rec + empty + hits * (7.0 + 7.003 + 36.0) + 8.0 + 8.0 + 21.0 * r["leaves"] / r["blen"] + 0.625 + sf * sg * c_surv
        work_sum += S * (p_all + sel_cells * p_keep)

    n_a_upper = va / n_draws
    population = n_a_upper * sum(vb * ((1.0 - r["qa"]) ** 3) * upper_cnt for r in rows)
    relcol_instr = n_a_upper * work_sum
    rels = n_a_upper * vb * lower_cnt * ysum * sel_yield
    yield_se = math.sqrt(max(0.0, y2 / n_draws - (ysum / n_draws) ** 2) / n_draws) / (ysum / n_draws) / LN2
    log2_U = math.log2(2.0**log2_li(lb) + 2.0**31)
    log2_need = log2_U + math.log2(excess)
    log2_rels = math.log2(rels)
    log2_relcol = math.log2(relcol_instr / UNIT)

    return {
        "n_draws": n_draws,
        "canonical_rows": len(rows),
        "primitive_rows": sum(1 for r in rows if r["primitive"]),
        "norm_bound_rows_verified": n_bound_verified,
        "nf_mean_bits": round(nf_m / wsum, 2),
        "ng_mean_bits": round(ng_m / wsum, 2),
        "nf_rel_weighted_bits": round(nf_r / ysum, 2),
        "ng_rel_weighted_bits": round(ng_r / ysum, 2),
        "log2_population_cells": round(math.log2(population), 4),
        "instructions_per_cell": round(relcol_instr / population, 2),
        "log2_U_orbits": round(log2_U, 4),
        "log2_need_1p3U": round(log2_need, 4),
        "log2_rels_obtained": round(log2_rels, 4),
        "margin_bits": round(log2_rels - log2_need, 4),
        "pool_se_log2_bits": round(yield_se, 4),
        "log2_relcol_units": round(log2_relcol, 4),
    }