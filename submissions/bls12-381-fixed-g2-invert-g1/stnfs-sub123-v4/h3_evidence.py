import math
import time
from array import array
from fractions import Fraction
from random import Random

T0 = time.time()

P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
U = 15132376222941642752
HH = [-1, -1, 1, 1, 0, 0, 1]
H_HI = [1, 0, 0, 1, 1, -1, -1]
H_DER = [6, 0, 0, 3, 2, -1]
FT = [(6, [1]), (5, [0, -2]), (3, [0, 0, 0, 2]), (1, [0, 0, 0, 0, 0, 1]), (0, [1, 1, -1, -1])]
TOY_U = 4

def say(*a):
    print(*a, flush=True)

def hred(c):
    c = list(c) + [0] * 14
    for i in range(len(c) - 1, 5, -1):
        v = c[i]
        if v:
            for j in range(6):
                c[i - 6 + j] -= v * HH[j]
        c[i] = 0
    return c[:6]

def pmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    c[i + j] += x * y
    return hred(c)

def pdet(M):
    M = [r[:] for r in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            r2 = next((j for j in range(k + 1, n) if M[j][k]), None)
            if r2 is None:
                return 0
            M[k], M[r2] = M[r2], M[k]; sign = -sign
        pk = M[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * pk - M[i][k] * M[k][j]) // prev
        prev = pk
    return sign * M[-1][-1]

def norm_h(w):
    m = [[pmul(w, [0] * i + [1])[j] for i in range(6)] for j in range(6)]
    return abs(pdet(m))

def snorm_h(w):
    m = [[pmul(w, [0] * i + [1])[j] for i in range(6)] for j in range(6)]
    return pdet(m)

def res_hi(pa, pb):
    da, db = len(pa) - 1, len(pb) - 1
    M = [[0] * (da + db) for _ in range(da + db)]
    for i in range(db):
        for j in range(da + 1):
            M[i][i + j] = pa[da - j]
    for i in range(da):
        for j in range(db + 1):
            M[db + i][i + j] = pb[db - j]
    return pdet(M)

def trace_elem(w, m):
    e = pmul(w, [0] * m + [1]) if m else w[:]
    mm = [[pmul(hred(e), [0] * i + [1])[j] for i in range(6)] for j in range(6)]
    return sum(mm[j][j] for j in range(6))

def htraces(mmax=13):
    return [trace_elem([1], m) for m in range(mmax + 1)]

def is_prime(n):
    if n < 2:
        return False
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % a == 0:
            return n == a
    d = n - 1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
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

def rho_factor(n):
    n = abs(n); fs = []
    d = 2
    while d * d <= n and d < 100000:
        while n % d == 0:
            fs.append(d); n //= d
        d += 1
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            fs.append(m); continue
        for seed in (2, 3, 5, 7, 11, 13, 17):
            x = y = seed; dd = 1; c = 1
            while dd == 1:
                x = (x * x + c) % m
                y = (y * y + c) % m
                y = (y * y + c) % m
                dd = math.gcd(abs(x - y), m)
            if dd != m:
                stack.append(dd); stack.append(m // dd)
                break
        else:
            raise ValueError("rho failed on %d" % m)
    return sorted(fs)

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def part_A():
    say("=" * 78)
    say("[A] EXACT DISCRIMINANT CONDITIONS MOD r — production tower")
    say("    K_h = Q[y]/(y^6+y^3+y^2-y-1);  F = z^6-2y z^5+2y^3 z^3+y^5 z+1+y-y^2-y^3 (z=x^2;")
    say("    f = F(x^2) degree 12 over K_h; the z-form Q[y,z]/(h,F) is the fixed field of x -> -x,")
    say("    degree 36 over Q - see [0d]-[0g]);  g = x^2+|u|y, |u| = 0x%x = %d" % (U, U))
    say("    r = 0x%x (%d bits);  p = %d bits" % (R, R.bit_length(), P.bit_length()))
    say("")
    dh_raw = res_hi(H_HI, H_DER)
    disc_h = (-1) ** 15 * dh_raw
    fs = rho_factor(disc_h)
    say("[A1] Res_y(h,h') = %d;  Disc(h) = (-1)^15 * Res = %d = %s" % (dh_raw, disc_h, fs))
    say("     squarefree: %s  =>  O_h = Z[y] and Disc(K_h) = Disc(h)" % (len(set(fs)) == len(fs)))
    say("     Disc(h) mod r = %d (nonzero: %s)" % (disc_h % R, disc_h % R != 0))
    assert disc_h == 58673 and fs == [23, 2551] and all(is_prime(x) for x in fs)
    say("")
    tr = htraces(13)
    A = [[2 * tr[a + b] for b in range(6)] for a in range(6)]
    Bg = [[-2 * U * tr[a + b + 1] for b in range(6)] for a in range(6)]
    G = [[0] * 12 for _ in range(12)]
    for a in range(6):
        for b in range(6):
            G[a][b] = A[a][b]; G[6 + a][6 + b] = Bg[a][b]
    dg = pdet(G)
    comp = disc_h ** 2 * (-(4 * U) ** 6)
    say("[A2] g side.  Order O_g' = Z[y][x]/(h, x^2+|u|y); trace Gram on the basis {y^a, y^a x}:")
    say("     Tr(y^a y^b) = 2 Tr_h(y^{a+b}), Tr(y^a x . y^b x) = -2|u| Tr_h(y^{a+b+1}), cross = 0.")
    say("     disc(O_g') = %d   (%d bits)" % (dg, abs(dg).bit_length()))
    say("     Disc(h)^2 * N_{K_h/Q}(-4|u|y) = %d ;  exact identity: %s" % (comp, dg == comp))
    uf = rho_factor(U >> 16)
    say("     N(y) = -h(0) = -1, so N(-4|u|y) = -(4|u|)^6 = -2^108 * %d^6 ;  %d = %s (prime factors)" %
        (U >> 16, U >> 16, uf))
    say("     disc(O_g') mod r = %d  (nonzero: %s)" % (dg % R, dg % R != 0))
    assert dg == comp and dg % R != 0
    say("     => r does not divide disc(O_g'); disc(K_g) = disc(O_g')/[O_Kg:O_g']^2 divides")
    say("        disc(O_g'), hence r does not divide Disc(K_g).")
    say("")
    def Fpoly_at(t):
        zc = [0] * 7
        for d, c in FT:
            zc[6 - d] += sum(ci * t ** i for i, ci in enumerate(c))
        return zc
    pts = [(t, res_hi(Fpoly_at(t), [(6 - i) * Fpoly_at(t)[i] for i in range(6)]))
           for t in range(-17, 17)]
    D = lagrange_int(pts)
    for t, v in pts[:3] + pts[-3:]:
        assert sum(c * t ** i for i, c in enumerate(D)) == v
    say("[A3] f side.  Order O_f' = Z[y][z]/(h, F) (degree 36).  D(y) := Res_z(F, dF/dz),")
    say("     interpolated from %d integer evaluations of the 11x11 Sylvester determinant:" % len(pts))
    say("     degree %d, max |coefficient| 2^%.1f, spot checks pass." %
        (len(D) - 1, math.log2(max(abs(c) for c in D))))
    Dr = hred(D)
    mD = [[pmul(Dr, [0] * i + [1])[j] for i in range(6)] for j in range(6)]
    ND = pdet(mD)
    rHD = res_hi(D[::-1], H_HI)
    say("     N_{K_h/Q}(D) = det(m_D in Z[y]/h) = %d (%d bits);  cross-check Res_y(h,D) = %d," %
        (ND, ND.bit_length(), rHD))
    say("     |equal|: %s" % (abs(rHD) == abs(ND)))
    assert abs(rHD) == abs(ND)
    Z6 = [0] * 6
    c = [Z6, [0, -2], Z6, [0, 0, 0, 2], Z6, [0, 0, 0, 0, 0, 1], [1, 1, -1, -1]]
    ps = [None] * 11
    ps[0] = [6]
    for k in range(1, 7):
        acc = [0] * 6
        for i in range(1, k):
            acc = [x + y for x, y in zip(acc, pmul(c[i], ps[k - i]))]
        acc = [x + k * v for x, v in zip(acc, c[k])]
        ps[k] = [-v for v in acc]
    for k in range(7, 11):
        acc = [0] * 6
        for i in range(1, 7):
            acc = [x + y for x, y in zip(acc, pmul(c[i], ps[k - i]))]
        ps[k] = [-v for v in acc]
    G36 = [[0] * 36 for _ in range(36)]
    for a in range(6):
        for b in range(6):
            for cc in range(6):
                for d in range(6):
                    G36[a + 6 * b][cc + 6 * d] = trace_elem(ps[b + d], a + cc)
    d36 = pdet(G36)
    comp = disc_h ** 6 * ND
    say("     disc(O_f') = det of the 36x36 trace Gram of {y^a z^b}, with power sums")
    say("     p_k = Tr_{K_f/K_h}(z^k) from Newton's identities = %d  (%d bits)" %
        (d36, abs(d36).bit_length()))
    say("     Disc(h)^6 * N_{K_h/Q}(disc_z F) = %d ;  exact identity: %s" % (comp, d36 == comp))
    assert d36 == comp
    say("     disc(O_f') mod r nonzero: %s;  N(D) mod r nonzero: %s;  Disc(h)^6 mod r nonzero: %s" %
        (d36 % R != 0, ND % R != 0, pow(disc_h, 6, R) != 0))
    assert d36 % R != 0 and ND % R != 0
    say("     => r does not divide disc(O_f'); disc(K_f) = disc(O_f')/[O_Kf:O_f']^2 divides")
    say("        disc(O_f'), hence r does not divide Disc(K_f).")
    say("")
    nmy = snorm_h([0, -U])
    say("[A4] K_g is a field: N_{K_h/Q}(-|u|y) = %d < 0, and squares have non-negative norm," % nmy)
    say("     so x^2+|u|y is irreducible over K_h.  (K_f being a field is the v21 §2 tower")
    say("     certificate [0][0b]; that irreducibility is assumed here, not recomputed.)")
    assert nmy < 0
    say("")
    return disc_h

def lagrange_int(pts):
    def polymul(a, b):
        c = [Fraction(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i + j] += x * y
        return c
    n = len(pts)
    coefs = [Fraction(0)] * n
    for i, (xi, yi) in enumerate(pts):
        num = [Fraction(1)]; den = Fraction(1)
        for j, (xj, _) in enumerate(pts):
            if i != j:
                num = polymul(num, [Fraction(-xj), Fraction(1)])
                den *= Fraction(xi - xj)
        term = [c * yi / den for c in num]
        coefs = [a + b for a, b in zip(coefs, term)]
    assert all(c.denominator == 1 for c in coefs)
    M = [int(c) for c in coefs]
    while len(M) > 1 and M[-1] == 0:
        M.pop()
    return M

def sturm_real_roots(p, lo=None, hi=None):
    def pp(q):
        g = 0
        for c in q:
            g = math.gcd(g, abs(c))
        return [c // g for c in q] if g > 1 else q
    def prem(a, b):
        a = a[:]
        if b[-1] < 0:
            b = [-v for v in b]
        db = len(b) - 1; lb = b[-1]
        while len(a) - 1 >= db and any(a):
            d = len(a) - 1 - db
            f = a[-1]
            nb = [0] * len(a)
            for i, c in enumerate(b):
                nb[d + i] = c * f
            a = [x * lb - y for x, y in zip(a, nb)]
            while a and a[-1] == 0:
                a.pop()
        return a
    chain = [p[:], [i * p[i] for i in range(1, len(p))] or [0]]
    while True:
        a, b = chain[-2], chain[-1]
        if not b or all(v == 0 for v in b):
            break
        r = pp([-v for v in prem(a, b)])
        if not r or all(v == 0 for v in r):
            break
        chain.append(r)
        if len(chain) > 80:
            break
    def sgn(q, x):
        v = sum(c * x ** i for i, c in enumerate(q))
        return (v > 0) - (v < 0)
    def variations(x):
        s = []
        for q in chain:
            v = sgn(q, x)
            if v:
                s.append(v)
        return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])
    def var_inf(pos):
        s = []
        for q in chain:
            lc = q[-1]; d = len(q) - 1
            v = lc if pos else lc * (-1) ** d
            if v:
                s.append((v > 0) - (v < 0))
        return sum(1 for i in range(len(s) - 1) if s[i] != s[i + 1])
    if lo is None:
        return var_inf(False) - var_inf(True)
    return variations(lo) - variations(hi)

def charpoly_of_sum(y_poly_hi_at_T, npoints, degmax):
    pts = []
    t = -npoints // 2
    while len(pts) < npoints:
        pts.append((t, res_hi(H_HI, y_poly_hi_at_T(t))))
        t += 1
    M = lagrange_int(pts)
    assert len(M) - 1 <= degmax
    for tt, v in pts[:3] + pts[-3:]:
        assert sum(c * tt ** i for i, c in enumerate(M)) == v
    return M

def part_B():
    say("=" * 78)
    say("[B] TOY-TOWER CERTIFICATES — real h; toy16.py f_toy = x^4+y x^2+y^2, g_toy = x^2+4y")
    disc_h = 58673
    n_y = snorm_h([0, 1])
    n_3y = snorm_h([0, -3])
    h45 = sum(HH[i] * Fraction(4, 5) ** i for i in range(7))
    h1 = sum(HH[i] for i in range(7))
    say("[B1] f_toy irreducible over K_h: x^4+y x^2+y^2 factors over K_h only as")
    say("     (x^2+px+q)(x^2-px+s) with p(s-q)=0, q+s-p^2=y, qs=y^2, which forces -3, y or -3y")
    say("     to be a square in K_h.  Computed: N(y) = %d < 0, N(-3y) = %d < 0 (squares have" % (n_y, n_3y))
    say("     non-negative norm), and h(4/5) = %s < 0 < h(1) = %s gives a real root y>0, where" % (h45, h1))
    say("     -3 < 0 is not a square.  K_f_toy is a field, degree 24.")
    assert n_y < 0 and n_3y < 0 and h45 < 0 < h1
    n_4y = snorm_h([0, -4])
    say("[B2] g_toy irreducible over K_h: N(-4y) = %d < 0, so -4y is not a square and" % n_4y)
    say("     x^2+4y has no root in K_h;  K_g_toy is a field, degree 12.")
    tr = htraces(13)
    A = [[2 * tr[a + b] for b in range(6)] for a in range(6)]
    Bg = [[-2 * TOY_U * tr[a + b + 1] for b in range(6)] for a in range(6)]
    G = [[0] * 12 for _ in range(12)]
    for a in range(6):
        for b in range(6):
            G[a][b] = A[a][b]; G[6 + a][6 + b] = Bg[a][b]
    dg_toy = pdet(G)
    comp = disc_h ** 2 * (-(4 * TOY_U) ** 6)
    say("[B3] disc(Z[y][x]/(h, x^2+4y)) = %d = Disc(h)^2 * N(-16y) = %d  (exact identity: %s);" %
        (dg_toy, comp, dg_toy == comp))
    say("     nonzero mod r: %s" % (dg_toy % R != 0))
    assert dg_toy == comp and dg_toy % R != 0
    Z6 = [0] * 6
    c = [Z6, Z6, [0, 1], Z6, [0, 0, 1]]
    ps = [None] * 7
    ps[0] = [4]
    for k in range(1, 5):
        acc = [0] * 6
        for i in range(1, k):
            acc = [x + y for x, y in zip(acc, pmul(c[i], ps[k - i]))]
        acc = [x + k * v for x, v in zip(acc, c[k])]
        ps[k] = [-v for v in acc]
    for k in range(5, 7):
        acc = [0] * 6
        for i in range(1, 5):
            acc = [x + y for x, y in zip(acc, pmul(c[i], ps[k - i]))]
        ps[k] = [-v for v in acc]
    G24 = [[0] * 24 for _ in range(24)]
    for a in range(6):
        for b in range(4):
            for cc in range(6):
                for d in range(4):
                    G24[a + 6 * b][cc + 6 * d] = trace_elem(ps[b + d], a + cc)
    d24 = pdet(G24)
    say("     disc(Z[y][x]/(h, f_toy)) = %d (%d bits); nonzero mod r: %s" %
        (d24, abs(d24).bit_length(), d24 % R != 0))
    assert d24 % R != 0
    say("     => r divides neither toy order discriminant, hence neither toy field discriminant")
    say("     — the premise the Schirokauer maps need at toy scale.")
    bad = [q for q in primes_upto(512) if d24 % q == 0 or dg_toy % q == 0]
    say("[B4] bad primes <= 512 (divide a toy order discriminant; excluded from the factor base")
    say("     and from relations): %s" % bad)
    n_real_h = sturm_real_roots(HH[:])
    n_neg_h = sturm_real_roots(HH[:], Fraction(-2), Fraction(0))
    say("[B5] exact Sturm: h has %d real roots (Cauchy bound 1+max|a_i| = 2), %d negative." %
        (n_real_h, n_neg_h))
    def f_toy_at(T):
        co = [T ** 4, -4 * T ** 3 + T * T, 6 * T * T - 2 * T + 1, -4 * T + 1, 1]
        return co[::-1]
    def g_toy_at(T):
        return [1, -2 * T + 4, T * T]
    Mf = charpoly_of_sum(f_toy_at, 34, 24)
    Mg = charpoly_of_sum(g_toy_at, 18, 12)
    n_real_f = sturm_real_roots(Mf)
    n_real_g = sturm_real_roots(Mg)
    say("     char. poly of x+y in K_f_toy: degree %d, %d real roots (K_f_toy is totally" % (len(Mf) - 1, n_real_f))
    say("       complex: x^4+y x^2+y^2 > 0 for real x, real y != 0).")
    say("     char. poly of x+y in K_g_toy: degree %d, %d real roots = 2 x (negative real roots of h)." %
        (len(Mg) - 1, n_real_g))
    assert n_real_f == 0 and n_real_g == 2 * n_neg_h
    c_f = 12 - 1
    c_g = n_real_g + (12 - n_real_g) // 2 - 1
    say("     unit ranks r1+r2-1: c_f = %d, c_g = %d  ->  %d Schirokauer columns." % (c_f, c_g, c_f + c_g))
    acc = [[0] * 9 for _ in range(5)]
    P1 = [[0, 4], [0], [1], [0], [0]]
    P2 = [[0, -3], [0], [1], [0], [0]]
    for i in range(5):
        for j in range(5):
            if i + j < 5:
                for ai, ac in enumerate(P1[i]):
                    for aj, bc in enumerate(P2[j]):
                        acc[i + j][ai + aj] += ac * bc
    ft = [[0] * 9 for _ in range(5)]
    ft[4][0] = 1; ft[2][1] = 1; ft[0][2] = 1
    ok = all(ft[k][m] - acc[k][m] == (13 if (k == 0 and m == 2) else 0)
             for k in range(5) for m in range(9))
    say("[B6] shared-root certificate: f_toy - (x^2+4y)(x^2-3y) = 13 y^2  (verified: %s)." % ok)
    assert ok
    say("     At every root of (h, x^2+4y) mod 13, f_toy vanishes: the toy target is the etale")
    say("     F_13-algebra F_13[y,x]/(h, x^2+4y), 12-dimensional over F_13, and 13 is a good prime")
    say("     (13 not in the bad list: %s).  This mirrors the production certificate [0]:" % (13 not in bad))
    say("     Res_y(h, f_toy(x^2 -> -4y)) = (13 y^2)^6 = 0 mod 13.")
    assert 13 not in bad
    comps = target_components(13)
    say("     target-algebra component degrees over F_13: %s (sum %d)" % (comps, sum(comps)))
    assert sum(comps) == 12
    say("")
    return dict(c_f=c_f, c_g=c_g, bad=bad, n_real_h=n_real_h, n_neg_h=n_neg_h)

def target_components(p):
    def pmod(a, b):
        a = [v % p for v in a]; b = [v % p for v in b]
        while a and a[-1] % p == 0:
            a.pop()
        if not a:
            return []
        db = len(b) - 1
        ib = pow(b[-1], p - 2, p)
        while a and len(a) - 1 >= db:
            d = len(a) - 1 - db
            f = a[-1] * ib % p
            for i, c in enumerate(b):
                a[d + i] = (a[d + i] - f * c) % p
            while a and a[-1] % p == 0:
                a.pop()
        return a
    def pgcd(a, b):
        while b and any(b):
            a, b = b, pmod(a, b)
        a = [v % p for v in a] if a else []
        while a and a[-1] % p == 0:
            a.pop()
        return a
    def pmulmod(a, b, m):
        c = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    c[i + j] = (c[i + j] + x * y) % p
        return pmod(c, m)
    def ppowmod(base, e, m):
        res = [1]
        base = pmod(base, m)
        while e:
            if e & 1:
                res = pmulmod(res, base, m)
            base = pmulmod(base, base, m)
            e >>= 1
        return res
    def pdiv(a, b):
        a = [v % p for v in a]
        while a and a[-1] % p == 0:
            a.pop()
        db = len(b) - 1; ib = pow(b[-1], p - 2, p)
        q = [0] * max(1, len(a) - db)
        while a and len(a) - 1 >= db:
            d = len(a) - 1 - db
            f = a[-1] * ib % p
            q[d] = f
            for i, c in enumerate(b):
                a[d + i] = (a[d + i] - f * c) % p
            while a and a[-1] % p == 0:
                a.pop()
        return q
    hbar = [v % p for v in HH]
    factors = []
    rem = hbar[:]
    for d in range(1, 7):
        while rem and len(rem) - 1 >= d:
            xpd = ppowmod([0, 1], p ** d, rem)
            t = [(xpd[i] if i < len(xpd) else 0) - (1 if i == 1 else 0)
                 for i in range(max(len(xpd), 2))]
            g = pgcd(t, rem)
            if not g or len(g) - 1 != d:
                break
            factors.append(g)
            rem = pdiv(rem, g)
    out = []
    for g in factors:
        e = len(g) - 1
        q = pmod([0, -TOY_U], g)
        val = ppowmod(q, (p ** e - 1) // 2, g)
        if val == [1]:
            out.extend([e, e])
        else:
            out.append(2 * e)
    return sorted(out)

def dickman_table(kmax=42, nt=600, prec=150):
    import decimal
    from decimal import Decimal as Dc
    with decimal.localcontext() as ctx:
        ctx.prec = prec
        zero, one = Dc(0), Dc(1)
        coefs = [None, [one] + [zero] * (nt - 1)]
        for k in range(2, kmax + 1):
            previous, c = coefs[k - 1], [zero] * nt
            for m in range(nt - 1):
                c[m + 1] = (previous[m] + m * c[m]) / (k * (m + 1))
            c[0] = previous[0] - sum(c[1:])
            coefs.append(c)
        return coefs

def dickman_eval(coefs, u):
    from decimal import Decimal as Dc
    import decimal
    decimal.getcontext().prec = 60
    if u <= 1:
        return 1.0
    k = int(u) + 1 if u > int(u) else int(u)
    if k >= len(coefs):
        return 0.0
    x = Dc(k) - Dc(repr(u))
    value = Dc(0)
    for c in reversed(coefs[k]):
        value = value * x + c
    return max(float(value), 0.0)

def part_D():
    say("=" * 78)
    say("[D] WHY THE REAL TOWER CANNOT BE TESTED AT TOY SCALE (measured)")
    rng = Random(11)
    lf, lg = [], []
    for _ in range(120):
        a = [rng.choice([-1, 0, 1]) for _ in range(6)]
        b = [rng.choice([-1, 0, 1]) for _ in range(6)]
        ap, bp = [[1]], [[1]]
        for _ in range(12):
            ap.append(pmul(ap[-1], a)); bp.append(pmul(bp[-1], b))
        f = [0] * 6
        for d, c in FT:
            t = pmul(c, pmul(ap[d], bp[12 - d]))
            f = [x + y for x, y in zip(f, t)]
        g = [x + y for x, y in zip(ap[2], pmul([0, U], bp[2]))]
        nf, ng = norm_h(f), norm_h(g)
        if nf:
            lf.append(math.log2(nf))
        if ng:
            lg.append(math.log2(ng))
    lf.sort(); lg.sort()
    say("     real tower, smallest possible relations ({−1,0,1}^12 box, 120 sampled pairs):")
    say("       log2|N_f|: min %.1f, median %.1f ;  log2|N_g|: min %.1f, median %.1f" %
        (lf[0], lf[len(lf) // 2], lg[0], lg[len(lg) // 2]))
    iv = dickman_table()
    def rho_at(u):
        return dickman_eval(iv, u)
    for u in (5.0, 8.0, 10.9, 20.0, lg[0] / 10.0):
        say("     Dickman rho(%.1f) = %.3e (the reference package's own interval power series)" % (u, rho_at(u)))
    ug = lg[0] / 10.0
    say("     at B = 2^10 the g-side smoothness exponent is u >= %.1f (measured min norm)," % ug)
    say("     so the expected number of relations in the whole %d-pair box is at most" % 3 ** 12)
    say("     %.1e : the real tower admits no toy-scale relation set.  (Production works at" %
        (3 ** 12 * rho_at(ug)))
    say("     B = 2^69 with log2|N_f|+log2|N_g| ~ 1285, u_f ~ 10.9 — the only regime where")
    say("     relations exist.  The toy rank test therefore runs on the toy16.py tower,")
    say("     disclosed as such.")
    say("")
