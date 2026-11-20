"""Exact O_H-basis of a special-q module L = {(a, b) in O_H^2 : a = rt b mod P}, P = (q, y - s) of degree one.
O_H = Z[zeta_15] is a PID: P = (pi). Given the generator pi (certificate: N(pi) = +-q and pi(s) = 0 mod q, checked) and
a primitive vector b1 of L (the first LLL row that is primitive), write b1 = x1 (pi, 0) + y1 (rt, 1); extended Euclid in
the norm-Euclidean ring O_H gives x, y with x1 y - y1 x = 1, so (b1, x (pi, 0) + y (rt, 1)) is an O_H-basis (determinant
pi, norm q); the second vector is size-reduced against the Z-lattice {y^i b1} by Babai rounding (a unimodular change)."""
import math


def check_generator(q, s, pi, Z):
    return abs(Z.norm_q15_z(pi)) == q and sum(c * pow(s, i, q) for i, c in enumerate(pi)) % q == 0


def div_round(x, y, Z):
    """round(x / y) coefficient-wise: x / y = x * prod(conj(y)) / N(y)."""
    sig = [Z.y_pow_z((7 * j) % 15) for j in range(8)]
    def act(v, k):  # y -> y^k
        out = [0] * 8
        for j, c in enumerate(v):
            if c:
                yk = Z.y_pow_z((j * k) % 15)
                for i in range(8): out[i] += c * yk[i]
        return out
    num = list(x)
    for k in (2, 4, 7, 8, 11, 13, 14):
        num = Z.mul_phi15_z(num, act(y, k))
    n = Z.norm_q15_z(y)
    return [(2 * c + n) // (2 * n) if n > 0 else -((2 * -c + -n) // (2 * -n)) for c in num] if n > 0 else [round(c / n) for c in num]
def _exact_div(x, y, Z):
    """x / y in O_H if exact, else None."""
    t = div_round(x, y, Z)
    return t if Z.mul_phi15_z(t, y) == list(x) else None


def ext_gcd(x, y, Z):
    """(g, s, t) with s x + t y = g, Euclid with nearest-coefficient rounding and a small neighbour search."""
    import itertools
    r0, r1 = list(x), list(y); s0, s1 = [1] + [0] * 7, [0] * 8; t0, t1 = [0] * 8, [1] + [0] * 7
    nb = [list(v) for v in itertools.product((-1, 0, 1), repeat=2)]
    while any(r1):
        n1 = abs(Z.norm_q15_z(r1)); qq = div_round(r0, r1, Z)
        best = None
        for d in [[0] * 8] + [[(c if j == k else 0) for j in range(8)] for k in range(8) for c in (-1, 1)]:
            qd = [a + b for a, b in zip(qq, d)]
            rr = Z.sub_phi15_z(r0, Z.mul_phi15_z(qd, r1)); nr = abs(Z.norm_q15_z(rr))
            if best is None or nr < best[0]: best = (nr, qd, rr)
        if best[0] >= n1: raise RuntimeError("no Euclidean step")
        _, qd, rr = best
        r0, r1 = r1, rr
        s0, s1 = s1, Z.sub_phi15_z(s0, Z.mul_phi15_z(qd, s1)); t0, t1 = t1, Z.sub_phi15_z(t0, Z.mul_phi15_z(qd, t1))
    return r0, s0, t0


def obasis2(q, s, rt, b1, Z, pi):
    """O_H-basis (b1, v2) of L with b1 a given primitive vector (e.g. the first LLL row)."""
    from fractions import Fraction as Fr
    a1v, b1v = b1
    y1 = list(b1v)
    x1 = _exact_div(Z.sub_phi15_z(a1v, Z.mul_phi15_z([rt] + [0] * 7, y1)), pi, Z)
    assert x1 is not None
    g, sx, ty = ext_gcd(x1, y1, Z)          # sx x1 + ty y1 = g, g a unit when b1 is primitive
    assert abs(Z.norm_q15_z(g)) == 1
    ginv = [0] * 8
    # g^-1 = g^(ord-1) up to sign: use exact division of 1 by g
    ginv = _exact_div([1] + [0] * 7, g, Z)
    # want x1 y - y1 x = 1: y = sx g^-1, x = -ty g^-1
    yv = Z.mul_phi15_z(sx, ginv); xv = [-c for c in Z.mul_phi15_z(ty, ginv)]
    wa = Z.add_phi15_z(Z.mul_phi15_z(xv, pi), Z.mul_phi15_z(yv, [rt] + [0] * 7)); wb = yv
    # Babai size reduction of w against the Z-lattice {y^i b1}
    rows = [Z.mul_phi15_z(Z.y_pow_z(i), a1v) + Z.mul_phi15_z(Z.y_pow_z(i), b1v) for i in range(8)]
    wv = wa + wb
    for _ in range(3):
        Bs = []; nb2 = []
        for r in rows:
            v = [Fr(c) for c in r]
            for bs, n2 in zip(Bs, nb2):
                m = sum(Fr(a) * b for a, b in zip(r, bs)) / n2; v = [a - m * b for a, b in zip(v, bs)]
            Bs.append(v); nb2.append(sum(c * c for c in v))
        c = [0] * 8; cur = [Fr(x) for x in wv]
        for i in range(7, -1, -1):
            m = sum(a * b for a, b in zip(cur, Bs[i])) / nb2[i]; k = round(m)
            if k: cur = [a - k * b for a, b in zip(cur, rows[i])]
        wv = [int(x) for x in cur]
    a2v, b2v = wv[:8], wv[8:]
    assert abs(Z.norm_q15_z(Z.sub_phi15_z(Z.mul_phi15_z(a1v, b2v), Z.mul_phi15_z(a2v, b1v)))) == q
    return a1v, b1v, a2v, b2v


def solve_int(A, b):
    """Integer solution u of A u = b (A: m x n integer, m <= n) by column Hermite reduction with a unimodular
    transform; returns None when no integer solution exists."""
    m, n = len(A), len(A[0])
    H = [row[:] for row in A]; U = [[int(i == j) for j in range(n)] for i in range(n)]
    def colop(j, k, a, b_, c, d):  # (col j, col k) <- (a col j + b col k, c col j + d col k)
        for r in range(m):
            x, y = H[r][j], H[r][k]; H[r][j], H[r][k] = a * x + b_ * y, c * x + d * y
        for r in range(n):
            x, y = U[r][j], U[r][k]; U[r][j], U[r][k] = a * x + b_ * y, c * x + d * y
    piv = 0
    for r in range(m):
        for k in range(piv + 1, n):
            if H[r][k]:
                x, y = H[r][piv], H[r][k]
                # extended gcd
                g0, s0, t0, g1, s1, t1 = x, 1, 0, y, 0, 1
                while g1:
                    qq = g0 // g1; g0, g1 = g1, g0 - qq * g1; s0, s1 = s1, s0 - qq * s1; t0, t1 = t1, t0 - qq * t1
                a, bb = s0, t0; c, d = -y // g0, x // g0
                colop(piv, k, a, bb, c, d)
        if H[r][piv]: piv += 1
        if piv == n: break
    # forward substitution on the lower-triangular part
    z = [0] * n; p = 0
    for r in range(m):
        acc = b[r] - sum(H[r][j] * z[j] for j in range(p))
        if p < n and H[r][p]:
            if acc % H[r][p]: return None
            z[p] = acc // H[r][p]; p += 1
        elif acc: return None
    return [sum(U[i][j] * z[j] for j in range(n)) for i in range(n)]


def complement(x1, y1, Z):
    """x, y in O_H with x1 y - y1 x = 1 (exact integer linear algebra; None if x1, y1 are not coprime)."""
    cols = [Z.mul_phi15_z(x1, Z.y_pow_z(i)) for i in range(8)] + [[-c for c in Z.mul_phi15_z(y1, Z.y_pow_z(i))] for i in range(8)]
    A = [[cols[j][r] for j in range(16)] for r in range(8)]
    u = solve_int(A, [1] + [0] * 7)
    if u is None: return None
    yv = [0] * 8; xv = [0] * 8
    for i in range(8):
        yi = Z.y_pow_z(i); yv = [a + u[i] * b for a, b in zip(yv, yi)]; xv = [a + u[8 + i] * b for a, b in zip(xv, yi)]
    assert Z.sub_phi15_z(Z.mul_phi15_z(x1, yv), Z.mul_phi15_z(y1, xv)) == [1] + [0] * 7
    return xv, yv


M0 = 0.028378          # min ray coordinate of the F_Vor cone (unitdomain.check, audit)
_PL = [cmath.exp(2j * math.pi * k / 15) for k in (1, 7, 4, 2)] if False else None


def _t2_gram():
    import cmath as _c
    pl = [_c.exp(2j * math.pi * k / 15) for k in (1, 7, 4, 2)]
    return [[2.0 * sum((z ** i * (z ** j).conjugate()).real for z in pl) for j in range(8)] for i in range(8)]


def ideal_generator(rows, N, Z):
    """A generator of the principal ideal with Z-basis rows (8 x 8) and norm N: complete Fincke-Pohst in the T2 metric up
    to T2 <= 2 N^(1/4)/M0. Some generator is T2-minimal in its unit orbit, hence in the F_Vor cone, where every coordinate
    X_k >= M0 sum X; with prod X_k = N this gives sum X <= N^(1/4)/M0 and T2 = 2 sum X <= 2 N^(1/4)/M0."""
    G = _t2_gram()
    Q = [[sum(rows[a][i] * G[i][j] * rows[b][j] for i in range(8) for j in range(8)) for b in range(8)] for a in range(8)]
    L = [[0.0] * 8 for _ in range(8)]
    for i in range(8):
        for j in range(i + 1):
            s = Q[i][j] - sum(L[i][k] * L[j][k] for k in range(j))
            L[i][j] = math.sqrt(s) if i == j else s / L[j][j]
    bound = 2.0 * N ** 0.25 / M0 * (1 + 1e-9)
    c = [0] * 8; found = []
    def rec(i, rem):
        if found: return
        if i < 0:
            if any(c):
                v = [sum(c[a] * rows[a][j] for a in range(8)) for j in range(8)]
                if abs(Z.norm_q15_z(v)) == N: found.append(v)
            return
        cc = sum(L[j][i] * c[j] for j in range(i + 1, 8)); r = math.sqrt(max(rem, 0.0)) / L[i][i]; ce = -cc / L[i][i]
        lo, hi = math.ceil(ce - r - 1e-9), math.floor(ce + r + 1e-9)
        order = sorted(range(lo, hi + 1), key=lambda t: abs(t - ce))
        for t in order:
            d = (L[i][i] * t + cc) ** 2
            if d <= rem + 1e-9:
                c[i] = t; rec(i - 1, rem - d)
                if found: break
        c[i] = 0
    rec(7, bound)
    assert found, "no generator within the proven bound"
    return found[0]


def _hnf_rows(gens):
    """Z-basis (8 rows) of the Z-span of the given integer vectors (square HNF)."""
    import copy
    M = [list(g) for g in gens]; out = []
    for col in range(8):
        piv = [r for r in M if r[col]]
        while len(piv) > 1:
            piv.sort(key=lambda r: abs(r[col])); p = piv[0]
            for r in piv[1:]:
                qq = r[col] // p[col]
                for j in range(8): r[j] -= qq * p[j]
            piv = [r for r in M if r[col]]
        if piv:
            p = piv[0]; M.remove(p); out.append(p)
    return out


def basis(q, s, rt, b1, Z, pi):
    """O_H-basis (v1, v2) of L: v1 = b1 / content(b1), v2 = exact complement, Babai-reduced."""
    from fractions import Fraction as Fr
    a1v, b1v = b1; y1 = list(b1v)
    x1 = _exact_div(Z.sub_phi15_z(a1v, Z.mul_phi15_z([rt] + [0] * 7, y1)), pi, Z)
    assert x1 is not None
    comp = complement(x1, y1, Z)
    if comp is None:
        # content ideal (x1, y1) = (g): generator by the proven enumeration, then divide
        gens = [Z.mul_phi15_z(Z.y_pow_z(i), x1) for i in range(8)] + [Z.mul_phi15_z(Z.y_pow_z(i), y1) for i in range(8)]
        rows = _hnf_rows(gens); Nn = abs(Z.norm_q15_z([1] + [0] * 7))
        det = 1
        for i, r in enumerate(rows): det *= r[i]
        g = ideal_generator(rows, abs(det), Z)
        x1 = _exact_div(x1, g, Z); y1 = _exact_div(y1, g, Z); assert x1 is not None and y1 is not None
        a1v = Z.add_phi15_z(Z.mul_phi15_z(x1, pi), Z.mul_phi15_z([rt] + [0] * 7, y1)); b1v = y1
        comp = complement(x1, y1, Z); assert comp is not None
    xv, yv = comp
    wv = Z.add_phi15_z(Z.mul_phi15_z(xv, pi), Z.mul_phi15_z(yv, [rt] + [0] * 7)) + list(yv)
    rows = [Z.mul_phi15_z(Z.y_pow_z(i), a1v) + Z.mul_phi15_z(Z.y_pow_z(i), b1v) for i in range(8)]
    Bs = []; nb2 = []
    for r in rows:
        v = [Fr(c) for c in r]
        for bs, n2 in zip(Bs, nb2):
            mm = sum(Fr(a) * b for a, b in zip(r, bs)) / n2; v = [a - mm * b for a, b in zip(v, bs)]
        Bs.append(v); nb2.append(sum(c * c for c in v))
    cur = [Fr(x) for x in wv]
    for i in range(7, -1, -1):
        k = round(sum(a * b for a, b in zip(cur, Bs[i])) / nb2[i])
        if k: cur = [a - k * b for a, b in zip(cur, rows[i])]
    a2v, b2v = [int(x) for x in cur[:8]], [int(x) for x in cur[8:]]
    assert abs(Z.norm_q15_z(Z.sub_phi15_z(Z.mul_phi15_z(a1v, b2v), Z.mul_phi15_z(a2v, b1v)))) == q
    return a1v, b1v, a2v, b2v


def prime_generator(q, s, Z):
    """Generator of P = (q, y - s) by the proven enumeration (Z-basis: q, y^i - s^i)."""
    rows = [[q] + [0] * 7] + [[(-pow(s, i, q)) % q if j == 0 else int(j == i) for j in range(8)] for i in range(1, 8)]
    # LLL-free: the enumeration is complete for any basis; reduce first for speed (exact integer LLL of dim 8)
    return ideal_generator(_lll8(rows), q, Z)


def _lll8(B):
    from fractions import Fraction as Fr
    B = [r[:] for r in B]; n = len(B)
    def gso():
        Bs = []; mu = [[Fr(0)] * n for _ in range(n)]; nb = []
        for i in range(n):
            v = [Fr(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = sum(a * b for a, b in zip(B[i], Bs[j])) / nb[j]; v = [a - mu[i][j] * b for a, b in zip(v, Bs[j])]
            Bs.append(v); nb.append(sum(a * a for a in v))
        return mu, nb
    k = 1; mu, nb = gso()
    while k < n:
        for j in range(k - 1, -1, -1):
            t = round(mu[k][j])
            if t: B[k] = [a - t * b for a, b in zip(B[k], B[j])]; mu, nb = gso()
        if nb[k] >= (Fr(99, 100) - mu[k][k - 1] ** 2) * nb[k - 1]: k += 1
        else: B[k], B[k - 1] = B[k - 1], B[k]; mu, nb = gso(); k = max(k - 1, 1)
    return B
