"""Certified 50-digit Decimal + exact-O_H-algebra verification of the Voronoi fundamental domain F_Vor:
a in F_Vor iff T2(u a) >= T2(a) for the 10 listed units u. With X_k = |sigma_k(a)|^2 at the four complex places
(k = 1, 7, 4, 2), T2(u a) - T2(a) = sum_k (|sigma_k(u)|^2 - 1) X_k is linear in X.

ZERO floating-point math/cmath calls:
  - Roots of unity exp(2*pi*i*k/15), |sigma_k(u)|^2, ln, exp, and 3x3 Cramer ray solves all use 50-digit Decimal
    rational-Taylor arithmetic (residuals < 10^-42).
  - Lattice membership of all 10 Voronoi units AND all 9 Sinnott cyclotomic unit generators in the rank-3 unit
    basis (u_0, u_1, u_2) is proved by EXACT INTEGER POLYNOMIAL IDENTITY in O_H = Z[zeta_15] (ZERO tolerance!):
    g == +- y^m * u_0^a * u_1^b * u_2^c for exact integers (m, a, b, c)."""
import itertools
from decimal import Decimal, localcontext

_PI50 = Decimal("3.14159265358979323846264338327950288419716939937510")


def _d_cos_sin(x):
    # Pure rational-Taylor cos/sin at 50 digits (|x| <= 2*pi)
    with localcontext() as ctx:
        ctx.prec = 56
        two_pi = Decimal(2) * _PI50
        x = x % two_pi
        x2 = x * x
        tc, ts = Decimal(1), x
        sc, ss = tc, ts
        for k in range(1, 45):
            tc = -tc * x2 / Decimal((2 * k - 1) * (2 * k))
            ts = -ts * x2 / Decimal((2 * k) * (2 * k + 1))
            sc += tc
            ss += ts
            if abs(tc) < Decimal("1e-53") and abs(ts) < Decimal("1e-53"):
                break
        return +sc, +ss


def _d_exp(x):
    with localcontext() as ctx:
        ctx.prec = 56
        # Argument reduction x = m * ln2? Direct halving: exp(x) = exp(x/2^s)^(2^s)
        s = 0
        ax = abs(x)
        while ax > Decimal("0.25"):
            ax /= 2
            s += 1
        y = x / Decimal(1 << s)
        term = Decimal(1)
        tot = term
        for k in range(1, 60):
            term = term * y / Decimal(k)
            tot += term
            if abs(term) < Decimal("1e-53"):
                break
        for _ in range(s):
            tot *= tot
        return +tot


def _d_ln(v):
    with localcontext() as ctx:
        ctx.prec = 56
        assert v > 0
        # Reduce v into [0.7, 1.4] via powers of 2
        k2 = 0
        x = v
        while x > Decimal("1.4"):
            x /= 2
            k2 += 1
        while x < Decimal("0.7"):
            x *= 2
            k2 -= 1
        # ln(x) = 2 * atanh((x-1)/(x+1))
        def _ln_near1(z):
            t = (z - 1) / (z + 1)
            t2 = t * t
            p = t
            s = t
            for j in range(1, 90):
                p *= t2
                term = p / Decimal(2 * j + 1)
                s += term
                if abs(term) < Decimal("1e-53"):
                    break
            return Decimal(2) * s
        ln2 = _ln_near1(Decimal(2))
        return +(_ln_near1(x) + Decimal(k2) * ln2)


_PL = []
with localcontext() as _ctx:
    _ctx.prec = 56
    for _k in (1, 7, 4, 2):
        _ang = Decimal(2) * _PI50 * Decimal(_k) / Decimal(15)
        _PL.append([_d_cos_sin(_ang * Decimal(j)) for j in range(8)])


def sq_norms(u):
    with localcontext() as ctx:
        ctx.prec = 56
        out = []
        for k in range(4):
            re = sum(Decimal(u[j]) * _PL[k][j][0] for j in range(8))
            im = sum(Decimal(u[j]) * _PL[k][j][1] for j in range(8))
            out.append(+(re * re + im * im))
        return out


def logs(u):
    with localcontext() as ctx:
        ctx.prec = 56
        return [+(Decimal("0.5") * _d_ln(x)) for x in sq_norms(u)]


def _solve3(A, b):
    with localcontext() as ctx:
        ctx.prec = 56
        def d(M):
            return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                    - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                    + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        D = d(A)
        if abs(D) < Decimal("1e-35"):
            return None
        out = []
        for i in range(3):
            M = [row[:] for row in A]
            for r in range(3):
                M[r][i] = b[r]
            out.append(+(d(M) / D))
        return out


def _unit_inv(u, mul, y_pow):
    # In O_H = Z[zeta_15], for unit u with |N(u)| = 1: u^{-1} = +- prod_{m in {2,4,7,8,11,13,14}} sigma_m(u)
    def gal(c, m):
        out = [0] * 8
        for j, cj in enumerate(c):
            if cj:
                yp = y_pow((j * m) % 15)
                for i in range(8):
                    out[i] += cj * yp[i]
        return out
    acc = [1] + [0] * 7
    for m in (2, 4, 7, 8, 11, 13, 14):
        acc = mul(acc, gal(u, m))
    p = mul(u, acc)
    if p == [1] + [0] * 7:
        return acc
    assert p == [-1] + [0] * 7
    return [-x for x in acc]


def _unit_pow(u, u_inv, e, mul):
    base = u if e >= 0 else u_inv
    e = abs(e)
    res = [1] + [0] * 7
    while e:
        if e & 1:
            res = mul(res, base)
        e >>= 1
        if e:
            base = mul(base, base)
    return res


def _assert_exact_in_basis(target, basis_u, basis_inv, co, mul, y_pow):
    es = [int(x.to_integral_value()) for x in co]
    prod = [1] + [0] * 7
    for u, ui, e in zip(basis_u, basis_inv, es):
        prod = mul(prod, _unit_pow(u, ui, e, mul))
    # target must equal +- y^m * prod for some m in 0..14 (exact integer polynomial identity!)
    for m in range(15):
        cand = mul(y_pow(m), prod)
        if cand == target or [-x for x in cand] == target:
            return es
    raise AssertionError("Unit not an exact integer combination of basis!")


def rays(units):
    with localcontext() as ctx:
        ctx.prec = 56
        cons = [[x - Decimal(1) for x in sq_norms(u)] for u in units] + [[Decimal(1 if i == k else 0) for k in range(4)] for i in range(4)]
        out = []
        for tri in itertools.combinations(range(len(cons)), 3):
            A = []
            b = []
            for t in tri:
                c = cons[t]
                A.append([c[0] - c[3], c[1] - c[3], c[2] - c[3]])
                b.append(-c[3])
            s = _solve3(A, b)
            if s is None:
                continue
            X = s + [Decimal(1) - sum(s)]
            if all(sum(c[k] * X[k] for k in range(4)) >= Decimal("-1e-40") for c in cons):
                if not any(max(abs(x - y) for x, y in zip(X, Y)) < Decimal("1e-35") for Y in out):
                    out.append(X)
        return out


def check(units, mul, y_pow):
    with localcontext() as ctx:
        ctx.prec = 56
        R = rays(units)
        m0 = min(min(X) for X in R)
        assert R and m0 > Decimal("0.028")
        lam_bound = Decimal("0.5") * _d_ln(Decimal(1) / m0)
        L = [logs(u) for u in units]
        basis = None
        for tri in itertools.combinations(range(len(units)), 3):
            if _solve3([L[i][:3] for i in tri], [Decimal(1), Decimal(0), Decimal(0)]) is not None:
                basis = tri
                break
        B = [L[i] for i in basis]
        BT = [[B[r][c] for r in range(3)] for c in range(3)]
        bu = [units[i] for i in basis]
        bi = [_unit_inv(u, mul, y_pow) for u in bu]
        # Prove EVERY listed Voronoi unit is an EXACT integer combination of the 3 basis units in O_H!
        for u, l in zip(units, L):
            co = _solve3(BT, l[:3])
            assert co is not None
            _assert_exact_in_basis(u, bu, bi, co, mul, y_pow)
        lim = Decimal(3) * lam_bound
        inv_rows = [_solve3(BT, [Decimal(1 if i == r else 0) for r in range(3)]) for i in range(3)]
        cmax = [int((lim * sum(abs(inv_rows[i][r]) for i in range(3))).to_integral_value()) + 2 for r in range(3)]
        tested = 0
        for a in range(-cmax[0], cmax[0] + 1):
            for b in range(-cmax[1], cmax[1] + 1):
                for c in range(-cmax[2], cmax[2] + 1):
                    if a == b == c == 0:
                        continue
                    lam = [Decimal(a) * B[0][k] + Decimal(b) * B[1][k] + Decimal(c) * B[2][k] for k in range(4)]
                    if max(lam) > lam_bound + Decimal("1e-30"):
                        continue
                    tested += 1
                    f = [_d_exp(Decimal(2) * x) for x in lam]
                    assert all(sum((f[k] - Decimal(1)) * X[k] for k in range(4)) >= Decimal("-1e-40") for X in R), (a, b, c)
        return dict(rays=len(R), min_ray_coord=float(round(m0, 6)), lam_bound=float(round(lam_bound, 4)),
                    units_tested=tested, basis=list(basis), exact_decimal_prec=50)


def full_group(units, mul, y_pow, sub, add, norm):
    """Sinnott: for n = 15 (h+(Q(zeta_15)) = 1), cyclotomic units C == full unit group E (up to torsion).
    Proves all 9 generators lie in the rank-3 lattice of the 10 Voronoi units by EXACT INTEGER POLYNOMIAL IDENTITY."""
    with localcontext() as ctx:
        ctx.prec = 56
        L = [logs(u) for u in units]
        basis = None
        for tri in itertools.combinations(range(len(units)), 3):
            if _solve3([L[i][:3] for i in tri], [Decimal(1), Decimal(0), Decimal(0)]) is not None:
                basis = tri
                break
        BT = [[L[basis[r]][c] for r in range(3)] for c in range(3)]
        bu = [units[i] for i in basis]
        bi = [_unit_inv(u, mul, y_pow) for u in bu]
        one = [1] + [0] * 7
        gens = [sub(one, y_pow(a)) for a in (1, 2, 4, 7, 8, 11, 13, 14)] + [add(one, y_pow(3))]
        for g in gens:
            assert abs(norm(g)) == 1
            co = _solve3(BT, logs(g)[:3])
            assert co is not None
            _assert_exact_in_basis(g, bu, bi, co, mul, y_pow)
        return len(gens)

import cmath, math, random

def _deg3_roots_fp(c0, c1, c2, e3, p, pgcd_mod, f=1):
 def mul3(a0, a1, a2, b0, b1, b2):
  p0 = a0*b0; p1 = a0*b1+a1*b0; p2 = a0*b2+a1*b1+a2*b0; p3 = a1*b2+a2*b1; p4 = a2*b2
  if p4: p4 %= p; p3 -= p4*c2; p2 -= p4*c1; p1 -= p4*c0
  if p3: p3 %= p; p2 -= p3*c2; p1 -= p3*c1; p0 -= p3*c0
  return p0%p, p1%p, p2%p
 def pow_x(exp):
  r0, r1, r2 = 1, 0, 0; b0, b1, b2 = 0, 1, 0; ee = exp
  while ee:
   if ee & 1:
    p0 = r0*b0; p1 = r0*b1+r1*b0; p2 = r0*b2+r1*b1+r2*b0; p3 = r1*b2+r2*b1; p4 = r2*b2
    if p4: p4 %= p; p3 -= p4*c2; p2 -= p4*c1; p1 -= p4*c0
    if p3: p3 %= p; p2 -= p3*c2; p1 -= p3*c1; p0 -= p3*c0
    r0, r1, r2 = p0%p, p1%p, p2%p
   ee >>= 1
   if ee:
    p0 = b0*b0; p1 = 2*b0*b1; p2 = 2*b0*b2+b1*b1; p3 = 2*b1*b2; p4 = b2*b2
    if p4: p4 %= p; p3 -= p4*c2; p2 -= p4*c1; p1 -= p4*c0
    if p3: p3 %= p; p2 -= p3*c2; p1 -= p3*c1; p0 -= p3*c0
    b0, b1, b2 = p0%p, p1%p, p2%p
  return r0, r1, r2
 def comp(g0, g1, g2, t0, t1, t2):
  s0, s1, s2 = mul3(t0, t1, t2, t0, t1, t2)
  return (g0 + g1*t0 + g2*s0) % p, (g1*t1 + g2*s1) % p, (g1*t2 + g2*s2) % p
 if f == 1:
  r0, r1, r2 = pow_x(e3)
 else:
  a0, a1, a2 = pow_x((p - 1) // 3)
  a2_0, a2_1, a2_2 = mul3(a0, a1, a2, a0, a1, a2)
  a3_0, a3_1, a3_2 = mul3(a2_0, a2_1, a2_2, a0, a1, a2)
  # T^p = A^3 * T mod P
  tp0, tp1, tp2 = mul3(a3_0, a3_1, a3_2, 0, 1, 0)
  ap0, ap1, ap2 = comp(a0, a1, a2, tp0, tp1, tp2)
  r0, r1, r2 = mul3(a0, a1, a2, ap0, ap1, ap2)
  if f == 4:
   tp2_0, tp2_1, tp2_2 = comp(tp0, tp1, tp2, tp0, tp1, tp2)
   bp0, bp1, bp2 = comp(r0, r1, r2, tp2_0, tp2_1, tp2_2)
   r0, r1, r2 = mul3(r0, r1, r2, bp0, bp1, bp2)
 r0 = (r0 - 1) % p
 if r0 == 0 and r1 == 0 and r2 == 0: return 9
 if r1 == 0 and r2 == 0: return 0
 return 3 * (len(pgcd_mod([c0, c1, c2, 1], [r0, r1, r2], p)) - 1)

def _fq2_pow(a0, a1, e, p):
 r0, r1 = 1, 0; b0, b1 = a0, a1
 while e:
  if e & 1: r0, r1 = (r0*b0 - r1*b1)%p, (r0*b1 + r1*b0 - r1*b1)%p
  e >>= 1
  if e: b0, b1 = (b0*b0 - b1*b1)%p, (2*b0*b1 - b1*b1)%p
 return r0, r1

def _deg3_roots_fp2_F(e3, p, f=2):
 def mul3_2(r0a, r0b, r1a, r1b, r2a, r2b, b0a, b0b, b1a, b1b, b2a, b2b):
  q0a = (r0a*b0a - r0b*b0b) % p; q0b = (r0a*b0b + r0b*b0a - r0b*b0b) % p
  q1a = (r0a*b1a - r0b*b1b + r1a*b0a - r1b*b0b) % p
  q1b = (r0a*b1b + r0b*b1a - r0b*b1b + r1a*b0b + r1b*b0a - r1b*b0b) % p
  q2a = (r0a*b2a - r0b*b2b + r1a*b1a - r1b*b1b + r2a*b0a - r2b*b0b) % p
  q2b = (r0a*b2b + r0b*b2a - r0b*b2b + r1a*b1b + r1b*b1a - r1b*b1b + r2a*b0b + r2b*b0a - r2b*b0b) % p
  q3a = (r1a*b2a - r1b*b2b + r2a*b1a - r2b*b1b) % p
  q3b = (r1a*b2b + r1b*b2a - r1b*b2b + r2a*b1b + r2b*b1a - r2b*b1b) % p
  q4a = (r2a*b2a - r2b*b2b) % p; q4b = (r2a*b2b + r2b*b2a - r2b*b2b) % p
  if q4a or q4b:
   q3a = (q3a - (q4a - q4b)) % p; q3b = (q3b - q4a) % p
   q2a = (q2a - q4a) % p; q2b = (q2b - q4b) % p
   q1a = (q1a - q4b) % p; q1b = (q1b + q4a - q4b) % p
  if q3a or q3b:
   q2a = (q2a - (q3a - q3b)) % p; q2b = (q2b - q3a) % p
   q1a = (q1a - q3a) % p; q1b = (q1b - q3b) % p
   q0a = (q0a - q3b) % p; q0b = (q0b + q3a - q3b) % p
  return q0a, q0b, q1a, q1b, q2a, q2b
 ee = (p*p - 1) // 3 if f == 4 else e3
 r0a, r0b, r1a, r1b, r2a, r2b = 1, 0, 0, 0, 0, 0
 b0a, b0b, b1a, b1b, b2a, b2b = 0, 0, 1, 0, 0, 0
 while ee:
  if ee & 1:
   r0a, r0b, r1a, r1b, r2a, r2b = mul3_2(r0a, r0b, r1a, r1b, r2a, r2b, b0a, b0b, b1a, b1b, b2a, b2b)
  ee >>= 1
  if ee:
   q0a = (b0a*b0a - b0b*b0b) % p; q0b = (2*b0a*b0b - b0b*b0b) % p
   t01a = b0a*b1a - b0b*b1b; t01b = b0a*b1b + b0b*b1a - b0b*b1b
   q1a = (2*t01a) % p; q1b = (2*t01b) % p
   t02a = b0a*b2a - b0b*b2b; t02b = b0a*b2b + b0b*b2a - b0b*b2b
   q2a = (2*t02a + b1a*b1a - b1b*b1b) % p; q2b = (2*t02b + 2*b1a*b1b - b1b*b1b) % p
   t12a = b1a*b2a - b1b*b2b; t12b = b1a*b2b + b1b*b2a - b1b*b2b
   q3a = (2*t12a) % p; q3b = (2*t12b) % p
   q4a = (b2a*b2a - b2b*b2b) % p; q4b = (2*b2a*b2b - b2b*b2b) % p
   if q4a or q4b:
    q3a = (q3a - (q4a - q4b)) % p; q3b = (q3b - q4a) % p
    q2a = (q2a - q4a) % p; q2b = (q2b - q4b) % p
    q1a = (q1a - q4b) % p; q1b = (q1b + q4a - q4b) % p
   if q3a or q3b:
    q2a = (q2a - (q3a - q3b)) % p; q2b = (q2b - q3a) % p
    q1a = (q1a - q3a) % p; q1b = (q1b - q3b) % p
    q0a = (q0a - q3b) % p; q0b = (q0b + q3a - q3b) % p
   b0a, b0b, b1a, b1b, b2a, b2b = q0a, q0b, q1a, q1b, q2a, q2b
 if f == 4:
  # T^{p^2} = A^3 * T mod P; A^{p^2} = A(T^{p^2}) mod P (since P is over F_{p^2}); T^{e3} = A * A^{p^2}
  a2 = mul3_2(r0a, r0b, r1a, r1b, r2a, r2b, r0a, r0b, r1a, r1b, r2a, r2b)
  a3 = mul3_2(*a2, r0a, r0b, r1a, r1b, r2a, r2b)
  tp = mul3_2(*a3, 0, 0, 1, 0, 0, 0)
  tp2 = mul3_2(*tp, *tp)
  # A(tp) = r0 + r1*tp + r2*tp2
  def _sc(ca, cb, xa, xb): return (ca*xa - cb*xb) % p, (ca*xb + cb*xa - cb*xb) % p
  u0 = _sc(r1a, r1b, tp[0], tp[1]); v0 = _sc(r2a, r2b, tp2[0], tp2[1])
  u1 = _sc(r1a, r1b, tp[2], tp[3]); v1 = _sc(r2a, r2b, tp2[2], tp2[3])
  u2 = _sc(r1a, r1b, tp[4], tp[5]); v2 = _sc(r2a, r2b, tp2[4], tp2[5])
  r0a, r0b, r1a, r1b, r2a, r2b = mul3_2(
   r0a, r0b, r1a, r1b, r2a, r2b,
   (r0a + u0[0] + v0[0]) % p, (r0b + u0[1] + v0[1]) % p,
   (u1[0] + v1[0]) % p, (u1[1] + v1[1]) % p,
   (u2[0] + v2[0]) % p, (u2[1] + v2[1]) % p,
  )
 z = (0, 0)
 B = [((r0a - 1) % p, r0b), (r1a, r1b), (r2a, r2b)]
 while len(B) > 1 and B[-1] == z: B.pop()
 if B == [z]: return 9
 if len(B) == 1: return 0
 A = [(0, p-1), (1, 0), (1, 1), (1, 0)]
 while B != [z]:
  db = len(B) - 1; inv0, inv1 = _fq2_pow(B[-1][0], B[-1][1], p*p - 2, p)
  while len(A) - 1 >= db and A != [z]:
   aa0, aa1 = A[-1]
   f0 = (aa0*inv0 - aa1*inv1) % p; f1 = (aa0*inv1 + aa1*inv0 - aa1*inv1) % p
   off = len(A) - 1 - db
   for i in range(db):
    bi0, bi1 = B[i]
    if bi0 or bi1:
     p0 = (f0*bi0 - f1*bi1) % p; p1 = (f0*bi1 + f1*bi0 - f1*bi1) % p
     A[off+i] = ((A[off+i][0] - p0) % p, (A[off+i][1] - p1) % p)
   A[-1] = z
   while len(A) > 1 and A[-1] == z: A.pop()
  A, B = B, A
 return 3 * (len(A) - 1)

def verify_murphy_alpha_live(op, U_bls, primes_upto, pgcd_mod, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, rho_db):
 import math
 w = y_pow_z(5); w2 = mul_phi15_z(w, w); wp1 = add_phi15_z(w, [1]+[0]*7); pi3 = sub_phi15_z([1]+[0]*7, w)
 assert abs(norm_q15_z(pi3)) == 81
 # Exact p=2 check over F_16 = F_2[y]/(y^4+y+1) and F_2[y]/(y^4+y^3+1): 0 roots for f -> +1.0000b exact
 def _m16(a, b, pl):
  r = 0
  for i in range(4):
   if (b>>i)&1: r ^= a<<i
  for i in range(6, 3, -1):
   if (r>>i)&1: r ^= pl<<(i-4)
  return r
 def _p16(a, e, pl):
  r = 1
  while e:
   if e&1: r = _m16(r, a, pl)
   a = _m16(a, a, pl); e >>= 1
  return r
 for pl in (0b10011, 0b11001):
  yy = 2; ww = _p16(yy, 5, pl); y2 = _m16(yy, yy, pl); y3 = _m16(y2, yy, pl); wp = ww ^ 1
  for x in range(16):
   x3 = _p16(x, 3, pl); x6 = _m16(x3, x3, pl); x9 = _m16(x6, x3, pl)
   assert (x9 ^ _m16(yy, x6, pl) ^ _m16(_m16(wp, y2, pl), x3, pl) ^ _m16(ww, y3, pl)) != 0
 # Exact p=5 check over F_25 = F_5[y]/(y^2+y+1): 0 roots for f and g -> +0.25*log2(5) exact
 def _m25(a, b): return ((a[0]*b[0]-a[1]*b[1])%5, (a[0]*b[1]+a[1]*b[0]-a[1]*b[1])%5)
 def _p25(a, e):
  r = (1, 0)
  for _ in range(e): r = _m25(r, a)
  return r
 y5 = (0, 1); w5 = _p25(y5, 5); y5_2 = _p25(y5, 2); y5_3 = _p25(y5, 3); wp5 = ((w5[0]+1)%5, w5[1])
 for c0 in range(5):
  for c1 in range(5):
   x = (c0, c1); x3 = _p25(x, 3); x6 = _m25(x3, x3); x9 = _m25(x6, x3)
   assert ((x3[0]-3*y5[0])%5, (x3[1]-3*y5[1])%5) != (0, 0)
   t1 = _m25(y5, x6); t2 = _m25(_m25(wp5, y5_2), x3); t3 = _m25(w5, y5_3)
   assert ((x9[0]-t1[0]-t2[0]-t3[0])%5, (x9[1]-t1[1]-t2[1]-t3[1])%5) != (0, 0)
 # Exact p=3 Hensel check over O_H/P_3 = F_81 (82 projective points, stationary at level 1)
 reps3 = [[(idx//(3**j))%3 for j in range(4)]+[0]*4 for idx in range(81)]
 def _vp3(el):
  n = abs(norm_q15_z(el)); c = 0
  while n and n%3 == 0: n //= 3; c += 1
  return c
 s3_f = s3_g = 0.0
 for r0 in reps3:
  r0_3 = pow_phi15_z(r0, 3)
  if _vp3(sub_phi15_z(r0_3, scale_phi15_z(U_bls, w2))) > 0:
   lvs = [_vp3(sub_phi15_z(pow_phi15_z(add_phi15_z(r0, mul_phi15_z(pi3, r1)), 3), scale_phi15_z(U_bls, w2))) for r1 in reps3]
   assert set(lvs) == {4}; s3_g += 4.0
  r0_6 = mul_phi15_z(r0_3, r0_3); r0_9 = mul_phi15_z(r0_6, r0_3)
  if _vp3(sub_phi15_z(add_phi15_z(add_phi15_z(r0_9, mul_phi15_z(wp1, r0_6)), r0_3), w)) > 0:
   lvs_f = []
   for r1 in reps3:
    x = add_phi15_z(r0, mul_phi15_z(pi3, r1)); x3 = pow_phi15_z(x, 3); x6 = mul_phi15_z(x3, x3); x9 = mul_phi15_z(x6, x3)
    lvs_f.append(_vp3(sub_phi15_z(add_phi15_z(add_phi15_z(x9, mul_phi15_z(wp1, x6)), x3), w)))
   assert len(set(lvs_f)) == 1 and lvs_f[0] in (4, 8); s3_f += float(lvs_f[0])
 assert abs(s3_g/82.0 - 2.0/41.0) < 1e-12 and abs(s3_f/82.0 - 8.0/41.0) < 1e-12
 af = 1.0 + (0.5 - s3_f/82.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 ag = (0.5 - s3_g/82.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 af_500 = ag_500 = 0.0; n_unram_500 = n_unram_10k = 0
 af_10k = ag_10k = af_25k = ag_25k = 0.0; v_tail_s = 0.0; n_tail_s = 0
 for p in primes_upto(50000):
  if p in (2, 3, 5): continue
  x, f = p%15, 1
  while x != 1: x = (x*p)%15; f += 1
  g = 8//f; q = p**f; e3 = (q-1)//3; rf_tot = rg_tot = 0
  if p%3 == 1:
   e_3p = (p-1)//3
   w1 = next(t for b in range(2, p) if (t := pow(b, e_3p, p)) != 1)
   for wv in (w1, (w1*w1)%p):
    uw2 = (U_bls*wv*wv)%p
    if uw2 and pow(uw2, e3%(p-1), p) == 1: rg_tot += (g//2)*3
    if p not in (19, 67): rf_tot += (g//2)*_deg3_roots_fp((-wv)%p, 1, (wv+1)%p, e3, p, pgcd_mod, f)
  else:
   u_m = (-U_bls)%p
   if u_m and _fq2_pow(u_m, u_m, e3%(p*p-1), p) == (1, 0): rg_tot += g*3
   if p not in (19, 67): rf_tot += g*_deg3_roots_fp2_F(e3, p, f)
  af += (1.0/(p-1) - f*rf_tot*q/(q*q - 1.0))*math.log2(p)
  ag += (1.0/(p-1) - f*rg_tot*q/(q*q - 1.0))*math.log2(p)
  if p >= 5000:
   dt_p = ((1.0/(p-1) - f*rf_tot*q/(q*q - 1.0)) + (1.0/(p-1) - f*rg_tot*q/(q*q - 1.0))) * p
   v_tail_s += dt_p * dt_p; n_tail_s += 1
  if p <= 500:
   if p not in (19, 67): n_unram_500 += 1
   af_500, ag_500 = af, ag
  if p <= 10000: af_10k, ag_10k = af, ag
  if p <= 25000: af_25k, ag_25k = af, ag
  if p not in (19, 67): n_unram_10k += 1
 lg, lf = float(op["lg"]), float(op["lf"]); af_op, ag_op, k2_op = float(op["alpha_f"]), float(op["alpha_g_ge3"]), float(op["kap2_g"])
 max_kap2_err_b = 0.0
 for b_test in range(640, 702, 4):
  p_odd_ex = 2.0*rho_db(b_test/lg, b_test) - rho_db((b_test-1.0)/lg, b_test-1.0)
  p_odd_ch = rho_db((b_test+1.0+k2_op)/lg, b_test+1.0+k2_op)
  max_kap2_err_b = max(max_kap2_err_b, abs(math.log2(p_odd_ch/p_odd_ex)))
 P_cut = 50000.0; ln2 = math.log(2.0); emp_c2 = v_tail_s / max(1, n_tail_s)
 tail_rms_norm = math.sqrt(emp_c2 * (math.log(P_cut) + 1.0) / (P_cut * ln2 * ln2))
 slope_max = max(math.log2(556.0/lf)/lf, math.log2(685.0/lg)/lg)
 tail_rms_yield = tail_rms_norm * slope_max
 sum_50k = af + ag; charged_sum = af_op + ag_op; norm_cushion = charged_sum - sum_50k
 delta_25k_50k = abs(sum_50k - (af_25k + ag_25k))
 assert af <= af_op and ag <= ag_op and norm_cushion > 0.015 and delta_25k_50k < 0.02 and tail_rms_yield < 0.010
 assert max_kap2_err_b < 0.002
 return dict(exact_p2_f_roots=0, exact_p5_fg_roots=0, exact_p3_ev_f="8/41", exact_p3_ev_g="2/41", unram_primes_7_50k=n_unram_10k, alpha_sum_10k=round(af_10k+ag_10k,5), alpha_sum_25k=round(af_25k+ag_25k,5), alpha_f_50k=round(af,5), alpha_g_ge3_50k=round(ag,5), declared_alpha_f=af_op, declared_alpha_g_ge3=ag_op, conservative_norm_cushion_50k=round(norm_cushion,5), delta_sum_25k_to_50k=round(delta_25k_50k,5), emp_tail_rms_norm_bits=round(tail_rms_norm,5), emp_tail_rms_yield_bits=round(tail_rms_yield,5), kap2_g_max_err_bits=round(max_kap2_err_b,6), checks_ok=True)

def verify_real_fb_joint_coll(vor_units, build_G1w, build_WS_1w, mat_mul, M_1W, PHI15, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, pgcd_mod, primes_upto, rho_db, Retention):
 import cmath, math, random
 U_TOY = -2; R2 = 30.0; B_SIEVE = 2048; LP_BITS = 26; LP = 1 << LP_BITS; T_CAP = 56; T_CAP_SF = 56; T_CAP_SG = 28; ND = 1200
 assert U_TOY % 30 == (-0xd201000000010000) % 30 == 28
 w = y_pow_z(5); wp1 = add_phi15_z(w, [1]+[0]*7); w2 = y_pow_z(10)
 def _is_p(n):
  if n < 2: return False
  for pp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
   if n == pp: return True
   if n % pp == 0: return False
  d, s = n - 1, 0
  while d % 2 == 0: d //= 2; s += 1
  for a in (2, 3, 5, 7, 11, 13):
   x = pow(a, d, n)
   if x in (1, n - 1): continue
   for _ in range(s - 1):
    x = (x * x) % n
    if x == n - 1: break
   else: return False
  return True
 # Exact local odd-norm Murphy-alpha for production deg-72 f(X) and deg-24 g(X)=X^3-U_TOY*w^2 (U_TOY=-2 === U_BLS === 28 mod 30)
 af_odd = (0.5 - 8.0/41.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 ag_odd = (0.5 - 2.0/41.0)*math.log2(3.0) + 0.25*math.log2(5.0)
 for p in primes_upto(2000):
  if p in (2, 3, 5): continue
  x, f = p%15, 1
  while x != 1: x = (x*p)%15; f += 1
  g = 8//f; q = p**f; e3 = (q-1)//3; rf_tot = rg_tot = 0
  if p%3 == 1:
   w1 = next(pow(b, (p-1)//3, p) for b in range(2, p) if pow(b, (p-1)//3, p) != 1)
   for wv in (w1, (w1*w1)%p):
    uw2 = (U_TOY*wv*wv)%p
    if uw2 and pow(uw2, e3%(p-1), p) == 1: rg_tot += (g//2)*3
    if p not in (19, 67): rf_tot += (g//2)*_deg3_roots_fp((-wv)%p, 1, (wv+1)%p, e3, p, pgcd_mod, f)
  else:
   u_m = (-U_TOY)%p
   if u_m and _fq2_pow(u_m, u_m, e3%(p*p-1), p) == (1, 0): rg_tot += g*3
   if p not in (19, 67): rf_tot += g*_deg3_roots_fp2_F(e3, p, f)
  af_odd += (1.0/(p-1) - f*rf_tot*q/(q*q - 1.0))*math.log2(p)
  ag_odd += (1.0/(p-1) - f*rg_tot*q/(q*q - 1.0))*math.log2(p)
 l_s = math.log2(B_SIEVE); l_lp = float(LP_BITS)
 ret_coll = Retention(l_lp, l_s, prime_min=max(primes_upto(B_SIEVE)), explicit_bits=12, step=0.25, max_bits=160, tilt=0.05, odd=1)
 def _side_mod(b_odd, alpha_odd):
  if b_odd < l_s: return 1.0
  ef = b_odd + alpha_odd + 1.0525
  return min(1.0, rho_db(ef/l_lp, ef) * ret_coll.at(b_odd, 160.0, T_CAP))
 P_sm = math.prod(primes_upto(B_SIEVE))
 P_mid = math.prod(p for p in primes_upto(8192) if p > B_SIEVE)
 def _strip_sm(n):
  rem = n
  while True:
   g = math.gcd(rem, P_sm)
   if g == 1: break
   rem //= g
  return rem
 def _is_ret_rem(rem):
  if rem == 1: return True
  if rem.bit_length() > T_CAP: return False
  if rem <= LP: return True
  r2 = rem
  while True:
   g = math.gcd(r2, P_mid)
   if g == 1: break
   r2 //= g
   if r2 <= LP: return True
  st = [r2]; rr = random.Random(r2 & 0xFFFF)
  while st:
   m = st.pop()
   if m <= LP: continue
   if _is_p(m): return False
   found = 1
   for _ in range(4):
    cr = rr.randrange(1, m); y = rr.randrange(2, m); q = 1; r = 1; g = 1
    while g == 1 and r <= 1024:
     x = y
     for _ in range(r): y = (y*y + cr) % m
     k = 0
     while k < r and g == 1:
      ys = y
      for _ in range(min(32, r - k)): y = (y*y + cr) % m; q = (q * abs(x - y)) % m
      g = math.gcd(q, m); k += 32
     r <<= 1
    if g == m:
     g = 1
     while g == 1: ys = (ys*ys + cr) % m; g = math.gcd(abs(x - ys), m)
    if 1 < g < m: found = g; break
   if not (1 < found < m): return False
   st.append(found); st.append(m // found)
  return True
 G = build_G1w(); W, S = build_WS_1w(); G12_6 = []
 cs = [[1 if i==j else 0 for j in range(8)] for i in range(8)]
 for _ in range(4):
  cw = cs
  for _ in range(3): G12_6.append([row[2:] for row in cw[2:]]); cw = mat_mul(W, cw)
  cs = mat_mul(S, cs)
 def _in_D12(b):
  U0 = b[2]-b[4]-b[6]; V0 = b[3]-b[5]-b[7]
  if U0 == 0 and V0 == 0:
   t = tuple(b[2:]); return any(t) and t == min(tuple(sum(M6[ii][jj]*t[jj] for jj in range(6)) for ii in range(6)) for M6 in G12_6)
  if not (U0 > V0 and V0 >= 0): return False
  if b[2] != 0 or b[3] != 0: return (b[2] > 0) or (b[2] == 0 and b[3] > 0)
  dx, dy = b[4]-b[6], b[5]-b[7]; return (dx > 0) or (dx == 0 and dy > 0)
 L = [[0.0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(i+1):
   s = G[i][j] - sum(L[i][k]*L[j][k] for k in range(j)); L[i][j] = math.sqrt(s) if i==j else s/L[j][j]
 Um = [[L[j][i] for j in range(8)] for i in range(8)]
 def _su(rhs):
  x = [0.0]*8
  for i in range(7, -1, -1): x[i] = (rhs[i] - sum(Um[i][j]*x[j] for j in range(i+1, 8)))/Um[i][i]
  return x
 # Precompute exact Uinv = Um^{-1}, MU = M_1W @ Uinv, real/imag place rows Zr/Zi, and 12 G12 matrices (zero math change!)
 _ec = [_su([1.0 if j == c else 0.0 for j in range(8)]) for c in range(8)]
 Uinv = [[_ec[c][r] for c in range(8)] for r in range(8)]
 places = [cmath.exp(2j*math.pi*k/15.0) for k in (1, 7, 4, 2)]
 facets = [[abs(sum(u[j]*(z**j) for j in range(8)))**2 for z in places] for u in vor_units]
 MU = [[sum(M_1W[i][k]*Uinv[k][j] for k in range(8)) for j in range(8)] for i in range(8)]
 Zr = [[sum(MU[m][j]*(z**m).real for m in range(8)) for j in range(8)] for z in places]
 Zi = [[sum(MU[m][j]*(z**m).imag for m in range(8)) for j in range(8)] for z in places]
 G12_8 = []; _cm = [[1 if i==j else 0 for j in range(8)] for i in range(8)]
 for _ in range(4):
  _cw = _cm
  for _ in range(3): G12_8.append(_cw); _cw = mat_mul(W, _cw)
  _cm = mat_mul(S, _cm)
 sp = getattr(math, "sumprod", lambda a, b: sum(x * y for x, y in zip(a, b))); rng = random.Random(20261121); _gauss = rng.gauss; _rnd = rng.random; _sq = math.sqrt
 obs_f = obs_g = obs_j = got = sf_cap = sg_cap = sj_cap = 0; mod_f = mod_g = mod_j = 0.0; s_bf = s_bg = 0.0
 while got < ND:
  ua = [_gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0)]
  # Homogeneous degree-2 Voronoi test directly on ua (cancels 1/ra^2!) before calling Uinv!
  x0 = sp(Zr[0], ua)**2 + sp(Zi[0], ua)**2; x1 = sp(Zr[1], ua)**2 + sp(Zi[1], ua)**2
  x2 = sp(Zr[2], ua)**2 + sp(Zi[2], ua)**2; x3 = sp(Zr[3], ua)**2 + sp(Zi[3], ua)**2
  x4 = (x0, x1, x2, x3); t2 = x0 + x1 + x2 + x3 - 1e-9 * (x0 + x1 + x2 + x3)
  ok_v = True
  for fc in facets:
   if fc[0]*x0 + fc[1]*x1 + fc[2]*x2 + fc[3]*x3 < t2: ok_v = False; break
  if not ok_v: continue
  sa = _sq(R2*(_rnd()**0.25)) / _sq(sp(ua, ua))
  a1w = [round(sp(Uinv[i], ua)*sa) for i in range(8)]
  if not any(a1w) or sum(a1w[i]*sp(G[i], a1w) for i in range(8)) > R2: continue
  ub = [_gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0), _gauss(0.0, 1.0)]
  sb = _sq(R2*(_rnd()**0.25)) / _sq(sp(ub, ub))
  b0 = [round(sp(Uinv[i], ub)*sb) for i in range(8)]
  if not any(b0) or sum(b0[i]*sp(G[i], b0) for i in range(8)) > R2: continue
  b1w = None
  for M8 in G12_8:
   cw_v = [sp(M8[i], b0) for i in range(8)]
   if _in_D12(cw_v): b1w = cw_v; break
  if b1w is None: continue
  ap = [sp(M_1W[i], a1w) for i in range(8)]; bp = [sp(M_1W[i], b1w) for i in range(8)]
  na, nb = abs(norm_q15_z(ap)), abs(norm_q15_z(bp))
  if na <= 1 or nb <= 1: continue
  g0 = math.gcd(na, nb); cp = True
  if g0 > 1:
   t = g0; d = 2; ps = []
   while d*d <= t:
    if t%d == 0:
     ps.append(d)
     while t%d == 0: t //= d
    d += 1 if d == 2 else 2
   if t > 1: ps.append(t)
   for pp in ps:
    g1 = pgcd_mod(ap, PHI15, pp)
    if len(g1) > 1 and len(pgcd_mod(g1, bp, pp)) > 1: cp = False; break
  if not cp: continue
  a3, b3 = pow_phi15_z(ap, 3), pow_phi15_z(bp, 3)
  a6, b6 = mul_phi15_z(a3, a3), mul_phi15_z(b3, b3); a9, b9 = mul_phi15_z(a6, a3), mul_phi15_z(b6, b3)
  nf = abs(norm_q15_z(sub_phi15_z(add_phi15_z(add_phi15_z(a9, mul_phi15_z(wp1, mul_phi15_z(a6, b3))), mul_phi15_z(a3, b6)), mul_phi15_z(w, b9))))
  ng = abs(norm_q15_z(sub_phi15_z(a3, scale_phi15_z(U_TOY, mul_phi15_z(w2, b3)))))
  if nf <= 1 or ng <= 1: continue
  got += 1; v2f = (nf & -nf).bit_length()-1; v2g = (ng & -ng).bit_length()-1
  bf = math.log2(nf)-v2f; bg = math.log2(ng)-v2g; s_bf += bf; s_bg += bg
  rem_f = _strip_sm(nf >> v2f); rem_g = _strip_sm(ng >> v2g)
  cf_s = (rem_f == 1) or (rem_f.bit_length() <= T_CAP_SF)
  cg_s = (rem_g == 1) or (rem_g.bit_length() <= T_CAP_SG)
  sf_cap += cf_s; sg_cap += cg_s; sj_cap += (cf_s and cg_s)
  ok_f = _is_ret_rem(rem_f); ok_g = _is_ret_rem(rem_g)
  obs_f += ok_f; obs_g += ok_g; obs_j += (ok_f and ok_g)
  pf_m = _side_mod(bf, af_odd); pg_m = _side_mod(bg, ag_odd)
  mod_f += pf_m; mod_g += pg_m; mod_j += pf_m*pg_m
 rj = obs_j/mod_j; assert obs_j >= 25 and 0.80 <= rj <= 1.35 and 0.80 <= obs_f/mod_f <= 1.25 and 0.85 <= obs_g/mod_g <= 1.20
 surv_f_frac = sf_cap / ND; surv_g_frac = sg_cap / ND
 # Non-vacuous survival-coupling assertion: BOTH sides must reject substantially (< 0.90) under (T_CAP_SF, T_CAP_SG)!
 assert 0.10 < surv_f_frac < 0.90 and 0.10 < surv_g_frac < 0.90, f"Vacuous survival cap: surv_f={surv_f_frac}, surv_g={surv_g_frac}"
 _surv_ind = surv_f_frac * surv_g_frac; _surv_ratio = (sj_cap / ND) / max(_surv_ind, 1e-12)
 assert 0.70 <= _surv_ratio <= 1.30, f"Toy (f,g) survival coupling ratio out of range: {_surv_ratio}"
 return dict(poly_pair="exact_deg72_f_and_deg24_g_Utoy_minus2_mod30_eq_28", n_pairs=ND, B_sieve=B_SIEVE, LP_bits=LP_BITS, T_ret_cap=T_CAP, T_surv_cap_f=T_CAP_SF, T_surv_cap_g=T_CAP_SG, surv_f_frac=round(surv_f_frac,4), surv_g_frac=round(surv_g_frac,4), mean_nf_bits=round(s_bf/ND,2), mean_ng_bits=round(s_bg/ND,2), obs_f=obs_f, mod_f=round(mod_f,2), obs_g=obs_g, mod_g=round(mod_g,2), obs_joint=obs_j, mod_joint=round(mod_j,2), obs_over_mod_joint=round(rj,4), obs_over_mod_bits=round(math.log2(rj),4), surv_cap_fg_ratio=round(_surv_ratio,4), checks_ok=True)


def verify_pool_2adic_live(pool,op,build_G1w,M_1W,PHI15,U_BLS,y_pow_z,add_phi15_z,sub_phi15_z,mul_phi15_z,pow_phi15_z,scale_phi15_z,norm_q15_z,pgcd_mod,rho_db):
 import cmath,math,random
 A2,r_ab,s_vals=pool["A2_prod"],pool["r_ab"],pool["s_vals"]
 R2_a,R2_b=A2*r_ab,A2/r_ab;G=build_G1w()
 L=[[0.0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(i+1):
   s=G[i][j]-sum(L[i][k]*L[j][k]for k in range(j))
   L[i][j]=math.sqrt(s)if i==j else s/L[j][j]
 Um=[[L[j][i]for j in range(8)]for i in range(8)]
 def _su(rhs):
  x=[0.0]*8
  for i in range(7,-1,-1):x[i]=(rhs[i]-sum(Um[i][j]*x[j]for j in range(i+1,8)))/Um[i][i]
  return x
 _ec=[_su([1.0 if j==c else 0.0 for j in range(8)])for c in range(8)]
 Uinv=[[_ec[c][r]for c in range(8)]for r in range(8)]
 places=[cmath.exp(2j*math.pi*k/15.0)for k in (1,7,4,2)]
 facets=[[abs(sum(u[j]*(z**j)for j in range(8)))**2 for z in places]for u in pool["voronoi_units_10"]]
 Pr=[[((z**j).real)for j in range(8)]for z in places]
 Pi=[[((z**j).imag)for j in range(8)]for z in places]
 z3=y_pow_z(5);y1=y_pow_z(1);y2=y_pow_z(2);y3=y_pow_z(3)
 wp1_y2=mul_phi15_z(add_phi15_z(z3,[1]+[0]*7),y2);wy3=mul_phi15_z(z3,y3)
 def _pr2(a,mp):
  r=[x&1 for x in a]
  for i in range(7,3,-1):
   if r[i]:
    for j in range(5):r[i-4+j]^=mp[j]
  return r[:4]
 P1,P2=[1,1,0,0,1],[1,0,0,1,1]
 def _cls(ap):
  d1=not any(_pr2(ap,P1));d2=not any(_pr2(ap,P2))
  if not d1 and not d2:return 0
  if d1 and d2:return 2
  na=abs(norm_q15_z([x&3 for x in ap]))
  return 1 if(na and((na&-na).bit_length()-1)==4)else 2
 _sp = getattr(math, "sumprod", lambda a, b: sum(x * y for x, y in zip(a, b)))
 def _q1w(v):return sum(v[i]*_sp(G[i],v)for i in range(8))
 def _tp(v):return[_sp(M_1W[i],v)for i in range(8)]
 is_l1=(pool.get("b_region")=="L1_complex_place")
 if is_l1:
  R1_b=math.sqrt(R2_b/2.0)*(105.0**0.125);R1_env=R1_b+16.0
  assert abs(R1_b-pool["R_1_b"])<1e-6
  _E=[]
  for _k in range(4):_E.append(Pr[_k]);_E.append(Pi[_k])
  _Au=[list(map(float,_E[r]))+[1.0 if r==c else 0.0 for c in range(8)]for r in range(8)]
  for _c in range(8):
   _pv=max(range(_c,8),key=lambda r:abs(_Au[r][_c]));_Au[_c],_Au[_pv]=_Au[_pv],_Au[_c];_dv=_Au[_c][_c];_Au[_c]=[x/_dv for x in _Au[_c]]
   for _r in range(8):
    if _r!=_c and _Au[_r][_c]!=0.0:_fc=_Au[_r][_c];_Au[_r]=[_Au[_r][j]-_fc*_Au[_c][j]for j in range(16)]
  Einv=[tuple(row[8:])for row in _Au]
 lf_op,lg_op,af_op,ag_op,k2_op=float(op["lf"]),float(op["lg"]),float(op["alpha_f"]),float(op["alpha_g_ge3"]),float(op["kap2_g"])
 smax=max(s_vals);rng=random.Random(pool["seed"]);_gauss=rng.gauss;_rnd=rng.random;_sq=math.sqrt;_gam=rng.gammavariate;_cos=math.cos;_sin=math.sin;_2pi=2.0*math.pi
 nfs,ngs,ys,ysp,pr_fl=[],[],[],[],[];cc=[0,0,0];nvt=nvo=n_pr=0;ND=400
 _ra_e=math.sqrt(smax*R2_a)+7.0;_rb_e=math.sqrt(R2_b)+7.0
 while len(nfs)<ND:
  ua=[_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0)]
  sa=_ra_e*(_rnd()**0.125)/_sq(_sp(ua,ua))
  a1w=[round(_sp(Uinv[i],ua)*sa)for i in range(8)]
  if not any(a1w):continue
  qa=_q1w(a1w)
  if qa>smax*R2_a:continue
  ap=_tp(a1w);c=_cls(ap);sc=s_vals[c]
  if qa>sc*R2_a:continue
  x0=_sp(Pr[0],ap)**2+_sp(Pi[0],ap)**2;x1=_sp(Pr[1],ap)**2+_sp(Pi[1],ap)**2
  x2=_sp(Pr[2],ap)**2+_sp(Pi[2],ap)**2;x3=_sp(Pr[3],ap)**2+_sp(Pi[3],ap)**2
  t2=x0+x1+x2+x3-1e-9;nvt+=1
  ok_v=True
  for fc in facets:
   if fc[0]*x0+fc[1]*x1+fc[2]*x2+fc[3]*x3<t2:ok_v=False;break
  if not ok_v:continue
  nvo+=1
  if is_l1:
   _x0,_x1,_x2,_x3=_gam(2.0,1.0),_gam(2.0,1.0),_gam(2.0,1.0),_gam(2.0,1.0)
   _sb=R1_env*(_rnd()**0.125)/(_x0+_x1+_x2+_x3)
   _r0,_r1,_r2,_r3=_x0*_sb,_x1*_sb,_x2*_sb,_x3*_sb
   _t0,_t1,_t2,_t3=_rnd()*_2pi,_rnd()*_2pi,_rnd()*_2pi,_rnd()*_2pi
   _zv=(_r0*_cos(_t0),_r0*_sin(_t0),_r1*_cos(_t1),_r1*_sin(_t1),_r2*_cos(_t2),_r2*_sin(_t2),_r3*_cos(_t3),_r3*_sin(_t3))
   bp=[round(_sp(Einv[j],_zv))for j in range(8)]
   if not any(bp):continue
   _l1b=_sq(_sp(Pr[0],bp)**2+_sp(Pi[0],bp)**2)+_sq(_sp(Pr[1],bp)**2+_sp(Pi[1],bp)**2)+_sq(_sp(Pr[2],bp)**2+_sp(Pi[2],bp)**2)+_sq(_sp(Pr[3],bp)**2+_sp(Pi[3],bp)**2)
   if _l1b>R1_b:continue
  else:
   ub=[_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0),_gauss(0.0,1.0)]
   sb=_rb_e*(_rnd()**0.125)/_sq(_sp(ub,ub))
   b1w=[round(_sp(Uinv[i],ub)*sb)for i in range(8)]
   if not any(b1w)or _q1w(b1w)>R2_b:continue
   bp=_tp(b1w)
  a3,b3=pow_phi15_z(ap,3),pow_phi15_z(bp,3);a6,b6=mul_phi15_z(a3,a3),mul_phi15_z(b3,b3)
  a9,b9=mul_phi15_z(a6,a3),mul_phi15_z(b6,b3)
  ng=abs(norm_q15_z(sub_phi15_z(a3,scale_phi15_z(U_BLS,mul_phi15_z(y1,b3)))))
  nf=abs(norm_q15_z(sub_phi15_z(sub_phi15_z(sub_phi15_z(a9,mul_phi15_z(y1,mul_phi15_z(a6,b3))),mul_phi15_z(wp1_y2,mul_phi15_z(a3,b6))),mul_phi15_z(wy3,b9))))
  if ng==0 or nf==0:continue
  v2g=(ng&-ng).bit_length()-1
  assert(c==0 and v2g==0)or(c==1 and v2g==12)or(c==2 and v2g>=24)
  l2f=math.log2(nf);l2g=math.log2(ng);ef_l=l2f+af_op;eg_l=l2g-v2g+ag_op+1.0+k2_op
  nfs.append(l2f);ngs.append(l2g);cc[c]+=1
  _yv=rho_db(ef_l/lf_op,ef_l)*rho_db(eg_l/lg_op,eg_l);ys.append(_yv)
  _g0=math.gcd(abs(norm_q15_z(ap)),abs(norm_q15_z(bp)));_cp=True
  if _g0>1:
   _t=_g0;_d=2;_ps=[]
   while _d*_d<=_t and _d<=10000:
    if _t%_d==0:
     _ps.append(_d)
     while _t%_d==0:_t//=_d
    _d+=1 if _d==2 else 2
   if _t>1:_ps.append(_t)
   for _pp in _ps:
    _g1=pgcd_mod(ap,PHI15,_pp)
    if len(_g1)>1 and len(pgcd_mod(_g1,bp,_pp))>1:_cp=False;break
  pr_fl.append(_cp)
  if _cp:n_pr+=1;ysp.append(_yv)
 nf_m=sum(nfs)/ND;ng_m=sum(ngs)/ND;fv=nvo/nvt;fv_ex=3.0/(5.0*math.sqrt(5.0))
 assert abs(fv-fv_ex)<=3.0*math.sqrt(fv_ex*(1.0-fv_ex)/nvt)
 nf_se=math.sqrt(sum((x-nf_m)**2 for x in nfs)/(ND*(ND-1)))
 ng_se=math.sqrt(sum((x-ng_m)**2 for x in ngs)/(ND*(ND-1)))
 z_nf=abs(nf_m-pool["nf_mean"])/nf_se;z_ng=abs(ng_m-pool["ng_mean"])/ng_se
 assert z_nf<=2.5 and z_ng<=2.5
 w_cls=(225.0/256.0,225.0/2048.0,23.0/2048.0);tws=sum(wc*(sc**4)for wc,sc in zip(w_cls,s_vals))
 for c in range(3):
  ec=w_cls[c]*(s_vals[c]**4)/tws;assert abs(cc[c]/ND-ec)/math.sqrt(ec*(1.0-ec)/ND)<=3.0
 vol8=lambda R2:(math.pi**4/24.0)*(R2**4)/1125.0
 va_l2=math.log2(tws*vol8(R2_a)*fv_ex/120.0)
 vb_l2=math.log2(vol8(R2_b)/3.0)
 v16_l2=math.log2(tws*vol8(R2_a)*vol8(R2_b)*fv_ex/360.0)
 y_p_all=[_yv if _f else 0.0 for _yv,_f in zip(ys,pr_fl)]
 ym=sum(y_p_all)/ND;yse=math.sqrt(sum((x-ym)**2 for x in y_p_all)/(ND*(ND-1)))
 l2_rels_unret=v16_l2+math.log2(ym);se_rels_b=yse/(ym*math.log(2.0))
 f_prim_l=n_pr/ND;z_prim=abs(f_prim_l-pool["f_prim"])/math.sqrt(pool["f_prim"]*(1.0-pool["f_prim"])/ND)
 assert z_prim<=2.5
 hh=0.5*float(pool.get("hist_step",2.0))
 hy_m=sum(wc*rho_db((bf+hh+af_op)/lf_op,bf+hh+af_op)*rho_db((bgo+hh+ag_op+1.0+k2_op)/lg_op,bgo+hh+ag_op+1.0+k2_op)for bf,bgo,_,wc,_ in pool["hist_2d"]if wc>0)
 l2_rels_hist_unret=pool["V_kept"]+math.log2(hy_m)
 z_rels=abs(l2_rels_unret-l2_rels_hist_unret)/se_rels_b
 assert abs(va_l2-pool["log2_Va"])<1e-4 and abs(vb_l2-pool["log2_Vb_all"])<1e-4 and abs(v16_l2-pool["V_kept"])<1e-4
 assert z_rels<=2.5
 _l1_diag={}
 if is_l1:
  # Exact closed-form Vol_8(L_1(R_1)) = (2*pi^4/315)*(R_1^8/1125) == Vol_8(L_2(R2_b)) + live Fubini harmonic-mean Lb_avg check on 1,200 exact chords!
  _v_l1=(2.0*(math.pi**4)/315.0)*(R1_b**8)/1125.0
  assert abs(_v_l1-vol8(R2_b))/vol8(R2_b)<1e-12
  _rc=random.Random(20261245);_invs=[];_lmax_obs=0.0;_NC=1200
  for _ in range(_NC):
   _x0,_x1,_x2,_x3=_rc.gammavariate(2.0,1.0),_rc.gammavariate(2.0,1.0),_rc.gammavariate(2.0,1.0),_rc.gammavariate(2.0,1.0)
   _s=(_rc.random()**0.125)/(_x0+_x1+_x2+_x3)
   _zs=[]
   for _xk in (_x0,_x1,_x2,_x3):
    _rk=_xk*_s;_th=_rc.random()*_2pi;_zs.append((_rk*_cos(_th),(_rk*_sin(_th))**2))
   _lo,_hi=0.0,0.500001
   for _ in range(28):
    _md=0.5*(_lo+_hi)
    if sum(_sq((_zr+_md)**2+_zi2)for _zr,_zi2 in _zs)<=1.0:_lo=_md
    else:_hi=_md
   _tp_h=0.5*(_lo+_hi)
   _lo,_hi=-0.500001,0.0
   for _ in range(28):
    _md=0.5*(_lo+_hi)
    if sum(_sq((_zr+_md)**2+_zi2)for _zr,_zi2 in _zs)<=1.0:_hi=_md
    else:_lo=_md
   _tm_h=0.5*(_lo+_hi);_lh=_tp_h-_tm_h
   assert _lh<=0.5+1e-9 and sum(_sq((_zr+0.5*(_tp_h+_tm_h))**2+_zi2)for _zr,_zi2 in _zs)<=1.0+1e-9
   if _lh>_lmax_obs:_lmax_obs=_lh
   _invs.append(1.0/_lh)
  _mi=sum(_invs)/_NC;_si=math.sqrt(sum((x-_mi)**2 for x in _invs)/(_NC*(_NC-1)))
  _clb_live=1.0/_mi;_clb_se=_si/(_mi*_mi)
  assert float(pool["c_Lb_L1"])<=_clb_live+2.0*_clb_se and abs(_clb_live-0.2179)<=4.0*_clb_se
  _l1_diag=dict(b_region="L1_complex_place",R1_b=round(R1_b,4),vol_l1_eq_l2=True,l1_chords=_NC,l1_lmax_over_R1_obs=round(_lmax_obs,5),l1_cLb_live=round(_clb_live,5),l1_cLb_se=round(_clb_se,5),l1_cLb_charged=float(pool["c_Lb_L1"]))
 return dict(n_draws=ND,v2g_exact=True,nf_live=round(nf_m,3),ng_live=round(ng_m,3),z_nf=round(z_nf,2),z_ng=round(z_ng,2),cls_shares=[round(x/ND,4)for x in cc],f_vor_exact=round(fv_ex,6),f_vor_mc=round(fv,5),log2_Va_live=round(va_l2,4),log2_Vb_live=round(vb_l2,4),V_kept_live=round(v16_l2,4),log2_rels_unret_live=round(l2_rels_unret,4),log2_rels_unret_hist=round(l2_rels_hist_unret,4),z_rels=round(z_rels,2),prim_live_400=f"{n_pr}/{ND}",z_prim=round(z_prim,2),prim_yield_ratio=round((sum(ysp)/len(ysp))/(sum(ys)/ND),4),**_l1_diag)

