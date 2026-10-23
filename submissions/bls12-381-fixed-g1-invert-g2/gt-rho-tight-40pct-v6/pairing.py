#!/usr/bin/env python3
"""Generic pairing on y^2 = x^3 + 4 over the Fp12 tower, for any Tower modulus.

The ate pairing is evaluated with a plain Miller loop over L = |x| (the BLS12
parameter magnitude) on the untwisted G2 point, then raised to (p^12-1)/r.  The
code is not specialised to BLS12-381, so the same routines check the real curve in
part_a.py and the toy towers in toy.py.
"""
from tower import ZERO2

Z6 = (ZERO2, ZERO2, ZERO2)


# ---------------- field inverses ----------------
def f2inv(T, a):
    p = T.p
    n = (a[0]*a[0] + a[1]*a[1]) % p
    if n == 0:
        raise ZeroDivisionError
    ni = pow(n, p-2, p)
    return (a[0]*ni % p, (-a[1])*ni % p)


def f6inv(T, a):
    a0, a1, a2 = a
    xi = T.xi
    t0 = T.f2sqr(a0); t1 = T.f2sqr(a1); t2 = T.f2sqr(a2)
    p01 = T.f2mul(a0, a1); p02 = T.f2mul(a0, a2); p12 = T.f2mul(a1, a2)
    n0 = T.f2sub(t0, T.f2mul(xi, p12))
    n1 = T.f2sub(T.f2mul(xi, t2), p01)
    n2 = T.f2sub(t1, p02)
    c0 = T.f2mul(t0, a0)
    c1 = T.f2mul(xi, T.f2mul(t1, a1))
    c2 = T.f2mul(T.f2mul(xi, xi), T.f2mul(t2, a2))
    t = T.f2mul(xi, T.f2mul(a0, p12))
    N = T.f2sub(T.f2add(T.f2add(c0, c1), c2), T.f2add(T.f2add(t, t), t))
    Ni = f2inv(T, N)
    return (T.f2mul(n0, Ni), T.f2mul(n1, Ni), T.f2mul(n2, Ni))


def f12inv(T, a):
    c0, c1 = a
    t = T.f6sub(T.f6mul(c0, c0), T.f6muxi(T.f6mul(c1, c1)))
    ti = f6inv(T, t)
    return (T.f6mul(c0, ti), T.f6neg(T.f6mul(c1, ti)))


# ---------------- Fp12 helpers ----------------
def f12add(T, a, b):
    return (T.f6add(a[0], b[0]), T.f6add(a[1], b[1]))


def f12sub(T, a, b):
    return (T.f6sub(a[0], b[0]), T.f6sub(a[1], b[1]))


def f12dbl(T, a):
    return (T.f6add(a[0], a[0]), T.f6add(a[1], a[1]))


def f12zero(T):
    return (Z6, Z6)


def f12_from_fp(T, x):
    return (((x % T.p, 0), ZERO2, ZERO2), Z6)


def f12_from_fp2(T, a):
    return ((a, ZERO2, ZERO2), Z6)


def f6_by_v(T, a):
    return (T.f2mul(T.xi, a[2]), a[0], a[1])


def f6_by_vinv(T, a):
    xinv = f2inv(T, T.xi)
    return (a[1], a[2], T.f2mul(a[0], xinv))


def mul_by_w(T, a, k):
    """a * w^k for k in {-3..3}."""
    if k == 0:
        return a
    if k > 0:
        b = mul_by_w(T, a, k-1)
        return (f6_by_v(T, b[1]), b[0])
    b = mul_by_w(T, a, k+1)
    return (b[1], f6_by_vinv(T, b[0]))


# ---------------- E: y^2 = x^3 + 4 over Fp12, affine ----------------
def e_dbl(T, P):
    if P is None:
        return None
    x1, y1 = P
    if y1[0] == Z6 and y1[1] == Z6:
        return None
    three_x2 = T.f6add(T.f6add(T.f6mul(x1[0], x1[0]), T.f6mul(x1[0], x1[0])),
                       T.f6mul(x1[0], x1[0]))
    m = T.mul((three_x2, Z6), f12inv(T, f12dbl(T, y1)))
    x3 = f12sub(T, T.mul(m, m), f12dbl(T, x1))
    y3 = f12sub(T, T.mul(m, f12sub(T, x1, x3)), y1)
    return (x3, y3)


def e_add(T, P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if y1[0] == y2[0] and y1[1] == y2[1]:
            return e_dbl(T, P)
        return None
    m = T.mul(f12sub(T, y2, y1), f12inv(T, f12sub(T, x2, x1)))
    x3 = f12sub(T, f12sub(T, T.mul(m, m), x1), x2)
    y3 = f12sub(T, T.mul(m, f12sub(T, x1, x3)), y1)
    return (x3, y3)


def line(T, A, B, Pt):
    xa, ya = A; xb, yb = B; xt, yt = Pt
    if xa != xb:
        m = T.mul(f12sub(T, yb, ya), f12inv(T, f12sub(T, xb, xa)))
        return f12sub(T, T.mul(m, f12sub(T, xt, xa)), f12sub(T, yt, ya))
    if ya[0] == yb[0] and ya[1] == yb[1]:
        three_x2 = T.f6add(T.f6add(T.f6mul(xa[0], xa[0]), T.f6mul(xa[0], xa[0])),
                           T.f6mul(xa[0], xa[0]))
        m = T.mul((three_x2, Z6), f12inv(T, f12dbl(T, ya)))
        return f12sub(T, T.mul(m, f12sub(T, xt, xa)), f12sub(T, yt, ya))
    return f12sub(T, xt, xa)


def untwist(T, Q):
    """E'(Fp2): y^2 = x^3 + 4*xi  ->  E(Fp12): (x*w^-2, y*w^-3)."""
    x, y = Q
    return (mul_by_w(T, f12_from_fp2(T, x), -2), mul_by_w(T, f12_from_fp2(T, y), -3))


def miller(T, P, Q, loop):
    """f_{loop, psi(Q)}(P) in Fp12.  P = (fp, fp) in G1, Q = (fp2, fp2) in G2."""
    Pp = (f12_from_fp(T, P[0]), f12_from_fp(T, P[1]))
    Qp = untwist(T, Q)
    f = T.one12()
    R = Qp
    for b in bin(loop)[3:]:
        f = T.mul(T.mul(f, f), line(T, R, R, Pp))
        R = e_dbl(T, R)
        if b == '1':
            f = T.mul(f, line(T, R, Qp, Pp))
            R = e_add(T, R, Qp)
    return f


def pairing(T, P, Q, loop, r):
    return T.pow(miller(T, P, Q, loop), (T.p**12 - 1)//r)


# ---------------- G1 / G2 point arithmetic ----------------
def g1_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        m = (3*x1*x1) * pow(2*y1, p-2, p) % p
    else:
        m = (y2-y1) * pow((x2-x1) % p, p-2, p) % p
    x3 = (m*m - x1 - x2) % p
    return (x3, (m*(x1-x3) - y1) % p)


def g1_mul(P, k, p):
    R = None
    while k:
        if k & 1:
            R = g1_add(R, P, p)
        P = g1_add(P, P, p)
        k >>= 1
    return R


def g2_add(T, P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and T.f2add(y1, y2) == ZERO2:
        return None
    if P == Q:
        m = T.f2mul((3, 0), T.f2mul(T.f2sqr(x1), f2inv(T, T.f2add(y1, y1))))
    else:
        m = T.f2mul(T.f2sub(y2, y1), f2inv(T, T.f2sub(x2, x1)))
    x3 = T.f2sub(T.f2sub(T.f2sqr(m), x1), x2)
    return (x3, T.f2sub(T.f2mul(m, T.f2sub(x1, x3)), y1))


def g2_mul(T, P, k):
    R = None
    while k:
        if k & 1:
            R = g2_add(T, R, P)
        P = g2_add(T, P, P)
        k >>= 1
    return R
