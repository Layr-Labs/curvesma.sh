#!/usr/bin/env python3
"""ec.py -- REAL prime-field elliptic-curve arithmetic for the v36 kernel, the
lifting stage and the ladder.  Two curves are instantiated, both genuine:

  FIXTURE  y^2 = x^3 + 4 over p37 = 68719499327 (37-bit prime, j=0 supersingular,
           group order p+1 = 192*357914059 verified by [order]G == O in
           census_ec, prime-order subgroup r37 = 357914059).
  BLS381   y^2 = x^3 + 4 over the real BLS12-381 base-field prime (381 bits,
           p % 4 == 3, asserted here).

Affine add/double (Fermat inversion) is the reference implementation.  The lifting
walk uses Jacobian coordinates plus a mixed add and an INVERSION-FREE equality test
(X == xT*Z^2 and Y == yT*Z^3), which is what makes 2^20 candidate checks affordable.
selfcheck() asserts the Jacobian path against the affine path on both curves before
anything is priced with it.
"""

P37 = 68719499327
R37 = 357914059
ORDER37 = P37 + 1
P381 = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
B = 4


def _mr(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _fix_p():
    """the constant above must be the 381-bit prime with p % 4 == 3; refuse otherwise"""
    global P381
    assert P381.bit_length() == 381, P381.bit_length()
    assert P381 % 4 == 3
    assert _mr(P381), "BLS12-381 base field prime is not prime"


_fix_p()


class Curve:
    def __init__(self, p, b=B):
        self.p, self.b = p, b

    # ---- affine (reference) ----
    def on(self, Pt):
        if Pt is None:
            return True
        x, y = Pt
        return (y * y - x * x * x - self.b) % self.p == 0

    def neg(self, Pt):
        return None if Pt is None else (Pt[0], (-Pt[1]) % self.p)

    def add(self, A, Bp):
        p = self.p
        if A is None:
            return Bp
        if Bp is None:
            return A
        x1, y1 = A
        x2, y2 = Bp
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = (3 * x1 * x1) * pow(2 * y1, p - 2, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def mul(self, k, A):
        R, Q = None, A
        while k:
            if k & 1:
                R = Q if R is None else self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def sqrt(self, v):
        p = self.p
        v %= p
        if v == 0:
            return 0
        if pow(v, (p - 1) // 2, p) != 1:
            return None
        if p % 4 == 3:
            y = pow(v, (p + 1) // 4, p)
            return y if y * y % p == v else None
        # Tonelli-Shanks, only needed if a curve prime were 1 mod 4
        q, s = p - 1, 0
        while q % 2 == 0:
            q //= 2
            s += 1
        z = 2
        while pow(z, (p - 1) // 2, p) != p - 1:
            z += 1
        m, c, t, r = s, pow(z, q, p), pow(v, q, p), pow(v, (q + 1) // 2, p)
        while t != 1:
            i, t2 = 0, t
            while t2 != 1:
                t2 = t2 * t2 % p
                i += 1
            bb = pow(c, 1 << (m - i - 1), p)
            m, c, t, r = i, bb * bb % p, t * bb * bb % p, r * bb % p
        return r

    def lift_x(self, x):
        y = self.sqrt((x * x * x + self.b) % self.p)
        return None if y is None else (x % self.p, y)

    def try_and_increment(self, start=1, tries=100):
        x = start
        for _ in range(tries):
            Pt = self.lift_x(x)
            if Pt is not None:
                return Pt
            x += 1
        return None

    # ---- Jacobian (fast path) ----
    def jdbl(self, P):
        p = self.p
        X, Y, Z = P
        if Y == 0 or Z == 0:
            return (0, 1, 0)
        XX = X * X % p
        YY = Y * Y % p
        YYYY = YY * YY % p
        ZZ = Z * Z % p
        S = 2 * ((X + YY) * (X + YY) - XX - YYYY) % p
        M = 3 * XX % p
        T = (M * M - 2 * S) % p
        X3 = T
        Y3 = (M * (S - T) - 8 * YYYY) % p
        Z3 = ((Y + Z) * (Y + Z) - YY - ZZ) % p
        return (X3, Y3, Z3)

    def jadd_mixed(self, P, Q):
        """P Jacobian, Q affine (Z2 = 1)"""
        p = self.p
        if Q is None:
            return P
        X1, Y1, Z1 = P
        if Z1 == 0:
            return (Q[0], Q[1], 1)
        x2, y2 = Q
        Z1Z1 = Z1 * Z1 % p
        U2 = x2 * Z1Z1 % p
        S2 = y2 * Z1 * Z1Z1 % p
        H = (U2 - X1) % p
        R = (S2 - Y1) % p
        if H == 0:
            if R == 0:
                return self.jdbl(P)
            return (0, 1, 0)
        HH = H * H % p
        HHH = H * HH % p
        V = X1 * HH % p
        X3 = (R * R - HHH - 2 * V) % p
        Y3 = (R * (V - X3) - Y1 * HHH) % p
        Z3 = Z1 * H % p
        return (X3, Y3, Z3)

    def jadd(self, P, Q):
        p = self.p
        X1, Y1, Z1 = P
        X2, Y2, Z2 = Q
        if Z1 == 0:
            return Q
        if Z2 == 0:
            return P
        Z1Z1 = Z1 * Z1 % p
        Z2Z2 = Z2 * Z2 % p
        U1 = X1 * Z2Z2 % p
        U2 = X2 * Z1Z1 % p
        S1 = Y1 * Z2 * Z2Z2 % p
        S2 = Y2 * Z1 * Z1Z1 % p
        H = (U2 - U1) % p
        R = (S2 - S1) % p
        if H == 0:
            if R == 0:
                return self.jdbl(P)
            return (0, 1, 0)
        HH = H * H % p
        HHH = H * HH % p
        V = U1 * HH % p
        X3 = (R * R - HHH - 2 * V) % p
        Y3 = (R * (V - X3) - S1 * HHH) % p
        Z3 = Z1 * Z2 * H % p
        return (X3, Y3, Z3)

    def jtoaff(self, P):
        X, Y, Z = P
        if Z == 0:
            return None
        p = self.p
        zi = pow(Z, p - 2, p)
        zi2 = zi * zi % p
        return (X * zi2 % p, Y * zi2 * zi % p)

    def jeq_affine(self, P, xT, yT):
        """inversion-free test  [P] == (xT,yT)"""
        p = self.p
        X, Y, Z = P
        if Z == 0:
            return False
        Z2 = Z * Z % p
        return X % p == xT * Z2 % p and Y % p == yT * Z2 * Z % p

    def jx(self, P):
        X, _, Z = P
        if Z == 0:
            return None
        p = self.p
        zi = pow(Z, p - 2, p)
        return X * zi % p * zi % p


def selfcheck():
    for cv, name in ((Curve(P37), "fixture37"), (Curve(P381), "bls381")):
        G = cv.try_and_increment(1)
        assert G is not None and cv.on(G), name
        A = cv.mul(12345, G)
        Bp = cv.mul(6789, G)
        S = cv.add(A, Bp)
        assert cv.mul(12345 + 6789, G) == S, name
        # Jacobian path must agree with affine
        Ja = (A[0], A[1], 1)
        Jb = (Bp[0], Bp[1], 1)
        assert cv.jtoaff(cv.jadd(Ja, Jb)) == S, name
        assert cv.jtoaff(cv.jadd_mixed(Ja, Bp)) == S, name
        assert cv.jtoaff(cv.jdbl(Ja)) == cv.mul(2, A), name
        assert cv.jeq_affine(Ja, A[0], A[1]) and not cv.jeq_affine(Ja, S[0], S[1]), name
        assert cv.jx(Ja) == A[0], name
    # the fixture order relation that the recovery range rests on
    cv = Curve(P37)
    G0 = cv.try_and_increment(1)
    assert cv.mul(ORDER37, G0) is None, "fixture group order"
    G = cv.mul(ORDER37 // R37, G0)
    assert G is not None and cv.mul(R37, G) is None and cv.on(G), "fixture r-subgroup"
    return G


G37 = selfcheck()

if __name__ == "__main__":
    print("fixture37 p=%d r=%d G=%r" % (P37, R37, G37))
    print("bls381  p_bits=%d on-curve=%s" % (P381.bit_length(),
                                             Curve(P381).on(Curve(P381).try_and_increment(1))))
