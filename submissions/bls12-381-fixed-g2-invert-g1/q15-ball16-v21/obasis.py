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
