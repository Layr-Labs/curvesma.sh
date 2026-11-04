"""Toy end-to-end check of the stored six-leaf birthday construction.

Runs the exact construction of proof.md on a small prime-order group:
  * verifies the automorphism phi(x,y) = (beta x, y) acts as multiplication by lambda,
  * verifies the one-inversion six-leaf x-formulas against direct affine addition,
  * runs the birthday table, resolves a collision by the two linear equations and
    recovers the discrete logarithm, checking it against [x]P = Q.

Input: JSON on stdin, {"seed": int, "trials": int}. Output: JSON on stdout.
This is evidence at toy scale only; it says nothing about production parameters.
"""
import json
import random
import sys

# ---------------------------------------------------------------- toy curve
P_MOD = 127             # p = 1 mod 3 and y^2 = x^3 + 7 has prime order 127
A_CURVE, B_CURVE = 0, 7


def inv(a, m):
    return pow(a % m, -1, m)


def add(P, Q, m=P_MOD):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % m == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + A_CURVE) * inv(2 * y1, m) % m
    else:
        lam = (y2 - y1) * inv(x2 - x1, m) % m
    x3 = (lam * lam - x1 - x2) % m
    return (x3, (lam * (x1 - x3) - y1) % m)


def mul(k, P, m=P_MOD):
    R, Q = None, P
    while k:
        if k & 1:
            R = add(R, Q, m)
        Q = add(Q, Q, m)
        k >>= 1
    return R


def on_curve(P, m=P_MOD):
    if P is None:
        return True
    x, y = P
    return (y * y - (x * x * x + A_CURVE * x + B_CURVE)) % m == 0


def curve_order(m=P_MOD):
    n = 1
    for x in range(m):
        rhs = (x * x * x + A_CURVE * x + B_CURVE) % m
        for y in range(m):
            if y * y % m == rhs:
                n += 1
    return n


def cube_root_unity(m=P_MOD):
    for b in range(2, m):
        if (b * b + b + 1) % m == 0:
            return b
    raise AssertionError("no cube root of unity")


def point_order(P, m=P_MOD):
    o = 1
    Q = P
    while Q is not None:
        Q = add(Q, P, m)
        o += 1
        if o > m * 3:
            break
    return o


def find_generator(m=P_MOD):
    """A point of the full group order, so that phi acts on <P> as a scalar."""
    order = curve_order(m)
    for x in range(1, m):
        rhs = (x * x * x + A_CURVE * x + B_CURVE) % m
        for y in range(m):
            if y * y % m == rhs and point_order((x, y), m) == order:
                return (x, y), order
    raise AssertionError("no point of full order; change the toy prime")


def main():
    cfg = json.load(sys.stdin)
    rng = random.Random(int(cfg.get("seed", 1)))
    trials = int(cfg.get("trials", 3))

    G, N = find_generator()
    beta = cube_root_unity()

    def phi(P):
        return None if P is None else ((beta * P[0]) % P_MOD, P[1])

    # phi acts as multiplication by some lambda; recover it
    lam = None
    for k in range(2, N):
        if mul(k, G) == phi(G):
            lam = k
            break
    assert lam is not None and pow(lam, 3, N) == 1 and lam != 1

    ok_formula = ok_search = 0
    for _ in range(trials):
        x = rng.randrange(1, N)
        Q = mul(x, G)

        # --- one record: six leaves, formula versus direct addition
        a, b, c = (rng.randrange(1, N) for _ in range(3))
        X = add(mul(a, G), mul(b, Q))
        Y = mul(c, G)
        xX, yX = X
        xY, yY = Y
        bt = [1, beta, beta * beta % P_MOD]
        d = [(xY - bt[t] * xX) % P_MOD for t in range(3)]
        if 0 in d:
            continue          # degenerate record (probability ~1/r on the real curve)
        D = d[0] * d[1] % P_MOD * d[2] % P_MOD
        Di = inv(D, P_MOD)
        u = [(yY - yX) % P_MOD, (yY + yX) % P_MOD]
        direct = []
        for t in range(3):
            Pt = ((bt[t] * xX) % P_MOD, yX)
            direct.append(add(Pt, Y))
            direct.append(add(Pt, (Y[0], (-Y[1]) % P_MOD)))
        for t in range(3):
            w = Di * d[(t + 1) % 3] % P_MOD * d[(t + 2) % 3] % P_MOD
            for s in range(2):
                xx = (u[s] * u[s] % P_MOD * w % P_MOD * w - xY - bt[t] * xX) % P_MOD
                want = direct[2 * t + s]
                assert want is None or want[0] == xx, (t, s, want, xx)
        ok_formula += 1

        # --- birthday table over a scaled-down leaf count
        leaves_needed = int(2.0 * (N ** 0.5)) + 8
        n = max(2, int((leaves_needed / 6.0) ** 0.5) + 1)
        table = {}
        found = None
        for k in range(n):
            ak, bk = rng.randrange(1, N), rng.randrange(1, N)
            Xk = add(mul(ak, G), mul(bk, Q))
            for l in range(n):
                cl = rng.randrange(1, N)
                Yl = mul(cl, G)
                xXk, yXk = Xk
                xYl, yYl = Yl
                dd = [(xYl - bt[t] * xXk) % P_MOD for t in range(3)]
                if 0 in dd:
                    continue      # degenerate record, charged as a failure row in the real ledger
                DD = dd[0] * dd[1] % P_MOD * dd[2] % P_MOD
                DDi = inv(DD, P_MOD)
                uu = [(yYl - yXk) % P_MOD, (yYl + yXk) % P_MOD]
                for t in range(3):
                    w = DDi * dd[(t + 1) % 3] % P_MOD * dd[(t + 2) % 3] % P_MOD
                    for s in range(2):
                        xx = (uu[s] * uu[s] % P_MOD * w % P_MOD * w - xYl - bt[t] * xXk) % P_MOD
                        if xx in table:
                            found = (table[xx], (t, s, k, l, ak, bk, cl))
                            break
                        table[xx] = (t, s, k, l, ak, bk, cl)
                    if found:
                        break
                if found:
                    break
            if found:
                break
        assert found, "toy run found no collision; increase the trial budget"
        (t1, s1, k1, l1, ak1, bk1, cl1), (t2, s2, k2, l2, ak2, bk2, cl2) = found
        # leaf i is [lam^t a + sg c]G + [lam^t b]Q  ->  m_i + n_i x
        def coeffs(t, s, ak, bk, cl):
            sgn = 1 if s == 0 else -1
            return ((pow(lam, t, N) * ak + sgn * cl) % N,
                    (pow(lam, t, N) * bk) % N)

        m1, n1 = coeffs(t1, s1, ak1, bk1, cl1)
        m2, n2 = coeffs(t2, s2, ak2, bk2, cl2)
        recovered = None
        for num, den in ((m2 - m1, n1 - n2), (-(m1 + m2), n1 + n2)):
            den %= N
            if den == 0:
                continue
            cand = (num % N) * inv(den, N) % N
            if mul(cand, G) == Q:
                recovered = cand
                break
        assert recovered == x, ("toy recovery failed", recovered, x)
        ok_search += 1

    print(json.dumps({
        "curve": "y^2 = x^3 + 7 over F_%d" % P_MOD,
        "group_order": N,
        "cube_root_beta": beta,
        "lambda": lam,
        "formula_checks": ok_formula,
        "toy_searches": ok_search,
        "result": "PASS",
    }, indent=1))


if __name__ == "__main__":
    main()





