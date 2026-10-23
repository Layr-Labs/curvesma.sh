"""Checks for the Karabina batched-orbit walk.

The BLS12-381 constants are the demontgomery values of blst coeffs2[] at
commit f7a88067ed98385ebc6cf41ccecf21ead706a5b7. The group checks use a small
tower and a small cyclic group, not GT.
"""

import json
import sys

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001

C1 = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAD
C2 = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAC
C3 = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAA


def check_constants():
    if pow(C1, 6, P) != 1:
        raise SystemExit("c1^6 != 1")
    if pow(C1, 3, P) == 1:
        raise SystemExit("c1 has order dividing 3")
    if pow(C1, 2, P) != C2:
        raise SystemExit("c2 != c1^2")
    if pow(C1, 3, P) != C3:
        raise SystemExit("c3 != c1^3")
    if pow(P, 6, R) != R - 1:
        raise SystemExit("p^6 != -1 mod r")
    if (R - 1) % 12 != 0:
        raise SystemExit("12 does not divide r-1")


def check_bound():
    # Certificate for log2((50400/60000)*S + 2^20) < 12572/100,
    # with S^2 >= (355/113)*n/2 >= pi*n/2.
    n = (R - 1) // 12
    ceil_arg = (355 * n + 225) // 226
    s = int(ceil_arg).bit_length()
    # integer square root
    root = 1 << ((ceil_arg.bit_length() + 1) // 2)
    while root * root > ceil_arg:
        root = (root + ceil_arg // root) // 2
    while root * root < ceil_arg:
        root += 1
    if root * root * 226 < 355 * n:
        raise SystemExit("sqrt bound failed")
    num = 50400 * root + 60000 * (1 << 20)
    den = 60000
    if num**100 >= (den**100) * (1 << 12572):
        raise SystemExit("log2 certificate failed")


def check_collision(trials):
    q = 13
    zeta = 2
    if pow(zeta, 12, q) != 1 or pow(zeta, 6, q) == 1:
        raise SystemExit("toy automorphism is not order 12")
    for i in range(trials):
        y = (2 + i) % q
        alpha1 = (5 + 3 * i) % q
        beta1 = (7 + i) % q
        delta = i % 12
        scale = pow(zeta, delta, q)
        beta2 = (beta1 * scale + 1 + (i % 5)) % q
        alpha2 = (alpha1 * scale + y * (beta1 * scale - beta2)) % q
        numer = (alpha1 * scale - alpha2) % q
        denom = (beta2 - beta1 * scale) % q
        if denom == 0:
            raise SystemExit("degenerate toy denominator")
        y_hat = (numer * pow(denom, -1, q)) % q
        if y_hat != y:
            raise SystemExit("collision equation failed")


def check_karabina(trials):
    """Compressed square matches a full cyclotomic square on a toy tower."""
    p = 19
    if p % 4 != 3:
        raise SystemExit("toy prime")

    def add2(x, y):
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)

    def sub2(x, y):
        return ((x[0] - y[0]) % p, (x[1] - y[1]) % p)

    def neg2(x):
        return ((-x[0]) % p, (-x[1]) % p)

    def mul2(x, y):
        return (
            (x[0] * y[0] - x[1] * y[1]) % p,
            (x[0] * y[1] + x[1] * y[0]) % p,
        )

    def sq2(x):
        return mul2(x, x)

    def nr(x):
        # (a+bu)(1+u) = (a-b)+(a+b)u with u^2 = -1
        return ((x[0] - x[1]) % p, (x[0] + x[1]) % p)

    def dbl(x):
        return add2(x, x)

    xi = (1, 1)

    def mul6(x, y):
        z0 = add2(mul2(x[0], y[0]), mul2(xi, add2(mul2(x[1], y[2]), mul2(x[2], y[1]))))
        z1 = add2(add2(mul2(x[0], y[1]), mul2(x[1], y[0])), mul2(xi, mul2(x[2], y[2])))
        z2 = add2(add2(mul2(x[0], y[2]), mul2(x[2], y[0])), mul2(x[1], y[1]))
        return (z0, z1, z2)

    def mul6_v(x):
        return (mul2(xi, x[2]), x[0], x[1])

    zero6 = ((0, 0), (0, 0), (0, 0))
    one6 = ((1, 0), (0, 0), (0, 0))

    def mul12(a, b):
        return (
            add6(mul6(a[0], b[0]), mul6_v(mul6(a[1], b[1]))),
            add6(mul6(a[0], b[1]), mul6(a[1], b[0])),
        )

    def add6(x, y):
        return tuple(add2(x[i], y[i]) for i in range(3))

    def sq12(a):
        return mul12(a, a)

    def pow12(a, e):
        r = (one6, zero6)
        b = a
        while e:
            if e & 1:
                r = mul12(r, b)
            b = sq12(b)
            e >>= 1
        return r

    def karabina(g1, g2, g3, g5):
        t0 = sq2(g1)
        t1 = sq2(g5)
        t5 = add2(g1, g5)
        t2 = sq2(t5)
        t3 = add2(t0, t1)
        t5 = sub2(t2, t3)
        t6 = add2(g3, g2)
        t3 = sq2(t6)
        t2 = sq2(g3)
        t6 = nr(t5)
        t5 = dbl(add2(t6, g3))
        z3 = add2(t5, t6)
        t4 = nr(t1)
        t5 = add2(t0, t4)
        t6 = sub2(t5, g2)
        t1 = sq2(g2)
        t6 = dbl(t6)
        z2 = add2(t6, t5)
        t4 = nr(t1)
        t5 = add2(t2, t4)
        t6 = sub2(t5, g1)
        t6 = dbl(t6)
        z1 = add2(t6, t5)
        t0 = add2(t2, t1)
        t5 = sub2(t3, t0)
        t6 = add2(t5, g5)
        t6 = dbl(t6)
        z5 = add2(t5, t6)
        return z1, z2, z3, z5

    phi = p**4 - p**2 + 1
    exp = (p**12 - 1) // phi
    seen = 0
    seed = 1
    while seen < trials:
        seed += 1
        a = tuple(((seed * (i + 3)) % p, (seed * (i + 5)) % p) for i in range(3))
        b = tuple(((seed * (i + 8)) % p, (seed * (i + 2)) % p) for i in range(3))
        x = pow12((a, b), exp)
        if x == (one6, zero6):
            continue
        xs = sq12(x)
        g1, g2, g3, g5 = x[0][1], x[0][2], x[1][0], x[1][2]
        z1, z2, z3, z5 = karabina(g1, g2, g3, g5)
        if (z1, z2, z3, z5) != (xs[0][1], xs[0][2], xs[1][0], xs[1][2]):
            raise SystemExit("karabina square mismatch")
        seen += 1


def main():
    raw = sys.stdin.read() if sys.stdin.isatty() is False else ""
    trials = 20
    if raw.strip():
        data = json.loads(raw)
        trials = int(data["trials"])
    if not 1 <= trials <= 50:
        raise SystemExit("trials out of range")
    check_constants()
    check_bound()
    check_collision(trials)
    check_karabina(trials)
    json.dump({"ok": True, "trials": trials}, sys.stdout)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
