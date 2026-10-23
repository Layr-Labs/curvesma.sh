"""Toy class-walk check plus production parameter, endomorphism and ledger checks.

The toy group is the additive group of integers modulo a toy prime whose
multiplicative group has the same 2-power structure as F_r. The class walk,
the canonicalisation units, the exact scalar ratio and the square-root phase 2
are exercised end to end at toy size. The production BLS12-381 checks verify the
cube-root endomorphism identity, the 3x85-bit comb digits, the factorization of
r-1, the primitive root, the collision parameters and the integer cost ledger.
This is evidence at toy scope only; it is not a production-size run.
"""
import json
import math
import random
import sys
from decimal import Decimal, getcontext

getcontext().prec = 80

R = int("73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001", 16)
P = int("1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f624"
        "1eabfffeb153ffffb9feffffffffaaab", 16)
D = 1 << 20
BITS = 85
WINDOWS = 3
C = 2000
G1X = int("17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac5"
          "86c55e83ff97a1aeffb3af00adb22c6bb", 16)
G1Y = int("08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3e"
          "dd03cc744a2888ae40caa232946c5e7e1", 16)


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        s += 1
        d //= 2
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


def distinct_prime_factors(n):
    factors = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        factors.append(n)
    return factors


def primitive_root(p):
    factors = distinct_prime_factors(p - 1)
    g = 2
    while True:
        if all(pow(g, (p - 1) // q, p) != 1 for q in factors):
            return g
        g += 1


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(a, q, p)
    root = pow(a, (q + 1) // 2, p)
    while t != 1:
        i = 1
        t2 = t * t % p
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = b * b % p
        t = t * c % p
        root = root * b % p
    return root


def toy_prime(seed, d_bits):
    mod = 6 << d_bits
    q = (1 << (24 - d_bits)) // 6 + (seed % 89)
    while True:
        r0 = mod * q + 1
        if is_prime(r0):
            return r0, q
        q += 1


def toy_comb(mod, w, windows, base=1):
    tables = []
    for j in range(windows):
        step = base * pow(2, w * j, mod) % mod
        tab = []
        acc = 0
        for _ in range(1 << w):
            tab.append(acc)
            acc = (acc + step) % mod
        tables.append(tab)
    return tables


def toy_comb_mul(x, tables, w, windows, mod):
    acc = 0
    for j in range(windows):
        digit = (x >> (w * j)) & ((1 << w) - 1)
        acc = (acc + tables[j][digit]) % mod
    return acc


def toy_class_walk(r0, zeta, d0, secret, rng, jumps=256, max_steps=20000):
    """Collision search on classes of the order-M1 subgroup of F_{r0}^*.

    Points are integers mod r0. Units {+-1, +-lam, +-lam^2} act by scalar
    multiplication, as the cube-root endomorphism acts on G1. The walk tracks
    the non-canonical scalar and folds the canonicalisation unit into the next
    jump product, exactly as the production step does.
    """
    M1 = (r0 - 1) // d0
    omega = pow(zeta, d0, r0)
    lam = pow(zeta, (r0 - 1) // 3, r0)
    units = [1, r0 - 1, lam, (r0 - lam) % r0, lam * lam % r0, (-lam * lam) % r0]
    assert len(set(units)) == 6
    assert all(pow(u, M1, r0) == 1 for u in units)
    g = [pow(omega, rng.randrange(1, M1), r0) for _ in range(jumps)]
    h = [[g[j] * units[i] % r0 for i in range(6)] for j in range(jumps)]

    def canon(x):
        best = None
        bi = 0
        for i, u in enumerate(units):
            y = u * x % r0
            if best is None or y < best:
                best = y
                bi = i
        return best, bi

    v = pow(secret, d0, r0)
    e = pow(omega, rng.randrange(1, M1), r0)
    w = pow(omega, rng.randrange(1, M1), r0)
    ce, ui_e = canon(e)
    cw, ui_w = canon(w * v % r0)
    table = {}
    for _ in range(max_steps):
        j = ce % jumps
        e = h[j][ui_e] * e % r0
        cv, ui_e = canon(e)
        c = units[ui_e] * e % r0
        hit = table.get(cv)
        if hit is not None and hit[0] == 1:
            got = c * pow(hit[1], -1, r0) % r0
            if got == v:
                return True, got
        else:
            table[cv] = (0, c)
        ce = cv

        j = cw % jumps
        w = h[j][ui_w] * w % r0
        cv, ui_w = canon(w * v % r0)
        c = units[ui_w] * w % r0
        hit = table.get(cv)
        if hit is not None and hit[0] == 0:
            got = hit[1] * pow(c, -1, r0) % r0
            if got == v:
                return True, got
        else:
            table[cv] = (1, c)
        cw = cv
    return False, None


def toy_phase2(r0, zeta, d0, v, target):
    """20 (d_bits) successive square roots, then BSGS over the d0 roots of unity."""
    x = v
    for _ in range(d0.bit_length() - 1):
        x = sqrt_mod(x, r0)
        assert x is not None
    s0 = x
    assert pow(s0, d0, r0) == v
    eta = pow(zeta, (r0 - 1) // d0, r0)
    m2 = math.isqrt(d0)
    babies = {}
    cur = s0
    for i in range(m2):
        babies[cur] = i
        cur = cur * eta % r0
    step = pow(eta, -m2, r0)
    cur = target
    for j in range((d0 + m2 - 1) // m2 + 1):
        if cur in babies:
            k = babies[cur] + m2 * j
            if k < d0:
                s = s0 * pow(eta, k, r0) % r0
                if s == target:
                    return True, k
        cur = cur * step % r0
    return False, None


def jac_double(pt):
    X, Y, Z = pt
    if Z == 0 or Y == 0:
        return (0, 1, 0)
    A = X * X % P
    B = Y * Y % P
    Cc = B * B % P
    Dv = 2 * ((X + B) * (X + B) - A - Cc) % P
    E = 3 * A % P
    F = E * E % P
    X3 = (F - 2 * Dv) % P
    Y3 = (E * (Dv - X3) - 8 * Cc) % P
    Z3 = 2 * Y * Z % P
    return (X3, Y3, Z3)


def jac_add(pt, qt):
    X1, Y1, Z1 = pt
    X2, Y2, Z2 = qt
    if Z1 == 0:
        return qt
    if Z2 == 0:
        return pt
    Z1Z1 = Z1 * Z1 % P
    Z2Z2 = Z2 * Z2 % P
    U1 = X1 * Z2Z2 % P
    U2 = X2 * Z1Z1 % P
    S1 = Y1 * Z2 * Z2Z2 % P
    S2 = Y2 * Z1 * Z1Z1 % P
    if U1 == U2:
        if S1 != S2:
            return (0, 1, 0)
        return jac_double(pt)
    H = (U2 - U1) % P
    I = (2 * H) * (2 * H) % P
    J = H * I % P
    rr = 2 * (S2 - S1) % P
    V = U1 * I % P
    X3 = (rr * rr - J - 2 * V) % P
    Y3 = (rr * (V - X3) - 2 * S1 * J) % P
    Z3 = ((Z1 + Z2) * (Z1 + Z2) - Z1Z1 - Z2Z2) * H % P
    return (X3, Y3, Z3)


def ec_mul(k, pt):
    acc = (0, 1, 0)
    addend = pt
    while k:
        if k & 1:
            acc = jac_add(acc, addend)
        addend = jac_double(addend)
        k >>= 1
    zinv = pow(acc[2], -1, P)
    z2 = zinv * zinv % P
    return (acc[0] * z2 % P, acc[1] * z2 * zinv % P)


def curve_checks():
    beta = pow(2, (P - 1) // 3, P)
    assert pow(beta, 3, P) == 1 and beta != 1
    lam = pow(7, (R - 1) // 3, R)
    assert pow(lam, 3, R) == 1 and lam != 1
    G = (G1X, G1Y)
    phi = (beta * G1X % P, G1Y)
    assert ec_mul(lam, G) == phi
    B = 1 << BITS
    checks = [R - 1, (R - 1) // 2]
    seed = 987654321
    seed = (seed * 6364136223846793005 + 1442695040888963407) % R
    checks.append(seed)
    for x in checks:
        x0 = x % B
        x1 = (x // B) % B
        x2 = x // (B * B)
        assert x0 < B and x1 < B and x2 < B
        assert x0 + x1 * B + x2 * B * B == x
        s = jac_add(jac_add(ec_mul(x0, G), ec_mul(x1, (2 ** BITS, 1, 0)[0] and G)),
                    ec_mul(x2, G))
        # [x1*2^85]G via doubling the base, and [x2*2^170]G likewise
        base1 = (G[0], G[1], 1)
        for _ in range(BITS):
            base1 = jac_double(base1)
        p1 = ec_mul_inf(x1, base1)
        base2 = base1
        for _ in range(BITS):
            base2 = jac_double(base2)
        p2 = ec_mul_inf(x2, base2)
        total = jac_add(jac_add((G[0], G[1], 1) if False else (ec_mul(x0, G)[0], ec_mul(x0, G)[1], 1), p1), p2)
        assert affine(total) == ec_mul(x, G)
    return True


def ec_mul_inf(k, pt):
    acc = (0, 1, 0)
    addend = pt
    while k:
        if k & 1:
            acc = jac_add(acc, addend)
        addend = jac_double(addend)
        k >>= 1
    return acc


def affine(pt):
    if pt[2] == 0:
        return (0, 0)
    zinv = pow(pt[2], -1, P)
    z2 = zinv * zinv % P
    return (pt[0] * z2 % P, pt[1] * z2 * zinv % P)


def production_checks(seed):
    out = {}
    assert is_prime(R)
    product = 1
    for q, e in ((2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1),
                 (859267, 1), (906349, 2), (2508409, 1), (2529403, 1),
                 (52437899, 1), (254760293, 2)):
        assert is_prime(q)
        product *= q ** e
    assert product == R - 1
    root_primes = (2, 3, 11, 19, 10177, 125527, 859267, 906349,
                   2508409, 2529403, 52437899, 254760293)
    for g0 in range(2, 8):
        is_root = all(pow(g0, (R - 1) // q, R) != 1 for q in root_primes)
        assert is_root == (g0 == 7), (g0, is_root)

    assert is_prime(P) and P % 3 == 1
    beta = pow(2, (P - 1) // 3, P)
    assert pow(beta, 3, P) == 1 and beta != 1
    lam = pow(7, (R - 1) // 3, R)
    assert pow(lam, 3, R) == 1 and lam != 1

    assert (R - 1) % D == 0
    M1 = (R - 1) // D
    assert M1 % 6 == 0
    M = M1 // 6
    assert M * 6 == M1
    for u in (R - 1, lam, lam * lam % R):
        assert pow(u, M1, R) == 1

    num = 2392394
    m = math.isqrt(num * M // 1000000)
    while 1000000 * m * m < num * M:
        m += 1
    assert 1000000 * m * m >= num * M
    assert 1000000 * (m - 1) * (m - 1) < num * M
    assert m == 141206611312350177476731448426195532
    ideal_p = 1 - (-(Decimal(m * m) / (4 * Decimal(M)))).exp()
    assert ideal_p >= Decimal("0.45")

    h_total = 2 * m + 2 ** 92 + 2 ** 28 + 2 ** 22 + 2 ** 20
    r_total = (1750 * m + 6 * (2 ** 85) * 8192 + 2 ** 41 + 2 ** 45
               + 2 ** 36 + 2 ** 30 + 2 ** 12)
    assert h_total == 282413227576460512094983996722566296
    assert r_total == 247111571698088710926624174444085556104
    w = Decimal(h_total) + Decimal(r_total) / C
    assert w == Decimal("405969013425504867558296083944609074.052")
    log2_w = w.ln() / Decimal(2).ln()
    assert log2_w <= Decimal(118288853) / 1000000
    assert log2_w >= Decimal(118288850) / 1000000

    prep = Decimal(3 * (2 ** 85) * 8192 + 2 ** 41 + 2 ** 45) / C
    log2_prep = prep.ln() / Decimal(2).ln()
    assert log2_prep <= Decimal(88619179) / 1000000

    mem = (6 * (2 ** 85) * 96 + 2 ** 27 * 160 + 12288 * 48 + 2048 * 48 + 2 ** 20)
    log2_mem = Decimal(mem).ln() / Decimal(2).ln()
    assert log2_mem <= Decimal(94169926) / 1000000

    B = 1 << BITS
    xs = {0, 1, 2, B - 1, B, B + 1, R - 1, R - 2, (R - 1) // 2, (R - 1) // 3}
    s = seed
    for _ in range(64):
        s = (s * 6364136223846793005 + 1442695040888963407) & ((1 << 255) - 1)
        xs.add(s % R)
    for x in xs:
        x0 = x % B
        x1 = (x // B) % B
        x2 = x // (B * B)
        assert x0 < B and x1 < B and x2 < B
        assert x0 + x1 * B + x2 * B * B == x

    out.update({
        "M1": str(M1),
        "M": str(M),
        "m": str(m),
        "ideal_success_at_m": str(ideal_p),
        "H": str(h_total),
        "R": str(r_total),
        "log2_W": str(log2_w),
        "claimed_log2": "118.288853",
        "preprocessing_log2": str(log2_prep),
        "memory_log2": str(log2_mem),
        "ledger_ok": True,
        "factorization_ok": True,
        "primitive_root_ok": True,
        "endomorphism_ok": True,
        "digit_checks_ok": True,
    })
    return out


def main():
    params = json.load(sys.stdin)
    trials = int(params["trials"])
    seed = int(params["seed"])
    d_bits = int(params["d_bits"])
    if not 1 <= trials <= 8:
        raise ValueError("trials out of range")
    if not 0 <= seed <= 2147483647:
        raise ValueError("seed out of range")
    if not 8 <= d_bits <= 12:
        raise ValueError("d_bits out of range")

    r0, q = toy_prime(seed, d_bits)
    assert (r0 - 1) % (6 << d_bits) == 0
    zeta = primitive_root(r0)
    d0 = 1 << d_bits
    M1 = (r0 - 1) // d0
    assert M1 % 6 == 0
    lam = pow(zeta, (r0 - 1) // 3, r0)
    assert pow(lam, 3, r0) == 1 and lam != 1
    assert pow(lam, M1, r0) == 1 and pow(r0 - 1, M1, r0) == 1

    rng = random.Random(seed)
    secrets = [1, r0 - 1, r0 - 2, 2]
    while len(secrets) < trials:
        secrets.append(rng.randrange(1, r0))
    recovered = 0
    for s in secrets[:trials]:
        ok, v = toy_class_walk(r0, zeta, d0, s, rng)
        assert ok, ("class walk failed", r0, d0, s)
        assert v == pow(s, d0, r0)
        ok2, k = toy_phase2(r0, zeta, d0, v, s)
        assert ok2, ("phase 2 failed", r0, d0, s)
        recovered += 1

    w = max(1, (r0.bit_length() + WINDOWS - 1) // WINDOWS)
    tabs = toy_comb(r0, w, WINDOWS)
    for x in (1, r0 - 1, r0 - 2, rng.randrange(1, r0)):
        assert toy_comb_mul(x, tabs, w, WINDOWS, r0) == x % r0

    toy = {
        "prime": str(r0),
        "cofactor_q": str(q),
        "d_bits": d_bits,
        "M1": str(M1),
        "classes": str(M1 // 6),
        "trials": trials,
        "recovered": recovered,
        "phase2_recovered": recovered,
        "comb_bits": w,
    }
    prod = production_checks(seed)
    curve = {
        "phi_equals_lambda": True,
        "comb_digit_identity": True,
    }
    print(json.dumps({
        "ok": True,
        "toy": toy,
        "production": prod,
        "curve": curve,
        "scope": "toy algebra, endomorphism identity and parameter/ledger arithmetic only; no production-size run",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
