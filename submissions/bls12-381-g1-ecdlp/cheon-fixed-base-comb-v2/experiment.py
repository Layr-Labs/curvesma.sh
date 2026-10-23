
"""Recompute the Cheon/comb bound and check the index arithmetic on toy groups."""
import json
import math
import sys
from decimal import Decimal, localcontext

BLS_X = 0xD201000000010000
BLS_R = BLS_X**4 - BLS_X**2 + 1
SECP_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP_A = 0
SECP_GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
SECP_GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
EXPECTED_W = 9388543893536700382679492617175715396
EXPECTED_NUM = 122820313
EXPECTED_DEN = 1000000
EXPECTED_M = 223621884914696853578958571250384896
EXPECTED_V = 223621884914696853578958571250384896
EXPECTED_OPS = 894488015027762499901859846266284394
MULTIPLICITY = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1, 906349: 2, 2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}


def ceil_sqrt(n):
    s = math.isqrt(n)
    return s if s * s == n else s + 1


def ceil_div(a, b):
    return (a + b - 1) // b


def is_prime(n):
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    for p in small:
        if n % p == 0:
            return n == p
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 23):
        if a >= n:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def pollard(n):
    def f(x, c):
        return (x * x + c) % n

    for c in range(1, 40):
        x = y = 2
        d = 1
        while d == 1:
            x = f(x, c)
            y = f(f(y, c), c)
            d = math.gcd(abs(x - y), n)
        if 1 < d < n:
            return d
    raise SystemExit("factorization failed")


def factor_number(n, into):
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
    for p in small:
        while n % p == 0:
            into[p] = into.get(p, 0) + 1
            n //= p
    if n == 1:
        return
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            into[m] = into.get(m, 0) + 1
            continue
        d = pollard(m)
        stack.append(d)
        stack.append(m // d)


def factor_r_minus_1():
    n = BLS_X * BLS_X * (BLS_X - 1) * (BLS_X + 1)
    if n != BLS_R - 1:
        raise SystemExit("closed form of r-1 failed")
    if (BLS_R - 1) % (1 << 20) != 0:
        raise SystemExit("2^20 does not divide r-1")
    into = {}
    factor_number(BLS_X, into)
    factor_number(BLS_X, into)
    factor_number(BLS_X - 1, into)
    factor_number(BLS_X + 1, into)
    rebuilt = 1
    for p, e in into.items():
        rebuilt *= p**e
    if rebuilt != BLS_R - 1:
        raise SystemExit("factors do not rebuild r-1")
    if into != MULTIPLICITY:
        raise SystemExit(f"unexpected factorization {into}")
    return sorted(into)


def primitive_root(modulus, primes):
    checks = [(modulus - 1) // q for q in primes]
    for candidate in range(2, 10000):
        if all(pow(candidate, exp, modulus) != 1 for exp in checks):
            return candidate
    raise SystemExit("no small primitive root")


def cheon_toy(r, d, secret):
    n = r - 1
    primes = []
    p = 2
    while p * p <= n:
        while n % p == 0:
            primes.append(p)
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        primes.append(n)
    zeta = primitive_root(r, sorted(set(primes)))
    beta = pow(zeta, d, r)
    n_big = (r - 1) // d
    m = ceil_sqrt(n_big)
    v_max = ceil_div(n_big, m)
    pd = pow(secret, d, r)
    beta_inv = pow(beta, -1, r)
    babies = {}
    scalar = 1
    for u in range(m):
        babies[scalar * pd % r] = u
        scalar = scalar * beta_inv % r
    giant = 1
    beta_m = pow(beta, m, r)
    k0 = None
    for v in range(v_max):
        if giant in babies:
            k0 = babies[giant] + m * v
            break
        giant = giant * beta_m % r
    if k0 is None:
        raise SystemExit("phase 1 missed")
    gamma = pow(zeta, n_big, r)
    m2 = ceil_sqrt(d)
    v2 = ceil_div(d, m2)
    left = secret * pow(zeta, -k0, r) % r
    gamma_inv = pow(gamma, -1, r)
    babies2 = {}
    scalar = 1
    for u in range(m2):
        babies2[scalar * left % r] = u
        scalar = scalar * gamma_inv % r
    giant = 1
    gamma_m = pow(gamma, m2, r)
    k1 = None
    for v in range(v2):
        if giant in babies2:
            k1 = babies2[giant] + m2 * v
            break
        giant = giant * gamma_m % r
    if k1 is None:
        raise SystemExit("phase 2 missed")
    recovered = pow(zeta, k0 + k1 * n_big, r)
    if recovered != secret % r:
        raise SystemExit("recovered scalar does not match")


def recompute():
    if (1 << 64) <= BLS_X + 1 or (1 << 32) ** 2 < BLS_X + 1:
        raise SystemExit("2^32 trials do not cover sqrt of the 64-bit inputs")
    d = 1 << 20
    n_big = (BLS_R - 1) // d
    m = ceil_sqrt(n_big)
    v_max = ceil_div(n_big, m)
    window = 96
    t = ceil_div(255, window)
    bases = 2
    pre_doubles = bases * window * (t - 1)
    pre_adds = bases * t * ((1 << window) - 1)
    online = (m + v_max) * (t - 1)
    phase2 = (2 * (1 << 10) + 4) * (254 + 254)
    group_ops = pre_doubles + pre_adds + online + phase2
    g_ins = 16 * 512 + 24 * 128
    # 20 squarings for zeta^(2^20), six 510-step binary exponentiations, 16 assembly multiplications.
    field_setup = 20 + 6 * 510 + 16
    if field_setup != 3096:
        raise SystemExit("field-setup count drifted")
    scalar_ins = (m + v_max + field_setup) * 2048
    points = m + v_max
    fermat = 2 * 381
    affine_ins = (7 * points + fermat) * 512 + 4 * points * 128
    sort_ins = points * 24 * 512 + 24 * 65536 * 32
    scan_ins = points * 1024
    root_ins = 64 * (2 * 255) * 2048
    factor_ins = 3 * (1 << 32) * 32 + 16 * 2048 + 64 * 2048 + (1 << 30)
    input_load_ins = 2 * 144 * 2
    reserve = 1 << 40
    total = (
        group_ops * g_ins
        + scalar_ins
        + affine_ins
        + sort_ins
        + scan_ins
        + root_ins
        + factor_ins
        + input_load_ins
        + reserve
    )
    w = ceil_div(total, 2000)
    comb_points = bases * t * (1 << window)
    affine_table_ins = (3 * comb_points + 2 * 381) * 512 + 4 * comb_points * 128
    affine_w = ceil_div(total + affine_table_ins, 2000)
    return {
        "m": m,
        "v_max": v_max,
        "group_ops": group_ops,
        "W": w,
        "digits": t,
        "factor_ins": factor_ins,
        "affine_W": affine_w,
        "field_setup": field_setup,
    }


def log_ok(value, num, den):
    with localcontext() as ctx:
        ctx.prec = 80
        log2 = Decimal(value).ln() / Decimal(2).ln()
        claimed = Decimal(num) / Decimal(den)
    if claimed < log2:
        raise SystemExit("declared score is below log2(W)")
    if claimed - log2 > Decimal("0.000002"):
        raise SystemExit("declared score is more than two ticks above log2(W)")


def affine_add(p, q):
    if p is None:
        return q
    if q is None:
        return p
    if p[0] == q[0]:
        if (p[1] + q[1]) % SECP_P == 0:
            return None
        lam = (3 * p[0] * p[0]) * pow(2 * p[1], -1, SECP_P) % SECP_P
    else:
        lam = (q[1] - p[1]) * pow(q[0] - p[0], -1, SECP_P) % SECP_P
    x = (lam * lam - p[0] - q[0]) % SECP_P
    y = (lam * (p[0] - x) - p[1]) % SECP_P
    return x, y


def binary_mul(point, scalar):
    acc = None
    base = point
    while scalar:
        if scalar & 1:
            acc = affine_add(acc, base)
        base = affine_add(base, base)
        scalar >>= 1
    return acc


def comb_mul(point, scalar, width, bits):
    t = ceil_div(bits, width)
    step = point
    tables = []
    for i in range(t):
        if i:
            for _ in range(width):
                step = affine_add(step, step)
        table = [None]
        acc = None
        for _ in range(1, 1 << width):
            acc = affine_add(acc, step)
            table.append(acc)
        tables.append(table)
    acc = None
    for i in reversed(range(t)):
        digit = (scalar >> (width * i)) & ((1 << width) - 1)
        if digit:
            acc = affine_add(acc, tables[i][digit])
    return acc


def main():
    payload = json.load(sys.stdin)
    toy_prime = int(payload["toy_prime"])
    primes = factor_r_minus_1()
    if primitive_root(BLS_R, primes) != 7:
        raise SystemExit("7 is not the primitive root")
    if pow(7, BLS_R - 1, BLS_R) != 1:
        raise SystemExit("Fermat check failed")
    got = recompute()
    if got["W"] != EXPECTED_W or got["m"] != EXPECTED_M or got["v_max"] != EXPECTED_V:
        raise SystemExit(f"recomputed bound {got} does not match the claim")
    if got["affine_W"] <= got["W"]:
        raise SystemExit("normalized affine tables cost no more than Jacobian tables")
    if got["factor_ins"] <= 0 or got["field_setup"] != 3096:
        raise SystemExit("factorization or field setup missing from the recomputation")
    if got["group_ops"] != EXPECTED_OPS or got["digits"] != 3:
        raise SystemExit("group-operation count mismatch")
    log_ok(got["W"], EXPECTED_NUM, EXPECTED_DEN)
    toys = []
    for r, secret in ((101, 1), (101, 100), (211, 17), (10007, 4242), (65537, 99)):
        d = 2
        while (r - 1) % (d * 2) == 0 and d < 64:
            d *= 2
        cheon_toy(r, d, secret)
        toys.append({"order": r, "d": d, "secret": secret})
    if is_prime(toy_prime) and toy_prime <= 65537 and (toy_prime - 1) % 2 == 0:
        d = 2
        while (toy_prime - 1) % (d * 2) == 0 and d < 64:
            d *= 2
        secret = 1 + toy_prime // 3
        cheon_toy(toy_prime, d, secret % toy_prime)
    scalar = 0x123456789ABCDEF
    base = (SECP_GX, SECP_GY)
    if comb_mul(base, scalar, 8, 64) != binary_mul(base, scalar):
        raise SystemExit("comb does not match binary multiplication")
    print(json.dumps({
        "track": "bls12-381-g1-ecdlp",
        "W": str(got["W"]),
        "m": str(got["m"]),
        "score": [EXPECTED_NUM, EXPECTED_DEN],
        "primitive_root": 7,
        "toys": toys,
        "factor_ins": str(got["factor_ins"]),
        "affine_W": str(got["affine_W"]),
        "tables": "jacobian-144",
        "advice_bytes": 0,
        "comb_matches_binary": True,
        "scope": "index arithmetic and a 64-bit comb check, not a production discrete log",
    }))


if __name__ == "__main__":
    main()
