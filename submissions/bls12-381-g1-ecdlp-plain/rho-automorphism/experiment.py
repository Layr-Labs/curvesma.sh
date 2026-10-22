"""Toy class-rho and GLV identities; not evidence of full-size costs."""
import json
import sys

ORDER = json.load(sys.stdin)["order"]

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
Z = -0xD201000000010000
LAMBDA = 0xAC45A4010001A40200000000FFFFFFFF
BETA = 0x1A0111EA397FE699EC02408663D4DE85AA0D857D89759AD4897D29650FB85F9B409427EB4F49FFFD8BFD00000000AAAC


def check_profile_identities():
    assert (R - 1) % 6 == 0
    assert (Z * Z - 1) % R == LAMBDA
    assert (LAMBDA * LAMBDA + LAMBDA + 1) % R == 0
    assert (BETA * BETA + BETA + 1) % P == 0
    return {"r_mod_6": 0, "lambda_cube_root": True, "beta_cube_root": True}


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def generator_of_automorphisms(n):
    """Multiplicative generator of a 2- or 6-element acting group on Z/nZ."""
    if is_prime(n) and (n - 1) % 6 == 0:
        for g in range(2, n):
            if pow(g, (n - 1) // 2, n) == n - 1 and pow(g, (n - 1) // 3, n) != 1:
                omega = pow(g, (n - 1) // 6, n)
                if pow(omega, 6, n) == 1 and pow(omega, 3, n) == n - 1:
                    return omega, 6
        raise AssertionError("no order-6 element")
    return n - 1, 2


def orbit(y, omega, n, size):
    out = []
    v = y % n
    for _ in range(size):
        out.append(v)
        v = (omega * v) % n
    return out


def canon(y, omega, n, size):
    if y % n == 0:
        return 0
    return min(orbit(y, omega, n, size))


def class_mu(y, y_canon, omega, n, size):
    v = y % n
    mu = 1
    for _ in range(size):
        if v == y_canon:
            return mu
        v = (omega * v) % n
        mu = (mu * omega) % n
    raise AssertionError("not in orbit")


def apply_canon_coeffs(a, b, secret, omega, n, size):
    raw = (a + b * secret) % n
    w = canon(raw, omega, n, size)
    mu = 1 if raw == 0 else class_mu(raw, w, omega, n, size)
    return (mu * a) % n, (mu * b) % n, w


def recover_from_relation(a1, b1, a2, b2, mu, n):
    left = (b1 - mu * b2) % n
    if left == 0:
        return None
    return ((mu * a2 - a1) * pow(left, -1, n)) % n


def check_solver(n, omega, size, trials=40):
    found = 0
    seed = 0
    while found < trials:
        seed += 1
        x = 1 + (seed * 5) % (n - 1)
        a1 = (3 * seed + 1) % n
        b1 = (7 * seed + 2) % n or 1
        mu = pow(omega, seed % size, n)
        b2 = (13 * seed + 4) % n or 1
        if (b1 - mu * b2) % n == 0:
            continue
        w = (a1 + b1 * x) % n
        a2_forced = (pow(mu, -1, n) * w - b2 * x) % n
        got = recover_from_relation(a1, b1, a2_forced, b2, mu, n)
        assert got == x
        found += 1
    return trials


def toy_rho_all(n, omega, size, r_parts=16, max_steps=500, retries=32):
    parts = [((j * 5 + 1) % n, (j * 9 + 3) % n) for j in range(r_parts)]

    def h(w):
        return w % r_parts

    recovered = 0
    for secret in range(n):
        ok = False
        for attempt in range(retries):
            table = {}
            starts = [
                ((1 + attempt) % n, (2 + 3 * attempt) % n),
                ((5 + 7 * attempt) % n, (11 + attempt) % n),
                ((13 + 17 * attempt) % n, (19 + 23 * attempt) % n),
            ]
            states = [list(apply_canon_coeffs(a, b, secret, omega, n, size)) + [0]
                      for a, b in starts]
            for _ in range(max_steps):
                for s in states:
                    a, b, w, ln = s
                    c, d = parts[h(w)]
                    a, b, w = apply_canon_coeffs(
                        (a + c) % n, (b + d) % n, secret, omega, n, size)
                    ln += 1
                    s[0], s[1], s[2], s[3] = a, b, w, ln
                    assert w == (a + b * secret) % n or w == 0
                    if w in table:
                        a1, b1 = table[w]
                        den = (b1 - b) % n
                        if den != 0:
                            x = ((a - a1) * pow(den, -1, n)) % n
                            if x == secret:
                                ok = True
                                break
                    else:
                        table[w] = (a, b)
                if ok:
                    break
            if ok:
                break
        if not ok:
            raise AssertionError("Discrete logarithm not recovered for secret %s" % secret)
        recovered += 1
    return recovered


def finite_bsgs_all(n):
    m = int(n ** 0.5) + 1
    babies = {j % n: j for j in range(m)}
    for secret in range(n):
        value = secret
        found = False
        for i in range((n + m - 1) // m):
            if value in babies:
                assert (i * m + babies[value]) % n == secret
                found = True
                break
            value = (value - m) % n
        if not found:
            raise AssertionError("BSGS failed for secret %s" % secret)
    return n


profile = check_profile_identities()
if is_prime(ORDER):
    omega, size = generator_of_automorphisms(ORDER)
    seen = set()
    orbits = 0
    for y in range(1, ORDER):
        if y in seen:
            continue
        o = orbit(y, omega, ORDER, size)
        assert len(set(o)) == size
        seen.update(o)
        orbits += 1
    assert orbits == (ORDER - 1) // size
    solver_trials = check_solver(ORDER, omega, size)
    recovered = toy_rho_all(ORDER, omega, size)
    method = "class-rho"
else:
    omega, size, orbits, solver_trials = None, 0, 0, 0
    recovered = finite_bsgs_all(ORDER)
    method = "bsgs-composite-fallback"

print(json.dumps({
    "group": "integers modulo order under addition",
    "order": ORDER,
    "method": method,
    "automorphism_order": size,
    "orbits": orbits,
    "profile_identities": profile,
    "solver_trials": solver_trials,
    "checked": recovered,
    "all_recovered": True,
    "scope": "toy group and profile identities only; does not validate production cost bounds",
}))
