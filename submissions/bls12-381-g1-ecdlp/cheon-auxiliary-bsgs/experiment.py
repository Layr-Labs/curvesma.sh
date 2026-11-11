"""Toy Cheon check; not evidence of full-size costs or security.

Implements the exact two-phase BSGS logic of proof.md on the additive
cyclic group Z/qZ (prime q, d | q-1), with direct scalar multiplications
in place of the production b-adic tables, and separately checks the
b-adic table evaluation identity on the same toy group. This checks
algorithm logic only; it does not validate production cost bounds,
memory bounds, or G1 arithmetic.
"""
import json
import math
import sys


def prime_factors(n):
    out = set()
    p = 2
    while p * p <= n:
        if n % p == 0:
            out.add(p)
            while n % p == 0:
                n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out.add(n)
    return out


def cheon(q, d, s):
    assert 2 <= q and (q - 1) % d == 0 and 1 <= s < q
    m = (q - 1) // d
    facs = prime_factors(q - 1)
    zeta = next(a for a in range(2, q)
                if all(pow(a, (q - 1) // ell, q) != 1 for ell in facs))
    eta = pow(zeta, d, q)
    theta = pow(zeta, m, q)
    a1 = s % q
    ad = pow(s, d, q)
    # Phase 1: find k0 in [0, m) with s^d = eta^k0.
    d1 = math.isqrt(m - 1) + 1
    baby = {pow(eta, j, q): j for j in range(d1)}
    k0 = None
    for i in range(d1):
        cand = ad * pow(eta, -i * d1, q) % q
        if cand in baby:
            k0 = i * d1 + baby[cand]
            break
    assert k0 is not None and pow(eta, k0, q) == ad
    # Phase 2: s = zeta^k0 * theta^k2; find k2 in [0, d).
    base = pow(zeta, k0, q)
    d2 = math.isqrt(d - 1) + 1
    baby2 = {pow(theta, j, q) * base % q: j for j in range(d2)}
    k2 = None
    for i in range(d2):
        cand = a1 * pow(theta, -i * d2, q) % q
        if cand in baby2:
            k2 = i * d2 + baby2[cand]
            break
    assert k2 is not None
    alpha = pow(zeta, k0, q) * pow(theta, k2, q) % q
    assert alpha == s
    return {"q": q, "d": d, "s": s, "zeta": zeta, "k0": k0, "k2": k2,
            "recovered": alpha}


def badic_table_check(q):
    # Table T(i, j) = j * b^i (mod q); sum of digit entries must rebuild delta.
    b = math.isqrt(q - 1) + 1
    c = 2
    assert b ** c >= q
    table = [[(j * (b ** i)) % q for j in range(b)] for i in range(c)]
    for delta in range(min(200, q)):
        digits = [(delta // (b ** i)) % b for i in range(c)]
        assert sum(table[i][digits[i]] for i in range(c)) % q == delta % q
    return {"b": b, "c": c, "checked": min(200, q)}


inp = json.load(sys.stdin)
q, d, s = inp["prime"], inp["d"], inp["secret"]
result = cheon(q, d, s)
table = badic_table_check(q)
print(json.dumps({"group": "integers modulo prime q under addition",
                  "cheon": result, "badic_table": table,
                  "scope": "toy group only; does not validate production cost bounds"}))
