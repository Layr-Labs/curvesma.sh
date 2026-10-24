"""Toy of the order-6 automorphism walk in proof.md.

The group is additive Z/n with n prime and n ≡ 1 (mod 6), so
X^2+X+1 splits and the six multipliers {±1, ±λ, ±λ²} are distinct.
Coefficients stay in the Z[λ] basis and are updated by the same
addition-only maps as the production claim. The observed point is the
expanded scalar; the secret is used only to build that point and to
check the recovered relation. This does not price secp256k1 arithmetic.
"""

from __future__ import annotations

import json
import math
import sys
from fractions import Fraction

WALKS = 8
BRANCHES = 16
PAD_NUM, PAD_DEN = 5, 4


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def inv(a: int, n: int) -> int:
    return pow(a % n, -1, n)


def ln2_upper() -> Fraction:
    """ln 2 < sum_{k=1..60} 1/(k 2^k) + 1/(61 2^60)."""
    s = Fraction(0)
    for k in range(1, 61):
        s += Fraction(1, k << k)
    return s + Fraction(1, 61 << 60)


def half_budget(n_orb: int) -> int:
    need = (2 * n_orb * ln2_upper()).__ceil__()
    t = math.isqrt(need)
    while t * (t + 1) < need:
        t += 1
    return t


def padded_budget(n_orb: int) -> int:
    t0 = half_budget(n_orb)
    return (t0 * PAD_NUM + PAD_DEN - 1) // PAD_DEN


def glv_lambda(n: int) -> int:
    disc = (-3) % n
    s, q = 0, n - 1
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (n - 1) // 2, n) != n - 1:
        z += 1
    m, c, t, r = s, pow(z, q, n), pow(disc, q, n), pow(disc, (q + 1) // 2, n)
    while t != 1:
        i, tt = 1, (t * t) % n
        while tt != 1:
            tt = (tt * tt) % n
            i += 1
        b = pow(c, 1 << (m - i - 1), n)
        m, c, t, r = i, (b * b) % n, (t * b * b) % n, (r * b) % n
    lam = ((-1 + r) * inv(2, n)) % n
    if (lam * lam + lam + 1) % n != 0:
        raise RuntimeError("lambda check failed")
    return lam


def multipliers(lam: int, n: int) -> list[int]:
    lam2 = (lam * lam) % n
    return [1, (n - 1) % n, lam, (n - lam) % n, lam2, (n - lam2) % n]


def mul_unit(u: int, v: int, kind: int, n: int) -> tuple[int, int]:
    """Multiply u+vλ by the kind-th unit. kind: 1,-1,λ,-λ,λ²,-λ²."""
    if kind == 0:
        return u % n, v % n
    if kind == 1:
        return (-u) % n, (-v) % n
    if kind == 2:
        return (-v) % n, (u - v) % n
    if kind == 3:
        return v % n, (v - u) % n
    if kind == 4:
        return (v - u) % n, (-u) % n
    if kind == 5:
        return (u - v) % n, u % n
    raise ValueError(kind)


def expand(u: int, v: int, lam: int, n: int) -> int:
    return (u + v * lam) % n


def branch_index(cv: int) -> int:
    mixed = ((cv & ((1 << 64) - 1)) * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
    return mixed >> 60


def branches(n: int) -> list[tuple[int, int]]:
    return [((2 * j + 1) % n, (3 * j + 5) % n) for j in range(BRANCHES)]


def canon(value: int, mults: list[int], n: int) -> tuple[int, int]:
    best_v, best_k = n, 0
    for kind, m in enumerate(mults):
        img = (m * value) % n
        if img < best_v:
            best_v, best_k = img, kind
    return best_v, best_k


def solve(a: tuple[int, int, int, int], b: tuple[int, int, int, int], lam: int, n: int, secret: int) -> bool:
    ae1 = expand(a[0], a[1], lam, n)
    be1 = expand(a[2], a[3], lam, n)
    ae2 = expand(b[0], b[1], lam, n)
    be2 = expand(b[2], b[3], lam, n)
    den = (be1 - be2) % n
    if den == 0:
        return False
    num = (ae2 - ae1) % n
    return (num * inv(den, n)) % n == secret


def recover(secret: int, n: int, lam: int, mults: list[int], br: list[tuple[int, int]], budget: int) -> tuple[bool, int]:
    seen: dict[int, tuple[int, int, int, int]] = {}
    states: list[list[int]] = []
    fruitless = 0

    def adopt(au: int, av: int, bu: int, bv: int, value: int) -> tuple[int, int, int, int, int] | None:
        nonlocal fruitless
        cv, kind = canon(value, mults, n)
        au, av = mul_unit(au, av, kind, n)
        bu, bv = mul_unit(bu, bv, kind, n)
        got = (expand(au, av, lam, n) + secret * expand(bu, bv, lam, n)) % n
        if got != cv:
            raise AssertionError("basis map drifted from the expanded scalar")
        coeffs = (au, av, bu, bv)
        if cv == 0:
            # α + β k ≡ 0
            den = expand(bu, bv, lam, n)
            if den != 0 and ((-expand(au, av, lam, n)) * inv(den, n)) % n == secret:
                return None
        prev = seen.get(cv)
        if prev is not None:
            if solve(prev, coeffs, lam, n, secret):
                return None
            fruitless += 1
        else:
            seen[cv] = coeffs
        return [au, av, bu, bv, cv]

    for w in range(WALKS):
        st = adopt((w + 1) % n, 0, 1, 0, (w + 1 + secret) % n)
        if st is None:
            return True, fruitless
        states.append(st)

    rounds = (budget + WALKS - 1) // WALKS
    for _ in range(rounds):
        for st in states:
            au, av, bu, bv, cv = st
            a, b = br[branch_index(cv)]
            au = (au + a) % n
            bu = (bu + b) % n
            value = (expand(au, av, lam, n) + secret * expand(bu, bv, lam, n)) % n
            nxt = adopt(au, av, bu, bv, value)
            if nxt is None:
                return True, fruitless
            st[:] = nxt
    return False, fruitless


def run(n: int) -> dict:
    if n < 67 or n > 2003 or n % 6 != 1 or not is_prime(n):
        raise ValueError("order must be a prime ≡ 1 (mod 6) in [67, 2003]")
    lam = glv_lambda(n)
    mults = multipliers(lam, n)
    if len(set(mults)) != 6:
        raise AssertionError("multipliers not distinct")
    # Every nonzero orbit has size 6 iff the multipliers are distinct mod n.
    orbit_bad = 0
    for v in range(1, min(n, 64)):
        if len({(m * v) % n for m in mults}) != 6:
            orbit_bad += 1
    n_orb = 1 + (n - 1) // 6
    t0 = half_budget(n_orb)
    budget = padded_budget(n_orb)
    br = branches(n)
    hits = 0
    fruitless = 0
    for secret in range(n):
        ok, fr = recover(secret, n, lam, mults, br, budget)
        hits += int(ok)
        fruitless += fr
    return {
        "group": "integers modulo order under addition, order-6 multiplier action",
        "order": n,
        "orbits": n_orb,
        "half_budget_steps": t0,
        "budget_steps": budget,
        "walks": WALKS,
        "branches": BRANCHES,
        "checked": n,
        "recovered": hits,
        "fruitless_relations": fruitless,
        "success": hits / n,
        "orbit_sample_failures": orbit_bad,
        "basis": "Z[lambda]/(lambda^2+lambda+1), addition-only unit maps",
        "scope": "toy group only; collision rate and coefficient identity, not production instruction costs",
    }


if __name__ == "__main__":
    print(json.dumps(run(json.load(sys.stdin)["order"])))
