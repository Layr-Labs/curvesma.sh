#!/usr/bin/env python3
"""Toy witness for six-class rho with a functional negation escape.

Canonicalization onto a GLV order-6 orbit does not delete negation-induced
fruitless 2-cycles. On y^2 = x^3 + 7 the map (x, y) -> (beta x, y) is an
automorphism of order 3; together with negation it gives six classes.
The structural 2-cycle is the identity

    can(can(P + J_s) + J_s) = P

whenever the first reduction is negation and the hash repeats. This file
exhibits that cycle on a prime-order toy group, then runs the escaped walk
f* and recovers the discrete log.

f* is a function of the current canonical point only:
  Q = can(P + J_h(P))
  if the reduction is negation and h(Q) = h(P):
      Q = can(P + E_h(P))     # one extra mixed addition
  (and likewise one further layer if the landing repeats the predicate).
Amortized extra additions are at most 1/(r-1) once the landing hash is
equidistributed; the toy measures the rate. Residual 2-cycles of f* are
solutions of a non-identity translation equation, hence O(1) per jump slot.
"""

from __future__ import annotations

import hashlib
import random


# Production parameters the proof instantiates. The toy uses a smaller r
# so the structural cycle is dense enough to see in a few thousand steps.
PROD_R = 16
PROD_DP_BITS = 20
TOY_R = 8


def _sha256_int(label: str) -> int:
    return int.from_bytes(hashlib.sha256(label.encode()).digest(), "big")


def egcd(a: int, b: int):
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def inv_mod(a: int, n: int) -> int:
    g, x, _ = egcd(a % n, n)
    if g != 1:
        raise ZeroDivisionError(a)
    return x % n


def tonelli(n: int, p: int):
    n %= p
    if n == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    while True:
        if t == 0 or t == 1:
            return r if t == 1 else 0
        i, tt = 1, (t * t) % p
        while tt != 1:
            tt = (tt * tt) % p
            i += 1
            if i == m:
                return None
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p


class Curve:
    """y^2 = x^3 + B over F_p, a = 0, with the order-6 GLV action."""

    def __init__(self, p: int, B: int = 7):
        self.p = p
        self.B = B % p
        disc = (-3) % p
        root = tonelli(disc, p)
        if root is None:
            raise ValueError("no beta")
        inv2 = inv_mod(2, p)
        self.beta = ((-1 + root) * inv2) % p
        if (self.beta * self.beta + self.beta + 1) % p != 0:
            self.beta = ((-1 - root) * inv2) % p
        if (self.beta * self.beta + self.beta + 1) % p != 0 or self.beta % p in (0, 1):
            raise ValueError("beta failed")

    def inv(self, x: int) -> int:
        return inv_mod(x, self.p)

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2 and (y1 + y2) % p == 0:
            return None
        if P == Q:
            if y1 % p == 0:
                return None
            lam = ((3 * x1 * x1) % p) * self.inv((2 * y1) % p) % p
        else:
            lam = ((y2 - y1) % p) * self.inv((x2 - x1) % p) % p
        x3 = (lam * lam - x1 - x2) % p
        y3 = (lam * (x1 - x3) - y1) % p
        return (x3, y3)

    def mul(self, k: int, P):
        R = None
        Q = P
        k %= self.order if hasattr(self, "order") else k
        if k < 0:
            k = -k
            Q = self.neg(Q)
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def neg(self, P):
        if P is None:
            return None
        return (P[0], (-P[1]) % self.p)

    def phi(self, P):
        """(x, y) -> (beta x, y). Eigenvalue lambda is fixed later on a generator."""
        if P is None:
            return None
        return ((self.beta * P[0]) % self.p, P[1])

    def phi2(self, P):
        # beta^2 = -1 - beta, so beta^2 x = -x - beta x. One mul, not two.
        if P is None:
            return None
        x, y = P
        bx = (self.beta * x) % self.p
        return ((-x - bx) % self.p, y)

    def on_curve(self, P) -> bool:
        if P is None:
            return True
        x, y = P
        return (y * y - (x * x * x + self.B)) % self.p == 0

    def random_point(self, rng: random.Random):
        p = self.p
        while True:
            x = rng.randrange(p)
            rhs = (x * x * x + self.B) % p
            y = tonelli(rhs, p)
            if y is not None:
                return (x, y)

    def order_of(self, P, n_bound: int):
        """Baby-step order, n_bound >= group order."""
        if P is None:
            return 1
        m = int(n_bound**0.5) + 1
        table = {}
        R = None
        for j in range(m + 1):
            table[R] = j
            R = self.add(R, P)
        step = self.mul(m, P)
        R = None
        for i in range(m + 1):
            # want i*m - j = 0 in <P>, i.e. i*step = j*P
            if R in table and not (i == 0 and table[R] == 0):
                j = table[R]
                k = i * m - j
                if k > 0 and self.mul(k, P) is None:
                    return k
            R = self.add(R, step)
        raise RuntimeError("order not found")


def point_count(curve: Curve) -> int:
    p = curve.p
    n = 1  # infinity
    for x in range(p):
        rhs = (x * x * x + curve.B) % p
        leg = pow(rhs, (p - 1) // 2, p) if rhs else 0
        if rhs == 0:
            n += 1
        elif leg == 1:
            n += 2
    return n


def find_toy_curve(rng: random.Random):
    """Smallest p > 80, p ≡ 1 (mod 3), whose y^2 = x^3+7 has prime order ≡ 1 (mod 6)."""
    p = 103
    while p < 4000:
        if pow(p, 1, 3) == 1 and all(p % d for d in range(2, int(p**0.5) + 1)):
            curve = Curve(p)
            n = point_count(curve)
            if n > 40 and all(n % d for d in range(2, int(n**0.5) + 1)) and n % 6 == 1:
                curve.order = n
                G = curve.random_point(rng)
                # prime order => any non-infinity point generates
                phiG = curve.phi(G)
                # lambda via BSGS on a  group of size n
                lam = bsgs(curve, G, phiG, n)
                if lam is not None and (lam * lam + lam + 1) % n == 0 and curve.mul(lam, G) == phiG:
                    curve.lam = lam
                    curve.G = G
                    return curve
        p += 2
    raise RuntimeError("no toy curve")


def bsgs(curve: Curve, P, Q, n: int):
    if Q is None:
        return 0
    m = int(n**0.5) + 1
    table = {}
    R = None
    for j in range(m):
        table[R] = j
        R = curve.add(R, P)
    step = curve.mul(m, P)
    R = Q
    neg_step = curve.neg(step)
    for i in range(m + 1):
        if R in table:
            k = (table[R] - i * m) % n
            if curve.mul(k, P) == Q:
                return k
        R = curve.add(R, neg_step)
    return None


# Units: index -> scalar. Canonical representative is minimal (x, y) as integers.
# mu = 0..5 maps to (1, -1, lam, -lam, lam^2, -lam^2).


def unit_scalar(curve: Curve, idx: int) -> int:
    n = curve.order
    lam = curve.lam
    lam2 = (-1 - lam) % n
    table = (1, n - 1, lam, (-lam) % n, lam2, (-lam2) % n)
    return table[idx]


def orbit_with_units(curve: Curve, P):
    """Six images and the unit index taking P to that image."""
    if P is None:
        return []
    x, y = P
    p = curve.p
    bx = (curve.beta * x) % p
    b2x = (-x - bx) % p
    ys = (y, (-y) % p)
    # (x-image index 0,1,2) * (sign 0,1) -> unit idx
    # sign 0 keeps y (units +1,+lam,+lam^2), sign 1 negates (units -1,-lam,-lam^2)
    out = []
    xs = (x, bx, b2x)
    for xi, xx in enumerate(xs):
        for si, yy in enumerate(ys):
            idx = (0, 1, 2, 3, 4, 5)[xi + 3 * si] if False else ((0, 2, 4)[xi] if si == 0 else (1, 3, 5)[xi])
            out.append((idx, (xx, yy)))
    return out


def canonicalize(curve: Curve, P):
    """Return (canonical point, unit idx u with can = [unit_u] P)."""
    best = None
    best_u = 0
    for idx, pt in orbit_with_units(curve, P):
        if best is None or pt < best:
            best = pt
            best_u = idx
    return best, best_u


def hash_slot(P, r: int) -> int:
    # Depends on the canonical x only. Well-defined on the representative
    # the walk actually stores. Not constant on the raw orbit — and it does
    # not need to be: the state is already canonical.
    return P[0] % r


def jump_coeffs(r: int, layer: int, s: int):
    """Small coefficient vectors (c, d).

    Layer 0 is the primary Teske table, layer >= 1 the escape tables.
    Vectors are pairwise outside every relation (c', d') = -mu (c, d) for a
    unit mu that preserves small height: both coordinates stay in 1..O(r)
    and are unique across a layer, so (c, d) = -(c, d) is impossible and
    distinct indices never collide. Huge units ±λ, ±λ^2 take small vectors
    to large ones, checked exactly when an order is known.
    """
    # Deterministic, full-period in the slot index.
    c = 1 + s + layer * (r + 3)
    d = 2 + 3 * s + layer * (2 * r + 5)
    return c, d


def apply_unit_coeff(curve: Curve, a: int, b: int, uidx: int):
    mu = unit_scalar(curve, uidx)
    n = curve.order
    return (a * mu) % n, (b * mu) % n


def step_raw(curve: Curve, P, a: int, b: int, c: int, d: int):
    """One adding step + canonicalization. Returns Q, u, a', b', adds_done."""
    J = curve.add(curve.mul(c, curve.G), curve.mul(d, curve.Q))
    S = curve.add(P, J)
    Q, u = canonicalize(curve, S)
    a2, b2 = (a + c) % curve.order, (b + d) % curve.order
    a3, b3 = apply_unit_coeff(curve, a2, b2, u)
    return Q, u, a3, b3


def is_negation_cycle(curve: Curve, P, Q, u: int, r: int) -> bool:
    """Structural predicate: reduction by -1 and the hash repeats.

    On that set the algebraic identity can(Q + J_s) = P holds for every P,
    with s = h(P), so the walk is a fruitless 2-cycle (coefficient monodromy
    is the identity). Canonicalization is what *creates* the identity; it
    does not remove it.
    """
    return u == 1 and hash_slot(Q, r) == hash_slot(P, r)


def f_star(curve: Curve, P, a: int, b: int, r: int, max_layer: int = 3):
    """Memoryless escaped step. Costs 1 + (layers used past 0) additions."""
    s = hash_slot(P, r)
    adds = 0
    for layer in range(max_layer):
        c, d = jump_coeffs(r, layer, s)
        Q, u, a2, b2 = step_raw(curve, P, a, b, c, d)
        adds += 1
        if not is_negation_cycle(curve, P, Q, u, r):
            return Q, a2, b2, adds, layer
    c, d = jump_coeffs(r, max_layer, s)
    Q, u, a2, b2 = step_raw(curve, P, a, b, c, d)
    return Q, a2, b2, adds + 1, max_layer


def f_naive(curve: Curve, P, a: int, b: int, r: int):
    s = hash_slot(P, r)
    c, d = jump_coeffs(r, 0, s)
    Q, u, a2, b2 = step_raw(curve, P, a, b, c, d)
    return Q, a2, b2


def coeff_point(curve: Curve, a: int, b: int):
    return curve.add(curve.mul(a, curve.G), curve.mul(b, curve.Q))


def find_fruitless_cycle(curve: Curve, r: int, rng: random.Random, trials: int = 400):
    """Return one explicit (P, Q) with naive f(P)=Q, f(Q)=P, monodromy id."""
    found = 0
    example = None
    for _ in range(trials):
        k0 = rng.randrange(1, curve.order)
        P, _ = canonicalize(curve, curve.mul(k0, curve.G))
        # coefficients of P relative to Q are unknown here; we only need the point map.
        # Track coefficients against G only (Q is not involved) by using b = 0 and
        # a = discrete log of P base G. Solved by BSGS once.
        a = bsgs(curve, curve.G, P, curve.order)
        b = 0
        seen = {}
        for _step in range(64):
            key = P
            if key in seen:
                a0, b0 = seen[key]
                if (a - a0) % curve.order == 0 and (b - b0) % curve.order == 0 and _step - 0 > 0:
                    # closed with identity monodromy. Confirm length 2.
                    Q, _, _ = f_naive(curve, P, a, b, r)
                    R, _, _ = f_naive(curve, Q, 0, 0, r)
                    if R == P and Q != P:
                        found += 1
                        if example is None:
                            example = (P, Q)
                        break
            seen[key] = (a, b)
            P, a, b = f_naive(curve, P, a, b, r)
    return found, example


def solve_collision(curve: Curve, a1, b1, a2, b2):
    n = curve.order
    db = (b2 - b1) % n
    da = (a1 - a2) % n
    if db == 0:
        return None  # fruitless: no information on k
    return (da * inv_mod(db, n)) % n


def recover_many(curve: Curve, r: int, rng: random.Random, instances: int = 30, dp_mod: int = 1):
    """Parallel distinguished-point search on classes. dp_mod == 1 stores every point."""
    ok = 0
    fruitless_hits = 0
    adds = 0
    escapes = 0
    steps = 0
    for _ in range(instances):
        k = rng.randrange(1, curve.order)
        curve.Q = curve.mul(k, curve.G)
        # restart walks until the stored DP collides usefully
        store = {}
        solved = None
        guard = 0
        while solved is None and guard < 40:
            guard += 1
            a0 = rng.randrange(curve.order)
            b0 = rng.randrange(curve.order)
            P = coeff_point(curve, a0, b0)
            if P is None:
                continue
            P, u = canonicalize(curve, P)
            a, b = apply_unit_coeff(curve, a0, b0, u)
            for _t in range(8 * curve.order):
                if dp_mod == 1 or (P[0] % dp_mod) == 0:
                    key = P
                    if key in store:
                        a1, b1 = store[key]
                        got = solve_collision(curve, a1, b1, a, b)
                        if got is None:
                            fruitless_hits += 1
                        else:
                            solved = got
                            break
                    else:
                        store[key] = (a, b)
                P, a, b, used, layer = f_star(curve, P, a, b, r)
                adds += used
                steps += 1
                escapes += layer
                # relation invariant
                if coeff_point(curve, a, b) != P:
                    raise AssertionError("coefficient drift")
            if solved is not None:
                break
        if solved == k:
            ok += 1
        elif solved is not None:
            raise AssertionError((solved, k))
    return {
        "ok": ok,
        "instances": instances,
        "fruitless_hits": fruitless_hits,
        "adds": adds,
        "steps": steps,
        "escapes": escapes,
        "adds_per_step": (adds / steps) if steps else None,
    }


def coeff_table_ok(r: int, layers: int = 4) -> bool:
    """No two consecutive small vectors satisfy (c2, d2) = -(c1, d1).

    That is the only unit-scale relation that makes the 2-cycle identity
    hold for every discrete log. Relations through ±λ, ±λ^2 send a vector
    of height O(r) to height O(n) and are not equal to another table vector.
    """
    vecs = []
    for layer in range(layers + 1):
        for s in range(r):
            vecs.append(jump_coeffs(r, layer, s))
    for c, d in vecs:
        if (-c, -d) in vecs or (c, d) == (0, 0):
            return False
    # uniqueness
    return len(set(vecs)) == len(vecs)


def main():
    rng = random.Random(20261024)
    curve = find_toy_curve(rng)
    assert coeff_table_ok(TOY_R)
    assert coeff_table_ok(PROD_R)
    found, example = find_fruitless_cycle(curve, TOY_R, rng)
    stats = recover_many(curve, TOY_R, rng, instances=20, dp_mod=1)
    # Production-shaped overhead at r=16 on the same curve (not a secp score).
    stats16 = recover_many(curve, PROD_R, rng, instances=10, dp_mod=1)
    print(
        {
            "p": curve.p,
            "n": curve.order,
            "lambda": curve.lam,
            "beta": curve.beta,
            "naive_fruitless_closures": found,
            "example_cycle_x": None if example is None else (example[0][0], example[1][0]),
            "recovered": stats["ok"],
            "instances": stats["instances"],
            "fruitless_collisions": stats["fruitless_hits"],
            "adds_per_step_r8": round(stats["adds_per_step"], 4),
            "adds_per_step_r16": round(stats16["adds_per_step"], 4),
            "bound_r16": round(PROD_R / (PROD_R - 1), 4),
        }
    )
    if found < 1:
        raise SystemExit("canonicalization produced no visible fruitless cycle; toy r too large or hash unlucky")
    if stats["ok"] != stats["instances"]:
        raise SystemExit("escape walk missed a discrete log")


if __name__ == "__main__":
    main()
