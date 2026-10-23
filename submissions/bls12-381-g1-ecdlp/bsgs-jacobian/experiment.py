"""Toy exhaustive check of the submitted BSGS construction.

Scope: small Weierstrass groups y^2 = x^3 + b over a small prime field.  For the
prime-order subgroup <P> of order r, every s in 1..r-1 is attacked with the same
construction as proof.md (baby walk, giant walk with D = -[m]P, sort both lists
by x, linear merge, infinity branch, +- false-match branch, verification) and the
recovered scalar is compared with s.  This is evidence at toy scope only: it does
not establish production curve behaviour, memory bounds or operation costs.
"""

import json
import sys
from math import isqrt


def inv(a, p):
    return pow(a, p - 2, p)


def add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1) * inv(2 * y1 % p, p) % p
    else:
        lam = (y2 - y1) * inv((x2 - x1) % p, p) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def mul(k, P, p):
    R = None
    Q = P
    while k:
        if k & 1:
            R = add(R, Q, p)
        Q = add(Q, Q, p)
        k >>= 1
    return R


def neg(P, p):
    if P is None:
        return None
    return (P[0], (-P[1]) % p)


def curve_order(p, b):
    total = 1
    for x in range(p):
        rhs = (x * x * x + b) % p
        if rhs == 0:
            total += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            total += 2
    return total


def largest_prime_factor(n):
    best = 1
    d = 2
    while d * d <= n:
        while n % d == 0:
            best = max(best, d)
            n //= d
        d += 1
    return max(best, n)


def find_point(p, b, r):
    sq = {}
    for y in range(p):
        sq.setdefault(y * y % p, y)
    for x in range(p):
        rhs = (x * x * x + b) % p
        y = sq.get(rhs)
        if y is None:
            continue
        Pt = (x, y)
        if Pt != (x, 0) and mul(r, Pt, p) is None:
            return Pt
    return None


def attack(P, T, r, p):
    m = isqrt(r)
    if m * m < r:
        m += 1
    h = m
    assert m * h >= r

    baby = []
    cur = P
    for j in range(1, m):
        if cur is not None:
            baby.append((cur[0], j))
        cur = add(cur, P, p)

    D = neg(mul(m, P, p), p)
    giant = []
    inf = []
    cur = T
    for i in range(h):
        if cur is None:
            inf.append(i)
        else:
            giant.append((cur[0], i))
        cur = add(cur, D, p)

    baby.sort()
    giant.sort()

    cands = []
    for i in inf:
        cands.append((i * m) % r)
    a = b = 0
    while a < len(baby) and b < len(giant):
        xa = baby[a][0]
        xb = giant[b][0]
        if xa < xb:
            a += 1
        elif xa > xb:
            b += 1
        else:
            a2 = a
            while a2 < len(baby) and baby[a2][0] == xa:
                a2 += 1
            b2 = b
            while b2 < len(giant) and giant[b2][0] == xa:
                b2 += 1
            for t in range(a, a2):
                for u in range(b, b2):
                    cands.append((giant[u][1] * m + baby[t][1]) % r)
            a = a2
            b = b2
    for cand in cands:
        if mul(cand, P, p) == T:
            return cand
    return None


def run(prime):
    p = prime
    best = None
    for b in (7, 11, 13, 17, 19, 23):
        n = curve_order(p, b)
        r = largest_prime_factor(n)
        if r > 16 and (best is None or r > best[1]):
            best = (b, r)
    if best is None:
        return {"prime": p, "order": 0, "checked": 0, "ok": False, "note": "no usable prime-order subgroup"}
    b, r = best
    P = find_point(p, b, r)
    if P is None:
        return {"prime": p, "order": r, "checked": 0, "ok": False, "note": "no point of order r"}
    checked = 0
    for s in range(1, r):
        T = mul(s, P, p)
        got = attack(P, T, r, p)
        if got != s:
            return {"prime": p, "order": r, "checked": checked, "ok": False, "first_failure": s}
        checked += 1
    return {"prime": p, "order": r, "checked": checked, "ok": True}


def main():
    raw = sys.stdin.read() or "{}"
    data = json.loads(raw)
    prime = int(data.get("prime", 1009))
    if not (101 <= prime <= 4096):
        raise SystemExit("prime out of range")
    print(json.dumps(run(prime), sort_keys=True))


if __name__ == "__main__":
    main()
