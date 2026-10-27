"""Toy-only arithmetic helpers adapted from the earlier submissions.
Executed only inside the pinned verifier container. No imported score claim.
"""
import math
import random

def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    e, s = n - 1, 0
    while e % 2 == 0:
        e //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, e, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
    mm, c, t, x = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (mm - i - 1), p)
        mm, c, t, x = i, b * b % p, t * b * b % p, x * b % p
    return x


def toy_prime(bits, dd, rng):
    while True:
        rr = rng.randrange(1 << (bits - 1), 1 << bits) // (6 * dd) * (6 * dd) + 1
        if is_prime(rr) and ((rr - 1) // dd) % 6 == 0 and pow(2, (rr - 1) // dd, rr) != 1:
            return rr


class Curve:
    def __init__(self, p, c):
        self.p, self.c = p, c
        self.inv = self.muls = 0

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
        else:
            lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def mul(self, k, P):
        R = None
        while k:
            if k & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            k >>= 1
        return R


def find_curve(bits, dd, rng):
    """p = 1 mod 3, y^2 = x^3 + c; prime group order r with 6d | r-1 and 2 not in H (as in production)."""
    while True:
        p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if p % 3 != 1 or not is_prime(p):
            continue
        for c in range(1, 8):
            E = Curve(p, c)
            while True:
                x = rng.randrange(p)
                y = sqrt_mod((x ** 3 + c) % p, p)
                if y is not None:
                    P = (x, y)
                    break
            lo, hi = p + 1 - 2 * math.isqrt(p) - 2, p + 1 + 2 * math.isqrt(p) + 2
            w = math.isqrt(hi - lo) + 1
            baby, R = {}, None
            for j in range(w):
                if R is not None:
                    baby.setdefault(R, j)
                R = E.add(R, P)
            step = E.mul(w, P)
            step = (step[0], (-step[1]) % p) if step else None
            Gk, cands = E.mul(lo, P), []
            for i in range(w + 1):
                if Gk is None:
                    cands.append(lo + i * w)
                elif (Gk[0], (-Gk[1]) % p) in baby:
                    cands.append(lo + i * w + baby[(Gk[0], (-Gk[1]) % p)])
                Gk = E.add(Gk, step)
            cands = sorted({N for N in cands if lo <= N <= hi and E.mul(N, P) is None})
            if len(cands) == 1 and is_prime(cands[0]) and (cands[0] - 1) % (6 * dd) == 0 \
                    and pow(2, (cands[0] - 1) // dd, cands[0]) != 1:
                return E, cands[0], P


def dth_root(t, dd, rr):
    k = dd.bit_length() - 1
    for j in range(k, 0, -1):
        y = sqrt_mod(t, rr)
        e = (rr - 1) >> (j - 1)
        t = next(z for z in (y, rr - y) if j == 1 or pow(z, e, rr) == 1)
    return t


