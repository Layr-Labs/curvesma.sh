#!/usr/bin/env python3
"""Standalone polynomial/tower arithmetic for supporting tests, not production code.
Fp12=Fp[w]/(w^12-2*w^6+2), xi=1+u=w^6, u^2=-1.
"""
class Tower:
    def __init__(self, p):
        assert p % 12 == 7
        self.p = p
        self.one = (1,) + (0,)*11
        self.xi = (1, 1)
        self.d = self.f2pow(self.xi, (p-1)//6)
        self.zeta = self.f2pow(self.xi, (p*p-1)//6)
        z, zi = self.zeta
        assert zi == 0 and (z*z-z+1) % p == 0
        self.factors = [[self.f2pow(self.xi, i*((p**k-1)//6))
                         for i in range(6)] for k in range(12)]

    def f2mul(self, a, b):
        p = self.p
        return ((a[0]*b[0]-a[1]*b[1]) % p,
                (a[0]*b[1]+a[1]*b[0]) % p)

    def f2pow(self, a, n):
        out = (1, 0)
        while n:
            if n & 1:
                out = self.f2mul(out, a)
            a = self.f2mul(a, a)
            n >>= 1
        return out

    def mul(self, a, b):
        t = [0]*23
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                t[i+j] += x*y
        for i in range(22, 11, -1):
            v = t[i]
            t[i-6] += 2*v
            t[i-12] -= 2*v
        p = self.p
        return tuple(v % p for v in t[:12])

    def pow(self, a, n):
        out = self.one
        while n:
            if n & 1:
                out = self.mul(out, a)
            n >>= 1
            if n:
                a = self.mul(a, a)
        return out

    def frob(self, a, k):
        p = self.p
        k %= 12
        out = [0]*12
        for i, f in enumerate(self.factors[k]):
            c = ((a[i]+a[i+6]) % p, a[i+6])
            if k & 1:
                c = (c[0], -c[1] % p)
            c = self.f2mul(c, f)
            out[i], out[i+6] = (c[0]-c[1]) % p, c[1]
        return tuple(out)

    def keys(self, y):
        # Four products: a*d0, b*d1, a*zeta, B*zeta.
        p = self.p
        a, b = (y[1]+y[7]) % p, y[7]
        B = (a*self.d[0]+b*self.d[1]) % p
        az, bz = a*self.zeta[0] % p, B*self.zeta[0] % p
        even = (a, az, (az-a) % p, -a % p, -az % p, (a-az) % p)
        odd = (B, bz, (bz-B) % p, -B % p, -bz % p, (B-bz) % p)
        return tuple(v for pair in zip(even, odd) for v in pair)

    def metadata(self, y):
        ks = self.keys(y)
        K = min(ks)
        ids = [k for k in range(12) if ks[k] == K]
        if len(ids) == 1:
            k = ids[0]
            return K, k, self.frob(y, k)
        # Toy fields need frequent exact ties. The production fast path forms
        # a full conjugate only on a tie or a class-equality candidate.
        rep, k = min((self.frob(y, k), k) for k in ids)
        return K, k, rep


def exact_order(a, r, order):
    return pow(a, order, r) == 1 and all(
        pow(a, order//q, r) != 1 for q in (2, 3) if order % q == 0)
