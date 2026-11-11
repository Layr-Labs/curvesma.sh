"""Executable mathematical reference for the Frobenius-class distinguished-point walk.

Python is NOT the priced machine (see machine.md). The same `solve` routine is
used for toy additive groups and for full-width BN254 GT arithmetic.
"""
import hashlib

SELECT = (0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11)  # Frobenius powers j != 6


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def digest(salt, inv):
    return hashlib.sha256(salt + inv).digest()


def selector(h):
    """Index j in SELECT from the first four digest bytes."""
    return SELECT[int.from_bytes(h[0:4], "big") % 11]


def distinguished(h, theta_bits):
    return int.from_bytes(h[4:20], "big") & ((1 << theta_bits) - 1) == 0


def step(G, salt, Y, C):
    """One walk step Y <- Y^2 * frob^j(Y); C <- C*(2+rho^j) mod r."""
    h = digest(salt, G.invariant(Y))
    j = selector(h)
    return G.mul(G.mul(Y, Y), G.frob(Y, j)), C * G.cj[j] % G.r


def solve(G, z, theta_bits, step_cap, walk_cap, walk_len_cap, coins):
    """Returns (x or None, steps, outcome). Never reads a hidden scalar."""
    if z == G.identity:
        return 0, 0, "identity"
    r = G.r
    salt = coins.getrandbits(128).to_bytes(16, "big")
    table = {}
    steps = 0
    for _ in range(walk_cap):
        a, b = sample_scalar(r, coins), sample_scalar(r, coins)
        Y = G.mul(G.power(G.g, a), G.power(z, b))
        if Y == G.identity:
            continue
        C, length, found = 1, 0, False
        while True:
            inv = G.invariant(Y)
            h = digest(salt, inv)
            if distinguished(h, theta_bits):
                found = True
                break
            if steps >= step_cap:
                return None, steps, "cap"
            if length >= walk_len_cap:
                break
            j = selector(h)
            Y = G.mul(G.mul(Y, Y), G.frob(Y, j))
            C = C * G.cj[j] % r
            steps += 1
            length += 1
        if not found:
            continue
        for (Y2, C2, a2, b2) in table.get(inv, ()):
            for k in range(12):
                if G.frob(Y2, k) == Y:
                    # C*(a+b x) = rho^k * C2*(a2+b2 x)  (mod r)
                    pk = pow(G.rho, k, r)
                    den = (C * b - pk * C2 * b2) % r
                    if den == 0:
                        return None, steps, "degenerate"
                    x = (pk * C2 * a2 - C * a) * pow(den, -1, r) % r
                    if G.power(G.g, x) != z:
                        return None, steps, "verification_failed"
                    return x, steps, "success"
        table.setdefault(inv, []).append((Y, C, a, b))
    return None, steps, "walk_cap"


class ToyGroup:
    """Additive group Z/q, q prime, q = 1 mod 12, with automorphism y -> rho*y
    of exact order 12 (the analogue of Frobenius on GT). Generator 1."""

    def __init__(self, q, rho):
        assert pow(rho, 12, q) == 1 and pow(rho, 6, q) == q - 1
        assert pow(rho, 4, q) != 1
        self.r, self.rho, self.g, self.identity = q, rho, 1, 0
        self.pw = [pow(rho, k, q) for k in range(12)]
        self.cj = [(2 + c) % q for c in self.pw]

    def mul(self, a, b):
        return (a + b) % self.r

    def power(self, base, e):
        return base * e % self.r

    def frob(self, y, j):
        return y * self.pw[j] % self.r

    def invariant(self, y):
        q = self.r
        return min(y * c % q for c in self.pw).to_bytes(32, "big")


class Fp12Group:
    """BN254 GT in the protected basis Fp[w]/(w^12-18w^6+82)."""

    def __init__(self, p, r, g):
        self.p, self.r, self.g = p, r, tuple(g)
        self.identity = (1,) + (0,) * 11
        self.rho = p % r
        self.cj = [(2 + pow(p, j, r)) % r for j in range(12)]
        w = (0, 1) + (0,) * 10
        # w^(p^j) = w * gamma_j with gamma_j in Fp[w^6]; images of w^i by powers.
        self.frob_basis = []
        wj = w
        for j in range(12):
            row, acc = [], self.identity
            for i in range(12):
                row.append(acc)
                acc = self.mul(acc, wj)
            self.frob_basis.append(row)
            wj = self.power(wj, p)
        assert wj == w

    def mul(self, a, b):
        p = self.p
        t = [0] * 23
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    t[i + j] += ai * bj
        for k in range(22, 11, -1):
            v = t[k] % p
            t[k - 6] += 18 * v
            t[k - 12] -= 82 * v
        return tuple(v % p for v in t[:12])

    def power(self, base, e):
        out = self.identity
        for bit in bin(e)[2:]:
            out = self.mul(out, out)
            if bit == "1":
                out = self.mul(out, base)
        return out

    def frob(self, y, j):
        p = self.p
        row = self.frob_basis[j % 12]
        acc = [0] * 12
        for ci, img in zip(y, row):
            if ci:
                for k in range(12):
                    acc[k] += ci * img[k]
        return tuple(v % p for v in acc)

    def trace_twelfth(self, y):
        """Tr_{Fp12/Fp}(y)/12 = c0 + 9*c6 = real part of the tower constant term."""
        return (y[0] + 9 * y[6]) % self.p

    def invariant(self, y):
        return self.trace_twelfth(y).to_bytes(32, "big")
