"""Executable reference for the tag-tracing 8-adding distinguished-point rho.

Python is NOT the priced machine (machine.md, ttstep.c, tower.c).  The same
`solve` routine runs on BN254 GT at full width (64-bit word arithmetic of the
tag sum emulated exactly) and on toy prime-field subgroups.

Element contexts provide: p (coordinate field), L (64-bit limbs per
coordinate), k (coordinates), r (group order), g, identity, mul, power,
coords(Y) -> k standard residues in the stored order, row(T) -> v with
coord0(Y*T) = sum_k coords(Y)[k]*v[k] mod p, and coord0(Y).
"""
import math

M64 = (1 << 64) - 1
RHO = 3                      # selector bits; R = 8 multipliers
R = 1 << RHO


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def step_cap(order):
    """S = floor(sqrt(1409*order/1000))+1, so S^2/(2*kappa*order) >= 0.6164375
    for kappa = 8/7: model collision probability 1-exp(-0.6164375) > 0.46."""
    return math.isqrt(1409 * order // 1000) + 1


class Entry:
    __slots__ = ("T", "F", "alpha", "beta", "depth", "children", "key")


class Table:
    """All multisets of at most `ell` multipliers, built lazily and memoised
    (production builds the whole table first; laziness changes no value)."""

    def __init__(self, G, mults, ell):
        self.G, self.mults, self.ell = G, mults, ell      # mults: [(M, c, d)]
        self.memo = {}
        self.rinv = pow(1 << (64 * G.L), -1, G.p)
        self.root = self._make((0,) * R, G.identity, 0, 0)

    def _make(self, key, T, alpha, beta):
        G = self.G
        e = Entry()
        e.key, e.T, e.alpha, e.beta = key, T, alpha, beta
        e.depth, e.children = sum(key), [None] * R
        p, L = G.p, G.L
        rinv = self.rinv
        F = []
        for vk in G.row(T):
            for i in range(L):
                u = vk * (1 << (64 * i)) * rinv % p
                f = (u << 128) // p                       # floor(2^128 u/p) < 2^128
                F.append((f & M64, f >> 64))              # (Flo, Fhi)
        e.F = F
        self.memo[key] = e
        return e

    def child(self, e, s):
        if e.depth == self.ell:
            return None
        c = e.children[s]
        if c is None:
            key = e.key[:s] + (e.key[s] + 1,) + e.key[s + 1:]
            c = self.memo.get(key)
            if c is None:
                M, cs, ds = self.mults[s]
                c = self._make(key, self.G.mul(e.T, M), (e.alpha + cs) % self.G.r,
                               (e.beta + ds) % self.G.r)
            e.children[s] = c
        return c


def mont_limbs(G, Y):
    """The stored window base: Montgomery limbs of the k coordinates."""
    out = []
    Rm = 1 << (64 * G.L)
    for c in G.coords(Y):
        m = c * Rm % G.p
        out.extend((m >> (64 * i)) & M64 for i in range(G.L))
    return out


def tag_sum(y, F):
    """A = sum_j (umulh(y_j, Flo_j) + y_j*Fhi_j) mod 2^64  -- the priced kernel."""
    A = 0
    for yj, (lo, hi) in zip(y, F):
        A += ((yj * lo) >> 64) + yj * hi
    return A & M64


def exact_selector(G, t):
    return (t << RHO) // G.p


def solve(G, z, ell, dbits, cbits, cap, walk_cap, walk_len_cap, coins, trace=None):
    """One bounded run.  dbits: a point is distinguished iff t*2^dbits < p.
    cbits: number of top bits of 2^64 t/p the certificate protects
    (production: cbits=3, dbits=80; toys: cbits=dbits<=16).
    Returns (x or None, steps, outcome, stats)."""
    stats = {"exact_paths": 0, "rebases": 0, "walks": 0}
    if z == G.identity:
        return 0, 0, "identity", stats
    r, p = G.r, G.p
    E = 2 * G.k * G.L
    cmask = (1 << (64 - cbits)) - 1
    mults = []
    for _ in range(R):
        c, d = sample_scalar(r, coins), sample_scalar(r, coins)
        mults.append((G.mul(G.power(G.g, c), G.power(z, d)), c, d))
    table = Table(G, mults, ell)
    store = {}
    steps = 0
    for _ in range(walk_cap):
        stats["walks"] += 1
        a, b = sample_scalar(r, coins), sample_scalar(r, coins)
        Y0 = G.mul(G.power(G.g, a), G.power(z, b))
        if Y0 == G.identity:
            continue
        y, e, length, found = mont_limbs(G, Y0), table.root, 0, False
        while True:
            A = tag_sum(y, e.F)
            if ((A + E) & M64) & cmask <= E:
                # not certified: re-base and read the coordinate exactly
                stats["exact_paths"] += 1
                Y0 = G.mul(Y0, e.T); a = (a + e.alpha) % r; b = (b + e.beta) % r
                y, e = mont_limbs(G, Y0), table.root
                t = G.coord0(Y0)
                s, dp = exact_selector(G, t), (t << dbits) < p
            else:
                s = A >> (64 - RHO)
                dp = dbits <= cbits and (A >> (64 - dbits)) == 0
                t = None
            if trace is not None:
                trace(table, Y0, e, s, dp, a, b)
            if dp:
                if t is None:                              # toy-only path (dbits <= cbits)
                    Y0 = G.mul(Y0, e.T); a = (a + e.alpha) % r; b = (b + e.beta) % r
                    y, e = mont_limbs(G, Y0), table.root
                    t = G.coord0(Y0)
                found = True
                break
            if steps >= cap:
                return None, steps, "cap", stats
            if length >= walk_len_cap:
                break
            nxt = table.child(e, s)
            if nxt is None:                                # window end: re-base
                stats["rebases"] += 1
                Y0 = G.mul(Y0, e.T); a = (a + e.alpha) % r; b = (b + e.beta) % r
                y = mont_limbs(G, Y0)
                nxt = table.child(table.root, s)
            e = nxt
            steps += 1
            length += 1
        if not found:
            continue
        hit = store.get(t)
        if hit is not None:
            a2, b2 = hit
            den = (b2 - b) % r
            if den == 0:
                return None, steps, "degenerate", stats
            x = (a - a2) * pow(den, -1, r) % r
            if G.power(G.g, x) == z:
                return x, steps, "success", stats
            continue                                       # equal tags, different elements
        store[t] = (a, b)
    return None, steps, "walk_cap", stats


class ToyField:
    """Order-n subgroup of F_q^*, one coordinate, one limb."""

    def __init__(self, q, n, g):
        assert (q - 1) % n == 0 and pow(g, n, q) == 1 and g != 1 and q < (1 << 62)
        self.p, self.L, self.k, self.r, self.g, self.identity = q, 1, 1, n, g, 1

    def mul(self, a, b):
        return a * b % self.p

    def power(self, base, e):
        return pow(base, e, self.p)

    def coords(self, Y):
        return [Y]

    def row(self, T):
        return [T]

    def coord0(self, Y):
        return Y


class BN254GT:
    """GT in Fp12.  Elements are coefficient tuples in the protected basis
    Fp[w]/(w^12-18w^6+82).  The stored coordinate order is the tower order of
    tower.c: Fp12 = Fp6 + Fp6*w, Fp6 = Fp2 + Fp2*v + Fp2*v^2 (v = w^2),
    Fp2 = Fp + Fp*u (u = w^6-9); coordinate 2*(3h+m)+part is the real
    (part 0) or u (part 1) component of the Fp2 coefficient a_i of w^i,
    i = 2m+h, with a_i = (c_i+9c_{i+6}) + c_{i+6}*u.  The tag coordinate is
    coord0 = c_0, the w^0 coefficient in the protected basis."""

    def __init__(self, p, r, g):
        self.p, self.r, self.g, self.L, self.k = p, r, tuple(g), 4, 12
        self.identity = (1,) + (0,) * 11
        self.order = [(2 * m + h, part) for h in range(2) for m in range(3) for part in range(2)]

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

    def coords(self, Y):
        p = self.p
        return [(Y[i] + 9 * Y[i + 6]) % p if part == 0 else Y[i + 6] for i, part in self.order]

    def coord0(self, Y):
        return Y[0]

    def row(self, T):
        """v with c_0(Y*T) = sum coords(Y)[k]*v[k].  In the protected basis
        c_0(Y*T) = y_0 T_0 - 82 sum_{i=1..11} y_i T_{12-i} - 1476 sum_{i=7..11} y_i T_{18-i}
        (w^12 = 18w^6-82, w^18 = 242w^6-1476); then y_i = Re a_i - 9 Im a_i,
        y_{i+6} = Im a_i for i < 6."""
        p = self.p
        vp = [T[0]] + [(-82 * T[12 - i] - (1476 * T[18 - i] if i >= 7 else 0)) % p
                       for i in range(1, 12)]
        return [vp[i] if part == 0 else (vp[i + 6] - 9 * vp[i]) % p for i, part in self.order]
