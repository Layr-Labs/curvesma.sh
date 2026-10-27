"""Toy-scale descent experiments for the v16 ball-region special-q descent.

Toy tower: h = y^6 + y^3 + y^2 - y - 1 (production), g_toy = x^2 + 4 y and
f_toy = y^2 P(x^2/y) with P(X) = X^2 + X + 1, i.e. f_toy = x^4 + y x^2 + y^2 (degree
4 in x instead of 12, so that norms are small enough for complete enumeration
and factorisation at toy scale).  The production descent's node structure is kept
exactly: a relative-degree-one special-q ideal (q, y - s, x - t) on the f or g side,
its lattice L_q = {(a, b) in Z[y]^2 : a(s) - t b(s) = 0 mod q} (determinant q),
an LLL-reduced basis, the c-space ball {c : ||c||^2 <= R^2, c != 0} with one
representative of each pair {c, -c}, exact norms, the relation predicate (all
prime factors <= 2^lb except at most k window primes in (2^lb, 2^n] per side),
distinct relations keyed by their complete factorisation (unit/sign duplicates
collapse), children = the exact degree-one ideals (q, y - s, x - t) of the window
primes containing the element, memoised recursive descent by ideal with backtracking.

Two experiments:
  fixed_balls : replicated fixed balls on random ideals: distinct relation counts
                against the semismooth model mean (quadrature at the actual norms)
                and their dispersion;
  tree        : complete toy descents from random root ideals to the factor base:
                node outcomes, relations, retries, work, model-predicted means.
Stdlib only.  Evidence at its toy scope only.
"""
import math, functools
from random import Random
from fractions import Fraction as F

H = [-1, -1, 1, 1, 0, 0, 1]
TOY_U = 4
F_TERMS = [(4, [1]), (2, [0, 1]), (0, [0, 0, 1])]          # x^4 + y x^2 + y^2
F_DEG = 4


H_NZ = [(j, H[j]) for j in range(6) if H[j]]  # nonzero low coefficients of h


def hmul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] += x * y
    for i in range(len(c) - 1, 5, -1):
        t = c[i]
        if t:
            for j, hj in H_NZ:
                c[i - 6 + j] -= t * hj
    return (c + [0] * 6)[:6]


def hnorm(a):
    # columns a y^i mod h: multiply the previous column by y and reduce
    cols = [hmul(a, [1])]
    for _ in range(5):
        c = [0] + cols[-1]
        cols.append([c[j] - c[6] * H[j] for j in range(6)])
    m = [[cols[j][i] for j in range(6)] for i in range(6)]
    prev, sign = 1, 1
    for k in range(5):
        if m[k][k] == 0:
            row = next((j for j in range(k + 1, 6) if m[j][k]), None)
            if row is None:
                return 0
            m[k], m[row] = m[row], m[k]; sign = -sign
        for i in range(k + 1, 6):
            for j in range(k + 1, 6):
                v = m[i][j] * m[k][k] - m[i][k] * m[k][j]
                assert v % prev == 0
                m[i][j] = v // prev
        prev = m[k][k]
    return abs(sign * m[-1][-1])


def toy_norms(a, b):
    ap, bp = [[1]], [[1]]
    for _ in range(F_DEG):
        ap.append(hmul(ap[-1], a)); bp.append(hmul(bp[-1], b))
    fs = [0] * 6
    for d, c in F_TERMS:
        term = hmul(c, hmul(ap[d], bp[F_DEG - d]))
        fs = [x + y for x, y in zip(fs, term)]
    gs = [x + y for x, y in zip(ap[2], hmul([0, TOY_U], bp[2]))]
    return hnorm(fs), hnorm(gs)


# ---------------- primes, roots, factoring ----------------
def mr_prime(n):
    if n < 2:
        return False
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % a == 0:
            return n == a
    d = n - 1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def roots_modp(poly, p, rng):
    poly = [c % p for c in poly]
    while poly and poly[-1] == 0:
        poly.pop()
    if len(poly) <= 1:
        return []

    def pmod(a, m):
        a = a[:]
        inv = pow(m[-1], -1, p)
        while len(a) >= len(m):
            if a[-1]:
                z = a[-1] * inv % p; s = len(a) - len(m)
                for j, v in enumerate(m):
                    a[s + j] = (a[s + j] - z * v) % p
            a.pop()
            while a and a[-1] == 0:
                a.pop()
        return a

    def pmul(a, b, m):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] = (r[i + j] + x * y) % p
        return pmod(r, m)

    def powmod(a, e, m):
        r = [1]
        while e:
            if e & 1:
                r = pmul(r, a, m)
            a = pmul(a, a, m); e >>= 1
        return r
    xp = powmod([0, 1], p, poly)
    xp = xp + [0] * max(0, 2 - len(xp)); xp[1] = (xp[1] - 1) % p
    a, b = poly, pmod(xp, poly)
    while b and any(b):
        a, b = b, pmod(a, b)
    if not a or not any(a):
        return []
    inv = pow(a[-1], -1, p); g = [x * inv % p for x in a]

    def split(g):
        if len(g) <= 1:
            return []
        if len(g) == 2:
            return [(-g[0] * pow(g[1], -1, p)) % p]
        while True:
            c = [rng.randrange(p) for _ in range(len(g) - 2)] + [rng.randrange(1, p), 1]
            t = powmod(c, (p - 1) // 2, g)
            t = t + [0] * max(0, 1 - len(t)); t = [(t[0] - 1) % p] + t[1:]
            h = pmod(t, g) if len(t) >= len(g) else t
            aa, bb = g, h
            while bb and any(bb):
                aa, bb = bb, pmod(aa, bb)
            if aa and 1 < len(aa) < len(g):
                iv = pow(aa[-1], -1, p); hh = [x * iv % p for x in aa]
                qq = g[:]; v = [0] * (len(g) - len(hh) + 1)
                for k in range(len(v) - 1, -1, -1):
                    v[k] = qq[k + len(hh) - 1]
                    for j, x in enumerate(hh):
                        qq[k + j] = (qq[k + j] - v[k] * x) % p
                assert not any(qq)
                return split(hh) + split(v)
    return sorted(split(g))


def factorize_window(n, small, lb, nbits, k):
    """exact relation decision for one side: n (special-q already divided out) is
    trial-divided by the primes <= 2^lb; the cofactor must then be a product of at
    most k primes in (2^lb, 2^nbits].  Returns (ok, factors) with the complete
    factorisation when ok (TD by all primes <= 2^nbits, deterministic), else
    (False, None).  Early rejections: cofactor > 2^(k*nbits), or a cofactor that
    is prime and > 2^nbits, or more than k window primes."""
    fs = []
    for p in small:
        if p * p > n:
            break
        while n % p == 0:
            fs.append(p); n //= p
    if n == 1:
        return True, fs
    LB = 1 << lb
    if n <= LB:
        fs.append(n); return True, fs
    if n >> (k * nbits):
        return False, None
    if mr_prime(n):
        if n <= (1 << nbits):
            fs.append(n); return True, fs
        return False, None
    # composite cofactor of primes > 2^lb, at most k*nbits bits: TD by the window primes
    win = 0
    for p in WINDOW_PRIMES(nbits):
        if p <= LB:
            continue
        if p * p > n:
            break
        while n % p == 0:
            fs.append(p); n //= p; win += 1
            if win > k:
                return False, None
    if n > 1:
        if n > (1 << nbits) or not mr_prime(n):
            return False, None
        fs.append(n); win += 1
    if win > k:
        return False, None
    return True, sorted(fs)


_WP = {}
def WINDOW_PRIMES(nbits):
    if nbits not in _WP:
        _WP[nbits] = primes_upto(1 << nbits)
    return _WP[nbits]


# ---------------- LLL (decimal proposals, exact verification) ----------------
def lll(B, prec=120):
    from decimal import Decimal as D, localcontext
    B = [v[:] for v in B]; n = len(B)
    with localcontext() as ctx:
        ctx.prec = prec
        mu = []; norm = []; star = []
        for v in B:
            s = list(map(D, v)); row = []
            for j, b in enumerate(star):
                x = sum(D(a) * c for a, c in zip(v, b)) / norm[j]; row.append(x)
                s = [a - x * c for a, c in zip(s, b)]
            star.append(s); mu.append(row); norm.append(sum(x * x for x in s))
        k = 1
        while k < n:
            for j in range(k - 1, -1, -1):
                q = int(mu[k][j].to_integral_value(rounding='ROUND_HALF_EVEN'))
                if q:
                    B[k] = [x - q * y for x, y in zip(B[k], B[j])]
                    for i in range(j):
                        mu[k][i] -= q * mu[j][i]
                    mu[k][j] -= q
            m = mu[k][k - 1]
            if norm[k] >= (D('0.99') - m * m) * norm[k - 1]:
                k += 1; continue
            A = norm[k] + m * m * norm[k - 1]; newm = m * norm[k - 1] / A
            norm[k] = norm[k] * norm[k - 1] / A; norm[k - 1] = A; mu[k][k - 1] = newm
            B[k], B[k - 1] = B[k - 1], B[k]
            for j in range(k - 1):
                mu[k][j], mu[k - 1][j] = mu[k - 1][j], mu[k][j]
            for i in range(k + 1, n):
                t = mu[i][k]; mu[i][k] = mu[i][k - 1] - m * t; mu[i][k - 1] = t + newm * mu[i][k]
            k = max(1, k - 1)
    return B


def gs_exact(B):
    stars = []; norm = []
    for v in B:
        s = list(map(F, v))
        for j, b in enumerate(stars):
            m = sum(x * y for x, y in zip(v, b)) / norm[j]
            s = [x - m * y for x, y in zip(s, b)]
        stars.append(s); norm.append(sum(x * x for x in s))
    return norm


# ---------------- special-q ideals and their lattices ----------------
def random_ideal(rng, bits, side):
    while True:
        q = rng.randrange(1 << (bits - 1), 1 << bits) | 1
        if not mr_prime(q):
            continue
        opts = []
        for s in roots_modp(H, q, rng):
            if side == 'g':
                poly = [TOY_U * s % q, 0, 1]
            else:
                poly = [s * s % q, 0, s, 0, 1]
            opts.extend((s, t) for t in roots_modp(poly, q, rng))
        if opts:
            s, t = rng.choice(opts)
            return q, s, t


def ideal_lattice(q, s, t):
    w = [pow(s, j, q) for j in range(6)] + [(-t * pow(s, j, q)) % q for j in range(6)]
    B0 = [[q] + [0] * 11] + [[(-w[j]) % q] + [int(k == j) for k in range(1, 12)] for j in range(1, 12)]
    B = lll(B0)
    assert all(sum(a * b for a, b in zip(v, w)) % q == 0 for v in B)
    assert math.prod(gs_exact(B)) == q * q
    return B


def ball_cells(R):
    """all integer c in Z^12 with 0 < ||c||^2 <= R^2 and last nonzero coordinate
    positive (one representative of each pair {c,-c}); lexicographic order."""
    out = []
    c = [0] * 12

    def rec(d, s, positive):
        if d < 0:
            if positive:
                out.append(tuple(c))
            return
        r = int(math.isqrt(s))
        lo = 0 if not positive else -r
        for x in range(lo, r + 1):
            if not positive and x == 0:
                c[d] = 0; rec(d - 1, s, False)
            else:
                if not positive and x < 0:
                    continue
                c[d] = x; rec(d - 1, s - x * x, True)
    rec(11, R * R, False)
    return out


# ---------------- the node: enumerate, factor, extract relations ----------------
def side_poly(s, q, side):
    return [TOY_U * s % q, 0, 1] if side == 'g' else [s * s % q, 0, s, 0, 1]


def child_ideals(a, b, q, side):
    """the degree-one prime ideals (q, y - s, x - t) on `side` containing a - b x:
    a(s) - t b(s) = 0 mod q (the exact ideals the element factors through)."""
    rng = Random(q)
    ev = lambda p, z: sum(c * pow(z, j, q) for j, c in enumerate(p)) % q
    out = []
    for s in sorted(set(roots_modp(H, q, rng))):
        for t in sorted(set(roots_modp(side_poly(s, q, side), q, rng))):
            if (ev(a, s) - t * ev(b, s)) % q == 0:
                out.append((q, s, t))
    return out


def node_relations(q, s, t, side, B, cells, lb, n, k, small, with_model=None):
    """returns dict: cells, raw relations, distinct relations (keyed by both
    factorisations), children per relation, unfactored, model mean (if with_model)."""
    LB = 1 << lb; NB = 1 << n
    raw = 0; distinct = {}; junk = 0
    model_sum = 0.0
    for c in cells:
        v = [sum(c[i] * B[i][j] for i in range(12)) for j in range(12)]
        nf, ng = toy_norms(v[:6], v[6:])
        if nf == 0 or ng == 0:
            continue
        if side == 'f':
            assert nf % q == 0; nf //= q
        else:
            assert ng % q == 0; ng //= q
        if with_model is not None:
            model_sum += with_model(nf, ng)
        okf, fs = factorize_window(nf, small, lb, n, k)
        if not okf:
            junk += 1; continue
        okg, gs = factorize_window(ng, small, lb, n, k)
        if not okg:
            junk += 1; continue
        wf = [x for x in fs if x > LB]; wg = [x for x in gs if x > LB]
        raw += 1
        key = (tuple(fs), tuple(gs))
        if key not in distinct:
            cf = [i for x in sorted(set(wf)) for i in child_ideals(v[:6], v[6:], x, 'f')]
            cg = [i for x in sorted(set(wg)) for i in child_ideals(v[:6], v[6:], x, 'g')]
            if len(cf) != len(set(wf)) or len(cg) != len(set(wg)):
                junk += 1; continue      # not one degree-one ideal per window prime
            distinct[key] = dict(children_f=cf, children_g=cg, cell=c)
    return dict(cells=len(cells), raw=raw, distinct=len(distinct), duplicates=raw - len(distinct),
                junk=junk, model_mean=model_sum, relations=distinct)


def fixed_balls(seed, bits, side, R, lb, n, k, ideals, model=None):
    """replicated fixed balls: one list per ideal with the distinct relation count
    and the model mean (quadrature at the actual norms when `model` is given)."""
    rng = Random(seed)
    small = primes_upto(1 << lb)
    cells = ball_cells(R)
    rows = []
    for i in range(ideals):
        q, s, t = random_ideal(rng, bits, side)
        B = ideal_lattice(q, s, t)
        r = node_relations(q, s, t, side, B, cells, lb, n, k, small, model)
        rows.append(dict(q=q, s=s, t=t, cells=r['cells'], raw=r['raw'], distinct=r['distinct'],
                         duplicates=r['duplicates'], junk=r['junk'], model_mean=r['model_mean']))
    return rows


class Tree:
    """memoised recursive descent with backtracking at toy scale."""

    def __init__(self, seed, lb, plan, small=None, max_nodes=400, order='random'):
        self.rng = Random(seed); self.lb = lb; self.plan = plan; self.order = order
        self.small = small or primes_upto(1 << lb)
        self.memo = {}; self.nodes = []; self.cells = 0; self.max_nodes = max_nodes
        self.cell_cache = {}

    def cells_for(self, R):
        if R not in self.cell_cache:
            self.cell_cache[R] = ball_cells(R)
        return self.cell_cache[R]

    def descend(self, ideal, side, depth=0):
        """descend the prime ideal (q, y - s, x - t) on `side` (memoised by ideal)."""
        q, s, t = ideal
        key = (q, s, t, side)
        if key in self.memo:
            return self.memo[key]
        if len(self.nodes) >= self.max_nodes:
            return False
        bits = q.bit_length()
        R, n, k = self.plan(bits)
        B = ideal_lattice(q, s, t)
        cells = self.cells_for(R)
        r = node_relations(q, s, t, side, B, cells, self.lb, n, k, self.small)
        self.cells += len(cells)
        rels = list(r['relations'].values())
        if self.order == 'fewest':
            rels.sort(key=lambda d: (len(d['children_f']) + len(d['children_g'])))
        elif self.order == 'children_first':
            # demonstration order: relations with children first, so the tree recurses
            self.rng.shuffle(rels)
            rels.sort(key=lambda d: not (d['children_f'] or d['children_g']))
        else:
            self.rng.shuffle(rels)
        tried = 0; ok = False
        for rel in rels:
            tried += 1
            good = True
            for cq in rel['children_f']:
                if not self.descend(cq, 'f', depth + 1):
                    good = False; break
            if good:
                for cq in rel['children_g']:
                    if not self.descend(cq, 'g', depth + 1):
                        good = False; break
            if good:
                ok = True; break
        zero_children = sum(1 for d in r['relations'].values() if not d['children_f'] and not d['children_g'])
        self.nodes.append(dict(q=q, side=side, bits=bits, depth=depth, R=R, n=n, k=k, cells=len(cells), zero_children=zero_children,
                               raw=r['raw'], distinct=r['distinct'], duplicates=r['duplicates'],
                               junk=r['junk'], tried=tried, success=ok))
        self.memo[key] = ok
        return ok


