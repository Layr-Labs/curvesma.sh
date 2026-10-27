"""Production-field boot and special-q lattice witnesses (stdlib only).
All arithmetic, membership and norm checks are exact; the decimal LLL proposes a
basis and an independent exact Gram-Schmidt verifies the LLL conditions and the
determinant.  v16: sample points are drawn from c-space BALLS (one of each +-pair,
last nonzero coordinate positive), matching the descent regions of DESCENT.md.
"""
import os, sys, json, math, time, random
from fractions import Fraction as F
from decimal import Decimal as D, localcontext
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'files') if os.path.isdir(os.path.join(HERE, 'files')) else HERE)
import common12 as C
P = 4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787
R = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
U = C.U
LOG2U6 = 6 * math.log2(U)
Z = [0] * 6; ONE = ([1] + [0] * 5, Z[:])


def add(a, b): return [(x + y) % P for x, y in zip(a, b)]
def neg(a): return [(-x) % P for x in a]
def mul(a, b): return [x % P for x in C.mul(a, b)]


def pow6(a, n):
    v = ONE[0]
    while n:
        if n & 1: v = mul(v, a)
        a = mul(a, a); n >>= 1
    return v


def mul12(a, b):
    ac = mul(a[0], b[0]); bd = mul(a[1], b[1]); cross = mul(add(a[0], a[1]), add(b[0], b[1]))
    return add(ac, neg(mul([0, U, 0, 0, 0, 0], bd))), add(cross, neg(add(ac, bd)))


def pow12(a, n):
    v = ONE
    while n:
        if n & 1: v = mul12(v, a)
        a = mul12(a, a); n >>= 1
    return v


def boot_basis(z):
    if not any(z[1]): raise ValueError('subfield/identity: solve separately')
    ratio = mul(neg(z[0]), pow6(z[1], P ** 6 - 2))
    return [[P if i == j else 0 for j in range(12)] for i in range(6)] + [mul(ratio, [0] * j + [1]) + [int(i == j) for i in range(6)] for j in range(6)]


def gs_exact(B):
    stars = []; mu = []; norm = []
    for v in B:
        row = []; s = list(map(F, v))
        for j, b in enumerate(stars):
            m = sum(x * y for x, y in zip(v, b)) / norm[j]; row.append(m)
            s = [x - m * y for x, y in zip(s, b)]
        stars.append(s); mu.append(row); norm.append(sum(x * x for x in s))
    return mu, norm


def lll(B):
    B = [v[:] for v in B]; n = len(B); swaps = 0
    with localcontext() as ctx:
        ctx.prec = 300
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
                    for i in range(j): mu[k][i] -= q * mu[j][i]
                    mu[k][j] -= q
            m = mu[k][k - 1]
            if norm[k] >= (D('0.99') - m * m) * norm[k - 1]: k += 1; continue
            A = norm[k] + m * m * norm[k - 1]; newm = m * norm[k - 1] / A
            norm[k] = norm[k] * norm[k - 1] / A; norm[k - 1] = A; mu[k][k - 1] = newm
            B[k], B[k - 1] = B[k - 1], B[k]
            for j in range(k - 1): mu[k][j], mu[k - 1][j] = mu[k - 1][j], mu[k][j]
            for i in range(k + 1, n):
                t = mu[i][k]; mu[i][k] = mu[i][k - 1] - m * t; mu[i][k - 1] = t + newm * mu[i][k]
            k = max(1, k - 1); swaps += 1
    return B, swaps


def certify_boot(z, B):
    for v in B:
        assert mul(v[:6], z[1]) == neg(mul(v[6:], z[0]))
    mu, norm = gs_exact(B)
    for i in range(12):
        assert all(abs(x) <= F(1, 2) for x in mu[i])
        if i: assert norm[i] >= (F(98, 100) - mu[i][-1] ** 2) * norm[i - 1]
    det2 = math.prod(norm); assert det2 == P ** 12
    return float(sum(norm)), mu, norm


def ball_point(rng, Rr, dim=12):
    """integer point of the ball ||c||^2 <= Rr^2 (uniform real point rounded; c != 0;
    last nonzero coordinate positive, i.e. one representative of {c, -c})."""
    while True:
        g = [rng.gauss(0, 1) for _ in range(dim)]
        nrm = math.sqrt(sum(x * x for x in g))
        r = Rr * rng.random() ** (1.0 / dim)
        c = [int(round(x * r / nrm)) for x in g]
        if any(c) and sum(x * x for x in c) <= Rr * Rr:
            for x in reversed(c):
                if x:
                    if x < 0: c = [-y for y in c]
                    break
            return c


def boot_ball_stats(z, B, seed, points=16, radii=(14, 96)):
    """exact g-side norms of ball points of the certified boot lattice: cg offsets
    cg = log2 N_g - 12 c - 6 log2|u|, c = (6 log2 p + S)/12, S = log2 of the half count."""
    import descent16 as DS
    cert, mu, norm = certify_boot(z, B); rng = random.Random(seed)
    out = []
    for Rr in radii:
        S = math.log2(DS.CELLS[Rr]); c = (6 * math.log2(P) + S) / 12
        cg = []
        for _ in range(points):
            k = ball_point(rng, Rr)
            v = [sum(k[i] * B[i][j] for i in range(12)) for j in range(12)]
            assert mul(v[:6], z[1]) == neg(mul(v[6:], z[0]))
            ng = C.norm([a + b for a, b in zip(C.mul(v[:6], v[:6]), C.mul([0, U], C.mul(v[6:], v[6:])))])
            assert ng > 0
            cg.append(math.log2(ng) - 12 * c - LOG2U6)
        out.append(dict(R=Rr, S=S, points=points, cg_mean=sum(cg) / points, cg_max=max(cg), cg_min=min(cg)))
    return out


# finite-field polynomial tools for actual relative-degree-one special-q ideals
def trim(a):
    while a and a[-1] == 0: a.pop()
    return a


def psub(a, b, q): return trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % q for i in range(max(len(a), len(b)))])


def pmod(a, b, q):
    a = trim([x % q for x in a]); iv = pow(b[-1], -1, q)
    while len(a) >= len(b):
        v = a[-1] * iv % q; k = len(a) - len(b)
        for i, x in enumerate(b): a[k + i] = (a[k + i] - v * x) % q
        trim(a)
    return a


def pmul(a, b, f, q):
    v = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): v[i + j] = (v[i + j] + x * y) % q
    return pmod(v, f, q)


def ppow(a, n, f, q):
    v = [1]
    while n:
        if n & 1: v = pmul(v, a, f, q)
        a = pmul(a, a, f, q); n >>= 1
    return v


def pgcd(a, b, q):
    while b: a, b = b, pmod(a, b, q)
    if not a: return []
    iv = pow(a[-1], -1, q); return [x * iv % q for x in a]


def roots(f, q, rng):
    f = trim([x % q for x in f]); f = pgcd(f, psub(ppow([0, 1], q, f, q), [0, 1], q), q)

    def split(g):
        if len(g) <= 1: return []
        if len(g) == 2: return [-g[0] * pow(g[1], -1, q) % q]
        while True:
            a = [rng.randrange(q) for _ in range(len(g) - 1)]
            h = pgcd(g, psub(ppow(a, (q - 1) // 2, g, q), [1], q), q)
            if 1 < len(h) < len(g):
                a = g[:]; v = [0] * (len(g) - len(h) + 1)
                for k in range(len(v) - 1, -1, -1):
                    v[k] = a[k + len(h) - 1]
                    for j, x in enumerate(h): a[k + j] = (a[k + j] - v[k] * x) % q
                assert not any(a)
                return split(h) + split(trim(v))
    return sorted(split(f))


def probable_prime(n):
    if n % 2 == 0: return False
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        if n == a: return True
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


def specialq_lattice(bits, side, rng):
    """an actual relative-degree-one special-q ideal (q, y - s, x - t) of `bits`
    bits on `side`, its lattice basis (certified: congruences and Gram determinant q^2)."""
    while True:
        q = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not probable_prime(q): continue
        ss = roots(C.H, q, rng)
        opts = []
        for s in ss:
            if side == 'g': ff = [U * s, 0, 1]
            else: ff = [1 + s - s * s - s ** 3, 0, s ** 5, 0, 0, 0, 2 * s ** 3, 0, 0, 0, -2 * s, 0, 1]
            opts.extend((s, t) for t in roots(ff, q, rng))
        if opts: break
    s, t = rng.choice(opts)
    w = [pow(s, j, q) for j in range(6)] + [(-t * pow(s, j, q)) % q for j in range(6)]
    B = [[q] + [0] * 11] + [[(-w[j]) % q] + [int(k == j) for k in range(1, 12)] for j in range(1, 12)]
    B, swaps = lll(B)
    assert all(sum(a * b for a, b in zip(v, w)) % q == 0 for v in B)
    mu, nn = gs_exact(B); assert math.prod(nn) == q * q
    return q, s, t, w, B, swaps


def specialq_ball_stats(bits, side, seed, points=16, radii=(16, 96)):
    """exact (cf, cg) offsets of ball points of an actual special-q lattice:
    cf = log2 N_f - 72 c, cg = log2 N_g - 12 c - 6 log2|u|, c = (log2 q + S)/12."""
    import descent16 as DS
    rng = random.Random(seed)
    q, s, t, w, B, swaps = specialq_lattice(bits, side, rng)
    rows = []
    lq = math.log2(q)
    for Rr in radii:
        S = math.log2(DS.CELLS[Rr]); c = (lq + S) / 12
        cf = []; cg = []
        for _ in range(points):
            k = ball_point(rng, Rr)
            v = [sum(k[i] * B[i][j] for i in range(12)) for j in range(12)]
            a, b = C.norms(v[:6], v[6:]); assert (a if side == 'f' else b) % q == 0
            cf.append(math.log2(a) - 72 * c); cg.append(math.log2(b) - 12 * c - LOG2U6)
        rows.append(dict(R=Rr, S=S, c=c, points=points, cf_mean=sum(cf) / points, cf_max=max(cf), cg_mean=sum(cg) / points, cg_max=max(cg)))
    return dict(bits=bits, side=side, seed=seed, q=q, s=s, t=t, swaps=swaps, stats=rows, basis=B)


def small_units(bound=2):
    """all units of O_h with coefficients in [-bound, bound]^6 (|norm| = 1), without +-1."""
    import itertools
    out = []
    for coeffs in itertools.product(range(-bound, bound + 1), repeat=6):
        if any(coeffs) and C.norm(list(coeffs)) == 1 and list(coeffs) not in ([1, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0, 0]):
            out.append(list(coeffs))
    return out


def unit_mates(B, seed, Rr, points, units):
    """for ball points c of the lattice with basis B: the number of units eps in
    `units` whose mate eps*(a - b x) (exact, in the lattice) lies in the same ball;
    returns the per-point counts and E[1/(1+m)] (distinct/raw under orbit counting)."""
    from fractions import Fraction as Fr
    n = 12
    a = [[Fr(x) for x in row] + [Fr(int(i == j)) for j in range(n)] for i, row in enumerate(B)]
    for k in range(n):
        piv = next(i for i in range(k, n) if a[i][k]); a[k], a[piv] = a[piv], a[k]
        p = a[k][k]; a[k] = [x / p for x in a[k]]
        for i in range(n):
            if i != k and a[i][k]:
                f = a[i][k]; a[i] = [x - f * y for x, y in zip(a[i], a[k])]
    Binv = [row[n:] for row in a]
    rng = random.Random(seed); counts = []
    for _ in range(points):
        c = ball_point(rng, Rr)
        v = [sum(c[i] * B[i][j] for i in range(12)) for j in range(12)]
        m = 0
        for eps in units:
            vm = C.mul(eps, v[:6]) + C.mul(eps, v[6:])
            cm = [sum(vm[i] * Binv[i][j] for i in range(12)) for j in range(12)]
            assert all(x.denominator == 1 for x in cm)
            if sum(int(x) * int(x) for x in cm) <= Rr * Rr: m += 1
        counts.append(m)
    return dict(points=points, units=len(units), counts=counts, mean=sum(counts) / points,
                distinct_over_raw=sum(1 / (1 + m) for m in counts) / points)
