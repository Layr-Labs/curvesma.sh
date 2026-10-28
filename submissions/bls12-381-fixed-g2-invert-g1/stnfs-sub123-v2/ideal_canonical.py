""
from fractions import Fraction as F
from collections import Counter
from math import gcd, isqrt, prod, lcm
from common12 import H, mul, norm, integer_interval
def determinant(a):
    a = [list(row) for row in a]
    n, previous, sign = len(a), 1, 1
    for k in range(n-1):
        j = next((j for j in range(k, n) if a[j][k]), None)
        if j is None:
            return 0
        if j != k:
            a[k], a[j] = a[j], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1, n):
            for j in range(k+1, n):
                value = a[i][j]*pivot-a[i][k]*a[k][j]
                if value % previous:
                    raise AssertionError('inexact Bareiss division')
                a[i][j] = value//previous
        previous = pivot
    return sign*a[-1][-1] if n else 1
def bezout(a, b, counts):
    x, y, u, v = 1, 0, 0, 1
    while b:
        q, r = divmod(a, b)
        a, b = b, r
        x, u, y, v = u, x-q*u, v, y-q*v
        counts['euclidean_divisions'] += 1
    if a < 0:
        return -a, -x, -y
    return a, x, y
def row_hnf(rows, determinant_multiple, counts=None):
    ""
    counts = Counter() if counts is None else counts
    d = determinant_multiple
    n = len(rows[0])
    if type(d) is not int or d < 1 or any(len(r) != n for r in rows):
        raise ValueError('positive index multiple and rectangular generators')
    pending = [[x % d for x in row] for row in rows]
    hnf = []
    for col in range(n):
        pending.append([d if j == col else 0 for j in range(n)])
        k = next(j for j, row in enumerate(pending) if row[col])
        pivot = pending.pop(k)
        for j, row in enumerate(pending):
            if row[col]:
                a, b = pivot[col], row[col]
                g, u, v = bezout(a, b, counts)
                pnew = [u*x+v*y for x, y in zip(pivot, row)]
                rnew = [(b//g)*x-(a//g)*y for x, y in zip(pivot, row)]
                assert pnew[col] == g and rnew[col] == 0
                pivot, pending[j] = pnew, rnew
                assert all(abs(x) <= 2*d*d for x in pivot+pending[j])
                counts['max_intermediate_bits'] = max(counts['max_intermediate_bits'],
                                                      max(abs(x).bit_length() for x in pivot+pending[j]))
                for k in range(col+1, n):
                    pivot[k] %= d
                    pending[j][k] %= d
                counts['gcd_row_pairs'] += 1
        if pivot[col] < 0:
            pivot = [-x for x in pivot]
        assert pivot[col] > 0
        hnf.append(pivot)
        pending = [row for row in pending if any(row)]
    assert not pending
    for col in range(n):
        for i in range(col):
            q = hnf[i][col]//hnf[col][col]
            hnf[i] = [x-q*y for x, y in zip(hnf[i], hnf[col])]
            assert all(abs(x) <= 2*d*d for x in hnf[i])
            counts['max_intermediate_bits'] = max(counts['max_intermediate_bits'],
                                                  max(abs(x).bit_length() for x in hnf[i]))
            for j in range(col+1, n):
                hnf[i][j] %= d
            counts['residue_row_reductions'] += 1
    return tuple(tuple(row) for row in hnf)
def inner(a, b, metric):
    return sum(a[i]*metric[i][j]*b[j] for i in range(len(a)) for j in range(len(b)))
def lll(basis, metric, counts, max_steps=10000):
    ""
    b = [list(row) for row in basis]
    n = len(b)
    denominator = lcm(*(F(x).denominator for row in metric for x in row))
    metric = [[int(F(x)*denominator) for x in row] for row in metric]
    def gs():
        gram = [[inner(u, v, metric) for v in b] for u in b]
        mu, lengths = [], []
        for i in range(n):
            row = []
            for j in range(i):
                row.append((F(gram[i][j])-sum(row[k]*mu[j][k]*lengths[k]
                                             for k in range(j)))/lengths[j])
            length = F(gram[i][i])-sum(row[k]**2*lengths[k] for k in range(i))
            if length <= 0:
                raise ValueError('dependent basis/indefinite metric')
            mu.append(row); lengths.append(length)
        counts['gram_schmidt_recomputations'] += 1
        return mu, lengths
    k = 1
    while k < n:
        counts['lll_steps'] += 1
        if counts['lll_steps'] > max_steps:
            raise RuntimeError('LLL work budget exhausted')
        mu, lengths = gs()
        for j in range(k-1, -1, -1):
            z = (mu[k][j]+F(1, 2)).numerator//(mu[k][j]+F(1, 2)).denominator
            if z:
                b[k] = [x-z*y for x, y in zip(b[k], b[j])]
                counts['lll_size_reductions'] += 1
                mu, lengths = gs()
        if lengths[k] >= (F(3, 4)-mu[k][k-1]**2)*lengths[k-1]:
            k += 1
        else:
            b[k], b[k-1] = b[k-1], b[k]
            counts['lll_swaps'] += 1
            k = max(1, k-1)
    return b
def ellipsoid_vectors(metric, radius2, counts, max_nodes=100000):
    ""
    q = [[F(x) for x in row] for row in metric]
    n = len(q)
    projections = [q]
    while len(q) > 1:
        if q[0][0] <= 0:
            raise ValueError('nonpositive projection pivot')
        q = [[q[i][j]-q[i][0]*q[0][j]/q[0][0]
              for j in range(1, len(q))] for i in range(1, len(q))]
        projections.append(q)
    if q[0][0] <= 0:
        raise ValueError('nonpositive projection pivot')
    x = [0]*n
    def walk(i):
        counts['enumeration_nodes'] += 1
        if counts['enumeration_nodes'] > max_nodes:
            raise RuntimeError('orbit enumeration work budget exhausted')
        q = projections[i]
        suffix = x[i+1:]
        t = sum(q[0][j+1]*v for j, v in enumerate(suffix))
        c = sum(q[j+1][k+1]*v*w for j, v in enumerate(suffix)
                for k, w in enumerate(suffix))-radius2
        for value in integer_interval(q[0][0], t, c):
            x[i] = value
            if i:
                yield from walk(i-1)
            else:
                counts['enumerated_vectors'] += 1
                yield tuple(x)
    yield from walk(n-1)
def canonical_a(a, metric, max_nodes=100000):
    ""
    if len(a) != 6 or a[-1] <= 0:
        raise ValueError('eligible a requires six coefficients and a5>0')
    counts = Counter()
    q = [[F(x) for x in row] for row in metric]
    qa, na = inner(a, a, q), norm(a)
    if not na:
        raise ValueError('zero a')
    basis = lll([mul(a, [0]*i+[1]) for i in range(6)], q, counts)
    assert abs(determinant(basis)) == na
    transformed = [[inner(u, v, q) for v in basis] for u in basis]
    for x in ellipsoid_vectors(transformed, qa, counts, max_nodes):
        v = tuple(sum(x[i]*basis[i][j] for i in range(6)) for j in range(6))
        if v[-1] <= 0 or v == tuple(a):
            continue
        counts['associate_norm_evaluations'] += 1
        if norm(v) == na:
            qv = inner(v, v, q)
            if (qv, v) < (qa, tuple(a)):
                return False, dict(counts), v
    return True, dict(counts), None
def field_certificate(p):
    ""
    n = 6
    derivative = [i*H[i] for i in range(1, len(H))]
    sylvester = []
    for i in range(5):
        sylvester.append([0]*i+list(reversed(H))+[0]*(4-i))
    for i in range(6):
        sylvester.append([0]*i+list(reversed(derivative))+[0]*(5-i))
    discriminant = (-1)**15*determinant(sylvester)
    d = abs(discriminant)
    factors = []
    for q in range(2, isqrt(d)+1):
        if d % q == 0:
            e = 0
            while d % q == 0:
                d //= q; e += 1
            factors.append((q, e))
    if d > 1:
        factors.append((d, 1))
    squarefree = all(e == 1 for _, e in factors)
    def trim(a):
        a = [v % p for v in a]
        while a and a[-1] == 0:
            a.pop()
        return a
    def remainder(a, b):
        a, b = trim(a), trim(b)
        inv = pow(b[-1], -1, p)
        while len(a) >= len(b):
            z, shift = a[-1]*inv % p, len(a)-len(b)
            for j, v in enumerate(b):
                a[j+shift] = (a[j+shift]-z*v) % p
            a = trim(a)
        return a
    def poly_gcd(a, b):
        while b:
            a, b = b, remainder(a, b)
        return trim(a)
    def power(a, e):
        out = [1]
        while e:
            if e & 1:
                out = [v % p for v in mul(out, a)]
            a = [v % p for v in mul(a, a)]
            e >>= 1
        return trim(out)
    x, tests = [0, 1], {}
    for k in range(1, 7):
        x = power(x, p)
        diff = trim((x+[0]*6)[:1]+[(x+[0]*6)[1]-1]+(x+[0]*6)[2:6])
        if k in (2, 3):
            tests[k] = len(poly_gcd(H, diff))-1
    irreducible = not diff and tests == {2: 0, 3: 0}
    return dict(discriminant=discriminant, discriminant_factors=factors,
                squarefree=squarefree, irreducible_mod_p=irreducible,
                maximal_power_basis=squarefree and irreducible,
                rabin_gcd_degrees=tests)
