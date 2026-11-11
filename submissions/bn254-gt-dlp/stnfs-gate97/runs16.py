"""Canonical a-runs along e_0 = 1: the maximal t-interval around 0 with a + t canonical and Q_a(a + t) <= lim."""
import cmath, math
import common12 as C
from fractions import Fraction as F
K_EXP, LOG_RADIUS = 40, 9.
def _roots():
    c = [complex(x) for x in C.H]; z = [cmath.exp(2j * math.pi * k / 6 + .4j) * 1.1 for k in range(6)]; ev = lambda x: sum(ck * x ** k for k, ck in enumerate(c))
    for _ in range(500): z = [zi - ev(zi) / math.prod(zi - zj for j, zj in enumerate(z) if j != i) for i, zi in enumerate(z)]
    return z
def _inv(u):
    rows = [C.mul([0] * j + [1], u) for j in range(6)]; M = [[F(rows[j][i]) for j in range(6)] + [F(int(i == 0))] for i in range(6)]
    for k in range(6):
        p = next(i for i in range(k, 6) if M[i][k]); M[k], M[p] = M[p], M[k]; pv = M[k][k]; M[k] = [x / pv for x in M[k]]
        for i in range(6):
            if i != k and M[i][k]: f = M[i][k]; M[i] = [x - f * y for x, y in zip(M[i], M[k])]
    inv = [int(M[i][6]) for i in range(6)]; assert C.mul(inv, u) == C.ONE; return inv
def _emb(u, z): return sum(c * z ** k for k, c in enumerate(u))
def unit_list():
    """All units +- u1^i u2^j u3^k (up to sign) whose logarithmic embeddings are at most LOG_RADIUS in absolute value."""
    r = _roots(); lg = [[math.log(abs(_emb(u, x))) for x in r] for u in C.UNITS]; pw = []
    for u in C.UNITS:
        d = {0: C.ONE}; ui = _inv(u)
        for e in range(1, K_EXP + 1): d[e] = C.mul(d[e - 1], u); d[-e] = C.mul(d[-e + 1], ui)
        pw.append(d)
    out = []; R = range(-K_EXP, K_EXP + 1)
    for i in R:
        for j in R:
            for k in R:
                if max(abs(i * a + j * b + k * c) for a, b, c in zip(*lg)) <= LOG_RADIUS: assert max(abs(i), abs(j), abs(k)) < K_EXP; out.append(C.mul(C.mul(pw[0][i], pw[1][j]), pw[2][k]))
    return out
_QU = {}
def run(a, units, q, lim):
    def form(x, y): return sum(x[i] * q[i][j] * y[j] for i in range(6) for j in range(6))
    E = [[int(i == j) for i in range(6)] for j in range(6)]; qa0, bad0, qd0 = form(a, a), form(a, E[0]), q[0][0]
    if id(units) not in _QU: _QU[id(units)] = [form(u, u) for u in units]
    coef = []
    for t in a, C.sigma(a):
        cols = [C.mul(t, E[j]) for j in range(6)]; G = [[form(cols[i], cols[j]) for j in range(6)] for i in range(6)]; Bd = [[form(cols[i], E[j]) for j in range(6)] for i in range(6)]
        for u, quu in zip(units, _QU[id(units)]):
            if t is a and u in (C.ONE, [-x for x in C.ONE]): continue
            qu = sum(u[i] * G[i][j] * u[j] for i in range(6) for j in range(6)); bu = sum(u[i] * Bd[i][j] * u[j] for i in range(6) for j in range(6)); coef.append((qu - qa0, 2 * (bu - bad0), quu - qd0))
    if min(c[0] for c in coef) < 0: return 0, 0
    def ok(t): x = [a[0] + t] + a[1:]; return form(x, x) <= lim and all(c0 + c1 * t + c2 * t * t >= 0 for c0, c1, c2 in coef)
    def first(sgn):
        best = math.inf
        for al, be, ga in coef:
            B_ = be * sgn
            if al < 0: return .0
            if ga == 0:
                if B_ < 0: best = min(best, -al / B_)
                continue
            d = B_ * B_ - 4 * ga * al
            if d < 0: continue
            s = math.sqrt(d)
            for rt in (-B_ - s) / (2 * ga), (-B_ + s) / (2 * ga):
                if rt > 0: best = min(best, rt)
        d = bad0 * bad0 - qd0 * (qa0 - lim); s = math.sqrt(max(d, 0)); return min(best, (-bad0 + s if sgn > 0 else bad0 + s) / qd0)
    ends = []
    for sgn in 1, -1:
        t = max(0, math.floor(first(sgn) + 1e-9))
        while t > 0 and not ok(sgn * t): t -= 1
        while ok(sgn * (t + 1)): t += 1
        ends.append(sgn * t)
    return ends[1], ends[0]
def lrs(a, t_lo, t_hi, q, den, A2, lam=1, kap=1):
    """log2 of sum over the run's a-values of (lam - kap Q_a(a + t)/A^2)^3 (exact integer sum)."""
    lam, kap = F(lam), F(kap); d = lam.denominator * kap.denominator; c = int(lam * d) * den * A2; k = int(kap * d)
    qa0 = sum(a[i] * q[i][j] * a[j] for i in range(6) for j in range(6)); b0 = sum(a[i] * q[i][0] for i in range(6)); q00 = q[0][0]
    s = sum((c - k * (qa0 + 2 * b0 * t + q00 * t * t)) ** 3 for t in range(t_lo, t_hi + 1)); return math.log2(s) - 3 * math.log2(d * den * A2)
