"""Murphy alpha of f and g over K_h for odd primes: Kummer-Dedekind root counts at regular primes, measured laws at 3, 5, 7.
alpha_p = log2(p) (1/(p-1) - E[v_p(N)]), E over primitive pairs; a simple root over a prime ideal of norm q gives q/(q^2-1)."""
import math
from random import Random
import common12 as C
FZ = [1, 6, 24, 36, 36]                      # F(z) = 36z^4 + 36z^3 + 24z^2 + 6z + 1, low degree first; f = w^4 F(x^2/w)
SING = (2, 3, 5, 7)                          # 2, 3 | 36 and disc F; 5 = N(w); 7 | disc h; disc h = 2^4 5^3 7^2, disc F = -2^8 3^5... (checked in _regular)
def _trim(a):
    while a and a[-1] == 0: a = a[:-1]
    return a
def pmulm(a, b, m, p):
    if not a or not b: return []
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): c[i + j] = (c[i + j] + x * y) % p
    return pmod(c, m, p)
def pmod(a, m, p):
    a = _trim([x % p for x in a]); dm = len(m) - 1
    if len(a) <= dm: return a
    inv = pow(m[-1], -1, p); a = list(a)
    for i in range(len(a) - 1, dm - 1, -1):
        if a[i]:
            f = a[i] * inv % p
            for j in range(dm + 1): a[i - dm + j] = (a[i - dm + j] - f * m[j]) % p
    return _trim(a[:dm])
def pgcd(a, b, p):
    a, b = _trim([x % p for x in a]), _trim([x % p for x in b])
    while b: a, b = b, pmod(a, b, p)
    if a: i = pow(a[-1], -1, p); a = [x * i % p for x in a]
    return a
def ppow(a, e, m, p):
    r = [1]
    while e:
        if e & 1: r = pmulm(r, a, m, p)
        a = pmulm(a, a, m, p); e >>= 1
    return r
def psub(a, b, p):
    n = max(len(a), len(b)); return _trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(n)])
def pdiv(a, b, p):
    a = list(a); q = [0] * max(1, len(a) - len(b) + 1); inv = pow(b[-1], -1, p)
    for i in range(len(a) - len(b), -1, -1):
        f = a[i + len(b) - 1] * inv % p; q[i] = f
        for j in range(len(b)): a[i + j] = (a[i + j] - f * b[j]) % p
    assert not _trim(a); return _trim(q)
def factor(m, p, rng):
    """Irreducible factors of a squarefree monic polynomial over F_p (distinct-degree, then equal-degree splitting)."""
    out = []; m = _trim([x % p for x in m]); i = pow(m[-1], -1, p); m = [x * i % p for x in m]; x = [0, 1]; xp = x; d = 0
    while len(m) - 1 >= 2 * (d + 1):
        d += 1; xp = ppow(xp, p, m, p); g = pgcd(psub(xp, x, p), m, p)
        if len(g) > 1:
            out += _eds(g, d, p, rng); m = pdiv(m, g, p); xp = pmod(xp, m, p)
    if len(m) > 1: out.append(m)
    return out
def _eds(g, d, p, rng):
    if len(g) - 1 == d: return [g]
    while True:
        t = [rng.randrange(p) for _ in range(len(g) - 1)]; s = psub(ppow(t, (p ** d - 1) // 2, g, p), [1], p); c = pgcd(s, g, p)
        if 1 < len(c) < len(g): return _eds(c, d, p, rng) + _eds(pdiv(g, c, p), d, p, rng)
def _regular(p):
    """p odd, h and F squarefree mod p with nonzero leading coefficients, p prime to N(w)."""
    if p in SING: return False
    der = lambda a: [i * c for i, c in enumerate(a)][1:]
    return len(pgcd(C.H, der(C.H), p)) == 1 and len(pgcd(FZ, der(FZ), p)) == 1
def local(p, rng=None):
    """(E[v_p N_f], E[v_p N_g]) at a regular prime p; g = x^2 - u w (a double root over a prime dividing u: 1/(q+1))."""
    rng = rng or Random(p); ef = eg = 0.
    for hi in factor(C.H, p, rng):
        d = len(hi) - 1; q = p ** d; w = pmod(C.W, hi, p); uw = pmod([C.U * x for x in C.W], hi, p)
        # g
        if not uw: eg += d / (q + 1.)
        elif ppow(uw, (q - 1) // 2, hi, p) == [1]: eg += d * 2 * q / (q * q - 1.)
        # f: roots z of F in F_q, x^2 = w z
        for fz in factor(FZ, p, rng):
            e = len(fz) - 1
            if d % e: continue
            for z in _roots_in(fz, hi, p, q, rng):
                s = pmulm(w, z, hi, p)
                if s and ppow(s, (q - 1) // 2, hi, p) == [1]: ef += d * 2 * q / (q * q - 1.)
    return ef, eg
def _roots_in(fz, hi, p, q, rng):
    """Roots in F_q = F_p[y]/hi of an irreducible fz over F_p whose degree divides [F_q : F_p]."""
    e = len(fz) - 1
    if e == 1: return [[(-fz[0] * pow(fz[1], -1, p)) % p]]
    # polynomials in z over F_q: lists of F_q elements (lists mod hi)
    K = lambda a, b: pmulm(a, b, hi, p)
    def zmulmod(A, B, M):
        c = [[] for _ in range(len(A) + len(B) - 1)]
        for i, x in enumerate(A):
            for j, y in enumerate(B): c[i + j] = psub(c[i + j], [(-t) % p for t in K(x, y)], p)
        return zmod(c, M)
    def zmod(A, M):
        A = [list(x) for x in A]
        while A and not A[-1]: A.pop()
        dm = len(M) - 1
        while len(A) > dm:
            f = A[-1]
            for j in range(dm): A[len(A) - 1 - dm + j] = psub(A[len(A) - 1 - dm + j], K(f, M[j]), p)
            A.pop()
            while A and not A[-1]: A.pop()
        return A
    def zgcd(A, B):
        while B:
            lead = B[-1]; inv = ppow(lead, q - 2, hi, p); B = [K(x, inv) for x in B]; A, B = B, zmod(A, B)
        return A
    M = [[c % p] if c % p else [] for c in fz]; inv = pow(fz[-1], -1, p); M = [[x * inv % p for x in c] for c in M]; todo = [M]; roots = []
    while todo:
        G = todo.pop()
        if len(G) == 2: roots.append([(-x) % p for x in G[0]]); continue
        while True:
            T = [[rng.randrange(p) for _ in range(len(hi) - 1)] for _ in range(len(G) - 1)]; R = [[1]]; e2 = (q - 1) // 2; Bq = T
            while e2:
                if e2 & 1: R = zmulmod(R, Bq, G)
                Bq = zmulmod(Bq, Bq, G); e2 >>= 1
            R = [list(x) for x in R] + [[]] * (1 if not R else 0); R[0] = psub(R[0], [1], p)
            while R and not R[-1]: R.pop()
            c = zgcd(G, R)
            if 1 < len(c) < len(G):
                # exact division G / c
                A = [list(x) for x in G]; Q = [[] for _ in range(len(G) - len(c) + 1)]
                for i in range(len(G) - len(c), -1, -1):
                    f = A[i + len(c) - 1]; Q[i] = f
                    for j in range(len(c)): A[i + j] = psub(A[i + j], K(f, c[j]), p)
                todo += [c, Q]; break
    return roots
def alpha(limit, sing, start=(3, 0., 0.)):
    """sing: {p: (E[v_p N_f], E[v_p N_g])} measured for the odd singular primes; the prime 2 is priced exactly per cell (model16.a2).
    start = (first prime, partial sums) continues a stored partial sum."""
    from threshold_model import primes_upto
    p0, af, ag = start
    for p in primes_upto(limit)[1:]:
        if p < p0: continue
        ef, eg = sing[p] if p in sing else local(p); assert p in sing or _regular(p); t = math.log2(p); af += t * (1 / (p - 1.) - ef); ag += t * (1 / (p - 1.) - eg)
    return af, ag
if __name__ == '__main__':
    import sys, time
    t = time.time(); irr = [p for p in __import__('threshold_model').primes_upto(3000) if p not in SING and not _regular(p)]; print('irregular odd primes < 3000 outside SING:', irr)
    for p in (11, 13, 17, 19, 23, 29, 31): print(p, local(p))
    print(time.time() - t)
