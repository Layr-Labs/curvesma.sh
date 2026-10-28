import math, random, functools, os, sys
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   

import retain16 as RT          
import descent16 as D          
from common12 import H, U, TERMS

GRID = RT.GRID                
BAD = RT.BAD                  
G_TERMS = RT.G_TERMS




OLD_HIST = {
    2:  dict(K=104, n=262144, f={0: 262144},
             g={0: 197164, 4: 37099, 8: 20852, 12: 4497, 16: 1900, 20: 435,
                24: 143, 28: 34, 32: 14, 34: 2, 36: 3, 40: 1}),
    3:  dict(K=66, n=262144, f={0: 262144}, g={0: 262144}),
    19: dict(K=25, n=262144, f={0: 262133, 4: 11}, g={0: 262138, 4: 6}),
    23: dict(K=23, n=262144, f={0: 262144},
             g={0: 220705, 1: 37963, 2: 3210, 3: 250, 4: 15, 5: 1}),
    67: dict(K=18, n=262144, f={0: 262144}, g={0: 262144}),
}







import json as _json
STRUCT_TABLE = _json.loads('{"101":[[1,0,2],[2,4,0],[3,0,0]],"103":[[1,0,0],[5,2,2]],"107":[[1,0,0],[5,0,2]],"109":[[2,6,2],[4,6,2]],"11":[[6,4,0]],"113":[[2,0,0],[4,0,0]],"127":[[1,4,2],[1,4,2],[4,0,0]],"13":[[1,0,0],[5,0,0]],"131":[[6,0,0]],"137":[[1,0,0],[1,0,2],[1,0,0],[3,0,2]],"139":[[1,2,0],[2,6,2],[3,0,2]],"149":[[1,0,2],[5,0,2]],"151":[[2,2,2],[4,0,0]],"157":[[1,0,2],[1,4,0],[4,0,0]],"163":[[1,4,2],[5,0,0]],"167":[[1,0,0],[5,0,2]],"17":[[6,12,2]],"173":[[1,0,0],[1,0,0],[4,4,2]],"179":[[6,0,0]],"181":[[1,0,0],[5,0,0]],"191":[[1,0,2],[5,0,0]],"193":[[2,6,2],[4,6,2]],"197":[[1,0,2],[5,0,2]],"199":[[6,4,0]],"211":[[6,4,0]],"223":[[1,0,2],[2,0,0],[3,0,2]],"227":[[1,0,0],[2,0,2],[3,0,2]],"229":[[1,2,2],[5,2,2]],"233":[[1,0,2],[1,0,2],[4,12,2]],"239":[[6,0,0]],"241":[[1,0,0],[1,0,0],[4,6,2]],"251":[[2,0,2],[4,8,0]],"257":[[6,0,2]],"263":[[1,0,2],[5,0,0]],"269":[[1,0,2],[2,0,2],[3,0,2]],"271":[[2,6,2],[2,6,2],[2,0,0]],"277":[[1,4,2],[1,4,2],[4,6,2]],"281":[[1,0,2],[5,0,2]],"283":[[2,0,2],[4,0,0]],"29":[[1,0,2],[5,0,2]],"293":[[1,0,0],[5,0,0]],"307":[[1,2,0],[5,0,2]],"31":[[1,2,2],[2,4,0],[3,2,2]],"311":[[1,0,0],[1,0,2],[4,4,2]],"313":[[1,0,2],[5,0,2]],"317":[[1,0,0],[1,0,2],[2,4,2],[2,0,0]],"331":[[6,0,0]],"337":[[6,4,2]],"347":[[2,4,0],[4,4,2]],"349":[[1,2,0],[1,2,2],[2,4,0],[2,8,2]],"353":[[1,0,0],[1,0,0],[4,12,2]],"359":[[1,0,0],[1,0,2],[4,0,2]],"367":[[1,2,2],[2,2,2],[3,6,0]],"37":[[3,6,2],[3,6,2]],"373":[[3,2,2],[3,2,2]],"379":[[1,2,2],[1,2,2],[1,2,2],[3,2,0]],"383":[[1,0,2],[1,0,0],[4,0,2]],"389":[[1,0,0],[5,0,0]],"397":[[1,2,0],[1,2,2],[1,2,0],[3,2,2]],"401":[[6,12,2]],"409":[[6,8,2]],"41":[[1,0,2],[5,0,2]],"419":[[2,4,2],[4,0,0]],"421":[[1,0,0],[1,0,0],[4,6,2]],"43":[[1,0,2],[5,2,0]],"431":[[1,0,0],[1,0,0],[2,0,0],[2,4,2]],"433":[[3,2,0],[3,2,0]],"439":[[1,2,2],[5,0,0]],"443":[[3,0,2],[3,0,0]],"449":[[1,0,0],[2,4,2],[3,0,0]],"457":[[1,4,0],[2,12,2],[3,4,0]],"461":[[3,0,2],[3,0,2]],"463":[[1,4,2],[2,0,0],[3,4,2]],"467":[[1,0,0],[2,8,0],[3,0,0]],"47":[[1,0,0],[5,0,2]],"479":[[1,0,2],[2,4,2],[3,0,0]],"487":[[6,0,0]],"491":[[1,0,2],[1,0,2],[1,0,2],[3,0,0]],"499":[[6,4,0]],"5":[[6,0,2]],"503":[[2,4,2],[4,0,0]],"509":[[1,0,2],[5,0,2]],"53":[[1,0,2],[1,0,2],[1,0,2],[1,0,0],[2,0,0]],"59":[[3,0,0],[3,0,2]],"61":[[1,4,2],[5,4,2]],"7":[[6,4,0]],"71":[[1,0,0],[1,0,0],[1,0,2],[3,0,0]],"73":[[1,0,2],[1,4,0],[2,4,2],[2,8,0]],"79":[[1,0,2],[5,0,0]],"83":[[1,0,2],[5,0,0]],"89":[[2,0,2],[4,4,2]],"97":[[1,2,0],[5,2,0]]}')

def structure_of17(p):
    t = STRUCT_TABLE.get(str(p))
    if t is not None:
        return [dict(d1=d1, f=[(1, nf)], g=[(1, ng)]) for d1, nf, ng in t]
    return RT.structure(p)

def verify_struct_table():
    primes = RT.primes(2, 512)
    checked = 0
    for i, p in enumerate(primes):
        if p in BAD:
            assert str(p) not in STRUCT_TABLE
            continue
        if i % 7 == 0:
            st = RT.structure(p)
            got = None if st is None else [[b['d1'], b['f'][0][1], b['g'][0][1]] for b in st]
            assert got == STRUCT_TABLE[str(p)], (p, got, STRUCT_TABLE[str(p)])
            checked += 1
    return checked, len(primes)

def _trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a

def pmod(a, m, p):
    a = _trim([x % p for x in a]); dm = len(m) - 1; inv = pow(m[-1], -1, p)
    while len(a) - 1 >= dm:
        c = a[-1] * inv % p; sh = len(a) - 1 - dm
        for i in range(len(m)):
            a[sh + i] = (a[sh + i] - c * m[i]) % p
        _trim(a)
    return a

def pmul(a, b, p):
    c = [0] * (len(a) + len(b) - 1) if a and b else []
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] = (c[i + j] + x * y) % p
    return _trim(c)

def padd(a, b, p):
    return _trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p
                  for i in range(max(len(a), len(b)))])

def pscal(a, c, p):
    return _trim([x * c % p for x in a])

def pdiv(a, b, p):
    a = _trim([x % p for x in a]); b = _trim([x % p for x in b])
    q = [0] * max(1, len(a) - len(b) + 1); inv = pow(b[-1], -1, p)
    while a and len(a) >= len(b):
        c = a[-1] * inv % p; sh = len(a) - len(b); q[sh] = c
        for i in range(len(b)):
            a[sh + i] = (a[sh + i] - c * b[i]) % p
        _trim(a)
    return _trim(q), a

def pgcd(a, b, p):
    a = _trim([x % p for x in a]); b = _trim([x % p for x in b])
    while b:
        a, b = b, pmod(a, b, p)
    if a:
        inv = pow(a[-1], -1, p); a = [x * inv % p for x in a]
    return a

def pderiv(a, p):
    return _trim([(i * a[i]) % p for i in range(1, len(a))])

def factor_Fp(a, p, seed=1):
    rng = random.Random(seed)
    def split(a, d):
        if len(a) - 1 == d:
            return [a]
        for _ in range(10000):
            r = _trim([rng.randrange(p) for _ in range(len(a) - 1)])
            if len(r) < 2:
                continue
            if p == 2:
                t = r; s = r
                for _ in range(d - 1):
                    t = pmod(pmul(t, t, p), a, p); s = padd(s, t, p)
                g = pgcd(a, s, p)
            else:
                g = pgcd(a, padd(RT.ppowmod(r, (p ** d - 1) // 2, a, p), [p - 1], p), p)
            if 1 < len(g) - 1 < len(a) - 1:
                q, _ = pdiv(a, g, p); return split(g, d) + split(q, d)
        raise RuntimeError('splitting failed')
    out = []; rest = _trim([x % p for x in a])
    for r in range(p):
        if len(rest) > 1 and sum(c * pow(r, i, p) for i, c in enumerate(rest)) % p == 0:
            out.append([(-r) % p, 1]); rest, _ = pdiv(rest, [(-r) % p, 1], p)
    d = 2; x = [0, 1]
    xp = RT.ppowmod(x, p, rest, p) if len(rest) > 1 else x
    while len(rest) - 1 >= 2 * d:
        xp = RT.ppowmod(xp, p, rest, p)
        g = pgcd(rest, padd(xp, pscal(x, p - 1, p), p), p)
        if len(g) > 1:
            out += split(g, d); rest, _ = pdiv(rest, g, p)
            if len(rest) > 1:
                xp = pmod(xp, rest, p)
        d += 1
    if len(rest) > 1:
        out.append(rest)
    return out

def h_factor_data(p):
    Hp = [x % p for x in H]
    g = pgcd(Hp, pderiv(Hp, p), p)
    sqf = pdiv(Hp, g, p)[0] if len(g) > 1 else Hp
    return factor_Fp(sqf, p), (factor_Fp(g, p) if len(g) > 1 else [])


def _canon(co):
    co = list(co)
    while co and co[-1] == 0:
        co.pop()
    return co

def _e_sub(A, B, p):
    n = max(len(A), len(B))
    out = []
    for i in range(n):
        a = A[i] if i < len(A) else []
        b = B[i] if i < len(B) else []
        out.append(padd(a, pscal(b, p - 1, p), p) if b else a)
    return RT.Etrim(out)

def _esub1(R, p):
    R = [list(c) for c in R]
    if not R:
        return [[p - 1]]
    R[0] = padd(R[0], [p - 1], p)
    return RT.Etrim(R)

def _e_div(A, B, m, p):
    A = [list(c) if c else [] for c in A]
    B = [list(c) if c else [] for c in B]
    inv = RT.einv(B[-1], m, p)
    Q = [[] for _ in range(len(A) - len(B) + 1)]
    while True:
        A = RT.Etrim(A)
        if not A or len(A) - 1 < len(B) - 1:
            break
        c = RT.emul(A[-1], inv, m, p)
        sh = len(A) - 1 - (len(B) - 1)
        Q[sh] = c
        for i in range(len(B)):
            A[sh + i] = padd(A[sh + i], pscal(RT.emul(c, B[i], m, p), p - 1, p), p)
    return RT.Etrim(Q)

def fq_roots(phi, hi, p, seed=1):
    q = p ** (len(hi) - 1)
    phi = RT.Etrim([_canon(list(c)) for c in phi])
    if not phi:
        return []
    XP = RT.Epowmod([[], [1]], q, phi, hi, p)
    Xp = [[1] if i == 1 else [] for i in range(max(len(XP), 2))]
    g = RT.Etrim(RT.Egcd(phi, _e_sub(XP, Xp, p), hi, p))
    if len(g) <= 1:
        return []
    rng = random.Random(seed)
    def split(g):
        if len(g) - 1 == 1:
            root = RT.emul(g[0] if g[0] else [], RT.einv(g[1], hi, p), hi, p)
            return [pscal(root, p - 1, p) if root else []]
        for _ in range(10000):
            r = RT.Etrim([[_canon([rng.randrange(p) for _ in range(len(hi) - 1)]) for _ in range(len(g) - 1)][0] for _ in range(len(g) - 1)])
            r = RT.Etrim([_canon([rng.randrange(p) for _ in range(len(hi) - 1)]) for _ in range(len(g) - 1)])
            if len(r) < 2:
                continue
            if p == 2:
                rr = RT.Emod(RT.Emul(r, r, hi, p), g, hi, p)
                hh = _e_sub(rr, r, p)
            else:
                R = RT.Epowmod(r, (q - 1) // 2, g, hi, p)
                hh = _esub1(R, p)
            h = RT.Etrim(RT.Egcd(g, hh, hi, p))
            if 0 < len(h) - 1 < len(g) - 1:
                return split(h) + split(_e_div(g, h, hi, p))
        raise RuntimeError('E-splitting failed')
    return split(g)


def _poly_mul_mod(a, b, mod):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                c[i + j] += x * y
    return [x % mod for x in c]

def _solve_modp(rows, n, p):
    A = [r[:] for r in rows]
    piv = []; r = 0
    for c in range(n):
        pr = None
        for i in range(r, len(A)):
            if A[i][c] % p:
                pr = i; break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        inv = pow(A[r][c] % p, -1, p)
        A[r] = [x * inv % p for x in A[r]]
        for i in range(len(A)):
            if i != r and A[i][c] % p:
                f = A[i][c]
                A[i] = [(A[i][m] - f * A[r][m]) % p for m in range(n + 1)]
        piv.append(c); r += 1
    sol = [0] * n
    for i, c in enumerate(piv):
        sol[c] = A[i][n] % p
    for row in A:
        if all(x % p == 0 for x in row[:n]) and row[n] % p:
            raise RuntimeError('inconsistent lifting system')
    return sol

def lift_coprime(facs, p, M):
    hh_full = list(H)
    n = len(hh_full) - 1
    F = [[x % p for x in f] for f in facs]
    degs = [len(f) - 1 for f in F]
    k = len(F)
    for j in range(1, M):
        pj = p ** j
        modj = p ** (j + 1)
        prod = [1]
        for f in F:
            prod = _poly_mul_mod(prod, [x % modj for x in f], modj)
            if len(prod) - 1 > n:
                prod = prod[:n + 1]
        err = [((hh_full[i] if i < len(hh_full) else 0) - (prod[i] if i < len(prod) else 0)) // pj % p
               for i in range(n + 1)]
        if not any(err):
            continue
        nunk = sum(degs)
        G = []
        for i in range(k):
            g = [1]
            for l in range(k):
                if l != i:
                    g = _poly_mul_mod(g, [x % p for x in F[l]], p)
            G.append(g)
        rows = []
        for deg in range(n + 1):
            row = [0] * nunk
            for i in range(k):
                for t in range(degs[i]):
                    if 0 <= deg - t < len(G[i]):
                        row[sum(degs[:i]) + t] = (row[sum(degs[:i]) + t] + G[i][deg - t]) % p
            rows.append(row + [err[deg] % p])
        sol = _solve_modp(rows, nunk, p)
        for i in range(k):
            for t in range(degs[i]):
                if sol[sum(degs[:i]) + t]:
                    F[i][t] = (F[i][t] + pj * sol[sum(degs[:i]) + t]) % p ** M
    return [([x % p ** M for x in f] + [0] * (d + 1 - len(f)))[:d + 1] for f, d in zip(F, degs)]


class UnramRing:
    kind = 'unram'
    def __init__(self, p, hi, M, HM):
        self.p, self.d, self.M, self.hi = p, len(hi) - 1, M, list(hi)
        self.HM = [x % (p ** M) for x in HM]
        self.mod = p ** M
        self.q = p ** self.d
    def zero(self):
        return (0,) * self.d
    def one(self):
        return tuple(1 if i == 0 else 0 for i in range(self.d))
    def lift_y(self, asc):
        c = [x % self.mod for x in asc]
        hm = self.HM
        while len(c) - 1 >= self.d:
            k = c[-1] % self.mod
            sh = len(c) - 1 - self.d
            for i in range(len(hm)):
                c[sh + i] = (c[sh + i] - k * hm[i]) % self.mod
            c.pop()
        return tuple((c + [0] * self.d)[:self.d])
    def add(self, a, b):
        return tuple((x + y) % self.mod for x, y in zip(a, b))
    def sub(self, a, b):
        return tuple((x - y) % self.mod for x, y in zip(a, b))
    def mul(self, a, b):
        c = [0] * (2 * self.d - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    c[i + j] += x * y
        hm = self.HM
        for i in range(len(c) - 1, self.d - 1, -1):
            k = c[i] % self.mod
            if k:
                sh = i - self.d
                for j in range(len(hm)):
                    c[sh + j] -= k * hm[j]
            c[i] = 0
        return tuple(x % self.mod for x in c[:self.d])
    def scale(self, a, k):
        return tuple(x * k % self.mod for x in a)
    def div_u(self, a):
        return tuple(x // self.p for x in a)
    def val(self, a):
        def vp(x):
            if x == 0:
                return self.M
            v = 0
            while v < self.M and x % self.p == 0:
                x //= self.p; v += 1
            return v
        return min(vp(x) for x in a)
    def res(self, a):
        return tuple(x % self.p for x in a)
    def lift_res(self, r):
        return tuple((r[i] if i < len(r) else 0) % self.mod for i in range(self.d))
    def upow(self, j):
        return self.scale(self.one(), self.p ** j) if j else self.one()

class RamRing23:
    kind = 'ram'
    def __init__(self, rho, F1, m):
        self.p = 23; self.M = m; self.mod = 23 ** m; self.d = 1; self.q = 23
        self.rho = rho % self.mod
        fr = 0; fpr = 0
        for i, c in enumerate(F1):
            fr = (fr + c * pow(self.rho, i, self.mod)) % self.mod
            if i:
                fpr = (fpr + i * c * pow(self.rho, i - 1, self.mod)) % self.mod
        assert fr % 23 == 0 and fpr % 23 == 0
        self.b23 = (fr // 23) % self.mod
        self.a23 = (fpr // 23) % self.mod
        assert self.b23 % 23 != 0, 'b23 must be a unit (this is v_P(y-rho) = 1)'
        self.binv = pow(self.b23, -1, self.mod)
        self.abinv = self.a23 * self.binv % self.mod
        self.hi = [(-rho) % 23, 1]
    def zero(self):
        return (0, 0)
    def one(self):
        return (1 % self.mod, 0)
    def lift_y(self, asc):
        pk = [(1 % self.mod, 0)]
        for _ in range(1, max(len(asc), 1)):
            a_, b_ = pk[-1]
            pk.append(((a_ * self.rho - 23 * b_ * self.b23) % self.mod,
                       (a_ + b_ * self.rho - 23 * b_ * self.a23) % self.mod))
        al = 0; be = 0
        for k, c in enumerate(asc):
            if c:
                al = (al + c * pk[k][0]) % self.mod
                be = (be + c * pk[k][1]) % self.mod
        return (al, be)
    def add(self, a, b):
        return ((a[0] + b[0]) % self.mod, (a[1] + b[1]) % self.mod)
    def sub(self, a, b):
        return ((a[0] - b[0]) % self.mod, (a[1] - b[1]) % self.mod)
    def mul(self, a, b):
        t = a[0] * b[0]
        u = a[0] * b[1] + a[1] * b[0]
        v = a[1] * b[1]
        return ((t - 23 * self.b23 * v) % self.mod, (u - 23 * self.a23 * v) % self.mod)
    def scale(self, a, k):
        return (a[0] * k % self.mod, a[1] * k % self.mod)
    def div_u(self, a):
        al, be = a
        assert al % 23 == 0, 'div_pi needs v_P >= 1'
        return ((be - self.abinv * al) % self.mod, ((-al // 23) * self.binv) % self.mod)
    def val(self, a):
        def v23(x):
            v = 0
            while x % 23 == 0 and v < self.M:
                x //= 23; v += 1
            return v
        return min(2 * v23(a[0]), 2 * v23(a[1]) + 1)
    def res(self, a):
        return (a[0] % 23,)
    def lift_res(self, r):
        return ((r[0] if r else 0) % self.mod, 0)
    def upow(self, j):
        if j == 0:
            return self.one()
        pi = self.sub(self.lift_y([0, 1]), (self.rho % self.mod, 0))
        out = self.one()
        for _ in range(j):
            out = self.mul(out, pi)
        return out

def local_primes(p, M):
    sqf, rep = h_factor_data(p)
    out = []
    if p == 23:
        assert len(rep) == 1 and len(rep[0]) == 2, 'expected one repeated linear factor at 23'
        rho = (-rep[0][0]) * pow(rep[0][1], -1, 23) % 23
        others = [f for f in sqf if not (len(f) == 2 and (-f[0]) * pow(f[1], -1, 23) % 23 == rho)]
        ram = pmul(rep[0], rep[0], 23)
        F = lift_coprime([ram] + others, 23, M)
        out.append((1, RamRing23(rho, F[0], 2 * M + 8)))
        for hi, HM in zip(others, F[1:]):
            out.append((len(hi) - 1, UnramRing(23, hi, M, HM)))
        return out
    F = lift_coprime(sqf, p, M)
    for hi, HM in zip(sqf, F):
        out.append((len(hi) - 1, UnramRing(p, hi, M, HM)))
    return out


def psi_coefs(side, ring):
    co = [ring.zero()] * 13
    if side == 'g':
        co[0] = ring.lift_y([0, U])          
        co[2] = ring.one()
    else:
        for deg, c in TERMS:
            co[deg] = ring.lift_y(c)
    return co

def child_poly(ring, phi, xi):
    xil = ring.lift_res(xi)
    xpow = [ring.one()]
    for _ in range(12):
        xpow.append(ring.mul(xpow[-1], xil))
    out = []
    for j in range(13):
        S = ring.zero()
        for k in range(j, 13):
            c = math.comb(k, j)
            S = ring.add(S, ring.mul(ring.scale(phi[k], c), xpow[k - j]))
        out.append(ring.div_u(ring.mul(S, ring.upow(j))))
    return out

def one_var_law(ring, psi, Vc, seed=1):
    q = ring.q
    nu = [Fr(0)] * Vc
    tail = Fr(0)
    stats = dict(nodes=0, max_run=0, root_calls=0, max_mult=0)
    stack = [(list(psi), 0, Fr(1))]
    while stack:
        phi, L, w = stack.pop()
        stats['nodes'] += 1
        run = 0
        while L < Vc and all(ring.val(c) >= 1 for c in phi):
            phi = [ring.div_u(c) for c in phi]
            L += 1; run += 1
        stats['max_run'] = max(stats['max_run'], run)
        if L >= Vc:
            tail += w
            continue
        assert any(ring.val(c) < ring.M for c in phi), 'precision exhausted'
        ephi = [_canon(list(ring.res(phi[k]))) for k in range(13)]
        roots = fq_roots(ephi, ring.hi, ring.p, seed)
        stats['root_calls'] += 1
        stats['max_mult'] = max(stats['max_mult'], len(roots))
        m = len(roots)
        assert m <= q
        nu[L] += w * Fr(q - m, q)
        for xi in roots:
            if L + 1 >= Vc:
                tail += w * Fr(1, q)
            else:
                stack.append((child_poly(ring, phi, xi), L + 1, w * Fr(1, q)))
    return nu, tail, stats

@functools.lru_cache(None)
def side_law(p, side, Vc):
    law = {0: Fr(1)}
    details = []
    for f_P, ring in local_primes(p, Vc + 12):
        psi = psi_coefs(side, ring)
        nu, tail, st = one_var_law(ring, psi, Vc // f_P + 2)
        q = ring.q
        qp = Fr(1, q)
        pl = {}
        for v in range(len(nu)):
            if v == 0:
                continue
            m = nu[v] / (1 + qp)
            if m:
                pl[f_P * v] = pl.get(f_P * v, Fr(0)) + m
        tot = sum(pl.values(), Fr(0))
        pl[0] = 1 - tot
        new = {}
        for x, px in law.items():
            for y, py in pl.items():
                new[x + y] = new.get(x + y, Fr(0)) + px * py
        law = new
        details.append(dict(f=f_P, q=q, nu=nu, tail=tail, stats=st))
    
    law = {v: m for v, m in law.items() if v < Vc}
    tail = 1 - sum(law.values(), Fr(0))
    return law, tail, details


VC = {2: 130, 3: 20, 19: 40, 23: 45, 67: 25}

@functools.lru_cache(None)
def exact_side_component(p, side):
    return side_law(p, side, VC[p])[:2]

def check_reduction_units():
    out = {}
    for p in BAD:
        for f_P, ring in local_primes(p, 8):
            if f_P not in out:
                out[f_P] = []
            c0 = ring.lift_y([1, 1, -1, -1])
            out[f_P].append((p, ring.kind, ring.val(c0)))
    assert all(v == 0 for lst in out.values() for (_, _, v) in lst), 'c0 not a unit!'
    return out

def bad_component(p, side, up, law=None, tail=None, inflate=None):
    lp = math.log2(p)
    rnd = math.ceil if up else math.floor
    if law is None:
        law, tail = exact_side_component(p, side)
    comp = {0: float(law.get(0, 0.0))}
    for v, m in law.items():
        if v == 0 or m == 0:
            continue
        mm = float(m) * inflate if inflate else float(m)
        if mm < 1e-33:
            continue
        x = rnd(v * lp * GRID)
        comp[x] = comp.get(x, 0.0) + mm
    if tail and float(tail) > 0.0:
        
        
        
        x = rnd(VC[p] * lp * GRID)
        comp[x] = comp.get(x, 0.0) + float(tail) * (inflate or 1.0)
    return comp

def bad_component_old(p, side, up):
    d = OLD_HIST[p]
    lp = math.log2(p)
    rnd = math.ceil if up else math.floor
    comp = {}
    for v, c in d[side].items():
        v = int(v)
        comp[rnd(v * lp * GRID)] = comp.get(rnd(v * lp * GRID), 0.0) + c / d['n']
    return comp

@functools.lru_cache(None)
def omitted_pdf17(side, up, old=False, inflate=None):
    pdf = {0: 1.0}; ordinary = 1.0
    for p in RT.primes(2, 512):
        if p in BAD:
            comp = bad_component_old(p, side, up) if old else bad_component(p, side, up, inflate=inflate)
        else:
            st = structure_of17(p)
            if st is None:
                raise ValueError('unexpected non-squarefree prime %d' % p)
            if p >= 256:
                for b in st:
                    d1 = b['d1']; n1 = sum(n for dd, n in b[side] if dd == 1)
                    ordinary -= n1 * p ** (-d1) if d1 > 1 else 0.0
            comp = RT.good_prime_pdf(p, st, side, p >= 256, up)
        pdf = RT._convolve(pdf, comp)
        pdf = {x: v for x, v in pdf.items() if v > 1e-30}
    j0 = int(math.log2(512))
    for j in range(j0, 69):
        pdf = RT._convolve(pdf, RT.large_block_pdf(j, up))
        pdf = {x: v for x, v in pdf.items() if v > 1e-30}
        ordinary -= (36.0 if side == 'f' else 6.0) * 1.26 * 2.0 ** (j + 1) / ((j + 1) * math.log(2)) * 2.0 ** (-2 * j)
    return pdf, ordinary

@functools.lru_cache(None)
def loss17(side, m, n, iA, kA, kB, alpha, slack, step=64, old=False, inflate=None):
    s_units = int(math.floor((slack - 1.0 - m / 128.0) * GRID))
    pdf_up, _ = omitted_pdf17(side, True, old=old, inflate=inflate)
    base = D.P.LX
    def hi(x_bits):
        xg = math.ceil(x_bits / step) * step
        return RT.ss_upper2(max(m - xg, 1.0), base, n, iA, kA, kB, D.P.MESH, alpha) if m - xg > 1 else 1.0
    def lo(x_bits):
        xg = math.floor(x_bits / step) * step
        return D.ss_lower2(max(m - xg, 1.0), base, n, iA, kA, kB, D.P.MESH, alpha, D.P.TRUNC) if m - xg > 1 else 1.0
    num = sum(pr * hi(x / GRID) for x, pr in pdf_up.items() if x > s_units)
    den = sum(pr * lo(x / GRID) for x, pr in pdf_up.items() if x <= s_units)
    return num / den if den > 0 else 1.0

def retention17(nf, ng, n, iA, kA, kB, slack, old=False, inflate=None):
    lf = loss17('f', math.ceil(nf), n, iA, kA, kB, D.P.ALPHA_F, slack, old=old, inflate=inflate)
    lg = loss17('g', math.ceil(ng), n, iA, kA, kB, D.P.ALPHA_G, slack, old=old, inflate=inflate)
    _, ordf = omitted_pdf17('f', True, old=old)
    _, ordg = omitted_pdf17('g', True, old=old)
    return max(0.0, (1 - lf - lg)) * max(0.0, ordf + ordg - 1) * RT.base_primitive(), lf, lg

def cell_loss17(side, slack, old=False):
    pdf, _ = omitted_pdf17(side, True, old=old)
    return sum(pr for x, pr in pdf.items() if x > slack * GRID)


def node_norms():
    bins = D.bins_of(D.P)
    cfmax = max(b[1] for b in D.P.QBINS); cgmax = max(b[2] for b in D.P.QBINS)
    rows = []
    for b in bins:
        for side in ('f', 'g'):
            n, iA, kA, kB, R, thin = D.PLAN16[(b[0], side)]
            S = math.log2(D.CELLS[R]); c = (b[1] + S) / 12
            nf = 72 * c + cfmax - (b[0] if side == 'f' else 0)
            ng = 12 * c + D.P.LOG2U6 + cgmax - (b[0] if side == 'g' else 0)
            rows.append((b[0], side, float(n), iA, kA, kB, R, nf, ng))
    return rows

def plan_retentions(slack, old=False, inflate=None):
    out = []
    for (bb, side, n, iA, kA, kB, R, nf, ng) in node_norms():
        ret, lf, lg = retention17(nf, ng, n, iA, kA, kB, slack, old=old, inflate=inflate)
        out.append((bb, side, ret, lf, lg))
    return out


def sample_side_hist(p, side, n, K, seed):
    rings = local_primes(p, K + 6)
    sqf, _ = h_factor_data(p)
    rng = random.Random(seed)
    M = p ** K
    hist = {}; rej = 0; done = 0
    while done < n:
        a = [rng.randrange(M) for _ in range(6)]
        b = [rng.randrange(M) for _ in range(6)]
        if any(not pmod(a, h, p) and not pmod(b, h, p) for h in sqf):
            rej += 1; continue
        v = 0
        for f_P, ring in rings:
            aP = ring.lift_y(a); bP = ring.lift_y(b)
            if side == 'g':
                phi = ring.add(ring.mul(aP, aP), ring.mul(ring.lift_y([0, U]), ring.mul(bP, bP)))
            else:
                ap_ = [ring.one()]; bp_ = [ring.one()]
                for _ in range(12):
                    ap_.append(ring.mul(ap_[-1], aP)); bp_.append(ring.mul(bp_[-1], bP))
                phi = ring.zero()
                for deg, c in TERMS:
                    phi = ring.add(phi, ring.mul(ring.lift_y(c), ring.mul(ap_[deg], bp_[12 - deg])))
            v += f_P * ring.val(phi)
        vv = min(v, K)
        hist[vv] = hist.get(vv, 0) + 1
        done += 1
    return hist, rej

def allowance_factor(p, side, n, K, seed):
    law, tail = exact_side_component(p, side)
    hist, rej = sample_side_hist(p, side, n, K, seed)
    worst_sigma = 0.0; infl = 1.0; rare = []
    for v, m in sorted(law.items()):
        if v == 0:
            continue
        qv = float(m)
        obs = hist.get(v, 0) / n
        sig = math.sqrt(qv * (1 - qv) / n)
        if qv * n >= 5:
            if sig > 0:
                worst_sigma = max(worst_sigma, abs(obs - qv) / sig)
            infl = max(infl, 1.0 + 3.0 * sig / qv)
        elif qv * n > 0 and obs > 0:
            rare.append((v, qv * n, hist[v]))
    qt = float(tail)
    if qt * n >= 5:
        ot = hist.get(K, 0) / n
        sigt = math.sqrt(qt * (1 - qt) / n)
        worst_sigma = max(worst_sigma, abs(ot - qt) / sigt)
        infl = max(infl, 1.0 + 3.0 * sigt / qt)
    return infl, worst_sigma, rare, hist, rej

def _lll(B):
    from decimal import Decimal as D, localcontext
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
            k = max(1, k - 1); swaps += 1
    return B, swaps

def _probable_prime(n):
    if n % 2 == 0:
        return False
    for _ in range(40):
        a = random.randrange(2, n - 1)
        d = n - 1; r = 0
        while d % 2 == 0:
            d //= 2; r += 1
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def mul6(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    for i in range(len(c) - 1, 5, -1):
        for j in range(6):
            c[i - 6 + j] -= c[i] * H[j]
    return (c + [0] * 6)[:6]

def norm6(a):
    cols = [mul6(a, [0] * i + [1]) for i in range(6)]
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
                m[i][j] = (m[i][j] * m[k][k] - m[i][k] * m[k][j]) // prev
        prev = m[k][k]
    return abs(sign * m[-1][-1])

def small_units(bound=2):
    import itertools
    out = []
    for coeffs in itertools.product(range(-bound, bound + 1), repeat=6):
        if any(coeffs) and norm6(list(coeffs)) == 1 and list(coeffs) not in ([1, 0, 0, 0, 0, 0], [-1, 0, 0, 0, 0, 0]):
            out.append(list(coeffs))
    return out

def _roots_mod(f, q, rng):
    f = _trim([x % q for x in f])
    if len(f) <= 1:
        return []
    if q > 2:
        f = pgcd(f, padd(RT.ppowmod([0, 1], q, f, q), pscal([0, 1], q - 1, q), q), q)
    if len(f) <= 1:
        return []
    out = []
    def split(a, rng):
        if len(a) - 1 == 1:
            out.append((-a[0]) * pow(a[1], -1, q) % q); return
        for _ in range(10000):
            r = _trim([rng.randrange(q) for _ in range(len(a) - 1)])
            if len(r) < 2:
                continue
            if q == 2:
                t = r; s = r
                for _ in range(len(a) - 2):
                    t = pmod(pmul(t, t, q), a, q); s = padd(s, t, q)
                h = pgcd(a, s, q)
            else:
                h = pgcd(a, padd(RT.ppowmod(r, (q - 1) // 2, a, q), [q - 1], q), q)
            if 0 < len(h) - 1 < len(a) - 1:
                split(h, rng); split(pdiv(a, h, q)[0], rng); return
    split(f, rng)
    return sorted(set(out))

def specialq_lattice(bits, side, seed):
    rng = random.Random(seed)
    while True:
        q = rng.getrandbits(bits) | (1 << (bits - 1)) | 1
        if not _probable_prime(q):
            continue
        ss = _roots_mod(H, q, rng)
        opts = []
        for s in ss:
            if side == 'g':
                ff = [U % q * s % q, 0, 1]
            else:
                ff = [(1 + s - s * s - s ** 3) % q, 0, s ** 5 % q, 0, 0, 0, 2 * s ** 3 % q, 0, 0, 0, (-2 * s) % q, 0, 1]
            opts.extend((s, t) for t in _roots_mod(ff, q, rng))
        if opts:
            break
    s, t = rng.choice(opts)
    w = [pow(s, j, q) for j in range(6)] + [(-t * pow(s, j, q)) % q for j in range(6)]
    B = [[q] + [0] * 11] + [[(-w[j]) % q] + [int(k == j) for k in range(1, 12)] for j in range(1, 12)]
    B, swaps = _lll(B)
    assert all(sum(a * b for a, b in zip(v, w)) % q == 0 for v in B)
    return q, s, t, w, B, swaps

@functools.lru_cache(None)
def _units(bound):
    return tuple(tuple(u) for u in small_units(bound))

def udup_measurement(bits, side, R, seed, points=120, bound=2):
    q, s, t, w, B, swaps = specialq_lattice(bits, side, seed)
    units = _units(bound)
    n = 12
    from fractions import Fraction as Fr_
    a = [[Fr_(x) for x in row] + [Fr_(int(i == j)) for j in range(n)] for i, row in enumerate(B)]
    for k in range(n):
        piv = next(i for i in range(k, n) if a[i][k])
        a[k], a[piv] = a[piv], a[k]
        p_ = a[k][k]
        a[k] = [x / p_ for x in a[k]]
        for i in range(n):
            if i != k and a[i][k]:
                f = a[i][k]; a[i] = [x - f * y for x, y in zip(a[i], a[k])]
    Binv = [r[n:] for r in a]
    den = q
    adj = [[int(row[j] * den) for j in range(n)] for row in Binv]
    
    
    
    
    
    Ks = []
    for eps in units:
        M_ = [[0] * 12 for _ in range(12)]   
        for col in range(6):
            ec = mul6([0] * col + [1], list(eps))
            for r in range(6):
                M_[r][col] = ec[r]
                M_[6 + r][6 + col] = ec[r]
        
        K = [[0] * 12 for _ in range(12)]
        for i in range(12):
            for j in range(12):
                s = 0
                for k in range(12):
                    if adj[k][i]:
                        for l in range(12):
                            if M_[k][l] and B[j][l]:
                                s += adj[k][i] * M_[k][l] * B[j][l]
                K[i][j] = s
        Ks.append(K)
    Rden2 = (R * den) ** 2
    rng = random.Random(seed + 1)
    counts = []
    checked = 0
    for _ in range(points):
        c = ball_point(rng, R)
        m = 0
        for K in Ks:
            cm = [sum(K[i][j] * c[j] for j in range(12)) for i in range(12)]
            if checked < 4:
                assert all(x % den == 0 for x in cm), 'mate not in the lattice'
                checked += 1
            if sum(x * x for x in cm) <= Rden2:
                m += 1
        counts.append(m)
    mean_m = sum(counts) / points
    d_r = sum(1 / (1 + m) for m in counts) / points
    var = sum((1 / (1 + m) - d_r) ** 2 for m in counts) / (points - 1)
    se = math.sqrt(var / points)
    return dict(bits=bits, side=side, R=R, q=q, points=points, units=len(units),
                mean_mates=mean_m, distinct_raw=d_r, se3=3 * se, counts=counts, swaps=swaps)

def ball_point(rng, Rr, dim=12):
    while True:
        g = [rng.gauss(0, 1) for _ in range(dim)]
        nrm = math.sqrt(sum(x * x for x in g))
        r = Rr * rng.random() ** (1.0 / dim)
        c = [int(round(x * r / nrm)) for x in g]
        if any(c) and sum(x * x for x in c) <= Rr * Rr:
            for x in reversed(c):
                if x:
                    if x < 0:
                        c = [-y for y in c]
                    break
            return c

