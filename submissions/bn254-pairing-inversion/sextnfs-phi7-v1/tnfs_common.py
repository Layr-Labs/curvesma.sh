"""Shared exact arithmetic for the SexTNFS (6,2) construction over Q(zeta_7).

Elements of Z[y]/(Phi_7) are 6-lists c_0..c_5 (sum c_i y^i, y^6 = -(1+y+...+y^5)).
A field context is a dict with p, u, EPS (twist), FCO (coefficients of f = eps^4 P(x^2/eps)).
Standard library only.
"""
import math
import random

ONE = [1, 0, 0, 0, 0, 0]
Y = [0, 1, 0, 0, 0, 0]
L2 = math.log2


def bn_P(x):
    return 36*x**4 + 36*x**3 + 24*x**2 + 6*x + 1


def bn_R(x):
    return 36*x**4 + 36*x**3 + 18*x**2 + 6*x + 1


def zy_mul(a, b, mod=None):
    c = [0]*11
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                c[i+j] += ai*bj
    for k in range(10, 6, -1):
        c[k-7] += c[k]
    c6 = c[6]
    out = [c[i] - c6 for i in range(6)]
    if mod is not None:
        out = [v % mod for v in out]
    return out


def zy_add(a, b, mod=None):
    out = [x+y for x, y in zip(a, b)]
    return [v % mod for v in out] if mod is not None else out


def zy_scale(a, s, mod=None):
    out = [x*s for x in a]
    return [v % mod for v in out] if mod is not None else out


def zy_pow(a, e, mod):
    result = ONE[:]
    base = [v % mod for v in a]
    while e:
        if e & 1:
            result = zy_mul(result, base, mod)
        base = zy_mul(base, base, mod)
        e >>= 1
    return result


def zy_conj(w, k):
    """sigma_k(w) = w(y^k), k in 1..6."""
    c7 = [0]*7
    for i, ci in enumerate(w):
        c7[(i*k) % 7] += ci
    return [c7[i] - c7[6] for i in range(6)]


def zy_norm(w):
    """N_{K_h/Q}(w) = product of the six conjugates (exact integer)."""
    prod = w[:]
    for k in range(2, 7):
        prod = zy_mul(prod, zy_conj(w, k))
    assert all(v == 0 for v in prod[1:]), "norm must be rational"
    return prod[0]


def zy_eval(c, s, q):
    return sum(ci*pow(s, i, q) for i, ci in enumerate(c)) % q


def make_context(p, u, eps):
    """f = eps^4 P(x^2/eps) = 36x^8 + 36 eps x^6 + 24 eps^2 x^4 + 6 eps^3 x^2 + eps^4 ; g = x^2 - u eps."""
    e2 = zy_mul(eps, eps); e3 = zy_mul(e2, eps); e4 = zy_mul(e3, eps)
    fco = {8: [36, 0, 0, 0, 0, 0], 6: zy_scale(eps, 36), 4: zy_scale(e2, 24), 2: zy_scale(e3, 6), 0: e4}
    return dict(p=p, u=u, EPS=eps, FCO=fco)


def check_context(ctx):
    """Phi_7 irreducible mod p (Frobenius of order 6 on y), u*eps a non-square, u a square, common root."""
    p, u, eps, fco = ctx["p"], ctx["u"], ctx["EPS"], ctx["FCO"]
    yp = zy_pow(Y, p, p); yp2 = zy_pow(yp, p, p); yp3 = zy_pow(yp2, p, p); yp6 = zy_pow(yp3, p**3, p)
    assert yp != Y and yp2 != Y and yp3 != Y and yp6 == Y
    half = (p**6 - 1)//2
    assert zy_pow(zy_scale(eps, u), half, p) == [(-1) % p, 0, 0, 0, 0, 0]
    assert zy_pow([u % p, 0, 0, 0, 0, 0], half, p) == ONE
    x2 = zy_scale(eps, u, p)
    acc = [0]*6; pw = ONE
    for k in range(0, 9, 2):
        acc = zy_add(acc, zy_mul(fco[k], pw, p), p)
        pw = zy_mul(pw, x2, p)
    assert acc == [0]*6
    return True


def norm_pair(ctx, a, b):
    """Exact |N_{K_f/Q}(a+bx)| and |N_{K_g/Q}(a+bx)|."""
    fco, eps, u = ctx["FCO"], ctx["EPS"], ctx["u"]
    na = zy_scale(a, -1)
    apow = [ONE]; bpow = [ONE]
    for _ in range(8):
        apow.append(zy_mul(apow[-1], na)); bpow.append(zy_mul(bpow[-1], b))
    F = [0]*6
    for k, c in fco.items():
        F = zy_add(F, zy_mul(c, zy_mul(apow[k], bpow[8-k])))
    nf = abs(zy_norm(F))
    G = zy_add(zy_mul(a, a), zy_scale(zy_mul(eps, zy_mul(b, b)), -u))
    ng = abs(zy_norm(G))
    return nf, ng


# ------------------------------------------------------------------ primes, smoothness
def primes_upto(n):
    sieve = bytearray([1])*(n+1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(n**0.5)+1):
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
    return [i for i in range(n+1) if sieve[i]]


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d = n-1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n-1):
            continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1:
                break
        else:
            return False
    return True


def primorial(bound):
    ps = primes_upto(bound)
    while len(ps) > 1:
        ps = [ps[i]*ps[i+1] if i+1 < len(ps) else ps[i] for i in range(0, len(ps), 2)]
    return ps[0]


def is_smooth(n, prim):
    """True iff every prime factor of n divides prim (n >= 1)."""
    if n <= 0:
        return False
    g = math.gcd(n, prim)
    while n > 1:
        if g == 1:
            return False
        n //= g
        g = math.gcd(n, g)
    return True


_SMALL_FOR_SMOOTH = None


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q-1)//2, q) == 1 else -1


def P_roots_mod(q):
    """roots theta of the BN polynomial P mod q via P = (6X^2+(3+w)X+1)(6X^2+(3-w)X+1), w^2 = -3."""
    w = sqrt_mod((-3) % q, q)
    if w is None:
        return []
    out = []
    inv12 = pow(12, -1, q)
    for sign in (1, -1):
        bq = (3 + sign*w) % q
        dr = sqrt_mod((bq*bq - 24) % q, q)
        if dr is None:
            continue
        for sg in (1, -1):
            theta = ((-bq + sg*dr)*inv12) % q
            out.append(theta)
    return sorted(set(th for th in out if (36*th**4 + 36*th**3 + 24*th**2 + 6*th + 1) % q == 0))


def max_prime_factor_bits(n, lim_primes):
    """Trial-divide |n| by the ascending primes `lim_primes` (bound LIM = lim_primes[-1]), early exit.
    Return (max_pf_bits, True) if n factors completely into primes <= LIM, where max_pf_bits is the
    bit length of the largest prime factor; else (None, False) meaning a prime factor > LIM remains.
    n is thus B-smooth for any B <= LIM iff the flag is True and max_pf_bits <= log2(B)+1 by value."""
    n = abs(n)
    if n <= 1:
        return 0, True
    LIM = lim_primes[-1]
    mx = 0
    for q in lim_primes:
        if q*q > n:
            break
        if n % q == 0:
            while n % q == 0:
                n //= q
            mx = q
    if n == 1:
        return mx.bit_length() if mx else 0, True
    # remaining n is prime (all smaller factors removed; q*q>n)
    if n <= LIM:
        return max(mx, n).bit_length(), True
    return None, False


def _pollard_rho(n, rnd):
    if n % 2 == 0:
        return 2
    while True:
        c = rnd.randrange(1, n)
        x = y = rnd.randrange(2, n)
        d = 1
        while d == 1:
            x = (x*x + c) % n
            y = (y*y + c) % n
            y = (y*y + c) % n
            d = math.gcd(abs(x - y), n)
        if d != n:
            return d


def smooth_bound(n, B, rnd):
    """True iff n (>=1) factors completely into primes <= B. Trial division to 2^16, then Pollard rho."""
    global _SMALL_FOR_SMOOTH
    if n <= 1:
        return True
    if _SMALL_FOR_SMOOTH is None:
        _SMALL_FOR_SMOOTH = primes_upto(1 << 16)
    lim = min(B, 1 << 16)
    for q in _SMALL_FOR_SMOOTH:
        if q > lim:
            break
        while n % q == 0:
            n //= q
    if n == 1:
        return True
    if (1 << 16) >= B:
        return False
    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            if m > B:
                return False
            continue
        if m > B*B*B*B:          # at least four prime factors each would have to be <= B; cheap reject by size
            pass
        d = _pollard_rho(m, rnd)
        stack.append(d)
        stack.append(m // d)
    return True


# ------------------------------------------------------------------ Dickman rho and semismooth
RHO_STEP = 1.0/400
RHO_MAX = 60.0


def build_rho():
    n1 = int(round(1/RHO_STEP))
    N = int(RHO_MAX/RHO_STEP)
    rho_ = [0.0]*(N+1)
    k = 0
    while k <= N:
        x = k*RHO_STEP
        if x <= 1:
            rho_[k] = 1.0
        elif x <= 2:
            rho_[k] = 1 - math.log(x)
        else:
            break
        k += 1
    window = sum(rho_[k-n1+1:k])
    for kk in range(k, N+1):
        x = kk*RHO_STEP
        val = RHO_STEP*(rho_[kk-n1]/2 + window)
        rho_[kk] = val/(x - RHO_STEP/2)
        window += rho_[kk] - rho_[kk-n1+1]
    return rho_


RHO = build_rho()


def rho(x):
    if x <= 1:
        return 1.0
    if x >= RHO_MAX:
        return 0.0
    k = x/RHO_STEP
    i = int(k); fr = k - i
    a, b = RHO[i], RHO[min(i+1, len(RHO)-1)]
    if a > 0 and b > 0:
        return math.exp((1-fr)*math.log(a) + fr*math.log(b))
    return a*(1-fr) + b*fr


def rho_check():
    assert abs(rho(2) - 0.306853) < 1e-3
    assert abs(rho(3)/0.0486084 - 1) < 1e-2
    assert abs(rho(5)/3.54725e-4 - 1) < 2e-2
    assert abs(rho(10)/2.77017e-11 - 1) < 5e-2
    return True


# ------------------------------------------------------------------ exact integer LLL (Cohen 2.6.7)
def lll_int(B, delta_num=3, delta_den=4):
    n = len(B); b = [list(v) for v in B]; m = len(b[0])

    def dot(x, y):
        return sum(p_*q_ for p_, q_ in zip(x, y))
    d = [1]*(n+1); lam = [[0]*n for _ in range(n)]

    def gs_init(k):
        for j in range(k+1):
            w = dot(b[k], b[j])
            for i in range(j):
                w = (d[i+1]*w - lam[k][i]*lam[j][i]) // d[i]
            if j < k:
                lam[k][j] = w
            else:
                d[k+1] = w

    def red(k, l):
        if 2*abs(lam[k][l]) > d[l+1]:
            q_ = (2*lam[k][l] + d[l+1]) // (2*d[l+1])
            bk, bl = b[k], b[l]
            for i in range(m):
                bk[i] -= q_*bl[i]
            lam[k][l] -= q_*d[l+1]
            for i in range(l):
                lam[k][i] -= q_*lam[l][i]

    def swap(k):
        b[k], b[k-1] = b[k-1], b[k]
        for j in range(k-1):
            lam[k][j], lam[k-1][j] = lam[k-1][j], lam[k][j]
        lmb = lam[k][k-1]
        Bv = (d[k-1]*d[k+1] + lmb*lmb) // d[k]
        for i in range(k+1, n):
            t_ = lam[i][k]
            lam[i][k] = (d[k+1]*lam[i][k-1] - lmb*t_) // d[k]
            lam[i][k-1] = (Bv*t_ + lmb*lam[i][k]) // d[k+1]
        d[k] = Bv
    for k in range(n):
        gs_init(k)
    k = 1
    while k < n:
        red(k, k-1)
        if delta_den*d[k+1]*d[k-1] < delta_num*d[k]*d[k] - delta_den*lam[k][k-1]*lam[k][k-1]:
            swap(k); k = max(k-1, 1)
        else:
            for l in range(k-2, -1, -1):
                red(k, l)
            k += 1
    return b


# ------------------------------------------------------------------ special-q ideals and lattices
def root_of_phi7(q):
    gen = 2
    while True:
        s = pow(gen, (q-1)//7, q)
        if s != 1:
            return s
        gen += 1


def sqrt_mod(v, q):
    if v == 0:
        return 0
    if pow(v, (q-1)//2, q) != 1:
        return None
    if q % 4 == 3:
        return pow(v, (q+1)//4, q)
    # Tonelli-Shanks
    Q = q-1; S = 0
    while Q % 2 == 0:
        Q //= 2; S += 1
    z = 2
    while pow(z, (q-1)//2, q) != q-1:
        z += 1
    M = S; c = pow(z, Q, q); t = pow(v, Q, q); R = pow(v, (Q+1)//2, q)
    while t != 1:
        i = 0; tt = t
        while tt != 1:
            tt = tt*tt % q; i += 1
        b = pow(c, 1 << (M-i-1), q)
        M = i; c = b*b % q; t = t*c % q; R = R*b % q
    return R


def ideal_roots(ctx, q, s, side):
    """roots t of g(s,x) (side 'g') or f(s,x) (side 'f') modulo q, for the degree-one prime (q, y-s) of K_h."""
    u, eps, fco = ctx["u"], ctx["EPS"], ctx["FCO"]
    if side == "g":
        v = u*zy_eval(eps, s, q) % q
        t = sqrt_mod(v, q)
        return [] if t is None else sorted({t, (q-t) % q})
    es = zy_eval(eps, s, q)
    roots = []
    # f(s,x) = es^4 P(x^2/es): x^2 = es*theta for roots theta of P mod q
    coeffs = [1, 6, 24, 36, 36]   # P in increasing degree
    for theta in range(q):
        if (36*theta**4 + 36*theta**3 + 24*theta**2 + 6*theta + 1) % q == 0:
            t = sqrt_mod(es*theta % q, q)
            if t is not None:
                roots += [t, (q-t) % q]
    return sorted(set(roots))


def special_q(ctx, bits, side, rnd):
    """A random degree-one prime ideal (q, s, t) with q = 1 mod 7 of the given bit size on the given side."""
    while True:
        q = rnd.getrandbits(bits) | (1 << (bits-1))
        q -= (q - 1) % 7
        if q.bit_length() != bits or not is_prime(q):
            continue
        s0 = root_of_phi7(q)
        for k in range(1, 7):
            s = pow(s0, k, q)
            roots = ideal_roots(ctx, q, s, side) if bits <= 20 or side == "g" else f_roots_big(ctx, q, s)
            if roots:
                return q, s, roots[rnd.randrange(len(roots))]


def f_roots_big(ctx, q, s):
    """roots of f(s,x) mod q for large q: roots of P mod q by gcd(x^q - x, P) are not needed; use
    the quadratic structure P = (6X^2+3X+1)^2 + 3X^2 = (6X^2+(3+w)X+1)(6X^2+(3-w)X+1) with w^2 = -3."""
    es = zy_eval(ctx["EPS"], s, q)
    w = sqrt_mod((-3) % q, q)
    if w is None:
        return []
    out = []
    for sign in (1, -1):
        bq = (3 + sign*w) % q
        disc = (bq*bq - 24) % q
        dr = sqrt_mod(disc, q)
        if dr is None:
            continue
        inv12 = pow(12, -1, q)
        for sg in (1, -1):
            theta = ((-bq + sg*dr) * inv12) % q
            assert (36*theta**4 + 36*theta**3 + 24*theta**2 + 6*theta + 1) % q == 0
            t = sqrt_mod(es*theta % q, q)
            if t is not None:
                out += [t, (q-t) % q]
    return sorted(set(out))


def spq_lattice(q, s, t):
    """rows of {(a,b) in Z^12 : a(s) + t b(s) = 0 mod q}; determinant q."""
    rows = [[q] + [0]*11]
    for i in range(1, 6):
        v = [0]*12; v[0] = (-pow(s, i, q)) % q; v[i] = 1; rows.append(v)
    for j in range(6):
        v = [0]*12; v[0] = (-t*pow(s, j, q)) % q; v[6+j] = 1; rows.append(v)
    return rows


def in_spq_lattice(q, s, t, a, b):
    return (zy_eval(a, s, q) + t*zy_eval(b, s, q)) % q == 0


def combine(R, cvec):
    w = [sum(cvec[i]*R[i][k] for i in range(12)) for k in range(12)]
    return w[:6], w[6:]


def boot_lattice(ctx, c):
    """{(a,b): a = b*c mod (p, Phi_7)}; determinant p^6."""
    p = ctx["p"]
    rows = []
    for i in range(6):
        v = [0]*12; v[i] = p; rows.append(v)
    cy = c
    for j in range(6):
        v = cy[:] + [0]*6; v[6+j] = 1; rows.append(v)
        cy = zy_mul(cy, Y, p)
    return rows


# ------------------------------------------------------------------ ECM stage 1 (Montgomery curves, Suyama)
def ecm_stage1(n, sigma, B1, prime_list):
    """One ECM curve, stage 1 only, with bound B1. Returns a nontrivial factor or None."""
    u_ = (sigma*sigma - 5) % n
    v_ = 4*sigma % n
    x0 = pow(u_, 3, n)
    z0 = pow(v_, 3, n)
    num = pow((v_ - u_) % n, 3, n) * (3*u_ + v_) % n
    den = 16*x0*v_ % n
    g = math.gcd(den, n)
    if 1 < g < n:
        return g
    if g == n:
        return None
    a24 = num * pow(den, -1, n) % n          # (A+2)/4

    def xdbl(x, z):
        s_ = (x + z) % n; d_ = (x - z) % n
        s2 = s_*s_ % n; d2 = d_*d_ % n
        t_ = (s2 - d2) % n
        return s2*d2 % n, t_*(d2 + a24*t_) % n

    def xadd(x1, z1, x2, z2, xd, zd):
        a_ = (x1 + z1)*(x2 - z2) % n
        b_ = (x1 - z1)*(x2 + z2) % n
        s_ = (a_ + b_) % n; d_ = (a_ - b_) % n
        return zd*s_*s_ % n, xd*d_*d_ % n

    def ladder(k, x, z):
        x1, z1 = x, z
        x2, z2 = xdbl(x, z)
        for bit in bin(k)[3:]:
            if bit == '1':
                x1, z1 = xadd(x1, z1, x2, z2, x, z)
                x2, z2 = xdbl(x2, z2)
            else:
                x2, z2 = xadd(x1, z1, x2, z2, x, z)
                x1, z1 = xdbl(x1, z1)
        return x1, z1
    x, z = x0, z0
    for pr in prime_list:
        if pr > B1:
            break
        k = pr
        while k*pr <= B1:
            k *= pr
        x, z = ladder(k, x, z)
    g = math.gcd(z, n)
    if 1 < g < n:
        return g
    return None


def ecm_schedule(n, B1, curves, prime_list, rnd):
    """Fixed schedule: returns (factor or None, curves used)."""
    for i in range(curves):
        sigma = rnd.randrange(6, n-1)
        g = ecm_stage1(n, sigma, B1, prime_list)
        if g is not None:
            return g, i+1
    return None, curves
