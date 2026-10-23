"""Toy machinery for the six-class (negation + GLV) rho at verifier scale.

This is the production mechanism, not a stand-in for it:

* the walk is a map on classes.  The twist word is the 3-bit alignment unit of
  the current point (the table stores all six unit multiples of every jump), so
  two trails in the same class add the same aligned jump and land in the same
  class;
* escapes come only from the production detectors - a two-step key compare (the
  2-cycle rule) and a 64-step checkpoint window with the traversal-minimum rule
  - and both double the canonical representative of the *minimum-fingerprint*
  class member, then clear the history;
* distinguished-point records carry the exact exponent pair (alpha, beta) with
  Q = [alpha]P0 + [beta]T, so a key match in the shared table yields a
  challenge logarithm through x = (u beta2 - beta1)/(alpha1 - u alpha2), and the
  candidate is verified as [x]P0 == T.

A full per-trail seen set exists nowhere in the walk.  The census harness owns a
global class->owner map, but that map is instrumentation for the all-arrival
statistic T_c: the walk never reads it, and every escape in this file is decided
by the two detectors above.
"""
import random
from math import isqrt


# ---------------------------------------------------------------- primality
def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in small:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small + (41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97):
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


def cube_root_unity(m):
    """Element of order 3 modulo prime m (b^2+b+1 == 0), or None."""
    if (m - 1) % 3:
        return None
    e = (m - 1) // 3
    for g in range(2, 40):
        c = pow(g, e, m)
        if c != 1 and (c * c + c + 1) % m == 0:
            return c
    return None


# ------------------------------------------------------------- group arithmetic
class Curve:
    __slots__ = ('p', 'r', 'b', 'beta', 'lam', 'lam2', 'U', 'UIDX', 'UINV')

    def __init__(self, p, r, b, beta, lam):
        self.p, self.r, self.b, self.beta, self.lam = p, r, b, beta, lam
        self.lam2 = lam * lam % r
        U = [1, lam, self.lam2, r - 1, (r - lam) % r, (r - self.lam2) % r]
        self.U = U
        self.UIDX = {u: i for i, u in enumerate(U)}
        inv = []
        for u in U:
            v = pow(u, r - 2, r)
            assert v in self.UIDX, (u, v)
            inv.append(self.UIDX[v])
        self.UINV = inv

    def is_on(self, P):
        if P is None:
            return True
        x, y = P
        return (y * y - x * x * x - self.b) % self.p == 0

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            m = (3 * x1 * x1) % p * pow(2 * y1, p - 2, p) % p
        else:
            m = (y2 - y1) * pow((x2 - x1) % p, p - 2, p) % p
        x3 = (m * m - x1 - x2) % p
        return (x3, (m * (x1 - x3) - y1) % p)

    def dbl(self, P):
        return self.add(P, P)

    def mul(self, k, P):
        R, Q = None, P
        k %= self.r
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.dbl(Q)
            k >>= 1
        return R

    def phi(self, P):
        return None if P is None else (self.beta * P[0] % self.p, P[1])

    # --- class key and alignment -------------------------------------------
    def class_of(self, Q):
        """(min x over the six-element orbit, index of u with Q = [u]c).

        c is the canonical representative: minimum x on the orbit, even y.
        phi(x, y) = (beta x, y) = [lam](x, y), so phi^i(Q) has x = beta^i x and
        the same y; the orbit representative with the minimum x is (min x, y),
        and the sign is fixed by the parity of y.
        """
        p, beta = self.p, self.beta
        x, y = Q
        x2 = beta * x % p
        x3 = (-x - x2) % p
        if x <= x2:
            if x <= x3:
                m, i = x, 0
            else:
                m, i = x3, 2
        elif x2 <= x3:
            m, i = x2, 1
        else:
            m, i = x3, 2
        j = (-i) % 3
        u = (1, self.lam, self.lam2)[j]
        if y & 1:
            u = self.r - u
        return m, self.UIDX[u]

    def canonical(self, m, Q):
        """Canonical representative (m, even y) of the class of Q."""
        y = Q[1]
        return (m, y if not (y & 1) else self.p - y)

    def exact_key(self, P):
        """The rule the ledger's metered fallback recomputes exactly."""
        p, beta = self.p, self.beta
        x, y = P
        x2 = beta * x % p
        x3 = (-x - x2) % p
        best = min(x, x2, x3)
        yy = y if y % 2 == 0 else p - y
        return (best, yy)


# ---------------------------------------------------------- toy curve search
def base_point(cur, rnd):
    while True:
        x = rnd.randrange(cur.p)
        rhs = (x * x * x + cur.b) % cur.p
        y = None
        if cur.p % 4 == 3:
            y = pow(rhs, (cur.p + 1) // 4, cur.p)
        else:
            for yy in range(2, 400):
                if yy * yy % cur.p == rhs:
                    y = yy
                    break
        if y is not None and (y * y - rhs) % cur.p == 0:
            return (x, y)


def find_curve(bits, start_seed):
    """j=0 curve y^2=x^3+7 over F_p with prime order r == 1 (mod 6), plus GLV."""
    rnd = random.Random(start_seed)
    p = (1 << (bits - 1)) + (rnd.randrange(1 << max(bits - 12, 1)) | 1)
    while p < (1 << bits):
        if p % 3 == 1 and is_prime(p):
            four_p = 4 * p
            v = 1
            while 3 * v * v < four_p:
                t = four_p - 3 * v * v
                a = isqrt(t)
                if a * a == t:
                    for cand in (p + 1 - a, p + 1 + a):
                        if cand % 6 == 1 and is_prime(cand):
                            beta = cube_root_unity(p)
                            lam = cube_root_unity(cand)
                            if beta and lam:
                                cur = Curve(p, cand, 7, beta, lam)
                                # orient lambda so that phi(P) == [lam]P exactly
                                G = base_point(cur, random.Random(7))
                                if cur.phi(G) != cur.mul(lam, G):
                                    cur = Curve(p, cand, 7, beta, lam * lam % cand)
                                return cur
                v += 1
        p += 2
    raise RuntimeError('no curve')


def check_glv(cur, rnd, trials=200):
    """phi(P) == [lam]P on random points, and beta^3 == lam^3 == 1."""
    assert pow(cur.beta, 3, cur.p) == 1 and cur.beta != 1
    assert pow(cur.lam, 3, cur.r) == 1 and cur.lam != 1
    G = base_point(cur, rnd)
    for _ in range(trials):
        P = cur.mul(rnd.randrange(1, cur.r), G)
        if cur.phi(P) != cur.mul(cur.lam, P):
            return False
    return True


def build_table(cur, R, base, rnd):
    """Twisted table: the six alignments [u*e_j]P0 of each of R jumps."""
    tbl = []
    for _ in range(R):
        e = rnd.randrange(1, cur.r)
        row = []
        for u in cur.U:
            k = u * e % cur.r
            row.append((cur.mul(k, base), k))
        tbl.append(row)
    return tbl


# ------------------------------------------------------- truncated index rule
def truncated_index(p, beta, x, gap=40, guard=8):
    """The production rule: 8 limb products, a carry fold and a gap test.

    beta_hat = floor(2^384 beta / p) has six 64-bit limbs; x has four.  The
    64-bit fraction of beta*x is limb 5 of x*beta_hat, computed from the four
    column-5 low halves and the four column-4 high halves (eight limb products).
    The fraction of x is the top limb of x, and the fraction of -x-beta*x is the
    two's-complement fold (v0+v1+v2 in {p, 2p}).

    Returns (argmin index, fallback flag, [f0,f1,f2], [u0,u1,u2], carries).
    """
    M = (1 << 64) - 1
    bh = (1 << 384) * beta // p
    xl = [(x >> (64 * i)) & M for i in range(4)]
    bl = [(bh >> (64 * j)) & M for j in range(6)]
    s = 0
    for i in range(4):
        for j in range(6):
            c = i + j
            if c == 5:
                s += (xl[i] * bl[j]) & M
            elif c == 4:
                s += (xl[i] * bl[j]) >> 64
    f0 = xl[3]
    f1 = s & M
    f2 = (-(f0 + f1)) & M
    # exact 64-bit fractions of the three orbit x-values
    x2 = beta * x % p
    x3 = (-x - x2) % p
    u = [(v << 64) // p for v in (x, x2, x3)]
    comp = [f0, f1, f2]
    g = sorted(comp)
    tie = g[1] - g[0] <= gap
    wrap = any(v <= guard or v >= (1 << 64) - guard for v in comp)
    fb = tie or wrap
    if fb:
        a = min(range(3), key=lambda i: (u[i], i))
    else:
        a = min(range(3), key=lambda i: (comp[i], i))
    return a, fb, comp, u, tie, wrap


def _exact_argmin(u):
    return min(range(3), key=lambda i: (u[i], i))


# --------------------------------------------------------------- the trail
def point_of(cur, P0, T, alpha, beta):
    """Q = [alpha]T + [beta]P0, the convention the recovery formula uses."""
    return cur.add(cur.mul(alpha, T), cur.mul(beta, P0))


def start_trail(cur, P0, T, alpha, beta):
    return Trail(cur, point_of(cur, P0, T, alpha, beta), alpha, beta)


class Trail:
    """Production walk state: point, class key, alignment, two-key history.

    Convention: Q = [alpha]T + [beta]P0, so the recovery formula is
    x = (u beta2 - beta1)/(alpha1 - u alpha2) and the u = 1 case is exactly
    x = (beta2 - beta1)/(alpha1 - alpha2).

    The history is O(1): the key, the key one step back, the key two steps back,
    the traversal minimum of the current 64-step window, a repeat flag and a
    window counter.  There is no set and no dictionary.
    """
    __slots__ = ('Q', 'K', 'u', 'k1', 'k2', 'zmin', 'zrep', 'ck', 'alpha',
                 'beta', 'esc2', 'escchk', 'steps')

    def __init__(self, cur, Q, alpha, beta):
        self.Q = Q
        self.alpha = alpha % cur.r
        self.beta = beta % cur.r
        self.K, self.u = cur.class_of(Q)
        self.k1 = None
        self.k2 = None
        self.zmin = self.K
        self.zrep = False
        self.ck = 0
        self.esc2 = 0
        self.escchk = 0
        self.steps = 0

    def escape(self, cur, kind):
        """Double the canonical representative of the current (minimum) class."""
        m = self.K
        c = cur.canonical(m, self.Q)
        uinv = cur.U[cur.UINV[self.u]]
        a = self.alpha * uinv % cur.r
        b = self.beta * uinv % cur.r
        Q = cur.dbl(c)
        self.alpha = 2 * a % cur.r
        self.beta = 2 * b % cur.r
        self.steps += 2                    # 7M + 6A + inversion share + control
        self.Q = Q
        self.K, self.u = cur.class_of(Q)
        self.k1 = self.k2 = None
        self.zmin = self.K
        self.zrep = False
        self.ck = 0
        if kind == 0:
            self.esc2 += 1
        else:
            self.escchk += 1

    def step(self, cur, tbl, R):
        """One production step; returns 1 if a checkpoint was reached."""
        j = self.K & (R - 1)
        pt, e = tbl[j][self.u]
        Q = cur.add(self.Q, pt)
        self.beta = (self.beta + e) % cur.r     # the jump is a multiple of P0
        self.steps += 1
        if Q is None:
            return 2                            # exceptional addition: restart
        self.Q = Q
        m, u = cur.class_of(self.Q)
        k2, k1 = self.k2, self.k1
        if k2 is not None and m == k2 and m < k1:
            self.K, self.u = m, u
            self.escape(cur, 0)            # same-order rule: K is the minimum
            return 0
        self.K, self.u = m, u
        self.k2, self.k1 = k1, m
        if m < self.zmin:
            self.zmin = m
            self.zrep = False
        elif m == self.zmin:
            self.zrep = True
        self.ck += 1
        if self.ck == 64:
            due = self.zrep and m == self.zmin
            self.ck = 0
            self.zmin = m
            self.zrep = False
            if due:
                self.escape(cur, 1)        # m is the window's minimum member
                return 0
            return 1
        return 0
