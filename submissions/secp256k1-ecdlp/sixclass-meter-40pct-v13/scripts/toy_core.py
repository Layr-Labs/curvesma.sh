"""Toy machinery for the six-class (negation + GLV) rho at verifier scale.

Everything is stdlib integer arithmetic. The walk is the production walk with
toy parameters: class key = min over the six-element unit orbit, twisted table,
2-cycle/checkpoint detection with a doubling escape, DP reports routed to
2^node_bits storage nodes with a per-node cap.
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
    def __init__(self, p, r, b, beta, lam):
        self.p, self.r, self.b, self.beta, self.lam = p, r, b, beta, lam

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

    # --- exact class key: minimum x over the orbit, y forced even
    def exact_key(self, P):
        p, beta = self.p, self.beta
        x, y = P
        x2 = beta * x % p
        x3 = (-x - x2) % p
        best = min(x, x2, x3)
        yy = y if y % 2 == 0 else p - y
        return (best, 0), (best, yy)

    # --- cheap rule: k-bit truncations of the three x-images with a metered
    #     exact fallback when the two smallest are within `gap` ulps
    def truncated_key(self, P, k, gap):
        p, beta = self.p, self.beta
        x, y = P
        shift = max(p.bit_length() - k, 0)
        x2 = beta * x % p
        x3 = (-x - x2) % p
        t = sorted((xx >> shift, i) for i, xx in enumerate((x, x2, x3)))
        if t[1][0] - t[0][0] <= gap:
            key, canon = self.exact_key(P)
            return key, canon, True
        xb = (x, x2, x3)[t[0][1]]
        yy = y if y % 2 == 0 else p - y
        return (xb, 0), (xb, yy), False


def base_point(cur, rnd):
    while True:
        x = rnd.randrange(cur.p)
        rhs = (x * x * x + cur.b) % cur.p
        y = pow(rhs, (cur.p + 1) // 4, cur.p) if cur.p % 4 == 3 else None
        if y is None:
            # Tonelli-Shanks is unnecessary: p == 1 mod 3 here, so trial y.
            for yy in range(2, 200):
                if yy * yy % cur.p == rhs:
                    y = yy
                    break
        if y is not None and (y * y - rhs) % cur.p == 0:
            P = cur.mul(cur.r // cur.r, (x, y))
            if P is not None and cur.is_on(P):
                return P


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
    """Twisted table: six alignments of each of R jumps."""
    U = [1, cur.lam, cur.lam * cur.lam % cur.r, cur.r - 1,
         (cur.r - cur.lam) % cur.r, (cur.r - cur.lam * cur.lam) % cur.r]
    tbl = []
    for _ in range(R):
        P = cur.mul(rnd.randrange(1, cur.r), base)
        tbl.append([cur.mul(u, P) for u in U])
    return tbl


def walk(cur, tbl, R, base, rnd, max_steps, escapes=True,
         escape_mode='doubling', additive=None, dp_bits=0, node_bits=0,
         node_cap=None, ring=8):
    """One trail of the production walk at toy scale."""
    logR = R.bit_length() - 1
    P = cur.mul(rnd.randrange(1, cur.r), base)
    key, P = cur.exact_key(P)
    seen = {key: 0}
    prev_key = None
    ringbuf = [key]
    steps = 0
    dbl_events = 0
    reports = []
    t_c = None
    nonfruitless = 0
    node_load = {}
    overflow = 0
    while steps < max_steps:
        xc = key[0]
        j = xc & (R - 1)
        twist = (xc >> (logR + 8)) % 6
        P = cur.add(P, tbl[j][twist])
        steps += 1
        nkey, P = cur.exact_key(P)
        if dp_bits and (nkey[0] & ((1 << dp_bits) - 1)) == 0:
            node = (nkey[0] >> 64) & ((1 << node_bits) - 1) if node_bits else 0
            reports.append(node)
            if node_cap is not None:
                node_load[node] = node_load.get(node, 0) + 1
                if node_load[node] > node_cap:
                    overflow += 1
        if nkey in seen and escapes:
            if t_c is None:
                t_c = steps
            if nkey == prev_key or nkey in ringbuf:
                if nkey != prev_key:
                    nonfruitless += 1
                if escape_mode == 'doubling':
                    P = cur.dbl(P)
                    steps += 2
                else:
                    P = cur.add(P, additive)
                    steps += 1
                dbl_events += 1
                nkey, P = cur.exact_key(P)
                seen = {nkey: steps}
                ringbuf = [nkey]
                prev_key = nkey
                continue
        elif nkey in seen:
            t_c = steps
            break
        seen[nkey] = steps
        ringbuf.append(nkey)
        if len(ringbuf) > ring:
            ringbuf.pop(0)
        prev_key = key
        key = nkey
    peak = max(node_load.values()) if node_load else 0
    return dict(steps=steps, t_c=t_c, doublings=dbl_events, reports=reports,
                nonfruitless=nonfruitless, overflow=overflow, peak=peak)
