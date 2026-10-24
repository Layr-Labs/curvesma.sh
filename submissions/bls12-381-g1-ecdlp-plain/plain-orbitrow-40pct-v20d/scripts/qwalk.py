"""qwalk.py -- the v20b quotient walk: the production algebra on a class universe.

The production walk lives on the classes of the mu_6 action (v20: equivariance
6,000,000/6,000,000, coalescence 1.000 at all 101 checkpoints).  This module is
that quotient walk with the field F_p replaced by Z_q, so it can be run at a
class universe of ~2^29.4 (the rival scale) instead of 2^379:

  * q prime, q = 1 mod 6, lambda^3 = 1 mod q, lambda != 1  (the mu_6 action:
    c -> (-1)^b lambda^k c on the additive group Z_q, exactly as phi = [lam] and
    [-1] act on the order-r subgroup of E).
  * class(c) = min{ +-lambda^k c }: the canonical representative; (k,b) is the
    orbit decomposition c = (-1)^b lambda^k canon.  This is the surrogate's index
    convention; the production index (argmin of the x-residues, y sign) is a
    different bijection orbit -> Z_3 x Z_2 and only the equivariance relation
    J(sigma c) = sigma J(c) is used by the walk, which both conventions satisfy.
  * one seed u[m] per column m; the slot at a canonical rep is J0 = [u[m]]T, and
    at c = sigma(canon) it is sigma(J0) = (c/canon) J0 -- BY CONSTRUCTION, the
    same table the ledger prices (one scalar multiple of T per column).
  * step: canon <- canon(canon + J) with the coefficient transform
    (a, c) <- e (a, c + u[m]),  e = canon/(canon + J) in mu_6,
    so that the invariant  a + c x == canon (mod q)  holds at EVERY step; the
    label of a point is its discrete log (Z_q with generator 1), so a record
    pair in one class solves x = (a1-a2)/(c2-c1) and is VERIFIED against the
    trial's x.  Recovery genuinely needs two walkers with different a0, exactly
    as in production, where every J is a multiple of T.

Two-step returns (the G2 trap channel) are counted on the same step function.
"""
import hashlib

import numpy as np

MASK64 = (1 << 64) - 1


def is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
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


def universe(q_min):
    """Smallest prime q >= q_min with q = 1 mod 6, plus a nontrivial cube root."""
    q = q_min + ((1 - q_min) % 6)
    while not is_prime(q):
        q += 6
    lam = None
    for a in range(2, 200):
        v = pow(a, (q - 1) // 3, q)
        if v != 1:
            lam = v
            break
    assert lam is not None and pow(lam, 3, q) == 1 and lam != 1
    return q, lam


class Universe:
    def __init__(self, q_min, R, seed=20250):
        self.q, self.lam = universe(q_min)
        self.q_log = self.q.bit_length() - 1
        self.R = R
        self.seed = seed
        self.lam2 = self.lam * self.lam % self.q
        self.LAMH = np.array([1, self.lam, self.lam2], dtype=np.uint64)
        self.classes = (self.q - 1) // 6 + 1
        u = np.empty(R, dtype=np.uint64)
        for m in range(R):
            d = hashlib.sha256(b"v20b-col:%d:%d:%d" % (seed, self.q, m)).digest()
            u[m] = (int.from_bytes(d[:8], "big") % (self.q - 1)) + 1
        self.u = u

    # ---- canon / orbit decomposition -------------------------------------
    def canon(self, c):
        q = self.q
        c = np.asarray(c, dtype=np.uint64)
        u1 = (self.lam * c) % q
        u2 = (self.lam2 * c) % q
        w = np.minimum(np.minimum(c, u1), u2)
        W = np.maximum(np.maximum(c, u1), u2)
        pos = w <= q - W
        canon = np.where(pos, w, q - W)
        hp = np.where(c <= u1, np.where(c <= u2, 0, 2), np.where(u1 <= u2, 1, 2))
        hq = np.where(c >= u1, np.where(c >= u2, 0, 2), np.where(u1 >= u2, 1, 2))
        e = np.where(pos, self.LAMH[hp], q - self.LAMH[hq])
        return canon.astype(np.uint64), e.astype(np.uint64), hp.astype(np.int64)

    def cell(self, canon):
        """Balanced R-way partition of the class set.

        Production uses m = (K*R) >> (bitlen(p)-1) with K = min(y, p-y), which is
        balanced because y is uniform on the curve.  The canonical class key here
        is a minimum of six residues and is NOT uniform, so the surrogate uses a
        keyed hash of the class: a balanced partition of the same R-way shape
        (the mechanism needs a balanced R-column table, not this particular map).
        """
        return (sm64(canon ^ np.uint64(0x9E3779B97F4A7C15))
                % np.uint64(self.R)).astype(np.int64)

    def dp(self, canon, d):
        """Class-function report gate (fingerprint of the canonical class key)."""
        thr = np.uint64((1 << 64) >> d if d < 64 else 0)
        return sm64(canon ^ np.uint64(0xD1B54A32D192ED03)) < thr


def sm64(a):
    a = np.asarray(a, dtype=np.uint64)
    z = (a + np.uint64(0x9E3779B97F4A7C15)) & np.uint64(MASK64)
    z = ((z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)) & np.uint64(MASK64)
    z = ((z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)) & np.uint64(MASK64)
    return z ^ (z >> np.uint64(31))


def step(U, canon, a, cc, x, d, want_dp=False):
    """One production step for every lane.  Returns (canon, a, cc, dp_hit)."""
    q = U.q
    m = U.cell(canon)
    u = U.u[m.astype(np.int64)]
    inc = (u * x) % q                     # the class value advances by u*x
    t = canon + inc
    t = np.where(t >= q, t - q, t)
    cn, e, _h = U.canon(t)
    a = (e * a) % q
    cc = (e * ((cc + u) % q)) % q         # the T-coefficient advances by u
    return cn, a.astype(np.uint64), cc.astype(np.uint64), (U.dp(cn, d) if want_dp else None)


def walk(U, start, x, a0, steps, d, rows=None):
    """Walk `steps` steps; return the class trace and the DP hit mask."""
    lanes = start.size
    canon = U.canon(start)[0]
    a = a0.copy()
    cc = np.zeros(lanes, dtype=np.uint64)
    trace = np.empty((steps, lanes), dtype=np.uint64)
    hits = np.zeros((steps, lanes), dtype=bool)
    for r in range(steps):
        trace[r] = canon
        canon, a, cc, hit = step(U, canon, a, cc, x, d, want_dp=True)
        hits[r] = hit
    return trace, hits
