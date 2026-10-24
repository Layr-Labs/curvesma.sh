"""qwalk.py -- v20e quotient walk: production algebra + COALESCENCE-PRESERVING TRAP ESCAPE.

The walk is the mu_6 quotient walk of v20b-d: q prime, q = 1 mod 6, lambda^3 = 1,
class(c) = min{+-lambda^k c}, one seed u[m] per column, and the step

    m = cell(canon);  t = canon + u[m]*x mod q;  canon' = canon(t);
    (a, c) <- e*(a, c + u[m]),   e = canon/(canon + u[m])

so the invariant  a + c*x == canon (mod q)  holds at every step and two records in
one class from two different walkers solve x = (a1-a2)/(c2-c1) exactly.

v20e ADDS THE COALESCENCE-PRESERVING TRAP ESCAPE.  A two-step return (the G2 trap
channel) is the state pattern s -> s1 -> s with s1 = canon(s + u[m(s)]x): the
successor of s1 is again s.  The escape is a deterministic function of the walk's
own classes:

    at class s1, if the successor would be the class s of the previous step,
    apply the ALTERNATE ORBIT ROW at s1 for this one step:
        column (m(s1) + 1) mod R instead of m(s1), same orbit construction,
    otherwise apply the normal row.

Why this preserves coalescence.  The alternate row is the same orbit-generated row
type one column over: at canonical class c the row is J = [e_c*u[m(c)]]T, so
J(sigma c) = sigma J(c) for all six sigma in mu_6 (construction checked in
scripts/medium_fixture.py), and the escape decision compares the successor CLASS
with the previous CLASS, both sigma-invariant.  The walk is never restarted, so two
walkers that occupy one class with the same predecessor keep sharing their future
exactly; the walk state for coalescence is the pair (class, previous class).

`prev` is the class one step back.  The escape fires at most once per closing
two-cycle and the walk continues.
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
        """Balanced R-way partition of the class set (class function)."""
        return (sm64(canon ^ np.uint64(0x9E3779B97F4A7C15))
                % np.uint64(self.R)).astype(np.int64)

    def dp(self, canon, d):
        """Class-function report gate (fingerprint of the canonical class key)."""
        thr = np.uint64((1 << 64) >> d if d < 64 else 0)
        return sm64(canon ^ np.uint64(0xD1B54A32D192ED03)) < thr

    # ---- rows ------------------------------------------------------------
    def row(self, canon, alt=False):
        """The orbit row applied at a canonical class.

        Returns (column, coefficient u).  The increment applied to the class value
        is u*x, i.e. J = [u]T; the representative's orbit multiplier e_c enters the
        coefficient update through `canon()`, exactly as in production where the
        slot stores sigma_{k,b}J0 and its multiple of T.
        """
        m = self.cell(canon)
        if alt:
            m = (m + 1) % self.R
        return m, self.u[m]


def sm64(a):
    a = np.asarray(a, dtype=np.uint64)
    z = (a + np.uint64(0x9E3779B97F4A7C15)) & np.uint64(MASK64)
    z = ((z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)) & np.uint64(MASK64)
    z = ((z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)) & np.uint64(MASK64)
    return z ^ (z >> np.uint64(31))


def _apply(U, canon, a, cc, x, m):
    """One step with an explicit column vector m.  Returns (canon, a, cc, e)."""
    q = U.q
    u = U.u[m]
    inc = (u * x) % q
    t = canon + inc
    t = np.where(t >= q, t - q, t)
    cn, e, _h = U.canon(t)
    a = (e * a) % q
    cc = (e * ((cc + u) % q)) % q
    return cn.astype(np.uint64), a.astype(np.uint64), cc.astype(np.uint64), e


def step(U, canon, a, cc, x, d, want_dp=False, prev=None, esc=False,
         force_every=0, tick=0, stats=None):
    """One production step for every lane.

    prev        : the class one step back (None disables the escape test)
    esc         : arm the two-cycle escape
    force_every : >0 forces the escape branch on every `force_every`-th lane/tick
                  (stress + equivariance runs only -- never in the census)
    Returns (canon, a, cc, dp_hit, prev_out, n_fired).
    """
    q = U.q
    m = U.cell(canon)
    cn, a2, cc2, _e = _apply(U, canon, a, cc, x, m)
    n_fired = 0
    if esc and prev is not None:
        fire = cn == prev
        if force_every:
            fire = fire | (((np.arange(canon.size) + tick) % force_every) == 0)
        if fire.any():
            idx = np.flatnonzero(fire)
            ma = (m[idx] + 1) % U.R
            cn2, a3, cc3, _e2 = _apply(U, canon[idx], a[idx], cc[idx], x[idx], ma)
            cn = cn.copy()
            a2 = a2.copy()
            cc2 = cc2.copy()
            cn[idx] = cn2
            a2[idx] = a3
            cc2[idx] = cc3
            n_fired = int(idx.size)
    if stats is not None:
        stats["fired"] = stats.get("fired", 0) + n_fired
        stats["steps"] = stats.get("steps", 0) + int(canon.size)
    hit = U.dp(cn, d) if want_dp else None
    return cn, a2.astype(np.uint64), cc2.astype(np.uint64), hit, canon.copy(), n_fired


def walk(U, start, x, a0, steps, d, esc=True):
    """Walk `steps` steps; return class trace, DP hit mask, escape-fire count."""
    lanes = start.size
    canon = U.canon(start)[0]
    a = a0.copy()
    cc = np.zeros(lanes, dtype=np.uint64)
    prev = canon.copy()
    trace = np.empty((steps, lanes), dtype=np.uint64)
    hits = np.zeros((steps, lanes), dtype=bool)
    fired = 0
    for r in range(steps):
        trace[r] = canon
        canon, a, cc, hit, prev, nf = step(U, canon, a, cc, x, d, want_dp=True,
                                           prev=prev, esc=esc)
        fired += nf
        hits[r] = hit
    return trace, hits, fired
