#!/usr/bin/env python3
"""dladder.py -- v24 CLASS-SPACE LADDER on the shipped walk, with the
PRODUCTION CLASS KEY, a REAL lag-3 detector, and PER-RUN pair recording.

WHAT v23 GOT WRONG (verdict 845b9bf0, points 1-3) AND WHAT THIS FILE DOES.

(1) CLASS KEY.  v23's klass(x) minimised over the six x-window residues
    +-beta^h*x.  Production's key is K = min(y, p-y) -- the orbit invariant that
    is UNIFORM on [0,p/2] for a random curve point -- and the h index is the
    argmin of the three x windows (selector.twist_index).  Two consequences v23
    got wrong, both fixed here:
      * P and -P share the x-key but sigma negates y, so v23 sent them to two
        DIFFERENT purported canonical points: its map was not a section of the
        mu6 action, so equal keys did not imply equal continuation.  Here the
        canonical point is (beta^h x, K) with K = min(y,p-y): P and -P land on
        the SAME point (asserted, control C1).
      * the x-side min-of-six is concentrated near 0 (density ~6(1-6x/p)^5), so
        v23's column m = (K*R)>>(bl-1) had the wrong law.  Here K is uniform on
        [0,p/2], so m is the production truncated quotient of a uniform key.
    The ladder key is byte-compared against the production key function on 20
    adversarial pairs, including P vs -P, the infinity-adjacent point and the
    p/2 sign boundary (key_diff + audit C1).

(2) LAG-3.  v23 tested k == prev1 or k == prev2 only (periods 1 and 2) while
    claiming a lag-3 escape.  Here prev1, prev2, prev3 are kept and all three are
    tested; every fire is attributed to the lag that caught it (fires_lag1/2/3),
    so the escape delta vs the v23 detector is MEASURED on the identical
    trajectories (esc_depth=2 reproduces v23 exactly), and a synthetic period-3
    cycle is a planted control that the old detector provably misses.

(3) d.  Production records only the 2^-d distinguished points.  v23 ran the law
    cells at d=0 and quoted production d=60.  The gate regime is measured here at
    every feasible density (gate_scan: d=0,4,8,12,16) and the collapse law is
    reported in RECORDS; the honest statement about d=60 is what that records
    axis says.

(4) w2.  Each paired trial runs TWO runs on the SAME challenge Q (fresh coins)
    and records both outcomes, so the two-run success is composed PER CHALLENGE
    (E[2p_Q - p_Q^2] measured directly), and the CI is bootstrapped over the
    paired challenge rows -- not built from single-run means.
"""
import json
import math
import sys
import time

# ---------------------------------------------------------------- arithmetic
def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
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


def sqrt_mod(a, p):
    a %= p
    if a == 0:
        return 0
    if pow(a, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(a, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c = i, b * b % p
        t, r = t * c % p, r * b % p
    return r


class Curve:
    """y^2 = x^3 + b over F_p, p = 7 mod 12, #E = N prime, gcd(N,6) = 1."""

    def __init__(self, p, b, N, beta):
        self.p, self.b, self.N, self.beta = p, b, N, beta
        self.beta2 = beta * beta % p
        self.N_CLASSES = (N - 1) // 6
        w = sqrt_mod(N - 3, N)
        for cand in ((w - 1) * pow(2, -1, N) % N, (-w - 1) * pow(2, -1, N) % N):
            if (cand * cand + cand + 1) % N == 0:
                self.omega = cand
                break
        else:
            raise AssertionError("no cube root of unity scalar")
        assert pow(self.omega, 3, N) == 1 and self.omega != 1

    def resolve_omega(self, G):
        p, N = self.p, self.N
        w = sqrt_mod(N - 3, N)
        roots = [(w - 1) * pow(2, -1, N) % N, (-w - 1) * pow(2, -1, N) % N]
        betas = [self.beta, self.beta * self.beta % p]
        for cand in roots:
            if pow(cand, 3, N) != 1 or cand == 1:
                continue
            for bt in betas:
                if self.mul(G, cand)[0] == bt * G[0] % p:
                    self.omega, self.beta = cand, bt
                    self.beta2 = bt * bt % p
                    return cand
        raise AssertionError("no cube-root eigenvalue matches a beta")

    def add_batch(self, lanes):
        p = self.p
        den = [(x2 - x1) % p for (x1, _y1, x2, _y2) in lanes]
        good = [i for i, d in enumerate(den) if d]
        out = [None] * len(lanes)
        pre = [1] * (len(good) + 1)
        for j, i in enumerate(good):
            pre[j + 1] = pre[j] * den[i] % p
        inv_all = pow(pre[-1], p - 2, p) if good else 1
        invs = {}
        for j in range(len(good) - 1, -1, -1):
            i = good[j]
            invs[i] = inv_all * pre[j] % p
            inv_all = inv_all * den[i] % p
        for i in good:
            x1, y1, x2, y2 = lanes[i]
            lam = (y2 - y1) * invs[i] % p
            x3 = (lam * lam - x1 - x2) % p
            out[i] = (x3, (lam * (x1 - x3) - y1) % p)
        return out

    def mul(self, P, k):
        R = None
        A = P
        while k:
            if k & 1:
                R = A if R is None else self.add(R, A)
            A = self.dbl(A)
            k >>= 1
        return R

    def add(self, P, Q):
        p = self.p
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if x1 == x2:
            if (y1 + y2) % p == 0:
                return None
            return self.dbl(P)
        lam = (y2 - y1) * pow((x2 - x1) % p, p - 2, p) % p
        x3 = (lam * lam - x1 - x2) % p
        return (x3, (lam * (x1 - x3) - y1) % p)

    def dbl(self, P):
        if P is None:
            return None
        p = self.p
        x, y = P
        if y == 0:
            return None
        lam = 3 * x * x % p * pow(2 * y % p, p - 2, p) % p
        x3 = (lam * lam - 2 * x) % p
        return (x3, (lam * (x - x3) - y) % p)

    def neg(self, P):
        return None if P is None else (P[0], (-P[1]) % self.p)

    def sigma(self, P, h, b):
        """mu6 action: x -> beta^h x, y -> (-1)^b y."""
        if P is None:
            return None
        x, y = P
        xb = x if h == 0 else (self.beta * x % self.p if h == 1 else self.beta2 * x % self.p)
        return (xb, (-y) % self.p if b else y)

    # ---- PRODUCTION KEY ------------------------------------------------
    def klass(self, P):
        """production key: K = min(y,p-y); h = argmin_h beta^h x (twist_index);
        b = [y > p-y].  Returns (K, h, b, eta).  klass(canon(P)) == klass(P) and
        canon(canon(P)) == canon(P): the representative is a section."""
        if P is None:
            raise ValueError("infinity has no class")
        x, y = P
        p = self.p
        K = y if y <= p - y else p - y
        b = 0 if y <= p - y else 1
        vals = (x, self.beta * x % p, self.beta2 * x % p)
        h = vals.index(min(vals))
        eta = pow(self.omega, h, self.N)
        if b:
            eta = (-eta) % self.N
        return K, h, b, eta

    def canon(self, P):
        """the unique representative of the mu6 orbit: (beta^h x, K)."""
        K, h, b, _e = self.klass(P)
        return self.sigma(P, h, b)

    def col(self, K, R):
        """PRODUCTION cell(K,R,p) = (K*R) >> (bitlen(p)-1), K = min(y,p-y) <= p/2."""
        return (K * R) >> (self.p.bit_length() - 1)


def prod_key(C, P, R):
    """The production key as the walk uses it: (column m, h, b, K)."""
    K, h, b, _e = C.klass(P)
    return (K * R) >> (C.p.bit_length() - 1), h, b, K


def find_fixture(target_classes):
    import numpy as np
    t = 6 * target_classes
    p = t + ((7 - t) % 12)
    while True:
        if is_prime(p) and p % 12 == 7:
            a = np.arange(p, dtype=np.int64)
            xs = ((a * a) % p * a) % p
            sq = np.zeros(p, dtype=bool)
            sq[(a * a) % p] = True
            for b in (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12):
                rhs = (xs + b) % p
                nzero = int(np.count_nonzero(rhs == 0))
                N = 1 + 2 * (int(np.count_nonzero(sq[rhs])) - nzero) + nzero
                if N > 1 and N % 2 and N % 3 and N % 3 == 1 and is_prime(N):
                    beta = None
                    for g in range(2, 60):
                        v = pow(g, (p - 1) // 3, p)
                        if v != 1:
                            beta = v
                            break
                    assert beta and pow(beta, 3, p) == 1
                    C = Curve(p, b, N, beta)
                    Gt = make_G(C)
                    if C.mul(Gt, N) is not None or C.mul(Gt, 1) is None:
                        continue
                    if any(C.mul(Gt, k) is None for k in (2, 3)):
                        continue
                    return C
        p += 12


FIXED_G = {}


class Table:
    """Orbit rows: J0(m) = [alpha_m]G + [d_m]Q, six mu6 slots per column."""

    def __init__(self, C, R, Q, seed):
        import hashlib
        self.C, self.R = C, R
        p, N = C.p, C.N
        self.J0, self.alpha, self.dcoef = [], [], []
        h = lambda tag, m: int.from_bytes(
            hashlib.sha256(b"%s:%d:%d" % (tag, seed, m)).digest()[:8], "big")
        G = FIXED_G[C.p]
        for m in range(R):
            # a row may not be the point at infinity (J0 = O means the row carries
            # no step); such a column is re-drawn -- the production table is
            # precomputed and simply omits degenerate rows.
            salt = 0
            while True:
                al = h(b"al%d" % salt, m) % N + 1
                dc = h(b"dc%d" % salt, m) % N + 1
                J = C.add(C.mul(G, al), C.mul(Q, dc))
                if J is not None:
                    break
                salt += 1
            self.alpha.append(al)
            self.dcoef.append(dc)
            self.J0.append(J)


def make_G(C, seed=7):
    import hashlib
    p = C.p
    i = 0
    while True:
        x = int.from_bytes(hashlib.sha256(b"G:%d:%d" % (seed, i)).digest()[:8], "big") % p
        v = (x * x % p * x + C.b) % p
        y = sqrt_mod(v, p)
        if y:
            return (x, y)
        i += 1


def sm64(v):
    M = (1 << 64) - 1
    z = (v + 0x9E3779B97F4A7C15) & M
    z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & M
    z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & M
    return (z ^ (z >> 31)) & M


def is_dp(K, d):
    if d <= 0:
        return True
    return sm64(K ^ 0xD1B54A32D192ED03) < ((1 << 64) >> d)


def run_trial(C, tab, G, Q, V, L, dgate=0, esc=True, esc_depth=3, rng=None, R=None,
              trace=None):
    """One production run: L lanes, V visits, escape armed, restart counter.

    State = (canonical point P, coefficient pair (a,c)) with P = [a]G + [c]Q.
    Canonical means klass(P) has b=0 and h=0; the class determines the next
    class, which is the quotient walk.  esc_depth=2 reproduces the v23 detector.
    """
    N = C.N
    R = R or tab.R
    a = [int(v) % N for v in rng.integers(0, N, size=L)]
    cc = [int(v) % N for v in rng.integers(0, N, size=L)]
    P, canon = [], []
    for i in range(L):
        Pt = C.add(C.mul(G, a[i]), C.mul(Q, cc[i]))
        while Pt is None:
            a[i] = int(rng.integers(0, N))
            cc[i] = int(rng.integers(0, N))
            Pt = C.add(C.mul(G, a[i]), C.mul(Q, cc[i]))
        K, h, b, e = C.klass(Pt)
        P.append(C.sigma(Pt, h, b))
        a[i] = e * a[i] % N
        cc[i] = e * cc[i] % N
        canon.append(K)
    prev1, prev2, prev3 = list(canon), list(canon), list(canon)
    esc_run = [0] * L
    seen = {}
    fires = {1: 0, 2: 0, 3: 0}
    restarts = visits = 0
    steps = max(1, V // L)
    order = list(range(L))
    for _t in range(steps):
        m = [C.col(canon[i], R) for i in range(L)]
        lanes = [(P[i][0], P[i][1], tab.J0[m[i]][0], tab.J0[m[i]][1]) for i in order]
        outs = C.add_batch(lanes)
        for pos, i in enumerate(order):
            Pn = outs[pos]
            if Pn is None:                     # J == +-P: absorbed, resample lane
                a0 = int(rng.integers(0, N)); c0 = int(rng.integers(0, N))
                P0 = C.add(C.mul(G, a0), C.mul(Q, c0))
                if P0 is None:
                    continue
                k0, h0, b0, e0 = C.klass(P0)
                P[i] = C.sigma(P0, h0, b0)
                a[i] = e0 * a0 % N; cc[i] = e0 * c0 % N
                canon[i] = k0
                esc_run[i] = 0
                restarts += 1
                prev1[i] = prev2[i] = prev3[i] = k0
                continue
            an = (a[i] + tab.alpha[m[i]]) % N
            cn = (cc[i] + tab.dcoef[m[i]]) % N
            k, h, b, e = C.klass(Pn)
            hit = []
            if esc:
                if k == prev1[i]:
                    hit.append(1)
                if k == prev2[i]:
                    hit.append(2)
                if esc_depth >= 3 and k == prev3[i]:
                    hit.append(3)
            fire = False
            if hit:
                m2 = (m[i] + 1) % R
                Pn2 = C.add(P[i], tab.J0[m2])
                k2, h2, b2, e2 = C.klass(Pn2)
                if k2 != prev1[i]:
                    Pn, an, cn, k, h, b, e = (Pn2, (a[i] + tab.alpha[m2]) % N,
                                             (cc[i] + tab.dcoef[m2]) % N,
                                             k2, h2, b2, e2)
                    fire = True
                    for lg in hit:
                        fires[lg] += 1
            P[i] = C.sigma(Pn, h, b)
            a[i] = e * an % N
            cc[i] = e * cn % N
            canon[i] = k
            esc_run[i] = esc_run[i] + 1 if fire else 0
            if esc_run[i] > 8:
                a0 = int(rng.integers(0, N)); c0 = int(rng.integers(0, N))
                P0 = C.add(C.mul(G, a0), C.mul(Q, c0))
                k0, h0, b0, e0 = C.klass(P0)
                P[i] = C.sigma(P0, h0, b0)
                a[i] = e0 * a0 % N; cc[i] = e0 * c0 % N
                canon[i] = k0
                esc_run[i] = 0
                restarts += 1
                prev1[i] = prev2[i] = prev3[i] = k0
                continue
            visits += 1
            if trace is not None:
                trace.append((i, k))
            if not is_dp(k, dgate):
                continue
            if k in seen and seen[k][0] != a[i]:
                a1, c1 = seen[k]
                for e6 in (1, C.omega, C.omega * C.omega % N,
                           N - 1, N - C.omega, N - C.omega * C.omega % N):
                    den = (cc[i] - e6 * c1) % N
                    if den == 0:
                        continue
                    xx = (e6 * a1 - a[i]) * pow(den, N - 2, N) % N
                    if C.mul(G, xx) == Q:
                        return dict(ok=True, visits=visits, fires=fires,
                                    restarts=restarts)
            if k not in seen:
                seen[k] = (a[i], cc[i])
        prev3 = list(prev2)
        prev2 = list(prev1)
        prev1 = list(canon)
    return dict(ok=False, visits=visits, fires=fires, restarts=restarts)


def cp_lower(k, n, alpha=0.05):
    """One-sided Clopper-Pearson lower bound, exact (bisection on the binomial)."""
    if n == 0 or k == 0:
        return 0.0
    lo, hi = 0.0, k / n
    for _ in range(90):
        mid = (lo + hi) / 2
        s = 0.0
        for j in range(k):
            s += math.exp(math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
                          + j * math.log(mid) + (n - j) * math.log(1 - mid))
        if 1 - s < alpha:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def cell(target_classes, trials, seed, esc=True, V_mult=1.0, dgate=0, L=16, R=32,
         esc_depth=3, pairs=False, boot=400):
    """One ladder cell.  pairs=True runs TWO runs on the SAME challenge Q and
    records the PAIRED outcomes, so w2 is composed per challenge."""
    import numpy as np
    C = find_fixture(target_classes)
    G = make_G(C)
    FIXED_G[C.p] = G
    C.resolve_omega(G)
    S_OPS = 1.009696 * math.sqrt(C.N_CLASSES)
    V = max(L, int(round(V_mult * S_OPS)))
    ok = visits_tot = rest_tot = 0
    fires_tot = {1: 0, 2: 0, 3: 0}
    rows = []
    t0 = time.time()
    rng = np.random.default_rng(seed)
    for tr in range(trials):
        x = int(rng.integers(1, C.N))
        Q = C.mul(G, x)
        tab1 = Table(C, R, Q, seed * 100003 + tr)
        r1 = run_trial(C, tab1, G, Q, V, L, dgate, esc=esc, esc_depth=esc_depth,
                       rng=rng, R=R)
        if pairs:
            tab2 = Table(C, R, Q, seed * 100003 + tr)   # same rows, fresh coins
            r2 = run_trial(C, tab2, G, Q, V, L, dgate, esc=esc, esc_depth=esc_depth,
                           rng=rng, R=R)
            rows.append((1 if r1["ok"] else 0, 1 if r2["ok"] else 0))
            ok += 1 if (r1["ok"] or r2["ok"]) else 0
            visits_tot += r1["visits"] + r2["visits"]
            rest_tot += r1["restarts"] + r2["restarts"]
            for lg in (1, 2, 3):
                fires_tot[lg] += r1["fires"][lg] + r2["fires"][lg]
        else:
            ok += 1 if r1["ok"] else 0
            visits_tot += r1["visits"]
            rest_tot += r1["restarts"]
            for lg in (1, 2, 3):
                fires_tot[lg] += r1["fires"][lg]
            rows.append((1 if r1["ok"] else 0, None))
    out = {
        "target_classes": target_classes,
        "p": C.p, "b": C.b, "N": C.N, "N_CLASSES": C.N_CLASSES,
        "classes_log2": round(math.log2(C.N_CLASSES), 4),
        "R": R, "L": L, "V": V, "trials": trials,
        "esc": esc, "esc_depth": esc_depth, "dgate": dgate, "V_mult": V_mult,
        "shape_V_over_sqrtN": V / math.sqrt(C.N_CLASSES),
        "successes": ok, "p_hat": ok / trials,
        "visits": visits_tot,
        "escape_fires": fires_tot[1] + fires_tot[2] + fires_tot[3],
        "fires_lag1": fires_tot[1], "fires_lag2": fires_tot[2],
        "fires_lag3": fires_tot[3],
        "restarts": rest_tot,
        "elapsed_s": round(time.time() - t0, 3),
        "cp95_run": cp_lower(ok, trials) if not pairs else None,
    }
    if pairs:
        n = len(rows)
        s1 = sum(r[0] for r in rows)
        s2 = sum(r[1] for r in rows)
        out["p_run1"] = s1 / n
        out["p_run2"] = s2 / n
        out["p_both"] = sum(1 for r in rows if r[0] and r[1]) / n
        out["p_either"] = sum(1 for r in rows if r[0] or r[1]) / n
        out["p_run_pooled"] = (s1 + s2) / (2 * n)
        out["w2_measured_per_challenge"] = out["p_either"]
        out["w2_model_from_means"] = 2 * out["p_run_pooled"] - out["p_run_pooled"] ** 2
        brng = np.random.default_rng(seed ^ 0xB007)
        arr = np.array(rows, dtype=np.int8)
        bs = []
        for _b in range(boot):
            idx = brng.integers(0, n, size=n)
            bs.append(float(np.mean(np.maximum(arr[idx][:, 0], arr[idx][:, 1]))))
        bs.sort()
        out["w2_boot_lo"] = bs[int(0.025 * boot)]
        out["w2_boot_hi"] = bs[int(0.975 * boot) - 1]
        out["w2_boot_mean"] = sum(bs) / boot
        out["w2_boot_se"] = float(np.std(bs, ddof=1))
        out["cp95_single_run_run1"] = cp_lower(s1, n)
        out["cp95_single_run_run2"] = cp_lower(s2, n)
        out["n_pairs"] = n
    return out


def adversarial_pairs(C, G, R=32):
    """20 pairs incl. P vs -P, the smallest-y point and the p/2 sign boundary."""
    p = C.p
    out = []
    for i, k in enumerate((1, 2, 3, 5, 8, 13, 21, 34, 55, 89)):
        P = C.mul(G, k)
        out.append(("P%d" % i, P))
        out.append(("negP%d" % i, C.neg(P)))
    yminP = None
    for x in range(1, 6000):
        v = (x * x % p * x + C.b) % p
        y = sqrt_mod(v, p)
        if y is None:
            continue
        K = min(y, p - y)
        if yminP is None or K < yminP[0]:
            yminP = (K, (x, y))
    if yminP:
        out.append(("min_K_point", yminP[1]))
        out.append(("min_K_neg", C.neg(yminP[1])))
    for x in range(1, 400000):
        v = (x * x % p * x + C.b) % p
        y = sqrt_mod(v, p)
        if y is None:
            continue
        if abs(y - (p // 2)) <= 1:
            out.append(("boundary_half", (x, y)))
            out.append(("boundary_half_neg", C.neg((x, y))))
            break
    i = 0
    while len(out) < 20:
        out.append(("extra_%d" % i, C.mul(G, 1000 + i)))
        i += 1
    return out[:20]


def key_diff(C, G, R=32):
    """C1: ladder key == production key on the 20 adversarial pairs."""
    rows = []
    for name, P in adversarial_pairs(C, G, R):
        m, h, b, K = prod_key(C, P, R)
        K2, h2, b2, _e = C.klass(P)
        m2 = (K2 * R) >> (C.p.bit_length() - 1)
        can = C.canon(P)
        can2 = C.canon(can)
        rows.append({"pair": name, "K": K, "K_ladder": K2, "h": h, "h_ladder": h2,
                     "b": b, "b_ladder": b2, "col": m, "col_ladder": m2,
                     "canon_x": can[0], "canon_y": can[1],
                     "idempotent": can2 == can,
                     "key_equal": (m == m2 and h == h2 and b == b2 and K == K2)})
    return rows


def gate_scan(target_classes, seed, trials, depths=(0, 4, 8, 12, 16), L=16, R=32):
    """(3) the gate law in RECORDS: p(d) at fixed V, and the collapse per bit."""
    rows = []
    for d in depths:
        r = cell(target_classes, trials, seed, dgate=d, L=L, R=R)
        r["records_per_run"] = round(r["visits"] / (trials * 2 ** d), 4) if d else None
        rows.append(r)
    return rows


def escape_delta(target_classes, seed, trials, L=16, R=32, pairs=True):
    """(2) v23 detector (depth 2) vs v24 (depth 3) on MATCHED seeds."""
    old = cell(target_classes, trials, seed, esc_depth=2, L=L, R=R, pairs=pairs)
    new = cell(target_classes, trials, seed, esc_depth=3, L=L, R=R, pairs=pairs)
    return old, new


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "law"
    T = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    out = {}
    if which == "law":
        out["law"] = [cell(2 ** 14, T, 11, L=16, R=32),
                      cell(2 ** 17, T, 12, L=16, R=32),
                      cell(2 ** 20, T, 13, L=16, R=32)]
    elif which == "w2":
        out["w2"] = [cell(2 ** 20, T, 13, L=16, R=32, pairs=True),
                     cell(2 ** 17, T, 12, L=16, R=32, pairs=True)]
    elif which == "gate":
        out["gate"] = gate_scan(2 ** 17, 12, T)
    elif which == "esc":
        o, n = escape_delta(2 ** 17, 12, T)
        out["esc"] = {"v23_depth2": o, "v24_depth3": n}
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
