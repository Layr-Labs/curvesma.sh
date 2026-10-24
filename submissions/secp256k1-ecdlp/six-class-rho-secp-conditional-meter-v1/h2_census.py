"""Toy census for premises H2 and H3 of proof.md (section 4).

Part "tables" (H2): G fresh tables x S trails of m arrivals with the walk rules of
upstream-proof.md section 2; one-way random-effects ANOVA of per-trail doubling and
traversal-step counts over tables. Part "collisions" (H3): inter-trail collision
count K up to the toy budget, against the random-mapping mean n^2/(2N').
Input: JSON on stdin (see main and submission.json). Output: JSON on stdout.
"""
import json
import math
import random
import sys
import time

TOYS = {"q25": (33556867, 33563419, 1), "q29": (536873251, 536839203, 3)}
CB = 7
STRICT_CYCLE = 8  # strict mode: a detected cycle longer than this is a non-fruitless return


def sqrt_mod(a, q):
    """Square root modulo a prime q (Tonelli-Shanks); None if a is a non-residue."""
    a %= q
    if a == 0:
        return 0
    if pow(a, (q - 1) // 2, q) != 1:
        return None
    if q % 4 == 3:
        return pow(a, (q + 1) // 4, q)
    s, d = 0, q - 1
    while d % 2 == 0:
        s, d = s + 1, d // 2
    z = 2
    while pow(z, (q - 1) // 2, q) != q - 1:
        z += 1
    m, c, t, x = s, pow(z, d, q), pow(a, d, q), pow(a, (d + 1) // 2, q)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2, i = t2 * t2 % q, i + 1
        b = pow(c, 1 << (m - i - 1), q)
        m, c, t, x = i, b * b % q, t * b * b % q, x * b % q
    return x


def cube_roots(m):
    """The two nontrivial cube roots of unity modulo a prime m = 1 mod 3."""
    for g in range(2, 1000):
        c = pow(g, (m - 1) // 3, m)
        if c != 1:
            return [c, c * c % m]


class Toy:
    """Affine arithmetic and the six-class canonical form on a toy subgroup."""

    def __init__(self, name, rng):
        """Build the toy curve.

        Args:
            name: key of TOYS.
            rng: random.Random used to pick the generator.
        """
        q, n, h = TOYS[name]
        self.q, self.r = q, n // h
        while True:
            x = rng.randrange(q)
            y = sqrt_mod(x ** 3 + CB, q)
            if y is None:
                continue
            G = self.mul((x, y), h)
            if G is not None and self.mul(G, self.r) is None:
                break
        self.G = G
        for b in cube_roots(q):
            for l in cube_roots(self.r):
                if self.mul(G, l) == (G[0] * b % q, G[1]):
                    self.beta, self.lam = b, l

    def add(self, P, Q):
        """Return P + Q (None is the point at infinity)."""
        if P is None:
            return Q
        if Q is None:
            return P
        q = self.q
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % q == 0:
                return None
            l = 3 * P[0] * P[0] * pow(2 * P[1], -1, q) % q
        else:
            l = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
        x3 = (l * l - P[0] - Q[0]) % q
        return (x3, (l * (P[0] - x3) - P[1]) % q)

    def mul(self, P, k):
        """Return [k]P by double-and-add, without reducing k."""
        acc = None
        while k:
            if k & 1:
                acc = self.add(acc, P)
            P = self.add(P, P)
            k >>= 1
        return acc

    def canon(self, P):
        """Return (K, t): the class key |y| and the twist index of the held member P."""
        q, x, y = self.q, P[0], P[1]
        bx = self.beta * x % q
        b2x = (-x - bx) % q
        j, cx = 0, x
        if bx < cx:
            cx, j = bx, 1
        if b2x < cx:
            j = 2
        if q - y < y:
            return q - y, 3 + j
        return y, j

    def table(self, R, rng):
        """Draw a fresh table: row h holds the six twists (-1)^s phi^-j (M_h).

        Args:
            R: number of partition classes (a power of two).
            rng: random.Random supplying the coefficients of M_h = [c]G.
        Returns:
            list of R rows of six affine points.
        """
        q, rows = self.q, []
        while len(rows) < R:
            M = self.mul(self.G, rng.randrange(1, self.r))
            row = []
            for t in range(6):
                s, j = divmod(t, 3)
                Mt = (M[0] * pow(self.beta, (-j) % 3, q) % q, M[1])
                if s:
                    Mt = (Mt[0], (q - Mt[1]) % q)
                row.append(Mt)
            rows.append(row)
        return rows


def trail(E, tab, R, P, m, on_arrival=None, dpmask=None, shared=None, owner=None):
    """Walk one trail for at most m class arrivals with the production escape rules.

    Args:
        E: the Toy curve.
        tab: the table from Toy.table.
        R: the number of partition classes.
        P: the start point.
        m: the arrival budget (every step, doubling and traversal step is one arrival).
        on_arrival: optional callback(K) -> bool; True ends the trail (collision census).
        dpmask: optional mask; a class with (K >> log2 R) & dpmask == 0 ends the trail.
        shared: optional dict class -> owner over all trails of the table. When given,
            the trail is strict: it ends with status "merge" at a class another trail
            visited, and with status "nfr" at a non-fruitless return (a class visited
            before its last escape, or a detected cycle longer than STRICT_CYCLE).
        owner: the trail's id in shared.
    Returns:
        (doublings, traversal_steps, arrivals, status).
    """
    rmask, rbits = R - 1, R.bit_length() - 1
    k1 = k2 = cp = None
    pending = False
    dbl = trav = 0
    length = 0
    old, seg = set(), set()

    def strict_stop(K):
        if K in old:
            return "nfr"
        o = shared.get(K)
        if o is not None and o != owner:
            return "merge"
        shared[K] = owner
        seg.add(K)
        return None

    while length < m:
        K, t = E.canon(P)
        if on_arrival is not None and on_arrival(K):
            return dbl, trav, length, "stop"
        if shared is not None:
            st = strict_stop(K)
            if st:
                return dbl, trav, length, st
        if dpmask is not None and ((K >> rbits) & dpmask) == 0:
            return dbl, trav, length, "dp"
        if pending:
            P = E.add(P, P)
            dbl += 1
            old |= seg
            seg = set()
            pending, k1, k2, cp = False, None, None, None
            length += 1
            continue
        if k2 is not None and K == k2:
            if K < k1:
                P = E.add(P, P)
                dbl += 1
                old |= seg
                seg = set()
                k1 = k2 = cp = None
                length += 1
                continue
            pending = True
        elif cp is not None and K == cp:
            best, Q, n = (K, P), P, 0
            while True:
                n += 1
                if n > 64:
                    return dbl, trav, length, "trav"
                KQ, tQ = E.canon(Q)
                Tp = tab[KQ & rmask][tQ]
                if Q[0] == Tp[0]:
                    return dbl, trav, length, "exc"
                Q = E.add(Q, Tp)
                trav += 1
                length += 1
                KQ = E.canon(Q)[0]
                if on_arrival is not None and on_arrival(KQ):
                    return dbl, trav, length, "stop"
                if KQ == K:
                    break
                if KQ < best[0]:
                    best = (KQ, Q)
            if shared is not None and n > STRICT_CYCLE:
                return dbl, trav, length, "nfr"
            P = E.add(best[1], best[1])
            dbl += 1
            old |= seg
            seg = set()
            k1 = k2 = cp = None
            length += 1
            continue
        if length % 64 == 0:
            cp = K
        Tp = tab[K & rmask][t]
        if P[0] == Tp[0]:
            return dbl, trav, length, "exc"
        P = E.add(P, Tp)
        k2, k1 = k1, K
        length += 1
    return dbl, trav, length, "cap"


def anova(groups):
    """One-way random-effects ANOVA over equal-size groups.

    Args:
        groups: list of G lists, each holding S per-trail counts.
    Returns:
        dict with the grand mean, the mean squares, F, the moment estimate of the
        between-group standard deviation of the group rate relative to the mean,
        and a 95% upper bound on that relative standard deviation.
    """
    G, S = len(groups), len(groups[0])
    means = [sum(g) / S for g in groups]
    grand = sum(means) / G
    msb = S * sum((x - grand) ** 2 for x in means) / (G - 1)
    msw = sum((v - mu) ** 2 for g, mu in zip(groups, means) for v in g) / (G * (S - 1))
    F = msb / msw if msw > 0 else float("inf")
    sig2 = max(0.0, (msb - msw) / S)
    # 95% upper bound for sigma_T^2: MS_B/chi2_{0.05}(G-1)*(G-1) - MS_W/S, with the
    # Wilson-Hilferty approximation of the lower chi-square quantile.
    nu, z = G - 1, 1.6448536
    chi_lo = nu * (1 - 2 / (9 * nu) - z * math.sqrt(2 / (9 * nu))) ** 3
    sig2_hi = max(0.0, msb * nu / chi_lo - msw) / S
    return dict(groups=G, per_group=S, mean_per_trail=grand, ms_between=msb,
                ms_within=msw, F=F,
                rel_sd_table_rate_estimate=math.sqrt(sig2) / grand if grand else None,
                rel_sd_table_rate_upper95=math.sqrt(sig2_hi) / grand if grand else None)


def tables_part(E, rng, R, G, S, m, deadline, strict=True):
    """Run part "tables" and return its summary (stops early at the deadline).

    Args:
        E, rng: the toy curve and the random source.
        R, G, S, m: table rows, tables, trails per table, arrivals per trail.
        deadline: time.time() value after which no new table is started.
        strict: if True, every trail is strict (see trail) and only trails that
            reach m arrivals without a merge or a non-fruitless return are kept;
            each table draws trails until it has S kept trails (at most 8*S tries).
            This is the regime of proof.md section 3: fresh arrivals only.
    Returns:
        dict with the trail statuses, the ANOVA of doublings and traversal steps
        over tables, and the pooled doubling rate.
    """
    dbl, trav, arrivals, statuses = [], [], 0, {}
    for _ in range(G):
        if time.time() > deadline:
            break
        tab = E.table(R, rng)
        shared = {} if strict else None
        gd, gt, tries, garr = [], [], 0, 0
        while len(gd) < S and tries < 8 * S:
            tries += 1
            P = E.mul(E.G, rng.randrange(1, E.r))
            d, v, n, st = trail(E, tab, R, P, m, shared=shared, owner=tries)
            statuses[st] = statuses.get(st, 0) + 1
            if strict and st != "cap":
                continue
            gd.append(d)
            gt.append(v)
            garr += n
        if len(gd) == S:
            dbl.append(gd)
            trav.append(gt)
            arrivals += garr
    out = dict(R=R, strict=strict, tables=len(dbl), per_table=S, arrivals_per_trail=m,
               arrivals=arrivals, trail_status=statuses)
    if len(dbl) >= 2:
        out["doublings"] = anova(dbl)
        out["traversal_steps"] = anova(trav)
        tot = sum(map(sum, dbl))
        out["doubling_rate"] = tot / arrivals if arrivals else None
        out["two_cycle_model_rate"] = 1 / (6 * R)
    return out


def collisions_part(E, rng, R, trials, M, dp_bits, kappa, p0, deadline):
    """Run part "collisions" and return its summary (stops early at the deadline)."""
    N1 = (E.r - 1) // 6 * R / (R - 1)
    n = int(kappa * math.sqrt(2 * N1 * math.log(1 / (1 - p0))))
    dpmask = (1 << dp_bits) - 1
    Ks, done = [], 0
    for _ in range(trials):
        if time.time() > deadline:
            break
        tab = E.table(R, rng)
        seen = {}
        total = [0]
        K = [0]
        walks = []
        for w in range(M):
            walks.append([E.mul(E.G, rng.randrange(1, E.r)), w, False])
        nxt = [M]
        while total[0] < n:
            for wk in walks:
                if total[0] >= n:
                    break
                owner = wk[1]

                def arrive(key, wk=wk, owner=owner):
                    total[0] += 1
                    prev = seen.get(key)
                    if prev is None:
                        seen[key] = owner
                    elif prev != owner and not wk[2]:
                        K[0] += 1
                        wk[2] = True
                    return total[0] >= n

                d, v, steps, st = trail(E, tab, R, wk[0], 4096, arrive, dpmask)
                wk[0] = E.mul(E.G, rng.randrange(1, E.r))
                wk[1] = nxt[0]
                wk[2] = False
                nxt[0] += 1
        Ks.append(K[0])
        done += 1
    mean = sum(Ks) / len(Ks) if Ks else None
    return dict(trials=done, walks=M, dp_bits=dp_bits, budget_arrivals=n,
                N_prime=N1, model_mean=n * n / (2 * N1), mean_K=mean,
                max_K=max(Ks) if Ks else None,
                se_mean_K=(math.sqrt(sum((k - mean) ** 2 for k in Ks) / (len(Ks) - 1) / len(Ks))
                           if len(Ks) > 1 else None))


def main():
    """Read the JSON configuration from stdin and print the census as JSON."""
    raw = sys.stdin.read().strip()
    cfg = json.loads(raw) if raw else {}
    rng = random.Random(cfg.get("seed", 1))
    t0 = time.time()
    budget = float(cfg.get("seconds", 12))
    E = Toy(cfg.get("curve", "q25"), rng)
    out = dict(config=cfg, curve=dict(q=E.q, r=E.r, beta=E.beta, lam=E.lam))
    split = float(cfg.get("tables_share", 0.6))
    out["tables"] = tables_part(E, rng, int(cfg.get("R", 64)), int(cfg.get("tables", 48)),
                                int(cfg.get("per_table", 8)), int(cfg.get("arrivals", 256)),
                                t0 + split * budget, bool(cfg.get("strict", True)))
    out["collisions"] = collisions_part(E, rng, int(cfg.get("coll_R", 2048)),
                                        int(cfg.get("coll_trials", 60)),
                                        int(cfg.get("coll_walks", 64)),
                                        int(cfg.get("dp_bits", 5)), 1.01, 0.39,
                                        t0 + budget)
    out["seconds"] = time.time() - t0
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
