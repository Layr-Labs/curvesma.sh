#!/usr/bin/env python3
"""v10b: v10's protocol at 4x the trials per density, both arms, cycle census.

Protocol is byte-for-byte v10's (same N_d = 784*4^d, same prime near N, same
R=4096 table, same four trails per group, same w_trail = max(S//4, 8*theta),
same exact exponent solve and planted-secret check).  The ONLY additions:

  * groups per density raised 400/400/400/300 -> 1600 everywhere (4x+);
  * fresh seeds (v10 used 500+d, this run uses 900+d) so the targets are fresh;
  * a CYCLE CENSUS: each trail records whether a class key repeats inside the
    trail.  Because the step is a function of the key alone, a repeated key is
    an exact trajectory cycle: without an escape rule the trail can produce no
    new distinguished point afterwards.  This is observation only -- the walk
    keeps stepping exactly as v10's protocol does, so recovery is unchanged.

Both arms are the SAME walk with the look-back window enabled/disabled
(escape.WindowRing vs None), so the paired comparison is exact.
"""
import json, random, sys, time
from math import sqrt, exp
import constants as K
import escape
from multiprocessing import Pool

MASK64 = (1 << 64) - 1
FIX = {}


def mix(x):
    x = (x + 0x9E3779B97F4A7C15) & MASK64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & MASK64
    return x ^ (x >> 31)


class Walk:
    """One trail of the class walk (verbatim v10, plus a cycle census)."""

    def __init__(self, N, R, r, table, secret, a, b, window=True, census=False):
        self.N = N; self.R = R; self.r = r; self.table = table
        self.secret = secret; self.a = a % r; self.b = b % r
        self.e = (self.a + self.b * secret) % r
        self.ring = escape.WindowRing() if window else None
        self.steps = 0; self.probes = 0
        self.key = mix(self.e)
        self.census = census
        self.seen = {self.key} if census else None
        self.cycle_at = -1
        self.repeats = 0

    def step(self):
        key = self.key
        row = (key >> 7) % self.R
        i0 = key % 12
        if self.ring is None:
            u, v = self.table[row]
            self._advance(u, v)
            return
        cands = []
        for p in range(escape.PROBES):
            i = (i0 + p) % 12
            u, v = self.table[(row + i) % self.R]
            e2 = (self.e + u + v * self.secret) % self.r
            cands.append((e2, (u, v)))
        for p in range(escape.PROBES):          # deterministic probe order
            e2, uv = cands[p]
            if not self.ring.contains(e2):
                if p:
                    self.probes += 1
                self.ring.push(e2)
                self._advance(*uv)
                return
        self.ring.blocked_all += 1
        self.ring.push(cands[0][0])
        self._advance(*cands[0][1])

    def _advance(self, u, v):
        self.a = (self.a + u) % self.r
        self.b = (self.b + v) % self.r
        self.e = (self.e + u + v * self.secret) % self.r
        self.key = mix(self.e)
        self.steps += 1
        if self.seen is not None and self.cycle_at < 0:
            if self.key in self.seen:
                self.cycle_at = self.steps
                self.repeats += 1
            else:
                self.seen.add(self.key)


def run_group(N, R, r, dp_bits, secret, rng, w_trail, n_trails, window, table, census):
    theta = 1 << dp_bits
    store = {}
    reps = hits = self_hits = cross_hits = probes = 0
    cyc_trails = cyc_steps = cyc_lost = 0
    got = False
    for seg in range(n_trails):
        w = Walk(N, R, r, table, secret, rng.randrange(1, r), rng.randrange(1, r),
                 window, census)
        for _ in range(w_trail):
            w.step()
            k = w.key
            if (k & (theta - 1)) == 0:
                reps += 1
                if k in store:
                    hits += 1
                    (a2, b2, seg2) = store[k]
                    if seg2 == seg:
                        self_hits += 1
                    else:
                        cross_hits += 1
                    den = (b2 - w.b) % r
                    if den:
                        s = (w.a - a2) * pow(den, r - 2, r) % r
                        if s == secret:
                            got = True
                            break
                else:
                    store[k] = (w.a, w.b, seg)
        probes += w.probes
        if w.cycle_at >= 0:
            cyc_trails += 1
            cyc_steps += w.cycle_at
            cyc_lost += (w_trail - w.cycle_at)
        if got:
            break
    return dict(recovered=int(got), reports=reps, hits=hits, self_hits=self_hits,
                cross_hits=cross_hits, probes=probes, store=len(store),
                cyc_trails=cyc_trails, cyc_steps=cyc_steps, cyc_lost=cyc_lost)


def _prime_r(N):
    r = N | 1
    while not _is_prime(r):
        r += 2
    return r


def _is_prime(n):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d = n - 1; s = 0
    while d % 2 == 0:
        d //= 2; s += 1
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


def setup(d):
    """Fixture per density: identical to v10's cell()."""
    N = 784 * (1 << (2 * d))          # sqrt(N) = 28*2^d, S = 32*2^d, trail = 8*2^d
    r = _prime_r(N)
    S = K.step_budget(n=r)
    n_trails = 4
    w_trail = max(S // n_trails, 8 * (1 << d))
    FIX[d] = dict(N=N, r=r, S=int(S), n_trails=n_trails, w_trail=w_trail,
                  theta=1 << d, budget=n_trails * w_trail)


def worker(arg):
    d, window, gi, seed, census = arg
    if d not in FIX:            # spawn start method: worker rebuilds its fixture
        setup(d)
    f = FIX[d]
    rng = random.Random(seed * 1000003 + gi)
    secret = rng.randrange(1, f['r'])
    table = [(rng.randrange(1, f['r']), rng.randrange(1, f['r'])) for _ in range(4096)]
    res = run_group(f['N'], 4096, f['r'], d, secret, rng, f['w_trail'],
                    f['n_trails'], window, table, census)
    res.update(dp_bits=d, N=f['N'], S=f['S'], group=gi, w_trail=f['w_trail'],
               budget=f['budget'], dp_intervals=round(f['w_trail'] / f['theta'], 2),
               budget_over_S=round(f['budget'] / f['S'], 4), seed=seed,
               arm=("window" if window else "control"))
    return res


def cholower(x, n, alpha=0.05):
    """Exact one-sided Clopper-Pearson LOWER bound, integer/rational arithmetic.

    Solves  sum_{k=x}^{n} C(n,k) p^k (1-p)^(n-k) = alpha  for p by bisection on
    exact Fractions evaluated at rationals with denominator 2^60; monotone in p,
    so bisection is exact to that grid.  x=0 gives 0.
    """
    from fractions import Fraction
    if x == 0:
        return Fraction(0)
    from math import comb

    def tail(p):
        q = 1 - p
        s = Fraction(0)
        # only terms with non-negligible mass, but keep it exact: n is <= 6400
        for k in range(x, n + 1):
            s += comb(n, k) * p ** k * q ** (n - k)
        return s

    lo, hi = Fraction(0), Fraction(x, n)
    for _ in range(60):
        mid = (lo + hi) / 2
        mid = Fraction(int(mid * (1 << 40)), 1 << 40)   # dyadic grid, exact
        if tail(mid) > alpha:
            lo = mid
        else:
            hi = mid
        if hi - lo < Fraction(1, 1 << 30):
            break
    return lo


def main():
    groups = int(sys.argv[1]) if len(sys.argv) > 1 else 1600
    census = len(sys.argv) < 3 or sys.argv[2] != 'nocensus'
    ds = (5, 7, 9, 11)
    for d in ds:
        setup(d)
    jobs = []
    for d in ds:
        for window in (True, False):
            for gi in range(groups):
                jobs.append((d, window, gi, 900 + d, census))
    # interleave arms so a mid-run crash still leaves a balanced partial sweep
    jobs.sort(key=lambda j: (j[0], j[2], not j[1]))
    t0 = time.time()
    rows = []
    with Pool(4) as pool:
        for i, res in enumerate(pool.imap_unordered(worker, jobs, chunksize=8)):
            rows.append(res)
            if (i + 1) % 400 == 0:
                print("  %d/%d  %.1fs" % (i + 1, len(jobs), time.time() - t0), flush=True)
    rows.sort(key=lambda x: (x['dp_bits'], x['arm'], x['group']))
    json.dump({'groups': groups, 'census': census, 'seed_base': 900,
               'seconds': round(time.time() - t0, 1), 'rows': rows},
              open('trials4x-out.json', 'w'))
    print("done %.1fs -> trials4x-out.json (%d rows)" % (time.time() - t0, len(rows)))


if __name__ == '__main__':
    main()
