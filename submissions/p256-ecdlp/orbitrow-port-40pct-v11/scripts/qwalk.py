"""qwalk.py -- v20f quotient walk: production algebra + two-step-return repair.

The walk is the mu_6 quotient walk of v20b-d: q prime, q = 1 mod 6, lambda^3 = 1,
class(c) = min{+-lambda^k c}, one seed u[m] per column, and the step

    m = col(canon);  t = canon + u[m]*x mod q;  canon' = canon(t);
    (a, c) <- e*(a, c + u[m]),   e = canon/(canon + u[m])

so the invariant  a + c*x == canon (mod q)  holds at every step and two records in
one class from two different walkers solve x = (a1-a2)/(c2-c1) exactly.

THE COLUMN (v20g, reviewer point 3).  v20e bucketed the canonical value with a
SplitMix hash; v20f replaced it with a truncation of the scalar minimum and
measured that this is the WRONG model (the toy's `canon` is the minimum of six
orbit residues and is concentrated near 0, so equal-width buckets of it are
skewed: column collision 0.065 vs 1/R = 0.0039, two-step return 18.8x the 1/(6R)
law).  v20g keeps the DISTRIBUTION right -- K is uniform on [0, p/2] in
production because it is min(y, p-y) of a random curve point, and the toy feeds a
uniform class-function K -- and imports the MAP verbatim from the ledger's own
selector module: `selector.cell(K, R, P) = (K*R) >> (bitlen(P)-1)`.  The report
gate `dp()` stays a fingerprint (SplitMix of the canonical class), which is the
production shape: production tests a fingerprint of the class key, not a range of
it.

THE TWO-STEP RETURN.  A two-step return is the state pattern

    c -> c1 -> c,      c1 = canon(c + u[m(c)]x),

and it is the only structural cycle of the map: writing e1, e2 for the two orbit
multipliers, the composition returns iff  c*(e1*e2 - 1) = v1 + e1*v2, and for
e1 = e2 = -1 (the mu_6 element -1, i.e. the canonical representative is the
negation) that identity degenerates to v2 = v1, true for EVERY c.  So the channel
has per-step density

    P(two-step return) = P(e1 = -1) * P(col(c1) = col(c)) = (1/6) * (1/R_eff),

with R_eff = selector.reachable_columns(R, P) the REACHABLE columns (the class
key is K <= P/2, so the shift is bitlen(P)-1 and not every bucket is reached),
measured with per-walk counters in scripts/medium_fixture.py and in the
production replay.  REPAIR (priced on every step, no abort, no restart): at class
c1, if the successor would be the class of the previous step, apply the ALTERNATE
ORBIT ROW at c1 for this one step -- column (m+1) mod R, same construction -- and
continue.  The row is orbit-generated, so J(sigma S) = sigma J(S) for all six
sigma in mu_6, and the escape decision compares CLASS fingerprints, which are
sigma-invariant.

WHAT COALESCENCE IS AND IS NOT USED FOR (v20f).  The success event is

    two walkers report the same CLASS somewhere in their own budgets,

and the record pool is keyed by class, so the two records may be produced at
DIFFERENT steps: recovery needs only that both trajectories CONTAIN the class,
not that the walkers stay together.  That is why the predecessor channel v20e
leaned on is not part of the model: the escape can change a walker's future after
a meeting, but it cannot remove the class both walkers already visited, and the
invariant a + c*x == canon holds at every step of both walks.  The escape's only
model cost is the case where exactly one of two walkers closing on a shared class
substitutes a row there, a per-meeting event of probability equal to the trap
density 1/(6*R_eff) = 2^-22.29 at production R_eff = 852105.

`prev` is the class one step back.  v20i: the escape fires up to K_ESC_MAX times in
a row and the walk continues -- there is NO raise and NO abort.  The run COST is
priced (ledger ROWS "escape run cost M/R", the deviation cap and the meter), so a
long run is a priced event, not an assertion.  `step()` takes the caller's per-lane
run-length array (`esc_run`); when a lane would exceed `kcap = K_ESC_MAX` the step
does not throw: it RESETS that lane's run counter and reports the restart through
`restart_out` (and the module counter `RESTARTS`), which is the same fresh-start
recovery the walker fleet already performs everywhere else.  Every walker loop in
this package passes the run array (medium_fixture.run, prod_sweep.shard, the run
fast path), so no path can silently exceed the charged run length.
"""
import hashlib

import numpy as np

import selector

K_ESC_MAX = 8                  # ALLOWED run length, == ledger.K_ESC_MAX
K_ESC = K_ESC_MAX               # legacy alias (callers pass kcap=qwalk.K_ESC)
MASK64 = (1 << 64) - 1
RESTARTS = 0                   # lanes restarted after an over-long escape run
COUNTER_ONLY = False           # CONTROL ARM ONLY (tests): the v20i defect, kept so
                               # the restart-faithfulness test can PROVE it fails.
FRESH_DRAWS = 0                # coin draws consumed by fresh starts (auditable)


class K_ESCExceeded(RuntimeError):
    """Unused in v20i (kept so an old caller cannot silently mis-handle it)."""


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
    """Smallest prime q >= q_min, and lam = 1.

    P-256 RETARGET (orbit size 2).  The plain/secp lineage walks the mu_6 orbit
    quotient because secp256k1 and BLS12-381 G1 have j = 0, so x -> beta*x is a
    curve automorphism and six points share one class key.  NIST P-256 has
    a = -3, b != 0 (j = 1728), so its automorphism group is {+1, -1} ONLY: the
    retarget's class map is canon(c) = min(c, q-c), i.e. the SAME code with
    lam = 1 (u1 = u2 = c, so the min/max over the orbit degenerates to the
    signed pair and e in {+1, -1}).  Nothing in the step, the escape or the
    recovery changes; only the class-space size does, N = (q-1)/2 instead of
    (q-1)/6, and that is charged in the ledger (+0.79 bits of S0 vs secp)."""
    q = q_min
    while not is_prime(q):
        q += 1
    lam = 1
    return q, lam


class Universe:
    def __init__(self, q_min, R, seed=20250, row_mode="orbit"):
        self.q, self.lam = universe(q_min)
        self.q_log = self.q.bit_length() - 1
        self.R = R
        self.seed = seed
        self.row_mode = row_mode
        self.lam2 = self.lam * self.lam % self.q
        self.LAMH = np.array([1, self.lam, self.lam2], dtype=np.uint64)
        self.Q6 = (self.q - 1) // 2 + 1
        self.classes = self.Q6
        # REACHABLE COLUMNS.  The production column is cell(K, R, p) with
        # K = min(y, p-y) <= p/2, so only R/2 of the R buckets are reachable --
        # the ledger books that as R_EFF = selector.reachable_columns(R, p) and
        # the trap density is 1/(6*R_EFF), not 1/(6*R).  The toy has the same
        # structure (K <= P//2), so it must use the same effective count.
        self.R_eff = int(selector.reachable_columns(self.R, 1 << 40))
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

    def col(self, canon):
        """PRODUCTION SELECTOR, imported verbatim (v20g, reviewer point 3).

        Production:  K = min(y, p-y)  (orbit-invariant, UNIFORM on [0, p/2] for a
        random point of the curve),  m = selector.cell(K, R, p) = (K*R) >> (bitlen(p)-1),
        i.e. a truncated quotient -- the same function this module calls.

        The toy's `canon` is the minimum of SIX orbit residues of a uniform scalar,
        which is concentrated near 0 (density ~6(1-6x/q)^5), so feeding it to the
        truncated quotient would model the wrong distribution: measured column
        collision 0.065 against 1/R = 3.9e-3 and a two-step return rate 18.8x the
        1/(6R) law.  What the toy CAN supply is the DISTRIBUTION: K = uniform on
        [0, 2^63), a class function (a SplitMix image of the class key).  So the
        map is the ledger's own cell() and the input is uniform, which is exactly
        the production pair (uniform K, truncated quotient).  The h-index of the
        production selector is a class function too; the toy collapses the six
        orbit slots of a column into one row, which is a coarsening of the row
        set, not of the column law.
        """
        K = (sm64(canon ^ np.uint64(0x243F6A8885A308D3)) >> np.uint64(25))
        K = K.astype(np.int64)                    # uniform on [0, 2^39)
        # exact arithmetic in int64: the production cell() multiplies K*R for a
        # ~381-bit K; here the uniform grid is 2^39 and R <= 2^20, so the product
        # is < 2^59 (asserted, never wrapped) and the shift is the same one.
        assert int(K.max(initial=0)) < (1 << 39)
        assert self.R.bit_length() + 39 < 63, (self.R, "K*R would overflow int64")
        return selector.cell_vec(K, self.R, 1 << 40).astype(np.int64)

    def dp(self, canon, d):
        """Class-function report gate: fingerprint of the canonical class key."""
        if d <= 0:
            return np.ones(np.asarray(canon).shape, dtype=bool)
        thr = np.uint64((1 << 64) >> d if d < 64 else 0)
        return sm64(canon ^ np.uint64(0xD1B54A32D192ED03)) < thr

    # ---- rows ------------------------------------------------------------
    def row(self, canon, alt=False):
        """The orbit row applied at a canonical class: (column, coefficient u)."""
        m = self.col(canon)
        if alt:
            m = (m + 1) % self.R
        if self.row_mode == "broken":
            # MUTATION ARM (tests only): a row not tied to the class's orbit, so
            # the equivariance test must FAIL.  Never used in a census.
            m = (m * 1103515245 + 12345) % self.R
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


def fresh_state(U, x, coins, n):
    """A FRESH lane start drawn from the walker's own coin stream (v20j).

    Mirrors the DP-triggered restart the fleet already performs: draw (a0, c0)
    uniform mod q, form the canonical representative of t = a0 + c0*x, and
    rescale so the invariant a + c*x == canon holds at the fresh state.  The
    draws come from `coins` (a numpy Generator owned by the caller's walker
    loop), so every restart is reproducible and its coin use is counted.
    """
    global FRESH_DRAWS
    a0 = coins.integers(0, U.q, size=n).astype(np.uint64)
    c0 = coins.integers(0, U.q, size=n).astype(np.uint64)
    FRESH_DRAWS += 2 * n
    t = np.array([(int(ai) + int(ci) * int(xi)) % U.q
                  for ai, ci, xi in zip(a0.tolist(), c0.tolist(), np.asarray(x).tolist())],
                 dtype=np.uint64)
    cn, e, _h = U.canon(t)
    return (cn.astype(np.uint64), (e * a0) % U.q, (e * c0) % U.q,
            a0.astype(np.uint64), c0.astype(np.uint64))


def step(U, canon, a, cc, x, d, want_dp=False, prev=None, esc=False,
         force_every=0, tick=0, esc_run=None, kcap=None, restart_out=None,
         coins=None):
    """One production step for every lane.

    prev        : the class one step back (None disables the escape test)
    esc         : arm the two-step-return repair
    force_every : >0 forces the escape branch on lanes (lane + tick) % force_every
                  == 0 (stress/equivariance runs only, never in the census)
    esc_run     : per-lane run of consecutive substitutions, UPDATED IN PLACE.  When
                  supplied with kcap, a lane that would exceed kcap is RESET and
                  flagged in restart_out (priced restart; no exception, no abort).
    Returns (canon, a, cc, dp_hit, prev_out, fired_mask).
    """
    q = U.q
    m = U.col(canon)
    cn, a2, cc2, _e = _apply(U, canon, a, cc, x, m)
    fired = np.zeros(canon.size, dtype=bool)
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
            fired = fire
    if esc_run is not None:
        esc_run[:] = np.where(fired, esc_run + 1, 0)
        if kcap is not None:
            global RESTARTS
            if restart_out is not None:
                restart_out[:] = False
            over = esc_run > kcap
            n_over = int(np.count_nonzero(over))
            if n_over:
                idx = np.flatnonzero(over)
                if COUNTER_ONLY:
                    # CONTROL ARM (v20i behaviour, tests only): the counter is
                    # cleared and RESTARTS is bumped -- the walk state and the
                    # coefficients are NOT touched.  `restart_faithful` in
                    # scripts/kcap_run.py asserts that this arm FAILS the
                    # state-continuity detector.
                    pass
                else:
                    if coins is None:
                        raise AssertionError(
                            "kcap breach must resample: pass coins=<Generator> "
                            "(a counter-only reset is not a restart)")
                    cn3, a4, cc4, _a0, _c0 = fresh_state(U, x[idx], coins, idx.size)
                    cn = cn.copy()
                    a2 = a2.copy()
                    cc2 = cc2.copy()
                    cn[idx] = cn3
                    a2[idx] = a4
                    cc2[idx] = cc4
                esc_run[idx] = 0
                RESTARTS += n_over
                if restart_out is not None:
                    restart_out[idx] = True
    hit = U.dp(cn, d) if want_dp else None
    prev_out = canon.copy()
    if restart_out is not None and restart_out.any():
        # a fresh lane has no predecessor: reset the escape memory too
        prev_out = prev_out.copy()
        prev_out[restart_out] = cn[restart_out]
    return cn, a2.astype(np.uint64), cc2.astype(np.uint64), hit, prev_out, fired


def walk(U, start, x, a0, steps, d, esc=True):
    """Walk `steps` steps; return per-step class trace, DP hit mask, fire mask."""
    lanes = start.size
    canon = U.canon(start)[0]
    a = a0.copy()
    cc = np.zeros(lanes, dtype=np.uint64)
    prev = canon.copy()
    trace = np.empty((steps, lanes), dtype=np.uint64)
    hits = np.zeros((steps, lanes), dtype=bool)
    fires = np.zeros((steps, lanes), dtype=bool)
    for r in range(steps):
        trace[r] = canon
        canon, a, cc, hit, prev, fired = step(U, canon, a, cc, x, d, want_dp=True,
                                              prev=prev, esc=esc, tick=r)
        fires[r] = fired
        hits[r] = hit
    return trace, hits, fires
