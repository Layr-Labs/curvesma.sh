#!/usr/bin/env python3
"""THE step of the v14 walk -- the single definition imported by the toy, the
evidence runners and (by name) the ledger description.  Nothing restates it.

Design (v14): EQUIVARIANT mu6-CLASS rho, quotient in the REPORT FILTER, no
canonicalisation, no state inversion, and ONE routing guard of depth D that
removes the absorbing all-multiplier cycle family structurally.

  state   y in GT (any representative), labels (a,b) with y = g^a z^b, plus the
          per-trail guard ring of the last D multiply-step (row, column) pairs.
          The state is NOT held canonical: it is whatever the walk reaches.

  key     keys12(y) = {A zeta^j} U {B zeta^j}, j = 0..5 (twelve Fp values).
              mkey  = min(keys12(y))   SIGMA-INVARIANT (class key)
              kappa = index of that minimum, kappa(sigma^j y) = kappa(y) - j

  step    t, dp, sw = routing.fields(mkey)
          sw < 3   (3 of 8 classes)  y <- y^2, (a,b) <- (2a,2b); the guard ring is
                                     not touched (a squaring forms no inverse)
          else     col <- guard.lookup(t) if the row is in the ring else kappa
                   m   <- TBL[t][col], y <- y*m, a <- a + U[t][col], b <- b
                   guard.push(t, col)
          where TBL[t][k] = sigma^{-k}(M_t) = g^{u_t p^{-k}},
                U[t][k]   = u_t p^{-k} mod r  (twelve pre-expanded columns)

  why it merges   for a sigma-partner y_j = sigma^j y: mkey is the same, so t and
                  sw are the same; kappa_j = kappa - j; both the recorded and the
                  looked-up column indices are shifted by -j, so the guard takes
                  the same decision and m_j = sigma^j m; hence y_j*m_j =
                  sigma^j(y*m).  The successor CLASS is a function of the class,
                  for every representative, with no canonicalisation.

  why no trap     v13's family was: same row t, kappa(y') = kappa(y)+6, so the
                  next multiplier was TBL[t][kappa+6] = TBL[t][kappa]^{-1} and
                  y returned exactly.  The guard removes it at the root: in a
                  cycle of length L <= 2D+1 every same-row pair is at
                  cycle-distance <= D, so the second visit reuses the FIRST
                  visit's column and the row contributes L*u_t*p^{-k}, L < r,
                  which is nonzero.  Longer cycles are bounded in argument.md
                  section 4 (>= L/2 row coincidences, rate (1/(12R))^(L/2)).

  report  the distinguished point is tested on mkey, a function of the class
          alone: two trails at a common class report at the same step with the
          same key and different labels.  That is the solvable event.  The guard
          never touches the DP predicate or the record key.

  solve   y1 = sigma^j y2 for exactly one j in Z/12, so
              a1 + b1 x = p^j (a2 + b2 x)  =>  x = (p^j a2 - a1)/(b1 - p^j b2).
          j = 6 is the INVOLUTION (p^6 = -1): x = -(a1+a2)/(b1+b2).  The twelve
          candidates are tried in order and the answer is VERIFIED against z.
"""
import constants as K
import routing


class Guard:
    """The depth-D ring of (row, column) pairs of the last D multiply steps.

    Most recent first.  lookup returns the recorded column of the most recent use
    of the row, or None.  The ring is per trail and 18 bits per slot at
    production width; it is the ONLY state the walk carries besides (y,a,b).
    """
    __slots__ = ('depth', 'slots')

    def __init__(self, depth=None):
        self.depth = K.GUARD_DEPTH if depth is None else depth
        self.slots = [(-1, 0)] * self.depth

    def lookup(self, t):
        for rt, col in self.slots:
            if rt == t:
                return col
        return None

    def push(self, t, col):
        s = self.slots
        for i in range(len(s) - 1, 0, -1):
            s[i] = s[i-1]
        s[0] = (t, col)

    def fires_at(self, t):
        """Did the guard rewrite the column at this step?  (evidence only)"""
        return self.lookup(t) is not None

    def copy(self):
        g = Guard(self.depth)
        g.slots = list(self.slots)
        return g


def keys12(T, y):
    return T.keys12(y)


def classkey(T, y):
    """(sigma-invariant class key, rotation index kappa).  The report filter."""
    ks = T.keys12(y)
    m = min(ks)
    return m, ks.index(m)


# name kept for the ledger/toy vocabulary of the family
kappa = classkey


def step(T, y, a, b, TBL, U, r, table_bits=None, guard=None):
    """One production step.  Returns (y', a', b', kind, info).

    kind is 'sq' (routing squaring channel) or 'mul'.  info = (kappa, row, col,
    guard_fired).  There is no escape branch and no restart: the guard rewrites
    one column and the trail continues.
    """
    mkey, kap = classkey(T, y)
    t, _dp, sw = routing.fields(mkey, table_bits)
    if sw < 3:                                     # 3 of 8 classes: squaring channel
        return T.cyclo_sqr_fast(y), (2 * a) % r, (2 * b) % r, 'sq', (kap, t, kap, 0)
    fired = 0
    if guard is None:
        col = kap
    else:
        col = guard.lookup(t)
        if col is None:
            col = kap
        else:
            fired = 1
        guard.push(t, col)
    u = U[t][col]
    return T.mul(y, TBL[t][col]), (a + u) % r, b, 'mul', (kap, t, col, fired)


def row_diff_mask(mkey, rows):
    """The v24 row-difference mask: 1 + (h mod (rows-1)), h from the SAME mix.

    rows is a power of two (the table size).  The mask is in [1, rows-1], so
    t XOR mask != t for every t in [0, rows), and m -> t XOR m is a bijection on
    [0, rows) for fixed t: the remapped row is uniform over the rows that are NOT
    the previous one.  Row marginal stays exactly uniform.
    """
    h = routing.mix(mkey)
    w = (h >> 128) | (h << 128)                      # bits the routing fields do not use
    return 1 + (w % (rows - 1)) if rows > 2 else 1


def step_v24(T, y, a, b, TBL, U, r, table_bits=None, rows=None, prev_row=-1):
    """THE v24 STEP: guard-free, and the same-row inverse-pair family REMOVED.

    Identical to step_v13 (column = kappa(y), no ring, no history in the class
    map) except that the multiply ROW may not repeat the previous multiply row.
    The exact-return family of v13 needs TWO CONSECUTIVE multiply steps on one
    row: multipliers TBL[t][k] then TBL[t][k+6] = TBL[t][k]^{-1}, so y comes back
    exactly and the trail is absorbing.  With rows forced to differ, the two
    multipliers live in different rows (different u_t), the return equation is
    w_t(L) = 0 for two distinct rows, and the family has NO preimage -- it is
    dead by construction, not debited.

    Returns (y', a', b', kind, info, prev_row).
    """
    mkey, kap = classkey(T, y)
    t, _dp, sw = routing.fields(mkey, table_bits)
    if sw < 3:                                     # squaring: forms no inverse
        return (T.cyclo_sqr_fast(y), (2 * a) % r, (2 * b) % r, 'sq',
                (kap, t, kap, 0), prev_row)
    if rows is not None and t == prev_row:
        t = t ^ row_diff_mask(mkey, rows)           # forced to a DIFFERENT row
    u = U[t][kap]
    return (T.mul(y, TBL[t][kap]), (a + u) % r, b, 'mul',
            (kap, t, kap, 0), t)


def step_v26(T, y, a, b, TBL, U, r, table_bits=None, rows=None):
    """THE v26 STEP: the SQUARED STEP -- family-free by the sigma-orbit theorem.

    y' = (y * m(kappa))^2.  Still class-equivariant: (s^j y * s^j m)^2 = s^j((y m)^2),
    so the successor class depends on the class alone (no prev_row, no ring, no
    history) and a class repeat IS a state repeat, so the cycle abort is sound.

    Why it is the MINIMAL family-free step (evidence/impossibility-theorem.md):
    class merging forces m(k) = sigma^(-k)(A_t), hence m(k+6) = m(k)^(-1) in every
    column -- measured 384/384 on the shipped fixture -- so ANY plain multiply step
    contains the exact-return pair (m, m^-1) at rate 1/(12R).  Squaring breaks it:
    a two-step return needs ((y m1)^2 m2)^2 = sigma^j y, i.e. y^3 = sigma^j/(m1^4 m2^2),
    which has EXACTLY ONE solution y per multiplier pair (cube is a bijection of the
    order-r group), a 2^-219 fraction of the space instead of a 2^-17.6 rate family.

    Returns (y', a', b', kind, info)."""
    mkey, kap = classkey(T, y)
    t, _dp, sw = routing.fields(mkey, table_bits)
    if sw < 3:                                     # squaring: forms no inverse
        return (T.cyclo_sqr_fast(y), (2 * a) % r, (2 * b) % r, 'sq', (kap, t))
    u = U[t][kap]
    return (T.cyclo_sqr_fast(T.mul(y, TBL[t][kap])), (2 * (a + u)) % r, (2 * b) % r,
            'mul', (kap, t))


def step_v13(T, y, a, b, TBL, U, r, table_bits=None):
    """CONTROL: the v13 rule exactly -- column = kappa(y), no ring.  This is the
    rule whose absorbing family the v13 verdict exhibits; T1/T1j must reproduce
    it, and the guard-free ring depth 0 is the same rule."""
    mkey, kap = classkey(T, y)
    t, _dp, sw = routing.fields(mkey, table_bits)
    if sw < 3:
        return T.cyclo_sqr_fast(y), (2 * a) % r, (2 * b) % r, 'sq', (kap, t)
    u = U[t][kap]
    return T.mul(y, TBL[t][kap]), (a + u) % r, b, 'mul', (kap, t)


def solve_candidates(p1, p2, a1, b1, a2, b2, r, order=None):
    """The twelve x candidates for one class collision, j = 0..11.

    j is the sigma-shift with y1 = sigma^j y2.  j = 6 is the involution case
    (the +-pair / label-bit case) and is solved by the same formula.
    """
    out = []
    js = range(12) if order is None else order
    for j in js:
        pj = p1[j] if isinstance(p1, (list, tuple)) else pow(p1, j, r)
        den = (b1 - pj * b2) % r
        if den == 0:
            continue
        out.append((j, (pj * a2 - a1) % r * pow(den, -1, r) % r))
    return out


def label_bit_solve(a1, b1, a2, b2, r):
    """The +-pair case alone (j = 6): y1 = y2^{-1}  =>  x = -(a1+a2)/(b1+b2)."""
    den = (b1 + b2) % r
    if den == 0:
        return None
    return (-(a1 + a2)) % r * pow(den, -1, r) % r


# ---------------------------------------------------------------------------
# CONTROL RULES -- deliberately broken predecessors.  Each property test runs on
# them too, so the test can be seen to FAIL on a mechanism that has the defect.
# ---------------------------------------------------------------------------

def kappa_raw(T, y):
    """v12c's key: c_0(y)^(1+p) on the RAW representative -- NOT sigma-invariant."""
    return T.kappa_key(y)


def step_classkeyed(T, y, a, b, tbl, us, r, table_bits=None):
    """v12c control: class-keyed multiplier on the raw representative."""
    key = T.kappa_key(y)
    t, _dp, sw = routing.fields(key, table_bits)
    if sw < 3:
        return T.cyclo_sqr_fast(y), (2 * a) % r, (2 * b) % r, 'sq', (0, t)
    return T.mul(y, tbl[t]), (a + us[t]) % r, b, 'mul', (0, t)


def step_v10f(T, y, a, b, tbl, us, r, frob_pow, ppow, kcache=None, table_bits=None):
    """v12b control: canonicalising step with NO escape."""
    _m, k = classkey(T, y)          # v12b routed on the ROTATION INDEX
    t, _dp, sw = routing.fields(k, table_bits)
    if sw < 3:
        raw = T.cyclo_sqr_fast(y)
        a2, b2 = (2 * a) % r, (2 * b) % r
    else:
        raw = T.mul(y, tbl[t])
        a2, b2 = (a + us[t]) % r, b
    _m2, kk = classkey(T, raw)
    y2 = frob_pow(raw, kk)
    pk = ppow[kk]
    return y2, a2 * pk % r, b2 * pk % r
