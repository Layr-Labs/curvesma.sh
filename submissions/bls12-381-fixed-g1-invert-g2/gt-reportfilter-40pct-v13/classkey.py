#!/usr/bin/env python3
"""THE step of the v13 walk -- the single definition imported by the toy, the
evidence runners and (by name) the ledger description.  Nothing restates it.

Design (v13): EQUIVARIANT mu6-CLASS rho.  The quotient lives in the REPORT
FILTER; the transition never canonicalises and never inverts.

  state   y in GT (any representative), labels (a,b) with y = g^a z^b.  The
          state is NOT held canonical: it is whatever the walk reaches.

  key     keys12(y) = {A zeta^j} U {B zeta^j}, j = 0..5 (twelve Fp values).
          The SET is a function of the Frobenius class: keys12(sigma^j y) is
          the same set, with the index of its minimum shifted by -j.  So

              mkey = min(keys12(y))          is SIGMA-INVARIANT (the class key)
              kappa = index of that minimum  satisfies kappa(sigma^j y) = kappa(y) - j

  step    t, dp, sw = routing.fields(mkey)
          sw < 3   (3 of 8 classes)  y <- y^2,          (a,b) <- (2a, 2b)
          else                       m <- TBL[t][kappa], y <- y*m,
                                     a <- a + U[t][kappa], b <- b
          where TBL[t][k] = sigma^{-k}(M_t) = g^{u_t p^{-k}} and
                U[t][k]   = u_t p^{-k} mod r     (twelve pre-expanded columns)

  why it merges   for a sigma-partner y_j = sigma^j y: mkey is the same, so t and
                  sw are the same; kappa_j = kappa - j, so
                  m_j = sigma^{-kappa+j}M_t = sigma^j m and
                  y_j * m_j = sigma^j (y*m).  The successor CLASS is a function of
                  the class, for EVERY representative, with no canonicalisation.
                  Squaring is sigma-equivariant on the nose.  T2 checks the
                  successor ELEMENTS, not just their keys, and the mutation
                  controls (v12c class-keyed, v12b canonicalising) both split.

  why no trap     no branch applies sigma^6 to a state, because no branch applies
                  sigma to a state at all.  The v12b absorbing cycle needs the
                  step to ADOPT an inverted product; here the inversion is never
                  formed.  The only state 2-cycle available is m*m' = 1 inside one
                  class, i.e. u_t p^{-k} + u_t' p^{-k'} = 0 mod r, a coincidence
                  with the table size, and the exhaustive fixture census (T1) and
                  the class-level census (T4-cycle mass) measure it.

  report  the distinguished point is tested on mkey, which is a function of the
          class alone: two trails at a common class report at the same step with
          the same key and different labels.  That is the solvable event.

  solve   y1 = sigma^j y2 for exactly one j in Z/12, so
              a1 + b1 x = p^j (a2 + b2 x)  =>  x = (p^j a2 - a1)/(b1 - p^j b2).
          j = 6 is the INVOLUTION (p^6 = -1): x = -(a1+a2)/(b1+b2).  The twelve
          candidates are tried in order and the answer is VERIFIED against z, so
          a wrong sign costs a verification, never a wrong output.
"""

import routing


def keys12(T, y):
    return T.keys12(y)


def classkey(T, y):
    """(sigma-invariant class key, rotation index kappa).  The report filter."""
    ks = T.keys12(y)
    m = min(ks)
    return m, ks.index(m)


# name kept for the ledger/toy vocabulary of the family
kappa = classkey


def step(T, y, a, b, TBL, U, r, table_bits=None):
    """One production step.  Returns (y', a', b', kind, info).

    kind is 'sq' (routing squaring channel) or 'mul'.  There is no escape branch:
    the transition has no canonicalisation to escape from.
    """
    mkey, kap = classkey(T, y)
    t, _dp, sw = routing.fields(mkey, table_bits)
    if sw < 3:                                     # 3 of 8 classes: squaring channel
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
    """v12c control: class-keyed multiplier on the raw representative, no sigma
    column, no canonicalisation.  Merge-unsafe by construction (the key of the
    successor is not a function of the class of the input); the mutation control
    for T2/T3."""
    key = T.kappa_key(y)
    t, _dp, sw = routing.fields(key, table_bits)
    if sw < 3:
        return T.cyclo_sqr_fast(y), (2 * a) % r, (2 * b) % r, 'sq', (0, t)
    return T.mul(y, tbl[t]), (a + us[t]) % r, b, 'mul', (0, t)


def step_v10f(T, y, a, b, tbl, us, r, frob_pow, ppow, kcache=None, table_bits=None):
    """v12b control: canonicalising step with NO escape -- the rule whose
    absorbing 2-cycle the v12b verdict exhibits and T1 must reproduce."""
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
