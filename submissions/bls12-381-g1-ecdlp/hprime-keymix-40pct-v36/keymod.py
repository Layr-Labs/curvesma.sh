#!/usr/bin/env python3
"""keymod.py v36 -- the recovery kernel AND the lifting stage the v35 verdict said
was missing (ticket 3dec3734, flaws 1).

WHAT CHANGED FROM v35, stated because it is the whole point of this file:

  * the key is the REAL map production uses, x([e]P0) mod 2^b, on a real curve
    (ec.Curve); v35 keyed on low bits of the exponent, which is not the production
    predicate and which is why its "collision" tests proved nothing about x equality.
  * `Pool` is built over the TRUE Cartesian product (bases x walks): every base
    walks the SAME list, and every (i,l) pair is inserted.  v35's census sampled one
    random walk offset per base; that is fixed here and reported per cell.
  * `lift()` implements the recovery from s0 = s mod m to s: with m = (r-1)/2^D the
    map s -> s mod m is 2^D-to-1 on [0,r), so

        s = s0 + k*m,  k in [0, 2^D)

    and the unique k is found by EXHAUSTIVE ENUMERATION with a full-point check
    [s0 + k*m]P0 == T (not a truncated-key check).  The walk is incremental:
    R_k = [s0]P0 + k*[m]P0, one mixed Jacobian add per candidate, and the equality
    test is inversion-free (X == xT*Z^2, Y == yT*Z^3).  A second, independent check
    [s^d]P0 == A_d is applied to the surviving candidate when A_d is supplied.

The truncated key is used ONLY as a filter to choose which records to test (its
false-positive rate is 2^-b per candidate, so over 2^D candidates the expected number
of spurious passes of the FILTER is 2^(D-b), which is charged as its own ledger row);
no equality of secrets is ever concluded from the truncated key.
"""


# ---------------------------------------------------------------- key predicate
class KeyFn:
    """key(e) = x([e]P0) mod 2^bits, on a real curve.  x(P) = x(-P), which is the
    one structural fact the collision step uses, and it is true by construction."""

    def __init__(self, curve, P0, bits, order):
        self.cv, self.P0, self.bits, self.order = curve, P0, bits, order
        self.mask = (1 << bits) - 1
        self._x = {}

    def point(self, e):
        e %= self.order
        if e not in self._x:
            self._x[e] = self.cv.mul(e, self.P0)
        return self._x[e]

    def key(self, e):
        Pt = self.point(e)
        if Pt is None:
            return None
        return Pt[0] & self.mask

    def key_of_point(self, Pt):
        return None if Pt is None else Pt[0] & self.mask


def key(e, m, bits):
    """legacy exponent projection -- kept ONLY so the v35 tests still run against a
    labelled stand-in.  Nothing in v36 charges or recovers with it."""
    return (e % m) & ((1 << bits) - 1)


def encode(c, k, s, m):
    return (c + k + s) % m


def decode(e, c, k, m):
    return (e - c - k) % m


# ---------------------------------------------------------------- the pool
class Pool:
    """One production pool: bases x walks, one challenge c, one shared walk list.

    records[(i,l)] = x([(c + k_i + s_l) mod m]P0) masked to `bits`.
    Built by walking: one scalar multiplication per base, then one point addition
    per walk step -- so the pool really is the Cartesian product, and the cost of
    building it that way is what the ledger's walk rows price.
    """

    def __init__(self, m, bits, c, bases, walks, keyfn=None):
        self.m, self.bits, self.c = m, bits, c
        self.bases, self.walks = list(bases), list(walks)
        self.keyfn = keyfn
        self.keys = {}
        self.records = {}
        if keyfn is not None:
            self._build_ec()
        else:
            for i, k in enumerate(self.bases):
                for l, s in enumerate(self.walks):
                    self.keys.setdefault(key(encode(c, k, s, m), m, bits), []).append((i, l))

    def _build_ec(self):
        cv, P0 = self.keyfn.cv, self.keyfn.P0
        for i, k in enumerate(self.bases):
            A = cv.mul((self.c + k) % self.m, P0)
            if A is None:
                continue
            Bp = cv.mul(self.walks[0] % self.m, P0)
            for l, s in enumerate(self.walks):
                if l:
                    Bp = cv.add(Bp, P0)
                R = cv.add(A, Bp)
                if R is None:
                    continue
                self.records[(i, l)] = R
                self.keys.setdefault(R[0] & ((1 << self.bits) - 1), []).append((i, l))

    # -- collision search against the target key
    def hits_for_target(self, target_pt):
        tk = self.keyfn.key_of_point(target_pt) if self.keyfn else None
        return list(self.keys.get(tk, [])), tk

    def find_collision(self, exclude=None):
        for kk, lst in self.keys.items():
            if len(lst) > 1:
                for a in range(len(lst)):
                    for b in range(a + 1, len(lst)):
                        if exclude and (lst[a], lst[b]) == exclude:
                            continue
                        return lst[a], lst[b]
        return None

    # -- recovery: collision -> s0 = s mod m (the lift seed)
    def recover(self, coll):
        (i, l), (i2, l2) = coll
        e1 = encode(self.c, self.bases[i], self.walks[l], self.m)
        e2 = encode(self.c, self.bases[i2], self.walks[l2], self.m)
        d = (self.bases[i] + self.walks[l]) - (self.bases[i2] + self.walks[l2])
        return dict(e1=e1, e2=e2, s0_mod_m=d % self.m, has_target=(e1 % self.m == (self.bases[i] + self.walks[l]) % self.m))


# ---------------------------------------------------------------- LIFTING
def lift(curve, P0, T, s0, m, D, ad=None, d=None, r=None):
    """Exhaustive lift from s0 = s mod m to s, over k in [0, 2^D).

    Returns dict(passes=[k...], s=..., checks=...).  Exactly one k passes when
    s0 really is s mod m and m*2^D >= r; zero pass otherwise (that is the control).
    """
    xT, yT = T
    R = curve.mul(s0, P0)                      # [s0]P0
    Rj = (R[0], R[1], 1) if R is not None else (0, 1, 0)
    Bm = curve.mul(m, P0)                      # [m]P0
    passes, checks = [], 0
    for k in range(1 << D):
        if curve.jeq_affine(Rj, xT, yT):
            passes.append(k)
            if len(passes) > 4:                # more than a handful = degenerate, stop
                break
        Rj = curve.jadd_mixed(Rj, Bm)
        checks += 1
    out = dict(passes=passes, candidates_checked=checks + 1, D=D, m=m,
               unique=len(passes) == 1, found=bool(passes))
    if passes:
        s = s0 + passes[0] * m
        out["s"] = s
        out["t_check"] = curve.mul(s, P0) == T
        if ad is not None and d is not None:
            assert r is not None, "the A_d check needs the subgroup order"
            out["ad_check"] = curve.mul(pow(s, d, r), P0) == ad
    else:
        out["s"] = None
        out["t_check"] = False
    return out


def cheon_root_check(s, d, r):
    return (r - 1) % d == 0, pow(s, d, r)
