#!/usr/bin/env python3
"""Toy v15: the PRODUCTION step rule and the PRODUCTION DP procedure at toy scale.

Differences from scripts/toy.py (v14), which the reviewer flagged:
  * branch/multiplier/DP bits are read from the class invariant T = Re(c0) with the
    production bit layout:  br = (T & 7) < 3,  h = (T >> 3) % R,
    dp = ((T >> (3 + log2 R)) & (2^dpbits - 1)) == 0.   v14 selected the branch from
    the table index h, which is a different function.
  * the store holds DP records only (one insert per DP step), not every visited
    state, exactly as the production spec: fingerprint -> (A, D, y).  A duplicate
    is replayed (twelve-Frobenius match, coefficient solve, g^x == z verification);
    a degenerate or failing replay is discarded, never counted as a collision.
  * a repeated class that is NOT a DP match is a structural 2-cycle: the trail
    escapes with one forced multiply step at index (h+1) % R, which is charged at
    the multiply price and adds no step, so escapes cost no extra budget.
  * every DP match is classified into the two collision channels:
        k == 0  same-index channel   (identical representative)
        k != 0  shifted-image channel (y1 = y2^(p^k)), the gap the reviewer named
    Both are counted and both are tested for usability, so the record shows which
    channel the toy's merges actually come from.
Toy-scale disclosure: the toy key space is one base-field element, so the exact
class identity is the sorted twelve-key set (v14 docstring); the production
fingerprint is 256 bits of the canonical minimum key.
"""
import json, random, sys
from math import isqrt

from toy import Tower, field_params, s_ops


class Toy15:
    def __init__(self, q, rp, R=128, dpbits=4, seed=1):
        self.q = q
        self.rp = rp
        self.R = R
        self.dpbits = dpbits
        self.nmul = (R - 1).bit_length()          # multiplier bits, production layout
        self.dpshift = 3 + self.nmul
        self.T = Tower(q)
        self.rng = random.Random(seed)
        self.lam = q % rp
        assert pow(self.lam, 12, rp) == 1 and pow(self.lam, 6, rp) == rp - 1
        self.N = (rp - 1) // 12
        self.linv = [pow(self.lam, (-k) % 12, rp) for k in range(12)]
        e = (q ** 12 - 1) // rp
        while True:
            x = tuple(tuple((self.rng.randrange(q), self.rng.randrange(q)) for _ in range(3)) for _ in range(2))
            g = self.T.pow(x, e)
            if g != self.T.one12() and self.T.pow(g, rp) == self.T.one12():
                break
        self.g = g

    # ---- production bit rule on the class invariant T = Re(c0) --------------
    # The BRANCH rule is production exactly: br = (T & 7) < 3.  At q=523 the
    # invariant has only 10 bits, so the multiplier index and the DP bits are read
    # from a fixed 32-bit hash L(T) instead of disjoint slices of a 381-bit T;
    # L is still a class function, which is the property the procedure uses.
    @staticmethod
    def label(t):
        return ((t * 0x9E3779B1) ^ (t >> 13)) & 0xFFFFFFFF

    def bits(self, y):
        t = y[0][0][0]
        L = self.label(t)
        br = (t & 7) < 3
        h = L % self.R
        dp = ((L >> 8) & ((1 << self.dpbits) - 1)) == 0
        return br, h, dp, t

    def build_table(self):
        rng, q, rp, R = self.rng, self.q, self.rp, self.R
        self.x = rng.randrange(1, rp)
        self.z = self.T.pow(self.g, self.x)
        self.ent = [None] * R
        for h in range(R):
            a = rng.randrange(1, rp)
            b = rng.randrange(1, rp)
            M = self.T.mul(self.T.pow(self.g, a), self.T.pow(self.z, b))
            cur, im = M, []
            for _ in range(12):
                im.append(cur)
                cur = self.T.frob(cur)
            self.ent[h] = (a, b, im)

    def ident(self, y):
        ks = self.T.keys12(y)
        m0 = min(ks)
        return tuple(sorted(ks)), m0, ks.index(m0)

    def inv_ok(self, y, A, D):
        return self.T.mul(self.T.pow(self.g, A % self.rp), self.T.pow(self.z, D % self.rp)) == y

    def step(self, y, A, D, h, kap, force_idx=None):
        T, rp = self.T, self.rp
        e = self.ent[h]
        j = (12 - kap) % 12
        y2 = T.mul(y, e[2][j])
        return y2, (A + self.linv[kap] * e[0]) % rp, (D + self.linv[kap] * e[1]) % rp, 0

    def sqr(self, y, A, D):
        return self.T.cyclo_sqr(y), (2 * A) % self.rp, (2 * D) % self.rp, 1

    def match_shift(self, y1, y2):
        cur = y1
        for k in range(12):
            if cur == y2:
                return k
            cur = self.T.frob(cur)
        return None

    def recover(self, A1, D1, A2, D2, k):
        rp = self.rp
        pj = pow(self.q, k, rp)
        den = (D2 - pj * D1) % rp
        if den == 0:
            return None
        return ((pj * A1 - A2) % rp) * pow(den, -1, rp) % rp

    def walk(self, max_steps, check_every=16):
        T, rp = self.T, self.rp
        A0 = self.rng.randrange(rp)
        D0 = self.rng.randrange(1, rp)
        y = T.mul(T.pow(self.g, A0), T.pow(self.z, D0))
        A, D = A0, D0
        assert self.inv_ok(y, A, D)
        K, m0, kap = self.ident(y)
        store = {}                                  # DP records only
        st = dict(steps=0, square=0, escapes=0, checks=0, bad=0, dp_hits=0,
                  dp_records=0, dp_matches=0, ch_same=0, ch_shift=0,
                  rec_same=0, rec_shift=0, degenerate=0, useful=False)
        prev = prev2 = None
        n = 0
        while n < max_steps:
            br, h, dp, t = self.bits(y)
            if prev is not None and K == prev or (prev2 is not None and K == prev2):
                # structural 2-cycle: forced multiply at a perturbed index, same
                # price as the multiply step it replaces, no extra step
                h2 = (h + 1) % self.R
                y, A, D, issq = self.step(y, A, D, h2, kap)
                st["escapes"] += 1
            elif br:
                y, A, D, issq = self.sqr(y, A, D)
            else:
                y, A, D, issq = self.step(y, A, D, h, kap)
            n += 1
            st["steps"] = n
            st["square"] += issq
            if n % check_every == 0:
                st["checks"] += 1
                if not self.inv_ok(y, A, D):
                    st["bad"] += 1
            K2, m02, kap = self.ident(y)
            if dp:
                st["dp_hits"] += 1
                rec = store.get(K2)
                if rec is None:
                    store[K2] = (A, D, y)
                    st["dp_records"] += 1
                else:
                    st["dp_matches"] += 1
                    k = self.match_shift(rec[2], y)
                    if k is None:
                        st["degenerate"] += 1
                    else:
                        if k == 0:
                            st["ch_same"] += 1
                        else:
                            st["ch_shift"] += 1
                        x = self.recover(rec[0], rec[1], A, D, k)
                        if x is not None and x != 0 and T.pow(self.g, x) == self.z:
                            if not self.inv_ok(y, A, D):
                                st["bad"] += 1
                            st["checks"] += 1
                            if k == 0:
                                st["rec_same"] += 1
                            else:
                                st["rec_shift"] += 1
                            st["useful"] = True
                            return st
                        st["degenerate"] += 1
                prev2, prev = prev, K2
                continue
            prev2, prev = prev, K2
        return st


if __name__ == "__main__":
    cfg = json.loads(sys.stdin.read() or "{}")
    q = cfg.get("q", 523)
    fp = field_params(q)
    if fp is None:
        raise SystemExit("bad q")
    rp, N, _o2, _o2o = fp
    R = cfg.get("R", 128)
    dbits = cfg.get("dpbits", 4)
    trials = cfg.get("trials", 180)
    S = s_ops(N, R)
    rows = []
    for i in range(trials):
        tw = Toy15(q, rp, R, dbits, seed=cfg.get("seed", 1) + i)
        tw.build_table()
        rows.append(tw.walk(max_steps=4 * S))
    keys = ("steps", "square", "escapes", "checks", "bad", "dp_hits", "dp_records",
            "dp_matches", "ch_same", "ch_shift", "rec_same", "rec_shift", "degenerate")
    agg = {k: sum(r[k] for r in rows) for k in keys}
    agg["useful"] = sum(1 for r in rows if r["useful"])
    agg["hits_at_S"] = sum(1 for r in rows if r["useful"] and r["steps"] <= S)
    out = dict(mode="v15", q=q, rp=rp, N=N, R=R, dpbits=dbits, trials=trials,
               S_ops=int(S), dp_shift=3 + (R - 1).bit_length(),
               branch_rule="(T & 7) < 3 with T = Re(c0); h and DP bits from L(T)", **agg)
    out["square_frac_num"] = agg["square"]
    out["square_frac_den"] = max(1, agg["steps"])
    print(json.dumps(out, indent=1))
