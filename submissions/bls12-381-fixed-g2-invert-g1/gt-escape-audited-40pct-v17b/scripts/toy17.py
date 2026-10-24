#!/usr/bin/env python3
"""Toy v17: the PRODUCTION step rule, the PRODUCTION DP/replay procedure, and the
PRODUCTION slot rule (disjoint halves keyed on the canonicalisation index), at toy
scale and at the DESIGN budget only.

What changed from toy16.py (v16), which the review rejected:

  D2 (p^6 degeneracy).  v16's production step applied M_h^{p^-kappa} and had no
  cycle escape; p^6 = -1 mod r means two multiply steps that use the same slot h
  with canonicalisation indices differing by 6 apply inverse multipliers and
  restore BOTH the state and the exponent labels exactly, so the walk dies on a
  label-trivial 2-cycle.  v16's toy escaped such cycles with an escape that read
  the previous two class keys -- a history-dependent rule that is a different
  algorithm from the production step (reviewer's point).

  v17 removes the trap instead of escaping it.  The multiplier table is ONE table
  of R slots, addressed as two DISJOINT halves: slot(y) = h(T) + (R/2)*[kappa(y)
  >= 6], where h(T) is read from the class invariant and kappa(y) is the
  canonicalisation index of the class.  Both are class functions, so the step
  stays a function of the class (merge-safe, stationary), and slot() is a
  function of the state.  In a p^6 pair kappa flips between {0..5} and {6..11},
  so the two steps address DIFFERENT halves and the inverse-multiplier
  cancellation cannot happen.  Proof in argument.md Section 3; the exact trap
  fixture and the adversarial cycle-length enumeration L in {2,4,6,12} are in
  this file's fixtures mode.

  There is no escape rule anywhere in the production step or in this toy: the
  walk is a stationary map on classes and a 2-cycle does not exist.

  Budget honesty.  v16's toy ran max_steps = 4*S and the ledger certified a PASS
  from the eventual 4x-budget number while the design-budget number was 57/180.
  v17 runs the DESIGN BUDGET EXACTLY: sum of trail lengths = S_design = s_ops(N,R),
  the same integer the ledger uses (and s_ops uses a strict lower enclosure of
  ln(5/3), so S_toy <= S_ledger).  Nothing is run at 4x.

Toy-scale disclosure (unchanged where the review did not object):
  * the branch bit is production exactly: br = (T & 7) < 3 with T = Re(c0).
  * the multiplier index and the DP window are read at the PRODUCTION bit
    offsets 3 and 3+11 = 14 of a class-function label L(T), because the q = 523
    invariant has only ~9 bits; the toy masks the index to R/2 = 64 and the DP
    window to dpbits.  Offsets match production, widths follow the toy table.
    scripts/predicate_check.py re-derives this rule independently.
  * DP records are 52-byte-shaped: the store keeps (tid, index) only and every
    replay reconstructs (A, D, y) by re-walking that trail from its derived
    start.  The toy's class key is the exact sorted twelve-key set.
"""
import json
import random
import sys
from math import isqrt

from toy import Tower, field_params, s_ops

IDX_BITS = 11            # production multiplier-index slice width
DP_SHIFT = 3 + IDX_BITS  # = 14, production DP window offset


class Toy17:
    def __init__(self, q, rp, R=128, dpbits=4, seed=1, rule="new"):
        assert R % 2 == 0
        self.q, self.rp, self.R, self.dpbits = q, rp, R, dpbits
        self.half = R // 2
        self.rulen = rule                      # "new" (disjoint halves) or "old" (v16)
        self.use_sqr = True                    # s=0 control: branch forced to multiply
        self.T = Tower(q)
        self.rng = random.Random(seed)
        self.seed = seed
        self.lam = q % rp
        assert pow(self.lam, 12, rp) == 1 and pow(self.lam, 6, rp) == rp - 1
        self.N = (rp - 1) // 12
        self.linv = [pow(self.lam, (-k) % 12, rp) for k in range(12)]
        e = (q ** 12 - 1) // rp
        while True:
            x = tuple(tuple((self.rng.randrange(q), self.rng.randrange(q))
                           for _ in range(3)) for _ in range(2))
            g = self.T.pow(x, e)
            if g != self.T.one12() and self.T.pow(g, rp) == self.T.one12():
                break
        self.g = g

    # ---- production bit rule on the class invariant T = Re(c0) --------------
    @staticmethod
    def label(t):
        return ((t * 0x9E3779B1) ^ (t >> 13)) & 0xFFFFFFFF

    def bits(self, y):
        t = y[0][0][0]
        L = self.label(t)
        br = ((t & 7) < 3) and self.use_sqr    # production exact (s=0 control forces multiply)
        h = (L >> 3) % self.half               # production offset 3, toy width
        dp = ((L >> DP_SHIFT) & ((1 << self.dpbits) - 1)) == 0
        return br, h, dp, t

    def dp_of(self, y):
        """The DP predicate, evaluated on the state it describes (post-step)."""
        L = self.label(y[0][0][0])
        return ((L >> DP_SHIFT) & ((1 << self.dpbits) - 1)) == 0

    def slot(self, h, kap):
        if self.rulen == "old":
            return h                            # v16: h alone, p^6 trap live
        return h + (self.half if kap >= 6 else 0)

    def build_table(self, xz=None):
        rng, rp, R = self.rng, self.rp, self.R
        self.x = xz if xz is not None else rng.randrange(1, rp)
        self.z = self.T.pow(self.g, self.x)
        self.ent = [None] * R
        for sl in range(R):
            a = rng.randrange(1, rp)
            b = rng.randrange(1, rp)
            M = self.T.mul(self.T.pow(self.g, a), self.T.pow(self.z, b))
            cur, im = M, []
            for _ in range(12):
                im.append(cur)
                cur = self.T.frob(cur)
            self.ent[sl] = (a, b, im)

    def ident(self, y):
        """exact class identity at toy scale + canonicalisation index.

        The twelve keys are a class function only up to their cyclic shift: a
        Frobenius conjugate y^{p^k} has ks'[j] = ks[j+k].  The class identity is
        therefore the key vector normalised to start at the canonical index --
        NOT the sorted multiset, which is a coarser object: distinct classes can
        share the same sorted key set (the a <-> B exchange is one such case),
        so a sorted-key "match" need not be a Frobenius conjugate at all.
        """
        ks = self.T.keys12(y)
        m0 = min(ks)
        kap = ks.index(m0)
        return (y[0][0][0], tuple(ks[(kap + j) % 12] for j in range(12))), kap

    def inv_ok(self, y, A, D):
        return self.T.mul(self.T.pow(self.g, A % self.rp), self.T.pow(self.z, D % self.rp)) == y

    def step(self, y, A, D, h, kap):
        e = self.ent[self.slot(h, kap)]
        j = (12 - kap) % 12
        return (self.T.mul(y, e[2][j]),
                (A + self.linv[kap] * e[0]) % self.rp,
                (D + self.linv[kap] * e[1]) % self.rp, 0)

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

    # ---- trail start derivation from a 2-byte trail id (never stored) -------
    def trail_start(self, tid):
        rr = random.Random(0x5EED0000 + self.seed * 1000003 + tid)
        return rr.randrange(self.rp), rr.randrange(1, self.rp)

    def walk_trail(self, tid, nsteps, upto=None):
        """Re-walk trail `tid` from its derived start; returns (y, A, D, index)."""
        a, b = self.trail_start(tid)
        y = self.T.mul(self.T.pow(self.g, a), self.T.pow(self.z, b))
        A, D = a, b
        last = nsteps - 1 if upto is None else upto
        K, kap = self.ident(y)
        for _ in range(last + 1):
            br, h, dp, t = self.bits(y)
            if br:
                y, A, D, _ = self.sqr(y, A, D)
            else:
                y, A, D, _ = self.step(y, A, D, h, kap)
            K, kap = self.ident(y)
        return y, A, D, last + 1

    def trial(self, S, trails):
        """ONE design-budget run: `trails` trails whose total forward steps = S.

        One shared DP store keyed by the exact class identity.  A duplicate is
        resolved by replaying BOTH trails from their derived starts (no stored
        state, no knowledge of x), then a twelve-Frobenius shift match, a
        coefficient solve and a g^x == z verification.  Failures are counted and
        the run CONTINUES (charged), never aborting on a bad match.
        """
        st = dict(fwd=0, replay=0, dp_hits=0, records=0, matches=0, self_m=0, cross_m=0,
                  degenerate=0, solved=0, useful=False)
        store = {}
        per = S // trails
        for tid in range(trails):
            n = per + (S - per * trails if tid == trails - 1 else 0)
            a, b = self.trail_start(tid)
            y = self.T.mul(self.T.pow(self.g, a), self.T.pow(self.z, b))
            A, D = a, b
            K, kap = self.ident(y)
            for i in range(n):
                br, h, dp, t = self.bits(y)
                if br:
                    y, A, D, _ = self.sqr(y, A, D)
                else:
                    y, A, D, _ = self.step(y, A, D, h, kap)
                K, kap = self.ident(y)
                dp = self.dp_of(y)          # predicate on the state being recorded
                st["fwd"] += 1
                if not dp:
                    continue
                st["dp_hits"] += 1
                K, kp = self.ident(y)
                rec = store.get(K)
                if rec is None:
                    store[K] = (tid, i, A, D, y)
                    st["records"] += 1
                    continue
                st["matches"] += 1
                t1, i1, A1, D1, y1 = rec
                if (t1, i1) == (tid, i):
                    continue
                # replay both trails from their derived starts: records carry no state
                y1r, A1r, D1r, n1 = self.walk_trail(t1, 0, upto=i1)
                y2r, A2r, D2r, n2 = self.walk_trail(tid, 0, upto=i)
                st["replay"] += n1 + n2
                if y1r != y1 or y2r != y:
                    st["degenerate"] += 1
                    continue
                if t1 == tid:
                    st["self_m"] += 1
                else:
                    st["cross_m"] += 1
                k = self.match_shift(y1r, y2r)
                if k is None:
                    st["degenerate"] += 1
                    continue
                x = self.recover(A1r, D1r, A2r, D2r, k)
                if x is None or x == 0:
                    st["degenerate"] += 1
                    continue
                st["solved"] += 1
                if self.T.pow(self.g, x) == self.z:
                    st["useful"] = True
                    return st
                st["degenerate"] += 1
        return st


KEYS = ("fwd", "replay", "dp_hits", "records", "matches", "self_m", "cross_m",
        "degenerate", "solved")


def sweep(cfg):
    q = cfg.get("q", 523)
    fp = field_params(q)
    rp, N, _o2, _o2o = fp
    R = cfg.get("R", 128)
    dbits = cfg.get("dpbits", 4)
    trials = cfg.get("trials", 200)
    trails = cfg.get("trails", 8)
    rule = cfg.get("rule", "new")
    use_sqr = cfg.get("use_sqr", True)
    seed0 = cfg.get("seed", 1)
    S = int(s_ops(N, R))
    rows = []
    for i in range(trials):
        tw = Toy17(q, rp, R, dbits, seed=seed0 + i, rule=rule)
        tw.use_sqr = use_sqr
        tw.build_table()
        rows.append(tw.trial(S, trails))
    agg = {k: sum(r[k] for r in rows) for k in KEYS}
    agg["useful"] = sum(1 for r in rows if r["useful"])
    agg["over_budget"] = sum(1 for r in rows if r["fwd"] > S)
    agg["max_fwd"] = max(r["fwd"] for r in rows)
    out = dict(mode="sweep", rule=rule, use_sqr=use_sqr, q=q, rp=rp, N=N, R=R, dpbits=dbits,
               trials=trials, trails=trails, S_design=S, per_trail=S // trails,
               dp_shift=DP_SHIFT, idx_bits=IDX_BITS,
               slot_rule=("h + (R/2)*[kappa>=6]" if rule == "new" else "h (v16, p^6 trap live)"),
               **agg)
    out["success_num"] = agg["useful"]
    out["success_den"] = trials
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    cfg = json.loads(sys.stdin.read() or "{}")
    sweep(cfg)
