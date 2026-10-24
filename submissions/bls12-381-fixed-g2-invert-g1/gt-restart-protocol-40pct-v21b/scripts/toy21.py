#!/usr/bin/env python3
"""toy21.py -- multi-target census under the RESTARTED-WORKER protocol and the
depth-d slot-history guard.  Reuses the v20 tower/step source (one implementation).

Modes:
  census : >=64 independent (target, table-seed) instances; per-instance success;
           pooled exact Clopper-Pearson lower end via the VERBATIM Decimal helper
           from secp256k1-ecdlp/sigma-table-40pct-v15k/scripts/gen_cert_v15k.py
  guard  : T3 coalescence with fully independent prev-slot history (depth d), plus a
           mutation control that must be able to collapse.
"""
import json
import math
import os
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction as F

getcontext().prec = 60
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toy import Tower, field_params          # noqa: E402  (SELF-CONTAINED: scripts/toy.py shipped)


# ---- verbatim Clopper-Pearson, copied from sigma-table-40pct-v15k/gen_cert_v15k.py
def cp95_lower(k, n, alpha=Decimal('0.05')):
    k = int(k); n = int(n)
    lnC = sum(Decimal(i).ln() for i in range(n - k + 1, n + 1)) - sum(Decimal(i).ln() for i in range(1, k + 1))
    def tail(p):
        p = Decimal(p); q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal('1e-40') * s: break
        return s
    lo, hi = Decimal('0.01'), Decimal('0.99')
    for _ in range(90):
        mid = (lo + hi) / 2
        if tail(mid) > alpha: hi = mid
        else: lo = mid
        if hi - lo < Decimal('1e-13'): break
    return lo


IDX_BITS = 12
DP_SHIFT = 3 + IDX_BITS
LAB_MASK = (1 << 128) - 1
LAB_MUL = 0x9E3779B97F4A7C15F39CC0605CEDC835
SQ_NUM, SQ_DEN = 3, 8


def mix128(x):
    return ((x * LAB_MUL) ^ (x >> 61)) & LAB_MASK


class T21:
    """v20 step + the depth-d slot-history guard (class-consistent), multi-instance."""

    def __init__(self, q, rp, R, dpbits, seed, depth):
        self.q, self.rp, self.R, self.dpbits = q, rp, R, dpbits
        self.depth = depth
        self.T = Tower(q)
        self.rng = random.Random(0xC0FFEE ^ (seed * 2654435761))
        self.lam = q % rp
        assert pow(self.lam, 12, rp) == 1 and pow(self.lam, 6, rp) == rp - 1
        self.N = (rp - 1) // 12
        self.linv = [pow(self.lam, (-k) % 12, rp) for k in range(12)]

    def new_instance(self, seed):
        """Fresh target z, fresh table, fresh g -- one independent challenge+table."""
        rng = random.Random(0x9E37 ^ (seed * 2246822519))
        e = (self.q ** 12 - 1) // self.rp
        while True:
            x = tuple(tuple((rng.randrange(self.q), rng.randrange(self.q))
                            for _ in range(3)) for _ in range(2))
            g = self.T.pow(x, e)
            if g != self.T.one12() and self.T.pow(g, self.rp) == self.T.one12():
                break
        self.g = g
        self.irng = rng
        self.x = rng.randrange(1, self.rp)
        self.z = self.T.pow(self.g, self.x)
        self.ent = []
        for _ in range(self.R):
            a, b = rng.randrange(1, self.rp), rng.randrange(1, self.rp)
            self.ent.append((a, b, self.T.pow(self.g, (a + self.x * b) % self.rp)))

    def keys(self, y):
        ks = self.T.keys12(y)
        best = None
        for i in range(12):
            cand = tuple(ks[(i + t) % 12] for t in range(12))
            if best is None or cand < best[0]:
                best = (cand, i)
        return ks, ks[best[1]], best[1]

    def bits(self, m0):
        L = ((m0 * LAB_MUL) ^ (m0 >> 61)) & LAB_MASK
        return (L & 7) < SQ_NUM, (L >> 3) % self.R, \
               ((L >> DP_SHIFT) & ((1 << self.dpbits) - 1)) == 0, L

    def decide(self, m0, hist):
        br, h, dp, L = self.bits(m0)
        if not br and h in hist:
            return True, -1, dp, L, True
        return br, (-1 if br else h), dp, L, False

    def _rot(self, M, j):
        cur = M
        for _ in range(j % 12):
            cur = self.T.frob(cur)
        return cur

    def fwd(self, y, A, D, kap, m0, hist):
        br, h, dp, L, g = self.decide(m0, hist)
        if br:
            y, A, D = self.T.cyclo_sqr(y), (2 * A) % self.rp, (2 * D) % self.rp
        else:
            a, b, M = self.ent[h]
            y = self.T.mul(y, self._rot(M, (12 - kap) % 12))
            A, D = (A + self.linv[kap] * a) % self.rp, (D + self.linv[kap] * b) % self.rp
            hist = ([h] + hist)[:self.depth]
        _ks, m0, kap = self.keys(y)
        return y, A, D, kap, m0, br, hist, dp

    def class_id(self, y):
        best = None
        cur = y
        for _ in range(12):
            ser = (cur[0][0], cur[0][1], cur[0][2], cur[1][0], cur[1][1], cur[1][2])
            if best is None or ser < best: best = ser
            cur = self.T.frob(cur)
        return best

    def conj_equal(self, y1, y2):
        cur = y1
        for _ in range(12):
            if cur == y2: return True
            cur = self.T.frob(cur)
        return False

    def match_shift(self, y1, y2):
        cur = y1
        for k in range(12):
            if cur == y2: return k
            cur = self.T.frob(cur)
        return None

    def start(self, a, b):
        return self.T.pow(self.g, (a + self.x * b) % self.rp), a % self.rp, b % self.rp

    def recover(self, A1, D1, A2, D2, k):
        pj = pow(self.q, k, self.rp)
        den = (D2 - pj * D1) % self.rp
        if den == 0: return None
        return ((pj * A1 - A2) % self.rp) * pow(den, -1, self.rp) % self.rp

    def ok(self, kk):
        return kk not in (None, 0) and self.T.pow(self.g, kk) == self.z


def replay(tw, start, steps, depth):
    a, b = start
    y, A, D = tw.start(a, b)
    _ks, m0, kap = tw.keys(y)
    hist = []
    for _ in range(steps):
        y, A, D, kap, m0, br, hist, dp = tw.fwd(y, A, D, kap, m0, hist)
    return y, A, D


def mode_census(cfg):
    """RESTARTED-WORKER protocol: each trail ends at its FIRST report, one record,
    restart from a fresh seed-derived start; success counted PER INSTANCE."""
    q, rp, N = cfg.get("q", 523), cfg["rp"], cfg["N"]
    R, dpbits, depth = cfg.get("R", 128), cfg.get("dpbits", 3), cfg.get("depth", 2)
    trails, trials, ninst = cfg.get("trails", 64), cfg.get("trials", 16), cfg.get("instances", 64)
    max_trail = cfg.get("max_trail", 64)
    tw = T21(q, rp, R, dpbits, 1, depth)
    per_inst = []
    tot_ok = tot_tr = 0
    fwd = rep = colls = degen = compares = recs = 0
    for inst in range(ninst):
        tw.new_instance(1000 + inst)
        rng = random.Random(0x5A21 ^ (inst * 7919))
        ok = 0
        for _t in range(trials):
            starts = [(rng.randrange(1, rp), rng.randrange(1, rp)) for _ in range(trails)]
            pool = {}
            for tid in range(trails):
                y, A, D = tw.start(*starts[tid])
                _ks, m0, kap = tw.keys(y)
                hist = []
                rec = None
                k = 0
                while True:
                    y, A, D, kap, m0, br, hist, dp = tw.fwd(y, A, D, kap, m0, hist)
                    k += 1; fwd += 1
                    if dp:
                        rec = (tid, k); recs += 1; break
                    if k >= max_trail: break
                if rec is None: continue
                fp = tw.class_id(y)
                o = pool.get(fp)
                if o is None:
                    pool[fp] = rec; continue
                t0, s0 = o
                if t0 == tid: continue
                colls += 1
                r0 = replay(tw, starts[t0], s0, depth)
                r1 = replay(tw, starts[tid], k, depth)
                rep += s0 + k
                y0, A0, D0 = r0
                y1, A1, D1 = r1
                sh = tw.match_shift(y0, y1)
                if sh is None: degen += 1; continue
                kk = tw.recover(A0, D0, A1, D1, sh)
                if not tw.ok(kk): degen += 1; continue
                ok += 1; break
        per_inst.append(ok)
        tot_ok += ok; tot_tr += trials
    cp = cp95_lower(tot_ok, tot_tr)
    out = {"mode": "census", "q": q, "rp": rp, "N": N, "R": R, "dpbits": dpbits,
           "depth": depth, "instances": ninst, "trials_per_instance": trials,
           "trails": trails, "max_trail": max_trail,
           "per_instance_success": per_inst,
           "per_instance_trials": trials,
           "instances_full": sum(1 for v in per_inst if v == trials),
           "instances_zero": sum(1 for v in per_inst if v == 0),
           "success": tot_ok, "trials": tot_tr, "rate": tot_ok / tot_tr,
           "cp95_lower": str(cp),
           "min_instance_rate": min(per_inst) / trials,
           "collisions": colls, "degenerate": degen,
           "forward_steps": fwd, "replay_steps": rep, "records": recs,
           "record_bytes": 50, "protocol": "restarted worker: trail ends at FIRST "
           "report, one 50 B record, restart from a fresh seed-derived start"}
    print(json.dumps(out, indent=1))
    return out


def mode_guard(cfg):
    """T3 coalescence, fully independent depth-d history on both sides."""
    q, rp, N = cfg.get("q", 523), cfg["rp"], cfg["N"]
    R, depth = cfg.get("R", 128), cfg.get("depth", 2)
    pairs, steps = cfg.get("pairs", 1200), cfg.get("steps", 60)
    tw = T21(q, rp, R, 5, 1, depth)
    tw.new_instance(7)
    rng = random.Random(0xA121)
    cps = [1, 2, 5, 10, steps]
    res = {}
    for arm in ("indep_hist", "prevparity"):
        same = dict((c, 0) for c in cps)
        for _ in range(pairs):
            a, b = rng.randrange(1, rp), rng.randrange(1, rp)
            j = rng.randrange(1, 12)
            y1, A1, D1 = tw.start(a, b)
            y2 = y1
            for _ in range(j): y2 = tw.T.frob(y2)
            _ks, m1, k1 = tw.keys(y1)
            _ks, m2, k2 = tw.keys(y2)
            assert m1 == m2
            h1 = [rng.randrange(R) for _ in range(depth)]
            h2 = [rng.randrange(R) for _ in range(depth)]
            A2, D2 = A1, D1
            ci = 0
            for i in range(steps + 1):
                if ci < len(cps) and i == cps[ci]:
                    if tw.conj_equal(y1, y2): same[i] += 1
                    ci += 1
                if i < steps:
                    if arm == "prevparity":
                        b1 = (h1[0] & 1) == 0
                        b2 = (h2[0] & 1) == 0
                        if b1 != b2:
                            y1 = tw.T.cyclo_sqr(y1); y2 = tw.T.mul(y2, tw.T.one12())
                            A1, D1 = 2 * A1 % rp, 2 * D1 % rp
                            _ks, m1, k1 = tw.keys(y1); _ks, m2, k2 = tw.keys(y2)
                            continue
                    y1, A1, D1, k1, m1, _b, h1, _d = tw.fwd(y1, A1, D1, k1, m1, h1)
                    y2, A2, D2, k2, m2, _b, h2, _d = tw.fwd(y2, A2, D2, k2, m2, h2)
        res[arm] = {"pairs": pairs, "steps": steps,
                    "coalesce_at_1": same[1] / pairs, "coalesce_at_last": same[steps] / pairs,
                    "by_step": dict(("%d" % c, same[c] / pairs) for c in cps)}
    out = {"mode": "guard", "q": q, "rp": rp, "N": N, "R": R, "depth": depth,
           "statistic": "fraction of conjugate pairs still in ONE class; both sides "
           "start with INDEPENDENT uniform depth-d history (arbitrary-merge test); "
           "`prevparity` is the mutation control", "rules": res}
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    c = json.loads(sys.stdin.read() or "{}")
    m = c.pop("mode", "census")
    {"census": mode_census, "guard": mode_guard}[m](c)
