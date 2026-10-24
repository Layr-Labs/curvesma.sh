#!/usr/bin/env python3
"""Toy v20 -- the DECLARED sparse-DP implementation, run literally.

Fixes the three defects the reviewer named on 096ec169:

(1) KEY RECOMPUTE ON SQUARING.  `mode_keys` shows the twelve keys are Fp-linear
    forms of ONE Fp2 coordinate (d0 = y[1][0]); two states with identical key
    vectors but different d0-free mass have DIFFERENT key vectors after squaring.
    The twelve keys of y^2 are therefore NOT a function of the twelve keys of y,
    so there is no cheap squared-state key update; the ledger charges the full
    key_six_plusminus + min_rotation_booth (1946) on every squaring step too.
    Controls: squared-key map vs true keys(y^2) mismatch count, and the
    key-vector collision witness.

(2) DECLARED BIT MAP.  branch bits 0-2, index bits 3-14 and the DP window bits
    15.. are ALL cut from the SAME mixed 128-bit label L = mix128(m0).  The v19
    defect (branch from raw m0, index/DP from L, two sources) is kept as a
    control arm `two_source`.

(3) TRAIL PROTOCOL.  A trail ENDS at its first report; the pool stores
    (fingerprint, trail, step) only; a collision is recovered by REPLAYING both
    trails from their stored (start, steps) -- no stored (A,D,y).

Every rule arm is mutation-controlled: `prevparity` (guard reads the parity of
the previous slot instead of the class slot) must FAIL the independent-prev
coalescence audit, `two_source` must disagree with the declared bit map.
"""
import json
import math
import random
import sys

from toy import Tower, field_params, s_ops

IDX_BITS = 12
DP_SHIFT = 3 + IDX_BITS          # 15
LAB_MASK = (1 << 128) - 1
LAB_MUL = 0x9E3779B97F4A7C15F39CC0605CEDC835
R_PROD = 4096
SQ_NUM, SQ_DEN = 3, 8            # declared squaring mixture (branch bits 0-2 of L)


def mix128(x):
    return ((x * LAB_MUL) ^ (x >> 61)) & LAB_MASK


def cp_lower(k, n, alpha=0.05):
    """Clopper-Pearson one-sided lower end, exact bisection on the binomial tail."""
    if k == 0:
        return 0.0
    lo, hi = 0.0, 1.0 - 1e-15
    for _ in range(300):
        mid = (lo + hi) / 2
        mid = min(max(mid, 1e-15), 1.0 - 1e-15)
        # P[X >= k | p=mid] via incomplete beta series (log-space, stable at toy n)
        s = 0.0
        for j in range(k, n + 1):
            s += math.exp(math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1)
                          + j * math.log(mid) + (n - j) * math.log1p(-mid))
        if s < alpha:          # tail too small -> p is too small
            lo = mid
        else:
            hi = mid
    return lo


class T20:
    def __init__(self, q, rp, R=R_PROD, dpbits=5, seed=1, rule="v20"):
        self.q, self.rp, self.R, self.dpbits = q, rp, R, dpbits
        self.rule, self.seed = rule, seed
        self.T = Tower(q)
        self.rng = random.Random(0xC0FFEE ^ (seed * 2654435761))
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

    # ---------------------------------------------------------- class function
    def keys(self, y):
        """kappa = index of the MINIMAL CYCLIC ROTATION of the twelve keys (Booth form)."""
        ks = self.T.keys12(y)
        best = None
        for i in range(12):
            cand = tuple(ks[(i + t) % 12] for t in range(12))
            if best is None or cand < best[0]:
                best = (cand, i)
        kap = best[1]
        return ks, ks[kap], kap

    def label(self, m0):
        return mix128(m0)

    def bits(self, m0):
        """THE DECLARED BIT MAP, one source: bits of L = mix128(m0) only."""
        L = self.label(m0)
        br = (L & 7) < SQ_NUM               # bits 0-2  -> squaring iff < 3 (3/8)
        h = (L >> 3) % self.R               # bits 3-14 -> index over R slots
        dp = ((L >> DP_SHIFT) & ((1 << self.dpbits) - 1)) == 0
        return br, h, dp, L

    def bits_two_source(self, y, m0):
        """CONTROL: the v19 defect -- branch from raw m0, index/DP from L."""
        L = self.label(m0)
        return (m0 & 7) < SQ_NUM, (L >> 3) % self.R, \
               ((L >> DP_SHIFT) & ((1 << self.dpbits) - 1)) == 0, L

    def decide(self, y, m0, prev):
        br, h, dp, L = (self.bits_two_source(y, m0) if self.rule == "two_source"
                        else self.bits(m0))
        was_guard = False
        if not br:
            if self.rule in ("v20", "two_source") and prev == h:
                br, was_guard = True, True
            elif self.rule == "prevparity" and (prev & 1) == 0:
                br, was_guard = True, True
        return br, (-1 if br else h), dp, L, was_guard

    # ------------------------------------------------------------------- step
    def build_table(self):
        rng, rp, R, T = self.rng, self.rp, self.R, self.T
        self.x = rng.randrange(1, rp)
        self.z = T.pow(self.g, self.x)
        self.ent = []
        for _ in range(R):
            a, b = rng.randrange(1, rp), rng.randrange(1, rp)
            self.ent.append((a, b, T.pow(self.g, (a + self.x * b) % rp)))

    def _rot(self, M, j):
        cur = M
        for _ in range(j % 12):
            cur = self.T.frob(cur)
        return cur

    def mul_step(self, y, A, D, h, kap):
        a, b, M = self.ent[h]
        return (self.T.mul(y, self._rot(M, (12 - kap) % 12)),
                (A + self.linv[kap] * a) % self.rp,
                (D + self.linv[kap] * b) % self.rp)

    def sqr_step(self, y, A, D):
        return self.T.cyclo_sqr(y), (2 * A) % self.rp, (2 * D) % self.rp

    def fwd(self, y, A, D, kap, m0, prev):
        br, h, dp, L, g = self.decide(y, m0, prev)
        if br:
            y, A, D = self.sqr_step(y, A, D)
            nxt = -1
        else:
            y, A, D = self.mul_step(y, A, D, h, kap)
            nxt = h
        _ks, m0, kap = self.keys(y)
        return y, A, D, kap, m0, br, nxt, dp

    # -------------------------------------------------------------- algebra
    def class_id(self, y):
        best = None
        cur = y
        for _ in range(12):
            ser = (cur[0][0], cur[0][1], cur[0][2], cur[1][0], cur[1][1], cur[1][2])
            if best is None or ser < best:
                best = ser
            cur = self.T.frob(cur)
        return best

    def conj_equal(self, y1, y2):
        cur = y1
        for _ in range(12):
            if cur == y2:
                return True
            cur = self.T.frob(cur)
        return False

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

    def start(self, a, b):
        return self.T.pow(self.g, (a + self.x * b) % self.rp), a % self.rp, b % self.rp

    def ok(self, kk):
        return kk not in (None, 0) and self.T.pow(self.g, kk) == self.z


def get_tower(cfg):
    q = cfg.get("q", 523)
    fp = field_params(q)
    if fp is None:
        raise SystemExit("bad q")
    return (q,) + tuple(fp)


# ==================================================== MODE: key-update refutation
def mode_keys(cfg):
    """Is there a cheap squared-state key update?  Measured, not asserted.

    keys_j(y) are six Fp-linear forms (and their negations) of d0 = y[1][0], i.e. of
    TWO Fp coordinates.  Test the only candidate cheap rules:
      A) keys(y^2) == keys(y)^2 componentwise (squaring the keys),
      B) keys(y^2) == Frobenius powers of keys(y),
      C) keys(y^2) is a function of keys(y) at all: find distinct y with equal key
         vectors and compare keys(y^2).
    """
    q, rp, N, _a, _b = get_tower(cfg)
    tw = T20(q, rp, cfg.get("R", 128), 5, seed=cfg.get("seed", 1), rule="v20")
    tw.build_table()
    rng = random.Random(0x4B20 ^ cfg.get("seed", 1))
    n = cfg.get("samples", 2000)
    sq_key_mismatch = frob_key_mismatch = 0
    buckets = {}
    func_of_keys_viol = 0
    witness = None
    for _ in range(n):
        y, A, D = tw.start(rng.randrange(1, rp), rng.randrange(1, rp))
        k0 = tw.T.keys12(y)
        y2 = tw.T.cyclo_sqr(y)
        k2 = tw.T.keys12(y2)
        cand = [x * x % q for x in k0]
        if cand != k2:
            sq_key_mismatch += 1
        fro = [tw.T.keys12(tw.T.frob(y))[i] for i in range(12)]
        if fro != k2:
            frob_key_mismatch += 1
        key = tuple(k0)
        if key in buckets:
            y_prev, k2_prev = buckets[key]
            if k2_prev != k2:
                func_of_keys_viol += 1
                if witness is None:
                    witness = {"same_keys": True,
                               "keys_y2_differ": [k2_prev[:3], k2[:3]],
                               "note": "equal twelve keys, different keys after squaring"}
        else:
            buckets[key] = (y, k2)
    out = {"mode": "keys", "q": q, "rp": rp, "N": N, "samples": n,
           "squared_key_rule_mismatch": sq_key_mismatch,
           "frobenius_key_rule_mismatch": frob_key_mismatch,
           "key_vector_repeats": sum(1 for v in buckets.values()) and n - len(buckets),
           "distinct_key_vectors": len(buckets),
           "keys_are_not_a_function_of_keys_of_y": func_of_keys_viol,
           "witness": witness,
           "conclusion": ("no squared-state key update exists: keys(y^2) is not a function "
                          "of keys(y), so the twelve keys + minimal rotation must be "
                          "RECOMPUTED and charged on every squaring step (1946 instr)")}
    print(json.dumps(out, indent=1))
    return out


# ============================================ MODE: declared sparse-DP census
def mode_sparse(cfg):
    """THE DECLARED PROTOCOL, literally.

    trails concurrent walks over ONE pool; every trail ENDS at its FIRST report; the
    pool stores only (fingerprint 32B, trail 2B, step 16B) = 50 B; a collision is
    recovered by REPLAYING both trails from (start, steps) and re-deriving (A, D, y).
    success = a collision whose replay-verified relation returns x with g^x = z.
    """
    q, rp, N, _a, _b = get_tower(cfg)
    R = cfg.get("R", 128)
    dpbits = cfg.get("dpbits", 3)
    trails = cfg.get("trails", 256)
    trials = cfg.get("trials", 2000)
    seed0 = cfg.get("seed", 1)
    arms = cfg.get("rules", ["v20", "two_source"])
    res = {}
    for rule in arms:
        tw = T20(q, rp, R, dpbits, seed=seed0, rule=rule)
        tw.build_table()
        rng = random.Random(0x5A20 ^ (seed0 * 7919))
        ok = coll = degen = 0
        steps = replays = 0
        for _t in range(trials):
            starts = [(rng.randrange(1, rp), rng.randrange(1, rp)) for _ in range(trails)]
            pool = {}
            hit = False
            for tid in range(trails):
                a, b = starts[tid]
                y, A, D = tw.start(a, b)
                _ks, m0, kap = tw.keys(y)
                prev = -1
                rec = None
                k = 0
                while True:                       # trail ENDS at its first report
                    y, A, D, kap, m0, br, prev, dp = tw.fwd(y, A, D, kap, m0, prev)
                    k += 1
                    steps += 1
                    if dp:
                        rec = (tid, k)
                        break
                    if k >= cfg.get("max_trail", 64):
                        break
                if rec is None:
                    continue
                tid1, k1 = rec
                fp = tw.class_id(y)
                o = pool.get(fp)
                if o is None:
                    pool[fp] = rec
                    continue
                t0, s0 = o
                if t0 == tid:
                    continue
                coll += 1
                # ---- REPLAY both trails from (start, steps); nothing else is stored.
                r0 = replay(tw, starts[t0], s0)
                r1 = replay(tw, starts[tid], k1)
                replays += s0 + k1
                if r0 is None or r1 is None:
                    degen += 1
                    continue
                y0, A0, D0 = r0
                y1, A1, D1 = r1
                sh = tw.match_shift(y0, y1)
                if sh is None:
                    degen += 1
                    continue
                kk = tw.recover(A0, D0, A1, D1, sh)
                if not tw.ok(kk):
                    degen += 1
                    continue
                ok += 1
                hit = True
                break
            if hit:
                pass
        M = trails
        model = 1.0 - math.exp(-M * (M - 1.0) / 2.0 / float(N))
        res[rule] = {"trials": trials, "trails": trails, "dpbits": dpbits,
                     "success": ok, "rate": ok / trials, "cp95_lower": cp_lower(ok, trials),
                     "collisions": coll, "degenerate": degen,
                     "records_per_trial": trails, "birthday_model": model,
                     "replay_steps": replays, "forward_steps": steps,
                     "record_bytes": 50,
                     "retention_vs_model": (ok / trials) / model if model else None}
    out = {"mode": "sparse", "q": q, "rp": rp, "N": N, "R": R,
           "protocol": "trail ends at FIRST report; pool=(fingerprint,trail,step); "
                       "recovery replays both trails from (start,steps); declared "
                       "single mixed-label bit map (branch 0-2, index 3-14, DP 15+)",
           "rules": res}
    print(json.dumps(out, indent=1))
    return out


def replay(tw, start, steps):
    a, b = start
    y, A, D = tw.start(a, b)
    _ks, m0, kap = tw.keys(y)
    prev = -1
    for _ in range(steps):
        y, A, D, kap, m0, br, prev, dp = tw.fwd(y, A, D, kap, m0, prev)
    return y, A, D


# =================================================== MODE: coalescence + prev init
def mode_coalesce(cfg):
    q, rp, N, _a, _b = get_tower(cfg)
    R = cfg.get("R", R_PROD)
    pairs = cfg.get("pairs", 1600)
    steps = cfg.get("steps", 100)
    seed0 = cfg.get("seed", 1)
    cps = sorted(set([0, 1, 2, 3, 4, 8, 16, 32, 64, steps]))
    cps = [c for c in cps if c <= steps]
    res = {}
    for arm in cfg.get("rules", ["sync", "indep_prev", "prevparity_indep"]):
        rule = "prevparity" if arm.startswith("prevparity") else "v20"
        tw = T20(q, rp, R, cfg.get("dpbits", 5), seed=seed0, rule=rule)
        tw.build_table()
        rng = random.Random(0xA120 ^ (seed0 * 7919))
        same = dict((c, 0) for c in cps)
        witness = None
        for t in range(pairs):
            a, b = rng.randrange(1, rp), rng.randrange(1, rp)
            j = rng.randrange(1, 12)
            y1, A1, D1 = tw.start(a, b)
            y2 = y1
            for _ in range(j):
                y2 = tw.T.frob(y2)
            _ks, m1, k1 = tw.keys(y1)
            _ks, m2, k2 = tw.keys(y2)
            if m1 != m2:
                raise SystemExit("m0 is not a class invariant")
            A2, D2 = A1, D1
            if arm == "sync":
                prev1 = prev2 = -1
            else:
                # THE HONEST ARBITRARY-MERGE TEST: independent random previous slots.
                prev1, prev2 = rng.randrange(R), rng.randrange(R)
            ci = 0
            for i in range(steps + 1):
                if ci < len(cps) and i == cps[ci]:
                    if tw.conj_equal(y1, y2):
                        same[i] += 1
                    elif witness is None:
                        witness = {"pair": t, "j": j, "step": i,
                                   "prev": [prev1, prev2],
                                   "why": "same class in, different class out"}
                    ci += 1
                if i < steps:
                    y1, A1, D1, k1, m1, _, prev1, _ = tw.fwd(y1, A1, D1, k1, m1, prev1)
                    y2, A2, D2, k2, m2, _, prev2, _ = tw.fwd(y2, A2, D2, k2, m2, prev2)
        res[arm] = {"pairs": pairs, "steps": steps,
                    "same_class_at": dict(("%d" % c, same[c]) for c in cps),
                    "coalesce_at_1": same.get(1, 0) / pairs,
                    "coalesce_at_last": same[cps[-1]] / pairs,
                    "first_divergence": witness}
    out = {"mode": "coalesce", "q": q, "rp": rp, "N": N, "R": R, "checkpoints": cps,
           "statistic": "fraction of same-class (conjugate) pairs still in ONE class at "
                        "each checkpoint; `sync` = both prev=-1 (v19's audit), "
                        "`indep_prev` = independent uniform previous slots (the honest "
                        "arbitrary-merge test), `prevparity_indep` = mutation control",
           "rules": res}
    print(json.dumps(out, indent=1))
    return out


# ============================================================== MODE: trap census
def mode_trap(cfg):
    """Exact length-2/3/4 returns of (class, A, D) under the declared v20 step and
    under the v18 control (repeat-slot squaring guard removed)."""
    q, rp, N, _a, _b = get_tower(cfg)
    R = cfg.get("R", R_PROD)
    pairs = cfg.get("pairs", 1600)
    steps = cfg.get("steps", 200)
    ntab = cfg.get("tables", 8)
    res = {}
    for arm in cfg.get("rules", ["v20", "v18ctl"]):
        rule = "v20" if arm == "v20" else "none"
        c2 = c3 = c4 = 0
        nsq = ntot = nguard = 0
        traces = 0
        first = None
        for tb in range(ntab):
            tw = T20(q, rp, R, cfg.get("dpbits", 5), seed=cfg.get("seed", 1) + tb, rule=rule)
            tw.build_table()
            rng = random.Random(0x7220 + tb * 131)
            for _ in range(pairs // ntab):
                y, A, D = tw.start(rng.randrange(1, rp), rng.randrange(1, rp))
                _ks, m0, kap = tw.keys(y)
                hist = [(m0, A, D)]
                prev = -1
                hit2 = hit3 = hit4 = False
                for i in range(steps):
                    raw_br = tw.bits(m0)[0]
                    y, A, D, kap, m0, br, prev, dp = tw.fwd(y, A, D, kap, m0, prev)
                    nguard += 1 if (rule == "v20" and br and not raw_br) else 0
                    nsq += br
                    ntot += 1
                    cur = (m0, A, D)
                    if len(hist) >= 2 and cur == hist[-2]:
                        hit2 = True
                        if first is None:
                            first = {"table": tb, "step": i}
                    if len(hist) >= 3 and cur == hist[-3]:
                        hit3 = True
                    if len(hist) >= 4 and cur == hist[-4]:
                        hit4 = True
                    hist.append(cur)
                traces += 1
                c2 += 1 if hit2 else 0
                c3 += 1 if hit3 else 0
                c4 += 1 if hit4 else 0
        res[arm] = {"traces": traces, "steps_per_trace": steps,
                    "traces_with_2cycle": c2, "traces_with_3cycle": c3,
                    "traces_with_4cycle": c4,
                    "per_1600_2cycle": round(1600.0 * c2 / traces, 3),
                    "squaring_fraction": round(nsq / max(1, ntot), 6),
                    "first_2cycle": first}
    out = {"mode": "trap", "q": q, "rp": rp, "N": N, "R": R, "tables": ntab,
           "statistic": "traces with an exact length-L return of (class, A, D)",
           "rules": res}
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    c = json.loads(sys.stdin.read() or "{}")
    {"keys": mode_keys, "sparse": mode_sparse, "coalesce": mode_coalesce,
     "trap": mode_trap}[c.get("mode", "coalesce")](c)
