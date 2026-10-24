#!/usr/bin/env python3
"""Toy v19 -- production step with a REPEAT-SLOT GUARD (no two consecutive multiply
steps may use the same table slot), the merge/coalescence audit, the p^6 fixture, the
second-order census and a control rule that hashes the sigma-invariant pair (t, n).

WHY v18 FAILED (our own trap census, 43,325/1600):  with ONE conjugacy-constant table,
every slot serves kappa and kappa+6 -- the exact pair with inverse multipliers.  If two
CONSECUTIVE multiply steps use the same slot h and the second has kappa' = kappa + 6,
then M_h^(p^-kappa) * M_h^(p^-kappa-6) = 1 and the (class, A, D) state returns exactly:
an absorbing 2-cycle.  The rate is (1/R)(1/12) per step, independent of what the slot
HASHES (only of how many slots exist), so no re-hashing can remove it.  At production
R=4096 that is 2.03e-5 per step against a trail of S/R ~ 2e34 steps: every trail traps.

THE v19 STEP keeps v18's Frobenius-equivariant map (minimal-rotation kappa, class-function
slot s = idx(H(m0)), multiplier M_s^(p^-kappa)) and adds ONE class-history guard:

    if the previous step was a multiply in slot h and the current slot is again h,
    take the cyclotomic squaring step instead.

The guard is a function of the CLASS TRAJECTORY only (slots and branch bits are functions
of m0, hence of the class), so two conjugate states -- which have the same class at every
step -- make the same decision at every step and their successors stay conjugate.  The
m^2 cancellation needs two consecutive multiplies in one slot, which the guard makes
unreachable, so the second-order channel is structurally empty at every scale, not just at
toy scale.

CONTROL RULE "tn" (the design seed, measured here rather than assumed): hash the pair
(t, n) = (z + z^sigma, z * z^sigma) for sigma = Frobenius^6.  n is Frobenius-invariant but
t is only sigma-invariant: t(y^(p^j)) = Frobenius^j(t(y)) != t(y) for j not divisible by 6,
so the slot is NOT a class function and the map is not equivariant.  The coalescence audit
measures the consequence.
"""
import json
import random
import sys

from toy import Tower, field_params, s_ops

IDX_BITS = 12            # production index field: bits 3..14  (4096 slots)
DP_SHIFT = 3 + IDX_BITS  # = 15, production DP window offset
LAB_MASK = (1 << 128) - 1
LAB_MUL = 0x9E3779B97F4A7C15F39CC0605CEDC835
R_PROD = 4096
RULES = ("v16", "v17", "v18", "v19", "tn")


def mix128(x):
    """Fixed universal multiply-shift hash, 128 bits (the fingerprint mixer)."""
    return ((x * LAB_MUL) ^ (x >> 61)) & LAB_MASK


class T19:
    def __init__(self, q, rp, R=R_PROD, dpbits=5, seed=1, rule="v19", use_sqr=True):
        assert rule in RULES and R % 2 == 0
        self.q, self.rp, self.R, self.dpbits = q, rp, R, dpbits
        self.rule, self.use_sqr, self.seed = rule, use_sqr, seed
        self.half = R // 2
        self.idxw = self.half if rule == "v17" else R
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

    # ------------------------------------------------------------------ key rule
    def keys(self, y):
        """m0 = min of the minimal cyclic rotation of the twelve Frobenius keys.

        kappa = argmin(ks) is a Frobenius torsor coordinate only when the minimum is
        unique; the minimal CYCLIC ROTATION is invariant under rotation of the key
        vector, so kappa(y^(p^j)) = kappa(y) - j exactly and m0 = ks[kappa] is a class
        invariant.  (v18's certified tie-break, kept verbatim.)
        """
        ks = self.T.keys12(y)
        if self.rule == "v18a":          # control: bare argmin, no tie-break
            m0 = min(ks)
            return ks, m0, ks.index(m0)
        best = None
        for i in range(12):
            cand = tuple(ks[(i + t) % 12] for t in range(12))
            if best is None or cand < best[0]:
                best = (cand, i)
        kap = best[1]
        return ks, ks[kap], kap

    # --------------------------------------------------------- sigma invariants
    def sigma(self, y):
        cur = y
        for _ in range(6):
            cur = self.T.frob(cur)
        return cur

    def tn_key(self, y):
        """(t, n) = (y + y^sigma, y * y^sigma): the seed's slot source.

        Complete invariant of the ORDER-TWO orbit {y, y^sigma} -- and therefore NOT a
        function of the order-twelve Frobenius class, which is the point of the audit.
        """
        s = self.sigma(y)
        # plain componentwise addition in the tower
        add = lambda u, v: tuple(tuple(((u[i][j][0] + v[i][j][0]) % self.q,
                                        (u[i][j][1] + v[i][j][1]) % self.q)
                                       for j in range(3)) for i in range(2))
        tt = add(y, s)
        nn = self.T.mul(y, s)
        ser = 0
        for i in range(2):
            for j in range(3):
                ser = (ser * self.q + tt[i][j][0]) % (1 << 128)
                ser = (ser * self.q + tt[i][j][1]) % (1 << 128)
        for i in range(2):
            for j in range(3):
                ser = (ser * self.q + nn[i][j][0]) % (1 << 128)
                ser = (ser * self.q + nn[i][j][1]) % (1 << 128)
        return ser

    def _label_from(self, L, t):
        """THE production bit map: disjoint ranges of the 128-bit class label.

        branch bits 0..2 | index bits 3..14 | DP window bits 15..15+dpbits.
        """
        br = ((t & 7) < 3) and self.use_sqr
        h = (L >> 3) % self.idxw
        dp = ((L >> DP_SHIFT) & ((1 << self.dpbits) - 1)) == 0
        return br, h, dp

    def bits(self, y, m0):
        """v16 reads a REPRESENTATIVE quantity; v17/v18/v19 read the class invariant."""
        if self.rule == "v16":
            return self._label_from(mix128(y[0][0][0]), y[0][0][0])
        if self.rule == "tn":
            return self._label_from(mix128(self.tn_key(y)), self.tn_key(y))
        return self._label_from(mix128(m0), m0)

    # ---------------------------------------------------------------- the step
    def build_table(self):
        rng, rp, R, T = self.rng, self.rp, self.R, self.T
        self.x = rng.randrange(1, rp)
        self.z = T.pow(self.g, self.x)
        self.ent = []
        for _ in range(R):
            a, b = rng.randrange(1, rp), rng.randrange(1, rp)
            # M_s = g^a * z^b = g^(a + x*b): one exponentiation per slot; the label
            # (a, b) is what the step adds to the exponent pair (A, D).
            self.ent.append((a, b, T.pow(self.g, (a + self.x * b) % rp)))

    def entry(self, h, kap):
        """Slot and p^-kappa rotation.  Only the SLOT rule changes between rules."""
        if self.rule == "v17":
            a, b, M = self.ent[h + (self.half if kap >= 6 else 0)]
        else:
            a, b, M = self.ent[h]
        return a, b, self._rot(M, (12 - kap) % 12)

    def _rot(self, M, j):
        if j == 0:
            return M
        cur = M
        for _ in range(j):
            cur = self.T.frob(cur)
        return cur

    def step(self, y, A, D, h, kap):
        a, b, M = self.entry(h, kap)
        return (self.T.mul(y, M),
                (A + self.linv[kap] * a) % self.rp,
                (D + self.linv[kap] * b) % self.rp)

    def sqr(self, y, A, D):
        return self.T.cyclo_sqr(y), (2 * A) % self.rp, (2 * D) % self.rp

    def decide(self, y, m0, prev):
        """(branch, slot).  v19: the REPEAT-SLOT GUARD.

        prev = slot of the immediately preceding MULTIPLY step, or -1 if the preceding
        step was a squaring (nothing to cancel against).  A multiply in a slot that
        repeats the previous multiply slot is replaced by a squaring, so the m^2
        cancellation M^(p^-kappa) * M^(p^-kappa-6) = 1 is unreachable in two steps.

        The guard reads only the class label and the previous class label, so it is a
        function of the class trajectory: conjugate states take identical decisions.
        """
        br, h, _dp = self.bits(y, m0)
        if self.rule == "v19" and not br and prev == h:
            br, h = True, -1
        return br, h

    def fwd(self, y, A, D, kap, m0, prev=-1):
        br, h = self.decide(y, m0, prev)
        if br:
            y, A, D = self.sqr(y, A, D)
            nxt = -1
        else:
            y, A, D = self.step(y, A, D, h, kap)
            nxt = h
        _ks, m0, kap = self.keys(y)
        return y, A, D, kap, m0, br, nxt

    # ------------------------------------------------------------ class algebra
    def conj_equal(self, y1, y2):
        """EXACT test: is y2 a Frobenius conjugate of y1 (same class)?"""
        cur = y1
        for _ in range(12):
            if cur == y2:
                return True
            cur = self.T.frob(cur)
        return False

    def class_id(self, y):
        """Exact class id at toy scale: lexicographic minimum coordinate tuple over
        the twelve conjugates.  Production stores m0 (381 bits), whose collision
        probability is negligible; the toy must be exact so a record match is a class
        match and never a fingerprint accident."""
        best = None
        cur = y
        for _ in range(12):
            ser = (cur[0][0], cur[0][1], cur[0][2], cur[1][0], cur[1][1], cur[1][2])
            if best is None or ser < best:
                best = ser
            cur = self.T.frob(cur)
        return best

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

    def inv_ok(self, y, A, D):
        return y == self.T.pow(self.g, (A + self.x * D) % self.rp)


def trail_start(tw, tid):
    rr = random.Random(0x5EED0000 + tw.seed * 1000003 + tid)
    return rr.randrange(1, tw.rp), rr.randrange(1, tw.rp)


def cfg_in():
    return json.loads(sys.stdin.read() or "{}")


def get_tower(cfg):
    q = cfg.get("q", 523)
    fp = field_params(q)
    if fp is None:
        raise SystemExit("bad q")
    return (q,) + tuple(fp)


# =============================================================== MODE: coalesce
def mode_coalesce(cfg):
    """MERGE AUDIT.  Start y and a non-trivial conjugate y^(p^j) -- the SAME class --
    walk the two independently with the production step, and test at fixed checkpoints
    whether the two states are still in one class (exact twelve-Frobenius test)."""
    q, rp, N, _a, _b = get_tower(cfg)
    R = cfg.get("R", R_PROD)
    pairs = cfg.get("pairs", 1600)
    steps = cfg.get("steps", 100)
    seed0 = cfg.get("seed", 1)
    jlo = cfg.get("jmin", 1)
    cps = sorted(set([0, 1, 2, 3, 4, 8, 16, 32, 64, steps]))
    cps = [c for c in cps if c <= steps]
    out_rules = {}
    for rule in cfg.get("rules", ["v19", "v18", "tn"]):
        tw = T19(q, rp, R, cfg.get("dpbits", 5), seed=seed0, rule=rule)
        tw.build_table()
        rng = random.Random(0xA11D ^ (seed0 * 7919))
        same = dict((c, 0) for c in cps)
        elem_eq = 0
        witness = None
        for t in range(pairs):
            a, b = rng.randrange(1, rp), rng.randrange(1, rp)
            j = rng.randrange(jlo, 12)
            y1, A1, D1 = tw.start(a, b)
            y2 = y1
            for _ in range(j):
                y2 = tw.T.frob(y2)
            assert tw.conj_equal(y1, y2)
            _ks, m1, k1 = tw.keys(y1)
            _ks, m2, k2 = tw.keys(y2)
            assert m1 == m2, "m0 must be a class invariant"
            A2, D2 = A1, D1
            prev1 = prev2 = -1
            ci = 0
            for i in range(steps + 1):
                if ci < len(cps) and i == cps[ci]:
                    if tw.conj_equal(y1, y2):
                        same[i] += 1
                    elif witness is None:
                        witness = {"pair": t, "j": j, "step": i, "m0_equal": m1 == m2,
                                   "slot": [tw.bits(y1, m1)[1], tw.bits(y2, m2)[1]],
                                   "kappa": [k1, k2],
                                   "why": "same class entering, different class leaving: "
                                          "the step is not a function of the class"}
                    ci += 1
                if i < steps:
                    y1, A1, D1, k1, m1, _, prev1 = tw.fwd(y1, A1, D1, k1, m1, prev1)
                    y2, A2, D2, k2, m2, _, prev2 = tw.fwd(y2, A2, D2, k2, m2, prev2)
                    if i == 0 and y1 == y2:
                        elem_eq += 1
        out_rules[rule] = {
            "pairs": pairs, "steps": steps, "R": R,
            "same_class_at": dict(("%d" % c, same[c]) for c in cps),
            "coalesce_within_1_step": same.get(1, 0) / pairs,
            "coalesce_within_2_steps": same.get(2, 0) / pairs,
            "same_class_at_last_step": same[cps[-1]] / pairs,
            "identical_elements_after_1_step": elem_eq,
            "first_divergence_witness": witness,
        }
    out = {"mode": "coalesce", "q": q, "rp": rp, "N": N, "R": R,
           "checkpoints": cps, "idx_bits": IDX_BITS, "dp_shift": DP_SHIFT,
           "statistic": "fraction of same-class (conjugate) start pairs still in ONE class "
                        "at each checkpoint, exact twelve-Frobenius test",
           "rules": out_rules}
    print(json.dumps(out, indent=1))
    return out


# =================================================================== MODE: trap
def mode_trap(cfg):
    """Degeneracy census under the step rules on identical starts.

    A trace is degenerate if the walk returns to a state whose exponent pair AND class
    equal the state two steps earlier: the p^6 signature (an exact 2-cycle, which
    restores the labels too, so no logarithm is recoverable and the trail is wasted).
    Also reports the squaring mixture and the guard firing rate."""
    q, rp, N, _a, _b = get_tower(cfg)
    R = cfg.get("R", R_PROD)
    pairs = cfg.get("pairs", 1600)
    steps = cfg.get("steps", 200)
    ntab = cfg.get("tables", 8)
    seed0 = cfg.get("seed", 1)
    res = {}
    for rule in cfg.get("rules", ["v16", "v17", "v18", "v19", "tn"]):
        c2 = c3 = 0
        nsq = ntot = nguard = 0
        traces = 0
        first = None
        for tb in range(ntab):
            tw = T19(q, rp, R, cfg.get("dpbits", 5), seed=seed0 + tb, rule=rule)
            tw.build_table()
            rng = random.Random(0x5150 + tb * 131)
            for _ in range(pairs // ntab):
                y, A, D = tw.start(rng.randrange(1, rp), rng.randrange(1, rp))
                _ks, m0, kap = tw.keys(y)
                hist = [(m0, A, D)]
                hit2 = hit3 = False
                prev = -1
                for i in range(steps):
                    raw = tw.bits(y, m0)[0]
                    br, h = tw.decide(y, m0, prev)
                    was_guard = bool(rule == "v19" and prev >= 0 and not raw and br)
                    y, A, D, kap, m0, br2, nxt = tw.fwd(y, A, D, kap, m0, prev)
                    prev = nxt
                    nguard += 1 if was_guard else 0
                    nsq += br2
                    ntot += 1
                    cur = (m0, A, D)
                    # hist holds s_0..s_i, so hist[-2] is TWO steps back: an exact
                    # length-2 return of (class, A, D) restores the labels and the
                    # element, i.e. the p^6 inverse-multiplier 2-cycle.
                    if len(hist) >= 2 and cur == hist[-2]:
                        hit2 = True
                        if first is None:
                            first = {"table": tb, "step": i}
                    if len(hist) >= 4 and cur == hist[-4]:
                        hit3 = True
                    hist.append(cur)
                traces += 1
                c2 += 1 if hit2 else 0
                c3 += 1 if hit3 else 0
        res[rule] = {"traces": traces, "steps_per_trace": steps,
                     "traces_with_2cycle": c2, "traces_with_3cycle": c3,
                     "per_1600_2cycle": round(1600.0 * c2 / traces, 3),
                     "rate_per_step": c2 / max(1, traces * steps),
                     "squaring_fraction": round(nsq / max(1, ntot), 6),
                     "guard_firings": nguard,
                     "first_2cycle": first}
    out = {"mode": "trap", "q": q, "rp": rp, "N": N, "R": R, "tables": ntab,
           "statistic": "traces with an exact length-2 return of (class, A, D); the p^6 "
                        "inverse-multiplier cancellation is structural whenever one slot "
                        "can serve kappa and kappa+6 in two consecutive multiply steps",
           "rules": res}
    print(json.dumps(out, indent=1))
    return out


# ==================================================================== MODE: tn
def mode_tn(cfg):
    """THE DESIGN-SEED TEST, measured.

    For random y and j = 1..11, the seed's slot source is the pair (t, n) of invariants
    of the ORDER-TWO orbit {y, y^sigma}, sigma = Frobenius^6.  Report:
      * whether n is Frobenius invariant (it must be: n = y^(1+p^6) in Fp6),
      * whether t is Frobenius invariant (it is only sigma-invariant),
      * the fraction of conjugate pairs whose (t, n) label AND slot agree.
    """
    q, rp, N, _a, _b = get_tower(cfg)
    R = cfg.get("R", R_PROD)
    tw = T19(q, rp, R, cfg.get("dpbits", 5), seed=cfg.get("seed", 1), rule="tn")
    tw.build_table()
    rng = random.Random(0x7E57 ^ (cfg.get("seed", 1) * 7919))
    n_inv = n_tot = t_inv = lab_eq = slot_eq = 0
    add = lambda u, v: tuple(tuple(((u[i][j][0] + v[i][j][0]) % q,
                                    (u[i][j][1] + v[i][j][1]) % q)
                                   for j in range(3)) for i in range(2))
    for _ in range(cfg.get("samples", 400)):
        y, _A, _D = tw.start(rng.randrange(1, rp), rng.randrange(1, rp))
        s = tw.sigma(y)
        n0 = tw.T.mul(y, s)
        t0 = add(y, s)
        lab0 = tw.tn_key(y)
        slot0 = tw.bits(y, min(tw.keys(y)[0]))[1]
        cur = y
        for j in range(1, 12):
            cur = tw.T.frob(cur)
            sj = tw.sigma(cur)
            nj = tw.T.mul(cur, sj)
            tj = add(cur, sj)
            n_tot += 1
            n_inv += 1 if nj == n0 else 0
            t_inv += 1 if tj == t0 else 0
            lab_eq += 1 if tw.tn_key(cur) == lab0 else 0
            slot_eq += 1 if tw.bits(cur, min(tw.keys(cur)[0]))[1] == slot0 else 0
    out = {"mode": "tn", "q": q, "rp": rp, "N": N, "R": R,
           "samples": cfg.get("samples", 400), "pairs": n_tot,
           "norm_frobenius_invariant": n_inv, "trace_frobenius_invariant": t_inv,
           "tn_label_agrees_on_conjugates": lab_eq,
           "tn_slot_agrees_on_conjugates": slot_eq,
           "verdict": ("(t,n) is a complete invariant of the ORDER-TWO orbit only; it is "
                       "NOT a class function on the order-twelve Frobenius class, so a "
                       "slot keyed on it splits conjugates and the map is not equivariant")
           if lab_eq < n_tot else "class function on this sample"}
    print(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    c = cfg_in()
    m = c.get("mode", "coalesce")
    {"coalesce": mode_coalesce, "trap": mode_trap, "tn": mode_tn}[m](c)
