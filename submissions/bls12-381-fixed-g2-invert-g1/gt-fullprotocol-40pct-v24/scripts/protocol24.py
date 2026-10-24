#!/usr/bin/env python3
"""protocol24.py -- THE SHIPPED PROTOCOL, METERED, v24.

v23 flaw 1 (verdict c3ec4294): measure23.one_fixture ran rp=521041, depth 2, 128
trails, dpbits=3 and charged FORWARD STEPS ONLY.  Table construction, starting
exponentiations, replay, lookup and verification were not charged to any meter, and
the two-scale test kept the same group.  This file replaces it.

PROTOCOL (one implementation, used by census and survival):
  * global meter in INSTRUCTIONS; before every action the worker asks the meter
    whether the action fits; if not it ABANDONS (trail dies, its record is kept).
  * every action is charged: table entry build, per-trail start exponentiation,
    walk step (multiply price / squaring price), record write, per-step table
    lookup, replay of a collided pair, and candidate verification.
  * the walk is the shipped transition: depth-5 slot-history guard, slot from a
    128-bit mix of the FULL class identity, kappa rotation (12-kappa) on multiply.
  * dpbits is the shipped record density.  v24 ships dpbits=0 (every visited class
    is recorded) because the merge-stay decay measured on the shipped map makes any
    spacing > 1 step unable to convert a coalescence into a record match: the
    collision must be REACHED by both workers at a recording step.  With dpbits=0
    the required stay-merge horizon is 1 step, which is inside the measured range.

CHARGES are declared in gen24.py from counted/derived rows; the meter uses them.

MODES
  census   : F fresh fixtures x 1 trial, meter-first, success = solved before
             exhaustion, counted directly.  Two groups (different rp).
  survival : coalescence survival at Rg=65536, arms + kappa-offset resolution.
  scan     : observer-only scan on the SHIPPED step function (initialsed from this
             module's step_msg -- same source symbol), full Rg multiplier index.
  cpctl    : planted control for the Clopper-Pearson helper.
"""
import json
import math
import os
import random
import sys
import time
from decimal import Decimal, getcontext
from fractions import Fraction as F

getcontext().prec = 80
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from toy import Tower, field_params                      # noqa: E402
from census23 import reduce_counts                       # noqa: E402

# ---------------------------------------------------------------- cost rows (Ir)
MUL_REF = 20000
SQR_REF = 27000
KEY = 1487
BOOTH = 405
TABLE = 150
EXPADD = 60
EXPDBL = 40
COMMON = 331
MUL_CHARGE = MUL_REF + KEY + BOOTH + TABLE + EXPADD + COMMON   # 22433
SQR_CHARGE = SQR_REF + EXPDBL + KEY + BOOTH + COMMON            # 29263
LOOKUP = 96          # per-step pool probe (row derivation in gen24.py)
REC_WRITE = 64       # record store: key write + label pair write (derivation)
VERIFY = 24002000    # e([x]P0,Q)=z : 1200 dense Fp12 muls + G1 scalar action
START_MULS = 383     # 2-scalar 255-bit exponentiation = 255 dbl + 128 add
START_IR = START_MULS * MUL_REF
ENTRY_MULS = 383     # one table entry: g^a z^b = 2 exps ~ 2*383 dense muls
ENTRY_IR = 2 * ENTRY_MULS * MUL_REF

M128 = (1 << 128) - 1


def mix128(x):
    for _ in range(3):
        x = (x * 0x9E3779B97F4A7C15F39CC0605CEDC835 + 0x0123456789ABCDEF) & M128
        x ^= x >> 67
    return x


class Group:
    """Production transition on a toy GT of order rp over F_q, q^12 tower."""

    def __init__(self, q, rp, Rg, dpbits, depth, seed):
        self.q, self.rp, self.Rg, self.dpbits, self.depth = q, rp, Rg, dpbits, depth
        self.T = Tower(q)
        self.lam = q % rp
        assert pow(self.lam, 12, rp) == 1 and pow(self.lam, 6, rp) == rp - 1
        self.N = (rp - 1) // 12
        self.linv = [pow(self.lam, (-k) % 12, rp) for k in range(12)]
        self.rng = random.Random(0xC0FFEE ^ (seed * 2654435761))

    # ---- instance: fresh target z, fresh g, fresh Rg-entry table (exact g^a z^b)
    def new_instance(self, seed):
        rng = random.Random(0x9E37 ^ (seed * 2246822519))
        e = (self.q ** 12 - 1) // self.rp
        while True:
            x = tuple(tuple((rng.randrange(self.q), rng.randrange(self.q))
                            for _ in range(3)) for _ in range(2))
            g = self.T.pow(x, e)
            if g != self.T.one12() and self.T.pow(g, self.rp) == self.T.one12():
                break
        self.g = g
        self.x = rng.randrange(1, self.rp)
        self.z = self.T.pow(self.g, self.x)
        # table: Rg entries, EXACT g^a z^b, built from two 256-limb sub-tables so the
        # fixture build is affordable; the shipped build price is charged separately.
        nb = 1 << 8
        def sub(base):
            out = []
            for _ in range(nb):
                e = rng.randrange(1, self.rp)
                out.append((e, self.T.pow(base, e)))
            return out
        self.ga = sub(self.g)
        self.zb = sub(self.z)
        self.ent = []
        for h in range(self.Rg):
            i1, i2 = h & (nb - 1), (h >> 8) & (nb - 1)
            j1, j2 = (h * 7 + 13) & (nb - 1), ((h >> 8) * 11 + 5) & (nb - 1)
            a1, M1 = self.ga[i1]
            a2, M2 = self.ga[i2]
            b1, M3 = self.zb[j1]
            b2, M4 = self.zb[j2]
            M = self.T.mul(self.T.mul(M1, M2), self.T.mul(M3, M4))
            self.ent.append((a1 + a2, b1 + b2, M, [None] * 12, 0))
        self.keys_cache = {}

    def rot(self, h, j):
        e = self.ent[h]
        j %= 12
        cache = e[3]
        if cache[j] is None:
            if cache[0] is None:
                cache[0] = e[2]
            cur = cache[j]
            if cur is None:
                k = j
                while cache[k - 1] is None:
                    k -= 1
                cur = cache[k - 1]
                for _ in range(k, j + 1):
                    cur = self.T.frob(cur)
                    cache[_] = cur
            cache[j] = cur
        return cache[j]

    def state(self, y):
        """(pool key, kappa, slot, branch, dp) from ONE keys12 evaluation."""
        ks = self.T.keys12(y)
        m = min(ks)
        kap = ks.index(m)
        v = 0
        for x in sorted(ks):
            v = v * self.q + x
        L = mix128(v)
        return v, kap, (L >> 3) % self.Rg, (L & 7) < 3, \
               ((L >> 3 + 0) & ((1 << self.dpbits) - 1)) == 0 if self.dpbits else True

    def step(self, y, A, D, kap, slot, branch, hist):
        """THE SHIPPED TRANSITION.  hist is updated ONLY on multiply steps."""
        guard = (not branch) and (slot in hist)
        if branch or guard:
            y, A, D = self.T.cyclo_sqr(y), (2 * A) % self.rp, (2 * D) % self.rp
        else:
            a, b, _M, _c, _n = self.ent[slot]
            y = self.T.mul(y, self.rot(slot, (12 - kap) % 12))
            A, D = (A + self.linv[kap] * a) % self.rp, (D + self.linv[kap] * b) % self.rp
            hist = ([slot] + hist)[:self.depth]
        return y, A, D, hist, guard

    def start(self, a, b):
        return self.T.pow(self.g, (a + self.x * b) % self.rp), a % self.rp, b % self.rp

    def shift(self, y1, y2):
        cur = y1
        for k in range(12):
            if cur == y2:
                return k
            cur = self.T.frob(cur)
        return None

    def recover(self, A1, D1, A2, D2, k):
        pj = pow(self.q, k, self.rp)
        den = (D2 - pj * D1) % self.rp
        if den == 0:
            return None
        return ((pj * A1 - A2) % self.rp) * pow(den, -1, self.rp) % self.rp

    def ok(self, xx):
        return xx not in (None, 0) and self.T.pow(self.g, xx) == self.z


# --------------------------------------------------------------------- CP helper
def cp95_lower(k, n, alpha=Decimal('0.05')):
    """EXACT one-sided Clopper-Pearson lower endpoint.  v24 repair: the search
    bracket is (0,1).  v23 searched (0.01,0.99) so 512/512 returned 0.99 although
    the exact endpoint is 0.05^(1/512) = 0.9941..."""
    k = int(k)
    n = int(n)
    if k == 0:
        return Decimal(0)
    if k == n:
        return alpha ** (Decimal(1) / Decimal(n))
    lnC = sum(Decimal(i).ln() for i in range(n - k + 1, n + 1)) - \
        sum(Decimal(i).ln() for i in range(1, k + 1))

    def tail(p):
        p = Decimal(p)
        q = Decimal(1) - p
        t = (lnC + k * p.ln() + (n - k) * q.ln()).exp()
        s = t
        for j in range(k, n):
            t *= (n - j) * p / ((j + 1) * q)
            s += t
            if t < Decimal('1e-45') * s:
                break
        return s

    lo, hi = Decimal(0), Decimal(1)
    for _ in range(140):
        mid = (lo + hi) / 2
        if tail(mid) > alpha:
            hi = mid
        else:
            lo = mid
        if hi - lo < Decimal('1e-15'):
            break
    return lo


# ------------------------------------------------------------------ metered walk
class Meter:
    def __init__(self, cap):
        self.cap = cap
        self.used = 0
        self.rows = {"fixed": 0, "walk": 0, "record": 0, "lookup": 0,
                     "replay": 0, "verify": 0}

    def fits(self, c):
        return self.used + c <= self.cap

    def charge(self, row, c):
        self.used += c
        self.rows[row] += c


def run_one(arg):
    gid, cfg = arg
    """One fresh fixture: fresh target, fresh table, fresh coins.  ONE Bernoulli.

    Meter units are INSTRUCTIONS.  Size-dependent components are charged at their
    true prices (walk, lookup, record, replay).  Group-size-independent components
    (table build, 1024 worker starts, code) are charged as their production
    PROPORTIONAL SHARE, cfg["fixed_frac"] = (table+starts+code)/W_prod: at fixture
    scale they are 2^99 times the whole fixture budget, so charging them at face
    value would make the fixture a measurement of the table build (fidelity-table
    row).  Verification is charged once in the ledger on success.
    """
    G = Group(cfg["q"], cfg["rp"], cfg["Rg"], cfg["dpbits"], cfg["depth"], 7)
    G.new_instance(cfg["seed0"] + gid)
    rng = random.Random(0x51CE ^ (gid * 2654435761))
    trail = cfg["trail"]
    nw = cfg["workers"]
    cap = cfg["cap_ir"]
    m = Meter(cap)
    m.charge("fixed", int(cap * cfg["fixed_frac"]))
    pool = {}
    steps = recs = matches = 0
    solved = False
    for wid in range(nw):
        a, b = rng.randrange(1, G.rp), rng.randrange(1, G.rp)
        y, A, D = G.start(a, b)
        hist = []
        key, kap, slot, br, dp = G.state(y)
        for _k in range(trail):
            c = MUL_CHARGE if not br else SQR_CHARGE
            if not m.fits(c + LOOKUP + REC_WRITE):
                break
            m.charge("walk", c)
            m.charge("lookup", LOOKUP)
            y, A, D, hist, _g = G.step(y, A, D, kap, slot, br, hist)
            steps += 1
            key, kap, slot, br, dp = G.state(y)
            if dp:
                m.charge("record", REC_WRITE)
                recs += 1
                old = pool.get(key)
                pool[key] = (A, D, y)
                if old is not None:
                    matches += 1
                    sh = G.shift(old[2], y)
                    if sh is not None:
                        rp_cost = 2 * (steps // max(1, nw)) * MUL_CHARGE
                        if m.fits(rp_cost):
                            m.charge("replay", rp_cost)
                        xx = G.recover(old[0], old[1], A, D, sh)
                        if xx is not None and G.ok(xx):
                            solved = True
                            break
        if solved:
            break
    return {"gid": gid, "solved": solved, "steps": steps, "records": recs,
            "matches": matches, "ir_used": m.used, "rows": m.rows, "cap_ir": cap,
            "exhausted": (not solved) and (m.used + MUL_CHARGE > cap)}


def mode_census(cfg):
    from multiprocessing import Pool
    t0 = time.time()
    out = {"mode": "census", "config": {k: v for k, v in cfg.items() if k != "seed0"},
           "trials_per_fixture": 1,
           "independence": "fresh target z, fresh g, fresh table, fresh coins per "
                           "fixture; one Bernoulli per fixture",
           "scales": []}
    nfix = cfg["fixtures"]
    with Pool(cfg["procs"]) as p:
        for sc in cfg["groups"]:
            gcfg = dict(cfg, **sc)
            N = (gcfg["rp"] - 1) // 12
            steps_budget = int(gcfg["ratio"] * math.sqrt(N))
            gcfg["trail"] = max(1, -(-steps_budget // gcfg["workers"]))
            gcfg["cap_ir"] = steps_budget * MUL_CHARGE
            res = p.map(run_one, [(i, gcfg) for i in range(nfix)], chunksize=4)
            ok = sum(1 for r in res if r["solved"])
            steps = sum(r["steps"] for r in res)
            cp = cp95_lower(ok, nfix)
            out["scales"].append({
                "q": gcfg["q"], "rp": gcfg["rp"], "N": N,
                "Rg": gcfg["Rg"], "depth": gcfg["depth"], "dpbits": gcfg["dpbits"],
                "workers": gcfg["workers"], "trail_cap": gcfg["trail"],
                "sqrtN": round(math.sqrt(N), 2),
                "steps_budget": steps_budget,
                "S_over_sqrtN": round(steps_budget / math.sqrt(N), 4),
                "fixtures": nfix, "success": ok, "rate": ok / nfix,
                "cp95_lower": str(cp),
                "steps_executed": steps, "records": sum(r["records"] for r in res),
                "matches": sum(r["matches"] for r in res),
                "mean_ir_used": sum(r["ir_used"] for r in res) // nfix,
                "exhausted": sum(1 for r in res if r["exhausted"]),
                "rows_ir": {k: sum(r["rows"][k] for r in res) for k in res[0]["rows"]}})
            print("scale", gcfg["rp"], ok, "/", nfix, "cp95", cp, flush=True)
    out["seconds"] = round(time.time() - t0, 1)
    return out


# --------------------------------------------------------------- survival (merge)
def survival_arm(args):
    rp, q, Rg, arm, pairs, gaps, depth = args
    G = Group(q, rp, Rg, 0, depth, 3)
    G.new_instance(11)
    rng = random.Random(0xB22 + Rg * 131 + len(arm))
    same = {g: 0 for g in gaps}
    by_j = {g: {} for g in gaps}
    gmax = max(gaps)
    for _ in range(pairs):
        a, b = rng.randrange(1, G.rp), rng.randrange(1, G.rp)
        j = rng.randrange(0, 12)
        y1, A1, D1 = G.start(a, b)
        y2 = y1
        for _i in range(j):
            y2 = G.T.frob(y2)
        h1 = [rng.randrange(Rg) for _ in range(depth)]
        h2 = list(h1) if arm == "true_merge" else [rng.randrange(Rg) for _ in range(depth)]
        A2, D2 = A1, D1
        k1 = G.state(y1)
        k2 = G.state(y2)
        for i in range(gmax + 1):
            if i in same:
                eq = (y1 == y2 if arm == "samerep" else G.shift(y1, y2) is not None)
                if eq:
                    same[i] += 1
                    by_j[i][j] = by_j[i].get(j, 0) + 1
            if i < gmax:
                y1, A1, D1, h1, _ = G.step(y1, A1, D1, k1[1], k1[2], k1[3], h1)
                y2, A2, D2, h2, _ = G.step(y2, A2, D2, k2[1], k2[2], k2[3], h2)
                k1 = G.state(y1)
                k2 = G.state(y2)
    return {"rp": rp, "Rg": Rg, "arm": arm, "pairs": pairs, "depth": depth,
            "survival": {str(g): same[g] / pairs for g in gaps},
            "survival_by_kappa_offset": {str(g): {str(k): v for k, v in sorted(by_j[g].items())}
                                         for g in gaps}}


def mode_survival(cfg):
    from multiprocessing import Pool
    groups = cfg["groups"]
    gaps = cfg["gaps"]
    arms = cfg.get("arms", ["true_merge", "indep_hist", "samerep"])
    jobs = [(g["rp"], g["q"], g["Rg"], a, cfg["pairs"], gaps, cfg["depth"])
            for g in groups for a in arms]
    with Pool(cfg["procs"]) as p:
        rs = p.map(survival_arm, jobs)
    out = {"mode": "survival", "pairs": cfg["pairs"], "gaps": gaps,
           "depth": cfg["depth"], "arms": {}}
    for r in rs:
        out["arms"]["%s_rp%d_Rg%d" % (r["arm"], r["rp"], r["Rg"])] = r
    # fit survival(gap) = rho^gap on the true_merge arm (log-linear, all gaps)
    fits = {}
    for k, r in out["arms"].items():
        sv = r["survival"]
        xs = [g for g in gaps if g >= 1 and sv[str(g)] > 0]
        if len(xs) >= 3:
            ys = [math.log(sv[str(g)]) for g in xs]
            n = len(xs)
            mx = sum(xs) / n
            my = sum(ys) / n
            b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / \
                sum((x - mx) ** 2 for x in xs)
            rho = math.exp(b)
            resid = [ys[i] - (my + b * (xs[i] - mx)) for i in range(n)]
            band = max(abs(e) for e in resid)
            fits[k] = {"rho": rho, "log_rho": b, "max_abs_residual_log": band,
                       "band_factor": math.exp(band), "n_points": n}
    out["fit"] = fits
    return out


# --------------------------------------------------------------------- scan (real)
def scan_segment(args):
    """Observer-only: the walk is advanced ONLY by step_msg (imported from this
    module's own transition), so history updates exactly as shipped."""
    wid, Rg, depth, nsteps, q, rp = args
    G = Group(q, rp, Rg, 0, depth, 5)
    G.new_instance(4400 + wid)
    rng = random.Random(0xA5C1 ^ (wid * 7919))
    a, b = rng.randrange(1, G.rp), rng.randrange(1, G.rp)
    y, A, D = G.start(a, b)
    hist = []
    key, kap, slot, br, dp = G.state(y)
    mult = brsq = guard = 0
    W = (4, 6, 8)
    adm = {m: 0 for m in W}
    alt = {m: 0 for m in W}
    delays = []
    slots = {}
    kaps = [0] * 12
    win = []
    stall = 0
    K_STALL = 32
    pending = None
    for _s in range(nsteps):
        y0, A0, D0, kap0, slot0, br0 = y, A, D, kap, slot, br
        y, A, D, hist, g = G.step(y, A, D, kap, slot, br, hist)
        if br0:
            brsq += 1
        elif g:
            guard += 1
        else:
            mult += 1
            slots[slot0] = slots.get(slot0, 0) + 1
            win = (win + [(slot0, kap0)])[-max(W):]
            for m in W:
                if len(win) >= m:
                    w = win[-m:]
                    if all(w[i][0] == w[i - 2][0] for i in range(2, m)):
                        alt[m] += 1
                        per = {}
                        for s, k in w:
                            per.setdefault(s, [0] * 12)[k] += 1
                        if all(reduce_counts(v) == [0, 0, 0, 0] for v in per.values()):
                            adm[m] += 1
                            if pending is None:
                                pending = 0
        # stall detector: K_STALL consecutive squaring steps -> trail abandoned
        stall = stall + 1 if (br0 or g) else 0
        if pending is not None:
            pending += 1
            if stall >= K_STALL:
                delays.append(pending)
                pending = None
        if stall >= 999999:
            break
        key, kap, slot, br, dp = G.state(y)
        kaps[kap] += 1
    return {"steps": nsteps, "multiply": mult, "branch_sqr": brsq, "guard": guard,
            "alt": alt, "adm": adm, "nslots": len(slots), "kaps": kaps,
            "delays": delays[:4000], "n_delays": len(delays)}


def mode_scan(cfg):
    from multiprocessing import Pool
    segs, nsteps = cfg["segs"], cfg["steps"]
    jobs = [(i, cfg["Rg"], cfg["depth"], nsteps, cfg["q"], cfg["rp"]) for i in range(segs)]
    with Pool(cfg["procs"]) as p:
        rs = p.map(scan_segment, jobs)
    num = lambda k: sum(r[k] for r in rs)
    dd = lambda k: {str(m): sum(r[k][m] for r in rs) for m in (4, 6, 8)}
    T = num("steps")
    kaps = [sum(r["kaps"][i] for r in rs) for i in range(12)]
    tot = sum(kaps)
    delays = [d for r in rs for d in r["delays"]]
    return {"mode": "scan", "Rg": cfg["Rg"], "depth": cfg["depth"], "q": cfg["q"],
            "rp": cfg["rp"], "segments": segs, "steps_per_segment": nsteps, "steps": T,
            "multiply": num("multiply"), "multiply_per_step": num("multiply") / T,
            "branch_sqr": num("branch_sqr"), "guard_conversions": num("guard"),
            "guard_per_step": num("guard") / T,
            "guard_per_multiply": num("guard") / num("multiply"),
            "word_alternating": dd("alt"), "word_admit": dd("adm"),
            "admit_per_step": {k: v / T for k, v in dd("adm").items()},
            "abandon_delays": delays, "abandon_events": len(delays),
            "abandon_delay_mean": (sum(delays) / len(delays)) if delays else None,
            "abandon_delay_max": max(delays) if delays else None,
            "distinct_slots_per_segment": num("nslots") / segs,
            "kappa_histogram": kaps,
            "kappa_max_rel_dev": max(abs(k * 12 - tot) / tot for k in kaps if tot)}


if __name__ == "__main__":
    cfg = json.loads(sys.argv[1] if len(sys.argv) > 1 else "{}")
    mode = cfg.pop("mode")
    if mode == "cpctl":
        k, n = cfg.get("k", 512), cfg.get("n", 512)
        got = cp95_lower(k, n)
        want = Decimal(0.05) ** (Decimal(1) / Decimal(n)) if k == n else None
        out = {"mode": "cpctl", "k": k, "n": n, "cp95_lower": str(got),
               "exact_endpoint": str(want) if want is not None else None,
               "matches_exact": (want is not None and abs(got - want) < Decimal('1e-18')),
               "v23_clipped_value": "0.99", "clipped_error": str(Decimal('0.99') - want)
               if want is not None else None}
        print(json.dumps(out, indent=1))
    else:
        out = {"census": mode_census, "survival": mode_survival,
               "scan": mode_scan}[mode](cfg)
    ev = cfg.get("out")
    if ev:
        p = os.path.join(os.path.dirname(HERE), "evidence", ev)
        json.dump(out, open(p, "w"), indent=1)
        print("wrote", p)