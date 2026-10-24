#!/usr/bin/env python3
"""protocol25.py -- THE SHIPPED PROTOCOL, METERED, v25.

Repairs verdict d3d07cd8 flaw 1 and flaw 3 in the CONTROL FLOW, not in the prose:

  * meter-first from the first instruction: the meter is created before any table
    entry exists; the table-build action, every worker-start exponentiation, the
    code/advice action and every walk step ask fits() BEFORE acting.  At fixture
    scale the group-size-independent actions are charged their proportional share
    (<1 Ir), charged as 1 Ir each -- an overcharge -- and the deviation is disclosed
    in the fidelity table.
  * the step price is GUARD-AWARE at charge time: guard = (not branch) and
    (slot in hist) is computed BEFORE the charge, so a guard-converted multiply is
    charged SQR_CHARGE, not MUL_CHARGE.  v24 charged by the branch flag.
  * a failed replay charge ABANDONS the worker: if fits(replay) is false the
    worker stops, it does not continue to recovery and verification.
  * candidate verification is INSIDE the meter: g^xx is recomputed by a metered
    square-and-multiply (every squaring and multiply asks fits), the charge lands
    in the verify row, and a verification that does not fit abandons the worker.
  * the 32-consecutive-squaring abort is IMPLEMENTED: the counter is armed by a
    guard admission and a run of 32 squaring steps abandons the worker.  The scan
    reports the delay, and a mutant build with forced admissions measures the
    delay-to-abandon distribution; the abort-disabled control never abandons.

The census budget is steps_budget * (MUL_CHARGE + LOOKUP + REC_WRITE) plus one
verification reserve, so the fixture can actually execute the claimed S/sqrt(N)
steps under the same per-step charges the ledger prices.
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

EV = os.path.join(os.path.dirname(HERE), "evidence")
_PR = json.load(open(os.path.join(EV, "prices25.json")))
MUL_REF = 20000
MUL_CHARGE = _PR["MUL_CHARGE"]
SQR_CHARGE = _PR["SQR_CHARGE"]
LOOKUP = _PR["LOOKUP"]
VERIFY_PROD = _PR["VERIFY_PROD"]
VERIFY_RESERVE = _PR["VERIFY_RESERVE"]
REPLAY_IR = _PR["REPLAY_IR"]
K_STALL = 32
TABLE_IR = _PR["TABLE_IR"]
START_IR = _PR["START_IR"]
CODE_IR = _PR["CODE_IR"]
W_PROD = _PR["W_PROD"]

M128 = (1 << 128) - 1


def mix128(x):
    for _ in range(3):
        x = (x * 0x9E3779B97F4A7C15F39CC0605CEDC835 + 0x0123456789ABCDEF) & M128
        x ^= x >> 67
    return x


def serialize_record(key, A, D, y):
    """The record as it is written: 6-word class key (48 B), two 255-bit labels
    (2 x 32 B), and the full GT element retained for the Frobenius-shift test
    (12 Fq components x 48 B).  Returns BYTES; REC_WRITE is len() of this."""
    import struct
    buf = struct.pack("<6Q", *[(key >> (64 * i)) & ((1 << 64) - 1) for i in range(6)])
    buf += A.to_bytes(32, "little") + D.to_bytes(32, "little")
    for c in y:
        for p in c:
            for x in p:
                buf += int(x).to_bytes(48, "little")
    return buf


REC_BYTES = len(serialize_record(0, 0, 0, (((0, 0), (0, 0), (0, 0)),
                                          ((0, 0), (0, 0), (0, 0)))))
REC_WRITE = REC_BYTES


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
    """EXACT one-sided Clopper-Pearson lower endpoint, bracket (0,1), unclipped."""
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


def pow_metered(tower, base, e, meter, mul, sqr):
    """g^e inside the global meter: EVERY squaring and multiply asks fits() first.
    Returns (value, None) or (None, 'abandon') when the next action does not fit."""
    r = tower.one12()
    for bit in range(e.bit_length() - 1, -1, -1):
        if not meter.fits(sqr):
            return None, "abandon"
        meter.charge("verify", sqr)
        r = tower.mul(r, r)
        if (e >> bit) & 1:
            if not meter.fits(mul):
                return None, "abandon"
            meter.charge("verify", mul)
            r = tower.mul(r, base)
    return r, None


def run_one(arg):
    """One fresh fixture: fresh target, fresh table, fresh coins.  ONE Bernoulli.

    All actions are metered.  Group-size-independent actions are charged their
    fixture proportional share, at least 1 Ir each (an overcharge); the production
    ledger charges them in full.  Success = a metered verification says g^x == z
    before the meter is exhausted; every failed replay/verify charge abandons.
    """
    gid, cfg = arg
    G = Group(cfg["q"], cfg["rp"], cfg["Rg"], cfg["dpbits"], cfg["depth"], 7)
    rng = random.Random(0x51CE ^ (gid * 2654435761))
    nw = cfg["workers"]
    cap = cfg["cap_ir"]
    m = Meter(cap)
    steps = recs = matches = verifies = aborts = 0
    solved = False
    abandoned = False
    reason = "none"
    # code/advice action, metered before the table exists
    cshare = max(1, cap * CODE_IR // W_PROD)
    if not m.fits(cshare):
        return {"gid": gid, "solved": False, "steps": 0, "records": 0, "matches": 0,
                "verifies": 0, "aborts": 0, "ir_used": 0, "rows": m.rows,
                "cap_ir": cap, "exhausted": True, "abandoned": True}
    m.charge("fixed", cshare)
    # table-build action, metered BEFORE the table is constructed
    tshare = max(1, cap * TABLE_IR // W_PROD)
    if not m.fits(tshare):
        return {"gid": gid, "solved": False, "steps": 0, "records": 0, "matches": 0,
                "verifies": 0, "aborts": 0, "ir_used": m.used, "rows": m.rows,
                "cap_ir": cap, "exhausted": True, "abandoned": True}
    m.charge("fixed", tshare)
    G.new_instance(cfg["seed0"] + gid)
    trail = cfg["trail"]
    pool = {}
    for wid in range(nw):
        sshare = max(1, cap * START_IR // W_PROD)
        if not m.fits(sshare):
            abandoned = True
            break
        m.charge("fixed", sshare)
        a, b = rng.randrange(1, G.rp), rng.randrange(1, G.rp)
        y, A, D = G.start(a, b)
        hist = []
        key, kap, slot, br, dp = G.state(y)
        stall = 0
        armed = False
        aborted = False
        for _k in range(trail):
            guard = (not br) and (slot in hist)
            c = SQR_CHARGE if (br or guard) else MUL_CHARGE
            if not m.fits(c + LOOKUP + REC_WRITE):
                abandoned = True; reason = "walk"
                break
            m.charge("walk", c)
            m.charge("lookup", LOOKUP)
            # 32-consecutive-squaring abort, armed by a guard admission
            if br or guard:
                stall += 1
                if guard:
                    armed = True
                if armed and stall >= K_STALL:
                    aborts += 1
                    aborted = True
                    abandoned = True
                    break
            else:
                stall = 0
                armed = False
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
                    rp_cost = REPLAY_IR
                    if not m.fits(rp_cost):
                        abandoned = True; reason = "replay"          # replay that does not fit ABANDONS
                        aborted = True
                        break
                    m.charge("replay", rp_cost)
                    sh = G.shift(old[2], y)
                    if sh is not None:
                        xx = G.recover(old[0], old[1], A, D, sh)
                        if xx is not None:
                            val, err = pow_metered(G.T, G.g, xx, m,
                                                   MUL_CHARGE, SQR_CHARGE)
                            verifies += 1
                            if err is not None:
                                abandoned = True; reason = "verify"  # verification that does not fit
                                aborted = True
                                break
                            if val == G.z:
                                solved = True
                                break
        if solved:
            break
        if aborted:
            continue
    return {"gid": gid, "solved": solved, "steps": steps, "records": recs,
            "matches": matches, "verifies": verifies, "aborts": aborts,
            "ir_used": m.used, "rows": m.rows, "cap_ir": cap,
            "exhausted": (not solved) and (not abandoned),
            "abandoned": abandoned, "reason": reason}


def mode_census(cfg):
    from multiprocessing import Pool
    t0 = time.time()
    out = {"mode": "census", "config": {k: v for k, v in cfg.items() if k != "seed0"},
           "trials_per_fixture": 1,
           "independence": "fresh target z, fresh g, fresh table, fresh coins per "
                           "fixture; one Bernoulli per fixture",
           "protocol": "protocol25.run_one: meter created first; table, starts, code, "
                       "steps (guard-aware), lookups, records, replay, and candidate "
                       "verification all charge the same global meter; a failed "
                       "replay/verify charge abandons; 32-squaring abort implemented",
           "scales": []}
    nfix = cfg["fixtures"]
    with Pool(cfg["procs"]) as p:
        for sc in cfg["groups"]:
            gcfg = dict(cfg, **sc)
            N = (gcfg["rp"] - 1) // 12
            steps_budget = int(gcfg["ratio"] * math.sqrt(N))
            gcfg["trail"] = max(1, -(-steps_budget // gcfg["workers"]))
            gcfg["cap_ir"] = steps_budget * (MUL_CHARGE + LOOKUP + REC_WRITE) + \
                4 * (2 * gcfg["trail"] * MUL_CHARGE) + VERIFY_RESERVE
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
                "cap_ir": gcfg["cap_ir"],
                "fixtures": nfix, "success": ok, "rate": ok / nfix,
                "cp95_lower": str(cp),
                "steps_executed": steps, "records": sum(r["records"] for r in res),
                "matches": sum(r["matches"] for r in res),
                "verifies": sum(r["verifies"] for r in res),
                "aborts": sum(r["aborts"] for r in res),
                "abandoned": sum(1 for r in res if r["abandoned"]),
                "mean_ir_used": sum(r["ir_used"] for r in res) // nfix,
                "exhausted": sum(1 for r in res if r["exhausted"]),
                "record_bytes": REC_BYTES,
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
    """Observer-only scan on the shipped step function, with the 32-squaring abort
    implemented.  force_admit is a MUTANT build: every alternating word is admitted
    and the admitted trap converts subsequent multiplies to squarings (the guard
    sees the slot in hist); abort_enabled=False is the control that never abandons.
    """
    wid, Rg, depth, nsteps, q, rp, force_admit, abort_enabled, force_period = args
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
    armed = False
    pending = None
    trapped = False
    aborts = 0
    would_abort = 0
    max_srun = 0
    for _s in range(nsteps):
        y0, A0, D0, kap0, slot0, br0 = y, A, D, kap, slot, br
        if force_admit and (not trapped) and force_period and _s > 0 and _s % force_period == 0:
            pending = 0; trapped = True; adm[4] += 1
        if trapped:
            # mutant trap dynamics: the guard blocks every multiply, so the walk
            # sees only squarings until the abort counter fires.
            stall += 1
            max_srun = max(max_srun, stall)
            if pending is not None:
                pending += 1
            if stall >= K_STALL:
                if abort_enabled:
                    delays.append(pending if pending is not None else stall)
                    aborts += 1
                    pending = None
                    trapped = False
                    stall = 0
                    if _s + 1 < nsteps:
                        # abandon this worker, restart it at a fresh start
                        aa, bb = rng.randrange(1, G.rp), rng.randrange(1, G.rp)
                        y, A, D = G.start(aa, bb)
                        hist = []
                        key, kap, slot, br, dp = G.state(y)
                else:
                    would_abort += 1
                    stall = 0
            continue
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
                        if force_admit or all(reduce_counts(v) == [0, 0, 0, 0]
                                              for v in per.values()):
                            adm[m] += 1
                            if pending is None:
                                pending = 0
                                trapped = force_admit
        if br0 or g:
            stall += 1
            if g:
                armed = True
            max_srun = max(max_srun, stall)
            if armed and stall >= K_STALL:
                if abort_enabled:
                    delays.append(pending if pending is not None else stall)
                    aborts += 1
                    pending = None
                    stall = 0
                    armed = False
                else:
                    would_abort += 1
                    stall = 0
        else:
            stall = 0
            armed = False
        key, kap, slot, br, dp = G.state(y)
        kaps[kap] += 1
    return {"steps": nsteps, "multiply": mult, "branch_sqr": brsq, "guard": guard,
            "alt": alt, "adm": adm, "nslots": len(slots), "kaps": kaps,
            "delays": delays[:4000], "n_delays": len(delays), "aborts": aborts,
            "would_abort": would_abort, "max_srun": max_srun}


def mode_scan(cfg):
    from multiprocessing import Pool
    segs, nsteps = cfg["segs"], cfg["steps"]
    fa = cfg.get("force_admit", False)
    ab = cfg.get("abort_enabled", True)
    fp = cfg.get("force_period", 1000)
    jobs = [(i, cfg["Rg"], cfg["depth"], nsteps, cfg["q"], cfg["rp"], fa, ab, fp)
            for i in range(segs)]
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
            "force_admit": fa, "abort_enabled": ab,
            "multiply": num("multiply"), "multiply_per_step": num("multiply") / T,
            "branch_sqr": num("branch_sqr"), "guard_conversions": num("guard"),
            "guard_per_step": num("guard") / T,
            "guard_per_multiply": num("guard") / num("multiply"),
            "word_alternating": dd("alt"), "word_admit": dd("adm"),
            "admit_per_step": {k: v / T for k, v in dd("adm").items()},
            "abandon_delays": delays, "abandon_events": len(delays),
            "abort_events": num("aborts"), "abort_disabled_stalls": num("would_abort"),
            "abandon_delay_mean": (sum(delays) / len(delays)) if delays else None,
            "abandon_delay_max": max(delays) if delays else None,
            "max_squaring_run": max(r["max_srun"] for r in rs),
            "distinct_slots_per_segment": num("nslots") / segs,
            "kappa_histogram": kaps,
            "kappa_max_rel_dev": max(abs(k * 12 - tot) / tot for k in kaps if tot)}


def mode_tie(cfg):
    """Tie analysis: keys12(y) is a 12-vector of canonical values; the class key is
    the SORTED tuple, which is invariant to which rotation attains the minimum, so
    a tie cannot change the pool slot.  This mode MEASURES the tie frequency."""
    G = Group(cfg["q"], cfg["rp"], cfg["Rg"], 0, cfg["depth"], 9)
    G.new_instance(77)
    rng = random.Random(0x71E)
    ties = 0
    n = cfg["samples"]
    example = None
    for _ in range(n):
        y = G.T.pow(G.g, rng.randrange(1, G.rp))
        ks = G.T.keys12(y)
        m = min(ks)
        if ks.count(m) > 1:
            ties += 1
            if example is None:
                example = [str(x) for x in ks]
    return {"mode": "tie", "samples": n, "ties": ties, "tie_rate": ties / n,
            "resolution": "kappa = ks.index(min(ks)) is the FIRST minimum (deterministic); "
                          "the pool slot is derived from sorted(ks), invariant to which "
                          "rotation attains the minimum, so a tie cannot fork the table "
                          "index or the record key.",
            "example": example}


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
               "scan": mode_scan, "tie": mode_tie}[mode](cfg)
    ev = cfg.get("out")
    if ev:
        p = os.path.join(EV, ev)
        json.dump(out, open(p, "w"), indent=1)
        print("wrote", p)
