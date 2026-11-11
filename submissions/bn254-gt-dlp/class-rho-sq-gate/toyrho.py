"""Toy-field run of proof.md's complete algorithm (imported by experiment.py, mode toy).
Python standard library only; evidence at toy scale only."""
import math
import random

from experiment import Tower, is_prime, khash, n_eff, nu_bound


def find_toy(rbits, qbits):
    rp = (1 << rbits) - ((1 << rbits) % 12) + 1
    while True:
        rp += 12
        if is_prime(rp): break
    for h in range(2, rp):
        lam = pow(h, (rp - 1) // 12, rp)
        if pow(lam, 6, rp) != 1 and pow(lam, 4, rp) != 1: break
    m = 12 * rp
    q0 = next(lam + k * rp for k in range(12) if (lam + k * rp) % 12 == 7)
    q = q0 + (((1 << qbits) - q0) // m + 1) * m
    while not is_prime(q):
        q += m
    for a0 in range(1, 1000):
        T = Tower.__new__(Tower); T.q = q
        xi = (a0, 1)
        if T.f2pow(xi, (q * q - 1) // 2) != (1, 0) and T.f2pow(xi, (q * q - 1) // 3) != (1, 0): break
    return q, rp, a0


def bn_toy(x):
    """BN-family toy: q = p(x), r' = r(x) from the BN polynomials (x odd and negative here), so q ~ r'
    as in production; xi = a0 + u with a0 the least value making xi neither a square nor a cube."""
    q = 36 * x ** 4 + 36 * x ** 3 + 24 * x ** 2 + 6 * x + 1
    rp = 36 * x ** 4 + 36 * x ** 3 + 18 * x ** 2 + 6 * x + 1
    assert is_prime(q) and is_prime(rp) and q % 12 == 7 and all(pow(q, k, rp) != 1 for k in range(1, 12))
    for a0 in range(1, 1000):
        T = Tower.__new__(Tower); T.q = q
        xi = (a0, 1)
        if T.f2pow(xi, (q * q - 1) // 2) != (1, 0) and T.f2pow(xi, (q * q - 1) // 3) != (1, 0): break
    return q, rp, a0


class ClassRho:
    """proof.md's algorithm with toy sizes: R indices, window w, A probes, dbits distinguished bits,
    trail cap 2^(dbits+7) slots. Runs trail after trail until it outputs x or has used `cap` slots."""

    def __init__(self, T, rp, g, R_idx, dbits, rng, w=24, A=8, sq=True):
        self.T, self.rp, self.g, self.R, self.dbits, self.rng, self.w, self.A = T, rp, g, R_idx, dbits, rng, w, A
        self.rb = R_idx.bit_length() - 1
        self.lam = T.q % rp
        self.sq = sq                                   # squaring branch on (bit 63 of hash word 2)

    def setup(self, z):
        T, rp = self.T, self.rp
        self.z, self.s1, self.s2 = z, self.rng.getrandbits(64), self.rng.getrandbits(64)
        li = pow(self.lam, -1, rp)
        self.TW, self.EX = [], []
        for h in range(self.R):
            u, v = self.rng.randrange(rp), self.rng.randrange(rp)
            Th = T.mul(T.pow(self.g, u), T.pow(z, v))
            self.TW.append([T.frob_k(Th, (12 - k) % 12) for k in range(12)])
            self.EX.append([(u * pow(li, k, rp) % rp, v * pow(li, k, rp) % rp) for k in range(12)])

    def run(self, cap):
        T, rp, lam = self.T, self.rp, self.lam
        tcap = 1 << (self.dbits + 7)
        store = {}
        S = slots = starts = rejects = kills = ties = capped = sqsteps = 0
        res = lambda x: dict(x=x, slots=slots, rejects=rejects, kills=kills, ties=ties, starts=starts,
                             capped=capped, steps=S, sq_steps=sqsteps)
        while slots < cap:
            starts += 1
            al, be = self.rng.randrange(rp), self.rng.randrange(rp)
            y = T.mul(T.pow(self.g, al), T.pow(self.z, be))
            sel = T.select(y)
            if sel is None:
                ties += 1; continue
            K, ks = sel
            win, n = [K], 0
            while True:
                w1 = khash(K, self.s1)
                h, dp = w1 & (self.R - 1), ((w1 >> self.rb) & ((1 << self.dbits) - 1)) == 0
                sqb = self.sq and (khash(K, self.s2) >> 63) == 1
                if dp:
                    f = pow(lam, ks, rp)
                    ac, bc = al * f % rp, be * f % rp
                    if K in store:
                        a2, b2 = store[K]
                        if b2 != bc:
                            x = (ac - a2) * pow((b2 - bc) % rp, -1, rp) % rp
                            if T.pow(self.g, x) == self.z:
                                return res(x)
                    else:
                        store[K] = (ac, bc)
                    break
                if n >= tcap:
                    capped += 1
                    break
                ok = False
                for a in range(self.A):
                    hh = (h + a) % self.R
                    square = a == 0 and sqb                 # probe 0 squares on a squaring class
                    if slots >= cap:                    # run slot cap: the run fails
                        return res(None)
                    y2 = T.mul(y, y) if square else T.mul(y, self.TW[hh][ks]); slots += 1; n += 1
                    s2 = T.select(y2)
                    if s2 is None:
                        ties += 1; break
                    if s2[0] in win:
                        rejects += 1; continue
                    ua, va = (al, be) if square else self.EX[hh][ks]
                    sqsteps += square
                    al, be, y, (K, ks) = (al + ua) % rp, (be + va) % rp, y2, s2
                    win.append(K)
                    if len(win) > self.w: win.pop(0)
                    S += 1; ok = True
                    break
                if not ok:
                    kills += 1; break
        return res(None)


def mean_slots(N, R_idx):
    """Model mean of the new states before the first class repeat: sqrt(pi N'/2)."""
    return math.sqrt(math.pi * float(n_eff(N, R_idx, True)) / 2)


def toy(seed, instances, R_idx, family="bn", rbits=22, qbits=40, bn_x=-15, dbits=4):
    q, rp, a0 = bn_toy(bn_x) if family == "bn" else find_toy(rbits, qbits)
    T = Tower(q, a0)
    rng = random.Random(10 ** 6)
    while True:
        y = tuple(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(3)) for _ in range(2))
        g = T.pow(y, (q ** 12 - 1) // rp)
        if g != T.one: break
    N = (rp - 1) // 12
    m0 = mean_slots(N, R_idx)
    cap = 64 * (math.isqrt(int(n_eff(N, R_idx, True))) + 1)          # toy slot cap: 64 sqrt(N')
    rows = []
    for i in range(instances):
        rr = random.Random((seed << 20) + (R_idx << 8) + i)
        x = rr.randrange(1, rp)
        alg = ClassRho(T, rp, g, R_idx, dbits, rr)
        alg.setup(T.pow(g, x))
        res = alg.run(cap)
        res["correct"] = res["x"] == x if res["x"] is not None else None
        rows.append(res)
    ratio = [o["slots"] / m0 for o in rows]
    mean = sum(ratio) / instances
    se = math.sqrt(sum((v - mean) ** 2 for v in ratio) / (instances - 1) / instances)
    slots = sum(o["slots"] for o in rows)
    return {"family": family, "q": q, "r_prime": rp, "xi": f"{a0}+u", "classes": N, "R": R_idx, "window": 24, "probes": 8,
            "dp_bits": dbits, "trail_cap": 1 << (dbits + 7), "slot_cap": cap, "instances": instances, "seed": seed,
            "model_mean_slots": round(m0, 1),
            "correct_outputs": sum(1 for o in rows if o["correct"]),
            "wrong_outputs": sum(1 for o in rows if o["correct"] is False),
            "runs_at_slot_cap": sum(1 for o in rows if o["x"] is None),
            "mean_slots_over_model": round(mean, 4), "standard_error": round(se, 4),
            "squaring_share_of_slots": round(sum(o["sq_steps"] for o in rows) / slots, 4),
            "reject_rate": round(sum(o["rejects"] for o in rows) / slots, 6),
            "nu_bound": round(nu_bound(R_idx, N), 6),
            "killed_trails": sum(o["kills"] for o in rows), "ties": sum(o["ties"] for o in rows),
            "trails_at_cap": sum(o["capped"] for o in rows),
            "note": "at 8 to 48 instances the mean has a standard error of 0.07 to 0.2; the offline runs are the evidence"}
