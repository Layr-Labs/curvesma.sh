#!/usr/bin/env python3
"""Special-q ball relation counts, toy tower-NFS shape (stdlib only).
h = y^2 - H1 y - H0, w unit, f = x^4 + C1 w x^2 + C0 w^2, g = x^2 - U w, cell (a, b) <-> a - b x.
Q = (q, y - s, x - t), q in (B, 8B]; L_Q = {a(s) = t b(s) mod q}, LLL, ball c in [-R, R]^4 up to
sign; RELATION = admissible cell whose norms (q removed on Q's side) are both B-smooth.
run(seed, balls=3, radius=None) | CLI: table TRIES SEED | part K SECONDS | merge"""
import json, math, os, random, sys, time
from math import gcd, log

H1, H0, W0, W1, C1, C0 = 1, 1, 0, 1, 1, 1
U = 62
B = 1 << 16
R = 6
HERE = os.path.dirname(os.path.abspath(__file__))
HD, PD = H1 * H1 + 4 * H0, C1 * C1 - 4 * C0

def is_prime(n):
    if n < 2:
        return False
    SP = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for q in SP:
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in SP:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def primes_upto(n):
    sv = bytearray([1]) * (n + 1)
    sv[0:2] = b"\0\0"
    for i in range(2, int(n ** 0.5) + 1):
        if sv[i]:
            sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return [i for i in range(n + 1) if sv[i]]

def sqrt_mod(a, q):
    a %= q
    if a == 0:
        return 0
    if pow(a, (q - 1) // 2, q) != 1:
        return None
    s, d = 0, q - 1
    while d % 2 == 0:
        d //= 2
        s += 1
    z = 2
    while pow(z, (q - 1) // 2, q) != q - 1:
        z += 1
    m, c, t, r = s, pow(z, d, q), pow(a, d, q), pow(a, (d + 1) // 2, q)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % q
            i += 1
        b = pow(c, 1 << (m - i - 1), q)
        m, c = i, b * b % q
        t, r = t * c % q, r * b % q
    return r

def quad_roots(c1, c0, q):
    sd = sqrt_mod(c1 * c1 - 4 * c0, q)
    if sd is None:
        return []
    i2 = (q + 1) // 2
    return sorted({(-c1 + sd) * i2 % q, (-c1 - sd) * i2 % q})

def pm_sqrt(z, q):
    t = sqrt_mod(z, q)
    return [] if t is None else sorted({t, (q - t) % q})

def ideals(side, q):
    out = []
    for s in quad_roots(-H1, -H0, q):
        ws = (W0 + W1 * s) % q
        if side == 0:
            for X in quad_roots(C1, C0, q):
                out += [(s, t) for t in pm_sqrt(ws * X, q)]
        else:
            out += [(s, t) for t in pm_sqrt(U * ws, q)]
    return out

def norms(a0, a1, b0, b1):
    A0 = a0 * a0 + H0 * a1 * a1
    A1 = 2 * a0 * a1 + H1 * a1 * a1
    B0 = b0 * b0 + H0 * b1 * b1
    B1 = 2 * b0 * b1 + H1 * b1 * b1
    V0 = W0 * B0 + H0 * W1 * B1
    V1 = W0 * B1 + W1 * B0 + H1 * W1 * B1
    g0 = A0 - U * V0
    g1 = A1 - U * V1
    f0 = (A0 * A0 + H0 * A1 * A1) + C1 * (A0 * V0 + H0 * A1 * V1) + C0 * (V0 * V0 + H0 * V1 * V1)
    f1 = ((2 * A0 * A1 + H1 * A1 * A1) + C1 * (A0 * V1 + A1 * V0 + H1 * A1 * V1)
          + C0 * (2 * V0 * V1 + H1 * V1 * V1))
    return abs(f0 * f0 + H1 * f0 * f1 - H0 * f1 * f1), abs(g0 * g0 + H1 * g0 * g1 - H0 * g1 * g1)

def lll(basis):
    M = [list(v) for v in basis]
    n = len(M)

    def gs():
        D = [1] * (n + 1)
        lam = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1):
                u = sum(x * y for x, y in zip(M[i], M[j]))
                for k in range(j):
                    u = (D[k + 1] * u - lam[i][k] * lam[j][k]) // D[k]
                if j < i:
                    lam[i][j] = u
                else:
                    D[i + 1] = u
        return D, lam
    k = 1
    while k < n:
        D, lam = gs()
        for j in range(k - 1, -1, -1):
            d = D[j + 1]
            qq = (2 * lam[k][j] + d) // (2 * d)
            if qq:
                M[k] = [x - qq * y for x, y in zip(M[k], M[j])]
                D, lam = gs()
        if 4 * D[k + 1] * D[k - 1] >= 3 * D[k] * D[k] - 4 * lam[k][k - 1] ** 2:
            k += 1
        else:
            M[k], M[k - 1] = M[k - 1], M[k]
            k = max(k - 1, 1)
    return M

def rho_table(h=0.002):
    m = int(round(1 / h))
    rho = [1.0] * (20 * m + 2)
    integ = 1.0
    for i in range(m + 1, len(rho)):
        rho[i] = (integ + 0.5 * h * (rho[i - 1] - rho[i - m] - rho[i - m - 1])) / (i * h - 0.5 * h)
        integ = rho[i] * i * h
    return rho

def cofactor(n, prod):
    while True:
        g = gcd(n, prod)
        if g == 1:
            return n
        n //= g

class Ctx:
    def __init__(self):
        self.bad = sorted(q for q in primes_upto(abs(2 * HD * U * C0 * PD)) if (2 * HD * U * C0 * PD) % q == 0)
        self.badprod = 1
        for q in self.bad:
            self.badprod *= q
        self.primes = primes_upto(8 * B)
        self.prod = 1
        chunk = 1
        for q in self.primes:
            if q > B:
                break
            chunk *= q
            if chunk.bit_length() > 4000:
                self.prod *= chunk
                chunk = 1
        self.prod *= chunk
        self.small = 1
        for q in self.primes:
            if q > 2048:
                break
            self.small *= q
        self.qs = [q for q in self.primes if q > B]
        self.rho = rho_table()
        self.rscale = 1.0 / (0.002 * log(B))
        self.table = None

    def cof(self, n):
        c = cofactor(n, self.small)
        if c == 1:
            return 1
        if c < 4194304:
            return 1 if c <= B else c
        if pow(2, c - 1, c) == 1 and is_prime(c):
            return c
        return cofactor(c, self.prod)

    def draw(self, rng, side, parent=False):
        while True:
            if parent:
                q = rng.randrange(8 * B + 1, 64 * B, 2)
                if not is_prime(q):
                    continue
            else:
                q = self.qs[rng.randrange(len(self.qs))]
            ids = ideals(side, q)
            if ids and rng.randrange(8) < len(ids):
                s, t = ids[rng.randrange(len(ids))]
                return q, s, t

    def sweep(self, side, q, s, t, reached=None, keep=None, radius=None):
        (p0, p1, p2, p3), (q0, q1, q2, q3), (r0, r1, r2, r3), (s0, s1, s2, s3) = lll(
            [(q, 0, 0, 0), (-s, 1, 0, 0), (t, 0, 1, 0), (t * s % q, 0, 0, 1)])
        cof, badprod, rho, rs = self.cof, self.badprod, self.rho, self.rscale
        tf, tg = self.table if self.table else (None, None)
        Rr = radius or R
        rg = range(-Rr, Rr + 1)
        nrel = nadm = 0
        pd = pt = 0.0
        B8 = 8 * B
        for c0 in range(0, Rr + 1):
            for c1 in rg:
                if c0 == 0 and c1 < 0:
                    continue
                z01 = c0 == 0 and c1 == 0
                for c2 in rg:
                    if z01 and c2 < 0:
                        continue
                    lo = 1 if (z01 and c2 == 0) else -Rr
                    a0 = c0 * p0 + c1 * q0 + c2 * r0 + (lo - 1) * s0
                    a1 = c0 * p1 + c1 * q1 + c2 * r1 + (lo - 1) * s1
                    b0 = c0 * p2 + c1 * q2 + c2 * r2 + (lo - 1) * s2
                    b1 = c0 * p3 + c1 * q3 + c2 * r3 + (lo - 1) * s3
                    for _ in range(lo, Rr + 1):
                        a0 += s0
                        a1 += s1
                        b0 += s2
                        b1 += s3
                        na = a0 * a0 + H1 * a0 * a1 - H0 * a1 * a1
                        nb = b0 * b0 + H1 * b0 * b1 - H0 * b1 * b1
                        if na == 0 or nb == 0 or gcd(na, nb) != 1:
                            continue
                        Nf, Ng = norms(a0, a1, b0, b1)
                        if gcd(Nf * Ng, badprod) != 1:
                            continue
                        nadm += 1
                        if side == 0:
                            Nf //= q
                        else:
                            Ng //= q
                        pd += rho[int(log(Nf) * rs + 0.5)] * rho[int(log(Ng) * rs + 0.5)]
                        if tf:
                            pt += tf[Nf.bit_length()] * tg[Ng.bit_length()]
                        mf, mg = Nf, Ng
                        if side == 0:
                            while mf % q == 0:
                                mf //= q
                        else:
                            while mg % q == 0:
                                mg //= q
                        cg = cof(mg)
                        if cg > B8:
                            if keep is not None and len(keep) < 4000:
                                keep.append((0, mf, mg))
                            continue
                        cf = cof(mf)
                        ok = cf == 1 and cg == 1
                        nrel += ok
                        if keep is not None and len(keep) < 4000:
                            keep.append((int(ok), mf, mg))
                        if reached is not None and not ok and cf * cg <= B8 and (cf == 1 or cg == 1):
                            l = cf * cg
                            sd = 0 if cf > 1 else 1
                            for s_ in quad_roots(-H1, -H0, l):
                                bs = (b0 + b1 * s_) % l
                                if bs == 0:
                                    continue
                                t_ = (a0 + a1 * s_) * pow(bs, -1, l) % l
                                if (s_, t_) in ideals(sd, l):
                                    reached.append((sd, l, s_, t_))
        return nrel, nadm, pd, pt

    def build_table(self, tries, rng):
        tal = [{}, {}]
        Ms = (4, 8, 16, 32, 64, 128, 256, 512)
        for M in Ms:
            for _ in range(tries // len(Ms)):
                a0, a1, b0, b1 = [rng.randint(-M, M) for _ in range(4)]
                na = a0 * a0 + H1 * a0 * a1 - H0 * a1 * a1
                nb = b0 * b0 + H1 * b0 * b1 - H0 * b1 * b1
                if na == 0 or nb == 0 or gcd(na, nb) != 1:
                    continue
                N2 = norms(a0, a1, b0, b1)
                if gcd(N2[0] * N2[1], self.badprod) != 1:
                    continue
                for sd in (0, 1):
                    e = tal[sd].setdefault(N2[sd].bit_length(), [0, 0])
                    e[1] += 1
                    e[0] += self.cof(N2[sd]) == 1
        out = []
        for sd in (0, 1):
            good = sorted(k for k, e in tal[sd].items() if e[1] >= 200)
            near = lambda k: tal[sd][min(good, key=lambda g_: abs(g_ - k))]
            out.append([round(near(k)[0] / near(k)[1], 5) for k in range(140)])
        return out

def full_check(n):
    big, d = 1, 2
    while d * d <= n:
        while n % d == 0:
            n //= d
            big = d
        d += 1 if d == 2 else 2
    return max(big, n)

def stats(rows):
    n = len(rows)
    if n == 0:
        return {"balls": 0}
    re = sorted(r[1] for r in rows)
    mean = sum(re) / n
    var = sum((x - mean) ** 2 for x in re) / n
    ratio = sorted(r[1] / r[3] for r in rows if r[3] > 0)
    qt = lambda v, p: v[min(len(v) - 1, int(p * len(v)))]
    width = max(1, -(-(re[-1] + 1) // 20))
    hist = [0] * 20
    for x in re:
        hist[x // width] += 1
    return {
        "balls": n, "empty": sum(1 for x in re if x == 0), "min": re[0], "q01": qt(re, 0.01),
        "q05": qt(re, 0.05), "median": qt(re, 0.5), "mean": round(mean, 3), "max": re[-1],
        "var_over_mean": round(var / mean, 3) if mean else None,
        "pred_dickman_mean": round(sum(r[2] for r in rows) / n, 3),
        "pred_table_mean": round(sum(r[3] for r in rows) / n, 3),
        "sum_real_over_dickman": round(sum(re) / sum(r[2] for r in rows), 4),
        "sum_real_over_table": round(sum(re) / sum(r[3] for r in rows), 4),
        "ratio_table_min": round(ratio[0], 4), "ratio_table_q01": round(qt(ratio, 0.01), 4),
        "ratio_table_q05": round(qt(ratio, 0.05), 4), 
        "P_below_mean_over_4": round(sum(1 for x in re if x < mean / 4) / n, 5),
        "P_below_mean_over_8": round(sum(1 for x in re if x < mean / 8) / n, 5),
        "hist_bin_width": width, "hist": hist,
    }

def load_summary():
    try:
        with open(os.path.join(HERE, "ballcount.json")) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None

def run(seed, balls=3, radius=None):
    t0 = time.time()
    C = Ctx()
    summ = load_summary()
    rng = random.Random(seed)
    C.table = summ["table"] if summ else C.build_table(4000, rng)
    out = []
    for side in (0, 1):
        for _ in range(balls):
            q, s, t = C.draw(rng, side)
            n, nadm, pd, pt = C.sweep(side, q, s, t, radius=radius)
            out.append({"side": "fg"[side], "q": q, "s": s, "t": t, "realized": n, "admissible": nadm,
                        "predicted_dickman": round(pd, 3), "predicted_table": round(pt, 3)})
    return {"seed": seed, "params": params(C), "balls": out, "seconds": round(time.time() - t0, 3),
            "radius": radius or R,
            "campaign": {k: v for k, v in summ.items() if k != "table"} if summ else None}

def params(C):
    return {"h": [H1, H0], "w": [W0, W1], "P": [C1, C0], "u": U, "p": U * U + C1 * U + C0, "B": B,
            "q_range": [B, 8 * B], "R": R, "bad_primes": C.bad}

def part(k, seconds):
    C = Ctx()
    with open(os.path.join(HERE, "ballcount_table.json")) as fh:
        C.table = json.load(fh)
    rng = random.Random(1000 + k)
    t0 = time.time()
    rows = {"rand_f": [], "rand_g": [], "reach_f": [], "reach_g": []}
    keep, parents, pending, seen = [], 0, [], set()
    while time.time() - t0 < seconds:
        for side in (0, 1):
            q, s, t = C.draw(rng, side)
            n, _, pd, pt = C.sweep(side, q, s, t, keep=keep if len(rows["rand_f"]) % 50 == 0 else None)
            rows["rand_" + "fg"[side]].append([q, n, round(pd, 3), round(pt, 3)])
        if len(rows["rand_f"]) % 2:
            continue
        if not pending:
            q, s, t = C.draw(rng, parents % 2, parent=True)
            parents += 1
            C.sweep((parents - 1) % 2, q, s, t, reached=pending, radius=4)
            rng.shuffle(pending)
            del pending[3:]
        if pending:
            sd, l, s, t = pending.pop()
            if (sd, l, s, t) not in seen:
                seen.add((sd, l, s, t))
                n, _, pd, pt = C.sweep(sd, l, s, t)
                rows["reach_" + "fg"[sd]].append([l, n, round(pd, 3), round(pt, 3)])
    ver = [0, 0]
    for ok, mf, mg in rng.sample(keep, min(150, len(keep))):
        ver[0] += 1
        ver[1] += (max(full_check(mf), full_check(mg)) <= B) == bool(ok)
    with open(os.path.join(HERE, "ballcount_part_%d.json" % k), "w") as fh:
        json.dump({"rows": rows, "parents": parents, "verify": ver, "seconds": time.time() - t0}, fh)

def merge():
    C = Ctx()
    rows, parents, ver, secs, k = {}, 0, [0, 0], 0.0, 0
    while os.path.exists(os.path.join(HERE, "ballcount_part_%d.json" % k)):
        with open(os.path.join(HERE, "ballcount_part_%d.json" % k)) as fh:
            d = json.load(fh)
        for key, v in d["rows"].items():
            rows.setdefault(key, []).extend(v)
        parents += d["parents"]
        ver = [ver[0] + d["verify"][0], ver[1] + d["verify"][1]]
        secs += d["seconds"]
        k += 1
    with open(os.path.join(HERE, "ballcount_table.json")) as fh:
        table = json.load(fh)
    doc = {"params": params(C), "workers": k, "cpu_seconds": round(secs), "parent_balls": parents,
           "verified_cells": ver[0], "verified_agree": ver[1],
           "populations": {key: stats(v) for key, v in sorted(rows.items())}, "table": table}
    with open(os.path.join(HERE, "ballcount.json"), "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    return doc

if __name__ == "__main__":
    a = sys.argv
    if len(a) > 1 and a[1] == "table":
        tb = Ctx().build_table(int(a[2]), random.Random(int(a[3])))
        with open(os.path.join(HERE, "ballcount_table.json"), "w") as fh:
            json.dump(tb, fh)
    elif len(a) > 1 and a[1] == "part":
        part(int(a[2]), float(a[3]))
    elif len(a) > 1 and a[1] == "merge":
        print(json.dumps({k: v for k, v in merge().items() if k != "table"}, indent=1))
    else:
        print(json.dumps(run(int(a[1]) if len(a) > 1 else 1, int(a[2]) if len(a) > 2 else 3), indent=1))
