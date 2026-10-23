"""Toy-scale runs of the proof's exact walk: parallel rho with distinguished points on classes under negation and the
order-3 endomorphism.

Everything the walk reads is the same as in the claim: the canonical class point (smallest of the three x values as
integers, y below p/2), the same 128-bit multiply-xorshift mixing of the six 64-bit limbs of that x, the same bit
fields (table index = hash bits 0..k-1, distinguished test = hash bits 16..16+d-1), the same fixed-period cycle rule
(compare before save, period ALPHA, counter reset only at a start or an escape) and the same cycle-determined escape
(walk the cycle once, restore the state at the member with the smallest canonical x, double it, reduce the scalar
accumulators). Group elements are stored as discrete logarithms of a fixed generator, an exact isomorphism: the walk
reads only the x and y coordinates of the point it holds, exactly as the production walk does. Toy groups only,
supporting evidence at toy scale, not a production complexity proof.

Part 1 measures, over several seeds, the number of class visits until the first collision between two different
walkers, for this walk and for a control that jumps to a uniformly random class at every step, and reports the mean
ratio to sqrt(pi N'/2) and to the control with 95% confidence intervals. Part 2 runs the full recovery of x through
distinguished points. Part 3 checks the two hash fields for uniformity and joint independence over every class of one
toy group.
"""
import json
import math
import random
import sys

MASK = (1 << 64) - 1
M1, M2, M3 = 0x9E3779B97F4A7C15, 0xBF58476D1CE4E5B9, 0x94D049BB133111EB
A0, B0 = 0x243F6A8885A308D3, 0x13198A2E03707344
ALPHA = 64  # checkpoint period, the production value
GROUPS = ((2251, 7, 2341), (30091, 2, 29881), (130003, 3, 130261))  # (p, b, n): y^2=x^3+b, prime order n = 1 mod 3
CONFIGS = ((0, 16, 4), (1, 16, 8), (1, 64, 8), (2, 64, 8))  # (group, table size, walkers)
SEEDS = (20261023, 20261024, 20261025, 20261026, 20261027)
DP_BITS = 4  # toy distinguished width; production uses 98
UNIFORM_INDEX_BITS = 4
UNIFORM_DP_BITS = 4


def hash128(x):
    """Mixing of the canonical x: six 64-bit limbs into two lanes, one cross mix, two finalisers. 128-bit output."""
    a, b = A0, B0
    for i in range(6):
        limb = (x >> (64 * i)) & MASK
        a = ((a ^ limb) * M1) & MASK
        a ^= a >> 32
        b = ((b + limb) * M2) & MASK
        b ^= b >> 29
    b ^= a
    out = []
    for z in (a, b):
        z ^= z >> 30
        z = (z * M2) & MASK
        z ^= z >> 27
        z = (z * M3) & MASK
        z ^= z >> 31
        out.append(z)
    return out[0], out[1]


def index_of(ab, nT):
    """Table index: low bits of the first hash lane (production: 16 bits, table 2^16)."""
    return ab[0] & (nT - 1)


def is_dp(ab, bits):
    """Distinguished test: bits 16..16+bits-1 of the 128-bit hash vanish (production: 98 bits, probability 2^-98)."""
    h0, h1 = ab
    if bits <= 48:
        return (h0 >> 16) & ((1 << bits) - 1) == 0
    return (h0 >> 16) == 0 and (h1 & ((1 << (bits - 48)) - 1)) == 0


def multiples(p, b, n):
    """All n multiples of a point of order n, in Jacobian coordinates, batch-normalised with one inversion."""
    for px in range(1, p):
        rhs = (px * px * px + b) % p
        if pow(rhs, (p - 1) // 2, p) != 1:
            continue
        py = next(y for y in range(p) if y * y % p == rhs)
        break
    X, Y, Z = [0] * n, [0] * n, [0] * n
    X[1], Y[1], Z[1] = px, py, 1
    yy = py * py % p
    s = 4 * px * yy % p
    m = 3 * px * px % p
    x3 = (m * m - 2 * s) % p
    y3 = (m * (s - x3) - 8 * yy * yy) % p
    z3 = 2 * py % p
    X[2], Y[2], Z[2] = x3, y3, z3
    x1, y1, z1 = x3, y3, z3
    for k in range(3, n):  # mixed addition of (px, py) to (x1, y1, z1)
        z1z1 = z1 * z1 % p
        u2 = px * z1z1 % p
        s2 = py * z1 % p * z1z1 % p
        h = (u2 - x1) % p
        r = 2 * (s2 - y1) % p
        hh = h * h % p
        i = 4 * hh % p
        j = h * i % p
        v = x1 * i % p
        x3 = (r * r - j - 2 * v) % p
        y3 = (r * (v - x3) - 2 * y1 * j) % p
        z3 = ((z1 + h) ** 2 - z1z1 - hh) % p
        x1, y1, z1 = x3, y3, z3
        X[k], Y[k], Z[k] = x3, y3, z3
    zs = Z[2:]
    pre, acc = [1] * len(zs), 1
    for idx, z in enumerate(zs):
        pre[idx] = acc
        acc = acc * z % p
    inv = pow(acc, -1, p)
    for idx in range(len(zs) - 1, -1, -1):
        zinv = inv * pre[idx] % p
        inv = inv * zs[idx] % p
        k = idx + 2
        X[k] = X[k] * zinv % p * zinv % p
        Y[k] = Y[k] * zinv % p * zinv % p * zinv % p
    return X, Y, px, py


class Group:
    def __init__(self, p, b, n):
        self.p, self.b, self.n = p, b, n
        self.beta = next(pow(g, (p - 1) // 3, p) for g in range(2, p) if pow(g, (p - 1) // 3, p) != 1)
        self.X, self.Y, self.px, self.py = multiples(p, b, n)
        self.lam = next(l for l in range(2, n)
                        if (l * l + l + 1) % n == 0 and self.X[l] == self.beta * self.px % p)
        self.mu = [(sg * pow(self.lam, k, n)) % n for k in range(3) for sg in (1, -1)]
        self.inv = [2 * ((-(i // 2)) % 3) + i % 2 for i in range(6)]
        self.classes = (n - 1) // 6
        self.cinfo = [None] * n  # s -> (x_C, unit index with unit*R = C)
        self.hcache = {}

    def canon(self, s):
        c = self.cinfo[s]
        if c is None:
            p, x = self.p, self.X[s]
            x1 = self.beta * x % p
            x2 = (-(x + x1)) % p
            if x <= x1 and x <= x2:
                k, xc = 0, x
            elif x1 <= x2:
                k, xc = 1, x1
            else:
                k, xc = 2, x2
            c = self.cinfo[s] = (xc, 2 * k + (1 if 2 * self.Y[s] > p else 0))
        return c

    def hashed(self, xc):
        h = self.hcache.get(xc)
        if h is None:
            h = self.hcache[xc] = hash128(xc)
        return h

    def sanity(self):
        n = self.n
        assert self.lam * self.lam % n == (-self.lam - 1) % n
        assert self.beta * self.beta % self.p == (-self.beta - 1) % self.p
        assert len({self.canon(s)[0] for s in range(1, n)}) == self.classes
        assert all((y * y - (x * x * x + self.b)) % self.p == 0 for x, y in zip(self.X[1:], self.Y[1:]))
        assert (self.X[n - 1], self.Y[n - 1]) == (self.px, self.p - self.py)


def step(g, s, a, c, ts, tab, nT):
    """One walk step: index from the hash of the canonical x, add the inverse image of that table point, update (a,c).
    Returns the new state and the canonical x of the class it started from."""
    xc, i = g.canon(s)
    h = index_of(g.hashed(xc), nT)
    j = g.mu[g.inv[i]]
    if tab is None:
        return (s + j * ts[h]) % g.n, 0, 0, xc
    return (s + j * ts[h]) % g.n, (a + j * tab[h][0]) % g.n, (c + j * tab[h][1]) % g.n, xc


def escape(g, s, a, c, ts, tab, nT):
    """The walker is on a cycle. Walk it once, restore the state at the member with the smallest canonical x, double it
    and reduce the scalars. The result is a function of the cycle alone, so every walker on the cycle leaves to the
    same class, whatever its entry point or checkpoint phase."""
    xc0 = g.canon(s)[0]
    best_xc, best_state, cur, length = None, None, (s, a, c), 0
    while True:
        ns, na, nc, _ = step(g, cur[0], cur[1], cur[2], ts, tab, nT)
        nxt = (ns, na, nc)
        xc = 0 if ns == 0 else g.canon(ns)[0]
        if best_xc is None or xc < best_xc:
            best_xc, best_state = xc, nxt
        cur, length = nxt, length + 1
        if ns == 0 or xc == xc0:
            bs, ba, bc = best_state
            return (2 * bs % g.n, 2 * ba % g.n, 2 * bc % g.n), length


def collide_walk(g, nT, m, rng, alpha, stats, cap):
    """Class visits until two different walkers report the same class, with the fixed-period cycle rule."""
    n = g.n
    ts = [rng.randrange(n) for _ in range(nT)]
    ws = [[rng.randrange(1, n), None, 0] for _ in range(m)]  # s, saved class, counter
    seen, visits = {}, 0
    while visits < cap:
        for wid, w in enumerate(ws):
            s = w[0]
            xc, i = g.canon(s)
            visits += 1
            if seen.setdefault(xc, wid) != wid:
                return visits
            if w[1] == xc:  # back at the saved class: this walker is on a cycle
                (ns, _, _), length = escape(g, s, 0, 0, ts, None, nT)
                stats["events"] += 1
                stats["max_len"] = max(stats["max_len"], length)
                stats["escape_steps"] += 2 * length
                ws[wid] = [ns, None, 0]
                continue
            if w[2] % alpha == 0:
                w[1] = xc
            w[2] += 1
            nxt = step(g, s, 0, 0, ts, None, nT)[0]
            w[0] = nxt if nxt else rng.randrange(1, n)
    raise RuntimeError("no collision within the visit cap")


def collide_uniform(classes, m, rng):
    """The same collision statistic for walkers that jump to a uniformly random class at every step."""
    cur = [rng.randrange(classes) for _ in range(m)]
    seen, visits = {}, 0
    while True:
        for wid in range(m):
            visits += 1
            if seen.setdefault(cur[wid], wid) != wid:
                return visits
            cur[wid] = rng.randrange(classes)


def mean_var(xs):
    mu = sum(xs) / len(xs)
    return mu, (sum((v - mu) ** 2 for v in xs) / (len(xs) - 1) if len(xs) > 1 else 0.0)


def ratio_ci(m1, v1, n1, m2, v2, n2):
    """Delta-method 95% interval for m1/m2 with independent samples."""
    r = m1 / m2
    se = math.sqrt(v1 / n1 / m2 ** 2 + m1 ** 2 * v2 / n2 / m2 ** 4)
    return [round(r - 1.96 * se, 4), round(r + 1.96 * se, 4)]


def run_config(g, nT, m, trials, seeds, alpha):
    stats = {"events": 0, "max_len": 0, "escape_steps": 0}
    walk, ctrl, per_seed = [], [], []
    per = max(1, trials // len(seeds))
    cap = 200 * g.classes
    for seed in seeds:
        rng = random.Random(seed)
        w = [collide_walk(g, nT, m, rng, alpha, stats, cap) for _ in range(per)]
        u = [collide_uniform(g.classes, m, rng) for _ in range(per)]
        walk += w
        ctrl += u
        per_seed.append(round(sum(w) / len(w) / (sum(u) / len(u)), 4))
    m1, v1 = mean_var(walk)
    m2, v2 = mean_var(ctrl)
    ref = math.sqrt(math.pi * g.classes / 2)
    return {"p": g.p, "classes": g.classes, "table": nT, "walkers": m, "trials": len(walk),
            "mean_visits_walk": round(m1, 2), "mean_visits_uniform_control": round(m2, 2),
            "sqrt_pi_n_over_2": round(ref, 2),
            "ratio_to_sqrt_pi_n_over_2": round(m1 / ref, 4),
            "ratio_to_sqrt_pi_n_over_2_ci95": ratio_ci(m1, v1, len(walk), ref, 0.0, 1),
            "ratio_to_uniform_control": round(m1 / m2, 4),
            "ratio_to_uniform_control_ci95": ratio_ci(m1, v1, len(walk), m2, v2, len(ctrl)),
            "per_seed_ratio_to_control": per_seed,
            "charged_penalty_sqrt_1_plus_1_over_table": round(math.sqrt(1 + 1.0 / nT), 4),
            "cycle_events": stats["events"], "longest_cycle": stats["max_len"],
            "escape_steps_per_visit": round(stats["escape_steps"] / sum(walk), 5)}


def solve(g, x, nT, bits, m, rng, alpha):
    """Full recovery through distinguished points: the same walk, with scalar accumulators and one final solve."""
    n = g.n
    tab = [(rng.randrange(n), rng.randrange(n)) for _ in range(nT)]
    ts = [(a + c * x) % n for a, c in tab]  # discrete logs of the table points a*P + c*T
    seen, visits = {}, 0

    def fresh():
        while True:
            a, c = rng.randrange(n), rng.randrange(n)
            if (a + c * x) % n:
                return [(a + c * x) % n, a, c, None, 0]

    ws = [fresh() for _ in range(m)]
    while visits < 200 * g.classes:
        for wi, w in enumerate(ws):
            s, a, c = w[0], w[1], w[2]
            xc, i = g.canon(s)
            ab = g.hashed(xc)
            visits += 1
            u = g.mu[i]
            if is_dp(ab, bits):
                ka, kc = u * a % n, u * c % n
                if xc in seen and seen[xc] != (ka, kc):
                    a2, c2 = seen[xc]
                    if (c2 - kc) % n:
                        return (ka - a2) * pow(c2 - kc, -1, n) % n, visits
                seen[xc] = (ka, kc)
                ws[wi] = fresh()
                continue
            if w[3] == xc:  # back at the saved class: cycle
                (ns, na, nc), _ = escape(g, s, a, c, ts, tab, nT)
                ws[wi] = [ns, na, nc, None, 0]
                continue
            if w[4] % alpha == 0:
                w[3] = xc
            w[4] += 1
            ns, na, nc, _ = step(g, s, a, c, ts, tab, nT)
            if ns == 0:
                ws[wi] = fresh()
                continue
            w[0], w[1], w[2] = ns, na, nc
    raise RuntimeError("no collision within the visit cap")


def uniformity(g, nT, bits):
    """The two hash fields over every class: index uniformity and independence of the distinguished test from it."""
    counts, dp = [0] * nT, [0] * nT
    seen = set()
    for s in range(1, g.n):
        xc = g.canon(s)[0]
        if xc in seen:  # the six representatives of a class share the canonical x; count the class once
            continue
        seen.add(xc)
        ab = g.hashed(xc)
        j = index_of(ab, nT)
        counts[j] += 1
        if is_dp(ab, bits):
            dp[j] += 1
    total = sum(counts)
    exp = total / nT
    chi2 = sum((v - exp) ** 2 / exp for v in counts)
    pdp = sum(dp) / total
    jchi = 0.0
    for j in range(nT):
        e = counts[j] * pdp
        if e > 0:
            jchi += (dp[j] - e) ** 2 / e
        e2 = counts[j] * (1 - pdp)
        if e2 > 0:
            jchi += (counts[j] - dp[j] - e2) ** 2 / e2
    return {"classes": total, "index_bits": int(round(math.log2(nT))), "dp_bits": bits,
            "index_chi2": round(chi2, 2), "index_dof": nT - 1,
            "dp_rate": round(pdp, 4), "dp_rate_expected": round(2.0 ** -bits, 4),
            "dp_rate_by_index": [round(d / c, 4) if c else None for d, c in zip(dp, counts)],
            "joint_chi2": round(jchi, 2), "joint_dof": nT - 1}


def main():
    cfg = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    trials = int(cfg.get("trials", 20))
    rng = random.Random(SEEDS[0])
    groups = [Group(*t) for t in GROUPS]
    for g in groups:
        g.sanity()
    out = [run_config(groups[gi], nT, m, trials, SEEDS, ALPHA) for gi, nT, m in CONFIGS]
    g = groups[1]
    tot, solves = 0, 6
    for _ in range(solves):
        x = rng.randrange(1, g.n)
        got, visits = solve(g, x, 64, DP_BITS, 8, rng, ALPHA)
        assert got == x, (x, got)
        tot += visits
    uni = uniformity(groups[2], 1 << UNIFORM_INDEX_BITS, UNIFORM_DP_BITS)
    ratios = [c["ratio_to_sqrt_pi_n_over_2"] for c in out]
    print(json.dumps({
        "walk": "class-6 canonical x, 128-bit hash index and distinguished bits, fixed-period cycle rule, "
                "cycle-minimum doubling escape",
        "collision_experiments": out,
        "mean_ratio_to_sqrt_pi_n_over_2_over_configs": round(sum(ratios) / len(ratios), 4),
        "uniformity": uni,
        "recovery": {"group_p": g.p, "solves": solves, "all_recovered": True,
                     "mean_visits": round(tot / solves, 1), "sqrt_class_count": round(math.sqrt(g.classes), 1)},
        "scope": "toy groups of 390 to 21710 classes with tables of 16 to 64 entries; supporting evidence only, "
                 "it does not validate production cost"}))


main()
