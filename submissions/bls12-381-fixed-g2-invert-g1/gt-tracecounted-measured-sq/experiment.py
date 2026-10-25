"""Frobenius-class rho with a cyclotomic-squaring branch read from a trace invariant: identities, ledger, toy runs.

Part A  exact identities on BLS12-381 itself (key formula, canonical conjugate,
        Granger-Scott squaring with operation counts, orders of p and 2 modulo r).
Part C  the production step rule (twisted-table multiply, squaring branch, branch and
        distinguished-point bits taken from K) executed in the real BLS12-381 GT for
        a few thousand steps, checked against explicit Frobenius arithmetic.
Part B  the complete algorithm (cycle escapes, distinguished points, collision
        solving) in a small field of the same tower shape, with first-collision
        statistics compared with the random-mapping law adjusted for the branch mix.
Part D  exact recomputation of the work ledger used in proof.md, with the squaring price read from the
        shipped measurement record (measure/: results.tsv and the raw Callgrind records).
Toy runs are evidence at their own scale only; they do not validate production costs.
"""
import json, math, os, random, sys, time
from fractions import Fraction

T0 = time.time()
cfg = json.load(sys.stdin)
Q_TOY = cfg.get("q", 67)
TRIALS = cfg.get("trials", 20)
R = cfg.get("R", 64)
SQ8 = cfg.get("sq8", 3)
DP_BITS = cfg.get("dp_bits", 3)
WALKS = cfg.get("walks", 4)
SEED = cfg.get("seed", 1)
PROD_CHECKS = cfg.get("production_checks", 2)
PROD_STEPS = cfg.get("production_steps", 3000)
BUDGET = cfg.get("seconds", 12)
PAIR_SAMPLES = cfg.get("pair_samples", 0)
rng = random.Random(SEED)

M64 = (1 << 64) - 1


def mix64(x):
    x &= M64
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
    return x ^ (x >> 31)


class Tower:
    """F_q2 = F_q[u]/(u^2+1), F_q12 = F_q2[w]/(w^6-xi), xi = 1+u (the blst/profile tower)."""

    def __init__(self, q):
        self.q = q
        self.count = None
        self.d = [self.f2pow((1, 1), m * (q - 1) // 6) for m in range(6)]
        z = self.f2pow((1, 1), (q * q - 1) // 6)
        assert z[1] == 0 and (z[0] * z[0] - z[0] + 1) % q == 0
        self.zeta = z[0]
        self.one = ((1, 0),) + ((0, 0),) * 5
        self._fk = {}

    def f2mul(self, a, b):
        q = self.q
        return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

    def f2pow(self, a, e):
        r = (1, 0)
        while e:
            if e & 1:
                r = self.f2mul(r, a)
            a = self.f2mul(a, a)
            e >>= 1
        return r

    def f2sqr(self, a):
        if self.count is not None:
            self.count["fp2_sqr"] += 1
        q = self.q
        return (((a[0] + a[1]) * (a[0] - a[1])) % q, (2 * a[0] * a[1]) % q)

    def f2add(self, a, b):
        if self.count is not None:
            self.count["fp2_addsub"] += 1
        q = self.q
        return ((a[0] + b[0]) % q, (a[1] + b[1]) % q)

    def f2sub(self, a, b):
        if self.count is not None:
            self.count["fp2_addsub"] += 1
        q = self.q
        return ((a[0] - b[0]) % q, (a[1] - b[1]) % q)

    def f2mulxi(self, a):
        if self.count is not None:
            self.count["fp2_addsub"] += 1
        q = self.q
        return ((a[0] - a[1]) % q, (a[0] + a[1]) % q)

    def mul(self, x, y):
        q = self.q
        acc = [[0, 0] for _ in range(11)]
        for i in range(6):
            a0, a1 = x[i]
            for j in range(6):
                b0, b1 = y[j]
                t = acc[i + j]
                t[0] += a0 * b0 - a1 * b1
                t[1] += a0 * b1 + a1 * b0
        out = []
        for m in range(6):
            c0, c1 = acc[m]
            if m + 6 <= 10:
                h0, h1 = acc[m + 6]
                c0 += h0 - h1
                c1 += h0 + h1
            out.append((c0 % q, c1 % q))
        return tuple(out)

    def pow(self, x, e):
        r = self.one
        while e:
            if e & 1:
                r = self.mul(r, x)
            x = self.mul(x, x)
            e >>= 1
        return r

    def frobk(self, y, k):
        k %= 12
        if k == 0:
            return y
        D = self._fk.get(k)
        if D is None:
            D = [self.f2pow((1, 1), m * (self.q ** k - 1) // 6) for m in range(6)]
            self._fk[k] = D
        q = self.q
        if k & 1:
            return tuple(self.f2mul((c[0], (-c[1]) % q), D[m]) for m, c in enumerate(y))
        return tuple(self.f2mul(c, D[m]) for m, c in enumerate(y))

    def keys(self, y):
        """Re(w^1 coefficient) of y^(q^k), k=0..11, from 4 multiplications in F_q."""
        q, z = self.q, self.zeta
        a = y[1][0]
        b = (y[1][0] * self.d[1][0] + y[1][1] * self.d[1][1]) % q
        za, zb = z * a % q, z * b % q
        A = [a, za, (za - a) % q, (-a) % q, (-za) % q, (a - za) % q]
        B = [b, zb, (zb - b) % q, (-b) % q, (-zb) % q, (b - zb) % q]
        out = [0] * 12
        for j in range(6):
            out[2 * j] = A[j]
            out[2 * j + 1] = B[j]
        return out

    def keys_direct(self, y):
        out = []
        for k in range(12):
            out.append(self.frobk(y, k)[1][0])
        return out

    def fp4sqr(self, a, b):
        a2, b2 = self.f2sqr(a), self.f2sqr(b)
        ab2 = self.f2sub(self.f2sub(self.f2sqr(self.f2add(a, b)), a2), b2)
        return self.f2add(a2, self.f2mulxi(b2)), ab2

    def cyc_sqr(self, y):
        """Granger-Scott squaring, valid on the cyclotomic subgroup (so on GT).
        With A=(c0,c3), B=(c1,c4), C=(c2,c5) in F_q4 = F_q2[v]/(v^2-xi), v=w^3:
          c0'=3A0-2c0, c3'=3A1+2c3, c1'=3(xi C1)+2c1, c4'=3C0-2c4, c2'=3B0-2c2, c5'=3B1+2c5
        where (A0,A1)=A^2, (B0,B1)=B^2, (C0,C1)=C^2."""
        c0, c1, c2, c3, c4, c5 = y
        A0, A1 = self.fp4sqr(c0, c3)
        B0, B1 = self.fp4sqr(c1, c4)
        C0, C1 = self.fp4sqr(c2, c5)
        xC1 = self.f2mulxi(C1)

        def minus(t, u):  # 3t - 2u = t + 2(t-u)
            d = self.f2sub(t, u)
            return self.f2add(t, self.f2add(d, d))

        def plus(t, u):  # 3t + 2u = t + 2(t+u)
            s = self.f2add(t, u)
            return self.f2add(t, self.f2add(s, s))

        return (minus(A0, c0), plus(xC1, c1), minus(B0, c2), plus(A1, c3), minus(C0, c4), plus(B1, c5))

    def rand_gt(self, rng, rp):
        hc = (self.q ** 12 - 1) // rp
        while True:
            x = tuple((rng.randrange(self.q), rng.randrange(self.q)) for _ in range(6))
            y = self.pow(x, hc)
            if y != self.one:
                return y


# ---------------------------------------------------------------- production constants
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
MR384 = pow(2, 384, p)      # Montgomery form of the pinned blst field: stored limbs = x*2^384 mod p
R_FACT = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1, 906349: 2,
          2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}


def mult_order(a, n, fact):
    o = n
    for pf, e in fact.items():
        for _ in range(e):
            if pow(a, o // pf, n) == 1:
                o //= pf
            else:
                break
    return o


def part_A(rng, checks):
    P = Tower(p)
    lam = p % r
    prod_fact = 1
    for k, e in R_FACT.items():
        prod_fact *= k ** e
    o2 = mult_order(2, r, {k: e for k, e in R_FACT.items()})
    ol = mult_order(lam, r, R_FACT)
    t_min = o2 // math.gcd(o2, 12)
    out = {
        "r_minus_1_factorisation_verified": prod_fact == r - 1,
        "order_of_p_mod_r": ol,
        "lambda6_is_minus1": pow(lam, 6, r) == r - 1,
        "one_plus_l4_plus_l8_is_0": (1 + pow(lam, 4, r) + pow(lam, 8, r)) % r == 0,
        "log2_order_of_2_mod_r": round(math.log2(o2), 3),
        "smallest_t_with_2^t_in_<p>_log2": round(math.log2(t_min), 3),
        "zeta_primitive_6th_root": (P.zeta * P.zeta - P.zeta + 1) % p == 0,
    }
    ok = dict.fromkeys(["gt", "frob", "keys", "distinct", "canon", "inv", "sq", "trace_inv", "trace_diff"], True)
    traces = []
    counts = []
    hcof = (p ** 12 - 1) // r
    for _ in range(checks):
        x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
        y = P.pow(x, hcof)
        ok["gt"] &= P.pow(y, r) == P.one and y != P.one
        ok["frob"] &= P.frobk(y, 1) == P.pow(y, p)
        kd, kf = P.keys_direct(y), P.keys(y)
        ok["keys"] &= kd == kf
        ok["distinct"] &= len(set(kf)) == 12
        conj = [P.frobk(y, k) for k in range(12)]
        canon = set()
        for c in conj:
            kk = P.keys(c)
            kap = min(range(12), key=lambda i: kk[i])
            canon.add(P.frobk(c, kap))
        ok["canon"] &= len(canon) == 1
        ok["inv"] &= P.mul(conj[6], y) == P.one
        ok["trace_inv"] &= len({c[0][0] for c in conj}) == 1       # Re(c_0) is the same for all 12 conjugates
        traces.append(y[0][0])
        P.count = {"fp2_sqr": 0, "fp2_addsub": 0}
        s = P.cyc_sqr(y)
        counts.append(dict(P.count))
        P.count = None
        ok["sq"] &= s == P.mul(y, y)
    ok["trace_diff"] &= len(set(traces)) == len(traces)
    out.update({"random_GT_elements": checks, **{k + "_ok": v for k, v in ok.items()}})
    out["granger_scott_op_counts"] = counts[0] if counts else None
    return out


def part_C(rng, steps, rprod=16, verify=150, stored=False):
    """Production step rule in the real GT: twisted-table multiply / squaring, bits from T = Re(c_0).

    stored=False reads T and the twelve keys as ordinary integers in [0,p) (what a
    canonical-residue implementation would see); stored=True reads the six limbs exactly as the
    pinned blst kernels store them, i.e. in Montgomery form x~ = x*2^384 mod p, which is the rule
    proof.md section 3 specifies and the one an implementation of this design runs at zero cost.
    Both modes run the same step and the same class checks; they differ only in which of the
    twelve keys is the minimum and in the bits taken from T.
    """
    P = Tower(p)
    hcof = (p ** 12 - 1) // r
    lam = [pow(p, k, r) for k in range(12)]
    g = P.pow(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6)), hcof)
    x0 = rng.randrange(1, 1 << 40)
    z0 = P.pow(g, x0)
    ah = [rng.randrange(1, 1 << 40) for _ in range(rprod)]
    bh = [rng.randrange(1, 1 << 40) for _ in range(rprod)]
    Mh = [P.mul(P.pow(g, ah[i]), P.pow(z0, bh[i])) for i in range(rprod)]
    Tab = [[P.frobk(Mh[i], k) for k in range(12)] for i in range(rprod)]
    # exponents are tracked exactly: the start element is g^a0 z0^b0
    a0, b0 = rng.randrange(1, r), rng.randrange(1, r)
    y = P.mul(P.pow(g, a0), P.pow(z0, b0))
    al, be = a0, b0
    n_sq = n_mul = bad_class = 0
    bins_sq = [0] * 8
    bins_h = [0] * 16
    two_cycles = 0
    idx_differs = 0
    hist = []
    for t in range(steps):
        # class invariant T = Re(c_0): no multiplication needed. In `stored` mode the limbs are read as
        # the kernel stores them (Montgomery), which costs nothing; in integer mode they are the residue.
        u = (y[0][0] * MR384) % p if stored else y[0][0]
        if len(hist) >= 2 and hist[-2][0] == u:          # fruitless 2-cycle: square the smaller-T class
            two_cycles += 1
            _, yy, aa, bb = min(hist[-1], (u, y, al, be), key=lambda s: s[0])
            y, al, be = P.cyc_sqr(yy), 2 * aa % r, 2 * bb % r
            hist = []
            continue
        hist.append((u, y, al, be))
        if len(hist) > 2:
            hist.pop(0)
        br = u & 7
        bins_sq[br] += 1
        if br < SQ8:
            y2 = P.cyc_sqr(y)                            # squaring step: no keys, no twist
            if t < verify:
                kap = min(range(12), key=lambda i: P.keys(y)[i])
                bad_class += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), P.frobk(y, kap))
            y, al, be = y2, 2 * al % r, 2 * be % r
            n_sq += 1
        else:
            kf = P.keys(y)                               # multiply step: 4 Fp multiplications give kappa
            kf_int = kf
            if stored:
                kf = [(v * MR384) % p for v in kf]       # the same keys as stored (Montgomery), no work
            kap = kf.index(min(kf))
            idx_differs += (kap != kf_int.index(min(kf_int)))
            h = (u >> 3) % rprod
            bins_h[(u >> 3) & 15] += 1
            idx = (12 - kap) % 12
            y2 = P.mul(y, Tab[h][idx])
            if t < verify:
                bad_class += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), Mh[h])
            y = y2
            al = (al + ah[h] * lam[idx]) % r
            be = (be + bh[h] * lam[idx]) % r
            n_mul += 1
    exp_ok = y == P.mul(P.pow(g, al), P.pow(z0, be))
    exp_sq = sum(bins_sq[:SQ8])
    n_br = sum(bins_h)
    n_all = sum(bins_sq)
    chi_sq = sum((c - n_all / 8) ** 2 / (n_all / 8) for c in bins_sq)
    chi_h = sum((c - n_br / 16) ** 2 / (n_br / 16) for c in bins_h) if n_br else None
    return {"steps": steps, "R_table": rprod, "squaring_steps": n_sq, "multiply_steps": n_mul,
            "key_rule": "stored Montgomery limbs (x*2^384 mod p)" if stored else "ordinary integers in [0,p)",
            "min_index_differs_from_integer_rule": idx_differs,
            "squaring_fraction": round(n_sq / max(1, n_all), 4), "target_fraction": SQ8 / 8,
            "explicit_class_check_failures_in_first_%d" % verify: int(bad_class),
            "exponent_invariant_y_equals_g^al_z^be": exp_ok,
            "chi2_low3_bits_of_Re_c0_7dof": round(chi_sq, 2),
            "chi2_next4_bits_of_Re_c0_15dof": round(chi_h, 2) if chi_h is not None else None,
            "two_cycle_events": two_cycles}


# ---------------------------------------------------------------- toy algorithm (Part B)
def toy_setup(q):
    phi = q ** 4 - q * q + 1
    rp = phi
    d = 2
    while d * d <= rp:
        while rp % d == 0 and rp > d:
            rp //= d
        d += 1
    return rp


def classinfo(T, y):
    """kappa, canonical conjugate, minimum key. Production never builds the canonical
    element; the toy does, and breaks equal-key ties by comparing full conjugates."""
    kf = T.keys(y)
    m = min(kf)
    ks = [k for k in range(12) if kf[k] == m]
    if len(ks) == 1:
        kap = ks[0]
        can = T.frobk(y, kap)
    else:
        can, kap = min((T.frobk(y, k), k) for k in ks)
    return kap, can, m


def toy_hash(can):
    h = 0x9E3779B97F4A7C15
    for c in can:
        h = mix64(h ^ (c[0] * 0x100000001B3 + c[1]))
    return h


def run_toy(T, q, rp, rng, walks, Rm, sq8, dp_bits, cap_factor=40):
    qk = [pow(q, k, rp) for k in range(12)]
    hc = (q ** 12 - 1) // rp
    g = T.rand_gt(rng, rp)
    xs = rng.randrange(1, rp)
    z0 = T.pow(g, xs)
    ah = [rng.randrange(rp) for _ in range(Rm)]
    bh = [rng.randrange(rp) for _ in range(Rm)]
    Mh = [T.mul(T.pow(g, ah[i]), T.pow(z0, bh[i])) for i in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    At = [[ah[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    Bt = [[bh[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    mask = (1 << dp_bits) - 1

    def new_walk():
        a, b = rng.randrange(rp), rng.randrange(rp)
        return {"y": T.mul(T.pow(g, a), T.pow(z0, b)), "a": a, "b": b, "hist": [], "ck": None, "n": 0, "since": 0}

    def advance(y, a, b):
        kap, can, K = classinfo(T, y)
        u = toy_hash(can)
        if (u & 7) < sq8:
            return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp, kap, can, K, u
        h = (u >> 3) % Rm
        idx = (12 - kap) % 12
        return T.mul(y, Tab[h][idx]), (a + At[h][idx]) % rp, (b + Bt[h][idx]) % rp, kap, can, K, u

    W = [new_walk() for _ in range(walks)]
    visited, dps = {}, {}
    steps, Tc, found = 0, None, None
    stats = {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0, "sq": 0, "ord": 0}
    while found is None:
        for w in W:
            y, a, b = w["y"], w["a"], w["b"]
            kap, can, K = classinfo(T, y)
            u = toy_hash(can)
            ca, cb = a * qk[kap] % rp, b * qk[kap] % rp
            if Tc is None:
                old = visited.get(can)
                if old is None:
                    visited[can] = (ca, cb)
                elif old != (ca, cb):
                    Tc = steps
            if (u >> 13) & mask == 0:
                old = dps.get(can)
                if old is not None and old[1] != cb:
                    xr = (old[0] - ca) * pow(cb - old[1], -1, rp) % rp
                    if T.pow(g, xr) != z0:
                        raise RuntimeError('recovered logarithm failed verification')
                    found = xr
                    break
                if old is not None:
                    stats["useless"] += 1
                dps[can] = (ca, cb)
                w["since"] = 0
            hist = w["hist"]
            key = (K, can)
            if len(hist) >= 2 and hist[-2][0][1] == can:
                stats["two_cycles"] += 1
                s1, s2 = hist[-1], (key, y, a, b)
                _, yy, aa, bb = min(s1, s2, key=lambda s: s[0])
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["ck"] is not None and w["ck"][0] == can and w["ck"][1] < w["n"]:
                stats["long_cycles"] += 1
                best = (key, y, a, b)
                yy, aa, bb = y, a, b
                for _ in range(w["n"] - w["ck"][1]):
                    yy, aa, bb, kp, cn, Kn, un = advance(yy, aa, bb)
                    steps += 1
                    k2, c2, K2 = classinfo(T, yy)
                    if (K2, c2) < best[0]:
                        best = ((K2, c2), yy, aa, bb)
                _, yy, aa, bb = best
                w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
                w["hist"], w["ck"] = [], None
                steps += 1
                continue
            if w["n"] % 64 == 0:
                w["ck"] = (can, w["n"])
            hist.append((key, y, a, b))
            if len(hist) > 2:
                hist.pop(0)
            stats["ord"] += 1
            stats["sq"] += (u & 7) < sq8
            w["y"], w["a"], w["b"], _, _, _, _ = advance(y, a, b)
            w["n"] += 1
            w["since"] += 1
            steps += 1
            if w["since"] > cap_factor << dp_bits:
                stats["restarts"] += 1
                w.update(new_walk())
    return Tc, steps, stats


def part_B(rng, q, trials, Rm, sq8, dp_bits, walks, budget):
    T = Tower(q)
    rp = toy_setup(q)
    if not (pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))):
        raise RuntimeError('toy field does not have embedding degree 12')
    N = (rp - 1) // 12
    s = sq8 / 8
    sumq2 = s * s + (1 - s) ** 2 / Rm
    Np = N / (1 - sumq2)
    Tcs, tot, ev, dev = [], [], {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0}, []
    t_start = time.time()
    while len(Tcs) < trials and time.time() - t_start < budget:
        Tc, steps, st = run_toy(T, q, rp, rng, walks, Rm, sq8, dp_bits)
        Tcs.append(Tc)
        tot.append(steps)
        for k in ev:
            ev[k] += st[k]
        sf = sq8 / 8
        # the toy's charge bookkeeping uses this package's price vector (section 5); the ratio
        # charge/nominal is invariant to which vector is used, both sides take the same one
        MUL_STEP, SQ_STEP = 1.11945, 0.58195
        nominal = (1 - sf) * MUL_STEP + sf * SQ_STEP
        charge = ((st["ord"] - st["sq"]) * MUL_STEP + st["sq"] * SQ_STEP) / max(1, st["ord"])
        dev.append((charge / nominal - 1, st["ord"], nominal))
    n = len(Tcs)
    sN = math.sqrt(N)
    pred = math.sqrt(math.pi * Np / 2)
    s40 = math.sqrt(2 * Np * math.log(1 / 0.60))
    s42 = math.sqrt(2 * Np * math.log(1 / 0.58))
    s43 = math.sqrt(2 * Np * math.log(1 / 0.57))
    allsteps = max(1, sum(tot))
    return {"q": q, "r_prime": rp, "classes_N": N, "R": Rm, "sq8": sq8, "squaring_prob": s,
            "sum_q2": round(sumq2, 6), "effective_N_over_N": round(1 / (1 - sumq2), 4),
            "trials_completed": n, "logs_recovered_and_verified": n,
            "mean_Tc_over_sqrtN": round(sum(Tcs) / n / sN, 4) if n else None,
            "model_mean_over_sqrtN": round(pred / sN, 4),
            "fraction_Tc_le_S_model_0.40": round(sum(t <= s40 for t in Tcs) / max(1, n), 3),
            "fraction_Tc_le_S_model_0.42": round(sum(t <= s42 for t in Tcs) / max(1, n), 3),
            "fraction_Tc_le_S_model_0.43": round(sum(t <= s43 for t in Tcs) / max(1, n), 3),
            "S_model_0.40_over_sqrtN": round(s40 / sN, 4), "S_model_0.42_over_sqrtN": round(s42 / sN, 4),
            "S_model_0.43_over_sqrtN": round(s43 / sN, 4),
            "mean_total_steps_incl_DP_delay_over_sqrtN": round(sum(tot) / max(1, n) / sN, 4),
            "charge_fluctuation_rms_of_(dev*sqrt(steps))_observed": round(math.sqrt(sum(d * d * k for d, k, _ in dev) / max(1, len(dev))), 4),
            "charge_fluctuation_iid_prediction": round(math.sqrt(s * (1 - s)) * 0.54 / dev[0][2], 4) if dev else None,
            "events": ev, "two_cycles_per_step": ev["two_cycles"] / allsteps,
            "Tc_over_sqrtN": [round(t / sN, 4) for t in Tcs] if cfg.get("dump") else None}


# ---------------------------------------------------------------- pair-collision law (Part E)
def part_E(rng, q, m, Rm, sq8):
    """Sample m distinct classes, apply one step F to each, and count image coincidences.
    H1 predicts N * P[F(c1) = F(c2)] = 1 - s^2 - (1-s)^2/(12R) for distinct classes c1, c2."""
    T = Tower(q)
    rp = toy_setup(q)
    N = (rp - 1) // 12
    g = T.rand_gt(rng, rp)
    Mh = [T.pow(g, rng.randrange(1, rp)) for _ in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    gp = [g]
    for _ in range(rp.bit_length()):
        gp.append(T.mul(gp[-1], gp[-1]))
    classes = {}
    while len(classes) < m:
        e = rng.randrange(1, rp)
        y, j = T.one, 0
        while e:
            if e & 1:
                y = T.mul(y, gp[j])
            e >>= 1
            j += 1
        kap, can, _ = classinfo(T, y)
        if can not in classes:
            classes[can] = (y, kap)
    images = {}
    for can, (y, kap) in classes.items():
        u = toy_hash(can)
        if (u & 7) < sq8:
            y2 = T.cyc_sqr(y)
        else:
            y2 = T.mul(y, Tab[(u >> 3) % Rm][(12 - kap) % 12])
        _, can2, _ = classinfo(T, y2)
        images[can2] = images.get(can2, 0) + 1
    pairs = sum(n * (n - 1) // 2 for n in images.values())
    tot = m * (m - 1) // 2
    s = sq8 / 8
    obs = pairs / tot * N
    return {"q": q, "classes_N": N, "sampled_classes": m, "R": Rm, "squaring_prob": s,
            "image_coincidence_pairs": pairs,
            "N_times_pair_collision_rate_observed": round(obs, 4),
            "relative_standard_error": round(1 / math.sqrt(max(1, pairs)), 4),
            "predicted_exact_1_minus_s2_minus_(1-s)^2/(12R)": round(1 - s * s - (1 - s) ** 2 / (12 * Rm), 4),
            "predicted_conservative_1_minus_s2_minus_(1-s)^2/R": round(1 - s * s - (1 - s) ** 2 / Rm, 4),
            "random_mapping": 1.0}


# ---------------------------------------------------------------- measurement record
OFFSETS = (0, 3, 101, 1009, 4001)


def locate_measure():
    """Directory holding results.tsv and callgrind/: /snapshot/files/measure under the grader's runner,
    measure/ beside this script otherwise; never the caller's working directory."""
    here = os.path.dirname(os.path.abspath(__file__))
    for base in (os.path.join(here, "files", "measure"), os.path.join(here, "measure"),
                 os.path.join(os.getcwd(), "files", "measure"), os.path.join(os.getcwd(), "measure")):
        if os.path.isfile(os.path.join(base, "results.tsv")):
            return base
    for root in (os.path.join(here, "files"), here):
        for dirpath, _, names in os.walk(root):
            if "results.tsv" in names and os.path.isdir(os.path.join(dirpath, "callgrind")):
                return dirpath
    raise SystemExit("experiment.py: measure/results.tsv not found under %s (looked in files/measure and "
                     "measure); the measurement record ships with the package." % here)


def load_measure():
    """Read results.tsv and the six raw Callgrind records; assert the three per-call marginals
    blst_fp_mul = 383 and blst_fp12_mul = 22,996 (the cost model's reference counts) and the
    blst_fp12_cyclotomic_sqr marginal of 9,803; return the squaring marginal."""
    base = locate_measure()
    rows = {}
    with open(os.path.join(base, "results.tsv")) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            routine, off, n, ir, gt = line.split()
            rows[(routine, int(off), int(n))] = (int(ir), gt)
    expect = {"fp-mul": 383, "fp12-mul": 22996, "cyc-sqr": 9803}
    marginals = {}
    for routine, want in expect.items():
        marg = {off: Fraction(rows[(routine, off, 65)][0] - rows[(routine, off, 1)][0], 64) for off in OFFSETS}
        one = {off: rows[(routine, off, 1)][0] for off in OFFSETS}
        assert len(set(marg.values())) == 1 and len(set(one.values())) == 1, (routine, marg, one)
        assert all(rows[(routine, off, 1)][1] == rows[(routine, off, 65)][1] for off in OFFSETS), routine
        m = marg[OFFSETS[0]]
        assert m.denominator == 1 and m == want, (routine, m, want)
        marginals[routine] = int(m)
    for routine in expect:
        for n in (1, 65):
            path = os.path.join(base, "callgrind", "%s-off0-n%d.callgrind" % (routine, n))
            with open(path) as fh:
                tot = [int(line.split()[1]) for line in fh if line.startswith("totals:")]
                fh.seek(0)
                summ = [int(line.split()[1]) for line in fh if line.startswith("summary:")]
            assert tot == summ == [rows[(routine, 0, n)][0]], (routine, n, tot, summ)
    assert all(rows[k][1] == "1" for k in rows if k[0] != "fp-mul"), "cyc2 must report GT membership 1"
    return {"measure_dir": base, "offsets": list(OFFSETS), "ir_per_call_marginal": marginals,
            "raw_callgrind_records_checked": 6,
            "ir_at_n1_and_n65_offset0": {k[0] + "@n" + str(k[2]): rows[k][0] for k in sorted(rows) if k[1] == 0}}


# ---------------------------------------------------------------- ledger (Part D)
def ledger():
    """Exact recomputation of the work ledger of proof.md section 5 (instruction counts / C, C = 20,000)."""
    from decimal import Decimal, getcontext, localcontext
    getcontext().prec = 60
    ln2 = Decimal(2).ln()
    lg = lambda v: float(Decimal(v).ln() / ln2)
    meas = load_measure()
    C = 20000
    FPMUL = meas["ir_per_call_marginal"]["fp-mul"]        # blst_fp_mul, measured (= the cost model's 383)
    FPADD = 45                                            # modular add/sub/neg as a called routine (analytic)
    SQ_MEASURED = meas["ir_per_call_marginal"]["cyc-sqr"]  # blst_fp12_cyclotomic_sqr, measured marginal
    FP2SQR = 2 * FPMUL + 3 * FPADD                        # analytic complex squaring (used for the comparison row)
    FP2ADD = 2 * FPADD
    GS_ANALYTIC = 9 * FP2SQR + 34 * FP2ADD + 200          # line-by-line count of the same formula (not charged)
    common_items = {"read Re(c_0), branch bits, DP test, table address": 50,
                    "(alpha, beta) update mod r (add or doubling)": 100,
                    "2-cycle and checkpoint comparisons on T, history": 60,
                    "loop control": 30,
                    "margin for loads and moves": 30}
    # The counters and the work meter of proof.md section 3.4 are multiword (192-bit meter, 2-word
    # stuck-walk counter); their per-step charge is itemised here and added to the 270, and the report store's
    # per-report charge is itemised below. The earlier revision priced neither.
    meter_items = {"work meter: 3-word (192-bit) accumulate, step charge in a register (1) and adds/adcs/adc (3), "
                   "then the exact 3-word compare with the cap held in three registers, cmp/sbcs/sbcs/b.hs (4), "
                   "executed on every step": 8,
                   "multiword counters: 2-word stuck-walk counter increment, adds/adc (2), and its "
                   "high-word threshold compare, cmp hi,#2560; b.hs (2), exact as 20/theta = 2560*2^64": 4}
    mult_items = {"4 Fp multiplications (a*zeta, a*d0, b*d1, B*zeta)": 4 * FPMUL,
                  "1 Fp addition (B) + 8 Fp sub/neg": 9 * FPADD,
                  "minimum of 12 keys (top limbs, tie test)": 150,
                  "load c_1, kappa, twist index": 32}
    common, extra = sum(common_items.values()), sum(mult_items.values())
    common_charged = common + sum(meter_items.values())
    assert (common, extra, GS_ANALYTIC, common_charged) == (270, 2119, 11369, 282)
    # ---------------------------------------------------------------------------------------------
    # Three price vectors (squaring routine, every-step work, multiply-only work) go through one code
    # path: the earliest ledger (rounded charges), the intermediate ledger (the itemised listing
    # without the meter and counter lines, squaring priced analytically) and this package (the meter
    # and counter lines charged, squaring priced by the measured marginal).
    PARENT = (11400, 300, 2200)
    RECORD = (GS_ANALYTIC, common, extra)
    OURS = (SQ_MEASURED, common_charged, extra)
    N = Fraction(r - 1, 12)
    s, Rt = Fraction(3, 8), 1024
    sumq2 = s * s + (1 - s) ** 2 / Rt
    Np = N / (1 - sumq2)
    Pm = Decimal("0.43")
    Npd = Decimal(Np.numerator) / Decimal(Np.denominator)
    S = (2 * Npd * (1 / (1 - Pm)).ln()).sqrt()
    eps = Decimal(5) / Decimal(12 * Rt) + Decimal(250) / Decimal(Rt * Rt)
    # Steps: the walk S (1 + eps) plus the detection window, 20/theta steps for each of the M walks,
    # plus the same allowance for each of the at most RESTART_CAP restarted walks (section 3.4): the
    # restart cap is charged here, so the cap, the meter and the success allowance rest on one number.
    # theta = 2^-THETA_BITS with THETA_BITS = 71 (proof.md section 3): the report store then fits one
    # 64-bit address space (section 0-ii). The earlier ledgers keep the DP rate 2^-60 they were priced at.
    THETA_BITS = 71
    E_OLD = Decimal(2 ** 20) * 20 * Decimal(2 ** 60)                      # detection window of the earlier revisions
    WALKS, RESTART_CAP = 2 ** 20, 2 ** 36     # 2^36 restarts vs 2^26.91 expected (proof.md 0-iii)
    STUCK = 20 * Decimal(2) ** THETA_BITS                                  # 20/theta, the stuck-walk threshold
    assert STUCK == 2560 * Decimal(2) ** 64      # so the counter's threshold test reads the high word only
    E = (WALKS + RESTART_CAP) * STUCK                                      # detection window + restart allowance
    # Expected restarts (proof.md 0-iii). The stuck counter is reset at every walk start, restart and report,
    # so the rule fires once per gap between consecutive reports, not once per walk: a gap of geometric length
    # with mean 1/theta exceeds 20/theta steps with probability (1-theta)^(20/theta) <= e^-20. The number of
    # gaps is the number of reports, expected theta times the steps. No tail bound is taken on these counts:
    # a run whose restarts exceed the cap is a cap failure, inside the 0.01 allowance of H1+ (proof.md 3.4).
    GAP_STEPS = S * (1 + eps)                                              # the steps of the walk
    EXP_GAPS = GAP_STEPS * Decimal(2) ** (-THETA_BITS)                     # expected reports = expected gap starts
    EXP_STUCK = EXP_GAPS * (-Decimal(20)).exp()                            # expected stuck-walk restarts
    EXP_LONG = GAP_STEPS * Decimal(2) ** (-100)                            # fruitless long cycles, rate below 2^-100 per step
    EXP_RESTARTS = EXP_STUCK + EXP_LONG
    assert RESTART_CAP == 2 ** 36
    assert Decimal("25.97") < lg(EXP_STUCK) < Decimal("25.98")             # quoted as 2^25.98 (rounded up)
    assert Decimal("25.83") < lg(EXP_LONG) < Decimal("25.84")              # quoted as 2^25.84
    assert Decimal("26.90") < lg(EXP_RESTARTS) < Decimal("26.91")          # quoted as 2^26.91
    assert lg(EXP_RESTARTS) - 36 < Decimal("-9.09") + Decimal("0.01") and lg(EXP_RESTARTS) < 36 - 9
    STUCK_FRAC = 20 * (-Decimal(20)).exp()                                 # steps wasted per step: 20/theta x e^-20 per gap of mean 1/theta
    LONG_FRAC = EXP_LONG * STUCK / S                                       # long-cycle events x one stuck allowance each, over S
    assert STUCK_FRAC < Decimal("4.2e-8") and LONG_FRAC < Decimal("3.8e-8")
    # Non-walk work (proof.md section 5): DP handling at the report cap, the shared meter's per-walk
    # updates, and setup. Each component is computed here so the bucket is a bound, not a constant.
    REPORT_CAP = Decimal(2) ** 55       # 1.1 S theta + 2^25 < 2^54.969, rounded up (proof.md 3.4)
    assert Decimal("1.1") * S * Decimal(2) ** (-THETA_BITS) + Decimal(2 ** 25) <= REPORT_CAP
    assert lg(Decimal("1.1") * S * Decimal(2) ** (-THETA_BITS) + Decimal(2 ** 25)) > 54.96
    SLOT_BYTES, SLOTS_PER_BUCKET, BUCKET_BITS = 80, 16, 52
    SLOTS = 2 ** BUCKET_BITS * SLOTS_PER_BUCKET
    report_items = {"12 keys for kappa (the multiply-step key list)": extra,
                    "two 255-bit modular multiplications by p^kappa mod r": 1600,
                    "report counter (one 64-bit word) increment and compare against the cap, branch": 5,
                    "bucket index (low 52 bits of F), bucket base, slot addresses (64-bit)": 10,
                    "16 slot probes (offset add, two 8-byte loads, two compares, branch)": 96,
                    "store the 80-byte record (ten 8-byte stores, address set-up)": 12,
                    "free-slot test and truncation branch": 4,
                    "reset of the 2-word stuck-walk counter at the report (two moves into the register pair)": 4,
                    "margin for loads, moves, call and return": 150}
    report_instr = sum(report_items.values())
    assert report_instr == 4000
    METER_UPDATE = 40        # two updates per walk (full 3-word compare on entry, hand-back on exit), 20 each
    dp_units = REPORT_CAP * report_instr / C
    start_units = Decimal(WALKS) * METER_UPDATE / C
    restart_units = Decimal(RESTART_CAP) * METER_UPDATE / C
    ONE_START = Decimal(2) ** Decimal("10.1")          # one fresh start g^alpha z^beta, < 1,021 units + loads
    restart_start_units = Decimal(RESTART_CAP) * ONE_START
    SETUP_OLD = Decimal(2) ** Decimal("63.8")           # non-walk bucket input of the two earlier price vectors in this recomputation (not the record's 2^64.3; the difference is below 2^-60 of the ledger)
    SETUP = Decimal(2) ** Decimal("52.7")               # non-walk bucket of this ledger (2^55 reports)
    assert dp_units + start_units + restart_units + restart_start_units + Decimal(2 ** 31) < SETUP
    # Truncation: a report is dropped when its bucket already holds SLOTS_PER_BUCKET records with other
    # fingerprints. The records already in the bucket number Bin(j, 2^-52) after j earlier reports (each
    # earlier report lands in this bucket with probability 2^-52), which grows with j, so the drop
    # probability of any report is at most the bucket-full tail at the cap, P(Bin(2^55, 2^-52) >= 16),
    # computed exactly below (every report is counted as a distinct record: an equal report is appended
    # like any other and occupies a slot, never replaced or shared). The mean is 8, so the
    # tail is that of a Poisson(8)-like count (the slots of one bucket are not independent coin flips).
    with localcontext() as _ctx:
        _ctx.prec = 80
        _n, _p = REPORT_CAP, Decimal(1) / Decimal(2 ** BUCKET_BITS)
        _lq = (1 - _p).ln()
        def _binom_pmf(k):
            c = Decimal(1)
            for i in range(k):
                c = c * (_n - i) / (i + 1)
            return c * _p ** k * (_lq * (_n - k)).exp()
        BUCKET_FULL = 1 - sum(_binom_pmf(k) for k in range(SLOTS_PER_BUCKET))
        _lam = _n * _p
        _pois = [(-_lam).exp() * _lam ** k / math.factorial(k) for k in range(SLOTS_PER_BUCKET)]
        BUCKET_FULL_POISSON = 1 - sum(_pois)
        # a useful collision is lost only if the first report of its class meets a full bucket (the
        # second then meets the same full bucket); charging both reports is the safe bound 2q.
        TRUNC_REL = 2 * BUCKET_FULL
        TRUNC_ABS = Pm * TRUNC_REL                 # absolute success lost, charged against the 0.43
    assert _lam == 8 and abs(BUCKET_FULL - BUCKET_FULL_POISSON) < Decimal("1e-15")
    assert Decimal("0.00823") < BUCKET_FULL < Decimal("0.00824")           # quoted as 0.00824, rounded up
    assert TRUNC_REL < Decimal("0.01648") and TRUNC_ABS < Decimal("0.00709")  # quoted rounded up
    assert Pm - TRUNC_ABS > Decimal("0.42291")                    # success after truncation, rounded down
    # The success budget (proof.md section 4): H1+ 0.01, truncation, and the small terms.
    SMALL_LOSSES = Decimal("4.2e-8") + Decimal("3.8e-8") + Decimal("2.1e-9") + Decimal("1e-6")   # stuck-walk steps, long-cycle steps, detection, 1/M
    LOSS_SUM = Decimal("0.01") + TRUNC_ABS + SMALL_LOSSES
    assert LOSS_SUM < Decimal("0.018")                                      # quoted as below 0.018
    assert Pm - Decimal("0.018") >= Decimal("0.39") + Decimal("0.022")      # 0.022 slack left for H1
    store = {"shape": "fixed-depth 16-slot buckets, direct bucket addressing from the low 52 bits of F, "
                      "record truncated when its bucket is full",
             "bucket_bits": BUCKET_BITS, "slots_per_bucket": SLOTS_PER_BUCKET, "slots": SLOTS,
             "record_bytes": SLOT_BYTES, "array_bytes_log2": lg(Decimal(SLOTS * SLOT_BYTES)),
             "load_factor_at_the_report_cap": float(REPORT_CAP / Decimal(SLOTS)),
             "worst_case_probes": SLOTS_PER_BUCKET,
             "mean_records_per_bucket_at_the_cap": float(_lam),
             "bucket_full_probability_at_the_cap_exact_binomial": float(BUCKET_FULL),
             "bucket_full_probability_at_the_cap_poisson8": float(BUCKET_FULL_POISSON),
             "worst_case_drop_fraction_per_report": float(BUCKET_FULL),
             "success_lost_by_truncation_relative": float(TRUNC_REL),
             "success_lost_by_truncation_absolute_of_0.43": float(TRUNC_ABS),
             "success_after_truncation": float(Pm - TRUNC_ABS),
             "success_loss_sum_h1plus_truncation_small_terms": float(LOSS_SUM),
             "success_after_every_loss": float(Pm - LOSS_SUM),
             "margin_over_0.39_after_every_loss": float(Pm - LOSS_SUM - Decimal("0.39"))}

    # Word widths (proof.md section 0-ii): every address, index and report counter fits one native
    # 64-bit word with the sign bit clear; the stuck-walk counter and the work meter are arithmetic
    # counters, not addresses, and are multiword with a priced carry chain (section 0-i).
    STRIDE = SLOTS_PER_BUCKET * SLOT_BYTES                       # 1,280 = 2^10 + 2^8
    assert STRIDE == (1 << 10) + (1 << 8)
    max_offset = (2 ** BUCKET_BITS - 1) * STRIDE + (SLOTS_PER_BUCKET - 1) * SLOT_BYTES + (SLOT_BYTES - 1)
    assert max_offset == SLOTS * SLOT_BYTES - 1 and max_offset < 2 ** 63        # one flat store, 64-bit offsets
    assert BUCKET_BITS <= 64 and 2 ** BUCKET_BITS * STRIDE < 2 ** 63            # i * 1,280 in one word
    assert int(REPORT_CAP) < 2 ** 63 and int(REPORT_CAP) * SLOT_BYTES < 2 ** 63  # report counter and its byte count
    T_BITS = p.bit_length()
    assert 3 + 10 + THETA_BITS + 128 <= T_BITS                   # branch, table, DP and fingerprint bits are disjoint bits of T
    F_SHIFT = 3 + 10 + THETA_BITS                                # F = (T >> 84) mod 2^128
    assert BUCKET_BITS <= 128                                    # bucket index = low 52 bits of F, in F's low word
    words = lambda v: -(-int(v).bit_length() // 64)
    assert words(STUCK) == 2 and words(Decimal(2) ** 140) == 3   # 20/theta and the meter need more than 64 bits
    widths = {"store_bytes_log2": lg(Decimal(SLOTS * SLOT_BYTES)), "max_byte_offset_bits": max_offset.bit_length(),
              "bucket_index_bits": BUCKET_BITS, "bucket_index_times_stride_bits": (2 ** BUCKET_BITS * STRIDE).bit_length(),
              "report_cap_bits": int(REPORT_CAP).bit_length(), "report_bytes_at_cap_bits": (int(REPORT_CAP) * SLOT_BYTES).bit_length(),
              "dp_bits_of_T": THETA_BITS, "dp_bit_range_of_T": [13, 13 + THETA_BITS - 1],
              "fingerprint_bit_range_of_T": [F_SHIFT, F_SHIFT + 127], "T_bits": T_BITS,
              "stuck_walk_threshold_log2": lg(STUCK), "stuck_walk_counter_words": words(STUCK),
              "work_meter_cap_bits_note": "2^139.97 instructions, 3 words"}
    FALSE_CANDIDATES = REPORT_CAP * REPORT_CAP / Decimal(2) ** 129
    assert FALSE_CANDIDATES == Decimal(2) ** -19                # expected false fingerprint candidates at the cap
    assert 25.2 < math.log2(2 ** 16 * 600) < 25.3

    def to_dec(fr):
        return Decimal(fr.numerator) / Decimal(fr.denominator)

    def walk_ledger(SQ, COM, MULT, scale=Decimal(1), window=E, setup=None):
        """cbar and W for one charge vector (identical algebra for every vector); `scale` multiplies the
        average step price (used for the drift rows), `window` is that revision's step bucket outside the
        walk itself (the earlier revisions charged M 20/theta only, this one adds the restart allowance)."""
        cbar = (1 - s) * (1 + Fraction(COM + MULT, C)) + s * Fraction(SQ + COM, C)
        cb = to_dec(cbar) * scale
        return cbar, cb, cb * S * (1 + eps) + cb * window + (SETUP if setup is None else setup)

    cbar_p, cb_p, W_p = walk_ledger(*PARENT, window=E_OLD, setup=SETUP_OLD)
    cbar_r, cb_r, W_r = walk_ledger(*RECORD, window=E_OLD, setup=SETUP_OLD)
    cbar, cb, W = walk_ledger(*OURS)
    # the two earlier ledgers must come out of this same code path unchanged
    if not (abs(lg(W_p) - 125.71472902377629) < 1e-12 and cbar_p == Fraction(9225, 10000)):
        raise RuntimeError('earliest ledger does not reproduce')
    if not (abs(lg(W_r) - 125.7074974545452) < 1e-12 and cbar_r == Fraction(9178875, 10000000)):
        raise RuntimeError('intermediate ledger does not reproduce')
    assert cbar == Fraction(889125, 1000000), cbar
    # step prices
    mul_step = 1 + Fraction(common_charged + extra, C)
    sq_step = Fraction(OURS[0] + common_charged, C)
    assert (mul_step, sq_step) == (Fraction(22401, 20000), Fraction(10085, 20000))
    gap = mul_step - sq_step                    # units a step gains when it is multiply instead of squaring
    # Total-work abort (proof.md 3.4). H1+ allows the squaring fraction to drop by up to DROP_POINTS
    # percentage points; a drop of D points raises the average step price by D/100 * gap, i.e. by the
    # factor 1 + delta(D) = 1 + D/100 * gap / cbar, applied to every step of the run (conservative: it
    # also scales the detection window). The cap is the least thousandth of a bit that covers that.
    DROP_POINTS = Decimal(2)

    def delta(D):
        return Decimal(str(D)) / 100 * to_dec(gap) / to_dec(cbar)

    _, _, W_drop = walk_ledger(*OURS, scale=1 + delta(DROP_POINTS))
    claim_thousandths = int((W_drop.ln() / ln2 * 1000).to_integral_value(rounding="ROUND_CEILING"))
    # exact integer inequality  W_drop <= 2^(claim/1000)  <=>  ceil(W_drop)^1000 <= 2^claim_thousandths
    # with W_drop rounded UP to 40 significant digits, as an exact rational
    W_up_digits = Decimal(W_drop).quantize(Decimal(1).scaleb(int(W_drop.adjusted()) - 40), rounding="ROUND_CEILING")
    fnum, fden = Fraction(W_up_digits).numerator, Fraction(W_up_digits).denominator
    assert fnum ** 1000 <= fden ** 1000 * 2 ** claim_thousandths, "claim is below the drift-covering cap"
    assert fnum ** 1000 > fden ** 1000 * 2 ** (claim_thousandths - 1) , "claim is not the least thousandth"
    CLAIM = Decimal(claim_thousandths) / 1000
    assert str(CLAIM) == "125.682", CLAIM
    Wcap = Decimal(2) ** CLAIM                  # the submitted bound, which is also the work meter's cap
    n = S * (1 + eps)
    cap_over_ledger = Wcap / W - 1
    slack = (Wcap - cb * E - SETUP) / n - cb    # per-step slack of the cap over the nominal ledger
    tol_sq_fraction = slack / to_dec(gap)       # tolerated drop of the squaring fraction below 3/8
    # drop of the squaring fraction at which the cap binds exactly
    assert tol_sq_fraction >= DROP_POINTS / 100

    def success_with_drift(drop_points):
        """Success probability if the squaring fraction drops by `drop_points` percentage points: the
        average step price rises by the factor 1 + delta, so the cap buys n/(1+delta) steps when it binds,
        and the collision law 1 - exp(-steps^2/(2N')) gives the success value."""
        d = delta(drop_points)
        steps = min(n, (Wcap - cb * E * (1 + d) - SETUP) / (cb * (1 + d)))
        return float(1 - (-(steps * steps) / (2 * Npd)).exp())

    drift_table = {p_: success_with_drift(p_) for p_ in ("0.5", "1", "1.5", "2", "2.5", "3", "4", "5")}

    def spent(D):
        return max(0.0, 0.43 - success_with_drift(D))

    lo_, hi_ = Decimal(2), Decimal(20)
    for _ in range(60):                          # drop at which the 0.01 budget of proof.md section 4 is spent
        mid = (lo_ + hi_) / 2
        if spent(mid) < 0.01:
            lo_ = mid
        else:
            hi_ = mid
    return {"measurement_record": meas,
            "fp2_sqr_instr_analytic": FP2SQR, "cyclotomic_square_instr_analytic_listing_not_charged": GS_ANALYTIC,
            "cyclotomic_square_instr_measured_charged": SQ_MEASURED,
            "common_items": common_items, "meter_and_counter_items": meter_items,
            "multiply_only_items": mult_items, "common_counted": common,
            "meter_and_counter_counted": common_charged - common, "common_charged": common_charged,
            "multiply_extra_counted": extra,
            "multiply_step_units": float(mul_step), "squaring_step_units": float(sq_step),
            "multiply_step_instr": common_charged + extra + C, "squaring_step_instr": SQ_MEASURED + common_charged,
            "step_price_gap_units": float(gap),
            "report_store": store, "word_widths": widths, "theta_bits": THETA_BITS,
            "expected_false_fingerprint_candidates_log2": lg(FALSE_CANDIDATES), "report_items": report_items, "report_instr": report_instr,
            "report_cap": float(REPORT_CAP), "report_cap_log2": lg(REPORT_CAP),
            "meter_updates_per_walk": 2, "meter_update_instr_each": METER_UPDATE // 2, "meter_update_instr_per_walk": METER_UPDATE,
            "meter_units_starts_log2": lg(start_units), "meter_units_restarts_log2": lg(restart_units),
            "restart_cap_log2": RESTART_CAP.bit_length() - 1,
            "expected_gaps_log2": lg(EXP_GAPS), "expected_stuck_restarts_log2": lg(EXP_STUCK),
            "expected_long_cycle_restarts_log2": lg(EXP_LONG), "expected_restarts_total_log2": lg(EXP_RESTARTS),
            "expected_restarts_over_cap_log2": lg(EXP_RESTARTS) - 36,
            "stuck_steps_fraction_of_S": float(STUCK_FRAC), "long_cycle_steps_fraction_of_S": float(LONG_FRAC),
            "dp_units_log2": lg(dp_units), "nonwalk_bucket_log2": lg(SETUP),
            "restart_start_units_log2": lg(restart_start_units),
            "earliest_ledger": {"charges": PARENT, "cbar": float(cbar_p), "log2_W": lg(W_p)},
            "intermediate_ledger": {"charges": RECORD, "cbar": float(cbar_r), "log2_W": lg(W_r)},
            "sum_q2": float(sumq2), "log2_N": lg(Decimal(N.numerator) / Decimal(N.denominator)),
            "log2_N_effective": lg(Npd), "log2_S_model_0.43": lg(S), "cbar_units_per_step": float(cb),
            "fruitless_cycle_eps": float(eps), "log2_detection_window_steps": lg(E), "log2_W": lg(W),
            "drift_design_points": float(DROP_POINTS), "drift_factor_at_design": float(delta(DROP_POINTS)),
            "log2_W_at_design_drift": lg(W_drop),
            "submitted_bound": str(CLAIM), "excess_of_cap_over_ledger": float(cap_over_ledger),
            "work_abort_slack_per_step": float(slack), "tolerated_squaring_fraction_drop_points": float(tol_sq_fraction * 100),
            "success_vs_squaring_fraction_drop_points": drift_table,
            "drop_points_where_0.01_is_spent": float(lo_),
            "success_at_model_0.43": float(1 - (-(S * S) / (2 * Npd)).exp()),
            "sensitivity_bits_per_100_instr_on_squaring": lg(walk_ledger(SQ_MEASURED + 100, common_charged, extra)[2]) - lg(W),
            "sensitivity_bits_per_100_instr_every_step": lg(walk_ledger(SQ_MEASURED, common_charged + 100, extra)[2]) - lg(W),
            "sensitivity_bits_at_full_multiplication_price": lg(walk_ledger(C - common_charged, common_charged, extra)[2]) - lg(W) }


out = {}
out["partA_production_BLS12_381"] = part_A(rng, PROD_CHECKS)
if cfg.get("production_steps", 3000) > 0 and PROD_CHECKS > 0:
    out["partC_production_walk"] = part_C(rng, PROD_STEPS, stored=True)
    out["partC_production_walk_integer_rule"] = part_C(rng, PROD_STEPS, stored=False)
out["partB_toy"] = part_B(rng, Q_TOY, TRIALS, R, SQ8, DP_BITS, WALKS, BUDGET)
if PAIR_SAMPLES > 0:
    out["partE_pair_collision_law"] = part_E(rng, Q_TOY, PAIR_SAMPLES, R, SQ8)
out["partD_ledger"] = ledger()
out["seconds"] = round(time.time() - T0, 2)
out["scope"] = ("Parts A and C check exact identities and the production step rule on real BLS12-381 "
                "GT elements; Part B is a toy field; none validates production costs.")
print(json.dumps(out, indent=1))
