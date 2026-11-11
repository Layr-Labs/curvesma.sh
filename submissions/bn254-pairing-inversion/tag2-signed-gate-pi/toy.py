"""The complete algorithm of proof.md in toy fields, and the two-level tag arithmetic shared with the
production checks (imported by experiment.py). Python standard library only; evidence at toy scale only.
The word size wb is 64 in production; toys use a smaller wb so that both tag levels occur at toy sizes."""
import math
import random

from kern import Tower, is_prime, mem_order, from_mem, M64

TSLOT = 6                  # the walk reads stored coordinate 6: the real part of the coefficient of w


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
    """BN-family toy: q = p(x), r' = r(x) from the BN polynomials, so q ~ r' as in production."""
    q = 36 * x ** 4 + 36 * x ** 3 + 24 * x ** 2 + 6 * x + 1
    rp = 36 * x ** 4 + 36 * x ** 3 + 18 * x ** 2 + 6 * x + 1
    assert is_prime(q) and is_prime(rp) and q % 12 == 7 and all(pow(q, k, rp) != 1 for k in range(1, 12))
    for a0 in range(1, 1000):
        T = Tower.__new__(Tower); T.q = q
        xi = (a0, 1)
        if T.f2pow(xi, (q * q - 1) // 2) != (1, 0) and T.f2pow(xi, (q * q - 1) // 3) != (1, 0): break
    return q, rp, a0


def recode(c, b, nd):
    """Balanced digits: d_i in [-2^(b-1), 2^(b-1)) for i < nd - 1, last digit >= 0, sum d_i 2^(b i) = c."""
    m, half, out, carry = (1 << b) - 1, 1 << (b - 1), [], 0
    for i in range(nd - 1):
        t = ((c >> (b * i)) & m) + carry
        carry = (t + half) >> b
        out.append(t - (carry << b))
    out.append((c >> (b * (nd - 1))) + carry)
    return out


def threshold(Rm, E, D, wb=64):
    """A tag A with interval width E is certified when (Rm A mod 2^wb) < this value."""
    return (1 << wb) - Rm * E - (1 << max(0, wb - D))


def exact_select(t, q, Rm, D):
    """Selector floor(Rm t / q) and the distinguished predicate, from the stored residue t."""
    s = Rm * t // q
    return s, Rm * t - s * q >= q - (q >> D)


class Tags:
    """Two-level tag arithmetic of proof.md section 2. scale: factor of the stored residues (2^256 mod p in
    production, 1 in the toys). b: digit bits; wb: word bits (64 in production)."""

    def __init__(self, T, scale=1, b=51, Rm=6, D=90, wb=64):
        self.T, self.q, self.scale, self.b, self.Rm, self.D, self.wb = T, T.q, scale % T.q, b, Rm, D, wb
        self.M = (1 << wb) - 1
        self.nd = -(-T.q.bit_length() // b)                # digits per coordinate
        self.ndig = 12 * self.nd
        self.top = ((T.q - 1) >> (b * (self.nd - 1))) + 1  # bound on the last digit of a coordinate
        self.Smax = 12 * ((self.nd - 1) * (1 << (b - 1)) + self.top)   # bound on the sum of |digits|
        self.E2 = 2 * self.ndig + 2                        # level-2 interval width
        self.thr2 = threshold(Rm, self.E2, D, wb)
        unit = [from_mem([int(i == k) for i in range(12)]) for k in range(12)]
        self.mat = [[mem_order(T.mul(unit[k], unit[i]))[TSLOT] for i in range(12)] for k in range(12)]

    def row(self, Te):
        """v with: coordinate TSLOT of Y * Te = sum_k y_k v_k mod q."""
        t = mem_order(Te)
        return [sum(m * x for m, x in zip(self.mat[k], t)) % self.q for k in range(12)]

    def words(self, Te):
        """G[k,j] = floor(2^wb u / q), V[k,j] = floor(2^(2wb-1) u / q) mod 2^(wb-1), u = v_k 2^(b j) mod q."""
        q, wb, G, V = self.q, self.wb, [], []
        for v in self.row(Te):
            for j in range(self.nd):
                u = (v << (self.b * j)) % q
                G.append((u << wb) // q)
                V.append(((u << (2 * wb - 1)) // q) & ((1 << (wb - 1)) - 1))
        return G, V

    def digits(self, Y):
        out = []
        for c in mem_order(Y):
            out += recode(c * self.scale % self.q, self.b, self.nd)
        return out

    def window(self, dg):
        """Per-window constants: N = sum of |negative digits|, E1 = sum |d| + 1, THR1."""
        N = sum(-d for d in dg if d < 0)
        E1 = N + sum(d for d in dg if d > 0) + 1
        return {"N": N, "E1": E1, "THR1": threshold(self.Rm, E1, self.D, self.wb)}

    def stored_t(self, Y):
        return mem_order(Y)[TSLOT] * self.scale % self.q

    def level1(self, dg, G, w):
        return (sum(d * g for d, g in zip(dg, G)) - w["N"]) & self.M

    def level2(self, dg, V, A1, w):
        return (A1 + 2 * sum((d * v) >> self.wb for d, v in zip(dg, V)) + w["N"] - 1) & self.M

    def select(self, dg, G, V, w):
        """(s, level): level 1 or 2 when certified, (None, 3) when the element must be evaluated exactly."""
        Rm, M, wb = self.Rm, self.M, self.wb
        A1 = self.level1(dg, G, w)
        if (Rm * A1) & M < w["THR1"]:
            return (Rm * A1) >> wb, 1
        A2 = self.level2(dg, V, A1, w)
        if (Rm * A2) & M < self.thr2:
            return (Rm * A2) >> wb, 2
        return None, 3

    def q1_bound(self):
        """Level-1 uncertified share for a uniform tag, with the deterministic bound on E1."""
        return (self.Rm * (self.Smax + 1) + (1 << max(0, self.wb - self.D)) + 1) / 2 ** self.wb


def degenerate(T, Ms):
    """The multiplier check of proof.md step 2: some M_a^(q^(2 i)) = M_b with (a, i) != (b, 0), i < 6."""
    imgs = {}
    for a, M in enumerate(Ms):
        Y = M
        for i in range(6):
            k_ = key(Y)
            if k_ in imgs:
                return True
            imgs[k_] = (a, i)
            Y = T.frob_k(Y, 2)
    return False


def key(Y):
    """Store key: the whole element."""
    return tuple(mem_order(Y))


class TagRho:
    """proof.md's algorithm with toy sizes: Rm multipliers, windows of L steps, distinguished density 2^-D.
    Table entries are built on first use (same contents as the full table). Runs walk after walk until it
    outputs x or has made `cap` exits from the step loop. The store is keyed by the whole element."""

    def __init__(self, T, rp, g, Rm, L, D, rng, b, wb):
        self.T, self.rp, self.g, self.Rm, self.L, self.D, self.rng = T, rp, g, Rm, L, D, rng
        self.tg = Tags(T, 1, b, Rm, D, wb)

    def setup(self, z):
        T, rp, rng = self.T, self.rp, self.rng
        self.z = z
        for draw in range(16):                         # proof.md step 2: redraw degenerate multiplier sets
            self.cd = [(rng.randrange(rp), rng.randrange(rp)) for _ in range(self.Rm)]
            self.M = [T.mul(T.pow(self.g, c), T.pow(z, d)) for c, d in self.cd]
            if not degenerate(T, self.M):
                break
        self.draws = draw + 1
        self.root = [T.one, 0, 0, self.tg.words(T.one), 0, [None] * self.Rm]
        self.gz = []
        for base in (self.g, z):
            tab = [base]
            for _ in range(rp.bit_length() - 1):
                tab.append(T.mul(tab[-1], tab[-1]))
            self.gz.append(tab)
        self.entries = 1

    def child(self, e, s):
        c = e[5][s]
        if c is None:
            Te = self.T.mul(e[0], self.M[s])
            c = [Te, (e[1] + self.cd[s][0]) % self.rp, (e[2] + self.cd[s][1]) % self.rp, self.tg.words(Te),
                 e[4] + 1, [None] * self.Rm]
            e[5][s] = c; self.entries += 1
        return c

    def start(self):
        T, rp = self.T, self.rp
        a, b = self.rng.randrange(rp), self.rng.randrange(rp)
        Y0 = T.one
        for tab, ex in zip(self.gz, (a, b)):
            for i in range(ex.bit_length()):
                if (ex >> i) & 1:
                    Y0 = T.mul(Y0, tab[i])
        return a, b, Y0

    def run(self, cap, check=0, trace=None):
        T, rp, q, Rm, tg = self.T, self.rp, self.T.q, self.Rm, self.tg
        store = {}
        st = dict(steps=0, exits=0, starts=0, reports=0, level1_fail=0, level2=0, level3=0, capped_walks=0,
                  checked=0, check_failures=0, entries=0)
        wcap = 128 << self.D
        while st["exits"] < cap:
            st["starts"] += 1
            a, b, Y0 = self.start()
            if Y0 == T.one:
                continue
            dg = tg.digits(Y0); w = tg.window(dg); e, n = self.root, 0
            while True:
                if n >= wcap:
                    st["capped_walks"] += 1; break
                if st["exits"] >= cap:
                    st["entries"] = self.entries
                    return dict(st, x=None)
                G, V = e[3]
                s, lev = tg.select(dg, G, V, w)
                st["level1_fail"] += lev > 1; st["level2"] += lev == 2; st["level3"] += lev == 3
                if lev < 3:
                    if st["checked"] < check:                  # against the exact definition
                        se, de = exact_select(tg.stored_t(T.mul(Y0, e[0])), q, Rm, self.D)
                        st["checked"] += 1; st["check_failures"] += se != s or de
                    if trace is not None:
                        trace.append(key(T.mul(Y0, e[0])))
                    if e[4] < self.L:
                        e = self.child(e, s)
                    else:                                      # window end: re-base
                        st["exits"] += 1
                        Y0 = T.mul(Y0, e[0]); a, b = (a + e[1]) % rp, (b + e[2]) % rp
                        dg = tg.digits(Y0); w = tg.window(dg); e = self.child(self.root, s)
                else:                                          # not certified: re-base, exact evaluation
                    st["exits"] += 1
                    if e is not self.root:
                        Y0 = T.mul(Y0, e[0]); a, b = (a + e[1]) % rp, (b + e[2]) % rp
                        dg = tg.digits(Y0); w = tg.window(dg)
                    if trace is not None:
                        trace.append(key(Y0))
                    s, dist = exact_select(tg.stored_t(Y0), q, Rm, self.D)
                    if dist:
                        st["reports"] += 1
                        k_ = key(Y0)
                        hit = store.get(k_)
                        if hit is None:
                            store[k_] = (a, b)
                        elif hit[1] != b:
                            x = (a - hit[0]) * pow((hit[1] - b) % rp, -1, rp) % rp
                            if T.pow(self.g, x) == self.z:
                                st["entries"] = self.entries
                                return dict(st, x=x)
                        break
                    e = self.child(self.root, s)
                st["steps"] += 1; n += 1
        st["entries"] = self.entries
        return dict(st, x=None)


class RefRho(TagRho):
    """Reference: the same walk evaluated by full products and the exact selector (no tags, no table)."""

    def run(self, cap=None, check=0, trace=None):
        T, rp, q, Rm = self.T, self.rp, self.T.q, self.Rm
        store, steps, starts, reports = {}, 0, 0, 0
        wcap = 128 << self.D
        while True:
            starts += 1
            a, b, Y = self.start()
            if Y == T.one:
                continue
            n = 0
            while n < wcap:
                if trace is not None:
                    trace.append(key(Y))
                s, dist = exact_select(mem_order(Y)[TSLOT] % q, q, Rm, self.D)
                if dist:
                    reports += 1
                    hit = store.get(key(Y))
                    if hit is None:
                        store[key(Y)] = (a, b)
                    elif hit[1] != b:
                        x = (a - hit[0]) * pow((hit[1] - b) % rp, -1, rp) % rp
                        if T.pow(self.g, x) == self.z:
                            return dict(x=x, steps=steps, starts=starts, reports=reports)
                    break
                Y = T.mul(Y, self.M[s]); a, b = (a + self.cd[s][0]) % rp, (b + self.cd[s][1]) % rp
                steps += 1; n += 1


def model_steps(rp, Rm):
    """Model mean of the steps before the first collision of an Rm-adding walk: sqrt(pi r Rm / (2 (Rm - 1)))."""
    return math.sqrt(math.pi * rp * Rm / (2 * (Rm - 1)))


def toy_field(family="bn", bn_x=-15, rbits=26, qbits=44):
    q, rp, a0 = bn_toy(bn_x) if family == "bn" else find_toy(rbits, qbits)
    T = Tower(q, a0)
    rng = random.Random(10 ** 6)
    while True:
        y = tuple(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(3)) for _ in range(2))
        g = T.pow(y, (q ** 12 - 1) // rp)
        if g != T.one: break
    return q, rp, a0, T, g


def toy_params(q):
    """Toy digit and word sizes: two digits per coordinate, and wb chosen so that the level-1 share
    6 Smax / 2^wb is between 1% and 2%, as in production (2.1%)."""
    b = -(-q.bit_length() // 2)
    nd = 2
    smax = 12 * ((nd - 1) * (1 << (b - 1)) + ((q - 1) >> b) + 1)
    wb = (6 * smax * 50).bit_length()
    return b, wb


def run_instance(T, rp, g, Rm, L, D, sd, check, b, wb, cap):
    rr = random.Random(sd)
    x = rr.randrange(1, rp)
    alg = TagRho(T, rp, g, Rm, L, D, rr, b, wb)
    alg.setup(T.pow(g, x))
    res = alg.run(cap, check)
    res["correct"] = res["x"] == x if res["x"] is not None else None
    return res


def toy(seed, instances, Rm, L, D=4, check=100, bn_x=-15, ref=4):
    """Mode `toy` of experiment.py: instances of the complete algorithm, and for the first `ref` instances
    the element sequence of the tag walk compared with the reference walk (same coins)."""
    q, rp, a0, T, g = toy_field("bn", bn_x)
    b, wb = toy_params(q)
    m0 = model_steps(rp, Rm)
    cap = 64 * (math.isqrt(rp) + 1)
    rows, ref_bad = [], 0
    for i in range(instances):
        sd = (seed << 20) + (Rm << 12) + (L << 4) + i
        rows.append(run_instance(T, rp, g, Rm, L, D, sd, check, b, wb, cap))
        if i < ref:
            tr1, tr2 = [], []
            for cls, tr in ((TagRho, tr1), (RefRho, tr2)):
                rr = random.Random(sd); x = rr.randrange(1, rp)
                alg = cls(T, rp, g, Rm, L, D, rr, b, wb); alg.setup(T.pow(g, x))
                res = alg.run(cap, 0, tr)
            ref_bad += tr1 != tr2
    ratio = [o["steps"] / m0 for o in rows]
    mean = sum(ratio) / instances
    se = math.sqrt(sum((v - mean) ** 2 for v in ratio) / (instances - 1) / instances)
    steps = sum(o["steps"] for o in rows) or 1
    return {"family": "bn", "q": q, "r_prime": rp, "xi": f"{a0}+u", "multipliers": Rm, "window": L, "dp_bits": D,
            "digit_bits": b, "word_bits": wb, "instances": instances, "seed": seed, "exit_cap": cap,
            "model_mean_steps": round(m0, 1),
            "correct_outputs": sum(1 for o in rows if o["correct"]),
            "wrong_outputs": sum(1 for o in rows if o["correct"] is False),
            "runs_at_exit_cap": sum(1 for o in rows if o["x"] is None),
            "mean_steps_over_model": round(mean, 4), "standard_error": round(se, 4),
            "level1_not_certified_per_step": round(sum(o["level1_fail"] for o in rows) / steps, 5),
            "level2_certified_per_step": round(sum(o["level2"] for o in rows) / steps, 5),
            "exact_evaluations_per_step": round(sum(o["level3"] for o in rows) / steps, 5),
            "reports_per_step": round(sum(o["reports"] for o in rows) / steps, 5),
            "selector_checks": sum(o["checked"] for o in rows),
            "selector_check_failures": sum(o["check_failures"] for o in rows),
            "reference_walk_comparisons": min(ref, instances), "reference_walk_mismatches": ref_bad,
            "capped_walks": sum(o["capped_walks"] for o in rows),
            "note": "at 8 to 32 instances the mean has a standard error of 0.08 to 0.2; the offline runs are the evidence"}


OFFLINE_SEED = 5 * 10 ** 7     # instance i of a row: random.Random(OFFLINE_SEED + 1000 i + 10 window + digits per coordinate)


def offline_instance(family, L, i, D=6, check=100):
    """One instance of the offline runs in toy_offline.json (the target x is drawn first from the seed)."""
    q, rp, a0, T, g = toy_field(family, -41)
    b, wb = toy_params(q)
    sd = OFFLINE_SEED + 1000 * i + 10 * L + 2
    res = run_instance(T, rp, g, 6, L, D, sd, check, b, wb, 64 * (math.isqrt(rp) + 1))
    return {k: res[k] for k in ("steps", "exits", "starts", "reports", "level1_fail", "level2", "level3",
                                "capped_walks", "checked", "check_failures", "correct")}


def census(bn_x=-15, slots=(1, TSLOT), dps=(4, 6, 8, 10)):
    """Every element of the order-r' subgroup of a BN-family toy field, read at each coordinate in slots:
    selector counts floor(6 t / q), distinguished counts for each D, and the number of distinct values of t."""
    q, rp, a0, T, g = toy_field("bn", bn_x)
    out = {"q": q, "r_prime": rp, "slots": {}}
    acc = {s: {"selector": [0] * 6, "distinguished": {str(D): 0 for D in dps}, "seen": bytearray(q // 8 + 1)}
           for s in slots}
    Y = T.one
    for k in range(rp):
        m = mem_order(Y)
        for s in slots:
            t, a = m[s], acc[s]
            sel = 6 * t // q
            a["selector"][sel] += 1
            for D in dps:
                if 6 * t - sel * q >= q - (q >> D):
                    a["distinguished"][str(D)] += 1
            a["seen"][t >> 3] |= 1 << (t & 7)
        Y = T.mul(Y, g)
    assert Y == T.one
    for s in slots:
        a = acc[s]
        out["slots"][str(s)] = {"selector": a["selector"], "distinguished": a["distinguished"],
                                "distinct_values": sum(bin(b).count("1") for b in a["seen"])}
    return out
