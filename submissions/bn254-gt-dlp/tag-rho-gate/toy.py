"""The complete algorithm of proof.md in toy fields, and the tag arithmetic shared with the production
checks (imported by experiment.py). Python standard library only; evidence at toy scale only."""
import math
import random

from kern import Tower, is_prime, mem_order, from_mem, M64

DIG = 51                   # bits per digit of a stored coordinate
TSLOT = 6                  # the walk reads stored coordinate 6 (the real part of e_1, the coefficient of w)


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


class Tags:
    """Tag arithmetic of proof.md section 2 for a tower field: rows, table words, digits, tags.
    scale is the factor of the stored residues (2^256 mod p in production, 1 in the toys)."""

    def __init__(self, T, scale=1):
        self.T, self.q, self.scale = T, T.q, scale % T.q
        self.nd = -(-T.q.bit_length() // DIG)              # digits per coordinate
        self.ndig = 12 * self.nd
        self.E = self.ndig << DIG                          # error bound of the first-level tag (lemma 1)
        self.E2 = 2 * self.ndig                            # error bound of the second-level tag (lemma 2)
        unit = [from_mem([int(i == k) for i in range(12)]) for k in range(12)]
        self.mat = [[mem_order(T.mul(unit[k], unit[i]))[TSLOT] for i in range(12)] for k in range(12)]

    def row(self, Te):
        """v with: coordinate TSLOT of Y * Te = sum_k y_k v_k mod q."""
        t = mem_order(Te)
        return [sum(m * x for m, x in zip(self.mat[k], t)) % self.q for k in range(12)]

    def words(self, Te):
        """Table words G[k][j] = floor(2^64 (v_k 2^(DIG j) mod q) / q)."""
        q = self.q
        return [(((v << (DIG * j)) % q) << 64) // q for v in self.row(Te) for j in range(self.nd)]

    def words2(self, Te):
        """Second table words: the next 64 fraction bits, floor(2^128 u / q) mod 2^64."""
        q = self.q
        return [((((v << (DIG * j)) % q) << 128) // q) & M64 for v in self.row(Te) for j in range(self.nd)]

    def digits(self, Y):
        """Digits of the stored residues of Y."""
        m = (1 << DIG) - 1
        return [((c * self.scale % self.q) >> (DIG * j)) & m for c in mem_order(Y) for j in range(self.nd)]

    def stored_t(self, Y):
        return mem_order(Y)[TSLOT] * self.scale % self.q

    @staticmethod
    def tag(dg, G):
        return sum(d * w for d, w in zip(dg, G)) & M64

    @staticmethod
    def tag2(dg, G2, A1):
        """Second-level tag: the first-level tag plus the high halves of digit x second word."""
        return (A1 + sum((d * w) >> 64 for d, w in zip(dg, G2))) & M64


def threshold(Rm, E, D):
    """A tag A is certified when (Rm A mod 2^64) < this value."""
    return (1 << 64) - Rm * E - (1 << max(0, 64 - D))


def exact_select(t, q, Rm, D):
    """Selector floor(Rm t / q) and the distinguished predicate, from the stored residue t."""
    s = Rm * t // q
    return s, Rm * t - s * q >= q - (q >> D)


class TagRho:
    """proof.md's algorithm with toy sizes: Rm multipliers, windows of L steps, distinguished density 2^-D.
    Table entries are built on first use (same contents as the full table). Runs walk after walk until it
    outputs x or has made `cap` exits from the step loop."""

    def __init__(self, T, rp, g, Rm, L, D, rng):
        self.T, self.rp, self.g, self.Rm, self.L, self.D, self.rng = T, rp, g, Rm, L, D, rng
        self.tg = Tags(T)
        self.thr = threshold(Rm, self.tg.E, D)
        self.thr2 = threshold(Rm, self.tg.E2, D)

    def setup(self, z):
        T, rp, rng = self.T, self.rp, self.rng
        self.z = z
        self.cd = [(rng.randrange(rp), rng.randrange(rp)) for _ in range(self.Rm)]
        self.M = [T.mul(T.pow(self.g, c), T.pow(z, d)) for c, d in self.cd]
        self.root = [T.one, 0, 0, self.tg.words(T.one), 0, [None] * self.Rm, self.tg.words2(T.one)]   # T_e, alpha, beta, words, size, children, second words
        self.gz = []                                           # g^(2^i) and z^(2^i) for the walk starts
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
                 e[4] + 1, [None] * self.Rm, self.tg.words2(Te)]
            e[5][s] = c; self.entries += 1
        return c

    def run(self, cap, check=0, xlog=None):
        """xlog (the target's logarithm, known to a test driver only) is used for one measurement: it labels
        every visited element by its logarithm so that revisited elements can be counted."""
        T, rp, q, Rm, tg = self.T, self.rp, self.T.q, self.Rm, self.tg
        store = {}
        steps = exits = starts = reports = unc = capped = checked = bad = revisits = lev2 = 0
        seen = set()
        wcap = 128 << self.D                                   # steps of one walk
        res = lambda x: dict(x=x, steps=steps, exits=exits, starts=starts, reports=reports, uncertified=unc,
                             capped_walks=capped, checked=checked, check_failures=bad, entries=self.entries,
                             revisits=revisits, second_level=lev2)
        while exits < cap:
            starts += 1
            a, b = self.rng.randrange(rp), self.rng.randrange(rp)
            Y0 = T.one                                         # g^a z^b from the stored squares
            for tab, ex in zip(self.gz, (a, b)):
                for i in range(ex.bit_length()):
                    if (ex >> i) & 1:
                        Y0 = T.mul(Y0, tab[i])
            if Y0 == T.one:
                continue
            dg, e, n = tg.digits(Y0), self.root, 0
            while True:
                if n >= wcap:
                    capped += 1; break
                if exits >= cap:
                    return res(None)
                if xlog is not None:                           # measurement only: was this element visited before?
                    lg = (a + e[1] + xlog * (b + e[2])) % rp
                    if lg in seen:
                        revisits += 1
                    else:
                        seen.add(lg)
                A = tg.tag(dg, e[3])
                cert = (Rm * A) & M64 < self.thr               # first level
                if not cert:                                   # second level: 64 more fraction bits
                    lev2 += 1
                    A = tg.tag2(dg, e[6], A)
                    cert = (Rm * A) & M64 < self.thr2
                if cert:                                       # certified: selector from the tag
                    s = (Rm * A) >> 64
                    if checked < check:                        # against the exact definition
                        se, de = exact_select(tg.stored_t(T.mul(Y0, e[0])), q, Rm, self.D)
                        checked += 1; bad += se != s or de
                    if e[4] < self.L:
                        e = self.child(e, s)
                    else:                                      # window end: re-base
                        exits += 1
                        Y0 = T.mul(Y0, e[0]); a, b = (a + e[1]) % rp, (b + e[2]) % rp
                        dg, e = tg.digits(Y0), self.child(self.root, s)
                else:                                          # not certified: re-base, exact evaluation
                    exits += 1; unc += 1
                    if e is not self.root:
                        Y0 = T.mul(Y0, e[0]); a, b = (a + e[1]) % rp, (b + e[2]) % rp
                        dg = tg.digits(Y0)
                    t = tg.stored_t(Y0)
                    s, dist = exact_select(t, q, Rm, self.D)
                    if dist:
                        reports += 1
                        hit = store.get(Y0)                    # the store key is the whole element
                        if hit is None:
                            store[Y0] = (a, b)
                        elif hit[1] != b:
                            x = (a - hit[0]) * pow((hit[1] - b) % rp, -1, rp) % rp
                            if T.pow(self.g, x) == self.z:
                                return res(x)
                        break
                    e = self.child(self.root, s)
                steps += 1; n += 1
        return res(None)


def model_steps(rp, Rm):
    """Model mean of the steps before the first collision of an Rm-adding walk: sqrt(pi r Rm / (2 (Rm - 1)))."""
    return math.sqrt(math.pi * rp * Rm / (2 * (Rm - 1)))


def toy_field(family="bn", bn_x=-15, rbits=22, qbits=40):
    q, rp, a0 = bn_toy(bn_x) if family == "bn" else find_toy(rbits, qbits)
    T = Tower(q, a0)
    rng = random.Random(10 ** 6)
    while True:
        y = tuple(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(3)) for _ in range(2))
        g = T.pow(y, (q ** 12 - 1) // rp)
        if g != T.one: break
    return q, rp, a0, T, g


def toy(seed, instances, Rm, L, D=4, check=100):
    q, rp, a0, T, g = toy_field()
    m0 = model_steps(rp, Rm)
    cap = 64 * (math.isqrt(rp) + 1)                           # toy exit cap
    rows = []
    for i in range(instances):
        rr = random.Random((seed << 20) + (Rm << 12) + (L << 4) + i)
        x = rr.randrange(1, rp)
        alg = TagRho(T, rp, g, Rm, L, D, rr)
        alg.setup(T.pow(g, x))
        res = alg.run(cap, check, x)
        res["correct"] = res["x"] == x if res["x"] is not None else None
        rows.append(res)
    ratio = [o["steps"] / m0 for o in rows]
    mean = sum(ratio) / instances
    se = math.sqrt(sum((v - mean) ** 2 for v in ratio) / (instances - 1) / instances)
    steps = sum(o["steps"] for o in rows)
    return {"family": "bn", "q": q, "r_prime": rp, "xi": f"{a0}+u", "multipliers": Rm, "window": L, "dp_bits": D,
            "digits_per_coordinate": alg.tg.nd, "instances": instances, "seed": seed, "exit_cap": cap,
            "model_mean_steps": round(m0, 1),
            "correct_outputs": sum(1 for o in rows if o["correct"]),
            "wrong_outputs": sum(1 for o in rows if o["correct"] is False),
            "runs_at_exit_cap": sum(1 for o in rows if o["x"] is None),
            "mean_steps_over_model": round(mean, 4), "standard_error": round(se, 4),
            "exits_per_step": round(sum(o["exits"] for o in rows) / steps, 4),
            "second_level_per_step": round(sum(o["second_level"] for o in rows) / steps, 5),
            "uncertified_per_step": round(sum(o["uncertified"] for o in rows) / steps, 5),
            "reports_per_step": round(sum(o["reports"] for o in rows) / steps, 5),
            "revisited_elements_per_run_in_walk_lengths": round(sum(o["revisits"] for o in rows) / instances / (1 << D), 3),
            "selector_checks": sum(o["checked"] for o in rows),
            "selector_check_failures": sum(o["check_failures"] for o in rows),
            "capped_walks": sum(o["capped_walks"] for o in rows),
            "note": "at 8 to 32 instances the mean has a standard error of 0.08 to 0.2; the offline runs are the evidence"}
