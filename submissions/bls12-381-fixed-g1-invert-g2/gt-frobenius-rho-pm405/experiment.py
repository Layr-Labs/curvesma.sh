"""Fixed-G1 pairing inversion into G2 by Frobenius-class rho with a squaring branch: identities, ledger, toy runs (parts A-F, S, D)."""
import json, math, os, random, runpy, sys, time
from fractions import Fraction

T0 = time.time()
cfg = json.load(sys.stdin)
Q_TOY = cfg.get("q", 67)
TRIALS = cfg.get("trials", 20)
R = cfg.get("R", 64)
BRANCH_BITS = cfg.get("branch_bits", 6)
SQ_CUT = cfg.get("sq_cut", 37)
DP_BITS = cfg.get("dp_bits", 3)
WALKS = cfg.get("walks", 4)
SEED = cfg.get("seed", 1)
PROD_CHECKS = cfg.get("production_checks", 2)
PROD_STEPS = cfg.get("production_steps", 3000)
PAIR_SAMPLES = cfg.get("pair_samples", 0)
rng = random.Random(SEED)

SAMPLE_R_MAX = 256
SAMPLE_BLOCK_INSTR = 2300
SAMPLE_BLOCK_ITEMS = {"20 rounds x 4 quarter-rounds x 12 instructions (4 add, 4 xor, 4 rotate)": 960,
                      "16-word feed-forward add": 16,
                      "load and store the 16-word state": 32,
                      "64-bit block counter, adds/adc": 2,
                      "address set-up and loop control": 10,
                      "register spills, prologue and return": 250}
assert sum(SAMPLE_BLOCK_ITEMS.values()) == 1270 < SAMPLE_BLOCK_INSTR
SAMPLE_ATTEMPT_INSTR = SAMPLE_BLOCK_INSTR + 20
assert SAMPLE_ATTEMPT_INSTR == 2320
SAMPLE_MAX_INSTR = SAMPLE_R_MAX * SAMPLE_ATTEMPT_INSTR


def sample_below(rng, m, r_max=SAMPLE_R_MAX):
    if m <= 1:
        return 0
    bits = max(1, m.bit_length())
    for _ in range(r_max):
        x = rng.getrandbits(bits)
        if x < m:
            return x
    return 1


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

        def minus(t, u):
            d = self.f2sub(t, u)
            return self.f2add(t, self.f2add(d, d))

        def plus(t, u):
            s = self.f2add(t, u)
            return self.f2add(t, self.f2add(s, s))

        return (minus(A0, c0), plus(xC1, c1), minus(B0, c2), plus(A1, c3), minus(C0, c4), plus(B1, c5))

    def rand_gt(self, rng, rp):
        hc = (self.q ** 12 - 1) // rp
        while True:
            x = tuple((sample_below(rng, self.q), sample_below(rng, self.q)) for _ in range(6))
            y = self.pow(x, hc)
            if y != self.one:
                return y


p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
MR384 = pow(2, 384, p)
BRANCH_BITS_PROD, SQ_CUT_PROD = 6, 36
MUL_CORE, SQ_CORE, COMMON_CHARGED, MULT_EXTRA = 14563, 7154, 262, 2119
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
        x = tuple((sample_below(rng, p), sample_below(rng, p)) for _ in range(6))
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
        ok["trace_inv"] &= len({c[0][0] for c in conj}) == 1
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


def part_S(checks=20000):
    """Exercise the bounded sampler of proof.md 0-iv: exact small-modulus reduction, R_MAX cap, residue 1."""
    class Src(object):
        def __init__(self, vals=None, reject=False):
            self.vals, self.reject, self.calls = list(vals or []), reject, 0

        def getrandbits(self, k):
            self.calls += 1
            return (1 << k) - 1 if self.reject or not self.vals else self.vals.pop(0)

    seq = Src([0, 1, 2, 3, 4, 5, 3])
    accepted = [sample_below(seq, 5) for _ in range(6)]
    exact = accepted == [0, 1, 2, 3, 4, 3] and seq.calls == 7
    ra = Src(reject=True)
    v = sample_below(ra, r)
    rng_s = random.Random(SEED ^ 0x5A5A)
    vals = [sample_below(rng_s, r) for _ in range(checks)]
    out = {"r_max_attempts": SAMPLE_R_MAX, "attempt_instr": SAMPLE_ATTEMPT_INSTR,
           "max_instr_per_sample": SAMPLE_MAX_INSTR, "max_units_per_sample": SAMPLE_MAX_INSTR / 20000,
           "acceptance_at_least": "r / 2^%d >= 1/2" % r.bit_length(), "exhaustion_rule": "return 1",
           "exhaustion_probability_log2_at_most": -SAMPLE_R_MAX, "exhaustion_allowance_log2_charged": -218.99,
           "chacha20_block_item_instr": sum(SAMPLE_BLOCK_ITEMS.values()),
           "chacha20_block_charged_instr": SAMPLE_BLOCK_INSTR,
           "exact_reduction_small_modulus": exact, "all_reject_returns": v, "all_reject_attempts": ra.calls,
           "seeded_draws_checked": checks, "seeded_all_in_range": all(0 <= x < r for x in vals)}
    out["all_checks_pass"] = bool(exact and v == 1 and ra.calls == SAMPLE_R_MAX and out["seeded_all_in_range"])
    assert out["all_checks_pass"]
    return out


def part_C(rng, steps, rprod=16, verify=150, stored=False):
    """Production step rule in the real GT: twisted-table multiply / squaring, bits from T = Re(c_0).

    stored=True reads the six limbs as the pinned blst kernels store them (Montgomery form, the rule of
    proof.md section 3); stored=False reads ordinary residues. Both run the same step and class checks.
    """
    P = Tower(p)
    hcof = (p ** 12 - 1) // r
    lam = [pow(p, k, r) for k in range(12)]
    g = P.pow(tuple((sample_below(rng, p), sample_below(rng, p)) for _ in range(6)), hcof)
    x0 = 1 + sample_below(rng, (1 << 40) - 1)
    z0 = P.pow(g, x0)
    ah = [1 + sample_below(rng, (1 << 40) - 1) for _ in range(rprod)]
    bh = [1 + sample_below(rng, (1 << 40) - 1) for _ in range(rprod)]
    Mh = [P.mul(P.pow(g, ah[i]), P.pow(z0, bh[i])) for i in range(rprod)]
    Tab = [[P.frobk(Mh[i], k) for k in range(12)] for i in range(rprod)]
    a0, b0 = sample_below(rng, r), sample_below(rng, r)
    y = P.mul(P.pow(g, a0), P.pow(z0, b0))
    al, be = a0, b0
    n_sq = n_mul = bad_class = 0
    bins_sq = [0] * (1 << BRANCH_BITS)
    bins_h = [0] * 16
    two_cycles = 0
    idx_differs = 0
    hist = []
    first_seen, last_new, tail = {}, 0, []
    for t in range(steps):
        u = (y[0][0] * MR384) % p if stored else y[0][0]
        if u not in first_seen:
            first_seen[u], last_new = t, t
        tail = (tail + [u])[-64:]
        if len(hist) >= 2 and hist[-2][0] == u:
            two_cycles += 1
            _, yy, aa, bb = min(hist[-1], (u, y, al, be), key=lambda s: s[0])
            y, al, be = P.cyc_sqr(yy), 2 * aa % r, 2 * bb % r
            hist = []
            continue
        hist.append((u, y, al, be))
        if len(hist) > 2:
            hist.pop(0)
        br = u & ((1 << BRANCH_BITS) - 1)
        bins_sq[br] += 1
        if br < SQ_CUT:
            y2 = P.cyc_sqr(y)
            if t < verify:
                kap = min(range(12), key=lambda i: P.keys(y)[i])
                bad_class += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), P.frobk(y, kap))
            y, al, be = y2, 2 * al % r, 2 * be % r
            n_sq += 1
        else:
            kf = P.keys(y)
            kf_int = kf
            if stored:
                kf = [(v * MR384) % p for v in kf]
            kap = kf.index(min(kf))
            idx_differs += (kap != kf_int.index(min(kf_int)))
            h = (u >> BRANCH_BITS) % rprod
            bins_h[(u >> BRANCH_BITS) & 15] += 1
            idx = (12 - kap) % 12
            y2 = P.mul(y, Tab[h][idx])
            if t < verify:
                bad_class += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), Mh[h])
            y = y2
            al = (al + ah[h] * lam[idx]) % r
            be = (be + bh[h] * lam[idx]) % r
            n_mul += 1
    exp_ok = y == P.mul(P.pow(g, al), P.pow(z0, be))
    exp_sq = sum(bins_sq[:SQ_CUT])
    n_br = sum(bins_h)
    n_all = sum(bins_sq)
    nb = 1 << BRANCH_BITS
    chi_sq = sum((c - n_all / nb) ** 2 / (n_all / nb) for c in bins_sq)
    chi_h = sum((c - n_br / 16) ** 2 / (n_br / 16) for c in bins_h) if n_br else None
    s_prod = SQ_CUT / nb
    return {"steps": steps, "R_table": rprod, "squaring_steps": n_sq, "multiply_steps": n_mul,
            "key_rule": "stored Montgomery limbs (x*2^384 mod p)" if stored else "ordinary integers in [0,p)",
            "min_index_differs_from_integer_rule": idx_differs,
            "squaring_fraction": round(n_sq / max(1, n_all), 4), "target_fraction": s_prod,
            "branch_rule_decisions": n_all,
            "binomial_sd_of_the_fraction_at_target": round(math.sqrt(s_prod * (1 - s_prod) / max(1, n_all)), 4),
            "z_score_of_the_fraction": round((n_sq / max(1, n_all) - s_prod) / math.sqrt(s_prod * (1 - s_prod) / max(1, n_all)), 2),
            "explicit_class_check_failures_in_first_%d" % verify: int(bad_class),
            "exponent_invariant_y_equals_g^al_z^be": exp_ok,
            "chi2_low%d_bits_of_Re_c0_%ddof" % (BRANCH_BITS, nb - 1): round(chi_sq, 2),
            "chi2_next4_bits_of_Re_c0_15dof": round(chi_h, 2) if chi_h is not None else None,
            "two_cycle_events": two_cycles,
            "distinct_T_values_visited": len(first_seen), "last_step_with_a_new_T": last_new,
            "period_of_T_over_the_last_64_steps": next((k for k in range(1, len(tail))
                                                        if all(tail[i] == tail[i - k] for i in range(k, len(tail)))), None)}


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


def heapsort_level_bound(n):
    """Worst-case sift-down levels of heapsort on n records (proof.md 0-ii): at most n for Floyd's
    construction plus sum_{k=2..n} floor(log2 k) for the n - 1 extractions, summed in power-of-two blocks."""
    tot, j = 0, 1
    while (1 << j) <= n:
        tot += j * (min(n, (1 << (j + 1)) - 1) - (1 << j) + 1)
        j += 1
    return n + tot


def heapsort(a):
    """In-place heapsort of the list a (Floyd construction, then n - 1 extractions); returns the number of
    sift-down levels executed (a level = one pass that compares a node with its children)."""
    n, levels = len(a), 0

    def sift(i, m):
        nonlocal levels
        while 2 * i + 1 < m:
            levels += 1
            c = 2 * i + 1
            if c + 1 < m and a[c + 1] > a[c]:
                c += 1
            if a[c] > a[i]:
                a[i], a[c] = a[c], a[i]
                i = c
            else:
                break
    for i in range(n // 2 - 1, -1, -1):
        sift(i, n)
    for m in range(n - 1, 0, -1):
        a[0], a[m] = a[m], a[0]
        sift(0, m)
    return levels


def match_log(log, verify, rp, pair_budget):
    """The matcher of proof.md 0-ii on an append-only log of (F, A, B) records: heapsort of the whole
    record, one adjacent comparison per record to drop exact duplicates, then every pair inside each
    run of equal F under a hard pair budget, halting at the first pair whose candidate verifies."""
    levels = heapsort(log)
    if levels > heapsort_level_bound(len(log)):
        raise RuntimeError('heapsort exceeded its worst-case level bound')
    ded = []
    for rec in log:
        if not ded or ded[-1] != rec:
            ded.append(rec)
    pairs, i = 0, 0
    while i < len(ded):
        j = i
        while j + 1 < len(ded) and ded[j + 1][0] == ded[i][0]:
            j += 1
        for u in range(i, j + 1):
            for v in range(u + 1, j + 1):
                pairs += 1
                if pairs > pair_budget:
                    return None, levels, pairs, len(log) - len(ded)
                (_, A1, B1), (_, A2, B2) = ded[u], ded[v]
                if B1 != B2:
                    x = (A1 - A2) * pow(B2 - B1, -1, rp) % rp
                    if verify(x):
                        return x, levels, pairs, len(log) - len(ded)
        i = j + 1
    return None, levels, pairs, len(log) - len(ded)


def run_toy(T, q, rp, rng, walks, Rm, branch_bits, sq_cut, dp_bits, cap_factor=40, pair_budget=2 ** 16):
    qk = [pow(q, k, rp) for k in range(12)]
    hc = (q ** 12 - 1) // rp
    g = T.rand_gt(rng, rp)
    xs = 1 + sample_below(rng, rp - 1)
    z0 = T.pow(g, xs)
    ah = [sample_below(rng, rp) for _ in range(Rm)]
    bh = [sample_below(rng, rp) for _ in range(Rm)]
    Mh = [T.mul(T.pow(g, ah[i]), T.pow(z0, bh[i])) for i in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    At = [[ah[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    Bt = [[bh[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
    mask = (1 << dp_bits) - 1
    window = walks * 20 << dp_bits

    def new_walk():
        a, b = sample_below(rng, rp), sample_below(rng, rp)
        return {"y": T.mul(T.pow(g, a), T.pow(z0, b)), "a": a, "b": b, "hist": [], "ck": None, "n": 0, "since": 0}

    def advance(y, a, b):
        kap, can, K = classinfo(T, y)
        u = toy_hash(can)
        if (u & ((1 << branch_bits) - 1)) < sq_cut:
            return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp, kap, can, K, u
        h = (u >> branch_bits) % Rm
        idx = (12 - kap) % 12
        return T.mul(y, Tab[h][idx]), (a + At[h][idx]) % rp, (b + Bt[h][idx]) % rp, kap, can, K, u

    W = [new_walk() for _ in range(walks)]
    TOY_K, TOY_BUDGET, TOY_AMAX = 1, 2 ** 40, 2 ** 25
    TOY_L = TOY_BUDGET - TOY_AMAX
    MUL_I, SQ_I, ESC_I, REP_I, RST_I, SW_I = (MUL_CORE + COMMON_CHARGED + MULT_EXTRA,
                                               SQ_CORE + COMMON_CHARGED,
                                               3 * 20000, 4000, int(2 ** 10.1 * 20000) + 43, 600)
    meter = [0, 0]

    def charge(x):
        if meter[0] >= TOY_L:
            return False
        assert x <= TOY_AMAX
        meter[0] += x
        meter[1] = max(meter[1], meter[0])
        if meter[0] > TOY_BUDGET:
            raise RuntimeError('toy meter exceeded its budget')
        return True
    visited, log = {}, []
    steps, Tc, found, check_at = 0, None, None, None
    stats = {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "duplicates_removed": 0, "sq": 0, "ord": 0,
             "records": 0, "heapsort_levels": 0, "heapsort_level_bound": 0, "pairs": 0, "matcher_runs": 0}
    stats["meter_stop"] = 0
    while found is None:
        for w in W:
            if not charge(SW_I):
                stats["meter_stop"] = 1
                break
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
                    check_at = Tc + window
                    stats["ord_at_Tc"], stats["sq_at_Tc"] = stats["ord"], stats["sq"]
            if not charge(SQ_I if (u & ((1 << branch_bits) - 1)) < sq_cut else MUL_I):
                stats["meter_stop"] = 1
                break
            if (u >> (branch_bits + 10)) & mask == 0:
                if not charge(REP_I):
                    stats["meter_stop"] = 1
                    break
                log.append((u >> (branch_bits + 10 + dp_bits), ca, cb))
                w["since"] = 0
            hist = w["hist"]
            key = (K, can)
            if len(hist) >= 2 and hist[-2][0][1] == can:
                stats["two_cycles"] += 1
                if not charge(ESC_I):
                    stats["meter_stop"] = 1
                    break
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
                if not charge(ESC_I + MUL_I * (w["n"] - w["ck"][1])):
                    stats["meter_stop"] = 1
                    break
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
            stats["sq"] += (u & ((1 << branch_bits) - 1)) < sq_cut
            w["y"], w["a"], w["b"], _, _, _, _ = advance(y, a, b)
            w["n"] += 1
            w["since"] += 1
            steps += 1
            if w["since"] > cap_factor << dp_bits:
                if not charge(RST_I):
                    stats["meter_stop"] = 1
                    break
                stats["restarts"] += 1
                w.update(new_walk())
        if stats["meter_stop"]:
            break
        if check_at is not None and steps >= check_at:
            stats["matcher_runs"] += 1
            snap = list(log)
            x, lv, pr, dup = match_log(snap, lambda v: T.pow(g, v) == z0, rp, pair_budget)
            stats["records"], stats["heapsort_levels"] = len(snap), lv
            stats["heapsort_level_bound"] = heapsort_level_bound(len(snap))
            stats["pairs"], stats["duplicates_removed"] = pr, dup
            if pr > pair_budget:
                raise RuntimeError('toy pair budget exceeded')
            if x is not None:
                if x != xs:
                    raise RuntimeError('recovered logarithm differs from the planted one')
                found = x
            else:
                check_at = steps + window
    stats["meter_max"], stats["meter_budget"] = meter[1], TOY_BUDGET
    return Tc, steps, stats


TOY_STEP_PRICES = (0.8472, 0.3708)


def part_B(rng, q, trials, Rm, branch_bits, sq_cut, dp_bits, walks):
    T = Tower(q)
    rp = toy_setup(q)
    if not (pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))):
        raise RuntimeError('toy field does not have embedding degree 12')
    N = (rp - 1) // 12
    s = sq_cut / (1 << branch_bits)
    sumq2 = s * s + (1 - s) ** 2 / Rm
    Np = N / (1 - sumq2)
    Tcs, tot, ev, dev = [], [], {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "duplicates_removed": 0,
                                 "records": 0, "heapsort_levels": 0, "heapsort_level_bound": 0, "matcher_runs": 0}, []
    max_pairs, max_level_ratio = 0, 0.0
    meter_stops, meter_ratio = 0, 0.0
    pool_ord, pool_sq = 0, 0
    while len(Tcs) < trials:
        Tc, steps, st = run_toy(T, q, rp, rng, walks, Rm, branch_bits, sq_cut, dp_bits)
        meter_ratio = max(meter_ratio, st["meter_max"] / st["meter_budget"])
        if st["meter_stop"]:
            meter_stops += 1
            continue
        Tcs.append(Tc)
        tot.append(steps)
        for k in ev:
            ev[k] += st[k]
        max_pairs = max(max_pairs, st["pairs"])
        max_level_ratio = max(max_level_ratio, st["heapsort_levels"] / max(1, st["heapsort_level_bound"]))
        pool_ord += st["ord"]
        pool_sq += st["sq"]
        sf = s
        MUL_STEP, SQ_STEP = TOY_STEP_PRICES
        nominal = (1 - sf) * MUL_STEP + sf * SQ_STEP
        o_, q_ = st["ord_at_Tc"], st["sq_at_Tc"]
        charge = ((o_ - q_) * MUL_STEP + q_ * SQ_STEP) / max(1, o_)
        dev.append((charge / nominal - 1, o_, nominal))
    n = len(Tcs)
    sN = math.sqrt(N)
    pred = math.sqrt(math.pi * Np / 2)
    s40 = math.sqrt(2 * Np * math.log(1 / 0.60))
    s402 = math.sqrt(2 * Np * math.log(1 / 0.598))
    s405 = math.sqrt(2 * Np * math.log(1 / 0.595))
    s42 = math.sqrt(2 * Np * math.log(1 / 0.58))
    s43 = math.sqrt(2 * Np * math.log(1 / 0.57))
    allsteps = max(1, sum(tot))
    return {"q": q, "r_prime": rp, "classes_N": N, "R": Rm, "branch_bits": branch_bits,
            "sq_cut": sq_cut, "squaring_prob": s,
            "sum_q2": round(sumq2, 6), "effective_N_over_N": round(1 / (1 - sumq2), 4),
            "trials_completed": n, "logs_recovered_and_verified": n,
            "mean_Tc_over_sqrtN": round(sum(Tcs) / n / sN, 4) if n else None,
            "model_mean_over_sqrtN": round(pred / sN, 4),
            "fraction_Tc_le_S_model_0.40": round(sum(t <= s40 for t in Tcs) / max(1, n), 3),
            "fraction_Tc_le_S_model_0.402": round(sum(t <= s402 for t in Tcs) / max(1, n), 3),
            "fraction_Tc_le_S_model_0.405": round(sum(t <= s405 for t in Tcs) / max(1, n), 3),
            "successes_within_S_model_0.402": sum(t <= s402 for t in Tcs),
            "binomial_standard_error_of_fraction_0.402": round(math.sqrt((sum(t <= s402 for t in Tcs) / max(1, n)) * (1 - sum(t <= s402 for t in Tcs) / max(1, n)) / max(1, n)), 4),
            "fraction_Tc_le_S_model_0.42": round(sum(t <= s42 for t in Tcs) / max(1, n), 3),
            "fraction_Tc_le_S_model_0.43": round(sum(t <= s43 for t in Tcs) / max(1, n), 3),
            "successes_within_S_model_0.43": sum(t <= s43 for t in Tcs),
            "binomial_standard_error_of_fraction_0.43": round(math.sqrt((sum(t <= s43 for t in Tcs) / max(1, n)) * (1 - sum(t <= s43 for t in Tcs) / max(1, n)) / max(1, n)), 4),
            "S_model_0.40_over_sqrtN": round(s40 / sN, 4), "S_model_0.402_over_sqrtN": round(s402 / sN, 4), "S_model_0.42_over_sqrtN": round(s42 / sN, 4),
            "mean_total_steps_incl_DP_delay_over_sqrtN": round(sum(tot) / max(1, n) / sN, 4),
            "charge_fluctuation_rms_of_(dev*sqrt(steps))_observed": round(math.sqrt(sum(d * d * k for d, k, _ in dev) / max(1, len(dev))), 4),
            "charge_fluctuation_iid_prediction": round(math.sqrt(s * (1 - s)) * (TOY_STEP_PRICES[0] - TOY_STEP_PRICES[1]) / dev[0][2], 4) if dev else None,
            "branch_rule_steps_pooled_over_trials": pool_ord,
            "squaring_fraction_pooled_over_trials": round(pool_sq / max(1, pool_ord), 4),
            "events": ev, "two_cycles_per_step": ev["two_cycles"] / allsteps,
            "matcher": "log, whole-record heapsort, dedup, equal-F pairs under a 2^16 budget, halt at first verified",
            "max_pairs_charged_in_a_trial": max_pairs,
            "schedule": "one processor, walks round-robin, chunk of 1 step iteration per turn, charged switch, one pre-debit meter",
            "meter_budget_instr_log2": 40, "max_meter_over_budget": round(meter_ratio, 8), "meter_stops": meter_stops,
            "max_heapsort_levels_over_worst_case_bound": round(max_level_ratio, 4),
            "Tc_over_sqrtN": [round(t / sN, 4) for t in Tcs] if cfg.get("dump") else None}


def part_E(rng, q, m, Rm, branch_bits, sq_cut):
    """Sample m distinct classes, apply one step F to each, and count image coincidences.
    H1 predicts N * P[F(c1) = F(c2)] = 1 - s^2 - (1-s)^2/(12R) for distinct classes c1, c2."""
    T = Tower(q)
    rp = toy_setup(q)
    N = (rp - 1) // 12
    g = T.rand_gt(rng, rp)
    Mh = [T.pow(g, 1 + sample_below(rng, rp - 1)) for _ in range(Rm)]
    Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
    gp = [g]
    for _ in range(rp.bit_length()):
        gp.append(T.mul(gp[-1], gp[-1]))
    classes = {}
    while len(classes) < m:
        e = 1 + sample_below(rng, rp - 1)
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
        if (u & ((1 << branch_bits) - 1)) < sq_cut:
            y2 = T.cyc_sqr(y)
        else:
            y2 = T.mul(y, Tab[(u >> branch_bits) % Rm][(12 - kap) % 12])
        _, can2, _ = classinfo(T, y2)
        images[can2] = images.get(can2, 0) + 1
    pairs = sum(n * (n - 1) // 2 for n in images.values())
    tot = m * (m - 1) // 2
    s = sq_cut / (1 << branch_bits)
    obs = pairs / tot * N
    return {"q": q, "classes_N": N, "sampled_classes": m, "R": Rm, "branch_bits": branch_bits,
            "sq_cut": sq_cut, "squaring_prob": s,
            "image_coincidence_pairs": pairs,
            "N_times_pair_collision_rate_observed": round(obs, 4),
            "relative_standard_error": round(1 / math.sqrt(max(1, pairs)), 4),
            "predicted_exact_1_minus_s2_minus_(1-s)^2/(12R)": round(1 - s * s - (1 - s) ** 2 / (12 * Rm), 4),
            "predicted_conservative_1_minus_s2_minus_(1-s)^2/R": round(1 - s * s - (1 - s) ** 2 / Rm, 4),
            "random_mapping": 1.0}


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
    """Read results.tsv and the six Callgrind records; assert the three blst marginals (383, 22,996,
    9,803) and each raw total; return them (the kernels' counts are charged instead)."""
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


def ledger(fc):
    """Exact recomputation of the work ledger of proof.md section 5 (instruction counts / C, C = 20,000);
    fc is part F (reduction.py): the counted pairing and the G2 output."""
    from decimal import Decimal, getcontext, localcontext
    getcontext().prec = 60
    ln2 = Decimal(2).ln()
    lg = lambda v: float(Decimal(v).ln() / ln2)
    meas = load_measure()
    C = 20000
    FPMUL = meas["ir_per_call_marginal"]["fp-mul"]
    FPADD = 45
    SQ_BLST = meas["ir_per_call_marginal"]["cyc-sqr"]
    SQ_MEASURED = SQ_CORE
    FP2SQR = 2 * FPMUL + 3 * FPADD
    FP2ADD = 2 * FPADD
    GS_ANALYTIC = 9 * FP2SQR + 34 * FP2ADD + 200
    common_items = {"read Re(c_0), branch bits, DP test, table address": 50,
                    "(alpha, beta) update mod r (add or doubling)": 100,
                    "2-cycle and checkpoint comparisons on T, history": 60,
                    "loop control": 11,
                    "margin for loads and moves": 28}
    meter_items = {"work meter: 3-word (192-bit) accumulate, step charge in a register (1) and adds/adcs/adc (3), "
                   "then the exact 3-word compare with L in three registers, cmp/sbcs/sbcs/b.hs, branch iff meter >= L (4), "
                   "executed on every step": 8,
                   "multiword counters: 2-word stuck-walk counter increment, adds/adc (2), and its "
                   "high-word threshold compare against 5120, exact as 20/theta = 5120*2^64; 5120 is not an "
                   "AArch64 compare immediate (imm12, optionally LSL #12), so it is materialised in a register "
                   "on every step, mov x16,#5120; cmp x9,x16; b.hs (3)": 5}
    mult_items = {"4 Fp multiplications (a*zeta, a*d0, b*d1, B*zeta)": 4 * FPMUL,
                  "1 Fp addition (B) + 8 Fp sub/neg": 9 * FPADD,
                  "minimum of 12 keys (top limbs, tie test)": 150,
                  "load c_1, kappa, twist index": 32}
    loop_items = {"back-branch to the step head (b)": 1,
                  "chunk counter, subs x12,x12,#1; b.eq switch": 2,
                  "branch select on the 6-bit field of T (and #63, cmp #36, b.lo; 2 of the old margin)": 3,
                  "distinguished-point test branch (b.eq)": 1,
                  "2-cycle and checkpoint test branches (b.ne x4)": 4}
    assert sum(loop_items.values()) == common_items["loop control"] == 11
    common, extra = sum(common_items.values()), sum(mult_items.values())
    common_charged = common + sum(meter_items.values())
    assert (common, extra, GS_ANALYTIC, common_charged) == (249, 2119, 11369, 262)
    assert (MUL_CORE, SQ_MEASURED, common_charged, extra) == (14563, 7154, COMMON_CHARGED, MULT_EXTRA)
    PARENT = (11400, 300, 2200)
    RECORD = (GS_ANALYTIC, 270, extra)
    OURS = (SQ_MEASURED, common_charged, extra)
    N = Fraction(r - 1, 12)
    Rt = 1024
    Pm = Decimal("0.402")
    eps = Decimal(5) / Decimal(12 * Rt) + Decimal(250) / Decimal(Rt * Rt)

    Pm_OLD = Decimal("0.43")

    def walk_length(s_, Pm_=Pm):
        """N' = N/(1 - s^2 - (1-s)^2/R), S = sqrt(2 N' ln(1/(1-Pm_))) for squaring fraction s_."""
        sumq2_ = s_ * s_ + (1 - s_) ** 2 / Rt
        Np_ = N / (1 - sumq2_)
        Npd_ = Decimal(Np_.numerator) / Decimal(Np_.denominator)
        return sumq2_, Npd_, (2 * Npd_ * (1 / (1 - Pm_)).ln()).sqrt()

    S_OLD = Fraction(3, 8)
    s = Fraction(36, 64)
    assert s == Fraction(SQ_CUT_PROD, 1 << BRANCH_BITS_PROD)
    sumq2, Npd, S = walk_length(s)
    THETA_BITS = 72
    E_OLD = Decimal(2 ** 20) * 20 * Decimal(2 ** 60)
    WALKS, RESTART_CAP = 2 ** 20, 2 ** 36
    STUCK = 20 * Decimal(2) ** THETA_BITS
    assert STUCK == 5120 * Decimal(2) ** 64
    E = (WALKS + RESTART_CAP) * STUCK
    N_MAX = Decimal(2) ** Decimal("125.162") * C / (SQ_MEASURED + common_charged)
    assert Decimal("126.58") < lg(N_MAX) < Decimal("126.61")
    GAP_STEPS = N_MAX
    EXP_GAPS = GAP_STEPS * Decimal(2) ** (-THETA_BITS)
    EXP_STUCK = EXP_GAPS * (-Decimal(20)).exp()
    EXP_LONG = GAP_STEPS * Decimal(2) ** (-100)
    EXP_RESTARTS = EXP_STUCK + EXP_LONG
    assert RESTART_CAP == 2 ** 36
    assert Decimal("25.73") < lg(EXP_STUCK) < Decimal("25.76")
    assert Decimal("26.58") < lg(EXP_LONG) < Decimal("26.61")
    assert Decimal("27.22") < lg(EXP_RESTARTS) < Decimal("27.25")
    assert lg(EXP_RESTARTS) - 36 < Decimal("-8.65") and lg(EXP_RESTARTS) < 36 - 8
    RESTART_STEPS = Decimal(RESTART_CAP) * STUCK
    assert lg(RESTART_STEPS) < Decimal("112.33")
    RESTART_FRAC = RESTART_STEPS / N_MAX
    assert RESTART_FRAC < Decimal("5.07e-5")
    xq = Decimal(1) - Decimal(s.numerator) / Decimal(s.denominator)
    assert xq == Decimal(28) / Decimal(64)
    tail = xq ** 148 * (148 - 147 * xq) / (1 - xq) ** 2
    EXP_DPCYC = N_MAX * Decimal(2) ** (-THETA_BITS) * (Decimal(2) ** (-288) * sum(range(65, 148)) + tail)
    DPCYC_LOSS = Decimal("2e-11")
    assert EXP_DPCYC < DPCYC_LOSS
    REPORT_CAP = 2 ** 55
    SLOT_BYTES = 80
    ZF_STORES, ZF_RECS = 2 ** 21 // 16, 2 ** 21 // 80
    ZF_INSTR = -(-(2 * ZF_STORES + 2 ** 16) // ZF_RECS)
    assert (ZF_STORES, ZF_RECS, ZF_INSTR) == (131072, 26214, 13)
    ZF_ITEM = "operating-system zero-fill of the log's 2 MiB pages on first write (131,072 stores and one fault per 26,214 records)"
    report_items = {"12 keys for kappa (the multiply-step key list)": extra,
                    "two 255-bit modular multiplications by p^kappa mod r": 1600,
                    "pre-debit: compare of the meter with the threshold (4), then the charge in a register and its 3-word add (4)": 8,
                    "report counter (one 64-bit word) compare against the cap, branch, increment": 5,
                    "slot address base + 80 c (two shifted adds)": 2,
                    "append the 80-byte record (ten 8-byte stores as five stp, address set-up)": 12,
                    "reset of the 2-word stuck-walk counter at the report (two moves into the register pair)": 4,
                    ZF_ITEM: ZF_INSTR,
                    "margin for loads, moves, call and return": 250 - ZF_INSTR}
    report_instr = sum(report_items.values())
    assert report_instr == 4000
    assert r.bit_length() == 255 and r < 2 ** 255
    SAMPLE_TABLE_DRAWS = 2 * 1024
    table_units = Decimal(1024 * (1021 + 11)) + Decimal(2 * 12288 * 800) / C + Decimal(SAMPLE_TABLE_DRAWS * SAMPLE_MAX_INSTR) / C + Decimal(12288 * 54 * 45) / C
    assert table_units < Decimal(2) ** Decimal("20.2")
    ONE_START = Decimal(2) ** Decimal("10.1")
    START_MAX = Decimal(1021) + Decimal(2 * SAMPLE_MAX_INSTR) / C + Decimal(43) / C
    assert START_MAX < ONE_START and ONE_START - START_MAX > 16
    SETUP_UNITS = Decimal(2 ** 31)
    SETUP_ITEMS = Decimal(2) ** Decimal("14.1") + 510 + Decimal(2) ** Decimal("20.2") + WALKS * ONE_START + Decimal(32) / C
    assert SETUP_ITEMS <= SETUP_UNITS
    PREPROC = SETUP_UNITS
    restart_count_units = Decimal(RESTART_CAP) * 3 / C
    assert lg(restart_count_units) < Decimal("23.30")
    METER_UPDATE = 40
    CHUNK_K = 2 ** 63
    assert STUCK % CHUNK_K == 0 and STUCK / CHUNK_K == 5 * 2 ** 11
    switch_items = {"store the 2 KB walk state (256 stores)": 256, "load the next walk's 2 KB state (256 loads)": 256,
                    "next walk index mod M and its state address": 4, "chunk counter set to K (movz)": 1,
                    "the switch's pre-debit compare and charge add": 8, "branch, call, return and margin": 75}
    SWITCH_INSTR = sum(switch_items.values())
    assert SWITCH_INSTR == 600
    SWITCHES = N_MAX / CHUNK_K + WALKS + RESTART_CAP + 1
    SWITCH_UNITS = SWITCHES * SWITCH_INSTR / C
    assert Decimal("58.52") < lg(SWITCH_UNITS) < Decimal("58.55")
    dp_units = Decimal(REPORT_CAP) * report_instr / C
    start_units = Decimal(WALKS) * METER_UPDATE / C
    restart_units = Decimal(RESTART_CAP) * METER_UPDATE / C
    restart_start_units = Decimal(RESTART_CAP) * ONE_START
    WALK_NONSTEP = dp_units + start_units + restart_units + restart_start_units + restart_count_units
    level_items = {"child compare (20 ldr, cmp and 9 sbcs, 11 csel; flags_check.py)": 41,
                   "parent compare (10 ldr, cmp and 9 sbcs, b.hs)": 21,
                   "move one 80-byte record on a swap (10 loads, 10 stores)": 20,
                   "index arithmetic and loop control": 50}
    LEVEL_INSTR = sum(level_items.values())
    assert LEVEL_INSTR == 132
    OUTER_INSTR = 40
    nrec = REPORT_CAP
    LEVELS = heapsort_level_bound(nrec)
    assert LEVELS == 54 * 2 ** 55 + 57
    assert sum(k.bit_length() - 1 for k in range(2, 4097)) + 4096 == heapsort_level_bound(4096)
    OUTER = nrec // 2 + nrec - 1
    SORT_INSTR = LEVEL_INSTR * LEVELS + OUTER_INSTR * OUTER
    DEDUP_INSTR = 100 * nrec
    PAIR_BUDGET, PAIR_UNITS = 2 ** 16, 600
    HALT_INSTR = 100
    PAIRING_UNITS = 2 ** 14
    assert max(fc["pairing_units_counted_g"], fc["pairing_units_counted_verification"]) <= PAIRING_UNITS
    G2_OUT = fc["g2_output_units_hard_max_255_bit_scalar"]
    OUTPUT_UNITS = G2_OUT + PAIRING_UNITS
    R_STOP_INSTR = SORT_INSTR + DEDUP_INSTR + HALT_INSTR + (PAIR_BUDGET * PAIR_UNITS + OUTPUT_UNITS) * C
    R_STOP = Decimal(R_STOP_INSTR) / C
    sort_units = Decimal(SORT_INSTR) / C
    assert lg(sort_units) < Decimal("53.53") and lg(R_STOP) < Decimal("53.55")
    SETUP_OLD = Decimal(2) ** Decimal("63.8")
    SETUP = Decimal(2) ** Decimal("54.2")
    NONWALK = WALK_NONSTEP + SETUP_UNITS + R_STOP
    assert NONWALK < SETUP and Decimal("54.1") < lg(NONWALK) < Decimal("54.181") and lg(SETUP_ITEMS) < Decimal("30.11")
    assert 2 * 20 / (Decimal(20).exp() - 1) < Decimal("8.25e-8")
    SAMPLER_LOSS = Decimal(2) ** Decimal("-218.99")
    TIE_LOSS = Decimal(2) ** -247
    assert N_MAX * 66 / p < TIE_LOSS
    SMALL_LOSSES = RESTART_FRAC + DPCYC_LOSS + Decimal("8.25e-8") + Decimal("1e-6") + SAMPLER_LOSS + TIE_LOSS
    LOSS_SUM = Decimal("0.01") + SMALL_LOSSES
    assert LOSS_SUM < Decimal("0.0101")
    assert Pm - Decimal("0.0101") >= Decimal("0.3919") >= Decimal("0.39")
    store = {"shape": "append-only log of 80-byte records addressed by the report counter; matched once after "
                      "the walk phase by an in-place heapsort of the whole record, adjacent dedup and a hard pair budget",
             "slots": REPORT_CAP, "record_bytes": SLOT_BYTES, "log_bytes_log2": lg(Decimal(REPORT_CAP * SLOT_BYTES)),
             "heapsort_levels_worst_case": "54 * 2^55 + 57", "heapsort_levels_per_record": float(Decimal(LEVELS) / nrec),
             "level_items": level_items, "level_instr": LEVEL_INSTR, "outer_iteration_instr": OUTER_INSTR,
             "sort_instr_per_record": float(Decimal(SORT_INSTR) / nrec), "sort_units_log2": lg(sort_units),
             "dedup_instr_per_record": 100, "pair_budget_log2": 16, "pair_units": PAIR_UNITS,
             "pair_units_log2": lg(Decimal(PAIR_BUDGET * PAIR_UNITS)),
             "stopping_allowance_units_log2": lg(R_STOP),
             "success_loss_sum_h1plus_small_terms": float(LOSS_SUM),
             "success_after_every_loss": float(Pm - LOSS_SUM),
             "margin_over_0.39_after_every_loss": float(Pm - LOSS_SUM - Decimal("0.39"))}

    max_offset = REPORT_CAP * SLOT_BYTES - 1
    assert max_offset < 2 ** 63 and 80 == (1 << 6) + (1 << 4)
    assert 2 * (REPORT_CAP - 1) + 2 < 2 ** 63
    assert int(REPORT_CAP) < 2 ** 63
    T_BITS = p.bit_length()
    assert BRANCH_BITS_PROD + 10 + THETA_BITS + 128 <= T_BITS
    F_SHIFT = BRANCH_BITS_PROD + 10 + THETA_BITS
    assert F_SHIFT == 88
    words = lambda v: -(-int(v).bit_length() // 64)
    assert words(STUCK) == 2 and words(Decimal(2) ** 140) == 3
    widths = {"log_bytes_log2": lg(Decimal(REPORT_CAP * SLOT_BYTES)), "max_byte_offset_bits": max_offset.bit_length(),
              "report_cap_bits": int(REPORT_CAP).bit_length(), "report_bytes_at_cap_bits": (int(REPORT_CAP) * SLOT_BYTES).bit_length(),
              "dp_bits_of_T": THETA_BITS, "dp_bit_range_of_T": [BRANCH_BITS_PROD + 10, BRANCH_BITS_PROD + 10 + THETA_BITS - 1],
              "fingerprint_bit_range_of_T": [F_SHIFT, F_SHIFT + 127], "T_bits": T_BITS,
              "stuck_walk_threshold_log2": lg(STUCK), "stuck_walk_counter_words": words(STUCK)}
    FALSE_CANDIDATES = Decimal(REPORT_CAP) * REPORT_CAP / Decimal(2) ** 129
    assert FALSE_CANDIDATES == Decimal(2) ** -19
    assert lg(FALSE_CANDIDATES / PAIR_BUDGET) == -35
    assert 25.2 < math.log2(2 ** 16 * 600) < 25.3

    def to_dec(fr):
        return Decimal(fr.numerator) / Decimal(fr.denominator)

    def walk_ledger(SQ, COM, MULT, s_=s, scale=Decimal(1), window=E, setup=None, mul_core=C, Pm_=Pm):
        """cbar and W for one charge vector; `scale` prices drift, `window` the outside bucket."""
        cbar = (1 - s_) * Fraction(mul_core + COM + MULT, C) + s_ * Fraction(SQ + COM, C)
        cb = to_dec(cbar) * scale
        return cbar, cb, cb * walk_length(s_, Pm_)[2] * (1 + eps) + cb * window + (SETUP + SWITCH_UNITS if setup is None else setup)

    cbar_p, cb_p, W_p = walk_ledger(*PARENT, s_=S_OLD, window=E_OLD, setup=SETUP_OLD, Pm_=Pm_OLD)
    cbar_r, cb_r, W_r = walk_ledger(*RECORD, s_=S_OLD, window=E_OLD, setup=SETUP_OLD, Pm_=Pm_OLD)
    cbar, cb, W = walk_ledger(*OURS, mul_core=MUL_CORE)
    if not (abs(lg(W_p) - 125.71472902377629) < 1e-12 and cbar_p == Fraction(9225, 10000)):
        raise RuntimeError('earliest ledger does not reproduce')
    if not (abs(lg(W_r) - 125.7074974545452) < 1e-12 and cbar_r == Fraction(9178875, 10000000)):
        raise RuntimeError('intermediate ledger does not reproduce')
    pass
    mul_step = Fraction(MUL_CORE + common_charged + extra, C)
    sq_step = Fraction(OURS[0] + common_charged, C)
    assert (mul_step, sq_step) == (Fraction(16944, 20000), Fraction(7416, 20000))
    assert (float(mul_step), float(sq_step)) == TOY_STEP_PRICES
    assert (MUL_CORE + common_charged + extra, OURS[0] + common_charged) == (16944, 7416)
    gap = mul_step - sq_step
    DROP_POINTS = Decimal("1.5")

    def delta(D):
        return Decimal(str(D)) / 100 * to_dec(gap) / to_dec(cbar)

    _, _, W_drop = walk_ledger(*OURS, scale=1 + delta(DROP_POINTS), mul_core=MUL_CORE)
    claim_thousandths = int((W_drop.ln() / ln2 * 1000).to_integral_value(rounding="ROUND_CEILING"))
    W_up_digits = Decimal(W_drop).quantize(Decimal(1).scaleb(int(W_drop.adjusted()) - 40), rounding="ROUND_CEILING")
    fnum, fden = Fraction(W_up_digits).numerator, Fraction(W_up_digits).denominator
    assert fnum ** 1000 <= fden ** 1000 * 2 ** claim_thousandths, "claim is below the drift-covering cap"
    assert fnum ** 1000 > fden ** 1000 * 2 ** (claim_thousandths - 1) , "claim is not the least thousandth"
    CLAIM = Decimal(claim_thousandths) / 1000
    assert CLAIM <= Decimal("125.162")
    assert CLAIM == Decimal("125.162"), CLAIM
    Wcap = Decimal(2) ** CLAIM
    widths["work_meter_cap_instr_log2"] = lg(Wcap * C)
    widths["work_meter_words"] = words(Wcap * C)
    assert lg(Wcap * C) < 139.55 and words(Wcap * C) == 3
    with localcontext() as _ctx:
        _ctx.prec = 90
        WCAP_INSTR = int((Decimal(2) ** CLAIM * C).to_integral_value(rounding="ROUND_FLOOR"))
    assert WCAP_INSTR ** 1000 <= 2 ** claim_thousandths * C ** 1000 < (WCAP_INSTR + 1) ** 1000
    A_MAX = 2 ** 25
    assert ONE_START * C + METER_UPDATE + 3 <= A_MAX
    assert SWITCH_INSTR <= A_MAX
    assert C + common_charged + extra <= A_MAX and report_instr <= A_MAX and 3 * C <= A_MAX
    assert C + common_charged + extra + report_instr + 3 * C <= A_MAX
    B_WALK = WCAP_INSTR - int(SETUP_UNITS) * C - R_STOP_INSTR
    L_THRESH = B_WALK - A_MAX
    assert int(SETUP_UNITS) * C + B_WALK + R_STOP_INSTR <= WCAP_INSTR
    assert 0 < L_THRESH < B_WALK < 2 ** 192
    walk_need_drop = cb * (1 + delta(DROP_POINTS)) * (S * (1 + eps) + E) + WALK_NONSTEP + SWITCH_UNITS
    assert walk_need_drop * C <= L_THRESH
    budget_units = Decimal(L_THRESH) / C
    n = S * (1 + eps)
    cap_over_ledger = Wcap / W - 1
    slack = (budget_units - WALK_NONSTEP - SWITCH_UNITS - cb * E) / n - cb
    tol_sq_fraction = slack * n / (to_dec(gap) * (n + E))
    assert tol_sq_fraction >= DROP_POINTS / 100
    assert tol_sq_fraction >= Decimal("0.015")
    assert (cb + to_dec(gap) * Decimal("0.015")) * (n + E) + WALK_NONSTEP + SWITCH_UNITS <= budget_units
    EXP_REPORTS = Decimal(B_WALK) / C / cb * Decimal(2) ** (-THETA_BITS)
    assert Decimal("53.94") < lg(EXP_REPORTS) < Decimal("53.97") and REPORT_CAP / EXP_REPORTS > Decimal("1.94")

    def success_with_drift(drop_points):
        """Success probability if the squaring fraction drops by `drop_points` percentage points: the
        average step price rises by the factor 1 + delta, so the cap buys n/(1+delta) steps when it binds,
        and the collision law 1 - exp(-steps^2/(2N')) gives the success value."""
        d = delta(drop_points)
        steps = min(n, (budget_units - WALK_NONSTEP - SWITCH_UNITS - cb * E * (1 + d)) / (cb * (1 + d)))
        return float(1 - (-(steps * steps) / (2 * Npd)).exp())

    drift_table = {p_: success_with_drift(p_) for p_ in ("0.5", "1", "1.5", "2", "2.5", "3", "4", "5")}

    def spent(D):
        return max(0.0, float(Pm) - success_with_drift(D))

    lo_, hi_ = Decimal(2), Decimal(20)
    for _ in range(60):
        mid = (lo_ + hi_) / 2
        if spent(mid) < 0.01:
            lo_ = mid
        else:
            hi_ = mid
    return {"measurement_record": meas,
            "fp2_sqr_instr_analytic": FP2SQR, "cyclotomic_square_instr_analytic_listing_not_charged": GS_ANALYTIC,
            "multiply_kernel_instr_charged": MUL_CORE, "cyclotomic_square_instr_charged": SQ_MEASURED,
            "cyclotomic_square_instr_blst_measured_baseline": SQ_BLST,
            "common_items": common_items, "loop_control_items": loop_items, "meter_and_counter_items": meter_items,
            "multiply_only_items": mult_items, "common_counted": common,
            "meter_and_counter_counted": common_charged - common, "common_charged": common_charged,
            "multiply_extra_counted": extra,
            "multiply_step_units": float(mul_step), "squaring_step_units": float(sq_step),
            "multiply_step_instr": MUL_CORE + common_charged + extra, "squaring_step_instr": SQ_MEASURED + common_charged,
            "step_price_gap_units": float(gap),
            "report_store": store, "word_widths": widths, "theta_bits": THETA_BITS,
            "expected_false_fingerprint_candidates_log2": lg(FALSE_CANDIDATES), "report_items": report_items, "report_instr": report_instr,
            "report_cap": float(REPORT_CAP), "report_cap_log2": lg(REPORT_CAP),
            "a_max_instr_log2": 25, "walk_budget_instr_log2": math.log2(B_WALK), "threshold_instr_log2": math.log2(L_THRESH),
            "setup_plus_walk_budget_plus_stopping_allowance_le_cap": True,
            "walk_need_at_design_drift_units_log2": lg(walk_need_drop),
            "expected_reports_log2": lg(EXP_REPORTS), "report_cap_over_expected_reports": float(REPORT_CAP / EXP_REPORTS),
            "expected_reports_at_theta_2^-71_over_report_cap": round(float(2 * EXP_REPORTS / REPORT_CAP), 4),
            "hard_max_walk_steps_log2": lg(N_MAX),
            "meter_bookkeeping_instr_per_start_or_restart": METER_UPDATE,
            "meter_units_starts_log2": lg(start_units), "meter_units_restarts_log2": lg(restart_units),
            "restart_cap_log2": RESTART_CAP.bit_length() - 1,
            "expected_gaps_log2": lg(EXP_GAPS), "expected_stuck_restarts_log2": lg(EXP_STUCK),
            "expected_long_cycle_restarts_log2": lg(EXP_LONG), "expected_restarts_total_log2": lg(EXP_RESTARTS),
            "expected_restarts_over_cap_log2": lg(EXP_RESTARTS) - 36,
            "restart_worst_case_step_fraction_of_hard_max": float(RESTART_FRAC),
            "restart_worst_case_steps_log2": lg(RESTART_STEPS),
            "expected_distinguished_long_cycles_log2": lg(EXP_DPCYC), "expected_distinguished_long_cycles": float(EXP_DPCYC),
            "distinguished_long_cycle_loss_charged": float(DPCYC_LOSS),
            "dp_units_log2": lg(dp_units), "nonwalk_bucket_log2": lg(SETUP),
            "chunk_steps_log2": 63, "switch_items": switch_items, "switch_instr": SWITCH_INSTR,
            "switches_hard_max_log2": lg(SWITCHES), "switch_units_log2": lg(SWITCH_UNITS),
            "preprocessing_units_log2": lg(PREPROC), "setup_items_units_log2": lg(SETUP_ITEMS),
            "sampler": {"r_max_attempts": SAMPLE_R_MAX, "attempt_instr": SAMPLE_ATTEMPT_INSTR,
                        "max_instr_per_sample": SAMPLE_MAX_INSTR, "max_units_per_sample": SAMPLE_MAX_INSTR / 20000,
                        "table_units_log2": lg(Decimal(SAMPLE_TABLE_DRAWS * SAMPLE_MAX_INSTR) / C),
                        "start_max_units": float(START_MAX), "one_start_bound_units": float(ONE_START),
                        "exhaustion_rule": "return 1", "exhaustion_probability_log2_at_most": -SAMPLE_R_MAX,
                        "exhaustion_allowance_log2_charged": -218.99,
                        "chacha20_block_item_instr": sum(SAMPLE_BLOCK_ITEMS.values()),
                        "chacha20_block_charged_instr": SAMPLE_BLOCK_INSTR},
            "sampler_loss_charged": float(SAMPLER_LOSS),
            "restart_counter_units_log2": lg(restart_count_units), "nonwalk_sum_log2": lg(NONWALK),
            "restart_start_units_log2": lg(restart_start_units),
            "earliest_ledger": {"charges": PARENT, "cbar": float(cbar_p), "log2_W": lg(W_p)},
            "intermediate_ledger": {"charges": RECORD, "cbar": float(cbar_r), "log2_W": lg(W_r)},
            "sum_q2": float(sumq2), "log2_N": lg(Decimal(N.numerator) / Decimal(N.denominator)),
            "log2_N_effective": lg(Npd), "log2_S_model_0.42": lg(S), "cbar_units_per_step": float(cb),
            "fruitless_cycle_eps": float(eps), "log2_detection_window_steps": lg(E), "log2_W": lg(W),
            "drift_design_points": float(DROP_POINTS), "drift_factor_at_design": float(delta(DROP_POINTS)),
            "log2_W_at_design_drift": lg(W_drop),
            "submitted_bound": str(CLAIM), "excess_of_cap_over_ledger": float(cap_over_ledger),
            "cap_minus_nominal_ledger_units_log2": lg(Wcap - W),
            "walk_budget_slack_per_step": float(slack), "tolerated_squaring_fraction_drop_points": float(tol_sq_fraction * 100),
            "success_vs_squaring_fraction_drop_points": drift_table,
            "drop_points_where_0.01_is_spent": float(lo_),
            "success_at_model_0.42": float(1 - (-(S * S) / (2 * Npd)).exp()),
            "sensitivity_bits_per_100_instr_on_squaring": lg(walk_ledger(SQ_MEASURED + 100, common_charged, extra, mul_core=MUL_CORE)[2]) - lg(W),
            "sensitivity_bits_per_100_instr_every_step": lg(walk_ledger(SQ_MEASURED, common_charged + 100, extra, mul_core=MUL_CORE)[2]) - lg(W),
            "sensitivity_bits_at_full_multiplication_price": lg(walk_ledger(C - common_charged, common_charged, extra, mul_core=MUL_CORE)[2]) - lg(W) }


out = {}
out["partA_production_BLS12_381"] = part_A(rng, PROD_CHECKS)
out["partS_bounded_sampler"] = part_S()
if cfg.get("production_steps", 3000) > 0 and PROD_CHECKS > 0:
    out["partC_production_walk"] = part_C(rng, PROD_STEPS, stored=True)
    out["partC_production_walk_integer_rule"] = part_C(rng, PROD_STEPS, stored=False)
out["partB_toy"] = part_B(rng, Q_TOY, TRIALS, R, BRANCH_BITS, SQ_CUT, DP_BITS, WALKS)
if PAIR_SAMPLES > 0:
    out["partE_pair_collision_law"] = part_E(rng, Q_TOY, PAIR_SAMPLES, R, BRANCH_BITS, SQ_CUT)
here = os.path.dirname(locate_measure())
out["partF_reduction_and_G2_output"] = runpy.run_path(os.path.join(here, "reduction.py"),
                                                      init_globals={"p": p, "r": r, "Tower": Tower})["run"](rng)
out["partD_ledger"] = ledger(out["partF_reduction_and_G2_output"])
out["nzcv_flags_check"] = runpy.run_path(os.path.join(here, "flags_check.py"))["run"]()
out["seconds"] = round(time.time() - T0, 2)
out["scope"] = ("Parts A and C check exact identities and the production branch rule (with a 16-entry table "
                "and no distinguished-point test) on real BLS12-381 GT elements; Part F checks the reduction "
                "and the G2 output on BLS12-381; Part B is a toy field; none validates production costs.")
print(json.dumps(out, indent=1))
