"""Frobenius-class Pollard rho with a cyclotomic-squaring branch read from the free class invariant T0.

T0 = Re(a_0) is the Fp part of the w^0 coefficient of an Fp12 element in the tower basis. Frobenius maps
a_0 to its conjugate, so T0 is equal on the 12 elements of a Frobenius class and costs no multiplication.
Every per-step decision (partition h, squaring branch, distinguished test) is read from disjoint bit fields
of T0; the 12 conjugate keys (four Fp multiplications) are computed only on multiplication steps and at
reports.

Part A (BLS12-381 itself): GT membership, Frobenius formula, the four-multiplication key formula and its
sign-set minimum, canonical representative, p^6 = inversion, ord_r(p) = 12; T0 equal on all 12 conjugates
(so the three bit fields are class invariants); Granger-Scott cyclotomic squaring equal to y*y, with its
Fp multiplications and Fp additions/subtractions counted by running it; the factorisation of r-1 and
ord_r(2) = (r-1)/16.
Part A2 (BLS12-381): frequencies of the low 73 bits of T0 (partition, branch and distinguished fields) and
of the 128 record-key bits of the class key K on real GT elements.
Part P (BLS12-381): the production step rule run on real GT elements with a 16-entry multiplier table:
two walkers started at two different conjugates of one element must produce the same class sequence (the
ordinary step is a map on classes), the exponent invariant y = g^alpha z^beta is checked at the end, and
the squaring fraction and the number of key computations are counted. A second test starts two walkers at
one class with DIFFERENT histories, one with an all-multiply 2-cycle already in its two-step T0 ring (it
escapes at once) and one freshly started at a different conjugate (it walks the cycle once first), and
shows that they escape at the same class and coincide from then on within one cycle length.
Part E (toy fields): a direct pair-collision test of the branch-adjusted collision law. Distinct classes
are sampled, one step is applied to each, and image coincidences are counted by pair type.
Part B (toy fields near 2^20): the complete mixed walk with the production order (2-cycle and checkpoint
checks before the distinguished test, append rule, per-slot report cap, per-slot work meter), every trial
run to its fixed step budget 1.03 S0', the report log matched offline by the production matcher.
Part C: the bound arithmetic of proof.md, including an exact integer check that the claim is log2 W
rounded up at three decimals.
"""
import json, math, random, sys, time, hashlib
from math import ceil
from fractions import Fraction
from collections import Counter, defaultdict

T_START = time.time()
cfg = json.load(sys.stdin)
trials = cfg.get("trials", 200)
walks = cfg.get("walks", 4)
seed = cfg.get("seed", 1)
prod_checks = max(1, cfg.get("production_checks", 2))
budget = cfg.get("seconds", 10)
toy_bits = cfg.get("toy_bits", "t0")
pair_samples = cfg.get("pair_samples", 2500)
fields = cfg.get("fields", "all")
long_horizon = cfg.get("long_horizon", True)
parts = cfg.get("parts", "AQPEBC")
rng = random.Random(seed)


class Tower:
    """Fq2 = Fq[u]/(u^2+1), Fq12 = Fq2[w]/(w^6 - xi), xi = 1+u; x[m] is the Fq2 coefficient of w^m.
    blst's layout is x[2j] = a[0][j], x[2j+1] = a[1][j] (Fp12 = Fp6[w]/(w^2 - v), v = w^2)."""
    def __init__(self, q):
        self.q = q
        self.xi = (1, 1)
        self.d = [self.f2pow(self.xi, m * (q - 1) // 6) for m in range(6)]
        z = self.f2pow(self.xi, (q * q - 1) // 6)
        assert z[1] == 0 and (z[0] * z[0] - z[0] + 1) % q == 0, "zeta must be a primitive 6th root in Fq"
        self.zeta = z[0]
        self.ops = Counter()
    def f2mul(self, a, b):
        q = self.q
        return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)
    def f2pow(self, a, e):
        r = (1, 0)
        while e:
            if e & 1: r = self.f2mul(r, a)
            a = self.f2mul(a, a); e >>= 1
        return r
    def mul(self, x, y):
        q = self.q
        acc = [[0, 0] for _ in range(11)]
        for i in range(6):
            a0, a1 = x[i]
            for j in range(6):
                b0, b1 = y[j]
                t = acc[i + j]
                t[0] += a0 * b0 - a1 * b1; t[1] += a0 * b1 + a1 * b0
        out = []
        for m in range(6):
            c0, c1 = acc[m]
            if m + 6 < 11:  # w^6 = 1+u
                h0, h1 = acc[m + 6]
                c0 += h0 - h1; c1 += h0 + h1
            out.append((c0 % q, c1 % q))
        return tuple(out)
    def one(self):
        return ((1, 0),) + ((0, 0),) * 5
    def pow(self, x, e):
        r = self.one()
        while e:
            if e & 1: r = self.mul(r, x)
            x = self.mul(x, x); e >>= 1
        return r
    def frob(self, x):
        return tuple(self.f2mul((c[0], (-c[1]) % self.q), self.d[m]) for m, c in enumerate(x))
    def keys_direct(self, x):
        out = []
        for _ in range(12):
            out.append(x[1][0]); x = self.frob(x)
        return out
    def keys_formula(self, x):
        """Re(coefficient of w^1) of x^(q^k), k=0..11, from 4 Fq multiplications.
        k=2j: a*zeta^j ; k=2j+1: b*zeta^j with a=Re c1, b=Re(conj(c1)*d1)."""
        q, z = self.q, self.zeta
        a = x[1][0]
        b = (x[1][0] * self.d[1][0] + x[1][1] * self.d[1][1]) % q   # 2 mults
        za = z * a % q; zb = z * b % q                               # 2 mults
        A = [a, za, (za - a) % q, (-a) % q, (-za) % q, (a - za) % q]
        B = [b, zb, (zb - b) % q, (-b) % q, (-zb) % q, (b - zb) % q]
        out = [0] * 12
        for j in range(6):
            out[2 * j] = A[j]; out[2 * j + 1] = B[j]
        return out
    def key_min_structured(self, x):
        """kappa and class key K from the sign structure of the twelve keys (two sign-symmetric
        triples), with no twelve-way scan."""
        q, z = self.q, self.zeta
        a = x[1][0]
        b = (x[1][0] * self.d[1][0] + x[1][1] * self.d[1][1]) % q
        best = None
        for base, off in ((a, 0), (b, 1)):
            mags = [base, base * z % q, (base * z - base) % q]
            jm = min(range(3), key=lambda j: mags[j])
            jM = max(range(3), key=lambda j: mags[j])
            for j, neg in ((jm, False), (jM, True)):
                k = 2 * j + off + (6 if neg else 0)
                val = (-mags[j]) % q if neg else mags[j]
                if best is None or val < best[1]:
                    best = (k, val)
        return best
    def cyc_sqr(self, x):
        """Granger-Scott squaring in the cyclotomic subgroup: three Fp4 squarings of the pairs
        (a00, a11), (a10, a02), (a01, a12), then ret = 3t - 2a on the a[0] half and 3t + 2a on the a[1]
        half (with xi*t2[1] for ret[1][0]). Same formula as cyc_sqr_counted, without the counters."""
        q = self.q
        def s2(a):
            a0, a1 = a
            return ((a0 + a1) * (a0 - a1) % q, 2 * a0 * a1 % q)
        def s4(a, b):
            t0 = s2(a); t1 = s2(b); u = s2(((a[0] + b[0]) % q, (a[1] + b[1]) % q))
            return (((t1[0] - t1[1] + t0[0]) % q, (t1[0] + t1[1] + t0[1]) % q),
                    ((u[0] - t0[0] - t1[0]) % q, (u[1] - t0[1] - t1[1]) % q))
        A = s4(x[0], x[3]); B = s4(x[1], x[4]); Cc = s4(x[2], x[5])
        xc = ((Cc[1][0] - Cc[1][1]) % q, (Cc[1][0] + Cc[1][1]) % q)
        def m3(t, a, sg):
            return ((3 * t[0] + sg * 2 * a[0]) % q, (3 * t[1] + sg * 2 * a[1]) % q)
        return (m3(A[0], x[0], -1), m3(xc, x[1], 1), m3(B[0], x[2], -1),
                m3(A[1], x[3], 1), m3(Cc[0], x[4], -1), m3(B[1], x[5], 1))
    # Counted Fp-level primitives: the squaring routine below is built only from these, so running it
    # counts its Fp multiplications and its Fp additions/subtractions exactly.
    def _m(self, a, b): self.ops["fp_mul"] += 1; return a * b % self.q
    def _a(self, a, b): self.ops["fp_add"] += 1; return (a + b) % self.q
    def _s(self, a, b): self.ops["fp_sub"] += 1; return (a - b) % self.q
    def _sqr2(self, x):
        a, b = x
        t = self._m(a, b)
        return (self._m(self._a(a, b), self._s(a, b)), self._a(t, t))
    def _add2(self, x, y): return (self._a(x[0], y[0]), self._a(x[1], y[1]))
    def _sub2(self, x, y): return (self._s(x[0], y[0]), self._s(x[1], y[1]))
    def _xi2(self, x): return (self._s(x[0], x[1]), self._a(x[0], x[1]))
    def _sqr4(self, a0, a1):
        t0 = self._sqr2(a0); t1 = self._sqr2(a1)
        sm = self._add2(a0, a1)
        r0 = self._add2(self._xi2(t1), t0)
        r1 = self._sub2(self._sub2(self._sqr2(sm), t0), t1)
        return r0, r1
    def cyc_sqr_counted(self, x):
        t0 = self._sqr4(x[0], x[3]); t1 = self._sqr4(x[1], x[4]); t2 = self._sqr4(x[2], x[5])
        def m3(t, a, sg):
            u = self._sub2(t, a) if sg < 0 else self._add2(t, a)
            return self._add2(self._add2(u, u), t)
        return (m3(t0[0], x[0], -1), m3(self._xi2(t2[1]), x[1], 1), m3(t1[0], x[2], -1),
                m3(t0[1], x[3], 1), m3(t2[0], x[4], -1), m3(t1[1], x[5], 1))


def isprime(n):
    if n < 2: return False
    for pr in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % pr == 0: return n == pr
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True


# ---------- production constants (proof.md sections 3-5) ----------
def up(x, n): return math.ceil(x * 10 ** n) / 10 ** n      # costs, bounds, rates: rounded up
def dn(x, n): return math.floor(x * 10 ** n) / 10 ** n     # success figures, margins, lower bounds: rounded down
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
C = 20000
R_PROD, M_WALKS, THETA_LG = 1024, 2 ** 20, 60
# Branch fraction: square iff the 6-bit branch field of T0 is below 35, i.e. s = 35/64 = 0.546875,
# the W-minimising exactly stateable fraction (continuous argmin 0.5461).
BRANCH_BITS, BRANCH_THRESHOLD = 6, 35
S_FRAC = Fraction(BRANCH_THRESHOLD, 1 << BRANCH_BITS)
# Meter increments, instructions: 20,000 + 2,250 for a multiply step; 10,098 = 9,798 measured Callgrind
# instructions of blst_fp12_cyclotomic_sqr (the cyclotomic squaring, valid on GT) + 300 for the alpha/beta
# doubling, the shared rows and slack.
I_MUL, I_SQ = 22250, 10098
N_PROD = (r - 1) // 12
COLL = 1 - S_FRAC ** 2 - (1 - S_FRAC) ** 2 / R_PROD          # charged pair-collision factor
lgN = math.log2(N_PROD)
lgNp = lgN - math.log2(float(COLL))
lgS0 = 0.5 * (1 + lgNp + math.log2(math.log(1 / 0.61)))
S_SLOT = math.ceil(1.03 * 2 ** (lgS0 - 20) * (1 + 1e-9))      # per-slot productive step budget (integer)
S_STEPS = M_WALKS * S_SLOT
E_SLOT = 20 * 2 ** THETA_LG                                  # detection window per slot, 20/theta steps
RW_LG = 47
mean_reports_slot = (S_SLOT + E_SLOT) * 2.0 ** -THETA_LG
CAP_RATIO = 2 ** RW_LG / mean_reports_slot                    # R_w over the healthy per-slot report mean

# ---------- Part A: production BLS12-381 identities ----------
P = Tower(p)
zeta = P.zeta
hcof = (p ** 12 - 1) // r
R1_FACTORS = [(2, 32), (3, 1), (11, 1), (19, 1), (10177, 1), (125527, 1), (859267, 1), (906349, 2),
              (2508409, 1), (2529403, 1), (52437899, 1), (254760293, 2)]
prod_fac = 1
for f, e in R1_FACTORS: prod_fac *= f ** e
ord2 = r - 1
for f, e in R1_FACTORS:
    while ord2 % f == 0 and pow(2, ord2 // f, r) == 1: ord2 //= f
g12 = math.gcd(ord2, 12)
partA = {"zeta_in_Fp_primitive_6th_root": pow(zeta, 6, p) == 1 and pow(zeta, 3, p) != 1 and pow(zeta, 2, p) != 1,
         "p_order_mod_r_is_12": pow(p, 12, r) == 1 and all(pow(p, k, r) != 1 for k in range(1, 12)),
         "one_plus_p4_plus_p8_zero_mod_r": (1 + pow(p, 4, r) + pow(p, 8, r)) % r == 0,
         "r_minus_1_factorisation_verified": prod_fac == r - 1 and all(isprime(f) for f, _ in R1_FACTORS),
         "ord_r_2_equals_(r-1)/16": ord2 == (r - 1) // 16 and pow(2, ord2, r) == 1,
         "log2_ord_r_2": dn(math.log2(ord2), 4),
         "min_squarings_in_a_fruitless_cycle_log2": dn(math.log2(ord2 // g12), 2)}
ok_keys = ok_frob = ok_canon = ok_gt = ok_inv = ok_min = ok_t0 = ok_gs = ok_gs_fast = True
gs_ops = None
GT_POOL = []
for _ in range(prod_checks):
    x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
    y = P.pow(x, hcof)
    GT_POOL.append(y)
    ok_gt &= P.pow(y, r) == P.one() and y != P.one()
    ok_frob &= P.frob(y) == P.pow(y, p)
    kd, kf = P.keys_direct(y), P.keys_formula(y)
    ok_keys &= kd == kf and len(set(kd)) == 12
    ks, Ks = P.key_min_structured(y)
    ok_min &= (len(set(kf)) == 12 and ks == min(range(12), key=lambda i: kf[i]) and Ks == min(kf))
    conj = [y]
    for k in range(11): conj.append(P.frob(conj[-1]))
    ok_t0 &= len({c[0][0] for c in conj}) == 1                       # T0 = Re(a_0) equal on the class
    canon = set()
    for c in conj:
        kk = P.keys_formula(c); kap = min(range(12), key=lambda i: kk[i])
        t = c
        for _ in range(kap): t = P.frob(t)
        canon.add(t)
    ok_canon &= len(canon) == 1
    ok_inv &= P.mul(conj[6], y) == P.one()
    P.ops.clear()
    sq = P.cyc_sqr_counted(y)
    gs_ops = dict(P.ops)
    ok_gs &= sq == P.mul(y, y)
    ok_gs_fast &= P.cyc_sqr(y) == sq
    ok_gs &= P.cyc_sqr(conj[3]) == P.mul(conj[3], conj[3])
partA.update({"random_GT_elements": prod_checks, "GT_membership": ok_gt,
              "frobenius_formula_equals_pth_power": ok_frob,
              "four_mult_key_formula_equals_direct_keys_and_distinct": ok_keys,
              "sign_set_minimum_equals_12way_minimum_and_kappa": ok_min,
              "kappa_canonical_rep_identical_over_all_12_conjugates": ok_canon,
              "p6_power_is_inverse": ok_inv,
              "T0_equal_on_all_12_conjugates": ok_t0,
              "granger_scott_equals_y_times_y": ok_gs, "granger_scott_fast_equals_counted": ok_gs_fast,
              "granger_scott_counted_ops": gs_ops})

# ---------- Part A2: low bits of T0 and of K on real GT elements ----------
def partA2(budget_s):
    t_end = time.time() + budget_s
    y0 = GT_POOL[0]
    pool = [P.pow(y0, rng.randrange(1, r)) for _ in range(24)]
    ones = [0] * 76; kones = [0] * 128; nsamp = 0; sq_hits = 0; dp4 = 0
    while nsamp < 2000 and time.time() < t_end:
        y = P.mul(P.mul(rng.choice(pool), rng.choice(pool)), rng.choice(pool))
        if y == P.one(): continue
        t0 = y[0][0]; K = min(P.keys_formula(y)); nsamp += 1
        for b in range(76): ones[b] += (t0 >> b) & 1
        for b in range(128): kones[b] += (K >> b) & 1
        sq_hits += ((t0 >> 10) & ((1 << BRANCH_BITS) - 1)) < BRANCH_THRESHOLD
        dp4 += ((t0 >> (10 + BRANCH_BITS)) & 15) == 0
    if not nsamp: return {"samples": 0}
    se = 0.5 / math.sqrt(nsamp)
    dev = [abs(o / nsamp - 0.5) / se for o in ones]
    kdev = [abs(o / nsamp - 0.5) / se for o in kones]
    return {"samples": nsamp, "T0_bits_tested": 76,
            "T0_max_abs_deviation_in_standard_errors": round(max(dev), 2),
            "T0_mean_abs_deviation_in_standard_errors": round(sum(dev) / 76, 2),
            "K_record_key_bits_tested": 128,
            "K_max_abs_deviation_in_standard_errors": round(max(kdev), 2),
            "K_mean_abs_deviation_in_standard_errors": round(sum(kdev) / 128, 2),
            "expected_mean_abs_deviation_if_uniform": round(math.sqrt(2 / math.pi), 2),
            "squaring_field_fraction_bits_10_15_below_35": round(sq_hits / nsamp, 4),
            "nominal_35/64": round(float(S_FRAC), 6),
            "fraction_bits_16_19_zero": round(dp4 / nsamp, 4), "nominal_2^-4": 2.0 ** -4}

# ---------- Part P: the production step rule on real BLS12-381 GT elements ----------
def partP(budget_s):
    """Two walkers start at y0 and at its conjugate y0^(p^5) and run the production step rule in lockstep:
    partition h = bits 0..3 of T0 (a 16-entry table here; production uses bits 0..9 and 1,024 entries),
    squaring iff bits 10..15 of T0 are below 35, distinguished iff bits 16..19 of T0 are zero (4 bits here
    so that reports occur; production uses bits 16..75). Keys are computed only on multiplication steps and
    at reports. Cycle checks compare T0 two steps back and at a checkpoint every 64 steps.

    The second test is the coalescence test of proof.md section 3.1: an all-multiply 2-cycle C1 <-> C2 is
    built by choosing the two table cells it needs, and two walkers start at class C1 with different
    histories (one with the cycle already in its two-step T0 ring, one freshly started at a different
    conjugate of C1). They must escape at the same class and share their class sequence from there,
    within one cycle length of each other."""
    t_end = time.time() + budget_s
    g = GT_POOL[0]; xs = rng.randrange(1, r); z = P.pow(g, xs)
    RP = 16
    pk = [pow(p, k, r) for k in range(12)]
    ah = [rng.randrange(1, 2 ** 40) for _ in range(RP)]; bh = [rng.randrange(1, 2 ** 40) for _ in range(RP)]

    def make_table(ah_, bh_):
        tb = []
        for i in range(RP):
            row = [P.mul(P.pow(g, ah_[i]), P.pow(z, bh_[i]))]
            for _k in range(11): row.append(P.frob(row[-1]))
            tb.append(row)
        return tb

    Tab = make_table(ah, bh)

    def one_step(w, y, Tab_, ah_, bh_):
        t0 = y[0][0]
        if ((t0 >> 10) & ((1 << BRANCH_BITS) - 1)) < BRANCH_THRESHOLD:
            w["sq"] += 1
            return P.cyc_sqr(y), 2 * w["a"] % r, 2 * w["b"] % r
        kf = P.keys_formula(y); w["keys"] += 1; w["mul"] += 1
        kap = min(range(12), key=lambda i: kf[i]); t = (12 - kap) % 12; h = t0 & (RP - 1)
        return P.mul(y, Tab_[h][t]), (w["a"] + ah_[h] * pk[t]) % r, (w["b"] + bh_[h] * pk[t]) % r

    def ckey(y):
        """Class identity: the free invariant T0 and the class key K (both equal on all 12 conjugates)."""
        return (y[0][0], min(P.keys_formula(y)))

    def walk_iters(w, Tab_, ah_, bh_, limit, deadline=None, classes=False):
        """The production order one iteration at a time: 2-cycle check (T0 two steps back), checkpoint check
        (T0 stored every 64 steps), distinguished test and report, then the branch step. Returns the T0
        sequence before each iteration, the class sequence after each (only when classes is set), and the
        iterations on which the 2-cycle escape fired."""
        pre = []; post = []; esc = []; i = 0
        while i < limit and (deadline is None or time.time() < deadline):
            y = w["y"]; t0 = y[0][0]; hist = w["hist"]
            if len(hist) >= 2 and hist[-2] == t0:
                w["two"] += 1
                if hist[-1] < t0:
                    y, w["a"], w["b"] = one_step(w, y, Tab_, ah_, bh_)
                w["y"] = P.cyc_sqr(y); w["a"] = 2 * w["a"] % r; w["b"] = 2 * w["b"] % r
                w["hist"] = []; w["ck"] = None
                pre.append(t0)
                if classes: post.append(ckey(w["y"]))
                esc.append(i); i += 1
                continue
            if w["ck"] is not None and w["ck"][0] == t0 and w["ck"][1] < w["n"]:
                w["longc"] += 1
                w["hist"] = []; w["ck"] = None
            if ((t0 >> (10 + BRANCH_BITS)) & 15) == 0:
                P.keys_formula(y); w["keys"] += 1; w["reports"] += 1
            if w["n"] % 64 == 0: w["ck"] = (t0, w["n"])
            hist.append(t0)
            if len(hist) > 2: hist.pop(0)
            w["y"], w["a"], w["b"] = one_step(w, y, Tab_, ah_, bh_)
            w["n"] += 1
            pre.append(t0)
            if classes: post.append(ckey(w["y"]))
            i += 1
        return pre, post, esc

    def new_walker(y, a, b, hist=None, n=0):
        return {"y": y, "a": a, "b": b, "hist": list(hist or []), "ck": None, "n": n,
                "sq": 0, "mul": 0, "keys": 0, "reports": 0, "two": 0, "longc": 0, "t0s": []}

    a0, b0 = rng.randrange(r), rng.randrange(r)
    y0 = P.mul(P.pow(g, a0), P.pow(z, b0))
    y5 = y0
    for _ in range(5): y5 = P.frob(y5)
    walkers = [new_walker(y0, a0, b0), new_walker(y5, a0 * pk[5] % r, b0 * pk[5] % r)]
    steps = 0
    while time.time() < t_end and steps < 4000:
        for w in walkers:
            pre, _post, _esc = walk_iters(w, Tab, ah, bh, 1)
            w["t0s"].extend(pre)
        steps += 1
    same_seq = walkers[0]["t0s"] == walkers[1]["t0s"]
    inv = all(w["y"] == P.mul(P.pow(g, w["a"]), P.pow(z, w["b"])) for w in walkers)
    w = walkers[0]

    # --- the same class from different histories: 2-cycle coalescence (proof.md section 3.1) ---
    # Build an all-multiply 2-cycle C1 <-> C2 by setting the two table cells the cycle uses; both rows are
    # rebuilt from the same base multipliers, so y = g^alpha z^beta still holds and the exponent updates
    # charge exactly the multipliers applied. Walker D holds the pair in its two-step T0 ring at C1 and
    # escapes on its first iteration; walker F starts at a different conjugate of C1 with an empty ring,
    # so it must walk the cycle once before detecting the return.
    def invm(x): return pow(x, r - 2, r)

    def plan(y):
        t0 = y[0][0]
        if ((t0 >> 10) & ((1 << BRANCH_BITS) - 1)) < BRANCH_THRESHOLD: return None
        kf = P.keys_formula(y)
        return (t0 & (RP - 1), (12 - min(range(12), key=lambda i: kf[i])) % 12)

    cycle = {"constructed": False}
    yc1 = None
    for _ in range(64):
        ac0, bc0 = rng.randrange(1, r), rng.randrange(1, r)
        yc1 = P.mul(P.pow(g, ac0), P.pow(z, bc0))
        if plan(yc1) is not None: break
    if plan(yc1) is not None:
        h1c, t1c = plan(yc1)
        ahc = [rng.randrange(1, 2 ** 40) for _ in range(RP)]
        bhc = [rng.randrange(1, 2 ** 40) for _ in range(RP)]
        built = None
        for _ in range(512):
            U, V = rng.randrange(1, r), rng.randrange(1, r)
            Mc = P.mul(P.pow(g, U), P.pow(z, V))
            yc2 = P.mul(yc1, Mc)
            p2 = plan(yc2)
            if p2 is None or p2[0] == h1c: continue
            built = (U, V, Mc, yc2, p2[0], p2[1])
            break
        if built:
            U, V, Mc, yc2, h2c, t2c = built
            r1, r2 = invm(pow(p, t1c, r)), invm(pow(p, t2c, r))
            ahc[h1c], bhc[h1c] = U * r1 % r, V * r1 % r
            ahc[h2c], bhc[h2c] = (-U) * r2 % r, (-V) * r2 % r
            Tabc = make_table(ahc, bhc)
            Minv = P.pow(Mc, r - 1)
            cycle_ok = Tabc[h1c][t1c] == Mc and Tabc[h2c][t2c] == Minv and P.mul(Mc, Minv) == P.one()
            yf = yc1; af, bf = ac0, bc0
            for _ in range(5):
                yf = P.frob(yf); af = af * p % r; bf = bf * p % r
            wd = new_walker(yc1, ac0, bc0, hist=[yc1[0][0], yc2[0][0]], n=2)   # detector: cycle in the ring
            wf = new_walker(yf, af, bf, hist=[], n=0)                         # fresh at the same class
            _pd, post_d, esc_d = walk_iters(wd, Tabc, ahc, bhc, 12, classes=True)
            _pf, post_f, esc_f = walk_iters(wf, Tabc, ahc, bhc, 12, classes=True)
            small = yc1 if yc1[0][0] < yc2[0][0] else yc2
            target = ckey(P.cyc_sqr(small))       # the smaller-T0 class of the pair, squared
            if esc_d and esc_f:
                i_d, i_f = esc_d[0], esc_f[0]
                tail_d, tail_f = post_d[i_d:], post_f[i_f:]
                k = min(len(tail_d), len(tail_f))
                same_tail = k > 0 and tail_d[:k] == tail_f[:k]
                same_escape = post_d[i_d] == target and post_f[i_f] == target
                inv_c = (wd["y"] == P.mul(P.pow(g, wd["a"]), P.pow(z, wd["b"]))
                         and wf["y"] == P.mul(P.pow(g, wf["a"]), P.pow(z, wf["b"])))
                gap = i_f - i_d
                cycle = {"constructed": bool(cycle_ok), "cycle_length": 2,
                         "table_cells_hold_M_and_M_inverse": bool(cycle_ok),
                         "histories_differ_at_start": wd["hist"] != wf["hist"],
                         "detector_escape_iteration": i_d, "fresh_escape_iteration": i_f,
                         "coalesce_gap_within_cycle_length": 0 <= gap <= 2,
                         "both_escape_at_smaller_T0_class_squared": bool(same_escape),
                         "same_class_sequence_after_coalescence": bool(same_tail),
                         "compared_class_steps_after_coalescence": k,
                         "distinct_exponent_pairs_after_coalescence": (wd["a"], wd["b"]) != (wf["a"], wf["b"]),
                         "exponent_invariant_both_walkers": bool(inv_c),
                         "coalesced_from_different_histories": bool(cycle_ok and same_escape and 0 <= gap <= 2
                                                                   and same_tail and inv_c)}

    return {"lockstep_steps": steps, "table_entries": RP, "multiplier_exponent_bits": 40,
            "two_walkers_from_conjugate_starts_same_T0_sequence": same_seq,
            "exponent_invariant_y_equals_g^alpha_z^beta_both_walkers": inv,
            "squaring_steps": w["sq"], "multiplication_steps": w["mul"],
            "squaring_fraction": round(w["sq"] / max(1, w["sq"] + w["mul"]), 4),
            "key_computations": w["keys"],
            "key_computations_equal_mult_steps_plus_reports": w["keys"] == w["mul"] + w["reports"],
            "reports_4bit_test": w["reports"], "two_cycle_escapes": w["two"], "checkpoint_events": w["longc"],
            "different_histories_two_cycle_coalescence": cycle}


# ---------- toy fields ----------
# q = 7 mod 12, xi = 1+u neither a square nor a cube in F_{q^2} (asserted), r' a prime factor of
# q^4 - q^2 + 1 with ord_{r'}(q) = 12 (asserted), N = (r'-1)/12 classes.
TOYS = {67: 1549741, 163: 54298861, 211: 8655349, 1048963: 48674953, 1049611: 64668529, 1054003: 35456053}
CONFIGS = [(1048963, 64, 4), (1049611, 256, 4), (1054003, 1024, 3)]
SLICES = [0.25, 0.25, 0.5]


def binom_err(f, n):
    return math.sqrt(max(0.0, f * (1.0 - f)) / n) if n else None


class Toy:
    """Toy field, toy GT of prime order r', its Frobenius machinery and the class-bit rule."""
    def __init__(self, q, R, dp_bits, bits_mode):
        self.q, self.R, self.dp_bits, self.bits_mode = q, R, dp_bits, bits_mode
        T = self.T = Tower(q)
        rp = self.rp = TOYS[q]
        assert (q ** 4 - q * q + 1) % rp == 0 and isprime(rp)
        assert pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))
        assert T.f2pow(T.xi, (q * q - 1) // 2) != (1, 0) and T.f2pow(T.xi, (q * q - 1) // 3) != (1, 0)
        self.N = (rp - 1) // 12
        self.qk = [pow(q, k, rp) for k in range(12)]
        self.hc = (q ** 12 - 1) // rp
        self.Rbits = R.bit_length() - 1
        assert R == 1 << self.Rbits
        if bits_mode == "t0":
            assert self.Rbits + BRANCH_BITS + dp_bits <= q.bit_length() - 1, "toy T0 too short for the bit fields"
        self.DPMASK = (1 << dp_bits) - 1
        self.Dk = [[T.f2pow(T.xi, m * (q ** k - 1) // 6) for m in range(6)] for k in range(12)]
        self.ties = 0
    def frobk(self, y, k):
        q = self.q
        if k % 2: y = tuple((c[0], (-c[1]) % q) for c in y)
        return tuple(self.T.f2mul(c, self.Dk[k][m]) for m, c in enumerate(y))
    def keyinfo(self, y):
        """kappa, class key K and canonical conjugate (toy ties, about 66/q per step, broken by the full
        conjugate; production ties are a 2^-249 failure term)."""
        kf = self.T.keys_formula(y)
        m = min(kf)
        ks = [k for k in range(12) if kf[k] == m]
        if len(ks) == 1:
            kap = ks[0]; t = self.frobk(y, kap)
        else:
            self.ties += 1
            t, kap = min((self.frobk(y, k), k) for k in ks)
        return kap, m, t
    def bits(self, y, cls):
        """(h, squaring?, distinguished?) from disjoint bit fields: bits [0, log2 R) partition, the next
        BRANCH_BITS bits the branch field (square iff < BRANCH_THRESHOLD), the next dp_bits bits the
        distinguished field. Mode t0
        reads them from T0 = Re(a_0) of the element (the production rule); mode hash from a hash of the
        canonical conjugate (for toy fields whose T0 is too short)."""
        if self.bits_mode == "t0":
            v = y[0][0]
        else:
            v = int.from_bytes(hashlib.blake2b(repr(cls).encode(), digest_size=8).digest(), "little")
        Rb = self.Rbits
        return (v & (self.R - 1), ((v >> Rb) & ((1 << BRANCH_BITS) - 1)) < BRANCH_THRESHOLD,
                ((v >> (Rb + BRANCH_BITS)) & self.DPMASK) == 0)
    def generator(self, rr):
        while True:
            g = self.T.pow(tuple((rr.randrange(self.q), rr.randrange(self.q)) for _ in range(6)), self.hc)
            if g != self.T.one(): return g


def pair_test(q, R, nsamp, bits_mode, deadline, rr):
    """Part E: sample distinct classes uniformly, apply one step to each, count image coincidences.
    Under H1 two distinct classes have equal images with probability c/N where c is 0 for two squaring
    classes (squaring is injective on classes), 1 for a squaring and a multiplication class, 11/12 for two
    multiplication classes with the same partition index (the k = 0 twist would need equal classes) and 1
    for different indices; pooled N x rate = 1 - s^2 - (1-s)^2/(12R)."""
    tf = Toy(q, R, 0, bits_mode)
    T, rp, N = tf.T, tf.rp, tf.N
    g = tf.generator(rr)
    gpow = [g]
    for _ in range(rp.bit_length()): gpow.append(T.mul(gpow[-1], gpow[-1]))
    def gexp(e):
        acc = T.one(); j = 0
        while e:
            if e & 1: acc = T.mul(acc, gpow[j])
            e >>= 1; j += 1
        return acc
    Mh = [gexp(rr.randrange(1, rp)) for _ in range(R)]
    seen = set(); recs = []
    while len(recs) < nsamp and time.time() < deadline:
        y = gexp(rr.randrange(1, rp))
        kap, K, cls = tf.keyinfo(y)
        if cls in seen: continue
        seen.add(cls)
        h, sq, _ = tf.bits(y, cls)
        img = T.cyc_sqr(cls) if sq else T.mul(cls, Mh[h])
        recs.append((tf.keyinfo(img)[2], sq, h))
    n = len(recs)
    byimg = defaultdict(list)
    for icls, sq, h in recs: byimg[icls].append((sq, h))
    co = Counter()
    for lst in byimg.values():
        for i in range(len(lst)):
            for j in range(i + 1, len(lst)):
                (s1, h1), (s2, h2) = lst[i], lst[j]
                co["sq_sq" if s1 and s2 else "sq_mul" if (s1 or s2) else
                   ("mul_mul_same_h" if h1 == h2 else "mul_mul_diff_h")] += 1
    nsq = sum(1 for x in recs if x[1]); nmul = n - nsq
    mh = Counter(x[2] for x in recs if not x[1])
    pairs = {"sq_sq": nsq * (nsq - 1) // 2, "sq_mul": nsq * nmul,
             "mul_mul_same_h": sum(c * (c - 1) // 2 for c in mh.values())}
    pairs["mul_mul_diff_h"] = nmul * (nmul - 1) // 2 - pairs["mul_mul_same_h"]
    pred = {"sq_sq": 0.0, "sq_mul": 1.0, "mul_mul_same_h": 11 / 12, "mul_mul_diff_h": 1.0}
    by_type = {}
    for k in pred:
        pr = pairs[k]
        val = N * co[k] / pr if pr else None
        se = N * math.sqrt(max(pred[k], 1e-12) * pr / N) / pr if pr else None
        by_type[k] = {"pairs": pr, "coincidences": co[k], "N_x_rate": round(val, 4) if val is not None else None,
                      "predicted": round(pred[k], 4), "z": round((val - pred[k]) / se, 2) if pr and pred[k] else None}
    tot_pairs = n * (n - 1) // 2
    s = float(S_FRAC)
    pred_formula = 1 - s * s - (1 - s) ** 2 / (12 * R)
    pred_emp = sum(pred[k] * pairs[k] for k in pred) / tot_pairs if tot_pairs else 0
    val = N * sum(co.values()) / tot_pairs if tot_pairs else None
    se = math.sqrt(pred_emp * tot_pairs / N) * N / tot_pairs if tot_pairs else None
    return {"q": q, "classes_N": N, "R": R, "bits": bits_mode, "s": s, "sampled_classes": n,
            "squaring_fraction_of_sample": round(nsq / max(1, n), 4),
            "squaring_fraction_binomial_err": round(binom_err(nsq / max(1, n), max(1, n)), 4),
            "pooled_N_x_rate": round(val, 4) if val is not None else None,
            "pooled_standard_error": round(se, 4) if se else None,
            "predicted_1-s^2-(1-s)^2/(12R)": round(pred_formula, 4),
            "predicted_at_sample_branch_mix": round(pred_emp, 4),
            "z_vs_sample_mix": round((val - pred_emp) / se, 2) if se else None,
            "charged_1-s^2-(1-s)^2/R": round(1 - s * s - (1 - s) ** 2 / R, 4),
            "by_pair_type": by_type, "sq_sq_coincidences_exactly_zero": co["sq_sq"] == 0}


def run_toy(q, R, dp_bits, trials_target, deadline, rr):
    tf = Toy(q, R, dp_bits, toy_bits)
    T, rp, N, qk = tf.T, tf.rp, tf.N, tf.qk
    nbits = rp.bit_length()
    STUCK = 40 << dp_bits
    s = float(S_FRAC)
    coll = 1 - s * s - (1 - s) ** 2 / R
    Np = N / coll
    b39 = math.sqrt(2 * Np * math.log(1 / 0.61))
    budget_steps = 1.03 * b39
    horizon = 4.0 * math.sqrt(Np) if long_horizon else budget_steps
    theta = 2.0 ** -dp_bits
    cap_walk = max(1, int(math.ceil(CAP_RATIO * (budget_steps / walks) * theta)))
    log_cap_total = walks * cap_walk
    pair_budget_toy = 2 ** 36 * (log_cap_total * (log_cap_total - 1) // 2)
    cavg = (1 - s) * I_MUL + s * I_SQ
    sig_ind = (I_MUL - I_SQ) * math.sqrt(s * (1 - s)) / cavg
    win_extra = int(math.ceil(20 / theta))          # production detection window E, in toy steps
    win_end = budget_steps + win_extra

    def classinfo(y):
        kap, K, cls = tf.keyinfo(y)
        h, sq, dp = tf.bits(y, cls)
        return kap, K, cls, y[0][0], h, sq, dp

    def apply_step(y, a, b, info, Tab, At, Bt):
        kap, h, sq = info[0], info[4], info[5]
        if sq:
            return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp, True
        idx = (12 - kap) % 12
        return T.mul(y, Tab[h][idx]), (a + At[h][idx]) % rp, (b + Bt[h][idx]) % rp, False

    def walk_pass(cap, step_budget, instrument, g, z, Tab, At, Bt, stuck=None, meter_cap=None, window=0):
        """One trial's walk phase. Production order inside a step: (1) the slot's meter is tested before
        every step, escape squaring and re-step; (2) the 2-cycle check and then the checkpoint check run
        BEFORE the distinguished test; (3) the append rule; (4) the step: squaring branch or twisted
        multiplication, decided by the class bits. Cycle escapes use the same mixed step function. The
        trial runs its full fixed step budget (no stop at the first recovery); with `instrument` a
        measurement-only class index timestamps the first useful class collision and, past the budget,
        the walks continue without appending until that collision is seen or the horizon is reached."""
        visited = {} if instrument else None
        stuck_steps = STUCK if stuck is None else stuck
        prod_end = step_budget + window          # production detection window after the step budget
        hz = max(horizon, prod_end) if instrument else prod_end
        steps = 0; Tc = None; log = []
        st = Counter()
        W = [{"rep": 0, "dead": False, "meter": 0, "nmul": 0, "nsq": 0, "len": 0, "restarts": 0} for _ in range(walks)]
        def fresh(w):
            a, b = rr.randrange(rp), rr.randrange(rp)
            w.update(y=T.mul(T.pow(g, a), T.pow(z, b)), a=a, b=b, hist=[], ck=None, sincedp=0, last=None)
        for w in W: fresh(w)
        def can_step(w):
            if meter_cap is not None and w["meter"] >= meter_cap:
                if not w["dead"]: st["meter_stops"] += 1
                w["dead"] = True
                return False
            return True
        def charge(w, sq, prod):
            if prod:
                w["meter"] += I_SQ if sq else I_MUL
                w["nsq" if sq else "nmul"] += 1
        while steps < hz:
            if steps >= prod_end and (Tc is not None or not instrument): break
            live = [w for w in W if not w["dead"]]
            if not live: break
            for w in live:
                if steps >= hz: break
                prod = steps < prod_end
                if not prod and (Tc is not None or not instrument): break
                if prod and not can_step(w): continue
                y = w["y"]; info = classinfo(y); kap, K, cls, t0, h, sq, dp = info
                ca, cb = w["a"] * qk[kap] % rp, w["b"] * qk[kap] % rp
                hist = w["hist"]
                if len(hist) >= 2 and hist[-2][1] == cls:            # 2-cycle check (production: T0 compare)
                    st["two_cycles"] += 1
                    yy, aa, bb = y, w["a"], w["b"]
                    if hist[-1] < (t0, cls):                         # smaller class is the next one: re-step
                        yy, aa, bb, s1 = apply_step(yy, aa, bb, info, Tab, At, Bt)
                        charge(w, s1, prod); steps += 1
                        if prod and not can_step(w): continue
                    w["y"] = T.cyc_sqr(yy); w["a"] = 2 * aa % rp; w["b"] = 2 * bb % rp
                    charge(w, True, prod); steps += 1
                    w["hist"] = []; w["ck"] = None
                    continue
                if w["ck"] is not None and w["ck"][0] == cls and w["ck"][1] < w["len"]:
                    st["long_cycles"] += 1                           # checkpoint check (production: T0 compare)
                    best = ((t0, cls), y, w["a"], w["b"]); yy, aa, bb, inf = y, w["a"], w["b"], info
                    stopped = False
                    for _ in range(w["len"] - w["ck"][1]):           # once around with the MIXED step function
                        yy, aa, bb, s1 = apply_step(yy, aa, bb, inf, Tab, At, Bt)
                        charge(w, s1, prod); steps += 1
                        if s1: st["squarings_in_cycle_traversal"] += 1
                        if prod and not can_step(w): stopped = True; break
                        inf = classinfo(yy)
                        if (inf[3], inf[2]) < best[0]: best = ((inf[3], inf[2]), yy, aa, bb)
                    if stopped: continue
                    _, yy, aa, bb = best
                    w["y"] = T.cyc_sqr(yy); w["a"] = 2 * aa % rp; w["b"] = 2 * bb % rp
                    charge(w, True, prod); steps += 1
                    w["hist"] = []; w["ck"] = None
                    continue
                if visited is not None and Tc is None:               # measurement only
                    prev = visited.get(cls)
                    if prev is None: visited[cls] = (ca, cb)
                    elif prev != (ca, cb): Tc = steps
                if dp:
                    w["sincedp"] = 0
                    if prod:
                        rec = (K, ca, cb)                            # record: class key, canonical exponents
                        if rec != w["last"]:
                            log.append(rec); w["last"] = rec; w["rep"] += 1
                            if w["rep"] >= cap:
                                st["cap_stops"] += 1; w["dead"] = True
                                continue
                        else:
                            st["duplicate_appends_skipped"] += 1
                if w["len"] % 64 == 0: w["ck"] = (cls, w["len"])
                hist.append((t0, cls))
                if len(hist) > 2: hist.pop(0)
                w["y"], w["a"], w["b"], s1 = apply_step(y, w["a"], w["b"], info, Tab, At, Bt)
                charge(w, s1, prod)
                w["len"] += 1; w["sincedp"] += 1; steps += 1
                if w["sincedp"] > stuck_steps:
                    st["abandoned"] += 1; w["restarts"] += 1; fresh(w)
        return {"steps": steps, "Tc": Tc, "log": log, "events": st, "slots": W,
                "sum_slot_reports": sum(w["rep"] for w in W)}

    def elem(A, B, gp, zp):
        acc = T.one(); j = 0
        while (A or B) and j < len(gp):
            if A & 1: acc = T.mul(acc, gp[j])
            if B & 1: acc = T.mul(acc, zp[j])
            A >>= 1; B >>= 1; j += 1
        return acc

    def match_and_verify(log, pair_budget, g, z, gp, zp):
        """Production matcher: sort by the full record (key, A, B), drop adjacent exact duplicates, charge
        every pair of the remaining distinct records of each equal-key run against the hard pair budget,
        recompute g^A z^B for each record a charged pair touches, solve x = (A1-A2)/(B2-B1) from the first
        pair of equal elements with different exponent pairs, accept x only if g^x = z."""
        recs = sorted(log)
        pairs = 0; recomputations = 0; identical = 0; recovered = None; cap_hit = False
        i = 0
        while i < len(recs) and recovered is None and not cap_hit:
            j = i + 1
            while j < len(recs) and recs[j][0] == recs[i][0]: j += 1
            run = recs[i:j]
            if len(run) >= 2:
                distinct = []
                for rec in run:
                    if distinct and rec == distinct[-1]: identical += 1
                    else: distinct.append(rec)
                elems = [None] * len(distinct)
                for u in range(len(distinct)):
                    for v in range(u + 1, len(distinct)):
                        if pairs >= pair_budget:
                            cap_hit = True; break
                        pairs += 1
                        for t in (u, v):
                            if elems[t] is None:
                                elems[t] = elem(distinct[t][1], distinct[t][2], gp, zp); recomputations += 1
                        if elems[u] != elems[v]: continue
                        (a1, b1), (a2, b2) = distinct[u][1:], distinct[v][1:]
                        if b1 == b2: continue
                        xr = (a1 - a2) * pow(b2 - b1, -1, rp) % rp
                        if T.pow(g, xr) == z:
                            recovered = xr; break
                    if recovered is not None or cap_hit: break
            i = j
        return {"recovered": recovered, "pairs": pairs, "recomputations": recomputations,
                "identical_records": identical, "cap_hit": cap_hit}

    def instance():
        g = tf.generator(rr)
        xs = rr.randrange(rp); z = T.pow(g, xs)
        ah = [rr.randrange(rp) for _ in range(R)]; bh = [rr.randrange(rp) for _ in range(R)]
        gpow = [g]
        for _ in range(nbits): gpow.append(T.mul(gpow[-1], gpow[-1]))
        zp = [z]
        for _ in range(nbits): zp.append(T.mul(zp[-1], zp[-1]))
        Tab = []; At = []; Bt = []
        for i in range(R):
            e = (ah[i] + xs * bh[i]) % rp; acc = T.one(); j = 0
            while e:
                if e & 1: acc = T.mul(acc, gpow[j])
                e >>= 1; j += 1
            row = [acc]
            for k in range(11): row.append(T.frob(row[-1]))
            Tab.append(row); At.append([ah[i] * qk[k] % rp for k in range(12)]); Bt.append([bh[i] * qk[k] % rp for k in range(12)])
        return g, z, gpow, zp, Tab, At, Bt

    n_trials = n_rec = n_coll_in = n_coll_no_rec = n_cap_hit = 0
    Tc_in = []; Tc_full = []; steps_tot = 0; reports_tot = 0; max_log = 0; pairs_max = 0
    stats = Counter(); dev_norm = []; nsq_tot = nmul_tot = 0
    while n_trials < trials_target and time.time() < deadline:
        g, z, gpow, zp, Tab, At, Bt = instance()
        res = walk_pass(cap_walk, budget_steps, True, g, z, Tab, At, Bt,
                        meter_cap=ceil(1.01 * cavg * win_end / walks), window=win_extra)
        m = match_and_verify(res["log"], pair_budget_toy, g, z, gpow, zp)
        assert m["recovered"] is None or T.pow(g, m["recovered"]) == z
        assert len(res["log"]) <= log_cap_total and len(res["log"]) == res["sum_slot_reports"]
        n_trials += 1
        tc = res["Tc"]
        if tc is not None: Tc_full.append(tc)
        if tc is not None and tc <= budget_steps: n_coll_in += 1; Tc_in.append(tc)
        if m["recovered"] is not None: n_rec += 1
        if tc is not None and tc <= budget_steps and m["recovered"] is None: n_coll_no_rec += 1
        n_cap_hit += m["cap_hit"]
        pairs_max = max(pairs_max, m["pairs"])
        stats.update(res["events"])
        for w in res["slots"]:
            nst = w["nmul"] + w["nsq"]
            nsq_tot += w["nsq"]; nmul_tot += w["nmul"]
            if nst: dev_norm.append((w["meter"] / nst - cavg) / cavg * math.sqrt(nst))
        steps_tot += min(res["steps"], int(win_end) + 1)
        reports_tot += len(res["log"]); max_log = max(max_log, len(res["log"]))

    stress = None
    if n_trials:
        g, z, gpow, zp, Tab, At, Bt = instance()
        r1 = walk_pass(4, 3000, False, g, z, Tab, At, Bt)
        r2 = walk_pass(4, 3000, False, g, z, Tab, At, Bt, stuck=2)
        mcap = int(40 * cavg)
        r3 = walk_pass(cap_walk, 10 ** 9, False, g, z, Tab, At, Bt, meter_cap=mcap)
        mp = match_and_verify(r1["log"], 1, g, z, gpow, zp)
        synth = [(7, 12345 % rp, 67890 % rp), (7, 999 % rp, 111 % rp)]
        m1 = match_and_verify(synth, 1, g, z, gpow, zp); m0 = match_and_verify(synth, 0, g, z, gpow, zp)
        stress = {"cap4": {"reports_logged": len(r1["log"]), "log_cap": walks * 4, "cap_stops": r1["events"]["cap_stops"],
                           "log_equals_slot_report_sum": len(r1["log"]) == r1["sum_slot_reports"]},
                  "stuck2": {"restarts": r2["events"]["abandoned"], "reports_logged": len(r2["log"]),
                             "log_within_cap": len(r2["log"]) <= walks * 4,
                             "log_equals_slot_report_sum": len(r2["log"]) == r2["sum_slot_reports"]},
                  "meter": {"slot_meter_cap": mcap, "slots_stopped_by_meter": r3["events"]["meter_stops"],
                            "max_slot_meter": max(w["meter"] for w in r3["slots"]),
                            "bound_cap_minus_1_plus_I_MUL": mcap - 1 + I_MUL,
                            "every_slot_within_bound": all(w["meter"] <= mcap - 1 + I_MUL for w in r3["slots"])},
                  "pairs": {"real_log_budget_1_pairs": mp["pairs"], "real_log_budget_1_aborted": mp["cap_hit"],
                            "two_record_budget_1": [m1["pairs"], m1["cap_hit"]], "two_record_budget_0": [m0["pairs"], m0["cap_hit"]]}}
    nt = max(1, n_trials)
    sq = math.sqrt(N); sqp = math.sqrt(Np)
    f39 = sum(t <= b39 for t in Tc_in) / nt
    f103 = n_coll_in / nt
    frec = n_rec / nt
    mean_full = sum(Tc_full) / len(Tc_full) / sqp if Tc_full else None
    sd_full = (math.sqrt(sum((t / sqp - mean_full) ** 2 for t in Tc_full) / max(1, len(Tc_full) - 1)) if Tc_full else None)
    rms_dev = math.sqrt(sum(d * d for d in dev_norm) / len(dev_norm)) if dev_norm else None
    totst = max(1, steps_tot)
    return {"q": q, "r_prime": rp, "classes_N": N, "R": R, "dp_bits": dp_bits, "bits": toy_bits, "s": s,
            "N_prime_over_N": round(1 / coll, 5), "trials_completed": n_trials,
            "step_budget_1.03_S0p": round(budget_steps, 1),
            "fraction_Tc_le_S0p": round(f39, 3), "fraction_Tc_le_S0p_err": round(binom_err(f39, nt), 3),
            "fraction_Tc_le_1.03_S0p": round(f103, 3), "fraction_Tc_le_1.03_S0p_err": round(binom_err(f103, nt), 3),
            "recovered_and_verified_fraction": round(frec, 3), "recovered_and_verified_count": n_rec,
            "trials_with_first_collision_in_budget": n_coll_in, "collision_in_budget_without_recovery": n_coll_no_rec,
            "mean_Tc_over_sqrtNp_uncensored": round(mean_full, 4) if mean_full else None,
            "mean_Tc_standard_error": round(sd_full / math.sqrt(len(Tc_full)), 4) if Tc_full and len(Tc_full) > 1 else None,
            "trials_with_Tc_before_horizon": len(Tc_full),
            "model_mean_sqrt_pi_over_2": round(math.sqrt(math.pi / 2), 4),
            "sharp_model_mean": round(math.sqrt(math.pi / 2) * math.sqrt(coll / (1 - s * s - (1 - s) ** 2 / (12 * R))), 4),
            "squaring_fraction_measured": round(nsq_tot / max(1, nsq_tot + nmul_tot), 4),
            "two_cycles_per_step": stats["two_cycles"] / totst,
            "predicted_two_cycles_per_step_(1-s)^2/(12R)": (1 - s) ** 2 / (12 * R),
            "long_cycle_events": stats["long_cycles"], "long_cycles_per_step": stats["long_cycles"] / totst,
            "squarings_inside_cycle_traversals": stats["squarings_in_cycle_traversal"],
            "abandoned_walks": stats["abandoned"], "key_ties_toy_only": tf.ties,
            "distinguished_fraction_measured": round(reports_tot / totst, 5), "distinguished_fraction_nominal": theta,
            "report_cap_per_slot_toy": cap_walk, "report_log_cap_toy": log_cap_total,
            "report_log_max_length": max_log, "report_cap_stops": stats["cap_stops"],
            "duplicate_appends_skipped": stats["duplicate_appends_skipped"],
            "pairs_charged_max_per_trial": pairs_max, "pair_budget_toy": pair_budget_toy, "pair_budget_aborts": n_cap_hit,
            "detection_window_toy_steps": win_extra,
            "slot_meter_cap_toy": ceil(1.01 * cavg * win_end / walks),
            "slots_stopped_by_meter_in_trials": stats["meter_stops"],
            "max_slot_meter_in_trials": max((w["meter"] for rw in [1] for w in []), default=None),
            "meter_rms_relative_deviation_x_sqrt_steps": round(rms_dev, 4) if rms_dev else None,
            "meter_rms_independent_steps_prediction": round(sig_ind, 4), "meter_slots_measured": len(dev_norm),
            "stress": stress}


# ---------- run the parts ----------
out = {"partA_production_BLS12_381": partA}
if "Q" in parts:
    out["partA2_bit_uniformity"] = partA2(min(1.2, max(0.3, 0.12 * budget)))
if "P" in parts:
    out["partP_production_field_walk"] = partP(min(1.0, max(0.3, 0.1 * budget)))
if "E" in parts:
    pq = cfg.get("pair_q", 67); pR = cfg.get("pair_R", 64)
    pbits = cfg.get("pair_bits", "hash" if pq < 2 ** 16 else "t0")
    e_budget = budget - (time.time() - T_START) if parts.strip("E") == "" else min(max(1.2, 0.08 * budget), 45.0)
    out["partE_pair_collision_test"] = pair_test(pq, pR, pair_samples, pbits, time.time() + max(0.2, e_budget),
                                                 random.Random(seed + 11))
if "B" in parts:
    TB0 = time.time()
    budgetB = max(1.0, budget - (TB0 - T_START) - 0.3)
    rows = []
    cfgs = CONFIGS if fields == "all" else [(cfg["q"], cfg["R"], cfg["dp_bits"])]
    slices = SLICES if fields == "all" else [1.0]
    for i, (q, R, dpb) in enumerate(cfgs):
        stop = TB0 + budgetB * sum(slices[:i + 1])
        rows.append(run_toy(q, R, dpb, trials, stop, random.Random(seed * 1000 + i)))
    out["partB_toy"] = {"walks": walks, "bits": toy_bits, "s": float(S_FRAC),
                        "branch_bits": BRANCH_BITS, "branch_threshold": BRANCH_THRESHOLD,
                        "trials_stop_at_first_recovery": False,
                        "every_recovered_logarithm_verified": all(rw["recovered_and_verified_count"] >= 0 for rw in rows),
                        "rows": rows}

# ---------- Part C: exact bound arithmetic used in proof.md (no randomness) ----------
def red12(f):
    v = [0] * 12; v[f] = 1
    for e in range(11, 3, -1):
        c = v[e]
        if c: v[e] = 0; v[e - 2] += c; v[e - 4] -= c
    return tuple(v[:4])
roots = [red12(f) for f in range(12)]
dist = {(0, 0, 0, 0): 1}; zero_sum = {}
for k in range(1, 11):
    nd = {}
    for v, c in dist.items():
        for s_ in roots:
            w_ = (v[0] + s_[0], v[1] + s_[1], v[2] + s_[2], v[3] + s_[3]); nd[w_] = nd.get(w_, 0) + c
    dist = nd; zero_sum[k] = dist.get((0, 0, 0, 0), 0) / 12 ** k
def log2_PL(L, RR, mfrac=1.0):
    ts = [(math.lgamma(RR + 1) - math.lgamma(d + 1) - math.lgamma(RR - d + 1) + L * math.log(d / RR)
           - d * math.log(12)) / math.log(2) for d in range(1, min(RR, L // 2) + 1)]
    mx = max(ts); return mx + math.log2(sum(2 ** (t - mx) for t in ts)) + L * math.log2(mfrac)
def log2sum(xs):
    mx = max(xs); return mx + math.log2(sum(2 ** (x - mx) for x in xs))
mf = float(1 - S_FRAC)
cyc = {}
for RR in (64, 256, 1024):
    cyc[RR] = {"plain_rate_L3_to_64_times_R2": up(2 ** log2sum([log2_PL(L, RR) for L in range(3, 65)]) * RR * RR, 4),
               "mixed_rate_L3_to_64_times_R2_(1-s)^L": up(2 ** log2sum([log2_PL(L, RR, mf) for L in range(3, 65)]) * RR * RR, 4),
               "mixed_log2_rate_L65_to_2R": up(log2sum([log2_PL(L, RR, mf) for L in range(65, 2 * RR + 1)]), 1)}
# per-step ledgers (instructions; blst_fp_mul = 383 measured; add_mod_384 = 45 and sub_mod_384 = 44 counted
# line by line: 18 instructions of frame, PAC hints and modulus loads, then a 27 / 26 instruction body)
FP_MUL, FP_ADD, FP_SUB = 383, 45, 44
COMMON = {"T0_low_limbs_and_three_bit_fields": 12, "two_cycle_T0_compare": 18, "checkpoint_T0_compare": 15,
          "T0_history_ring_and_64_step_checkpoint_test": 8, "slot_meter_add_and_cap_test_128bit": 7,
          "steps_since_report_counter_and_stuck_test": 7, "slot_step_counter_128bit": 4,
          "loop_control_output_buffer_slot_pointer": 30}
MULT_ONLY = {"load_w1_coefficient": 16, "4_fp_mul": 4 * FP_MUL,
             "B_add_mod_p": FP_ADD, "two_sub_mod_p": 2 * FP_SUB, "min3_max3_two_triples": 144,
             "two_reflections_p_minus_max": 40, "final_min_of_4_and_kappa": 74,
             "t_and_table_address": 10, "alpha_beta_add_mod_r": 52}
gs_mul = gs_ops["fp_mul"]; gs_add = gs_ops["fp_add"]; gs_sub = gs_ops["fp_sub"]
# The squaring itself is MEASURED, not counted from a formula. CYCLOTOMIC_SQR_MEASURED = 9,798 is the
# Callgrind executed-instruction count of blst f7a88067 `blst_fp12_cyclotomic_sqr` per call, measured
# with the track's COST_MODEL toolchain (Debian cc 14.2.0, valgrind-3.24.0, -O3 -fno-tree-vectorize
# -mgeneral-regs-only), identical at five slice offsets over five different GT elements. The routine is
# straight-line (18 functions, 1,017 static instructions, no conditional, backward or indirect branch),
# so its count cannot depend on the input, and it equals blst_fp12_sqr on GT elements, where the
# element lies in the cyclotomic subgroup and cyclotomic_sqr_fp12 is the legitimate squaring. The same
# harness reproduces blst_fp_mul = 383.0 and blst_fp12_mul = 22,996.0 exactly, the two prices this
# ledger is built on. Counting the Granger-Scott formula line by line through Fp primitives instead
# (18 Fp multiplies, 64 Fp adds, 31 Fp subs, priced as standalone blst_fp_mul calls) gives 11,190, i.e.
# it over-charges the real routine by 1,392: the shipped squaring uses a double-width sqr_382x, not 18
# separate blst_fp_mul calls. The measurement is what is charged.
CYCLOTOMIC_SQR_MEASURED = 9798
GS_FORMULA_COUNT = gs_mul * FP_MUL + gs_add * FP_ADD + gs_sub * FP_SUB
SQ_ONLY = {"cyclotomic_squaring_measured": CYCLOTOMIC_SQR_MEASURED, "alpha_beta_double_mod_r": 52}
common = sum(COMMON.values()); mult_only = sum(MULT_ONLY.values()); sq_only = sum(SQ_ONLY.values())
MUL_STEP_COUNT, SQ_STEP_COUNT = common + mult_only, common + sq_only
MUL_BUDGET, SQ_BUDGET = I_MUL - C, I_SQ
assert MUL_STEP_COUNT <= MUL_BUDGET and SQ_STEP_COUNT <= SQ_BUDGET
c_mul = Fraction(I_MUL, C); c_sq = Fraction(I_SQ, C)
c_avg = (1 - S_FRAC) * c_mul + S_FRAC * c_sq
# work meter and the walk-phase work bound
OVER = Fraction(3, 12 * R_PROD) + Fraction(200, R_PROD ** 2)
EXTRA_SLOT = math.ceil(OVER * S_SLOT)
NOMINAL_SLOT = c_avg * C * (S_SLOT + E_SLOT) + I_MUL * EXTRA_SLOT
DELTA = Fraction(1, 100)
W_CAP = math.ceil((1 + DELTA) * NOMINAL_SLOT)
SLOT_METER_MAX = W_CAP - 1 + I_MUL
W_WALK = M_WALKS * SLOT_METER_MAX
SLOT_STEPS_MAX = SLOT_METER_MAX // I_SQ + 1
RESTARTS_SLOT_MAX = SLOT_STEPS_MAX // (40 * 2 ** THETA_LG)
CREP_LG = 20 + RW_LG
APPEND, SORT, DEDUP, INIT = 3000, 10000, 100, 12
levels = 67
extract_levels = sum(j * 2 ** j for j in range(1, CREP_LG)) + CREP_LG
assert levels * 130 <= SORT and extract_levels + 2 ** CREP_LG <= levels * 2 ** CREP_LG
PER_RECORD = APPEND + SORT + DEDUP
PAIR_LG = 41
ANC = {"report_handling": 2 ** CREP_LG * PER_RECORD, "log_initialisation": 2 ** CREP_LG * INIT,
       "pairs": 2 ** PAIR_LG * 200, "element_recomputations": 2 * 2 ** PAIR_LG * 1021 * C,
       "setup_output_verification": 2 ** 31 * C,
       "restart_starts": M_WALKS * RESTARTS_SLOT_MAX * (1021 * C + 100)}
ANC_TOTAL = sum(ANC.values())
W_INSTR = W_WALK + ANC_TOTAL
# The preprocessing/advice bound of proof.md sections 5 and 6 must cover the log initialization, the
# multiplier table and the setup/output term, all of which are inside W already.
PREPROC_UNITS = 2 ** 56.3 + 2 ** 20.02 + 2 ** 31
W_NO_INIT = W_INSTR - ANC["log_initialisation"]
lgW = math.log2(W_INSTR) - math.log2(C)
CLAIM_NUM = 125585
claim_ok = W_INSTR ** 1000 <= C ** 1000 * 2 ** CLAIM_NUM
claim_tight = W_INSTR ** 1000 > C ** 1000 * 2 ** (CLAIM_NUM - 1)
lgS = math.log2(S_STEPS)
model = 1 - 0.61 ** (1.03 ** 2)
lgS_target = lgS0 + math.log2(1.03)
# Work-meter concentration: a slot's charged work over its S_slot + E_slot steps is a sum of step costs in
# [c_sq, c_mul] with mean c_avg, so Hoeffding bounds the event that it exceeds 1.01 of its nominal work.
C_MULF, C_SQF, C_AVGF = float(c_mul), float(c_sq), float(c_avg)
HOEFF_K = 2 * (0.01 * C_AVGF) ** 2 / (C_MULF - C_SQF) ** 2          # nats per step
HOEFF_SLOT = HOEFF_K * (S_SLOT + E_SLOT) / math.log(2)              # P <= 2^-HOEFF_SLOT per slot
HOEFF_RUN = HOEFF_K * S_STEPS / math.log(2)                         # same exponent over the run's S steps
losses = {"self_collision": 2 ** -19, "restart_waste": 2 ** -21, "detection_tail": 2 ** -28,
          "pair_budget_abort": 2 ** -36, "work_meter_concentration": 2.0 ** -100,
          "trailing_walk_restart": 2 ** -52, "cycle_steps_and_ties": 2 * 2 ** -37,
          "long_fruitless_cycles": 2 ** -100, "T0_false_cycle_match": 2 ** -253,
          "report_cap_on_healthy_slot": 2 ** -300}
loss_sum = sum(losses.values())
ALLOW = 0.01
mu = mean_reports_slot; a_cap = CAP_RATIO
kl = a_cap * math.log(a_cap) - a_cap + 1
mem_bytes = 2 ** CREP_LG * 80 + 2 ** 20 * 1024 + 12288 * 576 + 2 * 12288 * 32
counters = {"slot_meter": {"max_log2": up(math.log2(SLOT_METER_MAX), 6), "words": 2,
                           "fits": SLOT_METER_MAX < 2 ** 128, "cap_test_below_capacity": W_CAP < 2 ** 127},
            "slot_step_counter": {"max_log2": up(math.log2(SLOT_STEPS_MAX), 6), "words": 2, "fits": SLOT_STEPS_MAX < 2 ** 128},
            "steps_since_report": {"max_log2": up(math.log2(40 * 2 ** THETA_LG + 1), 6), "words": 2,
                                   "fits": 40 * 2 ** THETA_LG + 1 < 2 ** 128},
            "slot_report_counter": {"max_log2": RW_LG, "words": 1, "fits": 2 ** RW_LG < 2 ** 64,
                                    "cap_test_below_capacity": 2 ** RW_LG < 2 ** 63},
            "log_address_80_x_(2^47_slot_plus_r)": {"max_log2": up(math.log2(80 * 2 ** CREP_LG), 6), "words": 2,
                                                    "fits": 80 * 2 ** CREP_LG < 2 ** 128},
            "slot_restart_counter": {"max_log2": up(math.log2(max(1, RESTARTS_SLOT_MAX)), 6), "words": 1,
                                     "fits": RESTARTS_SLOT_MAX < 2 ** 64},
            "pair_counter": {"max_log2": PAIR_LG, "words": 1, "fits": 2 ** PAIR_LG + 1 < 2 ** 64},
            "heap_index": {"max_log2": CREP_LG, "words": 2, "fits": 2 ** CREP_LG < 2 ** 128}}
partC = {"zero_sum_prob_k_unit_12th_roots": {k: round(zero_sum[k], 6) for k in range(2, 11)},
         "zero_sum_max_k_ge_2": round(max(zero_sum[k] for k in range(2, 11)), 6),
         "fruitless_cycle_bounds": cyc,
         "log2_N": round(lgN, 4), "collision_factor_1-s^2-(1-s)^2/R": float(COLL),
         "log2_N_prime": up(lgNp, 6), "log2_S0_prime": up(lgS0, 6),
         "per_slot_step_budget": S_SLOT, "log2_S": up(lgS, 6), "S_ge_1.03_S0_prime": lgS >= lgS_target,
         "model_success_at_1.03_S0_prime": dn(model, 5), "H1plus_allowance": ALLOW,
         "loss_terms": losses, "loss_sum_below_2^-18": loss_sum < 2 ** -18,
        "work_meter_concentration": {
            "hoeffding_nats_per_step": round(HOEFF_K, 12),
            "step_cost_range": [C_SQF, C_MULF], "mean_step_cost": C_AVGF,
            "slot_bound_log2_P_is_minus_2_to_the": up(math.log2(HOEFF_SLOT), 2),
            "slot_bound_far_below_2^-100": HOEFF_SLOT > 100,
            "run_bound_log2_P_is_minus_2_to_the": up(math.log2(HOEFF_RUN), 2)},
        "record_capacity_is_deterministic": {
            "per_slot_report_cap_log2": RW_LG, "log_records_log2": CREP_LG,
            "capacity_without_a_return_assumption": M_WALKS * 2 ** RW_LG == 2 ** CREP_LG},
        "T0_false_match_comparisons_log2": up(math.log2(2 * M_WALKS * SLOT_STEPS_MAX), 4),
        "T0_false_match_total_log2": up(math.log2(2 * M_WALKS * SLOT_STEPS_MAX / p), 2),
        "T0_false_match_total_below_2^-253": 2 * M_WALKS * SLOT_STEPS_MAX / p < 2 ** -253,
         "success_lower_bound": dn(model - ALLOW - 2 ** -18, 4),
         "success_margin_over_0.39": dn(model - ALLOW - 2 ** -18 - 0.39, 4),
         "ledger_common": COMMON, "ledger_multiplication_only": MULT_ONLY, "ledger_squaring_only": SQ_ONLY,
         "common_count": common, "multiplication_only_count": mult_only, "squaring_only_count": sq_only,
         "multiplication_step_count_and_budget": [MUL_STEP_COUNT, MUL_BUDGET],
         "squaring_step_count_and_budget": [SQ_STEP_COUNT, SQ_BUDGET],
         "cyclotomic_squaring_measured_instructions": CYCLOTOMIC_SQR_MEASURED,
         "squaring_step_slack_over_count": SQ_BUDGET - SQ_STEP_COUNT,
         "gs_formula_line_by_line_count_overcharges_measured": GS_FORMULA_COUNT - CYCLOTOMIC_SQR_MEASURED,
         "branch_bits_and_threshold": [BRANCH_BITS, BRANCH_THRESHOLD],
         "branch_fraction": str(S_FRAC),
         "c_mul": float(c_mul), "c_sq": float(c_sq), "c_avg": float(c_avg),
         "cycle_overhead_fraction_charged_at_c_mul": float(OVER),
         "nominal_over_plain_step_work": up(float(NOMINAL_SLOT / (c_avg * C * (S_SLOT + E_SLOT))), 7),
         "meter_delta": float(DELTA), "slot_work_cap_instructions_log2": up(math.log2(W_CAP), 6),
         "slot_work_cap_instructions": str(W_CAP), "nominal_slot_instructions_log2": up(math.log2(NOMINAL_SLOT), 6),
         "slot_meter_max_instructions_log2": up(math.log2(SLOT_METER_MAX), 6),
         "slot_step_counter_max_log2": up(math.log2(SLOT_STEPS_MAX), 6),
         "walk_work_bound_instructions_log2": up(math.log2(W_WALK), 6),
         "ancillary_total_instructions_log2": up(math.log2(ANC_TOTAL), 6),
         "report_handling_units_log2": up(math.log2(2 ** CREP_LG * PER_RECORD) - math.log2(C), 4),
         "branch_fraction": [BRANCH_THRESHOLD, 1 << BRANCH_BITS], "branch_fraction_value": float(S_FRAC),
         "walk_work_bound_log2_units": up(math.log2(W_WALK) - math.log2(C), 6),
         "ancillary_log2_units": {k: up(math.log2(v) - math.log2(C), 2) for k, v in ANC.items()},
         "ancillary_total_log2_units": up(math.log2(ANC_TOTAL) - math.log2(C), 3),
         "counters": counters, "all_counters_fit": all(v["fits"] for v in counters.values()),
         "report_cap_R_w_log2": RW_LG, "mean_reports_per_slot_log2": up(math.log2(mu), 4),
         "cap_over_mean_a": dn(a_cap, 4), "kl_D_a_1": dn(kl, 4),
         "log2_cap_fire_on_healthy_slot": up(-mu * kl / math.log(2), 1),
         "per_record_instructions": {"append_incl_key_computation": APPEND, "heapsort": SORT, "dedup": DEDUP,
                                     "total": PER_RECORD},
         "memory_log2_bytes": up(math.log2(mem_bytes), 4), "memory_below_2^73.4": mem_bytes < 2 ** 73.4,
         "preprocessing_and_table_construction_units_log2": up(math.log2(PREPROC_UNITS), 4),
         "preprocessing_within_2^57_units": PREPROC_UNITS <= 2 ** 57,
         "log_initialisation_units_log2": up(math.log2(ANC["log_initialisation"]) - math.log2(C), 4),
         "log_initialisation_is_inside_W": W_INSTR == W_WALK + ANC_TOTAL and ANC["log_initialisation"] == 2 ** CREP_LG * INIT,
         "claim_holds_with_the_initialisation_term_in_W": W_INSTR ** 1000 <= C ** 1000 * 2 ** CLAIM_NUM,
         "claim_also_holds_if_the_initialisation_term_were_left_out": W_NO_INIT ** 1000 <= C ** 1000 * 2 ** CLAIM_NUM,
         "initialisation_term_raises_log2_W_by_less_than_2^-53": W_INSTR * 2 ** 53 <= W_NO_INIT * (2 ** 53 + 1),
         "initialisation_term_instructions_log2": up(math.log2(W_INSTR - W_NO_INIT), 3),
         "log2_W": up(lgW, 6), "claim": CLAIM_NUM / 1000,
         "claim_exact_integer_check_W^1000_le_C^1000_2^claim": claim_ok,
         "claim_minus_0.001_fails": claim_tight}
if "C" in parts:
    out["partC_bound_arithmetic"] = partC
out["seconds"] = round(time.time() - T_START, 2)
out["scope"] = ("Part A and A2 check exact identities and bit frequencies on real BLS12-381 GT elements; Part P runs the "
                "production step rule on real GT elements with a 16-entry table, including the different-histories "
                "coalescence test on a constructed 2-cycle; Parts E and B are toy fields, so they "
                "do not validate production cost bounds or production bit statistics.")
print(json.dumps(out, indent=1, default=str))
