"""Fixed-G1 pairing inversion into G2 by Frobenius-class rho with a cyclotomic-squaring branch: identities, ledger, toy runs.

Part A  exact identities on BLS12-381 itself (key formula, canonical conjugate,
        Granger-Scott squaring with operation counts, orders of p and 2 modulo r).
Part C  the production step rule (twisted-table multiply, squaring branch, branch and
        distinguished-point bits taken from stored T) executed in the real BLS12-381 GT for
        a few thousand steps, checked against explicit Frobenius arithmetic.
Part B  the complete algorithm (cycle escapes, distinguished points, collision
        solving) in a small field of the same tower shape, with first-collision
        statistics compared with the random-mapping law adjusted for the branch mix.
Part F  the instance reduction and output on BLS12-381: g = e(P,Q0) by a transcription of the pinned
        pairing(Q0,P), Q = [x]Q0 in G2 with its counted cost, the verification e(P,Q) = z and the
        identity target z = 1 (Q = infinity).
Part D  exact recomputation of the work ledger used in proof.md, with part F's counts.
Toy runs are evidence at their own scale only; they do not validate production costs.
"""
import json, math, random, sys, time
from fractions import Fraction

T0 = time.time()
cfg = json.load(sys.stdin)
Q_TOY = cfg.get("q", 67)
TRIALS = cfg.get("trials", 20)
R = cfg.get("R", 64)
SQ64 = cfg.get("sq64", 35)
DP_BITS = cfg.get("dp_bits", 3)
WALKS = cfg.get("walks", 4)
SEED = cfg.get("seed", 1)
PROD_CHECKS = cfg.get("production_checks", 2)
PROD_STEPS = cfg.get("production_steps", 3000)
BUDGET = cfg.get("seconds", 12)
PAIR_SAMPLES = cfg.get("pair_samples", 0)
rng = random.Random(SEED)

M64 = (1 << 64) - 1

# ---------------------------------------------------------------- price vector (proof.md section 5)
# One dense Fp12 multiplication is C = 20,000 instructions.  The per-step bookkeeping and key
# charges are analytic listings (proof.md section 5); the squaring branch is the itemisation of
# measured-kernel.md: the kernel's fixed operation counts times the static instruction length of
# each blst f7a88067 arm64 primitive body, plus the disassembled compiled glue body and the loop,
# cross-checked by Callgrind on the same compiled kernel.
CHARGE_C = 20000
FPMUL, FPADD = 383, 45                       # blst_fp_mul; modular add/sub/neg as a called routine
FP2SQR = 2 * FPMUL + 3 * FPADD               # complex squaring: 2 Fp mul + 3 Fp add/sub
FP2ADD = 2 * FPADD
GS_LISTING = 9 * FP2SQR + 34 * FP2ADD + 200  # old narrow 9 + 34 + 200 listing, kept for reference
COMMON_ITEMS = {"read Re(c_0), branch bits, DP test, table address": 50,
                "(alpha, beta) update mod r (add or doubling)": 100,
                "2-cycle and checkpoint comparisons on T, history": 60,
                "loop control": 30,
                "margin for loads and moves": 30}
MULT_ITEMS = {"4 Fp multiplications (a*zeta, a*d0, b*d1, B*zeta)": 4 * FPMUL,
              "1 Fp addition (B) + 8 Fp sub/neg": 9 * FPADD,
              "minimum of 12 keys (top limbs, tie test)": 150,
              "load c_1, kappa, twist index": 32}
COMMON = sum(COMMON_ITEMS.values())
MULT_EXTRA = sum(MULT_ITEMS.values())
# Static instruction lengths of the pinned blst f7a88067 arm64 primitive bodies (one line per
# encoded instruction, prologue/epilogue, call and return included, alignment nops excluded) and of
# the compiled glue body sqmeasure_kernel; measured-kernel.md lists each symbol.  Every line below is
# operation count x those lengths; the items sum to 9,708.  Callgrind measures the same kernel at
# 79,527,954 Ir / 8,192 calls = 9,708.0022 per call (the items plus 18 instructions of harness
# entry); the charge is the ceiling of the larger figure, 9,709.
GS_IR, GS_CALLS = 79527954, 8192             # Callgrind record, measured-kernel.md
SQ_ITEMS = {
    "9 wide Fp2 squarings (sqr_382x 86 + 2 x __mul_384 162)": 9 * (86 + 2 * 162),
    "24 wide 768-bit Fp add/sub, 8 per Fp4 squaring (9 add x 62, 15 sub x 61)":
        9 * (17 + 45) + 15 * (17 + 44),
    "12 Montgomery reductions (redc 23 + mul_by_1 154 + tail 27)": 12 * (23 + 154 + 27),
    "18 narrow Fp2 adds (add_mod_384x 24 + 2 x __add_mod_384 27)": 18 * (24 + 2 * 27),
    "3 narrow Fp2 subs (sub_mod_384x 24 + 2 x __sub_mod_384 26)": 3 * (24 + 2 * 26),
    "1 xi multiply (wrapper 29 + __sub_mod_384 26 + __add_mod_384_ab_are_loaded 21)": 29 + 26 + 21,
    "sqmeasure_kernel body, straight-line (270 mov/add, 67 bl, 28 ldr/str, 12 ldp/stp, sub, ret)":
        270 + 67 + 28 + 12 + 1 + 1,
    "loop control": 10,
}
SQ_SUM = sum(SQ_ITEMS.values())
SQ_PRIMITIVES = SQ_SUM - 379 - 10            # the six primitive lines, without glue body and loop
assert SQ_SUM == 9708 and SQ_PRIMITIVES == 9319
GS_CHARGE = max(SQ_SUM, -(-GS_IR // GS_CALLS))
assert GS_CHARGE == 9709 and GS_CHARGE * GS_CALLS >= GS_IR and GS_CHARGE >= SQ_SUM
MUL_STEP = 1 + (COMMON + MULT_EXTRA) / CHARGE_C
SQ_STEP = (GS_CHARGE + COMMON) / CHARGE_C


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
        if (u & 63) < SQ64:
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
            h = (u >> 6) % rprod
            bins_h[(u >> 6) & 15] += 1
            idx = (12 - kap) % 12
            y2 = P.mul(y, Tab[h][idx])
            if t < verify:
                bad_class += P.frobk(y2, kap) != P.mul(P.frobk(y, kap), Mh[h])
            y = y2
            al = (al + ah[h] * lam[idx]) % r
            be = (be + bh[h] * lam[idx]) % r
            n_mul += 1
    exp_ok = y == P.mul(P.pow(g, al), P.pow(z0, be))
    n_br = sum(bins_h)
    n_all = sum(bins_sq)
    chi_sq = sum((c - n_all / 8) ** 2 / (n_all / 8) for c in bins_sq)
    chi_h = sum((c - n_br / 16) ** 2 / (n_br / 16) for c in bins_h) if n_br else None
    return {"steps": steps, "R_table": rprod, "squaring_steps": n_sq, "multiply_steps": n_mul,
            "key_rule": "stored Montgomery limbs (x*2^384 mod p)" if stored else "ordinary integers in [0,p)",
            "min_index_differs_from_integer_rule": idx_differs,
            "squaring_fraction": round(n_sq / max(1, n_all), 4), "target_fraction": SQ64 / 64,
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


def run_toy(T, q, rp, rng, walks, Rm, sq64, dp_bits, cap_factor=40):
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
        if (u & 63) < sq64:
            return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp, kap, can, K, u
        h = (u >> 6) % Rm
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
            if (u >> 16) & mask == 0:
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
            stats["sq"] += (u & 63) < sq64
            w["y"], w["a"], w["b"], _, _, _, _ = advance(y, a, b)
            w["n"] += 1
            w["since"] += 1
            steps += 1
            if w["since"] > cap_factor << dp_bits:
                stats["restarts"] += 1
                w.update(new_walk())
    return Tc, steps, stats


def part_B(rng, q, trials, Rm, sq64, dp_bits, walks, budget):
    T = Tower(q)
    rp = toy_setup(q)
    if not (pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))):
        raise RuntimeError('toy field does not have embedding degree 12')
    N = (rp - 1) // 12
    s = sq64 / 64
    sumq2 = s * s + (1 - s) ** 2 / Rm
    Np = N / (1 - sumq2)
    Tcs, tot, ev, dev = [], [], {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "useless": 0}, []
    t_start = time.time()
    while len(Tcs) < trials and time.time() - t_start < budget:
        Tc, steps, st = run_toy(T, q, rp, rng, walks, Rm, sq64, dp_bits)
        Tcs.append(Tc)
        tot.append(steps)
        for k in ev:
            ev[k] += st[k]
        sf = sq64 / 64
        # the toy's charge bookkeeping uses this package's price vector (section 5); the ratio
        # charge/nominal is invariant to which vector is used, both sides take the same one
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
    return {"q": q, "r_prime": rp, "classes_N": N, "R": Rm, "sq64": sq64, "squaring_prob": s,
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
            "charge_fluctuation_iid_prediction": round(math.sqrt(s * (1 - s)) * 0.6205 / dev[0][2], 4) if dev else None,
            "events": ev, "two_cycles_per_step": ev["two_cycles"] / allsteps,
            "Tc_over_sqrtN": [round(t / sN, 4) for t in Tcs] if cfg.get("dump") else None}


# ---------------------------------------------------------------- pair-collision law (Part E)
def part_E(rng, q, m, Rm, sq64):
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
        if (u & 63) < sq64:
            y2 = T.cyc_sqr(y)
        else:
            y2 = T.mul(y, Tab[(u >> 6) % Rm][(12 - kap) % 12])
        _, can2, _ = classinfo(T, y2)
        images[can2] = images.get(can2, 0) + 1
    pairs = sum(n * (n - 1) // 2 for n in images.values())
    tot = m * (m - 1) // 2
    s = sq64 / 64
    obs = pairs / tot * N
    return {"q": q, "classes_N": N, "sampled_classes": m, "R": Rm, "squaring_prob": s,
            "image_coincidence_pairs": pairs,
            "N_times_pair_collision_rate_observed": round(obs, 4),
            "relative_standard_error": round(1 / math.sqrt(max(1, pairs)), 4),
            "predicted_exact_1_minus_s2_minus_(1-s)^2/(12R)": round(1 - s * s - (1 - s) ** 2 / (12 * Rm), 4),
            "predicted_conservative_1_minus_s2_minus_(1-s)^2/R": round(1 - s * s - (1 - s) ** 2 / Rm, 4),
            "random_mapping": 1.0}


# ---------------------------------------------------------------- report store (Part G)
def part_G(rng, records, groups):
    """The report store of section 3.2: append 128-bit fingerprints to a flat buffer, sort with
    eight stable 16-bit LSD radix passes, then scan adjacent records for equal fingerprints.
    Checks the fixed pass count, completeness against a dictionary reference and sorted order."""
    recs = []
    for gidx in range(groups):
        f = rng.getrandbits(128)
        for j in range(3):
            recs.append((f, gidx, j))
    while len(recs) < records:
        recs.append((rng.getrandbits(128), -1, len(recs)))
    rng.shuffle(recs)
    order = list(range(len(recs)))
    for pass_ in range(8):
        shift = 16 * pass_
        count = [0] * 65536
        for i in order:
            count[(recs[i][0] >> shift) & 0xFFFF] += 1
        start = 0
        for b in range(65536):
            c = count[b]
            count[b] = start
            start += c
        out = [0] * len(order)
        for i in order:
            d = (recs[i][0] >> shift) & 0xFFFF
            out[count[d]] = i
            count[d] += 1
        order = out
    found = set()
    for a, b in zip(order, order[1:]):
        if recs[a][0] == recs[b][0]:
            found.add(recs[a][0])
    ref = {}
    for rec in recs:
        ref.setdefault(rec[0], []).append(rec)
    expected = {f for f, lst in ref.items() if len(lst) > 1}
    return {"records": len(recs), "planted_duplicate_groups": groups, "fixed_passes": 8,
            "duplicate_groups_found": len(found), "duplicate_groups_expected": len(expected),
            "complete": found == expected,
            "sorted_ascending": all(recs[order[i]][0] <= recs[order[i + 1]][0] for i in range(len(order) - 1))}


# ---------------------------------------------------------------- reduction and output on BLS12-381 (Part F)
P0_XY = (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb,
         0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1)
Q0_XY = ((0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8,
          0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e),
         (0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801,
          0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be))
ATE_PBE = [0] * 16 + [1] + [0] * 31 + [1] + [0] * 8 + [1, 0, 0, 1, 0, 1, 1]   # the pinned pseudo_binary_encoding


class Fld:
    """One field of the pinned library's point and line formulas. Every call is one counted operation;
    the ledger charges each at one full unit (a dense Fp12 multiplication), whatever the field."""

    def __init__(self, kind):
        self.kind, self.n = kind, 0
        if kind == 1:
            self.zero, self.one = 0, 1
        elif kind == 2:
            self.zero, self.one = (0, 0), (1, 0)
        else:
            self.zero, self.one = (0,) * 12, (1,) + (0,) * 11

    def mul(self, a, b):
        self.n += 1
        if self.kind == 1:
            return a * b % p
        if self.kind == 2:
            return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)
        acc = [0] * 23
        for i in range(12):
            if a[i]:
                for j in range(12):
                    acc[i + j] += a[i] * b[j]
        for k in range(22, 11, -1):            # w^12 = 2 w^6 - 2 (the profile modulus)
            c = acc[k]
            acc[k - 6] += 2 * c
            acc[k - 12] -= 2 * c
        return tuple(v % p for v in acc[:12])

    def lin(self, a, b, sgn):
        self.n += 1
        if self.kind == 1:
            return (a + sgn * b) % p
        return tuple((u + sgn * v) % p for u, v in zip(a, b))

    def add(self, a, b):
        return self.lin(a, b, 1)

    def sub(self, a, b):
        return self.lin(a, b, -1)

    def smul(self, k, a):
        self.n += 1
        if self.kind == 1:
            return k * a % p
        return tuple(k * v % p for v in a)

    def eq(self, a, b):
        self.n += 1
        return a == b


def pt_double(F, pt):                          # pinned optimized_curve.double, operation for operation
    x, y, z = pt
    W = F.mul(F.smul(3, x), x)
    S = F.mul(y, z)
    B = F.mul(F.mul(x, y), S)
    H = F.sub(F.mul(W, W), F.smul(8, B))
    S2 = F.mul(S, S)
    nx = F.mul(F.smul(2, H), S)
    ny = F.sub(F.mul(W, F.sub(F.smul(4, B), H)), F.mul(F.mul(F.smul(8, y), y), S2))
    nz = F.mul(F.smul(8, S), S2)
    return (nx, ny, nz)


def pt_add(F, p1, p2):                         # pinned optimized_curve.add, operation for operation
    if F.eq(p1[2], F.zero) or F.eq(p2[2], F.zero):
        return p1 if p2[2] == F.zero else p2
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    U1, U2, V1, V2 = F.mul(y2, z1), F.mul(y1, z2), F.mul(x2, z1), F.mul(x1, z2)
    if F.eq(V1, V2) and F.eq(U1, U2):
        return pt_double(F, p1)
    elif F.eq(V1, V2):
        return (F.one, F.one, F.zero)
    U, V = F.sub(U1, U2), F.sub(V1, V2)
    V2s = F.mul(V, V)
    V2sV2 = F.mul(V2s, V2)
    V3 = F.mul(V, V2s)
    W = F.mul(z1, z2)
    A = F.sub(F.sub(F.mul(F.mul(U, U), W), V3), F.smul(2, V2sV2))
    return (F.mul(V, A), F.sub(F.mul(U, F.sub(V2sV2, A)), F.mul(V3, U2)), F.mul(V3, W))


def twist12(F, pt):                            # pinned optimized_curve.twist (coefficient placement)
    F.n += 1
    (x0, x1), (y0, y1), (z0, z1) = pt
    nx = [0] * 12; nx[1], nx[7] = (x0 - x1) % p, x1
    ny = [0] * 12; ny[0], ny[6] = (y0 - y1) % p, y1
    nz = [0] * 12; nz[3], nz[9] = (z0 - z1) % p, z1
    return (tuple(nx), tuple(ny), tuple(nz))


def linefunc(F, P1, P2, T):                    # pinned optimized_pairing.linefunc
    x1, y1, z1 = P1
    x2, y2, z2 = P2
    xt, yt, zt = T
    mn = F.sub(F.mul(y2, z1), F.mul(y1, z2))
    md = F.sub(F.mul(x2, z1), F.mul(x1, z2))
    if not F.eq(md, F.zero):
        pass
    elif F.eq(mn, F.zero):
        mn = F.mul(F.smul(3, x1), x1)
        md = F.mul(F.smul(2, y1), z1)
    else:
        return F.sub(F.mul(xt, z1), F.mul(x1, zt)), F.mul(z1, zt)
    return (F.sub(F.mul(mn, F.sub(F.mul(xt, z1), F.mul(x1, zt))), F.mul(md, F.sub(F.mul(yt, z1), F.mul(y1, zt)))),
            F.mul(F.mul(md, zt), z1))


def to_tower(c):                               # profile basis -> blst tower (u = w^6 - 1): a ring isomorphism
    return tuple(((c[m] + c[m + 6]) % p, c[m + 6]) for m in range(6))


def from_tower(y):
    return tuple((a - b) % p for a, b in y) + tuple(b for _, b in y)


def tpow(Tw, x, e, cnt):                       # square-and-multiply as Tower.pow, one unit per product
    out = Tw.one
    while e:
        if e & 1:
            out = Tw.mul(out, x); cnt[0] += 1
        x = Tw.mul(x, x); cnt[0] += 1
        e >>= 1
    return out


def pairing_pinned(Tw, Qpt, Ppt):
    """e(P,Q) as the pinned library defines it: pairing(Q,P) = miller_loop(Q,P) then f**((p^12-1)/r).
    Qpt is projective over Fp2, Ppt projective over Fp; None stands for infinity. Returns (value in the
    tower basis, counted units)."""
    if Qpt is None or Ppt is None:
        return Tw.one, 0
    F2, F12 = Fld(2), Fld(12)
    F12.n += 1                                 # cast_point_to_fq12
    cP = tuple((c,) + (0,) * 11 for c in Ppt)
    tR = tQ = twist12(F12, Qpt)
    R = Qpt
    fn, fd = F12.one, F12.one
    for v in ATE_PBE[62::-1]:
        n_, d_ = linefunc(F12, tR, tR, cP)
        fn = F12.mul(F12.mul(fn, fn), n_)
        fd = F12.mul(F12.mul(fd, fd), d_)
        R = pt_double(F2, R)
        tR = twist12(F12, R)
        if v == 1:
            n_, d_ = linefunc(F12, tR, tQ, cP)
            fn = F12.mul(fn, n_)
            fd = F12.mul(fd, d_)
            R = pt_add(F2, R, Qpt)
            tR = twist12(F12, R)
    cnt = [0]
    a, d = to_tower(fn), to_tower(fd)          # the change of basis, charged as two operations
    cnt[0] += 2
    conj = Tw.one                              # f_den^-1 = (prod of its 11 other conjugates) / norm, norm in Fp
    for k in range(1, 12):
        conj = Tw.mul(conj, Tw.frobk(d, k)); cnt[0] += 2
    nrm = Tw.mul(conj, d); cnt[0] += 1
    if any(nrm[m] != (0, 0) for m in range(1, 6)) or nrm[0][1] != 0:
        raise RuntimeError('norm of f_den is not in Fp')
    e, inv = p - 2, 1                          # Fermat inversion in Fp, one unit per Fp product
    base = nrm[0][0]
    while e:
        if e & 1:
            inv = inv * base % p; cnt[0] += 1
        base = base * base % p; cnt[0] += 1
        e >>= 1
    dinv = tuple(((u * inv) % p, (v * inv) % p) for u, v in conj); cnt[0] += 1
    f = Tw.mul(a, dinv); cnt[0] += 1
    g = tpow(Tw, f, (p ** 12 - 1) // r, cnt)
    return g, F2.n + F12.n + cnt[0]


def g2_mul_counted(x):
    """Q = [x]Q0 by left-to-right double-and-add with the pinned projective formulas over Fp2, then the
    affine encoding (one Fp2 inversion). Returns (affine Q or None for infinity, ladder units,
    encoding units, doublings, additions)."""
    F = Fld(2)
    Q0p = (Q0_XY[0], Q0_XY[1], (1, 0))
    if x == 0:
        return None, 0, 0, 0, 0
    Rp, nd, na = Q0p, 0, 0
    for bit in bin(x)[3:]:
        Rp = pt_double(F, Rp); nd += 1
        if bit == "1":
            Rp = pt_add(F, Rp, Q0p); na += 1
    ladder = F.n
    zx = Rp[2]
    if F.eq(zx, F.zero):
        return None, ladder, F.n - ladder, nd, na
    nrm = F.add(F.mul((zx[0], 0), (zx[0], 0)), F.mul((zx[1], 0), (zx[1], 0)))[0]
    e, inv, base = p - 2, 1, nrm
    while e:
        if e & 1:
            inv = inv * base % p; F.n += 1
        base = base * base % p; F.n += 1
        e >>= 1
    zi = F.mul((zx[0], (-zx[1]) % p), (inv, 0))
    return (F.mul(Rp[0], zi), F.mul(Rp[1], zi)), ladder, F.n - ladder, nd, na


def on_g2_curve(Qa):
    F = Fld(2)
    x, y = Qa
    return F.mul(y, y) == F.add(F.mul(F.mul(x, x), x), (4, 4))


def g2_branch_maxima():
    """Operation count of each branch of the pinned double/add on Fp2 points, so the hard maximum per
    ladder step is the largest branch, not the branch the sample happened to take."""
    F = Fld(2)
    Q0p = (Q0_XY[0], Q0_XY[1], (1, 0))
    negQ0 = (Q0_XY[0], (-Q0_XY[1][0] % p, -Q0_XY[1][1] % p), (1, 0))
    D = pt_double(F, Q0p)
    out = {}
    for name, a, b in (("generic", D, Q0p), ("equal_points", Q0p, Q0p), ("opposite_points", Q0p, negQ0),
                       ("infinity_operand", (F.one, F.one, F.zero), Q0p)):
        n0 = F.n; pt_add(F, a, b); out["add_" + name] = F.n - n0
    n0 = F.n; pt_double(F, D); out["double"] = F.n - n0
    return out


def part_F(rng):
    """Instance reduction and output on BLS12-381: g = e(P,Q0) per instance (pinned order pairing(Q0,P)),
    Q = [x]Q0 in G2 with its counted cost, verification e(P,Q) = z, and the identity target z = 1."""
    Tw = Tower(p)
    Fp1 = Fld(1)
    b = rng.randrange(1, r)
    Pp, acc = None, (P0_XY[0], P0_XY[1], 1)
    for bit in bin(b)[:1:-1]:                  # the challenge's own sampling of P = [b]P0 (not attack work)
        if bit == "1":
            Pp = acc if Pp is None else pt_add(Fp1, Pp, acc)
        acc = pt_double(Fp1, acc)
    zi = pow(Pp[2], p - 2, p)
    Pa = (Pp[0] * zi % p, Pp[1] * zi % p)
    P_ok = (Pa[1] * Pa[1] - Pa[0] ** 3 - 4) % p == 0
    g, units_g = pairing_pinned(Tw, (Q0_XY[0], Q0_XY[1], (1, 0)), (Pa[0], Pa[1], 1))
    g_ok = tpow(Tw, g, r, [0]) == Tw.one and g != Tw.one
    x = rng.randrange(1, r)                    # stands in for the logarithm the walk returns
    z = tpow(Tw, g, x, [0])
    Qa, lad_q, enc_q, nd, na = g2_mul_counted(x)
    Qr = g2_mul_counted(r)[0]                  # acceptance-side check [r]Q0 = O (not attack work)
    Qrx = g2_mul_counted(r * x)[0] if Qa is not None else None
    v, units_v = pairing_pinned(Tw, (Qa[0], Qa[1], (1, 0)), (Pa[0], Pa[1], 1))
    one_v, _ = pairing_pinned(Tw, None, (Pa[0], Pa[1], 1))
    br = g2_branch_maxima()
    per_add = max(v_ for k_, v_ in br.items() if k_.startswith("add_"))
    steps = (r - 1).bit_length() - 1           # x <= r-1 < 2^255: at most 254 doublings and 254 additions
    q_max = steps * (br["double"] + per_add) + enc_q
    return {
        "P_on_G1_curve_and_nonidentity": P_ok and Pp[2] != 0,
        "g_in_GT_and_not_1": g_ok,
        "pairing_units_counted_g": units_g,
        "pairing_units_counted_verification": units_v,
        "pairing_counts_equal": units_g == units_v,
        "final_exponent_bits": ((p ** 12 - 1) // r).bit_length(),
        "Q_on_G2_curve": on_g2_curve(Qa),
        "r_times_Q0_is_infinity": Qr is None,
        "r_times_Q_is_infinity": Qrx is None,
        "verification_e(P,Q)_equals_z": v == z,
        "g2_ladder_doublings_additions_this_x": [nd, na],
        "g2_ladder_units_this_x": lad_q,
        "g2_ops_per_double": br["double"],
        "g2_ops_per_add_by_branch": {k_[4:]: v_ for k_, v_ in br.items() if k_.startswith("add_")},
        "g2_affine_encoding_units": enc_q,
        "g2_ladder_steps_max": steps,
        "g2_output_units_hard_max_255_bit_scalar": q_max,
        "identity_target": {"z_equals_1_test": "before any walk", "x": 0, "Q": "infinity",
                            "e(P,O)_equals_1": one_v == Tw.one},
        "_g_profile": from_tower(g), "_P": Pa,
    }


# ---------------------------------------------------------------- ledger (Part D)
def ledger(fc):
    """Exact recomputation of the work ledger of proof.md section 5 (instruction counts / C, C = 20,000).
    fc carries the counted pairing and G2-output units of part F."""
    from decimal import Decimal, getcontext
    getcontext().prec = 60
    ln2 = Decimal(2).ln()
    lg = lambda v: float(Decimal(v).ln() / ln2)
    C = CHARGE_C
    # The squaring price is GS_CHARGE above: SQ_ITEMS sum to 9,708, Callgrind reads 79,527,954 Ir /
    # 8,192 calls = 9,708.0022 per call, and the charge is the ceiling of the larger, 9,709.
    # measured-kernel.md carries the itemisation, the kernel sources and the Callgrind record.
    GS = GS_CHARGE
    common_items, mult_items = COMMON_ITEMS, MULT_ITEMS
    common, extra = COMMON, MULT_EXTRA
    # ---------------------------------------------------------------------------------------------
    # Prices. Both charge vectors are recomputed from the same itemisation. The parent's ledger is
    # reproduced at its own squaring probability s = 3/8; this package uses the itemised squaring
    # charge and its re-tuned s = 35/64.
    PARENT = (11400, 300, 2200)          # the parent's rounded charges (submission 747ad6c9)
    OURS = (GS, common, extra)           # this package: itemised squaring charge and counted step charges
    if not (GS <= PARENT[0] and common <= PARENT[1] and extra <= PARENT[2]):
        raise RuntimeError('parent charges below derived counts')
    N = Fraction(r - 1, 12)
    s, Rt = Fraction(35, 64), 1024
    Pm = Decimal("0.43")
    eps = Decimal(5) / Decimal(12 * Rt) + Decimal(250) / Decimal(Rt * Rt)
    E = Decimal(2 ** 20) * 20 * Decimal(2 ** 60)
    SETUP = Decimal(2) ** Decimal("64.5")      # non-walk work, proof.md section 5

    def walk_ledger(SQ, COM, MULT, sq):
        """cbar, W, S, N' and sum q^2 for one charge vector and squaring probability."""
        sumq2 = sq * sq + (1 - sq) ** 2 / Rt
        Np = N / (1 - sumq2)
        Npd = Decimal(Np.numerator) / Decimal(Np.denominator)
        S = (2 * Npd * (1 / (1 - Pm)).ln()).sqrt()
        cbar = (1 - sq) * (1 + Fraction(COM + MULT, C)) + sq * (Fraction(SQ + COM, C))
        cb = Decimal(cbar.numerator) / Decimal(cbar.denominator)
        return cbar, cb, cb * S * (1 + eps) + cb * E + SETUP, S, Npd, sumq2

    cbar_p, cb_p, W_p, _, _, _ = walk_ledger(*PARENT, Fraction(3, 8))
    cbar, cb, W, S, Npd, sumq2 = walk_ledger(*OURS, s)
    SQ_CHARGE, COMMON_CHARGE, MULT_EXTRA_CHARGE = OURS
    # the parent's own ledger must come out of this same code path unchanged
    if not (abs(lg(W_p) - 125.71472902377629) < 1e-12 and cbar_p == Fraction(9225, 10000)):
        raise RuntimeError('parent ledger does not reproduce')
    # total-work abort (proof.md 3.4): plain arithmetic on the slack between the cap and the nominal ledger
    Wcap = Decimal(2) ** Decimal("125.62")     # the submitted bound
    if W > Wcap:
        raise RuntimeError("nominal ledger exceeds the submitted bound")
    n = S * (1 + eps)
    slack = (Wcap - cb * E - SETUP) / n - cb                 # per-step slack of the cap over the nominal ledger
    tol_sq_fraction = slack / Decimal("0.6205")                # tolerated drop of the squaring fraction below 35/64

    def success_with_drift(drop_points):
        """Success probability if the squaring fraction drops by `drop_points` percentage points.

        A drop by D raises the average step price by 0.6205*D units, so n stations of the cap buy
        n/(1+delta) steps; the collision law 1 - exp(-steps^2/(2N')) then gives the success value.
        """
        delta = Decimal(str(drop_points)) / 100 * Decimal("0.6205") / cb
        steps = min(n, Wcap / (cb * (1 + delta)))        # the cap binds once the drift is large
        return float(1 - (-(steps * steps) / (2 * Npd)).exp())

    drift_table = {p: success_with_drift(p) for p in ("0.5", "1", "1.5", "2", "3")}
    # setup and output (proof.md section 5): the per-instance pairing, the verification pairing and the
    # G2 output at its hard maximum come from part F's counts; every field operation there is one unit
    setup_items = {"pairing g = e(P,Q0) (counted)": fc["pairing_units_counted_g"],
                   "verification pairing e(P,Q) (counted)": fc["pairing_units_counted_verification"],
                   "check z^r = 1 (square-and-multiply)": r.bit_length() + bin(r).count("1"),
                   "table: 1024 x (2 exponentiations and their product <= 1,021, 11 Frobenius maps)": Rt * (1021 + 11),
                   "2^20 starts g^alpha z^beta at <= 1,021 each": 2 ** 20 * 1021,
                   "output Q = [x]Q0 in G2, hard maximum (counted)": fc["g2_output_units_hard_max_255_bit_scalar"]}
    setup_total = sum(setup_items.values())
    # report cap (proof.md 3.4): D <= 1.1 S theta + 2^25 reports. The store is a flat buffer of
    # 80-byte records appended as reports arise, then sorted by the 128-bit fingerprint with eight
    # stable 16-bit LSD radix passes and scanned for equal fingerprints. Each item below is counted
    # in proof.md section 3.2; the itemised 3,831 is charged at 4,000.
    report_items = {"12 keys for kappa": 2119,
                    "two 255-bit modular multiplications by p^kappa (2 x 512)": 1024,
                    "append of the 80-byte record (128-bit address and counters)": 32,
                    "eight 16-bit LSD radix passes over the 128-bit fingerprint (8 x 80)": 640,
                    "final equal-fingerprint scan": 16}
    report_itemised, report_charge = sum(report_items.values()), 4000
    dp_cap = Decimal(11) / 10 * S * Decimal(2) ** Decimal("-60") + Decimal(2) ** 25
    dp_handling = dp_cap * report_charge / C
    if not (report_itemised <= report_charge and setup_total <= 2 ** 31
            and fc["pairing_units_counted_g"] <= 2 ** 14
            and fc["pairing_units_counted_verification"] <= 2 ** 14 and dp_handling + 2 ** 31 < SETUP):
        raise RuntimeError('setup, output or non-walk bucket exceeds its stated bound')
    fp_exp = lg(dp_cap * dp_cap / Decimal(2) ** 129)
    sd6000 = 100 * math.sqrt(float(s * (1 - s)) / 6000)
    return {"squaring_items": SQ_ITEMS, "squaring_primitive_sum": SQ_PRIMITIVES, "squaring_items_sum": SQ_SUM, "squaring_charge": GS_CHARGE,
            "squaring_callgrind_cross_check": {"Ir": GS_IR, "calls": GS_CALLS,
                                              "per_call": GS_IR / GS_CALLS,
                                              "Ir_minus_items_x_calls": GS_IR - SQ_SUM * GS_CALLS},
            "kernel_measured_Ir": GS_IR, "kernel_measured_calls": GS_CALLS,
            "kernel_charge_rounding": "ceiling of max(item sum, Callgrind per call)",
            "log2_false_equal_F_pairs_expectation_at_cap": fp_exp,
            "log2_markov_more_than_2^16_false_runs": fp_exp - 16,
            "squaring_fraction_sampling_sd_points_6000_steps": sd6000,
            "nominal_ledger_within_cap": W <= Wcap,
            "fp2_sqr_instr": FP2SQR, "cyclotomic_square_instr_kernel": GS, "cyclotomic_square_instr_narrow_listing": GS_LISTING, "cyclotomic_square_charged": SQ_CHARGE,
            "common_items": common_items, "multiply_only_items": mult_items, "common_counted": common, "common_charged": COMMON_CHARGE,
            "multiply_extra_counted": extra, "multiply_extra_charged": MULT_EXTRA_CHARGE,
            "parent_charges": PARENT, "parent_cbar_units_per_step": float(cb_p), "parent_log2_W": lg(W_p),
            "sum_q2": float(sumq2), "log2_N": lg(Decimal(N.numerator) / Decimal(N.denominator)),
            "log2_N_effective": lg(Npd), "log2_S_model_0.43": lg(S), "cbar_units_per_step": float(cb),
            "fruitless_cycle_eps": float(eps), "log2_detection_window_steps": lg(E), "log2_W": lg(W),
            "submitted_bound": "125.62", "excess_of_cap_over_ledger": float(Wcap / W - 1),
            "work_abort_slack_per_step": float(slack), "tolerated_squaring_fraction_drop": float(tol_sq_fraction),
            "success_vs_squaring_fraction_drop_points": drift_table,
            "success_at_model_0.43": float(1 - (-(S * S) / (2 * Npd)).exp()),
            "setup_and_output_items_units": setup_items, "log2_setup_and_output_total": lg(setup_total),
            "setup_and_output_le_2^31": setup_total <= 2 ** 31,
            "log2_dp_report_cap": lg(dp_cap),
            "report_store": {"form": "flat buffer, eight stable 16-bit LSD radix passes on the 128-bit "
                                     "fingerprint, linear equal-fingerprint scan; no load factor",
                             "items": report_items, "itemised": report_itemised, "charged": report_charge},
            "log2_dp_handling": lg(dp_handling),
            "log2_report_memory_bytes_at_cap": lg(dp_cap * 80 + dp_cap * 32 + 2 ** 20),
            "dp_handling_plus_2^31_lt_2^64.5": dp_handling + 2 ** 31 < SETUP}


out = {}
out["partA_production_BLS12_381"] = part_A(rng, PROD_CHECKS)
if cfg.get("production_steps", 3000) > 0 and PROD_CHECKS > 0:
    out["partC_production_walk"] = part_C(rng, PROD_STEPS, stored=True)
    out["partC_production_walk_integer_rule"] = part_C(rng, PROD_STEPS, stored=False)
out["partB_toy"] = part_B(rng, Q_TOY, TRIALS, R, SQ64, DP_BITS, WALKS, BUDGET)
if PAIR_SAMPLES > 0:
    out["partE_pair_collision_law"] = part_E(rng, Q_TOY, PAIR_SAMPLES, R, SQ64)
out["partG_report_store"] = part_G(rng, 4096, 8)
pf = part_F(rng)
out["partF_reduction_and_output_BLS12_381"] = {k: v for k, v in pf.items() if not k.startswith("_")}
out["partD_ledger"] = ledger(pf)
out["seconds"] = round(time.time() - T0, 2)
out["scope"] = ("Parts A and C check exact identities and the production step rule on real BLS12-381 "
                "GT elements; Part F checks the reduction and the G2 output on BLS12-381; Part G checks "
                "the report store's fixed-pass sort and adjacent scan on data only; Part B is a toy "
                "field; none validates production costs.")
print(json.dumps(out, indent=1))
