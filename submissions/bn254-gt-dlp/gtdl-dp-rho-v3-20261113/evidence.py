"""bn254-gt-dlp: gated segmented distinguished-point Pollard rho.

Executable evidence.  A JSON request is read from stdin:

  {}                  default: profile checks, the production walk on the real
                      BN254 GT (invariant, partition histogram, DP rate), a real
                      recovery test with a known logarithm, the toy campaign and
                      the complete ledger.
  {"part":"toys"}     the toy campaign alone.
  {"part":"ledger"}   the ledger alone.
  {"part":"rebuild"}  recompute the rounded declarations.
  {"part":"kernel"}   replay the recorded kernel counts and the C walk statistics.

The toy campaign runs the *same* walk code as production (same increment
derivation, same partition rule, same distinguished-point rule, same segment
bookkeeping, same recovery), only the group and the encoding size differ.
Standard library only, deterministic.
"""
import hashlib
import json
import math
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path

HERE = Path("/snapshot/files") if Path("/snapshot/files").is_dir() else Path(__file__).resolve().parent


def load_profile():
    d = json.loads((HERE / "target-profile.json").read_text())
    return int(d["field_modulus"], 16), int(d["subgroup_order"], 16), \
        [int(x, 16) for x in d["generator_g"]["coefficients"]], d


P, R, GCOL, PROFILE = load_profile()
SQRT_R = math.isqrt(R)

# ------------------------------------------------------------- arithmetic ---
# Fp2 = Fp[u]/(u^2+1); Fp6 = Fp2[v]/(v^3-xi), xi = 9+u; Fp12 = Fp6[w]/(w^2-v).
# Profile: Fp[w]/(w^12-18w^6+82) with u = w^6-9, so a_m = (c_m+9c_{m+6}) + c_{m+6}u.

def f2mul(a, b):
    t0 = a[0] * b[0] % P
    t1 = a[1] * b[1] % P
    return ((t0 - t1) % P, ((a[0] + a[1]) * (b[0] + b[1]) - t0 - t1) % P)


def f2add(a, b): return ((a[0] + b[0]) % P, (a[1] + b[1]) % P)
def f2sub(a, b): return ((a[0] - b[0]) % P, (a[1] - b[1]) % P)
def mul_xi(a):   return ((9 * a[0] - a[1]) % P, (9 * a[1] + a[0]) % P)


def f6mul(a, b):
    v0 = f2mul(a[0], b[0]); v1 = f2mul(a[1], b[1]); v2 = f2mul(a[2], b[2])
    t0 = f2sub(f2sub(f2mul(f2add(a[1], a[2]), f2add(b[1], b[2])), v1), v2)
    t1 = f2sub(f2sub(f2mul(f2add(a[0], a[1]), f2add(b[0], b[1])), v0), v1)
    t2 = f2sub(f2sub(f2mul(f2add(a[0], a[2]), f2add(b[0], b[2])), v0), v2)
    return (f2add(v0, mul_xi(t0)), f2add(t1, mul_xi(v2)), f2add(t2, v1))


def f6add(a, b): return tuple(f2add(x, y) for x, y in zip(a, b))
def f6sub(a, b): return tuple(f2sub(x, y) for x, y in zip(a, b))
def mul_v(a):    return (mul_xi(a[2]), a[0], a[1])


def f12mul(A, B):
    a0, a1 = A; b0, b1 = B
    ac = f6mul(a0, b0); bd = f6mul(a1, b1)
    t = f6sub(f6sub(f6mul(f6add(a0, a1), f6add(b0, b1)), ac), bd)
    return (f6add(ac, mul_v(bd)), t)


ONE12 = (((1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)))


def f12pow(a, e):
    acc = ONE12
    while e:
        if e & 1:
            acc = f12mul(acc, a)
        a = f12mul(a, a)
        e >>= 1
    return acc


class Pow:
    """base^e by summing precomputed base^(2^i); ~127 multiplications per call."""

    def __init__(self, base, bits=254):
        self.p = [base]
        for _ in range(1, bits):
            self.p.append(f12mul(self.p[-1], self.p[-1]))

    def __call__(self, e):
        acc = ONE12
        i = 0
        while e:
            if e & 1:
                acc = f12mul(acc, self.p[i])
            e >>= 1
            i += 1
        return acc


def from_wbasis(c):
    c0 = [(0, 0)] * 3; c1 = [(0, 0)] * 3
    for m in range(6):
        j, k = m & 1, m >> 1
        a = ((c[m] + 9 * c[m + 6]) % P, c[m + 6] % P)
        if j == 0:
            c0[k] = a
        else:
            c1[k] = a
    return (tuple(c0), tuple(c1))


def to_wbasis(A):
    c0, c1 = A; out = [0] * 12
    for j, part in ((0, c0), (1, c1)):
        for k, a in enumerate(part):
            m = j + 2 * k
            out[m] = (out[m] + a[0] - 9 * a[1]) % P
            out[m + 6] = (out[m + 6] + a[1]) % P
    return out


def wbasis_mul(a, b):
    c = [0] * 23
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    c[i + j] = (c[i + j] + x * y) % P
    for k in range(22, 11, -1):
        v = c[k] % P
        if v:
            c[k] = 0
            c[k - 6] = (c[k - 6] + 18 * v) % P
            c[k - 12] = (c[k - 12] - 82 * v) % P
    return [x % P for x in c[:12]]


G12 = from_wbasis(GCOL)

# ------------------------------------------------------------- primality ----

def _mr(n, bases):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in bases:
        a %= n
        if a in (0, 1):
            continue
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


def prime(n):
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71]
    x = 0x9E3779B97F4A7C15
    for _ in range(8):
        x = (x * 0x2545F4914F6CDD1D + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        bases.append(x % n)
    return _mr(n, bases)


def order_of_p_mod_r():
    e = P % R
    k = 1
    while e != 1:
        e = e * (P % R) % R
        k += 1
        if k > 64:
            return None
    return k


# ===================================================================== walk ===
# One walk definition, used for the real GT and for every toy group.
#
# increments   c_j = int(SHA256("bn254gt-c"||j)) mod n,  d_j likewise
# partition    j = five bits from the middle of the state's canonical encoding
# step         Y <- Y * M_j,  M_j = g^{c_j} z^{d_j};  a,b <- a+c_j, b+d_j mod n
# segment      fresh uniform start, run to the first distinguished point
# DP rule      a fixed 128-bit window of the encoding compared with D
# table        (fingerprint, a, b), 80 bytes per entry
# recovery     x = (a - a')/(b' - b) mod n, then verify g^x == z

def increments(n, count=32):
    return [(int.from_bytes(hashlib.sha256(b"bn254gt-c" + j.to_bytes(4, "big")).digest(), "big") % n,
             int.from_bytes(hashlib.sha256(b"bn254gt-d" + j.to_bytes(4, "big")).digest(), "big") % n)
            for j in range(count)]


class GtWalk:
    """The production walk: states are Fp12 tower elements, x is unknown."""

    def __init__(self, z, d_threshold, count=32, gp=None, zp=None):
        self.z = z
        self.n = R
        self.gp = gp if gp is not None else Pow(G12)
        self.zp = zp if zp is not None else Pow(z)
        self.inc = increments(R, count)
        self.M = [f12mul(self.gp(c), self.zp(d)) for c, d in self.inc]
        self.C = [c for c, _ in self.inc]
        self.D = [d for _, d in self.inc]
        self.d_threshold = d_threshold

    def start(self, rng):
        a = rng.randrange(R); b = rng.randrange(R)
        return f12mul(self.gp(a), self.zp(b)), a, b

    def partition(self, y):
        # five bits from the middle of the canonical first coefficient
        return (y[0][0][0] >> 120) & 31

    def encode_window(self, y):
        """128-bit window: the low 128 bits of the same canonical first coefficient."""
        return y[0][0][0] & ((1 << 128) - 1)

    def is_dp(self, y):
        return self.encode_window(y) < self.d_threshold

    def fingerprint(self, y):
        h = hashlib.sha256(json.dumps([c for f2 in list(y[0]) + list(y[1]) for c in f2]).encode()).digest()
        return h[:16]

    def step(self, y, a, b):
        j = self.partition(y)
        return f12mul(y, self.M[j]), (a + self.C[j]) % R, (b + self.D[j]) % R, j

    def verify(self, x):
        return x is not None and 0 <= x < R and self.gp(x) == self.z


class ToyWalk:
    """The same walk on a prime-order subgroup of F_q^*, q = 2n+1 prime.

    The subgroup of squares is half of [0,q), so middle-bit partitions and
    trailing-bit distinguished points have their nominal probabilities to within
    the subgroup's density bias.  Increment derivation (SHA-256), partition
    placement, DP rule, segment bookkeeping and recovery are the production ones.
    """

    def __init__(self, q, n, g0, x, dp_threshold, count=32):
        self.q, self.n, self.g0 = q, n, g0
        self.bits = q.bit_length()
        self.inc = increments(n, count)
        self.M = [pow(g0, c, q) * pow(g0, d * x % n, q) % q for c, d in self.inc]
        self.C = [c for c, _ in self.inc]
        self.D = [d for _, d in self.inc]
        self.d_threshold = dp_threshold
        self.zp = pow(g0, x, q)

    def start(self, rng):
        a = rng.randrange(self.n); b = rng.randrange(self.n)
        return pow(self.g0, a, self.q) * pow(self.zp, b, self.q) % self.q, a, b

    def partition(self, y):
        return (y >> (self.bits // 2 - 2)) & 31

    def encode_window(self, y):
        return y & ((1 << (self.bits - 2)) - 1)

    def is_dp(self, y):
        return self.encode_window(y) < self.d_threshold

    def fingerprint(self, y):
        return y

    def step(self, y, a, b):
        j = self.partition(y)
        return y * self.M[j] % self.q, (a + self.C[j]) % self.n, (b + self.D[j]) % self.n, j

    def verify(self, x):
        return x is not None and pow(self.g0, x, self.q) == self.zp


def run_segments(walk, budget, w_seg, rng, seg_cap=64):
    """Segmented distinguished-point search; identical for production and toys."""
    table = {}
    mults = 0
    while mults < budget and len(table) < w_seg:
        y, a, b = walk.start(rng)
        limit = min(seg_cap * max(1, budget // w_seg), budget - mults)
        for _ in range(limit):
            y, a, b, _ = walk.step(y, a, b)
            mults += 1
            if walk.is_dp(y):
                break
        if mults >= budget:
            break
        if not walk.is_dp(y):
            continue
        key = walk.fingerprint(y)
        if key in table:
            a2, b2 = table[key]
            den = (b2 - b) % walk.n
            x = None if den == 0 else (a - a2) * pow(den, -1, walk.n) % walk.n
            return mults, walk.verify(x)
        table[key] = (a, b)
    return mults, False


# ------------------------------------------------------------------ toy part --

def toy_campaign(trials, bits_list, u_max, w_seg, seed=7):
    import random
    rng = random.Random(seed)
    out = []
    for bits in bits_list:
        n = (1 << bits) + 1
        q = 2 * n + 1
        while not (prime(n) and prime(q)):
            n += 2
            q = 2 * n + 1
        g0 = 2
        while pow(g0, n, q) != 1 or g0 == 1:
            g0 += 1
        budget = int(u_max * math.sqrt(n))
        L = max(1, budget // w_seg)
        thr = max(1, (1 << (q.bit_length() - 2)) // L)
        recs = []
        for _ in range(trials):
            x = rng.randrange(1, n)
            walk = ToyWalk(q, n, g0, x, thr)
            m, ok = run_segments(walk, budget, w_seg, rng)
            recs.append((m / math.sqrt(n), ok))
        succ = sorted(t for t, ok in recs if ok)
        out.append(dict(bits=bits, n_bits=n.bit_length(), q_bits=q.bit_length(), trials=trials,
                        successes=len(succ),
                        mean_u=round(sum(t for t, _ in recs) / len(recs), 4),
                        max_u=round(max(t for t, _ in recs), 4),
                        u_median_success=round(succ[len(succ) // 2], 4) if succ else None,
                        success_fraction_at_u={str(u): round(sum(1 for t in succ if t <= u) / trials, 3)
                                               for u in (1.2, 1.4, 1.6, 1.8, 2.0)}))
    return out


# ------------------------------------------------------- production-GT part --

def gt_recovery_test(seed=11):
    """Build genuine colliding states on the real GT and run the real recovery."""
    import random
    rng = random.Random(seed)
    walk = GtWalk(f12pow(G12, 12345), 1)          # threshold irrelevant here
    results = []
    for _ in range(4):
        x = rng.randrange(R)
        z = f12pow(G12, x)
        a, b = rng.randrange(R), rng.randrange(R)
        y = f12mul(f12pow(G12, a), f12pow(z, b))
        # a' + b'x = a + b x (mod r) with b' != b: pick b', solve for a'
        b2 = (b + 1 + rng.randrange(R - 1)) % R
        a2 = (a - x * (b2 - b)) % R
        y2 = f12mul(f12pow(G12, a2), f12pow(z, b2))
        den = (b2 - b) % R
        xs = (a - a2) * pow(den, -1, R) % R
        results.append(dict(equal=bool(y == y2), recovered=(xs == x), verified=bool(f12pow(G12, xs) == z)))
    return results


def gt_walk_stats(steps=1500, seed=5, d_threshold=1 << 118):
    """Partition histogram, DP rate and invariant checks on the real BN254 GT."""
    import random
    rng = random.Random(seed)
    x = 12345
    z = f12pow(G12, x)
    walk = GtWalk(z, d_threshold)
    y, a, b = walk.start(rng)
    hist = [0] * 32
    dps = 0
    bad = 0
    for i in range(steps):
        y, a, b, j = walk.step(y, a, b)
        hist[j] += 1
        if walk.is_dp(y):
            dps += 1
        if (i + 1) % 250 == 0 and y != f12pow(G12, (a + b * x) % R):
            bad += 1
    exp = steps / 32.0
    chi2 = sum((h - exp) ** 2 / exp for h in hist)
    return dict(steps=steps, partition_histogram=hist, chi2_31=round(chi2, 2),
                chi2_31_null_95=50.0, dp_hits=dps,
                dp_expected=steps / (1 << 10), dp_probability=2.0 ** -10,
                invariant_bad=bad)


# --------------------------------------------------------------- the ledger --

KERNEL = json.loads((HERE / "kernel-counts.json").read_text())
GTSTATS = json.loads((HERE / "gt-stats.json").read_text())
C_UNIT = 20000.0
# analytic AArch64 bound for one dense Fp12 multiplication in this tower
FP_MUL = 640          # 32 products + 32 carries + accumulator + select + memory
FP_ADD = 56           # add/sub/neg on four limbs including conditional fix-up
TOWER_OPS = 287       # counted in kernel.c: 90 Fp2 products, 156 Fp6 layers, 41 Fp12
TRAFFIC = 2 ** 13     # loads, stores, spills, call/return, loop control
FP12_BOUND = 54 * FP_MUL + TOWER_OPS * FP_ADD + TRAFFIC      # = 58,824
P_CONVERT = 736       # one Montgomery reduction to canonical form per step (bits are read canonically)
P_STEP = float(FP12_BOUND + P_CONVERT + 200)
U_STEPS = 1.6
P_SUCCESS_CHARGED = 0.66
ONE_TIME = 2 ** 40    # table build 2^30, 12,700 fresh-segment exponentiation pairs 2^39.4,
                      # DP bookkeeping 2^32, recovery/verification 2^22, randomness 2^20
DP_ENTRIES = 12700
DP_BYTES = 80
PEAK_BYTES = DP_ENTRIES * DP_BYTES + 32 * 448 + 4096
DECLARED = json.loads((HERE / "declared.json").read_text())


def ledger(u=U_STEPS, p_succ=P_SUCCESS_CHARGED):
    T = math.ceil(u * SQRT_R)
    instr = T * P_STEP + ONE_TIME
    w_run = instr / C_UNIT
    p_open = 0.39 / p_succ
    ew = p_open * w_run + (1.0 - p_open) * 1.0
    getcontext().prec = 60
    lg = Decimal(ew).ln() / Decimal(2).ln()
    up = lg.quantize(Decimal("0.000000001"), rounding=ROUND_CEILING)
    d_thr = (DP_ENTRIES << 128) // T
    return dict(T=str(T), log2_T=math.log2(T), per_step=P_STEP,
                units_per_step=P_STEP / C_UNIT, log2_w_run=math.log2(w_run),
                p_open=p_open, log2_E_W=float(lg), score_ceil_9dp=str(up),
                numerator=int(up * 10 ** 9), denominator=10 ** 9,
                dp_threshold=d_thr, dp_log2_probability=math.log2(d_thr / 2 ** 128),
                segment_length_log2=math.log2(T / DP_ENTRIES),
                peak_bytes=PEAK_BYTES)


# -------------------------------------------------------------------- parts --

def checks():
    out = {}
    out["p_prime"] = prime(P)
    out["r_prime"] = prime(R)
    out["p_congruent_1_mod_6"] = P % 6 == 1
    out["r_divides_p12_minus_1"] = (pow(P, 12, R) - 1) % R == 0
    out["ord_r_of_p"] = order_of_p_mod_r()
    out["g_is_not_identity"] = G12 != ONE12
    out["g_pow_r_is_identity"] = f12pow(G12, R) == ONE12
    out["bits"] = dict(p=P.bit_length(), r=R.bit_length())
    x = 0x1234567890ABCDEF
    same = []
    for _ in range(3):
        a, b = [], []
        for _ in range(12):
            x = (x * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
            a.append(x % P)
            x = (x * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
            b.append(x % P)
        same.append(to_wbasis(f12mul(from_wbasis(a), from_wbasis(b))) == wbasis_mul(a, b))
    out["tower_matches_wbasis_products"] = all(same)
    z = f12pow(G12, 12345)
    orbit = [z]
    for _ in range(12):
        orbit.append(f12pow(orbit[-1], P))
    out["frobenius_orbit_size"] = len(set(orbit[:12]))
    out["frobenius_returns_after_12"] = orbit[12] == orbit[0]
    out["inc_recovery"] = gt_recovery_test()
    out["inc_recovery_all_ok"] = all(all(v for v in r.values()) for r in out["inc_recovery"])
    out["gt_stats"] = gt_walk_stats()
    out["gt_stats_ok"] = (out["gt_stats"]["invariant_bad"] == 0 and
                          out["gt_stats"]["chi2_31"] < 60.0 and
                          abs(out["gt_stats"]["dp_hits"] - out["gt_stats"]["dp_expected"]) <= 6)
    return out


def main():
    req = json.load(sys.stdin) if not sys.stdin.isatty() else {}
    part = req.get("part", "default")
    res = {"part": part}
    if part in ("default", "all"):
        res["checks"] = checks()
        res["toy"] = toy_campaign(int(req.get("trials", 24)), [26, 30], 2.0, 256)
        L = ledger()
        res["ledger"] = L
        res["declared_match"] = (str(Decimal(L["score_ceil_9dp"])) == DECLARED["claim"]
                                 and L["numerator"] == DECLARED["score_numerator"])
        res["success_model"] = dict(charged=P_SUCCESS_CHARGED,
                                    random_map=round(1 - math.exp(-U_STEPS ** 2 / 2), 6),
                                    overall=round((0.39 / P_SUCCESS_CHARGED) * P_SUCCESS_CHARGED, 6))
    elif part == "toys":
        res["toy"] = toy_campaign(int(req.get("trials", 60)),
                                  [int(b) for b in req.get("bits", [26, 32, 38])],
                                  2.0, 256, seed=int(req.get("seed", 7)))
    elif part == "ledger":
        res["ledger"] = ledger()
    elif part == "rebuild":
        res["recomputed"] = ledger(); res["declared"] = DECLARED
    elif part == "kernel":
        res["kernel"] = KERNEL
        res["gt_stats"] = GTSTATS
        res["kernel_check"] = (KERNEL["ir_8192"] - KERNEL["ir_4096"]) / 4096.0
    elif part == "gt":
        res["checks"] = checks()
    else:
        res["error"] = "unknown part"
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
