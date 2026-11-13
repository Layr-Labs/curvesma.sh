"""bn254-gt-dlp: segmented distinguished-point Pollard rho on BN254 GT.

Executable evidence.  A JSON request is read from stdin:

  {}                  default: profile/field checks, real-GT walk invariant,
                      a short toy campaign, and the complete cost ledger.
  {"part":"toys"}     the full toy campaign.
  {"part":"ledger"}   the cost ledger only.
  {"part":"rebuild"}  recompute the rounded declarations from the constants.
  {"part":"kernel"}   replay the recorded kernel instruction counts.

Standard library only, deterministic.  The toy walks are evidence at their
actual size only: they are not a production DLP run and not a proof of the
collision-time heuristic at 2^254.
"""
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
# The protected profile's Fp12 is Fp[w]/(w^12-18w^6+82) with u = w^6-9, so the
# w-basis coefficients c_0..c_11 map to the tower as
#     a_m = (c_m + 9 c_{m+6}) + c_{m+6} u ,  m = 0..5,
# where m = j + 2k selects the Fp6 half j and the v^k slot.

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
    """Miller-Rabin, 20 fixed prime bases plus 8 derived bases (error < 4^-28)."""
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


# ------------------------------------------------------------- the walk -----
# Y_i = g^{a_i} z^{b_i}; the partition is the top five bits of Y's first
# coefficient; one step is one dense Fp12 multiplication by the table element
# M_j = g^{c_j} z^{d_j}; labels add (c_j, d_j) modulo r.

NTP = 32


def table_scalars(j):
    c = int.from_bytes(b"bn254gt-c" + j.to_bytes(4, "big"), "big") * 0x9E3779B97F4A7C15 % R
    d = int.from_bytes(b"bn254gt-d" + j.to_bytes(4, "big"), "big") * 0x2545F4914F6CDD1D % R
    return c, d


def partition(Y):
    return Y[0][0][0] >> 249


def build_table(z):
    return [(f12mul(f12pow(G12, c), f12pow(z, d)), c, d) for c, d in
            (table_scalars(j) for j in range(NTP))]


def gt_walk_check(z, x, steps=200, verify_every=40):
    """Walk the real BN254 GT and re-derive Y = g^{a+bx} independently."""
    tab = build_table(z)
    a, b = 12345 % R, 6789 % R
    Y = f12mul(f12pow(G12, a), f12pow(z, b))
    bad = 0
    for i in range(steps):
        m, c, d = tab[partition(Y)]
        Y = f12mul(Y, m)
        a = (a + c) % R
        b = (b + d) % R
        if (i + 1) % verify_every == 0:
            if Y != f12pow(G12, (a + b * x) % R):
                bad += 1
    return bad


# ------------------------------------------------------------- toy attack ---
# The same walk and the same segment/DP bookkeeping on prime-order subgroups of
# F_q^* small enough to enumerate.

def toy_group(bits):
    n = (1 << bits) + 1
    while not prime(n):
        n += 2
    k = 2
    while True:
        q = k * n + 1
        if prime(q):
            for h in range(2, 200):
                g0 = pow(h, k, q)
                if g0 != 1:
                    return q, n, g0
        k += 1


def toy_trial(q, n, g0, u_max, w_seg, rng):
    """Mirrors the production bookkeeping: fixed step budget T = u_max*sqrt(n),
    at most w_seg stored distinguished points, one DP per segment."""
    x = rng.randrange(1, n)
    zp = pow(g0, x, q)
    cs = [rng.randrange(n) for _ in range(32)]
    ds = [rng.randrange(n) for _ in range(32)]
    ms = [pow(g0, cs[j], q) * pow(zp, ds[j], q) % q for j in range(32)]
    budget = int(u_max * math.sqrt(n))
    seg = max(8, budget // w_seg)          # mean segment length = DP modulus
    table = {}
    mults = 0
    while mults < budget and len(table) < w_seg:
        a = rng.randrange(n); b = rng.randrange(n)
        y = pow(g0, a, q) * pow(zp, b, q) % q
        for _ in range(64 * seg):
            j = y % 32
            y = y * ms[j] % q
            a = (a + cs[j]) % n
            b = (b + ds[j]) % n
            mults += 1
            if y % seg == 0:
                break
        if y % seg:
            continue                       # cap hit without a DP: drop the segment
        if y in table:
            a2, b2 = table[y]
            den = (b2 - b) % n
            if den == 0:
                return mults, False
            xs = (a - a2) * pow(den, -1, n) % n
            return mults, pow(g0, xs, q) == zp
        table[y] = (a, b)
    return mults, False


def toy_campaign(trials, bits_list, u_max, w_seg, seed=7):
    import random
    rng = random.Random(seed)
    out = []
    for bits in bits_list:
        q, n, g0 = toy_group(bits)
        recs = []
        for _ in range(trials):
            m, ok = toy_trial(q, n, g0, u_max, w_seg, rng)
            recs.append((m / math.sqrt(n), ok))
        succ = sorted(t for t, ok in recs if ok)
        out.append(dict(
            bits=bits, n_bits=n.bit_length(), trials=trials,
            successes=len(succ),
            mean_u=round(sum(t for t, _ in recs) / len(recs), 4),
            max_u=round(max(t for t, _ in recs), 4),
            u_median_success=round(succ[len(succ) // 2], 4) if succ else None,
            u_90_success=round(succ[int(.9 * len(succ))], 4) if succ else None,
            success_fraction_at_u={
                str(u): round(sum(1 for t in succ if t <= u) / trials, 3)
                for u in (1.2, 1.4, 1.6, 1.8, 2.0)},
        ))
    return out


# --------------------------------------------------------------- the ledger --

KERNEL = json.loads((HERE / "kernel-counts.json").read_text())
C_UNIT = 20000.0
P_FP12 = 22996.0           # reference-measured dense Fp12 multiplication
P_STEP_OVERHEAD = 200.0    # partition, two 4-limb label updates, DP test, loop
P_STEP = P_FP12 + P_STEP_OVERHEAD
U_STEPS = 1.6
P_SUCCESS_CHARGED = 0.68   # theory 0.721963; toys measure 0.675/0.725
ONE_TIME = 2 ** 36         # table build, DP handling, recovery, verification
DP_ENTRIES = 12700
DP_BYTES = 80              # 16-byte fingerprint + two 32-byte labels
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
    return dict(T=str(T), log2_T=math.log2(T), log2_w_run=math.log2(w_run),
                p_open=p_open, log2_E_W=float(lg), score_ceil_9dp=str(up),
                numerator=int(up * 10 ** 9), denominator=10 ** 9,
                per_step=P_STEP, units_per_step=P_STEP / C_UNIT,
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
    out["gt_walk_bad_steps"] = gt_walk_check(z, 12345, steps=200, verify_every=50)
    out["gt_walk_invariant_ok"] = out["gt_walk_bad_steps"] == 0
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
        res["props"] = {
            "success_at_T_heuristic": round(1 - math.exp(-U_STEPS ** 2 / 2), 6),
            "success_charged": P_SUCCESS_CHARGED,
            "overall_success": round((0.39 / P_SUCCESS_CHARGED) * P_SUCCESS_CHARGED, 6),
        }
    elif part == "toys":
        res["toy"] = toy_campaign(int(req.get("trials", 60)),
                                  [int(b) for b in req.get("bits", [26, 32, 38])],
                                  2.0, 256, seed=int(req.get("seed", 7)))
    elif part == "ledger":
        res["ledger"] = ledger()
    elif part == "rebuild":
        L = ledger()
        res["recomputed"] = L
        res["declared"] = DECLARED
    elif part == "kernel":
        k = dict(KERNEL)
        k["per_mult_from_counters"] = (KERNEL["ir_8192"] - KERNEL["ir_4096"]) / 4096.0
        res["kernel"] = k
    else:
        res["error"] = "unknown part"
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
