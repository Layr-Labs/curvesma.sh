"""Frobenius-class Pollard rho: production identity checks, a toy-field run and
an exact ledger.

Part A checks, on BLS12-381 itself, the key identities that let the attack pick
a class-canonical Frobenius conjugate from four Fp multiplications.
Part B runs the complete walk (kappa-trick, 2-cycle escape, checkpoint cycle
escape, stuck-walk restart after 40/theta steps, distinguished points,
collision solving) in the order-r' subgroup of
F_{q^12}^* for a small q with the same tower shape as BLS12-381, and measures
the number of steps to the first useful class collision. It then replays each
trial's distinguished reports through the production store shape (a fixed LSD
radix sort of the report array followed by a linear scan of the runs of equal
class key) and checks the same logarithm comes out. Toy evidence only: it does
not validate production cost bounds, and it does not apply the production step
and report caps.
Part C recomputes the production ledger of proof.md from exact integers and
fractions and checks the final bound as an integer inequality.
"""
import json, math, random, sys, time
from fractions import Fraction

T0 = time.time()
cfg = json.load(sys.stdin)
q = cfg.get("q", 163)
trials = cfg.get("trials", 40)
R = cfg.get("R", 64)
dp_bits = cfg.get("dp_bits", 4)
walks = cfg.get("walks", 4)
seed = cfg.get("seed", 1)
prod_checks = cfg.get("production_checks", 2)
budget = cfg.get("seconds", 16)
rng = random.Random(seed)

# ---------- generic tower Fq2=Fq[u]/(u^2+1), Fq12=Fq2[w]/(w^6-xi), xi=1+u ----------
class Tower:
    def __init__(self, q):
        self.q = q
        self.xi = (1, 1)
        self.d = [self.f2pow(self.xi, m * (q - 1) // 6) for m in range(6)]
        z = self.f2pow(self.xi, (q * q - 1) // 6)
        assert z[1] == 0 and (z[0] * z[0] - z[0] + 1) % q == 0, "zeta must be a primitive 6th root in Fq"
        self.zeta = z[0]
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

# ---------- Part A: production BLS12-381 identities ----------
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P = Tower(p)
partA = {"zeta_in_Fp_primitive_6th_root": True,
         "p_order_mod_r_is_12": pow(p, 12, r) == 1 and all(pow(p, k, r) != 1 for k in range(1, 12)),
         "one_plus_p4_plus_p8_zero_mod_r": (1 + pow(p, 4, r) + pow(p, 8, r)) % r == 0}
ok_keys = ok_frob = ok_canon = ok_gt = ok_inv = True
hcof = (p ** 12 - 1) // r
for _ in range(prod_checks):
    x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
    y = P.pow(x, hcof)                                    # random GT element
    ok_gt &= P.pow(y, r) == P.one() and y != P.one()
    ok_frob &= P.frob(y) == P.pow(y, p)                   # Frobenius formula = p-th power
    kd, kf = P.keys_direct(y), P.keys_formula(y)
    ok_keys &= kd == kf and len(set(kd)) == 12
    conj = [y]
    for k in range(11): conj.append(P.frob(conj[-1]))
    canon = set()
    for c in conj:
        kk = P.keys_formula(c); kap = min(range(12), key=lambda i: kk[i])
        t = c
        for _ in range(kap): t = P.frob(t)
        canon.add(t)
    ok_canon &= len(canon) == 1
    ok_inv &= P.mul(conj[6], y) == P.one()
partA.update({"random_GT_elements": prod_checks, "GT_membership": ok_gt,
              "frobenius_formula_equals_pth_power": ok_frob,
              "four_mult_key_formula_equals_direct_keys_and_distinct": ok_keys,
              "kappa_canonical_rep_identical_over_all_12_conjugates": ok_canon,
              "p6_power_is_inverse": ok_inv})

# ---------- Part B: toy-field run of the complete algorithm ----------
T = Tower(q)
n = q ** 4 - q * q + 1
rp = n
for sp in range(2, 100000):
    while rp % sp == 0 and rp > sp: rp //= sp
N = (rp - 1) // 12
qk = [pow(q, k, rp) for k in range(12)]
hc = (q ** 12 - 1) // rp
assert pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))

Dk = [[T.f2pow(T.xi, m * (q ** k - 1) // 6) for m in range(6)] for k in range(12)]
def frobk(y, k):
    if k % 2: y = tuple((c[0], (-c[1]) % q) for c in y)
    return tuple(T.f2mul(c, Dk[k][m]) for m, c in enumerate(y))

def keyinfo(y):
    """kappa and class key. Production key: the minimal Fp value from keys_formula.
    q is tiny here, so equal values occur; the toy breaks such ties by comparing
    the full conjugates (in production a tie has probability about 2^-370 per step)."""
    kf = T.keys_formula(y)
    m = min(kf)
    ks = [k for k in range(12) if kf[k] == m]
    if len(ks) == 1:
        kap = ks[0]; t = frobk(y, kap)
    else:
        t, kap = min((frobk(y, k), k) for k in ks)
    return kap, (m, t)

def H(key):
    h = 0
    for c in key[1]:
        h = (h * 1000003 + c[0] * 7919 + c[1]) % (1 << 61)
    return h

# ---------- batched distinguished store: fixed LSD radix sort + linear run scan ----------
# The production store keeps every report in a flat array and merges after the walk:
# a fixed LSD radix sort over the whole class key, then one linear scan of the runs.
# The toy replays each trial through that same shape. It sorts by a 192-bit
# serialisation of the canonical representative (12 coefficients x 16 bits, q < 2^8)
# with 8-bit digits; production sorts the 381-bit scalar class key K with 16-bit digits.
STORE_KEY_BITS = 192
STORE_DIGIT_BITS = 8
STORE_PASSES = STORE_KEY_BITS // STORE_DIGIT_BITS
STORE_REPLAY_TRIALS = 20        # bound the replay so the runner's time limit keeps its margin

def serialise(cls):
    v = 0
    for c in cls:
        v = (v << 16) | (c[0] << 8) | c[1]
    return v

def store_solve(reports, g, z, rp):
    """Sort (class key, ca, cb) records by the class key, then scan runs linearly:
    keep the first record of a run and compare every later record with it and with
    its predecessor. Returns x on the first verified pair, else None."""
    keys = [r[0] for r in reports]
    idx = list(range(len(keys)))
    mask = (1 << STORE_DIGIT_BITS) - 1
    for p in range(STORE_PASSES):
        sh = p * STORE_DIGIT_BITS
        cnt = [0] * (1 << STORE_DIGIT_BITS)
        for i in idx:
            cnt[(keys[i] >> sh) & mask] += 1
        s = 0
        for d in range(1 << STORE_DIGIT_BITS):
            c = cnt[d]; cnt[d] = s; s += c
        out = [0] * len(idx)
        for i in idx:
            d = (keys[i] >> sh) & mask
            out[cnt[d]] = i; cnt[d] += 1
        idx = out
    run_first = None
    prev = None
    for i in idx:
        k = keys[i]
        if prev is None or k != keys[prev]:
            run_first = i; prev = i; continue
        for j in (run_first, prev):
            if j == i: continue
            ca1, cb1 = reports[j][1], reports[j][2]
            ca2, cb2 = reports[i][1], reports[i][2]
            if (ca1, cb1) == (ca2, cb2): continue
            xr = (ca1 - ca2) * pow(cb2 - cb1, -1, rp) % rp
            if T.pow(g, xr) == z:
                return xr
        prev = i
    return None

results = []; Tc_list = []; total_list = []; stats = {"two_cycles": 0, "long_cycles": 0, "useless": 0, "abandoned": 0}
store = {"trials": 0, "recovered": 0, "records": 0}
while len(results) < trials and time.time() - T0 < budget:
    while True:
        g = T.pow(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(6)), hc)
        if g != T.one(): break
    xs = rng.randrange(rp); z = T.pow(g, xs)
    ah = [rng.randrange(rp) for _ in range(R)]; bh = [rng.randrange(rp) for _ in range(R)]
    Mh = [T.mul(T.pow(g, ah[i]), T.pow(z, bh[i])) for i in range(R)]
    Tab = []; At = []; Bt = []
    for i in range(R):
        row = [Mh[i]]
        for k in range(11): row.append(T.frob(row[-1]))
        Tab.append(row); At.append([ah[i] * qk[k] % rp for k in range(12)]); Bt.append([bh[i] * qk[k] % rp for k in range(12)])
    dps = {}; visited = {}; steps = 0; Tc = None; found = None; reports = []
    def new_walk():
        a, b = rng.randrange(rp), rng.randrange(rp)
        return {"y": T.mul(T.pow(g, a), T.pow(z, b)), "a": a, "b": b, "hist": [], "ck": None, "len": 0, "sincedp": 0}
    W = [new_walk() for _ in range(walks)]
    steps += 0
    while found is None:
        for w in W:
            y = w["y"]; kap, key = keyinfo(y)
            ca, cb = w["a"] * qk[kap] % rp, w["b"] * qk[kap] % rp     # canonical exponents
            cls = key[1]
            if Tc is None:
                prev = visited.get(cls)
                if prev is None: visited[cls] = (ca, cb)
                elif prev != (ca, cb): Tc = steps
            hv = H(key)
            if hv % (1 << dp_bits) == 0:                               # distinguished class
                reports.append((serialise(cls), ca, cb))               # batched store record
                old = dps.get(cls)
                if old is not None and old[1] != cb:
                    xr = (old[0] - ca) * pow(cb - old[1], -1, rp) % rp
                    assert T.pow(g, xr) == z
                    found = xr; break
                if old is not None: stats["useless"] += 1
                dps[cls] = (ca, cb)
                w["sincedp"] = 0
            hist = w["hist"]
            if len(hist) >= 2 and hist[-2][0][1] == cls:                 # fruitless 2-cycle
                stats["two_cycles"] += 1
                (k1, y1, a1, b1), (k2, y2, a2, b2) = hist[-1], (key, y, w["a"], w["b"])
                (ky, yy, aa, bb) = min(((k1, y1, a1, b1), (k2, y2, a2, b2)), key=lambda t: t[0])
                w["y"] = T.mul(yy, yy); w["a"] = 2 * aa % rp; w["b"] = 2 * bb % rp
                w["hist"] = []; w["ck"] = None; steps += 1; continue
            if w["ck"] is not None and w["ck"][0] == cls and w["ck"][1] < w["len"]:   # longer cycle
                stats["long_cycles"] += 1
                # walk once more around the cycle to find its minimal class, then square it
                best = (key, y, w["a"], w["b"]); yy, aa, bb = y, w["a"], w["b"]
                for _ in range(w["len"] - w["ck"][1]):
                    kk, ke = keyinfo(yy); hh = H(ke) % R; idx = (12 - kk) % 12
                    yy = T.mul(yy, Tab[hh][idx]); aa = (aa + At[hh][idx]) % rp; bb = (bb + Bt[hh][idx]) % rp
                    steps += 1
                    k2, ke2 = keyinfo(yy)
                    if ke2 < best[0]: best = (ke2, yy, aa, bb)
                _, yy, aa, bb = best
                w["y"] = T.mul(yy, yy); w["a"] = 2 * aa % rp; w["b"] = 2 * bb % rp
                w["hist"] = []; w["ck"] = None; steps += 1; continue
            if w["len"] % 64 == 0: w["ck"] = (cls, w["len"])
            hist.append((key, y, w["a"], w["b"]))
            if len(hist) > 2: hist.pop(0)
            h = hv % R; idx = (12 - kap) % 12
            w["y"] = T.mul(y, Tab[h][idx]); w["a"] = (w["a"] + At[h][idx]) % rp; w["b"] = (w["b"] + Bt[h][idx]) % rp
            w["len"] += 1; w["sincedp"] += 1; steps += 1
            if w["sincedp"] > 40 << dp_bits:
                stats["abandoned"] += 1; w.update(new_walk())
    assert found == xs % rp or T.pow(g, found) == z
    results.append(True); Tc_list.append(Tc); total_list.append(steps)
    if store["trials"] < STORE_REPLAY_TRIALS:                          # batched store replay
        xb = store_solve(reports, g, z, rp)
        store["trials"] += 1; store["records"] += len(reports)
        if xb is not None and T.pow(g, xb) == z:
            store["recovered"] += 1

sq = math.sqrt(N)
b40 = math.sqrt(2 * N * math.log(1 / 0.6) * R / (R - 1))
b39 = math.sqrt(2 * N * math.log(1 / 0.61) * R / (R - 1))
out = {"partA_production_BLS12_381": partA,
       "partB_toy": {"q": q, "r_prime": rp, "classes_N": N, "R": R, "dp_bits": dp_bits, "walks": walks,
                     "trials_completed": len(results), "all_logs_recovered_and_verified": all(results),
                     "first_useful_class_collision_mean_over_sqrtN": round(sum(Tc_list) / len(Tc_list) / sq, 4) if Tc_list else None,
                     "random_mapping_prediction_mean_over_sqrtN": round(math.sqrt(math.pi / 2), 4),
                     "fraction_Tc_le_budget39": round(sum(t <= b39 for t in Tc_list) / max(1, len(Tc_list)), 3),
                     "budget39_over_sqrtN": round(b39 / sq, 4),
                     "fraction_Tc_le_budget40": round(sum(t <= b40 for t in Tc_list) / max(1, len(Tc_list)), 3),
                     "budget40_over_sqrtN": round(b40 / sq, 4),
                     "total_steps_incl_DP_delay_mean_over_sqrtN": round(sum(total_list) / max(1, len(total_list)) / sq, 4),
                     "events": stats,
                     "two_cycles_per_step": stats["two_cycles"] / max(1, sum(total_list)),
                     "predicted_two_cycles_per_step": 1 / (12 * R),
                     "batched_sorted_store": {"passes": STORE_PASSES, "digit_bits": STORE_DIGIT_BITS,
                                              "key_bits": STORE_KEY_BITS,
                                              "records": store["records"],
                                              "trials_replayed": store["trials"],
                                              "trials_recovered": store["recovered"]}},
       "Tc_over_sqrtN": [round(t / sq, 4) for t in Tc_list] if cfg.get("dump") else None,
       "scope": "Part A checks exact identities on real BLS12-381 GT elements. Part B is a toy field; it does not validate production cost bounds, and its store replay sorts a 192-bit toy class serialisation, not the production 381-bit class key. Part C recomputes the production ledger exactly."}

# ---------- Part C: exact ledger (integers and fractions only) ----------
def _log2(x):
    """float log2 of a positive int/Fraction, for printing only."""
    if isinstance(x, int): x = Fraction(x)
    n, d = x.numerator, x.denominator
    bn, bd = n.bit_length(), d.bit_length()
    return (bn - bd) + math.log2((n / (1 << bn)) / (d / (1 << bd)))

def _ln_upper(num, den, terms=40):
    """Exact rational upper bound on ln(num/den) for num > den > 0."""
    t = (Fraction(num, den) - 1) / (Fraction(num, den) + 1)   # ln = 2 atanh t
    s = sum((Fraction(t ** (2 * k + 1), 2 * k + 1) for k in range(terms)), Fraction(0))
    return 2 * (s + t ** (2 * terms + 1) / ((2 * terms + 1) * (1 - t * t)))

def _sqrt_upper(v, bits=400):
    """Rational upper bound on sqrt(v) for a positive Fraction v."""
    n, d = v.numerator, v.denominator
    return Fraction(math.isqrt(n << (2 * bits)) + 1, math.isqrt(d << (2 * bits)))

Rprod = 1024
Nprod = (r - 1) // 12
Np = Fraction(Nprod * Rprod, Rprod - 1)
M = 1 << 20
E = M * 40 * (1 << 60)                       # detection window M*40/theta
LN = _ln_upper(100, 61)                      # ln(1/0.61) upper bound
S0 = _sqrt_upper(2 * Np * LN)                # P = 0.39 step point
S = Fraction(103, 100) * S0                  # 1.03 slack
CAP = Fraction(1003, 1000)
Lmax = CAP * S + E
c = Fraction(1) + Fraction(2600, 20000)      # 1.13
# additive terms
setup = Fraction(1 << 31)                    # pairings, table, starts, output; itemised sum < 2^31
reports = (1 << 67) + 1                      # hard report cap
# batched sorted store: 112-byte records, fixed 24-pass 16-bit radix sort, linear run scan
dp_handle = Fraction(reports * 2076 + (1 << 25), 20000)   # 24 x 16 x 2^16 counter work, rounded up
restart_count = Lmax // (40 * (1 << 60)) + 1
restarts = Fraction(restart_count) * 1021
solves = Fraction(reports * 1040)            # at most two candidate solves per record
Tadd = setup + dp_handle + restarts + solves
W = c * (CAP * S + E) + Tadd
lw = _log2(W)
claim = Fraction(math.ceil(lw * 100), 100)
nW, dW = W.numerator, W.denominator
exact_ok = nW ** 100 <= dW ** 100 * (1 << 12586)
exact_tight = nW ** 100 > dW ** 100 * (1 << 12585)
mem_records = (1 << 67) * 112                # report records
mem_idx = 2 * (1 << 67) * 12                 # two radix-scatter index buffers
mem_walk = 1 << 31
mem_mult = 1024 * 12 * 576 + 1024 * 12 * 64  # twisted table + exponent tables
mem_misc = 1 << 20
partC = {"log2_p": round(_log2(p), 4), "log2_N": round(_log2(Nprod), 4),
         "log2_N_prime": round(_log2(Np), 4), "ln_1_over_0.61_upper_float": round(float(LN), 18),
         "log2_S0_P039": round(_log2(S0), 4), "log2_S": round(_log2(S), 4),
         "model_success_at_S": round(1 - math.exp(-float(Fraction(103, 100) ** 2 * LN)), 4),
         "log2_L_max": round(_log2(Lmax), 4),
         "log2_E": round(_log2(E), 4),
         "restart_count_le": restart_count, "log2_restart_count": round(_log2(restart_count), 4),
         "log2_setup_term": round(_log2(setup), 4),
         "instructions_per_report": 2076,
         "log2_dp_handle_term": round(_log2(dp_handle), 4),
         "log2_restarts_term": round(_log2(restarts), 4),
         "solves_per_record": 2,
         "log2_solves_term": round(_log2(solves), 4),
         "log2_additive_total": round(_log2(Tadd), 4),
         "c_per_step": float(c), "log2_W": round(lw, 6),
         "claim_log2": float(claim),
         "slack_bits": round(float(claim) - lw, 6),
         "exact_W_to_100_le_2_pow_12586": exact_ok,
         "exact_W_to_100_gt_2_pow_12585": exact_tight,
         "memory_bytes": mem_records + mem_idx + mem_walk + mem_mult + mem_misc,
         "memory_bytes_log2": round(_log2(mem_records + mem_idx + mem_walk + mem_mult + mem_misc), 4),
         "tie_probability_log2_over_run": round(_log2(Fraction(66, 1) * (Lmax + 1) / p), 2),
         "note": "All fields are exact-rational computations; only the printed log2 values are rounded, and the two exact_ fields are the integer certificate."}
out["partC_exact_ledger"] = partC
out["seconds"] = round(time.time() - T0, 2)
print(json.dumps(out, indent=1))
