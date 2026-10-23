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
radix sort of the report array on the FULL canonical class key, followed by a
linear scan of the runs of equal key) and checks the same logarithm comes out.
Toy evidence only: it does
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
    """kappa and the class pair (minimal key, full canonical representative).
    The step function uses the minimal key, exactly as production does; the store
    uses the full canonical representative, which is the exact class identity.
    q is tiny here, so equal minimal keys occur; the toy breaks such ties by
    comparing the full conjugates (in production a tie has probability about
    2^-370 per step)."""
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
# The class key is the FULL canonical representative of the class (every Fp
# coefficient, no truncation and no hash), so two records share a key exactly when
# they are in the same class. The toy replays each trial through that same shape:
# it sorts the same full canonical key, serialised exactly (12 coefficients, q < 2^8,
# 16 bits each), with 8-bit digits. Production sorts the 4,608-bit key with 288
# 16-bit digits.
STORE_KEY_BITS = 192            # 12 Fp coefficients x 16 bits, exact at q < 2^8
STORE_DIGIT_BITS = 8
STORE_PASSES = STORE_KEY_BITS // STORE_DIGIT_BITS
STORE_REPLAY_TRIALS = 20        # bound the replay so the runner's time limit keeps its margin

def serialise(cls):
    """Exact full-canonical-key serialisation: every coefficient, no hash. Two
    classes collide here only if they are the same class."""
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
                                              "key": "full canonical class key, exact, no hash",
                                              "records": store["records"],
                                              "trials_replayed": store["trials"],
                                              "trials_recovered": store["recovered"]}},
       "Tc_over_sqrtN": [round(t / sq, 4) for t in Tc_list] if cfg.get("dump") else None,
       "scope": "Part A checks exact identities on real BLS12-381 GT elements. Part B is a toy field; it does not validate production cost bounds, and its store replay sorts the exact full canonical class key of the toy field (12 coefficients, no hash), which is the same class-identity map production uses. Part C recomputes the production ledger exactly."}

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

def _exp_neg_upper(a, terms=200):
    """Exact rational upper bound on e^-a (a > 0 int): 1 / (partial Taylor sum of e^a)."""
    s, t = Fraction(0), Fraction(1)
    for k in range(terms):
        s += t; t = t * a / (k + 1)
    return 1 / s

def _dyadic(x, up, bits=160):
    """(c, e) with c * 2^-e >= x (up) or <= x (down), c an integer of about `bits` bits."""
    x = Fraction(x)
    e = bits - (x.numerator.bit_length() - x.denominator.bit_length())
    n, d = (x.numerator << e, x.denominator) if e >= 0 else (x.numerator, x.denominator << -e)
    return (-(-n // d) if up else n // d), e

def _vs_pow2(x, printed, up):
    """Exact test: x <= 2^printed (up) or x >= 2^printed (down); printed a decimal string."""
    k = Fraction(printed); num, den = k.numerator, k.denominator
    c, e = _dyadic(x, up)
    t = num + e * den                         # (c 2^-e)^den  vs  2^num  <=>  c^den vs 2^t
    lhs, rhs = c ** den, 1
    if t >= 0: rhs <<= t
    else: lhs <<= -t
    return lhs <= rhs if up else lhs >= rhs

Rprod = 4096
Nprod = (r - 1) // 12
Np = Fraction(Nprod * Rprod, Rprod - 1)
M = 1 << 20
THETA_INV = 1 << 80
E = M * 40 * THETA_INV                       # detection window M*40/theta = 40 * 2^100
LN = _ln_upper(100, 61)                      # ln(1/0.61) upper bound
S0 = _sqrt_upper(2 * Np * LN)                # P = 0.39 step point
S = Fraction(103, 100) * S0                  # 1.03 slack
CAP = Fraction(101, 100)                     # 0.008 S cycle waste + 0.002 S abandonment waste
Lmax = CAP * S + E
c = Fraction(1) + Fraction(2600, 20000)      # 1.13
# setup items (units), each an itemised worst case; their sum is below the declared 2^31
PAIRING = 8921                               # Miller loop <= 292.8 + final exponentiation <= 8,628
TABLE_ENTRY = 1020 + 1 + 11 + 1              # 2 exps, product, 11 Frobenius, 22 exponent-table mults mod r (22 x 512 <= 20,000)
setup_items = {"pairing_g": PAIRING, "table_4096_entries": Rprod * TABLE_ENTRY,
               "walk_starts_2^20_x_1021": M * 1021, "output_735_plus_inversion_15": 750,
               "verification_pairing": PAIRING}
setup_sum = sum(setup_items.values())
setup = Fraction(1 << 31)                    # declared setup cap, what enters W
reports = (1 << 54) + 1                      # hard report cap, inside one 64-bit address space
# report handling per record: append 80 words 112, canonical exponents 2 x 512,
# canonical representative 12 x 408 budgeted 5,120, identity index 12,
# 288 radix passes x 64, run scan 512. Every address, index and counter below is one 64-bit word.
MODMUL_R = 512                               # one multiplication mod r, itemised 408 (see proof.md)
CANON_REP = 5120                             # full canonical representative, 12 Fp values
RECORD_APPEND = 112                          # 80 word stores + 64-bit address, bound check, counter
SORT_PASSES = 288                            # 4,608-bit key in 16-bit digits
IDENTITY_INDEX = 12
PER_PASS = 64                                # single-word counting and scatter sweeps
RUN_SCAN = 512
PER_RECORD = RECORD_APPEND + 2 * MODMUL_R + CANON_REP + IDENTITY_INDEX + SORT_PASSES * PER_PASS + RUN_SCAN
RECORD_ITEMISED = 98 + 2 * 408 + 12 * 408 + 8 + SORT_PASSES * 52 + 480
COUNTER_WORK = SORT_PASSES * 16 * (1 << 16)  # zero + prefix-sum the 2^16 counters, 288 passes
dp_instr = reports * PER_RECORD + COUNTER_WORK
dp_handle = Fraction(dp_instr, 20000)
restart_count = Lmax // (40 * THETA_INV) + 1
restarts = Fraction(restart_count) * 1021
solves = Fraction(reports * 1040)            # two 520-unit solves per record, in units
# randomness (H3 of proof.md section 6): ONE fresh uniform 256-bit seed K from the platform
# source at start, charged as one 64-word draw; every scalar is a 384-bit Keccak-f[1600]
# counter-mode output under K, reduced mod r at 8,192 budgeted instructions per scalar.
KECCAK_ROUND = {"theta": 55, "rho_pi": 49, "chi": 100, "iota": 2, "loads_stores_spills": 90, "round_loop": 4}
KECCAK_PERM = 24 * sum(KECCAK_ROUND.values())
SAMPLE_ITEMS = {"keccak_f1600_24_rounds": KECCAK_PERM, "absorb_seed_pad_counter_squeeze": 64, "draw_loads": 48,
                "barrett_y_mu": 126, "qhat_r": 80, "y_minus_qhat_r": 21, "two_cond_sub": 18, "control": 24}
SAMPLE_ITEMISED = sum(SAMPLE_ITEMS.values())
SAMPLE_BUDGET = 8192
DRAW_BITS = 384
sample_sites = 2 * (Rprod + M + restart_count + 1)   # table entries, starts, restarts, +1 spare
sampler = Fraction(sample_sites * SAMPLE_BUDGET, 20000)
# the one fresh seed: a 256-bit key, drawn as one 64-word block and charged whole
SEED_BITS = 256
SEED_WORDS = 64
SEED_ITEMS = {"getrandom_512_bytes_syscall": 4096, "kernel_csprng_512_bytes_at_64_per_byte": 32768,
              "word_loads_stores_64": 768, "absorb_seed_counter_padding": 1920,
              "buffer_counters_control": 8192}
SEED_ITEMISED = sum(SEED_ITEMS.values())
SEED_BUDGET = 1 << 16
seed_draw = Fraction(SEED_BUDGET, 20000)
Tadd = setup + dp_handle + restarts + solves + sampler + seed_draw
W = c * (CAP * S + E) + Tadd
lw = _log2(W)
claim = Fraction(math.ceil(lw * 100), 100)
nW, dW = W.numerator, W.denominator
exact_ok = nW ** 100 <= dW ** 100 * (1 << 12587)
exact_tight = nW ** 100 > dW ** 100 * (1 << 12586)
# success-side quantities: Markov tails on the H1 means (no per-step independence is used)
mean_reports = Lmax / THETA_INV
report_tail = mean_reports / (1 << 54)                       # P(reports >= 2^54)
stretches = reports + M + restart_count
e40 = _exp_neg_upper(40)
abandon_mean = stretches * e40
abandon_thresh = Fraction(2, 1000) * S / (40 * THETA_INV)   # abandonments that strand 0.002 S steps
abandon_markov = abandon_mean / abandon_thresh
ties = 66 * (Lmax + 1) / p
waste_rate = Fraction(3, 12 * Rprod) + Fraction(200, Rprod * Rprod)
mean_waste = waste_rate * Lmax
cycle_tail = mean_waste / (Fraction(8, 1000) * S)            # P(waste > 0.008 S)
long_cycle = Fraction(1, 1 << 654)                           # cycle of length > 64 inside the run
self_collision = Fraction(1, M)
model = 1 - math.exp(-float(Fraction(103, 100) ** 2 * LN))
success = (Fraction(model) - cycle_tail - report_tail - abandon_markov - e40
           - ties - long_cycle - self_collision)
mem_records = (1 << 54) * 640                # exact class key 576 B + two 32-byte exponents
mem_idx = 2 * (1 << 54) * 8                  # two radix-scatter 8-byte index buffers
mem_walk = 1 << 31
mem_mult = Rprod * 12 * 576 + Rprod * 12 * 64 + 144 * 48  # twisted table, exponent tables, Frobenius constants
mem_digit = (1 << 16) * 8                    # digit counters, 2^16 x 8
mem_misc = 1 << 20                           # scalars, caps, code, stack, tape buffer
mem_peak = mem_records + mem_idx + mem_walk + mem_mult + mem_digit + mem_misc

# ---------- safe-rounding self-check: every quoted log2 figure against its exact value ----------
# kind "up": the text uses it as an upper bound, so the exact value must be <= 2^printed;
# kind "down": a lower bound, exact >= 2^printed.
FIGURES = [
    ("S0", S0, "125.6280", "up"), ("S", S, "125.6707", "up"), ("L_max", Lmax, "125.6851", "up"),
    ("E 40 x 2^100", E, "105.33", "up"), ("restart cap n_r", restart_count, "40.37", "up"),
    ("samples", sample_sites, "41.39", "up"), ("sampler term", sampler, "40.11", "up"),
    ("table setup 4,231,168", setup_items["table_4096_entries"], "22.02", "up"),
    ("setup itemised sum", setup_sum, "31", "up"),
    ("non-record arrays 2,180,520,704", mem_walk + mem_mult + mem_digit + mem_misc, "31.03", "up"),
    ("handling instructions", dp_instr, "68.63", "up"), ("handling term", dp_handle, "54.35", "up"),
    ("counter work", COUNTER_WORK, "28.17", "up"),
    ("restarts term", restarts, "50.37", "up"), ("solves term", solves, "64.03", "up"),
    ("additive total T", Tadd, "64.04", "up"), ("T exact 4 dp", Tadd, "64.0243", "up"),
    ("W", W, "125.8614", "up"), ("W claim", W, "125.87", "up"),
    ("W not at 125.86", W, "125.86", "down"),
    ("fresh seed draw", seed_draw, "1.72", "up"),
    ("fresh seed draw instructions 2^16", SEED_BUDGET, "16", "up"),
    ("mean reports", mean_reports, "45.69", "up"), ("report-cap tail", report_tail, "-8.31", "up"),
    ("cycle-waste tail", cycle_tail, "-6.76", "up"),
    ("stretches", stretches, "54.02", "up"),
    ("e^-40", e40, "-57.7", "up"), ("expected abandonments", abandon_mean, "-3.67", "up"),
    ("Markov abandonment tail", abandon_markov, "-35.05", "up"),
    ("restart samples 2 n_r", 2 * restart_count, "41.37", "up"),
    ("abandonment threshold", abandon_thresh, "31.38", "down"),
    ("ties over run", ties, "-248.9", "up"),
    ("2-cycle rate 1/(12R)", Fraction(1, 12 * Rprod), "-15.58", "up"),
    ("long-cycle case", long_cycle, "-654", "up"),
    ("memory peak", mem_peak, "63.36", "up"), ("memory < 2^64", mem_peak, "64", "up"),
    ("records 640 x 2^54", mem_records, "63.33", "up"),
]
VALUES = [("log2 p", p, 380.7007), ("log2 N", Nprod, 251.2721), ("log2 N'", Np, 251.2725)]
fig_out = []
for name, val, pr, kind in FIGURES:
    fig_out.append({"figure": name, "printed_log2": pr, "kind": kind,
                    "exact_log2": round(_log2(val), 6), "safe": _vs_pow2(val, pr, kind == "up")})
for name, val, pr in VALUES:
    fig_out.append({"figure": name, "printed_log2": pr, "kind": "value (4 dp)",
                    "exact_log2": round(_log2(val), 6), "safe": abs(_log2(val) - pr) <= 5e-5})
# plain-number figures quoted in proof.md
plain = {
    "per-step itemised 2,108 <= budget 2,600": 1532 + 24 + 192 + 192 + 40 + 64 + 48 + 16 == 2108 <= 2600,
    "setup itemised sum is 1,074,845,856 < 2^31": setup_sum == 1074845856 and setup_sum < (1 << 31),
    "table 4096 x 1,033 = 4,231,168": setup_items["table_4096_entries"] == 4231168,
    "walk starts 1,070,596,096 <= 2^30": M * 1021 == 1070596096 <= (1 << 30),
    "Keccak round 300, permutation 7,200": sum(KECCAK_ROUND.values()) == 300 and KECCAK_PERM == 7200,
    "sample itemised 7,581 <= budget 8,192": SAMPLE_ITEMISED == 7581 <= SAMPLE_BUDGET,
    "seed itemised 47,744 <= budget 65,536": SEED_ITEMISED == 47744 <= SEED_BUDGET,
    "one 256-bit seed in a 64-word draw": SEED_BITS == 256 and SEED_WORDS == 64 and SEED_BUDGET == (1 << 16),
    "per record itemised 21,274 <= budget 25,212": RECORD_ITEMISED == 21274 <= PER_RECORD,
    "per record 25,212": PER_RECORD == 25212,
    "radix key 4,608 bits = 288 x 16": SORT_PASSES * 16 == 4608,
    "modmul itemised 112 + 120 + 96 + 40 + 40 = 408 <= 512": 112 + 120 + 96 + 40 + 40 == 408 <= MODMUL_R,
    "22 exponent-table mults <= 1 unit": 22 * MODMUL_R <= 20000,
    "samples exact 2,828,231,406,414": sample_sites == 2828231406414,
    "restart count exact 1,414,114,650,534": restart_count == 1414114650534,
    "report cap 2^54 is 2^8.31 above the mean": Fraction(1 << 54) >= Fraction(2) ** Fraction(831, 100) * mean_reports,
    "waste rate 3/(12R) + 200/R^2 <= 7.30e-5": waste_rate <= Fraction(73, 1000000),
    "mean waste <= 7.4e-5 S": mean_waste / S <= Fraction(74, 1000000),
    "cycle tail = mean waste / 0.008 S <= 0.0093": cycle_tail == mean_waste / (Fraction(8, 1000) * S) and cycle_tail <= Fraction(93, 10000),
    "model at S >= 0.408": model >= 0.408,
    "success >= 0.3957": success >= Fraction(3957, 10000),
    "pairing 64 x 4.44 + 6 x 1.44 + 8,628 <= 8,921": Fraction(6400 * 444 + 600 * 144, 10000) + 8628 <= 8921,
    "memory peak = 5 x 2^61 + 2^58 + 2,180,520,704": mem_peak == 5 * (1 << 61) + (1 << 58) + 2180520704,
    "memory peak < 2^64": mem_peak < (1 << 64),
    "records + both index buffers = 656 x 2^54 < 2^64": 656 * (1 << 54) < (1 << 64),
    "slack 125.87 - log2 W >= 0.0087": _vs_pow2(W, "125.8613", True),
}
# printed-operand chains: recomputed from the rounded figures as the reader would (60-digit decimal)
from decimal import Decimal, getcontext
getcontext().prec = 60
def _p2(x): return Decimal(2) ** Decimal(x)
chains = {
    "2 x (1025 + 2^20 + 2^40.37) <= 2^41.39": 2 * (1025 + _p2(20) + _p2("40.37")) <= _p2("41.39"),
    "2^41.39 x 8,192 / 20,000 <= 2^40.11": _p2("41.39") * 8192 / 20000 <= _p2("40.11"),
    "2^40.37 x 1,021 <= 2^50.37": _p2("40.37") * 1021 <= _p2("50.37"),
    "(2^54 + 1) x 25,212 + 2^28.17 <= 2^68.63": (_p2(54) + 1) * 25212 + _p2("28.17") <= _p2("68.63"),
    "2^68.63 / 20,000 <= 2^54.35": _p2("68.63") / 20000 <= _p2("54.35"),
    "(2^54 + 1) x 1,040 <= 2^64.03": (_p2(54) + 1) * 1040 <= _p2("64.03"),
    "2^31 + 2^54.35 + 2^50.37 + 2^64.03 + 2^40.11 + 2^1.72 <= 2^64.04":
        _p2(31) + _p2("54.35") + _p2("50.37") + _p2("64.03") + _p2("40.11") + _p2("1.72") <= _p2("64.04"),
    "1.13 (1.01 x 2^125.6707 + 2^105.33) + 2^64.04 <= 2^125.8614":
        Decimal("1.13") * (Decimal("1.01") * _p2("125.6707") + _p2("105.33")) + _p2("64.04") <= _p2("125.8614"),
    "1.01 x 2^125.6707 + 2^105.33 <= 2^125.6851": Decimal("1.01") * _p2("125.6707") + _p2("105.33") <= _p2("125.6851"),
    "2^16 instructions / 20,000 <= 2^1.72": _p2(16) / 20000 <= _p2("1.72"),
    "2^54 + 2^20 + 2^40.37 + 1 <= 2^54.02": _p2(54) + _p2(20) + _p2("40.37") + 1 <= _p2("54.02"),
    "2^54.02 x 2^-57.7 <= 2^-3.67": _p2("54.02") * _p2("-57.7") <= _p2("-3.67"),
    "2^-3.67 / 2^31.38 <= 2^-35.05": _p2("-3.67") / _p2("31.38") <= _p2("-35.05"),
    "2^45.69 / 2^54 <= 2^-8.31": _p2("45.69") / _p2(54) <= _p2("-8.31"),
    "2^54 / 2^45.69 > 2^8.3": _p2(54) / _p2("45.69") > _p2("8.3"),
}
self_check = {"figures_checked": len(fig_out), "figures_safe": sum(f["safe"] for f in fig_out),
              "unsafe_figures": [f for f in fig_out if not f["safe"]],
              "plain_checks": len(plain), "plain_failed": [k for k, v in plain.items() if not v],
              "printed_operand_chains": len(chains), "chains_failed": [k for k, v in chains.items() if not v]}
self_check["all_safe"] = not (self_check["unsafe_figures"] or self_check["plain_failed"] or self_check["chains_failed"])

partC = {"log2_p": round(_log2(p), 4), "log2_N": round(_log2(Nprod), 4),
         "log2_N_prime": round(_log2(Np), 4), "ln_1_over_0.61_upper_float": round(float(LN), 18),
         "log2_S0_P039": round(_log2(S0), 6), "log2_S": round(_log2(S), 6),
         "model_success_at_S": round(model, 6),
         "log2_L_max": round(_log2(Lmax), 6),
         "log2_E": round(_log2(E), 6),
         "restart_count_le": restart_count, "log2_restart_count": round(_log2(restart_count), 6),
         "setup_items_units": setup_items, "setup_itemised_sum": setup_sum, "setup_declared_log2": 31,
         "instructions_per_record": PER_RECORD,
         "log2_dp_handle_term": round(_log2(dp_handle), 6),
         "log2_restarts_term": round(_log2(restarts), 6),
         "solves_per_record": 2,
         "log2_solves_term": round(_log2(solves), 6),
         "keccak_round_instructions": KECCAK_ROUND, "keccak_permutation_instructions": KECCAK_PERM,
         "sampler_items": SAMPLE_ITEMS, "sampler_itemised": SAMPLE_ITEMISED,
         "sampler_instructions_per_sample": SAMPLE_BUDGET,
         "sampler_draw_bits": DRAW_BITS,
         "samples_le": sample_sites, "log2_samples": round(_log2(sample_sites), 6),
         "log2_sampler_term": round(_log2(sampler), 6),
         "fresh_seed": {"bits": SEED_BITS, "words": SEED_WORDS, "items": SEED_ITEMS,
                        "itemised_instructions": SEED_ITEMISED, "budget_instructions": SEED_BUDGET,
                        "units": round(float(seed_draw), 6), "log2_units": round(_log2(seed_draw), 6),
                        "source": "platform randomness source at start, one draw, no rejection"},
         "randomness_assumption": "H3: one fresh uniform 256-bit seed; the Keccak-f[1600] counter-mode tape under it is computationally indistinguishable from independent uniform draws at the claimed cost, and the platform seed is uniform; the success bound is conditional on H3 (no statistical-distance claim is made)",
         "log2_mean_reports": round(_log2(mean_reports), 6),
         "log2_report_cap_tail_markov": round(_log2(report_tail), 6),
         "log2_cycle_waste_tail_markov": round(_log2(cycle_tail), 6),
         "log2_mean_cycle_waste": round(_log2(mean_waste), 6),
         "success_probability_bound": round(float(success), 6),
         "log2_expected_abandonments": round(_log2(abandon_mean), 6),
         "log2_markov_abandonment_tail": round(_log2(abandon_markov), 6),
         "log2_ties_over_run": round(_log2(ties), 6),
         "class_key_bits": 4608, "class_key_limbs64": 72, "record_bytes": 640,
         "sort_passes": SORT_PASSES, "sort_digit_bits": 16,
         "record_itemised_instructions": RECORD_ITEMISED,
         "memory_peak_lt_2_pow_64": mem_peak < (1 << 64),
         "report_cap": 1 << 54, "device_address_bits": 64,
         "log2_additive_total": round(_log2(Tadd), 6),
         "additive_total_units_ceil": math.ceil(Tadd),
         "c_per_step": float(c), "log2_W": round(lw, 6),
         "claim_log2": float(claim),
         "slack_bits": round(float(claim) - lw, 6),
         "exact_W_to_100_le_2_pow_12587": exact_ok,
         "exact_W_to_100_gt_2_pow_12586": exact_tight,
         "memory_bytes": mem_peak,
         "memory_bytes_exact": str(mem_peak),
         "memory_allowance_bytes": mem_walk + mem_mult + mem_digit + mem_misc,
         "memory_records_bytes": mem_records, "memory_index_buffers_bytes": mem_idx,
         "memory_bytes_log2": round(_log2(mem_peak), 6),
         "ledger_self_check": self_check,
         "ledger_figures": fig_out,
         "note": "All fields are exact-rational computations; only the printed log2 values are rounded. ledger_self_check tests every log2 figure quoted in proof.md against its exact value in the safe direction (exact integer powers), every plain-number claim, and every chain recomputed from the printed operands; the two exact_ fields are the integer certificate. No statistical-distance figure is computed or claimed: the randomness premise is H3, a scoped computational assumption."}
out["partC_exact_ledger"] = partC
out["seconds"] = round(time.time() - T0, 2)
print(json.dumps(out, indent=1))
