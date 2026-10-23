"""Mixed squaring/multiplication Frobenius-class Pollard rho: identity checks, a
toy-field run, and the cost ledger.

Part A checks, on BLS12-381 itself, the identities the attack uses: the
four-multiplication class key, the canonical conjugate, GT = the cyclotomic
subgroup (r | p^4 - p^2 + 1), and that squaring y -> y^2 commutes with Frobenius,
so it is a map on Frobenius classes that doubles (alpha, beta).
Part B runs the complete walk (kappa trick, class-determined step type, 2-cycle
escape, checkpoint cycle escape, distinguished points, collision solving) in the
order-r' subgroup of F_{q^12}^* for a small q with the same tower shape as
BLS12-381, and measures the number of steps to the first useful class collision,
for a control (mix_1024 = 0) and for mixed walks. Toy evidence only: it does not
validate production cost bounds.
Part C recomputes the fruitless-cycle bound and the per-step ledger with both
step kinds, and checks the submitted claim against log2 W.
"""
import json, math, random, sys, time

T0 = time.time()
cfg = json.load(sys.stdin)
q = cfg.get("q", 1048963)
trials = cfg.get("trials", 40)
R = cfg.get("R", 64)
dp_bits = cfg.get("dp_bits", 4)
walks = cfg.get("walks", 4)
seed = cfg.get("seed", 1)
mix_1024 = cfg.get("mix_1024", 0)    # fraction of squaring steps, in 1/1024 (0 = pure multiplication walk)
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
         "r_divides_Phi12_p_so_GT_is_in_the_cyclotomic_subgroup": (p ** 4 - p ** 2 + 1) % r == 0,
         "gt_order_odd_so_squaring_injective_on_the_class": r % 2 == 1,
         "p_order_mod_r_is_12": pow(p, 12, r) == 1 and all(pow(p, k, r) != 1 for k in range(1, 12)),
         "one_plus_p4_plus_p8_zero_mod_r": (1 + pow(p, 4, r) + pow(p, 8, r)) % r == 0}
ok_keys = ok_frob = ok_canon = ok_gt = ok_inv = ok_sq = True
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
    ok_sq &= P.frob(P.mul(y, y)) == P.mul(P.frob(y), P.frob(y))   # squaring commutes with Frobenius
    ok_sq &= P.mul(y, y) != y and P.mul(y, y) != P.one()          # and is not the identity on GT
partA.update({"random_GT_elements": prod_checks, "GT_membership": ok_gt,
              "frobenius_formula_equals_pth_power": ok_frob,
              "four_mult_key_formula_equals_direct_keys_and_distinct": ok_keys,
              "kappa_canonical_rep_identical_over_all_12_conjugates": ok_canon,
              "p6_power_is_inverse": ok_inv})

# ---------- Part B: toy-field run of the complete algorithm ----------
# Toy fields: q = 7 mod 12 near 2^20, xi = 1+u neither a square nor a cube in F_{q^2},
# r' a prime factor of q^4 - q^2 + 1 with ord_{r'}(q) = 12. q is large enough that the
# class key K (minimum of 12 values in [0, q)) has about 16 significant bits, so the toy
# walk uses the PRODUCTION predicates on the raw key: partition h = K mod R (low log2 R
# bits) and "distinguished" = the next dp_bits bits of K are all zero (production: the
# 60 bits above the 10 partition bits).
TOYS = {1048963: 48674953, 1049611: 64668529, 1054003: 35456053}
T = Tower(q)
n = q ** 4 - q * q + 1
rp = TOYS[q]
assert n % rp == 0 and all(pow(b, rp - 1, rp) == 1 for b in (2, 3, 5, 7, 11, 13))
N = (rp - 1) // 12
qk = [pow(q, k, rp) for k in range(12)]
hc = (q ** 12 - 1) // rp
assert pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))
assert T.f2pow(T.xi, (q * q - 1) // 2) != (1, 0) and T.f2pow(T.xi, (q * q - 1) // 3) != (1, 0)
Rbits = R.bit_length() - 1
assert R == 1 << Rbits, "R must be a power of two"
DPMASK = (1 << dp_bits) - 1
STUCK = 40 << dp_bits                      # production: 40/theta steps without a distinguished class
# Mixed walk. The step type is a function of the class key K only, never of a step or
# walk counter, so two walks at the same class always take the same step and stay
# merged. Production reads bits 198..207 of K (10 bits, above the 10 partition bits,
# the 60 distinguished bits and the 128 record-key bits). The toy key has about 16
# significant bits, so a raw slice is not uniform there and the step type comes from a
# 32-bit mix of K instead; the realised squaring fraction is measured and reported.
SQLEV = 64
SQTHR = (mix_1024 * SQLEV) // 1024
def is_sqr(K):
    return (((K * 2654435761) & 0xFFFFFFFF) >> 26) < SQTHR

Dk = [[T.f2pow(T.xi, m * (q ** k - 1) // 6) for m in range(6)] for k in range(12)]
def frobk(y, k):
    if k % 2: y = tuple((c[0], (-c[1]) % q) for c in y)
    return tuple(T.f2mul(c, Dk[k][m]) for m, c in enumerate(y))

ties = [0]
def keyinfo(y):
    """kappa, class key K and canonical conjugate. Production key: the minimal Fp value
    from keys_formula. In the toy a tie among the 12 keys has probability about 66/q per
    step, so the toy breaks it by comparing the full conjugates (production: ~2^-249 total,
    counted as failure). The canonical conjugate is computed here only for measurement and
    for confirming collisions; the walk itself never applies it."""
    kf = T.keys_formula(y)
    m = min(kf)
    ks = [k for k in range(12) if kf[k] == m]
    if len(ks) == 1:
        kap = ks[0]; t = frobk(y, kap)
    else:
        ties[0] += 1
        t, kap = min((frobk(y, k), k) for k in ks)
    return kap, m, t

results = []; Tc_list = []; total_list = []
stats = {"two_cycles": 0, "long_cycles": 0, "useless": 0, "abandoned": 0, "dp_key_only_matches": 0, "sqr_steps": 0, "mul_steps": 0}
dp_reports = 0
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
    dps = {}; visited = {}; steps = 0; Tc = None; found = None
    def new_walk():
        a, b = rng.randrange(rp), rng.randrange(rp)
        return {"y": T.mul(T.pow(g, a), T.pow(z, b)), "a": a, "b": b, "hist": [], "ck": None, "len": 0, "sincedp": 0}
    W = [new_walk() for _ in range(walks)]
    while found is None:
        for w in W:
            y = w["y"]; kap, K, cls = keyinfo(y)
            ca, cb = w["a"] * qk[kap] % rp, w["b"] * qk[kap] % rp     # canonical exponents
            if Tc is None:                                             # measurement only
                prev = visited.get(cls)
                if prev is None: visited[cls] = (ca, cb)
                elif prev != (ca, cb): Tc = steps
            if (K >> Rbits) & DPMASK == 0:                             # distinguished: production predicate on raw K
                dp_reports += 1
                lst = dps.setdefault(K, [])                            # stored by key, confirmed by canonical element
                for (ocls, oa, ob) in lst:
                    if ocls != cls:
                        stats["dp_key_only_matches"] += 1; continue
                    if ob != cb:
                        xr = (oa - ca) * pow(cb - ob, -1, rp) % rp
                        assert T.pow(g, xr) == z
                        found = xr
                    else:
                        stats["useless"] += 1
                    break
                if found is not None: break
                if all(o[0] != cls for o in lst): lst.append((cls, ca, cb))
                w["sincedp"] = 0
            hist = w["hist"]
            if len(hist) >= 2 and hist[-2][0] == K and hist[-2][1] == cls:   # fruitless 2-cycle (production: K equality)
                stats["two_cycles"] += 1
                (k1, c1, y1, a1, b1), (k2, c2, y2, a2, b2) = hist[-1], (K, cls, y, w["a"], w["b"])
                (_, _, yy, aa, bb) = min(((k1, c1, y1, a1, b1), (k2, c2, y2, a2, b2)), key=lambda t: (t[0], t[1]))
                w["y"] = T.mul(yy, yy); w["a"] = 2 * aa % rp; w["b"] = 2 * bb % rp
                w["hist"] = []; w["ck"] = None; steps += 1; continue
            if w["ck"] is not None and w["ck"][0] == cls and w["ck"][1] < w["len"]:   # longer cycle
                stats["long_cycles"] += 1
                # walk once more around the cycle to find its minimal class, then square it
                best = ((K, cls), y, w["a"], w["b"]); yy, aa, bb = y, w["a"], w["b"]
                for _ in range(w["len"] - w["ck"][1]):
                    kk, KK, _c = keyinfo(yy); hh = KK & (R - 1); idx = (12 - kk) % 12
                    yy = T.mul(yy, Tab[hh][idx]); aa = (aa + At[hh][idx]) % rp; bb = (bb + Bt[hh][idx]) % rp
                    steps += 1
                    _k2, K2, c2 = keyinfo(yy)
                    if (K2, c2) < best[0]: best = ((K2, c2), yy, aa, bb)
                _, yy, aa, bb = best
                w["y"] = T.mul(yy, yy); w["a"] = 2 * aa % rp; w["b"] = 2 * bb % rp
                w["hist"] = []; w["ck"] = None; steps += 1; continue
            if w["len"] % 64 == 0: w["ck"] = (cls, w["len"])
            hist.append((K, cls, y, w["a"], w["b"]))
            if len(hist) > 2: hist.pop(0)
            if is_sqr(K):                                              # class-determined squaring step
                stats["sqr_steps"] += 1
                w["y"] = T.mul(y, y); w["a"] = (2 * w["a"]) % rp; w["b"] = (2 * w["b"]) % rp
            else:                                                      # dense multiplication step
                stats["mul_steps"] += 1
                h = K & (R - 1); idx = (12 - kap) % 12                 # production partition on raw K
                w["y"] = T.mul(y, Tab[h][idx]); w["a"] = (w["a"] + At[h][idx]) % rp; w["b"] = (w["b"] + Bt[h][idx]) % rp
            w["len"] += 1; w["sincedp"] += 1; steps += 1
            if w["sincedp"] > STUCK:
                stats["abandoned"] += 1; w.update(new_walk())
    assert found == xs % rp or T.pow(g, found) == z
    results.append(True); Tc_list.append(Tc); total_list.append(steps)

# ---------- Part C: exact bound arithmetic used in proof.md (no randomness) ----------
def red12(f):  # x^f in Z[x]/(x^4 - x^2 + 1) (Phi_12), basis 1, x, x^2, x^3
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
        for s in roots:
            w = (v[0] + s[0], v[1] + s[1], v[2] + s[2], v[3] + s[3]); nd[w] = nd.get(w, 0) + c
    dist = nd; zero_sum[k] = dist.get((0, 0, 0, 0), 0) / 12 ** k
def log2_PL(L, RR):
    """log2 of sum_{d=1}^{min(RR, L//2)} C(RR,d) (d/RR)^L 12^-d (fruitless L-cycle bound, proof.md 3.2)."""
    ts = [(math.lgamma(RR + 1) - math.lgamma(d + 1) - math.lgamma(RR - d + 1) + L * math.log(d / RR) - d * math.log(12)) / math.log(2)
          for d in range(1, min(RR, L // 2) + 1)]
    mx = max(ts); return mx + math.log2(sum(2 ** (t - mx) for t in ts))
def log2sum(xs):
    mx = max(xs); return mx + math.log2(sum(2 ** (x - mx) for x in xs))
cyc = {}
for RR in (64, 256, 1024):
    short = log2sum([log2_PL(L, RR) for L in range(3, 65)])
    mid = log2sum([log2_PL(L, RR) for L in range(65, 2 * RR + 1)])
    cyc[RR] = {"rate_L3_to_64_times_R2": round(2 ** short * RR * RR, 4),
               "log2_rate_L65_to_2R": round(mid, 1)}
rr = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
lgN = math.log2((rr - 1) // 12)
Rp = 1024
lgS0 = 0.5 * (1 + lgN + math.log2(Rp / (Rp - 1)) + math.log2(math.log(1 / 0.61)))
lgS = lgS0 + math.log2(1.03)
model = 1 - math.exp(-(2 ** (2 * (lgS - lgS0))) * math.log(1 / 0.61))
over = 1 + 3 / (12 * Rp) + 200 / Rp ** 2

# ---- per-step ledger, proof.md section 5 (pinned blst f7a8806, C = 20,000) ----
IR_FP_MUL, IR_FP12_MUL, C = 383, 22996, 20000       # measured, tracks/docs/COST_MODEL.md
ADD = 40                                            # Fp / double-width add-sub upper bound
DENSE_ADD = 30 * ADD + 40 * ADD                     # mul_fp12 add block
CYC_ADD = 42 * ADD + 18 * ADD                       # cyclotomic_sqr_fp12 add block
KEY_COUNTED = 4 * IR_FP_MUL + 9 * ADD + 300 + 40 + 60 + 30
KEY_BUDGET = 2700                                   # 2,322 counted, the frontier's 16% margin
F_MIX = mix_1024 / 1024.0                           # squaring fraction the ledger charges
# Measured mixing inflation (proof.md section 4): the toy mixed walk's fraction within
# 1.03 S0 is 0.372 against 0.418 for the multiplication-only control in the same
# harness (600 trials each), so it needs about 8% more steps for the same success.
# The ledger charges 10% more steps (not a per-step cost).
KAPPA = 1.10 if F_MIX else 1.0
def step_costs(ir_redc=200):
    """(c_mul, c_sqr) in reference units for one redc_mont_384 instruction count.

    cyclotomic_sqr_fp12 = 3 x sqr_fp4 = 9 x sqr_fp2x2 + 12 x redc_mont_384 + add block;
    mul_fp12            = 24 x mul_fp2x2 + 12 x redc_mont_384 + add block, and the
    organizer measures that dense multiplication at 22,996 instructions. sqr_fp2x2 is
    charged at the mul_fp2x2 rate (a square is cheaper), and the squaring gets the same
    16% budget margin as the key block.
    """
    per_fp2 = (IR_FP12_MUL - 12 * ir_redc - DENSE_ADD) / 24.0
    gs = 9 * per_fp2 + 12 * ir_redc + CYC_ADD
    return 1 + 2800 / C, (KEY_BUDGET + 60 + gs * 1.156) / C
c_mul, c_sqr = step_costs()
def log2_W(f, kappa=None):
    c = (1 - f) * c_mul + f * c_sqr
    k = KAPPA if kappa is None else kappa
    return math.log2(c * k * over * 2 ** (lgS - 100) + c * 2 ** (84.3 - 100)
                     + 2 ** (62 - 100)) + 100
lgW = log2_W(F_MIX)
lgcap = math.log2(KAPPA * over * 2 ** lgS)       # walk-step cap, proof.md sections 3 and 5
mu_rep = 2 ** (lgcap - 60)                        # mean distinguished reports under H1
a_rep = 2 ** 66 / mu_rep                          # C_rep / mu
cher = a_rep * math.log(a_rep) - a_rep + 1
CLAIM = 125.767 if F_MIX else 125.861
partC = {"zero_sum_prob_k_unit_12th_roots": {k: round(zero_sum[k], 6) for k in range(2, 11)},
         "zero_sum_max_k_ge_2": round(max(zero_sum[k] for k in range(2, 11)), 6),
         "fruitless_cycle_bounds": cyc,
         "log2_N": round(lgN, 4), "log2_S0": round(lgS0, 4), "log2_S": round(lgS, 4),
         "model_success_at_S": round(model, 4),
         "walk_step_cap_log2": round(lgcap, 5),
         "walk_step_cap_over_S0": round(2 ** (lgcap - lgS0), 4),
         "walk_step_cap_over_S": round(2 ** (lgcap - lgS), 5),
         "report_mean_log2": round(lgcap - 60, 5),
         "report_cap_ratio_a": round(a_rep, 4),
         "chernoff_exponent_log2": round(math.log2(cher * mu_rep), 2),
         "per_step_ledger": {"key_block_counted": KEY_COUNTED, "key_block_budget": KEY_BUDGET,
                             "granger_scott_squaring_counted": round(step_costs(200)[1] * C - KEY_BUDGET - 60, 0),
                             "c_mul": round(c_mul, 5), "c_sqr": round(c_sqr, 5),
                             "c_sqr_at_redc_120": round(step_costs(120)[1], 5),
                             "c_sqr_at_redc_300": round(step_costs(300)[1], 5),
                             "fruitless_cycle_factor": round(over, 6)},
         "mean_step_cost_by_f": {str(f): round((1 - f) * c_mul + f * c_sqr, 5)
                                 for f in (0.0, 0.25, 0.3333, 0.5)},
         "log2_W_by_f_no_inflation": {str(f): round(log2_W(f, 1.0), 5) for f in (0.0, 0.25, 0.3333, 0.5)},
         "log2_W_by_f_charged": {str(f): round(log2_W(f, 1.0 if f == 0 else KAPPA), 5)
                                 for f in (0.0, 0.25, 0.3333, 0.5)},
         "mixing_step_inflation_charged": KAPPA,
         "mixing_step_inflation_measured": 1.08,
         "charged_squaring_fraction": round(F_MIX, 4),
         "log2_W": round(lgW, 5), "claim": CLAIM,
         "claim_ge_log2_W": CLAIM >= lgW}

sq = math.sqrt(N)
b39 = math.sqrt(2 * N * math.log(1 / 0.61) * R / (R - 1))
tot = max(1, sum(total_list))
out = {"partA_production_BLS12_381": partA, "partC_bound_arithmetic": partC,
       "partB_toy": {"q": q, "r_prime": rp, "classes_N": N, "R": R, "dp_bits": dp_bits, "walks": walks,
                     "partition_and_distinguished_predicates": "raw class key K: h = K mod R, distinguished iff bits log2(R)..log2(R)+dp_bits-1 of K are zero (production rule)",
                     "mix_1024": mix_1024,
                     "step_type_rule": "squaring iff the class key selects it (production: bits 198..207 of K below 1024*f; toy: 6 bits of a 32-bit mix of K, measured), so the step is a function of the class and merged walks stay merged",
                     "squaring_steps": stats["sqr_steps"], "multiply_steps": stats["mul_steps"],
                     "squaring_fraction_measured": round(stats["sqr_steps"] / max(1, stats["sqr_steps"] + stats["mul_steps"]), 4),
                     "trials_completed": len(results), "all_logs_recovered_and_verified": all(results),
                     "first_useful_class_collision_mean_over_sqrtN": round(sum(Tc_list) / len(Tc_list) / sq, 4) if Tc_list else None,
                     "random_mapping_prediction_mean_over_sqrtN": round(math.sqrt(math.pi / 2 * R / (R - 1)), 4),
                     "budget39_over_sqrtN": round(b39 / sq, 4),
                     "fraction_Tc_le_budget39": round(sum(t <= b39 for t in Tc_list) / max(1, len(Tc_list)), 3),
                     "fraction_Tc_le_1.03_budget39": round(sum(t <= 1.03 * b39 for t in Tc_list) / max(1, len(Tc_list)), 3),
                     "distinguished_fraction_measured": round(dp_reports / tot, 5),
                     "distinguished_fraction_nominal": 2.0 ** -dp_bits,
                     "total_steps_incl_DP_delay_mean_over_sqrtN": round(tot / max(1, len(total_list)) / sq, 4),
                     "events": stats, "key_ties_toy_only": ties[0],
                     "two_cycles_per_step": stats["two_cycles"] / tot,
                     "predicted_two_cycles_per_step": 1 / (12 * R),
                     "long_cycles_per_step": stats["long_cycles"] / tot,
                     "long_cycle_bound_per_step": 1 / R ** 2},
       "Tc_over_sqrtN": [round(t / sq, 4) for t in Tc_list] if cfg.get("dump") else None,
       "seconds": round(time.time() - T0, 2),
       "scope": "Part A checks exact identities on real BLS12-381 GT elements. Part B is a toy field; it does not validate production cost bounds."}
print(json.dumps(out, indent=1))
