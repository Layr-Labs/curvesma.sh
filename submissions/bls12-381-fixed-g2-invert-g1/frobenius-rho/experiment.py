"""Frobenius-class Pollard rho: production identity checks and a toy-field run.

Part A checks, on BLS12-381 itself, the key identities that let the attack pick
a class-canonical Frobenius conjugate from four Fp multiplications.
Part B runs the complete walk (kappa-trick, 2-cycle escape, checkpoint cycle
escape, distinguished points, collision solving) in the order-r' subgroup of
F_{q^12}^* for a small q with the same tower shape as BLS12-381, and measures
the number of steps to the first useful class collision. Toy evidence only:
it does not validate production cost bounds.
"""
import json, math, random, sys, time

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

results = []; Tc_list = []; total_list = []; stats = {"two_cycles": 0, "long_cycles": 0, "useless": 0, "abandoned": 0}
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

sq = math.sqrt(N)
b40 = math.sqrt(2 * N * math.log(1 / 0.6) * R / (R - 1))
out = {"partA_production_BLS12_381": partA,
       "partB_toy": {"q": q, "r_prime": rp, "classes_N": N, "R": R, "dp_bits": dp_bits, "walks": walks,
                     "trials_completed": len(results), "all_logs_recovered_and_verified": all(results),
                     "first_useful_class_collision_mean_over_sqrtN": round(sum(Tc_list) / len(Tc_list) / sq, 4) if Tc_list else None,
                     "random_mapping_prediction_mean_over_sqrtN": round(math.sqrt(math.pi / 2), 4),
                     "fraction_Tc_le_budget40": round(sum(t <= b40 for t in Tc_list) / max(1, len(Tc_list)), 3),
                     "budget40_over_sqrtN": round(b40 / sq, 4),
                     "total_steps_incl_DP_delay_mean_over_sqrtN": round(sum(total_list) / max(1, len(total_list)) / sq, 4),
                     "events": stats,
                     "two_cycles_per_step": stats["two_cycles"] / max(1, sum(total_list)),
                     "predicted_two_cycles_per_step": 1 / (12 * R)},
       "Tc_over_sqrtN": [round(t / sq, 4) for t in Tc_list] if cfg.get("dump") else None,
       "seconds": round(time.time() - T0, 2),
       "scope": "Part A checks exact identities on real BLS12-381 GT elements. Part B is a toy field; it does not validate production cost bounds."}
print(json.dumps(out, indent=1))
