"""Frobenius-class Pollard rho: production identity checks and a toy-field run.

Part A checks, on BLS12-381 itself, the key identities that let the attack pick
a class-canonical Frobenius conjugate from four Fp multiplications, including that
the two-subtraction sign-set minimum equals the twelve-way argmin and class key.
Every check is computed on real production elements; none is a pre-set flag.

Part B runs the complete walk (kappa-trick, 2-cycle escape, checkpoint cycle
escape, distinguished points, offline matching) in the order-r' subgroup of
F_{q^12}^* for three small q with the same tower shape as BLS12-381 (R = 64,
256, 1024). Every trial runs its full fixed step budget; the fraction of trials
with a first useful class collision within the budget and the fraction whose
logarithm the offline matcher recovers and verifies are reported separately, each
with its binomial standard error. One invocation prints one row per field. Toy
evidence only: it does not validate production cost bounds or production key
statistics.

Part A2 measures, on real BLS12-381 GT elements, the frequency of each of the low
70 bits of the class key K (partition bits and distinguished bits).

Part B keeps the production order inside a step: the 2-cycle check and then the
checkpoint-cycle check run BEFORE the distinguished test, so a walk that comes back
to a distinguished class in a fruitless cycle appends nothing, and a report that
equals this walk's last appended record (same key and same canonical exponents) is
not stored. The per-walk report cap is scaled from R_w = 2^47 by the same ratio to
the healthy per-walk mean; a slot that reaches it is stopped. The counter lives in the
slot record, so a stuck-walk restart carries it over instead of resetting it, and the
log bound M x cap is a bound on a slot's lifetime appends. The log length is
reported against its cap, a stress pass with a tiny cap exercises the cap rule, and a
second stress pass forces stuck-walk restarts to exercise the persistent counter.

Part B's recovery path is the production matcher, not an online index: the walk phase
runs every trial to its full fixed step budget, then the collected report log is
sorted by the FULL 80-byte record (128-bit key T, then the canonical exponents A and B),
so identical records are adjacent. Duplicate removal is one adjacent full-record
comparison per record. Each maximal equal-key run of d distinct records after that
dedup is the only place pairwise work occurs: every one of its d(d-1)/2 pairs is
enumerated and charged, the canonical element g^A z^B is recomputed for each record a
charged pair touches, and the logarithm is solved from the first pair whose elements are
equal and whose exponent pairs differ, accepted only when exponentiation verifies
g^x = z. A hard pair budget (production P = 2^41, scaled to toy size) aborts the match
as a failure when the total pair count would exceed it. Genuine repeats are charged as
pairs like any other: the matcher draws no line between a false key match and a repeat.
A measurement-only class index (a toy convenience, never read by the matcher or by a
success flag) timestamps each trial's first useful class collision, so the first-collision
fraction and the recovered fraction are reported separately.

Part C recomputes, without randomness, the zero-sum probabilities, the
fruitless-cycle bound, the deterministic report-log bound M * R_w, the probability
that the per-walk cap fires on a healthy walk, the memory bound, the worst-case
per-record report-handling charge (append, full-record heapsort, adjacent dedup), the
hard pair budget P and the pair charge it bounds, and the work ledger of proof.md
section 5.
"""
import json, math, random, sys, time

T0 = time.time()
cfg = json.load(sys.stdin)
trials = cfg.get("trials", 200)
walks = cfg.get("walks", 16)
dp_bits_default = cfg.get("dp_bits")
seed = cfg.get("seed", 1)
prod_checks = max(1, cfg.get("production_checks", 2))
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

    def key_min_structured(self, x):
        """kappa and class key K from the sign structure of the twelve keys, with no
        twelve-way scan: the keys are the two sign-symmetric triples +-(a, a*zeta,
        a*zeta-a) and +-(b, b*zeta, b*zeta-b), and the minimum over a triple's six
        residues is min(min of magnitudes, p - max of magnitudes)."""
        q, z = self.q, self.zeta
        a = x[1][0]
        b = (x[1][0] * self.d[1][0] + x[1][1] * self.d[1][1]) % q
        best = None
        for base, off in ((a, 0), (b, 1)):
            mags = [base, base * z % q, (base * z - base) % q]       # triple magnitudes
            jm = min(range(3), key=lambda j: mags[j])
            jM = max(range(3), key=lambda j: mags[j])
            for j, neg in ((jm, False), (jM, True)):
                k = 2 * j + off + (6 if neg else 0)
                val = (-mags[j]) % q if neg else mags[j]
                if best is None or val < best[1]:
                    best = (k, val)
        return best                                                 # (kappa, K)

# ---------- Part A: production BLS12-381 identities ----------
p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
P = Tower(p)
zeta = P.zeta
partA = {"zeta_in_Fp_primitive_6th_root":
         pow(zeta, 6, p) == 1 and pow(zeta, 3, p) != 1 and pow(zeta, 2, p) != 1,
         "p_order_mod_r_is_12": pow(p, 12, r) == 1 and all(pow(p, k, r) != 1 for k in range(1, 12)),
         "one_plus_p4_plus_p8_zero_mod_r": (1 + pow(p, 4, r) + pow(p, 8, r)) % r == 0}
ok_keys = ok_frob = ok_canon = ok_gt = ok_inv = ok_min = True
hcof = (p ** 12 - 1) // r
for _ in range(prod_checks):
    x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
    y = P.pow(x, hcof)                                    # random GT element
    ok_gt &= P.pow(y, r) == P.one() and y != P.one()
    ok_frob &= P.frob(y) == P.pow(y, p)                   # Frobenius formula = p-th power
    kd, kf = P.keys_direct(y), P.keys_formula(y)
    ok_keys &= kd == kf and len(set(kd)) == 12
    ks, Ks = P.key_min_structured(y)                          # two-subtraction minimum
    ok_min &= (len(set(kf)) == 12 and ks == min(range(12), key=lambda i: kf[i])
               and Ks == min(kf))
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
              "sign_set_minimum_equals_12way_minimum_and_kappa": ok_min,
              "kappa_canonical_rep_identical_over_all_12_conjugates": ok_canon,
              "p6_power_is_inverse": ok_inv})

# ---------- Part A2: low bits of the class key K on real GT elements ----------
# Pool of 64 random powers of one random GT element; each sample is a product of three pool
# members (y0 to a random exponent mod r, a random element of GT). K = min of the 12 keys from the four-multiplication formula.
def partA2(budget_s):
    t_end = time.time() + budget_s
    x = tuple((rng.randrange(p), rng.randrange(p)) for _ in range(6))
    y0 = P.pow(x, hcof)                                   # a random GT element
    pool = [P.pow(y0, rng.randrange(1, r)) for _ in range(64)]   # 64 random powers of it
    ones = [0] * 70; nsamp = 0; z4 = z6 = 0
    while nsamp < 2000 and len(pool) >= 3 and time.time() < t_end + 1.0:
        y = P.mul(P.mul(rng.choice(pool), rng.choice(pool)), rng.choice(pool))
        if y == P.one(): continue
        K = min(P.keys_formula(y)); nsamp += 1
        for b in range(70):
            ones[b] += (K >> b) & 1
        z4 += ((K >> 10) & 15) == 0; z6 += ((K >> 10) & 63) == 0
        if time.time() > t_end: break
    if not nsamp: return {"samples": 0}
    se = 0.5 / math.sqrt(nsamp)
    dev = [abs(o / nsamp - 0.5) / se for o in ones]
    return {"samples": nsamp, "pool_GT_elements": len(pool), "bits_tested": 70,
            "max_abs_deviation_in_standard_errors": round(max(dev), 2),
            "bit_at_max_deviation": dev.index(max(dev)),
            "mean_abs_deviation_in_standard_errors": round(sum(dev) / 70, 2),
            "expected_mean_abs_deviation_if_uniform": round(math.sqrt(2 / math.pi), 2),
            "fraction_bits_10_13_zero": round(z4 / nsamp, 4), "nominal_2^-4": 2.0 ** -4,
            "fraction_bits_10_15_zero": round(z6 / nsamp, 4), "nominal_2^-6": 2.0 ** -6}
partA2_out = partA2(2.2)
TB0 = time.time()                                  # Part B's time starts here
budgetB = max(1.0, budget - (TB0 - T0))

# ---------- Part B: toy-field run of the complete algorithm, one row per field ----------
# Toy fields: q = 7 mod 12 near 2^20, xi = 1+u neither a square nor a cube in F_{q^2},
# r' a prime factor of q^4 - q^2 + 1 with ord_{r'}(q) = 12. q is large enough that the
# class key K (minimum of 12 values in [0, q)) has about 16 significant bits, so the toy
# walk uses the PRODUCTION predicates on the raw key: partition h = K mod R (low log2 R
# bits) and "distinguished" = the next dp_bits bits of K are all zero (production: the
# 60 bits above the 10 partition bits).
TOYS = {1048963: 48674953, 1049611: 64668529, 1054003: 35456053}
CONFIGS = [(1048963, 64, 4), (1049611, 256, 4), (1054003, 1024, 3)]
# Time slices favour the slower R=1024 field so its sample stays useful.
SLICES = [0.25, 0.25, 0.5]


def binom_err(f, n):
    return math.sqrt(max(0.0, f * (1.0 - f)) / n) if n else None


def run_toy(q, R, dp_bits, trials_target, deadline):
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
    nbits = rp.bit_length()

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

    b39 = math.sqrt(2 * N * math.log(1 / 0.61) * R / (R - 1))
    budget_steps = 1.03 * b39
    theta = 2.0 ** -dp_bits
    # Per-walk report cap, scaled to the toy: production R_w = 2^47 is 2^1.328 = 2.511 times the
    # healthy per-walk mean report count 2^45.672 over S/M = 2^105.67 steps at theta = 2^-60. The toy
    # applies the same ratio to its own per-walk mean, then the hard rule "cap reports -> stop the
    # walk, its slot appends nothing more".
    CAP_RATIO = 2.510954918647871
    cap_walk = max(1, int(math.ceil(CAP_RATIO * (budget_steps / walks) * theta)))
    log_cap_total = walks * cap_walk
    # Hard pair budget of the matcher, scaled to toy size with the production margin. Production:
    # P = 2^41 pairs, 2^36 times the ~2^5 false-key pairs it expects over a 2^67-record log. The
    # toy's key is only a few bits wide, so its worst-case pair count (every pair in one run, about
    # log_cap^2/2) is far above its expected false-key pairs; taking 2^36 times THAT worst case is
    # the conservative version of the same margin and no healthy toy run can approach it. The abort
    # rule is therefore exercised by a separate tight-budget stress pass below, not by the rows.
    pair_budget_toy = 2 ** 36 * (log_cap_total * (log_cap_total - 1) // 2)

    def walk_pass(cap, step_budget, instrument, stuck=None):
        """One trial's walk phase, run to its full fixed step budget: no trial stops at its first
        recovery. The order inside a step is the production order of proof.md section 3: the
        2-cycle check and then the checkpoint-cycle check run BEFORE the distinguished test, so a
        walk that comes back to a distinguished class in a fruitless cycle appends nothing; a
        report is stored only when it differs from this walk's last appended record (same 128-bit
        key AND same two canonical exponents), which drops the identical records a fruitless
        return would otherwise store. `instrument` adds a MEASUREMENT-ONLY class index (a toy
        convenience, never read by the matcher or by a success flag) that timestamps the first
        useful class collision; recovery is the offline matcher below. `stuck` overrides the
        stuck-walk threshold (the restart-stress pass uses a tiny one)."""
        visited = {} if instrument else None
        stuck_steps = STUCK if stuck is None else stuck
        steps = 0; Tc = None; log = []
        st = {"two_cycles": 0, "long_cycles": 0, "abandoned": 0,
              "duplicate_appends_skipped": 0, "cap_stops": 0}

        def new_walk(rep=0):
            """A fresh trajectory in an existing slot. `rep` is the slot's cumulative report
            count, carried over by a restart (proof.md section 3: no restart resets it)."""
            a, b = rng.randrange(rp), rng.randrange(rp)
            return {"y": T.mul(T.pow(g, a), T.pow(z, b)), "a": a, "b": b, "hist": [], "ck": None,
                    "len": 0, "sincedp": 0, "rep": rep, "dead": False, "last": None}
        W = [new_walk() for _ in range(walks)]
        while steps < step_budget:
            live = [w for w in W if not w["dead"]]
            if not live: break
            for w in live:
                if steps >= step_budget: break
                y = w["y"]; kap, K, cls = keyinfo(y)
                ca, cb = w["a"] * qk[kap] % rp, w["b"] * qk[kap] % rp     # canonical exponents
                hist = w["hist"]
                if len(hist) >= 2 and hist[-2][0] == K and hist[-2][1] == cls:
                    st["two_cycles"] += 1                                  # 2-cycle check, BEFORE any append
                    (k1, c1, y1, a1, b1), (k2, c2, y2, a2, b2) = hist[-1], (K, cls, y, w["a"], w["b"])
                    (_, _, yy, aa, bb) = min(((k1, c1, y1, a1, b1), (k2, c2, y2, a2, b2)), key=lambda t: (t[0], t[1]))
                    w["y"] = T.mul(yy, yy); w["a"] = 2 * aa % rp; w["b"] = 2 * bb % rp
                    w["hist"] = []; w["ck"] = None; steps += 1; continue
                if w["ck"] is not None and w["ck"][0] == cls and w["ck"][1] < w["len"]:
                    st["long_cycles"] += 1                                 # cycle of length <= 64, BEFORE any append
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
                if visited is not None and Tc is None:                    # measurement only, not the recovery path
                    prev = visited.get(cls)
                    if prev is None: visited[cls] = (ca, cb)
                    elif prev != (ca, cb): Tc = steps
                if (K >> Rbits) & DPMASK == 0:                             # distinguished: production predicate on raw K
                    rec = (K >> (Rbits + dp_bits), ca, cb)                  # 80-byte record: truncated key, then both exponents
                    if rec != w["last"]:                                   # append rule: no exact repeat of this walk's last record
                        log.append(rec); w["last"] = rec; w["rep"] += 1
                        w["sincedp"] = 0
                        if w["rep"] >= cap:                                # per-slot report cap: stop this walk
                            st["cap_stops"] += 1; w["dead"] = True
                            continue
                    else:
                        st["duplicate_appends_skipped"] += 1
                        w["sincedp"] = 0
                if w["len"] % 64 == 0: w["ck"] = (cls, w["len"])
                hist.append((K, cls, y, w["a"], w["b"]))
                if len(hist) > 2: hist.pop(0)
                h = K & (R - 1); idx = (12 - kap) % 12                     # production partition on raw K
                w["y"] = T.mul(y, Tab[h][idx]); w["a"] = (w["a"] + At[h][idx]) % rp; w["b"] = (w["b"] + Bt[h][idx]) % rp
                w["len"] += 1; w["sincedp"] += 1; steps += 1
                if w["sincedp"] > stuck_steps:
                    st["abandoned"] += 1
                    w.update(new_walk(w["rep"]))      # restart: fresh trajectory, same slot counter
        # The log length must equal the sum of the slots' persistent counters: every append
        # increments exactly one slot counter, and no restart resets one.
        sum_rep = sum(w["rep"] for w in W); max_rep = max([w["rep"] for w in W] + [0])
        return {"steps": steps, "Tc": Tc, "log": log, "events": st,
                "sum_slot_reports": sum_rep, "max_slot_reports": max_rep}

    def elem(A, B, gp, zp):
        """g^A z^B from precomputed powers of g and z (one multiplication per set bit)."""
        acc = T.one(); j = 0
        while (A or B) and j < len(gp):
            if A & 1: acc = T.mul(acc, gp[j])
            if B & 1: acc = T.mul(acc, zp[j])
            A >>= 1; B >>= 1; j += 1
        return acc

    def match_and_verify(log, pair_budget, g, z, gp, zp):
        """Production-shaped offline matcher (proof.md section 3) on the collected log: sort by the
        FULL 80-byte record (key T first, then canonical A and B), so identical records are adjacent;
        scan maximal equal-key runs; drop an exact duplicate of the previous record as the same record
        (one adjacent full-record comparison per record); enumerate every pair of the remaining
        distinct records of the run, charging each pair against the hard pair budget and recomputing
        the canonical element g^A z^B of each record a charged pair touches; solve
        x = (A1-A2)/(B2-B1) from the first pair whose elements are equal and whose exponent pairs
        differ; accept x only when exponentiation verifies g^x = z. A total pair count above the
        budget aborts the match as a failure (proof.md section 3: Markov over the false-key pairs). A
        genuine repeat is charged as a pair like any other: no record is exempt."""
        recs = sorted(log)                 # full record: key T first, then canonical A, then B
        pairs = 0; recomputations = 0; identical = 0; equal_key_runs = 0
        equal_element_pairs = 0; useless_pairs = 0; recovered = None; cap_hit = False
        i = 0
        while i < len(recs) and recovered is None and not cap_hit:
            j = i + 1
            while j < len(recs) and recs[j][0] == recs[i][0]: j += 1
            run = recs[i:j]
            if len(run) >= 2:
                equal_key_runs += 1
                distinct = []
                for rec in run:
                    if distinct and rec == distinct[-1]: identical += 1   # adjacent after the full-record sort
                    else: distinct.append(rec)
                elems = [None] * len(distinct)
                for u in range(len(distinct)):
                    for v in range(u + 1, len(distinct)):
                        if pairs >= pair_budget:
                            cap_hit = True; break
                        pairs += 1
                        for w in (u, v):
                            if elems[w] is None:
                                elems[w] = elem(distinct[w][1], distinct[w][2], gp, zp)
                                recomputations += 1
                        if elems[u] != elems[v]: continue
                        equal_element_pairs += 1
                        (a1, b1), (a2, b2) = distinct[u][1:], distinct[v][1:]
                        if b1 == b2: useless_pairs += 1; continue
                        xr = (a1 - a2) * pow(b2 - b1, -1, rp) % rp
                        if T.pow(g, xr) == z:
                            recovered = xr
                            break
                    if recovered is not None or cap_hit: break
            i = j
        return {"recovered": recovered, "pairs": pairs, "recomputations": recomputations,
                "identical_records": identical, "equal_key_runs": equal_key_runs,
                "equal_element_pairs": equal_element_pairs, "useless_pairs": useless_pairs,
                "cap_hit": cap_hit}

    Tc_list = []; total_list = []; reports_list = []
    stats = {"two_cycles": 0, "long_cycles": 0, "abandoned": 0,
             "duplicate_appends_skipped": 0, "cap_stops": 0}
    dp_reports = 0; max_log = 0
    n_trials = 0; n_collided = 0; n_recovered = 0; n_collided_no_recovery = 0
    n_cap_hit = 0; identical_total = 0; pairs_max = 0; pairs_total = 0
    recomp_max = 0; recomp_total = 0
    while n_trials < trials_target and time.time() < deadline:
        while True:
            g = T.pow(tuple((rng.randrange(q), rng.randrange(q)) for _ in range(6)), hc)
            if g != T.one(): break
        xs = rng.randrange(rp); z = T.pow(g, xs)
        ah = [rng.randrange(rp) for _ in range(R)]; bh = [rng.randrange(rp) for _ in range(R)]
        # Mh = g^{ah} z^{bh} = g^{ah + xs bh} (z = g^xs). One exponentiation per
        # multiplier from shared powers of g, instead of two; identical values.
        ch = [(ah[i] + xs * bh[i]) % rp for i in range(R)]
        gpow = [g]
        for _ in range(nbits): gpow.append(T.mul(gpow[-1], gpow[-1]))
        zp = [z]
        for _ in range(nbits): zp.append(T.mul(zp[-1], zp[-1]))
        Mh = []
        for c in ch:
            e = c; acc = T.one(); j = 0
            while e:
                if e & 1: acc = T.mul(acc, gpow[j])
                e >>= 1; j += 1
            Mh.append(acc)
        Tab = []; At = []; Bt = []
        for i in range(R):
            row = [Mh[i]]
            for k in range(11): row.append(T.frob(row[-1]))
            Tab.append(row); At.append([ah[i] * qk[k] % rp for k in range(12)]); Bt.append([bh[i] * qk[k] % rp for k in range(12)])
        r = walk_pass(cap_walk, budget_steps, True)          # full fixed step budget, every trial
        m = match_and_verify(r["log"], pair_budget_toy, g, z, gpow, zp)
        # Recovery counts only when the exponentiation check inside the matcher passed.
        assert m["recovered"] is None or T.pow(g, m["recovered"]) == z
        assert len(r["log"]) <= log_cap_total
        n_trials += 1
        if r["Tc"] is not None: n_collided += 1
        if m["recovered"] is not None: n_recovered += 1
        if r["Tc"] is not None and m["recovered"] is None: n_collided_no_recovery += 1
        if m["cap_hit"]: n_cap_hit += 1
        identical_total += m["identical_records"]
        pairs_total += m["pairs"]; pairs_max = max(pairs_max, m["pairs"])
        recomp_total += m["recomputations"]; recomp_max = max(recomp_max, m["recomputations"])
        Tc_list.append(r["Tc"]); total_list.append(r["steps"])
        for key, val in r["events"].items(): stats[key] += val
        dp_reports += len(r["log"]); reports_list.append(len(r["log"]))
        max_log = max(max_log, len(r["log"]))

    # Stress pass: the same walk code with a tiny per-walk cap, so the cap rule and the log bound
    # are exercised on real walk data instead of being asserted.
    stress = None
    restart_stress = None
    pair_stress = None
    if n_trials:
        CAP_SMALL = 4
        r = walk_pass(CAP_SMALL, 3000, False)
        ok1 = len(r["log"]) == r["sum_slot_reports"]
        assert ok1
        stress = {"cap_per_walk": CAP_SMALL, "walks": walks, "log_cap_total": walks * CAP_SMALL,
                  "reports_logged": len(r["log"]), "cap_stops": r["events"]["cap_stops"], "steps": r["steps"],
                  "duplicate_appends_skipped": r["events"]["duplicate_appends_skipped"],
                  "slot_reports_sum": r["sum_slot_reports"],
                  "log_equals_slot_report_sum": ok1,
                  "log_within_cap": len(r["log"]) <= walks * CAP_SMALL}
        # Restart stress: a stuck-walk threshold of 2 steps forces restarts inside a slot. The log
        # length still equals the sum of the slots' counters and stays under M x cap, because the
        # report counter persists across a restart (a counter reset by a restart would let one slot
        # append past its cap and leave the sum short of the log length).
        r2 = walk_pass(CAP_SMALL, 3000, False, stuck=2)
        ok2 = (len(r2["log"]) == r2["sum_slot_reports"]) and len(r2["log"]) <= walks * CAP_SMALL
        assert ok2
        restart_stress = {"cap_per_walk": CAP_SMALL, "stuck_walk_threshold": 2, "walks": walks,
                          "log_cap_total": walks * CAP_SMALL, "restarts": r2["events"]["abandoned"],
                          "reports_logged": len(r2["log"]), "slot_reports_sum": r2["sum_slot_reports"],
                          "max_slot_reports": r2["max_slot_reports"], "cap_stops": r2["events"]["cap_stops"],
                          "log_equals_slot_report_sum": len(r2["log"]) == r2["sum_slot_reports"],
                          "log_within_cap": len(r2["log"]) <= walks * CAP_SMALL}
        # Pair-budget stress: the production matcher on a real collected log, and on a two-record log
        # with one equal-key pair, at pair budgets 1 and 0. The hard pair rule and its abort are
        # exercised on real records instead of being asserted.
        mp = match_and_verify(r["log"], 1, g, z, gpow, zp)
        synth = [(7, 12345 % rp, 67890 % rp), (7, 999 % rp, 111 % rp)]
        m1 = match_and_verify(synth, 1, g, z, gpow, zp)
        m0 = match_and_verify(synth, 0, g, z, gpow, zp)
        pair_stress = {"pair_budget": 1, "reports": len(r["log"]), "pairs_charged": mp["pairs"],
                       "aborted": mp["cap_hit"], "pairs_within_budget": mp["pairs"] <= 1,
                       "two_record_log_pairs_at_budget_1": m1["pairs"],
                       "two_record_log_aborts_at_budget_1": m1["cap_hit"],
                       "two_record_log_pairs_at_budget_0": m0["pairs"],
                       "two_record_log_aborts_at_budget_0": m0["cap_hit"]}

    sq = math.sqrt(N)
    ntrials = max(1, n_trials)
    hits = [t for t in Tc_list if t is not None]
    f39 = sum(t <= b39 for t in hits) / ntrials
    f103 = sum(t <= 1.03 * b39 for t in hits) / ntrials
    frec = n_recovered / ntrials
    tot = max(1, sum(total_list))
    return {"q": q, "r_prime": rp, "classes_N": N, "R": R, "dp_bits": dp_bits,
            "trials_completed": n_trials, "step_budget_1.03_S0": round(budget_steps, 1),
            "first_collision_fraction_within_1.03_S0": round(f103, 3),
            "recovered_and_verified_fraction": round(frec, 3),
            "recovered_and_verified_count": n_recovered,
            "recovered_and_verified_fraction_binomial_err": round(binom_err(frec, ntrials), 3),
            "trials_with_first_collision": n_collided,
            "collision_without_recovery": n_collided_no_recovery,
            "every_recovered_logarithm_verified": n_recovered > 0,
            "first_useful_class_collision_mean_over_sqrtN": round(sum(hits) / max(1, len(hits)) / sq, 4) if hits else None,
            "random_mapping_prediction_mean_over_sqrtN": round(math.sqrt(math.pi / 2 * R / (R - 1)), 4),
            "budget39_over_sqrtN": round(b39 / sq, 4),
            "fraction_Tc_le_budget39": round(f39, 3), "fraction_Tc_le_budget39_binomial_err": round(binom_err(f39, ntrials), 3),
            "fraction_Tc_le_1.03_budget39": round(f103, 3), "fraction_Tc_le_1.03_budget39_binomial_err": round(binom_err(f103, ntrials), 3),
            "distinguished_fraction_measured": round(dp_reports / tot, 5),
            "distinguished_fraction_nominal": 2.0 ** -dp_bits,
            "total_steps_incl_DP_delay_mean_over_sqrtN": round(tot / max(1, len(total_list)) / sq, 4),
            "events": stats, "key_ties_toy_only": ties[0],
            "two_cycles_per_step": stats["two_cycles"] / tot,
            "predicted_two_cycles_per_step": 1 / (12 * R),
            "long_cycles_per_step": stats["long_cycles"] / tot,
            "long_cycle_events": stats["long_cycles"],
            "long_cycle_bound_per_step": 1 / R ** 2,
            "report_cap_per_walk_toy": cap_walk, "report_log_cap_total_toy": log_cap_total,
            "report_log_max_length_toy": max_log, "report_log_within_cap": max_log <= log_cap_total,
            "report_cap_stops": stats["cap_stops"],
            "reports_per_trial_mean": round(sum(reports_list) / max(1, len(reports_list)), 1),
            "per_walk_mean_reports_toy": round((budget_steps / walks) * theta, 1),
            "duplicate_appends_skipped": stats["duplicate_appends_skipped"],
            "pairs_charged_total": pairs_total, "pairs_charged_max_per_trial": pairs_max,
            "pair_budget_toy": pair_budget_toy,
            "pairs_within_budget": pairs_max <= pair_budget_toy and n_cap_hit == 0,
            "pair_budget_aborts": n_cap_hit,
            "element_recomputations_total": recomp_total,
            "element_recomputations_max_per_trial": recomp_max,
            "identical_record_duplicates_dropped": identical_total,
            "pair_stress": pair_stress, "cap_stress": stress, "restart_stress": restart_stress}


partB_rows = []
for i, (q, R, dpb) in enumerate(CONFIGS):
    stop = TB0 + budgetB * (sum(SLICES[:i + 1]) if i < len(SLICES) else 1.0)
    row = run_toy(q, R, dpb, trials, min(stop, TB0 + budgetB))
    partB_rows.append(row)

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
C = 20000
# Ancillary work of section 5. Every count below is a worst-case (hard) count, never a mean.
# Report handling at the deterministic log bound C_rep = M * R_w = 2^20 * 2^47 = 2^67: append,
# full-record heapsort and adjacent dedup at 11,100 instructions per record; the pair charge bounded
# by the hard pair budget P = 2^41 (200 instructions per pair, plus one canonical-element
# recomputation per record a charged pair touches, at most 2P records); and setup/output. Together
# they are bounded by a single 2^66.2-unit term.
M_WALKS, RW, CREP_LG = 2 ** 20, 2 ** 47, 67
CREP = M_WALKS * RW
assert CREP == 2 ** CREP_LG
NSTEPS_MEAN_LG = lgS + math.log2(1.000435)                       # steps incl. cycle handling (E adds < 2^-41)
mu_walk = 2 ** (NSTEPS_MEAN_LG - 60 - 20)                        # healthy reports per walk over S/M steps
a_cap = RW / mu_walk                                             # cap over the healthy per-walk mean
kl_div = a_cap * math.log(a_cap) - a_cap + 1                     # D(a||1): Bernoulli tail exponent
lg_fire = -mu_walk * kl_div / math.log(2)                        # log2 P(healthy walk reaches the cap)
APPEND, DEDUP, SIFT_LEVEL = 1000, 100, 130
levels = 67                                                      # heap construction (<= n) plus extractions
# exact: sum_{k=2}^{2^67} floor(log2 k) = 65 * 2^67 + 69 <= 66 n, so the total is <= 67 n levels
extract_levels = sum(j * 2 ** j for j in range(1, CREP_LG)) + CREP_LG
SORT = 10000                                                     # 67 levels x 130 instructions = 8,710, charged 10,000
assert levels * SIFT_LEVEL <= SORT
assert extract_levels <= (CREP_LG - 1) * CREP and extract_levels + CREP <= levels * CREP
per_record = SORT + DEDUP + APPEND
lg_report = math.log2(CREP * per_record / C)
PAIR_BUDGET_LG, PAIR_INSTR, RECOMP_UNITS = 41, 200, 1021
lg_pair = math.log2(2 ** PAIR_BUDGET_LG * PAIR_INSTR / C)
lg_recomp = math.log2(2 * 2 ** PAIR_BUDGET_LG * RECOMP_UNITS)    # one per record a charged pair touches
false_pairs_lg = 2 * CREP_LG - 1 - 128                           # (2^67)^2/2 pairs x 2^-128 = 2^5 expected
lg_abort = false_pairs_lg - PAIR_BUDGET_LG                       # Markov: P(pairs > P) <= 2^5 / 2^41
assert lg_abort <= -36
SETUP = 2 ** 31
lg_anc = math.log2(2 ** lg_report + 2 ** lg_pair + 2 ** lg_recomp + SETUP)
ANC = 66.2
assert lg_anc <= ANC
# Per-step ledger of proof.md section 5 (instructions; blst_fp_mul = 383 measured).
STEP_LEDGER = {"4_fp_mul": 4 * 383, "B_add_mod_p": 20, "two_sub_mod_p": 40,
               "min3_max3_two_triples": 144, "two_reflections_p_minus_max": 40,
               "final_min_of_4_and_kappa": 74, "load_w1_coefficient": 16,
               "h_t_distinguished_test_table_address": 20, "alpha_beta_add_mod_r": 52,
               "cycle_compares_key_history": 40, "loop_control_output_buffer": 30,
               "three_128bit_step_counters": 12, "stuck_walk_and_step_cap_tests": 6}
STEP_COUNT = sum(STEP_LEDGER.values())
STEP_BUDGET = 2250
assert STEP_COUNT <= STEP_BUDGET
CSTEP = 1 + STEP_BUDGET / C                                      # c = 1 + 2,250 / 20,000
lgW = math.log2(CSTEP * over * 2 ** (lgS - 100) + CSTEP * 2 ** (84.3 - 100) + 2 ** (ANC - 100)) + 100
mem_bytes = CREP * 80 + 2 ** 20 * 1024 + 12288 * 576 + 2 * 12288 * 32
partC = {"zero_sum_prob_k_unit_12th_roots": {k: round(zero_sum[k], 6) for k in range(2, 11)},
         "zero_sum_max_k_ge_2": round(max(zero_sum[k] for k in range(2, 11)), 6),
         "fruitless_cycle_bounds": cyc,
         "log2_N": round(lgN, 4), "log2_S0": round(lgS0, 4), "log2_S": round(lgS, 4),
         "model_success_at_S": round(model, 4),
         "report_log_rule": "M * R_w = 2^20 * 2^47, deterministic: a walk slot that reaches R_w is stopped",
         "append_rule": ("2-cycle and checkpoint-cycle checks run before the report append; a report equal to "
                         "this walk's last appended record (same 128-bit key and same two canonical exponents) "
                         "is not stored; a class fixes the canonical element and not the exponent pair, so a "
                         "revisit with a different pair is the useful collision x = (A1-A2)/(B2-B1); an exact "
                         "duplicate left by a longer cycle is dropped by the adjacent full-record dedup as the "
                         "same record, and every remaining pair in its run is charged like any other"),
         "matcher_rule": ("sort the log by the full 80-byte record (key T, then A, then B), so identical records "
                          "are adjacent; dedup is one adjacent full-record comparison per record; a maximal "
                          "equal-key run of d distinct records after that is the only place pairwise work occurs, "
                          "and all d(d-1)/2 of its pairs are charged at <= 200 instructions each (element compare "
                          "plus the canonical recomputation check) against a hard pair budget P = 2^41; the run "
                          "aborts (failure) if the log's total pair count would exceed P; the first pair whose "
                          "elements are equal with B1 != B2 gives x = (A1-A2)/(B2-B1). A genuine repeat is charged "
                          "as a pair like any other: no record is exempt"),
         "pair_budget_P_log2": PAIR_BUDGET_LG,
         "pair_charge_log2_units": round(lg_pair, 3),
         "element_recomputation_log2_units": round(lg_recomp, 3),
         "expected_false_pairs_log2": false_pairs_lg,
         "abort_prob_markov_log2": lg_abort,
         "report_log_capacity_C_rep_log2": CREP_LG,
         "per_walk_cap_R_w_log2": 47, "mean_reports_per_walk_log2": round(math.log2(mu_walk), 4),
         "cap_over_walk_mean_a": round(a_cap, 4), "kl_divergence_D_a_1": round(kl_div, 4),
         "log2_cap_fire_on_healthy_walk": round(lg_fire, 1),
         "cap_fire_prob_below_2^-300": lg_fire < -300,
         "per_record_instructions": {"append": APPEND, "full_record_heapsort_charged": SORT,
                                     "adjacent_dedup": DEDUP, "total": per_record},
         "heapsort_levels_per_record": levels, "instructions_per_sift_level": SIFT_LEVEL,
         "heapsort_levels_x_instr": levels * SIFT_LEVEL,
         "report_handling_log2_units_at_cap": round(lg_report, 3),
         "distinct_stored_records_log2": CREP_LG,
         "distinct_key_pairs_log2": false_pairs_lg + 128,
         "setup_and_output_log2_units": 31,
         "ancillary_total_log2_units": round(lg_anc, 3), "ancillary_charged_log2_units": ANC,
         "memory_log2_bytes_hard": round(math.log2(mem_bytes), 4), "memory_below_2^73.4": mem_bytes < 2 ** 73.4,
         "counter_widths_ok_128bit": {"steps_since_report_log2": round(math.log2(40 * 2 ** 60), 2),
                                      "per_walk_steps_log2": round(NSTEPS_MEAN_LG - 20, 2),
                                      "global_steps_log2": round(NSTEPS_MEAN_LG, 2), "report_counter_log2": CREP_LG,
                                      "record_address_log2": round(math.log2(CREP * 80), 2),
                                      "per_walk_report_counter_log2": 47},
         "per_step_instructions_counted": STEP_COUNT, "per_step_instruction_budget": STEP_BUDGET,
         "step_cost_c": round(CSTEP, 6),
         "log2_W": round(lgW, 4), "claim": 125.826,
         "claim_ge_log2_W": 125.826 >= lgW}

out = {"partA_production_BLS12_381": partA, "partA2_key_bit_uniformity": partA2_out,
       "partC_bound_arithmetic": partC,
       "partB_toy": {"walks": walks, "seconds": round(budget, 2),
                     "partition_and_distinguished_predicates": "raw class key K: h = K mod R, distinguished iff bits log2(R)..log2(R)+dp_bits-1 of K are zero (production rule)",
                     "step_order_and_append_rule": "2-cycle and checkpoint-cycle checks before the distinguished test; a report equal to this walk's last appended record (same key, same canonical exponents) is not stored; the slot's report counter persists across a stuck-walk restart",
                     "match_method": "offline: sort the report log by the full 80-byte record (key T, then A, then B), so identical records are adjacent; dedup is one adjacent full-record comparison per record; enumerate every pair of the remaining distinct records of each equal-key run, charge each pair against the hard pair budget, recompute g^A z^B for each record a charged pair touches, solve x = (A1-A2)/(B2-B1) from the first pair whose elements are equal and whose exponent pairs differ, accept x only after exponentiation verifies g^x = z (a trial counts as recovered only when that check passes); a total pair count above the budget aborts the match as a failure, and a repeat is charged as a pair like any other",
                     "trials_stop_at_first_recovery": False,
                     "every_recovered_logarithm_verified": all(row["every_recovered_logarithm_verified"] for row in partB_rows if row["trials_completed"]),
                     "rows_with_trials": sum(1 for row in partB_rows if row["trials_completed"]),
                     "recovered_total": sum(row["recovered_and_verified_count"] for row in partB_rows),
                     "collision_without_recovery_total": sum(row["collision_without_recovery"] for row in partB_rows),
                     "rows": partB_rows},
       "seconds": round(time.time() - T0, 2),
       "scope": "Part A checks exact identities on real BLS12-381 GT elements (all computed, none pre-set). Part B is a toy field: its walk runs to a fixed toy step budget and its report log is tiny, so neither the production cost bounds nor the production 128-bit key statistics are validated there."}
print(json.dumps(out, indent=1))
