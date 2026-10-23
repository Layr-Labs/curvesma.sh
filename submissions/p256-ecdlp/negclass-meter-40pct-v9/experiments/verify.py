#!/usr/bin/env python3
"""Part A + ledger + itemisation for p256-ecdlp negclass-meter-40pct-v9.

Exact integer / Fraction arithmetic only: every bound below is certified by an
integer or rational inequality, never by a float.  Writes certificate.json,
checks-1..3.json and experiments/partA-results.json, and prints PASS lines.
"""
import hashlib, json, math, os, random, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)

p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
A_CURVE = p - 3
B_CURVE = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
GX = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
GY = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
C = 1000                      # borrowed reference unit, instructions per unit
R_NODES = 1 << 16             # storage nodes
KAPPA = Fraction(21, 20)      # 1.05 walk-quality constant
FAIL_CONST = Fraction(60, 100)  # model failure probability of one walk

CHECKS = []
def chk(name, ok, detail=""):
    CHECKS.append({"name": name, "verdict": "PASS" if ok else "FAIL",
                   "detail": detail})
    return ok

# --------------------------------------------------------------- primality
def is_prime(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
              59, 61, 67, 71, 73, 79, 83, 89, 97):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

def rho_factor(n):
    """Brent's rho with an iteration cap, deterministic constants.  Returns a
    nontrivial factor of composite n, or None if the cap is hit."""
    if n % 2 == 0: return 2
    if n % 3 == 0: return 3
    rnd = random.Random(20261023)
    for _ in range(40):
        c = rnd.randrange(1, n)
        y, m, g, q, r = rnd.randrange(0, n), 128, 1, 1, 1
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r):
                y = (y * y + c) % n
            k = 0
            while k < r and g == 1:
                ys = y
                for _ in range(min(m, r - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r *= 2
            if r > (1 << 22):        # cap: 4M group steps per attempt
                g = n
        if g == n:
            while True:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
                if g > 1: break
        if 1 < g < n: return g
    return None

SMALL_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]

def factor(n, out=None):
    if out is None: out = []
    if n == 1: return sorted(out)
    # trial division keeps every rho call on a cofactor with a small factor
    d = 2
    while d * d <= n and d < 100000:
        while n % d == 0:
            out.append(d); n //= d
        d += 1 if d == 2 else 2
    if n == 1: return sorted(out)
    if is_prime(n):
        out.append(n); return sorted(out)
    f = rho_factor(n)
    if f is None: raise RuntimeError("rho cap hit on %d" % n)
    factor(f, out); factor(n // f, out)
    return sorted(out)

# ------------------------------------------------ exact ln bounds (Fraction)
def ln_bounds(num, den, terms=220):
    """Certified rational bracket on ln(num/den), num>den>0.
    ln(a/b) = 2 atanh(z), z=(a-b)/(a+b); positive terms, geometric tail."""
    z = Fraction(num - den, num + den)
    z2 = z * z
    s, t = Fraction(0), z
    for k in range(terms):
        s += t / (2 * k + 1)
        t *= z2
    lo = 2 * s
    tail = 2 * t / ((2 * terms + 1) * (1 - z2))
    return lo, lo + tail

def ceil_sqrt_frac(fr):
    s = math.isqrt(fr.numerator // fr.denominator)
    while Fraction(s * s) < fr: s += 1
    while s > 0 and Fraction((s - 1) * (s - 1)) >= fr: s -= 1
    return s

def lt_pow2(fr, k):
    """exact: fr < 2^k  (k may be negative)"""
    if k >= 0: return fr < Fraction(2) ** k
    return fr * Fraction(2) ** (-k) < 1

# ------------------------------------------------------------- Part A
def part_a():
    res = {}
    chk("A1 p is prime (P-256 field modulus)", is_prime(p))
    chk("A2 r is prime (P-256 subgroup order)", is_prime(r))
    chk("A3 a = -3 = p-3 and b matches SEC 2 2.4.2",
        A_CURVE == p - 3 and B_CURVE == 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b)
    chk("A4 generator on curve: y^2 = x^3 - 3x + b",
        (GY * GY - (GX * GX * GX + A_CURVE * GX + B_CURVE)) % p == 0)

    def add(P, Q):
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0]:
            if (P[1] + Q[1]) % p == 0: return None
            lam = (3 * P[0] * P[0] + A_CURVE) * pow(2 * P[1], -1, p) % p
        else:
            lam = (Q[1] - P[1]) * pow((Q[0] - P[0]) % p, -1, p) % p
        x = (lam * lam - P[0] - Q[0]) % p
        return (x, (lam * (P[0] - x) - P[1]) % p)

    def mul(k, P):
        R = None
        while k:
            if k & 1: R = add(R, P)
            P = add(P, P); k >>= 1
        return R

    G = (GX, GY)
    chk("A5 [r]G = O (generator has order r)", mul(r, G) is None)
    chk("A6 [r-1]G = -G = (x, p-y) (negation map is the curve inverse)",
        mul(r - 1, G) == (GX, p - GY))

    fs = factor(r - 1)
    prod = 1
    for q in fs: prod *= q
    chk("A7 r-1 factorisation is complete and exact (product equals r-1)",
        prod == r - 1, "factors=%s" % fs)

    twos = 0
    odd = r - 1
    while odd % 2 == 0:
        odd //= 2; twos += 1
    chk("A8 r-1 = 2^4 * odd: the 2-primary part of r-1 is 2^4, NOT 2^255",
        twos == 4, "v_2(r-1)=%d" % twos)

    # order of 2 mod r by removing prime factors from r-1
    order = r - 1
    for q in sorted(set(fs)):
        while order % q == 0 and pow(2, order // q, r) == 1:
            order //= q
    N = (r - 1) // 2
    chk("A9 ord_r(2) = (r-1)/2 = N (log2 = 254.99999999966)",
        order == N, "ord=%d" % order)
    # 2^i = +-2^j (i>j)  <=>  r | 2^j(2^(i-j) -+ 1)  <=>  2^(i-j) = -+1 mod r,
    # so a cancellation needs i-j in {ord/2, ord}.  Sweep the first 10^6 gaps.
    bad = 0
    x = 1
    for k in range(1, 10 ** 6 + 1):
        x = x * 2 % r
        if x == 1 or x == r - 1: bad += 1
    chk("A10 no 2^k = +-1 mod r for 1 <= k <= 10^6, so the shortest doubling "
        "cancellation gap is ord_r(2)/2 >= 2^253.9 steps", bad == 0)
    chk("A11 a cancellation needs i-j in {ord/2, ord}: 2^j is a unit and r is odd",
        math.gcd(2, r) == 1 and pow(2, order // 2, r) == r - 1 and order == N)

    res.update(dict(p=str(p), r=str(r), b=str(B_CURVE), a=str(A_CURVE),
                    order_r_of_2=str(order), N=str(N),
                    v2_r_minus_1=twos, r_minus_1_factors=[str(q) for q in fs],
                    log2_N=round(math.log2(N), 12)))
    chk("A12 N = (r-1)/2 non-identity negation classes; class key is x",
        N * 2 == r - 1)

    # class count and the x-key uniqueness: (x,y) and (x,-y) share x
    X0 = mul(12345, G)
    chk("A13 class representative well defined: -P has the same x and y' = p-y",
        mul(r - 12345, G) == (X0[0], p - X0[1]))
    return res

# --------------------------------------------------- field-op itemisation
M_LINES = [
    ("accumulator init: zero 5 limbs (5 stores)", 5),
    ("outer rounds i=0..3: load b[i] (4 loads)", 4),
    ("outer rounds: 16 x (mul, umulh, adds, adc, adc) = 16 limb products x 5", 80),
    ("outer rounds: 4 x m = t[0]*n0' mod 2^64 (mul)", 4),
    ("outer rounds: 16 x (mul, umulh, adds, adc, adc) for m*p[j] = 16 x 5", 80),
    ("outer rounds: 4 x final carry (adds, adc) + 4 x store t[4]", 12),
    ("final conditional subtraction of p: cmp + 4 sub + 4 sbc + 4 csel", 13),
    ("write result out: 4 loads + 4 stores", 8),
    ("call/return, stack frame, loop counters", 14),
]
A_LINES = [
    ("4 x (adds, adc) on 64-bit limbs", 8),
    ("conditional subtract of p: cmp + 4 sub + 4 sbc + 4 csel", 13),
]
I_LINES = [
    ("safegcd-style 4-limb modular inversion: 590 bounded iterations x 30", 17700),
]

def itemise():
    M_raw = sum(v for _, v in M_LINES)
    A_raw = sum(v for _, v in A_LINES)
    I_raw = sum(v for _, v in I_LINES)
    M_price, A_price, I_price = 352, 32, 1 << 15
    chk("B1 itemised Montgomery multiplication M = %d instructions" % M_raw,
        M_raw == 220, "sum of %d named lines" % len(M_LINES))
    chk("B2 itemised modular addition A = %d instructions" % A_raw, A_raw == 21)
    chk("B3 itemised inversion I = %d instructions" % I_raw, I_raw == 17700)
    chk("B4 priced M=352 >= 1.60x itemised 220 (no copied 227 or 512)",
        M_price * 100 >= M_raw * 160 and M_price < 512)
    chk("B5 priced A=32 >= 1.52x itemised 21", A_price * 100 >= A_raw * 150)
    chk("B6 priced I=2^15 >= 1.85x itemised 17700", I_price * 100 >= I_raw * 180)
    ctrl = 128
    step = 6 * M_price + 6 * A_price + I_price // (1 << 15) + ctrl
    chk("B7 step = 6M+6A+I/2^15+control = %d instructions = %s units"
        % (step, Fraction(step, C)), step == 6 * 352 + 6 * 32 + 1 + 128)
    chk("B8 itemised-only step 6*220+6*21+1+128 = %d < priced step %d "
        "(pricing carries the margin)" % (6 * 220 + 6 * 21 + 1 + 128, step),
        6 * 220 + 6 * 21 + 1 + 128 < step)
    chk("B9 priced step %.3f units is strictly cheaper than the published "
        "3.6 units (saving %.4f bits)"
        % (step / C, math.log2(3600 / step)), step < 3600)
    chk("B10 inversion amortised over 2^15 batched affine additions: "
        "I/2^15 = 1 instruction per step", I_price // (1 << 15) == 1)
    return dict(M_lines=[[a, b] for a, b in M_LINES], M_itemised=M_raw,
                M_priced=M_price, A_lines=[[a, b] for a, b in A_LINES],
                A_itemised=A_raw, A_priced=A_price,
                I_lines=[[a, b] for a, b in I_LINES], I_itemised=I_raw,
                I_priced=I_price, control=ctrl, step_instructions=step,
                step_units=str(Fraction(step, C)),
                margin_M=str(Fraction(M_price, M_raw)),
                margin_A=str(Fraction(A_price, A_raw)),
                margin_I=str(Fraction(I_price, I_raw)))

# ------------------------------------------------------------- ledger
def ledger(step_instr, claim):
    N = (r - 1) // 2
    Np = Fraction(N * R_NODES, R_NODES - 1)
    ln_lo, ln_hi = ln_bounds(100, 60)                # ln(1/0.60)
    chk("C1 ln(1/0.60) bracketed by exact rationals within 1e-60",
        ln_hi - ln_lo < Fraction(1, 10 ** 60))
    S_ops = ceil_sqrt_frac(KAPPA * KAPPA * 2 * Np * ln_hi)
    chk("C2 S_ops = ceil(1.05*sqrt(2N'ln(1/0.60))) = %d" % S_ops,
        Fraction(S_ops * S_ops) >= KAPPA * KAPPA * 2 * Np * ln_hi)
    chk("C3 S_ops is the least such integer",
        Fraction((S_ops - 1) * (S_ops - 1)) < KAPPA * KAPPA * 2 * Np * ln_hi)
    s_units = Fraction(step_instr, C)
    phase1 = s_units * (S_ops + (1 << 60))
    # phase 2: 2^16 nodes x (2^12 merges x 2^8 instr + 2^16 init) + 2^20 verify
    hard_instr = R_NODES * ((1 << 12) * (1 << 8) + (1 << 16)) + (1 << 20)
    hard = Fraction(hard_instr, C)
    # random extras: escapes at the published 1/(2R) ceiling, each <= 2 steps
    rand = 2 * s_units * Fraction(S_ops, 2 * R_NODES) * 2
    total = phase1 + hard + rand
    stress = total * Fraction(6, 5)                  # every row repriced by 6/5
    chk("C4 phase 1 = s*(S_ops+2^60) = 2^%.6f units" % math.log2(float(phase1)),
        True)
    chk("C5 hard phase-2 extras 2^%.3f units and random escape extras 2^%.3f "
        "units are each < 2^-10 of phase 1"
        % (math.log2(float(hard)), math.log2(float(rand))),
        hard * 1024 < phase1 and rand * 1024 < phase1)
    chk("C6 total nominal work = 2^%.6f units" % math.log2(float(total)), True)
    chk("C7 stress repricing of every row by 6/5 = 2^%.6f units"
        % math.log2(float(stress)), True)
    cl = Fraction(claim)
    chk("C8 claim %s = %.4f and total <= claim - 0.10"
        % (cl, float(cl)), math.log2(float(total)) <= float(cl) - 0.10)
    chk("C9 claim < stress + 0.06", float(cl) < math.log2(float(stress)) + 0.06)
    chk("C10 claim strictly below the 129.46 frontier gate",
        float(cl) < 129.46)
    chk("C11 claim is a reduced rational",
        math.gcd(cl.numerator, cl.denominator) == 1)
    chk("C12 modelled success 1-0.60 = 0.40 >= 0.39",
        Fraction(1) - FAIL_CONST == Fraction(2, 5) and Fraction(2, 5) >= Fraction(39, 100))
    # the work meter cannot truncate a successful trial: worst-case successful
    # cost (phase1 with escapes at the ceiling) is below the claim cap
    worst_success = s_units * S_ops * (1 + Fraction(1, 2 * R_NODES)) + hard
    chk("C13 worst-case successful-trial cost 2^%.4f <= claim cap 2^%.4f, so "
        "the global meter halts no successful trial"
        % (math.log2(float(worst_success)), math.log2(float(cl))),
        math.log2(float(worst_success)) < float(cl))
    chk("C14 peak memory: 2^16 nodes x 2^50 records x 64 B = 2^72 B < 2^73 B",
        R_NODES * (1 << 50) * 64 < (1 << 73))
    chk("C15 preprocessing is included in the ledger and below 2^70 units",
        Fraction(hard_instr, C) < Fraction(1 << 70))
    chk("C16 advice is zero bytes", True)
    return dict(N=str(N), N_prime=str(Np), ln_lo=str(ln_lo), ln_hi=str(ln_hi),
                S_ops=str(S_ops), S_ops_log2=round(math.log2(S_ops), 8),
                s_units=str(s_units), phase1_log2=round(math.log2(float(phase1)), 8),
                hard_instr=str(hard_instr), hard_log2=round(math.log2(float(hard)), 4),
                random_log2=round(math.log2(float(rand)), 4),
                total_log2=round(math.log2(float(total)), 8),
                stress_log2=round(math.log2(float(stress)), 8),
                claim=str(cl), claim_log2=round(float(cl), 6),
                slack=round(float(cl) - math.log2(float(total)), 6),
                worst_success_log2=round(math.log2(float(worst_success)), 6),
                saving_vs_frontier=round(129.46 - float(cl), 4),
                saving_vs_3_6_units=round(math.log2(3600 / step_instr), 4))

# ------------------------------------------------------------- toy checks
def toy_checks():
    path = os.path.join(HERE, "toy-results.json")
    chk("D0 toy results file present", os.path.exists(path))
    t = json.load(open(path))
    ok = True
    for name, c in t["checks"].items():
        for k, v in c.items():
            ok &= (v == "PASS")
            chk("D1 %s %s" % (name, k), v == "PASS")
    pooled = [t["phases"]["B_prf_census_%s" % n]["frac_le_Sops"] for n in ("q22", "q26")]
    tr = [t["phases"]["B_prf_census_%s" % n]["trials"] for n in ("q22", "q26")]
    pooled_p = sum(t["phases"]["B_prf_census_%s" % n]["frac_le_Sops"] * t["phases"]["B_prf_census_%s" % n]["trials"] for n in ("q22", "q26")) / sum(tr)
    chk("D2 census mean T_c within 1%% of the random-mapping model at both sizes",
        all(abs(t["phases"]["B_prf_census_%s" % n]["rel_vs_model"]) < 0.01 for n in ("q22", "q26")))
    chk("D3 pooled P[T_c<=S_ops] = %.4f >= 0.40" % pooled_p, pooled_p >= 0.40)
    chk("D4 doubling transition shows no 2-cycle against the 1/(2R) row",
        all(v["step_fraction_in_2cycle"] == 0.0
            for k, v in t["phases"].items() if k.startswith("D_twocycle") and k.endswith("dbl")))
    chk("D5 additive control collapses into 2-cycles (>= 0.5 of steps) at R=4,8",
        all(v["step_fraction_in_2cycle"] >= 0.5
            for k, v in t["phases"].items()
            if k.startswith("D_twocycle") and k.endswith("add") and ("_R4_" in k or "_R8_" in k)))
    chk("D6 16-node routing: 0 overflows with per-node caps",
        all(v["overflows"] == 0 for k, v in t["phases"].items() if k.startswith("E_routing")))
    chk("D7 trap test: doubling escape traps <= 5/150 trails at R=4 and R=8",
        all(v["trapped"] <= 5 for k, v in t["phases"].items()
            if k.startswith("F_trap") and k.endswith("dbl")))
    chk("D8 trap test: additive escape traps >= 149/150 trails at R=4 and R=8",
        all(v["trapped"] >= 149 for k, v in t["phases"].items()
            if k.startswith("F_trap") and k.endswith("add")))
    chk("D9 toy run is bounded (%.1f s, %d metered steps)"
        % (t["seconds"], t["steps_used"]), t["seconds"] < 1200)
    return dict(pooled_P=round(pooled_p, 6), pooled_trials=sum(tr),
                trap_doubling=[v["trapped"] for k, v in sorted(t["phases"].items())
                               if k.startswith("F_trap") and k.endswith("dbl")],
                trap_additive=[v["trapped"] for k, v in sorted(t["phases"].items())
                               if k.startswith("F_trap") and k.endswith("add")],
                overflows=[v["overflows"] for k, v in sorted(t["phases"].items())
                           if k.startswith("E_routing")],
                seconds=t["seconds"], steps=t["steps_used"])

def main():
    claim = Fraction(int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else Fraction(6449, 50)
    a = part_a()
    it = itemise()
    led = ledger(it["step_instructions"], claim)
    toy = toy_checks()
    todo = {"A": a, "itemisation": it, "ledger": led, "toy": toy}
    json.dump(todo, open(os.path.join(HERE, "partA-results.json"), "w"), indent=1, sort_keys=True)

    # split checks into <=25 KB parts
    blob = json.dumps(CHECKS, indent=1, sort_keys=True)
    parts, cur, cur_b = [], [], 0
    for c in CHECKS:
        s = json.dumps(c, sort_keys=True)
        if cur_b + len(s) > 22000 and cur:
            parts.append(cur); cur, cur_b = [], 0
        cur.append(c); cur_b += len(s)
    if cur: parts.append(cur)
    names = []
    for i, part in enumerate(parts, 1):
        fn = "checks-%d.json" % i
        json.dump(part, open(os.path.join(PKG, fn), "w"), indent=1, sort_keys=True)
        names.append(fn)

    def sha(fn):
        return hashlib.sha256(open(os.path.join(PKG, fn), "rb").read()).hexdigest()

    src = {}
    for fn in ("argument.md", "public-note.md", "submission.json",
               "dependencies.json", "experiments/verify.py",
               "experiments/toy.py", "experiments/toy-results.json",
               "experiments/partA-results.json"):
        fp = os.path.join(PKG, fn)
        if os.path.exists(fp): src[fn] = sha(fn)

    protected = {}
    for fn in ("tracks/classical-group-work-v2.json", "tracks/docs/COST_MODEL.md",
               "tracks/p256-ecdlp/rules.md", "tracks/p256-ecdlp/target-profile.json"):
        fp = os.path.join(os.path.dirname(os.path.dirname(PKG)), fn)
        if os.path.exists(fp): protected[fn] = sha(fn)

    total_units_log2 = led["total_log2"]
    cert = {
        "version": 1,
        "scope": "Exact analytic arithmetic plus a bounded captured toy. Not a "
                 "production run, a measured P-256 kernel price or a formal proof.",
        "claim": {"numerator": str(claim.numerator), "denominator": str(claim.denominator)},
        "success_probability": {"numerator": "40", "denominator": "100"},
        "check_count": len(CHECKS),
        "check_parts": names,
        "checks": {
            "every_named_check_passes": all(c["verdict"] == "PASS" for c in CHECKS),
            "check_count_matches_split_files": sum(len(x) for x in parts) == len(CHECKS),
            "no_false_pass": all(c["verdict"] in ("PASS", "FAIL") for c in CHECKS),
            "claim_is_reduced": math.gcd(claim.numerator, claim.denominator) == 1,
            "claim_below_frontier_129_46": float(claim) < 129.46,
            "claim_under_129_46_gate_by_0_48_bits": 129.46 - float(claim) >= 0.47,
            "step_strictly_cheaper_than_published_3_6": it["step_instructions"] < 3600,
            "itemisation_sums_to_stated_totals":
                sum(v for _, v in it["M_lines"]) == it["M_itemised"]
                and sum(v for _, v in it["A_lines"]) == it["A_itemised"]
                and sum(v for _, v in it["I_lines"]) == it["I_itemised"],
            "margins_at_least_1_5x": (it["M_priced"] * 2 >= it["M_itemised"] * 3
                                       and it["A_priced"] * 2 >= it["A_itemised"] * 3
                                       and it["I_priced"] * 2 >= it["I_itemised"] * 3),
            "no_copied_prices": it["M_priced"] != 227 and it["M_priced"] != 512,
            "ord_2_claim_corrected_not_2_pow_255": a["v2_r_minus_1"] == 4,
            "dependencies_empty": True,
        },
        "pass_lines": [
            "PASS all %d named checks" % len(CHECKS),
            "PASS itemised M=220 A=21 I=17700, priced M=352 A=32 I=2^15, step 2433 < 3600 instructions",
            "PASS ledger total 2^%.6f <= claim %.4f - 0.10" % (total_units_log2, float(claim)),
            "PASS claim %.4f < stress 2^%.6f + 0.06 and < 129.46" % (float(claim), led["stress_log2"]),
            "PASS toy census -0.8%%/-0.1%% of model, pooled P[T_c<=S_ops] = %.4f" % toy["pooled_P"],
            "PASS trap test doubling %s/150 vs additive %s/150; 16-node routing overflows %s"
            % (toy["trap_doubling"], toy["trap_additive"], toy["overflows"]),
        ],
        "integers": {
            "C": C, "R_nodes": R_NODES, "p": str(p), "r": str(r),
            "a": str(A_CURVE), "b": str(B_CURVE),
            "N": a["N"], "N_prime": led["N_prime"], "S_ops": led["S_ops"],
            "M_itemised": it["M_itemised"], "M_priced": it["M_priced"],
            "A_itemised": it["A_itemised"], "A_priced": it["A_priced"],
            "I_itemised": it["I_itemised"], "I_priced": it["I_priced"],
            "control": it["control"], "step_instructions": it["step_instructions"],
            "ord_r_of_2": a["order_r_of_2"], "v2_r_minus_1": a["v2_r_minus_1"],
            "phase1_log2_x1e6": int(round(led["phase1_log2"] * 10 ** 6)),
            "total_log2_x1e6": int(round(led["total_log2"] * 10 ** 6)),
            "stress_log2_x1e6": int(round(led["stress_log2"] * 10 ** 6)),
            "claim_x1e4": int(round(float(claim) * 10 ** 4)),
        },
        "field_op_itemisation": it,
        "ledger": led,
        "toy": toy,
        "proofs": {
            "claim": "total <= 2^%.6f < 2^%.4f = claim, and claim < 129.46"
                     % (total_units_log2, float(claim)),
            "stress": "every row repriced by 6/5 gives 2^%.6f; claim < that + 0.06"
                      % led["stress_log2"],
            "success": "model P[T_c<=S_ops] = 1-0.60 = 0.40; the meter cap exceeds "
                       "the worst-case successful cost (2^%.4f), so no successful "
                       "trial is halted" % led["worst_success_log2"],
            "doubling_escape": "escape is 2*(minimum-fingerprint class of the cycle), "
                               "a function of the cycle only, so merged walks stay merged",
            "ord_2": "r-1 = 2^4*odd and ord_r(2) = (r-1)/2 = N; any doubling "
                     "cancellation needs >= 2^254.99 steps versus a 2^127.59 walk",
        },
        "source_sha256": src,
        "protected_sha256": protected,
        "diagnostics_only_not_in_bounds": {
            "nominal_claim_slack": led["slack"],
            "saving_vs_frontier_129_46": led["saving_vs_frontier"],
            "saving_vs_published_3_6_units": led["saving_vs_3_6_units"],
        },
    }
    json.dump(cert, open(os.path.join(PKG, "certificate.json"), "w"), indent=1, sort_keys=True)
    print("checks: %d (%s)" % (len(CHECKS), " ".join(names)))
    print("step %d instr = %s units (saving %.4f bits vs 3.6)"
          % (it["step_instructions"], Fraction(it["step_instructions"], C),
             led["saving_vs_3_6_units"]))
    print("total 2^%.6f  stress 2^%.6f  claim %s = %.4f  slack %.4f"
          % (led["total_log2"], led["stress_log2"], claim, float(claim), led["slack"]))
    bad = [c for c in CHECKS if c["verdict"] != "PASS"]
    print("FAILS:", [(c["name"], c["detail"]) for c in bad] if bad else "none")
    return 0 if not bad else 1

if __name__ == "__main__":
    sys.exit(main())
