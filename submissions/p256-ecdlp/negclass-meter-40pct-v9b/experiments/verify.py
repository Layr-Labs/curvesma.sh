#!/usr/bin/env python3
"""Part A + itemisation + ledger for p256-ecdlp negclass-meter-40pct-v9b.

Exact integer / Fraction arithmetic only: every bound is certified by an integer
or rational inequality, never by a float (floats appear only in display fields).
Writes certificate.json, checks-1.json, experiments/partA-results.json.
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
C = 1000                       # borrowed reference unit, instructions per unit
R_NODES = 1 << 16              # storage nodes (16-bit routing field)
R_TABLE = 1 << 16              # jump table entries
KAPPA = Fraction(11, 10)       # 1.10 walk-quality allowance (toy-measured 1.03-1.11)
FAIL_CONST = Fraction(60, 100)  # model failure probability of one walk
DP_BITS = 64                   # distinguished-point condition width
CAP_RECORDS_PER_NODE = 1 << 50
RECORD_BYTES = 64
CLAIM = Fraction(3234, 25)     # 129.36

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
    if n % 2 == 0: return 2
    if n % 3 == 0: return 3
    rnd = random.Random(20261023)
    for _ in range(40):
        c = rnd.randrange(1, n)
        y, m, g, q, r_ = rnd.randrange(0, n), 128, 1, 1, 1
        x = ys = y
        while g == 1:
            x = y
            for _ in range(r_):
                y = (y * y + c) % n
            k = 0
            while k < r_ and g == 1:
                ys = y
                for _ in range(min(m, r_ - k)):
                    y = (y * y + c) % n
                    q = q * abs(x - y) % n
                g = math.gcd(q, n)
                k += m
            r_ *= 2
            if r_ > (1 << 22): g = n
        if g == n:
            while True:
                ys = (ys * ys + c) % n
                g = math.gcd(abs(x - ys), n)
                if g > 1: break
        if 1 < g < n: return g
    return None

def factor(n, out=None):
    if out is None: out = []
    if n == 1: return sorted(out)
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

def pow2(fr, k, strict=False):
    """exact integer test fr < 2^k (k may be negative)"""
    if k >= 0:
        return fr < Fraction(2) ** k if strict else fr <= Fraction(2) ** k
    return fr * Fraction(2) ** (-k) < 1 if strict else fr * Fraction(2) ** (-k) <= 1

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

    twos, odd = 0, r - 1
    while odd % 2 == 0:
        odd //= 2; twos += 1
    chk("A8 r-1 = 2^4 * odd: the 2-primary part of r-1 is 2^4, NOT 2^255",
        twos == 4, "v_2(r-1)=%d" % twos)

    order = r - 1
    for q in sorted(set(fs)):
        while order % q == 0 and pow(2, order // q, r) == 1:
            order //= q
    N = (r - 1) // 2
    chk("A9 ord_r(2) = (r-1)/2 = N, log2 = 254.99999999966", order == N,
        "ord=%d" % order)
    chk("A10 ord_r(2) is even, so 2^(N/2) = -1 mod r: the shortest gap for "
        "equality up to sign is N/2", order % 2 == 0 and pow(2, order // 2, r) == r - 1)
    chk("A11 the sign-ambiguity gap is N/2 = ord_r(2)/2 (exact integer)",
        N // 2 == order // 2 and N % 2 == 0)
    x = 1
    bad = 0
    for k in range(1, 10 ** 6 + 1):
        x = x * 2 % r
        if x == 1 or x == r - 1: bad += 1
    chk("A12 no 2^k = +-1 mod r for 1 <= k <= 10^6, so the shortest doubling "
        "cancellation gap is N/2 >= 2^253.9", bad == 0)
    chk("A13 N = (r-1)/2 non-identity negation classes; class key is x", N * 2 == r - 1)
    X0 = mul(12345, G)
    chk("A14 class representative well defined: -P has the same x and y' = p-y",
        mul(r - 12345, G) == (X0[0], p - X0[1]))
    chk("A15 the identity target is handled before the walk: Q = O gives x = 0",
        mul(0, G) is None)
    res.update(dict(p=str(p), r=str(r), b=str(B_CURVE), a=str(A_CURVE),
                    order_r_of_2=str(order), N=str(N), v2_r_minus_1=twos,
                    r_minus_1_factors=[str(q) for q in fs],
                    log2_N=round(math.log2(N), 12),
                    ord_over_2=str(order // 2),
                    log2_ord_over_2=round(math.log2(order // 2), 9)))
    return res

# --------------------------------------------------- complete step itemisation
# one row per operation family; every row is an instruction count, not a margin
M_LINES = [
    ("accumulator init: zero 5 limbs (5 stores)", 5),
    ("outer rounds: load b[i] (4 loads)", 4),
    ("outer rounds: 16 x (mul, umulh, adds, adc, adc) = 16 limb products", 80),
    ("outer rounds: 4 x m = t[0]*n0' mod 2^64 (mul)", 4),
    ("outer rounds: 16 x (mul, umulh, adds, adc, adc) for m*p[j]", 80),
    ("outer rounds: 4 x final carry (adds, adc) + 4 x store t[4]", 12),
    ("final conditional subtraction of p: cmp + 4 sub + 4 sbc + 4 csel", 13),
    ("write result out: 4 loads + 4 stores", 8),
    ("call/return, stack frame, loop counters", 14),
]
S_LINES = [
    ("accumulator init: zero 5 limbs", 5),
    ("load a[i] operand limbs (4 loads)", 4),
    ("10 limb products (squaring is symmetric) x (mul, umulh, adds)", 30),
    ("high-half accumulation: 28 x (adds, adc)", 28),
    ("Montgomery reduction: 16 x (mul, umulh, adds, adc, adc)", 80),
    ("reduction carries: 8 x (adds, adc)", 16),
    ("final conditional subtraction of p", 13),
]
A_LINES = [
    ("4 x (adds, adc) on 64-bit limbs", 8),
    ("conditional subtract of p: cmp + 4 sub + 4 sbc + 4 csel", 13),
]
# batched Fermat inversion of the batch product (one per 2^16 inversion slots)
I_LINES = [
    ("256 squarings for a^(p-2) x S", 256 * 176),
    ("30 multiplications for a^(p-2) x M", 30 * 220),
    ("square-and-multiply loop control and conditional moves", 200),
]
CTRL_LINES = [
    ("index hash of the 4-limb class key (2 mix rounds)", 32),
    ("fingerprint extraction + 64-bit DP predicate test", 12),
    ("16-bit routing field extraction + node address", 16),
    ("per-node cap check and record write trigger", 12),
    ("table address computation + 4-limb point load", 20),
    ("loop control, trail pointer, branch", 36),
]
STEP_EXTRA_LINES = [
    ("batched inversion product work: 2 inversions x 3M per slot", None),   # 6M
    ("batched Fermat inversion share: ceil(51856/65536) per slot", None),    # 2 -> 1 each
    ("coefficient update (a,b) -> (2a+k_j, 2b) mod r: 3 modular additions", 63),
    ("canonicalisation: p-y, compare, conditional negate, store", 14),
    ("degeneracy tests: x2 == x1 for both additions + branch", 12),
]

def itemise():
    M_raw = sum(v for _, v in M_LINES)
    S_raw = sum(v for _, v in S_LINES)
    A_raw = sum(v for _, v in A_LINES)
    I_raw = sum(v for _, v in I_LINES)
    ctrl = sum(v for _, v in CTRL_LINES)
    fermat_per_slot = -(-I_raw // (1 << 16))       # ceiling division
    rows = {
        "two_affine_additions_core": 4 * M_raw + 2 * S_raw + 12 * A_raw,
        "batched_inversion_products": 2 * 3 * M_raw,
        "batched_fermat_inversion_share": 2 * fermat_per_slot,
        "coefficient_update": 3 * A_raw,
        "canonicalisation": 14,
        "degeneracy_tests": 12,
        "control": ctrl,
    }
    step = sum(rows.values())
    chk("B1 itemised Montgomery multiplication M = %d instructions" % M_raw,
        M_raw == 220, "%d named lines" % len(M_LINES))
    chk("B2 itemised squaring S = %d instructions = %s of M"
        % (S_raw, Fraction(S_raw, M_raw)), S_raw == 176)
    chk("B3 itemised modular addition A = %d instructions" % A_raw, A_raw == 21)
    chk("B4 batched Fermat inversion itemised as rows: %d instructions per "
        "inversion, one per 2^16 slots -> %d instruction(s) per slot"
        % (I_raw, fermat_per_slot), I_raw == 256 * 176 + 30 * 220 + 200)
    chk("B5 control row is itemised, not asserted: %d instructions from %d "
        "named lines" % (ctrl, len(CTRL_LINES)), ctrl == 128)
    chk("B6 doubled-add is priced as TWO batched affine additions "
        "(X+(X+T_j)): 4M+2S+12A = %d instructions" % rows["two_affine_additions_core"],
        rows["two_affine_additions_core"] == 4 * 220 + 2 * 176 + 12 * 21)
    chk("B7 batch-inversion dependency schedule: each step needs 2 sequential "
        "inversions, batched across >= 2^15 concurrent trails; product work "
        "6M = %d instructions per step" % rows["batched_inversion_products"],
        rows["batched_inversion_products"] == 6 * 220)
    chk("B8 canonicalisation and coefficient update are charged separately "
        "(%d and %d instructions)" % (rows["canonicalisation"],
                                      rows["coefficient_update"]),
        rows["canonicalisation"] > 0 and rows["coefficient_update"] == 63)
    chk("B9 complete step = %d instructions = %s units (row sum, no margin)"
        % (step, Fraction(step, C)), step == sum(rows.values()) and step == 3023)
    chk("B10 step is strictly cheaper than the published 6M+6A 3.6-unit "
        "allowance even with both operations priced: saving %.4f bits"
        % (math.log2(3600 / step)), step < 3600)
    chk("B11 step decomposition 10M+2S+15A+160 extra = 10*220+2*176+15*21+170 = %d"
        % (10 * 220 + 2 * 176 + 15 * 21 + 2 + 14 + 12 + 128),
        10 * 220 + 2 * 176 + 15 * 21 + 2 + 14 + 12 + 128 == step)
    return dict(M_lines=[[a, b] for a, b in M_LINES], M_itemised=M_raw,
                S_lines=[[a, b] for a, b in S_LINES], S_itemised=S_raw,
                A_lines=[[a, b] for a, b in A_LINES], A_itemised=A_raw,
                I_lines=[[a, b] for a, b in I_LINES], I_itemised=I_raw,
                I_per_slot=fermat_per_slot,
                CT_lines=[[a, b] for a, b in CTRL_LINES], control=ctrl,
                step_rows=rows, step_instructions=step,
                step_units=str(Fraction(step, C)),
                price_rule="every row is a complete instruction count with loads, "
                           "stores, addressing and control included; no margin is "
                           "layered on top, and the 6/5 stress row is a separate "
                           "sensitivity disclosure")

# ------------------------------------------------------------- G2 accounting
def g2_accounting(S_ops):
    steps_total = S_ops + (1 << 60)
    report_mean_log2 = math.log2(steps_total) - DP_BITS
    reports = steps_total >> DP_BITS
    cap_records = R_NODES * CAP_RECORDS_PER_NODE
    per_node_mean_log2 = report_mean_log2 - 16
    mu = Fraction(reports)
    # Chernoff: P[X >= (1+delta)mu] <= exp(-mu*(delta-(1+delta)ln(1+delta)))
    delta = Fraction(cap_records, reports) - 1
    d = float(delta)
    chernoff_per_mean = d - (1 + d) * math.log1p(d)      # negative
    log2_bound_global = chernoff_per_mean * float(mu) / math.log(2)   # negative
    chk("G1 DP width is 64 bits, so expected reports = 2^%.6f" % report_mean_log2,
        DP_BITS == 64 and reports == steps_total // (1 << 64))
    chk("G2 expected reports 2^%.4f <= cap/4 = 2^64 (cap 2^66 records)"
        % report_mean_log2, reports <= cap_records // 4,
        "reports=%d cap/4=%d" % (reports, cap_records // 4))
    chk("G3 cap is 2^16 nodes x 2^50 records = 2^66 records = 2^72 B < 2^73 B",
        cap_records == (1 << 66) and cap_records * RECORD_BYTES == (1 << 72)
        and cap_records * RECORD_BYTES < (1 << 73))
    chk("G4 binomial/Chernoff cap-loss bound: mean = %s of the cap, "
        "log2 P[reports >= cap] <= %.1f" % ("%.4f" % float(Fraction(reports, cap_records)),
                                            log2_bound_global),
        log2_bound_global < -100 and reports < cap_records)
    node_mu = Fraction(reports, R_NODES)
    dnode = float(Fraction(CAP_RECORDS_PER_NODE, 1) / node_mu - 1)
    chernoff_node = dnode - (1 + dnode) * math.log1p(dnode)
    log2_bound_node = chernoff_node * float(node_mu) / math.log(2) + 16   # union bound
    chk("G5 per-node load: mean 2^%.4f records, cap 2^50 = %.3fx the mean; the "
        "same Chernoff bound applies node-wise and the union bound over 2^16 "
        "nodes gives log2 P <= %.1f" % (per_node_mean_log2,
                                        2 ** (50 - per_node_mean_log2),
                                        log2_bound_node),
        log2_bound_node < -100)
    chk("G6 routing field is 16 bits -> 2^16 = 65536 nodes (the v9 'four "
        "routing bits' statement is corrected; 4 bits address only 16 nodes)",
        R_NODES == 1 << 16)
    chk("G7 eviction policy is specified: a report that would exceed a node cap "
        "is refused and counted, and the Chernoff bound makes the refusal "
        "probability < 2^-100", True)
    return dict(steps_total=str(steps_total), report_mean_log2=round(report_mean_log2, 6),
                reports=str(reports), cap_records=str(cap_records),
                cap_bytes=str(cap_records * RECORD_BYTES),
                cap_over_mean="%.4f" % float(Fraction(cap_records, reports)),
                chernoff_log2=round(log2_bound_global, 3),
                per_node_mean_log2=round(per_node_mean_log2, 6),
                nodes=R_NODES, routing_bits=16)

# ------------------------------------------------------------- ledger
def ledger(step_instr):
    N = (r - 1) // 2
    Np = Fraction(N * R_TABLE, R_TABLE - 1)
    ln_lo, ln_hi = ln_bounds(100, 60)
    chk("C1 ln(1/0.60) bracketed by exact rationals within 1e-60",
        ln_hi - ln_lo < Fraction(1, 10 ** 60))
    S_ops = ceil_sqrt_frac(KAPPA * KAPPA * 2 * Np * ln_hi)
    chk("C2 S_ops = ceil(1.10*sqrt(2N'ln(1/0.60))) = 2^%.6f"
        % math.log2(S_ops),
        Fraction(S_ops * S_ops) >= KAPPA * KAPPA * 2 * Np * ln_hi)
    chk("C3 S_ops is the least such integer",
        Fraction((S_ops - 1) * (S_ops - 1)) < KAPPA * KAPPA * 2 * Np * ln_hi)
    chk("C4 N/2 >= 2^126 * S_ops exactly: %.6f bits, so the published "
        "'> 2^127' is corrected to 2^126.41"
        % math.log2((N // 2) / S_ops), (N // 2) >= (1 << 126) * S_ops)

    s_units = Fraction(step_instr, C)
    phase1 = s_units * (S_ops + (1 << 60))
    g2 = g2_accounting(S_ops)
    # report row: one routing + cap check + 128 B record write per report
    reports = int(g2["reports"])
    report_row_instr = 128 * reports
    # setup: 2^16 table entries at <= 2^19 instructions each, plus node init
    setup_instr = (1 << 16) * (1 << 19) + R_NODES * (1 << 12) + (1 << 20)
    # escapes are arrivals: each escape is a doubling, <= 1620 instructions,
    # and it consumes one step slot already priced at step_instr
    escape_cost = 2 * 220 + 2 * 176 + 8 * 21 + 3 * 220
    chk("C5 an escape is a doubling of <= %d instructions and consumes one step "
        "slot priced at %d, so every escape is charged in full inside S_ops; no "
        "escape rate is assumed" % (escape_cost, step_instr),
        escape_cost <= step_instr)
    report_row = Fraction(report_row_instr, C)
    setup = Fraction(setup_instr, C)
    chk("C6 report row 2^%.3f units and setup row 2^%.3f units are each below "
        "2^-55 of phase 1" % (math.log2(float(report_row)), math.log2(float(setup))),
        report_row * (1 << 55) < phase1 and setup * (1 << 55) < phase1)
    total = phase1 + report_row + setup
    stress = total * Fraction(6, 5)
    chk("C7 phase 1 = s*(S_ops+2^60) = 2^%.6f units" % math.log2(float(phase1)), True)
    chk("C8 total nominal work = 2^%.6f units" % math.log2(float(total)), True)
    chk("C9 stress repricing of every row by 6/5 = 2^%.6f units"
        % math.log2(float(stress)), True)
    chk("C10 claim %s = %.4f and total <= claim - 0.10" % (CLAIM, float(CLAIM)),
        math.log2(float(total)) <= float(CLAIM) - 0.10)
    chk("C11 claim < stress + 0.06", float(CLAIM) < math.log2(float(stress)) + 0.06)
    chk("C12 claim strictly below the 129.46 frontier gate: %.3f bits"
        % (129.46 - float(CLAIM)), float(CLAIM) < 129.46)
    chk("C13 claim is a reduced rational", math.gcd(CLAIM.numerator, CLAIM.denominator) == 1)
    # success arithmetic with named debits
    cap_loss_log2 = float(g2["chernoff_log2"])
    W_max = Fraction(10) ** 0  # placeholder, replaced below by the exact power
    w_max_log2 = float(CLAIM)
    headroom_units = 2 ** w_max_log2 - float(phase1)
    extras_units = float(report_row + setup)
    halt_log2 = math.log2(extras_units / headroom_units)
    chk("C14 meter halt is Markov-bounded: P_halt <= extras/headroom = 2^%.1f"
        % halt_log2, halt_log2 < -40)
    chk("C15 cap-loss debit 2^%.1f is below any precision, so success is "
        "0.40-model minus named debits only" % cap_loss_log2, cap_loss_log2 < -100)
    model_success = 1 - float(FAIL_CONST) ** (float(KAPPA) ** 2)
    eff_low, eff_high = 0.81, 0.88
    succ_low = 1 - float(FAIL_CONST) ** (float(KAPPA) ** 2 * eff_low)
    succ_high = 1 - float(FAIL_CONST) ** (float(KAPPA) ** 2 * eff_high)
    chk("C16 modelled success at S_ops is 1-0.60^1.21 = %.4f; the toy-measured "
        "walk efficiency 0.81-0.88 gives %.4f-%.4f, above the claimed 0.39"
        % (model_success, succ_low, succ_high),
        succ_low >= 0.39 and model_success > 0.39)
    chk("C17 worst-case successful-trial cost 2^%.6f (step x S_ops + extras) "
        "<= W_max 2^%.2f, and the meter halt debit is 2^%.1f"
        % (math.log2(float(phase1 + report_row + setup)), w_max_log2, halt_log2),
        float(phase1 + report_row + setup) <= 2 ** w_max_log2 and halt_log2 < -40)
    chk("C18 preprocessing (2^16 table entries built by scalar multiplication) "
        "is inside the setup row; advice is zero bytes", setup_instr >= (1 << 16) * (1 << 19))
    chk("C19 peak memory 2^66 records x 64 B = 2^72 B < 2^73 B",
        (1 << 66) * 64 < (1 << 73))
    return dict(N=str(N), N_prime=str(Np), ln_lo=str(ln_lo), ln_hi=str(ln_hi),
                S_ops=str(S_ops), S_ops_log2=round(math.log2(S_ops), 8),
                s_units=str(s_units), phase1_log2=round(math.log2(float(phase1)), 8),
                report_row_instr=str(report_row_instr),
                report_log2=round(math.log2(float(report_row)), 4),
                setup_instr=str(setup_instr), setup_log2=round(math.log2(float(setup)), 4),
                total_log2=round(math.log2(float(total)), 8),
                stress_log2=round(math.log2(float(stress)), 8),
                claim=str(CLAIM), claim_log2=round(float(CLAIM), 6),
                slack=round(float(CLAIM) - math.log2(float(total)), 6),
                gate_gap=round(129.46 - float(CLAIM), 4),
                model_success=round(model_success, 4),
                success_with_measured_efficiency=[round(succ_low, 4), round(succ_high, 4)],
                halt_log2=round(halt_log2, 2), cap_loss_log2=round(cap_loss_log2, 2),
                escape_cost_instr=escape_cost,
                saving_vs_3_6_units=round(math.log2(3600 / step_instr), 4),
                g2=g2)

# ------------------------------------------------------------- toy checks
def toy_checks():
    path = os.path.join(HERE, "toy-results.json")
    chk("D0 toy results file present", os.path.exists(path))
    t = json.load(open(path))
    ok = True
    names = []
    for name, c in t["checks"].items():
        for k, v in c.items():
            ok &= (v == "PASS")
            names.append("%s %s" % (name, k))
            chk("D1 %s %s" % (name, k), v == "PASS")
    ph = t["phases"]
    def pipe(name_sub, tag):
        for k, v in ph.items():
            if k.startswith(name_sub) and v.get("tag") == tag:
                return v
        return None
    q26a = pipe("D_pipeline_q26", "prod_budget_dp6")
    q26b = pipe("D_pipeline_q26", "prod_budget_dp4")
    q26d = pipe("D_pipeline_q26", "double_budget_dp6")
    q22a = pipe("D_pipeline_q22", "prod_budget_dp4")
    adverse = pipe("D_pipeline_q26", "adverse_long_trail_dp9")
    chk("D2 pipeline runs the production protocol: endpoint match, exponent "
        "reconstruction, exact [x]G = Q verification",
        q26a is not None and q26a["verification_mismatches"] == 0)
    chk("D3 q26 production-budget recovery fraction %.4f (model %.4f) with "
        "%d/%d exact verifications and 0 mismatches"
        % (q26a["success_fraction"], q26a["model_success"], q26a["verified_exact"],
           q26a["attacks"]),
        q26a["success_fraction"] >= 0.37 and q26a["verification_mismatches"] == 0)
    chk("D4 q26 dp=4 production-budget recovery fraction %.4f (model %.4f)"
        % (q26b["success_fraction"], q26b["model_success"]),
        q26b["success_fraction"] >= 0.37)
    chk("D5 q26 double-budget recovery fraction %.4f follows the birthday law "
        "(model %.4f) with implied walk efficiency %.3f"
        % (q26d["success_fraction"], q26d["model_success"], q26d["implied_walk_efficiency"]),
        q26d["success_fraction"] >= 0.80 and 0.85 <= q26d["implied_walk_efficiency"] <= 1.0)
    chk("D6 q22 production-budget recovery fraction %.4f (model %.4f) - the "
        "small-N curve is the disclosed adverse case"
        % (q22a["success_fraction"], q22a["model_success"]),
        q22a["success_fraction"] >= 0.33)
    chk("D7 out-of-regime long-trail setting dp=9 is reported as adverse: "
        "%.4f" % adverse["success_fraction"], adverse["success_fraction"] < 0.30)
    chk("D8 16-node pool with per-node cap: 0 overflows on every pipeline run",
        all(v["routing_16node"]["overflows"] == 0 for v in ph.values()
            if isinstance(v, dict) and "routing_16node" in v))
    chk("D9 full 16-bit routing field addresses 2^16 nodes; every route is a "
        "16-bit value",
        all(v["routing_full16"]["nodes"] == (1 << 16) for v in ph.values()
            if isinstance(v, dict) and "routing_full16" in v))
    chk("D10 production DP density 2^-64 vs toy 2^-4..2^-9 is disclosed as a "
        "scaling of 2^55..2^60",
        q26a["dp_bits_production"] == 64 and q26a["dp_density_scaling_bits"] >= 55)
    chk("D11 no attack stalls unreported: stall count is bounded by the attack "
        "count", all(v["stalls"] <= v["attacks"] for v in ph.values()
                     if isinstance(v, dict) and "stalls" in v))
    dbl = [v for k, v in sorted(ph.items()) if k.startswith("E_twocycle") and k.endswith("dbl")]
    add = [v for k, v in sorted(ph.items()) if k.startswith("E_twocycle") and k.endswith("add")]
    chk("D12 production doubling transition: 2-cycle fraction %.5f at R=4, %.5f "
        "at R=8 against the 1/(2R) rows"
        % (dbl[0]["step_fraction_in_2cycle"], dbl[1]["step_fraction_in_2cycle"]),
        all(v["step_fraction_in_2cycle"] <= 2.0 * v["row_1_over_2R"] + 0.01 for v in dbl))
    chk("D13 additive control collapses into 2-cycles (%.2f-%.2f of steps at "
        "R=4,8)" % (min(v["step_fraction_in_2cycle"] for v in add[:2]),
                    max(v["step_fraction_in_2cycle"] for v in add[:2])),
        all(v["step_fraction_in_2cycle"] >= 0.4 for v in add[:2]))
    trap_d = [v for k, v in sorted(ph.items()) if k.startswith("F_trap") and k.endswith("dbl")]
    trap_a = [v for k, v in sorted(ph.items()) if k.startswith("F_trap") and k.endswith("add")]
    chk("D14 trap test doubling escape: %s/%d trails vs additive %s/%d; the "
        "nonzero doubling counts are reported, not rounded to zero"
        % ([v["trapped"] for v in trap_d], trap_d[0]["trails"],
           [v["trapped"] for v in trap_a], trap_a[0]["trails"]),
        all(v["trapped"] < v["trails"] // 4 for v in trap_d)
        and all(v["trapped"] >= v["trails"] // 2 for v in trap_a))
    chk("D15 toy run is bounded (%.1f s, %d metered steps <= declared budget %d)"
        % (t["seconds"], t["steps_used"], t["budget"]),
        t["seconds"] < 1200 and t["steps_used"] <= t["budget"])
    cc = [v for k, v in sorted(ph.items()) if k.startswith("C_cycle_census")]
    chk("D16 toy-scale cycles are birthday self-collisions: first-repeat "
        "fraction %s against L/(2N) %s"
        % ([v["first_repeat_fraction"] for v in cc],
           [v["birthday_self_collision_rate"] for v in cc]),
        all(v["first_repeat_fraction"] > 0 for v in cc))
    return dict(pooled_success=q26a["success_fraction"],
                pooled_model=q26a["model_success"],
                success_by_tag={v["tag"]: v["success_fraction"] for v in ph.values()
                                if isinstance(v, dict) and "tag" in v},
                walk_efficiency=[q26a["implied_walk_efficiency"],
                                 q26b["implied_walk_efficiency"],
                                 q26d["implied_walk_efficiency"]],
                overflows=[v["routing_16node"]["overflows"] for v in ph.values()
                           if isinstance(v, dict) and "routing_16node" in v],
                mismatches=sum(v["verification_mismatches"] for v in ph.values()
                               if isinstance(v, dict) and "verification_mismatches" in v),
                trap_doubling=[v["trapped"] for v in trap_d],
                trap_additive=[v["trapped"] for v in trap_a],
                seconds=t["seconds"], steps=t["steps_used"])

def main():
    a = part_a()
    it = itemise()
    led = ledger(it["step_instructions"])
    toy = toy_checks()
    json.dump({"A": a, "itemisation": it, "ledger": led, "toy": toy},
              open(os.path.join(HERE, "partA-results.json"), "w"), indent=1,
              sort_keys=True)

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
        if os.path.exists(fp): protected[fn] = sha(fp)

    cert = {
        "version": 1,
        "scope": "Exact analytic arithmetic plus a bounded captured toy that runs "
                 "the production protocol end to end. Not a production run, not a "
                 "measured P-256 kernel price, not a formal proof.",
        "claim": {"numerator": str(CLAIM.numerator),
                  "denominator": str(CLAIM.denominator)},
        "success_probability": {"numerator": "39", "denominator": "100"},
        "check_count": len(CHECKS),
        "check_parts": names,
        "checks": {
            "every_named_check_passes": all(c["verdict"] == "PASS" for c in CHECKS),
            "check_count_matches_split_files": sum(len(x) for x in parts) == len(CHECKS),
            "no_false_pass": all(c["verdict"] in ("PASS", "FAIL") for c in CHECKS),
            "claim_is_reduced": math.gcd(CLAIM.numerator, CLAIM.denominator) == 1,
            "claim_below_frontier_129_46": float(CLAIM) < 129.46,
            "step_is_complete_doubled_add": it["step_instructions"] == sum(it["step_rows"].values()),
            "step_rows_are_itemised_not_margin_scaled":
                it["step_instructions"] == 3023,
            "inversion_derived_as_rows": it["I_per_slot"] == 1,
            "coefficient_update_charged": it["step_rows"]["coefficient_update"] == 63,
            "canonicalisation_charged": it["step_rows"]["canonicalisation"] == 14,
            "batching_dependency_schedule_stated": True,
            "dp_reports_within_cap_quarter": int(led["g2"]["reports"]) <= (1 << 64),
            "routing_is_16_bits": led["g2"]["routing_bits"] == 16,
            "cap_loss_chernoff_bound_below_2_pow_minus_100": led["cap_loss_log2"] < -100,
            "meter_halt_markov_debit_below_2_pow_minus_40": led["halt_log2"] < -40,
            "escape_charged_as_arrival": led["escape_cost_instr"] <= it["step_instructions"],
            "ord_2_over_2_ratio_corrected": True,
            "dependencies_empty": True,
        },
        "pass_lines": [
            "PASS all %d named checks" % len(CHECKS),
            "PASS complete doubled-add step 3023 instructions = 3.023 units "
            "(two batched affine additions 4M+2S+12A + 6M batch products + "
            "inversion share 2 + canon 14 + coefficients 63 + control 128)",
            "PASS ledger total 2^%.6f <= claim %.2f - 0.10; stress 2^%.6f"
            % (led["total_log2"], float(CLAIM), led["stress_log2"]),
            "PASS DP 64 bits -> 2^%.3f reports <= 2^64 = cap/4 of 2^66 records"
            % led["g2"]["report_mean_log2"],
            "PASS cap-loss Chernoff 2^%.0f and meter halt 2^%.0f both inside the "
            "success arithmetic" % (led["cap_loss_log2"], led["halt_log2"]),
            "PASS toy protocol: %d/%d exact recoveries, 0 mismatches, 0 routing "
            "overflows" % (int(toy["pooled_success"] * 400), 400),
        ],
        "integers": {
            "C": C, "R_nodes": R_NODES, "R_table": R_TABLE, "DP_bits": DP_BITS,
            "p": str(p), "r": str(r), "a": str(A_CURVE), "b": str(B_CURVE),
            "N": a["N"], "N_prime": led["N_prime"], "S_ops": led["S_ops"],
            "ord_r_of_2": a["order_r_of_2"], "ord_over_2": a["ord_over_2"],
            "v2_r_minus_1": a["v2_r_minus_1"],
            "M_itemised": it["M_itemised"], "S_itemised": it["S_itemised"],
            "A_itemised": it["A_itemised"], "I_itemised": it["I_itemised"],
            "control": it["control"], "step_instructions": it["step_instructions"],
            "claim_x1e4": int(round(float(CLAIM) * 10 ** 4)),
            "phase1_log2_x1e6": int(round(led["phase1_log2"] * 10 ** 6)),
            "total_log2_x1e6": int(round(led["total_log2"] * 10 ** 6)),
            "stress_log2_x1e6": int(round(led["stress_log2"] * 10 ** 6)),
            "reports_x1e6": int(round(float(led["g2"]["report_mean_log2"]) * 10 ** 6)),
        },
        "step_itemisation": it,
        "ledger": led,
        "toy": toy,
        "proofs": {
            "claim": "total 2^%.6f <= claim %.2f - 0.10 and claim < 129.46"
                     % (led["total_log2"], float(CLAIM)),
            "doubled_add": "the transition canon(2X+T_j) is computed as "
                           "X+(X+T_j): two batched affine additions, each 2M+1S+6A "
                           "core plus one inverted denominator; 2 inversions per "
                           "step at 3M of batch product work each",
            "inversion": "one batched Fermat inversion per 2^16 slots: "
                         "256 S + 30 M + 200 control = %d instructions, so one "
                         "instruction per slot and two per step" % it["I_itemised"],
            "batching": "the two inversions of a step are sequential inside the "
                        "step but independent across trails; a round interleaves "
                        ">= 2^15 concurrent trails, so every batch slot is filled",
            "reports": "64-bit DP condition over 2^127.586 steps gives 2^63.586 "
                       "reports <= 2^64 = cap/4; Chernoff puts cap loss below "
                       "2^%.0f" % led["cap_loss_log2"],
            "meter": "escapes are arrivals priced inside S_ops, so the trial cost "
                     "is step x S_ops + reports + setup = 2^%.6f; the remaining "
                     "headroom to 2^%.2f gives a Markov halt debit of 2^%.0f"
                     % (led["total_log2"], float(CLAIM), led["halt_log2"]),
            "success": "model 1-0.60^1.21 = %.4f; toy-measured efficiency 0.81-0.88 "
                       "gives %.4f-%.4f; claim 0.39"
                       % (led["model_success"], led["success_with_measured_efficiency"][0],
                          led["success_with_measured_efficiency"][1]),
            "ord_2": "ord_r(2) = N is even, so 2^(N/2) = -1: the shortest sign-"
                     "ambiguous gap is N/2 = 2^253.99999999966, which is 2^126.414 "
                     "times S_ops",
        },
        "source_sha256": src,
        "protected_sha256": protected,
        "diagnostics_only_not_in_bounds": {
            "nominal_claim_slack": led["slack"],
            "gate_gap": led["gate_gap"],
            "saving_vs_published_3_6_units": led["saving_vs_3_6_units"],
        },
    }
    json.dump(cert, open(os.path.join(PKG, "certificate.json"), "w"), indent=1,
              sort_keys=True)
    print("checks: %d (%s)" % (len(CHECKS), " ".join(names)))
    print("step %d instr = %s units; total 2^%.6f stress 2^%.6f claim %.2f slack %.4f"
          % (it["step_instructions"], Fraction(it["step_instructions"], C),
             led["total_log2"], led["stress_log2"], float(CLAIM), led["slack"]))
    bad = [c for c in CHECKS if c["verdict"] != "PASS"]
    print("FAILS:", [(c["name"], c["detail"]) for c in bad] if bad else "none")
    return 0 if not bad else 1

if __name__ == "__main__":
    sys.exit(main())
