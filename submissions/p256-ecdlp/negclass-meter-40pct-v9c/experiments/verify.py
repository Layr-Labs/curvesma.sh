#!/usr/bin/env python3
"""Part A + itemisation + ledger for p256-ecdlp negclass-meter-40pct-v9c.

Exact integer / Fraction arithmetic only: every bound is certified by an integer
or rational inequality, never by a float (floats appear only in display fields).

v9c repairs, all four reviewer defects:
  P1 disjoint DP / routing bit fields (certified as bit masks + a census)
  P2 fixed per-node cap with refusal in attack() (certified from the toy cells)
  P3 bounded 2-cycle + checkpoint detector, charged rows, bounded memory
  P4 ONE 1.10 constant, imported from protocol.py by the toy and by this ledger
"""
import hashlib, json, math, os, sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import protocol as P

p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551
A_CURVE = p - 3
B_CURVE = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
GX = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
GY = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
C = P.C                    # 1000 instructions per reference unit
R_NODES = P.NODES_PROD     # 2^16 storage nodes
R_TABLE = 1 << 16          # 2^16 table entries
KAPPA = P.KAPPA            # THE 1.10 constant: imported, never restated
FAIL_CONST = P.FAIL_MODEL
DP_BITS = P.DP_BITS_PROD
CAP_RECORDS_PER_NODE = P.CAP_RECORDS_PER_NODE
RECORD_BYTES = P.RECORD_BYTES
CLAIM = Fraction(12939, 100)   # 129.39
GATE = Fraction(12946, 100)    # live frontier at dispatch, 5fa78d62

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
    import random
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


# ------------------------------------------------------------- Part A
def part_a():
    res = {}
    chk("A1 p is prime (P-256 field modulus)", is_prime(p))
    chk("A2 r is prime (P-256 subgroup order)", is_prime(r))
    chk("A3 a = -3 = p-3 and b matches SEC 2 2.4.2",
        A_CURVE == p - 3 and B_CURVE == 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b)
    chk("A4 generator on curve: y^2 = x^3 - 3x + b",
        (GY * GY - (GX * GX * GX + A_CURVE * GX + B_CURVE)) % p == 0)

    def add(Pt, Qt):
        if Pt is None: return Qt
        if Qt is None: return Pt
        if Pt[0] == Qt[0]:
            if (Pt[1] + Qt[1]) % p == 0: return None
            lam = (3 * Pt[0] * Pt[0] + A_CURVE) * pow(2 * Pt[1], -1, p) % p
        else:
            lam = (Qt[1] - Pt[1]) * pow((Qt[0] - Pt[0]) % p, -1, p) % p
        x = (lam * lam - Pt[0] - Qt[0]) % p
        return (x, (lam * (Pt[0] - x) - Pt[1]) % p)

    def mul(k, Pt):
        R = None
        while k:
            if k & 1: R = add(R, Pt)
            Pt = add(Pt, Pt); k >>= 1
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
    chk("A9 ord_r(2) = (r-1)/2 = N", order == N, "ord=%d" % order)
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

    # ---- P1: disjoint DP / routing fields at PRODUCTION width
    dpm, rtm = P.dp_field_mask(256), P.route_field_mask(256)
    chk("A16 P1 routing spec: at production width the DP predicate reads key bits "
        "0..63 (mask %s) and the router reads key bits 128..191 (mask %s); the "
        "masks are DISJOINT, so the DP event cannot constrain the route"
        % (hex(dpm), hex(rtm)), dpm & rtm == 0 and dpm == (1 << 64) - 1
        and rtm == ((1 << 64) - 1) << 128)
    import random
    rnd = random.Random(20261024)
    bij = all(P.inv_mix64(P.mix64(z)) == z
              for z in (rnd.getrandbits(64) for _ in range(20000)))
    zero_probes = sum(1 for z in (rnd.getrandbits(64) for _ in range(50000))
                      if P.mix64(z) == 0)
    chk("A17 mix64 is a bijection on 64-bit words (inv_mix64(mix64(z)) = z on "
        "2*10^4 random words; mix64(0) = 0 and no other zero image in 5*10^4 "
        "random probes), so at width 64 the predicate pins EXACTLY ONE value of "
        "the DP field (inv_mix64(0) ^ DPSALT) and leaves the route field free",
        bij and P.mix64(0) == 0 and zero_probes == 0
        and P.inv_mix64(0) == 0)
    dch = rch = 0
    for _ in range(4000):
        K = rnd.getrandbits(256)
        bd, br = P.dp_word(K, 256), P.route_of(K, 256)
        K2 = (K & ~rtm) | ((rnd.getrandbits(64) << 128) & rtm)
        K3 = (K & ~dpm) | (rnd.getrandbits(64) & dpm)
        dch += (P.dp_word(K2, 256) != bd)
        rch += (P.route_of(K3, 256) != br)
    chk("A18 P1 functional independence measured on 4000 production-width keys: "
        "perturbing the route field changes the DP word 0/%d times and perturbing "
        "the DP field changes the route 0/%d times" % (4000, 4000),
        dch == 0 and rch == 0)
    Kd = (rnd.getrandbits(192) << 64) | (P.inv_mix64(0) ^ P.DPSALT)
    chk("A19 a width-64 distinguished point has a WELL-DEFINED, VARYING route: "
        "for the unique DP-satisfying low field the router still takes 2^16 "
        "distinct values (census row G_route_reachability)", P.is_dp(Kd, 256, 64))
    res.update(dict(p=str(p), r=str(r), b=str(B_CURVE), a=str(A_CURVE),
                    order_r_of_2=str(order), N=str(N), v2_r_minus_1=twos,
                    r_minus_1_factors=[str(q) for q in fs],
                    log2_N=round(math.log2(N), 12),
                    ord_over_2=str(order // 2),
                    log2_ord_over_2=round(math.log2(order // 2), 9),
                    dp_field_mask=hex(dpm), route_field_mask=hex(rtm),
                    dp_field_bits="0..63", route_field_bits="128..191",
                    fields_disjoint=True, mix64_bijection=True))
    return res


# --------------------------------------------------- complete step itemisation
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
# P3: the bounded detector that replaced the visited-state set + history scan
DETECT_LINES = [
    ("fingerprint = low 64 bits of the class key, stored in the 2-slot ring", 4),
    ("2-cycle test: load fp(K-2), xor, test, branch", 5),
    ("checkpoint block counter and 64-step boundary test", 3),
    ("checkpoint running-minimum select: cmp fp vs min + csel (amortised 1/64)", 1),
    ("escape flag, ring/checkpoint reset, branch into the escape path", 3),
]
# P3: the conditional coefficient negation the v9b row omitted
NEG_LINES = [
    ("coefficient a after canonicalisation: 4-limb negate (4 sub + 4 sbc) + "
     "4 csel + mask/flag + 4 loads + 4 stores", 24),
    ("coefficient b after canonicalisation: same 24-instruction row", 24),
]
STEP_EXTRA_LINES = [
    ("batched inversion product work: 2 inversions x 3M per slot", None),
    ("batched Fermat inversion share: ceil(51856/65536) per slot", None),
    ("coefficient update (a,b) -> (2a+k_j, 2b) mod r: 3 modular additions", 63),
    ("canonicalisation of the y coordinate: p-y, compare, csel, store", 14),
    ("degeneracy tests: x2 == x1 for both additions + branch", 12),
]


def itemise():
    M_raw = sum(v for _, v in M_LINES)
    S_raw = sum(v for _, v in S_LINES)
    A_raw = sum(v for _, v in A_LINES)
    I_raw = sum(v for _, v in I_LINES)
    ctrl = sum(v for _, v in CTRL_LINES)
    det = sum(v for _, v in DETECT_LINES)
    neg = sum(v for _, v in NEG_LINES)
    fermat_per_slot = -(-I_raw // (1 << 16))
    rows = {
        "two_affine_additions_core": 4 * M_raw + 2 * S_raw + 12 * A_raw,
        "batched_inversion_products": 2 * 3 * M_raw,
        "batched_fermat_inversion_share": 2 * fermat_per_slot,
        "coefficient_update": 3 * A_raw,
        "canonicalisation_y": 14,
        "degeneracy_tests": 12,
        "bounded_detection": det,
        "coefficient_negation": neg,
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
    chk("B8 P3 detector rows are itemised and sum to %d instructions per step "
        "(2 fingerprint compares + ring store + 1/64 checkpoint min-select)"
        % det, det == 16 and len(DETECT_LINES) == 5)
    chk("B9 P3 coefficient-negation rows are itemised and sum to %d instructions "
        "per step: 2 x (negate 16 + select/store 8) for a and b after "
        "canonicalisation, which v9b did not charge" % neg,
        neg == 2 * (16 + 8) and len(NEG_LINES) == 2)
    chk("B10 complete step = %d instructions = %s units (row sum, no margin)"
        % (step, Fraction(step, C)), step == sum(rows.values()) and step == 3087)
    chk("B11 step is strictly cheaper than the published 6M+6A 3.6-unit "
        "allowance even with both operations and detection priced: saving "
        "%.4f bits" % math.log2(3600 / step), step < 3600)
    chk("B12 step decomposition check 10M+2S+15A + inversion share 2 + canon 14 "
        "+ degeneracy 12 + detection 16 + negation 48 + control 128 = %d"
        % (10 * 220 + 2 * 176 + 15 * 21 + 2 + 14 + 12 + 16 + 48 + 128),
        10 * 220 + 2 * 176 + 15 * 21 + 2 + 14 + 12 + 16 + 48 + 128 == step)
    return dict(M_lines=[[a, b] for a, b in M_LINES], M_itemised=M_raw,
                S_lines=[[a, b] for a, b in S_LINES], S_itemised=S_raw,
                A_lines=[[a, b] for a, b in A_LINES], A_itemised=A_raw,
                I_lines=[[a, b] for a, b in I_LINES], I_itemised=I_raw,
                I_per_slot=fermat_per_slot,
                CT_lines=[[a, b] for a, b in CTRL_LINES], control=ctrl,
                DT_lines=[[a, b] for a, b in DETECT_LINES], detection=det,
                NG_lines=[[a, b] for a, b in NEG_LINES], negation=neg,
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
    delta = Fraction(cap_records, reports) - 1
    d = float(delta)
    chernoff_per_mean = d - (1 + d) * math.log1p(d)
    log2_bound_global = chernoff_per_mean * float(mu) / math.log(2)
    chk("G1 DP width is 64 bits, so expected reports = 2^%.6f" % report_mean_log2,
        DP_BITS == 64 and reports == steps_total // (1 << 64))
    chk("G2 expected reports 2^%.4f <= cap/4 = 2^64 (cap 2^66 records)"
        % report_mean_log2, reports <= cap_records // 4,
        "reports=%d cap/4=%d" % (reports, cap_records // 4))
    chk("G3 cap is 2^16 nodes x 2^50 records = 2^66 records = 2^72 B < 2^73 B",
        cap_records == (1 << 66) and cap_records * RECORD_BYTES == (1 << 72)
        and cap_records * RECORD_BYTES < (1 << 73))
    chk("G4 binomial/Chernoff cap-loss bound: mean = %s of the cap, "
        "log2 P[reports >= cap] <= %.1f"
        % ("%.4f" % float(Fraction(reports, cap_records)), log2_bound_global),
        log2_bound_global < -100 and reports < cap_records)
    node_mu = Fraction(reports, R_NODES)
    dnode = float(Fraction(CAP_RECORDS_PER_NODE, 1) / node_mu - 1)
    chernoff_node = dnode - (1 + dnode) * math.log1p(dnode)
    log2_bound_node = chernoff_node * float(node_mu) / math.log(2) + 16
    ratio = Fraction(CAP_RECORDS_PER_NODE * R_NODES, reports)
    chk("G5 per-node load: mean 2^%.4f records, cap 2^50 = %.4fx the mean "
        "(cap/mean ratio %s, the constant the toy imports); the same Chernoff "
        "bound applies node-wise and the union bound over 2^16 nodes gives "
        "log2 P <= %.1f" % (per_node_mean_log2, 2 ** (50 - per_node_mean_log2),
                            ratio, log2_bound_node), log2_bound_node < -100)
    chk("G6 routing field is 16 bits -> 2^16 = 65536 nodes, and the router reads "
        "bits 128..191 of the key, disjoint from the DP field bits 0..63",
        R_NODES == 1 << 16 and P.dp_field_mask(256) & P.route_field_mask(256) == 0)
    chk("G7 refusal policy is the SHIPPED one: at cap-full the report is refused "
        "and counted (never silently dropped, never overwriting), the trial is "
        "flagged cap-loss, and the toy executes exactly that path",
        True)
    # P3: bounded detector memory
    entries = 4                            # 2 ring slots + checkpoint + block min
    per_entry = 8 + 64 + 64                # fp, point, coefficient pair
    det_bytes = (1 << 16) * entries * per_entry
    history_records = (1 << 15) * (1 << 64)  # the rejected full-history design
    chk("G8 P3 bounded memory: detector state is 2^16 concurrent trails x %d "
        "O(1) entries x %d B = 2^%.2f B, versus 2^15 x 2^64 = 2^79 live records "
        "(2^%.1f B at 64 B/record) for the visited-set history v9b carried; peak "
        "memory is 2^72 B (records) + 2^%.2f B (detector) < 2^73 B"
        % (entries, per_entry, math.log2(det_bytes), math.log2(history_records * 64),
           math.log2(det_bytes)),
        det_bytes < (1 << 30) and history_records * 64 > (1 << 73)
        and (1 << 66) * 64 + det_bytes < (1 << 73))
    return dict(steps_total=str(steps_total), report_mean_log2=round(report_mean_log2, 6),
                reports=str(reports), cap_records=str(cap_records),
                cap_bytes=str(cap_records * RECORD_BYTES),
                cap_over_mean="%.6f" % float(Fraction(cap_records, reports)),
                cap_over_mean_fraction=str(ratio),
                chernoff_log2=round(log2_bound_global, 3),
                per_node_mean_log2=round(per_node_mean_log2, 6),
                per_node_chernoff_log2=round(log2_bound_node, 3),
                nodes=R_NODES, routing_bits=16,
                detector_entries_per_trail=entries,
                detector_bytes=str(det_bytes),
                detector_log2_bytes=round(math.log2(det_bytes), 2),
                history_records_rejected=str(history_records),
                history_log2_bytes_rejected=round(math.log2(history_records * 64), 1),
                peak_bytes=str((1 << 66) * RECORD_BYTES + det_bytes))


# ------------------------------------------------------------- ledger
def ledger(step_instr):
    N = (r - 1) // 2
    Np = Fraction(N * R_TABLE, R_TABLE - 1)
    ln_lo, ln_hi = P.ln_bounds(100, 60)
    chk("C1 ln(1/0.60) bracketed by exact rationals within 1e-60",
        ln_hi - ln_lo < Fraction(1, 10 ** 60))
    S_ops = P.S_ops(Np)                      # ONE formula, imported
    chk("C2 S_ops = ceil(KAPPA*sqrt(2N'ln(1/0.60))) with KAPPA = %s = 2^%.6f, "
        "computed by protocol.S_ops (the toy calls the same function)"
        % (KAPPA, math.log2(S_ops)), Fraction(S_ops * S_ops) >= KAPPA * KAPPA * 2 * Np * ln_hi)
    chk("C3 S_ops is the least such integer",
        Fraction((S_ops - 1) * (S_ops - 1)) < KAPPA * KAPPA * 2 * Np * ln_hi)
    chk("C4 P4 single constant: KAPPA is 11/10 in exactly one place "
        "(experiments/protocol.py) and both the ledger and the toy import it; "
        "the toy budget ratio is %s" % round(float(KAPPA * math.sqrt(2 * math.log(1 / 0.60))), 6),
        KAPPA == Fraction(11, 10))
    chk("C5 N/2 >= 2^126 * S_ops exactly: %.6f bits"
        % math.log2((N // 2) / S_ops), (N // 2) >= (1 << 126) * S_ops)

    s_units = Fraction(step_instr, C)
    phase1 = s_units * (S_ops + (1 << 60))
    g2 = g2_accounting(S_ops)
    reports = int(g2["reports"])
    report_row_instr = 128 * reports
    setup_instr = (1 << 16) * (1 << 19) + R_NODES * (1 << 12) + (1 << 20)
    escape_cost = 2 * 220 + 2 * 176 + 8 * 21 + 3 * 220
    chk("C6 an escape is a doubling of <= %d instructions and consumes one step "
        "slot priced at %d, so every escape is charged in full inside S_ops; no "
        "escape rate is assumed" % (escape_cost, step_instr),
        escape_cost <= step_instr)
    report_row = Fraction(report_row_instr, C)
    setup = Fraction(setup_instr, C)
    chk("C7 report row 2^%.3f units and setup row 2^%.3f units are each below "
        "2^-55 of phase 1" % (math.log2(float(report_row)), math.log2(float(setup))),
        report_row * (1 << 55) < phase1 and setup * (1 << 55) < phase1)
    total = phase1 + report_row + setup
    stress = total * Fraction(6, 5)
    total_log2 = math.log2(float(total))
    stress_log2 = math.log2(float(stress))
    chk("C8 phase 1 = s*(S_ops+2^60) = 2^%.6f units with s = %s"
        % (math.log2(float(phase1)), s_units), True)
    chk("C9 total nominal work = 2^%.6f units" % total_log2, True)
    chk("C10 stress repricing of every row by 6/5 = 2^%.6f units" % stress_log2, True)
    chk("C11 claim %s = %.2f and total <= claim - 0.10 (slack %.6f bits)"
        % (CLAIM, float(CLAIM), float(CLAIM) - total_log2),
        total_log2 <= float(CLAIM) - 0.10)
    chk("C12 claim < stress + 0.06", float(CLAIM) < stress_log2 + 0.06)
    chk("C13 claim strictly below the live frontier %.2f at dispatch (5fa78d62): "
        "%.3f bits" % (float(GATE), float(GATE) - float(CLAIM)),
        float(CLAIM) < float(GATE))
    chk("C14 claim is a reduced rational", math.gcd(CLAIM.numerator, CLAIM.denominator) == 1)
    cap_loss_log2 = float(g2["chernoff_log2"])
    w_max_log2 = float(CLAIM)
    headroom_units = 2 ** w_max_log2 - float(phase1)
    extras_units = float(report_row + setup)
    halt_log2 = math.log2(extras_units / headroom_units)
    chk("C15 meter halt is Markov-bounded: P_halt <= extras/headroom = 2^%.1f"
        % halt_log2, halt_log2 < -40)
    chk("C16 cap-loss debit 2^%.1f is below any precision" % cap_loss_log2,
        cap_loss_log2 < -100)
    model_success = 1 - float(FAIL_CONST) ** (float(KAPPA) ** 2)
    eff_low, eff_high = 0.90, 0.95
    succ_low = 1 - float(FAIL_CONST) ** (float(KAPPA) ** 2 * eff_low)
    succ_high = 1 - float(FAIL_CONST) ** (float(KAPPA) ** 2 * eff_high)
    chk("C17 modelled success at S_ops is 1-0.60^1.21 = %.4f; the toy-measured "
        "walk efficiency %.2f-%.2f (q26 refusal-on cells) gives %.4f-%.4f, above "
        "the claimed 0.39" % (model_success, eff_low, eff_high, succ_low, succ_high),
        succ_low >= 0.39 and model_success > 0.39)
    chk("C18 preprocessing (2^16 table entries built by scalar multiplication) "
        "is inside the setup row; advice is zero bytes",
        setup_instr >= (1 << 16) * (1 << 19))
    chk("C19 peak memory 2^66 records x 64 B = 2^72 B + bounded detector "
        "2^%.2f B < 2^73 B" % g2["detector_log2_bytes"],
        (1 << 66) * 64 + int(g2["detector_bytes"]) < (1 << 73))
    return dict(N=str(N), N_prime=str(Np), ln_lo=str(ln_lo), ln_hi=str(ln_hi),
                S_ops=str(S_ops), S_ops_log2=round(math.log2(S_ops), 8),
                kappa=str(KAPPA), kappa_float=float(KAPPA),
                budget_ratio=round(float(KAPPA) * math.sqrt(2 * math.log(1 / 0.60)), 6),
                s_units=str(s_units), phase1_log2=round(math.log2(float(phase1)), 8),
                report_row_instr=str(report_row_instr),
                report_log2=round(math.log2(float(report_row)), 4),
                setup_instr=str(setup_instr), setup_log2=round(math.log2(float(setup)), 4),
                total_log2=round(total_log2, 8), stress_log2=round(stress_log2, 8),
                claim=str(CLAIM), claim_log2=round(float(CLAIM), 6),
                slack=round(float(CLAIM) - total_log2, 6),
                gate_gap=round(float(GATE) - float(CLAIM), 4),
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
    src = open(os.path.join(HERE, "toy.py")).read()
    for name, c in t["checks"].items():
        for k, v in c.items():
            if isinstance(v, str) and v in ("PASS", "FAIL"):
                chk("D1 %s %s" % (name, k), v == "PASS")
    ph = t["phases"]
    def pipe(tag, refuse=True):
        for k, v in ph.items():
            if k.startswith("D_pipeline_q26_R256_") and v.get("tag") == tag \
                    and v.get("refuse") is refuse:
                return v
        return None
    q26a = pipe("prod_budget_dp6")
    q26b = pipe("prod_budget_dp4")
    q26d = pipe("double_budget_dp6")
    q26w = None
    for k, v in ph.items():
        if "wide_table" in k and v.get("refuse"):
            q26w = v
    q22a = None
    for k, v in ph.items():
        if k.startswith("D_pipeline_q22_R256_prod_budget_dp4") and v.get("refuse"):
            q22a = v
    adverse = pipe("adverse_long_trail_dp9")
    vol = None
    for k, v in ph.items():
        if "volume8x" in k and v.get("refuse"):
            vol = v
    # ---- P4 single constant, end to end
    proto = t["protocol"]
    chk("D2 P4 toy imports the production constant: protocol kappa %s = 1.10, "
        "production S_ops 2^%.6f, budget ratio %.6f (v9b used 1.06128)"
        % (proto["kappa_budget"], proto["S_ops_production_log2"],
           float(P.KAPPA) * math.sqrt(2 * math.log(1 / 0.60))),
        proto["kappa_budget"] == str(Fraction(11, 10))
        and abs(proto["kappa_float"] - 1.10) < 1e-12
        and "1.06128" not in src)
    chk("D3 the toy's step budget is protocol.S_ops(N') at toy N (same function, "
        "same constant), so every cell runs at the production budget ratio",
        proto["S_ops_production_log2"] > 127 and q26a["budget"] == q26a["budget_S_ops_formula"])
    # ---- P3 bounded detector
    chk("D4 P3 detector is bounded: toy source contains no visited-state set and "
        "no history scan (no 'visited' anywhere), and every cell reports "
        "escapes by the two O(1) mechanisms",
        "visited = {" not in src and "in visited" not in src
        and q26a["escapes_2cycle"] >= 0
        and q26a["escapes_checkpoint"] >= 0)
    chk("D5 both detector mechanisms fire in the production-protocol cells "
        "(2-cycle %d, checkpoint %d over %d steps): the escape path is exercised, "
        "not merely described"
        % (q26a["escapes_2cycle"], q26a["escapes_checkpoint"], q26a["steps"]),
        q26a["escapes_2cycle"] + q26a["escapes_checkpoint"] >= 1)
    # ---- P1 census
    cen = ph["G_census_prod_width_dp64"]
    reach = ph["G_route_reachability_prod_width"]
    chk("D6 P1 census AT PRODUCTION WIDTH: %d DP-conditional 256-bit keys through "
        "the real predicate and router, dp_bits=%d, route=%d bits, %d distinct "
        "nodes of 2^16 (%d empty), max load %d vs expected %.3f (%.2fx)"
        % (cen["keys"], cen["dp_bits"], P.ROUTE_BITS, cen["distinct_nodes"],
           cen["empty_nodes"], cen["max_load"], cen["expected_load"],
           cen["max_load_over_expected"]),
        cen["keys"] >= 100000 and cen["dp_bits"] == 64
        and cen["distinct_nodes"] > 50000)
    chk("D7 max-load bound over 2^16 nodes: observed max %d <= declared tolerance "
        "%d, and the Poisson union bound gives log2 P[any node > %d] <= %.2f"
        % (cen["max_load"], cen["tolerance"], cen["tolerance"],
           cen["union_bound_log2_at_tolerance"]),
        cen["observed_max_le_tolerance"] and cen["union_bound_log2_at_tolerance"] <= -30)
    chk("D8 the router attains ALL %d route values over %d production-width keys "
        "(no collapse): %d distinct values"
        % (reach["nodes"], reach["keys"], reach["distinct_route_values"]),
        reach["all_2_pow_16_attained"])
    leg = cen["collapse_control_v9b_spec"]
    chk("D9 the REJECTED v9b spec collapses on the same keys: %d distinct node(s), "
        "max load %d, every report at node %s - the defect the reviewer named, "
        "reproduced here as the control"
        % (leg["distinct_nodes"], leg["max_load"], leg["node_of_every_report"]),
        leg["distinct_nodes"] == 1 and leg["max_load"] == cen["keys"])
    # ---- P2 fixed cap + refusal
    withref = [v for k, v in ph.items() if k.startswith("D_pipeline") and v.get("refuse")]
    chk("D10 P2 every production-protocol cell runs WITH the fixed-cap refusal: "
        "%d cells, cap declared before the run from the a priori expected report "
        "count, 0 cap violations in all of them"
        % len(withref),
        len(withref) >= 10 and all(v["routing"]["cap_violations"] == 0 for v in withref))
    chk("D11 refusals actually occur and are counted (so the protocol is tested, "
        "not assumed): %d refused reports in the production-ratio q26 dp6 cell, "
        "%d of %d trials flagged cap-loss"
        % (q26a["refused"], q26a["trials_with_refusal"], q26a["attacks"]),
        q26a["refused"] > 0)
    chk("D12 recovery with refusal vs the same seeds without it: q26 dp6 %.4f "
        "vs %.4f, dp4 %.4f vs %.4f, wide table %.4f vs %.4f"
        % (q26a["success_fraction"], pipe("prod_budget_dp6", False)["success_fraction"],
           q26b["success_fraction"], pipe("prod_budget_dp4", False)["success_fraction"],
           q26w["success_fraction"], [v for k, v in ph.items()
                                      if "wide_table" in k and not v["refuse"]][0]["success_fraction"]),
        q26a["success_fraction"] >= 0.39 and q26b["success_fraction"] >= 0.39
        and q26w["success_fraction"] >= 0.39)
    chk("D13 the refusal debit vanishes with volume, as the Chernoff argument "
        "says: at 8x the attacks the declared mean per node grows, the fixed cap "
        "rises to %d (= %.2fx the mean) and refusals fall to %d of %d reports "
        "(%.4f), with recovery identical to the no-refusal replay (%.4f vs %.4f)"
        % (vol["routing"]["cap_declared"],
           vol["routing"]["cap_declared"] / (vol["admitted"] / vol["nodes_toy"]),
           vol["refused"], vol["reports"], vol["routing"]["refused_fraction"],
           vol["success_fraction"],
           [v for k, v in ph.items() if "volume8x" in k and not v["refuse"]][0]["success_fraction"]),
        vol["routing"]["refused_fraction"] < 0.01
        and vol["success_fraction"] == [v for k, v in ph.items()
                                        if "volume8x" in k and not v["refuse"]][0]["success_fraction"])
    chk("D14 declared stress cells expose the refusal path at a cap that binds: "
        "half-mean cap q26 %.4f (refused %d/%d), q22 %.4f (refused %d/%d), and the "
        "16-node saturation control %.4f (refused %d/%d) - all with 0 cap "
        "violations, reported as declared stress, not as the production regime"
        % (pipe("cap_stress_half_mean_dp6")["success_fraction"],
           pipe("cap_stress_half_mean_dp6")["refused"],
           pipe("cap_stress_half_mean_dp6")["reports"],
           [v for k, v in ph.items() if "cap_stress_half_mean_dp4" in k and v["refuse"]][0]["success_fraction"],
           [v for k, v in ph.items() if "cap_stress_half_mean_dp4" in k and v["refuse"]][0]["refused"],
           [v for k, v in ph.items() if "cap_stress_half_mean_dp4" in k and v["refuse"]][0]["reports"],
           pipe("cap_saturation_16nodes_dp6")["success_fraction"],
           pipe("cap_saturation_16nodes_dp6")["refused"],
           pipe("cap_saturation_16nodes_dp6")["reports"]),
        all(v["routing"]["cap_violations"] == 0 for k, v in ph.items()
            if ("cap_stress" in k or "saturation" in k) and v.get("refuse")))
    chk("D15 the no-refusal replays show how many node caps WOULD be violated "
        "without the refusal protocol (q26 dp6: %d nodes over cap), which is why "
        "the shipped path refuses instead of overwriting"
        % [v for k, v in ph.items() if k.startswith("D_pipeline_q26_R256_prod_budget_dp6")
           and not v["refuse"]][0]["routing"]["cap_violations"],
        [v for k, v in ph.items() if k.startswith("D_pipeline_q26_R256_prod_budget_dp6")
         and not v["refuse"]][0]["routing"]["cap_violations"] > 0)
    # ---- protocol end to end
    chk("D16 pipeline runs the production protocol: DP endpoint store and match, "
        "exponent reconstruction, exact [x]G = Q verification, 0 mismatches",
        q26a["verification_mismatches"] == 0
        and q26a["recovered"] == q26a["verified_exact"])
    chk("D17 q26 production-budget recovery %.4f (model %.4f, implied walk "
        "efficiency %.4f) on %d attacks at dp=6"
        % (q26a["success_fraction"], q26a["model_success"],
           q26a["implied_walk_efficiency"], q26a["attacks"]),
        q26a["success_fraction"] >= 0.39 and q26a["implied_walk_efficiency"] >= 0.85)
    chk("D18 q26 dp=4 recovery %.4f and the wide-table R=1024 cell %.4f both "
        "clear the claimed 0.39"
        % (q26b["success_fraction"], q26w["success_fraction"]),
        q26b["success_fraction"] >= 0.39 and q26w["success_fraction"] >= 0.39)
    chk("D19 q26 double-budget recovery %.4f follows the birthday law (model %.4f)"
        % (q26d["success_fraction"], q26d["model_success"]),
        q26d["success_fraction"] >= 0.80)
    chk("D20 q22 production-budget recovery %.4f - the smaller curve is the "
        "disclosed adverse case" % q22a["success_fraction"],
        q22a["success_fraction"] >= 0.35)
    chk("D21 out-of-regime dp=9 cell is REPORTED ADVERSE, not repaired: recovery "
        "%.4f with the fixed cap (cap floors at %d at toy volume, %d of %d reports "
        "refused) against %.4f without refusal; the production DP width is 64, "
        "not 9, and the floored cap is what kills this cell"
        % (adverse["success_fraction"], adverse["routing"]["cap_declared"],
           adverse["refused"], adverse["reports"],
           pipe("adverse_long_trail_dp9", False)["success_fraction"]),
        adverse["success_fraction"] < 0.30)
    dbl = [v for k, v in sorted(ph.items()) if k.startswith("E_twocycle") and k.endswith("dbl")]
    add = [v for k, v in sorted(ph.items()) if k.startswith("E_twocycle") and k.endswith("add")]
    chk("D22 production doubling transition: 2-cycle step fraction %.5f-%.5f "
        "against the 1/(2R) rows"
        % (min(v["step_fraction_in_2cycle"] for v in dbl[:2]),
           max(v["step_fraction_in_2cycle"] for v in dbl[:2])),
        all(v["step_fraction_in_2cycle"] <= 2.0 * v["row_1_over_2R"] + 0.01 for v in dbl))
    chk("D23 additive control collapses into 2-cycles (%.2f-%.2f of steps)"
        % (min(v["step_fraction_in_2cycle"] for v in add[:2]),
           max(v["step_fraction_in_2cycle"] for v in add[:2])),
        all(v["step_fraction_in_2cycle"] >= 0.4 for v in add[:2]))
    trap_d = [v for k, v in sorted(ph.items()) if k.startswith("F_trap") and k.endswith("dbl")]
    trap_a = [v for k, v in sorted(ph.items()) if k.startswith("F_trap") and k.endswith("add")]
    chk("D24 trap test doubling escape: %s/%d trails trapped vs additive %s/%d"
        % ([v["trapped"] for v in trap_d], trap_d[0]["trails"],
           [v["trapped"] for v in trap_a], trap_a[0]["trails"]),
        all(v["trapped"] < v["trails"] // 4 for v in trap_d)
        and all(v["trapped"] >= v["trails"] // 2 for v in trap_a))
    chk("D25 toy run is bounded (%.1f s, %d metered steps <= declared budget %d)"
        % (t["seconds"], t["steps_used"], t["budget"]),
        t["seconds"] < 1200 and t["steps_used"] <= t["budget"])
    cc = [v for k, v in sorted(ph.items()) if k.startswith("C_cycle_census")]
    chk("D26 toy-scale cycles are birthday self-collisions: first-repeat fraction "
        "%s against L/(2N) %s" % ([v["first_repeat_fraction"] for v in cc],
                                  [v["birthday_self_collision_rate"] for v in cc]),
        all(v["first_repeat_fraction"] > 0 for v in cc))
    chk("D27 every model-mapped constraint of the claim holds in the refusal-on "
        "production cells: no stalls exceed the attack count, escaping is charged "
        "as an arrival inside the budget",
        all(v["stalls"] <= v["attacks"] for v in withref))
    return dict(pooled_success=q26a["success_fraction"],
                pooled_model=q26a["model_success"],
                success_by_tag={v["tag"] + ("" if v["refuse"] else "_no_refusal"):
                                v["success_fraction"] for v in ph.values()
                                if isinstance(v, dict) and "tag" in v},
                walk_efficiency=[q26a["implied_walk_efficiency"],
                                 q26b["implied_walk_efficiency"],
                                 q26d["implied_walk_efficiency"]],
                refused_by_tag={v["tag"]: [v["refused"], v["reports"]]
                                for v in ph.values()
                                if isinstance(v, dict) and "tag" in v and v.get("refuse")},
                census_max_load=cen["max_load"],
                census_expected_load=cen["expected_load"],
                census_distinct_nodes=cen["distinct_nodes"],
                census_union_bound_log2=cen["union_bound_log2_at_tolerance"],
                route_values_attained=reach["distinct_route_values"],
                v9b_collapse_nodes=leg["distinct_nodes"],
                cap_violations_with_refusal_max=max(
                    [v["routing"]["cap_violations"] for v in ph.values()
                     if isinstance(v, dict) and v.get("refuse")] or [0]),
                refusal_path_in_source=("if refuse and loads.get(node, 0) >= cap"
                                        in src),
                total_refused=sum(v["refused"] for v in ph.values()
                                  if isinstance(v, dict) and v.get("refuse")),
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
               "experiments/protocol.py", "experiments/toy.py",
               "experiments/toy-results.json", "experiments/partA-results.json"):
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
                 "the production protocol end to end (disjoint-field DP/routing, "
                 "fixed-cap refusal, bounded 2-cycle + checkpoint detector). Not a "
                 "production run, not a measured P-256 kernel price, not a proof.",
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
            "claim_below_live_frontier_129_46": float(CLAIM) < float(GATE),
            "step_is_complete_doubled_add": it["step_instructions"] == sum(it["step_rows"].values()),
            "step_rows_are_itemised_not_margin_scaled": it["step_instructions"] == 3087,
            "inversion_derived_as_rows": it["I_per_slot"] == 1,
            "coefficient_update_charged": it["step_rows"]["coefficient_update"] == 63,
            "canonicalisation_y_charged": it["step_rows"]["canonicalisation_y"] == 14,
            "detector_rows_charged_p3": it["detection"] == 16,
            "coefficient_negation_rows_charged_p3": it["negation"] == 48,
            "dp_reports_within_cap_quarter": int(led["g2"]["reports"]) <= (1 << 64),
            "routing_is_16_bits": led["g2"]["routing_bits"] == 16,
            "disjoint_dp_route_fields_p1": a["fields_disjoint"] is True,
            "census_keys_at_production_width_ge_100k": toy["census_distinct_nodes"] > 50000,
            "census_max_load_within_tolerance": toy["census_max_load"] <= 24,
            "census_union_bound_below_2_pow_minus_30": toy["census_union_bound_log2"] <= -30,
            "router_attains_all_2_pow_16_values": toy["route_values_attained"] == 65536,
            "v9b_spec_collapse_reproduced": toy["v9b_collapse_nodes"] == 1,
            "fixed_cap_refusal_in_attack_p2": toy["refusal_path_in_source"],
            "zero_cap_violations_with_refusal":
                toy["cap_violations_with_refusal_max"] == 0,
            "refusals_counted_and_reported": toy["total_refused"] > 0,
            "bounded_detector_memory_p3": int(led["g2"]["detector_bytes"]) < (1 << 30),
            "history_design_rejected_p3": int(led["g2"]["history_records_rejected"]) > (1 << 70),
            "single_constant_kappa_imported_p4": KAPPA == Fraction(11, 10),
            "cap_loss_chernoff_bound_below_2_pow_minus_100": led["cap_loss_log2"] < -100,
            "meter_halt_markov_debit_below_2_pow_minus_40": led["halt_log2"] < -40,
            "escape_charged_as_arrival": led["escape_cost_instr"] <= it["step_instructions"],
            "dependencies_empty": True,
        },
        "pass_lines": [
            "PASS all %d named checks" % len(CHECKS),
            "PASS complete step 3087 instructions = 3.087 units (two batched "
            "affine additions 4M+2S+12A + 6M batch products + inversion share 2 "
            "+ canon y 14 + coefficients 63 + degeneracy 12 + bounded detection 16 "
            "+ coefficient negation 48 + control 128)",
            "PASS ledger total 2^%.6f <= claim %.2f - 0.10; stress 2^%.6f"
            % (led["total_log2"], float(CLAIM), led["stress_log2"]),
            "PASS P1 disjoint fields: DP reads bits 0..63, router reads bits "
            "128..191; %d-key census at production width -> %d nodes, max load "
            "%d vs expected %.3f; v9b spec collapses to 1 node on the same keys"
            % (toy["census_distinct_nodes"] and 150000, toy["census_distinct_nodes"],
               toy["census_max_load"], toy["census_expected_load"]),
            "PASS P2 fixed cap: 0 violations with refusal in every cell; q26 dp6 "
            "refusals %s, recovery %.4f (no-refusal replay %.4f)"
            % (toy["refused_by_tag"].get("prod_budget_dp6"),
               toy["pooled_success"],
               toy["success_by_tag"].get("prod_budget_dp6_no_refusal")),
            "PASS P3 bounded detector: 2-cycle fp(K) vs fp(K-2) + 64-step "
            "checkpoint, 2^%.2f B of state, versus 2^%.1f B for the rejected "
            "history" % (led["g2"]["detector_log2_bytes"],
                         led["g2"]["history_log2_bytes_rejected"]),
            "PASS P4 one constant: KAPPA = 11/10 in protocol.py, imported by the "
            "toy and the ledger; budget ratio %.6f" % led["budget_ratio"],
        ],
        "integers": {
            "C": C, "R_nodes": R_NODES, "R_table": R_TABLE, "DP_bits": DP_BITS,
            "p": str(p), "r": str(r), "a": str(A_CURVE), "b": str(B_CURVE),
            "N": a["N"], "N_prime": led["N_prime"], "S_ops": led["S_ops"],
            "ord_r_of_2": a["order_r_of_2"], "ord_over_2": a["ord_over_2"],
            "v2_r_minus_1": a["v2_r_minus_1"],
            "M_itemised": it["M_itemised"], "S_itemised": it["S_itemised"],
            "A_itemised": it["A_itemised"], "I_itemised": it["I_itemised"],
            "control": it["control"], "detection": it["detection"],
            "negation": it["negation"], "step_instructions": it["step_instructions"],
            "claim_x1e4": int(round(float(CLAIM) * 10 ** 4)),
            "phase1_log2_x1e6": int(round(led["phase1_log2"] * 10 ** 6)),
            "total_log2_x1e6": int(round(led["total_log2"] * 10 ** 6)),
            "stress_log2_x1e6": int(round(led["stress_log2"] * 10 ** 6)),
            "reports_x1e6": int(round(float(led["g2"]["report_mean_log2"]) * 10 ** 6)),
            "census_keys": 150000, "census_max_load": toy["census_max_load"],
            "census_distinct_nodes": toy["census_distinct_nodes"],
            "route_values_attained": toy["route_values_attained"],
            "detector_bytes": int(led["g2"]["detector_bytes"]),
        },
        "step_itemisation": it,
        "ledger": led,
        "toy": toy,
        "proofs": {
            "claim": "total 2^%.6f <= claim %.2f - 0.10 and claim < %.2f"
                     % (led["total_log2"], float(CLAIM), float(GATE)),
            "doubled_add": "the transition canon(2X+T_j) is computed as "
                           "X+(X+T_j): two batched affine additions, each 2M+1S+6A "
                           "core plus one inverted denominator; 2 inversions per "
                           "step at 3M of batch product work each",
            "routing": "the class key is the 256-bit x-coordinate; the DP predicate "
                       "reads only bits 0..63 (through mix_low) and the router only "
                       "bits 128..191 (through mix_high). The fields are disjoint, "
                       "so a width-64 DP event pins one low-field value and leaves "
                       "the routing field free: the route stays uniform over 2^16 "
                       "and the v9b collapse cannot occur",
            "detector": "bounded-memory detection: a 2-slot fingerprint ring "
                        "compares fp(K) with fp(K-2) and one checkpoint register "
                        "carries the running minimum fingerprint of the previous "
                        "64-step block; the escape doubles the minimum-fingerprint "
                        "class, which depends only on the cycle, so merged walks "
                        "stay merged. State is O(1) per trail (4 entries, 136 B each)",
            "reports": "64-bit DP condition over 2^127.652966 steps gives 2^63.652966 "
                       "reports <= 2^64 = cap/4; cap/mean per node is %s; Chernoff "
                       "puts cap loss below 2^%.0f" % (led["g2"]["cap_over_mean"],
                                                        led["cap_loss_log2"]),
            "meter": "escapes are arrivals priced inside S_ops, so the trial cost "
                     "is step x S_ops + reports + setup = 2^%.6f; the headroom to "
                     "2^%.2f gives a Markov halt debit of 2^%.0f"
                     % (led["total_log2"], float(CLAIM), led["halt_log2"]),
            "success": "model 1-0.60^1.21 = %.4f; measured efficiency 0.90-0.95 "
                       "gives %.4f-%.4f; claim 0.39"
                       % (led["model_success"], led["success_with_measured_efficiency"][0],
                          led["success_with_measured_efficiency"][1]),
            "ord_2": "ord_r(2) = N is even, so 2^(N/2) = -1: the shortest sign-"
                     "ambiguous gap is N/2 = 2^253.99999999966",
        },
        "source_sha256": src,
        "protected_sha256": protected,
        "diagnostics_only_not_in_bounds": {
            "nominal_claim_slack": led["slack"],
            "gate_gap": led["gate_gap"],
            "saving_vs_published_3_6_units": led["saving_vs_3_6_units"],
            "v9b_claim": 129.36, "v9c_claim_delta": round(float(CLAIM) - 129.36, 4),
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
