"""Exact certificate for proof.md: parameters, success bound, work ledger, claim, memory (integers and Fractions only).
Run python3 certificate.py; exit status 0 and a final PASS mean every check holds."""
import math
import sys
from fractions import Fraction

FAILED = []


def check(name, ok, detail=""):
    """Record one named check.

    name: description; ok: truth value; detail: optional text.
    returns: None (prints the line, remembers failures).
    """
    print("%-66s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
    if not ok:
        FAILED.append(name)


def lg(x):
    """Floating log2 of a positive rational, for display only.

    x: int or Fraction.
    returns: float.
    """
    f = Fraction(x)
    return math.log2(f.numerator) - math.log2(f.denominator)


def cdiv(a, b):
    """Ceiling of a/b for integers, b > 0.

    a, b: integers.
    returns: integer.
    """
    return -((-a) // b)


def exp_interval(x, terms=60):
    """Rational enclosure of exp(x) for |x| <= 1.

    x: Fraction; terms: series length.
    returns: (lo, hi) Fractions with lo <= exp(x) <= hi.
    """
    assert abs(x) <= 1
    s, t = Fraction(1), Fraction(1)
    for k in range(1, terms + 1):
        t = t * x / k
        s += t
    rem = abs(t * x) / (terms + 1) * 3
    return s - rem, s + rem


def ln2_interval(bits=300):
    """Rational enclosure of ln 2 = 2 atanh(1/3).

    bits: target width 2^-bits.
    returns: (lo, hi).
    """
    z, z2 = Fraction(1, 3), Fraction(1, 9)
    tot, pw, k = Fraction(0), z, 0
    while True:
        tot += pw / (2 * k + 1)
        rest = pw * z2 / ((2 * k + 3) * (1 - z2))
        if rest < Fraction(1, 2 ** bits):
            return 2 * tot, 2 * (tot + rest)
        pw *= z2
        k += 1


LN2 = ln2_interval()


def pow2_interval(c):
    """Rational enclosure of 2^c for a rational c.

    c: Fraction.
    returns: (lo, hi).
    """
    n = c.numerator // c.denominator
    f = c - n
    lo, hi = exp_interval(f * LN2[0])[0], exp_interval(f * LN2[1])[1]
    sc = Fraction(2) ** n
    return lo * sc, hi * sc


def ceil_lg(x, scale=100):
    """Least multiple of 1/scale that is at least log2(x), certified by the interval bounds of pow2_interval.

    x: positive int or Fraction; scale: denominator of the result.
    returns: Fraction.
    """
    c = Fraction(math.ceil(lg(x) * scale), scale)
    while pow2_interval(c)[0] < x:
        c += Fraction(1, scale)
    while pow2_interval(c - Fraction(1, scale))[1] >= x:
        c -= Fraction(1, scale)
    return c


# ---------------------------------------------------------------- 1. constants
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
C0 = 2 ** 32 + 977
C = 1000                      # instructions per target-group-operation equivalent
NODES = 2 ** 72               # 2^36 x 2^36 array
SIDE = 2 ** 36
Q2 = 2 ** 27                  # Y generators per node
D1, D2, D3, D4 = 36, 36, 24, 48
KEYBITS = 128                 # stored key words A = K2, B = K3
THRESH = Fraction(39001, 100000)   # target with slack above the rules' 0.39

# executed counts (experiment.py parts A and B re-derive each of these by running the kernels)
FWD_INC, FWD_FIX = 105, 3         # leaf.s fwd per incidence; mov, movk, ret
BWD_INC, BWD_FIX = 523, 5         # leaf.s bwd; mov, movk, x30 save and restore, ret
KEY_INC, KEY_FIX = 66, 11         # leaf.s keys per incidence (two leaves) incl. the two emits; set-up, x30 save/restore, ret
EMIT_INC = 13                     # the two emits (6 + 7): the sender half of hop 1
KMUL_PATH = 83                    # extra on a leaf whose class index j is 1 or 2 (KMUL 81, counter, branch)
INV = 23646                       # secp_fe.s fe_inv, executed
EVENT_PATH = 3                    # LDR of the event, ADD of the window address, STR clearing the flag
H1RECV, H1RECV_DROP = 8 + EVENT_PATH, 9   # route.s h1recv per event: stored, or dropped at a full buffer
H2SEND = 9                        # route.s h2send per record
H2RECV, H2RECV_DROP = 17 + EVENT_PATH, 12  # route.s h2recv per event: stored, or dropped at a full d3 bucket
PASS4_REC, PASS4_BKT, SCAN_REC = 9, 7, 7
ZERO4_PER4, PREFIX4_PER4, INIT3_BKT, WCLR_PER4 = 4, 14, 5, 7
HOP_PRICE = 40                    # per record hop: the executed protocol (17.5 or 29) plus a margin
ITEM_PRICE = 20                   # per record in the local pass (executed 9)
HOP1_EXEC, HOP2_EXEC = Fraction(EMIT_INC, 2) + H1RECV, H2SEND + H2RECV
check("key pass: 66 = 2*(2 LDP + KEY 21 + B.LO + CBZ) + emits 13 + loop 3",
      KEY_INC == 2 * (2 + 21 + 1 + 1) + EMIT_INC + 3)
check("hop 1 executed 17.5 and hop 2 executed 29 (event path included) are below 40",
      HOP1_EXEC == Fraction(35, 2) and HOP2_EXEC == 29 and max(HOP1_EXEC, HOP2_EXEC) < HOP_PRICE)
DISPATCH, HW_IN_SOFTWARE = 4, 4   # handler loop per record; event append (2) and uncleared-window counters (2)
DISPATCH_PER_CALL = 7             # alternative: the kernel called once per event (wait, call, CBZ, RET, wrap)
check("with dispatch, the event append and the window counters charged in software both hops stay below 40",
      HOP1_EXEC + DISPATCH + HW_IN_SOFTWARE < HOP_PRICE and HOP2_EXEC + DISPATCH + HW_IN_SOFTWARE < HOP_PRICE,
      "%s and %s" % (HOP1_EXEC + DISPATCH + HW_IN_SOFTWARE, HOP2_EXEC + DISPATCH + HW_IN_SOFTWARE))
check("with the kernel called once per event (at most 7 per record of dispatch) both hops stay below 40",
      HOP1_EXEC + DISPATCH_PER_CALL < HOP_PRICE and HOP2_EXEC + DISPATCH_PER_CALL < HOP_PRICE,
      "%s and %s" % (HOP1_EXEC + DISPATCH_PER_CALL, HOP2_EXEC + DISPATCH_PER_CALL))
check("with per-call dispatch and the hardware work charged in software both hops are at most 40 (hop 2: 29 + 7 + 4)",
      HOP1_EXEC + DISPATCH_PER_CALL + HW_IN_SOFTWARE <= HOP_PRICE
      and HOP2_EXEC + DISPATCH_PER_CALL + HW_IN_SOFTWARE <= HOP_PRICE,
      "%s and %s" % (HOP1_EXEC + DISPATCH_PER_CALL + HW_IN_SOFTWARE, HOP2_EXEC + DISPATCH_PER_CALL + HW_IN_SOFTWARE))
check("local pass executed 9 per record is below its price 20", PASS4_REC < ITEM_PRICE)

# priced allowances (analytic, per item)
GEN_MUL = 2 ** 21          # one complete scalar multiplication [k]P or [k]Q with conversion and table write
COIN_SCALAR = 64           # at most two draws of 4 RNDR words (18 with one charged retry each) plus a range test (7): 50
BATCH_SETUP = 64           # copy X into H, x_X registers, acc = 1, pointers, index word, meter test, loop
INV_CALL = 16              # store acc, pointer set-up, BL, load R
NODE_INIT = 2 ** 20        # program <= 2^16 bytes = 2^11 records over 2 hops at 40 (2^17.3), plus port, buffer and pointer state
HEADER = 148               # G1..G5 and the two row tables (37 words at 4 per word)
CAND_CAP = 2 ** 20         # candidate pairs resolved per home node
RESOLVE = 2 ** 25          # per candidate: two requests and two replies of <= 8 records over 2 hops at 40, two exact
                           # additions, unit test, solve, [x]P = Q (554,530 itemized)
REQ_CAP = 2 ** 21          # requests answered per source node
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 22        # relayed records held per node in each phase-4 sub-phase
WINDOWS = 2 ** 37          # one window per node of the row and of the column, its own included (loopback)
BARRIERS, BARRIER = 16, 2 ** 82
OUTPUT = 2 ** 20
DELTA_J = 2 ** 42          # meter slack for the class-index products

# ---------------------------------------------------------------- 2. success
P1 = Fraction(6 * (r - 1) + 1, r * r)    # P(B in class(A)) for independent uniform A, B
TWO_BOUND = Fraction(2, r - 1)           # x-coordinate of a point uniform on at least r - 1 points: each value <= 2/(r-1)

# batch arithmetic sites (proof.md 3.4). Only row shifts (SHX), folds (FDX, also inside KMUL) and the fused-square folds
# (SQP, SQM) drop a carry, and only they can return a value >= p. Given X != O and every other generator, the fresh Y_l
# puts each operand value of a site at probability at most TWO_BOUND, so a site fails with probability <= BAD * TWO_BOUND.
SHX_BAD = 2 ** 64 * 2 ** 33     # rows B with B mod 2^192 > 2^192 - 2^33: needed for a row-shift carry
FOLD_BAD = 2 ** 66 * C0         # values in [C0, 2^66 C0): needed for the carry out of a fold of a sum below 2^322
SQX_BAD = C0 * C0 + 3 * C0      # x(X +- Y_l) below C0^2 + 3 C0: needed for a fused-square carry or a result >= p
KMUL_BAD = 2 * 2 ** 66 * C0     # x(X +- Y_l) with beta^j x below 2^66 C0 for j = 1 or 2: fold carry or K >= p
SITES = [  # (site, executions per batch of Q2 incidences, bad operand values)
    ("fwd: row shifts of acc_(l-1), l >= 1", 3 * (Q2 - 1), SHX_BAD),
    ("fwd: fold of acc_l, l >= 1", Q2 - 1, FOLD_BAD),
    ("bwd: row shifts of R_l", 3 * Q2, SHX_BAD),
    ("bwd: fold of inv_l = R_l acc_(l-1), l >= 1", Q2 - 1, FOLD_BAD),
    ("bwd: fold of R_(l-1) = R_l dx_l, l >= 1", Q2 - 1, FOLD_BAD),
    ("bwd: row shifts of inv_l", 3 * Q2, SHX_BAD),
    ("bwd: folds of lambda+ and lambda-", 2 * Q2, FOLD_BAD),
    ("bwd: fused squares x+ and x-", 2 * Q2, SQX_BAD),
    ("keys: KMUL folds of the two leaves", 2 * Q2, KMUL_BAD),
]
SITES_PER_INCIDENCE = dict(SHX=9, FDX=7, SQX=2)   # per pass of the fwd, bwd and keys loops (experiment.py counts them)
EXACT_AT_L0 = dict(SHX=3, FDX=3)   # rows of acc_(-1) = 1, folds of dx_0 * 1, R_0 * 1 and R_0 dx_0 = 1 (never carry)
D_SITES = sum(n * bad for _, n, bad in SITES) * TWO_BOUND
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed), r - 1)  # P(fixed key bits of a leaf take given values)


def success_terms(q1):
    """Success lower bound pieces for q1 X generators per node.

    q1: integer.
    returns: dict of exact Fractions.
    """
    n0 = 2 * q1 * Q2
    S = NODES * n0
    gamma = Fraction(S * S - NODES * n0 * n0, 2)
    mu = gamma * P1
    deg = 4 * (q1 + Q2) * S
    cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
    mu_lo = Fraction((mu.numerator << 200) // mu.denominator, 2 ** 200)
    e_hi = exp_interval(-mu_lo, 40)[1]
    # per-leaf / per-pair defect probabilities
    d_batch = D_SITES + Q2 * TWO_BOUND + Fraction(1, r)   # sites, dx_l = 0 (Y_l = +-X), X = O
    d_win = TWO_BOUND * 13 * 2 ** 197
    d_wj = TWO_BOUND * 2 ** 64 * C0
    var_j = Fraction(2, 9) * n0 + 1
    d_j = var_j / Fraction(DELTA_J - n0 // q1 - 1) ** 2 + q1 * d_batch
    m_c = SIDE * n0 * KB(D1) + 1
    cap_c = 6 * n0 + 2 ** 50
    d_c = (m_c + 1) / Fraction(cap_c - 2 - m_c) ** 2
    m_3 = S * KB(D1 + D2 + D3) + 1
    cap_3 = cdiv(m_3.numerator, m_3.denominator) + 2 ** 32
    d_3 = (m_3 + 1) / Fraction(cap_3 - 2 - m_3) ** 2
    lam_int = S * KB(D1 + D2 + D3 + D4)
    other_h = Fraction(S) * S * KB(D1 + D2) * KB(KEYBITS + D1 + D2) + 2 * S * KB(KEYBITS + D1 + D2)
    d_cl = other_h / (CAND_CAP - 4)
    other_g = Fraction(n0) * S * KB(KEYBITS + D1 + D2) * 2
    d_rq = other_g / (REQ_CAP - 8)
    other_v = (2 + 2 * 8) * Fraction(S, SIDE) * S * KB(D2) * KB(KEYBITS + D1 + D2)  # requests 1 record, replies <= 8
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_win + d_wj + d_j + d_c + d_3 + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)  # RNDR is the rules' ideal coin source: its failure branch is charged, never taken
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + NODES * Q2 * Fraction(2 * C0, r - 1) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_win=d_win, d_wj=d_wj, d_j=d_j,
                d_c=d_c, d_3=d_3, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def, rows=rows, glob=glob,
                lo=lo, cap_c=cap_c,
                cap_3=cap_3, m_c=m_c, m_3=m_3)


lo_q, hi_q = 2 ** 26, 2 ** 27
while hi_q - lo_q > 1:
    mid = (lo_q + hi_q) // 2
    if success_terms(mid)["lo"] >= THRESH:
        hi_q = mid
    else:
        lo_q = mid
Q1 = hi_q
T = success_terms(Q1)
N0, S = T["n0"], T["S"]
check("q1 = %d is the least with success bound >= 0.39001" % Q1, T["lo"] >= THRESH and success_terms(Q1 - 1)["lo"] < THRESH)
check("leaves S = 2^73*q1*q2 = 2^%.6f" % lg(S), S == 2 ** 73 * Q1 * Q2)
check("mu = |Gamma|*P1 = %.9f (cross-node pairs)" % float(T["mu"]), T["mu"] > Fraction(494, 1000))
check("Chen-Stein b1+b2 < 2^-96", T["cs"] < Fraction(1, 2 ** 96), "2^%.2f" % lg(T["cs"]))
check("intruder rate per colliding pair < 2^-14.6", T["lam_int"] < pow2_interval(Fraction(-146, 10))[0], "2^%.3f" % lg(T["lam_int"]))
check("window-drop rate per leaf < 2^-54.29", T["d_win"] < pow2_interval(Fraction(-5429, 100))[0], "2^%.4f" % lg(T["d_win"]))
for name, n, bad in SITES:
    k = -(-n // Q2)
    print("  site %-46s %-6s bad <= 2^%.2f, each <= 2^%.2f" % (name, ("%dq2" % k if k > 1 else "q2") + (
        "%+d" % (n - k * Q2) if n != k * Q2 else ""), ceil_lg(bad), ceil_lg(bad * TWO_BOUND)))
check("site table: 9 row shifts, 7 folds, 2 fused squares per incidence, less 3 + 3 exact at l = 0",
      sum(n for _, n, b in SITES if b == SHX_BAD) == SITES_PER_INCIDENCE["SHX"] * Q2 - EXACT_AT_L0["SHX"]
      and sum(n for _, n, b in SITES if b in (FOLD_BAD, KMUL_BAD)) == SITES_PER_INCIDENCE["FDX"] * Q2 - EXACT_AT_L0["FDX"]
      and sum(n for _, n, b in SITES if b == SQX_BAD) == SITES_PER_INCIDENCE["SQX"] * Q2)
check("bad sets: rows 2^97, folds below 2^98.01, fused squares below 2^64.01, KMUL below 2^99.01",
      SHX_BAD == 2 ** 97 and all(b < pow2_interval(Fraction(e, 100))[0] for b, e in ((FOLD_BAD, 9801), (SQX_BAD, 6401),
                                                                                     (KMUL_BAD, 9901))))
check("batch defect rate (site table, dx = 0, X = O) < 2^-125", T["d_batch"] < Fraction(1, 2 ** 125),
      "2^%.3f" % lg(T["d_batch"]))
check("class-index meter rate < 2^-31", T["d_j"] < Fraction(1, 2 ** 31), "2^%.2f" % lg(T["d_j"]))
check("hop-1 buffer overflow rate < 2^-42", T["d_c"] < Fraction(1, 2 ** 42), "2^%.2f" % lg(T["d_c"]))
check("d3 bucket overflow rate < 2^-30", T["d_3"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_3"]))
check("candidate cap rate < 2^-33.4, request cap rate < 2^-35.9",
      T["d_cl"] < pow2_interval(Fraction(-334, 10))[0] and T["d_rq"] < pow2_interval(Fraction(-359, 10))[0],
      "2^%.3f 2^%.3f" % (lg(T["d_cl"]), lg(T["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-31.2", T["d_rl"] < pow2_interval(Fraction(-312, 10))[0],
      "2^%.3f" % lg(T["d_rl"]))
check("whole-run row: a computational defect in any batch < 2^-26", T["g_def"] < Fraction(1, 2 ** 26), "2^%.3f" % lg(T["g_def"]))
check("whole-run row: RNDR failure is 0 (RNDR is the rules' ideal coin source; its retry is charged, never taken)",
      T["glob"] >= 0, "0")
check("whole-run rows (defective batch anywhere, coin rejection twice, x_Y >= p - C0) < 2^-26",
      T["glob"] < Fraction(1, 2 ** 26), "2^%.3f" % lg(T["glob"]))
check("success >= 1 - e^-mu - (b1+b2) - rows - global >= 0.39001", T["lo"] >= THRESH,
      "%.10f" % (math.floor(T["lo"] * 10 ** 10) / 10 ** 10))

# ---------------------------------------------------------------- 3. work ledger (instructions)
CAP_J = cdiv(2 * N0, 3) + DELTA_J
GENS = Q1 + Q2 + 1
rows = []


def row(name, value, pre=False):
    """Add one ledger row.

    name: label; value: exact integer instructions; pre: True if the row is preprocessing (independent of Q).
    returns: None.
    """
    assert isinstance(value, int) and value >= 0
    rows.append((name, value, pre))


row("program delivery, node and port state: 2^72 * 2^20", NODES * NODE_INIT, True)
row("coins: (q1+q2+1) scalars * 64 per node", NODES * GENS * COIN_SCALAR, True)
row("Y generators [c]P: q2 * 2^21 per node", NODES * Q2 * GEN_MUL, True)
row("R_g = [b]Q and X = R_g + [a]P: (q1+1) * 2^21 per node", NODES * (Q1 + 1) * GEN_MUL)
row("header constants: 148 per node", NODES * HEADER, True)
row("window flags cleared: 2^37/4 * 7 + 1 per node", NODES * (WINDOWS // 4 * WCLR_PER4 + 1), True)
row("Q broadcast: 2 records to every node over two hops at 40", 2 * (SIDE + NODES) * HOP_PRICE)
row("d4 counters zeroed: 2^48/4 * 4 per node", NODES * (2 ** D4 // 4 * ZERO4_PER4 + 1), True)
row("d3 cursor table: 5 * 2^24 + 1 per node", NODES * (INIT3_BKT * 2 ** D3 + 1), True)
row("forward: q1 * (105 q2 + 3)", NODES * Q1 * (FWD_INC * Q2 + FWD_FIX))
row("inversion: q1 * (23,646 + 16)", NODES * Q1 * (INV + INV_CALL))
row("backward: q1 * (523 q2 + 5)", NODES * Q1 * (BWD_INC * Q2 + BWD_FIX))
row("key pass without emits: q1 * (53 q2 + 11)", NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX))
row("class-index products: cap_j * 83 per node", NODES * CAP_J * KMUL_PATH)
row("batch set-up: q1 * 64", NODES * Q1 * BATCH_SETUP)
row("hop 1 (emit + receive): S * 40", S * HOP_PRICE)
row("hop 2 (send + receive into d3 buckets): S * 40", S * HOP_PRICE)
row("prefix sums of d4: 2^48/4 * 14 per node", NODES * (2 ** D4 // 4 * PREFIX4_PER4 + 1))
row("local pass: S * 20 + (7 * 2^24 + 1) per node", S * ITEM_PRICE + NODES * (PASS4_BKT * 2 ** D3 + 1))
row("scan: S * 7 + (4 + 4 * (2^20 + 1)) per node", S * SCAN_REC + NODES * (4 + 4 * (CAND_CAP + 1)))
row("candidates: 2^20 * 2^25 per node", NODES * CAND_CAP * RESOLVE)
row("requests served: 2^21 * 2^10 per node", NODES * REQ_CAP * REQ_SERVE)
row("barriers: 16 * 2^82", BARRIERS * BARRIER)
row("output", OUTPUT)
R_TOTAL = sum(v for _, v, _ in rows)
R_PRE = sum(v for _, v, pre in rows if pre)
W = Fraction(R_TOTAL, C)
print()
for name, v, pre in rows:
    print("  %-58s %s 2^%.4f" % (name, "P" if pre else " ", lg(v)))
print("  %-58s   2^%.6f instructions" % ("total", lg(R_TOTAL)))
print()


def claim_of(w):
    """Smallest multiple of 10^-6 whose power of two is at least w, verified by interval bounds.

    w: Fraction > 1.
    returns: Fraction claim.
    """
    c = Fraction(math.ceil(lg(w) * 10 ** 6), 10 ** 6)
    for _ in range(5):
        if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 6))[1] < w:
            return c
        c += Fraction(1, 10 ** 6) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 6)
    raise AssertionError("claim not certified")


CLAIM = claim_of(W)
check("claim = least 6-decimal c with 2^c >= R/1000", pow2_interval(CLAIM)[0] >= W
      and pow2_interval(CLAIM - Fraction(1, 10 ** 6))[1] < W, "%s" % format(float(CLAIM), ".6f"))
PER_LEAF = Fraction(R_TOTAL, S)
check("priced work per leaf = %.4f instructions" % float(PER_LEAF), 495 < PER_LEAF < 510)

# executed-count variant: hops, local pass at their executed counts
R_EXEC = R_TOTAL - 2 * S * HOP_PRICE - S * ITEM_PRICE + S * HOP1_EXEC + S * HOP2_EXEC + S * PASS4_REC
CLAIM_EXEC = claim_of(Fraction(R_EXEC, C))
check("claim at executed counts", True, "%.6f (%.4f per leaf)" % (float(CLAIM_EXEC), float(R_EXEC / S)))

# ---------------------------------------------------------------- 4. memory (bytes per node, explicit layout)
CAP_C, CAP_3 = T["cap_c"], T["cap_3"]
LAYOUT = [("code, constants, stack, port and flag state", 2 ** 20), ("header", 512), ("X table: x, y", Q1 * 64),
          ("Y table: x, y, x + 2C0", Q2 * 96), ("scalars a, b, c", GENS * 32), ("batch slots", Q2 * 64),
          ("hop-1 receive buffer (cap_c records)", CAP_C * 32), ("d3 buckets (2^24 * cap_3 records)", 2 ** D3 * CAP_3 * 32),
          ("d3 cursor table", 2 ** D3 * 16), ("d4 counters", 2 ** D4 * 8), ("pointer array", 2 ** D3 * CAP_3 * 8),
          ("candidate list", CAND_CAP * 16), ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64),
          ("relay buffer (2^22 entries)", RELAY_CAP * 64), ("windows: 2^37 x (32-byte slot, flag)", WINDOWS * 64),
          ("event queue (2^37 entries)", WINDOWS * 8), ("row and column ports (device registers)", 2 * 64)]
NODE_BYTES = sum(b for _, b in LAYOUT)
for name, b in LAYOUT:
    print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("every node's regions fit one 64-bit address space", NODE_BYTES < 2 ** 64, "2^%.4f" % lg(NODE_BYTES))
PEAK = NODES * NODE_BYTES
check("peak memory 2^72 * node bytes = 2^%.4f" % lg(PEAK), PEAK < 2 ** 136)
check("caps exceed crude means: cap_c > m_c, cap_3 > m_3", CAP_C > T["m_c"] + 2 ** 49 and CAP_3 > T["m_3"] + 2 ** 31)
check("preprocessing (rows marked P) = 2^%.4f units" % lg(Fraction(R_PRE, C)), R_PRE < R_TOTAL)
ANALYTIC = sum(v for n, v, _ in rows if n.split(":")[0] in ("program delivery, node and port state", "coins",
               "Y generators [c]P", "R_g = [b]Q and X = R_g + [a]P", "header constants", "batch set-up", "candidates",
               "requests served", "barriers", "output", "inversion")) - NODES * Q1 * INV
# itemized work behind the allowances (proof.md section 7)
# on the exact arithmetic of proof.md 3.4: fe_mul_f and fe_sqr_f as executed (115 and 82) plus their calls (pointers and
# BL: 4 and 3, as in fe_inv), and FS or FADDW (9) with both operands brought below p (2 x 9: ADDS chain of C0, CSELs)
MULX, SQRX, ADDX = 115 + 4, 82 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
check("generator (exact arithmetic): doubling 933, mixed addition 1,477, total < 2^19.4 (2^21 is > 3x)",
      DBL == 933 and MADD_ == 1477 and GEN_ITEM < pow2_interval(Fraction(194, 10))[0] and GEN_MUL > 3 * GEN_ITEM,
      "%d" % GEN_ITEM)
BARRIER_ITEM = NODES * (4 * HOP_PRICE + 2 * H1RECV)
check("barrier: 4 messages per node at 40 + 2 leader handler steps at 11 < 2^79.6 (2^82 is > 5.5x)",
      BARRIER_ITEM < pow2_interval(Fraction(796, 10))[0] and BARRIER > Fraction(55, 10) * BARRIER_ITEM,
      "2^%.3f" % lg(BARRIER_ITEM))
CAND_ITEM = 4 * 8 * 2 * HOP_PRICE + 2 * (INV + 600) + 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, two additions, unit test, solve mod r, [x]P = Q < 2^20 (2^25 is > 30x)",
      CAND_ITEM < 2 ** 20 and RESOLVE > 30 * CAND_ITEM, "%d" % CAND_ITEM)
PROGRAM_BYTES = 2 ** 16
check("program delivery: 2^16 bytes = 2^11 records over 2 hops at 40 = 2^17.3 < node set-up 2^20",
      PROGRAM_BYTES // 32 * 2 * HOP_PRICE < NODE_INIT)
check("analytic allowances are below 2^-14.6 of the total", ANALYTIC * pow2_interval(Fraction(146, 10))[1] < R_TOTAL,
      "2^%.3f" % lg(Fraction(ANALYTIC, R_TOTAL)))

print()
print("q1 %d  q2 2^27  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1, N0, lg(S), float(T["mu"]),
                                                                    math.floor(T["lo"] * 10 ** 10) / 10 ** 10))
print("R_total %d instructions  W 2^%.6f  claim %s  executed-count claim %s" % (
    R_TOTAL, lg(W), format(float(CLAIM), ".6f"), format(float(CLAIM_EXEC), ".6f")))
print("memory_log2_bytes %.4f  preprocessing_log2 %.4f" % (lg(PEAK), lg(Fraction(R_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
    sys.exit(1 if FAILED else 0)
