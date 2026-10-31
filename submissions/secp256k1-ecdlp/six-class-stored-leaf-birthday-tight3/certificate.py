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
THRESH = Fraction(3900001, 10000000)   # target, 10^-7 above the rules' 0.39

# executed counts (experiment.py parts A and B re-derive each of these by running the kernels)
FWD_INC, FWD_FIX = 105, 3         # leaf.s fwd per incidence; mov, movk, ret
BWD_INC, BWD_FIX = 514, 5         # leaf.s bwd; mov, movk, x30 save and restore, ret
EMIT = (12, 13)                   # keys' emit of the + and - leaf: the whole of hop 1's work for a record
EMIT_CHECK = 5                    # d1, cursor entry address, LDP of (cursor, slot end), CMP, B.HS: run before a drop
EMIT_INC = sum(EMIT)              # 25 per incidence
KEY_INC, KEY_FIX = 78, 11         # leaf.s keys per incidence (two leaves) incl. the two emits; set-up, x30 save/restore, ret
KMUL_PATH = 83                    # extra on a leaf whose class index j is 1 or 2 (KMUL 81, counter, branch)
INV = 23646                       # secp_fe.s fe_inv, executed
TINIT_ENTRY = 3                   # route.s tinit per cursor-table entry
FIN_DEST = 7                      # route.s fin per destination: LDR, three STPs, ADD, SUBS, B.NE
H2S_REC, H2S_DROP, H2S_SLOT = 12, 7, 11   # route.s h2s per record sent, per record dropped at a full slot, per slot
BKT_REC, BKT_DROP, BKT_SLOT = 11, 6, 12   # route.s bkt per record bucketed (STADD counts d4), dropped at a full bucket, per slot
PASS4_REC, PASS4_BKT, SCAN_REC = 4, 11, 7
LOOP8 = 2                                 # h2s, bkt, pass4: CMP and B.LO per 8 records of a slot's or bucket's main loop
ZERO4_PER4, PREFIX4_PER4, INIT3_BKT = 4, 10, 5   # prefix4: two LDPs, two STPs and four ADDs per block of four
check("key pass: 78 = 2*(2 LDP + KEY 21 + B.LO + CBZ) + emits 12 + 13 + loop 3",
      KEY_INC == 2 * (2 + 21 + 1 + 1) + EMIT_INC + 3)
check("emit: 5 up to the slot test, then d2, w3 (and the sign bit), three STPs, cursor ADD and STR",
      EMIT == (EMIT_CHECK + 7, EMIT_CHECK + 8))
DONE = 2   # per message: the port's in-flight count, up when the link takes it and down when it is written
DLV = 6    # per message: its delivery, route.s dlv (three LDPs of the message, two STPs at its address, RET)
MSG_LINK = DONE + DLV   # link-side work of one message at one hop, beside the sender's own instructions
L8 = Fraction(LOOP8, 8)   # main-loop control per record, an upper bound: a slot of n records runs floor(n/8) passes
HOP1, HOP2 = Fraction(EMIT_INC, 2) + MSG_LINK, H2S_REC + L8 + MSG_LINK
GROUP = BKT_REC + L8 + PASS4_REC + L8 + SCAN_REC
check("per record: hop 1 20.5, hop 2 20.25 (each with in-flight count 2 and delivery 6), grouping 22.5; 63.25 in all",
      HOP1 + HOP2 + GROUP == Fraction(253, 4))

# priced allowances (analytic, per item)
GEN_MUL = 2 ** 20          # one complete scalar multiplication [k]P or [k]Q with conversion and table write;
                           # the itemized total (665,910 arithmetic + 2,048 control) fits with 1.5699x headroom
COIN_SCALAR = 64           # at most two draws of 4 RNDR words (18 with one charged retry each) plus a range test (7): 50
BATCH_SETUP = 64           # copy X into H, x_X registers, acc = 1, pointers, index word, meter test, loop
INV_CALL = 16              # store acc, pointer set-up, BL, load R
NODE_INIT = 2 ** 20        # program <= 2^16 bytes = 2^11 messages over 2 hops at 40 (2^17.3), plus node state
HEADER = 148               # G1..G5 and the two row tables (37 words at 4 per word)
MSG = 40                   # a message outside the two hops (program, Q, phase 4), per hop, with its in-flight count
                           # and delivery (a broadcast message: 16 in all)
CAND_CAP = 2 ** 20         # candidate pairs resolved per home node
RESOLVE = 2 ** 21          # per candidate: two requests and two replies of <= 8 records over 2 hops at 40, two exact
                           # additions, unit test, solve, [x]P = Q (769,058 itemized, 2.727x headroom)
REQ_CAP = 2 ** 21          # requests answered per source node
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 22        # relayed records held per node in each phase-4 sub-phase
EXCH = 2 ** 41             # one phase-4 count exchange per node (4 per run)
BARRIERS, BARRIER = 32, 2 ** 82   # 2^10 per node and barrier
OUTPUT = 2 ** 20
DELTA_J = 2 ** 42          # meter slack for the class-index products
CAP_P = 2 ** 22 - 2 ** 12  # records per pair slot, at hop 1 and at hop 2

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
    m_p = n0 * KB(D1) + 1        # a node's leaves on d1 (hop 1), or a column's leaves on (d1, d2) (hop 2): equal
    d_p = (m_p + 1) / Fraction(CAP_P - 2 - m_p) ** 2
    m_3 = S * KB(D1 + D2 + D3) + 1
    cap_3 = cdiv(m_3.numerator, m_3.denominator) + 2 ** 32
    d_3 = (m_3 + 1) / Fraction(cap_3 - 2 - m_3) ** 2
    lam_int = S * KB(D1 + D2 + D3 + D4)
    # Exceptional pair types (proof.md 6.4): pairs of leaves of one node that share a generator,
    # which the independent-pair products above do not cover. Same-X pairs are S*q2, same-Y pairs
    # S*q1 and the +- partners S/2. Fixing the first leaf (and the generators) leaves the second
    # uniform on at least r - 2 points for the same-X and same-Y types; for a +- partner the map
    # Y -> x(X - Y) hits a fixed 200-bit key pattern at most 6*2^56 times. So each such pair matches
    # on the 200 compared key bits with probability at most 2*KB(200): a mass far below the
    # independent products.
    dep_pairs = S * Q2 + S * q1 + S // 2
    dep_pair_mass = dep_pairs * 2 * KB(KEYBITS + D1 + D2)
    other_h = (Fraction(S) * S * KB(D1 + D2) * KB(KEYBITS + D1 + D2) + 2 * S * KB(KEYBITS + D1 + D2)
               + dep_pair_mass)
    d_cl = other_h / (CAND_CAP - 4)
    other_g = Fraction(n0) * S * KB(KEYBITS + D1 + D2) * 2 + 2 * dep_pair_mass
    d_rq = other_g / (REQ_CAP - 8)
    other_v = ((2 + 2 * 8) * Fraction(S, SIDE) * S * KB(D2) * KB(KEYBITS + D1 + D2)  # requests 1 record, replies <= 8
               + (2 + 2 * 8) * dep_pair_mass)
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_win + d_wj + d_j + 2 * d_p + d_3 + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)  # RNDR is the rules' ideal coin source: its failure branch is charged, never taken
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + NODES * Q2 * Fraction(2 * C0, r - 1) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_win=d_win, d_wj=d_wj, d_j=d_j,
                d_p=d_p, d_3=d_3, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, dep=dep_pair_mass, other_h=other_h, other_g=other_g, other_v=other_v, g_def=g_def, rows=rows, glob=glob,
                lo=lo, cap_3=cap_3, m_p=m_p, m_3=m_3)


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
check("q1 = %d is the least with success bound >= 0.3900001" % Q1, T["lo"] >= THRESH and success_terms(Q1 - 1)["lo"] < THRESH)
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
check("pair-slot overflow rate, each hop, < 2^-19.9", T["d_p"] < pow2_interval(Fraction(-199, 10))[0],
      "2^%.3f" % lg(T["d_p"]))
check("d3 bucket overflow rate < 2^-30", T["d_3"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_3"]))
check("candidate cap rate < 2^-33.4, request cap rate < 2^-35.9",
      T["d_cl"] < pow2_interval(Fraction(-334, 10))[0] and T["d_rq"] < pow2_interval(Fraction(-359, 10))[0],
      "2^%.3f 2^%.3f" % (lg(T["d_cl"]), lg(T["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-31.2", T["d_rl"] < pow2_interval(Fraction(-312, 10))[0],
      "2^%.3f" % lg(T["d_rl"]))
check("exceptional same-X, same-Y and +- pairs add 2^%.2f, over 2^25 below the independent terms" % lg(T["dep"]),
      T["dep"] < pow2_interval(Fraction(-40, 1))[0] and T["dep"] * 2 ** 25 < T["other_h"]
      and 2 * T["dep"] * 2 ** 25 < T["other_g"] and 18 * T["dep"] * 2 ** 25 < T["other_v"])
check("whole-run row: a computational defect in any batch < 2^-26", T["g_def"] < Fraction(1, 2 ** 26), "2^%.3f" % lg(T["g_def"]))
check("whole-run row: RNDR failure is 0 (RNDR is the rules' ideal coin source; its retry is charged, never taken)",
      T["glob"] >= 0, "0")
check("whole-run rows (defective batch anywhere, coin rejection twice, x_Y >= p - C0) < 2^-26",
      T["glob"] < Fraction(1, 2 ** 26), "2^%.3f" % lg(T["glob"]))
check("success >= 1 - e^-mu - (b1+b2) - rows - global >= 0.3900001", T["lo"] >= THRESH,
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
row("Y generators [c]P: q2 * 2^20 per node", NODES * Q2 * GEN_MUL, True)
row("R_g = [b]Q and X = R_g + [a]P: (q1+1) * 2^20 per node", NODES * (Q1 + 1) * GEN_MUL)
row("header constants: 148 per node", NODES * HEADER, True)
row("cursor tables: 2 * (3 * 2^36 + 1) per node", NODES * 2 * (TINIT_ENTRY * SIDE + 1), True)
row("Q broadcast: 2 messages to every node over two hops at 40", 2 * (SIDE + NODES) * MSG)
row("d4 counters zeroed: 2^48/4 * 4 + 1 per node", NODES * (2 ** D4 // 4 * ZERO4_PER4 + 1), True)
row("d3 cursor table: 5 * 2^24 + 1 per node", NODES * (INIT3_BKT * 2 ** D3 + 1), True)
row("forward: q1 * (105 q2 + 3)", NODES * Q1 * (FWD_INC * Q2 + FWD_FIX))
row("inversion: q1 * (23,646 + 16)", NODES * Q1 * (INV + INV_CALL))
row("backward: q1 * (514 q2 + 5)", NODES * Q1 * (BWD_INC * Q2 + BWD_FIX))
row("key pass without emits: q1 * (53 q2 + 11)", NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX))
row("class-index products: cap_j * 83 per node", NODES * CAP_J * KMUL_PATH)
row("batch set-up: q1 * 64", NODES * Q1 * BATCH_SETUP)
row("hop 1, the emits, in-flight counts and deliveries: q1 * 41 q2", NODES * Q1 * (EMIT_INC + 2 * MSG_LINK) * Q2)
row("hop 2, h2s, in-flight counts and deliveries: S * 20.25 + 11 * 2^36 + 1 per node",
    S * (H2S_REC + MSG_LINK) + S * LOOP8 // 8 + NODES * (H2S_SLOT * SIDE + 1))
row("final cursors, in-flight counts and deliveries, both hops: 2 * (15 * 2^36 + 3) per node",
    NODES * 2 * ((FIN_DEST + MSG_LINK) * SIDE + 3))
row("bucket pass bkt: S * 11.25 + 12 * 2^36 + 2 per node", S * BKT_REC + S * LOOP8 // 8 + NODES * (BKT_SLOT * SIDE + 2))
row("prefix sums of d4: 2^48/4 * 10 + 1 per node", NODES * (2 ** D4 // 4 * PREFIX4_PER4 + 1))
row("local pass pass4: S * 4.25 + 11 * 2^24 + 1 per node", S * PASS4_REC + S * LOOP8 // 8 + NODES * (PASS4_BKT * 2 ** D3 + 1))
row("scan: S * 7 + 4 + 4 * (2^20 + 1) per node", S * SCAN_REC + NODES * (4 + 4 * (CAND_CAP + 1)))
row("candidates: 2^20 * 2^21 per node", NODES * CAND_CAP * RESOLVE)
row("requests served: 2^21 * 2^10 per node", NODES * REQ_CAP * REQ_SERVE)
row("phase-4 count exchanges: 4 * 2^41 per node", NODES * 4 * EXCH)
row("barriers: 32 * 2^82", BARRIERS * BARRIER)
row("phase control: kernel calls between barriers, 2^10 per node", NODES * 2 ** 10)
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
    """Smallest multiple of 10^-7 whose power of two is at least w, verified by interval bounds.

    w: Fraction > 1.
    returns: Fraction claim.
    """
    c = Fraction(math.ceil(lg(w) * 10 ** 7), 10 ** 7)
    for _ in range(5):
        if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 7))[1] < w:
            return c
        c += Fraction(1, 10 ** 7) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 7)
    raise AssertionError("claim not certified")


CLAIM = claim_of(W)
check("claim = least 7-decimal c with 2^c >= R/1000", pow2_interval(CLAIM)[0] >= W
      and pow2_interval(CLAIM - Fraction(1, 10 ** 7))[1] < W, "%s" % format(float(CLAIM), ".7f"))
PER_LEAF = Fraction(R_TOTAL, S)
check("work per leaf = %.7f instructions, at most the declared 454.6389" % float(PER_LEAF),
      PER_LEAF <= Fraction(4546389, 10 ** 4))

# ---------------------------------------------------------------- 4. memory (bytes per node, layout per phase)
CAP_3 = T["cap_3"]
SLOTS = SIDE * CAP_P * 32     # 2^36 pair slots of cap_p records: R1 (hop-1 receive) or R2 (hop-2 receive)
FIXED = [("code, constants, stack", 2 ** 20), ("barrier words: 5 granules of at most 2 KiB", 5 * 2048), ("header", 512), ("X table: x, y", Q1 * 64),
         ("Y table: x, y, x + 2C0", Q2 * 96), ("scalars a, b, c", GENS * 32), ("batch slots", Q2 * 64),
         ("column and row ports, each with its cursor table", 2 * (64 + SIDE * 16)),
         ("final-cursor tables E1 and E2", 2 * SIDE * 32), ("d3 cursor table", 2 ** D3 * 16),
         ("d4 counters", 2 ** D4 * 8), ("candidate list", CAND_CAP * 16),
         ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64), ("relay buffer", RELAY_CAP * 64),
         ("phase-4 count, offset and cursor tables", 3 * SIDE * 32)]
FIXED_BYTES = sum(b for _, b in FIXED)
PHASE3 = [("d3 buckets: 2^24 * cap_3 entries (address, d4)", 2 ** D3 * CAP_3 * 16),
          ("pointer array", 2 ** D3 * CAP_3 * 8)]
PHASE3_BYTES = sum(b for _, b in PHASE3)
for name, b in FIXED + [("R1 at 2^52: hop-1 pair slots (phases 1-2)", SLOTS),
                        ("R2 at 2^52 + |R1|: hop-2 pair slots (phases 2-4)", SLOTS)] + PHASE3:
    print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("fixed regions lie below 2^52", FIXED_BYTES <= 2 ** 52, "2^%.4f" % lg(FIXED_BYTES))
check("bwd's loop end: the q2 64-byte batch slots fill the lower 2^33 bytes of their 2^34-aligned block",
      Q2 * 64 == 2 ** 33)
check("R1 and R2 end below 2^64", 2 ** 52 + 2 * SLOTS < 2 ** 64, "2^%.6f" % lg(2 ** 52 + 2 * SLOTS))
check("phase 3: buckets and pointer array fit in R1's range (R1 is dead after hop 2)", PHASE3_BYTES <= SLOTS,
      "2^%.4f" % lg(PHASE3_BYTES))
NODE_BYTES = FIXED_BYTES + 2 * SLOTS   # phase 2, the largest: fixed regions, R1 and R2
PEAK = NODES * NODE_BYTES
check("peak memory 2^72 * (fixed + R1 + R2) = 2^%.4f bytes" % lg(PEAK), PEAK < 2 ** 136)
check("caps exceed crude means: cap_p > m_p + 2^20, cap_3 > m_3 + 2^31",
      CAP_P > T["m_p"] + 2 ** 20 and CAP_3 > T["m_3"] + 2 ** 31)
check("preprocessing (rows marked P) = 2^%.7f units, included in the total" % lg(Fraction(R_PRE, C)), R_PRE < R_TOTAL)
# the figures declared in proof.md and submission.json: the claim and each supplementary upper bound rounded up (the
# least multiple of 10^-7 or 10^-4 at or above the exact value), the success figure rounded down
DECLARED = dict(time_log2=Fraction(1255620843, 10 ** 7), memory_log2_bytes=Fraction(1359988, 10 ** 4),
                preprocessing_log2=Fraction(1106207, 10 ** 4), success=Fraction(3900001, 10 ** 7))


def least_at_or_above(d, x, scale):
    """Whether 2^d >= x and 2^(d - 1/scale) < x, both certified by the interval bounds of pow2_interval.

    d: Fraction; x: positive int or Fraction; scale: denominator of the decimal grid.
    returns: bool.
    """
    return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x


check("declared time_log2 125.5620843: the claim, least 7-decimal value with 2^c >= R/1000",
      DECLARED["time_log2"] == CLAIM and least_at_or_above(DECLARED["time_log2"], W, 10 ** 7))
check("declared memory 2^135.9988 bytes: least 4-decimal bound >= peak (2^%.7f)" % lg(PEAK),
      least_at_or_above(DECLARED["memory_log2_bytes"], PEAK, 10 ** 4))
check("declared preprocessing 2^110.6207 units: least 4-decimal bound >= rows marked P / 1000",
      least_at_or_above(DECLARED["preprocessing_log2"], Fraction(R_PRE, C), 10 ** 4))
check("declared success 0.3900001 <= the bound", DECLARED["success"] <= T["lo"])
ANALYTIC = sum(v for n, v, _ in rows if n.split(":")[0] in ("program delivery, node and port state", "coins",
               "Y generators [c]P", "R_g = [b]Q and X = R_g + [a]P", "header constants", "batch set-up",
               "Q broadcast", "candidates", "requests served", "phase-4 count exchanges", "barriers", "phase control",
               "output",
               "inversion")) - NODES * Q1 * INV
# itemized work behind the allowances (proof.md section 7)
# on the exact arithmetic of proof.md 3.4: fe_mul_f and fe_sqr_f as executed (115 and 82) plus their calls (pointers and
# BL: 4 and 3, as in fe_inv), and FS or FADDW (9) with both operands brought below p (2 x 9: ADDS chain of C0, CSELs)
MULX, SQRX, ADDX = 115 + 4, 82 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_CTRL = 2 ** 11   # 256 iterations: scalar byte load, bit test, branch, pointer step (6); set-up and return
check("generator (exact arithmetic): doubling 933, mixed addition 1,477, itemized %d + control %d = %d < 2^20"
      % (GEN_ITEM, GEN_CTRL, GEN_ITEM + GEN_CTRL),
      DBL == 933 and MADD_ == 1477 and GEN_ITEM < pow2_interval(Fraction(194, 10))[0]
      and GEN_ITEM + GEN_CTRL <= GEN_MUL and 3 * GEN_MUL > 2 * GEN_ITEM,
      "%d" % GEN_ITEM)
# a barrier at one node: at most 5 waits (4 children, the release), each LDAXR, CMP, B.EQ and at most one WFE, B and
# re-check (8); at most 5 messages (one up, four down), each 3 STPs and 3 for its index and address, plus its in-flight
# count and delivery; DSB; 20 set-up
BARRIER_ITEM = 5 * 8 + 5 * (6 + MSG_LINK) + 1 + 20
check("barrier: at most 131 instructions per node, in-flight counts and deliveries included (2^10 per node is > 7x)",
      BARRIER_ITEM == 131 and BARRIER > 7 * NODES * BARRIER_ITEM)
# a phase-4 count exchange at one node: at most 12 per entry of its 2^36-entry received-count table (LDP, epoch CMP,
# CSEL, prefix ADD, zero test, cap CMP and branch, loop), and at most 88 per message (its count, the count message,
# the offset reply and the write, each at most 22 with its in-flight count and delivery) for at most 2^24 messages
# (2^21 requests of up to 8 records)
EXCH_ITEM = 12 * SIDE + 4 * (14 + MSG_LINK) * 2 ** 24
check("phase-4 count exchange: < 2^39.6 per node (2^41 is > 2.5x)",
      EXCH_ITEM < pow2_interval(Fraction(396, 10))[0] and 2 * EXCH > 5 * EXCH_ITEM, "2^%.3f" % lg(EXCH_ITEM))
# a broadcast message (program, Q): two LDPs of its words, three STPs, ADD, SUBS, B.NE, its in-flight count and delivery
check("a broadcast message costs at most 16, in-flight count and delivery included, within its price of 40",
      2 + 3 + 3 + MSG_LINK <= MSG)
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * (INV + 600) + 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, two additions, unit test, solve mod r, [x]P = Q = %d < 2^20, within 2^21 (2.727x)"
      % CAND_ITEM, CAND_ITEM < 2 ** 20 and RESOLVE == 2 ** 21 and RESOLVE > 2 * CAND_ITEM, "%d" % CAND_ITEM)
PROGRAM_BYTES = 2 ** 16
check("program delivery: 2^16 bytes = 2^11 messages over 2 hops at 40 = 2^17.3 < node set-up 2^20",
      PROGRAM_BYTES // 32 * 2 * MSG < NODE_INIT)
check("analytic allowances are below 2^-14.5 of the total", ANALYTIC * pow2_interval(Fraction(145, 10))[1] < R_TOTAL,
      "2^%.3f" % lg(Fraction(ANALYTIC, R_TOTAL)))

print()
print("q1 %d  q2 2^27  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1, N0, lg(S), float(T["mu"]),
                                                                    math.floor(T["lo"] * 10 ** 10) / 10 ** 10))
print("R_total %d instructions  W 2^%.7f  claim %s  per leaf %.7f" % (R_TOTAL, lg(W), format(float(CLAIM), ".7f"),
                                                                      float(PER_LEAF)))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
    format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK), format(float(DECLARED["preprocessing_log2"]), ".4f"),
    lg(Fraction(R_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
    sys.exit(1 if FAILED else 0)
