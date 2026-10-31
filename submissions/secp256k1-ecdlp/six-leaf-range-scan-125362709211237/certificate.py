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
Q2 = 2 ** 25                  # Y generators per node
D1, D2, D3, D4 = 36, 36, 24, 48
KEYBITS = 112                 # stored key words A = K2, B = K3
THRESH = Fraction(3900001, 10 ** 7)   # target, 10^-9 above the rules' 0.39

# executed counts (experiment.py parts A and B re-derive each of these by running the kernels)
LEAVES = 6                        # leaves per incidence: x(phi^t(X) + Y) and x(phi^t(X) - Y), t = 0, 1, 2
FWD_INC, FWD_FIX = 105, 3         # leaf.s fwd per incidence; mov, movk, ret
BWD_INC, BWD_FIX = 1453, 5        # leaf.s bwd per incidence (six leaves); mov, movk, x30 save and restore, ret
EMIT = 10                         # keys' emit of a leaf: the whole of hop 1's work for a record
EMIT_CHECK = 5                    # d1, cursor entry address, LDADD of the cursor, CMP with the slot end, B.HS: before a drop
EMIT_INC = LEAVES * EMIT          # 60 per incidence
KEY_INC, KEY_FIX = 212, 11        # leaf.s keys per incidence (six leaves) incl. the six emits; set-up, x30 save/restore, ret
KMUL_PATH = 84                    # extra on a leaf whose class index j is 1 or 2 (KMUL 82 with its C0 load, counter, branch)
INV = 23646                       # secp_fe.s fe_inv, executed
TINIT_ENTRY = 3                   # route.s tinit per cursor-table entry
FIN_DEST = 9                      # route.s fin per destination: LDR, CMP and CSEL (clamp at the slot end), three STPs,
                                  # ADD, SUBS, B.NE
H2S_REC, H2S_DROP, H2S_SLOT = 10, 9, 11   # route.s h2s per record sent (LDADD advances the cursor), per record dropped
                                          # at a full slot (the end stored back into the cursor), per slot
H2S_FIX, BKT_FIX, SCAN_FIX = 3, 4, 40 * 2 ** D4 + 2     # h2s: port offset, x17 = 32, RET; bkt: three MOVs, RET; scan: prologue, at most
                                          # 7 single comparisons of 10, the first end test and RET
BKT_REC, BKT_DROP, BKT_SLOT = 10, 8, 12   # route.s bkt per record bucketed (LDADD cursor, TST of its low 38 bits, STADD
                                          # d4 count), per record dropped at a full bucket (LDADD undone), per slot
PASS4_REC, PASS4_BKT = 4, 11
SCAN_REC = Fraction(15, 4)                # scan per record: half a pointer LDP, key LDP, CMP, CCMP, B.EQ, a quarter of the
                                          # end test per eight
LOOP8 = 2                                 # h2s, bkt, pass4: CMP and B.LO per 8 records of a slot's or bucket's main loop
ZERO4_PER4, PREFIX4_PER4, INIT3_BKT = 4, 10, 4   # prefix4: two LDPs, two STPs and four ADDs per block of four
check("key pass: 212 = 6*(2 LDP + KEY 21 + B.LO + CBZ) + six emits of 10 + entry step and TBZ 2",
      KEY_INC == LEAVES * (2 + 21 + 1 + 1) + EMIT_INC + 2)
check("emit: 5 up to the slot test (LDADD advances the cursor), then d2, w3 (leaf type and index from the leaf "
      "pointer) and three STPs", EMIT == EMIT_CHECK + 5)
DONE = 2   # per message: the port's in-flight count, up when the link takes it and down when it is written
DLV = 6    # per message: its delivery, route.s dlv (three LDPs of the message, two STPs at its address, RET)
MSG_LINK = DONE + DLV   # link-side work of one message at one hop, beside the sender's own instructions
L8 = Fraction(LOOP8, 8)   # main-loop control per record, an upper bound: a slot of n records runs floor(n/8) passes
HOP1, HOP2 = EMIT + MSG_LINK, H2S_REC + L8 + MSG_LINK
GROUP = BKT_REC + L8 + PASS4_REC + L8 + SCAN_REC
check("per record: hop 1 18, hop 2 18.25 (each with in-flight count 2 and delivery 6), grouping 18.25; 54.5 in all",
      HOP1 + HOP2 + GROUP == Fraction(109, 2))

# priced allowances (analytic, per item)
GEN_MUL = 2 ** 20          # one complete scalar multiplication [k]P or [k]Q with conversion, the point's table entry
                           # (for Y: y^2, y, x, 3x^2, x + C0; for X: the eight header values) and its write
COIN_SCALAR = 64           # at most two draws of 4 RNDR words (18 with one charged retry each) plus a range test (7): 50
BATCH_SETUP = 128          # copy X's eight header values into H, y_X^2 registers, acc = 1, pointers, index word, meter
                           # test, loop
INV_CALL = 16              # store acc, pointer set-up, BL, load R
NODE_INIT = 2 ** 20        # program <= 2^16 bytes = 2^11 messages over 2 hops at 40 (2^17.3), plus node state
HEADER = 156               # G1..G5 and the two row tables, each with its C0 word (39 words at 4 per word)
MSG = 40                   # a message outside the two hops (program, Q, phase 4), per hop, with its in-flight count
                           # and delivery (a broadcast message: 16 in all)
CAND_CAP = 2 ** 28         # candidate pairs resolved per home node
RESOLVE = 2 ** 21          # per candidate: two requests and two replies of <= 8 records over 2 hops at 40, two exact
                           # additions, unit test, solve, [x]P = Q (769,058 itemized)
REQ_CAP = 2 ** 29          # requests answered per source node
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 30        # relayed records held per node in each phase-4 sub-phase
EXCH = 2 ** 41             # one phase-4 count exchange per node (4 per run)
BARRIERS, BARRIER = 32, 2 ** 82   # 2^10 per node and barrier
OUTPUT = 2 ** 20
DELTA_J = 2 ** 42          # meter slack for the class-index products
CAP_P = 2 ** 22 - 2 ** 12  # records per pair slot, at hop 1 and at hop 2

# ---------------------------------------------------------------- 2. success
P1 = Fraction(6 * (r - 1) + 1, r * r)    # P(B in class(A)) for independent uniform A, B
TWO_BOUND = Fraction(2, r - 1)           # x-coordinate of a point uniform on at least r - 1 points: each value <= 2/(r-1)

# batch arithmetic sites (proof.md 3.4). Only row shifts (SHX), folds (FDX, also inside KMUL and after the fused
# subtraction of LT), the square folds (F5) and FS on a fold output can drop a carry or a borrow, and only the folds can
# return a value >= p. Given X != O and every other generator, the fresh Y_l puts each operand value of a site at
# probability at most MULT/(r - 1), MULT being the number of points Y_l giving that value (two per x_Y, three per y_Y,
# and the degree of the map from x_Y or y_Y); a site fails with probability <= BAD * MULT / (r - 1).
SHX_BAD = 2 ** 64 * 2 ** 33     # rows B with B mod 2^192 > 2^192 - 2^33: needed for a row-shift carry
FOLD_BAD = 2 ** 66 * C0         # values in [C0, 2^66 C0): needed for the carry out of a fold of a sum below 2^322
SQ_BAD = C0 * C0 + 2 * C0       # squares in [C0, C0^2 + 2 C0): needed for the carry of a square fold (fifth word <= C0 + 1)
SQC_BAD = C0 * C0 + 3 * C0      # Q_0, Q_1 in [C0, C0^2 + 3 C0): the same with c_t added (fifth word <= C0 + 2)
FSQ_BAD = C0                    # Q_0 or Q_1 below C0: needed for the fold to return it as Q + p, and FS to wrap
LT_BAD = (2 ** 66 + 2) * C0     # leaves below (2^66 + 2) C0: fused fold carry (sum below 2^322 + 2^256) or a leaf >= p
KMUL_BAD = 2 * 2 ** 66 * C0     # leaves with beta^j x below 2^66 C0 for j = 1 or 2: fold carry or K >= p
SITES = [  # (site, executions per batch of Q2 incidences, bad operand values, points of Y_l per value)
    ("fwd: row shifts of acc_(l-1), l >= 1", 3 * (Q2 - 1), SHX_BAD, 6),
    ("fwd: fold of acc_l, l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: row shifts of R_l", 3 * Q2, SHX_BAD, 6),
    ("bwd: fold of inv_l = R_l acc_(l-1), l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: fold of R_(l-1) = R_l D_l, l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: row shifts of inv_l", 3 * Q2, SHX_BAD, 6),
    ("bwd: folds of u = inv_l (y_Y - y_X) and v", 2 * Q2, FOLD_BAD, 3),
    ("bwd: square folds of u^2 and v^2", 2 * Q2, SQ_BAD, 6),
    ("bwd: square folds of Q_0 and Q_1", 2 * Q2, SQC_BAD, 4),
    ("bwd: FS for Q_2 = 3x_Y^2 - Q_0 - Q_1", 2 * Q2, FSQ_BAD, 4),
    ("bwd: square folds of Q_t^2, t = 0, 1, 2", 3 * Q2, SQ_BAD, 8),
    ("bwd: row shifts of Q_t^2", 9 * Q2, SHX_BAD, 8),
    ("bwd: fused folds of the six leaves", 6 * Q2, LT_BAD, 2),
    ("keys: KMUL folds of the six leaves", 6 * Q2, KMUL_BAD, 2),
]
SITES_PER_INCIDENCE = dict(SHX=18, FDX=17, F5=7)   # per pass of the fwd, bwd and keys loops (experiment.py counts them)
EXACT_AT_L0 = dict(SHX=3, FDX=3)   # rows of acc_(-1) = 1, folds of D_0 * 1, R_0 * 1 and R_0 D_0 = 1 (never carry)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed), r - 1)  # P(fixed key bits of a leaf take given values)
RHO = Fraction(r, r - 1) ** 4  # slack for the exclusion Y != O when leaves are uniform on r - 1 points or pairs
FIXED_LEAVES = 2 * LEAVES     # the leaves of the pair's two incidences, fixed by its four generators
PARTNER_PAIRS = 2 * 10        # pairs of two partners of one of the pair's leaves (proof.md 6.4), counted without the
                              # home or relay-column condition


def success_terms(q1):
    """Success lower bound pieces for q1 X generators per node.

    q1: integer.
    returns: dict of exact Fractions.
    """
    n0 = LEAVES * q1 * Q2
    S = NODES * n0
    gamma = Fraction(S * S - NODES * n0 * n0, 2)
    mu = gamma * P1
    deg = 2 * LEAVES * (q1 + Q2) * S       # pairs sharing a generator with a given pair: an X is in 6 q2 leaves, a Y in 6 q1
    cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
    mu_lo = Fraction((mu.numerator << 200) // mu.denominator, 2 ** 200)
    e_hi = exp_interval(-mu_lo, 40)[1]
    # per-leaf / per-pair defect probabilities
    d_batch = D_SITES + Q2 * Fraction(6, r - 1) + Fraction(1, r)   # sites, D_l = 0 (Y_l = +-phi^t(X)), X = O
    d_win = TWO_BOUND * 13 * 2 ** 197
    d_wj = TWO_BOUND * 2 ** 64 * C0
    var_j = Fraction(2, 9) * n0 + 1
    d_j = var_j / Fraction(DELTA_J - n0 // q1 - 1) ** 2 + q1 * d_batch
    # caps (Chebyshev): count means of the other leaves given the pair's four generators; the variance adds, for the
    # 2(q1 + q2) incidences through those generators, 2 * 15 pairs of one incidence at most one leaf's mean each
    m_p = n0 * KB(D1) + 1        # a node's leaves on d1 (hop 1), or a column's leaves on (d1, d2) (hop 2): equal
    v_p = m_p + 1 + 60 * (q1 + Q2) * KB(D1)
    d_p = v_p / Fraction(CAP_P - FIXED_LEAVES - m_p) ** 2
    m_3 = S * KB(D1 + D2 + D3) + 1
    cap_3 = cdiv(m_3.numerator, m_3.denominator) + 2 ** 32
    v_3 = m_3 + 1 + 60 * (q1 + Q2) * KB(D1 + D2 + D3)
    d_3 = v_3 / Fraction(cap_3 - FIXED_LEAVES - m_3) ** 2
    # intruders and pair loads (Markov): every pair of leaves other than the colliding pair itself contributes at most
    # RHO times its independent-leaf probability (proof.md 6.4: conditioning, or the other node's independence)
    kb200 = KB(KEYBITS + D1 + D2)
    lam_int = RHO * (S + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
    other_h = RHO * (Fraction(S) * S * KB(D1 + D2) * kb200 + (2 * S + PARTNER_PAIRS) * kb200)
    d_cl = other_h / (CAND_CAP - 4)
    other_g = RHO * Fraction(n0) * S * kb200 * 2
    d_rq = other_g / (REQ_CAP - 8)
    other_v = RHO * (2 + 2 * 8) * (Fraction(S, SIDE) * S * KB(D2) * kb200 + (2 * S + PARTNER_PAIRS) * kb200)  # 1 or <= 8
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_win + d_wj + d_j + 2 * d_p + d_3 + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)  # RNDR is the rules' ideal coin source: its failure branch is charged, never taken
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_win=d_win, d_wj=d_wj, d_j=d_j,
                d_p=d_p, d_3=d_3, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def, rows=rows, glob=glob,
                lo=lo, cap_3=cap_3, m_p=m_p, m_3=m_3)


lo_q, hi_q = 2 ** 27, 2 ** 28
while hi_q - lo_q > 1:
    mid = (lo_q + hi_q) // 2
    if success_terms(mid)["lo"] >= THRESH:
        hi_q = mid
    else:
        lo_q = mid
Q1 = hi_q
T = success_terms(Q1)
N0, S = T["n0"], T["S"]
check("q1 = %d is the least with success bound >= 0.3900000001" % Q1, T["lo"] >= THRESH and success_terms(Q1 - 1)["lo"] < THRESH)
check("leaves S = 2^72*6*q1*q2 = 2^%.6f" % lg(S), S == NODES * 6 * Q1 * Q2)
check("mu = |Gamma|*P1 = %.9f (cross-node pairs)" % float(T["mu"]), T["mu"] > Fraction(494, 1000))
check("Chen-Stein b1+b2 < 2^-95", T["cs"] < Fraction(1, 2 ** 95), "2^%.2f" % lg(T["cs"]))
check("intruder rate per colliding pair < 2^-14.6", T["lam_int"] < pow2_interval(Fraction(-146, 10))[0], "2^%.3f" % lg(T["lam_int"]))
check("window-drop rate per leaf < 2^-54.29", T["d_win"] < pow2_interval(Fraction(-5429, 100))[0], "2^%.4f" % lg(T["d_win"]))
for name, n, bad, mult in SITES:
    k = -(-n // Q2)
    print("  site %-52s %-6s bad <= 2^%.2f, x%d/(r-1), each <= 2^%.2f" % (name, ("%dq2" % k if k > 1 else "q2") + (
        "%+d" % (n - k * Q2) if n != k * Q2 else ""), ceil_lg(bad), mult, ceil_lg(Fraction(bad * mult, r - 1))))
check("site table: 18 row shifts, 17 folds, 7 square folds per incidence, less 3 + 3 exact at l = 0; 2 FS on folds",
      sum(n for _, n, b, _ in SITES if b == SHX_BAD) == SITES_PER_INCIDENCE["SHX"] * Q2 - EXACT_AT_L0["SHX"]
      and sum(n for _, n, b, _ in SITES if b in (FOLD_BAD, LT_BAD, KMUL_BAD)) == SITES_PER_INCIDENCE["FDX"] * Q2
      - EXACT_AT_L0["FDX"]
      and sum(n for _, n, b, _ in SITES if b in (SQ_BAD, SQC_BAD)) == SITES_PER_INCIDENCE["F5"] * Q2
      and sum(n for _, n, b, _ in SITES if b == FSQ_BAD) == 2 * Q2)
check("bad sets: rows 2^97, folds below 2^98.02, square folds below 2^64.01, KMUL below 2^99.01",
      SHX_BAD == 2 ** 97 and all(b < pow2_interval(Fraction(e, 100))[0] for b, e in ((FOLD_BAD, 9802), (LT_BAD, 9802),
                                                                                     (SQ_BAD, 6401), (SQC_BAD, 6401),
                                                                                     (KMUL_BAD, 9901))))
check("batch defect rate (site table, D = 0, X = O) < 2^-126", T["d_batch"] < Fraction(1, 2 ** 126),
      "2^%.3f" % lg(T["d_batch"]))
check("class-index meter rate < 2^-31", T["d_j"] < Fraction(1, 2 ** 31), "2^%.2f" % lg(T["d_j"]))
check("pair-slot overflow rate, each hop, < 2^-19.9", T["d_p"] < pow2_interval(Fraction(-199, 10))[0],
      "2^%.3f" % lg(T["d_p"]))
check("d3 bucket overflow rate < 2^-30", T["d_3"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_3"]))
check("candidate cap rate < 2^-22, request cap rate < 2^-24",
      T["d_cl"] < pow2_interval(Fraction(-22, 1))[0] and T["d_rq"] < pow2_interval(Fraction(-24, 1))[0],
      "2^%.3f 2^%.3f" % (lg(T["d_cl"]), lg(T["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-20", T["d_rl"] < pow2_interval(Fraction(-20, 1))[0],
      "2^%.3f" % lg(T["d_rl"]))
check("whole-run row: a computational defect in any batch < 2^-26.9", T["g_def"] < pow2_interval(Fraction(-269, 10))[0],
      "2^%.3f" % lg(T["g_def"]))
check("whole-run rows (defective batch anywhere, coin rejection twice) < 2^-26.9",
      T["glob"] < pow2_interval(Fraction(-269, 10))[0], "2^%.3f" % lg(T["glob"]))
check("success >= 1 - e^-mu - (b1+b2) - rows - global >= 0.3900000001", T["lo"] >= THRESH,
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
row("Y generators [c]P with their table entries: q2 * 2^20 per node", NODES * Q2 * GEN_MUL, True)
row("R_g = [b]Q and X = R_g + [a]P with their header values: (q1+1) * 2^20 per node", NODES * (Q1 + 1) * GEN_MUL)
row("header constants: 156 per node", NODES * HEADER, True)
row("cursor tables: 2 * (3 * 2^36 + 1) per node", NODES * 2 * (TINIT_ENTRY * SIDE + 1), True)
row("Q broadcast: 2 messages to every node over two hops at 40", 2 * (SIDE + NODES) * MSG)
row("d4 counters zeroed: 2^48/4 * 4 + 1 per node", NODES * (2 ** D4 // 4 * ZERO4_PER4 + 1), True)
row("d3 cursor table: 4 * 2^24 + 1 per node", NODES * (INIT3_BKT * 2 ** D3 + 1), True)
row("forward: q1 * (105 q2 + 3)", NODES * Q1 * (FWD_INC * Q2 + FWD_FIX))
row("inversion: q1 * (23,646 + 16)", NODES * Q1 * (INV + INV_CALL))
row("backward: q1 * (1453 q2 + 5)", NODES * Q1 * (BWD_INC * Q2 + BWD_FIX))
row("key pass without emits: q1 * (152 q2 + 11)", NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX))
row("class-index products: cap_j * 84 per node", NODES * CAP_J * KMUL_PATH)
row("batch set-up: q1 * 128", NODES * Q1 * BATCH_SETUP)
row("hop 1, the emits, in-flight counts and deliveries: q1 * 108 q2", NODES * Q1 * (EMIT_INC + LEAVES * MSG_LINK) * Q2)
row("hop 2, h2s, in-flight counts and deliveries: S * 18.25 + 11 * 2^36 + 3 per node",
    S * (H2S_REC + MSG_LINK) + S * LOOP8 // 8 + NODES * (H2S_SLOT * SIDE + H2S_FIX))
row("final cursors, in-flight counts and deliveries, both hops: 2 * (17 * 2^36 + 3) per node",
    NODES * 2 * ((FIN_DEST + MSG_LINK) * SIDE + 3))
row("bucket pass bkt: S * 10.25 + 12 * 2^36 + 4 per node", S * BKT_REC + S * LOOP8 // 8 + NODES * (BKT_SLOT * SIDE + BKT_FIX))
row("prefix sums of d4: 2^48/4 * 10 + 1 per node", NODES * (2 ** D4 // 4 * PREFIX4_PER4 + 1))
row("local pass pass4: S * 4.25 + 11 * 2^24 + 1 per node", S * PASS4_REC + S * LOOP8 // 8 + NODES * (PASS4_BKT * 2 ** D3 + 1))
row("scan: S * 4.75 + 79 + 4 * (2^20 + 1) per node", S * SCAN_REC.numerator // SCAN_REC.denominator
    + NODES * (SCAN_FIX + 4 * (CAND_CAP + 1)))
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
    """Smallest multiple of 10^-9 whose power of two is at least w, verified by interval bounds.

    w: Fraction > 1.
    returns: Fraction claim.
    """
    c = Fraction(math.ceil(lg(w) * 10 ** 12), 10 ** 12)
    for _ in range(5):
        if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 12))[1] < w:
            return c
        c += Fraction(1, 10 ** 12) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 12)
    raise AssertionError("claim not certified")


CLAIM = claim_of(W)
check("claim = least 12-decimal c with 2^c >= R/1000", pow2_interval(CLAIM)[0] >= W
      and pow2_interval(CLAIM - Fraction(1, 10 ** 12))[1] < W, "%s" % format(float(CLAIM), ".12f"))
PER_LEAF = Fraction(R_TOTAL, S)
check("work per leaf = %.7f instructions" % float(PER_LEAF), 395 < PER_LEAF < 396)

# ---------------------------------------------------------------- 4. memory (bytes per node, layout per phase)
CAP_3 = T["cap_3"]
SLOTS = SIDE * CAP_P * 32     # 2^36 pair slots of cap_p records: R1 (hop-1 receive) or R2 (hop-2 receive)
FIXED = [("code, constants, stack", 2 ** 20), ("barrier words: 5 granules of at most 2 KiB", 5 * 2048),
         ("header", 1024), ("X table: the eight header values and y^2 per X", Q1 * 288),
         ("Y table: y^2, y, x, 3x^2, x + C0", Q2 * 160), ("scalars a, b, c", GENS * 32), ("batch entries", Q2 * 256),
         ("column and row ports, each with its cursor table", 2 * (64 + SIDE * 8)),
         ("final-cursor tables E1 and E2", 2 * SIDE * 32), ("d3 cursor table", 2 ** D3 * 8),
         ("d4 counters", 2 ** D4 * 8), ("candidate list", CAND_CAP * 16),
         ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64), ("relay buffer", RELAY_CAP * 64),
         ("phase-4 count, offset and cursor tables", 3 * SIDE * 32)]
FIXED_BYTES = sum(b for _, b in FIXED)
BUCKET_STRIDE = 2 ** 38      # bucket t's cap_3 entries of 16 bytes end at R1 + (t + 1)*2^38 (bkt's test of 38 low bits)
PHASE3 = [("d3 buckets: 2^24 strides of 2^38 bytes, cap_3 entries (address, d4) at the end of each",
           2 ** D3 * BUCKET_STRIDE),
          ("pointer array", 2 ** D3 * CAP_3 * 8)]
PHASE3_BYTES = sum(b for _, b in PHASE3)
for name, b in FIXED + [("R1 at 2^52: hop-1 pair slots (phases 1-2)", SLOTS),
                        ("R2 at 2^52 + |R1|: hop-2 pair slots (phases 2-4)", SLOTS)] + PHASE3:
    print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("fixed regions lie below 2^52", FIXED_BYTES <= 2 ** 52, "2^%.4f" % lg(FIXED_BYTES))
check("bwd's and keys' loop end: the q2 256-byte batch entries fill the lower 2^33 bytes of their 2^34-aligned block",
      Q2 * 256 == 2 ** 33)
check("record fields: k < 2^28 in w2 bits 36-63; leaf type 1..6 in w3 bits 36-38 and l < 2^25 in bits 39-63",
      Q1 <= 2 ** 28 and Q2 <= 2 ** 25)
check("R1 and R2 end below 2^64", 2 ** 52 + 2 * SLOTS < 2 ** 64, "2^%.6f" % lg(2 ** 52 + 2 * SLOTS))
check("phase 3: buckets and pointer array fit in R1's range (R1 is dead after hop 2)", PHASE3_BYTES <= SLOTS,
      "2^%.4f" % lg(PHASE3_BYTES))
check("d3 buckets: cap_3 entries of 16 bytes fit below a stride of 2^38, and R1 is a multiple of the stride",
      16 * CAP_3 < BUCKET_STRIDE and 2 ** 52 % BUCKET_STRIDE == 0, "16*cap_3 = 2^%.4f" % lg(16 * CAP_3))
check("hop-1 cursors cannot wrap: the last slot's base plus 32 per leaf of a node stays below 2^64 (hop-2 and bucket "
      "cursors stop at their end)", 2 ** 52 + 32 * CAP_P * (SIDE - 1) + 32 * N0 < 2 ** 64,
      "2^%.4f" % lg(2 ** 52 + 32 * CAP_P * (SIDE - 1) + 32 * N0))
NODE_BYTES = FIXED_BYTES + 2 * SLOTS   # phase 2, the largest: fixed regions, R1 and R2
PEAK = NODES * NODE_BYTES
check("peak memory 2^72 * (fixed + R1 + R2) = 2^%.4f bytes" % lg(PEAK), PEAK < 2 ** 136)
check("caps exceed crude means: cap_p > m_p + 2^20, cap_3 > m_3 + 2^31",
      CAP_P > T["m_p"] + 2 ** 20 and CAP_3 > T["m_3"] + 2 ** 31)
check("preprocessing (rows marked P) = 2^%.7f units, included in the total" % lg(Fraction(R_PRE, C)), R_PRE < R_TOTAL)
# the figures declared in proof.md and submission.json: the claim and each supplementary upper bound rounded up (the
# least multiple of 10^-9 or 10^-4 at or above the exact value), the success figure rounded down
DECLARED = dict(time_log2=Fraction(125362709211237, 10 ** 12), memory_log2_bytes=ceil_lg(PEAK,10**4),
                preprocessing_log2=ceil_lg(Fraction(R_PRE,C),10**4), success=THRESH,
                per_leaf=Fraction((PER_LEAF.numerator*10000+PER_LEAF.denominator-1)//PER_LEAF.denominator,10000))


def least_at_or_above(d, x, scale):
    """Whether 2^d >= x and 2^(d - 1/scale) < x, both certified by the interval bounds of pow2_interval.

    d: Fraction; x: positive int or Fraction; scale: denominator of the decimal grid.
    returns: bool.
    """
    return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x


check("declared time_log2 %s: claim, least 12-decimal value with 2^c >= R/1000" % format(float(DECLARED["time_log2"]), ".12f"),
      DECLARED["time_log2"] == CLAIM and least_at_or_above(DECLARED["time_log2"], W, 10 ** 12))
check("declared memory 2^135.9988 bytes: least 4-decimal bound >= peak (2^%.7f)" % lg(PEAK),
      least_at_or_above(DECLARED["memory_log2_bytes"], PEAK, 10 ** 4))
check("declared preprocessing 2^%s units: least 4-decimal bound >= rows marked P / 1000" % format(float(DECLARED["preprocessing_log2"]), ".4f"),
      least_at_or_above(DECLARED["preprocessing_log2"], Fraction(R_PRE, C), 10 ** 4))
check("declared success 0.3900001 <= the bound", DECLARED["success"] <= T["lo"])
check("declared work per leaf 395.9575: least 4-decimal bound >= R/S (%.7f)" % float(PER_LEAF),
      DECLARED["per_leaf"] >= PER_LEAF > DECLARED["per_leaf"] - Fraction(1, 10 ** 4))
ANALYTIC = sum(v for n, v, _ in rows if n.split(":")[0] in ("program delivery, node and port state", "coins",
               "Y generators [c]P with their table entries", "R_g = [b]Q and X = R_g + [a]P with their header values",
               "header constants", "batch set-up",
               "Q broadcast", "candidates", "requests served", "phase-4 count exchanges", "barriers", "phase control",
               "output",
               "inversion")) - NODES * Q1 * INV
# itemized work behind the allowances (proof.md section 7)
# on the exact arithmetic of proof.md 3.4: fe_mul_f and fe_sqr_f as executed (115 and 82) plus their calls (pointers and
# BL: 4 and 3, as in fe_inv), and FS or FADDW (9) with both operands brought below p (2 x 9: ADDS chain of C0, CSELs)
MULX, SQRX, ADDX = 115 + 4, 82 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_CTRL = 2 ** 11   # 256 iterations of scalar bit test, branch and pointer step (at most 6 each), set-up and return
Y_ENTRY = 2 * SQRX + 2 * ADDX + 5 + 10     # y^2, x^2, 3x^2 = x^2 + x^2 + x^2, x + C0 (ADDS chain), ten STPs
X_ENTRY = MADD_ + 6 * MULX + 2 * SQRX + 18  # X = R_g + [a]P; beta x, beta^2 x, x/2, beta x/2 and, from x^2, 3x^2/4 and
                                            # 3beta^2 x^2/4 (products by constants); x^2, y^2; eighteen STPs
check("generator (exact arithmetic): doubling 933, mixed addition 1,477; with control, the Y entry and the X header "
      "values at most %d < 2^20" % (GEN_ITEM + GEN_CTRL + max(Y_ENTRY, X_ENTRY)),
      DBL == 933 and MADD_ == 1477 and GEN_ITEM + GEN_CTRL + max(Y_ENTRY, X_ENTRY) <= GEN_MUL,
      "%.4fx" % (math.floor(10 ** 4 * GEN_MUL / (GEN_ITEM + GEN_CTRL + max(Y_ENTRY, X_ENTRY))) / 10 ** 4))
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
EXCH_ITEM = 12 * SIDE + 4 * (14 + MSG_LINK) * (8 * REQ_CAP)
check("phase-4 count exchange: below charged 2^41 per node",
      EXCH_ITEM < EXCH, "2^%.3f" % lg(EXCH_ITEM))
# a broadcast message (program, Q): two LDPs of its words, three STPs, ADD, SUBS, B.NE, its in-flight count and delivery
check("a broadcast message costs at most 16, in-flight count and delivery included, within its price of 40",
      2 + 3 + 3 + MSG_LINK <= MSG)
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * (INV + 600) + 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, two additions, unit test, solve mod r, [x]P = Q = %d, within 2^21" % CAND_ITEM,
      RESOLVE > 2 * CAND_ITEM, "%.3fx" % (math.floor(10 ** 3 * RESOLVE / CAND_ITEM) / 10 ** 3))
PROGRAM_BYTES = 2 ** 16
check("program delivery: 2^16 bytes = 2^11 messages over 2 hops at 40 = 2^17.3 < node set-up 2^20",
      PROGRAM_BYTES // 32 * 2 * MSG < NODE_INIT)
check("analytic allowances are below 2^-12 of the total", ANALYTIC * 2**12 < R_TOTAL,
      "2^%.3f" % lg(Fraction(ANALYTIC, R_TOTAL)))

print()
print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1, N0, lg(S), float(T["mu"]),
                                                                    math.floor(T["lo"] * 10 ** 10) / 10 ** 10))
print("R_total %d instructions  W 2^%.12f  claim %s  per leaf %.4f" % (R_TOTAL, lg(W), format(float(CLAIM), ".12f"),
                                                                      float(PER_LEAF)))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
    format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK), format(float(DECLARED["preprocessing_log2"]), ".4f"),
    lg(Fraction(R_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
    sys.exit(1 if FAILED else 0)
