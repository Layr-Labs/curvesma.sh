"""Exact certificate for proof.md: parameters, success bound, work ledger, claim, memory (integers and Fractions only).
Imported by probe_direct.py, which asserts every diagnostic, the exact score and complete resource bounds."""
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


def up2(x):
    """Least multiple of 1/100 at or above log2(x), exact (x^100 against powers of two), for the tables' display.

    x: positive int or Fraction.
    returns: float.
    """
    f = Fraction(x)
    a, b = f.numerator ** 100, f.denominator ** 100

    def le(k):
        """Whether x^100 <= 2^k. k: integer. returns: bool."""
        return a <= b << k if k >= 0 else a << -k <= b
    k = math.ceil(lg(f) * 100)
    while not le(k):
        k += 1
    while le(k - 1):
        k -= 1
    return k / 100


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
BETA_P = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)   # beta, a cube root of 1 mod p
C = 1000                      # instructions per target-group-operation equivalent
NODES = 2 ** 72               # 2^36 x 2^36 array
SIDE = 2 ** 36
Q2 = 2 ** 25                  # Y generators per node
D1, D2, D3, D4 = 36, 36, 24, 48
KEYBITS = 112                 # label bits a candidate matches beyond d1, d2: K2 (64), and d4 (48) by its range
THRESH = Fraction(390000001, 10 ** 9)   # target 0.390000001, 10^-9 above the rules' 0.39

# executed counts (experiment.py parts A and B re-derive each of these by running the kernels)
LEAVES = 6                        # leaves per incidence: x(phi^t(X) + Y) and x(phi^t(X) - Y), t = 0, 1, 2
FWD_INC, FWD_FIX = 99, 5          # leaf.s fwd per incidence; mov, movk, the fixup word's STR and LDR, ret
BWD_INC, BWD_FIX = 1403, 5        # leaf.s bwd per incidence (six leaves); mov, movk, STP and LDP of x30 and the fixup
                                  # word, ret
FIX_SITES = dict(SHX=18, FDX=11, F5=7, FS=4, FADDW=6)   # fixup sites per incidence in fwd and bwd (proof.md 3.4)
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4)     # a fixup's instructions: the rest of the chain and STADD
EMIT = 10                         # keys' emit of a leaf: the whole of hop 1's work for a record
EMIT_CHECK = 5                    # d1, cursor entry address, LDADD of the cursor, CMP with the slot end, B.HS: before a drop
EMIT_INC = LEAVES * EMIT          # 60 per incidence
KEY_INC, KEY_FIX = 266, 8         # leaf.s keys per incidence (six leaves, each on the j = 0 path) incl. the six emits;
                                  # set-up (port offset, x4, three LDPs of G0..G5, LSR of G3, MOV of the mask) and ret
F1_COL, COL4 = 11, 15             # KEY's top column F1 (4 UMULH, 3 ADD, 4 MADD) and the next column as a sum of words
                                  # shifted right by 3 (8 products, 7 ADDs with LSR #3) (proof.md 3.3)
JPATH = 5                         # extra on a leaf whose class index j is 1 or 2: the two-word sum and its complement, the
                                  # meter's ADD and the B back, beyond the j = 0 path's
INV = 23646                       # secp_fe.s fe_inv, executed
TINIT_ENTRY = 3                   # route.s tinit per cursor-table entry
FIN_DEST = 9                      # route.s fin per destination: LDR, CMP and CSEL (clamp at the slot end), three STPs,
                                  # ADD, SUBS, B.NE
H2S_REC,H2S_DROP,H2S_SLOT=7,0,11
H2S_FIX=2
SCAN_RANGE=26
SCAN_FIX=SCAN_RANGE*2**D4+2
BKT_REC=BKT_DROP=BKT_SLOT=BKT_FIX=0  # obsolete bkt is never executed
PASS4_REC,PASS4_BKT=6,11
SCAN_REC=Fraction(195,64)
LOOP8=2
ZERO4_PER4,PREFIX4_PER4,INIT3_BKT=4,10,5
check("key pass: 266 = 6*(2 LDP + F1 11 + column 15 + ADD, CMP, B.HS + ADDS, B.HS + ADD of c6 + emit 10) + entry step "
      "and TBZ 2", KEY_INC == LEAVES * (2 + F1_COL + COL4 + 3 + 2 + 1) + EMIT_INC + 2)
check("j = 1 path (ADDS, B.HS, ADD, ADDS, ADC, MVN; emit; meter ADD; B) and j = 2 path (ADDS, ADC, ADD, MVN; emit; ADD; B) "
      "each exceed the j = 0 path (ADDS, B.HS, ADD; emit) by 5", (6 + EMIT + 2) - (3 + EMIT) == JPATH
      and (2 + 4 + EMIT + 2) - (3 + EMIT) == JPATH)
check("emit: 5 up to the slot test (LDADD advances the cursor), then d2, w3 (leaf type and index from the leaf "
      "pointer) and three STPs", EMIT == EMIT_CHECK + 5)
DONE = 2   # per message: the port's in-flight count, up when the link takes it and down when it is written
DLV = 6    # per message: its delivery, route.s dlv (three LDPs of the message, two STPs at its address, RET)
MSG_LINK = DONE + DLV   # link-side work of one message at one hop, beside the sender's own instructions
L8 = Fraction(LOOP8, 8)   # main-loop control per record, an upper bound: a slot of n records runs floor(n/8) passes
HOP1,HOP2=EMIT+MSG_LINK,H2S_REC+L8+15+DONE
GROUP=PASS4_REC+L8+SCAN_REC
check('direct transport and grouping per-record charge',HOP1+HOP2+GROUP==Fraction(3299,64))

# priced allowances (analytic, per item)
GEN_MUL = 669889          # one complete scalar multiplication [k]P or [k]Q with conversion, the point's table entry
                           # (for Y: y^2, y, x, 3x^2, x + C0; for X: the eight header values) and its write
COIN_SCALAR = 64           # at most two draws of 4 RNDR words (18 with one charged retry each) plus a range test (7): 50
BATCH_SETUP = 128          # copy X's eight header values into H, y_X^2 registers, acc = 1, pointers, index word, both
                           # meters' updates and tests, loop (itemized below)
INV_CALL = 16              # store acc, pointer set-up, BL, load R
NODE_INIT = 2 ** 20        # program <= 2^16 bytes = 2^11 messages over 2 hops at 40 (2^17.3), plus node state
HEADER = 27                # G0..G5: 6 words (MOVZ, 3 MOVK each) and 3 STPs
MSG = 40                   # a message outside the two hops (program, Q, phase 4), per hop, with its in-flight count
                           # and delivery (a broadcast message: 16 in all)
CAND_CAP = 2 ** 28         # candidate pairs resolved per home node
RESOLVE = 1538117          # per candidate: two requests and two replies of <= 8 records over 2 hops at 40, two exact
                           # additions, unit test, solve, [x]P = Q: twice the 769,058 itemized, plus 1 (from 1fbd26e8)
REQ_CAP = 2 ** 29          # requests answered per source node
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 34        # relayed records held per node in each phase-4 sub-phase
EXCH = 2 ** 41             # one phase-4 count exchange per node (4 per run)
BARRIERS, BARRIER = 64, 2 ** 82   # 2^10 per node and barrier
OUTPUT = 2 ** 20
DELTA_J = 2 ** 42          # meter slack for the j = 1, 2 key paths (proof.md 3.2)
DELTA_F = 2 ** 45          # cap on a node's arithmetic fixups (proof.md 3.4)
FIX_BATCH = sum(FIX_SITES.values()) * Q2   # most fixups one batch can take: every site of every incidence
CAP_P = 7 * 2 ** 19  # records per pair slot, at hop 1 and at hop 2

# ---------------------------------------------------------------- 2. success
P1 = Fraction(6 * (r - 1) + 1, r * r)    # P(B in class(A)) for independent uniform A, B
TWO_BOUND = Fraction(2, r - 1)           # x-coordinate of a point uniform on at least r - 1 points: each value <= 2/(r-1)

# batch arithmetic sites (proof.md 3.4). Only row shifts (SHX), folds (FDX, also after the fused subtraction of LT), the
# square folds (F5) and FS on a fold output can drop a carry or a borrow or return a value other than the canonical one. Given X != O and every other generator, the fresh Y_l puts each operand value of a site at
# probability at most MULT/(r - 1), MULT being the number of points Y_l giving that value (two per x_Y, three per y_Y,
# and the degree of the map from x_Y or y_Y); a site fails with probability <= BAD * MULT / (r - 1).
SHX_BAD = 2 ** 64 * 2 ** 33     # rows B with B mod 2^192 > 2^192 - 2^33: needed for a row-shift carry
FOLD_BAD = 2 ** 66 * C0         # values in [C0, 2^66 C0): needed for the carry out of a fold of a sum below 2^322
SQ_BAD = C0 * C0 + 2 * C0       # squares in [C0, C0^2 + 2 C0): needed for the carry of a square fold (fifth word <= C0 + 1)
SQC_BAD = C0 * C0 + 3 * C0      # Q_0, Q_1 in [C0, C0^2 + 3 C0): the same with c_t added (fifth word <= C0 + 2)
FSQ_BAD = C0                    # Q_0 or Q_1 below C0: needed for the fold to return it as Q + p, and FS to wrap
LT_BAD = (2 ** 66 + 2) * C0     # leaves below (2^66 + 2) C0: fused fold carry (sum below 2^322 + 2^256) or a leaf >= p
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
]
SITES_PER_INCIDENCE = dict(SHX=18, FDX=11, F5=7)   # per pass of the fwd and bwd loops (experiment.py counts them)
EXACT_AT_L0 = dict(SHX=3, FDX=3)   # rows of acc_(-1) = 1, folds of D_0 * 1, R_0 * 1 and R_0 D_0 = 1 (never carry)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)
# arithmetic fixups (proof.md 3.4): a site takes its fixup only if its exact operand, a row it shifts, or its exact
# output, as a word below 2^256, lies in a set FIX_BAD of words, each congruent to one value mod p; with the same fresh
# generator and points per value as the site table, the mean number of fixups of a batch without a defect is E_FIX
FIX_ROW = C0 * 2 ** 192         # rows whose low word exceeds 2^64 - C0: needed for SHX's carry out of its first ADDS
FIX_FOLD = 2 ** 194 * C0        # outputs with low 128 bits below 2^66 C0: needed for FDX's carry into limb 2
FIX_SQ = 2 ** 193               # outputs with low 128 bits below 2^65: needed for F5's carry into limb 2
FIX_SUB = C0 * 2 ** 192         # outputs whose low word is at least 2^64 - C0 (FS) or below C0 (FADDW)
FIXUPS = [  # (site, macros, executions per batch of Q2 incidences, words that can take the fixup, points of Y_l per value)
    ("fwd: FS forming D_l", "FS", Q2, FIX_SUB, 6),
    ("fwd: row shifts of acc_(l-1), l >= 1", "SHX", 3 * (Q2 - 1), FIX_ROW, 6),
    ("fwd: fold of acc_l, l >= 1", "FDX", Q2 - 1, FIX_FOLD, 6),
    ("bwd: row shifts of R_l", "SHX", 3 * Q2, FIX_ROW, 6),
    ("bwd: folds of inv_l and R_(l-1), l >= 1", "FDX", 2 * (Q2 - 1), FIX_FOLD, 6),
    ("bwd: FS of y_Y - y_X", "FS", Q2, FIX_SUB, 3),
    ("bwd: FADDW of y_Y + y_X", "FADDW", Q2, FIX_SUB, 3),
    ("bwd: row shifts of inv_l", "SHX", 3 * Q2, FIX_ROW, 6),
    ("bwd: folds of u and v", "FDX", 2 * Q2, FIX_FOLD, 3),
    ("bwd: square folds of u^2 and v^2", "F5", 2 * Q2, FIX_SQ, 6),
    ("bwd: FADDW of x_Y + beta^t x_X/2, t = 0, 1", "FADDW", 2 * Q2, FIX_SUB, 2),
    ("bwd: square folds of Q_0 and Q_1", "F5", 2 * Q2, FIX_SQ, 4),
    ("bwd: two FS forming Q_2", "FS", 2 * Q2, FIX_SUB, 4),
    ("bwd: square folds of Q_t^2", "F5", 3 * Q2, FIX_SQ, 8),
    ("bwd: row shifts of Q_t^2", "SHX", 9 * Q2, FIX_ROW, 8),
    ("bwd: FADDW of s'' = beta^t x_X + (x_Y + C0)", "FADDW", 3 * Q2, FIX_SUB, 2),
    ("bwd: fused folds of the six leaves", "FDX", 6 * Q2, FIX_FOLD, 2),
]
E_FIX = sum(Fraction(n * bad * mult, r - 1) for _, _, n, bad, mult in FIXUPS)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (2 ** 19 + 1), (r - 1) * 2 ** 19)  # P(fixed label bits of a leaf
                              # take given values): three paths, two points per x, the lattice count within 1 + 2^-19
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
    d_wj = TWO_BOUND * (13 * 2 ** 197 + 2 ** 64 * C0)   # KEY's index may differ from argmin (window set and C0 case)
    var_j = (Fraction(2, 9) + 2 * d_wj) * n0 + 1
    d_j = var_j / (DELTA_J - n0 // q1 - 1 - 2 * n0 * d_wj) ** 2 + q1 * d_batch
    d_f = q1 * E_FIX / (DELTA_F - FIX_BATCH) + q1 * d_batch   # fixup meter (Markov), or a defective batch of the node
    d_lab = TWO_BOUND * ((70 + 72) * 2 ** 228 + 2 ** 34 * 2 ** 100)   # the label-error set D (proof.md 3.3): T1 mod
                                               # 2^28 below 70, T2 mod 2^28 in 72 residues, or x mod 2^156 > 2^156 - 2^34
    # caps (Chebyshev): count means of the other leaves given the pair's four generators; the variance adds, for the
    # 2(q1 + q2) incidences through those generators, 2 * 15 pairs of one incidence at most one leaf's mean each
    q_phase = cdiv(q1,2)
    m_p = 6*q_phase*Q2 * KB(D1) + 1        # a node's leaves on d1 (hop 1), or a column's leaves on (d1, d2) (hop 2): equal
    v_p = m_p + 1 + 60 * (q_phase + Q2) * KB(D1)
    assert CAP_P > FIXED_LEAVES+m_p
    d_p = v_p / Fraction(CAP_P - FIXED_LEAVES - m_p) ** 2
    m_3 = S * KB(D1 + D2 + D3) + 1
    cap_3 = cdiv(m_3.numerator, m_3.denominator) + 2 ** 32
    v_3 = m_3 + 1 + 60 * (q1 + Q2) * KB(D1 + D2 + D3)
    d_3 = v_3 / Fraction(cap_3 - FIXED_LEAVES - m_3) ** 2
    # intruders and pair loads (Markov): every pair of leaves other than the colliding pair itself contributes at most
    # RHO times its independent-leaf probability (proof.md 6.4: conditioning, or the other node's independence)
    kb_m = KB(KEYBITS + D1 + D2)   # 184 bits: d1, d2, K2 and d4
    lam_int = RHO * (S + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
    other_h = RHO * (Fraction(S) * S * KB(D1 + D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m)
    d_cl = other_h / (CAND_CAP - 4)
    other_g = RHO * Fraction(n0) * S * kb_m * 2
    d_rq = other_g / (REQ_CAP - 8)
    other_v = RHO * (2 + 2 * 8) * (Fraction(S, SIDE) * S * KB(D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m)  # 1 or <= 8
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_wj + d_j + d_f + d_lab + d_p + d_3 + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)  # RNDR is the rules' ideal coin source: its failure branch is charged, never taken
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_wj=d_wj, d_j=d_j, d_f=d_f, d_lab=d_lab,
                d_p=d_p, d_3=d_3, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def, rows=rows, glob=glob,
                lo=lo, cap_3=cap_3, m_p=m_p, m_3=m_3)


def ledger(Q1):
    T=success_terms(Q1)
    N0,S=T["n0"],T["S"]
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
    
    
    row("program delivery, node and port state", NODES * NODE_INIT, True)
    row("coins", NODES * GENS * COIN_SCALAR, True)
    row("Y generators [c]P with their table entries", NODES * Q2 * GEN_MUL, True)
    row("R_g = [b]Q and X = R_g + [a]P with their header values", NODES * (Q1 + 1) * GEN_MUL)
    row("header constants", NODES * HEADER, True)
    row("initial hop1 cursor tables", NODES * (TINIT_ENTRY * SIDE + 1), True)
    row("second phase cursor reset", NODES * (TINIT_ENTRY * SIDE + 1))
    row("Q broadcast", 2 * (SIDE + NODES) * MSG)
    row("d4 counters zeroed", NODES * (2 ** D4 // 4 * ZERO4_PER4 + 1), True)
    row("d3 cursor table", NODES * (5 * 2 ** D3 + 1), True)
    row("forward", NODES * Q1 * (FWD_INC * Q2 + FWD_FIX))
    row("inversion", NODES * Q1 * (INV + INV_CALL))
    row("backward", NODES * Q1 * (BWD_INC * Q2 + BWD_FIX))
    row("arithmetic fixups", NODES * DELTA_F * max(FIX_EXTRA.values()))
    row("key pass without emits", NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX))
    row("j = 1 and j = 2 key paths", NODES * CAP_J * JPATH)
    row("batch set-up", NODES * Q1 * BATCH_SETUP)
    row("hop 1, the emits, in-flight counts and deliveries", NODES * Q1 * (EMIT_INC + LEAVES * MSG_LINK) * Q2)
    row("hop 2, h2s, in-flight counts and deliveries",
        S * (7 + 15 + DONE) + S * LOOP8 // 8 + NODES * 2 * (H2S_SLOT * SIDE + 2))
    row("hop1 final cursors, two phases",
        NODES * 2 * ((FIN_DEST + MSG_LINK) * SIDE + 3))
    row("prefix sums of d4", NODES * (2 ** D4 // 4 * PREFIX4_PER4 + 1))
    row("local pass pass4", S * 6 + S * LOOP8 // 8 + NODES * (PASS4_BKT * 2 ** D3 + 1))
    row("scan_ranges", S * SCAN_REC.numerator // SCAN_REC.denominator
        + NODES * (SCAN_FIX + 8 * (CAND_CAP + 1)))
    row("candidates", NODES * CAND_CAP * (RESOLVE+64))
    row("requests served", NODES * REQ_CAP * REQ_SERVE)
    row("phase-4 count exchanges", NODES * 4 * EXCH)
    row("barriers", 64 * NODES * 2**10)
    row("phase control", NODES * 2 ** 14)
    row("output", OUTPUT)
    R_TOTAL = sum(v for _, v, _ in rows)
    R_PRE = sum(v for _, v, pre in rows if pre)
    W = Fraction(R_TOTAL, C)
    return {"total":R_TOTAL,"preprocessing":R_PRE,"rows":rows,"per_leaf":Fraction(R_TOTAL,S),"log2_display":lg(W)}

MULX, SQRX, ADDX = 115 + 4, 82 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_CTRL = 1600   # 256 iterations of scalar bit test, branch and pointer step (at most 6 each), set-up and return
Y_ENTRY = 2 * SQRX + 2 * ADDX + 5 + 10     # y^2, x^2, 3x^2 = x^2 + x^2 + x^2, x + C0 (ADDS chain), ten STPs
X_ENTRY = MADD_ + 6 * MULX + 2 * SQRX + 18  # X = R_g + [a]P; beta x, beta^2 x, x/2, beta x/2 and, from x^2, 3x^2/4 and
                                            # 3beta^2 x^2/4 (products by constants); x^2, y^2; eighteen STPs

assert GEN_ITEM+GEN_CTRL+max(Y_ENTRY,X_ENTRY)==669889

BATCH_ITEM = (2 + 4) + 32 + (1 + 4 + 4 + 4 + 2 + 1) + (3 * 4 + 1) + (3 * 4 + 1 + 2 + 1 + 1 + 1) + 2 * (2 + 1 + 3) \
    + 2 * (1 + 4 + 1 + 1) + (2 + 2 + 1 + 1) + 2
check("batch set-up: at most %d instructions, within its allowance of 128" % BATCH_ITEM, BATCH_ITEM <= BATCH_SETUP)

CAND_ITEM = 4 * 8 * 2 * MSG + 2 * (INV + 600) + 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, two additions, unit test, solve mod r, [x]P = Q = %d; charged above twice that" % CAND_ITEM,
      RESOLVE > 2 * CAND_ITEM, "%.3fx" % (math.floor(10 ** 3 * RESOLVE / CAND_ITEM) / 10 ** 3))
