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


# ---------------------------------------------------------------- 1. constants
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
C0 = 2 ** 32 + 977
BETA_P = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)   # beta, a cube root of 1 mod p
C = 1000                      # instructions per target-group-operation equivalent
NODES = 2 ** 72               # 2^36 x 2^36 array
SIDE = 2 ** 36
Q2 = 2 ** 25                  # Y generators per node
D1, D2, D3, D4 = 36, 36, 24, 46
KEYBITS = 110                 # label bits a candidate matches beyond d1, d2: K2 (64), and d4 (46) by its range
THRESH = Fraction(390000001, 10 ** 9)   # target 0.390000001, 10^-9 above the rules' 0.39

# executed counts (experiment.py parts A and B re-derive each of these by running the kernels)
LEAVES = 6                        # leaves per incidence: x(phi^t(X) + Y) and x(phi^t(X) - Y), t = 0, 1, 2
FWD_INC, FWD_FIX = Fraction(195, 2), 5          # leaf.s fwd per incidence (two per pass: TBZ x3, #33); mov, movk, the fixup word's
                                  # STR and LDR, ret
BWD_INC, BWD_FIX = Fraction(2413, 2), 5   # leaf.s bwd per incidence (six components; TBZ per two); mov, movk, STP and
                                  # LDP of x30 and the fixup word, ret
FIX_SITES = dict(SHX=18, FDX=14, F5=2, FS=2, FADDW=1)   # fixup sites per incidence in fwd and bwd (proof.md 3.4)
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4)     # a fixup's instructions: the rest of the chain and STADD
EMIT = 7                          # keys' emit of a leaf into its destination's record block (proof.md 3.2)
EMIT_INC = LEAVES * EMIT          # 48 per incidence
KEY_INC, KEY_FIX = Fraction(657, 2), 11   # leaf.s keys per incidence (six leaves on the j = 0 path, TBZ per two) incl.
                                  # the emits; set-up (STR of x30, the mask, x4 = 32, x0 = 32k, three LDPs and an LDR of
                                  # G0..G5 and G3 >> 3, MOV of 2^36 - 1), LDR of x30 and ret
K = 256                           # records per block message (proof.md 2.1)
PB1C, PB2C, D4BASE, BLKC = 1 << 50, 1 << 51, 1 << 49, 1 << 14   # column blocks, row blocks, d4 counters, stride
FLUSH = 8                         # after slot k: two SUBs (8160 = 2*4080), STR (LC back to slot 1), ORR (RC entry),
                                  # B.LO, STP of the header: a block sent, in keys and in h2s
KEY_BLOCK_DROP, H2S_BLOCK_DROP = 7, 9   # a block dropped: without the STP; in h2s the end stored back (STR) and B
F1_COL, COL4 = 11, 15             # KEY's top column F1 (4 UMULH, 3 ADD, 4 MADD) and the next column as a sum of words
                                  # shifted right by 3 (8 products, 7 ADDs with LSR #3) (proof.md 3.3)
JPATH = 5                         # extra on a leaf whose class index j is 1 or 2: the two-word sum and its complement, the
                                  # meter's ADD and the B back, beyond the j = 0 path's
INV = 23646                       # secp_fe.s fe_inv, executed
TINIT_ENTRY = 3                   # route.s tinit per RC-table entry
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 14, 3, 2, 2   # route.s fin per destination: LDR of LC and RC, SUBS, B.EQ,
                                  # ADD, CMP, CSEL, three STPs, two ADDs, SUBS, B.NE; a final block sent CMP, B.HS, STP
                                  # (held: CMP, B.HS); per call MOV and RET
H2S_REC, H2S_SLOT = 7, 11                 # route.s h2s per record (LDADD of its LC entry returns its slot), per slot
SCAN_RANGE = 25                           # scan_ranges: at most 25 per d4 range above 195/64 per record (proof.md 4)
H2S_FIX, BKT_FIX, SCAN_FIX = 4, 4, SCAN_RANGE * 2 ** D4 + 2   # h2s: x10 = 2^39, x17 = 32, x18 = 32k, RET; bkt: three MOVs,
                                          # RET; scan_ranges: its ranges, the entry CBZ and RET
BKT_REC, BKT_DROP, BKT_SLOT = 8, 7, 12    # route.s bkt per record bucketed (BFI addresses, LDADD cursor, TST of its low 38
                                          # bits, STADD d4 count), per record dropped at a full bucket (LDADD undone), per slot
PASS4_REC, PASS4_BKT, PASS4_FIX = 3, 11, 2   # pass4 per entry (LDP, LDADD of its d4 cursor, STR), per bucket, MOV and RET
SCAN_REC = Fraction(195, 64)              # scan_ranges per 64 records: 32 pointer LDPs, 64 LDRs of K2, 32 CMP and 32 CCMP,
                                          # 32 B.EQ and the end test (3) (from 892a3df8 and b00400fa)
LOOP8 = 2                                 # h2s, bkt, pass4: CMP and B.LO per 8 records of a slot's or bucket's main loop
ZERO4_PER4, PREFIX4_PER4, INIT3_BKT = 4, 10, 4   # prefix4: two LDPs, two STPs and four ADDs per block of four;
                                                 # init3 per entry (d3 cursor table, LC table)
check("key pass: 328.5 = 2 signs * (11 LDP + 3 F1C of 27 (F1 11, column 15, ADD of c6) + 8 corrections + 28 sums + 3 "
      "leaves * (index ADD, CMP, B.HS, CMN, B.HS + emit 7)) + TBZ per two incidences",
      KEY_INC == 2 * (11 + 3 * (F1_COL + COL4 + 1) + 8 + 28 + 3 * (5 + EMIT)) + Fraction(1, 2))
check("j = 1 and j = 2 paths (B.HS taken; CMN, B.HS; ADDS, ADC, MVN of N'; emit; meter ADD; B back) exceed the j = 0 path "
      "(ADD, CMP, B.HS, CMN, B.HS; emit) by 5", (3 + 2 + 3 + EMIT + 2) - (5 + EMIT) == JPATH)
check("emit: AND or BIC (8*d1), LDADD of its LC entry (the slot), d2 and w3 (LSR or EOR, and ADD of the entry pointer "
      "shifted left by 31), two STPs, TBZ of bit 13 for slot k", EMIT == 1 + 1 + 2 + 2 + 1)
DONE = 2   # per message: the port's in-flight count, up when the link takes it and down when it is written
DLV = 6    # per one-record message: its delivery, route.s dlv (three LDPs of the message, two STPs at its address, RET)
DLVK = 4 * K + 2   # per k-record message: route.s dlvk (LDP of the header, per record two LDPs and two STPs, RET)
MSG_LINK = DONE + DLV   # link-side work of one one-record message at one hop, beside the sender's own instructions
BLOCK = FLUSH + DONE + DLVK   # a full block at either hop: sent, its flush, in-flight count and delivery (dropped: 6, 8)
FINAL = FIN_DEST + FIN_SENT + DONE + DLVK + MSG_LINK   # fin per destination and hop, with both messages' link work
L8 = Fraction(LOOP8, 8)   # main-loop control per record, an upper bound: a slot of n records runs floor(n/8) passes
HOP1, HOP2 = EMIT + Fraction(BLOCK, K), H2S_REC + L8 + Fraction(BLOCK, K)
GROUP = BKT_REC + L8 + PASS4_REC + L8 + SCAN_REC
check("a full block: 1036 = flush 8 + in-flight count 2 + dlvk 1026 when sent, more than 7 or 9 when dropped; fin "
      "per destination at most 1053 (14 + 3, the final block's 1028 and the final cursor's 8)",
      BLOCK == 1036 and BLOCK > max(KEY_BLOCK_DROP, H2S_BLOCK_DROP) and FINAL == 1053
      and FINAL > FIN_DEST + FIN_HELD + MSG_LINK)
check("per record: hop 1 7 + 1036/256, hop 2 7.25 + 1036/256 (a block of k = 256 records), grouping 14.546875; "
      "36.890625 in all", HOP1 + HOP2 + GROUP == Fraction(4722, 128))

# priced allowances (analytic, per item)
GEN_MUL = 76839           # one complete scalar multiplication [k]P or [k]Q with conversion, the point's table entry
                           # (for Y: y^2, y, b, b^3, 3b^2, b^4 and two correction words; for X: its header rows and
                           # values) and its write, itemized below (as in 7e8d13c7 and ecc3ec1b)

PRECOMP = 52309152         # comb tables for P and Q at each node (2 bases, points and doubles)
COIN_SCALAR = 64           # at most two draws of 4 RNDR words (18 with one charged retry each) plus a range test (7): 50
BATCH_SETUP = 160          # copy X's header into H and its correction words into K, y_X^2 registers, acc = 1,
                           # pointers, sp, index word, both meters' updates and tests, loop (itemized below)
INV_CALL = 16              # store acc, pointer set-up, BL, load R
NODE_INIT = 2 ** 20        # program <= 2^16 bytes = 2^11 messages over 2 hops at 40 (2^17.3), plus node state
HEADER = 29                # keys' block K: G0..G5 (MOVZ, 3 MOVK each), LSR of G3 and 4 STPs
MSG = 40                   # a message outside the two hops (program, Q, phase 4), per hop, with its in-flight count
                           # and delivery (a broadcast message: 16 in all)
CAND_CAP = 2 ** 23         # candidate pairs resolved per home node
RESOLVE = 235411           # per candidate at its itemization (comb [x]P = Q against the node P table)
                           # record's three possible leaves, nine unit tests, one solve, [x]P = Q: twice the 778,872
                           # itemized, plus 1 (from 1fbd26e8)
REQ_CAP = 2 ** 29          # requests answered per source node
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 34        # relayed records held per node in each phase-4 sub-phase
EXCH = 1202590842880      # one phase-4 count exchange per node (4 per run): the itemization
BARRIERS, BARRIER = 32, 2 ** 82   # 2^10 per node and barrier
OUTPUT = 2 ** 20
DELTA_J = 2 ** 38          # meter slack for the j = 1, 2 key paths (proof.md 3.2)
DELTA_F = 2 ** 45          # cap on a node's arithmetic fixups (proof.md 3.4)
FIX_BATCH = sum(FIX_SITES.values()) * Q2   # most fixups one batch can take: every site of every incidence
CAP_P = 2 ** 22 - 2 ** 12  # records per pair slot, at hop 1 and at hop 2

# ---------------------------------------------------------------- 2. success
P1 = Fraction(6 * (r - 1) + 1, r * r)    # P(B in class(A)) for independent uniform A, B
TWO_BOUND = Fraction(2, r - 1)           # x-coordinate of a point uniform on at least r - 1 points: each value <= 2/(r-1)

# batch arithmetic sites (proof.md 3.4). Only row shifts (SHX), folds (FDX) and the square folds (F5) can drop a
# carry. Given X != O and every other generator, the fresh Y_l puts each operand value of a site at probability at most
# MULT/(r - 1), MULT being the number of points Y_l giving that value (two per x_Y, three per y_Y, and the degree of the
# map from x_Y or y_Y); a site fails with probability <= BAD * MULT / (r - 1).
SHX_BAD = 2 ** 64 * 2 ** 33     # rows B with B mod 2^192 > 2^192 - 2^33: needed for a row-shift carry
FOLD_BAD = 2 ** 66 * C0         # values in [C0, 2^66 C0): needed for the carry out of a fold of a sum below 2^322
CP_BAD = (2 ** 66 + 1) * C0     # values in [C0, (2^66 + 2) C0): the same for a sum below 2^322 + 2^256 (E, F, G)
SQ_BAD = C0 * C0 + 2 * C0       # squares in [C0, C0^2 + 2 C0): needed for the carry of a square fold (fifth word <= C0 + 1)
SITES = [  # (site, executions per batch of Q2 incidences, bad operand values, points of Y_l per value)
    ("fwd: row shifts of acc_(l-1), l >= 1", 3 * (Q2 - 1), SHX_BAD, 6),
    ("fwd: fold of acc_l, l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: row shifts of R_l", 3 * Q2, SHX_BAD, 6),
    ("bwd: fold of inv_l = R_l acc_(l-1), l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: fold of R_(l-1) = R_l D_l, l >= 1", Q2 - 1, FOLD_BAD, 6),
    ("bwd: row shifts of inv_l", 3 * Q2, SHX_BAD, 6),
    ("bwd: folds of u = inv_l (y_Y - y_X) and v", 2 * Q2, FOLD_BAD, 3),
    ("bwd: square folds of u^2 and v^2", 2 * Q2, SQ_BAD, 6),
    ("bwd: folds of E, F, G (tabled rows; b^4, a^4 added to E, F)", 3 * Q2, CP_BAD, 8),
    ("bwd: row shifts of E, F, G", 9 * Q2, SHX_BAD, 8),
    ("bwd: folds of the six components", 6 * Q2, FOLD_BAD, 8),
]
SITES_PER_INCIDENCE = dict(SHX=18, FDX=14, F5=2)   # per incidence in the fwd and bwd loops (experiment.py counts them)
EXACT_AT_L0 = dict(SHX=3, FDX=3)   # rows of acc_(-1) = 1, folds of D_0 * 1, R_0 * 1 and R_0 D_0 = 1 (never carry)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)
# arithmetic fixups (proof.md 3.4): a site takes its fixup only if its exact operand, a row it shifts, or its exact
# output, as a word below 2^256, lies in a set FIX_BAD of words, each congruent to one value mod p; with the same fresh
# generator and points per value as the site table, the mean number of fixups of a batch without a defect is E_FIX
FIX_ROW = C0 * 2 ** 192         # rows whose low word exceeds 2^64 - C0: needed for SHX's carry out of its first ADDS
FIX_FOLD = (2 ** 66 + 1) * 2 ** 128 * C0   # outputs with low 128 bits below (2^66 + 1) C0: FDX's carry into limb 2
                                           # (a sum below 2^322, or 2^322 + 2^256 for E and F)
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
    ("bwd: folds of E, F, G", "FDX", 3 * Q2, FIX_FOLD, 8),
    ("bwd: row shifts of E, F, G", "SHX", 9 * Q2, FIX_ROW, 8),
    ("bwd: folds of the six components", "FDX", 6 * Q2, FIX_FOLD, 8),
]
E_FIX = sum(Fraction(n * bad * mult, r - 1) for _, _, n, bad, mult in FIXUPS)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (2 ** 18 + 1), (r - 1) * 2 ** 18)  # P(fixed label bits of a leaf
                              # take given values): three paths, two points per x, the lattice count within 1 + 2^-18
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
    d_wj = TWO_BOUND * 23 * 2 ** 192   # KEY's index may differ from argmin: the window set W, (W1) 4, (W2) 5, (W3) 4,
                                       # (W4) 10, in units of 2^192 leaf values (proof.md 3.3)
    var_j = (Fraction(2, 9) + 2 * d_wj) * n0 + 1
    d_j = var_j / (DELTA_J - n0 // q1 - 1 - 2 * n0 * d_wj) ** 2 + q1 * d_batch
    d_f = q1 * E_FIX / (DELTA_F - FIX_BATCH) + q1 * d_batch   # fixup meter (Markov), or a defective batch of the node
    d_lab = TWO_BOUND * 783 * 2 ** 228   # the label-error set D' (proof.md 3.3): residues of T(v_0), T(v_1), T(v_2)
                                         # mod 2^28 within 73 + 73, 206 + 73 and 146 + 212 of a block boundary
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
    kb_m = KB(KEYBITS + D1 + D2)   # 184 bits: d1, d2, K2 and d4
    lam_int = RHO * (S + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
    sip = RHO * Fraction(5, 2) * (1 + 11 * (d_wj + d_lab)) * kb_m   # per leaf: the pairs of two leaves of one incidence
                                   # (5 per leaf), the second within W or D' at most 12 labels (6.4), no home condition
    other_h = RHO * (Fraction(S) * S * KB(D1 + D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) + sip * S
    d_cl = other_h / (CAND_CAP - 4)
    other_g = RHO * Fraction(n0) * S * kb_m * 2 + sip * n0
    d_rq = other_g / (REQ_CAP - 8)
    other_v = RHO * (2 + 2 * 8) * (Fraction(S, SIDE) * S * KB(D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) \
        + 18 * sip * S  # a reply counts <= 8
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_wj + d_j + d_f + d_lab + 2 * d_p + d_3 + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)  # RNDR is the rules' ideal coin source: its failure branch is charged, never taken
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_wj=d_wj, d_j=d_j, d_f=d_f, d_lab=d_lab,
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
FA, FB = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16), int("e4437ed6010e88286f547fa90abfe4c3", 16)
BOX = (p >> 100) + 1                      # a label component's value spans at most this many consecutive v (or u)
WIDE = 359 * (p >> 128) + 4               # a side widened by the largest error range (N': 359 units of T, proof.md 3.3)
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p, parallelogram 2A - B wide "
      "in u and A in v; box factor below 1 + 2^-18, one side's below 1 + 2^-19",
      FA * FA - FA * FB + FB * FB == p and (FB * BETA_P - (FA - FB)) % p == 0 and ((FA - FB) * BETA_P + FA) % p == 0
      and (BOX + WIDE + FA) * (BOX + WIDE + 2 * FA - FB) * 2 ** 200 < p * p * (1 + Fraction(1, 2 ** 18))
      and (BOX + WIDE) * 2 ** 100 < p * (1 + Fraction(1, 2 ** 19)))
check("q1 = %d is the least with success bound >= 0.390000001" % Q1, T["lo"] >= THRESH and success_terms(Q1 - 1)["lo"] < THRESH)
check("leaves S = 2^72*6*q1*q2 = 2^%.6f" % lg(S), S == NODES * 6 * Q1 * Q2)
check("mu = |Gamma|*P1 = %.9f (cross-node pairs)" % float(T["mu"]), T["mu"] > Fraction(494, 1000))
check("Chen-Stein b1+b2 < 2^-95", T["cs"] < Fraction(1, 2 ** 95), "2^%.2f" % lg(T["cs"]))
check("intruder rate per colliding pair < 2^-12.6", T["lam_int"] < pow2_interval(Fraction(-126, 10))[0], "2^%.3f" % lg(T["lam_int"]))
U3, UC3 = ((0, 1), (0, 68), (-71, 0)), ((-1, 1), (-1, 68), (-71, 2))   # T - U_m of a component; corrected (3.3)


def isum(*rs):
    """Interval sum. rs: (lo, hi) pairs. returns: (lo, hi)."""
    return sum(x[0] for x in rs), sum(x[1] for x in rs)


WR = [[isum(UC3[m], UC3[(m + t) % 3], U3[(m + 2 * t) % 3], (0, 2)) for t in range(3)] for m in (0, 1)]   # W_m(t) + kappa
NR = [(-WR[0][t][1] - WR[1][t][1] - 1, -WR[0][t][0] - WR[1][t][0] + 1) for t in range(3)]   # N' = -e0 - e1 - k + 1
LABEL_ERR = [(min(x[0] for x in rs), max(x[1] for x in rs)) for rs in (WR[0], WR[1], NR)]
check("3.3's error table from the component ranges: W0 per t [-2, 5], [-73, 71], [-72, 73]; W1 [-2, 206], [-72, 73], "
      "[-73, 71]; N' [-212, 5], [-145, 146], [-145, 146]; over t [-73, 73], [-73, 206], [-212, 146]: D' 783 residues, "
      "widening 359", WR == [[(-2, 5), (-73, 71), (-72, 73)], [(-2, 206), (-72, 73), (-73, 71)]]
      and NR == [(-212, 5), (-145, 146), (-145, 146)] and LABEL_ERR == [(-73, 73), (-73, 206), (-212, 146)]
      and sum(hi - lo for lo, hi in LABEL_ERR) == 783 and max(hi - lo for lo, hi in LABEL_ERR) + 1 == 359)
check("wrong class index rate per leaf < 2^-58.47", T["d_wj"] < pow2_interval(Fraction(-5847, 100))[0], "2^%.4f" % lg(T["d_wj"]))
for name, n, bad, mult in SITES:
    k = -(-n // Q2)
    print("  site %-52s %-6s bad <= 2^%.2f, x%d/(r-1), each <= 2^%.2f" % (name, ("%dq2" % k if k > 1 else "q2") + (
        "%+d" % (n - k * Q2) if n != k * Q2 else ""), up2(bad), mult, up2(Fraction(bad * mult, r - 1))))
check("site table: 18 row shifts, 14 folds, 2 square folds per incidence, less 3 + 3 exact at l = 0",
      sum(n for _, n, b, _ in SITES if b == SHX_BAD) == SITES_PER_INCIDENCE["SHX"] * Q2 - EXACT_AT_L0["SHX"]
      and sum(n for _, n, b, _ in SITES if b in (FOLD_BAD, CP_BAD)) == SITES_PER_INCIDENCE["FDX"] * Q2
      - EXACT_AT_L0["FDX"]
      and sum(n for _, n, b, _ in SITES if b == SQ_BAD) == SITES_PER_INCIDENCE["F5"] * Q2 and len(SITES) == 11)
check("bad sets: rows 2^97, folds below 2^98.02, square folds below 2^64.01",
      SHX_BAD == 2 ** 97 and all(b < pow2_interval(Fraction(e, 100))[0] for b, e in ((FOLD_BAD, 9802), (CP_BAD, 9802),
                                                                                     (SQ_BAD, 6401))))
for name, _, n, bad, mult in FIXUPS:
    k = -(-n // Q2)
    print("  fixup %-51s %-6s words <= 2^%.2f, x%d/(r-1), each <= 2^%.2f" % (name, ("%dq2" % k if k > 1 else "q2") + (
        "%+d" % (n - k * Q2) if n != k * Q2 else ""), up2(bad), mult, up2(Fraction(bad * mult, r - 1))))
FIX_RUNS = {m: sum(n for _, ms, n, _, _ in FIXUPS if ms == m) for m in FIX_SITES}
check("fixup table: per incidence 18 row shifts, 14 folds, 2 square folds, 2 FS and 1 FADDW, less 3 + 3 exact at l = 0",
      FIX_RUNS == {m: FIX_SITES[m] * Q2 - EXACT_AT_L0.get(m, 0) for m in FIX_SITES}
      and all(bad == dict(SHX=FIX_ROW, FDX=FIX_FOLD, F5=FIX_SQ, FS=FIX_SUB, FADDW=FIX_SUB)[m] for _, m, _, bad, _ in FIXUPS)
      and all(FIX_SITES[m] == SITES_PER_INCIDENCE[m] for m in ("SHX", "FDX", "F5")))
check("a call's fixup word (C0 per fixup, at most 37 q2 fixups) stays below 2^63", FIX_BATCH * C0 < 2 ** 63,
      "2^%.3f" % lg(FIX_BATCH * C0))
check("mean fixups per node (clean batches) q1*E_FIX = 2^%.3f < 2^29.15; cap 2^45 > 2^15 times it plus 37 q2" % lg(Q1 * E_FIX),
      Q1 * E_FIX < pow2_interval(Fraction(2915, 100))[0] and DELTA_F > 2 ** 15 * Q1 * E_FIX + FIX_BATCH)
check("fixup meter rate (Markov) < 2^-15.85", T["d_f"] < pow2_interval(Fraction(-1585, 100))[0], "2^%.3f" % lg(T["d_f"]))
check("label-error rate per leaf < 2^-17.38", T["d_lab"] < pow2_interval(Fraction(-1738, 100))[0],
      "2^%.3f" % lg(T["d_lab"]))
check("batch defect rate (site table, D = 0, X = O) < 2^-125.6", T["d_batch"] < pow2_interval(Fraction(-1256, 10))[0],
      "2^%.3f" % lg(T["d_batch"]))
check("j = 1, 2 meter rate < 2^-23.4", T["d_j"] < pow2_interval(Fraction(-234, 10))[0], "2^%.3f" % lg(T["d_j"]))
check("pair-slot overflow rate, each hop, < 2^-19.9", T["d_p"] < pow2_interval(Fraction(-199, 10))[0],
      "2^%.3f" % lg(T["d_p"]))
check("d3 bucket overflow rate < 2^-30", T["d_3"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_3"]))
check("candidate cap rate < 2^-18.4, request cap rate < 2^-25.9",
      T["d_cl"] < pow2_interval(Fraction(-184, 10))[0] and T["d_rq"] < pow2_interval(Fraction(-259, 10))[0],
      "2^%.3f 2^%.3f" % (lg(T["d_cl"]), lg(T["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-25.2", T["d_rl"] < pow2_interval(Fraction(-252, 10))[0],
      "2^%.3f" % lg(T["d_rl"]))
check("whole-run row: a computational defect in any batch < 2^-26.5", T["g_def"] < pow2_interval(Fraction(-265, 10))[0],
      "2^%.3f" % lg(T["g_def"]))
check("whole-run rows (defective batch anywhere, coin rejection twice) < 2^-26.5",
      T["glob"] < pow2_interval(Fraction(-265, 10))[0], "2^%.3f" % lg(T["glob"]))
check("success >= 1 - e^-mu - (b1+b2) - rows - global >= 0.390000001", T["lo"] >= THRESH,
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
row("Y generators [c]P with their table entries: q2 * 76839 per node", NODES * Q2 * GEN_MUL, True)
row("R_g = [b]Q and X = R_g + [a]P with their header values: (q1+1) * 76839 per node", NODES * (Q1 + 1) * GEN_MUL)
row("comb tables for P and Q: 52309152 per node", NODES * PRECOMP, True)
row("keys' block constants: 29 per node", NODES * HEADER, True)
row("LC tables (init3) and RC tables (tinit): 2 * (4 * 2^36 + 1) + 2 * (3 * 2^36 + 1) per node",
    NODES * 2 * (INIT3_BKT * SIDE + 1 + TINIT_ENTRY * SIDE + 1), True)
row("Q broadcast: 2 messages to every node over two hops at 40", 2 * (SIDE + NODES) * MSG)
row("d4 counters zeroed: 2^46/4 * 4 + 1 per node", NODES * (2 ** D4 // 4 * ZERO4_PER4 + 1), True)
row("d3 cursor table: 4 * 2^24 + 1 per node", NODES * (INIT3_BKT * 2 ** D3 + 1), True)
row("forward: q1 * (97.5 q2 + 5)", int(NODES * Q1 * (FWD_INC * Q2 + FWD_FIX)))
row("inversion: q1 * (23,646 + 16)", NODES * Q1 * (INV + INV_CALL))
row("backward: q1 * (1206.5 q2 + 5)", int(NODES * Q1 * (BWD_INC * Q2 + BWD_FIX)))
row("arithmetic fixups: 2^45 * 4 per node", NODES * DELTA_F * max(FIX_EXTRA.values()))
row("key pass without emits: q1 * (286.5 q2 + 11)", int(NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX)))
row("j = 1 and j = 2 key paths: cap_j * 5 per node", NODES * CAP_J * JPATH)
row("batch set-up: q1 * 160", NODES * Q1 * BATCH_SETUP)
row("hop 1, the emits and full blocks (flush, in-flight count, dlvk): S * (7 + 1036/256)", S * (EMIT * K + BLOCK) // K)
row("hop 2, h2s and full blocks: S * (7.25 + 1036/256) + 11 * 2^36 + 4 per node",
    S * (H2S_REC * K + LOOP8 * K // 8 + BLOCK) // K + NODES * (H2S_SLOT * SIDE + H2S_FIX))
row("final blocks and cursors, in-flight counts and deliveries, both hops: 2 * (1053 * 2^36 + 2) per node",
    NODES * 2 * (FINAL * SIDE + FIN_FIX))
row("bucket pass bkt: S * 8.25 + 12 * 2^36 + 4 per node", S * BKT_REC + S * LOOP8 // 8 + NODES * (BKT_SLOT * SIDE + BKT_FIX))
row("prefix sums of d4: 2^46/4 * 10 + 1 per node", NODES * (2 ** D4 // 4 * PREFIX4_PER4 + 1))
row("local pass pass4: S * 3.25 + 11 * 2^24 + 2 per node",
    S * PASS4_REC + S * LOOP8 // 8 + NODES * (PASS4_BKT * 2 ** D3 + PASS4_FIX))
row("scan_ranges: S * 195/64 + 25 * 2^46 + 2 + 8 * (2^23 + 1) per node", S * SCAN_REC.numerator // SCAN_REC.denominator
    + NODES * (SCAN_FIX + 8 * (CAND_CAP + 1)))
row("candidates: 2^23 * 235411 per node", NODES * CAND_CAP * RESOLVE)
row("requests served: 2^29 * 2^10 per node", NODES * REQ_CAP * REQ_SERVE)
row("phase-4 count exchanges: 4 * 1202590842880 per node", NODES * 4 * EXCH)
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
    """Smallest multiple of 10^-15 whose power of two is at least w, verified by interval bounds.

    w: Fraction > 1.
    returns: Fraction claim.
    """
    c = Fraction(math.ceil(lg(w) * 10 ** 15), 10 ** 15)
    for _ in range(10000):
        if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 15))[1] < w:
            return c
        c += Fraction(1, 10 ** 15) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 15)
    raise AssertionError("claim not certified")


CLAIM = claim_of(W)
check("claim = least 15-decimal c with 2^c >= R/1000", pow2_interval(CLAIM)[0] >= W
      and pow2_interval(CLAIM - Fraction(1, 10 ** 15))[1] < W, "%d/10^15" % (CLAIM * 10 ** 15))
PER_LEAF = Fraction(R_TOTAL, S)
check("work per leaf = %.7f instructions" % float(PER_LEAF), 305 < PER_LEAF < 307)

# ---------------------------------------------------------------- 4. memory (bytes per node, layout per phase)
CAP_3 = T["cap_3"]
SLOTS = SIDE * CAP_P * 32     # 2^36 pair slots of cap_p records: R1 (hop-1 receive) or R2 (hop-2 receive)
FIXED = [("code, constants, stack", 2 ** 20), ("barrier words: 5 granules of at most 2 KiB", 5 * 2048),
         ("header H and keys' block K", 1024), ("X table: H's rows and values, -T(a), -T(beta a), y^2 per X", Q1 * 512),
         ("Y table: y^2, y, b, b^3, 3b^2, b^4", Q2 * 192), ("scalars a, b, c", GENS * 32), ("batch entries", Q2 * 256), ("comb tables for P and Q", 4 * 32 * 256 * 64),
         ("the ports' message blocks", 2 * 64), ("LC tables at 0 and 2^39", 2 * SIDE * 8),
         ("RC tables at 2^43 and 2^43 + 2^39", 2 * SIDE * 8), ("record blocks at 2^50 and 2^51",
           2 * SIDE * 2 ** 14),
         ("final-cursor tables E1 and E2", 2 * SIDE * 32), ("d3 cursor table at 2^40", 2 ** D3 * 8),
         ("d4 counters at 2^49", 2 ** D4 * 8), ("candidate list", CAND_CAP * 16),
         ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64), ("relay buffer", RELAY_CAP * 64),
         ("phase-4 count, offset and cursor tables", 3 * SIDE * 32)]
FIXED_BYTES = sum(b for _, b in FIXED)
BUCKET_STRIDE = 2 ** 38      # bucket t's cap_3 entries of 16 bytes end at R1 + (t + 1)*2^38 (bkt's test of 38 low bits)
PHASE3 = [("d3 buckets: 2^24 strides of 2^38 bytes, cap_3 entries (address, d4 counter) at the end of each",
           2 ** D3 * BUCKET_STRIDE),
          ("pointer array", 2 ** D3 * CAP_3 * 8)]
PHASE3_BYTES = sum(b for _, b in PHASE3)
for name, b in FIXED + [("R1 at 2^52: hop-1 pair slots (phases 1-2)", SLOTS),
                        ("R2 at 2^52 + |R1|: hop-2 pair slots (phases 2-4)", SLOTS)] + PHASE3:
    print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("fixed regions lie below 2^52", FIXED_BYTES <= 2 ** 52, "2^%.4f" % lg(FIXED_BYTES))
ALIGNED = 4 * SIDE * 8 + 2 ** D3 * 8 + 2 ** D4 * 8 + 2 * SIDE * 2 ** 14
check("aligned tables: LC [0, 2^39), [2^39, 2^40), d3 [2^40, 2^40 + 2^27), RC [2^43, 2^44), d4 [2^49, 2^50), "
      "column blocks [2^50, 2^51), row blocks [2^51, 2^51 + 2^50); the rest fits between 2^44 and 2^48",
      SIDE * 8 == 2 ** 39 and 2 ** D3 * 8 == 2 ** 27 and 2 ** D4 * 8 == 2 ** 49 and SIDE * 2 ** 14 == 2 ** 50
      and D4BASE == 1 << 49 and D4BASE + 2 ** D4 * 8 == PB1C and PB1C == 1 << 50 and PB2C == 1 << 51
      and PB2C + (1 << 50) <= 1 << 52 and FIXED_BYTES - ALIGNED <= (1 << 48) - (1 << 44),
      "2^%.4f" % lg(FIXED_BYTES - ALIGNED))
check("record block: header and k slots within 2^14; bit 13 set exactly at slot k (TBZ #13); bases/stride keep "
      "low 13 bits zero; cap_p a multiple of k (keep lemma)",
      32 + 32 * K <= 2 ** 14 and all((((32 * n) >> 13) & 1 == 1) == (n == K) for n in range(1, K + 1))
      and PB1C & 0x1FFF == 0 and PB2C & 0x1FFF == 0 and BLKC & 0x1FFF == 0 and CAP_P % K == 0
      and S % K == 0)
check("bwd's and keys' loop end: the q2 256-byte batch entries fill the lower 2^33 bytes of their 2^34-aligned block",
      Q2 * 256 == 2 ** 33)
check("record fields: k < 2^28 in w2 bits 36-63; the sign's code 2 or 5 in w3 bits 36-38 and l < 2^25 (l - 1 mod 2^25 "
      "for s = -1) in bits 39-63", Q1 <= 2 ** 28 and Q2 <= 2 ** 25)
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
# least multiple of 10^-12 or 10^-4 at or above the exact value), the success figure rounded down
DECLARED = dict(time_log2=Fraction(124988192390288425, 10 ** 15), memory_log2_bytes=Fraction(1359989, 10 ** 4),
                preprocessing_log2=Fraction(1081053, 10 ** 4), success=Fraction(390000001, 10 ** 9),
                per_leaf=Fraction(3053866, 10 ** 4))


def least_at_or_above(d, x, scale):
    """Whether 2^d >= x and 2^(d - 1/scale) < x, both certified by the interval bounds of pow2_interval.

    d: Fraction; x: positive int or Fraction; scale: denominator of the decimal grid.
    returns: bool.
    """
    return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x


check("declared time_log2 %d/10^15: the claim, least 15-decimal value with 2^c >= R/1000" % (DECLARED["time_log2"] * 10 ** 15),
      DECLARED["time_log2"] == CLAIM and least_at_or_above(DECLARED["time_log2"], W, 10 ** 15))
check("declared memory 2^135.9989 bytes: least 4-decimal bound >= peak (2^%.7f)" % lg(PEAK),
      least_at_or_above(DECLARED["memory_log2_bytes"], PEAK, 10 ** 4))
check("declared preprocessing 2^%s units: least 4-decimal bound >= rows marked P / 1000" % format(float(DECLARED["preprocessing_log2"]), ".4f"),
      least_at_or_above(DECLARED["preprocessing_log2"], Fraction(R_PRE, C), 10 ** 4))
check("declared success 0.390000001 <= the bound", DECLARED["success"] <= T["lo"])
check("declared work per leaf 305.3866: least 4-decimal bound >= R/S (%.7f)" % float(PER_LEAF),
      DECLARED["per_leaf"] >= PER_LEAF > DECLARED["per_leaf"] - Fraction(1, 10 ** 4))
ANALYTIC = sum(v for n, v, _ in rows if n.split(":")[0] in ("program delivery, node and port state", "coins",
               "Y generators [c]P with their table entries", "R_g = [b]Q and X = R_g + [a]P with their header values", "comb tables for P and Q",
               "keys' block constants", "batch set-up",
               "Q broadcast", "candidates", "requests served", "phase-4 count exchanges", "barriers", "phase control",
               "output",
               "inversion")) - NODES * Q1 * INV
# itemized work behind the allowances (proof.md section 7)
# on the exact arithmetic of proof.md 3.4: fe_mul_f and fe_sqr_f as executed (115 and 82) plus their calls (pointers and
# BL: 4 and 3, as in fe_inv), and FS or FADDW (9) with both operands brought below p (2 x 9: ADDS chain of C0, CSELs)
MULX, SQRX, ADDX = 115 + 4, 82 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
# fixed-base comb, w = 8: 32 byte-adds against 32 sub-tables of 256 affine points; per add: address (2), x/y
# loads (4 LDPs), mixed add (1,477) with Z/H below-p and zero tests, control (6), path (16); one conversion; 64 set-up
COMB_ADD = 2 + 4 + MADD_ + 2 * (9 + 4) + 6 + 16
GEN_ITEM = 32 * COMB_ADD + (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_TEST = 0   # tests and paths are inside COMB_ADD
GEN_CTRL = 0   # byte extract and control are inside COMB_ADD
NEGT = 30   # a correction word -T(v) mod 2^128 of a canonical v (proof.md 3.3): T(v) = (v >> 128) + e, e = 1 iff
            # 2^128 (v mod 2^128) + (v >> 128) C0 >= p (two MULs and UMULHs by C0, an ADDS chain and the carry of adding
            # C0 to it), then e added and the two words negated
RED = 9     # a value brought below p before it is stored or keyed (an ADDS chain of C0 and CSELs, as in ADDX)
Y_ENTRY = 3 * SQRX + 2 * MULX + 2 * ADDX + 2 * NEGT + 7 * RED + 14   # y^2, b^2, b^4 = (b^2)^2; b^3, beta b; 3b^2 = b^2 +
                                            # b^2 + b^2; -T(b), -T(beta b) into the batch entry; 7 reductions; 14 STPs
X_ENTRY = MADD_ + 3 * SQRX + 11 * MULX + 2 * ADDX + 2 * NEGT + 17 * RED + 32   # X = R_g + [a]P; a^2, a^4 = (a^2)^2,
                                            # y^2; a^3, beta a and the rows c 2^(64k) mod p, k = 1..3, of c = 2a^3, 2a and
                                            # a^2 (nine products by constants); 2a^3, 2a; -T(a), -T(beta a); 17 reductions
                                            # (the 12 rows, y_X, y^2, a^4, a, beta a); 32 STPs
check("generator (exact arithmetic): doubling 933, mixed addition 1,477; with the exceptional-case tests, control, "
      "the Y entry and the X header values at most %d, the allowance" % (GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY,
      X_ENTRY)), DBL == 933 and MADD_ == 1477 and GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY, X_ENTRY) <= GEN_MUL,
      "%.4fx" % (math.floor(10 ** 4 * GEN_MUL / (GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY, X_ENTRY))) / 10 ** 4))
# a barrier at one node: at most 5 waits (4 children, the release), each LDAXR, CMP, B.EQ and at most one WFE, B and
# re-check (8); at most 5 messages (one up, four down), each 3 STPs and 3 for its destination word and address, plus its in-flight
# count and delivery; DSB; 20 set-up
BARRIER_ITEM = 5 * 8 + 5 * (6 + MSG_LINK) + 1 + 20
# batch set-up at one node, with nothing kept in registers across calls: the X-table pointer and k from the driver's
# frame (2 LDRs) and H's address (4 MOVZ/MOVK); X's header (rows of 2a^3, 2a, a^2, y_X, a^4: 28 LDPs, 28 STPs) into H
# and its correction words into K (2 LDPs, 2 STPs); fwd's arguments (q2, two pointers of at most 4 MOVZ/MOVK each,
# acc = 1 by 4 MOVs, y_X^2 by 2 LDPs) and BL; bwd's two pointers and BL (R is loaded by the inversion call); keys' entry
# pointer, slot end (LDR), index word (LDR, ADD), key-path count (LDR) and BL, and the count's STR after it; sp set to H
# for bwd and to K for keys and back to the driver's (4 MOVZ/MOVK and a MOV each); after fwd and after bwd, C0 (MOVZ,
# MOVK), the fixup word over C0 (UDIV) added to the node's fixup count (LDR, ADD, STR); before the batch, each meter's
# test (LDR of the count, its cap by 4 MOVZ/MOVK, CMP, B.HI); the loop (ADD and STR of k, q1 by MOVZ/MOVK, CMP, B.LO)
# and the X-table pointer step (ADD, STR)
BATCH_ITEM = (2 + 4) + 56 + 4 + (1 + 4 + 4 + 4 + 2 + 1) + (2 * 4 + 1) + (4 + 1 + 2 + 1 + 1 + 1) + 4 * (4 + 1) \
    + 2 * (2 + 1 + 3) + 2 * (1 + 4 + 1 + 1) + (2 + 2 + 1 + 1) + 2
check("batch set-up: at most %d instructions, within its allowance of 160" % BATCH_ITEM, BATCH_ITEM <= BATCH_SETUP)
check("barrier: at most 131 instructions per node, in-flight counts and deliveries included (2^10 per node is > 7x)",
      BARRIER_ITEM == 131 and BARRIER > 7 * NODES * BARRIER_ITEM)
# a phase-4 count exchange at one node: at most 12 per entry of its 2^36-entry received-count table (LDP, epoch CMP,
# CSEL, prefix ADD, zero test, cap CMP and branch, loop), and at most 88 per message (its count, the count message,
# the offset reply and the write, each at most 22 with its in-flight count and delivery) for at most 8 * 2^29 = 2^32
# messages (2^29 requests of up to 8 records)
EXCH_ITEM = 12 * SIDE + 4 * (14 + MSG_LINK) * (8 * REQ_CAP)
check("phase-4 count exchange: %d per node, the allowance (at most 8 * 2^29 messages)" % EXCH_ITEM,
      EXCH_ITEM == EXCH, "2^%.3f" % lg(EXCH_ITEM))
# a broadcast message (program, Q): two LDPs of its words, three STPs, ADD of 8 to its destination word, SUBS, B.NE,
# its in-flight count and delivery
check("a broadcast message costs at most 16, in-flight count and delivery included, within its price of 40",
      2 + 3 + 3 + MSG_LINK <= MSG)
# a record's three candidate leaves phi^t(X) + sY (t = 0, 1, 2, as w3 gives s but not t): D = x_Y^3 - x_X^3 (two
# squares, two products, a subtraction) and its inverse; Q_t = x_Y^2 + beta^t x_X x_Y + beta^(2t) x_X^2 (x_X x_Y and four
# products by beta, beta^2, six additions); 1/(x_Y - beta^t x_X) = Q_t/D (three products); beta x_X, beta^2 x_X; per t
# the slope, its square, x and y (two products, a square, four additions), and s y_Y - y_X once
PRECOMP_ITEM = 256 * DBL + 8160 * (MADD_ + 2 * (9 + 4)) + 32 * DBL + 8160 * DBL + 3 * 16351 * MULX + INV + 16351 * 9
check("comb precompute per base, twice that the allowance", 2 * PRECOMP_ITEM <= PRECOMP)
CAND_LEAVES = INV + 5 * SQRX + 18 * MULX + 20 * ADDX
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 9 * 600 + 384 * 200 + (GEN_ITEM + GEN_TEST + GEN_CTRL) + INV
check("candidate: comb [x]P = Q = %d; charged at its itemization" % CAND_ITEM, RESOLVE == CAND_ITEM)
PROGRAM_BYTES = 2 ** 16
check("program delivery: 2^16 bytes = 2^11 messages over 2 hops at 40 = 2^17.3 < node set-up 2^20",
      PROGRAM_BYTES // 32 * 2 * MSG < NODE_INIT)
check("analytic allowances are below 2^-13.7 of the total", ANALYTIC * pow2_interval(Fraction(137, 10))[1] < R_TOTAL,
      "2^%.3f" % lg(Fraction(ANALYTIC, R_TOTAL)))

print()
print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1, N0, lg(S), float(T["mu"]),
                                                                    math.floor(T["lo"] * 10 ** 10) / 10 ** 10))
print("R_total %d instructions  W 2^%.13f  claim %d/10^15  per leaf %.4f" % (R_TOTAL, lg(W), CLAIM * 10 ** 15,
                                                                      float(PER_LEAF)))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
    format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK), format(float(DECLARED["preprocessing_log2"]), ".4f"),
    lg(Fraction(R_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
    sys.exit(1 if FAILED else 0)
