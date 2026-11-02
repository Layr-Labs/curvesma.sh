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
D1, D2, M_CELL = 36, 36, 57   # the digits; a record's cell is K3's low 57 bits (proof.md 2.2, 4)
KEYBITS = 64 + M_CELL         # label bits a candidate matches beyond d1, d2: K2 (64) and its cell (57)
THRESH = Fraction(390000001, 10 ** 9)   # target 0.390000001, 10^-9 above the rules' 0.39

# executed counts (experiment.py parts A and B re-derive each of these by running the kernels)
LEAVES = 6                        # leaves per incidence: x(phi^t(X) + Y) and x(phi^t(X) - Y), t = 0, 1, 2
FWD_INC, FWD_FIX = 98, 5          # leaf.s fwd per incidence (two per pass: SUBS #2, B.NE); mov, movk, the fixup word's
                                  # STR and LDR, ret
BWD_INC, BWD_FIX = Fraction(2595, 2), 5   # leaf.s bwd per incidence (twelve fractions; TBZ per two); mov, movk, STP and
                                  # LDP of x30 and the fixup word, ret
FIX_SITES = dict(SHX=18, FDX=10, F5=2, FS=2, FADDW=1)   # fixup sites per incidence in fwd and bwd, SHX with the
                                                        # three SHXT (proof.md 3.4)
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4)     # a fixup's instructions: the rest of the chain and STADD
EMIT = 8                          # keys' emit of a leaf into its destination's record block (proof.md 3.2)
EMIT_INC = LEAVES * EMIT          # 48 per incidence
KEY_INC, KEY_FIX = Fraction(293, 2), 5   # leaf.s keys per incidence (six leaves on the j = 0 path, TBZ per two) incl.
                                  # the emits; set-up (the mask, x4 = 32, x0 = 32k, MOV of 2^36 - 1) and ret
K = 128                           # records per block message (proof.md 2.1)
FLUSH = 7                         # after slot k: SUB, STR (LC back to slot 1), ORR (RC entry), LDADD of 32k, CMP, B.HS or
                                  # B.LO, STP of the header: a block sent, in keys and in h2s
KEY_BLOCK_DROP, H2S_BLOCK_DROP = 6, 8   # a block dropped: without the STP; in h2s the end stored back (STR) and B
JPATH = 5                         # extra on a leaf whose class index j is 1 or 2: the two-word sum and its complement, the
                                  # meter's ADD and the B back, beyond the j = 0 path's
INV = 23646                       # secp_fe.s fe_inv, executed
TINIT_ENTRY = 3                   # route.s tinit per RC-table entry
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 14, 3, 2, 2   # route.s fin per destination: LDR of LC and RC, SUBS, B.EQ,
                                  # ADD, CMP, CSEL, three STPs, two ADDs, SUBS, B.NE; a final block sent CMP, B.HS, STP
                                  # (held: CMP, B.HS); per call MOV and RET
H2S_REC, H2S_SLOT = 8, 11                 # route.s h2s per record (LDADD of its LC entry returns its slot), per slot
H2S_FIX = 4                               # h2s: x10 = 2^39, x17 = 32, x18 = 32k, RET
H2G_REC, H2G_SLOT, H2G_FIX = 5, 12, 3     # route.s h2g per record (LDP, BFI of its cell, SWP into it, STR of the old
                                          # head over K3, CBZ), per slot, and per call (RET, or a stop's TBNZ, MOV, RET)
H2G_STEP, H2G_HIT = 5, 6                  # per walk step (LDP, SUBS, CCMP, B.EQ, CBZ or CBNZ), and at most 6 more per
                                          # candidate (TBNZ, CMP, B.HS, STP, B or CBNZ and B; 5 as executed)
ZERO_PASS, ZERO_BYTES = 1026, 2 ** 14     # h2zero: 1024 STPs of xzr, CMP and B.NE per 16 KiB of the head table
SLACK_W = 2 ** 40                         # the walk budget above its mean bound (proof.md 4, 6.4)
LOOP8 = 2                                 # h2s, h2g: CMP and B.LO per 8 records of a slot's main loop
INIT3_ENTRY = 4                           # init3 per LC-table entry
check("key pass: 146.5 = 2 signs * (6 LDP of the fractions + 28 sums + 3 leaves * (index ADD, CMP, B.HS, CMN, B.HS + "
      "emit 8)) + TBZ per two incidences", KEY_INC == 2 * (6 + 28 + 3 * (5 + EMIT)) + Fraction(1, 2))
check("j = 1 and j = 2 paths (B.HS taken; CMN, B.HS; ADDS, ADC, MVN of N'; emit; meter ADD; B back) exceed the j = 0 path "
      "(ADD, CMP, B.HS, CMN, B.HS; emit) by 5", (3 + 2 + 3 + EMIT + 2) - (5 + EMIT) == JPATH)
check("emit: AND or BIC (8*d1), LDADD of its LC entry (the slot), d2 and w3 (LSR or EOR, and ADD of the entry pointer "
      "shifted left by 31), two STPs, TST and B.NE of slot k", EMIT == 1 + 1 + 2 + 2 + 2)
DONE = 2   # per message: the port's in-flight count, up when the link takes it and down when it is written
DLV = 6    # per one-record message: its delivery, route.s dlv (three LDPs of the message, two STPs at its address, RET)
DLVK = 4 * K + 2   # per k-record message: route.s dlvk (LDP of the header, per record two LDPs and two STPs, RET)
MSG_LINK = DONE + DLV   # link-side work of one one-record message at one hop, beside the sender's own instructions
BLOCK = FLUSH + DONE + DLVK   # a full block at either hop: sent, its flush, in-flight count and delivery (dropped: 6, 8)
FINAL = FIN_DEST + FIN_SENT + DONE + DLVK + MSG_LINK   # fin per destination and hop, with both messages' link work
L8 = Fraction(LOOP8, 8)   # main-loop control per record, an upper bound: a slot of n records runs floor(n/8) passes
HOP1, HOP2 = EMIT + Fraction(BLOCK, K), H2S_REC + L8 + Fraction(BLOCK, K)
GROUP = H2G_REC + L8
check("a full block: 523 = flush 7 + in-flight count 2 + dlvk 514 when sent, more than 6 or 8 when dropped; fin per "
      "destination at most 541 (14 + 3, the final block's 516 and the final cursor's 8)",
      BLOCK == 523 and BLOCK > max(KEY_BLOCK_DROP, H2S_BLOCK_DROP) and FINAL == 541
      and FINAL > FIN_DEST + FIN_HELD + MSG_LINK)
check("per record: hop 1 8 + 523/128, hop 2 8.25 + 523/128 (a block of k = 128 records), grouping (h2g's inserts) "
      "5.25; 29.671875 in all", HOP1 + HOP2 + GROUP == Fraction(3798, 128))

# priced allowances (analytic, per item)
GEN_MUL = 677570           # one complete scalar multiplication [k]P or [k]Q with conversion, the point's table entry
                           # (for Y: y^2, y, b, b^3, 3b^2, b^4 and two correction words; for X: its header rows and
                           # values) and its write, itemized below (as in 7e8d13c7 and ecc3ec1b)
COIN_SCALAR = 64           # at most two draws of 4 RNDR words (18 with one charged retry each) plus a range test (7): 50
BATCH_SETUP = 160          # copy X's header with its correction words into H, y_X^2 registers, acc = 1,
                           # pointers, sp, index word, both meters' updates and tests, loop (itemized below)
INV_CALL = 16              # store acc, pointer set-up, BL, load R
NODE_INIT = 2 ** 20        # program <= 2^16 bytes = 2^11 messages over 2 hops at 40 (2^17.3), plus node state
HEADER = 72                # the rows of beta below H: 16 words (MOVZ, 3 MOVK each) and 8 STPs
MSG = 40                   # a message outside the two hops (program, Q, phase 4), per hop, with its in-flight count
                           # and delivery (a broadcast message: 16 in all)
CAND_CAP = 2 ** 20         # candidate pairs resolved per home node
RESOLVE = 1557745          # per candidate: two requests and two replies of <= 8 records over 2 hops at 40, each
                           # record's three possible leaves, nine unit tests, one solve, [x]P = Q: twice the 778,872
                           # itemized, plus 1 (from 1fbd26e8)
REQ_CAP = 2 ** 29          # requests answered per source node
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 34        # relayed records held per node in each phase-4 sub-phase
EXCH = 2 ** 41             # one phase-4 count exchange per node (4 per run)
BARRIERS, BARRIER = 32, 2 ** 82   # 2^10 per node and barrier
OUTPUT = 2 ** 20
DELTA_J = 2 ** 42          # meter slack for the j = 1, 2 key paths (proof.md 3.2)
DELTA_F = 2 ** 45          # cap on a node's arithmetic fixups (proof.md 3.4)
FIX_BATCH = sum(FIX_SITES.values()) * Q2   # most fixups one batch can take: every site of every incidence
CAP_P = 2 ** 22 - 2 ** 12  # records per pair slot, at hop 1 and at hop 2

# ---------------------------------------------------------------- 2. success
P1 = Fraction(6 * (r - 1) + 1, r * r)    # P(B in class(A)) for independent uniform A, B
TWO_BOUND = Fraction(2, r - 1)           # x-coordinate of a point uniform on at least r - 1 points: each value <= 2/(r-1)

# batch arithmetic sites (proof.md 3.4): only row shifts (SHX, SHXT), folds (FDX) and square folds (F5) drop a carry;
# given the other generators a site fails with probability <= BAD * MULT / (r - 1), MULT the points Y_l per value
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
    ("bwd: folds of beta u^2 and beta v^2 (tabled rows of beta)", 2 * Q2, FOLD_BAD, 6),
]
SITES_PER_INCIDENCE = dict(SHX=18, FDX=10, F5=2)   # per incidence in the fwd and bwd loops (experiment.py counts them)
EXACT_AT_L0 = dict(SHX=3, FDX=3)   # rows of acc_(-1) = 1, folds of D_0 * 1, R_0 * 1 and R_0 D_0 = 1 (never carry)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)
# arithmetic fixups (proof.md 3.4): a site takes its fixup only if its row or output lies in its set of words, each
# word one value mod p; with the site table's generators and points per value, a clean batch's mean is E_FIX
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
    ("bwd: folds of beta u^2 and beta v^2", "FDX", 2 * Q2, FIX_FOLD, 6),
]
E_FIX = sum(Fraction(n * bad * mult, r - 1) for _, _, n, bad, mult in FIXUPS)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (2 ** 17 + 1), (r - 1) * 2 ** 17)  # P(fixed label bits of a leaf
                              # take given values): three paths, two points per x, the lattice count within 1 + 2^-17
RHO = Fraction(r, r - 1) ** 4  # slack for the exclusion Y != O when leaves are uniform on r - 1 points or pairs
FIXED_LEAVES = 2 * LEAVES     # the leaves of the pair's two incidences, fixed by its four generators
PARTNER_PAIRS = 2 * 10        # pairs of two partners of one of the pair's leaves (proof.md 6.4), counted without the
                              # home or relay-column condition


FA, FB = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16), int("e4437ed6010e88286f547fa90abfe4c3", 16)
BOX = (p >> 100) + 1                      # a label component's value spans at most this many consecutive v (or u)
WIDE = 654 * (p >> 128) + 4               # a side widened by the largest error range (N': 654 units of T, proof.md 3.3)
HV, HU = BOX + WIDE + FA, BOX + WIDE + 2 * FA - FB   # an expanded box's sides: a widened block plus P's extent
DIAM = math.isqrt(max(FA * FA + FB * FB, (2 * FB - FA) ** 2 + (2 * FA - FB) ** 2)) + 1   # diameter of P
KAPPA = Fraction(HV * HU * 2 ** 200, p * p)   # an expanded box's area over (p/2^100)^2


def eps_u(a, b, edge):
    """Lemma U's excess (proof.md 3.3) for cells fixing a low bits of a label's first part and b of its second.

    a, b: digit bits; edge: whether the cell is an edge cell (W's horizontal strips reach it).
    returns: Fraction eps, with at most p 2^-(a+b) (1 + eps) values x whose widened cell set holds the cell.
    """
    f = a + b
    unit = Fraction(p, 2 ** f)
    box = Fraction(HV * HU, p)                   # lattice points of one expanded box, at most
    can = KAPPA * (1 + Fraction(7 * 2 ** b, 4 * 2 ** 100) + Fraction(6 * 2 ** a, 2 ** 100))
    edges = Fraction(3 * 2 ** (100 - b) * 4 * 3 * DIAM * HU * 2 ** f, p * p)   # a diam(P) strip along C's four edges
    top = Fraction(3 * (FA - FB) * ((p // 2 + FA) * HV * 2 ** 100 // (2 ** a * p) + HV + 1) * 2 ** f, p * p)
    rows = (1 if edge else 0) if b > 38 else (2 ** 36 // 2 ** b + 1) + ((2 ** 38 + 1) // 2 ** b + 1)
    cols = (2 ** 36 // 2 ** a + 1) + ((2 ** 38 + 1) // 2 ** a + 1)
    bands = 16 * (-((-5 * 2 ** 36) // 2 ** a) + 3)   # four slanted strips, four bands each, per row of the cell
    strips = 3 * (rows * 2 ** (100 - a) + (cols + bands) * 2 ** (100 - b)) * box / unit
    return can - 1 + edges + top + strips


EPS_GEN, EPS_EDGE, EPS_HOME = eps_u(D1, D2 + M_CELL, False), eps_u(D1, D2 + M_CELL, True), eps_u(D1, D2, False)
EDGE_LO = -(-(2 ** 93 - 2 ** 38 - 1 - (2 ** 36 - 1)) // 2 ** 36)   # least top-edge K3 mod 2^57 (some d2)
N_EDGE = (2 ** 36 // 2 ** 36 + 1) + (2 ** M_CELL - EDGE_LO)        # edge cells per home: K3 mod 2^57 in {0, 1} or top
KBU = lambda fixed, eps: Fraction(2 * 2 ** (256 - fixed), r - 1) * (1 + eps)   # Lemma U's mass bound KB'


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
    d_lab = TWO_BOUND * 1305 * 2 ** 228   # the label-error set D' (proof.md 3.3): residues of T(v_0), T(v_1), T(v_2)
                                          # mod 2^28 within 194 + 132, 194 + 132 and 264 + 389 of a block boundary
    # caps (Chebyshev): count means of the other leaves given the pair's four generators; the variance adds, for the
    # 2(q1 + q2) incidences through those generators, 2 * 15 pairs of one incidence at most one leaf's mean each
    m_p = n0 * KB(D1) + 1        # a node's leaves on d1 (hop 1), or a column's leaves on (d1, d2) (hop 2): equal
    v_p = m_p + 1 + 60 * (q1 + Q2) * KB(D1)
    d_p = v_p / Fraction(CAP_P - FIXED_LEAVES - m_p) ** 2
    # the walk budget (Chebyshev): steps at alpha's home are at most P', the pairs of leaves whose widened cell sets share
    # one of its 2^57 cells, a function of leaf values (proof.md 3.3, success.md 6.4); a leaf has at most 12 such cells
    kg, ke = KBU(D1 + D2 + M_CELL, EPS_GEN), KBU(D1 + D2 + M_CELL, EPS_EDGE)
    cellsum = 2 ** M_CELL * kg + N_EDGE * ke
    m_free = RHO * Fraction(S) * (S - 1) / 2 * (2 ** M_CELL * kg * kg + N_EDGE * ke * ke)   # free pairs
    m_ff, m_fF, m_inc = 66, 144 * S * ke, 30 * (q1 + Q2) * cellsum   # fixed-fixed, fixed-free, within an incidence
    m_W = m_free + m_ff + m_fF + m_inc
    v_leaf = 12 * RHO * Fraction(S) ** 3 * ke * ke * cellsum   # pairs sharing a free leaf (at most 12 widened cells)
    v_fix = 5184 * Fraction(S) ** 2 * ke * ke + 1584 * S * ke   # pairs sharing a leaf, one or two of them alpha's
    v_gen = 24 * Fraction(S) ** 3 * (q1 + Q2) * (12 * ke / (r - 1)) ** 2   # pairs sharing a generator
    v_W = m_W + v_leaf + v_fix + v_gen + 1
    cap_W = cdiv(m_W.numerator, m_W.denominator) + SLACK_W
    d_W = v_W / Fraction(cap_W - m_W) ** 2
    # pair loads (Markov): every pair of leaves other than the colliding pair itself contributes at most RHO times its
    # independent-leaf probability (proof.md 6.4: conditioning, or the other node's independence)
    kb_m = KB(KEYBITS + D1 + D2)   # 193 bits: d1, d2, K2 and the cell
    sip = RHO * Fraction(5, 2) * (1 + 11 * (d_wj + d_lab)) * kb_m   # per leaf: the pairs of two leaves of one incidence
                                   # (5 per leaf), the second within W or D' at most 12 labels (6.4), no home condition
    other_h = RHO * (Fraction(S) * S * KBU(D1 + D2, EPS_HOME) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) + sip * S
    d_cl = other_h / (CAND_CAP - 4)
    other_g = RHO * Fraction(n0) * S * kb_m * 2 + sip * n0
    d_rq = other_g / (REQ_CAP - 8)
    other_v = RHO * (2 + 2 * 8) * (Fraction(S, SIDE) * S * KB(D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) \
        + 18 * sip * S  # a reply counts <= 8
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_wj + d_j + d_f + d_lab + 2 * d_p + d_W + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)  # RNDR is the rules' ideal coin source: its failure branch is charged, never taken
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_wj=d_wj, d_j=d_j, d_f=d_f, d_lab=d_lab,
                d_p=d_p, d_W=d_W, m_W=m_W, v_W=v_W, cap_W=cap_W, m_free=m_free, m_ff=m_ff, m_fF=m_fF, m_inc=m_inc,
                v_leaf=v_leaf, v_fix=v_fix, v_gen=v_gen, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def,
                rows=rows, glob=glob, lo=lo, m_p=m_p)


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
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p, parallelogram 2A - B wide "
      "in u and A in v; box factor below 1 + 2^-17, one side's below 1 + 2^-18",
      FA * FA - FA * FB + FB * FB == p and (FB * BETA_P - (FA - FB)) % p == 0 and ((FA - FB) * BETA_P + FA) % p == 0
      and (BOX + WIDE + FA) * (BOX + WIDE + 2 * FA - FB) * 2 ** 200 < p * p * (1 + Fraction(1, 2 ** 17))
      and (BOX + WIDE) * 2 ** 100 < p * (1 + Fraction(1, 2 ** 18)))
check("q1 = %d is the least with success bound >= 0.390000001" % Q1, T["lo"] >= THRESH and success_terms(Q1 - 1)["lo"] < THRESH)
check("leaves S = 2^72*6*q1*q2 = 2^%.6f" % lg(S), S == NODES * 6 * Q1 * Q2)
check("mu = |Gamma|*P1 = %.9f (cross-node pairs)" % float(T["mu"]), T["mu"] > Fraction(494, 1000))
check("Chen-Stein b1+b2 < 2^-95", T["cs"] < Fraction(1, 2 ** 95), "2^%.2f" % lg(T["cs"]))


def isum(*rs):
    """Interval sum. rs: (lo, hi) pairs. returns: (lo, hi)."""
    return sum(x[0] for x in rs), sum(x[1] for x in rs)


TFP_TERMS = ((0, 1), (0, 7), (0, 56))   # a truncated product's T(z) - A (proof.md 3.3): e, c_low, d
RAW = isum(*TFP_TERMS)                  # a stored fraction of gamma
COR = isum(RAW, (-1, 0))                # plus a correction word: T(c - y) = T(c) - T(y) - {0, 1}
U3, UC3 = (RAW, RAW, (-2 - 2 * RAW[1], -2 * RAW[0])), (COR, COR, (-2 - 2 * COR[1], -2 * COR[0]))   # U_2 = -U_0 - U_1
check("3.3's fraction ranges from the terms e 0..1, c_low 0..7, d 0..56: T - estimate in [0, 64], corrected [-1, 64]; "
      "U_2 in [-130, 0], corrected [-130, 2]", RAW == (0, 64) and COR == (-1, 64) and U3[2] == (-130, 0)
      and UC3[2] == (-130, 2))
WR = [[isum(UC3[m], UC3[(m + t) % 3], U3[(m + 2 * t) % 3], (0, 2)) for t in range(3)] for m in (0, 1)]   # W_m(t) + kappa
NR = [(-WR[0][t][1] - WR[1][t][1] - 1, -WR[0][t][0] - WR[1][t][0] + 1) for t in range(3)]   # N' = -e0 - e1 - k + 1
LABEL_ERR = [(min(x[0] for x in rs), max(x[1] for x in rs)) for rs in (WR[0], WR[1], NR)]
check("3.3's error table from the fraction ranges: W0 per t [-2, 194], [-132, 130], [-131, 132]; W1 [-2, 194], "
      "[-131, 132], [-132, 130]; N' [-389, 5], [-263, 264], [-263, 264]; over t [-132, 194], [-132, 194], [-389, 264]: "
      "D' 1305 residues (d_lab's), widening 654 (WIDE's)",
      WR == [[(-2, 194), (-132, 130), (-131, 132)], [(-2, 194), (-131, 132), (-132, 130)]]
      and NR == [(-389, 5), (-263, 264), (-263, 264)] and LABEL_ERR == [(-132, 194), (-132, 194), (-389, 264)]
      and sum(hi - lo for lo, hi in LABEL_ERR) == 1305 and max(hi - lo for lo, hi in LABEL_ERR) + 1 == 654
      and WIDE == 654 * (p >> 128) + 4)
def chord(u):
    """Row chord of 3.3's canonical region C at height u (p = 1), from its inequalities. u: Fraction. returns: Fraction."""
    low = max(Fraction(0), min(u, (1 - u) / 2))                                   # v < u, 2v + u < 1
    high = max(Fraction(0), min(u, (2 - u) / 2) - max(Fraction(0), 1 - u))       # 1 - u <= v < u, 2v + u < 2
    return low + high


CUTS = (Fraction(0), Fraction(1, 3), Fraction(1, 2), Fraction(2, 3), Fraction(1))
PIECE = (lambda u: u, lambda u: (1 - u) / 2, lambda u: (3 * u - 1) / 2, lambda u: Fraction(1, 2))   # 3.3's quarters
GRID = [Fraction(n, 1200) for n in range(1200)]
CHORD_OK = all(chord(u) == PIECE[sum(u >= c for c in CUTS[1:4])](u) for u in GRID)
AREA = sum((CUTS[n + 1] - CUTS[n]) * (PIECE[n](CUTS[n]) + PIECE[n](CUTS[n + 1])) / 2 for n in range(4))
FALL = sum(max(Fraction(0), PIECE[n](CUTS[n]) - PIECE[n](CUTS[n + 1])) for n in range(4)) + sum(
    max(Fraction(0), PIECE[n](CUTS[n + 1]) - PIECE[(n + 1) % 4](CUTS[n + 1] % 1)) for n in range(4))
check("Lemma U step 1: C's chord from its inequalities equals 3.3's four pieces at 1,200 heights; area p^2/3, total fall "
      "7p/12 (the wrap at p included)", CHORD_OK and AREA == Fraction(1, 3) and FALL == Fraction(7, 12),
      "area %s p^2, fall %s p" % (AREA, FALL))
WIDTH = max(hi - lo for lo, hi in LABEL_ERR) + 1
check("Lemma U step 2: outside W a leaf emits path j*'s pair, within the widest error range of 3.3's table, %d units "
      "(WIDE); beta has order 3, three values of x per canonical v" % WIDTH, WIDTH == 654 and WIDE == WIDTH * (p >> 128) + 4
      and BETA_P != 1 and pow(BETA_P, 3, p) == 1)
B1, B2 = (FB, FA - FB), (FA - FB, -FA)
VERTS = [(0, 0), B1, B2, (B1[0] + B2[0], B1[1] + B2[1])]
D2MAX = max((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2 for x in VERTS for y in VERTS)
check("Lemma U step 3: kappa, an expanded box's area over (p/2^100)^2, below 1 + 2^-17.64; diam(P) from its vertices "
      "below 2^128.6", Fraction((BOX + WIDE + FA) * (BOX + WIDE + 2 * FA - FB) * 2 ** 200, p * p) == KAPPA
      and KAPPA - 1 < pow2_interval(Fraction(-1764, 100))[0] and math.isqrt(D2MAX) + 1 == DIAM
      and DIAM < pow2_interval(Fraction(1286, 10))[0], "kappa - 1 = 2^%.4f, diam(P) = 2^%.4f" % (lg(KAPPA - 1), lg(DIAM)))


def rows_first(a, b):
    """Lemma U step 4's factor less 1 and C's edge and top terms, for digits a, b. returns: (factor - 1, terms)."""
    f = a + b
    fac = KAPPA * (1 + Fraction(7, 4) * Fraction(2 ** b, 2 ** 100) + 6 * Fraction(2 ** a, 2 ** 100)) - 1
    terms = Fraction(3 * 2 ** (100 - b) * 4 * 3 * DIAM * HU * 2 ** f, p * p) + Fraction(
        3 * (FA - FB) * ((p // 2 + FA) * HV * 2 ** 100 // (2 ** a * p) + HV + 1) * 2 ** f, p * p)
    return fac, terms


RF93, RF36 = rows_first(D1, D2 + M_CELL), rows_first(D1, D2)
check("Lemma U step 4: rows first, kappa (1 + 1.75 2^(b-100) + 6 2^(a-100)) - 1 for the cells (36, 93) and a home (36, "
      "36), with C's edges and top below 2^-35", RF93[1] < Fraction(1, 2 ** 35) and RF36[1] < Fraction(1, 2 ** 35),
      "(36, 93) 2^%.4f, (36, 36) 2^%.4f; edges and top 2^%.2f" % (lg(RF93[0]), lg(RF36[0]), lg(RF93[1])))
EDGE_K = [k for k in list(range(0, 5)) + list(range(2 ** M_CELL - 9, 2 ** M_CELL))
          if k << 36 <= 2 ** 36 or (k << 36) + 2 ** 36 - 1 >= 2 ** 93 - 2 ** 38 - 1]   # some d2 reaches a strip
BOXL = Fraction(HV * HU, p)
STRIP = [3 * n * 2 ** (100 - m) * BOXL * 2 ** (D1 + D2 + M_CELL) / p for n, m in ((1, D1), (7, D2 + M_CELL), (128, D2 + M_CELL))]
check("Lemma U step 5: W's strips in path coordinates: the horizontal one reaches only edge cells, N_EDGE = %d per home "
      "(K3 mod 2^57 in {0, 1} or the top five), one row each; the vertical one 7 columns, the slanted ones 128 boxes per "
      "row" % len(EDGE_K), len(EDGE_K) == N_EDGE == 7 and EDGE_K[:2] == [0, 1],
      "per cell: horizontal 2^%.2f (edge cells), vertical 2^%.2f, slanted 2^%.2f" % tuple(lg(x) for x in STRIP))
check("Lemma U step 6: eps (eps_u) below 2^-6.19 for the cells (36, 93), 2^-4.75 for edge cells, 2^-17.6 for a home",
      EPS_GEN < pow2_interval(Fraction(-619, 100))[0] and EPS_EDGE < pow2_interval(Fraction(-475, 100))[0]
      and EPS_HOME < pow2_interval(Fraction(-176, 10))[0], "2^%.4f 2^%.4f 2^%.4f" % (lg(EPS_GEN), lg(EPS_EDGE),
                                                                                   lg(EPS_HOME)))
check("walk budget, mean m_W = m_free + m_ff + m_fF + m_inc (success.md 6.4)", T["m_W"] == T["m_free"] + T["m_ff"]
      + T["m_fF"] + T["m_inc"], "m_free 2^%.4f, m_ff %d, m_fF %.3f, m_inc 2^%.2f; m_W 2^%.4f (%.6f per record)" % (
          lg(T["m_free"]), T["m_ff"], float(T["m_fF"]), lg(T["m_inc"]), lg(T["m_W"]), T["m_W"] / N0))
check("walk budget, variance v_W = m_W + v_leaf + v_fix + v_gen + 1", T["v_W"] == T["m_W"] + T["v_leaf"] + T["v_fix"]
      + T["v_gen"] + 1, "v_leaf 2^%.2f, v_fix %.3f, v_gen 2^%.1f; v_W 2^%.2f" % (lg(T["v_leaf"]), float(T["v_fix"]),
                                                                                lg(T["v_gen"]), lg(T["v_W"])))
check("walk budget: cap_W = ceil(m_W) + 2^40 = %d, rate (v_W/(cap_W - m_W)^2) < 2^-23" % T["cap_W"],
      T["d_W"] < pow2_interval(Fraction(-23))[0] and T["cap_W"] == cdiv(T["m_W"].numerator, T["m_W"].denominator)
      + SLACK_W and 0 < T["cap_W"] < 2 ** 63, "d_W 2^%.3f" % lg(T["d_W"]))
check("wrong class index rate per leaf < 2^-58.47", T["d_wj"] < pow2_interval(Fraction(-5847, 100))[0], "2^%.4f" % lg(T["d_wj"]))
for name, n, bad, mult in SITES:
    k = -(-n // Q2)
    print("  site %-52s %-6s bad <= 2^%.2f, x%d/(r-1), each <= 2^%.2f" % (name, ("%dq2" % k if k > 1 else "q2") + (
        "%+d" % (n - k * Q2) if n != k * Q2 else ""), up2(bad), mult, up2(Fraction(bad * mult, r - 1))))
check("site table: 18 row shifts (3 by SHXT), 10 folds, 2 square folds per incidence, less 3 + 3 exact at l = 0",
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
check("fixup table: per incidence 18 row shifts, 10 folds, 2 square folds, 2 FS and 1 FADDW, less 3 + 3 exact at l = 0",
      FIX_RUNS == {m: FIX_SITES[m] * Q2 - EXACT_AT_L0.get(m, 0) for m in FIX_SITES}
      and all(bad == dict(SHX=FIX_ROW, FDX=FIX_FOLD, F5=FIX_SQ, FS=FIX_SUB, FADDW=FIX_SUB)[m] for _, m, _, bad, _ in FIXUPS)
      and all(FIX_SITES[m] == SITES_PER_INCIDENCE[m] for m in ("SHX", "FDX", "F5")))
check("a call's fixup word (C0 per fixup, at most 33 q2 fixups) stays below 2^63", FIX_BATCH * C0 < 2 ** 63,
      "2^%.3f" % lg(FIX_BATCH * C0))
check("mean fixups per node (clean batches) q1*E_FIX = 2^%.3f < 2^28.68; cap 2^45 > 2^16 times it plus 33 q2" % lg(Q1 * E_FIX),
      Q1 * E_FIX < pow2_interval(Fraction(2868, 100))[0] and DELTA_F > 2 ** 16 * Q1 * E_FIX + FIX_BATCH)
check("fixup meter rate (Markov) < 2^-16.32", T["d_f"] < pow2_interval(Fraction(-1632, 100))[0], "2^%.3f" % lg(T["d_f"]))
check("label-error rate per leaf < 2^-16.65", T["d_lab"] < pow2_interval(Fraction(-1665, 100))[0],
      "2^%.4f" % lg(T["d_lab"]))
check("batch defect rate (site table, D = 0, X = O) < 2^-126.05", T["d_batch"] < pow2_interval(Fraction(-12605, 100))[0],
      "2^%.3f" % lg(T["d_batch"]))
check("j = 1, 2 meter rate < 2^-30", T["d_j"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_j"]))
check("pair-slot overflow rate, each hop, < 2^-19.9", T["d_p"] < pow2_interval(Fraction(-199, 10))[0],
      "2^%.3f" % lg(T["d_p"]))
check("candidate cap rate < 2^-28, request cap rate < 2^-37",
      T["d_cl"] < pow2_interval(Fraction(-28))[0] and T["d_rq"] < pow2_interval(Fraction(-37))[0],
      "2^%.3f 2^%.3f" % (lg(T["d_cl"]), lg(T["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-36.2", T["d_rl"] < pow2_interval(Fraction(-362, 10))[0],
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
row("Y generators [c]P with their table entries: q2 * 677570 per node", NODES * Q2 * GEN_MUL, True)
row("R_g = [b]Q and X = R_g + [a]P with their header values: (q1+1) * 677570 per node", NODES * (Q1 + 1) * GEN_MUL)
row("rows of beta: 72 per node", NODES * HEADER, True)
row("LC tables (init3) and RC tables (tinit): 2 * (4 * 2^36 + 1) + 2 * (3 * 2^36 + 1) per node",
    NODES * 2 * (INIT3_ENTRY * SIDE + 1 + TINIT_ENTRY * SIDE + 1), True)
row("Q broadcast: 2 messages to every node over two hops at 40", 2 * (SIDE + NODES) * MSG)
row("forward: q1 * (98 q2 + 5)", NODES * Q1 * (FWD_INC * Q2 + FWD_FIX))
row("inversion: q1 * (23,646 + 16)", NODES * Q1 * (INV + INV_CALL))
row("backward: q1 * (1297.5 q2 + 5)", int(NODES * Q1 * (BWD_INC * Q2 + BWD_FIX)))
row("arithmetic fixups: 2^45 * 4 per node", NODES * DELTA_F * max(FIX_EXTRA.values()))
row("key pass without emits: q1 * (98.5 q2 + 5)", int(NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX)))
row("j = 1 and j = 2 key paths: cap_j * 5 per node", NODES * CAP_J * JPATH)
row("batch set-up: q1 * 160", NODES * Q1 * BATCH_SETUP)
row("hop 1, the emits and full blocks (flush, in-flight count, dlvk): S * (8 + 523/128)", S * (EMIT * K + BLOCK) // K)
row("hop 2, h2s and full blocks: S * (8.25 + 523/128) + 11 * 2^36 + 4 per node",
    S * (H2S_REC * K + LOOP8 * K // 8 + BLOCK) // K + NODES * (H2S_SLOT * SIDE + H2S_FIX))
row("final blocks and cursors, in-flight counts and deliveries, both hops: 2 * (541 * 2^36 + 2) per node",
    NODES * 2 * (FINAL * SIDE + FIN_FIX))
row("head table zeroed by h2zero: 1026 * 2^46 + 1 per node", NODES * (ZERO_PASS * (8 * 2 ** M_CELL // ZERO_BYTES) + 1))
row("chained grouping h2g: S * 5.25 + 12 * 2^36 + 3 per node",
    S * H2G_REC + S * LOOP8 // 8 + NODES * (H2G_SLOT * SIDE + H2G_FIX))
row("walk steps: 5 * cap_W per node", NODES * H2G_STEP * T["cap_W"])
row("candidate appends: 6 * (2^20 + 1) per node", NODES * H2G_HIT * (CAND_CAP + 1))
row("candidates: 2^20 * 1557745 per node", NODES * CAND_CAP * RESOLVE)
row("requests served: 2^29 * 2^10 per node", NODES * REQ_CAP * REQ_SERVE)
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
    """Smallest multiple of 10^-12 whose power of two is at least w, verified by interval bounds.

    w: Fraction > 1.
    returns: Fraction claim.
    """
    c = Fraction(math.ceil(lg(w) * 10 ** 12), 10 ** 12)
    for _ in range(50):
        if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 12))[1] < w:
            return c
        c += Fraction(1, 10 ** 12) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 12)
    raise AssertionError("claim not certified")


CLAIM = claim_of(W)
check("claim = least 12-decimal c with 2^c >= R/1000", pow2_interval(CLAIM)[0] >= W
      and pow2_interval(CLAIM - Fraction(1, 10 ** 12))[1] < W, "%d/10^12" % (CLAIM * 10 ** 12))
PER_LEAF = Fraction(R_TOTAL, S)
check("work per leaf = %.7f instructions" % float(PER_LEAF), 286 < PER_LEAF < 287)

# ---------------------------------------------------------------- 4. memory (bytes per node, layout per phase)
SLOTS = SIDE * CAP_P * 32     # 2^36 pair slots of cap_p records: R1 (hop-1 receive) or R2 (hop-2 receive)
FIXED = [("code, constants, stack", 2 ** 20), ("barrier words: 5 granules of at most 2 KiB", 5 * 2048),
         ("header H and the rows of beta below it", 1024), ("X table: H's rows and values, -T(a), -T(beta a), y^2 per X", Q1 * 512),
         ("Y table: y^2, y, b, b^3, 3b^2, b^4", Q2 * 192), ("scalars a, b, c", GENS * 32), ("batch entries", Q2 * 256),
         ("the ports' message blocks", 2 * 64), ("LC tables at 0 and 2^39", 2 * SIDE * 8),
         ("RC tables at 2^43 and 2^43 + 2^39", 2 * SIDE * 8), ("record blocks at 2^49 and 2^50", 2 * SIDE * 2 ** 13),
         ("final-cursor tables E1 and E2", 2 * SIDE * 32), ("candidate list", CAND_CAP * 16),
         ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64), ("relay buffer", RELAY_CAP * 64),
         ("phase-4 count, offset and cursor tables", 3 * SIDE * 32)]
FIXED_BYTES = sum(b for _, b in FIXED)
TABLE = 2 ** 60               # the head table: 2^57 words at 2^60, in R1's range (dead after hop 2)
PHASE3 = [("head table at 2^60 (h2zero, h2g)", 8 * 2 ** M_CELL)]
PHASE3_BYTES = sum(b for _, b in PHASE3)
for name, b in FIXED + [("R1 at 2^52: hop-1 pair slots (phases 1-2)", SLOTS),
                        ("R2 at 2^52 + |R1|: hop-2 pair slots (phases 2-4)", SLOTS)] + PHASE3:
    print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("fixed regions lie below 2^52", FIXED_BYTES <= 2 ** 52, "2^%.4f" % lg(FIXED_BYTES))
ALIGNED = 4 * SIDE * 8 + 2 * SIDE * 2 ** 13
check("aligned tables: LC tables [0, 2^39) and [2^39, 2^40), RC tables [2^43, 2^44), record blocks [2^49, 2^50) and "
      "[2^50, 2^51); the other fixed regions fit between 2^44 and 2^49", SIDE * 8 == 2 ** 39 and SIDE * 2 ** 13
      == 2 ** 49 and FIXED_BYTES - ALIGNED <= 2 ** 49 - 2 ** 44, "2^%.4f" % lg(FIXED_BYTES - ALIGNED))
check("record block: header and k slots of 32 bytes within 2^13; slot k, at +32k = 2^12, is its only slot whose bits "
      "5..11 are zero (TST #0xfe0); cap_p is a multiple of k (the keep lemma, proof.md 4)",
      32 + 32 * K <= 2 ** 13 and all(((32 * n) & 0xfe0 == 0) == (n == K) for n in range(1, K + 1)) and CAP_P % K == 0
      and S % K == 0)
check("bwd's and keys' loop end: the q2 256-byte batch entries fill the lower 2^33 bytes of their 2^34-aligned block",
      Q2 * 256 == 2 ** 33)
check("record fields: k < 2^28 in w2 bits 36-63; the sign's code 2 or 5 in w3 bits 36-38 and l < 2^25 (l - 1 mod 2^25 "
      "for s = -1) in bits 39-63", Q1 <= 2 ** 28 and Q2 <= 2 ** 25)
check("R1 and R2 end below 2^64", 2 ** 52 + 2 * SLOTS < 2 ** 64, "2^%.6f" % lg(2 ** 52 + 2 * SLOTS))
check("phase 3: the head table, 2^57 words at 2^60 (a multiple of 2^60, so one BFI forms a cell's address), lies in "
      "R1's range [2^52, 2^52 + |R1|), dead after hop 2", TABLE % 2 ** 60 == 0 and 2 ** 52 <= TABLE
      and TABLE + PHASE3_BYTES <= 2 ** 52 + SLOTS, "2^%.4f" % lg(PHASE3_BYTES))
check("hop-1 cursors cannot wrap: the last slot's base plus 32 per leaf of a node stays below 2^64 (hop-2 cursors "
      "stop at their end)", 2 ** 52 + 32 * CAP_P * (SIDE - 1) + 32 * N0 < 2 ** 64,
      "2^%.4f" % lg(2 ** 52 + 32 * CAP_P * (SIDE - 1) + 32 * N0))
NODE_BYTES = FIXED_BYTES + 2 * SLOTS   # phase 2, the largest: fixed regions, R1 and R2
PEAK = NODES * NODE_BYTES
check("peak memory 2^72 * (fixed + R1 + R2) = 2^%.4f bytes" % lg(PEAK), PEAK < 2 ** 136)
check("pair-slot cap exceeds its crude mean: cap_p > m_p + 2^20", CAP_P > T["m_p"] + 2 ** 20)
check("preprocessing (rows marked P) = 2^%.7f units, included in the total" % lg(Fraction(R_PRE, C)), R_PRE < R_TOTAL)
# the figures declared in proof.md and submission.json: the claim and each supplementary upper bound rounded up (the
# least multiple of 10^-12 or 10^-4 at or above the exact value), the success figure rounded down
DECLARED = dict(time_log2=Fraction(124896269329358, 10 ** 12), memory_log2_bytes=Fraction(1359987, 10 ** 4),
                preprocessing_log2=Fraction(1064648, 10 ** 4), success=Fraction(390000001, 10 ** 9),
                per_leaf=Fraction(2865721, 10 ** 4))


def least_at_or_above(d, x, scale):
    """Whether 2^d >= x and 2^(d - 1/scale) < x, both certified by the interval bounds of pow2_interval.

    d: Fraction; x: positive int or Fraction; scale: denominator of the decimal grid.
    returns: bool.
    """
    return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x


check("declared time_log2 %d/10^12: the claim, least 12-decimal value with 2^c >= R/1000" % (DECLARED["time_log2"] * 10 ** 12),
      DECLARED["time_log2"] == CLAIM and least_at_or_above(DECLARED["time_log2"], W, 10 ** 12))
check("declared memory 2^%s bytes: least 4-decimal bound >= peak (2^%.7f)" % (format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK)),
      least_at_or_above(DECLARED["memory_log2_bytes"], PEAK, 10 ** 4))
check("declared preprocessing 2^%s units: least 4-decimal bound >= rows marked P / 1000" % format(float(DECLARED["preprocessing_log2"]), ".4f"),
      least_at_or_above(DECLARED["preprocessing_log2"], Fraction(R_PRE, C), 10 ** 4))
check("declared success 0.390000001 <= the bound", DECLARED["success"] <= T["lo"])
check("declared work per leaf %s: least 4-decimal bound >= R/S (%.7f)" % (format(float(DECLARED["per_leaf"]), ".4f"), float(PER_LEAF)),
      DECLARED["per_leaf"] >= PER_LEAF > DECLARED["per_leaf"] - Fraction(1, 10 ** 4))
ANALYTIC = sum(v for n, v, _ in rows if n.split(":")[0] in ("program delivery, node and port state", "coins",
               "Y generators [c]P with their table entries", "R_g = [b]Q and X = R_g + [a]P with their header values",
               "rows of beta", "batch set-up",
               "Q broadcast", "candidates", "requests served", "phase-4 count exchanges", "barriers", "phase control",
               "output",
               "inversion")) - NODES * Q1 * INV
# itemized work behind the allowances (proof.md section 7)
# on the exact arithmetic of proof.md 3.4: fe_mul_f and fe_sqr_f as executed (115 and 82) plus their calls (pointers and
# BL: 4 and 3, as in fe_inv), and FS or FADDW (9) with both operands brought below p (2 x 9: ADDS chain of C0, CSELs)
MULX, SQRX, ADDX = 115 + 4, 82 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
# 256 doublings and mixed additions, two conversions (inversion, 3 products, 2 squares, 10 additions each), 64 set-up
GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_TEST = 256 * 2 * (9 + 4) + 64   # exceptional cases: per iteration Z and H brought below p (9, as in ADDX) and tested
                                    # for zero (ORR, ORR, ORR, CBZ); 64 for the paths they select (Q = O; H = 0)
GEN_CTRL = 1600   # 256 iterations of scalar bit test, branch and pointer step (at most 6 each: 1,536), set-up and return
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
# frame (2 LDRs) and H's address (4 MOVZ/MOVK); X's header (rows of 2a^3, 2a, a^2, y_X, a^4, -T(a), -T(beta a): 30
# LDPs, 30 STPs) into H; fwd's arguments (q2, two pointers of at most 4 MOVZ/MOVK each,
# acc = 1 by 4 MOVs, y_X^2 by 2 LDPs) and BL; bwd's two pointers and BL (R is loaded by the inversion call); keys' entry
# pointer, slot end (LDR), index word (LDR, ADD), key-path count (LDR) and BL, and the count's STR after it; sp set to H
# for bwd and back to the driver's (4 MOVZ/MOVK and a MOV each); after fwd and after bwd, C0 (MOVZ,
# MOVK), the fixup word over C0 (UDIV) added to the node's fixup count (LDR, ADD, STR); before the batch, each meter's
# test (LDR of the count, its cap by 4 MOVZ/MOVK, CMP, B.HI); the loop (ADD and STR of k, q1 by MOVZ/MOVK, CMP, B.LO)
# and the X-table pointer step (ADD, STR)
BATCH_ITEM = (2 + 4) + 60 + (1 + 4 + 4 + 4 + 2 + 1) + (2 * 4 + 1) + (4 + 1 + 2 + 1 + 1 + 1) + 2 * (4 + 1) \
    + 2 * (2 + 1 + 3) + 2 * (1 + 4 + 1 + 1) + (2 + 2 + 1 + 1) + 2
check("batch set-up: at most %d instructions, within its allowance of 160" % BATCH_ITEM, BATCH_ITEM <= BATCH_SETUP)
check("barrier: at most 131 instructions per node, in-flight counts and deliveries included (2^10 per node is > 7x)",
      BARRIER_ITEM == 131 and BARRIER > 7 * NODES * BARRIER_ITEM)
# a phase-4 count exchange at one node: at most 12 per entry of its 2^36-entry received-count table (LDP, epoch CMP,
# CSEL, prefix ADD, zero test, cap CMP and branch, loop), and at most 88 per message (its count, the count message,
# the offset reply and the write, each at most 22 with its in-flight count and delivery) for at most 8 * 2^29 = 2^32
# messages (2^29 requests of up to 8 records)
EXCH_ITEM = 12 * SIDE + 4 * (14 + MSG_LINK) * (8 * REQ_CAP)
check("phase-4 count exchange: < 2^40.5 per node, within 2^41 (at most 8 * 2^29 messages)",
      EXCH_ITEM < pow2_interval(Fraction(405, 10))[0] and EXCH_ITEM < EXCH, "2^%.3f" % lg(EXCH_ITEM))
# a broadcast message (program, Q): two LDPs of its words, three STPs, ADD of 8 to its destination word, SUBS, B.NE,
# its in-flight count and delivery
check("a broadcast message costs at most 16, in-flight count and delivery included, within its price of 40",
      2 + 3 + 3 + MSG_LINK <= MSG)
# a record's three candidate leaves phi^t(X) + sY (t = 0, 1, 2, as w3 gives s but not t): D = x_Y^3 - x_X^3 (two
# squares, two products, a subtraction) and its inverse; Q_t = x_Y^2 + beta^t x_X x_Y + beta^(2t) x_X^2 (x_X x_Y and four
# products by beta, beta^2, six additions); 1/(x_Y - beta^t x_X) = Q_t/D (three products); beta x_X, beta^2 x_X; per t
# the slope, its square, x and y (two products, a square, four additions), and s y_Y - y_X once
CAND_LEAVES = INV + 5 * SQRX + 18 * MULX + 20 * ADDX
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 9 * 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, both records' three candidate leaves, nine unit tests (600 each), solve mod r, [x]P = Q = "
      "%d; charged above twice that" % CAND_ITEM, RESOLVE == 2 * CAND_ITEM + 1,
      "%.3fx" % (math.floor(10 ** 3 * RESOLVE / CAND_ITEM) / 10 ** 3))
PROGRAM_BYTES = 2 ** 16
check("program delivery: 2^16 bytes = 2^11 messages over 2 hops at 40 = 2^17.3 < node set-up 2^20",
      PROGRAM_BYTES // 32 * 2 * MSG < NODE_INIT)
check("analytic allowances are below 2^-13.7 of the total", ANALYTIC * pow2_interval(Fraction(137, 10))[1] < R_TOTAL,
      "2^%.3f" % lg(Fraction(ANALYTIC, R_TOTAL)))

print()
print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1, N0, lg(S), float(T["mu"]),
                                                                    math.floor(T["lo"] * 10 ** 10) / 10 ** 10))
print("R_total %d instructions  W 2^%.10f  claim %d/10^12  per leaf %.4f" % (R_TOTAL, lg(W), CLAIM * 10 ** 12,
                                                                      float(PER_LEAF)))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
    format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK), format(float(DECLARED["preprocessing_log2"]), ".4f"),
    lg(Fraction(R_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
    sys.exit(1 if FAILED else 0)
