"""Exact certificate for proof.md (parameters, E-work ledger, claim, memory)."""
import math
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction

FAILED = []


def check(name, ok, detail=""):
    print("%-66s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
    if not ok:
        FAILED.append(name)


def lg(x):
    f = Fraction(x)
    return math.log2(f.numerator) - math.log2(f.denominator)


def up6(x):
    """The least multiple of 10^-6 at or above log2 x: the certified worst-case figure."""
    return math.ceil(lg(x) * 10 ** 6) / 10 ** 6


def cdiv(a, b):
    return -((-a) // b)


def exp_interval(x, terms=60):
    assert abs(x) <= 1
    s, t = Fraction(1), Fraction(1)
    for k in range(1, terms + 1):
        t = t * x / k
        s += t
    rem = abs(t * x) / (terms + 1) * 3
    return s - rem, s + rem


def ln2_interval(bits=300):
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
    n = c.numerator // c.denominator
    f = c - n
    lo, hi = exp_interval(f * LN2[0])[0], exp_interval(f * LN2[1])[1]
    sc = Fraction(2) ** n
    return lo * sc, hi * sc


# ---------------------------------------------------------------- 1. constants
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
C0 = 2 ** 32 + 977
# ---------------------------------------------------------------- 2. success
# ---------------------------------------------------------------- 3. work ledger (instructions)
# ---------------------------------------------------------------- 5. grouped early abort (proof.md 2.3, 4, 6.6): E-work
# ---------------------------------------------------------------- 4. memory (bytes per node, layout per phase)
BETA_P = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
C = 1000
NODES = 2 ** 72
SIDE = 2 ** 36
Q2 = 2 ** 25
D1, D2, D3, D4 = 36, 36, 24, 46
KEYBITS = 110
THRESH = Fraction(390000001, 10 ** 9)


LEAVES = 6
FWD_INC, FWD_FIX = 98, 5

BWD_INC, BWD_FIX = Fraction(2571, 2), 5

FIX_SITES = dict(SHX=18, FDX=10, F5=2, FS=2, FADDW=1)
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4)
EMIT = 7
EMIT_INC = LEAVES * EMIT
KEY_INC, KEY_FIX = Fraction(281, 2), 5


K = 256
PB1C, PB2C, D4BASE, BLKC = 1 << 50, 1 << 51, 1 << 49, 1 << 14
FLUSH = 8

KEY_BLOCK_DROP, H2S_BLOCK_DROP = 7, 9

JPATH = 4

INV = 22549
TINIT_ENTRY = 3
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 14, 3, 2, 2


H2S_REC, H2S_SLOT = 7, 11
SCAN_RANGE = 26
H2S_FIX, BKT_FIX, SCAN_FIX = 4, 4, SCAN_RANGE * 2 ** D4 + 2

BKT_REC, BKT_DROP, BKT_SLOT = 8, 7, 12

PASS4_REC, PASS4_BKT, PASS4_FIX = 3, 11, 2
SCAN_REC = Fraction(195, 64)

LOOP8 = 2
ZERO4_PER4, PREFIX4_PER4, INIT3_BKT = 4, 10, 4

check("key pass: 140.5 = 2 signs * (6 LDP of the fractions + 28 sums + 3 leaves * (index ADD, CMP, B.HS, CMN, B.HS + "
      "emit 7)) + TBZ per two incidences", KEY_INC == 2 * (6 + 28 + 3 * (5 + EMIT)) + Fraction(1, 2))
check("j = 1 and j = 2 paths (B.HS taken; CMN, B.HS; ADDS, ADC, MVN of N'; emit; B back) exceed the j = 0 path "
      "(ADD, CMP, B.HS, CMN, B.HS; emit) by 4", (3 + 2 + 3 + EMIT + 1) - (5 + EMIT) == JPATH)
check("emit: AND or BIC (8*d1), LDADD of its LC entry (the slot), d2 and w3 (LSR or EOR, and ADD of the entry pointer "
      "shifted left by 31), two STPs, TBZ of bit 13 for slot k", EMIT == 1 + 1 + 2 + 2 + 1)
DONE = 2
DLV = 6
DLVK = 4 * K + 2
MSG_LINK = DONE + DLV
BLOCK = FLUSH + DONE + DLVK
FINAL = FIN_DEST + FIN_SENT + DONE + DLVK + MSG_LINK
L8 = Fraction(LOOP8, 8)
HOP1, HOP2 = EMIT + Fraction(BLOCK, K), H2S_REC + L8 + Fraction(BLOCK, K)
GROUP = BKT_REC + L8 + PASS4_REC + L8 + SCAN_REC
check("a full block: 1036 = flush 8 + in-flight count 2 + dlvk 1026 when sent, more than 7 or 9 when dropped; fin per "
      "destination at most 1053",
      BLOCK == 1036 and BLOCK > max(KEY_BLOCK_DROP, H2S_BLOCK_DROP) and FINAL == 1053
      and FINAL > FIN_DEST + FIN_HELD + MSG_LINK)
check("per record: hop 1 7 + 1036/256, hop 2 7.25 + 1036/256, grouping 14.546875; 36.890625 in all",
      HOP1 + HOP2 + GROUP == Fraction(4722, 128))


GEN_MUL = 645410
# grouped X generators also store the X affine point (x, y): 64 bytes, four more STPs.
GEN_MULX = GEN_MUL + 4
# a Y item is [c]P with its table entry only (no mixed addition, no X header): its itemization
# GEN_ITEM + GEN_TEST + GEN_CTRL + Y_ENTRY (checked in section 7); the base ledger R keeps GEN_MUL.
GEN_MULY = 642930


COIN_SCALAR = 64
BATCH_SETUP = 160

INV_CALL = 16
NODE_INIT = 2 ** 20
HEADER = 72
MSG = 40

CAND_CAP = 2 ** 28
RESOLVE = 1491011


REQ_CAP = 2 ** 29
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 34
EXCH = 2 ** 41
BARRIERS, BARRIER = 32, 2 ** 82
OUTPUT = 2 ** 20
DELTA_F = 2 ** 45
FIX_BATCH = sum(FIX_SITES.values()) * Q2
CAP_P = 2 ** 22 - 2 ** 12


P1 = Fraction(6 * (r - 1) + 1, r * r)
TWO_BOUND = Fraction(2, r - 1)


SHX_BAD = 2 ** 64 * 2 ** 33
FOLD_BAD = 2 ** 66 * C0
CP_BAD = (2 ** 66 + 1) * C0
SQ_BAD = C0 * C0 + 2 * C0
SITES = [
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
SITES_PER_INCIDENCE = dict(SHX=18, FDX=10, F5=2)
EXACT_AT_L0 = dict(SHX=3, FDX=3)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)


FIX_ROW = C0 * 2 ** 192
FIX_FOLD = (2 ** 66 + 1) * 2 ** 128 * C0

FIX_SQ = 2 ** 193
FIX_SUB = C0 * 2 ** 192
FIXUPS = [
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
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (2 ** 17 + 1), (r - 1) * 2 ** 17)

RHO = Fraction(r, r - 1) ** 4
FIXED_LEAVES = 2 * LEAVES
PARTNER_PAIRS = 2 * 10


def success_terms(q1):
    n0 = LEAVES * q1 * Q2
    S = NODES * n0
    gamma = Fraction(S * S - NODES * n0 * n0, 2)
    mu = gamma * P1
    deg = 2 * LEAVES * (q1 + Q2) * S
    cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
    mu_lo = Fraction((mu.numerator << 200) // mu.denominator, 2 ** 200)
    e_hi = exp_interval(-mu_lo, 40)[1]

    d_batch = D_SITES + Q2 * Fraction(6, r - 1) + Fraction(1, r)
    d_wj = TWO_BOUND * 23 * 2 ** 192

    d_f = q1 * E_FIX / (DELTA_F - FIX_BATCH) + q1 * d_batch
    d_lab = TWO_BOUND * 1305 * 2 ** 228


    m_p = n0 * KB(D1) + 1
    v_p = m_p + 1 + 60 * (q1 + Q2) * KB(D1)
    d_p = v_p / Fraction(CAP_P - FIXED_LEAVES - m_p) ** 2
    m_3 = S * KB(D1 + D2 + D3) + 1
    cap_3 = cdiv(m_3.numerator, m_3.denominator) + 2 ** 32
    v_3 = m_3 + 1 + 60 * (q1 + Q2) * KB(D1 + D2 + D3)
    d_3 = v_3 / Fraction(cap_3 - FIXED_LEAVES - m_3) ** 2


    kb_m = KB(KEYBITS + D1 + D2)
    lam_int = RHO * (S + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
    sip = RHO * Fraction(5, 2) * (1 + 11 * (d_wj + d_lab)) * kb_m

    other_h = RHO * (Fraction(S) * S * KB(D1 + D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) + sip * S
    d_cl = other_h / (CAND_CAP - 4)
    other_g = RHO * Fraction(n0) * S * kb_m * 2 + sip * n0
    d_rq = other_g / (REQ_CAP - 8)
    other_v = RHO * (2 + 2 * 8) * (Fraction(S, SIDE) * S * KB(D2) * kb_m + (2 * S + PARTNER_PAIRS) * kb_m) \
        + 18 * sip * S
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_wj + d_f + d_lab + 2 * d_p + d_3 + d_cl + d_rq + d_rl
    rows = (mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf))
    g_def = NODES * q1 * d_batch
    g_rndr = Fraction(0)
    glob = NODES * (q1 + Q2 + 1) * Fraction(1, 2 ** 254) + g_def + g_rndr
    lo = 1 - e_hi - cs - rows - glob
    return dict(n0=n0, S=S, gamma=gamma, mu=mu, cs=cs, e_hi=e_hi, d_batch=d_batch, d_wj=d_wj, d_f=d_f, d_lab=d_lab,
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
BOX = (p >> 100) + 1
WIDE = 654 * (p >> 128) + 4
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p, parallelogram 2A - B wide "
      "in u and A in v; box factor below 1 + 2^-17, one side's below 1 + 2^-18",
      FA * FA - FA * FB + FB * FB == p and (FB * BETA_P - (FA - FB)) % p == 0 and ((FA - FB) * BETA_P + FA) % p == 0
      and (BOX + WIDE + FA) * (BOX + WIDE + 2 * FA - FB) * 2 ** 200 < p * p * (1 + Fraction(1, 2 ** 17))
      and (BOX + WIDE) * 2 ** 100 < p * (1 + Fraction(1, 2 ** 18)))
check("q1 = %d is the least with success bound >= 0.390000001" % Q1, T["lo"] >= THRESH and success_terms(Q1 - 1)["lo"] < THRESH)
check("leaves S = 2^72*6*q1*q2 = 2^%.6f" % lg(S), S == NODES * 6 * Q1 * Q2)
check("mu = |Gamma|*P1 = %.9f (cross-node pairs)" % float(T["mu"]), T["mu"] > Fraction(494, 1000))
check("Chen-Stein b1+b2 < 2^-95", T["cs"] < Fraction(1, 2 ** 95), "2^%.2f" % lg(T["cs"]))
check("intruder rate per colliding pair < 2^-12.7", T["lam_int"] < pow2_interval(Fraction(-127, 10))[0], "2^%.3f" % lg(T["lam_int"]))


def isum(*rs):
    return sum(x[0] for x in rs), sum(x[1] for x in rs)


TFP_TERMS = ((0, 1), (0, 7), (0, 56))
RAW = isum(*TFP_TERMS)
COR = isum(RAW, (-1, 0))
U3, UC3 = (RAW, RAW, (-2 - 2 * RAW[1], -2 * RAW[0])), (COR, COR, (-2 - 2 * COR[1], -2 * COR[0]))
check("3.3's fraction ranges from the terms e 0..1, c_low 0..7, d 0..56: T - estimate in [0, 64], corrected [-1, 64]; "
      "U_2 in [-130, 0], corrected [-130, 2]", RAW == (0, 64) and COR == (-1, 64) and U3[2] == (-130, 0)
      and UC3[2] == (-130, 2))
WR = [[isum(UC3[m], UC3[(m + t) % 3], U3[(m + 2 * t) % 3], (0, 2)) for t in range(3)] for m in (0, 1)]
NR = [(-WR[0][t][1] - WR[1][t][1] - 1, -WR[0][t][0] - WR[1][t][0] + 1) for t in range(3)]
LABEL_ERR = [(min(x[0] for x in rs), max(x[1] for x in rs)) for rs in (WR[0], WR[1], NR)]
check("3.3's error table from the fraction ranges: W0 per t [-2, 194], [-132, 130], [-131, 132]; W1 [-2, 194], "
      "[-131, 132], [-132, 130]; N' [-389, 5], [-263, 264], [-263, 264]; over t [-132, 194], [-132, 194], [-389, 264]: "
      "D' 1305 residues (d_lab's), widening 654 (WIDE's)",
      WR == [[(-2, 194), (-132, 130), (-131, 132)], [(-2, 194), (-131, 132), (-132, 130)]]
      and NR == [(-389, 5), (-263, 264), (-263, 264)] and LABEL_ERR == [(-132, 194), (-132, 194), (-389, 264)]
      and sum(hi - lo for lo, hi in LABEL_ERR) == 1305 and max(hi - lo for lo, hi in LABEL_ERR) + 1 == 654
      and WIDE == 654 * (p >> 128) + 4)
check("wrong class index rate per leaf < 2^-58.47", T["d_wj"] < pow2_interval(Fraction(-5847, 100))[0], "2^%.4f" % lg(T["d_wj"]))
check("site table: 18 row shifts (3 by SHXT), 10 folds, 2 square folds per incidence, less 3 + 3 exact at l = 0",
      sum(n for _, n, b, _ in SITES if b == SHX_BAD) == SITES_PER_INCIDENCE["SHX"] * Q2 - EXACT_AT_L0["SHX"]
      and sum(n for _, n, b, _ in SITES if b in (FOLD_BAD, CP_BAD)) == SITES_PER_INCIDENCE["FDX"] * Q2
      - EXACT_AT_L0["FDX"]
      and sum(n for _, n, b, _ in SITES if b == SQ_BAD) == SITES_PER_INCIDENCE["F5"] * Q2 and len(SITES) == 11)
check("bad sets: rows 2^97, folds below 2^98.02, square folds below 2^64.01",
      SHX_BAD == 2 ** 97 and all(b < pow2_interval(Fraction(e, 100))[0] for b, e in ((FOLD_BAD, 9802), (CP_BAD, 9802),
                                                                                     (SQ_BAD, 6401))))
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
      "2^%.3f" % lg(T["d_lab"]))
check("batch defect rate (site table, D = 0, X = O) < 2^-126.05", T["d_batch"] < pow2_interval(Fraction(-12605, 100))[0],
      "2^%.3f" % lg(T["d_batch"]))
check("pair-slot overflow rate, each hop, < 2^-19.9", T["d_p"] < pow2_interval(Fraction(-199, 10))[0],
      "2^%.3f" % lg(T["d_p"]))
check("d3 bucket overflow rate < 2^-30", T["d_3"] < Fraction(1, 2 ** 30), "2^%.2f" % lg(T["d_3"]))
check("candidate cap rate < 2^-23.4, request cap rate < 2^-25.9",
      T["d_cl"] < pow2_interval(Fraction(-234, 10))[0] and T["d_rq"] < pow2_interval(Fraction(-259, 10))[0],
      "2^%.3f 2^%.3f" % (lg(T["d_cl"]), lg(T["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-25.2", T["d_rl"] < pow2_interval(Fraction(-252, 10))[0],
      "2^%.3f" % lg(T["d_rl"]))
check("whole-run row: a computational defect in any batch < 2^-26.5", T["g_def"] < pow2_interval(Fraction(-265, 10))[0],
      "2^%.3f" % lg(T["g_def"]))
check("whole-run rows (defective batch anywhere, coin rejection twice) < 2^-26.5",
      T["glob"] < pow2_interval(Fraction(-265, 10))[0], "2^%.3f" % lg(T["glob"]))
check("success >= 1 - e^-mu - (b1+b2) - rows - global >= 0.390000001", T["lo"] >= THRESH,
      "%.10f" % (math.floor(T["lo"] * 10 ** 10) / 10 ** 10))


GENS = Q1 + Q2 + 1
rows = []


def row(name, value, pre=False):
    assert isinstance(value, int) and value >= 0
    rows.append((name, value, pre))


row("node set-up: 2^20 per node", NODES * NODE_INIT, True)
row("coins: (q1+q2+1) scalars * 64 per node", NODES * GENS * COIN_SCALAR, True)
row("Y generators [c]P with their table entries: q2 * 645410 per node", NODES * Q2 * GEN_MUL, True)
row("R_g and X generators with their header values: (q1+1) * 645410 per node", NODES * (Q1 + 1) * GEN_MUL)
row("rows of beta: 72 per node", NODES * HEADER, True)
row("LC tables (init3) and RC tables (tinit): 2 * (4 * 2^36 + 1) + 2 * (3 * 2^36 + 1) per node",
    NODES * 2 * (INIT3_BKT * SIDE + 1 + TINIT_ENTRY * SIDE + 1), True)
row("Q broadcast: 2 messages to every node over two hops at 40", 2 * (SIDE + NODES) * MSG)
row("d4 counters zeroed: 2^46/4 * 4 + 1 per node", NODES * (2 ** D4 // 4 * ZERO4_PER4 + 1), True)
row("d3 cursor table: 4 * 2^24 + 1 per node", NODES * (INIT3_BKT * 2 ** D3 + 1), True)
row("forward: q1 * (98 q2 + 5)", NODES * Q1 * (FWD_INC * Q2 + FWD_FIX))
row("inversion: q1 * (22,549 + 16)", NODES * Q1 * (INV + INV_CALL))
row("backward: q1 * (1285.5 q2 + 5)", int(NODES * Q1 * (BWD_INC * Q2 + BWD_FIX)))
row("arithmetic fixups: 2^45 * 4 per node", NODES * DELTA_F * max(FIX_EXTRA.values()))
row("key pass without emits: q1 * (98.5 q2 + 5)", int(NODES * Q1 * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX)))
row("j = 1 and j = 2 key paths: n0 * 4 per node (every leaf, no meter)", NODES * N0 * JPATH)
row("batch set-up: q1 * 160", NODES * Q1 * BATCH_SETUP)
row("hop 1, emits and full blocks: S * (7 + 1036/256)", S * (EMIT * K + BLOCK) // K)
row("hop 2, h2s and full blocks: S * (7.25 + 1036/256) + 11 * 2^36 + 4 per node",
    S * (H2S_REC * K + LOOP8 * K // 8 + BLOCK) // K + NODES * (H2S_SLOT * SIDE + H2S_FIX))
row("final blocks and cursors, both hops: 2 * (1053 * 2^36 + 2) per node",
    NODES * 2 * (FINAL * SIDE + FIN_FIX))
row("bucket pass bkt: S * 8.25 + 12 * 2^36 + 4 per node", S * BKT_REC + S * LOOP8 // 8 + NODES * (BKT_SLOT * SIDE + BKT_FIX))
row("prefix sums of d4: 2^46/4 * 10 + 1 per node", NODES * (2 ** D4 // 4 * PREFIX4_PER4 + 1))
row("local pass pass4: S * 3.25 + 11 * 2^24 + 2 per node",
    S * PASS4_REC + S * LOOP8 // 8 + NODES * (PASS4_BKT * 2 ** D3 + PASS4_FIX))
row("scan_ranges: S * 195/64 + 26 * 2^46 + 2 + 8 * (2^28 + 1) per node", S * SCAN_REC.numerator // SCAN_REC.denominator
    + NODES * (SCAN_FIX + 8 * (CAND_CAP + 1)))
row("candidates: 2^28 * 1491011 per node", NODES * CAND_CAP * RESOLVE)
row("requests served: 2^29 * 2^10 per node", NODES * REQ_CAP * REQ_SERVE)
row("phase-4 count exchanges: 4 * 2^41 per node", NODES * 4 * EXCH)
row("barriers: 32 * 2^82", BARRIERS * BARRIER)
row("phase control: 2^10 per node", NODES * 2 ** 10)
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
check("work per leaf = %.7f instructions" % float(PER_LEAF), 287 < PER_LEAF < 289)


G_GROUPS = 114
CAP_PG = CAP_P // 8
# Unequal early-abort ladder (proof.md 2.3, 6.6, changes.md): rung g starts at
# round(q1*(g/G)^0.65). Rung count, barriers per rung and every kernel are the
# base package's; only the rung sizes move, lowering sum_g PHI[g] W_g.
LADDER_GAMMA = 0.65


def group_sizes(q1):
    ms = [0]
    for _g in range(1, G_GROUPS):
        _m = int(round(q1 * (_g / G_GROUPS) ** LADDER_GAMMA))
        if _m <= ms[-1]:
            _m = ms[-1] + 1
        if _m > q1 - (G_GROUPS - _g):
            _m = q1 - (G_GROUPS - _g)
        ms.append(_m)
    ms.append(q1)
    assert len(set(ms)) == G_GROUPS + 1
    qs = [ms[_g + 1] - ms[_g] for _g in range(G_GROUPS)]
    return qs, ms
HB = 57
HCHG_NEW = 7
HCHG_SLOT = 13
HCHG_FIXG = 2
HCHG_VISIT = 6
HCHG_END = 1
HCHG_VER = 5
HCHG_LS0 = 2
HCHG_LS1 = 4
check("hchg per record (LDP, BFI, LDR head, ADD next address, STR head into it, CBZ, final STR): 7; per visit (LDR K2, "
      "ADD, LDR next, CMP, B.EQ, CBZ or CBNZ): 6; the B to the insert after an even number of visits: 1",
      HCHG_NEW == 7 and HCHG_VISIT == 6 and HCHG_END == 1)
ZW = 4096
ZERO_PASS = ZW // 2 + 2
HOME_MAX = SIDE * CAP_P

STOP_CHECK = 5
NBAR_COUNT = 4 + 17 * G_GROUPS + 1
NBAR_ALLOW = 2048
PHASE_CTRL_G = 2 ** 16


# Hash walks are uncapped (no overflow event or row); a chain holds one home's
# records of one fold, at most HOME_MAX, which bounds the worst-case walks.
KB129 = KB(D1 + D2 + HB)
KB64 = KB(64)


def success_prefix(q1p, qs):
    n0p = LEAVES * q1p * Q2
    Sp = NODES * n0p
    gamma = Fraction(Sp * Sp - NODES * n0p * n0p, 2)
    mu = gamma * P1
    deg = 2 * LEAVES * (q1p + Q2) * Sp
    cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
    mu_lo = Fraction((mu.numerator << 200) // mu.denominator, 2 ** 200) if mu else Fraction(0)
    e_hi = exp_interval(-mu_lo, 40)[1]
    d_batch = D_SITES + Q2 * Fraction(6, r - 1) + Fraction(1, r)
    d_wj = TWO_BOUND * 23 * 2 ** 192
    d_f = q1p * E_FIX / (DELTA_F - FIX_BATCH) + q1p * d_batch
    d_lab = TWO_BOUND * 1305 * 2 ** 228
    m_p2 = n0p * KB(D1) + 1
    v_p2 = m_p2 + 1 + 60 * (q1p + Q2) * KB(D1)
    d_p2 = v_p2 / Fraction(CAP_P - FIXED_LEAVES - m_p2) ** 2
    d_p1 = Fraction(0)
    for _q in qs:
        _n0 = LEAVES * _q * Q2
        _m1 = _n0 * KB(D1) + 1
        _v1 = _m1 + 1 + 60 * (_q + Q2) * KB(D1)
        d_p1 += _v1 / Fraction(CAP_PG - FIXED_LEAVES - _m1) ** 2
    kb_m = KB(KEYBITS + D1 + D2)
    lam_int = RHO * (Sp + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
    sip = RHO * Fraction(5, 2) * (1 + 11 * (d_wj + d_lab)) * kb_m
    other_h = RHO * (Fraction(Sp) * Sp * KB(D1 + D2) * kb_m + (2 * Sp + PARTNER_PAIRS) * kb_m) + sip * Sp
    d_cl = other_h / (CAND_CAP - 4)
    other_g = RHO * Fraction(n0p) * Sp * kb_m * 2 + sip * n0p
    d_rq = other_g / (REQ_CAP - 8)
    other_v = RHO * (2 + 2 * 8) * (Fraction(Sp, SIDE) * Sp * KB(D2) * kb_m + (2 * Sp + PARTNER_PAIRS) * kb_m) \
        + 18 * sip * Sp
    d_rl = other_v / (RELAY_CAP - 16)
    per_leaf = d_batch + d_wj + d_f + d_lab + d_p1 + d_p2 + d_cl + d_rq + d_rl
    rows_ = mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf)
    g_def = NODES * q1p * d_batch
    glob = NODES * (q1p + Q2 + 1) * Fraction(1, 2 ** 254) + g_def
    lo = 1 - e_hi - cs - rows_ - glob
    return dict(mu=mu, lo=lo, cs=cs, rows=rows_, glob=glob, Sp=Sp, n0p=n0p, d_p1=d_p1, d_p2=d_p2)


def grouped_full_lo(q):
    qs, _ = group_sizes(q)
    return success_prefix(q, qs)["lo"]


lo_q, hi_q = Q1, Q1
while grouped_full_lo(hi_q) < THRESH:
    hi_q += 4096
while hi_q - lo_q > 1:
    mid = (lo_q + hi_q) // 2
    if grouped_full_lo(mid) >= THRESH:
        hi_q = mid
    else:
        lo_q = mid
Q1G, MQ1 = group_sizes(hi_q)
Q1GD = hi_q
N0G = LEAVES * Q1GD * Q2
SG_FULL = NODES * N0G
check("grouped q1 = %d is least with the full grouped bound >= 0.390000001 (base %d)" % (Q1GD, Q1),
      grouped_full_lo(Q1GD) >= THRESH and (Q1GD == Q1 or grouped_full_lo(Q1GD - 1) < THRESH),
      "lo %.10f" % (math.floor(grouped_full_lo(Q1GD) * 10 ** 10) / 10 ** 10))
check("grouped bound locally increasing at the least q1",
      all(grouped_full_lo(q) <= grouped_full_lo(q + 1) for q in range(Q1GD - 2, Q1GD + 3)))
# ---- 5. Grouped early abort: PHI[g] = Pr[T > g - 1] upper bounds via the per-group
# Chen-Stein prefix bounds (success.md 6.6); WG[g] is E-priced, the rest worst-case.
TPRE = [success_prefix(m, Q1G[:g]) for g, m in enumerate(MQ1[:-1])]
PHI = [Fraction(1)] + [min(Fraction(1), max(Fraction(0), 1 - t["lo"])) for t in TPRE[1:]]
TFULL = success_prefix(Q1GD, Q1G)
check("full grouped mu matches the same-S base form; lo carries the per-group rows",
      TFULL["mu"] == Fraction(SG_FULL * SG_FULL - NODES * N0G * N0G, 2) * P1 and TFULL["lo"] >= THRESH)
check("non-stopping bounds: PHI[0] = 1, nonincreasing, within [0, 1]",
      PHI[0] == 1 and all(0 <= p <= 1 for p in PHI) and all(PHI[g] >= PHI[g + 1] for g in range(G_GROUPS - 1)),
      "PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f" % (float(PHI[1]), float(PHI[16]), float(PHI[31])))
check("per-group hop-1 slots %d are a multiple of k (keep lemma)" % CAP_PG, CAP_PG % K == 0)


check("a hash chain holds one home's records of one fold: at most one slot cap per slot (2^36 slots)",
      HOME_MAX == SIDE * CAP_P, "2^%.4f" % lg(HOME_MAX))
check("uncapped walks always terminate: every chain is finite (at most HOME_MAX entries)",
      HOME_MAX < 2 ** 59)


# ---- E-ledger (success.md 6.6): E_TOTAL = worst-case rest + sum PHI[g]*WG[g] + global
# same-X, K2-verify, listing and defective-run rows; WGW[g] is the true worst case for
# R_WORST (uncapped walks to HOME_MAX, every record verifying and listing). A walk of v
# visits costs 6v plus one B when v is even and no match ends it, so VISIT = 6.5 per
# expected visit; EK2 charges 5 per K2 hit, ELIST 4 per listing. Main-pair walks are
# bounded by history-conditional domination with a validity split (KB129 + 1/r per
# pair, uniformly in the history), same-X pairs are charged globally (SXV), defective
# inserts globally at HOME_MAX visits (DEFV), defective K2/listing pairs in DEFK2/DEFL.
PAIRM = KB129 + Fraction(1, r)
VISIT = HCHG_VISIT + Fraction(HCHG_END, 2)
# The j = 1, 2 key paths (JPATH more per leaf on JH or JT) are E-priced per group like the walks (success.md 6.6,
# key-path means): a valid batch's leaf is uniform given the history, so P(JH or JT, valid | .) <= PJ; defective
# batches take at most JPATH on every leaf, charged globally in DEFJ; the worst case charges every leaf JPATH
# (there is no key-path meter: JH and JT only compute N' and branch back).
PJ = Fraction(2, 3) + 2 * T["d_wj"] + Fraction(1, r - 1)
JROW_W = NODES * N0G * JPATH
check("E-charge per visit 6.5: 6 per visit plus the B after an even visit count, at most one per two visits",
      VISIT == Fraction(13, 2))
SG = [NODES * LEAVES * q * Q2 for q in Q1G]
WG, WGW, SAMEX, WGJ = [], [], [], []
for _g, _q in enumerate(Q1G):
    _s = SG[_g]
    _leaf = NODES * _q * (FWD_INC * Q2 + FWD_FIX + INV + INV_CALL + BATCH_SETUP + GEN_MULX) \
        + int(NODES * _q * (BWD_INC * Q2 + BWD_FIX)) + int(NODES * _q * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX))
    _route = _s * (EMIT * K + BLOCK) // K + _s * (H2S_REC * K + LOOP8 * K // 8 + BLOCK) // K \
        + NODES * (H2S_SLOT * SIDE + H2S_FIX) + NODES * 2 * (FINAL * SIDE + FIN_FIX) \
        + NODES * (INIT3_BKT * SIDE + 1 + TINIT_ENTRY * SIDE + 1)
    _gp = NODES * LEAVES * MQ1[_g] * Q2
    _samex = NODES * _q * (6 * Q2) * (6 * Q2 - 1) // 2
    _evis = (_s * _gp + _s * (_s - 1) // 2 - _samex) * PAIRM
    _gfix = _s * HCHG_NEW + _s * LOOP8 // 8 + NODES * (HCHG_SLOT * SIDE + HCHG_FIXG) + NODES * STOP_CHECK
    _group = _gfix + VISIT * _evis
    assert isinstance(_group, (int, Fraction)) and _group >= 0
    _jkey = JPATH * _s * PJ
    WGJ.append(_jkey)
    WG.append(_leaf + _route + _group + _jkey)
    SAMEX.append(_samex)
    _wgroup = _gfix + _s * HOME_MAX * (HCHG_VISIT + HCHG_END + HCHG_VER + HCHG_LS1)
    assert isinstance(_wgroup, int)
    WGW.append(_leaf + _route + _wgroup)
check("group sizes sum to the grouped q1, each >= 1; group records multiples of k",
      sum(Q1G) == Q1GD and min(Q1G) > 0 and all(s % K == 0 for s in SG))
# The ladder's shape is this package's change (changes.md).
_UNIF = [Q1GD // G_GROUPS + (1 if _g < Q1GD % G_GROUPS else 0) for _g in range(G_GROUPS)]
check("ladder: rung g starts at round(q1*(g/%d)^%.2f), %d rungs, each >= 1, rung records multiples of k; "
      "first rung %d > uniform %d and last rung %d < uniform %d"
      % (G_GROUPS, LADDER_GAMMA, G_GROUPS, Q1G[0], _UNIF[0], Q1G[-1], _UNIF[-1]),
      Q1G[0] > _UNIF[0] and Q1G[-1] < _UNIF[-1] and all(s >= 1 for s in Q1G)
      and all(_q == sum(Q1G[:i + 1]) - sum(Q1G[:i]) for i, _q in enumerate(Q1G)))
check("ladder boundaries strictly increasing and PHI still nonincreasing on it",
      all(MQ1[i] < MQ1[i + 1] for i in range(G_GROUPS))
      and all(PHI[i] >= PHI[i + 1] for i in range(G_GROUPS - 1)),
      "2^%.4f first rung, 2^%.4f last rung" % (lg(SG[0]), lg(SG[-1])))


ZERO_INIT = NODES * (2 ** HB // ZW * ZERO_PASS + 1)
check("hash-head zeroing: zero4096 zeroes 4096 words per pass with 2048 STPs, SUBS and B.NE (2050); 2^57 is a "
      "multiple of 4096, so the row is 2^57/4096 passes plus the RET", ZERO_PASS == 2050 and 2 ** HB % ZW == 0
      and ZERO_INIT == NODES * (2 ** (HB - 12) * 2050 + 1))
HR_INIT = NODES * (INIT3_BKT * SIDE + 1)
E_PRE = NODES * NODE_INIT + NODES * (Q1GD + Q2 + 1) * COIN_SCALAR + NODES * Q2 * GEN_MULY + NODES * HEADER \
    + NODES * (INIT3_BKT * SIDE + 1 + TINIT_ENTRY * SIDE + 1) + ZERO_INIT + HR_INIT
E_WORST_REST = E_PRE + NODES * GEN_MUL + 2 * (SIDE + NODES) * MSG \
    + NODES * DELTA_F * max(FIX_EXTRA.values()) \
    + NODES * CAND_CAP * RESOLVE + NODES * REQ_CAP * REQ_SERVE + NODES * G_GROUPS * 4 * EXCH \
    + NBAR_ALLOW * BARRIER + NODES * PHASE_CTRL_G + OUTPUT + G_GROUPS * 2 * (SIDE + NODES) * MSG
_KB182 = KB(KEYBITS + D1 + D2)
_SIPT = RHO * Fraction(5, 2) * (1 + 11 * (T["d_wj"] + T["d_lab"]))
# v6: global listings sum over homes; a 182-agreeing pair sits at exactly one home
# (the 182 bits include d1, d2), so the leading term takes _KB182 alone. The extra
# KB(D1 + D2) was a fixed-home restriction, legitimate only for 6.4's other_h.
ELIST = RHO * (Fraction(SG_FULL) * SG_FULL * _KB182
               + (2 * SG_FULL + PARTNER_PAIRS) * _KB182) + _SIPT * _KB182 * SG_FULL
EK2 = RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) * KB64
             + (2 * SG_FULL + PARTNER_PAIRS) * KB64) + _SIPT * KB64 * SG_FULL
SXV = VISIT * PAIRM * sum(SAMEX)
BRACKET = RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) + 2 * SG_FULL + PARTNER_PAIRS) \
    + _SIPT * SG_FULL
DEFV = VISIT * T["d_batch"] * (HOME_MAX * SG_FULL + sum(SAMEX))
# v6: a defective record's own walk also verifies/lists on K2/full matches
# (at most HOME_MAX each); the pair terms cover valid-side visits to defectives.
DEFK2 = HCHG_VER * 2 * T["d_batch"] * BRACKET + HCHG_VER * T["d_batch"] * HOME_MAX * SG_FULL
DEFL = HCHG_LS1 * 2 * T["d_batch"] * BRACKET + HCHG_LS1 * T["d_batch"] * HOME_MAX * SG_FULL
# key paths of defective batches, at most JPATH per leaf of the run (unconditional d_batch per batch)
DEFJ = JPATH * T["d_batch"] * SG_FULL
E_TOTAL = E_WORST_REST + sum(PHI[_g] * WG[_g] for _g in range(G_GROUPS)) + SXV + EK2 * HCHG_VER \
    + ELIST * HCHG_LS1 + DEFV + DEFK2 + DEFL + DEFJ
R_WORST = E_WORST_REST + JROW_W + sum(WGW)


def claim6_of(x):
    c = Fraction(math.ceil(lg(x) * 10 ** 6), 10 ** 6)
    for _ in range(50):
        if pow2_interval(c)[0] >= x and pow2_interval(c - Fraction(1, 10 ** 6))[1] < x:
            return c
        c += Fraction(1, 10 ** 6) if pow2_interval(c)[0] < x else -Fraction(1, 10 ** 6)
    raise AssertionError("worst claim not certified")


WORST6 = claim6_of(R_WORST)
check("disclosed worst case 2^%d/10^6: least 6-decimal bound >= R_WORST" % (WORST6 * 10 ** 6),
      pow2_interval(WORST6)[0] >= R_WORST and pow2_interval(WORST6 - Fraction(1, 10 ** 6))[1] < R_WORST)
W_E = E_TOTAL / C


def claim15_of(w):
    # the least multiple of 10^-15 with 2^c >= w, certified both ways by pow2_interval; the search starts from an
    # 80-digit decimal log2 (a double-precision log2 of these huge fractions is off by ~10^-11, too coarse here)
    getcontext().prec = 80
    start = (Decimal(w.numerator).ln() - Decimal(w.denominator).ln()) / Decimal(2).ln() * 10 ** 15
    c = Fraction(int(start.to_integral_value(rounding=ROUND_CEILING)), 10 ** 15)
    for _ in range(400):
        if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 15))[1] < w:
            return c
        c += Fraction(1, 10 ** 15) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 15)
    raise AssertionError("E-claim not certified")


CLAIM_E = claim15_of(W_E)
check("E-work below the base total, worst-case above it", E_TOTAL < R_TOTAL < R_WORST,
      "E 2^%.6f base 2^%.6f worst 2^%.6f (certified 6-decimal bound)" % (lg(E_TOTAL), lg(R_TOTAL), up6(R_WORST)))
check("E-claim = least 15-decimal c with 2^c >= E/1000", pow2_interval(CLAIM_E)[0] >= W_E
      and pow2_interval(CLAIM_E - Fraction(1, 10 ** 15))[1] < W_E, "%d/10^15" % (CLAIM_E * 10 ** 15))
JROW_E = sum(PHI[_g] * WGJ[_g] for _g in range(G_GROUPS)) + DEFJ
check("key paths E-priced: PJ = 2/3 + 2 d_wj + 1/(r - 1) per leaf, per group under PHI, plus DEFJ; below every-leaf row",
      PJ > Fraction(2, 3) and T["d_wj"] > Fraction(1, r - 1) and JROW_E < JROW_W
      and all(WGJ[_g] == JPATH * SG[_g] * PJ for _g in range(G_GROUPS)) and DEFJ == JPATH * T["d_batch"] * SG_FULL,
      "2^%.4f, every-leaf row 2^%.4f" % (lg(JROW_E), lg(JROW_W)))
check("worst case charges every leaf's key path: R_WORST = E_WORST_REST + n0 * JPATH per node + sum WGW",
      R_WORST - JROW_W == E_WORST_REST + sum(WGW) and JROW_W == NODES * N0G * JPATH)
check("barrier allowance covers the grouped run", NBAR_ALLOW >= NBAR_COUNT, "%d of %d" % (NBAR_COUNT, NBAR_ALLOW))
PER_LEAF_E = E_TOTAL / S


CAP_3 = T["cap_3"]
SLOTS = SIDE * CAP_P * 32
FIXED = [("code, constants, stack", 2 ** 20), ("barrier words: 5 granules of at most 2 KiB", 5 * 2048),
         ("header H and the rows of beta below it", 1024), ("X table: H's rows and values, -T(a), -T(beta a), y^2 per X", Q1 * 512),
         ("Y table: y^2, y, b, b^3, 3b^2, b^4", Q2 * 192), ("scalars a, b, c", GENS * 32), ("batch entries", Q2 * 256),
         ("the ports' message blocks", 2 * 64), ("LC tables at 0 and 2^39", 2 * SIDE * 8),
         ("RC tables at 2^43 and 2^43 + 2^39", 2 * SIDE * 8), ("record blocks at 2^50 and 2^51",
          2 * SIDE * 2 ** 14),
         ("final-cursor tables E1 and E2", 2 * SIDE * 32), ("d3 cursor table at 2^40", 2 ** D3 * 8),
         ("hash-read cursors at 2^48", SIDE * 8),
         ("d4 counters at 2^49", 2 ** D4 * 8), ("candidate list", CAND_CAP * 16),
         ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64), ("relay buffer", RELAY_CAP * 64),
         ("phase-4 count, offset and cursor tables", 3 * SIDE * 32)]
FIXED_BYTES = sum(b for _, b in FIXED)
BUCKET_STRIDE = 2 ** 38
PHASE3 = [("d3 buckets: 2^24 strides of 2^38 bytes, cap_3 entries (address, d4 counter) at the end of each",
           2 ** D3 * BUCKET_STRIDE),
          ("pointer array", 2 ** D3 * CAP_3 * 8)]


SLOTS_G = 2 ** D1 * CAP_PG * 32
PHASE3G = [("per-group hop-1 slots", SLOTS_G),
           ("hash heads: 2^57 zeroed words", 2 ** HB * 8),
           ("hash next-array: one 8-byte next word per R2 slot position", SIDE * CAP_P * 8)]
PHASE3G_BYTES = sum(b for _, b in PHASE3G)
# v6: each node retains all groups' X affine points (x, y) for request serving.
AFFINE_BYTES = Q1GD * 64
FIXEDG_BYTES = FIXED_BYTES - Q1 * 512 + max(Q1G) * 512 + AFFINE_BYTES
check("grouped R1 partition (slots + hash heads + hash next-array) fits in R1's range; the next-array "
      "covers every receive pattern up to the slot caps by construction", PHASE3G_BYTES <= SLOTS,
      "2^%.4f" % lg(PHASE3G_BYTES))
R1_BASE = 2 ** 52
R2_BASE = R1_BASE + SLOTS
NEXT_BASE, HEAD_BASE = 2 ** 61, 2 ** 62
NEXT_C = NEXT_BASE - (R2_BASE >> 2)
check("R1 layout: hop-1 slots at 2^52, the next-array at 2^61 (word (A>>2) + C per R2 address A, C = 2^61 - R2/4), "
      "the heads at 2^62 (2^60 bytes: the head address is 2^62 with the bucket in bits 3..59, one BFI); disjoint, in R1",
      R2_BASE % 4 == 0 and R1_BASE + SLOTS_G <= NEXT_BASE and NEXT_BASE + (SLOTS >> 2) <= HEAD_BASE
      and HEAD_BASE + 2 ** HB * 8 <= R2_BASE and HEAD_BASE % 2 ** 60 == 0 and 2 ** HB * 8 == 2 ** 60
      and (SLOTS >> 2) == PHASE3G[2][1], "next-array end 2^%.4f, heads end 2^%.4f, R1 end 2^%.4f"
      % (lg(NEXT_BASE + (SLOTS >> 2)), lg(HEAD_BASE + 2 ** HB * 8), lg(R2_BASE)))
check("grouped fixed regions (per-group X table, all-groups X affine store, no snapshot buffer) lie below 2^52",
      FIXEDG_BYTES <= 2 ** 52, "2^%.4f" % lg(FIXEDG_BYTES))
PHASE3_BYTES = sum(b for _, b in PHASE3)
for name, b in FIXED + [("R1 at 2^52: hop-1 pair slots (phases 1-2)", SLOTS),
                        ("R2 at 2^52 + |R1|: hop-2 pair slots (phases 2-4)", SLOTS)] + PHASE3:
    print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("fixed regions lie below 2^52", FIXED_BYTES <= 2 ** 52, "2^%.4f" % lg(FIXED_BYTES))
ALIGNED = 4 * SIDE * 8 + 2 ** D3 * 8 + 2 ** D4 * 8 + 2 * SIDE * 2 ** 14 + SIDE * 8
check("aligned tables: LC at 0 and 2^39, d3 at 2^40, RC at 2^43, hash-read cursors at 2^48, d4 at 2^49, column "
      "blocks at 2^50, row blocks at 2^51; the rest fits between 2^44 and 2^48",
      SIDE * 8 == 2 ** 39 and 2 ** D3 * 8 == 2 ** 27 and 2 ** D4 * 8 == 2 ** 49 and SIDE * 2 ** 14 == 2 ** 50
      and D4BASE == 1 << 49 and D4BASE + 2 ** D4 * 8 == PB1C and PB1C == 1 << 50 and PB2C == 1 << 51
      and PB2C + (1 << 50) <= 1 << 52 and FIXED_BYTES - ALIGNED <= (1 << 48) - (1 << 44),
      "2^%.4f" % lg(FIXED_BYTES - ALIGNED))
check("record block: header and k slots within 2^14; bit 13 set exactly at slot k; bases/stride keep low 13 bits zero; "
      "cap_p a multiple of k (keep lemma)",
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
check("hop-1 cursors cannot wrap: the last slot's base plus 32 per leaf of a node stays below 2^64", 2 ** 52 + 32 * CAP_P * (SIDE - 1) + 32 * N0 < 2 ** 64,
      "2^%.4f" % lg(2 ** 52 + 32 * CAP_P * (SIDE - 1) + 32 * N0))
NODE_BYTES = FIXED_BYTES + 2 * SLOTS
PEAK = NODES * NODE_BYTES
check("peak memory 2^72 * (fixed + R1 + R2) = 2^%.4f bytes" % lg(PEAK), PEAK < 2 ** 136)
check("hash-read cursors occupy [2^48, 2^48 + 2^39); d4 at 2^49",
      (1 << 48) + SIDE * 8 <= 1 << 49 and D4BASE == 1 << 49)
PEAK_G = NODES * (FIXEDG_BYTES + 2 * SLOTS)
check("grouped peak memory = 2^%.4f bytes" % lg(PEAK_G), PEAK_G < 2 ** 136)
check("caps exceed crude means: cap_p > m_p + 2^20, cap_3 > m_3 + 2^31",
      CAP_P > T["m_p"] + 2 ** 20 and CAP_3 > T["m_3"] + 2 ** 31)
check("preprocessing (rows marked P) = 2^%.7f units, included in the total" % lg(Fraction(R_PRE, C)), R_PRE < R_TOTAL)
check("grouped preprocessing (hash-head zeroing, read-cursor init) = 2^%.7f units, included"
      % lg(Fraction(E_PRE, C)), E_PRE < R_WORST)


DECLARED = dict(time_log2=CLAIM_E, memory_log2_bytes=None, preprocessing_log2=None,
                success=Fraction(390000001, 10 ** 9),
                per_leaf=Fraction(cdiv(PER_LEAF_E.numerator * 10 ** 4, PER_LEAF_E.denominator), 10 ** 4))
BASE_DECLARED_TIME = Fraction(124903563415332, 10 ** 12)


def least_at_or_above(d, x, scale):
    return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x


def least_grid_up(x, scale):
    k = math.ceil(lg(x) * scale)
    for _ in range(20):
        d = Fraction(k, scale)
        if least_at_or_above(d, x, scale):
            return d
        k += 1 if pow2_interval(d)[0] < x else -1
    raise AssertionError("grid bound not certified")


DECLARED["memory_log2_bytes"] = least_grid_up(PEAK_G, 10 ** 4)
DECLARED["preprocessing_log2"] = least_grid_up(Fraction(E_PRE, C), 10 ** 4)
check("base claim intact: %d/10^12 least with 2^c >= R/1000" % (BASE_DECLARED_TIME * 10 ** 12),
      BASE_DECLARED_TIME == CLAIM and least_at_or_above(BASE_DECLARED_TIME, W, 10 ** 12))
check("declared E time_log2 %d/10^15: least 15-decimal value with 2^c >= E/1000" % (DECLARED["time_log2"] * 10 ** 15),
      DECLARED["time_log2"] == CLAIM_E and least_at_or_above(DECLARED["time_log2"], W_E, 10 ** 15))
check("declared memory 2^%s bytes: least 4-decimal bound >= grouped peak (2^%.7f)"
      % (format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK_G)),
      least_at_or_above(DECLARED["memory_log2_bytes"], PEAK_G, 10 ** 4))
check("declared preprocessing 2^%s units: least 4-decimal bound >= grouped rows marked P / 1000"
      % format(float(DECLARED["preprocessing_log2"]), ".4f"),
      least_at_or_above(DECLARED["preprocessing_log2"], Fraction(E_PRE, C), 10 ** 4))
check("declared success 0.390000001 <= the bound", DECLARED["success"] <= TFULL["lo"])
check("declared E work per leaf %s: least 4-decimal bound >= E/S (%.7f)"
      % (format(float(DECLARED["per_leaf"]), ".4f"), float(PER_LEAF_E)),
      DECLARED["per_leaf"] >= PER_LEAF_E > DECLARED["per_leaf"] - Fraction(1, 10 ** 4))
ANALYTIC = sum(v for n, v, _ in rows if n.split(":")[0] in ("program delivery, node and port state", "coins",
               "Y generators [c]P with their table entries", "R_g = [b]Q and X = R_g + [a]P with their header values",
               "rows of beta", "batch set-up",
               "Q broadcast", "candidates", "requests served", "phase-4 count exchanges", "barriers", "phase control",
               "output",
               "inversion")) - NODES * Q1 * INV


MULX, SQRX, ADDX = 105 + 4, 80 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX

GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_TEST = 256 * 2 * (9 + 4) + 64

GEN_CTRL = 1600
NEGT = 30


RED = 9
Y_ENTRY = 3 * SQRX + 2 * MULX + 2 * ADDX + 2 * NEGT + 7 * RED + 14

X_ENTRY = MADD_ + 3 * SQRX + 11 * MULX + 2 * ADDX + 2 * NEGT + 17 * RED + 32


check("generator (exact arithmetic): doubling 933, mixed addition 1,477; with the exceptional-case tests, control, "
      "the Y entry and the X header values at most %d, the allowance" % (GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY,
      X_ENTRY)), DBL == 903 and MADD_ == 1391 and GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY, X_ENTRY) <= GEN_MUL,
      "%.4fx" % (math.floor(10 ** 4 * GEN_MUL / (GEN_ITEM + GEN_TEST + GEN_CTRL + max(Y_ENTRY, X_ENTRY))) / 10 ** 4))
check("v6 grouped X generator: header values plus the four affine-store STPs (64 bytes) fit GEN_MULX",
      GEN_ITEM + GEN_TEST + GEN_CTRL + X_ENTRY + 4 <= GEN_MULX, "%d <= %d"
      % (GEN_ITEM + GEN_TEST + GEN_CTRL + X_ENTRY + 4, GEN_MULX))
check("Y generator: [c]P with its table entry (no mixed addition, no X header) fits GEN_MULY",
      GEN_ITEM + GEN_TEST + GEN_CTRL + Y_ENTRY <= GEN_MULY and Y_ENTRY < X_ENTRY, "%d <= %d"
      % (GEN_ITEM + GEN_TEST + GEN_CTRL + Y_ENTRY, GEN_MULY))


BARRIER_ITEM = 5 * 8 + 5 * (6 + MSG_LINK) + 1 + 20


BATCH_ITEM = (2 + 4) + 60 + (1 + 4 + 4 + 4 + 2 + 1) + (2 * 4 + 1) + (4 + 1 + 2 + 1 + 1 + 1) + 2 * (4 + 1) \
    + 2 * (2 + 1 + 3) + 2 * (1 + 4 + 1 + 1) + (2 + 2 + 1 + 1) + 2
check("batch set-up: at most %d instructions, within its allowance of 160" % BATCH_ITEM, BATCH_ITEM <= BATCH_SETUP)
check("grouped phase control: 12 kernel calls per group at most 33 each, stop test and loop 8, over %d groups < 2^16"
      % G_GROUPS, G_GROUPS * (12 * 33 + 8) < PHASE_CTRL_G)
check("barrier: at most 131 instructions per node, in-flight counts and deliveries included (2^10 per node is > 7x)",
      BARRIER_ITEM == 131 and BARRIER > 7 * NODES * BARRIER_ITEM)


EXCH_ITEM = 12 * SIDE + 4 * (14 + MSG_LINK) * (8 * REQ_CAP)
check("phase-4 count exchange: < 2^40.5 per node, within 2^41",
      EXCH_ITEM < pow2_interval(Fraction(405, 10))[0] and EXCH_ITEM < EXCH, "2^%.3f" % lg(EXCH_ITEM))


check("a broadcast message costs at most 16, in-flight count and delivery included, within its price of 40",
      2 + 3 + 3 + MSG_LINK <= MSG)


CAND_LEAVES = INV + 5 * SQRX + 18 * MULX + 20 * ADDX
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 9 * 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, both records' three candidate leaves, nine unit tests, solve mod r, [x]P = Q = %d; "
      "charged above twice that" % CAND_ITEM, RESOLVE == 2 * CAND_ITEM + 1,
      "%.3fx" % (math.floor(10 ** 3 * RESOLVE / CAND_ITEM) / 10 ** 3))
PROGRAM_BYTES = 2 ** 16
check("program delivery: 2^16 bytes = 2^11 messages over 2 hops at 40 = 2^17.3 < node set-up 2^20",
      PROGRAM_BYTES // 32 * 2 * MSG < NODE_INIT)
check("analytic allowances are below 2^-13.7 of the total", ANALYTIC * pow2_interval(Fraction(137, 10))[1] < R_TOTAL,
      "2^%.3f" % lg(Fraction(ANALYTIC, R_TOTAL)))
ANALYTIC_E = ANALYTIC + NBAR_ALLOW * BARRIER + G_GROUPS * 2 * (SIDE + NODES) * MSG
check("grouped analytic allowances are below 2^-13 of E-work",
      ANALYTIC_E * pow2_interval(Fraction(13, 1))[1] < E_TOTAL, "2^%.3f" % lg(Fraction(ANALYTIC_E, E_TOTAL)))

print()
print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1, N0, lg(S), float(T["mu"]),
                                                                    math.floor(T["lo"] * 10 ** 10) / 10 ** 10))
print("R_total %d instructions  W 2^%.10f  claim %d/10^12  per leaf %.4f" % (R_TOTAL, lg(W), CLAIM * 10 ** 12,
                                                                      float(PER_LEAF)))
print("E_total 2^%.6f instructions  W_E 2^%.10f  E-claim %d/10^15  per leaf %.7f declared %s  R_worst 2^%.6f"
      % (lg(E_TOTAL), lg(W_E), CLAIM_E * 10 ** 15, float(PER_LEAF_E), format(float(DECLARED["per_leaf"]), ".4f"),
         up6(R_WORST)))
print("key paths: E-priced 2^%.6f (every-leaf row 2^%.6f in the worst case, JPATH %d); Y items at %d (was %d); hchg "
      "%d fixed, %s per visit" % (lg(JROW_E), lg(JROW_W), JPATH, GEN_MULY, GEN_MUL, HCHG_NEW, VISIT))
print("declared: memory_log2_bytes %s (grouped exact %.7f)  preprocessing_log2 %s (grouped exact %.7f)" % (
    format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK_G),
    format(float(DECLARED["preprocessing_log2"]), ".4f"), lg(Fraction(E_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
    sys.exit(1 if FAILED else 0)
