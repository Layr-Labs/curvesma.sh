"""Exact certificate for proof.md (parameters, E-work ledger, claim, memory)."""
import functools
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

@functools.lru_cache(maxsize=None)  # pure: exact interval for 2^c, memoized (the run checks some grid points twice)
def pow2_interval(c):
 n = c.numerator // c.denominator
 f = c - n
 lo, hi = exp_interval(f * LN2[0])[0], exp_interval(f * LN2[1])[1]
 sc = Fraction(2) ** n
 return lo * sc, hi * sc

# 1. constants
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
C0 = 2 ** 32 + 977
# 2. success (grouped q1); 5. grouped early abort (2.3, 4, 6.6): E-work ledger; 4. memory (bytes per node)
BETA_P = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
C = 1000
NODES = 2 ** 72
SIDE = 2 ** 36
Q2 = 2 ** 25
D1, D2, D4 = 36, 36, 46
KEYBITS = 110
THRESH = Fraction(39, 100)  # the rules' minimum (as aburan28's 7516608f)

LEAVES = 6
FWD_U = 32  # fwd's main loop: 32 incidences per pass
FWD_INC, FWD_FIX = 97 + Fraction(2, FWD_U), 9

def fwd_count(q2):  # 32/pass + 2/pass tail; 97 each + 2/pass; 9 fixed
 return FWD_FIX + 97 * q2 + 2 * (q2 // FWD_U) + (q2 % FWD_U)

BWD_INC, BWD_FIX = Fraction(4 * 1230 + 1, 4), 5  # four 1230-instruction bodies per TBZ of bit 34 (q2 a multiple of 4)

FIX_SITES = dict(SHX=18, FDX=8, F5=2, FS=2, FADDW=1, BLK=2)
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4, BLK=5)
EMIT = 6
EMIT_INC = LEAVES * EMIT
KEY_INC, KEY_FIX = Fraction(277, 2), 4
KEY_SIGN = 1  # per sign: LSL x18 = x3<<30, hoisted out of the three emits (j = 1 takes EON x18,S,LSR #28)

K = 256
PB1C, PB2C, BLKC = 1 << 50, 1 << 51, 1 << 14
FLUSH = 9

KEY_BLOCK_DROP, H2S_BLOCK_DROP = 7, 9
LOOP8 = 2

JPATH = 1
SJOIN = 1  # key-path tree (C2, C3): no branch back; a sign with a j = 1 or 2 leaf ends in one B, charged on every sign

INV = 22549
TINIT_ENTRY = 3
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 14, 4, 2, 2

# h2s (FEAT_LS64): per record 6 (H2B: BFI, LDADD, 2 STP, TBZ; slots in pairs by one LD64B and one ADD); a slot read
# alone (slot 1, a last even slot) takes two LDPs, +1
H2S_ARR, H2S_ONE = 6, 1
H2S_FULL = 6 + 2 * H2S_ONE  # a full block (count 8192): MOVN x1, LDR count, LDR p, MOV x17, TBNZ of bit 13, RET; slots 1, 256 alone, no loop
H2S_EMPTY = 7  # count 0: MOVN, LDR, LDR, MOV, TBNZ, CBZ, RET
H2S_ARR_FIX = 14 + 2 * H2S_ONE  # other blocks: + CBZ, SUB, ADD, TBZ x3, CMP, B.HS; slot 1 and at most one last slot alone

def h2s_count(n):  # executed h2s instructions on an n-record block, 0 <= n <= K
 if n == K:
  return H2S_FULL + H2S_ARR * K
 if n == 0:
  return H2S_EMPTY
 return 14 + H2S_ONE * (1 + (n % 2 == 0)) + LOOP8 * ((n - 1) // 8) + H2S_ARR * n
check("h2s: 6 per record (one LD64B and one ADD per two slots), +1 per slot read alone by two LDPs; a full block 8 fixed, any "
   "other at most 16 + 2 per 8-record pass (n - 1 records after slot 1, at most k/8 - 1 passes): <= 16 + 62",
   h2s_count(K) == 6 * K + 8 and max(h2s_count(n) - H2S_ARR * n for n in range(K)) == H2S_ARR_FIX + LOOP8 * (K // 8 - 1)
   and H2S_EMPTY <= H2S_ARR_FIX)

INIT3_BKT = 4

check("key pass: 138.5 = 2 signs * (6 LDP of the fractions + 26 sums (13 two-word) + the hoisted LSL + 3 leaves * (index ADDS, "
   "ADC, CMP, B.HS, CMN, B.HS + emit 6)) + TBZ per two incidences",
   KEY_INC == 2 * (6 + 26 + KEY_SIGN + 3 * (6 + EMIT)) + Fraction(1, 2))
check("bwd: four BWDI bodies (1230 each) per TBZ of bit 34, 5 fixed; q2 a multiple of 4",
   BWD_INC == 1230 + Fraction(1, 4) and Q2 % 4 == 0 and BWD_FIX == 5)
check("fwd: 97 each + 2 per 32-pass, 9 fixed; Q2 a multiple of 32",
   Q2 % FWD_U == 0 and fwd_count(Q2) == FWD_FIX + FWD_INC * Q2 and fwd_count(2) == 9 + 2 * 97 + 2)
check("j = 1 and j = 2 paths (ADDS, ADC, CMP, B.HS; CMN, B.HS; MVN of sigma for N'; emit) exceed the j = 0 path "
   "(ADDS, ADC, CMP, B.HS, CMN, B.HS; emit) by 1 (no branch back; tree copies end in one B per sign: SJOIN)",
   (4 + 2 + 1 + EMIT) - (6 + EMIT) == JPATH and SJOIN == 1)
check("emit: AND or BIC (8*d1), LDADD of its LC entry (the slot), w3 in one ADD x18,W,LSR #28 (j = 0, 2) or EON "
   "x18,S,LSR #28 (j = 1: d2 of N' = ~S, c and l' complemented, c in {2, 4} -> {5, 3}) on the per-sign x18 = x3<<30, two STPs, TBZ of bit 13 for slot k",
   EMIT == 1 + 1 + 1 + 2 + 1)
DONE = 2
DLV = 6
DLVK = 4 * K + 2  # part B's reference delivery only
MSG_LINK = DONE + DLV
BLOCK = FLUSH + DONE  # hop-1 block; arrival by h2s
FINAL = FIN_DEST + FIN_SENT + DONE + MSG_LINK  # hop-1 fin; final block by h2s
check("a hop-1 block: flush 9 + in-flight count 2 when sent, more than 7 when dropped; hop-1 fin per destination at most "
   "28 (14, send 4, in-flight 2, E1 message 8), its final block consumed by h2s (h2s rows), no dlvk",
   BLOCK == 11 and BLOCK > KEY_BLOCK_DROP and FINAL == 28 and FINAL > FIN_DEST + FIN_HELD + MSG_LINK)

GEN_MUL = 645410
# grouped X also stores (x, y): +64 bytes, +4 STPs.
GEN_MULX = GEN_MUL + 4
# Y item: [c]P + table entry only.
GEN_MULY = 642930

COIN_SCALAR = 64
BATCH_SETUP = 160

INV_CALL = 16
NODE_INIT = 2 ** 20
HEADER = 72
BT_BASE, BT_W = 0xfff8 << 48, 43
TBLD = 4 + 6 * 12 + 20 * (5 * 2 ** BT_W + 2 ** (256 - 5 * BT_W))  # tbld: 20 per entry, 6 slots (route.s)
MSG = 40

CAND_CAP = 2 ** 28
RESOLVE = 1491011

REQ_CAP = 2 ** 29
REQ_SERVE = 2 ** 10
RELAY_CAP = 2 ** 34
EXCH_G = 12 * SIDE + 4 * (14 + MSG_LINK) * (8 * REQ_CAP)  # charged at its itemization (7)
BARRIER = 2 ** 82
OUTPUT = 2 ** 20
DELTA_F = 2 ** 62
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
]
SITES_PER_INCIDENCE = dict(SHX=18, FDX=8, F5=2)
EXACT_AT_L0 = dict(SHX=3, FDX=3)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)

FIX_ROW = C0 * 2 ** 192
FIX_FOLD = (2 ** 66 + 1) * 2 ** 128 * C0

FIX_SQ = 2 ** 193
FIX_BLK = 30 * C0 * 2 ** 128  # fixup iff low 128 bits >= 2^128 - t C0; sum of six entries < 6p: carry iff its low 256 bits >= 2^256 - t C0 (t <= 5), <= 2 w per sum
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
 ("bwd: window-sum folds of beta u^2 and beta v^2", "BLK", 2 * Q2, FIX_BLK, 6),
]
E_FIX = sum(Fraction(n * bad * mult, r - 1) for _, _, n, bad, mult in FIXUPS)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (2 ** 17 + 1), (r - 1) * 2 ** 17)

RHO = Fraction(r, r - 1) ** 4
FIXED_LEAVES = 2 * LEAVES
PARTNER_PAIRS = 2 * 10

D_BATCH = D_SITES + Q2 * Fraction(6, r - 1) + Fraction(1, r)
D_WJ = TWO_BOUND * 23 * 2 ** 192
D_LAB = TWO_BOUND * 1305 * 2 ** 228

FA, FB = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16), int("e4437ed6010e88286f547fa90abfe4c3", 16)
BOX = (p >> 100) + 1
WIDE = 654 * (p >> 128) + 4
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p, parallelogram 2A - B wide "
   "in u and A in v; box factor below 1 + 2^-17, one side's below 1 + 2^-18",
   FA * FA - FA * FB + FB * FB == p and (FB * BETA_P - (FA - FB)) % p == 0 and ((FA - FB) * BETA_P + FA) % p == 0
   and (BOX + WIDE + FA) * (BOX + WIDE + 2 * FA - FB) * 2 ** 200 < p * p * (1 + Fraction(1, 2 ** 17))
   and (BOX + WIDE) * 2 ** 100 < p * (1 + Fraction(1, 2 ** 18)))
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
check("wrong class index rate per leaf < 2^-58.47", D_WJ < pow2_interval(Fraction(-5847, 100))[0], "2^%.4f" % lg(D_WJ))
check("site table: 18 row shifts (3 by SHXT), 8 folds, 2 square folds per incidence, less 3 + 3 exact at l = 0 (the beta window sums are exact)",
   sum(n for _, n, b, _ in SITES if b == SHX_BAD) == SITES_PER_INCIDENCE["SHX"] * Q2 - EXACT_AT_L0["SHX"]
   and sum(n for _, n, b, _ in SITES if b in (FOLD_BAD, CP_BAD)) == SITES_PER_INCIDENCE["FDX"] * Q2
   - EXACT_AT_L0["FDX"]
   and sum(n for _, n, b, _ in SITES if b == SQ_BAD) == SITES_PER_INCIDENCE["F5"] * Q2 and len(SITES) == 10)
check("bad sets: rows 2^97, folds below 2^98.02, square folds below 2^64.01",
   SHX_BAD == 2 ** 97 and all(b < pow2_interval(Fraction(e, 100))[0] for b, e in ((FOLD_BAD, 9802), (CP_BAD, 9802),
                      (SQ_BAD, 6401))))
FIX_RUNS = {m: sum(n for _, ms, n, _, _ in FIXUPS if ms == m) for m in FIX_SITES}
check("fixup table: per incidence 18 row shifts, 8 folds, 2 window-sum folds, 2 square folds, 2 FS and 1 FADDW, less 3 + 3 exact at l = 0",
   FIX_RUNS == {m: FIX_SITES[m] * Q2 - EXACT_AT_L0.get(m, 0) for m in FIX_SITES}
   and all(bad == dict(SHX=FIX_ROW, FDX=FIX_FOLD, F5=FIX_SQ, FS=FIX_SUB, FADDW=FIX_SUB, BLK=FIX_BLK)[m] for _, m, _, bad, _ in FIXUPS)
   and all(FIX_SITES[m] == SITES_PER_INCIDENCE[m] for m in ("SHX", "FDX", "F5")))
check("a call's fixup word (C0 per fixup, at most 33 q2 fixups) stays below 2^63", FIX_BATCH * C0 < 2 ** 63,
   "2^%.3f" % lg(FIX_BATCH * C0))
check("label-error rate per leaf < 2^-16.65", D_LAB < pow2_interval(Fraction(-1665, 100))[0], "2^%.3f" % lg(D_LAB))
check("batch defect rate (site table, D = 0, X = O) < 2^-126.05", D_BATCH < pow2_interval(Fraction(-12605, 100))[0],
   "2^%.3f" % lg(D_BATCH))
G_GROUPS, GAM = 160, (21, 32)  # group g ends at q1*(g/G)^(21/32)
CAP_PG = 2 ** 20 - 2 ** 11  # 32*cap_pg*2^36 = 2^61 - 2^52: R1 below the heads

def group_sizes(q1):  # m_g: the largest m with m^32 G^21 <= q1^32 g^21 (exact); m_0 = 0, m_G = q1
 a, b = GAM
 ms = [0]
 for g in range(1, G_GROUPS + 1):
  t, m = q1 ** b * g ** a, int(q1 * (g / G_GROUPS) ** (a / b)) + 2
  while m ** b * G_GROUPS ** a > t:
   m -= 1
  assert (m + 1) ** b * G_GROUPS ** a > t
  ms.append(m)
 return [ms[g + 1] - ms[g] for g in range(G_GROUPS)], ms
HB = 58
HCHG_NEW = 4
HCHG_VISIT = 4
HCHG_END = 1
HCHG_VER = 5
HCHG_LS0 = 2
HCHG_LSA = 1  # a listing record whose address is not x2 (the second of a pair, a slot alone) forms it: ADD/SUB x12, x2, #32
HCHG_LS1 = 4 + HCHG_LSA
check("join per record (BFI, SWP of the entry address x13, post-indexed STP, CBZ; K2, K3 loaded by dlvh's copy; x13 "
   "steps 16 with the 32-byte slot): 4; per visit on entry addresses (LDP, CMP, B.EQ, CBZ or CBNZ): 4; B after an even "
   "visit count: 1; a K2 hit's d4 check (ADD of A = 2E - 2C, LDR, EOR, TST, B.NE): 5",
   HCHG_NEW == 4 and HCHG_VISIT == 4 and HCHG_END == 1 and HCHG_VER == 5)
# dlvh's copy (FEAT_LS64): two slots per LD64B, ADD, ST64B, ADD (2 per slot); a slot alone (odd first, unpaired last)
# two LDPs and two STPs (4: +2)
DLVC, DLVS1 = 2, 2
DLVH_FAST = 15 + 2 * DLVS1  # LDP, 2 LDR, 2 MOV, LDP, LDP, ADD x1, MOV x2, CMP, B.LO, TBNZ, ADD of x13, STR, RET; slots 1, 256 alone
DLVH_SEG, DLVH_SLOW_FIX = 7, 24  # a non-empty segment: CBZ, TBZ of x1 bit 5, 3 TBZ, CMP, B.HS (+ SUB after an odd first slot); 12, CSEL, 8 bounds, ADD of x13, STR, RET

def dlvh_seg(n, m):  # copy of m slots from slot n (1-based): CBZ alone if empty
 if m == 0:
  return 1
 head = n % 2
 rest = m - head
 return DLVH_SEG + head + DLVC * m + DLVS1 * (head + rest % 2) + LOOP8 * (rest // 8)

def dlvh_copy(m1, m2):  # dlvh's executed copy and fixed instructions, joins excluded: m1 copied, m2 joined, the rest copied
 if m1 == 0 and m2 == K:
  return DLVH_FAST + DLVC * K
 return DLVH_SLOW_FIX - (m1 > 0) + dlvh_seg(1, m1) + dlvh_seg(1 + m1, m2) + dlvh_seg(1 + m1 + m2, K - m1 - m2)
DLVHB = dlvh_copy(0, K)
DLVHS = max(dlvh_copy(a, b) for a in range(K + 1) for b in range(K + 1 - a) if (a, b) != (0, K))
BLOCKH = FLUSH + DONE + DLVHB
FINALH = FIN_DEST + FIN_SENT + DONE + DLVHS + MSG_LINK + (DLVHS - DLVHB)
check("dlvh: copy 2 per slot (LD64B, ADD, ST64B, ADD per two), 4 per slot alone; fresh full block 15 fixed (LDP, 2 LDR, "
   "2 MOV, 2 LDP, ADD, MOV, CMP, B.LO, TBNZ, ADD of x13, STR, RET) and slots 1, 256 alone; else at most 24 (12, CSEL, "
   "8 bounds, ADD of x13, STR, RET) and per segment 1 if empty, else 7, + 1 after an odd first slot, + 2 per 8-slot "
   "pass; DLVHS the maximum over all segment splits; hop-2 final: slow + surcharge",
   DLVHB == 531 and DLVHS == 633 and BLOCKH == 542 and BLOCKH > H2S_BLOCK_DROP and FINALH == FINAL + DLVHS + DLVHS - DLVHB)
ZW = 16384  # words per zero16384 pass: 512 DC ZVA of a 256-byte block (DCZID_EL0.BS = 6), each with ADD #256
ZERO_PASS = ZW // 16 + 2
HOME_MAX = SIDE * CAP_P

STOP_CHECK = 5
NBAR_COUNT = 4 + 17 * G_GROUPS + 1
NBAR_ALLOW = 4096
PHASE_CTRL_G = 2 ** 17

# walks uncapped; a chain holds one home/fold, at most HOME_MAX.
KB130 = KB(D1 + D2 + HB)
KB64 = KB(64)

@functools.lru_cache(maxsize=None)  # a group's hop-1 row (6.4), rounded up to a 2^-256 multiple, and its mean
def p1_row(q):
 m = LEAVES * q * Q2 * KB(D1) + 1
 x = (m + 1 + 60 * (q + Q2) * KB(D1)) / Fraction(CAP_PG - FIXED_LEAVES - m) ** 2
 return Fraction(cdiv(x.numerator << 256, x.denominator), 2 ** 256), m

def success_prefix(q1p, qs):
 n0p = LEAVES * q1p * Q2
 Sp = NODES * n0p
 gamma = Fraction(Sp * Sp - NODES * n0p * n0p, 2)
 mu = gamma * P1
 deg = 2 * LEAVES * (q1p + Q2) * Sp
 cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
 mu_lo = Fraction((mu.numerator << 200) // mu.denominator, 2 ** 200) if mu else Fraction(0)
 e_hi = exp_interval(-mu_lo, 40)[1]
 d_batch, d_wj, d_lab = D_BATCH, D_WJ, D_LAB
 d_f = q1p * E_FIX / (DELTA_F - FIX_BATCH) + q1p * d_batch
 m_p2 = n0p * KB(D1) + 1
 v_p2 = m_p2 + 1 + 60 * (q1p + Q2) * KB(D1)
 d_p2 = v_p2 / Fraction(CAP_P - FIXED_LEAVES - m_p2) ** 2
 d_p1, m_p1 = Fraction(0), Fraction(0)
 for _q in qs:  # one row per group
  _row, _m1 = p1_row(_q)
  d_p1 += _row
  m_p1 = max(m_p1, _m1)
 kb_m = KB(KEYBITS + D1 + D2)
 lam_int = RHO * (Sp + FIXED_LEAVES) * KB(D1 + D2 + KEYBITS)  # the 182 bits a walk compares
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
 return dict(mu=mu, lo=lo, cs=cs, rows=rows_, glob=glob, Sp=Sp, n0p=n0p, d_p1=d_p1, d_p2=d_p2, m_p1=m_p1, m_p2=m_p2,
    d_f=d_f, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def)

def grouped_full_lo(q):
 qs, _ = group_sizes(q)
 return success_prefix(q, qs)["lo"]

lo_q, hi_q = 2 ** 27, 2 ** 28
check("grouped q1 search starts at 2^27, where the full grouped bound fails", grouped_full_lo(lo_q) < THRESH)
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
check("grouped q1 = %d is least with the full grouped bound >= 0.39" % Q1GD,
   grouped_full_lo(Q1GD) >= THRESH and grouped_full_lo(Q1GD - 1) < THRESH,
   "lo %.10f" % (math.floor(grouped_full_lo(Q1GD) * 10 ** 10) / 10 ** 10))
check("grouped bound locally increasing at the least q1",
   all(grouped_full_lo(q) <= grouped_full_lo(q + 1) for q in range(Q1GD - 2, Q1GD + 3)))
# 5. grouped abort: PHI[g] = Pr[T>g-1] bounds (6.6); WG[g] E-priced.
TPRE = [success_prefix(m, Q1G[:g]) for g, m in enumerate(MQ1[:-1])]
PHI = [Fraction(1)] + [min(Fraction(1), max(Fraction(0), 1 - t["lo"])) for t in TPRE[1:]]
TFULL = success_prefix(Q1GD, Q1G)
check("leaves S = 2^72*6*q1*q2 = 2^%.6f" % lg(SG_FULL), SG_FULL == NODES * 6 * Q1GD * Q2 == TFULL["Sp"])
check("mu = |Gamma|*P1 = %.9f (cross-node pairs); lo carries the per-group rows" % float(TFULL["mu"]),
   TFULL["mu"] == Fraction(SG_FULL * SG_FULL - NODES * N0G * N0G, 2) * P1 and TFULL["mu"] > Fraction(494, 1000)
   and TFULL["lo"] >= THRESH)
check("Chen-Stein b1+b2 < 2^-95", TFULL["cs"] < Fraction(1, 2 ** 95), "2^%.2f" % lg(TFULL["cs"]))
check("intruder rate per colliding pair (same 182 bits) < 2^-52.7", TFULL["lam_int"] < pow2_interval(Fraction(-527, 10))[0],
   "2^%.3f" % lg(TFULL["lam_int"]))
check("mean fixups per node (clean batches) q1*E_FIX = 2^%.3f < 2^28.68; cap 2^62 > 2^16 times it plus 33 q2"
   % lg(Q1GD * E_FIX), Q1GD * E_FIX < pow2_interval(Fraction(2868, 100))[0]
   and DELTA_F > 2 ** 16 * Q1GD * E_FIX + FIX_BATCH)
check("fixup meter rate (Markov) < 2^-33.5", TFULL["d_f"] < pow2_interval(Fraction(-335, 10))[0],
   "2^%.3f" % lg(TFULL["d_f"]))
check("pair-slot overflow rate: hop 1, the per-group rows summed, < 2^-18.6; hop 2 < 2^-19.9",
   TFULL["d_p1"] < pow2_interval(Fraction(-186, 10))[0] and TFULL["d_p2"] < pow2_interval(Fraction(-199, 10))[0],
   "2^%.3f 2^%.3f" % (lg(TFULL["d_p1"]), lg(TFULL["d_p2"])))
check("candidate cap rate < 2^-23.4, request cap rate < 2^-25.9",
   TFULL["d_cl"] < pow2_interval(Fraction(-234, 10))[0] and TFULL["d_rq"] < pow2_interval(Fraction(-259, 10))[0],
   "2^%.3f 2^%.3f" % (lg(TFULL["d_cl"]), lg(TFULL["d_rq"])))
check("relay cap rate (a reply counts 8 records) < 2^-25.2", TFULL["d_rl"] < pow2_interval(Fraction(-252, 10))[0],
   "2^%.3f" % lg(TFULL["d_rl"]))
check("whole-run row: a computational defect in any batch < 2^-26.5", TFULL["g_def"] < pow2_interval(Fraction(-265, 10))[0],
   "2^%.3f" % lg(TFULL["g_def"]))
check("whole-run rows (defective batch anywhere, coin rejection twice) < 2^-26.5",
   TFULL["glob"] < pow2_interval(Fraction(-265, 10))[0], "2^%.3f" % lg(TFULL["glob"]))
check("success >= 1 - e^-mu - (b1+b2) - rows - global >= 0.39", TFULL["lo"] >= THRESH,
   "%.10f" % (math.floor(TFULL["lo"] * 10 ** 10) / 10 ** 10))
check("non-stopping bounds: PHI[0] = 1, nonincreasing, within [0, 1]",
   PHI[0] == 1 and all(0 <= p <= 1 for p in PHI) and all(PHI[g] >= PHI[g + 1] for g in range(G_GROUPS - 1)),
   "PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f" % (float(PHI[1]), float(PHI[16]), float(PHI[31])))
check("per-group hop-1 slots %d are a multiple of k (keep lemma)" % CAP_PG, CAP_PG % K == 0)

check("a hash chain holds one home's records of one fold: at most one slot cap per slot (2^36 slots)",
   HOME_MAX == SIDE * CAP_P, "2^%.4f" % lg(HOME_MAX))
check("uncapped walks always terminate: every chain is finite (at most HOME_MAX entries)",
   HOME_MAX < 2 ** 59)

# E-ledger: E_TOTAL = rest + sum PHI[g]*WG[g] + globals; WGW[g] -> R_WORST.
# VISIT 4.5/visit; EK2 5/hit (ADD,LDR,EOR,TST,B.NE); ELIST 5 (4, + 1 address ADD/SUB off x2).
PAIRM = KB130 + Fraction(1, r)
VISIT = HCHG_VISIT + Fraction(HCHG_END, 2)
# j=1,2 key paths E-priced/group (6.6); defectives global.
PJ = Fraction(2, 3) + 2 * D_WJ + Fraction(1, r - 1)
JROW_W = NODES * N0G * JPATH + NODES * N0G // 3 * SJOIN
check("E-charge per visit 4.5: 4 per visit plus the B after an even visit count, at most one per two visits",
   VISIT == Fraction(9, 2))
SG = [NODES * LEAVES * q * Q2 for q in Q1G]
WG, WGW, SAMEX, WGJ = [], [], [], []
# phase-4 exchanges: per-group row under PHI, EXCH_G (itemized) per node and sub-phase.
EXCH_ROW = NODES * 4 * EXCH_G
for _g, _q in enumerate(Q1G):
 _s = SG[_g]
 _leaf = NODES * _q * (fwd_count(Q2) + INV + INV_CALL + BATCH_SETUP + GEN_MULX) \
  + int(NODES * _q * (BWD_INC * Q2 + BWD_FIX)) + int(NODES * _q * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX))
 # h2s: 6 per record; a full arriving block 8 fixed (unrolled); each other block (at most SIDE + 1 per node: one final
 # per column-mate sender and the empty call) at most 16 fixed and K/8 - 1 loop passes
 _route = _s * (EMIT * K + FLUSH + DONE) // K + _s * H2S_ARR + (_s // K) * H2S_FULL \
  + NODES * (SIDE + 1) * (H2S_ARR_FIX + LOOP8 * (K // 8 - 1)) \
  + _s * BLOCKH // K + NODES * (FINAL * SIDE + FIN_FIX) + NODES * (FINALH * SIDE + FIN_FIX) \
  + NODES * (INIT3_BKT * SIDE + 1 + TINIT_ENTRY * SIDE + 1)
 _gp = NODES * LEAVES * MQ1[_g] * Q2
 _samex = NODES * _q * (6 * Q2) * (6 * Q2 - 1) // 2
 _evis = (_s * _gp + _s * (_s - 1) // 2 - _samex) * PAIRM
 _gfix = _s * HCHG_NEW + NODES * STOP_CHECK
 _group = _gfix + VISIT * _evis
 assert isinstance(_group, (int, Fraction)) and _group >= 0
 _jkey = JPATH * _s * PJ + _s // 3 * SJOIN
 WGJ.append(_jkey)
 WG.append(_leaf + _route + _group + _jkey + EXCH_ROW)
 SAMEX.append(_samex)
 _wgroup = _gfix + _s * HOME_MAX * (HCHG_VISIT + HCHG_END + HCHG_VER + HCHG_LS1) + EXCH_ROW
 assert isinstance(_wgroup, int)
 WGW.append(_leaf + _route + _wgroup)
check("group sizes sum to the grouped q1, each >= 1; group records multiples of k",
   sum(Q1G) == Q1GD and min(Q1G) > 0 and all(s % K == 0 for s in SG))

ZERO_INIT = NODES * (2 ** HB // ZW * ZERO_PASS + 1)
check("hash-head zeroing: zero16384 zeroes 16384 words per pass with 512 (DC ZVA, ADD #256; a 256-byte block), SUBS and "
   "B.NE (1026); 2^58 is a multiple of 16384, so the row is 2^58/16384 passes plus the RET", ZERO_PASS == 1026
   and ZW == 512 * 256 // 8 and 2 ** HB % ZW == 0 and ZERO_INIT == NODES * (2 ** (HB - 14) * 1026 + 1))
HR_INIT = NODES * (INIT3_BKT * SIDE + 1)
E_PRE = NODES * NODE_INIT + NODES * (Q1GD + Q2 + 1) * COIN_SCALAR + NODES * Q2 * GEN_MULY + NODES * HEADER \
 + NODES * (INIT3_BKT * SIDE + 1 + TINIT_ENTRY * SIDE + 1) + ZERO_INIT + HR_INIT + NODES * TBLD
E_WORST_REST = E_PRE + NODES * GEN_MUL + 2 * (SIDE + NODES) * MSG \
 + NODES * DELTA_F * max(FIX_EXTRA.values()) \
 + NODES * CAND_CAP * RESOLVE + NODES * REQ_CAP * REQ_SERVE \
 + NBAR_ALLOW * BARRIER + NODES * PHASE_CTRL_G + OUTPUT + G_GROUPS * 2 * (SIDE + NODES) * MSG
_KB182 = KB(KEYBITS + D1 + D2)
_SIPT = RHO * Fraction(5, 2) * (1 + 11 * (D_WJ + D_LAB))
# v6: a 182-pair sits at one home: _KB182 alone.
ELIST = RHO * (Fraction(SG_FULL) * SG_FULL * _KB182
      + (2 * SG_FULL + PARTNER_PAIRS) * _KB182) + _SIPT * _KB182 * SG_FULL
# K2 hits (6.6): main pairs at KB(136) + 1/r per group under PHI; same-X at their visits; dependent terms global.
KB136 = KB(D1 + D2 + 64)
PAIRK2 = KB136 + Fraction(1, r)
EK2G = [RHO * (_s * _gp + _s * (_s - 1) // 2 - _sx) * PAIRK2 for (_s, _gp, _sx) in
    zip(SG, [NODES * LEAVES * m * Q2 for m in MQ1[:-1]], SAMEX)]
EK2SX = PAIRM * sum(SAMEX)
EK2DEP = RHO * (2 * SG_FULL + PARTNER_PAIRS) * KB64 + _SIPT * KB64 * SG_FULL
EK2 = sum(PHI[_g] * EK2G[_g] for _g in range(G_GROUPS)) + EK2SX + EK2DEP
SXV = VISIT * PAIRM * sum(SAMEX)
BRACKET = RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) + 2 * SG_FULL + PARTNER_PAIRS) \
 + _SIPT * SG_FULL
DEFV = VISIT * D_BATCH * (HOME_MAX * SG_FULL + sum(SAMEX))
# v6: defectives' own walks verify/list too (<= HOME_MAX each).
DEFK2 = HCHG_VER * 2 * D_BATCH * BRACKET + HCHG_VER * D_BATCH * HOME_MAX * SG_FULL
DEFL = HCHG_LS1 * 2 * D_BATCH * BRACKET + HCHG_LS1 * D_BATCH * HOME_MAX * SG_FULL
# DEFJ: JPATH/leaf at unconditional d_batch.
DEFJ = JPATH * D_BATCH * SG_FULL
# candidates and fixups at expected counts (7); their caps stay in R_WORST
E_CAND = RESOLVE * (ELIST + DEFL / HCHG_LS1)
E_FIXUP = NODES * max(FIX_EXTRA.values()) * Q1GD * (E_FIX + D_BATCH * FIX_BATCH)
CAPS_W = NODES * CAND_CAP * RESOLVE + NODES * DELTA_F * max(FIX_EXTRA.values())
E_TOTAL = E_WORST_REST - CAPS_W + E_CAND + E_FIXUP + sum(PHI[_g] * WG[_g] for _g in range(G_GROUPS)) + SXV \
 + EK2 * HCHG_VER + ELIST * HCHG_LS1 + DEFV + DEFK2 + DEFL + DEFJ
R_WORST = E_WORST_REST + JROW_W + sum(WGW)

def claim6_of(x):
 c = Fraction(math.ceil(lg(x) * 10 ** 6), 10 ** 6)
 for _ in range(50):
  if pow2_interval(c)[0] >= x and pow2_interval(c - Fraction(1, 10 ** 6))[1] < x:
   return c
  c += Fraction(1, 10 ** 6) if pow2_interval(c)[0] < x else -Fraction(1, 10 ** 6)
 raise AssertionError("worst claim not certified")

WORST6 = claim6_of(R_WORST)
# 6b: declared literal; checker compares (red line o).
DECLARED_WORST6 = Fraction(188604756, 10 ** 6)
check("disclosed worst case 2^%d/10^6: least 6-decimal bound >= R_WORST" % (WORST6 * 10 ** 6),
   pow2_interval(WORST6)[0] >= R_WORST and pow2_interval(WORST6 - Fraction(1, 10 ** 6))[1] < R_WORST)
check("declared worst case 2^188604756/10^6 equals the least 6-decimal bound (manifest match)",
   WORST6 == DECLARED_WORST6 and pow2_interval(DECLARED_WORST6)[0] >= R_WORST
   and pow2_interval(DECLARED_WORST6 - Fraction(1, 10 ** 6))[1] < R_WORST)
# 6b: E-total display bound.
E6 = claim6_of(E_TOTAL)
check("E-total display 2^%d/10^6: least 6-decimal bound >= E_TOTAL" % (E6 * 10 ** 6),
   pow2_interval(E6)[0] >= E_TOTAL and pow2_interval(E6 - Fraction(1, 10 ** 6))[1] < E_TOTAL)
W_E = E_TOTAL / C

def claim15_of(w):
 # least 10^-15 multiple c with 2^c >= w, both ways; 80-digit start.
 getcontext().prec = 80
 start = (Decimal(w.numerator).ln() - Decimal(w.denominator).ln()) / Decimal(2).ln() * 10 ** 15
 c = Fraction(int(start.to_integral_value(rounding=ROUND_CEILING)), 10 ** 15)
 for _ in range(400):
  if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 15))[1] < w:
   return c
  c += Fraction(1, 10 ** 15) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 15)
 raise AssertionError("E-claim not certified")

CLAIM_E = claim15_of(W_E)
check("E-work below the worst case", E_TOTAL < R_WORST,
   "E 2^%.6f worst declared 2^%.6f" % (float(E6), float(DECLARED_WORST6)))
check("E-claim = least 15-decimal c with 2^c >= E/1000", pow2_interval(CLAIM_E)[0] >= W_E
   and pow2_interval(CLAIM_E - Fraction(1, 10 ** 15))[1] < W_E, "%d/10^15" % (CLAIM_E * 10 ** 15))
JROW_E = sum(PHI[_g] * WGJ[_g] for _g in range(G_GROUPS)) + DEFJ
check("key paths E-priced: PJ = 2/3 + 2 d_wj + 1/(r - 1) per leaf + a B per sign, per group under PHI, plus DEFJ; below every-leaf row",
   PJ > Fraction(2, 3) and D_WJ > Fraction(1, r - 1) and JROW_E < JROW_W
   and all(WGJ[_g] == JPATH * SG[_g] * PJ + SG[_g] // 3 * SJOIN for _g in range(G_GROUPS)) and DEFJ == JPATH * D_BATCH * SG_FULL,
   "2^%.4f, every-leaf row 2^%.4f" % (lg(JROW_E), lg(JROW_W)))
check("worst case charges every leaf's key path and every sign's B: R_WORST = E_WORST_REST + JROW_W + sum WGW",
   R_WORST - JROW_W == E_WORST_REST + sum(WGW) and JROW_W == NODES * N0G * JPATH + NODES * N0G // 3 * SJOIN)
check("K2 verifications: main pairs (walks' pair counts, same-X excluded) at joint KB(136) + 1/r per group under PHI; "
   "same-X pairs global at most their visits (PAIRM); dependent terms global as before; below the old global product row",
   KB136 == KB(D1 + D2 + 64) and KB136 < KB(D1 + D2) * KB64 and PAIRK2 == KB136 + Fraction(1, r)
   and all(EK2G[_g] == RHO * (SG[_g] * (NODES * LEAVES * MQ1[_g] * Q2) + SG[_g] * (SG[_g] - 1) // 2 - SAMEX[_g])
      * PAIRK2 for _g in range(G_GROUPS)) and EK2SX == PAIRM * sum(SAMEX)
   and EK2DEP == RHO * (2 * SG_FULL + PARTNER_PAIRS) * KB64 + _SIPT * KB64 * SG_FULL
   and EK2 == sum(PHI[_g] * EK2G[_g] for _g in range(G_GROUPS)) + EK2SX + EK2DEP
   and EK2 < RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) * KB64 + (2 * SG_FULL + PARTNER_PAIRS) * KB64)
      + _SIPT * KB64 * SG_FULL,
   "2^%.4f" % lg(EK2 * HCHG_VER))
check("barrier allowance covers the grouped run", NBAR_ALLOW >= NBAR_COUNT, "%d of %d" % (NBAR_COUNT, NBAR_ALLOW))
PER_LEAF_E = E_TOTAL / SG_FULL

SLOTS = SIDE * CAP_P * 32
FIXED = [("code, constants, stack", 2 ** 20), ("barrier words: 5 granules of at most 2 KiB", 5 * 2048),
   ("header H and the rows of beta below it", 1024),
   ("X table (one group): H's rows and values, -T(a), -T(beta a), y^2 per X", max(Q1G) * 512),
   ("X affine store (all groups): (x, y) per X", Q1GD * 64),
   ("scalars a, b, c", (Q1GD + Q2 + 1) * 32), ("batch entries (Y words, corrections, D, acc, fractions)", Q2 * 512),
   ("the ports' message blocks", 2 * 64), ("LC tables at 0 and 2^39", 2 * SIDE * 8),
   ("RC tables at 2^43 and 2^43 + 2^39", 2 * SIDE * 8), ("record blocks at 2^50 and 2^51",
    2 * SIDE * 2 ** 14),
   ("final-cursor tables E1 and E2 (two, by group parity)", 3 * SIDE * 32),
   ("hash-read cursors at 2^48", SIDE * 8), ("join state at 2^48 - 160 (list cursor and end, C, -2C), h2s's p word at 2^48 - 168", 40),
   ("candidate list", CAND_CAP * 16),
   ("request and reply buffers", REQ_CAP * 64 + 8 * 2 * CAND_CAP * 64), ("relay buffer", RELAY_CAP * 64),
   ("phase-4 count, offset and cursor tables", 3 * SIDE * 32)]
FIXED_BYTES = sum(b for _, b in FIXED)

SLOTS_G = 2 ** D1 * CAP_PG * 32
PHASE3G = [("per-group hop-1 slots", SLOTS_G),
     ("hash heads: 2^58 zeroed words", 2 ** HB * 8),
     ("hash next-array: one 16-byte entry (next address, K2) per R2 slot position", SIDE * CAP_P * 16)]
PHASE3G_BYTES = sum(b for _, b in PHASE3G)
check("grouped R1 partition (slots + hash heads + hash next-array) fits in R1's range; the next-array "
   "covers every receive pattern up to the slot caps by construction", PHASE3G_BYTES <= SLOTS,
   "2^%.4f" % lg(PHASE3G_BYTES))
R1_BASE = 2 ** 52
R2_PAD = 32  # R2 starts 32 bytes past R1's end, so every pair slot base and cursor c is 32 mod 64 (LS64 pairs, section 4)
R2_BASE = R1_BASE + SLOTS + R2_PAD
check("LS64 alignment: window slot n at 2^14-aligned base + 32n and R2 slot c + 32(n - 1) with c = 32 mod 64 (R2 base 32 mod "
   "64, pair slots 32*cap_p apart, cursors moving by 32k or 32L only within a block): both 64-byte aligned exactly for even n",
   R2_BASE % 64 == 32 and (32 * CAP_P) % 64 == 0 and (32 * K) % 64 == 0 and BLKC % 64 == 0 and PB2C % 64 == 0
   and all(((BLKC + 32 * n) % 64 == 0) == ((R2_BASE + 32 * (n - 1)) % 64 == 0) == (n % 2 == 0) for n in range(1, K + 1)))
HEAD_BASE, NEXT_BASE = 2 ** 61, 2 ** 62
NEXT_C = NEXT_BASE - (R2_BASE >> 1)
check("R1 layout: hop-1 slots at 2^52, the heads at 2^61 (2^61 bytes: the head address is 2^61 with the bucket in bits "
   "3..60, one BFI), the next-array at 2^62 (entry (A>>1) + C per R2 address A, C = 2^62 - R2/2) up to R1's end; "
   "disjoint, in R1", R2_BASE % 2 == 0 and R1_BASE + SLOTS_G <= HEAD_BASE and HEAD_BASE + 2 ** HB * 8 <= NEXT_BASE
   and NEXT_BASE + (SLOTS >> 1) <= R2_BASE - R2_PAD and HEAD_BASE % 2 ** 61 == 0 and 2 ** HB * 8 == 2 ** 61
   and (SLOTS >> 1) == PHASE3G[2][1], "heads end 2^%.4f, next-array end 2^%.4f, R1 end 2^%.4f"
   % (lg(HEAD_BASE + 2 ** HB * 8), lg(NEXT_BASE + (SLOTS >> 1)), lg(R2_BASE - R2_PAD)))
for name, b in FIXED + [("R1 at 2^52: hop-1 pair slots (phases 1-2), then the partition below", SLOTS)] + PHASE3G \
  + [("R2 at 2^52 + |R1| + 32: hop-2 pair slots (phases 2-4)", SLOTS)]:
 print("  %-58s 2^%.4f bytes" % (name, lg(b)))
check("fixed regions (one group's X table, all groups' X affine store) lie below 2^52", FIXED_BYTES <= 2 ** 52,
   "2^%.4f" % lg(FIXED_BYTES))
ALIGNED = 4 * SIDE * 8 + 2 * SIDE * 2 ** 14 + SIDE * 8
check("aligned tables: LC at 0 and 2^39, RC at 2^43, hash-read cursors at 2^48, column blocks at 2^50, row blocks at "
   "2^51; the rest fits between 2^44 and 2^48",
   SIDE * 8 == 2 ** 39 and SIDE * 2 ** 14 == 2 ** 50 and PB1C == 1 << 50 and PB2C == 1 << 51
   and PB2C + (1 << 50) <= 1 << 52 and FIXED_BYTES - ALIGNED <= (1 << 48) - (1 << 44),
   "2^%.4f" % lg(FIXED_BYTES - ALIGNED))
check("record block: header and k slots within 2^14; bit 13 set exactly at slot k; bases/stride keep low 13 bits zero; "
   "cap_p and the group records multiples of k (keep lemma)",
   32 + 32 * K <= 2 ** 14 and all((((32 * n) >> 13) & 1 == 1) == (n == K) for n in range(1, K + 1))
   and PB1C & 0x1FFF == 0 and PB2C & 0x1FFF == 0 and BLKC & 0x1FFF == 0 and CAP_P % K == 0
   and SG_FULL % K == 0 and all(s % K == 0 for s in SG))
check("bwd's and keys' loop end: the q2 512-byte batch entries fill the lower 2^34 bytes of their 2^35-aligned block",
   Q2 * 512 == 2 ** 34)
check("record fields: k < 2^28 in w2 bits 36-63; the sign's code in w3 bits 36-38 is 2 (+, emit at entry + 128) or 4 "
   "(-, entry - 256), complemented on j = 1 (5, 3), the four distinct; l < 2^25 (l - 1 mod 2^25 for s = -1) in bits "
   "39-63, complemented on j = 1", Q1GD <= 2 ** 28 and Q2 <= 2 ** 25 and (128 >> 6, (2 ** 9 - 256) % 2 ** 9 >> 6) == (2, 4)
   and len({2, 4, 7 - 2, 7 - 4}) == 4)
check("R1 and R2 end below 2^64, below the beta tables at 2^64 - 2^51 (2^43 entries of 256 bytes, slot k at +32k)",
   R2_BASE + SLOTS <= BT_BASE and BT_BASE + (2 ** BT_W << 8) == 2 ** 64 and 6 * BT_W >= 256 > 5 * BT_W,
   "2^%.6f" % lg(R2_BASE + SLOTS))
CURSOR_END = 2 ** 52 + 32 * CAP_P * (SIDE - 1) + 32 * N0G
check("cursors cannot wrap: the last slot's base (hop-2 stride, at least hop 1's) plus 32 per leaf of a node over all "
   "groups stays below 2^64", CAP_PG <= CAP_P and CURSOR_END < 2 ** 64, "2^%.4f" % lg(CURSOR_END))
check("hash-read cursors occupy [2^48, 2^48 + 2^39), below the column blocks at 2^50", (1 << 48) + SIDE * 8 <= PB1C)
PEAK = NODES * (FIXED_BYTES + 2 * SLOTS + R2_PAD + (2 ** BT_W << 8))
check("peak memory 2^72 * (fixed + R1 + 32 + R2 + beta tables) = 2^%.4f bytes" % lg(PEAK), PEAK < 2 ** 136)
check("caps exceed crude means: cap_p > m_p2 + 2^20 at hop 2, cap_pg > m_p1 + 2^17 per group at hop 1",
   CAP_P > TFULL["m_p2"] + 2 ** 20 and CAP_PG > TFULL["m_p1"] + 2 ** 17,
   "m_p2 2^%.3f, m_p1 2^%.3f" % (lg(TFULL["m_p2"]), lg(TFULL["m_p1"])))
check("preprocessing E_PRE (set-up, coins, Y generators, beta rows, LC/RC tables, hash-head zeroing, read-cursor "
   "init) = 2^%.7f units, included"
   % lg(Fraction(E_PRE, C)), E_PRE < E_TOTAL < R_WORST)

DECLARED = dict(time_log2=CLAIM_E, worst_log2=DECLARED_WORST6, memory_log2_bytes=None, preprocessing_log2=None,
    success=Fraction(39, 100),
    per_leaf=Fraction(cdiv(PER_LEAF_E.numerator * 10 ** 4, PER_LEAF_E.denominator), 10 ** 4))

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

DECLARED["memory_log2_bytes"] = least_grid_up(PEAK, 10 ** 4)
DECLARED["preprocessing_log2"] = least_grid_up(Fraction(E_PRE, C), 10 ** 4)
check("declared E time_log2 %d/10^15: least 15-decimal value with 2^c >= E/1000" % (DECLARED["time_log2"] * 10 ** 15),
   DECLARED["time_log2"] == CLAIM_E and least_at_or_above(DECLARED["time_log2"], W_E, 10 ** 15))
check("declared memory 2^%s bytes: least 4-decimal bound >= peak (2^%.7f)"
   % (format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK)),
   least_at_or_above(DECLARED["memory_log2_bytes"], PEAK, 10 ** 4))
check("declared preprocessing 2^%s units: least 4-decimal bound >= E_PRE / 1000"
   % format(float(DECLARED["preprocessing_log2"]), ".4f"),
   least_at_or_above(DECLARED["preprocessing_log2"], Fraction(E_PRE, C), 10 ** 4))
check("declared success 0.39 <= the bound", DECLARED["success"] <= TFULL["lo"])
# the declared figures as literals, and as written in submission.json's claim statement (not only recomputed)
STATED = dict(time_log2=Fraction(124534495156618070, 10 ** 15), worst_log2=Fraction(188604756, 10 ** 6),
     memory_log2_bytes=Fraction(1359990, 10 ** 4), preprocessing_log2=Fraction(1161108, 10 ** 4),
     success=Fraction(39, 100), per_leaf=Fraction(2230163, 10 ** 4))
import json, os, re
_CL = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "submission.json")))["claims"][0]
_ST, _SC = _CL["statement"], _CL["score"]
check("stated figures (literals here, submission.json's statement and score) equal the certified ones",
   STATED == DECLARED and Fraction(int(_SC["numerator"]), int(_SC["denominator"])) == STATED["time_log2"] and all(re.search(re.escape(t) + r"(?![0-9])", _ST) for t in ("2^124.534495156618070", "2^188.604756", "2^135.9990", "2^116.1108",
              "223.0163", "0.39")))
check("declared E work per leaf %s: least 4-decimal bound >= E/S (%.7f)"
   % (format(float(DECLARED["per_leaf"]), ".4f"), float(PER_LEAF_E)),
   DECLARED["per_leaf"] >= PER_LEAF_E > DECLARED["per_leaf"] - Fraction(1, 10 ** 4))

MULX, SQRX, ADDX = 105 + 4, 80 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX

GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_TEST = 256 * 2 * (9 + 4) + 64

GEN_CTRL = 1600
NEGT = 30

RED = 9
Y_ENTRY = 3 * SQRX + 2 * MULX + 2 * ADDX + 2 * NEGT + 7 * RED + 14

X_ENTRY = MADD_ + 3 * SQRX + 11 * MULX + 2 * ADDX + 2 * NEGT + 17 * RED + 32

check("generator (exact arithmetic): doubling 903, mixed addition 1,391; with the exceptional-case tests, control, "
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
check("grouped phase control: 12 kernel calls per group at most 33 each, stop test and loop 8, h2s's p word 4 (CSEL, CSEL, ADD, STR), "
   "over %d groups < 2^17" % G_GROUPS, G_GROUPS * (12 * 33 + 8 + 4) < PHASE_CTRL_G)
check("barrier: at most 131 instructions per node, in-flight counts and deliveries included (2^10 per node is > 7x)",
   BARRIER_ITEM == 131 and BARRIER > 7 * NODES * BARRIER_ITEM)

EXCH_ITEM = 12 * SIDE + 4 * (14 + MSG_LINK) * (8 * REQ_CAP)
check("phase-4 count exchange charged at its itemization per node and sub-phase: 12 per entry of the 2^36-entry table, "
   "88 per message for 8*2^29 messages", EXCH_ITEM == EXCH_G == 12 * 2 ** 36 + 88 * 2 ** 32
   and EXCH_ROW == NODES * 4 * EXCH_G, "2^%.4f" % lg(EXCH_ITEM))

check("a broadcast message costs at most 16, in-flight count and delivery included, within its price of 40",
   2 + 3 + 3 + MSG_LINK <= MSG)

CAND_LEAVES = INV + 5 * SQRX + 18 * MULX + 20 * ADDX
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 9 * 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("candidate: messages, both records' three candidate leaves, nine unit tests, solve mod r, [x]P = Q = %d; "
   "charged above twice that" % CAND_ITEM, RESOLVE == 2 * CAND_ITEM + 1,
   "%.3fx" % (math.floor(10 ** 3 * RESOLVE / CAND_ITEM) / 10 ** 3))
PROGRAM_BYTES = 2 ** 18
check("program delivery: 2^18 bytes (29,132 kernel instructions expanded and the rest) = 2^13 messages over 2 hops at "
   "40 = 2^19.3 < node set-up 2^20", PROGRAM_BYTES // 32 * 2 * MSG < NODE_INIT)
# analytic (itemized, not kernel-counted) E-items, batch items over every group (an upper bound on their E-rows)
ANALYTIC_E = NODES * (Q1GD + Q2 + 1) * COIN_SCALAR + NODES * Q2 * GEN_MULY + NODES * HEADER + NODES * NODE_INIT \
 + NODES * Q1GD * GEN_MULX + NODES * GEN_MUL \
 + NODES * Q1GD * (BATCH_SETUP + INV_CALL) + (1 + G_GROUPS) * 2 * (SIDE + NODES) * MSG \
 + E_CAND + NODES * REQ_CAP * REQ_SERVE + sum(PHI[_g] * EXCH_ROW for _g in range(G_GROUPS)) \
 + NBAR_ALLOW * BARRIER + NODES * PHASE_CTRL_G + OUTPUT
check("analytic allowances (the E-averaged exchange rows included) are below 2^-12.4 of E-work",
   ANALYTIC_E * pow2_interval(Fraction(62, 5))[1] < E_TOTAL, "2^%.3f" % lg(Fraction(ANALYTIC_E, E_TOTAL)))

print()
print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  success >= %.10f" % (Q1GD, N0G, lg(SG_FULL), float(TFULL["mu"]),
                 math.floor(TFULL["lo"] * 10 ** 10) / 10 ** 10))
print("E_total 2^%.6f instructions  W_E 2^%.10f  E-claim %d/10^15  per leaf %.7f declared %s  R_worst 2^%.6f"
   % (float(E6), lg(W_E), CLAIM_E * 10 ** 15, float(PER_LEAF_E), format(float(DECLARED["per_leaf"]), ".4f"),
   float(DECLARED_WORST6)))
print("key paths: E-priced 2^%.6f (every-leaf row 2^%.6f in the worst case, JPATH %d); Y items at %d (was %d); join "
   "%d fixed, %s per visit" % (lg(JROW_E), lg(JROW_W), JPATH, GEN_MULY, GEN_MUL, HCHG_NEW, VISIT))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
 format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK),
 format(float(DECLARED["preprocessing_log2"]), ".4f"), lg(Fraction(E_PRE, C))))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
 sys.exit(1 if FAILED else 0)
