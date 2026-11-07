import functools, json, math, os, re, sys
from decimal import Decimal, ROUND_CEILING, getcontext
from fractions import Fraction
FAILED = []
def check(name, ok, detail=""):
 print("%-66s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
 if not ok: FAILED.append(name)
def lg(x):
 f = Fraction(x)
 return math.log2(f.numerator) - math.log2(f.denominator)
def cdiv(a, b): return -((-a) // b)
def exp_interval(x, terms=60):
 assert abs(x) <= 1
 n, d, S = x.numerator, x.denominator, 1 << 400
 lo = hi = 0; tl = th = S
 for k in range(terms + 1):
  if k: a, b = sorted((tl * n, th * n)); tl, th = a // (d * k), cdiv(b, d * k)
  lo, hi = lo + tl, hi + th
 rem = cdiv(max(-tl, th) * abs(n) * 3, d * (terms + 1))
 return Fraction(lo - rem, S), Fraction(hi + rem, S)
def ln2_interval(bits=300):
 z, z2 = Fraction(1, 3), Fraction(1, 9)
 tot, pw, k = Fraction(0), z, 0
 while True:
  tot += pw / (2 * k + 1)
  rest = pw * z2 / ((2 * k + 3) * (1 - z2))
  if rest < Fraction(1, 2 ** bits): return 2 * tot, 2 * (tot + rest)
  pw *= z2; k += 1
LN2 = ln2_interval()
@functools.lru_cache(maxsize=None)
def pow2_interval(c):
 n = c.numerator // c.denominator
 f = c - n
 lo, hi = exp_interval(f * LN2[0])[0], exp_interval(f * LN2[1])[1]
 sc = Fraction(2) ** n
 return lo * sc, hi * sc
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
C0 = 2 ** 32 + 977
BETA_P = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
C = 1000
NODES = 2 ** 68
M_G = 1024
N_CLUST = NODES // M_G
F1, F2 = 2 ** 42, 2 ** 42
N_REL, N_HOME = F1, F1 * F2
ALLN = NODES + N_REL + N_HOME
SIDE = F1
MX = 16
YPC = 33
YROUNDS = 63549
L_CHAIN = YROUNDS * YPC
Q2 = MX * L_CHAIN
D1, D2, D3, D4 = 42, 42, 64, 22
KEYBITS = 64 + D4
PART = 64 + D1
KBF = 3195
K1, K2B = 16, 16
THRESH = Fraction(3901, 10 ** 4)
LEAVES = 6
FWD_INC, FWD_FIX = Fraction(1081, 16), 7
def fwd_count(q2):
 assert q2 % MX == 0
 return FWD_FIX + (1081 * q2) // 16
BWD_INC, BWD_FIX, BWD_YSW = Fraction(16233, 16), 7, 0
def bwd_count(q2):
 assert q2 % MX == 0
 return BWD_FIX + (16233 * q2) // 16
FIX_SITES_QUAD = dict(SHX=54, FDX=12, F5=8, FS=6, FADDW=8, BLK=28)
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4, BLK=5)
EMIT = 5
EMIT_INC = LEAVES * EMIT
KEY_INC, KEY_FIX = Fraction(513, 4), 5
KEY_W3 = Fraction(2, 1)
K = K1
FLUSH = 10
H2S_BLOCK_DROP = 9
FLUSH2_X, FLUSH2_D = 10, 8
LOOP8 = 0
JPATH, SJOIN = 1, 1
INV = 22549
TINIT_ENTRY = 3
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 15, 4, 2, 2
FINH_DEST, FINH_SENT = 16, 5
H2S_ARR, H2S_ONE = 6, 1
H2S_FULL = -1
FLUSH1_X, FIN1_X = 8, 14
H2S_EMPTY = 6
H2S_ARR_FIX = 14 + H2S_ONE
def h2s_count(n):
 if n == K: return H2S_FULL + H2S_ARR * K
 if n == 0: return H2S_EMPTY
 return 14 + H2S_ONE * (n % 2) + H2S_ARR * n
check("h2s: full 16-slot block 95 (6 persistent port regs, 64B-aligned LINK+64); partial 14 + (n%2) + 6*n <= 15 + 6*n",
   h2s_count(K) == 95 and max(h2s_count(n) - H2S_ARR * n for n in range(K)) == H2S_ARR_FIX and H2S_EMPTY <= H2S_ARR_FIX)
INIT3_BKT = 4
check("key pass: 128.25 = 8 loads/steps + 52 sums + 2.0 w3 + 6*(6 index + emit 5) + 1/4 TBZ (2x KEYPAIR7 = 4 inc/loop)",
   KEY_INC == 8 + 52 + KEY_W3 + 6 * (6 + EMIT) + Fraction(1, 4) and KEY_FIX == 5)
check("bwd: EFG_QUAD_4 + 8 * BWDP16 + Y-switch = 16233 per 16 incidences (1014.5625/inc), 7 fixed; Q2 % 16 == 0",
   BWD_INC == Fraction(16233, 16) and Q2 % MX == 0 and BWD_FIX == 7 and bwd_count(Q2) == BWD_FIX + (16233 * Q2) // 16)
check("fwd: 2 hoisted loads + 4*(LKP7_PAIR 97 + 2 store + 2 P1X 84) + 1 add + 4 loop = 1081 per 16 incidences (67.5625/inc)",
   FWD_INC == Fraction(1081, 16) and Q2 % MX == 0 and FWD_FIX == 7 and fwd_count(Q2) == FWD_FIX + (1081 * Q2) // 16)
check("j = 1 and j = 2 paths (exact-carry ADDS, ADC; CMP, B.HS; CMN, B.HS; MVN of sigma; emit 5) exceed j = 0 by 1; SJOIN 1",
   (4 + 2 + 1 + EMIT) - (6 + EMIT) == JPATH and SJOIN == 1)
check("emit: AND or BIC (8*d1), LDADD, STP K2,K3, STP x2,x26, TBZ of bit 9 for slot k1 = 16",
   EMIT == 1 + 1 + 2 + 1)
DONE = 2
RCV, WAKE, EXPECT = 1, 8, 8
DONE_C = DONE + RCV + WAKE
RCV1, FIN_CNT, H2S_DONE = 4, 6, 4
EXPECTH = 16
RDONE = 1 + H2S_DONE
DONE_C1 = DONE + RCV1
DLV = 6
MSG_LINK = DONE_C + DLV
FINAL = FIN_DEST + FIN_SENT + DONE_C1 + H2S_ARR_FIX + LOOP8 * (K // 8 - 1) + FIN_CNT
check("hop-1 executed code within its charge: flush %d <= FLUSH + DONE, fin1 %d <= FIN_DEST + FIN_SENT + DONE" % (FLUSH1_X, FIN1_X),
   FLUSH1_X <= FLUSH + DONE and FIN1_X <= FIN_DEST + FIN_SENT + DONE and FLUSH2_X <= FLUSH and FLUSH2_D <= H2S_BLOCK_DROP)
GEN_MUL = 645410
GEN_MULX = GEN_MUL + 4
YR_FIX, YR_BLK, YB = 22556, 43258, 256
CHAIN_ADD = 2 ** 15
CHAIN_SET = 2 * GEN_MUL + 17 * CHAIN_ADD
YSND, YDLV, YDLVE = 24, 32, 23
YRING_TOT = YSND + (M_G - 2) * YDLV + YDLVE
def y_group(q):
 q_gen = cdiv(q, M_G)
 gen_y = YROUNDS * (cdiv(q_gen, YB) * YR_FIX + q_gen * YR_BLK) + q_gen * CHAIN_SET
 ring_y = cdiv(q * L_CHAIN * YRING_TOT, M_G)
 return gen_y + ring_y
XT_SLOTS = 7
XT_W, XT_WHI = 37, 36
XT_ENTRIES_SINGLE = 4 * 2 ** XT_W + 3 * 2 ** XT_WHI
XTBLD_PAIR44 = XT_SLOTS * (4 + 2 * 11) + 20 * (2 * XT_ENTRIES_SINGLE)
XTBLD_ALL = 40 * XTBLD_PAIR44
XSETUP = 2 ** 17
ENTRY_X = 16
COIN_SCALAR = 64
BATCH_SETUP = 1280
INV_CALL = 16
NODE_INIT = 2 ** 20
HEADER = 72
BT_BASE, BT_W = 0xfff8 << 48, 43
TBLD = 4 + 7 * 11 + 20 * (6 * 2 ** BT_W + 2 ** (256 - 5 * BT_W))
MSG = 40
CAND_CAP = 2 ** 27
RESOLVE = 1491011
REQ_CAP = 49 * 2 ** 33
REQ_SERVE = 2 ** 10 + GEN_MUL
RELAY_CAP = 2 ** 40
BARRIER = 2 ** 11 * ALLN
OUTPUT = 2 ** 20
DELTA_F = 2 ** 52
QQ = Q2 // 4
FIX_BATCH = sum(FIX_SITES_QUAD.values()) * QQ
P1 = Fraction(6 * (r - 1) + 1, r * r)
TWO_BOUND = Fraction(2, r - 1)
SHX_BAD = 2 ** 64 * 2 ** 33
FOLD_BAD = 2 ** 66 * C0
SQ_BAD = C0 * C0 + 2 * C0
FS_P1_BAD = C0
QP = Q2 // 2
SITES = [
 ("fwd: row shifts of acc_(j-1) (uniform on F_p* via rho)", 3 * QP, SHX_BAD, 1),
 ("fwd: fold of acc_j", QP, FOLD_BAD, 1),
 ("bwd: row shifts of R_j", 3 * QP, SHX_BAD, 1),
 ("bwd: fold of inv_e = R_j acc_(j-1)", QP, FOLD_BAD, 12),
 ("bwd: fold of R_(j-1) = R_j e_j", QP, FOLD_BAD, 1),
 ("bwd: row shifts of D_2 = B - A_2", 3 * QP, SHX_BAD, 6),
 ("bwd: fold of 1/D_1 = D_2 inv_e", QP, FOLD_BAD, 6),
 ("bwd: square fold of inv^2 = 1/D^2 (P1)", Q2, SQ_BAD, 12),
 ("bwd: fold of N = P y_Y (P1, deg-4 in y_Y)", Q2, FOLD_BAD, 12),
 ("bwd: FS of u^2 = M - N and FADDW of v^2 = M + N (P1)", 2 * Q2, FS_P1_BAD, 6),
 ("bwd: FADDW of 1/D_2 (C1BSPLIT_BDAY) and E_odd (EDELTA_ONE)", 4 * QQ, 2 ** 98, 8),
 ("bwd: row shifts of E, F, G (lookup sums fold exactly)", 9 * Q2, SHX_BAD, 8),
]
SITES_QUAD = dict(SHX=54, FDX=12, F5=4, FS_FA=8, BDAY_FA=4)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)
FIX_ROW = C0 * 2 ** 192
FIX_FOLD = (2 ** 66 + 1) * 2 ** 128 * C0
FIX_SQ = 2 ** 193
FIX_LK7 = 21 * C0 * 2 ** 128
FIX_LK8 = 28 * C0 * 2 ** 128
FIX_BLK = 30 * C0 * 2 ** 128
FIX_SUB = C0 * 2 ** 192
FIXUPS = [
 ("fwd: lookup-sum fold of e_j (8 terms)", "BLK", QP, FIX_LK8, 12),
 ("fwd: row shifts of acc_(j-1)", "SHX", 3 * QP, FIX_ROW, 1),
 ("fwd: fold of acc_j", "FDX", QP, FIX_FOLD, 1),
 ("bwd: row shifts of R_j", "SHX", 3 * QP, FIX_ROW, 1),
 ("bwd: fold of inv_e", "FDX", QP, FIX_FOLD, 12),
 ("bwd: fold of R_(j-1)", "FDX", QP, FIX_FOLD, 1),
 ("bwd: FS of D_2 = B - A_2", "FS", QP, FIX_SUB, 6),
 ("bwd: row shifts of D_2", "SHX", 3 * QP, FIX_ROW, 6),
 ("bwd: fold of 1/D_1", "FDX", QP, FIX_FOLD, 6),
 ("bwd: F5 of delta_A * inv_e (C1BSPLIT_BDAY)", "F5", QP, FIX_SUB, 6),
 ("bwd: FADDW of 1/D_2 = 1/D_1 + delta_A * inv_e", "FADDW", QP, FIX_SUB, 6),
 ("bwd: square fold of inv^2 (P1)", "F5", Q2, FIX_SQ, 12),
 ("bwd: lookup-sum fold of M = inv + 2A_X inv^2 (8 terms)", "BLK", Q2, FIX_LK8, 12),
 ("bwd: lookup-sum fold of P = 2y_X inv^2 (7 terms)", "BLK", Q2, FIX_LK7, 12),
 ("bwd: fold of N = P y_Y (P1)", "FDX", Q2, FIX_FOLD, 12),
 ("bwd: FS of u^2 = M - N", "FS", Q2, FIX_SUB, 6),
 ("bwd: FADDW of v^2 = M + N", "FADDW", Q2, FIX_SUB, 6),
 ("bwd: window-sum folds of beta u^2 and beta v^2", "BLK", 2 * Q2, FIX_BLK, 6),
 ("bwd: lookup-sum fold of E_even (8 terms, 2 per quad)", "BLK", 2 * QQ, FIX_LK8, 8),
 ("bwd: F5 of 2*delta_A * b (EDELTA_ONE)", "F5", 2 * QQ, FIX_SUB, 8),
 ("bwd: FADDW of E_odd = E_even + 2*delta_A * b", "FADDW", 2 * QQ, FIX_SUB, 8),
 ("bwd: lookup-sum folds of F, G (7 terms)", "BLK", 2 * Q2, FIX_LK7, 8),
 ("bwd: row shifts of E, F, G", "SHX", 9 * Q2, FIX_ROW, 8),
]
E_FIX = sum(Fraction(n * bad * mult, r - 1) for _, _, n, bad, mult in FIXUPS)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (KBF + 1), (r - 1) * KBF)
RHO = Fraction(r, r - 1) ** 2
CLUST_BATCH = M_G * LEAVES * Q2
FIXED_LEAVES = 2 * CLUST_BATCH
PARTNER_PAIRS = CLUST_BATCH * (CLUST_BATCH - 1)
D_YGEN = Fraction(YROUNDS * (YPC + 2 * 16 + 2), r - 1)
D_BATCH = D_SITES + Q2 * Fraction(6, r - 1) + D_YGEN
D_WJ = TWO_BOUND * 23 * 2 ** 192
D_LAB = TWO_BOUND * 1305 * 2 ** (128 + PART)
FA, FB = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16), int("e4437ed6010e88286f547fa90abfe4c3", 16)
BOX = (p >> PART) + 1
WIDE = 654 * (p >> 128) + 4
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p; box factor (106-bit) < 1 + 1/3195",
   FA * FA - FA * FB + FB * FB == p and (FB * BETA_P - (FA - FB)) % p == 0 and ((FA - FB) * BETA_P + FA) % p == 0
   and (BOX + WIDE + FA) * (BOX + WIDE + 2 * FA - FB) * 2 ** (2 * PART) < p * p * (1 + Fraction(1, KBF)))
def isum(*rs): return sum(x[0] for x in rs), sum(x[1] for x in rs)
TFP_TERMS = ((0, 1), (0, 7), (0, 56))
RAW = isum(*TFP_TERMS); COR = isum(RAW, (-1, 0))
U3, UC3 = (RAW, RAW, (-2 - 2 * RAW[1], -2 * RAW[0])), (COR, COR, (-2 - 2 * COR[1], -2 * COR[0]))
check("3.3's fraction ranges: RAW [0, 64], COR [-1, 64], U_2 [-130, 0], corrected [-130, 2]",
   RAW == (0, 64) and COR == (-1, 64) and U3[2] == (-130, 0) and UC3[2] == (-130, 2))
WR = [[isum(UC3[m], UC3[(m + t) % 3], U3[(m + 2 * t) % 3], (0, 2)) for t in range(3)] for m in (0, 1)]
NR = [(-WR[0][t][1] - WR[1][t][1] - 1, -WR[0][t][0] - WR[1][t][0] + 1) for t in range(3)]
LABEL_ERR = [(min(x[0] for x in rs), max(x[1] for x in rs)) for rs in (WR[0], WR[1], NR)]
check("3.3's error table: W0/W1/N' [-132, 194], [-132, 194], [-389, 264]; D' 1305 residues, widening 654",
   LABEL_ERR == [(-132, 194), (-132, 194), (-389, 264)] and sum(hi - lo for lo, hi in LABEL_ERR) == 1305
   and max(hi - lo for lo, hi in LABEL_ERR) + 1 == 654 and WIDE == 654 * (p >> 128) + 4)
check("wrong class index rate per leaf < 2^-58.47", D_WJ < pow2_interval(Fraction(-5847, 100))[0], "2^%.4f" % lg(D_WJ))
check("site table per X-quad: 54 SHX, 12 FDX, 4 SQ-F5, 8 P1 FS/FA, 4 BDAY FA (lookup/beta sums exact)",
   sum(n for _, n, b, _ in SITES if b == SHX_BAD) == SITES_QUAD["SHX"] * QQ
   and sum(n for _, n, b, _ in SITES if b == FOLD_BAD) == SITES_QUAD["FDX"] * QQ
   and sum(n for _, n, b, _ in SITES if b == SQ_BAD) == SITES_QUAD["F5"] * QQ
   and sum(n for _, n, b, _ in SITES if b == FS_P1_BAD) == SITES_QUAD["FS_FA"] * QQ
   and sum(n for _, n, b, _ in SITES if b == 2 ** 98) == SITES_QUAD["BDAY_FA"] * QQ and len(SITES) == 12)
FIX_RUNS = {m: sum(n for _, ms, n, _, _ in FIXUPS if ms == m) for m in FIX_SITES_QUAD}
check("fixup table per X-quad: 54 SHX, 12 FDX, 8 F5, 6 FS, 8 FADDW, 28 BLK (116 per X-quad = 29/inc)",
   FIX_RUNS == {m: FIX_SITES_QUAD[m] * QQ for m in FIX_SITES_QUAD} and sum(FIX_SITES_QUAD.values()) == 116)
check("a call's fixup word (C0 per fixup, at most 29 q2 fixups) stays below 2^63", FIX_BATCH * C0 < 2 ** 63,
   "2^%.3f" % lg(FIX_BATCH * C0))
check("label-error rate per leaf (106-bit parts) < 2^-10.65", D_LAB < pow2_interval(Fraction(-1065, 100))[0], "2^%.3f" % lg(D_LAB))
check("batch defect rate (site table, D = 0, degenerate chain) < 2^-126.05", D_BATCH < pow2_interval(Fraction(-12605, 100))[0],
   "2^%.3f" % lg(D_BATCH))
G_GROUPS = 1536
G_EXP = 28
def group_sizes(q1):
 ms = [0]
 g_pow = G_GROUPS ** G_EXP
 q32 = q1 ** 32
 for g in range(1, G_GROUPS):
  rhs = q32 * (g ** G_EXP)
  m = int(q1 * ((g / G_GROUPS) ** (G_EXP / 32.0))) + 2
  while (m ** 32) * g_pow > rhs: m -= 1
  while ((m + 1) ** 32) * g_pow <= rhs: m += 1
  ms.append(m)
 ms.append(q1)
 qs = [ms[_g + 1] - ms[_g] for _g in range(G_GROUPS)]
 return qs, ms
HB = 47
HCHG_NEW, HCHG_VISIT, HCHG_END, HCHG_VER = 4, 4, 0, 5
HCHG_LS0, HCHG_LSA = 2, 1
HCHG_LS1 = 5 + HCHG_LSA
check("join per record 4, visit 4, out-of-line hit blocks (.pushsection) HCHG_END 0, d4 check 5, listing 5 + 1 (incl str ST)",
   HCHG_NEW == 4 and HCHG_VISIT == 4 and HCHG_END == 0 and HCHG_VER == 5 and HCHG_LS1 == 6)
DLVC, DLVS1 = 2, 2
DLVHB = 31
DLV_RESTORE = 5
DLVH_SEG, DLVH_SLOW_FIX = 7, 17 + DLV_RESTORE
def dlvh_seg(n, m):
 if m == 0: return 1
 head = n % 2; rest = m - head
 return DLVH_SEG + head + DLVC * m + DLVS1 * (head + rest % 2)
def dlvh_copy(m1, m2, K_=K2B):
 if m1 == 0 and m2 == K_: return DLVHB
 return DLVH_SLOW_FIX - (m1 > 0) + dlvh_seg(0, m1) + dlvh_seg(m1, m2) + dlvh_seg(m1 + m2, K_ - m1 - m2)
DLVHS = max(dlvh_copy(a, b) for a in range(K2B + 1) for b in range(K2B + 1 - a) if (a, b) != (0, K2B))
BLOCKH = FLUSH + DONE_C + DLVHB
FINALH = FINH_DEST + FINH_SENT + DONE_C + DLVHS + MSG_LINK + (DLVHS - DLVHB)
check("dlvh: DLVHB 31 (8 persistent port regs, dlvh0-hoisted ST/HT/CC, K2B = 16), DLVHS 84, BLOCKH 52, FINALH 186",
   DLVHB == 31 and DLVHS == 84 and BLOCKH == FLUSH + DONE_C + DLVHB and FINALH == 186)
ZW = 16384
ZERO_PASS = 2 * (ZW // 32) + 2
STOP_CHECK = 5
NBAR_COUNT = 8 + 19 * G_GROUPS + 1
NBAR_ALLOW = 30000
PHASE_CTRL_G = 2 ** 20
KB_BUCKET = KB(D1 + D2 + HB)
KB64 = KB(64)
R2_UP_TIB = 682
CAP_H = (((R2_UP_TIB + 4095) * (2 ** 40)) // 512) * 16
def exp_neg_hi_top(m_val):
 if not m_val: return Fraction(1)
 k_h, mv = 0, m_val
 while mv > 1: mv /= 2; k_h += 1
 mlo = Fraction((mv.numerator << 200) // mv.denominator, 2 ** 200)
 eh = exp_interval(-mlo, 45)[1]
 for _ in range(k_h):
  eh = eh * eh
  eh = Fraction(cdiv(eh.numerator << 300, eh.denominator), 1 << 300)
 return eh
LOG2_POOLS = 15
N_POOLS = 1 << LOG2_POOLS
K_GEN = NODES // N_POOLS
B_MINE_PAIRS = 16384
PTS_MINE_CHAIN = 2 * B_MINE_PAIRS
NBATCH_MINE = (1 << 47) // PTS_MINE_CHAIN
PTS_PER_GEN_MINE = NBATCH_MINE * PTS_MINE_CHAIN
M_POOL = K_GEN * PTS_PER_GEN_MINE
N_CHAINS_POOL = K_GEN * NBATCH_MINE
PAIRS_PER_MINE_CHAIN = PTS_MINE_CHAIN * (PTS_MINE_CHAIN - 1) // 2
P_BUCKET_MAX = Fraction(6 * (1 << 63), r - 1)
P_SAME_CHAIN_POOL = N_CHAINS_POOL * PAIRS_PER_MINE_CHAIN * P_BUCKET_MAX
GAMMA_MINE_CROSS = M_POOL * (M_POOL - 1) // 2 - N_CHAINS_POOL * PAIRS_PER_MINE_CHAIN
P_BDAY_HIT_MIN = Fraction(1, 1 << 193) - Fraction(6, r - 1)
MU_POOL_LO = GAMMA_MINE_CROSS * P_BDAY_HIT_MIN
CS_POOL = GAMMA_MINE_CROSS * (2 * PTS_MINE_CHAIN * M_POOL) * (2 * P_BUCKET_MAX * P_BUCKET_MAX)
P_TRIPLE_POOL = Fraction(M_POOL * (M_POOL - 1) * (M_POOL - 2), 6) * (P_BUCKET_MAX * P_BUCKET_MAX)
POIS_LE7 = exp_neg_hi_top(MU_POOL_LO) * sum((MU_POOL_LO ** _k) / math.factorial(_k) for _k in range(8))
P_FAIL_POOL = POIS_LE7 + CS_POOL + P_TRIPLE_POOL + P_SAME_CHAIN_POOL
P_FAIL_ALL_POOLS = N_POOLS * P_FAIL_POOL
XMINE_PAIR_HOT = 963
XMINE_EMIT_PAIR = 2 * 16
XMINE_FIX_ALLOW = 20
GEN_MINE_NODE = GEN_MUL + NBATCH_MINE * (2 * GEN_MUL + 2 * COIN_SCALAR + INV + B_MINE_PAIRS * (XMINE_PAIR_HOT + XMINE_EMIT_PAIR + XMINE_FIX_ALLOW) + 64)
S_MINE_TOT = NODES * PTS_PER_GEN_MINE
ROUTE_MINE = S_MINE_TOT * 24 + N_HOME * ((1 << 32) // ZW * ZERO_PASS + 1) + 2 * N_POOLS * (MX // 2) * (2 * PTS_MINE_CHAIN * 1000 + NODES * MSG // M_G)
E_MINE = NODES * GEN_MINE_NODE + ROUTE_MINE
check("Phase-0 transversal birthday X-pair mining: N_POOLS >= M_G, P_FAIL < 2^-45, E_MINE itemized",
   N_POOLS >= M_G and ((p - 1) >> 63) + 1 == 1 << 193 and MU_POOL_LO > 63 and P_FAIL_ALL_POOLS < Fraction(1, 1 << 45),
   "p_fail 2^%.2f E_mine 2^%.4f" % (lg(P_FAIL_ALL_POOLS), lg(E_MINE)))
def success_prefix(q1p, qs=None):
 n0p = LEAVES * q1p * Q2
 Sp = NODES * n0p
 pool_n0 = K_GEN * n0p
 clust_n0 = M_G * n0p
 gamma = Fraction(Sp * Sp - N_POOLS * pool_n0 * pool_n0 - N_CLUST * clust_n0 * clust_n0, 2)
 D_LAB_SYM = TWO_BOUND * (3 * 653) * 2 ** (128 + PART)
 P1_GOOD = Fraction(6 * (r - 1), r * r) * (1 - D_WJ - D_LAB_SYM)
 mu_raw = gamma * P1
 mu = gamma * P1_GOOD
 n_y = q1p * L_CHAIN
 deg = 2 * LEAVES * (K_GEN * n_y + M_G * Q2) * Sp
 cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
 def exp_neg_hi(m_val):
  if not m_val: return Fraction(1)
  k_h, mv = 0, m_val
  while mv > 1: mv /= 2; k_h += 1
  mlo = Fraction((mv.numerator << 200) // mv.denominator, 2 ** 200)
  eh = exp_interval(-mlo, 45)[1]
  for _ in range(k_h):
   eh = eh * eh
   eh = Fraction(cdiv(eh.numerator << 300, eh.denominator), 1 << 300)
  return eh
 e_hi = exp_neg_hi(mu)
 d_batch, d_wj, d_lab = D_BATCH, D_WJ, D_LAB
 d_f = q1p * E_FIX / (DELTA_F - FIX_BATCH) + q1p * d_batch
 m_p2 = Sp * KB(D1 + D2) + 1
 v_p2 = m_p2 + 1 + 60 * M_G * n_y * KB(D1 + D2)
 d_p2 = v_p2 / Fraction(CAP_H - FIXED_LEAVES - m_p2) ** 2
 kb_m = KB(KEYBITS + D1 + D2)
 lam_int = RHO * (Sp + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
 sip = RHO * Fraction(5, 2) * (1 + 11 * (d_wj + d_lab)) * kb_m
 alpha_cross = 12 * Sp + PARTNER_PAIRS + CLUST_BATCH * CLUST_BATCH
 other_h = RHO * (Fraction(Sp) * Sp * KB(D1 + D2) * kb_m + alpha_cross * kb_m) + sip * Sp
 d_cl = other_h / (CAND_CAP - 4)
 other_g = RHO * Fraction(n0p) * Sp * kb_m * 2 + sip * n0p
 d_rq = other_g / (REQ_CAP - 8)
 other_v = RHO * (2 + 2 * 8) * (Fraction(Sp, SIDE) * Sp * KB(D2) * kb_m + alpha_cross * kb_m) + 18 * sip * Sp
 d_rl = other_v / (RELAY_CAP - 16)
 per_leaf = d_batch + d_f + d_p2 + d_cl + d_rq + d_rl
 rows_ = mu_raw * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf)
 g_def = NODES * q1p * d_batch
 glob = (P_FAIL_ALL_POOLS
         + NODES * (3 * q1p + 2 * NBATCH_MINE + 2) * Fraction(1, 2 ** 254)
         + N_POOLS * MX * Fraction(1, r)
         + N_CLUST * ((M_G * MX) * (M_G * MX - 1) // 2) * Fraction(6, r - 1)
         + g_def)
 lo = 1 - e_hi - cs - rows_ - glob
 return dict(mu=mu, lo=lo, cs=cs, rows=rows_, glob=glob, Sp=Sp, n0p=n0p, d_p2=d_p2, m_p2=m_p2,
    d_f=d_f, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def)
def grouped_full_lo(q):
 return success_prefix(q)["lo"]
Q1GD = 8700000000
Q1G, MQ1 = group_sizes(Q1GD)
N0G = LEAVES * Q1GD * Q2
SG_FULL = NODES * N0G
LO_RUN = grouped_full_lo(Q1GD)
GATE_BITS = 128
GATE_INSNS = 256
K_GATE = cdiv(THRESH.numerator * LO_RUN.denominator * (1 << GATE_BITS), THRESH.denominator * LO_RUN.numerator)
P_GATE = Fraction(K_GATE, 1 << GATE_BITS)
GATED_LO = P_GATE * LO_RUN
check("opened run q1 = %d has certified success lo_run >= 0.9979" % Q1GD,
   LO_RUN >= Fraction(9979, 10000), "lo_run %.12f" % (math.floor(LO_RUN * 10 ** 12) / 10 ** 12))
check("128-bit coordinator run gate K_GATE is least with (K_GATE / 2^128) * lo_run >= 0.3901",
   0 < K_GATE < (1 << GATE_BITS) and GATED_LO >= THRESH and Fraction(K_GATE - 1, 1 << GATE_BITS) * LO_RUN < THRESH
   and 2 * COIN_SCALAR + 4 + 8 <= GATE_INSNS,
   "p_gate %.9f gated_lo %.12f" % (float(P_GATE), math.floor(GATED_LO * 10 ** 12) / 10 ** 12))
check("power-law group boundaries m_g = floor(q1*(g/1536)^(28/32)) exact",
   all((MQ1[g] ** 32) * (G_GROUPS ** G_EXP) <= (Q1GD ** 32) * (g ** G_EXP) < ((MQ1[g] + 1) ** 32) * (G_GROUPS ** G_EXP) for g in range(1, G_GROUPS)))
TPRE = [success_prefix(m) for m in MQ1[:-1]]
def ceil_256(x):
 x = min(Fraction(1), max(Fraction(0), x))
 return Fraction(cdiv(x.numerator << 256, x.denominator), 1 << 256)
PHI = [Fraction(1)] + [ceil_256(1 - t["lo"]) for t in TPRE[1:]]
TFULL = success_prefix(Q1GD)
check("leaves S = 2^68*6*q1*q2 = 2^%.6f" % lg(SG_FULL), SG_FULL == NODES * 6 * Q1GD * Q2 == TFULL["Sp"])
check("mu = |Gamma|*P1 = %.9f (opened run cross-pool & cross-cluster pairs); gated_lo >= 0.3901" % float(TFULL["mu"]),
   TFULL["mu"] == Fraction(SG_FULL * SG_FULL - N_POOLS * (K_GEN * N0G) ** 2 - N_CLUST * (M_G * N0G) ** 2, 2) * Fraction(6 * (r - 1), r * r) * (1 - D_WJ - TWO_BOUND * 1959 * 2 ** (128 + PART)) and GATED_LO >= THRESH)
check("Chen-Stein b1+b2 < 2^-10 (transversal pools K_GEN = 2^53)", TFULL["cs"] < Fraction(1, 2 ** 10), "2^%.2f" % lg(TFULL["cs"]))
check("intruder rate per colliding pair < 2^-12.7", TFULL["lam_int"] < pow2_interval(Fraction(-127, 10))[0], "2^%.3f" % lg(TFULL["lam_int"]))
check("mean fixups per generator q1*E_FIX = 2^%.3f < 2^35; cap 2^52 > 2^17 times it plus 28 q2" % lg(Q1GD * E_FIX),
   Q1GD * E_FIX < pow2_interval(Fraction(35, 1))[0] and DELTA_F > 2 ** 17 * Q1GD * E_FIX + FIX_BATCH)
check("fixup meter rate (Markov, cap 2^52) < 2^-17.0", TFULL["d_f"] < pow2_interval(Fraction(-170, 10))[0], "2^%.3f" % lg(TFULL["d_f"]))
check("home-slot overflow rate < 2^-31.5", TFULL["d_p2"] < pow2_interval(Fraction(-315, 10))[0], "2^%.3f" % lg(TFULL["d_p2"]))
check("candidate cap < 2^-16, request cap < 2^-15.5, relay cap < 2^-20",
   TFULL["d_cl"] < pow2_interval(Fraction(-16, 1))[0] and TFULL["d_rq"] < pow2_interval(Fraction(-155, 10))[0]
   and TFULL["d_rl"] < pow2_interval(Fraction(-20, 1))[0],
   "2^%.3f 2^%.3f 2^%.3f" % (lg(TFULL["d_cl"]), lg(TFULL["d_rq"]), lg(TFULL["d_rl"])))
check("whole-run rows (incl. same-cluster X collisions) < 2^-25.5", TFULL["g_def"] < pow2_interval(Fraction(-255, 10))[0]
   and TFULL["glob"] < pow2_interval(Fraction(-255, 10))[0], "2^%.3f" % lg(TFULL["glob"]))
check("non-stopping bounds: PHI[0] = 1, nonincreasing, within [0, 1], >= exact 1 - lo",
   PHI[0] == 1 and all(0 <= x <= 1 for x in PHI) and all(PHI[g] >= PHI[g + 1] for g in range(G_GROUPS - 1))
   and all(PHI[g] >= min(Fraction(1), max(Fraction(0), 1 - TPRE[g]["lo"])) for g in range(1, G_GROUPS)),
   "PHI[1] %.4f PHI[16] %.4f PHI[31] %.4f" % (float(PHI[1]), float(PHI[16]), float(PHI[31])))
HOME_MAX = CAP_H
check("hash chain holds at most CAP_H < 2^59 entries", HOME_MAX == CAP_H < 2 ** 59, "2^%.4f" % lg(HOME_MAX))
PAIRM = KB_BUCKET + Fraction(1, r)
PAIRK2 = KB(D1 + D2 + 64) + Fraction(1, r)
VISIT = HCHG_VISIT + Fraction(HCHG_END, 2)
PJ = Fraction(2, 3) + 2 * D_WJ + Fraction(1, r - 1)
JROW_W = NODES * N0G * JPATH + NODES * N0G // 3 * SJOIN
check("E-charge per visit 4 (hit blocks out of line)", VISIT == 4)
SG = [NODES * LEAVES * q * Q2 for q in Q1G]
WG, WGW, SAMEX, WGJ, WGK2 = [], [], [], [], []
EXCH_ROW = ALLN * 4 * 2 ** 10
for _g, _q in enumerate(Q1G):
 _s = SG[_g]
 _leaf = NODES * _q * (fwd_count(Q2) + INV + INV_CALL + BATCH_SETUP) + NODES * y_group(_q) \
  + int(NODES * _q * bwd_count(Q2)) + int(NODES * _q * ((KEY_INC - EMIT_INC) * Q2 + KEY_FIX))
 _route = _s * (EMIT * K + FLUSH + DONE_C1) // K + _s * H2S_ARR + (_s // K) * H2S_FULL \
  + NODES * (FINAL * F1 + FIN_FIX) + _s * BLOCKH // K2B + N_REL * (FINALH * F2 + FIN_FIX) \
  + N_HOME * EXPECTH + N_REL * (EXPECT + RDONE) + NODES * (INIT3_BKT * F1 + 1)
 _gp = NODES * LEAVES * MQ1[_g] * Q2
 _samex = N_CLUST * _q * CLUST_BATCH * (CLUST_BATCH - 1) // 2
 _main = _s * _gp + _s * (_s - 1) // 2 - _samex
 _evis = _main * PAIRM
 _gfix = _s * HCHG_NEW + ALLN * STOP_CHECK
 WGK2.append(HCHG_VER * _main * PAIRK2)
 _group = _gfix + VISIT * _evis + WGK2[-1]
 _jkey = JPATH * _s * PJ + _s // 3 * SJOIN
 WGJ.append(_jkey)
 WG.append(_leaf + _route + _group + _jkey + EXCH_ROW)
 SAMEX.append(_samex)
 _wgroup = _gfix + _s * HOME_MAX * (HCHG_VISIT + HCHG_END + HCHG_VER + HCHG_LS1) + EXCH_ROW
 WGW.append(_leaf + _route + _wgroup)
check("group sizes sum to q1, each >= 1; records multiples of k1",
   sum(Q1G) == Q1GD and min(Q1G) > 0 and all(s % K == 0 for s in SG))
check("relay's 128-bit D: |finals' counts - messages| <= 2^118 per group, no wrap before last final",
   NODES * 2 ** 50 <= 2 ** 118 and 2 ** 126 + 2 ** 118 < 2 ** 127 and NODES < 2 ** 69
   and max(Q1G) * LEAVES * Q2 // K + 1 < 2 ** 50)
ZERO_INIT = N_HOME * (2 ** HB // ZW * ZERO_PASS + 1)
check("hash-head zeroing per home: zero16384 (BS = 6, 256-B blocks) over 2^47 words at 2^64 - 2^52",
   ZERO_PASS == 1026 and 2 ** HB % ZW == 0 and (2 ** 64 - 2 ** 52) % 256 == 0)
HR_INIT = N_HOME * (INIT3_BKT + 1) + N_HOME * 15
E_PRE = E_MINE + ALLN * NODE_INIT + NODES * (3 * cdiv(Q1GD, M_G) + 2) * COIN_SCALAR + N_POOLS * MX * (GEN_MULX + XSETUP) + NODES * MX * XSETUP \
 + NODES * XTBLD_ALL + NODES * Q2 * ENTRY_X + NODES * HEADER \
 + N_REL * (INIT3_BKT * F2 + 1 + TINIT_ENTRY * F2 + 1) + NODES * (INIT3_BKT * F1 + 2) + ZERO_INIT + HR_INIT + NODES * TBLD
_KB170 = KB(KEYBITS + D1 + D2)
_SIPT = RHO * Fraction(5, 2) * (1 + 11 * (D_WJ + D_LAB))
_ACROSS_FULL = 12 * SG_FULL + PARTNER_PAIRS + CLUST_BATCH * CLUST_BATCH
ELIST = RHO * (Fraction(SG_FULL) * SG_FULL * _KB170 + _ACROSS_FULL * _KB170) + _SIPT * _KB170 * SG_FULL
EK2_PROD = RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) * KB64 + _ACROSS_FULL * KB64) + _SIPT * KB64 * SG_FULL
EK2 = RHO * (_ACROSS_FULL * KB64) + _SIPT * KB64 * SG_FULL
SXK2 = PAIRM * sum(SAMEX)
SXV = VISIT * PAIRM * sum(SAMEX)
BRACKET = RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) + _ACROSS_FULL) + _SIPT * SG_FULL
DEFV = VISIT * D_BATCH * (HOME_MAX * SG_FULL + sum(SAMEX))
DEFK2 = HCHG_VER * 2 * D_BATCH * BRACKET + HCHG_VER * D_BATCH * HOME_MAX * SG_FULL
DEFL_PAIRS = 2 * D_BATCH * BRACKET + D_BATCH * HOME_MAX * SG_FULL
DEFL = HCHG_LS1 * DEFL_PAIRS
DEFJ = JPATH * D_BATCH * SG_FULL
E_FIX_WORK = NODES * Q1GD * E_FIX * max(FIX_EXTRA.values())
E_CAND_WORK = (ELIST + DEFL_PAIRS) * RESOLVE
E_REQ_WORK = (ELIST + DEFL_PAIRS) * (2 * REQ_SERVE)
E_SHARED_REST = E_PRE + NODES * GEN_MUL + 2 * ALLN * MSG + NBAR_ALLOW * BARRIER + ALLN * PHASE_CTRL_G + OUTPUT + G_GROUPS * 2 * ALLN * MSG
E_WORST_ONLY = NODES * DELTA_F * max(FIX_EXTRA.values()) + N_HOME * CAND_CAP * RESOLVE + NODES * REQ_CAP * REQ_SERVE
R_WORST_OPEN = E_SHARED_REST + E_WORST_ONLY + JROW_W + sum(WGW)
R_WORST = GATE_INSNS + R_WORST_OPEN
E_RUN = E_SHARED_REST + E_FIX_WORK + E_CAND_WORK + E_REQ_WORK + sum(PHI[_g] * WG[_g] for _g in range(G_GROUPS)) \
 + SXV + EK2 * HCHG_VER + ELIST * HCHG_LS1 + DEFV + DEFK2 + DEFL + DEFJ + SXK2 * HCHG_VER
E_TOTAL = GATE_INSNS + P_GATE * E_RUN
def claim6_of(x):
 c = Fraction(math.ceil(lg(x) * 10 ** 6), 10 ** 6)
 for _ in range(50):
  if pow2_interval(c)[0] >= x and pow2_interval(c - Fraction(1, 10 ** 6))[1] < x: return c
  c += Fraction(1, 10 ** 6) if pow2_interval(c)[0] < x else -Fraction(1, 10 ** 6)
 raise AssertionError("worst claim not certified")
WORST6 = claim6_of(R_WORST)
E6 = claim6_of(E_TOTAL)
check("E-total display 2^%d/10^6: least 6-decimal bound >= E_TOTAL" % (E6 * 10 ** 6),
   pow2_interval(E6)[0] >= E_TOTAL and pow2_interval(E6 - Fraction(1, 10 ** 6))[1] < E_TOTAL)
W_E = E_TOTAL / C
def claim15_of(w):
 getcontext().prec = 80
 start = (Decimal(w.numerator).ln() - Decimal(w.denominator).ln()) / Decimal(2).ln() * 10 ** 15
 c = Fraction(int(start.to_integral_value(rounding=ROUND_CEILING)), 10 ** 15)
 for _ in range(400):
  if pow2_interval(c)[0] >= w and pow2_interval(c - Fraction(1, 10 ** 15))[1] < w: return c
  c += Fraction(1, 10 ** 15) if pow2_interval(c)[0] < w else -Fraction(1, 10 ** 15)
 raise AssertionError("E-claim not certified")
CLAIM_E = claim15_of(W_E)
check("E-work below the worst case", E_TOTAL < R_WORST, "E 2^%.6f worst 2^%.6f" % (float(E6), float(WORST6)))
check("E-claim = least 15-decimal c with 2^c >= E/1000",
   pow2_interval(CLAIM_E)[0] >= W_E and pow2_interval(CLAIM_E - Fraction(1, 10 ** 15))[1] < W_E, "%d/10^15" % (CLAIM_E * 10 ** 15))
JROW_E = sum(PHI[_g] * WGJ[_g] for _g in range(G_GROUPS)) + DEFJ
check("key paths E-priced (JPATH 1, SJOIN 1) below every-leaf row",
   PJ > Fraction(2, 3) and D_WJ > Fraction(1, r - 1) and JROW_E < JROW_W, "2^%.4f < 2^%.4f" % (lg(JROW_E), lg(JROW_W)))
K2ROW_E = sum(PHI[_g] * WGK2[_g] for _g in range(G_GROUPS)) + SXK2 * HCHG_VER + EK2 * HCHG_VER
check("joint K3 hit counting (KB(148) + 1/r under PHI) below product row",
   PAIRK2 < KB(D1 + D2) * KB64 and K2ROW_E < EK2_PROD * HCHG_VER and PAIRK2 <= PAIRM,
   "2^%.4f < 2^%.4f" % (lg(K2ROW_E), lg(EK2_PROD * HCHG_VER)))
check("barrier allowance covers the grouped run", NBAR_ALLOW >= NBAR_COUNT, "%d of %d" % (NBAR_COUNT, NBAR_ALLOW))
PER_LEAF_OPEN = E_RUN / SG_FULL
PER_LEAF_E = PER_LEAF_OPEN
LOW, UP_BASE = 2 ** 52, 2 ** 64 - 2 ** 52
def valid_va(a): return 0 <= a < 2 ** 64 and (a >> 52) in (0, 2 ** 12 - 1)
TIB = 2 ** 40
AXYX_BASE = 64 * TIB
AXYX_HI_BASE = 576 * TIB
XT_PAIR44_BASE = 768 * TIB
XT_PAIR44_BYTES = 24 * 44 * TIB
XT_BASE = AXYX_BASE
XT_BYTES = 40 * 44 * TIB
check("beta tables in TTBR1 and 40 x 44-TiB X-tables in TTBR0 [64 TiB, 1,824 TiB) (< 2^51): valid 52-bit VAs, 2^43/2^42-aligned AXYX",
   BT_BASE == 2 ** 64 - 2 ** 51 and AXYX_BASE % (2 ** 43) == 0 and AXYX_HI_BASE == AXYX_BASE + MX * 32 * TIB and AXYX_HI_BASE % (2 ** 42) == 0
   and XT_PAIR44_BASE == AXYX_HI_BASE + MX * 12 * TIB and XT_PAIR44_BASE + XT_PAIR44_BYTES == 1824 * TIB < 2 ** 51
   and valid_va(AXYX_BASE) and valid_va(XT_PAIR44_BASE + XT_PAIR44_BYTES - 1)
   and all(valid_va(BT_BASE | ((w << 8) % 2 ** 64)) for w in (0, 1, 2 ** 43 - 1, 2 ** 64 - 1)))
YR_SLOT = 2 ** 31
YR_BASE = 32 * TIB + YR_SLOT
EB_END = 33 * TIB + 2 ** 34
EB_START = EB_END - Q2 * 512
GEN_LOW = [("LC table at 0 (8*d1)", F1 * 8),
   ("257 x 2-GiB Y-slots at 32 TiB + 2 GiB & chain blocks", (YB + 1) * YR_SLOT + YB * 2304),
   ("pair entries [EB_START, EB_END) at 33 TiB + 2^34 - 512*Q2", Q2 * 512),
   ("MMIO, code, constants, stack, granules, H, descriptors, X, scalars, phase-4 at 34 TiB", 16384 + 2 ** 20 + 8 * 2048 + MX * 6 * XT_SLOTS * 64 + 1024
    + MX * 64 + (MX + 2 * cdiv(Q1GD, M_G) + Q1GD + 1) * 32 + 1024 + REQ_CAP * 64 + 8 * 2 * 2 ** 20 * 64),
   ("40 x 44-TiB convex X-tables (1,760 TiB in [64 TiB, 1,824 TiB)) + Phase-0 64-KiB table", XT_BYTES + 2 * B_MINE_PAIRS * 64),
   ("record blocks, first half of d1, at 2^51", (F1 // 2) * 64 * K1)]
GEN_UP = [("record blocks, second half of d1, at 2^64 - 2^52", (F1 // 2) * 64 * K1), ("beta tables at 2^64 - 2^51", 2 ** BT_W << 8)]
REL_LOW = [("MMIO, code, stack, barriers, counters, blocks", 8192 + 2 ** 20 + 8 * 2048 + 4096), ("row LC & RC tables", 2 * F2 * 8),
   ("E2 and p tables", 2 * F2 * 8), ("phase-4 relay buffer", RELAY_CAP * 64)]
REL_UP = [("row record blocks at 2^64 - 2^52", F2 * 64 * K2B)]
R2_BYTES = CAP_H * 32
R2_UP_BYTES = R2_UP_TIB * TIB
R2_LOW_BYTES = R2_BYTES - R2_UP_BYTES
LINK_HR = 2 ** 52 - 2 ** 40
LINK_HR_MSG = LINK_HR + 0x1000
LINK_GEN = 34 * TIB
LINK_GEN_MSG = LINK_GEN + 0x1000
PORTY_GEN = LINK_GEN + 0xB000
NEXT_BASE = 2 ** 64 - 3 * 2 ** 50
HOME_LOW = [("R2 upper continuation [0, 4095 TiB)", R2_LOW_BYTES),
   ("MMIO link windows, code, stack, barriers, counters, join state, list, replies in [2^52 - 2^40, 2^52)", 8192 + 2 ** 20 + 8 * 2048 + 4096 + CAND_CAP * 528)]
HOME_UP = [("hash heads at 2^64 - 2^52 (HB = 47, 1024 TiB)", 2 ** HB * 8),
   ("hash next-array at 2^64 - 3*2^50 (2388.5 TiB)", CAP_H * 16),
   ("R2 initial 682-TiB segment at 2^64 - 682 TiB", R2_UP_BYTES)]
def disjoint_va(iv):
 s_iv = sorted(iv)
 return (all(0 <= a < b <= 2 ** 64 and valid_va(a) and valid_va(b - 1) and ((a < LOW) == (b - 1 < LOW)) for a, b in s_iv)
         and all(s_iv[i][1] <= s_iv[i + 1][0] for i in range(len(s_iv) - 1)))
GEN_IV = [(0, F1 * 8), (YR_BASE - 64, YR_BASE - 64 + GEN_LOW[1][1]), (EB_START, EB_END),
          (LINK_GEN, LINK_GEN + 4096), (LINK_GEN_MSG, LINK_GEN_MSG + 4096), (PORTY_GEN, PORTY_GEN + 512),
          (34 * TIB + 0x10000, 34 * TIB + 0x10000 + GEN_LOW[3][1]),
          (AXYX_BASE, XT_PAIR44_BASE + XT_PAIR44_BYTES), (2 ** 51, 2 ** 52), (UP_BASE, UP_BASE + 2 ** 51), (BT_BASE, 2 ** 64)]
REL_IV = [(0x10000, 0x10000 + REL_LOW[0][1]), (2 ** 45, 2 ** 45 + F2 * 8), (3 * 2 ** 45, 3 * 2 ** 45 + F2 * 8),
          (2 ** 47, 2 ** 47 + REL_LOW[2][1]), (2 ** 48 - 65536, 2 ** 48), (2 ** 48, 2 ** 48 + RELAY_CAP * 64),
          (LINK_HR, LINK_HR + 4096), (LINK_HR_MSG, LINK_HR_MSG + 4096), (UP_BASE, 2 ** 64)]
HOME_IV = [(0, R2_LOW_BYTES), (LINK_HR, LINK_HR + 4096), (LINK_HR_MSG, LINK_HR_MSG + 4096),
           (LINK_HR + 0x10000, 2 ** 52 - 16384), (2 ** 52 - 16384, 2 ** 52),
           (UP_BASE, UP_BASE + 2 ** HB * 8), (NEXT_BASE, NEXT_BASE + CAP_H * 16), (2 ** 64 - R2_UP_BYTES, 2 ** 64)]
check("generator: LC [0, 32T), Y-buf [32T+2G-64, 33T), entries [EB_START, EB_END), MMIO+misc [34T, 59T), X-tables [59T, 2^51), d1<2^41 [2^51, 2^52)",
   F1 * 8 == 32 * TIB and YR_BASE - 64 > F1 * 8 and YROUNDS * YPC * 512 + 64 <= YR_SLOT
   and YR_BASE - 64 + GEN_LOW[1][1] <= 33 * TIB and Q2 * 512 < 2 ** 34 and EB_END <= 34 * TIB
   and 34 * TIB + 0x10000 + GEN_LOW[3][1] <= AXYX_BASE and XT_PAIR44_BASE + XT_PAIR44_BYTES <= 2 ** 51
   and (F1 // 2) * 64 * K1 == 2 ** 51 and sum(b for _, b in GEN_UP) == 2 ** 52 and disjoint_va(GEN_IV))
check("relay and home pairwise VA/MMIO disjointness (R2 [-682T, 4095T) disjoint from LINK 4095T & CNT/ST, heads at 2^64 - 2^52)",
   sum(b for _, b in REL_LOW) < 512 * TIB and sum(b for _, b in REL_UP) == 2 ** 52
   and R2_LOW_BYTES == LINK_HR and HOME_LOW[1][1] <= 2 ** 40 and NEXT_BASE + CAP_H * 16 <= 2 ** 64 - R2_UP_BYTES
   and disjoint_va(REL_IV) and disjoint_va(HOME_IV))
check("physical addresses: every role fits 56-bit PA (FEAT_D128)",
   all(sum(b for _, b in lo + up) <= 2 ** 56 for lo, up in ((GEN_LOW, GEN_UP), (REL_LOW, REL_UP), (HOME_LOW, HOME_UP))))
check("tag: node 68 + (2*l + s) 26 (2*Q2 < 2^26) + (k - 1) 34 (q1 <= 2^34) = 128 <= 128 bits in w2, w3; YROUNDS*33 >= L_CHAIN",
   Q1GD <= 2 ** 34 and 2 * Q2 < 2 ** 26 and YROUNDS * YPC >= L_CHAIN)
PEAK = NODES * sum(b for _, b in GEN_LOW + GEN_UP) + N_REL * sum(b for _, b in REL_LOW + REL_UP) + N_HOME * sum(b for _, b in HOME_LOW + HOME_UP)
_KH = KB(D1 + D2)
COV = Fraction(SG_FULL * 6 * M_G * L_CHAIN, r - 1) * _KH ** 2
D_P2 = (TFULL["m_p2"] + COV + 60 * M_G * Q1GD * L_CHAIN * _KH) / Fraction(CAP_H - FIXED_LEAVES - TFULL["m_p2"]) ** 2
check("home slot: Chebyshev d_p2 < 2^-31.5 with FIXED_LEAVES = 2*M_G*6*Q2",
   COV < 1 and D_P2 <= TFULL["d_p2"] and CAP_H % K2B == 0 and CAP_H > TFULL["m_p2"] + FIXED_LEAVES + 2 ** 38)
def least_at_or_above(d, x, scale):
 return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x
def least_grid_up(x, scale):
 k = math.ceil(lg(x) * scale)
 for _ in range(20):
  d = Fraction(k, scale)
  if least_at_or_above(d, x, scale): return d
  k += 1 if pow2_interval(d)[0] < x else -1
 raise AssertionError("grid bound not certified")
DECLARED = dict(time_log2=CLAIM_E, worst_log2=WORST6, memory_log2_bytes=least_grid_up(PEAK, 10 ** 4),
    preprocessing_log2=least_grid_up(Fraction(E_PRE, C), 10 ** 4), success=THRESH,
    per_leaf=Fraction(cdiv(PER_LEAF_E.numerator * 10 ** 4, PER_LEAF_E.denominator), 10 ** 4))
MULX, SQRX, ADDX = 105 + 4, 80 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_TEST, GEN_CTRL, NEGT, RED = 256 * 2 * 13 + 64, 1600, 30, 9
X_ENTRY = MADD_ + 3 * SQRX + 11 * MULX + 2 * ADDX + 2 * NEGT + 17 * RED + 32
check("generator allowance GEN_MUL / GEN_MULX covers complete Jacobian double-and-add + conversions",
   GEN_ITEM + GEN_TEST + GEN_CTRL + X_ENTRY + 4 <= GEN_MULX)
check("chain allowance: affine add/dbl <= CHAIN_ADD; X descriptors (56 slots) <= XSETUP; entry X words <= ENTRY_X",
   INV + 3 * MULX + 2 * SQRX + 10 * ADDX + 4 < CHAIN_ADD and 56 * (2 * MULX + ADDX + 4) < XSETUP and 10 <= ENTRY_X)
BARRIER_MSGS = 2 * (N_HOME + NODES + 2 * (N_REL - 1))
BARRIER_MSG = 7 + DONE + DLV + RCV + WAKE
BARRIER_ITEM = ALLN * (3 * 8 + 20 + 1) + BARRIER_MSGS * BARRIER_MSG
check("batch set-up (meter, rho, kernel calls, yround/ring loop) <= 1280, phase control < 2^17, barrier > 7x, broadcast <= 40",
   64 + cdiv(4 * YROUNDS, YB) <= BATCH_SETUP and G_GROUPS * (12 * 33 + 16 + 4) < PHASE_CTRL_G and BARRIER >= 7 * BARRIER_ITEM and 9 + MSG_LINK <= MSG)
DLVQ, DLVQ_DROP, REQH_C, FWDQ, SENDR, FWDR = 10, 7, 17, 11, 10, 15
PH4_HOPS = {k: v + DONE + RCV + WAKE for k, v in dict(req_home_relay=cdiv(REQH_C, 2) + DLVQ, req_relay_gen=FWDQ + DLVQ,
    rep_gen_relay=SENDR + DLVQ, rep_relay_home=FWDR + DLV).items()}
_h = list(PH4_HOPS.values())
CAND_LEAVES = INV + 5 * SQRX + 18 * MULX + 20 * ADDX
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 9 * 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("phase-4 hops <= MSG, candidate resolve > 2x itemized, program delivery (<= 2^17 B) < NODE_INIT",
   max(_h) <= MSG and RESOLVE == 2 * CAND_ITEM + 1 and (2 ** 17 // 32) * 2 * MSG < NODE_INIT)
print("q1 %d  q2 %d  n0 %d  S 2^%.6f  mu %.9f  lo_run %.12f  p_gate %.9f  success >= %.12f" % (
 Q1GD, Q2, N0G, lg(SG_FULL), float(TFULL["mu"]), math.floor(LO_RUN * 10 ** 12) / 10 ** 12, float(P_GATE), math.floor(GATED_LO * 10 ** 12) / 10 ** 12))
print("E_total 2^%.6f instructions  W_E 2^%.15f  E-claim %d/10^15  per leaf %.7f declared %s  R_worst 2^%.6f" % (
 float(E6), lg(W_E), CLAIM_E * 10 ** 15, float(PER_LEAF_E), format(float(DECLARED["per_leaf"]), ".4f"), float(WORST6)))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
 format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK),
 format(float(DECLARED["preprocessing_log2"]), ".4f"), lg(Fraction(E_PRE, C))))
STATED = dict(time_log2=Fraction(123493382876856411, 1000000000000000), worst_log2=Fraction(179732087, 10 ** 6),
     memory_log2_bytes=Fraction(1369996, 10 ** 4), preprocessing_log2=Fraction(1172303, 10 ** 4),
     success=Fraction(3901, 10 ** 4), per_leaf=Fraction(740795, 10 ** 4))
_HERE = os.path.dirname(os.path.abspath(__file__))
_CL = json.load(open(os.path.join(_HERE, "submission.json")))["claims"][0]
_ST, _SC = _CL["statement"], _CL["score"]
_PF = open(os.path.join(_HERE, "proof.md")).read()
check("stated figures (literals here, K_GATE, proof.md §5, submission.json's statement and score) equal the certified ones",
   K_GATE == 133015161165470784583994612803444211703 and CAP_H == 164136470183936 and STATED == DECLARED
   and Fraction(int(_SC["numerator"]), int(_SC["denominator"])) == STATED["time_log2"]
   and format(CAP_H, ",") in _PF and all(f"2^{b}" in _PF for b in ("52.9588", "52.0661", "52.9996"))
    and all(s in _PF for s in ("[64 TiB, 1,824 TiB)", "16233/16 = 1,014.5625", "1081/16 = 67.5625", "513/4 = 128.25"))
    and "1,888" not in _PF and "48-TiB" not in _PF and "16217" not in _PF and "48-TiB" not in json.dumps(json.load(open(os.path.join(_HERE, "submission.json"))))
   and all(re.search(re.escape(t) + r"(?![0-9])", _ST) for t in (
    "2^123.493382876856411", "2^179.732087", "2^136.9996", "2^117.2303", "74.0795", "0.3901", "133015161165470784583994612803444211703")))
print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
 sys.exit(1 if FAILED else 0)
