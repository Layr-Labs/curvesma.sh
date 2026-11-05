"""Exact certificate for proof.md (52-bit-VA flipped-grid machine, E-work ledger, success bound, claim, memory)."""
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
F1, F2 = 2 ** 42, 2 ** 41
N_REL, N_HOME = F1, F1 * F2
ALLN = NODES + N_REL + N_HOME
SIDE = F1
Q2 = 2 ** 25
MX = 256
L_CHAIN = Q2 // MX
D1, D2, D3, D4 = 42, 41, 64, 23
KEYBITS = 64 + D4
PART = 64 + D1
KBF = 3000
K1, K2B = 16, 32
THRESH = Fraction(390000001, 10 ** 9)

LEAVES = 6
FWD_U = 32
FWD_INC, FWD_FIX = 78 + Fraction(4, FWD_U) + Fraction(1, MX), 9

def fwd_count(q2):
 return FWD_FIX + 78 * q2 + 4 * (q2 // FWD_U) + (q2 // MX) + (q2 % FWD_U)

BWD_INC, BWD_FIX, BWD_YSW = Fraction(2289, 2), 11, 6

def bwd_count(q2):
 return BWD_FIX + BWD_INC * q2 + BWD_YSW * (-(-q2 // MX) - 1)

# Sites and fixups per X-pair (2 incidences: C1FWD8 + BWDP; proof.md 3.4)
FIX_SITES_PAIR = dict(SHX=27, FDX=6, F5=4, FS=3, FADDW=2, BLK=14)
FIX_EXTRA = dict(SHX=3, FDX=3, F5=3, FS=4, FADDW=4, BLK=5)
EMIT = 5
EMIT_INC = LEAVES * EMIT
KEY_INC, KEY_FIX = Fraction(265, 2), 4
KEY_SIGN = 1

K = K1
FLUSH = 10
H2S_BLOCK_DROP = 9
FLUSH2_X, FLUSH2_D = 10, 8
LOOP8 = 2
JPATH, SJOIN = 1, 1
INV = 22549
TINIT_ENTRY = 3
FIN_DEST, FIN_SENT, FIN_HELD, FIN_FIX = 15, 4, 2, 2
FINH_DEST, FINH_SENT = 16, 5

H2S_ARR, H2S_ONE = 6, 1
H2S_FULL = 6 + 2 * H2S_ONE
FLUSH1_X, FIN1_X = 7, 14
H2S_EMPTY = 7
H2S_ARR_FIX = 14 + 2 * H2S_ONE

def h2s_count(n):
 if n == K: return H2S_FULL + H2S_ARR * K
 if n == 0: return H2S_EMPTY
 return 14 + H2S_ONE * (1 + (n % 2 == 0)) + LOOP8 * ((n - 1) // 8) + H2S_ARR * n

check("h2s: 6 per record (LD64B + ADD per two slots), +1 per lone slot; full block 8 fixed, other <= 16 + 2*(k1/8 - 1)",
   h2s_count(K) == 6 * K + 8 and max(h2s_count(n) - H2S_ARR * n for n in range(K)) == H2S_ARR_FIX + LOOP8 * (K // 8 - 1)
   and H2S_EMPTY <= H2S_ARR_FIX)

INIT3_BKT = 4
check("key pass: 132.5 = 2 signs * (6 LDP + 26 sums + hoisted ADD + 3 * (6 index + emit 5)) + TBZ per two incidences",
   KEY_INC == 2 * (6 + 26 + KEY_SIGN + 3 * (6 + EMIT)) + Fraction(1, 2))
check("bwd: two BWDP bodies (2288 per 2 incidences = 1144 each) per TST, B.NE (+1/2); per Y switch 6; 11 fixed; q2 % 4 == 0",
   BWD_INC == 1144 + Fraction(1, 2) and Q2 % 4 == 0 and BWD_FIX == 11 and BWD_YSW == 6
   and bwd_count(Q2) == BWD_FIX + BWD_INC * Q2 + BWD_YSW * (Q2 // MX - 1))
check("fwd: C1FWD8 156 per 2 incidences (78 each) + 4 per 32-pass + 1 per Y switch, 9 fixed; Q2 a multiple of 32 and MX",
   Q2 % FWD_U == 0 and Q2 % MX == 0 and fwd_count(Q2) == FWD_FIX + FWD_INC * Q2 and fwd_count(4) == 9 + 4 * 79)
check("j = 1 and j = 2 paths (exact-carry ADDS, ADC; CMP, B.HS; CMN, B.HS; MVN of sigma; emit 5) exceed j = 0 by 1; SJOIN 1",
   (4 + 2 + 1 + EMIT) - (6 + EMIT) == JPATH and SJOIN == 1)
check("emit: AND or BIC (8*d1), LDADD, STP K2,K3, STP w2,x18, TBZ of bit 9 for slot k1 = 16",
   EMIT == 1 + 1 + 2 + 1)

DONE = 2
RCV, WAKE, EXPECT = 1, 8, 8
DONE_C = DONE + RCV + WAKE
RCV1, FIN_CNT, H2S_DONE = 4, 6, 5
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
YR_FIX, YR_BLK, YB = 22556, 43093, 256
YROUNDS, YPC = 3972, 33
CHAIN_ADD = 2 ** 15
CHAIN_SET = 2 * GEN_MUL + 17 * CHAIN_ADD

def y_group(q):
 return YROUNDS * (-(-q // YB) * YR_FIX + q * YR_BLK) + q * CHAIN_SET

# Unified 8x32b X-tables (proof.md 3.1): 5 tables per X (2a^3, 2a [+a^4], a^2, y_X; and per X-pair -sigma [+pi-7sigma], A2-A1)
XT_SLOTS = 8
XT_W = 32
XT_PER_X = 5
XTBLD_ONE = 4 + XT_SLOTS * 10 + 20 * (XT_SLOTS * 2 ** XT_W)
XTBLD_ALL = MX * XT_PER_X * XTBLD_ONE
XSETUP = 2 ** 16
ENTRY_X = 16

COIN_SCALAR = 64
BATCH_SETUP = 160
INV_CALL = 16
NODE_INIT = 2 ** 20
HEADER = 72
BT_BASE, BT_W = 0xfff8 << 48, 43
TBLD = 4 + 7 * 10 + 20 * (6 * 2 ** BT_W + 2 ** (256 - 5 * BT_W))
MSG = 40

CAND_CAP = 2 ** 19
RESOLVE = 1491011
REQ_CAP = 2 ** 34
REQ_SERVE = 2 ** 10 + GEN_MUL
RELAY_CAP = 2 ** 38
BARRIER = 2 ** 11 * ALLN
OUTPUT = 2 ** 20
DELTA_F = 2 ** 49
FIX_BATCH = sum(FIX_SITES_PAIR.values()) * (Q2 // 2)

P1 = Fraction(6 * (r - 1) + 1, r * r)
TWO_BOUND = Fraction(2, r - 1)

SHX_BAD = 2 ** 64 * 2 ** 33
FOLD_BAD = 2 ** 66 * C0
SQ_BAD = C0 * C0 + 2 * C0
FS_C2_BAD = C0
QP = Q2 // 2
SITES = [
 ("fwd: row shifts of acc_(j-1) (uniform on F_p* via rho)", 3 * QP, SHX_BAD, 1),
 ("fwd: fold of acc_j", QP, FOLD_BAD, 1),
 ("bwd: row shifts of R_j", 3 * QP, SHX_BAD, 1),
 ("bwd: fold of inv_e = R_j acc_(j-1)", QP, FOLD_BAD, 12),
 ("bwd: fold of R_(j-1) = R_j e_j", QP, FOLD_BAD, 1),
 ("bwd: row shifts of D_2 = B - A_2", 3 * QP, SHX_BAD, 6),
 ("bwd: fold of 1/D_1 = D_2 inv_e", QP, FOLD_BAD, 6),
 ("bwd: fold of A = inv y_Y (C2)", Q2, FOLD_BAD, 6),
 ("bwd: FS of u = A - L and FADDW of v = A + L (C2)", 2 * Q2, FS_C2_BAD, 3),
 ("bwd: square folds of u^2 and v^2", 2 * Q2, SQ_BAD, 6),
 ("bwd: row shifts of E, F, G (lookup sums fold exactly)", 9 * Q2, SHX_BAD, 8),
]
SITES_PAIR = dict(SHX=27, FDX=6, F5=4, FS_FA=4)
D_SITES = sum(Fraction(n * bad * mult, r - 1) for _, n, bad, mult in SITES)

FIX_ROW = C0 * 2 ** 192
FIX_FOLD = (2 ** 66 + 1) * 2 ** 128 * C0
FIX_SQ = 2 ** 193
FIX_LK8 = 28 * C0 * 2 ** 128
FIX_LK9 = 36 * C0 * 2 ** 128
FIX_BLK = 30 * C0 * 2 ** 128
FIX_SUB = C0 * 2 ** 192
FIXUPS = [
 ("fwd: lookup-sum fold of e_j (9 terms)", "BLK", QP, FIX_LK9, 12),
 ("fwd: row shifts of acc_(j-1)", "SHX", 3 * QP, FIX_ROW, 1),
 ("fwd: fold of acc_j", "FDX", QP, FIX_FOLD, 1),
 ("bwd: row shifts of R_j", "SHX", 3 * QP, FIX_ROW, 1),
 ("bwd: fold of inv_e", "FDX", QP, FIX_FOLD, 12),
 ("bwd: fold of R_(j-1)", "FDX", QP, FIX_FOLD, 1),
 ("bwd: FS of D_2 = B - A_2", "FS", QP, FIX_SUB, 6),
 ("bwd: row shifts of D_2", "SHX", 3 * QP, FIX_ROW, 6),
 ("bwd: fold of 1/D_1", "FDX", QP, FIX_FOLD, 6),
 ("bwd: lookup-sum fold of 1/D_2 (9 terms)", "BLK", QP, FIX_LK9, 6),
 ("bwd: fold of A = inv y_Y", "FDX", Q2, FIX_FOLD, 6),
 ("bwd: lookup-sum fold of L = inv y_X (8 terms)", "BLK", Q2, FIX_LK8, 6),
 ("bwd: FS of u = A - L", "FS", Q2, FIX_SUB, 3),
 ("bwd: FADDW of v = A + L", "FADDW", Q2, FIX_SUB, 3),
 ("bwd: square folds of u^2 and v^2", "F5", 2 * Q2, FIX_SQ, 6),
 ("bwd: window-sum folds of beta u^2 and beta v^2", "BLK", 2 * Q2, FIX_BLK, 6),
 ("bwd: lookup-sum fold of E (9 terms)", "BLK", Q2, FIX_LK9, 8),
 ("bwd: lookup-sum folds of F, G (8 terms)", "BLK", 2 * Q2, FIX_LK8, 8),
 ("bwd: row shifts of E, F, G", "SHX", 9 * Q2, FIX_ROW, 8),
]
E_FIX = sum(Fraction(n * bad * mult, r - 1) for _, _, n, bad, mult in FIXUPS)
KB = lambda fixed: Fraction(6 * 2 ** (256 - fixed) * (KBF + 1), (r - 1) * KBF)

RHO = Fraction(r, r - 1) ** 2
FIXED_LEAVES = 2 * LEAVES * Q2
PARTNER_PAIRS = LEAVES * Q2 * (LEAVES * Q2 - 1)

D_YGEN = Fraction(YROUNDS * (YPC + 2 * 16 + 2), r - 1)
D_BATCH = D_SITES + Q2 * Fraction(6, r - 1) + D_YGEN
D_WJ = TWO_BOUND * 23 * 2 ** 192
D_LAB = TWO_BOUND * 1305 * 2 ** (128 + PART)

FA, FB = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16), int("e4437ed6010e88286f547fa90abfe4c3", 16)
BOX = (p >> PART) + 1
WIDE = 654 * (p >> 128) + 4
check("lattice lemma: basis (B, A - B), (A - B, -A) of {(v, u): u = beta v mod p}, |det| = p; box factor (106-bit) < 1 + 1/3000",
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
check("site table per X-pair: 27 SHX, 6 FDX, 4 F5, 4 C2 FS/FADDW (lookup and beta window sums exact)",
   sum(n for _, n, b, _ in SITES if b == SHX_BAD) == SITES_PAIR["SHX"] * QP
   and sum(n for _, n, b, _ in SITES if b == FOLD_BAD) == SITES_PAIR["FDX"] * QP
   and sum(n for _, n, b, _ in SITES if b == SQ_BAD) == SITES_PAIR["F5"] * QP
   and sum(n for _, n, b, _ in SITES if b == FS_C2_BAD) == SITES_PAIR["FS_FA"] * QP and len(SITES) == 11)
FIX_RUNS = {m: sum(n for _, ms, n, _, _ in FIXUPS if ms == m) for m in FIX_SITES_PAIR}
check("fixup table per X-pair: 27 SHX, 6 FDX, 4 F5, 3 FS, 2 FADDW, 14 BLK (28 per incidence)",
   FIX_RUNS == {m: FIX_SITES_PAIR[m] * QP for m in FIX_SITES_PAIR} and sum(FIX_SITES_PAIR.values()) == 56)
check("a call's fixup word (C0 per fixup, at most 28 q2 fixups) stays below 2^63", FIX_BATCH * C0 < 2 ** 63,
   "2^%.3f" % lg(FIX_BATCH * C0))
check("label-error rate per leaf (106-bit parts) < 2^-10.65", D_LAB < pow2_interval(Fraction(-1065, 100))[0], "2^%.3f" % lg(D_LAB))
check("batch defect rate (site table, D = 0, degenerate chain) < 2^-126.05", D_BATCH < pow2_interval(Fraction(-12605, 100))[0],
   "2^%.3f" % lg(D_BATCH))

G_GROUPS = 160
def group_sizes(q1):
 q, rem = divmod(q1, G_GROUPS)
 qs = [q + 1 if g < rem else q for g in range(G_GROUPS)]
 ms = [0] * (G_GROUPS + 1)
 for _g in range(G_GROUPS): ms[_g + 1] = ms[_g] + qs[_g]
 return qs, ms

HB = 47
HCHG_NEW, HCHG_VISIT, HCHG_END, HCHG_VER = 4, 4, 0, 5
HCHG_LS0, HCHG_LSA = 2, 1
HCHG_LS1 = 4 + HCHG_LSA
check("join per record 4, visit 4, out-of-line hit blocks (.pushsection) HCHG_END 0, d4 check 5, listing 4 + 1",
   HCHG_NEW == 4 and HCHG_VISIT == 4 and HCHG_END == 0 and HCHG_VER == 5 and HCHG_LS1 == 5)

DLVC, DLVS1 = 2, 2
DLVH_FAST = 15 + 2 * DLVS1
DLVH_SEG, DLVH_SLOW_FIX = 7, 24
def dlvh_seg(n, m):
 if m == 0: return 1
 head = n % 2; rest = m - head
 return DLVH_SEG + head + DLVC * m + DLVS1 * (head + rest % 2) + LOOP8 * (rest // 8)
def dlvh_copy(m1, m2, K_=K2B):
 if m1 == 0 and m2 == K_: return DLVH_FAST + DLVC * K_
 return DLVH_SLOW_FIX - (m1 > 0) + dlvh_seg(1, m1) + dlvh_seg(1 + m1, m2) + dlvh_seg(1 + m1 + m2, K_ - m1 - m2)
DLVHB = dlvh_copy(0, K2B)
DLVHS = max(dlvh_copy(a, b) for a in range(K2B + 1) for b in range(K2B + 1 - a) if (a, b) != (0, K2B))
BLOCKH = FLUSH + DONE_C + DLVHB
FINALH = FINH_DEST + FINH_SENT + DONE_C + DLVHS + MSG_LINK + (DLVHS - DLVHB)
check("dlvh: DLVHB 83, DLVHS 129, BLOCKH 104, FINALH 207",
   DLVHB == 19 + 2 * K2B and DLVHS == 129 and BLOCKH == FLUSH + DONE_C + DLVHB)

ZW = 16384
ZERO_PASS = 2 * (ZW // 32) + 2
STOP_CHECK = 5
NBAR_COUNT = 4 + 17 * G_GROUPS + 1
NBAR_ALLOW = 3072
PHASE_CTRL_G = 2 ** 16

KB_BUCKET = KB(D1 + D2 + HB)
KB64 = KB(64)
_MB = (2 ** 126 * KB(D1 + D2)); CAP_H = (int(_MB * 1.64) // (2 ** 36) + 1) * 2 ** 36

def success_prefix(q1p, qs):
 n0p = LEAVES * q1p * Q2
 Sp = NODES * n0p
 gamma = Fraction(Sp * Sp - NODES * n0p * n0p, 2)
 mu = gamma * P1
 n_y = q1p * L_CHAIN
 deg = 2 * LEAVES * (n_y + Q2) * Sp
 cs = gamma * deg * (P1 * P1 + Fraction(36, (r - 1) ** 2))
 mu_lo = Fraction((mu.numerator << 200) // mu.denominator, 2 ** 200) if mu else Fraction(0)
 e_hi = exp_interval(-mu_lo, 40)[1]
 d_batch, d_wj, d_lab = D_BATCH, D_WJ, D_LAB
 d_f = q1p * E_FIX / (DELTA_F - FIX_BATCH) + q1p * d_batch
 m_p2 = Sp * KB(D1 + D2) + 1
 v_p2 = m_p2 + 1 + 60 * n_y * KB(D1 + D2)
 d_p2 = v_p2 / Fraction(CAP_H - FIXED_LEAVES - m_p2) ** 2
 kb_m = KB(KEYBITS + D1 + D2)
 lam_int = RHO * (Sp + FIXED_LEAVES) * KB(D1 + D2 + D3 + D4)
 sip = RHO * Fraction(5, 2) * (1 + 11 * (d_wj + d_lab)) * kb_m
 other_h = RHO * (Fraction(Sp) * Sp * KB(D1 + D2) * kb_m + (2 * Sp + PARTNER_PAIRS) * kb_m) + sip * Sp
 d_cl = other_h / (CAND_CAP - 4)
 other_g = RHO * Fraction(n0p) * Sp * kb_m * 2 + sip * n0p
 d_rq = other_g / (REQ_CAP - 8)
 other_v = RHO * (2 + 2 * 8) * (Fraction(Sp, SIDE) * Sp * KB(D2) * kb_m + (2 * Sp + PARTNER_PAIRS) * kb_m) + 18 * sip * Sp
 d_rl = other_v / (RELAY_CAP - 16)
 per_leaf = d_batch + d_wj + d_f + d_lab + d_p2 + d_cl + d_rq + d_rl
 rows_ = mu * (1 + Fraction(1, r)) * (Fraction(1, r) + lam_int + 2 * per_leaf)
 g_def = NODES * q1p * d_batch
 glob = NODES * (3 * q1p + MX + 1) * Fraction(1, 2 ** 254) + NODES * MX * Fraction(1, r) + g_def
 lo = 1 - e_hi - cs - rows_ - glob
 return dict(mu=mu, lo=lo, cs=cs, rows=rows_, glob=glob, Sp=Sp, n0p=n0p, d_p2=d_p2, m_p2=m_p2,
    d_f=d_f, lam_int=lam_int, d_cl=d_cl, d_rq=d_rq, d_rl=d_rl, g_def=g_def)

def grouped_full_lo(q):
 qs, _ = group_sizes(q)
 return success_prefix(q, qs)["lo"]

lo_q, hi_q = 2 ** 31, 2 ** 32
check("grouped q1 search starts at 2^31, where the full grouped bound fails", grouped_full_lo(lo_q) < THRESH)
while hi_q - lo_q > 1:
 mid = (lo_q + hi_q) // 2
 if grouped_full_lo(mid) >= THRESH: hi_q = mid
 else: lo_q = mid
Q1G, MQ1 = group_sizes(hi_q)
Q1GD = hi_q
N0G = LEAVES * Q1GD * Q2
SG_FULL = NODES * N0G
check("grouped q1 = %d is least with the full grouped bound >= 0.39" % Q1GD,
   grouped_full_lo(Q1GD) >= THRESH and grouped_full_lo(Q1GD - 1) < THRESH,
   "lo %.12f" % (math.floor(grouped_full_lo(Q1GD) * 10 ** 12) / 10 ** 12))
check("grouped bound locally increasing at the least q1",
   all(grouped_full_lo(q) <= grouped_full_lo(q + 1) for q in range(Q1GD - 2, Q1GD + 3)))

TPRE = [success_prefix(m, Q1G[:g]) for g, m in enumerate(MQ1[:-1])]
PHI = [Fraction(1)] + [min(Fraction(1), max(Fraction(0), 1 - t["lo"])) for t in TPRE[1:]]
TFULL = success_prefix(Q1GD, Q1G)
check("leaves S = 2^68*6*q1*q2 = 2^%.6f" % lg(SG_FULL), SG_FULL == NODES * 6 * Q1GD * Q2 == TFULL["Sp"])
check("mu = |Gamma|*P1 = %.9f (cross-node pairs); lo >= 0.39" % float(TFULL["mu"]),
   TFULL["mu"] == Fraction(SG_FULL * SG_FULL - NODES * N0G * N0G, 2) * P1 and TFULL["lo"] >= THRESH)
check("Chen-Stein b1+b2 < 2^-75", TFULL["cs"] < Fraction(1, 2 ** 75), "2^%.2f" % lg(TFULL["cs"]))
check("intruder rate per colliding pair < 2^-12.7", TFULL["lam_int"] < pow2_interval(Fraction(-127, 10))[0], "2^%.3f" % lg(TFULL["lam_int"]))
check("mean fixups per generator q1*E_FIX = 2^%.3f < 2^32.68; cap 2^49 > 2^16 times it plus 28 q2" % lg(Q1GD * E_FIX),
   Q1GD * E_FIX < pow2_interval(Fraction(3268, 100))[0] and DELTA_F > 2 ** 16 * Q1GD * E_FIX + FIX_BATCH)
check("fixup meter rate (Markov, cap 2^49) < 2^-16.32", TFULL["d_f"] < pow2_interval(Fraction(-1632, 100))[0], "2^%.3f" % lg(TFULL["d_f"]))
check("home-slot overflow rate < 2^-32.6", TFULL["d_p2"] < pow2_interval(Fraction(-326, 10))[0], "2^%.3f" % lg(TFULL["d_p2"]))
check("candidate cap < 2^-12, request cap < 2^-14, relay cap < 2^-20",
   TFULL["d_cl"] < pow2_interval(Fraction(-12, 1))[0] and TFULL["d_rq"] < pow2_interval(Fraction(-14, 1))[0]
   and TFULL["d_rl"] < pow2_interval(Fraction(-20, 1))[0],
   "2^%.3f 2^%.3f 2^%.3f" % (lg(TFULL["d_cl"]), lg(TFULL["d_rq"]), lg(TFULL["d_rl"])))
check("whole-run rows < 2^-26.5", TFULL["g_def"] < pow2_interval(Fraction(-265, 10))[0]
   and TFULL["glob"] < pow2_interval(Fraction(-265, 10))[0], "2^%.3f" % lg(TFULL["glob"]))
check("non-stopping bounds: PHI[0] = 1, nonincreasing, within [0, 1]",
   PHI[0] == 1 and all(0 <= x <= 1 for x in PHI) and all(PHI[g] >= PHI[g + 1] for g in range(G_GROUPS - 1)),
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
 _samex = NODES * _q * (6 * Q2) * (6 * Q2 - 1) // 2
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
check("relay's 128-bit D: |finals' counts - messages| <= 2^117 per group, no wrap before last final",
   NODES * 2 ** 49 <= 2 ** 117 and 2 ** 126 + 2 ** 117 < 2 ** 127 and NODES < 2 ** 69
   and max(Q1G) * LEAVES * Q2 // K + 1 < 2 ** 49)

ZERO_INIT = N_HOME * (2 ** HB // ZW * ZERO_PASS + 1)
check("hash-head zeroing per home: zero16384 (BS = 6, 256-B blocks) over 2^47 words at 2^64 - 2^50",
   ZERO_PASS == 1026 and 2 ** HB % ZW == 0 and (2 ** 64 - 2 ** 50) % 256 == 0)
HR_INIT = N_HOME * (INIT3_BKT + 1)
E_PRE = ALLN * NODE_INIT + NODES * (3 * Q1GD + MX + 1) * COIN_SCALAR + NODES * MX * (GEN_MULX + XSETUP) \
 + NODES * XTBLD_ALL + NODES * Q2 * ENTRY_X + NODES * HEADER \
 + N_REL * (INIT3_BKT * F2 + 1 + TINIT_ENTRY * F2 + 1) + NODES * (INIT3_BKT * F1 + 2) + ZERO_INIT + HR_INIT + NODES * TBLD
E_WORST_REST = E_PRE + NODES * GEN_MUL + 2 * ALLN * MSG \
 + NODES * DELTA_F * max(FIX_EXTRA.values()) \
 + N_HOME * CAND_CAP * RESOLVE + NODES * REQ_CAP * REQ_SERVE \
 + NBAR_ALLOW * BARRIER + ALLN * PHASE_CTRL_G + OUTPUT + G_GROUPS * 2 * ALLN * MSG

_KB170 = KB(KEYBITS + D1 + D2)
_SIPT = RHO * Fraction(5, 2) * (1 + 11 * (D_WJ + D_LAB))
ELIST = RHO * (Fraction(SG_FULL) * SG_FULL * _KB170 + (2 * SG_FULL + PARTNER_PAIRS) * _KB170) + _SIPT * _KB170 * SG_FULL
EK2_PROD = RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) * KB64 + (2 * SG_FULL + PARTNER_PAIRS) * KB64) + _SIPT * KB64 * SG_FULL
EK2 = RHO * ((2 * SG_FULL + PARTNER_PAIRS) * KB64) + _SIPT * KB64 * SG_FULL
SXK2 = PAIRM * sum(SAMEX)
SXV = VISIT * PAIRM * sum(SAMEX)
BRACKET = RHO * (Fraction(SG_FULL) * SG_FULL * KB(D1 + D2) + 2 * SG_FULL + PARTNER_PAIRS) + _SIPT * SG_FULL
DEFV = VISIT * D_BATCH * (HOME_MAX * SG_FULL + sum(SAMEX))
DEFK2 = HCHG_VER * 2 * D_BATCH * BRACKET + HCHG_VER * D_BATCH * HOME_MAX * SG_FULL
DEFL = HCHG_LS1 * 2 * D_BATCH * BRACKET + HCHG_LS1 * D_BATCH * HOME_MAX * SG_FULL
DEFJ = JPATH * D_BATCH * SG_FULL

E_TOTAL = E_WORST_REST + sum(PHI[_g] * WG[_g] for _g in range(G_GROUPS)) + SXV + EK2 * HCHG_VER \
 + ELIST * HCHG_LS1 + DEFV + DEFK2 + DEFL + DEFJ + SXK2 * HCHG_VER
R_WORST = E_WORST_REST + JROW_W + sum(WGW)

def claim6_of(x):
 c = Fraction(math.ceil(lg(x) * 10 ** 6), 10 ** 6)
 for _ in range(50):
  if pow2_interval(c)[0] >= x and pow2_interval(c - Fraction(1, 10 ** 6))[1] < x: return c
  c += Fraction(1, 10 ** 6) if pow2_interval(c)[0] < x else -Fraction(1, 10 ** 6)
 raise AssertionError("worst claim not certified")

WORST6 = claim6_of(R_WORST)
DECLARED_WORST6 = Fraction(176807579, 10 ** 6)
check("disclosed worst case 2^%d/10^6: least 6-decimal bound >= R_WORST" % (WORST6 * 10 ** 6),
   WORST6 == DECLARED_WORST6 and pow2_interval(WORST6)[0] >= R_WORST and pow2_interval(WORST6 - Fraction(1, 10 ** 6))[1] < R_WORST)
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
check("E-work below the worst case", E_TOTAL < R_WORST, "E 2^%.6f worst 2^%.6f" % (float(E6), float(DECLARED_WORST6)))
check("E-claim = least 15-decimal c with 2^c >= E/1000",
   pow2_interval(CLAIM_E)[0] >= W_E and pow2_interval(CLAIM_E - Fraction(1, 10 ** 15))[1] < W_E, "%d/10^15" % (CLAIM_E * 10 ** 15))
JROW_E = sum(PHI[_g] * WGJ[_g] for _g in range(G_GROUPS)) + DEFJ
check("key paths E-priced (JPATH 1, SJOIN 1) below every-leaf row",
   PJ > Fraction(2, 3) and D_WJ > Fraction(1, r - 1) and JROW_E < JROW_W, "2^%.4f < 2^%.4f" % (lg(JROW_E), lg(JROW_W)))
K2ROW_E = sum(PHI[_g] * WGK2[_g] for _g in range(G_GROUPS)) + SXK2 * HCHG_VER + EK2 * HCHG_VER
check("joint K3 hit counting (KB(147) + 1/r under PHI) below product row",
   PAIRK2 < KB(D1 + D2) * KB64 and K2ROW_E < EK2_PROD * HCHG_VER and PAIRK2 <= PAIRM,
   "2^%.4f < 2^%.4f" % (lg(K2ROW_E), lg(EK2_PROD * HCHG_VER)))
check("barrier allowance covers the grouped run", NBAR_ALLOW >= NBAR_COUNT, "%d of %d" % (NBAR_COUNT, NBAR_ALLOW))
PER_LEAF_E = E_TOTAL / SG_FULL

# 52-bit VA memory layout per role
LOW, UP_BASE = 2 ** 52, 2 ** 64 - 2 ** 52
def valid_va(a): return 0 <= a < 2 ** 64 and (a >> 52) in (0, 2 ** 12 - 1)
XT_BYTES = MX * XT_PER_X * (2 ** XT_W << 8)
XT_BASE = 2 ** 48
check("beta tables in TTBR1 and X-tables [2^48, 2^48 + 1280*2^40) in TTBR0: valid 52-bit VAs, 64-byte aligned",
   BT_BASE == 2 ** 64 - 2 ** 51 and XT_BASE + XT_BYTES <= 2 ** 51 and valid_va(XT_BASE) and valid_va(XT_BASE + XT_BYTES - 1)
   and all(valid_va(BT_BASE | ((w << 8) % 2 ** 64)) for w in (0, 1, 2 ** 43 - 1, 2 ** 64 - 1)))

GEN_LOW = [("code, constants, stack", 2 ** 20), ("barrier & counter granules", 8 * 2048), ("header H, descriptors", MX * XT_PER_X * XT_SLOTS * 64 + 1024),
   ("X points (x, y)", MX * 64), ("scalars (b, a_j, chain c, s, rho)", (MX + 3 * Q1GD + 1) * 32),
   ("pair entries (2^25 * 512 B)", Q2 * 512), ("Y records & chain blocks (YB=256)", YB * (YROUNDS * YPC * 512 + 2304)),
   ("message blocks, request/reply buffers", 128 + REQ_CAP * 64 + 8 * 2 * 2 ** 20 * 64),
   ("LC table at 0 (8*d1)", F1 * 8), ("X-tables (1280 * 1 TiB) at 2^48", XT_BYTES),
   ("record blocks, first half of d1, at 2^51", (F1 // 2) * 64 * K1)]
GEN_UP = [("record blocks, second half of d1, at 2^64 - 2^52", (F1 // 2) * 64 * K1), ("beta tables at 2^64 - 2^51", 2 ** BT_W << 8)]
REL_LOW = [("code, stack, barriers, counters, blocks", 2 ** 20 + 8 * 2048 + 128), ("row LC & RC tables", 2 * F2 * 8),
   ("E2 and p tables", 2 * F2 * 8), ("phase-4 relay buffer", RELAY_CAP * 64)]
REL_UP = [("row record blocks at 2^64 - 2^52", F2 * 64 * K2B)]
R2_BYTES = CAP_H * 32
HOME_LOW = [("code, stack, barriers, counters, join state, list, replies", 2 ** 20 + 8 * 2048 + 64 + CAND_CAP * 528),
   ("R2 at 2^50 + 32 (CAP_H records)", R2_BYTES + 32)]
HOME_UP = [("hash next-array at 2^64 - 2^52", CAP_H * 16), ("hash heads at 2^64 - 2^50 (HB = 47)", 2 ** HB * 8)]

_GO = sum(b for n, b in GEN_LOW if not n.startswith(("LC table", "X-tables", "record blocks")))
check("generator: LC [0, 2^45), other regions from 2^46 below 2^48, X-tables [2^48, 2^51), d1<2^41 blocks [2^51, 2^52), upper half exact",
   F1 * 8 == 2 ** 45 and 2 ** 46 + _GO + 2 ** 35 <= XT_BASE and XT_BASE + XT_BYTES <= 2 ** 51
   and (F1 // 2) * 64 * K1 == 2 ** 51 and sum(b for _, b in GEN_UP) == 2 ** 52, "other 2^%.4f, XT 2^%.4f" % (lg(_GO), lg(XT_BYTES)))
check("relay and home fit their two 52-bit VA halves (home heads at 2^64 - 2^50 disjoint from next-array)",
   sum(b for _, b in REL_LOW) < 2 ** 48 - 8192 and sum(b for _, b in REL_UP) == 2 ** 52
   and 2 ** 50 + 32 + R2_BYTES <= 2 ** 52 and HOME_LOW[0][1] <= 2 ** 50
   and UP_BASE + CAP_H * 16 <= 2 ** 64 - 2 ** 50 and valid_va(2 ** 64 - 2 ** 50))
check("physical addresses: every role fits 56-bit PA (FEAT_D128)",
   all(sum(b for _, b in lo + up) <= 2 ** 56 for lo, up in ((GEN_LOW, GEN_UP), (REL_LOW, REL_UP), (HOME_LOW, HOME_UP))))
check("tag: node 68 + k 32 (q1 < 2^32) + l 25 + code 3 = 128 bits in w2, w3; YROUNDS*33 >= L_CHAIN",
   Q1GD < 2 ** 32 and Q2 == 2 ** 25 and YROUNDS * YPC >= L_CHAIN)

PEAK = NODES * sum(b for _, b in GEN_LOW + GEN_UP) + N_REL * sum(b for _, b in REL_LOW + REL_UP) + N_HOME * sum(b for _, b in HOME_LOW + HOME_UP)
_KH = KB(D1 + D2)
COV = Fraction(SG_FULL * 6 * L_CHAIN, r - 1) * _KH ** 2
D_P2 = (TFULL["m_p2"] + COV + 60 * Q1GD * L_CHAIN * _KH) / Fraction(CAP_H - FIXED_LEAVES - TFULL["m_p2"]) ** 2
check("home slot: Chebyshev d_p2 < 2^-32.6 with FIXED_LEAVES = 12*Q2",
   COV < 1 and D_P2 <= TFULL["d_p2"] and CAP_H % K2B == 0 and CAP_H > TFULL["m_p2"] + 2 ** 34)

def least_at_or_above(d, x, scale):
 return pow2_interval(d)[0] >= x and pow2_interval(d - Fraction(1, scale))[1] < x

def least_grid_up(x, scale):
 k = math.ceil(lg(x) * scale)
 for _ in range(20):
  d = Fraction(k, scale)
  if least_at_or_above(d, x, scale): return d
  k += 1 if pow2_interval(d)[0] < x else -1
 raise AssertionError("grid bound not certified")

DECLARED = dict(time_log2=CLAIM_E, worst_log2=DECLARED_WORST6, memory_log2_bytes=least_grid_up(PEAK, 10 ** 4),
    preprocessing_log2=least_grid_up(Fraction(E_PRE, C), 10 ** 4), success=THRESH,
    per_leaf=Fraction(cdiv(PER_LEAF_E.numerator * 10 ** 4, PER_LEAF_E.denominator), 10 ** 4))

MULX, SQRX, ADDX = 105 + 4, 80 + 3, 9 + 2 * 9
DBL, MADD_ = 2 * MULX + 5 * SQRX + 10 * ADDX, 8 * MULX + 3 * SQRX + 10 * ADDX
GEN_ITEM = 256 * (DBL + MADD_) + 2 * (INV + 3 * MULX + 2 * SQRX + 10 * ADDX) + 64
GEN_TEST, GEN_CTRL, NEGT, RED = 256 * 2 * 13 + 64, 1600, 30, 9
X_ENTRY = MADD_ + 3 * SQRX + 11 * MULX + 2 * ADDX + 2 * NEGT + 17 * RED + 32
check("generator allowance GEN_MUL / GEN_MULX covers complete Jacobian double-and-add + conversions",
   GEN_ITEM + GEN_TEST + GEN_CTRL + X_ENTRY + 4 <= GEN_MULX)
check("chain allowance: affine add/dbl <= CHAIN_ADD; X descriptors (40 slots) <= XSETUP; entry X words <= ENTRY_X",
   INV + 3 * MULX + 2 * SQRX + 10 * ADDX + 4 < CHAIN_ADD and 40 * (2 * MULX + ADDX + 4) < XSETUP and 10 <= ENTRY_X)

BARRIER_MSGS = 2 * (N_HOME + NODES + 2 * (N_REL - 1))
BARRIER_MSG = 7 + DONE + DLV + RCV + WAKE
BARRIER_ITEM = ALLN * (3 * 8 + 20 + 1) + BARRIER_MSGS * BARRIER_MSG
check("batch set-up (meter, rho, kernel calls, yround loop) <= 160, phase control < 2^16, barrier > 7x, broadcast <= 40",
   44 + cdiv(4 * YROUNDS, YB) <= BATCH_SETUP and G_GROUPS * (12 * 33 + 8 + 4) < PHASE_CTRL_G and BARRIER >= 7 * BARRIER_ITEM and 9 + MSG_LINK <= MSG)

DLVQ, DLVQ_DROP, REQH_C, FWDQ, SENDR, FWDR = 10, 7, 17, 11, 10, 15
PH4_HOPS = {k: v + DONE + RCV + WAKE for k, v in dict(req_home_relay=cdiv(REQH_C, 2) + DLVQ, req_relay_gen=FWDQ + DLVQ,
    rep_gen_relay=SENDR + DLVQ, rep_relay_home=FWDR + DLV).items()}
_h = list(PH4_HOPS.values())
CAND_LEAVES = INV + 5 * SQRX + 18 * MULX + 20 * ADDX
CAND_ITEM = 4 * 8 * 2 * MSG + 2 * CAND_LEAVES + 9 * 600 + 384 * 200 + 2 * 256 * (DBL + MADD_) // 2 + INV
check("phase-4 hops <= MSG, candidate resolve > 2x itemized, program delivery (<= 2^17 B) < NODE_INIT",
   max(_h) <= MSG and RESOLVE == 2 * CAND_ITEM + 1 and (2 ** 17 // 32) * 2 * MSG < NODE_INIT)

print("q1 %d  q2 2^25  n0 %d  S 2^%.6f  mu %.9f  success >= %.12f" % (
 Q1GD, N0G, lg(SG_FULL), float(TFULL["mu"]), math.floor(TFULL["lo"] * 10 ** 12) / 10 ** 12))
print("E_total 2^%.6f instructions  W_E 2^%.10f  E-claim %d/10^15  per leaf %.7f declared %s  R_worst 2^%.6f" % (
 float(E6), lg(W_E), CLAIM_E * 10 ** 15, float(PER_LEAF_E), format(float(DECLARED["per_leaf"]), ".4f"), float(DECLARED_WORST6)))
print("declared: memory_log2_bytes %s (exact %.7f)  preprocessing_log2 %s (exact %.7f)" % (
 format(float(DECLARED["memory_log2_bytes"]), ".4f"), lg(PEAK),
 format(float(DECLARED["preprocessing_log2"]), ".4f"), lg(Fraction(E_PRE, C))))

STATED = dict(time_log2=Fraction(124448260242728181, 10 ** 15), worst_log2=Fraction(176807579, 10 ** 6),
     memory_log2_bytes=Fraction(1352301, 10 ** 4), preprocessing_log2=Fraction(1160470, 10 ** 4),
     success=Fraction(390000001, 10 ** 9), per_leaf=Fraction(2098230, 10 ** 4))
_CL = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "submission.json")))["claims"][0]
_ST, _SC = _CL["statement"], _CL["score"]
check("stated figures (literals here, submission.json's statement and score) equal the certified ones",
   STATED == DECLARED and Fraction(int(_SC["numerator"]), int(_SC["denominator"])) == STATED["time_log2"]
   and all(re.search(re.escape(t) + r"(?![0-9])", _ST) for t in (
    "2^124.448260242728181", "2^176.807579", "2^135.2301", "2^116.0470", "209.8230", "0.390000001")))

print("FAIL: " + ", ".join(FAILED) if FAILED else "PASS")
if __name__ == "__main__":
 sys.exit(1 if FAILED else 0)
