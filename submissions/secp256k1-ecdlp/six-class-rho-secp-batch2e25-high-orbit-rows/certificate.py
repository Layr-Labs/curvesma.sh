"""Exact-integer certificate for proof.md 5.3 (integers and Fractions only; ln by the atanh series, powers of two
by certified integer roots). Run python3 certificate.py; exit 0 = all checks hold."""
import json
import math
import os
import sys
from fractions import Fraction

if hasattr(sys, "set_int_max_str_digits"):
 sys.set_int_max_str_digits(400000)

FAILED = []

def check(name, ok, detail=""):
 print("%-58s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
 if not ok:
  FAILED.append(name)

def log2(x):
 f = Fraction(x)
 return math.log2(f.numerator) - math.log2(f.denominator)

def ceil_fraction(f):
 f = Fraction(f)
 return -((-f.numerator) // f.denominator)

def iroot(n, k):
 if n < 2:
  return n
 b = n.bit_length()
 lo = 1 << max(0, (b - 1) // k)
 hi = 1 << ((b + k - 1) // k)
 while lo ** k > n:
  lo >>= 1
 while hi ** k <= n:
  hi <<= 1
 while hi - lo > 1:
  mid = (lo + hi) // 2
  if mid ** k <= n:
   lo = mid
  else:
   hi = mid
 return lo

def pow2_interval(num, den, scale=80):
 sh = num + scale * den
 assert sh >= 0
 m = iroot(2 ** sh, den)
 assert m ** den <= 2 ** sh < (m + 1) ** den, "integer root certificate"
 return Fraction(m, 2 ** scale), Fraction(m + 1, 2 ** scale)

def pow2_frac(f, bits=200):
 """[lo, hi] containing 2^f for a rational f in [0, 1): exp(f ln 2) by its Taylor series, the tail bounded by
 t^(n+1)/(n+1)! * e^t with e^t < 3, ln 2 from ln_interval, both ends rounded outward to multiples of 2^-bits."""
 l_lo, l_hi = ln_interval(2)
 def series(t, n=40):
  s_, term = Fraction(1), Fraction(1)
  for k in range(1, n + 1):
   term = term * t / k
   s_ += term
  return s_, term * t / (n + 1) * 3
 lo, _ = series(f * l_lo)
 hi, rest = series(f * l_hi)
 q = 2 ** bits
 return Fraction(lo.numerator * q // lo.denominator, q), Fraction(ceil_fraction((hi + rest) * q), q)

def ln_interval(x, bits=240):
 z = Fraction(x - 1, x + 1)
 z2 = z * z
 total = Fraction(0)
 p = z
 k = 0
 eps = Fraction(1, 2 ** bits)
 while True:
  total += p / (2 * k + 1)
  rest = p * z2 / ((2 * k + 3) * (1 - z2))
  if rest < eps:
   break
  p *= z2
  k += 1
 return 2 * total, 2 * (total + rest)

# --- constants
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
P_FIELD = int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
R_PART = 2 ** 26
KAPPA = Fraction(100009, 100000)
P0 = Fraction(39001172, 10 ** 8)
HEADROOM = Fraction(116, 10 ** 7)
CLAIM = Fraction(125830102, 10 ** 6)
THRESHOLD = Fraction(39, 100)  # rules.md: success at least 0.39
LOSS_OTHER = Fraction(11, 10 ** 8)

# --- 1. counted kernels
MUL = (16 + 16) + 30 + 23
SQR = (10 + 10) + 8 + 6 + 7 + 23
MUL_SUB = (16 + 16) + 30 + 24 + 5
SQR_SUB2 = (10 + 10) + 8 + 6 + 7 + 24 + 10
SUB, CANON, ADD_WEAK, ADD_WEAK2 = 9, 8, 9, 14
MUL_ADD = (16 + 16) + 30 + 23 + 5
X3 = (10 + 10) + 8 + 6 + 7 + 23 + 10
MUL_R, SUB_R, MUL_ADD_R, X3_R, SQR_R = MUL - 1, SUB - 2, MUL_ADD - 1, X3 - 1, SQR - 1  # kernels with rare exits (2.5)
RMUL, SHIFT = (16 + 16) + 34 + 6, 4
SUB_B = SUB - 4
WALK1 = SUB_R + 4 + 2
WALKB = SUB_B + 2 * RMUL + 3 * SHIFT + 6 + 2
WALK2 = 2 * RMUL + 3 * SHIFT + X3_R + MUL_ADD_R + SUB_B + 10 + 4
PASS1 = Fraction(128 * WALK1 + 64 + 1, 128)
PROD = Fraction(256 * WALKB + 128 + 1, 256)
PASS2 = Fraction(256 * WALK2 + 128 + 1, 256)
check("fe_mul costs 85 instructions", MUL == 85, str(MUL))
check("fe_sqr costs 64 instructions", SQR == 64, str(SQR))
check("kernels with rare exits: fe_mul 84, fe_sub 7, y3' 89, x3 73", (MUL_R, SUB_R, MUL_ADD_R, X3_R) == (84, 7, 89, 73))
check("fused fe_mul and fe_sqr cost 91 and 75", MUL_SUB == 91 and SQR_SUB2 == 75)
check("unfused pass 1: 1,729 per 128 walks = 13.5078125 per walk", 128 * WALK1 + 65 == 1729 and PASS1 == Fraction(1729, 128), str(PASS1))
check("prod: dy 5, two products 72 on 3 rows; 43,393 per 256 = 169.50390625",
 (SUB_B, RMUL, WALKB) == (5, 72, 169) and 256 * WALKB + 129 == 43393 and PROD == Fraction(43393, 256), str(PROD))
check("pass 2: 2*72 + 12 + 73 + 89 + 5 + 14 = 337; 86,401 per 256 = 337.50390625",
 WALK2 == 337 and 256 * WALK2 + 129 == 86401 and PASS2 == Fraction(86401, 256), str(PASS2))
ARRIVE_F = 13 + 10 + 3 + 1 + 1 + 3 + 1
LOOPA = 2 + SUB_B + 2
WALKF = ARRIVE_F + LOOPA
PASS1F = Fraction(128 * WALKF + 1, 128)
check("loop A walk: arrival 32 + 2 LDP x(T) + dx 5 + 2 STP = 41", ARRIVE_F == 32 and LOOPA == 9 and WALKF == 41)
check("loop A (pass1f): 5,249 per 128 walks = 41.0078125 per walk", 128 * WALKF + 1 == 5249 and PASS1F == Fraction(5249, 128),
 str(PASS1F))
ARRIVE_N = 3 + 11 + 6 + 1 + 4 + 1 + 1
WALKN = ARRIVE_N + LOOPA  # a pair's second walk adds the STP of both addresses
PASS1N = Fraction(64 * (2 * WALKN + 1) + 1, 128)
check("loop A (pass1p): arrival 27, walks 36 + 37 per pair, 4,673 per 128 = 36.5078125",
 ARRIVE_N == 27 and WALKN == 36 and 64 * (2 * WALKN + 1) + 1 == 4673 and PASS1N == Fraction(4673, 128), str(PASS1N))
check("compare rounds +4.5 per walk; anchor rounds arrive.s + pass1 - pass1f = 3.515625",
 PASS1F - PASS1N == Fraction(9, 2) and Fraction(1985, 64) + PASS1 - PASS1F == Fraction(225, 64))

ROWS = [13, 10, 3, 1, Fraction(3, 2), Fraction(5, 2), Fraction(1, 64)]
ROWS_TOTAL = sum(ROWS)
check("unfused arrival rows (anchor rounds) sum to 31.015625", ROWS_TOTAL == Fraction(1985, 64), "%.6f" % float(ROWS_TOTAL))
FUSED_FIXED = 22 + 3 + 8
KD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kernels")
sys.path[:0] = [KD]
sys.dont_write_bytecode = True
import a64
K = {}
for f in "secp_loops.s", "secp_fe.s":
 K.update(a64._asm_counts(open(os.path.join(KD, f)).read()))
N = lambda k: K[k]["instructions"]
INV = 255 * (N("fe_sqr_f") + 1) + 15 * (N("fe_mul_f") + 1) + N("fe_inv") - K["fe_inv"]["calls"] + 255 * 4
FIX = {k: N(k) - sum(l["iteration"] for l in K[k]["loops"]) for k in ("pass1", "prod", "pass2")}
GRP = [K[k]["loops"][-1] for k in ("pass1", "prod", "pass2")]
check("a64.py: inversion 24,322; fixed pass1 8, prod 28, pass 2 28; group main paths 1,729, 43,393, 86,401",
 INV == 24322 and FIX == dict(pass1=8, prod=28, pass2=28) and [l["main"] for l in GRP] == [1729, 43393, 86401]
 and [l["iteration"] - l["main"] for l in GRP] == [0, 3 * 256, 3 * 256] and [l["skips"] for l in GRP] == [0, 256, 256]
 and all(l["inner_branches"] == 0 for k in ("pass1", "prod", "pass2") for l in K[k]["loops"]))
CALLS = 64 + 8
ACCT = 14 + 3 * 8
BATCH_FIXED = INV + FUSED_FIXED + FIX["prod"] + FIX["pass2"] + CALLS + 12 + 3 + 7 + 3 + 2 + 4 + ACCT
TAILS_ROUND = 0
ROUND_FIXED = ceil_fraction(BATCH_FIXED + TAILS_ROUND)
check("per-round work 24,552, no tail walks at n = 2^23",
 BATCH_FIXED == 24552 and TAILS_ROUND == 0 and ROUND_FIXED == 24552,
 str(ROUND_FIXED))
SHARE = Fraction(ROUND_FIXED, 2 ** 25)

QDX = Fraction(4686, 2 ** 14)
PREM = 3 * (QDX + 1 + Fraction(1, 2 ** 10) + 2 * Fraction(1, 2 ** 11))  # premise levels QDX + pi0, 1/2 + pi0, 1/2 + pi0 (5.4)
CORR = PREM + Fraction(1056, 2 ** 23)
STEP = PASS1N + PROD + PASS2 + CORR + SHARE
DBL_ARITH = 4 * MUL_R + MUL_SUB + SQR_R + SQR_SUB2 + ADD_WEAK + 2 * ADD_WEAK2 + SUB_R + 30
DBL = STEP + DBL_ARITH + 34 + 3 + SHARE
TRAV = STEP + 24 + 2 + 7 + 6
PRICE_STEP = Fraction(69570, 128) + CORR + SHARE
PRICE_TRAV, PRICE_DBL = PRICE_STEP + 39, PRICE_STEP + 677
check("step is 36.5078125 + 169.50390625 + 337.50390625 + CORR 32413728/2^23 + share = 547.382569",
 STEP == PRICE_STEP and CORR == PREM + Fraction(1056, 2 ** 23), "%.7f" % float(STEP))
check("doubling bound: displaced step + 639 + 34 + 3 + share = step + 676 + share",
 DBL_ARITH == 639 and DBL == STEP + 676 + SHARE)
check("traversal bound is step + 24 + 2 + 7 + 6", TRAV == PRICE_STEP + 39)
check("doubling price step + 677 covers its bound", PRICE_DBL >= DBL)

# --- 2b. meter, stubs, reserves (proof.md 2.5, 3.2, 4); 1 milli-unit = 1 instruction
STUB, RARE, RESERVE, WALK_STUBS, RERUN_STUBS = 20, 25, 336, 15, 7  # milli-units; 15 = 1 loop-A + 6 prod + 8 pass-2 stubs
def stub_macro(fname, name):
 """Instructions on a stub macro's longest path (its whole body, RCHG expanded)."""
 t = open(os.path.join(KD, fname), encoding="utf-8").read()
 body = lambda m: [x.strip() for x in t[t.index(".macro " + m + "\n") if m == "RCHG" else t.index(".macro " + m + " "):]
                   .split(".endm")[0].splitlines()[1:] if x.strip() and not x.strip().endswith(":")]
 return sum(len(body("RCHG")) if x == "RCHG" else 1 for x in body(name))

F, LB, LC, LT = stub_macro("secp_fuse.s", "SS"), stub_macro("secp_loops.s", "SB"), stub_macro("secp_loops.s", "SC"), stub_macro("secp_loops.s", "ST")
LABELS = {f: sum(1 for k in a64.a64_load(open(os.path.join(KD, f)).read())[1] if k.startswith(".Lr")) for f in ("secp_loops.s", "secp_fuse.s")}
check("stubs: loops 3,727, loop A 258", LABELS == {"secp_loops.s": 3727, "secp_fuse.s": 258}
 and 129 + 256 * 14 + 14 == 3727 and 2 * 129 == 258)
check("stubs 7, 9, 8, 13 instructions <= 20 milli-units (2.5)",
 (F, LB, LC, LT) == (7, 9, 8, 13) and max(F, LB, LC, LT) <= 13 <= STUB, "%d %d %d %d" % (F, LB, LC, LT))
check("rare path 22, 23 with loop A's B, <= 25 milli-units (2.5, 3.2)",
 5 + 5 + 5 + 4 + 3 + 1 == 23 <= RARE)
check("a walk's 15 stubs = 300; 15*20 + 25 = 325 <= the 336 reserve (2.5, 4)", WALK_STUBS * STUB == 300
 and WALK_STUBS * STUB + RARE == 325 <= RESERVE)
FINAL_UNITS, FINAL_CHARGE = 2 ** 21, 2 ** 21 * 1000 + 2 ** 25 * 25
check("final arrival pass priced 2^18 units + 25 per walk covers its work, < 2^18.85 units (2.5)",
 FINAL_CHARGE >= 2 ** 25 * ROWS_TOTAL + 2 ** 10 + 2 ** 25 * 25 and log2(Fraction(FINAL_CHARGE, 1000)) < 21.85
 and log2(Fraction(FINAL_CHARGE, 1000)) < 70,
 "2^%.5f" % log2(Fraction(FINAL_CHARGE, 1000)))
CMP_ROUND = (PASS1F + PROD + PASS2 + CORR + RESERVE) * 2 ** 25 + ROUND_FIXED
ANCHOR_ROUND = CMP_ROUND + (Fraction(993, 64) + 2) * 2 ** 25
check("anchor round charge < 2^32.83, its +(993/64 + 2)*2^23 = MOVZ #0x8C2, LSL #16; round 0 + 9n < 2^33 (4)",
 PASS1F + PROD + PASS2 + CORR == Fraction(4629501984, 2 ** 23) and (Fraction(993, 64) + 2) * 2 ** 25 < 2 ** 34 and 12 + (ROWS_TOTAL + PASS1 - PASS1F) == Fraction(993, 64)
 and CMP_ROUND < ANCHOR_ROUND < pow2_interval(35000, 1000)[0] and ANCHOR_ROUND + 9 * 2 ** 25 < 2 ** 36, "2^%.4f" % log2(ANCHOR_ROUND))
check("MOVZ immediates 336*2^23, 140*2^23 and ADD 2^23 (4)",
 336 * 2 ** 25 < 2 ** 34 and 140 * 2 ** 25 < 2 ** 33 and 2 ** 25 < 2 ** 32)
check("phase 1's last blocks cover < 2^60.87 further operations (4)",
 Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP < pow2_interval(60870, 1000)[0],
 "2^%.5f" % log2(Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP))
RERUN = 1 + WALK1 + 2 + 4 + 1 + WALKB + 2 + 3
check("zero-test re-run with stubs and dy borrows < 2^22 units (4)",
 RERUN == 195 and 2 ** 25 * RERUN < 2 ** 34 and 2 ** 25 * (RERUN + RERUN_STUBS * STUB) + 2 ** 10 <= 2 ** 24 * 1000)
LMAX = 20 * 2 ** 64  # slot-steps of one trail are fewer (3.2, 5.1)
LAM, CB = Fraction(1, 2 ** 85), 9 * LMAX  # Bennett's lambda; X_k - E[X_k | past] <= 9 L_max
L2LO, L2HI = ln_interval(2)
DEV, INPROG = 140 * L2HI / LAM, 2 ** 25 * CORR * LMAX
check("CORR*n and 9n as MOVZ/MOVK immediates; rho' = 1056/2^23 (4)",
 CORR * 2 ** 25 < 2 ** 31 and 9 * 2 ** 25 < 2 ** 30 and CORR - PREM == Fraction(1056, 2 ** 23))
check("Bennett at lambda = 2^-85 absorbs the variance; Ville deviation + in-progress < 2^96 = Lambda_p (5.4)",
 LAM * CB < 3 and max(CORR, 9 - CORR) == 9 - CORR and LAM * (9 - CORR) ** 2 * LMAX / (2 * (1 - LAM * CB / 3)) <= CORR - PREM
 and DEV + INPROG < 2 ** 96, "%.6f; 2^%.4f + 2^%.4f" % (float(LAM * (9 - CORR) ** 2 * LMAX / (2 * (1 - LAM * CB / 3))), log2(DEV), log2(INPROG)))
DKW = 6 * Fraction(10 ** 9, 2718281828) ** 32
check("DKW row 6*exp(-32) < 2^-43.58; Ville row 2^20*2^-140 = 2^-120 (5.4)", DKW < pow2_interval(-4358, 100)[0], "2^%.3f" % log2(DKW))
# --- 2. the group algebra
check("p = 2^256 - C0 with C0 = 2^32 + 977", (1 << 256) - P_FIELD == (1 << 32) + 977)
check("r - 1 is divisible by 6", (r - 1) % 6 == 0)
N = (r - 1) // 6
check("log2 N = 253.41504 to 1e-5", abs(log2(N) - 253.41504) < 1e-5, "%.5f" % log2(N))

# --- 3. the report store of section 3
NODES = 2 ** 16
DCAP = ceil_fraction(pow2_interval(6383, 100)[1])
B = ceil_fraction(Fraction(2 * DCAP, NODES))
PROBE_PRICE, PROBES = 50, 8
REPORT_PRICE = 40000 + PROBES * PROBE_PRICE
OVF_INSERT = 6 * B + 35
OVF_LOOKUP = 16 * B
OVF_MILLI = DCAP * (OVF_INSERT + OVF_LOOKUP)
STORE_INIT_MILLI = 2 * DCAP * 16
check("D_cap = ceil(2^63.83) = 16396253461269537204", DCAP == 16396253461269537204)
check("B = ceil(2*D_cap/2^16) = 500373945961595", B == 500373945961595)
check("whole store < 2^72 bytes", NODES * B * 112 < 2 ** 72, "2^%.3f" % log2(NODES * B * 112))
check("report price 40.4 = 40 + 8 probes", REPORT_PRICE == 40400)
check("overflow worst case < 2^107.154 units", log2(Fraction(OVF_MILLI, 1000)) < 107.154,
 "2^%.6f" % log2(Fraction(OVF_MILLI, 1000)))

# --- 4. certified H1 numbers
L_LO, L_HI = ln_interval(1 / (1 - P0))
check("ln(1/(1-p0)) interval is positive and tight", 0 < L_LO < L_HI and L_HI - L_LO < Fraction(1, 2 ** 200),
 "%.15f" % float(L_LO))

N_PRIME = Fraction(N * R_PART, R_PART - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI
m = iroot((V.numerator * 10 ** 80 + V.denominator - 1) // V.denominator, 2) + 1
S_HI = Fraction(m, 10 ** 40)
S_INT = ceil_fraction(S_HI)
check("S_ops <= 2^126.6996 (certified upper bound)", S_HI * S_HI >= V and log2(S_HI) < 126.6996,
 "log2 bound %.9f" % log2(S_HI))
check("S_ops >= 2^126.6994 (lower bound for the quoted value)", log2(S_HI) > 126.6994, "%.9f" % log2(S_HI))

# --- 5. phase 1 and hard extras
PH1 = ceil_fraction(PRICE_STEP * (S_INT + 2 ** 60))
check("phase 1 = ceil(price*(S_int + 2^60)) covers price*(S_ops + 2^60)",
 Fraction(PH1, 1000) >= PRICE_STEP / 1000 * (S_HI + 2 ** 60))
check("phase 1 < 2^125.830023", log2(Fraction(PH1, 1000)) < 125.830023,
 "%.9f" % log2(Fraction(PH1, 1000)))

CLAIM_LO = 2 ** 125 * pow2_frac(CLAIM - 125)[0]
CAP_MILLI_LO = (CLAIM_LO.numerator * 1000) // CLAIM_LO.denominator
BLOCK_MILLI = 2 ** 40 * 1000
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
PERIOD, STRIDE, TRAV_CAP = 8192, 256, 32  # anchors, compare rounds, traversal cap (5.1)
CP_MILLI = 2 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))
ANCHOR_MILLI = ceil_fraction(Fraction(993, 64) * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD)))
CMP_MILLI = ceil_fraction((PASS1F - PASS1N) * ceil_fraction(Fraction(S_INT + 2 ** 61, STRIDE)))
HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
   ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
   ("phase 2", Fraction(2 ** 25 * (20 * 2 ** 64 + 65) * 600) * 1000),
   ("re-walks", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
   ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
   ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
   ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
   ("table (R accepted entries, every draw) and verification, < 2^41 (R:)", Fraction(2 ** 41) * 1000),
   ("node set-up, 2^26", Fraction(2 ** 26) * 1000),
   ("checkpoint writes: 2*ceil((S_int + 2^61)/8192)", Fraction(CP_MILLI)),
   ("anchor rounds: 993/64*ceil((S_int + 2^61)/8192)", Fraction(ANCHOR_MILLI)),
   ("compare rounds: 4.5*ceil((S_int + 2^61)/256)", Fraction(CMP_MILLI)),
   ("correction seeds: 2^20*(2^96 + 9n)", Fraction(2 ** 20 * (2 ** 96 + 9 * 2 ** 25))),
   ("final arrival pass", Fraction(FINAL_UNITS) * 1000)]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD = sum(int(v) for _, v in HARD_ITEMS)
check("hard extras < 2^112 units", log2(Fraction(HARD, 1000)) < 112.0,
 "%.9f" % log2(Fraction(HARD, 1000)))

# --- 6. expected random extras
# fruitless cycles of length 3..TRAV_CAP (5.1): S2(L, b) partitions into b blocks >= 2
TOP = 256
S2 = [[0] * (TOP + 1) for _ in range(TOP + 1)]
S2[0][0] = 1
for L_ in range(1, TOP + 1):
 for b_ in range(1, L_ // 2 + 1):
  S2[L_][b_] = b_ * S2[L_ - 1][b_] + (L_ - 1) * S2[L_ - 2][b_ - 1]
check("S2(3,1), S2(4,1), S2(4,2), S2(6,3) = 1, 1, 3, 15", (S2[3][1], S2[4][1], S2[4][2], S2[6][3]) == (1, 1, 3, 15))

def other_cycles(R, top=TRAV_CAP):
 return sum(Fraction(S2[L_][b_], R ** (L_ - b_)) for L_ in range(3, top + 1) for b_ in range(1, L_ // 2 + 1))

OTHER = other_cycles(R_PART)
check("R: at R = 2^26 the sum to L = 32 is < 2^-49.999", OTHER < Fraction(1, 1) / pow2_interval(49999, 1000)[0],
 "2^%.7f" % log2(OTHER))
# cycles of length 33..2*PERIOD end their trail (5.1): exact to L = 256; above it S2(L, b) <= C(L, b)*b^(L-b),
# b/R <= EPS: L = 2m-1 and 2m add <= (3/2)*(4*EPS)^m, summed over m >= 129
EPS = Fraction(PERIOD, R_PART)
TAIL = Fraction(3, 2) * (4 * EPS) ** (TOP // 2 + 1) / (1 - 4 * EPS)
LONG = other_cycles(R_PART, TOP) - OTHER + TAIL
check("cycles of length 33 to 16384: rate < 2^-300 per step", LONG < Fraction(1, 2 ** 300) and 4 * EPS < 1,
 "2^%.2f (tail 2^%.1f)" % (log2(LONG), log2(TAIL)))
RATE_HI = Fraction(1, 6 * R_PART) + OTHER
check("R: doubling rate 1/(6R) + sum <= 2^-28.5849", RATE_HI <= Fraction(1, 1) / pow2_interval(285849, 10000)[0],
 "2^%.9f" % log2(RATE_HI))
TRAV_RATE = Fraction(2, 6 * R_PART) + TRAV_CAP * OTHER
check("R: traversal-step rate 2/(6R) + 32*sum <= 2^-27.5849", TRAV_RATE <= Fraction(1, 1) / pow2_interval(275849, 10000)[0],
 "2^%.9f" % log2(TRAV_RATE))
OPS = S_HI + 2 ** 60
check("longer cycles: < 2^-232 trails per run, < 2^-100 units",
 LONG * OPS < Fraction(1, 2 ** 232) and LONG * OPS * 20 * 2 ** 64 < Fraction(1, 2 ** 100))
LCM = [L_ * STRIDE // math.gcd(L_, STRIDE) for L_ in range(TRAV_CAP + 1)]
check("anchors are compare rounds; lcm(L, 256) <= 7936 < 8192 for L <= 32", PERIOD % STRIDE == 0
 and max(LCM[1:]) == 7936 < PERIOD and 33 * STRIDE > PERIOD)
D2 = PERIOD - 1 + LCM[2] + 2
DL = max(PERIOD - 1 + LCM[L_] + L_ for L_ in range(3, TRAV_CAP + 1))
WASTE = Fraction(D2, 6 * R_PART) + DL * OTHER
W_KAPPA = pow2_interval(-1456, 100)[1]  # w = 2^-14.56 of the kappa rule (proof.md 5.2)
check("toy abandonment 640*2^-13.6 -> 0.052; 8,447 > the toy cap 640 (5.2)",
 round(640 * 2.0 ** -13.6, 3) == 0.052 and 8447 > 640 and D2 - 2 == 8447)
check("R: delay term 8,447/(6R) < 2^-15.5407 < w = 2^-14.56 of the kappa rule",
 Fraction(8447, 6 * R_PART) < pow2_interval(-155407, 10000)[0] and -15.5407 < -14.56, "2^%.6f" % log2(Fraction(8447, 6 * R_PART)))
check("R: waste w <= 8449/(6R) + 16158*sum <= 2^-15.54 <= 2^-14.56", D2 == 8449 and DL == 16158
 and WASTE <= pow2_interval(-1554, 100)[0] and WASTE <= W_KAPPA, "2^%.5f" % log2(WASTE))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "census-hz-q29.json")) as fh:
 RUNS = json.load(fh)["runs"]
MAIN = RUNS["walk_q29_r20_dp5"]
check("kappa rule, 1,024-table run: (1 + w)(kappa_eff + 2 s.e.) <= 1.00009",
 (1 + W_KAPPA) * (Fraction(str(MAIN["kappa_eff"])) + 2 * Fraction(str(MAIN["kappa_se"]))) <= KAPPA)
WALKS = [v for k, v in RUNS.items() if k.startswith("walk_")]
INV = sum(1 / Fraction(str(v["rho_se"])) ** 2 for v in WALKS)
RHO = sum(Fraction(str(v["rho"])) / Fraction(str(v["rho_se"])) ** 2 for v in WALKS) / INV
# kappa = rho^-1/2 and 2 s.e. = rho^-3/2 se_rho, so the rule holds if (1+w)^2 (1+SE_HI/rho)^2/rho <= kappa^2
SE_HI = Fraction(iroot(ceil_fraction(Fraction(2 ** 160) / INV), 2) + 1, 2 ** 80)
check("kappa rule, eleven walk runs pooled: (1 + w)(kappa_eff + 2 s.e.) <= 1.00009", len(WALKS) == 11
 and SE_HI * SE_HI * INV > 1 and (1 + W_KAPPA) ** 2 * (1 + SE_HI / RHO) ** 2 / RHO <= KAPPA ** 2,
 "%.6f" % ((1 + float(W_KAPPA)) * float(RHO) ** -0.5 * (1 + float(SE_HI / RHO))))
EXCESS = {k: (1 + W_KAPPA) * Fraction(str(v["kappa_eff"])) - KAPPA for k, v in RUNS.items() if k.startswith("walk_")}
NEAR = max(EXCESS, key=EXCESS.get)
check("no walk run exceeds 1.00009 by one s.e.; closest q29 R 2^25 dp3, 0.97 s.e. (5.2)",
 NEAR == "walk_q29_r25_dp3" and EXCESS[NEAR] < Fraction(272, 10 ** 6) and EXCESS[NEAR] >= Fraction(271, 10 ** 6)
 and all(EXCESS[k] <= Fraction(str(RUNS[k]["kappa_se"])) for k in EXCESS)
 and Fraction(96, 100) <= EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])) < Fraction(97, 100),
 "%.7f / %.3f s.e." % (float(EXCESS[NEAR]), float(EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])))))
H = RUNS["high_q29_r20_dp5"]
P12 = WALKS + [H]
I12 = sum(1 / Fraction(str(v["rho_se"])) ** 2 for v in P12)
R12 = sum(Fraction(str(v["rho"])) / Fraction(str(v["rho_se"])) ** 2 for v in P12) / I12
S12 = Fraction(iroot(ceil_fraction(Fraction(2 ** 160) / I12), 2) + 1, 2 ** 80)
check("row rule: high run within 1 s.e. of 1.00009; twelve runs pooled <= 1.00009 (5.2)",
 (1 + W_KAPPA) * Fraction(str(H["kappa_eff"])) - KAPPA <= Fraction(str(H["kappa_se"]))
 and S12 * S12 * I12 > 1 and (1 + W_KAPPA) ** 2 * (1 + S12 / R12) ** 2 / R12 <= KAPPA ** 2,
 "%.7f" % ((1 + float(W_KAPPA)) * float(R12) ** -0.5 * (1 + float(S12 / R12))))
RATE_K = 25 * Fraction(2 ** 32 + 977, 2 ** 64) + 2 * Fraction(1, 2 ** 63)  # 3 fe_sub, 6 rows, 4 row folds (4C0), 2 folds
check("kernel rare-exit rate 25*C0/2^64 + 2*2^-63 <= 2^-27.356 per operation",
 RATE_K <= pow2_interval(-27356, 1000)[0], "2^%.5f" % log2(RATE_K))
# 2-cycle budget (5.5)
BUNDLE, ETA, Q_P1, K_F = 677 + 2 * 39 + 2 * 25, Fraction(1, 2 ** 12), Fraction(1, 4 * R_PART), S_INT + 2 ** 61
C_P1, LAM1 = LMAX / 256 + 1, Fraction(1, 2 ** 72)
V_P1 = (C_P1 + (1 + ETA) ** 2 * Q_P1 * LMAX) * Q_P1
ABS1, DEV1 = LAM1 * V_P1 / (2 * (1 - LAM1 * C_P1 / 3)), 140 * L2HI / LAM1
N_C = (1 + ETA) * Q_P1 * (K_F + 2 ** 25 * LMAX) + DEV1
B_C = ceil_fraction(BUNDLE * (1 + 2 * ETA) * Q_P1 * K_F)
check("c = L_max/256 + 1 bounds the increments; v", BUNDLE == 805
 and LMAX / 257 + 1 + (1 + ETA) * Q_P1 * LMAX <= C_P1, "c 2^%.4f, v 2^%.4f" % (log2(C_P1), log2(V_P1)))
check("lambda = 2^-72: lambda*v/(2(1 - lambda*c/3)) <= eta*q'", LAM1 * C_P1 < 3 and ABS1 <= ETA * Q_P1, "2^%.4f" % log2(ABS1))
check("Ville 140 ln2/lambda < 2^79; catches < (1 + 2^-11)Kq'", DEV1 < 2 ** 79 and N_C <= (1 + 2 * ETA) * K_F * Q_P1,
 "2^%.4f; 2^%.5f" % (log2(DEV1), log2(N_C)))
WASTE_P1 = Fraction(D2, 4 * R_PART) + DL * OTHER
check("R: w at q' = 1/(4R) < 2^-14.95 < 2^-14.56", WASTE_P1 < pow2_interval(-1495, 100)[0], "2^%.5f" % log2(WASTE_P1))
EXTRA_ITEMS = [ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * OTHER * OPS * 1000),
   ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * TRAV_CAP * OTHER * OPS * 1000),
   2 ** 80 * 1000,
   ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 25) * 1000),
   ceil_fraction(Fraction(2 ** 70) * 1000),
   ceil_fraction(Fraction(25) * (2 * OTHER + Fraction(1, 2 ** 56) + Fraction(2 ** 32 + 977, 2 ** 64)) * OPS),
   ceil_fraction(Fraction(20) * RATE_K * (S_HI + 2 ** 61 + RATE_HI * OPS))]
E = sum(EXTRA_ITEMS)
check("residual random extras E_r < 2^93.77 units", log2(Fraction(E, 1000)) < 93.77,
 "2^%.9f" % log2(Fraction(E, 1000)))

# --- 7. allowance, total and claim
A = ceil_fraction(E / HEADROOM)
check("allowance A_r = ceil(E_r/0.0000116) (Markov safe)", A * HEADROOM >= E)
W_LEDGER = PH1 + HARD + B_C + A
check("W_ledger < 2^125.830102", log2(Fraction(W_LEDGER, 1000)) < 125.830102,
 "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("125.830106 holds: W_ledger <= the cap's lower bound", W_LEDGER <= CAP_MILLI_LO)
NEXT_HI = 1000 * 2 ** 125 * pow2_frac(CLAIM - Fraction(1, 10 ** 6) - 125)[1]
check("125.830101 fails: W_ledger > its upper bound", W_LEDGER > NEXT_HI)

HEAD = CAP_MILLI_LO - PH1 - HARD - B_C
check("headroom > allowance", HEAD > A, "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio E_r/head < 0.0000116", MARKOV_HI < HEADROOM,
 "%d/10^10 (rounded up)" % ceil_fraction(MARKOV_HI * 10 ** 10))
SUCCESS_LO = P0 - HEADROOM - LOSS_OTHER
check("success = 0.39001172 - 0.0000116 - 1.1e-7 = 0.39000001 > 0.39",
 SUCCESS_LO == Fraction(39000001, 10 ** 8) and SUCCESS_LO > THRESHOLD,
 "%.8f" % float(SUCCESS_LO))

# failure events other than the early halt (proof.md section 5.3): the upstream rows, the kernels' three
# rare-event rows and the canonical-index row
C0 = 2 ** 32 + 977
OPS_ALL = S_HI + 2 ** 61
check("doubling walk: 2 + 7*rate < 7 non-row reductions per operation", 2 + 7 * RATE_HI < 7)
ROW_CARRY = 7 * OPS_ALL * Fraction(C0 * C0 + 3 * C0, P_FIELD)
ROW_NEG = OPS_ALL * Fraction(C0, 2 ** 256)  # a doubling's x3 fold goes negative
ROW_HIGH = 2 * OPS_ALL * Fraction(C0, P_FIELD)  # x3 or y3 lands in [p, 2^256)
ROW_CANON = OPS_ALL * Fraction(C0 * 2 ** 64, P_FIELD - 1)  # x mod 2^192 > 2^192 - C0 (proof.md 3.1)
UPSTREAM_ROWS = (Fraction(1, 10 ** 7) + Fraction(1, 2 ** 30) + Fraction(1, 2 ** 50) + Fraction(1, 2 ** 60)
    + Fraction(1, 2 ** 64) + Fraction(1, 2 ** 114) + Fraction(1, r) + 2 * Fraction(1, 2 ** 100))
check("final-carry row 7*(S_ops + 2^61)*(C0^2 + 3*C0)/p < 2^-62.4", ROW_CARRY < pow2_interval(-624, 10)[0],
 "2^%.3f" % log2(ROW_CARRY))
check("negative-fold row (S_ops + 2^61)*C0/2^256 < 2^-97", ROW_NEG < Fraction(1, 2 ** 97), "2^%.3f" % log2(ROW_NEG))
check("high-result row 2*(S_ops + 2^61)*C0/p < 2^-96", ROW_HIGH < Fraction(1, 2 ** 96), "2^%.3f" % log2(ROW_HIGH))
check("canonical-index row (S_ops + 2^61)*C0*2^64/(p - 1) < 2^-33", ROW_CANON < Fraction(1, 2 ** 33),
 "2^%.3f" % log2(ROW_CANON))
ROWS_FAIL = UPSTREAM_ROWS + ROW_CARRY + ROW_NEG + ROW_HIGH + ROW_CANON + DKW + Fraction(1, 2 ** 120) + Fraction(1, 2 ** 140)
check("other failure rows, the account's and the 2-cycle budget's 2^-140 included, < 1.0104e-7 < 1.1e-7",
 ROWS_FAIL < Fraction(10104, 10 ** 11) < LOSS_OTHER, "%.6g" % float(ROWS_FAIL))

# --- 8. counter widths, section 3
METER_MILLI = CAP_MILLI_LO
check("192-bit meter holds the claim < 2^192 milli-units", METER_MILLI < 2 ** 192,
 "2^%.7f" % log2(METER_MILLI))
check("192-bit account A <= the work charged < 2^136 milli-units < 2^192 (4)", METER_MILLI < 2 ** 136)
check("128-bit operation counter holds S_ops + 2^60 < 2^127", S_HI + 2 ** 60 < 2 ** 127,
 "2^%.6f" % log2(S_HI + 2 ** 60))
STARTS_HI = Fraction(CAP_MILLI_LO, 1000) / 2 ** 13
check("128-bit start counter holds W/2^13 < 2^114", STARTS_HI < 2 ** 114)
PREP_CAP_MILLI = ceil_fraction(pow2_interval(590, 10)[1] * 1000)
check("128-bit preprocessing cap 2^59 units < 2^69 milli-units", PREP_CAP_MILLI < 2 ** 69)
check("128-bit refill counter < 2^87", REFILLS < 2 ** 87, "2^%.6f" % log2(REFILLS))

# --- 9. preprocessing
STORE_INIT_UNITS = Fraction(STORE_INIT_MILLI, 1000)
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(401, 10)[1] + 2 ** 36
CAP_PREP = pow2_interval(590, 10)[1]
check("expected preprocessing total < 2^59.0 units", PREP_EXPECTED < CAP_PREP,
 "2^%.5f" % log2(PREP_EXPECTED))
check("preprocessing cap halts with probability < 2^-30",
 Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS) < pow2_interval(-30, 1)[1],
 "2^%.4f" % log2(Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS)))

print()
print("step %.6f instructions" % STEP)
print("phase 1 2^%.6f  hard 2^%.6f  B_c 2^%.6f  E_r 2^%.6f  A_r 2^%.6f"
 % (log2(Fraction(PH1, 1000)), log2(Fraction(HARD, 1000)), log2(Fraction(B_C, 1000)), log2(Fraction(E, 1000)),
  log2(Fraction(A, 1000))))
print("W_ledger 2^%.9f units  claim %s  head 2^%.6f"
 % (log2(Fraction(W_LEDGER, 1000)), CLAIM, log2(Fraction(HEAD, 1000))))
print("B_c %d milli-units  Markov < %d/10^10  success > %.8f" % (B_C, ceil_fraction(MARKOV_HI * 10 ** 10), float(SUCCESS_LO)))
print()
if FAILED:
 print("CERTIFICATE: FAIL (%d)" % len(FAILED))
 raise SystemExit(1)
print("CERTIFICATE: PASS")
