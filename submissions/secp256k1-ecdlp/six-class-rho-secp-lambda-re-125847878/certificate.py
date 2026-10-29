"""Exact-integer certificate for W <= 2^125.847878 (proof.md 5.3): integers and Fractions only; ln(1/(1 - p0))
from the atanh series with its remainder, powers of two by certified integer roots.  Run python3 certificate.py
(exit 0 = all checks hold).  For kappa 1.0001, R 2^25, stride 128, p0 0.3900653, h 0.000065, step 554.0234375 +
26290/2^20: phase1 = ceil(price*(S_int + 2^60)), A = ceil(E/h), W = phase1 + hard + A, Markov
E/(cap - phase1 - hard) < h, and 1000*2^(claim - 10^-6) < W <= 1000*2^claim, claim = 125.847878.
"""
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
R_PART = 2 ** 25  # partition table rows
KAPPA = Fraction(10001, 10000)  # H1 slowdown allowance 1.0001 (proof.md section 5.2)
P0 = Fraction(3900653, 10000000)  # H1 floor p0 at S0 = S_ops / KAPPA
HEADROOM = Fraction(65, 1000000)  # 0.000065 Markov headroom h
CLAIM = Fraction(125847878, 10 ** 6)  # W = 2^125.847878
THRESHOLD = Fraction(39, 100)  # rules.md: success at least 0.39
LOSS_OTHER = Fraction(11, 10 ** 8)

# --- 1. counted kernels
MUL = (16 + 16) + 30 + 23  # products, six carry chains, reduction without a final fold
SQR = (10 + 10) + 8 + 6 + 7 + 23  # products, off-diagonal sum, doubling, diagonal chain, reduction
MUL_SUB = (16 + 16) + 30 + 24 + 5  # lambda*(x_W - x3) - y_W: fused reduction with one subtrahend
SQR_SUB2 = (10 + 10) + 8 + 6 + 7 + 24 + 10  # lambda^2 - x_W - x_T: fused reduction with two subtrahends
SUB, CANON, ADD_WEAK, ADD_WEAK2 = 9, 8, 9, 14
MUL_ADD = (16 + 16) + 30 + 23 + 5  # the loops' y3' = lambda*(x3 - x_W) + y_W: fused addition, unsigned fold
X3 = (10 + 10) + 8 + 6 + 7 + 23 + 10  # the loops' x3: subtrahends x_W and x_T + 2C0, both chains end ADC, unsigned fold
MUL_R, SUB_R, MUL_ADD_R, X3_R, SQR_R = MUL - 1, SUB - 2, MUL_ADD - 1, X3 - 1, SQR - 1  # kernels with rare exits (2.5)
RMUL, SHIFT = (16 + 16) + 34 + 6, 4  # a*b on b's shared rows B_i = b*2^(64i) mod p (2.5): products, sums, fold; one row
WALK1 = 2 * SUB_R + 8 + 4  # unfused pass 1 walk: dx, dy, 8 LDP, 4 STP (its table pointer comes with the group)
WALKB = 2 * RMUL + 3 * SHIFT + 4 + 2  # prod walk: acc's rows, e = acc*dy and acc*dx, 4 LDP, 2 STP
WALK2 = 2 * RMUL + 3 * SHIFT + X3_R + MUL_ADD_R + SUB_R + 10 + 4  # pass 2 walk: R's rows, lambda and R', x3, u, y3'; 10 LDP, 4 STP
PASS1 = Fraction(256 * WALK1 + 128 + 2, 256)  # 256 walks, 128 LDP of their pointers, SUBS + B.NE
PROD = Fraction(256 * WALKB + 2, 256)  # 256 walks, SUBS + B.NE
PASS2 = Fraction(256 * WALK2 + 128 + 2, 256)  # 256 walks, 128 LDP of their pointers, SUBS + B.NE
check("fe_mul costs 85 instructions", MUL == 85, str(MUL))
check("fe_sqr costs 64 instructions", SQR == 64, str(SQR))
check("kernels with rare exits: fe_mul 84, fe_sub 7, y3' 89, x3 73", (MUL_R, SUB_R, MUL_ADD_R, X3_R) == (84, 7, 89, 73))
check("fused fe_mul and fe_sqr cost 91 and 75", MUL_SUB == 91 and SQR_SUB2 == 75)
check("unfused pass 1: 6,786 per 256 walks = 26.5078125 per walk", 256 * WALK1 + 130 == 6786 and PASS1 == Fraction(3393, 128), str(PASS1))
check("prod: two products 72 on 3 rows of 4; 41,474 per 256 walks = 162.0078125 per walk",
 (RMUL, WALKB) == (72, 162) and 256 * WALKB + 2 == 41474 and PROD == Fraction(20737, 128), str(PROD))
check("pass 2: 2*72 + 12 + 73 + 89 + 7 + 14 = 339; 86,914 per 256 walks = 339.5078125 per walk",
 WALK2 == 339 and 256 * WALK2 + 130 == 86914 and PASS2 == Fraction(43457, 128), str(PASS2))
ARRIVE_F = 13 + 10 + 3 + 1 + 1 + 3 + 1
LOOPA = 2 + SUB_R + 2 + 2 + 2 + SUB_R + 2  # 2 LDP x(T), dx, 2 STP, 2 LDP y(W), 2 LDP y(T), dy, 2 STP
WALKF = ARRIVE_F + LOOPA
PASS1F = Fraction(256 * WALKF + 4, 256)
check("loop A walk: arrival 32 + 2 LDP x(T) + dx 7 + 2 STP + 4 LDP y(W), y(T) + dy 7 + 2 STP = 56", ARRIVE_F == 32 and LOOPA == 24 and WALKF == 56)
check("loop A (pass1f): 14,340 per 256 walks = 56.015625 per walk", 256 * WALKF + 4 == 14340 and PASS1F == Fraction(3585, 64),
 str(PASS1F))
ARRIVE_N = 3 + 11 + 6 + 1 + 1 + 4 + 1 + 1  # pass1p (2.5): 3 LDP, F1, select, TBNZ, one key test, 4 windows, B.HS, BFI
WALKN = ARRIVE_N + LOOPA  # a pair's second walk adds the STP of both addresses
PASS1N = Fraction(128 * (2 * WALKN + 1) + 2, 256)  # 128 pairs, SUBS + B.NE on x0
check("loop A (pass1p): arrival 28, walks 52 + 53 per pair, 13,442 per 256 walks = 52.5078125 per walk",
 ARRIVE_N == 28 and WALKN == 52 and 128 * (2 * WALKN + 1) + 2 == 13442 and PASS1N == Fraction(6721, 128), str(PASS1N))
check("compare rounds cost pass1f - pass1p = 3.5078125 per walk; anchor rounds arrive.s + pass1 - pass1f <= 1.508",
 PASS1F - PASS1N == Fraction(449, 128) and Fraction(1985, 64) + PASS1 - PASS1F == Fraction(193, 128))

ROWS = [("F1 truncated product (2 LDP, 4 UMULH, 3 ADD, 4 MADD; G1..G5 in registers)", 13),
  ("canonical index and slot: ADD, SUBS, CSEL, CSEL, ADDS, CSEL; CCMN, CCMN, CCMN, CCMP", 10),
  ("parity branch (TBNZ y0); K0 + 1 and K1 by SUBS, SBCS, or CMP y1 and ORR", 3),
  ("table address: BFI of K0 + 1 bits 0..25 into bits 9..34", 1),
  ("arrival traffic (LDP of y0, y1; one STP of a pair's two table addresses)", Fraction(3, 2)),
  ("checkpoint filter (one LDP of a pair's K_cp words, CCMP) and the rare branch (B.EQ)", Fraction(5, 2)),
  ("arrival loop of 32 pairs of walks (one TBZ on the pointer array's end bit per 64 walks)", Fraction(1, 64))]
ROWS_TOTAL = sum(v for _, v in ROWS)
check("unfused arrival rows (anchor rounds) sum to 31.015625", ROWS_TOTAL == Fraction(1985, 64), "%.6f" % float(ROWS_TOTAL))
# per-batch work (proof.md section 4)
FUSED_FIXED = 23 + 4 + 7  # loop A: prologue, the two group-count tests, epilogue and RET (executed: 34)
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
check("a64.py counts: inversion 24,322; fixed work pass1 9, prod 25, pass 2 26; group loops 6,786, 41,474, 86,914",
 INV == 24322 and FIX == dict(pass1=9, prod=25, pass2=26)
 and [K[k]["loop_iteration_max"] for k in ("pass1", "prod", "pass2")] == [6786, 41474, 86914])
CALLS = 64 + 8  # the passes' call and return; prod's call (MOV x0, MOV x3, BL: at most 8)
BATCH_FIXED = INV + FUSED_FIXED + FIX["prod"] + FIX["pass2"] + CALLS + 12 + 3 + 7 + 3 + 2 + 4
TAILS_ROUND = 255 * (WALKN + 3 - PASS1N) + 255 * (WALKB + 2 - PROD) + 255 * (1 + WALK2 + 2 - PASS2)  # at most 255 tail walks per loop
ROUND_FIXED = ceil_fraction(BATCH_FIXED + TAILS_ROUND)
check("per-round work 24,510 + 1,779.0234375 -> 26,290, per round of 2^20 walks",
 BATCH_FIXED == 24510 and TAILS_ROUND == Fraction(227715, 128) and ROUND_FIXED == 26290
 and 34 + 255 * (WALKN + 3 - PASS1N) <= FUSED_FIXED + 255 * (WALKF + 4 - PASS1F),  # pass1p's fixed 34 and tail walk 55 are below pass1f's
 str(ROUND_FIXED))
SHARE = Fraction(ROUND_FIXED, 2 ** 20)

STEP = PASS1N + PROD + PASS2 + SHARE
DBL_ARITH = 4 * MUL_R + MUL_SUB + SQR_R + SQR_SUB2 + ADD_WEAK + 2 * ADD_WEAK2 + SUB_R + 30
DBL = STEP + DBL_ARITH + 34 + 3 + SHARE
TRAV = STEP + 24 + 2 + 7 + 6
PRICE_STEP = Fraction(70915, 128) + SHARE
PRICE_TRAV, PRICE_DBL = PRICE_STEP + 39, PRICE_STEP + 677
check("step is 52.5078125 + 162.0078125 + 339.5078125 + share = 554.0234375 + 26,290/2^20", STEP == PRICE_STEP, "%.7f" % float(STEP))
check("doubling bound: displaced step + 639 + 34 + 3 + share = step + 676 + share",
 DBL_ARITH == 639 and DBL == STEP + 676 + SHARE)
check("traversal bound is step + 24 + 2 + 7 + 6", TRAV == PRICE_STEP + 39)
check("traversal price covers its bound", PRICE_TRAV >= TRAV)
check("doubling price step + 677 covers its bound", PRICE_DBL >= DBL)

# --- 2b. meter, stubs, reserves: every such bound proof.md states (2.5, 3.2, 4); 1 milli-unit = 1 instruction
STUB, RARE, RESERVE, WALK_STUBS, RERUN_STUBS = 20, 25, 336, 15, 7  # milli-units; 15 = 2 loop-A + 5 prod + 8 pass-2 stubs
def stub_macro(fname, name):
 """Instructions on a stub macro's longest path (its whole body, RCHG expanded)."""
 t = open(os.path.join(KD, fname), encoding="utf-8").read()
 body = lambda m: [x.strip() for x in t[t.index(".macro " + m + "\n") if m == "RCHG" else t.index(".macro " + m + " "):]
                   .split(".endm")[0].splitlines()[1:] if x.strip() and not x.strip().endswith(":")]
 return sum(len(body("RCHG")) if x == "RCHG" else 1 for x in body(name))


F, LB, LC, LT = stub_macro("secp_fuse.s", "SS"), stub_macro("secp_loops.s", "SB"), stub_macro("secp_loops.s", "SC"), stub_macro("secp_loops.s", "ST")
LABELS = {f: sum(1 for k in a64.a64_load(open(os.path.join(KD, f)).read())[1] if k.startswith(".Lr")) for f in ("secp_loops.s", "secp_fuse.s")}
check("stubs: loops 3,855 = 257*(2 + 5 + 8), loop A 1,028 = 2*257*2", LABELS == {"secp_loops.s": 3855, "secp_fuse.s": 1028}
 and 257 * (2 + 5 + 8) == 3855 and 2 * 257 * 2 == 1028)
check("shipped stubs loop A 7, loops SB 9, SC 8, row ST 13: all <= 13 instructions, 20 milli-units (2.5)",
 (F, LB, LC, LT) == (7, 9, 8, 13) and max(F, LB, LC, LT) <= 13 <= STUB, "%d %d %d %d" % (F, LB, LC, LT))
check("rare path 5 + 5 + 5 + 4 + 3 = 22 instructions <= its 25 milli-units (3.2)",
 5 + 5 + 5 + 4 + 3 == 22 <= RARE)
check("a walk's 15 stubs = 300; 15*20 + 25 = 325 <= the 336 reserve (2.5, 4)", WALK_STUBS * STUB == 300
 and WALK_STUBS * STUB + RARE == 325 <= RESERVE)
FINAL_UNITS, FINAL_CHARGE = 2 ** 15, 2 ** 15 * 1000 + 2 ** 20 * 25  # 2^15 units plus 25 milli-units per walk
check("final arrival pass (2.5): 2^15 units + 25 per walk cover 2^20 x 31.015625 + one entry each, under 2^15.85 units",
 FINAL_CHARGE >= 2 ** 20 * ROWS_TOTAL + 2 ** 20 * 25 and log2(Fraction(FINAL_CHARGE, 1000)) < 15.85
 and log2(Fraction(FINAL_CHARGE, 1000)) < 70,
 "2^%.5f" % log2(Fraction(FINAL_CHARGE, 1000)))
CMP_ROUND = (PASS1F + PROD + PASS2 + RESERVE + 2) * 2 ** 20 + ROUND_FIXED
check("a compare round's charge (557.53125 + 336 + 2)*2^20 + 26,290 on a 2^40-unit block is < 2^30 (4)",
 PASS1F + PROD + PASS2 == Fraction(17841, 32) and CMP_ROUND < 2 ** 30, "2^%.4f" % log2(CMP_ROUND))
check("phase 1's last blocks cover < 2^60.86 further operations (4)",
 Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP < pow2_interval(60860, 1000)[0],
 "2^%.5f" % log2(Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP))
RERUN = 1 + WALK1 + 2 + 4 + WALKB + 2  # tail walks of pass1 (with the per-walk dx = 0 test) and prod
check("the zero test's re-run 2^20*(197 + 7*20) + fixed < 2^19 units: its own action, stubs included (4)",
 RERUN == 197 and 2 ** 20 * RERUN < 2 ** 28 and 2 ** 20 * (RERUN + RERUN_STUBS * STUB) + 2 ** 10 <= 2 ** 19 * 1000)
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
L_LO, L_HI = ln_interval(1 / (1 - P0))  # ln(1/(1 - p0))
check("ln(1/(1-p0)) interval is positive and tight", 0 < L_LO < L_HI and L_HI - L_LO < Fraction(1, 2 ** 200),
 "%.15f" % float(L_LO))

N_PRIME = Fraction(N * R_PART, R_PART - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI  # upper bound on S_ops^2
m = iroot((V.numerator * 10 ** 80 + V.denominator - 1) // V.denominator, 2) + 1
S_HI = Fraction(m, 10 ** 40)
S_INT = ceil_fraction(S_HI)
check("S_ops <= 2^126.6996 (certified upper bound)", S_HI * S_HI >= V and log2(S_HI) < 126.6996,
 "log2 bound %.9f" % log2(S_HI))
check("S_ops >= 2^126.6995 (lower bound for the quoted value)", log2(S_HI) > 126.6995, "%.9f" % log2(S_HI))

# --- 5. phase 1 and hard extras
PH1 = ceil_fraction(PRICE_STEP * (S_INT + 2 ** 60))  # milli-units
check("phase 1 = ceil(price*(S_int + 2^60)) covers price*(S_ops + 2^60)",
 Fraction(PH1, 1000) >= PRICE_STEP / 1000 * (S_HI + 2 ** 60))
check("phase 1 < 2^125.84769", log2(Fraction(PH1, 1000)) < 125.84769,
 "%.9f" % log2(Fraction(PH1, 1000)))

CLAIM_LO, CLAIM_HI = (2 ** 125 * v for v in pow2_frac(CLAIM - 125))
P4_LO, P4_HI = pow2_interval(8928, 10000)
check("2^f by series agrees with the integer-root interval at f = 0.8928", pow2_frac(Fraction(8928, 10000))[0] <= P4_HI
 and P4_LO <= pow2_frac(Fraction(8928, 10000))[1] and CLAIM_HI - CLAIM_LO < Fraction(1, 2 ** 60))
CAP_MILLI_LO = (CLAIM_LO.numerator * 1000) // CLAIM_LO.denominator
BLOCK_MILLI = 2 ** 40 * 1000
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
PERIOD, STRIDE, TRAV_CAP = 4096, 128, 32  # anchors, compare rounds, traversal cap (5.1)
CP_MILLI = 2 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))  # <= 1 write (STP + STR) per anchor
ANCHOR_MILLI = 14 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))  # anchor work per walk: full K 7, counter 5, unfused 1.5
CMP_RATE = (PASS1F - PASS1N) + Fraction(380, 2 ** 20)
CMP_MILLI = ceil_fraction(CMP_RATE * ceil_fraction(Fraction(S_INT + 2 ** 61, STRIDE)))  # compare rounds run pass1f and compare-round fixed delta
HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
   ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
   ("phase 2", Fraction(2 ** 20 * (20 * 2 ** 64 + 65) * 600) * 1000),
   ("re-walks", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
   ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
   ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
   ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
   ("table (2^25 entries) and verification, < 2^39", Fraction(2 ** 39) * 1000),
   ("node set-up, 2^26", Fraction(2 ** 26) * 1000),
   ("checkpoint writes: 2*ceil((S_int + 2^61)/4096)", Fraction(CP_MILLI)),
   ("anchor rounds: 14*ceil((S_int + 2^61)/4096)", Fraction(ANCHOR_MILLI)),
   ("compare rounds: (449/128 + 380/2^20)*ceil((S_int + 2^61)/128)", Fraction(CMP_MILLI))]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD = sum(int(v) for _, v in HARD_ITEMS)
check("hard extras < 2^111.7984 units", log2(Fraction(HARD, 1000)) < 111.7984,
 "%.9f" % log2(Fraction(HARD, 1000)))
check("the compare-round row is the largest hard row", CMP_MILLI == max(int(v) for _, v in HARD_ITEMS))

# --- 6. expected random extras
# fruitless cycles of length 3..TRAV_CAP (proof.md section 5.1: S2(L, b) counts partitions into b blocks >= 2)
TOP = 256
S2 = [[0] * (TOP + 1) for _ in range(TOP + 1)]
S2[0][0] = 1
for L_ in range(1, TOP + 1):
 for b_ in range(1, L_ // 2 + 1):
  S2[L_][b_] = b_ * S2[L_ - 1][b_] + (L_ - 1) * S2[L_ - 2][b_ - 1]
check("S2(3,1), S2(4,1), S2(4,2), S2(6,3) = 1, 1, 3, 15", (S2[3][1], S2[4][1], S2[4][2], S2[6][3]) == (1, 1, 3, 15))


def other_cycles(R, top=TRAV_CAP):
 max_exp = max(L_ - b_ for L_ in range(3, top + 1) for b_ in range(1, L_ // 2 + 1))
 num = sum(S2[L_][b_] * (R ** (max_exp - (L_ - b_))) for L_ in range(3, top + 1) for b_ in range(1, L_ // 2 + 1))
 return Fraction(num, R ** max_exp)


check("at R = 2048 the sum to L = 64 is 2^-19.9954 (the public package's value)",
 abs(log2(other_cycles(2048, 64)) + 19.9954) < 5e-5, "2^%.6f" % log2(other_cycles(2048, 64)))
OTHER = other_cycles(R_PART)
check("at R = 2^25 the sum to L = 32 is < 2^-47.999", OTHER < Fraction(1, 1) / pow2_interval(47999, 1000)[0],
 "2^%.7f" % log2(OTHER))
# cycles longer than TRAV_CAP end their trail (proof.md 5.1): exact to L = 256, above it S2(L, b) <= (L/2)^L
TAIL = sum(Fraction((L_ // 2) * L_ ** L_, 2 ** L_ * R_PART ** ((L_ + 1) // 2)) for L_ in range(TOP + 1, 2 * PERIOD + 1))
LONG = other_cycles(R_PART, TOP) - OTHER + TAIL
check("cycles of length 33 to 8192: rate < 2^-300 per step", LONG < Fraction(1, 2 ** 300), "2^%.2f" % log2(LONG))
RATE_HI = Fraction(1, 6 * R_PART) + OTHER
check("doubling rate 1/(6R) + sum <= 2^-27.5849", RATE_HI <= Fraction(1, 1) / pow2_interval(275849, 10000)[0],
 "2^%.9f" % log2(RATE_HI))
TRAV_RATE = Fraction(2, 6 * R_PART) + TRAV_CAP * OTHER  # traversal steps per operation
check("traversal-step rate 2/(6R) + 32*sum <= 2^-26.5849", TRAV_RATE <= Fraction(1, 1) / pow2_interval(265849, 10000)[0],
 "2^%.9f" % log2(TRAV_RATE))
OPS = S_HI + 2 ** 60
check("trails ended by a longer cycle: LONG*ops < 2^-232 per run, LONG*ops*L_max < 2^-100 units",
 LONG * OPS < Fraction(1, 2 ** 232) and LONG * OPS * 20 * 2 ** 64 < Fraction(1, 2 ** 100))
# the strided compare (proof.md 5.1): the coverage, delay and traversal bounds below
LCM = [L_ * STRIDE // math.gcd(L_, STRIDE) for L_ in range(TRAV_CAP + 1)]
check("anchors are compare rounds and lcm(L, 128) <= 3968 < 4096 for L = 1..32, not for L = 33", PERIOD % STRIDE == 0
 and max(LCM[1:]) == 3968 < PERIOD and 33 * STRIDE > PERIOD)
D2 = PERIOD - 1 + LCM[2] + 2
DL = max(PERIOD - 1 + LCM[L_] + L_ for L_ in range(3, TRAV_CAP + 1))
WASTE = Fraction(D2, 6 * R_PART) + DL * OTHER
W_KAPPA = pow2_interval(-1456, 100)[1]  # w = 2^-14.56 of the kappa rule (proof.md 5.2)
# proof.md 5.2's toy arithmetic, both products stated there and the delay term they identify
check("toy abandonment 640*2^-13.6 = 0.0515 rounds to the 0.052 proof.md prints, and 4,223 > the toy cap 640",
 round(640 * 2.0 ** -13.6, 3) == 0.052 and 4223 > 640 and D2 - 2 == 4223)
check("the delay term 4,223*2^-27.585 = 2^-15.5409 lies below the w = 2^-14.56 the kappa rule allows",
 round(log2(4223 * 2.0 ** -27.585), 4) == -15.5409 and log2(4223 * 2.0 ** -27.585) < -14.56)
check("waste w <= 4225/(6R) + 8094*sum <= 2^-15.54 <= 2^-14.56", D2 == 4225 and DL == 8094
 and WASTE <= pow2_interval(-1554, 100)[0] and WASTE <= W_KAPPA, "2^%.5f" % log2(WASTE))
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "census-hz-q29.json")) as fh:
 RUNS = json.load(fh)["runs"]
MAIN = RUNS["walk_q29_r20_dp5"]
check("kappa rule, 1,024-table run: (1 + w)(kappa_eff + 2 s.e.) <= 1.0001",
 (1 + W_KAPPA) * (Fraction(str(MAIN["kappa_eff"])) + 2 * Fraction(str(MAIN["kappa_se"]))) <= KAPPA)
WALKS = [v for k, v in RUNS.items() if k.startswith("walk_")]
INV = sum(1 / Fraction(str(v["rho_se"])) ** 2 for v in WALKS)
RHO = sum(Fraction(str(v["rho"])) / Fraction(str(v["rho_se"])) ** 2 for v in WALKS) / INV
# kappa = rho^-1/2 and 2 s.e. = rho^-3/2 se_rho, so the rule holds if (1+w)^2 (1+SE_HI/rho)^2/rho <= kappa^2
SE_HI = Fraction(iroot(ceil_fraction(Fraction(2 ** 160) / INV), 2) + 1, 2 ** 80)
check("kappa rule, ten walk runs pooled: (1 + w)(kappa_eff + 2 s.e.) <= 1.0001", len(WALKS) == 10
 and SE_HI * SE_HI * INV > 1 and (1 + W_KAPPA) ** 2 * (1 + SE_HI / RHO) ** 2 / RHO <= KAPPA ** 2,
 "%.6f" % ((1 + float(W_KAPPA)) * float(RHO) ** -0.5 * (1 + float(SE_HI / RHO))))
# proof.md 5.2: no walk run exceeds kappa by one s.e.; the closest is q29 dp3
EXCESS = {k: (1 + W_KAPPA) * Fraction(str(v["kappa_eff"])) - KAPPA for k, v in RUNS.items() if k.startswith("walk_")}
NEAR = max(EXCESS, key=EXCESS.get)
check("no walk run exceeds 1.0001 by one s.e. (5.2); closest q29 R 2^25 dp3, 0.00026 = 0.93 s.e.",
 NEAR == "walk_q29_r25_dp3" and EXCESS[NEAR] < Fraction(27, 100000) and EXCESS[NEAR] >= Fraction(26, 100000)
 and all(EXCESS[k] <= Fraction(str(RUNS[k]["kappa_se"])) for k in EXCESS)
 and Fraction(93, 100) <= EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])) < Fraction(94, 100),
 "%.7f / %.3f s.e." % (float(EXCESS[NEAR]), float(EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])))))
# rare exits per operation (proof.md 5.3)
RATE_K = 25 * Fraction(2 ** 32 + 977, 2 ** 64) + 2 * Fraction(1, 2 ** 63)  # 3 fe_sub, 6 row shifts, 4 row folds (4C0 each), 2 folds
check("kernel rare-exit rate 25*C0/2^64 + 2*2^-63 <= 2^-27.356 per operation",
 RATE_K <= pow2_interval(-27356, 1000)[0], "2^%.5f" % log2(RATE_K))
HITS = 2 * RATE_HI  # true filter hits: the doubling's and the closing arrival's (proof.md 5.3)
check("true filter hits 2*(1/(6R) + sum) <= 2^-26.5849 per operation", HITS <= Fraction(1, 1) / pow2_interval(265849, 10000)[0],
 "2^%.9f" % log2(HITS))
EXTRA_ITEMS = [ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * RATE_HI * OPS * 1000),
   ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * TRAV_RATE * OPS * 1000),
   ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 20) * 1000),
   ceil_fraction(Fraction(2 ** 70) * 1000),
   ceil_fraction(Fraction(25) * (HITS + Fraction(1, 2 ** 56) + Fraction(2 ** 32 + 977, 2 ** 64)) * OPS),  # rare-path entries at 0.025 units
   ceil_fraction(Fraction(20) * RATE_K * (S_HI + 2 ** 61 + RATE_HI * OPS))]  # rare exits at 0.02 units; + a doubling's own
E = sum(EXTRA_ITEMS)
check("expected random extras < 2^98.8452 units", log2(Fraction(E, 1000)) < 98.8452,
 "2^%.9f" % log2(Fraction(E, 1000)))

# --- 7. allowance, total and claim
A = ceil_fraction(E / HEADROOM)
check("allowance A = ceil(E/0.000065) (Markov safe)", A * 65 >= E * 1000000)
W_LEDGER = PH1 + HARD + A
check("W_ledger < 2^125.847878", log2(Fraction(W_LEDGER, 1000)) < 125.847878,
 "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("the claim is 125847878/10^6", CLAIM == Fraction(125847878, 10 ** 6))
check("125.847878 holds: W_ledger <= floor(1000 * 2^125 * (lower bound of 2^0.847878))", W_LEDGER <= CAP_MILLI_LO)
NEXT_HI = 1000 * 2 ** 125 * pow2_frac(CLAIM - Fraction(1, 10 ** 6) - 125)[1]
check("125.847877 fails: W_ledger > 1000 * 2^125 * (upper bound of 2^0.847877)", W_LEDGER > NEXT_HI)

HEAD = CAP_MILLI_LO - PH1 - HARD
check("headroom > allowance A", HEAD > A, "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio < 0.000065", MARKOV_HI < HEADROOM and MARKOV_HI * 1000000 < 65,
 "%d/10^10 (rounded up)" % ceil_fraction(MARKOV_HI * 10 ** 10))
SUCCESS_LO = P0 - HEADROOM - LOSS_OTHER
check("success = 0.3900653 - 0.000065 - 1.1e-7 = 0.39000019 > 0.39",
 SUCCESS_LO == Fraction(3900001900, 10 ** 10) and SUCCESS_LO > THRESHOLD,
 "%.7f" % float(SUCCESS_LO))
check("success floor >= the 0.39 of rules.md", SUCCESS_LO > THRESHOLD, "%.7f" % float(SUCCESS_LO))

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
check("all failure rows other than the early halt < 1.1e-7",
 UPSTREAM_ROWS + ROW_CARRY + ROW_NEG + ROW_HIGH + ROW_CANON < LOSS_OTHER,
 "%.4g" % float(UPSTREAM_ROWS + ROW_CARRY + ROW_NEG + ROW_HIGH + ROW_CANON))

# --- 8. counter widths, section 3
METER_MILLI = CAP_MILLI_LO
check("192-bit meter holds the claim < 2^192 milli-units", METER_MILLI < 2 ** 192,
 "2^%.7f" % log2(METER_MILLI))
check("128-bit operation counter holds S_ops + 2^60 < 2^127", S_HI + 2 ** 60 < 2 ** 127,
 "2^%.6f" % log2(S_HI + 2 ** 60))
STARTS_HI = Fraction(CAP_MILLI_LO, 1000) / 2 ** 13
check("128-bit start counter holds W/2^13 < 2^114", STARTS_HI < 2 ** 114)
PREP_CAP_MILLI = ceil_fraction(pow2_interval(590, 10)[1] * 1000)
check("128-bit preprocessing cap 2^59 units < 2^69 milli-units", PREP_CAP_MILLI < 2 ** 69)
check("128-bit refill counter < 2^87", REFILLS < 2 ** 87, "2^%.6f" % log2(REFILLS))

# --- 9. preprocessing
STORE_INIT_UNITS = Fraction(STORE_INIT_MILLI, 1000)
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(381, 10)[1] + 2 ** 33  # table 2^38.1, first 2^20 starts 2^33
CAP_PREP = pow2_interval(590, 10)[1]
check("expected preprocessing total < 2^59.0 units", PREP_EXPECTED < CAP_PREP,
 "2^%.5f" % log2(PREP_EXPECTED))
check("preprocessing cap halts with probability < 2^-30",
 Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS) < pow2_interval(-30, 1)[1],
 "2^%.4f" % log2(Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS)))


print()
print("step      %.6f instructions = %.9f units" % (STEP, STEP / 1000.0))
print("doubling  %.6f instructions (with the displaced step) = %.9f units" % (DBL, DBL / 1000.0))
print("traversal %.6f instructions = %.9f units" % (TRAV, TRAV / 1000.0))
print("phase 1 2^%.6f  hard 2^%.6f  extras 2^%.6f  allowance 2^%.6f"
 % (log2(Fraction(PH1, 1000)), log2(Fraction(HARD, 1000)), log2(Fraction(E, 1000)),
  log2(Fraction(A, 1000))))
print("W_ledger 2^%.9f units  claim %s  head 2^%.6f"
 % (log2(Fraction(W_LEDGER, 1000)), CLAIM, log2(Fraction(HEAD, 1000))))
print("Markov < %d/10^10  success > %.7f" % (ceil_fraction(MARKOV_HI * 10 ** 10), float(SUCCESS_LO)))
print("milli-units: S_int %d\n phase 1 %d\n hard %d\n extras %d\n A %d\n W_ledger %d\n cap %d"
 % (S_INT, PH1, HARD, E, A, W_LEDGER, CAP_MILLI_LO))
print()
if FAILED:
 print("CERTIFICATE: FAIL (%d)" % len(FAILED))
 raise SystemExit(1)
print("CERTIFICATE: PASS")
