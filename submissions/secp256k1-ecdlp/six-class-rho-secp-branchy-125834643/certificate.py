"""Exact-integer certificate for proof.md 5.3 (integers and Fractions only; ln by the atanh series, powers of two
by certified integer roots): W = phase1 + hard + ceil(E/h) <= 1000*2^claim, the next lower six-decimal claim fails,
and the Markov ratio is below h. Run python3 certificate.py; exit 0 = all checks hold."""
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
R_PART = 2 ** 26  # partition table rows (the R-dependent bounds below are marked "R:")
KAPPA = Fraction(10001, 10000)  # H1 slowdown allowance 1.0001 (proof.md section 5.2)
P0 = Fraction(3900482, 10000000)  # H1 floor p0 at S0 = S_ops / KAPPA
HEADROOM = Fraction(48, 1000000)  # 0.000048 Markov headroom h
CLAIM = Fraction(125834643, 10 ** 6)  # W = 2^125.834643
THRESHOLD = Fraction(39, 100)  # rules.md: success at least 0.39
LOSS_OTHER = Fraction(11, 10 ** 8)

# --- 1. counted kernels
MUL = (16 + 16) + 30 + 23  # products, six carry chains, reduction without a final fold
SQR = (10 + 10) + 8 + 6 + 7 + 23  # products, off-diagonal sum, doubling, diagonal chain, reduction
MUL_SUB = (16 + 16) + 30 + 24 + 5  # lambda*(x_W - x3) - y_W: fused reduction with one subtrahend
SQR_SUB2 = (10 + 10) + 8 + 6 + 7 + 24 + 10  # lambda^2 - x_W - x_T: fused reduction with two subtrahends
SUB, CANON, ADD_WEAK, ADD_WEAK2 = 9, 8, 9, 14
MUL_ADD = (16 + 16) + 30 + 23 + 5  # the loops' y3' = lambda*(x3 - x_W) + y_W: fused addition, unsigned fold
X3 = (10 + 10) + 8 + 6 + 7 + 23 + 10  # the loops' x3: x_W and x_T + 2C0 off, chains end ADC, unsigned fold
MUL_R, SUB_R, MUL_ADD_R, X3_R, SQR_R = MUL - 1, SUB - 2, MUL_ADD - 1, X3 - 1, SQR - 1  # kernels with rare exits (2.5)
RMUL, SHIFT = (16 + 16) + 34 + 6, 4  # a*b on b's shared rows B_i = b*2^(64i) mod p (2.5); one row
SUB_B = SUB - 4  # branchy fe_sub (2.4) on its main path: SUBS, 3 SBCS, B.HS; a borrow adds ADD, SUBS, B.LO
WALK1 = SUB_R + 4 + 2  # unfused pass 1 walk: dx, 4 LDP, 2 STP (its table pointer comes with the group)
WALKB = SUB_B + 2 * RMUL + 3 * SHIFT + 6 + 2  # prod walk: dy, acc's rows, e = acc*dy and acc*dx, 6 LDP, 2 STP
WALK2 = 2 * RMUL + 3 * SHIFT + X3_R + MUL_ADD_R + SUB_B + 10 + 4  # pass 2: R's rows, lambda, R', x3, u, y3'; 10 LDP, 4 STP
PASS1 = Fraction(128 * WALK1 + 64 + 1, 128)  # each loop: 128 walks, 64 LDP of their pointers, TBZ or CBNZ
PROD = Fraction(128 * WALKB + 64 + 1, 128)
PASS2 = Fraction(128 * WALK2 + 64 + 1, 128)
check("fe_mul costs 85 instructions", MUL == 85, str(MUL))
check("fe_sqr costs 64 instructions", SQR == 64, str(SQR))
check("kernels with rare exits: fe_mul 84, fe_sub 7, y3' 89, x3 73", (MUL_R, SUB_R, MUL_ADD_R, X3_R) == (84, 7, 89, 73))
check("fused fe_mul and fe_sqr cost 91 and 75", MUL_SUB == 91 and SQR_SUB2 == 75)
check("unfused pass 1: 1,729 per 128 walks = 13.5078125 per walk", 128 * WALK1 + 65 == 1729 and PASS1 == Fraction(1729, 128), str(PASS1))
check("prod: dy 5, two products 72 on 3 rows; 21,697 per 128 = 169.5078125",
 (SUB_B, RMUL, WALKB) == (5, 72, 169) and 128 * WALKB + 65 == 21697 and PROD == Fraction(21697, 128), str(PROD))
check("pass 2: 2*72 + 12 + 73 + 89 + 5 + 14 = 337; 43,201 per 128 = 337.5078125",
 WALK2 == 337 and 128 * WALK2 + 65 == 43201 and PASS2 == Fraction(43201, 128), str(PASS2))
ARRIVE_F = 13 + 10 + 3 + 1 + 1 + 3 + 1
LOOPA = 2 + SUB_B + 2  # 2 LDP x(T), dx, 2 STP
WALKF = ARRIVE_F + LOOPA
PASS1F = Fraction(128 * WALKF + 1, 128)  # 128 walks, TBZ x2, #26
check("loop A walk: arrival 32 + 2 LDP x(T) + dx 5 + 2 STP = 41", ARRIVE_F == 32 and LOOPA == 9 and WALKF == 41)
check("loop A (pass1f): 5,249 per 128 walks = 41.0078125 per walk", 128 * WALKF + 1 == 5249 and PASS1F == Fraction(5249, 128),
 str(PASS1F))
ARRIVE_N = 3 + 11 + 6 + 1 + 1 + 4 + 1 + 1  # pass1p (2.5): 3 LDP, F1, select, TBNZ, one key test, 4 windows, B.HS, BFI
WALKN = ARRIVE_N + LOOPA  # a pair's second walk adds the STP of both addresses
PASS1N = Fraction(64 * (2 * WALKN + 1) + 1, 128)  # 64 pairs, TBZ x2, #26
check("loop A (pass1p): arrival 28, walks 37 + 38 per pair, 4,801 per 128 = 37.5078125",
 ARRIVE_N == 28 and WALKN == 37 and 64 * (2 * WALKN + 1) + 1 == 4801 and PASS1N == Fraction(4801, 128), str(PASS1N))
check("compare rounds cost pass1f - pass1p = 3.5 per walk; anchor rounds arrive.s + pass1 - pass1f = 3.515625",
 PASS1F - PASS1N == Fraction(7, 2) and Fraction(1985, 64) + PASS1 - PASS1F == Fraction(225, 64))

ROWS = [13, 10, 3, 1, Fraction(3, 2), Fraction(5, 2), Fraction(1, 64)]  # the arrival rows of proof.md section 3
ROWS_TOTAL = sum(ROWS)
check("unfused arrival rows (anchor rounds) sum to 31.015625", ROWS_TOTAL == Fraction(1985, 64), "%.6f" % float(ROWS_TOTAL))
# per-batch work (proof.md section 4)
FUSED_FIXED = 22 + 3 + 8  # loop A: prologue, tail-count and group tests, MOV x0, epilogue, RET (executed: 33)
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
GRP = [K[k]["loops"][-1] for k in ("pass1", "prod", "pass2")]  # the group loops; main = every correction skip taken
check("a64.py: inversion 24,322; fixed pass1 8, prod 28, pass 2 28; group main paths 1,729, 21,697, 43,201",
 INV == 24322 and FIX == dict(pass1=8, prod=28, pass2=28) and [l["main"] for l in GRP] == [1729, 21697, 43201]
 and [l["iteration"] - l["main"] for l in GRP] == [0, 3 * 128, 3 * 128] and [l["skips"] for l in GRP] == [0, 128, 128]
 and all(l["inner_branches"] == 0 for k in ("pass1", "prod", "pass2") for l in K[k]["loops"]))
CALLS = 64 + 8  # the passes' call and return; prod's call (MOV x0, x1, x2, x3, BL: at most 8)
ACCT = 14 + 3 * 8  # 192-bit account A (4): deposit and check, then 3k after loop A, prod and pass 2
BATCH_FIXED = INV + FUSED_FIXED + FIX["prod"] + FIX["pass2"] + CALLS + 12 + 3 + 7 + 3 + 2 + 4 + ACCT
TAILS_ROUND = 127 * (WALKF + 4 - PASS1F) + 127 * (1 + WALKB + 2 - PROD) + 127 * (1 + WALK2 + 2 - PASS2)  # <= 127 tail walks per loop
ROUND_FIXED = ceil_fraction(BATCH_FIXED + TAILS_ROUND)
check("per-round work 24,552 + 1,140.0234375 -> 25,693, per round of 2^23 walks",
 BATCH_FIXED == 24552 and TAILS_ROUND == Fraction(145923, 128) and ROUND_FIXED == 25693
 and FUSED_FIXED + 127 * (WALKN + 3 - PASS1N) <= FUSED_FIXED + 127 * (WALKF + 4 - PASS1F),  # pass1p's tail charge is below pass1f's
 str(ROUND_FIXED))
SHARE = Fraction(ROUND_FIXED, 2 ** 23)

CORR = 3 * (Fraction(3, 2) + Fraction(1, 2 ** 9) + 2 * Fraction(1, 2 ** 10))  # 3 per borrow x sum over dx, dy, u of (1/2 + pad)
PREM = 3 * (Fraction(3, 2) + Fraction(1, 2 ** 10) + 2 * Fraction(1, 2 ** 11))  # the same at the premise levels 1/2 + pi0 (5.4)
STEP = PASS1N + PROD + PASS2 + CORR + SHARE
DBL_ARITH = 4 * MUL_R + MUL_SUB + SQR_R + SQR_SUB2 + ADD_WEAK + 2 * ADD_WEAK2 + SUB_R + 30
DBL = STEP + DBL_ARITH + 34 + 3 + SHARE
TRAV = STEP + 24 + 2 + 7 + 6
PRICE_STEP = Fraction(69699, 128) + Fraction(1155, 256) + SHARE
PRICE_TRAV, PRICE_DBL = PRICE_STEP + 39, PRICE_STEP + 677
check("step is 37.5078125 + 169.5078125 + 337.5078125 + CORR 1155/256 + share = 549.03515625 + 25,693/2^23",
 STEP == PRICE_STEP and CORR == Fraction(1155, 256), "%.7f" % float(STEP))
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
check("stubs: loops 1,935 = 129*(1 + 6 + 8), loop A 258 = 2*129", LABELS == {"secp_loops.s": 1935, "secp_fuse.s": 258}
 and 129 * (1 + 6 + 8) == 1935 and 2 * 129 == 258)
check("stubs loop A 7, SB 9, SC 8, ST 13: all <= 13 instructions <= 20 milli-units (2.5)",
 (F, LB, LC, LT) == (7, 9, 8, 13) and max(F, LB, LC, LT) <= 13 <= STUB, "%d %d %d %d" % (F, LB, LC, LT))
check("rare path 5 + 5 + 5 + 4 + 3 = 22, 23 with loop A's B, <= 25 milli-units (2.5, 3.2)",
 5 + 5 + 5 + 4 + 3 + 1 == 23 <= RARE)
check("a walk's 15 stubs = 300; 15*20 + 25 = 325 <= the 336 reserve (2.5, 4)", WALK_STUBS * STUB == 300
 and WALK_STUBS * STUB + RARE == 325 <= RESERVE)
FINAL_UNITS, FINAL_CHARGE = 2 ** 18, 2 ** 18 * 1000 + 2 ** 23 * 25  # 2^18 units plus 25 milli-units per walk
check("final arrival pass: 2^18 units + 25 per walk cover 2^23 x 31.015625 + 2^10 + one entry each, < 2^18.85 units (2.5)",
 FINAL_CHARGE >= 2 ** 23 * ROWS_TOTAL + 2 ** 10 + 2 ** 23 * 25 and log2(Fraction(FINAL_CHARGE, 1000)) < 18.85
 and log2(Fraction(FINAL_CHARGE, 1000)) < 70,
 "2^%.5f" % log2(Fraction(FINAL_CHARGE, 1000)))
CMP_ROUND = (PASS1F + PROD + PASS2 + CORR + RESERVE) * 2 ** 23 + ROUND_FIXED
ANCHOR_ROUND = CMP_ROUND + 18 * 2 ** 23  # the anchor variant's 16 per walk (3.2) and each walk's checkpoint write, 2
check("anchor round charge (552.53515625 + 336 + 18)*2^23 + 25,693 < 2^32.83, +18*2^23 = MOVZ #0x900, LSL #16; round 0 + 9n < 2^33 (4)",
 PASS1F + PROD + PASS2 + CORR == Fraction(141449, 256) and 0x900 << 16 == 18 * 2 ** 23 and 12 + (ROWS_TOTAL + PASS1 - PASS1F) + Fraction(31, 64) == 16
 and CMP_ROUND < ANCHOR_ROUND < pow2_interval(32830, 1000)[0] and ANCHOR_ROUND + 9 * 2 ** 23 < 2 ** 33, "2^%.4f" % log2(ANCHOR_ROUND))
check("MOVZ #0xa800, #0x4600 (LSL #16) = 336*2^23, 140*2^23; ADD #0x800, LSL #12 = 2^23",
 336 * 2 ** 23 == 0xA800 << 16 and 140 * 2 ** 23 == 0x4600 << 16 and 2 ** 23 == 0x800 << 12)
check("phase 1's last blocks cover < 2^60.87 further operations (4)",
 Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP < pow2_interval(60870, 1000)[0],
 "2^%.5f" % log2(Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP))
RERUN = 1 + WALK1 + 2 + 4 + 1 + WALKB + 2 + 3  # tail walks of pass1 (with the per-walk dx = 0 test) and prod, a dy borrow
check("zero-test re-run 2^23*(195 + 7*20) + fixed < 2^22 units, stubs and dy borrows included (4)",
 RERUN == 195 and 2 ** 23 * RERUN < 2 ** 31 and 2 ** 23 * (RERUN + RERUN_STUBS * STUB) + 2 ** 10 <= 2 ** 22 * 1000)
LMAX = 20 * 2 ** 64  # slot-steps of one trail are fewer (3.2, 5.1)
LAM, CB = Fraction(1, 2 ** 80), 9 * LMAX  # Bennett's lambda; X_k - E[X_k | past] <= 9 L_max
L2LO, L2HI = ln_interval(2)
DEV, INPROG = 140 * L2HI / LAM, 2 ** 23 * CORR * LMAX  # Ville at 2^-140; trails in progress, each <= CORR*L_max
check("CORR*n = 1155*2^15 = MOVZ #0x8000 + MOVK #0x241, LSL #16; 9n = MOVZ #0x480, LSL #16; rho' = 3/512 (4)",
 CORR * 2 ** 23 == 0x8000 + (0x241 << 16) and 9 * 2 ** 23 == 0x480 << 16 and CORR - PREM == Fraction(3, 512))
check("Bennett at lambda = 2^-80 absorbs the variance (<= rho'); Ville 2^-140 deviation + in-progress < 2^94 = Lambda_p (5.4)",
 LAM * CB < 3 and max(CORR, 9 - CORR) == CORR and LAM * CORR ** 2 * LMAX / (2 * (1 - LAM * CB / 3)) <= CORR - PREM
 and DEV + INPROG < 2 ** 94, "%.6f; 2^%.4f + 2^%.4f" % (float(LAM * CORR ** 2 * LMAX / (2 * (1 - LAM * CB / 3))), log2(DEV), log2(INPROG)))
DKW = 6 * Fraction(10 ** 9, 2718281828) ** 32  # >= 6 exp(-2 R eps^2) at R = 2^26, eps = 2^-11 (2.718281828 < e)
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
check("phase 1 < 2^125.834481", log2(Fraction(PH1, 1000)) < 125.834481,
 "%.9f" % log2(Fraction(PH1, 1000)))

CLAIM_LO = 2 ** 125 * pow2_frac(CLAIM - 125)[0]
CAP_MILLI_LO = (CLAIM_LO.numerator * 1000) // CLAIM_LO.denominator
BLOCK_MILLI = 2 ** 40 * 1000
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
PERIOD, STRIDE, TRAV_CAP = 8192, 256, 32  # anchors, compare rounds, traversal cap (5.1)
CP_MILLI = 2 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))  # <= 1 write (STP + STR) per anchor
ANCHOR_MILLI = 16 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))  # anchor work per walk: full K 7, counter 5, unfused 3.515625
CMP_MILLI = ceil_fraction((PASS1F - PASS1N) * ceil_fraction(Fraction(S_INT + 2 ** 61, STRIDE)))
HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
   ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
   ("phase 2", Fraction(2 ** 23 * (20 * 2 ** 64 + 65) * 600) * 1000),
   ("re-walks", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
   ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
   ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
   ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
   ("table (R entries) and verification, < 2^40 (R:)", Fraction(2 ** 40) * 1000),
   ("node set-up, 2^26", Fraction(2 ** 26) * 1000),
   ("checkpoint writes: 2*ceil((S_int + 2^61)/8192)", Fraction(CP_MILLI)),
   ("anchor rounds: 16*ceil((S_int + 2^61)/8192)", Fraction(ANCHOR_MILLI)),
   ("compare rounds: 3.5*ceil((S_int + 2^61)/256)", Fraction(CMP_MILLI)),
   ("correction seeds: 2^20*(2^94 + 9n)", Fraction(2 ** 20 * (2 ** 94 + 9 * 2 ** 23)))]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD = sum(int(v) for _, v in HARD_ITEMS)
check("hard extras < 2^110.883919 units", log2(Fraction(HARD, 1000)) < 110.883919,
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
TRAV_RATE = Fraction(2, 6 * R_PART) + TRAV_CAP * OTHER  # traversal steps per operation
check("R: traversal-step rate 2/(6R) + 32*sum <= 2^-27.5849", TRAV_RATE <= Fraction(1, 1) / pow2_interval(275849, 10000)[0],
 "2^%.9f" % log2(TRAV_RATE))
OPS = S_HI + 2 ** 60
check("trails ended by a longer cycle: LONG*ops < 2^-232 per run, LONG*ops*L_max < 2^-100 units",
 LONG * OPS < Fraction(1, 2 ** 232) and LONG * OPS * 20 * 2 ** 64 < Fraction(1, 2 ** 100))
# the strided compare (proof.md 5.1)
LCM = [L_ * STRIDE // math.gcd(L_, STRIDE) for L_ in range(TRAV_CAP + 1)]
check("anchors are compare rounds and lcm(L, 256) <= 7936 < 8192 for L = 1..32, not for L = 33", PERIOD % STRIDE == 0
 and max(LCM[1:]) == 7936 < PERIOD and 33 * STRIDE > PERIOD)
D2 = PERIOD - 1 + LCM[2] + 2
DL = max(PERIOD - 1 + LCM[L_] + L_ for L_ in range(3, TRAV_CAP + 1))
WASTE = Fraction(D2, 6 * R_PART) + DL * OTHER
W_KAPPA = pow2_interval(-1456, 100)[1]  # w = 2^-14.56 of the kappa rule (proof.md 5.2)
# proof.md 5.2's toy arithmetic and the delay term
check("toy abandonment 640*2^-13.6 = 0.0515 -> 0.052 (5.2); 8,447 > the toy cap 640",
 round(640 * 2.0 ** -13.6, 3) == 0.052 and 8447 > 640 and D2 - 2 == 8447)
check("R: delay term 8,447/(6R) < 2^-15.5407 < w = 2^-14.56 of the kappa rule",
 Fraction(8447, 6 * R_PART) < pow2_interval(-155407, 10000)[0] and -15.5407 < -14.56, "2^%.6f" % log2(Fraction(8447, 6 * R_PART)))
check("R: waste w <= 8449/(6R) + 16158*sum <= 2^-15.54 <= 2^-14.56", D2 == 8449 and DL == 16158
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
check("kappa rule, eleven walk runs pooled: (1 + w)(kappa_eff + 2 s.e.) <= 1.0001", len(WALKS) == 11
 and SE_HI * SE_HI * INV > 1 and (1 + W_KAPPA) ** 2 * (1 + SE_HI / RHO) ** 2 / RHO <= KAPPA ** 2,
 "%.6f" % ((1 + float(W_KAPPA)) * float(RHO) ** -0.5 * (1 + float(SE_HI / RHO))))
# proof.md 5.2: no walk run exceeds kappa by one s.e.; the closest is q29 R 2^25 dp3
EXCESS = {k: (1 + W_KAPPA) * Fraction(str(v["kappa_eff"])) - KAPPA for k, v in RUNS.items() if k.startswith("walk_")}
NEAR = max(EXCESS, key=EXCESS.get)
check("no walk run exceeds 1.0001 by one s.e. (5.2); closest q29 R 2^25 dp3, 0.00026 = 0.93 s.e.",
 NEAR == "walk_q29_r25_dp3" and EXCESS[NEAR] < Fraction(27, 100000) and EXCESS[NEAR] >= Fraction(26, 100000)
 and all(EXCESS[k] <= Fraction(str(RUNS[k]["kappa_se"])) for k in EXCESS)
 and Fraction(93, 100) <= EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])) < Fraction(94, 100),
 "%.7f / %.3f s.e." % (float(EXCESS[NEAR]), float(EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])))))
# rare exits per operation (proof.md 5.3)
RATE_K = 25 * Fraction(2 ** 32 + 977, 2 ** 64) + 2 * Fraction(1, 2 ** 63)  # 3 fe_sub, 6 rows, 4 row folds (4C0), 2 folds
check("kernel rare-exit rate 25*C0/2^64 + 2*2^-63 <= 2^-27.356 per operation",
 RATE_K <= pow2_interval(-27356, 1000)[0], "2^%.5f" % log2(RATE_K))
HITS = 2 * RATE_HI  # true filter hits: the doubling's and the closing arrival's (proof.md 5.3)
check("R: true filter hits 2*(1/(6R) + sum) <= 2^-27.5849 per operation", HITS <= Fraction(1, 1) / pow2_interval(275849, 10000)[0],
 "2^%.9f" % log2(HITS))
EXTRA_ITEMS = [ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * RATE_HI * OPS * 1000),
   ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * TRAV_RATE * OPS * 1000),
   ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 23) * 1000),
   ceil_fraction(Fraction(2 ** 70) * 1000),
   ceil_fraction(Fraction(25) * (HITS + Fraction(1, 2 ** 56) + Fraction(2 ** 32 + 977, 2 ** 64)) * OPS),  # rare-path entries
   ceil_fraction(Fraction(20) * RATE_K * (S_HI + 2 ** 61 + RATE_HI * OPS))]  # rare exits, 0.02 units; + a doubling's
E = sum(EXTRA_ITEMS)
check("expected random extras < 2^97.8873 units", log2(Fraction(E, 1000)) < 97.8873,
 "2^%.9f" % log2(Fraction(E, 1000)))

# --- 7. allowance, total and claim
A = ceil_fraction(E / HEADROOM)
check("allowance = ceil(E/0.000048) (Markov safe)", A * 48 >= E * 1000000)
W_LEDGER = PH1 + HARD + A
check("W_ledger < 2^125.8346426", log2(Fraction(W_LEDGER, 1000)) < 125.8346426,
 "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("125.834643 holds: W_ledger <= floor(1000 * 2^125 * (lower bound of 2^0.834643))", W_LEDGER <= CAP_MILLI_LO)
NEXT_HI = 1000 * 2 ** 125 * pow2_frac(CLAIM - Fraction(1, 10 ** 6) - 125)[1]
check("125.834642 fails: W_ledger > 1000 * 2^125 * (upper bound of 2^0.834642)", W_LEDGER > NEXT_HI)

HEAD = CAP_MILLI_LO - PH1 - HARD
check("headroom > allowance", HEAD > A, "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio < 0.000048", MARKOV_HI < HEADROOM and MARKOV_HI * 1000000 < 48,
 "%d/10^10 (rounded up)" % ceil_fraction(MARKOV_HI * 10 ** 10))
SUCCESS_LO = P0 - HEADROOM - LOSS_OTHER
check("success = 0.3900482 - 0.000048 - 1.1e-7 = 0.39000009 > 0.39",
 SUCCESS_LO == Fraction(3900000900, 10 ** 10) and SUCCESS_LO > THRESHOLD,
 "%.7f" % float(SUCCESS_LO))

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
ROWS_FAIL = UPSTREAM_ROWS + ROW_CARRY + ROW_NEG + ROW_HIGH + ROW_CANON + DKW + Fraction(1, 2 ** 120)
check("failure rows other than the early halt, the account's included, < 1.0104e-7 < 1.1e-7",
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
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(391, 10)[1] + 2 ** 36  # table 2^39.1, first starts 2^36
CAP_PREP = pow2_interval(590, 10)[1]
check("expected preprocessing total < 2^59.0 units", PREP_EXPECTED < CAP_PREP,
 "2^%.5f" % log2(PREP_EXPECTED))
check("preprocessing cap halts with probability < 2^-30",
 Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS) < pow2_interval(-30, 1)[1],
 "2^%.4f" % log2(Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS)))


print()
print("step %.6f instructions" % STEP)
print("phase 1 2^%.6f  hard 2^%.6f  extras 2^%.6f  allowance 2^%.6f"
 % (log2(Fraction(PH1, 1000)), log2(Fraction(HARD, 1000)), log2(Fraction(E, 1000)),
  log2(Fraction(A, 1000))))
print("W_ledger 2^%.9f units  claim %s  head 2^%.6f"
 % (log2(Fraction(W_LEDGER, 1000)), CLAIM, log2(Fraction(HEAD, 1000))))
print("Markov < %d/10^10  success > %.7f" % (ceil_fraction(MARKOV_HI * 10 ** 10), float(SUCCESS_LO)))
print()
if FAILED:
 print("CERTIFICATE: FAIL (%d)" % len(FAILED))
 raise SystemExit(1)
print("CERTIFICATE: PASS")
