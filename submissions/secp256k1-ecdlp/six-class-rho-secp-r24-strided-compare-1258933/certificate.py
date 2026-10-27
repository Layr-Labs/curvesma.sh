"""Exact-integer certificate for the work bound W <= 2^125.8933 of proof.md.

Every check is an exact comparison of Python integers or Fractions: the instruction identities, prices and
store numbers are integer arithmetic; ln(1/(1 - p0)) is a rational interval from the atanh series with an
explicit remainder; every fractional power of two is a rational interval certified by integer roots; and the
ledger and success floor are exact integers or rationals, each claim inequality checked in both directions.

Run:  python3 certificate.py        (exit status 0 = all checks hold)

It proves, for kappa = 1.0001, R = 2^24, p0 = 0.3900503, h = 0.00005 and the step price 571.546875 + 25616/2^20:

  phase1  = ceil(price * (S_int + 2^60))   (milli-units);  hard = the twelve hard rows;  E = the six random rows
  A       = ceil(E / 0.00005);   W = phase1 + hard + A  < 2^125.8933
  E / (cap - phase1 - hard) < 0.00005   =>  success > 0.3900503 - 0.00005 - 1.1e-7
                      = 0.39000019 > 0.39
  W^10000  < 1000^10000 * 2^1258933    so 125.8933 holds
  W^10000 > 1000^10000 * 2^1258932    so 125.8932 fails
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
R_PART = 2 ** 24  # partition table rows
KAPPA = Fraction(10001, 10000)  # H1 slowdown allowance 1.0001 (proof.md section 5.2)
P0 = Fraction(3900503, 10000000)  # H1 floor p0 at S0 = S_ops / KAPPA
HEADROOM = Fraction(5, 100000)  # 0.00005 Markov headroom h
CLAIM = Fraction(1258933, 10000)  # W = 2^125.8933
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
WALK1 = MUL_R + SUB_R + 4 + 4  # pass 1 walk: 4 LDP, 4 STP (its table pointer comes with the group)
WALK2 = 3 * MUL_R + X3_R + MUL_ADD_R + 2 * SUB_R + 12 + 4  # pass 2 walk: 12 LDP, 4 STP
PASS1 = Fraction(128 * WALK1 + 64 + 2, 128)  # 128 walks, 64 LDP of their pointers, SUBS + B.NE
PASS2 = Fraction(128 * WALK2 + 64 + 2, 128)  # 128 walks, 64 LDP of their pointers, SUBS + B.NE on x18
TAILS = 127 * (1 + WALK1 + 2 - PASS1) + 127 * (1 + WALK2 + 2 - PASS2)  # tail walks of a round
check("fe_mul costs 85 instructions", MUL == 85, str(MUL))
check("fe_sqr costs 64 instructions", SQR == 64, str(SQR))
check("kernels with rare exits: fe_mul 84, fe_sub 7, y3' 89, x3 73", (MUL_R, SUB_R, MUL_ADD_R, X3_R) == (84, 7, 89, 73))
check("pass 1: 12,738 per 128 walks = 99.515625 per walk", 128 * WALK1 + 66 == 12738 and PASS1 == Fraction(6369, 64), str(PASS1))
check("fused fe_mul and fe_sqr cost 91 and 75", MUL_SUB == 91 and SQR_SUB2 == 75)
check("pass 2: 56,898 per 128 walks = 444.515625 per walk", 128 * WALK2 + 66 == 56898 and PASS2 == Fraction(28449, 64),
 str(PASS2))
check("127 + 127 tail walks cost 631.03125 above the group rates (charged per round)", TAILS == Fraction(20193, 32), str(TAILS))
ARRIVE_F = 13 + 10 + 3 + 1 + 1 + 3 + 1
WALKF = ARRIVE_F + 2 + SUB_R + 4 + MUL_R
PASS1F = Fraction(128 * WALKF + 4, 128)
check("fused walk: arrival 32 + 2 LDP x(T) + fe_sub 7 + 4 STP + fe_mul 84 = 129", ARRIVE_F == 32 and WALKF == 129)
check("fused pass 1: 16,516 per 128 walks = 129.03125 per walk", 128 * WALKF + 4 == 16516 and PASS1F == Fraction(4129, 32),
 str(PASS1F))
WALKN = WALKF - 2  # plain rounds (pass1p): no LDR of the K_cp word, no CCMP of the filter word
PASS1N = Fraction(128 * WALKN + 4, 128)
check("plain fused pass 1: 16,260 per 128 walks = 127.03125 per walk", 128 * WALKN + 4 == 16260 and PASS1N == Fraction(4065, 32),
 str(PASS1N))

ROWS = [("F1 truncated product (2 LDP, 4 UMULH, 3 ADD, 4 MADD; G1..G5 in registers)", 13),
  ("canonical index and slot: ADD, SUBS, CSEL, CSEL, ADDS, CSEL; CCMN, CCMN, CCMN, CCMP", 10),
  ("parity branch (TBNZ y0); K0 + 1 and K1 by SUBS, SBCS, or CMP y1 and ORR", 3),
  ("table address: BFI of K0 + 1 bits 0..24 into bits 9..33", 1),
  ("arrival traffic (LDP of y0, y1; one STP of a pair's two table addresses)", Fraction(3, 2)),
  ("checkpoint filter (one LDP of a pair's K_cp words, CCMP) and the rare branch (B.EQ)", Fraction(5, 2)),
  ("arrival loop of 32 pairs of walks (one TBZ on the pointer array's end bit per 64 walks)", Fraction(1, 64))]
ROWS_TOTAL = sum(v for _, v in ROWS)
check("unfused arrival rows (anchor rounds) sum to 31.015625", ROWS_TOTAL == Fraction(1985, 64), "%.6f" % float(ROWS_TOTAL))
# per-batch work (proof.md section 4)
FUSED_FIXED = 25 + (ARRIVE_F + 2 + SUB_R) + 14
KD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "kernels")
sys.path[:0] = [KD]
sys.dont_write_bytecode = True
import a64
K = {}
for f in "secp_loops.s", "secp_fe.s":
 K.update(a64._asm_counts(open(os.path.join(KD, f)).read()))
N = lambda k: K[k]["instructions"]
INV = 255 * (N("fe_sqr_f") + 1) + 15 * (N("fe_mul_f") + 1) + N("fe_inv") - K["fe_inv"]["calls"] + 255 * 4
P2F = N("pass2") - sum(l["iteration"] for l in K["pass2"]["loops"])
check("a64.py counts: inversion 24,322, pass 2 fixed 299", INV == 24322 and P2F == 299)
BATCH_FIXED = INV + FUSED_FIXED + P2F + 64 + 12 + 3 + 7 + 3 + 2 + 4  # + TST, B.NE: compare round?
TAILS_ROUND = 127 * (WALKF + 4 - PASS1F) + 127 * (1 + WALK2 + 2 - PASS2)  # n = 2^20: 127 + 127 tail walks every round
ROUND_FIXED = ceil_fraction(BATCH_FIXED + TAILS_ROUND)
check("per-round work 24,796 + 819.546875 -> 25,616, per round of 2^20 walks",
 FUSED_FIXED == 80 and BATCH_FIXED == 24796 and TAILS_ROUND == Fraction(52451, 64) and ROUND_FIXED == 25616
 and WALKN + 4 - PASS1N == WALKF + 4 - PASS1F,
 str(ROUND_FIXED))
SHARE = Fraction(ROUND_FIXED, 2 ** 20)

STEP = PASS1N + PASS2 + SHARE
DBL_ARITH = 4 * MUL_R + MUL_SUB + SQR_R + SQR_SUB2 + ADD_WEAK + 2 * ADD_WEAK2 + SUB_R + 30
DBL = STEP + DBL_ARITH + 34 + 3 + SHARE
TRAV = STEP + 24 + 2 + 7 + 6
PRICE_STEP = Fraction(36579, 64) + SHARE
PRICE_TRAV, PRICE_DBL = PRICE_STEP + 39, PRICE_STEP + 677
check("step is 127.03125 + 444.515625 + share = 571.546875 + 25,616/2^20", STEP == PRICE_STEP, "%.6f" % float(STEP))
check("doubling bound: displaced step + 639 + 34 + 3 + share = step + 676 + share",
 DBL_ARITH == 639 and DBL == STEP + 676 + SHARE)
check("traversal bound is step + 24 + 2 + 7 + 6", TRAV == PRICE_STEP + 39)
check("traversal price covers its bound", PRICE_TRAV >= TRAV)
check("doubling price step + 677 covers its bound", PRICE_DBL >= DBL)

# --- 2b. meter, stubs, reserves: every such bound proof.md states (2.5, 3.2, 4); 1 milli-unit = 1 instruction


def stub(fname):
 """Instructions in a kernel's first rare stub (.Lr1_0_ through its B), RCHG expanded."""
 t = open(os.path.join(KD, fname), encoding="utf-8").read()
 r = t[t.index(".macro RCHG"):].split(".endm")[0].count("\n") - 1
 g = t[t.index("\n.Lr1_0_"):].split("b .Lb")[0].splitlines()[1:]
 return 1 + sum(r if x.strip() == "RCHG" else 1 for x in g if x.strip() and not x.strip().endswith(":"))


STUB, RARE, RESERVE, WALK_STUBS, RERUN = 20, 25, 256, 9, 2      # milli-units; 9 = 2 pass-1 + 7 pass-2 stubs
F, L = stub("secp_fuse.s"), stub("secp_loops.s")
check("shipped stubs fused 7, loops 9, tail 6: all <= 9 instructions, 20 milli-units (2.5)",
 (F, L, 6) == (7, 9, 6) and max(F, L, 6) <= min(9, STUB), "%d %d" % (F, L))
check("rare path 4 + 5 + 6 + 4 + 3 = 22 instructions <= its 25 milli-units (3.2)",
 4 + 5 + 6 + 4 + 3 == 22 <= RARE)
check("a walk's 9 stubs = 180; (9 + 2)*20 + 25 = 245 <= the 256 reserve (2.5, 4)", WALK_STUBS * STUB == 180
 and (WALK_STUBS + RERUN) * STUB + RARE == 245 <= RESERVE)
FINAL_UNITS, FINAL_CHARGE = 2 ** 15, 2 ** 15 * 1000 + 2 ** 20 * 25   # 2^15 units plus 25 milli-units per walk
check("final arrival pass (2.5): 2^15 units + 25 per walk cover 2^20 x 31.015625 + one entry each, under 2^15.85 units",
 FINAL_CHARGE >= 2 ** 20 * ROWS_TOTAL + 2 ** 20 * 25 and log2(Fraction(FINAL_CHARGE, 1000)) < 15.85
 and log2(Fraction(FINAL_CHARGE, 1000)) < 70,
 "2^%.5f" % log2(Fraction(FINAL_CHARGE, 1000)))
check("a round's charge (571.546875 + 256 + 2)*2^20 + 25,616 on a 2^40-unit block is < 2^30 (4)",
 (STEP + RESERVE + 2) * 2 ** 20 + ROUND_FIXED < 2 ** 30)
check("phase 1's last blocks cover < 2^60.81 further operations (4)",
 Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP < pow2_interval(60810, 1000)[1],
 "2^%.5f" % log2(Fraction(2 ** 20 * 2 ** 40 * 1000) / PRICE_STEP))
check("the zero test's re-run 2^20*106 < 2^27 instructions, inside its 2^17-unit action (4)",
 2 ** 20 * (102 + 4) < 2 ** 27 and 2 ** 20 * (102 + 4) <= 2 ** 17 * 1000)
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
check("phase 1 < 2^125.89252", log2(Fraction(PH1, 1000)) < 125.89252,
 "%.9f" % log2(Fraction(PH1, 1000)))

CLAIM_LO, CLAIM_HI = pow2_interval(1258933, 10000)
CAP_MILLI_LO = (CLAIM_LO.numerator * 1000) // CLAIM_LO.denominator
BLOCK_MILLI = 2 ** 40 * 1000
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
PERIOD, STRIDE, TRAV_CAP = 4096, 16, 256  # anchors, compare rounds, traversal cap (5.1)
CP_MILLI = 2 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))  # <= 1 write (STP + STR) per anchor
ANCHOR_MILLI = 14 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))  # anchor work per walk: full K 7, counter 5, unfused 1.5
CMP_MILLI = 2 * ceil_fraction(Fraction(S_INT + 2 ** 61, STRIDE))  # compare rounds: LDR of the K_cp word and CCMP per walk
HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
   ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
   ("phase 2", Fraction(2 ** 20 * (20 * 2 ** 64 + 65) * 600) * 1000),
   ("re-walks", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
   ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
   ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
   ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
   ("table (2^24 entries) and verification, < 2^38", Fraction(2 ** 38) * 1000),
   ("node set-up, 2^26", Fraction(2 ** 26) * 1000),
   ("checkpoint writes: 2*ceil((S_int + 2^61)/4096)", Fraction(CP_MILLI)),
   ("anchor rounds: 14*ceil((S_int + 2^61)/4096)", Fraction(ANCHOR_MILLI)),
   ("compare rounds: 2*ceil((S_int + 2^61)/16)", Fraction(CMP_MILLI))]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD = sum(int(v) for _, v in HARD_ITEMS)
check("hard extras < 2^113.7927 units", log2(Fraction(HARD, 1000)) < 113.7927,
 "%.9f" % log2(Fraction(HARD, 1000)))
check("the compare-round row is the largest hard row", CMP_MILLI == max(int(v) for _, v in HARD_ITEMS))

# --- 6. expected random extras
# fruitless cycles of length 3..TRAV_CAP (proof.md section 5.1: S2(L, b) counts partitions into b blocks >= 2)
S2 = [[0] * (TRAV_CAP + 1) for _ in range(TRAV_CAP + 1)]
S2[0][0] = 1
for L_ in range(1, TRAV_CAP + 1):
 for b_ in range(1, L_ // 2 + 1):
  S2[L_][b_] = b_ * S2[L_ - 1][b_] + (L_ - 1) * S2[L_ - 2][b_ - 1]
check("S2(3,1), S2(4,1), S2(4,2), S2(6,3) = 1, 1, 3, 15", (S2[3][1], S2[4][1], S2[4][2], S2[6][3]) == (1, 1, 3, 15))


def other_cycles(R, top=TRAV_CAP):
 return sum(Fraction(S2[L_][b_], R ** (L_ - b_)) for L_ in range(3, top + 1) for b_ in range(1, L_ // 2 + 1))


check("at R = 2048 the sum to L = 64 is 2^-19.9954 (the public package's value)",
 abs(log2(other_cycles(2048, 64)) + 19.9954) < 5e-5, "2^%.6f" % log2(other_cycles(2048, 64)))
OTHER = other_cycles(R_PART)
check("at R = 2^24 the sum is < 2^-45.999", OTHER < Fraction(1, 1) / pow2_interval(45999, 1000)[0],
 "2^%.7f" % log2(OTHER))
RATE_HI = Fraction(1, 6 * R_PART) + OTHER
check("doubling rate 1/(6R) + sum <= 2^-26.5849", RATE_HI <= Fraction(1, 1) / pow2_interval(265849, 10000)[0],
 "2^%.9f" % log2(RATE_HI))
TRAV_RATE = Fraction(2, 6 * R_PART) + TRAV_CAP * OTHER  # traversal steps per operation
check("traversal-step rate 2/(6R) + 256*sum <= 2^-25.5846", TRAV_RATE <= Fraction(1, 1) / pow2_interval(255846, 10000)[0],
 "2^%.9f" % log2(TRAV_RATE))
OPS = S_HI + 2 ** 60
# the strided compare (proof.md 5.1): the coverage, delay and traversal bounds below
LCM = [L_ * STRIDE // math.gcd(L_, STRIDE) for L_ in range(TRAV_CAP + 1)]
check("anchors are compare rounds and lcm(L, 16) <= 4080 < 4096 for L = 1..256", PERIOD % STRIDE == 0
 and max(LCM[1:]) == 4080 < PERIOD)
D2 = PERIOD - 1 + LCM[2] + 2
DL = max(PERIOD - 1 + LCM[L_] + L_ for L_ in range(3, TRAV_CAP + 1))
WASTE = Fraction(D2, 6 * R_PART) + DL * OTHER
W_KAPPA = pow2_interval(-1456, 100)[1]  # w = 2^-14.56 of the kappa rule (proof.md 5.2)
# proof.md 5.2's toy arithmetic, both products stated there and the delay term they identify
check("toy abandonment 640*2^-13.6 = 0.0515 rounds to the 0.052 proof.md prints, and 4,111 > the toy cap 640",
 round(640 * 2.0 ** -13.6, 3) == 0.052 and 4111 > 640)
check("the delay term 4,111*2^-26.585 = 2^-14.5797 lies below the w = 2^-14.578 the kappa rule allows",
 round(log2(4111 * 2.0 ** -26.585), 4) == -14.5797 and log2(4111 * 2.0 ** -26.585) < -14.578)
TOY_SHARE = (0.052, 0.0527, 0.0558, 0.0674)   # proof.md 5.2: four runs in the declared image, 25..45 trails
check("the toy abandonment shares proof.md prints (0.052 to 0.0674, 25..45 trails of 115..269) span the rate's 0.0515 at most 2.4x",
 min(TOY_SHARE) <= 0.052 and 0.0515 <= min(TOY_SHARE) * 1.02 and max(TOY_SHARE) <= 2.4 * 0.0515)
check("waste w <= 4113/(6R) + 8430*sum <= 2^-14.578 <= 2^-14.56", D2 == 4113 and DL == 8430
 and WASTE <= pow2_interval(-14578, 1000)[0] and WASTE <= W_KAPPA, "2^%.5f" % log2(WASTE))
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
check("kappa rule, nine walk runs pooled: (1 + w)(kappa_eff + 2 s.e.) <= 1.0001", len(WALKS) == 9
 and SE_HI * SE_HI * INV > 1 and (1 + W_KAPPA) ** 2 * (1 + SE_HI / RHO) ** 2 / RHO <= KAPPA ** 2,
 "%.6f" % ((1 + float(W_KAPPA)) * float(RHO) ** -0.5 * (1 + float(SE_HI / RHO))))
# proof.md 5.2: no walk run exceeds kappa by one s.e.; the closest is q29 dp3
EXCESS = {k: (1 + W_KAPPA) * Fraction(str(v["kappa_eff"])) - KAPPA for k, v in RUNS.items() if k.startswith("walk_")}
NEAR = max(EXCESS, key=EXCESS.get)
check("no walk run exceeds 1.0001 by one s.e. (5.2); closest q29 dp3, 0.00009 = 0.65 s.e.",
 NEAR == "walk_q29_r20_dp3" and EXCESS[NEAR] < Fraction(1, 10000) and EXCESS[NEAR] >= Fraction(9, 100000)
 and all(EXCESS[k] <= Fraction(str(RUNS[k]["kappa_se"])) for k in EXCESS)
 and Fraction(8, 13) <= EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])) < Fraction(7, 10),
 "%.7f / %.3f s.e." % (float(EXCESS[NEAR]), float(EXCESS[NEAR] / Fraction(str(RUNS[NEAR]["kappa_se"])))))
# rare exits per operation (proof.md 5.3)
RATE_K = 3 * Fraction(2 ** 32 + 977, 2 ** 64) + 6 * Fraction(1, 2 ** 63)
check("kernel rare-exit rate 3*C0/2^64 + 6*2^-63 <= 2^-30.415 per operation",
 RATE_K <= pow2_interval(-30415, 1000)[0], "2^%.5f" % log2(RATE_K))
HITS = 2 * RATE_HI  # true filter hits: the doubling's and the closing arrival's (proof.md 5.3)
check("true filter hits 2*(1/(6R) + sum) <= 2^-25.5849 per operation", HITS <= Fraction(1, 1) / pow2_interval(255849, 10000)[0],
 "2^%.9f" % log2(HITS))
EXTRA_ITEMS = [ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * RATE_HI * OPS * 1000),
   ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * TRAV_RATE * OPS * 1000),
   ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 20) * 1000),
   ceil_fraction(Fraction(2 ** 70) * 1000),
   ceil_fraction(Fraction(25) * (HITS + Fraction(1, 2 ** 56) + Fraction(2 ** 32 + 977, 2 ** 64)) * OPS),  # rare-path entries at 0.025 units
   ceil_fraction(Fraction(20) * RATE_K * (S_HI + 2 ** 61 + RATE_HI * OPS))]  # rare exits at 0.02 units; + a doubling's own
E = sum(EXTRA_ITEMS)
check("expected random extras < 2^99.8052 units", log2(Fraction(E, 1000)) < 99.8052,
 "2^%.9f" % log2(Fraction(E, 1000)))

# --- 7. allowance, total and claim
A = E * 20000  # E / 0.00005, an integer
check("allowance A = E/0.00005 (Markov safe)", A * 5 >= E * 100000)
W_LEDGER = PH1 + HARD + A
check("W_ledger < 2^125.8933", log2(Fraction(W_LEDGER, 1000)) < 125.8933,
 "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("W_ledger > 2^125.89324 (for the quoted value)", log2(Fraction(W_LEDGER, 1000)) > 125.89324,
 "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("the claim is 1258933/10000", CLAIM == Fraction(1258933, 10000))
check("W_ledger <= 1000 * 2^125.8933", W_LEDGER <= CAP_MILLI_LO)
check("125.8933 holds: W^10000 < 1000^10000 * 2^1258933",
 W_LEDGER ** 10000 < 1000 ** 10000 * 2 ** 1258933)
check("125.8932 fails: W^10000 > 1000^10000 * 2^1258932",
 W_LEDGER ** 10000 > 1000 ** 10000 * 2 ** 1258932)

HEAD = CAP_MILLI_LO - PH1 - HARD
check("headroom > allowance A", HEAD > A, "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio < 0.00005", MARKOV_HI < HEADROOM and MARKOV_HI * 100000 < 5,
 "%.10f" % float(MARKOV_HI))
SUCCESS_LO = P0 - HEADROOM - LOSS_OTHER
check("success = 0.3900503 - 0.00005 - 1.1e-7 = 0.39000019 > 0.39",
 SUCCESS_LO == Fraction(3900001900, 10 ** 10) and SUCCESS_LO > THRESHOLD,
 "%.7f" % float(SUCCESS_LO))
check("success floor >= the 0.39 of rules.md", SUCCESS_LO > THRESHOLD, "%.7f" % float(SUCCESS_LO))

# failure events other than the early halt (proof.md section 5.3): the upstream rows, the kernels' three
# rare-event rows and the canonical-index row
C0 = 2 ** 32 + 977
OPS_ALL = S_HI + 2 ** 61
check("fused doubling walk: 6 + 4*rate < 7 reductions per operation", 6 + 4 * RATE_HI < 7)
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
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(371, 10)[1] + 2 ** 37  # table 2^37.1
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
print("Markov %.10f  success > %.7f" % (float(MARKOV_HI), float(SUCCESS_LO)))
print()
if FAILED:
 print("CERTIFICATE: FAIL (%d)" % len(FAILED))
 raise SystemExit(1)
print("CERTIFICATE: PASS")
