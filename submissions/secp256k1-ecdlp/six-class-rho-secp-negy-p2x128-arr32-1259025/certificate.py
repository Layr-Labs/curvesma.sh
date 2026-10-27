"""Exact-integer certificate for the work bound W <= 2^125.9019 of proof.md.

Every check is an exact comparison of Python integers or Fractions: the instruction identities, prices and
store numbers are integer arithmetic; ln(1/(1 - p0)) is a rational interval from the atanh series with an
explicit remainder; every fractional power of two is a rational interval certified by integer roots; and the
ledger and success floor are exact integers or rationals, each claim inequality checked in both directions.

Run:  python3 certificate.py        (exit status 0 = all checks hold)

It proves, for kappa = 1.00017, R = 2^22, p0 = 0.3901303, h = 0.00013 and the step price 575.046875 + 25391/2^21:

    phase1  = ceil(price * (S_int + 2^60))   (milli-units);  hard = the eleven hard rows;  E = the six random rows
    A       = ceil(E / 0.00013);  W = phase1 + hard + A  < 2^125.9019
    E / (cap - phase1 - hard) < 0.00013  =>  success > 0.3901303 - 0.00013 - 1.1e-7
                                            = 0.39000019 > 0.39
    W^10000  < 1000^10000 * 2^1259019    so 125.9019 holds
    W^10000 > 1000^10000 * 2^1259018    so 125.9018 fails
"""
import math
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


r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
P_FIELD = int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
R_PART = 2 ** 22
KAPPA = Fraction(100017, 100000)
P0 = Fraction(3901303, 10000000)
HEADROOM = Fraction(13, 100000)
CLAIM = Fraction(1259019, 10000)
THRESHOLD = Fraction(39, 100)
LOSS_OTHER = Fraction(11, 10 ** 8)

MUL = (16 + 16) + 30 + 23
SQR = (10 + 10) + 8 + 6 + 7 + 23
MUL_SUB = (16 + 16) + 30 + 24 + 5
SQR_SUB2 = (10 + 10) + 8 + 6 + 7 + 24 + 10
SUB, CANON, ADD_WEAK, ADD_WEAK2 = 9, 8, 9, 14
MUL_ADD = (16 + 16) + 30 + 23 + 5
X3 = (10 + 10) + 8 + 6 + 7 + 23 + 10
MUL_R, SUB_R, MUL_ADD_R, X3_R = MUL - 1, SUB - 2, MUL_ADD - 1, X3 - 1
SQR_R = SQR - 1
WALK1 = MUL_R + SUB_R + 4 + 4
WALK2 = 3 * MUL_R + X3_R + MUL_ADD_R + 2 * SUB_R + 12 + 4
PASS1 = Fraction(128 * WALK1 + 64 + 2, 128)
PASS2 = Fraction(128 * WALK2 + 64 + 2, 128)
TAILS = 127 * (1 + WALK1 + 2 - PASS1) + 127 * (1 + WALK2 + 2 - PASS2)
check("fe_mul costs 85 instructions", MUL == 85, str(MUL))
check("fe_sqr costs 64 instructions", SQR == 64, str(SQR))
check("kernels with rare exits: fe_mul 84, fe_sub 7, y3' 89, x3 73", (MUL_R, SUB_R, MUL_ADD_R, X3_R) == (84, 7, 89, 73))
check("pass 1: 12,738 per 128 walks = 99.515625 per walk", 128 * WALK1 + 66 == 12738 and PASS1 == Fraction(6369, 64), str(PASS1))
check("fused fe_mul and fe_sqr cost 91 and 75", MUL_SUB == 91 and SQR_SUB2 == 75)
check("pass 2: 56,898 per 128 walks = 444.515625 per walk", 128 * WALK2 + 66 == 56898 and PASS2 == Fraction(28449, 64),
   str(PASS2))
check("127 + 127 tail walks cost 631.03125 above the group rates (charged per round)", TAILS == Fraction(20193, 32), str(TAILS))

ROWS = [("F1 truncated product (2 LDP, 4 UMULH, 3 ADD, 4 MADD; G1..G5 in registers)", 13),
    ("canonical index and slot: ADD, SUBS, CSEL, CSEL, ADDS, CSEL; CCMN, CCMN, CCMN, CCMP", 10),
    ("parity branch (TBNZ y0); K0 + 1 and K1 by SUBS, SBCS, or CMP y1 and ORR", 3),
    ("table address: BFI of K0 + 1 bits 0..20 into bits 9..29", 1),
    ("arrival traffic (LDP of y0, y1; one STP of a pair's two table addresses)", Fraction(3, 2)),
    ("checkpoint filter (one LDP of a pair's K_cp words, CCMP) and the rare branch (B.EQ)", Fraction(5, 2)),
    ("arrival loop of 32 pairs of walks (one TBZ on the pointer array's end bit per 64 walks)", Fraction(1, 64)),
    ("share of the per-round work of the one batch of 2^21 walks", Fraction(25391, 2 ** 21))]
ROWS_TOTAL = sum(v for _, v in ROWS)
BATCH_FIXED = 24322 + (12875 - 102 - 12738 + 57644 - 447 - 56898) + 64 + 12 + 3 + 3 + 7 + 3 + 4 + 7
TAILS_ROUND = 127 * (1 + WALK1 + 2 - PASS1) + 127 * (1 + WALK2 + 2 - PASS2)
ROUND_FIXED = ceil_fraction(BATCH_FIXED + TAILS_ROUND)
check("per-round work is 24,759 + 631.03125 -> 25,391 instructions, charged per round of the 2^21 walks",
   BATCH_FIXED == 24759 and ROUND_FIXED == 25391 and ROWS[-1][1] == Fraction(ROUND_FIXED, 2 ** 21), str(ROUND_FIXED))
check("arrival rows sum to 31.015625 + 25,391/2^21", ROWS_TOTAL == Fraction(1985, 64) + Fraction(25391, 2 ** 21),
   "%.6f" % float(ROWS_TOTAL))

STEP = PASS1 + PASS2 + ROWS_TOTAL
DBL_ARITH = 4 * MUL_R + MUL_SUB + SQR_R + SQR_SUB2 + ADD_WEAK + 2 * ADD_WEAK2 + SUB_R + 30
DBL = DBL_ARITH + ROWS_TOTAL + 34 + 3
TRAV = STEP + 24 + 2 + 7
SHARE = ROWS[-1][1]
PRICE_STEP, PRICE_TRAV, PRICE_DBL = Fraction(36803, 64) + SHARE, Fraction(38915, 64) + SHARE, Fraction(45251, 64) + SHARE
check("step is 99.515625 + 444.515625 + 31.015625 + share = 575.046875 + 25,391/2^21", STEP == Fraction(36803, 64) + SHARE,
   "%.6f" % float(STEP))
check("doubling bound is 639 + 31.015625 + share + 34 + 3 = 707.015625 + share",
   DBL_ARITH == 639 and DBL == Fraction(45249, 64) + SHARE)
check("traversal bound is step + 24 + 2 + 7 = 608.046875 + share", TRAV == Fraction(38915, 64) + SHARE)
check("step price equals the table", PRICE_STEP == STEP, "%.6f" % float(PRICE_STEP))
check("traversal price covers its bound", PRICE_TRAV >= TRAV)
check("doubling price 707.046875 + share covers its bound", PRICE_DBL >= DBL)
check("surcharges 33 (traversal) and 132 (doubling) are whole milli-units",
   PRICE_TRAV - PRICE_STEP == 33 and PRICE_DBL - PRICE_STEP == 132)

check("p = 2^256 - C0 with C0 = 2^32 + 977", (1 << 256) - P_FIELD == (1 << 32) + 977)
check("r - 1 is divisible by 6", (r - 1) % 6 == 0)
N = (r - 1) // 6
check("log2 N = 253.41504 to 1e-5", abs(log2(N) - 253.41504) < 1e-5, "%.5f" % log2(N))

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

L_LO, L_HI = ln_interval(1 / (1 - P0))
check("ln(1/(1-p0)) interval is positive and tight", 0 < L_LO < L_HI and L_HI - L_LO < Fraction(1, 2 ** 200),
   "%.15f" % float(L_LO))

N_PRIME = Fraction(N * R_PART, R_PART - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI
m = iroot((V.numerator * 10 ** 80 + V.denominator - 1) // V.denominator, 2) + 1
S_HI = Fraction(m, 10 ** 40)
S_INT = ceil_fraction(S_HI)
check("S_ops <= 2^126.6999 (certified upper bound)", S_HI * S_HI >= V and log2(S_HI) < 126.6999,
   "log2 bound %.9f" % log2(S_HI))
check("S_ops >= 2^126.6997 (lower bound for the quoted value)", log2(S_HI) > 126.6997, "%.9f" % log2(S_HI))

PH1 = ceil_fraction(PRICE_STEP * (S_INT + 2 ** 60))
check("phase 1 = ceil(price*(S_int + 2^60)) covers price*(S_ops + 2^60)",
   Fraction(PH1, 1000) >= PRICE_STEP / 1000 * (S_HI + 2 ** 60))
check("phase 1 < 2^125.90182", log2(Fraction(PH1, 1000)) < 125.90182,
   "%.9f" % log2(Fraction(PH1, 1000)))

CLAIM_LO, CLAIM_HI = pow2_interval(1259019, 10000)
CAP_MILLI_LO = (CLAIM_LO.numerator * 1000) // CLAIM_LO.denominator
BLOCK_MILLI = 2 ** 40 * 1000
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
PERIOD = 256
CP_MILLI = 2 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))
ANCHOR_MILLI = 12 * ceil_fraction(Fraction(S_INT + 2 ** 61, PERIOD))
HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
       ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
       ("phase 2", Fraction(2 ** 20 * (20 * 2 ** 64 + 65) * 600) * 1000),
       ("re-walks", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
       ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
       ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
       ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
       ("table (2^22 entries) and verification, < 2^36", Fraction(2 ** 36) * 1000),
       ("node set-up, 2^26", Fraction(2 ** 26) * 1000),
       ("checkpoint writes: 2*ceil((S_int + 2^61)/256)", Fraction(CP_MILLI)),
       ("anchor rounds: 12*ceil((S_int + 2^61)/256)", Fraction(ANCHOR_MILLI))]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD = sum(int(v) for _, v in HARD_ITEMS)
check("hard extras < 2^112.5783 units", log2(Fraction(HARD, 1000)) < 112.5783,
   "%.9f" % log2(Fraction(HARD, 1000)))
check("the anchor-round row is the largest hard row", ANCHOR_MILLI == max(int(v) for _, v in HARD_ITEMS))

S2 = [[0] * (PERIOD + 1) for _ in range(PERIOD + 1)]
S2[0][0] = 1
for L_ in range(1, PERIOD + 1):
  for b_ in range(1, L_ // 2 + 1):
    S2[L_][b_] = b_ * S2[L_ - 1][b_] + (L_ - 1) * S2[L_ - 2][b_ - 1]
check("S2(3,1), S2(4,1), S2(4,2), S2(6,3) = 1, 1, 3, 15", (S2[3][1], S2[4][1], S2[4][2], S2[6][3]) == (1, 1, 3, 15))


def other_cycles(R, top=PERIOD):
  return sum(Fraction(S2[L_][b_], R ** (L_ - b_)) for L_ in range(3, top + 1) for b_ in range(1, L_ // 2 + 1))


check("at R = 2048 the sum to L = 64 is 2^-19.9954 (the public package's value)",
   abs(log2(other_cycles(2048, 64)) + 19.9954) < 5e-5, "2^%.6f" % log2(other_cycles(2048, 64)))
OTHER = other_cycles(R_PART)
check("at R = 2^22 the sum is < 2^-41.999", OTHER < Fraction(1, 1) / pow2_interval(41999, 1000)[0],
   "2^%.7f" % log2(OTHER))
RATE_HI = Fraction(1, 6 * R_PART) + OTHER
check("doubling rate 1/(6R) + sum <= 2^-24.5849", RATE_HI <= Fraction(1, 1) / pow2_interval(245849, 10000)[0],
   "2^%.9f" % log2(RATE_HI))
TRAV_RATE = Fraction(2, 6 * R_PART) + PERIOD * OTHER
check("traversal-step rate 2/(6R) + 256*sum <= 2^-23.5839", TRAV_RATE <= Fraction(1, 1) / pow2_interval(235839, 10000)[0],
   "2^%.9f" % log2(TRAV_RATE))
OPS = S_HI + 2 ** 60
RATE_K = 3 * Fraction(2 ** 32 + 977, 2 ** 64) + 6 * Fraction(1, 2 ** 63)
check("kernel rare-exit rate 3*C0/2^64 + 6*2^-63 <= 2^-30.415 per operation",
   RATE_K <= pow2_interval(-30415, 1000)[0], "2^%.5f" % log2(RATE_K))
EXTRA_ITEMS = [ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * RATE_HI * OPS * 1000),
       ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * TRAV_RATE * OPS * 1000),
       ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 21) * 1000),
       ceil_fraction(Fraction(2 ** 70) * 1000),
       ceil_fraction(Fraction(20) * (RATE_HI + Fraction(1, 2 ** 56)) * OPS),
       ceil_fraction(Fraction(20) * RATE_K * (S_HI + 2 ** 61))]
E = sum(EXTRA_ITEMS)
check("expected random extras < 2^100.1282 units", log2(Fraction(E, 1000)) < 100.1282,
   "2^%.9f" % log2(Fraction(E, 1000)))

A = -((-(E * 100000)) // 13)
check("allowance A = ceil(E/0.00013) (Markov safe)", A * 13 >= E * 100000)
W_LEDGER = PH1 + HARD + A
check("W_ledger < 2^125.9019", log2(Fraction(W_LEDGER, 1000)) < 125.9019,
   "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("W_ledger > 2^125.90184 (for the quoted value)", log2(Fraction(W_LEDGER, 1000)) > 125.90184,
   "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("the claim is 1259019/10000", CLAIM == Fraction(1259019, 10000))
check("W_ledger <= 1000 * 2^125.9019", W_LEDGER <= CAP_MILLI_LO)
check("125.9019 holds: W^10000 < 1000^10000 * 2^1259019",
   W_LEDGER ** 10000 < 1000 ** 10000 * 2 ** 1259019)
check("125.9018 fails: W^10000 > 1000^10000 * 2^1259018",
   W_LEDGER ** 10000 > 1000 ** 10000 * 2 ** 1259018)

HEAD = CAP_MILLI_LO - PH1 - HARD
check("headroom > allowance A", HEAD > A, "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio < 0.00013", MARKOV_HI < HEADROOM and MARKOV_HI * 100000 < 13,
   "%.10f" % float(MARKOV_HI))
SUCCESS_LO = P0 - HEADROOM - LOSS_OTHER
check("success = 0.3901303 - 0.00013 - 1.1e-7 = 0.39000019 > 0.39",
   SUCCESS_LO == Fraction(3900001900, 10 ** 10) and SUCCESS_LO > THRESHOLD,
   "%.7f" % float(SUCCESS_LO))
check("success floor >= the 0.39 of rules.md", SUCCESS_LO > THRESHOLD, "%.7f" % float(SUCCESS_LO))

C0 = 2 ** 32 + 977
OPS_ALL = S_HI + 2 ** 61
ROW_CARRY = 7 * OPS_ALL * Fraction(C0 * C0 + 3 * C0, P_FIELD)
ROW_NEG = OPS_ALL * Fraction(C0, 2 ** 256)
ROW_HIGH = 2 * OPS_ALL * Fraction(C0, P_FIELD)
ROW_CANON = OPS_ALL * Fraction(C0 * 2 ** 64, P_FIELD - 1)
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

STORE_INIT_UNITS = Fraction(STORE_INIT_MILLI, 1000)
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(331, 10)[1] + 2 ** 33
CAP_PREP = pow2_interval(590, 10)[1]
check("expected preprocessing total < 2^59.0 units", PREP_EXPECTED < CAP_PREP,
   "2^%.5f" % log2(PREP_EXPECTED))
check("preprocessing cap halts with probability < 2^-30",
   Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS) < pow2_interval(-30, 1)[1],
   "2^%.4f" % log2(Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS)))

print()
print("step      %.6f instructions = %.9f units" % (STEP, STEP / 1000.0))
print("doubling  %.6f instructions = %.9f units" % (DBL, DBL / 1000.0))
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
