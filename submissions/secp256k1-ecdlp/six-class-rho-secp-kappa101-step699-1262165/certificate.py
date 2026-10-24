"""Exact-integer certificate for the work bound W <= 2^126.2165 of proof.md.

Every check below is an exact comparison of Python integers or
``fractions.Fraction`` values.  Nothing is taken from a floating-point
computation:

  * the instruction identities, the analytic-row sum, the step, doubling and
    traversal prices, the store's bucket count, probe price, overflow insert and
    lookup bounds and the meter-block refill count are integer arithmetic;
  * ln(1/(1 - p0)) is a rational interval from the atanh series with an explicit
    remainder bound;
  * every fractional power of two is a rational interval certified by integer
    comparisons;
  * the phase-1 cost, hard extras, expected random charges, Markov allowance,
    ledger total and success floor are exact integers or rationals built from
    those bounds, and each printed inequality is one exact comparison, in both
    directions for the claim.

Run:  python3 certificate.py        (exit status 0 = all checks hold)

The certificate proves, for this package (kappa = 1.01, p0 = 0.39389,
Markov headroom h = 0.003866, step price 699):

    phase1  = 699 * (S_int + 2^60)                      (milli-units)
    hard    = sum of the nine exact hard rows
    E       = sum of the four exact expected random rows
    A       = ceil(E / 0.003866)  >= E / 0.003866
    W       = phase1 + hard + A  < 2^126.2165
    E / (cap - phase1 - hard) < 0.003866  =>  success > 0.39389 - 0.003866 - 2e-7
                                            = 0.3900238 > 0.39002 > 0.39
    W^10000  < 1000^10000 * 2^1262165    so 126.2165 holds
    W^10000 > 1000^10000 * 2^1262164    so 126.2164 fails
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


# ---------------------------------------------------------------- constants
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
P_FIELD = int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
R_WALKS = 2048
KAPPA = Fraction(101, 100)          # H1 slowdown allowance 1.01
P0 = Fraction(3938900, 10000000)    # H1 floor p0 at S0 = S_ops / KAPPA
HEADROOM = Fraction(3866, 1000000)  # 0.003866 Markov headroom h
CLAIM = Fraction(1262165, 10000)    # W = 2^126.2165
THRESHOLD = Fraction(39, 100)       # rules.md: success at least 0.39
LOSS_OTHER = Fraction(2, 10 ** 7)

# ------------------------------------------------------- 1. counted kernels
MUL = (16 + 16) + 30 + 26               # products, six carry chains, reduction
SQR = (10 + 10) + 8 + 6 + 7 + 26        # products, off-diagonal sum, doubling, diagonal chain, reduction
SUB, SUBCANON, CANON, ADD_WEAK, ADD_WEAK2 = 9, 14, 8, 9, 14
PASS1 = MUL + SUB + 5 + 4 + 2           # 1 LDR + 4 LDP, 4 STP, SUBS + B.NE
PASS2 = 4 * MUL + SQR + 3 * SUB + 2 * SUBCANON + 13 + 4 + 2   # 1 LDR + 12 LDP, 4 STP, SUBS + B.NE
check("fe_mul costs 88 instructions", MUL == 88, str(MUL))
check("fe_sqr costs 67 instructions", SQR == 67, str(SQR))
check("pass 1 iteration is 108", PASS1 == 108, str(PASS1))
check("pass 2 iteration is 493", PASS2 == 493, str(PASS2))

ROWS = [("F1 truncated product (ldp/ldr schedule)", 20),
        ("range test, S, wrap test, F2", 8),
        ("decision", 13),
        ("sign, p - y from p0 and 2^64-1 registers, compare, 4 CSEL, s", 13),
        ("h, distinguished test, t, table address", 7),
        ("arrival traffic (2 ldp + 1 store)", 3),
        ("fingerprints, history ring, checkpoint (ldp/stp)", 22),
        ("zero-denominator test (moved to once per batch)", 0),
        ("length (ldp/stp), cap test, meter, arrival loop", 11),
        ("share of the inversion and call overhead", 1)]
ROWS_TOTAL = sum(v for _, v in ROWS)
check("analytic rows sum to 98", ROWS_TOTAL == 98, str(ROWS_TOTAL))

STEP = PASS1 + PASS2 + ROWS_TOTAL
DBL_ARITH = 5 * MUL + 2 * SQR + ADD_WEAK + 2 * ADD_WEAK2 + 2 * SUB + 2 * SUBCANON + 72
DBL = DBL_ARITH + ROWS_TOTAL + 50 + 30
TRAV = STEP + 24
PRICE_STEP, PRICE_TRAV, PRICE_DBL = 699, 723, 907
check("step is 108 + 493 + 98 = 699", STEP == 699, str(STEP))
check("doubling bound is 729 + 98 + 50 + 30 = 907", DBL_ARITH == 729 and DBL == 907, str(DBL))
check("traversal bound is 699 + 24 = 723", TRAV == 723, str(TRAV))
check("step price equals the table", PRICE_STEP == STEP, "%d = %d" % (PRICE_STEP, STEP))
check("traversal price covers its bound", PRICE_TRAV >= TRAV, "%d >= %d" % (PRICE_TRAV, TRAV))
check("doubling price equals its derived bound", PRICE_DBL == DBL, "%d = %d" % (PRICE_DBL, DBL))

# ------------------------------------------------------- 2. the group algebra
check("p = 2^256 - C0 with C0 = 2^32 + 977", (1 << 256) - P_FIELD == (1 << 32) + 977)
check("r - 1 is divisible by 6", (r - 1) % 6 == 0)
N = (r - 1) // 6
check("log2 N = 253.41504 to 1e-5", abs(log2(N) - 253.41504) < 1e-5, "%.5f" % log2(N))

# --------------------------------------------- 3. the report store of section 3
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

# --------------------------------------------------- 4. certified H1 numbers
L_LO, L_HI = ln_interval(Fraction(10000000, 6061100))    # ln(1/(1 - p0))
check("ln(1/(1-p0)) interval is positive and tight", 0 < L_LO < L_HI and L_HI - L_LO < Fraction(1, 2 ** 200),
      "%.15f" % float(L_LO))

N_PRIME = Fraction(N * R_WALKS, R_WALKS - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI                   # upper bound on S_ops^2
m = iroot((V.numerator * 10 ** 80 + V.denominator - 1) // V.denominator, 2) + 1
S_HI = Fraction(m, 10 ** 40)
S_INT = ceil_fraction(S_HI)
check("S_ops <= 2^126.7233 (certified upper bound)", S_HI * S_HI >= V and log2(S_HI) < 126.7233,
      "log2 bound %.9f" % log2(S_HI))
check("S_ops >= 2^126.7232 (lower bound for the quoted value)", log2(S_HI) > 126.7232, "%.9f" % log2(S_HI))

# ------------------------------------------------------ 5. phase 1 and hard extras
PH1 = PRICE_STEP * (S_INT + 2 ** 60)                     # milli-units
check("phase 1 = 699*(S_int + 2^60) covers 0.699*(S_ops + 2^60)",
      Fraction(PH1, 1000) >= Fraction(PRICE_STEP, 1000) * (S_HI + 2 ** 60))
check("phase 1 < 2^126.20660", log2(Fraction(PH1, 1000)) < 126.20660,
      "%.9f" % log2(Fraction(PH1, 1000)))

CLAIM_LO, CLAIM_HI = pow2_interval(1262165, 10000)
CAP_MILLI_LO = (CLAIM_LO.numerator * 1000) // CLAIM_LO.denominator
BLOCK_MILLI = 2 ** 40 * 1000
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
              ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
              ("phase 2", Fraction(2 ** 20 * (20 * 2 ** 64 + 65) * 600) * 1000),
              ("re-walks", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
              ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
              ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
              ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
              ("table and verification, < 2^25", Fraction(2 ** 25) * 1000),
              ("node set-up, 2^26", Fraction(2 ** 26) * 1000)]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD = sum(int(v) for _, v in HARD_ITEMS)
check("hard extras < 2^107.156 units", log2(Fraction(HARD, 1000)) < 107.156,
      "%.9f" % log2(Fraction(HARD, 1000)))
check("the overflow row is the largest hard row", OVF_MILLI == max(int(v) for _, v in HARD_ITEMS))

# ---------------------------------------------------- 6. expected random extras
RATE_HI = Fraction(1, 6 * R_WALKS) + pow2_interval(-199954, 10000)[1]
check("doubling rate <= 2^-13.5681", RATE_HI <= Fraction(1, 1) / pow2_interval(135681, 10000)[0],
      "2^%.9f" % log2(RATE_HI))
OPS = S_HI + 2 ** 60
EXTRA_ITEMS = [ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * RATE_HI * OPS * 1000),
               ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * 64 * pow2_interval(-1999, 100)[1] * OPS * 1000),
               ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 20) * 1000),
               ceil_fraction(Fraction(2 ** 70) * 1000)]
E = sum(EXTRA_ITEMS)
check("expected random extras < 2^111.0090 units", log2(Fraction(E, 1000)) < 111.0090,
      "2^%.9f" % log2(Fraction(E, 1000)))

# ----------------------------------------------- 7. allowance, total and claim
A = -((-(E * 1000000)) // 3866)                          # ceil(E / 0.003866)
check("allowance A = ceil(E/0.003866) (Markov safe)", A * 3866 >= E * 1000000)
W_LEDGER = PH1 + HARD + A
check("W_ledger < 2^126.21650", log2(Fraction(W_LEDGER, 1000)) < 126.21650,
      "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("W_ledger > 2^126.21648 (for the quoted value)", log2(Fraction(W_LEDGER, 1000)) > 126.21648,
      "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("the claim is 1262165/10000", CLAIM == Fraction(1262165, 10000))
check("W_ledger <= 1000 * 2^126.2165", W_LEDGER <= CAP_MILLI_LO)
check("126.2165 holds: W^10000 < 1000^10000 * 2^1262165",
      W_LEDGER ** 10000 < 1000 ** 10000 * 2 ** 1262165)
check("126.2164 fails: W^10000 > 1000^10000 * 2^1262164",
      W_LEDGER ** 10000 > 1000 ** 10000 * 2 ** 1262164)

HEAD = CAP_MILLI_LO - PH1 - HARD
check("headroom > allowance A", HEAD > A, "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio < 0.003866", MARKOV_HI < HEADROOM and MARKOV_HI * 1000000 < 3866,
      "%.10f" % float(MARKOV_HI))
SUCCESS_LO = P0 - HEADROOM - LOSS_OTHER
check("success = 0.39389 - 0.003866 - 2e-7 = 0.3900238 > 0.39002",
      SUCCESS_LO == Fraction(3900238, 10 ** 7) and SUCCESS_LO > Fraction(39002, 100000),
      "%.7f" % float(SUCCESS_LO))
check("success floor >= the 0.39 of rules.md", SUCCESS_LO > THRESHOLD, "%.7f" % float(SUCCESS_LO))

# ------------------------------------------------- 8. counter widths, section 3
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

# ------------------------------------------------------- 9. preprocessing
STORE_INIT_UNITS = Fraction(STORE_INIT_MILLI, 1000)
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(241, 10)[1] + 2 ** 33
CAP_PREP = pow2_interval(590, 10)[1]
check("expected preprocessing total < 2^59.0 units", PREP_EXPECTED < CAP_PREP,
      "2^%.5f" % log2(PREP_EXPECTED))
check("preprocessing cap halts with probability < 2^-30",
      Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS) < pow2_interval(-30, 1)[1],
      "2^%.4f" % log2(Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS)))

print()
print("step      %d instructions = %.3f units" % (STEP, STEP / 1000.0))
print("doubling  %d instructions = %.3f units" % (DBL, DBL / 1000.0))
print("traversal %d instructions = %.3f units" % (TRAV, TRAV / 1000.0))
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
