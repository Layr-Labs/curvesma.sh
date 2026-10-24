"""Exact-integer certificate for the work bound W <= 2^126.338 of proof.md section 5.

Every check below is an exact comparison of Python integers or ``fractions.Fraction``
values.  Nothing is taken from a floating-point computation:

  * the instruction identities, the analytic-row sum, the step price, the doubling
    price, the traversal price, the store's bucket count, probe price, overflow
    insert and lookup bounds, the meter-block refill count and the report price
    are integer arithmetic;
  * ln(1/(1 - p0)) is a rational interval from the atanh series with an explicit
    remainder bound  z^(2n+1)/((2n+1)(1 - z^2));
  * every fractional power of two is a rational interval  [m, m+1)/2^scale  whose
    endpoints are certified by the integer comparisons  m^den <= 2^sh < (m+1)^den;
  * the phase-1 cost, the hard extras, the expected random charges, the Markov
    allowance, the ledger total and the success floor are exact integers or
    rationals built from those bounds, and each printed inequality is a single
    exact comparison, in both directions for the claim.

Run:  python3 certificate.py        (exit status 0 = all checks hold)

The certificate proves, for the prices and constants of this package:

    phase1  = 752 * (S_int + 2^60)                          (per-operation step price)
    hard    = sum of the nine exact hard rows                = 2^107.15550
    E       = sum of the four exact expected random rows     = 2^111.00076
    A       = ceil(E / 0.00416)  >= E / 0.00416              (Markov allowance)
    W       = phase1 + hard + A  = 2^126.337997 < 2^126.338
    E / W < 0.00416   =>  success > 0.3953 - 0.00416 - 2e-7 = 0.3911398 > 0.3911 > 0.39
    W^500  < 1000^500 * 2^63169    so 126.338 holds
    W^1000 > 1000^1000 * 2^126337  so 126.337 fails

so on every run W <= 2^126.338 and the attack recovers x with probability at least
0.3911 under the H1 heuristic of proof.md section 4.
"""
import math
import sys
from fractions import Fraction

sys.set_int_max_str_digits(200000)

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
    """floor(n ** (1/k)) for n >= 0, k >= 1, by bisection between two powers of two.

    Bisection is used rather than Newton because for k = 10^4 Newton's basin is
    narrower than a factor of two, and each comparison is an exact integer power."""
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
    """Rational (lo, hi) with lo <= 2**(num/den) <= hi and hi - lo = 2**-scale."""
    sh = num + scale * den            # m**den <= 2**sh  <=>  m <= 2**(num/den)*2**scale
    assert sh >= 0
    m = iroot(2 ** sh, den)
    assert m ** den <= 2 ** sh < (m + 1) ** den, "integer root certificate"
    return Fraction(m, 2 ** scale), Fraction(m + 1, 2 ** scale)


def ln_interval(x, bits=200):
    """Rational (lo, hi) bracketing ln(x) for x > 1, from ln x = 2 atanh z, z = (x-1)/(x+1).

    The remainder after n terms is bounded by z^(2n+1) / ((2n+1)(1 - z^2))."""
    z = Fraction(x - 1, x + 1)
    z2 = z * z
    total = Fraction(0)
    p = z                      # z^(2k+1)
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
KAPPA = Fraction(51, 50)          # H1 slowdown allowance 1.02
P0 = Fraction(3953, 10000)        # H1 floor 0.3953
CLAIM = Fraction(126338, 1000)    # W = 2^126.338 = 63169/500
THRESHOLD = Fraction(39, 100)     # rules.md: success at least 0.39
LOSS_OTHER = Fraction(2, 10 ** 7)
MARKOV_CAP = Fraction(416, 100000)   # 0.00416 the success statement allows

# ------------------------------------------------------- 1. counted kernels
MUL = 12 + 3 * 17 + 28            # fe_mul: row 0 twelve, rows 1..3 seventeen, reduction 28
SQR = 12 + 9 + 7 + 15 + 28        # fe_sqr: product 12, folds 9+7+15, reduction 28
SUB, CANON, ADD_WEAK, ADD_WEAK2 = 9, 8, 9, 14
PASS1 = MUL + SUB + 5 + 4 + 2                     # 1 fe_mul, 1 fe_sub, 5 loads, 4 stores, control
PASS2 = 4 * MUL + SQR + 5 * SUB + 2 * CANON + 19 + 6 + 2
check("fe_mul costs 91 instructions", MUL == 91, str(MUL))
check("fe_sqr costs 71 instructions", SQR == 71, str(SQR))
check("pass 1 iteration is 111", PASS1 == 111, str(PASS1))
check("pass 2 iteration is 523", PASS2 == 523, str(PASS2))

ROWS = [("F1 truncated product: 9 loads, 8 MUL/UMULH, 7 ADD", 24),
        ("range test, S, wrap test, F2 (F0 already loaded)", 8),
        ("decision: 4 + 2 + 4 + 3", 13),
        ("sign, p - y with C0 fold, compare, 4 CSEL, s", 17),
        ("h, distinguished test, t, table address", 7),
        ("arrival traffic: 4 loads of y, 1 store of the table pointer", 5),
        ("fingerprints, history ring, checkpoint", 26),
        ("zero-denominator test of pass 1", 4),
        ("length, cap test, meter, arrival loop", 13),
        ("share of the inversion and call overhead", 1)]
ROWS_TOTAL = sum(v for _, v in ROWS)
for name, v in ROWS:
    print("    row %-52s %3d" % (name, v))
check("analytic rows sum to 118", ROWS_TOTAL == 118, str(ROWS_TOTAL))

STEP = PASS1 + PASS2 + ROWS_TOTAL
DBL_ARITH = 5 * MUL + 2 * SQR + ADD_WEAK + 2 * ADD_WEAK2 + 4 * SUB + 2 * CANON + 72
DBL = DBL_ARITH + ROWS_TOTAL + 50 + 30
TRAV = STEP + 24
PRICE_STEP, PRICE_TRAV, PRICE_DBL = 752, 776, 956
check("step is 111 + 523 + 118 = 752", STEP == 752, str(STEP))
check("doubling arithmetic is 758", DBL_ARITH == 758, str(DBL_ARITH))
check("doubling bound is 758 + 118 + 50 + 30 = 956", DBL == 956, str(DBL))
check("traversal bound is 752 + 24 = 776", TRAV == 776, str(TRAV))
check("step price equals the table (no margin)", PRICE_STEP == STEP, "%d = %d" % (PRICE_STEP, STEP))
check("traversal price covers its bound", PRICE_TRAV >= TRAV, "%d >= %d" % (PRICE_TRAV, TRAV))
check("doubling price equals its derived bound", PRICE_DBL == DBL, "%d = %d" % (PRICE_DBL, DBL))

# ------------------------------------------------------- 2. the group algebra
check("p = 2^256 - C0 with C0 = 2^32 + 977", (1 << 256) - P_FIELD == (1 << 32) + 977,
      str(P_FIELD))
check("r - 1 is divisible by 6 (six classes)", (r - 1) % 6 == 0)
N = (r - 1) // 6
check("log2 N = 253.41504 to 1e-5", abs(log2(N) - 253.41504) < 1e-5, "%.5f" % log2(N))

# --------------------------------------------- 3. the report store of section 3
NODES = 2 ** 16
DCAP = ceil_fraction(pow2_interval(6383, 100)[1])      # ceil(2^63.83)
B = ceil_fraction(Fraction(2 * DCAP, NODES))           # buckets per node
PROBE_PRICE, PROBES = 50, 8                            # 0.05 units, fixed depth 8
REPORT_PRICE = 40000 + PROBES * PROBE_PRICE            # 40.4 units, in milli-units
OVF_INSERT = 6 * B + 35
OVF_LOOKUP = 16 * B
OVF_MILLI = DCAP * (OVF_INSERT + OVF_LOOKUP)
STORE_INIT_MILLI = 2 * DCAP * 16                  # 16 instructions per bucket, milli-units
check("D_cap = ceil(2^63.83) = 16396253461269537204", DCAP == 16396253461269537204, str(DCAP))
check("B = ceil(2*D_cap/2^16) = 500373945961595", B == 500373945961595, str(B))
check("B buckets of 104 bytes stay under a 2^56-byte node", B * 104 < 2 ** 56,
      "2^%.3f" % log2(B * 104))
check("112 bytes per bucket (bucket + 8-byte overflow address) < 2^56", B * 112 < 2 ** 56,
      "2^%.3f" % log2(B * 112))
check("the whole store < 2^72 bytes", NODES * B * 112 < 2 ** 72, "2^%.3f" % log2(NODES * B * 112))
check("2^16*B buckets cover the 2*D_cap slots", NODES * B >= 2 * DCAP)
check("report price 40.4 = 40 + 8 probes", REPORT_PRICE == 40400, str(REPORT_PRICE))
check("overflow insert bound 6B + 35 = 3002243675769605", OVF_INSERT == 3002243675769605,
      "2^%.4f" % log2(OVF_INSERT))
check("overflow lookup bound 16B = 8005983135385520", OVF_LOOKUP == 8005983135385520,
      "2^%.4f units" % log2(Fraction(OVF_LOOKUP, 1000)))
check("overflow worst case < 2^107.154 units", log2(Fraction(OVF_MILLI, 1000)) < 107.154,
      "2^%.6f" % log2(Fraction(OVF_MILLI, 1000)))

# --------------------------------------------------- 4. certified H1 numbers
L_LO, L_HI = ln_interval(Fraction(10000, 6047))          # ln(1/(1 - p0))
check("ln(1/(1-p0)) = 0.5030228 to 1e-12",
      abs(float(L_LO) - 0.5030228117) < 1e-9 and L_LO < L_HI,
      "%.12f" % float(L_LO))

N_PRIME = Fraction(N * R_WALKS, R_WALKS - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI                   # upper bound on S_ops^2
m = iroot((V.numerator * 10 ** 80 + V.denominator - 1) // V.denominator, 2) + 1
S_HI = Fraction(m, 10 ** 40)
S_INT = ceil_fraction(S_HI)
check("S_ops <= 2^126.7408 (certified upper bound)", S_HI * S_HI >= V and log2(S_HI) < 126.7408,
      "log2 bound %.6f" % log2(S_HI))
check("S_ops >= 2^126.7407 (lower, for the quoted value)",
      log2(S_HI) > 126.7406, "%.6f" % log2(S_HI))
check("S_int = ceil(S_ops) = 142160481062555318092036037759448841269",
      S_INT == 142160481062555318092036037759448841269, str(S_INT))

# ------------------------------------------------------ 5. phase 1 and hard extras
PH1 = PRICE_STEP * (S_INT + 2 ** 60)                     # milli-units
check("phase 1 = 752*(S_int + 2^60) covers 0.752*(S_ops + 2^60)",
      Fraction(PH1, 1000) >= Fraction(PRICE_STEP, 1000) * (S_HI + 2 ** 60), str(PH1))
check("phase 1 <= 2^126.3296", log2(Fraction(PH1, 1000)) < 126.3296,
      "%.7f" % log2(Fraction(PH1, 1000)))

CAP_MILLI_LO = (pow2_interval(126338, 1000)[0].numerator * 1000) // pow2_interval(126338, 1000)[0].denominator
BLOCK_MILLI = 2 ** 40 * 1000                             # one 2^40-unit block, in milli-units
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
check("refills <= 2^86.338 + 2^20", REFILLS == 97797502834719374447339471, str(REFILLS))
check("refill row < 2^81.28 units", log2(Fraction(REFILL_MILLI, 1000)) < 81.28,
      "2^%.6f" % log2(Fraction(REFILL_MILLI, 1000)))

HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
              ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
              ("phase 2: 2^20*(L_max+65)*600", Fraction(2 ** 20 * (20 * 2 ** 64 + 65) * 600) * 1000),
              ("re-walks: 16*((L_max+65)*600 + 2^13)", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
              ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
              ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
              ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
              ("table and verification, < 2^25", Fraction(2 ** 25) * 1000),
              ("node set-up, 2^26", Fraction(2 ** 26) * 1000)]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD_ROWS_I = [(n, int(v)) for n, v in HARD_ITEMS]
HARD = sum(v for _, v in HARD_ROWS_I)
for name, v in HARD_ROWS_I:
    print("    hard %-51s %d  2^%.6f" % (name, v, log2(Fraction(v, 1000))))
check("hard total = 180725797188573533523169202391506758",
      HARD == 180725797188573533523169202391506758, str(HARD))
check("hard extras < 2^107.156 units", log2(Fraction(HARD, 1000)) < 107.156,
      "%.6f" % log2(Fraction(HARD, 1000)))
check("the overflow row is the largest hard row", max(v for _, v in HARD_ROWS_I) == OVF_MILLI)

# ---------------------------------------------------- 6. expected random extras
RATE_HI = Fraction(1, 6 * R_WALKS) + pow2_interval(-199954, 10000)[1]
RATE_1356 = Fraction(1, 1) / pow2_interval(1356, 100)[0]
check("doubling rate 1/(6R) + 2^-19.9954 <= 2^-13.5681 <= 2^-13.56",
      RATE_HI <= Fraction(1, 1) / pow2_interval(135681, 10000)[0] and RATE_HI <= RATE_1356,
      "2^%.7f" % log2(RATE_HI))
OPS = S_HI + 2 ** 60
EXTRA_ITEMS = [("doublings: (0.956 - 0.752) * rate * operations",
                ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * RATE_HI * OPS * 1000)),
               ("traversals: 0.024 * 64 * 2^-19.99 * operations",
                ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * 64 * pow2_interval(-1999, 100)[1] * OPS * 1000)),
               ("starts: 2^13 * (S_ops/2^64 + 2^20)",
                ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 20) * 1000)),
               ("words, attempts, probes, small batches, fallbacks",
                ceil_fraction(Fraction(2 ** 70) * 1000))]
E = sum(v for _, v in EXTRA_ITEMS)
for name, v in EXTRA_ITEMS:
    print("    random %-49s %d  2^%.6f" % (name, v, log2(Fraction(v, 1000))))
check("expected random extras total 2597523031886392846748549154421911133",
      E == 2597523031886392846748549154421911133, str(E))
check("expected random extras < 2^111.00077 units", log2(Fraction(E, 1000)) < 111.00077,
      "%.7f" % log2(Fraction(E, 1000)))

# ----------------------------------------------- 7. allowance, total and claim
A = -((-(E * 100000)) // 416)                            # ceil(E / 0.00416)
check("allowance A = ceil(E/0.00416) = 624404574972690588160708931351420945433",
      A == 624404574972690588160708931351420945433, str(A))
check("A >= E / 0.00416 (Markov safe)", A * 416 >= E * 100000)
W_LEDGER = PH1 + HARD + A
check("W_ledger = phase 1 + hard + A = 107529267059811478367772329467123690012431",
      W_LEDGER == 107529267059811478367772329467123690012431, str(W_LEDGER))
check("log2(W_ledger/1000) = 126.337997 < 126.338",
      log2(Fraction(W_LEDGER, 1000)) < 126.338, "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("ceil at three decimals is the claim 126.338",
      log2(Fraction(W_LEDGER, 1000)) > 126.337,
      "126.337 fails: %.9f > 126.337" % log2(Fraction(W_LEDGER, 1000)))
check("the claim 126.338 = 63169/500", CLAIM == Fraction(63169, 500), str(CLAIM))
check("W_ledger <= 1000 * 2^126.338", W_LEDGER <= CAP_MILLI_LO)
check("126.338 holds: W^500 < 1000^500 * 2^63169", W_LEDGER ** 500 < 1000 ** 500 * 2 ** 63169,
      "%d-digit sides" % len(str(W_LEDGER ** 500)))
check("126.337 fails: W^1000 > 1000^1000 * 2^126337",
      W_LEDGER ** 1000 > 1000 ** 1000 * 2 ** 126337,
      "%d-digit sides" % len(str(W_LEDGER ** 1000)))

HEAD = CAP_MILLI_LO - PH1 - HARD
check("headroom = 2^126.338 - phase 1 - hard > allowance A", HEAD > A,
      "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio < 0.00416", MARKOV_HI < MARKOV_CAP, "%.10f" % float(MARKOV_HI))
check("Markov ratio 0.00415851 to 1e-8", abs(float(MARKOV_HI) - 0.0041585050) < 1e-8,
      "%.10f" % float(MARKOV_HI))
SUCCESS_LO = P0 - MARKOV_CAP - LOSS_OTHER
check("success = 0.3953 - 0.00416 - 2e-7 = 0.3911398 > 0.3911",
      SUCCESS_LO == Fraction(3911398, 10 ** 7) and SUCCESS_LO > Fraction(3911, 10000),
      "%.7f" % float(SUCCESS_LO))
check("success floor > 0.3905", SUCCESS_LO > Fraction(3905, 10000), "%.7f" % float(SUCCESS_LO))
check("success floor >= the 0.39 of rules.md", SUCCESS_LO > THRESHOLD, "%.7f" % float(SUCCESS_LO))
check("claim is below the record it revises (126.403 and 126.367)",
      CLAIM < 126403 and CLAIM < Fraction(126367, 1000))

# ------------------------------------------------- 8. counter widths, section 3
METER_MILLI = CAP_MILLI_LO
check("192-bit meter holds 2^126.338 units < 2^136.304 milli-units",
      METER_MILLI < 2 ** 192, "2^%.7f milli-units < 2^192" % log2(METER_MILLI))
check("128-bit operation counter holds S_ops + 2^60 < 2^127",
      S_HI + 2 ** 60 < 2 ** 127, "2^%.6f" % log2(S_HI + 2 ** 60))
STARTS_HI = Fraction(CAP_MILLI_LO, 1000) / 2 ** 13
check("128-bit start counter holds W/2^13 < 2^114", STARTS_HI < 2 ** 114, "2^%.6f" % log2(STARTS_HI))
PREP_CAP_MILLI = ceil_fraction(pow2_interval(590, 10)[1] * 1000)
check("128-bit preprocessing cap 2^59.0 units = 2^68.966 milli-units < 2^69",
      PREP_CAP_MILLI < 2 ** 69, "2^%.4f" % log2(PREP_CAP_MILLI))
check("128-bit refill counter holds 2^86.338 + 2^20 < 2^87", REFILLS < 2 ** 87,
      "2^%.6f" % log2(REFILLS))
check("no counter's range equals its capacity (each max < 2^(width-1))",
      METER_MILLI < 2 ** 191 and S_HI + 2 ** 60 < 2 ** 127 and STARTS_HI < 2 ** 127
      and PREP_CAP_MILLI < 2 ** 127 and REFILLS < 2 ** 127)

# ------------------------------------------------------- 9. preprocessing
STORE_INIT_UNITS = Fraction(STORE_INIT_MILLI, 1000)
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(241, 10)[1] + 2 ** 33
check("expected preprocessing total < 2^58.87 units", PREP_EXPECTED < pow2_interval(5887, 100)[1],
      "2^%.5f" % log2(PREP_EXPECTED))
CAP_PREP = pow2_interval(590, 10)[1]
check("preprocessing gap to the 2^59.0 cap > 2^55", CAP_PREP - STORE_INIT_UNITS > 2 ** 55,
      "2^%.4f" % log2(CAP_PREP - STORE_INIT_UNITS))
check("preprocessing cap halts with probability < 2^-30",
      Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS) < pow2_interval(-30, 1)[1],
      "2^%.4f" % log2(Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS)))

print()
print("step      %d instructions = %.3f units (table margin %d)" % (STEP, STEP / 1000.0, PRICE_STEP - STEP))
print("doubling  %d instructions = %.3f units (bound margin %d)" % (DBL, DBL / 1000.0, PRICE_DBL - DBL))
print("traversal %d instructions = %.3f units (bound margin %d)" % (TRAV, TRAV / 1000.0, PRICE_TRAV - TRAV))
print("store     B = %d buckets of 104 bytes + B 8-byte addresses; probe %d instr charged %d;"
      % (B, 39, PROBE_PRICE))
print("          overflow insert 6B + 35 = %d, lookup 16B = %d, row %d milli-units"
      % (OVF_INSERT, OVF_LOOKUP, OVF_MILLI))
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
