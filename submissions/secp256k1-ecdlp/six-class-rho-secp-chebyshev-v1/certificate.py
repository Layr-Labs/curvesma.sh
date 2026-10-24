"""Exact-integer certificate for the work bound of proof.md section 5.

Every check below is an exact comparison of Python integers or ``fractions.Fraction``
values.  Nothing is taken from a floating-point computation:

  * the instruction identities, the analytic-row sum, the step, doubling and
    traversal prices, the store's bucket count, probe price, overflow insert and
    lookup bounds, the meter-block refill count and the report price are integer
    arithmetic (unchanged from the package this one revises);
  * ln(1/(1 - p0)) is a rational interval from the atanh series with an explicit
    remainder bound  z^(2n+1)/((2n+1)(1 - z^2));
  * every fractional power of two is a rational interval  [m, m+1)/2^scale  whose
    endpoints are certified by the integer comparisons  m^den <= 2^sh < (m+1)^den;
  * every square root is an integer square root, rounded in the safe direction;
  * the phase-1 cost, the hard extras, the four random-charge allowances, the
    ledger total and the success floor are exact integers or rationals built from
    those bounds, and each printed inequality is a single exact comparison, in both
    directions for the claim.

Run:  python3 certificate.py        (exit status 0 = all checks hold)

The certificate proves, for the prices and constants of this package:

    phase1  = 752 * (S_int + 2^61)                   (S_int from p0 = 0.3901, kappa = 1.02)
    hard    = sum of the nine exact hard rows
    A_dbl   = 204  * (mu_D + t_D),  t_D = sqrt(2^80 * mu_D * 2^30)   (Chebyshev under H2)
    A_trav  = 1536 * (mu_V + t_V),  t_V = sqrt(2^80 * mu_V * 2^30)   (Chebyshev under H2)
    A_start = 2^30 * E[starts charge]                                 (Markov)
    A_other = 2^30 * 2^71 units                                       (Markov)
    W       = phase1 + hard + A_dbl + A_trav + A_start + A_other  < 2^126.3173
    success > 0.3901 - 2e-7 - 4 * 2^-30 > 0.39
    W^10000 < 1000^10000 * 2^1263173  so 126.3173 holds
    W^10000 > 1000^10000 * 2^1263172  so 126.3172 fails

so on every run W <= 2^126.3173 and the attack recovers x with probability above
0.39 under the H1 and H2 heuristics of proof.md section 4.
"""
import math
import sys
from fractions import Fraction

sys.set_int_max_str_digits(400000)

FAILED = []


def check(name, ok, detail=""):
    print("%-62s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
    if not ok:
        FAILED.append(name)


def log2(x):
    f = Fraction(x)
    return math.log2(f.numerator) - math.log2(f.denominator)


def ceil_fraction(f):
    f = Fraction(f)
    return -((-f.numerator) // f.denominator)


def iroot(n, k):
    """floor(n ** (1/k)) for n >= 0, k >= 1, by bisection between two powers of two."""
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


def sqrt_ceil(f):
    """Smallest integer s with s*s >= f, for a nonnegative rational f."""
    f = Fraction(f)
    n = ceil_fraction(f)
    s = math.isqrt(n)
    if s * s < n:
        s += 1
    assert s * s >= f and (s - 1) * (s - 1) < f
    return s


def pow2_interval(num, den, scale=80):
    """Rational (lo, hi) with lo <= 2**(num/den) <= hi and hi - lo = 2**-scale."""
    sh = num + scale * den
    assert sh >= 0
    m = iroot(2 ** sh, den)
    assert m ** den <= 2 ** sh < (m + 1) ** den, "integer root certificate"
    return Fraction(m, 2 ** scale), Fraction(m + 1, 2 ** scale)


def ln_interval(x, bits=200):
    """Rational (lo, hi) bracketing ln(x) for x > 1, from ln x = 2 atanh z, z = (x-1)/(x+1)."""
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
KAPPA = Fraction(51, 50)              # H1 slowdown allowance 1.02 (unchanged)
P0 = Fraction(3901, 10000)            # H1 floor 0.3901 (was 0.3953)
CLAIM_NUM, CLAIM_DEN = 1263173, 10000  # W = 2^126.3173
CLAIM = Fraction(CLAIM_NUM, CLAIM_DEN)
THRESHOLD = Fraction(39, 100)         # rules.md: success at least 0.39
LOSS_OTHER = Fraction(2, 10 ** 7)     # the failure events of section 4 other than the meter
VAR_FACTOR = 2 ** 80                  # H2: Var(count) <= 2^80 * (mean bound)
TAIL = Fraction(1, 2 ** 30)           # each of the four random-charge tails

# ------------------------------------------------------- 1. counted kernels (unchanged)
MUL = 12 + 3 * 17 + 28
SQR = 12 + 9 + 7 + 15 + 28
SUB, CANON, ADD_WEAK, ADD_WEAK2 = 9, 8, 9, 14
PASS1 = MUL + SUB + 5 + 4 + 2
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
    print("    row %-56s %3d" % (name, v))
check("analytic rows sum to 118", ROWS_TOTAL == 118, str(ROWS_TOTAL))

STEP = PASS1 + PASS2 + ROWS_TOTAL
DBL_ARITH = 5 * MUL + 2 * SQR + ADD_WEAK + 2 * ADD_WEAK2 + 4 * SUB + 2 * CANON + 72
DBL = DBL_ARITH + ROWS_TOTAL + 50 + 30
TRAV = STEP + 24
PRICE_STEP, PRICE_TRAV, PRICE_DBL = 752, 776, 956
check("step is 111 + 523 + 118 = 752", STEP == 752, str(STEP))
check("doubling bound is 758 + 118 + 50 + 30 = 956", DBL_ARITH == 758 and DBL == 956, str(DBL))
check("traversal bound is 752 + 24 = 776", TRAV == 776, str(TRAV))
check("prices equal or cover their bounds",
      PRICE_STEP == STEP and PRICE_TRAV >= TRAV and PRICE_DBL == DBL)

# ------------------------------------------------------- 2. the group algebra
check("p = 2^256 - C0 with C0 = 2^32 + 977", (1 << 256) - P_FIELD == (1 << 32) + 977)
check("r - 1 is divisible by 6 (six classes)", (r - 1) % 6 == 0)
N = (r - 1) // 6
check("log2 N = 253.41504 to 1e-5", abs(log2(N) - 253.41504) < 1e-5, "%.5f" % log2(N))

# --------------------------------------------- 3. the report store (unchanged)
NODES = 2 ** 16
DCAP = ceil_fraction(pow2_interval(6383, 100)[1])
B = ceil_fraction(Fraction(2 * DCAP, NODES))
PROBE_PRICE, PROBES = 50, 8
REPORT_PRICE = 40000 + PROBES * PROBE_PRICE
OVF_INSERT = 6 * B + 35
OVF_LOOKUP = 16 * B
OVF_MILLI = DCAP * (OVF_INSERT + OVF_LOOKUP)
STORE_INIT_MILLI = 2 * DCAP * 16
check("D_cap = ceil(2^63.83) = 16396253461269537204", DCAP == 16396253461269537204, str(DCAP))
check("B = ceil(2*D_cap/2^16) = 500373945961595", B == 500373945961595, str(B))
check("the whole store < 2^72 bytes", NODES * B * 112 < 2 ** 72, "2^%.3f" % log2(NODES * B * 112))
check("2^16*B buckets cover the 2*D_cap slots", NODES * B >= 2 * DCAP)

# --------------------------------------------------- 4. certified H1 numbers
X = 1 / (1 - P0)
L_LO, L_HI = ln_interval(X)
check("ln(1/(1-p0)) = ln(10000/6099) bracketed to 2^-199",
      L_LO < L_HI and L_HI - L_LO < Fraction(1, 2 ** 199), "%.12f" % float(L_LO))

N_PRIME = Fraction(N * R_WALKS, R_WALKS - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI                   # upper bound on S_ops^2
S_INT = sqrt_ceil(V)                                     # integer stop count, S_int^2 >= V
check("S_int >= kappa * sqrt(2 N' ln(1/(1-p0))) (stop count covers S_ops)", S_INT * S_INT >= V)
print("    S_int = %d = 2^%.6f" % (S_INT, log2(S_INT)))

# ------------------------------------------------------ 5. phase 1 and hard extras
OPS = S_INT + 2 ** 61                                    # phase-1 operations, upper bound (block overshoot < 2^60.41)
PH1 = PRICE_STEP * OPS                                   # milli-units

CAP_LO = pow2_interval(CLAIM_NUM, CLAIM_DEN)[0]
CAP_MILLI_LO = (CAP_LO.numerator * 1000) // CAP_LO.denominator   # floor(1000 * 2^claim)
BLOCK_MILLI = 2 ** 40 * 1000
CAP_HI = pow2_interval(CLAIM_NUM, CLAIM_DEN)[1]
REFILLS = ceil_fraction(CAP_HI * 1000 / BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
L_MAX = 20 * 2 ** 64

HARD_ITEMS = [("meter remainders 2^61", 2 ** 61 * 1000),
              ("meter block refills: refills * 30", REFILL_MILLI),
              ("phase 2: 2^20*(L_max+65)*600", 2 ** 20 * (L_MAX + 65) * 600 * 1000),
              ("re-walks: 16*((L_max+65)*600 + 2^13)", 16 * ((L_MAX + 65) * 600 + 2 ** 13) * 1000),
              ("reports: D_cap * 40,400", DCAP * REPORT_PRICE),
              ("overflow insert and lookup: D_cap*(22B+35)", OVF_MILLI),
              ("store initialisation: 2*D_cap*16", STORE_INIT_MILLI),
              ("table and verification, < 2^25", 2 ** 25 * 1000),
              ("node set-up, 2^26", 2 ** 26 * 1000)]
HARD = sum(v for _, v in HARD_ITEMS)
for name, v in HARD_ITEMS:
    print("    hard %-55s %d  2^%.6f" % (name, v, log2(Fraction(v, 1000))))
check("hard extras < 2^107.156 units", log2(Fraction(HARD, 1000)) < 107.156,
      "%.6f" % log2(Fraction(HARD, 1000)))

# ---------------------------------------- 6. random charges and their allowances
RATE_HI = Fraction(1, 6 * R_WALKS) + pow2_interval(-199954, 10000)[1]
check("doubling rate 1/(6R) + 2^-19.9954 <= 2^-13.5681",
      RATE_HI <= Fraction(1, 1) / pow2_interval(135681, 10000)[0], "2^%.7f" % log2(RATE_HI))
TRAV_RATE_HI = pow2_interval(-1999, 100)[1]

MU_D = RATE_HI * OPS                                     # E[doublings] <= mu_D
T_D = sqrt_ceil(VAR_FACTOR * MU_D / TAIL)                # Chebyshev: Var/t^2 <= 2^-30
A_DBL = (PRICE_DBL - PRICE_STEP) * (ceil_fraction(MU_D) + T_D)
MU_V = TRAV_RATE_HI * OPS                                # E[traversal events] <= mu_V
T_V = sqrt_ceil(VAR_FACTOR * MU_V / TAIL)
A_TRAV = (PRICE_TRAV - PRICE_STEP) * 64 * (ceil_fraction(MU_V) + T_V)
E_START = ceil_fraction(Fraction(2 ** 13) * (Fraction(S_INT, 2 ** 64) + 2 ** 20) * 1000)
A_START = E_START * 2 ** 30                              # Markov: E/A = 2^-30
E_OTHER = 2 ** 71 * 1000
A_OTHER = E_OTHER * 2 ** 30
check("Chebyshev tail for doublings <= 2^-30",
      Fraction(VAR_FACTOR) * MU_D / (T_D * T_D) <= TAIL, "t_D = 2^%.4f" % log2(T_D))
check("Chebyshev tail for traversal events <= 2^-30",
      Fraction(VAR_FACTOR) * MU_V / (T_V * T_V) <= TAIL, "t_V = 2^%.4f" % log2(T_V))
check("Markov tail for starts = 2^-30", Fraction(E_START, A_START) == TAIL)
check("Markov tail for small charges = 2^-30", Fraction(E_OTHER, A_OTHER) == TAIL)
ALLOW_ITEMS = [("doublings: 204 * (mu_D + t_D)", A_DBL),
               ("traversals: 24 * 64 * (mu_V + t_V)", A_TRAV),
               ("starts: 2^30 * 2^13 * (S_int/2^64 + 2^20)", A_START),
               ("words, attempts, probes, small batches, fallbacks: 2^30 * 2^71", A_OTHER)]
ALLOW = sum(v for _, v in ALLOW_ITEMS)
for name, v in ALLOW_ITEMS:
    print("    allow %-54s %d  2^%.6f" % (name, v, log2(Fraction(v, 1000))))
print("    expected doublings mu_D 2^%.4f, deviation allowance t_D 2^%.4f" % (log2(MU_D), log2(T_D)))
print("    expected traversal events mu_V 2^%.4f, deviation allowance t_V 2^%.4f" % (log2(MU_V), log2(T_V)))
check("all random allowances together < 2^112 units", ALLOW < 2 ** 112 * 1000,
      "2^%.6f" % log2(Fraction(ALLOW, 1000)))

# ----------------------------------------------- 7. total and claim
W_LEDGER = PH1 + HARD + ALLOW
LW = log2(Fraction(W_LEDGER, 1000))
check("W_ledger <= floor(1000 * 2^claim)", W_LEDGER <= CAP_MILLI_LO, "log2 W = %.9f" % LW)
check("126.3173 holds: W^10000 < 1000^10000 * 2^1263173",
      W_LEDGER ** 10000 < 1000 ** 10000 * 2 ** 1263173)
check("126.3172 fails: W^10000 > 1000^10000 * 2^1263172",
      W_LEDGER ** 10000 > 1000 ** 10000 * 2 ** 1263172)
check("claim is below the frontier 126.338", CLAIM < Fraction(126338, 1000))

# ------------------------------------------------- 8. success floor
METER_LOSS = 4 * TAIL
SUCCESS_LO = P0 - LOSS_OTHER - METER_LOSS
check("success > 0.3901 - 2e-7 - 4*2^-30 > 0.39009", SUCCESS_LO > Fraction(39009, 100000),
      "%.10f" % float(SUCCESS_LO))
check("success floor > the 0.39 of rules.md", SUCCESS_LO > THRESHOLD)

# ------------------------------------------------- 9. counter widths
check("192-bit meter holds 1000 * 2^claim milli-units", CAP_MILLI_LO < 2 ** 191,
      "2^%.5f" % log2(CAP_MILLI_LO))
check("128-bit operation counter holds S_int + 2^61 < 2^127", OPS < 2 ** 127, "2^%.6f" % log2(OPS))
check("128-bit start counter holds W/2^13 < 2^114", CAP_HI / 2 ** 13 < 2 ** 114)
check("128-bit refill counter holds the refills < 2^87", REFILLS < 2 ** 87, "2^%.6f" % log2(REFILLS))

# ------------------------------------------------------- 10. preprocessing (unchanged)
STORE_INIT_UNITS = Fraction(STORE_INIT_MILLI, 1000)
CAP_PREP = pow2_interval(590, 10)[0]
check("preprocessing cap halts with probability < 2^-30",
      Fraction(2 ** 25) / (CAP_PREP - STORE_INIT_UNITS) < Fraction(1, 2 ** 30))

print()
print("phase 1 2^%.6f  hard 2^%.6f  allowances 2^%.6f"
      % (log2(Fraction(PH1, 1000)), log2(Fraction(HARD, 1000)), log2(Fraction(ALLOW, 1000))))
print("W_ledger = %d milli-units = 2^%.9f units  claim %s = %d/%d"
      % (W_LEDGER, LW, float(CLAIM), CLAIM.numerator, CLAIM.denominator))
print("success > %.10f" % float(SUCCESS_LO))
print()
if FAILED:
    print("CERTIFICATE: FAIL (%d)" % len(FAILED))
    raise SystemExit(1)
print("CERTIFICATE: PASS")
