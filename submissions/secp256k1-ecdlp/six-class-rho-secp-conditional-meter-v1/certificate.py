"""Exact-integer certificate for the work bound W <= 2^126.3029 of proof.md.

Every inequality below is one exact comparison of Python integers or
``fractions.Fraction`` values; floating point is used only to print.

  * kernel identities, prices, analytic rows, store sizes and hard rows are
    integer arithmetic (unchanged from the inherited package);
  * ln(1/(1 - p0)) is a rational interval from the atanh series with an explicit
    remainder bound, and fractional powers of two are integer-root intervals;
  * the conditional-meter allowance of proof.md section 3 is built from rational
    upper bounds: twice the doubling-cycle mean and eight times the traversal-cycle
    mean over the trail mass Lambda, the retrace row for K <= 2^44 collisions, and the inherited starts and
    small-charge rows at 2^30 times their expectations;
  * the four tails (Kolmogorov for doublings and traversal cycles, Markov for K,
    starts and small charges) and the success floor are exact rationals;
  * the claim is checked in both directions with exact powers.

Run:  python3 certificate.py        (exit status 0 = all checks hold)
"""
import math
import sys
from fractions import Fraction

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(400000)

FAILED = []


def check(name, ok, detail=""):
    """Print one named check and record a failure.

    Args:
        name: label printed for the check.
        ok: the exact boolean result.
        detail: optional value printed after the verdict.
    """
    print("%-62s %s%s" % (name, "ok" if ok else "FAIL", ("  " + detail) if detail else ""))
    if not ok:
        FAILED.append(name)


def log2(x):
    """Return log2 of a positive rational as a float (printing only)."""
    f = Fraction(x)
    return math.log2(f.numerator) - math.log2(f.denominator)


def ceil_fraction(f):
    """Return the integer ceiling of a rational."""
    f = Fraction(f)
    return -((-f.numerator) // f.denominator)


def iroot(n, k):
    """Return floor(n^(1/k)) for integers n >= 0, k >= 1."""
    if n < 2:
        return n
    b = n.bit_length()
    lo, hi = 1 << max(0, (b - 1) // k), 1 << ((b + k - 1) // k)
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
    """Return rationals (lo, hi) with lo <= 2^(num/den) < hi and hi - lo = 2^-scale."""
    sh = num + scale * den
    assert sh >= 0
    m = iroot(2 ** sh, den)
    assert m ** den <= 2 ** sh < (m + 1) ** den, "integer root certificate"
    return Fraction(m, 2 ** scale), Fraction(m + 1, 2 ** scale)


def ln_interval(x, bits=240):
    """Return rationals (lo, hi) bracketing ln(x) for rational x > 1 (atanh series)."""
    z = Fraction(x - 1, x + 1)
    z2, total, p, k = z * z, Fraction(0), z, 0
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
R_WALKS = 2048                       # table rows R
M_WALKS = 2 ** 20                    # concurrent walks M
L_MAX = 20 * 2 ** 64                 # trail cap
KAPPA = Fraction(101, 100)           # H1 slowdown allowance
P0 = Fraction(3900003, 10 ** 7)      # H1 floor p0 at S0 = S_ops / KAPPA
CLAIM_NUM, CLAIM_DEN = 1263029, 10000
CLAIM = Fraction(CLAIM_NUM, CLAIM_DEN)
THRESHOLD = Fraction(39, 100)        # rules.md: success at least 0.39
LOSS_OTHER = Fraction(2, 10 ** 7)    # inherited failure table, events other than the meter
MEAN_MULT_D = 2                      # allowance = 2 x mean for doubling-cycle entries
MEAN_MULT_V = 8                      # allowance = 8 x mean for traversal cycles (sparser evidence)
K_EXP = 44                           # retrace row covers K <= 2^44 inter-trail collisions
EK_BOUND = 2 ** 14                   # H3: E[K] <= 2^14
J_TRAILS = 2 ** 114                  # starts <= W / 2^13 < 2^114 (meter-enforced)
WINDOW = 127                         # Y'' depends on a window of 64 labels: 2*63 + 1 neighbours
MARKOV_SMALL = 2 ** 30               # starts and small charges: 2^30 x expectation

# ------------------------------------------------------- 1. counted kernels
MUL = 12 + 3 * 17 + 28
SQR = 12 + 9 + 7 + 15 + 28
SUB, CANON, ADD_WEAK, ADD_WEAK2 = 9, 8, 9, 14
PASS1 = MUL + SUB + 5 + 4 + 2
PASS2 = 4 * MUL + SQR + 5 * SUB + 2 * CANON + 19 + 6 + 2
check("fe_mul costs 91 instructions", MUL == 91, str(MUL))
check("fe_sqr costs 71 instructions", SQR == 71, str(SQR))
check("pass 1 iteration is 111", PASS1 == 111, str(PASS1))
check("pass 2 iteration is 523", PASS2 == 523, str(PASS2))
ROWS = [24, 8, 13, 17, 7, 5, 26, 4, 13, 1]
check("analytic rows sum to 118", sum(ROWS) == 118, str(sum(ROWS)))
STEP = PASS1 + PASS2 + sum(ROWS)
DBL = 5 * MUL + 2 * SQR + ADD_WEAK + 2 * ADD_WEAK2 + 4 * SUB + 2 * CANON + 72 + sum(ROWS) + 50 + 30
TRAV = STEP + 24
PRICE_STEP, PRICE_TRAV, PRICE_DBL = 752, 776, 956
check("step is 111 + 523 + 118 = 752", STEP == PRICE_STEP, str(STEP))
check("doubling bound is 956", DBL == PRICE_DBL, str(DBL))
check("traversal bound is 776", TRAV == PRICE_TRAV, str(TRAV))
SUR_DBL, SUR_TRAV = PRICE_DBL - PRICE_STEP, PRICE_TRAV - PRICE_STEP     # milli-units
check("surcharges are 204 and 24 milli-units", (SUR_DBL, SUR_TRAV) == (204, 24))

# ------------------------------------------------------- 2. the group algebra
check("p = 2^256 - (2^32 + 977)", (1 << 256) - P_FIELD == (1 << 32) + 977)
check("r - 1 is divisible by 6", (r - 1) % 6 == 0)
N = (r - 1) // 6
check("log2 N = 253.41504 to 1e-5", abs(log2(N) - 253.41504) < 1e-5, "%.5f" % log2(N))

# ------------------------------------------------------- 3. the report store
NODES = 2 ** 16
DCAP = ceil_fraction(pow2_interval(6383, 100)[1])
B = ceil_fraction(Fraction(2 * DCAP, NODES))
REPORT_PRICE = 40000 + 8 * 50
OVF_MILLI = DCAP * ((6 * B + 35) + 16 * B)
STORE_INIT_MILLI = 2 * DCAP * 16
check("D_cap = ceil(2^63.83) = 16396253461269537204", DCAP == 16396253461269537204)
check("B = ceil(2*D_cap/2^16) = 500373945961595", B == 500373945961595)
check("whole store < 2^72 bytes", NODES * B * 112 < 2 ** 72, "2^%.3f" % log2(NODES * B * 112))

# ------------------------------------------------------- 4. H1 numbers
L_LO, L_HI = ln_interval(Fraction(10 ** 7, 10 ** 7 - 3900003))     # ln(1/(1 - p0))
check("ln(1/(1-p0)) interval is tight", 0 < L_LO < L_HI and L_HI - L_LO < Fraction(1, 2 ** 200),
      "%.15f" % float(L_LO))
N_PRIME = Fraction(N * R_WALKS, R_WALKS - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI                    # upper bound on S_ops^2
m = iroot((V.numerator * 10 ** 80 + V.denominator - 1) // V.denominator, 2) + 1
S_HI = Fraction(m, 10 ** 40)
S_INT = ceil_fraction(S_HI)
check("S_ops <= 2^126.71396 (certified upper bound)", S_HI * S_HI >= V and log2(S_HI) < 126.71396,
      "%.9f" % log2(S_HI))

# ------------------------------------------------------- 5. phase 1 and hard extras
N_OPS = S_INT + 2 ** 61                                   # covers the < 2^60.41 block overshoot
PH1 = PRICE_STEP * N_OPS                                  # milli-units
check("phase 1 = 752*(S_int + 2^61) < 2^126.30276", log2(Fraction(PH1, 1000)) < 126.30276,
      "%.9f" % log2(Fraction(PH1, 1000)))
CLAIM_LO, CLAIM_HI = pow2_interval(CLAIM_NUM, CLAIM_DEN)
CAP_MILLI_LO = (CLAIM_LO.numerator * 1000) // CLAIM_LO.denominator
REFILLS = -((-CAP_MILLI_LO) // (2 ** 40 * 1000)) + 2 ** 20
HARD_ITEMS = [2 ** 61 * 1000,                                       # meter remainders
              REFILLS * 30,                                         # block refills
              2 ** 20 * (L_MAX + 65) * 600 * 1000,                  # phase 2
              16 * ((L_MAX + 65) * 600 + 2 ** 13) * 1000,           # re-walks
              DCAP * REPORT_PRICE,                                  # reports
              OVF_MILLI,                                            # overflow insert and lookup
              STORE_INIT_MILLI,                                     # store initialisation
              2 ** 25 * 1000,                                       # table and verification
              2 ** 26 * 1000]                                       # node set-up
HARD = sum(HARD_ITEMS)
check("hard extras < 2^107.156 units", log2(Fraction(HARD, 1000)) < 107.156,
      "%.9f" % log2(Fraction(HARD, 1000)))

# ------------------------------------------------------- 6. conditional meter allowance
Q_D = Fraction(1, 6 * R_WALKS) + pow2_interval(-199954, 10000)[1]   # doubling-cycle entries / fresh arrival
Q_V = pow2_interval(-199954, 10000)[1]                              # cycles of length 3..64 / fresh arrival
check("q_D <= 2^-13.5681", Q_D <= 1 / pow2_interval(135681, 10000)[0], "2^%.9f" % log2(Q_D))

CAPPED_MAX = ceil_fraction(Fraction(N_OPS, L_MAX))                  # each capped trail used L_max operations
OPEN_AND_CAPPED = 65 * CAPPED_MAX + M_WALKS * (L_MAX + 65)
check("65*capped + M*(L_max+65) < 2^89", OPEN_AND_CAPPED < 2 ** 89, "2^%.4f" % log2(OPEN_AND_CAPPED))
LAMBDA = N_OPS + 2 ** 89 + 2 ** K_EXP * (L_MAX + 65)                # trail mass bound of section 3.3
check("Lambda < 2^126.71409", log2(LAMBDA) < 126.71409, "2^%.9f" % log2(LAMBDA))

A_D = ceil_fraction(SUR_DBL * MEAN_MULT_D * Q_D * LAMBDA)
A_V = ceil_fraction(SUR_TRAV * 64 * MEAN_MULT_V * Q_V * LAMBDA)
A_K = SUR_DBL * 2 ** K_EXP * (L_MAX + 65)
A_S = ceil_fraction(MARKOV_SMALL * 2 ** 13 * (S_HI / 2 ** 64 + 2 ** 20) * 1000)
A_O = MARKOV_SMALL * 2 ** 70 * 1000
ALLOW = A_D + A_V + A_K + A_S + A_O
for name, v in (("doublings 0.204*2*q_D*Lambda", A_D), ("traversals 0.024*64*8*q_V*Lambda", A_V),
                ("retrace 0.204*2^44*(L_max+65)", A_K), ("starts 2^30*E", A_S),
                ("small charges 2^30*2^70", A_O), ("allowance total", ALLOW)):
    print("  %-40s 2^%.7f units" % (name, log2(Fraction(v, 1000))))
check("allowance < 2^112.58 units", log2(Fraction(ALLOW, 1000)) < 112.58)

# Kolmogorov: P(max_{j<=J} S'_j >= (c-1)*Q*Lambda) <= J * WINDOW * Q * 2^64 / ((c-1)*Q*Lambda)^2
TAIL_D = Fraction(J_TRAILS * WINDOW * 2 ** 64) / (Q_D * (MEAN_MULT_D - 1) ** 2 * LAMBDA * LAMBDA)
TAIL_V = Fraction(J_TRAILS * WINDOW * 2 ** 64) / (Q_V * (MEAN_MULT_V - 1) ** 2 * LAMBDA * LAMBDA)
TAIL_K = Fraction(EK_BOUND, 2 ** K_EXP)
TAIL_S = Fraction(1, MARKOV_SMALL)
TAIL_O = Fraction(1, MARKOV_SMALL)
NFR_ALL = Fraction((L_MAX + 65) ** 2 * 2 ** 63, 2 * N)              # non-fruitless return, any phase-1 trail
check("doubling tail < 2^-54", TAIL_D < Fraction(1, 2 ** 54), "2^%.4f" % log2(TAIL_D))
check("traversal-cycle tail < 2^-54", TAIL_V < Fraction(1, 2 ** 54), "2^%.4f" % log2(TAIL_V))
check("collision tail = 2^14/2^44 = 2^-30", TAIL_K == Fraction(1, 2 ** 30))
check("non-fruitless return in any trail < 2^-54", NFR_ALL < Fraction(1, 2 ** 54), "2^%.4f" % log2(NFR_ALL))
METER_LOSS = TAIL_D + TAIL_V + TAIL_K + TAIL_S + TAIL_O + NFR_ALL
check("meter loss < 2.8e-9", METER_LOSS < Fraction(28, 10 ** 10), "%.4e" % float(METER_LOSS))

# ------------------------------------------------------- 7. total, success and claim
W_LEDGER = PH1 + HARD + ALLOW
check("W_ledger < 2^126.30287", log2(Fraction(W_LEDGER, 1000)) < 126.30287,
      "%.9f" % log2(Fraction(W_LEDGER, 1000)))
check("W_ledger <= 1000 * 2^126.3029", W_LEDGER <= CAP_MILLI_LO)
check("126.3029 holds: W^10000 < 1000^10000 * 2^1263029",
      W_LEDGER ** 10000 < 1000 ** 10000 * 2 ** 1263029)
check("126.3028 fails: W^10000 > 1000^10000 * 2^1263028",
      W_LEDGER ** 10000 > 1000 ** 10000 * 2 ** 1263028)
SUCCESS_LO = P0 - LOSS_OTHER - METER_LOSS
check("success >= 0.3900003 - 2e-7 - meter loss > 0.39", SUCCESS_LO > THRESHOLD,
      "%.10f" % float(SUCCESS_LO))

# ------------------------------------------------------- 8. counter widths
check("192-bit meter holds the claim", CAP_MILLI_LO < 2 ** 192, "2^%.6f" % log2(CAP_MILLI_LO))
check("128-bit operation counter holds S_ops + 2^61", S_HI + 2 ** 61 < 2 ** 127)
check("128-bit start counter holds W/2^13 < J = 2^114", Fraction(CAP_MILLI_LO, 1000) / 2 ** 13 < J_TRAILS)
check("128-bit refill counter < 2^87", REFILLS < 2 ** 87, "2^%.6f" % log2(REFILLS))

# ------------------------------------------------------- 9. preprocessing
CAP_PREP = pow2_interval(590, 10)[1]
PREP_EXPECTED = Fraction(STORE_INIT_MILLI, 1000) + 2 ** 26 + pow2_interval(241, 10)[1] + 2 ** 33
check("expected preprocessing < 2^59.0 units", PREP_EXPECTED < CAP_PREP, "2^%.5f" % log2(PREP_EXPECTED))
check("preprocessing cap halts with probability < 2^-30",
      Fraction(2 ** 25) / (CAP_PREP - Fraction(STORE_INIT_MILLI, 1000)) < pow2_interval(-30, 1)[1])

print()
print("phase 1 2^%.7f  hard 2^%.7f  allowance 2^%.7f"
      % (log2(Fraction(PH1, 1000)), log2(Fraction(HARD, 1000)), log2(Fraction(ALLOW, 1000))))
print("W_ledger 2^%.9f units  claim %s" % (log2(Fraction(W_LEDGER, 1000)), CLAIM))
print("success > %.10f" % float(SUCCESS_LO))
print()
if FAILED:
    print("CERTIFICATE: FAIL (%d)" % len(FAILED))
    raise SystemExit(1)
print("CERTIFICATE: PASS")
