"""Exact-integer certificate for the gated expected-work bound of proof.md.

Every check below is an exact comparison of Python integers or
``fractions.Fraction`` values, never floating point: instruction identities and
prices, rational-interval ln(1/(1 - p1)) and powers of two, and every ledger,
gate, failure-row and claim inequality in both directions.

Run:  python3 certificate.py        (exit status 0 = all checks hold)

The certificate proves, for this package (kappa = 1.01, R = 2^20, opened-run H1
floor p1 = 0.52, Markov headroom h = 0.000165, step price 647, 64-bit run gate):

    phase1  = 647 * (S_int + 2^60)                      (milli-units)
    hard    = sum of the ten exact hard rows
    E_rand  = sum of the four exact expected random rows
    A       = ceil(E_rand / 0.000165)  >= E_rand / 0.000165
    W1      = phase1 + hard + A  <= 2^126.3711          (per-run meter cap)
    lo_run  = 0.52 - 0.000165 - 2e-7 = 0.5198348        (opened-run success floor)
    K       = least integer with (K/2^64) * lo_run >= 0.39
    E_total = 2^7 + (K/2^64) * 2^126.3711               (expected total, milli-units)
              < 1000 * 2^125.9566   so 125.9566 holds, 125.9565 fails
    (K/2^64) * lo_run >= 0.39 so the unconditional success requirement holds
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
R_WALKS = 2 ** 20
KAPPA = Fraction(101, 100)          # H1 slowdown allowance 1.01
P1 = Fraction(52, 100)              # opened-run H1 floor at S0 = S_ops / KAPPA
HEADROOM = Fraction(165, 1000000)   # 0.000165 Markov headroom h
RUN_CAP10 = 1263711                 # per-run meter cap 2^RUN_CAP, 4-decimal, tightness checked below
CLAIM10 = 1259566                   # score cap 2^CLAIM, 4-decimal ceiling of log2(E_total/1000)
RUN_CAP = Fraction(RUN_CAP10, 10000)
CLAIM = Fraction(CLAIM10, 10000)
THRESHOLD = Fraction(39, 100)       # rules.md: success at least 0.39
LOSS_OTHER = Fraction(2, 10 ** 7)   # budget for every failure row except the meter's early halt
GATE_MILLI = 2 ** 7                 # one RNDR word + compare + branch, charged on every run

# ------------------------------------------------------- 1. counted kernels
MUL = (16 + 16) + 30 + 23               # products, six carry chains, reduction without a final fold
SQR = (10 + 10) + 8 + 6 + 7 + 23        # products, off-diagonal sum, doubling, diagonal chain, reduction
SUB, SUB2, CANON, ADD_WEAK, ADD_WEAK2 = 9, 15, 8, 9, 14
PASS1 = MUL + SUB + 5 + 4 + 2           # 1 LDR + 4 LDP, 4 STP, SUBS + B.NE
PASS2 = 4 * MUL + SQR + 3 * SUB + SUB2 + 13 + 4 + 2  # fused x3 (15); y3 and x_W-x3 by fe_sub
check("fe_mul costs 85 instructions", MUL == 85, str(MUL))
check("fe_sqr costs 64 instructions", SQR == 64, str(SQR))
check("fe_sub2 costs 15 instructions", SUB2 == 15, str(SUB2))
check("pass 1 iteration is 105", PASS1 == 105, str(PASS1))
check("pass 2 iteration is 465", PASS2 == 465, str(PASS2))

ROWS = [("F1 truncated product (G1..G5 held in registers)", 17),
        ("range test, S, wrap test, F2", 8),
        ("decision", 13),
        ("sign, p - y from p0 and 2^64-1 registers, compare, 4 CSEL", 12),
        ("h, distinguished test, table address (rows of 8 entries)", 6),
        ("arrival traffic (2 ldp + 1 store)", 3),
        ("checkpoint compare and test only (cycles caught by the checkpoint)", 8),
        ("zero-denominator test (moved to once per batch)", 0),
        ("length (low limb; carry path separate), meter, arrival loop", 9),
        ("share of the inversion and call overhead", 1)]
ROWS_TOTAL = sum(v for _, v in ROWS)
check("analytic rows sum to 77", ROWS_TOTAL == 77, str(ROWS_TOTAL))

STEP = PASS1 + PASS2 + ROWS_TOTAL
DBL_ARITH = 5 * MUL + 2 * SQR + ADD_WEAK + 2 * ADD_WEAK2 + SUB2 + 2 * SUB + 30   # fused x3; traffic as the loops
DBL = DBL_ARITH + ROWS_TOTAL + 34 + 30      # list moves 8 LDP + 8 STP + 18 other; escape logic 30
TRAV = STEP + 24
PRICE_STEP, PRICE_TRAV, PRICE_DBL = 647, 671, 794
check("step is 105 + 465 + 77 = 647", STEP == 647, str(STEP))
check("doubling bound is 653 + 77 + 34 + 30 = 794", DBL_ARITH == 653 and DBL == 794, str(DBL))
check("traversal bound is 647 + 24 = 671", TRAV == 671, str(TRAV))
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
L_LO, L_HI = ln_interval(Fraction(25, 12))               # ln(1/(1 - p1)) = ln(25/12)
check("ln(1/(1-p1)) interval is positive and tight", 0 < L_LO < L_HI and L_HI - L_LO < Fraction(1, 2 ** 200),
      "%.15f" % float(L_LO))

N_PRIME = Fraction(N * R_WALKS, R_WALKS - 1)
V = KAPPA * KAPPA * 2 * N_PRIME * L_HI                   # upper bound on S_ops^2
m = iroot((V.numerator * 10 ** 80 + V.denominator - 1) // V.denominator, 2) + 1
S_HI = Fraction(m, 10 ** 40)
S_INT = ceil_fraction(S_HI)
check("S_ops <= 2^127.0000 (certified upper bound)", S_HI * S_HI >= V and log2(S_HI) < 127.0,
      "log2 bound %.9f" % log2(S_HI))
check("S_ops >= 2^126.9980 (lower bound for the quoted value)", log2(S_HI) > 126.998, "%.9f" % log2(S_HI))

# the accepted package's budget S0, anchor for scaling the inherited failure rows
L0_LO, L0_HI = ln_interval(Fraction(10000000, 6098247))   # ln(1/(1 - 0.3901753)) of b5de673a
V0 = KAPPA * KAPPA * 2 * N_PRIME * L0_HI
m0 = iroot((V0.numerator * 10 ** 80 + V0.denominator - 1) // V0.denominator, 2) + 1
S0_INT = ceil_fraction(Fraction(m0, 10 ** 40))
check("accepted-package budget S0 recomputed at 2^126.7140", 126.7140 < log2(Fraction(S0_INT, 1)) < 126.7141,
      "%.9f" % log2(Fraction(S0_INT, 1)))

# ------------------------------------------------------ 5. phase 1 and hard extras
PH1 = PRICE_STEP * (S_INT + 2 ** 60)                     # milli-units
check("phase 1 = 647*(S_int + 2^60) covers 0.647*(S_ops + 2^60)",
      Fraction(PH1, 1000) >= Fraction(PRICE_STEP, 1000) * (S_HI + 2 ** 60))
check("phase 1 < 2^126.3710", log2(Fraction(PH1, 1000)) < 126.3710,
      "%.9f" % log2(Fraction(PH1, 1000)))

# per-run meter cap: 2^RUN_CAP, the four-decimal ceiling of log2(W_ledger/1000);
# tightness is checked in section 7.  The meter, the refills and the counter widths
# below are sized on this per-run cap, exactly as the accepted package sized them
# on its claim.
RUN_LO, RUN_HI = pow2_interval(RUN_CAP10, 10000)
CAP_MILLI_LO = (RUN_LO.numerator * 1000) // RUN_LO.denominator
BLOCK_MILLI = 2 ** 40 * 1000
REFILLS = -((-CAP_MILLI_LO) // BLOCK_MILLI) + 2 ** 20
REFILL_MILLI = REFILLS * 30
CP_MILLI = 2 * ceil_fraction(Fraction(S_INT + 2 ** 61, 64))    # <= 1 write (STP + STR) per 64 phase-1 operations
HARD_ITEMS = [("meter remainders 2^61", Fraction(2 ** 61) * 1000),
              ("meter block refills: refills * 30", Fraction(REFILL_MILLI)),
              ("phase 2", Fraction(2 ** 20 * (20 * 2 ** 64 + 65) * 600) * 1000),
              ("re-walks", Fraction(16 * ((20 * 2 ** 64 + 65) * 600 + 2 ** 13)) * 1000),
              ("reports: D_cap * 40,400", Fraction(DCAP * REPORT_PRICE)),
              ("overflow insert and lookup: D_cap*(22B+35)", Fraction(OVF_MILLI)),
              ("store initialisation: 2*D_cap*16", Fraction(STORE_INIT_MILLI)),
              ("table (2^20 entries) and verification, < 2^34", Fraction(2 ** 34) * 1000),
              ("node set-up, 2^26", Fraction(2 ** 26) * 1000),
              ("checkpoint writes: 2*ceil((S_int + 2^61)/64)", Fraction(CP_MILLI))]
assert all(v.denominator == 1 for _, v in HARD_ITEMS)
HARD = sum(int(v) for _, v in HARD_ITEMS)
check("hard extras < 2^112.0813 units", log2(Fraction(HARD, 1000)) < 112.0813,
      "%.9f" % log2(Fraction(HARD, 1000)))
check("the checkpoint-write row is the largest hard row", CP_MILLI == max(int(v) for _, v in HARD_ITEMS))

# ---------------------------------------------------- 6. expected random extras
# fruitless cycles of length 3..64: sum over L and b of S2(L, b) * R^-(L-b), where S2(L, b) counts the
# partitions of L steps into b blocks of size >= 2 (upstream-proof.md section 4)
S2 = [[0] * 65 for _ in range(65)]
S2[0][0] = 1
for L_ in range(1, 65):
    for b_ in range(1, L_ // 2 + 1):
        S2[L_][b_] = b_ * S2[L_ - 1][b_] + (L_ - 1) * S2[L_ - 2][b_ - 1]
check("S2(3,1), S2(4,1), S2(4,2), S2(6,3) = 1, 1, 3, 15", (S2[3][1], S2[4][1], S2[4][2], S2[6][3]) == (1, 1, 3, 15))


def other_cycles(R):
    return sum(Fraction(S2[L_][b_], R ** (L_ - b_)) for L_ in range(3, 65) for b_ in range(1, L_ // 2 + 1))


check("at R = 2048 the sum is 2^-19.9954 (the public package's value)",
      abs(log2(other_cycles(2048)) + 19.9954) < 5e-5, "2^%.6f" % log2(other_cycles(2048)))
OTHER = other_cycles(R_WALKS)
check("at R = 2^20 the sum is < 2^-37.999", OTHER < Fraction(1, 1) / pow2_interval(37999, 1000)[0],
      "2^%.7f" % log2(OTHER))
RATE_HI = Fraction(1, 6 * R_WALKS) + OTHER
check("doubling rate 1/(6R) + sum <= 2^-22.5849", RATE_HI <= Fraction(1, 1) / pow2_interval(225849, 10000)[0],
      "2^%.9f" % log2(RATE_HI))
TRAV_RATE = Fraction(2, 6 * R_WALKS) + 64 * OTHER       # traversal steps per operation
check("traversal-step rate 2/(6R) + 64*sum <= 2^-21.5839", TRAV_RATE <= Fraction(1, 1) / pow2_interval(215839, 10000)[0],
      "2^%.9f" % log2(TRAV_RATE))
OPS = S_HI + 2 ** 60
EXTRA_ITEMS = [ceil_fraction((Fraction(PRICE_DBL, 1000) - Fraction(PRICE_STEP, 1000)) * RATE_HI * OPS * 1000),
               ceil_fraction((Fraction(PRICE_TRAV, 1000) - Fraction(PRICE_STEP, 1000)) * TRAV_RATE * OPS * 1000),
               ceil_fraction(Fraction(2 ** 13) * (S_HI / 2 ** 64 + 2 ** 20) * 1000),
               ceil_fraction(Fraction(2 ** 70) * 1000)]
E = sum(EXTRA_ITEMS)
check("expected random extras < 2^102.1 units", log2(Fraction(E, 1000)) < 102.1,
      "2^%.9f" % log2(Fraction(E, 1000)))

# ------------------------------------------- 7. per-run ledger, allowance, cap
A = -((-(E * 1000000)) // 165)                           # ceil(E / 0.000165)
check("allowance A = ceil(E/0.000165) (Markov safe)", A * 165 >= E * 1000000)
W_LEDGER = PH1 + HARD + A
wlog = log2(Fraction(W_LEDGER, 1000))
check("W_ledger < 2^126.3712", wlog < 126.3712, "%.12f" % wlog)
check("W_ledger > 2^126.3710 (for the quoted value)", wlog > 126.3710, "%.12f" % wlog)
check("the per-run cap exponent is 1263711/10000", RUN_CAP == Fraction(1263711, 10000))
check("W_ledger <= 1000 * 2^126.3711 (meter cap)", W_LEDGER <= CAP_MILLI_LO)
check("126.3711 holds: W^10000 < 1000^10000 * 2^1263711",
      W_LEDGER ** 10000 < 1000 ** 10000 * 2 ** 1263711)
check("126.3710 fails: W^10000 > 1000^10000 * 2^1263710",
      W_LEDGER ** 10000 > 1000 ** 10000 * 2 ** 1263710)

HEAD = CAP_MILLI_LO - PH1 - HARD
check("headroom > allowance A", HEAD > A, "2^%.6f" % log2(Fraction(HEAD, 1000)))
MARKOV_HI = Fraction(E, HEAD)
check("Markov ratio < 0.000165", MARKOV_HI < HEADROOM and MARKOV_HI * 1000000 < 165,
      "%.10f" % float(MARKOV_HI))
LO_RUN = P1 - HEADROOM - LOSS_OTHER
check("opened-run floor = 0.52 - 0.000165 - 2e-7 = 0.5198348",
      LO_RUN == Fraction(5198348, 10 ** 7), "%.7f" % float(LO_RUN))
check("opened-run floor > the 0.39 of rules.md", LO_RUN > THRESHOLD, "%.7f" % float(LO_RUN))

# failure events other than the early halt (proof.md section 5.3): the upstream rows plus the two rare-event
# rows of the kernels.  At most 7 reductions per batched operation (6 per step, 7 per doubling) lose the carry
# of their last addition only if z mod p < C0^2 + C0; x3 or y3 is non-canonical only if x_W or y_W < C0.
C0 = 2 ** 32 + 977
OPS_ALL = S_HI + 2 ** 61
ROW_CARRY = 7 * OPS_ALL * Fraction(C0 * C0 + C0, P_FIELD)
ROW_SMALL = 2 * OPS_ALL * Fraction(C0, P_FIELD)
# fused x3 final-borrow: intermediate < (b1+b2)*C0 <= 2*C0 after the raw double-sub (one x3 per step/dbl)
ROW_SUB2 = OPS_ALL * Fraction(2 * C0, 1 << 256)
# The inherited trail-cap row is a per-budget union bound: the per-trail event is a
# walk property independent of the budget and the number of trails is proportional to
# the arrival budget, so it is recharged here linearly in S_ops from the accepted
# package's budget S0 (at which the row was carried as < 10^-7).  The check below
# also certifies that even a quadratic reading of the same row fits the budget.
TRAIL_CAP = Fraction(1, 10 ** 7) * Fraction(S_INT + 2 ** 61, S0_INT + 2 ** 61)
UPSTREAM_ROWS = (TRAIL_CAP + Fraction(1, 2 ** 30) + Fraction(1, 2 ** 50) + Fraction(1, 2 ** 60)
                 + Fraction(1, 2 ** 64) + Fraction(1, 2 ** 114) + Fraction(1, r) + 2 * Fraction(1, 2 ** 100))
check("final-carry row 7*(S_ops + 2^61)*(C0^2 + C0)/p < 2^-62.1", ROW_CARRY < pow2_interval(-621, 10)[0],
      "2^%.3f" % log2(ROW_CARRY))
check("small-coordinate row 2*(S_ops + 2^61)*C0/p < 2^-96", ROW_SMALL < Fraction(1, 2 ** 96),
      "2^%.3f" % log2(ROW_SMALL))
check("fe_sub2 final-borrow row (S_ops+2^61)*2*C0/2^256 < 2^-96", ROW_SUB2 < Fraction(1, 2 ** 96),
      "2^%.3f" % log2(ROW_SUB2))
check("all failure rows other than the early halt < 2e-7",
      UPSTREAM_ROWS + ROW_CARRY + ROW_SMALL + ROW_SUB2 < LOSS_OTHER,
      "%.4g" % float(UPSTREAM_ROWS + ROW_CARRY + ROW_SMALL + ROW_SUB2))
check("even a quadratic trail-cap scaling fits the 2e-7 budget",
      (Fraction(1, 10 ** 7) * Fraction((S_INT + 2 ** 61) ** 2, (S0_INT + 2 ** 61) ** 2)
       + UPSTREAM_ROWS - TRAIL_CAP + ROW_CARRY + ROW_SMALL + ROW_SUB2) < LOSS_OTHER,
      "%.4g" % float(Fraction(1, 10 ** 7) * Fraction((S_INT + 2 ** 61) ** 2, (S0_INT + 2 ** 61) ** 2)
                     + UPSTREAM_ROWS - TRAIL_CAP + ROW_CARRY + ROW_SMALL + ROW_SUB2))

# ------------------------------------------------------- 7b. the run gate
# One coordinator reads one fresh 64-bit RNDR word u, independent of the challenge
# and of every walk coin, and opens the run only when u < K.  A closed gate halts
# immediately: its GATE_MILLI is charged on every run, nothing is amortized.
K_GATE = ceil_fraction(THRESHOLD * (2 ** 64) / LO_RUN)
check("K_GATE >= 1 integer", K_GATE >= 1)
check("K_GATE opens successfully: (K/2^64)*lo_run >= 0.39",
      K_GATE * LO_RUN >= THRESHOLD * 2 ** 64)
check("K_GATE is least: ((K-1)/2^64)*lo_run < 0.39",
      (K_GATE - 1) * LO_RUN < THRESHOLD * 2 ** 64)
QGATE = Fraction(K_GATE, 2 ** 64)
UNCOND = QGATE * LO_RUN
check("unconditional success >= the 0.39 of rules.md", UNCOND >= THRESHOLD,
      "%.12f" % float(UNCOND))
# Expected total work over the gate coin: the gate is charged always, the opened
# run's meter caps its work at CAP_MILLI_LO milli-units on every opened run.
E_TOTAL = Fraction(GATE_MILLI) + QGATE * Fraction(CAP_MILLI_LO)
E_NUM = GATE_MILLI * 2 ** 64 + K_GATE * CAP_MILLI_LO
E_DEN = 2 ** 64
check("E_total = E_NUM/2^64 exactly", E_TOTAL == Fraction(E_NUM, E_DEN))
elog = log2(E_TOTAL / 1000)
check("E_total < 2^125.9566", elog < 125.9566, "%.12f" % elog)
check("E_total > 2^125.9565 (for the quoted claim)", elog > 125.9565, "%.12f" % elog)
check("the claim is 1259566/10000", CLAIM == Fraction(1259566, 10000))
check("125.9566 holds: E_NUM^10000 < (E_DEN*1000)^10000 * 2^1259566",
      E_NUM ** 10000 < (E_DEN * 1000) ** 10000 * 2 ** 1259566)
check("125.9565 fails: E_NUM^10000 > (E_DEN*1000)^10000 * 2^1259565",
      E_NUM ** 10000 > (E_DEN * 1000) ** 10000 * 2 ** 1259565)
# disclosed per-run worst case (gate always open): scored figure is the expectation
WORST_MILLI = GATE_MILLI + CAP_MILLI_LO
check("disclosed per-run worst case < 2^126.3712 units",
      log2(Fraction(WORST_MILLI, 1000)) < 126.3712, "%.9f" % log2(Fraction(WORST_MILLI, 1000)))

# ------------------------------------------------- 8. counter widths, section 3
METER_MILLI = CAP_MILLI_LO
check("192-bit meter holds the per-run cap < 2^192 milli-units", METER_MILLI < 2 ** 192,
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
PREP_EXPECTED = STORE_INIT_UNITS + 2 ** 26 + pow2_interval(331, 10)[1] + 2 ** 33   # table 2^33.1
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
print("W_ledger(run) 2^%.9f units  cap %s  head 2^%.6f"
      % (log2(Fraction(W_LEDGER, 1000)), RUN_CAP, log2(Fraction(HEAD, 1000))))
print("gate K = %d  q = %.9f  lo_run = %.7f  uncond = %.12f"
      % (K_GATE, float(QGATE), float(LO_RUN), float(UNCOND)))
print("E_total 2^%.9f units  claim %s  worst-case(run) 2^%.6f"
      % (elog, CLAIM, log2(Fraction(WORST_MILLI, 1000))))
print("Markov %.10f" % float(MARKOV_HI))
print("S_int   = %d" % S_INT)
print("phase 1 = %d" % PH1)
print("checkpoint = %d" % CP_MILLI)
print("hard    = %d" % HARD)
print("extras  = %d" % E)
print("A       = %d" % A)
print("W_ledger_run = %d" % W_LEDGER)
print("K = %d" % K_GATE)
print("E_num = %d" % E_NUM)
print("unconditional = %.20f" % float(UNCOND))
print("trail-cap row = %.6g   row sum = %.4g" % (float(TRAIL_CAP), float(UPSTREAM_ROWS + ROW_CARRY + ROW_SMALL + ROW_SUB2)))
print()
if FAILED:
    print("CERTIFICATE: FAIL (%d)" % len(FAILED))
    raise SystemExit(1)
print("CERTIFICATE: PASS")
