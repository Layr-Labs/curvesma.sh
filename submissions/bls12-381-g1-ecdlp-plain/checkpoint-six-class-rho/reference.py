"""Six-class checkpoint rho: exact work bound for plain BLS12-381 G1.

The walk is an r-adding map on canonical order-6 orbits, with a jump table
of 2**20 affine points. Checkpoints are the only stored states. A repeated
checkpoint with independent coefficient pairs solves for the scalar. A
dependent repeat, or a trail that runs out of budget, is a cycle escape.

Units follow tracks/classical-group-work-v2.json. One batched incomplete
affine G1 addition costs 1. Every other kernel is an instruction count
divided by C=2000. Field multiplication is the pinned blst_fp_mul price
of 383 instructions, and S is charged as M.
"""

from __future__ import annotations

from fractions import Fraction

GROUP_ORDER = int(
    "73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001",
    16,
)
# Smaller root of X^2+X+1 = 0 over F_r. lambda**2 + lambda + 1 == r.
LAMBDA = 0xAC45A4010001A40200000000FFFFFFFF

C = 2000
FP_MUL_INS = 383
JUMPS = 1 << 20
N_STARTS = 64
# Checkpoint predicate: bits [20, 117) of the canonical x are zero.
# theta = 2**(-97). Slots above that rate are not allocated.
THETA_BITS = 97
MAX_CHECKPOINTS = 1 << 30
# Eight mean gaps. A geometric wait longer than this has probability < exp(-8).
DETECTION_STEPS = 8 << THETA_BITS
WALK_SLACK = 4

# Montgomery CIOS, 4-limb value by 2-limb lambda, plus conditional subtract.
# A direct limb count is about 151 instructions; 256 is the charged cap.
LAMBDA_MUL_INS = 256
N_LAMBDA_MULS = 2
# phi and phi^2 on the affine x coordinate. Negation does not touch x.
N_BETA_MULS = 2
# Six-limb compares, jump-scalar adds, index and checkpoint tests, loop.
BOOKKEEPING_INS = 128
STEP_INS = N_LAMBDA_MULS * LAMBDA_MUL_INS + N_BETA_MULS * FP_MUL_INS + BOOKKEEPING_INS

# a=0 Jacobian double, EFD dbl-2009-l: 2M+5S. Mixed add: 7M+4S.
DBL_M = 2 + 5
MADD_M = 7 + 4
SCALAR_BITS = 255
DOUBLES = SCALAR_BITS - 1
ADDS = SCALAR_BITS - 1
# Fermat inversion in Fp, p.bit_length()==381, plus 6M to reach affine.
FERMAT_M = 381 + 381
AFFINE_M = 6
JUMP_INS = (
    DOUBLES * DBL_M * FP_MUL_INS
    + ADDS * MADD_M * FP_MUL_INS
    + MADD_M * FP_MUL_INS
    + (FERMAT_M + AFFINE_M) * FP_MUL_INS
)
START_INS = JUMP_INS
VERIFY_INS = DOUBLES * DBL_M * FP_MUL_INS + ADDS * MADD_M * FP_MUL_INS + AFFINE_M * FP_MUL_INS
# Open-address probe charged at 200 instructions, once per allocated slot.
LOOKUP_INS_EACH = 200
LOOKUP_INS_TOTAL = MAX_CHECKPOINTS * LOOKUP_INS_EACH

JUMP_ENTRY_BYTES = 96 + 64
CHECKPOINT_SLOT_BYTES = 256
WALK_STATE_BYTES = 160
CODE_SCRATCH_BYTES = 1 << 20


def class_count() -> int:
    if (GROUP_ORDER - 1) % 6 != 0:
        raise ValueError("r is not 1 mod 6")
    return (GROUP_ORDER - 1) // 6


def _arctan_bounds(z: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    acc = Fraction(0)
    zp = z
    z2 = z * z
    lower = Fraction(0)
    upper = Fraction(0)
    for n in range(terms):
        term = zp / (2 * n + 1)
        if n % 2 == 0:
            acc += term
            upper = acc
        else:
            acc -= term
            lower = acc
        zp *= z2
    if terms < 2:
        raise ValueError("arctan terms")
    return lower, upper


def pi_upper(terms: int = 12) -> Fraction:
    """Machin: pi/4 = 4*arctan(1/5) - arctan(1/239), truncated from above."""
    hi5, = (_arctan_bounds(Fraction(1, 5), terms)[1],)
    lo239 = _arctan_bounds(Fraction(1, 239), terms)[0]
    return 4 * (4 * hi5 - lo239)


def sigma_upper_int(n_classes: int | None = None) -> int:
    """Smallest integer s >= sqrt(pi_upper * N / 2)."""
    n = class_count() if n_classes is None else n_classes
    pi = pi_upper()
    # s*s >= pi*N/2
    numer = pi.numerator * n
    denom = pi.denominator * 2
    lo = 0
    hi = 1
    while hi * hi * denom < numer:
        hi <<= 1
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid * mid * denom >= numer:
            hi = mid
        else:
            lo = mid
    return hi


def bound_holds(n_classes: int | None = None) -> bool:
    """Birthday envelope: i.i.d. uniform samples satisfy E[T] <= s + WALK_SLACK.

    P(T>k) <= exp(-k(k-1)/(2N)). The sum of that envelope is at most
    2 + exp(1/(8N))*(sigma+1/2). exp(x) < 1+x+x^2 for x=1/(8N), and
    sigma <= s, so the integer comparison below is sufficient. The walk
    uses this same integer as its step budget; identifying the adding walk
    with that envelope is heuristic H1 in proof.md, not this inequality.
    """
    n = class_count() if n_classes is None else n_classes
    s = sigma_upper_int(n)
    # (1 + 1/(8N) + 1/(64 N^2))*(s + 1/2) + 2 <= s + WALK_SLACK
    # with WALK_SLACK >= 4 iff the s+2 form holds.
    if WALK_SLACK < 4:
        return False
    left = (64 * n * n + 8 * n + 1) * (2 * s + 1)
    right = (s + 2) * 128 * n * n
    return left <= right


def walk_steps() -> int:
    if not bound_holds():
        raise ValueError("expectation slack failed")
    return sigma_upper_int() + WALK_SLACK + DETECTION_STEPS


def instruction_numerator() -> int:
    steps = walk_steps()
    return (
        steps * (C + STEP_INS)
        + JUMPS * JUMP_INS
        + N_STARTS * START_INS
        + VERIFY_INS
        + LOOKUP_INS_TOTAL
    )


def work() -> Fraction:
    return Fraction(instruction_numerator(), C)


def preprocessing() -> Fraction:
    return Fraction(JUMPS * JUMP_INS + N_STARTS * START_INS, C)


def peak_bytes() -> int:
    return (
        JUMPS * JUMP_ENTRY_BYTES
        + MAX_CHECKPOINTS * CHECKPOINT_SLOT_BYTES
        + N_STARTS * WALK_STATE_BYTES
        + CODE_SCRATCH_BYTES
    )


def _ln2_bounds(terms: int = 80) -> tuple[Fraction, Fraction]:
    acc = Fraction(0)
    for k in range(1, terms + 1):
        acc += Fraction(1, k * (1 << k))
    rem = Fraction(1, (terms + 1) * (1 << terms))
    return acc, acc + rem


def _ln_ratio_bounds(numer: int, denom: int, terms: int = 48) -> tuple[Fraction, Fraction]:
    if numer <= 0 or denom <= 0:
        raise ValueError("ln domain")
    if numer == denom:
        return Fraction(0), Fraction(0)
    if numer == 2 * denom:
        lo, hi = _ln2_bounds()
        return lo, hi
    if numer < denom or numer > 2 * denom:
        raise ValueError("ln ratio not in [1, 2]")
    w = Fraction(numer - denom, numer + denom)
    w2 = w * w
    acc = Fraction(0)
    wpow = w
    for j in range(terms):
        acc += wpow / (2 * j + 1)
        wpow *= w2
    rem = wpow / ((2 * terms + 1) * (1 - w2))
    return 2 * acc, 2 * (acc + rem)


def log2_bounds(n: int) -> tuple[Fraction, Fraction]:
    """Rational bounds with lower <= log2(n) <= upper."""
    if n <= 0:
        raise ValueError(n)
    e = n.bit_length() - 1
    if n == 1 << e:
        exact = Fraction(e)
        return exact, exact
    k = 4096
    lower_mant = (n << k) >> e
    upper_mant = lower_mant + 1
    if upper_mant > (1 << (k + 1)):
        upper_mant = 1 << (k + 1)
    ln2_lo, ln2_hi = _ln2_bounds()
    if lower_mant <= 1 << k:
        frac_lo = Fraction(0)
    else:
        ln_lo, _ln_hi = _ln_ratio_bounds(lower_mant, 1 << k)
        frac_lo = ln_lo / ln2_hi
    ln_lo_u, ln_hi_u = _ln_ratio_bounds(upper_mant, 1 << k)
    frac_hi = ln_hi_u / ln2_lo
    return e + frac_lo, e + frac_hi


def log2_fraction_ceiling_micro(value: Fraction, places: int = 6) -> int:
    """Smallest integer d with d/10**places >= log2(value)."""
    if value <= 0:
        raise ValueError(value)
    upper = log2_bounds(value.numerator)[1] - log2_bounds(value.denominator)[0]
    scale = 10 ** places
    return (upper.numerator * scale + upper.denominator - 1) // upper.denominator


def micro_to_decimal(micro: int, places: int = 6) -> str:
    scale = 10 ** places
    whole, frac = divmod(micro, scale)
    return f"{whole}.{frac:0{places}d}"


def time_log2_decimal(places: int = 6) -> str:
    return micro_to_decimal(log2_fraction_ceiling_micro(work(), places), places)


def memory_log2_decimal(places: int = 6) -> str:
    return micro_to_decimal(log2_fraction_ceiling_micro(Fraction(peak_bytes()), places), places)


def preprocessing_log2_decimal(places: int = 6) -> str:
    return micro_to_decimal(log2_fraction_ceiling_micro(preprocessing(), places), places)


def self_check() -> None:
    if LAMBDA * LAMBDA + LAMBDA + 1 != GROUP_ORDER:
        raise ValueError("lambda polynomial")
    if LAMBDA.bit_length() > 128:
        raise ValueError("lambda wider than 128 bits")
    if not bound_holds():
        raise ValueError("slack")
    if STEP_INS != 1406:
        raise ValueError("step instruction total drifted")
    lo3, hi3 = log2_bounds(3)
    # 1.584962 < log2(3) < 1.584963
    if not (lo3 < Fraction(1584963, 1000000) and hi3 > Fraction(1584962, 1000000)):
        raise ValueError("log2(3) bounds")
    if log2_fraction_ceiling_micro(Fraction(8)) != 3_000_000:
        raise ValueError("log2 power of two")
    micro3 = log2_fraction_ceiling_micro(Fraction(3))
    if micro3 != 1_584_963:
        raise ValueError(f"ceil log2(3) {micro3}")
    if walk_steps() <= sigma_upper_int():
        raise ValueError("budget")
    if pi_upper() >= Fraction(355, 113):
        raise ValueError("machin pi")
    if pi_upper() <= 3:
        raise ValueError("pi lower")


if __name__ == "__main__":
    self_check()
    w = work()
    print(time_log2_decimal())
    print(f"{w.numerator}/{w.denominator}")
    print(memory_log2_decimal())
    print(peak_bytes())
    print(preprocessing_log2_decimal())
    print(f"{preprocessing().numerator}/{preprocessing().denominator}")
    print(class_count())
    print(sigma_upper_int())
    print(walk_steps())
    print(instruction_numerator())
