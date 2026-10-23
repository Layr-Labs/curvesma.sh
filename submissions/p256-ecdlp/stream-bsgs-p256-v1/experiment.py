"""Bounded stream-join/Jacobian/Montgomery evidence; execute only in the shared verifier.

One JSON request on stdin, one JSON result on stdout. The fixture harness knows
test scalars; recover() receives only group targets and generated records.
No event counter here is an AArch64 instruction count or production benchmark.
"""

import json
import sys


FIELD = int("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16)
ORDER = int("ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551", 16)
COEFF_A = int("ffffffff00000001000000000000000000000000fffffffffffffffffffffffc", 16)
COEFF_B = int("5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b", 16)
GENERATOR = (
    int("6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296", 16),
    int("4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5", 16),
)
# Public profile fixture, the repository's own P-256 doubling vector (tests/test_p256.py).
TWICE_GENERATOR = (
    int("7cf27b188d034f7e8a52380304b51ac3c08969e277f21b35a60b48fc47669978", 16),
    int("07775510db8ed040293d9ac69f7430dbba7dade63ce982299e04b79d227873d1", 16),
)
J_INFINITY = (0, 1, 0)
WORD = 1 << 64
MASK = WORD - 1
LIMBS = 4
R_MONT = 1 << 256
FIELD_LIMBS = [(FIELD >> (64 * k)) & MASK for k in range(LIMBS)]
# (Z/2^64)* has exponent 2^62, so FIELD^(2^62-1) is FIELD^-1 mod 2^64 (odd FIELD).
P_PRIME = (-pow(FIELD, (1 << 62) - 1, WORD)) % WORD
R2_MOD_P = R_MONT * R_MONT % FIELD
ONE_MONT = R_MONT % FIELD
COUNTS = {"mont_mul": 0, "limb_mul_pairs": 0, "limb_add_steps": 0}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def integer(request, name, default, lower, upper, power_of_two=False):
    value = request.get(name, default)
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError("%s must be an integer in [%d,%d]" % (name, lower, upper))
    if power_of_two and value & (value - 1):
        raise ValueError(name + " must be a power of two")
    return value


def new_merge_stats():
    return {"leaf_emissions": 0, "internal_emissions": 0,
            "key_comparisons": 0, "child_end_markers": 0}


def lazy_merge(records, stats, reverse_comparison=False):
    """Each leaf has one record; each internal node retains at most two heads.

    The record format is (ordered_key, scalar_index, represented_group_point).
    None denotes end-of-stream only; the infinity record has key -1 instead.
    Equal keys are taken from the left child first, preserving duplicates.
    """
    size = len(records)
    require(size > 0 and size & (size - 1) == 0, "merge tree needs a power-of-two leaf count")

    def subtree(begin, end):
        if end - begin == 1:
            stats["leaf_emissions"] += 1
            yield records[begin]
            return
        middle = (begin + end) // 2
        left = subtree(begin, middle)
        right = subtree(middle, end)
        left_head = next(left, None)
        right_head = next(right, None)
        while left_head is not None or right_head is not None:
            if left_head is None:
                take_left = False
            elif right_head is None:
                take_left = True
            else:
                stats["key_comparisons"] += 1
                take_left = (left_head[0] >= right_head[0] if reverse_comparison
                             else left_head[0] <= right_head[0])
            stats["internal_emissions"] += 1
            if take_left:
                yield left_head
                left_head = next(left, None)
                stats["child_end_markers"] += int(left_head is None)
            else:
                yield right_head
                right_head = next(right, None)
                stats["child_end_markers"] += int(right_head is None)

    return subtree(0, size)


def inspect_merge(records):
    """Independent list sorting is an evidence oracle, not attack machinery."""
    stats = new_merge_stats()
    merged = list(lazy_merge(records, stats))
    expected = sorted(records, key=lambda record: record[0])
    require(merged == expected, "lazy merge lost, reordered, or duplicated a record")
    levels = len(records).bit_length() - 1
    require(stats["leaf_emissions"] == len(records), "wrong leaf-emission count")
    require(stats["internal_emissions"] == len(records) * levels,
            "record traversals differ from full binary-tree bound")
    require(stats["key_comparisons"] <= len(records) * levels,
            "merge comparison count exceeds bound")
    require(stats["child_end_markers"] == 2 * (len(records) - 1),
            "end markers were lost or double-counted")
    return stats


def expected_coverage(order, m, h):
    """Scalar-set oracle independent of merging, signs, and point arithmetic."""
    d = 2 * m - 1
    return {(i * d + offset) % order
            for i in range(h) for offset in range(1 - m, m)}


def recover(target, order, m, babies, giants, sign_of_match, verify, omit_minus=False):
    """Merge-intersect ascending streams, advancing giants after a failed check.

    Baby keys are unique. Giant keys need not be: retain the baby head after a
    failed mutant verification so all giant records sharing that key are seen.
    Normal verified recovery returns at the first match.
    """
    baby_stats, giant_stats = new_merge_stats(), new_merge_stats()
    baby_stream = lazy_merge(babies, baby_stats)
    giant_stream = lazy_merge(giants, giant_stats)
    baby = next(baby_stream, None)
    giant = next(giant_stream, None)
    comparisons = verifications = negative_matches = 0
    d = 2 * m - 1
    while baby is not None and giant is not None:
        comparisons += 1
        if baby[0] < giant[0]:
            baby = next(baby_stream, None)
        elif giant[0] < baby[0]:
            giant = next(giant_stream, None)
        else:
            sign = sign_of_match(giant[2], baby[2])
            negative_matches += int(sign == -1)
            candidate = (giant[1] * d + (baby[1] if omit_minus else sign * baby[1])) % order
            verifications += 1
            if verify(candidate, target):
                return candidate, comparisons, verifications, negative_matches
            giant = next(giant_stream, None)
    return None, comparisons, verifications, negative_matches


def toy_case(order, m, h):
    require(order % 2 == 1, "toy order must be odd")
    require(m <= (order + 1) // 2, "toy baby quotient keys must be unique")
    d = 2 * m - 1
    babies = [(min(j, (-j) % order), j, j) for j in range(m)]
    require(len({row[0] for row in babies}) == m, "duplicate baby key")
    baby_merge = inspect_merge(babies)
    oracle = expected_coverage(order, m, h)
    require(len(oracle) == min(order, h * d), "contiguous-interval count mismatch")

    def sign_of_match(giant, baby):
        sign = 1 if giant == baby else -1
        require(giant == (sign * baby) % order, "toy sign mismatch")
        return sign

    def verify(candidate, target):
        return candidate == target

    successes = failures = comparisons = verifications = negatives = 0
    mutant_misses = duplicate_giant_targets = 0
    first_mutant = None
    for secret in range(order):
        values = [(secret - i * d) % order for i in range(h)]
        giants = [(min(value, (-value) % order), i, value)
                  for i, value in enumerate(values)]
        duplicate_giant_targets += int(len({row[0] for row in giants}) < h)
        inspect_merge(giants)
        recovered, compared, verified, signed = recover(
            secret, order, m, babies, giants, sign_of_match, verify)
        require((recovered is not None) == (secret in oracle), "toy coverage mismatch")
        require(recovered is None or recovered == secret, "wrong toy output")
        mutant, _, _, _ = recover(secret, order, m, babies, giants,
                                  sign_of_match, verify, omit_minus=True)
        require(mutant is None or mutant == secret, "mutant escaped scalar verification")
        if recovered is not None and mutant is None:
            mutant_misses += 1
            if first_mutant is None:
                first_mutant = secret
        successes += int(recovered is not None)
        failures += int(recovered is None)
        comparisons += compared
        verifications += verified
        negatives += signed
    return {
        "order": order, "m": m, "h": h, "checked_targets": order,
        "oracle_covered": len(oracle), "successful_targets": successes,
        "uncovered_failures": failures, "wrong_outputs": 0,
        "intersection_comparisons": comparisons, "verification_attempts": verifications,
        "negative_matches": negatives, "targets_with_duplicate_giant_keys": duplicate_giant_targets,
        "minus_omission_mutant_false_negatives": mutant_misses,
        "first_mutant_witness": first_mutant, "full_baby_merge_events": baby_merge,
        "scope": "exhaustive for this odd-order additive toy group only",
    }


def limbs_of(value):
    return [(value >> (64 * k)) & MASK for k in range(LIMBS)]


def value_of(limbs):
    return sum(limb << (64 * k) for k, limb in enumerate(limbs))


def mont_mul(left, right):
    """64-bit limb schoolbook product followed by four Montgomery REDC rounds.

    Inputs and output are canonical integers below FIELD. Each limb product is
    one MUL/UMULH pair; every add of a 64-bit word with carry is tallied. The
    returned value is left*right*2^-256 mod FIELD.
    """
    require(0 <= left < FIELD and 0 <= right < FIELD, "noncanonical Montgomery operand")
    COUNTS["mont_mul"] += 1
    a, b = limbs_of(left), limbs_of(right)
    acc = [0] * (2 * LIMBS + 1)
    for i in range(LIMBS):
        carry = 0
        for j in range(LIMBS):
            COUNTS["limb_mul_pairs"] += 1
            word = acc[i + j] + a[i] * b[j] + carry
            acc[i + j] = word & MASK
            carry = word >> 64
            COUNTS["limb_add_steps"] += 2
        acc[i + LIMBS] = carry
    for i in range(LIMBS):
        q = (acc[i] * P_PRIME) & MASK
        COUNTS["limb_mul_pairs"] += 1
        carry = 0
        for j in range(LIMBS):
            COUNTS["limb_mul_pairs"] += 1
            word = acc[i + j] + q * FIELD_LIMBS[j] + carry
            acc[i + j] = word & MASK
            carry = word >> 64
            COUNTS["limb_add_steps"] += 2
        k = i + LIMBS
        while carry:
            require(k < len(acc), "REDC accumulator overflow")
            word = acc[k] + carry
            acc[k] = word & MASK
            carry = word >> 64
            k += 1
            COUNTS["limb_add_steps"] += 1
        require(acc[i] == 0, "REDC round did not clear its low limb")
    result = value_of(acc[LIMBS:])
    require(result < 2 * FIELD, "REDC output not below 2p")
    if result >= FIELD:
        result -= FIELD
        COUNTS["limb_add_steps"] += LIMBS + 1
    return result


def to_mont(value):
    return mont_mul(value % FIELD, R2_MOD_P)


def from_mont(value):
    return mont_mul(value, 1)


def fadd(left, right):
    return (left + right) % FIELD


def fsub(left, right):
    return (left - right) % FIELD


def mont_inverse(value):
    """value^(p-2) by left-to-right binary exponentiation in Montgomery form."""
    require(value != 0, "zero inversion")
    start = COUNTS["mont_mul"]
    result = ONE_MONT
    for bit in bin(FIELD - 2)[2:]:
        result = mont_mul(result, result)
        if bit == "1":
            result = mont_mul(result, value)
    used = COUNTS["mont_mul"] - start
    require(used <= 512, "inversion exceeded 512 multiplications")
    return result


def montgomery_checks():
    require(P_PRIME == 1, "P-256 has p = -1 mod 2^64, so p' = 1")
    require(FIELD == (1 << 256) - (1 << 224) + (1 << 192) + (1 << 96) - 1, "field shape")
    values = [0, 1, 2, 3, FIELD - 1, FIELD - 2, (1 << 255), (1 << 224) - 1,
              GENERATOR[0], GENERATOR[1], TWICE_GENERATOR[0], TWICE_GENERATOR[1], COEFF_B]
    state = 0x243F6A8885A308D313198A2E03707344A4093822299F31D0082EFA98EC4E6C89
    for _ in range(24):
        state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 256)
        values.append(state % FIELD)
    pairs = 0
    before = dict(COUNTS)
    for left in values:
        for right in values:
            got = mont_mul(left, right)
            require(got == left * right * pow(R_MONT, FIELD - 2, FIELD) % FIELD, "REDC mismatch")
            require(from_mont(mont_mul(to_mont(left), to_mont(right))) == left * right % FIELD,
                    "Montgomery round trip mismatch")
            pairs += 1
    calls = COUNTS["mont_mul"] - before["mont_mul"]
    limb_pairs = COUNTS["limb_mul_pairs"] - before["limb_mul_pairs"]
    adds = COUNTS["limb_add_steps"] - before["limb_add_steps"]
    require(limb_pairs == calls * 36, "limb product count is not 16 + 4*(1+4) per multiply")
    # 32 product adds + 32 REDC adds + at most 5+4+3+2 carry ripples + 5 subtract steps = 83.
    require(adds <= calls * 83, "carry step tally exceeded 83 per multiply")
    for value in values:
        inverse = mont_inverse(to_mont(value)) if value else None
        if inverse is not None:
            require(from_mont(mont_mul(inverse, to_mont(value))) == 1, "inverse failed")
    return {"checked_products": pairs, "p_prime": P_PRIME,
            "limb_mul_pairs_per_multiply": 36,
            "max_carry_steps_per_multiply_bound": 83,
            "scope": "exact REDC algebra on fixed and pseudorandom canonical operands; "
                     "tallies are word-operation events, not AArch64 instruction counts"}


def affine_negate(point):
    return None if point is None else (point[0], (-point[1]) % FIELD)


def affine_add(left, right):
    """Independent plain-integer affine oracle; no Montgomery, Jacobian or batch calls."""
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if (y + v) % FIELD == 0:
            return None
        require(y == v, "invalid affine equal-x case")
        numerator, denominator = 3 * x * x + COEFF_A, 2 * y
    else:
        numerator, denominator = v - y, u - x
    slope = numerator * pow(denominator % FIELD, FIELD - 2, FIELD) % FIELD
    result_x = (slope * slope - x - u) % FIELD
    return result_x, (slope * (x - result_x) - y) % FIELD


def affine_multiply(scalar, point):
    require(scalar >= 0, "negative scalar")
    result = None
    running = point
    while scalar:
        if scalar & 1:
            result = affine_add(result, running)
        running = affine_add(running, running)
        scalar >>= 1
    return result


def jacobian_double(point):
    """dbl-2001-b for a=-3, on Montgomery-form coordinates."""
    x, y, z = point
    if z == 0 or y == 0:
        return J_INFINITY
    delta = mont_mul(z, z)
    gamma = mont_mul(y, y)
    beta = mont_mul(x, gamma)
    t = mont_mul(fsub(x, delta), fadd(x, delta))
    alpha = fadd(fadd(t, t), t)
    beta4 = fadd(beta, beta)
    beta4 = fadd(beta4, beta4)
    beta8 = fadd(beta4, beta4)
    out_x = fsub(mont_mul(alpha, alpha), beta8)
    out_z = fsub(fsub(mont_mul(fadd(y, z), fadd(y, z)), gamma), delta)
    gamma2 = mont_mul(gamma, gamma)
    gamma8 = fadd(gamma2, gamma2)
    gamma8 = fadd(gamma8, gamma8)
    gamma8 = fadd(gamma8, gamma8)
    out_y = fsub(mont_mul(alpha, fsub(beta4, out_x)), gamma8)
    return out_x, out_y, out_z


def jacobian_mixed_add(point, affine):
    """Complete wrapper for madd-2004-hmv; all values in Montgomery form."""
    if affine is None:
        return point
    x, y, z = point
    u, v = affine
    if z == 0:
        return u, v, ONE_MONT
    z2 = mont_mul(z, z)
    u2 = mont_mul(u, z2)
    s2 = mont_mul(v, mont_mul(z, z2))
    h = fsub(u2, x)
    delta = fsub(s2, y)
    if h == 0:
        return jacobian_double(point) if delta == 0 else J_INFINITY
    hh = mont_mul(h, h)
    hhh = mont_mul(h, hh)
    v1 = mont_mul(x, hh)
    out_x = fsub(fsub(mont_mul(delta, delta), hhh), fadd(v1, v1))
    out_y = fsub(mont_mul(delta, fsub(v1, out_x)), mont_mul(y, hhh))
    out_z = mont_mul(z, h)
    return out_x, out_y, out_z


def counted(function, *args):
    start = COUNTS["mont_mul"]
    result = function(*args)
    return result, COUNTS["mont_mul"] - start


def block_normalize(points):
    """Batch inversion of Montgomery Z factors, replacing zero factors with one.

    Returns canonical plain affine coordinates (converted out of Montgomery form).
    """
    if not points:
        return []
    start = COUNTS["mont_mul"]
    factors = [point[2] if point[2] != 0 else ONE_MONT for point in points]
    prefixes = []
    product = ONE_MONT
    for factor in factors:
        prefixes.append(product)
        product = mont_mul(product, factor)
    inverse_start = COUNTS["mont_mul"]
    inverse = mont_inverse(product)
    inversion_cost = COUNTS["mont_mul"] - inverse_start
    inverses = [0] * len(points)
    for index in range(len(points) - 1, -1, -1):
        inverses[index] = mont_mul(inverse, prefixes[index])
        inverse = mont_mul(inverse, factors[index])
    normalized = []
    for (x, y, z), zi in zip(points, inverses):
        if z == 0:
            normalized.append(None)
        else:
            zi2 = mont_mul(zi, zi)
            zi3 = mont_mul(zi2, zi)
            normalized.append((from_mont(mont_mul(x, zi2)), from_mont(mont_mul(y, zi3))))
    used = COUNTS["mont_mul"] - start - inversion_cost
    require(used <= 9 * len(points), "block normalization exceeded 9 multiplications per point")
    return normalized


def jacobian_of(point, scale=1):
    if point is None:
        return J_INFINITY
    s = to_mont(scale)
    s2 = mont_mul(s, s)
    return (mont_mul(to_mont(point[0]), s2),
            mont_mul(to_mont(point[1]), mont_mul(s2, s)), s)


def mont_affine(point):
    return None if point is None else (to_mont(point[0]), to_mont(point[1]))


def progression(start, step, count, block_size):
    """Generate projectively, then normalize blocks; indices stay in order."""
    current = jacobian_of(start)
    increment = mont_affine(step)
    output = []
    max_mults = 0
    while len(output) < count:
        block = []
        for _ in range(min(block_size, count - len(output))):
            block.append(current)
            current, used = counted(jacobian_mixed_add, current, increment)
            max_mults = max(max_mults, used)
        output.extend(block_normalize(block))
    require(max_mults <= 16, "a Jacobian step exceeded 16 multiplications")
    return output


def affine_progression(start, step, count):
    output = []
    current = start
    for _ in range(count):
        output.append(current)
        current = affine_add(current, step)
    return output


def arithmetic_checks():
    require((GENERATOR[1] ** 2 - GENERATOR[0] ** 3 - COEFF_A * GENERATOR[0] - COEFF_B) % FIELD == 0,
            "protected generator is off curve")
    require(COEFF_A == FIELD - 3, "profile coefficient a is not -3")
    require(affine_add(GENERATOR, GENERATOR) == TWICE_GENERATOR, "public 2P fixture failed")
    require(affine_multiply(ORDER, GENERATOR) is None, "rP is not infinity")
    montgomery = montgomery_checks()
    points = [None, GENERATOR, TWICE_GENERATOR, affine_negate(GENERATOR),
              affine_multiply(3, GENERATOR), affine_multiply(17, GENERATOR)]
    encodings = []
    expected = []
    mixed_cases = doubling_cases = 0
    max_mixed = max_double = 0
    for point in points:
        for scale in (1, 3, FIELD - 1):
            projective = jacobian_of(point, scale)
            encodings.append(projective)
            expected.append(point)
            doubled_point, used = counted(jacobian_double, projective)
            max_double = max(max_double, used)
            require(block_normalize([doubled_point])[0] == affine_add(point, point),
                    "Jacobian a=-3 doubling differs from affine")
            doubling_cases += 1
            for added in points:
                summed, used = counted(jacobian_mixed_add, projective, mont_affine(added))
                max_mixed = max(max_mixed, used)
                require(block_normalize([summed])[0] == affine_add(point, added),
                        "mixed addition differs from affine")
                mixed_cases += 1
    require(max_double <= 8, "doubling used more than 8 multiplications")
    require(max_mixed <= 12, "mixed branch used more than 4+8 multiplications")
    require(block_normalize(encodings) == expected, "mixed zero/nonzero-Z block failed")
    require(block_normalize([J_INFINITY] * 8) == [None] * 8, "all-infinity block failed")
    require(block_normalize([]) == [], "empty block failed")
    # The y=0 branch is unreachable for a nonidentity point in this odd-order
    # production group; this synthetic tuple checks the defensive branch only.
    require(jacobian_double((ONE_MONT, 0, ONE_MONT)) == J_INFINITY, "zero-Y defensive branch failed")
    return {"mixed_addition_cases": mixed_cases, "doubling_cases": doubling_cases,
            "max_doubling_multiplications": max_double,
            "max_mixed_branch_multiplications": max_mixed,
            "mixed_Z_block_size": len(encodings), "all_infinity_block_size": 8,
            "rescaling_factors": ["1", "3", "p-1"],
            "synthetic_zero_Y_branch_only": True,
            "montgomery": montgomery,
            "scope": "bounded full-field formula checks, not an instruction-cost measurement"}


def curve_key(point):
    return -1 if point is None else point[0]


def curve_checks(request):
    m = integer(request, "m", 16, 2, 64, power_of_two=True)
    h = integer(request, "h", 16, 1, 64, power_of_two=True)
    block_size = integer(request, "block_size", 4, 1, 16)
    arithmetic = arithmetic_checks()
    d = 2 * m - 1
    step = affine_negate(affine_multiply(d, GENERATOR))
    baby_points = progression(None, GENERATOR, m, block_size)
    require(baby_points == affine_progression(None, GENERATOR, m), "baby progression mismatch")
    babies = [(curve_key(point), j, point) for j, point in enumerate(baby_points)]
    require(len({row[0] for row in babies}) == m, "curve baby keys are not unique")
    baby_merge = inspect_merge(babies)
    oracle = expected_coverage(ORDER, m, h)
    require(len(oracle) == h * d, "small full-field coverage count mismatch")
    low, high = 1 - m, h * d - m
    raw = {0, 1, 2, m - 1, m, d - (m - 1), d, d + (m - 1),
           (h - 1) * d - (m - 1), (h - 1) * d,
           low - 1, low, low + 1, high - 1, high, high + 1,
           ORDER - 1, ORDER // 2}
    secrets = sorted({value % ORDER for value in raw})

    def sign_of_match(giant, baby):
        if giant is None:
            require(baby is None, "infinity matched a finite point")
            return 1
        require(baby is not None and giant[0] == baby[0], "invalid curve key match")
        if giant[1] == baby[1]:
            return 1
        require(giant[1] == (-baby[1]) % FIELD, "same-x points are not opposite")
        return -1

    def verify(candidate, target):
        return affine_multiply(candidate, GENERATOR) == target

    outcomes = []
    successes = failures = mutants = duplicate_giant_targets = 0
    for secret in secrets:
        target = affine_multiply(secret, GENERATOR)
        giant_points = progression(target, step, h, block_size)
        require(giant_points == affine_progression(target, step, h), "giant progression mismatch")
        giants = [(curve_key(point), i, point) for i, point in enumerate(giant_points)]
        duplicate_giant_targets += int(len({row[0] for row in giants}) < h)
        inspect_merge(giants)
        recovered, _, _, _ = recover(target, ORDER, m, babies, giants, sign_of_match, verify)
        require((recovered is not None) == (secret in oracle), "full-field coverage mismatch")
        require(recovered is None or recovered == secret, "wrong full-field scalar")
        mutant, _, _, _ = recover(target, ORDER, m, babies, giants, sign_of_match,
                                  verify, omit_minus=True)
        require(mutant is None or mutant == secret, "full-field mutant escaped verification")
        missed = recovered is not None and mutant is None
        mutants += int(missed)
        successes += int(recovered is not None)
        failures += int(recovered is None)
        outcomes.append({"scalar_hex": format(secret, "064x"), "covered": secret in oracle,
                         "recovered": recovered is not None, "mutant_false_negative": missed})
    require(mutants > 0, "full-field fixtures did not detect omitted minus sign")
    require(successes > 0 and failures > 0, "both recovery and failure need coverage")
    if h >= 4:
        require(duplicate_giant_targets > 0, "full-field giant-duplicate fixture missing")
    return {
        "group": "protected production P-256 (a=-3), small search bounds",
        "m": m, "h": h, "block_size": block_size, "checked_targets": len(secrets),
        "successful_targets": successes, "uncovered_failures": failures, "wrong_outputs": 0,
        "targets_with_duplicate_giant_keys": duplicate_giant_targets,
        "minus_omission_mutant_false_negatives": mutants,
        "arithmetic": arithmetic, "full_baby_merge_events": baby_merge, "outcomes": outcomes,
        "scope": "targeted full-field fixtures; not a production success-frequency or cost estimate",
    }


def merge_mutation_check():
    records = [(3, 0, None), (0, 1, None), (2, 2, None), (1, 3, None)]
    correct = list(lazy_merge(records, new_merge_stats()))
    mutant = list(lazy_merge(records, new_merge_stats(), reverse_comparison=True))
    expected = sorted(records, key=lambda row: row[0])
    require(correct == expected and mutant != expected, "merge-order mutation was not detected")
    duplicate_records = [(1, 0, None), (0, 1, None), (1, 2, None), (0, 3, None)]
    inspect_merge(duplicate_records)
    inspect_merge([(0, 0, None)])
    # In Z/17 with d=3, target16 gives giant0=-1 and giant5=+1.
    # Both keys equal1. The omit-minus mutant rejects the first but must retain
    # baby1 long enough to verify the second. There is no infinity giant here.
    babies = [(0, 0, 0), (1, 1, 1)]
    giants = [(min((16 - 3 * i) % 17, (-16 + 3 * i) % 17),
               i, (16 - 3 * i) % 17) for i in range(8)]
    recovered, _, attempts, _ = recover(
        16, 17, 2, babies, giants,
        lambda giant, baby: 1 if giant == baby else -1,
        lambda candidate, target: candidate == target, omit_minus=True)
    require(recovered == 16 and attempts == 2, "duplicate-key retry dropped the baby head")
    return {"reversed_comparison_detected": True, "duplicate_keys_preserved": True,
            "singleton_stream_checked": True, "end_markers_checked": True,
            "duplicate_key_retry_after_failed_verification_checked": True}


def accounting_checks():
    m = h = 1 << 127
    d = 2 * m - 1
    covered = h * d
    require(covered < ORDER and 2 * covered > ORDER, "production half-coverage failed")
    require(100 * covered >= 39 * ORDER, "production threshold failed")
    terms = {
        "fixed_setup_validation_output": 1 << 40,
        "all_tree_initialization": m * 4 * (1 << 18),
        "both_progressions_normalization_and_placement": m * 2 * ((1 << 18) + (1 << 17) + 4096),
        "both_complete_merge_streams": m * 2 * (1 << 18),
        "root_intersection": m * (1 << 14),
    }
    total = sum(terms.values())
    memory = 4 * m * (1 << 16) + (1 << 24)
    # Raising to the 10000th power checks score 1382191/10000 without floating logarithms.
    require(total ** 10000 < 1000 ** 10000 * (1 << 1382191), "submitted 1382191/10000 bound failed")
    require(memory < (1 << 146), "memory bound failed")
    require(2 * covered - ORDER == 26959946660873538059280334322932747074583865233571950339412849384111,
            "stated 2*COV-r integer mismatch")
    require(2 * m - 2 < ORDER and covered - m < ORDER, "index range does not fit below r")
    require(not total ** 10000 < 1000 ** 10000 * (1 << 1382190), "1382190/10000 would also hold; score not tight")
    return {
        "m": str(m), "h": str(h), "d": str(d), "r": str(ORDER),
        "coverage_count": str(covered), "twice_coverage_minus_order": str(2 * covered - ORDER),
        "success_at_least_half": True,
        "success_at_least_39_percent": True,
        "instruction_terms": {name: str(value) for name, value in terms.items()},
        "total_instructions": str(total), "reference_price": 1000,
        "score_numerator": 1382191, "score_denominator": 10000,
        "exact_10000_power_score_inequality": True,
        "preprocessing_bound": "at most total work, hence at most 2^(1382191/10000) equivalents",
        "peak_memory_bytes_bound": str(memory), "memory_below_2_to_146": True,
        "scope": "arithmetic checks of proposed ledger; does not establish per-kernel or communication costs",
    }


def main():
    request = json.load(sys.stdin)
    allowed = {"mode", "order", "m", "h", "block_size"}
    if not isinstance(request, dict) or set(request) - allowed:
        raise ValueError("expected object with only mode/order/m/h/block_size")
    mode = request.get("mode", "all")
    if mode not in ("all", "toy", "curve", "accounting"):
        raise ValueError("unknown mode")
    result = {"mode": mode, "heuristics": [],
              "event_count_scope": "bounded fixture behavior, not production instruction costs"}
    if mode == "all":
        configurations = [(17, 2, 1), (17, 2, 8), (17, 4, 4),
                          (101, 8, 4), (101, 8, 8), (121, 8, 4),
                          (257, 16, 8), (257, 16, 16)]
        result["toy"] = [toy_case(*parameters) for parameters in configurations]
        require(any(item["minus_omission_mutant_false_negatives"] > 0
                    for item in result["toy"]), "toy minus mutation was not detected")
        require(any(item["targets_with_duplicate_giant_keys"] > 0
                    for item in result["toy"]), "toy giant duplicates were not exercised")
        result["merge_controls"] = merge_mutation_check()
        result["curve"] = curve_checks(request)
        result["accounting"] = accounting_checks()
    elif mode == "toy":
        order = integer(request, "order", 101, 17, 257)
        m = integer(request, "m", 8, 2, 64, power_of_two=True)
        h = integer(request, "h", 8, 1, 64, power_of_two=True)
        result["toy"] = toy_case(order, m, h)
        witness = toy_case(17, 2, 1)
        require(witness["minus_omission_mutant_false_negatives"] > 0, "minus mutation not detected")
        result["mandatory_minus_mutation_witness"] = witness
        result["merge_controls"] = merge_mutation_check()
    elif mode == "curve":
        result["curve"] = curve_checks(request)
        result["merge_controls"] = merge_mutation_check()
    else:
        result["accounting"] = accounting_checks()
    result["all_requested_checks_passed"] = True
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
