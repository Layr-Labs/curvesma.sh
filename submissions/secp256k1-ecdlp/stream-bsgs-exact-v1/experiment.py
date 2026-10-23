"""Bounded stream-join/Jacobian evidence; execute only in the shared verifier.

One JSON request on stdin, one JSON result on stdout. The fixture harness knows
test scalars; recover() receives only group targets and generated records.
No event counter here is an AArch64 instruction count or production benchmark.
"""

import json
import sys


FIELD = int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
ORDER = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
GENERATOR = (
    int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
    int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16),
)
# Public profile fixture, reused from the repository's production-v1 vectors.
TWICE_GENERATOR = (
    int("c6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5", 16),
    int("1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52a", 16),
)
J_INFINITY = (0, 1, 0)


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


def affine_negate(point):
    return None if point is None else (point[0], (-point[1]) % FIELD)


def affine_add(left, right):
    """Independent affine oracle; no Jacobian or batch-inversion calls."""
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
        numerator, denominator = 3 * x * x, 2 * y
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
    x, y, z = point
    if z == 0 or y == 0:
        return J_INFINITY
    a = x * x % FIELD
    b = y * y % FIELD
    c = b * b % FIELD
    d = 2 * ((x + b) ** 2 - a - c) % FIELD
    e = 3 * a % FIELD
    f = e * e % FIELD
    out_x = (f - 2 * d) % FIELD
    out_y = (e * (d - out_x) - 8 * c) % FIELD
    out_z = 2 * y * z % FIELD
    return out_x, out_y, out_z


def jacobian_mixed_add(point, affine):
    """Complete wrapper for mixed Jacobian/affine addition on y^2=x^3+7."""
    if affine is None:
        return point
    x, y, z = point
    u, v = affine
    if z == 0:
        return u, v, 1
    z2 = z * z % FIELD
    u2 = u * z2 % FIELD
    s2 = v * z * z2 % FIELD
    h = (u2 - x) % FIELD
    delta = (s2 - y) % FIELD
    if h == 0:
        return jacobian_double(point) if delta == 0 else J_INFINITY
    hh = h * h % FIELD
    hhh = h * hh % FIELD
    v1 = x * hh % FIELD
    out_x = (delta * delta - hhh - 2 * v1) % FIELD
    out_y = (delta * (v1 - out_x) - y * hhh) % FIELD
    out_z = z * h % FIELD
    return out_x, out_y, out_z


def block_normalize(points):
    """Montgomery inversion of Z factors, replacing zero factors with one."""
    if not points:
        return []
    factors = [point[2] if point[2] != 0 else 1 for point in points]
    prefixes = []
    product = 1
    for factor in factors:
        prefixes.append(product)
        product = product * factor % FIELD
    inverse = pow(product, FIELD - 2, FIELD)
    inverses = [0] * len(points)
    for index in range(len(points) - 1, -1, -1):
        inverses[index] = inverse * prefixes[index] % FIELD
        inverse = inverse * factors[index] % FIELD
    normalized = []
    for (x, y, z), zi in zip(points, inverses):
        if z == 0:
            normalized.append(None)
        else:
            zi2 = zi * zi % FIELD
            normalized.append((x * zi2 % FIELD, y * zi2 * zi % FIELD))
    return normalized


def progression(start, step, count, block_size):
    """Generate projectively, then normalize blocks; indices stay in order."""
    current = J_INFINITY if start is None else (start[0], start[1], 1)
    output = []
    while len(output) < count:
        block = []
        for _ in range(min(block_size, count - len(output))):
            block.append(current)
            current = jacobian_mixed_add(current, step)
        output.extend(block_normalize(block))
    return output


def affine_progression(start, step, count):
    output = []
    current = start
    for _ in range(count):
        output.append(current)
        current = affine_add(current, step)
    return output


def arithmetic_checks():
    require((GENERATOR[1] ** 2 - GENERATOR[0] ** 3 - 7) % FIELD == 0,
            "protected generator is off curve")
    require(affine_add(GENERATOR, GENERATOR) == TWICE_GENERATOR, "public 2P fixture failed")
    require(affine_multiply(ORDER, GENERATOR) is None, "rP is not infinity")
    points = [None, GENERATOR, TWICE_GENERATOR, affine_negate(GENERATOR),
              affine_multiply(3, GENERATOR), affine_multiply(17, GENERATOR)]
    encodings = []
    expected = []
    mixed_cases = doubling_cases = 0
    for point in points:
        for scale in (1, 3, FIELD - 1):
            if point is None:
                projective = J_INFINITY
            else:
                projective = (point[0] * scale * scale % FIELD,
                              point[1] * scale * scale * scale % FIELD, scale)
            encodings.append(projective)
            expected.append(point)
            doubled = block_normalize([jacobian_double(projective)])[0]
            require(doubled == affine_add(point, point), "Jacobian doubling differs from affine")
            doubling_cases += 1
            for added in points:
                result = block_normalize([jacobian_mixed_add(projective, added)])[0]
                require(result == affine_add(point, added), "mixed addition differs from affine")
                mixed_cases += 1
    require(block_normalize(encodings) == expected, "mixed zero/nonzero-Z block failed")
    require(block_normalize([J_INFINITY] * 8) == [None] * 8, "all-infinity block failed")
    require(block_normalize([]) == [], "empty block failed")
    # The y=0 branch is unreachable for a nonidentity point in this odd-order
    # production group; this synthetic tuple checks the defensive branch only.
    require(jacobian_double((1, 0, 1)) == J_INFINITY, "zero-Y defensive branch failed")
    return {"mixed_addition_cases": mixed_cases, "doubling_cases": doubling_cases,
            "mixed_Z_block_size": len(encodings), "all_infinity_block_size": 8,
            "rescaling_factors": ["1", "3", "p-1"],
            "synthetic_zero_Y_branch_only": True,
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
        "group": "protected production secp256k1, small search bounds",
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
    return {
        "m": str(m), "h": str(h), "d": str(d), "r": str(ORDER),
        "coverage_count": str(covered), "success_at_least_half": True,
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
