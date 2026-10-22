"""Bounded signed-BSGS checks for the shared, credential-free verifier only.

Read one JSON object on stdin; write one JSON object on stdout. Modes are
all, toy, curve, accounting. This is correctness evidence, not the production
implementation or a measurement of the analytic instruction/memory ledger.
The test harness knows fixture scalars; recover_* receive only group targets.
"""

import bisect
import json
import sys
from typing import Optional, Tuple


FIELD = int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16)
ORDER = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
GENERATOR = (
    int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
    int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16),
)
# Public 2P fixture from tests/fixtures/production-v1-vectors.json.
TWICE_GENERATOR = (
    int("c6047f9441ed7d6d3045406e95c07cd85c778e4b8cef3ca7abac09b95c709ee5", 16),
    int("1ae168fea63dc339a3c58419466ceaeef7f632653266d0e1236431a950cfe52a", 16),
)
Point = Optional[Tuple[int, int]]


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def bounded_integer(request, name, default, lower, upper):
    value = request.get(name, default)
    if type(value) is not int or not lower <= value <= upper:
        raise ValueError("%s must be an integer in [%d,%d]" % (name, lower, upper))
    return value


def expected_coverage(order, m, h):
    """Independent scalar-set oracle; not used by either recovery algorithm."""
    d = 2 * m - 1
    return {(i * d + offset) % order
            for i in range(h) for offset in range(1 - m, m)}


def recover_toy(target, order, m, h, table, omit_minus=False):
    keys, records = table
    d = 2 * m - 1
    gamma = target
    probes = verifications = negative_matches = 0
    for i in range(h):
        probes += 1
        key = min(gamma, (-gamma) % order)
        position = bisect.bisect_left(keys, key)
        if position < len(keys) and keys[position] == key:
            j = records[position]
            sign = 1 if gamma == j else -1
            require(gamma == (sign * j) % order, "invalid toy quotient match")
            negative_matches += int(sign == -1)
            candidate = (i * d + (j if omit_minus else sign * j)) % order
            verifications += 1
            # Independent attack-side verification rejects mutant false outputs.
            if candidate == target:
                return candidate, probes, verifications, negative_matches
        gamma = (gamma - d) % order
    return None, probes, verifications, negative_matches


def toy_case(order, m, h):
    require(order % 2 == 1, "toy order must be odd")
    require(m <= (order + 1) // 2, "baby quotient representatives must be unique")
    pairs = sorted((min(j, (-j) % order), j) for j in range(m))
    table = ([key for key, _ in pairs], [j for _, j in pairs])
    require(len(set(table[0])) == m, "duplicate toy baby key")
    oracle = expected_coverage(order, m, h)
    require(len(oracle) == min(order, h * (2 * m - 1)), "toy interval count mismatch")
    successes = failures = probes = verifications = negative_matches = 0
    mutant_misses = 0
    first_mutant_witness = None
    for secret in range(order):
        recovered, count, verified, negatives = recover_toy(secret, order, m, h, table)
        require((recovered is not None) == (secret in oracle), "toy coverage mismatch")
        require(recovered is None or recovered == secret, "incorrect toy scalar")
        probes += count
        verifications += verified
        negative_matches += negatives
        successes += int(recovered is not None)
        failures += int(recovered is None)
        mutant, _, _, _ = recover_toy(secret, order, m, h, table, omit_minus=True)
        require(mutant is None or mutant == secret, "mutant escaped verification")
        if recovered is not None and mutant is None:
            mutant_misses += 1
            if first_mutant_witness is None:
                first_mutant_witness = secret
    return {
        "group": "additive integers modulo odd order",
        "order": order, "m": m, "h": h, "d": 2 * m - 1,
        "checked_targets": order, "oracle_covered": len(oracle),
        "successful_targets": successes, "uncovered_failures": failures,
        "wrong_outputs": 0, "giant_probes": probes,
        "verification_attempts": verifications, "negative_matches": negative_matches,
        "minus_omission_mutant_false_negatives": mutant_misses,
        "first_mutant_witness": first_mutant_witness,
        "scope": "exhaustive for this toy group only; no production-cost inference",
    }


def toy_checks(request):
    order = bounded_integer(request, "order", 101, 17, 257)
    m = bounded_integer(request, "m", 8, 2, 64)
    h = bounded_integer(request, "h", 8, 1, 64)
    configured = toy_case(order, m, h)
    witness = toy_case(17, 2, 1)
    require(witness["minus_omission_mutant_false_negatives"] > 0,
            "negative-control test did not detect omitted minus sign")
    return {"configured": configured, "mandatory_mutation_witness": witness}


def negate(point: Point) -> Point:
    return None if point is None else (point[0], (-point[1]) % FIELD)


def add(left: Point, right: Point) -> Point:
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if (y + v) % FIELD == 0:
            return None
        require(y == v, "invalid equal-x curve points")
        numerator, denominator = 3 * x * x, 2 * y
    else:
        numerator, denominator = v - y, u - x
    slope = numerator * pow(denominator % FIELD, FIELD - 2, FIELD) % FIELD
    result_x = (slope * slope - x - u) % FIELD
    return result_x, (slope * (x - result_x) - y) % FIELD


def multiply(scalar: int, point: Point) -> Point:
    require(scalar >= 0, "negative scalar multiplication request")
    result = None
    running = point
    while scalar:
        if scalar & 1:
            result = add(result, running)
        running = add(running, running)
        scalar >>= 1
    return result


def curve_key(point: Point) -> int:
    return -1 if point is None else point[0]


def make_curve_table(m):
    rows = []
    running = None
    for j in range(m):
        rows.append((curve_key(running), j, None if running is None else running[1]))
        running = add(running, GENERATOR)
    rows.sort()
    keys = [row[0] for row in rows]
    require(len(set(keys)) == m, "duplicate full-field baby quotient key")
    return keys, rows


def recover_curve(target, m, h, table, step, omit_minus=False):
    keys, rows = table
    gamma = target
    d = 2 * m - 1
    probes = verifications = negative_matches = 0
    for i in range(h):
        probes += 1
        key = curve_key(gamma)
        position = bisect.bisect_left(keys, key)
        if position < len(keys) and keys[position] == key:
            _, j, baby_y = rows[position]
            if gamma is None:
                require(j == 0 and baby_y is None, "infinity table mismatch")
                sign = 1
            elif gamma[1] == baby_y:
                sign = 1
            else:
                require(baby_y is not None and gamma[1] == (-baby_y) % FIELD,
                        "same-x match is not a point or its negative")
                sign = -1
            negative_matches += int(sign == -1)
            candidate = (i * d + (j if omit_minus else sign * j)) % ORDER
            verifications += 1
            if multiply(candidate, GENERATOR) == target:
                return candidate, probes, verifications, negative_matches
        gamma = add(gamma, step)
    return None, probes, verifications, negative_matches


def curve_checks(request):
    m = bounded_integer(request, "m", 16, 2, 64)
    h = bounded_integer(request, "h", 16, 1, 64)
    d = 2 * m - 1
    require((GENERATOR[1] ** 2 - GENERATOR[0] ** 3 - 7) % FIELD == 0,
            "generator is off curve")
    require(add(GENERATOR, GENERATOR) == TWICE_GENERATOR, "public 2P fixture mismatch")
    require(multiply(ORDER, GENERATOR) is None, "generator order mismatch")
    require(add(GENERATOR, negate(GENERATOR)) is None, "opposite-point mismatch")
    table = make_curve_table(m)
    step = negate(multiply(d, GENERATOR))
    oracle = expected_coverage(ORDER, m, h)
    require(len(oracle) == h * d, "small full-field interval coverage mismatch")
    low, high = 1 - m, h * d - m
    # Both ends, block boundaries, infinity, both signs, and outside the interval.
    raw_scalars = {0, 1, 2, m - 1, m, d - (m - 1), d, d + (m - 1),
                   (h - 1) * d - (m - 1), (h - 1) * d,
                   low - 1, low, low + 1, high - 1, high, high + 1,
                   ORDER - 1, ORDER // 2}
    secrets = sorted(value % ORDER for value in raw_scalars)
    secrets = sorted(set(secrets))
    successes = failures = probes = verifications = negative_matches = 0
    mutant_misses = 0
    outcomes = []
    for secret in secrets:
        target = multiply(secret, GENERATOR)
        recovered, count, verified, negatives = recover_curve(target, m, h, table, step)
        require((recovered is not None) == (secret in oracle), "curve coverage mismatch")
        require(recovered is None or recovered == secret, "incorrect curve scalar")
        mutant, _, _, _ = recover_curve(target, m, h, table, step, omit_minus=True)
        require(mutant is None or mutant == secret, "curve mutant escaped verification")
        missed = recovered is not None and mutant is None
        mutant_misses += int(missed)
        successes += int(recovered is not None)
        failures += int(recovered is None)
        probes += count
        verifications += verified
        negative_matches += negatives
        outcomes.append({"scalar_hex": format(secret, "064x"),
                         "covered": secret in oracle,
                         "recovered": recovered is not None,
                         "mutant_false_negative": missed})
    require(mutant_misses > 0, "full-field negative control missed omitted minus sign")
    require(successes > 0 and failures > 0, "curve fixtures lack both outcomes")
    return {
        "group": "protected production secp256k1 with small table and search bounds",
        "m": m, "h": h, "d": d, "checked_targets": len(secrets),
        "oracle_covered_scalars": len(oracle), "successful_targets": successes,
        "uncovered_failures": failures, "wrong_outputs": 0,
        "giant_probes": probes, "verification_attempts": verifications,
        "negative_matches": negative_matches,
        "minus_omission_mutant_false_negatives": mutant_misses,
        "public_fixture_checks": ["P on curve", "2P", "rP=infinity", "P+(-P)=infinity"],
        "outcomes": outcomes,
        "scope": "full-size field and group; deliberately small search only; targeted fixtures are not a success-frequency or cost estimate",
    }


def accounting_checks():
    """Recompute the proof's exact integers, not the validity of kernel budgets."""
    m = h = 1 << 127
    d = 2 * m - 1
    covered = h * d
    require(covered < ORDER, "production coverage would overlap modulo r")
    require(2 * covered >= ORDER, "production success is below one half")
    require(100 * covered >= 39 * ORDER, "production success is below 39 percent")
    components = {
        "fixed_setup_validation_output": 1 << 40,
        "storage_construction": (1 << 22) * m,
        "baby_generation_and_storage": ((1 << 23) + (1 << 17) + (1 << 12)) * m,
        "merge_sort": 127 * ((1 << 20) + (1 << 12)) * m + (1 << 30),
        "giant_search": ((1 << 23) + 128 * ((1 << 17) + (1 << 12)) + (1 << 16)) * h,
    }
    total = sum(components.values())
    preprocessing = total - components["giant_search"]
    memory = (1 << 18) * m + (1 << 24)
    # Pinned totals detect transcription drift against the accompanying proof.
    require(total == 29292028818269973496407165951123623938690646016,
            "total ledger differs from proof")
    require(preprocessing == 24909932386758706049095397310548220902765494272,
            "preprocessing ledger differs from proof")
    require(memory == 44601490397061246283071436545296723028738048,
            "memory ledger differs from proof")
    require(total < 1000 * (1 << 145), "total exceeds claimed work bound")
    require(preprocessing < 1000 * (1 << 145), "preprocessing exceeds stated bound")
    require(memory < (1 << 146), "memory exceeds stated bound")
    return {
        "m": str(m), "h": str(h), "d": str(d), "r": str(ORDER),
        "covered_scalar_count": str(covered),
        "success_exact_numerator": str(covered), "success_exact_denominator": str(ORDER),
        "success_at_least_half": True, "success_at_least_39_percent": True,
        "instruction_components": {name: str(value) for name, value in components.items()},
        "total_instructions": str(total), "preprocessing_instructions": str(preprocessing),
        "reference_instructions_per_group_equivalent": 1000,
        "total_below_2_to_145_group_equivalents": True,
        "preprocessing_below_2_to_145_group_equivalents": True,
        "peak_memory_bytes_bound": str(memory), "peak_memory_below_2_to_146": True,
        "nonuniform_advice_bytes": 0,
        "scope": "exact arithmetic and threshold checks only; instruction and storage budgets require the analytic argument and are not measured here",
    }


def main():
    request = json.load(sys.stdin)
    if not isinstance(request, dict) or set(request) - {"mode", "order", "m", "h"}:
        raise ValueError("expected an object with only mode/order/m/h")
    mode = request.get("mode", "all")
    if mode not in ("all", "toy", "curve", "accounting"):
        raise ValueError("unsupported mode")
    result = {"mode": mode, "heuristics": [],
              "runtime_counter_units": "fixture-algorithm events, not charged AArch64 instructions"}
    if mode == "all":
        configurations = [(17, 2, 1), (17, 2, 8), (17, 4, 3),
                          (101, 8, 3), (101, 8, 7), (121, 7, 3),
                          (257, 16, 8), (257, 16, 9)]
        result["toy"] = [toy_case(*parameters) for parameters in configurations]
        require(any(case["minus_omission_mutant_false_negatives"] > 0
                    for case in result["toy"]), "toy mutation control failed")
        result["curve"] = curve_checks(request)
        result["accounting"] = accounting_checks()
    elif mode == "toy":
        result["toy"] = toy_checks(request)
    elif mode == "curve":
        result["curve"] = curve_checks(request)
    else:
        result["accounting"] = accounting_checks()
    result["all_requested_checks_passed"] = True
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
