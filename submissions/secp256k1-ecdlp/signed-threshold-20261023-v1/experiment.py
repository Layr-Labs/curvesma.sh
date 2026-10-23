"""Bounded specification evidence; run only in the shared pinned verifier.

This is NOT a production attack, AArch64 instruction trace, distributed router
implementation, or empirical production-success estimate. Python big integers
are used for finite checks only and are not assigned unit cost by the argument.
The predecessor's published affine batching/radix/fold construction is restated
here; no predecessor execution result is imported as this candidate's evidence.
"""

import json
import math
import sys
from collections import Counter


P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
R = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
GENERATOR = (
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)


class Curve:
    def __init__(self, p, a, b, generator, order):
        self.p, self.a, self.b = p, a, b
        self.generator, self.order = generator, order
        self.bits = p.bit_length()
        assert self.on_curve(generator)

    def on_curve(self, point):
        if point is None:
            return True
        x, y = point
        return 0 <= x < self.p and 0 <= y < self.p and (
            y * y - x * x * x - self.a * x - self.b
        ) % self.p == 0

    def neg(self, point):
        return None if point is None else (point[0], -point[1] % self.p)

    def add(self, left, right):
        if left is None:
            return right
        if right is None:
            return left
        x1, y1 = left
        x2, y2 = right
        if x1 == x2:
            if (y1 + y2) % self.p == 0:
                return None
            slope = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, self.p) % self.p
        else:
            slope = (y2 - y1) * pow(x2 - x1, -1, self.p) % self.p
        x3 = (slope * slope - x1 - x2) % self.p
        return (x3, (slope * (x1 - x3) - y1) % self.p)

    def multiply(self, scalar, point=None):
        if point is None:
            point = self.generator
        if scalar < 0:
            return self.multiply(-scalar, self.neg(point))
        result = None
        while scalar:
            if scalar & 1:
                result = self.add(result, point)
            point = self.add(point, point)
            scalar >>= 1
        return result


def record(point, label):
    # Native layout: x[4 words], tag[1 word], index[2 words] = 56 bytes.
    # Toy fields shorten only the test key width, not the production ledger.
    return (0, 0, label) if point is None else (point[0], 2 | (point[1] & 1), label)


def orbit_key(rec, field_bits):
    return 0 if rec[1] == 0 else (1 << field_bits) | rec[0]


def radix_sort(records, field_bits, digit_bits):
    """Actual stable LSD counting/scatter sort; bounded small test histograms."""
    assert 1 <= digit_bits <= 8
    radix = 1 << digit_bits
    result = list(records)
    for shift in range(0, field_bits + 1, digit_bits):
        counts = [0] * radix
        for rec in result:
            counts[(orbit_key(rec, field_bits) >> shift) & (radix - 1)] += 1
        cursors = []
        total = 0
        for count in counts:
            cursors.append(total)
            total += count
        destination = [None] * len(result)
        for rec in result:
            digit = (orbit_key(rec, field_bits) >> shift) & (radix - 1)
            destination[cursors[digit]] = rec
            cursors[digit] += 1
        assert all(rec is not None for rec in destination)
        result = destination
    assert result == sorted(records, key=lambda rec: orbit_key(rec, field_bits))
    return result


def batch_advance(curve, states, step, active, counts):
    """Independent lanes, one real batch inversion, and complete exceptions."""
    assert step is not None and len(states) == len(active)
    p = curve.p
    denominators, ordinary = [], []
    for state, live in zip(states, active):
        use_formula = live and state is not None and state[0] != step[0]
        ordinary.append(use_formula)
        denominators.append((step[0] - state[0]) % p if use_formula else 1)
    prefixes, product = [], 1
    for denominator in denominators:
        assert denominator != 0
        prefixes.append(product)
        product = product * denominator % p
    inverse = pow(product, p - 2, p)
    inverses = [None] * len(states)
    for lane in range(len(states) - 1, -1, -1):
        inverses[lane] = inverse * prefixes[lane] % p
        inverse = inverse * denominators[lane] % p
        assert denominators[lane] * inverses[lane] % p == 1
    counts["batch_inversions"] += 1
    counts["padded_slots"] += len(states)
    result = []
    for state, live, use_formula, inverse in zip(states, active, ordinary, inverses):
        if not live:
            counts["masked_slots"] += 1
            updated = state
        elif use_formula:
            counts["ordinary_slots"] += 1
            slope = (step[1] - state[1]) * inverse % p
            xnew = (slope * slope - state[0] - step[0]) % p
            updated = (xnew, (slope * (state[0] - xnew) - state[1]) % p)
        else:
            if state is None:
                counts["infinity_exceptions"] += 1
            elif state == step:
                counts["doubling_exceptions"] += 1
            else:
                assert state == curve.neg(step)
                counts["opposite_exceptions"] += 1
            updated = curve.add(state, step)
        assert curve.on_curve(updated)
        assert updated == (curve.add(state, step) if live else state)
        result.append(updated)
    return result


def interleaved_stream(curve, offset, step, m, batch_width, counts):
    length = (m + batch_width - 1) // batch_width
    states = [
        curve.add(offset, curve.multiply(lane * length, step))
        if lane * length < m else None
        for lane in range(batch_width)
    ]
    result = []
    for layer in range(length):
        active = [lane * length + layer < m for lane in range(batch_width)]
        for lane, live in enumerate(active):
            if live:
                label = lane * length + layer
                assert states[lane] == curve.add(offset, curve.multiply(label, step))
                result.append(record(states[lane], label))
        # Includes every final unused update and every padded batch slot.
        states = batch_advance(curve, states, step, active, counts)
    assert sorted(rec[2] for rec in result) == list(range(m))
    return result


def signed_match(curve, target, babies, giants, m):
    baby_head = giant_head = 0
    while baby_head < len(babies) and giant_head < len(giants):
        baby, giant = babies[baby_head], giants[giant_head]
        bk, gk = orbit_key(baby, curve.bits), orbit_key(giant, curve.bits)
        if bk < gk:
            baby_head += 1
        elif bk > gk:
            giant_head += 1
        else:
            center = (m - 1) + giant[2] * (2 * m - 1)
            delta = baby[2] if (baby[1] & 1) == (giant[1] & 1) else -baby[2]
            candidate = center + delta
            assert 0 <= candidate < m * (2 * m - 1) < curve.order
            assert curve.multiply(candidate) == target
            return candidate
    return None


def attack(curve, target, m, batch_width, counts, return_streams=False):
    assert 1 <= m and m * (2 * m - 1) < curve.order
    n = 2 * m - 1
    baby_raw = interleaved_stream(curve, None, curve.generator, m, batch_width, counts)
    offset = curve.add(target, curve.neg(curve.multiply(m - 1)))
    step = curve.neg(curve.multiply(n))
    giant_raw = interleaved_stream(curve, offset, step, m, batch_width, counts)
    babies = radix_sort(baby_raw, curve.bits, min(8, curve.bits))
    giants = radix_sort(giant_raw, curve.bits, min(8, curve.bits))
    assert len({orbit_key(rec, curve.bits) for rec in babies}) == m
    answer = signed_match(curve, target, babies, giants, m)
    return (answer, babies, giants) if return_streams else answer


def minimal_grid(order):
    threshold = (39 * order + 99) // 100
    m = max(1, (1 + math.isqrt(1 + 8 * threshold)) // 4)
    while m * (2 * m - 1) < threshold:
        m += 1
    assert m * (2 * m - 1) >= threshold
    assert m == 1 or (m - 1) * (2 * (m - 1) - 1) < threshold
    return threshold, m


def production_ledger():
    threshold, m = minimal_grid(R)
    n = 2 * m - 1
    B, K = 1 << 16, 1 << 43
    L, N = (m + B - 1) // B, 2 * m
    M, A, I, G, F = 4096, 256, 2130176, 78037248, 8448
    assert m == 150264624583687887066834814587324391195
    assert n == 300529249167375774133669629174648782389
    assert L == 2292856210078245347089154275319281
    assert threshold <= m * n < R
    assert 100 * m * n - 39 * R == 781642910497772634128346537537968206357
    assert I == 512 * (M + 64) + 256
    slow_M = 256 * (2 * A + 64) + 256
    slow_I = 512 * (slow_M + 64) + 256
    assert G == slow_I + 16 * slow_M + 32 * A + 4096
    rows = {
        "padded_affine_slots": 2 * B * L * (6 * M + 8 * A + 1024),
        "batch_inversions": 2 * L * I,
        "complete_seeds_and_steps": B * 1025 * G + 1536 * G,
        "record_writes": N * (F + 1024),
        "both_six_pass_sorts": N * 6 * (3 * F + 512),
        "bucket_clearing_and_prefixes": 12 * K * 192,
        "complete_merge": N * (F + 1024),
        "exceptions": 16 * G,
        "reconstruction_verification": 4096 * G,
        "routing_nodes": (1 << 85) * 8192,
        "fixed_input_setup_output": 1 << 44,
    }
    total = sum(rows.values())
    memory = 4 * m * 56 + (1 << 85) * 8192 + K * 64 + (1 << 28)
    assert 56 * (1 << 46) < (1 << 52)
    assert 92 * 84 + 512 == 8240 < F
    assert total == 60634932750577016349807926009615425840505344
    assert total ** 4 < 1000 ** 4 * (1 << 543)
    assert memory == 33659275907062999353028056380885061285792
    assert memory < 1 << 135
    assert 2 * B * L - N == 16842
    # Exact decomposition checks six 43-bit production digits without
    # attempting to allocate a 2^43-bucket test histogram.
    for key in (0, 1 << 256, (1 << 256) | (P - 1), (1 << 257) - 1):
        digits = [(key >> (43 * i)) & (K - 1) for i in range(6)]
        assert sum(digit << (43 * i) for i, digit in enumerate(digits)) == key
    return {
        "m": m, "n": n, "L": L, "success_numerator": m * n,
        "success_denominator": R, "rows_instructions": rows,
        "total_instructions": total, "padded_extra_slots": 2 * B * L - N,
        "score_upper_bound": {"numerator": 543, "denominator": 4},
        "score_integer_inequality_passed": True, "peak_bytes_bound": memory,
        "memory_below_2_to_135": True,
    }


def folding_checks():
    c0, mask = (1 << 32) + 977, (1 << 256) - 1
    edges = [0, 1, 2, c0 - 1, c0, c0 + 1, P - 2, P - 1]
    for bit in (31, 32, 63, 64, 127, 128, 191, 192, 255):
        edges.extend(((1 << bit) - 1, 1 << bit, (1 << bit) + 1))
    pairs = [(a, b) for a in edges for b in edges]
    # Fixed deterministic arithmetic examples, not random sampling evidence.
    value = 0x6A09E667F3BCC909
    for _ in range(256):
        value = (value * 6364136223846793005 + 1442695040888963407) & mask
        left = value % P
        value = (value * 6364136223846793005 + 1442695040888963407) & mask
        pairs.append((left, value % P))
    fifth_limb_cases = 0
    for left, right in pairs:
        z = left * right
        z = (z & mask) + c0 * (z >> 256)
        assert z < 1 << 290
        z = (z & mask) + c0 * (z >> 256)
        assert z < (1 << 256) + (1 << 67)
        z = (z & mask) + c0 * (z >> 256)
        assert z < (1 << 256) + c0 < 2 * P
        fifth_limb_cases += int(z >> 256 != 0)
        reduced = z - P if z >= P else z
        assert 0 <= reduced < P and reduced == left * right % P
    # Explicitly exercise fifth-limb-preserving correction on admissible
    # post-fold intermediate values; no claim that these all occur as products.
    for z in (P - 1, P, P + 1, 1 << 256, (1 << 256) + c0 - 1):
        words = [(z >> (64 * i)) & ((1 << 64) - 1) for i in range(5)]
        restored = sum(word << (64 * i) for i, word in enumerate(words))
        reduced = restored - P if restored >= P else restored
        assert restored == z and reduced == z % P and 0 <= reduced < P
    return {"multiplication_cases": len(pairs), "product_fifth_limb_cases": fifth_limb_cases,
            "explicit_correction_cases": 5}


def mechanism_checks():
    curves = [
        Curve(211, 0, 7, (3, 33), 199),
        Curve(17, 2, 2, (5, 1), 19),
        Curve(7, 0, 3, (1, 2), 13),
        Curve(97, 2, 3, (3, 6), 5),
    ]
    counts, exhaustive_cases, duplicate_cases = Counter(), 0, 0
    per_curve = []
    for curve in curves:
        points = [None]
        for _ in range(1, curve.order):
            points.append(curve.add(points[-1], curve.generator))
        assert len(set(points)) == curve.order
        assert curve.add(points[-1], curve.generator) is None
        curve_cases = 0
        m = 1
        while m * (2 * m - 1) < curve.order:
            for batch_width in (1, 3, 16):
                for secret, target in enumerate(points):
                    before = counts["padded_slots"]
                    answer, babies, giants = attack(
                        curve, target, m, batch_width, counts, return_streams=True)
                    assert answer == (secret if secret < m * (2 * m - 1) else None)
                    length = (m + batch_width - 1) // batch_width
                    assert counts["padded_slots"] - before == 2 * batch_width * length
                    if len({orbit_key(rec, curve.bits) for rec in giants}) != m:
                        duplicate_cases += 1
                    exhaustive_cases += 1
                    curve_cases += 1
            m += 1
        per_curve.append({"field": curve.p, "order": curve.order,
                          "exhaustive_target_grid_batch_cases": curve_cases})
    toy = curves[0]
    answer, babies, giants = attack(toy, toy.multiply(45), 7, 3, counts, True)
    assert answer == 45
    assert len({orbit_key(rec, toy.bits) for rec in giants}) < len(giants)
    # Negative controls: parity is payload, and its sign must be used.
    answer, babies, giants = attack(toy, None, 7, 16, counts, True)
    assert answer == 0
    matches = [(b, g) for b in babies for g in giants
               if orbit_key(b, toy.bits) == orbit_key(g, toy.bits)]
    assert len(matches) == 1
    baby, giant = matches[0]
    assert baby[2] == 6 and giant[2] == 0
    assert (baby[1] & 1) != (giant[1] & 1)
    incorrect_unsigned = 6 + giant[2] * 13 + baby[2]
    assert incorrect_unsigned == 12 and toy.multiply(incorrect_unsigned) is not None
    assert not {(b[0], b[1]) for b in babies}.intersection(
        {(g[0], g[1]) for g in giants})
    # Equal orbit keys with distinct parity/index payloads must remain stable.
    records = [(4, 3, 9), (4, 2, 7), (0, 0, 11), (4, 3, 2), (0, 2, 5)]
    for digit_bits in (1, 2, 3, 7):
        assert radix_sort(records, 8, digit_bits) == [
            records[2], records[4], records[0], records[1], records[3]]
    # All exceptional branches, including a separate order-two point.
    direct = batch_advance(toy, [None, toy.generator, toy.neg(toy.generator), None],
                           toy.generator, [True, True, True, False], counts)
    assert direct == [toy.generator, toy.multiply(2), None, None]
    order_two = Curve(17, 0, 7, (3, 0), 2)
    assert order_two.add(order_two.generator, order_two.generator) is None
    assert batch_advance(order_two, [order_two.generator], order_two.generator,
                         [True], counts) == [None]
    assert all(counts[key] > 0 for key in (
        "masked_slots", "ordinary_slots", "infinity_exceptions",
        "doubling_exceptions", "opposite_exceptions"))
    assert duplicate_cases > 0
    return {"curves": per_curve, "exhaustive_cases": exhaustive_cases,
            "cases_with_duplicate_giant_orbits": duplicate_cases,
            "batch_counters": dict(counts), "sign_negative_controls_passed": True,
            "stable_payload_order_passed": True, "zero_y_exception_passed": True}


def production_field_small_grids():
    curve, counts, cases = Curve(P, 0, 7, GENERATOR, R), Counter(), 0
    assert curve.multiply(R) is None
    for m in (7, 11):
        n, coverage = 2 * m - 1, m * (2 * m - 1)
        secrets = sorted({0, 1, m - 2, m - 1, m, n - 1, n, 2 * n - 1,
                          coverage - 1, coverage, coverage + 1, R - 1})
        for batch_width in (4, 16):
            for secret in secrets:
                answer = attack(curve, curve.multiply(secret), m, batch_width, counts)
                assert answer == (secret if secret < coverage else None)
                cases += 1
    return {"planted_small_grid_cases": cases, "batch_counters": dict(counts)}


def main():
    raw = sys.stdin.buffer.read(4097)
    if len(raw) > 4096 or json.loads(raw) != {}:
        raise ValueError("The only permitted bounded experiment input is {}")
    result = {"status": "passed", "production_ledger": production_ledger(),
              "field_folding": folding_checks(), "small_curve_mechanisms": mechanism_checks(),
              "secp256k1_small_grids": production_field_small_grids(),
              "limits": [
                  "No production-size search was executed.",
                  "No compiled AArch64 instruction counts were measured.",
                  "No physical distributed routing realization was tested.",
                  "Toy coverage and planted full-field cases do not estimate production success.",
                  "The success argument and instruction ledger are analytic claims.",
              ]}
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
