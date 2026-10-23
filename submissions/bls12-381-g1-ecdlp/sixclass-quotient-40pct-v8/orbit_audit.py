"""Exhaustive tiny-curve subgroup and 6x6 record-digit audit."""
from hashlib import sha256
import json


def audit_orbits(curve, inputs, shifts, table):
    assert len(set(inputs)) == curve.r
    assert curve.add(inputs[-1], curve.P) is None
    assert curve.lam != 1 and pow(curve.lam, 3, curve.r) == 1
    assert curve.beta != 1 and pow(curve.beta, 3, curve.p) == 1
    zeta = pow(11, 64, curve.r)
    Q = ((curve.r - 1) // 64) // 6
    multipliers = [pow(zeta, shift * Q, curve.r) for shift in shifts]
    classes, fixed = {}, [0] * 5
    exceptions = alignment_errors = 0
    for scalar, point in enumerate(inputs[1:], 1):
        orbit = curve.orbit(point)
        exceptions += len(set(orbit)) != 6
        alignment_errors += curve.tau(point) != inputs[curve.lam * scalar % curve.r]
        for j in range(1, 6):
            fixed[j - 1] += orbit[j] == point
        assert all((y * y - x * x * x - 4) % curve.p == 0 for x, y in orbit)
        assert all(inputs[m * scalar % curve.r] == orbit[j] for j, m in enumerate(multipliers))
        rec = curve.canonical_record(point, 0, 0)
        digit = int.from_bytes(rec[120:], 'big') >> 1
        x, y = int.from_bytes(rec[:48], 'big'), int.from_bytes(rec[48:96], 'big')
        assert orbit[digit] == (x, y)
        members = classes.setdefault(rec[:104], {})
        assert digit not in members
        members[digit] = scalar
    wrong = count = 0
    pair_histogram = [0] * 36
    transcript = sha256()
    for key, members in sorted(classes.items()):
        assert set(members) == set(range(6))
        for u in range(6):
            for v in range(6):
                correction = pow(zeta, table[u][v] * Q, curve.r)
                wrong += members[v] != correction * members[u] % curve.r
                count += 1
                pair_histogram[6 * u + v] += 1
                transcript.update((json.dumps([members[u], members[v], u, v, correction], separators=(',', ':')) + '\n').encode())
    assert exceptions == alignment_errors == wrong == sum(fixed) == 0
    assert count == (curve.r - 1) // 6 * 36
    return dict(nonzero_points=curve.r - 1, orbit_classes=len(classes),
                orbit_size_exceptions=exceptions, fixed_nonidentity_counts=fixed,
                scalar_field_alignment_errors=alignment_errors, digit_pair_tests=count,
                wrong_digit_relations=wrong, digit_pair_histogram=pair_histogram,
                digit_shifts=shifts, digit_table=table, transcript_sha256=transcript.hexdigest())
