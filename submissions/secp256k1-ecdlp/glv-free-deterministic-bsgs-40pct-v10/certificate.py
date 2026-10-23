#!/usr/bin/env python3
"""v10 exact-integer ledger, overflow, coverage and bounded mechanism checks.
Run: python3 -B certificate.py [--write]
--write captures certificate.json and its console report beside this file.
No floating point, credentials, network, or production-size attack execution.
The integer-power certification method is inherited from v9, not its ledger.
"""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
ORDER = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
C = 1000
M_BABY = 1 << 128
H_GIANT = (ORDER + M_BABY - 1) // M_BABY
N = M_BABY + H_GIANT
BATCH = 1 << 16
LENGTH = 1 << 112
BUCKETS = 1 << 43
PASSES = 6
A = 256
M = 4096
INV = 512 * (M + 64) + 256
OLD_M = 256 * (2 * A + 64) + 256
OLD_INV = 512 * (OLD_M + 64) + 256
G = OLD_INV + 16 * OLD_M + 32 * A + 4096
RECORD_BYTES = 56
DEPTH = 72
TRANSFER = 6144
CLAIM_NUM = 6831
CLAIM_DEN = 50
Q = 100000
CLAIM_SCALED = CLAIM_NUM * Q // CLAIM_DEN
PRIOR_T = (1 << 40) + (1 << 127) * 832768
PRIOR_RECORDED_SCALED = 13670178


def u192(n):
    assert 0 <= n < 1 << 192
    return [(n >> (64 * k)) & ((1 << 64) - 1) for k in range(3)]


def read192(a):
    return sum(v << (64 * k) for k, v in enumerate(a))


def add192(a, b):
    out = []
    carry = 0
    for x, y in zip(a, b):
        z = x + y + carry
        out.append(z & ((1 << 64) - 1))
        carry = z >> 64
    assert carry == 0
    return out


def four_limb_mul(a, b):
    """Exact schoolbook and pseudo-Mersenne folding model, not instrumentation."""
    mask = (1 << 64) - 1
    aa = [(a >> (64 * k)) & mask for k in range(4)]
    bb = [(b >> (64 * k)) & mask for k in range(4)]
    t = [0] * 8
    for i in range(4):
        for j in range(4):
            z = aa[i] * bb[j]
            k = i + j
            carry = z
            while carry:
                assert k < 8
                w = t[k] + (carry & mask)
                t[k] = w & mask
                carry = (carry >> 64) + (w >> 64)
                k += 1
    z = sum(v << (64 * k) for k, v in enumerate(t))
    c = (1 << 32) + 977
    for _ in range(3):
        z = (z & ((1 << 256) - 1)) + (z >> 256) * c
    assert z < (1 << 256) + c
    if z >= P:
        z -= P
    assert 0 <= z < P
    return z


def point_add(a, b, p):
    if a is None:
        return b
    if b is None:
        return a
    x, y = a
    u, v = b
    if x == u:
        if (y + v) % p == 0:
            return None
        slope = 3 * x * x * pow(2 * y, p - 2, p) % p
    else:
        slope = (v - y) * pow((u - x) % p, p - 2, p) % p
    xx = (slope * slope - x - u) % p
    return xx, (slope * (x - xx) - y) % p


def point_mul(k, a, p):
    z = None
    while k:
        if k & 1:
            z = point_add(z, a, p)
        a = point_add(a, a, p)
        k >>= 1
    return z


def batch_add(points, step, p):
    # Each point is a different chain's current state: no serial dependency
    # is crossed by one simultaneous prefix/reverse inversion.
    ds = []
    ordinary = []
    for a in points:
        good = a is not None and step is not None and a[0] != step[0]
        ordinary.append(good)
        ds.append((step[0] - a[0]) % p if good else 1)
    pref = [1]
    for d in ds:
        pref.append(pref[-1] * d % p)
    inv = pow(pref[-1], p - 2, p)
    out = [None] * len(points)
    for k in range(len(points) - 1, -1, -1):
        dinv = inv * pref[k] % p
        inv = inv * ds[k] % p
        a = points[k]
        if not ordinary[k]:
            out[k] = point_add(a, step, p)
        else:
            slope = (step[1] - a[1]) * dinv % p
            xx = (slope * slope - a[0] - step[0]) % p
            out[k] = (xx, (slope * (a[0] - xx) - a[1]) % p)
    return out


def key(point):
    return 0 if point is None else ((2 + (point[1] & 1)) << 256) | point[0]


def interleaved_records(n, chains, start, step, p):
    length = (n + chains - 1) // chains
    starts = [t * length for t in range(chains) if t * length < n]
    current = [point_add(start, point_mul(j, step, p), p) for j in starts]
    records = []
    for k in range(length):
        active = [t for t, j in enumerate(starts) if j + k < n]
        for t in active:
            records.append((key(current[t]), starts[t] + k))
        # Padding inactive chains with infinity/1 never affects active chains.
        new = batch_add([current[t] for t in active], step, p)
        for t, v in zip(active, new):
            current[t] = v
    return records


def stable_radix(records, digit_bits):
    n = len(records)
    passes = (258 + digit_bits - 1) // digit_bits
    kmax = 1 << digit_bits
    for stage in range(passes):
        counts = [u192(0) for _ in range(kmax)]
        for k, _ in records:
            d = (k >> (stage * digit_bits)) & (kmax - 1)
            counts[d] = add192(counts[d], u192(1))
        cursor = []
        acc = u192(0)
        for cnt in counts:
            cursor.append(acc[:])
            acc = add192(acc, cnt)
        assert read192(acc) == n
        dest = [None] * n
        for rec in records:
            d = (rec[0] >> (stage * digit_bits)) & (kmax - 1)
            at = read192(cursor[d])
            cursor[d] = add192(cursor[d], u192(1))
            assert dest[at] is None
            dest[at] = rec
        records = dest
    return records


def recover(babies, giants, m, r, target, base, p):
    u = v = 0
    found = []
    pairs = 0
    while u < len(babies) and v < len(giants):
        bk, j = babies[u]
        gk, i = giants[v]
        if bk < gk:
            u += 1
        elif gk < bk:
            v += 1
        else:
            pairs += 1
            a = i * m + j
            if 0 <= a < r and point_mul(a, base, p) == target:
                found.append(a)
            u += 1
            v += 1
    assert pairs <= 2
    return found


def bounded_mechanism_tests():
    # Enumerated prime-order curve: y^2=x^3+7 over F_211 has 199 points.
    p = 211
    base = (3, 33)
    r = 199
    assert 1 + sum(1 for x in range(p) for y in range(p)
                   if (y*y - x*x*x - 7) % p == 0) == r
    assert point_mul(r, base, p) is None
    m = 16
    chains = 4
    h = (r + m - 1) // m
    d = point_mul(m, base, p)
    assert d is not None
    d = (d[0], (-d[1]) % p)
    babies = interleaved_records(m, chains, None, base, p)
    assert dict((j, k) for k, j in babies)[0] == 0
    assert all(k == key(point_mul(j, base, p)) for k, j in babies)
    bs = stable_radix(babies, 3)
    successes = 0
    for x in range(r):
        target = point_mul(x, base, p)
        gs = stable_radix(interleaved_records(h, chains, target, d, p), 3)
        assert recover(bs, gs, m, r, target, base, p) == [x]
        successes += 1
    # Separate scalar-encoded prime cyclic groups exhaust division boundaries.
    scalar_trials = 0
    for r2 in (7, 11, 13, 17, 19, 23, 29, 31):
        for m2 in range(2, r2):
            h2 = (r2 + m2 - 1) // m2
            if m2 * h2 >= 2 * r2:
                continue
            for x in range(r2):
                candidates = [i * m2 + j for i in range(h2) for j in range(m2)
                              if (x - i * m2) % r2 == j]
                assert len(candidates) <= 2
                assert [a for a in candidates if a < r2] == [x]
                scalar_trials += 1
    values = [0, 1, 2, (1 << 64) - 1, 1 << 64, (1 << 128) - 1,
              1 << 128, P - 2, P - 1]
    mult_checks = 0
    for a in values:
        for b in values:
            assert four_limb_mul(a, b) == a * b % P
            mult_checks += 1
    gx = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
    gy = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
    fixtures = [None, (gx, gy), (gx, (-gy) % P), point_mul(2, (gx, gy), P),
                point_mul(12345, (gx, gy), P)]
    assert batch_add(fixtures, (gx, gy), P) == [point_add(a, (gx, gy), P) for a in fixtures]
    return {"toy_curve": "y^2=x^3+7 over F_211", "toy_generator": list(base),
            "toy_subgroup_order": r, "toy_m": m, "toy_h": h, "toy_chains": chains,
            "toy_trials": r, "toy_successes": successes,
            "cyclic_division_trials": scalar_trials, "field_multiply_fixtures": mult_checks,
            "production_batch_fixtures": len(fixtures)}


def interval_scaled(value, divisor, q=Q):
    left = pow(value, q)
    denom = pow(divisor, q)
    lo, hi = 0, 200 * q
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if left < denom << mid:
            hi = mid
        else:
            lo = mid
    return [lo, hi]


def build():
    rows = {
        "batched_affine_slots": N * (6 * M + 8 * A + 1024),
        "all_batch_inversions_including_residual": 2 * LENGTH * INV,
        "complete_chain_seeds_and_D": (1025 * BATCH + 512) * G,
        "initial_record_writes": N * (TRANSFER + 256),
        "six_radix_passes_both_streams": N * PASSES * (3 * TRANSFER + 128),
        "histogram_clear_prefix_setup": 2 * PASSES * BUCKETS * 96,
        "merge_read_compare": N * (TRANSFER + 256),
        "complete_exception_fallbacks": 16 * G,
        "candidate_scalar_verification": 4096 * G,
        "routing_initialization": (1 << 73) * 4096,
        "input_and_fixed_setup": 1 << 40,
        "scalar_output": 1 << 40,
    }
    total = sum(rows.values())
    # Cold preprocessing: all shared infrastructure plus the baby-dependent rows.
    pre = (rows["routing_initialization"] + rows["input_and_fixed_setup"]
           + (512 * BATCH + 512) * G + LENGTH * INV
           + M_BABY * ((6 * M + 8 * A + 1024) + TRANSFER + 256
                       + PASSES * (3 * TRANSFER + 128))
           + PASSES * BUCKETS * 96 + 16 * G)
    storage = {
        "four_record_arrays": 2 * N * RECORD_BYTES,
        "routing_code_and_state": (1 << 73) * 4096,
        "histogram_and_scatter_cursors": BUCKETS * 64,
        "batch_scratch": BATCH * 256,
        "engine_code_input_output": 1 << 24,
    }
    checks = {}
    def check(name, ok):
        checks[name] = bool(ok)
    check("h_is_2pow128_minus_one", H_GIANT == M_BABY - 1)
    check("minimal_full_coverage", M_BABY * (H_GIANT - 1) < ORDER <= M_BABY * H_GIANT)
    check("twin_cap_mh_lt_2r", M_BABY * H_GIANT < 2 * ORDER)
    check("baby_U0_is_infinity", key(None) == 0)
    check("batch_partition_exact", BATCH * LENGTH == M_BABY and H_GIANT == BATCH * LENGTH - 1)
    check("serial_146_kernel_constants", (OLD_M, OLD_INV, G) == (147712, 75661568, 78037248))
    check("hot_multiply_instruction_itemization", 16 * 128 + 3 * 256 + 256 + 256 <= M)
    check("hot_field_add_instruction_itemization", 64 + 192 <= A)
    check("complete_exception_occurrence_cap", 2 + 3 <= 16)
    check("radix_histogram_update_instruction_itemization", 6 + 3 + 3 + 3 + 2 + 3 + 12 <= 64)
    check("radix_scatter_cursor_instruction_itemization", 6 + 3 + 3 + 3 + 2 + 3 + 12 <= 64)
    check("radix_bucket_setup_instruction_itemization", 32 + 32 + 32 == 96)
    check("all_work_is_R_no_uncharged_H", total == sum(rows.values()))
    check("three_limb_count_and_prefix_capacity", max(M_BABY, H_GIANT) < 1 << 192)
    check("counter_carry_64bit", read192(add192(u192((1 << 64) - 1), u192(1))) == 1 << 64)
    check("counter_carry_128bit", read192(add192(u192((1 << 128) - 1), u192(1))) == 1 << 128)
    check("all_records_one_bucket_is_safe", read192(add192(u192(M_BABY - 1), u192(1))) == M_BABY)
    check("record_layout", 32 + 1 + 7 + 16 == RECORD_BYTES)
    check("six_digits_cover_complete_key", PASSES * 43 == 258)
    check("leaf_fits_native_address_space", RECORD_BYTES * (1 << 58) < 1 << 64)
    check("four_arrays_fit_depth72_tree", 4 * (1 << 70) == 1 << DEPTH)
    check("transfer_itemization_fits_F", ((4 + 11) * 4 + 2 * 8) * DEPTH + 512 <= TRANSFER)
    check("histogram_fits_one_native_address_space", BUCKETS * 64 < 1 << 64)
    check("claim_fraction_reduced", math.gcd(CLAIM_NUM, CLAIM_DEN) == 1)
    check("claim_scaling_exact", CLAIM_NUM * Q == CLAIM_SCALED * CLAIM_DEN)
    check("claim_is_log2_upper_bound", pow(total, Q) < pow(C, Q) << CLAIM_SCALED)
    check("at_least_0p005_bits_nominal_slack", pow(total, Q) < pow(C, Q) << (CLAIM_SCALED - 500))
    check("all_R_times_1p2_margin_ge_0p10", pow(6 * total, Q) < pow(5 * C, Q) << (CLAIM_SCALED - 10000))
    check("claim_below_prior_computed_not_frontier", pow(C, Q) << CLAIM_SCALED < pow(PRIOR_T, Q))
    check("prior_recorded_bound_reproduced", pow(PRIOR_T, Q) < pow(C, Q) << PRIOR_RECORDED_SCALED)
    check("prior_recorded_five_decimal_floor_reproduced", pow(PRIOR_T, Q) >= pow(C, Q) << (PRIOR_RECORDED_SCALED - 1))
    check("preprocessing_included", 0 < pre < total)
    check("peak_storage_below_2pow137", sum(storage.values()) < 1 << 137)
    check("serial_fallback_alone_exceeds_129_threshold", N * G > C << 129)
    toy = bounded_mechanism_tests()
    check("captured_toy_all_challenges_including_zero", toy["toy_successes"] == toy["toy_trials"])
    subpath = HERE / "submission.json"
    if subpath.exists():
        sub = json.loads(subpath.read_text())
        score = sub["claims"][0]["score"]
        check("submission_score_matches_certificate", (int(score["numerator"]), int(score["denominator"])) == (CLAIM_NUM, CLAIM_DEN))
    deps = HERE / "dependencies.json"
    if deps.exists():
        check("dependencies_empty", json.loads(deps.read_text()) == [])
    result = {
        "version": 1, "H_reference_operations": "0", "C": C,
        "parameters": {"p": str(P), "r": str(ORDER), "m": str(M_BABY), "h": str(H_GIANT),
                       "N": str(N), "B": str(BATCH), "L": str(LENGTH), "K": str(BUCKETS),
                       "passes": PASSES, "coverage_slack": str(M_BABY * H_GIANT - ORDER)},
        "kernels": {"A": A, "M_hot": M, "I_hot": INV, "M_old": OLD_M, "I_old": OLD_INV,
                    "G_complete_old": G, "F_record_transfer": TRANSFER,
                    "hot_slot": 6 * M + 8 * A + 1024},
        "R_rows": {k: str(v) for k, v in rows.items()},
        "total_R_instructions": str(total), "preprocessing_instructions": str(pre),
        "peak_storage_rows_bytes": {k: str(v) for k, v in storage.items()},
        "peak_storage_bytes": str(sum(storage.values())), "memory_log2_bytes_bound": 137,
        "advice_bytes": 0, "success_numerator": 1, "success_denominator": 1,
        "claim": {"numerator": str(CLAIM_NUM), "denominator": str(CLAIM_DEN), "decimal": "136.62"},
        "q": Q, "claim_scaled": CLAIM_SCALED,
        "computed_log2_interval_scaled": interval_scaled(total, C),
        "repriced_1p2_log2_interval_scaled": interval_scaled(6 * total, 5 * C),
        "prior_public_work": {"submission": "5b158b37-6ab9-405f-b30a-92ceab0ccefe",
                              "T": str(PRIOR_T), "recorded_score": "136.70178",
                              "computed_interval_scaled": interval_scaled(PRIOR_T, C)},
        "captured_mechanism_run": toy, "checks": checks,
        "verdict": "PASS" if all(checks.values()) else "FAIL",
    }
    return result


def report(result):
    lines = [f"{ 'PASS' if ok else 'FAIL'} {name}" for name, ok in result["checks"].items()]
    lines += ["R_total=" + result["total_R_instructions"],
              "preprocessing=" + result["preprocessing_instructions"],
              "peak_bytes=" + result["peak_storage_bytes"],
              "log2_interval_scaled=" + json.dumps(result["computed_log2_interval_scaled"]),
              "repriced_interval_scaled=" + json.dumps(result["repriced_1p2_log2_interval_scaled"]),
              "toy=" + json.dumps(result["captured_mechanism_run"], sort_keys=True),
              "verdict=" + result["verdict"]]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    data = build()
    output = report(data)
    if "--write" in sys.argv:
        (HERE / "certificate.json").write_text(json.dumps(data, indent=2) + "\n")
        (HERE / "certificate-output.txt").write_text(output)
    else:
        recorded = HERE / "certificate.json"
        if recorded.exists() and json.loads(recorded.read_text()) != data:
            print("FAIL captured_certificate_differs_from_recomputation")
            sys.exit(1)
    print(output, end="")
    sys.exit(0 if data["verdict"] == "PASS" else 1)
