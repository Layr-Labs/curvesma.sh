"""Bounded supporting evidence for the class Pollard-rho DP candidate (v3).

Standard library only. Reads JSON on stdin. It checks the production secp256k1
endomorphism constants, runs the complete class walk with distinguished points
and coefficient re-run extraction on a small additive group, and recomputes the
analytic production ledger. It does not run a production attack or measure
AArch64 instructions.
"""

import json
import math
import random
import sys

P_FIELD = int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16)
R_ORDER = int("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141", 16)
GEN = (
    int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
    int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16),
)
BETA = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
LAM = int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16)
MASK64 = (1 << 64) - 1
C_INSTRUCTIONS = 1000
E_INSTRUCTIONS = 630
DP_BITS = 56
WALKERS = 1 << 10
NU = 1.03


def inv_mod(a, m):
    return pow(a % m, m - 2, m)


def ec_add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if (y1 + y2) % P_FIELD == 0:
            return None
        slope = (3 * x1 * x1) * inv_mod(2 * y1, P_FIELD) % P_FIELD
    else:
        slope = (y2 - y1) * inv_mod(x2 - x1, P_FIELD) % P_FIELD
    x3 = (slope * slope - x1 - x2) % P_FIELD
    y3 = (slope * (x1 - x3) - y1) % P_FIELD
    return (x3, y3)


def ec_mul(k, point):
    result = None
    addend = point
    while k:
        if k & 1:
            result = ec_add(result, addend)
        addend = ec_add(addend, addend)
        k >>= 1
    return result


def profile_check():
    phi_p = ((BETA * GEN[0]) % P_FIELD, GEN[1])
    orbit_sizes = []
    for k in (1, 2, 3, 5, 123456789):
        point = ec_mul(k, GEN)
        orbit = set()
        for e in range(3):
            x_e = pow(BETA, e, P_FIELD) * point[0] % P_FIELD
            orbit.add((x_e, point[1]))
            orbit.add((x_e, (-point[1]) % P_FIELD))
        orbit_sizes.append(len(orbit))
    return {
        "beta_cubed_is_one": pow(BETA, 3, P_FIELD) == 1,
        "beta_is_nontrivial": BETA % P_FIELD != 1,
        "lambda_cubed_is_one": pow(LAM, 3, R_ORDER) == 1,
        "lambda_is_nontrivial": LAM % R_ORDER != 1,
        "phi_P_equals_lambda_P": phi_p == ec_mul(LAM, GEN),
        "x_zero_is_off_curve": pow(7, (P_FIELD - 1) // 2, P_FIELD) == P_FIELD - 1,
        "generic_orbit_sizes": orbit_sizes,
        "all_generic_orbits_are_six": all(size == 6 for size in orbit_sizes),
    }


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def find_lambda(n):
    g = 2
    while g < n:
        lam = pow(g, (n - 1) // 3, n)
        if lam != 1 and pow(lam, 3, n) == 1:
            return lam
        g += 1
    raise ValueError("no nontrivial cube root of unity")


def automorphisms(n, lam):
    return sorted({
        pow(lam, e, n) if sign == 0 else (-pow(lam, e, n)) % n
        for e in range(3)
        for sign in (0, 1)
    })


def canonical(value, automorphism_list, n):
    values = [(value * h) % n for h in automorphism_list]
    best = min(values)
    return best, automorphism_list[values.index(best)]


def step_index(canonical_x):
    return ((canonical_x * 6364136223846793005 + 1442695040888963407) & MASK64)


def walk_step(point, a, b, table, alphas, betas, autos, n, increments):
    canonical_point, _ = canonical(point, autos, n)
    k = step_index(canonical_point) % increments
    raw = (point + table[k]) % n
    a = (a + alphas[k]) % n
    b = (b + betas[k]) % n
    canonical_point, unit = canonical(raw, autos, n)
    return canonical_point, (unit * a) % n, (unit * b) % n


def run_seed(seed, table, alphas, betas, autos, n, increments, dp_bits, cap, track=False):
    rng = random.Random(seed)
    a = rng.randrange(n)
    b = 0
    point, unit = canonical(a, autos, n)
    a = (unit * a) % n
    b = (unit * b) % n
    for length in range(cap + 1):
        canonical_point, _ = canonical(point, autos, n)
        if canonical_point % (1 << dp_bits) == 0:
            if track:
                return canonical_point, seed, length, a, b
            return canonical_point, seed, length
        if length == cap:
            break
        point, a, b = walk_step(point, a, b, table, alphas, betas, autos, n, increments)
    return None


def toy_trial(order, trials, seed, increments, dp_bits, walkers):
    if order < 7 or not is_prime(order) or (order - 1) % 6 != 0:
        raise ValueError("order must be a prime congruent to 1 mod 6")
    lam = find_lambda(order)
    autos = automorphisms(order, lam)
    if len(autos) != 6:
        raise ValueError("automorphism group is not of order 6")
    rng = random.Random(seed)
    classes = (order - 1) // 6 + 1
    cap = 40 * (1 << dp_bits)
    recovered = 0
    failed = 0
    collisions = 0
    rechecks = 0
    for _ in range(trials):
        x = rng.randrange(1, order)
        alphas = [rng.randrange(order) for _ in range(increments)]
        betas = [rng.randrange(order) for _ in range(increments)]
        table = [(alphas[k] + betas[k] * x) % order for k in range(increments)]
        seen = {}
        solved = False
        for _ in range(walkers):
            walk_seed = rng.randrange(1 << 31)
            record = run_seed(walk_seed, table, alphas, betas, autos, order,
                              increments, dp_bits, cap)
            if record is None:
                continue
            canonical_point, old_seed, old_length = record
            if canonical_point in seen:
                collisions += 1
                first_seed, first_length = seen[canonical_point]
                run1 = run_seed(walk_seed, table, alphas, betas, autos, order,
                                increments, dp_bits, cap, track=True)
                run2 = run_seed(first_seed, table, alphas, betas, autos, order,
                                increments, dp_bits, cap, track=True)
                if run1 and run2:
                    c1, _, l1, a1, b1 = run1
                    c2, _, l2, a2, b2 = run2
                    if c1 == c2 and l1 == old_length and l2 == first_length:
                        rechecks += 1
                        denominator = (b1 - b2) % order
                        if denominator != 0:
                            candidate = ((a2 - a1) * inv_mod(denominator, order)) % order
                            if candidate == x % order:
                                solved = True
                                recovered += 1
                                break
            else:
                seen[canonical_point] = (walk_seed, old_length)
        if not solved:
            failed += 1
    return {
        "order": order,
        "lambda": str(lam),
        "automorphisms": [str(h) for h in autos],
        "classes": classes,
        "trials": trials,
        "walkers": walkers,
        "increments": increments,
        "dp_bits": dp_bits,
        "recovered": recovered,
        "failed": failed,
        "collisions": collisions,
        "rechecks": rechecks,
    }


def production_ledger():
    classes = (R_ORDER - 1) // 6
    k0 = math.sqrt(0.98859 * classes)
    theta = 1 << DP_BITS
    lag = WALKERS * 40 * theta
    walk_steps = math.ceil(NU * k0) + lag
    total_instructions = walk_steps * (C_INSTRUCTIONS + E_INSTRUCTIONS) + (1 << 83) * C_INSTRUCTIONS
    time_log2 = math.log2(total_instructions) - math.log2(C_INSTRUCTIONS)
    ideal_success = 1 - math.exp(-(NU * k0) ** 2 / (2 * classes))
    required_quality = -math.log(0.61) / ((NU * NU) * 0.98859 / 2)
    return {
        "r": str(R_ORDER),
        "p": str(P_FIELD),
        "classes_N_prime": str(classes),
        "k0": k0,
        "k0_log2": math.log2(k0),
        "nu": NU,
        "walk_steps": str(walk_steps),
        "walk_steps_log2": math.log2(walk_steps),
        "step_instructions": C_INSTRUCTIONS + E_INSTRUCTIONS,
        "step_reference_equivalents": (C_INSTRUCTIONS + E_INSTRUCTIONS) / C_INSTRUCTIONS,
        "other_ref_ops": str(1 << 83),
        "total_instructions": str(total_instructions),
        "time_log2": time_log2,
        "submitted_time_log2": math.ceil(time_log2 * 1_000_000) / 1_000_000,
        "ideal_success_rho_1": ideal_success,
        "required_walk_quality_for_0_39": required_quality,
        "memory_log2_bytes": 76.6,
        "preprocessing_log2": 32.01,
        "nonuniform_advice_bytes": "0",
    }


def main():
    try:
        config = json.load(sys.stdin)
        order = int(config["order"])
        trials = int(config["trials"])
        seed = int(config["seed"])
        increments = int(config["increments"])
        dp_bits = int(config["dp_bits"])
        walkers = int(config["walkers"])
        if not (7 <= order <= 1_000_003):
            raise ValueError("order out of range")
        if not (1 <= trials <= 100):
            raise ValueError("trials out of range")
        if not (0 <= seed <= (1 << 31) - 1):
            raise ValueError("seed out of range")
        if not (2 <= increments <= 256):
            raise ValueError("increments out of range")
        if not (1 <= dp_bits <= 12):
            raise ValueError("dp_bits out of range")
        if not (1 <= walkers <= 512):
            raise ValueError("walkers out of range")
    except Exception as exc:
        print(json.dumps({"error": "invalid input", "detail": str(exc)}, sort_keys=True))
        return

    toy = toy_trial(order, trials, seed, increments, dp_bits, walkers)
    profile = profile_check()
    ledger = production_ledger()
    checks = {
        "profile_constants_ok": all([
            profile["beta_cubed_is_one"],
            profile["beta_is_nontrivial"],
            profile["lambda_cubed_is_one"],
            profile["lambda_is_nontrivial"],
            profile["phi_P_equals_lambda_P"],
            profile["x_zero_is_off_curve"],
            profile["all_generic_orbits_are_six"],
        ]),
        "toy_recovered_all": toy["recovered"] == trials,
        "toy_rechecks_match": toy["rechecks"] >= 1,
        "ledger_total_below_2^128_ref_ops": int(ledger["total_instructions"]) < (1 << 128) * C_INSTRUCTIONS,
        "ledger_success_above_0.39": ledger["ideal_success_rho_1"] > 0.39,
        "ledger_required_quality_at_most_0.943": ledger["required_walk_quality_for_0_39"] <= 0.943,
        "submitted_score_below_127.47": ledger["submitted_time_log2"] < 127.47,
    }
    output = {
        "scope": (
            "Bounded supporting evidence only. The toy group runs the complete DP walk, "
            "merging, re-run extraction and scalar recovery at small order; the profile "
            "check verifies the secp256k1 endomorphism constants; the ledger recomputes "
            "the analytic integers. No production attack is run and no AArch64 "
            "instruction count is measured."
        ),
        "profile_check": profile,
        "toy": toy,
        "ledger": ledger,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
