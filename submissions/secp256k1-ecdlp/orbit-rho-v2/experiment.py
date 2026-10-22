"""Bounded supporting evidence for the orbit-collision r-adding walk.

This program is deliberately small and uses only the Python standard library. It
does three things:

1. Checks the production secp256k1 endomorphism constants and orbit sizes.
2. Runs the same orbit-class walk on a small additive cyclic group, where the
   full walk can be simulated exhaustively, and compares the observed first-repeat
   time with sqrt(pi*N/2).
3. Recomputes the exact production ledger integers and checks the inequalities.

It does not run the production attack, measure AArch64 instructions, or certify a
production cost bound. Those are analytic claims in proof.md. The toy check is
evidence only at its small order.
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
STEP_INSTRUCTIONS = 1814
INCREMENTS = 1 << 20
SPEED_CONSTANT = 1.05


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


def orbit_walk(order, trials, seed, increments):
    if order < 7 or not is_prime(order) or (order - 1) % 6 != 0:
        raise ValueError("order must be a prime congruent to 1 mod 6")
    lam = find_lambda(order)
    automorphisms = sorted({
        value
        for e in range(3)
        for value in (pow(lam, e, order), (-pow(lam, e, order)) % order)
    })
    if len(automorphisms) != 6:
        raise ValueError("automorphism group is not of order 6")
    classes = (order - 1) // 6 + 1
    rng = random.Random(seed)
    ratios = []
    recovered = 0
    failed = 0
    step_counts = []
    for _ in range(trials):
        x = rng.randrange(1, order)
        alpha = [rng.randrange(order) for _ in range(increments)]
        beta = [rng.randrange(order) for _ in range(increments)]
        increments_points = [(alpha[k] + beta[k] * x) % order
                             for k in range(increments)]
        a = rng.randrange(order)
        b = 0
        point = a
        seen = {}
        steps = 0
        limit = max(100, 10 * math.isqrt(classes) + 100)
        success = False
        for _ in range(limit):
            steps += 1
            canonical = min((point * h) % order for h in automorphisms)
            if canonical in seen:
                old_point, old_a, old_b = seen[canonical]
                w = None
                for h in automorphisms:
                    if (h * old_point) % order == point:
                        w = h
                        break
                if w is None:
                    failed += 1
                    break
                denominator = (b - w * old_b) % order
                if denominator == 0:
                    failed += 1
                    break
                candidate = ((w * old_a - a) * inv_mod(denominator, order)) % order
                if candidate == x % order:
                    success = True
                    recovered += 1
                else:
                    failed += 1
                break
            seen[canonical] = (point, a, b)
            k = ((canonical * 6364136223846793005 + 1442695040888963407)
                 & MASK64) % increments
            point = (point + increments_points[k]) % order
            a = (a + alpha[k]) % order
            b = (b + beta[k]) % order
        else:
            failed += 1
        ratios.append(steps / math.sqrt(classes))
        step_counts.append(steps)
    return {
        "order": order,
        "lambda": str(lam),
        "automorphisms": [str(h) for h in automorphisms],
        "classes": classes,
        "trials": trials,
        "increments": increments,
        "recovered": recovered,
        "failed": failed,
        "mean_steps": sum(step_counts) / len(step_counts),
        "mean_ratio_to_sqrt_classes": sum(ratios) / len(ratios),
        "ideal_ratio": math.sqrt(math.pi / 2),
        "min_ratio": min(ratios),
        "max_ratio": max(ratios),
        "budget_ratio": SPEED_CONSTANT,
        "success_fraction_at_budget": sum(
            1 for count in step_counts if count <= SPEED_CONSTANT * math.sqrt(classes)
        ) / trials,
        "ideal_success_at_budget": 1 - math.exp(-(SPEED_CONSTANT ** 2) / 2),
    }


def production_ledger():
    classes = (R_ORDER - 1) // 6 + 1
    steps = math.ceil(SPEED_CONSTANT * math.sqrt(classes))
    preprocessing = INCREMENTS * 1024 + INCREMENTS * 512
    bitmap = 1 << 120
    verification = 1 << 10
    total_instructions = (
        steps * STEP_INSTRUCTIONS + (preprocessing + bitmap + verification) * C_INSTRUCTIONS
    )
    time_log2 = math.log2(total_instructions) - math.log2(C_INSTRUCTIONS)
    memory = (1 << 128) * 128 + (1 << 125) + (1 << 27) + (1 << 24)
    return {
        "r": str(R_ORDER),
        "p": str(P_FIELD),
        "classes": str(classes),
        "speed_constant": SPEED_CONSTANT,
        "steps": str(steps),
        "step_instructions": STEP_INSTRUCTIONS,
        "preprocessing_ref_ops": str(preprocessing),
        "bitmap_ref_ops": str(bitmap),
        "verification_ref_ops": str(verification),
        "total_instructions": str(total_instructions),
        "time_log2": time_log2,
        "submitted_time_log2": math.ceil(time_log2 * 1_000_000) / 1_000_000,
        "success_probability_heuristic": 1 - math.exp(-(steps * steps) / (2 * classes)),
        "memory_bytes": str(memory),
        "memory_log2_bytes": math.log2(memory),
        "preprocessing_log2": math.log2(preprocessing),
        "nonuniform_advice_bytes": "0",
    }


def main():
    try:
        config = json.load(sys.stdin)
        order = int(config["order"])
        trials = int(config["trials"])
        seed = int(config["seed"])
        increments = int(config["increments"])
        if not (7 <= order <= 1_000_003):
            raise ValueError("order out of range")
        if not (1 <= trials <= 500):
            raise ValueError("trials out of range")
        if not (0 <= seed <= (1 << 31) - 1):
            raise ValueError("seed out of range")
        if not (2 <= increments <= 4096):
            raise ValueError("increments out of range")
    except Exception as exc:
        print(json.dumps({"error": "invalid input", "detail": str(exc)}, sort_keys=True))
        return

    toy = orbit_walk(order, trials, seed, increments)
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
        "toy_ratio_matches_ideal": abs(toy["mean_ratio_to_sqrt_classes"] - toy["ideal_ratio"]) < 0.15,
        "toy_success_fraction_within_sampling_error": abs(
            toy["success_fraction_at_budget"] - toy["ideal_success_at_budget"]
        ) < 0.08,
        "ledger_total_below_2^128_ref_ops": int(ledger["total_instructions"]) < (1 << 128) * C_INSTRUCTIONS,
        "ledger_memory_below_2^136_bytes": int(ledger["memory_bytes"]) < (1 << 136),
        "success_probability_above_0.39": ledger["success_probability_heuristic"] > 0.39,
    }
    output = {
        "scope": (
            "Bounded supporting evidence only. The toy group checks the class-walk algebra "
            "and birthday scaling at small order; the profile check verifies the secp256k1 "
            "endomorphism constants; the ledger recomputes the analytic integers. No "
            "production attack is run and no AArch64 instruction count is measured."
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
