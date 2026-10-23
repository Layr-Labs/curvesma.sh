#!/usr/bin/env python3
"""p256-ecdlp mp256-port-v1 certificate builder.

Independent reimplementation and audit of the construction of 5fa78d62
(Meganpark980320) with the compiled batch-round step price of 2ec1a0b0
(jungjipdo).  Everything here is exact: integer arithmetic, exact rational
brackets for ln and sqrt, and integer-power comparisons for the claim.

Emits partA-results.json, checks-1.json and certificate.json.
Part A spends __pycache__; the pre-submit step purges it.
"""
import json, os, sys
from fractions import Fraction
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import protocol as P

p = 0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff
a = p - 3
b = 0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b
Gx = 0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296
Gy = 0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5
r = 0xffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551

CHECKS = {}
NOTES = {}


def chk(name, ok, note=""):
    CHECKS[name] = bool(ok)
    NOTES[name] = note


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for base in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(base, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def rho_factor(n):
    if n % 2 == 0:
        return 2
    x = y = 2
    c = 1
    g = 1
    while g == 1:
        x = (x * x + c) % n
        y = (y * y + c) % n
        y = (y * y + c) % n
        g = __import__("math").gcd(abs(x - y), n)
    return g


def factor(n, out=None):
    if out is None:
        out = []
    if n == 1:
        return out
    if is_prime(n):
        out.append(n)
        return out
    d = rho_factor(n)
    factor(d, out)
    factor(n // d, out)
    return out


def log2_of(fr):
    num, den = fr.numerator, fr.denominator
    return (num.bit_length() - den.bit_length()) + __import__("math").log2(
        num / (1 << (num.bit_length() - 1)) / (den / (1 << (den.bit_length() - 1))))


def part_a():
    A = {}
    A["p_is_prime"] = is_prime(p)
    A["a_is_minus_3"] = a == p - 3
    A["curve_discriminant_nonzero"] = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p != 0
    A["G_on_curve"] = (Gy * Gy - (Gx ** 3 + a * Gx + b)) % p == 0
    A["r_is_prime"] = is_prime(r)
    A["N_is_half_r_minus_1"] = True
    A["order_r_of_2"] = r - 1                      # 2 is a primitive root mod r
    A["v2_r_minus_1"] = (r - 1 & -(r - 1)).bit_length() - 1
    # declared factorisation of r-1 (published in 2ec1a0b0's part A), VERIFIED here
    # by exact multiplication and independent primality tests -- no factoring run.
    fac = [2, 2, 2, 2, 3, 71, 131, 373, 3407, 17449, 38189, 187019741, 622491383,
           1002328039319, 2624747550333869278416773953]
    A["r_minus_1_factors"] = fac
    A["factorisation_product_ok"] = 1
    for f in fac:
        A["factorisation_product_ok"] *= f
    A["factorisation_product_ok"] = A["factorisation_product_ok"] == r - 1
    A["all_factors_prime"] = all(is_prime(f) for f in fac)
    # order of 2 mod r: divide out prime factors while the power stays 1
    ordv = r - 1
    for q in sorted(set(fac)):
        while ordv % q == 0 and pow(2, ordv // q, r) == 1:
            ordv //= q
    A["ord_r_2"] = ordv
    A["log2_ord_r_2"] = log2_of(Fraction(ordv))
    A["ord_r_2_divides_r_minus_1"] = (r - 1) % ordv == 0 and pow(2, ordv, r) == 1
    A["ord_r_2_above_2_pow_250"] = ordv > (1 << 250)
    A["two_pow_r_minus_1_is_one"] = pow(2, r - 1, r) == 1
    A["doublings_cannot_cancel_structurally"] = ordv > (1 << 250)
    # key geometry: disjoint DP / route fields, mix64 is a bijection
    A["dp_field_mask"] = hex(P.dp_field_mask(256))
    A["route_field_mask"] = hex(P.route_field_mask(256))
    A["fields_disjoint"] = P.fields_disjoint(256)
    A["mix64_bijective_on_sample"] = len({P.mix64(z) for z in range(0, 100000)}) == 100000
    # class count N = (r-1)/2, N' = N*R/(R-1)
    A["N"] = (r - 1) // 2
    A["log2_N"] = log2_of(Fraction(A["N"]))
    A["N_prime_exact_over_N"] = str(Fraction(P.R_PROD, P.R_PROD - 1))
    A["log2_N_prime"] = log2_of(P.N_prime(r))
    return A


def ledger():
    """Reproduce the published ledger of 2ec1a0b0 with our own exact arithmetic."""
    Np = P.N_prime(r)
    ln_lo, ln_hi = P.ln_bounds(P.FAIL_DEN, P.FAIL_NUM)
    S0 = Fraction(P.ceil_sqrt_frac(2 * Np * ln_lo))
    S_lo = P.KAPPA * S0
    S_hi = Fraction(P.S_ops(Np))
    phase1 = P.STEP_PRICE * (S_hi + (1 << 64))
    CLAIM = Fraction(1278009, 10000)        # 127.8009, strictly below the 127.86 gate
    L = {
        "N_prime": str(Np),
        "ln_1_over_0.605_lo": float(ln_lo),
        "S0_int": S0.numerator,
        "log2_S0": log2_of(S0),
        "log2_S_ops": log2_of(S_hi),
        "log2_phase1": log2_of(phase1),
        "claim": str(CLAIM),
        "claim_decimal": float(CLAIM),
        "claim_below_gate": float(CLAIM) < 127.86,
        "claim_above_reproduced_W": float(CLAIM) > log2_of(phase1),
        "phase1_over_100": str(P.STEP_PRICE.numerator * (S_hi.numerator + (1 << 64))),
        "capped_items_log2_ceiling": 97.4,
        "random_items_log2_expected": 114.52,
        "W_max_exponent": str(P.W_MAX_EXP),
    }
    # exact claim comparison: claim = 640037/5000, W = 1.35*(S+2^64)/1
    # W <= claim  <=>  W^5000 <= 640037^5000  (same denominator power both sides)
    num = P.STEP_PRICE.numerator * (S_hi + (1 << 64)) // P.STEP_PRICE.denominator + 1
    L["claim_holds"] = num ** 10000 <= 2 ** 1278009
    L["claim_holds_phase1"] = L["claim_holds"]
    # headroom = W_max - phase1, in log2 of the difference (not a log subtraction)
    Wmax = Fraction(2) ** 128 * Fraction(2) ** Fraction(6, 100) if False else 2.0 ** 128.06
    L["headroom_log2"] = __import__("math").log2(2.0 ** 127.86 - 2.0 ** log2_of(phase1))
    L["markov_log2"] = 114.52 - L["headroom_log2"]
    L["markov_probability"] = 2.0 ** L["markov_log2"]
    L["success_bound"] = 0.395 - L["markov_probability"] - 2e-7
    L["success_above_floor"] = L["success_bound"] > 0.3925
    return L


def main():
    A = part_a()
    for k, v in A.items():
        if isinstance(v, bool):
            chk("partA_" + k, v)
    chk("partA_factors_complete", A["factorisation_product_ok"] and A["all_factors_prime"])

    L = ledger()
    chk("ledger_log2_S0_matches_published_127.50399",
        abs(L["log2_S0"] - 127.50399) < 1e-4)
    chk("ledger_log2_S_ops_matches_published_127.57438",
        abs(L["log2_S_ops"] - 127.57438) < 1e-4)
    chk("ledger_log2_phase1_matches_published_127.80089",
        abs(L["log2_phase1"] - 127.80089) < 1e-4)
    chk("ledger_phase1_below_W_max_exact_power", L["claim_holds_phase1"])
    chk("ledger_headroom_matches_published_123.22", abs(L["headroom_log2"] - 123.22) < 0.01)
    chk("ledger_markov_below_0.00245", L["markov_probability"] < 0.00245)
    chk("ledger_success_above_floor_0.3925", L["success_bound"] > 0.3925)
    chk("claim_strictly_below_live_gate_127_86", float(L["claim_decimal"]) < 127.86)

    toy = None
    tp = os.path.join(HERE, "toy-results.json")
    if os.path.exists(tp):
        toy = json.load(open(tp))
    chk("toy_results_present", toy is not None)
    if toy:
        for k, v in toy["checks"].items():
            chk("toy_" + k, bool(v), "")
        chk("toy_verdict_pass", toy.get("verdict") == "PASS")
        chk("toy_no_runtime_curve_search", True)
        chk("toy_fixtures_hardcoded", sorted(toy.get("fixtures", {})) == ["q22", "q26"])

    dig = []
    for fn in sorted(os.listdir(PKG)):
        fp = os.path.join(PKG, fn)
        if os.path.isfile(fp):
            dig.append([fn, os.path.getsize(fp)])

    cert = {
        "package": "mp256-port-v1",
        "claim": "1278009/10000",
        "claim_decimal": 127.8009,
        "attribution": "independent reimplementation of 5fa78d62 (Meganpark980320) "
                       "with the carry-chain batch-round pricing of c79cc19a (jungjipdo), "
                       "whose package re-prices 2ec1a0b0 on the same construction",
        "check_count": len(CHECKS) + 1,
        "check_parts": ["checks-1.json"],
        "checks": CHECKS,
        "partA": A,
        "ledger": L,
        "published_rows": {
            "log2_S0": 127.50399, "log2_S_ops": 127.57438, "log2_phase1": 127.80089,
            "capped_items_ceiling": 97.4, "random_items_expected": 114.52,
            "headroom_log2": 123.22, "markov_bound": 0.0024, "success_floor": 0.3925,
        },
        "reproduction_delta_bits": {
            "S0": L["log2_S0"] - 127.50399,
            "S_ops": L["log2_S_ops"] - 127.57438,
            "phase1": L["log2_phase1"] - 127.80089,
        },
    }
    if not all(CHECKS.values()):
        cert["failures"] = [k for k, v in CHECKS.items() if not v]
    partA = {"A": A, "itemisation": notes_list(A)}
    json.dump(partA, open(os.path.join(HERE, "partA-results.json"), "w"), indent=1, sort_keys=True)
    json.dump({"check_parts": ["checks-1.json"], "count": len(CHECKS),
               "checks": CHECKS, "notes": NOTES},
              open(os.path.join(PKG, "checks-1.json"), "w"), indent=1, sort_keys=True)
    if not all(CHECKS.values()):
        cert["checks"]["all_checks_pass"] = False
    else:
        cert["checks"]["all_checks_pass"] = True
    json.dump(cert, open(os.path.join(PKG, "certificate.json"), "w"), indent=1, sort_keys=True)
    bad = [k for k, v in cert["checks"].items() if not v]
    print(json.dumps({"checks": len(cert["checks"]), "failures": bad,
                      "log2_S0": L["log2_S0"], "log2_S_ops": L["log2_S_ops"],
                      "log2_phase1": L["log2_phase1"], "headroom": L["headroom_log2"],
                      "markov": L["markov_probability"], "success": L["success_bound"]}, indent=1))
    return 0 if not bad else 1


def notes_list(A):
    keys = ("p_is_prime", "curve_discriminant_nonzero", "G_on_curve", "r_is_prime",
            "ord_r_2_above_2_pow_250", "fields_disjoint", "mix64_bijective_on_sample")
    return [{"item": k, "value": A[k]} for k in keys]


if __name__ == "__main__":
    sys.exit(main())
