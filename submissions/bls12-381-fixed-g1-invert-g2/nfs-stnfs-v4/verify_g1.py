"""End-to-end check of the fixed-G1 reduction with the pinned library.

This is *not* an attack run and it makes no claim about STNFS.  It checks the
algebra that the submission's correctness statement rests on, at production
field size, using the pinned `py_ecc.optimized_bls12_381` library and the
encoding/acceptance procedure of `tracks/bls12-381-fixed-g1-invert-g2/
target-profile.json`:

    P = [b]P0,  z = e(P,Q0)^a,  solver recovers a = log_{e(P,Q0)} z
    (here a is *given* to the checker, not computed),  Q = [a]Q0,
    acceptance: e(P,Q) == z with all 12 GT coefficients compared.

It also covers the two edge cases the claim names, `z = 1` (a = 0, Q = O) and
`b` near the top of the scalar range.  Deterministic: seed 20261023.

Run:  .venv/bin/python verify_g1.py   (writes toy-g1.json next to itself)
"""
import hashlib
import json
import os
import random

from py_ecc.fields import optimized_bls12_381_FQ2 as FQ2
from py_ecc.optimized_bls12_381 import (FQ12, G1, G2, b2, curve_order, is_inf,
                                        is_on_curve, multiply, normalize, pairing)

R = curve_order
SEED = 20261023
P_HEX = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab


def coeffs(gt):
    return [int(c.n) if hasattr(c, "n") else int(c) for c in gt.coeffs]


def enc_g2(pt):
    """target-profile.json G2 encoding: affine coordinates, Fp2 as [c0, c1]."""
    if is_inf(pt):
        return None
    xa, ya = normalize(pt)
    return {"x": [format(int(c), "096x") for c in xa.coeffs],
            "y": [format(int(c), "096x") for c in ya.coeffs]}


def dec_g2(obj):
    if obj is None:
        return None
    return (FQ2([int(s, 16) for s in obj["x"]]),
            FQ2([int(s, 16) for s in obj["y"]]),
            FQ2.one())


def on_g2(pt):
    """y^2 = x^3 + 4*(1+u) over Fp2, via the pinned library's own check."""
    if pt is None:
        return True
    return is_on_curve(pt, b2)


def one_instance(b, a, label):
    P = multiply(G1, b)
    g = pairing(G2, P, final_exponentiate=True)          # e(P,Q0): G2 argument first
    z = g ** a                                            # e(P,Q0)^a
    Q = multiply(G2, a)                                   # [a]Q0, O when a = 0
    # --- acceptance procedure of target-profile.json ---------------------
    enc = enc_g2(Q)
    Q_dec = dec_g2(enc)
    ok_curve = on_g2(Q_dec)
    ok_order = (Q_dec is None) if a == 0 else (Q_dec is not None and is_inf(multiply(Q_dec, R)))
    lhs = pairing(Q_dec, P, final_exponentiate=True) if Q_dec is not None else FQ12.one()
    ok_pair = coeffs(lhs) == coeffs(z)
    ok_12 = len(coeffs(lhs)) == 12 and len(coeffs(z)) == 12
    ok_identity = (a != 0) or (enc is None and coeffs(z) == coeffs(FQ12.one()))
    return {
        "label": label,
        "b_hex": format(b, "064x"),
        "a_hex": format(a, "064x"),
        "z_sha256_12coeff": hashlib.sha256(
            b"".join(int(c).to_bytes(48, "big") for c in coeffs(z))).hexdigest(),
        "Q_is_infinity": enc is None,
        "pairing_matches_target_12coefficients": ok_pair,
        "returned_Q_on_curve": ok_curve,
        "returned_Q_in_subgroup_r": ok_order,
        "identity_case_consistent": ok_identity,
        "gt_element_has_12_coefficients": ok_12,
    }


def main():
    rng = random.Random(SEED)
    cases = [("random", rng.randrange(1, R), rng.randrange(0, R)) for _ in range(5)]
    cases.append(("b=1", 1, rng.randrange(0, R)))
    cases.append(("b=r-1", R - 1, rng.randrange(0, R)))
    cases.append(("z=1 (a=0)", rng.randrange(1, R), 0))

    rows = [one_instance(b, a, label) for label, b, a in cases]

    # the supplied P must be the second pairing argument in this library
    P0m = multiply(G1, 12345)
    order_check = None
    try:
        pairing(P0m, multiply(G2, 7), final_exponentiate=True)
        order_check = "accepted G1 first (unexpected)"
    except Exception as exc:                       # noqa: BLE001 - pinning check
        order_check = type(exc).__name__

    checks = {
        "all_pairings_match_target": all(r["pairing_matches_target_12coefficients"] for r in rows),
        "all_returned_points_on_curve": all(r["returned_Q_on_curve"] for r in rows),
        "all_returned_points_in_subgroup": all(r["returned_Q_in_subgroup_r"] for r in rows),
        "all_identity_cases_consistent": all(r["identity_case_consistent"] for r in rows),
        "all_gt_elements_12_coefficients": all(r["gt_element_has_12_coefficients"] for r in rows),
        "pinned_argument_order_is_pairing_G2_first": order_check in
        ("TypeError", "ValueError", "AssertionError", "IndexError"),
    }
    out = {
        "scope": ("algebraic end-to-end check of the reduction at production field "
                  "size with the pinned py_ecc library; a is supplied to the "
                  "checker, no discrete logarithm is computed, and no complexity "
                  "claim is supported by this file"),
        "library": "py_ecc.optimized_bls12_381 (py_ecc 8.0.0)",
        "curve_order_bits": R.bit_length(),
        "field_prime_hex": format(P_HEX, "0x"),
        "seed": SEED,
        "instances": rows,
        "pinned_argument_order_probe": order_check,
        "checks": checks,
    }
    failed = [k for k, v in checks.items() if not v]
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "toy-g1.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")
    print(json.dumps({"checks": checks, "instances": len(rows)}, indent=1))
    if failed:
        raise SystemExit("failed checks: " + ", ".join(failed))


if __name__ == "__main__":
    main()
