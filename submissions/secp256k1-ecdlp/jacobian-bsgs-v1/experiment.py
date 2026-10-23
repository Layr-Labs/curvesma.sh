"""Exact word-operation ledger for the Jacobian-coordinate BSGS on secp256k1.

Two independent pieces of evidence are emitted:

  * a toy exhaustive check of scalar recovery on a small additive cyclic
    group, bound to the claim's input schema, and
  * the exact integer ledger of the production construction, recomputed
    from the formulas proof.md derives.

The production search is analytic; this module never executes it.
"""
import json
import math
import sys

N = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141

# classical-group-work-v2 primitives
K = 256
S_LIMBS = 1
A = 128 * (S_LIMBS + 1)
M = K * (2 * A + 64) + 128 * (S_LIMBS + 1)
I = 2 * K * (M + 64) + 128 * (S_LIMBS + 1)
C = I + 16 * M + 32 * A + 4096

# projective kernels
J = 12 * M + 16 * A + 4096      # Jacobian doubling / mixed addition
Q = 6 * M + 4 * A + 4096        # Jacobian-vs-affine equality probe
BN = 3 * M + 64                 # Montgomery batch normalisation per point

W = 4
B = 16 * 3 + 64
D = 16 * W + 64
L = 128
V = (L + 1) * (B + 64) + D + 256

SETUP = 2 ** 24
VERIFY = 4 * K * C + 2 ** 24
PRICE = C


def extents():
    """Baby-table width m and giant-walk length h covering the whole order."""
    m = 1 << L
    h = -(-N // m)
    return m, h


def ledger():
    m, h = extents()
    baby = m * (J + BN + D) + I
    sort = m * L * (B + 2 * D + 128) + 4096 * L
    giant = h * (J + Q + V)
    total = baby + sort + giant + SETUP + VERIFY
    return {
        "m": m,
        "h": h,
        "coverage_num": m * h,
        "coverage_den": N,
        "baby_table": baby,
        "merge_sort": sort,
        "giant_search": giant,
        "setup": SETUP,
        "verification": VERIFY,
        "preprocessing": baby + sort,
        "total": total,
        "price": PRICE,
        "time_log2": math.log2(total / PRICE),
        "preprocessing_log2": math.log2((baby + sort) / PRICE),
        "peak_bytes": 2 * m * W * 32 + 2 ** 20,
    }


def toy(order):
    """Exhaustive BSGS scalar recovery in the additive group modulo order."""
    m = math.isqrt(order) + 1
    babies = {j % order: j for j in range(m)}
    for secret in range(order):
        value = secret
        for i in range((order + m - 1) // m):
            if value in babies:
                recovered = (i * m + babies[value]) % order
                assert recovered == secret
                break
            value = (value - m) % order
        else:
            raise AssertionError("Discrete logarithm not recovered")


def main():
    req = json.load(sys.stdin)
    order = req["order"]
    toy(order)
    row = ledger()
    assert row["coverage_num"] >= N
    assert row["time_log2"] <= 124.1
    assert row["preprocessing_log2"] <= 122.0
    out = {
        "group": "integers modulo order under addition",
        "order": order,
        "checked": order,
        "all_recovered": True,
        "scope": "toy group only; production bounds are analytic",
        "ledger": {
            key: (str(value) if isinstance(value, int) and value > 2 ** 53 else value)
            for key, value in row.items()
        },
    }
    print(json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
