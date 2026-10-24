"""Toy witness and integer ledger for fixed-Q half-range BSGS.

Production group is BLS12-381 GT. Frobenius^6 negates odd monomial
coefficients and, on the order-r subgroup, equals inversion because
p^6 ≡ -1 (mod r). This process never builds Fp12. It checks that
constant, checks the coverage integers, and exhausts the same control
flow on Z/nZ with negation standing in for π^6.
"""

import json
import math
import sys

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
C = 20000
# Concrete AArch64 schedules times two. See proof.md.
B = 512
D = 320
ADDR = 64
CONJ = 576
CTRL = 96
CANON = CONJ + B + D
PASS = B + 2 * D + CTRL
L = 126
M = 1 << L
REC = 592


def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


def poly_mul(a, b):
    c = [0] * 24
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj:
                c[i + j] = (c[i + j] + ai * bj) % P
    for deg in range(22, 11, -1):
        coef = c[deg]
        if coef == 0:
            continue
        c[deg] = 0
        c[deg - 6] = (c[deg - 6] + 2 * coef) % P
        c[deg - 12] = (c[deg - 12] - 2 * coef) % P
    return c[:12]


def poly_pow_w(exp: int):
    result = [1] + [0] * 11
    base = [0, 1] + [0] * 10
    while exp:
        if exp & 1:
            result = poly_mul(result, base)
        base = poly_mul(base, base)
        exp >>= 1
    return result


def frobenius6_is_odd_negation() -> bool:
    # w^12 = 2 w^6 - 2, so (w^6)^2 - 2 w^6 + 2 = 0.
    w6 = poly_pow_w(6)
    w6_sq = poly_mul(w6, w6)
    expect = [(-2) % P, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0]
    if w6 != [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]:
        return False
    if w6_sq != expect:
        return False
    image = poly_pow_w(P**6)
    want = [0, (P - 1) % P] + [0] * 10
    return image == want


def production_ledger():
    if pow(P, 6, R) != R - 1:
        raise SystemExit("p^6 != -1 mod r")
    stride = 2 * (M - 1)
    n = ceil_div(R - 1, stride)
    if 2 * n * (M - 1) + 1 < R:
        raise SystemExit("coverage short")
    if 2 * (M - 1) >= R:
        raise SystemExit("baby pairs wrap")
    h_mul = M + n - 1
    h_misc = 9216  # two pairings at 4096 and z^r at 1024
    sort_i = L * M * PASS
    baby_canon_i = M * CANON
    store_i = M * D
    v = CANON + (L + 1) * (B + ADDR)
    look_i = n * v
    r_instr = sort_i + baby_canon_i + store_i + look_i + CONJ + (1 << 26)
    w = h_mul + h_misc + (r_instr + C - 1) // C
    r_pre = sort_i + baby_canon_i + store_i + CONJ + (1 << 20)
    w_pre = M + 4096 + (r_pre + C - 1) // C
    mem = 2 * M * REC + (1 << 20)
    return {
        "M": M,
        "stride": stride,
        "n_giants": n,
        "interval_length": 2 * n * (M - 1) + 1,
        "H_mul": h_mul,
        "H_misc": h_misc,
        "sort_instructions": sort_i,
        "canon_instructions": baby_canon_i,
        "store_instructions": store_i,
        "lookup_instructions": look_i,
        "lookup_each": v,
        "R": r_instr,
        "W": w,
        "W_pre": w_pre,
        "memory_bytes": mem,
    }


def choose_width(order: int) -> int:
    m = max(2, math.isqrt(order))
    while m > 2 and (m - 1) + (m - 2) >= order:
        m -= 1
    return m


def recover(secret: int, order: int, width: int) -> int:
    """Same control flow as production, additive generator 1."""
    stride = 2 * (width - 1)
    n = ceil_div(order - 1, stride)
    table = {}
    for j in range(width):
        neg = (-j) % order
        if (j % order) <= neg:
            key = j % order
            signed = j
        else:
            key = neg
            signed = -j
        if key in table:
            raise SystemExit("canonical key collision")
        table[key] = signed
    y = secret % order
    step = stride % order
    for i in range(n):
        neg = (-y) % order
        key = y if y <= neg else neg
        if key in table:
            signed = table[key]
            if y == signed % order:
                found = i * stride + signed
            else:
                found = i * stride - signed
            return found % order
        y -= step
        if y < 0:
            y += order
    raise SystemExit("miss")


def exhaust(order: int) -> None:
    if not isinstance(order, int) or order < 2 or order > 10000:
        raise SystemExit("order out of range")
    width = choose_width(order)
    for secret in range(order):
        got = recover(secret, order, width)
        if got != secret:
            raise SystemExit(f"mismatch {secret}->{got}")


def main() -> None:
    payload = json.load(sys.stdin)
    order = payload["order"]
    exhaust(order)
    frob = frobenius6_is_odd_negation()
    if not frob:
        raise SystemExit("frobenius")
    ledger = production_ledger()
    json.dump(
        {
            "group": "integers modulo order under addition",
            "order": order,
            "checked": order,
            "all_recovered": True,
            "frobenius6_odd_negation": frob,
            "p6_is_minus_one_mod_r": True,
            "production": {
                "M": str(ledger["M"]),
                "stride": str(ledger["stride"]),
                "n_giants": str(ledger["n_giants"]),
                "interval_length": str(ledger["interval_length"]),
                "W": str(ledger["W"]),
                "W_pre": str(ledger["W_pre"]),
                "memory_bytes": str(ledger["memory_bytes"]),
            },
            "scope": "toy group plus exact public constants; not a production run",
        },
        sys.stdout,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
