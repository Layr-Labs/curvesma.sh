"""Toy witness for half-range BSGS with a cyclotomic odd-part key.

Production group is BLS12-381 GT in Fp[w]/(w^12-2*w^6+2). Frobenius-6
negates odd monomial coefficients and, on the order-r subgroup, equals
inversion because p^6 ≡ -1 (mod r). The lookup key keeps those six odd
coefficients, canonicalised under negation, plus one bit of the even part.
That bit separates the two roots of X^2 - c in a field. This process never
builds Fp12. It checks the public identities, checks the coverage integers,
exhausts the stride on Z/nZ, and exhausts the same key on the norm-1 set
of a small prime field.
"""

import json
import math
import sys

P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
C = 20000
# Schedules in proof.md, each already doubled.
B = 274
D = 184
ADDR = 64
CONJ = 576
SIGN = 224
CTRL = 96
CANON = CONJ + SIGN + B + D
PASS = B + 2 * D + CTRL
L = 126
M = 1 << L
REC = 320
FIELD_TOY = 1009


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
    h_misc = 9216
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


def norm1_keys(prime: int) -> int:
    """Norm-1 keys in F_p. X^2 - c has at most two roots; the low bit separates them."""
    if prime < 3 or prime % 2 == 0 or prime > 10000:
        raise SystemExit("field toy out of range")
    for a in range(2, int(math.isqrt(prime)) + 1):
        if prime % a == 0:
            raise SystemExit("field toy not prime")

    def neg(x: int) -> int:
        return (prime - x) % prime

    def sign(a: int) -> int:
        if a == 0:
            return 0
        return a & 1

    for a in range(prime):
        other = neg(a)
        if a == other:
            if a != 0:
                raise SystemExit("sign fixed point")
            continue
        if sign(a) == sign(other):
            raise SystemExit("sign does not separate roots")

    pairs = {}
    counted = 0
    for b in range(prime):
        c = (1 + (b * b) % prime) % prime
        roots = [a for a in range(prime) if (a * a) % prime == c]
        if len(roots) > 2:
            raise SystemExit("degree-2 polynomial with 3 roots")
        canon = b if b <= neg(b) else neg(b)
        for a in roots:
            key = (canon, sign(a))
            pair = frozenset(((a, b), (a, neg(b))))
            previous = pairs.get(key)
            if previous is None:
                pairs[key] = pair
            elif previous != pair:
                raise SystemExit("norm-1 key collision")
            counted += 1
    # a^2 - b^2 = 1 has p - 1 solutions: a - b runs through F_p^x.
    if counted != prime - 1:
        raise SystemExit("norm-1 cardinality")
    return counted


def main() -> None:
    payload = json.load(sys.stdin)
    order = payload["order"]
    exhaust(order)
    frob = frobenius6_is_odd_negation()
    if not frob:
        raise SystemExit("frobenius")
    counted = norm1_keys(FIELD_TOY)
    ledger = production_ledger()
    json.dump(
        {
            "group": "integers modulo order under addition",
            "order": order,
            "checked": order,
            "all_recovered": True,
            "frobenius6_odd_negation": frob,
            "p6_is_minus_one_mod_r": True,
            "norm1_field": FIELD_TOY,
            "norm1_elements": counted,
            "norm1_keys_injective": True,
            "production": {
                "M": str(ledger["M"]),
                "stride": str(ledger["stride"]),
                "n_giants": str(ledger["n_giants"]),
                "interval_length": str(ledger["interval_length"]),
                "lookup_each": str(ledger["lookup_each"]),
                "W": str(ledger["W"]),
                "W_pre": str(ledger["W_pre"]),
                "memory_bytes": str(ledger["memory_bytes"]),
            },
            "scope": "toy group, small-field norm-1 keys, exact public constants; not a production run",
        },
        sys.stdout,
    )
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
