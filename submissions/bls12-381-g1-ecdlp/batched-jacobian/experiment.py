"""Toy Cheon identity check plus the production integer bounds.

Stdin JSON: {"modulus": <prime>, "divisor": <d dividing modulus-1>}.
The multiplicative group mod `modulus` stands in for a cyclic group of that
order. This does not price BLS12-381 additions and does not recover a
production scalar. It also verifies the width-5 wNAF nonzero-digit bound for
every scalar the production recurrence multiplies by.
"""

import json
import sys

R = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
P = 0x1A0111EA397FE69A4B1BA7B6434BACD764774B84F38512BF6730D2A0F6B0F6241EABFFFEB153FFFFB9FEFFFFFFFFAAAB
D_SCORE = 1 << 20
M_SCORE = 223621884914696853578958571250384896
LIMB = 1 << 64
RECORD_BYTES = 112
INDEX_OFFSET = 96
WNAF_W = 5
WNAF_NONZERO_BOUND = 64


def ceil_sqrt(n: int) -> int:
    if n < 0:
        raise ValueError("negative")
    if n <= 1:
        return n
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            break
        x = y
    return x if x * x == n else x + 1


def prime_factors(n: int) -> list[int]:
    factors = []
    x = n
    d = 2
    while d * d <= x:
        if x % d == 0:
            factors.append(d)
            while x % d == 0:
                x //= d
        d += 1 if d == 2 else 2
    if x > 1:
        factors.append(x)
    return factors


def is_primitive_root(g: int, modulus: int, factors: list[int]) -> bool:
    return all(pow(g, (modulus - 1) // q, modulus) != 1 for q in factors)


def pack_index(j: int) -> bytes:
    if j < 0 or j >= LIMB * LIMB:
        raise ValueError("index does not fit in two limbs")
    return (j & (LIMB - 1)).to_bytes(8, "little") + (j >> 64).to_bytes(8, "little")


def unpack_index(raw: bytes) -> int:
    if len(raw) != 16:
        raise ValueError("index width")
    return int.from_bytes(raw[:8], "little") | (int.from_bytes(raw[8:], "little") << 64)


def wnaf_digits(k: int, w: int) -> list[int]:
    """Width-w non-adjacent form of k as a little-endian digit list."""
    digits = []
    mod = 1 << w
    half = 1 << (w - 1)
    while k >= 1:
        if k & 1:
            digit = k % mod
            if digit >= half + 1:
                digit -= mod
            k -= digit
        else:
            digit = 0
        digits.append(digit)
        k >>= 1
    return digits


def check_wnaf(k: int) -> None:
    digits = wnaf_digits(k, WNAF_W)
    nonzero = sum(1 for d in digits if d != 0)
    if nonzero > WNAF_NONZERO_BOUND:
        raise AssertionError("wnaf density")
    if sum(d << i for i, d in enumerate(digits)) != k:
        raise AssertionError("wnaf reconstruction")
    for i, d in enumerate(digits):
        if d != 0:
            if d % 2 == 0 or abs(d) >= (1 << (WNAF_W - 1)) + 1:
                raise AssertionError("wnaf digit range")
            if any(digits[i + j] != 0 for j in range(1, WNAF_W) if i + j < len(digits)):
                raise AssertionError("wnaf nonadjacency")


def recover_all(modulus: int, divisor: int) -> None:
    if modulus < 3 or (modulus - 1) % divisor != 0:
        raise ValueError("divisor must divide modulus-1")
    order = modulus - 1
    m = order // divisor
    factors = prime_factors(order)
    zeta = 2
    while zeta < modulus and not is_primitive_root(zeta, modulus, factors):
        zeta += 1
    if zeta >= modulus:
        raise AssertionError("no primitive root")
    width = ceil_sqrt(m)
    giants = (m + width - 1) // width
    width2 = ceil_sqrt(divisor)
    giants2 = (divisor + width2 - 1) // width2
    lam = pow(zeta, -divisor, modulus)
    rho = pow(zeta, divisor * width, modulus)
    base = pow(zeta, m, modulus)
    lam2 = pow(base, -1, modulus)
    rho2 = pow(base, width2, modulus)
    for secret in range(1, modulus):
        babies = {}
        acc = pow(secret, divisor, modulus)
        for u in range(width):
            babies[acc] = u
            acc = acc * lam % modulus
        acc = 1
        hit = None
        for v in range(giants):
            if acc in babies:
                hit = (babies[acc], v)
                break
            acc = acc * rho % modulus
        if hit is None:
            raise AssertionError("k0 miss")
        k0 = (hit[0] + width * hit[1]) % m
        t2 = secret * pow(zeta, -k0, modulus) % modulus
        babies = {}
        acc = t2
        for u in range(width2):
            babies[acc] = u
            acc = acc * lam2 % modulus
        acc = 1
        hit = None
        for v in range(giants2):
            if acc in babies:
                hit = (babies[acc], v)
                break
            acc = acc * rho2 % modulus
        if hit is None:
            raise AssertionError("k1 miss")
        k1 = (hit[0] + width2 * hit[1]) % divisor
        if pow(zeta, k0 + m * k1, modulus) != secret:
            raise AssertionError("reconstruction")


def production_integers() -> None:
    m = (R - 1) // D_SCORE
    if (M_SCORE - 1) * (M_SCORE - 1) >= m or M_SCORE * M_SCORE < m:
        raise AssertionError("search width")
    if (m + M_SCORE - 1) // M_SCORE != M_SCORE:
        raise AssertionError("giant count")
    if not (2**117 < M_SCORE < 2**118):
        raise AssertionError("index width")
    raw = pack_index(M_SCORE - 1)
    if unpack_index(raw) != M_SCORE - 1:
        raise AssertionError("limb roundtrip")
    if len(raw) != 16 or INDEX_OFFSET + len(raw) != RECORD_BYTES:
        raise AssertionError("record layout")
    factors = [2, 3, 11, 19, 10177, 125527, 859267, 906349, 2508409, 2529403, 52437899, 254760293]
    prod = (
        (1 << 32)
        * 3
        * 11
        * 19
        * 10177
        * 125527
        * 859267
        * 906349
        * 906349
        * 2508409
        * 2529403
        * 52437899
        * 254760293
        * 254760293
    )
    if prod != R - 1:
        raise AssertionError("factor product")
    if not is_primitive_root(7, R, factors):
        raise AssertionError("primitive root 7")
    # Algebraic split r-1 = X^2 (X-1) (X+1), X = -parameter.
    x = 0xD201000000010000
    if x * x * (x - 1) * (x + 1) != R - 1:
        raise AssertionError("parameter factorization")
    if P.bit_length() != 381:
        raise AssertionError("field width")
    # Width-5 wNAF density of every scalar a production step multiplies by.
    # The ledger charges 64 nonzero digits per 255-bit scalar.
    stage_scalars = [
        pow(7, -D_SCORE, R),          # lambda  = 7^{-d}
        pow(7, D_SCORE * M_SCORE, R), # rho     = 7^{dM}
        pow(7, -m, R),                #           7^{-m}
        pow(7, m * (1 << 10), R),     #           7^{m M_2}
        R - 1,                        # worst-case runtime scalar shape
    ]
    for scalar in stage_scalars:
        if scalar.bit_length() > 255:
            raise AssertionError("scalar length")
        check_wnaf(scalar)


def main() -> None:
    spec = json.load(sys.stdin)
    recover_all(int(spec["modulus"]), int(spec["divisor"]))
    production_integers()
    print(
        json.dumps(
            {
                "checked_modulus": spec["modulus"],
                "checked_divisor": spec["divisor"],
                "all_recovered": True,
                "production_integers_ok": True,
                "wnaf5_nonzero_bound": WNAF_NONZERO_BOUND,
                "scope": "toy multiplicative group plus integer and recoding bounds; not a curve timing",
            }
        )
    )


if __name__ == "__main__":
    main()
