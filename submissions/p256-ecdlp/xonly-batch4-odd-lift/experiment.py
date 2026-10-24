"""Toy odd-lift BSGS, plus the integer arithmetic of the production bound.

The curve search and the discrete-log checks use a prime-order curve over a
17-element field. They do not price P-256. The production section only checks
the same integer inequalities written in proof.md: instruction prices times
operation counts, divided by C = 1000, sit under 2^133.
"""

from __future__ import annotations

import json
import sys

P256_P = 2**256 - 2**224 + 2**192 + 2**96 - 1
M_INSTR = 672
A_INSTR = 32
CMP_INSTR = 48
MOVE_INSTR = 32
BLOCK_LOG = 12  # giant x-conversion blocks of 2^12


def egcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return a, 1, 0
    g, x, y = egcd(b, a % b)
    return g, y, x - (a // b) * y


def modinv(a: int, p: int) -> int:
    g, x, _ = egcd(a % p, p)
    if g != 1:
        raise ZeroDivisionError(a)
    return x % p


Point = tuple[int, int] | None
Jac = tuple[int, int, int] | None


def affine_add(p: int, a: int, p1: Point, p2: Point) -> Point:
    if p1 is None:
        return p2
    if p2 is None:
        return p1
    x1, y1 = p1
    x2, y2 = p2
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if p1 == p2:
        if y1 % p == 0:
            return None
        lam = ((3 * x1 * x1 + a) * modinv(2 * y1, p)) % p
    else:
        lam = ((y2 - y1) * modinv(x2 - x1, p)) % p
    x3 = (lam * lam - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return x3, y3


def affine_mul(p: int, a: int, k: int, pt: Point) -> Point:
    r: Point = None
    q = pt
    while k:
        if k & 1:
            r = affine_add(p, a, r, q)
        q = affine_add(p, a, q, q)
        k >>= 1
    return r


def jac_double(p: int, a: int, pt: Jac) -> Jac:
    if pt is None:
        return None
    x1, y1, z1 = pt
    if z1 % p == 0 or y1 % p == 0:
        return None
    xx = (x1 * x1) % p
    yy = (y1 * y1) % p
    yyyy = (yy * yy) % p
    zz = (z1 * z1) % p
    s = (2 * ((x1 + yy) ** 2 - xx - yyyy)) % p
    m = (3 * xx + a * (zz * zz % p)) % p
    x3 = (m * m - 2 * s) % p
    y3 = (m * (s - x3) - 8 * yyyy) % p
    z3 = ((y1 + z1) ** 2 - yy - zz) % p
    return x3, y3, z3


def jac_add_affine(p: int, a: int, pt: Jac, q: Point) -> Jac:
    if q is None:
        return pt
    if pt is None or pt[2] % p == 0:
        return q[0], q[1], 1
    x1, y1, z1 = pt
    x2, y2 = q
    z1z1 = (z1 * z1) % p
    u2 = (x2 * z1z1) % p
    s2 = (y2 * z1 * z1z1) % p
    h = (u2 - x1) % p
    if h == 0:
        if s2 % p == y1 % p:
            return jac_double(p, a, pt)
        return None
    hh = (h * h) % p
    i = (4 * hh) % p
    j = (h * i) % p
    r = (2 * (s2 - y1)) % p
    v = (x1 * i) % p
    x3 = (r * r - j - 2 * v) % p
    y3 = (r * (v - x3) - 2 * y1 * j) % p
    z3 = ((z1 + h) ** 2 - z1z1 - hh) % p
    return x3, y3, z3


def batch_x(p: int, pts: list[Jac]) -> list[int | None]:
    out: list[int | None] = [None] * len(pts)
    finite = [(i, pt) for i, pt in enumerate(pts) if pt is not None and pt[2] % p != 0]
    k = len(finite)
    if k == 0:
        return out
    zs = [pt[2] % p for _, pt in finite]
    prefix = [1]
    for t in range(k - 1):
        prefix.append(prefix[-1] * zs[t] % p)
    acc = modinv(prefix[-1] * zs[-1] % p, p)
    invz = [0] * k
    for t in range(k - 1, 0, -1):
        invz[t] = acc * prefix[t] % p
        acc = acc * zs[t] % p
    invz[0] = acc
    for slot, (i, pt) in enumerate(finite):
        zinv2 = invz[slot] * invz[slot] % p
        out[i] = pt[0] * zinv2 % p
    return out


def odd_representative(k: int, n: int) -> int:
    if k % 2 == 1:
        return k
    return k + n


def unique_odds(m: int) -> bool:
    limit = 2 ** (2 * m + 1)
    seen: set[int] = set()
    for i in range(2**m):
        base = (2 * i + 1) << m
        for j in range(2 ** (m - 1)):
            step = 2 * j + 1
            for sigma in (1, -1):
                z = base + sigma * step
                if z <= 0 or z >= limit or z in seen or z % 2 == 0:
                    return False
                seen.add(z)
    return seen == set(range(1, limit, 2))


def mergesort(rows: list[tuple[int, int]]) -> tuple[list[tuple[int, int]], int]:
    comps = 0

    def merge(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> list[tuple[int, int]]:
        nonlocal comps
        i = j = 0
        out: list[tuple[int, int]] = []
        while i < len(a) and j < len(b):
            comps += 1
            if a[i][0] <= b[j][0]:
                out.append(a[i])
                i += 1
            else:
                out.append(b[j])
                j += 1
        out.extend(a[i:])
        out.extend(b[j:])
        return out

    def rec(xs: list[tuple[int, int]]) -> list[tuple[int, int]]:
        if len(xs) <= 1:
            return list(xs)
        mid = len(xs) // 2
        return merge(rec(xs[:mid]), rec(xs[mid:]))

    return rec(rows), comps


def mergesort_bound(n: int) -> int:
    if n <= 1:
        return 0
    ceil_log = (n - 1).bit_length()
    return n * ceil_log - (1 << ceil_log) + 1


def binary_search(table: list[tuple[int, int]], x: int) -> tuple[int | None, int]:
    lo, hi, comps = 0, len(table), 0
    while lo < hi:
        mid = (lo + hi) // 2
        comps += 1
        y = table[mid][0]
        if y == x:
            return table[mid][1], comps
        if y < x:
            lo = mid + 1
        else:
            hi = mid
    return None, comps


def build_babies(p: int, a: int, gen: Point, m: int) -> list[tuple[int, int]]:
    assert gen is not None
    two = affine_add(p, a, gen, gen)
    acc: Jac = (gen[0], gen[1], 1)
    babies: list[Jac] = []
    for _ in range(2 ** (m - 1)):
        babies.append(acc)
        acc = jac_add_affine(p, a, acc, two)
    xs = batch_x(p, babies)
    rows = [(xs[j], j) for j in range(len(xs))]
    if any(x is None for x, _ in rows):
        raise RuntimeError("baby at infinity")
    table, comps = mergesort(rows)  # type: ignore[arg-type]
    # Data-dependent merges can exhaust a run early. The proof charges the
    # worst-case recurrence, which adds len-1 at every merge.
    if comps > mergesort_bound(len(rows)):
        raise RuntimeError(f"mergesort {comps} > {mergesort_bound(len(rows))}")
    if len({x for x, _ in table}) != len(table):
        raise RuntimeError("x collision")
    return table


def giant_xs(p: int, a: int, gen: Point, q: Point, m: int) -> list[int | None]:
    assert gen is not None and q is not None
    s = affine_mul(p, a, 2**m, gen)
    assert s is not None
    stride = affine_add(p, a, s, s)
    assert stride is not None
    neg_s = (s[0], (-s[1]) % p)
    neg_stride = (stride[0], (-stride[1]) % p)
    acc = jac_add_affine(p, a, (q[0], q[1], 1), neg_s)
    pts: list[Jac] = []
    for _ in range(2**m):
        pts.append(acc)
        acc = jac_add_affine(p, a, acc, neg_stride)
    return batch_x(p, pts)


def recover(p: int, a: int, n: int, gen: Point, q: Point, m: int) -> int:
    if q is None:
        return 0
    if 2 * n > 2 ** (2 * m + 1):
        raise RuntimeError("partition does not cover 2n")
    table = build_babies(p, a, gen, m)
    xs = giant_xs(p, a, gen, q, m)
    for i, x in enumerate(xs):
        if x is None:
            continue
        j, comps = binary_search(table, x)
        if comps > (len(table) - 1).bit_length() + 1:
            raise RuntimeError("binary probe")
        if j is None:
            continue
        base = (2 * i + 1) << m
        step = 2 * j + 1
        for sigma in (1, -1):
            z = base + sigma * step
            if not 0 < z < 2 * n:
                continue
            k = z % n
            if affine_mul(p, a, k, gen) == q:
                return k
    raise RuntimeError("miss")


def point_count(p: int, a: int, b: int) -> list[Point]:
    squares = {pow(y, 2, p) for y in range(p)}
    pts: list[Point] = [None]
    for x in range(p):
        rhs = (x * x * x + a * x + b) % p
        if rhs not in squares:
            continue
        for y in range(p):
            if (y * y - rhs) % p == 0:
                pts.append((x, y))
    return pts


def find_toy() -> tuple[int, int, int, int, Point]:
    p = 17
    for a in range(p):
        for b in range(p):
            pts = point_count(p, a, b)
            n = len(pts)
            if n not in (11, 13, 17, 19, 23):
                continue
            gen = next(pt for pt in pts if pt is not None)
            if affine_mul(p, a, n, gen) is not None:
                continue
            if affine_mul(p, a, 1, gen) != gen:
                continue
            return p, a, b, n, gen
    raise RuntimeError("no prime-order toy curve")


def cios(a: int, b: int, p: int, n0: int) -> int:
    """Four-limb CIOS. Returns a*b*2^{-256} mod p.

    Limb schedule matches Koc-Acar-Kaliski: the modulus product is written one
    limb lower as it is formed, which is the division by 2^64. Two conditional
    subtractions of p cover a raw result in [0, 2p).
    """
    w = 64
    s = 4
    word = 1 << w
    t = [0] * (s + 2)
    a_limbs = [(a >> (w * j)) & (word - 1) for j in range(s)]
    b_limbs = [(b >> (w * i)) & (word - 1) for i in range(s)]
    n_limbs = [(p >> (w * j)) & (word - 1) for j in range(s)]
    for i in range(s):
        carry = 0
        for j in range(s):
            uv = t[j] + a_limbs[j] * b_limbs[i] + carry
            t[j] = uv % word
            carry = uv // word
        uv = t[s] + carry
        t[s] = uv % word
        t[s + 1] = uv // word
        m = (t[0] * n0) % word
        uv = t[0] + m * n_limbs[0]
        if uv % word != 0:
            raise RuntimeError("CIOS low limb")
        carry = uv // word
        for j in range(1, s):
            uv = t[j] + m * n_limbs[j] + carry
            t[j - 1] = uv % word
            carry = uv // word
        uv = t[s] + carry
        t[s - 1] = uv % word
        carry = uv // word
        t[s] = t[s + 1] + carry
        t[s + 1] = 0
    r = 0
    for j in range(s):
        r |= t[j] << (w * j)
    # The top limb is 0 or 1. Fold 2^256 into the integer, then at most two
    # subtractions of p. A third borrow means the schedule's bound is wrong.
    if t[s] > 1:
        raise RuntimeError("CIOS extra limb")
    if t[s] == 1:
        r += 1 << 256
    if r >= p:
        r -= p
    if r >= p:
        r -= p
    if r >= p:
        raise RuntimeError("CIOS not reduced")
    return r


def check_cios() -> None:
    p = P256_P
    n0 = (-modinv(p, 1 << 64)) % (1 << 64)
    r_inv = modinv(1 << 256, p)
    samples = [1, 2, 3, p - 1, p - 2, (1 << 255) - 19, 0x123456789ABCDEF]
    for a in samples:
        for b in samples:
            got = cios(a, b, p, n0)
            expect = (a * b * r_inv) % p
            if got != expect:
                raise RuntimeError("CIOS mismatch")


def production_coeffs() -> dict[str, int]:
    mixed_normal = 11 * M_INSTR + 16 * A_INSTR
    double = 9 * M_INSTR + 24 * A_INSTR
    # Baby adds never hit the equal-summand branch; giant adds hit a special
    # branch at most four times per input, charged in the additive 2^25 term.
    mixed = mixed_normal
    # Coefficients of 2^128 in the instruction total.
    attack = (3 * mixed) // 2
    attack += 6 * M_INSTR
    attack += (3 * M_INSTR) // 2
    attack += (382 * M_INSTR) // 4096 + 1
    attack += 128 * CMP_INSTR
    attack += (126 * CMP_INSTR) // 2
    attack += (127 * MOVE_INSTR) // 2
    attack += MOVE_INSTR // 2
    attack += 64 + 32 + 128
    pre = mixed // 2
    pre += 2 * M_INSTR
    pre += M_INSTR // 2
    pre += (126 * CMP_INSTR) // 2
    pre += (127 * MOVE_INSTR) // 2
    pre += 16 + 32 + 64 + 1
    return {"attack": attack, "preprocessing": pre, "mixed": mixed, "double": double}


def check_production_integers() -> dict[str, int]:
    e = P256_P - 2
    if e.bit_count() != 128 or e.bit_length() != 256:
        raise RuntimeError("Fermat weight")
    inv_muls = (e.bit_length() - 1) + (e.bit_count() - 1)
    if inv_muls != 382:
        raise RuntimeError("Fermat length")
    coeffs = production_coeffs()
    if not coeffs["attack"] < 32000:
        raise RuntimeError(f"attack coeff {coeffs['attack']}")
    if not coeffs["preprocessing"] < 16000:
        raise RuntimeError(f"pre coeff {coeffs['preprocessing']}")
    # T < attack * 2^128 + 2^25, and 1000 * 2^133 = 32000 * 2^128.
    if coeffs["attack"] * 2**128 + 2**25 >= 1000 * 2**133:
        raise RuntimeError("time bound")
    if coeffs["preprocessing"] * 2**128 + 2**25 >= 1000 * 2**132:
        raise RuntimeError("preprocessing bound")
    if 3 * 2**132 + 2**20 >= 2**134:
        raise RuntimeError("memory bound")
    # Block inversions: 2^116 + 1, each 382 field multiplications.
    inv_coeff_num = 382 * M_INSTR
    # inv_coeff_num * 2^116 / 2^128 = inv_coeff_num / 4096 < 63
    if inv_coeff_num // 4096 >= 63:
        raise RuntimeError("inv coeff")
    return coeffs


def main() -> None:
    spec = json.load(sys.stdin)
    width = spec["width"]
    if not isinstance(width, int) or not 3 <= width <= 4:
        raise SystemExit("width")
    for m in (3, 4, 5):
        if not unique_odds(m):
            raise SystemExit(f"partition {m}")
        n = 2 ** (m - 1)
        if mergesort_bound(n) != _merge_rec(n):
            raise SystemExit("merge closed form")
    p, a, b, n, gen = find_toy()
    for k in range(n):
        pt = affine_mul(p, a, k, gen)
        if pt is None:
            continue
        # Jacobian mixed chain of k additions of G must match affine k*G.
        acc: Jac = None
        for _ in range(k):
            acc = jac_add_affine(p, a, acc, gen)
        xs = batch_x(p, [acc])
        if xs[0] != pt[0]:
            raise SystemExit(f"jacobian x at {k}")
    if 2 * n > 2 ** (2 * width + 1):
        raise SystemExit("width too small for toy order")
    # Production babies sit below the group order. The same separation is
    # required on the toy curve, otherwise an odd multiple lands on infinity.
    if 2**width >= n:
        raise SystemExit("baby window meets the group order")
    checked = 0
    for secret in range(n):
        target = affine_mul(p, a, secret, gen)
        got = recover(p, a, n, gen, target, width)
        if got != secret:
            raise SystemExit(f"{secret} -> {got}")
        checked += 1
    # Odd representative agrees with the point.
    for secret in range(1, n):
        z = odd_representative(secret, n)
        if z % 2 == 0 or affine_mul(p, a, z % n, gen) != affine_mul(p, a, secret, gen):
            raise SystemExit("odd representative")
    check_cios()
    coeffs = check_production_integers()
    print(
        json.dumps(
            {
                "toy_field": p,
                "toy_a": a,
                "toy_b": b,
                "toy_order": n,
                "width": width,
                "checked": checked,
                "all_recovered": True,
                "cios_samples": 49,
                "attack_coeff": coeffs["attack"],
                "preprocessing_coeff": coeffs["preprocessing"],
                "scope": "toy recovery and integer bound arithmetic only",
            },
            sort_keys=True,
        )
    )


def _merge_rec(n: int) -> int:
    memo = {0: 0, 1: 0}

    def rec(length: int) -> int:
        if length in memo:
            return memo[length]
        mid = length // 2
        memo[length] = rec(mid) + rec(length - mid) + length - 1
        return memo[length]

    return rec(n)


if __name__ == "__main__":
    main()
