"""Executable mathematical reference. Python runtime is NOT the scored machine.

The priced scalar implementation and its finite instruction bounds are in
machine.md and arithmetic.c. Test coins are reproducible, not a success proof.

Attribution: the 32-adding walk, the salted SHA-256 partition, the bounded
rejection sampling, the Montgomery word reference, and the polynomial-basis
Fp12 routines follow the track's rho-baseline-v2-20261110 package. The
distinguished-point detection, the tower representation, and the Karatsuba
formulas are new to this package; the walk itself is not claimed as original.
"""
import hashlib
import math


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def rho(r, target, identity, mul, power, encode, coins, dp_bits=110,
        record_cap=1 << 20):
    """One bounded 32-adding walk, distinguished-point detection, no restart.

    The walk, table construction, partition, and cap are the baseline's. The
    only detection change is Floyd tortoise-and-hare replaced by a
    distinguished-point table: a state is distinguished when the low dp_bits
    of the big-endian integer formed by SHA-256 digest bytes 1..14 are zero.
    The digest is the same one already computed for the partition index.
    Each visited state is tested before stepping, so exactly one digest is
    computed per step. A repeated distinguished point triggers recovery; a
    degenerate coefficient difference, the step cap, or a full record table
    reports failure within the cap. No hidden scalar is used.
    """
    if target == identity:
        return 0, 0, "identity"
    salt = coins.getrandbits(256).to_bytes(32, "big")
    table = []
    for _ in range(32):
        a, b = sample_scalar(r, coins), sample_scalar(r, coins)
        table.append((mul(power(a), power(b, target)), a, b))
    a, b = sample_scalar(r, coins), sample_scalar(r, coins)
    state = (mul(power(a), power(b, target)), a, b)
    stored = {}
    cap = 2 * (math.isqrt(r) + (math.isqrt(r) ** 2 < r))
    mask = (1 << dp_bits) - 1
    for step in range(1, cap + 1):
        y, a, b = state
        digest = hashlib.sha256(salt + encode(y)).digest()
        j = digest[0] & 31
        if int.from_bytes(digest[1:15], "big") & mask == 0:
            if y in stored:
                a1, b1 = stored[y]
                denominator = (b - b1) % r
                if denominator == 0:
                    return None, step, "degenerate"
                x = ((a1 - a) * pow(denominator, r - 2, r)) % r
                if power(x) != target:
                    return None, step, "verification_failed"
                return x, step, "success"
            if len(stored) >= record_cap:
                return None, step, "table_full"
            stored[y] = (a, b)
        m, c, d = table[j]
        state = (mul(y, m), (a + c) % r, (b + d) % r)
    return None, cap, "cap"


def fp12_mul(a, b, p):
    """Polynomial-basis Fp12 multiply (baseline formula), for cross-checks."""
    c = [0] * 23
    for i in range(12):
        for j in range(12):
            c[i + j] = (c[i + j] + a[i] * b[j]) % p
    for k in range(22, 11, -1):
        c[k - 6] = (c[k - 6] + 18 * c[k]) % p
        c[k - 12] = (c[k - 12] - 82 * c[k]) % p
    return tuple(c[:12])


def fp12_pow(a, n, p):
    out = (1,) + (0,) * 11
    while n:
        if n & 1:
            out = fp12_mul(out, a, p)
        a = fp12_mul(a, a, p)
        n >>= 1
    return out


# ---- Tower representation ----
# Fp2 = Fp[u]/(u^2+1); Fp6 = Fp2[v]/(v^3-xi), xi = 9+u; Fp12 = Fp6[t]/(t^2-v).
# Fp2 is (c0, c1); Fp6 is (b0, b1, b2) of Fp2s; Fp12 is (a0, a1) of Fp6s.
# The embedding into the protected polynomial basis Fp[w]/(w^12-18w^6+82) is
# u = w^6-9, v = w^2, t = w. All tower operations below use the Karatsuba
# formulas priced in machine.md. The optional count dict records the exact
# numbers of Fp multiplications and Fp add/sub operations, so experiment.py
# can check the priced operation counts instead of asserting them.

def _fp_add(a, b, p, count):
    if count is not None:
        count["adds"] += 1
    return (a + b) % p


def _fp_sub(a, b, p, count):
    if count is not None:
        count["adds"] += 1
    return (a - b) % p


def _fp_mul(a, b, p, count):
    if count is not None:
        count["muls"] += 1
    return (a * b) % p


def fp2_add(a, b, p, count=None):
    return (_fp_add(a[0], b[0], p, count), _fp_add(a[1], b[1], p, count))


def fp2_sub(a, b, p, count=None):
    return (_fp_sub(a[0], b[0], p, count), _fp_sub(a[1], b[1], p, count))


def fp2_mul(a, b, p, count=None):
    """Complex Karatsuba: 3 Fp muls + 5 Fp add/subs."""
    ac = _fp_mul(a[0], b[0], p, count)
    bd = _fp_mul(a[1], b[1], p, count)
    s = _fp_mul(_fp_add(a[0], a[1], p, count),
                _fp_add(b[0], b[1], p, count), p, count)
    return (_fp_sub(ac, bd, p, count),
            _fp_sub(_fp_sub(s, ac, p, count), bd, p, count))


def fp2_mul_xi(a, p, count=None):
    """Multiply by xi = 9+u: (9a0-a1) + (a0+9a1)u. 9x = 8x+x via 3 doubles."""
    def mul9(x):
        d = _fp_add(x, x, p, count)
        d = _fp_add(d, d, p, count)
        d = _fp_add(d, d, p, count)
        return _fp_add(d, x, p, count)
    return (_fp_sub(mul9(a[0]), a[1], p, count),
            _fp_add(a[0], mul9(a[1]), p, count))


def fp6_add(a, b, p, count=None):
    return tuple(fp2_add(a[i], b[i], p, count) for i in range(3))


def fp6_sub(a, b, p, count=None):
    return tuple(fp2_sub(a[i], b[i], p, count) for i in range(3))


def fp6_mul(a, b, p, count=None):
    """6 Fp2 muls. v^3 = xi; c0 = v0+xi*v5, c1 = v3+xi*v2, c2 = v4+v1."""
    v0 = fp2_mul(a[0], b[0], p, count)
    v1 = fp2_mul(a[1], b[1], p, count)
    v2 = fp2_mul(a[2], b[2], p, count)
    v3 = fp2_sub(fp2_sub(fp2_mul(fp2_add(a[0], a[1], p, count),
                                 fp2_add(b[0], b[1], p, count), p, count),
                         v0, p, count), v1, p, count)
    v4 = fp2_sub(fp2_sub(fp2_mul(fp2_add(a[0], a[2], p, count),
                                 fp2_add(b[0], b[2], p, count), p, count),
                         v0, p, count), v2, p, count)
    v5 = fp2_sub(fp2_sub(fp2_mul(fp2_add(a[1], a[2], p, count),
                                 fp2_add(b[1], b[2], p, count), p, count),
                         v1, p, count), v2, p, count)
    return (fp2_add(v0, fp2_mul_xi(v5, p, count), p, count),
            fp2_add(v3, fp2_mul_xi(v2, p, count), p, count),
            fp2_add(v4, v1, p, count))


def fp6_mul_v(a, p, count=None):
    """Multiply by v: (xi*a2) + a0*v + a1*v^2."""
    return (fp2_mul_xi(a[2], p, count), a[0], a[1])


def fp12_tower_mul(a, b, p, count=None):
    """Karatsuba over Fp6: 3 Fp6 muls, t^2 = v."""
    v0 = fp6_mul(a[0], b[0], p, count)
    v1 = fp6_mul(a[1], b[1], p, count)
    v2 = fp6_sub(fp6_sub(fp6_mul(fp6_add(a[0], a[1], p, count),
                                 fp6_add(b[0], b[1], p, count), p, count),
                         v0, p, count), v1, p, count)
    return (fp6_add(v0, fp6_mul_v(v1, p, count), p, count), v2)


def fp12_tower_pow(a, n, p):
    one = (((1, 0), (0, 0), (0, 0)), ((0, 0), (0, 0), (0, 0)))
    out = one
    while n:
        if n & 1:
            out = fp12_tower_mul(out, a, p)
        a = fp12_tower_mul(a, a, p)
        n >>= 1
    return out


def to_tower(c, p):
    """Protected polynomial coefficients [c0..c11] -> tower form.

    With u = w^6-9, v = w^2, t = w: the even part is
    (c0+9c6) + c6*u + ((c2+9c8) + c8*u)*v + ((c4+9c10) + c10*u)*v^2,
    and the odd part is the same pattern on c1, c3, c5, c7, c9, c11 times t.
    """
    def f2(x0, x1):
        return (x0 % p, x1 % p)
    a0 = (f2(c[0] + 9 * c[6], c[6]), f2(c[2] + 9 * c[8], c[8]),
          f2(c[4] + 9 * c[10], c[10]))
    a1 = (f2(c[1] + 9 * c[7], c[7]), f2(c[3] + 9 * c[9], c[9]),
          f2(c[5] + 9 * c[11], c[11]))
    return (a0, a1)


def to_poly(t, p):
    """Tower form -> protected polynomial coefficients [c0..c11]."""
    a0, a1 = t
    c = [0] * 12
    c[0] = (a0[0][0] - 9 * a0[0][1]) % p
    c[6] = a0[0][1] % p
    c[2] = (a0[1][0] - 9 * a0[1][1]) % p
    c[8] = a0[1][1] % p
    c[4] = (a0[2][0] - 9 * a0[2][1]) % p
    c[10] = a0[2][1] % p
    c[1] = (a1[0][0] - 9 * a1[0][1]) % p
    c[7] = a1[0][1] % p
    c[3] = (a1[1][0] - 9 * a1[1][1]) % p
    c[9] = a1[1][1] % p
    c[5] = (a1[2][0] - 9 * a1[2][1]) % p
    c[11] = a1[2][1] % p
    return tuple(c)


def montgomery_word_reference(a, b, modulus):
    """Exact fixed-loop word algorithm matching arithmetic.c, including carries.

    Reused from the rho-baseline-v2-20261110 package; the priced bound in
    machine.md is a tighter cell count over this same loop structure.
    """
    mask = (1 << 64) - 1
    aa = [(a >> (64 * i)) & mask for i in range(4)]
    bb = [(b >> (64 * i)) & mask for i in range(4)]
    mm = [(modulus >> (64 * i)) & mask for i in range(4)]
    n0 = (-pow(mm[0], -1, 1 << 64)) & mask
    t = [0] * 9
    for i in range(4):
        carry = 0
        for j in range(4):
            v = aa[i] * bb[j] + t[i + j] + carry
            assert v < 1 << 128
            t[i + j], carry = v & mask, v >> 64
        for k in range(i + 4, 9):
            v = t[k] + carry
            t[k], carry = v & mask, v >> 64
        assert carry == 0
    for i in range(4):
        q = (t[i] * n0) & mask
        carry = 0
        for j in range(4):
            v = q * mm[j] + t[i + j] + carry
            assert v < 1 << 128
            t[i + j], carry = v & mask, v >> 64
        for k in range(i + 4, 9):
            v = t[k] + carry
            t[k], carry = v & mask, v >> 64
        assert carry == 0 and t[i] == 0
    u = sum(t[i + 4] << (64 * i) for i in range(5))
    assert 0 <= u < 2 * modulus < 1 << 255
    return u - modulus if u >= modulus else u
