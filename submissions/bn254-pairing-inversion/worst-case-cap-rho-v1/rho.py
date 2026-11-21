"""Executable mathematical reference. Python runtime is NOT the scored machine.

The priced scalar implementation and its finite instruction bounds are in
machine.md and arithmetic.c. Test coins are reproducible, not a success proof.
"""
import hashlib
import math


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def rho(r, target, identity, mul, power, encode, coins):
    """One bounded 32-adding walk, Floyd collision detection, no hidden scalar."""
    if target == identity:
        return 0, 0, "identity"
    salt = coins.getrandbits(256).to_bytes(32, "big")
    table = []
    for _ in range(32):
        a, b = sample_scalar(r, coins), sample_scalar(r, coins)
        table.append((mul(power(a), power(b, target)), a, b))
    a, b = sample_scalar(r, coins), sample_scalar(r, coins)
    initial = (mul(power(a), power(b, target)), a, b)

    def step(state):
        y, a, b = state
        j = hashlib.sha256(salt + encode(y)).digest()[0] & 31
        m, c, d = table[j]
        return mul(y, m), (a + c) % r, (b + d) % r

    slow = fast = initial
    cap = math.isqrt(121 * r) // 10 + 1
    while 100 * cap * cap < 121 * r:
        cap += 1
    for k in range(1, cap + 1):
        slow = step(slow)
        fast = step(step(fast))
        if slow[0] == fast[0]:
            denominator = (fast[2] - slow[2]) % r
            if denominator == 0:
                return None, k, "degenerate"
            x = ((slow[1] - fast[1]) * pow(denominator, r - 2, r)) % r
            if power(x) != target:
                return None, k, "verification_failed"
            return x, k, "success"
    return None, cap, "cap"


def fp12_mul(a, b, p):
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


def montgomery_word_reference(a, b, modulus):
    """Exact fixed-loop word algorithm matching arithmetic.c, including carries."""
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
