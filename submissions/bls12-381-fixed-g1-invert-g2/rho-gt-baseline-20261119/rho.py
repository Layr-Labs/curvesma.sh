"""Executable mathematical reference. Python runtime is NOT the scored machine.

The priced scalar implementation and its finite instruction bounds are in
machine.md and arithmetic.c. Test coins are reproducible, not a success proof.

In the attack, the base is g=e(P,Q0) (protected py_ecc pairing, computed
once), the target is z, and the output is Q=[x]Q0 for the recovered x.
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
    """One bounded 32-adding walk, Floyd collision detection, no hidden scalar.

    power(a) = g^a, power(b, target) = target^b. Returns (x, k, outcome).
    """
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
    cap = 2 * (math.isqrt(r) + (math.isqrt(r) ** 2 < r))
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
    """Protected basis Fp[w]/(w^12 - 2 w^6 + 2): w^12 = 2 w^6 - 2."""
    c = [0] * 23
    for i in range(12):
        for j in range(12):
            c[i + j] = (c[i + j] + a[i] * b[j]) % p
    for k in range(22, 11, -1):
        t = (2 * c[k]) % p
        c[k - 6] = (c[k - 6] + t) % p
        c[k - 12] = (c[k - 12] - t) % p
    return tuple(c[:12])


def fp12_pow(a, n, p):
    out = (1,) + (0,) * 11
    while n:
        if n & 1:
            out = fp12_mul(out, a, p)
        a = fp12_mul(a, a, p)
        n >>= 1
    return out


def montgomery_word_reference(a, b, modulus, limbs):
    """Exact fixed-loop word algorithm matching arithmetic.c, including carries."""
    mask = (1 << 64) - 1
    n = limbs
    aa = [(a >> (64 * i)) & mask for i in range(n)]
    bb = [(b >> (64 * i)) & mask for i in range(n)]
    mm = [(modulus >> (64 * i)) & mask for i in range(n)]
    n0 = (-pow(mm[0], -1, 1 << 64)) & mask
    t = [0] * (2 * n + 1)
    for i in range(n):
        carry = 0
        for j in range(n):
            v = aa[i] * bb[j] + t[i + j] + carry
            assert v < 1 << 128
            t[i + j], carry = v & mask, v >> 64
        for k in range(i + n, 2 * n + 1):
            v = t[k] + carry
            t[k], carry = v & mask, v >> 64
        assert carry == 0
    for i in range(n):
        q = (t[i] * n0) & mask
        carry = 0
        for j in range(n):
            v = q * mm[j] + t[i + j] + carry
            assert v < 1 << 128
            t[i + j], carry = v & mask, v >> 64
        for k in range(i + n, 2 * n + 1):
            v = t[k] + carry
            t[k], carry = v & mask, v >> 64
        assert carry == 0 and t[i] == 0
    u = sum(t[i + n] << (64 * i) for i in range(n + 1))
    assert 0 <= u < 2 * modulus < 1 << (64 * n)
    return u - modulus if u >= modulus else u
