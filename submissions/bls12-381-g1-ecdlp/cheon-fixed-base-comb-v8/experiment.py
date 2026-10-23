from __future__ import annotations
import os
import sys
import hashlib
from collections import Counter
MASK = (1 << 64) - 1
_ROT = (
    (0, 36, 3, 41, 18),
    (1, 44, 10, 45, 2),
    (62, 6, 43, 15, 61),
    (28, 55, 25, 21, 56),
    (27, 20, 39, 8, 14),
)
_RC = (
    0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000,
    0x000000000000808B, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
    0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
    0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003,
    0x8000000000008002, 0x8000000000000080, 0x000000000000800A, 0x800000008000000A,
    0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
)
class Cost:
    def __init__(self) -> None:
        self.ops: Counter[str] = Counter()
    def op(self, name: str, k: int = 1) -> None:
        self.ops[name] += k
    @property
    def total(self) -> int:
        return sum(self.ops.values())
    def snapshot(self) -> int:
        return self.total
    def copy(self) -> dict[str, int]:
        return dict(self.ops)
def _rotl(x: int, n: int, cost: Cost) -> int:
    cost.op("ror")
    n &= 63
    if n == 0:
        return x & MASK
    return ((x << n) | (x >> (64 - n))) & MASK
def keccak_f(state: list[int], cost: Cost) -> None:
    a = state
    for rc in _RC:
        c = [0] * 5
        for x in range(5):
            v = a[x]
            for y in range(1, 5):
                cost.op("eor")
                v ^= a[x + 5 * y]
            c[x] = v
        d = [0] * 5
        for x in range(5):
            d[x] = c[(x - 1) % 5] ^ _rotl(c[(x + 1) % 5], 1, cost)
            cost.op("eor")
        for x in range(5):
            for y in range(5):
                cost.op("eor")
                a[x + 5 * y] ^= d[x]
        b = [0] * 25
        for x in range(5):
            for y in range(5):
                cost.op("str")
                rot = _ROT[x][y]
                val = a[x + 5 * y]
                if rot:
                    val = _rotl(val, rot, cost)
                nx, ny = y, (2 * x + 3 * y) % 5
                b[nx + 5 * ny] = val
                cost.op("ldr")
        for y in range(5):
            row = [b[x + 5 * y] for x in range(5)]
            for x in range(5):
                cost.op("ldr")
            for x in range(5):
                cost.op("mvn")
                cost.op("and")
                cost.op("eor")
                cost.op("str")
                a[x + 5 * y] = row[x] ^ ((~row[(x + 1) % 5]) & row[(x + 2) % 5] & MASK)
        cost.op("ldr")
        cost.op("eor")
        a[0] ^= rc
def shake256(data: bytes, out_len: int, cost: Cost | None = None) -> bytes:
    local = cost is None
    if cost is None:
        cost = Cost()
    rate = 136
    st = [0] * 25
    msg = bytearray(data)
    msg.append(0x1F)
    offset = 0
    while offset + rate <= len(msg):
        block = msg[offset : offset + rate]
        _xor_block(st, block, cost)
        keccak_f(st, cost)
        offset += rate
    tail = bytearray(rate)
    tail[: len(msg) - offset] = msg[offset:]
    tail[rate - 1] ^= 0x80
    _xor_block(st, tail, cost)
    keccak_f(st, cost)
    out = bytearray()
    while len(out) < out_len:
        out.extend(_store_block(st, rate, cost))
        if len(out) < out_len:
            keccak_f(st, cost)
    if local:
        shake256.last_cost = cost.total  # type: ignore[attr-defined]
    return bytes(out[:out_len])
def _xor_block(st: list[int], block: bytes, cost: Cost) -> None:
    for i in range(0, len(block), 8):
        chunk = block[i : i + 8]
        word = int.from_bytes(chunk.ljust(8, b"\x00"), "little")
        cost.op("ldr")
        cost.op("eor")
        cost.op("str")
        lane = i // 8
        st[lane] ^= word
def _store_block(st: list[int], rate: int, cost: Cost) -> bytes:
    raw = bytearray()
    for i in range((rate + 7) // 8):
        cost.op("ldr")
        cost.op("str")
        raw.extend(st[i].to_bytes(8, "little"))
    return bytes(raw[:rate])
def _mac(lo_in: int, a: int, b: int, carry: int, cost: Cost) -> tuple[int, int]:
    cost.op("mul")
    cost.op("umulh")
    cost.op("adds")
    cost.op("adc")
    cost.op("adds")
    cost.op("adc")
    total = lo_in + a * b + carry
    return total >> 64, total & MASK
class Field:
    def __init__(self, prime: int, limbs: int, cost: Cost | None = None) -> None:
        if prime.bit_length() > 64 * limbs:
            raise ValueError("prime does not fit")
        self.p = prime
        self.n = limbs
        self.cost = cost if cost is not None else Cost()
        self.mod = [(prime >> (64 * i)) & MASK for i in range(limbs)]
        self.n0inv = _neg_inv64(self.mod[0], self.cost)
        r = (1 << (64 * limbs)) % prime
        self.r2 = self._limbs((r * r) % prime)
        self.one = self._limbs(r)
        self.zero = [0] * limbs
        self._setup_cost = self.cost.total
    def _limbs(self, x: int) -> list[int]:
        return [(x >> (64 * i)) & MASK for i in range(self.n)]
    def _int(self, limbs: list[int]) -> int:
        v = 0
        for i, w in enumerate(limbs[: self.n]):
            v |= w << (64 * i)
        return v
    def to_mont(self, x: int) -> list[int]:
        return self._cios(self._limbs(x % self.p), self.r2)
    def from_mont(self, a: list[int]) -> int:
        return self._int(self._cios(a, self._limbs(1))) % self.p
    def mul(self, a: list[int], b: list[int]) -> list[int]:
        return self._cios(a, b)
    def sqr(self, a: list[int]) -> list[int]:
        return self._cios(a, a)
    def add(self, a: list[int], b: list[int]) -> list[int]:
        return self._add(a, b)
    def sub(self, a: list[int], b: list[int]) -> list[int]:
        return self._sub(a, b)
    def mul2(self, a: list[int]) -> list[int]:
        return self._add(a, a)
    def mul3(self, a: list[int]) -> list[int]:
        return self._add(self._add(a, a), a)
    def mul4(self, a: list[int]) -> list[int]:
        return self.mul2(self.mul2(a))
    def mul8(self, a: list[int]) -> list[int]:
        return self.mul2(self.mul4(a))
    def inv(self, a: list[int]) -> list[int]:
        e = self.p - 2
        result = self.to_mont(1)
        base = a
        bits = e.bit_length()
        for i in range(bits - 1, -1, -1):
            self.cost.op("lsr")
            self.cost.op("tst")
            self.cost.op("b")
            self.cost.op("add")
            result = self.sqr(result)
            if (e >> i) & 1:
                result = self.mul(result, base)
        return result
    def _cios(self, a: list[int], b: list[int]) -> list[int]:
        n = self.n
        cost = self.cost
        t = [0] * (n + 2)
        cost.op("ldr")  # n0inv
        for i in range(n):
            cost.op("ldr")  # b[i], then held
            carry = 0
            for j in range(n):
                cost.op("ldr")  # a[j]
                cost.op("ldr")  # t[j]
                carry, t[j] = _mac(t[j], a[j], b[i], carry, cost)
                cost.op("str")
            cost.op("ldr")
            cost.op("adds")
            cost.op("adc")
            cost.op("str")
            cost.op("str")
            wide = t[n] + carry
            t[n] = wide & MASK
            t[n + 1] = wide >> 64
            cost.op("ldr")  # t[0]
            cost.op("mul")
            m = (t[0] * self.n0inv) & MASK
            carry = 0
            cost.op("ldr")  # modulus limb 0
            carry, _low = _mac(t[0], m, self.mod[0], carry, cost)
            for j in range(1, n):
                cost.op("ldr")  # t[j]
                cost.op("ldr")  # modulus[j]
                carry, t[j - 1] = _mac(t[j], m, self.mod[j], carry, cost)
                cost.op("str")
            cost.op("ldr")
            cost.op("adds")
            cost.op("adc")
            cost.op("str")
            wide = t[n] + carry
            t[n - 1] = wide & MASK
            cost.op("adds")
            cost.op("str")
            t[n] = t[n + 1] + (wide >> 64)
        return self._cond_sub(t[:n], t[n])
    def _cond_sub(self, t: list[int], extra: int) -> list[int]:
        cost = self.cost
        n = self.n
        diff = [0] * n
        borrow = 0
        for i in range(n):
            cost.op("ldr")
            cost.op("ldr")
            cost.op("subs")
            cost.op("sbc")
            cost.op("str")
            wide = t[i] - self.mod[i] - borrow
            diff[i] = wide & MASK
            borrow = 1 if wide < 0 else 0
        cost.op("cmp")
        cost.op("csel")
        take = extra != 0 or borrow == 0
        out = [0] * n
        for i in range(n):
            cost.op("ldr")
            cost.op("ldr")
            cost.op("csel")
            cost.op("str")
            out[i] = diff[i] if take else t[i]
        return out
    def _add(self, a: list[int], b: list[int]) -> list[int]:
        cost = self.cost
        n = self.n
        t = [0] * n
        carry = 0
        for i in range(n):
            cost.op("ldr")
            cost.op("ldr")
            cost.op("adds")
            cost.op("adc")
            cost.op("str")
            wide = a[i] + b[i] + carry
            t[i] = wide & MASK
            carry = wide >> 64
        return self._cond_sub(t, carry)
    def _sub(self, a: list[int], b: list[int]) -> list[int]:
        cost = self.cost
        n = self.n
        t = [0] * n
        borrow = 0
        for i in range(n):
            cost.op("ldr")
            cost.op("ldr")
            cost.op("subs")
            cost.op("sbc")
            cost.op("str")
            wide = a[i] - b[i] - borrow
            t[i] = wide & MASK
            borrow = 1 if wide < 0 else 0
        added = [0] * n
        carry = 0
        for i in range(n):
            cost.op("ldr")
            cost.op("ldr")
            cost.op("adds")
            cost.op("adc")
            cost.op("str")
            wide = t[i] + self.mod[i] + carry
            added[i] = wide & MASK
            carry = wide >> 64
        cost.op("cmp")
        out = [0] * n
        for i in range(n):
            cost.op("ldr")
            cost.op("ldr")
            cost.op("csel")
            cost.op("str")
            out[i] = added[i] if borrow else t[i]
        return out
def _neg_inv64(a0: int, cost: Cost) -> int:
    cost.op("mov")
    x = 1
    two = 2
    for _ in range(6):
        cost.op("mul")
        cost.op("mul")
        cost.op("sub")
        x = (x * (two - (a0 * x & MASK))) & MASK
    cost.op("neg")
    return (-x) & MASK
def reduce_hash(digest: bytes, modulus: int, cost: Cost) -> int:
    words = [int.from_bytes(digest[i : i + 8], "little") for i in range(0, 64, 8)]
    acc = 0
    bits = 512
    limb_n = (modulus.bit_length() + 63) // 64
    for i in range(bits):
        byte_i = (bits - 1 - i) // 64
        bit = (bits - 1 - i) % 64
        bit_val = (words[byte_i] >> bit) & 1
        cost.op("lsl")
        cost.op("adc")
        cost.op("adds", limb_n)  # shift the accumulator limbs
        acc = ((acc << 1) | bit_val)
        cost.op("cmp", limb_n)
        cost.op("subs", limb_n)
        cost.op("csel", limb_n)
        if acc >= modulus:
            acc -= modulus
        cost.op("add")
        cost.op("cmp")
        cost.op("b")
    return acc
INF = None
def _elem_zero(field: Field, limbs) -> bool:
    field.cost.op("ldr", field.n)
    field.cost.op("cmp", field.n)
    field.cost.op("b")
    return all(w == 0 for w in limbs)
def _z_is_zero(field: Field, p) -> bool:
    if p is INF:
        field.cost.op("ldr", field.n)
        field.cost.op("cmp", field.n)
        field.cost.op("b")
        return True
    return _elem_zero(field, p[2])
def _flag_inf(field: Field, point) -> bool:
    field.cost.op("cmp")
    field.cost.op("b")
    return point is INF
def _write_inf(field: Field):
    field.cost.op("str", field.n)
    return INF
def dbl_2009_l(field: Field, p):
    if _z_is_zero(field, p):
        return INF
    x1, y1, z1 = p
    a = field.sqr(x1)
    b = field.sqr(y1)
    c = field.sqr(b)
    d = field.mul2(field.sub(field.sub(field.sqr(field.add(x1, b)), a), c))
    e = field.mul3(a)
    f = field.sqr(e)
    x3 = field.sub(f, field.mul2(d))
    y3 = field.sub(field.mul(e, field.sub(d, x3)), field.mul8(c))
    z3 = field.mul(field.mul2(y1), z1)
    return x3, y3, z3
def dbl_2001_b(field: Field, p):
    if _z_is_zero(field, p):
        return INF
    x1, y1, z1 = p
    delta = field.sqr(z1)
    gamma = field.sqr(y1)
    beta = field.mul(x1, gamma)
    alpha = field.mul3(field.mul(field.sub(x1, delta), field.add(x1, delta)))
    x3 = field.sub(field.sqr(alpha), field.mul8(beta))
    z3 = field.sub(field.sub(field.sqr(field.add(y1, z1)), gamma), delta)
    y3 = field.sub(field.mul(alpha, field.sub(field.mul4(beta), x3)), field.mul8(field.sqr(gamma)))
    return x3, y3, z3
def _copy_affine(field: Field, q_affine):
    field.cost.op("ldr", field.n)
    field.cost.op("str", field.n)
    field.cost.op("ldr", field.n)
    field.cost.op("str", field.n)
    field.cost.op("ldr", field.n)
    field.cost.op("str", field.n)
    x2, y2 = q_affine
    return x2, y2, list(field.one)
def add_2007_bl(field: Field, p, q, dbl=None):
    if dbl is None:
        dbl = dbl_2009_l
    if _z_is_zero(field, p):
        return q
    if _z_is_zero(field, q):
        return p
    x1, y1, z1 = p
    x2, y2, z2 = q
    z1z1 = field.sqr(z1)
    z2z2 = field.sqr(z2)
    u1 = field.mul(x1, z2z2)
    u2 = field.mul(x2, z1z1)
    s1 = field.mul(field.mul(y1, z2), z2z2)
    s2 = field.mul(field.mul(y2, z1), z1z1)
    h = field.sub(u2, u1)
    r = field.mul2(field.sub(s2, s1))
    h_zero = _elem_zero(field, h)
    r_zero = _elem_zero(field, r)
    if h_zero and r_zero:
        return dbl(field, p)
    if h_zero:
        return _write_inf(field)
    i = field.sqr(field.mul2(h))
    j = field.mul(h, i)
    v = field.mul(u1, i)
    x3 = field.sub(field.sub(field.sqr(r), j), field.mul2(v))
    y3 = field.sub(field.mul(r, field.sub(v, x3)), field.mul(field.mul2(s1), j))
    z3 = field.mul(field.sub(field.sub(field.sqr(field.add(z1, z2)), z1z1), z2z2), h)
    return x3, y3, z3
def madd_2007_bl(field: Field, p, q_affine, dbl=None):
    if dbl is None:
        dbl = dbl_2009_l
    if _z_is_zero(field, p):
        if _flag_inf(field, q_affine):
            return INF
        return _copy_affine(field, q_affine)
    if _flag_inf(field, q_affine):
        return p
    x1, y1, z1 = p
    x2, y2 = q_affine
    z1z1 = field.sqr(z1)
    u2 = field.mul(x2, z1z1)
    s2 = field.mul(field.mul(y2, z1), z1z1)
    h = field.sub(u2, x1)
    r = field.mul2(field.sub(s2, y1))
    h_zero = _elem_zero(field, h)
    r_zero = _elem_zero(field, r)
    if h_zero and r_zero:
        return dbl(field, p)
    if h_zero:
        return _write_inf(field)
    hh = field.sqr(h)
    i = field.mul4(hh)
    j = field.mul(h, i)
    v = field.mul(x1, i)
    x3 = field.sub(field.sub(field.sqr(r), j), field.mul2(v))
    y3 = field.sub(field.mul(r, field.sub(v, x3)), field.mul(field.mul2(y1), j))
    z3 = field.sub(field.sub(field.sqr(field.add(z1, h)), z1z1), hh)
    return x3, y3, z3
def affine_add_mont(field: Field, p, q):
    x1, y1 = p
    x2, y2 = q
    field.cost.op("ldr", field.n)
    field.cost.op("cmp", field.n)
    field.cost.op("b")
    lam = field.mul(field.sub(y2, y1), field.inv(field.sub(x2, x1)))
    x3 = field.sub(field.sub(field.sqr(lam), x1), x2)
    y3 = field.sub(field.mul(lam, field.sub(x1, x3)), y1)
    return x3, y3
def affine_dbl_mont(field: Field, p, a_mont):
    x1, y1 = p
    field.cost.op("ldr", field.n)
    field.cost.op("cmp", field.n)
    field.cost.op("b")
    num = field.add(field.mul3(field.sqr(x1)), a_mont)
    den = field.mul2(y1)
    lam = field.mul(num, field.inv(den))
    x3 = field.sub(field.sub(field.sqr(lam), x1), x1)
    y3 = field.sub(field.mul(lam, field.sub(x1, x3)), y1)
    return x3, y3
def affine_step(field: Field, p, q, a_mont):
    p_inf = _flag_inf(field, p)
    q_inf = _flag_inf(field, q)
    if p_inf:
        return q
    if q_inf:
        return p
    x1, y1 = p
    x2, y2 = q
    field.cost.op("ldr", field.n)
    field.cost.op("cmp", field.n)
    field.cost.op("b")
    if x1 == x2:
        summed = field.add(y1, y2)
        if _elem_zero(field, summed):
            return INF
        num = field.add(field.mul3(field.sqr(x1)), a_mont)
        den = field.mul2(y1)
        lam = field.mul(num, field.inv(den))
        x3 = field.sub(field.sub(field.sqr(lam), x1), x1)
        y3 = field.sub(field.mul(lam, field.sub(x1, x3)), y1)
        return x3, y3
    lam = field.mul(field.sub(y2, y1), field.inv(field.sub(x2, x1)))
    x3 = field.sub(field.sub(field.sqr(lam), x1), x2)
    y3 = field.sub(field.mul(lam, field.sub(x1, x3)), y1)
    return x3, y3
def scale_z(field: Field, p):
    if p is INF:
        return None
    x1, y1, z1 = p
    a = field.inv(z1)
    aa = field.sqr(a)
    x3 = field.mul(x1, aa)
    y3 = field.mul(field.mul(y1, aa), a)
    return field.from_mont(x3), field.from_mont(y3)
def jacobian_mul(field: Field, k: int, point, dbl, add) -> object:
    result = INF
    base = point
    for bit in range(max(k.bit_length(), 1) - 1, -1, -1):
        if result is not INF:
            result = dbl(field, result)
        if (k >> bit) & 1:
            result = base if result is INF else add(field, result, base, dbl)
    if k == 0:
        return INF
    return result
def measure_call(field: Field, fn) -> int:
    start = field.cost.total
    fn()
    return field.cost.total - start
import hashlib
import json
import math
import time
from math import isqrt
BLS_X = 0xD201000000010000
BLS_P = int(
    "1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab",
    16,
)
BLS_R = BLS_X**4 - BLS_X**2 + 1
BLS_GX = int(
    "17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb",
    16,
)
BLS_GY = int(
    "08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1",
    16,
)
SECP_P = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
SECP_R = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
SECP_GX = 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
SECP_GY = 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8
P256_P = 0xFFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF
P256_R = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551
P256_A = P256_P - 3
P256_GX = 0x6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296
P256_GY = 0x4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5
C_PRICE = {
    "secp256k1-ecdlp": 1000,
    "p256-ecdlp": 1000,
    "bls12-381-g1-ecdlp": 2000,
    "bls12-381-g1-ecdlp-plain": 2000,
    "bls12-381-fixed-g1-invert-g2": 20000,
    "bls12-381-fixed-g2-invert-g1": 20000,
}
PROFILE_SHA = {
    "secp256k1-ecdlp": "18c3338de676f1909715814325bffdb1513f1621d678b942e8a7a9563b8fd086",
    "p256-ecdlp": "802d2350d1804b19fcf0e1f6d80d816862e966625710abfa2b4244ca764aa605",
    "bls12-381-g1-ecdlp": "7efc64c4327d983a9dcd2e26eee449aba20039000d368a451b7e58f61192b602",
    "bls12-381-g1-ecdlp-plain": "5094e8348d457a4f9eaf7f08564a7a9f5aa5b9860cc46cce877296a886629514",
    "bls12-381-fixed-g1-invert-g2": "af94a59fd63a014f6dd199b84e54e5a8f8ab3e5b220d432610ba785d0abef888",
    "bls12-381-fixed-g2-invert-g1": "c18402d08303b68d536a4c1ea2241420f04455e7e90b970cf9310d12374788a9",
}
COST_SHA = "e0edc25a3e8b4a4b116442b2c2eb37f81f4aa0068fb6a1906312b040d5e17294"
BASE_COMMIT = "6423c2728dc03d5f2d647c7f22adf1dac0643a66"
FLOYD_MULT = 4
RECORD_BYTES = 128
ROOT_CANDIDATES = 64
RADIX_PASSES = 6 * 8
RADIX_BUCKETS = 256
WORD_BITS = 64
def limb_words(n: int) -> int:
    if n < 0:
        raise SystemExit("negative width")
    return max(1, (max(n, 1).bit_length() + WORD_BITS - 1) // WORD_BITS)
def _winc(cost: Cost, words: int) -> None:
    for i in range(words):
        cost.op("ldr")
        cost.op("adds" if i == 0 else "adcs")
        cost.op("str")
def _wadd(cost: Cost, words: int) -> None:
    for i in range(words):
        cost.op("ldr")
        cost.op("ldr")
        cost.op("adds" if i == 0 else "adcs")
        cost.op("str")
def _wcmp(cost: Cost, words: int) -> None:
    for i in range(words):
        cost.op("ldr")
        cost.op("ldr")
        cost.op("cmp" if i == 0 else "sbcs")
    cost.op("b")
def _wshift7(cost: Cost, src_words: int, dst_words: int) -> None:
    if dst_words < src_words:
        raise SystemExit("address width shorter than the index")
    for i in range(dst_words):
        if i < src_words:
            cost.op("ldr")
        cost.op("lsl")
        cost.op("lsr")
        cost.op("orr")
        cost.op("str")
def wide_increment(words: int) -> int:
    c = Cost()
    _winc(c, words)
    return c.total
def wide_add(words: int) -> int:
    c = Cost()
    _wadd(c, words)
    return c.total
def wide_compare_branch(words: int) -> int:
    c = Cost()
    _wcmp(c, words)
    return c.total
def wide_shift7(src_words: int, dst_words: int) -> int:
    c = Cost()
    _wshift7(c, src_words, dst_words)
    return c.total
def _waddr(cost: Cost, index_words: int, addr_words: int) -> None:
    _wshift7(cost, index_words, addr_words)
    _wadd(cost, addr_words)
def wide_address(index_words: int, addr_words: int) -> int:
    c = Cost()
    _waddr(c, index_words, addr_words)
    return c.total
def _slot(cost: Cost) -> None:
    cost.op("ldrb")
    cost.op("add")
    cost.op("add")
def _copy_record(cost: Cost, addr_words: int) -> None:
    for _ in range(RECORD_BYTES // 8):
        _winc(cost, addr_words)
        _winc(cost, addr_words)
        cost.op("ldr")
        cost.op("str")
def _loop_step(cost: Cost, trips: int) -> None:
    if trips <= 0:
        return
    words = limb_words(trips)
    if words == 1:
        cost.op("add")
        cost.op("cmp")
        cost.op("b")
        return
    _winc(cost, words)
    _wcmp(cost, words)
def loop_control(trips: int) -> int:
    c = Cost()
    _loop_step(c, trips)
    return c.total
def loop_cost(trips: int, body: int) -> int:
    if trips <= 0:
        return 0
    return trips * (body + loop_control(trips))
def _trial_body() -> int:
    c = Cost()
    for name in ("udiv", "msub", "cmp", "b.ne", "add", "cmp", "b.lo"):
        c.op(name)
    if limb_words(1 << 32) != 1:
        raise SystemExit("trial index does not fit in one word")
    return c.total
TRIAL_BODY = _trial_body()
def radix_pass_parts(points: int) -> tuple[int, int, int]:
    if points <= 0:
        raise SystemExit("radix width of an empty list")
    cw = limb_words(points)
    aw = limb_words(points * RECORD_BYTES)
    addr = wide_address(cw, aw)
    slot = 3
    hist = addr + wide_increment(aw) + slot + wide_increment(cw)
    copy = (RECORD_BYTES // 8) * (2 * wide_increment(aw) + 2)
    scatter = addr + wide_increment(aw) + slot + addr + copy + wide_increment(cw)
    per = hist + scatter + 2 * loop_control(points)
    prefix = loop_cost(RADIX_BUCKETS, wide_add(cw))
    clear = loop_cost(RADIX_BUCKETS, slot + cw)
    return per, prefix, clear
def radix_sort_ins(points: int) -> int:
    if points <= 0:
        return 0
    per, prefix, clear = radix_pass_parts(points)
    one_pass = points * per + prefix + clear
    return loop_cost(RADIX_PASSES, one_pass)
def execute_radix_sort(records: list[bytes]) -> tuple[list[bytes], int]:
    n = len(records)
    if n <= 0:
        raise SystemExit("empty radix sort")
    if any(len(rec) != RECORD_BYTES for rec in records):
        raise SystemExit("record width")
    cw = limb_words(n)
    aw = limb_words(n * RECORD_BYTES)
    cost = Cost()
    src = list(records)
    for pass_i in range(RADIX_PASSES):
        counts = [0] * RADIX_BUCKETS
        for _ in range(RADIX_BUCKETS):
            _slot(cost)
            cost.op("str", cw)
            _loop_step(cost, RADIX_BUCKETS)
        for i in range(n):
            _waddr(cost, cw, aw)
            _winc(cost, aw)
            _slot(cost)
            _winc(cost, cw)
            _loop_step(cost, n)
            counts[src[i][pass_i]] += 1
        offsets = [0] * RADIX_BUCKETS
        running = 0
        for b in range(RADIX_BUCKETS):
            _wadd(cost, cw)
            _loop_step(cost, RADIX_BUCKETS)
            offsets[b] = running
            running += counts[b]
        if running != n:
            raise SystemExit("histogram does not sum to the record count")
        dest = [b""] * n
        for i in range(n):
            _waddr(cost, cw, aw)
            _winc(cost, aw)
            _slot(cost)
            _waddr(cost, cw, aw)
            _copy_record(cost, aw)
            _winc(cost, cw)
            _loop_step(cost, n)
            byte = src[i][pass_i]
            dest[offsets[byte]] = src[i]
            offsets[byte] += 1
        src = dest
        _loop_step(cost, RADIX_PASSES)
    return src, cost.total
def sort_self_check() -> int:
    recs = []
    for i in range(20):
        key = (i * 13 % 17).to_bytes(48, "little")
        body = i.to_bytes(8, "little") + bytes(RECORD_BYTES - 56)
        recs.append(key + body)
    recs.append(recs[2][:48] + b"\xab" + bytes(RECORD_BYTES - 49))
    out, counted = execute_radix_sort(recs)
    expect = sorted(recs, key=lambda rec: int.from_bytes(rec[:48], "little"))
    if out != expect:
        raise SystemExit("radix sort is not the sorted permutation")
    if sorted(out) != sorted(recs):
        raise SystemExit("radix sort dropped or invented a record")
    formula = radix_sort_ins(len(recs))
    if counted != formula:
        raise SystemExit(f"radix counter {counted} != formula {formula}")
    return counted
def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b
def ceil_sqrt(n: int) -> int:
    s = isqrt(n)
    return s if s * s == n else s + 1
def ceil_log2_rational(value: int, den: int = 1_000_000) -> tuple[int, int]:
    if value <= 0:
        raise ValueError("log2 of a non-positive integer")
    if value == 1:
        return 0, den
    if value & (value - 1) == 0:
        return (value.bit_length() - 1) * den, den
    L = value.bit_length() - 1
    yn = value - (1 << L)
    yd = value + (1 << L)
    scale_bits = 4096
    scale = 1 << scale_bits
    y_scaled = (yn * scale) // yd
    y2 = (y_scaled * y_scaled) // scale
    term = y_scaled
    acc = 0
    k = 0
    while term > 0 and k < 100000:
        acc += term // (2 * k + 1)
        term = (term * y2) // scale
        k += 1
    ln_m_scaled = 2 * acc + 4 * scale_bits
    ln2_lower_num, ln2_lower_den = 69314718055, 10**11
    logm_num = ln_m_scaled * ln2_lower_den
    logm_den = scale * ln2_lower_num
    numerator_unreduced = (L * logm_den + logm_num) * den
    num = (numerator_unreduced + logm_den - 1) // logm_den
    return num, den
def format_rational(num: int, den: int) -> str:
    whole, rem = divmod(num, den)
    return f"{whole}.{rem:0{len(str(den - 1))}d}"
def score_of(value: int) -> dict:
    num, den = ceil_log2_rational(value)
    approx = math.log2(value)
    gap = num / den - approx
    if gap < -1e-9 or gap >= 0.000002:
        raise SystemExit(f"log ceiling gap {gap} for {value}")
    return {"numerator": num, "denominator": den, "text": format_rational(num, den), "gap": gap}
def delta(field: Field, fn) -> int:
    start = field.cost.total
    out = fn()
    return field.cost.total - start, out
def measure_field(prime: int, limbs: int) -> dict:
    field = Field(prime, limbs, Cost())
    a = field.to_mont(3)
    b = field.to_mont(11)
    mul, _ = delta(field, lambda: field.mul(a, b))
    sqr, _ = delta(field, lambda: field.sqr(a))
    add, _ = delta(field, lambda: field.add(a, b))
    sub, _ = delta(field, lambda: field.sub(a, b))
    inv, got = delta(field, lambda: field.inv(a))
    if field.from_mont(field.mul(a, got)) != 1:
        raise SystemExit("measured inverse is wrong")
    e = prime - 2
    return {
        "mul": mul,
        "sqr": sqr,
        "add": add,
        "sub": sub,
        "inv": inv,
        "fermat_muls": (e.bit_length() - 1) + (e.bit_count() - 1),
        "setup": field._setup_cost,
        "limbs": limbs,
    }
def _affine_of(field: Field, point):
    if point is INF:
        return None
    return scale_z(field, point)
def measure_group(prime: int, limbs: int, g, dbl) -> dict:
    field = Field(prime, limbs, Cost())
    base = (field.to_mont(g[0]), field.to_mont(g[1]), field.to_mont(1))
    two = dbl(field, base)
    dbl_ins, _ = delta(field, lambda: dbl(field, base))
    add_ins, three = delta(field, lambda: add_2007_bl(field, two, base, dbl))
    aff = (field.to_mont(g[0]), field.to_mont(g[1]))
    madd_ins, _ = delta(field, lambda: madd_2007_bl(field, two, aff, dbl))
    neg_aff = (aff[0], field.sub(field.zero, aff[1]))
    neg_jac = (base[0], field.sub(field.zero, base[1]), base[2])
    madd_eq, eq_pt = delta(field, lambda: madd_2007_bl(field, base, aff, dbl))
    if _affine_of(field, eq_pt) != _affine_of(field, dbl(field, base)):
        raise SystemExit("mixed addition of a point to itself is not a double")
    madd_neg, neg_pt = delta(field, lambda: madd_2007_bl(field, base, neg_aff, dbl))
    if neg_pt is not INF:
        raise SystemExit("mixed addition of a point and its negative is not infinity")
    madd_inf, inf_pt = delta(field, lambda: madd_2007_bl(field, INF, aff, dbl))
    if _affine_of(field, inf_pt) != (g[0], g[1]):
        raise SystemExit("mixed addition of infinity did not return the affine point")
    madd_ainf, ainf_pt = delta(field, lambda: madd_2007_bl(field, base, INF, dbl))
    if _affine_of(field, ainf_pt) != (g[0], g[1]):
        raise SystemExit("mixed addition by an affine infinity changed the point")
    add_eq, add_eq_pt = delta(field, lambda: add_2007_bl(field, base, base, dbl))
    if _affine_of(field, add_eq_pt) != _affine_of(field, two):
        raise SystemExit("addition of a point to itself is not a double")
    add_neg, add_neg_pt = delta(field, lambda: add_2007_bl(field, base, neg_jac, dbl))
    if add_neg_pt is not INF:
        raise SystemExit("addition of a point and its negative is not infinity")
    add_inf, add_inf_pt = delta(field, lambda: add_2007_bl(field, INF, base, dbl))
    if add_inf_pt is not base:
        raise SystemExit("adding infinity on the left did not return the other point")
    two_aff = _affine_of(field, two)
    p_m = (field.to_mont(two_aff[0]), field.to_mont(two_aff[1]))
    q_m = aff
    a_mont = field.to_mont(0 if prime != P256_P else P256_A)
    add_aff, got = delta(field, lambda: affine_step(field, p_m, q_m, a_mont))
    three_aff = (field.from_mont(got[0]), field.from_mont(got[1]))
    if three_aff != _affine_of(field, three):
        raise SystemExit("affine addition does not match jacobian addition")
    dbl_aff, got_d = delta(field, lambda: affine_step(field, q_m, q_m, a_mont))
    if (field.from_mont(got_d[0]), field.from_mont(got_d[1])) != two_aff:
        raise SystemExit("affine doubling does not match jacobian doubling")
    neg_step, got_n = delta(field, lambda: affine_step(field, q_m, neg_aff, a_mont))
    if got_n is not INF:
        raise SystemExit("affine negation is not infinity")
    inf_step, got_i = delta(field, lambda: affine_step(field, INF, q_m, a_mont))
    if (field.from_mont(got_i[0]), field.from_mont(got_i[1])) != (g[0], g[1]):
        raise SystemExit("affine step from infinity did not return the point")
    both_inf, got_b = delta(field, lambda: affine_step(field, INF, INF, a_mont))
    if got_b is not INF:
        raise SystemExit("infinity plus infinity is not infinity")
    slope, slope_pt = delta(field, lambda: affine_add_mont(field, p_m, q_m))
    if (field.from_mont(slope_pt[0]), field.from_mont(slope_pt[1])) != three_aff:
        raise SystemExit("slope addition diverged")
    slope_d, slope_d_pt = delta(field, lambda: affine_dbl_mont(field, q_m, a_mont))
    if (field.from_mont(slope_d_pt[0]), field.from_mont(slope_d_pt[1])) != two_aff:
        raise SystemExit("slope doubling diverged")
    madd_price = max(madd_ins, madd_eq, madd_neg, madd_inf, madd_ainf)
    add_price = max(add_ins, add_eq, add_neg, add_inf)
    step = max(add_aff, dbl_aff, neg_step, inf_step, both_inf)
    return {
        "dbl": dbl_ins,
        "add": add_price,
        "add_unequal": add_ins,
        "madd": madd_price,
        "madd_unequal": madd_ins,
        "madd_equal": madd_eq,
        "affine_add": add_aff,
        "affine_dbl": dbl_aff,
        "step": step,
        "slope_add": slope,
        "slope_dbl": slope_d,
    }
def shake_cost() -> int:
    c = Cost()
    data = b"\x11" * 32 + (7).to_bytes(8, "big") + b"v3"
    digest = shake256(data, 64, c)
    if digest != hashlib.shake_256(data).digest(64):
        raise SystemExit("SHAKE256 mismatch")
    c0 = Cost()
    shake256(b"\x00" * 42, 64, c0)
    if c.total != c0.total:
        raise SystemExit("SHAKE cost depends on the message bytes")
    return c.total
def reduce_cost(modulus: int) -> int:
    c = Cost()
    digest = hashlib.shake_256(b"reduce").digest(64)
    got = reduce_hash(digest, modulus, c)
    if got != int.from_bytes(digest, "little") % modulus:
        raise SystemExit("reduction mismatch")
    return c.total
def partition_value(words: list[int]) -> int:
    h = 0
    for w in words:
        h ^= w & MASK
        h = (h * 0xBF58476D1CE4E5B9) & MASK
        h ^= h >> 30
    h ^= h >> 27
    h = (h * 0x94D049BB133111EB) & MASK
    h ^= h >> 31
    return h & 31
def partition_words(words: list[int], cost: Cost) -> int:
    h = 0
    cost.op("cmp")
    cost.op("b")
    cost.op("mov")
    for w in words:
        cost.op("eor")
        h ^= w & MASK
        cost.op("mul")
        h = (h * 0xBF58476D1CE4E5B9) & MASK
        cost.op("eor")
        cost.op("lsr")
        h ^= h >> 30
        cost.op("add")
        cost.op("cmp")
        cost.op("b")
    cost.op("eor")
    cost.op("lsr")
    h ^= h >> 27
    cost.op("mul")
    h = (h * 0x94D049BB133111EB) & MASK
    cost.op("eor")
    cost.op("lsr")
    h ^= h >> 31
    cost.op("and")
    return h & 31
def partition_cost(nlimbs: int) -> int:
    c = Cost()
    charged = partition_words([1] * nlimbs, c)
    if charged != partition_value([1] * nlimbs):
        raise SystemExit("partition schedule diverges from the value")
    return c.total
def equality_cost(nlimbs: int) -> int:
    c = Cost()
    c.op("ldr")
    c.op("ldr")
    c.op("cmp")
    c.op("b")
    for _ in range(nlimbs):
        c.op("ldr")
        c.op("ldr")
        c.op("cmp")
        c.op("b")
        c.op("add")
        c.op("cmp")
        c.op("b")
    return c.total
def factor_integer(n: int) -> list[int]:
    primes: list[int] = []
    def absorb(m: int) -> None:
        if m <= 1:
            return
        for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
            while m % p == 0:
                primes.append(p)
                m //= p
        if m == 1:
            return
        if _prime(m):
            primes.append(m)
            return
        d = _pollard(m)
        absorb(d)
        absorb(m // d)
    absorb(n)
    primes.sort()
    return primes
def _prime(n: int) -> bool:
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    for p in small:
        if n == p:
            return True
        if n % p == 0:
            return False
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in small:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True
def _pollard(n: int) -> int:
    def f(x, c):
        return (x * x + c) % n
    for c in range(1, 80):
        x = y = 2
        d = 1
        while d == 1:
            x = f(x, c)
            y = f(f(y, c), c)
            d = math.gcd(abs(x - y), n)
        if 1 < d < n:
            return d
    raise SystemExit(f"Pollard rho failed for {n}")
def exp_cost(bitlen: int, hamming: int, mul_ins: int) -> int:
    return mul_ins * (1 + bitlen + hamming) + 4 * bitlen
def batch_field_ins(npoints: int, mul_ins: int, inv_ins: int) -> int:
    if npoints <= 0:
        return 0
    return inv_ins + (7 * npoints - 3) * mul_ins
def _mask_body(limbs: int) -> int:
    return (2 * limbs + 1) + 4 * limbs
def batch_invert_ins(npoints: int, mul_ins: int, inv_ins: int, limbs: int) -> int:
    if npoints <= 0:
        return 0
    trips = npoints - 1
    sanitize = loop_cost(npoints, _mask_body(limbs))
    prefix = loop_cost(trips, mul_ins)
    backward = loop_cost(trips, 2 * mul_ins)
    affine = loop_cost(npoints, (2 * limbs + 1) + 4 * mul_ins)
    return inv_ins + sanitize + prefix + backward + affine
def execute_batch_affine(field: Field, npoints: int, zero_at: int = -1):
    zs = [field.zero if i == zero_at else field.to_mont(i + 3) for i in range(npoints)]
    xs = [field.to_mont(i + 11) for i in range(npoints)]
    ys = [field.to_mont(i + 19) for i in range(npoints)]
    start = field.cost.total
    one = field.one
    sanitized = []
    flags = []
    for i in range(npoints):
        field.cost.op("ldr", field.n)
        field.cost.op("cmp", field.n)
        field.cost.op("b")
        zero = all(w == 0 for w in zs[i])
        chosen = []
        for limb in range(field.n):
            field.cost.op("ldr")
            field.cost.op("ldr")
            field.cost.op("csel")
            field.cost.op("str")
            chosen.append(one[limb] if zero else zs[i][limb])
        sanitized.append(chosen)
        flags.append(zero)
        _loop_step(field.cost, npoints)
    products = [sanitized[0]]
    for i in range(1, npoints):
        products.append(field.mul(products[-1], sanitized[i]))
        _loop_step(field.cost, npoints - 1)
    acc = field.inv(products[-1])
    inverses = [None] * npoints
    for i in range(npoints - 1, 0, -1):
        inverses[i] = field.mul(acc, products[i - 1])
        acc = field.mul(acc, sanitized[i])
        _loop_step(field.cost, npoints - 1)
    inverses[0] = acc
    out = []
    for i in range(npoints):
        field.cost.op("ldr", field.n)
        field.cost.op("cmp", field.n)
        field.cost.op("b")
        z2 = field.sqr(inverses[i])
        z3 = field.mul(z2, inverses[i])
        x3 = field.mul(xs[i], z2)
        y3 = field.mul(ys[i], z3)
        out.append(None if flags[i] else (x3, y3))
        _loop_step(field.cost, npoints)
    return field.cost.total - start, out
def poly_check_ins(scalar_mul: int) -> int:
    field = Field(BLS_R, 4, Cost())
    start = field.cost.total
    x = field.to_mont(BLS_X)
    x2 = field.sqr(x)
    x4 = field.sqr(x2)
    got = field.add(field.sub(x4, x2), field.to_mont(1))
    if field.from_mont(got) != BLS_R % BLS_R:
        if field.from_mont(got) != 0:
            raise SystemExit("parameter polynomial does not match the order")
    return field.cost.total - start
def digit_extract_ins(window: int) -> int:
    span = 2 if window <= 65 else 3
    c = Cost()
    c.op("ldr", span)
    for _ in range(max(span - 1, 0)):
        c.op("lsr")
        c.op("lsl")
        c.op("orr")
    c.op("and", span)
    c.op("add")
    c.op("adc")
    return c.total
def recoded_digit(scalar: int, index: int, window: int) -> int:
    return ((scalar >> (index * window)) & ((1 << window) - 1)) + 1
def recoding_identity(scalar: int, window: int, bits: int) -> int:
    t = ceil_div(bits, window)
    base = 1 << window
    corr = (base**t - 1) // (base - 1)
    acc = 0
    for i in range(t):
        digit = recoded_digit(scalar, i, window)
        if digit < 1 or digit > base:
            raise SystemExit("recoded digit left 1..2^w")
        acc += digit * (base**i)
    if acc - corr != scalar:
        raise SystemExit("recoding does not reproduce the scalar")
    return corr
def cheon_at(window: int, g: dict, scalar: dict, shake: int, red: int) -> dict:
    if (BLS_R - 1) % (1 << 20) != 0:
        raise SystemExit("2^20 does not divide r-1")
    d = 1 << 20
    n_big = (BLS_R - 1) // d
    m = ceil_sqrt(n_big)
    v_max = ceil_div(n_big, m)
    if (m - 1) + m * (v_max - 1) < n_big - 1:
        raise SystemExit("BSGS range misses an index")
    points = m + v_max
    bits = 255
    t = ceil_div(bits, window)
    bases = 2
    base_pow = 1 << window
    if base_pow >= BLS_R:
        raise SystemExit("a table index meets the group order")
    corr = 0
    for sample in (0, 1, m, 1 << (window + 3), BLS_R - 1):
        corr = recoding_identity(sample, window, bits)
    if corr % BLS_R == 0:
        raise SystemExit("the correction multiple is the identity")
    comb_points = bases * t * base_pow
    tables = bases * t
    pre_doubles = bases * window * (t - 1) + tables
    pre_adds = tables * (base_pow - 2)
    online = points * t
    phase2_scalars = 2 * (1 << 10) + 4
    cw = limb_words(points)
    aw = limb_words(points * RECORD_BYTES)
    if cw != 2 or aw != 2:
        raise SystemExit("search-list counters or addresses are not two words")
    list_addr = wide_address(cw, aw)
    ccw = limb_words(comb_points)
    caw = limb_words(max(comb_points, 1) * RECORD_BYTES)
    table_addr = wide_address(ccw, caw)
    corr_addr = wide_address(1, 1)
    copy_ins = 3 * 2 * 6
    extract = digit_extract_ins(window)
    construct = (
        loop_cost(pre_doubles, g["dbl"] + table_addr)
        + loop_cost(pre_adds, g["add"] + table_addr)
        + loop_cost(tables, copy_ins + table_addr)
        + loop_cost(bases, loop_cost(t - 1, g["add"] + table_addr) + FIELD["bls"]["sub"] + corr_addr)
    )
    first = extract + table_addr + copy_ins
    rest = loop_cost(t - 1, g["madd"] + table_addr + extract)
    jac_store = 3 * 2 * 6
    corr_madd = g["madd"] + corr_addr + 2
    online_ins = loop_cost(points, first + rest + corr_madd + list_addr + jac_store)
    phase2_ins = loop_cost(phase2_scalars, loop_cost(254, g["dbl"] + g["madd"] + 2))
    group_ins = construct + online_ins + phase2_ins
    comb_norm = batch_invert_ins(comb_points, FIELD["bls"]["mul"], FIELD["bls"]["inv"], 6)
    list_norm = batch_invert_ins(points, FIELD["bls"]["mul"], FIELD["bls"]["inv"], 6)
    scalar_index = loop_cost(points, scalar["mul"] + list_addr)
    setup = {}
    zd = 1 << 20
    setup["zeta_d"] = exp_cost(zd.bit_length(), zd.bit_count(), scalar["mul"])
    setup["beta_m"] = exp_cost(m.bit_length(), m.bit_count(), scalar["mul"])
    setup["zeta_N"] = exp_cost(n_big.bit_length(), n_big.bit_count(), scalar["mul"])
    setup["s_unknown"] = exp_cost(255, 255, scalar["mul"])
    setup["k_combine"] = loop_cost(16, scalar["mul"])
    setup["three_inverses"] = loop_cost(3, scalar["inv"])
    setup_ins = sum(setup.values())
    factors = FACTORS
    distinct = sorted(set(factors))
    root_one = 0
    for q in distinct:
        exponent = (BLS_R - 1) // q
        root_one += exp_cost(exponent.bit_length(), exponent.bit_count(), scalar["mul"])
    root_ins = loop_cost(ROOT_CANDIDATES, root_one)
    factor_ins = 3 * (1 << 32) * TRIAL_BODY
    rebuild = loop_cost(len(factors), scalar["mul"]) + 4 * 4
    poly = POLY_INS
    loads = 2 * 2 * 6  # two affine points, two coordinates, six limbs
    to_mont = 4 * FIELD["bls"]["mul"]
    per_record, radix_prefix, radix_clear = radix_pass_parts(points)
    sort_ins = radix_sort_ins(points)
    scan_ins = loop_cost(points, equality_cost(12) + 2 * list_addr)
    found_ins = 2 * 2 * cw + 4
    field_setup = 2 * scalar["setup"] + FIELD["bls"]["setup"]
    total = (
        group_ins
        + comb_norm
        + list_norm
        + scalar_index
        + setup_ins
        + root_ins
        + factor_ins
        + rebuild
        + poly
        + loads
        + to_mont
        + sort_ins
        + scan_ins
        + found_ins
        + field_setup
    )
    c = C_PRICE["bls12-381-g1-ecdlp"]
    w = ceil_div(total, c)
    pre_ins = (
        construct
        + comb_norm
        + root_ins
        + factor_ins
        + rebuild
        + poly
        + field_setup
    )
    pre_w = ceil_div(pre_ins, c)
    bucket_bytes = RADIX_BUCKETS * cw * 8
    memory = comb_points * (RECORD_BYTES + 48) + points * RECORD_BYTES * 2 + 2 * bucket_bytes + 65536
    return {
        "window": window,
        "digits": t,
        "d": d,
        "N": n_big,
        "m": m,
        "v_max": v_max,
        "points": points,
        "comb_points": comb_points,
        "pre_doubles": pre_doubles,
        "pre_adds": pre_adds,
        "online": online,
        "group_ins": group_ins,
        "comb_norm": comb_norm,
        "list_norm": list_norm,
        "scalar_index": scalar_index,
        "setup_ins": setup_ins,
        "setup": setup,
        "root_ins": root_ins,
        "root_primes": len(distinct),
        "factor_ins": factor_ins,
        "rebuild": rebuild,
        "poly": poly,
        "loads": loads + to_mont,
        "sort_ins": sort_ins,
        "scan_ins": scan_ins,
        "field_setup": field_setup,
        "total_ins": total,
        "W": w,
        "pre_W": pre_w,
        "memory_bytes": memory,
        "phase2_scalars": phase2_scalars,
        "dbl": g["dbl"],
        "add": g["add"],
        "madd": g["madd"],
        "mul": FIELD["bls"]["mul"],
        "add_f": FIELD["bls"]["add"],
        "sub_f": FIELD["bls"]["sub"],
        "inv": FIELD["bls"]["inv"],
        "scalar_mul": scalar["mul"],
        "scalar_inv": scalar["inv"],
        "count_words": cw,
        "addr_words": aw,
        "construct": construct,
        "online_ins": online_ins,
        "phase2_ins": phase2_ins,
        "record_ins": per_record,
        "radix_prefix": radix_prefix,
        "radix_clear": radix_clear,
        "digit_ins": extract,
        "copy_ins": copy_ins,
        "found_ins": found_ins,
        "correction": corr,
    }
FIELD: dict = {}
FACTORS: list[int] = []
POLY_INS = 0
def best_cheon(g: dict, scalar: dict, shake: int, red: int) -> dict:
    best = None
    rows = []
    for window in (16, 24, 32, 48, 64, 80, 85, 96, 100, 104, 112, 127):
        row = cheon_at(window, g, scalar, shake, red)
        rows.append({"window": window, "digits": row["digits"], "W": row["W"]})
        if best is None or row["W"] < best["W"]:
            best = row
    best["window_scan"] = rows
    return best
def attach_scores(row: dict) -> dict:
    row["time"] = score_of(row["W"])
    row["memory"] = score_of(row["memory_bytes"])
    row["pre"] = score_of(row["pre_W"])
    return row
WALK_ORDERS = (
    (4294967291, 4),
    (68719476731, 2),
    (1099511627689, 1),
    (17592186044399, 1),
)
def cheon_algebra_check() -> int:
    r, d = 65537, 256
    if (r - 1) % d != 0:
        raise SystemExit("toy divisor")
    primes = factor_integer(r - 1)
    distinct = sorted(set(primes))
    zeta = 3
    if any(pow(zeta, (r - 1) // q, r) == 1 for q in distinct):
        raise SystemExit("toy root")
    beta = pow(zeta, d, r)
    n_big = (r - 1) // d
    if pow(beta, n_big, r) != 1 or any(pow(beta, n_big // q, r) == 1 for q in distinct):
        raise SystemExit("toy beta order")
    m = ceil_sqrt(n_big)
    v_max = ceil_div(n_big, m)
    if (m - 1) + m * (v_max - 1) < n_big - 1:
        raise SystemExit("toy phase-1 rectangle")
    m2 = ceil_sqrt(d)
    v2 = ceil_div(d, m2)
    if (m2 - 1) + m2 * (v2 - 1) < d - 1:
        raise SystemExit("toy phase-2 rectangle")
    beta_inv = pow(beta, -1, r)
    beta_m = pow(beta, m, r)
    gamma = pow(zeta, n_big, r)
    ginv = pow(gamma, -1, r)
    gm = pow(gamma, m2, r)
    wrong = 0
    for secret in range(1, r):
        pd = pow(secret, d, r)
        babies = {}
        scalar = pd
        for u in range(m):
            babies[scalar] = u
            scalar = scalar * beta_inv % r
        giant = 1
        k0 = None
        for v in range(v_max):
            if giant in babies:
                k0 = babies[giant] + m * v
                break
            giant = giant * beta_m % r
        if k0 is None:
            raise SystemExit("phase 1 missed a scalar")
        left = secret * pow(zeta, -k0, r) % r
        babies2 = {}
        scalar = left
        for u in range(m2):
            babies2[scalar] = u
            scalar = scalar * ginv % r
        giant = 1
        k1 = None
        for v in range(v2):
            if giant in babies2:
                k1 = babies2[giant] + m2 * v
                break
            giant = giant * gm % r
        if k1 is None or pow(zeta, k0 + k1 * n_big, r) != secret:
            raise SystemExit("phase 2 missed a scalar")
        if pow(secret, n_big, r) == 1:
            wrong += 1
    if wrong != n_big:
        raise SystemExit("A_1 phase-1 image is not the index-d subgroup")
    prod_d = 1 << 20
    prod_n = (BLS_R - 1) // prod_d
    prod_beta = pow(7, prod_d, BLS_R)
    if pow(prod_beta, prod_n, BLS_R) != 1:
        raise SystemExit("production beta^N")
    twos = 0
    rest = []
    for prime in FACTORS:
        if prime == 2 and twos < 20:
            twos += 1
        else:
            rest.append(prime)
    if twos != 20:
        raise SystemExit("production power of two")
    seen = set()
    for prime in rest:
        if prime in seen:
            continue
        seen.add(prime)
        if pow(prod_beta, prod_n // prime, BLS_R) == 1:
            raise SystemExit("production beta order collapses")
    prod_m = ceil_sqrt(prod_n)
    prod_v = ceil_div(prod_n, prod_m)
    if (prod_m - 1) + prod_m * (prod_v - 1) < prod_n - 1:
        raise SystemExit("production phase-1 rectangle")
    prod_m2 = ceil_sqrt(prod_d)
    prod_v2 = ceil_div(prod_d, prod_m2)
    if (prod_m2 - 1) + prod_m2 * (prod_v2 - 1) < prod_d - 1:
        raise SystemExit("production phase-2 rectangle")
    return r - 1
def bls_order_check() -> None:
    field = Field(BLS_P, 6, Cost())
    base = (field.to_mont(BLS_GX), field.to_mont(BLS_GY), field.to_mont(1))
    acc = INF
    k = BLS_R
    for i in range(k.bit_length() - 1, -1, -1):
        if acc is not INF:
            acc = dbl_2009_l(field, acc)
        if (k >> i) & 1:
            acc = base if acc is INF else add_2007_bl(field, acc, base)
    if acc is not INF and any(w != 0 for w in acc[2]):
        raise SystemExit("[r]G is not infinity")
def _points_match(field: Field, got, expect) -> bool:
    if (got is INF) != (expect is INF):
        return False
    if got is INF:
        return True
    return _affine_of(field, got) == _affine_of(field, expect)
def comb_recoding_check() -> None:
    field = Field(SECP_P, 4, Cost())
    dbl = dbl_2009_l
    window, bits = 4, 20
    base_pow = 1 << window
    t = ceil_div(bits, window)
    if ((1 >> window) & (base_pow - 1)) != 0 or recoded_digit(1, 1, window) != 1:
        raise SystemExit("scalar 1 was expected to have a zero ordinary digit")
    if any(recoded_digit(1, i, window) < 1 for i in range(t)):
        raise SystemExit("scalar 1 recoding still has a zero digit")
    gen = (field.to_mont(SECP_GX), field.to_mont(SECP_GY), field.to_mont(1))
    tables = []
    cursor = gen
    for _digit in range(t):
        table = [None] * (base_pow + 1)
        table[1] = cursor
        table[2] = dbl(field, cursor)
        acc = table[2]
        for e in range(3, base_pow + 1):
            acc = add_2007_bl(field, acc, cursor, dbl)
            table[e] = acc
        tables.append(table)
        for _shift in range(window):
            cursor = dbl(field, cursor)
    corr = tables[0][1]
    for i in range(1, t):
        corr = add_2007_bl(field, corr, tables[i][1], dbl)
    neg = _affine_of(field, (corr[0], field.sub(field.zero, corr[1]), corr[2]))
    if neg is None:
        raise SystemExit("correction point is infinity")
    neg_m = (field.to_mont(neg[0]), field.to_mont(neg[1]))
    affine_tables = []
    for table in tables:
        row = [None]
        for e in range(1, base_pow + 1):
            xy = _affine_of(field, table[e])
            if xy is None:
                raise SystemExit("recoded comb table contains the point at infinity")
            row.append((field.to_mont(xy[0]), field.to_mont(xy[1])))
        affine_tables.append(row)
    samples = (0, 1, base_pow, (1 << bits) - 1, 0xABCDE & ((1 << bits) - 1))
    for k in samples:
        recoding_identity(k, window, bits)
        acc = tables[0][recoded_digit(k, 0, window)]
        for i in range(1, t):
            acc = madd_2007_bl(field, acc, affine_tables[i][recoded_digit(k, i, window)], dbl)
        acc = madd_2007_bl(field, acc, neg_m, dbl)
        expect = jacobian_mul(field, k, gen, dbl, add_2007_bl)
        if not _points_match(field, acc, expect):
            raise SystemExit(f"recoded comb diverged on {k}")
def build(
    with_sample: bool = True,
    with_cheon: bool = True,
    with_walks: bool = True,
    with_ec: bool = True,
    with_gt: bool = True,
    with_sort: bool = True,
) -> dict:
    global FIELD, FACTORS, POLY_INS, FIELD_FOR
    print("measuring fields", flush=True)
    FIELD = {
        "secp": measure_field(SECP_P, 4),
        "p256": measure_field(P256_P, 4),
        "bls": measure_field(BLS_P, 6),
        "secp_r": measure_field(SECP_R, 4),
        "p256_r": measure_field(P256_R, 4),
        "bls_r": measure_field(BLS_R, 4),
    }
    for label, prime, key in (
        ("secp", SECP_P, "secp"),
        ("p256", P256_P, "p256"),
        ("bls", BLS_P, "bls"),
        ("bls_r", BLS_R, "bls_r"),
    ):
        e = prime - 2
        got = exp_cost(e.bit_length(), e.bit_count(), FIELD[key]["mul"])
        if got != FIELD[key]["inv"]:
            raise SystemExit(f"{label} exp_cost {got} != measured inv {FIELD[key]['inv']}")
    print("measuring groups", flush=True)
    g_secp = measure_group(SECP_P, 4, (SECP_GX, SECP_GY), dbl_2009_l)
    g_p256 = measure_group(P256_P, 4, (P256_GX, P256_GY), dbl_2001_b)
    g_bls = measure_group(BLS_P, 6, (BLS_GX, BLS_GY), dbl_2009_l)
    shake = shake_cost()
    red = reduce_cost(BLS_R)
    if reduce_cost(SECP_R) != red or reduce_cost(P256_R) != red:
        raise SystemExit("reduction cost changed with the order")
    parts = []
    if with_cheon:
        print("factoring parameter", flush=True)
        for label, value in ("X", BLS_X), ("X-1", BLS_X - 1), ("X+1", BLS_X + 1):
            fac = factor_integer(value)
            print(label, fac, flush=True)
            parts.append(fac)
        FACTORS = parts[0] + parts[0] + parts[1] + parts[2]
        product = 1
        for p in FACTORS:
            product *= p
        if product != BLS_R - 1:
            raise SystemExit("factor product is not r-1")
        distinct = sorted(set(FACTORS))
        if any(pow(7, (BLS_R - 1) // q, BLS_R) == 1 for q in distinct):
            raise SystemExit("7 is not a primitive root for these primes")
        POLY_INS = poly_check_ins(FIELD["bls_r"]["mul"])
    else:
        FACTORS = []
        distinct = []
        POLY_INS = 0
    for key in ("bls", "secp", "p256"):
        row = FIELD[key]
        if row["sqr"] != row["mul"]:
            raise SystemExit(f"{key} squaring is not the multiplication routine")
        prime = {"bls": BLS_P, "secp": SECP_P, "p256": P256_P}[key]
        limbs = 6 if key == "bls" else 4
        probe = Field(prime, limbs, Cost())
        for npoints in (1, 4):
            measured, outputs = execute_batch_affine(probe, npoints)
            if any(point is None for point in outputs):
                raise SystemExit(f"{key} nonzero batch wrote infinity")
            full = batch_invert_ins(npoints, row["mul"], row["inv"], limbs)
            if measured != full:
                raise SystemExit(f"{key} batch L={npoints} measured {measured} formula {full}")
            zprobe = Field(prime, limbs, Cost())
            zmeasured, zout = execute_batch_affine(zprobe, npoints, zero_at=npoints - 1)
            if zmeasured != full or zout[-1] is not None:
                raise SystemExit(f"{key} zero denominator was not masked")
            for i, point in enumerate(zout):
                if point is None:
                    continue
                z = i + 3
                got_x = zprobe.from_mont(point[0])
                got_y = zprobe.from_mont(point[1])
                if got_x != (i + 11) * pow(z, -2, prime) % prime:
                    raise SystemExit(f"{key} batch x poisoned by a zero denominator")
                if got_y != (i + 19) * pow(z, -3, prime) % prime:
                    raise SystemExit(f"{key} batch y poisoned by a zero denominator")
    if TRIAL_BODY != 7 or RADIX_PASSES != 48 or RADIX_BUCKETS != 256:
        raise SystemExit("trial body or radix pass count changed")
    if (wide_increment(2), wide_compare_branch(2), wide_increment(3), wide_compare_branch(3)) != (6, 7, 9, 10):
        raise SystemExit("wide primitive width changed")
    if wide_address(1, 1) != 9 or wide_address(2, 2) != 18:
        raise SystemExit("address schedule changed")
    if loop_control(1) != 3 or limb_words(1 << 32) != 1:
        raise SystemExit("one-word induction changed")
    if with_sort:
        print("radix sort", flush=True)
        sort_self_check()
    cheon = None
    if with_cheon:
        print("comb recoding", flush=True)
        comb_recoding_check()
        print("cheon windows", flush=True)
        cheon = best_cheon(g_bls, FIELD["bls_r"], shake, red)
        if cheon["count_words"] != 2 or cheon["addr_words"] != 2:
            raise SystemExit("Cheon search list is not two words")
        per, prefix, clear = radix_pass_parts(cheon["points"])
        if cheon["record_ins"] != per or cheon["radix_prefix"] != prefix or cheon["radix_clear"] != clear:
            raise SystemExit("Cheon sort schedule drifted from the radix formula")
        if cheon["sort_ins"] != radix_sort_ins(cheon["points"]):
            raise SystemExit("Cheon sort total drifted")
        if per - wide_increment(2) != 328:
            raise SystemExit(f"two-word record schedule {per} is not 328 plus a two-word increment")
        if cheon["window"] != 85:
            print("window moved", cheon["window"], [(r["window"], str(r["W"])) for r in cheon["window_scan"]], flush=True)
    FIELD_FOR = {
        "secp256k1-ecdlp": FIELD["secp"],
        "p256-ecdlp": FIELD["p256"],
        "bls12-381-g1-ecdlp-plain": FIELD["bls"],
    }
    walks = {}
    gt = {}
    fp2 = {"dbl": 0, "madd": 0, "add": 0, "prices": {}, "leaves_dbl": {}, "leaves_madd": {}}
    if with_walks and with_ec:
        print("walks", flush=True)
        walks = {
            "secp256k1-ecdlp": ec_walk(
                "secp256k1-ecdlp", SECP_R, SECP_P, g_secp, FIELD["secp_r"], shake, red, 8
            ),
            "p256-ecdlp": ec_walk("p256-ecdlp", P256_R, P256_P, g_p256, FIELD["p256_r"], shake, red, 8),
            "bls12-381-g1-ecdlp-plain": ec_walk(
                "bls12-381-g1-ecdlp-plain", BLS_R, BLS_P, g_bls, FIELD["bls_r"], shake, red, 12
            ),
        }
    if with_walks and with_gt:
        print("pairing walks", flush=True)
        fp2 = measure_fp2()
        g2_dbl = fp2["dbl"]
        g2_madd = fp2["madd"]
        gt = {
            "bls12-381-fixed-g1-invert-g2": gt_walk(
                "bls12-381-fixed-g1-invert-g2", g_bls, g2_dbl, g2_madd, FIELD["bls_r"], shake, red
            ),
            "bls12-381-fixed-g2-invert-g1": gt_walk(
                "bls12-381-fixed-g2-invert-g1", g_bls, g2_dbl, g2_madd, FIELD["bls_r"], shake, red
            ),
        }
    if with_walks:
        for row in [*walks.values(), *gt.values()]:
            if row["iter_words"] != 3:
                raise SystemExit(f"{row['track']} iteration index is not three words")
        h6 = 2747865721831794936778242923524729653278
        for row in gt.values():
            if row["H"] != h6:
                raise SystemExit(f"{row['track']} H moved to {row['H']}")
    for row in [*walks.values(), *gt.values()] + ([cheon] if cheon is not None else []):
        attach_scores(row)
    sample = None
    if with_sample:
        print("large-order sample", flush=True)
        sample = run_walk_sample(hashlib.sha256(b"curvesmash-v3-walk").digest())
        print(
            "sample",
            sample["hits"],
            "/",
            sample["trials"],
            "max",
            round(sample["max_ratio"], 4),
            "sec",
            round(sample["seconds"], 2),
            flush=True,
        )
        if sample["hits"] != sample["trials"]:
            raise SystemExit("a large-order trial missed the charged budget")
        if sample["max_ratio"] > FLOYD_MULT:
            raise SystemExit("ratio exceeded the charged multiplier")
    return {
        "field": {k: v for k, v in FIELD.items()},
        "group": {"secp": g_secp, "p256": g_p256, "bls": g_bls},
        "shake": shake,
        "reduce": red,
        "factors": (
            {"X": parts[0], "Xm": parts[1], "Xp": parts[2], "distinct": distinct}
            if parts
            else {"X": [], "Xm": [], "Xp": [], "distinct": []}
        ),
        "fp2": fp2,
        "cheon": cheon,
        "walks": walks,
        "gt": gt,
        "sample": sample,
        "pairing_one": pairing_bound() if gt else 0,
    }
def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, int):
        return str(obj) if obj.bit_length() > 60 else obj
    return obj
TRACK = 'bls12-381-g1-ecdlp'
CLAIM_W = 11740148960245930949733787335886895195
CLAIM_NUM = 123142791
CLAIM_DEN = 1000000
def entropy_check() -> None:
    seed = os.urandom(32)
    if len(seed) != 32:
        raise SystemExit("seed is not 32 bytes")
    seen = []
    for i in range(8):
        raw = seed + i.to_bytes(8, "big") + b"v3"
        digest = shake256(raw, 64)
        if digest != hashlib.shake_256(raw).digest(64):
            raise SystemExit("SHAKE256 does not match the library")
        value = int.from_bytes(digest, "little")
        if value.bit_length() <= 400:
            raise SystemExit("coefficient digest is narrower than the reduction argument")
        seen.append(value)
    if len(set(seen)) != len(seen):
        raise SystemExit("distinct indices produced one digest")
def main() -> None:
    raw = sys.stdin.buffer.read()
    if raw.strip() not in (b"", b"{}"):
        raise SystemExit("this program takes an empty object")
    real_out = sys.stdout
    sys.stdout = sys.stderr
    try:
        entropy_check()
        ledger = build(with_sample=False, with_walks=False)
        row = ledger["cheon"]
        if row["W"] != CLAIM_W or row["time"]["numerator"] != CLAIM_NUM or row["time"]["denominator"] != CLAIM_DEN:
            raise SystemExit("recomputed W does not match the claim")
        if row["count_words"] != 2 or row["addr_words"] != 2 or row["record_ins"] != 334:
            raise SystemExit("search-list width")
        covered = cheon_algebra_check()
        bls_order_check()
        if covered != 65536:
            raise SystemExit("coverage count")
        result = {
            "ok": True,
            "track": TRACK,
            "W": str(row["W"]),
            "time_log2": row["time"]["text"],
            "score_numerator": row["time"]["numerator"],
            "score_denominator": row["time"]["denominator"],
            "covered": covered,
            "phase1_bases": ["P0", "A_d"],
            "phase2_bases": ["A_1", "P0"],
            "counter_words": row["count_words"],
            "address_words": row["addr_words"],
            "record_ins": row["record_ins"],
        }
        payload = json.dumps(result, separators=(",", ":")) + "\n"
    finally:
        sys.stdout = real_out
    real_out.write(payload)
if __name__ == "__main__":
    main()
