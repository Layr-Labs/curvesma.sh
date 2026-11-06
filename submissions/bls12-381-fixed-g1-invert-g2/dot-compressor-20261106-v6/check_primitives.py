"""Hosted-only semantic fixtures. Does not measure A64 cost or attack success."""
import hashlib
import json
import math
import sys

MASK = (1 << 32) - 1


def primes(count):
    out = []
    n = 2
    while len(out) < count:
        if all(n % p for p in out if p * p <= n):
            out.append(n)
        n += 1
    return out


def cube_root(n):
    lo, hi = 0, 1 << ((n.bit_length() + 2) // 3)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    assert lo**3 <= n < (lo + 1)**3
    return lo


P = primes(64)
IV = tuple(math.isqrt(p << 64) & MASK for p in P[:8])
K = tuple(cube_root(p << 96) & MASK for p in P)
assert IV[0] == 0x6A09E667
assert K[0] == 0x428A2F98


def rotr(x, n):
    return ((x >> n) | (x << (32 - n))) & MASK


def hash_schedule(message):
    bit_length = len(message) * 8
    padded = message + b'\x80'
    padded += b'\x00' * ((56 - len(padded)) % 64)
    padded += bit_length.to_bytes(8, 'big')
    state = IV
    for offset in range(0, len(padded), 64):
        block = padded[offset:offset + 64]
        w = [int.from_bytes(block[i:i+4], 'big') for i in range(0, 64, 4)]
        for i in range(16, 64):
            t0 = w[i-15]
            t1 = rotr(t0, 7)
            t1 ^= rotr(t0, 18)
            t1 ^= t0 >> 3
            t0 = w[i-2]
            t2 = rotr(t0, 17)
            t2 ^= rotr(t0, 19)
            t2 ^= t0 >> 10
            t0 = w[i-16]
            t1 = (t1 + t0) & MASK
            t0 = w[i-7]
            t1 = (t1 + t0) & MASK
            t1 = (t1 + t2) & MASK
            w.append(t1)
        work = state
        for i in range(64):
            a, b, c, d, e, f, g, h = work
            t0 = rotr(e, 6)
            t0 ^= rotr(e, 11)
            t0 ^= rotr(e, 25)
            t1 = f ^ g
            t1 &= e
            t1 ^= g
            t2 = (h + t0) & MASK
            t2 = (t2 + t1) & MASK
            t3 = w[i]
            t2 = (t2 + t3) & MASK
            t3 = K[i]
            t2 = (t2 + t3) & MASK
            t0 = rotr(a, 2)
            t0 ^= rotr(a, 13)
            t0 ^= rotr(a, 22)
            t1 = a | b
            t1 &= c
            t3 = a & b
            t1 |= t3
            t0 = (t0 + t1) & MASK
            d = (d + t2) & MASK
            h = (t2 + t0) & MASK
            work = (h, a, b, c, d, e, f, g)
        state = tuple((a+b) & MASK for a, b in zip(state, work))
    return b''.join(x.to_bytes(4, 'big') for x in state)


params = json.load(sys.stdin)
trials = params['trials']
assert isinstance(trials, int) and 1 <= trials <= 128
messages = [b'', b'abc']
for length in [1, 39, 40, 55, 56, 63, 64, 65, 104, 127, 128, 1000]:
    messages.append(bytes((i * 29 + length) % 256 for i in range(length)))
for i in range(trials):
    length = (i * 137 + 23) % 1025
    messages.append(hashlib.shake_256(str(i).encode()).digest(length))
for i, message in enumerate(messages):
    assert hash_schedule(message) == hashlib.sha256(message).digest(), i

Q = (1 << 61) - 1
MASK64 = (1 << 64) - 1
MASK58 = (1 << 58) - 1
DIGIT_MASK = (1 << 60) - 1


def dot_evaluate(digits, key):
    lo, hi = digits[0], 0
    power = key
    for digit in digits[1:]:
        product = digit * power
        product_lo, product_hi = product & MASK64, product >> 64
        low_sum = lo + product_lo
        lo = low_sum & MASK64
        hi = hi + product_hi + (low_sum >> 64)
        assert hi < 1 << 63
        power = (power * key) % Q
    t0 = lo & Q
    t1 = lo >> 61
    t2 = hi & MASK58
    t1 |= t2 << 3
    t2 = hi >> 58
    t0 += t1
    t0 += t2
    assert t0 < 1 << 63
    t1 = t0 >> 61
    t0 &= Q
    t0 += t1
    if t0 >= Q:
        t0 -= Q
    assert 0 <= t0 < Q
    return t0


cases = [([0] * 39, [0, 1, 2, Q-2, Q-1]),
         ([DIGIT_MASK] * 39, [0, 1, 2, Q-2, Q-1]),
         ([DIGIT_MASK if i % 2 else 0 for i in range(39)],
          [0, 1, 2, Q-2, Q-1])]
for i in range(trials):
    raw = hashlib.shake_256(('polynomial' + str(i)).encode()).digest(44 * 8)
    words = [int.from_bytes(raw[j:j+8], 'big') for j in range(0, len(raw), 8)]
    cases.append(([v & DIGIT_MASK for v in words[:39]], [v % Q for v in words[39:]]))
checks = 0
for digits, keys in cases:
    for key in keys:
        horner = 0
        for digit in reversed(digits):
            horner = (horner * key + digit) % Q
        assert dot_evaluate(digits, key) == horner
        checks += 1
print(json.dumps({'all_equal': True, 'sha256_messages': len(messages),
                  'polynomial_evaluations': checks,
                  'scope': 'finite SHA and polynomial semantics only; no A64 trace or attack execution'}))
