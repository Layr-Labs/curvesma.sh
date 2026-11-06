"""Hosted-only semantic fixtures. Does not measure A64 cost or attack success."""
import hashlib
import json
import math
import sys

MASK = (1 << 64) - 1


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


P = primes(80)
IV = tuple(math.isqrt(p << 128) & MASK for p in P[:8])
K = tuple(cube_root(p << 192) & MASK for p in P)
assert IV[0] == 0x6A09E667F3BCC908
assert K[0] == 0x428A2F98D728AE22


def rotr(x, n):
    return ((x >> n) | (x << (64 - n))) & MASK


def hash_schedule(message):
    bit_length = len(message) * 8
    padded = message + b'\x80'
    padded += b'\x00' * ((112 - len(padded)) % 128)
    padded += bit_length.to_bytes(16, 'big')
    state = IV
    for offset in range(0, len(padded), 128):
        block = padded[offset:offset + 128]
        w = [int.from_bytes(block[i:i+8], 'big') for i in range(0, 128, 8)]
        for i in range(16, 80):
            t0 = w[i-15]
            t1 = rotr(t0, 1)
            t1 ^= rotr(t0, 8)
            t1 ^= t0 >> 7
            t0 = w[i-2]
            t2 = rotr(t0, 19)
            t2 ^= rotr(t0, 61)
            t2 ^= t0 >> 6
            t0 = w[i-16]
            t1 = (t1 + t0) & MASK
            t0 = w[i-7]
            t1 = (t1 + t0) & MASK
            t1 = (t1 + t2) & MASK
            w.append(t1)
        work = state
        for i in range(80):
            a, b, c, d, e, f, g, h = work
            t0 = rotr(e, 14)
            t0 ^= rotr(e, 18)
            t0 ^= rotr(e, 41)
            t1 = f ^ g
            t1 &= e
            t1 ^= g
            t2 = (h + t0) & MASK
            t2 = (t2 + t1) & MASK
            t3 = w[i]
            t2 = (t2 + t3) & MASK
            t3 = K[i]
            t2 = (t2 + t3) & MASK
            t0 = rotr(a, 28)
            t0 ^= rotr(a, 34)
            t0 ^= rotr(a, 39)
            t1 = a | b
            t1 &= c
            t3 = a & b
            t1 |= t3
            t0 = (t0 + t1) & MASK
            d = (d + t2) & MASK
            h = (t2 + t0) & MASK
            work = (h, a, b, c, d, e, f, g)
        state = tuple((a+b) & MASK for a, b in zip(state, work))
    return b''.join(x.to_bytes(8, 'big') for x in state)


params = json.load(sys.stdin)
trials = params['trials']
assert isinstance(trials, int) and 1 <= trials <= 128
messages = [b'', b'abc']
for length in [1, 39, 40, 111, 112, 127, 128, 129, 168, 255, 256, 1000]:
    messages.append(bytes((i * 29 + length) % 256 for i in range(length)))
for i in range(trials):
    length = (i * 137 + 23) % 1025
    messages.append(hashlib.shake_256(str(i).encode()).digest(length))
for i, message in enumerate(messages):
    assert hash_schedule(message) == hashlib.sha512(message).digest(), i
print(json.dumps({'all_equal': True, 'messages': len(messages),
                  'scope': 'SHA-512 semantics only; no instruction trace or attack execution'}))
