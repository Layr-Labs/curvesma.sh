"""Hosted-only original semantic fixtures for standard RFC7693 BLAKE2s-256.
Algorithm constants are standard numerical data; no reference C code is copied.
No instruction trace, hash-model proof or production attack is asserted.
"""
import hashlib
import json
import math
import sys

MASK = (1 << 32) - 1
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19]
IV = tuple(math.isqrt(p << 64) & MASK for p in PRIMES)
SIGMA_HEX = ['000102030405060708090a0b0c0d0e0f', '0e0a0408090f0d06010c00020b070503', '0b080c0005020f0d0a0e030607010904', '070903010d0c0b0e0206050a04000f08', '0900050702040a0f0e010b0c0608030d', '020c060a000b0803040d07050f0e0109', '0c05010f0e0d040a000706030902080b', '0d0b070e0c01030905000f040806020a', '060f0e090b0300080c020d0701040a05', '0a020804070601050f0b090e030c0d00']
SIGMA = [tuple(bytes.fromhex(row)) for row in SIGMA_HEX]
assert all(sorted(row) == list(range(16)) for row in SIGMA)
QUADS = [(0,4,8,12),(1,5,9,13),(2,6,10,14),(3,7,11,15),
         (0,5,10,15),(1,6,11,12),(2,7,8,13),(3,4,9,14)]


def rotr(x, n):
    return ((x >> n) | (x << (32-n))) & MASK


def initial_state():
    h = list(IV)
    h[0] ^= 0x01010020
    return h


def compress(h, block, counter, final):
    assert len(block) == 64 and 0 <= counter < 1 << 64
    words = [int.from_bytes(block[i:i+4], 'little') for i in range(0,64,4)]
    v = list(h) + list(IV)
    v[12] ^= counter & MASK
    v[13] ^= counter >> 32
    if final:
        v[14] ^= MASK
    for row in SIGMA:
        for j,(a,b,c,d) in enumerate(QUADS):
            x,y = words[row[2*j]],words[row[2*j+1]]
            v[a] = (v[a]+v[b]+x) & MASK
            v[d] = rotr(v[d]^v[a],16)
            v[c] = (v[c]+v[d]) & MASK
            v[b] = rotr(v[b]^v[c],12)
            v[a] = (v[a]+v[b]+y) & MASK
            v[d] = rotr(v[d]^v[a],8)
            v[c] = (v[c]+v[d]) & MASK
            v[b] = rotr(v[b]^v[c],7)
    return [(h[i]^v[i]^v[i+8]) & MASK for i in range(8)]


def hash_schedule(message):
    h = initial_state()
    total = 0
    for offset in range(0,max(1,len(message)),64):
        chunk = message[offset:offset+64]
        total += len(chunk)
        block = chunk + bytes(64-len(chunk))
        h = compress(h,block,total,offset+64 >= len(message))
    return b''.join(word.to_bytes(4,'little') for word in h)


params = json.load(sys.stdin)
trials = params['trials']
assert isinstance(trials, int) and 1 <= trials <= 128
messages = [b'', b'abc']
for length in [1, 39, 40, 56, 63, 64, 65, 104, 127, 128, 688, 1000]:
    messages.append(bytes((i * 29 + length) % 256 for i in range(length)))
for i in range(trials):
    length = (i * 137 + 23) % 1025
    messages.append(hashlib.shake_256(str(i).encode()).digest(length))
for i, message in enumerate(messages):
    assert hash_schedule(message) == hashlib.blake2s(message).digest(), i

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

encodings = [bytes(289), bytes([255]) * 288 + bytes([1]), bytes([128]) + bytes(288)]
for i in range(trials):
    encodings.append(hashlib.shake_256(('encoding' + str(i)).encode()).digest(289))
extraction_checks = 0
for encoding in encodings:
    b_words = [int.from_bytes(encoding[j:j+8], 'big') for j in range(0, 288, 8)]
    sign_byte = encoding[-1]
    limbs = [((b_words[35] << 8) | sign_byte) & MASK64]
    for j in range(1, 36):
        limbs.append(((b_words[35-j] << 8) | (b_words[36-j] >> 56)) & MASK64)
    limbs.append(b_words[0] >> 56)
    integer = int.from_bytes(encoding, 'big')
    assert sum(word << (64*j) for j, word in enumerate(limbs)) == integer
    for l in range(39):
        j, t = divmod(60*l, 64)
        lower, upper = limbs[j], limbs[j+1]
        digit = (((upper << 64) | lower) >> t) & MASK64
        digit &= DIGIT_MASK
        assert digit == (integer >> (60*l)) & DIGIT_MASK
        extraction_checks += 1

walk_cases = 0
for i in range(trials + 3):
    seed = hashlib.shake_256(('seed' + str(i)).encode()).digest(32)
    payload = hashlib.shake_256(('payload' + str(i)).encode()).digest(40)
    prefix = b'FAPI-RHO-WALK-v1' + seed + bytes(16)
    assert len(prefix) == 64
    message = prefix + payload
    assert len(message) == 104
    prepared = payload + bytes(24)
    assert len(prepared) == 64
    mid = compress(initial_state(),prefix,64,False)
    final_state = compress(mid,prepared,104,True)
    digest = b''.join(word.to_bytes(4,'little') for word in final_state)
    assert digest == hash_schedule(message) == hashlib.blake2s(message).digest()
    words = [int.from_bytes(digest[j:j+4], 'little') for j in range(0,32,4)]
    rev7 = int.from_bytes(words[7].to_bytes(4,'little'),'big')
    label = ((words[6] >> 24) << 32) | rev7
    assert label == int.from_bytes(digest[-5:],'big')
    walk_cases += 1

FIELD_P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
FIELD_R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
assert pow(FIELD_P, 6, FIELD_R) == FIELD_R-1
HALF = (FIELD_P + 1) // 2


def sign_fast(coefficients):
    for x in coefficients:
        if x:
            return int(x >= HALF)
    return 0


def compare_keys(left, right):
    for pair in range(18):
        for j in (2*pair, 2*pair+1):
            if left[j] != right[j]:
                return -1 if left[j] < right[j] else 1
    if left[36] != right[36]:
        return -1 if left[36] < right[36] else 1
    return 0


assert pow(2, 50, 101) == 100
norm_keys = {}
for a in range(101):
    for b in range(101):
        if (a*a - 2*b*b) % 101 == 1:
            sign = int(a > (-a) % 101)
            key = (b, sign)
            assert key not in norm_keys
            norm_keys[key] = (a, b)
assert len(norm_keys) == 102

vectors = [([0]*6, [0]*6), ([HALF]+[0]*5, [FIELD_P-1]*6),
           ([FIELD_P-1]*6, [1]*6)]
for i in range(trials):
    raw = hashlib.shake_256(('canonical' + str(i)).encode()).digest(12*48)
    c = [int.from_bytes(raw[j:j+48], 'big') % FIELD_P for j in range(0, len(raw), 48)]
    vectors.append((c[:6], c[6:]))
key_cases = []
for a, b in vectors:
    neg = [(-x) % FIELD_P for x in a]
    sign = sign_fast(a)
    assert sign == int(tuple(a) > tuple(neg))
    neg_sign = sign_fast(neg)
    assert neg_sign == int(tuple(neg) > tuple(a))
    encoded_b = b''.join(x.to_bytes(48, 'big') for x in b)
    physical = b''.join(x.to_bytes(48, 'little') for x in b)
    words = []
    for i in range(36):
        offset = 48*(i//6)+8*(5-(i%6))
        value = int.from_bytes(physical[offset:offset+8], 'little')
        assert value == int.from_bytes(encoded_b[8*i:8*i+8], 'big')
        words.append(value)
    key = words + [sign]
    neg_key = words + [neg_sign]
    assert (key == neg_key) == (a == neg)
    for candidate in (key, neg_key):
        stored = b''.join(w.to_bytes(8, 'little') for w in candidate) + bytes(8)
        assert len(stored) == 304
        assert [int.from_bytes(stored[j:j+8], 'little') for j in range(0, 296, 8)] == candidate
        key_cases.append(candidate)
key_comparisons = 0
for left in key_cases:
    for right in key_cases:
        expected = (tuple(left) > tuple(right)) - (tuple(left) < tuple(right))
        assert compare_keys(left, right) == expected
        key_comparisons += 1
print(json.dumps({'all_equal': True, 'blake2s_messages': len(messages),
                  'polynomial_evaluations': checks,
                  'digit_extractions': extraction_checks,
                  'walk_block_and_label_cases': walk_cases,
                  'public_frobenius_parameter_check': True,
                  'toy_norm_one_points': len(norm_keys),
                  'canonical_sign_vectors': len(vectors),
                  'exact_key_comparisons': key_comparisons,
                  'scope': 'finite primitive/toy encoding checks only; not production GT success, A64 trace or hash-model proof'}))
