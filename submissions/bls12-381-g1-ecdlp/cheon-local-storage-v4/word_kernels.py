"""Scalar word algorithms and conservative instruction certificates.

Python models word semantics for verification; these counters are analytical
ARM instruction bounds, not timings or measurements of Python/library calls.
"""
import hashlib
import random

MASK = (1 << 64)-1
MASK32 = (1 << 32)-1


def montgomery(a, b, modulus, limbs):
    assert limbs in (4,6)
    assert 0 <= a < modulus and 0 <= b < modulus and modulus < 1 << (64*limbs-1)
    words = lambda x: [(x >> (64*i)) & MASK for i in range(limbs)]
    aa, bb, nn = words(a), words(b), words(modulus)
    t = [0]*(2*limbs+1)
    n0 = -pow(modulus, -1, 1 << 64) & MASK
    instructions = 96  # loads, result stores, constants, prologue/epilogue allowance

    def addmul(x, y, k):
        nonlocal instructions
        product = x*y
        lo, hi = product & MASK, product >> 64
        instructions += 2  # MUL, UMULH
        total = t[k]+lo
        t[k], carry = total & MASK, total >> 64
        instructions += 1  # ADDS
        total = t[k+1]+hi+carry
        t[k+1], carry = total & MASK, total >> 64
        instructions += 1  # ADCS
        for z in range(k+2, 2*limbs+1):
            total = t[z]+carry
            t[z], carry = total & MASK, total >> 64
            instructions += 1  # ADCS with zero; unrolled even for carry=0
        assert carry == 0

    for i in range(limbs):
        for j in range(limbs):
            addmul(aa[i], bb[j], i+j)
    for i in range(limbs):
        q = t[i]*n0 & MASK
        instructions += 1
        for j in range(limbs):
            addmul(q, nn[j], i+j)
        assert t[i] == 0
    value = sum(t[i+limbs] << (64*i) for i in range(limbs+1))
    assert value < 2*modulus
    if value >= modulus:
        value -= modulus
    instructions += 2*limbs  # SUBS/SBCS then CSEL, charged on every call
    assert instructions <= (512 if limbs==4 else 1024)
    return value, instructions


def montgomery4(a,b,modulus):
    return montgomery(a,b,modulus,4)


K = [
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2]
INITIAL = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]


def sha256_48(message):
    assert len(message) == 48
    block = message+b'\x80'+bytes(7)+(384).to_bytes(8, 'big')
    assert len(block) == 64
    w = [int.from_bytes(block[4*i:4*i+4], 'big') for i in range(16)]
    rot = lambda x, n: ((x >> n) | (x << (32-n))) & MASK32
    for i in range(16, 64):
        s0 = rot(w[i-15], 7)^rot(w[i-15], 18)^(w[i-15] >> 3)
        s1 = rot(w[i-2], 17)^rot(w[i-2], 19)^(w[i-2] >> 10)
        w.append((w[i-16]+s0+w[i-7]+s1) & MASK32)
    a,b,c,d,e,f,g,h = INITIAL
    for i in range(64):
        s1 = rot(e,6)^rot(e,11)^rot(e,25)
        choice = (e & f)^((~e) & g)
        t1 = (h+s1+choice+K[i]+w[i]) & MASK32
        s0 = rot(a,2)^rot(a,13)^rot(a,22)
        majority = (a & b)^(a & c)^(b & c)
        t2 = (s0+majority) & MASK32
        h,g,f,e,d,c,b,a = g,f,e,(d+t1)&MASK32,c,b,a,(t1+t2)&MASK32
    digest = b''.join(((x+y)&MASK32).to_bytes(4,'big')
                      for x,y in zip(INITIAL,(a,b,c,d,e,f,g,h)))
    # Native W-register operations truncate to 32 bits without explicit ANDs.
    bound = 48*24 + 64*48 + 512
    assert bound <= 8192
    return digest, bound


def check(modulus, seed, samples=512, limbs=4):
    rng = random.Random(seed)
    rinv = pow(1 << (64*limbs), -1, modulus)
    maximum = 0
    values = [(0,0),(1,1),(modulus-1,modulus-1),(1,modulus-1)]
    values += [(rng.randrange(modulus),rng.randrange(modulus)) for _ in range(samples)]
    for a,b in values:
        out,cost = montgomery(a,b,modulus,limbs)
        assert out == a*b*rinv % modulus
        maximum = max(maximum,cost)
    for _ in range(samples):
        message = rng.getrandbits(384).to_bytes(48,'big')
        digest,cost = sha256_48(message)
        assert digest == hashlib.sha256(message).digest()
    return {'montgomery_vectors':len(values),'maximum_analytical_montgomery_instructions':maximum,
            'montgomery_limbs':limbs,'montgomery_budget':512 if limbs==4 else 1024,'sha256_vectors':samples,
            'analytical_sha256_bound':cost,'sha256_budget':8192,
            'scope':'Word semantics checked against independent big-integer/hashlib results; analytical scalar instruction bounds, not measured Python work'}
