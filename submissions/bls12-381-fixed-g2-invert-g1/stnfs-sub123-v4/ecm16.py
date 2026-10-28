import math
from functools import lru_cache
import descent16 as D
MONT5 = D.MONT5
TORSION = 12
WHEELS = (30030,)
def mulmod_instr(nbits):
    L = max(5, math.ceil(nbits / 64))
    return MONT5 * (L / 5.0) ** 2, L
def gcd_instr(nbits):
    L = max(5, math.ceil(nbits / 64))
    return 2.0 * nbits * (12 * L + 16)
def primes_upto(n):
    a = bytearray(b'\1') * (n + 1)
    a[0:2] = b'\0\0'
    for p in range(2, math.isqrt(n) + 1):
        if a[p]:
            a[p * p::p] = b'\0' * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if a[p]]
def pi_upper(x):
    return 1.26 * x / math.log(x) if x >= 17 else 7.0
@lru_cache(None)
def semismooth(fbits, b1, b2, steps=256):
    m = fbits - math.log2(TORSION)
    base = math.log2(b1)
    top = math.log2(b2)
    acc = D.rho_lo(m / base)
    if top > base:
        for i in range(steps):
            lo = base + (top - base) * i / steps
            hi = base + (top - base) * (i + 1) / steps
            acc += D.rho_lo((m - lo) / base) * math.log(hi / lo)
    return min(acc, 1 - 1e-15)
@lru_cache(None)
def power_correction(fbits, b1, b2):
    m = fbits - math.log2(TORSION)
    u = m / math.log2(b1)
    tilt = D.rho_hi(max(u - 2.0, 1.0)) / max(D.rho_lo(u), 1e-300)
    y = math.sqrt(b1)
    return max(0.0, 1.0 - tilt * 2.0 / (y * math.log(y)))
def stage_ops(b1, b2, wheel):
    s1 = math.ceil(1.5 * b1) * 11
    phi = wheel
    for q in (2, 3, 5, 7, 11, 13):
        if wheel % q == 0:
            phi = phi // q * (q - 1)
    J = phi // 2
    G = math.ceil(b2 / wheel) + 1
    s2 = 6 * (wheel // 2 + G) + 3 * J * G if b2 > b1 else 0
    return 1.1 * (s1 + s2), s1, s2, J, G
@lru_cache(None)
def schedule(fbits, nbits, success=0.999):
    nbits = int(math.ceil(nbits))
    mul, L = mulmod_instr(nbits)
    gcd = gcd_instr(nbits)
    best = None
    for k in range(10, 37, 2):
        b1 = 2 ** k
        for ratio in (1, 10, 30, 100, 300):
            b2 = b1 * ratio
            p1 = semismooth(fbits, b1, b2) * power_correction(fbits, b1, b2)
            if p1 <= 1e-12:
                continue
            curves = math.ceil(math.log1p(-success) / math.log1p(-p1))
            for wheel in (WHEELS if ratio > 1 else (0,)):
                if wheel and wheel > b2:
                    continue
                ops, s1, s2, J, G = stage_ops(b1, b2, wheel) if wheel else (1.1 * math.ceil(1.5 * b1) * 11, math.ceil(1.5 * b1) * 11, 0, 0, 0)
                per_curve = ops * mul + 2 * gcd
                redo = p1 * (ops * mul + (pi_upper(b1) + G) * gcd)
                cost = curves * (per_curve + redo)
                if best is None or cost < best['instructions']:
                    best = dict(fbits=fbits, nbits=nbits, B1=b1, B2=b2, wheel=wheel, curves=curves, p_curve=p1,
                                stage1_M=s1, stage2_M=s2, baby=J, giant=G, mul_instr=mul, limbs=L, gcd_instr=gcd,
                                instructions=cost, log2_instructions=math.log2(cost))
    return best
def ecm_instr(fbits, nbits, success):
    return schedule(fbits, nbits, success)['instructions']
def xdbl(P_, a24, n):
    x, z = P_
    u = (x + z) ** 2 % n; v = (x - z) ** 2 % n; w = (u - v) % n
    return u * v % n, w * (v + a24 * w) % n
def xadd(P_, Q, diff, n):
    x, z = P_; y, w = Q; dx, dz = diff
    u = (x + z) * (y - w) % n; v = (x - z) * (y + w) % n
    return dz * (u + v) ** 2 % n, dx * (u - v) ** 2 % n
def ladder(P_, k, a24, n):
    if k == 0:
        return (1, 0)
    R0 = P_; R1 = xdbl(P_, a24, n)
    for bit in bin(k)[3:]:
        if bit == '0':
            R1 = xadd(R0, R1, P_, n); R0 = xdbl(R0, a24, n)
        else:
            R0 = xadd(R0, R1, P_, n); R1 = xdbl(R1, a24, n)
    return R0
def curve(p, sigma):
    u = (sigma * sigma - 5) % p; v = 4 * sigma % p
    den = 16 * pow(u, 3, p) * v % p
    if not den:
        return None
    a24 = pow(v - u, 3, p) * (3 * u + v) * pow(den, -1, p) % p
    if (4 * a24 - 2) ** 2 % p == 4:
        return None
    return (pow(u, 3, p), pow(v, 3, p)), a24
def stage1(p, sigma, b1):
    c = curve(p, sigma)
    if c is None:
        return None
    Q, a24 = c
    for q in primes_upto(b1):
        power = q
        while power * q <= b1:
            power *= q
        Q = ladder(Q, power, a24, p)
    return Q, a24
def attempt(p, sigma, b1, b2, wheel=210):
    st = stage1(p, sigma, b1)
    if st is None:
        return None
    Q, a24 = st
    if Q[1] == 0:
        return 1
    baby = {j: ladder(Q, j, a24, p) for j in range(1, wheel // 2 + 1) if math.gcd(j, wheel) == 1}
    step = ladder(Q, wheel, a24, p)
    prev = (1, 0); cur = step
    for i in range(1, math.ceil(b2 / wheel) + 1):
        x, z = cur
        for j, (y, w) in baby.items():
            if not (b1 < i * wheel - j <= b2 or b1 < i * wheel + j <= b2):
                continue
            if (x * w - y * z) % p == 0:
                return 2
        nxt = xdbl(step, a24, p) if i == 1 else xadd(cur, step, prev, p)
        prev, cur = cur, nxt
    return 0
def toy(trials, seed, sizes=((36, 1024, 50), (40, 2048, 50), (44, 4096, 50))):
    import random
    rng = random.Random(seed); out = []
    for bits, b1, ratio in sizes:
        b2 = b1 * ratio; counts = [0, 0, 0]; pred = 0.0; predc = 0.0; n = 0
        for _ in range(trials):
            p = rng.randrange(2 ** (bits - 1), 2 ** bits) | 1
            while not D_probable(p):
                p += 2
            r = attempt(p, rng.randrange(6, p), b1, b2)
            if r is None:
                continue
            counts[r] += 1; n += 1
            pred += semismooth(math.log2(p), b1, b2)
            predc += semismooth(math.log2(p), b1, b2) * power_correction(math.log2(p), b1, b2)
        out.append(dict(bits=bits, B1=b1, B2=b2, trials=n, stage1=counts[1], stage2=counts[2],
                        empirical=(counts[1] + counts[2]) / n, model=pred / n, model_corrected=predc / n))
    return out
def D_probable(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d = n - 1; s = 0
    while d % 2 == 0:
        s += 1; d //= 2
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        x = pow(a % n, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True
