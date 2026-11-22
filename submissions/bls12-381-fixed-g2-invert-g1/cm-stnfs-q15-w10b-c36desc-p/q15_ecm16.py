"""q15_ecm16.py — Exact port of c36p/ecm16.py for Q(zeta_15) x c36p-DESCENT."""
import math
from math import ceil, isqrt, log, log1p, log2, sqrt
from functools import lru_cache
import q15_descent16 as D

MONT5 = D.MONT5
TORSION = 12
WHEELS = (30030,)


def mulmod_instr(nbits):
    L = max(5, ceil(nbits / 64))
    return MONT5 * (L / 5.0) ** 2, L


def gcd_instr(nbits):
    L = max(5, ceil(nbits / 64))
    return 2.0 * nbits * (12 * L + 16)


def primes_upto(n):
    a = bytearray(b"\x01") * (n + 1)
    a[0:2] = b"\x00\x00"
    for p in range(2, isqrt(n) + 1):
        if a[p]:
            a[p * p :: p] = b"\x00" * ((n - p * p) // p + 1)
    return [p for p in range(2, n + 1) if a[p]]


def pi_upper(x):
    return 1.26 * x / log(x) if x >= 17 else 7.0


@lru_cache(None)
def semismooth(fbits, b1, b2, steps=256):
    m = fbits - log2(TORSION)
    base = log2(b1)
    top = log2(b2)
    acc = D.rho_lo(m / base)
    if top > base:
        for i in range(steps):
            lo = base + (top - base) * i / steps
            hi = base + (top - base) * (i + 1) / steps
            acc += D.rho_lo((m - lo) / base) * log(hi / lo)
    return min(acc, 1.0 - 1e-15)


@lru_cache(None)
def power_correction(fbits, b1, b2):
    m = fbits - log2(TORSION)
    u = m / log2(b1)
    tilt = D.rho_hi(max(u - 2.0, 1.0)) / max(D.rho_lo(u), 1e-300)
    y = sqrt(b1)
    return max(0.0, 1.0 - tilt * 2.0 / (y * log(y)))


def stage_ops(b1, b2, wheel):
    s1 = ceil(1.5 * b1) * 11
    phi = wheel
    for q in (2, 3, 5, 7, 11, 13):
        if wheel % q == 0:
            phi = phi // q * (q - 1)
    J = phi // 2
    G = ceil(b2 / wheel) + 1
    s2 = 6 * (wheel // 2 + G) + 3 * J * G if b2 > b1 else 0
    return 1.1 * (s1 + s2), s1, s2, J, G


@lru_cache(None)
def schedule(fbits, nbits, success=0.999):
    nbits = int(ceil(nbits))
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
            curves = ceil(log1p(-success) / log1p(-p1))
            for wheel in WHEELS if ratio > 1 else (0,):
                if wheel and wheel > b2:
                    continue
                ops, s1, s2, J, G = (
                    stage_ops(b1, b2, wheel)
                    if wheel
                    else (1.1 * ceil(1.5 * b1) * 11, ceil(1.5 * b1) * 11, 0, 0, 0)
                )
                per_curve = ops * mul + 2 * gcd
                redo = p1 * (ops * mul + (pi_upper(b1) + G) * gcd)
                cost = curves * (per_curve + redo)
                if best is None or cost < best["instructions"]:
                    best = dict(
                        fbits=fbits,
                        nbits=nbits,
                        B1=b1,
                        B2=b2,
                        wheel=wheel,
                        curves=curves,
                        p_curve=p1,
                        stage1_M=s1,
                        stage2_M=s2,
                        baby=J,
                        giant=G,
                        mul_instr=mul,
                        limbs=L,
                        gcd_instr=gcd,
                        instructions=cost,
                        log2_instructions=log2(cost),
                    )
    return best


def ecm_instr(fbits, nbits, success):
    return schedule(fbits, nbits, success)["instructions"]
