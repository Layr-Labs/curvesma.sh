"""v15 shared kit: toy curve arithmetic (stdlib only), record layout, measured prices.

v15 fixes the v14 cls conflation: the production class key is the canonical
x-coordinate (negation maps (x,y) to (x,-y), so P and -P share a class, while
x and p-x are DIFFERENT points in general and are NOT merged). cls() now returns
the x-coordinate, matching record_pack()'s 32-byte key and the shipped sign rule.
Inverse is the extended-Euclid inverse (Python pow(a,-1,n)), valid for every
unit of any modulus; no Fermat claim is made anywhere.
"""
import hashlib
import math
import os

P_TOY = 1000003
B_TOY = 314159
N_TOY = 1000839            # fixture group order, COMPOSITE: 1000839 = 3 * 333613.


def is_probable_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def inv(a, n=None):
    """Extended-Euclid inverse; raises ValueError when gcd(a,n) != 1."""
    return pow(a % (N_TOY if n is None else n), -1, N_TOY if n is None else n)


def add(Pt, Qt, p=None):
    p = P_TOY if p is None else p
    if Pt is None:
        return Qt
    if Qt is None:
        return Pt
    x1, y1 = Pt
    x2, y2 = Qt
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if Pt == Qt:
        lam = (3 * x1 * x1 - 3) * inv(2 * y1, p) % p
    else:
        lam = (y2 - y1) * inv((x2 - x1) % p, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


def neg(Pt, p=None):
    if Pt is None:
        return None
    p = P_TOY if p is None else p
    return (Pt[0], (-Pt[1]) % p)


def mul(k, Pt, p=None):
    R = None
    while k:
        if k & 1:
            R = add(R, Pt, p)
        Pt = add(Pt, Pt, p)
        k >>= 1
    return R


def canon(Pt, p=None):
    """Shipped sign rule: negate iff y is odd, so the representative has even y."""
    if Pt is None:
        return None
    p = P_TOY if p is None else p
    return neg(Pt, p) if (Pt[1] & 1) else Pt


def cls(Pt):
    """Production class key = canonical x. Negation preserves x, so P and -P
    collide; x and p-x are different keys and are never merged."""
    if Pt is None:
        return None
    return Pt[0]


RECORD_BYTES = 48
TRAIL_BITS = 20
STEP_BITS = 108
IDENTITY_BITS = TRAIL_BITS + STEP_BITS


def identity(trail, step):
    assert 0 <= trail < (1 << TRAIL_BITS)
    assert 0 <= step < (1 << STEP_BITS)
    return (trail << STEP_BITS) | step


def unidentity(v):
    return v >> STEP_BITS, v & ((1 << STEP_BITS) - 1)


def record_pack(x, trail, step):
    return x.to_bytes(32, "big") + identity(trail, step).to_bytes(16, "big")


def record_unpack(buf):
    assert len(buf) == RECORD_BYTES
    return int.from_bytes(buf[:32], "big"), unidentity(int.from_bytes(buf[32:], "big"))


# ---- prices measured in evidence/insn_price.json (objdump of the shipped build) ----
FE_MUL_IR = 148            # fe_mul_obj span, arm64 -O2 -fPIC, p256_fe_v15.o
FE_ADD_IR = 27             # fe_add_obj span (out-of-line upper bound on the inlined form)
FE_SUB_IR = 27             # fe_sub_obj span
BATCH_INV_MULS = 385       # fe_inv chain: 256 squarings + popcount(p-2)=128 muls + 1 fixup
ADD_MULS = 3               # 2M+1S for a=-3: num*inv, lam^2, lam*(x1-x3)
NORM_C = 1000.0


def units(ir):
    return ir / NORM_C


def log2(x):
    return math.log2(x)


def base_point():
    for x in range(1, P_TOY):
        rhs = (x * x * x - 3 * x + B_TOY) % P_TOY
        if pow(rhs, (P_TOY - 1) // 2, P_TOY) == 1:
            y = pow(rhs, (P_TOY + 1) // 4, P_TOY)
            if (y * y - rhs) % P_TOY == 0:
                return (x, y)
    raise SystemExit("no base point")


def row_table(P, Q, R=64):
    tab = []
    for r in range(R):
        a = 1 + (r * 7 + 3) % R
        b = 1 + (r * 11 + 5) % R
        tab.append((a, b, add(mul(a, P), mul(b, Q))))
    return tab
