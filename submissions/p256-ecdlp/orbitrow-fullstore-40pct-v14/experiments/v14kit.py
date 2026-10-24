"""v14 shared kit: toy curve arithmetic (stdlib only), record layout, cost constants.

No numpy, no third-party imports. Everything here is deterministic and replayable.
"""
import hashlib
import math
import os

# ---------------------------------------------------------------- toy curve
# y^2 = x^3 - 3x + b over F_p with p prime; a = -3 matches NIST P-256's curve shape.
P_TOY = 1000003
B_TOY = 314159


def pt_order_count():
    """Count points by Legendre-sum (used once at build time to pick b)."""
    p = P_TOY
    total = 1
    for x in range(p):
        rhs = (x * x * x - 3 * x + B_TOY) % p
        if rhs == 0:
            total += 1
        else:
            total += 2 if pow(rhs, (p - 1) // 2, p) == 1 else 0
    return total


def inv(a, p=None):
    p = P_TOY if p is None else p
    return pow(a % p, p - 2, p)


def add(Pt, Qt, p=None):
    """Affine addition on y^2=x^3-3x+b. Complete enough for the toy fixtures."""
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
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


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
    """Class representative: negate iff y is odd (the shipped sign rule)."""
    if Pt is None:
        return None
    p = P_TOY if p is None else p
    return neg(Pt, p) if (Pt[1] & 1) else Pt


def cls(Pt):
    """Negation-invariant class key: (min(x, p-x), y^2 parity tag is unnecessary for x)."""
    if Pt is None:
        return None
    return min(Pt[0], P_TOY - Pt[0])


# ---------------------------------------------------------------- record layout
# The 40-byte v13 record could not identify a visit at production scale.
# v14 sizes the record from the requirement: 20-bit trail id + 108-bit step id = 128 bits.
RECORD_BYTES = 48           # 32-byte class key (x, 256 bit) + 16-byte visit identity
TRAIL_BITS = 20
STEP_BITS = 108
IDENTITY_BITS = TRAIL_BITS + STEP_BITS


def identity(trail, step):
    assert 0 <= trail < (1 << TRAIL_BITS)
    assert 0 <= step < (1 << STEP_BITS)
    return (trail << STEP_BITS) | step


def unidentity(v):
    return v >> STEP_BITS, v & ((1 << STEP_BITS) - 1)


def record_pack(x, trail, step, p=None):
    """48-byte record: 32-byte big-endian x || 16-byte big-endian identity."""
    p = P_TOY if p is None else p
    return x.to_bytes(32, "big") + identity(trail, step).to_bytes(16, "big")


def record_unpack(buf):
    assert len(buf) == RECORD_BYTES
    x = int.from_bytes(buf[:32], "big")
    return x, unidentity(int.from_bytes(buf[32:], "big"))


# ---------------------------------------------------------------- cost rows
# P-256 span prices, pinned in the track's own accepted accounting lineage
# (v10p evidence/stepspans_core.json, object sha256 96b736ef..., TEXT-only nm spans).
FE_MUL_IR = 147            # static span of the shipped p256 fe_mul
FE_SQR_IR = 147            # charged as a multiplication (conservative, S <= M)
FE_ADD_IR = 8              # 4 limbs x (add + carry)
STEP_FORWARD_IR = 239      # shipped step_forward span
STEP_CLOSE_IR = 694        # shipped step_close span
BATCH_INV_ITER_IR = 293    # per-walker back-substitution cost
INV_MULS = 385             # Fermat chain for p-2: 256 squarings + 129 multiplications
NORM_C = 1000.0            # protected track normalization C (p256-ecdlp = 1000)


def units(ir):
    return ir / NORM_C


def log2(x):
    return math.log2(x)
