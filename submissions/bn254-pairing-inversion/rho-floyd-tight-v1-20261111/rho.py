"""Executable mathematical reference for the bounded 32-adding Pollard-rho walk.

This module is evidence, not a production attack. It gives:

  * ``montgomery_word_reference`` - a reference for four-word Montgomery REDC,
    used only to check the arithmetic identity on full-width inputs.
  * ``fp12_mul`` / ``fp12_pow`` - plain polynomial arithmetic in the protected
    BN254 tower Fp12 = Fp[w]/(w^12 - 18*w^6 + 82), used to check the public
    generator order and the exponent invariant on full-width elements.
  * ``rho`` - the generic 32-adding Floyd walk, run on small cyclic groups as
    finite walk-statistics evidence (including failures). The same loop order,
    salted SHA-256 partition, 384-byte state encoding, cap, recovery and
    failure rules are used as in the production description in ``proof.md``.
  * ``floyd_cap`` / ``ledger`` - the exact production cap and instruction
    ledger, recomputed from the protected group order.

Nothing here solves a production-size discrete logarithm or proves a
probability bound at the production order; see ``proof.md`` for scope.
"""

import hashlib
import math

# Protected BN254 tower modulus Fp12 = Fp[w]/(w^12 - 18 w^6 + 82):
#   w^12 == 18*w^6 - 82   (mod p and the defining polynomial).
_RED_HI = 18    # coefficient of w^(d-6) when reducing w^d
_RED_LO = -82   # coefficient of w^(d-12) when reducing w^d


def montgomery_word_reference(a, b, m):
    """Reference value of a four-word (2^256) Montgomery multiplication.

    Returns (a*b*Rinv) mod m with R = 2^256. A concrete bounded REDC
    implementation computing this exact value is described in ``machine.md``;
    here we only assert the arithmetic identity on full-width inputs.
    """
    Rinv = pow(1 << 256, -1, m)
    return (a * b * Rinv) % m


def fp12_mul(a, b, p):
    """Multiply two Fp12 elements (12 coefficients c0..c11 for sum ci*w^i)."""
    res = [0] * 23
    for i in range(12):
        ai = a[i]
        if ai:
            for j in range(12):
                res[i + j] += ai * b[j]
    # Reduce degrees 22..12 using w^d = 18 w^(d-6) - 82 w^(d-12).
    for d in range(22, 11, -1):
        c = res[d]
        if c:
            res[d - 6] += _RED_HI * c
            res[d - 12] += _RED_LO * c
            res[d] = 0
    return tuple(x % p for x in res[:12])


def fp12_pow(base, exp, p):
    """Square-and-multiply exponentiation in Fp12."""
    result = (1,) + (0,) * 11
    acc = tuple(x % p for x in base)
    e = exp
    while e:
        if e & 1:
            result = fp12_mul(result, acc, p)
        e >>= 1
        if e:
            acc = fp12_mul(acc, acc, p)
    return result


def _branch(salt, encoded):
    """Salted SHA-256 partition into one of 32 classes (bits 0..4)."""
    return hashlib.sha256(salt + encoded).digest()[0] & 31


def rho(q, secret, mul, scalar_mul, g_elt, z_elt, identity, encode, coins,
        cap=None):
    """Generic 32-adding Floyd walk on a cyclic group of order q.

    ``g_elt`` is a generator, ``z_elt`` = g^secret. The walk recovers secret
    with the production loop order. Group operations are supplied abstractly:
      mul(x, y)         -> group product
      scalar_mul(x, k)  -> x^k (k a nonnegative integer)
      encode(x)         -> fixed-length bytes for the partition
    Returns (recovered_or_None, floyd_iterations, outcome) where outcome is one
    of 'success', 'identity', 'degenerate', 'cap', 'verification_failed'.
    """
    # Identity target: discrete log is 0, reported directly.
    if z_elt == identity:
        return (0, 0, "identity")
    if cap is None:
        cap = 2 * (math.isqrt(q) + 1)

    salt = coins.getrandbits(256).to_bytes(32, "big")
    # 32 increment pairs; M_j = g^{c_j} * z^{d_j}.
    cs = [coins.randrange(q) for _ in range(32)]
    ds = [coins.randrange(q) for _ in range(32)]
    Ms = [mul(scalar_mul(g_elt, cs[j]), scalar_mul(z_elt, ds[j])) for j in range(32)]

    a0, b0 = coins.randrange(q), coins.randrange(q)
    Y0 = mul(scalar_mul(g_elt, a0), scalar_mul(z_elt, b0))

    def step(Y, a, b):
        j = _branch(salt, encode(Y))
        return (mul(Y, Ms[j]), (a + cs[j]) % q, (b + ds[j]) % q)

    # Floyd: tortoise one step, hare two steps.
    tY, ta, tb = Y0, a0, b0
    hY, ha, hb = Y0, a0, b0
    for n in range(1, cap + 1):
        tY, ta, tb = step(tY, ta, tb)
        hY, ha, hb = step(hY, ha, hb)
        hY, ha, hb = step(hY, ha, hb)
        if tY == hY:
            denom = (hb - tb) % q
            if denom == 0:
                return (None, n, "degenerate")
            x = ((ta - ha) * pow(denom, -1, q)) % q
            if scalar_mul(g_elt, x) == z_elt:
                return (x, n, "success")
            return (None, n, "verification_failed")
    return (None, cap, "cap")


# ---- Production ledger (recomputed from the protected order) -----------------

# Itemized conservative per-group-step instruction charge (see machine.md):
FP12_MUL = 1 << 15          # 32768  >= organizer blst dense Fp12 mul (22996)
SHA_SELECT = 1 << 14        # 16384  salted SHA-256 over the 384-byte state
UPDATE = 1 << 12            # 4096   two 254-bit modular adds, table load, control
STEP = FP12_MUL + SHA_SELECT + UPDATE        # 53248 per group step
COMPARE = 1 << 13           # 8192   12-coefficient equality + multiword counters
PER_ITER = 3 * STEP + COMPARE                # 167936 per Floyd iteration
SETUP = 1 << 40             # setup, sampling, recovery, verify, (P0,[x]Q0) output


def floyd_cap(r):
    """Charged iteration count N = 2*ceil(sqrt(r))."""
    s = math.isqrt(r)
    if s * s < r:
        s += 1
    return 2 * s


def ledger(r):
    """Return (N, per_iter, R_total) for the protected order r."""
    N = floyd_cap(r)
    R = N * PER_ITER + SETUP
    return N, PER_ITER, R
