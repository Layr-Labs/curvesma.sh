"""Executable mathematical reference. Python runtime is NOT the scored machine.

Single 32-adding walk with distinguished points (DPs). In the attack the base
is g=e(P,Q0) (protected py_ecc pairing, computed once), the target is z, all
walk elements are kept in the blst tower representation with canonical
Montgomery coefficients, and the output is Q=[x]Q0 for the recovered x.

key(y) returns (index_word, dp_word, fingerprint): in production these are
word 0 of canonical coefficient t0, the low 128 bits of t1 and the low 128 bits
of t2 (machine.md). Test coins are reproducible, not a success proof.
"""
import math


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def walk_cap(r, gap):
    """S = floor(3*(isqrt(r)+1)/2) + gap >= 1.5*sqrt(r) + gap."""
    return (3 * (math.isqrt(r) + 1)) // 2 + gap


def rho_dp(r, target, identity, mul, power, key, coins, dbits, gap,
           max_entries=1 << 21, max_verify=64):
    """Returns (x, steps, outcome)."""
    if target == identity:
        return 0, 0, "identity"
    table = []
    for _ in range(32):
        a, b = sample_scalar(r, coins), sample_scalar(r, coins)
        table.append((mul(power(a), power(b, target)), a, b))
    a, b = sample_scalar(r, coins), sample_scalar(r, coins)
    y = mul(power(a), power(b, target))
    dmask = (1 << dbits) - 1
    cap = walk_cap(r, gap)
    stored = {}
    verifications = 0
    kk = key(y)
    for k in range(1, cap + 1):
        m, c, d = table[kk[0] & 31]
        y = mul(y, m)
        a, b = (a + c) % r, (b + d) % r
        kk = key(y)
        if kk[1] & dmask:
            continue
        prev = stored.get(kk[2])
        if prev is None:
            if len(stored) >= max_entries:
                return None, k, "table_full"
            stored[kk[2]] = (a, b)
            continue
        A, B = prev
        den = (b - B) % r
        if den == 0:
            return None, k, "degenerate"
        x = ((A - a) * pow(den, r - 2, r)) % r
        verifications += 1
        if power(x) == target:
            return x, k, "success"
        if verifications >= max_verify:
            return None, k, "verification_failed"
        stored[kk[2]] = (a, b)
    return None, cap, "cap"


# ---- Protected polynomial basis Fp[w]/(w^12 - 2 w^6 + 2) ----
def fp12_mul(a, b, p):
    c = [0] * 23
    for i in range(12):
        for j in range(12):
            c[i + j] = (c[i + j] + a[i] * b[j]) % p
    for k in range(22, 11, -1):
        t = (2 * c[k]) % p
        c[k - 6] = (c[k - 6] + t) % p
        c[k - 12] = (c[k - 12] - t) % p
    return tuple(c[:12])


def fp12_pow(a, n, p, mul=None):
    mul = mul or (lambda x, y: fp12_mul(x, y, p))
    out = None
    while n:
        if n & 1:
            out = a if out is None else mul(out, a)
        a = mul(a, a)
        n >>= 1
    return out


# ---- blst tower: Fp2=Fp[u]/(u^2+1), Fp6=Fp2[v]/(v^3-(1+u)), Fp12=Fp6[w]/(w^2-v)
# Flat tower order: index 6*e + 2*k + i is coefficient c_i of the Fp2 value
# multiplying v^k w^e.  Since w^2=v and w^6=1+u, this is the same w as the
# protected basis, so u = w^6 - 1.
def tower_to_poly(t, p):
    out = [0] * 12
    for e in range(2):
        for k in range(3):
            c0, c1 = t[6 * e + 2 * k], t[6 * e + 2 * k + 1]
            n = 2 * k + e
            out[n] = (out[n] + c0 - c1) % p
            out[n + 6] = (out[n + 6] + c1) % p
    return tuple(out)


def poly_to_tower(f, p):
    out = [0] * 12
    for e in range(2):
        for k in range(3):
            n = 2 * k + e
            out[6 * e + 2 * k] = (f[n] + f[n + 6]) % p
            out[6 * e + 2 * k + 1] = f[n + 6] % p
    return tuple(out)


def _f2mul(a, b, p):
    return ((a[0] * b[0] - a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)


def _f2add(a, b, p):
    return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)


def _f2xi(a, p):  # multiply by 1+u
    return ((a[0] - a[1]) % p, (a[0] + a[1]) % p)


def _f6mul(a, b, p):
    m = [[_f2mul(a[i], b[j], p) for j in range(3)] for i in range(3)]
    c0 = _f2add(m[0][0], _f2xi(_f2add(m[1][2], m[2][1], p), p), p)
    c1 = _f2add(_f2add(m[0][1], m[1][0], p), _f2xi(m[2][2], p), p)
    c2 = _f2add(_f2add(m[0][2], m[1][1], p), m[2][0], p)
    return (c0, c1, c2)


def _f6add(a, b, p):
    return tuple(_f2add(a[i], b[i], p) for i in range(3))


def _f6v(a, p):  # multiply by v
    return (_f2xi(a[2], p), a[0], a[1])


def tower_mul(x, y, p):
    def split(t):
        return [tuple((t[6 * e + 2 * k], t[6 * e + 2 * k + 1]) for k in range(3))
                for e in range(2)]
    A0, A1 = split(x)
    B0, B1 = split(y)
    c0 = _f6add(_f6mul(A0, B0, p), _f6v(_f6mul(A1, B1, p), p), p)
    c1 = _f6add(_f6mul(A0, B1, p), _f6mul(A1, B0, p), p)
    return tuple(v for half in (c0, c1) for pair in half for v in pair)


def production_key(t, p):
    """Canonical Montgomery residues (R=2^384) of the 12 tower coefficients."""
    mont = [(c << 384) % p for c in t]
    mask128 = (1 << 128) - 1
    return (mont[0] & ((1 << 64) - 1), mont[1] & mask128, mont[2] & mask128)
