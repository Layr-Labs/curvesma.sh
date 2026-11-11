"""Executable mathematical reference.  Python is NOT the priced machine.

The priced scalar implementation is arithmetic.c; its executed-instruction
counts are derived in machine.md.  This file mirrors that code on integers.
"""
import math

MASK = (1 << 64) - 1
K1 = 0x9E3779B97F4A7C15
K2 = 0xD6E8FEB86659FD93
TABLE = 1024


# ---------------- word-level Montgomery reference (mirrors fp_mul) ---------
def mont_mul_words(a, b, m):
    aw = [(a >> (64 * i)) & MASK for i in range(4)]
    bw = [(b >> (64 * i)) & MASK for i in range(4)]
    mw = [(m >> (64 * i)) & MASK for i in range(4)]
    n0 = (-pow(mw[0], -1, 1 << 64)) & MASK
    t = [0] * 6
    for i in range(4):
        c = 0
        for j in range(4):
            v = aw[j] * bw[i] + t[j] + c
            assert v < 1 << 128
            t[j], c = v & MASK, v >> 64
        v = t[4] + c
        t[4], t[5] = v & MASK, v >> 64
        q = (t[0] * n0) & MASK
        v = q * mw[0] + t[0]
        assert v & MASK == 0
        c = v >> 64
        for j in range(1, 4):
            v = q * mw[j] + t[j] + c
            assert v < 1 << 128
            t[j - 1], c = v & MASK, v >> 64
        v = t[4] + c
        t[3], t[4] = v & MASK, t[5] + (v >> 64)
        assert t[4] <= 1
    u = sum(t[i] << (64 * i) for i in range(5))
    assert u < 2 * m
    return u - m if u >= m else u


# ---------------- tower arithmetic mirroring arithmetic.c ------------------
# Fp2 = (c0, c1), Fp6 = (Fp2, Fp2, Fp2), Fp12 = (Fp6, Fp6)
def fp2_mul(a, b, p):
    v0, v1 = a[0] * b[0] % p, a[1] * b[1] % p
    m = (a[0] + a[1]) * (b[0] + b[1]) % p
    return ((v0 - v1) % p, (m - v0 - v1) % p)


def fp2_add(a, b, p):
    return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)


def fp2_sub(a, b, p):
    return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)


def fp2_mul_xi(a, p):
    return ((9 * a[0] - a[1]) % p, (9 * a[1] + a[0]) % p)


def fp6_mul(a, b, p):
    v0, v1, v2 = fp2_mul(a[0], b[0], p), fp2_mul(a[1], b[1], p), fp2_mul(a[2], b[2], p)
    m = fp2_mul(fp2_add(a[1], a[2], p), fp2_add(b[1], b[2], p), p)
    c0 = fp2_add(v0, fp2_mul_xi(fp2_sub(fp2_sub(m, v1, p), v2, p), p), p)
    m = fp2_mul(fp2_add(a[0], a[1], p), fp2_add(b[0], b[1], p), p)
    c1 = fp2_add(fp2_sub(fp2_sub(m, v0, p), v1, p), fp2_mul_xi(v2, p), p)
    m = fp2_mul(fp2_add(a[0], a[2], p), fp2_add(b[0], b[2], p), p)
    c2 = fp2_add(fp2_sub(fp2_sub(m, v0, p), v2, p), v1, p)
    return (c0, c1, c2)


def fp6_add(a, b, p):
    return tuple(fp2_add(x, y, p) for x, y in zip(a, b))


def fp6_sub(a, b, p):
    return tuple(fp2_sub(x, y, p) for x, y in zip(a, b))


def fp12_mul(a, b, p):
    v0, v1 = fp6_mul(a[0], b[0], p), fp6_mul(a[1], b[1], p)
    m = fp6_mul(fp6_add(a[0], a[1], p), fp6_add(b[0], b[1], p), p)
    c1 = fp6_sub(fp6_sub(m, v0, p), v1, p)
    c0 = (fp2_add(v0[0], fp2_mul_xi(v1[2], p), p), fp2_add(v0[1], v1[0], p),
          fp2_add(v0[2], v1[1], p))
    return (c0, c1)


def tower_one():
    z = (0, 0)
    return (((1, 0), z, z), (z, z, z))


def fp12_pow(a, n, p):
    out = tower_one()
    while n:
        if n & 1:
            out = fp12_mul(out, a, p)
        a = fp12_mul(a, a, p)
        n >>= 1
    return out


def poly_to_tower(c, p):
    """Profile basis sum c_k w^k, w^12=18w^6-82, to the tower: w^6 = 9+u."""
    a = [((c[i] + 9 * c[i + 6]) % p, c[i + 6] % p) for i in range(6)]
    return ((a[0], a[2], a[4]), (a[1], a[3], a[5]))


def tower_to_poly(t, p):
    a = [t[0][0], t[1][0], t[0][1], t[1][1], t[0][2], t[1][2]]
    return tuple([(x - 9 * y) % p for x, y in a] + [y for _, y in a])


def poly_mul(a, b, p):
    c = [0] * 23
    for i in range(12):
        for j in range(12):
            c[i + j] = (c[i + j] + a[i] * b[j]) % p
    for k in range(22, 11, -1):
        c[k - 6] = (c[k - 6] + 18 * c[k]) % p
        c[k - 12] = (c[k - 12] - 82 * c[k]) % p
    return tuple(c[:12])


# ---------------- partition and distinguished-point predicate --------------
def mix(word_a, word_b, salt):
    h = ((word_a ^ salt) * K1) & MASK
    h ^= h >> 32
    h = (h + word_b) & MASK
    h = (h * K2) & MASK
    h ^= h >> 32
    return h


def bn_words(y, p):
    """The four words read by walk_step, from Montgomery limbs (R=2^256)."""
    mont = lambda v: (v << 256) % p
    wa = mont(y[0][0][0]) & MASK                 # c0.c0.c0.l[0]
    wb = (mont(y[1][1][1]) >> 64) & MASK         # c1.c1.c1.l[1]
    wc = mont(y[0][1][0]) & MASK                 # c0.c1.c0.l[0]
    wd = mont(y[1][0][1]) & MASK                 # c1.c0.c1.l[0]
    return wa, wb, wc, wd


def sample_scalar(r, coins):
    for _ in range(1024):
        v = coins.getrandbits(r.bit_length())
        if v < r:
            return v
    raise RuntimeError("bounded randomness rejection")


def cap_steps(r):
    """N = ceil(13*ceil(sqrt(r))/10)."""
    s = math.isqrt(r)
    s += s * s < r
    return -(-13 * s // 10)


def rho(r, target, identity, mul, power, words, dp_bits, coins, table_cap):
    """One bounded 1024-adding walk with distinguished points; no restarts.

    words(Y) -> (wa, wb, wc, wd); partition index = mix(wa,wb,salt)>>54;
    Y is distinguished iff the low dp_bits bits of wc | (wd << 64) vanish.
    Returns (x or None, steps, outcome).
    """
    if target == identity:
        return 0, 0, "identity"
    salt = coins.getrandbits(64)
    table = []
    for _ in range(TABLE):
        c, d = sample_scalar(r, coins), sample_scalar(r, coins)
        table.append((mul(power(c), power(d, target)), c, d))
    a, b = sample_scalar(r, coins), sample_scalar(r, coins)
    y = mul(power(a), power(b, target))
    dp_mask = (1 << dp_bits) - 1
    seen = []                      # append-only table, linear scan as priced
    index = {}                     # Python-only accelerator for the same scan
    cap = cap_steps(r)
    for k in range(1, cap + 1):
        wa, wb, wc, wd = words(y)
        m, c, d = table[mix(wa, wb, salt) >> 54]
        y, a, b = mul(y, m), (a + c) % r, (b + d) % r
        wa, wb, wc, wd = words(y)
        if (wc | (wd << 64)) & dp_mask:
            continue
        if y in index:
            a0, b0 = seen[index[y]][1:]
            den = (b - b0) % r
            if den == 0:
                return None, k, "degenerate"
            x = ((a0 - a) * pow(den, r - 2, r)) % r
            if power(x) != target:
                return None, k, "verification_failed"
            return x, k, "success"
        if len(seen) >= table_cap:
            return None, k, "table_full"
        index[y] = len(seen)
        seen.append((y, a, b))
    return None, cap, "cap"
