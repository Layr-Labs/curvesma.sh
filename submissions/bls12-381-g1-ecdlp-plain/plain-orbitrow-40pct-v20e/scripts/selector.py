"""selector.py -- v20: ONE 12-product window + one fixed-point division + subtraction.

v19d's selector called TWO 12-product windows (c1 and c2) while its ledger priced
one, and priced neither the (x<<128)//p division nor the near-tie fallback.  This
file is the actual implementation, and `ledger.py` prices exactly these functions:

    v1 = top128 of the low 448-bit window of x*floor(2^448 beta/p)
       = 2^128 * frac(beta*x/p) to relative 2^-64          (12 partial products)
    v0 = (x*floor(2^512/p)) >> 384
       = 2^128 * frac(x/p) to < 1 unit                      (3 partial products)
    v2 = (1<<128) - v0 - v1  mod 2^128
       = 2^128 * frac(beta^2*x/p), because frac(x/p)+frac(beta x/p)
         +frac(beta^2 x/p) is an integer for x in F_p^*      (2 subtracts)

h = argmin(v0, v1, v2); a top-2 gap below BAND is a near-tie and takes the exact
fallback, which is METERED (383 instructions) and measured in the toy.
"""
MASKL = (1 << 64) - 1
M128 = (1 << 128) - 1
BAND = 1 << 66            # 4x the measured 2^64.34 window error


def limbs(v, n):
    return [(v >> (64 * i)) & MASKL for i in range(n)]


def spec_window12(x, c):
    """VERBATIM production window: 12 partial products, (i,j) with i+j in {5,6}."""
    xl, cl = limbs(x, 6), limbs(c, 7)
    s5 = s6 = 0
    n = 0
    for i in range(6):
        for j in range(7):
            if i + j == 5:
                s5 += xl[i] * cl[j]
                n += 1
            elif i + j == 6:
                s6 += xl[i] * cl[j]
                n += 1
    return (((s6 + (s5 >> 64)) & MASKL) << 64) | (s5 & MASKL), n


def spec_recip3(x, c):
    """v0: 3 partial products, (i,j) with i+j in {6,7} of x*floor(2^512/p)."""
    xl, cl = limbs(x, 6), limbs(c, 3)
    s6 = s7 = 0
    n = 0
    for i in range(6):
        for j in range(3):
            if i + j == 6:
                s6 += xl[i] * cl[j]
                n += 1
            elif i + j == 7:
                s7 += xl[i] * cl[j]
                n += 1
    return (((s7 + (s6 >> 64)) & MASKL) << 64) | (s6 & MASKL), n


def window_generic(x, c, n_limbs, band):
    xl, cl = limbs(x, n_limbs), limbs(c, n_limbs + 1)
    acc = {}
    for i in range(n_limbs):
        for j in range(n_limbs + 1):
            k = i + j
            if k in band:
                acc[k] = acc.get(k, 0) + xl[i] * cl[j]
    top = max(band)
    s_top = acc.get(top, 0)
    s_below = acc.get(top - 1, 0)
    return (((s_top + (s_below >> 64)) & MASKL) << 64) | (s_below & MASKL)


def twist_index(x, p, beta, shift=None, n_limbs=None, band=None):
    """Exact-integer toy instantiation (full windows, no truncation)."""
    if n_limbs is None:
        n_limbs = max(1, (p.bit_length() + 63) // 64)
    if shift is None:
        shift = 64 * n_limbs
    if band is None:
        band = {n_limbs - 1, n_limbs}
    M = 1 << shift
    c1 = (M * beta) // p
    beta2 = (beta * beta) % p
    c2 = (M * beta2) // p
    v0 = (x << 128) // p
    v1 = window_generic(x, c1, n_limbs, band)
    v2 = window_generic(x, c2, n_limbs, band)
    vals = [v0, v1, v2]
    return vals.index(min(vals)), vals


def twist_index_exact(x, p, beta):
    beta2 = (beta * beta) % p
    vals = [x, (beta * x) % p, (beta2 * x) % p]
    return vals.index(min(vals)), vals


def twist_index_production(x, p, beta, band=None):
    """v20 production selector: 1 window (12 products) + v0 (3 products) + 2 subs.

    Returns (h, vals, near_tie).  On near_tie the caller MUST take the exact
    fallback (383 instructions, metered).
    """
    band = BAND if band is None else band
    M = 1 << 448
    c1 = (M * beta) // p
    c0 = (1 << 512) // p
    v1, _ = spec_window12(x, c1)
    v0, _ = spec_recip3(x, c0)
    v2 = ((1 << 128) - v0 - v1) & M128
    vals = [v0, v1, v2]
    ss = sorted(vals)
    return vals.index(ss[0]), vals, (ss[1] - ss[0]) < band


def twist_index_production_fallback(x, p, beta, band=None):
    """What the metered fallback returns: the exact argmin, plus the flag."""
    h, vals, nt = twist_index_production(x, p, beta, band)
    if nt:
        return twist_index_exact(x, p, beta)[0], True
    return h, False


def cell(K, R, P):
    """Table cell: (K*R) >> (bitlen(P)-1), exact shift (v17's map, unchanged)."""
    if not (0 <= K <= P // 2):
        raise ValueError("class key out of range")
    return (K * R) >> (P.bit_length() - 1)


def reachable_columns(R, P):
    return ((P // 2) * R >> (P.bit_length() - 1)) + 1


def window_error(p, beta, samples=4000):
    """max |12-product v1 - exact 2^128 frac(beta x/p)| (the near-tie band basis)."""
    M = 1 << 448
    c1 = (M * beta) // p
    worst = 0
    x = 7
    for _ in range(samples):
        x = (x * 6364136223846793005 + 1442695040888963407) % p
        if x == 0:
            continue
        v1, _ = spec_window12(x, c1)
        exact = (((beta * x) % p) << 128) // p
        worst = max(worst, abs(v1 - exact))
    return worst
