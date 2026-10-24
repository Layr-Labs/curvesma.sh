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
    return cell_vec(K, R, P)


def cell_vec(K, R, P):
    """The SAME map, array-valued (qwalk.col imports this: one implementation).

    m = (K*R) >> (bitlen(P)-1), the truncated-quotient column.  The scalar `cell`
    above is this function plus the range guard, so the toy's walk and the ledger
    cannot drift apart.
    """
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


# ---------------------------------------------------------------------------
# secp256k1 retarget (v16, product count corrected in v17): 4x64-bit limbs,
# 320-bit v1 window, 384-bit v0 window.  The bls versions took 6-limb x and
# 7-limb c1 with bands i+j in {5,6}; at 256 bits the same construction is 4-limb x
# and 5-limb c1 with bands i+j in {3,4}.
#
# THE PRODUCT COUNT IS 8, NOT 7 (verdict d717195a point 3, confirmed on the
# listing): the band i+j == 3 has the FOUR pairs (0,3),(1,2),(2,1),(3,0) and the
# band i+j == 4 has the FOUR pairs (0,4),(1,3),(2,2),(3,1) -- i is 0..3, so
# (4,0) does not exist and the second band is not short by one.  8 = 4 + 4 is
# what the loops below execute, what ledger.ROWS now prices (35 insn, not 31),
# and what the assertion n1 == PRODUCTS_V1 below enforces at import time.  v0 is
# 3 products of the 384-bit window at shift 256 (bands {4,5}), exact (error 0).
# ---------------------------------------------------------------------------
PRODUCTS_V1 = 8
PRODUCTS_V0 = 3
NL_SECP = 4


def spec_window_secp(x, c):
    """v1: 8 partial products (4 at i+j=3, 4 at i+j=4) of the 4x5 limb 320-bit
    window; limbs 3..4 of the product are the top 128 bits of
    frac(beta x/p)*2^320."""
    xl, cl = limbs(x, 4), limbs(c, 5)
    s3 = s4 = 0
    n = 0
    for i in range(4):
        for j in range(5):
            if i + j == 3:
                s3 += xl[i] * cl[j]
                n += 1
            elif i + j == 4:
                s4 += xl[i] * cl[j]
                n += 1
    return (((s4 + (s3 >> 64)) & MASKL) << 64) | (s3 & MASKL), n


def spec_recip3_secp(x, c):
    """v0: 3 partial products (i,j) with i+j in {4,5} of x*floor(2^384/p)."""
    xl, cl = limbs(x, 4), limbs(c, 3)
    s4 = s5 = 0
    n = 0
    for i in range(4):
        for j in range(3):
            if i + j == 4:
                s4 += xl[i] * cl[j]
                n += 1
            elif i + j == 5:
                s5 += xl[i] * cl[j]
                n += 1
    return (((s5 + (s4 >> 64)) & MASKL) << 64) | (s4 & MASKL), n


def twist_index_production_secp(x, p, beta, band=None):
    """secp256k1 production selector: 1 eight-product window + v0 (3) + 2 subs."""
    band = BAND if band is None else band
    M = 1 << 320
    c1 = (M * beta) // p
    c0 = (1 << 384) // p
    v1, n1 = spec_window_secp(x, c1)
    v0, n0 = spec_recip3_secp(x, c0)
    assert n1 == PRODUCTS_V1 and n0 == PRODUCTS_V0, (n1, n0)
    v2 = ((1 << 128) - v0 - v1) & M128
    vals = [v0, v1, v2]
    ss = sorted(vals)
    return vals.index(ss[0]), vals, (ss[1] - ss[0]) < band


def window_error_secp(p, beta, samples=4000):
    """max |8-product v1 - exact 2^128 frac(beta x/p)| and the same for v0."""
    M = 1 << 320
    c1 = (M * beta) // p
    c0 = (1 << 384) // p
    worst1 = worst0 = 0
    x = 7
    for _ in range(samples):
        x = (x * 6364136223846793005 + 1442695040888963407) % p
        if x == 0:
            continue
        v1, _ = spec_window_secp(x, c1)
        v0, _ = spec_recip3_secp(x, c0)
        worst1 = max(worst1, abs(v1 - (((beta * x) % p) << 128) // p))
        worst0 = max(worst0, abs(v0 - ((x << 128) // p)))
    return worst1, worst0
