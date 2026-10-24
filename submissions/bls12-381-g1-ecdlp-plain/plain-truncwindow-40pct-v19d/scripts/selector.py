"""selector.py -- v19b: the ONE implementation of the twist-index predicate and the cell map.

Production and toy share this module; ledger.py prices exactly these functions, and
scripts/toy.py calls the same `twist_index` that the production ledger prices.

Spec (argument.md section 2).  For a class representative S = (x, y) on
y^2 = x^3 + b with j = 0 and a cube root of unity beta mod p, the three twists
x, beta*x, beta^2*x are the x-coordinates of S, phi(S), phi^2(S).  The canonical
twist index is

    h = argmin_{k in {0,1,2}}  ( beta^k * x  mod p )

which is a CLASS FUNCTION: multiplying x by beta permutes the three residues
cyclically, so the minimum is the same element for all six maps of the orbit.

The production estimate of the residues is a 12-partial-product truncated
product, verbatim from the probe that validated it:

    c1 = floor(2^448 * beta / p),  c2 = floor(2^448 * beta^2 / p)
    t  = (x * c1) mod 2^448                      # low 448-bit window
    v1 = t >> 320                                # top 128 bits of that window

`v1` is a relative 2^-64 approximation of 2^128 * (beta*x mod p) / p, because the
low window of x*c1 is the Montgomery residue beta*x*2^-448: the product's low
limbs are kept and the 12 partial products (i+j in {5,6}) drop only carries from
limbs 0..4.  Ranking v1, v2 against (x << 128)//p therefore ranks the three
residues.  beta^2*x needs no second product: beta^2 = -1-beta gives
B = (2p - x - A) mod p, two six-limb subtractions (priced in the decision row).

`spec_window12` is the probe's function, kept bit-identical; `twist_index` is the
generic form used by both the ledger and the toy, and the certificate asserts
they agree on 100 production keys (TOY-SPEC FIDELITY).
"""
MASKL = (1 << 64) - 1


def limbs(v, n):
    return [(v >> (64 * i)) & MASKL for i in range(n)]


def spec_window12(x, c):
    """VERBATIM production window: 12 partial products, (i,j) with i+j in {5,6}.

    Returns the top 128 bits of the low 448-bit window of x*c, and the number of
    partial products used (always 12 for a 6-limb x and a 7-limb c).
    """
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


def window_generic(x, c, n_limbs, band):
    """The same rule with the limb band as a parameter (the toy instantiates this).

    n_limbs = number of 64-bit limbs of x; band = the set of i+j kept.  For the
    production constants (n_limbs=6, band={5,6}) this is bit-identical to
    spec_window12; verify_certificate.py asserts that on 100 keys.
    """
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
    """h = argmin of the three residues, decided from the truncated windows.

    Production: shift=448, n_limbs=6, band={5,6}  (12 partial products).
    Toy:        shift=2*64*n_limbs (the full product), band={n_limbs-1, n_limbs}.
    Exact integers throughout; no reciprocal, no correction step.
    """
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
    v1 = window_generic(x, c1, n_limbs, band)
    v2 = window_generic(x, c2, n_limbs, band)
    v0 = (x << 128) // p
    vals = [v0, v1, v2]
    return vals.index(min(vals)), vals


def twist_index_exact(x, p, beta):
    """Exact argmin, integer arithmetic: the reference the window is measured against."""
    beta2 = (beta * beta) % p
    vals = [x, (beta * x) % p, (beta2 * x) % p]
    return vals.index(min(vals)), vals


def twist_index_production(x, p, beta):
    """The production instantiation the ledger prices: 448-bit window, 12 products."""
    return twist_index(x, p, beta, shift=448, n_limbs=6, band={5, 6})


def cell(K, R, P):
    """Table cell: (K*R) >> (bitlen(P)-1), exact shift (v17's map, unchanged)."""
    if not (0 <= K <= P // 2):
        raise ValueError("class key out of range")
    return (K * R) >> (P.bit_length() - 1)


def reachable_columns(R, P):
    return ((P // 2) * R >> (P.bit_length() - 1)) + 1
