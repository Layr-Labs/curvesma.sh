#!/usr/bin/env python3
"""The v12b class-step rule: CANONICALISING (no held-offset correction).

Why the rule changed
--------------------
v10f held an arbitrary representative y = sigma^(-k)(r) of its class and made the
step a class function by applying sigma^(-k)(m) instead of m, i.e. the table index
was i = (tw - k) mod 12 with k = kappa(y).  That is correct as a class function --
and it is exactly what created the rejected two-cycle: sigma^6 is inversion
(p^6 = -1 mod r), so sigma^(-(k+6))(m) = sigma^(-k)(m^-1).  As k spans Z/12 the
applied multiplier set {sigma^(-k)(m)} is closed under inversion for EVERY table,
so a step of m followed by a step of sigma^(-6)(m) returns the element, and the
label sum vanishes (delta_i + delta_{i+6} = 0).  The walk is absorbed with no
distinguished point emitted.  Restricting the routing field's modulus does not
help: k still spans Z/12, so the derived index still spans all twelve residues.

The rule
--------
Keep the state at the canononical representative and renormalise the labels:

    y  <-  sigma^k( y * m )          m = tbl[t][col],   k = kappa(y * m)
    (a, b) <- ( (a + delta) p^k , b p^k ) mod r

Now the multiplier applied to the element is always the table entry itself, never
sigma^(-k) of one, so the entire (m, m^-1) family is out of reach as soon as no
table entry's inverse is an entry.  That is a table property, made deterministic
by taking independent lower-half exponents 1 <= u <= (r-1)/2 (see build_table).

The squaring step renormalises the same way: y <- sigma^k(y^2), (a,b) ->
(a p^k, b p^k) doubled, i.e. (2a p^k, 2b p^k).

Solve
-----
Every stored record is at the canonical representative of its class, so two
records of one class have equal elements and labels (a1,b1), (a2,b2) satisfying
a1 + b1 x = a2 + b2 x; hence x = (a1 - a2) / (b2 - b1) mod r, with no Frobenius
k-search at all.  (v10f needed the k-search precisely because its records could
carry any of the twelve representatives.)
"""


def renorm_index(k, p_pow, r):
    """p^k mod r: the label renormalisation factor for a step ending at offset k."""
    return p_pow[k % 12] % r


def solve_x(a1, b1, a2, b2, r):
    """x = (a1 - a2)/(b2 - b1) mod r for two canonical records of one class.

    Returns None when the denominator vanishes, which is the (a1,b1)==(a2,b2)
    case: the same state, not a solvable pair.
    """
    den = (b2 - b1) % r
    if den == 0:
        return None
    return ((a1 - a2) % r) * pow(den, -1, r) % r


def table_is_inverse_free(us, r):
    """Deterministic check: no entry's exponent is the negation of another's.

    With every u in [1, (r-1)/2] the inverse exponents r-u lie in
    [(r+1)/2, r-1], which is disjoint from the used range.  Checked on the actual
    list, not argued from the range.
    """
    s = set(us)
    return len([u for u in us if (r - u) % r in s]) == 0


def class_step_desc():
    return ('y <- sigma^kappa(y*m)(y*m); (a,b) <- ((a+delta) p^kappa, b p^kappa); '
            'table exponents in [1,(r-1)/2] so no entry inverse is an entry; '
            'solve x = (a1-a2)/(b2-b1) mod r with no k-search')
