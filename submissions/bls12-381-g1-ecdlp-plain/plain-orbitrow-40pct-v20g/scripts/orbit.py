"""orbit.py -- v20: the six-element orbit structure as executable code.

A point of E: y^2 = x^3 + b with j = 0 has automorphism group mu_6 = <-1, phi>,
phi(x, y) = (beta x, y), beta^3 = 1, phi^3 = 1.  On the group, phi acts as
[lam] with lam^3 = 1 (mod the order), so the six maps are

    sigma_{k,b}(P) = (-1)^b * phi^k(P) = [-+lam^k] P ,  k in {0,1,2}, b in {0,1}

THE CLASS.  class(P) = orbit under those six maps.  canon(P) is the unique orbit
element with the minimum x-residue and y <= p - y.  Every class function in this
package is a function of canon(P).

THE INDEX.  h(P) = argmin_{i in {0,1,2}} (beta^i x_P mod p).  The three residues
are permuted cyclically by phi, so the SET is orbit-invariant and

    P = sigma_{k,b}(canon(P)),   k = (-h) mod 3,   b = [y > p - y]

exactly (asserted in the toy on 1e6 points).  Likewise K = min(y, p - y) is
orbit-invariant, and m = cell(K, R) is the table column.

THE EQUIVARIANT TABLE.  One seed increment J0[m] per column.  The six slots are

    Tab[m][h][b] = sigma_{k(h), b} J0[m]        (k(h) = (-h) mod 3)

so that, for every sigma in mu_6,

    Tab[m(sigma P)][h(sigma P)][b(sigma P)] = sigma Tab[m(P)][h(P)][b(P)]

i.e. J(sigma P) = sigma J(P) BY CONSTRUCTION, hence
class(P + J(P)) == class(sigma P + J(sigma P)): the induced map on the
(r-1)/6 classes is well defined and two walkers that meet in a class stay
together for ever.
"""


def canon(P, p, beta, beta2):
    x, y = P
    r0, r1, r2 = x, beta * x % p, beta2 * x % p
    if r0 <= r1 and r0 <= r2:
        cx = r0
    elif r1 <= r2:
        cx = r1
    else:
        cx = r2
    cy = y if y <= p - y else p - y
    return (cx, cy)


def orbit_index(P, p, beta, beta2):
    """(k, b) with P = sigma_{k,b}(canon(P))."""
    x, y = P
    vals = [x, beta * x % p, beta2 * x % p]
    h = vals.index(min(vals))
    k = (-h) % 3
    b = 1 if y > p - y else 0
    return k, b


def sig_point(P, p, beta, beta2, k, b):
    """(-1)^b phi^k(P)."""
    x, y = P
    for _ in range(k):
        x = beta * x % p
    return (x, (p - y) % p) if b else (x, y)


def sig_scalar(s, Q, lam, k, b):
    """The mu_6 action on a scalar coefficient: s -> (-1)^b lam^k s mod Q."""
    v = s * pow(lam, k, Q) % Q
    return (Q - v) % Q if b else v


def table_slot(h, b):
    """The six table slots per column, ordered by (h, b)."""
    return 2 * h + b


def seed_scalar(seed_bytes, Q):
    return int.from_bytes(seed_bytes, "big") % Q


def build_column(Q, lam, T, a0, m, s_val, mul, neg):
    """Build the six orbit slots of one column.

    J0 = [s_val]T (a multiple of T: this is the one-coefficient update the ledger
    prices).  Slot (h, b) holds (sigma_{k,b} J0, sigma_{k,b} s_val) with
    k = (-h) mod 3, so that adding slot (h(P), b(P)) to P makes the successor's
    class a function of the class of P alone.
    """
    J0 = mul(s_val, T)
    slots = []
    for h in range(3):
        k = (-h) % 3
        for b in range(2):
            J = sig_point(J0, jglob["p"], jglob["beta"], jglob["beta2"], k, b) if J0 is not None else None
            slots.append((J, sig_scalar(s_val, Q, lam, k, b)))
    return slots


jglob = {}
