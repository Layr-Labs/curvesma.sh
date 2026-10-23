"""Exact affine EC arithmetic. Production pricing uses the explicit Fermat path.

The reference uses Python integers, not a claim that a Python operation costs one
instruction. The instruction ledger is in argument.md and certificate.json.
"""

def add(P, Q, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    X, Y = Q
    if x == X:
        if (y + Y) % p == 0:
            return None
        numerator, denominator = 3*x*x % p, 2*y % p
    else:
        numerator, denominator = (Y-y) % p, (X-x) % p
    slope = numerator * pow(denominator, p-2, p) % p
    xx = (slope*slope-x-X) % p
    return xx, (slope*(x-xx)-y) % p


def mul(k, P, p):
    Q = None
    for i in range(k.bit_length()-1, -1, -1):
        Q = add(Q, Q, p)
        if (k >> i) & 1:
            Q = add(Q, P, p)
    return Q


def canonical(P, c):
    """Return the 49-byte complete class key, canonical point and eigenvalue."""
    if P is None:
        return b'\x00' + bytes(48), None, 1
    p, r = c['p'], c['r']
    x, y = P
    bx = c['beta']*x % p
    xs = [x, bx, (-x-bx) % p]
    j = xs.index(min(xs))
    sign = 1 if 2*y < p else -1
    yy = sign*y % p
    u = sign * c['units'][j] % r
    return b'\x01' + yy.to_bytes(48, 'big'), (xs[j], yy), u


def toy_context(p, r, cofactor):
    beta = next(t for g in range(2, p) if (t := pow(g, (p-1)//3, p)) != 1)
    G = None
    for x in range(p):
        ys = [y for y in range(1, p) if (y*y-x*x*x-4) % p == 0]
        if ys:
            G = mul(cofactor, (x, ys[0]), p)
            if G is not None:
                break
    assert G is not None and mul(r, G, p) is None
    lam = next(k for k in range(2, r) if (k*k+k+1) % r == 0
               and mul(k, G, p) == (beta*G[0] % p, G[1]))
    return {'p': p, 'r': r, 'beta': beta, 'lam': lam, 'G': G,
            'units': (1, lam, lam*lam % r)}


def budget(r):
    from math import isqrt
    N = (r-1)//6
    B = isqrt(169*N//100)
    return B + (100*B*B < 169*N)


def solve(T, c, seed, check=False):
    """One bounded attempt. No restarts, table, checkpoints or unbounded sampling.

    The public reference expands a reproducible seed for initial toy coins.
    Its mathematical ideal model treats the keyed SHAKE values as independent.
    """
    from hashlib import shake_256
    p, r, G = c['p'], c['r'], c['G']
    B = budget(r)
    count = 0
    def done(x, reason, mu=None, period=None):
        if x is not None:
            assert 0 <= x < r and mul(x, G, p) == T
        return {'x': x, 'reason': reason, 'evaluations': count,
                'mu': mu, 'period': period, 'budget': B}
    if T is None:
        return done(0, 'identity target')
    assert len(seed) == 32
    raw = shake_256(b'initial-v3'+seed).digest(128)
    a, b = int.from_bytes(raw[:64], 'big') % r, int.from_bytes(raw[64:], 'big') % r
    if b == 0:
        return done(None, 'zero initial target coefficient')
    P = add(mul(a, G, p), mul(b, T, p), p)
    _, P, u = canonical(P, c)
    initial = (P, u*a % r, u*b % r)
    if P is None:
        return done(-a*pow(b, r-2, r) % r, 'initial identity relation')
    def key(state):
        Q = state[0]
        return bytes(49) if Q is None else b'\x01'+Q[1].to_bytes(48, 'big')
    def step(state):
        nonlocal count
        assert count < 6*B
        count += 1
        Q, aa, bb = state
        h = int.from_bytes(shake_256(b'orbit6-v3'+seed+key(state)).digest(64), 'big') % r
        Q = add(Q, mul(h, G, p), p)
        _, Q, u = canonical(Q, c)
        new = (Q, u*(aa+h) % r, u*bb % r)
        if check:
            assert add(mul(new[1], G, p), mul(new[2], T, p), p) == Q
        return new
    slow = fast = initial
    met = False
    for _ in range(B):
        slow = step(slow)
        fast = step(step(fast))
        if key(slow) == key(fast):
            met = True
            break
    if not met:
        return done(None, 'Floyd iteration cap')
    # Replay locates the entry by geometry; meeting-state labels are not used
    # as a relation because the gap may be several full laps of the cycle.
    left, right, mu = initial, slow, 0
    while key(left) != key(right):
        if mu >= B:
            return done(None, 'entry cap')
        left, right = step(left), step(right)
        mu += 1
    entry = left
    lap, period = step(entry), 1
    while key(lap) != key(entry):
        if period >= B:
            return done(None, 'period cap')
        lap = step(lap)
        period += 1
    denominator = (lap[2]-entry[2]) % r
    if denominator == 0:
        return done(None, 'one-lap orbit product is identity', mu, period)
    x = (entry[1]-lap[1])*pow(denominator, r-2, r) % r
    return done(x, 'verified one-lap relation', mu, period)


def wire_output(x, r):
    """The exact protected output object; never emit a point or a relation."""
    if not isinstance(x, int) or not 0 <= x < r:
        raise ValueError('noncanonical scalar')
    return {'x': format(x, '064x')}
