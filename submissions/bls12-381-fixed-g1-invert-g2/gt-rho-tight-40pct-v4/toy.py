#!/usr/bin/env python3
"""End-to-end toy: the production walk on toy towers of exactly the production shape.

Toy scale only changes three numbers: the field modulus, the class count N and the
distinguished-point density (2^-5 instead of 2^-80, because a toy trail must be
shorter than sqrt(N)).  The class rule, the fold, the table with twelve
Frobenius twists, the step mixture (routing.fields), the ordered store, the replay
with exponent bookkeeping, the Frobenius alignment, the exponent solve and the
verification against z are the production algorithm.

Memory is bounded: the store refuses new reports at a fixed capacity and sets a
loss flag, as the production cap does.
"""
import json, random, sys, time
from fractions import Fraction
import tower
import routing
import constants as K

TOY_DP_BITS = 7           # toy distinguished density 2^-7 (production: 2^-80)
STORE_CAP = 4096          # fixed-cap refusal
MAX_SEG = 4096            # hard segment cap
R_TABLE = 16


def primes(n):
    s = [True]*(n+1); s[0] = s[1] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]:
            for j in range(i*i, n+1, i):
                s[j] = False
    return [i for i in range(n+1) if s[i]]


def is_prime(n):
    return tower.is_prime(n)


SMALL = primes(1 << 16)


TOY_X = -56          # BLS12 parameter; the same parameterisation as BLS12-381


def find_toy():
    """A BLS12 tower of production shape: r = x^4-x^2+1, p = (x^6-2x^5+2x^3+x+1)/3.

    Every property the algorithm relies on is re-checked here, not assumed: p and r
    prime, p = 7 mod 12, ord_r(p) = 12, and a random final exponentiation landing in
    the order-r subgroup of Fp12*.
    """
    x = TOY_X
    r = x**4 - x**2 + 1
    p = (x**6 - 2*x**5 + 2*x**3 + x + 1)//3
    assert is_prime(p) and is_prime(r), (p, r)
    assert p % 12 == 7, p % 12
    assert (r - 1) % 12 == 0
    assert all(pow(p, k, r) != 1 for k in (1, 2, 3, 4, 6)) and pow(p, 12, r) == 1
    return p, r, (r-1)//12


def build_table(T, r, g, z, rng, R=R_TABLE):
    muls = []
    for _ in range(R):
        u = rng.randrange(1, r); v = rng.randrange(1, r)
        el = T.mul(T.pow(g, u), T.pow(z, v))
        twists = []
        for i in range(12):
            twists.append((el, (u*pow(T.p, i, r)) % r, (v*pow(T.p, i, r)) % r))
            el = T.frob(el)
        muls.append(twists)
    return muls


def key_of(T, y):
    return min(T.keys12(y))


def fold(T, y, r):
    """Canonical representative of the class and the Frobenius power used."""
    ks = T.keys12(y)
    k = ks.index(min(ks))
    for _ in range(k):
        y = T.frob(y)
    return y, k


def step(T, y, muls, r, a, b):
    """One production step.  Returns (y', a', b')."""
    k = key_of(T, y)
    tbl, tw, dp, st = routing.fields(k)
    if st < 3:
        y = T.cyclo_sqr_fast(y)
        a = (2*a) % r; b = (2*b) % r
        i = 0
    else:
        j = tbl % R_TABLE
        i = tw % 12
        el, uu, vv = muls[j][i]
        y = T.mul(y, el)
        a = (a + uu) % r; b = (b + vv) % r
    y, kap = fold(T, y, r)
    a = a*pow(T.p, kap, r) % r; b = b*pow(T.p, kap, r) % r
    return y, a, b


def run_group(T, r, g, z, muls, rng, cap=STORE_CAP, segments=4):
    """Independent trails, each ending at one report; an ordered store detects a repeat."""
    store = {}
    loss = False
    trails = []
    reports = 0
    for seg in range(segments):
        a0 = rng.randrange(1, r); b0 = rng.randrange(1, r)
        y = T.mul(T.pow(g, a0), T.pow(z, b0))
        trails.append((a0, b0, y))
        a, b = a0, b0
        for idx in range(MAX_SEG):
            y, a, b = step(T, y, muls, r, a, b)
            k = key_of(T, y)
            if (routing.fields(k)[2] & ((1 << TOY_DP_BITS)-1)) != 0:
                continue
            reports += 1
            if k in store:
                return dict(hit=True, pair=(store[k], (seg, idx+1)), loss=loss,
                            reports=reports, trails=trails)
            if len(store) >= cap:
                loss = True
                break
            store[k] = (seg, idx+1)
            break
    return dict(hit=False, loss=loss, reports=reports, trails=trails, size=len(store))


def replay(T, r, g, z, muls, trails, tid, steps):
    a0, b0, y = trails[tid]
    a, b = a0, b0
    for _ in range(steps):
        y, a, b = step(T, y, muls, r, a, b)
    return y, a, b


def group(T, r, g, z, muls, rng, secret, cap=STORE_CAP, segments=4):
    out = run_group(T, r, g, z, muls, rng, cap, segments)
    if not out['hit']:
        return dict(recovered=False, hit=False, loss=out['loss'])
    (t1, s1), (t2, s2) = out['pair']
    y1, a1, b1 = replay(T, r, g, z, muls, out['trails'], t1, s1)
    y2, a2, b2 = replay(T, r, g, z, muls, out['trails'], t2, s2)
    j = None
    for k in range(12):
        if y1 == T.pow(y2, pow(T.p, k, r)):
            j = k
            break
    if j is None:
        return dict(recovered=False, hit=True, aligned=False)
    pj = pow(T.p, j, r)
    num = (a2*pj - a1) % r
    den = (b1 - b2*pj) % r
    if den == 0:
        return dict(recovered=False, hit=True, aligned=True, singular=True)
    x = num*pow(den, r-2, r) % r
    ok = (T.pow(g, x) == z)
    return dict(recovered=True, hit=True, aligned=True, verified=ok, secret_match=(x == secret))


def main():
    t0 = time.time()
    p, r, N = find_toy()
    T = tower.Tower(p)
    rng = random.Random(424242)
    e = (p**12 - 1)//r
    g = None
    while g is None:
        x = tuple(tuple((rng.randrange(p), rng.randrange(p)) for _ in range(3)) for _ in range(2))
        h = T.pow(x, e)
        if h != T.one12() and T.pow(h, r) == T.one12():
            g = h
    groups, segments = 40, 16
    res = dict(groups=0, hits=0, recovered=0, verified=0, secret_match=0, cap_losses=0,
               aligned_fail=0, singular=0, p=p, r=r, N=N, dp_bits=TOY_DP_BITS,
               segments=segments, store_cap=STORE_CAP, toy_x=TOY_X)
    for gi in range(groups):
        secret = rng.randrange(1, r)
        z = T.pow(g, secret)
        muls = build_table(T, r, g, z, rng)
        o = group(T, r, g, z, muls, rng, secret, segments=segments)
        res['groups'] += 1
        if o.get('hit'):
            res['hits'] += 1
        if o.get('loss'):
            res['cap_losses'] += 1
        if o.get('recovered'):
            res['recovered'] += 1
            if o.get('verified'):
                res['verified'] += 1
            if o.get('secret_match'):
                res['secret_match'] += 1
        if o.get('hit') and o.get('aligned') is False:
            res['aligned_fail'] += 1
        if o.get('singular'):
            res['singular'] += 1
    # fixed-cap refusal: a store that is too small must refuse and set the loss flag
    secret = rng.randrange(1, r)
    z = T.pow(g, secret)
    muls = build_table(T, r, g, z, rng)
    small = run_group(T, r, g, z, muls, rng, cap=1, segments=segments)
    res['cap_probe_refused'] = bool(small['loss'])
    res['cap_probe_hit'] = bool(small['hit'])
    res['seconds'] = round(time.time()-t0, 1)
    res['all_verified'] = (res['recovered'] == res['verified'] == res['secret_match']
                           and res['recovered'] > 0 and res['cap_probe_refused'])
    json.dump(res, open("toy-out.json", "w"), indent=1)
    print("toy p=%d r=%d N=%d dp_bits=%d segments=%d" % (p, r, N, TOY_DP_BITS, segments))
    print("groups %d hits %d recovered %d verified %d secret_match %d cap_loss %d aligned_fail %d"
          % (res['groups'], res['hits'], res['recovered'], res['verified'],
             res['secret_match'], res['cap_losses'], res['aligned_fail']))
    print("cap probe: refused %s" % res['cap_probe_refused'])
    print("TOY", "PASS" if res['all_verified'] else "FAIL", "%.1fs" % res['seconds'])
    return 0 if res['all_verified'] else 1


if __name__ == "__main__":
    sys.exit(main())
