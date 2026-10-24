#!/usr/bin/env python3
"""make_scale_fixtures.py — random CM search (D=-3) for the same curve construction at
three class-space scales, for the p-transfer stability census.

Construction identical to evidence/census_v10f41.c (the 2^39.000003 fixture):
    p prime,  p = (t^2 + 3 v^2)/4,  E_b: y^2 = x^3 + b / F_p,
    #E = q = p + 1 - t prime (cofactor 1),  N' = (q+1)/2 x-classes,
    T = sqrt(2*Lambda*N') with Lambda fixed -> T/sqrt(N') preserved at every scale.

usage: make_scale_fixtures.py out.json
"""
import json
import math
import random
import sys

def is_prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, r = n-1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(r-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def ec_add(P, Q, p):
    if P is None: return Q
    if Q is None: return P
    if P[0] == Q[0]:
        if (P[1] + Q[1]) % p == 0: return None
        l = (3*P[0]*P[0]) * pow(2*P[1], -1, p) % p
    else:
        l = (Q[1]-P[1]) * pow(Q[0]-P[0], -1, p) % p
    x = (l*l - P[0] - Q[0]) % p
    return (x, (l*(P[0]-x) - P[1]) % p)

def ec_mul(k, P, p):
    R = None
    while k:
        if k & 1: R = ec_add(R, P, p)
        P = ec_add(P, P, p); k >>= 1
    return R

def tonelli(n, p):
    if p == 2: return n
    if pow(n, (p-1)//2, p) != 1: return None
    q, s = p-1, 0
    while q % 2 == 0: q //= 2; s += 1
    if s == 1: return pow(n, (p+1)//4, p)
    z = 2
    while pow(z, (p-1)//2, p) != p-1: z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q+1)//2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1: t2 = t2*t2 % p; i += 1
        bb = pow(c, 1 << (m-i-1), p)
        m, c, t, r = i, bb*bb % p, t*bb*bb % p, r*bb % p
    return r

def gen(k, tries=4_000_000, rng=None):
    """p ~ 2^(k+1), N' = (q+1)/2 ~ 2^k."""
    if rng is None:
        rng = random.Random(1234+k)
    tgt = 1 << (k+1)
    for _ in range(tries):
        # t^2 + 3v^2 = 4p ~ 2^(k+3); pick v uniform in a band, t from the circle equation
        v = rng.randrange(1, 1 << ((k+3)//2))
        rem = 4*tgt - 3*v*v
        if rem <= 0: continue
        t0 = math.isqrt(rem)
        t = t0 - rng.randrange(0, max(1, t0 >> 6))
        if t < 2: continue
        if (t*t + 3*v*v) % 4: continue
        p = (t*t + 3*v*v)//4
        if not (0.94*tgt < p < 1.06*tgt): continue
        if p % 3 != 1 or not is_prime(p): continue
        q = p + 1 - t
        if not is_prime(q): continue
        if abs(t) > 2*math.isqrt(p): continue
        N = (q+1)//2
        if not (0.9*(1 << k) < N < 1.1*(1 << k)): continue
        # find b and a generator on the order-q sextic twist
        for b in (2,3,5,6,7,10,11,13,14,15,17,19):
            for gx in range(1, 60):
                rhs = (gx*gx*gx + b) % p
                gy = tonelli(rhs, p)
                if gy is None or gy*gy % p != rhs: continue
                if ec_mul(q, (gx, gy), p) is None:
                    return dict(p=p, q=q, t=t, v=v, b=b, gx=gx, gy=gy, classes=N,
                                classes_log2=math.log2(N), p_log2=math.log2(p),
                                target_log2=k)
    return None

if __name__ == "__main__":
    outp = sys.argv[1] if len(sys.argv) > 1 else "fixtures_scales.json"
    out = {}
    for k in (25, 32):
        r = gen(k)
        if not r:
            print(f"scale 2^{k}: NOT FOUND", file=sys.stderr); continue
        out[k] = r
        print(f"2^{k}: p={r['p']} q={r['q']} t={r['t']} v={r['v']} b={r['b']} "
              f"G=({r['gx']},{r['gy']}) classes={r['classes']} log2={r['classes_log2']:.6f} "
              f"plog2={r['p_log2']:.6f}", flush=True)
    json.dump(out, open(outp, "w"), indent=1)
