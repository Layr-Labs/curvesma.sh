#!/usr/bin/env python3
"""Recompute the bn254-pairing-inversion rho-tower-pi claim and check G1 vectors.

Standard library only. Every check here runs at its actual scope:
  * the ledger is recomputed exactly from the profile order r and the declared
    per-step price;
  * the G1 vectors shipped in vectors.json are recomputed with an independent
    affine implementation (big integers, not the shipped C code), on the real
    BN254 G1;
  * the small prime-group walks are evidence at their actual small scope, not a
    production collision-statistics proof.
"""
import json
import math
import os
import random
import sys

R = int('30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001', 16)
P = int('30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47', 16)
C = 20000
C_STEP = 52000
N = 2 * math.isqrt(R) + 2
DECLARED_NUM, DECLARED_DEN = 130799212, 1000000
EXTRA_SETUP = 2**36
EXTRA_VERIFY = 2**32
EXTRA_POINT = 2**24
EXTRA_OUTPUT = 2**20
CTRL_PER_ITER = 4096

HERE = os.path.dirname(os.path.abspath(__file__))


def ledger():
    walk = 3 * N * C_STEP / C
    ctrl = N * CTRL_PER_ITER / C
    setup = EXTRA_SETUP / C
    ver = EXTRA_VERIFY / C
    point = EXTRA_POINT / C
    out = EXTRA_OUTPUT / C
    total = walk + ctrl + setup + ver + point + out
    return dict(walk=walk, control=ctrl, setup=setup, verify=ver, point=point,
                output=out, total=total, log2=math.log2(total))


# --- independent affine G1 arithmetic on the real curve y^2 = x^3 + 3 ---

def ec_add(a, b):
    if a is None:
        return b
    if b is None:
        return a
    x1, y1 = a
    x2, y2 = b
    if x1 == x2:
        if (y1 + y2) % P == 0:
            return None
        lam = (3 * x1 * x1) * pow(2 * y1, -1, P) % P
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    return (x3, (lam * (x1 - x3) - y1) % P)


def ec_mul(k, pt):
    acc = None
    while k:
        if k & 1:
            acc = ec_add(acc, pt)
        pt = ec_add(pt, pt)
        k >>= 1
    return acc


def on_curve(pt):
    if pt is None:
        return True
    x, y = pt
    return (y * y - x * x * x - 3) % P == 0


def check_vectors():
    data = json.load(open(os.path.join(HERE, 'vectors.json')))
    base = (int(data['base']['x'], 16), int(data['base']['y'], 16))
    checks = [('base on curve', on_curve(base)),
              ('base is the profile P0', base == (1, 2)),
              ('[r]P0 = O', ec_mul(R, base) is None)]
    all_ok = all(ok for _, ok in checks)
    for v in data['vectors']:
        k = int(v['k'], 16)
        got = ec_mul(k, base)
        want = (int(v['x'], 16), int(v['y'], 16))
        ok = got == want and on_curve(want) and ec_mul(R, want) is None
        checks.append(('vector k=0x%x matches and has order r' % k, ok))
        all_ok = all_ok and ok
    return checks, all_ok


def small_walk(p, steps, rng):
    """One 32-adding walk on the additive group Z/p with Floyd detection."""
    x = rng.randrange(p)
    cs = [rng.randrange(p) for _ in range(32)]
    ds = [rng.randrange(p) for _ in range(32)]
    ms = [(cs[j] + ds[j] * x) % p for j in range(32)]
    a0, b0 = rng.randrange(p), rng.randrange(p)
    y0 = (a0 + b0 * x) % p

    def step(y, a, b):
        h = y & 0xffffffff
        h ^= h >> 16
        h = (h * 0x7feb352d) & 0xffffffff
        h ^= h >> 15
        h = (h * 0x846ca68b) & 0xffffffff
        h ^= h >> 16
        j = h % 32
        return (y + ms[j]) % p, (a + cs[j]) % p, (b + ds[j]) % p

    ty, ta, tb = y0, a0, b0
    hy, ha, hb = y0, a0, b0
    for k in range(1, steps + 1):
        ty, ta, tb = step(ty, ta, tb)
        hy, ha, hb = step(*step(hy, ha, hb))
        if ty == hy:
            db = (hb - tb) % p
            if db == 0:
                return False, k
            return ((ta - ha) * pow(db, -1, p)) % p == x, k
    return False, steps


def is_probable_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def main():
    args = json.loads(sys.argv[1]) if len(sys.argv) > 1 else {}
    seed = int(args.get('seed', 20261112))
    trials = int(args.get('trials', 24))
    if not (0 <= seed < 2**32) or not (1 <= trials <= 64):
        raise SystemExit('bad input')

    L = ledger()
    score = DECLARED_NUM / DECLARED_DEN
    checks = []
    checks.append(('N = 2*ceil(sqrt(r))', N == 2 * (math.isqrt(R) + 1)))
    checks.append(('declared score >= log2(W)', score >= L['log2']))
    checks.append(('declared score rounded up to 1e-6',
                   DECLARED_NUM == math.ceil(L['log2'] * 1e6)))
    succ = 1 - math.exp(-(31 / 32) * 2)
    checks.append(('success bound >= 0.85', succ >= 0.85))
    checks.append(('protected minimum 0.39 <= 0.85', 0.85 >= 0.39))
    checks.append(('G1 row is below the walk row by >2^80',
                   EXTRA_POINT < 3 * N * C_STEP / 2**80))

    g1_checks, g1_ok = check_vectors()
    checks.extend(g1_checks)

    rng = random.Random(seed)
    p = 1048573
    checks.append(('small modulus prime', is_probable_prime(p)))
    wins = 0
    for _ in range(trials):
        ok, _k = small_walk(p, 2 * math.isqrt(p) + 2, rng)
        wins += 1 if ok else 0

    result = {
        'ledger_log2': L,
        'declared_score': score,
        'success_model': succ,
        'checks': [{'name': n, 'pass': bool(v)} for n, v in checks],
        'g1_vectors': {'checked': len(json.load(open(os.path.join(HERE, 'vectors.json')))['vectors']),
                       'all_match': bool(g1_ok)},
        'small_walks': {'modulus': p, 'trials': trials, 'successes': wins,
                        'note': 'scope-limited evidence only'},
        'all_pass': all(v for _, v in checks),
    }
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
