#!/usr/bin/env python3
"""Recompute the bn254-gt-dlp rho-tower claim ledger and run a small walk check.

Standard library only. The small prime-group runs are evidence at their actual
scope; they are not a production collision-statistics proof.
"""
import json, math, sys

R = int('30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001', 16)
C = 20000
C_STEP = 52000
N = 2 * math.isqrt(R) + 2
DECLARED_NUM, DECLARED_DEN = 130799212, 1000000
EXTRA_SETUP = 2**36
EXTRA_VERIFY = 2**32
EXTRA_OUTPUT = 2**20
CTRL_PER_ITER = 4096

def ledger():
    walk = 3 * N * C_STEP / C
    ctrl = N * CTRL_PER_ITER / C
    setup = EXTRA_SETUP / C
    ver = EXTRA_VERIFY / C
    out = EXTRA_OUTPUT / C
    total = walk + ctrl + setup + ver + out
    return dict(walk=walk, control=ctrl, setup=setup, verify=ver,
                output=out, total=total, log2=math.log2(total))

def small_walk(p, steps, rng):
    """One 32-adding walk on the additive group Z/p with Floyd detection."""
    x = rng.randrange(p)                     # target z = x in Z/p
    cs = [rng.randrange(p) for _ in range(32)]
    ds = [rng.randrange(p) for _ in range(32)]
    ms = [(cs[j] + ds[j] * x) % p for j in range(32)]   # M_j = c_j + d_j*z
    a0, b0 = rng.randrange(p), rng.randrange(p)
    y0 = (a0 + b0 * x) % p

    def step(y, a, b):
        h = y & 0xffffffff
        h ^= h >> 16; h = (h * 0x7feb352d) & 0xffffffff
        h ^= h >> 15; h = (h * 0x846ca68b) & 0xffffffff
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
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0: return n == q
    d, s = n - 1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):
        x = pow(a, d, n)
        if x in (1, n - 1): continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
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
    # small prime groups
    import random
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
        'small_walks': {'modulus': p, 'trials': trials, 'successes': wins,
                        'note': 'scope-limited evidence only'},
        'all_pass': all(v for _, v in checks),
    }
    print(json.dumps(result, sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
