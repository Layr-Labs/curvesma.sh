"""Toy Cheon p-1 BSGS on (Z/rZ, +) as a stand-in cyclic group. Not G1 cost evidence."""
import json, math, random, sys

def is_prime(n):
    if n < 2: return False
    for q in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % q == 0: return n == q
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,3,5,7,11,13,17,19,23):
        x = pow(a, d, n)
        if x in (1, n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def primitive_root(r, rng):
    # r prime => need generator of F_r^*
    fac = []
    m = r-1
    # trial factor for small toy r
    x=2
    while x*x <= m:
        if m%x==0:
            fac.append(x)
            while m%x==0: m//=x
        x+=1
    if m>1: fac.append(m)
    while True:
        g = rng.randrange(2, r)
        if all(pow(g, (r-1)//q, r) != 1 for q in fac):
            return g

def cheon(r, d, s, rng):
    # additive group: "g"=1, "g^a"=a mod r; scalar mul = field mul
    # We work directly in F_r for the DL of the secret, using Cheon on (Z/r)*
    # For the toy we just recover s from (s, s^d) via Cheon in F_r^* structure,
    # mirroring the two BSGS phases on exponents.
    xi = primitive_root(r, rng)
    # phase1: k2 from s^d = xi^{d k2}
    target = pow(s, d, r)
    zeta = pow(xi, d, r)
    m1 = math.isqrt((r-1)//d) + 1
    table = {}
    cur = target
    inv_zeta = pow(zeta, -1, r)
    for v in range(m1):
        table[cur] = v
        cur = cur * inv_zeta % r
    step = pow(zeta, m1, r)
    cur = 1
    k2 = None
    for u in range(m1):
        if cur in table:
            k2 = (m1*u + table[cur]) % ((r-1)//d)
            break
        cur = cur * step % r
    assert k2 is not None
    # phase2
    # s = xi^{((r-1)/d)*k1 + k2} => s * xi^{-k2} = xi^{((r-1)/d)*k1}
    t2 = s * pow(xi, -k2, r) % r
    delta = pow(xi, (r-1)//d, r)
    m2 = math.isqrt(d) + 1
    table = {}
    cur = t2
    inv_delta = pow(delta, -1, r)
    for v in range(m2):
        table[cur] = v
        cur = cur * inv_delta % r
    step = pow(delta, m2, r)
    cur = 1
    k1 = None
    for u in range(m2):
        if cur in table:
            k1 = (m2*u + table[cur]) % d
            break
        cur = cur * step % r
    assert k1 is not None
    rec = pow(xi, ((r-1)//d)*k1 + k2, r)
    return rec

def run(bits, trials, seed):
    rng = random.Random(seed)
    # find toy prime r with 2^k || r-1 for k>=8, take d=2^8
    while True:
        r = rng.getrandbits(bits) | (1<<(bits-1)) | 1
        if not is_prime(r): continue
        v=0; m=r-1
        while m%2==0: m//=2; v+=1
        if v >= 8: break
    d = 2**8
    ok = 0
    for _ in range(trials):
        s = rng.randrange(1, r)
        rec = cheon(r, d, s, rng)
        if rec == s: ok += 1
    return {"r": r, "d": d, "bits": bits, "trials": trials, "recovered": ok,
            "success_rate": ok/trials,
            "scope": "toy multiplicative group of prime field only; not G1 cost"}

if __name__ == '__main__':
    args = json.load(sys.stdin)
    print(json.dumps(run(args["bits"], args["trials"], args.get("seed", 1))))
