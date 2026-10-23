"""Own bounded coordinate-faithful GLV subgroup oracle at THREE quotient scales.
Stdlib only. The oracle accelerates toy additions only; every partition/DP
predicate uses actual x. No production input or production-size search.
Small-scale build searches CM candidates and accepts only checked prime-order
points with matching GLV actions. Each Group stores its own subgroup order.
"""
import math
import random
from array import array

Q = 523969  # default large toy subgroup order; instances keep their own Q

def prime(n):
    if n < 2: return False
    for p in (2,3,5,7,11,13,17,19,23,29,31,37):
        if n % p == 0: return n == p
    d, s = n-1, 0
    while d % 2 == 0: d //= 2; s += 1
    for a in (2,325,9375,28178,450775,9780504,1795265022):
        x = pow(a, d, n)
        if x in (1,n-1): continue
        for _ in range(s-1):
            x = x*x % n
            if x == n-1: break
        else: return False
    return True

def add(P, a, b):
    if a is None: return b
    if b is None: return a
    x,y = a; X,Y = b
    if x == X and (y+Y) % P == 0: return None
    s = (3*x*x*pow(2*y,-1,P) if x == X else (Y-y)*pow((X-x)%P,-1,P)) % P
    z = (s*s-x-X) % P
    return z, (s*(x-z)-y) % P

def find_curve(q, bits):
    """Prime p of `bits` bits, p = lam mod q (lam = fixed cube root of unity
    in F_q), p = 19 mod 24 (so B=1 works: 2 has order 3 in F_p, curve
    y^2=x^3+1 is supersingular-ish enough), point of exact order q found
    by try-and-increment."""
    lam = next(u for u in range(2,q) if u != 1 and pow(u,3,q) == 1)
    step = q*3
    start = ((1<<(bits-1)) - 19)//step*step
    p = start + ((19 - start) % step)
    while True:
        p += step if p >= (1<<(bits-1)) and p != start else 0
        p = start + ((19 - start) % step) if p == start else p
        break
    # walk candidates p0 + t*lcm(3q) congruent 19 mod 24
    t = 0
    while True:
        p = start + ((19 - start) % step) + t*step
        t += 1
        if p < (1<<(bits-1)): continue
        if p >= (1<<bits): return None
        if not prime(p): continue
        # find G of order q: random x, point on y^2=x^3+1, multiply by (p+1)/q? order divides p+1 for supersingular; brute-force scalar search instead
        rng = random.Random(1000+bits)
        found = False
        for _ in range(4000):
            x = rng.randrange(1,p)
            rhs = (x*x*x+1) % p
            y = tonelli(rhs,p)
            if y is None: continue
            G = (x,y)
            h = (p+1)//q
            if (p+1) % q: break  # q must divide p+1 for cheap cofactor map
            g = scalar_mul(p, G, h)
            if g is not None and scalar_mul(p, g, q) is None:
                return p, 1, g, lam
        return p, 1, None, lam  # fallthrough caller retries

def tonelli(n, p):
    if n == 0: return 0
    if pow(n,(p-1)//2,p) != 1: return None
    if p % 4 == 3: return pow(n,(p+1)//4,p)
    q_, s = p-1, 0
    while q_ % 2 == 0: q_ //= 2; s += 1
    z = 2
    while pow(z,(p-1)//2,p) != p-1: z += 1
    m, c, t, r = s, pow(z,q_,p), pow(n,q_,p), pow(n,(q_+1)//2,p)
    while t != 1:
        i = 0; tt = t
        while tt != 1: tt = tt*tt % p; i += 1
        b = pow(c,1<<(m-i-1),p); m, c, t, r = i, b*b % p, t*b*b % p, r*b % p
    return r

def scalar_mul(P, G, k):
    R = None; A = G
    while k:
        if k & 1: R = add(P, R, A)
        A = add(P, A, A); k >>= 1
    return R

class Group:
    """One quotient scale. P: field prime; B: curve B; Q: order; G: generator."""
    def __init__(self, P, B, G, order=Q):
        Q = self.Q = order
        self.P, self.B, self.G = P, B, G
        assert prime(P) and prime(Q) and (G[1]*G[1]-G[0]**3-B) % P == 0
        assert scalar_mul(P, G, Q) is None and scalar_mul(P, G, Q//2) is not None
        N = (Q-1)//6
        self.N = N
        points = [None]*Q
        logs = {}
        cur = None
        for k in range(1,Q):
            cur = add(P, cur, G)
            assert cur is not None and cur not in logs
            points[k] = cur; logs[cur] = k
        assert add(P, cur, G) is None
        self.points = points
        beta = cube_root(P)
        lam = logs[(beta*G[0] % P, G[1])]
        assert lam != 1 and pow(lam,3,Q) == 1
        self.beta, self.lam = beta, lam
        units = (1, lam, lam*lam % Q)
        canon_log = array('I',[0])*Q
        unit = array('I',[0])*Q
        reps = {}
        for k in range(1,Q):
            x,y = points[k]
            vals = (x, x*beta % P, x*beta*beta % P)
            j = min(range(3), key=vals.__getitem__)
            u = units[j]
            if y & 1: u = -u % Q
            ck = k*u % Q
            assert points[ck] == (vals[j], P-y if y & 1 else y)
            canon_log[k] = ck; unit[k] = u
            reps[vals[j]] = ck
        self.canon_log, self.unit = canon_log, unit
        self.xs = sorted(reps)
        assert len(self.xs) == N
        self.kl = [reps[x] for x in self.xs]
        ix = {k:i for i,k in enumerate(self.kl)}
        canon_index = array('i',[-1])*Q
        for k in range(1,Q): canon_index[k] = ix[canon_log[k]]
        self.canon_index = canon_index
        self.double = [canon_index[2*k%Q] for k in self.kl]
        self.double_unit = [unit[2*k%Q] for k in self.kl]
        rng = random.Random(92017)
        for _ in range(2000):
            a,b = rng.randrange(Q), rng.randrange(Q)
            assert add(P, points[a], points[b]) == points[(a+b)%Q]
        self.validation = dict(subgroup_points_enumerated=Q-1,
            subgroup_order_verified=True, matching_glv_all_nonidentity=True,
            independent_group_addition_fixtures=2000, classes=N,
            p=P, b=B, order=Q, generator=G, beta=beta, lambda_scalar=lam)

    def transitions(self, width, seed, shift):
        Q = self.Q
        rng = random.Random(seed)
        jumps = [rng.randrange(Q) for _ in range(width)]
        branch = [(x >> shift) & (width-1) for x in self.xs]
        hist = [branch.count(i) for i in range(width)]
        assert min(hist) > 0
        dest, unit = [], []
        for i,k in enumerate(self.kl):
            z = (k + jumps[branch[i]]) % Q
            dest.append(self.canon_index[z]); unit.append(self.unit[z])
        for i in rng.sample(range(self.N), min(1000,self.N)):
            z = add(self.P, self.points[self.kl[i]], self.points[jumps[branch[i]]])
            zz = (self.kl[i] + jumps[branch[i]]) % Q
            assert z == self.points[zz]
        return jumps, branch, dest, unit, hist

def cube_root(mod):
    g = 2
    while pow(g,(mod-1)//3,mod) == 1: g += 1
    return pow(g,(mod-1)//3,mod)

def build(bits, seed=7, order=None):
    """Search checked j=0 curves, using CM candidates p=a^2+3b^2.

    (a-1)^2+3b^2=0 mod q makes q divide p+1-2a. This is only
    a candidate curve order, not an assertion about the B=1 twist: accept
    ONLY a nonzero point killed by q and the full Group GLV assertions.
    Never use p-1 or p+1 as an elliptic-curve order without checking.
    """
    q = Q if order is None else order
    assert prime(q) and q % 24 == 1
    root = tonelli((-3) % q, q)
    assert root is not None
    rng = random.Random(seed)
    limit = math.isqrt((1 << bits)-1)
    for _ in range(100000):
        b = rng.randrange(1, max(2,limit//2))
        residue = (1+root*b) % q
        a = residue + q*rng.randrange(max(1,(limit-residue)//q+1))
        p = a*a+3*b*b
        if p.bit_length() != bits or p % 3 != 1 or not prime(p): continue
        candidate_order = p+1-2*a
        assert candidate_order % q == 0
        h = candidate_order//q
        for _ in range(12):
            x = rng.randrange(1,p)
            y = tonelli((x*x*x+1) % p,p)
            if y is None: continue
            g = scalar_mul(p,(x,y),h)
            if g is not None and scalar_mul(p,g,q) is None:
                return Group(p,1,g,order=q)
    raise RuntimeError('no fully checked GLV curve for this seed/bits/order')

# Production-scale v9 curve kept verbatim as the largest toy scale
P9, B9 = 6780512759743, 1
G9 = (1942568123522, 36770172157)

def scale_groups():
    return dict(
        large=Group(P9, B9, G9),
        medium=build(32, order=5449),
        # q=601 gives delta=0.64 even at s=1; q=1297 is the first
        # prime ==1 mod24 with delta<0.30. No assertion is weakened.
        small=build(24, order=1297),
    )

if __name__ == '__main__':
    import sys, json
    gs = scale_groups()
    for name,g in gs.items():
        print(json.dumps(dict(scale=name, p_bits=g.P.bit_length(),
            classes=g.N, lam=g.lam)), flush=True)
