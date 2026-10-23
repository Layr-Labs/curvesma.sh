"""Own bounded coordinate-faithful GLV subgroup oracle; stdlib only.
The oracle accelerates toy additions only. Every partition/DP uses actual x.
No production input or production-size search is accepted.
"""
import math
import random
from array import array
P, B, Q = 6780512759743, 1, 523969
G = (1942568123522, 36770172157)
N = (Q-1)//6

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

def add(a, b):
    if a is None: return b
    if b is None: return a
    x,y = a; X,Y = b
    if x == X and (y+Y) % P == 0: return None
    s = (3*x*x*pow(2*y,-1,P) if x == X else (Y-y)*pow((X-x)%P,-1,P)) % P
    z = (s*s-x-X) % P
    return z, (s*(x-z)-y) % P

def cube_root(mod):
    g = 2
    while pow(g,(mod-1)//3,mod) == 1: g += 1
    return pow(g,(mod-1)//3,mod)

class Group:
    def __init__(self):
        assert prime(P) and prime(Q) and (G[1]*G[1]-G[0]**3-B)%P == 0
        self.points: list[tuple[int,int] | None] = [None]*Q
        logs = {}
        cur = None
        for k in range(1,Q):
            cur = add(cur,G)
            assert cur is not None and cur not in logs
            self.points[k] = cur; logs[cur] = k
        assert add(cur,G) is None
        self.beta = cube_root(P)
        self.lam = logs[(self.beta*G[0]%P,G[1])]
        assert self.lam != 1 and pow(self.lam,3,Q) == 1
        units = (1,self.lam,self.lam*self.lam%Q)
        self.canon_log = array('I',[0])*Q
        self.unit = array('I',[0])*Q
        reps = {}
        for k in range(1,Q):
            x,y = self.points[k]
            vals = (x,x*self.beta%P,x*self.beta*self.beta%P)
            j = min(range(3),key=vals.__getitem__)
            u = units[j]
            if y&1: u = -u % Q
            ck = k*u % Q
            assert self.points[ck] == (vals[j],P-y if y&1 else y)
            self.canon_log[k] = ck; self.unit[k] = u
            reps[vals[j]] = ck
        self.xs = sorted(reps)
        assert len(self.xs) == N
        self.kl = [reps[x] for x in self.xs]
        ix = {k:i for i,k in enumerate(self.kl)}
        self.canon_index = array('i',[-1])*Q
        for k in range(1,Q): self.canon_index[k] = ix[self.canon_log[k]]
        self.double = [self.canon_index[2*k%Q] for k in self.kl]
        self.double_unit = [self.unit[2*k%Q] for k in self.kl]
        rng = random.Random(92017)
        for _ in range(2000):
            a,b = rng.randrange(Q),rng.randrange(Q)
            assert add(self.points[a],self.points[b]) == self.points[(a+b)%Q]
        self.validation = dict(subgroup_points_enumerated=Q-1,
            subgroup_order_verified=True,matching_glv_all_nonidentity=True,
            independent_group_addition_fixtures=2000,classes=N,
            p=P,b=B,order=Q,generator=G,beta=self.beta,lambda_scalar=self.lam)

    def transitions(self,width,seed,shift=14):
        rng = random.Random(seed)
        jumps = [rng.randrange(Q) for _ in range(width)]
        branch = [(x>>shift)&(width-1) for x in self.xs]
        hist = [branch.count(i) for i in range(width)]
        assert min(hist)>0
        dest,unit = [],[]
        for i,k in enumerate(self.kl):
            z = (k+jumps[branch[i]]) % Q
            dest.append(self.canon_index[z]); unit.append(self.unit[z])
        for i in rng.sample(range(N),1000):
            z = add(self.points[self.kl[i]],self.points[jumps[branch[i]]])
            zz = (self.kl[i]+jumps[branch[i]]) % Q
            assert z == self.points[zz]
        return jumps,branch,dest,unit,hist
