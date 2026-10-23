#!/usr/bin/env python3
"""Captured, exact small-curve doubled-add/Floyd/DP mechanism experiment."""
import math
import random
MASK = (1 << 64) - 1

def mix64(z):
    z &= MASK
    z = ((z ^ (z >> 30)) * 0xbf58476d1ce4e5b9) & MASK
    z = ((z ^ (z >> 27)) * 0x94d049bb133111eb) & MASK
    return z ^ (z >> 31)

def point_hash(point, key):
    x = 0 if point is None else point[0]
    for i in range(4):
        key = mix64(key ^ ((x >> (64 * i)) & MASK))
    return key

class Curve:
    def __init__(self, p, a, b):
        self.p, self.a, self.b = p, a, b
    def add(self, x, y):
        if x is None: return y
        if y is None: return x
        p = self.p
        if x[0] == y[0]:
            if (x[1] + y[1]) % p == 0: return None
            slope = (3*x[0]*x[0]+self.a)*pow(2*x[1],-1,p) % p
        else:
            slope = (y[1]-x[1])*pow((y[0]-x[0])%p,-1,p) % p
        u = (slope*slope-x[0]-y[0]) % p
        return u, (slope*(x[0]-u)-x[1]) % p
    def mul(self, a, point):
        result = None
        while a:
            if a & 1: result = self.add(result, point)
            point = self.add(point, point)
            a >>= 1
        return result
    def canon(self, point):
        if point is None or point[1] <= self.p//2: return point, 1
        return (point[0], self.p-point[1]), -1

def prime_small(n):
    return n >= 2 and all(n%d for d in range(2,math.isqrt(n)+1))

def run():
    p = 65519
    roots = {}
    for y in range(p): roots.setdefault(y*y%p,y)
    for b in range(1,100):
        if (4*(-3)**3 + 27*b*b)%p == 0: continue
        order, first = 1, None
        for x in range(p):
            rhs = (x*x*x-3*x+b)%p
            if rhs in roots:
                y = roots[rhs]
                order += 1 if y == 0 else 2
                if first is None and y: first = (x,y)
        if prime_small(order): break
    else: raise AssertionError('prime-order toy curve unavailable')
    curve = Curve(p,-3,b)
    assert prime_small(p) and curve.mul(order,first) is None
    points, point = [], None
    for _ in range(order):
        points.append(point)
        point = curve.add(point,first)
    assert point is None and len(set(points)) == order
    k = math.isqrt(1211*((order+1)//2)//1000)+2
    trials, walkers, table_size, dpbits = 400, 16, 256, 4
    rounds, tail = (k+walkers-1)//walkers, 24<<dpbits
    rng = random.Random(20261023)
    birthday = recovered = wrong = starts = restarts = events = 0
    meetings = maximum = 0
    for _ in range(trials):
        secret = rng.randrange(order)
        cs = [rng.randrange(order) for _ in range(table_size)]
        ds = [rng.randrange(order) for _ in range(table_size)]
        offsets = [(a+secret*d)%order for a,d in zip(cs,ds)]
        ik,dk = rng.getrandbits(64),rng.getrandbits(64)
        def canon(z,a,bb):
            _,s = curve.canon(points[z])
            return (z,a,bb) if s>0 else ((-z)%order,(-a)%order,(-bb)%order)
        def fresh():
            a,bb = rng.randrange(order),rng.randrange(order)
            return canon((a+secret*bb)%order,a,bb)
        def advance(state):
            z,a,bb = state
            j = point_hash(points[z],ik)&(table_size-1)
            out = canon((2*z+offsets[j])%order,(2*a+cs[j])%order,(2*bb+ds[j])%order)
            assert (out[1]+secret*out[2])%order == out[0]
            return out
        tort = [fresh() for _ in range(walkers)]
        hare,age = tort.copy(),[0]*walkers
        starts += walkers
        seen,store = {},{}
        collision = solved = bad = False
        nr = 0
        for tick in range(rounds+tail):
            for w in range(walkers):
                tort[w] = advance(tort[w])
                hare[w] = advance(advance(hare[w]))
                age[w] += 1
                z,a,bb = tort[w]
                if tick<rounds:
                    if z in seen and seen[z][0]!=w: collision=True
                    seen.setdefault(z,(w,a,bb))
                if point_hash(points[z],dk)&((1<<dpbits)-1)==0:
                    events += 1
                    if z in store:
                        aa,bbb = store[z]
                        den = (bb-bbb)%order
                        if den:
                            answer = (aa-a)*pow(den,-1,order)%order
                            good = points[answer]==points[secret]
                            recovered += int(good)
                            wrong += int(not good)
                            solved,bad = good,not good
                            break
                    else: store[z] = a,bb
                if tort[w][0]==hare[w][0]:
                    meetings += 1
                    maximum = max(maximum,age[w])
                    restarts += 1
                    nr += 1
                    if nr>walkers:
                        bad=True
                        break
                    tort[w] = fresh()
                    hare[w] = tort[w]
                    age[w] = 0
                    starts += 1
            if solved or bad: break
        birthday += int(collision)
    examples = []
    for mu,lam in ((0,1),(0,257),(100,193),(513,1021)):
        def f(i): return i+1 if i+1<mu+lam else mu
        t=h=0
        for meet in range(1,mu+lam+1):
            t,h = f(t),f(f(h))
            if t==h: break
        assert t==h and meet<=mu+lam
        examples.append({'mu':mu,'lambda':lam,'meeting_tick':meet})
    return {'seed':20261023,'curve':{'p':p,'a':-3,'b':b,'order':order,'generator':list(first),
            'order_method':'all affine x counted; trial-division primality'},'trials':trials,
            'walkers':walkers,'table_entries':table_size,'theta':'1/16','first_phase_rounds':rounds,
            'continuation_rounds':tail,'first_phase_cross_collisions_observed_before_stop':birthday,
            'recovered_and_verified':recovered,'wrong_verified_answers':wrong,
            'starts_including_restarts':starts,'floyd_restarts':restarts,'dp_reports':events,
            'floyd_meetings':meetings,'max_observed_meeting_tick':maximum,
            'long_cycle_floyd_examples':examples,
            'scope':'Mechanism and bounded toy pipeline only. Stop censors birthday count. Not production partition, probability, cost, or random-mapping certification.'}
