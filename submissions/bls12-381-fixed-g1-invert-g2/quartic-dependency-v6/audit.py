"""Only the shared verifier's credential-free isolated runner may execute this file.
Scope: arbitrary full-production field coefficients and exact work/success bounds.
This does not execute native assembly, invert a production instance, or test H/storage.
"""
from fractions import Fraction as F
from math import comb, factorial
import json
import random
import sys

P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
R = int('73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001', 16)
XI = (1, 1)
ZERO = (0, 0)

def add(x, y):
    return ((x[0]+y[0]) % P, (x[1]+y[1]) % P)

def sub(x, y):
    return ((x[0]-y[0]) % P, (x[1]-y[1]) % P)

def mul(x, y):
    return ((x[0]*y[0]-x[1]*y[1]) % P,
            (x[0]*y[1]+x[1]*y[0]) % P)

def quartic(x, y):
    ac, bd = mul(x[0], y[0]), mul(x[1], y[1])
    cross = mul(add(x[0], x[1]), add(y[0], y[1]))
    return (add(ac, mul(XI, bd)), sub(sub(cross, ac), bd))

def projection(x, y):
    # x,y are increasing powers of w; t=w^3 and t^2=XI.
    u = [(x[i], x[i+3]) for i in range(3)]
    v = [(y[i], y[i+3]) for i in range(3)]
    a, b, c = quartic(u[0], v[0]), quartic(u[1], v[2]), quartic(u[2], v[1])
    return (add(a[0], mul(XI, add(b[1], c[1]))),
            add(a[1], add(b[0], c[0])))

def dense(x, y):
    # Independent schoolbook convolution in Fp2[w]/(w^6-XI).
    z = [ZERO for _ in range(11)]
    for i in range(6):
        for j in range(6):
            q = (x[i][0]*y[j][0]-x[i][1]*y[j][1],
                 x[i][0]*y[j][1]+x[i][1]*y[j][0])
            z[i+j] = ((z[i+j][0]+q[0]) % P, (z[i+j][1]+q[1]) % P)
    for i in range(10, 5, -1):
        a, b = z[i]
        z[i-6] = ((z[i-6][0]+a-b) % P, (z[i-6][1]+a+b) % P)
    return z[0], z[3]

def certificate():
    w, m = 1 << 92, (1 << 61)-1
    length = 65726586350650495195887117980710461107
    lam = F(comb(length,2)*12*(R-1),R*R)
    bstar = 144*(F(comb(length,2)+12*comb(length,3),R*R)
                 +F(96*comb(length,4),R*R*w)
                 +F(6*comb(length,3),R*w*w)+F(54*comb(length,4),R*w**4))
    epsilon = 2*bstar+lam/(w*w)
    assert epsilon < F(7,10**17)
    cap = comb(length,2)*(F(1,w*w)+F(12,R)+F(1 << 23,m**3))/(length//160+1)
    t6 = sum((-lam)**k/factorial(k) for k in range(7))
    success = 1-t6-epsilon-lam/R-cap-F(length,R)-F(1,1 << 970)
    assert success > F(39001,100000)
    total = 27840*length+(1 << 124)
    assert total**200 < 20000**200*(1 << 25221)
    assert 9*1280+52*40+256 == 13856 < 13952
    return {'success_inequality': True, 'work_inequality': True,
            'success_diagnostic': float(success), 'score_cap': '25221/200'}

request = json.loads(sys.stdin.read())
count, seed = request['count'], request['seed']
assert isinstance(count, int) and 1 <= count <= 256
assert isinstance(seed, int) and 0 <= seed <= 4294967295
rng = random.Random(seed)
for i in range(count):
    x = [(rng.randrange(P), rng.randrange(P)) for _ in range(6)]
    y = [(rng.randrange(P), rng.randrange(P)) for _ in range(6)]
    if i == 0:
        x, y = [ZERO]*6, [ZERO]*6
    if i == 1:
        x, y = [(P-1,P-1)]*6, [(P-1,0)]*6
    assert projection(x,y) == dense(x,y)
print(json.dumps({'status':'passed','field_projection_cases':count,
                  'certificate':certificate(),
                  'scope':'full-width coefficient formulas only; no native execution, production inversion, hash or storage certification'}, sort_keys=True))
