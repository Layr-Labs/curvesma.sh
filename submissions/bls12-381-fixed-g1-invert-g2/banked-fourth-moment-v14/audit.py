"""Only the shared verifier's credential-free isolated runner may execute this file.
Scope: arbitrary full-production field coefficients and exact work/success bounds.
This does not execute native assembly, invert a production instance, or test H/storage.
"""
from fractions import Fraction as F
from math import comb
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

def grid_rows():
    phases = [pow(P,i,R) for i in range(12)]
    assert len(set(phases)) == 12 and phases[6] == R-1
    totals = []
    for width,height in [(2,2),(2,4),(3,3),(4,3)]:
        types = [(j,k) for k in range(height) for j in range(width)]
        seen = set()
        for s,(js,ks) in enumerate(types):
            for jt,kt in types[s+1:]:
                for i,lam in enumerate(phases):
                    if i == 0 and (js == jt or ks == kt):
                        continue
                    if i == 6 and js != jt and ks != kt and js > jt:
                        continue
                    row = [0]*(width+height)
                    row[jt] += 1
                    row[js] -= lam
                    row[width+kt] += 1
                    row[width+ks] -= lam
                    row = [v % R for v in row]
                    pivot = next(v for v in row if v)
                    inv = pow(pivot,-1,R)
                    canonical = tuple(v*inv % R for v in row)
                    assert canonical not in seen
                    seen.add(canonical)
        expected = 11*comb(width*height,2)+comb(width,2)*comb(height,2)
        assert len(seen) == expected
        totals.append({'width':width,'height':height,'distinct_projective_rows':len(seen)})
    return totals

def certificate():
    width,height,m = 1 << 92,13568405678,(1 << 61)-1
    length = width*height
    assert length == 67187490632253187137321911992595775488
    assert 23*width**2*10000**2*height**2 >= 2*R*9950**2
    assert 23*width**2*10000**2*(height-1)**2 < 2*R*9950**2
    events = 11*comb(length,2)+comb(width,2)*comb(height,2)
    degree = 24*length*(width+height)
    p0 = F(R-1,R*R)
    mu = events*p0
    a2 = F(comb(events,2),R*R)
    assert comb(events,3) > events**2*degree
    a3 = (comb(events,3)-events**2*degree)*p0**3
    k2,k3 = 12*comb(16,2),12*comb(36,2)
    a4 = F(comb(events,4),R**4)+F(3*events**2*degree*k3,R**3)+F(events*degree*k2*k2,R**2)
    cap = comb(length,2)*(F(12,R)+F(1 << 23,m**3))/(length//160+1)
    success = mu-a2+a3-a4-mu/R-cap-F(length,R)-F(1,1 << 970)
    assert success > F(39004,100000)
    total = 19200*length+(1 << 120)
    assert total**100 < 20000**100*(1 << 12561)
    assert 1024*(32*255+81*61) == 13415424 < (1 << 24)
    assert 9*1280+52*40+256 == 13856 < 13952
    return {'success_inequality':True,'work_inequality':True,
            'success_diagnostic':float(success),'score_cap':'12561/100',
            'rank2_error_diagnostic':float(F(events*degree*k2*k2,R**2))}

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
                  'certificate':certificate(),'small_grid_row_checks':grid_rows(),
                  'scope':'full-width coefficient formulas only; no native execution, production inversion, hash or storage certification'}, sort_keys=True))
