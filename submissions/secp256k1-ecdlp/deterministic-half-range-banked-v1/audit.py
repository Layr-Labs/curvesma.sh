"""Bounded supporting checks. Execute only in the isolated verifier runner."""
import bisect
import json
import sys

assert json.load(sys.stdin) == {"mode": "audit"}
p = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
r = 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
P = (0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
     0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8)
m, h, G, T, F0 = 2**127, 2**128, 2**23, 2**15, 2**104
total = F0 + (m+h)*G + m*T + 127*m*(6*T+4096) + h*(128*(T+512)+4096)
pre = F0 + m*G + m*T + 127*m*(6*T+4096)
assert total == F0 + 59215872*m
assert pre == F0 + 33910784*m
assert pre < total < 1000*2**143
assert 2**255 < r < 2**256
assert 2**135 + 2**108 + 2**32 < 2**136
assert 150*128+8192 < T
assert 256*(2*256+64)+512 < 2**18
assert 512*2**18+256*64+512 < 2**28
assert 16*2**18+64*256+4096+7*2**18+2**28//4096+8192 < G

def partial_bsgs(order, width, height, target):
    records = sorted((j % order, j) for j in range(width))
    value = target
    for i in range(height):
        at = bisect.bisect_left(records, (value, -1))
        if at < width and records[at][0] == value:
            answer = i*width+records[at][1]
            assert 0 <= answer < width*height < order
            assert answer % order == target
            return answer
        value = (value-width) % order
    return None

groups_checked = 0
targets_checked = 0
for order in (17, 31, 67, 127, 257):
    for width, height in ((2, 4), (3, 5), (4, 4)):
        if width*height >= order:
            continue
        for target in range(order):
            answer = partial_bsgs(order, width, height, target)
            assert answer == (target if target < width*height else None)
            targets_checked += 1
        groups_checked += 1

def affine_add(a, b):
    if a is None:
        return b
    if b is None:
        return a
    x, y = a
    u, v = b
    if x == u:
        if (y+v) % p == 0:
            return None
        slope = 3*x*x*pow(2*y, -1, p) % p
    else:
        slope = (v-y)*pow(u-x, -1, p) % p
    xx = (slope*slope-x-u) % p
    return xx, (slope*(x-xx)-y) % p

def double(a):
    X, Y, Z = a
    if Z == 0 or Y == 0:
        return (0, 1, 0)
    A, C1 = X*X % p, Y*Y % p
    C2 = C1*C1 % p
    D = 2*((X+C1)**2-A-C2) % p
    E = 3*A % p
    F = E*E % p
    xx = (F-2*D) % p
    return xx, (E*(D-xx)-8*C2) % p, 2*Y*Z % p

def mixed(a, b):
    if b is None:
        return a
    X, Y, Z = a
    u, v = b
    if Z == 0:
        return u, v, 1
    z2 = Z*Z % p
    U = u*z2 % p
    S = v*Z*z2 % p
    H, t = (U-X) % p, (S-Y) % p
    if H == 0:
        return double(a) if t == 0 else (0, 1, 0)
    HH = H*H % p
    J, V, q = 4*H*HH % p, 4*X*HH % p, 2*t % p
    xx = (q*q-J-2*V) % p
    return xx, (q*(V-xx)-2*Y*J) % p, ((Z+H)**2-z2-HH) % p

def normalize(block):
    prefixes = []
    acc = 1
    for _, _, z in block:
        prefixes.append(acc)
        acc = acc*(z if z else 1) % p
    inv = pow(acc, -1, p)
    answer = [None]*len(block)
    for i in range(len(block)-1, -1, -1):
        X, Y, Z = block[i]
        zi = inv*prefixes[i] % p
        inv = inv*(Z if Z else 1) % p
        if Z:
            z2 = zi*zi % p
            answer[i] = X*z2 % p, Y*z2*zi % p
    return answer

points = [None]
for _ in range(32):
    points.append(affine_add(points[-1], P))
checks = 0
all_jacobian = []
expected = []
for i, a in enumerate(points):
    z = i+2
    ja = (0, 1, 0) if a is None else (a[0]*z*z % p, a[1]*z*z*z % p, z)
    for b in (None, P, a, None if a is None else (a[0], -a[1] % p)):
        got = mixed(ja, b)
        want = affine_add(a, b)
        assert normalize([got])[0] == want
        all_jacobian.append(got)
        expected.append(want)
        checks += 1
assert normalize(all_jacobian) == expected
assert normalize([(0, 1, 0)]*8) == [None]*8

print(json.dumps({
    "exact_accounting": True,
    "claimed_time_log2": 143,
    "instruction_upper_bound": str(total),
    "preprocessing_instruction_upper_bound": str(pre),
    "covered_probability": {"numerator": str(2**255), "denominator": str(r)},
    "toy_rectangles_checked": groups_checked,
    "toy_targets_checked": targets_checked,
    "production_curve_identity_checks": checks,
    "batch_normalization_checks": 2,
    "scope": "finite identities and integer budgets only; no production attack or instruction measurement"
}, sort_keys=True))
