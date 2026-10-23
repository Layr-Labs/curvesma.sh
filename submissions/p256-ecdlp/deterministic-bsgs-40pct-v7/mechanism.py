"""Bounded mechanism checks, not an instruction measurement or full-size run.
All affine coordinates here are ordinary integers. Jacobian tests use CIOS
Montgomery field multiplication; the independent affine oracle uses pow().
"""
from dataclasses import dataclass
from typing import Any

WORD = 1 << 64
MASK = WORD - 1


class Field:
    def __init__(self, p):
        self.p = p
        self.mod = [p >> (64 * i) & MASK for i in range(4)]
        self.radix = 1 << 256
        self.radix_inv = pow(self.radix, -1, p)
        self.pprime = -pow(p, -1, WORD) % WORD
        self.calls = 0
        self.pairs = 0

    def cios(self, a, b):
        """Four CIOS rounds, exact nine-word accumulator, canonical output."""
        assert 0 <= a < self.p and 0 <= b < self.p
        aa = [a >> (64 * i) & MASK for i in range(4)]
        bb = [b >> (64 * i) & MASK for i in range(4)]
        t = [0] * 9

        def add_pair(k, x, y):
            # MUL/UMULH, ADDS/ADCS and a bounded carry ripple.
            product = x * y
            lo, hi = product & MASK, product >> 64
            s = t[k] + lo
            t[k], carry = s & MASK, s >> 64
            s = t[k + 1] + hi + carry
            t[k + 1], carry = s & MASK, s >> 64
            for pos in range(k + 2, 9):
                s = t[pos] + carry
                t[pos], carry = s & MASK, s >> 64
            assert carry == 0
            self.pairs += 1

        for i in range(4):
            for j in range(4):
                add_pair(j, aa[j], bb[i])
            q = t[0] * self.pprime & MASK
            for j in range(4):
                add_pair(j, q, self.mod[j])
            assert t[0] == 0
            t = t[1:] + [0]
        value = sum(x << (64 * i) for i, x in enumerate(t))
        assert value < 2 * self.p
        if value >= self.p:
            value -= self.p
        self.calls += 1
        return value

    def enc(self, a):
        return a * self.radix % self.p

    def dec(self, a):
        return self.cios(a, 1)

    def mul(self, a, b):
        return self.cios(a, b)

    def add(self, a, b):
        return (a + b) % self.p

    def sub(self, a, b):
        return (a - b) % self.p

    def inv(self, a):
        assert a != 0
        result = self.enc(1)
        for i in range(255, -1, -1):
            result = self.mul(result, result)
            if (self.p - 2) >> i & 1:
                result = self.mul(result, a)
        return result


def affine_add(a, b, p):
    if a is None:
        return b
    if b is None:
        return a
    x, y = a
    u, v = b
    if x == u:
        if (y + v) % p == 0:
            return None
        slope = (3 * x * x - 3) * pow(2 * y, -1, p) % p
    else:
        slope = (v - y) * pow((u - x) % p, -1, p) % p
    xx = (slope * slope - x - u) % p
    return xx, (slope * (x - xx) - y) % p


def affine_scalar(k, point, p):
    if k < 0:
        return affine_scalar(-k, (point[0], -point[1] % p), p)
    out = None
    while k:
        if k & 1:
            out = affine_add(out, point, p)
        point = affine_add(point, point, p)
        k >>= 1
    return out


def jac(point, f):
    return (f.enc(1), f.enc(1), 0) if point is None else (
        f.enc(point[0]), f.enc(point[1]), f.enc(1))


def twice(point, f):
    x, y, z = point
    if not z or not y:
        return jac(None, f)
    mul, add, sub = f.mul, f.add, f.sub
    delta, gamma = mul(z, z), mul(y, y)
    beta = mul(x, gamma)
    alpha = mul(sub(x, delta), add(x, delta))
    alpha = add(add(alpha, alpha), alpha)
    two_beta = add(beta, beta)
    four_beta = add(two_beta, two_beta)
    eight_beta = add(four_beta, four_beta)
    xx = sub(mul(alpha, alpha), eight_beta)
    zz = sub(sub(mul(add(y, z), add(y, z)), gamma), delta)
    gamma2 = mul(gamma, gamma)
    two_gamma2 = add(gamma2, gamma2)
    four_gamma2 = add(two_gamma2, two_gamma2)
    eight_gamma2 = add(four_gamma2, four_gamma2)
    yy = sub(mul(alpha, sub(four_beta, xx)), eight_gamma2)
    return xx, yy, zz


def mixed(point, addend, f):
    if addend is None:
        return point
    if not point[2]:
        return jac(addend, f)
    x, y, z = point
    u, v = f.enc(addend[0]), f.enc(addend[1])
    mul, add, sub = f.mul, f.add, f.sub
    z2 = mul(z, z)
    u2 = mul(u, z2)
    s2 = mul(mul(v, z), z2)
    h, d = sub(u2, x), sub(s2, y)
    if not h:
        return twice(point, f) if not d else jac(None, f)
    hh = mul(h, h)
    ii = add(add(hh, hh), add(hh, hh))
    jj = mul(h, ii)
    rr = add(d, d)
    vv = mul(x, ii)
    xx = sub(sub(mul(rr, rr), jj), add(vv, vv))
    yj = mul(y, jj)
    yy = sub(mul(rr, sub(vv, xx)), add(yj, yj))
    zz = sub(sub(mul(add(z, h), add(z, h)), z2), hh)
    return xx, yy, zz


def normalize(points, f):
    """3B multiplication batch inverse, 2M x + 3M y + 2M decoding."""
    one = f.enc(1)
    zz = [z if z else one for _, _, z in points]
    prefixes, acc = [], one
    for z in zz:
        prefixes.append(acc)
        acc = f.mul(acc, z)
    inv = f.inv(acc)
    zi = [0] * len(zz)
    for i in range(len(zz) - 1, -1, -1):
        zi[i] = f.mul(inv, prefixes[i])
        inv = f.mul(inv, zz[i])
    result = []
    for (x, y, z), q in zip(points, zi):
        # Recompute q^2 for y: the ledger intentionally charges five scalings.
        ax = f.mul(x, f.mul(q, q))
        ay = f.mul(y, f.mul(f.mul(q, q), q))
        ox, oy = f.dec(ax), f.dec(ay)
        result.append((ox, oy) if z else None)
    return result


def progression(start, addend, n, f, batch=1024):
    out, state, block = [], jac(start, f), []
    for _ in range(n):
        block.append(state)                 # start at O for the baby j=0
        state = mixed(state, addend, f)      # dependency-respecting Jacobian chain
        if len(block) == batch:
            out.extend(normalize(block, f))
            block = []
    if block:
        out.extend(normalize(block, f))
    return out


def key(record):
    point, _ = record
    return (0, 0) if point is None else (1, point[0])


@dataclass
class Stream:
    left: Any = None
    right: Any = None
    leaf_record: Any = None
    leaf: bool = False
    used: bool = False
    started: bool = False
    a: Any = None
    b: Any = None

    def pull(self):
        if self.leaf:
            if self.used:
                return None
            self.used = True
            return self.leaf_record
        if not self.started:
            self.a, self.b = self.left.pull(), self.right.pull()
            self.started = True
        if self.a is None and self.b is None:
            return None
        if self.b is None or (self.a is not None and key(self.a) <= key(self.b)):
            result, self.a = self.a, self.left.pull()
        else:
            result, self.b = self.b, self.right.pull()
        return result


def make_tree(records, capacity):
    assert capacity > 0 and capacity & (capacity - 1) == 0
    layer = [Stream(leaf=True, leaf_record=records[i] if i < len(records) else None)
             for i in range(capacity)]
    while len(layer) != 1:
        layer = [Stream(left=layer[i], right=layer[i + 1])
                 for i in range(0, len(layer), 2)]
    return layer[0]


def recover(q, gen, p, r, m, h, f):
    babies = progression(None, gen, m, f, batch=4)
    dm = affine_scalar(m, gen, p)
    assert dm is not None
    d = (dm[0], -dm[1] % p)
    giants = progression(q, d, h, f, batch=4)
    cap = 1 << (max(m, h) - 1).bit_length()
    bt = make_tree([(point, j) for j, point in enumerate(babies)], cap)
    gt = make_tree([(point, i) for i, point in enumerate(giants)], cap)
    b, g = bt.pull(), gt.pull()
    verifications = 0
    while b is not None and g is not None:
        if key(b) < key(g):
            b = bt.pull()
        elif key(g) < key(b):
            g = gt.pull()
        else:
            j, i = b[1], g[1]
            # X equality means Q-imP=+/-jP, not c versus c+r.
            # Full y comparison selects the branch; both candidates are still
            # verified here to exercise the conservative two-action bound.
            same = b[0] == g[0]
            candidates = [(i * m + j) % r, (i * m - j) % r]
            if not same:
                candidates.reverse()
            for c in candidates:
                verifications += 1
                if affine_scalar(c, gen, p) == q:
                    return c, verifications
            raise AssertionError('a true X match must resolve as im +/- j')
    return None, verifications
