from fractions import Fraction as F
from math import comb, cos, sin, pi, nextafter, inf, isfinite, frexp, isqrt
from bisect import bisect_left, bisect_right
from collections import Counter
from common12 import H, U, norms
def cadd(a, b):
    return a[0]+b[0], a[1]+b[1]
def cmul(a, b):
    return a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]
def cdiv(a, b):
    d = b[0]*b[0]+b[1]*b[1]
    return (a[0]*b[0]+a[1]*b[1])/d, (a[1]*b[0]-a[0]*b[1])/d
def ceval(coefficients, z):
    out = (F(0), F(0))
    for c in reversed(coefficients):
        out = cadd(cmul(out, z), (F(c), F(0)))
    return out
def root_disks():
    roots = [1.3*complex(cos(2*pi*(i+.2)/6), sin(2*pi*(i+.2)/6)) for i in range(6)]
    for _ in range(150):
        new = []
        for i, z in enumerate(roots):
            numerator = sum(c*z**j for j, c in enumerate(H))
            denominator = 1
            for j, w in enumerate(roots):
                if i != j:
                    denominator *= z-w
            new.append(z-numerator/denominator)
        if max(abs(x-y) for x, y in zip(new, roots)) < 1e-15:
            roots = new
            break
        roots = new
    quantum, radius = 1 << 128, F(1, 1 << 100)
    derivative = [j*H[j] for j in range(1, len(H))]
    disks = []
    for z in roots:
        centre = (F(z.real), F(0) if abs(z.imag) < 1e-12 else F(z.imag))
        for _ in range(8):
            correction = cdiv(ceval(H, centre), ceval(derivative, centre))
            updated = (centre[0]-correction[0], centre[1]-correction[1])
            centre = tuple(F(round(x*quantum), quantum) for x in updated)
        coefficients = []
        for k in range(7):
            coefficients.append(ceval([comb(j, k)*H[j] for j in range(k, 7)], centre))
        lower_linear = max(abs(x) for x in coefficients[1])*radius
        upper_remainder = sum(abs(x) for x in coefficients[0])
        upper_remainder += sum(sum(abs(x) for x in coefficients[k])*radius**k for k in range(2, 7))
        if not lower_linear > upper_remainder:
            raise AssertionError('uncertified root disk')
        disks.append((centre, radius))
    assert all(max(abs(c[i]-d[i]) for i in range(2)) > r+s
               for j, (c, r) in enumerate(disks) for d, s in disks[j+1:])
    assert sum(c[1] == 0 for c, _ in disks) == 2
    assert sum(c[1] > r for c, r in disks) == 2
    return sorted(disks, key=lambda x: (x[0][1] != 0, x[0][0], x[0][1]))
class Interval:
    def __init__(self, lo, hi=None, counts=None):
        self.lo, self.hi = float(lo), float(lo if hi is None else hi)
        self.counts = Counter() if counts is None else counts
        if self.lo > self.hi or not isfinite(self.lo) or not isfinite(self.hi):
            raise ArithmeticError('nonfinite/empty floating interval')
    @classmethod
    def rational(cls, lo, hi=None, counts=None):
        lo, hi = F(lo), F(lo if hi is None else hi)
        return cls(nextafter(float(lo), -inf), nextafter(float(hi), inf), counts)
    def __add__(self, other):
        if not isinstance(other, Interval):
            other = Interval.rational(other, counts=self.counts)
        self.counts['interval_additions'] += 1
        return Interval(nextafter(self.lo+other.lo, -inf),
                        nextafter(self.hi+other.hi, inf), self.counts)
    __radd__ = __add__
    def __neg__(self):
        return Interval(-self.hi, -self.lo, self.counts)
    def __sub__(self, other):
        return self+-other
    def __mul__(self, other):
        if not isinstance(other, Interval):
            other = Interval.rational(other, counts=self.counts)
        self.counts['interval_multiplications'] += 1
        candidates = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(nextafter(min(candidates), -inf), nextafter(max(candidates), inf), self.counts)
    __rmul__ = __mul__
    def absolute(self):
        lo = 0.0 if self.lo <= 0 <= self.hi else min(abs(self.lo), abs(self.hi))
        return Interval(lo, max(abs(self.lo), abs(self.hi)), self.counts)
    def square(self):
        a = self.absolute()
        self.counts['interval_squares'] += 1
        return Interval(max(0.0, nextafter(a.lo*a.lo, -inf)),
                        nextafter(a.hi*a.hi, inf), self.counts)
class ComplexInterval:
    def __init__(self, real, imag):
        self.real, self.imag = real, imag
    @classmethod
    def constant(cls, value, counts):
        return cls(Interval.rational(value, counts=counts), Interval(0, counts=counts))
    def __add__(self, other):
        if not isinstance(other, ComplexInterval):
            other = ComplexInterval.constant(other, self.real.counts)
        return ComplexInterval(self.real+other.real, self.imag+other.imag)
    __radd__ = __add__
    def __neg__(self):
        return ComplexInterval(-self.real, -self.imag)
    def __sub__(self, other):
        return self+-other
    def __mul__(self, other):
        if not isinstance(other, ComplexInterval):
            other = ComplexInterval.constant(other, self.real.counts)
        return ComplexInterval(self.real*other.real-self.imag*other.imag,
                               self.real*other.imag+self.imag*other.real)
    __rmul__ = __mul__
    def absolute_square(self):
        value = self.real.square()+self.imag.square()
        return Interval(max(0.0, value.lo), value.hi, self.real.counts)
def log_table():
    lower, upper = [], []
    for j in range(17):
        value = 1 << (j+16*52)
        k = value
        for _ in range(4):
            k = isqrt(k)
        assert k**16 <= value < (k+1)**16
        lower.append(k/(1 << 52))
        upper.append((k+int(k**16 != value))/(1 << 52))
    return lower, upper
LOG_LOWER, LOG_UPPER = log_table()
def scaled_log_float_bounds(lo, hi):
    if lo <= 0 or not isfinite(lo) or not isfinite(hi) or hi < lo:
        return None
    m, e = frexp(lo)
    lower = 16*(e-1)+bisect_right(LOG_UPPER, 2*m)-1
    m, e = frexp(hi)
    upper = 16*(e-1)+bisect_left(LOG_LOWER, 2*m)
    return lower, upper
class NormBounds:
    def __init__(self, a, disks=None):
        if len(a) != 6 or any(type(x) is not int for x in a):
            raise ValueError('six integer a coefficients required')
        self.a = tuple(a)
        self.counts, self.places = Counter(), []
        for (r, i), radius in (root_disks() if disks is None else disks):
            if i < -radius:
                continue
            real = i == 0
            root = (Interval.rational(r-radius, r+radius, self.counts) if real else
                    ComplexInterval(Interval.rational(r-radius, r+radius, self.counts),
                                    Interval.rational(i-radius, i+radius, self.counts)))
            one = Interval(1, counts=self.counts) if real else ComplexInterval.constant(1, self.counts)
            av = one*0
            for coefficient in reversed(a):
                av = av*root+coefficient
            ap = [one]
            yp = [one]
            for _ in range(12):
                ap.append(ap[-1]*av)
            for _ in range(6):
                yp.append(yp[-1]*root)
            coefficients = [ap[12], -2*root*ap[10], one*0, 2*yp[3]*ap[6],
                            one*0, yp[5]*ap[2], yp[6]]
            self.places.append((real, root, one, ap[2], U*root, coefficients))
        assert len(self.places) == 4
    def row(self, suffix):
        contexts = []
        for real, root, one, a2, uy, coefficients in self.places:
            bv = one*0
            for coefficient in reversed(suffix):
                bv = bv*root+coefficient
            contexts.append((real, bv*root, a2, uy, coefficients))
        return contexts
    def block(self, contexts, lo, hi):
        self.counts['block_evaluations'] += 1
        logs = [[0, 0], [0, 0]]
        for real, base, a2, uy, coefficients in contexts:
            b0 = Interval(lo, hi, self.counts)
            bv = base+b0 if real else base+ComplexInterval(b0, Interval(0, counts=self.counts))
            w = bv*bv
            fv = coefficients[6]
            for k in range(5, -1, -1):
                fv = fv*w
                if k not in (2, 4):
                    fv = fv+coefficients[k]
            gv = a2+uy*w
            for side, value in enumerate((fv, gv)):
                magnitude = value.absolute() if real else value.absolute_square()
                bound = scaled_log_float_bounds(magnitude.lo, magnitude.hi)
                self.counts['magnitude_log_bounds'] += 1
                if bound is None:
                    logs[side] = None
                elif logs[side] is not None:
                    logs[side][0] += bound[0]
                    logs[side][1] += bound[1]
        return logs
    def adaptive_row(self, a, suffix, lo, hi, max_nodes=10000, singleton_policy='exact', width=32):
        if tuple(a) != self.a or len(suffix) != 5:
            raise ValueError('a/context or suffix dimension mismatch')
        if type(lo) is not int or type(hi) is not int or lo > hi:
            raise ValueError('nonempty integer row interval required')
        if singleton_policy not in ('exact', 'reject'):
            raise ValueError('unknown singleton policy')
        contexts, pending, leaves = self.row(suffix), [(lo, hi)], []
        start_nodes = self.counts['block_evaluations']
        while pending:
            if self.counts['block_evaluations']-start_nodes >= max_nodes:
                raise RuntimeError('embedding-bound node budget exhausted')
            lo, hi = pending.pop()
            bounds = self.block(contexts, lo, hi)
            if all(v is not None and v[1]-v[0] <= width for v in bounds):
                leaves.append((lo, hi, bounds))
                self.counts['accepted_blocks'] += 1
            elif lo == hi:
                if singleton_policy == 'reject':
                    leaves.append((lo, hi, {'rejected': True}))
                    self.counts['rejected_unresolved_points'] += 1
                else:
                    ns = norms(a, [lo]+list(suffix))
                    leaves.append((lo, hi, {'exact_norms': ns}))
                    self.counts['exact_norm_fallbacks'] += 1
            else:
                middle = (lo+hi)//2
                pending.extend(((middle+1, hi), (lo, middle)))
                self.counts['block_splits'] += 1
        return leaves
