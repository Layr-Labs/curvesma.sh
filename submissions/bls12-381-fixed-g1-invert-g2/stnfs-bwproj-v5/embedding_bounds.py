""
from fractions import Fraction as F
from math import comb, cos, sin, pi, nextafter, inf, isfinite, frexp, isqrt, sqrt
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
  ""
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
  ""
  if lo <= 0 or not isfinite(lo) or not isfinite(hi) or hi < lo:
    return None
  m, e = frexp(lo)
  lower = 16*(e-1)+bisect_right(LOG_UPPER, 2*m)-1
  m, e = frexp(hi)
  upper = 16*(e-1)+bisect_left(LOG_LOWER, 2*m)
  return lower, upper
EPS50 = 2.0**-50
KAPPA = 2.0**-32
UP, DN, UP44 = 1+2.0**-50, 1-2.0**-50, 1+2.0**-44
def _iv_mid_rad(v):
  m = 0.5*(v.lo+v.hi)
  return m, nextafter(max(m-v.lo, v.hi-m), inf)
def _mid_rad(v, real):
  """centre (float or complex pair) and radius r with |x - centre| <= r for x in v."""
  if real:
    return _iv_mid_rad(v)
  (mr, rr), (mi, ri) = _iv_mid_rad(v.real), _iv_mid_rad(v.imag)
  return (mr, mi), nextafter(sqrt(rr*rr+ri*ri), inf)*UP
def _fast_place(real, root, one, a2, uy, coefficients):
  """centres cm and majorants ch >= |c| of the coefficients of F_v(w) = sum c_k w^k and of
  G_v = a2 + uy w."""
  def cen_maj(v):
    c, r = _mid_rad(v, real)
    a = abs(c) if real else sqrt(c[0]*c[0]+c[1]*c[1])*UP
    h = (a+r)*UP
    if r > h*2.0**-34 and r > 2.0**-1000:
      raise ArithmeticError('coefficient interval too wide')
    return c, h
  cm, ch = zip(*[cen_maj(c) for c in coefficients])
  (a2m, a2h), (uym, uyh) = cen_maj(a2), cen_maj(uy)
  return real, cm, ch, a2m, a2h, uym, uyh
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
    self.fast = [_fast_place(*p) for p in self.places]
  def row(self, suffix):
    contexts = []
    for (real, root, one, a2, uy, coefficients), fp in zip(self.places, self.fast):
      bv = one*0
      for coefficient in reversed(suffix):
        bv = bv*root+coefficient
      contexts.append((fp, _mid_rad(bv*root, real)))
    return contexts
  def block(self, contexts, lo, hi):
    self.counts['block_evaluations'] += 1
    m, r0 = 0.5*(lo+hi), 0.5*(hi-lo)
    prod = [[1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0]]
    ok = [True, True]
    for (real, cm, ch, a2m, a2h, uym, uyh), (bm, bd) in contexts:
      if real:
        x = bm+m
        s = abs(x)
        w = x*x
        sw = w
      else:
        xr, xi = bm[0]+m, bm[1]
        s = sqrt(xr*xr+xi*xi)*UP
        w = (xr*xr-xi*xi, 2*xr*xi)
        sw = s*s
      rho = (r0+bd+s*EPS50)*UP
      rw = (rho*(s+s+rho)+sw*EPS50)*UP
      t = (sw+rw)*UP
      ph = ch[6]
      for k in range(5, -1, -1):
        ph = ph*t+ch[k]
      eg = (rw*uyh+KAPPA*(a2h+uyh*t))*UP44
      if real:
        T = list(cm)
        for i in range(6):
          for k in range(5, i-1, -1):
            T[k] = T[k]+T[k+1]*w
        e1 = 0.0
        for k in range(6, 0, -1):
          e1 = (e1+abs(T[k]))*rw
        g = a2m+uym*w
        vals = ((abs(T[0]), (e1+KAPPA*ph)*UP44), (abs(g), eg))
      else:
        Tr, Ti = [c[0] for c in cm], [c[1] for c in cm]
        for i in range(6):
          for k in range(5, i-1, -1):
            Tr[k], Ti[k] = Tr[k]+(Tr[k+1]*w[0]-Ti[k+1]*w[1]), Ti[k]+(Tr[k+1]*w[1]+Ti[k+1]*w[0])
        e1 = 0.0
        for k in range(6, 0, -1):
          e1 = (e1+sqrt(Tr[k]*Tr[k]+Ti[k]*Ti[k]))*rw
        gr, gi = a2m[0]+(uym[0]*w[0]-uym[1]*w[1]), a2m[1]+(uym[0]*w[1]+uym[1]*w[0])
        vals = ((sqrt(Tr[0]*Tr[0]+Ti[0]*Ti[0]), (e1*UP+KAPPA*ph)*UP44), (sqrt(gr*gr+gi*gi), eg))
      for side, (v, e) in enumerate(vals):
        if real:
          vlo, vhi = (v-e)*DN, (v+e)*UP
        else:
          vlo, vhi = (v*DN-e)*DN, (v*UP+e)*UP
          vlo, vhi = (vlo*vlo*DN if vlo > 0 else vlo), vhi*vhi*UP
        k = 0 if real else 2
        if not vlo > 0:
          ok[side] = False
        prod[side][k] *= vlo
        prod[side][k+1] *= vhi
    logs = [None, None]
    for side in range(2):
      if ok[side]:
        b0 = scaled_log_float_bounds(prod[side][0]*DN, prod[side][1]*UP)
        b1 = scaled_log_float_bounds(prod[side][2]*DN, prod[side][3]*UP)
        if b0 is not None and b1 is not None:
          logs[side] = [b0[0]+b1[0], b0[1]+b1[1]]
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
