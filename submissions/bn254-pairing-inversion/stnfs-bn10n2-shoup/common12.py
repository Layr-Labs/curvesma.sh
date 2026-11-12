"""BN254 tower arithmetic: K_h = Q[y]/h, O_h = Z[y] (power basis e_j = y^j), f, g, sigma, norms, sampling."""
from fractions import Fraction as F
import sys as _sys, os as _os
def d95(s):
 o=[]
 for i in range(0,len(s),39):
  t=s[i:i+39]; n=0
  for c in t:n=n*95+ord(c)-32
  o.append(n.to_bytes(((95**len(t)).bit_length()-1)//8,'big'))
 return b''.join(o)
def _install_lib():
 """lib.b95 (base95 lzma text, '@@ name' sections): the byte-identical modules of the accepted BLS12-381 package and descent_witness.py, served to import."""
 import lzma, base64, re, importlib.abc, importlib.util
 here = _os.path.dirname(_os.path.abspath(__file__)); p = _os.path.join(here, 'lib.b95')
 if not _os.path.exists(p) or any(getattr(f, '_bnlib', 0) for f in _sys.meta_path): return
 t = lzma.decompress(d95("".join(open(_os.path.join(here,n)).read() for n in ("lib.b95","libextra.b95")))).decode(); parts = re.split(r"(?m)^@@ (.*)\n", t); src = {parts[i][:-3]: parts[i + 1] for i in range(1, len(parts), 2)}
 class Lib(importlib.abc.MetaPathFinder, importlib.abc.Loader):
  _bnlib = 1
  def find_spec(self, name, path, target=None): return importlib.util.spec_from_loader(name, self, origin=_os.path.join(here, name + '.py')) if name in src else None
  def create_module(self, spec): return None
  def exec_module(self, m): m.__file__ = _os.path.join(here, m.__name__ + '.py'); exec(compile(src[m.__name__], m.__file__, 'exec'), m.__dict__)
 _sys.meta_path.append(Lib())
_install_lib()
from math import isqrt, lcm
H = [1, -1, -2, 5, -2, -1, 1]            # h = y^6 - y^5 - 2y^4 + 5y^3 - 2y^2 - y + 1 (low degree first)
U = 4965661367192848881                  # BN parameter u: p = 36u^4 + 36u^3 + 24u^2 + 6u + 1
W = [-2, 2, 2, -3, 0, 1]                 # w = y^5 - 3y^3 + 2y^2 + 2y - 2, N(w) = 5
SIGP = [1, 2, -5, 2, 1, -1]              # sigma(y) = -y^5 + y^4 + 2y^3 - 5y^2 + 2y + 1
CU = [1, 2, -2, 0, 1, 0]                 # c = y^4 - 2y^2 + 2y + 1, c sigma(c) = -1, sigma(w) = c^2 w
UNITS = [[0, -1, 0, 0, 0, 0], [0, -2, 2, 0, -1, 0], [-1, 1, 0, 0, 0, 0]]
FC = [36, 36, 24, 6, 1]                  # f = sum_k FC[k] w^k x^(8-2k) = w^4 F(x^2/w); g = x^2 - u w
ONE = [1, 0, 0, 0, 0, 0]
def pmul(a, b):
 c = [0] * (len(a) + len(b) - 1)
 for i, x in enumerate(a):
  if x:
   for j, y in enumerate(b): c[i + j] += x * y
 for i in range(len(c) - 1, 5, -1):
  if c[i]:
   for j in range(6): c[i - 6 + j] -= c[i] * H[j]
 return (c + [0] * 6)[:6]
mul = pmul
def inverse(matrix):
 n = len(matrix); a = [[F(x) for x in row] + [F(i == j) for j in range(n)] for i, row in enumerate(matrix)]
 for k in range(n):
  j = next(j for j in range(k, n) if a[j][k]); a[k], a[j] = a[j], a[k]; p = a[k][k]; a[k] = [x / p for x in a[k]]
  for j in range(n):
   if j != k: p = a[j][k]; a[j] = [x - p * y for x, y in zip(a[j], a[k])]
 return [row[n:] for row in a]
def _pw(x, n):
 r = ONE
 for _ in range(n): r = mul(r, x)
 return r
WPOW = [_pw(W, k) for k in range(5)]
UW = [U * x for x in W]
def _sig_rows():
 out = []; s = ONE
 for _ in range(6): out.append(s); s = pmul(s, SIGP)
 return out
S = _sig_rows()
def sigma(a): return [sum(a[i] * S[i][j] for i in range(6)) for j in range(6)]
SD = [[1, 0, 0, -1, -1, -2], [0, 0, -1, -1, -2, 0], [0, 1, 0, 1, -2, 1], [0, 0, 0, 1, -1, 2], [0, 0, 0, 1, 0, 2], [0, 0, -1, 0, -2, 3]]
SB = [[int(x) for x in r] for r in zip(*inverse(SD))]
"""Sieve basis for b: b = sum_j c_j SB[j] with c = SD b (unimodular; SB[0] = 1; SD[5] is the shortest dual vector of Q_b: lines c_5 = k)."""
assert SB[0] == ONE and all(sum(SD[i][k] * SB[j][k] for k in range(6)) == (i == j) for i in range(6) for j in range(6))
def sieve_coords(b): return [sum(SD[j][k] * b[k] for k in range(6)) for j in range(6)]
def sieve_metric(q): return [[sum(SB[i][k] * q[k][l] * SB[j][l] for k in range(6) for l in range(6)) for j in range(6)] for i in range(6)]
def norm(a):
 cols = [mul(a, [0] * i + [1] + [0] * (5 - i)) for i in range(6)]; m = [[cols[j][i] for j in range(6)] for i in range(6)]; prev, sign = 1, 1
 for k in range(5):
  if m[k][k] == 0:
   row = next((j for j in range(k + 1, 6) if m[j][k]), None)
   if row is None: return 0
   m[k], m[row] = m[row], m[k]; sign = -sign
  for i in range(k + 1, 6):
   for j in range(k + 1, 6): value = m[i][j] * m[k][k] - m[i][k] * m[k][j]; assert value % prev == 0; m[i][j] = value // prev
  prev = m[k][k]
 return abs(sign * m[-1][-1])
def fg(a, b):
 """The K_h elements b^8 f(a/b) and b^2 g(a/b) for the candidate a - b x."""
 A = mul(a, a); V = mul(W, mul(b, b)); ap = [ONE, A]; vp = [ONE, V]
 for _ in range(3): ap.append(mul(ap[-1], A)); vp.append(mul(vp[-1], V))
 f = [0] * 6
 for k, c in enumerate(FC): f = [x + c * y for x, y in zip(f, mul(ap[4 - k], vp[k]))]
 return f, [x - U * y for x, y in zip(A, V)]
def norms(a, b):
 f, g = fg(a, b); return norm(f), norm(g)
def integer_interval(a, t, c):
 a, t, c = F(a), F(t), F(c)
 if a <= 0: raise ValueError('non-positive pivot')
 centre = -t / a; radius2 = (t * t - a * c) / (a * a)
 if radius2 < 0: return range(0)
 n, d = centre.numerator, centre.denominator; z = isqrt(radius2.numerator * d * d // radius2.denominator); lo = -(-(n - z) // d); hi = (n + z) // d; return range(lo, hi + 1)
def scaled_log(n, scale=16, upper=False):
 if n < 1 or scale < 1: raise ValueError('positive norm and scale required')
 x = n ** scale; floor = x.bit_length() - 1; return floor + int(upper and x != 1 << floor)
def integer_metric(matrix): den = lcm(*(F(x).denominator for row in matrix for x in row)); q = [[int(F(x) * den) for x in row] for row in matrix]; return q, den, inverse(q)
def quadratic(q, a): return sum(a[i] * q[i][j] * a[j] for i in range(6) for j in range(6))
def rho(u):
 import descent16; return 1. if u <= 1 else descent16.rho_series(u)
_FR = {}
def _frame(q):
 key = str(q)
 if key not in _FR:
  from collections import Counter; import ideal_canonical as IC
  V = IC.lll([[int(i == j) for j in range(6)] for i in range(6)], [[F(x) for x in r] for r in q], Counter()); assert abs(IC.determinant(V)) == 1
  qq = [[sum(V[i][k] * q[k][l] * V[j][l] for k in range(6) for l in range(6)) for j in range(6)] for i in range(6)]; _FR[key] = V, inverse(qq)
 return _FR[key]
def uniform_half_ellipsoid(q, inverse_q, radius2, rng):
 """Uniform integer point of {x5 >= 1, q(x) <= radius2}: rejection from the box of an LLL-reduced frame."""
 V, iq = _frame(q); cap = []
 for j in range(6): c = F(radius2) * iq[j][j]; cap.append(isqrt(c.numerator // c.denominator))
 while True:
  y = [rng.randrange(-c, c + 1) for c in cap]; a = [sum(y[i] * V[i][j] for i in range(6)) for j in range(6)]
  if a[5] >= 1 and quadratic(q, a) <= radius2: return a, quadratic(q, a)
def draw_a(cfg, seed, index):
 """a: uniform integer point of {a5 >= 1, Q_a(a) <= theta A^2} (power coordinates)."""
 from random import Random
 qa, da, ia = integer_metric(cfg['matrix_a']); A2, theta = cfg['production_A_squared'], F(cfg['production_theta'])
 return uniform_half_ellipsoid(qa, ia, int(A2 * theta * da), Random(f"{cfg['sample_tag']}:{seed}:{index}"))
_SM = {}
def draw_b(cfg, seed, index, aq):
 """b: uniform integer point of the half-slice {c5 >= 1, Q_b(b) <= lam A^2 - Q_a(a)}, c = sieve coordinates of b (own random stream)."""
 from random import Random
 key = str(cfg['matrix_b'])
 if key not in _SM: qb, db, _ = integer_metric(cfg['matrix_b']); qs = sieve_metric(qb); _SM[key] = qs, db, inverse(qs), integer_metric(cfg['matrix_a'])[1]
 qs, db, iqs, da = _SM[key]; A2 = cfg['production_A_squared']
 c, bq = uniform_half_ellipsoid(qs, iqs, (F(cfg['production_region_lam']) * A2 - F(aq, da)) * db, Random(f"{cfg['sample_tag']}:b:{seed}:{index}"))
 return [sum(c[j] * SB[j][k] for j in range(6)) for k in range(6)], bq
def regenerate_pair(cfg, seed, index):
 a, aq = draw_a(cfg, seed, index); b, bq = draw_b(cfg, seed, index, aq); return a, b, aq, bq
def stats():
 import json, lzma, base64, os
 return json.loads(lzma.decompress(base64.b85decode(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'stats16.b85')).read())))
