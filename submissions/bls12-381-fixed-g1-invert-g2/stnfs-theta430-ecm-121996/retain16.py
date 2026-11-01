"""Relation-conditioned survivor retention of the descent sieve (DESCENT.md
"Retention"). The sieve credits one log2 p per hit of a degree-one root pair
with p >= 2^8. A relation cell is rejected only when its omitted factor-base"""
import math, functools
import descent16 as D
from common12 import H, U, TERMS
GRID = 4                     # masses in 1/4 bits
G_TERMS = [(2, [1]), (0, [0, U])]        # g(x, y) = x^2 + U y
def _trim(a):
  while a and a[-1] == 0:
    a.pop()
  return a
def pmod(a, m, p):
  a = _trim([x % p for x in a]); dm = len(m) - 1; inv = pow(m[-1], -1, p)
  while len(a) - 1 >= dm:
    c = a[-1] * inv % p; sh = len(a) - 1 - dm
    for i in range(len(m)):
      a[sh + i] = (a[sh + i] - c * m[i]) % p
    _trim(a)
  return a
def pmul(a, b, p):
  c = [0] * (len(a) + len(b) - 1) if a and b else []
  for i, x in enumerate(a):
    if x:
      for j, y in enumerate(b):
        c[i + j] = (c[i + j] + x * y) % p
  return _trim(c)
def padd(a, b, p):
  return _trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(max(len(a), len(b)))])
def pscal(a, c, p):
  return _trim([x * c % p for x in a])
def pdiv(a, b, p):
  a = _trim([x % p for x in a]); b = _trim([x % p for x in b]); q = [0] * max(1, len(a) - len(b) + 1); inv = pow(b[-1], -1, p)
  while a and len(a) >= len(b):
    c = a[-1] * inv % p; sh = len(a) - len(b); q[sh] = c
    for i in range(len(b)):
      a[sh + i] = (a[sh + i] - c * b[i]) % p
    _trim(a)
  return _trim(q), a
def pgcd(a, b, p):
  a = _trim([x % p for x in a]); b = _trim([x % p for x in b])
  while b:
    a, b = b, pmod(a, b, p)
  if a:
    inv = pow(a[-1], -1, p); a = [x * inv % p for x in a]
  return a
def pderiv(a, p):
  return _trim([(i * a[i]) % p for i in range(1, len(a))])
def ppowmod(a, e, m, p):
  r = [1]; a = pmod(a, m, p)
  while e:
    if e & 1:
      r = pmod(pmul(r, a, p), m, p)
    a = pmod(pmul(a, a, p), m, p); e >>= 1
  return r
def roots_Fp(a, p):
  return [r for r in range(p) if sum(c * pow(r, i, p) for i, c in enumerate(a)) % p == 0]
def factor_Fp(a, p, seed=1):
  """irreducible factors of a squarefree polynomial over F_p: roots by search, then distinct-degree and Cantor-Zassenhaus splitting."""
  import random
  rng = random.Random(seed)
  def split(a, d):
    if len(a) - 1 == d:
      return [a]
    for _ in range(10000):
      r = _trim([rng.randrange(p) for _ in range(len(a) - 1)])
      if len(r) < 2:
        continue
      if p == 2:
        t = r; s = r
        for _ in range(d - 1):
          t = pmod(pmul(t, t, p), a, p); s = padd(s, t, p)
        g = pgcd(a, s, p)
      else:
        g = pgcd(a, padd(ppowmod(r, (p ** d - 1) // 2, a, p), [p - 1], p), p)
      if 1 < len(g) - 1 < len(a) - 1:
        q, _ = pdiv(a, g, p); return split(g, d) + split(q, d)
    raise RuntimeError('splitting failed')
  out = []; rest = _trim([x % p for x in a])
  for r in roots_Fp(rest, p):
    out.append([(-r) % p, 1]); rest, _ = pdiv(rest, [(-r) % p, 1], p)
  d = 2; x = [0, 1]
  xp = ppowmod(x, p, rest, p) if len(rest) > 1 else x
  while len(rest) - 1 >= 2 * d:
    xp = ppowmod(xp, p, rest, p)
    g = pgcd(rest, padd(xp, pscal(x, p - 1, p), p), p)
    if len(g) > 1:
      out += split(g, d); rest, _ = pdiv(rest, g, p)
      if len(rest) > 1:
        xp = pmod(xp, rest, p)
    d += 1
  if len(rest) > 1:
    out.append(rest)
  return out
def emul(a, b, m, p):
  return pmod(pmul(a, b, p), m, p) if a and b else []
def einv(a, m, p):
  r0, r1 = list(m), _trim(list(a)); s0, s1 = [], [1]
  while r1:
    q, r = pdiv(r0, r1, p); r0, r1 = r1, r; s0, s1 = s1, padd(s0, pscal(pmul(q, s1, p), p - 1, p), p)
  return pmod(pscal(s0, pow(r0[0], -1, p), p), m, p)
def Etrim(A):
  while A and not A[-1]:
    A.pop()
  return A
def Emul(A, B, m, p):
  C = [[] for _ in range(len(A) + len(B) - 1)] if A and B else []
  for i, x in enumerate(A):
    if x:
      for j, y in enumerate(B):
        if y:
          C[i + j] = padd(C[i + j], emul(x, y, m, p), p)
  return Etrim(C)
def Emod(A, M, m, p):
  A = Etrim([list(x) for x in A]); dm = len(M) - 1; inv = einv(M[-1], m, p)
  while len(A) - 1 >= dm:
    c = emul(A[-1], inv, m, p); sh = len(A) - 1 - dm
    for i in range(len(M)):
      A[sh + i] = padd(A[sh + i], pscal(emul(c, M[i], m, p), p - 1, p), p)
    Etrim(A)
  return A
def Egcd(A, B, m, p):
  A = Etrim([list(x) for x in A]); B = Etrim([list(x) for x in B])
  while B:
    A, B = B, Emod(A, B, m, p)
  if A:
    inv = einv(A[-1], m, p); A = [emul(x, inv, m, p) for x in A]
  return A
def Ederiv(A, p):
  return Etrim([pscal(A[i], i % p, p) for i in range(1, len(A))])
def Epowmod(A, e, M, m, p):
  R = [[1]]; A = Emod(A, M, m, p)
  while e:
    if e & 1:
      R = Emod(Emul(R, A, m, p), M, m, p)
    A = Emod(Emul(A, A, m, p), M, m, p); e >>= 1
  return R
def Equo(A, M, m, p):
  """quotient of A by the monic M over F_p[t]/(m)."""
  A = Etrim([list(x) for x in A]); dm = len(M) - 1; Q = [[] for _ in range(max(0, len(A) - dm))]
  while len(A) - 1 >= dm:
    c = A[-1]; sh = len(A) - 1 - dm; Q[sh] = c
    for i in range(len(M)):
      A[sh + i] = padd(A[sh + i], pscal(emul(c, M[i], m, p), p - 1, p), p)
    Etrim(A)
  return Etrim(Q)
def root_multiplicities(A, q, m, p):
  """multiplicities of the roots of A in F_q = F_p[t]/(m), sorted: R_1 = gcd(A, X^q - X) has the distinct roots, R_(k+1) = gcd(A / (R_1 ... R_k), R_k) the roots of multiplicity > k."""
  A = Etrim([list(x) for x in A])
  XP = Epowmod([[], [1]], q, A, m, p)
  diff = Etrim([padd(XP[i] if i < len(XP) else [], pscal([1] if i == 1 else [], p - 1, p), p) for i in range(max(len(XP), 2))])
  R = Egcd(A, diff, m, p); counts = []
  while len(R) > 1:
    counts.append(len(R) - 1); A = Equo(A, R, m, p); R = Egcd(A, R, m, p)
  counts.append(0)
  return [k + 1 for k in range(len(counts) - 1) for _ in range(counts[k] - counts[k + 1])]

def poly_over_ext(terms, m, p):
  """f(x, eta) with eta = t in F_p[t]/(m): coefficients in x."""
  deg = max(d for d, _ in terms); co = [[] for _ in range(deg + 1)]
  for d, c in terms:
    v = []
    for i, ci in enumerate(c):
      if ci:
        v = padd(v, pscal(ppowmod([0, 1], i, m, p), ci % p, p), p)
    co[d] = pmod(v, m, p) if v else []
  return co
@functools.lru_cache(None)
def structure(p):
  """prime ideals (p, h_i) of the base field: for every irreducible factor h_i of H mod p (Kummer-Dedekind, O_K = Z[y]) its degree d1, exponent e and the root multiplicities of f and g in F_q, q = p^d1."""
  Hp = [x % p for x in H]; sq = pgcd(Hp, pderiv(Hp, p), p)
  out = []
  for hi in factor_Fp(pdiv(Hp, sq, p)[0] if len(sq) > 1 else Hp, p):
    e, rest = 0, Hp
    while True:
      quo, rem = pdiv(rest, hi, p)
      if rem:
        break
      rest, e = quo, e + 1
    d1 = len(hi) - 1
    out.append([d1, e] + [root_multiplicities(poly_over_ext(t, hi, p), p ** d1, hi, p) for t in (TERMS, G_TERMS)])
  return out
def primes(lo, hi):
  return [p for p in range(lo, hi) if all(p % q for q in range(2, math.isqrt(p) + 1))]
def _pdf_add(pdf, mass_units, prob):
  pdf[mass_units] = pdf.get(mass_units, 0.0) + prob
def local_pdf(p, st, side, sieved, up):
  """P(X_p = x) of a base-primitive cell (DESCENT.md): at a prime ideal with residue field F_q, b is a unit with probability q/(q+1) and a root of multiplicity mu of the side polynomial then has v >= j with probability at most min(1/q, mu q^-ceil(j/mu)), exactly q^-j when mu = 1 (Hensel); the ideals above p are independent. One unit of v weighs d1 log2 p (the first power is credited at sieved primes, whose degree > 1 ideals are excluded instead)."""
  lp = math.log2(p); rnd = math.ceil if up else math.floor
  pdf = {0: 1.0}
  for d1, e, mf, mg in st:
    mults = mf if side == 'f' else mg
    if not mults or (sieved and d1 != 1):
      continue
    q = p ** d1
    tail = lambda j: q / (q + 1) * sum(min(1.0 / q, mu * q ** (-((j + mu - 1) // mu))) for mu in mults)
    comp = {}; prev = tail(1); comp[0] = 1 - prev
    for j in range(1, 1000):
      nxt = tail(j + 1)
      _pdf_add(comp, rnd(d1 * (j - 1 if sieved else j) * lp * GRID), prev - nxt)
      prev = nxt
      if nxt < 1e-40:
        _pdf_add(comp, rnd(d1 * (j + (0 if sieved else 1)) * lp * GRID), nxt); break
    pdf = _convolve(pdf, comp)
  return pdf
@functools.lru_cache(None)
def base_primitive(pmax_exact=512):
  """lower bound of the base-primitive fraction: prod (1 - p^(-2 d1)) over the prime ideals above p < 2^9 exactly, (1 - 6 p^-2) per prime above."""
  f = 1.0
  for p in primes(2, pmax_exact):
    for d1, e, mf, mg in structure_of(p):
      f *= 1 - p ** (-2.0 * d1)
  j0 = int(math.log2(pmax_exact))
  for j in range(j0, 200):
    f *= 1 - 6.0 * 1.26 * 2.0 ** (j + 1) / ((j + 1) * math.log(2)) * 2.0 ** (-2 * j)
  return f
def _convolve(a, b):
  out = {}
  for x, px in a.items():
    for y, py in b.items():
      out[x + y] = out.get(x + y, 0.0) + px * py
  return out
def _prune(pdf):
  """entries below 1e-30 are moved to the largest kept mass (a larger mass only raises the loss)."""
  keep = {x: v for x, v in pdf.items() if v > 1e-30}
  lost = sum(v for x, v in pdf.items() if x not in keep)
  if lost:
    keep[max(keep)] += lost
  return keep

def large_block_pdf(j, up):
  """excess powers of sieved primes in [2^j, 2^(j+1)): <= 72 p^-2 per level over <= 1.26 x/ln x primes; one and two events."""
  lam = 72.0 * 1.26 * 2.0 ** (j + 1) / ((j + 1) * math.log(2)) * 2.0 ** (-2 * j)
  pdf = {0: 1.0}
  if lam < 1e-30:
    return pdf
  comp = {}; tot = 0.0
  mass = (j + 1) if up else j
  for t in range(1, 12):
    pr = lam * 2.0 ** (-j * (t - 1)) * (1 - 2.0 ** -j)
    _pdf_add(comp, math.ceil(mass * t * GRID) if up else math.floor(mass * t * GRID), pr); tot += pr
  for t in range(2, 8):
    pr = lam * lam / 2 * 2.0 ** (-j * (t - 2)) * (t - 1)
    _pdf_add(comp, math.ceil(mass * t * GRID) if up else math.floor(mass * t * GRID), pr); tot += pr
  comp[0] = max(0.0, 1 - tot)
  return comp
@functools.lru_cache(None)
def omitted_pdf(side, up, pmax_exact=512):
  """distribution of the omitted mass X (GRID units) of a base-primitive cell on one side, and a lower bound of the ordinary fraction."""
  pdf = {0: 1.0}; ordinary = 1.0
  for p in primes(2, pmax_exact):
    st = structure_of(p)
    if p >= 256:
      for d1, e, mf, mg in st:
        if d1 > 1:
          ordinary -= len(mf if side == 'f' else mg) * p ** (-d1)
    pdf = _prune(_convolve(pdf, local_pdf(p, st, side, p >= 256, up)))
  j0 = int(math.log2(pmax_exact))
  for j in range(j0, 69):
    pdf = _prune(_convolve(pdf, large_block_pdf(j, up)))
    ordinary -= (36.0 if side == 'f' else 6.0) * 1.26 * 2.0 ** (j + 1) / ((j + 1) * math.log(2)) * 2.0 ** (-2 * j)
  return pdf, ordinary
STRUCT_TABLE = None
def load_table(stats):
  global STRUCT_TABLE
  STRUCT_TABLE = stats['prime_structure']['table']
  omitted_pdf.cache_clear(); base_primitive.cache_clear()
def structure_of(p):
  """the stored table (stats16.b85 prime_structure, recomputed by [11b] at the five primes where H, f or g is not squarefree and spot-checked elsewhere) or structure(p)."""
  return STRUCT_TABLE[str(p)] if STRUCT_TABLE is not None and str(p) in STRUCT_TABLE else structure(p)

@functools.lru_cache(None)
def ss_upper2(m, base, n, iA, kA, kB, mesh, alpha):
  """upper sum companion of descent16.ss_lower2: rho at the right endpoints, upper bracket, every cell meeting the region kept, floor(m)."""
  m = math.floor(m)
  if n <= base:
    return D.rho_hi((m + alpha) / base)
  rows, h = D.window_rows(n, base, iA, mesh)
  tot = 0.0
  for (jA, jB), row in rows.items():
    if jA > kA or (jB and not kB):
      continue
    j = jA + jB
    for i, weight in enumerate(row):
      if not weight:
        continue
      left = j * base + i * h
      if left > m:
        continue
      u = (m - left - j * h + alpha) / base
      tot += weight * (1.0 if u <= 1 else D.rho_hi(u))
  return min(tot, 1.0)
@functools.lru_cache(None)
def loss(side, m, n, iA, kA, kB, alpha, slack, step=64):
  """P(rejected | relation) = sum_{x > s} P(X = x) Pss_hi(m - x) / sum_{x <= s} P(X = x) Pss_lo(m - x) (DESCENT.md)."""
  Pm = D.P
  s_units = int(math.floor((slack - 1.0 - m / 128.0) * GRID))
  pdf_up, _ = omitted_pdf(side, True)
  base = Pm.LX
  def hi(x_bits):
    xg = math.ceil(x_bits / step) * step
    return ss_upper2(max(m - xg, 1.0), base, n, iA, kA, kB, Pm.MESH, alpha) if m - xg > 1 else 1.0
  def lo(x_bits):
    xg = math.floor(x_bits / step) * step
    return D.ss_lower2(max(m - xg, 1.0), base, n, iA, kA, kB, Pm.MESH, alpha, Pm.TRUNC) if m - xg > 1 else 1.0
  num = sum(pr * hi(x / GRID) for x, pr in pdf_up.items() if x > s_units)
  den = sum(pr * lo(x / GRID) for x, pr in pdf_up.items() if x <= s_units)
  return num / den if den > 0 else 1.0
def counted_fraction():
  """base-primitive cells minus the two (disjoint) non-factor-base exclusions."""
  _, ordf = omitted_pdf('f', True); _, ordg = omitted_pdf('g', True)
  return max(0.0, ordf + ordg - 1) * base_primitive()


def retention(nf, ng, n, iA, kA, kB, slack):
  """mu factor of a sieve node: the base-primitive and ordinary fractions times 1 - loss_f - loss_g."""
  lf = loss('f', math.ceil(nf), n, iA, kA, kB, D.P.ALPHA_F, slack)
  lg = loss('g', math.ceil(ng), n, iA, kA, kB, D.P.ALPHA_G, slack)
  _, ordf = omitted_pdf('f', True); _, ordg = omitted_pdf('g', True)
  return max(0.0, (1 - lf - lg)) * max(0.0, ordf + ordg - 1) * base_primitive(), lf, lg
def cell_loss(side, slack):
  """unconditional P(X > slack): the fraction of all counted cells rejected by the omitted mass alone."""
  pdf, _ = omitted_pdf(side, True)
  return sum(pr for x, pr in pdf.items() if x > slack * GRID)
