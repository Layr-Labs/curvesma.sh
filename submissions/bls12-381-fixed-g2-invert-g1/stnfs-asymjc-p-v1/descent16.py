"""stnfs descent model (DESCENT.md): nodes, lower sums, survivors, the work recurrence."""
import math, functools
from common12 import ALPHA_F_CHARGED
from decimal import Decimal, localcontext
LN2 = math.log(2); UNIT = 20000.0
KMAX = 64
RSTEP = 1.0 / 1024
def _rho_coefficients(kmax=KMAX, nt=600):
  with localcontext() as ctx:
    ctx.prec = 200
    zero, one = Decimal(0), Decimal(1)
    coefs = [None, [one] + [zero] * (nt - 1)]
    for k in range(2, kmax + 1):
      previous, c = coefs[k - 1], [zero] * nt
      for m in range(nt - 1):
        c[m + 1] = (previous[m] + m * c[m]) / (k * (m + 1))
      c[0] = previous[0] - sum(c[1:])
      coefs.append(c)
    out = [None]
    for c in coefs[1:]:
      big = max(abs(x) for x in c)
      while len(c) > 1 and abs(c[-1]) < big * Decimal('1e-25'):
        c = c[:-1]
      out.append([float(x) for x in c])
    return out
RHO_COEFS = _rho_coefficients()
def rho_series(u):
  if u <= 1:
    return 1.0
  k = math.ceil(u)
  if k > KMAX:
    return 0.0
  x, value = k - u, 0.0
  for c in reversed(RHO_COEFS[k]):
    value = value * x + c
  return max(value, 0.0)
RHO_TAB = [rho_series(i * RSTEP) for i in range(int(KMAX / RSTEP) + 2)]
def rho_lo(u):
  """lower bracket of rho(u): the table value at the next grid point."""
  if u <= 1:
    return 1.0
  i = math.ceil(u / RSTEP)
  return RHO_TAB[i] if i < len(RHO_TAB) else 0.0
def rho_hi(u):
  """upper bracket: the previous grid point;"""
  if u <= 1:
    return 1.0
  i = math.floor(u / RSTEP)
  if i + 1 >= len(RHO_TAB):
    return 1.0 / math.factorial(math.floor(u))
  return RHO_TAB[i]
def kappa(u, t):
  return rho_series(u - t) / rho_series(u) if rho_series(u) > 0 else 0.0
CALIB = ((1.0, 2.0, 0.17), (2.0, 3.0, 0.11), (3.0, 5.0, 0.052), (5.0, 8.0, 0.022))   # (lambda lo, hi, c_r): cover of the measured fraction of OCCURRING ideals still without a logarithm after purge and back-substitution, 1.5 x CALIB_MAX rounded up ([6b], [7c], DESCENT.md)
CALIB_MAX = {(1.0, 2.0): 0.1129, (2.0, 3.0): 0.0732, (3.0, 5.0): 0.0326, (5.0, 8.0): 0.0082}   # largest measured unknown fraction among occurring ideals over the 16 filtsim.missing_recon runs (seeds 1-3 at 2^12..2^16, seed 1 at 2^17; means 741.8/498.8, excess 1.1)
ZERO_EPS = 2.0 ** -8   # relative allowance on the zero-occurrence missing mass for the realized main-phase table (Chebyshev, [7m]; DESCENT.md)
def p_missing(lam):
  """probability that a factor-base ideal of expected occurrence lambda has no logarithm: the ideals that never occur
  ((1 + ZERO_EPS) e^-lambda, unrecoverable) plus the occurring ones still unknown after purge and back-substitution,
  charged at the band's cover c_r of the measured conditional fraction; lambda < 1 and lambda >= 8 keep the
  uncalibrated P(Poisson <= 2) (DESCENT.md, Factor-base ideals without a logarithm)."""
  p0 = math.exp(-lam)
  for lo, hi, c in CALIB:
    if lo <= lam < hi:
      return (1.0 + ZERO_EPS) * p0 + (1.0 - p0) * c
  return p0 * (1 + lam + lam * lam / 2)
@functools.lru_cache(None)
def miss_side(u4, log2_rels, u_main, L_main=69.0):
  """missing-log ideals on one m-bit side (DESCENT.md): sum over 1/10-bit bands b < 69 of kappa(u, b/69) ln((b+.1)/b) P_miss(lambda)."""
  u = u4 / 4.0; tot = 0.0
  for b10 in range(400, 690):
    b = b10 / 10.0
    lam = 2.0 ** (log2_rels - b) * 2 * kappa(u_main, b / L_main)
    tot += kappa(u, b / 69.0) * math.log((b + 0.1) / b) * p_missing(lam)
  return tot
U_GRID = (7.0, 8.5, 10.0, 12.0, 14.0, 17.0, 20.0, 23.0, 26.0, 29.0, 32.0, 35.0, 38.0, 40.0)   # relation norm sizes u = m/69 spanning the shipped nodes and the boot
BAND_MAX = 2.0 ** 0.1   # within a 1/10-bit band the occurrence weight varies as 1/q: the largest ideal's weight is at most 2^0.1 x the band average
def band_weights(um, log2_rels, L_main=69.0):
  """per 1/10-bit band b < 69 on a side with main-phase norm size um: (b, n_b, P_miss, e^-lambda) with n_b = (li(2^(b+0.1)) -
  li(2^b))/2 sigma-orbits of degree-one ideals (NBAR = 2 ideals per prime, quotiented by sigma; model16.log2_U) and lambda_b
  the occurrence model of miss_side."""
  out = []
  for b10 in range(400, 690):
    b = b10 / 10.0
    lam = 2.0 ** (log2_rels - b) * 2 * kappa(um, b / L_main)
    n = (2 ** log2li(b + 0.1) - 2 ** log2li(b)) / 2.0
    out.append((b, n, p_missing(lam), math.exp(-lam)))
  return out
def band_w(u, b, n):
  """band-average probability that a relation of norm size u contains a given orbit of the band: n_b w_b = kappa(u, b/69)
  ln((b + 0.1)/b), the per-band count of miss_side, so sum_I w_I m_I is the recurrence's missing count exactly."""
  return kappa(u, b / 69.0) * math.log((b + 0.1) / b) / n
def miss_shared(Pm):
  """shared missing-set quantities (DESCENT.md): (a) log2 of the probability that two missing-log child visits are the same
  ideal, sum_I m_I wh_I^2 with wh_I = max over u of BAND_MAX w_I(u) / E Z(u) (the largest share of a relation's missing children
  one ideal can take; E Z(u) = sum_I w_I(u) m_I = miss_side); (b) log2 of Var S_0 / (ZERO_EPS E S_0)^2 at the worst u with the
  within-band maximum weights, S_0 = sum_I w_I 1[I never occurs] (Chebyshev, pairwise negatively correlated occupancies).
  Returns the worse side of each."""
  pair = -1e9; cheb = -1e9
  for um, L in ((Pm.MAIN_UF, Pm.MAIN_LF), (Pm.MAIN_UG, Pm.MAIN_LG)):
    bands = band_weights(um, Pm.LOG2_RELS, L)
    wh = [0.0] * len(bands)
    for u in U_GRID:
      ws = [band_w(u, b, n) for b, n, pm, p0 in bands]
      EZ = sum(n * w * pm for (b, n, pm, p0), w in zip(bands, ws))
      ES0 = sum(n * w * p0 for (b, n, pm, p0), w in zip(bands, ws))
      VS0 = sum(n * (BAND_MAX * w) ** 2 * p0 for (b, n, pm, p0), w in zip(bands, ws))
      cheb = max(cheb, math.log2(VS0 / (ZERO_EPS * ES0) ** 2))
      for i, w in enumerate(ws):
        wh[i] = max(wh[i], BAND_MAX * w / EZ)
    pair = max(pair, math.log2(sum(n * pm * w * w for (b, n, pm, p0), w in zip(bands, wh))))
  return pair, cheb
def miss_consistency(Pm):
  """largest relative difference between sum_I w_I(u) m_I (the shared-set population) and miss_side's per-relation count over
  the plan's norm sizes and both sides ([7u] checks it is below 1e-9: one ideal universe)."""
  worst = 0.0
  for um, L in ((Pm.MAIN_UF, Pm.MAIN_LF), (Pm.MAIN_UG, Pm.MAIN_LG)):
    bands = band_weights(um, Pm.LOG2_RELS, L)
    for u in U_GRID:
      EZ = sum(n * band_w(u, b, n) * pm for b, n, pm, p0 in bands)
      ms = miss_side(int(round(4 * u)), Pm.LOG2_RELS, um, L)
      worst = max(worst, abs(EZ - ms) / ms)
  return worst
def missing(Pm, nf, ng):
  """missing-log ideals per relation: f side at nf bits, g side at ng bits."""
  f = lambda m, um, L: miss_side(math.ceil(4 * m / 69.0), Pm.LOG2_RELS, um, L)
  return (f(nf, Pm.MAIN_UF, Pm.MAIN_LF) if nf else 0.0) + (f(ng, Pm.MAIN_UG, Pm.MAIN_LG) if ng else 0.0)
KA_MAX = 16
def _conv(a, b):
  """out[k] = sum_i a[i] b[k - i] over the nonzero supports, one dot product per k."""
  from operator import mul as _mul
  la, lb = len(a), len(b)
  out = [0.0] * (la + lb - 1)
  A = la
  while A > 0 and a[A - 1] == 0.0:
    A -= 1
  B0 = 0
  while B0 < lb and b[B0] == 0.0:
    B0 += 1
  B1 = lb
  while B1 > B0 and b[B1 - 1] == 0.0:
    B1 -= 1
  if A == 0 or B0 >= B1:
    return out
  bb = b[B0:B1][::-1]; L = B1 - B0
  for k in range(B0, A + B1 - 1):
    kk = k - B0
    i0 = max(0, kk - L + 1); i1 = min(kk, A - 1)
    out[k] = sum(map(_mul, a[i0:i1 + 1], bb[L - 1 - kk + i0:L - kk + i1]))
  return out
@functools.lru_cache(None)
def window_rows(n, base, iA, mesh, kmax=KA_MAX, kB=1):
  """rows[(jA, jB)]: measure/(jA! jB!) of jA A-type and jB B-type interval products by index sum; h = (n - base)/mesh; jB rows only when the B window is used."""
  h = (n - base) / mesh
  wA = [math.log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
  wB = [math.log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
  PA = [[1.0]]
  for j in range(1, kmax + 1):
    PA.append([x / j for x in _conv(PA[-1], wA)])
  PB = [[1.0], wB]
  rows = {}
  for jA in range(kmax + 1):
    for jB in range(2 if kB else 1):
      rows[(jA, jB)] = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
  return rows, h
@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha, trunc, kmax=KA_MAX, kB=1):
  """lower sums L[jA][jB] (cumulative in jA) at integer m."""
  rows, h = window_rows(n, base, iA, mesh, kmax, kB)
  out = [[0.0, 0.0] for _ in range(kmax + 1)]
  for (jA, jB), row in rows.items():
    j = jA + jB
    imax = len(row) - 1
    ifeas = min(imax, int(math.floor((m - j * h - j * base) / h))) if h > 0 else imax
    if ifeas < 0:
      continue
    imin = 0 if not trunc else max(0, ifeas - int(trunc * base / h) - 1)
    tab, ntab, ceil = RHO_TAB, len(RHO_TAB), math.ceil
    k0 = (m - j * base + alpha) / (base * RSTEP); hk = h / (base * RSTEP)
    s = 0.0
    for i in range(imin, ifeas + 1):
      weight = row[i]
      if weight:
        k = ceil(k0 - i * hk)
        s += weight * (1.0 if k <= 0 else (tab[k] if k < ntab else 0.0))
    out[jA][jB] = s
  cum = [[0.0, 0.0] for _ in range(kmax + 1)]
  for jB in range(2):
    acc = 0.0
    for jA in range(kmax + 1):
      acc += out[jA][jB]
      cum[jA][jB] = acc
  return tuple(tuple(x) for x in cum)
def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=3.0):
  """P(all primes <= 2^n, <= kA in (2^base, 2^nA], <= kB in (2^nA, 2^n]) lower sum."""
  m = math.ceil(m)
  if n <= base:
    return rho_lo((m + alpha) / base)
  cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc, kA, 1 if kB else 0)
  return cum[kA][0] + (cum[kA][1] if kB else 0.0)
def surv_upper(m, base, T):
  """P(product of the prime factors above 2^base of an m-bit integer <= 2^T) <= (T/base) rho((m - T)/base) for base <= T < m."""
  if T <= 0:
    return 0.0
  if T >= m:
    return 1.0
  return min(1.0, max(1.0, T / base) * rho_hi((math.floor(m) - T) / base))
MONT5 = 367.76               # measured five-limb Montgomery multiplication (ecm12.c, KERNELS.md)
def ecm_instr(fbits, nbits, success):
  """instructions of the ecm16 schedule finding every prime of <= fbits bits in an nbits-bit integer with probability >= success per prime."""
  import ecm16
  return ecm16.ecm_instr(fbits, nbits, success)
def log2li(l):
  x = l * LN2
  return l - math.log2(x) + math.log2(1 + 1 / x + 2 / x / x)
CELLS = {4:13315520,6:1573359464,8:47982017060,10:687193400292,12:6073556677752,14:38414493101628,
     16:190058478539564,20:2754504316057176,24:24505531366261200,28:155610653675916932,
     32:771909602591773000,40:11221331971811657968,48:99995402971660463340,
     56:635616373364584622852,64:3155035383949694198996,80:45900097950907999286464,
     96:409193506988359583551012,112:2601677688965738307448548,128:12916156661164168796080328,160:187942845067274337246511324,192:1675659512501707896389090868,
     224:10654623198171176305341831664,256:52897551248440914796084637156,288:217401149762189752879411886484,
     320:769748406668707279529221545272,384:6863095027560050329435265332068}
class P:
  """the parameter set (proof.md / DESCENT.md)."""
  LB = 69.0; LX = 69.0; PMIN = 8.0; SLACK = 96.0
  ALPHA_F = ALPHA_F_CHARGED; ALPHA_G = 1.0; LOG2U6 = 6 * 63.714   # Murphy alpha charged above the p < 50,000 sums (5.085, 0.732; [14])
  QBINS = [(1 / 6, -26.4 + 23, 14.5 + 6), (1 / 6, -11.2 + 23, 14.5 + 6), (1 / 6, -0.1 + 23, 14.5 + 6),
       (1 / 6, 10.5 + 23, 14.5 + 6), (1 / 6, 23.2 + 23, 14.5 + 6), (1 / 6, 86.6 + 23, 14.5 + 6)]
  BOOT_BINS = [(0.5, -191.5 + 4), (0.4, -185.9 + 4), (0.1, -177.0)]   # sqrt(U)-weighted boot basis: worse-radius p50/p90, ceiled to 1 decimal, + 4; the top bin covers the sample maximum (-178.33) by 1.3 bits; 9,216 exact norms on 48 certified lattices (DESCENT.md, [11])
  QBINS_LO = [(-140.05 - 23, -54.56 - 6), (-26.4 - 23, -54.56 - 6), (-11.2 - 23, -54.56 - 6),
        (-0.1 - 23, -54.56 - 6), (10.5 - 23, -54.56 - 6), (23.2 - 23, -54.56 - 6)]
  LOG2_RELS = 63.80; MAIN_LF = 69.25; MAIN_LG = 69.0; MAIN_UF = 740.5 / MAIN_LF; MAIN_UG = 498.6 / MAIN_LG   # pool estimate - 2 SE (floor) at theta 633/2000 and at or below its relation-weighted mean norms; [4] checks both
  THIN_MAX = 90               # nodes with q <= 2^95 discard relations touching a missing log
  UDUP = 0.25                 # distinct relations per raw relation (unit duplicates; stratified relation-conditioned minimum 0.283, DESCENT.md)
  MIN_RELS = 16.0; PSI_MAX = 0.8; PHI_FLOOR = 0.02
  TABLE_SLACK = 2.0 ** -16    # relative margin on the boot failure for the realized missing-log table (DESCENT.md, Recurrence and success)
  RET_FLOOR = 0.8592          # retention floor of every sieve node (retain16.retention; part=losses recomputes the exact per-node values)
  COUNTED = 0.8592             # base-primitive x ordinary fractions (retain16.counted_fraction; part=losses), the boot's counted fraction
  NB_R1 = 0.5; NB_R2 = 1.0    # negative-binomial shapes: windows of <= 1 / >= 2 primes (toy fits 0.60 / 1.86-12.2)
  BOUND_CELL = 2000.0         # certified per-cell norm bounds, both sides (84 forms; ~1500 measured, KERNELS.md)
  HITS_PER_CELL = 2 * math.log(69.0 / 8.0)
  SIEVE_HIT = 850.0; SCAN = 10.0; ZERO = 2.0   # separate f and g score arrays (757 enumeration >= 2x the largest value measured at R >= 192 + 84 transfer + 9 apply; KERNELS.md)
  PAIR_SETUP = 2.0 ** 25.1; LLL_NODE = 2.0 ** 36; NORMS = 2.0 ** 26; TDIV = 2.0 ** 24.5; IDEAL = 2.0 ** 34; PRIME = 2.0 ** 24
  ECM_SUCCESS = 0.999
  RADII = (4, 6, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 288, 320, 384)
  MESH = 384; TRUNC = 0.0
  LOG2P6 = 6 * 380.99700
  LMAX = 225; BIN = 5
def pairs_count(Pm):
  return 2.0 * 2 ** log2li(Pm.LX)
def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R, thin=None):
  """own cost (instructions) and statistics of one node plan;"""
  pairs = pairs_count(Pm)
  cfmax = max(b[1] for b in Pm.QBINS); cgmax = max(b[2] for b in Pm.QBINS)
  S = math.log2(CELLS[R]); c = (l_hi + S) / 12
  h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
  nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
  T = kA * nA + kB * n + Pm.SLACK
  c_lo = (l_lo + S) / 12
  pr = 0.0; sv = 0.0
  for (w, cf, cg), (cf_lo, cg_lo) in zip(Pm.QBINS, Pm.QBINS_LO):
    nf = 72 * c + cf - (l_lo if side == 'f' else 0)
    ng = 12 * c + Pm.LOG2U6 + cg - (l_lo if side == 'g' else 0)
    pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
    pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
    pr += w * pf * pg
    nf_lo = 72 * c_lo + cf_lo - (l_hi if side == 'f' else 0)
    ng_lo = 12 * c_lo + Pm.LOG2U6 + cg_lo - (l_hi if side == 'g' else 0)
    sv += w * surv_upper(nf_lo, Pm.LX, T) * surv_upper(ng_lo, Pm.LX, T)
  if pr <= 0:
    return None
  nf_max = 72 * c + cfmax - (l_lo if side == 'f' else 0)
  ng_max = 12 * c + Pm.LOG2U6 + cgmax - (l_lo if side == 'g' else 0)
  det = 1.0 - (1 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 20.0 + 2 * (kA + kB))
  miss = missing(Pm, nf_max, ng_max)
  if thin is None or l_lo <= Pm.THIN_MAX:
    thin = thin or l_lo <= Pm.THIN_MAX
  mu = CELLS[R] * pr * Pm.UDUP * det * (math.exp(-miss) if thin else 1.0)
  if mu < Pm.MIN_RELS:
    return None
  surv_cost = (Pm.NORMS + Pm.TDIV + ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS) + ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS)
         + 2 * (ecm_instr(nA, T, Pm.ECM_SUCCESS) + Pm.PRIME) + Pm.IDEAL)
  survivors = CELLS[R] * sv + CELLS[R] * pr
  sieve = CELLS[R] * (Pm.HITS_PER_CELL * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
  trials = CELLS[R] * surv_cost
  if trials < sieve + survivors * surv_cost:
    own = Pm.LLL_NODE + trials; how = 'ecm'; survivors = CELLS[R]
    ret, loss_f, loss_g = 1.0, 0.0, 0.0
  else:
    own = Pm.LLL_NODE + sieve + survivors * surv_cost; how = 'sieve'
    ret, loss_f, loss_g = Pm.RET_FLOOR, None, None
    mu *= ret                 # relation-conditioned retention floor (DESCENT.md, retain16.py; part=losses recomputes every node's exact value)
    if mu < Pm.MIN_RELS:
      return None
  return dict(own=own, n=n, iA=iA, nA=nA, kA=kA, kB=kB, R=R, S=S, mu=mu, how=how, miss=miss, ret=ret, loss_f=loss_f, loss_g=loss_g,
        miss_children=0.0 if thin else miss, thin=bool(thin), pr_log2=math.log2(pr),
        surv_log2=math.log2(survivors), own_log2=math.log2(own), nf=nf_max, ng=ng_max, det=det, T=T,
        sieve_log2=math.log2(sieve), survcost_log2=math.log2(survivors * surv_cost))
def bins_of(Pm):
  bins = [(int(Pm.LB), int(Pm.LB) + Pm.BIN + 1)]
  l = int(Pm.LB) + Pm.BIN + 1
  while l < Pm.LMAX:
    bins.append((l, l + Pm.BIN)); l += Pm.BIN
  return bins
def phi_nb(mean, Pm, window):
  """P(no distinct relation) for a negative-binomial count of the given mean."""
  r = Pm.NB_R1 if window <= 1 else Pm.NB_R2
  return max((1.0 + mean / r) ** (-r), Pm.PHI_FLOOR)
def evaluate(Pm, PLAN, BOOT):
  """evaluate the shipped plan: per-bin W (instructions), phi, details, and the two-log expected work;"""
  shared = ('LX', 'MESH', 'TRUNC', 'ALPHA_F', 'ALPHA_G')
  if any(getattr(Pm, name) != getattr(P, name) for name in shared):
    raise ValueError('descent and retention parameters must agree')
  bins = bins_of(Pm)
  W = {}; PHI = {}; DET = {}
  def wmax(side, n):
    return max([W[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
  def phimax(side, n):
    return max([PHI[(b, side)] for b in bins if b[0] < n - 1e-9] + [0.0])
  def miss_terms():
    """failure probability and work of a missing-log child Q (prime special-q node in its own lattice L_Q): the larger base-bin values over both sides (DESCENT.md)."""
    base = [k for k in W if k[0][0] == bins[0][0]]
    return (max([PHI[k] for k in base] + [0.0]), max([W[k] for k in base] + [0.0]))
  for b in bins:
    for side in ('f', 'g'):
      n, iA, kA, kB, R, thin = PLAN[(b[0], side)]
      d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R, bool(thin))
      if d is None:
        raise ValueError('infeasible plan entry %s %s' % (b, side))
      nA = d['nA']
      pm, wm = miss_terms()
      psi = kA * (phimax('f', nA) + phimax('g', nA)) + kB * (phimax('f', n) + phimax('g', n)) + d['miss_children'] * pm
      if psi > Pm.PSI_MAX:
        raise ValueError('psi too large at %s %s' % (b, side))
      retry = 1.0 / (1.0 - psi)
      sub = kA * (wmax('f', nA) + wmax('g', nA)) + kB * (wmax('f', n) + wmax('g', n)) + d['miss_children'] * wm
      W[(b, side)] = d['own'] + retry * sub
      window = (kA if nA > Pm.LX + 1e-9 else 0) + (kB if n > nA + 1e-9 else 0)
      PHI[(b, side)] = phi_nb(d['mu'] * (1 - psi), Pm, window)
      DET[(b, side)] = dict(d, psi=psi, retry=retry, sub_log2=math.log2(sub) if sub > 0 else None,
                 W_log2=math.log2(W[(b, side)]), phi=PHI[(b, side)])
  n1, iA, kA, Rb, nrels = BOOT   # nrels distinct boot relations are kept and tried in turn
  if n1 > Pm.LMAX:
    raise ValueError('boot window above the last bin')
  Sb = math.log2(CELLS[Rb]); cb = (Pm.LOG2P6 + Sb) / 12
  h = (n1 - Pm.LX) / Pm.MESH; nA = Pm.LX + iA * h
  m_max = 12 * cb + Pm.LOG2U6 + max(x[1] for x in Pm.BOOT_BINS)
  pr = sum(w * ss_lower2(12 * cb + Pm.LOG2U6 + cg, Pm.LX, float(n1), iA, kA, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
       for w, cg in Pm.BOOT_BINS)
  det = 1.0 - (1 - Pm.ECM_SUCCESS) * (m_max / 20.0 + kA + 1)
  trials = nrels / (pr * det * Pm.UDUP * Pm.COUNTED)
  balls = 1.0 + trials / CELLS[Rb]
  trial = Pm.NORMS + ecm_instr(nA, m_max, Pm.ECM_SUCCESS) + Pm.PRIME + Pm.IDEAL
  pm, wm = miss_terms(); mb = missing(Pm, 0, m_max)
  psi = kA * phimax('g', nA) + phimax('g', float(n1)) + mb * pm
  if psi > Pm.PSI_MAX:
    raise ValueError('boot psi too large')
  retry = 1 / (1 - psi)
  sub = kA * wmax('g', nA) + wmax('g', float(n1)) + mb * wm
  one_log = Pm.LLL_NODE * balls + trials * trial + retry * sub
  phi_boot = psi ** nrels * (1 + Pm.TABLE_SLACK)   # all nrels kept relations fail (each with probability <= psi conditional on the main-phase table; the relative margin covers the table's sampling fluctuation, DESCENT.md)
  boot = dict(n1=n1, nA=nA, kA=kA, R=Rb, nrels=nrels, S=Sb, m=m_max, pr_log2=math.log2(pr), trials_log2=math.log2(trials),
        trial_log2=math.log2(trial), balls=balls, miss=mb, own_log2=math.log2(trials * trial), sub_log2=math.log2(sub),
        psi=psi, retry=retry, phi=phi_boot, one_log_log2=math.log2(one_log))
  return dict(W=W, PHI=PHI, DET=DET, bins=bins, boot=boot, two_log_units=math.log2(2 * one_log / UNIT),
        one_log_instr=one_log)
PLAN16 = {(int(k[:-1]), k[-1]): tuple(float(x) if "." in x else int(x) for x in v.split(","))   # (bin, side): (n, iA x 8, kA, kB, R, discard); searched at this revision's inputs (DESCENT.md)
          for k, v in (e.split(":") for e in (
 "69f:69,368,0,0,288,1 69g:69,384,0,0,320,1 75f:74,268,1,0,256,1 75g:75,336,1,0,256,1 80f:75,328,1,0,320,0 80g:75,260,2,0,192,0 85f:80,188,2,0,192,1 "
 "85g:75,364,2,0,224,1 90f:77,276,2,0,256,1 90g:76,320,2,0,320,0 95f:81,192,3,0,288,1 95g:80,192,2,1,256,1 100f:85,88,1,1,320,1 100g:85,132,2,1,288,1 "
 "105f:84,136,2,1,320,1 105g:85,140,3,1,320,0 110f:90,148,2,1,288,1 110g:94,152,2,1,320,1 115f:92,180,2,1,320,1 115g:100,196,2,1,288,1 120f:94,228,2,1,320,0 "
 "120g:94,372,3,0,384,0 125f:100,188,2,1,384,1 125g:100,192,3,1,384,1 130f:100,192,3,1,320,1 130g:100,192,4,1,384,1 135f:105,212,3,1,288,1 135g:110,236,3,1,320,0 "
 "140f:103,272,3,1,320,0 140g:110,276,3,1,384,0 145f:112,232,3,1,320,1 145g:114,260,4,1,288,1 150f:110,272,3,1,384,0 150g:110,336,4,1,384,0 155f:114,220,4,1,288,0 "
 "155g:123,244,4,1,384,1 160f:115,216,4,1,384,0 160g:128,236,4,1,384,0 165f:124,196,4,1,384,0 165g:130,268,4,1,384,0 170f:124,224,4,1,384,0 170g:135,288,4,1,384,1 "
 "175f:136,240,4,1,288,1 175g:150,240,4,1,384,0 180f:130,288,4,1,320,0 180g:130,272,6,1,320,0 185f:140,240,4,1,384,0 185g:140,244,6,1,320,0 190f:145,240,4,1,384,0 "
 "190g:157,188,6,1,384,0 195f:150,240,4,1,384,0 195g:161,316,4,1,384,0 200f:150,288,4,1,320,0 200g:150,280,6,1,288,0 205f:165,240,4,1,384,1 205g:159,244,6,1,384,0 "
 "210f:160,288,4,1,320,0 210g:172,220,6,1,384,0 215f:165,288,4,1,320,0 215g:170,288,6,1,224,0 220f:170,288,4,1,384,1 220g:170,288,6,1,288,0 ").split())}
BOOT16 = (210, 356, 13, 4, 16)
def report(res):
  """compact per-bin table for the experiment report."""
  rows = []
  for b in res['bins']:
    for side in ('f', 'g'):
      d = res['DET'][(b, side)]
      rows.append(dict(bin=b[0], side=side, n=int(d['n']), nA=round(d['nA'], 1), kA=d['kA'], kB=d['kB'], R=d['R'], ret=round(d['ret'], 4),
               loss_f=d['loss_f'], loss_g=d['loss_g'], how=d['how'], mu=round(d['mu'], 2), own_log2=round(d['own_log2'], 3),
               W_log2=round(d['W_log2'], 3), psi=round(d['psi'], 4), phi=round(d['phi'], 6), nf=round(d['nf'], 1), ng=round(d['ng'], 1), T=round(d['T'], 1)))
  return rows
