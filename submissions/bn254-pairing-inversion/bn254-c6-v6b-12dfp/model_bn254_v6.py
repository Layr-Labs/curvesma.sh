"""
model_bn254_v6.py
Pure-Python-stdlib verification & ledger engine for BN254 v6 (C6-m184_w216-rec):
  - Tower: C6-m184_w216-rec (content-reduced reciprocal rank-2 CM-sextic, order-6 tau = Frob^10,
  exact rational sigma(w) == w twist invariance)
  - 100% Pure Baseline Scalar AArch64 (-mgeneral-regs-only, simd = 0, zero FEAT_LSE):
      * Band 1 [2^14, 2^44): bn254_apply_8b_u16_u4 (4.25 flat / 6.25 own64) + bn254_b1_early4_find
      * Band 2 [2^44, 2^54): 8-lane 64-bit SWAR sb_u7 per-SB tight q_min_j gate (W_SB=256, max_sb=5<=8,
    C_SWAR_PROBE=18.0, C_SB_ENTER=14.0, C_CAND=13.0, C_HIT_HI=9.0 flat / 14.5 own64)
      * Band 3 [2^54, 2^59): sigma-orbit-folded 12D Optimal Bounding Ellipsoid E* Fincke-Pohst walk
    (1/3 as many ideals; exact VR_12=2.238988x provable ratio over active R={u<=1, v+kap*u<=1};
    MULT_E12=5.85 > VR_12/f_FD=5.7543; C_FP12_NODE=26.0, C_FP12_LEAF_GATE=14.0 compiled 12 Ir,
    C_FP12_LEAF_FD=36.0 on MULT_FD_ACT=2.62 > 1/f_FD=2.5700 active leaves; C_BKT_B3=16.0)
      * Scan+Zero: uint16_t 4-cell/word 64-bit scalar SWAR (C_SCAN_ZERO=2.15, step=0.25)
      * Stage B: 100% score-only (ZERO sieve prime recording!); C_SB=80,000 strips exact 2^14-smooth
    part m of N_g via 380 64-bit word-product chunks and tests O = 16*log2(N_g) - score - 16*log2(m) <= 16*S_B
  - Primary model: Full de Bruijn local-density correction rho_dB (live norm-histogram-weighted Retention
  & Stage B) + fixed-instance H3 4-parameter hard-cap Z[omega] LA (W_ROW=161.0, NU_SW_CAP=0.00160,
  NU_L_CAP=0.00098, K5_CAP=12.00, TAB_R_CAP=0.85, lazy=0) + un-discounted Schirokauer maps +
  live 62-row LB=59 pure QBINS/BOOT_BINS 12D descent + strict FIX-U + exact rational gate K.
"""
import json, math, random
from math import log, log2, exp, sqrt, ceil, floor, pi, log1p, gamma
from fractions import Fraction as F
from decimal import Decimal, localcontext
from collections import Counter
from array import array
from pathlib import Path
import descent_bn254_v6 as DESC
LN2 = math.log(2.0)
GAMMA = 0.5772156649015329
UNIT = 20000.0
U_BN254 = 4965661367192848881
P_BN254 = 36 * U_BN254**4 + 36 * U_BN254**3 + 24 * U_BN254**2 + 6 * U_BN254 + 1
R_BN254 = 36 * U_BN254**4 + 36 * U_BN254**3 + 18 * U_BN254**2 + 6 * U_BN254 + 1
H_COEFFS = [3, -4, 0, 2, 1, -2, 1]
W_COEFFS = [0, -2, -1, 0, 1, -1]
SIGMA_NUM = [9, -4, -5, 1, 3, -2]
ZK_NUM = [
  [3,  0,  0,  0,  0,  0],
  [0,  3,  0,  0,  0,  0],
  [9, -4, -5,  1,  3, -2],
  [0,  2,  1, -2,  0,  1],
  [6,  1, -1, -1,  3, -1],
  [-6, 6,  6,  0, -3,  3],
]
F_REC_COEFFS = {
  8: [3, -4, -2, 0, 2, -2],
  6: [0, 6, 3, 0, -3, 3],
  4: [-12, 0, 0, 0, 0, 0],
  2: [12, -6, -3, 0, 3, -3],
  0: [-5, 4, 2, 0, -2, 2],
}
CONTENT_BITS = 12.0 + 6.0 * log2(3.0)
EXACT_PRIM_T4 = 0.683542
ALPHA_F_RAW = -20.010 + 0.08
ALPHA_G_RAW = -0.119 + 0.08
_DATA = None
def load_data():
  global _DATA
  if _DATA is None:
    _DATA = json.loads((Path(__file__).resolve().parent / "c6_data.json").read_text())
  return _DATA
# Verified pure baseline scalar AArch64 kernel prices (-mgeneral-regs-only, simd=0, zero FEAT_LSE)
C_HIT_LO = 4.25       # bn254_apply_8b_u16_u4: 17 instr / 4 hits = 4.25 Ir/hit (6.25 if own64 buffered)
C_SWAR_PROBE = 18.0   # bn254_swar_sb_u7_probe_batch: 17.75 <= 18.0 Ir per (p, sign, b-row) covering <= 8 SBs
C_SB_ENTER = 14.0     # Band-2 entered SB directory load + scan setup + terminator
C_CAND = 13.0         # bn254_b2_scan_direct_u16: 13 Ir/candidate check
C_HIT_HI = 9.0        # bn254_b2_scan_direct_u16: +7 <= 9.0 Ir/hit direct single-owner (14.5 if own64 buffered)
C_FP12_NODE = 26.0     # Band-3 12D Fincke-Pohst internal node (12D shift + fsqrt + bounds)
C_FP12_LEAF_GATE = 14.0 # Band-3 E*->R finite-diff gate u<=1 & (1-cb)(1-u)<=nr (compiled 12<=14 Ir, simd=0)
C_FP12_LEAF_FD = 36.0  # Band-3 active R-leaf exact integer re-check + 6-facet Voronoi unit-FD + sigma canon
F_FD_UNITS = 0.3891    # Exact Voronoi unit-FD angular fraction inside Q_a(a) <= 1 (1/0.3891 = 2.5700)
MULT_FD_ACT = 2.62     # Conservative active-R leaf-per-canonical-hit multiplier (+1.9% above exact 2.5700)
MULT_E12 = 5.85        # Conservative E* leaf-per-canonical-hit multiplier (+1.7% above exact VR_12/f_FD=5.7543)
C_LLL12 = 2.0 ** 22   # Generous 12D LLL basis reduction budget per sigma-orbit (1/3 as many ideals!)
C_BKT_B3 = 16.0       # Band-3 8B log_p-grouped bucket: 6 write + 5 xfer + 4.25 apply = 15.25 <= 16.0 Ir/hit
C_SCAN_ZERO = 2.15    # bn254_scan_zero_u16: 1.91 <= 2.15 Ir/cell (8x unrolled 4-cell/word 64-bit scalar SWAR)
C_LIST_A = 35.0       # 3-pass pure-C a-table build per entry
C_LBUCKET = 4.0       # Conservative per-SB build cost (1-byte sb_u7 + offset)
C_UNIT_5 = 7.8        # Unit-entry Z[omega] cost per 5-lane BW block
C_GEN_1 = 33.0        # General non-unit Z[omega] entry cost per 5-lane BW block
# Score-only Stage B (ZERO Band-1/2/3 prime recording!):
C_NORMS = 27000.0     # Analytical scalar norm evaluation (12-coord degree-8/2 f_rec, g_rec norms)
C_SB = 80000.0        # Exact P_14 score-only strip (380 64-bit prime-product chunks on 8-limb N_g ~18k <= 80k ops!)
C_TDIV = 2.0 ** 23.45
C_BATCH = 2.0 ** 22.12
C_BATCH2 = 2.0 ** 23.13
C_GBATCH = 2.0 ** 23.0
LIST_NS = 256
LI_PAD = 1.001
Q_F_SAFETY = 100.0
KMAX = 4
FFT_A, FFT_B, FFT_XFER, FFT_LOCAL, FFT_COPY, FFT_ROOT, FFT_LEAF = 8.5, 1200.0, 504.0, 67108864, 3.0, 6.0, 32768.0
NTT3 = (69.75, 119.0, 90.0, 151.5, 371.0, 139.25, 271.5)
NTT3_XFER = 3.0 * FFT_XFER + 128.0
SPLIT_Q, SPLIT_REC, SPLIT_C, SPLIT_IDEAL, NORM_BITS = 1099511627776.0, 56, 33024.0, 67108864.0, 576
SPLIT_CAP, SPLIT_RES = 112.0, 84.0
# Fixed-instance H3 4-parameter post-filtering hard caps (deterministic pre-merge W_avg <= 161.0 guard, lazy=0!)
W_ROW = 161.0         # Enforced deterministically by pre-merge weight guard W_target = 160
NU_SW_CAP = 0.00160   # +73% above U'=2^22 mean (0.000923), +55% above U'=2^21 mean (0.001031)
NU_L_CAP = 0.00098    # +72% above U'=2^22 mean (0.000570), +56% above U'=2^21 mean (0.000630)
K5_CAP = 12.00        # +14.3% above U'=2^22 mean (10.497)
TAB_R_CAP = 0.85      # +14.7% above U'=2^22 mean (0.7411)
C_PROJ = 28.0
C_FULL = 97.0
C_Z64 = 40.0
C_RED9 = 270.0
BW_SIGMA = 0.0078125
PREP_OWNER = 2.0 ** 80
PREP_GUARD = 2.0 ** 80
PREP_JOINT = 2.0 ** 80
LA_RESERVE = 2.0 ** 78
G2_OUTPUT_BOUND = 69456896.0
OP = dict(
  lf=61.0, lg=60.0, exc=0.82, W=160,
  p_min=14, p_cut=44.0, p_mid=54.0, lx2=59.0,
  slack=114.0, slack_b=119.0, w_sb=256, c_tail_z6=420.0,
  margin=0.0666, work_allow=1.0033,
)
def fp12_gauss_nph(h):
  """Exact parameter-free 12D Gaussian-volume-ratio internal nodes per ellipsoid leaf at density h = V_12/(2p)."""
  v12 = (pi ** 6.0) / gamma(7.0)
  tot = 1.0
  for d in range(1, 12):
    vd = (pi ** (0.5 * d)) / gamma(0.5 * d + 1.0)
    tot += (vd / (v12 ** (d / 12.0))) * (h ** (d / 12.0))
  return tot / max(h, 1e-12)
# 1. EXACT DICKMAN RHO & DE BRUIJN RHO_DB (PURE STDLIB)
_RHO_COEFS = None
def _build_rho_coefs(kmax=33, nt=140):
  with localcontext() as ctx:
    ctx.prec = 70
    zero, one = Decimal(0), Decimal(1)
    coefs = [None, [one] + [zero] * (nt - 1)]
    for k in range(2, kmax + 1):
      prev, c = coefs[k - 1], [zero] * nt
      for m in range(nt - 1):
        c[m + 1] = (prev[m] + m * c[m]) / (k * (m + 1))
      c[0] = prev[0] - sum(c[1:])
      coefs.append(c)
    out = [None]
    for c in coefs[1:]:
      big = max(abs(x) for x in c)
      while len(c) > 1 and abs(c[-1]) < big * Decimal("1e-18"):
        c = c[:-1]
      out.append([float(x) for x in c])
    return out
def rho(u):
  global _RHO_COEFS
  if u <= 1.0:
    return 1.0
  if _RHO_COEFS is None:
    _RHO_COEFS = _build_rho_coefs()
  k = ceil(u)
  if k >= len(_RHO_COEFS):
    return 0.0
  x = k - u
  v = 0.0
  for c in reversed(_RHO_COEFS[k]):
    v = v * x + c
  return max(0.0, v)
def rho_db(u, log2_N):
  if u <= 1.0:
    return 1.0
  return max(1e-300, rho(u) - GAMMA * rho(u - 1.0) / max(1.0, log2_N * LN2))
def kappa_rho(u, t):
  ru = rho(u)
  return rho(u - t) / ru if ru > 0.0 else 0.0
# 2. LIVE EULER-PRODUCT RETENTION & SURVIVOR MODEL
def _primes_upto(n):
  sieve = bytearray(b"\x01") * (n + 1)
  sieve[:2] = b"\x00\x00"
  for p in range(2, int(sqrt(n)) + 1):
    if sieve[p]:
      sieve[p * p :: p] = b"\x00" * ((n - p * p) // p + 1)
  return [p for p in range(2, n + 1) if sieve[p]]
def _series_coeffs(weights, normalizer, length):
  terms = [(j, j * v) for j, v in enumerate(weights) if j and abs(v) >= 1e-24]
  out = [exp(-normalizer)] + [0.0] * length
  for n in range(1, length + 1):
    val = 0.0
    for j, w in terms:
      if j > n: break
      val += w * out[n - j]
    out[n] = val / n
  return out
class Retention:
  def __init__(self, B_bits, X_bits, prime_min=2**14, step=0.25, max_bits=600, tilt=0.05, explicit_bits=14, power_terms=8, odd=1):
    self.B, self.X, self.minimum, self.step = B_bits, X_bits, prime_min, step
    self.tilt, self.length = tilt, ceil(max_bits / step)
    self.names = ("small", "sieve_geo", "large_geo", "sieve_sf", "large_sf")
    weights = {k: [0.0] * (self.length + 1) for k in self.names}
    totals = Counter()
    cutoff = min(float(explicit_bits), B_bits)
    for p in _primes_upto(int(2 ** cutoff) - 1)[odd:]:
      b = log2(p)
      bi = round(b / step)
      q = exp(tilt * bi * step - LN2 * b)
      groups = [("small", False)] if p < prime_min else [(("sieve" if b < X_bits else "large") + "_geo", False), (("sieve" if b < X_bits else "large") + "_sf", True)]
      for name, sf in groups:
        totals[name] += log1p(q) if sf else -log1p(-q)
        qk = q
        for k in range(1, self.length + 1):
          idx = k * bi
          if idx > self.length: break
          weights[name][idx] += (-1 if sf and k % 2 == 0 else 1) * qk / k
          qk *= q
          if qk < 1e-25: break
    def continuum(left, right, pref):
      if left >= right: return
      w = step / 4.0
      cells = ceil((right - left) / w)
      w = (right - left) / cells
      for i in range(cells):
        b = left + (i + 0.5) * w
        bi = round(b / step)
        for k in range(1, power_terms + 1):
          amt = exp(k * tilt * bi * step - (k - 1) * LN2 * b) * w / (k * b)
          idx = k * bi
          sign = -1 if k % 2 == 0 else 1
          grps = ((pref, amt),) if pref == "small" else ((pref + "_geo", amt), (pref + "_sf", sign * amt))
          for g, x in grps:
            totals[g] += x
            if idx <= self.length: weights[g][idx] += x
    lm = log2(prime_min)
    continuum(cutoff, min(lm, B_bits), "small")
    continuum(max(cutoff, lm), min(X_bits, B_bits), "sieve")
    continuum(max(cutoff, X_bits), B_bits, "large")
    self.weights, self.normalizers = weights, {k: totals[k] for k in self.names}
    all_w = [sum(weights[k][j] for k in ("small", "sieve_geo", "large_geo")) for j in range(self.length + 1)]
    self.denominator = _series_coeffs(all_w, sum(totals[k] for k in ("small", "sieve_geo", "large_geo")), self.length)
    self.selection_normalizer = exp(totals["sieve_sf"] + totals["large_sf"] - totals["sieve_geo"] - totals["large_geo"])
    self._c = {}
  def co(self, k):
    if k not in self._c:
      self._c[k] = _series_coeffs(self.weights[k], self.normalizers[k], self.length)
    return self._c[k]
  def two(self, n, a, b=9**9, win=None, sf=False):
    a = min(a, n)
    if a < 0 or b < 0: return 0.0
    key = (n, b, win, sf)
    if key not in self._c:
      if key[1:] not in self._c:
        st = self.step; sm = self.co("small"); D = min(self.length, 2048); m = [0.0] * (D + 1)
        for l, v in enumerate(self.co("large_sf" if sf else "large_geo")[: min(b, D) + 1]):
          if v and (win is None or l * st <= win[2] or any(self.X * k - win[2] <= l * st <= win[1] * k + win[2] for k in range(1, win[0] + 1))):
            for t in range(D + 1 - l): m[l + t] += v * sm[t]
        self._c[key[1:]] = m
      m = self._c[key[1:]]; sv = self.co("sieve_sf" if sf else "sieve_geo"); c = array("d"); t = 0.0
      for x in range(min(n, len(m) - 1) + 1):
        t += m[x] * sv[n - x]; c.append(t)
      self._c[key] = c
    c = self._c[key]
    return c[min(a, len(c) - 1)]
  def cap(self, nb): return floor(ceil(nb) / log2(self.minimum))
  def at(self, nb, sl, sb=None):
    n = round(nb / self.step); rb = (self.cap(nb) + 1) / 16.0
    r = self.selection_normalizer * self.two(n, floor((sl - rb) / self.step), 9**9 if sb is None else floor((sb - rb - 0.125) / self.step), None, True) / self.denominator[n]
    return max(0.0, min(1.0, r))
  def survivor_upper(self, nb, sl, ex=0.0, sb=None, lb=None, kmax=4):
    n = round(nb / self.step); a = floor((sl + ex) / self.step)
    return (self.two(n, a) if sb is None else self.two(n, a, floor((sb + 0.25) / self.step), (kmax, lb, (self.cap(nb) + 2) / 16.0 + 0.25))) / self.denominator[n]
# 3. LIVE ALGEBRA + EXACT RATIONAL sigma(w) == w + MOD-r CRT FACTORIZATION
def _ptrim(a):
  while len(a) > 1 and a[-1] == 0: a.pop()
  return a
def _pmul(a, b, mod):
  c = [0] * (len(a) + len(b) - 1)
  for i, x in enumerate(a):
    if x:
      for j, y in enumerate(b):
        if y: c[i + j] = (c[i + j] + x * y) % mod
  return _ptrim(c)
def _pmod(a, b, mod):
  a = [x % mod for x in a]; b = _ptrim([x % mod for x in b])
  inv = pow(b[-1], -1, mod)
  while len(a) >= len(b) and not (len(a) == 1 and a[0] == 0):
    d = len(a) - len(b); fac = (a[-1] * inv) % mod
    for i in range(len(b)): a[d + i] = (a[d + i] - fac * b[i]) % mod
    _ptrim(a)
  return a
def _pdiv(a, b, mod):
  a = _ptrim([x % mod for x in a]); b = _ptrim([x % mod for x in b])
  q = [0] * max(1, len(a) - len(b) + 1); inv = pow(b[-1], -1, mod)
  while len(a) >= len(b) and not (len(a) == 1 and a[0] == 0):
    d = len(a) - len(b); fac = (a[-1] * inv) % mod; q[d] = fac
    for i in range(len(b)): a[d + i] = (a[d + i] - fac * b[i]) % mod
    _ptrim(a)
  return _ptrim(q), a
def _pgcd(a, b, mod):
  a = _ptrim([x % mod for x in a]); b = _ptrim([x % mod for x in b])
  while not (len(b) == 1 and b[0] == 0):
    a, b = b, _pmod(a, b, mod)
  inv = pow(a[-1], -1, mod)
  return [(x * inv) % mod for x in a]
def _ppowmod(base, e, m, mod):
  res = [1]; b = _pmod(base, m, mod)
  while e > 0:
    if e & 1: res = _pmod(_pmul(res, b, mod), m, mod)
    b = _pmod(_pmul(b, b, mod), m, mod); e >>= 1
  return res
def _ddf_degrees(f, mod):
  f = _ptrim([x % mod for x in f]); inv = pow(f[-1], -1, mod); f = [(x * inv) % mod for x in f]
  df = [(i * f[i]) % mod for i in range(1, len(f))]
  if len(_pgcd(f, df, mod)) != 1: return None
  degs = []; h = [0, 1]; rem = f[:]; d = 1
  while len(rem) - 1 >= 2 * d:
    h = _ppowmod(h, mod, rem, mod)
    hx = h[:] + [0] * max(0, 2 - len(h)); hx[1] = (hx[1] - 1) % mod
    g = _pgcd(hx, rem, mod); gd = len(g) - 1
    if gd > 0:
      degs.extend([d] * (gd // d))
      rem, _ = _pdiv(rem, g, mod)
      h = _pmod(h, rem, mod)
    d += 1
  if len(rem) - 1 > 0: degs.append(len(rem) - 1)
  return sorted(degs)
def _roots_mod(f, mod):
  rng = random.Random(12345)
  f = _ptrim([x % mod for x in f]); inv = pow(f[-1], -1, mod); f = [(x * inv) % mod for x in f]
  out, stack = [], [f]
  while stack:
    cur = stack.pop(); deg = len(cur) - 1
    if deg == 0: continue
    if deg == 1: out.append((-cur[0]) % mod); continue
    while True:
      a = rng.randrange(mod)
      pw = _ppowmod([a, 1], (mod - 1) // 2, cur, mod)
      pw[0] = (pw[0] - 1) % mod
      g = _pgcd(pw, cur, mod)
      if 0 < len(g) - 1 < deg:
        q, _ = _pdiv(cur, g, mod)
        stack.append(g); stack.append(q); break
  return sorted(out)
def _pmul_h(a, b, mod):
  r = _pmod(_pmul(a, b, mod), H_COEFFS, mod)
  return (r + [0] * 6)[:6]
def _ppow_h(a, e, mod):
  r = _ppowmod(a, e, H_COEFFS, mod)
  return (r + [0] * 6)[:6]
def _pinv_h(a, mod):
  r0, r1 = [x % mod for x in H_COEFFS], _ptrim([x % mod for x in a])
  t0, t1 = [0], [1]
  while not (len(r1) == 1 and r1[0] == 0):
    q, rem = _pdiv(r0, r1, mod)
    qt = _pmul(q, t1, mod)
    L = max(len(t0), len(qt))
    t2 = _ptrim([((t0[i] if i < len(t0) else 0) - (qt[i] if i < len(qt) else 0)) % mod for i in range(L)])
    r0, r1, t0, t1 = r1, rem, t1, t2
  inv = pow(r0[0], -1, mod)
  return ([(x * inv) % mod for x in t0] + [0] * 6)[:6]
def _norm_Kh_exact(coeffs_int):
  M = [[F(0)] * 6 for _ in range(6)]; cur = [F(c) for c in coeffs_int]
  for col in range(6):
    for r in range(6): M[r][col] = cur[r]
    top = cur[5]; nxt = [F(0)] + cur[:5]
    for j in range(6): nxt[j] -= top * H_COEFFS[j]
    cur = nxt
  det = F(1)
  for col in range(6):
    piv = next(r for r in range(col, 6) if M[r][col] != 0)
    if piv != col: M[col], M[piv] = M[piv], M[col]; det = -det
    pv = M[col][col]; det *= pv
    for r in range(col + 1, 6):
      if M[r][col] != 0:
        fac = M[r][col] / pv
        for j in range(col, 6): M[r][j] -= fac * M[col][j]
  return int(det)
def _pmul_h_q(a, b):
  res = [F(0)] * 11
  for i in range(6):
    for j in range(6): res[i + j] += a[i] * b[j]
  for k in range(10, 5, -1):
    c = res[k]
    if c:
      res[k] = F(0)
      for j in range(6): res[k - 6 + j] -= c * H_COEFFS[j]
  return res[:6]
def verify_tower_algebra_and_crt():
  """Live big-int tower algebra, exact rational sigma(w) == w invariance, and live mod-r CRT factorization!"""
  ord_ok = (pow(P_BN254, 12, R_BN254) == 1) and (pow(P_BN254, 6, R_BN254) != 1) and (pow(P_BN254, 4, R_BN254) != 1)
  y = [0, 1, 0, 0, 0, 0]
  inv3 = pow(3, -1, P_BN254)
  sig_p = [(c * inv3) % P_BN254 for c in SIGMA_NUM]
  y_p1 = _ppow_h(y, P_BN254, P_BN254)
  y_p2 = _ppow_h(y_p1, P_BN254, P_BN254)
  y_p3 = _ppow_h(y_p2, P_BN254, P_BN254)
  y_p4 = _ppow_h(y_p3, P_BN254, P_BN254)
  y_p6 = _ppow_h(_ppow_h(y_p4, P_BN254, P_BN254), P_BN254, P_BN254)
  h_irred = (y_p6 == y) and (y_p2 != y) and (y_p3 != y)
  sig_frob_ok = (y_p4 == sig_p)
  # Exact rational sigma(w) == w check in K = Q[y]/(h(y)) (foundation of sigma-orbit-folded 12D Band 3!)
  sig_q = [F(c, 3) for c in SIGMA_NUM]
  pows = [[F(1)] + [F(0)] * 5, sig_q]
  for _ in range(4): pows.append(_pmul_h_q(pows[-1], sig_q))
  w_sig = [sum(W_COEFFS[k] * pows[k][j] for k in range(6)) for j in range(6)]
  sigma_fixes_w = (w_sig == [F(c) for c in W_COEFFS])
  uw = [(U_BN254 * c) % P_BN254 for c in W_COEFFS]
  uw_half = _ppow_h(uw, (P_BN254 - 1) // 2, P_BN254)
  cur, acc = uw_half, [1, 0, 0, 0, 0, 0]
  val_p6 = val_p10 = None
  for j in range(10):
    acc = _pmul_h(acc, cur, P_BN254)
    if j == 5: val_p6 = acc[:]
    if j == 9: val_p10 = acc[:]
    cur = _ppow_h(cur, P_BN254, P_BN254)
  minus_one = [P_BN254 - 1, 0, 0, 0, 0, 0]
  g_irred = (val_p6 == minus_one)
  tau_frob10 = sig_frob_ok and (val_p10 == minus_one)
  chi = pow(P_BN254, 10, R_BN254)
  chi_expr = 36 * U_BN254**3 + 18 * U_BN254**2 + 6 * U_BN254 + 2
  chi_ok = (chi == chi_expr) and ((chi * chi - chi + 1) % R_BN254 == 0) and ((pow(chi, 3, R_BN254) + 1) % R_BN254 == 0)
  inv_uw = _pinv_h(uw, P_BN254)
  f_eval = [0] * 6; x2_pow = [1, 0, 0, 0, 0, 0]
  for d in (0, 2, 4, 6, 8):
    term = _pmul_h([c % P_BN254 for c in F_REC_COEFFS[d]], x2_pow, P_BN254)
    f_eval = [(a + b) % P_BN254 for a, b in zip(f_eval, term)]
    x2_pow = _pmul_h(x2_pow, inv_uw, P_BN254)
  f_vanishes = all(c == 0 for c in f_eval)
  nw, n_lead, n_const = abs(_norm_Kh_exact(W_COEFFS)), abs(_norm_Kh_exact(F_REC_COEFFS[8])), abs(_norm_Kh_exact(F_REC_COEFFS[0]))
  h_mod_r_degs = _ddf_degrees(H_COEFFS, R_BN254)
  rts_r = _roots_mod(H_COEFFS, R_BN254)
  g_mod_r_degs, f_mod_r_degs = [], []
  for s in rts_r:
    ws = sum(W_COEFFS[j] * pow(s, j, R_BN254) for j in range(6)) % R_BN254
    g_mod_r_degs.append(_ddf_degrees([(-1) % R_BN254, 0, (U_BN254 * ws) % R_BN254], R_BN254))
    fp = [0] * 9
    for d, cs in F_REC_COEFFS.items():
      fp[d] = sum(cs[j] * pow(s, j, R_BN254) for j in range(6)) % R_BN254
    f_mod_r_degs.append(_ddf_degrees(fp, R_BN254))
  crt_ok = (h_mod_r_degs == [1, 1, 1, 1, 1, 1] and len(rts_r) == 6
              and all(gd == [2] for gd in g_mod_r_degs)
              and all(fd == [2, 2, 2, 2] for fd in f_mod_r_degs))
  T4 = load_data()["t4_geom"]
  S = T4["SMAT"]
  def iso_err(G):
    return max(abs(sum(S[a][i] * G[a][b] * S[b][j] for a in range(6) for b in range(6)) - G[i][j]) for i in range(6) for j in range(6))
  iso_a, iso_b = iso_err(T4["GA"]), iso_err(T4["GB"])
  # Live field-map / boot-lattice algebraic identity check: gcd(r, p^6 - 1) == 1 & a + (R*a)x == (a/c0)*(c0 + c1*x)
  gcd_r_p6_ok = (math.gcd(R_BN254, pow(P_BN254, 6, R_BN254) - 1) == 1)
  c0_t = [(17 * k + 3) % P_BN254 for k in range(6)]
  c1_t = [(29 * k + 11) % P_BN254 for k in range(6)]
  a_t = [(7 * k + 5) % P_BN254 for k in range(6)]
  inv_c0 = _pinv_h(c0_t, P_BN254)
  R_boot = _pmul_h(c1_t, inv_c0, P_BN254)
  b_t = _pmul_h(R_boot, a_t, P_BN254)
  cof_t = _pmul_h(a_t, inv_c0, P_BN254)
  boot_id_ok = (_pmul_h(cof_t, c0_t, P_BN254) == a_t and _pmul_h(cof_t, c1_t, P_BN254) == b_t)
  ok = (ord_ok and h_irred and g_irred and tau_frob10 and sigma_fixes_w and chi_ok and f_vanishes
          and crt_ok and gcd_r_p6_ok and boot_id_ok and (nw == 216) and (n_lead == 729) and (n_const == 729)
          and (iso_a < 1e-10) and (iso_b < 1e-10))
  return dict(
    ok=ok, ord_r_p_12=ord_ok, h_irred_mod_p=h_irred, g_rec_irred_mod_p6=g_irred,
    tau_is_frob10=tau_frob10, sigma_fixes_w_exact=sigma_fixes_w, chi_192bit_ok=chi_ok,
    f_rec_vanishes_at_inv_uw=f_vanishes, gcd_r_p6_minus_1_is_1=gcd_r_p6_ok, boot_identity_ok=boot_id_ok,
    norm_w=nw, norm_f_lead=n_lead, norm_f_const=n_const,
    mod_r_h_degs=h_mod_r_degs, mod_r_g_degs=g_mod_r_degs, mod_r_f_degs=f_mod_r_degs,
    mod_r_crt_all_deg2_ok=crt_ok, sigma_iso_err_Ga=iso_a, sigma_iso_err_Gb=iso_b,
  )
def verify_exhaustive_fd(R2=50.0):
  """Exhaustive Fincke-Pohst enumeration of ALL nonzero v with T2(v) <= R2: asserts 0 missed & 0 dups!"""
  T4 = load_data()["t4_geom"]
  EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
  FAC6, UMATS, S = T4["FAC6"], T4["UMATS"], T4["SMAT"]
  GT2 = [[sum(2.0 * (EMB_ZK[j][a] * EMB_ZK[j][b].conjugate()).real for j in range(3)) for b in range(6)] for a in range(6)]
  def mv(M, v): return tuple(sum(M[r][c] * v[c] for c in range(6)) for r in range(6))
  def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
  S2 = mm(S, S)
  def qf(v):
    s = 0.0
    for r in range(6):
      vr = v[r]
      if vr:
        s += GT2[r][r] * vr * vr
        for c in range(r + 1, 6):
          if v[c]: s += 2.0 * GT2[r][c] * vr * v[c]
    return s
  L = [[0.0] * 6 for _ in range(6)]
  for r in range(6):
    for c in range(r + 1):
      s = GT2[r][c] - sum(L[r][k] * L[c][k] for k in range(c))
      L[r][c] = sqrt(s) if r == c else s / L[c][c]
  pts, v = [], [0] * 6
  def enum(idx, rem):
    shift = sum(L[j][idx] * v[j] for j in range(idx + 1, 6))
    diag = L[idx][idx]; rad = sqrt(max(0.0, rem))
    for x in range(ceil((-rad - shift) / diag - 1e-12), floor((rad - shift) / diag + 1e-12) + 1):
      d = diag * x + shift; nrem = rem - d * d
      if nrem >= -1e-9:
        v[idx] = x
        if idx == 0:
          if any(v): pts.append(tuple(v))
        else:
          enum(idx - 1, max(0.0, nrem))
  enum(5, R2 + 1e-9)
  def xv_of(p):
    return [2.0 * ((z := sum(EMB_ZK[j][k] * p[k] for k in range(6))).real**2 + z.imag**2) for j in range(3)]
  def is_can(p):
    if tuple(-x for x in p) < p: return False
    xv = xv_of(p)
    for fi in range(6):
      dot = sum(FAC6[fi][j] * xv[j] for j in range(3))
      if dot < -1e-9: return False
      if abs(dot) <= 1e-9:
        uv = mv(UMATS[fi], p); nuv = tuple(-x for x in uv)
        if (nuv if nuv < uv else uv) < p: return False
    for Sm in (S, S2):
      vk = mv(Sm, p); nvk = tuple(-x for x in vk)
      if (nvk if nvk < vk else vk) < p: return False
    return True
  def orbit_rep(p):
    cur = p
    for _ in range(50):
      xv = xv_of(cur); bf = -1; bd = -1e-9
      for fi in range(6):
        dot = sum(FAC6[fi][j] * xv[j] for j in range(3))
        if dot < bd: bd = dot; bf = fi
      if bf < 0: break
      cur = mv(UMATS[bf], cur)
    mq = qf(cur); comp = {cur}; q = [cur]
    while q:
      u = q.pop()
      for w in [tuple(-x for x in u), mv(S, u), mv(S2, u)] + [mv(U, u) for U in UMATS]:
        if w not in comp and abs(qf(w) - mq) <= 1e-7:
          comp.add(w); q.append(w)
    return min(comp)
  orbits = {}
  for p in pts: orbits.setdefault(orbit_rep(p), []).append(p)
  exact1 = missed = dups = 0
  for members in orbits.values():
    c = sum(1 for m in members if is_can(m))
    if c == 1: exact1 += 1
    elif c == 0: missed += 1
    else: dups += c - 1
  return dict(
    ok=(exact1 == len(orbits) and missed == 0 and dups == 0 and len(orbits) == 107),
    R2=R2, nonzero_pts=len(pts), total_orbits=len(orbits),
    exact1=exact1, missed_orbits=missed, duplicate_orbits=dups,
  )
# 4. PURE-STDLIB 6D ELLIPSOID sb_u7 PROFILE & FULL 16-ROW LEDGER
def log2_li(log2x):
  x = log2x * LN2
  return log2x - log2(x) + log2(1.0 + 1.0 / x + 2.0 / (x * x))
def log2_U_z6(lf, lg):
  uf = (2.0 ** log2_li(lf)) / 6.0
  ug = (2.0 ** log2_li(lg)) / 6.0
  deg2 = (2.0 ** (0.5 * lf) / (0.5 * lf * LN2) + 2.0 ** (0.5 * lg) / (0.5 * lg * LN2)) / 6.0
  return log2(uf + ug + deg2)
def fft_F(N):
  return FFT_A * N * (log2(max(N, 2.0)) ** 2) + FFT_B + (FFT_XFER * N if N > FFT_LOCAL else 0.0)
def ntt_F(N):
  k = max(2, ceil(log2(N))); n = 2.0 ** k; o = k % 2
  return NTT3[0]*n*k + NTT3[5 if o else 1]*n + NTT3[2]*N + NTT3[3]*k + NTT3[6 if o else 4] + (NTT3_XFER*n if n > FFT_LOCAL else 0.0)
def mul_F(N):
  return fft_F(N) if max(2, ceil(log2(N))) > 55 else min(fft_F(N), ntt_F(N))
def pbits(lo, hi):
  return 1.000028 * (2.0 ** hi - (2.0 ** lo if lo else 0)) / LN2 if lo == 0 else 1.000028 * (2.0 ** hi) / LN2
def gbatch_cost(pb, leaf=9, L=None, c_leaf=FFT_LEAF):
  SP = ceil(pb / 64.0)
  L = L if L is not None else ceil(log2((SP + 1) / leaf))
  G = L + 2; n = lambda j: leaf * (2 ** j)
  PT = sum((2 ** (L - j)) * (mul_F(2 * n(j - 1)) + FFT_COPY * 2 * n(j - 1)) for j in range(1, L + 1))
  RT = sum((2 ** (L - j + 1)) * (mul_F(n(j) + G + n(j - 1)) + FFT_COPY * (n(j - 1) + G)) for j in range(1, L + 1))
  return (PT + RT + FFT_ROOT * mul_F(2 * (max(SP, n(L)) + G))) / (2 ** L) + c_leaf
def split_parts(R, lb):
  np_ = 1.25506 * (2.0 ** lb) / (lb * LN2); h = ceil(log2(np_)); D = h + 1
  SP = 2.0 ** ceil(log2(np_ * (lb / 64.0 + 1.0)))
  L = ceil(log2(SPLIT_REC * SPLIT_Q)); G = L + 2; C = SPLIT_C
  W = R / SPLIT_Q + 1.0; B = SPLIT_REC * R * D; occ = 2 * NORM_BITS
  return dict(
    trees=B * (2 * (C * (36 * L + 2 * G) + 3 * FFT_B) + 2**15),
    roots=6 * C * (2 * SP * W * D + 136 * B) + 6 * FFT_B * B,
    gcd=B * 2.0**24, small=R * 2.0**26, output=occ * R * 2.0**20,
    ideal=occ * R * SPLIT_IDEAL, control=B * 2.0**16,
    ptree=D * SP * C + 2 * FFT_B * (2.0**h), debit=128 * (W + 2),
  )
def _det_mat(M):
  n = len(M); A = [row[:] for row in M]; det = 1.0
  for col in range(n):
    piv = max(range(col, n), key=lambda r: abs(A[r][col]))
    if piv != col: A[col], A[piv] = A[piv], A[col]; det = -det
    pv = A[col][col]; det *= pv
    for r in range(col + 1, n):
      fac = A[r][col] / pv
      for j in range(col, n): A[r][j] -= fac * A[col][j]
  return det
def _build_sb_profile_pure(L_max, w_sb, kap, nx=120, nm=24, nphase=8, n_ph_tbl=8, q_slack=1.0/126.0):
  """Pure-stdlib 6D ellipsoid integral for per-SB tight q_min_j + 8-phase 7-bit sb_u7 gate."""
  h_tot = c_tot = sb_row = tot_len = 0.0
  # Bin fj_gate and fb1 into 200 quantile bins for instant pure-Python prime integration!
  nb = 200
  wj_hist = [0.0] * (nb + 1)
  wb1_hist = [0.0] * (nb + 1)
  ph_pad = w_sb / float(n_ph_tbl)
  max_sb = 0
  for ix in range(nx):
    r2 = (ix + 0.5) / nx
    w5 = 2.5 * (r2 ** 1.5) / nx
    L = L_max * sqrt(1.0 - r2)
    tot_len += w5 * L
    du = L / nm
    for im in range(nm):
      u = ((im + 0.5) / nm - 0.5) * L
      qb = r2 + (2.0 * u / L_max) ** 2
      fa = min(1.0, max(0.0, (1.0 - qb) / kap)) ** 3
      h_tot += w5 * du * fa
      wb1_hist[min(nb, int(fa * nb))] += w5 * du
    u_lo, u_hi = -0.5 * L, 0.5 * L
    for ip in range(nphase):
      off = (ip + 0.5) / nphase * w_sb
      pos = -off; cnt = 0
      while pos < L:
        t0 = max(0.0, pos); t1 = min(L, pos + w_sb)
        if t1 > t0:
          cnt += 1
          a_u = u_lo + t0; b_u = min(u_hi, u_lo + t1)
          u_min = 0.0 if (a_u <= 0.0 <= b_u) else (a_u if abs(a_u) < abs(b_u) else b_u)
          q_min = r2 + (2.0 * u_min / L_max) ** 2
          fa_ex = min(1.0, max(0.0, (1.0 - q_min) / kap)) ** 3
          a_ug = max(u_lo, a_u - ph_pad); b_ug = min(u_hi, b_u + ph_pad)
          u_mg = 0.0 if (a_ug <= 0.0 <= b_ug) else (a_ug if abs(a_ug) < abs(b_ug) else b_ug)
          q_mg = max(0.0, r2 + (2.0 * u_mg / L_max) ** 2 - q_slack * kap)
          fa_g = min(1.0, max(0.0, (1.0 - q_mg) / kap)) ** 3
          wt = w5 / nphase
          c_tot += wt * w_sb * fa_ex; sb_row += wt
          wj_hist[min(nb, int(fa_g * nb))] += wt
        pos += w_sb
      if cnt > max_sb: max_sb = cnt
  ws = sum(wj_hist); wb1s = sum(wb1_hist)
  fj_bins = [(min(1.0, (i + 1.0) / nb), w / ws) for i, w in enumerate(wj_hist) if w > 0.0]
  fb1_bins = [(min(1.0, (i + 0.5) / nb), w / wb1s) for i, w in enumerate(wb1_hist) if w > 0.0]
  return dict(cr2=c_tot / h_tot, sb_row=sb_row, max_sb=max_sb,
        frac_reg_analytic=h_tot / tot_len, fj_bins=fj_bins, fb1_bins=fb1_bins)
def _peel_fixed_point(lf, lg, exc, nf=403.5, ng=468.2):
  lU = log2_U_z6(lf, lg)
  terms = []
  for m, lb in ((nf, lf), (ng, lg)):
    u = m / lb; steps = int(round((lb - 30.0) * 10))
    for i in range(steps):
      b = 30.0 + i * 0.1; c = b + 0.1
      k = kappa_rho(u, c / lb)
      terms.append((2.0 ** (lU - c) * 6.0 * k, k * log(c / b)))
  def run(e_val):
    def step(h):
      tot = 0.0; omh = 1.0 - h
      for lam0, wt in terms:
        lam = lam0 * e_val
        if lam < 100.0: tot += wt * exp(-lam * omh)
      return 1.0 - exp(-tot)
    h = 1.0
    for _ in range(100):
      hn = step(h)
      if abs(hn - h) < 1e-10: return hn
      h = hn
    return h
  h_val = run(exc)
  lo, hi = 0.50, 1.00
  for _ in range(24):
    mid = 0.5 * (lo + hi)
    if run(mid) > 0.75: lo = mid
    else: hi = mid
  return dict(h=h_val, exc_crit=hi, stalls=(h_val > 0.75))
_RET_CACHE = {}
_PEEL_CACHE = {}
_DESC_CACHE = {}
_PROF_CACHE = {}
def compute_ledger(track="bn254-gt-dlp", use_db=True, c_tail_z6=None, cap_scale=1.0,
                   strict_per_owner=False, bkt_b3=None):
  D = load_data()
  PQ = D["pool_quad"]; T4 = D["t4_geom"]; FILT = D["filtsim_raw"]
  lf, lg, exc = OP["lf"], OP["lg"], OP["exc"]
  p_min, p_cut, pm, lx2 = OP["p_min"], OP["p_cut"], OP["p_mid"], OP["lx2"]
  sl, sb, w_sb = OP["slack"], OP["slack_b"], OP["w_sb"]
  if c_tail_z6 is None: c_tail_z6 = OP["c_tail_z6"]
  if bkt_b3 is None: bkt_b3 = C_BKT_B3
  # 1. Live Euler-product Retention & Survivor model (histogram-weighted across ALL 200k-pool norm bins!)
  rkey = (lg, lx2, p_min)
  if rkey not in _RET_CACHE:
    _RET_CACHE[rkey] = (
      Retention(lg, lx2, prime_min=2**p_min, explicit_bits=min(p_min, 14), step=0.25, max_bits=600, tilt=0.05, odd=1),
      Retention(600.0, lx2, prime_min=2**p_min, explicit_bits=min(p_min, 14), step=0.25, max_bits=601, tilt=0.005, odd=1),
    )
  ret_obj, cost_obj = _RET_CACHE[rkey]
  # 2. Live 200k-pool self-consistent yield (PRIMARY = FULL DE BRUIJN rho_db + bin-weighted Retention!)
  prim_tail_corr = EXACT_PRIM_T4 / PQ["p_prim"]
  lU = log2_U_z6(lf, lg)
  target_lr = lU + log2(exc) + OP["margin"]
  log2_V_sub0 = PQ["log2_V0"] + log2(PQ["frac_reg"])
  n_in = PQ["n_in"]; bins = PQ["bins"]
  def eval_d(d):
    ys = ys2 = nf_w = ng_w = rg_w = 0.0
    for lnf0, lng0, cnt in bins:
      nf_raw = lnf0 + 48.0 * d; ng_raw = lng0 + 12.0 * d
      nf = nf_raw + ALPHA_F_RAW; ng = ng_raw + ALPHA_G_RAW
      ng_q = round(ng_raw * 2.0) / 2.0
      r_g_bin = (8.0 / 9.0) * ret_obj.at(ng_q, sl, sb) + (1.0 / 9.0) * ret_obj.at(ng_q - 3.0, sl - 3.0, sb)
      y = (rho_db(nf / lf, nf) * rho_db(ng / lg, ng) if use_db else rho(nf / lf) * rho(ng / lg)) * r_g_bin * (1.0 - 0.0014) * prim_tail_corr
      wy = cnt * y; ys += wy; ys2 += cnt * (y * y)
      nf_w += wy * nf_raw; ng_w += wy * ng_raw; rg_w += wy * r_g_bin
    my = ys / n_in
    return log2_V_sub0 + 12.0 * d + log2(max(my, 1e-300)), my, ys2 / n_in, nf_w / ys - CONTENT_BITS, ng_w / ys, rg_w / ys
  lo, hi = -2.5, +2.5
  for _ in range(32):
    mid = 0.5 * (lo + hi)
    if eval_d(mid)[0] < target_lr: lo = mid
    else: hi = mid
  d_opt = 0.5 * (lo + hi)
  lr, my, my2, nf_rel, ng_rel, r_g = eval_d(d_opt)
  se_bits = sqrt(max(my2 - my * my, 0.0) / n_in) / max(my, 1e-300) / LN2
  # 3. Exact C6 geometry + pure-scalar 3-band collection pricing
  log2_V = log2_V_sub0 + 12.0 * d_opt
  log2_Vb = PQ["log2_Vb"] + 6.0 * d_opt
  log2_Na = PQ["log2_Va"] + 6.0 * d_opt
  cells = 2.0 ** log2_V; Vb = 2.0 ** log2_Vb; Na = 2.0 ** log2_Na
  GB = T4["GB"]; g00 = GB[0][0]; det_GB = _det_mat(GB)
  S5 = [[GB[r][c] - GB[r][0] * GB[0][c] / g00 for c in range(1, 6)] for r in range(1, 6)]
  det_S5 = _det_mat(S5)
  Rb = (Vb * sqrt(det_GB) / ((pi**3) / 12.0)) ** (1.0 / 6.0)
  L_max = 2.0 * Rb / sqrt(g00)
  n5 = ((4.0 * pi**2) / 15.0) * (Rb ** 5) / sqrt(det_S5)
  L_b_avg = Vb / n5
  score_arena_log2 = log2((L_max + 2.0) * Na * 2.0 + Na * 2.0)
  vb_score_log2 = log2(Vb * 2.0)
  b2_1p_ram_log2 = log2(16.0 * Na + 9.0 * (2.0 ** pm) / w_sb + (L_max + 2.0) * Na * 2.0 + Na * 2.0)
  pkey = (round(L_max, 2), w_sb)
  if pkey not in _PROF_CACHE:
    _PROF_CACHE[pkey] = _build_sb_profile_pure(L_max, w_sb, PQ["kap"])
  prof = _PROF_CACHE[pkey]
  cr2_chg = ceil(prof["cr2"] * 1.005 * 1000.0) / 1000.0
  # Score-only Stage B cofactorization (ZERO Band-1/2/3 prime recording!):
  c_gb = max(C_GBATCH, gbatch_cost(pbits(0, lf), leaf=9))
  c_b2 = max(C_BATCH2, gbatch_cost(pbits(0, lg), c_leaf=2**19))
  cc = (c_gb, C_TDIV + C_BATCH + c_b2)
  cof_sum = sf_sum = s0_sum = 0.0
  for lnf0, lng0, cnt in bins:
    nf_raw = lnf0 + 48.0 * d_opt - CONTENT_BITS; ng_raw = lng0 + 12.0 * d_opt
    ng_q = round(ng_raw * 2.0) / 2.0
    sf_b = (8.0 / 9.0) * cost_obj.survivor_upper(ng_q, sl, 0.25) + (1.0 / 9.0) * cost_obj.survivor_upper(ng_q - 3.0, sl - 3.0, 0.25)
    s0_b = (8.0 / 9.0) * cost_obj.survivor_upper(ng_q, sl, 0.25, sb, lg, KMAX) + (1.0 / 9.0) * cost_obj.survivor_upper(ng_q - 3.0, sl - 3.0, 0.25, sb, lg, KMAX)
    q_f_b = min(1.0, Q_F_SAFETY * rho((nf_raw + 1.58) / lf))
    cof_sum += cnt * (sf_b * (C_NORMS + C_SB) + s0_b * (cc[0] + q_f_b * cc[1]))
    sf_sum += cnt * sf_b; s0_sum += cnt * s0_b
  cof_b = cof_sum / n_in; sf = sf_sum / n_in; s0 = s0_sum / n_in
  c_hit_lo = 6.25 if strict_per_owner else C_HIT_LO
  c_hit_hi = 14.5 if strict_per_owner else C_HIT_HI
  # Band 1 [2^14, 2^44): bn254_apply_8b_u16_u4 + exact min(bsearch, bn254_b1_early4_find)
  tot_lk = 0.0; dl1 = (p_cut - p_min) / 100.0
  fb1_bins = prof["fb1_bins"]
  for i in range(100):
    l0 = p_min + i * dl1; l1 = l0 + dl1; lm = 0.5 * (l0 + l1)
    dn = LI_PAD * (2.0**log2_li(l1) - 2.0**log2_li(l0))
    mu = Na / (2.0 ** lm)
    c_bs = 24.0 + 8.0 * (log2(mu + 1.0) + 1.0)
    c_e4 = sum(wt * (11.0 + (1.0 - exp(-mu * fa)) * (14.0 + mu * fa)) for fa, wt in fb1_bins)
    tot_lk += dn * min(c_bs, c_e4)
  n1_tot = LI_PAD * (2.0**log2_li(p_cut) - 2.0**log2_li(p_min))
  rep_lo = log(p_cut / p_min) * (c_hit_lo * 1.004 + 0.06) + (tot_lk * Vb + n1_tot * n5 * 24.0) / cells
  # Band 2 [2^44, 2^54): 8-lane SWAR sb_u7 per-SB tight q_min_j gate (W_SB=256, max_sb=5 <= 8 guaranteed!)
  assert prof["max_sb"] <= 8, f"SWAR overflow: max_sb={prof['max_sb']} > 8"
  n_words_row = 1.0
  n2_tot = n2_occ = n2_ent = 0.0; dl2 = (pm - p_cut) / 100.0
  fj_bins = prof["fj_bins"]
  for i in range(100):
    l0 = p_cut + i * dl2; l1 = l0 + dl2; lm = 0.5 * (l0 + l1)
    dn = LI_PAD * (2.0**log2_li(l1) - 2.0**log2_li(l0))
    mu_sb = (w_sb * Na) / (2.0 ** lm)
    p_occ = 1.0 - exp(-mu_sb) if mu_sb < 40.0 else 1.0
    p_ent = sum(wt * (1.0 - exp(-mu_sb * fa)) for fa, wt in fj_bins)
    n2_tot += dn; n2_occ += dn * p_occ; n2_ent += dn * p_ent
  probe = (n2_tot * n5 * n_words_row * C_SWAR_PROBE) / cells
  occ_sb = (n2_ent * n5 * prof["sb_row"] * C_SB_ENTER) / cells
  hi_ap = log(pm / p_cut) * (C_CAND * cr2_chg + c_hit_hi)
  # Band 3 [2^54, 2^59): sigma-orbit-folded 12D Optimal Bounding Ellipsoid E* Fincke-Pohst
  # Active R = {u <= 1, v + kap*u <= 1}, norm vol I(kap) = 1 - 9/4*kap + 9/5*kap^2 - 1/2*kap^3 = 0.445438.
  # Optimal bounding 12D ellipsoid E* = {0.5*u + v/(2*(1-kap)) <= 1} contains R by vertex convexity,
  # with exact provable volume ratio VR_12 = 3.2*(1-kap)^3 / I(kap) = 2.238988x!
  kap_r = PQ["kap"]
  i_kap = 1.0 - 2.25 * kap_r + 1.8 * (kap_r ** 2) - 0.5 * (kap_r ** 3)
  vr_12 = (3.2 * ((1.0 - kap_r) ** 3)) / i_kap
  assert MULT_E12 >= vr_12 / F_FD_UNITS and MULT_FD_ACT >= 1.0 / F_FD_UNITS
  hits_top = log(lx2 / pm)
  dp_top = LI_PAD * 0.5 * (2.0**log2_li(lx2) - 2.0**log2_li(pm))
  V_full12_half = 3.0 * cells * (vr_12 / F_FD_UNITS)
  w_s = npl_s = 0.0
  for ii in range(200):
    lp_b3 = pm + (lx2 - pm) * (ii + 0.5) / 200.0
    wt_b3 = 1.0 / lp_b3
    w_s += wt_b3; npl_s += wt_b3 * fp12_gauss_nph(V_full12_half / (2.0 ** lp_b3))
  fp12_nodes_int = npl_s / w_s
  c_fp_walk = hits_top * (MULT_E12 * (C_FP12_NODE * fp12_nodes_int + C_FP12_LEAF_GATE) + MULT_FD_ACT * C_FP12_LEAF_FD)
  c_babai_init = (dp_top / 3.0) * 500.0 / cells  # 1 root per sigma-orbit (~10^-40 ipc!)
  bkt_top = hits_top * bkt_b3
  per_cell = rep_lo + probe + occ_sb + hi_ap + c_fp_walk + c_babai_init + bkt_top + C_SCAN_ZERO + cof_b
  coll_u = log2(cells * per_cell * OP["work_allow"] / UNIT)
  # Lists build
  n_lo = LI_PAD * 0.5 * (2.0**log2_li(p_cut) - 2.0**log2_li(p_min))
  sp_lo = LI_PAD * 0.5 * (2.0**(2.0 * p_cut)) / (2.0 * p_cut * LN2 - 1.0)
  ins_lo = (n_lo + 1.0) * (C_LIST_A * Na) + C_LBUCKET * LIST_NS * sp_lo
  byt_lo = 16.0 * Na * n_lo + 8.0 * LIST_NS * sp_lo
  n_hi = LI_PAD * 0.5 * (2.0**log2_li(pm) - 2.0**log2_li(p_cut))
  sp_hi = LI_PAD * 0.5 * (2.0**(2.0 * pm)) / (2.0 * pm * LN2 - 1.0)
  ins_hi = n_hi * (C_LIST_A * Na) + C_LBUCKET * (sp_hi / w_sb)
  byt_hi = 16.0 * Na * n_hi + 9.0 * (sp_hi / w_sb)  # 8B offset + 1B sb_u7
  ins_pre = Na * (1000.0 + 12.0 * log2_Na)
  lists_u = log2((ins_lo + ins_hi + ins_pre) / UNIT)
  ram_log2 = log2(byt_lo + byt_hi + 8.0 * cells * hits_top + (L_max + 2.0) * Na * 2.0)
  # 4. 3-Scale W=160 Real-C Ladder (2^20 -> 2^21 -> 2^22) & Fixed-Instance H3 4-Cap LA!
  def stats_of(key):
    rs = FILT[key]; n = len(rs)
    def ms(idx):
      vs = [x[idx] for x in rs]; m = sum(vs) / n
      s2 = sum((v - m) ** 2 for v in vs) / (n - 1) if n > 1 else 0.0
      return m, sqrt(s2 / n), max(vs)
    sw, sw_se, sw_mx = ms(1); nl, nl_se, nl_mx = ms(2); k5, _, k5_mx = ms(3); tr, _, tr_mx = ms(4); uf, _, _ = ms(5); eu, _, _ = ms(6)
    k5_exp_vs = [x[6] * C_UNIT_5 + (1.0 - x[6]) * C_GEN_1 for x in rs]
    k5_exp_mean = sum(k5_exp_vs) / n
    return dict(n=n, N_sw_U=sw, N_sw_se=sw_se, N_sw_max=sw_mx, N_l_U=nl, N_l_se=nl_se, N_l_max=nl_mx,
          k5_ev=k5, k5_exp=k5_exp_mean, k5_exp_max=max(k5_exp_vs), tab_r=tr, tab_r_max=tr_mx, unit_frac=uf, exp_unit=eu)
  cur20 = stats_of("m184_z6_61.0_60.0_e0.82_W160_20")
  cur21 = stats_of("m184_z6_61.0_60.0_e0.82_W160_21")
  cur22 = stats_of("m184_z6_61.0_60.0_e0.82_W160_22")
  assert cur20["n"] == 4 and cur21["n"] == 8 and cur22["n"] == 4
  nu_sw_cap = NU_SW_CAP * cap_scale
  nu_l_cap = NU_L_CAP * cap_scale
  k5_cap = K5_CAP; tab_r_cap = TAB_R_CAP
  pf_rl_cap = min(1.0, nu_l_cap / nu_sw_cap)
  log2_Ns_cap = lU + log2(nu_sw_cap); log2_N_cap = lU + log2(nu_l_cap)
  n_sm = 5
  c_prod = (W_ROW * k5_cap * tab_r_cap + c_tail_z6) / pf_rl_cap
  def bw_steps(lN, m, n): return (1.0 + BW_SIGMA) * (2.0**lN + n_sm) * (1.0 / m + 1.0 / n) + 6.0
  def gen_at_cap(lN_cap, m, n):
    L_cap = bw_steps(lN_cap, m, n); k = m + n
    return log2(4.0 * (k**3) * log2(L_cap) * (fft_F(18.0 * L_cap) + 640.0 * L_cap) / UNIT)
  best_la = None
  for m in (256, 512, 768, 1024):
    for n in (20, 25, 30, 35, 40):
      N = 2.0 ** log2_N_cap; L_st = bw_steps(log2_N_cap, m, n); s = N * (1.0 / m + 1.0 / n)
      kap_bw = L_st / s; hor = max(s, (1.0 + BW_SIGMA) * (N + n_sm) / n + 4.0) / s
      per_n2 = (c_prod + C_PROJ) * (n / m + 1.0) * kap_bw + (1.0 / m + 1.0 / n) * (c_prod + (n_sm * C_FULL + (n - n_sm) * C_Z64) + C_RED9) * hor
      per_n1 = (m + n) * C_RED9
      la_u = log2((N * N * per_n2 + N * per_n1 + PREP_OWNER + PREP_GUARD + PREP_JOINT + LA_RESERVE + 2.0**20) / UNIT)
      gen_u = gen_at_cap(log2_N_cap, m, n)
      both = log2(2.0**la_u + 2.0**gen_u)
      if best_la is None or both < best_la["both"]:
        ncap = 2.0 ** log2_Ns_cap; nslots = ceil(n / 5.0)
        mem_1own = log2(ncap * (160.0 * nslots + 160.0 * nslots + W_ROW * 8.0) + (2.0**55.04) * nslots + 65.0 * (2.0**53))
        best_la = dict(la=la_u, gen=gen_u, both=both, m=m, n=n, mem_1own=mem_1own, log2_N_cap=log2_N_cap, log2_Ns_cap=log2_Ns_cap, c_prod=c_prod)
  pkey_peel = (lf, lg, exc)
  if pkey_peel not in _PEEL_CACHE:
    _PEEL_CACHE[pkey_peel] = _peel_fixed_point(lf, lg, exc, nf=nf_rel, ng=ng_rel)
  peel = _PEEL_CACHE[pkey_peel]
  # 5. LIVE 62-ROW LB=59 PURE QBINS/BOOT_BINS DESCENT RECURRENCE!
  rels_low = lr - 2.25 * se_bits
  dkey = (round(rels_low, 3), round(peel["h"], 4))
  if dkey not in _DESC_CACHE:
    _DESC_CACHE[dkey] = DESC.evaluate_descent(
      D["descent"], rho, fft_F, gbatch_cost, pbits,
      log2_rels=rels_low, peel_h=peel["h"], nf_main=nf_rel, ng_main=ng_rel,
    )
  desc = _DESC_CACHE[dkey]
  # 6. Remaining ledger rows (FULL un-discounted Schirokauer sm_u & 1/3 sigma-orbit 12D LLL bases_u!)
  rels_up = lr + 2.25 * se_bits; u_log = log2(UNIT); lb_max = max(lf, lg)
  canon_u = log2_Na + 42.8 - u_log
  r_bits = R_BN254.bit_length()
  bit_f = 48*48*512 + 7*4*36*512 + 2**17 + 3*48*1024 + 2*48*6*8*40 + 48*1024
  bit_g = 12*12*512 + 1*6*6*512 + 2**15 + 3*12*1024 + 2*12*6*8*40 + 12*1024
  sm_u = rels_up + log2(48 * r_bits * bit_f + 12 * r_bits * bit_g) - u_log
  misc_u = log2(2.0**(rels_up + 30.4) + 2.0**(rels_up + 14.0)) - u_log
  pi_b = log2(1.25506 * (2.0**lb_max) / (lb_max * LN2))
  fb_u = pi_b + 23.4 - u_log
  bases_u = log2((dp_top / 3.0) * C_LLL12 / UNIT)
  split_u = log2(sum(split_parts(2.0**rels_up, lb=lb_max).values())) - u_log + 1.0
  L_gt = ceil(lb_max - 6.0 + 0.53)
  pg = sum((2.0**(lb_max - 6.0 + 0.53 - j - 1)) * (fft_F(2.0**(j + 1)) + 3.0 * (2.0**(j + 1))) for j in range(L_gt))
  gtree_u = log2(2.0 * pg) - u_log
  rows = dict(
    collection=coll_u, chi_setup=65.712288, la=best_la["la"], gen=best_la["gen"],
    lists=lists_u, two_logs=desc["two_log"], canon=canon_u, sm=sm_u, misc=misc_u,
    fb=fb_u, bases=bases_u, split=split_u, gtree=gtree_u,
    cert=50.000001, verify=11.000001, wrapper=53.712288,
  )
  if track == "bn254-pairing-inversion":
    rows["output"] = log2(G2_OUTPUT_BOUND / UNIT)
  exp_run = log2(sum(2.0**v for v in rows.values()))
  # 7. Strict Additive Union Bound (FIX-U) under Fixed-Instance H3 (P(cap abort | H3) = 0)!
  ncap = ceil(2.0 ** best_la["log2_Ns_cap"])
  fail_terms = dict(
    main_matrix=4.0 * ncap / (2.0**64),
    filt_cap_h3_given_H3=0.0,
    split_abort=(2.0**split_u) / (2.0**SPLIT_CAP - (2.0**SPLIT_RES) / UNIT),
    log_phase_cap=2.0 * (2.0**(desc["two_log"] - 1.0 - 112.0)),
    node_abort=2.0 ** desc["abort_log2"],
    repeat_ideal=2.0 ** desc["repeat_log2"],
    group_order_and_sm_comp=7.0 / R_BN254,  # 2/r group order + 5/r random SM compression W_f, W_g
    numerical=2.0 ** (-48.0),
  )
  succ_exact = 1.0 - sum(fail_terms.values())
  succ_decl = f"{floor(succ_exact * 1e6) / 1e6:.6f}"
  s_frac = F(succ_decl)
  K_gate = int(-(-(F(39, 100) * (2**64)) // s_frac))
  gate_passes = (F(K_gate, 2**64) * s_frac >= F(39, 100))
  gate_prev_fails = (F(K_gate - 1, 2**64) * s_frac < F(39, 100))
  exp_run_decl = ceil(exp_run * 1e9) / 1e9
  claim_exact = log2(1.0 + float(F(K_gate, 2**64)) * (2.0**exp_run_decl))
  claim_num = ceil(claim_exact * 1e9)
  claim_decl = f"{claim_num / 1e9:.9f}"
  all4_mono = (
    cur20["N_sw_U"] > cur21["N_sw_U"] > cur22["N_sw_U"] and
    cur20["N_l_U"] > cur21["N_l_U"] > cur22["N_l_U"] and
    cur20["k5_exp"] > cur21["k5_exp"] > cur22["k5_exp"] and
    cur20["tab_r"] > cur21["tab_r"] > cur22["tab_r"]
  )
  all_below_4caps = (
    max(cur20["N_sw_max"], cur21["N_sw_max"], cur22["N_sw_max"]) < nu_sw_cap and
    max(cur20["N_l_max"], cur21["N_l_max"], cur22["N_l_max"]) < nu_l_cap and
    max(cur20["k5_exp_max"], cur21["k5_exp_max"], cur22["k5_exp_max"]) < k5_cap and
    max(cur20["tab_r_max"], cur21["tab_r_max"], cur22["tab_r_max"]) < tab_r_cap
  )
  return dict(
    track=track, use_db=use_db, c_tail_z6=c_tail_z6, cap_scale=cap_scale,
    strict_per_owner=strict_per_owner, bkt_b3=bkt_b3,
    op=OP,
    pool=dict(n_total=PQ["n_total"], n_in=n_in, n_prim_in=PQ["n_prim_in"],
                  frac_reg=PQ["frac_reg"], frac_reg_analytic=prof["frac_reg_analytic"]),
    geometry=dict(log2_U=lU, log2_V=log2_V, se_bits=se_bits, log2_Vb=log2_Vb, log2_Na=log2_Na,
                      log2_N5=log2(n5), L_b_avg=L_b_avg, L_max=L_max, sb_row=prof["sb_row"],
                      max_sb=prof["max_sb"], vr_12=vr_12,
                      cr2_tight=cr2_chg, p_occ_b2=n2_occ/max(n2_tot,1.0), p_ent_b2=n2_ent/max(n2_tot,1.0),
                      nf_rel=nf_rel, ng_rel=ng_rel, r_g=r_g, sf=sf, s0=s0,
                      score_arena_log2=score_arena_log2, vb_score_log2=vb_score_log2,
                      b2_1p_ram_log2=b2_1p_ram_log2, ram_log2=ram_log2),
    per_cell=dict(total=per_cell, rep_lo=rep_lo, probe=probe, occ_sb=occ_sb, hi_ap=hi_ap,
                      fp12_walk=c_fp_walk, fp12_nodes_int=fp12_nodes_int, babai_init=c_babai_init,
                      bkt_top=bkt_top, scan_zero=C_SCAN_ZERO, cof=cof_b),
    la=dict(bw=[best_la["m"], best_la["n"]], log2_N_cap=best_la["log2_N_cap"], log2_Ns_cap=best_la["log2_Ns_cap"],
        W_row=W_ROW, nu_sw_cap=nu_sw_cap, nu_l_cap=nu_l_cap, k5_cap=k5_cap, tab_r_cap=tab_r_cap,
        pf_rl_cap=pf_rl_cap, c_prod=c_prod,
        mem_1own_log2=best_la["mem_1own"], fits_1own=(best_la["mem_1own"] < 64.0),
        h3_mode="fixed_instance_deterministic_4cap",
        ladder_all4_monotonic=all4_mono, all_runs_below_4caps=all_below_4caps,
        filtsim_cur20=cur20, filtsim_cur21=cur21, filtsim_cur22=cur22, peel=peel),
    descent=desc,
    rows=rows,
    exp_run=exp_run, exp_run_decl=f"{exp_run_decl:.9f}",
    fail_terms=fail_terms,
    succ_exact=succ_exact, succ_decl=succ_decl,
    K_gate=K_gate, gate_passes=gate_passes, gate_prev_fails=gate_prev_fails,
    claim_exact=claim_exact, claim_decl=claim_decl, score_numerator=str(claim_num),
  )
