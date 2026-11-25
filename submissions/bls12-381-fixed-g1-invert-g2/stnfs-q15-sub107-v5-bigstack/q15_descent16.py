"""q15_descent16.py — Pure Mode-A k=1 Zero-Retry Boole + Unbounded i.i.d. Fresh-e Wald Boot
16D Q(zeta_15) Descent Evaluator.

Pure Mode-A architecture (zero internal retries, single special-q lattice L_Q, unconditional Boole):
  - Single special-q lattice L_Q (index q in Z^16) per slot; LB = LX, LMAX = 260, 5-bit bins.
  - Pure k = 1 zero-retry node failure bound:
      p_0(b, s) = max(P0_FLOOR, (1.0 + NB_ALPHA * mu(b, s) / NB_R) ** (-NB_R))
 with NB_R = 1.0, NB_ALPHA = 0.5, P0_FLOOR = 0.001, generous MU_MIN >= 500.
  - Bottom-up unconditional Boole tree recursion over c_det deterministic slots:
      c_det(b, s) = 1 + kA*(cdet_max(f, nA) + cdet_max(g, nA)) + kB*(cdet_max(f, n) + cdet_max(g, n))
      phi(b, s)   = p_0(b, s) + kA*(phi_max(f, nA) + phi_max(g, nA)) + kB*(phi_max(f, n) + phi_max(g, n))
      W(b, s)     = own_1(b, s) + kA*(w_max(f, nA) + w_max(g, nA)) + kB*(w_max(f, n) + w_max(g, n))
  - Unbounded i.i.d. fresh-e Wald boot (BOOT16 = (230, 264, 13, 2), ZERO B_max/T_boot_max):
      c_det_one   = kA_b * cdet_max(g, nA_b) + cdet_max(g, n1)
      Phi_tree    = kA_b * phi_max(g, nA_b) + phi_max(g, n1) <= MAX_PHI_TREE = 0.60
      e_boot_runs = 1.0 / (1.0 - Phi_tree),  boot_phi = 0.0,  p_descent_fail = 0.0
      one_log     = e_boot_runs * (boot_own_1run + ch_boot)
  - Explicit BOUND_CELL charged on BOTH sieve and direct-ECM code paths, linear union-bound missing-leaf
 retention max(0.0, 1.0 - miss), and branch-tree-capped L_min = 2, Torsion-16 Edwards (a=+1) 1-rung Theorem T2 ECM.
"""
import functools
import json
import math
from decimal import Decimal, localcontext
from math import exp, factorial, floor, log, log1p, log2
from pathlib import Path

LN2 = log(2.0)
UNIT = 20000.0
KMAX = 64
RSTEP = 1.0 / 1024.0
MONT5 = 367.76


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
   while len(c) > 1 and abs(c[-1]) < big * Decimal("1e-25"):
    c = c[:-1]
   out.append([float(x) for x in c])
  return out


RHO_COEFS = _rho_coefficients()
_RHO_COEFS_REV = tuple(tuple(reversed(c)) if c is not None else None for c in RHO_COEFS)


def rho_series(u):
 if u <= 1.0:
  return 1.0
 k = math.ceil(u)
 if k > KMAX:
  return 0.0
 x, value = k - u, 0.0
 for c in _RHO_COEFS_REV[k]:
  value = value * x + c
 return value if value > 0.0 else 0.0


RHO_TAB = [rho_series(i * RSTEP) for i in range(int(KMAX / RSTEP) + 2)]
_ONE_IDX = int(round(1.0 / RSTEP))


import struct as _ST_ULP
def _step_ulp(v, target, steps=16):
 if v >= 0.0 and math.isfinite(v):
  u = _ST_ULP.unpack("<Q", _ST_ULP.pack("<d", v))[0]
  if target > v:
   return _ST_ULP.unpack("<d", _ST_ULP.pack("<Q", u + steps))[0]
  if target < v:
   return _ST_ULP.unpack("<d", _ST_ULP.pack("<Q", u - steps))[0] if u >= steps else 0.0
 for _ in range(steps):
  v = math.nextafter(v, target)
 return v


RHO_TAB_LO = [_step_ulp(v, 0.0, 16) if i > _ONE_IDX and v > 0.0 else v for i, v in enumerate(RHO_TAB)]
RHO_TAB_HI = [_step_ulp(v, math.inf, 16) if i > _ONE_IDX else v for i, v in enumerate(RHO_TAB)]
import operator as _OP
_SUMPROD = getattr(math, "sumprod", lambda a, b, _m=_OP.mul: sum(map(_m, a, b)))


def rho_lo(u):
 if u <= 1.0:
  return 1.0
 i = math.ceil(u / RSTEP)
 return RHO_TAB_LO[i] if i < len(RHO_TAB_LO) else 0.0


def rho_hi(u):
 if u <= 1.0:
  return 1.0
 i = floor(u / RSTEP)
 if i + 1 >= len(RHO_TAB_HI):
  return 1.0 / factorial(floor(u))
 return RHO_TAB_HI[i]


def kappa(u, t):
 ru = rho_series(u)
 return rho_series(u - t) / ru if ru > 0 else 0.0


@functools.lru_cache(None)
def _miss_side_table(log2_rels, u_main, lx=69.0, main_lb=71.56, orbit=12, p_sieve=69.0, ret_mult=1.0):
 """Precompute strictly-decreasing collection-cut-conditioned lam(b) & right-endpoint weights on [40, lx]."""
 steps = 50
 h = 1.0 / steps
 lo_i = int(round(40.0 * steps))
 hi_i = int(round(lx * steps))
 prev_lam = None
 entries = []
 for i in range(lo_i, hi_i):
  b_lo = i * h
  b_hi = (i + 1) * h
  rm = ret_mult if b_hi > p_sieve else 1.0
  lam = (2.0 ** (log2_rels - b_hi)) * orbit * kappa(u_main, b_hi / main_lb) * rm
  if prev_lam is not None and lam >= prev_lam:
   raise ValueError(f"miss_side lam(b) not strictly decreasing at b={b_hi}: {lam} >= {prev_lam}")
  prev_lam = lam
  lp = 0.75 * lam
  p_thin = exp(-lp) * (1.0 + lp) if lp < 700.0 else 1e-290
  p_p3 = exp(-lam) * (1.0 + lam + 0.5 * lam * lam + (lam * lam * lam) / 6.0) if lam < 700.0 else 1e-290
  p_miss = max(p_p3, p_thin)
  entries.append((b_hi / lx, log(b_hi / b_lo) * p_miss))
 lam40 = (2.0 ** (log2_rels - 40.0)) * orbit * kappa(u_main, 40.0 / main_lb) * (ret_mult if 40.0 > p_sieve else 1.0)
 if lam40 * 0.75 < 700.0:
  raise ValueError(f"miss_side b<40 tail assumption violated: 0.75*lam(40)={0.75*lam40}")
 entries.append((40.0 / lx, log(40.0 / 2.0) * 1e-290))
 return tuple(entries)


@functools.lru_cache(None)
def miss_side(u4, log2_rels, u_main, lx=69.0, main_lb=71.56, orbit=12, p_sieve=69.0, ret_mult=1.0):
 """Rigorous right-endpoint upper Riemann sum (0.02-bit step) on [40, lx], plus explicit b < 40 analytic tail."""
 u = u4 / 4.0
 ru = rho_series(u)
 if ru <= 0.0:
  return 0.0
 entries = _miss_side_table(log2_rels, u_main, lx, main_lb, orbit, p_sieve, ret_mult)
 tot = 0.0
 for t, w in entries:
  tot += rho_series(u - t) * w
 return tot / ru


def missing(Pm, nf, ng):
 f = lambda m, um, lb, ps, rm: miss_side(
 math.ceil(4.0 * m / Pm.LX), Pm.LOG2_RELS, um, Pm.LX, lb, Pm.ORBIT, ps, rm
 )
 return (
 (f(nf, Pm.MAIN_UF, Pm.MAIN_LB_F, Pm.P_SIEVE_F, Pm.RET_MULT_F) if nf else 0.0)
 + (f(ng, Pm.MAIN_UG, Pm.MAIN_LB, Pm.P_SIEVE_G, Pm.RET_MULT_G) if ng else 0.0)
 )


KA_MAX = 20


def _conv(a, b):
 out = [0.0] * (len(a) + len(b) - 1)
 nz_a = [i for i, x in enumerate(a) if x]
 nz_b = [j for j, y in enumerate(b) if y]
 if not nz_a or not nz_b:
  return out
 a0, a1 = nz_a[0], nz_a[-1] + 1
 b0, b1 = nz_b[0], nz_b[-1] + 1
 aa = a[a0:a1]
 bb_rev = b[b0:b1][::-1]
 La, Lb = a1 - a0, b1 - b0
 pad = [0.0] * (Lb - 1)
 ap = pad + aa + pad
 sp = _SUMPROD
 base = a0 + b0
 out[base : base + La + Lb - 1] = [sp(ap[k : k + Lb], bb_rev) for k in range(La + Lb - 1)]
 return out


@functools.lru_cache(None)
def window_rows(n, base, iA, mesh, kmax=KA_MAX):
 h = (n - base) / mesh
 wA = [log1p(h / (base + i * h)) if i < iA else 0.0 for i in range(mesh)]
 wB = [log1p(h / (base + i * h)) if i >= iA else 0.0 for i in range(mesh)]
 PA = [[1.0]]
 for j in range(1, kmax + 1):
  PA.append([x / j for x in _conv(PA[-1], wA)])
 PB = [[1.0], wB]
 rows = {}
 for jA in range(kmax + 1):
  for jB in range(2):
   r = PA[jA] if jB == 0 else _conv(PA[jA], PB[1])
   i0 = next((idx for idx, v in enumerate(r) if v), len(r))
   i1 = len(r) - next((idx for idx, v in enumerate(reversed(r)) if v), len(r))
   rows[jA, jB] = (r, i0, i1)
 return rows, h


@functools.lru_cache(None)
def ss_cum2(m, base, n, iA, mesh, alpha, trunc, kmax=KA_MAX, need_b=1):
 rows, h = window_rows(n, base, iA, mesh, kmax)
 out = [[0.0, 0.0] for _ in range(kmax + 1)]
 tab, ntab, ceil = RHO_TAB_LO, len(RHO_TAB_LO), math.ceil
 sp = _SUMPROD
 inv_br = 1.0 / (base * RSTEP)
 hk = h * inv_br
 tr_span = int(trunc * base / h) + 1 if (trunc and h > 0) else 0
 jBs = (0, 1) if need_b else (0,)
 j_lim = kmax + (2 if need_b else 1)
 # base * inv_br == 1 / RSTEP == 1024 exactly; precompute j=0 ceiling indices once across all j.
 step_j = int(round(1.0 / RSTEP))
 k0_0 = (m + alpha) * inv_br
 plan = []
 g_lo, g_hi = 10**9, -1
 for j in range(j_lim):
  items = []
  ifeas_j = int(floor((m - j * h - j * base) / h)) if h > 0 else 10**9
  k0 = k0_0 - j * step_j
  i_ntab = int((k0 - (ntab - 1)) / hk) if hk > 0.0 else 0
  for jB in jBs:
   jA = j - jB
   if not (0 <= jA <= kmax):
    continue
   row, r_lo, r_hi = rows[jA, jB]
   if r_lo >= r_hi:
    continue
   ifeas = r_hi - 1 if (r_hi - 1) < ifeas_j else ifeas_j
   if ifeas < r_lo:
    continue
   imin = r_lo if not trunc else (r_lo if r_lo > ifeas - tr_span else ifeas - tr_span)
   if i_ntab > imin:
    imin = i_ntab
   if imin <= ifeas:
    items.append((jA, jB, row, imin, ifeas))
  if not items:
   continue
  t0 = min(it[3] for it in items)
  t1 = max(it[4] for it in items)
  if t0 < g_lo:
   g_lo = t0
  if t1 > g_hi:
   g_hi = t1
  plan.append((j, k0, t0, t1, items))
 if plan:
  idx0 = [ceil(k0_0 - i * hk) for i in range(g_lo, g_hi + 1)]
  for j, k0, t0, t1, items in plan:
   sh = j * step_j
   sl = idx0[t0 - g_lo : t1 - g_lo + 1]
   if k0 - t1 * hk <= 0.0 or k0 - t0 * hk >= ntab:
    arr = [1.0 if (k := x - sh) <= 0 else (tab[k] if k < ntab else 0.0) for x in sl]
   else:
    arr = [tab[x - sh] for x in sl]
   for jA, jB, row, imin, ifeas in items:
    out[jA][jB] = sp(row[imin : ifeas + 1], arr[imin - t0 : ifeas - t0 + 1])
 cum = [[0.0, 0.0] for _ in range(kmax + 1)]
 for jB in jBs:
  acc = 0.0
  for jA in range(kmax + 1):
   acc += out[jA][jB]
   cum[jA][jB] = acc
 return tuple(tuple(x) for x in cum)


def ss_lower2(m, base, n, iA, kA, kB, mesh, alpha, trunc=0.0):
 m = math.ceil(m)
 if n <= base:
  return rho_lo((m + alpha) / base)
 nb = 1 if kB else 0
 cum = ss_cum2(m, base, n, iA, mesh, alpha, trunc, kA, nb)
 return cum[kA][0] + (cum[kA][1] if nb else 0.0)


@functools.lru_cache(None)
def _ss_lower2_tau_c(m, base, n, iA, kA, kB, mesh, alpha, tau):
 rows, h = window_rows(n, base, iA, mesh, kA)
 tab, ntab, ceil = RHO_TAB_LO, len(RHO_TAB_LO), math.ceil
 inv_br = 1.0 / (base * RSTEP)
 hk = h * inv_br
 sp = _SUMPROD
 step_j = int(round(1.0 / RSTEP))
 k0_0 = (m + alpha) * inv_br
 cap = m if m < tau else tau
 jBs = (0, 1) if kB else (0,)
 j_lim = kA + (2 if kB else 1)
 plan = []
 g_lo, g_hi = 10**9, -1
 for j in range(j_lim):
  items = []
  ifeas_j = int(floor((cap - j * h - j * base) / h)) if h > 0 else 10**9
  k0 = k0_0 - j * step_j
  i_ntab = int((k0 - (ntab - 1)) / hk) if hk > 0.0 else 0
  for jB in jBs:
   jA = j - jB
   if not (0 <= jA <= kA):
    continue
   row, r_lo, r_hi = rows[jA, jB]
   if r_lo >= r_hi:
    continue
   ifeas = r_hi - 1 if (r_hi - 1) < ifeas_j else ifeas_j
   if ifeas < r_lo:
    continue
   imin = r_lo if r_lo > i_ntab else i_ntab
   if imin <= ifeas:
    items.append((row, imin, ifeas))
  if not items:
   continue
  t0 = min(it[1] for it in items)
  t1 = max(it[2] for it in items)
  if t0 < g_lo:
   g_lo = t0
  if t1 > g_hi:
   g_hi = t1
  plan.append((j, k0, t0, t1, items))
 tot = 0.0
 if plan:
  idx0 = [ceil(k0_0 - i * hk) for i in range(g_lo, g_hi + 1)]
  for j, k0, t0, t1, items in plan:
   sh = j * step_j
   sl = idx0[t0 - g_lo : t1 - g_lo + 1]
   if k0 - t1 * hk <= 0.0 or k0 - t0 * hk >= ntab:
    arr = [1.0 if (k := x - sh) <= 0 else (tab[k] if k < ntab else 0.0) for x in sl]
   else:
    arr = [tab[x - sh] for x in sl]
   for row, imin, ifeas in items:
    tot += sp(row[imin : ifeas + 1], arr[imin - t0 : ifeas - t0 + 1])
 return tot


def ss_lower2_tau(m, base, n, iA, kA, kB, mesh, alpha, tau):
 """Exact ss_lower2 restricted to window-prime sum j*base + (i+j)*h <= tau (2-stage boot cofactor cap)."""
 m = math.ceil(m)
 if n <= base:
  return rho_lo((m + alpha) / base)
 return _ss_lower2_tau_c(m, base, n, iA, kA, kB, mesh, alpha, tau)


def surv_upper(m, base, T):
 """Clamped uniform 2.0x H-EA upper bound min(1.0, 2.0 * psurv512_hi(m, base, T)) for all m >= 0, T >= 0.
 Uses the identical 512-mesh IEEE-754 outward-rounded interval psurv_iv(m, base, T) across ALL T >= 0
 including T == 0 (pure base-smooth check R == 1), guaranteeing surv_upper >= min(1.0, 2.0 * P_BD(m/base, T/base))."""
 if T < 0.0:
  return 0.0
 if T >= m:
  return 1.0
 return min(1.0, 2.0 * psurv(m, base, T))


def smooth_part_upper(c_bits, base, mesh=512):
 """Computer-checked IEEE-754 outward-rounded uniform-over-all-norms upper bound on the
 constant-threshold sieve-survival event S[cell] >= c_bits (with 1/16-bit rounded-down fixed-point
 weights w_p = floor(16*log2(p))/16 <= log2(p) on sieved primes p >= 2^8 = 256, giving both
 S[cell] <= log2(s_smooth) deterministically and rounding deficit < (m/8)*(1/16) = m/128 matching
 q15_retain16.py::loss's slack - 1.0 - m/128.0 reservation):
      - For 0 < m < c_bits: smooth part <= N < 2^c_bits, so sieve survival is 0 deterministically.
      - For m >= c_bits (u = m/base >= c = c_bits/base): since 0 <= omega(t) <= 1 and rho >= 0
  is non-increasing, P_BD(u, u - c) <= rho(c) + int_c^inf rho(w) dw <= r_hi[k_c] + tail_hi[k_c]
  with k_c = floor(nextafter(c_bits/base, 0) * mesh) <= exact(c_bits/base) * mesh,
  yielding min(1.0, 2.0 * (r_hi[k_c] + tail_hi[k_c])) >= min(1.0, 2.0 * P_BD(u, u - c)).
 Used to price NORMS across ALL sieve survivors (both sub-floor A^c and above-floor A) before floor admission."""
 if c_bits <= 0.0:
  return 1.0
 c_lo = math.nextafter(c_bits / base, 0.0)
 um = max(45.0, math.ceil(c_lo) + 2.0)
 m, _, _, r_hi, _, _, _, tail_hi = psurv_iv_tab(um, mesh)
 k_c = max(0, int(floor(c_lo * m)))
 ub = math.nextafter(r_hi[k_c] + tail_hi[k_c], math.inf)
 return min(1.0, math.nextafter(2.0 * ub, math.inf))


@functools.lru_cache(None)
def psurv_iv_tab(u_max=45.0, m=512):
 """IEEE-754 outward-rounded interval Buchstab-Dickman tables (rho_lo <= rho <= rho_hi,
 om_lo <= omega <= om_hi, plus upper tail int_w^inf rho(s)ds <= tail_hi[floor(w*m)]) on grid 1/m."""
 na = math.nextafter
 N = int(math.ceil(u_max * m)) + 2
 h_hi = na(1.0 / m, math.inf)
 h_lo = na(1.0 / m, 0.0)
 r_hi = [1.0] * (m + 1)
 r_lo = [1.0] * (m + 1)
 B = 32
 assert m % B == 0
 blk_hi = []; blk_lo = []
 suf_hi = [0.0] * (N + 1); suf_lo = [0.0] * (N + 1)
 def _close_blk(bi):
  s0 = bi * B; s1 = s0 + B; ah = 0.0; al = 0.0
  for idx in range(s1 - 1, s0 - 1, -1):
   ah = na(ah + r_hi[idx], math.inf); al = na(al + r_lo[idx], 0.0)
   suf_hi[idx] = ah; suf_lo[idx] = al
  blk_hi.append(ah); blk_lo.append(al)
 for bi in range(m // B):
  _close_blk(bi)
 pref_hi = r_hi[m]; pref_lo = r_lo[m]
 for k in range(m + 1, N + 1):
  u = k / m
  if k <= 2 * m:
   v = 1.0 - log(u)
   vh = na(na(v, math.inf), math.inf)
   vl = max(0.0, na(na(v, 0.0), 0.0))
  else:
   b1 = k // B; b0 = (k - m) // B
   if k == b1 * B:
    sh = math.fsum(blk_hi[b0:b1])
    sl = math.fsum(blk_lo[b0 + 1 : b1]) + suf_lo[k - m + 1]
   else:
    sh = pref_hi + suf_hi[k - m] + math.fsum(blk_hi[b0 + 1 : b1])
    s1_idx = k - m + 1
    sl = (pref_lo + math.fsum(blk_lo[b0 + 1 : b1])) if s1_idx == (b0 + 1) * B else (pref_lo + suf_lo[s1_idx] + math.fsum(blk_lo[b0 + 1 : b1]))
   vh = sh / k; vl = sl / (k - 1)
   for _ in range(16):
    vh = na(vh, math.inf); vl = na(vl, 0.0)
   vl = max(0.0, vl)
  r_hi.append(vh); r_lo.append(vl)
  if (k + 1) % B == 0:
   _close_blk((k + 1) // B - 1)
   pref_hi = 0.0; pref_lo = 0.0
  else:
   pref_hi = na(pref_hi + vh, math.inf); pref_lo = na(pref_lo + vl, 0.0)
 om_hi = [0.0] * N
 om_lo = [0.0] * N
 I_hi = [1.0] * (N + 1)
 I_lo = [1.0] * (N + 1)
 for j in range(m, N):
  if j < 2 * m:
   vh = m / float(j)
   vl = m / float(j + 1)
  else:
   vh = I_hi[j + 1 - m] * m / float(j)
   vl = I_lo[j - m] * m / float(j + 1)
  for _ in range(4):
   vh = na(vh, math.inf)
   vl = na(vl, 0.0)
  om_hi[j] = vh
  om_lo[j] = max(0.0, vl)
  ih = I_hi[j] + h_hi * vh
  il = I_lo[j] + h_lo * vl
  for _ in range(4):
   ih = na(ih, math.inf)
   il = na(il, 0.0)
  I_hi[j + 1] = ih
  I_lo[j + 1] = il
 tail_hi = [0.0] * (N + 1)
 iu_max = int(u_max)
 acc_t = math.nextafter(2.0 * math.exp(-math.lgamma(iu_max + 1)), math.inf) if iu_max > 170 else (2.0 / factorial(iu_max))
 tail_hi[N] = acc_t
 for j in range(N - 1, -1, -1):
  acc_t = na(acc_t + h_hi * r_hi[j], math.inf)
  tail_hi[j] = acc_t
 return m, h_hi, h_lo, tuple(r_hi), tuple(r_lo), tuple(om_hi), tuple(om_lo), tuple(tail_hi)


@functools.lru_cache(None)
def _psurv_iv_rev(u_max=45.0, m=512):
 _, _, _, r_hi, r_lo, _, _, _ = psurv_iv_tab(u_max, m)
 return r_hi[::-1], r_lo[::-1]


@functools.lru_cache(None)
def psurv_iv(m_bits, B, T, mesh=512):
 """Computer-checked outward-rounded interval enclosure [ps_lo, ps_hi] containing exact P_BD(m_bits/B, T/B)."""
 if T < 0.0:
  return 0.0, 0.0
 if m_bits <= T:
  return 1.0, 1.0
 u = m_bits / B
 v = T / B
 um = max(45.0, math.ceil(u) + 2.0)
 m, h_hi, h_lo, r_hi, r_lo, om_hi, om_lo, _ = psurv_iv_tab(um, mesh)
 r_hi_rev, r_lo_rev = _psurv_iv_rev(um, mesh)
 u_lo_i = int(floor(u * m))
 u_hi_i = int(math.ceil(u * m))
 hi = min(u, v)
 nr = len(r_hi)
 ans_hi = r_hi[u_lo_i]
 ans_lo = r_lo[u_hi_i] if u_hi_i < nr else 0.0
 if hi <= 1.0:
  return min(1.0, ans_lo), min(1.0, ans_hi)
 j_end_hi = int(math.ceil(hi * m))
 j_end_lo = int(floor(hi * m))
 sp = _SUMPROD
 j_cut_hi = u_lo_i - 1 if (u_lo_i - 1) < j_end_hi else j_end_hi
 s_hi = 0.0
 if j_cut_hi > m:
  s_hi += sp(r_hi_rev[nr - u_lo_i + m : nr - u_lo_i + j_cut_hi], om_hi[m : j_cut_hi])
 if j_end_hi > (j_cut_hi if j_cut_hi > m else m):
  s_hi += math.fsum(om_hi[max(m, j_cut_hi) : j_end_hi])
 j_start_lo = m if (u_hi_i - m) < nr else (u_hi_i - nr + 1)
 j_cut_lo = u_hi_i if u_hi_i < j_end_lo else j_end_lo
 s_lo = 0.0
 if j_cut_lo > j_start_lo:
  s_lo += sp(r_lo_rev[nr - 1 - u_hi_i + j_start_lo : nr - 1 - u_hi_i + j_cut_lo], om_lo[j_start_lo : j_cut_lo])
 if j_end_lo > (j_cut_lo if j_cut_lo > j_start_lo else j_start_lo):
  s_lo += math.fsum(om_lo[max(j_start_lo, j_cut_lo) : j_end_lo])
 na = math.nextafter
 for _ in range(8):
  s_hi = na(s_hi, math.inf); s_lo = na(s_lo, 0.0)
 v_hi = min(1.0, na(ans_hi + na(h_hi * s_hi, math.inf), math.inf))
 v_lo = min(1.0, max(0.0, na(ans_lo + na(h_lo * max(0.0, s_lo), 0.0), 0.0)))
 return v_lo, v_hi


def psurv(m_bits, B, T, mesh=512):
 """Strict IEEE-754 outward-rounded upper bound on P_BD(m_bits/B, T/B)."""
 return psurv_iv(m_bits, B, T, mesh)[1]


def ecm_instr(fbits, nbits, success):
 import q15_ecm16
 return q15_ecm16.ecm_instr(fbits, nbits, success)


def log2li(l):
 x = l * LN2
 return l - log2(x) + log2(1.0 + 1.0 / x + 2.0 / (x * x))


def log2li_upper(l):
 """Explicit upper bracket on log2(pi_K(2^l)) for l >= 40 (Dusart 2.334/x^2 + 1% finite-range margin)."""
 x = l * LN2
 return l - log2(x) + log2(1.01 * (1.0 + 1.0 / x + 2.334 / (x * x)))


CELLS = {}


def load_cells():
 global CELLS
 if CELLS_16D:
  for k, v in CELLS_16D.items():
   if k in CELLS:
    assert CELLS[k] == v, f"CELLS mismatch at R={k}"
   else:
    CELLS[k] = v
 elif not CELLS:
  cells, _, _ = compute_cells_16d(RADII_16D)
  CELLS.update(cells)
 return CELLS


class P:
 # Mode-A pure k=1 zero-retry Boole + unbounded i.i.d. fresh-e Wald boot parameters
 LB = 69.0
 LX = 69.0
 NB_R = 1.0
 NB_ALPHA = 0.5
 P0_FLOOR = 0.001
 MU_MIN = 800.0
 MAX_PHI_TREE = 0.60
 PMIN = 8.0
 SLACK = 96.0
 ALPHA_F = 2.74
 ALPHA_G = 0.77
 LOG2U8 = 8.0 * log2(0xd201000000010000)
 LOG2P4 = 4.0 * 380.99700
 QBINS = [
 (1.0 / 6.0, 33.37, 46.5),
 (1.0 / 6.0, 47.14, 46.5),
 (1.0 / 6.0, 57.50, 46.5),
 (1.0 / 6.0, 67.20, 46.5),
 (1.0 / 6.0, 79.76, 46.5),
 (1.0 / 6.0, 131.91, 46.5),
 ]
 QBINS_LO = [
 (-64.83, -67.35),
 (-4.63, -67.35),
 (9.14, -67.35),
 (19.50, -67.35),
 (29.20, -67.35),
 (41.76, -67.35),
 ]
 QBINS_SUB = []
 BOOT_BINS = [(0.5, -260.0), (0.4, -250.0), (0.1, -225.22)]
 BOOT_SUB = []
 BOOT_KEEP_FRAC = 0.90
 TAU_BOOT = 2095
 BOOT_LO_OFF = -300.0
 BOOT_SURV_X = 2.0
 PRCAP_NODE = 4.0
 PRCAP_BOOT = 1.5
 LOG2_RELS = 63.235
 MAIN_LB_F = 71.20
 MAIN_LB = 71.56
 MAIN_UF = 545.0 / MAIN_LB_F
 MAIN_UG = 665.0 / MAIN_LB
 P_SIEVE_F = 60.4
 P_SIEVE_G = 62.4
 RET_MULT_F = 0.916
 RET_MULT_G = 0.82
 ORBIT = 12
 UDUP = 0.40
 UDUP_BOOT = 0.40
 UDUP_NODE = 0.55
 E16_COUNT_LO = 0.95
 E16_COUNT_HI = 1.05
 BALL_PAIR = 0.0
 BOUND_CELL = 4400.0
 SIEVE_HIT = 1400.0
 SCAN = 10.0
 ZERO = 2.0
 PAIR_SETUP = 2.0 ** 35.2
 LLL_NODE = 2.0 ** 40.0
 LLL_BOOT = 2.0 ** 47.0
 NORMS = 2.0 ** 27.0
 TDIV = 2.0 ** 25.0
 IDEAL = 2.0 ** 35.0
 PRIME = 2.0 ** 24.0
 ECM_SUCCESS = 0.999
 MESH = 384
 TRUNC = 0.0
 LMAX = 260
 BIN = 5
 FALLBACK_RET = None


def p0_modeA(mu, Pm=P):
 """Pure Mode-A k = 1 zero-retry slot failure probability: max(P0_FLOOR, (1 + NB_ALPHA * mu / NB_R)^(-NB_R))."""
 return max(Pm.P0_FLOOR, (1.0 + Pm.NB_ALPHA * mu / Pm.NB_R) ** (-Pm.NB_R))


def load_stats(stats_or_path, Pm=None):
 """Load QBINS, QBINS_LO, QBINS_SUB, BOOT_BINS, BOOT_SUB, CELLS, UDUP, tower, and retain16 tables from q15_stats."""
 st = json.loads(Path(stats_or_path).read_text()) if isinstance(stats_or_path, (str, Path)) else stats_or_path
 tgt = Pm if Pm is not None else P
 if "cells_16d" in st:
  CELLS.update({int(k): int(v) for k, v in st["cells_16d"].items()})
 load_cells()
 tw = st.get("tower", {})
 if "log2_U8" in tw:
  tgt.LOG2U8 = float(tw["log2_U8"])
 if "log2_P4" in tw:
  tgt.LOG2P4 = float(tw["log2_P4"])
 if "alpha_f" in tw:
  tgt.ALPHA_F = float(tw["alpha_f"])
 if "alpha_g" in tw:
  tgt.ALPHA_G = float(tw["alpha_g"])
 bnp = st.get("ball_norm_pool", st)
 for key in ("QBINS", "QBINS_LO", "QBINS_SUB"):
  if key in bnp:
   setattr(tgt, key, [tuple(x) for x in bnp[key]])
 btp = st.get("boot_norm_pool", st)
 if "BOOT_BINS" in btp:
  tgt.BOOT_BINS = [tuple(x) for x in btp["BOOT_BINS"]]
 if "BOOT_SUB" in btp:
  tgt.BOOT_SUB = list(btp["BOOT_SUB"])
 if "BOOT_KEEP_FRAC" in btp:
  tgt.BOOT_KEEP_FRAC = float(btp["BOOT_KEEP_FRAC"])
 ud = st.get("unit_orbit_udup", {})
 if "recommended_UDUP" in ud:
  tgt.UDUP = float(ud["recommended_UDUP"])
 if "udup_boot_charged" in ud:
  tgt.UDUP_BOOT = float(ud["udup_boot_charged"])
 if "udup_node_charged" in ud:
  tgt.UDUP_NODE = float(ud["udup_node_charged"])
 if "e16_count_lo" in ud:
  tgt.E16_COUNT_LO = float(ud["e16_count_lo"])
 if "e16_count_hi" in ud:
  tgt.E16_COUNT_HI = float(ud["e16_count_hi"])
 if "prime_structure" in st and "omitted_law" in st:
  try:
   import q15_retain16
   q15_retain16.load_table(st)
  except ImportError:
   pass
 return st


def pairs_count(Pm):
 return 2.0 * (2.0 ** log2li_upper(Pm.LX))


def _get_retention(Pm, nf_max, ng_max, n, iA, kA, kB):
 if Pm.FALLBACK_RET is not None:
  return Pm.FALLBACK_RET, 0.0, 0.0
 import q15_retain16
 return q15_retain16.retention(
 nf_max, ng_max, n, iA, kA, kB,
 slack=Pm.SLACK, lx=Pm.LX, mesh=Pm.MESH,
 alpha_f=Pm.ALPHA_F, alpha_g=Pm.ALPHA_G,
 )


def _get_counted_fraction(Pm):
 if Pm.FALLBACK_RET is not None:
  return Pm.FALLBACK_RET
 import q15_retain16
 return q15_retain16.counted_fraction()


def bins_of(Pm):
 """Bins from LX up to LMAX: first bin [LX, 75), then 5-bit bins [75, 80), ..., [255, 260)."""
 first_hi = 75.0 if Pm.LX < 75.0 else float(int(Pm.LX // Pm.BIN) * Pm.BIN + Pm.BIN)
 bins = [(float(Pm.LX), first_hi)]
 l = first_hi
 while l < Pm.LMAX - 1e-9:
  bins.append((float(l), float(min(Pm.LMAX, l + Pm.BIN))))
  l += Pm.BIN
 return bins


def node_cost(Pm, l_lo, l_hi, side, n, iA, kA, kB, R):
 """Evaluate ONE pure k=1 zero-retry single-special-q (L_Q) descent slot with linear missing-leaf retention
 max(0.0, 1.0 - miss), 100% eager cofactorization, BOUND_CELL on BOTH sieve and direct-ECM branches,
 and branch-tree-capped L_min=2, Torsion-16 Edwards (a=+1) 1-rung Theorem T2 ECM."""
 cells_map = load_cells()
 pairs = pairs_count(Pm)
 cfmax = max(b[1] for b in Pm.QBINS)
 cgmax = max(b[2] for b in Pm.QBINS)
 C_R = cells_map[R]
 S = log2(C_R)
 c_hi = (l_hi + S) / 16.0
 c_lo = (l_lo + S) / 16.0
 s_f_hi = l_lo if side == "f" else 0.0
 s_g_hi = l_lo if side == "g" else 0.0
 s_f_lo = l_hi if side == "f" else 0.0
 s_g_lo = l_hi if side == "g" else 0.0

 h = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
 nA = Pm.LX + iA * h if n > Pm.LX else Pm.LX
 T = kA * nA + kB * n + Pm.SLACK
 T0 = kA * nA + kB * n  # Exact post-ecm(LX) zero-slack cofactor bound (ZERO sieve slack!)

 nf_floor = 72.0 * c_lo + min(x[0] for x in Pm.QBINS_LO) - s_f_lo
 ng_floor = 24.0 * c_lo + Pm.LOG2U8 + min(x[1] for x in Pm.QBINS_LO) - s_g_lo
 ps_f0 = surv_upper(nf_floor, Pm.LX, T0)
 ps_g0 = surv_upper(ng_floor, Pm.LX, T0)
 sv0 = ps_f0 * ps_g0
 # Global exact norm-floor Stage-1 bounds (ZERO sextile-frequency averaging in node_cost!):
 ps_fT = surv_upper(nf_floor, Pm.LX, T)
 ps_gT = surv_upper(ng_floor, Pm.LX, T)
 sv = ps_fT * ps_gT
 sv_fT_g0 = ps_fT * ps_g0
 sv_f0_gT = ps_f0 * ps_gT
 # Pre-NORMS uniform-over-all-norms (0 < m <= nf_max, including ALL sub-floor cells m < nf_floor!)
 # sieve-survivor bound under constant node-wide 1/16-bit fixed-point sieve threshold S >= m_floor - T:
 sv_norms = smooth_part_upper(nf_floor - T, Pm.LX) * smooth_part_upper(ng_floor - T, Pm.LX)
 pr = 0.0
 for w, cf, cg in Pm.QBINS:
  nf = 72.0 * c_hi + cf - s_f_hi
  ng = 24.0 * c_hi + Pm.LOG2U8 + cg - s_g_hi
  pf = ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
  pg = ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
  pr += w * pf * pg
 if Pm.QBINS_SUB:
  pri = 0.0
  for cf, cg in Pm.QBINS_SUB:
   nf = 72.0 * c_hi + cf - s_f_hi
   ng = 24.0 * c_hi + Pm.LOG2U8 + cg - s_g_hi
   pri += (
   ss_lower2(nf, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_F, Pm.TRUNC)
   * ss_lower2(ng, Pm.LX, n, iA, kA, kB, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
   )
  pr = min(pri / len(Pm.QBINS_SUB), Pm.PRCAP_NODE * pr)
 if pr <= 0.0:
  return None

 nf_max = 72.0 * c_hi + cfmax - s_f_hi
 ng_max = 24.0 * c_hi + Pm.LOG2U8 + cgmax - s_g_hi
 nf_fac = math.ceil(nf_max / 16.0 + kA + kB)
 ng_fac = math.ceil(ng_max / 16.0 + kA + kB)
 det = 1.0 - (1.0 - Pm.ECM_SUCCESS) * ((nf_max + ng_max) / 16.0 + 2.0 * (kA + kB))
 miss = missing(Pm, nf_max, ng_max)
 lin_keep = max(0.0, 1.0 - miss)
 if lin_keep <= 0.0:
  return None
 C_R_lo = C_R * Pm.E16_COUNT_LO
 C_R_hi = C_R * Pm.E16_COUNT_HI
 mu = C_R_lo * pr * Pm.UDUP_NODE * det * lin_keep

 n_lp_max = 2 * (kA + kB)
 nt = Pm.NORMS + Pm.TDIV  # Exact integer norms |N_f|, |N_g| enforcing lower floors (nf_floor, ng_floor) before TDIV/ECM
 ef1 = ecm_instr(Pm.LX, nf_max, Pm.ECM_SUCCESS)
 eg1 = ecm_instr(Pm.LX, ng_max, Pm.ECM_SUCCESS)
 e_s1 = nt + ef1 + eg1
 if kA + kB == 0:
  e_s2 = 0.0
 else:
  e_s2 = (
  2.0 * ecm_instr(nA, T0, Pm.ECM_SUCCESS)
  + max(2, n_lp_max) * Pm.PRIME
  + max(1, n_lp_max) * Pm.IDEAL
  )
 # Pre-admission NORMS is priced separately at sv_norms (uniform over ALL m > 0, covering BOTH sub-floor A^c
 # and above-floor A sieve survivors). Inside NORMS, exact |N_f|, |N_g| deterministically reject any cell
 # with log2(|N|/q) < m_floor or S < log2(|N|/q) - T BEFORE TDIV and Stage-1 ECM run.
 # Unconditional sequential Stage-1 ECM (g-first vs f-first) runs strictly AFTER admission (m >= m_floor):
 # {sieve_f(T) & sieve_g(T) & R_g <= 2^T0} subset {sieve_f(T) & R_g <= 2^T0} (rate <= sv_fT_g0 + pr),
 # and symmetrically {sieve_f(T) & sieve_g(T) & R_f <= 2^T0} subset {R_f <= 2^T0 & sieve_g(T)} (rate <= sv_f0_gT + pr).
 # Every per-cell work term is conservatively scaled by E16_COUNT_HI = 1.05 (covering E_16 cap C_30 lattice-count variance).
 s2_tot = C_R_hi * (sv0 + pr) * e_s2
 c_norms_bc = C_R_hi * (sv_norms + pr) * Pm.NORMS
 c_s1_bc = c_norms_bc + C_R_hi * min(
 (sv + pr) * (Pm.TDIV + eg1) + (sv_fT_g0 + pr) * ef1,
 (sv + pr) * (Pm.TDIV + ef1) + (sv_f0_gT + pr) * eg1,
 )
 hits_pc = 2.0 * log(Pm.LX / Pm.PMIN)
 sieve = C_R_hi * (hits_pc * Pm.SIEVE_HIT + Pm.SCAN + Pm.ZERO + Pm.BOUND_CELL) + pairs * Pm.PAIR_SETUP
 own_sieve_bc = Pm.LLL_NODE + sieve + c_s1_bc + s2_tot
 trials_direct = C_R_hi * (Pm.BOUND_CELL + e_s1 + (sv0 + pr) * e_s2)
 own_direct = Pm.LLL_NODE + trials_direct
 if own_direct < own_sieve_bc:
  own_1 = own_direct
  how = "ecm"
  survivors = C_R_hi
  surv_cost_tot = C_R_hi * (e_s1 + (sv0 + pr) * e_s2)
  ret, loss_f, loss_g = _get_counted_fraction(Pm), 0.0, 0.0
  mu *= ret
 else:
  own_1 = own_sieve_bc
  how = "sieve"
  survivors = C_R_hi * (sv + pr)
  surv_cost_tot = c_s1_bc + s2_tot
  ret, loss_f, loss_g = _get_retention(Pm, nf_max, ng_max, n, iA, kA, kB)
  mu *= ret
 if mu < Pm.MU_MIN:
  return None
 p0 = p0_modeA(mu, Pm)
 return dict(
 own_1=own_1,
 own_1_log2=log2(own_1),
 n=n,
 iA=iA,
 nA=nA,
 kA=kA,
 kB=kB,
 R=R,
 S=S,
 mu=mu,
 p0=p0,
 how=how,
 miss=miss,
 lin_keep=lin_keep,
 ret=ret,
 loss_f=loss_f,
 loss_g=loss_g,
 pr_log2=log2(pr),
 surv_log2=log2(survivors),
 sieve_log2=log2(sieve),
 survcost_log2=log2(max(1.0, surv_cost_tot)),
 nf=nf_max,
 ng=ng_max,
 nf_floor=nf_floor,
 ng_floor=ng_floor,
 det=det,
 T=T,
 T0=T0,
 sv0=sv0,
 sv_s1=sv,
 sv_norms=sv_norms,
 sv_fT_g0=sv_fT_g0,
 sv_f0_gT=sv_f0_gT,
 global_floor_s1=True,
 cofac_frac=1.0,
 eager=True,
 )


def evaluate(Pm, PLAN, BOOT):
 """Pure Mode-A k=1 zero-retry unconditional Boole DAG + Unbounded i.i.d. Fresh-e Wald Boot."""
 cells_map = load_cells()
 bins = bins_of(Pm)
 W_NODE, PHI_NODE, C_DET, DET = {}, {}, {}, {}

 def w_max(side, lim):
  return max([W_NODE[b, side] for b in bins if (b, side) in W_NODE and b[0] < lim - 1e-9] + [0.0])

 def phi_max(side, lim):
  return max([PHI_NODE[b, side] for b in bins if (b, side) in PHI_NODE and b[0] < lim - 1e-9] + [0.0])

 def cdet_max(side, lim):
  return max([C_DET[b, side] for b in bins if (b, side) in C_DET and b[0] < lim - 1e-9] + [0])

 for b in bins:
  for side in ("f", "g"):
   n, iA, kA, kB, R = PLAN[b[0], side]
   h_tmp = (n - Pm.LX) / Pm.MESH if n > Pm.LX else 0.0
   nA_tmp = Pm.LX + iA * h_tmp if n > Pm.LX else Pm.LX
   if not ((kA == 0 or nA_tmp < b[0] - 1e-9) and (kB == 0 or n < b[0] - 1e-9)):
    raise ValueError(f"non-decreasing plan at {b} {side}: n={n}, nA={nA_tmp}")
   d = node_cost(Pm, b[0], b[1], side, float(n), iA, kA, kB, R)
   if d is None:
    raise ValueError(f"infeasible plan entry {b} {side}")
   nA = d["nA"]
   c_det = (
   1
   + (kA * (cdet_max("f", nA) + cdet_max("g", nA)) if kA else 0)
   + (kB * (cdet_max("f", float(n)) + cdet_max("g", float(n))) if kB else 0)
   )
   phi_n = (
   d["p0"]
   + (kA * (phi_max("f", nA) + phi_max("g", nA)) if kA else 0.0)
   + (kB * (phi_max("f", float(n)) + phi_max("g", float(n))) if kB else 0.0)
   )
   w_n = (
   d["own_1"]
   + (kA * (w_max("f", nA) + w_max("g", nA)) if kA else 0.0)
   + (kB * (w_max("f", float(n)) + w_max("g", float(n))) if kB else 0.0)
   )
   C_DET[b, side] = c_det
   PHI_NODE[b, side] = phi_n
   W_NODE[b, side] = w_n
   DET[b, side] = dict(
   d,
   c_det=c_det,
   phi_node=phi_n,
   w_node_log2=log2(w_n / UNIT),
   )

 n1, iA_b, kA_b, Rb = BOOT
 if n1 > Pm.LMAX:
  raise ValueError(f"boot window n1={n1} above LMAX={Pm.LMAX}")
 C_Rb = cells_map[Rb]
 Sb = log2(C_Rb)
 cb = (Pm.LOG2P4 + Sb) / 16.0
 hb = (n1 - Pm.LX) / Pm.MESH
 nA_b = Pm.LX + iA_b * hb
 m_max = 24.0 * cb + Pm.LOG2U8 + max(x[1] for x in Pm.BOOT_BINS)
 tau_b = getattr(Pm, "TAU_BOOT", None)
 boot_lo_off = getattr(Pm, "BOOT_LO_OFF", -300.0)
 m_min = 24.0 * cb + Pm.LOG2U8 + boot_lo_off
 if tau_b is None:
  pr_b = sum(
  w * ss_lower2(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
  for w, cg in Pm.BOOT_BINS
  )
  if Pm.BOOT_SUB:
   pri_b = sum(
   ss_lower2(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, Pm.TRUNC)
   for cg in Pm.BOOT_SUB
   )
   pr_b = min(pri_b / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * pr_b)
 else:
  pr_b = sum(
  w * ss_lower2_tau(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, tau_b)
  for w, cg in Pm.BOOT_BINS
  )
  if Pm.BOOT_SUB:
   pri_b = sum(
   ss_lower2_tau(24.0 * cb + Pm.LOG2U8 + cg, Pm.LX, float(n1), iA_b, kA_b, 1, Pm.MESH, Pm.ALPHA_G, tau_b)
   for cg in Pm.BOOT_SUB
   )
   pr_b = min(pri_b / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * pr_b)
 if pr_b <= 0.0:
  raise ValueError("infeasible boot pr_b <= 0")
 boot_n_factors = math.ceil(m_max / 16.0 + kA_b + 1.0)
 det_b = 1.0 - (1.0 - Pm.ECM_SUCCESS) * (m_max / 16.0 + kA_b + 1.0)
 mb = missing(Pm, 0.0, m_max)
 lin_keep_b = max(0.0, 1.0 - mb)
 if lin_keep_b <= 0.0:
  raise ValueError("infeasible boot lin_keep_b <= 0")
 prob_b = pr_b * Pm.BOOT_KEEP_FRAC * det_b * Pm.UDUP_BOOT * _get_counted_fraction(Pm) * lin_keep_b
 if tau_b is None:
  sv_b = 1.0
  trial_b = (
  Pm.LLL_BOOT
  + Pm.NORMS
  + Pm.TDIV
  + ecm_instr(nA_b, m_max, Pm.ECM_SUCCESS)
  + (kA_b + 1) * (Pm.PRIME + Pm.IDEAL)
  )
 else:
  # 2-stage tau-capped K=1 (BALL_PAIR=0.0 exact) boot with explicit two-sided pre-ECM norm window [m_min, m_max]:
  # Before Stage-1 ECM, reject any cell with log2(|N_g|) < m_min or log2(|N_g|) > m_max
  # (acceptance probability = 0.0 identically for all m < m_min).
  # On admitted cells m in [m_min, m_max], Stage-1 ecm(LX, m_max) strips <=69.0b factors,
  # and Stage-2 runs only if residual cofactor log2(R_g) <= tau_b = TAU_BOOT (2095); on admitted norms,
  # Stage-1 survival is monotonically bounded under H1 / H-EA by BOOT_SURV_X * psurv(m_min, LX, tau_b)
  # (IEEE-754 outward-rounded interval Buchstab-Dickman integral * 2.0x safety margin, >2.04x exact P_BD).
  ps_raw_b = psurv(m_min, Pm.LX, tau_b)
  sv_b = min(1.0, float(getattr(Pm, "BOOT_SURV_X", 2.0)) * ps_raw_b)
  ecm_s1 = ecm_instr(Pm.LX, m_max, Pm.ECM_SUCCESS)
  ecm_s2 = ecm_instr(nA_b, tau_b, Pm.ECM_SUCCESS)
  trial_b = (
  Pm.LLL_BOOT
  + Pm.NORMS
  + Pm.TDIV
  + ecm_s1
  + (sv_b + pr_b) * (ecm_s2 + (kA_b + 1) * (Pm.PRIME + Pm.IDEAL))
  )
 boot_own_1run = trial_b / prob_b
 ch_boot = kA_b * w_max("g", nA_b) + w_max("g", float(n1))
 phi_tree = kA_b * phi_max("g", nA_b) + phi_max("g", float(n1))
 c_det_one = kA_b * cdet_max("g", nA_b) + cdet_max("g", float(n1))
 if phi_tree > Pm.MAX_PHI_TREE + 1e-12:
  raise ValueError(f"Phi_tree={phi_tree:.6f} exceeds MAX_PHI_TREE={Pm.MAX_PHI_TREE}")
 e_boot_runs = 1.0 / (1.0 - phi_tree)
 one_log = e_boot_runs * (boot_own_1run + ch_boot)
 two_log_u = log2(2.0 * one_log / UNIT)
 boot_u = log2(2.0 * e_boot_runs * boot_own_1run / UNIT)
 sub_u = log2(2.0 * e_boot_runs * ch_boot / UNIT) if ch_boot > 0 else 0.0
 boot = dict(
 n1=n1,
 nA=nA_b,
 kA=kA_b,
 R=Rb,
 S=Sb,
 m=m_max,
 m_min=m_min,
 tau_boot=tau_b,
 boot_lo_off=boot_lo_off if tau_b is not None else None,
 surv_b=sv_b,
 boot_psurv_raw=ps_raw_b if tau_b is not None else 1.0,
 boot_keep_frac=Pm.BOOT_KEEP_FRAC,
 boot_n_factors=boot_n_factors,
 pr_log2=log2(pr_b),
 trials_1run_log2=log2(1.0 / prob_b),
 trial_log2=log2(trial_b),
 e_trials_log=e_boot_runs / prob_b,
 miss=mb,
 lin_keep=lin_keep_b,
 own_log2=log2(e_boot_runs * boot_own_1run),
 sub_log2=log2(e_boot_runs * ch_boot),
 boot_u=boot_u,
 sub_u=sub_u,
 phi_tree=phi_tree,
 c_det_one=c_det_one,
 e_boot_runs=e_boot_runs,
 p0_max=max(d["p0"] for d in DET.values()),
 mu_min=min(d["mu"] for d in DET.values()),
 phi=0.0,
 one_log_log2=log2(one_log),
 )
 return dict(
 W_NODE=W_NODE,
 PHI_NODE=PHI_NODE,
 C_DET=C_DET,
 DET=DET,
 bins=bins,
 boot=boot,
 two_log_units=two_log_u,
 one_log_instr=one_log,
 )


def report(res):
 rows = []
 for b in res["bins"]:
  for side in ("f", "g"):
   d = res["DET"][b, side]
   rows.append(
   dict(
   bin=b[0],
   side=side,
   n=d["n"],
   nA=round(d["nA"], 2),
   kA=d["kA"],
   kB=d["kB"],
   R=d["R"],
   ret=round(d["ret"], 4),
   loss_f=d["loss_f"],
   loss_g=d["loss_g"],
   how=d["how"],
   mu=round(d["mu"], 2),
   p0=round(d["p0"], 6),
   c_det=d["c_det"],
   phi_node=round(d["phi_node"], 6),
   lin_keep=round(d["lin_keep"], 4),
   own_1_log2=round(d["own_1_log2"], 3),
   w_node_log2=round(d["w_node_log2"], 3),
   nf=round(d["nf"], 1),
   ng=round(d["ng"], 1),
   T=round(d["T"], 1),
   )
   )
 return rows


# Shipped Mode-A pure k=1 zero-retry Boole + E_16 cap C_30 cap F_Vor(a,b) Gram-ellipsoid + Superadditive Kronecker-NTT3 BLS Stage-2 ECM + Lagrangian-DP (Phi<=0.48, tau=2095) Wald boot plan:
PLAN_TAG = "MODEA_SUB107_V5_BIGSTACK_E16_NTT3_PHI48_K1"
_PLAN = {
69.0: ((69.0, 384, 0, 0, 44), (69.0, 384, 0, 0, 44)),
75.0: ((69.0, 384, 0, 0, 47), (69.0, 384, 0, 0, 47)),
80.0: ((69.0, 384, 0, 0, 55), (69.0, 384, 0, 0, 55)),
85.0: ((75.0, 0, 0, 1, 44), (75.0, 0, 0, 1, 45)),
90.0: ((80.0, 0, 0, 1, 45), (80.0, 0, 0, 1, 46)),
95.0: ((80.0, 0, 0, 1, 52), (80.0, 0, 0, 1, 53)),
100.0: ((85.0, 24, 0, 1, 56), (85.0, 384, 1, 0, 56)),
105.0: ((85.0, 384, 1, 0, 61), (80.0, 384, 2, 0, 52)),
110.0: ((80.0, 384, 2, 0, 58), (80.0, 384, 2, 0, 59)),
115.0: ((90.0, 292, 2, 0, 56), (90.0, 292, 2, 0, 58)),
120.0: ((85.0, 264, 2, 1, 54), (85.0, 264, 2, 1, 57)),
125.0: ((85.0, 264, 2, 1, 60), (85.0, 264, 2, 1, 61)),
130.0: ((100.0, 136, 2, 1, 51), (100.0, 136, 2, 1, 55)),
135.0: ((100.0, 136, 2, 1, 59), (100.0, 136, 2, 1, 61)),
140.0: ((100.0, 198, 2, 1, 58), (105.0, 170, 2, 1, 59)),
145.0: ((100.0, 198, 2, 1, 64), (100.0, 198, 3, 1, 55)),
150.0: ((100.0, 198, 3, 1, 55), (100.0, 198, 3, 1, 60)),
155.0: ((100.0, 198, 3, 1, 61), (105.0, 170, 3, 1, 63)),
160.0: ((105.0, 170, 3, 1, 65), (100.0, 384, 4, 0, 57)),
165.0: ((100.0, 384, 4, 0, 55), (105.0, 330, 3, 1, 60)),
170.0: ((105.0, 330, 3, 1, 59), (110.0, 290, 3, 1, 60)),
175.0: ((110.0, 290, 3, 1, 58), (115.0, 258, 3, 1, 61)),
180.0: ((115.0, 258, 3, 1, 59), (120.0, 233, 3, 1, 65)),
185.0: ((120.0, 233, 3, 1, 62), (120.0, 271, 3, 1, 67)),
190.0: ((120.0, 271, 3, 1, 63), (120.0, 233, 4, 1, 59)),
195.0: ((130.0, 226, 3, 1, 62), (120.0, 271, 4, 1, 57)),
200.0: ((140.0, 194, 3, 1, 63), (125.0, 246, 4, 1, 64)),
205.0: ((145.0, 181, 3, 1, 69), (130.0, 226, 4, 1, 67)),
210.0: ((140.0, 248, 3, 1, 63), (145.0, 181, 4, 1, 66)),
215.0: ((140.0, 275, 3, 1, 64), (155.0, 160, 4, 1, 69)),
220.0: ((150.0, 241, 3, 1, 64), (145.0, 232, 4, 1, 61)),
225.0: ((160.0, 215, 3, 1, 65), (140.0, 275, 4, 1, 65)),
230.0: ((160.0, 236, 3, 1, 66), (155.0, 227, 4, 1, 63)),
235.0: ((160.0, 257, 3, 1, 67), (160.0, 215, 4, 1, 69)),
240.0: ((160.0, 278, 3, 1, 67), (160.0, 236, 4, 1, 68)),
245.0: ((160.0, 299, 3, 1, 69), (165.0, 244, 4, 1, 65)),
250.0: ((165.0, 304, 3, 1, 65), (160.0, 215, 5, 1, 68)),
255.0: ((185.0, 235, 3, 1, 67), (185.0, 168, 5, 1, 63)),
}
PLAN16 = {}
for _b0, (_pf, _pg) in _PLAN.items():
 PLAN16[float(_b0), "f"] = _pf
 PLAN16[float(_b0), "g"] = _pg
BOOT16 = (235, 279, 13, 2)


def expected_nodes(Pm, dres):
 bt = dres['boot']
 return 2.0 * bt['e_boot_runs'] * bt['c_det_one'], dres['C_DET']

# --- Exact 16D Euclidean ball half-counts theta_3(q)^16 (merged from q15_cells16) ---
RADII_16D = [
2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20, 24,
28, 32, 40, 44, 45, 46, 47, 48, 49, 50, 51, 52,
53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88,
89, 90, 91, 92, 93, 94, 95, 96, 104, 112, 128, 160,
192, 224, 256,
]


def _poly_sqr_kronecker(coeffs, max_deg):
 """Exact square of non-negative integer polynomial truncated to degree <= max_deg."""
 n = min(len(coeffs), max_deg + 1)
 c = coeffs[:n]
 mx = max(c) if c else 0
 if mx == 0:
  return [0] * (max_deg + 1)
    # Max output coeff <= n * mx^2
 bits = (n * mx * mx).bit_length() + 4
    # Pack via divide-and-conquer binary joining (O(n log n) bit-ops instead of O(n^2))
 def pack(arr):
  m = len(arr)
  if m == 1:
   return arr[0]
  h = m // 2
  return pack(arr[:h]) | (pack(arr[h:]) << (h * bits))

 P = pack(c)
 P2 = P * P
    # Unpack first max_deg + 1 coefficients via divide-and-conquer
 out = [0] * (max_deg + 1)
 def unpack(val, lo, hi):
  m = hi - lo
  if m == 1:
   out[lo] = val
   return
  h = m // 2
  sh = h * bits
  mask = (1 << sh) - 1
  unpack(val & mask, lo, lo + h)
  unpack(val >> sh, lo + h, hi)

 mask_tot = (1 << ((max_deg + 1) * bits)) - 1
 unpack(P2 & mask_tot, 0, max_deg + 1)
 return out


def compute_cells_16d(radii=RADII_16D):
 """Return dict {R: exact half-count of nonzero integer vectors in B_16(R)}."""
 rmax = max(radii)
 nmax = rmax * rmax
    # theta_3(q) = 1 + 2 * sum_{k>=1} q^{k^2}
 t1 = [0] * (nmax + 1)
 t1[0] = 1
 k = 1
 while k * k <= nmax:
  t1[k * k] = 2
  k += 1
 t2 = _poly_sqr_kronecker(t1, nmax)
 t4 = _poly_sqr_kronecker(t2, nmax)
 t8 = _poly_sqr_kronecker(t4, nmax)
    # Only need prefix sums of t16 = t8^2! Still compute exact t16 via Kronecker square.
 t16 = _poly_sqr_kronecker(t8, nmax)

    # Prefix sums (excluding 0-vector -> (sum_{n=0}^{R^2} r_16(n) - 1) // 2)
 out = {}
 acc = 0
 r_set = {r * r: r for r in radii}
 for n, cnt in enumerate(t16):
  acc += cnt
  if n in r_set:
   assert (acc - 1) % 2 == 0
   out[r_set[n]] = (acc - 1) // 2
 return out, t8, t16


# Exact precomputed 16D Euclidean ball half-counts (verified live by verify_cells_16d())
CELLS_16D = {
2: 17072,
3: 7610768,
4: 642328192,
5: 20983317024,
6: 370380040992,
8: 35245626789232,
10: 1224659017256368,
12: 22368828372683552,
14: 261572690513523664,
16: 2204885775774351824,
20: 77892709501111553824,
24: 1435703378669541869120,
28: 16881378095768665289072,
32: 142809330198162449751264,
40: 5066430304652070286222016,
44: 23269941169404213064198608,
45: 33335978184921058854566528,
46: 47380626373836043862806816,
47: 66835388445263268301947808,
48: 93598338771342552134554192,
49: 130171139947638217131397584,
50: 179832662152985340438430416,
51: 246856409932819017251941136,
52: 336782308453559835527614784,
53: 456757285178657459704636224,
54: 615955018736661825271247424,
55: 826095266149915577867513520,
56: 1102084932864321560794312816,
57: 1462800356008087947444783536,
58: 1932042325232195220017376624,
59: 2539706340907815842931956864,
60: 3323187630288119357324959904,
61: 4329086373471401391716596448,
62: 5615272264880597345250781936,
63: 7253346400709206862581807216,
64: 9331585100118065448589077424,
65: 11958494099092118238763422272,
66: 15266988354513029968737665600,
67: 19419375543081442133396096288,
68: 24613288501933780179515846688,
69: 31088638942526637085060714512,
70: 39135809669454454203068585456,
71: 49105406324342637387321483120,
72: 61419488442265210044703617280,
73: 76584824527332961951020493696,
74: 95208479101175351992739707904,
75: 118015754739057731370033840000,
76: 145871021853259038435964024976,
77: 179802333609094776550937275632,
78: 221029325394756883118369094096,
79: 270995717319705887250635226656,
80: 331407089644239582815695476352,
81: 404273867108651477775029192608,
82: 491960716062717087460259314928,
83: 597244362800894441818385164048,
84: 723377881979613983521413539504,
85: 874165336126779543442006244784,
86: 1054047930076009515411124285408,
87: 1268200336846180102996640758816,
88: 1522640172260996956164628542368,
89: 1824355987131552241003996067760,
90: 2181447619495065257460434538256,
91: 2603287600764568839049745195600,
92: 3100705632070257969313444986160,
93: 3686193354248248506731826264256,
94: 4374135475387317443402687358432,
95: 5181078015351411277917626479008,
96: 6126017222460374568589259782640,
104: 22046708351612438386783201872208,
112: 72156655045774710398864896025072,
128: 611103589253768021332109984071136,
160: 21708837938222727070858553213329808,
192: 401342941836214191182349314922514608,
224: 4727718937138462024086156615685585216,
256: 40041895073292555867202006947082587760,
}



# --- Exact outward-rounded integer BOUND_CELL (+6) kernel (merged from q15_bound_cell) ---
U_BLS = -0xd201000000010000
LIM60 = 1 << 60


def _norm_ball(re, im, err, sh):
 """Renormalize integer ball (re, im, err, sh) so |re|, |im|, err < 2^60 (shifting right with outward err rounding if needed)."""
 mx = max(abs(re), abs(im), err)
 if mx == 0:
  return 0, 0, 0, sh
 bl = mx.bit_length()
 if bl >= 60:
  rsh = bl - 59
  # Truncating re and im drops < 2^rsh each => Euclidean error <= sqrt(2)*2^rsh <= 2 * 2^rsh => +2 after shift!
  re2 = re >> rsh if re >= 0 else -((-re) >> rsh)
  im2 = im >> rsh if im >= 0 else -((-im) >> rsh)
  err2 = ((err + (1 << rsh) - 1) >> rsh) + 2
  if err2 >= LIM60:
   re2 = re2 >> 1 if re2 >= 0 else -((-re2) >> 1)
   im2 = im2 >> 1 if im2 >= 0 else -((-im2) >> 1)
   err2 = ((err2 + 1) >> 1) + 2
   rsh += 1
  return re2, im2, err2, sh + rsh
 if bl < 60:
  lsh = 60 - bl
  return re << lsh, im << lsh, err << lsh, sh - lsh
 return re, im, err, sh


def bmul_6(a, b):
 """Exact integer complex ball multiply with +6 truncation/rounding margin (ceil(sqrt(5) + 3) = 6)."""
 x0, x1, x2, s1 = a
 x3, x4, x5, s2 = b
 x8 = ((x0 * x3) >> 64) - ((x1 * x4) >> 64)
 x10 = ((x0 * x4) >> 64) + ((x1 * x3) >> 64)
 x12 = abs(x0) + abs(x1)
 x13 = abs(x3) + abs(x4)
 x14 = ((x12 * x5) >> 64) + ((x13 * x2) >> 64) + ((x2 * x5) >> 64) + 6
 return _norm_ball(x8, x10, x14, s1 + s2 + 64)


def badd(a, b):
 """Exact outward-rounded integer complex ball addition."""
 x0, x1, x2, s1 = a
 x3, x4, x5, s2 = b
 if x0 == 0 and x1 == 0 and x2 == 0:
  return b
 if x3 == 0 and x4 == 0 and x5 == 0:
  return a
 if s1 < s2:
  x0, x1, x2, s1, x3, x4, x5, s2 = x3, x4, x5, s2, x0, x1, x2, s1
 d = s1 - s2
 if d >= 120:
  # |b| <= (|x3|+|x4|+x5) * 2^s2 <= 3 * 2^60 * 2^(s1 - 120) < 2^s1
  return _norm_ball(x0, x1, x2 + 1, s1)
 # Shift b right by d with outward rounding (+2 for re/im truncation + ceil(x5/2^d))
 b_re = x3 >> d if x3 >= 0 else -((-x3) >> d)
 b_im = x4 >> d if x4 >= 0 else -((-x4) >> d)
 b_err = ((x5 + (1 << d) - 1) >> d) + (2 if d > 0 else 0)
 return _norm_ball(x0 + b_re, x1 + b_im, x2 + b_err, s1)


def bsub(a, b):
 return badd(a, (-b[0], -b[1], b[2], b[3]))


def _precompute_q15_balls():
 """High-precision (50-digit Decimal) certified integer balls for zeta_15^(k*j), w_k, w_k+1, U_BLS*w_k^2."""
 with localcontext() as ctx:
  ctx.prec = 50
  # Compute pi to 50 digits via Machin's formula
  def atan_dec(inv_x):
   x2 = Decimal(inv_x * inv_x)
   term = Decimal(1) / Decimal(inv_x)
   s = term
   k = 1
   while True:
    term = -term / x2
    add = term / Decimal(2 * k + 1)
    if abs(add) < Decimal("1e-48"):
     break
    s += add
    k += 1
   return s

  pi = Decimal(16) * atan_dec(5) - Decimal(4) * atan_dec(239)

  def cos_sin(theta):
   t2 = theta * theta
   c = Decimal(1); term = Decimal(1); k = 1
   while True:
    term = -term * t2 / Decimal((2 * k - 1) * (2 * k))
    if abs(term) < Decimal("1e-48"):
     break
    c += term; k += 1
   s = theta; term = theta; k = 1
   while True:
    term = -term * t2 / Decimal((2 * k) * (2 * k + 1))
    if abs(term) < Decimal("1e-48"):
     break
    s += term; k += 1
   return c, s

  scale = Decimal(1 << 58)
  emb = []
  for k in (1, 2, 4, 7):
   pows = []
   for j in range(15):
    ang = Decimal(2) * pi * Decimal((k * j) % 15) / Decimal(15)
    cs, sn = cos_sin(ang)
    re = int((cs * scale).to_integral_value())
    im = int((sn * scale).to_integral_value())
    # Rounded to nearest integer => |err_re|<=0.5, |err_im|<=0.5 => Euclidean err <= 0.708 < 1!
    pows.append((re, im, 1, -58))
   w = pows[5]
   wp1 = badd(w, (1 << 58, 0, 0, -58))
   w2 = pows[10]
   # Exact integer ball for U_BLS * w^2:
   # U_BLS is a 64-bit integer; represent as normalized ball (re, 0, 0, sh)
   u_ball = _norm_ball(U_BLS, 0, 0, 0)
   uw2 = bmul_6(u_ball, w2)
   emb.append((pows[:8], w, wp1, uw2))
  return emb


_EMB4 = _precompute_q15_balls()


def _ball_abs_ub(b):
 """Exact integer upper bound (rad, sh) such that |z| <= rad * 2^sh with rad > 0.
    Uses exact math.isqrt ceiling; hardware ASM uses 9-op SQ_LUT (<= +0.5343b < +0.625b envelope added by callers)."""
 re, im, err, sh = b
 if re == 0 and im == 0 and err == 0:
  return 0, sh
 m2 = re * re + im * im
 sr = math.isqrt(m2)
 if sr * sr < m2:
  sr += 1
 return sr + err, sh


def _eval_lin_ball(vec8, pows):
 """Evaluate sum_{j=0..7} vec8[j] * zeta_k^j as a strictly outward-rounded integer ball."""
 mx = max(abs(x) for x in vec8)
 if mx == 0:
  return (0, 0, 0, 0)
 # Shift vec8 right only if any coeff >= 2^57 so sum_j |c_j| * 2^58 stays exact!
 # Actually in Python arbitrary-precision ints, we compute the exact integer sum at scale 2^-58
 # and normalize once — equivalent to hardware multi-limb accumulation before 1 renormalization!
 re = 0
 im = 0
 err = 0
 for j in range(8):
  cj = vec8[j]
  if cj:
   pr, pi_, pe, _ = pows[j]
   re += cj * pr
   im += cj * pi_
   err += abs(cj) * pe
 return _norm_ball(re, im, err, -58)


def bound_cell_int_fg(aa, bb):
 """Return (nf_ub_bits, ng_ub_bits) using ONLY outward-rounded integer ball arithmetic (+6 margin).

    Strictly guarantees:
      nf_ub_bits >= log2(nf)  and  ng_ub_bits >= log2(ng)  (whenever nf, ng > 0).
    For zero norm (rad == 0), returns (-1.0, ...) so callers detect exact zero.
    Also returns exact integer ceil_log2 upper bounds (nf_ceil, ng_ceil) via integer bit_length!
    """
 prod_f = 1
 sh_f = 0
 prod_g = 1
 sh_g = 0
 for pows, w, wp1, uw2 in _EMB4:
  A = _eval_lin_ball(aa, pows)
  B = _eval_lin_ball(bb, pows)
  A2 = bmul_6(A, A)
  A3 = bmul_6(A2, A)
  B2 = bmul_6(B, B)
  B3 = bmul_6(B2, B)
  # g_3: A^3 - U_BLS * w^2 * B^3
  G = bsub(A3, bmul_6(uw2, B3))
  rg, sg = _ball_abs_ub(G)
  if rg == 0:
   prod_g = 0
  else:
   prod_g *= rg * rg
   sh_g += 2 * sg
  # f_9: A^9 + (w+1) A^6 B^3 + A^3 B^6 - w B^9
  A6 = bmul_6(A3, A3)
  B6 = bmul_6(B3, B3)
  A9 = bmul_6(A6, A3)
  B9 = bmul_6(B6, B3)
  t1 = bmul_6(wp1, bmul_6(A6, B3))
  t2 = bmul_6(A3, B6)
  t3 = bmul_6(w, B9)
  F = bsub(badd(badd(A9, t1), t2), t3)
  rf, sf = _ball_abs_ub(F)
  if rf == 0:
   prod_f = 0
  else:
   prod_f *= rf * rf
   sh_f += 2 * sf

 # Exact integer ceiling of log2(prod * 2^sh):
 # For any integer P > 0: 2^(P.bit_length()-1) <= P <= 2^P.bit_length(), and strictly P == 2^(bl-1) iff (P & (P-1)) == 0.
 def _ub_bits(P, sh):
  if P == 0:
   return -1.0, 0
  bl = P.bit_length()
  c_int = (bl - 1 if (P & (P - 1)) == 0 else bl) + sh
  # Tight continuous upper bound (outward-rounded by +1e-12 for IEEE-754 log2):
  f_ub = math.log2(P) + sh + 1e-12
  return f_ub, c_int

 nf_ub, nf_ceil = _ub_bits(prod_f, sh_f)
 ng_ub, ng_ceil = _ub_bits(prod_g, sh_g)
 return nf_ub, ng_ub, nf_ceil, ng_ceil


def bound_cell_check_and_cap(aa, bb, nf_max, ng_max, lq=0.0, side=None, hw_env=0.625, nf_floor=None, ng_floor=None):
 """Deterministic per-cell BOUND_CELL (+6) evaluation WITH hard (nf_max, ng_max) norm-cap discard
    and optional pre-sieve/pre-NORMS lower-floor screen (since nf_rem >= log2(|N_f|/q_f) unconditionally,
    any cell with nf_rem < nf_floor or ng_rem < ng_floor is guaranteed sub-floor and discarded before NORMS).

    Guarantees ZERO cells exceeding (nf_max, ng_max) ever enter cofactorization or ECM (so ECM limb
    counts never exceed ceil(nf_max / 64) or ceil(ng_max / 64) on ANY input lattice!).
    Returns (admitted: bool, nf_rem_ub: float, ng_rem_ub: float, nf_ub: float, ng_ub: float).
    """
 nf_ub, ng_ub, _, _ = bound_cell_int_fg(aa, bb)
 if nf_ub < 0.0 or ng_ub < 0.0:
  return False, -1.0, -1.0, nf_ub, ng_ub
 nf_rem = (nf_ub + hw_env) - (lq if side == "f" else 0.0)
 ng_rem = (ng_ub + hw_env) - (lq if side == "g" else 0.0)
 if nf_rem > nf_max or ng_rem > ng_max:
  return False, nf_rem, ng_rem, nf_ub, ng_ub
 if (nf_floor is not None and nf_rem < nf_floor) or (ng_floor is not None and ng_rem < ng_floor):
  return False, nf_rem, ng_rem, nf_ub, ng_ub
 return True, nf_rem, ng_rem, nf_ub, ng_ub


def screen_special_q_basis(B_red, q, side, norms_fn):
 """Deterministic post-LLL basis screening guard: reject any special-q lattice whose reduced basis
    contains a rank-deficient q^2-or-higher short-subspace basis vector (v_q(N_side(B_red[i])) != 1, e.g.
    q^4 short subspace when t = a/b mod q with small a, b). Costs 16 norm checks (< 2^31 Ir << LLL_NODE = 2^40 Ir);
    rejected bases are discarded under DISC. Individual surviving combinations with v_q(N_side) >= 2 are
    deterministically rejected during exact norm/cofactor verification via (N_side // q) % q == 0."""
 q2 = q * q
 for row in B_red:
  nf, ng = norms_fn(row[:8], row[8:])
  ns = nf if side == "f" else ng
  if ns == 0 or (ns % q) != 0 or (ns % q2) == 0:
   return False
 return True


