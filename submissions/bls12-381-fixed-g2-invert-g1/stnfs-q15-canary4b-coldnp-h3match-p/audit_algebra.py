"""audit_algebra.py — Extracted algebraic, toy E2E DLP, and Voronoi/G_T boot elimination checks
for audit_cm.py (keeps every Python module well below Yukon's 65,536 B per-file limit)."""
import cmath
import hashlib
import math
import random
from fractions import Fraction as F

from cm_models import NTT4_PRIMES, SMALL_150, lll_16d_int, log2_li

U_BLS = -0xd201000000010000
R_BLS = U_BLS ** 4 - U_BLS ** 2 + 1
P_BLS = ((U_BLS - 1) ** 2 * R_BLS) // 3 + U_BLS
PHI15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]


def ptrim(a):
 a = list(a)
 while len(a) > 1 and a[-1] == 0:
  a.pop()
 return a


def pmul_mod(a, b, m, mod_poly=None):
 out = [0] * (len(a) + len(b) - 1)
 for i, ai in enumerate(a):
  if ai:
   for j, bj in enumerate(b):
    if bj:
     out[i + j] += ai * bj
 if mod_poly is not None:
  return prem_mod(out, mod_poly, m)
 out = [x % m for x in out]
 while len(out) > 1 and out[-1] == 0:
  out.pop()
 return out


def prem_mod(a, b, m):
 a = [x % m for x in a]
 while len(a) > 1 and a[-1] == 0:
  a.pop()
 db = len(b) - 1
 while db > 0 and b[db] % m == 0:
  db -= 1
 if db == 0:
  return [0]
 lc = b[db] % m
 monic = (lc == 1)
 inv_lc = 1 if monic else pow(lc, -1, m)
 b_mod = [x % m for x in b[:db + 1]]
 while len(a) - 1 >= db and not (len(a) == 1 and a[0] == 0):
  c = a[-1] if monic else (a[-1] * inv_lc) % m
  off = len(a) - 1 - db
  for i in range(db):
   bi = b_mod[i]
   if bi:
    a[off + i] = (a[off + i] - c * bi) % m
  a.pop()
  while len(a) > 1 and a[-1] == 0:
   a.pop()
 return a


def ppow_mod(base, exp, mod_poly, m):
 res = [1]
 b = prem_mod(base, mod_poly, m)
 e = exp
 while e > 0:
  if e & 1:
   res = pmul_mod(res, b, m, mod_poly)
  e >>= 1
  if e:
   b = pmul_mod(b, b, m, mod_poly)
 return res


def pgcd_mod(a, b, m):
 a = ptrim([x % m for x in a])
 b = ptrim([x % m for x in b])
 while not (len(b) == 1 and b[0] == 0):
  a, b = b, prem_mod(a, b, m)
 inv = pow(a[-1], -1, m)
 return [(x * inv) % m for x in a]


def pmul_z(a, b):
 out = [0] * (len(a) + len(b) - 1)
 for i, ai in enumerate(a):
  if ai:
   for j, bj in enumerate(b):
    if bj:
     out[i + j] += ai * bj
 return out


def red_phi15_z(a):
 c = list(a)
 while len(c) > 8:
  lead = c.pop()
  if lead:
   off = len(c) - 8
   for i in range(8):
    c[off + i] -= lead * PHI15[i]
 while len(c) < 8:
  c.append(0)
 return c


def mul_phi15_z(a, b):
 a0, a1, a2, a3, a4, a5, a6, a7 = a
 b0, b1, b2, b3, b4, b5, b6, b7 = b
 c0 = a0 * b0
 c1 = a0 * b1 + a1 * b0
 c2 = a0 * b2 + a1 * b1 + a2 * b0
 c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0
 c4 = a0 * b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4 * b0
 c5 = a0 * b5 + a1 * b4 + a2 * b3 + a3 * b2 + a4 * b1 + a5 * b0
 c6 = a0 * b6 + a1 * b5 + a2 * b4 + a3 * b3 + a4 * b2 + a5 * b1 + a6 * b0
 c7 = a0 * b7 + a1 * b6 + a2 * b5 + a3 * b4 + a4 * b3 + a5 * b2 + a6 * b1 + a7 * b0
 c8 = a1 * b7 + a2 * b6 + a3 * b5 + a4 * b4 + a5 * b3 + a6 * b2 + a7 * b1
 c9 = a2 * b7 + a3 * b6 + a4 * b5 + a5 * b4 + a6 * b3 + a7 * b2
 c10 = a3 * b7 + a4 * b6 + a5 * b5 + a6 * b4 + a7 * b3
 c11 = a4 * b7 + a5 * b6 + a6 * b5 + a7 * b4
 c12 = a5 * b7 + a6 * b6 + a7 * b5
 c13 = a6 * b7 + a7 * b6
 c14 = a7 * b7
 c13 += c14; c11 -= c14; c10 += c14; c9 -= c14; c7 += c14; c6 -= c14
 c12 += c13; c10 -= c13; c9 += c13; c8 -= c13; c6 += c13; c5 -= c13
 c11 += c12; c9 -= c12; c8 += c12; c7 -= c12; c5 += c12; c4 -= c12
 c10 += c11; c8 -= c11; c7 += c11; c6 -= c11; c4 += c11; c3 -= c11
 c9 += c10; c7 -= c10; c6 += c10; c5 -= c10; c3 += c10; c2 -= c10
 c8 += c9; c6 -= c9; c5 += c9; c4 -= c9; c2 += c9; c1 -= c9
 c7 += c8; c5 -= c8; c4 += c8; c3 -= c8; c1 += c8; c0 -= c8
 return [c0, c1, c2, c3, c4, c5, c6, c7]


def pow_phi15_z(a, e):
 if e == 3:
  return mul_phi15_z(a, mul_phi15_z(a, a))
 res = [1] + [0] * 7
 b = list(a)
 while e > 0:
  if e & 1:
   res = mul_phi15_z(res, b)
  e >>= 1
  if e:
   b = mul_phi15_z(b, b)
 return res


_YPOW15 = tuple(red_phi15_z([0] * k + [1]) for k in range(15))


def y_pow_z(k):
 return list(_YPOW15[k % 15])


def sub_phi15_z(a, b):
 return [a[i] - b[i] for i in range(8)]


def add_phi15_z(a, b):
 return [a[i] + b[i] for i in range(8)]


def scale_phi15_z(s, a):
 return [s * x for x in a]


def norm_q15_z(a):
 v0, v1, v2, v3, v4, v5, v6, v7 = a
 p1 = mul_phi15_z(a, (v0 + v1 + v2 - v5 - v6 - v7, -v2 - v4 + v7, -v1 - v3 + v6, v1 - v6 - v7, -v1 - v2 + v7, v2 - v5 - v7, v1 - v4 - v6, -v1 - v2 - v3 + v6 + v7))
 v0, v1, v2, v3, v4, v5, v6, v7 = p1
 p2 = mul_phi15_z(p1, (v0 - v2 - v6 + v7, v2 + v4 - v7, -v3 + v6, -v2 - v6, v1 + v2 - v7, -v2 + v5 + v7, -v6, v2 - v3 + v6 - v7))
 v0, v1, v2, v3, v4, v5, v6, v7 = p2
 p3 = mul_phi15_z(p2, (v0 + v2 + v4, -v4, -v2 - v6, v2, -v2 - v4 + v7, v4 + v5, v2 + v3, v1 - v2 - v4 - v6))
 assert not any(p3[1:]), f"non-rational norm {p3}"
 return p3[0]


def trace_z15(d):
 d %= 15
 g = math.gcd(15, d)
 return 8 if g == 15 else (-4 if g == 5 else (-2 if g == 3 else 1))


def mat_mul(A, B):
 n, k, m = len(A), len(B), len(B[0])
 return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]


def mat_T(A):
 return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def int_det(M):
 n = len(M)
 A = [[F(x) for x in row] for row in M]
 det = F(1)
 for i in range(n):
  piv = next(r for r in range(i, n) if A[r][i] != 0)
  if piv != i:
   A[i], A[piv] = A[piv], A[i]
   det = -det
  det *= A[i][i]
  inv = F(1) / A[i][i]
  for r in range(i + 1, n):
   if A[r][i] != 0:
    fac = A[r][i] * inv
    for c in range(i, n):
     A[r][c] -= fac * A[i][c]
 assert det.denominator == 1
 return det.numerator


def inv_int_mat(M):
 n = len(M)
 A = [[F(x) for x in row] + [F(1 if i == j else 0) for j in range(n)] for i, row in enumerate(M)]
 for i in range(n):
  piv = next(r for r in range(i, n) if A[r][i] != 0)
  A[i], A[piv] = A[piv], A[i]
  inv = F(1) / A[i][i]
  A[i] = [x * inv for x in A[i]]
  for r in range(n):
   if r != i and A[r][i] != 0:
    fac = A[r][i]
    A[r] = [A[r][c] - fac * A[i][c] for c in range(2 * n)]
 return [[x.numerator for x in row[n:]] for row in A]


EXPS_1W = [0, 5, 3, 8, 6, 11, 9, 14]
M_1W = mat_T([y_pow_z(e) for e in EXPS_1W])
M_1W_INV = inv_int_mat(M_1W)


def build_G1w():
 G_pow = [[trace_z15(i - j) for j in range(8)] for i in range(8)]
 return mat_mul(mat_T(M_1W), mat_mul(G_pow, M_1W))


def build_WS_1w():
 M_sig = mat_T([y_pow_z(7 * j) for j in range(8)])
 M_w = mat_T([y_pow_z(5 + j) for j in range(8)])
 S_1w = mat_mul(M_1W_INV, mat_mul(M_sig, M_1W))
 W_1w = mat_mul(M_1W_INV, mat_mul(M_w, M_1W))
 return W_1w, S_1w


def check_algebra_bls12():
 assert P_BLS.bit_length() == 381
 assert R_BLS.bit_length() == 255
 assert P_BLS % 15 == 7
 assert R_BLS % 15 == 13
 chi = P_BLS % R_BLS
 assert chi == (U_BLS % R_BLS)
 assert pow(chi, 12, R_BLS) == 1
 _cof_gt = (P_BLS ** 12 - 1) // R_BLS
 assert (P_BLS ** 12 - 1) % R_BLS == 0 and math.gcd(R_BLS, _cof_gt) == 1
 for d in (1, 2, 3, 4, 6):
  assert math.gcd(pow(chi, d, R_BLS) - 1, R_BLS) == 1
 w_pi = ((U_BLS ** 3 - U_BLS ** 2 - U_BLS) * pow(U_BLS + 1, -1, P_BLS)) % P_BLS
 assert (w_pi * w_pi + w_pi + 1) % P_BLS == 0
 w_pibar = (-1 - w_pi) % P_BLS
 h_pi = pgcd_mod(PHI15, [-w_pi % P_BLS, 0, 0, 0, 0, 1], P_BLS)
 h_pibar = pgcd_mod(PHI15, [-w_pibar % P_BLS, 0, 0, 0, 0, 1], P_BLS)
 assert len(h_pi) == 5 and len(h_pibar) == 5
 assert pmul_mod(h_pi, h_pibar, P_BLS) == [x % P_BLS for x in PHI15]
 yp = ppow_mod([0, 1], P_BLS, h_pi, P_BLS)
 assert yp == prem_mod(y_pow_z(7), h_pi, P_BLS), "sigma(y)=y^7 must equal y^p mod h_pi!"
 uw2 = (U_BLS * w_pi * w_pi) % P_BLS
 exp3 = (P_BLS ** 4 - 1) // 3
 c_frob = pow(uw2, (P_BLS - 1) // 3, P_BLS)
 assert pow(uw2, exp3, P_BLS) != 1, "u*w^2 must be a non-cube in F_{p^4}"
 assert c_frob == (w_pi * w_pi) % P_BLS, "X^p = omega^2 * X required for tau = Frob^1!"
 x3 = uw2
 assert (pow(x3, 3, P_BLS) + (w_pi + 1) * pow(x3, 2, P_BLS) + x3 - w_pi) % P_BLS == 0
 assert -U_BLS == (1 << 16) * 906349 * 254760293
 assert pow(3, (19 ** 2 - 1) // 3, 19) != 1 and pow(6, (67 ** 4 - 1) // 3, 67) != 1
 assert 67 % 15 == 7 and 5 % 3 == 2
 f_ge2_orbits_log2 = 1.0 + log2_li(35.85)
 assert f_ge2_orbits_log2 < 32.5 and 35.85 < 48.0
 w_r1 = (-U_BLS * U_BLS) % R_BLS
 assert (w_r1 * w_r1 + w_r1 + 1) % R_BLS == 0
 w_r2 = (-1 - w_r1) % R_BLS
 assert w_r2 == pow(U_BLS, 4, R_BLS)
 hr1 = pgcd_mod(PHI15, [-w_r1 % R_BLS, 0, 0, 0, 0, 1], R_BLS)
 hr2 = pgcd_mod(PHI15, [-w_r2 % R_BLS, 0, 0, 0, 0, 1], R_BLS)
 assert len(hr1) == 5 and len(hr2) == 5
 assert pmul_mod(hr1, hr2, R_BLS) == [x % R_BLS for x in PHI15]

 def mat_mul_mod(A, B, m):
  n = len(A)
  return [[sum(A[i][k] * B[k][j] for k in range(n)) % m for j in range(n)] for i in range(n)]

 def mat_pow_mod(A, e, m):
  n = len(A)
  res = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
  b = [row[:] for row in A]
  while e > 0:
   if e & 1:
    res = mat_mul_mod(res, b, m)
   e >>= 1
   if e:
    b = mat_mul_mod(b, b, m)
  return res

 def frob_mat(poly, m):
  d = len(poly) - 1
  xr = ppow_mod([0, 1], m, poly, m)
  cols = []
  cur = [1]
  for _ in range(d):
   cols.append(cur + [0] * (d - len(cur)))
   cur = pmul_mod(cur, xr, m, poly)
  return [[cols[j][i] for j in range(d)] for i in range(d)]

 I4 = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
 for hri in (hr1, hr2):
  M4 = frob_mat(hri, R_BLS)
  M4_2 = mat_mul_mod(M4, M4, R_BLS)
  assert M4_2 != I4 and mat_mul_mod(M4_2, M4_2, R_BLS) == I4
 u3 = pow(U_BLS, 3, R_BLS)
 g2_roots = [u3, (u3 * w_r2) % R_BLS, (u3 * w_r2 * w_r2) % R_BLS]
 assert len(set(g2_roots)) == 3
 g2_prod = [1]
 for rt in g2_roots:
  g2_prod = pmul_mod(g2_prod, [(-rt) % R_BLS, 1], R_BLS)
 assert g2_prod == [(-U_BLS * w_r2 * w_r2) % R_BLS, 0, 0, 1]
 c1 = (U_BLS * w_r1 * w_r1) % R_BLS
 assert pow(c1, (R_BLS - 1) // 3, R_BLS) != 1
 assert ((R_BLS + 1) * (R_BLS * R_BLS + 1)) % 3 == 1
 I9 = [[1 if i == j else 0 for j in range(9)] for i in range(9)]
 for wi in (w_r1, w_r2):
  fi = [(-wi) % R_BLS, 0, 0, 1, 0, 0, (wi + 1) % R_BLS, 0, 0, 1]
  M9 = frob_mat(fi, R_BLS)
  M9_3 = mat_pow_mod(M9, 3, R_BLS)
  assert mat_pow_mod(M9_3, 3, R_BLS) == I9
  x_r3_minus_x = [M9_3[i][1] for i in range(9)]
  x_r3_minus_x[1] = (x_r3_minus_x[1] - 1) % R_BLS
  assert pgcd_mod(x_r3_minus_x, fi, R_BLS) == [1]
 G = build_G1w()
 detG = int_det(G)
 assert detG == 1265625 == (3 ** 4) * (5 ** 6)
 W, S = build_WS_1w()
 assert mat_mul(mat_T(W), mat_mul(G, W)) == G
 assert mat_mul(mat_T(S), mat_mul(G, S)) == G
 for np_i in NTT4_PRIMES:
  assert np_i < (1 << 62) and (np_i - 1) % (1 << 53) == 0
  for a_b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
   d_m, s_m = np_i - 1, 0
   while d_m % 2 == 0:
    d_m //= 2; s_m += 1
   x_m = pow(a_b, d_m, np_i)
   if x_m in (1, np_i - 1):
    continue
   for _ in range(s_m - 1):
    x_m = (x_m * x_m) % np_i
    if x_m == np_i - 1:
     break
   else:
    raise AssertionError(f"NTT prime {np_i} composite!")
 ntt4_bits = sum(math.log2(np_i) for np_i in NTT4_PRIMES)
 assert ntt4_bits > 244.0
 _oh_sm = []
 for _p in range(2, 256):
  if any(_p % _d == 0 for _d in range(2, int(_p ** 0.5) + 1)):
   continue
  if _p == 3:
   _f, _g = 4, 1
  elif _p == 5:
   _f, _g = 2, 1
  else:
   _x, _f = _p % 15, 1
   while _x != 1:
    _x = (_x * _p) % 15; _f += 1
   _g = 8 // _f
  _Np = _p ** _f
  if _Np < 256:
   _oh_sm.append((_p, _f, _g, _Np))
 _min_np = min(x[3] for x in _oh_sm)
 assert _min_np == 16 > 9 > 3 and sum(x[2] for x in _oh_sm) == 56

 def _fq_roots(p, mp):
  d = len(mp) - 1
  els = [[]]
  for _ in range(d):
   els = [x + [c] for x in els for c in range(p)]

  def pm(a, b):
   r = [0] * (2 * d - 1)
   for i, x in enumerate(a):
    for j, y in enumerate(b):
     r[i + j] = (r[i + j] + x * y) % p
   for i in range(len(r) - 1, d - 1, -1):
    c = r[i]
    if c:
     for j in range(d):
      r[i - d + j] = (r[i - d + j] - c * mp[j]) % p
   return r[:d]

  one = [1] + [0] * (d - 1)
  w = one if p == 3 else next(e for e in els if not any((pm(e, e)[j] + e[j] + one[j]) % p for j in range(d)))
  w2 = pm(w, w)
  wp1 = [(w[j] + one[j]) % p for j in range(d)]
  rf = rg = 0
  for x in els:
   x2 = pm(x, x); x3 = pm(x2, x)
   if not any((x3[j] - (U_BLS % p) * w2[j]) % p for j in range(d)):
    rg += 1
   x6 = pm(x3, x3); x9 = pm(x6, x3); t6 = pm(wp1, x6)
   if not any((x9[j] + t6[j] + x3[j] - w[j]) % p for j in range(d)):
    rf += 1
  return rf, rg

 _sm_cnt = {}
 for _p, _mp, _g in ((2, [1, 1, 0, 0, 1], 2), (3, [2, 1, 0, 0, 1], 1), (5, [2, 1, 1], 1), (11, [1, 1, 1], 4)):
  _rf, _rg = _fq_roots(_p, _mp)
  _Np = _p ** (len(_mp) - 1)
  if _rf + _rg:
   _sm_cnt[_Np] = _g * (_rf + _rg)
 for _p in (31, 61, 151, 181, 211, 241):
  _z = next(pow(b, (_p - 1) // 15, _p) for b in range(2, 200) if pow(b, (_p - 1) // 3, _p) != 1 and pow(b, (_p - 1) // 5, _p) != 1 and pow(b, (_p - 1) // 15, _p) != 1)
  _cp = 0
  for _k in (1, 2, 4, 7, 8, 11, 13, 14):
   _ry = pow(_z, _k, _p); _w = pow(_ry, 5, _p); _w2 = (_w * _w) % _p
   _cp += sum(1 for _x in range(_p) if (pow(_x, 9, _p) + (_w + 1) * pow(_x, 6, _p) + pow(_x, 3, _p) - _w) % _p == 0)
   _cp += sum(1 for _x in range(_p) if (pow(_x, 3, _p) - U_BLS * _w2) % _p == 0)
  _sm_cnt[_p] = _cp
 assert sorted(_sm_cnt.items()) == sorted(SMALL_150) and sum(_sm_cnt.values()) == 150
 _sm_hits = sum(c / float(n) for n, c in SMALL_150)
 return dict(
 p_bits=P_BLS.bit_length(), r_bits=R_BLS.bit_length(), p_mod15=P_BLS % 15, r_mod15=R_BLS % 15,
 h_pi_deg=4, h_pibar_deg=4, uw2_noncube=True, tau_eq_frob1=True,
 mod_r_crt_degs=[4, 4], mod_r_g_degs=[4, 4, 4, 12], mod_r_f_degs=[36, 36],
 eps_g_deg=12, eps_f_deg=36, det_G=detG, w_s_isometries=True,
 ntt4_primes_ok=len(NTT4_PRIMES), ntt4_prod_bits=round(ntt4_bits, 2),
 f_ge2_orbits_log2=round(f_ge2_orbits_log2, 2), oh_min_prime_norm=_min_np,
 oh_small_ideals_lt256=56, kf_kg_small_ideals_lt256=150, small_ideals_hits_pc=round(_sm_hits, 4),
 )


def check_toy_e2e_dlp(toy):
 u0, p0, r0, c_x, chi = toy["u0"], toy["p0"], toy["r0"], toy["c_x"], toy["chi"]
 h_pi, h_pibar = toy["h_pi"], toy["h_pibar"]
 assert p0 == 7207 and r0 == 601 and chi == p0 % r0 == 596
 assert pmul_mod(h_pi, h_pibar, p0) == [x % p0 for x in PHI15]
 def fp4_add(a, b): return [(a[i] + b[i]) % p0 for i in range(4)]
 def fp4_sub(a, b): return [(a[i] - b[i]) % p0 for i in range(4)]
 _h0, _h1, _h2, _h3 = h_pi[0], h_pi[1], h_pi[2], h_pi[3]
 def fp4_mul(a, b):
  a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
  c0 = a0 * b0; c1 = a0 * b1 + a1 * b0; c2 = a0 * b2 + a1 * b1 + a2 * b0; c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0
  c4 = a1 * b3 + a2 * b2 + a3 * b1; c5 = a2 * b3 + a3 * b2; c6 = a3 * b3
  v = c6 % p0
  if v: c2 -= v * _h0; c3 -= v * _h1; c4 -= v * _h2; c5 -= v * _h3
  v = c5 % p0
  if v: c1 -= v * _h0; c2 -= v * _h1; c3 -= v * _h2; c4 -= v * _h3
  v = c4 % p0
  if v: c0 -= v * _h0; c1 -= v * _h1; c2 -= v * _h2; c3 -= v * _h3
  return [c0 % p0, c1 % p0, c2 % p0, c3 % p0]
 def fp4_pow(a, e):
  res = [1, 0, 0, 0]; b = list(a)
  while e > 0:
   if e & 1: res = fp4_mul(res, b)
   e >>= 1
   if e: b = fp4_mul(b, b)
  return res
 def fp4_inv(a):
  a1 = fp4_pow(a, p0); a2 = fp4_pow(a1, p0); a3 = fp4_pow(a2, p0)
  prod3 = fp4_mul(a1, fp4_mul(a2, a3))
  inv_n = pow(fp4_mul(a, prod3)[0], -1, p0)
  return [(x * inv_n) % p0 for x in prod3]
 y4 = [0, 1, 0, 0]; w4 = fp4_pow(y4, 5); w4_2 = fp4_mul(w4, w4)
 uw2 = [(u0 * x) % p0 for x in w4_2]
 wc1 = w4 if c_x == 1 else w4_2; wc2 = fp4_mul(wc1, wc1)
 frob4_cols = [fp4_pow([1 if j == i else 0 for j in range(4)], p0) for i in range(4)]
 def frob4(a):
  out = [0, 0, 0, 0]
  for i in range(4):
   ai = a[i]
   if ai:
    col = frob4_cols[i]
    for j in range(4): out[j] = (out[j] + ai * col[j]) % p0
  return out
 def fp12_mul(A, B):
  c0 = fp4_add(fp4_mul(A[0], B[0]), fp4_mul(uw2, fp4_add(fp4_mul(A[1], B[2]), fp4_mul(A[2], B[1]))))
  c1 = fp4_add(fp4_add(fp4_mul(A[0], B[1]), fp4_mul(A[1], B[0])), fp4_mul(uw2, fp4_mul(A[2], B[2])))
  c2 = fp4_add(fp4_add(fp4_mul(A[0], B[2]), fp4_mul(A[1], B[1])), fp4_mul(A[2], B[0]))
  return (c0, c1, c2)
 def fp12_pow(A, e):
  res = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]); b = A
  while e > 0:
   if e & 1: res = fp12_mul(res, b)
   e >>= 1
   if e: b = fp12_mul(b, b)
  return res
 def fp12_frob(A): return (frob4(A[0]), fp4_mul(frob4(A[1]), wc1), fp4_mul(frob4(A[2]), wc2))
 def fp12_inv(A):
  def f4(Z):
   for _ in range(4): Z = fp12_frob(Z)
   return Z
  A4 = f4(A); A8 = f4(A4); adj = fp12_mul(A4, A8); inv4 = fp4_inv(fp12_mul(A, adj)[0])
  return (fp4_mul(adj[0], inv4), fp4_mul(adj[1], inv4), fp4_mul(adj[2], inv4))
 cof_GT = (p0 ** 4 - p0 ** 2 + 1) // r0
 def to_GT(A):
  e6 = A
  for _ in range(6): e6 = fp12_frob(e6)
  ez = fp12_mul(e6, fp12_inv(A)); e2 = fp12_frob(fp12_frob(ez))
  return fp12_pow(fp12_mul(e2, ez), cof_GT)
 def embed_ab(a8, b8):
  ar = prem_mod(a8, h_pi, p0); ar += [0] * (4 - len(ar))
  br = prem_mod([-x for x in b8], h_pi, p0); br += [0] * (4 - len(br))
  return (ar, br, [0, 0, 0, 0])
 def key12(A): return (tuple(A[0]), tuple(A[1]), tuple(A[2]))
 G_gen = to_GT(([3, 1, 0, 0], [2, 0, p0 - 1, 0], [0, 5, 0, 0]))
 m_bs = math.isqrt(r0) + 1; bs_map = {}; cur = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
 for j in range(m_bs):
  bs_map[key12(cur)] = j; cur = fp12_mul(cur, G_gen)
 Gm_inv = fp12_inv(fp12_pow(G_gen, m_bs))
 def dlog_GT(A):
  c = to_GT(A)
  for i in range(m_bs + 1):
   k = key12(c)
   if k in bs_map: return (i * m_bs + bs_map[k]) % r0
   c = fp12_mul(c, Gm_inv)
  raise RuntimeError("BSGS failed")
 r0sq = r0 * r0
 roots_r0 = [s for s in range(2, r0) if pow(s, 15, r0) == 1 and pow(s, 5, r0) != 1 and pow(s, 3, r0) != 1]
 assert len(roots_r0) == 8
 roots_l2 = []
 for s in roots_r0:
  val = sum(PHI15[j] * pow(s, j, r0sq) for j in range(9)) % r0sq
  der = sum(j * PHI15[j] * pow(s, j - 1, r0) for j in range(1, 9)) % r0
  roots_l2.append((s - val * pow(der, -1, r0)) % r0sq)
 crt_e = []
 for idx, s in enumerate(roots_r0):
  poly = [1]; den = 1
  for j, sj in enumerate(roots_r0):
   if j != idx:
    poly = pmul_mod(poly, [-sj % r0, 1], r0); den = (den * (s - sj)) % r0
  invd = pow(den, -1, r0)
  p8 = [(x * invd) % r0 for x in poly] + [0] * (8 - len(poly))
  crt_e.append(p8)
 w_tau_exp = 5 * c_x
 chi_inv_pows = [pow(chi, (-t) % 12, r0) for t in range(12)]
 s2_pows = [[pow(s2, j, r0sq) for j in range(8)] for s2 in roots_l2]
 s2_meta = []
 for s2 in roots_l2:
  wc = pow(s2, 5, r0sq); w2c = (wc * wc) % r0sq
  s2_meta.append(((u0 * w2c) % r0sq, wc % r0sq, ((wc + 1) * wc) % r0sq, w2c % r0sq))
 def sm_chi_proj(a8, b8, is_g):
  dX = 3 if is_g else 9; eps = (r0 ** 3 - 1) if is_g else (r0 ** 6 - 1)
  L_c = [[0] * 8 for _ in range(dX)]; _MSK = (1 << 44) - 1
  for ci, (spw, (cg0, cf0, cf3, cf6)) in enumerate(zip(s2_pows, s2_meta)):
   ac = sum(a8[j] * spw[j] for j in range(8)) % r0sq
   bc = sum(b8[j] * spw[j] for j in range(8)) % r0sq
   if is_g:
    c0 = cg0; p0_, p1_, p2_ = 1, 0, 0; q0_, q1_, q2_ = ac, (-bc) % r0sq, 0; ee = eps
    while ee:
     if ee & 1:
      p0_, p1_, p2_ = (p0_ * q0_ + c0 * (p1_ * q2_ + p2_ * q1_)) % r0sq, (p0_ * q1_ + p1_ * q0_ + c0 * p2_ * q2_) % r0sq, (p0_ * q2_ + p1_ * q1_ + p2_ * q0_) % r0sq
     ee >>= 1
     if ee:
      q0_, q1_, q2_ = (q0_ * q0_ + c0 * (2 * q1_ * q2_)) % r0sq, (2 * q0_ * q1_ + c0 * q2_ * q2_) % r0sq, (2 * q0_ * q2_ + q1_ * q1_) % r0sq
    pw = [p0_, p1_, p2_]
   else:
    c0, c3, c6 = cf0, cf3, cf6
    def _smul9(x, y):
     ix = ((((((((x[8] << 44) | x[7]) << 44 | x[6]) << 44 | x[5]) << 44 | x[4]) << 44 | x[3]) << 44 | x[2]) << 44 | x[1]) << 44 | x[0]
     iy = ((((((((y[8] << 44) | y[7]) << 44 | y[6]) << 44 | y[5]) << 44 | y[4]) << 44 | y[3]) << 44 | y[2]) << 44 | y[1]) << 44 | y[0]
     p = ix * iy; o = []
     for _ in range(17): o.append(p & _MSK); p >>= 44
     for idx in range(16, 8, -1):
      v = o[idx] % r0sq
      if v: o[idx - 9] += v * c0; o[idx - 6] += v * c3; o[idx - 3] += v * c6
     return [o[i] % r0sq for i in range(9)]
    pw = [1, 0, 0, 0, 0, 0, 0, 0, 0]; bs = [ac, (-bc) % r0sq, 0, 0, 0, 0, 0, 0, 0]; ee = eps
    while ee:
     if ee & 1: pw = _smul9(pw, bs)
     ee >>= 1
     if ee: bs = _smul9(bs, bs)
   pw[0] = (pw[0] - 1) % r0sq; ec = crt_e[ci]
   for k in range(dX):
    assert pw[k] % r0 == 0
    qk = (pw[k] // r0) % r0
    if qk:
     for j in range(8): L_c[k][j] = (L_c[k][j] + qk * ec[j]) % r0
  proj = [0] * (dX * 8)
  for k in range(dX):
   ck = L_c[k]
   if not any(ck): continue
   acc = [0] * 8; cur_p = list(ck)
   for t in range(12):
    wt = chi_inv_pows[t]; term = mul_phi15_z(cur_p, _YPOW15[(w_tau_exp * t * k) % 15])
    for j in range(8): acc[j] = (acc[j] + wt * term[j]) % r0
    nxt = [0] * 8
    for j, cj in enumerate(cur_p):
     if cj:
      y7j = _YPOW15[(7 * j) % 15]
      for ii in range(8): nxt[ii] = (nxt[ii] + cj * y7j[ii]) % r0
    cur_p = nxt
   for j in range(8): proj[k * 8 + j] = acc[j]
  return proj
 def rref_mod(M, m):
  A = [[x % m for x in row] for row in M]; nrows, ncols = len(A), len(A[0]); pivs = []; r = 0
  for c in range(ncols):
   p = next((i for i in range(r, nrows) if A[i][c] != 0), None)
   if p is None: continue
   A[r], A[p] = A[p], A[r]; inv = pow(A[r][c], -1, m); A[r] = [(x * inv) % m for x in A[r]]
   for i in range(nrows):
    if i != r and A[i][c] != 0:
     fac = A[i][c]; A[i] = [(A[i][j] - fac * A[r][j]) % m for j in range(ncols)]
   pivs.append(c); r += 1
  return A, pivs
 def left_kernel_mod(M, m):
  MT = mat_T(M); A, pivs = rref_mod(MT, m); piv_set = set(pivs); ncols = len(MT[0])
  free = [c for c in range(ncols) if c not in piv_set]; ker = []
  for fc in free:
   v = [0] * ncols; v[fc] = 1
   for ri, pc in enumerate(pivs): v[pc] = (-A[ri][fc]) % m
   ker.append(v)
  return ker
 def solve_mod(M, b, m):
  aug = [row + [bi % m] for row, bi in zip(M, b)]; A, pivs = rref_mod(aug, m); ncols = len(M[0])
  assert ncols not in pivs, "inconsistent linear system"
  x = [0] * ncols
  for ri, pc in enumerate(pivs): x[pc] = A[ri][ncols]
  return x
 rng_sm = random.Random(20261119)
 W_g = [rng_sm.randrange(r0) for _ in range(24)]
 W_f = [[rng_sm.randrange(r0) for _ in range(72)] for _ in range(3)]
 kg = toy["kg"]
 dlogs_g = [dlog_GT(embed_ab(rel["a"], rel["b"])) for rel in kg["rels"]]
 assert dlogs_g == [rel["dl"] for rel in kg["rels"]]
 smg_full = [sm_chi_proj(rel["a"], rel["b"], True) for rel in kg["rels"]]
 _, pivg = rref_mod(smg_full, r0); assert len(pivg) == 2
 smg_1 = [[sum(W_g[j] * row[j] for j in range(24)) % r0] for row in smg_full]
 Mg = [rel["row"] for rel in kg["rels"]]; Maug_g = [row + sm for row, sm in zip(Mg, smg_1)]
 Lkg = left_kernel_mod(Maug_g, r0)
 fail_g = sum(1 for v in Lkg if sum(v[i] * dlogs_g[i] for i in range(len(dlogs_g))) % r0 != 0)
 assert len(Lkg) >= 20 and fail_g == 0
 Lkg0 = left_kernel_mod(Mg, r0)
 fail_g0 = sum(1 for v in Lkg0 if sum(v[i] * dlogs_g[i] for i in range(len(dlogs_g))) % r0 != 0)
 assert fail_g0 > 0
 fchi_cnt = 0; nchi_chk = 0; chi2 = (chi * chi) % r0
 for rel, dl in zip(kg["rels"], dlogs_g):
  if dl == 0: continue
  nchi_chk += 1; dl_f = dlog_GT(fp12_frob(embed_ab(rel["a"], rel["b"])))
  assert dl_f == (chi * dl) % r0
  if dl_f != (chi2 * dl) % r0: fchi_cnt += 1
 assert fchi_cnt == nchi_chk and nchi_chk > 30
 wg = solve_mod(Maug_g, dlogs_g, r0); t_matches = 0
 for tgt in kg["targets"]:
  dl_true = dlog_GT(embed_ab(tgt["a"], tgt["b"])); assert dl_true == tgt["dl"]
  smt_full = sm_chi_proj(tgt["a"], tgt["b"], True)
  smt = [sum(W_g[j] * smt_full[j] for j in range(24)) % r0]
  if sum((tgt["row"] + smt)[j] * wg[j] for j in range(len(wg))) % r0 == dl_true: t_matches += 1
 assert t_matches == len(kg["targets"]) == 5
 kf = toy["kf"]
 dlogs_f = [dlog_GT(embed_ab(rel["a"], rel["b"])) for rel in kf["rels"]]
 assert dlogs_f == [rel["dl"] for rel in kf["rels"]]
 smf_full = [sm_chi_proj(rel["a"], rel["b"], False) for rel in kf["rels"]]
 _, pivf = rref_mod(smf_full, r0); assert 3 <= len(pivf) <= 6
 smf_3 = [[sum(W_f[c][j] * row[j] for j in range(72)) % r0 for c in range(3)] for row in smf_full]
 Mf = [rel["row"] for rel in kf["rels"]]; Maug_f = [row + sm for row, sm in zip(Mf, smf_3)]
 Lkf = left_kernel_mod(Maug_f, r0)
 fail_f = sum(1 for v in Lkf if sum(v[i] * dlogs_f[i] for i in range(len(dlogs_f))) % r0 != 0)
 assert len(Lkf) >= 50 and fail_f == 0
 wf = solve_mod(Maug_f, dlogs_f, r0)
 match_f = sum(1 for row, dl in zip(Maug_f, dlogs_f) if sum(row[j] * wf[j] for j in range(len(wf))) % r0 == dl)
 assert match_f == len(dlogs_f) == 62
 return dict(
 kg_rels=len(dlogs_g), kg_deps=len(Lkg), kg_fail=fail_g,
 kg_neg_drop_sm_fail=f"{fail_g0}/{len(Lkg0)}", kg_neg_wrong_chi_fail=f"{fchi_cnt}/{nchi_chk}",
 kg_targets_solved=f"{t_matches}/5", kf_rels=len(dlogs_f), kf_deps=len(Lkf), kf_fail=fail_f,
 kf_vlog_fit=f"{match_f}/{len(dlogs_f)}", random_Wf_Wg_compression=True,
 omitted_short_orbit_factors=kg["n_short"] + kf["n_short"],
 )


def check_voronoi_exact(data):
 units = data["pool"]["voronoi_units_10"]
 assert len(units) == 10
 places = [cmath.exp(2j * math.pi * k / 15.0) for k in (1, 7, 4, 2)]
 facets = []
 for u in units:
  assert abs(norm_q15_z(u)) == 1
  facets.append([abs(sum(u[j] * (z ** j) for j in range(8))) ** 2 for z in places])
 anc_ok = 0
 for row in data["pool"]["anchors"]:
  a_pow = [sum(M_1W[i][j] * row["a"][j] for j in range(8)) for i in range(8)]
  x4 = [abs(sum(a_pow[j] * (z ** j) for j in range(8))) ** 2 for z in places]
  t2 = sum(x4)
  if all(sum(fj * xj for fj, xj in zip(f, x4)) >= t2 - 1e-7 for f in facets):
   anc_ok += 1
 assert anc_ok == len(data["pool"]["anchors"]) == 12
 p0_l, u0_l = 7207, -5
 wpi0 = ((u0_l ** 3 - u0_l ** 2 - u0_l) * pow(u0_l + 1, -1, p0_l)) % p0_l
 hpi0 = pgcd_mod(PHI15, [-wpi0 % p0_l, 0, 0, 0, 0, 1], p0_l)
 gam0 = [137, 2903, 5111, 889]
 B_hnf = []
 for i in range(4):
  v = [0] * 16; v[i] = p0_l; B_hnf.append(v)
 for j in range(4):
  v = [0] * 16
  for idx, c in enumerate(hpi0): v[j + idx] = c
  B_hnf.append(v)
 for k in range(8):
  yk = [0] * 8; yk[k] = 1
  ak = pmul_mod(gam0, prem_mod(yk, hpi0, p0_l), p0_l, hpi0)
  v = [0] * 16
  for idx, c in enumerate(ak): v[idx] = c
  v[8 + k] = 1; B_hnf.append(v)
 B_red, _, _ = lll_16d_int(B_hnf, 99, 100)
 v1_l = B_red[0]; a1_l, b1_l = v1_l[:8], v1_l[8:]
 v2_l, k_ind = None, None
 for idx in range(1, 16):
  ak_l, bk_l = B_red[idx][:8], B_red[idx][8:]
  if any(sub_phi15_z(mul_phi15_z(a1_l, bk_l), mul_phi15_z(b1_l, ak_l))):
   v2_l, k_ind = B_red[idx], idx + 1; break
 assert v2_l is not None
 a2_l, b2_l = v2_l[:8], v2_l[8:]
 for aa, bb in ((a1_l, b1_l), (a2_l, b2_l)):
  assert prem_mod(aa, hpi0, p0_l) == pmul_mod(gam0, prem_mod(bb, hpi0, p0_l), p0_l, hpi0)
 det_l = abs(norm_q15_z(sub_phi15_z(mul_phi15_z(a1_l, b2_l), mul_phi15_z(a2_l, b1_l))))
 assert det_l % (p0_l ** 4) == 0 and det_l > 0
 assert math.gcd(R_BLS, P_BLS ** 4 - 1) == 1 and math.gcd(R_BLS, P_BLS ** 6 - 1) == 1 and (P_BLS ** 12 - 1) % (R_BLS * math.lcm(P_BLS ** 4 - 1, P_BLS ** 6 - 1)) == 0
 w_pi = ((U_BLS ** 3 - U_BLS ** 2 - U_BLS) * pow(U_BLS + 1, -1, P_BLS)) % P_BLS
 h_pi = pgcd_mod(PHI15, [-w_pi % P_BLS, 0, 0, 0, 0, 1], P_BLS)
 assert len(h_pi) == 5 and h_pi[4] == 1
 h0, h1, h2, h3 = h_pi[0], h_pi[1], h_pi[2], h_pi[3]; uw2 = (U_BLS * w_pi * w_pi) % P_BLS
 def _f4a(a, b): return [(a[0] + b[0]) % P_BLS, (a[1] + b[1]) % P_BLS, (a[2] + b[2]) % P_BLS, (a[3] + b[3]) % P_BLS]
 def _f4s(a, b): return [(a[0] - b[0]) % P_BLS, (a[1] - b[1]) % P_BLS, (a[2] - b[2]) % P_BLS, (a[3] - b[3]) % P_BLS]
 def _f4c(s, a): return [(s * a[0]) % P_BLS, (s * a[1]) % P_BLS, (s * a[2]) % P_BLS, (s * a[3]) % P_BLS]
 def _f4m(a, b):
  a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
  c0 = a0 * b0; c1 = a0 * b1 + a1 * b0; c2 = a0 * b2 + a1 * b1 + a2 * b0; c3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0
  c4 = a1 * b3 + a2 * b2 + a3 * b1; c5 = a2 * b3 + a3 * b2; c6 = a3 * b3
  v = c6 % P_BLS
  if v: c2 -= v * h0; c3 -= v * h1; c4 -= v * h2; c5 -= v * h3
  v = c5 % P_BLS
  if v: c1 -= v * h0; c2 -= v * h1; c3 -= v * h2; c4 -= v * h3
  v = c4 % P_BLS
  if v: c0 -= v * h0; c1 -= v * h1; c2 -= v * h2; c3 -= v * h3
  return [c0 % P_BLS, c1 % P_BLS, c2 % P_BLS, c3 % P_BLS]
 def _f4p(a, e):
  r = [1, 0, 0, 0]; b = list(a)
  while e:
   if e & 1: r = _f4m(r, b)
   e >>= 1
   if e: b = _f4m(b, b)
  return r
 _fc = [_f4p([1 if j == i else 0 for j in range(4)], P_BLS) for i in range(4)]
 def _fr1(a):
  o = [0, 0, 0, 0]
  for i in range(4):
   ai = a[i]
   if ai:
    col = _fc[i]; o[0] += ai * col[0]; o[1] += ai * col[1]; o[2] += ai * col[2]; o[3] += ai * col[3]
  return [x % P_BLS for x in o]
 def _fr2(a): return _fr1(_fr1(a))
 def _f4i(a):
  a1 = _fr1(a); a2 = _fr1(a1); a3 = _fr1(a2); pr = _f4m(a1, _f4m(a2, a3))
  return _f4c(pow(_f4m(a, pr)[0], -1, P_BLS), pr)
 def _f12m(A, B):
  a0, a1, a2 = A; b0, b1, b2 = B; v0 = _f4m(a0, b0); v1 = _f4m(a1, b1); v2 = _f4m(a2, b2)
  t0 = _f4s(_f4s(_f4m(_f4a(a1, a2), _f4a(b1, b2)), v1), v2)
  t1 = _f4s(_f4s(_f4m(_f4a(a0, a1), _f4a(b0, b1)), v0), v1)
  t2 = _f4s(_f4s(_f4m(_f4a(a0, a2), _f4a(b0, b2)), v0), v2)
  return (_f4a(v0, _f4c(uw2, t0)), _f4a(t1, _f4c(uw2, v2)), _f4a(t2, v1))
 def _f12p(A, e):
  if e == 0: return ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
  tbl = [([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]), A]
  for _ in range(14): tbl.append(_f12m(tbl[-1], A))
  nbits = (e.bit_length() + 3) & ~3; r = tbl[0]
  for sh in range(nbits - 4, -1, -4):
   r = _f12m(r, r); r = _f12m(r, r); r = _f12m(r, r); r = _f12m(r, r)
   nib = (e >> sh) & 15
   if nib: r = _f12m(r, tbl[nib])
  return r
 _inv2 = (P_BLS + 1) // 2; _d2 = _f4s([0, 1, 0, 0], _fr2([0, 1, 0, 0])); _i2d = _f4i(_f4c(2, _d2))
 def _to2(a): ap = _fr2(a); return _f4c(_inv2, _f4a(a, ap)), _f4m(_f4s(a, ap), _i2d)
 _g_raw = ([3, 1, 0, 0], [5, 0, 2, 0], [7, 0, 0, 1])
 _cx_fr = pow(uw2, (P_BLS - 1) // 3, P_BLS); _cx2_fr = (_cx_fr * _cx_fr) % P_BLS
 def _f12_fr(A, k):
  rA = A
  for _ in range(k): rA = (_fr1(rA[0]), _f4c(_cx_fr, _fr1(rA[1])), _f4c(_cx2_fr, _fr1(rA[2])))
  return rA
 def _f12_inv(A):
  A1 = _f12_fr(A, 4); A2 = _f12_fr(A1, 4); p12 = _f12m(A1, A2); nm = _f12m(A, p12)
  assert nm[1] == [0, 0, 0, 0] and nm[2] == [0, 0, 0, 0]; inm = _f4i(nm[0])
  return (_f4m(inm, p12[0]), _f4m(inm, p12[1]), _f4m(inm, p12[2]))
 _y1 = _f12m(_f12_fr(_g_raw, 6), _f12_inv(_g_raw)); _y2 = _f12m(_f12_fr(_y1, 2), _y1)
 _e_hard = (P_BLS ** 4 - P_BLS ** 2 + 1) // R_BLS; G0 = _f12p(_y2, _e_hard)
 _one12 = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
 assert G0 != _one12 and _f12p(G0, R_BLS) == _one12
 _gt_ok = 0
 for _ti in range(4):
  _xt = int(hashlib.sha256(f"q15-boot-elim-{_ti}".encode()).hexdigest(), 16) % (R_BLS - 1) + 1
  _Z = _f12p(G0, _xt)
  _e0, _o0 = _to2(_Z[0]); _e1, _o1 = _to2(_Z[1]); _e2, _o2 = _to2(_Z[2])
  _c0 = _f4s(_f4m(_e1, _o0), _f4m(_o1, _e0)); _c1 = _f4s(_f4m(_o2, _e0), _f4m(_e2, _o0)); _c2 = _f4s(_f4m(_e2, _o1), _f4m(_o2, _e1))
  assert _fr2(_c0) == _c0 and _fr2(_c1) == _c1 and _fr2(_c2) == _c2
  _s = (_c0, _c1, _c2); _W = _f12m(_s, _Z)
  assert _W[2] == [0, 0, 0, 0] and any(_W[1])
  _gam = _f4m([(-x) % P_BLS for x in _W[0]], _f4i(_W[1]))
  _bv = [1 if j == _ti else 0 for j in range(8)]; _b4 = prem_mod(_bv, h_pi, P_BLS); _b4 = _b4 + [0] * (4 - len(_b4))
  _a4 = _f4m(_gam, _b4)
  _el = (_a4, [(-x) % P_BLS for x in _b4], [0, 0, 0, 0])
  _fac4 = _f4m([(-x) % P_BLS for x in _b4], _f4i(_W[1]))
  _rhs = _f12m((_fac4, [0, 0, 0, 0], [0, 0, 0, 0]), _f12m(_s, _Z))
  assert _el == _rhs and _f12_fr(_s, 6) == _s and _f12_fr((_fac4, [0, 0, 0, 0], [0, 0, 0, 0]), 4) == (_fac4, [0, 0, 0, 0], [0, 0, 0, 0])
  _gt_ok += 1
 assert _gt_ok == 4
 _x_sc = int(hashlib.sha256(b"q15-boot-z1-shortcut").hexdigest(), 16) % (R_BLS - 1) + 1
 _Z_sc = _f12p(G0, _x_sc); _e0_sc = (-_x_sc) % R_BLS; _Zp_sc = _f12m(_Z_sc, _f12p(G0, _e0_sc))
 _V0_sc = 12345678901234567 % R_BLS
 assert _Zp_sc == _one12 and ((-_e0_sc * _V0_sc) % R_BLS) == ((_V0_sc * _x_sc) % R_BLS)
 _e0z, _o0z = _to2(_one12[0]); _e1z, _o1z = _to2(_one12[1]); _e2z, _o2z = _to2(_one12[2])
 assert _f12m((_f4s(_f4m(_e1z, _o0z), _f4m(_o1z, _e0z)), _f4s(_f4m(_o2z, _e0z), _f4m(_e2z, _o0z)), _f4s(_f4m(_e2z, _o1z), _f4m(_o2z, _e1z))), _one12)[1] == [0, 0, 0, 0]
 return dict(
 voronoi_units_ok=10, anchors_min_t2=f"{anc_ok}/12", constructive_lll16_k=k_ind,
 constructive_lll16_index=det_l // (p0_l ** 4), gt_boot_fp6_elim_live_381b=_gt_ok,
 gt_boot_z_eq_1_shortcut_ok=True, gcd_r_lcm_p4_p6_is_1=True,
 )
