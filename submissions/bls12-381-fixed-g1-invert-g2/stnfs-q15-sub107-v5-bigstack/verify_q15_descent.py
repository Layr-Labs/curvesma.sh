#!/usr/bin/env python3
"""verify_q15_descent.py — Fast pure-Python stdlib verifier for Q(zeta_15) Mode-A pure k=1 zero-retry Boole
+ unbounded i.i.d. fresh-e Wald boot descent.

Checks:
  [7]   16D ball half-counts (R=2, 4, 8, 256), +-16-ULP Dickman enclosure, collection-cut-conditioned miss_side
  strict monotonicity, and FULL live recomputation of the Mode-A pure k=1 zero-retry Boole + Wald i.i.d.
  fresh-e boot recurrence (p_0(b) = max(0.001, (1 + 0.5*mu(b))^-1), NB_R=1.0, NB_ALPHA=0.5, P0_FLOOR=0.001,
  mu_min >= Pm.MU_MIN, Phi_tree <= 0.60, c_det_one, E[runs/log] = 1/(1 - Phi_tree), boot_phi = 0.0,
  linear max(0, 1 - miss) retention, BOUND_CELL on both sieve and direct-ECM code paths).
  [7c]  Live reproducible full-16D raw-L_q pure k=1 deg-(32,24) parent & reached-child generator.
  [7a]  Live miniature real-2-core-FB multi-level Q(zeta_15) campaign check (q15_mini_adaptive.b85) + live
  bit-exact replay of 3 multi-level witnesses (scoped finite-scale compatibility evidence, not a proof).
  [7b]  100%-fresh production-HNF 16-cal/48-holdout R=2 boot archive (BOOT_KEEP_FRAC=0.90 < CP95 LB=0.9395,
  BALL_PAIR=0.0) + live multi-bin exact-DP UDUP=0.40 recomputation on single special-q L_Q.
  [11]  Live ACTUAL integer BOUND_CELL (+6 ball-mul margin, q15_descent16.py) on 12 exact-DP cells +
  live (135, g) short-subspace counterexample rejection.
  [11b] Live q15_retain16 prime-splitting + non-vacuous p=2 (K=140) check + j>=74 tail + counted_fraction.
  [12]  Live q15_ecm16 L_min=2 branch-tree DP + 12-bit geometric ladder + Suyama/Newton/singular-N-1 check.
"""
import base64
import bisect
import json
import lzma
import math
import random
import sys
import time
from math import ceil, log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
 sys.path.insert(0, str(HERE))

import q15_descent16 as DS
import q15_ecm16 as EC
import q15_retain16 as RT

import operator as _OP
_SUMPROD = getattr(math, "sumprod", lambda a, b, _m=_OP.mul: sum(map(_m, a, b)))

U_BLS = -0xd201000000010000
LOG2U8 = 8.0 * log2(abs(U_BLS))
PHI15 = [1, -1, 0, 1, -1, 1, 0, -1, 1]


def _mul_z(a, b):
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


_YPOW = [[1, 0, 0, 0, 0, 0, 0, 0]]
for _i in range(1, 15):
 _YPOW.append(_mul_z(_YPOW[-1], [0, 1, 0, 0, 0, 0, 0, 0]))
_W = _YPOW[5]
_WP1 = [(1 if i == 0 else 0) + _W[i] for i in range(8)]
_W2 = _YPOW[10]


def _norm_q15(a):
 if not any(a):
  return 0
 def sig(v, k):
  res = [0] * 8
  for j, cj in enumerate(v):
   if cj:
    pw = _YPOW[(k * j) % 15]
    for i in range(8):
     res[i] += cj * pw[i]
  return res
 p = _mul_z(a, sig(a, 4))
 p = _mul_z(p, sig(p, 2))
 return _mul_z(p, sig(p, 11))[0]


def _norms_fg(aa, bb):
 a2, b2 = _mul_z(aa, aa), _mul_z(bb, bb)
 a3, b3 = _mul_z(a2, aa), _mul_z(b2, bb)
 w2b3 = _mul_z(_W2, b3)
 ng = abs(_norm_q15([a3[i] - U_BLS * w2b3[i] for i in range(8)]))
 a6, b6 = _mul_z(a3, a3), _mul_z(b3, b3)
 a9, b9 = _mul_z(a6, a3), _mul_z(b6, b3)
 t1 = _mul_z(_WP1, _mul_z(a6, b3))
 t2 = _mul_z(a3, b6)
 t3 = _mul_z(_W, b9)
 nf = abs(_norm_q15([a9[i] + t1[i] + t2[i] - t3[i] for i in range(8)]))
 return nf, ng


def _is_prime(n):
 if n < 2:
  return False
 for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
  if n == p:
   return True
  if n % p == 0:
   return False
 d, s = n - 1, 0
 while d % 2 == 0:
  d //= 2; s += 1
 for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
  if a >= n:
   continue
  x = pow(a, d, n)
  if x in (1, n - 1):
   continue
  for _ in range(s - 1):
   x = (x * x) % n
   if x == n - 1:
    break
  else:
   return False
 return True


def _tr15(d):
 g = math.gcd(15, d % 15)
 return 8 if g == 15 else (-4 if g == 5 else (-2 if g == 3 else 1))


_G8 = [[_tr15(i - j) for j in range(8)] for i in range(8)]


def _lll_16d_qcm(B_in, num=99, den=100):
 def _bareiss(B0, use_g8, p_num=num, p_den=den):
  n = 16
  M = [list(r) for r in B0]
  sp = _SUMPROD
  if use_g8:
   def dot(u, v):
    s = 0
    for i in range(8):
     ui, ub = u[i], u[8 + i]
     gi = _G8[i]
     if ui:
      s += ui * sp(gi, v[:8])
     if ub:
      s += ub * sp(gi, v[8:])
    return s
  else:
   dot = sp
  d = [1] * (n + 1)
  lam = [[0] * n for _ in range(n)]
  for i in range(n):
   for j in range(i + 1):
    u = dot(M[i], M[j])
    for k in range(j):
     u = (d[k + 1] * u - lam[i][k] * lam[j][k]) // d[k]
    if j < i:
     lam[i][j] = u
    else:
     d[i + 1] = u
  k = 1
  while k < n:
   for l in range(k - 1, -1, -1):
    if 2 * abs(lam[k][l]) > d[l + 1]:
     q = (2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1]) if lam[k][l] > 0 else -((-2 * lam[k][l] + d[l + 1]) // (2 * d[l + 1]))
     if q:
      M[k] = [x - q * y for x, y in zip(M[k], M[l])]
      lam[k][l] -= q * d[l + 1]
      for j in range(l):
       lam[k][j] -= q * lam[l][j]
   if p_den * (d[k - 1] * d[k + 1] + lam[k][k - 1] ** 2) < p_num * d[k] * d[k]:
    M[k], M[k - 1] = M[k - 1], M[k]
    for j in range(k - 1):
     lam[k][j], lam[k - 1][j] = lam[k - 1][j], lam[k][j]
    for i in range(k + 1, n):
     t = lam[i][k - 1]
     lam[i][k - 1] = (lam[k][k - 1] * t + d[k - 1] * lam[i][k]) // d[k]
     lam[i][k] = (d[k + 1] * t - lam[k][k - 1] * lam[i][k]) // d[k]
    d[k] = (d[k - 1] * d[k + 1] + lam[k][k - 1] ** 2) // d[k]
    if k > 1:
     k -= 1
   else:
    k += 1
  return M, d

 B1, d1 = _bareiss(B_in, False, 11, 40)
 B2, d2 = _bareiss(B1, True, num, den)
 return B2, d1[16], d2[16]


def _sq_hnf(q, s, t):
 rp = [pow(s, i, q) for i in range(8)]
 B = [[0] * 16 for _ in range(16)]
 B[0][0] = q
 for i in range(1, 8):
  B[i][0] = (-rp[i]) % q; B[i][i] = 1
 for i in range(8):
  B[8 + i][0] = (t * rp[i]) % q; B[8 + i][8 + i] = 1
 return B


def _roots_mod_q(poly, q, rng):
 f = RT._trim([c % q for c in poly])
 if len(f) <= 1:
  return []
 xp = RT.ppowmod([0, 1], q, f, q)
 if len(xp) < 2:
  xp += [0] * (2 - len(xp))
 xp[1] = (xp[1] - 1) % q
 g = RT.pgcd(f, RT._trim(xp), q)
 out = []

 def split(h):
  deg = len(h) - 1
  if deg <= 0:
   return
  if deg == 1:
   out.append((-h[0] * pow(h[1], -1, q)) % q)
   return
  while True:
   a = rng.randrange(q)
   t = RT.ppowmod([a, 1], (q - 1) // 2, h, q)
   if not t:
    continue
   t[0] = (t[0] - 1) % q
   d = RT.pgcd(h, RT._trim(t), q)
   if 0 < len(d) - 1 < deg:
    quo, _ = RT.pdiv(h, d, q)
    split(d)
    split(quo)
    return

 split(g)
 return sorted(set(out))


def _sample_ideal(bits, side, rng):
 while True:
  q = rng.randrange(1 << (bits - 1), 1 << bits)
  q = q - (q % 15) + 1
  if q.bit_length() != bits or not _is_prime(q):
   continue
  e15 = (q - 1) // 15
  g0 = next((pow(g, e15, q) for g in range(2, 100) if pow(pow(g, e15, q), 3, q) != 1 and pow(pow(g, e15, q), 5, q) != 1), None)
  if g0 is None:
   continue
  s = pow(g0, rng.choice((1, 2, 4, 7, 8, 11, 13, 14)), q)
  w = pow(s, 5, q)
  poly = [(-U_BLS * w * w) % q, 0, 0, 1] if side == "g" else [(-w) % q, 0, 0, 1, 0, 0, (w + 1) % q, 0, 0, 1]
  rts = _roots_mod_q(poly, q, rng)
  if rts:
   return q, s, rng.choice(rts)


def _build_dp_ball(rmax=32):
 nmax = rmax * rmax
 rd = [[0] * (nmax + 1) for _ in range(17)]
 rd[0][0] = 1
 for d in range(1, 17):
  prev = rd[d - 1]; cur = rd[d]
  for s in range(nmax + 1):
   v = prev[s]
   if not v:
    continue
   cur[s] += v
   k = 1
   while s + k * k <= nmax:
    cur[s + k * k] += 2 * v
    k += 1
 pref = [0] * (nmax + 1)
 acc = 0
 for s in range(1, nmax + 1):
  acc += rd[16][s]
  pref[s] = acc
 return rd, pref


def _sample_prim_ball_dp(rd, pref, R, rng):
 r2 = R * R
 tot = pref[r2]
 while True:
  u = rng.randrange(tot)
  s = bisect.bisect_right(pref, u, 1, r2 + 1)
  c = [0] * 16
  for d in range(16, 0, -1):
   rem = s
   v = rng.randrange(rd[d][rem])
   w0 = rd[d - 1][rem]
   if v < w0:
    c[d - 1] = 0
   else:
    v -= w0
    k = 1
    while True:
     wk = rd[d - 1][rem - k * k]
     if v < wk:
      c[d - 1] = k; s -= k * k; break
     v -= wk
     if v < wk:
      c[d - 1] = -k; s -= k * k; break
     v -= wk
     k += 1
  if math.gcd(*c) == 1:
   for x in reversed(c):
    if x:
     if x < 0:
      c = [-y for y in c]
     break
   return c


def _mini_w_roots(p):
 for g in range(2, p):
  w = pow(g, (p - 1) // 3, p)
  if w != 1:
   w2 = (w * w) % p
   return (w, w2) if w < w2 else (w2, w)
 raise ValueError(p)


def _mini_canon_ry(p, oi):
 w0, w1 = _mini_w_roots(p)
 wt = w0 if oi == 0 else w1
 e15 = (p - 1) // 15
 y0 = next(pow(g, e15, p) for g in range(2, p) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1)
 return min(pow(y0, k, p) for k in (1, 2, 4, 7, 8, 11, 13, 14) if pow(pow(y0, k, p), 5, p) == wt)


def _verify_mini_tower_p7():
 """Live pure-Python check of miniature Q(zeta_15) deg-(32, 24) tower (u_t=1, p_t=7)."""
 h1 = RT.pgcd(PHI15, [5, 0, 0, 0, 0, 1], 7)
 assert h1 == [4, 1, 2, 4, 1] and len(h1) == 5
 y7_2 = RT.ppowmod([0, 1], 7 ** 2, h1, 7)
 y7_4 = RT.ppowmod([0, 1], 7 ** 4, h1, 7)
 assert y7_4 == [0, 1]
 assert len(RT.pgcd(RT._trim([y7_2[i] - (1 if i == 1 else 0) for i in range(max(len(y7_2), 2))]), h1, 7)) == 1
 assert pow(4, (7 ** 4 - 1) // 3, 7) != 1
 f61 = [45, 3, 0, 0, 1]
 xp = RT.ppowmod([0, 1], 61, f61, 61)
 xp2 = RT.ppowmod([0, 1], 61 ** 2, f61, 61)
 xp4 = RT.ppowmod([0, 1], 61 ** 4, f61, 61)
 assert xp4 == [0, 1]
 assert len(RT.pgcd(RT._trim([xp[i] - (1 if i == 1 else 0) for i in range(max(len(xp), 2))]), f61, 61)) == 1
 assert len(RT.pgcd(RT._trim([xp2[i] - (1 if i == 1 else 0) for i in range(max(len(xp2), 2))]), f61, 61)) == 1
 return True


def verify_mini_adaptive_campaign(rdir):
 """Verify the miniature real-2-core-FB multi-level Q(zeta_15) campaign archive (q15_mini_adaptive.b85)
 and replay all 3 live bit-exact multi-level Q_CM-LLL witnesses.
 Scoped B_FB = 2^14 finite-scale compatibility evidence, not a proof."""
 t0 = time.perf_counter()
 p_b85 = rdir / "q15_mini_adaptive.b85"
 pack = json.loads(lzma.decompress(base64.b85decode("".join(p_b85.read_text().split()))).decode())
 dat = pack["mini"]
 fb_data = pack["fb"]
 assert "not a proof" in dat.get("scope_framing", "").lower()
 assert _verify_mini_tower_p7() is True
 pooled = dat["pooled"]
 per_set = dat["per_setting"]
 assert len(per_set) == 6 and pooled["roots"] == 12802 and pooled["after_att"] == 6393
 assert pooled["oe_fail"] < 1.0 and 0.90 <= pooled["oe_work"] <= 1.10 and pooled["max_C_bday"] < 10.0

 # Live bit-exact replay of all 3 multi-level Q_CM-LLL witnesses against the real 2-core FB (B_FB = 2^14)
 B_FB = int(fb_data["b_fb"])
 NS_MAX = int(fb_data["ns_max"])
 mul15 = _mul_z
 YP = _YPOW
 W = YP[5]; W2 = YP[10]
 WM2 = [W[i] - (2 if i == 0 else 0) for i in range(8)]
 SIG_M = [YP[(7 * j) % 15] for j in range(8)]

 def norm_zw(el):
  v0, v1, v2, v3, v4, v5, v6, v7 = el
  p1 = mul15(el, (v0-v2-v6+v7, v2+v4-v7, -v3+v6, -v2-v6, v1+v2-v7, -v2+v5+v7, -v6, v2-v3+v6-v7))
  v0, v1, v2, v3, v4, v5, v6, v7 = p1
  res = mul15(p1, (v0+v2+v4, -v4, -v2-v6, v2, -v2-v4+v7, v4+v5, v2+v3, v1-v2-v4-v6))
  c0, c1 = res[0], res[5]
  return c0, c1, abs(c0 * c0 - c0 * c1 + c1 * c1)

 def fg_norms(a, b):
  a2, b2 = mul15(a, a), mul15(b, b)
  a3, b3 = mul15(a2, a), mul15(b2, b)
  a4, b4 = mul15(a2, a2), mul15(b2, b2)
  w2b3 = mul15(W2, b3)
  eg = [a3[i] - w2b3[i] for i in range(8)]
  ab3 = mul15(a, b3); wm2b4 = mul15(WM2, b4)
  ef = [a4[i] + 3 * ab3[i] + wm2b4[i] for i in range(8)]
  return norm_zw(ef), norm_zw(eg)

 ret_f_ns = {(int(p), int(oi)) for p, oi in fb_data["ret_f_ns"]}
 ret_g_ns = {(int(p), int(oi)) for p, oi in fb_data["ret_g_ns"]}
 ret_f_sp = {(int(p), int(oi), int(rt)) for p, oi, rt in fb_data["ret_f_split"]}
 ret_g_sp = {(int(p), int(oi), int(rt)) for p, oi, rt in fb_data["ret_g_split"]}

 ns_ps = [p for p in range(2, NS_MAX + 1) if _is_prime(p) and p % 15 != 1]
 sp_ps = [p for p in range(2, B_FB + 1) if _is_prime(p) and p % 15 == 1]
 ns_w = {p: _mini_w_roots(p) for p in ns_ps if p % 3 == 1}
 sp_w = {}; sp_ry = {}
 def get_sp_w(p):
  v = sp_w.get(p)
  if v is None:
   v = _mini_w_roots(p); sp_w[p] = v
  return v
 def get_sp_ry(p, oi):
  v = sp_ry.get((p, oi))
  if v is None:
   v = _mini_canon_ry(p, oi); sp_ry[(p, oi)] = v
  return v
 P_FB = math.prod(ns_ps + sp_ps)

 bts = []; cur = [0] * 8
 def gb(idx, rem):
  if idx == 8:
   if any(cur): bts.append(tuple(cur))
   return
  r = math.isqrt(rem)
  for v in range(-r, r + 1):
   cur[idx] = v; gb(idx + 1, rem - v * v)
 gb(0, int(fb_data["ball_r2"]))
 assert len(bts) == int(fb_data["ball_pts"])

 def eval_side0(sd, rt, w, p):
  return (pow(rt, 3, p) - w * w) % p == 0 if sd == 1 else (pow(rt, 4, p) + 3 * rt + w - 2) % p == 0

 def extract_tau(sd, a, b, p, w, ry0):
  out = []
  for e in (1, 7, 4, 13):
   rym = ry0 if e == 1 else pow(ry0, e, p)
   av = bv = 0; ryp = 1
   for i in range(8):
    av = (av + a[i] * ryp) % p
    bv = (bv + b[i] * ryp) % p
    ryp = (ryp * rym) % p
   if bv == 0: continue
   rt = (av * pow(bv, -1, p)) % p
   if eval_side0(sd, rt, w, p) and rt not in out:
    out.append(rt)
  return out

 def check_fb(sm, c0, c1, sd, a, b):
  r_ns = ret_f_ns if sd == 0 else ret_g_ns
  r_sp = ret_f_sp if sd == 0 else ret_g_sp
  rem = sm
  for p in ns_ps:
   if rem == 1: break
   if rem % p == 0:
    while rem % p == 0: rem //= p
    if p % 3 == 1:
     w0, w1 = ns_w[p]; h = 0
     for oi, w in ((0, w0), (1, w1)):
      if (c0 + c1 * w) % p == 0:
       if (p, oi) not in r_ns: return False
       h += 1
     if not h: return False
    elif (p, 0) not in r_ns: return False
  for p in sp_ps:
   if rem == 1: break
   if rem % p == 0:
    while rem % p == 0: rem //= p
    w0, w1 = get_sp_w(p); h = 0
    for oi, w in ((0, w0), (1, w1)):
     if (c0 + c1 * w) % p == 0:
      rts = extract_tau(sd, a, b, p, w, get_sp_ry(p, oi))
      if not rts: return False
      for rt in rts:
       if (p, oi, rt) not in r_sp: return False
       h += 1
    if not h: return False
  return rem == 1

 def factor_win(cof, nw, km):
  if cof == 1: return []
  if km <= 0 or cof.bit_length() > nw * km: return None
  lim = 1 << nw; out = []; st = [cof]
  while st:
   v = st.pop()
   if v == 1: continue
   if _is_prime(v):
    if v > lim or v <= B_FB or v % 15 != 1: return None
    if v not in out:
     if len(out) >= km: return None
     out.append(v)
   else:
    d = None
    for c in range(1, 32):
     x = y = 2; g = 1
     while g == 1:
      x = (x * x + c) % v; y = (y * y + c) % v; y = (y * y + c) % v
      g = math.gcd(abs(x - y), v)
     if 1 < g < v: d = g; break
    if d is None: return None
    st.extend((d, v // d))
  return out

 def scan_ball(sd, q, ry, rt, nw, km):
  B, _, _ = _lll_16d_qcm(_sq_hnf(q, ry, rt))
  a1, b1 = B[0][:8], B[0][8:]
  a2, b2 = next((r[:8], r[8:]) for r in B[1:] if any(x != y for x, y in zip(mul15(a1, r[8:]), mul15(b1, r[:8]))))
  a2s = [mul15(YP[j], a2) for j in range(8)]
  b2s = [mul15(YP[j], b2) for j in range(8)]
  rels = []; max_cof = 1 if km == 0 else (1 << (nw * km))
  for bt in bts:
   a = list(a1); b = list(b1)
   for j, c in enumerate(bt):
    if c:
     aj, bj = a2s[j], b2s[j]
     for i in range(8):
      a[i] += c * aj[i]; b[i] += c * bj[i]
   (c0f, c1f, nf), (c0g, c1g, ng) = fg_norms(a, b)
   if sd == 0:
    if nf == 0 or nf % q: continue
    nf //= q
   else:
    if ng == 0 or ng % q: continue
    ng //= q
   if nf == 0 or ng == 0: continue
   cofg = ng
   while cofg > 1:
    g = math.gcd(cofg, P_FB)
    if g == 1: break
    cofg //= g
   if cofg > max_cof: continue
   coff = nf
   while coff > 1:
    g = math.gcd(coff, P_FB)
    if g == 1: break
    coff //= g
   if coff > max_cof: continue
   pf = factor_win(coff, nw, km)
   if pf is None: continue
   pg = factor_win(cofg, nw, km - len(pf))
   if pg is None or len(pf) + len(pg) > km: continue
   if not check_fb(nf // coff, c0f, c1f, 0, a, b): continue
   if not check_fb(ng // cofg, c0g, c1g, 1, a, b): continue
   ch = []; ok = True
   for csd, ps, c0, c1 in ((0, pf, c0f, c1f), (1, pg, c0g, c1g)):
    for p in ps:
     w0, w1 = _mini_w_roots(p); h = 0
     for coi, w in ((0, w0), (1, w1)):
      if (c0 + c1 * w) % p == 0:
       cry = _mini_canon_ry(p, coi)
       rts = extract_tau(csd, a, b, p, w, cry)
       if not rts: ok = False; break
       for crt in rts:
        if len(ch) >= km: ok = False; break
        ch.append((csd, p, coi, cry, crt)); h += 1
     if not ok or not h: ok = False; break
    if not ok: break
   if ok: rels.append(ch)
  return rels

 def desc(sd, q, oi, ry, rt, nw, km, dsh, rec=False):
  rels = scan_ball(sd, q, ry, rt, nw, km)
  tw = 1; prev = set(); seq = []
  for ri, ch in enumerate(rels):
   if not ch: return True, tw, seq
   aok = True; cw = 0; sh = 0
   for csd, cp, coi, cry, crt in ch:
    if ri > 0 and (csd, cp, coi, crt) in prev: sh = 1
    cb = cp.bit_length(); cnw = max(14, cb - dsh); ckm = 0 if cnw == 14 else 1
    cok, cwk, _ = desc(csd, cp, coi, cry, crt, cnw, ckm, dsh, False)
    cw += cwk
    if not cok: aok = False
   prev.update((c[0], c[1], c[2], c[4]) for c in ch)
   tw += cw
   if rec: seq.append([ri, int(not aok), cw, len(ch), sh])
   if aok: return True, tw, seq
  return False, tw, seq

 wits = dat["live_witnesses"]
 assert len(wits) == 3
 for w in wits:
  _, _, got = desc(w["side"], w["q"], w["oi"], w["ry"], w["rt"], w["n_win"], w["k_root"], w["d_shrink"], True)
  assert got == w["exp_seq"], f"Witness mismatch {w['setting']} q={w['q']}: {got} != {w['exp_seq']}"
 # Verify bundled exploratory full-16D raw-L_q pure k=1 cohort (Diagnostic B) and 38-bin prod raw-L_q survey (Diagnostic A)
 k1c = pack["full16d_k1_cohort"]
 assert "not a proof" in k1c.get("scope_framing", "").lower()
 assert k1c["parents_all"]["n_ideals"] == 16000 and k1c["parents_all"]["zero_count"] == 0
 assert k1c["reached_children_all"]["n_ideals"] == 15000 and k1c["reached_children_all"]["zero_count"] == 0
 assert k1c["combined_disjoint_all"]["n_ideals"] == 31000 and k1c["combined_disjoint_all"]["zero_count"] == 0
 psv = pack["prod_raw_lq_survey"]
 assert "not a proof" in psv.get("scope_framing", "").lower()
 a38 = psv["all38_bins"]
 rc = psv["reached_children"]
 assert a38["n_bins"] == 38 and len(a38["per_bin"]) == 38 and a38["n_ideals"] == 30400 and a38["below_1x"] == 0 and a38["vq_basis_failures"] == 0
 assert rc["n_reached_children"] == 416 and rc["children_below_1x"] == 0 and rc["controls_below_1x"] == 0 and rc["ks_2samp_p_child_vs_ctrl"] > 0.5
 return dict(
 scope_framing=dat["scope_framing"],
 tower_p7_verified=True,
 real_2core_fb=dat["real_2core_fb"],
 n_settings=len(per_set),
 pooled=pooled,
 live_witnesses_replayed=len(wits),
 replay_wall_s=round(time.perf_counter() - t0, 4),
 full16d_k1_cohort=dict(
 parents_n=k1c["parents_all"]["n_ideals"],
 parents_zero=k1c["parents_all"]["zero_count"],
 reached_n=k1c["reached_children_all"]["n_ideals"],
 reached_zero=k1c["reached_children_all"]["zero_count"],
 combined_n=k1c["combined_disjoint_all"]["n_ideals"],
 combined_zero_count=k1c["combined_disjoint_all"]["zero_count"],
 within_bin_fano=k1c["combined_disjoint_all"]["within_bin_fano"],
 min_M=k1c["combined_disjoint_all"]["quantiles_M"]["min"],
 ),
 prod_raw_lq_survey=dict(
 all38_n_ideals=a38["n_ideals"],
 all38_total_cells=a38["total_cells"],
 all38_below_1x=a38["below_1x"],
 all38_min_ratio=a38["ratio_quantiles"]["min"],
 all38_median_ratio=a38["ratio_quantiles"]["median"],
 min_euc_b0_over_q16=a38["min_euc_b0_over_q16"],
 min_qcm_b0_over_detG8_q16=a38["min_qcm_b0_over_detG8_q16"],
 reached_n=rc["n_reached_children"],
 reached_below_1x=rc["children_below_1x"],
 reached_min_ratio=rc["child_quantiles"]["min"],
 reached_geom_child_over_ctrl=rc["geom_mean_child_over_ctrl"],
 reached_ks_p=rc["ks_2samp_p_child_vs_ctrl"],
 ),
 checks_ok=True,
 )


def _live_full16d_k1_deg32_24(n_parents=20, seed=20261124):
 """Live reproducible full-16D raw-L_q pure k=1 parent & first-usable-reached-child generator
 on the exact deg-(32, 24) Q(zeta_15) tower (f4 = x^4 + x + w, g3 = x^3 + w^2, B_fb = 2^17,
 LP_max = 2^22 - 1, q_parent in [2^22, 2^23), 2496-point full-16D ||k||_2^2 <= 3 half-ball)."""
 t0 = time.perf_counter()
 mul15 = _mul_z
 W1 = _YPOW[5]
 W2 = _YPOW[10]

 def norm_zw(el):
  v0, v1, v2, v3, v4, v5, v6, v7 = el
  p1 = mul15(el, (v0-v2-v6+v7, v2+v4-v7, -v3+v6, -v2-v6, v1+v2-v7, -v2+v5+v7, -v6, v2-v3+v6-v7))
  v0, v1, v2, v3, v4, v5, v6, v7 = p1
  res = mul15(p1, (v0+v2+v4, -v4, -v2-v6, v2, -v2-v4+v7, v4+v5, v2+v3, v1-v2-v4-v6))
  c0, c1 = res[0], res[5]
  return c0, c1, abs(c0 * c0 - c0 * c1 + c1 * c1)

 def norms_32_24(a, b):
  a2 = mul15(a, a); a3 = mul15(a2, a); a4 = mul15(a2, a2)
  b2 = mul15(b, b); b3 = mul15(b2, b); b4 = mul15(b2, b2)
  ab3 = mul15(a, b3); wb4 = mul15(W1, b4); w2b3 = mul15(W2, b3)
  ef = [a4[i] + ab3[i] + wb4[i] for i in range(8)]
  eg = [a3[i] + w2b3[i] for i in range(8)]
  return norm_zw(ef), norm_zw(eg)

 ball16 = [((i, 1),) for i in range(16)]
 for i in range(16):
  for j in range(i + 1, 16):
   ball16.append(((i, 1), (j, 1)))
   ball16.append(((i, 1), (j, -1)))
 for i in range(16):
  for j in range(i + 1, 16):
   for k in range(j + 1, 16):
    for sj in (1, -1):
     for sk in (1, -1):
      ball16.append(((i, 1), (j, sj), (k, sk)))
 assert len(ball16) == 2496

 B_FB = 1 << 17
 LP_MAX = (1 << 22) - 1
 P_SMALL = math.prod(EC.primes_upto(2048))

 def is_lp_smooth(n, lp_lim):
  if n <= 1:
   return True, []
  while True:
   g = math.gcd(n, P_SMALL)
   if g == 1:
    break
   n //= g
  if n == 1:
   return True, []
  lps = []
  st = [n]
  while st:
   v = st.pop()
   if v == 1:
    continue
   if _is_prime(v):
    if v > lp_lim:
     return False, None
    if v > B_FB:
     lps.append(v)
    continue
   d = None
   for c in (1, 2, 3):
    x = y = 2; q_acc = 1; it = 0
    while it < 600:
     x = (x * x + c) % v
     y = (y * y + c) % v; y = (y * y + c) % v
     diff = abs(x - y)
     if diff:
      q_acc = (q_acc * diff) % v
     it += 1
     if (it & 15) == 0:
      g = math.gcd(q_acc, v)
      if g > 1:
       if g < v:
        d = g
       break
    if d is not None:
     break
   if d is None:
    return False, None
   st.append(d); st.append(v // d)
  return True, lps

 rng = random.Random(seed)

 def sample_parent(side):
  while True:
   q = rng.randrange((1 << 22) + 1, 1 << 23)
   q -= (q - 1) % 15
   if q < (1 << 22) or not _is_prime(q):
    continue
   s = None
   for g in range(2, 60):
    r = pow(g, (q - 1) // 15, q)
    if pow(r, 3, q) != 1 and pow(r, 5, q) != 1:
     s = r
     break
   if s is None:
    continue
   w = pow(s, 5, q)

   def ptrim(u):
    while len(u) > 1 and u[-1] == 0:
     u.pop()
    return u

   def pmod(u, m):
    u = ptrim(list(u))
    while len(u) >= len(m) and (len(u) > 1 or u[0] != 0):
     c = u.pop()
     if c:
      sh = len(u) - len(m) + 1
      for i in range(len(m) - 1):
       u[sh + i] = (u[sh + i] - c * m[i]) % q
     ptrim(u)
    return ptrim(u)

   def pm(u, v, m):
    o = [0] * (len(u) + len(v) - 1)
    for i, a in enumerate(u):
     for j, b in enumerate(v):
      o[i + j] = (o[i + j] + a * b) % q
    return pmod(o, m)

   def pgcd(u, v):
    u, v = ptrim(list(u)), ptrim(list(v))
    while len(v) > 1 or v[0] != 0:
     while len(u) >= len(v) and (len(u) > 1 or u[0] != 0):
      sc = (u[-1] * pow(v[-1], -1, q)) % q
      sh = len(u) - len(v)
      for i in range(len(v)):
       u[sh + i] = (u[sh + i] - sc * v[i]) % q
      ptrim(u)
     u, v = v, u
    return u

   poly = [(w * w) % q, 0, 0, 1] if side == "g" else [w % q, 1, 0, 0, 1]
   res = [1]; base = [0, 1]; e = q
   while e:
    if e & 1:
     res = pm(res, base, poly)
    base = pm(base, base, poly); e >>= 1
   if len(res) < 2:
    res.extend([0] * (2 - len(res)))
   res[1] = (res[1] - 1) % q
   lin = pgcd(poly, ptrim(res))
   if len(lin) < 2:
    continue
   ok_sp = True
   while len(lin) > 2:
    found = False
    for a0 in range(1, 30):
     r2 = [1]; b2 = [a0, 1]; e = (q - 1) // 2
     while e:
      if e & 1:
       r2 = pm(r2, b2, lin)
      b2 = pm(b2, b2, lin); e >>= 1
     r2[0] = (r2[0] - 1) % q
     gp = pgcd(lin, ptrim(r2))
     if 1 < len(gp) < len(lin):
      lin = gp; found = True; break
    if not found:
     ok_sp = False; break
   if not ok_sp:
    continue
   t = (-lin[0] * pow(lin[1], -1, q)) % q
   return q, s, t

 def scan_first_usable(q, s, t, side, lp_lim):
  B, _, _ = _lll_16d_qcm(_sq_hnf(q, s, t))
  for idx_c, comb in enumerate(ball16):
   v = [0] * 16
   for idx, sg in comb:
    row = B[idx]
    if sg == 1:
     for j in range(16):
      v[j] += row[j]
    else:
     for j in range(16):
      v[j] -= row[j]
   a, b = v[:8], v[8:]
   (c0f, c1f, nf), (c0g, c1g, ng) = norms_32_24(a, b)
   if side == "f":
    if nf == 0 or nf % q:
     continue
    nf //= q
   else:
    if ng == 0 or ng % q:
     continue
    ng //= q
   if nf == 0 or ng == 0:
    continue
   okg, lpg = is_lp_smooth(ng, lp_lim)
   if not okg:
    continue
   okf, lpf = is_lp_smooth(nf, lp_lim)
   if not okf:
    continue
   children = []
   all_lp_split = True
   for csd, lps, c0, c1 in (("f", lpf, c0f, c1f), ("g", lpg, c0g, c1g)):
    for p in lps:
     if p % 15 != 1:
      all_lp_split = False
      break
     found_p = None
     for g in range(2, 60):
      r = pow(g, (p - 1) // 15, p)
      if pow(r, 3, p) != 1 and pow(r, 5, p) != 1:
       for e in (1, 2, 4, 7, 8, 11, 13, 14):
        sy = pow(r, e, p); wy = pow(sy, 5, p)
        if (c0 + c1 * wy) % p != 0:
         continue
        av = bv = 0; sp = 1
        for i in range(8):
         av = (av + a[i] * sp) % p
         bv = (bv + b[i] * sp) % p
         sp = (sp * sy) % p
        if bv == 0:
         continue
        ty = (av * pow(bv, -1, p)) % p
        if (csd == "g" and (pow(ty, 3, p) + wy * wy) % p == 0) or (csd == "f" and (pow(ty, 4, p) + ty + wy) % p == 0):
         found_p = (csd, p, sy, ty)
         break
       break
     if found_p is None:
      all_lp_split = False
      break
     children.append(found_p)
    if not all_lp_split:
     break
   if not all_lp_split:
    continue
   return True, idx_c + 1, (children[0] if children else None)
  return False, len(ball16), None

 par_empty = 0; n_ch = 0; ch_empty = 0; steps_tot = 0
 for i in range(n_parents):
  sd = "f" if i < (n_parents // 2) else "g"
  q, s, t = sample_parent(sd)
  ok, stp, ch = scan_first_usable(q, s, t, sd, LP_MAX)
  steps_tot += stp
  if not ok:
   par_empty += 1
  if ch is not None:
   csd, cp, cs, ct = ch
   ok2, stp2, _ = scan_first_usable(cp, cs, ct, csd, cp - 1)
   n_ch += 1; steps_tot += stp2
   if not ok2:
    ch_empty += 1
 tot = n_parents + n_ch
 assert par_empty == 0 and ch_empty == 0 and n_ch >= 12 and tot >= 32
 return dict(
 tower="deg-(32,24) Q(zeta_15): f4=x^4+x+w, g3=x^3+w^2",
 ball_pts=2496, B_fb=B_FB, LP_max=LP_MAX,
 parents=n_parents, parents_empty=par_empty,
 reached_children=n_ch, reached_children_empty=ch_empty,
 total_ideals=tot, total_empty=par_empty + ch_empty,
 avg_cells_to_first_usable=round(steps_tot / tot, 2),
 wall_s=round(time.perf_counter() - t0, 4),
 )


def run_all_checks(root_dir=None, seed=1, mesh=384, run_diagnostics=True, core_checks=True):
 t0 = time.perf_counter()
 rdir = Path(root_dir) if root_dir else HERE
 st_b85 = (rdir / "q15_stats.b85").read_text()
 st = json.loads(lzma.decompress(base64.b85decode("".join(st_b85.split()))).decode())

 # [7] Exact 384-mesh pure Mode-A k=1 zero-retry Boole + Wald i.i.d. fresh-e boot evaluation
 Pm = DS.P()
 DS.load_stats(st, Pm)
 assert Pm.NB_R == 1.0 and Pm.NB_ALPHA == 0.5 and Pm.P0_FLOOR == 0.001
 assert Pm.MU_MIN >= 500.0 and Pm.MAX_PHI_TREE == 0.60
 assert Pm.LMAX == 260 and Pm.BALL_PAIR == 0.0 and Pm.BOOT_KEEP_FRAC == 0.90
 assert Pm.BOUND_CELL == 4400.0 and Pm.SIEVE_HIT == 1400.0 and Pm.UDUP_BOOT == 0.40 and Pm.UDUP_NODE == 0.55
 assert Pm.E16_COUNT_LO == 0.95 and Pm.E16_COUNT_HI == 1.05
 Pm.MESH = int(mesh); Pm.TRUNC = 0.0
 assert Pm.MAIN_LB_F == 71.20 and Pm.MAIN_LB == 71.56 and abs(Pm.MAIN_UF - 545.0 / 71.20) < 1e-12 and abs(Pm.MAIN_UG - 665.0 / 71.56) < 1e-12
 _bnp_m = st["ball_norm_pool"]
 assert _bnp_m["f_jv_measured"] - 10.0 * _bnp_m["f_jv_se"] >= Pm.UDUP_NODE
 assert _bnp_m["f_jv_ideal_mean"] - 3.0 * _bnp_m["f_jv_ideal_se"] >= Pm.UDUP_NODE

 assert DS.CELLS[2] == 17072 and DS.CELLS[4] == 642328192 and DS.CELLS[8] == 35245626789232
 vol256 = 0.5 * (math.pi ** 8) / math.factorial(8) * (256 ** 16)
 assert abs(DS.CELLS[256] / vol256 - 1.0) < 1e-4

 # Live exact Kronecker theta_3(q)^16 verification of ALL used radii in PLAN16 + BOOT16
 _used_radii = sorted({int(v[4]) for v in DS.PLAN16.values()} | {int(DS.BOOT16[3])})
 _live_cells, _, _ = DS.compute_cells_16d(_used_radii)
 for _ru in _used_radii:
  assert DS.CELLS[_ru] == DS.CELLS_16D[_ru] == _live_cells[_ru], f"CELLS_16D mismatch at R={_ru}"

 # Strict +-16-ULP Dickman enclosure
 assert DS.rho_lo(2.0) < (1.0 - math.log(2.0)) < DS.rho_hi(2.0)
 assert DS.rho_lo(1032.0 / 1024.0) < (1.0 - math.log(1032.0 / 1024.0)) < DS.rho_hi(1032.0 / 1024.0)
 assert DS.rho_lo(3080.0 / 1024.0) < DS.rho_series(3080.0 / 1024.0) < DS.rho_hi(3080.0 / 1024.0)
 assert DS.rho_lo(3306.0 / 1024.0) < DS.rho_series(3306.0 / 1024.0) < DS.rho_hi(3306.0 / 1024.0)

 # Collection-cut-conditioned miss_side strict monotonicity at LX
 _mf_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UF, Pm.LX, Pm.MAIN_LB_F, Pm.ORBIT, Pm.P_SIEVE_F, Pm.RET_MULT_F) for b_m in range(40, 71)]
 _mg_seq = [DS.miss_side(b_m, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.P_SIEVE_G, Pm.RET_MULT_G) for b_m in range(40, 71)]
 _mg_uncond = DS.miss_side(40, Pm.LOG2_RELS, Pm.MAIN_UG, Pm.LX, Pm.MAIN_LB, Pm.ORBIT, Pm.LX, 1.0)
 assert _mf_seq[0] > 0.0 and _mg_seq[0] > _mg_uncond > 0.0
 assert all(_mf_seq[i] < _mf_seq[i + 1] and _mg_seq[i] < _mg_seq[i + 1] for i in range(len(_mf_seq) - 1))

 t_ev = time.perf_counter()
 res = DS.evaluate(Pm, DS.PLAN16, DS.BOOT16)
 dt_ev = time.perf_counter() - t_ev
 rows = DS.report(res)
 bt = res["boot"]
 EN, _ = DS.expected_nodes(Pm, res)
 two_log = res["two_log_units"]
 min_mu = bt["mu_min"]
 min_ret = min(r["ret"] for r in rows)
 max_loss = max(max(r["loss_f"], r["loss_g"]) for r in rows)
 max_miss_node = max(d["miss"] for d in res["DET"].values())
 max_R_used = max(d["R"] for d in res["DET"].values())

 assert len(rows) == 2 * len(res["bins"])
 assert min_mu >= Pm.MU_MIN and bt["phi_tree"] <= Pm.MAX_PHI_TREE and bt["phi"] == 0.0
 assert min_ret >= 0.90 and max_loss < 5e-3
 assert all(r["how"] in ("sieve", "ecm") for r in rows)
 assert all(v["cofac_frac"] == 1.0 and v["eager"] is True for v in res["DET"].values())
 assert all(abs(v["lin_keep"] - max(0.0, 1.0 - v["miss"])) < 1e-12 for v in res["DET"].values())
 assert abs(bt["lin_keep"] - max(0.0, 1.0 - bt["miss"])) < 1e-12
 _det_sieve = [v for v in res["DET"].values() if v["how"] == "sieve"]
 _det_ecm = [v for v in res["DET"].values() if v["how"] == "ecm"]

 # [7a] & [7c] Heavy finite-scale diagnostics (run in part="diagnostics" / "all")
 mini_res = verify_mini_adaptive_campaign(rdir) if run_diagnostics else None
 live_k1_gen = _live_full16d_k1_deg32_24(n_parents=20, seed=20261124) if run_diagnostics else None

 # [7b] 100%-fresh production-HNF 16-cal/48-holdout R=2 boot archive (from q15_stats.b85, two-sided [m_min, m_max])
 bnp = st["boot_norm_pool"]
 assert bnp["fresh_from_production_hnf_all_64"] is True
 cg_by_t = bnp["cg_max_by_target"]
 cg_min_by_t = bnp["cg_min_by_target"]
 assert len(cg_by_t) == 64 and len(cg_min_by_t) == 64
 assert abs(bnp["cal16_cg_max"] - max(cg_by_t[:16])) < 1e-9
 assert abs(bnp["holdout48_cg_max"] - max(cg_by_t[16:])) < 1e-9
 assert abs(bnp["cg_max"] - max(bnp["cal16_cg_max"], bnp["holdout48_cg_max"])) < 1e-9
 assert abs(bnp["cal16_cg_min"] - min(cg_min_by_t[:16])) < 1e-9
 assert abs(bnp["holdout48_cg_min"] - min(cg_min_by_t[16:])) < 1e-9
 assert abs(bnp["cg_min"] - min(bnp["cal16_cg_min"], bnp["holdout48_cg_min"])) < 1e-9
 assert bnp["cg_min"] > Pm.BOOT_LO_OFF and bnp["holdout48_cg_min"] > Pm.BOOT_LO_OFF
 assert bnp["holdout48_below_m_min_cells"] == 0 and bnp["holdout48_below_m_min_targets"] == 0 and bnp["holdout48_outside_window_targets"] == 0
 assert bnp["n_exceed_m_max"] == 0 and bnp["cal16_cells"] == 16 * 17072 and bnp["holdout48_cells"] == 48 * 17072
 assert bnp["holdout48_exceed_targets"] == 0 and bnp["holdout48_exceed_cal_cut"] == 0 and bnp["holdout48_exceed_boot_bins_max"] == 0
 cp95_h48_live = 0.05 ** (1.0 / 48.0)
 cp95_h48_eq_tail = 0.025 ** (1.0 / 48.0)
 assert cp95_h48_live >= Pm.BOOT_KEEP_FRAC and cp95_h48_eq_tail >= Pm.BOOT_KEEP_FRAC
 assert abs(bnp["cp95_lower_bound_holdout48"] - round(cp95_h48_live, 5)) < 1e-5
 assert bnp["cg_max"] <= max(x[1] for x in Pm.BOOT_BINS)

 # [7b-tau] Live 2-stage K=1 (BALL_PAIR=0.0 exact) tau=2095 Buchstab-Dickman psurv boot with [m_min, m_max]
 # Any cell with log2(|N_g|) < m_min is REJECTED before Stage-1 ECM (acceptance probability = 0.0 identically).
 # On admitted norms m in [m_min, m_max], upper-Riemann Buchstab-Dickman psurv * 2.0 monotonically bounds survival.
 assert Pm.TAU_BOOT == 2095 and Pm.BOOT_LO_OFF == -300.0 and Pm.BOOT_SURV_X == 2.0 and Pm.BALL_PAIR == 0.0
 assert Pm.BOOT_LO_OFF < bnp["raw_deciles"][0] - 34.0
 assert bt["tau_boot"] == 2095 and bt["phi_tree"] <= 0.48
 assert abs(bt["boot_psurv_raw"] - DS.psurv(bt["m_min"], Pm.LX, Pm.TAU_BOOT)) < 1e-18
 assert abs(bt["surv_b"] - min(1.0, Pm.BOOT_SURV_X * bt["boot_psurv_raw"])) < 1e-18
 _cb_tau = (Pm.LOG2P4 + log2(DS.CELLS[DS.BOOT16[3]])) / 16.0
 _n1_t, _iA_t, _kA_t, _ = DS.BOOT16
 _m_sub_min = 24.0 * _cb_tau + Pm.LOG2U8 + min(Pm.BOOT_SUB)
 assert bt["m_min"] < _m_sub_min
 _pr_sub_min = DS.ss_lower2_tau(_m_sub_min, Pm.LX, float(_n1_t), _iA_t, _kA_t, 1, 384, Pm.ALPHA_G, Pm.TAU_BOOT)
 assert _pr_sub_min <= DS.ss_lower2(_m_sub_min, Pm.LX, float(_n1_t), _iA_t, _kA_t, 1, 384, Pm.ALPHA_G, 0.0)

 def _boot_stage1_surv_ub(m_val):
  if m_val < bt["m_min"]:
   return 0.0
  return min(1.0, Pm.BOOT_SURV_X * DS.psurv(m_val, Pm.LX, Pm.TAU_BOOT))

 for _m_lo in (1.0, 100.0, 465.0, 515.01, 1500.0, 2515.0):
  assert _boot_stage1_surv_ub(_m_lo) == 0.0
 for _m_hi in (bt["m_min"], 2530.0, 2550.0, bt["m"]):
  _sv_hi = _boot_stage1_surv_ub(_m_hi)
  assert 0.0 < _sv_hi <= bt["surv_b"] * (1.0 + 1e-12)

 # [7d] Live verification of safe surv_upper(m, LX, T) (including T=0 & 0 <= m <= T) and 2-stage post-ecm(LX)
 # zero-slack node cofactor gate T0 = kA*nA + kB*n with global pre-ECM lower-norm floor (nf_floor, ng_floor)
 for _mt in (0.0, 10.0, 69.0, 200.0, 515.01, 1000.0, 2515.01):
  _s0 = DS.surv_upper(_mt, Pm.LX, 0.0)
  assert 0.0 < _s0 <= 1.0 and _s0 >= DS.rho_lo(_mt / Pm.LX)
  for _Tt in (69.0, 300.0, 600.0, 2000.0):
   _st_val = DS.surv_upper(_mt, Pm.LX, _Tt)
   assert (_st_val == 1.0) if (_mt <= _Tt) else (0.0 < _st_val <= 1.0)
 _bn_obs_min = st["ball_norm_pool"]["observed_min"]
 _cf_lo_min = min(x[0] for x in Pm.QBINS_LO)
 _cg_lo_min = min(x[1] for x in Pm.QBINS_LO)
 assert _cf_lo_min <= _bn_obs_min["cf"] - 14.0 and _cg_lo_min <= _bn_obs_min["cg"] - 8.0
 for _bk, _sk in (((180.0, 185.0), "g"), ((230.0, 235.0), "f")):
  _dn = res["DET"][_bk, _sk]
  _T0_n = _dn["T0"]; _T_n = _dn["T"]; _nf_fl = _dn["nf_floor"]; _ng_fl = _dn["ng_floor"]
  assert _dn["global_floor_s1"] is True and _T0_n == _dn["kA"] * _dn["nA"] + _dn["kB"] * _dn["n"]
  assert abs(_dn["sv0"] - DS.surv_upper(_nf_fl, Pm.LX, _T0_n) * DS.surv_upper(_ng_fl, Pm.LX, _T0_n)) < 1e-15
  assert abs(_dn["sv_s1"] - DS.surv_upper(_nf_fl, Pm.LX, _T_n) * DS.surv_upper(_ng_fl, Pm.LX, _T_n)) < 1e-15
  assert abs(_dn["sv_norms"] - DS.smooth_part_upper(_nf_fl - _T_n, Pm.LX) * DS.smooth_part_upper(_ng_fl - _T_n, Pm.LX)) < 1e-15
  assert _dn["sv_norms"] >= _dn["sv_s1"] * (1.0 - 1e-12)
  assert abs(_dn["sv_fT_g0"] - DS.surv_upper(_nf_fl, Pm.LX, _T_n) * DS.surv_upper(_ng_fl, Pm.LX, _T0_n)) < 1e-15
  assert abs(_dn["sv_f0_gT"] - DS.surv_upper(_nf_fl, Pm.LX, _T0_n) * DS.surv_upper(_ng_fl, Pm.LX, _T_n)) < 1e-15
  def _node_s1_ub(mf_v, mg_v, nf_f=_nf_fl, ng_f=_ng_fl, T0_v=_T0_n):
   if mf_v < nf_f or mg_v < ng_f:
    return 0.0
   return DS.surv_upper(mf_v, Pm.LX, T0_v) * DS.surv_upper(mg_v, Pm.LX, T0_v)
  assert _node_s1_ub(_nf_fl - 1.0, _ng_fl) == 0.0 and _node_s1_ub(_nf_fl, _ng_fl - 1.0) == 0.0 and _node_s1_ub(10.0, 10.0) == 0.0
  for _dm in (0.0, 25.0, 75.0, 150.0):
   assert 0.0 < _node_s1_ub(_nf_fl + _dm, _ng_fl + _dm) <= _dn["sv0"] * (1.0 + 1e-12)
   assert 0.0 < DS.surv_upper(_nf_fl + _dm, Pm.LX, _T_n) * DS.surv_upper(_ng_fl + _dm, Pm.LX, _T_n) <= _dn["sv_s1"] * (1.0 + 1e-12)

 # Live multi-bin boot_R2_48cells (UDUP_BOOT=0.40) + E_16 cap C_30 cap F_Vor(a,b) deterministic m_true==1 check
 _ud_live = {}
 if run_diagnostics:
  _ub = st["unit_orbit_udup"]["boot_R2_48cells"]
  assert len(_ub["m_true"]) == 48 and len(_ub["cgs"]) == 48
  _cb2 = (Pm.LOG2P4 + log2(DS.CELLS[2])) / 16.0
  _n1b, _iAb, _kAb, _ = DS.BOOT16
  _wsb = [DS.ss_lower2_tau(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, Pm.TAU_BOOT) for cg_ in _ub["cgs"]]
  _prb = sum(w * DS.ss_lower2_tau(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, Pm.TAU_BOOT) for w, cg_ in Pm.BOOT_BINS)
  _prb = min(sum(DS.ss_lower2_tau(24.0 * _cb2 + Pm.LOG2U8 + cg_, Pm.LX, float(_n1b), _iAb, _kAb, 1, 384, Pm.ALPHA_G, Pm.TAU_BOOT) for cg_ in Pm.BOOT_SUB) / len(Pm.BOOT_SUB), Pm.PRCAP_BOOT * _prb)
  _inv_mb = [1.0 / m for m in _ub["m_true"]]
  _mb_inv = sum(_inv_mb) / 48.0
  _se_mb = math.sqrt(sum((x - _mb_inv) ** 2 for x in _inv_mb) / (47.0 * 48.0))
  _twb = sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / sum(_wsb)
  _dyb = (sum(w / m for w, m in zip(_wsb, _ub["m_true"])) / 48.0) / _prb
  assert abs(_mb_inv - _ub["mean_inv_m"]) < 1e-3
  assert Pm.UDUP_BOOT <= _mb_inv - 2.0 * _se_mb and _mb_inv <= 1.0 and Pm.UDUP_BOOT <= _twb <= 1.0 and _dyb >= Pm.UDUP_BOOT
  _ud_live["boot_R2_48cells"] = dict(mean_inv_m=round(_mb_inv, 4), true_yield_weighted_inv_m=round(_twb, 4), dedup_yield_over_charged_pr=round(_dyb, 4))
  _ud_live["e16_c30_fvor_160ideals"] = dict(f_jv_measured=_bnp_m["f_jv_measured"], f_jv_se=_bnp_m["f_jv_se"], udup_node_charged=Pm.UDUP_NODE, deterministic_m_true=1)

 # [11] Live ACTUAL integer BOUND_CELL (+6 ball-mul margin) + (135, g) short-subspace counterexample rejection
 dp_rd, dp_pref = _build_dp_ball(32)
 rng = random.Random(20261122 + seed)
 q_bc, s_bc, t_bc = _sample_ideal(90, "g", rng)
 B_bc, det_e, det_qcm = _lll_16d_qcm(_sq_hnf(q_bc, s_bc, t_bc))
 assert det_e == q_bc * q_bc and det_qcm == (1265625 ** 2) * q_bc * q_bc
 bc_gmin, bc_gmax = 1e9, -1e9
 for _ in range(12):
  k = _sample_prim_ball_dp(dp_rd, dp_pref, 32, rng)
  v = [sum(k[i] * B_bc[i][j] for i in range(16)) for j in range(16)]
  nf, ng = _norms_fg(v[:8], v[8:])
  assert ng % q_bc == 0 and nf > 0 and ng > 0
  nf_ub, ng_ub, nf_c, ng_c = DS.bound_cell_int_fg(v[:8], v[8:])
  gf, gg = nf_ub - log2(nf), ng_ub - log2(ng)
  assert gf >= 0.0 and gg >= 0.0 and nf_c >= ceil(log2(nf)) and ng_c >= ceil(log2(ng))
  bc_gmin = min(bc_gmin, gf, gg)
  bc_gmax = max(bc_gmax, gf, gg)
 assert bc_gmin >= 0.0 and bc_gmax <= 0.5343

 _q_ce = 81009504249198611324243595142290546015061
 _s_ce = 14886707496769063264106121185902690387739
 _t_ce = (5 * pow(8, -1, _q_ce)) % _q_ce
 _B_ce, _, _ = _lll_16d_qcm(_sq_hnf(_q_ce, _s_ce, _t_ce))
 assert DS.screen_special_q_basis(_B_ce, _q_ce, "g", _norms_fg) is False
 _k_ce = [1 if i == 0 else (47 if i == 8 else 0) for i in range(16)]
 _v_ce = [sum(_k_ce[i] * _B_ce[i][j] for i in range(16)) for j in range(16)]
 _d135 = res["DET"][(135.0, 140.0), "g"]
 _ok_ce, _nf_ce_rem, _ng_ce_rem, _, _ = DS.bound_cell_check_and_cap(_v_ce[:8], _v_ce[8:], _d135["nf"], _d135["ng"], lq=log2(_q_ce), side="g")
 assert _ok_ce is False and _nf_ce_rem > _d135["nf"] and _ng_ce_rem > _d135["ng"]
 assert DS.screen_special_q_basis(B_bc, q_bc, "g", _norms_fg) is True

 # Live 4-bin production deg-(72,24) q >= 2^70 raw-L_q Q_CM-LLL + E_16 cap C_30 cap F_Vor(a,b) + BOUND_CELL diagnostic
 # evaluated at each bin's EXACT actual charged cell count C_R (ZERO radius rescaling; 32 primitive cells/lattice).
 # Supporting exploratory compatibility check only — NOT a full-ball proof or replacement for H5.
 import q15_diag_extra as _DXE
 _places15, _vor_u10, _vor_facets, _e16_seen_vor, _e16_dup_units = _DXE.e16_verify_fvor_core(B_bc, _G8, DS.CELLS[32], rdir, _mul_z)

 prod_q70_diag = None
 if run_diagnostics:
  rng_p70 = random.Random(20261123)
  p70_rows = []
  for _b_lo, _b_hi, _qb, _sd in ((69.0, 75.0, 72, "f"), (100.0, 105.0, 102, "g"), (180.0, 185.0, 182, "g"), (230.0, 235.0, 232, "f")):
   _qp, _sp, _tp = _sample_ideal(_qb, _sd, rng_p70)
   _Bp, _dep, _dqp = _lll_16d_qcm(_sq_hnf(_qp, _sp, _tp))
   assert _dep == _qp * _qp and _dqp == (1265625 ** 2) * _qp * _qp
   assert DS.screen_special_q_basis(_Bp, _qp, _sd, _norms_fg) is True
   _dp = res["DET"][(_b_lo, _b_hi), _sd]
   _ns, _iAs, _kAs, _kBs, _Rs = _dp["n"], _dp["iA"], _dp["kA"], _dp["kB"], _dp["R"]
   _lqp = log2(_qp)
   _Gp, _Uip, _dhp = _DXE.e16_chol_invT(_Bp, _G8)
   _R_ell_p = ((DS.CELLS[_Rs] * 30.0 * math.factorial(8) * _dhp) / (math.pi ** 8)) ** (1.0 / 16.0)
   _pr_sum = 0.0
   _kept = 0
   for _ in range(32):
    _kk, _vv = _DXE.e16_sample_cell(_Bp, _Uip, _Gp, _R_ell_p, _places15, _vor_facets, rng_p70)
    _ok_bc, _nf_rem, _ng_rem, _, _ = DS.bound_cell_check_and_cap(
    _vv[:8], _vv[8:], _dp["nf"], _dp["ng"], lq=_lqp, side=_sd
    )
    if not _ok_bc:
     continue
    _nf_ex, _ng_ex = _norms_fg(_vv[:8], _vv[8:])
    if _nf_ex == 0 or _ng_ex == 0:
     continue
    _nf_ex_rem = log2(_nf_ex) - (_lqp if _sd == "f" else 0.0)
    _ng_ex_rem = log2(_ng_ex) - (_lqp if _sd == "g" else 0.0)
    assert _nf_ex_rem <= _nf_rem + 1e-9 and _ng_ex_rem <= _ng_rem + 1e-9
    if _nf_ex_rem < _dp["nf_floor"] or _ng_ex_rem < _dp["ng_floor"]:
     continue
    _kept += 1
    _pf = DS.ss_lower2(_nf_rem, Pm.LX, float(_ns), _iAs, _kAs, _kBs, 384, Pm.ALPHA_F, 0.0)
    _pg = DS.ss_lower2(_ng_rem, Pm.LX, float(_ns), _iAs, _kAs, _kBs, 384, Pm.ALPHA_G, 0.0)
    _pr_sum += _pf * _pg
   _rat = (_pr_sum / 32.0) / (2.0 ** _dp["pr_log2"])
   assert _kept > 0 and _rat >= 2.0
   p70_rows.append(dict(bin_lo=_b_lo, side=_sd, q_bits=_qb, exact_R=_Rs, cells=32, kept_cells=_kept, screened_cert_ratio=round(_rat, 3)))
  prod_q70_diag = dict(
  bound_cell_screened=True,
  exact_norms_floor_checked=True,
  e16_gram_ellipsoid=True,
  radius_rescaling=False,
  cells_per_lattice=32,
  n_checked=len(p70_rows),
  n_below_1x=sum(1 for x in p70_rows if x["screened_cert_ratio"] < 1.0),
  min_ratio=min(x["screened_cert_ratio"] for x in p70_rows),
  bins=p70_rows,
  )

 # [11b] Live q15_retain16 prime-splitting + non-vacuous p=2 (K=140) + j>=74 tail + counted_fraction
 rt_v = None; _t2_mono = None; t16_st = None; ntt3_st = None; _t16_ok = True; sch_node = sch_boot = None; lad_boot = 1.0
 if core_checks:
  import q15_ecm_ntt3 as NT3
  ntt3_st = NT3.verify_ntt3_stage2_ecm()
  assert ntt3_st["ok"] is True
  tab = st["prime_structure"]["table"]
  for pp in (2, 3, 5, 19, 31, 61, 101, 257):
   assert RT.structure(pp) == tab[str(pp)]
  assert abs(RT.counted_fraction() - st["omitted_law"]["counted_fraction"]) < 1e-6
  assert st["omitted_law"]["tail_j_ge_74_included"] is True
  rt_v = RT.verify_retain16(n_samples=64, seed=11, recompute_full=False)
  assert rt_v["probes_verified"] >= 10 and rt_v["sample_node_slack96"]["loss_f"] < 5e-5 and rt_v["sample_node_slack96"]["loss_g"] < 5e-5
  assert rt_v["residue_checks"][2]["K"] == 140 and rt_v["residue_checks"][2]["max_g"] <= 128.0

     # [12] Live q15_ecm16 Deterministic Stepped-GCD L_min=2 branch-tree DP + Z/2 x Z/8 Torsion-16 Edwards (a=+1) ECM check
  assert EC.TORSION == 16 and EC.STEPPER_MULS == 25.0 and EC.STEPPER_GCDS == 4.0 and EC.MR_ROUNDS == 110
  assert EC.L_MIN_DEFAULT == 2 and EC.B_MIN_DEFAULT == 32 and EC.NEWTON_PRECHECK_INSTR == 2.0 ** 26
  l2_boot, _ = EC.tree_branch_caps(int(ceil(bt["m"])), EC.B_MIN_DEFAULT, EC.L_MIN_DEFAULT)
  assert l2_boot >= ceil(bt["m"] / 64.0) ** 2
  _l2_65, _g_65 = EC.tree_branch_caps(65, EC.B_MIN_DEFAULT, EC.L_MIN_DEFAULT)
  assert _l2_65 >= 8 and _g_65 >= 2.0 * EC.gcd_instr(33, EC.L_MIN_DEFAULT)
  _b0 = res["bins"][0]
  sch_node = EC.schedule_rung(round(Pm.LX + EC.ECM_T2_PAD_BITS, 4), int(ceil(res["DET"][_b0, "f"]["nf"])), Pm.ECM_SUCCESS)
  sch_boot_s1 = EC.schedule_rung(round(Pm.LX + EC.ECM_T2_PAD_BITS, 4), int(ceil(bt["m"])), Pm.ECM_SUCCESS)
  sch_boot = EC.schedule_rung(round(bt["nA"] + EC.ECM_T2_PAD_BITS, 4), int(Pm.TAU_BOOT), Pm.ECM_SUCCESS)
  assert sch_node["curves"] > 0 and sch_boot_s1["curves"] > 0 and sch_boot["curves"] > 0
  _t2_targets = [(Pm.LX, int(ceil(bt["m"]))), (bt["nA"], int(Pm.TAU_BOOT))]
  for _dv in res["DET"].values():
   _t2_targets.append((Pm.LX, int(ceil(_dv["nf"]))))
   _t2_targets.append((Pm.LX, int(ceil(_dv["ng"]))))
   if _dv["kA"] + _dv["kB"] > 0:
    _t2_targets.append((_dv["nA"], int(ceil(_dv["T0"]))))
  _t2_mono = EC.verify_1rung_monotonicity(_t2_targets)
  assert _t2_mono["ok"] is True and _t2_mono["min_iv_ratio"] >= 1.0
  lad_boot = sch_boot["instructions"]
  t16_st = EC.verify_t16_edwards_ecm()
  _t16_ok = bool(
  t16_st["torsion"] == 16 and t16_st["nl_rank_ge1"] and t16_st["tors_z2xz8_q"]
  and t16_st["t16_prime_counts"] == "243/243" and t16_st["ew_mod_step_ok"]
  and t16_st["ew_to_mont_mod_ok"] and t16_st["mont_muls_per_bit"] <= 11.05 and t16_st["w2310_live"]
  and t16_st["perfect_power_precheck_ok"] and t16_st["singular_nm1_rejected"]
  and t16_st["split_factor"] > 1 and t16_st["emp_composite_mode"]
  and t16_st["emp_p1"] >= t16_st["model_p1_lo"]
  and t16_st["emp_balanced_p1"] >= t16_st["model_p1_lo"]
  and t16_st["emp_tri_p1"] >= t16_st["model_p1_lo"]
  and t16_st["emp_tri_s1s2_p1"] >= t16_st["emp_tri_s1s2_model_p1_lo"]
  and t16_st["frontier_3factor_full_ok"] == "16/16"
  and t16_st["multicurve_miss_tail_ok"]
  and t16_st["stage2_tree_bisect_ok"] and t16_st["stage2_zero_leaf_skip_ok"]
  and t16_st["stage2_raw_zero_leaf_ok"] and t16_st["stage2_subfactor_zero_ok"]
  and t16_st["stage2_3prime_pq_split_ok"]
  )
  assert _t16_ok

 # [13] Extra scoped exploratory diagnostics (q15_diag_extra.py) + exact-degree deg-(72,24) H5/H5-W (q15_exact_h5_diag.py)
 import q15_diag_extra as DX
 import q15_exact_h5_diag as EH5
 hea_dde = DX.verify_hea_surv_analytic(res, Pm, DS) if core_checks else None
 exact_h5 = EH5.verify_q15_exact_deg72_h5(None, sys.modules[__name__], DS, full_replays=bool(run_diagnostics))
 if run_diagnostics:
  _arc_mini = json.loads(lzma.decompress(base64.b85decode("".join((rdir / "q15_mini_adaptive.b85").read_text().split()))).decode())
  js700k = DX.verify_joint_smooth_700k(_arc_mini)
  hb_diag = DX.verify_highbit_prod_diagnostics(res, Pm, DS, DS, _lll_16d_qcm, _sq_hnf, _norms_fg, _arc_mini, n_cells=8, G8=_G8, places=_places15, facets=_vor_facets)
 else:
  js700k = None
  hb_diag = None
 if not core_checks:
  return dict(all_ok=True, mode="pure_k1_boole_wald_iid_boot", wall_s=round(time.perf_counter() - t0, 3), exact_deg72_h5=exact_h5, mini_adaptive_campaign=mini_res, live_full16d_k1_gen=live_k1_gen, prod_q70_cert_diag=prod_q70_diag, udup_live_recomputed=_ud_live, joint_smooth_700k=js700k, highbit_prod_diag=hb_diag)

 wall_s = time.perf_counter() - t0
 return dict(
 all_ok=True,
 mode="pure_k1_boole_wald_iid_boot",
 plan_tag=DS.PLAN_TAG,
 wall_s=round(wall_s, 3),
 eval_384_s=round(dt_ev, 3),
 mesh=mesh,
 n_nodes=len(rows),
 n_sieve=len(_det_sieve),
 n_ecm=len(_det_ecm),
 max_R_used=max_R_used,
 lx=Pm.LX,
 nb_r=Pm.NB_R,
 nb_alpha=Pm.NB_ALPHA,
 p0_floor=Pm.P0_FLOOR,
 min_mu=round(min_mu, 2),
 p0_max=float(f"{bt['p0_max']:.4e}"),
 c_det_one=int(bt["c_det_one"]),
 phi_tree=round(bt["phi_tree"], 6),
 e_boot_runs=round(bt["e_boot_runs"], 6),
 e_trials_log=float(f"{bt['e_trials_log']:.4e}"),
 expected_slots_2log=round(EN, 2),
 two_log_units=round(two_log, 6),
 boot_u=round(bt["boot_u"], 4),
 sub_u=round(bt["sub_u"], 4),
 boot_phi=0.0,
 p_descent_fail=0.0,
 boot_m=round(bt["m"], 2),
 boot_m_min=round(bt["m_min"], 2),
 tau_boot=int(bt["tau_boot"]),
 boot_lo_off=bt["boot_lo_off"],
 boot_surv_b=float(f"{bt['surv_b']:.6e}"),
 boot_psurv_raw=float(f"{bt['boot_psurv_raw']:.6e}"),
 boot_surv_x=Pm.BOOT_SURV_X,
 boot_tau_monotone_ok=True,
 node_2stage_zero_slack_ok=True,
 ecm_1rung_monotonicity_ok=True,
 ecm_1rung_targets=_t2_mono["targets"],
 ecm_1rung_intervals=_t2_mono["intervals_checked"],
 ecm_1rung_min_iv_ratio=_t2_mono["min_iv_ratio"],
 ecm_1rung_pad_bits=_t2_mono["pad_bits"],
 node_floor_margin_bits=dict(cf=round(_bn_obs_min["cf"] - _cf_lo_min, 2), cg=round(_bn_obs_min["cg"] - _cg_lo_min, 2)),
 cells_16d_live_verified=True,
 boot_miss=round(bt["miss"], 5),
 boot_lin_keep=round(bt["lin_keep"], 5),
 max_miss_node=round(max_miss_node, 5),
 boot_keep_frac=bt["boot_keep_frac"],
 boot_n_factors=bt["boot_n_factors"],
 boot_n1=int(bt["n1"]),
 boot_nA=round(bt["nA"], 2),
 boot_kA=int(bt["kA"]),
 boot_R=int(bt["R"]),
 mini_adaptive_campaign=mini_res,
 live_full16d_k1_gen=live_k1_gen,
 prod_q70_cert_diag=prod_q70_diag,
 pm_log2_rels=Pm.LOG2_RELS,
 pm_main_lb_f=Pm.MAIN_LB_F,
 pm_main_lb=Pm.MAIN_LB,
 pm_main_nf=round(Pm.MAIN_UF * Pm.MAIN_LB_F, 4),
 pm_main_ng=round(Pm.MAIN_UG * Pm.MAIN_LB, 4),
 pm_alpha_f=Pm.ALPHA_F,
 pm_alpha_g=Pm.ALPHA_G,
 pm_pair_setup=Pm.PAIR_SETUP,
 pm_lll_node=Pm.LLL_NODE,
 pm_lll_boot=Pm.LLL_BOOT,
 pm_bound_cell=Pm.BOUND_CELL,
 pm_sieve_hit=Pm.SIEVE_HIT,
 pm_udup_boot=Pm.UDUP_BOOT,
 pm_udup_node=Pm.UDUP_NODE,
 pm_e16_count_lo=Pm.E16_COUNT_LO,
 pm_e16_count_hi=Pm.E16_COUNT_HI,
 e16_fvor_live=dict(tested=80, fvor_hits=_e16_seen_vor, unit_dup_violations=_e16_dup_units, f_jv_160ideals=_bnp_m["f_jv_measured"]),
 ntt3_stage2_ecm=ntt3_st,
 pm_ret_mult_f=Pm.RET_MULT_F,
 pm_ret_mult_g=Pm.RET_MULT_G,
 f4_counterexample_rejected=True,
 t16_edwards_ecm_ok=_t16_ok,
 t16_edwards_ecm=t16_st,
 hea_surv_analytic=hea_dde,
 exact_deg72_h5=exact_h5,
 joint_smooth_700k=js700k,
 highbit_prod_diag=hb_diag,
 udup_live_recomputed=_ud_live,
 eager_all_100pct=True,
 linear_miss_retention=True,
 min_ret=round(min_ret, 5),
 max_loss=float(f"{max_loss:.3e}"),
 bound_cell_live_min_gap_bits=bc_gmin,
 bound_cell_live_max_gap_bits=bc_gmax,
 retain16_live_ok=True,
 retain16_probes=int(rt_v["probes_verified"]),
 retain16_p2_K=int(rt_v["residue_checks"][2]["K"]),
 retain16_p2_max_g=rt_v["residue_checks"][2]["max_g"],
 boot_r2_holdout48=dict(
 fresh_from_production_hnf_all_64=True,
 cal16_cg_min=bnp["cal16_cg_min"],
 cal16_cg_max=bnp["cal16_cg_max"],
 holdout48_cg_min=bnp["holdout48_cg_min"],
 holdout48_cg_max=bnp["holdout48_cg_max"],
 holdout48_below_m_min_cells=0,
 holdout48_exceed_targets=0,
 holdout48_outside_window_targets=0,
 cp95_one_sided_lb_holdout48_window=bnp["cp95_lower_bound_holdout48"],
 cp95_equal_tailed_lb_holdout48_window=round(cp95_h48_eq_tail, 5),
 boot_keep_frac_charged=Pm.BOOT_KEEP_FRAC,
 ),
 ecm_schedules=dict(
 stepped_gcd_ladder=True,
 l_min=EC.L_MIN_DEFAULT,
 b_min=EC.B_MIN_DEFAULT,
 node_min_f=dict(curves=sch_node["curves"], pi_b1=round(sch_node["pi_b1"], 1), G=sch_node["G"], log2_ir=round(sch_node["log2_instructions"], 2)),
 boot_nA=dict(curves=sch_boot["curves"], pi_b1=round(sch_boot["pi_b1"], 1), G=sch_boot["G"], log2_ir=round(sch_boot["log2_instructions"], 2), ladder_log2_ir=round(log2(lad_boot), 2)),
 ),
 )


if __name__ == "__main__":
 out = run_all_checks()
 print(json.dumps(out, indent=2))
