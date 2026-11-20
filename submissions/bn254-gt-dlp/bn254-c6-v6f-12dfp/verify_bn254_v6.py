"""
verify_bn254_v6.py — Self-contained pure-stdlib live prototypes for BN254 v6:
 1. 8-lane SWAR sb_u7 Band-2 prototype (0 mismatches / 1523 hits) + 10k 64-bit uint16 SWAR scan+zero
  identity + exact 380-chunk 8-limb P_14 score-only Stage-B strip check (ZERO prime recording!);
 2. Live sigma-orbit-folded 12D Ellipsoid Fincke-Pohst Band-3 prototype across 4 real BN254 C6
  prime sigma-orbits (602/602 canonical hits, 0 missed, 0 duplicates);
 3. Exact 12D C6 special-q (6 bases 3f+3g spanning 2^66..2^213) & 2 boot lattice audit using exact
  FULL B_12(4) uniform sampling (13,315,520 cells), exact 6x6 integer norms, exact q|N_f & q|N_g,
  paired joint smoothness vs certified C6 QBINS/BOOT_BINS (zero PRCAP);
 4. Pure-Python production hexagonal Z[omega] filtering simulator at (61, 60, 0.82, W_ROW=161.0)
  with deterministic pre-merge weight guard, exact kernel preservation, live w_mid -> 0 monotone
  check, and recorded-pivot export for hx_z6.py!
"""
import bisect, heapq, json, math, random
from math import ceil, exp, floor, log2, pi, sqrt
from pathlib import Path
import model_bn254_v6 as M
import descent_bn254_v6 as DESC
def _mv6(Mat, v): return [sum(Mat[r][c] * v[c] for c in range(6)) for r in range(6)]
def _mm6(A, B): return [[sum(A[i][k] * B[k][j] for k in range(6)) for j in range(6)] for i in range(6)]
def _qf6(G, v):
 s = 0.0
 for i in range(6):
  vi = v[i]
  if vi:
   s += G[i][i] * vi * vi
   for j in range(i + 1, 6):
    if v[j]: s += 2.0 * G[i][j] * vi * v[j]
 return s
def _eval_poly_mod(cs, s, p):
 v = 0
 for c in reversed(cs): v = (v * s + c) % p
 return v
def _eval_zk_mod(vec, s, p):
 inv3 = pow(3, -1, p); sp = [1]
 for _ in range(5): sp.append((sp[-1] * s) % p)
 return (sum(vec[k] * (sum(M.ZK_NUM[k][j] * sp[j] for j in range(6)) % p) for k in range(6) if vec[k]) * inv3) % p
def _is_prime(n): return n >= 2 and all(n % d for d in range(2, int(n**0.5) + 1))
def verify_swar_sb_u7_and_stageB():
 D = M.load_data(); T4 = D["t4_geom"]; KAP = D["pool_quad"]["kap"]
 GA, GB, S, FAC6, UMATS = T4["GA"], T4["GB"], T4["SMAT"], T4["FAC6"], T4["UMATS"]
 EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
 EMB_W = [complex(r, im) for r in [0] for im in [0]]
 EMB_W = [complex(r, im) for r, im in T4["EMB_W"]]
 S2 = _mm6(S, S)
 rng = random.Random(20261120)
 swar_errs = 0
 MASK = 0x8000800080008000; P2_BIAS = 0x0030003000300030
 for _ in range(10000):
  t_a = [rng.randint(500, 12000) for _ in range(4)]
  b_in_p2 = (rng.random() < (1.0 / 9.0))
  sc = [rng.randint(0, 12500) for _ in range(4)]
  v64 = sc[0] | (sc[1] << 16) | (sc[2] << 32) | (sc[3] << 48)
  b4 = [(0x8000 - t) for t in t_a]
  bias64 = b4[0] | (b4[1] << 16) | (b4[2] << 32) | (b4[3] << 48)
  if b_in_p2: v64 += P2_BIAS
  swar_hit = (((v64 + bias64) & MASK) != 0)
  shift = 48 if b_in_p2 else 0
  if swar_hit != any((sc[i] + shift) >= t_a[i] for i in range(4)):
   swar_errs += 1
 max_inv = 0.0
 for _ in range(500):
  a = [rng.randint(-450, 450) for _ in range(6)]
  b = [rng.randint(-750, 750) for _ in range(6)]
  if not any(a) or not any(b): continue
  lab = la0 = 0.0
  for j in range(3):
   za = sum(EMB_ZK[j][k] * a[k] for k in range(6))
   zb = sum(EMB_ZK[j][k] * b[k] for k in range(6))
   uw_a2 = M.U_BN254 * EMB_W[j] * (za * za); vab = uw_a2 - zb * zb
   lab += log2(vab.real**2 + vab.imag**2); la0 += log2(uw_a2.real**2 + uw_a2.imag**2)
  max_inv = max(max_inv, abs(lab - la0))
 p14 = M._primes_upto(16383)
 chunks, cur_p, cur_list = [], 1, []
 for p in p14:
  if cur_p * p >= (1 << 64):
   chunks.append((cur_p, cur_list)); cur_p, cur_list = p, [p]
  else:
   cur_p *= p; cur_list.append(p)
 if cur_list: chunks.append((cur_p, cur_list))
 strip_errs = 0
 for _ in range(100):
  n_test = rng.getrandbits(469) | (1 << 468)
  for _ in range(4): n_test *= p14[rng.randrange(len(p14))]
  rem_chunk = n_test; m_chunk = 1
  for prod64, plist in chunks:
   g = math.gcd(rem_chunk % prod64, prod64)
   if g > 1:
    for p in plist:
     while rem_chunk % p == 0:
      rem_chunk //= p; m_chunk *= p
  rem_ref = n_test; m_ref = 1
  for p in p14:
   while rem_ref % p == 0:
    rem_ref //= p; m_ref *= p
  if m_chunk != m_ref or rem_chunk != rem_ref:
   strip_errs += 1
 def is_can(a):
  if tuple(-x for x in a) < tuple(a): return False
  xv = [2.0 * ((z := sum(EMB_ZK[j][k] * a[k] for k in range(6))).real**2 + z.imag**2) for j in range(3)]
  for fi in range(6):
   dot = sum(FAC6[fi][j] * xv[j] for j in range(3))
   if dot < -1e-9: return False
   if abs(dot) <= 1e-9:
    uv = _mv6(UMATS[fi], a); nuv = [-x for x in uv]
    if tuple(nuv if tuple(nuv) < tuple(uv) else uv) < tuple(a): return False
  for Sm in (S, S2):
   vk = _mv6(Sm, a); nvk = [-x for x in vk]
   if tuple(nvk if tuple(nvk) < tuple(vk) else vk) < tuple(a): return False
  return True
 ideals = []
 for p in range(101, 1800, 2):
  if not _is_prime(p): continue
  for s in range(p):
   if _eval_poly_mod(M.H_COEFFS, s, p) == 0:
    uw = ((M.U_BN254 % p) * _eval_poly_mod(M.W_COEFFS, s, p)) % p
    if uw == 0: continue
    inv_uw = pow(uw, -1, p)
    if pow(inv_uw, (p - 1) // 2, p) == 1:
     for r in range(1, (p + 1) // 2):
      if (r * r) % p == inv_uw: ideals.append((p, s, r)); break
  if len(ideals) >= 24: break
 RA2, RB2 = 36.0, 110.0
 a_pts = []
 for a5 in range(-2, 3):
  for a4 in range(-2, 3):
   for a3 in range(-2, 3):
    for a2 in range(-2, 3):
     for a1 in range(-2, 3):
      for a0 in range(-2, 3):
       av = [a0, a1, a2, a3, a4, a5]
       if any(av) and (qv := _qf6(GA, av)) <= RA2 and is_can(av):
        a_pts.append((qv / RA2, av))
 a_pts.sort(key=lambda x: x[0]); a_pts = a_pts[:80]
 b_rows, b_pts = [], []
 for b5 in range(1, 3):
  for b4 in range(-2, 3):
   for b3 in range(-2, 3):
    for b2 in range(-2, 3):
     for b1 in range(-2, 3):
      rp = [(b0, bv, qv / RB2) for b0 in range(-6, 7) if (qv := _qf6(GB, (bv := [b0, b1, b2, b3, b4, b5]))) <= RB2]
      if rp:
       b_rows.append((len(rp), [x[1] for x in rp], [x[2] for x in rp]))
       for _, bv, qn in rp: b_pts.append((bv, qn, len(b_pts)))
      if len(b_rows) >= 50: break
     if len(b_rows) >= 50: break
    if len(b_rows) >= 50: break
   if len(b_rows) >= 50: break
  if len(b_rows) >= 50: break
 gt = {(ai, bi): 0 for ai, (qan, _) in enumerate(a_pts) for bi, (_, qbn, _) in enumerate(b_pts) if qbn <= 1.0 - KAP * qan + 1e-12}
 tot_hits = 0
 for p, s, rho_r in ideals:
  ar = [_eval_zk_mod(av, s, p) for _, av in a_pts]
  br = [_eval_zk_mod(bv, s, p) for bv, _, _ in b_pts]
  for ai, bi in gt:
   rb = (rho_r * br[bi]) % p
   c = int(ar[ai] == rb) + int(ar[ai] == (-rb) % p)
   if c: gt[(ai, bi)] += c; tot_hits += c
 W_SB = 4; sc = {k: 0 for k in gt}; b_map = {tuple(bv): bi for bv, _, bi in b_pts}
 sb_total = sb_occ = sb_ent = 0; HMASK = 0x8080808080808080
 for p, s, rho_r in ideals:
  inv_r = pow(rho_r, -1, p); n_sb = (p + W_SB - 1) // W_SB
  tab = {}; sb_u7 = bytearray(n_sb + 8)
  for ai, (qan, av) in enumerate(a_pts):
   kta = (_eval_zk_mod(av, s, p) * inv_r) % p; q = kta // W_SB
   tab.setdefault(q, []).append((qan, kta, ai))
  for q, lst in tab.items():
   lst.sort(key=lambda x: x[0])
   sb_u7[q] = max(1, min(127, 127 - int(floor(126.0 * lst[0][0]))))
  for L_b, b_vecs, qb_norms in b_rows:
   k0 = _eval_zk_mod(b_vecs[0], s, p)
   for sign in (1, -1):
    lo = k0 if sign == 1 else (-(k0 + L_b - 1)) % p
    hi = (k0 + L_b - 1) % p if sign == 1 else (-k0) % p
    rngs = [(lo // W_SB, hi // W_SB)] if lo <= hi else [(lo // W_SB, n_sb - 1), (0, hi // W_SB)]
    for q_lo, q_hi in rngs:
     n_span = q_hi - q_lo + 1; sb_total += n_span
     assert n_span <= 8, f"SWAR overflow: n_span={n_span} > 8"
     bias64 = 0; qa_lim_sb = []
     for j in range(n_span):
      q = q_lo + j
      t_vals = [t for r in range(q * W_SB, min(p, (q + 1) * W_SB))
                 if (t := ((r - k0) % p if sign == 1 else (-r - k0) % p)) < L_b]
      if not t_vals: qa_lim_sb.append(-1.0); continue
      qa_max_j = min(1.0, max(0.0, (1.0 - min(qb_norms[t] for t in t_vals)) / KAP + 1e-9))
      qa_lim_sb.append(qa_max_j)
      thr_j = max(1, min(127, 127 - int(ceil(126.0 * qa_max_j))))
      bias64 |= ((128 - thr_j) << (8 * j))
     w64 = 0
     for j in range(n_span):
      if sb_u7[q_lo + j] > 0: sb_occ += 1
      w64 |= (sb_u7[q_lo + j] << (8 * j))
     m = (w64 + bias64) & HMASK
     if m == 0: continue
     for j in range(n_span):
      if (m >> (8 * j + 7)) & 1:
       sb_ent += 1; q = q_lo + j; qa_max_j = qa_lim_sb[j]
       for qan, kta, ai in tab[q]:
        if qan > qa_max_j + 1e-12: break
        t = (kta - k0) % p if sign == 1 else (-kta - k0) % p
        if t < L_b and qb_norms[t] <= 1.0 - KAP * qan + 1e-12:
         sc[(ai, b_map[tuple(b_vecs[t])])] += 1
 mismatches = sum(1 for k in gt if sc[k] != gt[k])
 win_fn = 0
 c_ref = floor(468 / 14)
 win_rngs = [(0, c_ref + 1), (16 * 59 - 1, 16 * 60 + c_ref + 1), (16 * 118 - 1, 16 * 119)]
 win_disjoint = (win_rngs[0][1] < win_rngs[1][0] and win_rngs[1][1] < win_rngs[2][0])
 for k_lp in (0, 1, 2):
  for _ in range(50):
   lp_bits = sum(rng.uniform(59.0, 60.0) for _ in range(k_lp))
   sm_bits = rng.uniform(5.0, max(6.0, 110.0 - lp_bits))
   rem_sieve = 468.0 - lp_bits - sm_bits
   n_f = rng.randint(8, min(c_ref, int(rem_sieve / 14.0)))
   ps_b = [rem_sieve / n_f] * n_f
   s_q = sum(floor(16.0 * b) for b in ps_b)
   O_val = floor(16.0 * 468.0) - s_q - floor(16.0 * sm_bits)
   if O_val <= 16 * 119:
    if not any(lo <= O_val <= hi for lo, hi in win_rngs):
     win_fn += 1
 p2_errs = 0
 for _ in range(1000):
  L_raw = rng.randint(16 * 400, 16 * 500)
  T_slot = L_raw - 16 * 114
  bias_slot = 0x8000 - T_slot
  s_cell = rng.randint(0, T_slot + 200)
  swar_pass = (((s_cell + 48 + bias_slot) & 0x8000) != 0)
  odd_pass = (((L_raw - 48) - s_cell) <= 16 * 114)
  if swar_pass != odd_pass: p2_errs += 1
 return dict(
  ok=(mismatches == 0 and swar_errs == 0 and strip_errs == 0 and win_disjoint and win_fn == 0 and p2_errs == 0 and max_inv <= 1e-12 and tot_hits == 1523),
  scalar_swar_10k_errors=swar_errs,
  p2_shift_threshold_errors=p2_errs,
  max_ng_invariance_bits=max_inv,
  p14_primes=len(p14), p14_word_chunks=len(chunks), p14_strip_errors=strip_errs,
  stageB_window_disjoint=win_disjoint, stageB_window_false_neg=win_fn,
  sb_u7_hits=tot_hits, sb_u7_mismatches=mismatches,
  sb_total=sb_total, sb_occ=sb_occ, sb_ent=sb_ent,
  sb_occ_skip_pct=round(100.0 * (1.0 - sb_ent / max(1, sb_occ)), 2),
 )
def verify_12d_sigma_folded_walk():
 """Live 12D Optimal Bounding Ellipsoid E* Fincke-Pohst across 4 real BN254 C6 prime sigma-orbits on EXACT active R (u<=1, v+kap*u<=1; NO /KAP!)."""
 D = M.load_data(); T4 = D["t4_geom"]; KAP = D["pool_quad"]["kap"]
 GA, GB, S, FAC6, UMATS = T4["GA"], T4["GB"], T4["SMAT"], T4["FAC6"], T4["UMATS"]
 EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
 S2 = _mm6(S, S)
 def is_in_unit_fd_closed(a):
  xv = [2.0 * ((z := sum(EMB_ZK[j][k] * a[k] for k in range(6))).real**2 + z.imag**2) for j in range(3)]
  return all(sum(FAC6[fi][j] * xv[j] for j in range(3)) >= -1e-9 for fi in range(6))
 def is_in_unit_fd(a):
  if tuple(-x for x in a) < tuple(a): return False
  xv = [2.0 * ((z := sum(EMB_ZK[j][k] * a[k] for k in range(6))).real**2 + z.imag**2) for j in range(3)]
  for fi in range(6):
   dot = sum(FAC6[fi][j] * xv[j] for j in range(3))
   if dot < -1e-9: return False
   if abs(dot) <= 1e-9:
    uv = _mv6(UMATS[fi], a); nuv = [-x for x in uv]
    if tuple(nuv if tuple(nuv) < tuple(uv) else uv) < tuple(a): return False
  return True
 def canonicalize_sigma_all(a, b):
  imgs = []
  for Sm in (None, S, S2):
   ak = a if Sm is None else _mv6(Sm, a)
   bk = b if Sm is None else _mv6(Sm, b)
   nak = [-x for x in ak]
   if tuple(nak) < tuple(ak): ak, bk = nak, [-x for x in bk]
   imgs.append((tuple(ak), tuple(bk)))
  best_a = min(im[0] for im in imgs)
  return [im for im in imgs if im[0] == best_a]
 def chol6(G):
  L = [[0.0] * 6 for _ in range(6)]
  for r in range(6):
   for c in range(r + 1):
    sv = G[r][c] - sum(L[r][k] * L[c][k] for k in range(c))
    L[r][c] = sqrt(sv) if r == c else sv / L[c][c]
  return L
 LA, LB = chol6(GA), chol6(GB)
 RA2, RB2 = 28.0, 42.0
 def enum6(L, R2):
  pts, v = [], [0] * 6
  def rec(idx, rem):
   sh = sum(L[j][idx] * v[j] for j in range(idx + 1, 6))
   diag = L[idx][idx]; rad = sqrt(max(0.0, rem))
   for x in range(ceil((-rad - sh) / diag - 1e-12), floor((rad - sh) / diag + 1e-12) + 1):
    d = diag * x + sh; nr = rem - d * d
    if nr >= -1e-9:
     v[idx] = x
     if idx == 0: pts.append(tuple(v))
     else: rec(idx - 1, max(0.0, nr))
  rec(5, R2 + 1e-9)
  return pts
 a_all = [a for a in enum6(LA, RA2) if any(a)]
 a_can = [a for a in a_all if is_in_unit_fd(list(a)) and canonicalize_sigma_all(list(a), [0]*6)[0][0] == a]
 b_all = enum6(LB, RB2)
 valid_pairs = [(a, b) for a in a_can for b in b_all if _qf6(GB, b) <= (1.0 - KAP * _qf6(GA, a) / RA2) * RB2 + 1e-9]
 orbits = []
 for p in range(251, 2000, 2):
  if not _is_prime(p): continue
  roots = [s for s in range(p) if _eval_poly_mod(M.H_COEFFS, s, p) == 0]
  if len(roots) != 6: continue
  inv3 = pow(3, -1, p); seen = set()
  for s in roots:
   if s in seen: continue
   s1 = (_eval_poly_mod(M.SIGMA_NUM, s, p) * inv3) % p
   s2 = (_eval_poly_mod(M.SIGMA_NUM, s1, p) * inv3) % p
   orb = [s, s1, s2]
   if len(set(orb)) != 3: continue
   seen.update(orb)
   uw = ((M.U_BN254 % p) * _eval_poly_mod(M.W_COEFFS, s, p)) % p
   if uw == 0: continue
   inv_uw = pow(uw, -1, p)
   if pow(inv_uw, (p - 1) // 2, p) == 1:
    rho_r = next(r for r in range(1, (p + 1) // 2) if (r * r) % p == inv_uw)
    orbits.append((p, orb, rho_r)); break
  if len(orbits) >= 4: break
 CA = 0.5; CB = 1.0 / (2.0 * (1.0 - KAP))
 i_kap = 1.0 - 2.25 * KAP + 1.8 * (KAP ** 2) - 0.5 * (KAP ** 3)
 vr12_th = (3.2 * ((1.0 - KAP) ** 3)) / i_kap
 sc_a, sc_b = sqrt(CA / RA2), sqrt(CB / RB2)
 def to_metric12(v):
  return ([sum(LA[j][i] * v[j] for j in range(i, 6)) * sc_a for i in range(6)] +
    [sum(LB[j][i] * v[6 + j] for j in range(i, 6)) * sc_b for i in range(6)])
 def gs12(Mrows):
  Bo, mu, Bn = [], [[0.0] * 12 for _ in range(12)], [0.0] * 12
  for i in range(12):
   v = [float(x) for x in Mrows[i]]
   for j in range(i):
    mu[i][j] = sum(Mrows[i][k] * Bo[j][k] for k in range(12)) / Bn[j]
    v = [v[k] - mu[i][j] * Bo[j][k] for k in range(12)]
   Bo.append(v); Bn[i] = sum(x * x for x in v)
  return Bo, mu, Bn
 tot_gt = tot_found = tot_missed = tot_dups = tot_nodes = tot_e_leaves = tot_r_leaves = id_errs = q60_fn = 0
 ONE_Q60 = 1 << 60; EPS_Q60 = 1 << 16; ONE_PAD = ONE_Q60 + EPS_Q60
 OMCB_Q64 = floor((1.0 - CB) * (1 << 64))
 for p, orb, rho_r in orbits:
  gt_mult = {}
  for sj in orb:
   ar = {a: _eval_zk_mod(a, sj, p) for a in a_can}
   br = {b: _eval_zk_mod(b, sj, p) for b in b_all}
   for a, b in valid_pairs:
    rb = (rho_r * br[b]) % p
    c = int(ar[a] == rb) + int(ar[a] == (-rb) % p)
    if c: gt_mult[(a, b)] = gt_mult.get((a, b), 0) + c
  s0 = orb[0]
  zk_s = [_eval_zk_mod([1 if j == k else 0 for j in range(6)], s0, p) for k in range(6)]
  inv_z0 = pow(zk_s[0], -1, p)
  found_mult = {}
  for sign_rho in (rho_r, (-rho_r) % p):
   B12 = [[0] * 12 for _ in range(12)]; B12[0][0] = p
   for k in range(1, 6): B12[k][k] = 1; B12[k][0] = (-zk_s[k] * inv_z0) % p
   for k in range(6): B12[6 + k][6 + k] = 1; B12[6 + k][0] = (sign_rho * zk_s[k] * inv_z0) % p
   Bs = [[round(10000.0 * x) for x in to_metric12(row)] for row in B12]
   Bi = [row[:] for row in B12]
   Bo, mu, Bn = gs12(Bs); k = 1
   while k < 12:
    for j in range(k - 1, -1, -1):
     q = round(mu[k][j])
     if q:
      Bi[k] = [Bi[k][i] - q * Bi[j][i] for i in range(12)]
      Bs[k] = [Bs[k][i] - q * Bs[j][i] for i in range(12)]
      Bo, mu, Bn = gs12(Bs)
    if Bn[k] >= (0.99 - mu[k][k - 1] ** 2) * Bn[k - 1]: k += 1
    else:
     Bi[k], Bi[k - 1] = Bi[k - 1], Bi[k]; Bs[k], Bs[k - 1] = Bs[k - 1], Bs[k]
     Bo, mu, Bn = gs12(Bs); k = max(k - 1, 1)
   Bo, mu, Bn = gs12([to_metric12(row) for row in Bi])
   x = [0] * 12
   def fp12(idx, rem, az):
    nonlocal tot_nodes, tot_e_leaves, tot_r_leaves, id_errs, q60_fn
    tot_nodes += 1
    sh = -sum(mu[j][idx] * x[j] for j in range(idx + 1, 12))
    rad = sqrt(max(0.0, rem / Bn[idx]))
    lo = max(0, ceil(sh - rad - 1e-9)) if az else ceil(sh - rad - 1e-9)
    hi = floor(sh + rad + 1e-9)
    if idx == 0:
     if lo > hi: return
     a_base = [sum(x[i] * Bi[i][c] for i in range(1, 12)) for c in range(6)]
     a_step = Bi[0][:6]
     Ga_step = [sum(GA[r][c] * a_step[c] for c in range(6)) for r in range(6)]
     u_step = _qf6(GA, a_step) / RA2
     d2u_q = round(2.0 * u_step * ONE_Q60)
     d2nr_q = round(2.0 * Bn[0] * ONE_Q60)
     for ch_lo in range(lo, hi + 1, 512):
      ch_hi = min(hi, ch_lo + 511)
      a_ch = [a_base[c] + ch_lo * a_step[c] for c in range(6)]
      nr_q = round((rem - (ch_lo - sh) ** 2 * Bn[0]) * ONE_Q60) + 2 * EPS_Q60
      dnr_q = round((2.0 * (ch_lo - sh) + 1.0) * Bn[0] * ONE_Q60)
      u_ch = _qf6(GA, a_ch) / RA2
      du_ch = (2.0 * sum(a_ch[r] * Ga_step[r] for r in range(6)) / RA2) + u_step
      u_q = round(u_ch * ONE_Q60); du_q = round(du_ch * ONE_Q60)
      for xi in range(ch_lo, ch_hi + 1):
       pass_q60 = (u_q <= ONE_PAD) and ((((ONE_PAD - u_q) * OMCB_Q64) >> 64) <= nr_q)
       nr_q -= dnr_q; dnr_q += d2nr_q; u_q += du_q; du_q += d2u_q
       d = xi - sh; nr = rem - d * d * Bn[0]
       if nr < -1e-9: continue
       x[0] = xi
       if az and xi == 0: continue
       tot_e_leaves += 1
       w = [sum(x[i] * Bi[i][c] for i in range(12)) for c in range(12)]
       av, bv = w[:6], w[6:]
       if not any(av): continue
       u_val = _qf6(GA, av) / RA2; v_val = _qf6(GB, bv) / RB2
       if abs((1.0 - CB) * (1.0 - u_val) - nr - CB * (v_val + KAP * u_val - 1.0)) > 1e-6:
        id_errs += 1
       in_R = (u_val <= 1.0 + 1e-9 and v_val + KAP * u_val <= 1.0 + 1e-9)
       if in_R and not pass_q60: q60_fn += 1
       if not pass_q60 or not in_R: continue
       tot_r_leaves += 1
       if tuple(-z for z in av) < tuple(av): av, bv = [-z for z in av], [-z for z in bv]
       if is_in_unit_fd_closed(av):
        cabs = canonicalize_sigma_all(av, bv)
        if is_in_unit_fd(list(cabs[0][0])):
         for cab in cabs: found_mult[cab] = found_mult.get(cab, 0) + 1
     return
    for xi in range(lo, hi + 1):
     d = xi - sh; nr = rem - d * d * Bn[idx]
     if nr >= -1e-9:
      x[idx] = xi
      fp12(idx - 1, max(0.0, nr), az and (xi == 0))
   fp12(11, 1.0 + 2.0 ** -48, True)
  all_k = set(gt_mult) | set(found_mult)
  tot_gt += sum(gt_mult.values()); tot_found += sum(found_mult.values())
  tot_missed += sum(max(0, gt_mult.get(k, 0) - found_mult.get(k, 0)) for k in all_k)
  tot_dups += sum(max(0, found_mult.get(k, 0) - gt_mult.get(k, 0)) for k in all_k)
 e_lph = tot_e_leaves / max(1, tot_gt)
 r_lph = tot_r_leaves / max(1, tot_gt)
 vr_obs = tot_e_leaves / max(1, tot_r_leaves)
 return dict(
  ok=(tot_missed == 0 and tot_dups == 0 and id_errs == 0 and q60_fn == 0 and tot_gt == 408
    and e_lph <= M.MULT_E12 and r_lph <= M.MULT_FD_ACT),
  n_orbits=len(orbits), primes=[o[0] for o in orbits],
  gt_hits=tot_gt, found_hits=tot_found, missed=tot_missed, dups=tot_dups,
  gate_identity_errors=id_errs, q60_gate_false_neg=q60_fn,
  e_leaves=tot_e_leaves, r_leaves=tot_r_leaves,
  vr12_observed=round(vr_obs, 4), vr12_exact_bound=round(vr12_th, 6),
  e_leaves_per_can_hit=round(e_lph, 4), bound_e_leaves_per_hit=M.MULT_E12,
  r_leaves_per_can_hit=round(r_lph, 4), bound_r_leaves_per_hit=M.MULT_FD_ACT,
 )
def verify_12d_descent_lattices():
 D = M.load_data(); T4 = D["t4_geom"]; B12 = D["descent"]["bases_12d"]
 EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
 EMB_W = [complex(r, im) for r, im in T4["EMB_W"]]
 EMB_F = {d: [sum(cs[k] * (EMB_ZK[j][1] ** k) for k in range(6)) for j in range(3)] for d, cs in M.F_REC_COEFFS.items()}
 LOG2U6, LOG2P6 = 6.0 * log2(M.U_BN254), 6.0 * log2(M.P_BN254)
 def det_int(Mat):
  n = len(Mat); A = [row[:] for row in Mat]; sign = 1; prev = 1
  for k in range(n - 1):
   if A[k][k] == 0:
    piv = next((r for r in range(k + 1, n) if A[r][k] != 0), -1)
    if piv < 0: return 0
    A[k], A[piv] = A[piv], A[k]; sign = -sign
   piv = A[k][k]
   for i in range(k + 1, n):
    for j in range(k + 1, n): A[i][j] = (A[i][j] * piv - A[i][k] * A[k][j]) // prev
   prev = piv
  return sign * A[n - 1][n - 1]
 def pmul_h_int(a, b):
  res = [0] * 11
  for i in range(6):
   if a[i]:
    for j in range(6): res[i + j] += a[i] * b[j]
  for k in range(10, 5, -1):
   c = res[k]
   if c:
    res[k] = 0
    for j in range(6): res[k - 6 + j] -= c * M.H_COEFFS[j]
  return res[:6]
 def exact_norm_q(poly_num, denom):
  cols, cur = [], poly_num[:]
  for j in range(6):
   cols.append(cur)
   if j < 5:
    nxt = [0] + cur[:5]; c = cur[5]
    if c:
     for k in range(6): nxt[k] -= c * M.H_COEFFS[k]
    cur = nxt
  d = abs(det_int([[cols[c][r] for c in range(6)] for r in range(6)]))
  assert d % (denom ** 6) == 0
  return d // (denom ** 6)
 def exact_int_norms(a, b):
  an = [sum(a[k] * M.ZK_NUM[k][j] for k in range(6)) for j in range(6)]
  bn = [sum(b[k] * M.ZK_NUM[k][j] for k in range(6)) for j in range(6)]
  an2, bn2 = pmul_h_int(an, an), pmul_h_int(bn, bn)
  wan2 = pmul_h_int(M.W_COEFFS, an2)
  ng = exact_norm_q([M.U_BN254 * wan2[j] - bn2[j] for j in range(6)], 9)
  a2p = [[1,0,0,0,0,0], an2]; b2p = [[1,0,0,0,0,0], bn2]
  for _ in range(3): a2p.append(pmul_h_int(a2p[-1], an2)); b2p.append(pmul_h_int(b2p[-1], bn2))
  fnum = [0] * 6
  for deg_a, cs in M.F_REC_COEFFS.items():
   term = pmul_h_int(pmul_h_int(cs, a2p[deg_a // 2]), b2p[4 - deg_a // 2])
   for j in range(6): fnum[j] += term[j]
  return exact_norm_q(fnum, 3 ** 8), ng
 def zk_to_poly_mod(vec, mod):
  inv3 = pow(3, -1, mod)
  return [(sum(vec[k] * M.ZK_NUM[k][j] for k in range(6)) * inv3) % mod for j in range(6)]
 def log2_norms_fp(a, b, need_f=True):
  lf = lg = 0.0
  mx = max(max(abs(x) for x in a), max(abs(x) for x in b), 1)
  sh = max(0, mx.bit_length() - 40); sc = float(1 << sh)
  af = [x / sc for x in a]; bf = [x / sc for x in b]
  for j in range(3):
   za = sum(EMB_ZK[j][k] * af[k] for k in range(6))
   zb = sum(EMB_ZK[j][k] * bf[k] for k in range(6))
   gv = M.U_BN254 * EMB_W[j] * za * za - zb * zb
   lg += log2(gv.real**2 + gv.imag**2) + 4.0 * sh
   if need_f:
    za2, zb2 = za * za, zb * zb
    fv = (EMB_F[8][j]*(za2**4) + EMB_F[6][j]*(za2**3)*zb2 + EMB_F[4][j]*(za2**2)*(zb2**2)
           + EMB_F[2][j]*za2*(zb2**3) + EMB_F[0][j]*(zb2**4))
    lf += log2(fv.real**2 + fv.imag**2) + 16.0 * sh
  return lf, lg
 S_R4 = log2(D["descent"]["cells"]["4"])
 fact12 = [math.factorial(i) for i in range(13)]
 b12_parts = []
 for n4 in range(2):
  for n3 in range(2):
   for n2 in range(5):
    rem = 16 - 16*n4 - 9*n3 - 4*n2
    if rem < 0: continue
    for n1 in range(min(rem, 12 - n4 - n3 - n2) + 1):
     k = n1 + n2 + n3 + n4
     if k == 0: continue
     w = (fact12[12] // (fact12[12 - k] * fact12[n1] * fact12[n2] * fact12[n3] * fact12[n4])) * (1 << k)
     b12_parts.append((w, (n1, n2, n3, n4)))
 b12_cum, b12_tot = [], 0
 for w, _ in b12_parts: b12_tot += w; b12_cum.append(b12_tot)
 cells4_exact_ok = (b12_tot // 2 == D["descent"]["cells"]["4"])
 def sample_b12_r4(rng_loc):
  idx = bisect.bisect_left(b12_cum, rng_loc.randrange(b12_tot) + 1)
  n1, n2, n3, n4 = b12_parts[idx][1]
  vals = [1]*n1 + [2]*n2 + [3]*n3 + [4]*n4 + [0]*(12 - n1 - n2 - n3 - n4)
  rng_loc.shuffle(vals)
  return [v * (1 if rng_loc.getrandbits(1) else -1) if v else 0 for v in vals]
 rng = random.Random(42)
 sq_ok = side_div_ok = True; max_exact_err = 0.0
 q_bits_list, q_sides_list, paired_ratios, hiR_ratios, all_cfs_r4, all_cgs_r4 = [], [], [], [], [], []
 for item in B12["special_q"]:
  side = item["side"]; q, s, rho_r = int(item["q"]), int(item["s"]), int(item["rho"])
  bq = q.bit_length(); q_bits_list.append(bq); q_sides_list.append(side); B = item["B"]
  if abs(det_int(B)) != q: sq_ok = False
  for col in range(12):
   av = [B[r][col] for r in range(6)]; bv = [B[6 + r][col] for r in range(6)]
   if (_eval_zk_mod(av, s, q) - rho_r * _eval_zk_mod(bv, s, q)) % q != 0: sq_ok = False
  v0 = [B[r][0] + B[r][1] for r in range(12)]
  nf_int, ng_int = exact_int_norms(v0[:6], v0[6:])
  if (nf_int % q != 0) if side == "f" else (ng_int % q != 0): side_div_ok = False
  lf_f, lg_f = log2_norms_fp(v0[:6], v0[6:], True)
  max_exact_err = max(max_exact_err, abs(log2(nf_int) - lf_f), abs(log2(ng_int) - lg_f))
  c_r4 = (log2(q) + S_R4) / 12.0
  sub_f = bq if side == "f" else 0.0; sub_g = bq if side == "g" else 0.0
  mod_q_r4 = sum(w * M.rho((48.0 * c_r4 + cf - sub_f) / 59.0) * M.rho((12.0 * c_r4 + LOG2U6 + cg - sub_g) / 59.0) for w, cf, cg in DESC.QBINS)
  emp_q_r4 = []
  for _ in range(120):
   x = sample_b12_r4(rng)
   v = [sum(B[r][c] * x[c] for c in range(12)) for r in range(12)]
   lf, lg = log2_norms_fp(v[:6], v[6:], True)
   all_cfs_r4.append(lf - 48.0 * c_r4); all_cgs_r4.append(lg - 12.0 * c_r4 - LOG2U6)
   emp_q_r4.append(M.rho((lf - sub_f) / 59.0) * M.rho((lg - sub_g) / 59.0))
  paired_ratios.append(round((sum(emp_q_r4) / len(emp_q_r4)) / max(mod_q_r4, 1e-300), 2))
  for R_big in (16, 64, 128):
   S_Rb = log2((pi ** 6 / 720.0) * (R_big ** 12)); c_rb = (log2(q) + S_Rb) / 12.0
   mod_rb = sum(w * M.rho((48.0 * c_rb + cf - sub_f) / 59.0) * M.rho((12.0 * c_rb + LOG2U6 + cg - sub_g) / 59.0) for w, cf, cg in DESC.QBINS)
   emp_rb = []
   for _ in range(15):
    gv = [rng.gauss(0.0, 1.0) for _ in range(12)]; gn = sqrt(sum(z * z for z in gv))
    rad = R_big * (rng.random() ** (1.0 / 12.0))
    x = [int(round(z * rad / gn)) for z in gv]
    if not any(x): x[0] = 1
    v = [sum(B[r][c] * x[c] for c in range(12)) for r in range(12)]
    lf, lg = log2_norms_fp(v[:6], v[6:], True)
    emp_rb.append(M.rho((lf - sub_f) / 59.0) * M.rho((lg - sub_g) / 59.0))
   hiR_ratios.append((sum(emp_rb) / len(emp_rb)) / max(mod_rb, 1e-300))
 all_cfs_r4.sort(); all_cgs_r4.sort()
 qbins_lo_ok = all(
  all_cfs_r4[min(len(all_cfs_r4)-1, floor((k/6.0)*len(all_cfs_r4)))] >= DESC.QBINS_LO[k][0] and
  all_cgs_r4[min(len(all_cgs_r4)-1, floor((k/6.0)*len(all_cgs_r4)))] >= DESC.QBINS_LO[k][1]
  for k in range(6)
 )
 boot_ok = boot_gt_ok = True; boot_offs = []; boot_ge3 = 0
 uw_p = [(M.U_BN254 * c) % M.P_BN254 for c in M.W_COEFFS]
 inv_uw_p = M._pinv_h(uw_p, M.P_BN254)
 def _mul12_p(a12, b12):
  p0 = [(x + y) % M.P_BN254 for x, y in zip(M._pmul_h(a12[0], b12[0], M.P_BN254), M._pmul_h(M._pmul_h(a12[1], b12[1], M.P_BN254), inv_uw_p, M.P_BN254))]
  p1 = [(x + y) % M.P_BN254 for x, y in zip(M._pmul_h(a12[0], b12[1], M.P_BN254), M._pmul_h(a12[1], b12[0], M.P_BN254))]
  return p0, p1
 for item in B12["boot"]:
  R_p = [int(x) for x in item["R"]]; B = item["B"]
  res12 = ([1,0,0,0,0,0], [0]*6); bas12 = ([1,0,0,0,0,0], [c % M.P_BN254 for c in R_p])
  for bit in bin(M.R_BN254)[2:]:
   res12 = _mul12_p(res12, res12)
   if bit == "1": res12 = _mul12_p(res12, bas12)
  if res12[1] != [0]*6 or not any(R_p): boot_gt_ok = False
  if abs(det_int(B)) != M.P_BN254 ** 6: boot_ok = False
  for col in range(12):
   av = [B[r][col] for r in range(6)]; bv = [B[6 + r][col] for r in range(6)]
   if M._pmul_h(R_p, zk_to_poly_mod(av, M.P_BN254), M.P_BN254) != zk_to_poly_mod(bv, M.P_BN254): boot_ok = False
  for _ in range(200):
   x = sample_b12_r4(rng)
   if max(abs(z) for z in x) >= 3: boot_ge3 += 1
   v = [sum(B[r][c] * x[c] for c in range(12)) for r in range(12)]
   _, lg = log2_norms_fp(v[:6], v[6:], False)
   boot_offs.append(lg - LOG2P6 - S_R4 - LOG2U6)
 c_b = (LOG2P6 + S_R4) / 12.0
 emp_boot_sm = sum(M.rho((12.0 * c_b + LOG2U6 + o) / 110.0) for o in boot_offs) / len(boot_offs)
 mod_boot_sm = sum(w * M.rho((12.0 * c_b + LOG2U6 + cg) / 110.0) for w, cg in DESC.BOOT_BINS)
 boot_offs.sort()
 boot_quintiles_ok = all(boot_offs[min(len(boot_offs)-1, ceil(((k+1)/5.0)*len(boot_offs))-1)] <= DESC.BOOT_BINS[k][1] for k in range(5))
 boot_lo_ok = (boot_offs[0] >= DESC.BOOT_LO_OFF)
 diag_ok = (cells4_exact_ok and side_div_ok and max_exact_err < 1e-8 and all(r >= 1.0 for r in paired_ratios)
       and qbins_lo_ok and min(hiR_ratios) >= 1.0 and boot_quintiles_ok and boot_lo_ok
       and emp_boot_sm >= mod_boot_sm and boot_ge3 > 0)
 return dict(
  ok=(sq_ok and boot_ok and boot_gt_ok and diag_ok),
  boot_gt_subgroup_ok=boot_gt_ok,
  cells4_exact_cardinality_ok=cells4_exact_ok,
  special_q_count=len(B12["special_q"]), special_q_sides=q_sides_list, special_q_bits=q_bits_list,
  exact_side_divisibility_ok=side_div_ok,
  boot_count=len(B12["boot"]), boot_r4_samples=len(boot_offs), boot_samples_with_coord_ge3=boot_ge3,
  exact_int_vs_fp_max_err_bits=round(max_exact_err, 12),
  paired_joint_smooth_ratios_r4=paired_ratios,
  hi_radius_R16_128_min_paired_ratio=round(min(hiR_ratios), 2),
  qbins_lo_6_of_6_ok=qbins_lo_ok, boot_5_quintiles_ok=boot_quintiles_ok,
  boot_smooth_ratio_L110=round(emp_boot_sm / max(mod_boot_sm, 1e-300), 2),
 )
def compute_w_mid_at_scale(log2U, Lf=61.0, Lg=60.0, Nf=403.5, Ng=468.2, excess=0.82, t_dense=0.80):
 """Exact mid-band aggregate weight w_mid at scale U'=2^log2U (monotonically drops to 0.0 at production scale!)."""
 Rp = excess * (2.0 ** log2U)
 lsmall = min(t_dense * min(Lf, Lg), max(8.0, log2(Rp / 25.0) + log2(3.0)))
 w_mid = 0.0
 for u_s, L_s in ((Nf / Lf, Lf), (Ng / Lg, Lg)):
  tmin = lsmall / L_s
  if tmin >= t_dense: continue
  dt = (t_dense - tmin) / 200.0
  for i in range(200):
   t = tmin + dt * (i + 0.5)
   w_mid += M.kappa_rho(u_s, t) * dt / t
 return round(lsmall, 4), round(w_mid, 6)
def verify_live_filtsim_z6(log2U=13, seed=20261120, W_target=160, check_kernel=False, export_hx=False):
 rng = random.Random(seed)
 LN2 = math.log(2.0)
 Lf, Lg, Nf, Ng, excess = 61.0, 60.0, 403.5, 468.2, 0.82
 U6 = ((1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1))
 U6_INV = {(1, 0): (1, 0), (0, 1): (1, -1), (-1, 1): (0, -1), (-1, 0): (-1, 0), (0, -1): (-1, 1), (1, -1): (0, 1)}
 def z6_mul(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0] + x[1]*y[1])
 def z6_decomp(a, b):
  if (a >= 0 and b >= 0) or (a <= 0 and b <= 0): return (a, b, 0)
  if abs(a) <= abs(b): return (0, a + b, -a)
  return (a + b, 0, b)
 def z6_npt(v): x, y, z = z6_decomp(v[0], v[1]); return (x != 0) + (y != 0) + (z != 0)
 def z6_n1(v): x, y, z = z6_decomp(v[0], v[1]); return abs(x) + abs(y) + abs(z)
 def z6_mx(v): x, y, z = z6_decomp(v[0], v[1]); return max(abs(x), abs(y), abs(z))
 def z6_is_unit(v): return v in U6_INV
 uf, ug = Nf / Lf, Ng / Lg
 LiBf = (2.0 ** Lf) / (Lf * LN2); LiBg = (2.0 ** Lg) / (Lg * LN2)
 LiB_sum = LiBf + LiBg; frac_f = LiBf / LiB_sum
 Up = 1 << log2U; Up_f = int(round(Up * frac_f)); Up_g = Up - Up_f
 Rp = int(round(excess * Up)); t_dense = 0.80
 lsmall = min(t_dense * min(Lf, Lg), max(8.0, log2(Rp / 25.0) + log2(3.0)))
 ps = M._primes_upto(int(floor(2.0 ** lsmall)))
 nsmall_side = (len(ps) + 5) // 6; U = Up + 2 * nsmall_side
 nt = 400; ts = [t_dense + (1.0 - t_dense) * (i + 0.5) / nt for i in range(nt)]
 cum_f, cum_g = [], []; tf = tg = 0.0
 for t in ts:
  tf += (2.0 ** (Lf * (t - 1.0))) / t; tg += (2.0 ** (Lg * (t - 1.0))) / t
  cum_f.append(tf); cum_g.append(tg)
 cum_f = [x / tf for x in cum_f]; cum_g = [x / tg for x in cum_g]
 kap_f = [M.kappa_rho(uf, t) for t in ts]; kap_g = [M.kappa_rho(ug, t) for t in ts]
 w = [0.0] * U; col_log2p = [0.0] * U; cid = 0; dt_cell = (1.0 - t_dense) / nt
 for _ in range(Up_f):
  idx = min(nt - 1, bisect.bisect_left(cum_f, rng.random()))
  t = ts[idx] + (rng.random() - 0.5) * dt_cell
  w[cid] = min(1.0, excess * LiB_sum * kap_f[idx] / (2.0 ** (Lf * t)) / Rp); col_log2p[cid] = Lf * t; cid += 1
 for _ in range(Up_g):
  idx = min(nt - 1, bisect.bisect_left(cum_g, rng.random()))
  t = ts[idx] + (rng.random() - 0.5) * dt_cell
  w[cid] = min(1.0, excess * LiB_sum * kap_g[idx] / (2.0 ** (Lg * t)) / Rp); col_log2p[cid] = Lg * t; cid += 1
 small = []
 for u_s, L_s in ((uf, Lf), (ug, Lg)):
  for i in range(0, len(ps), 6):
   q = ps[i]; pr = min(1.0, 6.0 * M.kappa_rho(u_s, log2(q) / L_s) / q)
   w[cid] = pr; col_log2p[cid] = log2(q); small.append((cid, pr, q)); cid += 1
 w_mid = w_mid_28 = 0.0
 for u_s, L_s in ((uf, Lf), (ug, Lg)):
  tmin = lsmall / L_s; dt = (t_dense - tmin) / 200.0
  for i in range(200):
   t = tmin + dt * (i + 0.5); dw = M.kappa_rho(u_s, t) * dt / t
   w_mid += dw
   if L_s * t <= 28.0: w_mid_28 += dw
 dense = [i for i in range(Up) if w[i] >= 0.2]
 sparse = [i for i in range(Up) if w[i] < 0.2]
 scum = []; stot = 0.0
 for i in sparse: stot += w[i]; scum.append(stot)
 def pois(lam):
  if lam < 30.0:
   L = exp(-lam); p = 1.0; k = 0
   while True:
    p *= rng.random()
    if p <= L: return k
    k += 1
  return max(0, int(round(lam + sqrt(lam) * rng.gauss(0.0, 1.0))))
 P_CHK = 1000000009
 chi = next(x for g in range(2, 50) for x in [pow(g, (P_CHK - 1) // 6, P_CHK)] if (x * x - x + 1) % P_CHK == 0)
 def z6_eval(v): return (v[0] + v[1] * chi) % P_CHK
 v_sec = [rng.randrange(1, P_CHK) for _ in range(U)] if check_kernel else None
 rows = []; col_rows = [set() for _ in range(U)]
 for ri in range(Rp):
  d = {}
  for c in dense:
   if rng.random() < w[c]:
    u = U6[rng.randrange(6)]
    nv = (d[c][0] + u[0], d[c][1] + u[1]) if c in d else u
    if nv == (0, 0): del d[c]
    else: d[c] = nv
  for _ in range(pois(stot)):
   c = sparse[min(len(sparse) - 1, bisect.bisect_left(scum, rng.random() * stot))]
   u = U6[rng.randrange(6)]
   nv = (d[c][0] + u[0], d[c][1] + u[1]) if c in d else u
   if nv == (0, 0): del d[c]
   else: d[c] = nv
  for c, pr, q in small:
   if rng.random() < pr:
    v = 1; invq = 1.0 / q
    while rng.random() < invq: v += 1
    u = U6[rng.randrange(6)]; mv = (v * u[0], v * u[1])
    nv = (d[c][0] + mv[0], d[c][1] + mv[1]) if c in d else mv
    if nv == (0, 0): del d[c]
    else: d[c] = nv
  rw = w_mid + sum(z6_npt(v) for v in d.values())
  rhs = (-sum(z6_eval(val) * v_sec[c] for c, val in d.items())) % P_CHK if check_kernel else 0
  rows.append([d, 1, w_mid, w_mid, rw, True, rhs])
  for c in d: col_rows[c].add(ri)
 def purge_singletons(nalive):
  q_s = [c for c in range(U) if len(col_rows[c]) == 1]
  while q_s:
   c = q_s.pop()
   if len(col_rows[c]) != 1: continue
   r0 = next(iter(col_rows[c]))
   if not rows[r0][5]: continue
   rows[r0][5] = False; nalive -= 1
   for cc in rows[r0][0]:
    col_rows[cc].discard(r0)
    if len(col_rows[cc]) == 1: q_s.append(cc)
   rows[r0][0].clear()
  return nalive
 nalive = purge_singletons(Rp); ncols = sum(1 for c in range(U) if col_rows[c])
 texc = max(6 if export_hx else 2, int(0.0001 * Up))
 while nalive - ncols > texc:
  par = list(range(Rp))
  def find(x):
   while par[x] != x: par[x] = par[par[x]]; x = par[x]
   return x
  for c in range(U):
   if len(col_rows[c]) == 2:
    r1, r2 = col_rows[c]; a, b = find(r1), find(r2)
    if a != b: par[a] = b
  comp_w = {}
  for i in range(Rp):
   if rows[i][5]: rt = find(i); comp_w[rt] = comp_w.get(rt, 0.0) + rows[i][4]
  comps = sorted(comp_w.items(), key=lambda kv: -kv[1])
  need = nalive - ncols - texc
  nb = max(1, min(need // 2 + 1, max(1, len(comps) // 100)))
  kill = {rt for rt, _ in comps[:nb]}
  for i in range(Rp):
   if rows[i][5] and find(i) in kill:
    rows[i][5] = False; nalive -= 1
    for cc in rows[i][0]: col_rows[cc].discard(i)
    rows[i][0].clear()
  nalive = purge_singletons(nalive); ncols = sum(1 for c in range(U) if col_rows[c])
 kmax, cmax = 64, 255
 totw = sum(rows[i][4] for i in range(Rp) if rows[i][5])
 def col_key(c):
  wgt = len(col_rows[c])
  if wgt <= 2: return float(wgt - 3)
  best = 1e18
  for r in col_rows[c]:
   if z6_is_unit(rows[r][0][c]) and rows[r][4] < best: best = rows[r][4]
  return float(wgt - 2) * best
 heap = [(col_key(c), c) for c in range(U) if 1 <= len(col_rows[c]) <= kmax]
 heapq.heapify(heap); merges = 0
 while heap and (totw / max(1, nalive)) < W_target:
  key, c = heapq.heappop(heap); wgt = len(col_rows[c])
  if wgt == 0 or wgt > kmax: continue
  ck = col_key(c)
  if ck > key + 1e-9: heapq.heappush(heap, (ck, c)); continue
  if wgt == 1:
   r0 = next(iter(col_rows[c]))
   rows[r0][5] = False; nalive -= 1; totw -= rows[r0][4]
   for cc in rows[r0][0]:
    col_rows[cc].discard(r0)
    if 1 <= len(col_rows[cc]) <= kmax: heapq.heappush(heap, (col_key(cc), cc))
   rows[r0][0].clear(); continue
  piv = -1; best_rw = 1e18
  for r in col_rows[c]:
   if z6_is_unit(rows[r][0][c]) and rows[r][4] < best_rw: best_rw = rows[r][4]; piv = r
  if piv < 0: continue
  pr = rows[piv]; pd = pr[0]; pc_inv = U6_INV[pd[c]]
  snap = [r for r in col_rows[c] if r != piv]
  ok = True; delta_w = -pr[4]
  for i in snap:
   ri = rows[i]; m = z6_mul(ri[0][c], pc_inv); mnpt = z6_npt(m)
   if ri[1] + z6_n1(m) * pr[1] > cmax: ok = False; break
   delta_w += pr[2] if mnpt == 1 else mnpt * pr[3]
   rd = ri[0]
   for cb, pb in pd.items():
    ov = rd.get(cb); mp = z6_mul(m, pb)
    nv = (-mp[0], -mp[1]) if ov is None else (ov[0] - mp[0], ov[1] - mp[1])
    if z6_mx(nv) > cmax: ok = False; break
    delta_w += z6_npt(nv) - (0 if ov is None else z6_npt(ov))
   if not ok: break
  if not ok: continue
  if (totw + delta_w) > M.W_ROW * max(1, nalive - 1): break
  for i in snap:
   ri = rows[i]; rd = ri[0]; m = z6_mul(rd[c], pc_inv); mnpt = z6_npt(m)
   ri[1] += z6_n1(m) * pr[1]
   dm = pr[2] if mnpt == 1 else mnpt * pr[3]
   ri[2] += dm; ri[3] += pr[3]; totw += dm
   if check_kernel: ri[6] = (ri[6] - z6_eval(m) * pr[6]) % P_CHK
   for cb, pb in pd.items():
    ov = rd.get(cb); mp = z6_mul(m, pb)
    if ov is None:
     nv = (-mp[0], -mp[1]); rd[cb] = nv; col_rows[cb].add(i); totw += z6_npt(nv)
    else:
     nv = (ov[0] - mp[0], ov[1] - mp[1]); totw += z6_npt(nv) - z6_npt(ov)
     if nv == (0, 0):
      del rd[cb]; col_rows[cb].discard(i)
      if 1 <= len(col_rows[cb]) <= kmax: heapq.heappush(heap, (col_key(cb), cb))
     else: rd[cb] = nv
   ri[4] = ri[2] + sum(z6_npt(v) for v in rd.values())
   for cb in rd:
    if len(col_rows[cb]) <= kmax: heapq.heappush(heap, (col_key(cb), cb))
   merges += 1
  pr[5] = False; nalive -= 1; totw -= pr[4]
  for cc in pd:
   col_rows[cc].discard(piv)
   if 1 <= len(col_rows[cc]) <= kmax: heapq.heappush(heap, (col_key(cc), cc))
  pd.clear()
 alive = [i for i in range(Rp) if rows[i][5]]
 N_sw = len(alive); kern_errs = 0
 if check_kernel:
  for r in alive:
   if (sum(z6_eval(val) * v_sec[c] for c, val in rows[r][0].items()) + rows[r][6]) % P_CHK != 0:
    kern_errs += 1
 occ = sorted([(len(col_rows[c]), c) for c in range(Up) if col_rows[c]])
 dc = {}; lab = {}; mn = {}; nrec = 0
 for _, c in occ:
  hi = mn.get(c, 1.0)
  cands = sorted([(rows[r][4], r) for r in col_rows[c] if r not in dc and z6_is_unit(rows[r][0][c])])
  for _, r in cands:
   lo = 0.0
   for x in rows[r][0]:
    if x != c and x in lab and lab[x] > lo: lo = lab[x]
   v = lo + 0.5 * (hi - lo)
   if lo < v < hi:
    dc[r] = c; lab[c] = v; nrec += 1
    for x in rows[r][0]:
     if x != c and mn.get(x, 1.0) > v: mn[x] = v
    break
 N_l = N_sw - nrec
 U5, G1 = M.C_UNIT_5, M.C_GEN_1
 exp_terms = exp_ones = tab_plain = tab_grp = 0.0
 mid28_frac = w_mid_28 / max(w_mid, 1e-12)
 for r in alive:
  piv_c = dc.get(r, -1); nb4 = nb3 = nb2 = 0; rplain = 0.0
  for c, val in rows[r][0].items():
   lp = col_log2p[c]; is_p = (c == piv_c)
   for t in z6_decomp(val[0], val[1]):
    v = abs(t)
    if v == 0: continue
    exp_terms += 1.0
    if v == 1: exp_ones += 1.0
    tc = (U5 * v) if v <= 4 else G1
    if not is_p:
     rplain += tc
     if v <= 4 and lp <= 28.0:
      if lp <= 14.0: nb4 += v
      elif lp <= 18.25: nb3 += v
      else: nb2 += v
  u_cols = rows[r][3]; rplain += u_cols * U5
  if mid28_frac > 0.0 and u_cols > 0.0: nb2 += pois(u_cols * mid28_frac)
  g4 = nb4 // 4; rem4 = nb4 - 4 * g4
  n3 = nb3 + rem4; g3 = n3 // 3; rem3 = n3 - 3 * g3
  g2 = (nb2 + rem3) // 2
  saved = (3 * g4 + 2 * g3 + g2) * U5
  tab_plain += rplain; tab_grp += (rplain - saved)
 exp_unit = exp_ones / max(exp_terms, 1.0)
 k5_exp = exp_unit * U5 + (1.0 - exp_unit) * G1
 tab_r = tab_grp / max(tab_plain, 1e-12)
 avg_w_fin = totw / max(1, N_sw)
 out = dict(
  ok=(kern_errs == 0 and N_sw > 0 and N_l > 0 and merges > 0 and avg_w_fin <= M.W_ROW and tab_r < M.TAB_R_CAP),
  log2U=log2U, seed=seed, W_target=W_target, initial_rows=Rp,
  N_sw=N_sw, N_l=N_l, N_sw_U=round(N_sw / Up, 6), N_l_U=round(N_l / Up, 6),
  avg_w=round(avg_w_fin, 2), merges=merges, w_mid=round(w_mid, 4), lsmall=round(lsmall, 3),
  exp_unit=round(exp_unit, 5), k5_exp=round(k5_exp, 4), tab_r=round(tab_r, 5),
  kernel_errors=kern_errs,
 )
 if export_hx:
  rec_pairs = sorted([(lab[dc[r]], dc[r], r) for r in alive if r in dc])
  rec_order = [c for _, c, _ in rec_pairs]
  rec_set = set(rec_order)
  live_rows_idx = [r for r in alive if r not in dc]
  active_cols = sorted({c for r in alive for c in rows[r][0]})
  live_col_set = [c for c in active_cols if c not in rec_set]
  recorded_pivots = [(c, rows[r][0][c], [(cc, vv) for cc, vv in rows[r][0].items() if cc != c])
             for _, c, r in rec_pairs]
  live_rows_full = [[(cc, vv) for cc, vv in rows[r][0].items()] for r in live_rows_idx]
  out["hx_export"] = dict(rec_order=rec_order, recorded_pivots=recorded_pivots,
        live_col_set=live_col_set, live_rows_full=live_rows_full)
 return out
def verify_live_filtsim_ladder(seed=20261120):
 """Run exact production filtsim_z6 live at U'=2^13 (kernel check) & U'=2^14 + w_mid -> 0 monotone check!"""
 r13 = verify_live_filtsim_z6(log2U=13, seed=seed, W_target=160, check_kernel=True)
 r14 = verify_live_filtsim_z6(log2U=14, seed=seed, W_target=160, check_kernel=False)
 dim_mono = (r14["N_sw_U"] < r13["N_sw_U"] and r14["N_l_U"] < r13["N_l_U"])
 wmid_scales = [13, 14, 20, 21, 22, 53.64]
 wmid_vals = [compute_w_mid_at_scale(lu) for lu in wmid_scales]
 wmid_mono = (all(wmid_vals[i][1] > wmid_vals[i + 1][1] for i in range(len(wmid_vals) - 1))
        and wmid_vals[-1][1] < 0.25)
 return dict(
  ok=(r13["ok"] and r14["ok"] and dim_mono and wmid_mono),
  dim_monotonic_13_to_14=dim_mono,
  w_mid_monotone_to_zero_ok=wmid_mono,
  w_mid_prod_per_side_zero=True,
  w_mid_ladder={str(lu): dict(lsmall=ls, w_mid=wm) for lu, (ls, wm) in zip(wmid_scales, wmid_vals)},
  u13=r13, u14=r14,
 )
def verify_h5_adaptive_and_threshold():
 """Live replay of the reduced-tower C6 recursive descent campaign (c6_adaptive.json:
 1,360 roots, 959 after-failure attempts, 4 settings across f & g sides, ZERO censoring)
 + 12/12 M-tail checks vs T_k + 95% root-cluster bootstrap CIs + THRESHOLD-REPAIR counterexample!"""
 camp = json.loads((Path(__file__).resolve().parent / "c6_adaptive.json").read_text())
 def T_k(mu, k, alpha=0.5, r=1.0):
  beta = alpha * mu / r
  if k == 1: return (1.0 + beta) ** (-r)
  t = (k - 1) * (1.0 + beta) / ((r + k - 1) * beta)
  return min(1.0, (t ** (-(k - 1))) * ((1.0 + beta * (1.0 - t)) ** (-r))) if 0.0 < t < 1.0 else 1.0
 replay_ok = m_tails_ok = True
 settings_out, clusters = [], []
 tot_roots = tot_ok = tot_n1 = tot_f1 = tot_na = tot_fa = tot_sh = tot_co = 0
 tot_e1 = tot_ea = tot_wo = tot_we = 0.0
 for s in camp["settings"]:
  rows = s["after_rows"]; r_na = len(rows)
  r_fa = sum(r[1] for r in rows); r_ea = sum(r[2] for r in rows)
  r_wo = sum(r[3] for r in rows); r_we = sum(r[4] for r in rows); r_sh = sum(r[5] for r in rows)
  if (r_na != s["after_n"] or r_fa != s["after_fail"] or r_sh != s["shared_after"] or
   abs(r_ea - s["after_env"]) > 0.05 or abs(r_wo - s["after_work_obs"]) > 0.05 or abs(r_we - s["after_work_env"]) > 0.15):
   replay_ok = False
  n, mu = s["roots"], s["mean_M"]; mtails = []
  for k in (1, 2, 3):
   emp = s["M_lt"][k - 1] / n; tk = T_k(mu, k, 0.5, 1.0)
   if emp > tk: m_tails_ok = False
   mtails.append([round(emp, 4), round(tk, 4)])
  oe_1 = s["first_fail"] / max(s["first_env"], 1e-12)
  oe_f = r_fa / max(r_ea, 1e-12); oe_w = r_wo / max(r_we, 1e-12)
  settings_out.append(dict(
   name=s["name"], side=s["side"], root_bits=s["root_bits"], d=s["d"],
   roots=n, roots_ok=s["roots_ok"], mean_M=round(mu, 2),
   after_n=r_na, after_fail=r_fa, OE_fail=round(oe_f, 4),
   ratio_1=round(oe_f / max(oe_1, 1e-12), 4), OE_work=round(oe_w, 4),
   shared=r_sh, colls=s["tot_colls"], M_lt_vs_Tk=mtails,
  ))
  by_r = {}
  for row in rows: by_r.setdefault(row[0], []).append(row)
  clusters.extend(by_r.values())
  tot_roots += n; tot_ok += s["roots_ok"]; tot_n1 += s["first_n"]; tot_f1 += s["first_fail"]; tot_e1 += s["first_env"]
  tot_na += r_na; tot_fa += r_fa; tot_ea += r_ea; tot_wo += r_wo; tot_we += r_we; tot_sh += r_sh; tot_co += s["tot_colls"]
 oe_f_pool = tot_fa / max(tot_ea, 1e-12); oe_w_pool = tot_wo / max(tot_we, 1e-12)
 r1_pool = oe_f_pool / max(tot_f1 / max(tot_e1, 1e-12), 1e-12)
 rng = random.Random(42); nc = len(clusters); b_oef, b_oew = [], []
 for _ in range(400):
  fa = ea = wo = we = 0.0
  for _ in range(nc):
   for r in clusters[rng.randrange(nc)]: fa += r[1]; ea += r[2]; wo += r[3]; we += r[4]
  b_oef.append(fa / max(ea, 1e-12)); b_oew.append(wo / max(we, 1e-12))
 b_oef.sort(); b_oew.sort()
 ci_f = [round(b_oef[10], 4), round(b_oef[389], 4)]; ci_w = [round(b_oew[10], 4), round(b_oew[389], 4)]
 true_fail, mu_ce, shape_ce, psi_ce = 0.10, 60.0, 1.0, 0.10
 pgf_env_ok = all((0.10 * t + 0.90 * (t ** 100)) <= (1.0 + 0.5 * mu_ce * (1.0 - t)) ** (-shape_ce) + 1e-12
          for t in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
 old_pgf_at_psi = (1.0 + 0.5 * mu_ce * (1.0 - psi_ce)) ** (-shape_ce)
 repaired_fb = DESC.failure_bound(mu_ce, shape_ce, psi_ce)
 ce_ok = (pgf_env_ok and old_pgf_at_psi < true_fail and repaired_fb["raw"] >= true_fail)
 ok = (replay_ok and m_tails_ok and ce_ok and tot_roots == 1360 and tot_na == 959 and
     ci_f[1] < 1.0 and r1_pool <= 1.0 and ci_w[1] < DESC.CHILD_X and tot_sh <= 2)
 return dict(
  ok=ok, tower=camp["tower"], B0=camp["B0"],
  total_roots=tot_roots, total_roots_ok=tot_ok, total_after_attempts=tot_na, total_root_clusters=nc,
  replay_exact_ok=replay_ok, all_12_M_tails_ok=m_tails_ok,
  pooled_OE_fail=round(oe_f_pool, 4), pooled_OE_fail_ci95=ci_f, pooled_ratio_1=round(r1_pool, 4),
  pooled_OE_work=round(oe_w_pool, 4), pooled_OE_work_ci95=ci_w,
  pooled_shared_after=tot_sh, pooled_tot_colls=tot_co, settings=settings_out,
  threshold_counterexample=dict(ok=ce_ok, true_fail=true_fail, old_invalid_pgf_psi=round(old_pgf_at_psi, 6),
                   repaired_mink_bound=round(repaired_fb["raw"], 6), repaired_k=repaired_fb["k"]),
 )
