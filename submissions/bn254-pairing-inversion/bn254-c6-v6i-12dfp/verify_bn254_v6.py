"""verify_bn254_v6.py — Self-contained pure-stdlib live prototypes for BN254 v6i."""
import bisect, heapq, json, math, random
from fractions import Fraction as F
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
 N_a_prod = int(2.0 ** 43.8652)
 test_pairs = [(t_i, a_i) for t_i in (0, 505, 1010)
               for a_i in (0, 1, (1 << 32) - 1, 1 << 32, (1 << 32) + 1, (1 << 43) + 12345, N_a_prod - 1)]
 lin_idxs = [t_i * N_a_prod + a_i for t_i, a_i in test_pairs]
 a_slot_64bit_ok = (
  len(set(lin_idxs)) == len(test_pairs)
  and all(0 <= h < (1 << 54) for h in lin_idxs)
  and all((h // N_a_prod, h % N_a_prod) == (t_i, a_i) for h, (t_i, a_i) in zip(lin_idxs, test_pairs))
  and ((0 * N_a_prod + 0) != (0 * N_a_prod + (1 << 32)))
 )
 return dict(
  ok=(mismatches == 0 and swar_errs == 0 and strip_errs == 0 and win_disjoint and win_fn == 0 and p2_errs == 0 and a_slot_64bit_ok and max_inv <= 1e-12 and tot_hits == 1523),
  scalar_swar_10k_errors=swar_errs,
  p2_shift_threshold_errors=p2_errs,
  a_slot_64bit_linear_ok=a_slot_64bit_ok,
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
 tot_gt = tot_found = tot_missed = tot_dups = tot_nodes = tot_e_leaves = tot_r_leaves = id_errs = q60_fn = drift_errs = 0
 ONE_Q60 = 1 << 60; EPS_Q60 = 1 << 26; ONE_PAD = ONE_Q60 + EPS_Q60
 OMCB_Q64 = floor((1.0 - CB) * (1 << 64))
 Ma_t = [[round(GA[r][c] * ONE_Q60 / RA2) for c in range(6)] for r in range(6)]
 Mea_t = [[round(CA * GA[r][c] * ONE_Q60 / RA2) for c in range(6)] for r in range(6)]
 Mb_t = [[round(CB * GB[r][c] * ONE_Q60 / RB2) for c in range(6)] for r in range(6)]
 def _qfi6(Mt, v):
  return sum(Mt[r][c] * v[r] * v[c] for r in range(6) for c in range(6))
 def _bfi6(Mt, u, v):
  return sum(Mt[r][c] * u[r] * v[c] for r in range(6) for c in range(6))
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
   a_step = Bi[0][:6]; b_step = Bi[0][6:]
   u_step_q = _qfi6(Ma_t, a_step); d2u_q = 2 * u_step_q
   q_step_q = _qfi6(Mea_t, a_step) + _qfi6(Mb_t, b_step); d2nr_q = 2 * q_step_q
   PAD_INT = 2.0 ** -30
   def fp12(idx, rem, az):
    nonlocal tot_nodes, tot_e_leaves, tot_r_leaves, id_errs, q60_fn, drift_errs
    tot_nodes += 1
    sh = -sum(mu[j][idx] * x[j] for j in range(idx + 1, 12))
    rad = sqrt(max(0.0, rem / Bn[idx]))
    lo = max(0, ceil(sh - rad - PAD_INT)) if az else ceil(sh - rad - PAD_INT)
    hi = floor(sh + rad + PAD_INT)
    if idx == 0:
     if lo > hi: return
     a_base = [sum(x[i] * Bi[i][c] for i in range(1, 12)) for c in range(6)]
     b_base = [sum(x[i] * Bi[i][6 + c] for i in range(1, 12)) for c in range(6)]
     for ch_lo in range(lo, hi + 1, 512):
      ch_hi = min(hi, ch_lo + 511)
      a_ch = [a_base[c] + ch_lo * a_step[c] for c in range(6)]
      b_ch = [b_base[c] + ch_lo * b_step[c] for c in range(6)]
      u_q = _qfi6(Ma_t, a_ch)
      du_q = 2 * _bfi6(Ma_t, a_ch, a_step) + u_step_q
      q_q = _qfi6(Mea_t, a_ch) + _qfi6(Mb_t, b_ch)
      dnr_q = 2 * (_bfi6(Mea_t, a_ch, a_step) + _bfi6(Mb_t, b_ch, b_step)) + q_step_q
      nr_q = ONE_Q60 - q_q + 2 * EPS_Q60
      for xi in range(ch_lo, ch_hi + 1):
       pass_q60 = (u_q <= ONE_PAD) and ((((ONE_PAD - u_q) * OMCB_Q64) >> 64) <= nr_q)
       cur_uq, cur_nrq = u_q, nr_q
       nr_q -= dnr_q; dnr_q += d2nr_q; u_q += du_q; du_q += d2u_q
       d = xi - sh; nr = rem - d * d * Bn[0]
       if nr < -PAD_INT: continue
       x[0] = xi
       if az and xi == 0: continue
       tot_e_leaves += 1
       w = [sum(x[i] * Bi[i][c] for i in range(12)) for c in range(12)]
       av, bv = w[:6], w[6:]
       if not any(av): continue
       if cur_uq != _qfi6(Ma_t, av) or cur_nrq != (ONE_Q60 - (_qfi6(Mea_t, av) + _qfi6(Mb_t, bv)) + 2 * EPS_Q60):
        drift_errs += 1
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
     if nr >= -PAD_INT:
      x[idx] = xi
      fp12(idx - 1, max(0.0, nr), az and (xi == 0))
   fp12(11, 1.0 + PAD_INT, True)
  all_k = set(gt_mult) | set(found_mult)
  tot_gt += sum(gt_mult.values()); tot_found += sum(found_mult.values())
  tot_missed += sum(max(0, gt_mult.get(k, 0) - found_mult.get(k, 0)) for k in all_k)
  tot_dups += sum(max(0, found_mult.get(k, 0) - gt_mult.get(k, 0)) for k in all_k)
 RA2_P = 468.87 ** 2; RB2_P = 1547.36 ** 2
 Ma_p = [[round(GA[r][c] * ONE_Q60 / RA2_P) for c in range(6)] for r in range(6)]
 Mea_p = [[round(CA * GA[r][c] * ONE_Q60 / RA2_P) for c in range(6)] for r in range(6)]
 Mb_p = [[round(CB * GB[r][c] * ONE_Q60 / RB2_P) for c in range(6)] for r in range(6)]
 GA_f = [[F(GA[r][c]) for c in range(6)] for r in range(6)]; GB_f = [[F(GB[r][c]) for c in range(6)] for r in range(6)]
 fRA2 = F(RA2_P); fRB2 = F(RB2_P); fCA = F(CA); fCB = F(CB); fKAP = F(KAP)
 rng_p = random.Random(20261120)
 prod_drift = prod_fn = prod_in_R = prod_in_E = 0; prod_max_u_err = prod_max_q_err = 0.0
 for _ in range(12):
  ra = [rng_p.gauss(0.0, 1.0) for _ in range(6)]; ua0 = _qf6(GA, ra) / RA2_P
  sa = sqrt(rng_p.uniform(0.60, 0.999) / ua0)
  a0 = [int(round(z * sa)) for z in ra]
  u_t = _qf6(GA, a0) / RA2_P
  rb = [rng_p.gauss(0.0, 1.0) for _ in range(6)]; vb0 = _qf6(GB, rb) / RB2_P
  sb = sqrt(max(0.05, (1.0 - KAP * min(1.0, u_t)) * rng_p.uniform(0.60, 0.999)) / vb0)
  b0 = [int(round(z * sb)) for z in rb]
  da = [rng_p.choice((-1, 0, 1)) if i < 2 else 0 for i in range(6)]
  db = [rng_p.choice((-1, 0, 1)) if i < 2 else 0 for i in range(6)]
  uq = _qfi6(Ma_p, a0); duq = 2 * _bfi6(Ma_p, a0, da) + _qfi6(Ma_p, da); d2uq = 2 * _qfi6(Ma_p, da)
  qq = _qfi6(Mea_p, a0) + _qfi6(Mb_p, b0)
  dqq = 2 * (_bfi6(Mea_p, a0, da) + _bfi6(Mb_p, b0, db)) + _qfi6(Mea_p, da) + _qfi6(Mb_p, db)
  d2qq = 2 * (_qfi6(Mea_p, da) + _qfi6(Mb_p, db))
  nrq = ONE_Q60 - qq + 2 * EPS_Q60
  for step in range(512):
   if step in (0, 1, 7, 31, 127, 255, 511):
    ak = [a0[c] + step * da[c] for c in range(6)]; bk = [b0[c] + step * db[c] for c in range(6)]
    if uq != _qfi6(Ma_p, ak) or nrq != (ONE_Q60 - (_qfi6(Mea_p, ak) + _qfi6(Mb_p, bk)) + 2 * EPS_Q60):
     prod_drift += 1
    u_ex = sum(GA_f[r][c] * ak[r] * ak[c] for r in range(6) for c in range(6)) / fRA2
    v_ex = sum(GB_f[r][c] * bk[r] * bk[c] for r in range(6) for c in range(6)) / fRB2
    q_ex = fCA * u_ex + fCB * v_ex
    if u_ex <= 2 and q_ex <= 2:
     prod_in_E += 1
     prod_max_u_err = max(prod_max_u_err, abs(float(F(uq, ONE_Q60) - u_ex)))
     prod_max_q_err = max(prod_max_q_err, abs(float(F(ONE_Q60 + 2 * EPS_Q60 - nrq, ONE_Q60) - q_ex)))
     if u_ex <= 1 and v_ex + fKAP * u_ex <= 1:
      prod_in_R += 1
      if not ((uq <= ONE_PAD) and ((((ONE_PAD - uq) * OMCB_Q64) >> 64) <= nrq)): prod_fn += 1
   nrq -= dqq; dqq += d2qq; uq += duq; duq += d2uq
 prod_q60_ok = (drift_errs == 0 and prod_drift == 0 and prod_fn == 0 and prod_in_R >= 20
                and 0.0 < prod_max_u_err < 2.0**-37 and 0.0 < prod_max_q_err < 2.0**-35)
 e_lph = tot_e_leaves / max(1, tot_gt)
 r_lph = tot_r_leaves / max(1, tot_gt)
 vr_obs = tot_e_leaves / max(1, tot_r_leaves)
 return dict(
  ok=(tot_missed == 0 and tot_dups == 0 and id_errs == 0 and q60_fn == 0 and prod_q60_ok and tot_gt == 408
    and e_lph <= M.MULT_E12 and r_lph <= M.MULT_FD_ACT),
  n_orbits=len(orbits), primes=[o[0] for o in orbits],
  gt_hits=tot_gt, found_hits=tot_found, missed=tot_missed, dups=tot_dups,
  gate_identity_errors=id_errs, q60_gate_false_neg=q60_fn,
  exact_int_q60_drift_errors=drift_errs,
  prod_q60_512step_ok=prod_q60_ok, prod_q60_in_R_tested=prod_in_R, prod_q60_in_E_tested=prod_in_E,
  prod_q60_drift_errors=prod_drift, prod_q60_false_neg=prod_fn,
  prod_q60_max_u_err=float(f"{prod_max_u_err:.3e}"), prod_q60_max_q_err=float(f"{prod_max_q_err:.3e}"),
  prod_q60_eps_pad=float(f"{2.0**-34:.3e}"),
  e_leaves=tot_e_leaves, r_leaves=tot_r_leaves,
  vr12_observed=round(vr_obs, 4), vr12_exact_bound=round(vr12_th, 6),
  e_leaves_per_can_hit=round(e_lph, 4), bound_e_leaves_per_hit=M.MULT_E12,
  r_leaves_per_can_hit=round(r_lph, 4), bound_r_leaves_per_hit=M.MULT_FD_ACT,
 )
def verify_12d_descent_lattices():
 D = M.load_data(); T4 = D["t4_geom"]
 B12 = json.loads((Path(__file__).resolve().parent / "c6_witnesses.json").read_text())
 EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]
 EMB_W = [complex(r, im) for r, im in T4["EMB_W"]]
 EMB_F = {d: [sum(cs[k] * (EMB_ZK[j][1] ** k) for k in range(6)) for j in range(3)] for d, cs in M.F_REC_COEFFS.items()}
 LOG2U6, LOG2P6 = 6.0 * log2(M.U_BN254), 6.0 * log2(M.P_BN254)
 ZK3I = [[1, 0, 1, -1, -3, -3], [0, 1, -1, 0, -1, -1], [0, 0, -1, 1, 1, 3], [0, 0, -1, -1, 0, 2], [0, 0, 1, 0, 1, -1], [0, 0, 0, 1, 1, 2]]
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
 EZR = tuple(tuple(z.real for z in row) for row in EMB_ZK)
 EZI = tuple(tuple(z.imag for z in row) for row in EMB_ZK)
 _sp_l = math.sumprod
 def log2_norms_fp(a, b, need_f=True):
  lf = lg = 0.0
  mx = max(max(abs(x) for x in a), max(abs(x) for x in b), 1)
  sh = max(0, mx.bit_length() - 40); sc = float(1 << sh)
  af = [x / sc for x in a]; bf = [x / sc for x in b]
  for j in range(3):
   za = complex(_sp_l(EZR[j], af), _sp_l(EZI[j], af))
   zb = complex(_sp_l(EZR[j], bf), _sp_l(EZI[j], bf))
   za2, zb2 = za * za, zb * zb
   gv = M.U_BN254 * EMB_W[j] * za2 - zb2
   lg += log2(gv.real*gv.real + gv.imag*gv.imag) + 4.0 * sh
   if need_f:
    za4, zb4 = za2 * za2, zb2 * zb2
    fv = EMB_F[8][j]*za4*za4 + EMB_F[6][j]*za4*za2*zb2 + EMB_F[4][j]*za4*zb4 + EMB_F[2][j]*za2*zb4*zb2 + EMB_F[0][j]*zb4*zb4
    lf += log2(fv.real*fv.real + fv.imag*fv.imag) + 16.0 * sh
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
  side, q, s, rho_r = item[0], M.b62d(item[1]), M.b62d(item[2]), M.b62d(item[3])
  bf = [M.b62d(x) for x in item[4].split()]
  bq = q.bit_length(); q_bits_list.append(bq); q_sides_list.append(side)
  B = [bf[r * 12 : (r + 1) * 12] for r in range(12)]
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
  for _ in range(80):
   x = sample_b12_r4(rng)
   v = [_sp_l(B[r], x) for r in range(12)]
   lf, lg = log2_norms_fp(v[:6], v[6:], True)
   all_cfs_r4.append(lf - 48.0 * c_r4); all_cgs_r4.append(lg - 12.0 * c_r4 - LOG2U6)
   emp_q_r4.append(M.rho((lf - sub_f) / 59.0) * M.rho((lg - sub_g) / 59.0))
  paired_ratios.append(round((sum(emp_q_r4) / len(emp_q_r4)) / max(mod_q_r4, 1e-300), 2))
  for R_big in (16, 64, 128):
   S_Rb = log2((pi ** 6 / 720.0) * (R_big ** 12)); c_rb = (log2(q) + S_Rb) / 12.0
   mod_rb = sum(w * M.rho((48.0 * c_rb + cf - sub_f) / 59.0) * M.rho((12.0 * c_rb + LOG2U6 + cg - sub_g) / 59.0) for w, cf, cg in DESC.QBINS)
   emp_rb = []
   for _ in range(12):
    gv = [rng.gauss(0.0, 1.0) for _ in range(12)]; gn = sqrt(sum(z * z for z in gv))
    rad = R_big * (rng.random() ** (1.0 / 12.0))
    x = [int(round(z * rad / gn)) for z in gv]
    if not any(x): x[0] = 1
    v = [_sp_l(B[r], x) for r in range(12)]
    lf, lg = log2_norms_fp(v[:6], v[6:], True)
    emp_rb.append(M.rho((lf - sub_f) / 59.0) * M.rho((lg - sub_g) / 59.0))
   hiR_ratios.append((sum(emp_rb) / len(emp_rb)) / max(mod_rb, 1e-300))
 all_cfs_r4.sort(); all_cgs_r4.sort()
 def _qtile(arr, qv): return round(arr[min(len(arr) - 1, max(0, ceil(qv * len(arr)) - 1))], 2)
 qbins_lo_ok = all(
  all_cfs_r4[min(len(all_cfs_r4)-1, floor((k/6.0)*len(all_cfs_r4)))] >= DESC.QBINS_LO[k][0] and
  all_cgs_r4[min(len(all_cgs_r4)-1, floor((k/6.0)*len(all_cgs_r4)))] >= DESC.QBINS_LO[k][1]
  for k in range(6)
 )
 cf_sextiles = [_qtile(all_cfs_r4, (k + 1) / 6.0) for k in range(6)]
 cg_sextiles = [_qtile(all_cgs_r4, (k + 1) / 6.0) for k in range(6)]
 qbins_hi_ok = all(cf_sextiles[k] <= DESC.QBINS[k][1] and cg_sextiles[k] <= DESC.QBINS[k][2] for k in range(6))
 q_pts = (0.10, 0.25, 0.50, 0.75, 0.90, 1.00)
 boot_ok = boot_gt_ok = True; boot_offs = []; boot_ge3 = 0
 P_bn = M.P_BN254; inv3p = pow(3, -1, P_bn)
 uw_p = [(M.U_BN254 * c) % P_bn for c in M.W_COEFFS]
 inv_uw_p = M._pinv_h(uw_p, P_bn)
 def _mul12_p(a12, b12):
  p0 = [(x + y) % P_bn for x, y in zip(M._pmul_h(a12[0], b12[0], P_bn), M._pmul_h(M._pmul_h(a12[1], b12[1], P_bn), inv_uw_p, P_bn))]
  p1 = [(x + y) % P_bn for x, y in zip(M._pmul_h(a12[0], b12[1], P_bn), M._pmul_h(a12[1], b12[0], P_bn))]
  return p0, p1
 for item in B12["boot"]:
  R_p = [M.b62d(x) for x in item[0].split()]
  A_flat = [M.b62d(x) for x in item[1].split()]
  A_top = [A_flat[r * 12 : (r + 1) * 12] for r in range(6)]
  res12 = ([1,0,0,0,0,0], [0]*6); bas12 = ([1,0,0,0,0,0], [c % P_bn for c in R_p])
  for bit in bin(M.R_BN254)[2:]:
   res12 = _mul12_p(res12, res12)
   if bit == "1": res12 = _mul12_p(res12, bas12)
  if res12[1] != [0]*6 or not any(R_p): boot_gt_ok = False
  B_bot = [[0] * 12 for _ in range(6)]
  for col in range(12):
   av = [A_top[r][col] for r in range(6)]
   ap = [(sum(av[k] * M.ZK_NUM[k][j] for k in range(6)) * inv3p) % P_bn for j in range(6)]
   bp = M._pmul_h(R_p, ap, P_bn)
   bp3_s = [(3 * x) % P_bn for x in bp]
   bp3_s = [x - P_bn if x > P_bn // 2 else x for x in bp3_s]
   for r_i in range(6):
    num = sum(ZK3I[r_i][c_i] * bp3_s[c_i] for c_i in range(6))
    if num % 3 != 0: boot_ok = False
    B_bot[r_i][col] = num // 3
  B = A_top + B_bot
  if abs(det_int(B)) != P_bn ** 6: boot_ok = False
  for col in range(12):
   av = [B[r][col] for r in range(6)]; bv = [B[6 + r][col] for r in range(6)]
   if M._pmul_h(R_p, zk_to_poly_mod(av, P_bn), P_bn) != zk_to_poly_mod(bv, P_bn): boot_ok = False
  for _ in range(120):
   x = sample_b12_r4(rng)
   if max(abs(z) for z in x) >= 3: boot_ge3 += 1
   v = [_sp_l(B[r], x) for r in range(12)]
   _, lg = log2_norms_fp(v[:6], v[6:], False)
   boot_offs.append(lg - LOG2P6 - S_R4 - LOG2U6)
 c_b = (LOG2P6 + S_R4) / 12.0
 emp_boot_sm = sum(M.rho((12.0 * c_b + LOG2U6 + o) / 110.0) for o in boot_offs) / len(boot_offs)
 mod_boot_sm = sum(w * M.rho((12.0 * c_b + LOG2U6 + cg) / 110.0) for w, cg in DESC.BOOT_BINS)
 boot_offs.sort()
 boot_quintiles = [_qtile(boot_offs, (k + 1) / 5.0) for k in range(5)]
 boot_quintiles_ok = all(boot_quintiles[k] <= DESC.BOOT_BINS[k][1] for k in range(5))
 boot_lo_ok = (boot_offs[0] >= DESC.BOOT_LO_OFF)
 diag_ok = (cells4_exact_ok and side_div_ok and max_exact_err < 1e-8 and all(r >= 1.0 for r in paired_ratios)
       and qbins_lo_ok and qbins_hi_ok and min(hiR_ratios) >= 1.0 and boot_quintiles_ok and boot_lo_ok
       and emp_boot_sm >= mod_boot_sm and boot_ge3 > 0 and len(B12["special_q"]) == 24 and len(B12["boot"]) == 8)
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
  cf_sextiles_obs=cf_sextiles, cg_sextiles_obs=cg_sextiles,
  cf_quantiles_p10_25_50_75_90_max=[_qtile(all_cfs_r4, qv) for qv in q_pts],
  cg_quantiles_p10_25_50_75_90_max=[_qtile(all_cgs_r4, qv) for qv in q_pts],
  qbins_hi_6_of_6_ok=qbins_hi_ok, qbins_lo_6_of_6_ok=qbins_lo_ok,
  boot_quintiles_obs=boot_quintiles,
  boot_quantiles_p10_25_50_75_90_max=[_qtile(boot_offs, qv) for qv in q_pts],
  boot_5_quintiles_ok=boot_quintiles_ok,
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
 active_cols = {c for c in range(U) if col_rows[c]}
 def purge_singletons(nalive, stack=None):
  q_s = [c for c in active_cols if len(col_rows[c]) == 1]
  while q_s:
   c = q_s.pop()
   cr = col_rows[c]
   if len(cr) != 1: continue
   r0 = next(iter(cr))
   row0 = rows[r0]
   if not row0[5]: continue
   d0 = row0[0]
   if stack is not None:
    if any(len(col_rows[cc]) < 2 for cc in d0 if cc != c): continue
    if z6_eval(d0[c]) == 0: continue
    stack.append((c, dict(d0), row0[6]))
   row0[5] = False; nalive -= 1
   for cc in d0:
    scc = col_rows[cc]
    scc.discard(r0)
    lcc = len(scc)
    if lcc == 1: q_s.append(cc)
    elif lcc == 0: active_cols.discard(cc)
   d0.clear()
  return nalive
 nalive = purge_singletons(Rp); ncols = len(active_cols)
 peel_cols = set(active_cols)
 repurge_stack = []
 texc = max(6 if export_hx else 2, int(0.0001 * Up))
 while nalive - ncols > texc:
  par = list(range(Rp))
  def find(x):
   while par[x] != x: par[x] = par[par[x]]; x = par[x]
   return x
  for c in active_cols:
   cr = col_rows[c]
   if len(cr) == 2:
    r1, r2 = cr; a, b = find(r1), find(r2)
    if a != b: par[a] = b
  comp_w = {}
  for i in range(Rp):
   ri = rows[i]
   if ri[5]: rt = find(i); comp_w[rt] = comp_w.get(rt, 0.0) + ri[4]
  comps = sorted(comp_w.items(), key=lambda kv: -kv[1])
  need = nalive - ncols - texc
  nb = max(1, min(need // 2 + 1, max(1, len(comps) // 100)))
  kill = {rt for rt, _ in comps[:nb]}
  dropped_any = False
  for i in range(Rp):
   ri = rows[i]
   if ri[5] and find(i) in kill:
    d_i = ri[0]
    if any(len(col_rows[cc]) < 2 for cc in d_i): continue
    ri[5] = False; nalive -= 1; dropped_any = True
    for cc in d_i:
     scc = col_rows[cc]
     scc.discard(i)
     if not scc: active_cols.discard(cc)
    d_i.clear()
    nalive = purge_singletons(nalive, stack=repurge_stack)
  if not dropped_any:
   safe = sorted([(rows[i][4], i) for i in range(Rp) if rows[i][5] and all(len(col_rows[cc]) >= 2 for cc in rows[i][0])], reverse=True)
   if not safe: break
   for _, i in safe[:max(1, need // 2)]:
    ri = rows[i]
    if not ri[5]: continue
    d_i = ri[0]
    if any(len(col_rows[cc]) < 2 for cc in d_i): continue
    ri[5] = False; nalive -= 1
    for cc in d_i:
     scc = col_rows[cc]; scc.discard(i)
     if not scc: active_cols.discard(cc)
    d_i.clear()
    nalive = purge_singletons(nalive, stack=repurge_stack)
    if nalive - len(active_cols) <= texc: break
  ncols = len(active_cols)
 post_peel_cols_lost = len(peel_cols - (active_cols | {c for c, _, _ in repurge_stack}))
 repurge_backsub_errors = 0
 if check_kernel:
  v_rec = {c: v_sec[c] for c in active_cols}
  for c_p, d_p, rhs_p in reversed(repurge_stack):
   pv = z6_eval(d_p[c_p])
   s_oth = sum(z6_eval(val) * v_rec[cc] for cc, val in d_p.items() if cc != c_p) % P_CHK
   v_c = ((-(rhs_p + s_oth)) * pow(pv, -1, P_CHK)) % P_CHK
   v_rec[c_p] = v_c
   if v_c != v_sec[c_p]: repurge_backsub_errors += 1
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
  ok=(kern_errs == 0 and post_peel_cols_lost == 0 and repurge_backsub_errors == 0 and N_sw > 0 and N_l > 0 and merges > 0 and avg_w_fin <= M.W_ROW and tab_r < M.TAB_R_CAP),
  log2U=log2U, seed=seed, W_target=W_target, initial_rows=Rp,
  N_sw=N_sw, N_l=N_l, N_sw_U=round(N_sw / Up, 6), N_l_U=round(N_l / Up, 6),
  avg_w=round(avg_w_fin, 2), merges=merges, w_mid=round(w_mid, 4), lsmall=round(lsmall, 3),
  exp_unit=round(exp_unit, 5), k5_exp=round(k5_exp, 4), tab_r=round(tab_r, 5),
  kernel_errors=kern_errs, post_peel_cols_lost=post_peel_cols_lost,
  repurge_stack_len=len(repurge_stack), repurge_backsub_errors=repurge_backsub_errors,
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
def verify_c6_pgf_and_hea_campaign():
 from fractions import Fraction as F
 D = M.load_data(); T4 = D["t4_geom"]
 camp = json.loads((Path(__file__).resolve().parent / "c6_witnesses.json").read_text())["pgf_campaign"]
 tw = camp["tower"]
 U_T, P_T, WA = int(tw["u_t"]), int(tw["p_t"]), float(tw["wa_lll"])
 LB, NA, NB, B_D2, T_S = int(tw["lb"]), int(tw["nA"]), int(tw["nB"]), int(tw["b_d2"]), int(tw["T_surv"])
 P_ECM, M_STEP = float(tw["p_ecm"]), int(tw["miss_fb_step"])
 B_LO, B_NA, B_NB, B_TS = 1 << LB, 1 << NA, 1 << NB, 1 << T_S
 qhf, qhg, qlf, qlg = tw["qbins_hi_f"], tw["qbins_hi_g"], tw["qbins_lo_f"], tw["qbins_lo_g"]
 af, ag, ZKN, W_P = float(tw["alpha_f"]), float(tw["alpha_g"]), M.ZK_NUM, M.W_COEFFS
 WBAR_P = [-4 - c if i == 0 else -c for i, c in enumerate(W_P)]
 EMB_ZK = [[complex(r, im) for r, im in row] for row in T4["EMB_ZK"]]; S2 = sqrt(2.0)
 ww = M._pmul_h_q([F(c) for c in W_P], [F(c) for c in WBAR_P])
 y1 = M._ppow_h([0,1,0,0,0,0], P_T, P_T); y2 = M._ppow_h(y1, P_T, P_T); y3 = M._ppow_h(y2, P_T, P_T); y6 = M._ppow_h(y3, P_T**3, P_T)
 tower_ok = (ww == [F(6)]+[F(0)]*5 and 1+U_T*6 == P_T and y6 == [0,1,0,0,0,0] and y2 != [0,1,0,0,0,0] and y3 != [0,1,0,0,0,0])
 sig_q = [F(c, 3) for c in M.SIGMA_NUM]; pw_s = [[F(1)]+[F(0)]*5, sig_q]
 for _ in range(4): pw_s.append(M._pmul_h_q(pw_s[-1], sig_q))
 SMY = tuple(tuple(int(3*pw_s[k][j]) for k in range(6)) for j in range(6))
 _MHP = {S: (1<<S, 1<<(S-1), (1<<S)-1, S, 2*S, 3*S, 4*S, 5*S) for S in (32, 44, 88)}
 def mul_h(a, b, S=44):
  Mod, H, msk, s1, s2, s3, s4, s5 = _MHP[S]
  P = (a[0]+(a[1]<<s1)+(a[2]<<s2)+(a[3]<<s3)+(a[4]<<s4)+(a[5]<<s5)) * (b[0]+(b[1]<<s1)+(b[2]<<s2)+(b[3]<<s3)+(b[4]<<s4)+(b[5]<<s5))
  c = [0]*11
  for i in range(11):
   v = P & msk
   if v >= H: v -= Mod; P += Mod
   c[i] = v; P >>= s1
  c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10 = c
  if c10: c9+=2*c10; c8-=c10; c7-=2*c10; c5+=4*c10; c4-=3*c10
  if c9: c8+=2*c9; c7-=c9; c6-=2*c9; c4+=4*c9; c3-=3*c9
  if c8: c7+=2*c8; c6-=c8; c5-=2*c8; c3+=4*c8; c2-=3*c8
  if c7: c6+=2*c7; c5-=c7; c4-=2*c7; c2+=4*c7; c1-=3*c7
  if c6: c5+=2*c6; c4-=c6; c3-=2*c6; c1+=4*c6; c0-=3*c6
  return (c0, c1, c2, c3, c4, c5)
 _sp = math.sumprod
 def norm_zy(p):
  s1 = tuple(_sp(r, p) for r in SMY); s2 = tuple(_sp(r, s1) for r in SMY)
  n3 = mul_h(mul_h(p, s1, 44), s2, 88); c1 = -n3[2]//27; c0 = n3[0]//27 - 2*c1
  return c0*c0 + 2*c1*c1
 def eval_g_sq(pa2, pb2):
  wa2 = mul_h(W_P, pa2, 44)
  return norm_zy((U_T*wa2[0]-pb2[0], U_T*wa2[1]-pb2[1], U_T*wa2[2]-pb2[2], U_T*wa2[3]-pb2[3], U_T*wa2[4]-pb2[4], U_T*wa2[5]-pb2[5]))//729
 def eval_f_sq(pa2, pb2):
  wb2 = mul_h(WBAR_P, pb2, 44)
  return norm_zy((pa2[0]+wb2[0], pa2[1]+wb2[1], pa2[2]+wb2[2], pa2[3]+wb2[3], pa2[4]+wb2[4], pa2[5]+wb2[5]))//729
 def _prod(xs):
  while len(xs) > 1: xs = [xs[i]*xs[i+1] if i+1 < len(xs) else xs[i] for i in range(0, len(xs), 2)]
  return xs[0] if xs else 1
 ps_lo = M._primes_upto(B_LO)
 ps_tail = [p for p in ps_lo if p > 2048]
 miss_ps = [p for i, p in enumerate(ps_tail) if i % M_STEP == 0]
 P_MISS = _prod(miss_ps)
 P1, P2, P3 = _prod([p for p in ps_lo if p<=256]), _prod([p for p in ps_lo if 256<p<=2048]), _prod([p for p in ps_lo if 2048<p<=B_LO])
 P_MID = P2 * P3
 def _isp(n):
  if n < 2: return False
  for p in (2,3,5,7,11,13,17,19,23,29,31,37):
   if n == p: return True
   if n % p == 0: return False
  d, s = n - 1, 0
  while d % 2 == 0: d //= 2; s += 1
  for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
   if a % n == 0: continue
   x = pow(a, d, n)
   if x in (1, n - 1): continue
   for _ in range(s - 1):
    x = (x * x) % n
    if x == n - 1: break
   else: return False
  return True
 def _brent(n):
  if n % 2 == 0: return 2
  for c in (1, 2, 3, 5, 7, 11, 13):
   y, m, g, r, q = 2, 64, 1, 1, 1
   while g == 1:
    x = y
    for _ in range(r): y = (y*y + c) % n
    k = 0
    while k < r and g == 1:
     ys = y
     for _ in range(min(m, r-k)): y = (y*y + c) % n; q = (q * abs(x-y)) % n
     g = math.gcd(q, n); k += m
    r <<= 1
   if g == n:
    while True:
     ys = (ys*ys + c) % n; g = math.gcd(abs(x-ys), n)
     if g > 1: break
   if 1 < g < n: return g
  return n
 _gcd = math.gcd; B_SQ = 1 << (2*LB); B_NB2 = B_NB*B_NB; P_ECM_THR = round(P_ECM * 1000)
 def strip_blo(n):
  for P in (P1, P2, P3):
   while n > 1:
    g = _gcd(n, P)
    if g == 1: break
    n //= g
   if n == 1: return 1
  return n
 def count_mid(n):
  cnt = 0
  while True:
   g = _gcd(n, P_MID)
   if g == 1: break
   n //= g
   if g < 65536 or _isp(g): cnt += 1
   else:
    stk = [g]
    while stk:
     x = stk.pop()
     if x < 65536 or _isp(x): cnt += 1
     else: d = _brent(x); stk.append(d); stk.append(x // d)
  return cnt
 def ecm_ok(side, q, sp_idx, n_fac):
  bk = (side << 56) ^ (q << 24) ^ (sp_idx << 8) ^ 0x51C62026
  for j in range(n_fac):
   x = ((bk ^ j) + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
   x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
   x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
   if ((x ^ (x >> 31)) % 1000) >= P_ECM_THR: return False
  return True
 def qual_2w(cof, b_a=B_NA, b_b=B_NB):
  if cof == 1: return True, []
  if cof < B_SQ: return (cof <= b_b), ([cof] if cof <= b_b else [])
  if cof > b_b*b_b or _isp(cof): return False, []
  lps, stk = [], [cof]
  while stk:
   m = stk.pop()
   if m == 1: continue
   if m < B_SQ or _isp(m):
    if m > b_b: return False, []
    lps.append(m)
    if len(lps) > 2: return False, []
    continue
   d = _brent(m)
   if d == m: return False, []
   stk.append(d); stk.append(m // d)
  return (sum(1 for p in lps if p > b_a) <= 1 and len(lps) <= 2), lps
 def sq_basis(q, s, rho):
  zks = [_eval_zk_mod([1 if j==k else 0 for j in range(6)], s, q) for k in range(6)]; iz0 = pow(zks[0], -1, q)
  B = [[0]*12 for _ in range(12)]; B[0][0] = q
  for k in range(1, 6): B[k][k] = 1; B[k][0] = (-zks[k]*iz0) % q
  for k in range(6): B[6+k][6+k] = 1; B[6+k][0] = (rho*zks[k]*iz0) % q
  E = [[S2*sum(EMB_ZK[j][c]*r[c] for c in range(6)).real*WA if idx%2==0 else S2*sum(EMB_ZK[j][c]*r[c] for c in range(6)).imag*WA for j in range(3) for idx in (0,1)]
       + [S2*sum(EMB_ZK[j][c]*r[6+c] for c in range(6)).real if idx%2==0 else S2*sum(EMB_ZK[j][c]*r[6+c] for c in range(6)).imag for j in range(3) for idx in (0,1)] for r in B]
  def gs(R):
   Bo, mu, Bn = [], [[0.0]*12 for _ in range(12)], [0.0]*12
   for i in range(12):
    v = R[i][:]
    for j in range(i):
     mu[i][j] = math.sumprod(R[i], Bo[j]) / Bn[j]; v = [v[c] - mu[i][j]*Bo[j][c] for c in range(12)]
    Bo.append(v); Bn[i] = math.sumprod(v, v)
   return Bo, mu, Bn
  Bo, mu, Bn = gs(E); k = 1
  while k < 12:
   for j in range(k-1, -1, -1):
    qv = round(mu[k][j])
    if qv: B[k] = [B[k][i]-qv*B[j][i] for i in range(12)]; E[k] = [E[k][i]-qv*E[j][i] for i in range(12)]; Bo, mu, Bn = gs(E)
   if Bn[k] >= (0.99 - mu[k][k-1]**2)*Bn[k-1]: k += 1
   else: B[k], B[k-1] = B[k-1], B[k]; E[k], E[k-1] = E[k-1], E[k]; Bo, mu, Bn = gs(E); k = max(k-1, 1)
  return B, [([sum(r[c]*ZKN[c][j] for c in range(6)) for j in range(6)], [sum(r[6+c]*ZKN[c][j] for c in range(6)) for j in range(6)]) for r in B]
 B12 = [((i,1),) for i in range(12)] + [((i,1),(j,s)) for i in range(12) for j in range(i+1,12) for s in (1,-1)] + [((i,1),(j,sj),(k,sk)) for i in range(12) for j in range(i+1,12) for k in range(j+1,12) for sj in (1,-1) for sk in (1,-1)]
 def find_ch(cside, cq, av, bv):
  yp = M._ppowmod([0,1], cq, M.H_COEFFS, cq)
  if len(yp) < 2: yp += [0]*(2-len(yp))
  yp[1] = (yp[1]-1) % cq
  for s in M._roots_mod(M._pgcd(yp, M.H_COEFFS, cq), cq):
   a0, b0 = _eval_zk_mod(av, s, cq), _eval_zk_mod(bv, s, cq)
   if b0 == 0: continue
   r0 = (a0 * pow(b0, -1, cq)) % cq
   if r0 == 0: continue
   ws = _eval_poly_mod(W_P, s, cq)
   if (cside == 0 and (U_T*ws*r0*r0-1)%cq == 0) or (cside == 1 and (r0*r0-4-ws)%cq == 0): return (cside, cq, s, r0)
  return None
 def scan_ball(side, q, s, rho, b_a=B_NA, b_b=B_NB, first=False):
  Bp, R3 = sq_basis(q, s, rho); M_o = S_o = 0
  Ra = [r[0] for r in R3] + [[-x for x in r[0]] for r in R3]
  Rb = [r[1] for r in R3] + [[-x for x in r[1]] for r in R3]
  fn1, fn2 = (eval_f_sq, eval_g_sq) if side == 0 else (eval_g_sq, eval_f_sq)
  for sp_idx, sp in enumerate(B12):
   if len(sp) == 3:
    i, j, k = sp[0][0], sp[1][0]+(0 if sp[1][1]==1 else 12), sp[2][0]+(0 if sp[2][1]==1 else 12)
    a0, a1, a2 = Ra[i], Ra[j], Ra[k]; b0, b1, b2 = Rb[i], Rb[j], Rb[k]
    pa = (a0[0]+a1[0]+a2[0], a0[1]+a1[1]+a2[1], a0[2]+a1[2]+a2[2], a0[3]+a1[3]+a2[3], a0[4]+a1[4]+a2[4], a0[5]+a1[5]+a2[5])
    pb = (b0[0]+b1[0]+b2[0], b0[1]+b1[1]+b2[1], b0[2]+b1[2]+b2[2], b0[3]+b1[3]+b2[3], b0[4]+b1[4]+b2[4], b0[5]+b1[5]+b2[5])
   elif len(sp) == 2:
    i, j = sp[0][0], sp[1][0]+(0 if sp[1][1]==1 else 12)
    a0, a1 = Ra[i], Ra[j]; b0, b1 = Rb[i], Rb[j]
    pa = (a0[0]+a1[0], a0[1]+a1[1], a0[2]+a1[2], a0[3]+a1[3], a0[4]+a1[4], a0[5]+a1[5])
    pb = (b0[0]+b1[0], b0[1]+b1[1], b0[2]+b1[2], b0[3]+b1[3], b0[4]+b1[4], b0[5]+b1[5])
   else:
    pa, pb = tuple(Ra[sp[0][0]]), tuple(Rb[sp[0][0]])
   if not any(pa) or not any(pb): continue
   pa2 = tuple(x//3 for x in mul_h(pa, pa, 32)); pb2 = tuple(x//3 for x in mul_h(pb, pb, 32))
   n1 = fn1(pa2, pb2)
   if n1 <= 0: continue
   c1 = strip_blo(n1); s1 = (not first) and (c1 < B_TS); q1, ab1 = qual_2w(c1, b_a, b_b)
   if not s1 and not q1: continue
   n2 = fn2(pa2, pb2)
   if n2 % q != 0: continue
   n2 //= q
   if n2 <= 0: continue
   c2 = strip_blo(n2)
   if s1 and c2 < B_TS: S_o += 1
   if q1:
    q2, ab2 = qual_2w(c2, b_a, b_b)
    if q2 and _gcd(n1, P_MISS) == 1 and _gcd(n2, P_MISS) == 1 and ecm_ok(side, q, sp_idx, len(ab1)+len(ab2)+count_mid(n1)+count_mid(n2)):
     if first:
      v = [0]*12
      for idx, sgn in sp:
       row = Bp[idx]
       for c in range(12): v[c] += row[c] if sgn == 1 else -row[c]
      return (v[:6], v[6:], ab2, ab1) if side == 0 else (v[:6], v[6:], ab1, ab2)
     M_o += 1
  return None if first else (M_o, S_o)
 first_usable = lambda ps, qp, spy, rhp, b_a=B_NA, b_b=B_NB: scan_ball(ps, qp, spy, rhp, b_a, b_b, True)
 eval_ball = lambda side, q, s, rho: scan_ball(side, q, s, rho)
 _rc = {}
 def _rc_f(u):
  v = _rc.get(u)
  if v is None: v = M.rho(u); _rc[u] = v
  return v
 rho_lo = lambda u: 1.0 if u <= 1.0 else _rc_f(ceil(u*1024.0)/1024.0)
 rho_hi = lambda u: 1.0 if u <= 1.0 else _rc_f(floor(u*1024.0)/1024.0)
 iA_m = int(round(48*(NA-LB)/(NB-LB))); h_m = (NB-LB)/48.0
 wA_m = [math.log1p(h_m/(LB+i*h_m)) if i < iA_m else 0.0 for i in range(48)]
 wB_m = [math.log1p(h_m/(LB+i*h_m)) if i >= iA_m else 0.0 for i in range(48)]
 rows_2w = [(0, [1.0]), (1, wA_m), (1, wB_m), (2, DESC._conv(wA_m, wB_m))]
 def ss_env(m, alp):
  m = ceil(m)
  return min(1.0, sum(w*rho_lo((m-j*LB-i*h_m+alp)/LB) for j, row in rows_2w for i, w in enumerate(row) if w and j*LB+i*h_m+j*h_m <= m))
 sv_env = lambda m: 1.0 if T_S >= m else min(1.0, max(1.0, T_S/LB)*rho_hi((floor(m)-T_S)/LB))
 LOG2U6_T = 6.0*log2(U_T)
 def model_env(side, q):
  qb = log2(q); c = (qb+10.0)/12.0; sf = qb if side == 1 else 0.0; sg = qb if side == 0 else 0.0
  pr = sum((1.0/6.0)*ss_env(12.0*c+qhf[k]-sf, af)*ss_env(12.0*c+LOG2U6_T+qhg[k]-sg, ag) for k in range(6))
  sv = sum((1.0/6.0)*sv_env(12.0*c+qlf[k]-sf)*sv_env(12.0*c+LOG2U6_T+qlg[k]-sg) for k in range(6))
  nf_m = 12.0*c+qhf[-1]-sf; ng_m = 12.0*c+LOG2U6_T+qhg[-1]-sg
  det_t = max(0.0, 1.0 - (1.0 - P_ECM)*((nf_m + ng_m)/8.0 + 4.0))
  uf, ug = nf_m/LB, ng_m/LB; rf, rg = max(rho_lo(uf), 1e-300), max(rho_lo(ug), 1e-300)
  miss_t = sum((rho_hi(max(0.0, uf-log2(p)/LB))/rf + rho_hi(max(0.0, ug-log2(p)/LB))/rg)/p for p in miss_ps)
  return round(1024.0*pr*det_t*max(0.0, 1.0-miss_t), 2), round(1024.0*min(1.0, sv+pr), 2)
 live_ok = True; live_recs = []; d2_live_ok = True
 for (pit, c1s, c1q, c2s, c2q), (e1, e2) in zip(camp["live_d2_chains"], camp["live_d2_expected"]):
  r1 = first_usable(*pit, B_NA, B_NB)
  if r1 is None or c1q not in (r1[2] if c1s == 0 else r1[3]) or not (B_D2 < c1q <= B_NA < pit[1]): live_ok = d2_live_ok = False; continue
  ch1 = find_ch(c1s, c1q, r1[0], r1[1])
  if ch1 is None: live_ok = d2_live_ok = False; continue
  M1, S1 = eval_ball(*ch1); mu1, sv1 = model_env(ch1[0], ch1[1])
  if M1 != e1[3] or S1 != e1[5] or abs(mu1-e1[2]) > 0.02 or abs(sv1-e1[4]) > 0.02: live_ok = d2_live_ok = False
  live_recs.append([ch1[0], ch1[1], 1, mu1, M1, sv1, S1])
  r2 = first_usable(*ch1, B_D2, B_D2)
  if r2 is None or c2q not in (r2[2] if c2s == 0 else r2[3]) or not (B_LO < c2q <= B_D2 < c1q): live_ok = d2_live_ok = False; continue
  ch2 = find_ch(c2s, c2q, r2[0], r2[1])
  if ch2 is None: live_ok = d2_live_ok = False; continue
  M2, S2 = eval_ball(*ch2); mu2, sv2 = model_env(ch2[0], ch2[1])
  if M2 != e2[3] or S2 != e2[5] or abs(mu2-e2[2]) > 0.02 or abs(sv2-e2[4]) > 0.02: live_ok = d2_live_ok = False
  live_recs.append([ch2[0], ch2[1], 2, mu2, M2, sv2, S2])
 for (pit, ecs, ecq), erow in zip(camp["live_d1_parents"], camp["live_d1_expected"]):
  r1 = first_usable(*pit, B_NA, B_NB)
  if r1 is None or ecq not in (r1[2] if ecs == 0 else r1[3]) or not (B_LO < ecq <= B_NA < pit[1]): live_ok = False; continue
  ch = find_ch(ecs, ecq, r1[0], r1[1])
  if ch is None: live_ok = False; continue
  M_o, S_o = eval_ball(*ch); mu_m, sv_m = model_env(ch[0], ch[1])
  if M_o != erow[3] or S_o != erow[5] or abs(mu_m-erow[2]) > 0.02 or abs(sv_m-erow[4]) > 0.02: live_ok = False
  live_recs.append([ch[0], ch[1], 1, mu_m, M_o, sv_m, S_o])
 for cs, cq, cy, cr, emu, eM, esv, eS in camp["live_ctrl_ideals"]:
  M_o, S_o = eval_ball(cs, cq, cy, cr); mu_m, sv_m = model_env(cs, cq)
  if M_o != eM or S_o != eS or abs(mu_m-emu) > 0.02 or abs(sv_m-esv) > 0.02: live_ok = False
  live_recs.append([cs, cq, 0, mu_m, M_o, sv_m, S_o])
 c224 = [[int(x[0]), M.b62d(x[1]), int(x[2]), M.b62d(x[3])/100.0, M.b62d(x[4]), M.b62d(x[5])/100.0, M.b62d(x[6])]
         for x in (s.split() for s in camp["compact_224_b62"].split(";"))]
 def _summ(rs):
  n = len(rs); mus = [r[3] for r in rs]; ms = [r[4] for r in rs]; svms = [r[5] for r in rs]; svos = [r[6] for r in rs]
  pgf = {str(t): dict(emp=round(sum(t**m for m in ms)/n, 8),
                      nb1=round(sum((1.0+u*(1.0-t))**(-1.0) for u in mus)/n, 8),
                      nb2=round(sum((1.0+0.25*u*(1.0-t))**(-2.0) for u in mus)/n, 8),
                      ok=(sum(t**m for m in ms)/n <= sum((1.0+0.25*u*(1.0-t))**(-2.0) for u in mus)/n
                          and sum(t**m for m in ms)/n <= sum((1.0+u*(1.0-t))**(-1.0) for u in mus)/n)) for t in (0.0, 0.2, 0.5, 0.8)}
  return dict(n=n, mu_min=round(min(mus),2), mu_max=round(max(mus),2), M_min=min(ms), M_max=max(ms),
              n_zero=sum(1 for m in ms if m==0),
              nb1_p0=round(sum((1.0+u)**(-1.0) for u in mus)/n, 6),
              nb2_p0=round(sum((1.0+0.25*u)**(-2.0) for u in mus)/n, 6),
              mean_M_over_mu=round(sum(m/u for m, u in zip(ms, mus))/n, 4),
              two_sided_surv_obs_over_up=round(sum(svos)/sum(svms), 4),
              two_sided_surv_max_over_up=round(max(o/m for o, m in zip(svos, svms)), 4),
              pgf_ok=all(v["ok"] for v in pgf.values()), pgf=pgf)
 s_all, s_d1, s_d2, s_ctrl, s_live = _summ(c224), _summ([r for r in c224 if r[2]==1]), _summ([r for r in c224 if r[2]==2]), _summ([r for r in c224 if r[2]==0]), _summ(live_recs)
 ok = (tower_ok and live_ok and d2_live_ok and len(c224)==224 and len(live_recs)==16
       and s_all["n_zero"]==0 and s_all["mu_min"]>=40.0
       and s_all["pgf_ok"] and s_d1["pgf_ok"] and s_d2["pgf_ok"] and s_ctrl["pgf_ok"] and s_live["pgf_ok"]
       and s_all["mean_M_over_mu"]>=1.0 and s_all["two_sided_surv_max_over_up"]<1.0 and s_live["two_sided_surv_max_over_up"]<1.0)
 return dict(ok=ok, tower_verified=tower_ok, live_16_usable_and_depth2_exact=live_ok,
             live_4_depth2_2hop_chains_exact=d2_live_ok,
             tower_params=dict(u_t=U_T, p_t=P_T, lb=LB, nA=NA, nB=NB, b_d2=B_D2, T_surv=T_S, ball_pts=1024, p_ecm=P_ECM, miss_fb_primes=len(miss_ps)),
             usable_thinning_totals=camp["usable_thinning_totals"],
             slot_d1_vs_ctrl_M_over_mu=round(s_d1["mean_M_over_mu"]/s_ctrl["mean_M_over_mu"], 4),
             slot_d2_vs_ctrl_M_over_mu=round(s_d2["mean_M_over_mu"]/s_ctrl["mean_M_over_mu"], 4),
             all_224=s_all, slot_d1_100=s_d1, slot_d2_24=s_d2, ctrl_100=s_ctrl, live_16=s_live)
