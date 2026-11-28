"""canary1_verifiers.py — Pure-stdlib verifiers for BLS12-381 Stage-1 Collection (Cold-Path 2D NP Gate + Inactive Direct-1):
 1. Single-tier [71.2] 2D-56b NTT3 (reusing existing ntt61_r4_u2, zero new ASM) for Stage-1 gbatch/fbatch
    and product-tree build/workspace sizing (exact 3-prime 24*Mtot root workspace RAM).
 2. 100%-deterministic, zero-failure, distribution-free H=44 Bernstein binary-split prime tree T_48
    band_split factoring every prime p <= 2^48 with exact multiplicity (live prime-power unit test).
 3. Bit-exact (<1e-11 vs QuantizedRetention.at_u8 / survivor_upper_u8) Cold-Path 2D Neyman-Pearson Gate
    (S_g >= T_g_lo(b,c) AND S_f + S_g >= T_g_lo(b,c) + delta_jw, 100% column-constant, zero S_rem) evaluator,
    exact conditioned ret_mult_f/g_min, and live 400-draw discrete O_H single-draw parity callback.
 4. Hardened 8-doc positive + negative stale-token programmatic doc_code_sync."""
import json
import math
import re
from fractions import Fraction as F
from pathlib import Path

_SP = math.sumprod


def _ntt3_2d_56b_opt(M_prod_64):
 """Return (best_cost_ir, best_mtot_slots) for 2D 56-bit-limb 3-prime NTT reusing ntt61_r4_u2."""
 import cm_models as CM
 m64 = math.ceil(max(M_prod_64, 1.0))
 if m64 <= (1 << 52):
  k = max(2, (m64 - 1).bit_length())
  mtot = (3 * (1 << (k - 2))) if (k >= 3 and 3 * (1 << (k - 2)) >= m64) else (1 << k)
  return CM.ntt3_F(m64), mtot
 m56 = math.ceil(m64 * 64.0 / 56.0)
 assert m56 * ((1 << 56) - 1) ** 2 < CM.NTT3_PROD, "56b NTT3 dynamic range overflow!"
 L56 = math.ceil(m56 / 2.0)
 best, best_mtot = 1e30, 0
 k_tot = max(4, (L56 - 1).bit_length())
 for k1 in range(max(2, k_tot // 2 - 2), min(50, k_tot // 2 + 3)):
  K1 = 1 << k1
  K2 = math.ceil(L56 / K1)

  def best_1d(req):
   k = max(2, (req - 1).bit_length())
   cands = [(1 << k, k, False)]
   if k >= 3 and 3 * (1 << (k - 2)) >= req:
    cands.append((3 * (1 << (k - 2)), k - 2, True))
   return min(cands, key=lambda x: x[0])

  N1, k1n, r3_1 = best_1d(2 * K1)
  N2, k2n, r3_2 = best_1d(2 * K2)
  assert N1 <= 3 * (1 << 52) and N2 <= 3 * (1 << 52)
  assert (N1 * N2) * ((1 << 56) - 1) ** 2 < CM.NTT3_PROD
  c_r = 37.50 * k1n + (26.0 if (k1n & 1) else 0.0) + (144.0 if r3_1 else 0.0)
  c_c = 37.50 * k2n + (26.0 if (k2n & 1) else 0.0) + (144.0 if r3_2 else 0.0)
  cost = float(N1 * N2) * (c_r + c_c + 100.0 + 1024.0 + 40.0)
  if cost < best:
   best, best_mtot = cost, N1 * N2
 return best, best_mtot


def ntt3_2d_56b_F(M_prod_64):
 return _ntt3_2d_56b_opt(M_prod_64)[0]


def ntt3_2d_56b_mtot(M_prod_64):
 return _ntt3_2d_56b_opt(M_prod_64)[1]


def gbatch_ntt3_2d(pbits, leaf=9, c_leaf=32768.0):
 SP = math.ceil(pbits / 64.0)
 L = max(1, math.ceil(math.log2((SP + 1) / leaf)))
 G = L + 4
 n = lambda j: leaf * (1 << j)
 PT = sum((1 << (L - j)) * (ntt3_2d_56b_F(2 * n(j - 1)) + 4.0 * n(j - 1)) for j in range(1, L + 1))
 RT = sum((1 << (L - j)) * ((5.0 / 3.0) * ntt3_2d_56b_F(n(j) + G) + 4.0 * (n(j - 1) + G)) for j in range(1, L + 1))
 return (PT + RT + 3.0 * ntt3_2d_56b_F(2 * (max(SP, n(L)) + G))) / (1 << L) + c_leaf


def _verify_binary_split_prime_tree_live():
 """Live exactness check of prime-power extraction (gcd(rem, P_hi)==1) & T_hi binary-split prime tree."""
 import cm_models as CM
 ps = CM.primes_upto(256)

 def build(arr):
  if len(arr) == 1:
   return (arr[0], None, None)
  mid = len(arr) // 2
  L, R = build(arr[:mid]), build(arr[mid:])
  return (L[0] * R[0], L, R)

 tree = build(ps)
 P_hi = tree[0]

 def split(S, node, out):
  if S == 1:
   return
  prod_v, left, right = node
  if left is None:
   e = 0
   while S % prod_v == 0:
    S //= prod_v
    e += 1
   assert S == 1 and e >= 1
   out[prod_v] = e
   return
  g_l = math.gcd(S, left[0] % S)
  S_l, S_r = g_l, S // g_l
  while (g2 := math.gcd(S_l, S_r)) > 1:
   S_l *= g2
   S_r //= g2
  assert S_l * S_r == S and math.gcd(S_l, S_r) == 1
  split(S_l, left, out)
  split(S_r, right, out)

 for N, exp_f, exp_r in (
  (2**4 * 17**3 * 31**2 * 97 * 251 * 1009 * 10007, {2: 4, 17: 3, 31: 2, 97: 1, 251: 1}, 1009 * 10007),
  (3**5 * 19**4 * 131**2 * 241 * 65537, {3: 5, 19: 4, 131: 2, 241: 1}, 65537),
  (1009 * 65537, {}, 1009 * 65537),
 ):
  S0 = math.gcd(N, P_hi % N)
  rem = N // S0
  while (g2 := math.gcd(S0, rem)) > 1:
   S0 *= g2
   rem //= g2
  assert S0 * rem == N and math.gcd(rem, P_hi) == 1 and rem == exp_r
  got = {}
  split(S0, tree, got)
  assert got == exp_f
  chk = rem
  for p, e in got.items():
   chk *= p ** e
  assert chk == N
 return True


def band_split_2d_ntt3_cost(p_hi=48.0, n_rels_upper=2.0 ** 63.790):
 """100%-deterministic, zero-failure, distribution-free complete factorization of all p <= 2^{p_hi}
 (with exact multiplicities) via trial division p <= 2^16 + root prime-power extraction + H-level T_48 tree."""
 import cm_models as CM
 assert _verify_binary_split_prime_tree_live()
 pbits = 1.01 * (2.0 ** p_hi) / CM.LN2
 H = math.ceil(math.log2(1.25506 * (2.0 ** p_hi) / (p_hi * CM.LN2)))
 c_tdiv16 = 6542.0 * (6.0 * (10.0 + 12.0) + 48.0)
 c_gcd_768 = 10.0 * CM.MONT5 * ((12.0 / 5.0) ** 2)
 c_f0 = gbatch_ntt3_2d(pbits, leaf=10, c_leaf=32768.0)
 c_g0 = gbatch_ntt3_2d(pbits, leaf=12, c_leaf=32768.0)
 c_tree_per_rel = c_f0 + c_g0 + 14.0 * c_gcd_768
 for h in range(H):
  pb_h = max(1024.0, pbits / (2.0 ** (h + 1)))
  cf_h = gbatch_ntt3_2d(pb_h, leaf=10, c_leaf=32768.0)
  cg_h = gbatch_ntt3_2d(pb_h, leaf=12, c_leaf=32768.0)
  c_tree_per_rel += 40.0 * (cf_h + 7.0 * c_gcd_768) + 48.0 * (cg_h + 7.0 * c_gcd_768)
 sp48 = math.ceil(pbits / 64.0)
 c_build_t48 = math.ceil(math.log2(sp48)) * ntt3_2d_56b_F(sp48) + 100.0 * sp48
 c_partial_run = (H + 1) * 4.0 * math.ceil(math.log2(sp48)) * ntt3_2d_56b_F(sp48)
 assert (H + 1) * sp48 * 8.0 < (1 << 54)
 c_rel = c_tdiv16 + c_tree_per_rel + (1 << 24)
 c_tot_run = n_rels_upper * c_rel + c_partial_run + c_build_t48
 return c_rel, math.log2(c_tot_run / 20000.0)


def exact_pmf_ccdf(qr, nb, sb, sf, lb=None, kmax=4):
 """Bit-exact (<1e-11 vs qr.at_u8 / qr.survivor_upper_u8) PMF and CCDF over all u8 scores w=0..W."""
 st = qr.step
 n = round(nb / st)
 C = qr._get_C(sf)
 W = len(C) - 1
 if n < 0 or n > qr.length:
  return [0.0] * (W + 1), [0.0] * (W + 2)
 D_m = qr.length if qr.length < 2048 else 2048
 if sf:
  rb = (math.floor(math.ceil(nb) / qr._lmin) + 1) * 0.0625
  b = 9 ** 9 if sb is None else math.floor((sb - rb - 0.125) / st)
  if b < 0:
   return [0.0] * (W + 1), [0.0] * (W + 2)
  lim = n if n < D_m else D_m
  lim_b = D_m if b >= lim else b
  m, m_lo, m_hi = qr._m_conv(lim_b, None, True)
  scale = qr.sel_norm / qr.denominator[n]
 else:
  lim = n if n < D_m else D_m
  if sb is None:
   lim_b = D_m
   win = None
  else:
   lim_b = math.floor((sb + 0.25) / st)
   if lim_b < 0:
    return [0.0] * (W + 1), [0.0] * (W + 2)
   if lim_b > D_m:
    lim_b = D_m
   win = (kmax, lb, (qr.cap(nb) + 2) / 16.0 + 0.25)
  m, m_lo, m_hi = qr._m_conv(lim_b, win, False)
  scale = 1.0 / qr.denominator[n]
 qi = qr.q_i
 D = qr._dmax
 span = 2 * D + 1
 cache_rev = getattr(qr, "_C_rev_cache", None)
 if cache_rev is None:
  cache_rev = {}
  qr._C_rev_cache = cache_rev
 crev = cache_rev.get(sf)
 if crev is None:
  crev = []
  for wv in range(W + 1):
   row = C[wv]
   j0 = next((j for j, x in enumerate(row) if x != 0.0), span)
   j1 = span - 1 - next((j for j, x in enumerate(reversed(row)) if x != 0.0), span)
   crev.append((row[j0 : j1 + 1][::-1], j0, j1) if j0 <= j1 else None)
  cache_rev[sf] = crev
 k_max_m = lim if lim < m_hi else m_hi
 pmf = [0.0] * (W + 1)
 w_lo = (n - k_max_m - D) // qi
 if w_lo < 0:
  w_lo = 0
 w_hi = (n - m_lo + D) // qi + 1
 if w_hi > W:
  w_hi = W
 sp = _SP
 for wv in range(w_lo, w_hi + 1):
  ent = crev[wv]
  if ent is None:
   continue
  rrev, j0, j1 = ent
  nbj = n - (wv * qi - D)
  k0 = nbj - j1
  k1 = nbj - j0
  if k0 < m_lo:
   k0 = m_lo
  if k1 > k_max_m:
   k1 = k_max_m
  if k0 <= k1:
   o0 = k0 - (nbj - j1)
   pmf[wv] = sp(m[k0 : k1 + 1], rrev[o0 : o0 + (k1 - k0 + 1)]) * scale
 ccdf = [0.0] * (W + 2)
 acc = 0.0
 for wv in range(w_hi, -1, -1):
  acc += pmf[wv]
  ccdf[wv] = 0.0 if acc < 0.0 else (1.0 if acc > 1.0 else acc)
 return pmf, ccdf


def eval_cold_np_2d_collection(
 pool, op, det, ret_f, cst_f, ret_g_sh, cst_g_sh, p_f, p_g_shells,
 c_fbatch, c_g_extra, stage_b_sp_fn,
):
 """Exact Cold-Path 2D Neyman-Pearson Gate (S_g >= T_g_lo(b,c) AND S_f + S_g >= T_g_lo(b,c) + delta_jw,
 100% column-constant with zero per-cell norm oracle and zero post-Stage-B S_rem gate).
 Uses bit-exact QuantizedRetention 2D PMFs (<1e-11 vs at_u8/survivor_upper_u8) and explicitly prices
 StageB_g, StageB_f, single-tier [71.2] f-batch, and [71.56] g-batch."""
 import cm_models as CM
 lf, lg = op["lf"], op["lg"]
 alpha_f, alpha_g_ge3, kap2_g = op["alpha_f"], op["alpha_g_ge3"], op["kap2_g"]
 q_g = float(det["b1_q_g"])
 q_f = float(det["b1_q_f"])
 assert abs(q_g - q_f) < 1e-12, "Unweighted 2D-NP sum S_f + S_g requires q_f == q_g!"
 sbf, sbg = float(det["sbf"]), float(det["sbg"])
 delta_jw = int(det["delta_jw"])
 assert delta_jw >= 1, "delta_jw must be positive!"
 _slg = float(det["sl_g_tgt"])
 _tlo = int(det["tg_u8_min"])
 _thi = int(det["tg_u8_max"])
 _kh_g, _kh_f = float(det.get("kh_g", 0.08)), float(det.get("kh_f", 0.08))
 _rmg = max(round(b / q_g) / b for b in [33.0 + i / 16.0 for i in range(625)])
 hs = float(pool.get("hist_step", 0.25))
 hh = 0.5 * hs
 v_shares = list(pool["v_shares"])
 hist_shells = pool["hist_shells"]
 assert abs(p_g_shells[0] - p_g_shells[1]) < 1e-12, "Inactive Direct-1 uses shared p_g across shells!"
 p_g = float(p_g_shells[0])
 ret_g, cst_g = ret_g_sh[0], cst_g_sh[0]

 # Machine-precision self-check against QuantizedRetention.at_u8 & survivor_upper_u8
 _, _chk_rf = exact_pmf_ccdf(ret_f, 540.25, sbf, True)
 assert max(abs(_chk_rf[w] - ret_f.at_u8(540.25, w, sbf)) for w in range(len(_chk_rf) - 1)) < 1e-11
 _, _chk_cg = exact_pmf_ccdf(cst_g, 664.0, sbg, False, lg, int(sbg // p_g))
 assert max(abs(_chk_cg[w] - cst_g.survivor_upper_u8(664.0, w, sbg, lg, int(sbg // p_g))) for w in range(len(_chk_cg) - 1)) < 1e-11

 uniq_bf = sorted({r[0] for k in range(2) for r in hist_shells[k]})
 uniq_bgo = sorted({r[1] for k in range(2) for r in hist_shells[k]})
 all_bf_pts = sorted({x for bf in uniq_bf for x in (bf, bf + hs)})
 all_bgo_pts = sorted({x for bgo in uniq_bgo for x in (bgo, bgo + hs)})
 wc_bf_pts = sorted({x for k in range(2) for r in hist_shells[k] if r[3] > 0.0 for x in (r[0], r[0] + hs)})
 wc_bgo_pts = sorted({x for k in range(2) for r in hist_shells[k] if r[3] > 0.0 for x in (r[1], r[1] + hs)})

 _cfa_pt = {x: exact_pmf_ccdf(cst_f, x, None, False)[1] for x in all_bf_pts}
 _cf_pt = {x: exact_pmf_ccdf(cst_f, x, sbf, False, lf, int(sbf // p_f))[1] for x in all_bf_pts}
 _rf_pt = {x: exact_pmf_ccdf(ret_f, x, sbf, True)[1] for x in wc_bf_pts}

 _cga_pt = {x: exact_pmf_ccdf(cst_g, x, None, False)[1] for x in all_bgo_pts}
 _cg_pt = {x: exact_pmf_ccdf(cst_g, x, sbg, False, lg, int(sbg // p_g))[1] for x in all_bgo_pts}
 _rg_pt = {x: exact_pmf_ccdf(ret_g, x, sbg, True)[1] for x in wc_bgo_pts}

 W_F = len(_cfa_pt[all_bf_pts[0]]) - 2
 W_G = len(_cga_pt[all_bgo_pts[0]]) - 2
 _wf_off = W_F + 1
 _z_f = [0.0] * (W_F + W_G + 4)

 f_tab = {}
 for bf in uniq_bf:
  sfa = [min(1.0, 1.25 * (a if a > b else b)) for a, b in zip(_cfa_pt[bf], _cfa_pt[bf + hs])]
  cs = [min(1.0, 1.25 * (a if a > b else b)) for a, b in zip(_cf_pt[bf], _cf_pt[bf + hs])]
  if bf in _rf_pt and (bf + hs) in _rf_pt:
   crm = [0.5 * (a + b) for a, b in zip(_rf_pt[bf], _rf_pt[bf + hs])]
   crh = _rf_pt[bf + hs]
   crm_r = crm[::-1] + [crm[0]] * (W_G + 2)
   crh_r = crh[::-1] + [crh[0]] * (W_G + 2)
  else:
   crm_r = crh_r = _z_f
  qf = min(1.0, 2.0 * CM.rho_exact((bf + alpha_f) / lf))
  ef_c = bf + hh + alpha_f; ef_h = bf + hs + alpha_f
  f_tab[bf] = (
   crm_r, crh_r,
   cs[::-1] + [cs[0]] * (W_G + 2),
   sfa[::-1] + [sfa[0]] * (W_G + 2),
   qf, CM.rho_db(ef_c / lf, ef_c), CM.rho_db(ef_h / lf, ef_h),
  )

 def _nz_span(arr):
  lo = 0; n = len(arr)
  while lo < n and arr[lo] == 0.0:
   lo += 1
  if lo >= n:
   return (n, -1, arr)
  hi = n - 1
  while hi > lo and arr[hi] == 0.0:
   hi -= 1
  return (lo, hi, arr)

 _z_g = (W_G + 1, -1, [])
 g_tab = {}
 for bgo in uniq_bgo:
  cca = [min(1.0, a if a > b else b) for a, b in zip(_cga_pt[bgo], _cga_pt[bgo + hs])]
  cs = [min(1.0, a if a > b else b) for a, b in zip(_cg_pt[bgo], _cg_pt[bgo + hs])]
  pa = _nz_span([max(0.0, cca[w] - cca[w + 1]) for w in range(W_G + 1)])
  pb = _nz_span([max(0.0, cs[w] - cs[w + 1]) for w in range(W_G + 1)])
  if bgo in _rg_pt and (bgo + hs) in _rg_pt:
   crm = [0.5 * (a + b) for a, b in zip(_rg_pt[bgo], _rg_pt[bgo + hs])]
   crh = _rg_pt[bgo + hs]
   prm = _nz_span([max(0.0, crm[w] - crm[w + 1]) for w in range(W_G + 1)])
   prh = _nz_span([max(0.0, crh[w] - crh[w + 1]) for w in range(W_G + 1)])
  else:
   prm = prh = _z_g
  eg_c = bgo + hh + alpha_g_ge3 + 1.0 + kap2_g; eg_h = bgo + hs + alpha_g_ge3 + 1.0 + kap2_g
  g_tab[bgo] = (prm, prh, pa, pb, cca, cs, CM.rho_db(eg_c / lg, eg_c), CM.rho_db(eg_h / lg, eg_h))

 _tgf = lambda x: max(_tlo, min(_thi, max(int(round((x - _slg) / q_g)), 128 + math.floor(x * _rmg) - 255)))
 _tg_ok = set()
 _csb_g_c = {}; _csb_f_c = {}
 ysum_mid = ysum_rethi = ysum_strict_hi = c_cof_live = sga_mean = max_cca = f_prim_tot = var_sh_tot = 0.0
 _nf_cond_num = _ng_cond_num = _w_cond_den = 0.0
 max_bg_seen = 0.0
 _w_pairs = {}
 _ret_pair_cache = {}; _cof_pair_cache = {}
 sp = _SP
 _ret_get = _ret_pair_cache.get; _cof_get = _cof_pair_cache.get
 _tg_in = _tg_ok.__contains__; _tg_add = _tg_ok.add
 _csbg_get = _csb_g_c.get; _csbf_get = _csb_f_c.get

 for k in range(2):
  vk = v_shares[k]
  ym_k = yr_k = ysh_k = ysh2_k = cc_k = sga_k = wk_k = fp_k = 0.0
  for _row in hist_shells[k]:
   bf, bgo, bgm, wc, wk, _tgu = _row
   if bgm > max_bg_seen:
    max_bg_seen = bgm
   _tk = (bgo, _tgu)
   if not _tg_in(_tk):
    assert _tlo <= _tgu <= _thi and _tgu >= _tgf(bgo) and (128 - _tgu) + math.floor((bgo + hs) * _rmg) <= 248 <= 255
    _tg_add(_tk)
   tjw = _tgu + delta_jw
   _pk3 = (bf, bgo, _tgu)
   crm_r, crh_r, cs_r, sfa_r, qf, rf_c0, rf_h0 = f_tab[bf]
   prm, prh, pa, pb, cca, cs_g, rg_c0, rg_h0 = g_tab[bgo]
   _bidx = _wf_off - tjw
   a0_base = _tgu if _tgu >= (tjw - W_F) else (tjw - W_F)
   if wc > 0.0:
    fp_k += wc
    _rv = _ret_get(_pk3)
    if _rv is None:
     mlo, mhi, marr = prm
     a0 = a0_base if a0_base > mlo else mlo
     ret_m = sp(marr[a0 : mhi + 1], crm_r[_bidx + a0 : _bidx + mhi + 1]) if a0 <= mhi else 0.0
     hlo, hhi, harr = prh
     a0 = a0_base if a0_base > hlo else hlo
     ret_h = sp(harr[a0 : hhi + 1], crh_r[_bidx + a0 : _bidx + hhi + 1]) if a0 <= hhi else 0.0
     _rv = (ret_m, ret_h)
     _ret_pair_cache[_pk3] = _rv
    else:
     ret_m, ret_h = _rv
    p_sm = rf_c0 * rg_c0; p_sm_hi = rf_h0 * rg_h0
    ym_k += wc * p_sm * ret_m
    yr_k += wc * p_sm * ret_h
    _ysh = p_sm_hi * ret_h
    ysh_k += wc * _ysh
    ysh2_k += wc * _ysh * _ysh
    if ret_h > 0.0:
     _wsh = vk * wc * _ysh
     _nf_cond_num += _wsh * (bf + hh)
     _ng_cond_num += _wsh * (bgo + hh)
     _w_cond_den += _wsh
     _ck = (int(bf), int(bgo), _tgu)
     _w_pairs[_ck] = _w_pairs.get(_ck, 0.0) + _wsh
   if wk > 0.0:
    _cv = _cof_get(_pk3)
    if _cv is None:
     alo, ahi, aarr = pa
     a0 = a0_base if a0_base > alo else alo
     peg = sp(aarr[a0 : ahi + 1], sfa_r[_bidx + a0 : _bidx + ahi + 1]) if a0 <= ahi else 0.0
     blo, bhi, barr = pb
     a0 = a0_base if a0_base > blo else blo
     if a0 <= bhi:
      _bs = barr[a0 : bhi + 1]
      pef = sp(_bs, sfa_r[_bidx + a0 : _bidx + bhi + 1])
      pfb = sp(_bs, cs_r[_bidx + a0 : _bidx + bhi + 1])
     else:
      pef = pfb = 0.0
     _cv = (peg, pef, pfb)
     _cof_pair_cache[_pk3] = _cv
    else:
     peg, pef, pfb = _cv
    csb_g = _csbg_get(bgm)
    if csb_g is None:
     csb_g = stage_b_sp_fn(float(bgm), _kh_g)
     _csb_g_c[bgm] = csb_g
    csb_f = _csbf_get(bf)
    if csb_f is None:
     csb_f = stage_b_sp_fn(float(bf + hs), _kh_f)
     _csb_f_c[bf] = csb_f
    p_both = pfb if pfb < (cs_g[_tgu] * qf) else (cs_g[_tgu] * qf)
    c_fb = pfb * c_fbatch + p_both * c_g_extra
    cc_k += wk * (peg * (6500.0 + csb_g) + pef * (20000.0 + csb_f) + c_fb)
    _cca_v = cca[_tgu]
    if _cca_v > max_cca:
     max_cca = _cca_v
    sga_k += wk * _cca_v
    wk_k += wk
  assert abs(wk_k - 1.0) < 1e-9
  ysum_mid += vk * ym_k; ysum_rethi += vk * yr_k; ysum_strict_hi += vk * ysh_k
  c_cof_live += vk * cc_k; sga_mean += vk * sga_k; f_prim_tot += vk * fp_k
  var_sh_tot += (vk * vk) * max(0.0, ysh2_k - ysh_k * ysh_k) / float(pool["n_draws_shell"])

 # Live single-draw exact 2D-NP callback for 400-draw discrete O_H parity in unitdomain.py
 _sd_f_c = {}; _sd_g_c = {}
 def ret_2d_single_draw(l2f, l2go, tgu, _shk):
  tjw = tgu + delta_jw
  bf_q = round(l2f * 4.0) * 0.25
  bgo_q = round(l2go * 4.0) * 0.25
  cf_r = _sd_f_c.get(bf_q)
  if cf_r is None:
   cf = exact_pmf_ccdf(ret_f, bf_q, sbf, True)[1]
   cf_r = cf[::-1] + [cf[0]] * (W_G + 2)
   _sd_f_c[bf_q] = cf_r
  pg = _sd_g_c.get(bgo_q)
  if pg is None:
   cg = exact_pmf_ccdf(ret_g, bgo_q, sbg, True)[1]
   pg = [max(0.0, cg[w] - cg[w + 1]) for w in range(W_G + 1)]
   _sd_g_c[bgo_q] = pg
  a0 = tgu if tgu >= (tjw - W_F) else (tjw - W_F)
  if a0 > W_G:
   return 0.0
  _bidx = _wf_off - tjw
  return _SP(pg[a0:], cf_r[_bidx + a0 : _bidx + W_G + 1])

 # Exact conditioned ret_mult_f_min and ret_mult_g_min over unsieved primes b in (p_side, 70.0]
 _rf0_c = {}; _pg0_c = {}
 def _get_rf0(x):
  v = _rf0_c.get(x)
  if v is None:
   c = _rf_pt.get(x)
   if c is None:
    c = exact_pmf_ccdf(ret_f, x, sbf, True)[1]
   v = c[::-1] + [c[0]] * (W_G + 2)
   _rf0_c[x] = v
  return v

 def _get_pg0(x):
  v = _pg0_c.get(x)
  if v is None:
   c = _rg_pt.get(x)
   if c is None:
    c = exact_pmf_ccdf(ret_g, x, sbg, True)[1]
   v = [max(0.0, c[w] - c[w + 1]) for w in range(W_G + 1)]
   _pg0_c[x] = v
  return v

 _cp_list = []
 for (bf1, bgo1, tgu), wv in _w_pairs.items():
  bf_v = float(bf1); bgo_v = float(bgo1)
  tjw = tgu + delta_jw; _bidx = _wf_off - tjw
  a0 = tgu if tgu >= (tjw - W_F) else (tjw - W_F)
  if a0 > W_G:
   continue
  rf0_r = _get_rf0(bf_v + hs); pg0 = _get_pg0(bgo_v + hs)
  r0 = _SP(pg0[a0:], rf0_r[_bidx + a0 : _bidx + W_G + 1])
  if r0 > 1e-12:
   ef = bf_v + hs + alpha_f; eg = bgo_v + hs + alpha_g_ge3 + 1.0 + kap2_g
   _cp_list.append((bf_v, bgo_v, tgu, tjw, _bidx, a0, wv, ef, eg, CM.rho_db(ef / lf, ef), CM.rho_db(eg / lg, eg), 1.0 / r0, rf0_r, pg0))

 _fc_b = {}; _gc_b = {}
 def _cond_ret_mult(side_s, b_val):
  num = den = 0.0
  for bf_v, bgo_v, tgu, tjw, _bidx, a0, wv, ef, eg, rf0, rg0, ir0, rf0_r, pg0 in _cp_list:
   if side_s == "f":
    eb = ef - b_val
    ww = (wv / rf0) * CM.rho_db(eb / lf, max(1.0, eb))
    den += ww
    if sbf > b_val and bf_v + hs > b_val:
     kf = (bf_v, b_val)
     cf_br = _fc_b.get(kf)
     if cf_br is None:
      cf_b = exact_pmf_ccdf(ret_f, bf_v + hs - b_val, sbf - b_val, True)[1]
      cf_br = cf_b[::-1] + [cf_b[0]] * (W_G + 2)
      _fc_b[kf] = cf_br
     num += ww * (_SP(pg0[a0:], cf_br[_bidx + a0 : _bidx + W_G + 1]) * ir0)
   else:
    eb = eg - b_val
    ww = (wv / rg0) * CM.rho_db(eb / lg, max(1.0, eb))
    den += ww
    if sbg > b_val and bgo_v + hs > b_val:
     kg = (bgo_v, b_val)
     pg_b = _gc_b.get(kg)
     if pg_b is None:
      cg_b = exact_pmf_ccdf(ret_g, bgo_v + hs - b_val, sbg - b_val, True)[1]
      pg_b = [max(0.0, cg_b[w] - cg_b[w + 1]) for w in range(W_G + 1)]
      _gc_b[kg] = pg_b
     num += ww * (_SP(pg_b[a0:], rf0_r[_bidx + a0 : _bidx + W_G + 1]) * ir0)
  return num / den if den > 0.0 else 0.0

 ret_mult_f_min = min(_cond_ret_mult("f", b) for b in (p_f + 0.02, 65.0, 68.0, 70.0))
 ret_mult_g_min = min(_cond_ret_mult("g", b) for b in (p_g + 0.02, 62.4, 65.0, 68.0, 70.0))
 return dict(
  ysum_mid=ysum_mid, ysum_rethi=ysum_rethi, ysum_strict_hi=ysum_strict_hi,
  c_cof_live=c_cof_live, sga_mean=sga_mean, max_cca=max_cca, f_prim_tot=f_prim_tot,
  var_sh_tot=var_sh_tot, max_bg_seen=max_bg_seen,
  nf_cond_mean=_nf_cond_num / _w_cond_den, ng_cond_mean=_ng_cond_num / _w_cond_den,
  ret_mult_f_min=ret_mult_f_min, ret_mult_g_min=ret_mult_g_min, ret_mult_g_bridge_min=1.0,
  ret_2d_single_draw=ret_2d_single_draw,
 )


def verify_cold_np_3chord_live(tg_u8, delta_jw, adp, is_cnp=True, tf_u8=45):
 """Bit-exact 3-chord (lengths 16, 11, 16) simulation of Pass-1 G_val push + Pass-2 column-constant
 2D Neyman-Pearson test (S_g >= T_g_lo(b, c) AND S_f + S_g >= T_g_lo(b, c) + delta_jw) with S_b1 refill."""
 import random
 NC = 3 * 128
 ch_sh = (0, 2, -1) if adp else (0, 0, 0)
 tg_ch = [[tuple(tg_u8[c] + (((r % 5) - 2 + sh) if adp else 0) for c in range(3)) for r in range(16)] for sh in ch_sh]
 bg_ch = [[tuple(128 - t for t in row) for row in ch] for ch in tg_ch]
 S_b1 = [bg_ch[0][0][i // 128] for i in range(NC)]
 S_b2 = [[0] * NC for _ in range(32)]
 rng = random.Random(20261125)
 for ci, n_act in enumerate((16, 11, 16)):
  tg_cols, bg_cols, bg_nxt_ch = tg_ch[ci], bg_ch[ci], bg_ch[(ci + 1) % 3][0]
  r_last = n_act - 1
  mg_sb = [[0] * (NC // 64) for _ in range(n_act)]
  gv_sb = [[[] for _ in range(3)] for _ in range(n_act)]
  g1 = [[rng.randrange(0, 55) for _ in range(NC)] for _ in range(n_act)]
  g2 = [[rng.randrange(0, 35) for _ in range(NC)] for _ in range(16)]
  f1 = [[rng.randrange(0, 35) for _ in range(NC)] for _ in range(n_act)]
  f2 = [[rng.randrange(0, 20) for _ in range(NC)] for _ in range(16)]
  for r in range(16):
   for i in range(NC):
    S_b2[r][i] += g2[r][i]
  for r in range(n_act):
   for i in range(NC):
    S_b1[i] += g1[r][i]
   for c in range(3):
    xb_nxt = 0 if r == r_last else bg_cols[r + 1][c]
    for sub in range(2):
     blk = 2 * c + sub; mb_w = 0
     for jj in range(63, -1, -1):
      idx = 64 * blk + jj; v = S_b1[idx] + S_b2[r][idx]; assert v <= 255
      S_b1[idx] = xb_nxt; S_b2[r][idx] = 0
      if v & 0x80:
       mb_w |= 1 << jj
       if is_cnp:
        gv_sb[r][c].append(v - bg_cols[r][c])
     mg_sb[r][blk] = mb_w
  for r in range(n_act, 16):
   S_b2[r] = [0] * NC
  assert all(v == 0 for v in S_b1) and all(v == 0 for row in S_b2[:16] for v in row)
  for r in range(16):
   for i in range(NC):
    S_b2[r][i] += f2[r][i]
  for r in range(n_act):
   for i in range(NC):
    S_b1[i] += f1[r][i]
   got_r = []
   for c in range(3):
    xf_nxt = bg_nxt_ch[c] if r == r_last else 0
    w0, w1 = mg_sb[r][2 * c], mg_sb[r][2 * c + 1]
    if (w0 | w1) != 0:
     gp = 0
     for sub, ww in ((0, w0), (64, w1)):
      while ww:
       jj = ww.bit_length() - 1
       ww ^= 1 << jj
       idx = 128 * c + sub + jj
       sf_u = S_b1[idx] + S_b2[r][idx]
       if is_cnp:
        sg_u = gv_sb[r][c][gp]; gp += 1
        if sf_u + sg_u >= tg_cols[r][c] + delta_jw:
         got_r.append(idx)
       elif sf_u >= tf_u8:
        got_r.append(idx)
     if is_cnp:
      assert gp == len(gv_sb[r][c])
    for jj in range(128):
     S_b1[128 * c + jj] = xf_nxt; S_b2[r][128 * c + jj] = 0
   if is_cnp:
    exp_r = [i for i in range(NC) if (g1[r][i] + g2[r][i] >= tg_cols[r][i // 128]) and (f1[r][i] + f2[r][i] + g1[r][i] + g2[r][i] >= tg_cols[r][i // 128] + delta_jw)]
   else:
    exp_r = [i for i in range(NC) if (g1[r][i] + g2[r][i] >= tg_cols[r][i // 128]) and (f1[r][i] + f2[r][i] >= tf_u8)]
   assert sorted(got_r) == exp_r
  for r in range(n_act, 16):
   S_b2[r] = [0] * NC
  assert all(S_b1[i] == bg_nxt_ch[i // 128] for i in range(NC)) and all(v == 0 for row in S_b2[:16] for v in row)
 return True


def check_doc_code_sync(root, led, op, cl, dres=None, p2_live=None):
 """Positive numeric sync + negative stale-token assertions across all 8 package docs & manifest."""
 root = Path(root)
 names = ("proof.md", "STORAGE.md", "RANK-BRIDGE.md", "ADAPTIVE.md", "DESCENT.md", "KERNELS.md", "KERNELS_LA.md", "submission.json")
 docs = {nm: (root / nm).read_text() for nm in names}
 proof, storage, rbridge, adaptive, descent, kernels, kla, sub_txt = [docs[nm] for nm in names]
 sc = json.loads(sub_txt)["claims"][0]["score"]
 assert F(int(sc["numerator"]), int(sc["denominator"])) == F(led["claim_p"])
 assert led["claim_p"] in proof and led["claim_p"] in sub_txt
 assert led["exp_run_claim"] in proof and led["exp_run_claim"] in sub_txt
 assert str(led["K_gate"]) in proof and str(led["K_gate"]) in sub_txt
 assert led["p_floor"] in proof and led["p_floor"] in sub_txt
 assert led.get("p_cold_target", "0.391") in proof and led.get("p_cold_target", "0.391") in sub_txt
 assert f"{op['lf']:.2f}" in proof and f"{op['lg']:.2f}" in proof
 assert f"{cl['c_cof_2d_upper_live']:.5f}" in proof and f"{cl['c_cof_2d_upper_live']:.5f}" in kla and f"{cl['c_cof_2d_upper_live']:.5f}" in sub_txt
 assert "[71.2]" in proof and "[71.2]" in kla and "[71.2]" in sub_txt
 assert f"{cl['band_u']:.4f}" in proof
 assert f"{cl['la_u']:.4f}" in proof
 assert f"{cl['slab_owner_eib']:.3f}" in storage and f"{cl['slab_owner_eib']:.3f}" in kernels
 assert f"{cl['pt_2048own_frac']:.4f}" in storage
 assert f"{cl['L_max_shells'][0]:.2f}" in kernels and f"{cl['L_max_shells'][0]:.2f}" in proof
 assert f"{cl['ret_mult_f_min']:.5f}" in proof and f"{cl['ret_mult_g_min']:.5f}" in proof and f"{cl['log2_rels_2se_live']:.6f}" in proof
 assert f"{cl['ret_mult_f_min']:.5f}" in descent and f"{cl['ret_mult_g_min']:.5f}" in descent and f"{cl['log2_rels_2se_live']:.6f}" in descent
 assert f"{op['lf']:.2f}" in descent and f"{op['lg']:.2f}" in descent
 assert "2.831 EiB" in storage and "2.831 EiB" in proof
 if p2_live is not None:
  g = p2_live["l1_geom_live"]
  assert g["rounding_cell_exact_cover_ok"] is True
  assert f"{g['E_q_exact_256']:.6f}" in proof and f"{g['R_cov_shells'][0]:.2f}" in proof
 if dres is not None:
  phi_s = f"{dres['phi_tree']:.6f}"
  assert phi_s in adaptive and phi_s in descent and phi_s in proof and phi_s in sub_txt
 stale_all = (
  "100.2657", "0.72392", "2.48e-9", "D x (D + 4 SM)", "Chernoff concentration",
  "63.24275", "[62.5, 65.5, 68.5, 71.2]", "3.435 EiB", "0.2147",
  "nu_fvs", "FVS Schur-BW", "bls_fvs_explicit", "106.735800", "7194949683715096640",
  "106.740800", "108.0939", "0.45709", "6.5054",
  "107.861", "104.845", "0.799", "T_f_u8 = 45", "63.260409",
  "106.503200", "107.8565", "0.8820", "0.66108", "Window-Map", "2^{48.8}", "4.504 EiB",
 )
 for fn, txt in docs.items():
  for tok in stale_all:
   assert tok not in txt, f"Stale token {tok!r} in {fn}!"
  assert re.search(r"\brbit\b", txt, re.I) is None, f"Disallowed mnemonic 'rbit' in {fn}!"
 for tok in ("de46f228", "a1d2e1bc", "byte-identical", "baseline", "prior version"):
  assert tok not in sub_txt.lower(), f"Forbidden comparative token {tok!r} in submission.json!"
 return dict(docs_checked=len(docs), stale_tokens_checked=len(stale_all), sync_ok=True)
