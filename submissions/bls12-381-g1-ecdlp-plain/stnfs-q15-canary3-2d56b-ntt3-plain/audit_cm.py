"""audit_cm.py — Self-contained pure-Python-stdlib hosted Docker audit for BLS12-381 CM-STNFS over Q(zeta_15)
(K=2 F_Vor(a) x L_{3/4}(b), M=15 slabs, per-column odd-norm adaptive T_g(b,c), split p_rot=56.5 4-rot W_B=1,
single-tier [71.2] 2D-56b NTT3 f-batch/g-batch & deterministic T_48 binary-split band_split,
E_16 cap C_30 cap F_Vor(a,b) Gram-ellipsoid descent, and superadditive Kronecker-NTT3 + BLS Stage-2 ECM)."""
import base64
import json
import lzma
import math
import sys
from decimal import Decimal, ROUND_CEILING, localcontext
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE / "files" if (HERE / "files").is_dir() else HERE
if str(ROOT) not in sys.path:
 sys.path.insert(0, str(ROOT))
if str(HERE) not in sys.path:
 sys.path.insert(0, str(HERE))

from audit_algebra import (
 EXPS_1W, M_1W, M_1W_INV, P_BLS, PHI15, R_BLS, U_BLS,
 add_phi15_z, build_G1w, build_WS_1w, check_algebra_bls12, check_toy_e2e_dlp,
 check_voronoi_exact, int_det, inv_int_mat, mat_T, mat_mul, mul_phi15_z,
 norm_q15_z, pgcd_mod, pmul_mod, pmul_z, pow_phi15_z, ppow_mod, prem_mod,
 ptrim, red_phi15_z, scale_phi15_z, sub_phi15_z, trace_z15, y_pow_z,
)
from cm_models import (
 GAMMA, LN2, MONT5, NTT3_PRIMES, NTT3_PROD, NTT4_PRIMES, NTT4_PROD,
 QuantizedRetention, Retention, ecm_instr, ecm_mean_var, fft_F, gbatch_ntt3,
 lll_16d_int, log2_li, missing_recon_z12, mr_instr, ntt3_F, pbits_split,
 _ret_coeffs, primes_upto, rho_db, rho_exact, srt_toy_check, verify_kernels_symbolic,
)

# --- Exact 4-band QuantizedRetention DP for deterministic 4-tier f-batch (merged from fbatch_tiers) ---
import operator as _OP
_SUMPROD = getattr(math, "sumprod", lambda a, b, _m=_OP.mul: sum(map(_m, a, b)))


def _conv_sp(a, b, D):
 """Exact discrete convolution c[s] = sum_{i+j=s} a[i]*b[j] for 0 <= s <= D via math.sumprod."""
 a_lo = next((i for i, x in enumerate(a) if x != 0.0), len(a))
 b_lo = next((i for i, x in enumerate(b) if x != 0.0), len(b))
 if a_lo >= len(a) or b_lo >= len(b):
  return [0.0] * (D + 1)
 a_hi = len(a) - 1 - next(i for i, x in enumerate(reversed(a)) if x != 0.0)
 b_hi = len(b) - 1 - next(i for i, x in enumerate(reversed(b)) if x != 0.0)
 lb = b_hi + 1
 b_rev = b[:lb][::-1]
 out = [0.0] * (D + 1)
 sp = _SUMPROD
 s_end = min(D, a_hi + b_hi)
 for s in range(a_lo + b_lo, s_end + 1):
  i0 = s - b_hi if s - b_hi > a_lo else a_lo
  i1 = s - b_lo if s - b_lo < a_hi else a_hi
  if i0 <= i1:
   off = lb - 1 - s
   out[s] = sp(a[i0 : i1 + 1], b_rev[off + i0 : off + i1 + 1])
 return out


def build_fbatch_tier_m(cst_f, p_f, b_mid, sbf, lf):
 """Precompute 4-band convolution table (small * mid[p_f..b_mid] * large_all[b_mid..cst_f.B=660])
    on the exact QuantizedRetention grid of cst_f (cst_f.B = 660 >= 660 integrates ALL primes up to
    660 bits — including every single prime R_t in [125, 362] bits — NOT capped at lf = 71.2!),
    asserting algebraic splitting identity mid * large_all == large_geo to < 1e-6!"""
 assert cst_f.B >= 660, f"cst_f.B={cst_f.B} < 660: must integrate all primes up to full f-norm!"
 st = cst_f.step
 D = cst_f.length if cst_f.length < 2048 else 2048
 tilt = cst_f.tilt
 w_mid = [0.0] * (D + 1)
 w_lg = [0.0] * (D + 1)
 t_mid = 0.0
 t_lg = 0.0
 w0 = st / 4.0
 for left, right, is_mid in ((p_f, b_mid, True), (b_mid, cst_f.B, False)):
  if left >= right:
   continue
  cells = math.ceil((right - left) / w0)
  w = (right - left) / cells
  for i in range(cells):
   b = left + (i + 0.5) * w
   bi = round(b / st)
   idx = bi
   for k in range(1, 9):
    ex = k * tilt * bi * st - (k - 1) * LN2 * b
    if ex < -75.0 and idx > D:
     break
    amt = math.exp(ex) * w / (k * b)
    if is_mid:
     t_mid += amt
     if idx <= D:
      w_mid[idx] += amt
    else:
     t_lg += amt
     if idx <= D:
      w_lg[idx] += amt
    idx += bi
 c_mid = _ret_coeffs(w_mid, t_mid, D)
 c_lg = _ret_coeffs(w_lg, t_lg, D)
 c_sm = cst_f.co("small")
 c_lgeo = cst_f.co("large_geo")
 # Verify exact algebraic splitting identity: (mid * lg_all)[s] == large_geo[s]!
 mid_lg_all = _conv_sp(c_mid, c_lg, D)
 max_id_err = max(abs(mid_lg_all[s] - c_lgeo[s]) for s in range(D + 1))
 assert max_id_err < 1e-6, f"4-band f-batch splitting identity error {max_id_err:.2e} >= 1e-6!"

 # Gate ONLY the exact known post-Tier-t integer R_t: drop R_t == 1 (already smooth),
 # R_t > 2^{sbf} (impossible for Stage-B R_0 <= 2^{sbf}), or R_t outside union_{k>=2} [b_mid*k, lf*k].
 # ZERO truncation on i + j (sum of stripped mid primes + R_t runs across full D = 2048 >= u_lim)!
 w2 = (cst_f.cap(660.0) + 2) / 16.0 + 0.25
 lim_lg = min(D, math.floor((sbf + w2) / st))
 kmax = int(sbf // b_mid)
 lg_win = [
  v if (1 <= l <= lim_lg and v and any(b_mid * k - w2 <= l * st <= lf * k + w2 for k in range(2, kmax + 1))) else 0.0
  for l, v in enumerate(c_lg[: lim_lg + 1])
 ]
 m_ml_f = _conv_sp(c_mid, lg_win, D)
 m_tot = _conv_sp(m_ml_f, c_sm, D)
 m_lo = next((i for i, x in enumerate(m_tot) if x != 0.0), D + 1)
 m_hi = D - next((i for i, x in enumerate(reversed(m_tot)) if x != 0.0), D + 1)
 return m_tot, m_lo, m_hi, max_id_err


def p_fbatch_tier_u8(cst_f, m_info, nb, T_u8):
 """Exact 4-band survivor probability for entering the next f-batch tier (untruncated i+j <= D)."""
 m_tot, m_lo, m_hi, _ = m_info
 n = round(nb / cst_f.step)
 if n < 0 or n > cst_f.length:
  return 0.0
 jmin = max(0, int(T_u8) * cst_f.q_i - cst_f._dmax)
 u_lim = n - jmin
 if u_lim < 0:
  return 0.0
 D = cst_f.length if cst_f.length < 2048 else 2048
 assert u_lim <= D, f"p_fbatch_tier_u8 truncation: u_lim={u_lim} > D={D} at nb={nb}"
 lim = min(n, D, u_lim)
 sv_r, s_lo, s_hi = cst_f._spass_rev(T_u8, False)
 off = cst_f.length - n
 lo_k = max(m_lo, n - s_hi)
 hi_k = min(lim, m_hi, n - s_lo)
 if lo_k > hi_k:
  return 0.0
 r = _SUMPROD(m_tot[lo_k : hi_k + 1], sv_r[off + lo_k : off + hi_k + 1]) / cst_f.denominator[n]
 return 0.0 if r < 0.0 else (1.0 if r > 1.0 else r)
from filtsim_z12 import verify_h3_ladder_and_fit
from hx_z12 import run_hx_on_filtered_z12, unplanted_two_sided_bw_check
from w8_verifiers import check_sigmaB_exactness_live, verify_genuine_q15_filter_hx_live
from canary1_verifiers import (
 ntt3_2d_56b_F, ntt3_2d_56b_mtot, gbatch_ntt3_2d, band_split_2d_ntt3_cost, check_doc_code_sync,
)

UNIT = 20000.0
_RCTX = []


def check_sigmaB_exactness():
 return check_sigmaB_exactness_live(build_G1w, build_WS_1w, mat_mul, M_1W, norm_q15_z, U_BLS)


def check_anchor_norms(anchors):
 import unitdomain as _UD
 return _UD.check_anchor_norms(
  anchors, M_1W, U_BLS, y_pow_z, add_phi15_z, sub_phi15_z,
  mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z,
 )


def verify_pool_2adic_live(pool, op, run_ht_mc=True, ht_only=False, ret_ctx=None):
 import unitdomain as _UD
 return _UD.verify_pool_2adic_live(
  pool, op, build_G1w, build_WS_1w, mat_mul, M_1W, PHI15, U_BLS,
  y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z,
  scale_phi15_z, norm_q15_z, pgcd_mod, rho_db,
  run_ht_mc=run_ht_mc, ht_only=ht_only, ret_ctx=ret_ctx,
 )


def check_coll_la_core(data):
 op = data["operating_point"]
 pool = data["pool"]
 coll = data["collection"]
 la = data["la"]
 led = data["ledger"]
 lf, lg = op["lf"], op["lg"]
 p_g, p_f, p_lo = coll["p_g"], coll["p_f"], coll["p_lo"]
 p_g_shells = list(coll["p_g_shells"])
 assert len(p_g_shells) == 2 and abs(p_g - max(p_g_shells)) < 1e-9
 det = coll["det"]
 log2_Up = math.log2((2.0 ** log2_li(lf) + 2.0 ** log2_li(lg)) / 12.0)
 assert abs(log2_Up - op["log2_Up"]) < 1e-9
 p_min = int(coll.get("p_min", 33))
 q_g = float(det["b1_q_g"])
 q_f = float(det["b1_q_f"])
 tg_u8 = tuple(det["tg_u8"])
 tf_u8 = int(det["tf_u8"])
 _tgm = det.get("tg_mode", "static")
 _slg = float(det.get("sl_g_tgt", 376.0))
 _tlo = int(det.get("tg_u8_min", 60))
 _thi = int(det.get("tg_u8_max", 127))
 _rmg = max(round(b / q_g) / b for b in [33.0 + i / 16.0 for i in range(625)])
 sbg, sbf = float(det["sbg"]), float(det["sbf"])
 ret_f = QuantizedRetention(round(lf, 2), p_f, q_step=q_f, prime_min=2 ** p_min, explicit_bits=14, step=0.25, max_bits=660, tilt=0.05, odd=1)
 cst_f = QuantizedRetention(660, p_f, q_step=q_f, prime_min=2 ** p_min, explicit_bits=14, step=0.25, max_bits=661, tilt=0.005, odd=1)
 ret_g_sh = [QuantizedRetention(round(lg, 2), pgk, q_step=q_g, prime_min=2 ** p_min, explicit_bits=14, step=0.25, max_bits=740, tilt=0.05, odd=1) for pgk in p_g_shells]
 cst_g_sh = [QuantizedRetention(740, pgk, q_step=q_g, prime_min=2 ** p_min, explicit_bits=14, step=0.25, max_bits=741, tilt=0.005, odd=1) for pgk in p_g_shells]
 kmf = int(sbf // p_f)
 kmg_sh = [int(sbg // pgk) for pgk in p_g_shells]
 _kh_g = float(det.get("kh_g", 0.08))
 _kh_f = float(det.get("kh_f", 0.08))
 # Exact deterministic g-batch / f-batch (2D-56b NTT3 reusing existing ntt61_r4_u2, ZERO new ASM!)
 fbatch_tiers = [float(x) for x in det.get("fbatch_tiers", [lf])]
 _use_2d56 = (det.get("fbatch_ntt_mode") == "2d_56b_ntt3")
 _gb_fn = gbatch_ntt3_2d if _use_2d56 else gbatch_ntt3
 _ntt_pt_fn = ntt3_2d_56b_F if _use_2d56 else ntt3_F
 c_gbatch = _gb_fn(pbits_split(lg, _kh_g), leaf=11)
 if len(fbatch_tiers) >= 2:
  assert abs(fbatch_tiers[-1] - lf) < 1e-9 and p_f < fbatch_tiers[0]
  assert all(fbatch_tiers[i] < fbatch_tiers[i + 1] for i in range(len(fbatch_tiers) - 1))
  c_ftiers = [_gb_fn(pbits_split(bt, _kh_f), leaf=9) for bt in fbatch_tiers]
  c_mr2 = 0.0
  m_ftiers = [build_fbatch_tier_m(cst_f, p_f, bt, sbf, lf) for bt in fbatch_tiers[:-1]]
  fbatch_id_err = max(mt[3] for mt in m_ftiers)
 else:
  c_fbatch = _gb_fn(pbits_split(lf, _kh_f), leaf=9)
  m_ftiers = []
  fbatch_id_err = 0.0

 # Empirical prime-splitting density verification up to p <= 16384
 _ps14 = primes_upto(16384)
 _tb14 = _gb14 = _fb14 = 0.0
 for _p14 in _ps14:
  _lb14 = math.log2(_p14)
  _tb14 += _lb14
  if _p14 in (2, 3, 5):
   _gb14 += _lb14; _fb14 += _lb14; continue
  if _p14 % 15 != 1:
   continue
  _e14 = (_p14 - 1) // 15
  _s14 = next(t for g in range(2, 80) for t in (pow(g, _e14, _p14),) if pow(t, 3, _p14) != 1 and pow(t, 5, _p14) != 1)
  _w14 = pow(_s14, 5, _p14); _w24 = (_w14 * _w14) % _p14
  if pow((U_BLS * _w24) % _p14, (_p14 - 1) // 3, _p14) == 1 or pow((U_BLS * _w14) % _p14, (_p14 - 1) // 3, _p14) == 1:
   _gb14 += _lb14
  for _ww in (_w14, _w24):
   _fp = [(-_ww) % _p14, 0, 0, 1, 0, 0, (_ww + 1) % _p14, 0, 0, 1]
   _xp = ppow_mod([0, 1], _p14, _fp, _p14)
   _xm = list(_xp) + [0] * max(0, 2 - len(_xp))
   _xm[1] = (_xm[1] - 1) % _p14
   if len(pgcd_mod(_fp, _xm, _p14)) > 1:
    _fb14 += _lb14; break
 assert _gb14 / _tb14 <= _kh_g and _fb14 / _tb14 <= _kh_f

 _sb_cache = {}
 def _stage_b_sp(nbits, kh):
  _k = (int(round(nbits)), kh)
  if _k not in _sb_cache:
   L = max(5, math.ceil(nbits / 64.0))
   _sb_cache[_k] = gbatch_ntt3(pbits_split(p_min, kh), leaf=L, c_leaf=10.0 * MONT5 * ((L / 5.0) ** 2) + (nbits / 24.0) * (64.0 * L + 544.0))
  return _sb_cache[_k]

 # Product-tree memory & build cost (includes P_g, all fbatch_tiers trees, candidate tree, and exact root NTT workspace)
 _n_pt_own = 2048
 _sp_lg = math.ceil(pbits_split(lg, _kh_g) / 64.0)
 _L_pt_g = math.ceil(math.log2(_sp_lg / 11.0))
 _cand_tree_g = (_L_pt_g + 1.0) * (11.0 * (1 << _L_pt_g)) * 8.0
 _ntt_root_ws = (24.0 * ntt3_2d_56b_mtot(2 * _sp_lg)) if _use_2d56 else (4.0 * _sp_lg * 8.0)
 _pt_mem_tot = (_L_pt_g + 1.0) * _sp_lg * 8.0 + _cand_tree_g + _ntt_root_ws
 _pt_cost = math.ceil(math.log2(_sp_lg)) * _ntt_pt_fn(_sp_lg) + 100.0 * _sp_lg
 _sub_l11_max = max(
  (_L_pt_g - 11.0 + 1.0) * (_sp_lg / (1 << 11)) * 8.0,
  (_L_pt_g - 11.0 + 1.0) * (11.0 * (1 << (_L_pt_g - 11))) * 8.0,
 )
 for _bf_t in fbatch_tiers:
  _sp_t = math.ceil(pbits_split(float(_bf_t), _kh_f) / 64.0)
  _L_t = math.ceil(math.log2(_sp_t / 9.0))
  _pt_mem_tot += (_L_t + 1.0) * _sp_t * 8.0
  _pt_cost += math.ceil(math.log2(_sp_t)) * _ntt_pt_fn(_sp_t) + 100.0 * _sp_t
  _sub_l11_max = max(_sub_l11_max, (_L_t - 11.0 + 1.0) * (_sp_t / (1 << 11)) * 8.0)
 _pt_mem_own = _pt_mem_tot / _n_pt_own
 assert _pt_mem_own < 2.0 ** 64 and _sub_l11_max < 2.0 ** 62
 assert math.log2(_pt_cost / UNIT) <= dict(led["rows"])["pg_pf_product_trees"]

 alpha_f, alpha_g_ge3, kap2_g = op["alpha_f"], op["alpha_g_ge3"], op["kap2_g"]
 rf_c, sfa_c, sfb_c, qf_c = {}, {}, {}, {}
 sft_c = [{} for _ in m_ftiers]
 _rf_sm, _csb_g_c, _csb_f_c = {}, {}, {}
 v_shares = list(pool["v_shares"])
 hist_shells = pool["hist_shells"]
 assert len(v_shares) == 2 and len(hist_shells) == 2 and abs(sum(v_shares) - 1.0) < 1e-12
 hs = float(pool.get("hist_step", 0.25))
 hh = 0.5 * hs
 ysum_mid = ysum_rethi = ysum_strict_hi = c_cof_live = sga_mean = 0.0
 _nf_cond_num = _ng_cond_num = _w_cond_den = 0.0
 var_sh_tot = 0.0
 max_bg_seen = 0
 f_prim_tot = 0.0
 _f1_f = {}
 _f1_g_sh = [{}, {}]
 _rg_c_sh = [{}, {}]
 _tgf = lambda x: max(_tlo, min(_thi, max(int(round((x - _slg) / q_g)), 128 + math.floor(x * _rmg) - 255)))
 _tg_ok_c = set()
 _has_tgu = (pool.get("hist_enc") == "cls_tgu_v1")
 assert (_tgm != "per_col_bgmax_adaptive") or _has_tgu, "Adaptive T_g_u8 requires draw-level recorded tgu!"

 for k in range(2):
  vk = v_shares[k]
  rg_c = _rg_c_sh[k]
  sga_c, sgb_c, _rg_sm = {}, {}, {}
  ret_g = ret_g_sh[k]
  cst_g = cst_g_sh[k]
  kmg = kmg_sh[k]
  ym_k = yr_k = ysh_k = ysh2_k = cc_k = sga_k = wk_k = fp_k = 0.0
  _csb_g_get = _csb_g_c.get; _csb_f_get = _csb_f_c.get
  _rf_get = rf_c.get; _rg_get = rg_c.get; _rg_sm_get = _rg_sm.get
  _f1_f_get = _f1_f.get; _f1_gk = _f1_g_sh[k]; _f1_gk_get = _f1_gk.get
  _tg_ok_in = _tg_ok_c.__contains__; _tg_ok_add = _tg_ok_c.add
  _is_adp = (_tgm == "per_col_bgmax_adaptive")
  _has_mft = bool(m_ftiers); _c_g_extra = (2.0 ** 24.5) + (2.0 ** 24.0) + c_gbatch
  for _row in hist_shells[k]:
   bf, bgo, bgm, wc, wk = _row[:5]
   if bgm > max_bg_seen:
    max_bg_seen = bgm
   if _is_adp:
    _tgu = _row[5]
    _tk = (bgo, _tgu)
    if not _tg_ok_in(_tk):
     assert _tlo <= _tgu <= _thi and _tgu >= _tgf(bgo) and (128 - _tgu) + math.floor((bgo + hs) * _rmg) <= 246 <= 255
     _tg_ok_add(_tk)
   else:
    _cls = 0 if abs((bgm - bgo) - hs) < 1e-6 else (1 if abs((bgm - bgo) - (12.0 + hs)) < 1e-6 else 2)
    _tgu = tg_u8[_cls]
   _rf_v = _rf_get(bf)
   if _rf_v is None:
    _rf_v = (ret_f.at_u8(bf + hh, tf_u8, sbf), ret_f.at_u8(bf + hs, tf_u8, sbf))
    rf_c[bf] = _rf_v
    sfa_c[bf] = min(1.0, 1.25 * max(cst_f.survivor_upper_u8(float(bf), tf_u8), cst_f.survivor_upper_u8(float(bf + hs), tf_u8)))
    sfb_c[bf] = min(1.0, 1.25 * max(cst_f.survivor_upper_u8(float(bf), tf_u8, sbf, lf, kmf), cst_f.survivor_upper_u8(float(bf + hs), tf_u8, sbf, lf, kmf)))
    _pv_t = sfa_c[bf]
    for _ti, _mt in enumerate(m_ftiers):
     _vt = min(_pv_t, 1.25 * max(p_fbatch_tier_u8(cst_f, _mt, float(bf), tf_u8), p_fbatch_tier_u8(cst_f, _mt, float(bf + hs), tf_u8)))
     sft_c[_ti][bf] = _vt
     _pv_t = _vt
    qf_c[bf] = min(1.0, 2.0 * rho_exact((bf + alpha_f) / lf))
    ef_c = bf + hh + alpha_f; ef_h = bf + hs + alpha_f
    _rf_sm[bf] = (rho_db(ef_c / lf, ef_c), rho_db(ef_h / lf, ef_h))
   _gk = (bgo, _tgu)
   _rg_v = _rg_get(_gk)
   if _rg_v is None:
    _rg_v = (ret_g.at_u8(bgo + hh, _tgu, sbg), ret_g.at_u8(bgo + hs, _tgu, sbg))
    rg_c[_gk] = _rg_v
    sga_c[_gk] = min(1.0, max(cst_g.survivor_upper_u8(float(bgo), _tgu), cst_g.survivor_upper_u8(float(bgo + hs), _tgu)))
    sgb_c[_gk] = min(1.0, max(cst_g.survivor_upper_u8(float(bgo), _tgu, sbg, lg, kmg), cst_g.survivor_upper_u8(float(bgo + hs), _tgu, sbg, lg, kmg)))
   _rg_sm_v = _rg_sm_get(bgo)
   if _rg_sm_v is None:
    eg_c = bgo + hh + alpha_g_ge3 + 1.0 + kap2_g; eg_h = bgo + hs + alpha_g_ge3 + 1.0 + kap2_g
    _rg_sm_v = (rho_db(eg_c / lg, eg_c), rho_db(eg_h / lg, eg_h))
    _rg_sm[bgo] = _rg_sm_v
   if wc > 0.0:
    fp_k += wc
    _rf0, _rf1 = _rf_sm[bf]; _rg0, _rg1 = _rg_sm_v
    p_sm = _rf0 * _rg0; p_sm_hi = _rf1 * _rg1
    ym_k += wc * p_sm * _rf_v[0] * _rg_v[0]
    yr_k += wc * p_sm * _rf_v[1] * _rg_v[1]
    _ysh = p_sm_hi * _rf_v[1] * _rg_v[1]
    ysh_k += wc * _ysh; ysh2_k += wc * _ysh * _ysh
    if _rf_v[1] > 0.0 and _rg_v[1] > 0.0:
     _wsh = vk * wc * _ysh
     _f1_f[bf] = _f1_f_get(bf, 0.0) + _wsh
     _f1_gk[_gk] = _f1_gk_get(_gk, 0.0) + _wsh
     _nf_cond_num += _wsh * (bf + hh); _ng_cond_num += _wsh * (bgo + hh); _w_cond_den += _wsh
   if wk > 0.0:
    sga, sgb = sga_c[_gk], sgb_c[_gk]
    sfa, sfb, qf = sfa_c[bf], sfb_c[bf], qf_c[bf]
    csb_g = _csb_g_get(bgm)
    if csb_g is None:
     csb_g = _stage_b_sp(float(bgm), _kh_g); _csb_g_c[bgm] = csb_g
    csb_f = _csb_f_get(bf)
    if csb_f is None:
     csb_f = _stage_b_sp(float(bf + hs), _kh_f); _csb_f_c[bf] = csb_f
    if _has_mft:
     c_fb = (sgb * sfa) * c_ftiers[0]
     for _ti in range(len(m_ftiers)):
      c_fb += (sgb * sft_c[_ti][bf]) * c_ftiers[_ti + 1]
     c_fb += (sgb * (sfa if sfa < qf else qf)) * _c_g_extra
    else:
     q_f_cond = qf / sfb if sfb > 1e-12 else 1.0
     if q_f_cond > 1.0:
      q_f_cond = 1.0
     c_fb = (sgb * sfb) * (c_fbatch + q_f_cond * _c_g_extra)
    cc_k += wk * ((sga * sfa) * (6500.0 + csb_g) + (sgb * sfa) * (20000.0 + csb_f) + c_fb)
    sga_k += wk * sga
    wk_k += wk
  assert abs(wk_k - 1.0) < 1e-9
  ysum_mid += vk * ym_k; ysum_rethi += vk * yr_k; ysum_strict_hi += vk * ysh_k
  c_cof_live += vk * cc_k; sga_mean += vk * sga_k; f_prim_tot += vk * fp_k
  var_sh_tot += (vk * vk) * max(0.0, ysh2_k - ysh_k * ysh_k) / float(pool["n_draws_shell"])

 nf_cond_mean = _nf_cond_num / _w_cond_den
 ng_cond_mean = _ng_cond_num / _w_cond_den

 # Conditioned ret_mult scan
 _f_it = [(_bf, _wr / rho_db((_bf + hs + alpha_f) / lf, _bf + hs + alpha_f), _bf + hs + alpha_f, 1.0 / rf_c[_bf][1]) for _bf, _wr in _f1_f.items()]
 _g_it = [(_k, _bgo, _tgu, _wr / rho_db((_bgo + hs + alpha_g_ge3 + 1.0 + kap2_g) / lg, _bgo + hs + alpha_g_ge3 + 1.0 + kap2_g), _bgo + hs + alpha_g_ge3 + 1.0 + kap2_g, 1.0 / _rg_c_sh[_k][(_bgo, _tgu)][1]) for _k in range(2) for (_bgo, _tgu), _wr in _f1_g_sh[_k].items()]

 def _min_ret_mult(side_s, b_lo_100):
  _rd = rho_db; _rdc = {}; _rac = {}; _sf = (side_s == "f")
  if _sf:
   _ra = ret_f.at_u8; _il = 1.0 / lf; _it = [(int(round(_bf * 100)), _bf, _w0, _ef, _irf) for _bf, _w0, _ef, _irf in _f_it]
   def _ev(_bi):
    _b = _bi * 0.01; _num = _den = 0.0; _ok_s = (sbf > _b)
    for _bf100, _bf, _w0, _ef, _irf in _it:
     _ek = _bf100 - _bi; _rv = _rdc.get(_ek)
     if _rv is None:
      _eb = _ef - _b; _rv = _rd(_eb * _il, _eb if _eb > 1.0 else 1.0); _rdc[_ek] = _rv
     _ww = _w0 * _rv; _den += _ww
     if _ok_s and _bf + hs > _b:
      _rk = (_ek, _bi); _rav = _rac.get(_rk)
      if _rav is None:
       _rav = _ra(_bf + hs - _b, tf_u8, sbf - _b); _rac[_rk] = _rav
      _num += _ww * (_rav * _irf)
    return _num / _den
  else:
   _il = 1.0 / lg; _it = [(int(round(_bgo * 100)), _k, _bgo, _tgu, _w0, _eg, _irg) for _k, _bgo, _tgu, _w0, _eg, _irg in _g_it]
   def _ev(_bi):
    _b = _bi * 0.01; _num = _den = 0.0; _ok_s = (sbg > _b)
    for _bg100, _k, _bgo, _tgu, _w0, _eg, _irg in _it:
     _ek = _bg100 - _bi; _rv = _rdc.get(_ek)
     if _rv is None:
      _eb = _eg - _b; _rv = _rd(_eb * _il, _eb if _eb > 1.0 else 1.0); _rdc[_ek] = _rv
     _ww = _w0 * _rv; _den += _ww
     if _ok_s and _bgo + hs > _b:
      _rk = (_k, _ek, _tgu, _bi); _rav = _rac.get(_rk)
      if _rav is None:
       _rav = ret_g_sh[_k].at_u8(_bgo + hs - _b, _tgu, sbg - _b); _rac[_rk] = _rav
      _num += _ww * (_rav * _irg)
    return _num / _den
  _cgrid = list(range(7000, b_lo_100 - 1, -100)) + [b_lo_100]
  _rm, _best = min((_ev(_bi), _bi) for _bi in _cgrid)
  for _bi in range(max(b_lo_100, _best - 75), min(7000, _best + 75) + 1, 25):
   if _bi not in _cgrid:
    _r = _ev(_bi)
    if _r < _rm:
     _rm = _r
  return _rm

 ret_mult_f_min = _min_ret_mult("f", int(round(p_f * 100)) + 2)
 ret_mult_g_min = _min_ret_mult("g", int(round(min(p_g_shells) * 100)) + 2)
 _rm_g_br = 1e9; _qg_b = float(det["b1_q_g"])
 for _bi in (int(round(min(p_g_shells) * 100)) + 2, 6300, 6340, int(round(max(p_g_shells) * 100))):
  _b = _bi * 0.01; _num = _den = 0.0; _wt = int(round(_b / _qg_b))
  for _k, _bgo, _tgu, _w0, _eg, _irg in _g_it:
   _eb = _eg - _b; _rv = rho_db(_eb / lg, _eb if _eb > 1.0 else 1.0); _ww = _w0 * _rv; _den += _ww
   if sbg > _b and _bgo + hs > _b:
    _rav = ret_g_sh[_k].at_u8(_bgo + hs - _b, max(0, _tgu - _wt), sbg) if _b <= p_g_shells[_k] else ret_g_sh[_k].at_u8(_bgo + hs - _b, _tgu, sbg - _b)
    _num += _ww * (_rav * _irg)
  _rm_g_br = min(_rm_g_br, _num / _den)
 assert ret_mult_f_min >= 0.90 and ret_mult_g_min >= 0.82 and _rm_g_br >= 1.0, "Conditioned ret_mult below (0.90, 0.82) or bridge < 1.0!"
 assert abs(c_cof_live - det["c_cof_2d_upper"]) < 1e-4, f"c_cof_live={c_cof_live:.5f} != {det['c_cof_2d_upper']}"
 assert max_bg_seen >= 684, "2D histogram must cover high-norm tail!"
 assert abs(f_prim_tot - pool["f_prim"]) < 1e-6 and 0.94 < f_prim_tot < 0.97

 log2_rels_h = (pool["log2_base"] + pool["dV"]) + math.log2(ysum_mid)
 log2_rels_rethi = (pool["log2_base"] + pool["dV"]) + math.log2(ysum_rethi)
 log2_rels_strict_hi = (pool["log2_base"] + pool["dV"]) + math.log2(ysum_strict_hi)
 ecm_coll_succ = 0.999999
 ecm_rel_keep = 1.0 - 29.0 * (1.0 - ecm_coll_succ)
 det_cred_pen = -math.log2(float(coll.get("det_credited", 0.99)) * ecm_rel_keep)
 margin_net = (log2_rels_rethi - det_cred_pen) - (log2_Up + math.log2(op["exc"]))
 margin_strict_hi = (log2_rels_strict_hi - det_cred_pen) - (log2_Up + math.log2(op["exc"]))
 assert margin_net >= 0.025 and margin_strict_hi >= 0.012 and margin_strict_hi >= 4.5 * pool["se_y_bits"]
 _se_sh_b = math.sqrt(var_sh_tot) / (ysum_strict_hi * LN2)
 log2_rels_2se_live = log2_rels_strict_hi - det_cred_pen - 2.0 * max(_se_sh_b, pool["se_y_bits"])
 assert log2_rels_2se_live >= pool["log2_rels"]
 assert coll["h2_allow"] >= 1.0 + 8.0 * max(pool["se_c_rel"], pool["se_w_rel"])

 # Stage-1 shell geometry & pricing (K=2 L_{3/4}(b), M=15 slabs, split p_rot=56.5 4-rot W_B=1)
 V_lin = 2.0 ** pool["V_kept"]
 p2_mult = pool["p2_mult"]
 N_a_k = 4.0 * (2.0 ** pool["log2_Va"])
 R1_shells = list(pool["R1_shells"])
 Lb_shells = list(pool["Lb_shells"])
 L_max_shells = list(pool["L_max_shells"])
 N7_shells = [2.0 ** x for x in pool["log2_N7_shells"]]
 lb_r1 = float(pool["Lb_over_R1_charged"])
 lmax_c = float(pool.get("L_max_coeff", 0.5))
 for k in range(2):
  assert abs(L_max_shells[k] - lmax_c * R1_shells[k]) < 1e-4 and abs(Lb_shells[k] - lb_r1 * R1_shells[k]) < 1e-4 and L_max_shells[k] < coll["w_tile"]
  Vk = V_lin * v_shares[k]
  assert abs(N_a_k * N7_shells[k] * Lb_shells[k] - Vk) / Vk < 1e-6, f"Shell {k} volume mismatch!"
 w_tile, w_sb, p_m4, p_m4_s = coll["w_tile"], coll["w_sb"], coll["p_m4"], coll["p_m4_s"]
 bp = int(coll["bp"])
 n_slabs = int(coll.get("n_slabs", 15))
 wb2 = int(coll.get("w_b2_mid", 8))
 assert bp == 1 and w_sb == 32 and n_slabs in (13, 15) and w_tile == n_slabs * 32 and max(L_max_shells) + 31.0 < w_tile
 p_lo_g = float(coll.get("p_lo_g", p_lo))
 p_rot = float(coll.get("p_rot", p_lo))
 c_hit_rot = float(coll.get("c_hit_hi_rot", coll["c_hit_hi_u4"]))
 p_w1 = float(coll["p_w1"])
 p_dir = float(coll["p_dir"])
 c_list_w1 = float(coll["c_list_a_w1"])
 c_list_b2 = float(coll["c_list_a_b2"])
 c_mset = float(coll["c_memset_dir"])
 c_lbkt = float(coll["c_lbkt"])
 bld_ovh_sh = [float(x) for x in pool["bld_ovh_shells"]]
 assert len(bld_ovh_sh) == 2

 def price_shell_side(p_hi, side, k):
  N7_k = N7_shells[k]; Lb_k = Lb_shells[k]; Vk = V_lin * v_shares[k]; ovh_k = bld_ovh_sh[k]
  pl = min(p_lo_g if side == "g" else p_lo, p_hi)
  pm = min(p_m4, pl); pms = min(p_m4_s, pm)
  n1al = 1.01 * (2.0 ** log2_li(pms) - 2.0 ** log2_li(float(p_min)))
  n1as = 1.01 * (2.0 ** log2_li(pm) - 2.0 ** log2_li(pms)) if pm > pms else 0.0
  n1b = 1.01 * (2.0 ** log2_li(pl) - 2.0 ** log2_li(pm)) if pl > pm else 0.0
  b1_lk = n1al * float(coll["b1a_long_lk"]) + n1as * float(coll["b1a_short_lk"]) + n1b * float(coll["b1bl_lk"])
  c1al = float(coll["c_hit_lo_m4_" + side]); c1as = float(coll["c_hit_lo_m4_s_" + side]); c1bl = float(coll["c_hit_1bl_" + side])
  ib1 = (math.log(pms / float(p_min)) * c1al + (math.log(pm / pms) * c1as if pm > pms else 0.0) + (math.log(pl / pm) * c1bl if pl > pm else 0.0)) * p2_mult
  ib1 += (b1_lk * N7_k * Lb_k + bp * (n1al + n1as + n1b) * N7_k * 24.0) / Vk
  sp1a = 4.0 * 1.01 * (2.0 ** (2.0 * pm)) / (2.0 * pm * LN2 - 1.0) / 4.0
  sp1b = (bp * 1.01 * (2.0 ** (2.0 * pl) - 2.0 ** (2.0 * pm)) / (2.0 * pl * LN2 - 1.0) / 4.0) if pl > pm else 0.0
  bld = (n1al + n1as + n1b) * (c_list_w1 + ovh_k) * N_a_k + c_lbkt * (sp1a + sp1b) + (4.0 * (n1al + n1as) + bp * n1b) / 4.0
  mem = (n1al + n1as + n1b) * 8.0 * N_a_k + 8.0 * (sp1a + sp1b)
  ipr = ihi = 0.0
  if p_hi > pl:
   sb32 = Lb_k / 32.0 + 1.5
   p1h = min(p_w1, p_hi)
   if p1h > pl:
    _n_wb = len({int(round((pl + 0.05 * i) / float(det["b1_q_" + side]))) for i in range(int(round((p_hi - pl) * 20)) + 1)})
    assert _n_wb <= 4 and _n_wb * n_slabs * 109 <= 23 * 15 * 60
    # Split W_B=1 at p_rot: [pl, pr_hi) uses 4 pre-rotated tables (rev_b2_emit_id_u16 on 4/4 conjugates),
    # [pr_hi, p1h) uses 1 shared table (rev_b2_emit_id/c123).
    pr_hi = min(max(p_rot, pl), p1h)
    if pr_hi > pl:
     n2_r = 1.01 * (2.0 ** log2_li(pr_hi) - 2.0 ** log2_li(pl))
     ipr += bp * n2_r * N7_k * (float(coll["c_tile_hdr_w1"]) + coll["c_b2_row_k0"] + sb32 * float(coll["c_sb_hdr_u4"])) / Vk
     ihi += math.log(pr_hi / pl) * c_hit_rot * p2_mult
     sp2_r = 4.0 * bp * 1.01 * (2.0 ** (2.0 * pr_hi) - 2.0 ** (2.0 * pl)) / (2.0 * pr_hi * LN2 - 1.0) / 4.0
     bld += n2_r * (c_list_w1 + ovh_k) * N_a_k + c_lbkt * sp2_r + bp * n2_r
     mem += n2_r * 8.0 * N_a_k + 8.0 * sp2_r
    if p1h > pr_hi:
     n2_s = 1.01 * (2.0 ** log2_li(p1h) - 2.0 ** log2_li(pr_hi))
     ipr += bp * n2_s * N7_k * (float(coll["c_tile_hdr_w1"]) + coll["c_b2_row_k0"] + sb32 * float(coll["c_sb_hdr_u4"])) / Vk
     ihi += math.log(p1h / pr_hi) * float(coll["c_hit_hi_u4"]) * p2_mult
     sp2_s = bp * 1.01 * (2.0 ** (2.0 * p1h) - 2.0 ** (2.0 * pr_hi)) / (2.0 * p1h * LN2 - 1.0) / 4.0
     bld += (n2_s / 4.0) * (c_list_w1 + ovh_k) * N_a_k + c_lbkt * sp2_s + bp * (n2_s / 4.0)
     mem += (n2_s / 4.0) * 8.0 * N_a_k + 8.0 * sp2_s
    ihi += bp * (23.0 * 15.0 * 60.0) / (N_a_k * Lb_k)
   p2h = min(p_dir, p_hi)
   if p2h > p_w1:
    n2_8 = 1.01 * (2.0 ** log2_li(p2h) - 2.0 ** log2_li(p_w1))
    cand_bd = 1.0 + wb2 / Lb_k
    ipr += bp * n2_8 * N7_k * (float(coll["c_tile_hdr_w8"]) + coll["c_b2_row_k0"] + sb32 * float(coll["c_sb_hdr_u2"])) / Vk
    ihi += math.log(p2h / p_w1) * float(coll["c_hit_hi_u2"]) * cand_bd * p2_mult
    sp2_8 = bp * 1.01 * (2.0 ** (2.0 * p2h) - 2.0 ** (2.0 * p_w1)) / (2.0 * p2h * LN2 - 1.0) / 4.0
    bld += (n2_8 / 4.0) * (c_list_b2 + ovh_k) * N_a_k + c_lbkt * (sp2_8 / wb2) + bp * (n2_8 / 4.0)
    mem += (n2_8 / 4.0) * 8.0 * N_a_k + 8.0 * (sp2_8 / wb2)
   if p_hi > p_dir:
    np_pan = 64; dx = (p_hi - p_dir) / np_pan
    for pi in range(np_pan):
     x0 = p_dir + pi * dx; x1 = x0 + dx; xm = 0.5 * (x0 + x1); pm_v = 2.0 ** xm
     n_p = 1.01 * (2.0 ** log2_li(x1) - 2.0 ** log2_li(x0))
     lam = w_sb * N_a_k / pm_v; p_ne = 1.0 - math.exp(-lam); ovf_m = max(0.0, lam - p_ne)
     bld += (n_p / 4.0) * ((pm_v / w_sb) * c_mset + N_a_k * (11.25 + ovh_k + 3.0 * (ovf_m / lam)) + 1.0)
     mem += (n_p / 4.0) * (8.0 * (pm_v / w_sb) + 16.0 * N_a_k * p_ne)
     c_sb_d = 3.5 + 13.0 * p_ne + 11.0 * ovf_m + 1.54 * lam
     ipr += bp * n_p * N7_k * (coll["c_tile_hdr"] + coll["c_b2_row_k0"] + 18.0 + sb32 * c_sb_d) / Vk
  return bld, mem, ib1 + ipr + ihi

 bg_tot = bf_tot = ipc_sieve = mem_tbl_peak = 0.0
 for k in range(2):
  bgk, mgk, ipgk = price_shell_side(p_g_shells[k], "g", k)
  bfk, mfk, ipfk = price_shell_side(p_f, "f", k)
  bg_tot += bgk; bf_tot += bfk; ipc_sieve += v_shares[k] * (ipgk + ipfk)
  mem_tbl_peak = max(mem_tbl_peak, mgk + mfk)
 ipc_tot = ipc_sieve + coll["c_scan_zero"] + c_cof_live
 coll_u = pool["V_kept"] + math.log2(ipc_tot * coll["h2_allow"] / UNIT)
 build_u = math.log2((bg_tot + bf_tot) / UNIT)
 rows_d = dict(led["rows"])
 assert coll_u <= rows_d["relation_collection"] + 1e-4, f"coll_u={coll_u:.5f} > {rows_d['relation_collection']}"
 assert build_u <= rows_d["replay_a_tables"] + 1e-4, f"build_u={build_u:.5f} > {rows_d['replay_a_tables']}"

 _sc_mem = 0.055 * float(1 << 60)
 col_own_mem = (1.0 + 32.0 + 4.0) * float(1 << 58)
 m1a_bld = 16.0 * N_a_k + 32.0 * (2.0 ** p_m4) + w_tile * (1 << 20) + _sc_mem
 m2_rot_bld = 16.0 * N_a_k + 32.0 * (2.0 ** p_rot) + w_tile * (1 << 20) + _sc_mem
 m2_w1_bld = 16.0 * N_a_k + 8.0 * (2.0 ** p_w1) + w_tile * (1 << 20) + _sc_mem
 m2_w8_bld = 16.0 * N_a_k + 8.0 * ((2.0 ** p_dir) / wb2) + w_tile * (1 << 20) + _sc_mem
 _lam_dir = w_sb * N_a_k / (2.0 ** p_g)
 m2_dir_bld = 8.0 * ((2.0 ** p_g) / w_sb) + 16.0 * N_a_k * (1.0 - math.exp(-_lam_dir)) + w_tile * (1 << 20) + _sc_mem
 assert 2.0 ** pool["log2_Va"] < (1 << 56) * 0.95 and col_own_mem < 2.0 ** 64, "Fused u8 slab owner (9.25 EiB) must fit < 2^64 B!"
 assert max(m1a_bld, m2_rot_bld, m2_w1_bld, m2_w8_bld, m2_dir_bld) < 2.0 ** 64, "All 1-owner build peaks must fit < 2^64 B!"
 assert mem_tbl_peak + n_slabs * col_own_mem < 2.0 ** 120, "Sequential-shell peak tables + slabs must fit < 2^120 B!"
 mem_tbl_tot = mem_tbl_peak

 n_rels_upper_log2 = pool["log2_rels"] + 0.55
 assert n_rels_upper_log2 <= 63.95
 n_rels_upper = 2.0 ** n_rels_upper_log2
 if _use_2d56:
  c_rel_band, band_u_live = band_split_2d_ntt3_cost(48.0, n_rels_upper)
 else:
  np_band = 1.25506 * (2.0 ** 48) / (48.0 * LN2)
  c_rel_band = np_band * (6 * (10 + 12) + 44 + 4) + (13 + 16) * 1024 + (1 << 26)
  band_u_live = math.log2(n_rels_upper * c_rel_band / UNIT)
 assert band_u_live <= rows_d["band_split"]

 ef_e, ef_v = ecm_mean_var(lf, 640, success=ecm_coll_succ)
 eg_e, eg_v = ecm_mean_var(lg, 768, success=ecm_coll_succ)
 assert ef_v == 0.0 and eg_v == 0.0
 c_ecm_rel = 13.0 * ef_e + 14.0 * mr_instr(640) + 13.0 * (1 << 26) + 16.0 * eg_e + 17.0 * mr_instr(768) + 16.0 * (1 << 26)
 ew_ecm_live = (c_ecm_rel / UNIT) * n_rels_upper
 assert (2.0 ** rows_d["ecm_split_expected"]) - ew_ecm_live > 0
 assert n_rels_upper >= (2.0 ** pool["log2_rels"]) * 1.46 and ecm_rel_keep < 1.0

 _h3 = verify_h3_ladder_and_fit(la, log2_Up, run_u12=False, op=op, nf_cond=nf_cond_mean, ng_cond=ng_cond_mean)
 k5_live = _h3["k5_live"]
 assert led["cap_abort"] == 0.0
 assert la["tab_ratio_charged"] == 1.0 and la["n_owners"] == 1 and la["te"] == 0
 N = 2.0 ** (log2_Up + math.log2(la["nl_cap"]))
 N_sw = 2.0 ** (log2_Up + math.log2(la["sw_cap"]))
 m, n, n_sm = la["m"], la["n"], la["n_sm"]
 la_mem_full = 2.0 * N_sw * float(la.get("bytes_per_slot", 128.0))
 horner_1own_mem = 2.0 * N_sw * 32.0
 _n_proj_own = int(la.get("n_proj_owners", 16))
 la_proj_mem = (124.0 * 32.0 * N + 128.0 * 32.0 * N) / _n_proj_own + 16.0 * N
 _L_bm_mem = (1.0 + 2.0 ** -7) * (N + n_sm) * (1.0 / m + 1.0 / n) + 6.0
 _n_seq_own = 2048
 la_seq_mem = (_L_bm_mem * m * n * 32.0) / _n_seq_own
 _n_pm_own = 1 << 22
 _M_pm_sub = float(1 << math.ceil(math.log2(_L_bm_mem)))
 la_pm_mem = (2.0 * _L_bm_mem * ((m + n) ** 2) * 32.0 + 3.0 * 4.0 * _M_pm_sub * ((m + n) ** 2) * 8.0) / _n_pm_own
 la_gen_vec_mem = _L_bm_mem * n * 32.0
 _nu_ub = float(la.get("nonunit_ub", 0.25))
 _n_filt_own = 8192
 filt_mut_mem = (n_rels_upper * (64.0 + 201.0 * (8.0 + _nu_ub) + 7.0 + 201.0 * 10.0 + 160.0)) / _n_filt_own
 _n_purge_own = 2048
 purge_mut_mem = (n_rels_upper * 1040.0) / _n_purge_own
 assert filt_mut_mem <= 8.0 * (2.0 ** 60) and purge_mut_mem <= 8.0 * (2.0 ** 60)
 assert all((1.0 - x["exp_unit"]) <= _nu_ub for x in la["ladder"][2:])
 n_mat_own = int(la.get("n_mat_owners", 16))
 assert 1.0 - (la["eu_26"] - 2.0 * la["se_eu_26"]) <= _nu_ub
 assert la["k5_charged"] >= (1.0 - _nu_ub) * 7.65 + _nu_ub * 22.0 and float(la.get("k5_horner_charged", 14.0)) >= (1.0 - _nu_ub) * 9.60 + _nu_ub * 24.0
 _k5_h_live = (1.0 - _nu_ub) * 9.60 + _nu_ub * 24.0
 _k5_h_chg = float(la.get("k5_horner_charged", 14.0))
 assert all(_k5_h_chg > float(x["k5"]) + 1.95 for x in la["ladder"])
 _w_max_row = 4096
 assert _w_max_row < 0xFFFF
 mat_shard_mem = (N_sw * (64.0 + 201.0 * (8.0 + _nu_ub) + 7.0)) / n_mat_own + (64.0 + _w_max_row * 9.0)
 _w_row_avg = 8.0 + 201.0 + (201.0 * _nu_ub) / 8.0 + 7.0 / 8.0
 _c_sm4 = float(la["c_stream_mat"]); _c_sm1 = float(la["c_stream_mat_1vec"])
 _ops_row_agg = 1.0625 * _w_row_avg + 24.0
 assert _c_sm4 >= _ops_row_agg / 4.0 and _c_sm1 >= _ops_row_agg
 assert max(la_mem_full, horner_1own_mem, la_proj_mem, mat_shard_mem, la_seq_mem, la_pm_mem, la_gen_vec_mem, filt_mut_mem, purge_mut_mem) < 2.0 ** 64
 _c_xproj = float(la.get("c_xproj", 70.0))
 assert _c_xproj >= 64.5
 c_prod = (la["sw_cap"] / la["nl_cap"]) * (201.0 * la["k5_charged"] * la["tab_ratio_charged"] + la["c_row_chi"] + _c_sm4)
 c_prod_horner = (la["sw_cap"] / la["nl_cap"]) * (201.0 * _k5_h_chg * la["tab_ratio_charged"] + la["c_row_chi"] + _c_sm1)
 bw_sig = 2.0 ** -7
 kap = ((1.0 + bw_sig) * (N + n_sm) + 6.0 / (1.0 / m + 1.0 / n)) / N
 per_n2 = (c_prod + _c_xproj) * (n / m + 1.0) * kap + (1.0 / m + (1.0 + 1.0 / 128.0) / n) * (c_prod_horner + n * (97.0 + 4.5) + MONT5)
 la_u = math.log2((N * N * per_n2 + N * (m + n) * MONT5) / UNIT)
 L_bm = (1.0 + bw_sig) * (N + n_sm) * (1.0 / m + 1.0 / n) + 6.0
 t_bm = m + n; depth_bm = math.ceil(math.log2(L_bm)); k_bm = math.ceil(math.log2(4.0 * L_bm))
 assert (1 << k_bm) <= (1 << 53) and t_bm * L_bm * (1 << 128) < NTT4_PROD
 w_base_bm = L_bm * m * (t_bm ** 2) * MONT5
 w_rec_bm = 2.0 * depth_bm * (8.0 * L_bm) * (1536.0 * (t_bm ** 3) + (t_bm ** 2) * (1536.0 * k_bm + 6220.0))
 gen_u = math.log2((w_base_bm + w_rec_bm) / UNIT)
 assert la_u <= rows_d["linear_algebra"] and gen_u <= rows_d["matrix_generator"]
 _RCTX[:] = [(ret_f, ret_g_sh, tf_u8, sbf, sbg, log2_rels_h)]
 return dict(
  c_cof_2d_upper_live=round(c_cof_live, 5), fbatch_4band_id_err=fbatch_id_err, sga_mean=sga_mean,
  ret_mult_f_min=round(ret_mult_f_min, 5), ret_mult_g_min=round(ret_mult_g_min, 5), ret_mult_g_bridge_min=round(_rm_g_br, 5),
  nf_cond_mean=round(nf_cond_mean, 3), ng_cond_mean=round(ng_cond_mean, 3), max_bg_bin_covered=max_bg_seen,
  log2_rels_2d_mid=round(log2_rels_h, 6), log2_rels_2d_rethi=round(log2_rels_rethi, 6),
  log2_rels_2d_strict_hi=round(log2_rels_strict_hi, 6), log2_rels_2se_live=round(log2_rels_2se_live, 6),
  margin_net_rethi_bits=round(margin_net, 5), margin_strict_hi_bits=round(margin_strict_hi, 5),
  L_max_shells=[round(x, 2) for x in L_max_shells], Lb_shells=[round(x, 2) for x in Lb_shells],
  ipc_tot=round(ipc_tot, 4), coll_u=round(coll_u, 6), build_u=round(build_u, 4),
  pt_2048own_frac=round(_pt_mem_own / (2.0 ** 64), 4), pt_l11_sub_max_log2_B=round(math.log2(_sub_l11_max), 3),
  slab_owner_eib=round(col_own_mem / (2.0 ** 60), 3),
  max_1own_bld_eib=round(max(m1a_bld, m2_rot_bld, m2_w1_bld, m2_w8_bld, m2_dir_bld) / (2.0 ** 60), 3),
  mem_tables_log2_B=round(math.log2(mem_tbl_tot), 3),
  band_u=round(band_u_live, 4), ecm_u=round(math.log2(ew_ecm_live), 4), ecm_var_live=0.0, ecm_rel_keep=round(ecm_rel_keep, 6),
  k5_live=round(k5_live, 4), k5_horner_live=round(_k5_h_live, 4), la_u=round(la_u, 6), gen_u=round(gen_u, 4),
  la_mem_full_log2=round(math.log2(la_mem_full), 3), la_proj_mem_log2=round(math.log2(la_proj_mem), 3), n_proj_owners=_n_proj_own,
  la_seq_mem_log2=round(math.log2(la_seq_mem), 3), n_seq_owners=_n_seq_own,
  la_pm_mem_log2=round(math.log2(la_pm_mem), 3), n_pm_owners=_n_pm_own,
  la_gen_vec_mem_log2=round(math.log2(la_gen_vec_mem), 3), n_filt_owners=_n_filt_own,
  filt_mut_mem_log2=round(math.log2(filt_mut_mem), 3), purge_mut_mem_log2=round(math.log2(purge_mut_mem), 3),
  mat_shard_mem_log2=round(math.log2(mat_shard_mem), 3), n_mat_owners=n_mat_own,
  c_stream_mat=float(la["c_stream_mat"]), N_la=N, N_sw=N_sw, f_prim_mc_2m=round(f_prim_tot, 6),
 )


def check_hx_filtered_live():
 res = run_hx_on_filtered_z12(log2U=12, seed=5, m=12, n=8, s_sm=4)
 assert (
  res["checks_ok"] and res["l_rec_strict_lower_tri"] and res["prod_eq_schur"]
  and res["pm_basis_lead_sm_nonzero"] and res["pm_basis_zero_sm_ctrl_zero_vec"]
  and res["pm_basis_horner_passes"] == 1 and res["rank_M"] == res["N"]
  and res["rank_A"] == res["D"] - 1 and res["kernel_eq_zero"] and res["recovered_true_kernel"]
  and res["post_peel_cols_lost"] == 0 and res["post_merge_cols_lost"] == 0
  and res["repurge_backsub_errors"] == 0 and res["n_repurge_backsub"] == res["repurge_stack_len"]
  and res["purge_backsub_errors"] == 0 and res["n_purge_backsub"] > 0
  and res["n_reconstruct_unit"] == 0 and res["all_2core_rels_zero"]
  and res["recovered_cols_match_true"] and res["cols_unrecovered_tri"] == 0
  and res["backsub_2core_exact"]
 )
 return res


def main():
 raw = sys.stdin.read().strip()
 req = json.loads(raw) if raw else {}
 allowed_keys = {"part", "descent_mesh"}
 if not set(req).issubset(allowed_keys):
  raise ValueError("audit_cm.py accepts {}, {'part':'core'|'diagnostics'|'all'}, or {'descent_mesh':384}")
 part = req.get("part", "core")
 if part not in ("core", "diagnostics", "all"):
  raise ValueError("part must be 'core', 'diagnostics', or 'all'")
 import unitdomain as _UD
 import verify_q15_descent as _VD
 data = _UD.load_pkg_tables(ROOT)
 _dm = req.get("descent_mesh", 384)
 if _dm not in (384,):
  raise ValueError("descent_mesh must be 384")
 out = {}
 ck = {}
 if part in ("diagnostics", "all"):
  toy = json.loads(lzma.decompress(base64.b85decode((ROOT / "toy_e2e.b85").read_bytes())).decode())
  out["murphy_alpha_exact_live"] = _UD.verify_murphy_alpha_live(data["operating_point"], U_BLS, primes_upto, pgcd_mod, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, rho_db); ck["murphy_alpha_exact_live"] = True
  out["real_fb_joint_coll_live"] = _UD.verify_real_fb_joint_coll(data["pool"]["voronoi_units_10"], build_G1w, build_WS_1w, mat_mul, M_1W, PHI15, y_pow_z, add_phi15_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, pgcd_mod, primes_upto, rho_db, Retention); ck["real_fb_joint_coll_live"] = True
  out["toy_e2e_dlp"] = check_toy_e2e_dlp(toy); ck["toy_e2e_dlp"] = True
  out["filtsim_and_hx_live"] = check_hx_filtered_live(); ck["filtsim_and_hx_live"] = True
  out["missing_recon_6b_live"] = missing_recon_z12(14, 0, lf=71.20, lg=71.56, nf=545.0, ng=665.0, excess=1.02, p_f=60.4, p_g=62.4, ret_f=0.916, ret_g=0.82, sbf=362.0, sbg=315.0); ck["missing_recon_6b_live"] = True
  _h3r = verify_h3_ladder_and_fit(data["la"], data["operating_point"]["log2_Up"], run_u12=True, op=data["operating_point"], nf_cond=545.242, ng_cond=665.034)
  _u12p = _h3r.get("live_u12_seed0")
  out["h3_ladder_and_genuine_q15_hx"] = dict(h3=_h3r, genuine_q15_hx=verify_genuine_q15_filter_hx_live(toy, primes_upto, PHI15, y_pow_z, sub_phi15_z, mul_phi15_z, pow_phi15_z, scale_phi15_z, norm_q15_z, pgcd_mod, pmul_mod, prem_mod, u12_proxy=_u12p))
  ck["h3_ladder_and_genuine_q15_hx"] = True
  out["unplanted_two_sided_bw"] = unplanted_two_sided_bw_check()
  assert out["unplanted_two_sided_bw"]["pm_basis_lead_sm_nonzero"] and out["unplanted_two_sided_bw"]["pm_basis_zero_sm_ctrl_zero_vec"] and out["unplanted_two_sided_bw"]["pm_basis_horner_passes"] == 1
  ck["unplanted_two_sided_bw"] = True
 if part == "diagnostics":
  out["lq_g120_ht_live"] = verify_pool_2adic_live(data["pool"], data["operating_point"], run_ht_mc=True, ht_only=True); ck["lq_g120_ht_live"] = True
  dres = _VD.run_all_checks(root_dir=ROOT, seed=1, mesh=_dm, run_diagnostics=True, core_checks=False)
  assert dres["all_ok"] is True and dres["mini_adaptive_campaign"]["checks_ok"] is True and dres["mini_adaptive_campaign"]["live_witnesses_replayed"] >= 3
  assert dres["mini_adaptive_campaign"]["full16d_k1_cohort"]["combined_zero_count"] == 0 and dres["mini_adaptive_campaign"]["prod_raw_lq_survey"]["all38_below_1x"] == 0 and dres["mini_adaptive_campaign"]["prod_raw_lq_survey"]["reached_below_1x"] == 0
  assert dres["live_full16d_k1_gen"]["total_ideals"] >= 32 and dres["live_full16d_k1_gen"]["total_empty"] == 0 and dres["live_full16d_k1_gen"]["reached_children"] >= 12
  assert dres["prod_q70_cert_diag"]["radius_rescaling"] is False and dres["prod_q70_cert_diag"]["bound_cell_screened"] is True and dres["prod_q70_cert_diag"]["exact_norms_floor_checked"] is True and dres["prod_q70_cert_diag"]["cells_per_lattice"] == 32 and dres["prod_q70_cert_diag"]["n_checked"] == 4 and dres["prod_q70_cert_diag"]["n_below_1x"] == 0 and dres["prod_q70_cert_diag"]["min_ratio"] >= 2.0
  assert dres["joint_smooth_700k"]["ok"] is True and dres["highbit_prod_diag"]["ok"] is True
  _eh5 = dres["exact_deg72_h5"]
  assert _eh5["usable_implies_surv"] is True and _eh5["disjoint_c12_orbits"] == 100 and _eh5["live_replays"] == 6 and _eh5["pooled_obs_over_mu_pred_2core"] >= 1.0
  out["q15_modeA_descent_diagnostics"] = dres; ck["q15_modeA_descent_diagnostics"] = True
  res = dict(checks_ok=all(ck.values()), part=part, checks=ck, summary=out)
  print(json.dumps(res, separators=(",", ":")))
  if not res["checks_ok"]:
   raise SystemExit(1)
  return

 out["algebra_bls12"] = check_algebra_bls12(); ck["algebra_bls12"] = True
 out["sigmaB_exactness"] = check_sigmaB_exactness(); ck["sigmaB_exactness"] = True
 out["coll_la_core"] = check_coll_la_core(data); ck["coll_la_core"] = True
 out["anchor_norms"] = check_anchor_norms(data["pool"]["anchors"])
 out["pool_2adic_live"] = verify_pool_2adic_live(data["pool"], data["operating_point"], run_ht_mc=(part == "all"), ret_ctx=_RCTX[0])
 ck["anchor_norms"] = True
 out["voronoi_exact"] = check_voronoi_exact(data); ck["voronoi_exact"] = True
 out["srt_toy"] = srt_toy_check(); ck["srt_toy"] = True
 _dp_kern = data.get("descent_kernels", {"lll_node": 2.0 ** 40.0, "lll_boot": 2.0 ** 47.0, "pair_setup": 2.0 ** 35.2, "bound_cell": 4400.0})
 out["kernels_symbolic"] = verify_kernels_symbolic(
  (ROOT / "KERNELS.md").read_text() + "\n" + (ROOT / "KERNELS_LA.md").read_text(),
  data["collection"], data["la"], R_BLS,
  float(_dp_kern["lll_node"]), float(_dp_kern["lll_boot"]), _dp_kern,
  sga_mean=out["coll_la_core"]["sga_mean"], lb_live=min(data["pool"]["Lb_shells"]),
 )
 ck["kernels_symbolic"] = True
 dres = _VD.run_all_checks(root_dir=ROOT, seed=1, mesh=_dm, run_diagnostics=(part == "all"))
 out["q15_modeA_descent"] = dres; ck["q15_modeA_descent"] = True

 _cl = out["coll_la_core"]; _op = data["operating_point"]; _led = data["ledger"]
 assert dres["all_ok"] is True and dres["mode"] == "pure_k1_boole_wald_iid_boot"
 assert dres["nb_r"] == 1.0 and dres["nb_alpha"] == 0.5 and dres["p0_floor"] == 0.001
 assert dres["phi_tree"] <= 0.48 and dres["boot_phi"] == 0.0 and dres["p_descent_fail"] == 0.0
 if part == "all":
  assert dres["mini_adaptive_campaign"]["checks_ok"] is True and dres["mini_adaptive_campaign"]["live_witnesses_replayed"] >= 3
  assert dres["mini_adaptive_campaign"]["full16d_k1_cohort"]["combined_zero_count"] == 0 and dres["mini_adaptive_campaign"]["prod_raw_lq_survey"]["all38_below_1x"] == 0 and dres["mini_adaptive_campaign"]["prod_raw_lq_survey"]["reached_below_1x"] == 0
  assert dres["live_full16d_k1_gen"]["total_ideals"] >= 32 and dres["live_full16d_k1_gen"]["total_empty"] == 0 and dres["live_full16d_k1_gen"]["reached_children"] >= 12
  assert dres["prod_q70_cert_diag"]["radius_rescaling"] is False and dres["prod_q70_cert_diag"]["bound_cell_screened"] is True and dres["prod_q70_cert_diag"]["exact_norms_floor_checked"] is True and dres["prod_q70_cert_diag"]["cells_per_lattice"] == 32 and dres["prod_q70_cert_diag"]["n_checked"] == 4 and dres["prod_q70_cert_diag"]["n_below_1x"] == 0 and dres["prod_q70_cert_diag"]["min_ratio"] >= 2.0
  assert dres["joint_smooth_700k"]["ok"] is True and dres["highbit_prod_diag"]["ok"] is True
 assert dres["hea_surv_analytic"]["ok"] is True
 _eh5 = dres["exact_deg72_h5"]
 assert _eh5["usable_implies_surv"] is True and _eh5["disjoint_c12_orbits"] == 100 and _eh5["live_replays"] >= (6 if part == "all" else 2) and _eh5["pooled_obs_over_mu_pred_2core"] >= 1.0
 assert dres["pm_log2_rels"] <= _op["log2_Up"] + math.log2(_op["exc"]) and dres["pm_log2_rels"] <= _cl["log2_rels_2se_live"]
 assert abs(dres["pm_main_lb_f"] - _op["lf"]) < 1e-9 and abs(dres["pm_main_lb"] - _op["lg"]) < 1e-9
 assert dres["pm_main_nf"] <= _cl["nf_cond_mean"] and dres["pm_main_ng"] <= _cl["ng_cond_mean"]
 assert abs(dres["pm_alpha_f"] - _op["alpha_f"]) < 1e-9 and abs(dres["pm_alpha_g"] - _op["alpha_g_total"]) < 1e-9
 assert dres["pm_pair_setup"] == float(_dp_kern["pair_setup"]) == 2.0 ** 35.2
 assert dres["pm_lll_node"] == float(_dp_kern["lll_node"]) == 2.0 ** 40.0
 assert dres["pm_lll_boot"] == float(_dp_kern["lll_boot"]) == 2.0 ** 47.0
 assert dres["pm_bound_cell"] == float(_dp_kern["bound_cell"]) == 4400.0 and dres["pm_sieve_hit"] == 1400.0
 assert dres["pm_ret_mult_f"] <= _cl["ret_mult_f_min"] and dres["pm_ret_mult_g"] <= _cl["ret_mult_g_min"]
 assert dres["ntt3_stage2_ecm"]["ok"] is True and dres["e16_fvor_live"]["unit_dup_violations"] == 0
 _t16_ecm = dres["t16_edwards_ecm"]
 assert (
  dres["eager_all_100pct"] is True and dres["linear_miss_retention"] is True
  and dres["f4_counterexample_rejected"] is True and dres["t16_edwards_ecm_ok"] is True
  and _t16_ecm["torsion"] == 16 and _t16_ecm["nl_rank_ge1"] is True and _t16_ecm["tors_z2xz8_q"] is True
  and _t16_ecm["emp_composite_mode"] is True and _t16_ecm["stage2_tree_bisect_ok"] is True
  and _t16_ecm["stage2_zero_leaf_skip_ok"] is True and _t16_ecm["stage2_raw_zero_leaf_ok"] is True
  and _t16_ecm["stage2_subfactor_zero_ok"] is True and _t16_ecm["multicurve_miss_tail_ok"] is True
  and _t16_ecm["emp_p1"] >= _t16_ecm["model_p1_lo"]
  and _t16_ecm["emp_balanced_p1"] >= _t16_ecm["model_p1_lo"]
  and _t16_ecm["emp_tri_p1"] >= _t16_ecm["model_p1_lo"]
 )
 import q15_descent16 as _DS
 for _u_chk in (1032.0 / 1024.0, 2.0, 3080.0 / 1024.0, 3306.0 / 1024.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 14.0 + 1.0 / 1024.0, 14.5, 15.0):
  _re_u = rho_exact(_u_chk)
  assert _DS.rho_lo(_u_chk) < _re_u < _DS.rho_hi(_u_chk), f"Dickman enclosure failed at u={_u_chk}"
 _st_q15 = json.loads(lzma.decompress(base64.b85decode("".join((ROOT / "q15_stats.b85").read_text().split()))).decode())
 _Pm_c = _DS.P()
 _DS.load_stats(_st_q15, _Pm_c)
 _Pm_c.MAIN_LB_F = _op["lf"]
 _Pm_c.MAIN_LB = _op["lg"]
 _Pm_c.MAIN_UF = 545.0 / _Pm_c.MAIN_LB_F
 _Pm_c.MAIN_UG = 665.0 / _Pm_c.MAIN_LB
 for _mf in (600, 900, 1200):
  assert _DS.miss_side(math.ceil(4 * _mf / _Pm_c.LX), _Pm_c.LOG2_RELS, _Pm_c.MAIN_UF, _Pm_c.LX, _Pm_c.MAIN_LB_F, _Pm_c.ORBIT, _Pm_c.P_SIEVE_F, _Pm_c.RET_MULT_F) >= _DS.miss_side(math.ceil(4 * _mf / _Pm_c.LX), _cl["log2_rels_2se_live"], _cl["nf_cond_mean"] / _op["lf"], _Pm_c.LX, _op["lf"], _Pm_c.ORBIT, data["collection"]["p_f"], _cl["ret_mult_f_min"])
 for _mg in (700, 1200, 2590):
  assert _DS.miss_side(math.ceil(4 * _mg / _Pm_c.LX), _Pm_c.LOG2_RELS, _Pm_c.MAIN_UG, _Pm_c.LX, _Pm_c.MAIN_LB, _Pm_c.ORBIT, _Pm_c.P_SIEVE_G, _Pm_c.RET_MULT_G) >= _DS.miss_side(math.ceil(4 * _mg / _Pm_c.LX), _cl["log2_rels_2se_live"], _cl["ng_cond_mean"] / _op["lg"], _Pm_c.LX, _op["lg"], _Pm_c.ORBIT, min(data["collection"]["p_g_shells"]), _cl["ret_mult_g_min"])
 assert _Pm_c.P_SIEVE_G <= min(data["collection"]["p_g_shells"])
 for _lm in (0.1, 1.0, 2.0, 5.0, 10.0, 18.46, 32.0, 100.0, 250.0, 600.0, 700.0, 933.0, 1000.0):
  _lp = 0.75 * _lm
  _p_thin = math.exp(-_lp) * (1.0 + _lp) if _lp < 700.0 else 1e-290
  _p_p3 = math.exp(-_lm) * (1.0 + _lm + 0.5 * _lm * _lm + (_lm ** 3) / 6.0) if _lm < 700.0 else 1e-290
  _ex_thin = math.exp(-_lp) * (1.0 + _lp) if _lp < 709.0 else 0.0
  _ex_p3 = math.exp(-_lm) * (1.0 + _lm + 0.5 * _lm * _lm + (_lm ** 3) / 6.0) if _lm < 709.0 else 0.0
  assert max(_p_p3, _p_thin) >= _ex_thin and max(_p_p3, _p_thin) >= _ex_p3
 dres["missing_exact_thin_dom_rho25"] = True

 _rows_d = dict(_led["rows"])
 assert dres["two_log_units"] <= _rows_d["both_individual_logs_exp"]
 with localcontext() as ctx:
  ctx.prec = 80
  two = Decimal(2); ln2 = two.ln()
  p2 = lambda v: (Decimal(str(v)) * ln2).exp()
  lg2 = lambda v: v.ln() / ln2
  exp_core = sum(p2(v) for _, v in _led["rows"])
  exp_core_l2 = lg2(exp_core)
  assert exp_core_l2 <= Decimal(_led["exp_run_claim"])
  t_run_l2 = Decimal(str(_led["t_run_log2"]))
  assert t_run_l2 == Decimal("124.0")
  p_meter = p2(_led["exp_run_claim"]) / p2(t_run_l2)
  p_meter_110 = Decimal("1.1") * p_meter
  assert p_meter <= Decimal(_led["p_meter_upper"])
  order_d = Decimal(str(R_BLS))
  assert data.get("track") == "bls12-381-g1-ecdlp-plain" and data.get("unit_ir") == UNIT and data.get("g1_unit_ir") == 2000.0
  _g1u = int(data["g1_unit_ir"])
  assert int(UNIT) % _g1u == 0
  g1s = Decimal(int(UNIT) // _g1u)
  assert g1s == Decimal(10)
  fail_main = two ** -66 + two ** -80 + Decimal(2) / order_d
  fail_desc = p_meter + Decimal(4) / order_d
  fail_desc_110 = p_meter_110 + Decimal(4) / order_d
  succ_core = (Decimal(1) - fail_main) * (Decimal(1) - fail_desc)
  succ_head = (Decimal(1) - fail_main) * (Decimal(1) - fail_desc_110) - two ** -48
  lo = F(_led["p_floor"])
  p_tgt = F(_led.get("p_cold_target", "0.391"))
  assert F(str(succ_head)) >= lo > p_tgt >= F(39, 100)
  K = int(_led["K_gate"])
  assert 0 < K < (1 << 64)
  assert K * lo >= p_tgt * (1 << 64) > (K - 1) * lo
  _wr = p2(_led["exp_run_claim"]) + p2(_led["wrapper_log2"])
  ew_p = lg2(Decimal(1) + Decimal(K) * _wr / (two ** 64))
  assert ew_p < Decimal(_led["claim_p"])
  ew_g1 = lg2(Decimal(1) + Decimal(K) * g1s * _wr / (two ** 64))
  assert ew_g1 < Decimal(_led["claim_g1"])
  _wh = p2(t_run_l2) + p2(_led["wrapper_log2"])
  w_hard_l2 = lg2(_wh)
  w_hard_g1_l2 = lg2(g1s * _wh)
  manifest_ok = False; track_claim = None
  for mp in (ROOT / "submission.json", HERE / "submission.json"):
   if mp.is_file():
    mf = json.loads(mp.read_text())
    sc = mf["claims"][0]["score"]
    sn, sd = int(sc["numerator"]), int(sc["denominator"])
    assert math.gcd(sn, sd) == 1
    declared = F(sn, sd)
    ev_files = mf["claims"][0].get("evidence", [])
    ex_files = [ex["source"] for ex in mf.get("experiments", []) if "source" in ex]
    if declared == F(_led["claim_g1"]) and ev_files and all((mp.parent / ef).is_file() for ef in ev_files + ex_files):
     track_claim = _led["claim_g1"]
     assert ("2^" + track_claim) in mf["title"]
     assert ("2^" + _led["exp_run_claim"]) in mf["claims"][0]["statement"]
     assert _led["claim_p"] in (mp.parent / "proof.md").read_text()
     manifest_ok = True
    break
  assert manifest_ok and track_claim is not None
  manifest_subbounds_ok = bool(
   manifest_ok
   and _cl["coll_u"] <= _rows_d["relation_collection"]
   and _cl["build_u"] <= _rows_d["replay_a_tables"]
   and _cl["la_u"] <= _rows_d["linear_algebra"]
   and _cl["gen_u"] <= _rows_d["matrix_generator"]
   and dres["two_log_units"] <= _rows_d["both_individual_logs_exp"]
   and _cl["band_u"] <= _rows_d["band_split"]
   and _cl["ecm_u"] <= _rows_d["ecm_split_expected"]
   and exp_core_l2 <= Decimal(_led["exp_run_claim"])
   and ew_p < Decimal(_led["claim_p"])
   and ew_g1 < Decimal(_led["claim_g1"])
  )
  assert manifest_subbounds_ok
  _q10 = Decimal("1e-10"); _up = lambda d: str(d.quantize(_q10, rounding=ROUND_CEILING))
 out["ledger_gate_claim"] = dict(
  exp_run_log2=str(exp_core_l2)[:18], exp_run_claim=_led["exp_run_claim"],
  t_run_log2=str(t_run_l2), w_hard_log2=str(w_hard_l2)[:18], w_hard_g1_log2=str(w_hard_g1_l2)[:18],
  p_meter=_up(p_meter), p_meter_110=_up(p_meter_110),
  p_descent_fail="0.0",
  succ_core=str(succ_core)[:10], succ_head_110_wrap=str(succ_head)[:10],
  p_floor=_led["p_floor"], p_cold_target=str(p_tgt), K_gate=K, ew_p_log2=str(ew_p)[:18],
  claim_p=_led["claim_p"], ew_g1_log2=str(ew_g1)[:18], claim_g1=track_claim,
  track_claim=track_claim, manifest_subbounds_ok=manifest_subbounds_ok,
 )
 out["doc_code_sync"] = check_doc_code_sync(ROOT, dict(_led, claim_p=track_claim), _op, _cl, dres=dres, p2_live=out["pool_2adic_live"])
 ck["doc_code_sync"] = True
 ck["ledger_gate_claim"] = True
 res = dict(checks_ok=all(ck.values()), checks=ck, summary=out)
 print(json.dumps(res, separators=(",", ":")))
 if not res["checks_ok"]:
  raise SystemExit(1)


if __name__ == "__main__":
 main()
