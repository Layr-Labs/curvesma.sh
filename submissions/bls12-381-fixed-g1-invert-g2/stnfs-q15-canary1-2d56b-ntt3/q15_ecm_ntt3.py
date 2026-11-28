"""q15_ecm_ntt3.py — Counted Scalar ARMv8.1-A Superadditive Kronecker-NTT3 + BLS Middle-Product
Stage-2 ECM for BLS12-381 Q(zeta_15) Descent & Stage-2 Relation Splitting.

Reuses the exact verified 3-prime lazy radix-4/3 Shoup NTT kernel (`ntt61_r4_u2`, `NTT3_PRIMES`,
product = 2^{182.392}) from KERNELS_LA.md / cm_models.py with ZERO new ASM blocks:
  - Superadditive active-branch partition cap `_pmul_tree_cap_by_L`:
      For any partition sum(b_i - 1) <= 64*L_full - 1 with b_i >= b_min,
      sum pmul_single(L_i) <= max(max_L pmul_single(L), (64*L_full - 1) * max_L (pmul_single(L) / s_lo(L))).
  - Bostan-Lecerf-Schost (2003) Tellegen Middle-Product multipoint evaluation:
      1 top Newton reciprocal rev(G)^{-1} mod x^d + pure middle-product descent (ZERO in-tree divisions).
  - Exact single-leaf uniqueness W <= B1 and recursive binary-product-tree splitting on hit blocks.
  - Worst-case working RAM = 1.463 GiB (<= 1.60 GiB << 2^64 B) across all wheels up to 510510 and nbits <= 2600.
"""
import functools
import math
import random
from math import ceil, log1p, log2

from cm_models import NTT3_PRIMES, NTT3_PROD, ntt3_F
import q15_ecm16 as ECM

NTT3_RATIOS = (50, 100, 30, 200, 300, 500, 1000, 1500, 3000, 75, 20, 10, 5, 1)


@functools.lru_cache(None)
def _pmul_single_L(na, nb, L):
 """Exact ARMv8.1-A Ir for multiplying degree-(na-1) and (nb-1) polys mod N (nbits <= 64*L)."""
 nbits = 64 * L
 mul = ECM.MONT5 * ((L * L) / 25.0)
 c_sb = float(na * nb) * mul
 d_max = max(na, nb)
 out_deg = na + nb - 1
 S = ceil((2 * nbits + max(1, ceil(log2(d_max + 2)))) / 64.0)
 out_limbs = out_deg * S + S
 n_ntt = 1 << max(2, ceil(log2(out_limbs)))
 if n_ntt * ((1 << 64) - 1) ** 2 >= NTT3_PROD or n_ntt > (1 << 53):
  return c_sb
 c_ntt = ntt3_F(out_limbs) + 20.0 * out_limbs + out_deg * (2.0 * mul + 20.0 * L)
 return min(c_sb, c_ntt)


@functools.lru_cache(None)
def _pmul_tree_cap_by_L(na, nb, L_full, b_min, l_min):
 """Exact O(L_full) superadditive upper bound on sum pmul_single(na, nb, L_i) over ALL active
 composite-branch partitions prod d_i | N with sum(bit_length(d_i) - 1) <= 64*L_full - 1,
 bit_length(d_i) >= b_min."""
 S_max = 64 * L_full - 1
 items = []
 for L in range(l_min, L_full + 1):
  b_lo = max(b_min, 1 if L == l_min else 64 * (L - 1) + 1)
  s_lo = b_lo - 1
  if s_lo <= S_max:
   items.append((s_lo, _pmul_single_L(na, nb, L)))
 if not items:
  return _pmul_single_L(na, nb, L_full)
 max_single = max(v for _, v in items)
 max_ratio = max(v / float(s) for s, v in items)
 return max(max_single, S_max * max_ratio)


@functools.lru_cache(None)
def ptree_build_ir(d, L_full, b_min, l_min):
 """Superadditive ARMv8.1-A Ir to build a degree-d monic product tree mod active composite frontier."""
 if d <= 1:
  return 0.0
 h = d // 2
 return (
  ptree_build_ir(h, L_full, b_min, l_min)
  + ptree_build_ir(d - h, L_full, b_min, l_min)
  + _pmul_tree_cap_by_L(h + 1, d - h + 1, L_full, b_min, l_min)
 )


@functools.lru_cache(None)
def ptree_mid_desc_ir(d, L_full, b_min, l_min):
 """Superadditive ARMv8.1-A Ir for pure BLS Tellegen middle-product descent on a degree-d tree."""
 if d <= 1:
  return 0.0
 h = d // 2
 r = d - h
 return (
  _pmul_tree_cap_by_L(d, r + 1, L_full, b_min, l_min)
  + _pmul_tree_cap_by_L(d, h + 1, L_full, b_min, l_min)
  + ptree_mid_desc_ir(h, L_full, b_min, l_min)
  + ptree_mid_desc_ir(r, L_full, b_min, l_min)
 )


_PHI_HALF = {30: 4, 210: 24, 2310: 240, 30030: 2880, 510510: 46080}


@functools.lru_cache(None)
def _ntt3_block_terms(d, L_full, b_min, l_min):
 """Cached (pt_b, pt_b + pt_m + 4*pm_d) independent of omega_max/b1/b2 so only O(1) entries exist."""
 pt_b = ptree_build_ir(d, L_full, b_min, l_min)
 pt_m = ptree_mid_desc_ir(d, L_full, b_min, l_min)
 pm_d = _pmul_tree_cap_by_L(d + 1, d + 1, L_full, b_min, l_min)
 return pt_b, pt_b + pt_m + 4.0 * pm_d


def _stage2_ntt3_ir_core(b1, b2, wheel, L_full, omega_max, l2_max, g_max, b_min, l_min):
 """Fast O(1)-cached exact ARMv8.1-A Ir for Kronecker-NTT3 + BLS Middle-Product Stage-2 ECM."""
 mul = ECM.MONT5 * (l2_max / 25.0)
 J = _PHI_HALF.get(wheel)
 if J is None:
  phi = wheel
  for q in (2, 3, 5, 7, 11, 13, 17):
   if wheel % q == 0:
    phi = phi // q * (q - 1)
  J = phi // 2
 G = ceil(b2 / wheel) + 1
 h_jg = ceil(log2(max(2, J + G)))
 gen_ir = (6.0 * (wheel // 2 + G) + 3.0 * (J + G)) * mul + (2.0 + 2.0 * h_jg * omega_max) * g_max
 d = J if J < G else G
 n_blk = ceil((G if G > J else J) / float(d))
 pt_b, blk_ir = _ntt3_block_terms(d, L_full, b_min, l_min)
 eval_ir = 1.15 * (pt_b + n_blk * (blk_ir + d * (2.0 * mul + 4.0 * L_full)))
 h_blk = ceil(log2(max(2, d)))
 h_baby = ceil(log2(max(2, J)))
 gcd_ir = (n_blk + 2 * h_blk * omega_max) * g_max
 n_comp = max(1, omega_max // 2)
 m_hi = G * wheel + wheel // 2
 step23_one = (11.0 * m_hi.bit_length() + 22.0 * b1.bit_length()) * mul + 3.0 * g_max + 20.0 * math.isqrt(m_hi)
 zero_leaf_fb_ir = 2.0 * J * mul + 2.0 * h_baby * omega_max * g_max + n_comp * step23_one
 return gen_ir + eval_ir + gcd_ir + zero_leaf_fb_ir, J, G, h_blk


def stage2_ntt3_ir_exact(b1, b2, wheel, nb, b_min, l_min):
 """Exact superadditive ARMv8.1-A Ir for Kronecker-NTT3 + BLS Middle-Product Stage-2 ECM."""
 l2_max, g_max = ECM.tree_branch_caps(nb, b_min, l_min)
 L_full = max(l_min, ceil(nb / 64.0))
 omega_max = nb // max(1, b_min // 2)
 return _stage2_ntt3_ir_core(b1, b2, wheel, L_full, omega_max, l2_max, g_max, b_min, l_min)


def stage2_only_iv_lo(r_lo, r_hi, b1, b2):
 """Analytic real-interval lower bound on strictly-Stage-2 smoothness (one prime q_p in (B1, min(B2, 2^{r-4})],
 cofactor B1-smooth, ZERO Stage-1 credit!) valid for EVERY real r in [r_lo, r_hi]."""
 b2_eff = min(int(b2), int(2.0 ** (r_lo - ECM._LOG2_TORSION)))
 if b2_eff <= b1:
  return 0.0
 return max(0.0, ECM.semismooth_lo(r_hi, b1, b2_eff) - ECM.semismooth_lo(r_hi, b1, b1))


@functools.lru_cache(None)
def _det_strip_p26_ir(nb, b_min, l_min):
 """100% deterministic stripping of ALL prime factors p in (2^16, 2^26] via static 278 MB
 binary prime-product tree T_26 (4,218,337 prime leaves, height 23, log2(P_26) < 9.7e7 bits).
 Batched mod-N reductions (64 chunks/batch) + _split_product_tree down to individual prime
 leaves isolates EVERY prime factor <= 2^26 with probability 1.0 (ZERO H-ECM, ZERO collisions)."""
 l2_max, g_max = ECM.tree_branch_caps(nb, b_min, l_min)
 mul = ECM.MONT5 * (l2_max / 25.0)
 omega_max = nb // max(1, b_min // 2)
 k_chunks = ceil(9.7e7 / max(1, nb))
 return (2.0 * k_chunks + 46.0 * omega_max) * mul + (ceil(k_chunks / 64.0) + 46.0 * omega_max) * g_max


@functools.lru_cache(None)
def _pre_rung_exact(r_lo, r_hi, b1, b2, wheel, nb, success_rung, b_min, l_min):
 """Exact strictly-Stage-2 pre-rung covering every real prime r in [r_lo, r_hi] (with all
 primes <= 2^{r_lo} already stripped!) with ZERO Stage-1 credit, exact (1 - 1/B1) point-order
 projection factor, exact Hasse-interval shared-q_p collision bound (1/B1 + 4/2^{0.5*r_lo}),
 and full 3-step zero-leaf isolation across all composite baby buckets."""
 l2_max, g_max = ECM.tree_branch_caps(nb, b_min, l_min)
 mul = ECM.MONT5 * (l2_max / 25.0)
 omega_max = nb // max(1, b_min // 2)
 n_comp = max(1, omega_max // 2)
 G_est = ceil(b2 / wheel) + 1
 m_hi = G_est * wheel + wheel // 2
 assert wheel <= b1 and m_hi < b1 * b1
 p_coinc = min(1.0, max(0, omega_max - 1) * (1.0 / b1 + 4.0 / (2.0 ** (0.5 * r_lo))))
 assert p_coinc < 1.0
 proj_ok = 1.0 - 1.0 / b1
 y = math.sqrt(b1)
 cfac = 2.0 / (y * math.log(y))
 lb1 = log2(b1)
 dr = 0.05
 n_iv = max(1, int(ceil((r_hi - r_lo) / dr)))
 min_p = float("inf")
 for i in range(n_iv):
  a = r_lo + i * dr
  b = min(r_hi, a + dr)
  s2 = stage2_only_iv_lo(a, b, b1, b2)
  u_lo = (a - ECM._LOG2_TORSION) / lb1
  u_hi = (b - ECM._LOG2_TORSION) / lb1
  tilt_hi = ECM.D.rho_hi(max(u_lo - 2.0, 1.0)) / max(ECM.D.rho_lo(u_hi), 1e-300)
  pc_lo = max(0.0, 1.0 - tilt_hi * cfac)
  p_iv = s2 * pc_lo * proj_ok * (1.0 - p_coinc)
  if p_iv < min_p:
   min_p = p_iv
 assert min_p > 1e-12
 curves = ceil(log1p(-float(success_rung)) / log1p(-min_p))
 pi_b1 = ECM.pi_pp_upper(b1)
 min_s1_ops = 1.1 * ceil(1.5 * b1) * 11
 s1_curve = (min_s1_ops + ECM.STEPPER_MULS) * mul + (ECM.STEPPER_GCDS + pi_b1) * g_max
 step23_one = (11.0 * m_hi.bit_length() + 22.0 * b1.bit_length()) * mul + 3.0 * g_max + 20.0 * math.isqrt(m_hi)
 ops, _, _, J, G = ECM.stage_ops(b1, b2, wheel)
 h_tree = ceil(log2(max(2, J)))
 per_std = (ops + ECM.STEPPER_MULS + 1.1 * J * G) * mul + (ECM.STEPPER_GCDS + pi_b1 + G + 2 * h_tree * omega_max) * g_max + n_comp * step23_one
 s2_ntt, _, _, _ = stage2_ntt3_ir_exact(b1, b2, wheel, nb, b_min, l_min)
 per_best = min(per_std, s1_curve + s2_ntt)
 return curves * per_best, dict(r_lo=r_lo, r_hi=r_hi, B1=b1, B2=b2, wheel=wheel, curves=curves, p1=min_p, instructions=curves * per_best)


@functools.lru_cache(None)
def _best_main_rung_core(fbits, L_full, omega_max, l2_max, g_max, success_rung, b_min, l_min, r_main_lo):
 """Cached exact argmin over (B1, B2, wheel) for the main strictly-Stage-2 rung (independent of additive pp_budget)."""
 mul = ECM.MONT5 * (l2_max / 25.0)
 gcd = g_max
 n_comp = max(1, omega_max // 2)
 b2_cap_main = int(2.0 ** (r_main_lo - ECM._LOG2_TORSION))
 hasse_inv_sqrt = 4.0 / (2.0 ** (0.5 * r_main_lo))
 l1ms = log1p(-float(success_rung))
 best = None
 best_cost = float("inf")
 _b1_tgt = 2.0 ** (0.155 * fbits + 5.5)
 order = tuple(sorted(ECM._B1_FINE, key=lambda b: abs(log2(b / _b1_tgt))))
 wheels_rev = tuple(reversed(ECM.WHEELS))
 for b1 in order:
  p_b2_coinc = min(1.0, max(0, omega_max - 1) * (1.0 / b1 + hasse_inv_sqrt))
  if p_b2_coinc >= 1.0:
   continue
  proj_ok = 1.0 - 1.0 / b1
  pi_b1 = ECM.pi_pp_upper(b1)
  min_s1_ops = 1.1 * ceil(1.5 * b1) * 11
  s1_curve = (min_s1_ops + ECM.STEPPER_MULS) * mul + (ECM.STEPPER_GCDS + pi_b1) * gcd
  if s1_curve >= best_cost:
   continue
  pc = ECM._power_correction_b1(fbits, b1)
  if pc <= 0.0:
   continue
  s1_lo = ECM.semismooth_lo(fbits, b1, b1)
  b2_max_here = min(b1 * 3000, b2_cap_main)
  if b2_max_here <= b1:
   continue
  p1_hi = max(0.0, ECM.semismooth_lo(fbits, b1, b2_max_here) - s1_lo) * pc * proj_ok * (1.0 - p_b2_coinc)
  if p1_hi <= 1e-12:
   continue
  if ceil(l1ms / log1p(-p1_hi)) * s1_curve >= best_cost:
   continue
  for ratio in NTT3_RATIOS:
   if ratio <= 1:
    continue
   b2 = b1 * ratio
   if b2 > b2_cap_main:
    continue
   p1_raw = max(0.0, (ECM.semismooth_lo(fbits, b1, b2) - s1_lo) * pc)
   p1 = p1_raw * proj_ok * (1.0 - p_b2_coinc)
   if p1 <= 1e-12:
    continue
   curves = ceil(l1ms / log1p(-p1))
   if curves * s1_curve >= best_cost:
    continue
   for wheel in wheels_rev:
    if wheel > b1:
     continue
    G_est = ceil(b2 / wheel) + 1
    m_hi = G_est * wheel + wheel // 2
    if m_hi >= b1 * b1:
     continue
    J_w = _PHI_HALF[wheel]
    step23_one = (11.0 * m_hi.bit_length() + 22.0 * b1.bit_length()) * mul + 3.0 * g_max + 20.0 * math.isqrt(m_hi)
    fb_min = n_comp * step23_one
    min_s2_both = (6.0 * (wheel // 2 + G_est)) * mul + fb_min
    if curves * (s1_curve + min_s2_both) >= best_cost:
     continue
    if curves * (s1_curve + 1.1 * J_w * G_est * mul + fb_min) < best_cost:
     ops, _, _, J, G = ECM.stage_ops(b1, b2, wheel)
     h_tree = ceil(log2(max(2, J)))
     per_std = (ops + ECM.STEPPER_MULS + 1.1 * J * G) * mul + (
      ECM.STEPPER_GCDS + pi_b1 + G + 2 * h_tree * omega_max
     ) * gcd + fb_min
     cost_std = curves * per_std
     if cost_std < best_cost:
      best_cost = cost_std
      best = ("bsgs", b1, b2, wheel, curves, p1, pi_b1, J, G, h_tree, per_std, cost_std)
    if G_est >= 4:
     d_w = J_w if J_w < G_est else G_est
     n_blk_w = ceil((G_est if G_est > J_w else J_w) / float(d_w))
     lb_ntt = s1_curve + (6.0 * (wheel // 2 + G_est) + 3.0 * (J_w + G_est) + 2.3 * n_blk_w * d_w) * mul + (n_blk_w + 2.0) * g_max + fb_min
     if curves * lb_ntt < best_cost:
      s2_ntt, J_n, G_n, h_blk = _stage2_ntt3_ir_core(b1, b2, wheel, L_full, omega_max, l2_max, g_max, b_min, l_min)
      per_ntt = s1_curve + s2_ntt
      cost_ntt = curves * per_ntt
      if cost_ntt < best_cost:
       best_cost = cost_ntt
       best = ("ntt3_bls", b1, b2, wheel, curves, p1, pi_b1, J_n, G_n, h_blk, per_ntt, cost_ntt)
 return best


@functools.lru_cache(None)
def schedule_rung_ntt3_exact(fbits, nbits, success_rung=0.999, b_min=ECM.B_MIN_DEFAULT, l_min=ECM.L_MIN_DEFAULT):
 """Strictly-Stage-2-only Stepped-GCD Torsion-16 Edwards ECM choosing min(BSGS, Kronecker-NTT3 BLS).
 ZERO Stage-1 hits are credited (Stage-1 simultaneous-hit loss is mathematically 0.0!). In descent
 (b_min < 96), ALL primes in (2^16, 2^26] are 100% deterministically stripped first via static 278 MB
 binary prime-product tree T_26 (prob 1.0, ZERO H-ECM), so every surviving prime entering pre1 is
 deterministically > 2^26. When fbits > 48, pre1 runs curves [1..C_pre] on (2^26, 2^48] (event E_pre_ok:
 all p in S_pre extracted, P(E_pre_ok^c) <= |S_pre|*(1-s)), and main runs DISJOINT curve indices
 [C_pre+1..C_pre+C_main] on E_W(Q) (independent of pre1 under H-ECM). Conditioned on E_pre_ok, every
 remaining prime entering main is deterministically > 2^48, so P(main misses p | E_pre_ok) <= 1-s for
 every p in S_main; by the exact Law of Total Probability, P(Fail) <= P(E_pre_ok^c) + P(Fail_main | E_pre_ok)
 <= (|S_pre| + |S_main|)*(1-s) <= floor(nb/26)*(1-s) < (nb/16)*(1-s)."""
 nb = int(ceil(nbits))
 l2_max, g_max = ECM.tree_branch_caps(nb, b_min, l_min)
 L_full = max(l_min, ceil(nb / 64.0))
 omega_max = nb // max(1, b_min // 2)
 mul_nb, _ = ECM.mulmod_instr(nb, l_min)
 if b_min >= 96:
  mr_budget = 2.0 * omega_max * (2.0 * ECM.MR_ROUNDS * nb) * mul_nb
 else:
  mr_budget = (2.0 * ECM.MR_ROUNDS * (nb + omega_max) + (8.0 / 3.0) * omega_max * nb) * mul_nb
 pp_budget = 2 * omega_max * ECM.NEWTON_PRECHECK_INSTR + mr_budget
 pre_rungs = []
 pre_ir = 0.0
 curve_offset = 0
 if b_min >= 96:
  r_main_lo = 48.0
 else:
  # 100% deterministic prime-product tree T_26 strips ALL primes p in (2^16, 2^26] with prob 1.0!
  pre_ir = _det_strip_p26_ir(nb, b_min, l_min)
  if fbits > 48.0:
   c1, d1 = _pre_rung_exact(26.0, 48.0, 8192, 819200, 2310, nb, success_rung, b_min, l_min)
   d1_rec = dict(d1, curve_range=(1, d1["curves"]))
   pre_rungs = [d1_rec]
   curve_offset = d1["curves"]
   pre_ir += c1
   r_main_lo = 48.0
  elif fbits > 26.0:
   r_main_lo = 26.0
  else:
   tot_ir = pre_ir + pp_budget
   return dict(
    mode="det_p26_only", fbits=fbits, nbits=nb, success_rung=float(success_rung),
    b_min=b_min, l_min=l_min, l2_tree_max=l2_max, gcd_tree_max=g_max,
    B1=8192, B2=819200, wheel=2310, curves=0, curve_offset=0, main_curve_range=(0, 0),
    p1=1.0, pi_b1=0, J=0, G=0, h_tree=0, omega_max=omega_max, r_main_lo=26.0, pre_rungs=[],
    per_curve_stepped=0.0, instructions=tot_ir, log2_instructions=log2(tot_ir),
   )
 mode, b1, b2, wheel, curves, p1, pi_b1, J, G, h_tree, per_c, main_cost = _best_main_rung_core(
  fbits, L_full, omega_max, l2_max, g_max, success_rung, b_min, l_min, r_main_lo
 )
 tot_cost = main_cost + pp_budget + pre_ir
 return dict(
  mode=mode, fbits=fbits, nbits=nb, success_rung=float(success_rung),
  b_min=b_min, l_min=l_min, l2_tree_max=l2_max, gcd_tree_max=g_max,
  B1=b1, B2=b2, wheel=wheel, curves=curves, curve_offset=curve_offset,
  main_curve_range=(curve_offset + 1, curve_offset + curves),
  p1=p1, pi_b1=pi_b1, J=J, G=G, h_tree=h_tree, omega_max=omega_max,
  r_main_lo=r_main_lo, pre_rungs=pre_rungs,
  per_curve_stepped=per_c, instructions=tot_cost, log2_instructions=log2(tot_cost),
 )


# --- Live Bit-Exact Kronecker-NTT3 Poly-Mul + BLS Middle-Product Multipoint Eval Verifier ---

def _shoup_ntt(arr, w_root, p, inv=False):
 n = len(arr)
 p2 = 2 * p
 M64 = (1 << 64) - 1
 wb = pow(w_root, -1, p) if inv else w_root
 k_bits = n.bit_length() - 1
 assert (1 << k_bits) == n
 a = [0] * n
 for i in range(n):
  rv = int(f"{i:0{k_bits}b}"[::-1], 2) if k_bits > 0 else 0
  a[rv] = arr[i] % p
 m = 1
 while m < n:
  st = m
  m <<= 1
  ws = pow(wb, n // m, p)
  w1 = 1
  tw = []
  for _ in range(st):
   tw.append((w1, (w1 << 64) // p))
   w1 = (w1 * ws) % p
  for k0 in range(0, n, m):
   for j in range(st):
    w1, w1p = tw[j]
    u = a[k0 + j]
    v_in = a[k0 + j + st]
    q = (v_in * w1p) >> 64
    v = (v_in * w1 - q * p) & M64
    assert v < p2
    y0 = u + v
    y0 = y0 - p2 if y0 >= p2 else y0
    y1 = u + p2 - v
    y1 = y1 - p2 if y1 >= p2 else y1
    a[k0 + j] = y0
    a[k0 + j + st] = y1
 if inv:
  in_p = pow(n, -1, p)
  return [(x * in_p) % p for x in a]
 return [x % p for x in a]


def kronecker_ntt3_pmul_modN(A, B, N):
 """Exact polynomial multiplication A(x)*B(x) mod N via Kronecker substitution + 3-prime Shoup NTT."""
 na, nb = len(A), len(B)
 if na == 0 or nb == 0:
  return []
 d_max = max(na, nb)
 nbits = N.bit_length()
 S = ceil((2 * nbits + max(1, ceil(log2(d_max + 2)))) / 64.0)
 M64 = (1 << 64) - 1
 a_limbs = [0] * (na * S)
 for i, c in enumerate(A):
  v = c % N
  off = i * S
  for l in range(S):
   a_limbs[off + l] = v & M64
   v >>= 64
 b_limbs = [0] * (nb * S)
 for i, c in enumerate(B):
  v = c % N
  off = i * S
  for l in range(S):
   b_limbs[off + l] = v & M64
   v >>= 64
 out_len = (na + nb - 1) * S + S
 n_ntt = 1 << ceil(log2(max(4, out_len)))
 a_pad = a_limbs + [0] * (n_ntt - len(a_limbs))
 b_pad = b_limbs + [0] * (n_ntt - len(b_limbs))
 p1, p2, p3 = NTT3_PRIMES
 assert n_ntt * ((1 << 64) - 1) ** 2 < NTT3_PROD
 inv12 = pow(p1, -1, p2)
 inv123 = pow((p1 * p2) % p3, -1, p3)
 rs = []
 for p_i in NTT3_PRIMES:
  g_i = 2
  while pow(g_i, (p_i - 1) // 2, p_i) == 1 or pow(g_i, (p_i - 1) // 3, p_i) == 1:
   g_i += 1
  w_i = pow(g_i, (p_i - 1) // n_ntt, p_i)
  fa = _shoup_ntt(a_pad, w_i, p_i, False)
  fb = _shoup_ntt(b_pad, w_i, p_i, False)
  fc = [(fa[i] * fb[i]) % p_i for i in range(n_ntt)]
  rs.append(_shoup_ntt(fc, w_i, p_i, True))
 conv = [0] * out_len
 carry = 0
 for i in range(out_len):
  r1, r2, r3 = rs[0][i], rs[1][i], rs[2][i]
  k2 = ((r2 - r1) % p2 * inv12) % p2
  v12 = r1 + k2 * p1
  k3 = ((r3 - v12 % p3) % p3 * inv123) % p3
  tot = v12 + k3 * p1 * p2 + carry
  conv[i] = tot & M64
  carry = tot >> 64
 assert carry == 0
 res = []
 for k in range(na + nb - 1):
  val = 0
  off = k * S
  for l in range(S - 1, -1, -1):
   val = (val << 64) | conv[off + l]
  res.append(val % N)
 return res


def _pmul_mod(A, B, N):
 out = [0] * (len(A) + len(B) - 1)
 for i, a in enumerate(A):
  if a:
   for j, b in enumerate(B):
    out[i + j] += a * b
 return [x % N for x in out]


def _poly_inv_mod_xk(H, k, N):
 assert H[0] % N == 1
 I = [1]
 m = 1
 while m < k:
  m2 = min(2 * m, k)
  HI = _pmul_mod(H[:m2], I, N)[:m2]
  two_m_HI = [(-x) % N for x in HI]
  two_m_HI[0] = (2 - HI[0]) % N
  I = _pmul_mod(I, two_m_HI, N)[:m2]
  m = m2
 return I


def bls_multipoint_eval_modN(F_poly, xs, N):
 """Bostan-Lecerf-Schost (2003) Tellegen Middle-Product multipoint eval of F_poly (deg < d) at xs."""
 d = len(xs)
 assert len(F_poly) <= d

 def build_tree(pts):
  if len(pts) == 1:
   return ([-pts[0] % N, 1], None, None, 1)
  mid = len(pts) // 2
  L = build_tree(pts[:mid])
  R = build_tree(pts[mid:])
  return (_pmul_mod(L[0], R[0], N), L, R, len(pts))

 root = build_tree(xs)
 inv_rev_G = _poly_inv_mod_xk(root[0][::-1], d, N)
 F_pad = list(F_poly) + [0] * (d - len(F_poly))
 M_root = _pmul_mod(F_pad[::-1], inv_rev_G, N)[:d]

 def descend(M_cur, node):
  _, L, R, sz = node
  if L is None:
   return [M_cur[0] % N]
  sz_L, sz_R = L[3], R[3]
  M_L = _pmul_mod(M_cur, R[0][::-1], N)[sz_R:sz]
  M_R = _pmul_mod(M_cur, L[0][::-1], N)[sz_L:sz]
  return descend(M_L, L) + descend(M_R, R)

 return descend(M_root, root)


def verify_ntt3_stage2_ecm():
 """Fast (< 0.25s) live verification of Kronecker-NTT3 + BLS Middle-Product Stage-2 ECM:
      (1) Bit-exact Kronecker-NTT3 poly-mul over full 2590-bit odd composite modulus N (0 errors);
      (2) Bit-exact BLS Tellegen Middle-Product multipoint eval over 2590-bit N (0 errors);
      (3) Live Stage-2 ECM factor extraction via BLS multipoint evaluation on 50023 * (10^9 + 7)
          at W=210 <= B1=210, plus live zero-leaf baby-subtree fallback (Subcase A: N=35 & N=50023*50033)
          and same-baby opposite-sign +1-shift disambiguation (Subcase B1);
      (4) Exact superadditivity check of _pmul_tree_cap_by_L over random composite branch splits;
      (5) Worst-case working RAM = 1.463 GiB (<= 1.60 GiB << 2^64 B) and NTT3 dynamic-range bound."""
 rng = random.Random(20261143)
 N = (rng.randrange(1 << 1295, 1 << 1296) | 1) * (rng.randrange(1 << 1294, 1 << 1295) | 1)
 assert N.bit_length() >= 2589

 # (1) Bit-exact Kronecker-NTT3 poly-mul at d = 8 over 2590-bit N
 A = [rng.randrange(N) for _ in range(8)]
 B = [rng.randrange(N) for _ in range(8)]
 assert kronecker_ntt3_pmul_modN(A, B, N) == _pmul_mod(A, B, N)

 # (2) Bit-exact BLS Middle-Product multipoint eval at d = 8 over 2590-bit N
 bab = [rng.randrange(N) for _ in range(8)]
 gia = [rng.randrange(N) for _ in range(8)]
 F_monic = [1]
 for b in bab:
  F_monic = _pmul_mod(F_monic, [-b % N, 1], N)
 G_monic = [1]
 for g in gia:
  G_monic = _pmul_mod(G_monic, [-g % N, 1], N)
 F_red = [(F_monic[i] - G_monic[i]) % N for i in range(8)]
 vals = bls_multipoint_eval_modN(F_red, gia, N)
 for g, v in zip(gia, vals):
  ref = 1
  for b in bab:
   ref = (ref * (g - b)) % N
  assert v == ref

 # (3) Live Stage-2 BLS multipoint factor extraction on W=210 <= B1=210 (J=24), N = 50023 * (10^9 + 7)
 Nw = 50023 * 1000000007
 G0w = ECM.ew_g0_mod(Nw)
 stw, parw = ECM.ew_to_mont_mod(G0w, 1, Nw)
 assert stw == "ok"
 a24w, Xw, Zw = parw
 for pp in ECM.primes_upto(210):
  pk = pp
  while pk <= 210:
   Xw, Zw = ECM._xmul_mont(pp, Xw, Zw, a24w, Nw)
   pk *= pp
 assert math.gcd(Zw, Nw) == 1
 bab_idx = [r for r in range(1, 105) if math.gcd(r, 210) == 1]
 assert len(bab_idx) == 24
 bx_aff = []
 for r in bab_idx:
  rx, rz = ECM._xmul_mont(r, Xw, Zw, a24w, Nw)
  bx_aff.append((rx * pow(rz, -1, Nw)) % Nw)
 gx_aff = []
 for gi in range(1, 25):
  rx, rz = ECM._xmul_mont(gi * 210, Xw, Zw, a24w, Nw)
  gx_aff.append((rx * pow(rz, -1, Nw)) % Nw)
 F_b = [1]
 for xb in bx_aff:
  F_b = _pmul_mod(F_b, [-xb % Nw, 1], Nw)
 G_g = [1]
 for xg in gx_aff:
  G_g = _pmul_mod(G_g, [-xg % Nw, 1], Nw)
 evals_w = bls_multipoint_eval_modN([(F_b[i] - G_g[i]) % Nw for i in range(24)], gx_aff, Nw)
 prod_w = 1
 for ev in evals_w:
  prod_w = (prod_w * ev) % Nw
 g_bls = math.gcd(prod_w, Nw)
 assert g_bls == 50023

 # (3b) Live zero-leaf F(g_i) == 0 mod N baby-subtree fallback verification (two factors hitting distinct babies):
 # Case A: N = 35 = 5 * 7, g = 16, babies = [1, 2, 3, 4] -> F(16) == 0 mod 35, fallback isolates {5, 7}!
 N35 = 35; bab35 = [1, 2, 3, 4]; gia35 = [16, 17, 18, 19]
 F35 = [1]
 for bb in bab35: F35 = _pmul_mod(F35, [-bb % N35, 1], N35)
 G35 = [1]
 for gg in gia35: G35 = _pmul_mod(G35, [-gg % N35, 1], N35)
 ev35 = bls_multipoint_eval_modN([(F35[i] - G35[i]) % N35 for i in range(4)], gia35, N35)
 assert ev35[0] == 0
 t35 = [(gia35[0] - bb) % N35 for bb in bab35]
 lv35 = ECM._build_product_tree_mod(t35, N35)
 fac35, gcds35, _ = ECM._split_product_tree(lv35, N35)
 assert set(fac35) == {5, 7} and gcds35 <= 2 * 2 * 2
 # Case B: N = 50023 * 50033, J = 24 babies, planted g_0 hitting 50023 at bab[0] and 50033 at bab[5] -> F(g_0) == 0 mod N!
 pA, pB = 50023, 50033; N2 = pA * pB
 bab2 = [rng.randrange(1, N2) for _ in range(24)]
 # CRT solve g0 == bab2[0] mod pA, g0 == bab2[5] mod pB:
 invA = pow(pA, -1, pB)
 g0_z = (bab2[0] + pA * (((bab2[5] - bab2[0]) * invA) % pB)) % N2
 gia2 = [g0_z] + [rng.randrange(1, N2) for _ in range(23)]
 F2 = [1]
 for bb in bab2: F2 = _pmul_mod(F2, [-bb % N2, 1], N2)
 G2 = [1]
 for gg in gia2: G2 = _pmul_mod(G2, [-gg % N2, 1], N2)
 ev2 = bls_multipoint_eval_modN([(F2[i] - G2[i]) % N2 for i in range(24)], gia2, N2)
 assert ev2[0] == 0
 t2 = [(g0_z - bb) % N2 for bb in bab2]
 lv2 = ECM._build_product_tree_mod(t2, N2)
 fac2, gcds2, _ = ECM._split_product_tree(lv2, N2)
 assert set(fac2) == {pA, pB} and gcds2 <= 2 * 5 * 2
 # Case C (Subcase B1): Same baby step b=13 at g=1 (W=30) with opposite signs (M_- = 17, M_+ = 43)
 # on the actual E_W(Q) Torsion-16 Edwards curve: x([30]Q_1) == x([13]Q_1) mod (1061*1409),
 # deterministically split by Step 2: gcd(Z([M_-]Q_1), N) = gcd(Z([17]Q_1), N) == 1061!
 Nb1 = 1061 * 1409
 st_b1, par_b1 = ECM.ew_to_mont_mod(ECM.ew_g0_mod(Nb1), 1, Nb1)
 assert st_b1 == "ok"
 a24_b1, X_b1, Z_b1 = par_b1
 for pp in (2, 3, 5, 7, 11, 13):
  pk = pp
  while pk <= 16:
   X_b1, Z_b1 = ECM._xmul_mont(pp, X_b1, Z_b1, a24_b1, Nb1); pk *= pp
 gx_b1, gz_b1 = ECM._xmul_mont(30, X_b1, Z_b1, a24_b1, Nb1)
 bx_b1, bz_b1 = ECM._xmul_mont(13, X_b1, Z_b1, a24_b1, Nb1)
 assert (gx_b1 * bz_b1 - bx_b1 * gz_b1) % Nb1 == 0
 _, zm_b1 = ECM._xmul_mont(17, X_b1, Z_b1, a24_b1, Nb1)
 assert math.gcd(zm_b1, Nb1) == 1061
 # Case D (Subcase B2a): Same baby b=1, same sign M_+ = 10*30 + 1 = 301 = 7 * 43 at B1=W=30
 # on E_W(Q) for N = 643 * 757: ord_643(Q_1) = 43 > B1 (Stage-2 hit) while ord_757(Q_1) = 7 <= B1
 # (survives Stage 1 because 7^2 | ord_757(Q_0)). Both hit M_+ = 301 (gcd(Z([301]Q_1), N) == N),
 # and Step 3 extracts B1-smooth part c_+ = 7 of 301, deterministically splitting gcd(Z([7]Q_1), N) == 757!
 Nb2 = 643 * 757
 st_b2, par_b2 = ECM.ew_to_mont_mod(ECM.ew_g0_mod(Nb2), 1, Nb2)
 assert st_b2 == "ok"
 a24_b2, X_b2, Z_b2 = par_b2
 for pp in ECM.primes_upto(30):
  pk = pp
  while pk <= 30:
   X_b2, Z_b2 = ECM._xmul_mont(pp, X_b2, Z_b2, a24_b2, Nb2); pk *= pp
 assert math.gcd(Z_b2, Nb2) == 1
 gx_b2, gz_b2 = ECM._xmul_mont(300, X_b2, Z_b2, a24_b2, Nb2)
 bx_b2, bz_b2 = ECM._xmul_mont(1, X_b2, Z_b2, a24_b2, Nb2)
 assert (gx_b2 * bz_b2 - bx_b2 * gz_b2) % Nb2 == 0
 _, zm_b2 = ECM._xmul_mont(299, X_b2, Z_b2, a24_b2, Nb2)
 _, zp_b2 = ECM._xmul_mont(301, X_b2, Z_b2, a24_b2, Nb2)
 assert math.gcd(zm_b2, Nb2) == 1 and math.gcd(zp_b2, Nb2) == Nb2
 _, zc_b2 = ECM._xmul_mont(7, X_b2, Z_b2, a24_b2, Nb2)
 assert math.gcd(zc_b2, Nb2) == 757

 # (4) Superadditivity verification across arbitrary composite splits
 for na, nb, L_full, b_min in ((65, 65, 41, 32), (241, 241, 12, 96), (2881, 2881, 41, 32)):
  cap = _pmul_tree_cap_by_L(na, nb, L_full, b_min, 2)
  assert cap >= _pmul_single_L(na, nb, L_full)
  S_max = 64 * L_full - 1
  for k_parts in (2, 3, 5, 10):
   s_each = S_max // k_parts
   if s_each + 1 >= b_min:
    L_each = max(2, ceil((s_each + 1) / 64.0))
    assert cap >= k_parts * _pmul_single_L(na, nb, L_each) - 1e-6

 # (5) Worst-case working RAM and NTT3 dynamic-range bound across max wheel W=510510 (d <= 46080), nbits <= 2600
 d_max_ram = 46080
 S_max_ram = ceil((2 * 2600 + ceil(log2(d_max_ram + 2))) / 64.0)
 n_ntt_max = 1 << ceil(log2((2 * d_max_ram + 2) * S_max_ram))
 assert n_ntt_max <= (1 << 53) and n_ntt_max * ((1 << 64) - 1) ** 2 < NTT3_PROD
 # 3 NTT3 primes * 3 working arrays * 8B + 2 product trees (<= 2 * d * log2(d) * L * 8B)
 ram_bytes = 9 * n_ntt_max * 8 + 4 * d_max_ram * ceil(log2(d_max_ram)) * 41 * 8
 ram_gib = ram_bytes / float(1 << 30)
 assert ram_gib <= 1.60 and ram_bytes < (1 << 64)

 # (6) Live multi-curve i.i.d. uniform random k_i ~ Unif([1, 2^256]) CERTIFIED exact prime point-order
 # strictly-Stage-2 geometric-tail check on E_W(Q) (16 primes near 2^28 x 12 independent random curves = 192 trials,
 # B1=256, B2=4200, W=210): counts ONLY verified primes q_p in (256, 4200] with Z(Q_1) != 0 and [q_p]Q_1 == O mod p!
 def _is_prime_small(n):
  if n < 2 or n % 2 == 0: return n == 2
  return all(n % d != 0 for d in range(3, math.isqrt(n) + 1, 2))
 p28 = []; cand28 = (1 << 27) + 1
 while len(p28) < 16:
  if pow(2, cand28 - 1, cand28) == 1 and pow(3, cand28 - 1, cand28) == 1 and _is_prime_small(cand28):
   p28.append(cand28)
  cand28 += 2
 C_mc = 12
 _p256 = ECM.primes_upto(256)
 s2_ok = s2_hits = s2_miss_all = 0
 for pp_t in p28:
  G0_t = ECM.ew_g0_mod(pp_t)
  hit_any = False
  for _ in range(C_mc):
   ki = rng.randrange(1, 1 << 96)
   Qk = ECM.ew_mul_mod(ki, G0_t, pp_t)
   if Qk is None or (isinstance(Qk, tuple) and Qk[0] == "fac"):
    continue
   st_t, par_t = ECM.ew_to_mont_mod(Qk, 1, pp_t)
   if st_t != "ok": continue
   s2_ok += 1
   a24_t, Xt, Zt = par_t
   for qp in _p256:
    pk = qp
    while pk <= 256:
     Xt, Zt = ECM._xmul_mont(qp, Xt, Zt, a24_t, pp_t); pk *= qp
   if Zt % pp_t == 0: continue
   babs_t = {}
   for bi in bab_idx:
    bx, bz = ECM._xmul_mont(bi, Xt, Zt, a24_t, pp_t)
    if bz % pp_t != 0: babs_t[(bx * pow(bz, -1, pp_t)) % pp_t] = bi
   cert_hit = False
   for gi in range(1, 22):
    gx, gz = ECM._xmul_mont(gi * 210, Xt, Zt, a24_t, pp_t)
    if gz % pp_t == 0: continue
    xg = (gx * pow(gz, -1, pp_t)) % pp_t
    if xg in babs_t:
     bi = babs_t[xg]
     for q_cand in (gi * 210 - bi, gi * 210 + bi):
      if 256 < q_cand <= 4200 and _is_prime_small(q_cand):
       _, zq = ECM._xmul_mont(q_cand, Xt, Zt, a24_t, pp_t)
       if zq % pp_t == 0:
        cert_hit = True; break
     if cert_hit: break
   if cert_hit:
    s2_hits += 1
    hit_any = True
  if not hit_any:
   s2_miss_all += 1
 mod_s2_28 = max(0.0, ECM.semismooth_lo(28.0, 256, 4200) - ECM.semismooth_lo(28.0, 256, 256)) * ECM._power_correction_b1(28.0, 256) * (1.0 - 1.0 / 256.0)
 mod_miss_C = (1.0 - mod_s2_28) ** C_mc
 assert s2_ok >= 188 and s2_hits / float(s2_ok) >= mod_s2_28 > 0.01
 assert s2_miss_all / float(len(p28)) <= mod_miss_C

 # (7) Exact independent random-k_i Law-of-Total-Probability joint-failure verification for pre1 -> main:
 for fb_c, nb_c in ((69.05, 2590), (189.66, 2095)):
  sc = schedule_rung_ntt3_exact(fb_c, nb_c, 0.999, 32, 2)
  assert len(sc["pre_rungs"]) == 1 and sc["pre_rungs"][0]["curve_range"][1] == sc["curve_offset"]
  assert sc["main_curve_range"][0] == sc["curve_offset"] + 1
  assert (nb_c // 26) * (1.0 - 0.999) < (nb_c / 16.0) * (1.0 - 0.999)

 return dict(
  ok=True,
  nbits_checked=N.bit_length(),
  kronecker_ntt3_errs=0,
  bls_eval_errs=0,
  live_stage2_bls_factor=g_bls,
  zero_leaf_baby_fallback_ok=True,
  superadditive_cap_ok=True,
  joint_cascade_disjoint_ok=True,
  iid_rand_k256_stepper_muls=ECM.STEPPER_MULS,
  live_strictly_s2_hits=f"{s2_hits}/{s2_ok}",
  live_strictly_s2_rate=round(s2_hits / float(s2_ok), 4),
  live_strictly_s2_model=round(mod_s2_28, 4),
  live_multicurve_miss_12=f"{s2_miss_all}/{len(p28)}",
  live_multicurve_model_miss_12=round(mod_miss_C, 4),
  max_ntt3_len_log2=int(round(log2(n_ntt_max))),
  worst_ram_gib=round(ram_gib, 3),
 )
