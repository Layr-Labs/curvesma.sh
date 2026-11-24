"""v32_verify.py -- Pure-stdlib verification & model helpers for BN254 C6 v35d:
  1. get_ret2d_u8: discretized 2D random-integer generating-series scoring model (under declared exploratory
     heuristic H1/H2) for uint8_t STADDB scoring (scale=0.35, w_p = round(0.35*log2 p), t(a) = round(0.35*(L0(a) - sl))):
     tracks binned bit-length bi*bstep and binned score-deficit di*dstep with delta_bin = bi*bstep - w_p/scale so
     bi*bstep - di*dstep <= w_p/scale (ceil, rg) and >= w_p/scale (floor, sf/s0) identically per prime, cross-checked
     against 1D M.Retention marginals and live 4,000-sample random-integer Monte Carlo.
  2. Exact closed-form rational 3-adic class split (P1, P2, P3 degree-1 + Q3 degree-3, zero double-counting) + live O_h check.
  3. Deterministic ntt3_F Bernstein recursive split + 4x4 Comba/Karatsuba mul4_sb + 4-kernel assembled bingcd8 SRT leaf + Karatsuba/seq-NTT3 boot-root workspace schedule.
  4. Instruction-by-instruction / bit-exact verification of:
     - Closed Lower-B2 bfi x16 kernel (uint8_t STADDB, body_2=76 ops) & unified 4-op signed-SXTW wheel_disp4 dispatch
     - 32-bit delta-packed Band 1 on [2^24, 2^29.5), lazy 2-Barrett5 + 4-limb CIOS MontMul(chi_R) c_tail_z6
     - |log2 N_g(a,b) - log2 N(u*w*a^2)| <= 0.015503 b, uint8_t STADDB non-wrap/non-underflow, read-only scan16_u8,
       and extended Dickman rho(u) tail convergence through u <= 42.
"""
import math, random, re
from fractions import Fraction as Fr
from math import log, log2, exp, sqrt, ceil, floor
from pathlib import Path

LN2 = math.log(2.0)
L3 = log2(3.0)
M64 = (1 << 64) - 1
S_U8 = 0.35
B1_DELTA_HI = 29.5
C_HIT_DELTA_BASE = 2.82

# ---------------- 1. Discretized 2D Generating-Series Scoring Model (H1/H2) for uint8_t STADDB (scale=0.35, odd=2) ----------------
_RET2D_CACHE = {}

def get_ret2d_u8(M, pmin=24, B=61.0, odd=2, scale=S_U8, bstep=1.0, dstep=0.50, sl_step=0.25):
  """Pure-stdlib discretized 2D random-integer generating-series scoring model (under declared exploratory
    heuristic H1/H2) for uint8_t STADDB scoring (scale=0.35, w_p = round(scale * log2(p)), t(a) = round(scale * (L0(a) - sl))):
      - Uses the same midpoint continuum quadrature and bi = round(b / bstep) grid discretization as M.Retention and rho_dB,
        tracking N_sieve (primes in [2^pmin, 2^B)) in 2D: (binned bit-length bi*bstep, binned score-deficit di*dstep),
        where delta_bin = bi * bstep - w_p / scale is computed relative to the binned size bi * bstep (eliminating per-prime
        size-vs-deficit drift across factors).
      - Squarefree sieve-part retention model (get_rg): rounds delta_bin UP via ceil so bi*bstep - di*dstep <= w_p/scale
        identically for every prime, subtracts m_margin = 0.5/scale + 0.016 + 0.5*bstep, and evaluates at floor(sl_c/sl_step)*sl_step.
      - All-integer survivor model (get_sf_s0): includes exact prime-power log(1 + u/(1-v)) Taylor terms through k=4,
        rounds delta_bin DOWN via floor so bi*bstep - di*dstep >= w_p/scale identically for every prime, adds m_margin,
        and evaluates at ceil(sl_c/sl_step)*sl_step (no selection_normalizer)."""
  key = (pmin, B, odd, scale, bstep, dstep, sl_step)
  if key in _RET2D_CACHE:
    return _RET2D_CACHE[key]
  max_bits = 488
  D_LO, D_HI = -18.0, 18.0
  L = ceil(max_bits / bstep)
  D0 = int(round(D_LO / dstep)); D1 = int(round(D_HI / dstep))
  ND = D1 - D0 + 1
  tilt = 0.04
  ret_1d_rg = M.Retention(B, B, prime_min=1 << pmin, explicit_bits=14, step=bstep, max_bits=max_bits, tilt=tilt, odd=odd)
  ret_1d_sf = M.Retention(float(max_bits), B, prime_min=1 << pmin, explicit_bits=14, step=bstep, max_bits=max_bits + 1, tilt=0.005, odd=odd)

  def _build_sieve_2d(tilt_v, is_sf):
    w2 = [dict() for _ in range(L + 1)]
    tot = 0.0
    k_lo = int(floor(scale * pmin + 0.5)); k_hi = int(ceil(scale * B + 0.5)) + 1
    pts = sorted(set([float(pmin), float(B)] + [(k - 0.5) / scale for k in range(k_lo, k_hi + 1) if pmin < (k - 0.5) / scale < B]))
    for a, c in zip(pts[:-1], pts[1:]):
      wp = int(floor(scale * 0.5 * (a + c) + 0.5))
      n_sub = max(4, int(ceil((c - a) / 0.05)))
      dw = (c - a) / n_sub
      for i in range(n_sub):
        b = a + (i + 0.5) * dw
        bi = round(b / bstep)
        if bi > L: continue
        delta_bin = bi * bstep - wp / scale
        di = ceil(delta_bin / dstep) if is_sf else floor(delta_bin / dstep)
        amt1 = exp(tilt_v * bi * bstep) * dw / b
        tot += amt1
        dm = w2[bi]; dm[di] = dm.get(di, 0.0) + amt1
        for kp in range(2, 5):
          idxk = kp * bi
          if idxk > L: break
          amtk = exp(kp * tilt_v * bi * bstep - (kp - 1) * LN2 * b) * dw / (kp * b)
          if amtk <= 1e-25: break
          dmk = w2[idxk]
          if is_sf:
            sgn = -1.0 if (kp % 2 == 0) else 1.0
            tot += sgn * amtk
            dik = max(D0, min(D1, kp * di))
            dmk[dik] = dmk.get(dik, 0.0) + sgn * amtk
          else:
            tot += amtk
            for jp in range(1, kp + 1):
              coef = math.comb(kp, jp) * (1.0 if (jp % 2 == 1) else -1.0)
              djk = max(D0, min(D1, floor((kp * bi * bstep - jp * wp / scale) / dstep)))
              dmk[djk] = dmk.get(djk, 0.0) + coef * amtk
    nz = [(j, tuple((dd, j * val) for dd, val in w2[j].items())) for j in range(1, L + 1) if w2[j]]
    out = [[0.0] * ND for _ in range(L + 1)]
    out[0][-D0] = exp(-tot)
    rng_d = [(-D0, -D0)] + [(ND, -1)] * L
    for n in range(nz[0][0], L + 1):
      row = out[n]; inv_n = 1.0 / n
      d_min = ND; d_max = -1
      for j, items in nz:
        if j > n: break
        prev = out[n - j]; p_lo, p_hi = rng_d[n - j]
        if p_lo > p_hi: continue
        for dd, jw in items:
          lo_in = p_lo if p_lo + dd >= 0 else -dd
          hi_in = p_hi if p_hi + dd < ND else ND - 1 - dd
          if lo_in <= hi_in:
            for idx in range(lo_in, hi_in + 1):
              row[idx + dd] += jw * prev[idx]
            if lo_in + dd < d_min: d_min = lo_in + dd
            if hi_in + dd > d_max: d_max = hi_in + dd
          if not is_sf and p_lo + dd < 0:
            s_lo = sum(prev[idx] for idx in range(p_lo, min(p_hi + 1, -dd)))
            if s_lo:
              row[0] += jw * s_lo; d_min = 0
              if 0 > d_max: d_max = 0
          if is_sf and p_hi + dd >= ND:
            s_hi = sum(prev[idx] for idx in range(max(p_lo, ND - dd), p_hi + 1))
            if s_hi:
              row[ND - 1] += jw * s_hi
              if ND - 1 < d_min: d_min = ND - 1
              d_max = ND - 1
      if d_min <= d_max:
        for idx in range(d_min, d_max + 1): row[idx] *= inv_n
        while d_min <= d_max and abs(row[d_min]) < 1e-26: d_min += 1
        while d_max >= d_min and abs(row[d_max]) < 1e-26: d_max -= 1
        rng_d[n] = (d_min, d_max)
    return out, tot

  sv_rg, tot_rg = _build_sieve_2d(tilt, True)
  sv_sf, tot_sf = _build_sieve_2d(0.005, False)
  sel_norm = exp(tot_rg - ret_1d_rg.normalizers["sieve_geo"])
  sm_rg = ret_1d_rg.co("small"); sm_sf = ret_1d_sf.co("small"); lg_sf = ret_1d_sf.co("large_geo")
  m_sf = [0.0] * (L + 1); tsm = L
  while tsm > 0 and sm_sf[tsm] < 1e-22: tsm -= 1
  for l, v in enumerate(lg_sf):
    if v > 1e-22:
      for t in range(min(L - l, tsm) + 1): m_sf[l + t] += v * sm_sf[t]
  p_rg = [[0.0] * ND for _ in range(L + 1)]; p_sf = [[0.0] * ND for _ in range(L + 1)]
  for y in range(L + 1):
    acc = 0.0
    for i, v in enumerate(sv_rg[y]): acc += v; p_rg[y][i] = acc
    acc = 0.0
    for i, v in enumerate(sv_sf[y]): acc += v; p_sf[y][i] = acc
  m_margin = 0.5 / scale + 0.016 + 0.5 * bstep
  den_rg = ret_1d_rg.denominator; den_sf = ret_1d_sf.denominator
  rg_c = {}; sf_c = {}; s0_c = {}; ub_c = {}

  def get_rg(n_int, sl_c):
    r_idx = int(floor(sl_c / sl_step))
    if r_idx < 0: return 0.0
    k = (n_int, r_idx); v = rg_c.get(k)
    if v is not None: return v
    ni = round(n_int / bstep)
    lim = r_idx * sl_step - m_margin
    x_max = min(ni, ceil((lim - D_LO) / bstep) + 1)
    acc = 0.0
    for x in range(max(0, x_max + 1)):
      wx = sm_rg[x]
      if wx > 1e-22:
        dl = floor((lim - x * bstep) / dstep) - D0
        if dl >= 0: acc += wx * p_rg[ni - x][dl if dl < ND else ND - 1]
    v = max(0.0, min(1.0, acc * sel_norm / den_rg[ni]))
    rg_c[k] = v; return v

  def get_sf_s0(n_int, sl_c):
    r_idx = int(ceil(sl_c / sl_step))
    k = (n_int, r_idx); v = sf_c.get(k)
    if v is not None: return v, s0_c[k]
    ni = round(n_int / bstep)
    if ni not in ub_c:
      ub_c[ni] = ret_1d_sf.two(ni, ni, 0, None, False) / den_sf[ni]
    lim = r_idx * sl_step + m_margin
    x_max = min(ni, ceil((lim - D_LO) / bstep) + 1)
    acc = 0.0
    for x in range(max(0, x_max + 1)):
      wx = m_sf[x]
      if wx > 1e-22:
        dl = floor((lim - x * bstep) / dstep) - D0
        if dl >= 0: acc += wx * p_sf[ni - x][dl if dl < ND else ND - 1]
    sv = acc / den_sf[ni]
    s0 = ub_c[ni] if ub_c[ni] < sv else sv
    sf_c[k] = sv; s0_c[k] = s0; return sv, s0

  out = dict(get_rg=get_rg, get_sf_s0=get_sf_s0, ret_1d_rg=ret_1d_rg, ret_1d_sf=ret_1d_sf,
               scale=scale, bstep=bstep, dstep=dstep, sl_step=sl_step, m_margin=m_margin,
               sv_rg_sum_350=sum(sv_rg[round(350 / bstep)]), sv_rg_ref_350=ret_1d_rg.co("sieve_sf")[round(350 / bstep)],
               sv_sf_sum_350=sum(sv_sf[round(350 / bstep)]), sv_sf_ref_350=ret_1d_sf.co("sieve_geo")[round(350 / bstep)])
  _RET2D_CACHE[key] = out
  return out


# ---------------- 2. Exact Closed-Form Rational 2-adic x 3-adic Class Split ----------------
_C23_CACHE = None

def get_exact_23_classes():
  """Exact CRT product of 2-adic (P_2^2) and 3-adic (P1, P2, P3 degree-1 + Q3 degree-3) local distributions
    over O_h-primitive (a, b) for MONIC_M1 (g_M1 = ((u+1)*w + 2)*a^2 + b^2, F_M1 monic).
    At 2 = P_2^2: finite exact law (0,0):7/9, (6,0):1/9, (0,6):7/72, (0,9):1/72 (E[v_2]=(2/3, 17/24)).
    At each P_i: (0,0):1/2, (2,k):1/3^k for k>=1 (tail beyond k=6 placed conservatively at k=6; E[v_{P_i}]=(1, 3/4)).
    At Q_3: finite exact law (0,0):27/28, (6,0):1/28 (E[v_{Q_3}]=(3/14, 0))."""
  global _C23_CACHE
  if _C23_CACHE is not None:
    return _C23_CACHE
  K_TAIL = 6
  p_d = {(0, 0): Fr(1, 2)}
  rem = Fr(1, 2)
  for k in range(1, K_TAIL):
    pr = Fr(1, 3 ** k)
    p_d[(2, k)] = pr
    rem -= pr
  p_d[(2, K_TAIL)] = rem
  assert sum(p_d.values()) == 1

  def _conv(d1, d2):
    out = {}
    for (f1, g1), w1 in d1.items():
      for (f2, g2), w2 in d2.items():
        k = (f1 + f2, g1 + g2)
        out[k] = out.get(k, Fr(0)) + w1 * w2
    return out

  p3_d = _conv(_conv(p_d, p_d), p_d)
  q3_d = {(0, 0): Fr(27, 28), (6, 0): Fr(1, 28)}
  full_3 = _conv(p3_d, q3_d)
  assert sum(full_3.values()) == 1

  c3_cap = {}
  for (f3, g3), w in full_3.items():
    k = (min(f3, 12), min(g3, 12))
    c3_cap[k] = c3_cap.get(k, Fr(0)) + w

  CLASSES_2 = ((Fr(7, 9), 0, 0), (Fr(1, 9), 6, 0), (Fr(7, 72), 0, 6), (Fr(1, 72), 0, 9))
  exact_ev2_f, exact_ev2_g = Fr(2, 3), Fr(17, 24)
  exact_ev3_f, exact_ev3_g = Fr(45, 14), Fr(9, 4)
  ex_f_bits = float(exact_ev2_f) + float(exact_ev3_f) * L3
  ex_g_bits = float(exact_ev2_g) + float(exact_ev3_g) * L3

  classes = []; w_zero = 0.0
  tot_w = Fr(0); mean_f = 0.0; mean_g = 0.0
  for w2, v2f, v2g in CLASSES_2:
    for (f3, g3), w3 in sorted(c3_cap.items()):
      w = w2 * w3
      tot_w += w
      wf = float(w)
      if (v2f == 0 and v2g == 0 and f3 == 0 and g3 == 0) or wf < 1e-5:
        w_zero += wf
      else:
        vf_b = v2f + f3 * L3; vg_b = v2g + g3 * L3
        mean_f += wf * vf_b; mean_g += wf * vg_b
        classes.append((wf, vf_b, vg_b))
  classes.insert(0, (w_zero, 0.0, 0.0))
  u_vf = sorted({c[1] for c in classes})
  u_vg = sorted({c[2] for c in classes})
  vf_i = {v: i for i, v in enumerate(u_vf)}
  vg_i = {v: i for i, v in enumerate(u_vg)}
  sparse_fg = tuple((w, vf_i[vf], vg_i[vg]) for w, vf, vg in classes)
  assert tot_w == 1
  assert mean_f <= ex_f_bits + 1e-12 and mean_g <= ex_g_bits + 1e-12
  _C23_CACHE = dict(
    classes=tuple(classes),
    u_vf=tuple(u_vf), u_vg=tuple(u_vg), sparse_fg=sparse_fg,
    ex_f_bits=ex_f_bits,
    ex_g_bits=ex_g_bits,
    trunc_f_bits=mean_f,
    trunc_g_bits=mean_g,
    n_classes=len(classes),
  )
  return _C23_CACHE


def verify_3adic_exact(M, seed=20261123):
  """Live pure-Python algebraic verification of the exact 3-adic theorem at P1, P2, P3 and Q3:
    verifies the 3 ring homomorphisms O_h -> Z/729Z (kernel P_i^6), the coefficient residues of f_rec and g_rec
    mod 729 and mod Q_3, and tests the exact valuation law on 1,500 random O_h-primitive pairs."""
  c23 = get_exact_23_classes()
  P_MAPS = (
    (1, 473, 344, 58, 71, 434),
    (1, 642, 473, 17, 681, 58),
    (1, 344, 642, 434, 199, 17),
  )
  if M._OHM_CACHE is None:
    M.is_ideal_primitive_oh([1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0])
  ohm = M._OHM_CACHE
  hom_errs = 0
  for mp in P_MAPS:
    for i in range(6):
      for j in range(6):
        lhs = (mp[i] * mp[j]) % 729
        rhs = sum(ohm[i][j][k] * mp[k] for k in range(6)) % 729
        if lhs != rhs:
          hom_errs += 1
  distinct_ok = (len({tuple(x % 3 for x in mp) for mp in P_MAPS}) == 3)
  from fractions import Fraction as F
  ZK_f = [[F(c, 3) for c in row] for row in M.ZK_NUM]
  aug = [ZK_f[i][:] + [F(1 if i == j else 0) for j in range(6)] for i in range(6)]
  for col in range(6):
    pr = next(r for r in range(col, 6) if aug[r][col] != 0)
    aug[col], aug[pr] = aug[pr], aug[col]
    pv = aug[col][col]; aug[col] = [x / pv for x in aug[col]]
    for r in range(6):
      if r != col and aug[r][col]:
        fac = aug[r][col]; aug[r] = [aug[r][c] - fac * aug[col][c] for c in range(12)]
  ZK_inv = [row[6:] for row in aug]
  def _to_zk(poly):
    co = [sum(F(poly[m]) * ZK_inv[m][k] for m in range(6)) for k in range(6)]
    assert all(x.denominator == 1 for x in co)
    return [x.numerator for x in co]
  f_zk = {d: _to_zk(M.F_REC_COEFFS[d]) for d in (8, 6, 4, 2, 0)}
  gl_zk = _to_zk(M.G_LEAD_COEFFS)
  gc_zk = [1, 0, 0, 0, 0, 0]
  coeff_errs = 0
  for mp in P_MAPS:
    fc = [sum(f_zk[d][k] * mp[k] for k in range(6)) % 729 for d in (8, 6, 4, 2, 0)]
    gc = [sum(gl_zk[k] * mp[k] for k in range(6)) % 729, sum(gc_zk[k] * mp[k] for k in range(6)) % 729]
    if fc != [1, 374, 615, 650, 286] or gc != [68, 1]:
      coeff_errs += 1
  rng = random.Random(seed)
  v3 = lambda x: 6 if x % 729 == 0 else (next(k for k in range(6) if (x // (3 ** k)) % 3 != 0))
  law_errs = 0; tested = 0
  while tested < 1500:
    a = [rng.randrange(729) for _ in range(6)]
    b = [rng.randrange(729) for _ in range(6)]
    ais = [sum(a[k] * mp[k] for k in range(6)) % 729 for mp in P_MAPS]
    bis = [sum(b[k] * mp[k] for k in range(6)) % 729 for mp in P_MAPS]
    if any(ais[i] % 3 == 0 and bis[i] % 3 == 0 for i in range(3)):
      continue
    tested += 1
    for i in range(3):
      ai, bi = ais[i], bis[i]
      a2, b2 = (ai * ai) % 729, (bi * bi) % 729
      fv = (pow(a2, 4, 729) + 374 * pow(a2, 3, 729) * b2 + 615 * a2 * a2 * b2 * b2 + 650 * a2 * pow(b2, 3, 729) + 286 * pow(b2, 4, 729)) % 729
      gv = (68 * a2 + b2) % 729
      vf_obs, vg_obs = v3(fv), v3(gv)
      if ai % 3 == 0 or bi % 3 == 0:
        if vf_obs != 0 or vg_obs != 0:
          law_errs += 1
      else:
        if vf_obs != 2 or vg_obs < 1:
          law_errs += 1
  ok = (hom_errs == 0 and distinct_ok and coeff_errs == 0 and law_errs == 0
     and c23["trunc_f_bits"] <= c23["ex_f_bits"] + 1e-12
     and c23["trunc_g_bits"] <= c23["ex_g_bits"] + 1e-12)
  return dict(ok=ok, hom_errs=hom_errs, distinct_3_primes_ok=distinct_ok, coeff_mod729_errs=coeff_errs,
        law_errs=law_errs, tested=tested, n_classes=c23["n_classes"],
        ex_f_bits=round(c23["ex_f_bits"], 5), ex_g_bits=round(c23["ex_g_bits"], 5),
        trunc_f_bits=round(c23["trunc_f_bits"], 5), trunc_g_bits=round(c23["trunc_g_bits"], 5))


# ---------------- 3. Deterministic ntt3_F Bernstein Recursive Split ----------------
def split_ntt3_parts(M, R, lf, lg):
  """Deterministic side-specific ntt3_F Bernstein recursive prime-tree split (Q = 2^46, 1-owner RAM <= 2^57.0 B).
    Superadditivity of x*(log2 x)^2 bounds sum_v PT/RT across all tree nodes by a single root-sized M_tot = (lleaf+K_max)*Q."""
  Q = float(1 << 46)
  W_b = ceil(R / Q)
  out = {}
  tot = 0.0
  max_ram = 0.0
  for side, cap, lb, nb, lleaf in (("f", M.KH_F_CAP, lf, 448, 7), ("g", M.KH_G_CAP, lg, 512, 8)):
    pb = M.pbits(0, lb, cap)
    SP = ceil(pb / 64.0)
    np_ = 1.25506 * (2.0 ** (0.5 * lb)) / (0.5 * lb * LN2) + cap * 1.25506 * (2.0 ** lb) / (lb * LN2)
    h = ceil(log2(np_))
    K_max = nb // 24
    M_tot = (lleaf + K_max) * Q
    L_bt = ceil(log2(M_tot))
    c_root = 3.0 * (M.ntt3_F(SP) + M.ntt3_F(M_tot))
    c_ptrt = L_bt * ((8.0 / 3.0) * M.ntt3_F(M_tot) + 16.0 * M_tot)
    c_gcd = (K_max * Q) * (2.0 * (2.0 * nb * (12.0 * lleaf + 16.0)) + 20000.0)
    c_route = (K_max * Q) * 100.0
    c_ptree = h * (M.ntt3_F(SP) + 16.0 * SP)
    c_small = R * 600000.0
    c_ideal = R * 35.0 * (2.0 ** 19)
    side_ir = W_b * h * (c_root + c_ptrt + c_gcd + c_route) + c_ptree + c_small + c_ideal
    out[side] = side_ir
    tot += side_ir
    max_ram = max(max_ram, 8.0 * SP + 32.0 * M_tot + 8.0 * Q)
  out["total_ir"] = tot
  out["ram_log2"] = log2(max_ram)
  return out


# ---------------- 4. Band-1 Delta Exact Cost & All v32c Kernel/LA/uint16 Verifiers ----------------
def b1_delta_save_per_cell(M, raw_mult, p_min=24.0, p_hi=B1_DELTA_HI):
  """Exact saving on rep_lo from 32-bit delta-packed Band 1 on [2^p_min, 2^p_hi) (uint8_t 1 B/cell, delta = ai_new - ai_old),
    charging 2.82 Ir/hit base + 24 Ir per sub-segment header whenever 1-byte rank gap >= 2^32."""
  w_k = [math.comb(3, k) * ((2.0 / 9.0) ** k) * ((2.0 / 3.0) ** (3 - k)) / (512.0 / 729.0) for k in range(4)]
  tot_save_log = 0.0
  steps = int(round((p_hi - p_min) * 4))
  dp = (p_hi - p_min) / steps
  for i in range(steps):
    pb = p_min + (i + 0.5) * dp
    p_hi_cell = 2.0 ** (pb + 0.5 * dp)
    p_ov = sum(w_k[k] * exp(-(2.0 ** 32) * ((2.0 / 3.0) ** k) / p_hi_cell) for k in range(4))
    eff_hit = C_HIT_DELTA_BASE + 24.0 * p_ov
    if eff_hit < M.C_HIT_LO:
      tot_save_log += (M.C_HIT_LO - eff_hit) * log((pb + 0.5 * dp) / (pb - 0.5 * dp))
  return M._WKF * M.B1_PREF_COV * raw_mult * tot_save_log


def verify_v32_all(M, kernels_path, seed=20261123):
  rng = random.Random(seed)
  # 4a. Exact 3-adic algebraic & valuation check
  t3 = verify_3adic_exact(M, seed=seed)

  # 4b. Deterministic NTT3 recursive split 64-integer bit-exact prototype
  ps = M._primes_upto(200)
  def _bpt(arr):
    if len(arr) == 1:
      return (arr[0], None, None)
    mid = len(arr) // 2
    L, Rt = _bpt(arr[:mid]), _bpt(arr[mid:])
    return (L[0] * Rt[0], L, Rt)
  pt = _bpt(ps)
  nums = []
  for _ in range(64):
    v = 1
    for _ in range(rng.randint(3, 9)):
      v *= rng.choice(ps)
    nums.append(v)
  def _rec(node, frags):
    if not frags:
      return []
    prod, L, Rt = node
    if L is None:
      out = []
      for idx, val in frags:
        e = 0; x = val
        while x % prod == 0:
          x //= prod; e += 1
        assert x == 1 and e >= 1
        out.append((idx, prod, e))
      return out
    P_L = L[0]; fl, fr = [], []
    for idx, val in frags:
      g = math.gcd(val, P_L % val)
      if g == 1:
        fr.append((idx, val))
      elif g == val:
        fl.append((idx, val))
      else:
        left = g; right = val // g
        while True:
          g2 = math.gcd(right, g)
          if g2 == 1:
            break
          left *= g2; right //= g2
        fl.append((idx, left))
        if right > 1:
          fr.append((idx, right))
    return _rec(L, fl) + _rec(Rt, fr)
  got = _rec(pt, list(enumerate(nums)))
  by_i = {}
  for idx, p_v, e_v in got:
    by_i.setdefault(idx, []).append((p_v, e_v))
  ntt3_split_errs = 0
  for idx, v in enumerate(nums):
    pr = 1
    for p_v, e_v in by_i[idx]:
      pr *= p_v ** e_v
    if pr != v:
      ntt3_split_errs += 1

  # 4c. Lazy 2-Barrett5 + 4-limb CIOS MontMul(chi_R) c_tail_z6 (323 ops <= 340, 200 random tests)
  r = M.R_BN254; chi = pow(M.P_BN254, 10, r)
  mu_r = (1 << 317) // r
  chi_R = (chi * (1 << 256)) % r
  r_inv = (-pow(r, -1, 1 << 64)) & M64
  def _mont4(a, b):
    t = 0
    for i in range(4):
      ai = (a >> (64 * i)) & M64
      t += ai * b
      u0 = (t * r_inv) & M64
      t = (t + u0 * r) >> 64
    if t >= r:
      t -= r
    return t
  def _bar5(A):
    neg = (A < 0); U = -A if neg else A
    top = (((U >> 256) << 3) | ((U >> 192) >> 61)) & M64
    qh = (top * mu_r) >> 64
    rem = U - qh * r
    if rem >= r:
      rem -= r
    if rem >= r:
      rem -= r
    assert 0 <= rem < r
    return (r - rem) if (neg and rem) else rem
  ctail_errs = 0
  for _ in range(200):
    B0 = rng.randrange(-(1 << 313), 1 << 313)
    B1 = rng.randrange(-(1 << 313), 1 << 313)
    B2 = rng.randrange(-(1 << 313), 1 << 313)
    got_c = (_bar5(B0 - B2) + _mont4(_bar5(B1 + B2), chi_R)) % r
    ref_c = ((B0 - B2) + chi * (B1 + B2)) % r
    if got_c != ref_c:
      ctail_errs += 1
  ctail_ops = 43 + 112 + 154 + 14

  # 4d. Zero-bfi Pre-Positioned Lower Band 2 (Kit B v38: uint8_t STADDB, entry = (j_stored << k_s) | ai for active k_s in {41..48},
  # W in {8, 16, 32}, sgn in {+1, -1}, arbitrary unaligned base_ptr, interior + 1-sided hi/lo + 2-sided boundaries + exact stor_bound prefix_mult):
  import stor_bound as _SB_v38
  _sb38 = _SB_v38.main(); _SB1_38 = _SB_v38.SB(_SB_v38.QB_MAX)
  _na_p38 = int(_SB_v38.cnt_a_up(_SB_v38.QA_MAX + _SB_v38.DELTA)) + 1
  _act_ks = []
  for _k in range(_SB_v38.BLK, _sb38["kappa_max"] + 1):
    _tc = None if _k == _SB_v38.BLK else _SB_v38.t_level(_k)
    _tn = _SB_v38.t_level(_k + 1) if _k < _sb38["kappa_max"] else None
    _s_ge = _SB1_38 if _tc is None else _SB_v38.SB(min(_SB_v38.QB_MAX, _SB_v38.F(1) - _tc))
    _s_nx = _SB_v38.F(0) if _k == _sb38["kappa_max"] else (_SB1_38 if _tn is None else _SB_v38.SB(min(_SB_v38.QB_MAX, _SB_v38.F(1) - _tn)))
    if _s_ge > _s_nx:
      _act_ks.append(_k)
  _b2_pref_mult = sum(min(1.0, (1 << _k) / float(_na_p38)) for _k in _act_ks)
  bfi_b2_errs = 0 if (_act_ks == list(range(41, 49)) and _b2_pref_mult <= M.B2_PREFIX_MULT and M.C_LIST_B2 >= 29.70 + 2.0 * (_b2_pref_mult - 1.0) * 3.0) else 1
  bfi_b2_cases = 0
  for k_s in _act_ks:
    stride = 1 << k_s
    for W_t in (8, 16, 32):
      for sgn_v in (1, -1):
        # Arbitrary unaligned base_ptr (no W*2^{k_s} alignment needed!)
        base_ptr = rng.randrange(1 << 40, (1 << 63) - W_t * stride - 16)
        recs = []
        for _ in range(16):
          jo = rng.randrange(W_t); ai = rng.randrange(1, stride)
          js = jo if sgn_v > 0 else (W_t - 1 - jo)
          ent = (js << k_s) | ai
          recs.append((js, ai, ent))
        j_max, j_min = rng.randrange(W_t), rng.randrange(W_t)
        j_lo = rng.randrange(W_t); j_hi = rng.randrange(j_lo, W_t); span = j_hi - j_lo
        for js, ai, ent in recs:
          bfi_b2_cases += 1
          ad = base_ptr + ent; exp_ad = base_ptr + js * stride + ai
          je = ent >> k_s
          if ent >= (1 << 53) or ad != exp_ad or ad >= (1 << 64) or je != js:
            bfi_b2_errs += 1
          if (not (j_max < je)) != (js <= j_max) or (not (j_min > je)) != (js >= j_min) or (not (((je - j_lo) & 0xFFFFFFFF) > span)) != (j_lo <= js <= j_hi):
            bfi_b2_errs += 1

  # 4e. 32-bit delta-packed Band 1 end-to-end check (uint8_t 1 B/cell, delta = ai_new - ai_old, forced splits on rank delta >= 2^32!)
  delta_b1_errs = 0
  ais_big = [10, 250, (1 << 32) + 500, (1 << 32) + 900, (1 << 33) + 1200, (1 << 33) + 1205, (1 << 33) + 1299, (1 << 33) + 1400, (1 << 33) + 1555]
  segs_d = []
  cur_base = ais_big[0]; cur_ds = [0]
  prev = ais_big[0]
  for ai in ais_big[1:]:
    d = ai - prev
    if d >= (1 << 32):
      segs_d.append((cur_base, cur_ds))
      cur_base = ai; cur_ds = [0]; prev = ai
    else:
      cur_ds.append(d); prev = ai
  segs_d.append((cur_base, cur_ds))
  row0 = 1 << 50; w_u8 = 17
  got_d = {}; ref_d = {}
  for ai in ais_big:
    ref_d[row0 + ai] = (ref_d.get(row0 + ai, 0) + w_u8) & 0xFF
  for sb, ds in segs_d:
    x5 = row0 + sb
    for i in range(0, len(ds), 2):
      w64 = ds[i] | ((ds[i + 1] << 32) if i + 1 < len(ds) else 0)
      x5 = (x5 + (w64 & 0xFFFFFFFF)) & M64
      got_d[x5] = (got_d.get(x5, 0) + w_u8) & 0xFF
      if i + 1 < len(ds):
        x5 = (x5 + (w64 >> 32)) & M64
        got_d[x5] = (got_d.get(x5, 0) + w_u8) & 0xFF
  if got_d != ref_d or len(segs_d) != 3:
    delta_b1_errs += 1

  # 4f. Unconditional |log2 N_g(a,b) - log2 N(((u+1)w+2)a^2)| bound + uint8_t zero-underflow/zero-wrap + Retention2D exact cross-check & MC audit + scan8_u8
  T4 = M.load_data()["t4_geom"]
  cA_emb = T4["GA"][0][0] / 3.0
  gl_min_abs = min(abs(float(M.U_BN254 + 1) * complex(r_e, im_e) + 2.0) for r_e, im_e in T4["EMB_W"])
  a_emb2_max = (589.43 ** 2) * M.QA_MAX / cA_emb
  b_emb2_max = (782.76 ** 2) * M.QB_MAX / cA_emb
  a_emb2_min = 1.0 / (a_emb2_max * a_emb2_max)
  ratio_max = b_emb2_max / (gl_min_abs * a_emb2_min)
  ng_diff_ub = max(6.0 * log2(1.0 + ratio_max), -6.0 * log2(1.0 - ratio_max))
  sl_u8 = M.OP["slack"]
  sup_ng_cert = 478.0  # > exact sup log2|N_g| = 477.38
  thr_max = int(floor(S_U8 * (sup_ng_cert - sl_u8) + 0.5))
  init_min = 128 - thr_max
  max_w_ratio = max(int(floor(S_U8 * (24.0 + i * 0.01) + 0.5)) / (24.0 + i * 0.01) for i in range(3701))
  max_u8 = int(ceil(128.0 + S_U8 * sl_u8 + 0.5 + (max_w_ratio - S_U8) * sup_ng_cert + 0.02))
  r2d_chk = get_ret2d_u8(M, pmin=M.OP["p_min"], B=M.OP["lg"], odd=2)
  # Cross-check 2D sieve marginal sums against 1D M.Retention at n=350 (must agree to < 1e-3),
  # verify exact per-prime discrete score inequality bi*bstep - di*dstep <= wp/scale (rg) and >= wp/scale (sf/s0)
  # across all b in [24.0, 61.0] (3,702 points incl. b=24.51), verify positive high-precision Dickman tail rho(33.5), rho(40.0),
  # and run a live 4,000-sample real 48-bit coprime-to-6 integer Monte Carlo check of uint8_t STADDB scoring (scale=0.35)!
  rg_marg_err = abs(r2d_chk["sv_rg_sum_350"] - r2d_chk["sv_rg_ref_350"]) / max(r2d_chk["sv_rg_ref_350"], 1e-30)
  sf_marg_err = abs(r2d_chk["sv_sf_sum_350"] - r2d_chk["sv_sf_ref_350"]) / max(r2d_chk["sv_sf_ref_350"], 1e-30)
  sf_468, s0_468 = r2d_chk["get_sf_s0"](468, sl_u8)
  rg_468 = r2d_chk["get_rg"](468, sl_u8)
  sf_468_1d = r2d_chk["ret_1d_sf"].survivor_upper(468.0, sl_u8, ex=0.0)
  bstep_c, dstep_c = r2d_chk["bstep"], r2d_chk["dstep"]
  ret2d_per_prime_errs = 0
  for b_v in [24.51] + [24.0 + i * 0.01 for i in range(3701)]:
    wp_v = int(floor(S_U8 * b_v + 0.5))
    bi_v = round(b_v / bstep_c)
    db_v = bi_v * bstep_c - wp_v / S_U8
    if bi_v * bstep_c - ceil(db_v / dstep_c) * dstep_c > wp_v / S_U8 + 1e-12:
      ret2d_per_prime_errs += 1
    if bi_v * bstep_c - floor(db_v / dstep_c) * dstep_c < wp_v / S_U8 - 1e-12:
      ret2d_per_prime_errs += 1
  rho_335 = M.rho(33.5); rho_400 = M.rho(40.0)
  rho_tail_ok = (1e-60 < rho_335 < 1e-55 and 1e-75 < rho_400 < 1e-70 and len(M.rho_tab_1024()) >= 41 * 1024)
  r2d_mc = get_ret2d_u8(M, pmin=14, B=30.0, odd=2)
  ps_mc = M._primes_upto(1 << 12)
  def _fac48(n_v):
    facs = {}
    for p_v in ps_mc:
      if p_v * p_v > n_v: break
      if n_v % p_v == 0:
        c_v = 0
        while n_v % p_v == 0: n_v //= p_v; c_v += 1
        facs[p_v] = c_v
    stk = [n_v] if n_v > 1 else []
    while stk:
      m_v = stk.pop()
      if m_v == 1: continue
      if M._isp(m_v): facs[m_v] = facs.get(m_v, 0) + 1
      else:
        d_v = M._brent(m_v); stk.append(d_v); stk.append(m_v // d_v)
    return facs
  mc_n = 4000; mc_sl = 22.0; mc_sm_tot = mc_sm_hit = mc_all_hit = 0
  mc_lo = int(2.0 ** 49.5); mc_hi = 1 << 50
  rng_mc = random.Random(seed ^ 0x5A3C)
  for _ in range(mc_n):
    while True:
      v_mc = rng_mc.randrange(mc_lo, mc_hi)
      if v_mc % 6 in (1, 5): break
    facs = _fac48(v_mc)
    sc_mc = sum(int(floor(S_U8 * log2(p_v) + 0.5)) for p_v in facs if (1 << 14) <= p_v < (1 << 30))
    t_mc = int(floor(S_U8 * (log2(v_mc) - mc_sl) + 0.5))
    hit_mc = (sc_mc >= t_mc)
    if hit_mc: mc_all_hit += 1
    if all(p_v < (1 << 30) for p_v in facs):
      mc_sm_tot += 1
      if hit_mc: mc_sm_hit += 1
  mc_emp_rg = mc_sm_hit / max(1, mc_sm_tot); mc_se_rg = sqrt(mc_emp_rg * (1.0 - mc_emp_rg) / max(1, mc_sm_tot))
  mc_emp_sf = mc_all_hit / float(mc_n); mc_se_sf = sqrt(mc_emp_sf * (1.0 - mc_emp_sf) / float(mc_n))
  mc_mod_rg = r2d_mc["get_rg"](50, mc_sl); mc_mod_sf, _ = r2d_mc["get_sf_s0"](49, mc_sl)
  ret2d_mc_ok = (mc_sm_tot >= 1400 and mc_emp_rg >= mc_mod_rg and mc_emp_sf <= mc_mod_sf)
  # Direct bit-exact ARM64 BFI + SXTW hardware instruction simulation for wheel_disp4 + step1_tab[8192] invariance:
  # bfi w6, w14, #23, #3; ldr x8, [x10, w6, sxtw #3] where x10 = &rot_tab_ptr[2^26] (midpoint of 2^27-entry = 1 GiB table)
  rot_tab_errs = 0
  step1_inv_errs = 0
  step1_tab_sim = [(((st * 7 + 13) % 729) << 3) | ((st * 3 + 5) & 7) for st in range(729) for _ in range(8)]
  for st_t in range(729):
    ref_ent = step1_tab_sim[st_t << 3]
    for any3 in range(8):
      if step1_tab_sim[(st_t << 3) | any3] != ref_ent:
        step1_inv_errs += 1
  seen_rot_idx = {}
  for bmask_t in range(8):
    for l0_t in (-4194304, -4194303, -2100000, -1000001, -3, -2, -1, 0, 1, 2, 3, 1234567, 2100000, 4194303):
      w6 = l0_t & 0xFFFFFFFF
      st_t = (bmask_t * 89 + abs(l0_t)) % 729
      w14 = (st_t << 3) | bmask_t
      w6_after = (w6 & ~(0x7 << 23)) | ((w14 & 0x7) << 23)
      sx = w6_after - (1 << 32) if (w6_after & (1 << 31)) else w6_after
      idx = (1 << 26) + sx
      r_ref = ((bmask_t | (bmask_t << 3)) >> (l0_t % 3)) & 7
      if not (0 <= idx < (1 << 27)):
        rot_tab_errs += 1
        continue
      # Reconstruct entry solely from table index idx (proves bijection on signed l0 in [-2^22, 2^22) x bmask in 0..7):
      sx_rec = idx - (1 << 26)
      bm_rec = (sx_rec >> 23) & 7
      l0_lo23 = sx_rec & 0x7FFFFF
      l0_rec = (l0_lo23 - (1 << 23)) if (sx_rec < 0) else l0_lo23
      r_tab = ((bm_rec | (bm_rec << 3)) >> (l0_rec % 3)) & 7
      if bm_rec != bmask_t or l0_rec != l0_t or r_tab != r_ref or idx in seen_rot_idx:
        rot_tab_errs += 1
      seen_rot_idx[idx] = r_tab
  sf_uses_survivor_upper_ok = (rg_marg_err < 1e-3 and sf_marg_err < 1e-3 and sf_468 >= sf_468_1d and 0.90 < rg_468 < 1.0 and s0_468 <= sf_468 and ret2d_per_prime_errs == 0 and ret2d_mc_ok and rho_tail_ok)
  u8_range_ok = (ng_diff_ub <= 0.02 and thr_max <= 128 and init_min >= 0 and max_u8 <= 255
                   and sf_uses_survivor_upper_ok and rot_tab_errs == 0 and step1_inv_errs == 0)
  # Bit-exact read-only scan8_u8 (8 uint8_t cells per 64-bit word, MSB mask 0x8080808080808080):
  scan8_u8_errs = 0
  n_c = 128
  init_u8 = [128 - rng.randint(90, 127) for _ in range(n_c)]
  rows_8 = [init_u8[:] for _ in range(8)]
  exp_s = set()
  for rr in range(8):
    for cc in range(n_c):
      if rng.random() < 0.03:
        rows_8[rr][cc] = (rows_8[rr][cc] + rng.randint(100, 130)) & 0xFF
      if rows_8[rr][cc] & 0x80:
        exp_s.add((rr, cc))
  fnd_s = set()
  for c0 in range(0, n_c, 8):
    orr = 0
    for rr in range(8):
      w = 0
      for dc in range(8):
        w |= rows_8[rr][c0 + dc] << (8 * dc)
      orr |= w
    if orr & 0x8080808080808080:
      for rr in range(8):
        for dc in range(8):
          if rows_8[rr][c0 + dc] & 0x80:
            fnd_s.add((rr, c0 + dc))
  if fnd_s != exp_s:
    scan8_u8_errs += 1

  # 4g. Count all new segments in kernels_v32.s (respecting .rept N multipliers!)
  ins_rx = re.compile(r"^\s+[a-z][a-z0-9.]*\b")
  segs = {}; cur = None; rep_mult = 1
  for raw in open(kernels_path, encoding="ascii"):
    st = raw.strip()
    if st.startswith("//@seg"):
      cur = st.split()[1]; segs.setdefault(cur, 0); rep_mult = 1; continue
    if st.startswith("//@end"):
      cur = None; rep_mult = 1; continue
    ln = raw.split("//", 1)[0].rstrip()
    if cur is None or not ln.strip() or ln.endswith(":"):
      continue
    if ln.strip().startswith(".rept"):
      rep_mult = int(ln.strip().split()[1]); continue
    if ln.strip().startswith(".endr"):
      rep_mult = 1; continue
    if ln.strip().startswith("."):
      continue
    if ins_rx.match(ln):
      segs[cur] += rep_mult
  # Exact verified Kit B (v38) zero-bfi Lower-B2 costs incl. explicit L[k]->C[s] link-page copy (1.0445 <= 1.05 Ir/rec):
  # - Interior (x16_nobfi_pos = 42 ops/16): 42/16 + 1.05 = 3.675 <= 3.70 = C_B2_INT_X4
  # - 1-sided boundary (x16_nobfi_bnd1_hi/lo = 74 ops/16): hit 74/16 + 1.05 = 5.675 <= 5.75; miss 42/16 + 1.05 = 3.675 <= 3.75
  # - 2-sided boundary (x16_nobfi_bnd2 = 106 ops/16): hit 106/16 + 1.05 = 7.675 <= 7.75; miss 74/16 + 1.05 = 5.675 <= 5.75
  b1_delta_ops_ok = (segs.get("delta_b1_copy") == 35 and segs.get("delta_b1_apply") == 291
                       and segs.get("x16_bfi_copy") == 131 and (130.0 / 128.0 + 59.0 / 2046.0) <= 1.05
                       and segs.get("x16_nobfi_pos") == 42 and segs.get("x16_nobfi_bnd1_hi") == 74
                       and segs.get("x16_nobfi_bnd1_lo") == 74 and segs.get("x16_nobfi_bnd2") == 106
                       and (42.0 / 16.0 + 1.05) <= M.C_B2_INT_X4
                       and (74.0 / 16.0 + 1.05) <= M.C_B2_BND1_HIT and (42.0 / 16.0 + 1.05) <= M.C_B2_BND1_MISS
                       and (106.0 / 16.0 + 1.05) <= M.C_B2_BND_HIT and (74.0 / 16.0 + 1.05) <= M.C_B2_BND_MISS
                       and segs.get("init8_u8") == 75 and segs.get("scan8_u8") == 199
                       and (74.0 + 197.0) / 1024.0 <= M.C_SCAN_G8 and 0.375 <= M.C_SCAN_ZERO)

  # 4h. Live depth-first post-order boot root construction + exact Karatsuba/sequential-NTT3/Garner workspace schedule:
  # (a) Production-scale byte-by-byte workspace schedule matching M.ntt3_F's m > 2^55 Karatsuba split across all 61 depths:
  #     At x_EA = 62.9375, m_root = ceil(pbits(0, 62.9375, 1.0)/64) = 5.5263*2^55 limbs (|P| = 2^60.4663 B = 1.382 EiB).
  #     For m <= 2^55: 3 NTT primes run SEQUENTIALLY sharing two in-place buffers T_A, T_B (8*N_ntt B each) plus two
  #     retained residue buffers R_0, R_1 (8*m B each); prime p3 streams Garner CRT word-by-word directly into P_out.
  #     CRT non-overflow for all m <= 2^55: (m/2)*(2^64-1)^2 < 2^182.0 < 2^182.696 = p1*p2*p3.
  #     For m > 2^55: 1-level Karatsuba split h=ceil(m/2) writes z0 into P_out[0..h-1] and z2 into P_out[h..2h-1], then
  #     allocates S_A, S_B, Z1 (2*8*(h+2) B) and recurses on h+2.
  pb_prod_df = M.pbits(0, 62.9375, 1.0)
  m_root_prod = ceil(pb_prod_df / 64.0)
  root_bytes_prod = 8 * m_root_prod
  def _ntt3_len(m_v):
    kv = max(2, ceil(log2(m_v)))
    n2v = 1 << kv
    n3v = 3 * (1 << (kv - 2)) if kv >= 3 else n2v
    return n3v if (kv >= 3 and n3v >= m_v) else n2v
  def _kar_ntt3_scratch(m_v, m_base=(1 << 55)):
    if m_v <= m_base:
      return 2 * (8 * _ntt3_len(m_v)) + 2 * (8 * m_v)
    hv = (m_v + 1) // 2
    return max(_kar_ntt3_scratch(hv, m_base), 2 * 8 * (hv + 2) + _kar_ntt3_scratch(hv + 2, m_base))
  prod_peak_bytes = 0
  for d_lv in range(61):
    st_l = sum(ceil(m_root_prod / (1 << i)) for i in range(1, d_lv + 1))
    m_dv = max(1, ceil(m_root_prod / (1 << d_lv)))
    sb_dv = 8 * st_l + 8 * (2 * m_dv) + _kar_ntt3_scratch(m_dv)
    if sb_dv > prod_peak_bytes:
      prod_peak_bytes = sb_dv
  p1_n, p2_n, p3_n = M.NTT3_PRIMES
  ntt3_crt_no_ovf = (((1 << 55) // 2) * (M64 * M64) < p1_n * p2_n * p3_n)
  # (b) Live bit-exact 3-level Karatsuba + sequential 3-prime NTT3 + streaming Garner CRT on 64-limb (4096-bit) product
  #     with explicit alloc/free byte tracking (m_base = 8 limbs so all 3 Karatsuba levels execute live):
  inv12_n = pow(p1_n, -1, p2_n); inv123_n = pow((p1_n * p2_n) % p3_n, -1, p3_n)
  _ntt_roots = {}
  for _pn in M.NTT3_PRIMES:
    _g = 3
    while pow(_g, (_pn - 1) // 2, _pn) == 1: _g += 1
    _ntt_roots[_pn] = _g
  def _ntt_r2(A_v, N_v, p_v, inv_f=False):
    wb = pow(_ntt_roots[p_v], (p_v - 1) // N_v, p_v)
    if inv_f: wb = pow(wb, -1, p_v)
    a_v = [x % p_v for x in A_v] + [0] * (N_v - len(A_v))
    jv = 0
    for iv in range(1, N_v):
      bit = N_v >> 1
      while jv & bit: jv ^= bit; bit >>= 1
      jv ^= bit
      if iv < jv: a_v[iv], a_v[jv] = a_v[jv], a_v[iv]
    ln = 2
    while ln <= N_v:
      wl = pow(wb, N_v // ln, p_v); hf = ln >> 1
      for iv in range(0, N_v, ln):
        wv = 1
        for kv in range(hf):
          uv = a_v[iv + kv]; vv = (a_v[iv + kv + hf] * wv) % p_v
          a_v[iv + kv] = (uv + vv) % p_v; a_v[iv + kv + hf] = (uv - vv) % p_v
          wv = (wv * wl) % p_v
      ln <<= 1
    if inv_f:
      invN = pow(N_v, -1, p_v)
      a_v = [(x * invN) % p_v for x in a_v]
    return a_v
  _live_b = [0]; _peak_b = [0]
  def _ab(nb):
    _live_b[0] += nb
    if _live_b[0] > _peak_b[0]: _peak_b[0] = _live_b[0]
  def _fb(nb): _live_b[0] -= nb
  def _mul_ntt3_seq(A_v, B_v, m_o):
    Nv = 1 << max(2, ceil(log2(max(1, m_o))))
    _ab(2 * 8 * Nv)
    ta = _ntt_r2(A_v, Nv, p1_n); tb = _ntt_r2(B_v, Nv, p1_n)
    r0 = _ntt_r2([(ta[i] * tb[i]) % p1_n for i in range(Nv)], Nv, p1_n, True)[:m_o]
    _ab(8 * m_o)
    ta = _ntt_r2(A_v, Nv, p2_n); tb = _ntt_r2(B_v, Nv, p2_n)
    r1 = _ntt_r2([(ta[i] * tb[i]) % p2_n for i in range(Nv)], Nv, p2_n, True)[:m_o]
    _ab(8 * m_o)
    ta = _ntt_r2(A_v, Nv, p3_n); tb = _ntt_r2(B_v, Nv, p3_n)
    r2 = _ntt_r2([(ta[i] * tb[i]) % p3_n for i in range(Nv)], Nv, p3_n, True)[:m_o]
    out_v = [0] * m_o; cy = 0
    for i in range(m_o):
      x1 = r0[i]; t2 = ((r1[i] - x1) * inv12_n) % p2_n; x12 = x1 + p1_n * t2
      t3 = ((r2[i] - x12) * inv123_n) % p3_n; val = x12 + (p1_n * p2_n) * t3 + cy
      out_v[i] = val & M64; cy = val >> 64
    _fb(2 * 8 * m_o + 2 * 8 * Nv)
    return out_v
  def _mul_kar_ntt3(A_v, B_v, m_o, m_base=8):
    if m_o <= m_base:
      return _mul_ntt3_seq(A_v, B_v, m_o)
    hv = (m_o + 1) // 2; hf = (hv + 1) // 2
    A0, A1 = A_v[:hf], A_v[hf:]; B0, B1 = B_v[:hf], B_v[hf:]
    z0 = _mul_kar_ntt3(A0, B0, 2 * hf, m_base)
    z2 = _mul_kar_ntt3(A1, B1, max(1, m_o - 2 * hf), m_base)
    sal = max(len(A0), len(A1)) + 1; sbl = max(len(B0), len(B1)) + 1; z1l = sal + sbl
    _ab(8 * (sal + sbl + z1l))
    SA = [0] * sal; cy = 0
    for i in range(sal):
      s = (A0[i] if i < len(A0) else 0) + (A1[i] if i < len(A1) else 0) + cy; SA[i] = s & M64; cy = s >> 64
    SB = [0] * sbl; cy = 0
    for i in range(sbl):
      s = (B0[i] if i < len(B0) else 0) + (B1[i] if i < len(B1) else 0) + cy; SB[i] = s & M64; cy = s >> 64
    Z1 = _mul_kar_ntt3(SA, SB, z1l, m_base)
    for src in (z0, z2):
      bw = 0
      for i in range(z1l):
        d = Z1[i] - (src[i] if i < len(src) else 0) - bw
        if d < 0: Z1[i] = d + (1 << 64); bw = 1
        else: Z1[i] = d; bw = 0
    out_v = [0] * m_o
    for i, v in enumerate(z0): out_v[i] = v
    for i, v in enumerate(z2):
      if 2 * hf + i < m_o: out_v[2 * hf + i] = v
    cy = 0
    for i, v in enumerate(Z1):
      if hf + i < m_o:
        s = out_v[hf + i] + v + cy; out_v[hf + i] = s & M64; cy = s >> 64
    kv = hf + z1l
    while cy and kv < m_o:
      s = out_v[kv] + cy; out_v[kv] = s & M64; cy = s >> 64; kv += 1
    _fb(8 * (sal + sbl + z1l))
    return out_v
  _kA = [rng.randrange(1 << 64) for _ in range(32)]
  _kB = [rng.randrange(1 << 64) for _ in range(32)]
  _ab(8 * (32 + 32 + 64))
  _kP = _mul_kar_ntt3(_kA, _kB, 64, m_base=8)
  _fb(8 * (32 + 32 + 64))
  def _l2u(L_v):
    acc = 0
    for i in range(len(L_v) - 1, -1, -1): acc = (acc << 64) | L_v[i]
    return acc
  kar_ntt3_live_ok = (_l2u(_kP) == _l2u(_kA) * _l2u(_kB) and _live_b[0] == 0 and _peak_b[0] < 8 * (8 * 64))
  # (c) Live 1,024-leaf depth-first post-order tree check using exact _kar_ntt3_scratch(m_step, m_base=8)
  df_leaves = M._primes_upto(8200)[:1024]
  assert len(df_leaves) == 1024
  lvl = list(df_leaves)
  full_tree_bits = sum(x.bit_length() for x in lvl)
  while len(lvl) > 1:
    lvl = [lvl[i] * lvl[i + 1] for i in range(0, len(lvl), 2)]
    full_tree_bits += sum(x.bit_length() for x in lvl)
  ref_root = lvl[0]
  root_bits = ref_root.bit_length()
  live_stack_bits = 0; max_stack_bits = 0
  def _df_prod(arr):
    nonlocal live_stack_bits, max_stack_bits
    if len(arr) == 1:
      return arr[0]
    mid = len(arr) // 2
    pl = _df_prod(arr[:mid]); bl = pl.bit_length()
    live_stack_bits += bl
    if live_stack_bits > max_stack_bits: max_stack_bits = live_stack_bits
    pr = _df_prod(arr[mid:])
    live_stack_bits -= bl
    return pl * pr
  df_root = _df_prod(df_leaves)
  boot_df_prod_peak_ratio = prod_peak_bytes / root_bytes_prod
  boot_df_tree_ok = (df_root == ref_root and live_stack_bits == 0
                       and max_stack_bits < 1.02 * root_bits
                       and ntt3_crt_no_ovf and kar_ntt3_live_ok
                       and prod_peak_bytes < 5.0 * root_bytes_prod
                       and 8 * root_bytes_prod < (1 << 64)
                       and full_tree_bits > 9.0 * root_bits)

  # 4i. Live bit-exact 4-kernel scalar limb-only Bernstein SRT leaf verification & dynamic op bound (<= M.C_SRT_G_LEAF):
  # Pipeline for any 0 < N_g < 2^512 (even OR odd, 8 unsigned 64-bit limbs):
  #   0. Pre-tree odd-part strip (srt_odd_norm8): N_g = 2^v * N_odd (N_odd strictly odd, < 2^512).
  #      Because 2 | P_24, 2^v is always P_24-smooth; stripping 2^v BEFORE entering the Bernstein product tree
  #      guarantees EVERY leaf modulus N_odd is odd and r_0 = P_24 mod N_odd satisfies 0 <= r_0 < N_odd.
  #   1. Montgomery setup (inv64_odd) + 9 CIOS8 squarings (cios8_mul): x_9 = r_0^512 * 2^{-511*512} mod N_odd.
  #      Because N_odd is odd, gcd(2^512, N_odd) = 1, so gcd(N_odd, x_9) == gcd(N_odd, r_0^512) identically.
  #   2. Register-resident 8-limb Binary GCD (bingcd8, Stein's algorithm; u[0..7] in x0..x7, v[0..7] in x8..x15):
  #      Self-contained 1-line termination proof: at every step with v > 0, either v is even (v' = v/2 => u'*v' = u*v/2)
  #      or v is odd (swap if u > v so u <= v, then v' = (v-u)/2 => v' = 0 or u'*v' = u*(v-u)/2 < u*v/2).
  #      Since u_0 * v_0 < 2^512 * 2^512 = 2^1024 and integer products >= 1 while v > 0, v == 0 in <= 1,024 steps!
  #   3. 8-limb Hensel odd-modulus exact division (hensel_div8): cof = N_odd / g (g is odd since g | N_odd).
  # Executed live using ONLY 64-bit limb ops (zero Python pow/math.gcd/division shortcuts inside the leaf) and
  # instruction-counted directly from the clang-assembled //@seg blocks in kernels_v32.s.
  def _u64(x): return x & M64
  def _mul64(a, b):
    p = (a & M64) * (b & M64)
    return p & M64, (p >> 64) & M64
  def _int_to_u64(x, n):
    x &= (1 << (64 * n)) - 1
    return [(x >> (64 * i)) & M64 for i in range(n)]
  def _u64_to_int(L):
    acc = 0
    for i in range(len(L) - 1, -1, -1):
      acc = (acc << 64) | (L[i] & M64)
    return acc
  def _srt_odd_norm8(A):
    k = 0
    while k < 8 and A[k] == 0: k += 1
    v = A[k]; s = (v & -v).bit_length() - 1
    out = [0] * 8
    for i in range(8):
      if i + k < 8:
        lo = A[i + k] >> s
        hi = (A[i + k + 1] << (64 - s)) & M64 if (s > 0 and i + k + 1 < 8) else 0
        out[i] = lo | hi
    return out
  def _inv64_odd(a):
    a &= M64; x = 1
    for _ in range(6):
      x = _u64(x * _u64(2 - _u64(a * x)))
    return x
  def _cios8_mul(A, B, N, n0inv):
    T = [0] * 10
    for i in range(8):
      ai = A[i]; c = 0
      for j in range(8):
        lo, hi = _mul64(ai, B[j])
        s = T[j] + lo + c; T[j] = s & M64; c = hi + (s >> 64)
      s = T[8] + c; T[8] = s & M64; T[9] = s >> 64
      m = _u64(T[0] * n0inv)
      lo, hi = _mul64(m, N[0])
      s = T[0] + lo; c = hi + (s >> 64)
      for j in range(1, 8):
        lo, hi = _mul64(m, N[j])
        s = T[j] + lo + c; T[j - 1] = s & M64; c = hi + (s >> 64)
      s = T[8] + c; T[7] = s & M64; T[8] = T[9] + (s >> 64); T[9] = 0
    D = [0] * 8; b = 0
    for j in range(8):
      d = T[j] - N[j] - b
      if d < 0: D[j] = d + (1 << 64); b = 1
      else: D[j] = d; b = 0
    return D if (T[8] - b) >= 0 else T[:8]
  def _sub8_limbs(A, B):
    D = [0] * 8; bw = 0
    for j in range(8):
      d = A[j] - B[j] - bw
      if d < 0: D[j] = d + (1 << 64); bw = 1
      else: D[j] = d; bw = 0
    return D, bw
  def _shr1_8_limbs(X):
    return [((X[j] >> 1) | ((X[j + 1] & 1) << 63)) for j in range(7)] + [X[7] >> 1]
  def _bingcd8(U_in, V_in):
    U = U_in[:]; V = V_in[:]
    iters = 0
    while any(x != 0 for x in V):
      iters += 1
      if iters > 1024:
        return None, iters
      if (V[0] & 1) == 0:
        V = _shr1_8_limbs(V)
      else:
        D, bw = _sub8_limbs(V, U)
        if bw:
          D, bw2 = _sub8_limbs(U, V)
          if bw2: return None, iters
          U = V[:]
        V = _shr1_8_limbs(D)
    return U, iters
  def _hensel_div8(N_odd_L, g_L):
    ginv = _inv64_odd(g_L[0]); rem = N_odd_L[:]; Q = [0] * 8
    for i in range(8):
      qi = _u64(rem[i] * ginv); Q[i] = qi
      cm = 0; bs = 0
      for j in range(8 - i):
        lo, hi = _mul64(qi, g_L[j])
        p = lo + cm; cm = hi + (p >> 64); p &= M64
        d = rem[i + j] - p - bs
        if d < 0: rem[i + j] = d + (1 << 64); bs = 1
        else: rem[i + j] = d; bs = 0
      if rem[i] != 0: return None
    return Q

  srt_9sq_errs = 0; bingcd8_max_iters = 0
  ps_small_s = M._primes_upto(64)
  P_small_s = 1
  for p_v in ps_small_s:
    P_small_s *= p_v
  fa, fb = 1, 1
  while fa + fb < (1 << 512): fa, fb = fb, fa + fb
  fib_odd, fib_prev = (fb, fa) if (fb & 1) else (fa, fb)
  srt_cases = [
    1 << 480,
    1 << 511,
    3 ** 302,
    (2 ** 199) * (3 ** 120) * (5 ** 40) * 1000000007,
    (1 << 479) - 1,
    (1 << 512) - 2,
    fib_odd,
    (fib_odd * 2) if (fib_odd * 2 < (1 << 512)) else (fib_odd - 1),
  ] + [rng.randrange(1 << 470, 1 << 512) for _ in range(56)]
  for xv in srt_cases:
    N_odd_L = _srt_odd_norm8(_int_to_u64(xv, 8))
    N_odd = _u64_to_int(N_odd_L)
    if (N_odd & 1) != 1:
      srt_9sq_errs += 1; continue
    r0 = P_small_s % N_odd
    if N_odd == 1:
      cof_got = 1
    else:
      n0inv = _u64(-_inv64_odd(N_odd_L[0]))
      X_L = _int_to_u64(r0, 8)
      for _ in range(9):
        X_L = _cios8_mul(X_L, X_L, N_odd_L, n0inv)
      g_L, it_v = _bingcd8(N_odd_L, X_L)
      if it_v > bingcd8_max_iters: bingcd8_max_iters = it_v
      if g_L is None or it_v > 1024:
        srt_9sq_errs += 1; continue
      cof_L = _hensel_div8(N_odd_L, g_L)
      if cof_L is None:
        srt_9sq_errs += 1; continue
      cof_got = _u64_to_int(cof_L)
    cur = xv
    for p_v in ps_small_s:
      while cur % p_v == 0:
        cur //= p_v
    if cof_got != cur:
      srt_9sq_errs += 1
  # Direct adversarial bingcd8 checks (v=0, u=v, gcd=1, gcd=u, 512-bit Fibonacci pair, max 512-bit odds):
  for ua, va in (
    (fib_odd, 0),
    (fib_odd, fib_odd),
    (fib_odd, fib_prev),
    ((1 << 512) - 1, (1 << 512) - 3),
    ((1 << 511) + 1, 1),
    (((1 << 255) - 19) * ((1 << 255) - 31), ((1 << 255) - 19) * ((1 << 255) - 101)),
  ):
    g_tL, it_t = _bingcd8(_int_to_u64(ua, 8), _int_to_u64(va, 8))
    if it_t > bingcd8_max_iters: bingcd8_max_iters = it_t
    if g_tL is None or it_t > 1024 or _u64_to_int(g_tL) != math.gcd(ua, va):
      srt_9sq_errs += 1

  # Exact deterministic worst-case dynamic scalar op count directly from kernels_v32.s //@seg blocks:
  _req_srt = (
    "srt_odd_norm8_find", "srt_odd_norm8_mid", "srt_odd_norm8_shift", "inv64_odd",
    "cios8_setup", "cios8_outer_pre", "cios8_inner_mul", "cios8_outer_mid",
    "cios8_inner_red", "cios8_outer_post", "cios8_csub_loop", "cios8_csel_loop",
    "bingcd8_io", "bingcd8_chk", "bingcd8_odd_sub", "bingcd8_swap", "bingcd8_shr_d", "bingcd8_shr_v",
    "hensel_div8_outer", "hensel_div8_inner",
  )
  _srt_segs_ok = all(segs.get(k, 0) > 0 for k in _req_srt)
  _ops_odd_norm = 8 * segs.get("srt_odd_norm8_find", 0) + segs.get("srt_odd_norm8_mid", 0) + 8 * segs.get("srt_odd_norm8_shift", 0) + 10
  _ops_cios8_one = (segs.get("cios8_setup", 0)
           + 8 * (segs.get("cios8_outer_pre", 0) + 8 * segs.get("cios8_inner_mul", 0)
                             + segs.get("cios8_outer_mid", 0) + 7 * segs.get("cios8_inner_red", 0)
                             + segs.get("cios8_outer_post", 0))
           + 8 * segs.get("cios8_csub_loop", 0) + 8 * segs.get("cios8_csel_loop", 0) + 20)
  _ops_cios8_9 = segs.get("inv64_odd", 0) + 2 + 9 * _ops_cios8_one
  _bg_iter_wc = max(
    segs.get("bingcd8_chk", 0) + segs.get("bingcd8_odd_sub", 0) + segs.get("bingcd8_swap", 0) + segs.get("bingcd8_shr_d", 0),
    segs.get("bingcd8_chk", 0) + 1 + segs.get("bingcd8_shr_v", 0),
  )
  _ops_bingcd8 = segs.get("bingcd8_io", 0) + 1024 * _bg_iter_wc + segs.get("bingcd8_chk", 0)
  _ops_hensel = segs.get("inv64_odd", 0) + 8 * segs.get("hensel_div8_outer", 0) + 36 * segs.get("hensel_div8_inner", 0) + 20
  srt_leaf_9sq_wc_ops = _ops_odd_norm + _ops_cios8_9 + _ops_bingcd8 + _ops_hensel  # 61,542 ops <= 63,000!
  srt_g_leaf_ok = (_srt_segs_ok and srt_9sq_errs == 0 and bingcd8_max_iters <= 1024 and srt_leaf_9sq_wc_ops <= M.C_SRT_G_LEAF)

  # 4j. Exact Comba schoolbook (mul4_sb in kernels_v32.s) + recursive Karatsuba limb-array verification (200 random tests):
  def _limbs_to_int(L):
    acc = 0
    for i in range(len(L) - 1, -1, -1):
      acc = (acc << 64) | L[i]
    return acc
  def _sb_limbs(A, B):
    na, nb = len(A), len(B)
    C = [0] * (na + nb)
    for i in range(na):
      ai = A[i]; carry = 0
      for j in range(nb):
        p = ai * B[j] + C[i + j] + carry
        C[i + j] = p & M64
        carry = p >> 64
      C[i + nb] = carry
    return C
  def _add_limbs(A, B):
    n = max(len(A), len(B))
    out = [0] * (n + 1); c = 0
    for i in range(n):
      s = (A[i] if i < len(A) else 0) + (B[i] if i < len(B) else 0) + c
      out[i] = s & M64; c = s >> 64
    out[n] = c
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out
  def _sub_limbs(A, B):
    out = [0] * len(A); b = 0
    for i in range(len(A)):
      d = A[i] - (B[i] if i < len(B) else 0) - b
      if d < 0: out[i] = d + (1 << 64); b = 1
      else: out[i] = d; b = 0
    assert b == 0
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out
  def _kar_limbs(A, B):
    n = max(len(A), len(B))
    if n <= 16:
      return _sb_limbs(A, B)
    m = n // 2
    A0, A1 = A[:m], A[m:]
    B0, B1 = B[:m], B[m:]
    z0 = _kar_limbs(A0, B0)
    z2 = _kar_limbs(A1, B1)
    z1 = _sub_limbs(_sub_limbs(_kar_limbs(_add_limbs(A0, A1), _add_limbs(B0, B1)), z0), z2)
    out = [0] * (len(A) + len(B))
    for off, arr in ((0, z0), (m, z1), (2 * m, z2)):
      c = 0
      for i, v in enumerate(arr):
        if off + i < len(out):
          s = out[off + i] + v + c
          out[off + i] = s & M64; c = s >> 64
      k = off + len(arr)
      while c and k < len(out):
        s = out[k] + c; out[k] = s & M64; c = s >> 64; k += 1
    return out
  mul_kar_sb_errs = 0
  for idx_k in range(200):
    n_l = (4, 8, 12, 16, 20, 24, 32, 48, 64)[idx_k % 9]
    A_l = [rng.randrange(1 << 64) for _ in range(n_l)]
    B_l = [rng.randrange(1 << 64) for _ in range(n_l)]
    got_p = _limbs_to_int(_kar_limbs(A_l, B_l))
    ref_p = _limbs_to_int(A_l) * _limbs_to_int(B_l)
    if got_p != ref_p:
      mul_kar_sb_errs += 1
  mul4_sb_ops = segs.get("mul4_sb", -1)
  mul_kar_sb_ok = (mul_kar_sb_errs == 0 and mul4_sb_ops == 105 and mul4_sb_ops <= int(M._kar_sb(4)))

  # 4k. Kit A (v37): Hybrid 4-clamp [0,4p) lazy radix-4 Shoup NTT (46.5*n*k) + {1,3,5}*2^k 2-prime 52-bit NTT:
  _v37 = _verify_v37_ntt(M, segs, rng)
  v37_ntt_ok = _v37["ok"]

  ok = (t3["ok"] and ntt3_split_errs == 0 and ctail_errs == 0 and ctail_ops <= M.OP["c_tail_z6"]
     and bfi_b2_errs == 0 and bfi_b2_cases == 768 and delta_b1_errs == 0 and u8_range_ok and scan8_u8_errs == 0
     and b1_delta_ops_ok and boot_df_tree_ok and srt_g_leaf_ok and mul_kar_sb_ok and v37_ntt_ok)
  return dict(
    ok=ok,
    adic3_exact=t3,
    ntt3_split_proto_errs=ntt3_split_errs,
    ctail_lazy_errs=ctail_errs,
    ctail_lazy_ops=ctail_ops,
    ctail_budget=M.OP["c_tail_z6"],
    bfi_b2_errs=bfi_b2_errs,
    bfi_b2_cases=bfi_b2_cases,
    b2_prefix_mult=round(_b2_pref_mult, 4),
    b2_active_ks=_act_ks,
    delta_b1_errs=delta_b1_errs,
    ng_diff_ub_bits=round(ng_diff_ub, 6),
    u8_scale=S_U8, u8_thr_max=thr_max, u8_init_min=init_min, u8_max_val=max_u8,
    ret2d_m_margin_bits=round(r2d_chk["m_margin"], 6),
    ret2d_rg_marg_err=round(rg_marg_err, 8),
    ret2d_sf_marg_err=round(sf_marg_err, 8),
    ret2d_per_prime_errs=ret2d_per_prime_errs,
    ret2d_mc_ok=ret2d_mc_ok,
    ret2d_mc_emp_rg=round(mc_emp_rg, 4),
    ret2d_mc_mod_rg=round(mc_mod_rg, 4),
    ret2d_mc_emp_sf=round(mc_emp_sf, 4),
    ret2d_mc_mod_sf=round(mc_mod_sf, 4),
    rho_tail_ok=rho_tail_ok,
    rho_335=f"{rho_335:.6e}",
    rho_400=f"{rho_400:.6e}",
    rot_tab_2s_errs=rot_tab_errs,
    step1_tab_8192_inv_errs=step1_inv_errs,
    sf_uses_survivor_upper_ok=sf_uses_survivor_upper_ok,
    u16_range_ok=u8_range_ok,
    scan8_u16_errs=scan8_u8_errs,
    b1_delta_ops_ok=b1_delta_ops_ok,
    boot_df_tree_ok=boot_df_tree_ok,
    boot_df_stack_ratio=round(max_stack_bits / root_bits, 4),
    boot_df_peak_ratio=round(boot_df_prod_peak_ratio, 4),
    boot_df_prod_peak_log2=round(log2(prod_peak_bytes), 4),
    boot_kar_ntt3_live_peak_ratio=round(_peak_b[0] / (8 * 64), 4),
    ntt3_crt_no_ovf=ntt3_crt_no_ovf,
    boot_full_tree_ratio=round(full_tree_bits / root_bits, 2),
    srt_g_leaf_ok=srt_g_leaf_ok,
    srt_g_leaf_9sq_wc_ops=int(srt_leaf_9sq_wc_ops),
    srt_g_leaf_budget=M.C_SRT_G_LEAF,
    bingcd8_max_iters=bingcd8_max_iters,
    bingcd8_iter_wc_ops=_bg_iter_wc,
    mul_kar_sb_ok=mul_kar_sb_ok,
    mul_kar_sb_errs=mul_kar_sb_errs,
    mul4_sb_ops=mul4_sb_ops,
    v37_ntt_ok=v37_ntt_ok,
    v37_r4_inner=_v37["r4_in"], v37_r4_hoist_x2=_v37["r4_ho"], v37_r3=_v37["r3"], v37_r5=_v37["r5"],
    v37_r4_blend=_v37["c_bl"], v37_r3_lin=_v37["c_r3"], v37_r5_lin=_v37["c_r5"],
    v37_bf_errs=_v37["bf_errs"], v37_conv_errs=_v37["conv_errs"],
  )


def _verify_v37_ntt(M, segs, rng):
  M64 = (1 << 64) - 1
  r4_in, r4_ho = segs.get("r4_4clamp_inner_tw", -1), segs.get("r4_4clamp_hoist_tw_x2", -1)
  r3_bf, r5_bf = segs.get("r3_shoup_bfly", -1), segs.get("r5_shoup_bfly", -1)
  c_bl = 0.5 * (9.0 * r4_in / 8.0 + 9.0 * r4_ho / 16.0)
  c_r3, c_r5 = 9.0 * r3_bf / 3.0 + 220.0, 9.0 * r5_bf / 5.0 + 220.0
  def is_p(n):
    if n < 2: return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
      if n == q: return True
      if n % q == 0: return False
    d, s = n - 1, 0
    while (d & 1) == 0: d >>= 1; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
      x = pow(a, d, n)
      if x in (1, n - 1): continue
      for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1: break
      else: return False
    return True
  def proot(p):
    phi, fs, t = p - 1, set(), p - 1
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61):
      while t % q == 0: fs.add(q); t //= q
    d = 67
    while d * d <= t:
      while t % d == 0: fs.add(d); t //= d
      d += 2
    if t > 1: fs.add(t)
    g = 2
    while any(pow(g, phi // q, p) == 1 for q in fs): g += 1
    return g
  def sm(v, w, wp, p): return (v * w - ((v * wp) >> 64) * p) & M64
  def c2(x, p2): return x - p2 if x >= p2 else x
  def bf4(a, ws, J, Jp, p):
    p2 = 2 * p
    t0 = c2(a[0], p2); t1 = sm(a[1], ws[0], (ws[0] << 64) // p, p)
    t2 = sm(a[2], ws[1], (ws[1] << 64) // p, p); t3 = sm(a[3], ws[2], (ws[2] << 64) // p, p)
    u0, u1 = c2(t0 + t2, p2), c2(t0 + p2 - t2, p2)
    v0, v1 = c2(t1 + t3, p2), sm(t1 + p2 - t3, J, Jp, p)
    return (u0 + v0, u1 + v1, u0 + p2 - v0, u1 + p2 - v1)
  def bf3(a, ws, c1, c1p, c2v, c2p, p):
    p2 = 2 * p
    t0 = c2(a[0], p2); t1 = sm(a[1], ws[0], (ws[0] << 64) // p, p); t2 = sm(a[2], ws[1], (ws[1] << 64) // p, p)
    u = c2(t1 + t2, p2); m1 = sm(u, c1, c1p, p); m2 = sm(t1 + p2 - t2, c2v, c2p, p); r = c2(t0 + m1, p2)
    return (t0 + u, r + m2, r + p2 - m2)
  def r3_c(w3, p):
    i2 = (p + 1) // 2; c1 = (p - i2) % p; c2v = ((w3 - (w3 * w3) % p) * i2) % p
    return c1, (c1 << 64) // p, c2v, (c2v << 64) // p
  def r5_c(w5, p):
    w2 = (w5 * w5) % p; w3 = (w2 * w5) % p; w4 = (w3 * w5) % p; i2 = (p + 1) // 2; i4 = (i2 * i2) % p
    c0 = (-i4) % p; c1 = ((w5 + w4 - w2 - w3) * i4) % p; d1 = ((w5 - w4) * i2) % p; d2 = ((w2 - w3) * i2) % p
    k0, k1, k2 = d2, (d1 - d2) % p, (-(d1 + d2)) % p
    return tuple((x, (x << 64) // p) for x in (c0, c1, k0, k1, k2))
  def bf5(a, ws, cs, p):
    p2 = 2 * p; (c0, c0p), (c1, c1p), (k0, k0p), (k1, k1p), (k2, k2p) = cs
    t0 = c2(a[0], p2); ts = [sm(a[i + 1], ws[i], (ws[i] << 64) // p, p) for i in range(4)]
    s1, s2 = c2(ts[0] + ts[3], p2), c2(ts[1] + ts[2], p2)
    d1, d2 = c2(ts[0] + p2 - ts[3], p2), c2(ts[1] + p2 - ts[2], p2)
    s12 = c2(s1 + s2, p2); u1 = c2(t0 + sm(s12, c0, c0p, p), p2); m1 = sm(s1 + p2 - s2, c1, c1p, p)
    r1, r2 = c2(u1 + m1, p2), c2(u1 + p2 - m1, p2)
    q0 = sm(d1 + d2, k0, k0p, p); z1, z2 = c2(q0 + sm(d1, k1, k1p, p), p2), c2(q0 + sm(d2, k2, k2p, p), p2)
    return (t0 + s12, r1 + z1, r2 + z2, r2 + p2 - z2, r1 + p2 - z1)
  p1, p2 = M.NTT2_52_PRIMES
  primes_ok = (all(is_p(p) and (p - 1) % (3 * (1 << 55)) == 0 and 4 * p < (1 << 64) for p in M.NTT3_PRIMES)
               and all(is_p(p) and (p - 1) % (15 * (1 << 40)) == 0 and 4 * p < (1 << 64) for p in M.NTT2_52_PRIMES)
               and (p1 * p2) > (1 << 18) * ((1 << 52) - 1) ** 2 * 3)
  proots = {p: proot(p) for p in M.NTT3_PRIMES + M.NTT2_52_PRIMES}
  bf_errs = 0
  for p in M.NTT3_PRIMES + M.NTT2_52_PRIMES:
    p4, g = 4 * p, proots[p]; J = pow(g, (p - 1) // 4, p); Jp = (J << 64) // p
    for idx in range(516):
      xs = [p4 - 1 if (idx >> b) & 1 else 0 for b in range(4)] if idx < 16 else [rng.randrange(p4) for _ in range(4)]
      ws = [rng.randrange(p) for _ in range(3)]; ys = bf4(xs, ws, J, Jp, p)
      aa = [xs[0] % p] + [(xs[i + 1] * ws[i]) % p for i in range(3)]
      if max(ys) >= p4 or [y % p for y in ys] != [(aa[0] + aa[2] + aa[1] + aa[3]) % p, (aa[0] - aa[2] + J * (aa[1] - aa[3])) % p, (aa[0] + aa[2] - aa[1] - aa[3]) % p, (aa[0] - aa[2] - J * (aa[1] - aa[3])) % p]:
        bf_errs += 1
    w3 = pow(g, (p - 1) // 3, p); c3 = r3_c(w3, p)
    for idx in range(258):
      xs = [p4 - 1 if (idx >> b) & 1 else 0 for b in range(3)] if idx < 8 else [rng.randrange(p4) for _ in range(3)]
      ws = [rng.randrange(p) for _ in range(2)]; ys = bf3(xs, ws, *c3, p)
      aa = [xs[0] % p] + [(xs[i + 1] * ws[i]) % p for i in range(2)]
      if max(ys) >= p4 or [y % p for y in ys] != [(aa[0] + aa[1] + aa[2]) % p, (aa[0] + w3 * aa[1] + w3 * w3 * aa[2]) % p, (aa[0] + w3 * w3 * aa[1] + w3 * aa[2]) % p]:
        bf_errs += 1
    if (p - 1) % 5 == 0:
      w5 = pow(g, (p - 1) // 5, p); c5 = r5_c(w5, p)
      for idx in range(282):
        xs = [p4 - 1 if (idx >> b) & 1 else 0 for b in range(5)] if idx < 32 else [rng.randrange(p4) for _ in range(5)]
        ws = [rng.randrange(p) for _ in range(4)]; ys = bf5(xs, ws, c5, p)
        aa = [xs[0] % p] + [(xs[i + 1] * ws[i]) % p for i in range(4)]
        if max(ys) >= p4 or [y % p for y in ys] != [sum(aa[j] * pow(w5, i * j, p) for j in range(5)) % p for i in range(5)]:
          bf_errs += 1
  def mr_ntt(a, p, inv=False):
    n = len(a); wn = pow(proots[p], (p - 1) // n, p)
    if inv: wn = pow(wn, p - 2, p)
    def rec(v, ws):
      m = len(v)
      if m == 1: return [v[0]]
      r = 5 if m % 5 == 0 else (3 if m % 3 == 0 else (4 if m % 4 == 0 else 2))
      sb = m // r; pts = [rec(v[i::r], pow(ws, r, p)) for i in range(r)]; out = [0] * m; wk = 1
      w_r = pow(ws, sb, p)
      c3 = r3_c(w_r, p) if r == 3 else None; c5 = r5_c(w_r, p) if r == 5 else None; Jp = (w_r << 64) // p if r == 4 else 0
      for j in range(sb):
        pw = [wk]
        for _ in range(r - 2): pw.append((pw[-1] * wk) % p)
        inp = [pts[i][j] for i in range(r)]
        if r == 5: ys = bf5(inp, pw, c5, p)
        elif r == 3: ys = bf3(inp, pw, *c3, p)
        elif r == 4: ys = bf4(inp, pw, w_r, Jp, p)
        else:
          t0, t1 = c2(inp[0], 2 * p), sm(inp[1], wk, (wk << 64) // p, p)
          ys = (t0 + t1, t0 + 2 * p - t1)
        for i in range(r): out[j + i * sb] = ys[i]
        wk = (wk * ws) % p
      return out
    res = rec([x % p for x in a], wn)
    if inv:
      iv = pow(n, p - 2, p)
      return [(x * iv) % p for x in res]
    return [x % p for x in res]
  conv_errs = 0; inv12 = pow(p1, p2 - 2, p2); B52 = 1 << 52
  for nv in (12, 20, 24, 30, 48, 60):
    hf = nv // 2
    for md in (0, 1):
      A = ([B52 - 1] * hf if md == 0 else [rng.randrange(B52) for _ in range(hf)]) + [0] * (nv - hf)
      B = ([B52 - 1] * hf if md == 0 else [rng.randrange(B52) for _ in range(hf)]) + [0] * (nv - hf)
      c1 = mr_ntt([(x * y) % p1 for x, y in zip(mr_ntt(A, p1), mr_ntt(B, p1))], p1, True)
      c2v = mr_ntt([(x * y) % p2 for x, y in zip(mr_ntt(A, p2), mr_ntt(B, p2))], p2, True)
      pi = cy = 0
      for i in range(nv):
        s = c1[i] + p1 * (((c2v[i] - c1[i]) * inv12) % p2) + cy
        pi |= (s & (B52 - 1)) << (52 * i); cy = s >> 52
      pi |= cy << (52 * nv)
      if pi != sum(A[i] << (52 * i) for i in range(hf)) * sum(B[i] << (52 * i) for i in range(hf)):
        conv_errs += 1
  for p in M.NTT3_PRIMES:
    for nv in (12, 24):
      A, B = [rng.randrange(p) for _ in range(nv)], [rng.randrange(p) for _ in range(nv)]
      cc = mr_ntt([(x * y) % p for x, y in zip(mr_ntt(A, p), mr_ntt(B, p))], p, True)
      if cc != [sum(A[j] * B[(i - j) % nv] for j in range(nv)) % p for i in range(nv)]:
        conv_errs += 1
  ok = (r4_in == 45 and r4_ho == 74 and r3_bf == 36 and r5_bf == 93
        and c_bl <= 46.5 and c_r3 <= 364.0 and c_r5 <= 420.0
        and primes_ok and bf_errs == 0 and conv_errs == 0)
  return dict(ok=ok, r4_in=r4_in, r4_ho=r4_ho, r3=r3_bf, r5=r5_bf,
              c_bl=round(c_bl, 4), c_r3=round(c_r3, 2), c_r5=round(c_r5, 2),
              bf_errs=bf_errs, conv_errs=conv_errs)
