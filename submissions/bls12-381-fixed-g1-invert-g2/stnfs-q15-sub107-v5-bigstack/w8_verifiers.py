"""w8_verifiers.py — Live verifiers for BLS12-381 Q(zeta_15) CM-STNFS:
  1. verify_threshold_counterexample(): exact algebraic witness that naive PGF substitution G(psi) fails.
  2. verify_genuine_q15_filter_hx_live():
     - 100% UNSCREENED, ZERO-ORACLE, TWO-SIDED Q(zeta_15) Z[zeta_12] filtering + UNPLANTED Block Wiedemann H-X
       + post-hoc independent BSGS dlog_GT verification (34/34 exact on BOTH f and g sides) + algebraic gamma_0 blinding.
     - Matched-scale genuine-vs-synthetic Z[zeta_12] calibration (exact same R0, act_cols, row-wt dist, 1/p col wts).
     - TRUE same-underlying-matrix full-incidence (0 omitted cols) vs truncated-proxy (middle band -> w_mid)
       comparison on ONE concrete Z[zeta_12] matrix M_full.
"""
import bisect
import hashlib
import heapq
import math
import random
from filtsim_z12 import (
U12, U12_SET, U12_INV, z12_add, z12_sub, z12_mul, z12_is_zero, z12_maxabs,
z12_eval_mod, _vinfo,
)
from hx_z12 import _right_nullspace, _rank_rows
def _isp_det(n):
 if n < 2:
  return False
 for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
  if n == p:
   return True
  if n % p == 0:
   return False
 d, s = n - 1, 0
 while d % 2 == 0:
  d //= 2
  s += 1
 for a in (2, 3, 5, 7, 11, 13):
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
def verify_threshold_counterexample():
 """Exact algebraic witness that naive PGF substitution G(psi) fails for general non-i.i.d. laws."""
 mu, rv, psi = 60.0, 1.0, 0.10
 bt = 0.5 * mu / rv
 G = lambda t: (1.0 + bt * (1.0 - t)) ** (-rv)
 for i in range(81):
  t = i / 100.0
  assert 0.10 * t + 0.90 * (t ** 100) <= G(t) + 1e-12
 assert 0.10 > G(psi) + 0.05
 best = G(0.0) + psi
 for k in range(2, 65):
  tk = min(0.8, (k - 1) * (1.0 + 1.0 / bt) / (k - 1 + rv))
  best = min(best, min(1.0, G(tk) / (tk ** (k - 1))) + (psi ** k))
 assert best >= 0.10 - 1e-12
 return dict(pgf_psi_invalid=round(G(psi), 5), true_fail=0.10, tr_bound=round(best, 5), tr_valid=True)
def _same_matrix_full_vs_proxy_z12(U=1024, seed=20261121, W_lim=60):
 """TRUE same-underlying-matrix comparison:
    Generate ONE concrete Z[zeta_12] matrix M_full with ALL light + middle + heavy columns explicitly
    retained (NO 0.90 incidence cap, NO omitted middle band), then run the exact same zero-column-loss
    _dc-guarded filter on:
      (A) full_explicit: 100% of M_full (0 omitted columns, w_mid = 0),
      (B) proxy_same_matrix: EXACT SAME M_full with its top-15% densest middle-band columns removed from
          the sparse graph and charged as static per-row w_mid weight (stopping at N_sw / N_l_eff >= W_lim).
    """
 rng = random.Random(seed)
 R_raw = int(1.65 * U)
 wts = [1.0 / (c * math.log(c + 3.0) + 7.0) for c in range(U)]
 sw = sum(wts)
 cum = []; acc = 0.0
 for w in wts:
  acc += w / sw; cum.append(acc)
 M_full = []
 for _ in range(R_raw):
  rd = {}
  k_ent = rng.randint(10, 18)
  for _ in range(k_ent):
   c = min(U - 1, bisect.bisect_left(cum, rng.random()))
   u = U12[rng.randrange(12)]
   nv = z12_add(rd.get(c, (0, 0, 0, 0)), u)
   if z12_is_zero(nv):
    rd.pop(c, None)
   else:
    rd[c] = nv
  if len(rd) >= 2:
   M_full.append(rd)
 deg0 = [0] * U
 for rd in M_full:
  for c in rd:
   deg0[c] += 1
 s_degs = sorted(d for d in deg0 if d >= 2)
 d_cut = s_degs[int(0.85 * len(s_degs))]
 mid_cols = {c for c in range(U) if deg0[c] >= d_cut}
 n_mid = len(mid_cols)
 def _run_z12_filter(rows_in, w_mid_in, extra_dense_cols):
  rows = [dict(r) for r in rows_in]
  w_mid = list(w_mid_in)
  al = [True] * len(rows)
  cr = {}
  for ri, rd in enumerate(rows):
   for c in rd:
    cr.setdefault(c, set()).add(ri)
  q1 = [c for c, s in cr.items() if len(s) == 1]
  while q1:
   c = q1.pop()
   if len(cr[c]) != 1:
    continue
   ri = next(iter(cr[c]))
   if not al[ri]:
    continue
   al[ri] = False
   for cc in list(rows[ri]):
    cr[cc].discard(ri)
    if len(cr[cc]) == 1:
     q1.append(cc)
   rows[ri].clear()
  cols_2core = {c for c, s in cr.items() if len(s) >= 2}
  n_2core = len(cols_2core)
  repurge = []
  alive_idx = [ri for ri in range(len(rows)) if al[ri] and rows[ri]]
  alive_idx.sort(key=lambda ri: len(rows[ri]) + w_mid[ri], reverse=True)
  n_al = len(alive_idx); targ = int(1.25 * n_2core)
  for ri in alive_idx:
   if n_al <= targ:
    break
   if not al[ri]:
    continue
   if all(len(cr[c]) >= 2 for c in rows[ri]):
    al[ri] = False; n_al -= 1
    q_rp = []
    for c in list(rows[ri]):
     cr[c].discard(ri)
     if len(cr[c]) == 1:
      q_rp.append(c)
    rows[ri].clear()
    while q_rp:
     cc = q_rp.pop()
     if len(cr[cc]) == 1:
      r2 = next(iter(cr[cc]))
      if al[r2] and all(len(cr[ccc]) >= 2 for ccc in rows[r2] if ccc != cc):
       al[r2] = False; n_al -= 1; repurge.append(cc)
       for ccc in list(rows[r2]):
        cr[ccc].discard(r2)
        if len(cr[ccc]) == 1:
         q_rp.append(ccc)
       rows[r2].clear()
  pq = [(len(cr[c]), c) for c in cols_2core if 1 <= len(cr[c]) <= 12]
  heapq.heapify(pq)
  merged = []
  n_sw = sum(len(rows[ri]) + w_mid[ri] for ri in range(len(rows)) if al[ri])
  n_l = sum(1 for c in cols_2core if cr[c])
  while pq and (n_sw / max(1, n_l + extra_dense_cols)) < W_lim:
   d0, c = heapq.heappop(pq); d = len(cr[c])
   if d == 0 or d != d0 or d > 12:
    continue
   if d == 1:
    ri = next(iter(cr[c]))
    if al[ri] and rows[ri].get(c) in U12_SET and all(len(cr[cc]) >= 2 for cc in rows[ri] if cc != c):
     merged.append(c); al[ri] = False
     n_sw -= (len(rows[ri]) + w_mid[ri]); n_l -= 1
     for cc in list(rows[ri]):
      cr[cc].discard(ri)
      if 1 <= len(cr[cc]) <= 12:
       heapq.heappush(pq, (len(cr[cc]), cc))
     rows[ri].clear()
    continue
   piv = None; pw = 10**9
   for ri in cr[c]:
    if rows[ri][c] in U12_SET:
     w = len(rows[ri]) + w_mid[ri]
     if w < pw:
      pw = w; piv = ri
   if piv is None:
    continue
   pr = rows[piv]; invp = U12_INV[pr[c]]; p_wm = w_mid[piv]
   oth = [ri for ri in cr[c] if ri != piv]
   ok = True; up = []; _dc = {cc: 1 for cc in pr if cc != c}
   for ri in oth:
    fc = z12_mul(rows[ri][c], invp); nr = dict(rows[ri])
    for cc, pv in pr.items():
     had = cc in rows[ri]
     nv = z12_sub(nr.get(cc, (0, 0, 0, 0)), z12_mul(fc, pv))
     if z12_is_zero(nv):
      nr.pop(cc, None)
      if cc != c and had:
       _dc[cc] = _dc.get(cc, 0) + 1
     elif z12_maxabs(nv) > 255:
      ok = False; break
     else:
      nr[cc] = nv
      if cc != c and not had:
       _dc[cc] = _dc.get(cc, 0) - 1
    if not ok:
     break
    nwm = min(n_mid, w_mid[ri] + p_wm) if (p_wm or w_mid[ri]) else 0
    up.append((ri, nr, nwm))
   if ok and any(len(cr[cc]) <= d1 for cc, d1 in _dc.items() if d1 > 0):
    ok = False
   if not ok:
    continue
   merged.append(c); al[piv] = False; n_l -= 1
   n_sw -= (len(pr) + p_wm)
   for cc in list(pr):
    cr[cc].discard(piv)
   rows[piv].clear()
   for ri, nr, nwm in up:
    n_sw += (len(nr) + nwm) - (len(rows[ri]) + w_mid[ri])
    old_c, new_c = set(rows[ri]), set(nr)
    for cc in old_c - new_c:
     cr[cc].discard(ri)
    for cc in new_c - old_c:
     cr[cc].add(ri)
    rows[ri] = nr; w_mid[ri] = nwm
   for cc in set(pr) | {k for _, nr, _ in up for k in nr}:
    if 1 <= len(cr[cc]) <= 12:
     heapq.heappush(pq, (len(cr[cc]), cc))
  live_c = {c for c in cols_2core if cr[c]}
  lost = n_2core - (len(live_c) + len(repurge) + len(merged))
  tot_e = sum(len(rows[ri]) for ri in range(len(rows)) if al[ri])
  unit_e = sum(1 for ri in range(len(rows)) if al[ri] for v in rows[ri].values() if v in U12_SET)
  from filtsim_z12 import _vinfo
  et = eo = ec = 0.0
  for ri in range(len(rows)):
   if not al[ri]:
    continue
   for v in rows[ri].values():
    for tv in _vinfo(v)[3]:
     if tv:
      av = -tv if tv < 0 else tv
      et += 1.0
      if av == 1:
       eo += 1.0
      ec += (7.65 * av) if av <= 4 else 22.0
  exp_u = eo / max(1.0, et)
  k5_ex = ec / max(1.0, et)
  k5_ub = exp_u * 7.65 + (1.0 - exp_u) * 22.0
  tot_wm = sum(w_mid[ri] for ri in range(len(rows)) if al[ri])
  return dict(
  n_2core=n_2core, N_l_eff=len(live_c) + extra_dense_cols, N_l_sparse=len(live_c),
  N_sw_rows=sum(1 for x in al if x), N_sw_entries=tot_e + tot_wm,
  unit_frac_ring=round(unit_e / max(1, tot_e), 4),
  exp_unit_hex6=round(exp_u, 5), k5_exact=round(k5_ex, 4), k5_ub=round(k5_ub, 4),
  lost=lost, merged=len(merged),
  )
 r_full = _run_z12_filter(M_full, [0] * len(M_full), 0)
 M_prox = [{c: v for c, v in rd.items() if c not in mid_cols} for rd in M_full]
 W_prox = [sum(1 for c in rd if c in mid_cols) for rd in M_full]
 r_prox = _run_z12_filter(M_prox, W_prox, n_mid)
 assert r_full["lost"] == 0 and r_prox["lost"] == 0
 assert r_full["N_l_eff"] <= r_prox["N_l_eff"]
 return dict(
 same_underlying_matrix=True, U=U, R_raw=len(M_full), mid_cols_omitted_in_proxy=n_mid,
 d_cut_mid=d_cut, W_lim=W_lim, full_explicit=r_full, proxy_same_matrix=r_prox,
 full_nl_le_proxy_nl=(r_full["N_l_eff"] <= r_prox["N_l_eff"]),
 full_entries_over_proxy_entries=round(r_full["N_sw_entries"] / max(1, r_prox["N_sw_entries"]), 4),
 )
def verify_genuine_q15_filter_hx_live(
toy, primes_upto, PHI15, y_pow_z, sub_phi15_z, mul_phi15_z, pow_phi15_z,
scale_phi15_z, norm_q15_z, pgcd_mod, pmul_mod, prem_mod, u12_proxy=None,
):
 """100% UNSCREENED, ZERO-ORACLE, TWO-SIDED Q(zeta_15) Z[zeta_12] filtering + UNPLANTED Block Wiedemann H-X
    + matched-scale genuine-vs-synthetic calibration + TRUE same-underlying-matrix full-vs-proxy check."""
 P0, R0, CHI = 7, 181, 7
 W_EL = y_pow_z(5); W2 = y_pow_z(10); ONE = [1, 0, 0, 0, 0, 0, 0, 0]
 CG_EL = W2; CF_EL = mul_phi15_z(sub_phi15_z(W_EL, ONE), W2)
 w_pi = next(w for w in range(2, P0) if (w * w + w + 1) % P0 == 0)
 h_pi = pgcd_mod(PHI15, [-w_pi % P0, 0, 0, 0, 0, 1], P0)
 assert len(h_pi) == 5 and not any(prem_mod(sub_phi15_z(CF_EL, CG_EL), h_pi, P0))
 assert pow((w_pi * w_pi) % P0, (P0 - 1) // 3, P0) == w_pi
 w_tau_exp = 5
 B_FB = 2000
 ps_1 = [p for p in primes_upto(B_FB) if p % 15 == 1 and p != R0]
 ps_van = [2, 3, 5] + [p for p in primes_upto(B_FB) if p % 15 in (4, 7, 13)]
 P_van = math.prod(ps_van); P_1 = math.prod(ps_1)
 def strip_all(n):
  if n <= 0 or n % R0 == 0:
   return None
  t = n
  while True:
   g = math.gcd(t, P_van)
   if g == 1:
    break
   t //= g
  if t == 1:
   return []
  r = t; hits = []
  g1 = math.gcd(r, P_1)
  if g1 == 1:
   return None
  for p in ps_1:
   if g1 % p == 0:
    e = 0
    while r % p == 0:
     r //= p; e += 1
    hits.append((p, e))
    if r == 1:
     return hits
  return hits if r == 1 else None
 pool = [[1 if j == i else 0 for j in range(8)] for i in range(8)]
 for i in range(8):
  for j in range(i + 1, 8):
   for s in (1, -1):
    v = [0] * 8; v[i] = 1; v[j] = s; pool.append(v)
 for i in range(8):
  for j in range(8):
   if i != j:
    for s in (1, -1):
     v = [0] * 8; v[i] = 2; v[j] = s; pool.append(v)
 for i in range(8):
  for j in range(i + 1, 8):
   for k in range(j + 1, 8):
    for s1 in (1, -1):
     for s2 in (1, -1):
      v = [0] * 8; v[i] = 1; v[j] = s1; v[k] = s2
      if abs(norm_q15_z(v)) <= 220:
       pool.append(v)
 seen_v = set(); uniq = []
 for v in pool:
  t = tuple(v) if v >= [-x for x in v] else tuple(-x for x in v)
  if t not in seen_v:
   seen_v.add(t); uniq.append(list(t))
 def canon_y(v):
  best = None; cur = list(v)
  for _ in range(15):
   for s in (1, -1):
    t = tuple(s * x for x in cur)
    if best is None or t < best:
     best = t
   cur = mul_phi15_z(cur, [0, 1, 0, 0, 0, 0, 0, 0])
  return best
 seen_b = set(); b_uniq = []
 for v in uniq:
  c = canon_y(v)
  if c not in seen_b:
   seen_b.add(c); b_uniq.append(v)
 pre_a = [(v, pow_phi15_z(v, 3), abs(norm_q15_z(v))) for v in uniq]
 pre_b = [(v, b3, abs(norm_q15_z(v)), mul_phi15_z(CG_EL, b3), mul_phi15_z(CF_EL, b3)) for v in b_uniq for b3 in (pow_phi15_z(v, 3),)]
 raw = []; used_p = set()
 for a8, a3, na in pre_a:
  for b8, b3, nb, cgb3, cfb3 in pre_b:
   if a8 == b8:
    continue
   ng = abs(norm_q15_z(sub_phi15_z(a3, cgb3)))
   if ng <= 1:
    continue
   fg = strip_all(ng)
   if fg is None:
    continue
   nf = abs(norm_q15_z(sub_phi15_z(a3, cfb3)))
   if nf <= 1:
    continue
   ff = strip_all(nf)
   if ff is None or (not fg and not ff):
    continue
   ar = prem_mod(a8, h_pi, P0); br = prem_mod(b8, h_pi, P0)
   if not any(ar) and not any(br):
    continue
   g0 = math.gcd(na, nb)
   if g0 > 1:
    ok_cp = True; t = g0; d = 2
    while d * d <= t:
     if t % d == 0:
      while t % d == 0:
       t //= d
      g1 = pgcd_mod(a8, PHI15, d)
      if len(g1) > 1 and len(pgcd_mod(g1, b8, d)) > 1:
       ok_cp = False; break
     d += 1 if d == 2 else 2
    if ok_cp and t > 1:
     g1 = pgcd_mod(a8, PHI15, t)
     if len(g1) > 1 and len(pgcd_mod(g1, b8, t)) > 1:
      ok_cp = False
    if not ok_cp:
     continue
   raw.append((a8, b8, fg, ff))
   for p, _ in fg:
    used_p.add(p)
   for p, _ in ff:
    used_p.add(p)
 orbs = []; by_sp = {}
 for p in sorted(used_p):
  e15 = (p - 1) // 15
  g0 = next(pow(g, e15, p) for g in range(2, p) if pow(pow(g, e15, p), 3, p) != 1 and pow(pow(g, e15, p), 5, p) != 1)
  yrs = [pow(g0, k, p) for k in (1, 2, 4, 7, 8, 11, 13, 14)]
  for sd, C_el in ((0, CG_EL), (1, CF_EL)):
   idp = set()
   for ry in yrs:
    rhs = sum(C_el[j] * pow(ry, j, p) for j in range(8)) % p
    wp = pow(ry, 5, p); w2p = (wp * wp) % p
    if pow(rhs, (p - 1) // 3, p) == 1:
     for x in range(1, p):
      if pow(x, 3, p) == rhs:
       idp.update(((ry, x), (ry, (x * wp) % p), (ry, (x * w2p) % p))); break
   sn = set()
   for ry0, rx0 in sorted(idp):
    if (ry0, rx0) in sn:
     continue
    ob = []; cy, cx = ry0, rx0
    for _ in range(12):
     pw = (1, cy, (cy*cy)%p, pow(cy,3,p), pow(cy,4,p), pow(cy,5,p), pow(cy,6,p), pow(cy,7,p))
     ob.append((cy, cx, pw)); sn.add((cy, cx))
     cy, cx = pow(cy, 7, p), (cx * pow(cy, w_tau_exp, p)) % p
    oid = len(orbs); orbs.append((sd, p, ob)); by_sp.setdefault((sd, p), []).append(oid)
 def lift_root(ry, rx, C_el, p, e):
  mod = p; cy, cx = ry, rx
  for _ in range(1, e):
   nmod = mod * p
   fy = sum(PHI15[j] * pow(cy, j, nmod) for j in range(9)) % nmod
   dfy = sum(j * PHI15[j] * pow(cy, j - 1, p) for j in range(1, 9)) % p
   cy = (cy - fy * pow(dfy, -1, p)) % nmod
   rhs = sum(C_el[j] * pow(cy, j, nmod) for j in range(8)) % nmod
   fx = (pow(cx, 3, nmod) - rhs) % nmod
   dfx = (3 * cx * cx) % p
   cx = (cx - fx * pow(dfx, -1, p)) % nmod
   mod = nmod
  return cy, cx, mod
 def fac_z12(a8, b8, fs, sd, C_el):
  row = {}
  a0,a1,a2,a3,a4,a5,a6,a7 = a8; b0,b1,b2,b3,b4,b5,b6,b7 = b8
  for p, ev in fs:
   rem = ev
   for oid in by_sp.get((sd, p), ()):
    for k, (ry, rx, pw) in enumerate(orbs[oid][2]):
     av = (a0 + a1*pw[1] + a2*pw[2] + a3*pw[3] + a4*pw[4] + a5*pw[5] + a6*pw[6] + a7*pw[7]) % p
     bv = (b0 + b1*pw[1] + b2*pw[2] + b3*pw[3] + b4*pw[4] + b5*pw[5] + b6*pw[6] + b7*pw[7]) % p
     if (av - bv * rx) % p == 0:
      vp = 1
      if ev > 1:
       rye, rxe, pmod = lift_root(ry, rx, C_el, p, ev)
       val = (sum(a8[j] * pow(rye, j, pmod) for j in range(8)) - sum(b8[j] * pow(rye, j, pmod) for j in range(8)) * rxe) % pmod
       vp = 0; t = p
       while vp < ev and val % t == 0:
        vp += 1; t *= p
      u = U12[(-k) % 12]
      nv = row.get(oid, (0, 0, 0, 0))
      for _ in range(vp):
       nv = z12_add(nv, u)
      if z12_is_zero(nv):
       row.pop(oid, None)
      else:
       row[oid] = nv
      rem -= vp
   if rem != 0:
    return None
  return row
 rels = []; seen_sig = set()
 for a8, b8, fg, ff in raw:
  zg = fac_z12(a8, b8, fg, 0, CG_EL); zf = fac_z12(a8, b8, ff, 1, CF_EL)
  if zg is None or zf is None:
   continue
  zrow = dict(zf)
  for o, zv in zg.items():
   zrow[o] = tuple(-x for x in zv)
  if not zrow:
   continue
  first_o = min(zrow)
  inv0 = pow(z12_eval_mod(zrow[first_o], CHI, R0), -1, R0)
  sig = tuple(sorted((o, (z12_eval_mod(zv, CHI, R0) * inv0) % R0) for o, zv in zrow.items()))
  if sig in seen_sig:
   continue
  seen_sig.add(sig); rels.append((a8, b8, zrow, zf, zg))
 while True:
  deg = {}
  for _, _, zr, _, _ in rels:
   for o in zr:
    deg[o] = deg.get(o, 0) + 1
  sing = {o for o, d in deg.items() if d < 2}
  if not sing:
   break
  rels = [r for r in rels if not (set(r[2]) & sing)]
 act = sorted({o for _, _, zr, _, _ in rels for o in zr})
 nf_o = sum(1 for o in act if orbs[o][0] == 1)
 ng_o = sum(1 for o in act if orbs[o][0] == 0)
 assert len(rels) >= 30 and len(act) >= 24 and nf_o >= 10 and ng_o >= 10
 r0sq = R0 * R0
 roots_r0 = [s for s in range(2, R0) if pow(s, 15, R0) == 1 and pow(s, 5, R0) != 1 and pow(s, 3, R0) != 1]
 roots_l2 = [
 (s - (sum(PHI15[j]*pow(s, j, r0sq) for j in range(9)) % r0sq)
 * pow(sum(j*PHI15[j]*pow(s, j-1, R0) for j in range(1, 9)) % R0, -1, R0)) % r0sq
 for s in roots_r0
 ]
 crt_e = []
 for idx, s in enumerate(roots_r0):
  poly = [1]; den = 1
  for j, sj in enumerate(roots_r0):
   if j != idx:
    poly = pmul_mod(poly, [-sj % R0, 1], R0); den = (den * (s - sj)) % R0
  invd = pow(den, -1, R0); p8 = [(x * invd) % R0 for x in poly]; p8 += [0]*(8 - len(p8)); crt_e.append(p8)
 chi_inv_pows = [pow(CHI, (-t) % 12, R0) for t in range(12)]
 s2_pows = [[pow(s2, j, r0sq) for j in range(8)] for s2 in roots_l2]
 yp15 = [y_pow_z(k) for k in range(15)]
 def sm_3(a8, b8, C_el):
  eps = R0**3 - 1; L_c = [[0]*8 for _ in range(3)]
  for ci, spw in enumerate(s2_pows):
   ac = sum(a8[j] * spw[j] for j in range(8)) % r0sq
   bc = sum(b8[j] * spw[j] for j in range(8)) % r0sq
   c0 = sum(C_el[j] * spw[j] for j in range(8)) % r0sq
   p0, p1, p2 = 1, 0, 0
   q0, q1, q2 = ac, (-bc) % r0sq, 0
   ee = eps
   while ee:
    if ee & 1:
     p0, p1, p2 = (
     (p0*q0 + c0*(p1*q2 + p2*q1)) % r0sq,
     (p0*q1 + p1*q0 + c0*p2*q2) % r0sq,
     (p0*q2 + p1*q1 + p2*q0) % r0sq,
     )
    ee >>= 1
    if ee:
     q0, q1, q2 = (
     (q0*q0 + c0*(2*q1*q2)) % r0sq,
     (2*q0*q1 + c0*q2*q2) % r0sq,
     (2*q0*q2 + q1*q1) % r0sq,
     )
   pw = [(p0 - 1) % r0sq, p1, p2]; ec = crt_e[ci]
   for k in range(3):
    qk = (pw[k] // R0) % R0
    if qk:
     for j in range(8):
      L_c[k][j] = (L_c[k][j] + qk * ec[j]) % R0
  proj = [0] * 24
  for k in range(3):
   ck = L_c[k]
   if not any(ck):
    continue
   acc = [0]*8; cur_p = list(ck)
   for t in range(12):
    wt = chi_inv_pows[t]; shift = yp15[(w_tau_exp * t * k) % 15]; term = mul_phi15_z(cur_p, shift)
    for j in range(8):
     acc[j] = (acc[j] + wt * term[j]) % R0
    nxt = [0]*8
    for j, cj in enumerate(cur_p):
     if cj:
      y7j = yp15[(7 * j) % 15]
      for ii in range(8):
       nxt[ii] = (nxt[ii] + cj * y7j[ii]) % R0
    cur_p = nxt
   for j in range(8):
    proj[k*8 + j] = acc[j]
  return proj
 smf_raw = [sm_3(a8, b8, CF_EL) for a8, b8, _, _, _ in rels]
 smg_raw = [sm_3(a8, b8, CG_EL) for a8, b8, _, _, _ in rels]
 rng_sm = random.Random(20261121)
 Wf = [rng_sm.randrange(R0) for _ in range(24)]
 Wg = [rng_sm.randrange(R0) for _ in range(24)]
 s_sm = 2
 sms_all = [
 [sum(Wf[k] * smf_raw[ri][k] for k in range(24)) % R0,
 (-sum(Wg[k] * smg_raw[ri][k] for k in range(24))) % R0]
 for ri in range(len(rels))
 ]
 om = {o: i for i, o in enumerate(act)}
 nc_full = len(act) + s_sm
 M_full = []
 for ri, (_, _, zr, _, _) in enumerate(rels):
  r = [0] * nc_full
  for o, zv in zr.items():
   r[om[o]] = z12_eval_mod(zv, CHI, R0)
  r[-2] = sms_all[ri][0]; r[-1] = sms_all[ri][1]
  M_full.append(r)
 rk_M_full = _rank_rows([r[:len(act)] for r in M_full], R0)
 rk_A_full, ns_full = _right_nullspace(M_full, nc_full, R0)
 assert rk_M_full == len(act) and rk_A_full == nc_full - 1 and len(ns_full) == 1
 rows = [dict(r[2]) for r in rels]
 sms = [list(s) for s in sms_all]
 al = [True] * len(rows)
 cr = {c: set() for c in act}
 for ri, rd in enumerate(rows):
  for c in rd:
   cr[c].add(ri)
 ncols_2core = len(act); merge_stack = []
 target_live = max(12, ncols_2core - 10)
 pq = [(len(cr[c]), c) for c in act if 1 <= len(cr[c]) <= 14]
 heapq.heapify(pq)
 while pq and sum(1 for c in act if cr[c]) > target_live:
  d0, c = heapq.heappop(pq); d = len(cr[c])
  if d == 0 or d != d0 or d > 14:
   continue
  if d == 1:
   ri = next(iter(cr[c]))
   if al[ri] and rows[ri].get(c) in U12_SET and all(len(cr[cc]) >= 2 for cc in rows[ri] if cc != c):
    merge_stack.append((c, dict(rows[ri]), list(sms[ri])))
    al[ri] = False
    for cc in list(rows[ri]):
     cr[cc].discard(ri)
     if 1 <= len(cr[cc]) <= 14:
      heapq.heappush(pq, (len(cr[cc]), cc))
    rows[ri].clear()
   continue
  piv = None; pw = 10**9
  for ri in sorted(cr[c]):
   if rows[ri][c] in U12_SET:
    rw = sum(_vinfo(v)[0] for v in rows[ri].values())
    if rw < pw:
     pw = rw; piv = ri
  if piv is None:
   continue
  pr = rows[piv]; invp = U12_INV[pr[c]]; psm = sms[piv]
  oth = [ri for ri in sorted(cr[c]) if ri != piv]
  ok = True; up = []; _dc = {cc: 1 for cc in pr if cc != c}
  for ri in oth:
   fc = z12_mul(rows[ri][c], invp); fc_q = z12_eval_mod(fc, CHI, R0)
   nr = dict(rows[ri])
   for cc, pv in pr.items():
    had = cc in rows[ri]
    nv = z12_sub(nr.get(cc, (0, 0, 0, 0)), z12_mul(fc, pv))
    if z12_is_zero(nv):
     nr.pop(cc, None)
     if cc != c and had:
      _dc[cc] = _dc.get(cc, 0) + 1
    elif z12_maxabs(nv) > 255:
     ok = False; break
    else:
     nr[cc] = nv
     if cc != c and not had:
      _dc[cc] = _dc.get(cc, 0) - 1
   if not ok:
    break
   nsm = [(sms[ri][s] - fc_q * psm[s]) % R0 for s in range(s_sm)]
   up.append((ri, nr, nsm))
  if ok and any(len(cr[cc]) <= d1 for cc, d1 in _dc.items() if d1 > 0):
   ok = False
  if not ok:
   continue
  merge_stack.append((c, dict(pr), list(psm)))
  al[piv] = False
  for cc in list(pr):
   cr[cc].discard(piv)
  rows[piv].clear()
  for ri, nr, nsm in up:
   old_cs = set(rows[ri]); new_cs = set(nr)
   for cc in old_cs - new_cs:
    cr[cc].discard(ri)
   for cc in new_cs - old_cs:
    cr[cc].add(ri)
   rows[ri] = nr; sms[ri] = nsm
  for cc in set(pr) | {k for _, nr, _ in up for k in nr}:
   if 1 <= len(cr[cc]) <= 14:
    heapq.heappush(pq, (len(cr[cc]), cc))
 live_cols = sorted(c for c in act if cr[c])
 cols_lost = ncols_2core - (len(live_cols) + len(merge_stack))
 assert cols_lost == 0 and len(merge_stack) >= 6
 live_idx = [ri for ri in range(len(rows)) if al[ri] and rows[ri]]
 tot_entries = sum(len(rows[ri]) for ri in live_idx)
 unit_entries = sum(1 for ri in live_idx for v in rows[ri].values() if v in U12_SET)
 gen_exp_unit = unit_entries / max(1, tot_entries)
 gen_k5_hex6 = gen_exp_unit * 7.65 + (1.0 - gen_exp_unit) * 22.0
 assert gen_k5_hex6 < 10.5
 N = len(live_cols); D = N + s_sm
 assert len(live_idx) >= D
 lc_ix = {c: j for j, c in enumerate(live_cols)}
 _ldeg = {}
 for ri in live_idx:
  for c in rows[ri]:
   _ldeg[c] = _ldeg.get(c, 0) + 1
 _keep_d = set(live_idx)
 for ri in sorted(live_idx, key=lambda i: (-len(rows[i]), -i)):
  if len(_keep_d) <= D:
   break
  if all(_ldeg[c] >= 2 for c in rows[ri]):
   _keep_d.remove(ri)
   for c in rows[ri]:
    _ldeg[c] -= 1
 sel_idx = [ri for ri in live_idx if ri in _keep_d]
 assert len(sel_idx) == D and all(_ldeg[c] >= 1 for c in live_cols)
 A_live = []
 for ri in sel_idx:
  r = [0] * D
  for c, zv in rows[ri].items():
   r[lc_ix[c]] = z12_eval_mod(zv, CHI, R0)
  for s in range(s_sm):
   r[N + s] = sms[ri][s]
  A_live.append(r)
 Mc = [[A_live[i][j] for i in range(D)] for j in range(N)]
 S = [[A_live[i][N + t] for i in range(D)] for t in range(s_sm)]
 rM = _rank_rows(Mc, R0); rA = _rank_rows(A_live, R0)
 assert rM == N and rA == D - 1
 def op(x):
  return [sum(A_live[i][j] * x[j] for j in range(N)) % R0 for i in range(D)]
 g = random.Random(20261121); m, n = 6, 4; d = -(-D // n) + 1; e = -(-D // m) + 1
 Z_cols = [[g.randrange(R0) for _ in range(N)] for _ in range(n - s_sm)]
 Y = [S[t] for t in range(s_sm)] + [op(zc) for zc in Z_cols]
 Bv = lambda y: op(y[:N])
 K_rows = []; R_blocks = [Y]; cur = [Bv(y) for y in Y]
 for _ in range(d):
  K_rows.extend(cur); R_blocks.append(cur); cur = [Bv(y) for y in cur]
 rC = _rank_rows(K_rows, R0)
 def Xof(onepad):
  idx = list(range(N)); g.shuffle(idx); blk = [[] for _ in range(m)]
  for j, c in enumerate(idx):
   blk[j % m].append(c)
  for t in range(s_sm):
   blk[0 if onepad else t % m].append(N + t)
  return [[(g.randrange(R0) if c in set(b) else 0) for c in range(D)] for b in blk]
 def build_Oe(X):
  Ok = []; Xi = [r[:] for r in X]
  for _ in range(e):
   Ok.extend(Xi)
   Xi = [[sum(r[i]*A_live[i][k] for i in range(D)) % R0 if k < N else 0 for k in range(D)] for r in Xi]
  return Ok
 Oe_good = build_Oe(Xof(False)); Oe_bad = build_Oe(Xof(True))
 rOe_good = _rank_rows(Oe_good, R0); rOe_bad = _rank_rows(Oe_bad, R0)
 rO = _rank_rows([[sum(r[i]*A_live[i][k] for i in range(D)) % R0 for k in range(D)] for r in Oe_good], R0)
 R_cols = [col for blk in R_blocks for col in blk]
 H_rows = [[sum(r[i]*R_cols[c][i] for i in range(D)) % R0 for c in range(len(R_cols))] for r in Oe_good]
 rH, ker_H = _right_nullspace(H_rows, len(R_cols), R0)
 assert rC == N and rOe_good == D and rOe_bad < D and rO == D - 1 and rH == D - 1
 f_sel = next((fv for fv in ker_H if any(fv[:s_sm])), None)
 assert f_sel is not None
 c_rec = f_sel[:s_sm]; e_rec = f_sel[s_sm:n]
 acc = [0] * D
 for i in range(d, 0, -1):
  fi = f_sel[i*n:(i+1)*n]
  yi = [sum(Y[t][r]*fi[t] for t in range(n)) % R0 for r in range(D)]
  acc = [(Bv(acc)[r] + yi[r]) % R0 for r in range(D)] if i < d else yi
 w_rec = [(acc[j] + sum(Z_cols[t][j]*e_rec[t] for t in range(n - s_sm))) % R0 for j in range(N)]
 vlog = {live_cols[j]: w_rec[j] for j in range(N)}
 for piv_c, ents, psm in reversed(merge_stack):
  piv_q = z12_eval_mod(ents[piv_c], CHI, R0)
  acc_v = sum(psm[s] * c_rec[s] for s in range(s_sm)) % R0
  for cc, zv in ents.items():
   if cc != piv_c:
    acc_v = (acc_v + z12_eval_mod(zv, CHI, R0) * vlog[cc]) % R0
  vlog[piv_c] = (-pow(piv_q, -1, R0) * acc_v) % R0
 _h0, _h1, _h2, _h3 = h_pi[:4]
 def fp4_add(a, b): return [(a[i] + b[i]) % P0 for i in range(4)]
 def fp4_mul(a, b):
  a0, a1, a2, a3 = a; b0, b1, b2, b3 = b
  c0 = a0*b0; c1 = a0*b1+a1*b0; c2 = a0*b2+a1*b1+a2*b0; c3 = a0*b3+a1*b2+a2*b1+a3*b0
  c4 = a1*b3+a2*b2+a3*b1; c5 = a2*b3+a3*b2; c6 = a3*b3
  v = c6 % P0
  if v: c2 -= v*_h0; c3 -= v*_h1; c4 -= v*_h2; c5 -= v*_h3
  v = c5 % P0
  if v: c1 -= v*_h0; c2 -= v*_h1; c3 -= v*_h2; c4 -= v*_h3
  v = c4 % P0
  if v: c0 -= v*_h0; c1 -= v*_h1; c2 -= v*_h2; c3 -= v*_h3
  return [c0 % P0, c1 % P0, c2 % P0, c3 % P0]
 def fp4_pow(a, e):
  res, b = [1, 0, 0, 0], list(a)
  while e:
   if e & 1: res = fp4_mul(res, b)
   e >>= 1
   if e: b = fp4_mul(b, b)
  return res
 def fp4_inv(a):
  a1 = fp4_pow(a, P0); a2 = fp4_pow(a1, P0); a3 = fp4_pow(a2, P0)
  p3 = fp4_mul(a1, fp4_mul(a2, a3)); inv = pow(fp4_mul(a, p3)[0], -1, P0)
  return [(x * inv) % P0 for x in p3]
 w4 = fp4_pow([0, 1, 0, 0], 5); w4_2 = fp4_mul(w4, w4); uw2 = list(w4_2)
 frob4_cols = [fp4_pow([1 if j == i else 0 for j in range(4)], P0) for i in range(4)]
 def frob4(a):
  out = [0, 0, 0, 0]
  for i, ai in enumerate(a):
   if ai:
    col = frob4_cols[i]
    for j in range(4): out[j] = (out[j] + ai * col[j]) % P0
  return out
 def fp12_mul(A, B):
  return (fp4_add(fp4_mul(A[0], B[0]), fp4_mul(uw2, fp4_add(fp4_mul(A[1], B[2]), fp4_mul(A[2], B[1])))),
  fp4_add(fp4_add(fp4_mul(A[0], B[1]), fp4_mul(A[1], B[0])), fp4_mul(uw2, fp4_mul(A[2], B[2]))),
  fp4_add(fp4_add(fp4_mul(A[0], B[2]), fp4_mul(A[1], B[1])), fp4_mul(A[2], B[0])))
 def fp12_pow(A, e):
  res, b = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]), A
  while e:
   if e & 1: res = fp12_mul(res, b)
   e >>= 1
   if e: b = fp12_mul(b, b)
  return res
 def fp12_frob(A): return (frob4(A[0]), fp4_mul(frob4(A[1]), w4), fp4_mul(frob4(A[2]), w4_2))
 def fp12_inv(A):
  def f4(Z):
   for _ in range(4): Z = fp12_frob(Z)
   return Z
  A4 = f4(A); A8 = f4(A4); adj = fp12_mul(A4, A8); inv4 = fp4_inv(fp12_mul(A, adj)[0])
  return (fp4_mul(adj[0], inv4), fp4_mul(adj[1], inv4), fp4_mul(adj[2], inv4))
 cof = (P0**4 - P0**2 + 1) // R0
 def to_GT(A):
  e6 = A
  for _ in range(6): e6 = fp12_frob(e6)
  ez = fp12_mul(e6, fp12_inv(A)); e2 = fp12_frob(fp12_frob(ez))
  return fp12_pow(fp12_mul(e2, ez), cof)
 G_gen = to_GT(([3, 1, 0, 0], [2, 0, P0 - 1, 0], [0, 5, 0, 0]))
 m_bs = math.isqrt(R0) + 1; bs_map = {}; cur = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
 for j in range(m_bs):
  bs_map[(tuple(cur[0]), tuple(cur[1]), tuple(cur[2]))] = j
  cur = fp12_mul(cur, G_gen)
 Gm_inv = fp12_inv(fp12_pow(G_gen, m_bs))
 def dlog_GT(a8, b8):
  ar = prem_mod(a8, h_pi, P0); ar += [0] * (4 - len(ar))
  br = prem_mod([-x for x in b8], h_pi, P0); br += [0] * (4 - len(br))
  c = to_GT((ar, br, [0, 0, 0, 0]))
  for i in range(m_bs + 1):
   k = (tuple(c[0]), tuple(c[1]), tuple(c[2]))
   if k in bs_map: return (i * m_bs + bs_map[k]) % R0
   c = fp12_mul(c, Gm_inv)
  raise RuntimeError("BSGS failed")
 lam = None; exact_f = exact_g = nz_dl = 0
 for ri, (a8, b8, zr, zf, zg) in enumerate(rels):
  vf = (sum(z12_eval_mod(zv, CHI, R0) * vlog[c] for c, zv in zf.items()) + sms_all[ri][0] * c_rec[0]) % R0
  vg = (sum(z12_eval_mod(zv, CHI, R0) * vlog[c] for c, zv in zg.items()) - sms_all[ri][1] * c_rec[1]) % R0
  assert vf == vg
  dl = dlog_GT(a8, b8)
  if dl != 0:
   nz_dl += 1
   if lam is None:
    lam = (vf * pow(dl, -1, R0)) % R0
  if lam is not None:
   if vf == (lam * dl) % R0: exact_f += 1
   if vg == (lam * dl) % R0: exact_g += 1
 assert lam is not None and lam != 0 and exact_f == len(rels) and exact_g == len(rels)
 assert (P0 ** 12 - 1) % R0 == 0
 _ecof = (P0 ** 12 - 1) // R0
 assert math.gcd(R0, _ecof) == 1
 def _rel_to_GT(a8_, b8_):
  ar_ = prem_mod(a8_, h_pi, P0); ar_ += [0] * (4 - len(ar_))
  br_ = prem_mod([-x for x in b8_], h_pi, P0); br_ += [0] * (4 - len(br_))
  return to_GT((ar_, br_, [0, 0, 0, 0]))
 _vgs_all = [(sum(z12_eval_mod(zv, CHI, R0) * vlog[c] for c, zv in zg.items()) - sms_all[ri][1] * c_rec[1]) % R0 for ri, (_, _, _, _, zg) in enumerate(rels)]
 _nz_idx = [i for i, v in enumerate(_vgs_all) if v != 0]
 assert len(_nz_idx) >= 4
 _i0, _ir1, _ir2 = _nz_idx[0], _nz_idx[1], _nz_idx[-1]
 gamma_0 = _rel_to_GT(rels[_i0][0], rels[_i0][1])
 V_0 = (_ecof * _vgs_all[_i0]) % R0
 _id12 = ([1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0])
 assert gamma_0 != _id12 and fp12_pow(gamma_0, R0) == _id12 and V_0 != 0
 _e1, _e2 = 37, 113
 _r1_gt = _rel_to_GT(rels[_ir1][0], rels[_ir1][1]); _r2_gt = _rel_to_GT(rels[_ir2][0], rels[_ir2][1])
 _h_el = fp12_mul(_r1_gt, fp12_pow(gamma_0, (R0 - _e1) % R0))
 _T_el = fp12_mul(_r2_gt, fp12_pow(gamma_0, (R0 - _e2) % R0))
 assert fp12_mul(_h_el, fp12_pow(gamma_0, _e1)) == _r1_gt and fp12_mul(_T_el, fp12_pow(gamma_0, _e2)) == _r2_gt
 _V_hp = (_ecof * _vgs_all[_ir1]) % R0; _V_Tp = (_ecof * _vgs_all[_ir2]) % R0
 _V_h_unbl = (_V_hp - _e1 * V_0) % R0; _V_T_unbl = (_V_Tp - _e2 * V_0) % R0
 assert _V_h_unbl != 0
 _ratio_unbl = (_V_T_unbl * pow(_V_h_unbl, -1, R0)) % R0
 assert fp12_pow(_h_el, _ratio_unbl) == _T_el
 wts_m = [1.0 / orbs[o][1] for o in act]; sw_m = sum(wts_m); cum_m = []; acc_m = 0.0
 for w in wts_m:
  acc_m += w / sw_m; cum_m.append(acc_m)
 syn_uf = []; syn_k5 = []
 for sd_i in range(4):
  sr = random.Random(20261121 + 100 + sd_i); s_rows = []
  for _, _, zr, _, _ in rels:
   d_r = {}
   for _ in range(len(zr)):
    ci = act[min(len(act) - 1, bisect.bisect_left(cum_m, sr.random()))]
    nv = z12_add(d_r.get(ci, (0, 0, 0, 0)), U12[sr.randrange(12)])
    if z12_is_zero(nv):
     d_r.pop(ci, None)
    else:
     d_r[ci] = nv
   if len(d_r) >= 2:
    s_rows.append(d_r)
  sc = {c: set() for c in act}
  sal = [True] * len(s_rows)
  for ri, rd in enumerate(s_rows):
   for c in rd:
    sc[c].add(ri)
  for c in sorted(act, key=lambda cc: len(sc[cc])):
   if len(sc[c]) != 2:
    continue
   r1, r2 = sorted(sc[c])
   if not (sal[r1] and sal[r2]) or s_rows[r1].get(c) not in U12_SET:
    continue
   fc = z12_mul(s_rows[r2][c], U12_INV[s_rows[r1][c]])
   for cc, pv in list(s_rows[r1].items()):
    nv = z12_sub(s_rows[r2].get(cc, (0, 0, 0, 0)), z12_mul(fc, pv))
    if z12_is_zero(nv):
     s_rows[r2].pop(cc, None); sc[cc].discard(r2)
    else:
     s_rows[r2][cc] = nv; sc[cc].add(r2)
   sal[r1] = False
   for cc in list(s_rows[r1]):
    sc[cc].discard(r1)
  st = sum(len(s_rows[ri]) for ri in range(len(s_rows)) if sal[ri])
  su = sum(1 for ri in range(len(s_rows)) if sal[ri] for v in s_rows[ri].values() if v in U12_SET)
  uf_s = su / max(1, st)
  syn_uf.append(uf_s); syn_k5.append(uf_s * 7.65 + (1.0 - uf_s) * 22.0)
 same_mat = _same_matrix_full_vs_proxy_z12(U=1024, seed=20261121, W_lim=60)
 same_mat_ladder = {"U1024_W60": same_mat}
 _prev_rat = same_mat["full_entries_over_proxy_entries"]
 for _u_sm, _w_sm in ((4096, 200), (8192, 200), (16384, 200), (32768, 200)):
  _sm_r = _same_matrix_full_vs_proxy_z12(U=_u_sm, seed=20261121, W_lim=_w_sm)
  assert _sm_r["full_nl_le_proxy_nl"] and _sm_r["full_explicit"]["lost"] == 0 and _sm_r["proxy_same_matrix"]["lost"] == 0
  assert _sm_r["full_entries_over_proxy_entries"] < _prev_rat
  _prev_rat = _sm_r["full_entries_over_proxy_entries"]
  same_mat_ladder[f"U{_u_sm}_W{_w_sm}"] = _sm_r
 return dict(
 tower="unscreened_zero_oracle_2sided_Q15_p7_r181_chi7",
 polys="g=X^3-w^2, f=X^3-(w-1)*w^2 (congruent mod pi=(7,h_pi), non-isomorphic N(w-1)=3)",
 oracle_free_matrix_filter_solver=True,
 bad_cols_screened=0,
 rref_row_selection=False,
 rels_2core=len(rels),
 orbits_2core=ncols_2core,
 orbits_f=nf_o,
 orbits_g=ng_o,
 rk_M_2core=rk_M_full,
 rk_A_2core=rk_A_full,
 nullity_2core=len(ns_full),
 merged_pivots=len(merge_stack),
 live_cols=N,
 cols_2core_lost=cols_lost,
 s_sm=s_sm,
 held_out_surviving_rows=len(live_idx) - D,
 gen_exp_unit=round(gen_exp_unit, 4),
 gen_k5_hex6=round(gen_k5_hex6, 3),
 syn_matched_exp_unit=round(sum(syn_uf) / 4.0, 4),
 syn_matched_k5_hex6=round(sum(syn_k5) / 4.0, 3),
 bw_rank_Cd=rC,
 bw_rank_Oe_good=rOe_good,
 bw_rank_Oe_bad=rOe_bad,
 bw_rank_OeA=rO,
 bw_rank_H=rH,
 posthoc_bsgs_lam=lam,
 posthoc_bsgs_nz_dlogs=nz_dl,
 posthoc_exact_f=exact_f, gamma0_blind_unblind_lambda_cancel_ok=True,
 posthoc_exact_g=exact_g,
 posthoc_total_rels=len(rels),
 same_matrix_full_vs_proxy=same_mat,
 same_matrix_full_vs_proxy_ladder=same_mat_ladder,
 checks_ok=True,
 )
def check_sigmaB_exactness_live(build_G1w, build_WS_1w, mat_mul, M_1W, norm_q15_z, U_BLS):
 import struct
 G = build_G1w()
 W, S = build_WS_1w()
 def mv(M, v): return [sum(M[i][j] * v[j] for j in range(8)) for i in range(8)]
 def q_1w(v): return sum(v[i] * G[i][j] * v[j] for i in range(8) for j in range(8))
 I8 = [[1 if i == j else 0 for j in range(8)] for i in range(8)]
 G12_6 = []
 cur_s = I8
 for _ in range(4):
  cur_w = cur_s
  for _ in range(3):
   assert all(cur_w[i][j] == 0 for i in range(2, 8) for j in range(2))
   G12_6.append([row[2:] for row in cur_w[2:]])
   cur_w = mat_mul(W, cur_w)
  cur_s = mat_mul(S, cur_s)
 assert len(G12_6) == 12
 def in_D12(b):
  U = b[2] - b[4] - b[6]
  V = b[3] - b[5] - b[7]
  if U == 0 and V == 0:
   t = tuple(b[2:])
   if not any(t): return False
   return t == min(tuple(sum(M6[i][j] * t[j] for j in range(6)) for i in range(6)) for M6 in G12_6)
  if not (U > V and V >= 0): return False
  if b[2] != 0 or b[3] != 0: return (b[2] > 0) or (b[2] == 0 and b[3] > 0)
  d1x, d1y = b[4] - b[6], b[5] - b[7]
  if d1x != 0 or d1y != 0: return (d1x > 0) or (d1x == 0 and d1y > 0)
  return False
 R = [[0.0] * 8 for _ in range(8)]
 for i in range(8):
  for j in range(i + 1):
   s = sum(R[k][i] * R[k][j] for k in range(j))
   R[j][i] = math.sqrt(G[i][i] - s) if i == j else (G[i][j] - s) / R[j][j]
 pts = []
 cur = [0] * 8
 def fp_rec(k, rem):
  s = sum(R[k][j] * cur[j] for j in range(k + 1, 8))
  d = R[k][k]; rad = math.sqrt(max(0.0, rem)) / d; ctr = -s / d
  lo = math.ceil(ctr - rad - 1e-11); hi = math.floor(ctr + rad + 1e-11)
  if k == 0:
   for v in range(lo, hi + 1):
    cur[0] = v
    q = q_1w(cur)
    if 0 < q <= 30: pts.append(tuple(cur))
  else:
   for v in range(lo, hi + 1):
    cur[k] = v; diff = d * v + s; nrem = rem - diff * diff
    if nrem >= -1e-9: fp_rec(k - 1, nrem)
 fp_rec(7, 30.0)
 seen = set(); n_orb12 = exact1 = dup = missed = 0
 for b in pts:
  if b in seen: continue
  orb = set(); cs = list(b)
  for _ in range(4):
   cw = list(cs)
   for _ in range(3):
    orb.add(tuple(cw)); cw = mv(W, cw)
   cs = mv(S, cs)
  seen |= orb
  if len(orb) == 12:
   n_orb12 += 1
   h = sum(1 for x in orb if in_D12(x))
   if h == 1: exact1 += 1
   elif h > 1: dup += 1
   else: missed += 1
 assert dup == 0 and missed == 0 and len(pts) == 2940 and n_orb12 == 239 and exact1 == 239
 def _pows_1w(s, p):
  s3 = pow(s, 3, p); s5 = pow(s, 5, p); s6 = (s3 * s3) % p; s9 = (s6 * s3) % p
  return (1, s5, s3, (s3 * s5) % p, s6, (s6 * s5) % p, s9, (s9 * s5) % p)
 def eval_1w_pows(v, pw, p):
  return (v[0] * pw[0] + v[1] * pw[1] + v[2] * pw[2] + v[3] * pw[3] + v[4] * pw[4] + v[5] * pw[5] + v[6] * pw[6] + v[7] * pw[7]) % p
 rng = random.Random(20261119)
 R2_a, R2_b = 80.0, 70.0; SCALES = (1.00, 1.50, 2.00)
 a_orbits = []; seen_a = set()
 while len(a_orbits) < 40:
  v0 = [rng.randint(-2, 2) for _ in range(8)]
  if not any(v0): continue
  v1 = mv(S, v0); v2 = mv(S, v1); v3 = mv(S, v2)
  orb = [tuple(x) for x in (v0, v1, v2, v3)]
  if len(set(orb)) != 4 or any(x in seen_a for x in orb): continue
  qa = q_1w(v0) / R2_a
  a0_p = [sum(M_1W[i][j] * v0[j] for j in range(8)) for i in range(8)]
  _na0 = abs(norm_q15_z(a0_p)); _v2 = (_na0 & -_na0).bit_length() - 1 if _na0 else 8
  cls = 0 if _v2 == 0 else (1 if _v2 == 4 else 2)
  for _vk in (v1, v2, v3):
   _akp = [sum(M_1W[i][j] * _vk[j] for j in range(8)) for i in range(8)]
   _nak = abs(norm_q15_z(_akp)); _v2k = (_nak & -_nak).bit_length() - 1 if _nak else 8
   assert (0 if _v2k == 0 else (1 if _v2k == 4 else 2)) == cls
  if qa <= SCALES[cls]:
   seen_a.update(orb)
   a_orbits.append(([v0, v1, v2, v3], cls))
 b_rows = []; seen_t = set()
 while len(b_rows) < 35:
  tail = [0] + [rng.randint(-1, 1) for _ in range(7)]
  if not in_D12(tail): continue
  tt = tuple(tail)
  if tt in seen_t: continue
  pts_r = []
  for b0 in range(-6, 7):
   bv = list(tail); bv[0] = b0
   assert in_D12(bv)
   if q_1w(bv) <= R2_b: pts_r.append((b0, bv))
  if len(pts_r) >= 4: seen_t.add(tt); b_rows.append((tail, pts_r))
 def is_prime(n): return n >= 2 and all(n % d != 0 for d in range(2, math.isqrt(n) + 1))
 orbits4 = []
 for p in range(31, 350):
  if p % 15 != 1 or not is_prime(p): continue
  s_roots = [s for s in range(2, p) if pow(s, 15, p) == 1 and pow(s, 5, p) != 1 and pow(s, 3, p) != 1]
  seen_s = set()
  for s0 in s_roots:
   if s0 in seen_s: continue
   sk = [s0, pow(s0, 7, p), pow(s0, 4, p), pow(s0, 13, p)]
   seen_s.update(sk)
   w = pow(s0, 5, p); w2 = (w * w) % p
   rhs_g = (U_BLS * w2) % p
   for r in range(1, p):
    if pow(r, 3, p) == rhs_g: orbits4.append(("g", p, sk, r)); break
   c6 = (-w2) % p; c3 = ((-w - 1) * w) % p; c0 = (-w) % p
   for r in range(1, p):
    r3 = pow(r, 3, p)
    if (pow(r3, 3, p) + c6 * r3 * r3 + c3 * r3 + c0) % p == 0: orbits4.append(("f", p, sk, r)); break
  if len(orbits4) >= 12: break
 tot_hits = mis = 0; MASK = (1 << 64) - 1; OID_M = (1 << 56) - 1; BIT58 = 1 << 58; CLR58 = MASK ^ BIT58
 long_b0s = list(range(-80, 80))
 for idx_o, (side, p, sk, r) in enumerate(orbits4):
  tier = idx_o % 3
  s0 = sk[0]; r_inv = pow(r, -1, p); pw0 = _pows_1w(s0, p)
  if tier == 0:
   hsb = [0] * p; kbuf = []
   for m in range(4):
    for vs, _ in a_orbits:
     ka = (eval_1w_pows(vs[m], pw0, p) * r_inv) % p; hsb[ka] += 1; kbuf.append(ka)
   pos_buf = [0] * (p + 1); cur = 0
   for qi in range(p): pos_buf[qi + 1] = cur; cur += hsb[qi] * 8
   fba = bytearray(cur); idx_kb = 0
   for m in range(4):
    x_om = m << 56
    for _ in range(len(a_orbits)):
     ka = kbuf[idx_kb]; idx_kb += 1
     struct.pack_into("<Q", fba, pos_buf[ka + 1], x_om | ((ka & 31) << 59))
     pos_buf[ka + 1] += 8; x_om += 1
   off = pos_buf; flat = bytes(fba)
  elif tier == 1:
   n_sb = (p + 7) // 8; hsb = [0] * n_sb; kbuf = []
   for m in range(4):
    for vs, _ in a_orbits:
     ka = (eval_1w_pows(vs[m], pw0, p) * r_inv) % p; hsb[ka >> 3] += 1; kbuf.append(ka)
   off = [0]
   for c in hsb: off.append(off[-1] + c * 8)
   fba = bytearray(off[-1]); pos = list(off[:-1]); idx_kb = 0
   for m in range(4):
    x_om = m << 56
    for _ in range(len(a_orbits)):
     ka = kbuf[idx_kb]; idx_kb += 1; sbi = ka >> 3
     struct.pack_into("<Q", fba, pos[sbi], x_om | ((ka & 31) << 59)); pos[sbi] += 8; x_om += 1
   flat = bytes(fba)
  else:
   n_sb = (p + 31) // 32; slot = [MASK] * n_sb; ovf = {}; ovf_base = 0x3800000000000000; ovf_adj = (ovf_base - (1 << 62)) & MASK; nxt_tag = BIT58
   d_a0 = (4 * r_inv) % p
   for m in range(4):
    x_om = m << 56; prev_a = None; prev_c = None; ka = 0
    for vs, cls_a in a_orbits:
     av = vs[m]; ka_dir = (eval_1w_pows(av, pw0, p) * r_inv) % p
     if prev_a is not None and av[1:] == prev_a[1:] and av[0] == prev_a[0] + 4:
      ka = (ka + d_a0) % p; assert ka == ka_dir and cls_a == prev_c
     else:
      ka = ka_dir
     prev_a = av; prev_c = cls_a; sbi = ka >> 5
     rec = x_om | ((ka & 31) << 59); x_om += 1; cv = slot[sbi]
     if cv == MASK: slot[sbi] = rec
     else: ovf[ovf_base + 16 * (nxt_tag & (BIT58 - 1))] = (rec, cv); slot[sbi] = nxt_tag; nxt_tag += 1
  for k in range(4):
   s_k = sk[k]; pwk = _pows_1w(s_k, p); c_rot = (-k) & 3
   a_evals = [(eval_1w_pows(vs[0], pwk, p), eval_1w_pows(vs[1], pwk, p), eval_1w_pows(vs[2], pwk, p), eval_1w_pows(vs[3], pwk, p)) for (vs, cls) in a_orbits]
   test_rows = list(b_rows) + [(b_rows[0][0], [(b0, [b0] + b_rows[0][0][1:]) for b0 in long_b0s])]
   for tail, pts_r in test_rows:
    L_b = len(pts_r); k0 = eval_1w_pows(pts_r[0][1], pwk, p)
    slabs = [[{} for _ in range(32)] for _ in range(15)]
    tx_bufs = [[] for _ in range(15)]
    def _drain_tx(sid):
     for ptr in tx_bufs[sid]:
      cl = ptr >> 59; key = (ptr & OID_M, (ptr >> 56) & 3)
      slabs[sid][cl][key] = slabs[sid][cl].get(key, 0) + 1
     tx_bufs[sid].clear()
    def _emit_to(sid, ptr):
     assert -1 <= sid <= 15 and ((ptr >> 58) & 1) == 0
     if 0 <= sid < 15:
      tx_bufs[sid].append(ptr)
      if len(tx_bufs[sid]) >= 16: _drain_tx(sid)
    col0 = 0; rem_L = L_b; ks = k0
    while rem_L > 0:
     Lr = min(rem_L, p - ks); ke = ks + Lr
     if tier == 0:
      xt0 = ((((col0 - ks) & 31) << 59) | (c_rot << 56)) & MASK
      sid = col0 >> 5; oc = ks; oe = ke
      if (col0 & 31) != 0:
       on = min(oc + (32 - (col0 & 31)), oe); assert 0 <= oc <= on <= p
       for ofs in range(off[oc], off[on], 8):
        rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(sid, ((rec + xt0) & CLR58) & MASK)
       oc = on; sid += 1
      x21 = (oe - oc) >> 5
      for _ in range(x21):
       assert 0 <= oc and oc + 32 <= oe <= p
       for ofs in range(off[oc], off[oc + 32], 8):
        rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(sid, ((rec + xt0) & CLR58) & MASK)
       oc += 32; sid += 1
      if oc < oe:
       assert 0 <= oc < oe <= p
       for ofs in range(off[oc], off[oe], 8):
        rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(sid, ((rec + xt0) & CLR58) & MASK)
     elif tier == 1:
      delta = (col0 - ks) & 31; r_bnd = (-delta) & 31
      xt0 = ((delta << 59) | (c_rot << 56)) & MASK
      q = ks >> 3; q_end = (ke + 7) >> 3; assert 0 <= q < q_end <= len(off) - 1
      pref_n = ((r_bnd >> 3) - (q & 3)) & 3
      if pref_n > 0:
       s_cur = (col0 + 8 * q - ks) >> 5; qp = min(q + pref_n, q_end)
       for ofs in range(off[q], off[qp], 8):
        rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(s_cur, ((rec + xt0) & CLR58) & MASK)
       q = qp
      if q < q_end:
       if delta == 0:
        s_cur = (col0 + 8 * q - ks) >> 5; x21 = (q_end - q) >> 2
        for _ in range(x21):
         assert q + 4 <= q_end <= len(off) - 1
         for ofs in range(off[q], off[q + 4], 8):
          rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(s_cur, ((rec + xt0) & CLR58) & MASK)
         q += 4; s_cur += 1
        if q < q_end:
         for ofs in range(off[q], off[q_end], 8):
          rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(s_cur, ((rec + xt0) & CLR58) & MASK)
       else:
        s_lo = (col0 + 8 * q - ks - 1) >> 5; x21 = (q_end - q) >> 2
        for _ in range(x21):
         assert q + 4 <= q_end <= len(off) - 1
         for ofs in range(off[q], off[q + 1], 8):
          rec, = struct.unpack_from("<Q", flat, ofs); s_raw = rec + xt0
          _emit_to(s_lo + (1 if s_raw > MASK else 0), (s_raw & CLR58) & MASK)
         for ofs in range(off[q + 1], off[q + 4], 8):
          rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(s_lo + 1, ((rec + xt0) & CLR58) & MASK)
         q += 4; s_lo += 1
        if q < q_end:
         assert q + 1 <= q_end <= len(off) - 1
         for ofs in range(off[q], off[q + 1], 8):
          rec, = struct.unpack_from("<Q", flat, ofs); s_raw = rec + xt0
          _emit_to(s_lo + (1 if s_raw > MASK else 0), (s_raw & CLR58) & MASK)
         for ofs in range(off[q + 1], off[q_end], 8):
          rec, = struct.unpack_from("<Q", flat, ofs); _emit_to(s_lo + 1, ((rec + xt0) & CLR58) & MASK)
     else:
      cb0 = col0 - (ks & 31); s_base = cb0 >> 5; delta = cb0 & 31
      xt0 = ((delta << 59) | (c_rot << 56)) & MASK
      x4 = ks >> 5; x6 = (ke - 1) >> 5; x21 = (x6 - x4) & MASK
      cf_w = (x21 >= 1); x21 = (x21 - 1) & MASK
      si = x4; st = 0
      def _emit_bkt(bi, s_off):
       cur_v = slot[bi]
       if cur_v == MASK: return
       while cur_v & BIT58:
        val, cur_v = ovf[(ovf_adj + ((cur_v << 4) & MASK)) & MASK]
        s_raw = val + xt0; _emit_to(s_base + s_off + (1 if s_raw > MASK else 0), (s_raw & CLR58) & MASK)
       s_raw = cur_v + xt0; _emit_to(s_base + s_off + (1 if s_raw > MASK else 0), (s_raw & CLR58) & MASK)
      if cf_w:
       while True:
        assert si + 1 < len(slot)
        x1_s, x2_s = slot[si], slot[si + 1]
        if (x1_s & x2_s) != MASK:
         _emit_bkt(si, st); _emit_bkt(si + 1, st + 1)
        si += 2; st += 2
        cf_l = (x21 >= 2); x21 = (x21 - 2) & MASK
        if not cf_l: break
      x21 = (x21 + 2) & MASK
      if x21 != 0:
       assert x21 == 1 and si < len(slot)
       _emit_bkt(si, st)
     col0 += Lr; rem_L -= Lr; ks = 0
    for _sid in range(15): _drain_tx(_sid)
    rep_hits = {}
    for t in range(L_b):
     for key, cnt in slabs[t >> 5][t & 31].items(): rep_hits[(key[0], key[1], t)] = cnt
    if tier == 0:
     for t in range(L_b, 15 * 32): assert not slabs[t >> 5][t & 31]
    gt_hits = {}
    for t, (b0, bv) in enumerate(pts_r):
     rbs = (r * eval_1w_pows(bv, pwk, p)) % p
     for oid, ev4 in enumerate(a_evals):
      for mm in range(4):
       if ev4[mm] == rbs: kk = (oid, mm, t); gt_hits[kk] = gt_hits.get(kk, 0) + 1
    if rep_hits != gt_hits: mis += len(set(rep_hits.items()) ^ set(gt_hits.items()))
    tot_hits += sum(gt_hits.values())
 assert mis == 0 and tot_hits > 5000
 return dict(d12_pts=len(pts), d12_orb12=n_orb12, d12_exact1=exact1, d12_dups=dup, d12_missed=missed, tables=len(orbits4), ideals=4 * len(orbits4), n_slabs=15, long_chord_cols=len(long_b0s), gt_hits=tot_hits, arch="Arch3-K2-Lq75-M15-W1W8Dir32-Bit58-U8Fused", mismatches=mis)
def verify_v06_kernels_extra(parsed, coll, la, r_bls, ntt3_primes):
 import math, random
 rng = random.Random(20261124)
 M64 = (1 << 64) - 1; OID_M = (1 << 56) - 1; BIT58 = 1 << 58; CLR58 = M64 ^ BIT58
 pb = 509; nsb1 = ((pb + 15) // 16) * 16; db = 137; pmsb = pb - db; xk = 41
 for _pt in ((1 << 64) - 59, (1 << 64) - 1, (1 << 63) + 1, pb):
  for _aa in (0, 1, _pt // 2, _pt - 2, _pt - 1):
   for _bb in (0, 1, _pt // 2, _pt - 2, _pt - 1, (_pt - _aa) % _pt):
    _x3 = (_pt - _bb) & M64; _x1 = (_aa - _x3) & M64; _x2 = (_aa + _bb) & M64
    _res = _x1 if _aa >= _x3 else _x2
    assert _res == (_aa + _bb) % _pt
 _b1a_cnt = 4000000
 assert 8 * _b1a_cnt < (1 << 64) and 8 * int(0.9351 * (1 << 58)) < (1 << 61) < (1 << 64)
 fbase = 0x10000000; pos_buf = [0] * (nsb1 + 1); pos_buf[0] = fbase
 kbuf1 = []; ref1 = [[] for _ in range(nsb1)]
 for _ in range(256):
  for _j in range(8):
   xk = (xk - pmsb) if xk >= pmsb else (xk + db)
   kbuf1.append(xk); pos_buf[xk + 1] += 8
 xcur = fbase
 for s in range(0, nsb1, 8):
  for j in range(8):
   _cb = pos_buf[s + j + 1]; pos_buf[s + j + 1] = xcur; xcur += _cb
 flat1 = {}; xom = 0
 for kv in kbuf1:
  rec = ((kv & 31) << 59) | xom; oid = xom & OID_M; mb = (xom >> 56) & 3; xom += 1
  ptr = pos_buf[kv + 1]; flat1[ptr] = rec; pos_buf[kv + 1] = ptr + 8
  ref1[kv].append(((kv & 31) << 59) | (mb << 56) | oid)
 for s in range(nsb1): assert [flat1[a] for a in range(pos_buf[s], pos_buf[s + 1], 8)] == ref1[s]
 nsb8 = (((pb + 7) // 8) + 15) // 16 * 16; pos8 = [0] * (nsb8 + 1); pos8[0] = fbase; ref8 = [[] for _ in range(nsb8)]
 for kv in kbuf1: pos8[(kv >> 3) + 1] += 8
 xcur = fbase
 for s in range(0, nsb8, 8):
  for j in range(8):
   _cb = pos8[s + j + 1]; pos8[s + j + 1] = xcur; xcur += _cb
 flat8 = {}; xom = 0
 for kv in kbuf1:
  s = kv >> 3; rec = ((kv & 31) << 59) | xom; oid = xom & OID_M; mb = (xom >> 56) & 3; xom += 1
  ptr = pos8[s + 1]; flat8[ptr] = rec; pos8[s + 1] = ptr + 8
  ref8[s].append(((kv & 31) << 59) | (mb << 56) | oid)
 for s in range(nsb8): assert [flat8[a] for a in range(pos8[s], pos8[s + 1], 8)] == ref8[s]
 B8 = (pb + 7) // 8
 for ks_w8, col0_w8 in ((pb - 1, 0), (pb - 5, 0), (pb - 9, 0), (pb - 32, 0), (0, 0), (0, 5), (0, 8), (0, 9), (8, 8)):
  ke_w8 = min(pb, ks_w8 + 160); q_w = ks_w8 >> 3; qe_w = (ke_w8 + 7) >> 3; assert 0 <= q_w < qe_w <= B8
  delta_w = (col0_w8 - ks_w8) & 31; r_bnd_w = (-delta_w) & 31; pref_w = ((r_bnd_w >> 3) - (q_w & 3)) & 3
  q_after_pre = min(q_w + pref_w, qe_w); x21_w = (qe_w - q_after_pre) >> 2; q_after = q_after_pre + 4 * x21_w
  assert q_after <= qe_w <= B8 and (q_after == qe_w or q_after + 1 <= qe_w <= B8)
  if ks_w8 == 0 and col0_w8 == 9:
   assert pref_w == 2 and (pos8[2] - pos8[0]) // 8 >= 4 and q_after_pre < qe_w and ((col0_w8 + 8 * q_after_pre - ks_w8 - 1) >> 5) == 0
 nsb32 = (pb + 31) // 32; slot = [M64] * nsb32; ovf_nodes = {}; ovf_base = 0x3800000000000000; ovf_adj = (ovf_base - (1 << 62)) & M64; nxt_tag = BIT58; dir_ref = [[] for _ in range(nsb32)]
 for oid in range(600):
  kv = 32 if oid == 0 else rng.randrange(pb); s = kv >> 5; m = rng.randrange(4)
  rec = ((kv & 31) << 59) | (m << 56) | oid; dir_ref[s].append(rec)
  if oid == 0: assert (rec >> 59) == 0
  cur = slot[s]
  if cur == M64: slot[s] = rec
  else: ovf_nodes[ovf_base + 16 * (nxt_tag & (BIT58 - 1))] = (rec, cur); slot[s] = nxt_tag; nxt_tag += 1
 for _ob in (0, 0x0123456789abcde0, 0x3800000000000000, 0x8000000000000000):
  _oadj = (_ob - (1 << 62)) & M64
  for _ni in (0, 1, (1 << 54) - 1, 1 << 54, (1 << 54) + 7, 1 << 55, 1 << 56, int(0.9351 * (1 << 58)) - 1, (1 << 58) - 1):
   _tg = BIT58 | _ni
   assert (_tg & BIT58) != 0 and ((_oadj + ((_tg << 4) & M64)) & M64) == ((_ob + 16 * _ni) & M64)
 for s in range(0, nsb32, 2):
  x1 = slot[s]; x2 = slot[s + 1] if (s + 1) < nsb32 else M64
  if (x1 & x2) == M64:
   assert not dir_ref[s] and ((s + 1) >= nsb32 or not dir_ref[s + 1]); continue
  for bi, cur in ((s, x1), (s + 1, x2)):
   if bi >= nsb32: continue
   got = []
   if cur != M64:
    while cur & BIT58: val, cur = ovf_nodes[(ovf_adj + ((cur << 4) & M64)) & M64]; got.append(val)
    got.append(cur)
   assert sorted(got) == sorted(dir_ref[bi])
 for nb_test in range(1, 8):
  b0 = nsb32 - nb_test; ks_e = b0 * 32 + (pb - 1 - b0 * 32 if nb_test == 1 else 0); ke_e = pb
  x4 = ks_e >> 5; x6 = (ke_e - 1) >> 5; assert x4 == b0 and x6 == nsb32 - 1
  x21 = (x6 - x4) & M64; cf_w = (x21 >= 1); x21 = (x21 - 1) & M64
  visited = []; si = x4
  if cf_w:
   while True:
    assert si + 1 < nsb32; visited.extend((si, si + 1)); si += 2
    cf_l = (x21 >= 2); x21 = (x21 - 2) & M64
    if not cf_l: break
  x21 = (x21 + 2) & M64
  if x21 != 0: assert x21 == 1 and si < nsb32; visited.append(si)
  assert visited == list(range(b0, nsb32))
 for k0 in range(9):
  for k1 in range(9):
   if k0 == 0 and k1 == 0: continue
   n_ne = (1 if k0 > 0 else 0) + (1 if k1 > 0 else 0); m_ovf = max(0, k0 - 1) + max(0, k1 - 1)
   ops = (2 if k0 == 0 else (10 if k0 == 1 else 13 + 9 * (k0 - 1))) + (2 if k1 == 0 else (11 if k1 == 1 else 13 + 9 * (k1 - 1)))
   assert ops <= 13.0 * n_ne + 11.0 * m_ovf + 1e-12
 _mx_b_t = 0.0
 for n in range(1, 201):
  q8, r8 = divmod(n, 8)
  for rt_b, b8_ops, k4, k2, k1 in ((5.90, 47, 22, 11, 6), (6.90, 55, 26, 13, 7), (11.25, 90, 44, 22, 11)):
   op8 = 4 + q8 * b8_ops + ((k4 + 1 if (r8 & 4) else 1) + (k2 + 1 if (r8 & 2) else 1) + (k1 + 1 if (r8 & 1) else 1) if r8 else 0)
   _mx_b_t = max(_mx_b_t, op8 - rt_b * n)
 assert 12.0 + _mx_b_t <= 19.25 + 1e-12
 # Exhaustive exact 2-adic class invariance: a -> a + 2 preserves cls==0 (2 in P1*P2); a -> a + 4 preserves ALL cls in {0,1,2} (4 in P1^2*P2^2):
 for _c in range(256):
  _v = [(_c >> i) & 1 for i in range(8)]
  _v2 = [_v[0] + 2] + _v[1:]; _v4 = [_v[0] + 4] + _v[1:]
  assert all((_v2[i] - _v[i]) % 2 == 0 for i in range(8)) and all((_v4[i] - _v[i]) % 4 == 0 for i in range(8))
 # Adversarial TX-buffer occupancy stress test across W1, W8 (split+interior), and Direct1 handoff (starting at 128 KiB - 8 B):
 _B17 = 1 << 17; _B18 = 1 << 18; _tx_o = [131064, 131064]; _max_o = 0
 for _n_adv in (20000, 16385, 16384, 31, 15, 1):
  if _tx_o[0] & _B17: _tx_o[0] = 0
  for _ in range(_n_adv // 16):
   _tx_o[0] += 128; _max_o = max(_max_o, _tx_o[0])
   if _tx_o[0] & _B17: _tx_o[0] = 0
  _tx_o[0] += 8 * (_n_adv % 16); _max_o = max(_max_o, _tx_o[0])
  assert _tx_o[0] < _B18
 for _n_sp, _n_in in ((20000, 20000), (16385, 16385), (7, 15)):
  if _tx_o[0] & _B17: _tx_o[0] = 0
  if _tx_o[1] & _B17: _tx_o[1] = 0
  for _ in range(_n_sp // 4):
   for _j in range(4): _tx_o[_j & 1] += 8; _max_o = max(_max_o, _tx_o[0], _tx_o[1])
   if _tx_o[0] & _B17: _tx_o[0] = 0
   if _tx_o[1] & _B17: _tx_o[1] = 0
  for _j in range(_n_sp % 4): _tx_o[0] += 8; _max_o = max(_max_o, _tx_o[0])
  for _ in range(_n_in // 16):
   _tx_o[0] += 128; _max_o = max(_max_o, _tx_o[0])
   if _tx_o[0] & _B17: _tx_o[0] = 0
  _tx_o[0] += 8 * (_n_in % 16); _max_o = max(_max_o, _tx_o[0])
  assert _tx_o[0] < _B18 and _tx_o[1] < _B18
  for _d in range(32):
   _tx_o[0] += 8; _max_o = max(_max_o, _tx_o[0])
   if _tx_o[0] & _B17: _tx_o[0] = 0
 assert _max_o <= 131216 < _B18
 mx_u16 = mx_spl = -1e9
 for n in range(1, 201):
  q16, r16 = divmod(n, 16)
  for base16, b16_ops, k8, k4, k2, k1 in ((35.0/16.0, 35, 17, 9, 5, 4), (51.0/16.0, 51, 25, 13, 7, 5)):
   op16 = 4 + q16 * b16_ops + ((k8 if (r16 & 8) else 1) + (k4 if (r16 & 4) else 1) + (k2 if (r16 & 2) else 1) + (k1 if (r16 & 1) else 1) if r16 else 0)
   mx_u16 = max(mx_u16, op16 - base16 * n)
  q4, r4 = divmod(n, 4)
  for base4, b4_ops, k2, k1 in ((5.5, 22, 10, 6), (6.5, 26, 12, 7)):
   op4 = 4 + q4 * b4_ops + ((k2 if (r4 & 2) else 1) + (k1 if (r4 & 1) else 1) if r4 else 0)
   mx_spl = max(mx_spl, op4 - base4 * n)
 assert 11.0 + mx_u16 <= float(coll["c_sb_hdr_u4"]) + 1e-12
 assert 16.0 + mx_u16 + mx_spl <= float(coll["c_sb_hdr_u2"]) + 1e-12
 mx1al = -1e9
 for kl in range(256):
  oids = [rng.randrange(64) for _ in range(kl)]; ref_a = [0] * 64; wlp = 13
  for o in oids: ref_a[o] += wlp
  si = so = 0; sd = []; xc = 8 * kl
  if xc:
   so += 2; cf = (xc >= 512); xc = (xc - 512) & M64
   if cf:
    while True:
     sd.extend(oids[si:si+64]); si += 64; so += 66; cf = (xc >= 512); xc = (xc - 512) & M64
     if not cf: break
   xc = (xc + 512) & M64; so += 2
   if xc:
    for sh_b, cnt_b, op_b in ((256, 32, 32), (128, 16, 16), (64, 8, 8), (32, 4, 4), (16, 2, 2), (8, 1, 2)):
     so += 1
     if xc & sh_b: sd.extend(oids[si:si+cnt_b]); si += cnt_b; so += op_b
  assert si == kl and sd == oids
  ai = ao = 0; ac = [0] * 64; xc = 8 * kl
  if xc:
   ao += 2; cf = (xc >= 512); xc = (xc - 512) & M64
   if cf:
    while True:
     for t in range(64): ac[sd[ai+t]] += wlp
     ai += 64; ao += 98; cf = (xc >= 512); xc = (xc - 512) & M64
     if not cf: break
   xc = (xc + 512) & M64; ao += 2
   if xc:
    for sh_b, cnt_b, op_b in ((256, 32, 48), (128, 16, 24), (64, 8, 12), (32, 4, 6), (16, 2, 3), (8, 1, 2)):
     ao += 1
     if xc & sh_b:
      for t in range(cnt_b): ac[sd[ai+t]] += wlp
      ai += cnt_b; ao += op_b
  assert ai == kl and ac == ref_a
  mx1al = max(mx1al, (so + ao) - 2.5625 * kl)
 assert 6.0 + 1.0 * (4.0 + mx1al) <= float(coll["b1a_long_lk"])
 for p in ntt3_primes:
  assert 8 * p < (1 << 64)
  g = 2
  while pow(g, (p - 1) // 2, p) == 1: g += 1
  jv = pow(g, (p - 1) // 4, p); jp = (jv << 64) // p; p2v = 2 * p; p4v = 4 * p
  for t in (8 * p - 1, 7 * p, 4 * p - 1, rng.randrange(8 * p)):
   r = (t * jv - ((t * jp) >> 64) * p) & M64; assert r < p2v and r % p == (t * jv) % p
  for _ in range(32):
   a0, a1, a2, a3 = [rng.randrange(8 * p) for _ in range(4)]
   w1 = rng.randrange(p); w2 = (w1 * w1) % p; w3 = (w2 * w1) % p
   x1 = (a1 * w1 - ((a1 * ((w1 << 64) // p)) >> 64) * p) & M64
   x2 = (a2 * w2 - ((a2 * ((w2 << 64) // p)) >> 64) * p) & M64
   x3 = (a3 * w3 - ((a3 * ((w3 << 64) // p)) >> 64) * p) & M64
   x0 = a0 - p4v if a0 >= p4v else a0; x0 = x0 - p2v if x0 >= p2v else x0
   assert 0 <= x0 < p2v and 0 <= x1 < p2v and 0 <= x2 < p2v and 0 <= x3 < p2v
   x8 = x0 + x2; x10 = x0 + p2v - x2; x11 = x1 + x3; x12 = x1 + p2v - x3
   assert 0 <= x8 < p4v and 0 <= x10 < p4v and 0 <= x11 < p4v and 0 <= x12 < p4v
   x12 = (x12 * jv - ((x12 * jp) >> 64) * p) & M64; assert x12 < p2v
   y0 = x8 + x11; y2 = x8 + p4v - x11; y1 = x10 + x12; y3 = x10 + p2v - x12
   assert all(0 <= y < 8 * p for y in (y0, y1, y2, y3))
   assert (y0 % p, y1 % p, y2 % p, y3 % p) == (
   (a0 + a1*w1 + a2*w2 + a3*w3) % p, (a0 + jv*a1*w1 - a2*w2 - jv*a3*w3) % p,
   (a0 - a1*w1 + a2*w2 - a3*w3) % p, (a0 - jv*a1*w1 - a2*w2 + jv*a3*w3) % p)
 def _leh_gcd(U, V):
  if U < V: U, V = V, U
  ops = 0
  while V > (1 << 64):
   bl = U.bit_length(); L = (bl + 63) // 64; sh = bl - 64
   x0, x1, x2, x3, x4, x5 = U >> sh, V >> sh, 1, 0, 0, 1; ops += 32
   while x1:
    q = x0 // x1; r = x0 - q * x1; nx2 = x2 + q * x4; nx3 = x3 + q * x5; ops += 12
    if r < nx3 or (x1 - r) < (x4 + nx2) or r == 0: break
    q2 = x1 // r; x14 = x1 - q2 * r; x15 = x4 + q2 * nx2; x16 = x5 + q2 * nx3; ops += 17
    if x14 < x15 or (r - x14) < (nx3 + x16): break
    x0, x1, x2, x3, x4, x5 = r, x14, nx2, nx3, x15, x16
   if x3 == 0: U, V = V, U % V; ops += 32 * L
   else:
    Un = x2 * U - x3 * V; Vn = x5 * V - x4 * U; assert Un >= 0 and Vn >= 0
    U, V = (Un, Vn) if Un >= Vn else (Vn, Un); ops += 64 * L
  while V: U, V = V, U % V; ops += 16
  return U, ops
 nb = 680; L_g = math.ceil(nb / 64); bnd_leh = (nb / 24.0) * (64.0 * L_g + 544.0); max_lo = 0
 for _ in range(48):
  a = rng.getrandbits(nb) | (1 << (nb - 1)); b = rng.getrandbits(nb) | (1 << (nb - 1))
  g, o = _leh_gcd(a, b); assert g == math.gcd(a, b); max_lo = max(max_lo, o)
 f0, f1 = 1, 1
 while f1.bit_length() < nb: f0, f1 = f1, f0 + f1
 gf, of = _leh_gcd(f1, f0); assert gf == 1 and max(max_lo, of) < bnd_leh
 mu_r = (1 << 318) // r_bls
 for Av in [0, 1, -1, r_bls, -r_bls, 2*r_bls - 1, -(2*r_bls - 1), (1 << 290), -(1 << 290), 4096*255*(r_bls - 1)] + [rng.randint(-(1 << 290), 1 << 290) for _ in range(64)]:
  U = abs(Av); top = (U >> 254) & M64; qh = (top * mu_r) >> 64; rem = U - qh * r_bls
  assert 0 <= rem < 2 * r_bls
  if rem >= r_bls: rem -= r_bls
  if Av < 0 and rem != 0: rem = r_bls - rem
  assert rem == Av % r_bls
 return dict(r1_build_u16_ok=True, direct1_sparse_b2_ok=True, b1a_long_lk_sim=round(6.0 + 1.0 * (4.0 + mx1al), 2), ntt61_lazy8p_ok=True, lehmer_max_ops=max(max_lo, of), lehmer_bound=round(bnd_leh, 1), row_red5_57_ok=True)

LN2 = math.log(2.0)
def missing_recon_z12(log2U=14,seed=1,lf=71.56,lg=71.56,nf=545.0,ng=665.0,excess=1.02,p_f=60.0,p_g=62.6,ret_f=0.917,ret_g=0.82,sbf=362.0,sbg=315.0):
 import bisect as _bs,random as _rd,filtsim_z12 as _FS
 rg=_rd.Random(seed);Us=1<<log2U;uf,ug=nf/lf,ng/lg
 LiBf=(2.0**lf)/(lf*LN2);LiBg=(2.0**lg)/(lg*LN2);LiBs=LiBf+LiBg
 Upf=int(round(Us*LiBf/LiBs));Upg=Us-Upf;nr=int(round(excess*Us))
 rf0,rg0=_FS.rho_exact(uf),_FS.rho_exact(ug)
 kf=lambda t:_FS.rho_exact(uf-t)/rf0 if rf0>0 else 0.0
 kg=lambda t:_FS.rho_exact(ug-t)/rg0 if rg0>0 else 0.0
 nt=1000;ts=[0.80+0.20*(i+0.5)/nt for i in range(nt)];cf,cg=[],[];tf=tg=0.0
 for t in ts:tf+=(2.0**(lf*(t-1.0)))/t;tg+=(2.0**(lg*(t-1.0)))/t;cf.append(tf);cg.append(tg)
 cf=[x/tf for x in cf];cg=[x/tg for x in cg];cl_raw=[];cl=[];c_meta=[];ds=0.20/nt
 for side_id,cnt,cm,L,kp,ps,rm in ((0,Upf,cf,lf,kf,p_f,ret_f),(1,Upg,cg,lg,kg,p_g,ret_g)):
  for _ in range(cnt):
   ix=min(nt-1,_bs.bisect_left(cm,rg.random()));t=ts[ix]+(rg.random()-0.5)*ds
   lb=L*t;base_l=excess*LiBs*kp(t)/(2.0**lb)
   cl_raw.append(min(float(nr),base_l))
   cl.append(min(float(nr),base_l*(rm if lb>ps else 1.0)))
   c_meta.append((side_id,lb if lb>ps else 0.0))
 nc=len(cl);p_keep_est=0.688
 while True:
  nr_raw=max(nr,int(math.ceil(nr/p_keep_est)));sc_raw=nr_raw/nr
  row_bf=[0.0]*nr_raw;row_bg=[0.0]*nr_raw;uns_buf=[]
  for c in range(nc):
   sd_c,lb_c=c_meta[c]
   if lb_c<=0.0:continue
   k=min(nr_raw,_FS._poisson(rg,min(0.90*nr_raw,cl_raw[c]*sc_raw)))
   if not k:continue
   rb_arr=row_bf if sd_c==0 else row_bg
   if k==1:
    r=rg._randbelow(nr_raw);rb_arr[r]+=lb_c;uns_buf.append((r,c))
   else:
    for r in rg.sample(range(nr_raw),k):rb_arr[r]+=lb_c;uns_buf.append((r,c))
  surv=[r for r in range(nr_raw)if row_bf[r]<=sbf and row_bg[r]<=sbg]
  if len(surv)>=nr:
   raw_to_surv=[-1]*nr_raw
   for idx_s,r in enumerate(surv[:nr]):raw_to_surv[r]=idx_s
   assert all(row_bf[r]<=sbf+1e-9 and row_bg[r]<=sbg+1e-9 for r in surv[:nr]),"Stage-B residual cap violated in missing_recon_z12!"
   break
  p_keep_est*=0.5;assert p_keep_est>1e-6
 rc=[[]for _ in range(nr)];cr=[[]for _ in range(nc)];rd=[0]*nc
 for r,c in uns_buf:
  sr=raw_to_surv[r]
  if sr>=0:rc[sr].append(c);cr[c].append(sr);rd[c]+=1
 for c in range(nc):
  if c_meta[c][1]>0.0:continue
  k=min(nr,_FS._poisson(rg,cl_raw[c]))
  if not k:continue
  if k==1:
   sr=rg._randbelow(nr);rc[sr].append(c);cr[c].append(sr);rd[c]+=1
  else:
   for sr in rg.sample(range(nr),k):rc[sr].append(c);cr[c].append(sr);rd[c]+=1
 dg=list(rd);al=[True]*nr;st=[];qq=[c for c in range(nc)if dg[c]==1]
 while qq:
  c=qq.pop()
  if dg[c]!=1:continue
  r=next((rr for rr in cr[c]if al[rr]),-1)
  if r<0:continue
  al[r]=False;st.append((r,c))
  for c2 in rc[r]:
   dg[c2]-=1
   if dg[c2]==1:qq.append(c2)
 kn=[dg[c]>=2 for c in range(nc)]
 for r,c in reversed(st):
  if all(kn[c2]or c2==c for c2 in rc[r]):kn[c]=True
 out_b=[]
 for lo,hi in ((1.5,3.0),(3.0,5.0),(5.0,8.0),(8.0,14.0),(14.0,30.0)):
  cs=[c for c,lm in enumerate(cl)if lo<=lm<hi]
  if not cs:continue
  n_=len(cs);mf=sum(1 for c in cs if not kn[c])/n_
  md=sum(max(math.exp(-cl[c])*(1.0+cl[c]+0.5*cl[c]**2+(cl[c]**3)/6.0),math.exp(-0.75*cl[c])*(1.0+0.75*cl[c]))for c in cs)/n_
  l3=sum(1 for c in cs if(not kn[c])and rd[c]>=3)
  assert mf<=md,f"missing_recon band [{lo},{hi}): meas {mf} > mod {md}"
  out_b.append(dict(band=[lo,hi],n=n_,meas_miss=round(mf,5),mod_miss=round(md,5),model_over_meas=round(md/mf,2)if mf>0 else 999.0,lost_ge3=l3))
 return out_b
