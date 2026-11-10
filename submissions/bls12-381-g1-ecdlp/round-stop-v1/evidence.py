#!/usr/bin/env python3 -I
"""Paired class-walk evidence: ledger/planner/ec/proto/census. Stdlib only."""
import json, math, random, sys, time
from fractions import Fraction
R_PRODUCTION = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D_PRODUCTION = 1 << 20
M_PRODUCTION = (R_PRODUCTION - 1) // D_PRODUCTION
def is_prime(n):
 if n < 2:
  return False
 for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
  if n % q == 0:
   return n == q
 e, s = n - 1, 0
 while e % 2 == 0:
  e //= 2
  s += 1
 for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
  x = pow(a, e, n)
  if x in (1, n - 1):
   continue
  for _ in range(s - 1):
   x = x * x % n
   if x == n - 1:
    break
  else:
   return False
 return True
def primes_in(lo, hi):
 return [q for q in range(lo, hi + 1) if is_prime(q)]
def sqrt_mod(a, p):
 a %= p
 if a == 0:
  return 0
 if pow(a, (p - 1) // 2, p) != 1:
  return None
 q, s = p - 1, 0
 while q % 2 == 0:
  q //= 2
  s += 1
 z = next(z for z in range(2, p) if pow(z, (p - 1) // 2, p) == p - 1)
 mm, c, t, x = s, pow(z, q, p), pow(a, q, p), pow(a, (q + 1) // 2, p)
 while t != 1:
  i, t2 = 0, t
  while t2 != 1:
   t2 = t2 * t2 % p
   i += 1
  b = pow(c, 1 << (mm - i - 1), p)
  mm, c, t, x = i, b * b % p, t * b * b % p, x * b % p
 return x
def toy_prime(bits, dd, rng):
 while True:
  rr = rng.randrange(1 << (bits - 1), 1 << bits) // (6 * dd) * (6 * dd) + 1
  if is_prime(rr) and ((rr - 1) // dd) % 6 == 0 and pow(2, (rr - 1) // dd, rr) != 1:
   return rr
def character(rr, dd):
 """chi: integers -> Z/d."""
 mm = (rr - 1) // dd
 gen = next(g for g in range(2, 10 ** 6)
         if all(pow(g, (rr - 1) // q, rr) != 1 for q in (2, 3)))
 zeta = pow(gen, mm, rr)
 zlog = {pow(zeta, i, rr): i for i in range(dd)}
 return (lambda x: zlog[pow(x, mm, rr)]), gen
def class_lists(rr, dd, X2):
 chi, _ = character(rr, dd); f, c = [0] * (X2 + 1), [0] * (X2 + 1)
 for q in range(2, math.isqrt(X2) + 1): f[q * q::q] = [q] * ((X2 - q * q) // q + 1)
 for b in range(2, X2 + 1): c[b] = (c[f[b]] + c[b // f[b]]) % dd if f[b] else chi(b)
 assert all(c[b] == chi(b) for b in range(1, X2 + 1, 97))
 classes = [[] for _ in range(dd)]
 for b in range(1, X2 + 1): classes[c[b]].append(b)
 return classes
def h2_sizes(seed=3):
 """s35: toy class-length spread."""
 rng = random.Random(seed); out = []
 for lg, lx in ((4, 14), (4, 18), (6, 18), (8, 18), (8, 20)):
  dd = 1 << lg; L = [len(c) for c in class_lists(toy_prime(40, dd, rng), dd, 1 << lx)]
  out.append([dd, lx - lg, round(min(L) * dd / 2 ** lx, 4), round(max(L) * dd / 2 ** lx, 4)])
 return out
from planner import REL32_NONE, REL32_MAX, rel32, plan, OPS_BOUND, CC, CCY
from smx_audit import SMX, audit_smx
import batch_audit, kernel_check, smx_audit, re
from batch_audit import BATCH, ROUND_FIXED, RMAXC, RC_FIX, audit_round_ctl, eps_lit
class Curve:
 def __init__(self, p, c):
  self.p, self.c = p, c
  self.inv = self.muls = 0
 def add(self, P, Q):
  if P is None:
   return Q
  if Q is None:
   return P
  p = self.p
  (x1, y1), (x2, y2) = P, Q
  if x1 == x2:
   if (y1 + y2) % p == 0:
    return None
   lam = 3 * x1 * x1 * pow(2 * y1, -1, p) % p
  else:
   lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
  x3 = (lam * lam - x1 - x2) % p
  return (x3, (lam * (x1 - x3) - y1) % p)
 def single(self, P, T):
  """Root step P+T."""
  assert P is not None and T is not None and P[0] != T[0]
  self.inv += 1
  self.muls += 6
  return self.add(P, T)
 def pair(self, P, T):
  """P+T,P-T from one inversion."""
  p = self.p
  (x1, y1), (x2, y2) = P, T
  assert x1 != x2
  inv = pow(x2 - x1, -1, p)
  self.inv += 1
  self.muls += 9
  lp = (y2 - y1) * inv % p
  lm = (-y2 - y1) * inv % p
  x3p = (lp * lp - x1 - x2) % p
  x3m = (lm * lm - x1 - x2) % p
  return ((x3p, (lp * (x1 - x3p) - y1) % p),
    (x3m, (lm * (x1 - x3m) - y1) % p))
 def mul(self, k, P):
  R = None
  while k:
   if k & 1:
    R = self.add(R, P)
   P = self.add(P, P)
   k >>= 1
  return R
def cm_ok(p, dd):
 """Speed only: skip p if no j=0 order p+1-t can pass."""
 x, y = p, sqrt_mod(-3, p)
 while y * y > p: x, y = y, x % y
 b = math.isqrt((p - y * y) // 3)
 return y * y + 3 * b * b != p or any(is_prime(N) and (N - 1) % (6 * dd) == 0 and pow(2, (N - 1) // dd, N) != 1 for t in (2 * y, y + 3 * b, y - 3 * b) for N in (p + 1 - t, p + 1 + t))
def find_curve(bits, dd, rng):
 """Toy curve: 6d | r-1, 2 not in H."""
 while True:
  p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
  if p % 3 != 1 or not is_prime(p):
   continue
  ok = cm_ok(p, dd)
  for c in range(1, 8):
   E = Curve(p, c)
   while True:
    x = rng.randrange(p)
    y = sqrt_mod((x ** 3 + c) % p, p)
    if y is not None:
     P = (x, y)
     break
   if not ok:
    continue
   lo, hi = p + 1 - 2 * math.isqrt(p) - 2, p + 1 + 2 * math.isqrt(p) + 2
   w = math.isqrt(hi - lo) + 1
   baby, R = {}, None
   for j in range(w):
    if R is not None:
     baby.setdefault(R, j)
    R = E.add(R, P)
   step = E.mul(w, P)
   step = (step[0], (-step[1]) % p) if step else None
   Gk, cands = E.mul(lo, P), []
   for i in range(w + 1):
    if Gk is None:
     cands.append(lo + i * w)
    elif (Gk[0], (-Gk[1]) % p) in baby:
     cands.append(lo + i * w + baby[(Gk[0], (-Gk[1]) % p)])
    Gk = E.add(Gk, step)
   cands = sorted({N for N in cands if lo <= N <= hi and E.mul(N, P) is None})
   if len(cands) == 1 and is_prime(cands[0]) and (cands[0] - 1) % (6 * dd) == 0 \
     and pow(2, (cands[0] - 1) // dd, cands[0]) != 1:
    return E, cands[0], P
def dth_root(t, dd, rr):
 k = dd.bit_length() - 1
 for j in range(k, 0, -1):
  y = sqrt_mod(t, rr)
  e = (rr - 1) >> (j - 1)
  t = next(z for z in (y, rr - y) if j == 1 or pow(z, e, rr) == 1)
 return t
def ec_attack(cfg):
 bits, dd, K = cfg["bits"], 1 << cfg["logd"], float(cfg["K"])
 trials, seed = cfg["trials"], cfg["seed"]
 KEYB, NODEB, HEADB = cfg["keybits"], cfg["nodebits"], cfg["headbits"]
 STOREB = KEYB - NODEB - HEADB
 rng = random.Random(seed)
 E, rr, P = find_curve(bits, dd, rng)
 p, mm = E.p, (rr - 1) // dd
 beta = next(z for z in (pow(g, (p - 1) // 3, p) for g in range(2, 100)) if z != 1)
 lamg = next(L for L in (pow(g, (rr - 1) // 3, rr) for g in range(2, 100))
    if L != 1 and E.mul(L, P) == (beta * P[0] % p, P[1]))
 U = [1, rr - 1, lamg, rr - lamg, lamg * lamg % rr, rr - lamg * lamg % rr]
 assert all(pow(u, mm, rr) == 1 for u in U)
 chi, gen = character(rr, dd)
 size = math.isqrt(int(K * mm / 6)) + 1
 Qr = Fraction(cfg.get('q', '1')); size_b = size_g = size
 if Qr != 1:
  size_g = 1
  while 6 * (-(-(Qr.numerator * size_g) // Qr.denominator)) * size_g < K * mm: size_g += 1
  size_b = -(-(Qr.numerator * size_g) // Qr.denominator); size = max(size_b, size_g)
 X2, GMAX = 1 << cfg["logX2"], 1 << cfg["logG"]
 classes = [[] for _ in range(dd)]
 for b in range(1, X2 + 1):
  classes[chi(b)].append(b)
 plans = [plan(L, GMAX) if L else (None, None) for L in classes]
 pair_frac = (sum(2 * st["pairs"] for _, st in plans if st) /
         sum(st["n"] for _, st in plans if st))
 X1 = X2 + int(X2 ** 0.921)
 while True:
  base_primes = primes_in(X2 + 1, X1)
  bases, cnt, baby_n, giant_n = [], 0, None, None
  for a in base_primes:
   bases.append(a)
   cnt += len(classes[(-chi(a)) % dd])
   if baby_n is None:
    if cnt >= size_b:
     baby_n, cnt = len(bases), 0
   elif cnt >= size_g:
    giant_n = len(bases)
    break
  if baby_n and giant_n:
   break
  X1 *= 2
 assert baby_n and giant_n
 def expected(a, limit):
  vals, c = [], (-chi(a)) % dd
  L, (ops, _) = classes[c], plans[c]
  for op in ops:
   positions = (op[1],) if op[0] in (0, 1) else (op[3], op[2])
   for position in positions:
    vals.append(a * L[position])
    if len(vals) == limit:
     return vals
  return vals
 def walk(base_point, a, emit, limit, checks):
  produced, c = 0, (-chi(a)) % dd
  L, (ops, _) = classes[c], plans[c]
  Qa = E.mul(a, base_point)
  tab = {}
  def T(g):
   if g not in tab:
    tab[g] = E.mul(g, Qa)
   return tab[g]
  pts, root = {}, None
  for op in ops:
   if op[0] == 0:
    root = E.mul(L[op[1]], Qa)
    pts[op[1]] = root
    emit(a * L[op[1]], root)
    produced += 1
   elif op[0] == 1:
    root = E.single(root, T(op[2]))
    pts[op[1]] = root
    emit(a * L[op[1]], root)
    produced += 1
   else:
    _, x, i, j, g = op
    Pp, Pm = E.pair(pts[x], T(g))
    pts[j], pts[i] = Pp, Pm
    emit(a * L[j], Pp)
    produced += 1
    if produced < limit:
     emit(a * L[i], Pm)
     produced += 1
   if checks[0] > 0 and produced % 41 == 0:
    v = a * L[op[1] if op[0] in (0, 1) else op[3]]
    assert E.mul(v, base_point) == pts[op[1] if op[0] in (0, 1) else op[3]]
    checks[0] -= 1
   if produced >= limit:
    return produced
  return produced
 NNODES, NHEADS = 1 << NODEB, 1 << HEADB
 HEAD_CAP, KEPT_CAP, CAND_CAP = cfg["head_cap"], 4, cfg["cand_cap"]
 mask = (1 << KEYB) - 1
 keyf = lambda Pt: (Pt[0] >> 1) & mask
 keys3 = lambda Pt: [(Pt[0] >> 1) & mask, (beta * Pt[0] % p >> 1) & mask, (beta * beta * Pt[0] % p >> 1) & mask]
 split = lambda k: (k >> (KEYB - NODEB), (k >> STOREB) & (NHEADS - 1),
             k & ((1 << STOREB) - 1))
 succ = inv_t = mul_t = cand_total = cmp_total = rw_mul_t = 0
 rewalk_fails = cand_over = 0
 occ_chi2_sum = 0.0
 spot = [cfg["spotchecks"]]
 for _ in range(trials):
  s = rng.randrange(1, rr)
  t_secret = pow(s, dd, rr)
  A1, Ad = E.mul(s, P), E.mul(t_secret, P)
  E.inv = E.muls = 0
  heads = {}
  candidates = []
  emitted_b, emitted_g = [], []
  def baby_emit(idv):
   def emit(v, Pt):
    k = keyf(Pt)
    emitted_b.append(v)
    node, head, st = split(k)
    h = heads.setdefault((node, head), [])
    h.append((st, idv))
    assert len(h) <= HEAD_CAP
   return emit
  giant_state = {}
  def giant_emit(v, Pt):
   emitted_g.append(v)
   for k in keys3(Pt):
    node, head, st = split(k)
    chain = heads.get((node, head), ())
    for (st_b, idv) in chain:
     if st_b == st:
      candidates.append((k, idv, giant_state["id"], v))
    giant_state["cmp"] += len(chain)
  remaining = size_b
  for bi in range(baby_n):
   if remaining <= 0:
    break
   remaining -= walk(P, bases[bi], baby_emit(bi), remaining, spot)
  assert len(emitted_b) == size_b
  lam_occ = size_b / (NNODES * NHEADS)
  cells = NNODES * NHEADS
  obs = [len(h) for h in heads.values()]
  chi2 = sum((o - lam_occ) ** 2 for o in obs)
  occ_chi2_sum += chi2 / (size_b * (1 - 1 / cells))
  giant_state["id"], giant_state["cmp"] = baby_n, 0
  remaining = size_g
  for gi in range(baby_n, giant_n):
   if remaining <= 0:
    break
   remaining -= walk(Ad, bases[gi], giant_emit, remaining, spot)
   giant_state["id"] = gi + 1
  assert len(emitted_g) == size_g
  walk_inv, walk_mul = E.inv, E.muls
  exp_b = []
  for bi in range(baby_n):
   exp_b += expected(bases[bi], size_b)
   if len(exp_b) >= size_b:
    break
  exp_b = exp_b[:size_b]
  exp_g = []
  for gi in range(baby_n, giant_n):
   exp_g += expected(bases[gi], size_g)
   if len(exp_g) >= size_g:
    break
  exp_g = exp_g[:size_g]
  assert exp_b == emitted_b and exp_g == emitted_g
  inverses = [pow(g, -1, rr) for g in exp_g]
  exact_ratios = {u * b * inv % rr for u in U for b in exp_b for inv in inverses}
  hit = t_secret in exact_ratios
  cand_total += len(candidates)
  if len(candidates) > CAND_CAP:
   cand_over += 1
  solved = False
  for (k, idb, idg, gv) in candidates[:CAND_CAP]:
   def collect(base_point, a, key, kf, lim):
    kept = []
    def emit(v, Pt):
     if key in kf(Pt):
      kept.append(v)
    walk(base_point, a, emit, lim, [0])
    return kept
   keptb = collect(P, bases[idb], k, lambda Pt: [keyf(Pt)], size_b)
   keptg = collect(Ad, bases[idg], k, keys3, size_g)
   if len(keptb) > KEPT_CAP or len(keptg) > KEPT_CAP:
    rewalk_fails += 1
    continue
   done = False
   for bv in keptb:
    for gvv in keptg:
     for u in U:
      z = u * bv * pow(gvv, -1, rr) % rr
      if E.mul(z, P) == Ad:
       s0 = dth_root(z, dd, rr)
       w = pow(gen, (rr - 1) // dd, rr)
       sol = next(s0 * pow(w, i, rr) % rr for i in range(dd)
                     if E.mul(s0 * pow(w, i, rr) % rr, P) == A1)
       assert sol == s, "wrong scalar"
       done = True
     if done:
      break
    if done:
     break
   if done:
    solved = True
    break
  assert solved == hit, f"resolution mismatch solved={solved} exact={hit}"
  cmp_total += giant_state["cmp"]
  inv_t += walk_inv
  mul_t += walk_mul
  rw_mul_t += E.muls - walk_mul
  if solved:
   succ += 1
 Ke = 6 * size_b * size_g / mm
 return {
  "p_bits": round(math.log2(p), 2), "r": rr, "d": dd, "X2": X2, "GMAX": GMAX,
  "GMAX_over_d2": round(GMAX / dd / dd, 1), "X2_over_GMAX": round(X2 / GMAX, 1),
  "bases_baby_giant": [baby_n, giant_n - baby_n], "sizes_baby_giant": [size_b, size_g], "q": str(Qr), "keys": "x; giants x, beta x, beta^2 x",
  "plan_pair_fraction": round(pair_frac, 4), "K_eff": round(Ke, 4),
  "predicted": round(1 - math.exp(-Ke), 4),
  "keybits": KEYB, "nodes": NNODES, "heads_per_node": NHEADS,
  "trials": trials, "successes": succ, "success_rate": round(succ / trials, 4),
  "coverage_events_match_exact_each_trial": True,
  "candidates_per_trial": round(cand_total / trials, 3),
  "false_match_model": round(3 * size_b * size_g / (1 << KEYB), 3),
  "comparisons_per_trial": round(cmp_total / trials, 1),
  "comparison_model": round(3 * size_b * size_g / (NNODES * NHEADS), 1),
  "occupancy_chi2_scaled_mean": round(occ_chi2_sum / trials, 4),
  "rewalk_cap_fails": rewalk_fails, "candidate_cap_over": cand_over,
  "all_recovered_scalars_correct": True,
  "mults_per_element": round(mul_t / trials / (size_b + size_g), 3),
  "rewalk_mults_per_element_at_cand_rate": round(rw_mul_t / trials / (size_b + size_g), 3),
  "reference_mults_per_element": 6.0,
 }
SMOOTH_PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
         67, 71, 73, 79, 83, 89, 97)
def mark(marks, sbd, u, M):
 for t in sbd: marks[(t - u) % M] = 1
def census(cfg):
 m_bits, k = cfg["m_bits"], cfg["k"]
 K, seed = float(cfg["K"]), cfg["seed"]
 m_target = 1 << m_bits
 c_lo = m_target // 6
 smooths = []
 def gen(idx, val):
  if idx == len(SMOOTH_PRIMES):
   smooths.append(val)
   return
  v = val
  while v <= 4 * c_lo:
   gen(idx + 1, v)
   v *= SMOOTH_PRIMES[idx]
 gen(0, 1)
 smooths = sorted({s for s in smooths if s >= c_lo})
 count = 0
 r = None
 for c in smooths:
  rr_try = 6 * c * (1 << k) + 1
  if is_prime(rr_try) and ((rr_try - 1) // (1 << k)) % 6 == 0 \
    and pow(2, (rr_try - 1) // (1 << k), rr_try) != 1:
   if count == seed:
    r = rr_try
    break
   count += 1
 assert r is not None, "no toy prime for this seed"
 d = 1 << k
 m = (r - 1) // d
 n_all = r - 1
 fac = {}
 for q in (2, 3):
  while n_all % q == 0:
   fac[q] = fac.get(q, 0) + 1
   n_all //= q
 rest, q = n_all, 5
 while q * q <= rest:
  while rest % q == 0:
   fac[q] = fac.get(q, 0) + 1
   rest //= q
  q += 2
 if rest > 1:
  fac[rest] = fac.get(rest, 0) + 1
 gen = None
 for cand in range(2, 100000):
  if all(pow(cand, (r - 1) // q, r) != 1 for q in fac) and pow(cand, r - 1, r) == 1:
   gen = cand
   break
 assert gen is not None
 parts = []
 for q, e in fac.items():
  pe = q ** e
  beta = pow(gen, (r - 1) // pe, r)
  gam = pow(beta, q ** (e - 1), r)
  b = int(math.isqrt(q)) + 1
  baby = {}
  cur = 1
  for j in range(b):
   baby[cur] = j
   cur = cur * gam % r
  gamminv = pow(pow(gam, b, r), -1, r)
  parts.append((q, e, pe, beta, gam, b, baby, gamminv))
 def dlog_q(t, b, baby, gamminv):
  cur = t
  for i in range(b + 1):
   j = baby.get(cur)
   if j is not None:
    return i * b + j
   cur = cur * gamminv % r
  raise ValueError("dlog_q")
 def dlog(x):
  crt, mods = [], []
  for q, e, pe, beta, gam, b, baby, gamminv in parts:
   xq = pow(x, (r - 1) // pe, r)
   digit = 0
   for i in range(e):
    t = pow(xq, q ** (e - 1 - i), r)
    di = dlog_q(t, b, baby, gamminv)
    digit += di * q ** i
    xq = xq * pow(beta, (pe - (di * q ** i) % pe) % pe, r) % r
   crt.append(digit)
   mods.append(pe)
  X, Nm = 0, 1
  for a, mo in zip(crt, mods):
   t = (a - X) % mo
   X = X + Nm * (t * pow(Nm, -1, mo) % mo)
   Nm *= mo
  return X % (r - 1)
 assert dlog(2) % d != 0, "2 in H_t"
 Qr = Fraction(cfg.get('q', '1')); NB = NG = math.ceil(math.sqrt(K * m / 6))
 if Qr != 1:
  NG = math.isqrt(int(K * m / 6 / Qr)) // 2
  while 6 * (-(-(Qr.numerator * NG) // Qr.denominator)) * NG < K * m: NG += 1
  NB = -(-(Qr.numerator * NG) // Qr.denominator)
 N = max(NB, NG)
 nbar = cfg.get('nbar')
 e_t = max(2, round(NG / (nbar * d))) if nbar else max(2, round(N ** float(cfg.get('e_exp', '0.4811'))))
 X2 = d * e_t
 if cfg.get("window") == "top":
  I0 = r // X2
  assert I0 > X2, "toy too small for top-window control"
 else:
  I0 = X2
 Ilen = int(I0 ** 0.921)
 n_side = math.ceil(N / e_t)
 primes_I = primes_in(I0 + 1, I0 + Ilen)
 need = (NB + NG) * 13 // (10 * e_t) + 50 if nbar else 2 * (n_side + 8) + 10
 while len(primes_I) < need:
  Ilen *= 2
  primes_I = primes_in(I0 + 1, I0 + Ilen)
 spf = list(range(X2 + 1))
 for q in range(2, int(math.isqrt(X2)) + 1):
  if spf[q] == q:
   for mm2 in range(q * q, X2 + 1, q):
    if spf[mm2] == mm2:
     spf[mm2] = q
 Ep = {}
 for q in range(2, X2 + 1):
  if spf[q] == q:
   Ep[q] = dlog(q)
 modn = r - 1
 E = [0] * (X2 + 1)
 for b in range(2, X2 + 1):
  q = spf[b]
  E[b] = (E[b // q] + Ep[q]) % modn if q != b else Ep[q] % modn
 classes = {c2: [] for c2 in range(d)}
 for b in range(1, X2 + 1):
  classes[E[b] % d].append(b)
 modS = (r - 1) // 6
 M = m // 6
 def emitted(bs, N):
  vals, cnt, used = [], 0, 0
  for a in bs:
   if cnt >= N:
    break
   used += 1
   Ea = dlog(a)
   L = classes.get((-Ea) % d) or []
   take = min(len(L), N - cnt)
   vals.extend((Ea + E[b]) % modS for b in L[:take])
   cnt += take
  return vals, cnt, used
 sb, cb, nb = emitted(primes_I, NB)
 if nbar:
  gpc, gb, tot, i = cfg.get("gpc", 1), [], 0, nb
  while tot < NG:
   Ea = dlog(primes_I[i]); i += 1; L = classes[(-Ea) % d]
   gb.append((Ea % d, [(Ea + E[b]) % modS // d for b in L])); tot += len(L)
  byc = {}
  for j, g in enumerate(gb): byc.setdefault(g[0], []).append(j)
  R = -(-max(map(len, byc.values())) // gpc)
  rr = [[byc[c][x] for c in sorted(byc) for x in range(j * gpc, min(len(byc[c]), (j + 1) * gpc))] for j in range(R)]
  inc, pos = [], 0
  for o in rr: inc.append(list(range(pos, pos + len(o)))); pos += len(o)
  sbd = [t // d for t in sb]
  def sweep(order):
   marks, out, cum = bytearray(M), [], 0
   for grp in order:
    for j in grp:
     for u in gb[j][1]: mark(marks, sbd, u, M)
     cum += len(gb[j][1])
    Kj = 6 * NB * cum / m
    out.append([round(Kj, 4), round(marks.count(1) / M - 1 + math.exp(-Kj), 6)])
   return out
  rng, tp, dp = random.Random(1), 0, 0
  bch = [(-dlog(a)) % d for a in primes_I[:nb]]
  for _ in range(cfg.get("lit", 500)):
   A, B = classes[rng.choice(bch)], classes[(-rng.choice(gb)[0]) % d]
   tp += len(A) * len(B); dp += len({(x // math.gcd(x, y), y // math.gcd(x, y)) for x in A for y in B})
  eps, a, b = 1 - dp / tp, sweep(rr), sweep(inc)
  return {"r": r, "d": d, "m": m, "X2": X2, "N": [NB, tot], "rounds": R, "eps_lit": round(eps, 6), "production_order": a, "increasing_order": b,
    "lowest": [min(x[1] for x in a), min(x[1] for x in b)], "lowest_minus_literal": [round(min(x[1] + eps * x[0] * math.exp(-x[0]) for x in o), 6) for o in (a, b)]}
 sg, cg, ng = emitted(primes_I[nb:], NG)
 assert (cb, cg) == (NB, NG)
 if cfg.get("multiplicity"):
  nu = [0] * M
  for ti in sb:
   for uj in sg:
    nu[((ti - uj) % modS) // d] += 1
  marks = bytes(1 if v else 0 for v in nu)
  e1 = sum(nu)
  e2 = sum(v * v for v in nu)
  energy = e2 / (e1 + e1 * e1 / M)
  mult = {"energy_ratio": round(energy, 4),
    "mean_nu_marked": round(e1 / max(1, sum(marks)), 4),
    "max_nu": max(nu)}
 else:
  marks = bytearray(M)
  for ti in sb:
   for uj in sg:
    marks[((ti - uj) % modS) // d] = 1
  mult = None
 cov = sum(marks) / M
 Kr = 6 * NB * NG / m
 pred = 1 - math.exp(-Kr)
 repfrac = 1 - sum(marks) / (NB * NG)
 return {
  "r": r, "d": d, "m": m, "window": cfg.get("window", "low"), "K": round(Kr, 6),
  "N": [NB, NG], "q": str(Qr), "n_side": nb,
  "e_t": e_t, "X2": X2, "coverage": round(cov, 12), "predicted": round(pred, 12),
  "diff": round(cov - pred, 12), "repeat_fraction": round(repfrac, 6),
  "uniform_repeat_fraction": round(1 - (1 - math.exp(-Kr)) / Kr, 5),
  **({"multiplicity": mult} if mult else {}),
 }
def planner_mode(cfg):
 out = []
 for spec in cfg["specs"]:
  dd, logX2, seed = spec["d"], spec["logX2"], spec["seed"]
  rng = random.Random(seed)
  rr = toy_prime(28 + int(math.log2(dd)) + 4, dd, rng)
  X2 = 1 << logX2
  cl = class_lists(rr, dd, X2)
  GMAX = 16 * spec.get("kappa", 1) * dd * dd
  rows = []
  for c in range(min(spec.get("classes", 8), dd)):
   L = cl[c]
   if len(L) < 50:
    continue
   _, st = plan(L, GMAX, max_tree=spec.get("max_tree", 1 << 20))
   rows.append(st)
  if rows:
   out.append({
    "r": rr, "d": dd, "X2": X2, "GMAX": GMAX, "GMAX_over_d2": GMAX // (dd * dd),
    "span_over_2GMAX": round(X2 / (2 * GMAX), 1),
    "classes": len(rows),
    "list_lengths": [min(st["n"] for st in rows), max(st["n"] for st in rows)],
    "min_pair_fraction": round(min(st["pair_fraction"] for st in rows), 4),
    "mean_pair_fraction": round(
     sum(st["pair_fraction"] for st in rows) / len(rows), 4),
    "max_cycle_rejections": max(st["cycle_rejections"] for st in rows),
    "max_depth_seen": max(st["max_depth"] for st in rows),
    "depth_cap": 1 << 14, "max_tree_size_seen": max(st["max_tree_size"] for st in rows),
    "size_rejections": sum(st["size_rejections"] for st in rows),
    "tree_cap": spec.get("max_tree", 1 << 20),
    "roots_per_class": [st["roots"] for st in rows],
    "min_leaf_fraction": round(min(st["leaf_fraction"] for st in rows), 5),
    "mean_leaf_fraction": round(sum(st["leaf_fraction"] for st in rows) / len(rows), 5),
    "centres_hosting": [st["centres_hosting"] for st in rows], "fallback_new_centres": [st["fallback_new_centres"] for st in rows],
    "max_ops_per_element": max(st["max_ops_per_element"] for st in rows), "ops_bound": OPS_BOUND,
   })
 return out
def barrier_tree(nodes, fanin):
 """Fan-in tree: (internal nodes, sum children^2)."""
 internal = sq = 0
 level = nodes
 while level > 1:
  parents = -(-level // fanin)
  full, rest = divmod(level, fanin)
  sq += full * fanin * fanin + (rest * rest if rest else 0)
  internal += parents
  level = parents
 return internal, sq
WAKE = 256
def barrier_instr(nodes, fanin, count):
 """Barriers: arrive/wake/sweep/release (worst)."""
 internal, sq = barrier_tree(nodes, fanin)
 arrivals = nodes + internal - 1
 return count * (arrivals * 4 + arrivals * WAKE + sq * 2 + arrivals * (WAKE + 2)), internal
LINK_WAVE = 16
UNT, UNK = 64, 6
REC_GEN = 12; REC_RELAY, REC_INS, REC_LOOK = (Fraction(UNT * b + 2, UNT) for b in (12, 10, 4))
CMP_EACH, CAP_EACH, EXT_EACH = 7, 2, 3
GEN12 = ["ubfx","add","ldp","ldp","cmp","b.hs","stp","stp","str","str","ldp","ldr"]
RELAY_BODY = ["ldp","ldp","ldr","ubfx","add","ldp","cmp","b.hs","stp","stp","str","str"]
RELAY_U = RELAY_BODY * UNT + ["subs","b.ne"]
UNROLL_SETUP = ["and","lsr","cbz"]
UNROLL_FIX = len(UNROLL_SETUP) + UNK + 1
def unroll(name, body, U=UNT):
 """s29: U-fold loop (one subs + b.ne per U records) + binary tail ladder; label SKIP renamed per copy."""
 k = U.bit_length() - 1; assert U == 1 << k >= 2; c = iter(range(4 * U)); b = lambda: body.replace("SKIP", "S%d" % next(c))
 L = [name + ":", " and x5, x23, #%d" % (U - 1), " lsr x23, x23, #%d" % k, " cbz x23, tail", "pair:"] + [b() for _ in range(U)] + [" subs x23, x23, #1", " b.ne pair", "tail:"]
 for j in range(k - 1, -1, -1): L += [" tbz x5, #%d, T%d" % (j, j)] + [b() for _ in range(1 << j)] + ["T%d:" % j]
 return "\n".join(L + [" ret"])
INS_BODY = """ ldp x0, x1, [x20], #16
 ldp x2, x3, [x20], #16
 ldr x4, [x20], #8
 lsr x9, x1, #5
 ldr x11, [x21, x9, lsl #3]
 extr x10, x1, x0, #60
 bfi x4, x10, #55, #9
 stp x3, x4, [x22, #16]
 str x22, [x21, x9, lsl #3]
 stp x2, x11, [x22], #32"""
INS_U = unroll("ins2", INS_BODY)
def audit_ins12():
 """s24/s29: INS_U loop, executed, oracle memory."""
 import a64
 from collections import defaultdict
 def run(n, seed):
  r = random.Random(seed); P = (1 << 62) + 32; mem = defaultdict(int); want = {}
  for i in range(n):
   w = [r.getrandbits(64), r.randrange(40) << 6 | r.getrandbits(6), r.getrandbits(64), r.getrandbits(64), r.getrandbits(55)]
   for j in range(5): mem[2**40 + 40 * i + 8 * j] = w[j]
   a = P + 32 * i; h = 8 * (w[1] >> 5)
   want.update({a: w[2], a + 8: want.get(h, 0), a + 16: w[3], a + 24: w[4] | ((w[1] << 64 | w[0]) >> 60 & 511) << 55, h: a})
  want.update(mem); o = {}
  st = a64.a64_run(*a64.a64_load(INS_U), "ins2", {"x20": 2**40, "x21": 0, "x22": P, "x23": n}, mem, out_regs=o)
  assert st == 10 * n + 2 * (n // UNT) + UNROLL_FIX and o["x22"] == P + 32 * n and mem == want
 for n, sd in ((600, 22), (601, 23), (1, 24), (0, 25), (639, 26), (38, 27)): run(n, sd)
LOOK_BODY = """ ldp x0, x1, [x20], #40
 lsr x9, x1, #5
 ldr x11, [x21, x9, lsl #3]
 cbz x11, SKIP
 subs x24, x24, #1
 extr x10, x1, x0, #60
 ldur x2, [x20, #-24]
SKIP:"""
LOOK_U = unroll("look2", LOOK_BODY)
EXT3 = ["subs","extr","ldur"]
def audit_look6():
 """s24/s29: LOOK_U loop, executed."""
 import a64
 from collections import defaultdict
 C2, C10 = 0xC0FFEE, 0x5EED
 def chk(txt, N, seed):
  r = random.Random(seed); B = 2**40; H = 40
  rec = [[r.getrandbits(64), r.randrange(H) << 6 | r.getrandbits(6), r.getrandbits(64)] for _ in range(N)]
  hid = [w[1] >> 5 for w in rec]
  head = [0 if r.random() < 0.7 else (1 << 62) + 32 + 32 * r.randrange(1000) for _ in range((max(hid) + 1) if hid else 1)]
  ne = sum(1 for h in hid if head[h])
  pg, lb = a64.a64_load(txt)
  m = defaultdict(int)
  for i, w in enumerate(rec):
   for j in range(3): m[B + 40*i + 8*j] = w[j]
  for h in range(len(head)): m[8*h] = head[h]
  sn = dict(m); o = {}
  assert a64.a64_run(pg, lb, "look2", {"x20": B, "x21": 0, "x23": N, "x24": N}, m, out_regs=o) == 4*N + 2*(N // UNT) + 3*ne + UNROLL_FIX
  assert o["x20"] == B + 40*N and o["x24"] == N - ne and dict(m) == sn
  nz = next((i for i in range(N) if head[hid[i]]), None)
  if nz is not None:
   wz = rec[nz]; o3 = {}
   a64.a64_run(pg, lb, "look2", {"x20": B + 40*nz, "x21": 0, "x23": 1, "x24": 5, "x2": C2, "x10": C10}, m, out_regs=o3)
   assert o3["x2"] == wz[2]
  return rec, hid, head, pg, lb, B
 rec, hid, head, pg, lb, B = chk(LOOK_U, 600, 23)
 for n, sd in ((601, 24), (1, 25), (0, 26), (639, 27), (38, 28)): chk(LOOK_U, n, sd)
 N = 600
 i1 = next(i for i in range(N) if head[hid[i]]); i0 = next(i for i in range(N) if not head[hid[i]])
 m = defaultdict(int)
 for i, w in enumerate(rec):
  for j in range(3): m[B + 40*i + 8*j] = w[j]
 for h in range(len(head)): m[8*h] = head[h]
 for idx, want in ((i1, True), (i0, False)):
  w = rec[idx]; o1 = {}
  st = a64.a64_run(pg, lb, "look2", {"x20": B + 40*idx, "x21": 0, "x23": 1, "x24": 5, "x2": C2, "x10": C10}, m, out_regs=o1)
  h = w[1] >> 5
  assert st == UNROLL_FIX + (7 if want else 4) and o1["x20"] == B + 40*idx + 40 and o1["x9"] == h and o1["x11"] == head[h]
  if want:
   assert o1["x2"] == w[2] and o1["x10"] & 0x1FF == ((w[1] << 64 | w[0]) >> 60) & 0x1FF and o1["x24"] == 4
  else:
   assert o1["x2"] == C2 and o1["x10"] == C10 and o1["x24"] == 5
 for nm, v in (("adv", LOOK_U.replace("#40", "#16")), ("off", LOOK_U.replace("#-24", "#-16")), ("lsr", LOOK_U.replace("x1, #5", "x1, #4")), ("drop", LOOK_U.replace(" ldur x2, [x20, #-24]\n", ""))):
  try: chk(v, 600, 23)
  except AssertionError: continue
  raise AssertionError("negative passed: " + nm)
K10_LO, HEAD_LSR, INS_BFI, CMP_AND_POS = 60, 5, (55, 9), (55, 9)
CMP_AND_MASK = 0xFF80000000000000
M64 = (1 << 64) - 1
def cmp_field(c4, x10):
 """The displayed field compare: and b,c4,#CMP_AND_MASK ; cmp b,x10,lsl #55 (x10 = extr of the giant key, unmasked)."""
 return (c4 & CMP_AND_MASK) == (x10 << CMP_AND_POS[0]) & M64
def cmp_copy(base, link, lab):
 """One chain comparison (s25): cap pair 2 + CMP 7; the str opens candidate formation (own row)."""
 return """SKIP%s:
 subs x25, x25, #1
 b.lo fail
 ldp x12, %s, [%s]
 ldr x14, [%s, #24]
 and x14, x14, #0xFF80000000000000
 cmp x14, x10, lsl #55
 ccmp x12, x2, #0, eq
 b.ne SKIP%sN
 str %s, [x26], #8
SKIP%sN:""" % (lab, link, base, base, lab, base, lab)
CMP_A = cmp_copy("x11", "x13", "A") + "\n cbz x13, SKIP"
CMP_B = cmp_copy("x13", "x11", "B") + "\n cbnz x11, SKIPA"
CMP7 = [l.split()[0] for l in CMP_A.split("\n") if l[0] == " "][2:8] + ["cbz"]
LOOKC_BODY = LOOK_BODY.replace("SKIP:", CMP_A + "\n" + CMP_B + "\nSKIP:")
LOOKC_U = unroll("look2", LOOKC_BODY) + "\nfail:\n mov x27, #1\n ret"
def audit_cmp(ins=None, look=None):
 """s25: executed insert and lookup loops with the inline comparator against the full-key oracle."""
 import a64
 from collections import defaultdict
 ins, look = ins or INS_U, look or LOOKC_U
 r = random.Random(20261109); P, T, G, CL = (1 << 62) + 32, 1 << 40, 1 << 41, 1 << 42
 def key(f9, head, k64, own=0): return own | f9 << 60 | head << 69 | k64 << 128
 baby = [key(0, 1, 0), key(0, 1, 7), key(511, 1, 7), key(256, 2, 5), key(1, 3, 1 << 63), key(3, 4, 9, 5), key(3, 4, 9, 6)]
 baby += [key(r.getrandbits(9), r.randrange(8), r.getrandbits(64), r.getrandbits(60)) for _ in range(120)]
 giant = [1 << 69, key(0, 0, 0), key(0, 1, 7), key(1, 1, 7), key(511, 1, 7), key(256, 2, 5), key(0, 2, 5), key(1, 3, 1 << 63), key(1, 3, 0), key(3, 4, 9, 77)]
 for k in baby[7:60]:
  giant += [k, k ^ 1 << r.randrange(60, 69), k ^ 1 << r.randrange(128, 192), k ^ 1 << 69, k ^ r.getrandbits(60)]
 giant += [key(r.getrandbits(9), r.randrange(10), r.getrandbits(64)) for _ in range(40)]
 mem = defaultdict(int); chains = defaultdict(list)
 for i, k in enumerate(baby):
  w = [k & M64, k >> 64 & M64, k >> 128, r.getrandbits(64), r.getrandbits(55)]
  for j in range(5): mem[T + 40 * i + 8 * j] = w[j]
  chains[k >> 69 & (1 << 59) - 1].insert(0, (k, P + 32 * i))
 a64.a64_run(*a64.a64_load(ins), "ins2", {"x20": T, "x21": 0, "x22": P, "x23": len(baby)}, mem)
 for i, k in enumerate(giant):
  for j in range(3): mem[G + 40 * i + 8 * j] = k >> 64 * j & M64
 want, C, ne = [], 0, 0
 for k in giant:
  ch = chains.get(k >> 69 & (1 << 59) - 1, [])
  ne += bool(ch); C += len(ch); want += [a for (bk, a) in ch if bk >> 60 == k >> 60]
 odd = sum(1 for k in giant for (bk, a) in chains.get(k >> 69 & (1 << 59) - 1, []) if bk >> 60 == k >> 60 and k >> 69 & 1)
 assert len(want) >= 60 and (1 << 62) + 32 in want and C > 2 * len(giant) and 10 < odd < len(want) - 10
 snap = dict(mem); o = {}; n = len(giant)
 st = a64.a64_run(*a64.a64_load(look), "look2", {"x20": G, "x21": 0, "x23": n, "x24": n, "x25": C, "x26": CL}, mem, out_regs=o)
 got = [mem[CL + 8 * i] for i in range((o["x26"] - CL) // 8)]
 assert got == want and o["x27"] == 0 and o["x25"] == 0 and o["x24"] == n - ne
 assert st == 4 * n + 2 * (n // UNT) + 3 * ne + (CMP_EACH + CAP_EACH) * C + len(want) + UNROLL_FIX
 assert {a: v for a, v in mem.items() if not CL <= a < CL + 8 * len(want) and (v or a in snap)} == snap
 o2 = {}; a64.a64_run(*a64.a64_load(look), "look2", {"x20": G, "x21": 0, "x23": n, "x24": n, "x25": C - 1, "x26": CL}, mem, out_regs=o2)
 assert o2["x27"] == 1
 return {"baby": len(baby), "giant": n, "comparisons": C, "matches": len(want), "matches_with_key_bit_69_set": odd, "nonempty": ne, "steps": st}
def audit_packed_match():
 assert (K10_LO, HEAD_LSR, INS_BFI, CMP_AND_POS) == (60, 5, (55, 9), (55, 9))
 assert CMP_AND_MASK == (((1 << 9) - 1) << 55) == 0xFF80000000000000 and K10_LO + 9 == 69 == 64 + HEAD_LSR and 128 - 69 == 59
 assert CMP7 == ["ldp","ldr","and","cmp","ccmp","b.ne","cbz"] and CMP_A.count("\n ") == CMP_B.count("\n ") == CMP_EACH + CAP_EACH + 1
 def words(k): return (k & M64, (k >> 64) & M64, (k >> 128) & M64)
 def ex(k): return (k >> K10_LO) & M64
 def f9(k): return (k >> K10_LO) & 0x1FF
 def pk(k, junk=0): return (junk & ((1 << 55) - 1)) | f9(k) << INS_BFI[0]
 def slot(k): return (k >> 69) & ((1 << 59) - 1)
 r = random.Random(20261103)
 assert slot(1 << 69) == 1 and f9(1 << 69) == 0 and ex(1 << 69) == 512 and cmp_field(pk(1 << 69), ex(1 << 69))
 assert not (pk(1 << 69) & 0x7FC0000000000000) == (ex(1 << 69) << 54) & M64
 kk = [0, (1 << 192) - 1, 1 << 60, ((1 << 192) - 1) ^ (1 << 60), 0x3FF << 60, ((0x3FF << 60) | 0xDEAD), 1 << 63, 1 << 64, (1 << 60) - 1, (1 << 70) - 1, 1 << 68, 1 << 69, 3 << 68]
 nb = [x for x in range(192) if x < 60 or x > 68]
 for k in kk + [r.getrandbits(192) for _ in range(3000)]:
  assert (words(k)[1] >> HEAD_LSR) == slot(k) and ex(k) == ((words(k)[1] << 64 | words(k)[0]) >> 60) & M64
  c4 = pk(k, r.getrandbits(64))
  assert cmp_field(c4, ex(k))
  for bit in range(60, 69): assert not cmp_field(c4, ex(k ^ 1 << bit))
  ko = k ^ (1 << r.choice(nb)); assert cmp_field(c4, ex(ko))
  assert slot(k ^ (1 << r.randrange(69, 128))) != slot(k)
def audit_record_path(rows, N, per_point, baby_pt, giant_pt):
 assert (baby_pt, giant_pt) == (BABY_PIN, GIANT_PIN)
 assert len(GEN12) == REC_GEN and len(RELAY_BODY) == 12 and REC_RELAY == Fraction(len(RELAY_U), UNT) == Fraction(770, 64)
 assert INS_BODY.count("\n") + 1 == 10 and INS_U.count(INS_BODY) == 2 * UNT - 1 and REC_INS == Fraction(642, 64)
 assert LOOK_BODY.count("\n") + 1 == 8 and LOOK_U.count("ldp x0, x1, [x20], #40") == 2 * UNT - 1 and REC_LOOK == Fraction(258, 64)
 assert len(EXT3) == EXT_EACH and len(CMP7) == CMP_EACH and len(UNROLL_SETUP) == 3 and UNROLL_FIX == 10
 assert "stp x2, x11, [x22], #32" in INS_BODY
 assert "x5" not in INS_BODY and "x5" not in LOOK_BODY
 assert 60 + 9 + 59 + 64 == 192 and (127 - 69 + 1) == 59 and (68 - 60 + 1) + (191 - 128 + 1) == 9 + 64
 for lsb, w in [(0,30),(30,30),(55,9),(55,9)]:
  assert 0 <= lsb < 64 and 1 <= w <= 64 and lsb + w <= 64
 assert K10_LO < 64 and HEAD_LSR < 64 and CMP_AND_POS[0] < 64
 assert 16 + 16 + 8 == 40 and 16 + 24 == 40
 assert 24 <= 255 and 40 <= 504 and 40 % 8 == 0
 assert 64 + 64 + 9 + 62 + 35 <= 256 and 62 + 2 == 64 and 29 + 35 == 64 and 35 <= 55
 assert (1 << 57) * 32 < 1 << 63
 m = 0xCCCCCCCD
 assert m * 5 == 2 ** 34 + 1 and m < 2 ** 32 and 2 ** 24 * 5 < 2 ** 32
 assert 4 * (2 ** 34 + 1) + 5 * (2 ** 24 - 1) < 5 * 2 ** 34
 assert all((x * 2341) >> 14 == x // 7 for x in range(1792))
 assert 2341 < 2 ** 16 and 1791 * 2341 < 2 ** 32
 assert 8 + 3 + 5 == LINK_WAVE and 4 <= 8 and 8 <= 9 and 11 <= 16
 assert "owner pool-cap compares" not in rows and "head-index word combine" not in rows
 assert REC_GEN + REC_RELAY + REC_INS == Fraction(2180, 64)
 assert EXT3 == ["subs","extr","ldur"] and K10_LO == 60
 audit_packed_match()
 audit_ins12()
 audit_look6()
 cmp_run = audit_cmp()
 for bad in (INS_U.replace("#55, #9", "#54, #9"), None):
  try: audit_cmp(bad, None if bad else LOOKC_U.replace("0xFF80000000000000", "0x7FC0000000000000").replace("lsl #55", "lsl #54"))
  except AssertionError: continue
  raise AssertionError("comparator negative control passed")
 audit_smx()
 b = batch_audit.run({}); assert (b['batch_per_lane'], b['round_fixed']) == (str(BATCH), ROUND_FIXED)
 b['comparator'] = cmp_run
 return b
def proto(cfg):
 from collections import Counter, defaultdict
 rng = random.Random(cfg.get("seed", 2026))
 GB, SB, TB, UB = cfg.get("gen_bits", 3), cfg.get("s_bits", 1), cfg.get("t_bits", 2), cfg.get("u_bits", 2)
 LANES, WAVES, STEPS = cfg.get("lanes", 4), cfg.get("waves", 3), cfg.get("steps", 6)
 CAP1, CAP2, CAPC, FANIN = cfg.get("cap1", 96), cfg.get("cap2", 96), cfg.get("capc", 8), cfg.get("fanin", 4)
 pattern = cfg.get("pattern", "uniform")
 assert GB >= SB
 G, R, O = 1 << GB, 1 << (TB + SB), 1 << (TB + UB)
 NODES = G + R + O
 cost = Counter()
 buf = {}
 tag = {}
 cnt = {}
 ptr = {}
 single = cfg.get("single_count", False)
 emitted, delivered = [], defaultdict(list)
 heads, candidates = defaultdict(list), []
 overflow = []
 def open_buf(sender, receiver, parity, wave):
  key = (sender, receiver, parity)
  lst = buf.setdefault(key, [])
  w, rd = tag.get(key, (None, None))
  if w != wave:
   assert w is None or rd == w + 1, "buffer reused before read"
   lst.clear(); tag[key] = (wave, None)
  return lst
 def put(sender, receiver, parity, rec, cap, wave):
  lst = open_buf(sender, receiver, parity, wave)
  pw, pv = ptr.get((sender, receiver), (wave, 0))
  if pw != wave: assert pv == 0; pw, pv = wave, 0
  if pv >= cap * 40:
   overflow.append((sender, receiver, wave)); return False
  lst.append(rec); ptr[(sender, receiver)] = (wave, pv + 40); return True
 def publish(sender, receiver, parity, wave):
  lst = open_buf(sender, receiver, parity, wave)
  pw, pv = ptr.get((sender, receiver), (wave, 0))
  val = pv if pw == wave else 0
  if pw != wave: assert pv == 0
  assert val % 40 == 0 and val // 40 == len(lst)
  ptr[(sender, receiver)] = (wave, 0)
  cnt[(sender, receiver) if single else (sender, receiver, parity)] = (len(lst), wave)
 def take(sender, receiver, parity, wave):
  key = (sender, receiver, parity)
  w, rd = tag.get(key, (None, None))
  if w is None: return []
  assert w == wave - 1 and rd is None, "read in the wrong wave"
  tag[key] = (w, wave)
  n, cw = cnt[(sender, receiver) if single else key]
  if (n, cw) != (len(buf[key]), w):
   raise RuntimeError("count word mismatch")
  return list(buf[key][:n])
 def relay_of(key, g): return ((key & (R // (1 << SB) - 1)) << SB) | (g & ((1 << SB) - 1))
 def owner_of(key): return key & (O - 1)
 def t_of(x): return x & ((1 << TB) - 1)
 def pack_of(key): return ((key >> K10_LO) & 0x1FF) << INS_BFI[0]
 bk, planted = [], [0]
 def barrier():
  instr, _ = barrier_instr(NODES, FANIN, 1)
  cost["barrier"] += instr
 wave_count = Counter()
 seq = 0
 for phase in (0, 1):
  for wl in range(WAVES + 2):
   w = phase * (WAVES + 2) + wl
   if wl < WAVES:
    for g in range(G):
     for lane in range(LANES):
      for st in range(STEPS):
       n_out = 2 if rng.random() < 0.98 else 1
       for _ in range(n_out):
        if pattern == "one_owner": key = 0 | (rng.getrandbits(180) << 12)
        elif pattern == "two_owners": key = (rng.randrange(2) * (O // 2)) | (rng.getrandbits(180) << 12)
        elif phase and rng.random() < 0.25: key = rng.choice(bk); planted[0] += 1
        else: key = rng.getrandbits(192)
        if not phase: bk.append(key)
        rec = (key, phase, g, lane, seq); seq += 1
        emitted.append(rec)
        cost["rec_gen"] += 1
        if not put(g, G + relay_of(key, g), w & 1, rec, CAP1, w):
         return {"pattern": pattern, "failed": True, "reason": "link cap overflow (run fails)",
           "overflow_at": overflow[0], "records_before_failure": len(emitted)}
     for rl in range(g & ((1 << SB) - 1), R, 1 << SB):
      publish(g, G + rl, w & 1, w)
     cost["link"] += R >> SB
     wave_count["gen"] += 1
   if 1 <= wl <= WAVES:
    for rl in range(R):
     s = rl & ((1 << SB) - 1)
     for g in range(s, G, 1 << SB):
      for rec in take(g, G + rl, (w - 1) & 1, w):
       cost["rec_relay"] += 1
       assert relay_of(rec[0], g) == rl
       if not put(G + rl, G + R + owner_of(rec[0]), (w - 1) & 1, rec, CAP2, w - 1):
        return {"pattern": pattern, "failed": True, "reason": "link cap overflow (run fails)",
          "overflow_at": overflow[0], "records_before_failure": len(emitted)}
     for o in range(O):
      if t_of(o) == (rl >> SB): publish(G + rl, G + R + o, (w - 1) & 1, w - 1)
     cost["link"] += (G >> SB) + O // (1 << TB)
     wave_count["relay"] += 1
   if 2 <= wl <= WAVES + 1:
    for o in range(O):
     t = t_of(o)
     for s in range(1 << SB):
      rl = (t << SB) | s
      for rec in take(G + rl, G + R + o, w & 1, w - 1):
       assert owner_of(rec[0]) == o
       delivered[o].append(rec)
       slot = (rec[0] >> (TB + UB)) & 63
       rest = rec[0] >> (TB + UB + 6)
       if phase == 0:
        heads[(o, slot)].append((rest, pack_of(rec[0]), rec)); cost["insert"] += 1
       else:
        cost["lookup"] += 1
        if heads.get((o, slot)): cost["nonempty"] += 1
        for (rb, chb, bb) in heads[(o, slot)]:
         cost["compare"] += 1
         pkm = cmp_field(chb, (rec[0] >> K10_LO) & M64) and rb == rest
         assert pkm == (rb == rest)
         if pkm: candidates.append((o, rec[0], bb, rec))
     cost["link"] += 1 << SB
     wave_count["owner"] += 1
   barrier(); cost["barriers"] += 1
 got = [rec for o in delivered for rec in delivered[o]]
 assert Counter(got) == Counter(emitted) and len(got) == len(emitted)
 n_b = sum(1 for r in emitted if r[1] == 0); n_g = len(emitted) - n_b
 assert cost["rec_gen"] == cost["rec_relay"] == len(emitted) and cost["insert"] == n_b and cost["lookup"] == n_g
 cand_cost = 0; resolved = 0
 for (o, key, bb, gg) in candidates:
  t = t_of(o); replies = []
  for side in (bb, gg):
   g = side[2]; rl = (t << SB) | (g & ((1 << SB) - 1))
   cand_cost += 2 * REC_RELAY
   assert side in emitted
   rl2 = relay_of(key, g); assert owner_of(key) == o
   cand_cost += 2 * REC_RELAY; replies.append(side[4])
  resolved += 1
 assert resolved == len(candidates) >= planted[0] and (planted[0] > 0 or pattern != "uniform")
 links = G * (R >> SB) + R * (O // (1 << TB))
 expected_link_events = (R >> SB) * wave_count["gen"] + ((G >> SB) + O // (1 << TB)) * wave_count["relay"] + (1 << SB) * wave_count["owner"]
 assert cost["link"] == expected_link_events
 internal, _ = barrier_tree(NODES, FANIN)
 total = (cost["rec_gen"] * REC_GEN + cost["rec_relay"] * REC_RELAY + cost["insert"] * REC_INS + cost["lookup"] * REC_LOOK
       + cost["compare"] * CMP_EACH + cost["link"] * LINK_WAVE + cost["barrier"] + cand_cost)
 formula = (len(emitted) * (REC_GEN + REC_RELAY) + n_b * REC_INS + n_g * REC_LOOK + cost["compare"] * CMP_EACH
         + expected_link_events * LINK_WAVE + barrier_instr(NODES, FANIN, cost["barriers"])[0] + cand_cost)
 assert total == formula
 return {"pattern": pattern, "failed": False, "gens": G, "relays": R, "owners": O, "lanes": LANES, "waves_per_phase": WAVES,
   "records": len(emitted), "baby": n_b, "giant": n_g, "delivered_exactly_once": True,
   "every_buffer_read_in_the_wave_after_it_was_written": True, "barriers": cost["barriers"],
   "barrier_tree_internal_nodes": internal, "candidates": len(candidates), "candidates_resolved": resolved,
   "compares": cost["compare"], "links": links, "link_wave_events": cost["link"],
   "nonempty_lookups": cost["nonempty"], "nonempty_fraction": cost["nonempty"] / max(1, cost["lookup"]),
   "instructions_by_item": {"gen": cost["rec_gen"] * REC_GEN, "relay": str(cost["rec_relay"] * REC_RELAY),
                   "insert": str(cost["insert"] * REC_INS), "lookup": str(cost["lookup"] * REC_LOOK),
                   "compare": cost["compare"] * CMP_EACH, "link_words": cost["link"] * LINK_WAVE,
                   "barriers": cost["barrier"], "candidate_path": str(cand_cost)},
   "instructions_total": str(total), "closed_form_equal": True,
   "scope": "sequential simulation; listing applied to counted events"}
def proto_mode(cfg):
 out = []
 c = dict(cfg); c["pattern"] = "uniform"; out.append(proto(c))
 c = dict(cfg); c.update({"pattern": "two_owners", "cap2": 8 * cfg.get("cap2", 96)}); out.append(proto(c))
 c = dict(cfg); c.update({"pattern": "one_owner"}); out.append(proto(c))
 assert out[2]["failed"] and not out[0]["failed"] and not out[1]["failed"]
 c = dict(cfg); c.update({"pattern": "uniform", "single_count": True})
 try:
  proto(c); clobbered = False
 except RuntimeError:
  clobbered = True
 assert clobbered
 out.append({"pattern": "uniform", "single_count_word_control": "one count word per link is clobbered by the next wave (detected); per-parity words never are"})
 out.append(lane_schedule(cfg.get("schedule", {})))
 return out
def lane_schedule(cfg):
 """Segment-to-lane assignment: padding lemma."""
 import heapq
 rng = random.Random(cfg.get("seed", 7))
 LANES, SEG, TM = cfg.get("lanes", 64), cfg.get("seg", 32), cfg.get("tm", 16)
 res = []
 for kind in ("uniform", "extremes", "long_first", "long_last", "alternating"):
  for S in (LANES, LANES + 1, 3 * LANES + 5, 50 * LANES + 17):
   if kind == "uniform": L = [rng.randint(SEG, SEG + TM) for _ in range(S)]
   elif kind == "extremes": L = [rng.choice((SEG, SEG + TM)) for _ in range(S)]
   elif kind == "long_first": L = sorted(rng.randint(SEG, SEG + TM) for _ in range(S))[::-1]
   elif kind == "long_last": L = sorted(rng.randint(SEG, SEG + TM) for _ in range(S))
   else: L = [SEG + TM if (i // LANES) % 2 == i % 2 else SEG for i in range(S)]
   free = [(0, j) for j in range(LANES)]
   heapq.heapify(free)
   fin = [0] * LANES
   for x in L:
    t, j = heapq.heappop(free)
    fin[j] = t + x
    heapq.heappush(free, (t + x, j))
   span, ops = max(fin), sum(L)
   pad = LANES * span - ops
   assert span - min(fin) <= SEG + TM and pad <= LANES * (SEG + TM) and span <= ops / LANES + SEG + TM
   st = [sum(L[j::LANES]) for j in range(LANES)]
   res.append((kind, S, pad, LANES * max(st) - ops))
 return {"lanes": LANES, "seg": SEG, "treemax": TM, "cases": len(res), "padding_lemma_holds": True,
   "padding_bound": LANES * (SEG + TM),
   "max_greedy_padding": max(r[2] for r in res), "max_stride_padding": max(r[3] for r in res)}
KERNEL = 916
BABY_PIN, GIANT_PIN = Fraction(7230079554123, 6553600000), Fraction(8384332354123, 6553600000)
def exp_neg_bounds(x, terms=40):
 """Rational enclosure of exp(-x), 0<x<1: alternating Taylor, decreasing terms."""
 assert Fraction(0) < x < Fraction(1)
 t, s, lo, up = Fraction(1), Fraction(0), None, None
 for n in range(terms):
  s += t
  if n % 2 == 0: up = s
  else: lo = s
  t = t * (-x) / (n + 1)
 assert lo < up and up - lo < Fraction(1, 10 ** 50)
 return lo, up
def exp4(K):
 """Enclosure of exp(-K), 1 < K < 4, as exp(-K/4)^4."""
 lo, up = exp_neg_bounds(K / 4); assert 0 < lo and 1 < K < 4 and up ** 4 - lo ** 4 < Fraction(1, 10 ** 49)
 return lo ** 4, up ** 4
def success_floor_lower(K):
 """Lower bound of an opened run's success (s2)."""
 lo, up = exp4(K)
 x = 1 - up - DELTA_P - FAILP
 return Fraction(x.numerator * 10 ** 60 // x.denominator, 10 ** 60)
TARGET, GATE_C, KT_RUN, DELTA_P = Fraction(39, 100), 64, Fraction(12, 5), Fraction(1, 2000)
FAILP = Fraction(1, 1 << 20) + Fraction(7, 1 << 40) + Fraction(1, 1 << 110)
def gate(lo):
 Kg = -(-(TARGET * (1 << 64)).numerator * lo.denominator // ((TARGET * (1 << 64)).denominator * lo.numerator))
 assert Fraction(Kg, 1 << 64) * lo >= TARGET > Fraction(Kg - 1, 1 << 64) * lo and 0 < Kg < 1 << 64
 return Kg
def ledger(cfg):
 from decimal import Decimal, localcontext, ROUND_CEILING
 ceil = lambda x: -(-x.numerator // x.denominator) if isinstance(x, Fraction) else x
 r, d, m = R_PRODUCTION, D_PRODUCTION, M_PRODUCTION
 assert (r - 1) % d == 0 and m % 6 == 0
 v2, t_ = 0, r - 1
 while t_ % 2 == 0: t_ //= 2; v2 += 1
 assert v2 == 32
 lam = (-0xd201000000010000) ** 2 - 1
 assert (lam * lam + lam + 1) % r == 0 and pow(lam, m, r) == 1 and pow(2, m, r) != 1
 KT, F = Fraction(cfg.get('KT', KT_RUN)), Fraction(998, 1000)
 ELL, Q = Fraction(193, 200), Fraction(cfg.get('ratio', '23/50'))
 KX, BETA, ROUTE = kernel_check.leaf_worst(), batch_audit.BETA_LEN, 6
 assert (KX, BETA, len(batch_audit.BUILD)) == (590, 128, 24)
 F_PRE = Fraction(39, 40)
 assert F_PRE <= F
 NB_ = lambda g: -(-(Q.numerator * g) // Q.denominator)
 lo_, hi_ = 1, 1 << 120
 while lo_ < hi_:
  mid_ = (lo_ + hi_) // 2
  if 6 * NB_(mid_) * mid_ >= KT * m: hi_ = mid_
  else: lo_ = mid_ + 1
 NG = lo_; NB = NB_(NG); N = NB
 K = Fraction(6 * NB * NG, m)
 assert 6 * NB_(NG - 1) * (NG - 1) < KT * m <= 6 * NB * NG and K >= KT
 X2, GMAX, EMIN, EMAX = 1 << 76, 1 << 47, (1 << 56) - (1 << 52), 1 << 57
 LANES, SEG, TREEMAX, WSTEP, DEPTH = 1 << 19, 1 << 22, 1 << 20, 1 << 33, 1 << 17
 NGT, NG = NG, NG + EMAX
 GENS, SBITS, TBITS, UBITS = 1 << 47, 20, 30, 30
 RELAYS, OWNERS = 1 << (TBITS + SBITS), 1 << (TBITS + UBITS)
 HEADS, POOL, CAP1, CAP2, CAPC, CAND = 1 << 59, 1 << 57, 1 << 24, 1 << 21, 1 << 8, 1 << 44
 RECB, FANIN = 40, 1 << 10
 MARG = Fraction(1, 1)
 SM = 8192000
 SM77 = -(-77 * SM // 255)
 assert SM77 * 255 >= 77 * SM and SM77 < SM and 4 * 77 > 255
 inv_share = Fraction(ROUND_FIXED + 5, LANES)
 assert inv_share < 1 and (ROUND_FIXED, BATCH) == (172585, Fraction(29411, 32))
 tbl_share = ceil(Fraction(232864, LANES))
 PAIR_BK = 8 + 1 + 6 + 1 + 2 + 4 + inv_share
 ROOT_BK = 5 + 6 + 2 + 1 + 2 + inv_share + 5
 TABLE_BK = 8 + 7 + 6 + tbl_share + 5 + 21
 assert TABLE_BK <= 50
 REC_AVG = Fraction(2 * (REC_GEN + REC_RELAY) + REC_INS + REC_LOOK, 2)
 PAIR = 2 * KERNEL + BATCH + MARG * (PAIR_BK + ROUTE)
 ROOT = KERNEL + BATCH + MARG * ROOT_BK
 CAPF = Fraction(7, 32)
 KXC = 580 + (KX - 580) * CAPF
 PAIR_PT = Fraction(BATCH + PAIR_BK + ROUTE, 2) + ELL * KXC + (1 - ELL) * KERNEL
 BABY_PT = F * PAIR_PT + (1 - F) * ROOT + MARG * (REC_GEN + REC_RELAY + REC_INS)
 GIANT_PT = F * PAIR_PT + (1 - F) * ROOT + MARG * (3 * (REC_GEN + REC_RELAY + REC_LOOK) - 2) + BETA
 per_point = Fraction(NB * BABY_PT + NG * GIANT_PT, NB + NG)
 LOOKS = 3 * NG
 BASES_B = ceil(Fraction(NB, EMIN)) + 1; BASES = BASES_B + ceil(Fraction(NG, EMIN)) + 1
 BCLASS = 1 << 42
 assert BASES < 1 << 62 and BASES < d * BCLASS
 assert 6 * BASES_B * (BASES - BASES_B) * (1 << 110) < m
 mu_f = Fraction(3 * NB * NG, 1 << 192)
 assert mu_f < CAND - (1 << 20) and mu_f / (CAND - (1 << 20) - mu_f) ** 2 < Fraction(1, 1 << 40)
 GPC = GENS // d
 BPG = ceil(Fraction(BCLASS, GPC))
 SEGLEN = SEG + TREEMAX
 LANE_OPS = (EMAX // SEG // LANES + 1) * SEGLEN
 assert 2 * WSTEP * LANES == (1 << 23) * (1 << TBITS)
 WAVES_SIDE = ceil(Fraction(BPG * LANE_OPS, WSTEP)) + 2
 PLAN_ROUNDS = ceil(Fraction(EMAX * 16, 1 << 48))
 RMAX, RW = BPG, ceil(Fraction(LANE_OPS, WSTEP // 4)) + 10
 WAVES_G = RMAX * RW
 WAVES = WAVES_SIDE + WAVES_G + PLAN_ROUNDS + 64
 PADOPS = BASES * LANES * SEGLEN + 2 * LANES * LANE_OPS
 fixed_call_bound = (BASES + 2) * (8 * 48 * LANE_OPS + 64)
 assert fixed_call_bound % 2 == 0
 assert 32 + 1 + 5 + 2 == 40 and 5 < 64
 ALLOW = 1
 rows = {}
 rows["records: N_b + N_g points"] = ceil(NB * BABY_PT + NG * GIANT_PT)
 rows["padding ops at the pair price"] = ceil(PADOPS * (PAIR + 2 * BETA))
 rows["side-output fixed calls"] = fixed_call_bound
 rows["leaf count cap (s37): one-wave slack, barrier sums"] = GENS * WSTEP * LANES * 2 * 10 + WAVES * GENS * 16
 rows["listing allowance (1/point)"] = (NB + NG) * ALLOW
 rows["owner pool-cap wave checks"] = OWNERS * (5 * WAVES + 6)
 THRESH = POOL - (1 << SBITS) * CAP2
 assert (1 << SBITS) * CAP2 == (1 << 41) and THRESH + (1 << 41) == POOL
 assert THRESH == sum(((THRESH >> s) & 0xFFFF) << s for s in (0, 16, 32, 48)) and THRESH < (1 << 57)
 assert HEADS * 8 > 0 and HEADS * 8 + POOL * 32 <= (1 << 63)
 PB = HEADS * 8; TP = PB + 32 * THRESH
 assert all(sum(1 for sh in (0, 16, 32, 48) if (v >> sh) & 0xFFFF) <= 2 for v in (PB + 32, TP)) and HEADS == 1 << (128 - 69)
 assert TP + 32 * (1 << 41) == PB + 32 * POOL == 1 << 63 and PB == HEADS * 8
 mu_own, dOwn = Fraction(NB, OWNERS), Fraction(THRESH * OWNERS - NB, NB)
 assert 0 < dOwn < 1 and mu_own * dOwn * dOwn / 3 > 70.5
 assert SMX == 630722 and SM77 == 2473663 and SMX < SM77
 rows["segment restarts"] = (ceil(Fraction(NB + NG, SEG)) + BASES * LANES) * (SMX + 16)
 TL, Mt = [], 3
 while Mt < GMAX:
  pr = min(Mt - 2, -(-(GMAX - Mt - 1) // 2)); TL.append(pr); Mt = Mt + 2 * pr + 1
 assert len(TL) == 30 and Mt >= GMAX and all(1 <= x for x in TL)
 TROUNDS = sum(-(-x // LANES) + 1 for x in TL)
 rows["gap tables"] = BASES * (ceil(TROUNDS * (LANES * PAIR + 128)) + GMAX * 21)
 rows["table strand starts and base points"] = BASES * (LANES + 1) * SM
 G_o, E_o = Fraction(LOOKS, OWNERS), Fraction(LOOKS * NB, OWNERS * OWNERS * HEADS)
 dC, dH = Fraction(1, 10**7), Fraction(97, 10**7)
 lam2 = NG * 3 / Fraction(1 << 120)
 assert lam2 < 1 and 300 * (math.log(300) - 1 - math.log(lam2)) > 2000 * math.log(2)
 assert G_o * dC * dC / 3 > 180 and 2 * (dH * E_o) ** 2 / ((G_o + 1800) * 64 * 64) > 76 and 2 * dC + dH < Fraction(1, 10**5) and 600 < dC * G_o / 100
 assert 76 > 70.5
 cmp_cap = ceil(Fraction(100001 * LOOKS * NB, 100000 * HEADS * OWNERS * OWNERS))
 assert (HEADS * 8) % 16 == 0
 rows["comparisons at the per-owner cap"] = OWNERS * cmp_cap * CMP_EACH
 rows["comparison cap enforcement"] = OWNERS * cmp_cap * CAP_EACH
 OC = OWNERS * HEADS
 y = Fraction(NB, OC - 1)
 assert y < 1 and y ** 3 / 6 < y ** 2 / 2 < y
 p_occ = y - y ** 2 / 2 + y ** 3 / 6
 mu_L = Fraction(LOOKS, OWNERS)
 assert p_occ < Fraction(322, 1000) and mu_L > 10 ** 16 and p_occ > Fraction(9, 100)
 lam = Fraction(NB, OC)
 assert lam - lam * lam / 2 > Fraction(9, 100) and HEADS * 9 // 100 > 10 ** 16
 cap_ne = ceil(Fraction(1000031 * mu_L * p_occ, 1000000))
 dS = Fraction(1, 10**5)
 assert (1 + dS) ** 3 < Fraction(1000031, 1000000) and cap_ne * 3 * OWNERS < LOOKS
 assert mu_L * dS * dS / 3 > 70.5 and HEADS * Fraction(9, 100) * dS * dS / 3 > 70.5
 assert Fraction(cap_ne - CAP2) >= (1 + dS) ** 3 * mu_L * p_occ
 assert CAP2 == (0x200 << 12) and 0x200 < (1 << 12)
 assert 8 + 8 + 8 <= 32 and cap_ne == sum(((cap_ne >> s) & 0xFFFF) << s for s in (0, 16, 32, 48))
 rows["nonempty-lookup extracts and cap"] = OWNERS * cap_ne * EXT_EACH
 rows["owner EXT-cap buffer checks"] = OWNERS * (WAVES * (2 + 2 * (1 << SBITS)) + 5)
 rows["owner initialisation"] = OWNERS * (HEADS * 8 // 64 * 2 + 64)
 LINKS1, LINKS2 = GENS * (1 << TBITS), RELAYS * (1 << UBITS)
 COPY = d * ((1 << 6) + (1 << 16) + (1 << 26))
 NODES = OWNERS + RELAYS + GENS + 2 * d + COPY + (1 << 14)
 bar, BNODES = barrier_instr(NODES, FANIN, WAVES)
 rows["per-wave link words, memory barriers"] = WAVES * ((LINKS1 + LINKS2) * LINK_WAVE + NODES * 24)
 rows["global barriers"] = bar
 rows["plan distribution (copy tree)"] = d * EMAX * ((1 << 6) + (1 << 16) + (1 << 26) + (1 << 27)) * 4
 rows["plan distribution (copy tree)"] += GENS * (EMAX // SEG + 2) * 4
 WIN = 2 * (2 * GMAX // d)
 KAPPA = GMAX // (16 * d * d); assert KAPPA == 8
 rows["planner (scans, matching, emission)"] = d * EMAX * (6 * WIN * 16 + 5 * 64 * KAPPA * (2 * DEPTH * 3 + 20) + 64 + 16)
 rows["planner descendant traversals (adjacency lists)"] = d * EMAX * (64 * TREEMAX + 8 * DEPTH)
 AUG_IN = d * EMAX
 assert (CC, CCY) == (32, 8) and OPS_BOUND == 2 * CC + CC * (6 + CC + CCY * (5 + CC)) == 10752
 rows["planner augmentation: partner cache"] = 2 * AUG_IN * WIN * 16
 rows["planner augmentation: tentative operations"] = 2 * AUG_IN * OPS_BOUND * (16 * DEPTH + 64 + 16)
 rows["planner subtree-size maintenance and pair-fraction check"] = (d * EMAX // 2) * (4 * DEPTH + 4) + d * EMAX * 4
 rows["partial-plan remainders at the root price"] = ceil(2 * EMAX * (ROOT - PAIR_PT))
 rows["character sieve, class lists, their transfer"] = X2 * (400000 + 16 + 20) + d * (1 << 33) * 400000
 rows["base interval sieve, characters, assignment"] = ((1 << 70) + (1 << 43)) * 16 + (1 << 65) * (400000 + 64) + BASES * 64
 rows["candidate delivery and replies"] = ceil(CAND * (8 * REC_RELAY + 200))
 rows["candidate re-walks"] = CAND * 2 * (SEGLEN * (3300 + 232864) + (GMAX - 1) * (2000 + 50) + 2 * SM)
 rows["candidate tests"] = CAND * (96 * SM + (1 << 16) * 20)
 rows["window-table initialisation"] = (LINKS1 + LINKS2) * 8
 rows["chain-id formation"] = (ceil(Fraction(NB + NG, SEG)) + BASES * LANES) * 4
 rows["candidate-link wave counts"] = 4 * (LINKS1 + LINKS2) * 9 * RMAX
 assert RMAX == RMAXC and audit_round_ctl() == 6 and 8 + 10 <= RC_FIX
 rows["round control and stop flags"] = (NODES + BNODES) * RC_FIX * RMAX
 rows["candidate formation"] = CAND * 16
 rows["transport unroll setup and tails"] = UNROLL_FIX * ((LINKS1 + LINKS2) * WAVES + 2 * OWNERS * WAVES)
 rows["beta tables: build (s38)"] = GENS * (batch_audit.TABLE_ENTRIES * len(batch_audit.BUILD) + 7 * 2000 + 64)
 rows["root recovery and input distribution"] = (1 << 24) * 2000 + NODES * 64 + (1 << 40)
 batch = audit_record_path(rows, N, per_point, BABY_PT, GIANT_PT)
 total = sum(rows.values())
 lo_run = success_floor_lower(K); Kg = gate(lo_run)
 MAXR = GENS * EMAX
 EG = ceil(NGT * ((1 - exp4(KT)[0]) / KT + DELTA_P + FAILP)) + MAXR + 2 * EMAX
 GIANT_WC = GIANT_PT + F * ELL * (KX - KXC)
 assert GIANT_WC > GIANT_PT and KX == 590
 assert MAXR + 2 * EMAX < EG < NGT and MAXR * 5000 < NGT and RMAX * RW + WAVES_SIDE < 1 << 23
 e_run = total - rows["records: N_b + N_g points"] + ceil(NB * BABY_PT + EG * GIANT_WC)
 E_num = GATE_C * (1 << 64) + Kg * e_run
 W = Fraction(E_num, 2000 << 64)
 W_run = Fraction(total + GATE_C, 2000)
 baby_point = F_PRE * PAIR_PT + (1 - F_PRE) * ROOT + MARG * (REC_GEN + REC_RELAY + REC_INS)
 QB = max(Fraction(NB, NB + NG), Fraction(BASES_B, BASES), Fraction(BASES_B + 1, BASES + 2))
 pre = (ceil(NB * baby_point) + ceil(QB * rows["padding ops at the pair price"])
       + ceil(QB * rows["side-output fixed calls"])
       + ceil(QB * rows["listing allowance (1/point)"])
       + rows["owner pool-cap wave checks"]
       + rows["owner EXT-cap buffer checks"]
       + ceil(QB * rows["segment restarts"]) + ceil(QB * rows["gap tables"])
       + ceil(QB * rows["table strand starts and base points"]) + rows["owner initialisation"]
       + ceil(Fraction(rows["per-wave link words, memory barriers"], 2)) + ceil(Fraction(rows["global barriers"], 2))
       + rows["plan distribution (copy tree)"] + rows["planner (scans, matching, emission)"]
       + rows["planner descendant traversals (adjacency lists)"]
       + rows["planner augmentation: partner cache"] + rows["planner augmentation: tentative operations"]
       + rows["planner subtree-size maintenance and pair-fraction check"]
       + rows["partial-plan remainders at the root price"] // 2
       + rows["character sieve, class lists, their transfer"] + rows["base interval sieve, characters, assignment"]
       + rows["window-table initialisation"] + ceil(QB * rows["chain-id formation"])
       + rows["beta tables: build (s38)"])
 assert pre < total and 0 < QB < Fraction(3, 5) and 2 * (WAVES_SIDE + PLAN_ROUNDS + 64) < WAVES
 gen_mem = EMAX * 16 + GMAX * 192 + batch_audit.TABLE_BYTES + 512 + LANES * (2 * DEPTH + 3) * 192 + LANES * 16 + (1 << TBITS) * (16 + 2 * CAP1 * RECB) + (1 << TBITS) * CAPC * 64 + BPG * 16 + (EMAX // SEG + 2) * 8 + (LANES + 1) * 272 + (1 << 20) + (1 << 30)
 relay_mem = (GENS >> SBITS) * (16 + 2 * CAP1 * RECB) + (1 << UBITS) * (16 + 2 * CAP2 * RECB) + (1 << UBITS) * CAPC * 64 + (GENS >> SBITS) * CAPC * 64 + (1 << 30)
 owner_mem = HEADS * 8 + POOL * 32 + (1 << SBITS) * (16 + 2 * CAP2 * RECB) + (1 << 16) * 64 + (1 << SBITS) * CAPC * 64 + (1 << 30)
 plan_mem = EMAX * 43 + (1 << 40) + (1 << 36) + (1 << 58) * 8 + (1 << 6) * (1 << 48) * 2 + (1 << 30)
 sieve_mem = (1 << 56) * 4 + (1 << 33) * 4 + d * (1 << 40) + (1 << 30)
 copy_mem = 2 * (1 << 48) + (1 << 10) * (1 << 48) + (1 << 30)
 assert gen_mem < (1 << 64)
 tx = '\n'.join(set('\n'.join(batch_audit.LISTING + batch_audit.beta_listing() + batch_audit.BUILD + smx_audit.LISTING + GEN12 + RELAY_BODY + EXT3 + [INS_U, LOOK_U, LOOKC_U] + [t for q in (0, 1) for t in (kernel_check.kernel_text(q), kernel_check.leaf_text(kernel_check.kernel_text(), q))]).split('\n')))
 assert not re.search(r'\bb\.(lt|ge|gt|le|mi|pl|vs|vc)\b|,\s*(lt|ge|gt|le|mi|pl|vs|vc)\b|\b(ldrs|sxt|asr)', tx) and set(re.findall(r'tbn?z\s+(x\d+),\s*#63,\s*(\.\w+)', tx)) == {('x16', '.Lx_nonneg'), ('x19', '.Lfound'), ('x22', '.Lskip')}
 for v in (relay_mem, plan_mem, sieve_mem, copy_mem):
  assert v < (1 << 63)
 assert (1 << 63) < owner_mem < (1 << 64)
 mem = GENS * gen_mem + RELAYS * relay_mem + OWNERS * owner_mem + d * plan_mem + d * sieve_mem + COPY * copy_mem + (BNODES + (1 << 14)) * (1 << 20)
 with localcontext() as ctx:
  ctx.prec = 90
  log2 = lambda x: (Decimal(x.numerator) / Decimal(x.denominator)).ln() / Decimal(2).ln()
  s = log2(W)
  grid = (s * 10 ** 15).to_integral_value(rounding=ROUND_CEILING)
  assert grid - s * 10 ** 15 > Decimal(10) ** -14
  score = int(grid)
  bound = lambda x, scale: int((log2(x) * scale).to_integral_value(rounding=ROUND_CEILING)) / scale
  succ_exact = Fraction(Kg, 1 << 64) * lo_run
  assert succ_exact >= TARGET and Kg * e_run < E_num and e_run < total
  succ = float(succ_exact)
  return {"K": str(K), "N": N, "N_b": NB, "N_g": NGT, "N_g_run": NG, "expected_giant_points": EG, "expected_run_instructions": e_run, "rounds_max": RMAX, "round_waves": RW, "giant_point_worst_leaf": str(GIANT_WC), "round_mass_max": MAXR, "expected_giant_share": float(Fraction(EG, NGT)), "giant_bases": BASES - BASES_B, "prime_margin": float(Fraction(1 << 70, 53) / BASES), "eps_lit_bound": float(eps_lit()), "baby_point": str(BABY_PT), "giant_point": str(GIANT_PT), "ell": str(ELL), "log2_2N": float(log2(Fraction(2 * N))), "f": str(F), "kernel": KERNEL, "batch_round": batch, "fixed_call_bound": str(fixed_call_bound),
    "pair_op_instructions_two_points": str(PAIR), "root_op_instructions": str(ROOT),
    "record_instructions_listed_avg": str(REC_AVG), "per_point_instructions": str(per_point),
    "per_point_float": float(per_point), "allowance_per_point": ALLOW, "lane_ops": LANE_OPS, "bases": BASES, "waves": WAVES, "nodes_log2": math.log2(NODES),
    "barrier_internal_nodes": BNODES, "rows_log2_units": {k: float(log2(Fraction(v, 2000))) for k, v in rows.items()},
    "rows_instructions": rows,
    "run_total_instructions": total, "expected_numerator_x2pow64": E_num, "pre_instructions_opened_run": pre, "expected_log2": str(s), "score_numerator": score, "score_denominator": 10 ** 15,
    "score": score / 10 ** 15, "opened_run_preprocessing_log2_ceil_1e4": bound(Fraction(pre, 2000), 10 ** 4), "false_match_mean_log2": float(log2(mu_f)), "cmp_cap": cmp_cap, "cap_ne": cap_ne, "h5_conditional_mean_nu": float(Fraction(6 * NB * NG, m) / (1 - exp4(K)[1] - Fraction(5, 10000))), "bases_log2": float(log2(Fraction(BASES))),
    "memory_log2_ceil_1e4": bound(Fraction(mem), 10 ** 4),
    "per_node_bytes_log2": {"gen": math.log2(gen_mem), "relay": math.log2(relay_mem), "owner": math.log2(owner_mem),
          "planner": math.log2(plan_mem), "sieve": math.log2(sieve_mem), "copy": math.log2(copy_mem)},
    "success_floor": succ, "success_floor_exact_lower": str(succ_exact), "success_claimed": 0.3900, "KT": str(KT), "gate_K": Kg, "gate_cost": GATE_C, "E_numerator_x2pow64": E_num, "run_success_lower": str(lo_run), "run_success_float": float(lo_run), "run_instructions": total + GATE_C, "run_log2_units": str(log2(W_run)), "run_score_numerator": int((log2(W_run) * 10 ** 15).to_integral_value(rounding=ROUND_CEILING)), "all_per_node_structures_below_2pow64": True, "generating_node_bytes": gen_mem}
def main():
 cfg = json.loads(sys.stdin.read())
 mode = cfg.get("mode", "all")
 t0 = time.time()
 out = {"mode": mode}
 if mode in ("ledger", "all"):
  out["ledger"] = ledger(cfg)
 if mode in ("planner", "all"):
  out["planner"] = planner_mode(cfg.get("planner", {
   "specs": [{"d": 16, "logX2": 15, "seed": 0, "classes": 8}]}))
 if mode in ("ec", "all"):
  out["ec"] = ec_attack(cfg.get("ec", {
   "bits": 23, "logd": 5, "K": "2.5", "trials": 60, "seed": 2026,
   "keybits": 18, "nodebits": 5, "headbits": 6, "head_cap": 64,
   "cand_cap": 4096, "logX2": 13, "logG": 12, "spotchecks": 6}))
 if mode in ("proto", "all"):
  out["proto"] = proto_mode(cfg.get("proto", {}))
 if mode in ("census", "all"):
  out["census"] = [census(dict(cfg.get("census",
    {"m_bits": 21, "k": 6, "K": "2.5", "seed": 0, "q": "23/50", "nbar": 6})))]
 if mode == "all":
  out["h2"] = h2_sizes()
 out["seconds"] = round(time.time() - t0, 2)
 print(json.dumps(out, indent=1))
if __name__ == "__main__":
 main()
