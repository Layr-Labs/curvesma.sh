#!/usr/bin/env python3 -I
"""Paired class-walk Cheon evidence. Stdin {"mode":...}: ledger, planner, ec, proto, census, all. Stdlib only."""
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
  chi, _ = character(rr, dd)
  classes = [[] for _ in range(dd)]
  for b in range(1, X2 + 1):
    classes[chi(b)].append(b)
  return classes
# planner: planner.py (greedy + s16 augmentation)
from planner import REL32_NONE, REL32_MAX, rel32, plan, OPS_BOUND, CC, CCY
from smx_audit import SMX, audit_smx
# curve
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
def find_curve(bits, dd, rng):
  """Toy curve: 6d | r-1, 2 not in H."""
  while True:
    p = rng.randrange(1 << (bits - 1), 1 << bits) | 1
    if p % 3 != 1 or not is_prime(p):
      continue
    for c in range(1, 8):
      E = Curve(p, c)
      while True:
        x = rng.randrange(p)
        y = sqrt_mod((x ** 3 + c) % p, p)
        if y is not None:
          P = (x, y)
          break
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
# ec toy attack
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
        if cnt >= size:
          baby_n, cnt = len(bases), 0
      elif cnt >= size:
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
  # orbit key: even one of y,p-y, bit 0 dropped
  keyf = lambda Pt: ((Pt[1] if Pt[1] % 2 == 0 else p - Pt[1]) >> 1) & mask
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
      k = keyf(Pt)
      emitted_g.append(v)
      node, head, st = split(k)
      chain = heads.get((node, head), ())
      for (st_b, idv) in chain:
        if st_b == st:
          candidates.append((k, idv, giant_state["id"], v))
      giant_state["cmp"] += len(chain)
    remaining = size
    for bi in range(baby_n):
      if remaining <= 0:
        break
      remaining -= walk(P, bases[bi], baby_emit(bi), remaining, spot)
    assert len(emitted_b) == size
    # H3: occupancy chi-square over node x head cells
    lam_occ = size / (NNODES * NHEADS)
    cells = NNODES * NHEADS
    obs = [len(h) for h in heads.values()]
    chi2 = sum((o - lam_occ) ** 2 for o in obs)
    occ_chi2_sum += chi2 / (size * (1 - 1 / cells))   # scaled by the binomial expectation
    giant_state["id"], giant_state["cmp"] = baby_n, 0
    remaining = size
    for gi in range(baby_n, giant_n):
      if remaining <= 0:
        break
      remaining -= walk(Ad, bases[gi], giant_emit, remaining, spot)
      giant_state["id"] = gi + 1
    assert len(emitted_g) == size
    walk_inv, walk_mul = E.inv, E.muls
    exp_b = []
    for bi in range(baby_n):
      exp_b += expected(bases[bi], size)
      if len(exp_b) >= size:
        break
    exp_b = exp_b[:size]
    exp_g = []
    for gi in range(baby_n, giant_n):
      exp_g += expected(bases[gi], size)
      if len(exp_g) >= size:
        break
    exp_g = exp_g[:size]
    assert exp_b == emitted_b and exp_g == emitted_g
    inverses = [pow(g, -1, rr) for g in exp_g]
    exact_ratios = {u * b * inv % rr for u in U for b in exp_b for inv in inverses}
    hit = t_secret in exact_ratios
    cand_total += len(candidates)
    if len(candidates) > CAND_CAP:
      cand_over += 1
    solved = False
    for (k, idb, idg, gv) in candidates[:CAND_CAP]:
      def collect(base_point, a, key):
        kept = []
        def emit(v, Pt):
          if keyf(Pt) == key:
            kept.append(v)
        walk(base_point, a, emit, size, [0])
        return kept
      keptb = collect(P, bases[idb], k)
      keptg = collect(Ad, bases[idg], k)
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
  Ke = 6 * size * size / mm
  return {
    "p_bits": round(math.log2(p), 2), "r": rr, "d": dd, "X2": X2, "GMAX": GMAX,
    "GMAX_over_d2": round(GMAX / dd / dd, 1), "X2_over_GMAX": round(X2 / GMAX, 1),
    "bases_baby_giant": [baby_n, giant_n - baby_n], "size_per_side": size,
    "plan_pair_fraction": round(pair_frac, 4), "K_eff": round(Ke, 4),
    "predicted": round(1 - math.exp(-Ke), 4),
    "keybits": KEYB, "nodes": NNODES, "heads_per_node": NHEADS,
    "trials": trials, "successes": succ, "success_rate": round(succ / trials, 4),
    "coverage_events_match_exact_each_trial": True,
    "candidates_per_trial": round(cand_total / trials, 3),
    "false_match_model": round(size * size / (1 << KEYB), 3),
    "comparisons_per_trial": round(cmp_total / trials, 1),
    "comparison_model": round(size * size / (NNODES * NHEADS), 1),
    "occupancy_chi2_scaled_mean": round(occ_chi2_sum / trials, 4),
    "rewalk_cap_fails": rewalk_fails, "candidate_cap_over": cand_over,
    "all_recovered_scalars_correct": True,
    "mults_per_element": round(mul_t / trials / (2 * size), 3),
    "rewalk_mults_per_element_at_cand_rate": round(rw_mul_t / trials / (2 * size), 3),
    "reference_mults_per_element": 6.0,
  }
# census
SMOOTH_PRIMES = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
         67, 71, 73, 79, 83, 89, 97)
def census(cfg):
  m_bits, k = cfg["m_bits"], cfg["k"]
  K, seed = float(cfg["K"]), cfg["seed"]
  m_target = 1 << m_bits
  # toy prime r=6c2^k+1, c smooth (primes<=97), (seed+1)-th candidate
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
  N = math.ceil(math.sqrt(K * m / 6))
  e_t = max(2, round(N ** 0.4843))
  X2 = d * e_t
  if cfg.get("window") == "top":
    # control: base window near r//X2 (max wraparound); needs r//X2 > X2
    I0 = r // X2
    assert I0 > X2, "toy too small for top-window control"
  else:
    I0 = X2
  Ilen = int(I0 ** 0.921)
  n_side = math.ceil(N / e_t)
  primes_I = primes_in(I0 + 1, I0 + Ilen)
  need = 2 * (n_side + 8) + 10
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
  def emitted(bs):
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
  sb, cb, nb = emitted(primes_I)
  sg, cg, ng = emitted(primes_I[nb:])
  assert cb == N and cg == N, f"emitted {cb},{cg} != {N}"
  if cfg.get("multiplicity"):
    # H5 diagnostic: per-coset multiplicity nu
    nu = [0] * M
    for ti in sb:
      for uj in sg:
        nu[((ti - uj) % modS) // d] += 1
    marks = bytes(1 if v else 0 for v in nu)
    e1 = sum(nu)
    e2 = sum(v * v for v in nu)
    energy = e2 / (e1 + e1 * e1 / M)   # 1.0, independent uniform pairs
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
  Kr = 6 * N * N / m
  pred = 1 - math.exp(-Kr)
  repfrac = 1 - sum(marks) / (N * N)
  return {
    "r": r, "d": d, "m": m, "window": cfg.get("window", "low"), "K": round(Kr, 6),
    "N": N, "n_side": nb,
    "e_t": e_t, "X2": X2, "coverage": round(cov, 12), "predicted": round(pred, 12),
    "diff": round(cov - pred, 12), "repeat_fraction": round(repfrac, 6),
    "uniform_repeat_fraction": round(1 - (1 - math.exp(-Kr)) / Kr, 5),
    **({"multiplicity": mult} if mult else {}),
  }
# planner mode
def planner_mode(cfg):
  out = []
  for spec in cfg["specs"]:
    dd, logX2, seed = spec["d"], spec["logX2"], spec["seed"]
    rng = random.Random(seed)
    rr = toy_prime(28 + int(math.log2(dd)) + 4, dd, rng)
    X2 = 1 << logX2
    cl = class_lists(rr, dd, X2)
    GMAX = 16 * dd * dd
    rows = []
    for c in range(min(spec.get("classes", 8), dd)):
      L = cl[c]
      if len(L) < 50:
        continue
      _, st = plan(L, GMAX, max_tree=spec.get("max_tree", 1 << 20))
      rows.append(st)
    if rows:
      out.append({
        "r": rr, "d": dd, "X2": X2, "GMAX": GMAX, "GMAX_over_d2": 16,
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
        "min_greedy_fraction": round(min((st["n"] - st["aug_unmatched_in"]) / st["n"] for st in rows), 4),
        "max_unmatched_after_greedy": round(max(st["aug_unmatched_in"] / st["n"] for st in rows), 5),
        "aug_len3_total": sum(st["aug_len3"] for st in rows), "aug_len5_total": sum(st["aug_len5"] for st in rows), "max_ops_per_element": max(st["max_ops_per_element"] for st in rows), "aug_direct_total": sum(st["aug_direct"] for st in rows),
        "aug_candidates_total": sum(st["aug_cands"] for st in rows),
        "c_est": round((1 - sum(st["pair_fraction"] for st in rows) / len(rows)) * dd, 4)
              if spec.get("max_tree", 1 << 20) == dd else None,
      })
  return out
# shared transport formulas
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
WAKE = 256            # interrupt entry, handler dispatch, return (per wake)
def barrier_instr(nodes, fanin, count):
  """Barriers: arrive/wake/sweep/release (worst)."""
  internal, sq = barrier_tree(nodes, fanin)
  arrivals = nodes + internal - 1
  return count * (arrivals * 4 + arrivals * WAKE + sq * 2 + arrivals * (WAKE + 2)), internal
LINK_WAVE = 16        # per link per wave (s9: sender 8 + receiver 3 + 5 margin)
REC_GEN, REC_RELAY, REC_INS, REC_LOOK = 12, 14, 14, 7   # per record; k10 contiguous
CMP_EACH, CAP_EACH, EXT_EACH = 8, 2, 2   # compare (s19: CMP8), cap counter + branch, EXT2 (s21)
GEN12 = ["ubfx","add","ldp","ldp","cmp","b.hs","stp","stp","str","str","ldp","ldr"]
RELAY14 = ["ldp","ldp","ldr","ubfx","add","ldp","cmp","b.hs","stp","stp","str","str","subs","b.ne"]
INS14 = ["ldp","ldp","ldr","lsr","ldr","add","extr","bfi","stp","stp","str","add","subs","b.ne"]
LOOK7 = ["ldp","ldr","lsr","ldr","cbz","subs","b.ne"]
EXT2 = ["subs","extr"]
CMP8 = ["ldr","ldr","ldr","and","cmp","ccmp","b.eq","cbnz"]
K10_LO, HEAD_LSR, INS_BFI, CMP_AND_POS = 60, 6, (54, 10), (54, 10)   # k10 lsb; head lsr; bfi/and pos (s17)
CMP_AND_MASK = 0xFFC0000000000000   # w4 bits 54..63 (s17)
def audit_packed_match():
  M = (1 << 64) - 1
  assert (K10_LO, HEAD_LSR, INS_BFI, CMP_AND_POS) == (60, 6, (54, 10), (54, 10))
  assert CMP_AND_MASK == (((1 << 10) - 1) << 54) == 0xFFC0000000000000   # 10 contiguous ones: ARM64 logical immediate
  def words(k): return (k & M, (k >> 64) & M, (k >> 128) & M)   # w0, w1, w2 of the 192-bit key
  def extr(w1, w0): return (((w1 << 64) | w0) >> K10_LO) & M    # extr full word
  def k10(k): return (k >> K10_LO) & 0x3FF                      # key bits 60..69
  def pk(k): return k10(k) << INS_BFI[0]                        # INS bfi into zeroed-high w4
  def slot(k): return (k >> 70) & ((1 << 58) - 1)               # head bits 70..127
  def ex(k):
    w0, w1, _ = words(k); return extr(w1, w0)                 # EXT full word
  r = random.Random(20261103)
  assert k10(1 << 60) == (ex(1 << 60) & 0x3FF) == 1             # 2^60: k10=1 both sides
  hk = (r.getrandbits(58) << 70) | r.getrandbits(60) | (r.getrandbits(64) << 128)
  assert (words(hk)[1] >> HEAD_LSR) == slot(hk)                # head index is the single lsr (s17)
  kk = [0, (1 << 192) - 1, 1 << 60, ((1 << 192) - 1) ^ (1 << 60), 0x3FF << 60, ((0x3FF << 60) | 0xDEAD),
      1 << 63, 1 << 64, (1 << 60) - 1, (1 << 70) - 1]         # straddles 60/63/64/69/70
  for k in kk + [r.getrandbits(192) for _ in range(1500)]:
    assert (ex(k) & 0x3FF) == k10(k)                         # extr low 10 = oracle k10
    assert pk(k) < 1 << 64 and pk(k) & ((1 << 54) - 1) == 0  # spare-only write
    assert (words(k)[1] >> HEAD_LSR) == slot(k)             # head lsr = oracle
  nb = [x for x in range(192) if x < 60 or x > 69]
  for _ in range(1500):
    k = r.getrandbits(192); g = k10(k)
    chain = (r.getrandbits(31) | pk(k)) & M                  # baby chain word (s17)
    b = chain & CMP_AND_MASK                                 # CMP and-extract
    assert ((b >> 54) & 0x3FF) == g                          # extracted k10 = oracle
    assert b == ((ex(k) << 54) & M)                          # shifted-cmp equal on match
    kn = k ^ (1 << r.randrange(60, 70))
    assert k10(kn) != g and (pk(kn) & CMP_AND_MASK) != b     # k10 flip mismatches
    ko = k ^ (1 << r.choice(nb))
    assert k10(ko) == g and (pk(ko) & CMP_AND_MASK) == b     # other bits inert
    kh = k ^ (1 << r.randrange(70, 128))
    assert slot(kh) != slot(k)                               # head flip changes the chain
def audit_record_path(rows, N, per_point):
  assert per_point == Fraction(1529177, 1000)   # 997/1000 x 2982/2 + 3/1000 x 2050 + 73/2
  assert len(GEN12) == REC_GEN and len(RELAY14) == REC_RELAY and len(INS14) == REC_INS and len(LOOK7) == REC_LOOK and len(EXT2) == EXT_EACH and len(CMP8) == CMP_EACH
  assert 60 + 10 + 58 + 64 == 192 and (127 - 70 + 1) == 58 and (69 - 60 + 1) + (191 - 128 + 1) == 10 + 64   # partition; head 58; cmp 74 (s17)
  for lsb, w in [(0,30),(30,30),(54,10),(54,10)]:   # bfi/and position in range (s17)
    assert 0 <= lsb < 64 and 1 <= w <= 64 and lsb + w <= 64
  assert K10_LO < 64 and HEAD_LSR < 64 and CMP_AND_POS[0] < 64   # shift immediates < 64 (s17)
  assert 16 + 16 + 8 == 40 and 16 + 24 == 40  # = one record
  assert 24 < 256 and 16 <= 504  # imm9/imm7 fit
  assert 64 + 64 + 10 + 62 + 33 <= 256 and 62 + 2 == 64 and 31 + 33 == 64  # pool fits (233/256); chain verbatim (s19)
  assert 33 - 10 >= 0 and (1 << 56) * 32 < 1 << 63       # k10 fits spare [54..63]; pool < 2^63 (s17/s19)
  m = 0xCCCCCCCD                                                          # /5 magic: m*5 == 2^34+1
  assert m * 5 == 2 ** 34 + 1 and m < 2 ** 32 and 2 ** 24 * 5 < 2 ** 32
  assert 4 * (2 ** 34 + 1) + 5 * (2 ** 24 - 1) < 5 * 2 ** 34               # (x*m)>>34 == x//5 for x < 2^24*5
  assert all((x * 2341) >> 14 == x // 7 for x in range(1792))   # candidate /56 magic
  assert 2341 < 2 ** 16 and 1791 * 2341 < 2 ** 32                    # single-movz magic; umull exact
  assert 8 + 3 + 5 == LINK_WAVE and 4 <= 8 and 8 <= 9 and 11 <= 16        # 8+3+margin; fits
  assert "owner pool-cap compares" not in rows and "head-index word combine" not in rows   # deleted by s17/s18
  assert REC_GEN + REC_RELAY + REC_INS == 40   # baby-phase listed
  assert EXT2 == ["subs","extr"] and K10_LO == 60   # giant extr bits 60..69 (s21: b.lo at buffers)
  audit_packed_match()   # packed INS/EXT/CMP semantics
  audit_smx()            # exact restart price: the generated restart listing executed end to end (s15)
# proto (sync wave transport)
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
  buf = {}                       # (sender,receiver,parity) -> records
  tag = {}                       # (sender,receiver,parity) -> (written_wave,read_wave)
  cnt = {}                       # count words: (sender,receiver,parity) -> (count,written_wave)
  ptr = {}  # s9 fill pointers
  single = cfg.get("single_count", False)   # neg control
  emitted, delivered = [], defaultdict(list)
  heads, candidates = defaultdict(list), []
  overflow = []
  def open_buf(sender, receiver, parity, wave):
    key = (sender, receiver, parity)
    lst = buf.setdefault(key, [])
    w, rd = tag.get(key, (None, None))
    if w != wave:              # first use this wave: previous contents read in wave w+1
      assert w is None or rd == w + 1, "buffer reused before read"
      lst.clear(); tag[key] = (wave, None)
    return lst
  def put(sender, receiver, parity, rec, cap, wave):
    lst = open_buf(sender, receiver, parity, wave)
    pw, pv = ptr.get((sender, receiver), (wave, 0))
    if pw != wave: assert pv == 0; pw, pv = wave, 0
    if pv >= cap * 40:                                # pre-write cap check on fill pointer
      overflow.append((sender, receiver, wave)); return False
    lst.append(rec); ptr[(sender, receiver)] = (wave, pv + 40); return True
  def publish(sender, receiver, parity, wave):   # wave end: derive count, store it, reset
    lst = open_buf(sender, receiver, parity, wave)
    pw, pv = ptr.get((sender, receiver), (wave, 0))
    val = pv if pw == wave else 0
    if pw != wave: assert pv == 0
    assert val % 40 == 0 and val // 40 == len(lst)
    ptr[(sender, receiver)] = (wave, 0)
    cnt[(sender, receiver) if single else (sender, receiver, parity)] = (len(lst), wave)
  def take(sender, receiver, parity, wave):      # read count, then records
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
  def owner_of(key): return key & (O - 1)           # t = low TB bits, u = next UB
  def t_of(x): return x & ((1 << TB) - 1)
  def k10_of(key): return (key >> K10_LO) & 0x3FF                       # EXT bits 60..69
  def pack_of(key): return k10_of(key) << INS_BFI[0]                    # INS bfi to w4[54..63]
  def barrier():
    instr, _ = barrier_instr(NODES, FANIN, 1)
    cost["barrier"] += instr
  wave_count = Counter()
  seq = 0
  for phase in (0, 1):                               # 0 = baby, 1 = giant
    for wl in range(WAVES + 2):
      w = phase * (WAVES + 2) + wl               # global wave; wl = phase index
      # generating nodes emit wave w
      if wl < WAVES:
        for g in range(G):
          for lane in range(LANES):
            for st in range(STEPS):
              n_out = 2 if rng.random() < 0.98 else 1
              for _ in range(n_out):
                if pattern == "one_owner": key = 0 | (rng.getrandbits(180) << 12)
                elif pattern == "two_owners": key = (rng.randrange(2) * (O // 2)) | (rng.getrandbits(180) << 12)
                else: key = rng.getrandbits(192)
                rec = (key, phase, g, lane, seq); seq += 1
                emitted.append(rec)
                cost["rec_gen"] += 1
                if not put(g, G + relay_of(key, g), w & 1, rec, CAP1, w):
                  return {"pattern": pattern, "failed": True, "reason": "link cap overflow (run fails)",
                      "overflow_at": overflow[0], "records_before_failure": len(emitted)}
          for rl in range(g & ((1 << SB) - 1), R, 1 << SB):
            publish(g, G + rl, w & 1, w)       # count word written per link, empty or not
          cost["link"] += R >> SB            # count words to 2^TB relays
          wave_count["gen"] += 1
      # relays: forward wave w-1
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
          cost["link"] += (G >> SB) + O // (1 << TB)   # inbound read, outbound write
          wave_count["relay"] += 1
      # owners: consume wave w-2
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
                if heads.get((o, slot)): cost["nonempty"] += 1   # s12: extracts execute here only
                for (rb, chb, bb) in heads[(o, slot)]:
                  cost["compare"] += 1
                  pkm = (((chb >> CMP_AND_POS[0]) & ((1 << CMP_AND_POS[1]) - 1)) == k10_of(rec[0]) and rb == rest)
                  assert pkm == (rb == rest)   # packed decides as full compare
                  if pkm: candidates.append((o, rec[0], bb, rec))
          cost["link"] += 1 << SB
          wave_count["owner"] += 1
      barrier(); cost["barriers"] += 1
  # exact delivery
  got = [rec for o in delivered for rec in delivered[o]]
  assert Counter(got) == Counter(emitted) and len(got) == len(emitted)
  n_b = sum(1 for r in emitted if r[1] == 0); n_g = len(emitted) - n_b
  assert cost["rec_gen"] == cost["rec_relay"] == len(emitted) and cost["insert"] == n_b and cost["lookup"] == n_g
  # candidate path owner->relay->gen, both chains, then replies
  cand_cost = 0; resolved = 0
  # tags: write at w', read at w'+1, rewrite at w'+2
  for (o, key, bb, gg) in candidates:
    t = t_of(o); replies = []
    for side in (bb, gg):
      g = side[2]; rl = (t << SB) | (g & ((1 << SB) - 1))
      cand_cost += 2 * REC_RELAY                 # owner->relay->gen (two hops)
      assert side in emitted                     # gen node finds chain (re-walk charged)
      # reply: gen -> relay R(t_key,s) -> owner(key) via forward links
      rl2 = relay_of(key, g); assert owner_of(key) == o
      cand_cost += 2 * REC_RELAY; replies.append(side[4])
    resolved += 1
  assert resolved == len(candidates)
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
      "instructions_by_item": {"gen": cost["rec_gen"] * REC_GEN, "relay": cost["rec_relay"] * REC_RELAY,
                   "insert": cost["insert"] * REC_INS, "lookup": cost["lookup"] * REC_LOOK,
                   "compare": cost["compare"] * CMP_EACH, "link_words": cost["link"] * LINK_WAVE,
                   "barriers": cost["barrier"], "candidate_path": cand_cost},
      "instructions_total": total, "closed_form_equal": True,
      "scope": "sequential simulation; listing applied to counted events"}
def proto_mode(cfg):
  out = []
  c = dict(cfg); c["pattern"] = "uniform"; out.append(proto(c))
  c = dict(cfg); c.update({"pattern": "two_owners", "cap2": 8 * cfg.get("cap2", 96)}); out.append(proto(c))   # skewed keys, wide caps
  c = dict(cfg); c.update({"pattern": "one_owner"}); out.append(proto(c))   # one-owner overflow fails
  assert out[2]["failed"] and not out[0]["failed"] and not out[1]["failed"]
  c = dict(cfg); c.update({"pattern": "uniform", "single_count": True})   # neg control
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
      # stride assignment (v20): lane j takes segments j, j+LANES, ...
      st = [sum(L[j::LANES]) for j in range(LANES)]
      res.append((kind, S, pad, LANES * max(st) - ops))
  return {"lanes": LANES, "seg": SEG, "treemax": TM, "cases": len(res), "padding_lemma_holds": True,
      "padding_bound": LANES * (SEG + TM),
      "max_greedy_padding": max(r[2] for r in res), "max_stride_padding": max(r[3] for r in res)}
# ledger
KERNEL = 929  # worst lane path of measure/side_lanes.S (kernel_check.py; K1/K2 of proof.md s1; Callgrind mean 922.925 below it)
def exp_neg_bounds(x, terms=40):
  """Exact rational enclosure of exp(-x), 0<x<1: alternating Taylor series with decreasing terms (even/odd partial sums bracket)."""
  assert Fraction(0) < x < Fraction(1)
  t, s, lo, up = Fraction(1), Fraction(0), None, None
  for n in range(terms):
    s += t
    if n % 2 == 0: up = s
    else: lo = s
    t = t * (-x) / (n + 1)
  assert lo < up and up - lo < Fraction(1, 10 ** 50)
  return lo, up
def success_floor_lower(K):
  """Exact lower bound of 1 - exp(-K) - 0.0005 - 2^-20 - 7 x 2^-40 - 2^-111 (s2)."""
  return 1 - exp_neg_bounds(K)[1] - Fraction(5, 10000) - Fraction(1, 1 << 20) - Fraction(7, 1 << 40) - Fraction(1, 1 << 111)
FLOOR_TARGET = Fraction(39, 100) + Fraction(1, 10 ** 10)   # declared 0.39 plus a 1e-10 arithmetic margin
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
  KT, F = Fraction(12377947373, 25 * 10**9), Fraction(997, 1000)    # KT: least 1e-12 multiple with floor >= 0.39+1e-10 (s2); F = per-class floor (s16), runtime-checked
  assert success_floor_lower(KT) >= FLOOR_TARGET > success_floor_lower(KT - Fraction(1, 10**12))
  F_PRE = Fraction(39, 40)  # pre-only floor (below F; preprocessing is not scored)
  assert F_PRE <= F
  N = math.isqrt((KT * m) // 6)                          # least N with 6 N^2 >= KT m
  while 6 * N * N < KT * m: N += 1
  K = Fraction(6 * N * N, m)                             # realized K >= 0.4952
  assert 6 * (N - 1) ** 2 < KT * m <= 6 * N * N and K >= KT
  X2, GMAX, EMIN, EMAX = 1 << 76, 1 << 44, 1 << 55, 1 << 57
  LANES, SEG, TREEMAX, WSTEP, DEPTH = 1 << 18, 1 << 24, 1 << 20, 1 << 34, 1 << 14
  GENS, SBITS, TBITS, UBITS = 1 << 47, 20, 30, 30
  RELAYS, OWNERS = 1 << (TBITS + SBITS), 1 << (TBITS + UBITS)
  HEADS, POOL, CAP1, CAP2, CAPC, CAND = 1 << 58, 1 << 56, 1 << 24, 1 << 21, 1 << 8, 1 << 44
  RECB, FANIN = 40, 1 << 10
  MARG = Fraction(1, 1)  # charged exactly
  SM = 8192000                                          # one scalar mult, instructions
  SM77 = -(-77 * SM // 255)                             # 77-bit price, same per-bit rate (s13)
  assert SM77 * 255 >= 77 * SM and SM77 < SM and 4 * 77 > 255   # ceil exact; below SM; slack covers setup
  inv_share = ceil(Fraction(232864, LANES))             # batched inversion share/lane op
  PAIR_BK = 8 + 7 + 12 + 11 + 7 + inv_share + 0 + 7    # branchy-7 key; control in KERNEL
  ROOT_BK = 5 + 12 + 14 + 6 + inv_share + 5 + 7
  TABLE_BK = 8 + 7 + 6 + inv_share + 5 + 21   # +21: v = p - x, n = p - y per entry (3 ldp, 2 subs + 10 sbcs, 6 stp)
  assert TABLE_BK <= 50                         # re-walk rebuilds tables at 2000 + 50
  REC_AVG = Fraction(2 * (REC_GEN + REC_RELAY) + REC_INS + REC_LOOK, 2)
  PAIR = 2000 + KERNEL + MARG * PAIR_BK                 # two points
  ROOT = 2000 + MARG * ROOT_BK
  per_point = F * PAIR / 2 + (1 - F) * ROOT + MARG * REC_AVG
  # sizes (H2 bounds, sieve-checked at runtime; violation stops run = failure)
  BASES_SIDE = ceil(Fraction(N, EMIN)) + 1
  BASES = 2 * BASES_SIDE
  BCLASS = 1 << 42                                      # one-character bases, both sides
  GPC = GENS // d                                       # gen nodes/class
  BPG = ceil(Fraction(BCLASS, GPC))                     # bases/gen node
  SEGLEN = SEG + TREEMAX                                # entries/plan segment, at most
  LANE_OPS = (EMAX // SEG // LANES + 1) * SEGLEN         # ops/lane/base (ops/LANES + one segment)
  WAVES_SIDE = ceil(Fraction(BPG * LANE_OPS, WSTEP)) + 2
  PLAN_ROUNDS = ceil(Fraction(EMAX * 16, 1 << 48))
  WAVES = 2 * WAVES_SIDE + PLAN_ROUNDS + 64
  PADOPS = BASES * LANES * SEGLEN + 2 * LANES * LANE_OPS   # per-lane pad bound
  fixed_call_bound = 40 * (BASES + 2) * LANE_OPS + 64 * (BASES + 2)
  assert fixed_call_bound == (BASES + 2) * (40 * LANE_OPS + 64) and fixed_call_bound % 2 == 0   # form; baby/online split exact
  assert 32 + 1 + 5 + 2 == 40 and 5 < 64   # s20: 32+1+5+2=40; terminator lsl #5
  ALLOW = 1  # per-point unlisted allowance
  rows = {}
  rows["records: 2N points"] = ceil(2 * N * per_point)
  rows["padding ops at the pair price"] = ceil(PADOPS * PAIR)
  rows["side-output fixed calls"] = fixed_call_bound        # 40/invocation + 64/base (s20)
  rows["listing allowance (1/point)"] = 2 * N * ALLOW
  rows["owner pool-cap wave checks"] = OWNERS * (5 * WAVES + 6)   # s18: wave-start check + wave-end spill; +6/owner init
  THRESH = POOL - (1 << SBITS) * CAP2
  assert (1 << SBITS) * CAP2 == (1 << 41) and THRESH + (1 << 41) == POOL   # one wave <= 2^41/owner; start + wave = cap (s18)
  assert THRESH == sum(((THRESH >> s) & 0xFFFF) << s for s in (0, 16, 32, 48)) and THRESH < (1 << 56)   # movz+3xmovk init (s18)
  assert HEADS * 8 > 0 and HEADS * 8 + POOL * 32 < (1 << 63)   # s19: pool after heads, nonzero base; addrs < 2^63
  mu_own, dOwn = Fraction(N, OWNERS), Fraction(THRESH * OWNERS - N, N)
  assert 0 < dOwn < 1 and mu_own * dOwn * dOwn / 3 > 70.5   # s18 Chernoff vs 70.5 nats (same H3 union member)
  assert SMX == 630722 and SM77 == 2473663 and SMX < SM77
  rows["segment restarts"] = (ceil(Fraction(2 * N, SEG)) + BASES * LANES) * (SMX + 16)   # s15 exact price from the shipped kernel bodies; +dir/cursor/handoff
  rows["gap tables"] = BASES * (GMAX - 1) * (2000 + ceil(MARG * TABLE_BK))
  rows["table strand starts and base points"] = BASES * (LANES + 1) * SM
  # s11 1.00001: counts within 1e-7 (Chernoff >= 180 nats) + Hoeffding on chains <= 64 (2^-200 tail) at 9.7e-6 (>= 94 nats)
  G_o, E_o = Fraction(N, OWNERS), Fraction(N * N, OWNERS * OWNERS * HEADS)
  dC, dH = Fraction(1, 10**7), Fraction(97, 10**7)
  assert G_o * dC * dC / 3 > 180 and 2 * (dH * E_o) ** 2 / (G_o * 64 * 64) > 94 and 2 * dC + dH < Fraction(1, 10**5)
  assert 94 > 70.5   # 2^-100 per owner over 2^60 owners stays below 2^-40 (70.5 nats); both exponents exceed it
  cmp_cap = ceil(Fraction(100001 * N * N, 100000 * HEADS * OWNERS * OWNERS))
  rows["comparisons at the per-owner cap"] = OWNERS * cmp_cap * CMP_EACH
  rows["comparison cap enforcement"] = OWNERS * cmp_cap * CAP_EACH  # online only (s9)
  OC = OWNERS * HEADS                                    # occupancy cells; keys uniform (H3)
  y = Fraction(N, OC - 1)
  assert y < 1 and y ** 3 / 6 < y ** 2 / 2 < y            # Taylor terms decrease (s12 bound valid)
  p_occ = y - y ** 2 / 2 + y ** 3 / 6                    # P(nonempty) <= p_occ (s12)
  mu_L = Fraction(N, OWNERS)                             # exact mean lookups/inserts/owner
  assert p_occ < Fraction(322, 1000) and mu_L > 10 ** 16 and p_occ > Fraction(17, 100)   # value + floors (s12)
  lam = Fraction(N, OC)
  assert lam - lam * lam / 2 > Fraction(17, 100) and HEADS * 17 // 100 > 10 ** 16   # occupied floor (s12)
  cap_ne = ceil(Fraction(1000031 * mu_L * p_occ, 1000000))       # 1.000031 covers three 1.00001 stages (s12)
  dS = Fraction(1, 10**5)
  assert (1 + dS) ** 3 < Fraction(1000031, 1000000) and cap_ne * 2 * OWNERS < N   # covered; EXT beats 2/lookup
  assert mu_L * dS * dS / 3 > 70.5 and HEADS * Fraction(17, 100) * dS * dS / 3 > 70.5   # per-stage Chernoff exponents at the asserted floors
  assert Fraction(cap_ne - CAP2) >= (1 + dS) ** 3 * mu_L * p_occ   # s21: buffer threshold keeps the 3-stage tail
  assert CAP2 == (0x200 << 12) and 0x200 < (1 << 12)   # s21: cmp-imm 2^21 (12-bit << 12)
  assert 8 + 8 + 8 <= 32 and cap_ne == sum(((cap_ne >> s) & 0xFFFF) << s for s in (0, 16, 32, 48))   # s21: slot-0 fit; movz+3xmovk init
  rows["nonempty-lookup extracts and cap"] = OWNERS * cap_ne * EXT_EACH   # online only (s21: EXT2)
  rows["owner EXT-cap buffer checks"] = OWNERS * (WAVES * (2 + 2 * (1 << SBITS)) + 5)   # s21: cmp+b.lo/buffer; ldr+str/wave; +5/owner
  rows["owner initialisation"] = OWNERS * (HEADS * 8 // 16 + 64)
  LINKS1, LINKS2 = GENS * (1 << TBITS), RELAYS * (1 << UBITS)
  COPY = d * ((1 << 6) + (1 << 16) + (1 << 26))
  NODES = OWNERS + RELAYS + GENS + 2 * d + COPY + (1 << 14)     # 2^14 spare nodes
  bar, BNODES = barrier_instr(NODES, FANIN, WAVES)
  rows["per-wave link words, memory barriers"] = WAVES * ((LINKS1 + LINKS2) * LINK_WAVE + NODES * 24)
  rows["global barriers"] = bar
  rows["plan distribution (copy tree)"] = d * EMAX * ((1 << 6) + (1 << 16) + (1 << 26) + (1 << 27)) * 4
  rows["plan distribution (copy tree)"] += GENS * (EMAX // SEG + 2) * 4         # per-node segment directory
  WIN = 2 * (2 * GMAX // d)
  rows["planner (scans, matching, emission)"] = d * EMAX * (3 * WIN * 16 + 64 * (2 * DEPTH * 3 + 20) + 64 + 16)
  # traversals: <=4l max_tree visits at 16 + l DEPTH steps at 8
  rows["planner descendant traversals (adjacency lists)"] = d * EMAX * (64 * TREEMAX + 8 * DEPTH)
  # s16 augmentation: one list-order pass (lengths 1/3/5) charged for every element: one window scan fills the partner
  # cache, <= OPS_BOUND tentative ops at one ancestor walk each (depth cap) + candidate reads.
  AUG_IN = d * EMAX
  assert (CC, CCY) == (32, 8) and OPS_BOUND == 2 * CC + CC * (6 + CC + CCY * (5 + CC)) == 10752
  rows["planner augmentation: partner cache"] = AUG_IN * WIN * 16
  rows["planner augmentation: tentative operations"] = AUG_IN * OPS_BOUND * (16 * DEPTH + 64 + 16)
  rows["planner subtree-size maintenance and pair-fraction check"] = (d * EMAX // 2) * (4 * DEPTH + 4) + d * EMAX * 2
  rows["partial-plan remainders at the root price"] = ceil(2 * EMAX * (ROOT - PAIR / 2))   # the one partial class per side, every point as a root
  rows["character sieve, class lists, their transfer"] = X2 * (400000 + 16 + 20) + d * (1 << 33) * 400000
  rows["base interval sieve, characters, assignment"] = ((1 << 70) + (1 << 43)) * 16 + (1 << 65) * (400000 + 64) + BASES * 64
  rows["candidate delivery and replies"] = CAND * (8 * REC_RELAY + 200)
  rows["candidate re-walks"] = CAND * 2 * (SEGLEN * (3300 + 232864) + (GMAX - 1) * (2000 + 50) + 2 * SM)   # rebuild at 2000 + 50 >= 2000 + TABLE_BK
  rows["candidate tests"] = CAND * (96 * SM + (1 << 16) * 20)
  rows["window-table initialisation"] = (LINKS1 + LINKS2) * 8   # s9: 4+4 margin
  rows["chain-id formation"] = (ceil(Fraction(2 * N, SEG)) + BASES * LANES) * 4   # s10: 3+stp
  rows["candidate-link wave counts"] = 4 * (LINKS1 + LINKS2) * 9   # 9 on 4 traversals (s9)
  rows["candidate formation"] = CAND * 16                          # 10 used/match (proof.md s9)
  rows["root recovery and input distribution"] = (1 << 24) * 2000 + NODES * 64 + (1 << 40)
  audit_record_path(rows, N, per_point)
  total = sum(rows.values())
  W = Fraction(total, 2000)
  # Preprocessing = full baby phase: baby 40/record (giant 33),
  # walk part at F_PRE (s4).
  baby_point = F_PRE * PAIR / 2 + (1 - F_PRE) * ROOT + MARG * (REC_GEN + REC_RELAY + REC_INS)
  pre = (ceil(N * baby_point) + ceil(Fraction(rows["padding ops at the pair price"], 2))
       + rows["side-output fixed calls"] // 2
       + ceil(Fraction(rows["listing allowance (1/point)"], 2))
       + rows["owner pool-cap wave checks"]   # full baby phase (giant checks conservatively incl.)
       + rows["owner EXT-cap buffer checks"]   # s21: full (baby-wave checks conservatively incl.)
       + ceil(Fraction(rows["segment restarts"], 2)) + ceil(Fraction(rows["gap tables"], 2))
       + ceil(Fraction(rows["table strand starts and base points"], 2)) + rows["owner initialisation"]
       + ceil(Fraction(rows["per-wave link words, memory barriers"], 2)) + ceil(Fraction(rows["global barriers"], 2))
       + rows["plan distribution (copy tree)"] + rows["planner (scans, matching, emission)"]
       + rows["planner descendant traversals (adjacency lists)"]
       + rows["planner augmentation: partner cache"] + rows["planner augmentation: tentative operations"]
       + rows["planner subtree-size maintenance and pair-fraction check"]
       + rows["partial-plan remainders at the root price"] // 2
       + rows["character sieve, class lists, their transfer"] + rows["base interval sieve, characters, assignment"]
       + rows["window-table initialisation"] + ceil(Fraction(rows["chain-id formation"], 2)))
  assert pre < total
  # memory: every per-node structure + total below 2^63
  gen_mem = EMAX * 16 + GMAX * 192 + LANES * (2 * DEPTH + 1) * 128 + LANES * 16 + (1 << TBITS) * (16 + 2 * CAP1 * RECB) + (1 << TBITS) * CAPC * 64 + BPG * 16 + (EMAX // SEG + 2) * 8 + (1 << 30)   # s8 slots; gap 192 (v, n stored); s10 chain slots
  relay_mem = (GENS >> SBITS) * (16 + 2 * CAP1 * RECB) + (1 << UBITS) * (16 + 2 * CAP2 * RECB) + (1 << UBITS) * CAPC * 64 + (GENS >> SBITS) * CAPC * 64 + (1 << 30)
  owner_mem = HEADS * 8 + POOL * 32 + (1 << SBITS) * (16 + 2 * CAP2 * RECB) + (1 << 16) * 64 + (1 << SBITS) * CAPC * 64 + (1 << 30)
  plan_mem = EMAX * 32 + EMAX * 4 + (1 << 37) + (1 << 58) * 8 + (1 << 6) * (1 << 48) * 2 + (1 << 30)   # arrays+adjacency; s16 child prev link, partner cache
  sieve_mem = (1 << 56) * 4 + (1 << 33) * 4 + d * (1 << 40) + (1 << 30)
  copy_mem = 2 * (1 << 48) + (1 << 10) * (1 << 48) + (1 << 30)
  for v in (gen_mem, relay_mem, owner_mem, plan_mem, sieve_mem, copy_mem):
    assert v < (1 << 63)
  mem = GENS * gen_mem + RELAYS * relay_mem + OWNERS * owner_mem + d * plan_mem + d * sieve_mem + COPY * copy_mem + (BNODES + (1 << 14)) * (1 << 20)
  # success: 1-e^-K minus H1 allowance, modelling budget, cap tails, diagonal
  with localcontext() as ctx:
    ctx.prec = 90
    log2 = lambda x: (Decimal(x.numerator) / Decimal(x.denominator)).ln() / Decimal(2).ln()
    s = log2(W)
    grid = (s * 10 ** 15).to_integral_value(rounding=ROUND_CEILING)
    assert grid - s * 10 ** 15 > Decimal(10) ** -14
    score = int(grid)
    bound = lambda x, scale: int((log2(x) * scale).to_integral_value(rounding=ROUND_CEILING)) / scale
    succ_exact = success_floor_lower(K)               # exact enclosure at the realized K >= KT (s2)
    assert succ_exact >= FLOOR_TARGET
    succ = float(succ_exact)
    return {"K": str(K), "N": N, "log2_2N": float(log2(Fraction(2 * N))), "f": str(F), "kernel": KERNEL, "fixed_call_bound": str(fixed_call_bound),
        "pair_op_instructions_two_points": str(PAIR), "root_op_instructions": str(ROOT),
        "record_instructions_listed_avg": str(REC_AVG), "per_point_instructions": str(per_point),
        "per_point_float": float(per_point), "allowance_per_point": ALLOW, "lane_ops": LANE_OPS, "bases": BASES, "waves": WAVES, "nodes_log2": math.log2(NODES),
        "barrier_internal_nodes": BNODES, "rows_log2_units": {k: float(log2(Fraction(v, 2000))) for k, v in rows.items()},
        "rows_instructions": rows,
        "W_instructions": total, "pre_instructions": pre, "W_log2": str(s), "score_numerator": score, "score_denominator": 10 ** 15,
        "score": score / 10 ** 15, "preprocessing_log2_ceil_1e4": bound(Fraction(pre, 2000), 10 ** 4),
        "memory_log2_ceil_1e4": bound(Fraction(mem), 10 ** 4),
        "per_node_bytes_log2": {"gen": math.log2(gen_mem), "relay": math.log2(relay_mem), "owner": math.log2(owner_mem),
                    "planner": math.log2(plan_mem), "sieve": math.log2(sieve_mem), "copy": math.log2(copy_mem)},
        "success_floor": succ, "success_floor_exact_lower": str(succ_exact), "success_claimed": 0.3900, "KT": str(KT), "all_per_node_structures_below_2pow63": True}
# main
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
      "bits": 23, "logd": 5, "K": "0.50", "trials": 60, "seed": 2026,
      "keybits": 18, "nodebits": 5, "headbits": 6, "head_cap": 64,
      "cand_cap": 4096, "logX2": 13, "logG": 12, "spotchecks": 6}))
  if mode in ("proto", "all"):
    out["proto"] = proto_mode(cfg.get("proto", {}))
  if mode in ("census", "all"):
    out["census"] = [census(dict(cfg.get("census",
        {"m_bits": 21, "k": 4, "K": "0.50", "seed": 0, "multiplicity": True})))]
  out["seconds"] = round(time.time() - t0, 2)
  print(json.dumps(out, indent=1))
if __name__ == "__main__":
  main()
