#!/usr/bin/env python3 -I
"""Paired class-walk evidence: ledger/planner/ec/proto/census. Stdlib only."""
import json,math,random,sys,time
from fractions import Fraction
R_PRODUCTION = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
D_PRODUCTION = 1 << 20
M_PRODUCTION = (R_PRODUCTION - 1) // D_PRODUCTION
def is_prime(n):
 if n < 2:
  return False
 for q in (2,3,5,7,11,13,17,19,23,29,31,37):
  if n % q == 0:
   return n == q
 e,s = n - 1,0
 while e % 2 == 0:
  e //= 2
  s += 1
 for a in (2,3,5,7,11,13,17,19,23,29,31,37):
  x = pow(a,e,n)
  if x in (1,n - 1):
   continue
  for _ in range(s - 1):
   x = x * x % n
   if x == n - 1:
    break
  else:
   return False
 return True
def primes_in(lo,hi):
 return [q for q in range(lo,hi + 1) if is_prime(q)]
def sqrt_mod(a,p):
 a %= p
 if a == 0:
  return 0
 if pow(a,(p - 1) // 2,p) != 1:
  return None
 q,s = p - 1,0
 while q % 2 == 0:
  q //= 2
  s += 1
 z = next(z for z in range(2,p) if pow(z,(p - 1) // 2,p) == p - 1)
 mm,c,t,x = s,pow(z,q,p),pow(a,q,p),pow(a,(q + 1) // 2,p)
 while t != 1:
  i,t2 = 0,t
  while t2 != 1:
   t2 = t2 * t2 % p
   i += 1
  b = pow(c,1 << (mm - i - 1),p)
  mm,c,t,x = i,b * b % p,t * b * b % p,x * b % p
 return x
def toy_prime(bits,dd,rng):
 while True:
  rr = rng.randrange(1 << (bits - 1),1 << bits) // (6 * dd) * (6 * dd) + 1
  if is_prime(rr) and ((rr - 1) // dd) % 6 == 0 and pow(2,(rr - 1) // dd,rr) != 1:
   return rr
def character(rr,dd):
 """chi: integers -> Z/d."""
 mm = (rr - 1) // dd
 gen = next(g for g in range(2,10 ** 6)
 if all(pow(g,(rr - 1) // q,rr) != 1 for q in (2,3)))
 zeta = pow(gen,mm,rr)
 zlog = {pow(zeta,i,rr): i for i in range(dd)}
 return (lambda x: zlog[pow(x,mm,rr)]),gen
def class_lists(rr,dd,X2):
 chi,_ = character(rr,dd); f,c = [0] * (X2 + 1),[0] * (X2 + 1)
 for q in range(2,math.isqrt(X2) + 1): f[q * q::q] = [q] * ((X2 - q * q) // q + 1)
 for b in range(2,X2 + 1): c[b] = (c[f[b]] + c[b // f[b]]) % dd if f[b] else chi(b)
 assert all(c[b] == chi(b) for b in range(1,X2 + 1,97))
 classes = [[] for _ in range(dd)]
 for b in range(1,X2 + 1): classes[c[b]].append(b)
 return classes
def h2_sizes(seed=3):
 """s35: toy class-length spread."""
 rng = random.Random(seed); out = []
 for lg,lx in ((4,14),(4,18),(6,18),(8,18),(8,20)):
  dd = 1 << lg; L = [len(c) for c in class_lists(toy_prime(40,dd,rng),dd,1 << lx)]
  out.append([dd,lx - lg,round(min(L) * dd / 2 ** lx,4),round(max(L) * dd / 2 ** lx,4)])
 return out
from planner import REL32_NONE,REL32_MAX,rel32,plan,OPS_BOUND,CC,CCY
from smx_audit import SMX,audit_smx
import batch_audit,kernel_check,smx_audit,re
from batch_audit import BATCH,ROUND_FIXED,RMAXC,RC_FIX,audit_round_ctl,eps_lit
class Curve:
 def __init__(self,p,c):
  self.p,self.c = p,c
  self.inv = self.muls = 0
 def add(self,P,Q):
  if P is None:
   return Q
  if Q is None:
   return P
  p = self.p
  (x1,y1),(x2,y2) = P,Q
  if x1 == x2:
   if (y1 + y2) % p == 0:
    return None
   lam = 3 * x1 * x1 * pow(2 * y1,-1,p) % p
  else:
   lam = (y2 - y1) * pow(x2 - x1,-1,p) % p
  x3 = (lam * lam - x1 - x2) % p
  return (x3,(lam * (x1 - x3) - y1) % p)
 def single(self,P,T):
  """Root step P+T."""
  assert P is not None and T is not None and P[0] != T[0]
  self.inv += 1
  self.muls += 6
  return self.add(P,T)
 def pair(self,P,T):
  """P+T,P-T from one inversion."""
  p = self.p
  (x1,y1),(x2,y2) = P,T
  assert x1 != x2
  inv = pow(x2 - x1,-1,p)
  self.inv += 1
  self.muls += 9
  lp = (y2 - y1) * inv % p
  lm = (-y2 - y1) * inv % p
  x3p = (lp * lp - x1 - x2) % p
  x3m = (lm * lm - x1 - x2) % p
  return ((x3p,(lp * (x1 - x3p) - y1) % p),
 (x3m,(lm * (x1 - x3m) - y1) % p))
 def mul(self,k,P):
  R = None
  while k:
   if k & 1:
    R = self.add(R,P)
   P = self.add(P,P)
   k >>= 1
  return R
def cm_ok(p,dd):
 """Speed only: skip p if no j=0 order p+1-t can pass."""
 x,y = p,sqrt_mod(-3,p)
 while y * y > p: x,y = y,x % y
 b = math.isqrt((p - y * y) // 3)
 return y * y + 3 * b * b != p or any(is_prime(N) and (N - 1) % (6 * dd) == 0 and pow(2,(N - 1) // dd,N) != 1 for t in (2 * y,y + 3 * b,y - 3 * b) for N in (p + 1 - t,p + 1 + t))
def find_curve(bits,dd,rng):
 """Toy curve: 6d | r-1, 2 not in H."""
 while True:
  p = rng.randrange(1 << (bits - 1),1 << bits) | 1
  if p % 3 != 1 or not is_prime(p):
   continue
  ok = cm_ok(p,dd)
  for c in range(1,8):
   E = Curve(p,c)
   while True:
    x = rng.randrange(p)
    y = sqrt_mod((x ** 3 + c) % p,p)
    if y is not None:
     P = (x,y)
     break
   if not ok:
    continue
   lo,hi = p + 1 - 2 * math.isqrt(p) - 2,p + 1 + 2 * math.isqrt(p) + 2
   w = math.isqrt(hi - lo) + 1
   baby,R = {},None
   for j in range(w):
    if R is not None:
     baby.setdefault(R,j)
    R = E.add(R,P)
   step = E.mul(w,P)
   step = (step[0],(-step[1]) % p) if step else None
   Gk,cands = E.mul(lo,P),[]
   for i in range(w + 1):
    if Gk is None:
     cands.append(lo + i * w)
    elif (Gk[0],(-Gk[1]) % p) in baby:
     cands.append(lo + i * w + baby[(Gk[0],(-Gk[1]) % p)])
    Gk = E.add(Gk,step)
   cands = sorted({N for N in cands if lo <= N <= hi and E.mul(N,P) is None})
   if len(cands) == 1 and is_prime(cands[0]) and (cands[0] - 1) % (6 * dd) == 0 \
     and pow(2,(cands[0] - 1) // dd,cands[0]) != 1:
    return E,cands[0],P
def dth_root(t,dd,rr):
 k = dd.bit_length() - 1
 for j in range(k,0,-1):
  y = sqrt_mod(t,rr)
  e = (rr - 1) >> (j - 1)
  t = next(z for z in (y,rr - y) if j == 1 or pow(z,e,rr) == 1)
 return t
def ec_attack(cfg):
 bits,dd,K = cfg["bits"],1 << cfg["logd"],float(cfg["K"])
 trials,seed = cfg["trials"],cfg["seed"]
 KEYB,NODEB,HEADB = cfg["keybits"],cfg["nodebits"],cfg["headbits"]
 STOREB = KEYB - NODEB - HEADB
 rng = random.Random(seed)
 E,rr,P = find_curve(bits,dd,rng)
 p,mm = E.p,(rr - 1) // dd
 beta = next(z for z in (pow(g,(p - 1) // 3,p) for g in range(2,100)) if z != 1)
 lamg = next(L for L in (pow(g,(rr - 1) // 3,rr) for g in range(2,100))
 if L != 1 and E.mul(L,P) == (beta * P[0] % p,P[1]))
 U = [1,rr - 1,lamg,rr - lamg,lamg * lamg % rr,rr - lamg * lamg % rr]
 assert all(pow(u,mm,rr) == 1 for u in U)
 chi,gen = character(rr,dd)
 size = math.isqrt(int(K * mm / 6)) + 1
 Qr = Fraction(cfg.get('q','1')); size_b = size_g = size
 if Qr != 1:
  size_g = 1
  while 6 * (-(-(Qr.numerator * size_g) // Qr.denominator)) * size_g < K * mm: size_g += 1
  size_b = -(-(Qr.numerator * size_g) // Qr.denominator); size = max(size_b,size_g)
 X2,GMAX = 1 << cfg["logX2"],1 << cfg["logG"]
 classes = [[] for _ in range(dd)]
 for b in range(1,X2 + 1):
  classes[chi(b)].append(b)
 plans = [plan(L,GMAX) if L else (None,None) for L in classes]
 pair_frac = (sum(2 * st["pairs"] for _,st in plans if st) /
 sum(st["n"] for _,st in plans if st))
 X1 = X2 + int(X2 ** 0.921)
 while True:
  base_primes = primes_in(X2 + 1,X1)
  bases,cnt,baby_n,giant_n = [],0,None,None
  for a in base_primes:
   bases.append(a)
   cnt += len(classes[(-chi(a)) % dd])
   if baby_n is None:
    if cnt >= size_b:
     baby_n,cnt = len(bases),0
   elif cnt >= size_g:
    giant_n = len(bases)
    break
  if baby_n and giant_n:
   break
  X1 *= 2
 assert baby_n and giant_n
 def expected(a,limit):
  vals,c = [],(-chi(a)) % dd
  L,(ops,_) = classes[c],plans[c]
  for op in ops:
   positions = (op[1],) if op[0] in (0,1) else (op[3],op[2])
   for position in positions:
    vals.append(a * L[position])
    if len(vals) == limit:
     return vals
  return vals
 def walk(base_point,a,emit,limit,checks):
  produced,c = 0,(-chi(a)) % dd
  L,(ops,_) = classes[c],plans[c]
  Qa = E.mul(a,base_point)
  tab = {}
  def T(g):
   if g not in tab:
    tab[g] = E.mul(g,Qa)
   return tab[g]
  pts,root = {},None
  for op in ops:
   if op[0] == 0:
    root = E.mul(L[op[1]],Qa)
    pts[op[1]] = root
    emit(a * L[op[1]],root)
    produced += 1
   elif op[0] == 1:
    root = E.single(root,T(op[2]))
    pts[op[1]] = root
    emit(a * L[op[1]],root)
    produced += 1
   else:
    _,x,i,j,g = op
    Pp,Pm = E.pair(pts[x],T(g))
    pts[j],pts[i] = Pp,Pm
    emit(a * L[j],Pp)
    produced += 1
    if produced < limit:
     emit(a * L[i],Pm)
     produced += 1
   if checks[0] > 0 and produced % 41 == 0:
    v = a * L[op[1] if op[0] in (0,1) else op[3]]
    assert E.mul(v,base_point) == pts[op[1] if op[0] in (0,1) else op[3]]
    checks[0] -= 1
   if produced >= limit:
    return produced
  return produced
 NNODES,NHEADS = 1 << NODEB,1 << HEADB
 HEAD_CAP,KEPT_CAP,CAND_CAP = cfg["head_cap"],4,cfg["cand_cap"]
 mask = (1 << KEYB) - 1
 keyf = lambda Pt: (Pt[0] >> 1) & mask
 keys3 = lambda Pt: [(Pt[0] >> 1) & mask,(beta * Pt[0] % p >> 1) & mask,(beta * beta * Pt[0] % p >> 1) & mask]
 split = lambda k: (k >> (KEYB - NODEB),(k >> STOREB) & (NHEADS - 1),
 k & ((1 << STOREB) - 1))
 succ = inv_t = mul_t = cand_total = cmp_total = rw_mul_t = 0
 rewalk_fails = cand_over = 0
 occ_chi2_sum = 0.0
 spot = [cfg["spotchecks"]]
 order,cb_,cg_,ib,ig = [],0,0,0,baby_n
 while ib < baby_n or ig < giant_n:
  if ig >= giant_n or (ib < baby_n and cb_ * size_g <= cg_ * size_b):
   order.append((0,ib)); cb_ += len(expected(bases[ib],size_b)); ib += 1
  else:
   order.append((1,ig)); cg_ += len(expected(bases[ig],size_g)); ig += 1
 stops = early = 0
 for _ in range(trials):
  s = rng.randrange(1,rr)
  t_secret = pow(s,dd,rr)
  A1,Ad = E.mul(s,P),E.mul(t_secret,P)
  E.inv = E.muls = 0
  heads = {}
  candidates = []
  emitted_b,emitted_g = [],[]
  state = {"cmp": 0}
  def emit_side(side,idv):
   def emit(v,Pt):
    (emitted_g if side else emitted_b).append(v)
    for k in (keys3(Pt) if side else [keyf(Pt)]):
     node,head,st = split(k)
     chain = heads.get((1 - side,node,head),())
     for (st_o,ido) in chain:
      if st_o == st:
       candidates.append((k,ido,idv,v) if side else (k,idv,ido,v))
     state["cmp"] += len(chain)
     h = heads.setdefault((side,node,head),[])
     h.append((st,idv))
     assert len(h) <= HEAD_CAP
   return emit
  def resolve(cands):
   nonlocal rewalk_fails
   for (k,idb,idg,gv) in cands:
    def collect(base_point,a,key,kf,lim):
     kept = []
     def emit(v,Pt):
      if key in kf(Pt):
       kept.append(v)
     walk(base_point,a,emit,lim,[0])
     return kept
    keptb = collect(P,bases[idb],k,lambda Pt: [keyf(Pt)],size_b)
    keptg = collect(Ad,bases[idg],k,keys3,size_g)
    if len(keptb) > KEPT_CAP or len(keptg) > KEPT_CAP:
     rewalk_fails += 1
     continue
    for bv in keptb:
     for gvv in keptg:
      for u in U:
       z = u * bv * pow(gvv,-1,rr) % rr
       if E.mul(z,P) == Ad:
        s0 = dth_root(z,dd,rr)
        w = pow(gen,(rr - 1) // dd,rr)
        sol = next(s0 * pow(w,i,rr) % rr for i in range(dd)
 if E.mul(s0 * pow(w,i,rr) % rr,P) == A1)
        assert sol == s,"wrong scalar"
        return True
   return False
  rem,solved,stop,seen = [size_b,size_g],False,None,0
  for rnd,(side,bi) in enumerate(order):
   if rem[side] > 0:
    rem[side] -= walk(Ad if side else P,bases[bi],emit_side(side,bi),rem[side],spot)
   new = candidates[seen:]; seen = len(candidates)
   if len(candidates) <= CAND_CAP and resolve(new):
    solved,stop = True,rnd
    break
  walk_inv,walk_mul = E.inv,E.muls
  vb,vg,val_r = {},[],{}
  for rnd,(side,bi) in enumerate(order):
   for v in expected(bases[bi],size_b if side == 0 else size_g):
    if side == 0: vb.setdefault(v % rr,rnd)
    else: vg.append((v,rnd))
  first = None
  for (gv,rg) in vg:
   for u in U:
    rb = vb.get(t_secret * gv * pow(u,-1,rr) % rr)
    if rb is not None and (first is None or max(rb,rg) < first): first = max(rb,rg)
  pre_b = [v for (side,bi) in order[:len(order) if stop is None else stop + 1] if side == 0 for v in expected(bases[bi],size_b)]
  assert emitted_b == pre_b[:size_b]
  hit = first is not None
  cand_total += len(candidates)
  if len(candidates) > CAND_CAP:
   cand_over += 1
  assert solved == hit and stop == first,f"resolution mismatch solved={solved} stop={stop} exact={first}"
  stops += solved; early += solved and stop < len(order) - 1
  lam_occ = len(emitted_b) / (NNODES * NHEADS)
  cells = NNODES * NHEADS
  obs = [len(h) for (sd,_,_),h in heads.items() if sd == 0]
  occ_chi2_sum += sum((o - lam_occ) ** 2 for o in obs) / (max(1,len(emitted_b)) * (1 - 1 / cells))
  cmp_total += state["cmp"]
  inv_t += walk_inv
  mul_t += walk_mul
  rw_mul_t += E.muls - walk_mul
  if solved:
   succ += 1
 Ke = 6 * size_b * size_g / mm
 return {
 "p_bits": round(math.log2(p),2),"r": rr,"d": dd,"X2": X2,"GMAX": GMAX,
 "GMAX_over_d2": round(GMAX / dd / dd,1),"X2_over_GMAX": round(X2 / GMAX,1),
 "bases_baby_giant": [baby_n,giant_n - baby_n],"sizes_baby_giant": [size_b,size_g],"q": str(Qr),"keys": "x; giants x, beta x, beta^2 x",
 "plan_pair_fraction": round(pair_frac,4),"K_eff": round(Ke,4),
 "predicted": round(1 - math.exp(-Ke),4),
 "keybits": KEYB,"nodes": NNODES,"heads_per_node": NHEADS,
 "trials": trials,"successes": succ,"success_rate": round(succ / trials,4),
 "coverage_events_match_exact_each_trial": True,"rounds": len(order),"stopped_before_last_round": early,"stop_round_is_first_covering_round_each_trial": True,
 "candidates_per_trial": round(cand_total / trials,3),
 "false_match_model": round(3 * size_b * size_g / (1 << KEYB),3),
 "comparisons_per_trial": round(cmp_total / trials,1),
 "comparison_model": round(3 * size_b * size_g / (NNODES * NHEADS),1),
 "occupancy_chi2_scaled_mean": round(occ_chi2_sum / trials,4),
 "rewalk_cap_fails": rewalk_fails,"candidate_cap_over": cand_over,
 "all_recovered_scalars_correct": True,
 "mults_per_element": round(mul_t / trials / (size_b + size_g),3),
 "rewalk_mults_per_element_at_cand_rate": round(rw_mul_t / trials / (size_b + size_g),3),
 "reference_mults_per_element": 6.0,
 }
SMOOTH_PRIMES = (5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,
 67,71,73,79,83,89,97)
def mark(marks,sbd,u,M):
 for t in sbd: marks[(t - u) % M] = 1
def census(cfg):
 m_bits,k = cfg["m_bits"],cfg["k"]
 K,seed = float(cfg["K"]),cfg["seed"]
 m_target = 1 << m_bits
 c_lo = m_target // 6
 smooths = []
 def gen(idx,val):
  if idx == len(SMOOTH_PRIMES):
   smooths.append(val)
   return
  v = val
  while v <= 4 * c_lo:
   gen(idx + 1,v)
   v *= SMOOTH_PRIMES[idx]
 gen(0,1)
 smooths = sorted({s for s in smooths if s >= c_lo})
 count = 0
 r = None
 for c in smooths:
  rr_try = 6 * c * (1 << k) + 1
  if is_prime(rr_try) and ((rr_try - 1) // (1 << k)) % 6 == 0 \
    and pow(2,(rr_try - 1) // (1 << k),rr_try) != 1:
   if count == seed:
    r = rr_try
    break
   count += 1
 assert r is not None,"no toy prime for this seed"
 d = 1 << k
 m = (r - 1) // d
 n_all = r - 1
 fac = {}
 for q in (2,3):
  while n_all % q == 0:
   fac[q] = fac.get(q,0) + 1
   n_all //= q
 rest,q = n_all,5
 while q * q <= rest:
  while rest % q == 0:
   fac[q] = fac.get(q,0) + 1
   rest //= q
  q += 2
 if rest > 1:
  fac[rest] = fac.get(rest,0) + 1
 gen = None
 for cand in range(2,100000):
  if all(pow(cand,(r - 1) // q,r) != 1 for q in fac) and pow(cand,r - 1,r) == 1:
   gen = cand
   break
 assert gen is not None
 parts = []
 for q,e in fac.items():
  pe = q ** e
  beta = pow(gen,(r - 1) // pe,r)
  gam = pow(beta,q ** (e - 1),r)
  b = int(math.isqrt(q)) + 1
  baby = {}
  cur = 1
  for j in range(b):
   baby[cur] = j
   cur = cur * gam % r
  gamminv = pow(pow(gam,b,r),-1,r)
  parts.append((q,e,pe,beta,gam,b,baby,gamminv))
 def dlog_q(t,b,baby,gamminv):
  cur = t
  for i in range(b + 1):
   j = baby.get(cur)
   if j is not None:
    return i * b + j
   cur = cur * gamminv % r
  raise ValueError("dlog_q")
 def dlog(x):
  crt,mods = [],[]
  for q,e,pe,beta,gam,b,baby,gamminv in parts:
   xq = pow(x,(r - 1) // pe,r)
   digit = 0
   for i in range(e):
    t = pow(xq,q ** (e - 1 - i),r)
    di = dlog_q(t,b,baby,gamminv)
    digit += di * q ** i
    xq = xq * pow(beta,(pe - (di * q ** i) % pe) % pe,r) % r
   crt.append(digit)
   mods.append(pe)
  X,Nm = 0,1
  for a,mo in zip(crt,mods):
   t = (a - X) % mo
   X = X + Nm * (t * pow(Nm,-1,mo) % mo)
   Nm *= mo
  return X % (r - 1)
 assert dlog(2) % d != 0,"2 in H_t"
 Qr = Fraction(cfg.get('q','1')); NB = NG = math.ceil(math.sqrt(K * m / 6))
 if Qr != 1:
  NG = math.isqrt(int(K * m / 6 / Qr)) // 2
  while 6 * (-(-(Qr.numerator * NG) // Qr.denominator)) * NG < K * m: NG += 1
  NB = -(-(Qr.numerator * NG) // Qr.denominator)
 N = max(NB,NG)
 nbar = cfg.get('nbar')
 e_t = max(2,round(NG / (nbar * d))) if nbar else max(2,round(N ** float(cfg.get('e_exp','0.4811'))))
 X2 = d * e_t
 if cfg.get("window") == "top":
  I0 = r // X2
  assert I0 > X2,"toy too small for top-window control"
 else:
  I0 = X2
 Ilen = int(I0 ** 0.921)
 n_side = math.ceil(N / e_t)
 primes_I = primes_in(I0 + 1,I0 + Ilen)
 need = (NB + NG) * 13 // (10 * e_t) + 50 if nbar else 2 * (n_side + 8) + 10
 while len(primes_I) < need:
  Ilen *= 2
  primes_I = primes_in(I0 + 1,I0 + Ilen)
 spf = list(range(X2 + 1))
 for q in range(2,int(math.isqrt(X2)) + 1):
  if spf[q] == q:
   for mm2 in range(q * q,X2 + 1,q):
    if spf[mm2] == mm2:
     spf[mm2] = q
 Ep = {}
 for q in range(2,X2 + 1):
  if spf[q] == q:
   Ep[q] = dlog(q)
 modn = r - 1
 E = [0] * (X2 + 1)
 for b in range(2,X2 + 1):
  q = spf[b]
  E[b] = (E[b // q] + Ep[q]) % modn if q != b else Ep[q] % modn
 classes = {c2: [] for c2 in range(d)}
 for b in range(1,X2 + 1):
  classes[E[b] % d].append(b)
 modS = (r - 1) // 6
 M = m // 6
 def emitted(bs,N):
  vals,cnt,used = [],0,0
  for a in bs:
   if cnt >= N:
    break
   used += 1
   Ea = dlog(a)
   L = classes.get((-Ea) % d) or []
   take = min(len(L),N - cnt)
   vals.extend((Ea + E[b]) % modS for b in L[:take])
   cnt += take
  return vals,cnt,used
 sb,cb,nb = emitted(primes_I,NB)
 if nbar:
  gpc,i = cfg.get("gpc",1),0
  def side(N):
   nonlocal i
   out,tot = [],0
   while tot < N:
    Ea = dlog(primes_I[i]); i += 1; L = classes[(-Ea) % d]
    out.append((Ea % d,[(Ea + E[b]) % modS // d for b in L])); tot += len(L)
   return out,tot
  bb,totb = side(NB); gb,totg = side(NG)
  def rounds(bs):
   byc = {}
   for j,g in enumerate(bs): byc.setdefault(g[0],[]).append(j)
   R = -(-max(map(len,byc.values())) // gpc)
   return [[byc[c][x] for c in sorted(byc) for x in range(j * gpc,min(len(byc[c]),(j + 1) * gpc))] for j in range(R)]
  rb,rg = rounds(bb),rounds(gb)
  def sweep(kinds):
   marks,out,B,G,ib,ig = bytearray(M),[],[],[],0,0
   for kd in kinds:
    if kd == 0:
     new = [u for j in rb[ib] for u in bb[j][1]]; ib += 1
     ng_neg = [-g for g in G]
     for u in new: mark(marks,ng_neg,-u,M)
     B += new
    else:
     new = [u for j in rg[ig] for u in gb[j][1]]; ig += 1
     for u in new: mark(marks,B,u,M)
     G += new
    Kj = 6 * len(B) * len(G) / m
    out.append([kd,round(Kj,4),round(marks.count(1) / M - 1 + math.exp(-Kj),6)])
   return out
  kinds,ib,ig,cb_,cg_ = [],0,0,0,0
  while ib < len(rb) or ig < len(rg):
   if ig >= len(rg) or (ib < len(rb) and cb_ * totg <= cg_ * totb):
    kinds.append(0); cb_ += sum(len(bb[j][1]) for j in rb[ib]); ib += 1
   else:
    kinds.append(1); cg_ += sum(len(gb[j][1]) for j in rg[ig]); ig += 1
  rng = random.Random(seed)
  ctl = [0] * len(rb) + [1] * len(rg); rng.shuffle(ctl)
  tp,dp = 0,0
  for _ in range(cfg.get("lit",500)):
   A,Bc = classes[(-rng.choice(bb)[0]) % d],classes[(-rng.choice(gb)[0]) % d]
   tp += len(A) * len(Bc); dp += len({(x // math.gcd(x,y),y // math.gcd(x,y)) for x in A for y in Bc})
  eps,a,b = 1 - dp / tp,sweep(kinds),sweep(ctl)
  lw = lambda o: [round(min(x[2] for x in o if x[1] > 0),6),round(min(x[2] + eps * x[1] * math.exp(-x[1]) for x in o if x[1] > 0),6)]
  return {"r": r,"d": d,"m": m,"X2": X2,"N": [totb,totg],"rounds": len(kinds),"eps_lit": round(eps,6),"production_order": a,"shuffled_order": b,
 "lowest": [lw(a)[0],lw(b)[0]],"lowest_plus_literal": [lw(a)[1],lw(b)[1]]}
 sg,cg,ng = emitted(primes_I[nb:],NG)
 assert (cb,cg) == (NB,NG)
 if cfg.get("multiplicity"):
  nu = [0] * M
  for ti in sb:
   for uj in sg:
    nu[((ti - uj) % modS) // d] += 1
  marks = bytes(1 if v else 0 for v in nu)
  e1 = sum(nu)
  e2 = sum(v * v for v in nu)
  energy = e2 / (e1 + e1 * e1 / M)
  mult = {"energy_ratio": round(energy,4),
 "mean_nu_marked": round(e1 / max(1,sum(marks)),4),
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
 "r": r,"d": d,"m": m,"window": cfg.get("window","low"),"K": round(Kr,6),
 "N": [NB,NG],"q": str(Qr),"n_side": nb,
 "e_t": e_t,"X2": X2,"coverage": round(cov,12),"predicted": round(pred,12),
 "diff": round(cov - pred,12),"repeat_fraction": round(repfrac,6),
 "uniform_repeat_fraction": round(1 - (1 - math.exp(-Kr)) / Kr,5),
 **({"multiplicity": mult} if mult else {}),
 }
def planner_mode(cfg):
 out = []
 for spec in cfg["specs"]:
  dd,logX2,seed = spec["d"],spec["logX2"],spec["seed"]
  rng = random.Random(seed)
  rr = toy_prime(28 + int(math.log2(dd)) + 4,dd,rng)
  X2 = 1 << logX2
  cl = class_lists(rr,dd,X2)
  GMAX = 16 * spec.get("kappa",1) * dd * dd
  rows = []
  for c in range(min(spec.get("classes",8),dd)):
   L = cl[c]
   if len(L) < 50:
    continue
   _,st = plan(L,GMAX,max_tree=spec.get("max_tree",1 << 20))
   rows.append(st)
  if rows:
   out.append({
 "r": rr,"d": dd,"X2": X2,"GMAX": GMAX,"GMAX_over_d2": GMAX // (dd * dd),
 "span_over_2GMAX": round(X2 / (2 * GMAX),1),
 "classes": len(rows),
 "list_lengths": [min(st["n"] for st in rows),max(st["n"] for st in rows)],
 "min_pair_fraction": round(min(st["pair_fraction"] for st in rows),4),
 "mean_pair_fraction": round(
 sum(st["pair_fraction"] for st in rows) / len(rows),4),
 "max_cycle_rejections": max(st["cycle_rejections"] for st in rows),
 "max_depth_seen": max(st["max_depth"] for st in rows),
 "depth_cap": 1 << 14,"max_tree_size_seen": max(st["max_tree_size"] for st in rows),
 "size_rejections": sum(st["size_rejections"] for st in rows),
 "tree_cap": spec.get("max_tree",1 << 20),
 "roots_per_class": [st["roots"] for st in rows],
 "min_leaf_fraction": round(min(st["leaf_fraction"] for st in rows),5),
 "mean_leaf_fraction": round(sum(st["leaf_fraction"] for st in rows) / len(rows),5),
 "centres_hosting": [st["centres_hosting"] for st in rows],"fallback_new_centres": [st["fallback_new_centres"] for st in rows],
 "max_ops_per_element": max(st["max_ops_per_element"] for st in rows),"ops_bound": OPS_BOUND,
 })
 return out
def barrier_tree(nodes,fanin):
 """Fan-in tree: (internal nodes, sum children^2)."""
 internal = sq = 0
 level = nodes
 while level > 1:
  parents = -(-level // fanin)
  full,rest = divmod(level,fanin)
  sq += full * fanin * fanin + (rest * rest if rest else 0)
  internal += parents
  level = parents
 return internal,sq
WAKE = 256
def barrier_instr(nodes,fanin,count):
 """Barriers: arrive/wake/sweep/release (worst)."""
 internal,sq = barrier_tree(nodes,fanin)
 arrivals = nodes + internal - 1
 return count * (arrivals * 4 + arrivals * WAKE + sq * 2 + arrivals * (WAKE + 2)),internal
LINK_WAVE = 16
UNT,UNK = 64,6
REC_GEN = 12; REC_RELAY,REC_IL = (Fraction(UNT * b + 2,UNT) for b in (12,12))
CMP_EACH,CAP_EACH = 7,2
GEN12 = ["ubfx","add","ldp","ldp","cmp","b.hs","stp","stp","str","str","ldp","ldr"]
RELAY_BODY = ["ldp","ldp","ldr","ubfx","add","ldp","cmp","b.hs","stp","stp","str","str"]
RELAY_U = RELAY_BODY * UNT + ["subs","b.ne"]
UNROLL_SETUP = ["and","lsr","cbz"]
UNROLL_FIX = len(UNROLL_SETUP) + UNK + 1
def unroll(name,body,U=UNT):
 """s29: U-fold loop (one subs + b.ne per U records) + binary tail ladder; label SKIP renamed per copy."""
 k = U.bit_length() - 1; assert U == 1 << k >= 2; c = iter(range(4 * U)); b = lambda: body.replace("SKIP","S%d" % next(c))
 L = [name + ":"," and x5, x23, #%d" % (U - 1)," lsr x23, x23, #%d" % k," cbz x23, tail","pair:"] + [b() for _ in range(U)] + [" subs x23, x23, #1"," b.ne pair","tail:"]
 for j in range(k - 1,-1,-1): L += [" tbz x5, #%d, T%d" % (j,j)] + [b() for _ in range(1 << j)] + ["T%d:" % j]
 return "\n".join(L + [" ret"])
K10_LO,HEAD_LSR,INS_BFI,CMP_AND_POS = 62,7,(55,9),(55,9)
HEAD_BITS = 128 - 64 - HEAD_LSR
CMP_AND_MASK = 0xFF80000000000000
M64 = (1 << 64) - 1
def cmp_field(c4,x10):
 """The displayed field compare: and b,c4,#CMP_AND_MASK ; cmp b,x10,lsl #55 (x10 = extr of the record key, unmasked)."""
 return (c4 & CMP_AND_MASK) == (x10 << CMP_AND_POS[0]) & M64
def cmp_copy(base,link,lab):
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
SKIP%sN:""" % (lab,link,base,base,lab,base,lab)
CMP_A = cmp_copy("x11","x13","A") + "\n cbz x13, SKIP"
CMP_B = cmp_copy("x13","x11","B") + "\n cbnz x11, SKIPA"
CMP7 = [l.split()[0] for l in CMP_A.split("\n") if l[0] == " "][2:8] + ["cbz"]
def il_body(own):
 """s45 owner body: insert into the record's own table (head word +own of the head pair), look up the other table."""
 a,b = ("x13","x11") if own == 0 else ("x11","x13")
 return """ ldp x0, x1, [x20], #16
 ldp x2, x3, [x20], #16
 ldr x4, [x20], #8
 lsr x9, x1, #%d
 add x12, x21, x9, lsl #4
 ldp %s, %s, [x12]
 extr x10, x1, x0, #%d
 bfi x4, x10, #55, #9
 stp x3, x4, [x22, #16]
 str x22, [x12, #%d]
 stp x2, x13, [x22], #32
 cbz x11, SKIP
""" % (HEAD_LSR,a,b,K10_LO,own) + CMP_A + "\n" + CMP_B + "\nSKIP:"
IL_B,IL_G = il_body(0),il_body(8)
IL12 = [l.split()[0] for l in IL_B.split("SKIPA:")[0].split("\n") if l.startswith(" ")]
FAIL = "\nfail:\n mov x27, #1\n ret"
IL_BU,IL_GU = unroll("ilb",IL_B) + FAIL,unroll("ilg",IL_G) + FAIL
def key_parts(k):
 return k & M64,k >> 64 & M64,k >> 128
def audit_il(bodies=None,seed=20261116):
 """s45: executed kind-homogeneous rounds of the owner loops against the full-key oracle (heads, pools, candidate lists);
 every equal-key (baby, giant key) pair is reported exactly once, by whichever record arrives later."""
 import a64
 from collections import defaultdict
 BU,GU = bodies or (IL_BU,IL_GU)
 pB,pG = a64.a64_load(BU),a64.a64_load(GU)
 r = random.Random(seed); HB = HEAD_LSR + 64; PB,PG,T,CL = (1 << 61) + 32,(1 << 62) + 32,1 << 40,1 << 39
 def key(f9,head,k64,own=0): return own | f9 << K10_LO | head << HB | k64 << 128
 pool = [key(0,1,0),key(0,1,7),key(511,1,7),key(256,2,5),key(1,3,1 << 63),key(3,4,9,5)]
 pool += [key(r.getrandbits(9),r.randrange(10),r.getrandbits(64),r.getrandbits(62)) for _ in range(40)]
 def variant(k):
  return r.choice([k,k,k ^ 1 << r.randrange(K10_LO,K10_LO + 9),k ^ 1 << r.randrange(128,192),k ^ 1 << HB,k ^ r.getrandbits(K10_LO)])
 kinds = [0,1,0,0,1,1,0,1,1,0,1,0]
 mem = defaultdict(int); chains = defaultdict(list); ptr = [PB,PG]; seen = [[],[]]; found = []; steps = 0; C_all = 0
 rec_of = {}
 for rnd,kind in enumerate(kinds):
  n = r.choice([1,0,63,64,65,70]) if rnd < 6 else r.randrange(20,140)
  ks = [variant(r.choice(pool)) for _ in range(n)]
  base = T + (rnd << 20)
  want,C = [],0
  for i,k in enumerate(ks):
   w0,w1,w2 = key_parts(k); w3,w4 = r.getrandbits(64),r.getrandbits(55)
   for j,w in enumerate((w0,w1,w2,w3,w4)): mem[base + 40 * i + 8 * j] = w
   h = w1 >> HEAD_LSR
   other = chains[(h,1 - kind)]
   C += len(other)
   for (bk,a) in other:
    if bk >> K10_LO == k >> K10_LO: want.append(a); found.append((rec_of[a],(kind,rnd,i)) if kind else ((kind,rnd,i),rec_of[a]))
   a = ptr[kind] + 32 * i; rec_of[a] = (kind,rnd,i)
   chains[(h,kind)].insert(0,(k,a)); seen[kind].append(((kind,rnd,i),k))
  snap = dict(mem); o = {}
  st = a64.a64_run(*(pB if kind == 0 else pG),"ilb" if kind == 0 else "ilg",{"x20": base,"x21": 0,"x22": ptr[kind],"x23": n,"x25": C,"x26": CL},mem,out_regs=o)
  got = [mem[CL + 8 * i] for i in range((o["x26"] - CL) // 8)]
  assert got == want and o["x27"] == 0 and o["x25"] == 0 and o["x22"] == ptr[kind] + 32 * n
  assert st == 12 * n + 2 * (n // UNT) + (CMP_EACH + CAP_EACH) * C + len(want) + UNROLL_FIX
  for i,k in enumerate(ks):
   w0,w1,w2 = key_parts(k); a = ptr[kind] + 32 * i; h = w1 >> HEAD_LSR
   assert mem[a] == w2 and mem[a + 16] == snap[base + 40 * i + 24] and mem[a + 24] == snap[base + 40 * i + 32] | (k >> K10_LO & 511) << 55
  for (h,kd),ch in chains.items(): assert not ch or mem[16 * h + 8 * kd] == ch[0][1] and all(mem[ch[x][1] + 8] == (ch[x + 1][1] if x + 1 < len(ch) else 0) for x in range(len(ch)))
  if C:
   o2 = {}; m2 = defaultdict(int,snap)
   a64.a64_run(*(pB if kind == 0 else pG),"ilb" if kind == 0 else "ilg",{"x20": base,"x21": 0,"x22": ptr[kind],"x23": n,"x25": C - 1,"x26": CL},m2,out_regs=o2)
   assert o2["x27"] == 1
  ptr[kind] += 32 * n; steps += st; C_all += C
 pairs = sorted(((b,g) for (b,kb) in seen[0] for (g,kg) in seen[1] if kb >> K10_LO == kg >> K10_LO))
 assert sorted(found) == pairs and len(set(found)) == len(found) and len(pairs) > 40
 late_b = sum(1 for (b,g) in pairs if b[1] > g[1])
 assert 10 < late_b < len(pairs) - 10
 return {"rounds": len(kinds),"records": len(seen[0]) + len(seen[1]),"equal_key_pairs": len(pairs),"found_by_baby_record": late_b,
 "comparisons": C_all,"steps": steps,"each_pair_reported_once": True}
def audit_packed_match():
 assert (K10_LO,HEAD_LSR,INS_BFI,CMP_AND_POS,HEAD_BITS) == (62,7,(55,9),(55,9),57)
 assert CMP_AND_MASK == (((1 << 9) - 1) << 55) == 0xFF80000000000000 and K10_LO + 9 == 71 == 64 + HEAD_LSR and 128 - 71 == 57
 assert CMP7 == ["ldp","ldr","and","cmp","ccmp","b.ne","cbz"] and CMP_A.count("\n ") == CMP_B.count("\n ") == CMP_EACH + CAP_EACH + 1
 def ex(k): return (k >> K10_LO) & M64
 def f9(k): return (k >> K10_LO) & 0x1FF
 def pk(k,junk=0): return (junk & ((1 << 55) - 1)) | f9(k) << INS_BFI[0]
 def slot(k): return (k >> 71) & ((1 << 57) - 1)
 r = random.Random(20261103)
 assert slot(1 << 71) == 1 and f9(1 << 71) == 0 and ex(1 << 71) == 512 and cmp_field(pk(1 << 71),ex(1 << 71))
 kk = [0,(1 << 192) - 1,1 << 62,((1 << 192) - 1) ^ (1 << 62),0x3FF << 62,1 << 63,1 << 64,(1 << 62) - 1,(1 << 72) - 1,1 << 70,1 << 71,3 << 70]
 nb = [x for x in range(192) if x < 62 or x > 70]
 for k in kk + [r.getrandbits(192) for _ in range(3000)]:
  w0,w1,_ = key_parts(k)
  assert (w1 >> HEAD_LSR) == slot(k) and ex(k) == ((w1 << 64 | w0) >> K10_LO) & M64
  c4 = pk(k,r.getrandbits(64))
  assert cmp_field(c4,ex(k))
  for bit in range(62,71): assert not cmp_field(c4,ex(k ^ 1 << bit))
  assert cmp_field(c4,ex(k ^ (1 << r.choice(nb))))
  assert slot(k ^ (1 << r.randrange(71,128))) != slot(k)
def audit_record_path(rows,baby_pt,giant_pt):
 assert (baby_pt,giant_pt) == (BABY_PIN,GIANT_PIN)
 assert len(GEN12) == REC_GEN and len(RELAY_BODY) == 12 and REC_RELAY == Fraction(len(RELAY_U),UNT) == Fraction(770,64)
 assert IL12 == ["ldp","ldp","ldr","lsr","add","ldp","extr","bfi","stp","str","stp","cbz"] and REC_IL == Fraction(770,64)
 assert IL_BU.count("ldp x0, x1, [x20], #16") == IL_GU.count("ldp x0, x1, [x20], #16") == 2 * UNT - 1
 assert "x5" not in IL_B and "x5" not in IL_G and IL_B.replace("x13, x11, [x12]","x11, x13, [x12]").replace("[x12, #0]","[x12, #8]") == IL_G
 assert len(CMP7) == CMP_EACH and len(UNROLL_SETUP) == 3 and UNROLL_FIX == 10
 assert 62 + 9 + 57 + 64 == 192 and 16 + 16 + 8 == 40 and 64 + 64 + 9 + 63 + 35 <= 256 and 29 + 35 == 64 and 35 <= 55
 assert (1 << 56) * 32 == 1 << 61 and (1 << 57) * 32 == 1 << 62 and (1 << 57) * 16 == 1 << 61
 m = 0xCCCCCCCD
 assert m * 5 == 2 ** 34 + 1 and m < 2 ** 32 and 2 ** 24 * 5 < 2 ** 32 and 4 * (2 ** 34 + 1) + 5 * (2 ** 24 - 1) < 5 * 2 ** 34
 assert all((x * 2341) >> 14 == x // 7 for x in range(1792)) and 2341 < 2 ** 16 and 1791 * 2341 < 2 ** 32
 assert 8 + 3 + 5 == LINK_WAVE
 assert REC_GEN + REC_RELAY + REC_IL == Fraction(2308,64)
 audit_packed_match()
 il = audit_il()
 for bad in ((IL_BU.replace("x1, #7","x1, #6"),IL_GU),(IL_BU,IL_GU.replace("[x12, #8]","[x12, #0]")),(IL_BU.replace("#55, #9","#54, #9"),IL_GU),
  (IL_BU.replace("lsl #55","lsl #54"),IL_GU.replace("lsl #55","lsl #54"))):
  try: audit_il(bad)
  except (AssertionError,KeyError): continue
  raise AssertionError("owner-loop negative control passed")
 audit_smx()
 b = batch_audit.run({}); assert (b['batch_per_lane'],b['round_fixed']) == (str(BATCH),ROUND_FIXED)
 b['owner_loops'] = il
 return b
def proto(cfg):
 from collections import Counter,defaultdict
 rng = random.Random(cfg.get("seed",2026))
 GB,SB,TB,UB = cfg.get("gen_bits",3),cfg.get("s_bits",1),cfg.get("t_bits",2),cfg.get("u_bits",2)
 LANES,WAVES,STEPS = cfg.get("lanes",4),cfg.get("waves",3),cfg.get("steps",6)
 CAP1,CAP2,CAPC,FANIN = cfg.get("cap1",96),cfg.get("cap2",96),cfg.get("capc",8),cfg.get("fanin",4)
 pattern = cfg.get("pattern","uniform")
 assert GB >= SB
 G,R,O = 1 << GB,1 << (TB + SB),1 << (TB + UB)
 NODES = G + R + O
 cost = Counter()
 buf = {}
 tag = {}
 cnt = {}
 ptr = {}
 single = cfg.get("single_count",False)
 emitted,delivered = [],defaultdict(list)
 heads,candidates = defaultdict(list),[]
 overflow = []
 def open_buf(sender,receiver,parity,wave):
  key = (sender,receiver,parity)
  lst = buf.setdefault(key,[])
  w,rd = tag.get(key,(None,None))
  if w != wave:
   assert w is None or rd == w + 1,"buffer reused before read"
   lst.clear(); tag[key] = (wave,None)
  return lst
 def put(sender,receiver,parity,rec,cap,wave):
  lst = open_buf(sender,receiver,parity,wave)
  pw,pv = ptr.get((sender,receiver),(wave,0))
  if pw != wave: assert pv == 0; pw,pv = wave,0
  if pv >= cap * 40:
   overflow.append((sender,receiver,wave)); return False
  lst.append(rec); ptr[(sender,receiver)] = (wave,pv + 40); return True
 def publish(sender,receiver,parity,wave):
  lst = open_buf(sender,receiver,parity,wave)
  pw,pv = ptr.get((sender,receiver),(wave,0))
  val = pv if pw == wave else 0
  if pw != wave: assert pv == 0
  assert val % 40 == 0 and val // 40 == len(lst)
  ptr[(sender,receiver)] = (wave,0)
  cnt[(sender,receiver) if single else (sender,receiver,parity)] = (len(lst),wave)
 def take(sender,receiver,parity,wave):
  key = (sender,receiver,parity)
  w,rd = tag.get(key,(None,None))
  if w is None: return []
  assert w == wave - 1 and rd is None,"read in the wrong wave"
  tag[key] = (w,wave)
  n,cw = cnt[(sender,receiver) if single else key]
  if (n,cw) != (len(buf[key]),w):
   raise RuntimeError("count word mismatch")
  return list(buf[key][:n])
 def relay_of(key,g): return ((key & (R // (1 << SB) - 1)) << SB) | (g & ((1 << SB) - 1))
 def owner_of(key): return key & (O - 1)
 def t_of(x): return x & ((1 << TB) - 1)
 def pack_of(key): return ((key >> K10_LO) & 0x1FF) << INS_BFI[0]
 bk,planted = [[],[]],[0]
 kinds = cfg.get("kinds",[0,1,0,1,1,0])
 def barrier():
  instr,_ = barrier_instr(NODES,FANIN,1)
  cost["barrier"] += instr
 wave_count = Counter()
 seq = 0
 for rnd,phase in enumerate(kinds):
  for wl in range(WAVES + 2):
   w = rnd * (WAVES + 2) + wl
   if wl < WAVES:
    for g in range(G):
     for lane in range(LANES):
      for st in range(STEPS):
       n_out = 2 if rng.random() < 0.98 else 1
       for _ in range(n_out):
        if pattern == "one_owner": key = 0 | (rng.getrandbits(180) << 12)
        elif pattern == "two_owners": key = (rng.randrange(2) * (O // 2)) | (rng.getrandbits(180) << 12)
        elif bk[1 - phase] and rng.random() < 0.25: key = rng.choice(bk[1 - phase]); planted[0] += 1
        else: key = rng.getrandbits(192)
        bk[phase].append(key)
        rec = (key,phase,g,lane,seq); seq += 1
        emitted.append(rec)
        cost["rec_gen"] += 1
        if not put(g,G + relay_of(key,g),w & 1,rec,CAP1,w):
         return {"pattern": pattern,"failed": True,"reason": "link cap overflow (run fails)",
 "overflow_at": overflow[0],"records_before_failure": len(emitted)}
     for rl in range(g & ((1 << SB) - 1),R,1 << SB):
      publish(g,G + rl,w & 1,w)
     cost["link"] += R >> SB
     wave_count["gen"] += 1
   if 1 <= wl <= WAVES:
    for rl in range(R):
     s = rl & ((1 << SB) - 1)
     for g in range(s,G,1 << SB):
      for rec in take(g,G + rl,(w - 1) & 1,w):
       cost["rec_relay"] += 1
       assert relay_of(rec[0],g) == rl
       if not put(G + rl,G + R + owner_of(rec[0]),(w - 1) & 1,rec,CAP2,w - 1):
        return {"pattern": pattern,"failed": True,"reason": "link cap overflow (run fails)",
 "overflow_at": overflow[0],"records_before_failure": len(emitted)}
     for o in range(O):
      if t_of(o) == (rl >> SB): publish(G + rl,G + R + o,(w - 1) & 1,w - 1)
     cost["link"] += (G >> SB) + O // (1 << TB)
     wave_count["relay"] += 1
   if 2 <= wl <= WAVES + 1:
    for o in range(O):
     t = t_of(o)
     for s in range(1 << SB):
      rl = (t << SB) | s
      for rec in take(G + rl,G + R + o,w & 1,w - 1):
       assert owner_of(rec[0]) == o
       delivered[o].append(rec)
       slot = (rec[0] >> (TB + UB)) & 63
       rest = rec[0] >> (TB + UB + 6)
       cost["il"] += 1
       if heads.get((1 - phase,o,slot)): cost["nonempty"] += 1
       for (rb,chb,bb) in heads[(1 - phase,o,slot)]:
        cost["compare"] += 1
        pkm = cmp_field(chb,(rec[0] >> K10_LO) & M64) and rb == rest
        assert pkm == (rb == rest)
        if pkm: candidates.append((o,rec[0]) + ((bb,rec) if phase else (rec,bb)))
       heads[(phase,o,slot)].append((rest,pack_of(rec[0]),rec))
     cost["link"] += 1 << SB
     wave_count["owner"] += 1
   barrier(); cost["barriers"] += 1
 got = [rec for o in delivered for rec in delivered[o]]
 assert Counter(got) == Counter(emitted) and len(got) == len(emitted)
 n_b = sum(1 for r in emitted if r[1] == 0); n_g = len(emitted) - n_b
 assert cost["rec_gen"] == cost["rec_relay"] == cost["il"] == len(emitted)
 byk = defaultdict(list)
 for g in emitted:
  if g[1] == 1: byk[g[0]].append(g)
 pairs = Counter((b,g) for b in emitted if b[1] == 0 for g in byk.get(b[0],()))
 assert Counter((c[2],c[3]) for c in candidates) == pairs
 cand_cost = 0; resolved = 0
 for (o,key,bb,gg) in candidates:
  t = t_of(o); replies = []
  for side in (bb,gg):
   g = side[2]; rl = (t << SB) | (g & ((1 << SB) - 1))
   cand_cost += 2 * REC_RELAY
   assert side in emitted
   rl2 = relay_of(key,g); assert owner_of(key) == o
   cand_cost += 2 * REC_RELAY; replies.append(side[4])
  resolved += 1
 assert resolved == len(candidates) >= planted[0] and (planted[0] > 0 or pattern != "uniform")
 links = G * (R >> SB) + R * (O // (1 << TB))
 expected_link_events = (R >> SB) * wave_count["gen"] + ((G >> SB) + O // (1 << TB)) * wave_count["relay"] + (1 << SB) * wave_count["owner"]
 assert cost["link"] == expected_link_events
 internal,_ = barrier_tree(NODES,FANIN)
 total = (cost["rec_gen"] * REC_GEN + cost["rec_relay"] * REC_RELAY + cost["il"] * REC_IL
 + cost["compare"] * CMP_EACH + cost["link"] * LINK_WAVE + cost["barrier"] + cand_cost)
 formula = (len(emitted) * (REC_GEN + REC_RELAY + REC_IL) + cost["compare"] * CMP_EACH
 + expected_link_events * LINK_WAVE + barrier_instr(NODES,FANIN,cost["barriers"])[0] + cand_cost)
 assert total == formula
 return {"pattern": pattern,"failed": False,"gens": G,"relays": R,"owners": O,"lanes": LANES,"waves_per_phase": WAVES,
 "records": len(emitted),"baby": n_b,"giant": n_g,"delivered_exactly_once": True,
 "every_buffer_read_in_the_wave_after_it_was_written": True,"barriers": cost["barriers"],
 "barrier_tree_internal_nodes": internal,"candidates": len(candidates),"candidates_resolved": resolved,
 "compares": cost["compare"],"links": links,"link_wave_events": cost["link"],
 "nonempty_lookups": cost["nonempty"],"nonempty_fraction": cost["nonempty"] / max(1,cost["il"]),
 "instructions_by_item": {"gen": cost["rec_gen"] * REC_GEN,"relay": str(cost["rec_relay"] * REC_RELAY),
 "owner_insert_lookup": str(cost["il"] * REC_IL),
 "compare": cost["compare"] * CMP_EACH,"link_words": cost["link"] * LINK_WAVE,
 "barriers": cost["barrier"],"candidate_path": str(cand_cost)},
 "instructions_total": str(total),"closed_form_equal": True,"rounds": len(kinds),"each_equal_key_pair_once": True,
 "scope": "sequential simulation of kind-homogeneous rounds; listing applied to counted events"}
def proto_mode(cfg):
 out = []
 c = dict(cfg); c["pattern"] = "uniform"; out.append(proto(c))
 c = dict(cfg); c.update({"pattern": "two_owners","cap2": 8 * cfg.get("cap2",96)}); out.append(proto(c))
 c = dict(cfg); c.update({"pattern": "one_owner"}); out.append(proto(c))
 assert out[2]["failed"] and not out[0]["failed"] and not out[1]["failed"]
 c = dict(cfg); c.update({"pattern": "uniform","single_count": True})
 try:
  proto(c); clobbered = False
 except RuntimeError:
  clobbered = True
 assert clobbered
 out.append({"pattern": "uniform","single_count_word_control": "one count word per link is clobbered by the next wave (detected); per-parity words never are"})
 out.append(lane_schedule(cfg.get("schedule",{})))
 return out
def lane_schedule(cfg):
 """Segment-to-lane assignment: padding lemma."""
 import heapq
 rng = random.Random(cfg.get("seed",7))
 LANES,SEG,TM = cfg.get("lanes",64),cfg.get("seg",32),cfg.get("tm",16)
 res = []
 for kind in ("uniform","extremes","long_first","long_last","alternating"):
  for S in (LANES,LANES + 1,3 * LANES + 5,50 * LANES + 17):
   if kind == "uniform": L = [rng.randint(SEG,SEG + TM) for _ in range(S)]
   elif kind == "extremes": L = [rng.choice((SEG,SEG + TM)) for _ in range(S)]
   elif kind == "long_first": L = sorted(rng.randint(SEG,SEG + TM) for _ in range(S))[::-1]
   elif kind == "long_last": L = sorted(rng.randint(SEG,SEG + TM) for _ in range(S))
   else: L = [SEG + TM if (i // LANES) % 2 == i % 2 else SEG for i in range(S)]
   free = [(0,j) for j in range(LANES)]
   heapq.heapify(free)
   fin = [0] * LANES
   for x in L:
    t,j = heapq.heappop(free)
    fin[j] = t + x
    heapq.heappush(free,(t + x,j))
   span,ops = max(fin),sum(L)
   pad = LANES * span - ops
   assert span - min(fin) <= SEG + TM and pad <= LANES * (SEG + TM) and span <= ops / LANES + SEG + TM
   st = [sum(L[j::LANES]) for j in range(LANES)]
   res.append((kind,S,pad,LANES * max(st) - ops))
 return {"lanes": LANES,"seg": SEG,"treemax": TM,"cases": len(res),"padding_lemma_holds": True,
 "padding_bound": LANES * (SEG + TM),
 "max_greedy_padding": max(r[2] for r in res),"max_stride_padding": max(r[3] for r in res)}
KERNEL = 894
ROWS_PT, VROWS = 283, batch_audit.VROWS
BABY_PIN,GIANT_PIN = Fraction(5936963848153,6553600000),Fraction(7235395848153,6553600000)
def exp_neg_bounds(x,terms=60):
 """Rational enclosure of exp(-x), 0<x<1: alternating Taylor, decreasing terms."""
 assert Fraction(0) < x < Fraction(1)
 t,s,lo,up = Fraction(1),Fraction(0),None,None
 for n in range(terms):
  s += t
  if n % 2 == 0: up = s
  else: lo = s
  t = t * (-x) / (n + 1)
 assert lo < up and up - lo < Fraction(1,10 ** 50)
 return lo,up
def exp4(K):
 """Enclosure of exp(-K), 1 < K < 8, as exp(-K/8)^8."""
 lo,up = exp_neg_bounds(K / 8); assert 0 < lo and 1 < K < 8 and up ** 8 - lo ** 8 < Fraction(1,10 ** 48)
 return lo ** 8,up ** 8
def gauss_int_up(K):
 """Upper enclosure of J(K) = int_0^1 exp(-K y^2) dy = sum (-K)^n/(n!(2n+1)); stops on a decreasing negative term."""
 s,t,n = Fraction(0),Fraction(1),0
 while True:
  s += t / (2 * n + 1); nt = t * (-K) / (n + 1)
  if n > K and nt < 0 and -nt / (2 * n + 3) < Fraction(1,10 ** 50): return s
  t,n = nt,n + 1
def success_floor_lower(K):
 """Lower bound of an opened run's success (s2)."""
 lo,up = exp4(K)
 x = 1 - up - DELTA_P - FAILP
 return Fraction(x.numerator * 10 ** 60 // x.denominator,10 ** 60)
TARGET,GATE_C,KT_RUN,DELTA_P = Fraction(39,100),64,Fraction(10463,2000),Fraction(1,2000)
FAILP = Fraction(1,1 << 20) + Fraction(7,1 << 40) + Fraction(1,1 << 108)
def gate(lo):
 Kg = -(-(TARGET * (1 << 64)).numerator * lo.denominator // ((TARGET * (1 << 64)).denominator * lo.numerator))
 assert Fraction(Kg,1 << 64) * lo >= TARGET > Fraction(Kg - 1,1 << 64) * lo and 0 < Kg < 1 << 64
 return Kg
def ledger(cfg):
 from decimal import Decimal,localcontext,ROUND_CEILING
 ceil = lambda x: -(-x.numerator // x.denominator) if isinstance(x,Fraction) else x
 r,d,m = R_PRODUCTION,D_PRODUCTION,M_PRODUCTION
 assert (r - 1) % d == 0 and m % 6 == 0
 v2,t_ = 0,r - 1
 while t_ % 2 == 0: t_ //= 2; v2 += 1
 assert v2 == 32
 lam = (-0xd201000000010000) ** 2 - 1
 assert (lam * lam + lam + 1) % r == 0 and pow(lam,m,r) == 1 and pow(2,m,r) != 1
 KT,F = Fraction(cfg.get('KT',KT_RUN)),Fraction(998,1000)
 ELL,Q = Fraction(193,200),Fraction(cfg.get('ratio','609/500'))
 KX,BETA,ROUTE = kernel_check.leaf_worst(),batch_audit.BETA_LEN,6
 assert (KX,BETA,len(batch_audit.BUILD)) == (233,128,24) and max(n for n,_ in kernel_check.lane_paths(kernel_check.kernel_text())) == KERNEL
 NB_ = lambda g: -(-(Q.numerator * g) // Q.denominator)
 lo_,hi_ = 1,1 << 120
 while lo_ < hi_:
  mid_ = (lo_ + hi_) // 2
  if 6 * NB_(mid_) * mid_ >= KT * m: hi_ = mid_
  else: lo_ = mid_ + 1
 NG = lo_; NB = NB_(NG)
 K = Fraction(6 * NB * NG,m)
 assert 6 * NB_(NG - 1) * (NG - 1) < KT * m <= 6 * NB * NG and K >= KT
 X2,GMAX,EMIN,EMAX = 1 << 76,1 << 47,(1 << 56) - (1 << 52),1 << 57
 LANES,SEG,TREEMAX,WSTEP,DEPTH = 1 << 19,1 << 22,1 << 20,1 << 33,1 << 17
 NBT,NGT = NB,NG
 NB,NG = NB + EMAX,NG + EMAX
 GENS,SBITS,TBITS,UBITS = 1 << 47,20,30,32
 RELAYS,OWNERS = 1 << (TBITS + SBITS),1 << (TBITS + UBITS)
 HEADS,POOLB,POOLG,CAP1,CAP2,CAPC,CAND = 1 << HEAD_BITS,1 << 56,1 << 57,1 << 24,1 << 21,1 << 8,1 << 45
 assert TBITS + UBITS == K10_LO and HEADS * OWNERS == 1 << 119
 RECB,FANIN = 40,1 << 10
 SM = 8192000
 SM77 = -(-77 * SM // 255)
 assert SM77 * 255 >= 77 * SM and SM77 < SM and 4 * 77 > 255
 inv_share = Fraction(ROUND_FIXED + 5,LANES)
 assert inv_share < 1 and (ROUND_FIXED,BATCH) == (172597,Fraction(9553,8))
 tbl_share = ceil(Fraction(232864,LANES))
 PAIR_BK = 8 + 1 + 6 + 1 + 2 + 4 + inv_share
 ROOT_BK = 5 + 6 + 2 + 1 + 2 + inv_share + 5
 TABLE_BK = 8 + 7 + 6 + tbl_share + 5 + 21
 assert TABLE_BK <= 50
 assert 2 * KERNEL + BATCH + PAIR_BK + ROUTE < 3300
 PAIR = 2 * KERNEL + BATCH + PAIR_BK + ROUTE
 ROOT = KERNEL + BATCH + ROOT_BK
 CAPF = Fraction(7,32)
 KXC = KX
 PAIR_PT = Fraction(BATCH + PAIR_BK + ROUTE,2) + ELL * KXC + (1 - ELL) * KERNEL
 BABY_PT = F * PAIR_PT + (1 - F) * ROOT + REC_GEN + REC_RELAY + REC_IL
 GIANT_PT = F * PAIR_PT + (1 - F) * ROOT + 3 * (REC_GEN + REC_RELAY + REC_IL) - 2 + BETA
 LOOKS = 3 * NG
 BASES_B = ceil(Fraction(NB,EMIN)) + 1; BASES = BASES_B + ceil(Fraction(NG,EMIN)) + 1
 BSIDE = 1 << 42
 assert BASES < 1 << 63 and BASES_B < d * BSIDE and BASES - BASES_B < d * BSIDE
 assert 6 * BASES_B * (BASES - BASES_B) * (1 << 108) < m
 mu_f = Fraction(3 * NB * NG,1 << 192)
 assert mu_f < CAND - (1 << 20) and mu_f / (CAND - (1 << 20) - mu_f) ** 2 < Fraction(1,1 << 40)
 GPC = GENS // d
 BPG = ceil(Fraction(BSIDE,GPC))
 SEGLEN = SEG + TREEMAX
 LANE_OPS = (EMAX // SEG // LANES + 1) * SEGLEN
 assert 2 * WSTEP * LANES == (1 << 23) * (1 << TBITS)
 PLAN_ROUNDS = ceil(Fraction(EMAX * 16,1 << 48))
 RMAX,RW = 2 * BPG,ceil(Fraction(LANE_OPS,WSTEP // 4)) + 10
 WAVES = RMAX * RW + PLAN_ROUNDS + 64
 PADOPS = BASES * LANES * SEGLEN + 2 * LANES * LANE_OPS
 fixed_call_bound = (BASES + 2) * (8 * 48 * LANE_OPS + 64)
 assert fixed_call_bound % 2 == 0 and 32 + 1 + 5 + 2 == 40
 rows = {}
 rows["records: N_b + N_g points"] = ceil(NB * BABY_PT + NG * GIANT_PT)
 rows["padding ops at the pair price"] = ceil(PADOPS * (PAIR + 2 * BETA))
 rows["side-output fixed calls"] = fixed_call_bound
 rows["listing allowance (1/point)"] = NB + NG
 rows["owner pool-cap wave checks"] = 2 * OWNERS * (5 * WAVES + 6)
 for POOL,cnt in ((POOLB,NB),(POOLG,LOOKS)):
  THRESH = POOL - (1 << SBITS) * CAP2
  assert (1 << SBITS) * CAP2 == (1 << 41) and THRESH + (1 << 41) == POOL
  dOwn = Fraction(THRESH * OWNERS - cnt,cnt); assert 0 < dOwn < 1 and Fraction(cnt,OWNERS) * dOwn * dOwn / 3 > 82
 PBB,PBG = HEADS * 16,HEADS * 16 + POOLB * 32
 TPB,TPG = PBB + 32 * (POOLB - (1 << 41)),PBG + 32 * (POOLG - (1 << 41))
 assert PBB == 1 << 61 and PBG == 1 << 62 and PBG + 32 * POOLG == 1 << 63 and TPB + 32 * (1 << 41) == PBG
 assert all(sum(1 for sh in (0,16,32,48) if (v >> sh) & 0xFFFF) <= 2 for v in (PBB + 32,PBG + 32,TPB,TPG))
 assert SMX == 630722 and SM77 == 2473663 and SMX < SM77
 rows["segment restarts"] = (ceil(Fraction(NB + NG,SEG)) + BASES * LANES) * (SMX + ROWS_PT + 16)
 TL,Mt = [],3
 while Mt < GMAX:
  pr = min(Mt - 2,-(-(GMAX - Mt - 1) // 2)); TL.append(pr); Mt = Mt + 2 * pr + 1
 assert len(TL) == 30 and Mt >= GMAX and all(1 <= x for x in TL)
 TROUNDS = sum(-(-x // LANES) + 1 for x in TL)
 rows["gap tables"] = BASES * (ceil(TROUNDS * (LANES * PAIR + 128)) + GMAX * (21 + VROWS))
 rows["table strand starts and base points"] = BASES * (LANES + 1) * (SM + ROWS_PT)
 # s45: per-round cumulative comparison budgets, cap_j = ceil(CSC * nb_j * ng_j) + SLACK per owner
 dS,SLACK = Fraction(1,10 ** 5),1 << 40
 CSC = (1 + dS) * Fraction(3,OWNERS * OWNERS * HEADS)
 cmp_cap = ceil(CSC * NB * NG) + SLACK
 for e in range(0,41,4):
  for nb_,ng_ in ((max(NB >> e,1 << 100),NG),(NB,max(NG >> e,1 << 100)),(max(NB >> e,1 << 100),max(NG >> e,1 << 100))):
   Bm,Gm = Fraction(nb_,OWNERS),Fraction(3 * ng_,OWNERS)
   dB,dG = Bm / 10 ** 7 + (1 << 30),Gm / 10 ** 7 + (1 << 30)
   assert all(min(dx * dx / (3 * mx),dx / 3) > 120 for dx,mx in ((dB,Bm),(dG,Gm)))
   assert 2 * Fraction(1 << 39) ** 2 / ((Gm + dG + 1800) * 64 * 64) > 120
   assert (Bm + dB) * (Gm + dG) / HEADS + (1 << 39) <= ceil(CSC * nb_ * ng_) + SLACK
 rows["comparisons at the per-owner cap"] = OWNERS * cmp_cap * CMP_EACH
 rows["comparison cap enforcement"] = OWNERS * cmp_cap * CAP_EACH
 rows["per-round comparison budgets"] = OWNERS * RMAX * 8
 rows["owner initialisation"] = OWNERS * (HEADS * 16 // 64 * 2 + 64)
 LINKS1,LINKS2 = GENS * (1 << TBITS),RELAYS * (1 << UBITS)
 COPY = d * ((1 << 6) + (1 << 16) + (1 << 26))
 NODES = OWNERS + RELAYS + GENS + 2 * d + COPY + (1 << 14)
 bar,BNODES = barrier_instr(NODES,FANIN,WAVES)
 rows["per-wave link words, memory barriers"] = WAVES * ((LINKS1 + LINKS2) * LINK_WAVE + NODES * 24)
 rows["global barriers"] = bar
 rows["plan distribution (copy tree)"] = d * EMAX * ((1 << 6) + (1 << 16) + (1 << 26) + (1 << 27)) * 4
 rows["plan distribution (copy tree)"] += GENS * (EMAX // SEG + 2) * 4
 WIN = 2 * (2 * GMAX // d)
 KAPPA = GMAX // (16 * d * d); assert KAPPA == 8
 rows["planner (scans, matching, emission)"] = d * EMAX * (6 * WIN * 16 + 5 * 64 * KAPPA * (2 * DEPTH * 3 + 20) + 64 + 16)
 rows["planner descendant traversals (adjacency lists)"] = d * EMAX * (64 * TREEMAX + 8 * DEPTH)
 AUG_IN = d * EMAX
 assert (CC,CCY) == (32,8) and OPS_BOUND == 2 * CC + CC * (6 + CC + CCY * (5 + CC)) == 10752
 rows["planner augmentation: partner cache"] = 2 * AUG_IN * WIN * 16
 rows["planner augmentation: tentative operations"] = 2 * AUG_IN * OPS_BOUND * (16 * DEPTH + 64 + 16)
 rows["planner subtree-size maintenance and pair-fraction check"] = (d * EMAX // 2) * (4 * DEPTH + 4) + d * EMAX * 4
 rows["partial-plan remainders at the root price"] = ceil(2 * EMAX * (ROOT - PAIR_PT))
 rows["character sieve, class lists, their transfer"] = X2 * (400000 + 16 + 20) + d * (1 << 33) * 400000
 rows["base interval sieve, characters, assignment"] = ((1 << 70) + (1 << 43)) * 16 + (1 << 65) * (400000 + 64) + BASES * 64
 rows["candidate delivery and replies"] = ceil(CAND * (8 * REC_RELAY + 200))
 rows["candidate re-walks"] = CAND * 2 * (SEGLEN * (3300 + 232864) + (GMAX - 1) * (2000 + 50 + VROWS + ROWS_PT) + 2 * (SM + ROWS_PT))
 rows["candidate tests"] = CAND * (96 * SM + (1 << 16) * 20)
 rows["window-table initialisation"] = (LINKS1 + LINKS2) * 8
 rows["chain-id formation"] = (ceil(Fraction(NB + NG,SEG)) + BASES * LANES) * 4
 rows["candidate-link wave counts"] = 4 * (LINKS1 + LINKS2) * 9 * RMAX
 assert RMAX == RMAXC and audit_round_ctl() == 6 and 8 + 10 <= RC_FIX
 rows["round control and stop flags"] = (NODES + BNODES) * RC_FIX * RMAX
 rows["candidate formation"] = CAND * 16
 rows["transport unroll setup and tails"] = UNROLL_FIX * ((LINKS1 + LINKS2) * WAVES + 2 * OWNERS * WAVES)
 rows["beta tables: build (s38)"] = GENS * (batch_audit.TABLE_ENTRIES * len(batch_audit.BUILD) + 7 * 2000 + 64)
 rows["row-product tables (s46)"] = GENS * (1 << 15)
 rows["reduction tables (s47)"] = GENS * (8 * (1 << 48) * 24 + 512 * 16)
 rows["root recovery and input distribution"] = (1 << 24) * 2000 + NODES * 64 + (1 << 40)
 batch = audit_record_path(rows,BABY_PT,GIANT_PT)
 total = sum(rows.values())
 lo_run = success_floor_lower(K); Kg = gate(lo_run)
 MAXR = GENS * EMAX
 JUP = gauss_int_up(KT)
 EB = ceil(NBT * (JUP + DELTA_P + FAILP)) + MAXR + 2 * EMAX
 EG = ceil(NGT * (JUP + DELTA_P + FAILP)) + MAXR + 2 * EMAX
 assert MAXR * 5000 < min(NBT,NGT) and WAVES < 1 << 24 and EB < NBT and EG < NGT
 AW = rows["listing allowance (1/point)"]; AL_E = (EB + EG) + (1 << 116); assert AL_E < AW
 TB = rows["gap tables"] + rows["table strand starts and base points"]; BE = ceil(Fraction(EB + EG,EMIN)) + 2; assert TB % BASES == 0 and BE < BASES
 SG = ceil(Fraction(EB + EG,SEG)) + BE * LANES
 PE = ceil((BE * LANES * SEGLEN + 2 * LANES * LANE_OPS) * (PAIR + 2 * BETA)) + (BE + 2) * (8 * 48 * LANE_OPS + 64) + SG * (SMX + ROWS_PT + 20)
 PW = rows["padding ops at the pair price"] + rows["side-output fixed calls"] + rows["segment restarts"] + rows["chain-id formation"]; assert PE < PW
 DU = Fraction(2 * MAXR,min(NBT,NGT))
 UEND = Fraction(NB * NG,NBT * NGT)
 CW = Fraction(1) / KT + DU + (DELTA_P + FAILP) * UEND
 CB = OWNERS * (CMP_EACH + CAP_EACH) * (ceil(CSC * NBT * NGT * CW) + SLACK + RMAX + 1)
 CR = rows["comparisons at the per-owner cap"] + rows["comparison cap enforcement"]; assert CB < CR
 pre = sum(rows[k] for k in ("owner initialisation","plan distribution (copy tree)","planner (scans, matching, emission)","planner descendant traversals (adjacency lists)",
 "planner augmentation: partner cache","planner augmentation: tentative operations","planner subtree-size maintenance and pair-fraction check",
 "character sieve, class lists, their transfer","base interval sieve, characters, assignment","window-table initialisation","beta tables: build (s38)","row-product tables (s46)","reduction tables (s47)"))
 e_run = (total - rows["records: N_b + N_g points"] - AW + ceil(EB * BABY_PT + EG * GIANT_PT) + AL_E - TB + TB // BASES * BE - PW + PE - CR + CB)
 E_num = GATE_C * (1 << 64) + Kg * e_run
 W = Fraction(E_num,2000 << 64)
 W_run = Fraction(total + GATE_C,2000)
 assert 2 * (PLAN_ROUNDS + 64) < WAVES
 gen_mem = EMAX * 16 + (1 << 48) * 1024 + GMAX * 912 + batch_audit.TABLE_BYTES + 512 + LANES * (2 * DEPTH + 3) * 672 + 2 * 512 * 64 + LANES * 16 + (1 << TBITS) * (16 + 2 * CAP1 * RECB) + (1 << TBITS) * CAPC * 64 + 2 * BPG * 16 + (EMAX // SEG + 2) * 8 + (LANES + 1) * (batch_audit.SB + 400) + (1 << 20) + (1 << 30)
 relay_mem = (GENS >> SBITS) * (16 + 2 * CAP1 * RECB) + (1 << UBITS) * (16 + 2 * CAP2 * RECB) + (1 << UBITS) * CAPC * 64 + (GENS >> SBITS) * CAPC * 64 + (1 << 30)
 owner_mem = HEADS * 16 + POOLB * 32 + POOLG * 32 + (1 << SBITS) * (16 + 2 * CAP2 * RECB) + (1 << 16) * 64 + (1 << SBITS) * CAPC * 64 + (1 << 30)
 plan_mem = EMAX * 43 + (1 << 40) + (1 << 36) + (1 << 58) * 8 + (1 << 6) * (1 << 48) * 2 + (1 << 30)
 sieve_mem = (1 << 56) * 4 + (1 << 33) * 4 + d * (1 << 40) + (1 << 30)
 copy_mem = 2 * (1 << 48) + (1 << 10) * (1 << 48) + (1 << 30)
 assert gen_mem < (1 << 64)
 tx = '\n'.join(set('\n'.join(batch_audit.LISTING + batch_audit.beta_listing() + batch_audit.BUILD + smx_audit.LISTING + GEN12 + RELAY_BODY + [IL_BU,IL_GU] + [t for q in (0,1) for t in (kernel_check.kernel_text(q),kernel_check.leaf_text(kernel_check.kernel_text(),q))]).split('\n')))
 assert not re.search(r'\bb\.(lt|ge|gt|le|mi|pl|vs|vc)\b|,\s*(lt|ge|gt|le|mi|pl|vs|vc)\b|\b(ldrs|sxt|asr)',tx) and set(re.findall(r'tbn?z\s+(x\d+),\s*#63,\s*(\.\w+)',tx)) == {('x16','.Lx_nonneg'),('x19','.Lfound'),('x22','.Lskip')}
 for v in (relay_mem,plan_mem,sieve_mem,copy_mem):
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
  bound = lambda x,scale: int((log2(x) * scale).to_integral_value(rounding=ROUND_CEILING)) / scale
  succ_exact = Fraction(Kg,1 << 64) * lo_run
  assert succ_exact >= TARGET and Kg * e_run < E_num and e_run < total
  return {"K": str(K),"N_b": NBT,"N_g": NGT,"Q": str(Q),"expected_baby_points": EB,"expected_giant_points": EG,"gauss_integral_upper": float(JUP),
 "expected_run_instructions": e_run,"expected_allowance": AL_E,"expected_bases": BE,"expected_segments": SG,"walk_rows_expected": PE,
 "comparison_rows_expected": CB,"comparison_weight": float(CW),"cmp_cap_full": cmp_cap,"rounds_max": RMAX,"round_waves": RW,"round_mass_max": MAXR,
 "expected_baby_share": float(Fraction(EB,NBT)),"expected_giant_share": float(Fraction(EG,NGT)),"bases": BASES,"baby_bases": BASES_B,
 "prime_margin": float(Fraction(1 << 70,53) / BASES),"eps_lit_bound": float(eps_lit()),"baby_point": str(BABY_PT),"giant_point": str(GIANT_PT),
 "ell": str(ELL),"f": str(F),"kernel": KERNEL,"batch_round": batch,"pair_op_instructions_two_points": str(PAIR),"root_op_instructions": str(ROOT),
 "lane_ops": LANE_OPS,"waves": WAVES,"nodes_log2": math.log2(NODES),"barrier_internal_nodes": BNODES,
 "rows_log2_units": {k: float(log2(Fraction(v,2000))) for k,v in rows.items()},"rows_instructions": rows,
 "run_total_instructions": total,"preprocessing_instructions": pre,"preprocessing_log2_ceil_1e4": bound(Fraction(pre,2000),10 ** 4),"expected_numerator_x2pow64": E_num,"expected_log2": str(s),"score_numerator": score,"score_denominator": 10 ** 15,
 "score": score / 10 ** 15,"false_match_mean_log2": float(log2(mu_f)),"h5_conditional_mean_nu": float(K / (1 - exp4(K)[1] - Fraction(5,10000))),
 "bases_log2": float(log2(Fraction(BASES))),"memory_log2_ceil_1e4": bound(Fraction(mem),10 ** 4),
 "per_node_bytes_log2": {"gen": math.log2(gen_mem),"relay": math.log2(relay_mem),"owner": math.log2(owner_mem),"planner": math.log2(plan_mem),"sieve": math.log2(sieve_mem),"copy": math.log2(copy_mem)},
 "success_floor": float(succ_exact),"success_floor_exact_lower": str(succ_exact),"success_claimed": 0.3900,"KT": str(KT),"gate_K": Kg,"gate_cost": GATE_C,
 "E_numerator_x2pow64": E_num,"run_success_lower": str(lo_run),"run_success_float": float(lo_run),"run_instructions": total + GATE_C,"run_log2_units": str(log2(W_run)),
 "run_score_numerator": int((log2(W_run) * 10 ** 15).to_integral_value(rounding=ROUND_CEILING)),"all_per_node_structures_below_2pow64": True,"generating_node_bytes": gen_mem}
def main():
 cfg = json.loads(sys.stdin.read())
 mode = cfg.get("mode","all")
 t0 = time.time()
 out = {"mode": mode}
 if mode in ("ledger","all","declared"):
  out["ledger"] = ledger(cfg)
  import a64,gc; a64._a64_parse.cache_clear(); gc.collect()
 if mode in ("planner","all"):
  out["planner"] = planner_mode(cfg.get("planner",{
 "specs": [{"d": 16,"logX2": 15,"seed": 0,"classes": 8}]}))
 if mode in ("ec","all","declared"):
  out["ec"] = ec_attack(cfg.get("ec",{
 "bits": 23,"logd": 5,"K": "5.2315","trials": 60,"seed": 2026,"q": "59/50",
 "keybits": 18,"nodebits": 5,"headbits": 6,"head_cap": 64,
 "cand_cap": 4096,"logX2": 12,"logG": 11,"spotchecks": 6}))
 if mode in ("proto","all","declared"):
  out["proto"] = proto_mode(cfg.get("proto",{}))
 if mode in ("census","all","declared"):
  out["census"] = [census(dict(cfg.get("census",
 {"m_bits": 23,"k": 8,"K": "5.2315","seed": 0,"q": "59/50","nbar": 6})))]
 if mode == "all":
  out["h2"] = h2_sizes()
 out["seconds"] = round(time.time() - t0,2)
 print(json.dumps(out,indent=1))
if __name__ == "__main__":
 main()
