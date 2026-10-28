"""Six-class parallel Pollard rho on secp256k1 (proof.md).  Part A: production identities and the canonical
index.  Part B: the walk of proof.md 5.1 on a toy curve.  Kernels mode: proof.md sections 2-5, the shipped
kernels run by an AArch64 interpreter, the counts and the ledger.
"""
import json, math, random, re, sys, time
from fractions import Fraction

KAPPA = 1.0001  # slowdown allowance of H1 (proof.md section 5.2)
P0 = 0.3900653  # probability floor of H1 at S0 = S_ops / KAPPA
PERIOD = 4096  # anchor period of the walk (proof.md section 5.1)
STRIDE = 128  # checkpoint compares only in rounds (length + phase) % STRIDE == 0
TRAV_CAP = 32  # traversal cap: the escape handles cycles of length <= 32

T0 = time.time()
cfg = json.load(sys.stdin)
curve = cfg.get("curve", "q25")
trials = cfg.get("trials", 100)
R = cfg.get("R", 256)
dp_bits = cfg.get("dp_bits", 4)
seed = cfg.get("seed", 1)
prod_checks = cfg.get("production_checks", 2)
budget = cfg.get("seconds", 14)
dump = cfg.get("dump", False)
pipeline = cfg.get("pipeline", False)
long_trails = cfg.get("long_trails", 0)
long_len = cfg.get("long_len", 4096)
canon_random = cfg.get("canon_random", None)
canon_stress = cfg.get("canon_stress", None)
run_toy = cfg.get("toy", True)
CLS = cfg.get("classes", 6)
kernels = cfg.get("kernels", False)
rng = random.Random(seed)

# -- curves
PROD = dict(
 q=int("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16),
 r=int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16),
 h=1,
 gx=int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16),
 gy=int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16),
 beta=int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16),
 lam=int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16),
)
CB = 7  # curve coefficient b (y^2 = x^3 + 7) for production and toys
TOYS = {
 "q25": (33556867, 33563419, 1),
 "q29": (536873251, 536839203, 3),
}


class Curve:
 """Affine arithmetic on y^2 = x^3 + 7 over F_q (points are tuples, None = infinity)."""

 def __init__(self, q, r, h):
  self.q, self.r, self.h = q, r, h

 def on_curve(self, P):
  return P is None or (P[1] * P[1] - P[0] ** 3 - CB) % self.q == 0

 def add(self, P, Q):
  q = self.q
  if P is None:
   return Q
  if Q is None:
   return P
  x1, y1 = P
  x2, y2 = Q
  if x1 == x2:
   if (y1 + y2) % q == 0:
    return None
   l = 3 * x1 * x1 * pow(2 * y1, -1, q) % q
  else:
   l = (y2 - y1) * pow(x2 - x1, -1, q) % q
  x3 = (l * l - x1 - x2) % q
  return (x3, (l * (x1 - x3) - y1) % q)

 def neg(self, P):
  return None if P is None else (P[0], (-P[1]) % self.q)

 def mul(self, k, P):
  k %= self.r
  return self.mul_raw(k, P)

 def mul_raw(self, k, P):
  """Scalar multiplication without reducing by a presumed point order."""
  Rr = None
  while k:
   if k & 1:
    Rr = self.add(Rr, P)
   P = self.add(P, P)
   k >>= 1
  return Rr

 def fixed_base(self, P):
  """Precompute j * 2^(8i) * P for a fixed base P (toy curves only)."""
  self.fb = []
  Q = P
  for i in range((self.r.bit_length() + 7) // 8):
   row = [None]
   for j in range(1, 256):
    row.append(self.add(row[-1], Q))
   self.fb.append(row)
   Q = self.add(row[-1], Q)

 def fmul(self, k):
  """[k]G via the fixed-base table: at most 4 additions on the toy curves."""
  k %= self.r
  Rr, i = None, 0
  while k:
   Rr = self.add(Rr, self.fb[i][k & 255])
   k >>= 8
   i += 1
  return Rr

 def phi(self, P, j=1):
  if P is None:
   return None
  return (P[0] * pow(self.beta, j % 3, self.q) % self.q, P[1])


def sqrt_mod(a, q):
 a %= q
 if a == 0:
  return 0
 if pow(a, (q - 1) // 2, q) != 1:
  return None
 s, t = 0, q - 1
 while t % 2 == 0:
  t //= 2
  s += 1
 z = 2
 while pow(z, (q - 1) // 2, q) != q - 1:
  z += 1
 m, c, tt, x = s, pow(z, t, q), pow(a, t, q), pow(a, (t + 1) // 2, q)
 while tt != 1:
  i, t2 = 1, tt * tt % q
  while t2 != 1:
   t2 = t2 * t2 % q
   i += 1
  b = pow(c, 1 << (m - i - 1), q)
  m, c, tt, x = i, b * b % q, tt * b * b % q, x * b % q
 return x


def cube_roots_of_unity(m):
 for g in range(2, 1000):
  c = pow(g, (m - 1) // 3, m)
  if c != 1:
   return [c, c * c % m]


def make_toy(name):
 q, n, h = TOYS[name]
 r = n // h
 E = Curve(q, r, h)
 while True:
  x = rng.randrange(q)
  y = sqrt_mod(x ** 3 + CB, q)
  if y is None:
   continue
  G = None  # cofactor multiple [h](x, y), computed without reducing h mod r
  P, k = (x, y), h
  while k:
   if k & 1:
    G = E.add(G, P)
   P = E.add(P, P)
   k >>= 1
  if G is not None and E.mul_raw(r, G) is None:
   break
 E.G = G
 E.fixed_base(G)
 for b in cube_roots_of_unity(q):
  for l in cube_roots_of_unity(r):
   if E.mul(l, G) == (G[0] * b % q, G[1]):
    E.beta, E.lam = b, l
 return E


# -- the class walk
class Walker:
 """The walk of proof.md section 5.1 (upstream-proof.md section 2 with checkpoint-only cycle detection), on curve E."""

 def __init__(self, E, R, dp_bits, T, rng, start_point=None):
  self.E, self.R, self.dp_bits = E, R, dp_bits
  self.rbits = R.bit_length() - 1
  assert 1 << self.rbits == R
  self.dpmask = (1 << dp_bits) - 1
  self.Lmax = 20 << dp_bits
  r = E.r
  lam = E.lam
  self.units = [(-1) ** s * pow(lam, j, r) % r for s in (0, 1) for j in (0, 1, 2)]
  self.uinv = [pow(u, -1, r) for u in self.units]
  self.T = T
  self.table = []  # table[h][t] = (point u_t^{-1} M_h, coeff c, coeff d)
  while len(self.table) < R:
   c, d = rng.randrange(r), rng.randrange(r)
   M = start_point(c, d) if start_point else E.add(E.mul(c, E.G), E.mul(d, T))
   if M is None:  # M_h = O: redraw (probability 1/r)
    continue
   row = []
   for t in range(6):
    s, j = divmod(t, 3)
    Mt = E.phi(M, -j)
    if s:
     Mt = E.neg(Mt)
    ui = self.uinv[t]
    row.append((Mt, ui * c % r, ui * d % r))
   self.table.append(row)
  self.stats = dict(steps=0, doublings=0, esccp=0, abandoned=0, exceptional=0)

 def escape(self, P, a, b, track):
  """Escape from the class of P by doubling (the attack)."""
  r = self.E.r
  self.stats["doublings"] += 1
  Q = self.E.add(P, P)  # doubling of a G1 point never hits y = 0 (odd order)
  if track:
   a, b = 2 * a % r, 2 * b % r
  return Q, a, b

 def add_checked(self, P, Q):
  """Incomplete affine addition: None signals an exceptional operand pair
        (equal x-coordinates, i.e. P = +-Q), on which the trail is abandoned."""
  if P[0] == Q[0]:
   return None
  return self.E.add(P, Q)

 def canon(self, P):
  """Return (K, t, canonical point): canonical = units[t] * P."""
  q, beta = self.E.q, self.E.beta
  x, y = P
  if CLS == 1:  # control: point itself, key injective on points
   return x + q * y, 0, P
  if CLS == 2:  # control: negation classes, key x
   if q - y < y:
    return x, 3, (x, q - y)
   return x, 0, P
  bx = beta * x % q
  b2x = (-x - bx) % q  # beta^2 x = -(1+beta) x
  j = 0
  cx = x
  if bx < cx:
   cx, j = bx, 1
  if b2x < cx:
   cx, j = b2x, 2
  if y & 1:  # the class key is the even one of y and q - y
   return q - y, 3 + j, (cx, q - y)
  return y, j, (cx, y)

 def step_target(self, P, K, t):
  h = (K >> 1) & (self.R - 1)  # bit 0 of K is always 0
  return self.table[h][t]

 def is_dp(self, K):
  return ((K >> (self.rbits + 1)) & self.dpmask) == 0

 def run_trail(self, P, coeffs=None, on_class=None, phase=0, rule=(STRIDE, PERIOD)):
  """Walk from P to a distinguished class or the trail cap (proof.md
    section 5.1).  Returns (status, K, length, P_end, coeffs_end); coeffs=(a, b) tracks P = aG + bT for
    re-walks; on_class(K) is called on every class arrival; rule = (stride, anchor period)."""
  E, r = self.E, self.E.r
  stride, period = rule
  cp = None  # checkpoint key
  length = 0
  a = b = None
  if coeffs is not None:
   a, b = coeffs
  while True:
   K, t, c = self.canon(P)
   if on_class is not None:
    on_class(K)
   if self.is_dp(K):
    if coeffs is not None:
     u = self.units[t]
     return "dp", K, length, c, (u * a % r, u * b % r)
    return "dp", K, length, c, None
   if length >= self.Lmax:
    self.stats["abandoned"] += 1
    return "cap", K, length, c, None
   if cp is not None and (length + phase) % stride == 0 and K == cp:
    self.stats["esccp"] += 1
    best = (K, P, a, b)
    Q, qa, qb = P, a, b
    trav = 0
    while True:
     trav += 1
     if trav > TRAV_CAP:  # not a cycle of length <= TRAV_CAP: abandon the trail
      self.stats["abandoned"] += 1
      return "trav", K, length, None, None
     TP, tc, td = self.step_target(Q, *self.canon(Q)[:2])
     Q = self.add_checked(Q, TP)
     if Q is None:
      self.stats["exceptional"] += 1
      return "exc", K, length, None, None
     if coeffs is not None:
      qa, qb = (qa + tc) % r, (qb + td) % r
     length += 1
     self.stats["steps"] += 1
     KQ = self.canon(Q)[0]
     if on_class is not None:
      on_class(KQ)  # traversal arrivals are counted visits too
     if KQ == K:
      break
     if KQ < best[0]:
      best = (KQ, Q, qa, qb)
    _, P, a, b = best
    P, a, b = self.escape(P, a, b, coeffs is not None)
    if P is None:
     self.stats["exceptional"] += 1
     return "exc", K, length, None, None
    cp = None
    length += 1
    continue
   if (length + phase) % period == 0:  # anchor rounds: the trail's phase is its start round mod period
    cp = K
   TP, tc, td = self.step_target(P, K, t)
   P = self.add_checked(P, TP)
   if P is None:
    self.stats["exceptional"] += 1
    return "exc", K, length, None, None
   if coeffs is not None:
    a, b = (a + tc) % r, (b + td) % r
   self.stats["steps"] += 1
   length += 1


# -- fast canonical index (production)
M64 = (1 << 64) - 1
GQ = (PROD["beta"] << 384) // PROD["q"]  # floor(beta/p * 2^384): 6 limbs
GL = [(GQ >> (64 * j)) & M64 for j in range(6)]


def exact_j(x):
 q, b = PROD["q"], PROD["beta"]
 v = [x, b * x % q, b * b % q * x % q]
 return min(range(3), key=lambda j: v[j])


def fast_j(x):
 """j by the rule of proof.md 3.1 (bit-exact model of the arrival kernels), or None to fall back."""
 xl = [(x >> (64 * i)) & M64 for i in range(4)]
 F1 = sum(((xl[i] * GL[5 - i]) & M64) + ((xl[i] * GL[4 - i]) >> 64) for i in range(4)) & M64
 F0 = xl[3]
 S, d = (F0 + F1) & M64, (F0 - F1) & M64
 m = F1 if F0 >= F1 else F0
 e = m + S  # e >= 2^64 iff m >= F2*
 if F1 + 16 > M64 or S + 17 > M64 or (e & M64) + 33 > M64 or d < 16:
  return None
 return 2 if e > M64 else 1 if F0 >= F1 else 0


def stress_x(rr):
 """Near ties (two of v0, v1, v2 within 2^k, k = 150..249), a value near 0, or x + beta*x near p."""
 q, b = PROD["q"], PROD["beta"]
 mlt = rr.choice([b - 1, b * b - 1, b - b * b, 1 + b, b, 1, b * b])
 k = rr.randrange(150, 250)
 d = rr.randrange(-(1 << k), 1 << k)
 return (q - d if mlt == 1 + b else d) * pow(mlt % q, -1, q) % q or 1


def low_r_x(rr):
 """x with x mod 2^192 > 2^192 - C0: the only x with 2^64*x/p >= F0 + 1 (proof.md section 3.1)."""
 return (rr.randrange(PROD["q"] >> 192) << 192) + 2 ** 192 - rr.randrange(1, 2 ** 32 + 977)


def check_fast_j(n_random, n_stress):
 rr = random.Random(seed + 10 ** 7)  # separate stream: Part B is unaffected by this check
 res = {}
 for kind, n, gen in (("random", n_random, lambda: rr.randrange(1, PROD["q"])),
      ("stress", n_stress, lambda: stress_x(rr)), ("low_r", n_stress // 10, lambda: low_r_x(rr))):
  res.update({kind + "_requested": n, kind + "_done": 0, kind + "_fallbacks": 0, kind + "_mismatches": 0})
  for _ in range(n):
   if time.time() - T0 > budget:  # stay inside the runner's time limit
    break
   x = gen()
   res[kind + "_done"] += 1
   f = fast_j(x)
   if f is None:
    res[kind + "_fallbacks"] += 1
   elif f != exact_j(x):
    res[kind + "_mismatches"] += 1
 return res


# -- Part A: production checks
def part_a(nchecks):
 out = {}
 E = Curve(PROD["q"], PROD["r"], PROD["h"])
 q, r = E.q, E.r
 E.beta, E.lam = PROD["beta"], PROD["lam"]
 G = (PROD["gx"], PROD["gy"])
 E.G = G
 out["P0_on_curve"] = E.on_curve(G)
 out["r_P0_is_infinity"] = E.mul(r - 1, G) == E.neg(G)
 out["q_mod_3"] = q % 3
 out["beta_cube_root_of_unity"] = (E.beta * E.beta + E.beta + 1) % q == 0
 out["lambda_cube_root_of_unity"] = (E.lam * E.lam + E.lam + 1) % r == 0
 out["phi_P0_equals_lambda_P0"] = E.mul(E.lam, G) == E.phi(G)
 out["cofactor_1_r_prime_order_group"] = PROD["h"] == 1
 out["log2_r"] = math.log2(r)
 out["log2_N_classes"] = math.log2((r - 1) // 6)
 rm1, fac = r - 1, {}
 for pr in [2, 3, 149, 631, 107361793816595537, 174723607534414371449, 341948486974166000522343609283189]:
  while rm1 % pr == 0:
   fac[pr] = fac.get(pr, 0) + 1
   rm1 //= pr
 order = r - 1
 for pr, ex in fac.items():
  for _ in range(ex):
   if pow(2, order // pr, r) == 1:
    order //= pr
   else:
    break
 out["r_minus_1_fully_factored"] = rm1 == 1 and all(pow(2, order // pr, r) != 1 for pr in fac if order % pr == 0)
 out["log2_order_of_2_mod_r"] = round(math.log2(order), 4)
 xs = rng.randrange(1, r)
 T = E.mul(xs, G)
 W = Walker(E, 4, 0, T, rng)  # tiny table, for identity checks only
 ok_orbit = ok_key = ok_canon = ok_twist = ok_double = True
 for _ in range(nchecks * 4):
  k = rng.randrange(1, r)
  P = E.mul(k, G)
  members = [E.mul(u * k % r, G) for u in W.units]
  for idx, u in enumerate(W.units):
   s, j = divmod(idx, 3)
   Pj = E.phi(P, j)
   if s:
    Pj = E.neg(Pj)
   if members[idx] != Pj:
    ok_orbit = False
  if len(set(members)) != 6:
   ok_orbit = False
  cans = [W.canon(Mm) for Mm in members]
  if len({c[0] for c in cans}) != 1 or cans[0][0] != (q - P[1] if P[1] & 1 else P[1]) or cans[0][0] & 1:
   ok_key = False
  if len({c[2] for c in cans}) != 1:
   ok_canon = False
  for Mm, (K, t, c) in zip(members, cans):
   if E.mul(W.units[t], Mm) != c:
    ok_canon = False
   h = (K >> 1) & (W.R - 1)
   TP = W.table[h][t][0]
   Mh = W.table[h][0][0]
   n1 = W.canon(E.add(Mm, TP))
   n2 = W.canon(E.add(c, Mh))
   if n1[2] != n2[2]:
    ok_twist = False
   if W.canon(E.add(Mm, Mm))[2] != W.canon(E.add(c, c))[2]:
    ok_double = False
 out["orbit_size_6_and_unit_action"] = ok_orbit
 out["key_even_y_is_class_invariant"] = ok_key
 out["canonical_member_unique"] = ok_canon
 out["twisted_step_is_class_map"] = ok_twist
 out["doubling_is_class_map"] = ok_double
 ok_tab = all(W.table[h][t][0] == E.mul(W.uinv[t], W.table[h][0][0]) for h in range(W.R) for t in range(6))
 out["twisted_table_entries_equal_u_inverse_M"] = ok_tab

 Wp = Walker(E, 16, 64, T, rng)  # small table (R does not affect these identities)
 a, b = rng.randrange(r), rng.randrange(r)
 P = E.add(E.mul(a, G), E.mul(b, T))
 steps = 6 * nchecks
 trailA = []
 for i in range(steps):
  K, t, c = Wp.canon(P)
  TP, tc, td = Wp.step_target(P, K, t)
  P = E.add(P, TP)
  a, b = (a + tc) % r, (b + td) % r
  trailA.append((P, a, b))
 out["production_steps_coefficients_ok"] = E.add(E.mul(a, G), E.mul(b, T)) == P
 i0 = steps // 3
 PA, aA, bA = trailA[i0]
 v = rng.randrange(6)
 u = Wp.units[v]
 bB = rng.randrange(r)
 aB = (u * (aA + bA * xs) - bB * xs) % r
 PB = E.mul(u, PA)
 assert E.add(E.mul(aB, G), E.mul(bB, T)) == PB
 for i in range(i0 + 1, steps):
  K, t, c = Wp.canon(PB)
  TP, tc, td = Wp.step_target(PB, K, t)
  PB = E.add(PB, TP)
  aB, bB = (aB + tc) % r, (bB + td) % r
 PA, aA, bA = trailA[-1]
 KA, tA, cA = Wp.canon(PA)
 KB, tB, cB = Wp.canon(PB)
 out["merged_trails_same_canonical"] = cA == cB
 uA, uB = Wp.units[tA], Wp.units[tB]
 a1, b1 = uA * aA % r, uA * bA % r  # canonical coefficients of trail A
 a2, b2 = uB * aB % r, uB * bB % r
 xsol = (a2 - a1) * pow(b1 - b2, -1, r) % r
 out["solve_x_equals_(aB-aA)/(bA-bB)"] = xsol == xs and E.mul(xsol, G) == T
 nr = 15000 * nchecks if canon_random is None else canon_random
 ns = 5000 * nchecks if canon_stress is None else canon_stress
 out["fast_canonical_index"] = check_fast_j(nr, ns)
 out["part_a_seconds"] = round(time.time() - T0, 2)
 return out


# -- Part B: toy run
def toy_trial(E, rng, R, dp_bits):
 r = E.r
 xs = rng.randrange(1, r)
 T = E.fmul(xs)

 def point(a, b):
  return E.fmul(a + b * xs)

 W = Walker(E, R, dp_bits, T, rng, point)
 N = (r - 1) // CLS
 first_seen = {}  # class key -> trail id (instrumentation only)
 state = {"visits": 0, "tc": None, "trail": None}

 def on_class(K):
  state["visits"] += 1
  tid = state["trail"]
  prev = first_seen.get(K)
  if prev is None:
   first_seen[K] = tid
  elif prev != tid and state["tc"] is None:
   state["tc"] = state["visits"]

 dp_table = {}
 trail_id = 0
 solved = None
 Np = N * R / (R - 1)
 # toy only: + the mean own-trail age (1 - theta)/theta, which delays T_c on a toy (proof.md 5.2)
 B = KAPPA * math.sqrt(2 * Np * math.log(1 / (1 - P0))) + (1 << dp_bits) - 1
 D_cap = int(2 * (B + W.Lmax) / (1 << dp_bits)) + 1
 NODES = cfg.get("nodes", 16)
 node_shift = W.rbits + dp_bits + 1
 node_cap = -(-2 * D_cap // NODES)
 node_load = [0] * NODES
 matches = 0
 while solved is None:
  if pipeline and state["visits"] >= B:
   solved = False  # phase-1 budget exhausted and the open trail is finished
   break
  a, b = rng.randrange(r), rng.randrange(r)
  P = point(a, b)
  if P is None:
   continue
  trail_id += 1
  state["trail"] = trail_id
  ph = rng.randrange(PERIOD)
  v0 = state["visits"]
  status, K, length, c, _ = W.run_trail(P, None, on_class, ph)
  if status != "dp":
   W.stats["abandoned_visits"] = W.stats.get("abandoned_visits", 0) + state["visits"] - v0
   continue
  if K in dp_table:
   matches += 1
   if pipeline and matches > 8:
    solved = False
    break
   (a2, b2, l2, ph2) = dp_table[K]
   s1 = W.stats.copy()
   _, K1, _, c1, (A1, B1) = W.run_trail(point(a, b), (a, b), None, ph)
   _, K2, _, c2, (A2, B2) = W.run_trail(point(a2, b2), (a2, b2), None, ph2)
   W.stats.update(s1)  # re-walk cost is accounted separately
   if c1 == c2 and (B1 - B2) % r:
    x = (A2 - A1) * pow(B1 - B2, -1, r) % r
    if E.fmul(x) == T:
     solved = (x == xs)
     break
  else:
   if pipeline and len(dp_table) >= D_cap:
    solved = False
    break
   if pipeline:
    nd = (K >> node_shift) % NODES
    if node_load[nd] >= node_cap:
     solved = False
     W.stats["node_overflow"] = W.stats.get("node_overflow", 0) + 1
     break
    node_load[nd] += 1
   dp_table[K] = (a, b, length, ph)
  if time.time() - T0 > budget:
   return None
 nl = dict(max_node_load_over_mean=(max(node_load) * NODES / max(1, sum(node_load))) if pipeline and sum(node_load) else None,
   max_node_load_over_cap=(max(node_load) / node_cap) if pipeline else None)
 return dict(tc=state["tc"], solve_visits=state["visits"], ok=solved, trails=trail_id, **nl, **W.stats)


def long_trail_test(E, rng, R, L, ntrails):
 """Trails of L operations without distinguished points; counts returns other than fruitless-cycle repeats."""
 r = E.r
 xs = rng.randrange(1, r)
 T = E.fmul(xs)
 W = Walker(E, R, 0, T, rng, lambda c, d: E.fmul(c + d * xs))
 W.is_dp = lambda K: False
 W.Lmax = L
 out = dict(trails=0, operations=0, escapes=0, nonfruitless_returns=0, fruitless_repeats=0,
   trails_with_nonfruitless_returns=0)
 N = (r - 1) // 6
 for _ in range(ntrails):
  if time.time() - T0 > budget:
   break
  P = E.fmul(rng.randrange(1, r))
  seen = {}
  cnt = [0, 0]  # [non-fruitless returns, fruitless-cycle repeats]
  arr = [0]

  def on_class(K):
   ep = W.stats["doublings"]
   arr[0] += 1
   if K in seen:
    pep, parr = seen[K]
    if pep < ep or arr[0] - parr > TRAV_CAP:
     cnt[0] += 1
    else:
     cnt[1] += 1
   seen[K] = (ep, arr[0])

  s0 = dict(W.stats)
  W.run_trail(P, None, on_class)
  out["trails"] += 1
  out["operations"] += W.stats["steps"] + W.stats["doublings"] - s0["steps"] - s0["doublings"]
  out["escapes"] += W.stats["doublings"] - s0["doublings"]
  out["nonfruitless_returns"] += cnt[0]
  out["fruitless_repeats"] += cnt[1]
  out["trails_with_nonfruitless_returns"] += cnt[0] > 0
 out["expected_self_collisions_per_trail_L2_over_2N"] = round(L * L / (2 * N), 4)
 return out


def rule_check(E, rng, R, dp_bits, ntrails):
 """proof.md 5.1: each start walked under the rule (1, 256) (compare at every arrival, anchors every 256) and
  under the production rule (STRIDE, PERIOD), with the same phase and a trail cap of 2^15. Whenever both end at a
  distinguished class it must be the same class and canonical point after the same number of escapes; the
  operations are summed over those trails."""
 xs = rng.randrange(1, E.r)
 W = Walker(E, R, dp_bits, E.fmul(xs), rng, lambda c, d: E.fmul(c + d * xs))
 W.Lmax = 1 << 15
 out = dict(R=R, trails=0, both_dp=0, same_end=0, same_escapes=0, escapes=0, ops_1_256=0, ops_production=0)
 for _ in range(ntrails):
  P, ph, ends = E.fmul(rng.randrange(1, E.r)), rng.randrange(PERIOD), []
  for rule in ((1, 256), (STRIDE, PERIOD)):
   e0 = W.stats["esccp"]
   st, K, length, c, _ = W.run_trail(P, None, None, ph, rule)
   ends.append((st, K, c, W.stats["esccp"] - e0, length))
  out["trails"] += 1
  if ends[0][0] == ends[1][0] == "dp":
   out["both_dp"] += 1
   out["same_end"] += ends[0][:3] == ends[1][:3]
   out["same_escapes"] += ends[0][3] == ends[1][3]
   out["escapes"] += ends[0][3]
   out["ops_1_256"] += ends[0][4]
   out["ops_production"] += ends[1][4]
 return out

# -- kernels (proof.md section 2)
M64 = (1 << 64) - 1
C0 = (1 << 32) + 977  # 2^256 = C0 (mod p)
KST = {}  # statistics of the limb ports (extreme values met)


def _kmax(key, v):
 if v > KST.get(key, -1):
  KST[key] = v


def _addc(a, b, c):
 s = a + b + c
 return s & M64, s >> 64


def _subc(a, b, c):
 d = a - b - c
 return d & M64, 1 if d < 0 else 0


def _inc(a, c):
 """a + c where the kernel adds a carry bit without a carry out (asserted)."""
 assert a + c <= M64
 return a + c


def _limbs(v):
 return [(v >> (64 * i)) & M64 for i in range(4)]


def _val(l):
 return sum(x << (64 * i) for i, x in enumerate(l))


def _lohi(x, y):
 p = x * y
 return p & M64, p >> 64


# the six carry chains of fe_mul (kernels/gen_loops.py MUL_PLAN): (op, destination, sources)
MUL_PLAN = [
 ("adds", "s1", "H00", "L01"), ("adcs", "s2", "H01", "L02"), ("adcs", "s3", "H02", "L03"),
 ("adcs", "s4", "H03", "L13"), ("adc", "c5a", "H13"),
 ("adds", "z1", "s1", "L10"), ("adcs", "s2", "s2", "H10"), ("adcs", "s3", "s3", "H20"),
 ("adcs", "s4", "s4", "H30"), ("adc", "c5b", "H22"),
 ("adds", "s2", "s2", "L20"), ("adcs", "s3", "s3", "L30"), ("adcs", "s4", "s4", "L22"),
 ("adcs", "s5", "c5a", "c5b"), ("adc", "c6a", "H23"),
 ("adds", "z2", "s2", "L11"), ("adcs", "s3", "s3", "H11"), ("adcs", "s4", "s4", "H12"),
 ("adcs", "s5", "s5", "H31"), ("adc", "c6b", "H32"),
 ("adds", "s3", "s3", "L12"), ("adcs", "s4", "s4", "H21"), ("adcs", "s5", "s5", "L23"),
 ("adcs", "s6", "c6a", "c6b"), ("adc", "c7", "H33"),
 ("adds", "z3", "s3", "L21"), ("adcs", "z4", "s4", "L31"), ("adcs", "z5", "s5", "L32"),
 ("adcs", "z6", "s6", "L33"), ("adc", "z7", "c7"),
]


def p_mul(a, b, subs=None, add=None):
 """Port of the generated fe_mul: 32 products, six carry chains, then the reduction (fused if subs or add)."""
 v = {}
 for i in range(4):
  for j in range(4):
   v["L%d%d" % (i, j)], v["H%d%d" % (i, j)] = _lohi(a[i], b[j])
 v["z0"] = v.pop("L00")
 c = 0
 for st in MUL_PLAN:
  op, dst, src = st[0], st[1], [v.pop(w) for w in st[2:]]
  if op == "adc":
   v[dst] = _inc(src[0], c)  # a fresh high word (<= 2^64 - 2) or the top column
  else:
   v[dst], c = _addc(src[0], src[1], c if op == "adcs" else 0)
 assert sorted(v) == ["z%d" % k for k in range(8)]
 z = [v["z%d" % k] for k in range(8)]
 assert _val(z) == _val(a) * _val(b)
 return _reduce_sub(z, subs or [], add) if subs or add else _reduce(z)


def _reduce(z):
 """Port of the generated reduction (23 instructions); returns t mod 2^256, t = (s mod 2^256) + s4*C0
    (proof.md 2.2)."""
 g, k = zip(*[_lohi(z[4 + j], C0) for j in range(4)])
 q1, c = _addc(g[1], k[0], 0)
 q2, c = _addc(g[2], k[1], c)
 q3, c = _addc(g[3], k[2], c)
 q4 = _inc(k[3], c)
 s0, c = _addc(z[0], g[0], 0)
 s1, c = _addc(z[1], q1, c)
 s2, c = _addc(z[2], q2, c)
 s3, c = _addc(z[3], q3, c)
 s4 = _inc(q4, c)
 _kmax("s4", s4)
 assert s4 <= C0
 w0, w1 = _lohi(s4, C0)
 _kmax("w1", w1)
 assert w1 <= 1
 t0, c = _addc(s0, w0, 0)
 t1, c = _addc(s1, w1, c)
 t2, c = _addc(s2, 0, c)
 t3, c = _addc(s3, 0, c)
 if c:
  # t >= 2^256, so z mod p = t - p < C0^2 + C0 and the result t - 2^256 is z mod p - C0
  KST["final_carry"] = KST.get("final_carry", 0) + 1
  zp = _val(z) % PROD["q"]
  assert zp < C0 * C0 + C0 and _val([t0, t1, t2, t3]) == zp - C0
 return [t0, t1, t2, t3]


def _red_ok(z, r):
 """r equals the port's reduction of z, which is z mod p up to weak form or, on the rare event, z mod p - C0."""
 m = _val(_reduce(_limbs(z % 2**256) + _limbs(z >> 256)))
 zp = z % PROD["q"]
 return r == m and (m % PROD["q"] == zp or (zp < C0 * C0 + C0 and m == zp - C0))


def p_sqr(a, subs=None, bias=False):
 """Port of the generated fe_sqr: off-diagonal chains, doubling, diagonal chain, reduction."""
 a0, a1, a2, a3 = a
 l01, h01 = _lohi(a0, a1); l02, h02 = _lohi(a0, a2); l03, h03 = _lohi(a0, a3)
 l12, h12 = _lohi(a1, a2); l13, h13 = _lohi(a1, a3); l23, h23 = _lohi(a2, a3)
 y3, c = _addc(l03, l12, 0)
 y4, c = _addc(h03, h12, c)
 h13 = _inc(h13, c)
 o1 = l01
 o2, c = _addc(h01, l02, 0)
 o3, c = _addc(y3, h02, c)
 o4, c = _addc(y4, l13, c)
 o5, c = _addc(h13, l23, c)
 o6 = _inc(h23, c)
 d7 = o6 >> 63
 d6 = ((o6 << 1) & M64) | (o5 >> 63)
 d5 = ((o5 << 1) & M64) | (o4 >> 63)
 d4 = ((o4 << 1) & M64) | (o3 >> 63)
 d3 = ((o3 << 1) & M64) | (o2 >> 63)
 d2 = ((o2 << 1) & M64) | (o1 >> 63)
 z0, g0 = _lohi(a0, a0); y2, g1 = _lohi(a1, a1); y4d, g2 = _lohi(a2, a2); y6, g3 = _lohi(a3, a3)
 z1, c = _addc(g0, (o1 << 1) & M64, 0)
 z2, c = _addc(d2, y2, c)
 z3, c = _addc(d3, g1, c)
 z4, c = _addc(d4, y4d, c)
 z5, c = _addc(d5, g2, c)
 z6, c = _addc(d6, y6, c)
 z7 = d7 + g3 + c
 assert z7 <= M64  # the top word of a^2 < 2^512
 z = [z0, z1, z2, z3, z4, z5, z6, z7]
 assert _val(z) == _val(a) ** 2
 return _reduce_sub(z, subs, None, bias) if subs else _reduce(z)


def p_sub(a, b):
 """fe_sub: a any 256-bit value, b canonical."""
 assert _val(b) < PROD["q"]
 d0, bw = _subc(a[0], b[0], 0)
 d1, bw = _subc(a[1], b[1], bw)
 d2, bw = _subc(a[2], b[2], bw)
 d3, bw = _subc(a[3], b[3], bw)
 m = C0 if bw else 0
 r0, bw2 = _subc(d0, m, 0)
 r1, bw2 = _subc(d1, 0, bw2)
 r2, bw2 = _subc(d2, 0, bw2)
 r3, bw2 = _subc(d3, 0, bw2)
 assert bw2 == 0  # no second borrow: the wrapped value is above C0 whenever b is canonical
 return [r0, r1, r2, r3]


def _rx(v, c, k):
 """Row stub (proof.md 2.5): a carry c out of limb 1 runs into limbs 2, 3; one out of 2^256 adds C0."""
 if c:
  KST[k] = KST.get(k, 0) + 1
  v = _val(v) + 2**128
  if v >> 256:
   KST[k + "w"] = KST.get(k + "w", 0) + 1
   v += C0 - 2**256
  v = _limbs(v)
 return v


def p_rows(R):
 """Port of pass 2's rows B_i = R*2^(64i) mod p: B_(i+1) = [lo, b0 + hi, b1, b2], (lo, hi) = b3*C0."""
 B = [R]
 for _ in range(3):
  b = B[-1]
  lo, hi = _lohi(b[3], C0)
  b0, c = _addc(b[0], hi, 0)
  B.append(_rx([lo, b0] + b[1:3], c, "row_shift"))
 return B


def p_rmul(B, a):
 """Port of a row product: Z = sum a_i*B_i (rows 4 + 6 additions, z5 <= 3), then (Z mod 2^256) +
    (z4 + 2^64*z5)*C0 by MUL, UMULH, MADD, ADDS, ADCS and the row stub."""
 z = [0] * 6
 for i in range(4):
  l, h = zip(*[_lohi(a[i], B[i][j]) for j in range(4)])
  r, c = [l[0]], 0
  for j in range(3):
   x, c = _addc(h[j], l[j + 1], c)
   r.append(x)
  r.append(_inc(h[3], c))
  c = 0
  for k in range(5):
   z[k], c = _addc(z[k], r[k], c)
  z[5] = _inc(z[5], c)
 assert _val(z) == sum(a[i] * _val(B[i]) for i in range(4)) and z[5] <= 3
 lo, hi = _lohi(z[4], C0)
 t0, c = _addc(z[0], lo, 0)
 t1, c = _addc(z[1], _inc(hi, z[5] * C0), c)
 return _rx([t0, t1] + z[2:4], c, "row_fold")


def _reduce_sub(z, bs, add=None, bias=False):
 """Port of the fused reduction (24 + 5k instructions): s = L + H*C0, the subtrahends bs off all five limbs,
    then the signed fold.  add: y3' (23 + 5, unsigned fold); bias: x3 (23 + 10, chains end ADC, unsigned fold)."""
 g, k = zip(*[_lohi(z[4 + j], C0) for j in range(4)])
 q1, c = _addc(g[1], k[0], 0)
 q2, c = _addc(g[2], k[1], c)
 q3, c = _addc(g[3], k[2], c)
 s4 = _inc(k[3], c)
 s0, c = _addc(z[0], g[0], 0)
 s1, c = _addc(z[1], q1, c)
 s2, c = _addc(z[2], q2, c)
 s3, c = _addc(z[3], q3, c)
 s4 = _inc(s4, c)
 assert s4 <= C0
 for b in bs:
  assert _val(b) < (2**256 if bias else PROD["q"])
  s0, w = _subc(s0, b[0], 0)
  s1, w = _subc(s1, b[1], w)
  s2, w = _subc(s2, b[2], w)
  s3, w = _subc(s3, b[3], w)
  s4 += bias - w
 if add:
  assert _val(add) < PROD["q"]
  s0, c = _addc(s0, add[0], 0)
  s1, c = _addc(s1, add[1], c)
  s2, c = _addc(s2, add[2], c)
  s3, c = _addc(s3, add[3], c)
  s4 += c
  assert (s4 * C0) >> 64 <= 1
 assert (0 if bias else -2) <= s4 <= C0 + (2 if bias else 1 if add else 0)
 assert (s4 * C0) >> 64 <= 1 or s4 < 0
 t = _val([s0, s1, s2, s3]) + s4 * C0
 q = PROD["q"]
 v = (_val(z) - sum(_val(b) for b in bs) + _val(add or [0]) + 2 * C0 * bias) % q
 if t >= 2**256:  # as in _reduce: v < C0^2 + 3*C0 and the result is v - C0
  KST["fused_carry"] = KST.get("fused_carry", 0) + 1
  assert v < C0 * C0 + 3 * C0 and t - 2**256 == v - C0
 elif t < 0:  # two subtrahends, s4' = -1 and s' mod 2^256 < C0: the result is v + C0 (mod p)
  KST["fused_negative"] = KST.get("fused_negative", 0) + 1
  assert len(bs) == 2 and s4 == -1 and (t + 2**256 - v - C0) % q == 0
 return _limbs(t % 2**256)


def _fused_ok(z, bs, r, add=0, bias=False):
 m = _val(_reduce_sub(_limbs(z % 2**256) + _limbs(z >> 256), [_limbs(b) for b in bs], add and _limbs(add), bias))
 q = PROD["q"]
 v = (z - sum(bs) + add + 2 * C0 * bias) % q
 return r == m and (m % q == v or m == v - C0 or m % q == (v + C0) % q)


def _negative_inputs(rr, n):
 """(a, b, c) with a = 2^128, so a^2 = 2^256 gives s = C0, and b + c in (2^256, 2^256 + C0]: the fold goes
    negative (proof.md section 2.2)."""
 out = []
 for _ in range(n):
  e, d = rr.randrange(2**32), rr.randrange(1, C0 + 1)
  out.append((2**128, PROD["q"] - 1 - e, C0 + 1 + e + d))
 return out


def _edge_values():
 q = PROD["q"]
 vals = [0, 1, 2, C0, q - 1, q, q + 1, 2**256 - 1, 2**256 - 2, 2**256 - C0, 2**256 - C0 - 1,
   2**255, 2**255 - 1, 2**64 - 1, 2**64, 2**128 - 1, 2**192 - 1, 2**192, 2**33, q - 2**33]
 for k in (32, 40, 64, 128, 192, 250):
  vals += [q - 2**k, (q + 2**k) % 2**256, 2**256 - 2**k]
 for pat in range(16):
  vals.append(sum(((M64 if pat >> i & 1 else 0) << (64 * i)) for i in range(4)))
 return sorted(set(v % 2**256 for v in vals))


def _second_carry_inputs(rr, count):
 """512-bit z whose reduction's last addition carries out of 2^256."""
 out = []
 while len(out) < count:
  H = rr.getrandbits(256)
  if rr.random() < 0.5:
   H |= (2**256 - 1) ^ (2**rr.randrange(1, 64) - 1)
  HC = H * C0
  t = 2**256 - 1 - rr.randrange(0, 2**rr.randrange(1, 66))
  L = (t - HC) % 2**256
  s4 = (L + HC) >> 256
  if t + s4 * C0 >= 2**256:
   out.append(L + (H << 256))
  if len(out) < count and rr.random() < 0.25:  # exact boundary t + s4*C0 == 2^256
   t = 2**256 - s4 * C0
   L = (t - HC) % 2**256
   if (L + HC) >> 256 == s4:
    out.append(L + (H << 256))
 return out


def _kernel_port_check(nvec):
 q = PROD["q"]
 rr = random.Random(777)
 edge = _edge_values()
 vals = edge + [rr.getrandbits(256) for _ in range(nvec)] + [q - rr.getrandbits(40) for _ in range(nvec // 8)] \
  + [2**256 - 1 - rr.getrandbits(40) for _ in range(nvec // 8)] + [rr.getrandbits(64) << (64 * rr.randrange(4)) for _ in range(nvec // 8)]
 bad = 0
 n = 0
 for a in edge:
  for b in edge:
   r = _val(p_mul(_limbs(a), _limbs(b)))
   n += 1
   bad += not _red_ok(a * b, r)
 for _ in range(nvec):
  a, b = rr.choice(vals), rr.choice(vals)
  r = _val(p_mul(_limbs(a), _limbs(b)))
  n += 1
  bad += not _red_ok(a * b, r)
 ns = 0
 for a in vals:
  r = _val(p_sqr(_limbs(a)))
  ns += 1
  bad += not _red_ok(a * a, r)
 nsub = 0
 canon_edge = [v for v in edge if v < q]
 pairs = [(a, b) for a in edge for b in canon_edge] + [(rr.choice(vals), rr.choice(vals) % q) for _ in range(nvec)]
 for a, b in pairs:
  r = _val(p_sub(_limbs(a), _limbs(b)))
  nsub += 1
  bad += not (r < 2**256 and r % q == (a - b) % q)
  bad += not (r < q or (a >= q and b < C0))  # canonical unless a >= p and b < C0
  KST["sub_noncanonical"] = KST.get("sub_noncanonical", 0) + (r >= q)
 trip = [(a, b, c) for a in edge[::3] for b in canon_edge[::2] for c in canon_edge[::3]] + _negative_inputs(rr, 50)
 trip += [(rr.choice(vals), rr.randrange(q), rr.randrange(q)) for _ in range(nvec // 4)]
 trip += [(rr.choice(vals), rr.choice([0, q - 1, rr.randrange(q)]), q - C0 - 1) for _ in range(40)]  # x_T + 2C0 = 2^256 - 1
 nf = 0
 for a, b, c in trip:
  if c + 2 * C0 < 2**256:  # x3 of the loops: subtrahends x_W and x_T + 2C0 (proof.md 2.2)
   nf += 1
   bad += not _fused_ok(a * a, [b, c + 2 * C0], _val(p_sqr(_limbs(a), [_limbs(b), _limbs(c + 2 * C0)], True)), 0, True)
  bad += not _fused_ok(a * a, [b, c], _val(p_sqr(_limbs(a), [_limbs(b), _limbs(c)])))
  bad += not _fused_ok(a * b, [c], _val(p_mul(_limbs(a), _limbs(b), [_limbs(c)])))
  bad += not _fused_ok(a * b, [], _val(p_mul(_limbs(a), _limbs(b), add=_limbs(c))), c)
 zs = [rr.getrandbits(512) for _ in range(nvec)] + _second_carry_inputs(rr, 200) + [2**512 - 1, 0, 2**256 - 1]
 for z in zs:
  r = _val(_reduce(_limbs(z % 2**256) + _limbs(z >> 256)))
  bad += not _red_ok(z, r)
 a0 = [rr.randrange(2**63, 2**64) for _ in range(8)]  # a0*R in [T, T + a0): the folds wrap 2^256
 rows = [(R, a) for R in edge[::2] + vals[99:139] for a in edge[1::3] + vals[-9:]]
 rows += [(-(-(2**318 + 2**256 - 2**65) // a), a) for a in a0] + [(2**192 - 1 + (a << 192), a) for a in a0]
 for R, a in rows:
  r = _val(p_rmul(p_rows(_limbs(R)), _limbs(a)))
  bad += not (r < 2**256 and r % q == R * a % q)
 return dict(mul_checked=n, rows_checked=len(rows), sqr_checked=ns, sub_checked=nsub, fused_checked=3 * len(trip) + nf, reduce_checked=len(zs),
    mismatches=bad, extremes=dict(KST))


def aadd(A, B, q=PROD["q"]):
 """Affine addition (doubling when the x agree) on production secp256k1."""
 if A[0] == B[0]:
  lam = 3 * A[0] * A[0] * pow(2 * A[1], -1, q) % q
 else:
  lam = (B[1] - A[1]) * pow(B[0] - A[0], -1, q) % q
 x = (lam * lam - A[0] - B[0]) % q
 return (x, (lam * (A[0] - x) - A[1]) % q)


def amul(k, A):
 R = None
 while k:
  if k & 1:
   R = A if R is None else aadd(R, A)
  A = aadd(A, A)
  k >>= 1
 return R


def _batch_port_check(nwalks):
 """Port of pass1/pass2 (shared inversion by pow) on production secp256k1 points."""
 q = PROD["q"]
 rr = random.Random(4242)
 Gp = (PROD["gx"], PROD["gy"])

 tab = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(4)]
 walks = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(nwalks)]
 tps = []
 for w in walks:
  t = tab[rr.randrange(4)]
  if t[0] == w[0]:
   t = tab[(tab.index(t) + 1) % 4]
  if rr.random() < 0.5:
   t = (t[0], q - t[1])
  tps.append(t)
 W = [(_limbs(x), _limbs(y)) for x, y in walks]
 T = [(_limbs(x), _limbs(y)) for x, y in tps]
 n = len(W)
 acc = [None] * n
 acc[0] = p_sub(T[0][0], W[0][0])
 for i in range(1, n):
  dx = p_sub(T[i][0], W[i][0])
  acc[i] = p_mul(dx, acc[i - 1])
 R = _limbs(pow(_val(acc[n - 1]), -1, q))

 def finish(w, t, inv):
  dy = p_sub(t[1], w[1])
  lam = p_mul(dy, inv)
  x3 = p_sqr(lam, [w[0], _limbs(_val(t[0]) + 2 * C0)], True)
  return (x3, p_mul(lam, p_sub(x3, w[0]), add=w[1]))  # -(W + T): proof.md 2.5

 out = [None] * n
 for i in range(n - 1, 0, -1):
  B = p_rows(R)
  inv, R = p_rmul(B, acc[i - 1]), p_rmul(B, p_sub(T[i][0], W[i][0]))
  out[i] = finish(W[i], T[i], inv)
 out[0] = finish(W[0], T[0], R)
 bad = 0
 for i in range(n):
  want = aadd(walks[i], tps[i])
  want = (want[0], q - want[1])
  got = (_val(out[i][0]), _val(out[i][1]))
  bad += got != want
 return dict(walks=n, mismatches=bad)


def asm_exec_check(text, sizes=(1, 3, 6, 129)):
 """Run pass1/pass2 of the shipped loops on production points against affine arithmetic; 1 and 129 walks
    give the per-walk count, 3 and 6 run the tails."""
 q = PROD["q"]
 rr = random.Random(2468)
 Gp = (PROD["gx"], PROD["gy"])

 prog, labels = a64_load(text)
 tab = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(3)]
 tab += [(x, q - y) for x, y in tab]
 res = dict(walks_checked=0, mismatches=0, executed={})
 for n in sizes:
  walks = [amul(rr.randrange(1, PROD["r"]), Gp) for _ in range(n)]
  tps = [tab[rr.randrange(len(tab))] for _ in range(n)]
  mem = {}
  WB, TB, PB, RB, SB = 0x10000, 0x20000, 0x30000, 0x50000, 0x90000
  AB = (1 << 27) + (1 << 26) - 64  # the slot array
  for i in range(n):
   for k, v in enumerate(_limbs(walks[i][0]) + _limbs(walks[i][1])):
    mem[WB + 64 * i + 8 * k] = v
   for k, v in enumerate(_limbs(tps[i][0]) + _limbs(tps[i][1]) + [0] * 16 + _limbs(tps[i][0] + 2 * C0)):
    mem[TB + 256 * i + 8 * k] = v  # x, y and, at +192, x + 2C0
   mem[PB + 8 * i] = TB + 256 * i
  for a in range(0, 80 + 64, 8):
   mem[SB - 112 + a] = 0
  o1, o2 = {}, {}
  s1 = a64_run(prog, labels, "pass1", {"x0": n, "x1": WB, "x2": PB, "x3": AB, "x19": 10 ** 9, "sp": SB}, mem,
     out_regs=o1)
  last = sum(mem[AB + 64 * n + 32 + 8 * k] << (64 * k) for k in range(4))
  rinv = pow(last, -1, q)
  for k, v in enumerate(_limbs(rinv)):
   mem[RB + 8 * k] = v
  s2 = a64_run(prog, labels, "pass2", {"x0": n, "x1": WB, "x2": PB, "x3": AB, "x4": RB, "x19": 10 ** 9,
           "sp": SB}, mem, out_regs=o2)
  res["exits_taken"] = res.get("exits_taken", 0) + (2 * 10 ** 9 - o1["x19"] - o2["x19"]) // 20
  for i in range(n):
   got = (sum(mem[WB + 64 * i + 8 * k] << (64 * k) for k in range(4)),
    sum(mem[WB + 64 * i + 32 + 8 * k] << (64 * k) for k in range(4)))
   res["walks_checked"] += 1
   want = aadd(walks[i], tps[i])
   res["mismatches"] += got != (want[0], q - want[1])  # the loops store -(W + T)
  res["executed"][n] = dict(pass1=s1, pass2=s2)
 ns = [n for n in sorted(res["executed"]) if n % 128 == 1]
 res["per_extra_walk"] = dict(
  pass1=[res["executed"][b]["pass1"] - res["executed"][a]["pass1"] for a, b in zip(ns, ns[1:])],
  pass2=[res["executed"][b]["pass2"] - res["executed"][a]["pass2"] for a, b in zip(ns, ns[1:])])
 res["per_extra_walk_divided"] = dict(
  pass1=[(res["executed"][b]["pass1"] - res["executed"][a]["pass1"]) / (b - a) for a, b in zip(ns, ns[1:])],
  pass2=[(res["executed"][b]["pass2"] - res["executed"][a]["pass2"]) / (b - a) for a, b in zip(ns, ns[1:])])
 return res


def arrival_check(text, nb=8):
 """Run kernels/arrive.s on production inputs (random, near ties, low r, forced distinguished keys and filter
    hits) against the bit-exact model and the exact rule; rare-free batches of 64 to 193 walks give the count."""
 prog, labels = a64_load(text)
 q, rr = PROD["q"], random.Random(seed + 31)
 TP, B0, WK, KC, LG, CS = 1 << 25, 1 << 35, 1 << 20, 1 << 21, 1 << 22, 1 << 23
 consts = GL[1:6] + [(q & M64) + 1, M64, 33]

 def model(x, y, kcp):
  K, j = (y if y % 2 == 0 else q - y), fast_j(x)
  code = 1 if (K >> 64) & M64 == 0 else 2 if j is None else 3 if (K & M64) + 1 == kcp else 0
  return code, B0 + 64 * (j or 0) + 512 * (y & 1) + 1024 * ((K >> 1) & (2 ** 25 - 1))

 def run(ws):
  mem = {CS + 8 * i: c for i, c in enumerate(consts)}
  for i, (x, y, kcp) in enumerate(ws):
   for k in range(4):
    mem[WK + 64 * i + 8 * k], mem[WK + 64 * i + 32 + 8 * k] = (x >> (64 * k)) & M64, (y >> (64 * k)) & M64
   mem[KC + 8 * i], mem[LG + 8 * i] = kcp, 0
  tp0 = TP + 2 ** 16 - 8 * len(ws)
  mem[1 << 24] = LG  # the test stub's log base, at the caller's [sp]
  n = a64_run(prog, labels, "arrive", dict(x0=len(ws), x1=WK, x2=tp0, x3=KC, x4=B0, x5=CS, sp=1 << 24), mem)
  return n, [(mem[LG + 8 * i], mem.get(tp0 + 8 * i)) for i in range(len(ws))]

 out = dict(walks=0, plain=0, distinguished=0, window=0, hit=0, stress=0, low_r=0, mismatches_vs_model=0,
   fast_indices=0, mismatches_vs_exact=0)
 for _ in range(nb):
  ws = []
  for _ in range(rr.choice([1, 2, 3, 64, 65, 250, 251])):
   kind = rr.choice("rrrrsssl")
   x = rr.randrange(1, q) if kind == "r" else stress_x(rr) if kind == "s" else low_r_x(rr)
   y = rr.randrange(1, q)
   if rr.random() < 0.05:  # a distinguished key: limb 1 of K zero
    K = (rr.getrandbits(256) & ~(M64 << 64)) % q & ~1 or 2
    y = K if rr.random() < 0.5 else q - K
   K = y if y % 2 == 0 else q - y
   ws.append((x, y, (K & M64) + 1 if rr.random() < 0.15 else rr.getrandbits(64) & ~1))
   out["stress"] += kind == "s"
   out["low_r"] += kind == "l"
  for (x, y, kcp), (code, addr) in zip(ws, run(ws)[1]):
   mc, ma = model(x, y, kcp)
   out["walks"] += 1
   out[("plain", "distinguished", "window", "hit")[mc]] += 1
   out["mismatches_vs_model"] += code != mc or (mc in (0, 3) and addr != ma)
   j = fast_j(x)
   if j is not None:
    out["fast_indices"] += 1
    out["mismatches_vs_exact"] += j != exact_j(x)
 counts = {}
 for n in (64, 65, 128, 129, 192, 193):
  ws = []
  while len(ws) < n:
   x, y = rr.randrange(1, q), rr.randrange(1, q) | (len(ws) in (0, n - 2, n - 1))
   if model(x, y, 0)[0] == 0:
    ws.append((x, y, 0))
  counts[n] = run(ws)[0]
 out["per_walk"] = sorted({(counts[b] - counts[a]) / (b - a) for a, b in ((64, 128), (128, 192), (65, 129), (129, 193))})
 return out


def fused_check(text):
 """Run both copies of kernels/secp_fuse.s (proof.md 2.6 item 7), pass1f and pass1p, on every batch against the
   address model, big-integer dx and prefix products, the rare list and the meter; rare-free batches of 1, 129
   and 130 give the counts.  pass1p (plain rounds) has no filter test, so a filter hit there is not rare.

  :param text: the source of kernels/secp_fuse.s
  :return: dict of counts and mismatches
  """
 prog, labels = a64_load(text)
 q, rr = PROD["q"], random.Random(seed + 47)
 WK, PB, AB, KC, RB, SB, CS, B0 = 1 << 20, 1 << 22, 1 << 23, 1 << 24, 1 << 25, 1 << 26, 1 << 27, 1 << 35
 ctx = GL[1:6] + [(q & M64) + 1, M64]
 q0 = q & M64

 def model(x, y, kcp, cmp=True):
  K, j = (y if y % 2 == 0 else q - y), fast_j(x)
  y0 = y & M64
  rare = (K >> 64) & M64 == 0 or j is None or (y % 2 == 1 and y0 > q0 + 1) or (cmp and (K & M64) + 1 == kcp)
  if not rare:
   return False, B0 + 64 * j + 512 * (y & 1) + 1024 * ((K >> 1) & (2 ** 25 - 1))
  xl = _limbs(x)
  F1 = sum(((xl[i] * GL[5 - i]) & M64) + ((xl[i] * GL[4 - i]) >> 64) for i in range(4)) & M64
  F0 = xl[3]
  jk = 2 if min(F0, F1) + ((F0 + F1) & M64) > M64 else 1 if F0 >= F1 else 0
  w = y0 if y % 2 == 0 else (q0 + 1 - y0) & M64
  return True, B0 + 64 * jk + ((w & (2 ** 26 - 1)) << 9)

 def sub(a, b):
  d = a - b
  return (d, 0) if d >= 0 else ((d + 2 ** 256 - C0) % 2 ** 256, int((d + 2 ** 256) % 2 ** 64 < C0))

 def red(z):
  t = (z % 2 ** 256) + (z >> 256) * C0
  s4, sm = t >> 256, t % 2 ** 256
  return (sm + s4 * C0) % 2 ** 256, int(sm % 2 ** 128 + s4 * C0 >= 2 ** 128)

 def run(walks, kcp, xts, fn):
  n, mem = len(walks), {CS + 8 * k: c for k, c in enumerate(ctx)}
  want = [model(x, y, c, fn == "pass1f") for (x, y), c in zip(walks, kcp)]
  for i, (x, y) in enumerate(walks):
   for k, v in enumerate(_limbs(x) + _limbs(y)):
    mem[WK + 64 * i + 8 * k] = v
   mem[KC + 8 * i] = kcp[i]
  tab = {a: xts[i] for i, (_, a) in enumerate(want)}
  for a, xt in tab.items():
   for k, v in enumerate(_limbs(xt)):
    mem[a + 8 * k] = v
  regs, meter = {}, 10 ** 9
  cnt = a64_run(prog, labels, fn, dict(x0=n, x1=WK, x2=PB, x3=AB, x4=B0, x5=KC, x6=CS, x7=RB, x19=meter,
                        sp=SB), mem, out_regs=regs)
  bad, exits, F = 0, [0] * n, None
  for i, ((x, y), (_, a)) in enumerate(zip(walks, want)):
   bad += mem[PB + 8 * i] != a
   dx, e = sub(tab[a], x)
   exits[i] += e
   if i:
    bad += _val([mem[AB + 64 * i + 8 * k] for k in range(4)]) != dx
    bad += _val([mem[AB + 64 * i + 32 + 8 * k] for k in range(4)]) != F
    F, e = red(dx * F)
    exits[i] += e
   else:
    F = dx
  bad += _val([mem[AB + 64 * n + 32 + 8 * k] for k in range(4)]) != F
  listed, p = [], RB
  while mem.get(p):
   listed.append((mem[p + 8] - WK) // 64 - 1)
   p += 32
  rare = [i for i, (r_, _) in enumerate(want) if r_]
  bad += listed != rare
  bad += meter - regs["x19"] != 20 * sum(exits) + 25 * len(rare)
  return cnt, bad, len(rare), exits

 def plain(n):
  ws = []
  while len(ws) < n:
   x, y = rr.randrange(1, q), rr.randrange(1, q)
   if not model(x, y, 0)[0]:
    ws.append((x, y))
  return ws

 out = dict(walks=0, mismatches=0, rare_entries=0, exits=0, walk0_sub_exits=0, product_exits=0,
       plain_rare_entries=0, plain_exits=0, plain_filter_hits=0)
 FNS = (("pass1f", ""), ("pass1p", "plain_"))
 counts = {}
 for n in (1, 129, 130):
  ws = plain(n)
  kc, xt = [rr.getrandbits(64) & ~1 for _ in ws], [rr.randrange(q) for _ in ws]
  for fn, v in FNS:
   c, b, nr, ex = run(ws, kc, xt, fn)
   counts[fn, n] = c
   out["mismatches"] += b + nr + sum(ex)
  out["walks"] += n
 # rare kinds by position: 1 distinguished key, 2 near-tie x, 3 odd y with y0 > q0 + 1, 4 filter hit, 5 fe_sub exit
 for n in (2, 3, 7, 130):
  ws, kc, xt = plain(n), [], []
  for i, (x, y) in enumerate(ws):
   kind = (i - n) % 6
   if kind == 1:
    K = (rr.getrandbits(256) & ~(M64 << 64)) % q & ~1 or 2
    ws[i] = (x, K if rr.random() < 0.5 else q - K)
   elif kind == 2:
    ws[i] = (stress_x(rr), y)
   elif kind == 3:
    ws[i] = (x, ((rr.randrange(q >> 64) << 64) | q0 + 2 + 2 * rr.randrange(1, C0 // 2 - 2)) | 1)
   K = ws[i][1] if ws[i][1] % 2 == 0 else q - ws[i][1]
   kc.append((K & M64) + 1 if kind == 4 else rr.getrandbits(64) & ~1)
   if kind == 5:
    a = rr.randrange(q - (1 << 106))
    xt.append(a)
    ws[i] = (a + (rr.randrange(1, 1 << 40) << 64) - rr.randrange(C0), ws[i][1])
   else:
    xt.append(rr.randrange(q))
  for fn, v in FNS:
   c, b, nr, ex = run(ws, kc, xt, fn)
   out["mismatches"] += b
   out[v + "rare_entries"] += nr
   out[v + "exits"] += sum(ex)
   out["walk0_sub_exits"] += ex[0] and not v
  out["plain_filter_hits"] += sum((i - n) % 6 == 4 for i in range(n))
  out["walks"] += n
 # forced product carries: walk i - 1 makes the running product exactly 2^128, walk i then carries out of limb 1
 for n, at in ((3, (2,)), (12, (2, 6, 11)), (131, (2, 64, 128, 130))):
  ws, xt, F = plain(n), [], None
  for i, (x, y) in enumerate(ws):
   if i + 1 in at:
    d = 2 ** 128 * pow(F, -1, q) % q
   elif i in at:
    d = (2 ** 128 - 1) + ((rr.randrange(2, 1 << 20) << 128) // C0 << 128)
   else:
    d = rr.randrange(1, q)
   xt.append((x + d) % q)
   dx = sub(xt[-1], x)[0]
   F = dx if F is None else red(dx * F)[0]
  kc = [rr.getrandbits(64) & ~1 for _ in ws]
  for fn, v in FNS:
   c, b, nr, ex = run(ws, kc, xt, fn)
   out["mismatches"] += b + any(not ex[i] for i in at)
   out[v + "exits"] += sum(ex)
  out["walks"] += n
  out["product_exits"] += len(at)
 for fn, v in FNS:
  out.update({v + "fixed": counts[fn, 1], v + "group_per_walk": (counts[fn, 129] - counts[fn, 1]) / 128,
        v + "tail_walk": counts[fn, 130] - counts[fn, 129]})
 return out

def asm_unit_check(text, nvec):
 """Run the emitter wrappers of kernels/secp_tests.s on edge, random and carrying inputs."""
 q = PROD["q"]
 rr = random.Random(1357)
 prog, labels = a64_load(text)
 RP, AP, BP, CP, SB = 0x1000, 0x2000, 0x3000, 0x4000, 0x90000
 edge = _edge_values()
 canon_edge = [v for v in edge if v < q]

 def call(name, a=None, b=None, z=None, c=None):
  mem = {}
  if c is not None:
   for k, v in enumerate(_limbs(c)):
    mem[CP + 8 * k] = v
  if z is not None:
   for k, v in enumerate(_limbs(z % 2**256) + _limbs(z >> 256)):
    mem[AP + 8 * k] = v
  else:
   for k, v in enumerate(_limbs(a)):
    mem[AP + 8 * k] = v
   if b is not None:
    for k, v in enumerate(_limbs(b)):
     mem[BP + 8 * k] = v
  a64_run(prog, labels, name, {"x0": RP, "x1": AP, "x2": BP, "x3": CP, "sp": SB}, mem)
  return sum(mem[RP + 8 * k] << (64 * k) for k in range(4))

 res = dict(checked={}, mismatches=0)

 def tally(name, ok):
  res["checked"][name] = res["checked"].get(name, 0) + 1
  res["mismatches"] += not ok

 mul_pairs = [(a, b) for a in edge[::2] for b in edge[::3]] + [(rr.getrandbits(256), rr.getrandbits(256)) for _ in range(nvec)]
 for a, b in mul_pairs:
  r = call("t_mul", a, b)
  tally("t_mul", _red_ok(a * b, r))
 for a in edge + [rr.getrandbits(256) for _ in range(nvec)]:
  r = call("t_sqr", a)
  tally("t_sqr", _red_ok(a * a, r))
 sub_pairs = [(a, b) for a in edge[::2] for b in canon_edge[::2]] + [(rr.getrandbits(256), rr.randrange(q)) for _ in range(nvec)]
 for a, b in sub_pairs:
  r = call("t_sub", a, b)
  tally("t_sub", r % q == (a - b) % q and (r < q or (a >= q and b < C0)))
 trip = [(a, b, c) for a in edge[::4] for b in canon_edge[::3] for c in canon_edge[::4]] + _negative_inputs(rr, 40)
 for a, b, c in trip + [(rr.getrandbits(256), rr.randrange(q), rr.randrange(q)) for _ in range(nvec // 2)]:
  tally("t_sqr_sub2", _fused_ok(a * a, [b, c], call("t_sqr_sub2", a, b, c=c)))
  tally("t_mul_sub", _fused_ok(a * b, [c], call("t_mul_sub", a, b, c=c)))
 zs = _second_carry_inputs(rr, nvec // 2) + [rr.getrandbits(512) for _ in range(nvec // 2)] + [2**512 - 1, 0]
 for z in zs:
  r = call("t_reduce", z=z)
  tally("t_reduce", _red_ok(z, r))
 res["reduce_final_carry_inputs"] = nvec // 2
 return res


PRICE_SHARE = Fraction(25616, 2 ** 20)
PRICE_STEP = Fraction(35811, 64) + PRICE_SHARE
PRICE_TRAVERSAL, PRICE_DOUBLING = PRICE_STEP + 39, PRICE_STEP + 677
UNROLL = dict(pass1=128, pass2=128)
ROWS = dict(F1_truncated_product=13, canonical_index_and_slot=10, parity_branch_K_distinguished=3,
   bfi_address=1, arrival_traffic=1.5, checkpoint_filter_and_rare_branch=2.5, arrival_loop_tbz=1 / 64)
BATCH_EXTRAS = 12 + 3 + 7 + 3 + 2 + 4


def kernels_check():
 import os
 from decimal import Decimal as D, getcontext
 getcontext().prec = 60
 here = os.path.dirname(os.path.abspath(__file__))
 res = {}
 res["python_port_vs_bigint"] = _kernel_port_check(2500)
 res["python_batch_port_vs_affine"] = _batch_port_check(24)
 kdir = next((b for b in (os.path.join(here, "kernels"), "/snapshot/files/kernels", os.path.join(here, "files", "kernels"))
    if os.path.exists(os.path.join(b, "a64.py"))), None)
 if kdir is None:
  res["asm"] = "kernels/ not found"
  return res
 sys.dont_write_bytecode = True
 sys.path.insert(0, kdir)
 import a64  # the interpreter and counter (kernels/a64.py)
 globals().update(a64_load=a64.a64_load, a64_run=a64.a64_run)
 texts = {}
 for fname in ("secp_loops.s", "secp_fe.s", "secp_tests.s", "arrive.s", "secp_fuse.s"):
  with open(os.path.join(kdir, fname), encoding="utf-8") as fh:
   texts[fname] = fh.read()
 asm = {}
 for fname in ("secp_loops.s", "secp_fe.s"):  # the wrappers and the arrival rows are not loops
  asm.update(a64._asm_counts(texts[fname]))
 res["asm_executed_vs_affine"] = asm_exec_check(texts["secp_loops.s"])
 res["arrival_rows_executed"] = arr = arrival_check(texts["arrive.s"])
 res["fused_pass1_executed"] = fz = fused_check(texts["secp_fuse.s"])
 res["rare_exits_executed"] = rex = a64.exit_check(texts["secp_loops.s"], dict(
  p_sub=p_sub, p_mul=p_mul, p_sqr=p_sqr, p_rows=p_rows, p_rmul=p_rmul, limbs=_limbs, val=_val, C0=C0, q=PROD["q"], seed=seed))
 try:
  res["asm_emitters_vs_bigint"] = asm_unit_check(texts["secp_tests.s"], 400)
 except Exception as exc:  # keep the rest of the evidence if the wrapper run fails
  res["asm_emitters_vs_bigint"] = "error: %r" % (exc,)
 res["asm_per_function"] = {k: dict(instructions=v["instructions"], loop_iteration_max=v["loop_iteration_max"],
         inner_branches=[l["inner_branches"] for l in v["loops"]], calls=v["calls"])
       for k, v in sorted(asm.items())}
 loops_single_block = all(l["inner_branches"] == 0 for k in ("pass1", "pass2") for l in asm[k]["loops"]) \
  and len(asm["pass1"]["loops"]) == 2 and len(asm["pass2"]["loops"]) == 2  # tail loop, main loop
 per_walk = {k: asm[k]["loop_iteration_max"] / UNROLL[k] for k in UNROLL}
 rows_plain = Fraction(sum(ROWS.values()))
 fused = Fraction(int(fz["group_per_walk"] * 128), 128)
 arith = Fraction(int(fz["plain_group_per_walk"] * 128), 128) + Fraction(per_walk["pass2"])  # plain rounds
 mulf, sqrf = asm["fe_mul_f"]["instructions"], asm["fe_sqr_f"]["instructions"]
 glue = asm["fe_inv"]["instructions"] - asm["fe_inv"]["calls"]
 inv_bound = 255 * (sqrf + 1) + 15 * (mulf + 1) + glue + 255 * 4
 fixed_per_batch = (fz["fixed"] + asm["pass2"]["instructions"] - sum(l["iteration"] for l in asm["pass2"]["loops"])
     + 64 + BATCH_EXTRAS)
 batch_overhead = inv_bound + fixed_per_batch
 # n = 2^20 leaves 127 + 127 tail walks every round
 tails_round = 127 * (fz["tail_walk"] - fused) + 127 * (asm["pass2"]["loops"][0]["iteration"] - Fraction(per_walk["pass2"]))
 round_charged = math.ceil(batch_overhead + tails_round)
 step_exact = arith + Fraction(round_charged, 2 ** 20)
 step = float(step_exact)
 res.update(
  loops_are_single_block=loops_single_block,
  fused_pass1_per_walk=float(fused), unfused_arrival_rows=ROWS, unfused_rows_total=float(rows_plain),
  step_instructions_bound=step,
  fe_mul_f_instructions=mulf, fe_sqr_f_instructions=sqrf,
  inversion_instructions_bound=inv_bound, per_batch_fixed_instructions=batch_overhead,
  per_batch_overhead_below_2_pow_15=batch_overhead < 2**15,
  priced=dict(step=float(PRICE_STEP), traversal=float(PRICE_TRAVERSAL), doubling=float(PRICE_DOUBLING)),
  tail_walks_per_round=float(tails_round), per_round_work_charged=round_charged,
  step_price_exact="%s + %d/2^20" % (arith, round_charged),
  traversal_instructions_bound=step + 24 + 2 + 7 + 6,
  loop_iteration_per_walk=per_walk,
  asm_executed_matches_counts=bool(res["asm_executed_vs_affine"]["mismatches"] == 0 and all(set(res["asm_executed_vs_affine"]["per_extra_walk_divided"][k]) == {per_walk[k]} for k in UNROLL)),
  price_covers_bound=bool(step_exact == PRICE_STEP and step_exact + 24 + 2 + 7 + 6 <= PRICE_TRAVERSAL
        and batch_overhead == 24796 and round_charged == 25616 and batch_overhead < 2**15),
  margin_step_instructions=float(PRICE_STEP - step_exact),
  row_port_paths=all(KST.get(k) for k in ("row_shift", "row_shiftw", "row_fold", "row_foldw")),
  rare_exits_match=bool(rex["mismatches"] == 0 and rex["stubs"] == 1554 and rex["stubs_entered"] == 1554
       and rex["exits"] == rex["stub_entries"]
       and res["asm_executed_vs_affine"].get("exits_taken", 0) == 0),
  arrival_rows_match=bool(arr["mismatches_vs_model"] == 0 and arr["mismatches_vs_exact"] == 0
        and arr["per_walk"] == [float(rows_plain)]),
  fused_pass1_match=bool(fz["mismatches"] == 0 and fz["group_per_walk"] == 129.03125 and fz["tail_walk"] == 133
        and fz["fixed"] == 80 and fz["rare_entries"] > 0 and fz["exits"] > 0 and fz["walk0_sub_exits"] > 0),
  plain_pass1_match=bool(fz["mismatches"] == 0 and fz["plain_group_per_walk"] == 127.03125 and fz["plain_tail_walk"] == 131
        and fz["plain_fixed"] == 78 and fz["plain_filter_hits"] > 0
        and fz["plain_rare_entries"] == fz["rare_entries"] - fz["plain_filter_hits"] and fz["plain_exits"] == fz["exits"]))

 # operation costs (proof.md 2) and the identities that tie them to the loop counts
 OPC = dict(mul=32 + 30 + 23, sqr=20 + 8 + 6 + 7 + 23, mul_sub=62 + 24 + 5, sqr_sub2=41 + 24 + 10, sub=9, canon=8,
   add=9, add_weak_twice=14, mul_add=62 + 23 + 5, sqr_sub2_bias=41 + 23 + 10, rmul=32 + 34 + 6, shift=4)
 # the loops' kernels with rare exits (proof.md 2.5): fe_sub 7, plain fe_mul 84, y3' = lambda*u + y_W 89, x3 73
 mul_r, sub_r, mul_add_r, x3_r = OPC["mul"] - 1, OPC["sub"] - 2, OPC["mul_add"] - 1, OPC["sqr_sub2_bias"] - 1
 walk1 = mul_r + sub_r + 4 + 4  # 99
 walk2 = 2 * OPC["rmul"] + 3 * OPC["shift"] + mul_r + x3_r + mul_add_r + 2 * sub_r + 12 + 4  # 432
 ident1, ident2 = 128 * walk1 + 64 + 2, 128 * walk2 + 64 + 2  # + pointer LDPs, SUBS + B.NE
 tails = 127 * (1 + walk1 + 2 - ident1 / 128) + 127 * (1 + walk2 + 2 - ident2 / 128)
 # doubling (proof.md section 4): its arithmetic with the loops' rare-exit fe_mul, fe_sqr and fe_sub, and <= 30
 # loads/stores/pointers/loop control
 dbl_arith = 4 * mul_r + OPC["mul_sub"] + (OPC["sqr"] - 1) + OPC["sqr_sub2"] + OPC["add"] + 2 * OPC["add_weak_twice"] + sub_r + 30
 dbl_total = step_exact + dbl_arith + Fraction(round_charged, 2 ** 20) + 34 + 3  # the displaced step, moves, escape
 res["doubling"] = dict(op_costs=OPC, identity_pass1=ident1, identity_pass2=ident2, tail_walks_bound=tails,
      identities_hold=bool(ident1 == asm["pass1"]["loop_iteration_max"] and ident2 == asm["pass2"]["loop_iteration_max"]),
      arithmetic_bound=dbl_arith, total_bound=float(dbl_total), priced=float(PRICE_DOUBLING),
      price_covers_bound=bool(dbl_total <= PRICE_DOUBLING))
 two = D(2)
 lg = lambda x: x.ln() / two.ln()
 r = PROD["r"]
 N = D((r - 1) // 6)
 Rr = D(2 ** 25)
 Np = N * Rr / (Rr - 1)
 kappa, p0 = D(str(KAPPA)), D(str(P0))
 S0 = (2 * Np * (D(1) / (D(1) - p0)).ln()).sqrt()
 Sops = kappa * S0
 dec = lambda f: D(f.numerator) / D(f.denominator)
 step_u, dbl_u, trav_u = dec(PRICE_STEP) / 1000, dec(PRICE_DOUBLING) / 1000, dec(PRICE_TRAVERSAL) / 1000
 ph1 = step_u * (Sops + two ** 60)
 # report store of upstream-proof.md section 3: integer D_cap, bucket count B, probe depth 8,
 # overflow list of B 8-byte addresses, report price 40 + 8*0.05
 from decimal import ROUND_CEILING
 Dcap = int((two ** D("63.83")).to_integral_value(rounding=ROUND_CEILING))
 B = -((-(2 * Dcap)) // 2 ** 16)
 report_price = D(40) + 8 * D("0.05")
 ovf_insert = 6 * B + 35
 ovf_lookup = 16 * B
 claim = D("125.862153")
 Wmax = two ** claim
 refills = int((Wmax / two ** 40).to_integral_value(rounding=ROUND_CEILING)) + 2 ** 20
 Lmax = 20 * two ** 64
 capped = dict(
  meter_remainders=two ** 61,
  meter_block_refills=D(refills) * 30 / 1000,
  phase2=two ** 20 * (Lmax + 65) * 600,
  rewalks=8 * 2 * ((Lmax + 65) * 600 + two ** 13),
  reports=D(Dcap) * report_price,
  overflow_worst_case=D(Dcap) * (ovf_insert + ovf_lookup) / 1000,
  store_init=D(Dcap) * 2 * 16 / 1000,
  table_verification=two ** 39,
  node_setup=two ** 26,
  checkpoint_writes=(Sops + two ** 61) * 2 / PERIOD / 1000,  # 2 stores per anchor
  anchor_rounds=(Sops + two ** 61) * 14 / PERIOD / 1000,  # full K 7 + counter 5 + unfused pass 1 1.5, rounded
  compare_rounds=(Sops + two ** 61) * 2 / STRIDE / 1000,  # LDR of the K_cp word + CCMP per walk
 )
 hard = sum(capped.values())
 rate = two ** D("-27.5849")  # upper bound on 1/(6R) + 2^-47.999 = 2^-27.58496 (R = 2^25)
 random_items = dict(
  doubling=rate * (Sops + two ** 60) * (dbl_u - step_u),
  traversal=(trav_u - step_u) * two ** D("-26.5849") * (Sops + two ** 60),  # >= 2/(6R) + 32*2^-47.999 (TRAV_CAP)
  starts=two ** 13 * (Sops / two ** 64 + two ** 20),
  other=two ** 70,
  rare_entries=D("0.025") * (2 * rate + two ** -56 + D(C0) / two ** 64) * (Sops + two ** 60),  # hits (2 per doubling), windows, distinguished, odd-y borrows
  rare_exits=D("0.02") * (14 * D(C0) / two ** 64 + 4 * two ** -63) * (Sops + two ** 61 + rate * (Sops + two ** 60)),
 )
 extras = sum(random_items.values())
 head = Wmax - ph1 - hard
 markov = extras / head
 allowance = (extras / D("0.000065")).to_integral_value(rounding=ROUND_CEILING)
 w_ledger = ph1 + hard + allowance
 success = p0 - D("0.000065") - D("1.1e-7")
 res["ledger"] = dict(
  log2_S0=float(lg(S0)), log2_Sops=float(lg(Sops)), log2_phase1=float(lg(ph1)),
  log2_hard_extras={k: round(float(lg(v)), 2) for k, v in capped.items()},
  log2_hard_total=float(lg(hard)),
  log2_random_items={k: round(float(lg(v)), 2) for k, v in random_items.items()},
  log2_expected_random_total=float(lg(extras)),
  log2_headroom=float(lg(head)), markov_early_halt_probability=float(markov),
  success_lower_bound=float(success),
  success_with_unrounded_ratio=float(p0 - markov - D("1.1e-7")),
  claimed_time_log2=float(claim),
  store=dict(buckets_per_node=B, bucket_bytes=104, overflow_address_bytes=8,
    probe_instructions=39, probe_price=float(D("0.05")),
    report_price=float(report_price), overflow_insert=ovf_insert,
    overflow_lookup=ovf_lookup, blocks_per_refill=2 ** 40, refills=refills),
  log2_allowance=float(lg(allowance)), log2_w_ledger=float(lg(w_ledger)),
  w_ledger_below_claim=bool(w_ledger <= Wmax),
  claim_is_six_decimal_ceiling=bool(w_ledger <= Wmax and lg(w_ledger) > claim - D("0.000001")),
  claim_above_phase1_plus_hard=bool(Wmax > ph1 + hard),
  success_at_least_0_39=bool(success > D("0.39")),
  success_at_least_0_3900001=bool(success > D("0.3900001")),
  model_probability_at_S_ops=float(1 - (1 - p0) ** (kappa * kappa)),
  hard_below_2_pow_111_1493=bool(hard < two ** D("111.1493")),
  random_below_2_pow_98_8272=bool(extras < two ** D("98.8272")),
  kappa=float(kappa), p0=float(p0),
  markov_below_0_000065=bool(markov < D("0.000065")),
 )
 return res



def done(out):
 out["seconds"] = round(time.time() - T0, 2)
 print(json.dumps(out, indent=1))
 s = json.dumps(dict(out, config=0))
 if "false" in s:
  sys.exit("false in the result: %s" % re.findall(r'"([^"]*)": false', s))


def main():
 out = {"config": cfg}
 if kernels:
  out["kernels"] = kernels_check()
 if prod_checks:
  out["part_a"] = part_a(prod_checks)
 if long_trails:
  out["long_trail_test"] = dict(curve=curve, R=R, long_len=long_len, escape="doubling",
         **long_trail_test(make_toy(curve), rng, R, long_len, long_trails))
 if run_toy:
  out["rule_check"] = rule_check(make_toy(curve), random.Random(seed + 61), 8, dp_bits, 600)
 if not run_toy:
  return done(out)
 E = make_toy(curve)
 r = E.r
 N = (r - 1) // CLS
 Np = N * R / (R - 1)
 S0 = math.sqrt(2 * Np * math.log(1 / (1 - P0)))  # pure model: probability P0
 S1 = KAPPA * S0 + (1 << dp_bits) - 1  # the budget of H1 plus the toy's own-trail delay (proof.md 5.2)
 out["toy"] = dict(q=E.q, r=r, h=E.h, beta=E.beta, lam=E.lam, N=N, R=R, dp_bits=dp_bits,
     mean_trail_over_sqrtN=round((1 << dp_bits) / math.sqrt(N), 4),
     S0=round(S0, 1), S1=round(S1, 1),
     model_mean_tc=round(math.sqrt(math.pi * Np / 2), 1))
 res = []
 for _ in range(trials):
  tr = toy_trial(E, rng, R, dp_bits)
  if tr is None:
   break
  res.append(tr)
 n = len(res)
 if n:
  tcs = [x["tc"] for x in res if x["tc"] is not None]
  steps = sum(x["steps"] for x in res)
  cens = pipeline  # in pipeline mode T_c is censored at the budget
  summ = dict(
   trials_completed=n,
   mode="pipeline (scaled budget, phase 2, report cap, 8-match cap)" if pipeline else "run until solved",
   pipeline_success_rate=round(sum(1 for x in res if x["ok"]) / n, 4) if pipeline else None,
   all_solved_correctly=None if pipeline else all(x["ok"] for x in res),
   mean_tc_over_sqrtN=None if cens else round(sum(tcs) / len(tcs) / math.sqrt(N), 4),
   stderr_mean_tc_over_sqrtN=None if cens else round(math.sqrt(sum((t / math.sqrt(N) - sum(tcs) / len(tcs) / math.sqrt(N)) ** 2 for t in tcs) / max(1, len(tcs) - 1) / len(tcs)), 4),
   model_mean_over_sqrtN=round(math.sqrt(math.pi * Np / 2) / math.sqrt(N), 4),
   frac_tc_le_S0=round(sum(t <= S0 for t in tcs) / n, 4),
   frac_tc_le_S1=round(sum(t <= S1 for t in tcs) / n, 4),
   binomial_stderr_at_0p39=round(math.sqrt(0.39 * 0.61 / n), 4),
   model_P_S0=P0,
   H1_floor_P_S1=P0,
   model_P_S1=round(1 - math.exp(-S1 * S1 / (2 * Np)), 4),
   mean_solve_over_tc=None if cens else round(sum(x["solve_visits"] / x["tc"] for x in res if x["tc"]) / len(tcs), 3),
   visits_per_step=round(sum(x["solve_visits"] for x in res) / max(1, sum(x["steps"] + x["doublings"] for x in res)), 4),
   total_steps=steps,
   doublings_per_step=sum(x["doublings"] for x in res) / max(1, steps),
   cycle_model_1_over_6R=1 / (6 * R),
   checkpoint_escapes_per_step=sum(x["esccp"] for x in res) / max(1, steps),
   abandoned_trails=sum(x["abandoned"] for x in res),
   abandoned_visit_fraction=round(sum(x.get("abandoned_visits", 0) for x in res) / max(1, sum(x["solve_visits"] for x in res)), 4),
   exceptional_additions=sum(x["exceptional"] for x in res),
   node_overflows=sum(x.get("node_overflow", 0) for x in res),
   max_node_load_over_cap=max((x["max_node_load_over_cap"] or 0) for x in res) if pipeline else None,
   mean_max_node_load_over_mean=(sum(x["max_node_load_over_mean"] or 0 for x in res) / n) if pipeline else None,
  )
  out["toy_summary"] = summ
  if dump:
   out["tc_list"] = tcs
 done(out)


main()
