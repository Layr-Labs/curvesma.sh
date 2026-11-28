"""q15_ecm16.py — Deterministic Stepped-GCD + Stage-2 Binary-Product-Tree Positive-Rank Z/2 x Z/8
Torsion-16 Edwards (a = +1) ECM (BN254 v40 parametrization on E_W: y^2 = x^3 + 2x^2 - 8x, G_0 = (9/4, 15/8),
Montgomery x-only 11.0 muls/bit, explicit E_W(Z/NZ) i.i.d. uniform k_i ~ Unif([1, 2^256]) pricing +3,800*mul + 4*gcd, exact L_min=2 branch-tree
DP, W <= B1 single-leaf-unique BSGS wheels (30, 210, 2310, 30030, 510510), recursive Stage-2 product-tree
splitting, test-once Newton + MR_ROUNDS=110 primality budget, continuous 0.05b real-interval Theorem T2)."""
from fractions import Fraction as F
from functools import lru_cache
import math
import operator as _OP
import random
from math import ceil, isqrt, log, log1p, log2, sqrt
import q15_descent16 as D

_SUMPROD = getattr(math, "sumprod", lambda a, b, _m=_OP.mul: sum(map(_m, a, b)))

MONT5 = D.MONT5
TORSION = 16
_LOG2_TORSION = log2(TORSION)  # exact 4.0 (Z/2 x Z/8 subset E_d(Q) => 16 | #E_d(F_p) for all good p)
STEPPER_MULS = 3800.0  # Explicit E_W(Z/NZ) fresh i.i.d. uniform random k_i ~ Unif([1, 2^256]) scalar-mul budget (projective v=x/2 10M/bit ladder + Okeya-Sakurai + ew_to_mont + 4 Mont conv <= 2596 muls; honest CIOS + spills + 5L+2 adds + cswap <= 0.9599 * STEPPER_MULS * mul at L=2, <= 0.9109 at L=3, <= 0.669 at L>=5)
STEPPER_GCDS = 4.0   # Explicit E_W(Z/NZ) curve-mapping budget (2 extended Binary Stein GCDs: 1 in Okeya-Sakurai recovery + 1 in ew_to_mont_mod <= 4 plain GCDs)
WHEELS = (30, 210, 2310, 30030, 510510)
B_MIN_DEFAULT = 32   # Conservative 2 * 16 bits (all primes <= 2^16 stripped by TDIV)
L_MIN_DEFAULT = 2    # 2-limb (<=128b) mont2_cios kernel (47 <= 48 ARMv8 ops <= 58.84 Ir charged)
NEWTON_PRECHECK_INSTR = 2.0 ** 26  # Honest schoolbook upper bound (> 32.06M Ir for 2590b input)
MR_ROUNDS = 110      # Unconditional Rabin error <= 4^-110 = 2^-220 per test (< 2^-96 << 2^-80 over < 2^124 tests in T_run=2^124 u)


def primes_upto(n):
 a = bytearray(b"\x01") * (n + 1)
 a[0:2] = b"\x00\x00"
 for p in range(2, isqrt(n) + 1):
  if a[p]:
   a[p * p :: p] = b"\x00" * ((n - p * p) // p + 1)
 return [p for p in range(2, n + 1) if a[p]]


def _is_perfect_power(n):
 """Exact prime-exponent Newton integer perfect-power precheck before ECM (eliminates N = p^k;
    explicitly charged at NEWTON_PRECHECK_INSTR = 2^26 Ir)."""
 if n < 4:
  return None
 r2 = isqrt(n)
 if r2 * r2 == n:
  return (r2, 2)
 bl = n.bit_length()
 for k in primes_upto(bl):
  if k == 2:
   continue
  x = 1 << ((bl + k - 1) // k)
  km1 = k - 1
  while True:
   nx = (km1 * x + n // (x ** km1)) // k
   if nx >= x:
    break
   x = nx
  if x ** k == n:
   return (x, k)
 return None


# --- Positive-rank Z/2 x Z/8 Torsion-16 Edwards (a = +1) curve family on E_W: y^2 = x^3 + 2x^2 - 8x ---

def ew_add_q(P, Q):
 """Exact rational chord-and-tangent addition on E_W(Q): y^2 = x^3 + 2x^2 - 8x."""
 x1, y1 = P
 x2, y2 = Q
 lm = (3 * x1 * x1 + 4 * x1 - 8) / (2 * y1) if x1 == x2 else (y2 - y1) / (x2 - x1)
 x3 = lm * lm - 2 - x1 - x2
 return (x3, lm * (x1 - x3) - y1)


def ed_add_q(A, B, d):
 """Exact rational Edwards (a = +1) addition on x^2 + y^2 = 1 + d*x^2*y^2."""
 x1, y1 = A
 x2, y2 = B
 t = d * x1 * x2 * y1 * y2
 return ((x1 * y2 + y1 * x2) / (1 + t), (y1 * y2 - x1 * x2) / (1 - t))


def ew_g0_mod(N):
 """Base point G_0 = (9/4, 15/8) on E_W(Z/NZ) for odd N (exact shifts, zero GCDs)."""
 inv2 = (N + 1) >> 1
 inv4 = (inv2 * inv2) % N
 inv8 = (inv4 * inv2) % N
 return (9 * inv4) % N, (15 * inv8) % N


def ew_mod(P, Q, N):
 """Modular chord-and-tangent addition on E_W(Z/NZ): y^2 = x^3 + 2x^2 - 8x (<= 4 muls + 1 ext-GCD)."""
 if P is None:
  return Q
 if Q is None:
  return P
 x1, y1 = P
 x2, y2 = Q
 if x1 == x2:
  gy = math.gcd((y1 + y2) % N, N)
  if 1 < gy < N:
   return ("fac", gy)
  if gy == N:
   return None
  den = (2 * y1) % N
  num = (3 * x1 * x1 + 4 * x1 - 8) % N
 else:
  den = (x2 - x1) % N
  num = (y2 - y1) % N
 gd = math.gcd(den, N)
 if 1 < gd < N:
  return ("fac", gd)
 if gd == N:
  return None
 lm = (num * pow(den, -1, N)) % N
 x3 = (lm * lm - 2 - x1 - x2) % N
 return (x3, (lm * (x1 - x3) - y1) % N)


def ew_mul_aff_ref(k, G0, N):
 """Reference affine double-and-add scalar multiple [k]G_0 on E_W(Z/NZ) (used only to cross-check ew_mul_proj_counted)."""
 res = None
 base = G0
 e = k
 while e > 0:
  if e & 1:
   res = ew_mod(res, base, N)
   if isinstance(res, tuple) and len(res) == 2 and res[0] == "fac":
    return res
  e >>= 1
  if e:
   base = ew_mod(base, base, N)
   if isinstance(base, tuple) and len(base) == 2 and base[0] == "fac":
    return base
 return res


def ew_mul_proj_counted(k, N):
 """Projective v=x/2 10M/bit Montgomery ladder + Okeya-Sakurai (x,y) recovery on E_W: y^2=x^3+2x^2-8x, G0=(9/4,15/8).
 Exact match vs affine ew_mul_aff_ref for all k != 0, -1 (mod ord_p(G_0)) and gcd(60, N) == 1 (always true post-TDIV);
 captures any composite split whenever k == 0 or -1 (mod ord_p(G_0)) for a proper factor p | N."""
 if k <= 0:
  return None, 0, 0, 0
 if k == 1:
  return ew_g0_mod(N), 2, 1, 0
 muls = adds = 0
 def _xdbl(V, Z):
  nonlocal muls, adds
  V2 = (V * V) % N; Z2 = (Z * Z) % N; S = (V + Z) % N; twoVZ = ((S * S) % N - V2 - Z2) % N; t2 = (Z2 + Z2) % N
  V3 = pow((V2 + t2) % N, 2, N); Z3 = (twoVZ * ((V2 + V2 + twoVZ - t2 - t2) % N)) % N
  muls += 5; adds += 9
  return V3, Z3
 def _xadd(V1, Z1, V2, Z2):
  nonlocal muls, adds
  A = (V1 * V2) % N; B = (Z1 * Z2) % N; E = ((V1 + Z1) * (V2 - Z2)) % N
  s2 = pow((A + B + B) % N, 2, N); d2 = pow((A - B - E) % N, 2, N)
  s2_2 = (s2 + s2) % N; s2_4 = (s2_2 + s2_2) % N; V3 = (s2_4 + s2_4) % N
  d4 = ((d2 + d2) % N) * 2 % N; Z3 = ((d4 + d4) % N + d2) % N
  muls += 5; adds += 13
  return V3, Z3
 V0, Z0 = 9 % N, 8 % N
 V1, Z1 = _xdbl(V0, Z0)
 for bit in bin(k)[3:]:
  if bit == "1":
   V0, Z0 = _xadd(V0, Z0, V1, Z1); V1, Z1 = _xdbl(V1, Z1)
  else:
   V1, Z1 = _xadd(V0, Z0, V1, Z1); V0, Z0 = _xdbl(V0, Z0)
 X1 = (V0 + V0) % N; X2 = (V1 + V1) % N; X1_4 = (4 * X1) % N; Z0_9 = (9 * Z0) % N
 p1 = (((9 * X1 - 32 * Z0) % N) * ((Z0_9 + X1_4) % N)) % N; p2 = pow((Z0_9 - X1_4) % N, 2, N)
 X1Z0 = (X1 * Z0) % N; Z0Z1 = (Z0 * Z1) % N; t60 = (60 * Z0Z1) % N; den = (Z0 * t60) % N
 num = (((p1 + 144 * X1Z0) % N) * Z1 - X2 * p2) % N
 muls += 7; adds += 32
 g = math.gcd(den, N)
 if 1 < g < N:
  return ("fac", g), muls, adds, 1
 if g == N or g == 0:
  return None, muls, adds, 1
 inv_d = pow(den, -1, N); muls += 3
 return ((X1 * ((t60 * inv_d) % N)) % N, (num * inv_d) % N), muls, adds, 1


def ew_mul_mod(k, G0, N):
 """Scalar multiple [k]G_0 on E_W(Z/NZ) via projective v=x/2 10M/bit Montgomery ladder + Okeya-Sakurai recovery."""
 res, _, _, _ = ew_mul_proj_counted(k, N)
 return res


def ew_to_mont_mod(Qk, sg, N):
 """Map Qk = (x, u) in E_W(Z/NZ) and sign sg in {+1, -1} to Montgomery parameters
    (a24 = 1/(1-d), X0 = dn^2 + x8^2, Z0 = dn^2 - x8^2) of the Z/2 x Z/8 torsion-16 Edwards (a=+1)
    curve E_d(Z/NZ) with base point u_0 = (1 + y_0)/(1 - y_0) = X0/Z0 (17 muls + 1 ext-GCD).
    Catches d in {0, 1} (singular A^2 = 4), y_0 in {+1, -1} (2-torsion base point), and E_W poles
    in a single modular GCD check."""
 x, u = Qk
 dx = (x * x - 4) % N
 vn = (2 * (x - 2 + sg * u)) % N
 dx2 = (dx * dx) % N
 vn2 = (vn * vn) % N
 vndx = (vn * dx) % N
 dn = (vn2 - 2 * vndx - dx2) % N
 x8n = (vn2 + dx2) % N
 wn = (vn2 + 2 * vndx - dx2) % N
 x82 = (x8n * x8n) % N
 x84 = (x82 * x82) % N
 dn2 = (dn * dn) % N
 wn2 = (wn * wn) % N
 den_1md = (x84 - wn2 * dn2) % N
 X0 = (dn2 + x82) % N
 Z0 = (dn2 - x82) % N
 p_rest = (((((dx * dn) % N) * x8n % N) * wn % N) * ((X0 * Z0) % N)) % N
 chk = (p_rest * den_1md) % N
 g = math.gcd(chk, N)
 if 1 < g < N:
  return ("fac", g)
 if g == N or g == 0:
  return ("sing", N)
 a24 = (x84 * ((pow(chk, -1, N) * p_rest) % N)) % N
 return ("ok", (a24, X0, Z0))


def _xmul_mont(k, X, Z, a24, M):
 """Montgomery x-only ladder (5M per xdbl + 6M per xadd = 11.0 M/bit)."""
 if k == 1:
  return X, Z

 def _xd(X0, Z0):
  a = (X0 + Z0) % M
  b = (X0 - Z0) % M
  a2 = (a * a) % M
  b2 = (b * b) % M
  c = (a2 - b2) % M
  return (a2 * b2) % M, (c * (b2 + a24 * c)) % M

 def _xa(X1, Z1, X2, Z2, Xd, Zd):
  u1 = ((X1 - Z1) * (X2 + Z2)) % M
  v1 = ((X1 + Z1) * (X2 - Z2)) % M
  return (Zd * pow(u1 + v1, 2, M)) % M, (Xd * pow(u1 - v1, 2, M)) % M

 X0, Z0 = X, Z
 X1, Z1 = _xd(X, Z)
 for bit in bin(k)[3:]:
  if bit == "1":
   X0, Z0 = _xa(X0, Z0, X1, Z1, X, Z)
   X1, Z1 = _xd(X1, Z1)
  else:
   X1, Z1 = _xa(X0, Z0, X1, Z1, X, Z)
   X0, Z0 = _xd(X0, Z0)
 return X0, Z0


def _run_t16_stepped(N, Qk, sg, b1, primes_b1):
 """Run one Torsion-16 Edwards (a = +1) Stage-1 curve mod N taking gcd after EVERY prime power p^k <= B1
    (deterministic stepped-gcd schedule; zero curve replay)."""
 st, par = ew_to_mont_mod(Qk, sg, N)
 if st == "fac":
  return par
 if st != "ok":
  return "SINGULAR"
 a24, x, z = par
 for p in primes_b1:
  pk = p
  while pk <= b1:
   x, z = _xmul_mont(p, x, z, a24, N)
   g = math.gcd(z, N)
   if g > 1:
    return g
   pk *= p
 return 1


def _build_product_tree_mod(raw_leaves, N):
 """Build exact Stage-2 binary product tree mod N over J leaves, replacing any 0 mod N leaf
    with 1 BEFORE multiplying into the tree. Because every remaining leaf has gcd(t_j, N) < N,
    whenever gcd(root, N) == N, binary-tree bisection is unconditionally guaranteed to isolate a
    proper factor 1 < gcd(node, N) < N in <= 1 + 2*ceil(log2(J)) gcds (0 extra muls)."""
 leaves = [t % N for t in raw_leaves]
 leaves = [t if t != 0 else 1 for t in leaves]
 levels = [leaves]
 while len(levels[-1]) > 1:
  cur = levels[-1]
  nxt = [(cur[i] * cur[i + 1]) % N if i + 1 < len(cur) else cur[i] for i in range(0, len(cur), 2)]
  levels.append(nxt)
 return levels


def _split_product_tree(levels, N):
 """Recursive Stage-2 binary-product-tree splitting + test-once prime-power extraction (<= 2*h_tree*Omega(N) GCDs)."""
 g_root = math.gcd(levels[-1][0], N)
 if g_root == 1:
  return [N], 1, 0
 raw_out = []
 if g_root < N:
  raw_out.append(N // g_root)
 desc_gcds = 0

 def _rec(lv, idx, m):
  nonlocal desc_gcds
  if lv == 0:
   raw_out.append(m)
   return
  cur = levels[lv - 1]
  li = 2 * idx
  ri = li + 1
  if ri >= len(cur):
   _rec(lv - 1, li, m)
   return
  gl = math.gcd(cur[li], m)
  desc_gcds += 1
  if gl == m and lv == 1:
   # At leaf level ONLY (zero extra recursion): if left leaf cur[li] == 0 mod m cannot
   # separate m (e.g. a sub-factor collision [pq, p, q, r] or unsanitized [0, p, q, 1]),
   # check right sibling leaf cur[ri] in 1 GCD (<= 2 GCDs at this node, preserving the
   # unconditional deterministic bound desc_gcds <= 2 * h_tree * Omega(g_root) for ALL trees!).
   gr = math.gcd(cur[ri], m)
   desc_gcds += 1
   if 1 < gr < m:
    raw_out.append(gr)
    raw_out.append(m // gr)
   else:
    raw_out.append(m)
   return
  mr = m // gl
  if gl > 1:
   _rec(lv - 1, li, gl)
  if mr > 1:
   gr = math.gcd(cur[ri], mr)
   desc_gcds += 1
   if gr > 1:
    _rec(lv - 1, ri, gr)
   rem = mr // gr
   if rem > 1:
    raw_out.append(rem)

 _rec(len(levels) - 1, 0, g_root)
 out = []
 for m in raw_out:
  if m < N:
   pp = _is_perfect_power(m)
   if pp is not None:
    out.extend([pp[0]] * pp[1])
    continue
  out.append(m)
 return out, 1 + desc_gcds, desc_gcds


def _bisect_product_tree(levels, N):
 """Compatibility wrapper returning one proper factor (if any) and total GCD count via _split_product_tree."""
 parts, total_gcds, _ = _split_product_tree(levels, N)
 for d in parts:
  if 1 < d < N:
   return d, total_gcds
 return 1, total_gcds


def _run_t16_frontier_curve(N, Qk, sg, b1, primes_b1):
 """Run ONE Torsion-16 Edwards (a = +1) Stage-1 curve across the active composite frontier of N:
    whenever 1 < g < d splits an active composite d at curve setup or a prime-power step, immediately
    reduce (x, z, a24) mod (d // g) and CONTINUE THE SAME CURVE on the cofactor d // g (zero extra curve
    count; total mul/gcd across all active branches remains bounded by tree_branch_caps(nbits))."""

 def _init_branch(d):
  if d <= 1:
   return []
  pp = _is_perfect_power(d)
  if pp is not None:
   return _init_branch(pp[0])
  st, par = ew_to_mont_mod((Qk[0] % d, Qk[1] % d), sg, d)
  if st == "fac":
   return _init_branch(par) + _init_branch(d // par)
  if st != "ok":
   return [(d, None, None, None)]
  a24, bx, bz = par
  return [(d, bx, bz, a24)]

 branches = _init_branch(N)
 for p in primes_b1:
  pk = p
  while pk <= b1:
   nxt = []
   for d, bx, bz, a24 in branches:
    if bx is None:
     nxt.append((d, None, None, None))
     continue
    bx, bz = _xmul_mont(p, bx, bz, a24, d)
    g = math.gcd(bz, d)
    if 1 < g < d:
     nxt.append((g, None, None, None))
     cof = d // g
     nxt.append((cof, bx % cof, bz % cof, a24 % cof))
    elif g == d:
     nxt.append((d, None, None, None))
    else:
     nxt.append((d, bx, bz, a24))
   branches = nxt
   pk *= p
 return [d for d, _, _, _ in branches]


def _run_t16_frontier_s1s2(N, Qk, sg, b1, b2, wheel, primes_b1):
 """Run ONE full Torsion-16 Edwards Stage-1 stepped-GCD + Stage-2 binary-product-tree parallel-active-frontier
    curve on N. Returns (parts, total_tree_descent_gcds)."""

 def _init(d):
  if d <= 1:
   return []
  pp = _is_perfect_power(d)
  if pp is not None:
   return _init(pp[0])
  st, par = ew_to_mont_mod((Qk[0] % d, Qk[1] % d), sg, d)
  if st == "fac":
   return _init(par) + _init(d // par)
  if st != "ok":
   return [(d, None, None, None)]
  a24, bx, bz = par
  return [(d, bx, bz, a24)]

 branches = _init(N)
 for p in primes_b1:
  pk = p
  while pk <= b1:
   nxt = []
   for d, bx, bz, a24 in branches:
    if bx is None:
     nxt.append((d, None, None, None))
     continue
    bx, bz = _xmul_mont(p, bx, bz, a24, d)
    g = math.gcd(bz, d)
    if 1 < g < d:
     nxt.append((g, None, None, None))
     cof = d // g
     nxt.append((cof, bx % cof, bz % cof, a24 % cof))
    elif g == d:
     nxt.append((d, None, None, None))
    else:
     nxt.append((d, bx, bz, a24))
   branches = nxt
   pk *= p

 done = [d for d, bx, _, _ in branches if bx is None]
 active = [(d, bx, bz, a24) for d, bx, bz, a24 in branches if bx is not None]
 babies = [r for r in range(1, wheel // 2) if math.gcd(r, wheel) == 1]
 G = ceil(b2 / wheel) + 1
 total_desc_gcds = 0
 for d, bx, bz, a24 in active:
  cur = d
  cur_x, cur_z, cur_a24 = bx, bz, a24
  bpts = [_xmul_mont(r, cur_x, cur_z, cur_a24, cur) for r in babies]
  for g_idx in range(1, G + 1):
   if cur <= 1:
    break
   xg, zg = _xmul_mont(g_idx * wheel, cur_x, cur_z, cur_a24, cur)
   raw = [(xg * zb - xb * zg) % cur for xb, zb in bpts]
   levels = _build_product_tree_mod(raw, cur)
   g_root = math.gcd(levels[-1][0], cur)
   if g_root == 1:
    continue
   parts, _, dg = _split_product_tree(levels, cur)
   total_desc_gcds += dg
   if g_root == cur:
    done.extend(parts)
    cur = 1
    break
   rem = cur // g_root
   pp_r = _is_perfect_power(rem)
   if pp_r is not None:
    rem_base, _ = pp_r
    skipped = 0
    for pt in parts:
     if pt == rem_base and skipped == 0:
      skipped = 1
      continue
     done.append(pt)
    cur = rem_base
   else:
    skipped = 0
    for pt in parts:
     if pt == rem and skipped == 0:
      skipped = 1
      continue
     done.append(pt)
    cur = rem
   cur_x %= cur
   cur_z %= cur
   cur_a24 %= cur
   bpts = [(xb % cur, zb % cur) for xb, zb in bpts]
  if cur > 1:
   done.append(cur)
 return done, total_desc_gcds


def _is_prime_small(n):
 if n < 2:
  return False
 for d in range(2, isqrt(n) + 1):
  if n % d == 0:
   return False
 return True


def verify_t16_edwards_ecm():
 """Live verification of Positive-Rank Z/2xZ/8 Torsion-16 Edwards ECM, projective E_W ladder, stepped-GCD & Stage-2 tree."""
 assert TORSION == 16 and _LOG2_TORSION == 4.0 and STEPPER_MULS == 3800.0 and STEPPER_GCDS == 4.0
 wh_j = {w: stage_ops(1 << 20, 1 << 26, w)[3] for w in WHEELS}
 assert WHEELS == (30, 210, 2310, 30030, 510510)
 assert wh_j[30] == 4 and wh_j[210] == 24 and wh_j[2310] == 240 and wh_j[30030] == 2880 and wh_j[510510] == 46080

 # (1) Nagell-Lutz rank >= 1 on E_W(Q): y^2 = x^3 + 2x^2 - 8x
 G0 = (F(9, 4), F(15, 8))
 nl_ok = (G0[1] * G0[1] == G0[0] ** 3 + 2 * G0[0] * G0[0] - 8 * G0[0] and G0[0].denominator != 1)
 assert nl_ok

 # (2) Algebraic Q proof of Z/2 x Z/8 subset E_d(Q) across 8 distinct curves
 curves = []
 seen = set()
 Q_cur = G0
 tors_q_ok = True
 for _ in range(50):
  x, u = Q_cur
  if x not in (2, -2):
   for sg in (1, -1):
    v = 2 * (x - 2 + sg * u) / (x * x - 4)
    Y = (v * v * x - 2 - 2 * v) / 2
    dn = v * v - 2 * v - 1
    if dn == 0:
     continue
    x8 = (v * v + 1) / dn
    w = (v * v + 2 * v - 1) / dn
    dv = (w * w) / (x8 ** 4)
    if dv in (0, 1) or dv in seen:
     continue
    x0, y0 = Y / dn, x8 * x8
    if x0 * x0 + y0 * y0 != 1 + dv * x0 * x0 * y0 * y0:
     continue
    P8 = (x8, x8)
    P2 = ed_add_q(P8, P8, dv)
    P4 = ed_add_q(P2, P2, dv)
    P8_ = ed_add_q(P4, P4, dv)
    sd = w / (x8 * x8)
    Am = 2 * (1 + dv) / (1 - dv)
    up = -(1 + sd) / (1 - sd)
    if not (
     2 * x8 * x8 == 1 + dv * (x8 ** 4)
     and P2 == (1, 0)
     and P4 == (0, -1)
     and P8_ == (0, 1)
     and up != 0
     and up * up + Am * up + 1 == 0
    ):
     tors_q_ok = False
    seen.add(dv)
    curves.append((dv, x0, y0))
    if len(curves) >= 8:
     break
  if len(curves) >= 8:
   break
  Q_cur = ew_add_q(Q_cur, G0)
 assert tors_q_ok and len(curves) == 8

 # (3) Direct point count: 16 | #E_d(F_p) across all good primes p in (101, 450)
 ps_cnt = [p for p in range(101, 451) if _is_prime_small(p)]
 tot16 = 0
 ok16 = 0
 for dq, x0q, y0q in curves[:4]:
  for p in ps_cnt:
   if dq.denominator % p == 0 or x0q.denominator % p == 0 or y0q.denominator % p == 0:
    continue
   d = (dq.numerator * pow(dq.denominator, -1, p)) % p
   if d in (0, 1):
    continue
   iv = pow((1 - d) % p, -1, p)
   A = (2 * (1 + d) * iv) % p
   B = (4 * iv) % p
   c = 1
   for uu in range(p):
    rh = (B * (uu * uu * uu + A * uu * uu + uu)) % p
    if rh == 0:
     c += 1
    elif pow(rh, (p - 1) // 2, p) == 1:
     c += 2
   tot16 += 1
   if c % 16 == 0:
    ok16 += 1
 assert tot16 > 200 and ok16 == tot16

 # (4) Modular E_W stepper + ew_to_mont_mod + 11.0 muls/bit Montgomery vs affine Edwards
 pt = 1000000007
 G0_23 = ew_g0_mod(23)
 G0_pt = ew_g0_mod(pt)
 step_ok = (
  ew_mod((4, 1), (4, 6), 7) is None
  and ew_mod((60, 71), (60, 50), 77) == ("fac", 11)
  and ew_mod((4, 8), (4, 13), 35) == ("fac", 7)
  and ew_to_mont_mod(G0_23, 1, 23) == ("sing", 23)
 )
 Qm = G0_pt
 Qq = G0
 for _ in range(10):
  if (
   (Qq[0].numerator * pow(Qq[0].denominator, -1, pt)) % pt != Qm[0]
   or (Qq[1].numerator * pow(Qq[1].denominator, -1, pt)) % pt != Qm[1]
  ):
   step_ok = False
  Qq = ew_add_q(Qq, G0)
  Qm = ew_mod(Qm, G0_pt, pt)
 assert step_ok

 # (4b) Live verification of projective v=x/2 10M/bit Montgomery ladder + Okeya-Sakurai (x,y) recovery vs affine ew_mul_aff_ref
 _rng_ew = random.Random(20261128)
 _ew_max_m = _ew_max_a = _ew_max_g = 0
 for _p_ew in (pt, (1 << 127) - 1):
  _G0_ew = ew_g0_mod(_p_ew)
  _ks_ew = [1, 2, 3, 4, 7, 15, (1 << 256) - 1, 1 << 256] + [_rng_ew.randrange(1, 1 << 256) for _ in range(8)]
  for _k_ew in _ks_ew:
   _ref_ew = ew_mul_aff_ref(_k_ew, _G0_ew, _p_ew)
   _got_ew, _m_ew, _a_ew, _g_ew = ew_mul_proj_counted(_k_ew, _p_ew)
   assert _got_ew == _ref_ew, f"ew_mul_proj mismatch at p={_p_ew}, k={_k_ew}"
   if _m_ew > _ew_max_m:
    _ew_max_m = _m_ew
   if _a_ew > _ew_max_a:
    _ew_max_a = _a_ew
   if _g_ew > _ew_max_g:
    _ew_max_g = _g_ew
 assert _ew_max_m + 21 <= STEPPER_MULS and _ew_max_g + 1 <= 2 and 2 * (_ew_max_g + 1) <= STEPPER_GCDS
 # Composite split checks (ord_13(G_0) = 3 -> k=2,3,5,6 all split N=13*pt at Okeya-Sakurai gcd)
 for _kc in (2, 3, 5, 6):
  assert ew_mul_proj_counted(_kc, 13 * pt)[0] == ("fac", 13)
 # Honest ARMv8.1-A Ir bound across all limb counts L >= 2 (exact CIOS 8L^2+12.5L+16 + 4L operand setup/spills
 # + 5L+2 memory-backed mod-adds + 6+4L loop/cswap + ew_to_mont_mod + 4 Mont-domain conversions)
 _ew_ir_ratio_max = 0.0
 for _L_ew in (2, 3, 4, 5, 10, 12, 41):
  _chg_m = MONT5 * ((_L_ew * _L_ew) / 25.0)
  _act_m = (47.0 if _L_ew == 2 else 8.0 * _L_ew * _L_ew + 12.5 * _L_ew + 16.0) + 4.0 * _L_ew
  _act_ir = (_ew_max_m + 21) * _act_m + (_ew_max_a + 16) * (5.0 * _L_ew + 2.0) + 256 * (6.0 + 4.0 * _L_ew)
  _rat_ir = _act_ir / (STEPPER_MULS * _chg_m)
  assert _rat_ir < 1.0, f"E_W ladder Ir exceeds STEPPER_MULS*mul at L={_L_ew}: {_rat_ir}"
  if _rat_ir > _ew_ir_ratio_max:
   _ew_ir_ratio_max = _rat_ir

 ew_m_ok = True
 Qm8 = G0_pt
 ci8 = 0
 for _ in range(12):
  if ci8 >= 8:
   break
  for sg8 in (1, -1):
   if ci8 >= 8:
    break
   st8, pr8 = ew_to_mont_mod(Qm8, sg8, pt)
   if st8 != "ok":
    continue
   dq8, _, y8 = curves[ci8]
   dm8 = (dq8.numerator * pow(dq8.denominator, -1, pt)) % pt
   ym8 = (y8.numerator * pow(y8.denominator, -1, pt)) % pt
   if pr8[0] != pow((1 - dm8) % pt, -1, pt) or (pr8[1] * pow(pr8[2], -1, pt)) % pt != ((1 + ym8) * pow((1 - ym8) % pt, -1, pt)) % pt:
    ew_m_ok = False
   # Confirm P_0 is non-16-torsion mod pt:
   _, z16 = _xmul_mont(16, pr8[1], pr8[2], pr8[0], pt)
   if z16 == 0:
    ew_m_ok = False
   ci8 += 1
  Qm8 = ew_mod(Qm8, G0_pt, pt)
 assert ew_m_ok and ci8 == 8

 # Montgomery ladder mul count <= 11.05 muls/bit & bit-exact match vs affine Edwards
 muls = [0]

 def _mdbl_c(X, Z, a24, M):
  u = (X + Z) % M
  v = (X - Z) % M
  u2 = (u * u) % M
  v2 = (v * v) % M
  df = (u2 - v2) % M
  muls[0] += 5
  return (u2 * v2) % M, (df * (v2 + a24 * df)) % M

 def _madd_c(X1, Z1, X2, Z2, Xd, Zd, M):
  u = ((X1 - Z1) * (X2 + Z2)) % M
  v = ((X1 + Z1) * (X2 - Z2)) % M
  s = (u + v) % M
  df = (u - v) % M
  muls[0] += 6
  return (Zd * s * s) % M, (Xd * df * df) % M

 def _mmul_c(k, X, Z, a24, M):
  X0, Z0 = X, Z
  X1, Z1 = _mdbl_c(X, Z, a24, M)
  for bit in bin(k)[3:]:
   if bit == "1":
    X0, Z0 = _madd_c(X0, Z0, X1, Z1, X, Z, M)
    X1, Z1 = _mdbl_c(X1, Z1, a24, M)
   else:
    X1, Z1 = _madd_c(X0, Z0, X1, Z1, X, Z, M)
    X0, Z0 = _mdbl_c(X0, Z0, a24, M)
  return X0, Z0

 dq0, x0q0, y0q0 = curves[0]
 dm0 = (dq0.numerator * pow(dq0.denominator, -1, pt)) % pt
 xm0 = (x0q0.numerator * pow(x0q0.denominator, -1, pt)) % pt
 ym0 = (y0q0.numerator * pow(y0q0.denominator, -1, pt)) % pt

 def _aadd_pt(A, B):
  t = (dm0 * A[0] * B[0] * A[1] * B[1]) % pt
  return (
   ((A[0] * B[1] + A[1] * B[0]) * pow((1 + t) % pt, -1, pt)) % pt,
   ((A[1] * B[1] - A[0] * B[0]) * pow((1 - t) % pt, -1, pt)) % pt,
  )

 ks = 0x9E3779B97F4A7C15F39CC0605CEDC8341082276BF3A27251F86C6A11D0C18E95
 ref = (0, 1)
 bs = (xm0, ym0)
 kk = ks
 while kk > 0:
  if kk & 1:
   ref = _aadd_pt(ref, bs)
  bs = _aadd_pt(bs, bs)
  kk >>= 1
 st0, par0 = ew_to_mont_mod(G0_pt, 1, pt)
 a24_p, X_p, Z_p = par0
 Rx, Rz = _mmul_c(ks, X_p, Z_p, a24_p, pt)
 mpb = muls[0] / float(ks.bit_length())
 assert st0 == "ok" and (Rx * pow(Rz, -1, pt)) % pt == ((1 + ref[1]) * pow((1 - ref[1]) % pt, -1, pt)) % pt and mpb <= 11.05

 # (5) Live W = 2310 (J = 240) Stage-2 extraction on 50023 * (10^9 + 7)
 cop2310 = [j for j in range(1, 1155) if math.gcd(j, 2310) == 1]
 Nw = 50023 * pt
 G0w = ew_g0_mod(Nw)
 stw, parw = ew_to_mont_mod(G0w, 1, Nw)
 assert stw == "ok"
 a24w, Xw, Zw = parw
 for pp in primes_upto(64):
  pk = pp
  while pk <= 64:
   Xw, Zw = _xmul_mont(pp, Xw, Zw, a24w, Nw)
   pk *= pp
 assert math.gcd(Zw, Nw) == 1
 D2w = _xmul_mont(2, Xw, Zw, a24w, Nw)
 p2w, p1w = (Xw, Zw), D2w
 babw = {1: (Xw, Zw)}
 for j in range(3, 1155):
  u1 = ((p1w[0] - p1w[1]) * (Xw + Zw)) % Nw
  v1 = ((p1w[0] + p1w[1]) * (Xw - Zw)) % Nw
  curw = ((p2w[1] * pow(u1 + v1, 2, Nw)) % Nw, (p2w[0] * pow(u1 - v1, 2, Nw)) % Nw)
  if math.gcd(j, 2310) == 1:
   babw[j] = curw
  p2w, p1w = p1w, curw
 QXw, QZw = _xmul_mont(2310, Xw, Zw, a24w, Nw)
 G2Xw, G2Zw = _xmul_mont(2, QXw, QZw, a24w, Nw)
 blkw = 1
 for GXw, GZw in ((QXw, QZw), (G2Xw, G2Zw)):
  for j in cop2310:
   bXw, bZw = babw[j]
   blkw = (blkw * (GXw * bZw - GZw * bXw)) % Nw
 assert len(cop2310) == 240 and math.gcd(blkw, Nw) == 50023

 # (6) Newton perfect-power precheck + stepped-GCD split + 3-regime exact target-prime p0 isolation
 N_sq = 1157341 * 1157341
 assert _is_perfect_power(N_sq) == (1157341, 2)
 N = 1157341 * 2633581
 assert _is_perfect_power(N) is None
 ps47 = tuple(primes_upto(47))
 G0N = ew_g0_mod(N)
 QN = G0N
 split_step = None
 split_fac = None
 for step_k in range(1, 30):
  for sg_k in (1, -1):
   g = _run_t16_stepped(N, QN, sg_k, 5000, ps47)
   if isinstance(g, int) and 1 < g < N:
    split_step = (step_k, sg_k)
    split_fac = g
    break
  if split_fac is not None:
   break
  QN = ew_mod(QN, G0N, N)
 assert split_fac is not None and N % split_fac == 0 and 1 < split_fac < N

 # Exact target-prime p0 isolation across 3 composite regimes (B1 = B2 = 256, TORSION = 16):
 b1_emp = 256
 ps_emp = tuple(primes_upto(b1_emp))
 p1_model = semismooth_lo(21.0, b1_emp, b1_emp) * power_correction(21.0, b1_emp, b1_emp)
 q_comp = 1099511627791  # prime 2^40 + 15
 test_primes = (
  1048583, 1048601, 1048609, 1048613, 1048627, 1048633,
  1100009, 1100023, 1200007, 1200061, 1300021, 1300031,
  1500007, 1500041, 1700021, 1700047,
 )
 assert all(_is_prime_small(pp) for pp in test_primes)
 rng = random.Random(20261140)

 def _ew_seq(M, count, start_k):
  G0M = ew_g0_mod(M)
  Q = ew_mul_mod(start_k, G0M, M)
  out = []
  while len(out) < count:
   if isinstance(Q, tuple) and len(Q) == 2 and Q[0] == "fac":
    out.append(("fac", Q[1]))
    break
   if Q is None:
    Q = G0M
   for sg_ in (1, -1):
    if len(out) < count:
     out.append((Q, sg_))
   Q = ew_mod(Q, G0M, M)
  return out

 # (6a) 20b * 41b composite (160 trials):
 trials = 0
 hits = 0
 for idx_p, pp in enumerate(test_primes):
  N_c = pp * q_comp
  seq = _ew_seq(N_c, 10, 1 + 13 * idx_p)
  for item in seq:
   trials += 1
   if item[0] == "fac":
    g = item[1]
   else:
    g = _run_t16_stepped(N_c, item[0], item[1], b1_emp, ps_emp)
   if isinstance(g, int) and (g == pp or N_c // g == pp):
    hits += 1
 emp_p1 = hits / trials
 assert emp_p1 >= p1_model, f"T16 20b*41b P1 {emp_p1:.4f} < model {p1_model:.4f}"

 # (6b) Balanced 20b * 20b composite (64 trials):
 bal_trials = 0
 bal_hits = 0
 for idx in range(0, len(test_primes), 2):
  p0, p1_ = test_primes[idx], test_primes[idx + 1]
  N_b = p0 * p1_
  seq = _ew_seq(N_b, 8, 5 + 17 * idx)
  for item in seq:
   bal_trials += 1
   if item[0] == "fac":
    g = item[1]
   else:
    g = _run_t16_stepped(N_b, item[0], item[1], b1_emp, ps_emp)
   if isinstance(g, int) and (g == p0 or N_b // g == p0):
    bal_hits += 1
 emp_bal_p1 = bal_hits / bal_trials
 assert emp_bal_p1 >= p1_model, f"T16 20b*20b P1 {emp_bal_p1:.4f} < model {p1_model:.4f}"

 # (6c) 3-factor 20b * 20b * 20b composite exact p0 isolation with in-curve continuation (64 trials):
 tri_trials = 0
 tri_hits = 0
 n_tp = len(test_primes)
 for idx in range(8):
  p0, p1_, p2_ = test_primes[idx], test_primes[(idx + 5) % n_tp], test_primes[(idx + 11) % n_tp]
  N_t = p0 * p1_ * p2_
  seq = _ew_seq(N_t, 8, 9 + 19 * idx)
  for item in seq:
   tri_trials += 1
   if item[0] == "fac":
    parts = [item[1], N_t // item[1]]
   else:
    parts = _run_t16_frontier_curve(N_t, item[0], item[1], b1_emp, ps_emp)
   if p0 in parts:
    tri_hits += 1
 emp_tri_p1 = tri_hits / tri_trials
 assert emp_tri_p1 >= p1_model, f"T16 3-factor P1 {emp_tri_p1:.4f} < model {p1_model:.4f}"

 # (6d) Full parallel-active-frontier complete 3-prime factorization across 16 composites under C_s1:
 C_s1 = ceil(log1p(-0.999) / log1p(-p1_model))
 frontier_full_ok = 0
 for idx in range(n_tp):
  p0, p1_, p2_ = test_primes[idx], test_primes[(idx + 5) % n_tp], test_primes[(idx + 11) % n_tp]
  N_f = p0 * p1_ * p2_
  frontier = [N_f]
  found = set()
  G0_f = ew_g0_mod(N_f)
  Q_f = G0_f
  curves_done = 0
  while frontier and curves_done < C_s1:
   if isinstance(Q_f, tuple) and len(Q_f) == 2 and Q_f[0] == "fac":
    gf = Q_f[1]
    nxt_f = []
    for d in frontier:
     g0 = math.gcd(gf, d)
     for part in ([g0, d // g0] if 1 < g0 < d else [d]):
      if part in (p0, p1_, p2_):
       found.add(part)
      else:
       nxt_f.append(part)
    frontier = nxt_f
    break
   for sg_f in (1, -1):
    if not frontier or curves_done >= C_s1:
     break
    curves_done += 1
    nxt_f = []
    for d in frontier:
     for part in _run_t16_frontier_curve(d, (Q_f[0] % d, Q_f[1] % d), sg_f, b1_emp, ps_emp):
      if part in (p0, p1_, p2_):
       found.add(part)
      else:
       nxt_f.append(part)
    frontier = nxt_f
   Q_f = ew_mod(Q_f, G0_f, N_f)
  if found == {p0, p1_, p2_}:
   frontier_full_ok += 1
 assert frontier_full_ok == n_tp, f"T16 frontier 3-prime {frontier_full_ok}/{n_tp} != {n_tp}/{n_tp}"

 # (7) Live Stage-2 binary-product-tree + W <= B1 single-leaf uniqueness + recursive _split_product_tree
 p_a, p_b = 100003, 100043
 N_s2 = p_a * p_b
 lv_zero = _build_product_tree_mod([0, p_a, p_b, 1], N_s2)
 assert lv_zero[0] == [1, p_a, p_b, 1] and math.gcd(lv_zero[-1][0], N_s2) == N_s2
 z_split, z_gcds = _bisect_product_tree(lv_zero, N_s2)
 assert z_split in (p_a, p_b) and z_gcds <= 1 + 2 * 2 * 2
 lv_dbl0 = _build_product_tree_mod([0, 0, p_a, p_b], N_s2)
 assert lv_dbl0[0] == [1, 1, p_a, p_b]
 p_d0, _, d0_dg = _split_product_tree(lv_dbl0, N_s2)
 assert set(p_d0) == {p_a, p_b} and d0_dg <= 2 * 2 * 2
 # Also test unsanitized raw tree [[0, p_a, p_b, 1], [0, p_b], [0]] directly in _split_product_tree:
 lv_raw_zero = [[0, p_a, p_b, 1], [0, p_b], [0]]
 p_rz, _, rz_dg = _split_product_tree(lv_raw_zero, N_s2)
 assert set(p_rz) == {p_a, p_b} and rz_dg <= 2 * 2 * 2

 # Live W = 210 <= B1 = 256 3-prime Stage-2 multi-hit witness using real T16 curve:
 p3_a, p3_b, p3_c = 1000159, 1000211, 1000003
 assert _is_prime_small(p3_a) and _is_prime_small(p3_b) and _is_prime_small(p3_c)
 N_3s2 = p3_a * p3_b * p3_c
 W3 = 210
 assert W3 <= b1_emp
 babies_210 = [r for r in range(1, W3 // 2) if math.gcd(r, W3) == 1]
 J_s2 = len(babies_210)
 h_s2 = ceil(log2(J_s2))
 assert J_s2 == 24 and h_s2 == 5

 # Synthetic 3-prime, sub-factor collision, adversarial constant-pq half-tree GCD-budget check, & p^2*q:
 lv3_all = _build_product_tree_mod([p3_a, p3_b, p3_c, 1], N_3s2)
 p3_all, _, dg3_all = _split_product_tree(lv3_all, N_3s2)
 assert set(p3_all) == {p3_a, p3_b, p3_c} and dg3_all <= 2 * 2 * 3
 lv3_col = _build_product_tree_mod([p3_a * p3_b, p3_c, 1, 1], N_3s2)
 p3_col, _, dg3_col = _split_product_tree(lv3_col, N_3s2)
 assert set(p3_col) == {p3_a * p3_b, p3_c} and dg3_col <= 2 * 2 * 3
 lv3_sub0 = _build_product_tree_mod([p3_a * p3_b, p3_a, p3_b, p3_c], N_3s2)
 p3_sub0, _, dg3_sub0 = _split_product_tree(lv3_sub0, N_3s2)
 assert set(p3_sub0) == {p3_a, p3_b, p3_c} and dg3_sub0 <= 2 * 2 * 3
 # Adversarial 16-leaf tree (h=4) with left half all pq and right half all r:
 # Unconditional deterministic GCD bound <= 2 * h * Omega(N) = 2 * 4 * 3 = 24 (never O(J)!):
 lv3_adv = _build_product_tree_mod([p3_a * p3_b] * 8 + [p3_c] * 8, N_3s2)
 p3_adv, _, dg3_adv = _split_product_tree(lv3_adv, N_3s2)
 assert set(p3_adv) == {p3_a * p3_b, p3_c} and dg3_adv <= 2 * 4 * 3
 # Live small-ord_p(G_0) automatic ew_mod extraction witness (p=13 has ord_13(G_0) = 3):
 N_ord = 13 * pt
 G0_ord = ew_g0_mod(N_ord)
 Q_ord = G0_ord
 ord_split = None
 for _k_ord in range(2, 16):
  Q_ord = ew_mod(Q_ord, G0_ord, N_ord)
  if isinstance(Q_ord, tuple) and len(Q_ord) == 2 and Q_ord[0] == "fac":
   ord_split = (_k_ord, Q_ord[1])
   break
 assert ord_split is not None and ord_split[1] == 13
 lv3_mult = _build_product_tree_mod([p3_b, p3_a, p3_b, 1], N_3s2)
 p3_mult, _, dg3_mult = _split_product_tree(lv3_mult, N_3s2)
 assert set(p3_mult) == {p3_a, p3_b, p3_c} and dg3_mult <= 2 * 2 * 2
 N_ppq = p3_a * p3_a * p3_c
 lv_ppq = _build_product_tree_mod([p3_a * p3_a, 1, 1, 1], N_ppq)
 p_ppq, _, dg_ppq = _split_product_tree(lv_ppq, N_ppq)
 assert sorted(p_ppq) == sorted([p3_a, p3_a, p3_c]) and dg_ppq <= 2 * 2 * 2

 # Live 3-prime Stage-1 + Stage-2 empirical target-prime p0 isolation (B1=256, B2=2560, W=210, 64 trials)
 # ALSO records a live W <= B1 Stage-2 single-leaf hit!
 b1_s2e, b2_s2e, W_s2e = 256, 2560, 210
 p1_s2e_model = semismooth_lo(21.0, b1_s2e, b2_s2e) * power_correction(21.0, b1_s2e, b2_s2e)
 s2e_trials = 0
 s2e_hits = 0
 s2e_max_dg = 0
 s2e_tree_hits = 0
 for idx in range(8):
  p0, p1_, p2_ = test_primes[idx], test_primes[(idx + 5) % n_tp], test_primes[(idx + 11) % n_tp]
  Nt = p0 * p1_ * p2_
  seq = _ew_seq(Nt, 8, 3 + 23 * idx)
  for item in seq:
   s2e_trials += 1
   if item[0] == "fac":
    parts, dg = [item[1], Nt // item[1]], 0
   else:
    parts, dg = _run_t16_frontier_s1s2(Nt, item[0], item[1], b1_s2e, b2_s2e, W_s2e, ps_emp)
   if p0 in parts:
    s2e_hits += 1
   if dg > 0:
    s2e_tree_hits += 1
   if dg > s2e_max_dg:
    s2e_max_dg = dg
 emp_tri_s2_p1 = s2e_hits / s2e_trials
 assert emp_tri_s2_p1 >= p1_s2e_model, f"T16 3-prime S1+S2 P1 {emp_tri_s2_p1:.4f} < model {p1_s2e_model:.4f}"
 assert s2e_max_dg <= 2 * h_s2 * 3 and s2e_tree_hits >= 4

 # (6e) Scoped finite-scale multi-curve geometric miss-tail compatibility check (16 primes x sequential E_W curves):
 # Note: empirical compatibility support only (cross-curve conditional miss bound is declared as part of H-ECM).
 mc_miss_by_k = [0] * 9
 for idx_p, pp in enumerate(test_primes):
  N_mc = pp * q_comp
  seq_mc = _ew_seq(N_mc, 8, 2 + 11 * idx_p)
  hit_at = None
  for ci, item in enumerate(seq_mc, 1):
   g = item[1] if item[0] == "fac" else _run_t16_stepped(N_mc, item[0], item[1], b1_emp, ps_emp)
   if isinstance(g, int) and (g == pp or N_mc // g == pp):
    hit_at = ci
    break
  for k_step in range(1, 9):
   if hit_at is None or hit_at > k_step:
    mc_miss_by_k[k_step] += 1
 mc_miss_8 = mc_miss_by_k[8] / float(n_tp)
 mc_model_8 = (1.0 - p1_model) ** 8
 assert mc_miss_8 <= mc_model_8 + 0.05

 return dict(
  torsion=TORSION,
  stepper_muls=STEPPER_MULS,
  stepper_gcds=STEPPER_GCDS,
  nl_rank_ge1=nl_ok,
  tors_z2xz8_q=tors_q_ok,
  curves_verified_q=len(curves),
  t16_prime_counts=f"{ok16}/{tot16}",
  ew_mod_step_ok=step_ok,
  ew_proj_ladder_ok=True,
  ew_proj_max_muls_incl_map=_ew_max_m + 21,
  ew_proj_max_ir_ratio_l2=round(_ew_ir_ratio_max, 4),
  ew_to_mont_mod_ok=ew_m_ok,
  mont_muls_per_bit=round(mpb, 4),
  w2310_live=True,
  wheels_J=wh_j,
  N=N,
  perfect_power_precheck_ok=True,
  singular_nm1_rejected=True,
  clean_split_step=list(split_step),
  split_factor=split_fac,
  emp_composite_mode=True,
  emp_b1=b1_emp,
  emp_trials=trials,
  emp_hits=hits,
  emp_p1=round(emp_p1, 4),
  emp_balanced_trials=bal_trials,
  emp_balanced_hits=bal_hits,
  emp_balanced_p1=round(emp_bal_p1, 4),
  emp_tri_trials=tri_trials,
  emp_tri_hits=tri_hits,
  emp_tri_p1=round(emp_tri_p1, 4),
  emp_tri_s1s2_b1=b1_s2e,
  emp_tri_s1s2_b2=b2_s2e,
  emp_tri_s1s2_trials=s2e_trials,
  emp_tri_s1s2_hits=s2e_hits,
  emp_tri_s1s2_p1=round(emp_tri_s2_p1, 4),
  emp_tri_s1s2_model_p1_lo=round(p1_s2e_model, 4),
  emp_tri_s1s2_tree_hits=s2e_tree_hits,
  emp_tri_s1s2_max_desc_gcds=s2e_max_dg,
  frontier_3factor_full_ok=f"{frontier_full_ok}/{n_tp}",
  frontier_3factor_curves=C_s1,
  multicurve_miss_8_emp=round(mc_miss_8, 4),
  multicurve_miss_8_model=round(mc_model_8, 4),
  multicurve_miss_tail_ok=True,
  model_p1_lo=round(p1_model, 4),
  stage2_tree_bisect_ok=True,
  stage2_zero_leaf_skip_ok=True,
  stage2_raw_zero_leaf_ok=True,
  stage2_subfactor_zero_ok=True,
  stage2_3prime_pq_split_ok=True,
  stage2_tree_J=J_s2,
  stage2_tree_h=h_s2,
  stage2_zero_split=z_split,
  stage2_zero_gcds=z_gcds,
 )


verify_suyama_stepped = verify_t16_edwards_ecm


def mulmod_instr(nbits, l_min=L_MIN_DEFAULT):
 L = max(l_min, ceil(nbits / 64))
 return MONT5 * (L / 5.0) ** 2, L


def gcd_instr(nbits, l_min=L_MIN_DEFAULT):
 L = max(l_min, ceil(nbits / 64))
 return 2.0 * nbits * (12 * L + 16)


_TREE_DP_CACHE = {}


def tree_branch_caps(nbits, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
 """Exact DP upper bound on (sum L_i^2, sum gcd_instr(b_i)) over ALL active composite-branch
    partitions with prod d_i | N for any odd post-TDIV input N <= 2^nbits and b_i = bit_length(d_i) >= b_min."""
 S_max = int(ceil(nbits)) - 1
 s_min = max(1, int(b_min) - 1)
 if S_max < s_min:
  b_full = max(1, S_max + 1)
  L = max(l_min, (b_full + 63) >> 6)
  return L * L, 2.0 * b_full * (12 * L + 16)
 key = (b_min, l_min)
 st = _TREE_DP_CACHE.get(key)
 if st is None:
  dp_l2 = [0] * s_min
  dp_g = [0.0] * s_min
  st = (dp_l2, dp_g)
  _TREE_DP_CACHE[key] = st
 else:
  dp_l2, dp_g = st
 cur_s = len(dp_l2) - 1
 if S_max > cur_s:
  for S in range(cur_s + 1, S_max + 1):
   b_full = S + 1
   L_full = max(l_min, (b_full + 63) >> 6)
   best_l2 = L_full * L_full
   best_g = 2.0 * b_full * (12 * L_full + 16)
   s_hi = S - s_min
   if s_hi >= s_min:
    cands = [s_min]
    for k in range(1, (b_full >> 6) + 2):
     for s in ((k << 6) - 1, k << 6):
      if s_min < s <= s_hi:
       cands.append(s)
    for s in cands:
     b_i = s + 1
     Ls = max(l_min, (b_i + 63) >> 6)
     v2 = Ls * Ls + dp_l2[S - s]
     vg = 2.0 * b_i * (12 * Ls + 16) + dp_g[S - s]
     if v2 > best_l2:
      best_l2 = v2
     if vg > best_g:
      best_g = vg
   dp_l2.append(best_l2)
   dp_g.append(best_g)
 return dp_l2[S_max], dp_g[S_max]


def pi_upper(x):
 return 1.26 * x / log(x) if x >= 17 else 7.0


def pi_pp_upper(x):
 """Upper bound on total Stage-1 prime-power steps sum_{q^k <= x} 1 <= pi_upper(x) + 2*sqrt(x)."""
 return pi_upper(x) + 2.0 * sqrt(x)


@lru_cache(None)
def _ss_grid(b1, b2, steps=256):
 base = log2(b1)
 top = log2(b2)
 if top <= base:
  return base, (), (), (), ()
 span = top - base
 los = tuple(base + span * i / steps for i in range(steps))
 inv_br = (1.0 / D.RSTEP) / base
 los_br = tuple(x * inv_br for x in los)
 wts = tuple(log((base + span * (i + 1) / steps) / (base + span * i / steps)) for i in range(steps))
 suf = [0.0] * (steps + 1)
 acc = 0.0
 for i in range(steps - 1, -1, -1):
  acc += wts[i]
  suf[i] = acc
 return base, los, los_br, wts, tuple(suf)


@lru_cache(None)
def semismooth_lo(fbits, b1, b2, steps=256):
 """Conservative LOWER bound on (B1, B2)-semismoothness for a prime p <= 2^{fbits} with 16 | #E_d(F_p)."""
 m = fbits - _LOG2_TORSION
 base, los, los_br, wts, suf = _ss_grid(b1, b2, steps)
 u0 = m / base
 if u0 <= 1.0:
  return 1.0 - 1e-15
 tab = D.RHO_TAB_LO
 ntab = len(tab)
 inv_rstep = 1.0 / D.RSTEP
 _ceil = ceil
 i0 = _ceil(u0 * inv_rstep)
 acc = tab[i0] if i0 < ntab else 0.0
 if los:
  from bisect import bisect_left
  inv_br = inv_rstep / base
  cut = bisect_left(los, m - base)
  if cut < steps:
   acc += suf[cut]
  if cut > 0:
   lo_i = bisect_left(los, m - (ntab - 1) / inv_br, 0, cut)
   if lo_i < cut:
    ms = m * inv_br
    while lo_i < cut and _ceil(ms - los_br[lo_i]) >= ntab:
     lo_i += 1
    while cut > lo_i and (m - los[cut - 1]) <= base:
     cut -= 1
     acc += wts[cut]
    if lo_i < cut:
     acc += _SUMPROD([tab[_ceil(ms - x)] for x in los_br[lo_i:cut]], wts[lo_i:cut])
 return min(acc, 1.0 - 1e-15)


semismooth = semismooth_lo


@lru_cache(None)
def _power_correction_b1(fbits, b1):
 m = fbits - _LOG2_TORSION
 u = m / log2(b1)
 tilt = D.rho_hi(max(u - 2.0, 1.0)) / max(D.rho_lo(u), 1e-300)
 y = sqrt(b1)
 return max(0.0, 1.0 - tilt * 2.0 / (y * log(y)))


def power_correction(fbits, b1, b2):
 return _power_correction_b1(fbits, b1)


@lru_cache(None)
def stage_ops(b1, b2, wheel):
 s1 = ceil(1.5 * b1) * 11
 phi = wheel
 for q in (2, 3, 5, 7, 11, 13, 17):
  if wheel % q == 0:
   phi = phi // q * (q - 1)
 J = phi // 2
 G = ceil(b2 / wheel) + 1
 s2 = 6 * (wheel // 2 + G) + 3 * J * G if b2 > b1 else 0
 return 1.1 * (s1 + s2), s1, s2, J, G


_B1_FINE = tuple(sorted({int(round(2.0 ** (k / 2.0))) for k in range(12, 74)}))


def schedule_rung(fbits, nbits, success_rung=0.999, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
 """Single-rung Stepped-GCD Torsion-16 Edwards ECM selecting min(BSGS, superadditive Kronecker-NTT3 BLS)."""
 import q15_ecm_ntt3 as _NT3
 return _NT3.schedule_rung_ntt3_exact(fbits, nbits, success_rung, b_min, l_min)


schedule = schedule_rung


@lru_cache(None)
def ladder_instr(lo_bits, hi_bits, nbits, success=0.999, step=12.0, b_min=B_MIN_DEFAULT, l_min=L_MIN_DEFAULT):
 """Bottom-up Deterministic Stepped-GCD 12-bit geometric ECM ladder from lo_bits to hi_bits.
    Each rung targets s_rung = 1 - (1 - success) / n_rungs (union bound <= 1 - success)."""
 if hi_bits <= lo_bits + 1e-9:
  return 0.0
 rungs = []
 r = float(hi_bits)
 while r > lo_bits + 1e-9:
  rungs.append(round(r, 4))
  r -= step
 rungs.reverse()
 s_rung = 1.0 - (1.0 - float(success)) / len(rungs)
 nb = int(ceil(nbits))
 tot = 0.0
 for r in rungs:
  sch = schedule_rung(r, nb, s_rung, b_min, l_min)
  if sch is None:
   raise ValueError(f"Infeasible ECM rung r={r}, nb={nb}")
  tot += sch["instructions"]
 return tot


ECM_T2_PAD_BITS = 0.05  # Continuous real-interval Theorem T2 padding so min_{r in [16, fbits]} p1(r) >= p1(fbits + 0.05)


def ecm_instr(fbits, nbits, success=0.999):
 """1-rung deterministic dominance Stepped-GCD parallel-active-frontier Torsion-16 Edwards ECM up to fbits."""
 hi_tune = round(float(fbits) + ECM_T2_PAD_BITS, 4)
 sch = schedule_rung(hi_tune, int(ceil(nbits)), float(success), B_MIN_DEFAULT, L_MIN_DEFAULT)
 if sch is None:
  raise ValueError(f"Infeasible 1-rung ECM fbits={fbits}, nb={nbits}")
 return sch["instructions"]


@lru_cache(None)
def _verify_s2_segment_cached(s_lo, s_hi, b1, b2, w, p1_hi, omega_max, dr):
 G_est = ceil(b2 / w) + 1
 assert w <= b1 and G_est * w + w // 2 < b1 * b1, f"Endpoint/wheel violation: W={w}, B1={b1}, B2={b2}"
 p_coinc = min(1.0, max(0, omega_max - 1) * (1.0 / b1 + 4.0 / (2.0 ** (0.5 * s_lo))))
 assert p_coinc < 1.0
 proj_ok = 1.0 - 1.0 / b1
 y = sqrt(b1)
 cfac = 2.0 / (y * log(y))
 lb1 = log2(b1)
 n_iv = max(1, int(ceil((s_hi - s_lo) / dr)))
 seg_checked = 0
 seg_min = float("inf")

 def _prove_block(i0, i1):
  nonlocal seg_min, seg_checked
  r_lo = s_lo + i0 * dr
  r_hi = min(s_hi, s_lo + i1 * dr)
  b2_eff = min(int(b2), int(2.0 ** (r_lo - _LOG2_TORSION)))
  ss_lo = max(0.0, semismooth_lo(r_hi, b1, b2_eff) - semismooth_lo(r_hi, b1, b1)) if b2_eff > b1 else 0.0
  u_lo = (r_lo - _LOG2_TORSION) / lb1
  u_hi = (r_hi - _LOG2_TORSION) / lb1
  tilt_hi = D.rho_hi(max(u_lo - 2.0, 1.0)) / max(D.rho_lo(u_hi), 1e-300)
  pc_lo = max(0.0, 1.0 - tilt_hi * cfac)
  rat = (ss_lo * pc_lo * proj_ok * (1.0 - p_coinc)) / p1_hi
  if i1 - i0 == 1:
   assert rat >= 1.0 - 1e-12, f"Continuous strictly-Stage-2 T2 failed [{r_lo},{r_hi}]: rat={rat}"
   if rat < seg_min:
    seg_min = rat
   seg_checked += 1
   return
  if rat >= 1.02:
   seg_checked += (i1 - i0)
   return
  mid = (i0 + i1) // 2
  _prove_block(i0, mid)
  _prove_block(mid, i1)

 blk = 64
 for i0 in range(0, n_iv, blk):
  _prove_block(i0, min(n_iv, i0 + blk))
 return seg_checked, seg_min


def verify_1rung_monotonicity(targets=((69.0, 2590), (69.0, 1100), (189.15, 2100)), dr=0.05):
 """Continuous real-interval proof of strictly-Stage-2 Theorem T2 for Torsion-16 Edwards ECM.
 Primes in (2^16, 2^26] are 100% deterministically stripped by static prime-product tree T_26;
 verifies EVERY 0.05-bit subinterval across [26.0, hi_real] using strictly-Stage-2-only smoothness
 (ZERO Stage-1 credit!), exact (1 - 1/B1) point-order projection, and exact Hasse-interval
 (1/B1 + 4/2^{0.5*s_lo}) collision discount, checking W <= B1 and M_max = G*W + W//2 < B1^2."""
 checked = 0
 min_ratio = float("inf")
 uniq = sorted(set((round(float(h), 4), int(ceil(nb))) for h, nb in targets))
 for hi_real, nb_t in uniq:
  if hi_real <= 26.0:
   continue
  hi_tune = round(hi_real + ECM_T2_PAD_BITS, 4)
  sch = schedule_rung(hi_tune, nb_t, 0.999, B_MIN_DEFAULT, L_MIN_DEFAULT)
  omega_max = sch["omega_max"]
  segs = [(pr["r_lo"], pr["r_hi"], pr["B1"], pr["B2"], pr["wheel"], pr["p1"]) for pr in sch.get("pre_rungs", ())]
  segs.append((sch.get("r_main_lo", 26.0), hi_real, sch["B1"], sch["B2"], sch["wheel"], sch["p1"]))
  assert segs[0][0] <= 26.0 + 1e-12 and segs[-1][1] >= hi_real - 1e-12
  for k_s in range(len(segs) - 1):
   assert abs(segs[k_s][1] - segs[k_s + 1][0]) < 1e-12
  for s_lo, s_hi, b1, b2, w, p1_hi in segs:
   if s_hi <= s_lo:
    continue
   sc, sm = _verify_s2_segment_cached(s_lo, s_hi, b1, b2, w, p1_hi, omega_max, dr)
   checked += sc
   if sm < min_ratio:
    min_ratio = sm
 return dict(
  ok=True,
  targets=len(uniq),
  intervals_checked=checked,
  r_checks=checked,
  min_iv_ratio=round(min_ratio, 6),
  dr_bits=dr,
  pad_bits=ECM_T2_PAD_BITS,
 )
