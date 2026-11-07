import contextlib, functools, io, json, os, random, sys, time
from collections import defaultdict
from fractions import Fraction
T0 = time.time()
cfg = json.load(sys.stdin)
SEED = cfg.get("seed", 1)
BUDGET = cfg.get("seconds", 12)
DO_KERNELS = cfg.get("kernels", True)
RESERVE = 0.5
T_END = T0 + BUDGET - RESERVE
MIN_RANDOM, MIN_NEAR, MIN_CLASSES = 2000, 600, 3000
def package_root():
 here = os.path.dirname(os.path.abspath(__file__))
 for base in (here, os.path.join(here, "files"), "/snapshot/files"):
  if os.path.isfile(os.path.join(base, "kernels", "a64.py")): return base
 raise FileNotFoundError("kernels/a64.py not found")
HERE = package_root()
sys.path.insert(0, os.path.join(HERE, "kernels"))
import a64
sys.path.insert(0, HERE)
import sites, machine
M64 = 2 ** 64 - 1
M41, M42 = 2 ** 41 - 1, 2 ** 42 - 1
PART = 106
M128 = 2 ** 128
p = 2 ** 256 - 2 ** 32 - 977
r = int("fffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 16)
GX = int("79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 16)
GY = int("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16)
BETA = int("7ae96a2b657c07106e64479eac3434e99cf0497512f58995c1396c28719501ee", 16)
LAM = int("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72", 16)
C0 = 2 ** 32 + 977
FB = int("e4437ed6010e88286f547fa90abfe4c3", 16)
FA = int("114ca50f7a8e2f3f657c1108d9d44cfd9", 16)
FRAC_OFF = ((64, 80, 128, 144, 192, 208), (96, 112, 160, 176, 224, 240))
ENTRY, ACC, KEYS0, MX, Q2M = 512, 32, 64, 12, 2 ** 25 - 1
YR_B4, YR_B, YR_B2, YR_YROW, YR_NT = 256, 288, 64, 128, 320
W7_SHIFTS = (0, 37, 74, 111, 148, 185, 222)
W7_MASKS = ((1 << 37) - 1,) * 6 + ((1 << 34) - 1,)
BT = (1, BETA, BETA * BETA % p)
TERMS = dict(e=1, c_low=7, d=56, total=64)
ERR = ((-132, 194), (-132, 194), (-389, 264))
D_EDGE = ((194, 132), (194, 132), (264, 389))
OUT, CERT_LINES = {}, []
class Curve:
 def __init__(self, q, order): self.q, self.r = q, order
 def add(self, P, Q):
  q = self.q
  if P is None: return Q
  if Q is None: return P
  if P[0] == Q[0]:
   if (P[1] + Q[1]) % q == 0: return None
   lam = 3 * P[0] * P[0] * pow(2 * P[1], -1, q) % q
  else: lam = (Q[1] - P[1]) * pow(Q[0] - P[0], -1, q) % q
  x = (lam * lam - P[0] - Q[0]) % q
  return (x, (lam * (P[0] - x) - P[1]) % q)
 def neg(self, P): return None if P is None else (P[0], (-P[1]) % self.q)
 def mul(self, k, P):
  k %= self.r; acc = None
  while k:
   if k & 1: acc = self.add(acc, P)
   P = self.add(P, P); k >>= 1
  return acc
PROD = Curve(p, r); GP = (GX, GY)
def limbs(v, n=4): return [(v >> (64 * k)) & M64 for k in range(n)]
def put(mem, a, v, n=4):
 for k, w in enumerate(limbs(v, n)): mem[a + 8 * k] = w
def get(mem, a, n=4): return sum(mem[a + 8 * k] << (64 * k) for k in range(n))
def reg_int(regs, names): return sum(regs[n] << (64 * k) for k, n in enumerate(names))
def as_regs(v, names): return {n: w for n, w in zip(names, limbs(v, len(names)))}
def exact_key(x):
 v = [x, BETA * x % p, BETA * BETA % p * x % p]
 j = min(range(3), key=lambda t: v[t])
 return v[j], j
def t128(v): return (v << 128) // p
def comp(v): return (v << PART) // p
def label(v):
 c1, c2 = comp(v), comp(BETA * v % p)
 return c1 & M42, c1 >> 42 & M42, c1 >> 42, c2 >> 42
def class_label(x): return label(exact_key(x)[0])
def neg_t(v): return -t128(v % p) % M128
def corr(v): return neg_t(v), neg_t(BETA * v)
def fold(z): return (z % 2 ** 256 + (z >> 256) * C0) % 2 ** 256
def lkfold(z):
 z = z % 2 ** 256 + (z >> 256) * C0
 return z - 2 ** 256 + C0 if z >> 256 else z
def fs_op(a, b): return (a - b if a >= b else a - b + p) % 2 ** 256
def faddw_op(a, b): return (a + b if a + b < 2 ** 256 else a + b - p) % 2 ** 256
def shx_rows(c):
 R = [c]
 for _ in range(3): R.append((R[-1] % 2 ** 192 * 2 ** 64 + (R[-1] >> 192) * C0) % 2 ** 256)
 return R
def tfp(R, M, tm=None, sform=False):
 f = s = lo = d = 0
 for k in range(4):
  m, w = M >> 64 * k & M64, limbs(R[k])
  f += (m * w[3] & M64) + (m * w[2] >> 64)
  s += ((m * w[2] & M64) >> 3) + ((m * w[1] >> 64) >> 3)
  lo += m * w[0] + ((m * w[1] & M64) << 64)
  d += (m * w[2] & M64) % 8 + (m * w[1] >> 64) % 8
 A = ((f << 64) + 8 * s) % M128
 if tm is not None:
  Z = sum((M >> 64 * k & M64) * R[k] for k in range(4))
  e, tot = (t128(Z % p) - (Z >> 128)) % M128, (t128(Z % p) - A) % M128
  for k, v in (("e", e), ("c_low", lo >> 128), ("d", d), ("total", tot)): tm[k] = max(tm.get(k, 0), v)
  tm["sum_ok"] = tm.get("sum_ok", True) and tot == e + (lo >> 128) + d and s < 2 ** 64
 return A >> 64 << 64 | s if sform else A
def fracs(mem, e, s):
 v = [get(mem, e + o, 2) for o in FRAC_OFF[s]]
 for i in (4, 5): v[i] = v[i] >> 64 << 64 | (v[i] << 3 & M64)
 return v
def components(X, Y):
 a, b = X[0], Y[0]
 E, F, G = (b ** 4 + 2 * a ** 3 * b) % p, (a ** 4 + 2 * a * b ** 3) % p, 3 * a * a * b * b % p
 out = []
 for w in (pow(Y[1] + X[1], -2, p), pow(Y[1] - X[1], -2, p)): out += [E * w % p, F * w % p, G * w % p]
 return out
def frac_errs(mem, e, X, Y):
 cs, out = components(X, Y), []
 for s in (0, 1):
  ex = ((cs[3 * s] - Y[0]) % p, (cs[3 * s + 1] - X[0]) % p, cs[3 * s + 2])
  out += [(t128(BT[i % 2] * ex[i // 2] % p) - v + 2 ** 127) % M128 - 2 ** 127 for i, v in enumerate(fracs(mem, e, s))]
 return out
def leaf_w(f6, t):
 u = [[f6[2 * i], f6[2 * i + 1], -(f6[2 * i] + f6[2 * i + 1]) % M128] for i in range(3)]
 return [(u[0][m] + u[1][(m + t) % 3] + u[2][(m + 2 * t) % 3]) % M128 for m in (0, 1)]
def record_fields(w):
 f0, f1 = w[0] >> 64, w[1] >> 64
 s = ((w[0] + w[1]) >> 64) & M64
 j = (2 if f1 + s > M64 else 1) if f0 >= f1 else (2 if f0 + s > M64 else 0)
 n = ~(w[0] + w[1]) % M128
 a, b = ((w[0], w[1]), (w[1], n), (n, w[0]))[j]
 return j, a >> 22 & M42, a >> 64 & M42, a >> 64, b >> 64
def in_W(x):
 v = [x << 64, (BETA * x % p) << 64, (BETA * BETA % p * x % p) << 64]; top = p << 64
 return (any(w < p or w > top - p for w in v[:2]) or v[2] < p or v[2] > top - 4 * p or abs(v[0] - v[1]) < 2 * p
   or -p < v[0] - v[2] < 4 * p or -p < v[1] - v[2] < 4 * p)
def in_D(x):
 v = (x, BETA * x % p, BETA * BETA % p * x % p)
 return any(t128(w) % 2 ** 22 < lo or t128(w) % 2 ** 22 >= 2 ** 22 - hi for w, (lo, hi) in zip(v, D_EDGE))
def dft_leaf(c3, t, a, b):
 return (c3[0] - b + pow(BETA, t, p) * (c3[1] - a) + pow(BETA, 2 * t, p) * c3[2]) % p
def inverse_dft(L):
 i3, b2 = pow(3, -1, p), BETA * BETA % p
 return ((L[0] + L[1] + L[2]) * i3 % p, (L[0] + b2 * L[1] + BETA * L[2]) * i3 % p, (L[0] + BETA * L[1] + b2 * L[2]) * i3 % p)
def weak(rg): return rg.randrange(p) if rg.random() < 0.85 else rg.randrange(C0) + p
def represent(x, t, a, b, rg, tm=None):
 w, c, k = weak(rg) or 1, [rg.randrange(p) for _ in range(3)], rg.randrange(3); c[k] = 0
 v = (x + b + BT[t] * a) * pow(w, -1, p) - sum(BT[t * i % 3] * c[i] for i in range(3))
 c[k] = v * pow(BT[t * k % 3], -1, p) % p
 bw = BETA * w % p; bw += p if bw < C0 and rg.random() < 0.5 else 0
 out = []
 for ci, y in zip(c, (b, a, None)):
  R, cw = shx_rows(ci), corr(y) if y else (0, 0)
  out += [(tfp(R, w, tm) + cw[0]) % M128, (tfp(R, bw, tm) + cw[1]) % M128]
 return out, ((c[0] * w - b) % p, (c[1] * w - a) % p, c[2] * w % p)
def edge_leaf(rg, kind):
 m = rg.randrange(3)
 if kind == "D":
  r22 = rg.choice((rg.randrange(420), 2 ** 22 - 1 - rg.randrange(420), rg.randrange(2 ** 22)))
  v = -(-(rg.getrandbits(106) << 22 | r22) * p // M128)
  return v * pow(BETA, -m, p) % p or 1
 off = rg.randrange(6 << 192); w = off if rg.random() < 0.5 else p - 1 - off
 if kind == "Xm": return w * pow(BETA, -m, p) % p or 1
 return w * pow((1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA)[m] % p, -1, p) % p or 1
def deadline(frac): now = time.time(); return now + frac * max(0.0, T_END - now)
def near_tie_x(rg):
 mlt = rg.choice([1, BETA, BETA * BETA, 1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA])
 dd = rg.randrange(-(1 << rg.randrange(150, 250)), 1 << rg.randrange(150, 250))
 return dd * pow(mlt % p, -1, p) % p or 1
def index_check(x, rg):
 t, a, b = rg.randrange(3), rg.randrange(1, p), rg.randrange(1, p)
 j = record_fields(leaf_w(represent(x, t, a, b, rg)[0], t))[0]
 return in_W(x), j != exact_key(x)[1]
def class_labels(rg):
 out = dict(classes=0, edge_classes=0, members=[0, 0, 0], in_W=0, in_D=0, wrong_index=0, mismatches=0,
      errors=[[99, -99] for _ in range(3)], tfp_terms={})
 dl = deadline(0.05)
 while out["classes"] < MIN_CLASSES or (out["classes"] < 4 * MIN_CLASSES and time.time() < dl):
  edge = out["classes"] % 3 == 2
  x0 = edge_leaf(rg, ("D", "Xm", "diff")[out["classes"] // 3 % 3]) if edge else rg.randrange(1, p)
  v = exact_key(x0)[0]; out["classes"] += 1; out["edge_classes"] += edge; lab = label(v)
  for m in range(3):
   x = v * pow(BETA, (3 - m) % 3, p) % p; t, a, b = rg.randrange(3), rg.randrange(1, p), rg.randrange(1, p)
   w = leaf_w(represent(x, t, a, b, rg, out["tfp_terms"])[0], t)
   for k, c in enumerate((w[0], w[1], ~(w[0] + w[1]) % M128)):
    e = (t128(pow(BETA, k, p) * x % p) - c + 2 ** 127) % M128 - 2 ** 127
    out["errors"][k] = [min(out["errors"][k][0], e), max(out["errors"][k][1], e)]
   f = record_fields(w)
   if in_W(x): out["in_W"] += 1; continue
   if f[0] != m: out["wrong_index"] += 1; continue
   out["members"][m] += 1
   if in_D(x): out["in_D"] += 1
   elif f[1:] != lab: out["mismatches"] += 1
 tm = out["tfp_terms"]
 out["ok"] = (out["mismatches"] == 0 and out["wrong_index"] == 0 and out["in_D"] > 0 and out["in_W"] > 0
     and all(lo <= e[0] and e[1] <= hi for e, (lo, hi) in zip(out["errors"], ERR))
     and tm["sum_ok"] and all(tm[k] <= v for k, v in TERMS.items()))
 return out
def part_a_identities(rng):
 res = {}
 res["beta_cube_root_of_1_mod_p"] = pow(BETA, 3, p) == 1 and BETA != 1
 res["lambda_cube_root_of_1_mod_r"] = pow(LAM, 3, r) == 1 and LAM != 1
 res["lambda_G_equals_beta_x"] = PROD.mul(LAM, GP) == (BETA * GX % p, GY)
 res["seven_nonresidue_no_x0_point"] = pow(7, (p - 1) // 2, p) == p - 1
 res["window_forms_nonzero"] = all(c % p for c in (1 - BETA, 1 - BETA * BETA, BETA - BETA * BETA, 3))
 res["label_constants"] = (FA * FA - FA * FB + FB * FB == p and FB * BETA % p == FA - FB and (FA - FB) * BETA % p == p - FA)
 bad = drops = wrong = done = 0; kinds = {"random": 0, "near_tie": 0}; dl = deadline(0.10)
 while done < 26000 and (done < MIN_RANDOM or time.time() < dl):
  x, kind = (rng.randrange(1, p), "random") if done < 20000 else (near_tie_x(rng), "near_tie")
  kinds[kind] += 1; done += 1; win, err = index_check(x, rng)
  drops += win; wrong += win and err; bad += err and not win
 extra = max(0, MIN_NEAR - kinds["near_tie"]); rng2 = random.Random(SEED + 7919)
 for _ in range(extra):
  win, err = index_check(near_tie_x(rng2), rng2); kinds["near_tie"] += 1; drops += win; wrong += win and err; bad += err and not win
 res["class_rule"] = dict(checked=kinds, in_window_set=drops, wrong_index_in_window_set=wrong, mismatches_outside_window_set=bad, near_ties_from_second_stream=extra)
 res["class_labels"] = class_labels(random.Random(SEED + 104729))
 return res
SENTINEL = 0x5E5E5E5E5E5E5E5E
TIB = 1 << 40
A_SP = (34 << 40) + 0x7000000
EB_BASE = (1 << 45) + (1 << 40)
EB_END = EB_BASE + (1 << 34)
YR0, YR_STRIDE, XT0, BT0 = (1 << 45) + (1 << 31), 1 << 31, 59 * TIB, 0xfff8 << 48
DEL0 = 1088 * TIB
AXYX0 = 1280 * TIB
E0 = EB_END
SLOTW = (0, 1, 3, 2, 4, 5, 3)
def offs7(v): return [(k << 43) + (((v >> sh) & m) << 6) for k, (sh, m) in enumerate(zip(W7_SHIFTS, W7_MASKS))]
def x_consts(a): return 2 * a ** 3 % p, 2 * a % p, a * a % p, pow(a, 4, p)
def y_words(b): return b, pow(b, 3, p), 3 * b * b % p
def bt_addr(w, i):
 t = w >> 43 * i & (1 << 43) - 1
 return BT0 + (t << 8) + 32 * SLOTW.index(i) + (128 * (w >> 128 & 1) if i == 3 else 0)
def blook(w): return lkfold(sum(BETA * (w >> 43 * i & (1 << 43) - 1) * pow(2, 43 * i, p) % p for i in range(6)))
def lookup7(cval, w, add0=0):
 return lkfold(sum(((cval * ((w >> sh) & m) * pow(2, sh, p) + (add0 if k == 0 else 0)) % p) for k, (sh, m) in enumerate(zip(W7_SHIFTS, W7_MASKS))))
class LazyTable(dict):
 def __init__(self, tables, xpts=()):
  super().__init__()
  self.tables = tables
  self.pair49 = {}
  self.del32 = {}
  self.axyx = {}
  for j, (a, yx) in enumerate(xpts):
   AX = yx * yx % p
   self.axyx[AXYX0 + (j << 46)] = (2 * AX % p, 2 * yx % p)
  for pr in range(len(xpts) // 2):
   (a1, y1), (a2, y2) = xpts[2 * pr], xpts[2 * pr + 1]
   c10, c11, c12, a14 = x_consts(a1); c20, c21, c22, a24 = x_consts(a2)
   self.pair49[XT0 + (3 * pr + 0) * 49 * TIB] = (c10, 0, c20, 0)
   self.pair49[XT0 + (3 * pr + 1) * 49 * TIB] = (c11, a14, c21, a24)
   self.pair49[XT0 + (3 * pr + 2) * 49 * TIB] = (c12, 0, c22, 0)
   self.del32[DEL0 + (pr << 45)] = ((y2 * y2 - y1 * y1) % p, 0)
  for q in range(len(xpts) // 4):
   A1, A2 = xpts[4 * q][1] ** 2 % p, xpts[4 * q + 1][1] ** 2 % p
   A3, A4 = xpts[4 * q + 2][1] ** 2 % p, xpts[4 * q + 3][1] ** 2 % p
   sA, pA = (A1 + A2) % p, (A1 * A2) % p; sB, pB = (A3 + A4) % p, (A3 * A4) % p
   self.pair49[XT0 + (18 + q) * 49 * TIB] = ((-sB) % p, (pB - 7 * sB) % p, (-sA) % p, (pA - 7 * sA) % p)
 def __missing__(self, a):
  POISON = 0xDEADBEEFCAFEF00D
  if a >= BT0:
   o = a - BT0; k, t, w = (o >> 5) & 7, o >> 8, (o >> 3) & 3
   v = POISON if k == 7 or (k == 5 and t >= (1 << 41)) else (BETA * t * pow(2, 43 * SLOTW[k], p) % p) >> (64 * w) & M64
   self[a] = v; return v
  if AXYX0 <= a < (1 << 51):
   ab = a & ~((1 << 46) - 1)
   if ab in self.axyx:
    c_ax, c_yx = self.axyx[ab]; o = a - ab; t, k, half, w = o >> 9, (o >> 6) & 7, (o >> 5) & 1, (o >> 3) & 3
    v = POISON if k >= 7 or t > W7_MASKS[k] else (((c_ax if half == 0 else c_yx) * t * pow(2, W7_SHIFTS[k], p) % p) >> (64 * w) & M64)
    self[a] = v; return v
  if DEL0 <= a < AXYX0:
   tb = a & ~((1 << 45) - 1)
   if tb in self.del32:
    cval, add0 = self.del32[tb]; o = a - tb; t, k, w = o >> 8, (o >> 5) & 7, (o >> 3) & 3
    v = POISON if k >= 7 or t > W7_MASKS[k] else (((cval * t * pow(2, W7_SHIFTS[k], p) + (add0 if k == 0 else 0)) % p) >> (64 * w) & M64)
    self[a] = v; return v
  if XT0 <= a < DEL0:
   idx = (a - XT0) // (49 * TIB); tb = XT0 + idx * 49 * TIB
   if tb in self.pair49:
    c1, ad1, c2, ad2 = self.pair49[tb]; o = a - tb
    k, rem = o >> 43, o & ((1 << 43) - 1); t, half, w = rem >> 6, (rem >> 5) & 1, (rem >> 3) & 3
    cv, ad = (c1, ad1) if half == 0 else (c2, ad2)
    v = POISON if k >= 7 or t > W7_MASKS[k] else (((cv * t * pow(2, W7_SHIFTS[k], p) + (ad if k == 0 else 0)) % p) >> (64 * w)) & M64
    self[a] = v; return v
  for lo, hi, first, stride in self.tables:
   if lo <= a < hi: return first + (a - lo) // 8 * stride
  raise KeyError(hex(a))
def x_words(mem, e, Xs, l):
 j = l % MX; pr = j >> 1; X = Xs[j]
 for k, w in enumerate(corr(X[0])): put(mem, e + 256 + 16 * k, w, 2)
 put(mem, e + 288, X[1] * X[1] % p)
 mem[e + 320] = AXYX0 + (j << 46)
 if (j & 1) == 0:
  for k in range(3): mem[e + 328 + 8 * k] = XT0 + (3 * pr + k) * 49 * TIB
 else:
  mem[e + 328] = DEL0 + (pr << 45)
 if j == 0:
  mem[e + 384] = XT0 + 18 * 49 * TIB; mem[e + 392] = XT0 + 19 * 49 * TIB
 elif j == 8:
  mem[e + 384] = XT0 + 20 * 49 * TIB
def put_yrec(mem, e, Y):
 b, y = Y; B = y * y % p
 mem[e] = 0; put(mem, e + YR_B4, pow(b, 4, p)); put(mem, e + YR_B, B); put(mem, e + YR_B2, B * B % p)
 for k in range(4): put(mem, e + YR_YROW + 32 * k, y * pow(2, 64 * k, p) % p)
 put(mem, e + YR_NT, neg_t(b), 2); put(mem, e + YR_NT + 16, neg_t(BETA * b), 2)
 ob, ob3, o3b2 = offs7(b), offs7(pow(b, 3, p)), offs7(3 * b * b % p)
 for k in range(4): mem[e + 96 + 8 * k] = ob3[k]
 for k in range(3): mem[e + 352 + 8 * k] = ob3[4 + k]
 mem[e + 376] = 0
 for k in range(7): mem[e + 384 + 8 * k] = ob[k]
 mem[e + 440] = o3b2[0]
 for k in range(6): mem[e + 448 + 8 * k] = o3b2[1 + k]
 mem[e + 496] = mem[e + 504] = 0
def bwd_model(Xs, Ys, mem, acc_pairs, e_pairs, R, tm=None, ybase=None):
 if ybase is None: ybase = YR0
 out, prod = {}, lambda m, rr: sum((m >> 64 * k & M64) * rr[k] for k in range(4))
 for jp in range(len(acc_pairs) - 1, -1, -1):
  l1, l2 = 2 * jp, 2 * jp + 1
  (a1, yx1), (a2, yx2) = Xs[l1 % MX], Xs[l2 % MX]
  b, yy = Ys[l2 // MX]
  B, A1, A2 = yy * yy % p, yx1 * yx1 % p, yx2 * yx2 % p
  RR = shx_rows(R)
  inv_e, R = fold(prod(acc_pairs[jp], RR)), fold(prod(e_pairs[jp], RR))
  D2 = fs_op(B, A2)
  inv1 = fold(prod(inv_e, shx_rows(D2)))
  inv2 = lkfold(inv1 + sum((((A2 - A1) * ((inv_e >> sh) & m) * pow(2, sh, p)) % p) for sh, m in zip(W7_SHIFTS, W7_MASKS)))
  for l, (a, yx), inv in ((l2, (a2, yx2), inv2), (l1, (a1, yx1), inv1)):
   b_l, yy_l = Ys[l // MX]; AX = yx * yx % p
   yr = ybase + 512 * (l // MX)
   yrows = [get(mem, yr + YR_YROW + 32 * k) for k in range(4)]
   inv_sq = fold(inv * inv % 2 ** 256 + (inv * inv >> 256) * C0)
   Mval = lkfold(inv + sum((((2 * AX % p) * ((inv_sq >> sh) & m) * pow(2, sh, p)) % p) for sh, m in zip(W7_SHIFTS, W7_MASKS)))
   Pval = lookup7(2 * yx % p, inv_sq)
   Nval = fold(prod(Pval, yrows))
   w = [fs_op(Mval, Nval), faddw_op(Mval, Nval)]
   bw = [blook(v) for v in w]
   c0, c1, c2, a4 = x_consts(a); b4 = get(mem, yr + YR_B4)
   yw = y_words(b_l)
   cs = [lkfold(sum(((ci_v * ((wv >> sh) & m) * pow(2, sh, p) + (a4 if ci == 1 and k == 0 else 0)) % p) for k, (sh, m) in enumerate(zip(W7_SHIFTS, W7_MASKS))) + (b4 if ci == 0 else 0))
      for ci, (ci_v, wv) in enumerate(zip((c0, c1, c2), yw))]
   fr = {}
   for i, (c, yc) in enumerate(zip(cs, (b_l, a, None))):
    Rc, cw = shx_rows(c), corr(yc) if yc else (0, 0)
    for s in (0, 1):
     fr[s, 2 * i] = tfp(Rc, w[s], tm, True) if i == 2 else (tfp(Rc, w[s], tm) + cw[0]) % M128
     fr[s, 2 * i + 1] = tfp(Rc, bw[s], tm, True) if i == 2 else (tfp(Rc, bw[s], tm) + cw[1]) % M128
   out[l] = [fr[s, i] for s in (0, 1) for i in range(6)]
 return out, R
def six_leaves(X, Y):
 out = []
 for t in range(3):
  PX = (pow(BETA, t, p) * X[0] % p, X[1])
  out += [PROD.add(PX, Y)[0], PROD.add(PX, PROD.neg(Y))[0]]
 return out
def fixup_extra(prog):
 out, last = {}, -1
 for i, (mn, _) in enumerate(prog):
  if mn in ("b.lo", "b.hs", "cbnz"): last = i
  elif mn == "stadd": out[i] = i - last
 return out
def flags_dead_after_fixups(prog, labels):
 reads = ("adc", "adcs", "sbc", "sbcs", "csel", "cinc", "cset", "ccmp", "ccmn", "b.eq", "b.ne", "b.lo", "b.hs")
 def target(t, pc):
  if t[:-1].isdigit() and t[-1:] in ("f", "b"):
   c_ = labels[t[:-1]]
   return [q for q in c_ if q < pc][-1] if t[-1] == "b" else [q for q in c_ if q >= pc][0]
  return labels[t][0]
 for i in fixup_extra(prog):
  todo, seen = [i + 1], set()
  while todo:
   j = todo.pop()
   if j in seen: continue
   seen.add(j); mn, o = prog[j]
   if mn in reads: return False
   if mn in ("adds", "subs", "cmp", "cmn", "tst", "ret"): continue
   if mn in ("b", "tbz", "tbnz", "cbz", "cbnz"):
    todo.append(target(o[-1], j + 1))
    if mn == "b": continue
   todo.append(j + 1)
 return True
CBS, CB0, YSP = 2304, (34 << 40) + 0x2000000, (34 << 40) + 0x4000000
def chain_n(i):
 rho, q = divmod(i, 33)
 return 33 * rho + (0 if q == 32 else (16 - q // 2) * (1 - 2 * (q % 2)))
@functools.lru_cache(maxsize=None)
def chain_y(ch, i): return PROD.add(ch[0], PROD.mul(chain_n(i) % r, ch[1]))
def y_generator(fprog, flabels, rng, cert):
 mem, chains = {}, []
 for c in range(2):
  C, S = PROD.mul(rng.randrange(1, r), GP), PROD.mul(rng.randrange(1, r), GP); chains.append((C, S)); cb = CB0 + CBS * c
  for i, P in [(0, C)] + [(i, PROD.mul(i, S)) for i in list(range(1, 17)) + [33]]:
   put(mem, cb + (1088 if i == 33 else 64 * i), P[0]); put(mem, cb + (1088 if i == 33 else 64 * i) + 32, P[1])
  mem[cb + 2240] = YR0 + c * YR_STRIDE
 put(mem, CB0 - 128, 1)
 for a in range(YSP - 4096, YSP + 1024, 8): mem[a] = 0
 put(mem, YSP + 352, BETA); mem[YSP + 384], mem[YSP + 424] = C0, CB0
 counts = []
 for nch in (2, 1):
  mem[YSP + 416] = CB0 + CBS * nch
  counts.append(a64.a64_run(fprog, flabels, "yround", dict(sp=YSP, x30=SENTINEL), mem, limit=10 ** 7))
 bad = nrec = 0
 for c, rd in ((0, 0), (1, 0), (0, 1)):
  for q in range(33):
   i = 33 * rd + q; (b, y), e = chain_y(chains[c], i), YR0 + c * YR_STRIDE + 512 * i; B = y * y % p
   ok = (get(mem, e + YR_B4) % p == pow(b, 4, p) and get(mem, e + YR_B) == B and get(mem, e + YR_B2) == B * B % p)
   ok &= all(get(mem, e + YR_YROW + 32 * k) == y * pow(2, 64 * k, p) % p for k in range(4))
   ok &= get(mem, e + YR_NT, 2) == neg_t(b) and get(mem, e + YR_NT + 16, 2) == neg_t(BETA * b)
   ob, ob3, o3b2 = offs7(b), offs7(pow(b, 3, p)), offs7(3 * b * b % p)
   ok &= ([mem[e + 96 + 8 * k] for k in range(4)] + [mem[e + 352 + 8 * k] for k in range(3)]) == ob3
   ok &= [mem[e + 384 + 8 * k] for k in range(7)] == ob
   ok &= ([mem[e + 440]] + [mem[e + 448 + 8 * k] for k in range(6)]) == o3b2
   bad += not ok; nrec += 1
 cen = [(get(mem, CB0 + CBS * c), get(mem, CB0 + CBS * c + 32)) for c in range(2)]
 ok = bad == 0 and cen == [chain_y(chains[0], 33 * 2 + 32), chain_y(chains[1], 33 + 32)]
 cok = counts == [cert.YR_FIX + 2 * cert.YR_BLK, cert.YR_FIX + cert.YR_BLK]
 return dict(records=nrec, mismatches=bad, centres_ok=ok, counts=counts, counts_ok=cok), mem, chains[0], chains[1]
def leaf_batch(prog, labels, Xs, Ys, q2, ymem, cert, W2=0x123456789, KN=5 << 30, override=None, chain=None, rho=None, ybase=YR0):
 assert q2 % MX == 0 and ybase % YR_STRIDE == 0
 fx = fixup_extra(prog); vis = {k: dict.fromkeys(fx, 0) for k in ("fwd", "bwd")}
 nY, SP = q2 // MX, A_SP
 sb = EB_END - ENTRY * q2
 LS = ([(XT0, 1 << 51), (BT0, 1 << 64), (ybase, ybase + 512 * nY), (sb, EB_END)], [(sb, EB_END)])
 mem = LazyTable(machine.LC1_TABLES, Xs)
 if ymem is not None:
  for a in range(ybase - 64, ybase + 512 * nY, 8): mem[a] = ymem.get(a, 0)
 else:
  for a in range(ybase - 64, ybase, 8): mem[a] = 0
  for yi in range(nY): put_yrec(mem, ybase + 512 * yi, Ys[yi])
 for l in range(q2): x_words(mem, sb + ENTRY * l, Xs, l)
 rho = rho or random.Random(q2).randrange(1, p)
 regs = dict(x0=q2, x1=ybase + 64, x3=sb, sp=SP, x30=SENTINEL); regs.update(as_regs(rho, ("x24", "x25", "x26", "x27")))
 out = {}; nf = a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out, visits=vis["fwd"], ls64=LS)
 words, acc = [out["x0"]], reg_int(out, ("x24", "x25", "x26", "x27"))
 accs = [get(mem, sb + ENTRY * (2 * j) + ACC) for j in range(q2 // 2)]
 es = [get(mem, sb + ENTRY * (2 * j)) for j in range(q2 // 2)]
 R0 = pow(acc, p - 2, p); top = ybase + 512 * (nY - 1)
 regs = dict(x3=EB_END, x1=ybase - 32, sp=top, x30=SENTINEL); regs.update(as_regs(R0, ("x0", "x2", "x14", "x15")))
 out = {}; nb = a64.a64_run(prog, labels, "bwd", regs, mem, out_regs=out, visits=vis["bwd"], ls64=LS)
 words.append(out["x5"])
 want = rho
 for l in range(q2): want = want * (Ys[l // MX][1] ** 2 - Xs[l % MX][1] ** 2) % p
 final_R = reg_int(out, ("x0", "x2", "x14", "x15")); sp_kept = out["sp"] == top
 tm = {}; model, mR = bwd_model(Xs, Ys, mem, accs, es, R0, tm, ybase)
 bad_model, bad_bound, bad_leaf, err, leaves = int(mR != final_R), 0, 0, [[99, -99], [99, -99]], []
 for l in range(q2):
  X, Y = Xs[l % MX], Ys[l // MX]; e, cs = sb + ENTRY * l, components(X, Y)
  bad_model += sum(get(mem, e + o, 2) != w for o, w in zip(FRAC_OFF[0] + FRAC_OFF[1], model[l]))
  for i, d in enumerate(frac_errs(mem, e, X, Y)):
   g = i % 6 // 4; err[g] = [min(err[g][0], d), max(err[g][1], d)]
   bad_bound += not g - 1 <= d <= TERMS["total"]
  six = six_leaves(X, Y)
  bad_leaf += sum(dft_leaf(cs[3 * s:3 * s + 3], t, X[0], Y[0]) != six[2 * t + s] for s in (0, 1) for t in range(3))
  leaves.append(six)
 for (l, s), (x, t, rg) in (override or {}).items():
  f6, c3 = represent(x, t, Xs[l % MX][0], Ys[l // MX][0], rg)
  for i in (4, 5): f6[i] = f6[i] >> 64 << 64 | (f6[i] & M64) >> 3
  for o, v in zip(FRAC_OFF[s], f6): put(mem, sb + ENTRY * l + o, v, 2)
  for t_ in range(3): leaves[l][2 * t_ + s] = (c3[0] + BT[t_] * c3[1] + BT[2 * t_ % 3] * c3[2]) % p
 port = machine.Port1(mem)
 regs = dict(x2=W2, x3=sb + KEYS0, x26=KN, sp=SP, x30=SENTINEL)
 out = {}; nk = a64.a64_run(prog, labels, "keys", regs, mem, out_regs=out, dev=port.dev, ls64=([(sb, EB_END)], []))
 sp_kept &= out["sp"] == SP
 snd = machine.Sender1(cert)
 meter, fields = machine.keys_reference(sys.modules[__name__], mem, sb, q2, W2, KN, snd, cert)
 flat = [leaves[l][2 * t + s] for l in range(q2) for s in (0, 1) for t in range(3)]
 win = [in_W(x) for x in flat]; dset = [in_D(x) for x in flat]
 wrong = sum(not w and f[0] != exact_key(x)[1] for w, f, x in zip(win, fields, flat))
 canon = sum(f[1:] != class_label(x) for w, d, f, x in zip(win, dset, fields, flat) if not w and not d)
 dec = 0
 for k, (d, rec) in enumerate(snd.seq if chain else ()):
  v = (rec[3] >> 4) & ((1 << 26) - 1); s = v & 1; l = v >> 1
  Y, X = chain_y(chain, l // MX), Xs[l % MX]
  cand = [class_label(PROD.add((pow(BETA, t, p) * X[0] % p, X[1]), PROD.neg(Y) if s else Y)[0]) for t in range(3)]
  dec += (l, s) != (k // 6, k // 3 % 2) or Y != Ys[l // MX] or \
   (not win[k] and not dset[k] and (d, rec[0] & M42, rec[0], rec[1]) not in cand)
 nfix = {k: [sum(v for i, v in vis[k].items() if fx[i] == e) for e in (3, 4, 5)] for k in vis}
 ext = {k: 3 * nfix[k][0] + 4 * nfix[k][1] + 5 * nfix[k][2] for k in vis}
 return dict(fwd=nf, bwd=nb, keys=nk, acc_ok=acc % p == want, fraction_mismatches=bad_model,
    fractions_outside_bounds=bad_bound, fraction_errors=err, tfp_terms=tm, identity_mismatches=bad_leaf,
    fixups_3_4=nfix, fixup_extra=ext, record_decode_mismatches=dec,
    fixup_words_ok=words == [C0 * sum(nfix[k]) for k in ("fwd", "bwd")],
    fwd_bwd_formulas=nf == cert.fwd_count(q2) + ext["fwd"] and nb == cert.bwd_count(q2) + ext["bwd"],
    records_ok=port.blocks == snd.msgs and machine.state1_ok(mem, snd),
    records=len(snd.seq), key_paths=meter, in_W=sum(win),
    in_D=sum(d and not w for w, d in zip(win, dset)), wrong_index=wrong, keys_not_canonical=canon,
    keys_formula=nk == cert.KEY_FIX + cert.KEY_INC * q2 + meter + snd.extra,
    final_R_times_rho_is_one=final_R * rho % p == 1, sp_kept=sp_kept)
def macro_counts(src, entry, names):
 macros, cur, lines = {}, None, [ln.split("//")[0].strip() for ln in src.split("\n")]
 for s in lines:
  if s.startswith(".macro"): cur = s.split()[1].rstrip(","); macros[cur] = []
  elif s.startswith(".endm"): cur = None
  elif cur is not None: macros[cur].append(s)
 out = dict.fromkeys(names, 0)
 def walk(seq):
  for s in seq:
   if not s: continue
   w = s.split()[0].rstrip(",")
   if w in out: out[w] += 1
   if w in macros: walk(macros[w])
 start = lines.index(entry + ":")
 walk(lines[start + 1:lines.index("ret", start)])
 return out
def sqrt_mod_p(v):
 s_ = pow(v, (p + 1) // 4, p)
 return s_ if s_ * s_ % p == v % p else None
def cube_root(v):
 c_ = pow(v, (p + 2) // 9, p)
 return c_ if c_ * c_ * c_ % p == v % p else None
def degenerate_batch(prog, labels, Xs, Ys):
 sb = EB_END - ENTRY * MX
 LS = ([(XT0, 1 << 51), (BT0, 1 << 64), (YR0, YR0 + 512), (sb, EB_END)], [(sb, EB_END)])
 mem = LazyTable(machine.LC1_TABLES, Xs)
 for a in range(YR0 - 64, YR0, 8): mem[a] = 0
 put_yrec(mem, YR0, Ys[0])
 for l in range(MX): x_words(mem, sb + ENTRY * l, Xs, l)
 regs = dict(x0=MX, x1=YR0 + 64, x3=sb, x24=7, x25=0, x26=0, x27=0, sp=A_SP, x30=SENTINEL)
 out = {}; a64.a64_run(prog, labels, "fwd", regs, mem, out_regs=out, ls64=LS)
 acc = reg_int(out, ("x24", "x25", "x26", "x27"))
 regs = dict(x3=EB_END, x1=YR0 - 32, sp=YR0, x30=SENTINEL)
 regs.update(as_regs(pow(acc, p - 2, p), ("x0", "x2", "x14", "x15")))
 a64.a64_run(prog, labels, "bwd", regs, mem, ls64=LS)
 good = 0
 for l, X in enumerate(Xs):
  if (Ys[0][1] ** 2 - X[1] ** 2) % p:
   good += sum(-1 <= d <= TERMS["total"] for d in frac_errs(mem, sb + ENTRY * l, X, Ys[0]))
 return acc, good
def edge_batch(prog, labels, Xs, Ys, ymem, cert, rg, kinds, need):
 over = {}
 for l in range(MX):
  for s in (0, 1):
   while True:
    x = edge_leaf(rg, kinds[(2 * l + s) % len(kinds)])
    if len(over) >= 4 or (in_W(x) if need == "W" else in_D(x) and not in_W(x)): break
   over[(l, s)] = (x, rg.randrange(3), rg)
 return leaf_batch(prog, labels, Xs, Ys, MX, ymem, cert, override=over)
def tbld_production_tables(leaf_src, rng):
 rsrc = open(os.path.join(HERE, "kernels", "route.s")).read()
 rprog, rlab = a64.a64_load(rsrc)
 mrp, mrl = a64.a64_load(rsrc.replace("\ttst x6,#1\n\tcsel x18,x3,x2,ne\n\tlsr x6,x6,#1", "\tcmp x8,#2\n\tcsel x18,x3,x2,eq"))
 hprog, hlab = a64.a64_load(leaf_src + sites.SITE_KERNELS + """
t_lkp:
\tENTER
\tldp x6,x7,[x1,#288]
\tldr x30,[x1,#304]
\tLKP7_PAIR 384
\tLEAVE
t_c1b:
\tENTER
\tC1BSPLIT7 x25,x26,x27,x28
\tLEAVE
t_efg:
\tENTER
\tEFG_PAIR_6
\tLEAVE
""")
 PSN, DSC, sb = 0xDEADBEEFCAFEF00D, 0x300000, EB_END - ENTRY * MX
 LS_H = ([(XT0, 1 << 51), (BT0, 1 << 64), (YR0, YR0 + 512), (sb, EB_END)], [(sb, EB_END)])
 wval = lambda ts: sum(t << sh for t, sh in zip(ts, W7_SHIFTS))
 def psn(m, b, nw):
  for i in range(nw): m[b + 8 * i] = PSN
 def bp49(m, tb, c1, a1, c2, a2, rp_=rprog, rl_=rlab):
  tot = 0
  for k, sh in enumerate(W7_SHIFTS):
   psn(m, tb + (k << 43), 64); dk = DSC + 128 * k
   put(m, dk, c1 * pow(2, sh, p) % p); put(m, dk + 32, a1 if k == 0 else 0)
   put(m, dk + 64, c2 * pow(2, sh, p) % p); put(m, dk + 96, a2 if k == 0 else 0)
   tot += a64.a64_run(rp_, rl_, "tbld", dict(x0=tb + (k << 43), x1=dk, x2=8, x3=2, x4=2, x5=64, x6=3 if k == 6 else 0, x30=SENTINEL), m)
  return tot
 def run_h(ent, regs, m):
  o = {}; a64.a64_run(hprog, hlab, ent, dict(regs, sp=YR0, x30=SENTINEL), m, out_regs=o, ls64=LS_H); return o
 A1, A2, Bv = rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p)
 cd, inv_e, ep = (A2 - A1) % p, wval([3, 5, 2, 6, 4, 7, 1]), sb + 576
 def try_del(rp_, rl_):
  m = {}; psn(m, DEL0, 256)
  for k, sh in enumerate(W7_SHIFTS): put(m, DSC + 64 * k, cd * pow(2, sh, p) % p); put(m, DSC + 64 * k + 32, 0)
  n = a64.a64_run(rp_, rl_, "tbld", dict(x0=DEL0, x1=DSC, x2=8, x3=2, x4=7, x5=256, x6=1 << 6, x30=SENTINEL), m)
  put(m, YR0 + 288, Bv); put(m, ep + 224, A2); m[ep + 264] = DEL0
  o = run_h("t_c1b", dict(x3=ep, **as_regs(inv_e, ("x25", "x26", "x27", "x28"))), m)
  return n, get(m, ep - 64) % p, reg_int(o, ("x8", "x9", "x10", "x11")) % p
 nd, i1, i2 = try_del(rprog, rlab); _, _, i2_mut = try_del(mrp, mrl)
 w1 = inv_e * (Bv - A2) % p; w2 = (w1 + cd * inv_e) % p
 del_ok, mut_del = (nd == 4 + 7 * 11 + 20 * 50 and (i1, i2) == (w1, w2)), (i2_mut != w2)
 ma, cax, cyx, yYv = {}, rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p)
 inv_p1 = (1 << 111) + (2 << 74) + (1 << 37) + 1; isq = inv_p1 * inv_p1; mam = {}; psn(ma, AXYX0, 512); psn(mam, AXYX0, 512)
 for k, sh in enumerate(W7_SHIFTS):
  put(ma, DSC + 128 * k, cax * pow(2, sh, p) % p); put(ma, DSC + 128 * k + 32, 0)
  put(ma, DSC + 128 * k + 64, cyx * pow(2, sh, p) % p); put(ma, DSC + 128 * k + 96, 0)
 for m_ in (ma, mam):
  for a_, v_ in list(ma.items()):
   if DSC <= a_ < DSC + 896: m_[a_] = v_
 na = a64.a64_run(rprog, rlab, "tbld", dict(x0=AXYX0, x1=DSC, x2=8, x3=2, x4=14, x5=512, x6=3 << 12, x30=SENTINEL), ma)
 a64.a64_run(mrp, mrl, "tbld", dict(x0=AXYX0, x1=DSC, x2=8, x3=2, x4=14, x5=512, x6=3 << 12, x30=SENTINEL), mam)
 for m_ in (ma, mam):
  m_[ep + 256] = AXYX0
  for k in range(4): put(m_, YR0 + YR_YROW + 32 * k, yYv * pow(2, 64 * k, p) % p)
 oa = run_h("t_p1", dict(x3=ep, **as_regs(inv_p1, ("x8", "x9", "x10", "x11"))), ma)
 Mv, Nv = (inv_p1 + cax * isq) % p, cyx * isq % p * yYv % p
 u2_got = oa["x28"] | oa["x29"] << 64 | ma[YR0 + 496] << 128 | oa["x1"] << 192
 axyx_ok = (na == 4 + 14 * 11 + 20 * 100 and u2_got % p == (Mv - Nv) % p
  and reg_int(oa, ("x0", "x2", "x5", "x6")) % p == (Mv + Nv) % p
  and ma[AXYX0 + 2 * 512 + 13 * 32] == PSN and mam[AXYX0 + 2 * 512 + 13 * 32] != PSN)
 mp, c1s, ad1, c2s, ad2, B2v, wb3 = {}, rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p), wval([5, 3, 6, 2, 7, 7, 1])
 np49, ob3 = bp49(mp, XT0, c1s, ad1, c2s, ad2), offs7(wb3)
 put(mp, YR0 + 64, B2v); put(mp, YR0 + 96, B2v); mp[sb + 384] = XT0
 for k in range(4): mp[YR0 + 96 + 8 * k] = ob3[k]
 for k in range(3): mp[YR0 + 352 + 8 * k] = ob3[4 + k]
 ol = run_h("t_lkp", dict(x1=YR0 + 64, x3=sb), mp)
 mpm = dict(mp); bp49(mpm, XT0, c1s, ad1, c2s, ad2, mrp, mrl); olm = run_h("t_lkp", dict(x1=YR0 + 64, x3=sb), mpm)
 p49_ok = (np49 == 7 * (4 + 2 * 11) + 20 * 100 and reg_int(ol, ("x16", "x17", "x18", "x19")) % p == (c1s * wb3 + ad1 + B2v) % p
  and reg_int(ol, ("x20", "x21", "x22", "x23")) % p == (c2s * wb3 + ad2 + B2v) % p
  and reg_int(olm, ("x16", "x17", "x18", "x19")) % p != (c1s * wb3 + ad1 + B2v) % p)
 me, wb, w3b2, b4v = {}, wval([2, 6, 4, 5, 3, 7, 1]), wval([3, 4, 7, 2, 6, 6, 1]), rng.randrange(1, p)
 put(me, YR0 + 256, b4v); put(me, YR0 + 288, 0)
 for k, v in enumerate(offs7(wb3)[:4]): me[YR0 + 96 + 8 * k] = v
 for k, v in enumerate(offs7(wb3)[4:]): me[YR0 + 352 + 8 * k] = v
 for k, v in enumerate(offs7(wb)): me[YR0 + 384 + 8 * k] = v
 me[YR0 + 440] = offs7(w3b2)[0]
 for k, v in enumerate(offs7(w3b2)[1:]): me[YR0 + 448 + 8 * k] = v
 ec = [[(rng.randrange(1, p), rng.randrange(1, p) if mi == 1 else 0, rng.randrange(1, p), rng.randrange(1, p) if mi == 1 else 0) for mi in range(3)] for _ in range(6)]
 for pr in range(6):
  for mi, cs_ in enumerate(ec[pr]):
   tb = XT0 + (3 * pr + mi) * 49 * TIB; bp49(me, tb, *cs_); me[sb + ENTRY * (2 * pr) + 328 + 8 * mi] = tb
 run_h("t_efg", dict(x3=EB_END - 448, x0=1, x2=2, x14=3, x15=4), me)
 efg_ok = all(get(me, sb + ENTRY * (2 * pr) + 64 + 64 * mi) % p == (ec[pr][mi][0] * wi + ec[pr][mi][1] + (b4v if mi == 0 else 0)) % p
  and get(me, sb + ENTRY * (2 * pr) + 96 + 64 * mi) % p == (ec[pr][mi][2] * wi + ec[pr][mi][3] + (b4v if mi == 0 else 0)) % p
  for pr in range(6) for mi, wi in enumerate((wb, wb3, w3b2)))
 mb, mbm = {}, {}; psn(mb, BT0, 256); psn(mbm, BT0, 256)
 for idx, wi in enumerate(SLOTW):
  for m_ in (mb, mbm): put(m_, DSC + 64 * idx, BETA * pow(2, 43 * wi, p) % p); put(m_, DSC + 64 * idx + 32, 0)
 nbt = a64.a64_run(rprog, rlab, "tbld", dict(x0=BT0, x1=DSC, x2=8, x3=2, x4=7, x5=256, x6=1 << 5, x30=SENTINEL), mb)
 a64.a64_run(rprog, rlab, "tbld", dict(x0=BT0, x1=DSC, x2=8, x3=2, x4=7, x5=256, x6=1 << 4, x30=SENTINEL), mbm)
 u2 = sum(t << (43 * i) for i, t in enumerate([6, 5, 7, 5, 4, 1]))
 v2 = sum(t << (43 * i) for i, t in enumerate([4, 7, 6, 6, 7, 1]))
 ob = run_h("t_blk", dict(**as_regs(u2, ("x28", "x29", "x11", "x1")), **as_regs(v2, ("x0", "x2", "x5", "x6"))), mb)
 obm = run_h("t_blk", dict(**as_regs(u2, ("x28", "x29", "x11", "x1")), **as_regs(v2, ("x0", "x2", "x5", "x6"))), mbm)
 bt_ok = (nbt == 4 + 7 * 11 + 20 * 50 and reg_int(ob, ("x24", "x25", "x26", "x27")) % p == BETA * u2 % p
  and reg_int(ob, ("x20", "x21", "x22", "x23")) % p == BETA * v2 % p
  and reg_int(obm, ("x24", "x25", "x26", "x27")) % p != BETA * u2 % p)
 ok = del_ok and mut_del and axyx_ok and p49_ok and efg_ok and bt_ok
 return dict(ok=ok, del_ok=del_ok, mut_del_rejected=mut_del, axyx_ok=axyx_ok, pair49_ok=p49_ok and efg_ok, bt_ok=bt_ok)
def part_a_kernels(rng, cert):
 res = {}
 src = open(os.path.join(HERE, "kernels", "leaf.s")).read()
 prog, labels = a64.a64_load(src)
 fprog, flabels = a64.a64_load(open(os.path.join(HERE, "kernels", "secp_fe.s")).read())
 yg, ymem, ch0, ch1 = y_generator(fprog, flabels, rng, cert)
 res["y_generator"] = yg
 Xs = [PROD.mul(rng.randrange(1, r), GP) for _ in range(MX)]
 Ys = [chain_y(ch0, i) for i in range(12)]
 runs = {}
 for q2 in (12, 24, 36, 60, 96):
  runs[q2] = leaf_batch(prog, labels, Xs, Ys, q2, ymem, cert, chain=ch0)
 ov_top = {(0, 0): (planted_x(rng, (1 << 42) - 4, 0)[0], 0, rng), (0, 1): (planted_x(rng, (1 << 42) - 3, 1)[0], 1, rng)}
 c1_run = leaf_batch(prog, labels, Xs, [chain_y(ch1, i) for i in range(2)], 24, ymem, cert, override=ov_top, chain=ch1, ybase=YR0 + YR_STRIDE)
 res["yslot1_and_lc1_top_ok"] = (c1_run["records_ok"] and c1_run["keys_formula"] and c1_run["fwd_bwd_formulas"] and c1_run["sp_kept"])
 res["runs"] = runs
 net = {q: {k: d[k] - d["fixup_extra"][k] for k in ("fwd", "bwd")} for q, d in runs.items()}
 res["fwd_bwd_formulas_hold"] = all(d["fwd_bwd_formulas"] and d["fixup_words_ok"] for d in runs.values())
 res["counts_ok"] = all(net[q]["fwd"] == cert.fwd_count(q) and net[q]["bwd"] == cert.bwd_count(q) for q in runs)
 res["fixups"] = dict(sites=len(fixup_extra(prog)), flags_dead_at_every_join=flags_dead_after_fixups(prog, labels))
 tbzs = [o for mn, o in prog if mn in ("tbz", "tbnz")]
 res["loop_ends_on_bit_34"] = [o for o in tbzs if o[0] == "x3"] == [["x3", "#34", "1b"]] and all(o[1] == "#9" for o in tbzs if o[0] == "x5")
 res["keys_formula_holds"] = all(d["keys_formula"] for d in runs.values())
 res["tfp_terms"] = dict(largest={k: max(d["tfp_terms"][k] for d in runs.values()) for k in TERMS}, bounds=TERMS,
       sums_exact=all(d["tfp_terms"]["sum_ok"] for d in runs.values()))
 res["fraction_errors"] = [[min(d["fraction_errors"][g][0] for d in runs.values()), max(d["fraction_errors"][g][1] for d in runs.values())] for g in (0, 1)]
 res["all_correct"] = (yg["mismatches"] == 0 and yg["centres_ok"] and yg["counts_ok"] and res["yslot1_and_lc1_top_ok"] and all(
  d["acc_ok"] and d["fraction_mismatches"] == 0 and d["fractions_outside_bounds"] == 0 and d["identity_mismatches"] == 0
  and d["record_decode_mismatches"] == 0 and d["records_ok"] and d["keys_not_canonical"] == 0 and d["wrong_index"] == 0
  and d["final_R_times_rho_is_one"] and d["sp_kept"] for d in runs.values()))
 erng = random.Random(SEED + 49979687)
 for name, kinds, need in (("window_batch", ("Xm", "diff"), "W"), ("boundary_batch", ("D",), "D")):
  w = edge_batch(prog, labels, Xs, Ys, ymem, cert, erng, kinds, need)
  res[name] = dict(leaves=6 * MX, in_W=w["in_W"], in_D_outside_W=w["in_D"], wrong_index=w["wrong_index"],
       keys_not_canonical=w["keys_not_canonical"], records=w["records"], records_ok=w["records_ok"],
       ok=(w[{"W": "in_W", "D": "in_D"}[need]] >= 4 and w["wrong_index"] == 0 and w["keys_not_canonical"] == 0
        and w["records"] == 6 * MX and w["records_ok"] and w["keys_formula"] and w["fwd_bwd_formulas"] and w["fixup_words_ok"]))
 gp, gl = a64.a64_load(src.replace("\tb.lo 9f\n", "\tcbnz xzr,9f\n").replace("\tb.hs 9f\n", "\tcbnz xzr,9f\n"))
 g = leaf_batch(gp, gl, Xs, Ys, MX, ymem, cert, chain=ch0)
 tot_pair = sum(cert.FIX_SITES_PAIR.values())
 ext_pair = sum(cert.FIX_EXTRA[m] * cert.FIX_SITES_PAIR[m] for m in cert.FIX_SITES_PAIR)
 res["forced_fixups"] = dict(incidences=MX, fixups_3_4=g["fixups_3_4"], fixup_extra=g["fixup_extra"],
        fraction_mismatches=g["fraction_mismatches"],
        ok=(g["fraction_mismatches"] == 0 and g["fractions_outside_bounds"] == 0 and g["identity_mismatches"] == 0
         and g["fwd_bwd_formulas"] and g["fixup_words_ok"] and g["records_ok"]
         and sum(map(sum, g["fixups_3_4"].values())) == (MX // 2) * tot_pair and sum(g["fixup_extra"].values()) == (MX // 2) * ext_pair))
 cnts, ok = set(), True
 for _ in range(3):
  a = rng.randrange(1, p); mem = {}; put(mem, 0x1000, a)
  cnts.add(a64.a64_run(fprog, flabels, "fe_inv", dict(x0=0x2000, x1=0x1000, sp=0x100000, x30=SENTINEL), mem, limit=10 ** 6))
  ok &= get(mem, 0x2000) % p == pow(a, p - 2, p)
 res["fe_inv"] = dict(counts=sorted(cnts), correct=ok and cnts == {cert.INV})
 names = ("SHX", "SHXT", "FDX", "F5", "FS", "FADDW", "LKFOLD4")
 counts = {k: macro_counts(src, k, names + ("LKP7_PAIR", "BWDPAIR")) for k in ("fwd", "bwd", "keys")}
 copies_f, copies_b = 2 * counts["fwd"]["LKP7_PAIR"], counts["bwd"]["BWDPAIR"]
 per_pair = {}
 for n in names: per_pair[n] = counts["fwd"][n] // copies_f + counts["bwd"][n] // copies_b
 per_pair["SHX"] += per_pair.pop("SHXT")
 per_pair["BLK"] = per_pair.pop("LKFOLD4")
 res["site_counts"] = dict(per_x_pair=per_pair, pair_copies=dict(fwd=copies_f, bwd=copies_b))
 res["sites_match_certificate"] = (per_pair == cert.FIX_SITES_PAIR and copies_f == 6 and copies_b == 6 and sum(counts["keys"].values()) == 0)
 srng = random.Random(SEED + 15485863)
 res["site_checks"] = sites.site_checks(sys.modules[__name__], srng, fprog, flabels)
 res["tbld_production"] = tbld_production_tables(src, random.Random(SEED + 27182818))
 phY = (BETA * Ys[0][0] % p, Ys[0][1])
 acc0, good = degenerate_batch(prog, labels, [Xs[0], phY, Xs[1], PROD.neg(Ys[0])] + Xs[2:10], [Ys[0]])
 res["d_zero_batch"] = dict(acc_zero=acc0 % p == 0, fractions_within_bounds=good, all_other_fractions_wrong=good == 0)
 prog_r, _ = a64.a64_load(open(os.path.join(HERE, "kernels", "route.s")).read())
 tot_ins = len(prog) + len(prog_r) + len(fprog)
 res["program_instructions"] = dict(leaf=len(prog), route=len(prog_r), secp_fe=len(fprog), total=tot_ins, fits_2_17_bytes=4 * tot_ins < 2 ** 17)
 return res
def lattice_point(tv, tu):
 c1 = (2 * (FA * tv + (FA - FB) * tu) + p) // (2 * p)
 c2 = (2 * ((FA - FB) * tv - FB * tu) + p) // (2 * p)
 return c1 * FB + c2 * (FA - FB), c1 * (FA - FB) - c2 * FA
def planted_x(rg, d1, d2, v=None, j=None):
 while True:
  vc = v
  if vc is None:
   c1 = (rg.getrandbits(22) << 42 | d2) << 42 | d1
   c2 = (rg.getrandbits(17) << 47 | rg.randrange(16)) << 42 | rg.getrandbits(42)
   vv, uu = lattice_point((4 * c1 + 1) * p >> PART + 2, (4 * c2 + 1) * p >> PART + 2)
   if not (0 < vv < p and uu == BETA * vv % p and comp(vv) == c1 and comp(uu) == c2 and exact_key(vv)[0] == vv): continue
   vc = vv
  m = rg.randrange(3) if j is None else j
  x = vc * pow(BETA, (3 - m) % 3, p) % p
  if exact_key(x)[1] == m and not in_W(x) and not in_D(x): return x, vc
def load_certificate():
 buf = io.StringIO()
 with contextlib.redirect_stdout(buf): import certificate as cert
 lines = [" ".join(ln.split()) for ln in buf.getvalue().splitlines()]
 return cert, [ln for ln in lines if " ok" in ln or "FAIL" in ln or ln.startswith(("q1 ", "E_total", "declared:"))]
def part_d(cert):
 names = ("LEAVES", "FWD_INC", "FWD_FIX", "BWD_INC", "BWD_FIX", "KEY_INC", "KEY_FIX", "EMIT", "EMIT_INC",
    "JPATH", "SJOIN", "FIX_SITES_PAIR", "FIX_EXTRA", "DELTA_F", "INV", "TINIT_ENTRY", "K", "K2B", "FLUSH",
    "FLUSH1_X", "FIN1_X", "FLUSH2_X", "FLUSH2_D", "RCV", "RCV1", "FIN_CNT", "H2S_DONE", "WAKE", "EXPECT", "EXPECTH",
    "H2S_ARR", "H2S_FULL", "DLV", "DONE", "G_GROUPS", "HB", "HCHG_NEW", "DLVHB", "DLVHS", "HCHG_VISIT", "HCHG_VER",
    "VISIT", "HCHG_END", "ZW", "ZERO_PASS", "STOP_CHECK", "NBAR_ALLOW", "HOME_MAX", "GATE_INSNS", "K_GATE")
 return dict(passed=not cert.FAILED, failed=cert.FAILED, checks=CERT_LINES,
    claim="%d.%015d" % divmod(int(cert.CLAIM_E * 10 ** 15), 10 ** 15), q1=cert.Q1GD,
    per_leaf_declared=format(float(cert.DECLARED["per_leaf"]), ".4f"),
    mu=float(cert.TFULL["mu"]), opened_success_lower_bound=float(cert.LO_RUN), gate_prob=float(cert.P_GATE), success_lower_bound=float(cert.GATED_LO), per_leaf=float(cert.PER_LEAF_OPEN),
    groups=cert.G_GROUPS, counts={k: str(getattr(cert, k)) if not isinstance(getattr(cert, k), (int, tuple)) else getattr(cert, k) for k in names})
def main():
 global CERT_LINES
 rng = random.Random(SEED)
 cert, CERT_LINES = load_certificate()
 OUT["input"] = cfg
 OUT["budget"] = dict(seconds=BUDGET, reserve=RESERVE, minimums=dict(class_rule_random=MIN_RANDOM, class_rule_near_tie=MIN_NEAR, class_labels=MIN_CLASSES))
 A_id = OUT["A_identities"] = part_a_identities(rng)
 id_ok = (A_id["beta_cube_root_of_1_mod_p"] and A_id["lambda_cube_root_of_1_mod_r"] and A_id["lambda_G_equals_beta_x"]
    and A_id["seven_nonresidue_no_x0_point"] and A_id["window_forms_nonzero"] and A_id["label_constants"]
    and A_id["class_rule"]["mismatches_outside_window_set"] == 0 and A_id["class_labels"]["ok"])
 if DO_KERNELS:
  A_k = OUT["A_kernels"] = part_a_kernels(rng, cert)
  B_m = OUT["B_machine"] = machine.part_b(sys.modules[__name__], cert)
  OUT["counts_match_certificate"] = cm = (
   A_k["fwd_bwd_formulas_hold"] and A_k["counts_ok"] and A_k["fixups"]["flags_dead_at_every_join"]
   and A_k["loop_ends_on_bit_34"] and A_k["keys_formula_holds"] and A_k["tfp_terms"]["sums_exact"]
   and A_k["all_correct"] and A_k["window_batch"]["ok"] and A_k["boundary_batch"]["ok"]
   and A_k["forced_fixups"]["ok"] and A_k["fe_inv"]["correct"] and A_k["sites_match_certificate"]
   and A_k["site_checks"]["all_true"] and A_k["tbld_production"]["ok"] and A_k["d_zero_batch"]["acc_zero"] and A_k["d_zero_batch"]["all_other_fractions_wrong"]
   and A_k["program_instructions"]["fits_2_17_bytes"]
   and B_m["formulas_ok"] and B_m["block_paths"]["ok"] and B_m["grouped_paths"]["ok"]
   and B_m["ls64"]["ok"] and B_m["va52"]["ok"] and B_m["encoding"]["ok"] and B_m["phase4_round"]["ok"] and B_m["y_ring"]["ok"] and B_m["gate"]["ok"])
 else:
  cm = True
 D_c = OUT["D_certificate"] = part_d(cert)
 OUT["all_checks_pass"] = bool(id_ok and cm and D_c["passed"])
 OUT["seconds"] = round(time.time() - T0, 2)
 json.dump(OUT, sys.stdout, indent=1)
 print()
 if not OUT["all_checks_pass"]:
  sys.exit(1)
if __name__ == "__main__":
 main()
