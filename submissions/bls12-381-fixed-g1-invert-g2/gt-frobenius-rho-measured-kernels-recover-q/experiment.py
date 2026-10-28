"""Lazily evaluated (prefix-product trie) Frobenius-class rho in BLS12-381 GT: identities, the shipped kernels
run in the emulator, the production lazy walk driven by those kernels against an eager reference, toy-field
rho with trails and replay, the pair-collision law, the fruitless-cycle bound, the instance reduction and G2
output (reduction.py) and the exact ledger (parts A, S, K, C, G, B, E, F, R, D)."""
import importlib, json, math, os, random, runpy, sys, time
from fractions import Fraction

sys.dont_write_bytecode = True
T0 = time.time()
cfg = json.load(sys.stdin)
Q_TOY = cfg.get("q", 67)
TRIALS = cfg.get("trials", 20)
R = cfg.get("R", 32)               # toy census table size, the production R_PROD below unless overridden
BRANCH_BITS = cfg.get("branch_bits", 6)
SQ_CUT = cfg.get("sq_cut", 0)
DP_BITS = cfg.get("dp_bits", 3)
WALKS = cfg.get("walks", 4)
SEED = cfg.get("seed", 1)
PROD_CHECKS = cfg.get("production_checks", 2)
PROD_STEPS = cfg.get("production_steps", 240)
KERNEL_CASES = cfg.get("kernel_cases", 3)
PROD_DP_BITS = cfg.get("production_dp_bits", 3)
STAT_CLASSES = cfg.get("stat_classes", 3000)
PAIR_SAMPLES = cfg.get("pair_samples", 0)
CYCLE_M0 = 12
rng = random.Random(SEED)

HERE = os.path.dirname(os.path.abspath(__file__))

def locate(name):
  """A shipped file or directory beside this script (or under files/ for the grader's snapshot layout)."""
  for base in (HERE, os.path.join(HERE, "files"), os.getcwd(), os.path.join(os.getcwd(), "files")):
    cand = os.path.join(base, name)
    if os.path.exists(cand):
      return cand
  for dirpath, dirs, names in os.walk(HERE):
    if name in names or name in dirs:
      return os.path.join(dirpath, name)
  raise SystemExit("experiment.py: %s not found beside %s" % (name, HERE))

# Production parameters (proof.md sections 3 and 5).
R_PROD, L_PROD, THETA_BITS, STUCK_M = 32, 6, 86, 28
C_UNIT = 20000
LAZY_CORE, ROOT_CORE, MUL_CORE, SQ_CORE = 3457, 3074, 12735, 6786
EVERY_ITEMS = {
  "K reads and the multiplier index: ldr K0, mov #32, umulh (h = floor(32 K0 / 2^64))": 3,
  "distinguished-point test on bits 64..149 of K: ldp K1,K2; and x4, K2, #0x3fffff; orr x4, x4, K1; cbz": 4,
  "kappa load and twist t = (12 - kappa) mod 12: ldr, mov #12, sub, cmp, csel": 5,
  "child word index 384 nu + 1 + 12 h + t: two movs, two madd, add": 5,
  "depth increment, compare with L = 6, branch to the flush": 3,
  "node address base + 3,552 nu' (mov, madd) or leaf address (sub, madd, csel)": 5,
  "remember the previous node, base-buffer index and K slot for escapes (str, stp)": 2,
  "2-cycle and checkpoint comparisons on K, history ring rotation, checkpoint store test": 60,
  "loop control: back-branch 1, chunk counter subs/b.eq 2, root/lazy select 1, DP branch 1, 2-cycle and "
  "checkpoint branches 4, kernel call bl + 4 argument moves + caller frame link 3 (8)": 17,
  "margin for loads and moves": 53,
  "192-bit work meter: price in a register and adds/adcs/adc (4), exact 3-word compare cmp/sbcs/sbcs/b.hs (4)": 8,
  "2-word stuck-walk counter: adds/adc (2), mov x16,#0x7000000; cmp x9,x16; b.hs (3)": 5,
}
EVERY = sum(EVERY_ITEMS.values())
FLUSH_ITEMS = {"mul.s kernel (y <- y * Pi_leaf)": MUL_CORE,
       "call: bl, 4 argument moves, caller frame link": 8,
       "swap of the two base buffers, node and depth reset (3 moves)": 3,
       "margin": 7}
FLUSH = sum(FLUSH_ITEMS.values())
ESC_ITEMS = {"choose the smaller of the two cycle keys (6-limb compare and selects)": 14,
      "address of that class's node and base buffer": 4,
      "flush its representative: mul.s (skipped at the root) and its call": MUL_CORE + 8,
      "cyclotomic square sqr.s and its call": SQ_CORE + 8,
      "reset histories, checkpoint, node and depth": 20,
      "pre-debit compare and charge add": 8,
      "margin": 24000 - (14 + 4 + MUL_CORE + 8 + SQ_CORE + 8 + 20 + 8)}
ESC_INSTR = sum(ESC_ITEMS.values())
STEP_ROOT_I = ROOT_CORE + EVERY
STEP_LAZY_I = LAZY_CORE + EVERY
STEP_FLUSH_I = LAZY_CORE + EVERY + FLUSH
BLOCK_I = STEP_ROOT_I + (L_PROD - 2) * STEP_LAZY_I + STEP_FLUSH_I
assert (EVERY, FLUSH, ESC_INSTR) == (170, 12753, 24000)
assert (STEP_ROOT_I, STEP_LAZY_I, STEP_FLUSH_I, BLOCK_I) == (3244, 3627, 16380, 34132)

SAMPLE_R_MAX = 128
SAMPLE_BLOCK_INSTR = 2300
SAMPLE_BLOCK_ITEMS = {"20 rounds x 4 quarter-rounds x 12 instructions (4 add, 4 xor, 4 rotate)": 960,
           "16-word feed-forward add": 16,
           "load and store the 16-word state": 32,
           "64-bit block counter, adds/adc": 2,
           "address set-up and loop control": 10,
           "register spills, prologue and return": 250}
assert sum(SAMPLE_BLOCK_ITEMS.values()) == 1270 < SAMPLE_BLOCK_INSTR
SAMPLE_ATTEMPT_INSTR = SAMPLE_BLOCK_INSTR + 20
assert SAMPLE_ATTEMPT_INSTR == 2320
SAMPLE_MAX_INSTR = SAMPLE_R_MAX * SAMPLE_ATTEMPT_INSTR

def sample_below(rng, m, r_max=SAMPLE_R_MAX):
  if m <= 1:
    return 0
  bits = max(1, m.bit_length())
  for _ in range(r_max):
    x = rng.getrandbits(bits)
    if x < m:
      return x
  return 1

M64 = (1 << 64) - 1

def mix64(x):
  x &= M64
  x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9) & M64
  x = ((x ^ (x >> 27)) * 0x94D049BB133111EB) & M64
  return x ^ (x >> 31)

class Tower:
  """F_q2 = F_q[u]/(u^2+1), F_q12 = F_q2[w]/(w^6-xi), xi = 1+u (the blst/profile tower)."""

  def __init__(self, q):
    self.q = q
    self.count = None
    self.d = [self.f2pow((1, 1), m * (q - 1) // 6) for m in range(6)]
    z = self.f2pow((1, 1), (q * q - 1) // 6)
    assert z[1] == 0 and (z[0] * z[0] - z[0] + 1) % q == 0
    self.zeta = z[0]
    self.one = ((1, 0),) + ((0, 0),) * 5
    self._fk = {}

  def f2mul(self, a, b):
    q = self.q
    return ((a[0] * b[0] - a[1] * b[1]) % q, (a[0] * b[1] + a[1] * b[0]) % q)

  def f2pow(self, a, e):
    r = (1, 0)
    while e:
      if e & 1:
        r = self.f2mul(r, a)
      a = self.f2mul(a, a)
      e >>= 1
    return r

  def f2sqr(self, a):
    if self.count is not None:
      self.count["fp2_sqr"] += 1
    q = self.q
    return (((a[0] + a[1]) * (a[0] - a[1])) % q, (2 * a[0] * a[1]) % q)

  def f2add(self, a, b):
    if self.count is not None:
      self.count["fp2_addsub"] += 1
    q = self.q
    return ((a[0] + b[0]) % q, (a[1] + b[1]) % q)

  def f2sub(self, a, b):
    if self.count is not None:
      self.count["fp2_addsub"] += 1
    q = self.q
    return ((a[0] - b[0]) % q, (a[1] - b[1]) % q)

  def f2mulxi(self, a):
    if self.count is not None:
      self.count["fp2_addsub"] += 1
    q = self.q
    return ((a[0] - a[1]) % q, (a[0] + a[1]) % q)

  def mul(self, x, y):
    q = self.q
    acc = [[0, 0] for _ in range(11)]
    for i in range(6):
      a0, a1 = x[i]
      for j in range(6):
        b0, b1 = y[j]
        t = acc[i + j]
        t[0] += a0 * b0 - a1 * b1
        t[1] += a0 * b1 + a1 * b0
    out = []
    for m in range(6):
      c0, c1 = acc[m]
      if m + 6 <= 10:
        h0, h1 = acc[m + 6]
        c0 += h0 - h1
        c1 += h0 + h1
      out.append((c0 % q, c1 % q))
    return tuple(out)

  def pow(self, x, e):
    r = self.one
    while e:
      if e & 1:
        r = self.mul(r, x)
      x = self.mul(x, x)
      e >>= 1
    return r

  def frobk(self, y, k):
    k %= 12
    if k == 0:
      return y
    D = self._fk.get(k)
    if D is None:
      D = [self.f2pow((1, 1), m * (self.q ** k - 1) // 6) for m in range(6)]
      self._fk[k] = D
    q = self.q
    if k & 1:
      return tuple(self.f2mul((c[0], (-c[1]) % q), D[m]) for m, c in enumerate(y))
    return tuple(self.f2mul(c, D[m]) for m, c in enumerate(y))

  def keys(self, y):
    """Re(w^1 coefficient) of y^(q^k), k=0..11, from 4 multiplications in F_q."""
    q, z = self.q, self.zeta
    a = y[1][0]
    b = (y[1][0] * self.d[1][0] + y[1][1] * self.d[1][1]) % q
    za, zb = z * a % q, z * b % q
    A = [a, za, (za - a) % q, (-a) % q, (-za) % q, (a - za) % q]
    B = [b, zb, (zb - b) % q, (-b) % q, (-zb) % q, (b - zb) % q]
    out = [0] * 12
    for j in range(6):
      out[2 * j] = A[j]
      out[2 * j + 1] = B[j]
    return out

  def keys_direct(self, y):
    out = []
    for k in range(12):
      out.append(self.frobk(y, k)[1][0])
    return out

  def fp4sqr(self, a, b):
    a2, b2 = self.f2sqr(a), self.f2sqr(b)
    ab2 = self.f2sub(self.f2sub(self.f2sqr(self.f2add(a, b)), a2), b2)
    return self.f2add(a2, self.f2mulxi(b2)), ab2

  def cyc_sqr(self, y):
    """Granger-Scott squaring, valid on the cyclotomic subgroup (so on GT).
        With A=(c0,c3), B=(c1,c4), C=(c2,c5) in F_q4 = F_q2[v]/(v^2-xi), v=w^3:
          c0'=3A0-2c0, c3'=3A1+2c3, c1'=3(xi C1)+2c1, c4'=3C0-2c4, c2'=3B0-2c2, c5'=3B1+2c5
        where (A0,A1)=A^2, (B0,B1)=B^2, (C0,C1)=C^2."""
    c0, c1, c2, c3, c4, c5 = y
    A0, A1 = self.fp4sqr(c0, c3)
    B0, B1 = self.fp4sqr(c1, c4)
    C0, C1 = self.fp4sqr(c2, c5)
    xC1 = self.f2mulxi(C1)

    def minus(t, u):
      d = self.f2sub(t, u)
      return self.f2add(t, self.f2add(d, d))

    def plus(t, u):
      s = self.f2add(t, u)
      return self.f2add(t, self.f2add(s, s))

    return (minus(A0, c0), plus(xC1, c1), minus(B0, c2), plus(A1, c3), minus(C0, c4), plus(B1, c5))

  def rand_gt(self, rng, rp):
    hc = (self.q ** 12 - 1) // rp
    while True:
      x = tuple((sample_below(rng, self.q), sample_below(rng, self.q)) for _ in range(6))
      y = self.pow(x, hc)
      if y != self.one:
        return y

p = 0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab
r = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
MR384 = pow(2, 384, p)
R_FACT = {2: 32, 3: 1, 11: 1, 19: 1, 10177: 1, 125527: 1, 859267: 1, 906349: 2,
     2508409: 1, 2529403: 1, 52437899: 1, 254760293: 2}

def mult_order(a, n, fact):
  o = n
  for pf, e in fact.items():
    for _ in range(e):
      if pow(a, o // pf, n) == 1:
        o //= pf
      else:
        break
  return o

def part_A(rng, checks):
  P = Tower(p)
  lam = p % r
  prod_fact = 1
  for k, e in R_FACT.items():
    prod_fact *= k ** e
  o2 = mult_order(2, r, {k: e for k, e in R_FACT.items()})
  ol = mult_order(lam, r, R_FACT)
  t_min = o2 // math.gcd(o2, 12)
  out = {
    "r_minus_1_factorisation_verified": prod_fact == r - 1,
    "order_of_p_mod_r": ol,
    "lambda6_is_minus1": pow(lam, 6, r) == r - 1,
    "one_plus_l4_plus_l8_is_0": (1 + pow(lam, 4, r) + pow(lam, 8, r)) % r == 0,
    "log2_order_of_2_mod_r": round(math.log2(o2), 3),
    "smallest_t_with_2^t_in_<p>_log2": round(math.log2(t_min), 3),
    "zeta_primitive_6th_root": (P.zeta * P.zeta - P.zeta + 1) % p == 0,
  }
  ok = dict.fromkeys(["gt", "frob", "keys", "distinct", "canon", "inv", "sq", "trace_inv", "trace_diff"], True)
  traces = []
  counts = []
  hcof = (p ** 12 - 1) // r
  for _ in range(checks):
    x = tuple((sample_below(rng, p), sample_below(rng, p)) for _ in range(6))
    y = P.pow(x, hcof)
    ok["gt"] &= P.pow(y, r) == P.one and y != P.one
    ok["frob"] &= P.frobk(y, 1) == P.pow(y, p)
    kd, kf = P.keys_direct(y), P.keys(y)
    ok["keys"] &= kd == kf
    ok["distinct"] &= len(set(kf)) == 12
    conj = [P.frobk(y, k) for k in range(12)]
    canon = set()
    for c in conj:
      kk = P.keys(c)
      kap = min(range(12), key=lambda i: kk[i])
      canon.add(P.frobk(c, kap))
    ok["canon"] &= len(canon) == 1
    ok["inv"] &= P.mul(conj[6], y) == P.one
    ok["trace_inv"] &= len({c[0][0] for c in conj}) == 1
    traces.append(y[0][0])
    P.count = {"fp2_sqr": 0, "fp2_addsub": 0}
    s = P.cyc_sqr(y)
    counts.append(dict(P.count))
    P.count = None
    ok["sq"] &= s == P.mul(y, y)
  ok["trace_diff"] &= len(set(traces)) == len(traces)
  out.update({"random_GT_elements": checks, **{k + "_ok": v for k, v in ok.items()}})
  out["granger_scott_op_counts"] = counts[0] if counts else None
  return out

def part_S(checks=20000):
  """Exercise the bounded sampler of proof.md 0-iv: exact small-modulus reduction, R_MAX cap, residue 1."""
  class Src(object):
    def __init__(self, vals=None, reject=False):
      self.vals, self.reject, self.calls = list(vals or []), reject, 0

    def getrandbits(self, k):
      self.calls += 1
      return (1 << k) - 1 if self.reject or not self.vals else self.vals.pop(0)

  seq = Src([0, 1, 2, 3, 4, 5, 3])
  accepted = [sample_below(seq, 5) for _ in range(6)]
  exact = accepted == [0, 1, 2, 3, 4, 3] and seq.calls == 7
  ra = Src(reject=True)
  v = sample_below(ra, r)
  rng_s = random.Random(SEED ^ 0x5A5A)
  vals = [sample_below(rng_s, r) for _ in range(checks)]
  out = {"r_max_attempts": SAMPLE_R_MAX, "attempt_instr": SAMPLE_ATTEMPT_INSTR,
     "max_instr_per_sample": SAMPLE_MAX_INSTR, "max_units_per_sample": SAMPLE_MAX_INSTR / 20000,
     "acceptance_at_least": "r / 2^%d >= 1/2" % r.bit_length(), "exhaustion_rule": "return 1",
     "exhaustion_probability_log2_at_most": -SAMPLE_R_MAX, "exhaustion_allowance_log2_charged": -71.99,
     "chacha20_block_item_instr": sum(SAMPLE_BLOCK_ITEMS.values()),
     "chacha20_block_charged_instr": SAMPLE_BLOCK_INSTR,
     "exact_reduction_small_modulus": exact, "all_reject_returns": v, "all_reject_attempts": ra.calls,
     "seeded_draws_checked": checks, "seeded_all_in_range": all(0 <= x < r for x in vals)}
  out["all_checks_pass"] = bool(exact and v == 1 and ra.calls == SAMPLE_R_MAX and out["seeded_all_in_range"])
  assert out["all_checks_pass"]
  return out

def toy_setup(q):
  phi = q ** 4 - q * q + 1
  rp = phi
  d = 2
  while d * d <= rp:
    while rp % d == 0 and rp > d:
      rp //= d
    d += 1
  return rp

def classinfo(T, y):
  """kappa, canonical conjugate, minimum key. Production never builds the canonical
    element; the toy does, and breaks equal-key ties by comparing full conjugates."""
  kf = T.keys(y)
  m = min(kf)
  ks = [k for k in range(12) if kf[k] == m]
  if len(ks) == 1:
    kap = ks[0]
    can = T.frobk(y, kap)
  else:
    can, kap = min((T.frobk(y, k), k) for k in ks)
  return kap, can, m

def toy_hash(can):
  h = 0x9E3779B97F4A7C15
  for c in can:
    h = mix64(h ^ (c[0] * 0x100000001B3 + c[1]))
  return h

def heapsort_level_bound(n):
  """Worst-case sift-down levels of heapsort on n records (proof.md 0-ii): at most n for Floyd's
    construction plus sum_{k=2..n} floor(log2 k) for the n - 1 extractions, summed in power-of-two blocks."""
  tot, j = 0, 1
  while (1 << j) <= n:
    tot += j * (min(n, (1 << (j + 1)) - 1) - (1 << j) + 1)
    j += 1
  return n + tot

def heapsort(a):
  """In-place heapsort of the list a (Floyd construction, then n - 1 extractions); returns the number of
    sift-down levels executed (a level = one pass that compares a node with its children)."""
  n, levels = len(a), 0

  def sift(i, m):
    nonlocal levels
    while 2 * i + 1 < m:
      levels += 1
      c = 2 * i + 1
      if c + 1 < m and a[c + 1] > a[c]:
        c += 1
      if a[c] > a[i]:
        a[i], a[c] = a[c], a[i]
        i = c
      else:
        break
  for i in range(n // 2 - 1, -1, -1):
    sift(i, n)
  for m in range(n - 1, 0, -1):
    a[0], a[m] = a[m], a[0]
    sift(0, m)
  return levels

def match_log(log, replay, verify, rp, pair_budget):
  """The matcher of proof.md 0-ii on an append-only log of (F, alpha0, beta0) records, one per trail:
    heapsort of the whole record, one adjacent comparison per record to drop exact duplicates, then every
    pair inside each run of equal F under a hard pair budget. A pair's exponents come only from replaying
    its two trails from their start pairs to their distinguished points (`replay`); the matcher halts at
    the first pair whose candidate verifies."""
  levels = heapsort(log)
  if levels > heapsort_level_bound(len(log)):
    raise RuntimeError('heapsort exceeded its worst-case level bound')
  ded = []
  for rec in log:
    if not ded or ded[-1] != rec:
      ded.append(rec)
  pairs, replays, i = 0, [0, 0], 0
  while i < len(ded):
    j = i
    while j + 1 < len(ded) and ded[j + 1][0] == ded[i][0]:
      j += 1
    for u in range(i, j + 1):
      for v in range(u + 1, j + 1):
        pairs += 1
        if pairs > pair_budget:
          return None, levels, pairs, len(log) - len(ded), replays
        (F1, a01, b01), (F2, a02, b02) = ded[u], ded[v]
        A1, B1, n1 = replay(F1, a01, b01)
        A2, B2, n2 = replay(F2, a02, b02)
        replays[0] += 2
        replays[1] = max(replays[1], n1, n2)
        if B1 != B2:
          x = (A1 - A2) * pow(B2 - B1, -1, rp) % rp
          if verify(x):
            return x, levels, pairs, len(log) - len(ded), replays
    i = j + 1
  return None, levels, pairs, len(log) - len(ded), replays

def run_toy(T, q, rp, rng, walks, Rm, branch_bits, sq_cut, dp_bits, cap_factor=40, pair_budget=2 ** 16):
  qk = [pow(q, k, rp) for k in range(12)]
  g = T.rand_gt(rng, rp)
  xs = 1 + sample_below(rng, rp - 1)
  z0 = T.pow(g, xs)
  ah = [sample_below(rng, rp) for _ in range(Rm)]
  bh = [sample_below(rng, rp) for _ in range(Rm)]
  Mh = [T.mul(T.pow(g, ah[i]), T.pow(z0, bh[i])) for i in range(Rm)]
  Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
  At = [[ah[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
  Bt = [[bh[i] * qk[k] % rp for k in range(12)] for i in range(Rm)]
  mask = (1 << dp_bits) - 1
  window = walks * 20 << dp_bits
  gp, zp = [g], [z0]
  for _ in range(rp.bit_length()):
    gp.append(T.mul(gp[-1], gp[-1]))
    zp.append(T.mul(zp[-1], zp[-1]))

  def fpow(tab, e):
    """Fixed-base power from the squares table: the same reduced element as T.pow (only the toy's
        running time differs; part B's outputs are unchanged)."""
    out_, i = T.one, 0
    while e:
      if e & 1:
        out_ = T.mul(out_, tab[i])
      e >>= 1
      i += 1
    return out_

  def trail(a0, b0):
    """A trail's state from its start pair. (a, b) is carried for the toy's instrumentation (T_c) and
        for replay; the records written by the walk hold only (F, a0, b0)."""
    return {"y": T.mul(fpow(gp, a0), fpow(zp, b0)), "a0": a0, "b0": b0, "a": a0, "b": b0,
        "hist": [], "ck": None, "n": 0, "since": 0}

  def new_walk():
    return trail(sample_below(rng, rp), sample_below(rng, rp))

  def advance(y, a, b):
    kap, can, K = classinfo(T, y)
    u = toy_hash(can)
    if (u & ((1 << branch_bits) - 1)) < sq_cut:
      return T.cyc_sqr(y), 2 * a % rp, 2 * b % rp, kap, can, K, u
    h = (u >> branch_bits) % Rm
    idx = (12 - kap) % 12
    return T.mul(y, Tab[h][idx]), (a + At[h][idx]) % rp, (b + Bt[h][idx]) % rp, kap, can, K, u

  TOY_K, TOY_BUDGET, TOY_AMAX = 1, 2 ** 40, 2 ** 25
  TOY_L = TOY_BUDGET - TOY_AMAX
  MUL_I, SQ_I, ESC_I, REP_I, RST_I, SW_I = (STEP_FLUSH_I, STEP_FLUSH_I, ESC_INSTR, 4000,
                       int(2 ** 10.1 * 20000) + 43, 600)
  meter = [0, 0]

  def charge(x):
    if meter[0] >= TOY_L:
      return False
    assert x <= TOY_AMAX
    meter[0] += x
    meter[1] = max(meter[1], meter[0])
    if meter[0] > TOY_BUDGET:
      raise RuntimeError('toy meter exceeded its budget')
    return True

  stats = {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "dp_restarts": 0, "duplicates_removed": 0,
      "sq": 0, "ord": 0, "records": 0, "heapsort_levels": 0, "heapsort_level_bound": 0, "pairs": 0,
      "matcher_runs": 0, "replays": 0, "replay_iterations_max": 0}
  stats["meter_stop"] = 0

  def iterate(w, pay, on_visit=None, replaying=False):
    """One iteration of a walk: the same code for the walk phase and for replay. Returns 'stop' (the
        pre-debit compare refused a charge), 'dp' (the class is distinguished: the trail ends here) or
        'ok'. Replay passes a charge function that always accepts and never restarts a trail."""
    y, a, b = w["y"], w["a"], w["b"]
    kap, can, K = classinfo(T, y)
    u = toy_hash(can)
    if on_visit is not None:
      on_visit(can, a * qk[kap] % rp, b * qk[kap] % rp)
    if (u >> (branch_bits + 10)) & mask == 0:
      return "dp", u >> (branch_bits + 10 + dp_bits), kap
    if not pay(SQ_I if (u & ((1 << branch_bits) - 1)) < sq_cut else MUL_I):
      return "stop", None, None
    hist = w["hist"]
    key = (K, can)
    if len(hist) >= 2 and hist[-2][0][1] == can:
      if not replaying:
        stats["two_cycles"] += 1
      if not pay(ESC_I):
        return "stop", None, None
      s1, s2 = hist[-1], (key, y, a, b)
      _, yy, aa, bb = min(s1, s2, key=lambda s: s[0])
      w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
      w["hist"], w["ck"] = [], None
      return "esc", None, None
    if w["ck"] is not None and w["ck"][0] == can and w["ck"][1] < w["n"]:
      if not replaying:
        stats["long_cycles"] += 1
      best = (key, y, a, b)
      yy, aa, bb = y, a, b
      if not pay(ESC_I + MUL_I * (w["n"] - w["ck"][1])):
        return "stop", None, None
      for _ in range(w["n"] - w["ck"][1]):
        yy, aa, bb, kp, cn, Kn, un = advance(yy, aa, bb)
        k2, c2, K2 = classinfo(T, yy)
        if (K2, c2) < best[0]:
          best = ((K2, c2), yy, aa, bb)
      _, yy, aa, bb = best
      w["y"], w["a"], w["b"] = T.cyc_sqr(yy), 2 * aa % rp, 2 * bb % rp
      w["hist"], w["ck"] = [], None
      return "esc", None, None
    if w["n"] % 64 == 0:
      w["ck"] = (can, w["n"])
    hist.append((key, y, a, b))
    if len(hist) > 2:
      hist.pop(0)
    if not replaying:
      stats["ord"] += 1
      stats["sq"] += (u & ((1 << branch_bits) - 1)) < sq_cut
    w["y"], w["a"], w["b"], _, _, _, _ = advance(y, a, b)
    w["n"] += 1
    w["since"] += 1
    return "ok", None, None

  def replay(F, a0, b0):
    """Re-run the trail that starts at (a0, b0) to its distinguished point, tracking exponents, and
        return the canonical exponent pair there and the number of iterations replayed."""
    w = trail(a0, b0)
    its = 0
    while True:
      st, Fd, kap = iterate(w, lambda x: True, replaying=True)
      its += 1
      if st == "dp":
        if Fd != F:
          raise RuntimeError('replayed trail ends at a different distinguished class')
        return w["a"] * qk[kap] % rp, w["b"] * qk[kap] % rp, its
      if its > 2 * ((cap_factor << dp_bits) + 1) + 3:
        raise RuntimeError('replay exceeded the trail bound')

  W = [new_walk() for _ in range(walks)]
  visited, log = {}, []
  steps, Tc, found, check_at = 0, None, None, None
  state = {"Tc": None}

  def on_visit(can, ca, cb):
    if state["Tc"] is None:
      old = visited.get(can)
      if old is None:
        visited[can] = (ca, cb)
      elif old != (ca, cb):
        state["Tc"] = steps
        stats["ord_at_Tc"], stats["sq_at_Tc"] = stats["ord"], stats["sq"]
        stats["esc2_at_Tc"], stats["escL_at_Tc"] = stats["two_cycles"], stats["long_cycles"]

  while found is None:
    for w in W:
      if not charge(SW_I):
        stats["meter_stop"] = 1
        break
      st, Fd, _ = iterate(w, charge, on_visit)
      if Tc is None and state["Tc"] is not None:
        Tc = state["Tc"]
        check_at = Tc + window
      steps += 1
      if st == "stop":
        stats["meter_stop"] = 1
        break
      if st == "dp":
        if not charge(REP_I):
          stats["meter_stop"] = 1
          break
        log.append((Fd, w["a0"], w["b0"]))
        if not charge(RST_I):
          stats["meter_stop"] = 1
          break
        stats["dp_restarts"] += 1
        w.clear()
        w.update(new_walk())
        continue
      if w["since"] > cap_factor << dp_bits:
        if not charge(RST_I):
          stats["meter_stop"] = 1
          break
        stats["restarts"] += 1
        w.clear()
        w.update(new_walk())
    if stats["meter_stop"]:
      break
    if check_at is not None and steps >= check_at:
      stats["matcher_runs"] += 1
      snap = list(log)
      x, lv, pr, dup, rpl = match_log(snap, replay, lambda v: T.pow(g, v) == z0, rp, pair_budget)
      stats["records"], stats["heapsort_levels"] = len(snap), lv
      stats["heapsort_level_bound"] = heapsort_level_bound(len(snap))
      stats["pairs"], stats["duplicates_removed"] = pr, dup
      stats["replays"] += rpl[0]
      stats["replay_iterations_max"] = max(stats["replay_iterations_max"], rpl[1])
      if pr > pair_budget:
        raise RuntimeError('toy pair budget exceeded')
      if x is not None:
        if x != xs:
          raise RuntimeError('recovered logarithm differs from the planted one')
        found = x
      else:
        check_at = steps + window
  stats["meter_max"], stats["meter_budget"] = meter[1], TOY_BUDGET
  return Tc, steps, stats

# Toy success metric in the ledger's units: walk steps to the first useful collision (wasted 2-cycle steps included),
# each 2-cycle escape at 24,000 instructions in steps of c, each checkpoint escape at its worst case (64 steps at
# 2.8794 c, and the escape), against S (1 + eps), eps the production eps rounded up, 0.02379221. DP and fresh-start
# iterations are not counted (production has one per 2^86 steps).
TOY_ESC_STEPS = 24000 / 20000 / (34132 / 120000)
TOY_ESCL_STEPS = 64 * 16380 / 34132 * 6 + TOY_ESC_STEPS
TOY_EPS = 0.02379221
TOY_PM = 0.4005

def part_B(rng, q, trials, Rm, branch_bits, sq_cut, dp_bits, walks):
  T = Tower(q)
  rp = toy_setup(q)
  if not (pow(q, 12, rp) == 1 and all(pow(q, k, rp) != 1 for k in range(1, 12))):
    raise RuntimeError('toy field does not have embedding degree 12')
  N = (rp - 1) // 12
  s = sq_cut / (1 << branch_bits)
  sumq2 = s * s + (1 - s) ** 2 / Rm
  Np = N / (1 - sumq2)
  Tcs, tot, ev, dev = [], [], {"two_cycles": 0, "long_cycles": 0, "restarts": 0, "dp_restarts": 0,
                "duplicates_removed": 0, "records": 0, "heapsort_levels": 0,
                "heapsort_level_bound": 0, "matcher_runs": 0, "replays": 0}, []
  max_pairs, max_level_ratio, max_replay = 0, 0.0, 0
  meter_stops, meter_ratio = 0, 0.0
  pool_ord, pool_sq = 0, 0
  Teq, pre = [], [0, 0, 0]
  while len(Tcs) < trials:
    Tc, steps, st = run_toy(T, q, rp, rng, walks, Rm, branch_bits, sq_cut, dp_bits)
    meter_ratio = max(meter_ratio, st["meter_max"] / st["meter_budget"])
    if st["meter_stop"]:
      meter_stops += 1
      continue
    Tcs.append(Tc)
    Teq.append(st["ord_at_Tc"] + TOY_ESC_STEPS * st["esc2_at_Tc"] + TOY_ESCL_STEPS * st["escL_at_Tc"])
    pre[0] += st["ord_at_Tc"]
    pre[1] += st["esc2_at_Tc"]
    pre[2] += st["escL_at_Tc"]
    tot.append(steps)
    for k in ev:
      ev[k] += st[k]
    max_pairs = max(max_pairs, st["pairs"])
    max_replay = max(max_replay, st["replay_iterations_max"])
    max_level_ratio = max(max_level_ratio, st["heapsort_levels"] / max(1, st["heapsort_level_bound"]))
    pool_ord += st["ord"]
  n = len(Tcs)
  sN = math.sqrt(N)
  pred = math.sqrt(math.pi * Np / 2)
  s40 = math.sqrt(2 * Np * math.log(1 / 0.60))
  s405 = math.sqrt(2 * Np * math.log(1 / 0.595))
  s4005 = math.sqrt(2 * Np * math.log(1 / 0.5995))
  sPm = math.sqrt(2 * Np * math.log(1 / (1 - TOY_PM)))
  s42 = math.sqrt(2 * Np * math.log(1 / 0.58))
  s43 = math.sqrt(2 * Np * math.log(1 / 0.57))
  allsteps = max(1, sum(tot))
  return {"q": q, "r_prime": rp, "classes_N": N, "R": Rm, "branch_bits": branch_bits,
      "sq_cut": sq_cut, "squaring_prob": s,
      "sum_q2": round(sumq2, 6), "effective_N_over_N": round(1 / (1 - sumq2), 4),
      "toy_model_point": TOY_PM, "toy_eps": TOY_EPS,
      "trials_completed": n, "logs_recovered_and_verified": n,
      "mean_Tc_over_sqrtN": round(sum(Tcs) / n / sN, 4) if n else None,
      "model_mean_over_sqrtN": round(pred / sN, 4),
      "fraction_Tc_le_S_model_0.40": round(sum(t <= s40 for t in Tcs) / max(1, n), 3),
      "fraction_Tc_le_S_model_0.405": round(sum(t <= s405 for t in Tcs) / max(1, n), 3),
      "fraction_Tc_le_S_model_0.4005": round(sum(t <= s4005 for t in Tcs) / max(1, n), 3),
      "fraction_charged_steps_le_S_model_Pm_times_1_plus_eps": round(sum(t <= sPm * (1 + TOY_EPS) for t in Teq) / max(1, n), 4),
      "successes_charged_steps_le_S_model_Pm_times_1_plus_eps": sum(t <= sPm * (1 + TOY_EPS) for t in Teq),
      "fraction_charged_steps_le_S_model_0.4005_times_1_plus_eps": round(sum(t <= s4005 * (1 + TOY_EPS) for t in Teq) / max(1, n), 4),
      "mean_charged_steps_at_Tc_over_sqrtN": round(sum(Teq) / max(1, n) / sN, 4),
      "steps_before_Tc": pre[0], "two_cycle_escapes_before_Tc": pre[1], "checkpoint_escapes_before_Tc": pre[2],
      "two_cycles_per_step_before_Tc": round(pre[1] / max(1, pre[0]), 6),
      "fraction_Tc_le_S_model_0.42": round(sum(t <= s42 for t in Tcs) / max(1, n), 3),
      "fraction_Tc_le_S_model_0.43": round(sum(t <= s43 for t in Tcs) / max(1, n), 3),
      "successes_within_S_model_0.43": sum(t <= s43 for t in Tcs),
      "binomial_standard_error_of_fraction_0.43": round(math.sqrt((sum(t <= s43 for t in Tcs) / max(1, n)) * (1 - sum(t <= s43 for t in Tcs) / max(1, n)) / max(1, n)), 4),
      "S_model_0.40_over_sqrtN": round(s40 / sN, 4), "S_model_0.405_over_sqrtN": round(s405 / sN, 4), "S_model_0.4005_over_sqrtN": round(s4005 / sN, 4), "S_model_0.42_over_sqrtN": round(s42 / sN, 4),
      "S_model_0.43_over_sqrtN": round(s43 / sN, 4),
      "mean_total_steps_incl_DP_delay_over_sqrtN": round(sum(tot) / max(1, n) / sN, 4),
      "branch_rule_steps_pooled_over_trials": pool_ord,
      "squaring_fraction_pooled_over_trials": round(pool_sq / max(1, pool_ord), 4),
      "events": ev, "two_cycles_per_step": ev["two_cycles"] / allsteps,
      "matcher": "append-only log of (F, alpha0, beta0) trail records, in-place heapsort of the whole record, "
           "adjacent dedup, pairs inside equal-F runs under a hard pair budget of 2^16, each pair's "
           "exponents recovered only by replaying its two trails, halt at the first verified pair",
      "max_replay_iterations_in_a_trial": max_replay,
      "replay_iteration_bound": 2 * ((40 << dp_bits) + 1) + 3,
      "max_pairs_charged_in_a_trial": max_pairs,
      "schedule": "one processor, walks round-robin, chunk of 1 step iteration per turn, charged switch, one pre-debit meter; "
            "every trail ends at its first distinguished point and the walk restarts from a fresh start",
      "meter_budget_instr_log2": 40, "max_meter_over_budget": round(meter_ratio, 8), "meter_stops": meter_stops,
      "max_heapsort_levels_over_worst_case_bound": round(max_level_ratio, 4),
      "Tc_over_sqrtN": [round(t / sN, 4) for t in Tcs] if cfg.get("dump") else None}

def part_E(rng, q, m, Rm, branch_bits, sq_cut):
  """Sample m distinct classes, apply one step F to each, and count image coincidences.
    H1 predicts N * P[F(c1) = F(c2)] = 1 - s^2 - (1-s)^2/(12R) for distinct classes c1, c2."""
  T = Tower(q)
  rp = toy_setup(q)
  N = (rp - 1) // 12
  g = T.rand_gt(rng, rp)
  Mh = [T.pow(g, 1 + sample_below(rng, rp - 1)) for _ in range(Rm)]
  Tab = [[T.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
  gp = [g]
  for _ in range(rp.bit_length()):
    gp.append(T.mul(gp[-1], gp[-1]))
  classes = {}
  while len(classes) < m:
    e = 1 + sample_below(rng, rp - 1)
    y, j = T.one, 0
    while e:
      if e & 1:
        y = T.mul(y, gp[j])
      e >>= 1
      j += 1
    kap, can, _ = classinfo(T, y)
    if can not in classes:
      classes[can] = (y, kap)
  images = {}
  for can, (y, kap) in classes.items():
    u = toy_hash(can)
    if (u & ((1 << branch_bits) - 1)) < sq_cut:
      y2 = T.cyc_sqr(y)
    else:
      y2 = T.mul(y, Tab[(u >> branch_bits) % Rm][(12 - kap) % 12])
    _, can2, _ = classinfo(T, y2)
    images[can2] = images.get(can2, 0) + 1
  pairs = sum(n * (n - 1) // 2 for n in images.values())
  tot = m * (m - 1) // 2
  s = sq_cut / (1 << branch_bits)
  obs = pairs / tot * N
  return {"q": q, "classes_N": N, "sampled_classes": m, "R": Rm, "branch_bits": branch_bits,
      "sq_cut": sq_cut, "squaring_prob": s,
      "image_coincidence_pairs": pairs,
      "N_times_pair_collision_rate_observed": round(obs, 4),
      "relative_standard_error": round(1 / math.sqrt(max(1, pairs)), 4),
      "predicted_exact_1_minus_s2_minus_(1-s)^2/(12R)": round(1 - s * s - (1 - s) ** 2 / (12 * Rm), 4),
      "predicted_conservative_1_minus_s2_minus_(1-s)^2/R": round(1 - s * s - (1 - s) ** 2 / Rm, 4),
      "random_mapping": 1.0}

# Shipped kernels in the emulator (part K; kernels/harness.py) and the production lazy walk driven by them (part C)
sys.path.insert(0, locate("kernels"))
harness = importlib.import_module("harness")
assert (harness.p, harness.MR384) == (p, MR384)
Kernels = harness.Kernels

def stored_keys(P, y):
  return [v * MR384 % p for v in P.keys(y)]

def ref_key(P, y):
  ks = stored_keys(P, y)
  K = min(ks)
  return K, ks.index(K)

def rand_gt(P, rng):
  hcof = (p ** 12 - 1) // r
  return P.pow(tuple((sample_below(rng, p), sample_below(rng, p)) for _ in range(6)), hcof)

def part_K(rng, KN, P, gens, cases):
  """Each shipped kernel on real GT values against the big-integer reference: mul.s (x * m), sqr.s (x^2),
    root.s and lazy.s (the stored key minimum K, its index kappa and the stored a, b of y' = x or x * Pi);
    one lazy case is also checked against keys from explicit Frobenius maps (the definition)."""
  fails = {"mul": 0, "node_build": 0, "sqr": 0, "root": 0, "lazy": 0, "per_base_data": 0, "lazy_direct_frobenius": 0}
  x, m = gens[0], gens[1]
  for i in range(cases):
    x, m = P.mul(x, gens[2]), P.mul(m, gens[3])
    fails["mul"] += KN.mul(x, KN.mul_entry(m)) != P.mul(x, m)
    fails["node_build"] += KN.mul(m, KN.mul_entry(x)) != P.mul(x, m)
    fails["sqr"] += KN.sqr(x) != P.cyc_sqr(x)
    K, kap, a, b, ps = KN.key(x)
    fails["root"] += (K, kap) != ref_key(P, x) or (a, b) != (x[1][0] * MR384 % p, x[1][1] * MR384 % p)
    fails["per_base_data"] += ps != KN.base_data(x)
    y = P.mul(x, m)
    K, kap, a, b, _ = KN.key(x, KN.gamma(m))
    fails["lazy"] += (K, kap) != ref_key(P, y) or (a, b) != (y[1][0] * MR384 % p, y[1][1] * MR384 % p)
    if i == 0:
      direct = [v * MR384 % p for v in P.keys_direct(y)]
      fails["lazy_direct_frobenius"] += K != min(direct) or direct[kap] != K
  executed = {k: sorted(v) for k, v in KN.executed.items()}
  want = {"mul": MUL_CORE, "sqr": SQ_CORE, "lazy": LAZY_CORE, "root": ROOT_CORE}
  ok = (all(v == 0 for v in fails.values()) and KN.static == want
     and all(executed[k] == [want[k]] for k in want))
  assert ok, (fails, KN.static, executed)
  return {"cases": cases, "static_counts_emitted": KN.static, "executed_counts_every_call": executed,
      "mismatches": fails, "all_checks_pass": ok}

DP_SHIFT, F_SHIFT, F_BITS = 64, 64 + THETA_BITS, 128
assert F_SHIFT == 150 and F_SHIFT + F_BITS == 278

def class_bits(K, Rm, dp_bits):
  """(h, distinguished, F) from the stored key minimum K: h = floor(Rm K_0 / 2^64) from bits 0..63, the DP
    predicate on bits 64..64+dp_bits-1 (production: 86 bits, 64..149), F = bits 150..277."""
  return ((K & M64) * Rm) >> 64, (K >> DP_SHIFT) & ((1 << dp_bits) - 1) == 0, (K >> F_SHIFT) & ((1 << F_BITS) - 1)

def part_C(rng, KN, P, g, steps, dp_bits, Rm=R_PROD, depth=L_PROD):
  """The production lazy walk in the real GT: R = 32 twisted multipliers, trie depth L = 6, the 2-cycle escape,
    trails ending at distinguished classes (test density 2^-dp_bits on the production bit field, bits 64.. of K),
    a record (F, alpha0, beta0) per trail and a fresh start. Every key comes from the emulated root.s / lazy.s,
    every flush and escape from the emulated mul.s / sqr.s; an eager walk (explicit element, keys from F2,
    exponents carried) runs beside it and every step is compared. Afterwards every recorded trail is replayed
    eagerly from its start pair alone; the replay must reach a distinguished class with the recorded F, and the
    canonical pair p^kappa (alpha, beta) must give the canonical element g^alpha^ z^beta^."""
  lam = [pow(p, k, r) for k in range(12)]
  x0 = 1 + sample_below(rng, (1 << 40) - 1)
  z0 = P.pow(g, x0)
  ah = [1 + sample_below(rng, (1 << 40) - 1) for _ in range(Rm)]
  bh = [1 + sample_below(rng, (1 << 40) - 1) for _ in range(Rm)]
  Mh = [P.mul(P.pow(g, ah[i]), P.pow(z0, bh[i])) for i in range(Rm)]
  Tab = [[P.frobk(Mh[i], k) for k in range(12)] for i in range(Rm)]
  trie = {(): P.one}

  def Pi(word):
    if word not in trie:
      h, t = word[-1]
      trie[word] = P.mul(Pi(word[:-1]), Tab[h][t])
    return trie[word]

  def start():
    a0, b0 = 1 + sample_below(rng, (1 << 40) - 1), 1 + sample_below(rng, (1 << 40) - 1)
    return a0, b0, P.mul(P.pow(g, a0), P.pow(z0, b0))

  def eager_trail(a0, b0, y, limit):
    """The walk rule with explicit elements and exponents (what a replay executes)."""
    al, be, hist, n = a0, b0, [], 0
    while n < limit:
      K, kap = ref_key(P, y)
      h, dp, F = class_bits(K, Rm, dp_bits)
      if dp:
        return F, kap, y, al, be, n
      if len(hist) >= 2 and hist[-2][0] == K:
        _, yy, aa, bb = min(hist[-1], (K, y, al, be), key=lambda s: s[0])
        y, al, be, hist = P.cyc_sqr(yy), 2 * aa % r, 2 * bb % r, []
        n += 1
        continue
      hist = (hist + [(K, y, al, be)])[-2:]
      t = (12 - kap) % 12
      y = P.mul(y, Tab[h][t])
      al, be = (al + ah[h] * lam[t]) % r, (be + bh[h] * lam[t]) % r
      n += 1
    return None

  a0, b0, ye = start()
  al, be = a0, b0
  base, word, hist = ye, (), []
  log = []
  cnt = {"root_steps": 0, "lazy_steps": 0, "flushes": 0, "two_cycle_escapes": 0, "escape_flushes": 0,
     "reports": 0, "key_mismatches": 0, "element_mismatches": 0, "kernel_vs_python_mismatches": 0}
  hcount = [0] * Rm
  for _ in range(steps):
    if P.mul(base, Pi(word)) != ye:
      cnt["element_mismatches"] += 1
    if word:
      K, kap, _, _, _ = KN.key(base, KN.gamma(Pi(word)))
      cnt["lazy_steps"] += 1
    else:
      K, kap, _, _, _ = KN.key(base)
      cnt["root_steps"] += 1
    cnt["key_mismatches"] += (K, kap) != ref_key(P, ye)
    h, dp, F = class_bits(K, Rm, dp_bits)
    if dp:
      log.append((F, a0, b0))
      cnt["reports"] += 1
      a0, b0, ye = start()
      al, be = a0, b0
      base, word, hist = ye, (), []
      continue
    if len(hist) >= 2 and hist[-2][0] == K:
      cnt["two_cycle_escapes"] += 1
      Ks, bs, ws, yes, als, bes = min(hist[-1], (K, base, word, ye, al, be), key=lambda s: s[0])
      if ws:
        expl = KN.mul(bs, KN.mul_entry(Pi(ws)))
        cnt["escape_flushes"] += 1
      else:
        expl = bs
      cnt["kernel_vs_python_mismatches"] += expl != yes
      base, word, hist = KN.sqr(expl), (), []
      ye, al, be = P.cyc_sqr(yes), 2 * als % r, 2 * bes % r
      cnt["kernel_vs_python_mismatches"] += base != ye
      continue
    hist = (hist + [(K, base, word, ye, al, be)])[-2:]
    t = (12 - kap) % 12
    hcount[h] += 1
    ye = P.mul(ye, Tab[h][t])
    al, be = (al + ah[h] * lam[t]) % r, (be + bh[h] * lam[t]) % r
    word = word + ((h, t),)
    if len(word) == depth:
      base = KN.mul(base, KN.mul_entry(Pi(word)))
      word = ()
      cnt["flushes"] += 1
      cnt["kernel_vs_python_mismatches"] += base != ye
  exp_ok = ye == P.mul(P.pow(g, al), P.pow(z0, be))
  replay_ok, replay_steps = 0, 0
  for F, a0_, b0_ in log:
    res = eager_trail(a0_, b0_, P.mul(P.pow(g, a0_), P.pow(z0, b0_)), steps + 1)
    if res is None:
      continue
    F2, kap, y, al_, be_, n = res
    canon = P.frobk(y, kap)
    ahat, bhat = al_ * lam[kap] % r, be_ * lam[kap] % r
    replay_ok += F2 == F and canon == P.mul(P.pow(g, ahat), P.pow(z0, bhat))
    replay_steps += n
  ok = (exp_ok and cnt["key_mismatches"] == cnt["element_mismatches"] == cnt["kernel_vs_python_mismatches"] == 0
     and replay_ok == len(log))
  assert ok, (cnt, replay_ok, len(log))
  n = sum(hcount)
  chi = sum((c - n / Rm) ** 2 / (n / Rm) for c in hcount)
  return {"steps": steps, "R": Rm, "trie_depth_L": depth, "test_dp_bits": dp_bits, "trie_nodes_built": len(trie) - 1,
      **cnt, "records": len(log), "records_replayed_to_the_same_F_with_canonical_pair_verified": replay_ok,
      "replayed_steps": replay_steps, "exponent_invariant_of_the_last_trail": exp_ok,
      "multiplier_index_chi2_%ddof" % (Rm - 1): round(chi, 2), "all_checks_pass": ok,
      "note": "keys by emulated root.s/lazy.s, flushes and escapes by emulated mul.s/sqr.s, compared step by step "
          "with the eager walk; the production DP field read at a test density; exponents carried only for "
          "this check and recovered for the records only by replay"}

def part_G(rng, P, g, n):
  """Statistics of the bits of K that the rule reads, on n production classes (the elements g^(c + i d), i < n,
    for random c, d: distinct classes, one multiplication each): the multiplier index h (32 bins), the low 8 bits
    of the DP field (bits 64..71, 256 bins), the low 8 bits of F (bits 150..157, 256 bins), bit balance over all 278
    bits read, and K < p/2 (keys come in negation pairs)."""
  y = P.pow(g, 1 + sample_below(rng, r - 1))
  d = P.pow(g, 1 + sample_below(rng, r - 1))
  hb, db, fb = [0] * R_PROD, [0] * 256, [0] * 256
  ones = [0] * (F_SHIFT + F_BITS)
  kmax = 0
  for _ in range(n):
    K, _ = ref_key(P, y)
    kmax = max(kmax, K)
    h, _, F = class_bits(K, R_PROD, THETA_BITS)
    hb[h] += 1
    db[(K >> DP_SHIFT) & 255] += 1
    fb[F & 255] += 1
    for j in range(F_SHIFT + F_BITS):
      ones[j] += (K >> j) & 1
    y = P.mul(y, d)

  def chi(b):
    e = n / len(b)
    return round(sum((c - e) ** 2 / e for c in b), 2)
  zmax = max(abs(o - n / 2) / math.sqrt(n / 4) for o in ones)
  return {"classes": n, "chi2_h_31dof": chi(hb), "chi2_dp_field_low8_255dof": chi(db), "chi2_F_low8_255dof": chi(fb),
      "max_abs_z_bit_balance_over_278_bits": round(zmax, 2),
      "all_K_below_half_p": kmax < p // 2 + 1}

# Ledger (part D)
CLAIM_TXT = "123.896"
P_MODEL_TXT = "0.4005"

def ledger(cyc, red):
  """Exact recomputation of the work ledger of proof.md section 5 (instruction counts / C, C = 20,000)."""
  from decimal import Decimal, getcontext, localcontext
  getcontext().prec = 60
  ln2 = Decimal(2).ln()
  lg = lambda v: float(Decimal(v).ln() / ln2)
  D = lambda fr: Decimal(fr.numerator) / Decimal(fr.denominator)
  C = C_UNIT
  Rt, Lt = R_PROD, L_PROD
  # step prices (instructions) and the per-step average over a complete block
  c = Fraction(BLOCK_I, Lt * C)
  cmax = Fraction(STEP_FLUSH_I, C)
  assert c == Fraction(34132, 120000)
  assert all(Fraction(k * STEP_LAZY_I, 1) <= k * Lt * c * C for k in range(1, Lt))
  assert STEP_ROOT_I < STEP_LAZY_I < c * C < STEP_FLUSH_I
  # partial blocks: k < L steps cost at most STEP_ROOT_I + (k-1) STEP_LAZY_I < k c C
  assert all(STEP_ROOT_I + (k - 1) * STEP_LAZY_I < k * c * C for k in range(1, Lt))
  # collision law
  N = Fraction(r - 1, 12)
  Pm = Decimal(P_MODEL_TXT)
  sumq2 = Fraction(1, Rt)
  Np = N / (1 - sumq2)
  Npd = D(Np)
  S = (2 * Npd * (1 / (1 - Pm)).ln()).sqrt()
  # fruitless closed walks (part F): eps in step equivalents, 1.25 x the H2 model
  esc_steps = Fraction(ESC_INSTR, C) / c
  P2 = Fraction(1, 12 * Rt)
  eps2 = P2 * (2 + esc_steps)
  eps_long = sum(Fraction(pl).limit_denominator(10 ** 30) * (64 + l + l * cmax / c + esc_steps)
         for l, pl in zip(range(3, 65), cyc["P_l_by_length_3_to_64"]))
  EPS_MARGIN = Fraction(5, 4)
  eps = D(EPS_MARGIN * (eps2 + eps_long))
  assert cyc["P_2_equals_1_over_12R"]
  LONG_RATE = Decimal(2) ** -121
  assert cyc["rate_useless_walks_longer_than_64_log2"] < -122.6
  # trails, restarts, detection window
  WALKS, RESTART_CAP = 2 ** 20, 2 ** 20
  STUCK = STUCK_M * Decimal(2) ** THETA_BITS
  assert STUCK == 0x7000000 * Decimal(2) ** 64
  E = (WALKS + RESTART_CAP) * STUCK
  # Iteration maximum: a block of k <= 6 iterations with the escape, or report and fresh start, ending it is charged
  # at least k c C, so n <= B_walk / (c C) + 6 M (each walk's last block), B_walk < W_cap C.
  for k in range(1, Lt + 1):
    blk = STEP_ROOT_I + (k - 1) * STEP_LAZY_I + (FLUSH if k == Lt else 0)
    assert blk + ESC_INSTR >= k * c * C and blk + 4000 + int(2 ** 10) * C >= k * c * C
    assert k < Lt or blk == Lt * c * C
  N_MAX = Decimal(2) ** Decimal(CLAIM_TXT) / D(c) + 6 * 2 ** 20
  EXP_STUCK = N_MAX * Decimal(2) ** (-THETA_BITS) * (-Decimal(STUCK_M)).exp()
  EXP_LONG = N_MAX * LONG_RATE
  EXP_RESTARTS = EXP_STUCK + EXP_LONG
  RESTART_FRAC = Decimal(RESTART_CAP) * STUCK / N_MAX
  # caps: each count's expectation under H2's per-step rates, P(cap reached) <= E / cap (Markov)
  REPORT_CAP, SLOT_BYTES, PAIR_BUDGET = 2 ** 55, 80, 2 ** 4
  MU_REPORTS = N_MAX * Decimal(2) ** (-THETA_BITS)
  # pairs charged before the halting one: false-key (2^-F_BITS each) or beta^1 = beta^2 (1/r < 2^-254 each)
  assert r > 2 ** 254
  EXP_PAIRS = Decimal(REPORT_CAP) * (REPORT_CAP - 1) / 2 * (Decimal(2) ** -F_BITS + Decimal(2) ** -254)
  CAPS = {"report cap 2^55 reached": (MU_REPORTS, REPORT_CAP), "restart cap 2^20 reached": (EXP_RESTARTS, RESTART_CAP),
      "pair budget 2^4 exceeded": (EXP_PAIRS, PAIR_BUDGET)}
  markov_ratio_report, markov_ratio_restart, markov_ratio_pairs = (cap / e for e, cap in CAPS.values())
  REPORT_ITEMS = {"F from K": 5, "start pair": 4, "pre-debit compare and add": 8, "report counter": 5,
          "slot address": 2, "80-byte record": 12, "stuck counter reset": 4, "page zero-fill": 13,
          "margin": 3947}
  report_instr = sum(REPORT_ITEMS.values())
  assert report_instr == 4000 and Fraction(2 * 131072 + 65536, 2 ** 21 // 80) <= Fraction("12.5002") < 13
  # sampler
  SAMPLE_TABLE_DRAWS = 2 * Rt
  SAMPLES_MAX = Decimal(SAMPLE_TABLE_DRAWS + 2 * (WALKS + REPORT_CAP + RESTART_CAP))
  assert lg(SAMPLES_MAX) < Decimal("56.0001") and SAMPLES_MAX * SAMPLE_R_MAX < Decimal(2) ** 64
  # setup: g, z check, table, trie (each node 1.6 units), zero-fill of the trie pages, walk starts
  NODES_INT = sum((12 * Rt) ** d for d in range(1, Lt))
  NODES_LEAF = (12 * Rt) ** Lt
  INT_BYTES, LEAF_BYTES = 20 * 48 + 54 * 48, 54 * 48
  TRIE_BYTES = NODES_INT * INT_BYTES + NODES_LEAF * LEAF_BYTES
  node_items = {"argument set-up: addresses of T[h][t]'s 12 raw coordinates and of the parent's entry, margin": 72,
         "mul.s: Pi_child = T[h][t] * Pi_parent (table entry's raw coordinates against the parent's 54 forms)": MUL_CORE + 8,
         "the child's 54 mul.s table-side forms reduced mod p: 270 modular additions at 60": 270 * 60,
         "node entry gamma: four xi products (8 modular add/sub at 60)": 8 * 60,
         "its six pre-sums mod p and Ea, Ed: root.s's per-base sequence on gamma (1,753 executed), two Montgomery "
         "reductions with q-table subtraction, store and constant reloads (2 x 188), the offset p (2 x 12)": 1753 + 376 + 24,
         "store the 74 vectors (444 words, 222 stp) and address/loop control": 222 + 50}
  NODE_INSTR = sum(node_items.values())
  NODE_UNITS = Fraction(8, 5)
  assert NODE_INSTR == 31920 and NODE_INSTR <= NODE_UNITS * C
  ZERO_FILL_PER_BYTE = Fraction(2 * 131072 + 65536, 2 ** 21)
  table_units = (Decimal(Rt * (1021 + 11)) + Decimal(SAMPLE_TABLE_DRAWS * SAMPLE_MAX_INSTR) / C
         + Decimal(2 * 12 * Rt * 800 + 12 * Rt * 270 * 60) / C)
  trie_units = D(NODE_UNITS) * (NODES_INT + NODES_LEAF) + D(ZERO_FILL_PER_BYTE * TRIE_BYTES) / C
  ONE_START = Decimal(2) ** Decimal("10.1")
  START_MAX = Decimal(1021) + Decimal(2 * SAMPLE_MAX_INSTR) / C + Decimal(43) / C
  assert START_MAX < ONE_START
  SETUP_ITEMS = (Decimal(2) ** Decimal("14.1") + 510 + table_units + trie_units + WALKS * ONE_START
         + Decimal(32 + 2 * 1940) / C)
  assert red["pairing_units_counted_g"] <= 2 ** 14 and red["pairing_units_counted_verification"] <= 2 ** 14
  Q_OUT = red["g2_output_units_hard_max_255_bit_scalar"]
  assert Q_OUT <= 12809
  # The trie has 12R = 384 children per level: 2^52.21 units, the least hundredth above SETUP_ITEMS.
  SETUP_UNITS = Decimal(2) ** Decimal("52.21")
  assert SETUP_ITEMS <= SETUP_UNITS
  # walk-phase non-step work
  METER_UPDATE = 40
  restart_count_units = Decimal(RESTART_CAP) * 3 / C
  dp_units = Decimal(REPORT_CAP) * report_instr / C
  start_units = Decimal(WALKS) * METER_UPDATE / C
  restart_units = Decimal(RESTART_CAP) * METER_UPDATE / C
  restart_start_units = Decimal(RESTART_CAP) * ONE_START
  dp_restart_units = Decimal(REPORT_CAP) * (ONE_START + Decimal(METER_UPDATE) / C)
  WALK_NONSTEP = dp_units + start_units + restart_units + restart_start_units + restart_count_units + dp_restart_units
  CHUNK_K = 2 ** 63
  SWITCH_INSTR = 600
  SWITCHES = N_MAX / CHUNK_K + WALKS + RESTART_CAP + 1
  SWITCH_UNITS = SWITCHES * SWITCH_INSTR / C
  # stopping allowance: heapsort, dedup, pairs (two replays each), output, verification
  LEVEL_INSTR, OUTER_INSTR = 132, 40
  LEVELS = heapsort_level_bound(REPORT_CAP)
  assert LEVELS == 54 * 2 ** 55 + 57
  SORT_INSTR = LEVEL_INSTR * LEVELS + OUTER_INSTR * (REPORT_CAP // 2 + REPORT_CAP - 1)
  DEDUP_INSTR = 100 * REPORT_CAP
  A_MAX = 2 ** 25
  ESC_LONG_INSTR = ESC_INSTR + 64 * STEP_FLUSH_I + 65 * 64
  assert ESC_LONG_INSTR == 1076480 <= A_MAX
  assert START_MAX <= Decimal("1050.7") and ONE_START - START_MAX > 46
  START_I = ONE_START * C + 43
  assert all(a < A_MAX for a in (STEP_FLUSH_I, ESC_INSTR, ESC_LONG_INSTR, report_instr, SWITCH_INSTR, START_I,
                                 STEP_FLUSH_I + report_instr + ESC_LONG_INSTR, report_instr + START_I))
  REPLAY_ITERS = 2 * STUCK + 3
  REPLAY_ITER_INSTR = 2 * A_MAX
  REPLAY_UNITS = REPLAY_ITERS * REPLAY_ITER_INSTR / C
  PAIR_UNITS = 2 * REPLAY_UNITS + 2 * Decimal(4000) / C + 600
  OUTPUT_UNITS = 12809 + 2 ** 14
  R_STOP_INSTR = SORT_INSTR + DEDUP_INSTR + 100 + int(((PAIR_BUDGET * PAIR_UNITS + OUTPUT_UNITS) * C)
                             .to_integral_value(rounding="ROUND_CEILING"))
  R_STOP = Decimal(R_STOP_INSTR) / C
  NONWALK = WALK_NONSTEP + SETUP_UNITS + R_STOP
  cd = D(c)
  W = cd * S * (1 + eps) + cd * E + NONWALK + SWITCH_UNITS
  claim_thousandths = int((W.ln() / ln2 * 1000).to_integral_value(rounding="ROUND_CEILING"))
  CLAIM = Decimal(claim_thousandths) / 1000
  assert CLAIM == Decimal(CLAIM_TXT), (CLAIM, lg(W))
  W_up = W.quantize(Decimal(1).scaleb(int(W.adjusted()) - 40), rounding="ROUND_CEILING")
  fnum, fden = Fraction(W_up).numerator, Fraction(W_up).denominator
  assert fnum ** 1000 <= fden ** 1000 * 2 ** claim_thousandths, "ledger exceeds the claim"
  assert fnum ** 1000 > fden ** 1000 * 2 ** (claim_thousandths - 1), "claim is not the least thousandth"
  with localcontext() as ctx:
    ctx.prec = 90
    WCAP_INSTR = int((Decimal(2) ** CLAIM * C).to_integral_value(rounding="ROUND_FLOOR"))
  assert WCAP_INSTR ** 1000 <= 2 ** claim_thousandths * C ** 1000 < (WCAP_INSTR + 1) ** 1000
  B_WALK = WCAP_INSTR - int(SETUP_UNITS) * C - R_STOP_INSTR
  L_THRESH = B_WALK - A_MAX
  assert 0 < L_THRESH < B_WALK < 2 ** 192
  walk_need = (cd * (S * (1 + eps) + E) + WALK_NONSTEP + SWITCH_UNITS) * C
  assert walk_need <= L_THRESH
  budget_units = Decimal(L_THRESH) / C
  slack_frac = (budget_units - WALK_NONSTEP - SWITCH_UNITS) / (cd * (S * (1 + eps) + E)) - 1
  assert slack_frac > Decimal("0.000125") and WCAP_INSTR < 2 ** Decimal("138.184")
  assert lg(SETUP_ITEMS) < 52.2098 and lg(table_units) < 15.07 and eps <= Decimal("0.023793")
  EXP_REPORTS = Decimal(B_WALK) / C / cd * Decimal(2) ** (-THETA_BITS)
  assert markov_ratio_report >= 40000 and markov_ratio_restart >= 39000 and markov_ratio_pairs >= 8 * 10 ** 6
  assert REPORT_CAP / EXP_REPORTS > 2 and lg(MU_REPORTS) < 39.71 and EXP_PAIRS < Decimal(2) ** -19
  assert lg(EXP_RESTARTS) < 4.75 and lg(EXP_STUCK) < -0.68 and lg(EXP_LONG) < 4.71
  assert cyc["rate_useless_walks_longer_than_64_log2"] < -122.6 and cyc["sum_P_l_3_to_64"] < 3.45e-5
  assert lg(REPLAY_UNITS) < 103.52 and lg(R_STOP) < 108.52 and lg(WALK_NONSTEP) < 65.11 and lg(SWITCH_UNITS) < 57.651
  assert lg(E) < 111.81 and lg(D(ZERO_FILL_PER_BYTE * TRIE_BYTES) / C) < 45.9
  assert lg(RESTART_CAP * STUCK) < 110.81 and lg(cd * RESTART_CAP * STUCK) < 109 and lg(restart_units) < 11.04
  assert (-Decimal(STUCK_M)).exp() < Decimal("6.92e-13") and lg(Decimal(RESTART_CAP * METER_UPDATE)) < 25.33
  assert lg(SWITCHES) < 62.71 and lg(PAIR_UNITS) < 104.52
  assert SETUP_ITEMS + Decimal(4101 * 327680) / C < SETUP_UNITS and lg(SETUP_UNITS - SETUP_ITEMS) > 39.5
  # success chain
  SAMPLER_LOSS = Decimal(2) ** Decimal("-71.99")
  DET = Decimal("1.94e-11")
  assert STUCK_M * (-Decimal(STUCK_M)).exp() / (1 - (-Decimal(STUCK_M)).exp()) < DET
  assert lg(N_MAX * REPLAY_ITERS / D(N)) < -33.7 and 2 ** -33.7 < 1e-6
  assert SAMPLES_MAX * Decimal(2) ** -SAMPLE_R_MAX <= SAMPLER_LOSS
  TIE = Decimal(2) ** Decimal("-240")
  assert 70 * N_MAX / p < TIE
  IDENT = Decimal(2) ** -128
  assert (N_MAX + SAMPLES_MAX) / r < IDENT
  LOSS_ITEMS = {"restart cap worst case, as a fraction of the iteration maximum": RESTART_FRAC,
         **{k + " (Markov: expectation / cap)": e / cap for k, (e, cap) in CAPS.items()},
         "detection failure on the two colliding trails": 2 * DET,
         "collision inside a single trail": Decimal("1e-6"),
         "sampler exhaustion allowance": SAMPLER_LOSS, "key tie or zero key": TIE,
         "a start or iterate equal to the identity (zero keys), 1/r each": IDENT,
         "beta^1 = beta^2": Decimal(2) ** -254,
         "the work threshold ending the walk phase early (H1+)": Decimal("0.01")}
  LOSS_SUM = sum(LOSS_ITEMS.values())
  MARKOV_SUM = sum(e / cap for e, cap in CAPS.values())
  REST_LOSS = LOSS_SUM - Decimal("0.01") - RESTART_FRAC - MARKOV_SUM
  assert lg(N_MAX) < 125.710 and RESTART_FRAC < Decimal("3.27e-5") and LOSS_SUM < Decimal("0.0100843")
  assert MARKOV_SUM < Decimal("5.07e-5") and REST_LOSS < Decimal("1.1e-6")
  SUCCESS_LB = (Pm - LOSS_SUM).quantize(Decimal("0.00001"), rounding="ROUND_FLOOR")
  assert SUCCESS_LB >= Decimal("0.3904") and Pm - LOSS_SUM >= SUCCESS_LB
  # memory
  LOG_BYTES = REPORT_CAP * SLOT_BYTES
  MEM = LOG_BYTES + TRIE_BYTES + WALKS * 4096 + 2 ** 26
  # Memory is unscored; what binds is addressability: every array and every byte offset inside it fits one
  # unsigned 64-bit register (each array below 2^63, the total below 2^64).
  assert LOG_BYTES < 2 ** 63 and NODES_LEAF * LEAF_BYTES < 2 ** 63 and TRIE_BYTES < 2 ** 63 and MEM < 2 ** 64
  assert math.log2(NODES_INT) < 42.929 and math.log2(NODES_LEAF) < 51.51 and math.log2(MEM) < 63.284
  assert math.log2(TRIE_BYTES) < 62.855 and math.log2(NODES_LEAF * LEAF_BYTES) < 62.85
  # bit layout of K (stored, below p/2 < 2^380)
  assert F_SHIFT + F_BITS <= 380
  # sensitivity
  def W_with(cc=cd, ee=eps, ss=S):
    return cc * ss * (1 + ee) + cc * E + NONWALK + SWITCH_UNITS
  return {
    "every_step_items": EVERY_ITEMS, "every_step_instr": EVERY, "flush_items": FLUSH_ITEMS, "flush_instr": FLUSH,
    "escape_items": ESC_ITEMS, "escape_instr": ESC_INSTR,
    "step_prices_instr": {"root": STEP_ROOT_I, "lazy": STEP_LAZY_I, "lazy_with_flush": STEP_FLUSH_I},
    "block_instr": BLOCK_I, "c_units_per_step": float(c),
    "R": Rt, "trie_depth_L": Lt, "sum_q2_conservative": float(sumq2),
    "log2_N": lg(D(N)), "log2_N_effective": lg(Npd), "P_model": float(Pm), "log2_S": lg(S),
    "eps_two_cycles_model": float(eps2), "eps_cycles_3_to_64_model": float(eps_long),
    "eps_charged_1.25x": float(eps), "escape_in_step_equivalents": float(esc_steps),
    "useless_walks_longer_than_64_rate_charged_log2": lg(LONG_RATE),
    "hard_max_walk_steps_log2": lg(N_MAX), "expected_stuck_restarts_log2": lg(EXP_STUCK),
    "expected_long_cycle_restarts_log2": lg(EXP_LONG), "restart_cap_log2": lg(RESTART_CAP),
    "restart_worst_case_step_fraction": float(RESTART_FRAC),
    "trie_nodes_internal_log2": math.log2(NODES_INT), "trie_leaves_log2": math.log2(NODES_LEAF),
    "trie_bytes_log2": math.log2(TRIE_BYTES), "node_items": node_items, "node_instr": NODE_INSTR,
    "trie_units_log2": lg(trie_units), "setup_items_units_log2": lg(SETUP_ITEMS),
    "preprocessing_bound_units_log2": float(lg(SETUP_UNITS)),
    "report_log_bytes_log2": math.log2(LOG_BYTES), "peak_memory_bytes_log2": math.log2(MEM),
    "expected_reports_log2": lg(EXP_REPORTS), "report_cap_over_expected": float(REPORT_CAP / EXP_REPORTS),
    "caps_markov": {k: {"expectation_log2": lg(e), "cap_log2": lg(cap), "markov_ratio": float(cap / e),
              "markov_failure_bound": float(e / cap), "loss_item": k + " (Markov: expectation / cap)"}
            for k, (e, cap) in CAPS.items()},
    "walk_nonstep_units_log2": lg(WALK_NONSTEP), "switch_units_log2": lg(SWITCH_UNITS),
    "replay_units_per_trail_log2": lg(REPLAY_UNITS), "stopping_allowance_units_log2": lg(R_STOP),
    "nonwalk_units_log2": lg(NONWALK), "detection_window_steps_log2": lg(E),
    "log2_W": lg(W), "submitted_bound": str(CLAIM),
    "walk_budget_slack_fraction_over_need": float(slack_frac),
    "walk_budget_instr_log2": math.log2(B_WALK), "threshold_instr_log2": math.log2(L_THRESH),
    "success_loss_items": {k: float(v) for k, v in LOSS_ITEMS.items()},
    "success_loss_sum": float(LOSS_SUM), "success_after_every_loss": float(Pm - LOSS_SUM),
    "success_loss_markov_shares_sum": float(MARKOV_SUM), "success_loss_rest_after_h1plus_restart_markov": float(REST_LOSS),
    "success_lower_bound_floor_5_digits": str(SUCCESS_LB),
    "expected_reports_at_iteration_max_log2": lg(MU_REPORTS), "expected_pairs_log2": lg(EXP_PAIRS),
    "sensitivity_bits_eps_margin_1_instead_of_1.25": lg(W_with(ee=eps / Decimal("1.25"))) - lg(W),
    "sensitivity_bits_per_100_instr_on_every_step": lg(W_with(cc=cd + Decimal(100) / C)) - lg(W),
    "sensitivity_bits_if_lazy_kernel_doubled": lg(W_with(cc=cd + Decimal(4 * LAZY_CORE) / Decimal(6 * C))) - lg(W),
  }

out = {}
_P = Tower(p)
_gens = [rand_gt(_P, rng) for _ in range(2)]
_gens += [_P.mul(_gens[0], _gens[1]), _P.cyc_sqr(_gens[1])]
out["partA_production_BLS12_381"] = part_A(rng, PROD_CHECKS)
out["partS_bounded_sampler"] = part_S()
_KN = Kernels()
out["partK_kernels_in_emulator"] = part_K(rng, _KN, _P, _gens, KERNEL_CASES)
if PROD_STEPS > 0:
  out["partC_production_lazy_walk"] = part_C(rng, _KN, _P, _gens[0], PROD_STEPS, PROD_DP_BITS)
if STAT_CLASSES > 0:
  out["partG_key_bit_statistics"] = part_G(rng, _P, _gens[0], STAT_CLASSES)
out["partB_toy"] = part_B(rng, Q_TOY, TRIALS, R, BRANCH_BITS, SQ_CUT, DP_BITS, WALKS)
if PAIR_SAMPLES > 0:
  out["partE_pair_collision_law"] = part_E(rng, Q_TOY, PAIR_SAMPLES, R, BRANCH_BITS, SQ_CUT)
_cyc = runpy.run_path(locate("cycles.py"))["part_F"](CYCLE_M0)
out["partF_fruitless_closed_walks"] = {k: v for k, v in _cyc.items() if k != "P_l_by_length_3_to_64"}
out["partR_reduction_and_G2_output"] = runpy.run_path(locate("reduction.py"),
                          init_globals={"p": p, "r": r, "Tower": Tower})["run"](rng)
assert all(v for k, v in out["partR_reduction_and_G2_output"].items() if k.endswith(("_curve", "_infinity", "equals_z", "not_1",
                                          "nonidentity", "counts_equal")))
out["partD_ledger"] = ledger(_cyc, out["partR_reduction_and_G2_output"])
out["nzcv_flags_check"] = runpy.run_path(locate("flags_check.py"))["run"]()
out["seconds"] = round(time.time() - T0, 2)
out["scope"] = ("Parts A, K, C and G run on real BLS12-381 GT elements (C: the production rule with R = 32, L = 6, "
        "DP field at a test density); part R on BLS12-381 G1, G2 and GT; part B is a toy field; none "
        "validates production costs or production success.")
print(json.dumps(out, indent=1))
