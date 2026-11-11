"""Parts H and I of lazy.py (executed by lazy.py in its own namespace: rng, Tower, profile_to_tower, fwords,
row, f_off, build_ttloop, run_tt, lazy_env, lazy_out, run and the constants come from lazy.py)."""

# =============================================================== Part H: toy mean collision time (8-adding walk)
import math as _m


def _is_prime(n):
 if n < 2:
  return False
 for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
  if n % q == 0:
   return n == q
 d, s = n - 1, 0
 while d % 2 == 0:
  d //= 2
  s += 1
 for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
  x = pow(a, d, n)
  if x in (1, n - 1):
   continue
  for _ in range(s - 1):
   x = x * x % n
   if x == n - 1:
    break
  else:
   return False
 return True


def toy_group(bits, coins):
 """a fresh 61-bit prime q = 2 n h + 1 with prime n of the given size; generator of the order-n subgroup"""
 while True:
  n = coins.getrandbits(bits) | (1 << (bits - 1)) | 1
  if not _is_prime(n):
   continue
  for _ in range(200):
   h = coins.getrandbits(60 - bits) | (1 << (59 - bits))
   q = 2 * n * h + 1
   if _is_prime(q):
    while True:
     g = pow(coins.randrange(2, q - 1), 2 * h, q)
     if g != 1:
      return q, n, g


def toy_run(bits, dbits, coins, m=8):
 """one target, the walk of proof.md section 3 at toy scale: s = floor(m t / q), distinguished iff
    t 2^dbits < q, sequential walks to their distinguished point, run to the first recovery.
    Returns (touch steps, total steps to recovery, walks, n)."""
 q, n, g = toy_group(bits, coins)
 x = coins.randrange(1, n)
 z = pow(g, x, q)
 cs = [coins.randrange(n) for _ in range(m)]
 ds = [coins.randrange(n) for _ in range(m)]
 Ms = [pow(g, cs[i], q) * pow(z, ds[i], q) % q for i in range(m)]
 lim = q >> dbits
 store = {}
 seen = set()
 steps, touch, walks = 0, None, 0
 len_cap = 64 << dbits
 while True:
  walks += 1
  a, b = coins.randrange(n), coins.randrange(n)
  y = pow(g, a, q) * pow(z, b, q) % q
  trail = []
  k = 0
  while y >= lim and k < len_cap:
   if touch is None:
    trail.append(y)
    if y in seen:
     touch = steps
   s = y * m // q
   y = y * Ms[s] % q
   a = (a + cs[s]) % n
   b = (b + ds[s]) % n
   steps += 1
   k += 1
  if touch is None:
   seen.update(trail)
   if y in seen:
    touch = steps
   seen.add(y)
  if k >= len_cap:
   continue
  old = store.get(y)
  if old is not None and old[1] != b:
   xr = (old[0] - a) * pow(b - old[1], -1, n) % n
   assert pow(g, xr, q) == z
   return (touch if touch is not None else steps), steps, walks, n
  store[y] = (a, b)


def part_H(settings, coins, budget):
 """mean of touch/sqrt(pi kappa n/2) and of recovery/sqrt(pi kappa n/2), kappa = 8/7, per setting"""
 t0 = time.time()
 out = []
 kappa = 8 / 7
 for idx, (bits, dbits, trials) in enumerate(settings):
  RT, RS, Wk = [], [], []
  deadline = t0 + budget * (idx + 1) / len(settings)       # equal share of the budget per setting
  while len(RT) < trials and time.time() < deadline:
   a, b, c, n = toy_run(bits, dbits, coins)
   unit = (_m.pi * kappa * n / 2) ** 0.5
   RT.append(a / unit)
   RS.append(b / unit)
   Wk.append(c)
  if not RT:
   continue
  k = len(RT)
  mt, ms = sum(RT) / k, sum(RS) / k
  st = (sum((v - mt) ** 2 for v in RT) / max(1, k - 1)) ** 0.5 / k ** 0.5
  ss = (sum((v - ms) ** 2 for v in RS) / max(1, k - 1)) ** 0.5 / k ** 0.5
  out.append({"n_bits": bits, "dp_bits": dbits, "runs": k,
     "mean_touch_over_model": round(mt, 4), "se": round(st, 4),
     "mean_recovery_steps_over_model": round(ms, 4), "se_recovery": round(ss, 4),
     "mean_walks": round(sum(Wk) / k, 2)})
 return out

# =============================================================== Part S: generated tag-tracing step loop
# Entry layout (2,624 bytes; the pointer to an entry is P = entry + 512):
#   P - 512 .. P - 1   : F words of limbs j = 0..31   (Flo_j at P - 512 + 16 j, Fhi_j at +8)
#   P + 0   .. P + 63  : the eight child pointers (each again of the form entry' + 512; 0 at window depth)
#   P + 64  .. P + 575 : F words of limbs j = 32..47  (P + 64 + 16 (j - 32))
#   P + 576 ..         : the 54 window-product operands of T_e and the labels (not read by the loop)
# Registers inside tt_loop: x0 = P (current entry), x1 = y (48 limbs of the window base, canonical
# Montgomery words), x2 = trip budget, x3 = A, x4 = umulh temporary, x5/x6 = F pair, x7/x8 = reloaded limb
# pair, x9 = selector / scratch, x10 = child pointer, x11..x17 and x19..x30 = 20 resident limbs.
NRES = 20
RES = [11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 9]
MASK61 = (1 << 61) - 1


def f_off(j):
 return -512 + 16 * j if j < 32 else 64 + 16 * (j - 32)


def build_ttloop():
 """returns (prologue, body, exits) as instruction lists; body is one trip, ending in b.ne to its start"""
 res = RES[:NRES - 1] + [None]          # x9 is needed as scratch: keep 19 resident limbs in registers
 nres = NRES - 1
 pro = [("addi", SP, SP, -112), ("str", 3, SP, 96)]          # x3 = cursor pointer on entry
 for k in range(6):
  pro.append(("stp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k))
 for j in range(0, nres - 1, 2):
  pro.append(("ldp", res[j], res[j + 1], 1, 8 * j))
 if nres % 2:
  pro.append(("ldr", res[nres - 1], 1, 8 * (nres - 1)))
 body = []
 first = True
 j = 0
 while j < 48:
  if j < nres:
   ys = [res[j]]
  else:
   if j + 1 < 48 and j + 1 >= nres:
    body.append(("ldp", 7, 8, 1, 8 * j))
    ys = [7, 8]
   else:
    body.append(("ldr", 7, 1, 8 * j))
    ys = [7]
  for k, yr in enumerate(ys):
   jj = j + k
   body.append(("ldp", 5, 6, 0, f_off(jj)))
   body.append(("umulh", 4, yr, 5))
   if first:
    body.append(("mul", 3, yr, 6))
    first = False
   else:
    body.append(("madd", 3, yr, 6, 3))
   body.append(("add", 3, 3, 4))
  j += len(ys)
 body += [("addi", 9, 3, 96),                 # t = A + 96
    ("andi", 9, 9, MASK61),             # t mod 2^61
    ("subs", ZR, 9, ("#", 96)),         # certified iff t > 96
    ("b.ls", "EXACT"),
    ("lsr", 9, 3, 61),                  # s = A >> 61
    ("ldrx", 10, 0, 9),                 # n = [P + 8 s]
    ("cbz", 10, "REBASE"),
    ("mov", 0, 10),
    ("subsi", 2, 2, 1),
    ("b.ne", "LOOP")]
 # every exit: store the cursor (e, A, budget), status, restore, return
 def exit_path(code):
  return [("ldr", 9, SP, 96), ("stp", 0, 3, 9, 0), ("str", 2, 9, 16), ("movi", 0, code)] + \
    [("ldp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k) for k in range(6)] + [("addi", SP, SP, 112), ("ret",)]
 exits = {"BUDGET": exit_path(0), "EXACT": exit_path(1), "REBASE": exit_path(2)}
 return pro, body, exits


def validate_tt(code):
 bad = []
 for ins in code:
  op = ins[0]
  if op == "madd":
   ok = all(x in GPR for x in ins[1:5])
  elif op == "andi":
   ok = ins[1] in GPR and ins[2] in GPR and ins[3] == MASK61   # 61 ones: a valid bitmask immediate
  elif op == "ldrx":
   ok = all(x in GPR for x in ins[1:4])                          # LDR Xt, [Xn, Xm, LSL #3]
  elif op in ("b.ls", "b.ne", "cbz"):
   ok = True
  elif op == "mov":
   ok = ins[1] in GPR and ins[2] in GPR
  elif op == "subsi":
   ok = ins[1] in GPR and ins[2] in GPR and 0 <= ins[3] <= 4095
  else:
   ok = not validate([ins])
  if ok and any(isinstance(x, int) and x == 18 for x in ins[1:] if not isinstance(x, tuple)) and op not in ("stp", "str"):
   ok = False
  if not ok:
   bad.append(" ".join(str(x) for x in ins))
 return bad


def run_tt(pro, body, exits, mem, regs):
 """Execute prologue, trips and one exit. Returns (status, trips, instructions per trip list, total)."""
 C = Z = 0
 total = 0
 trips = []

 def ex(ins):
  nonlocal C, Z
  op = ins[0]
  if op == "ldp":
   a = regs[ins[3]] + ins[4]
   regs[ins[1]], regs[ins[2]] = mem[a], mem[a + 8]
  elif op == "ldr":
   regs[ins[1]] = mem[regs[ins[2]] + ins[3]]
  elif op == "ldrx":
   regs[ins[1]] = mem[regs[ins[2]] + 8 * regs[ins[3]]]
  elif op == "stp":
   a = regs[ins[3]] + ins[4]
   mem[a], mem[a + 8] = regs[ins[1]], regs[ins[2]]
  elif op == "str":
   mem[regs[ins[2]] + ins[3]] = regs[ins[1]]
  elif op == "umulh":
   regs[ins[1]] = (regs[ins[2]] * regs[ins[3]]) >> 64
  elif op == "mul":
   regs[ins[1]] = (regs[ins[2]] * regs[ins[3]]) & M64
  elif op == "madd":
   regs[ins[1]] = (regs[ins[2]] * regs[ins[3]] + regs[ins[4]]) & M64
  elif op == "add":
   regs[ins[1]] = (regs[ins[2]] + regs[ins[3]]) & M64
  elif op == "addi":
   regs[ins[1]] = (regs[ins[2]] + ins[3]) & M64
  elif op == "andi":
   regs[ins[1]] = regs[ins[2]] & ins[3]
  elif op == "subs":
   d = regs[ins[2]] - ins[3][1]
   C, Z = int(d >= 0), int(d == 0)
   if ins[1] != ZR:
    regs[ins[1]] = d & M64
  elif op == "subsi":
   d = regs[ins[2]] - ins[3]
   C, Z = int(d >= 0), int(d == 0)
   regs[ins[1]] = d & M64
  elif op == "lsr":
   regs[ins[1]] = regs[ins[2]] >> ins[3]
  elif op == "mov":
   regs[ins[1]] = regs[ins[2]]
  elif op == "movi":
   regs[ins[1]] = ins[2]
  elif op == "ret":
   pass
  else:
   raise RuntimeError(op)
  regs[ZR] = 0

 for ins in pro:
  ex(ins)
  total += 1
 while True:
  n = 0
  target = None
  for ins in body:
   n += 1
   op = ins[0]
   if op == "b.ls":
    if (not C) or Z:
     target = ins[1]
     break
   elif op == "cbz":
    if regs[ins[1]] == 0:
     target = ins[2]
     break
   elif op == "b.ne":
    if Z:
     target = "BUDGET"
     break
   else:
    ex(ins)
  total += n
  trips.append(n)
  if target is not None:
   for ins in exits[target]:
    ex(ins)
    total += 1
   return target, trips, total


def part_S(cases_entries, steps):
 """random 128-bit F words, random child structure: compare A, the selector, the exits and the
    certificate test with the definition of proof.md section 2 at every trip"""
 pro, body, exits = build_ttloop()
 viol = validate_tt(pro + body + [i for e in exits.values() for i in e])
 nent = cases_entries
 base_addr = 0x10000000
 P = [base_addr + 4096 * i + 512 for i in range(nent)]
 mem = {}
 Fw = []
 for i in range(nent):
  F = [(rng.getrandbits(64), rng.getrandbits(64)) for _ in range(48)]
  if i % 7 == 3:
   F[0] = (F[0][0], 0)
  Fw.append(F)
  for j, (lo, hi) in enumerate(F):
   mem[P[i] + f_off(j)] = lo
   mem[P[i] + f_off(j) + 8] = hi
  for s in range(8):
   mem[P[i] + 8 * s] = 0 if i % 11 == 10 else P[rng.randrange(nent)]     # i = 10 mod 11: window depth
 ylimbs = [rng.getrandbits(64) for _ in range(48)]
 ylimbs[47] |= 1                            # odd: lets the test force a chosen A below
 for i in range(nent):
  if i % 13 == 5:                        # force an uncertified tag: (A + 96) mod 2^61 <= 96
   F = Fw[i]
   A = 0
   for j in range(48):
    A = (A + ((ylimbs[j] * F[j][0]) >> 64) + ylimbs[j] * F[j][1]) & M64
   target = (rng.getrandbits(3) << 61) | ((1 << 61) - 96 + rng.randrange(0, 96))   # low 61 bits in [2^61-96, 2^61-1]: (A+96) mod 2^61 <= 95
   lo, hi = F[47]
   hi = (hi + (target - A) * pow(ylimbs[47], -1, 1 << 64)) & M64
   F[47] = (lo, hi)
   mem[P[i] + f_off(47) + 8] = hi
 Y = 0x20000000
 for j in range(48):
  mem[Y + 8 * j] = ylimbs[j]
 CUR = 0x30000000
 mism = 0
 trip_lengths = set()
 statuses = {}
 checked = 0
 starts = list(range(min(nent, 40))) + [i for i in range(nent) if i % 13 == 5]
 for start in starts:
  regs = [rng.getrandbits(64) for _ in range(33)]
  regs[0], regs[1], regs[2], regs[3], regs[SP] = P[start], Y, steps, CUR, 0x40000000
  status, trips, total = run_tt(pro, body, exits, mem, regs)
  statuses[status] = statuses.get(status, 0) + 1
  # replay with the definition
  e = start
  for k in range(len(trips)):
   A = 0
   for j in range(48):
    lo, hi = Fw[e][j]
    A = (A + ((ylimbs[j] * lo) >> 64) + ylimbs[j] * hi) & M64
   cert = ((A + 96) & M64) & MASK61 > 96
   nxt = mem[P[e] + 8 * (A >> 61)] if cert else None
   last = k == len(trips) - 1
   checked += 1
   if not cert:
    mism += not (last and status == "EXACT")
    break
   if nxt == 0:
    mism += not (last and status == "REBASE")
    break
   e = P.index(nxt)
   if last:
    mism += status != "BUDGET"
   else:
    trip_lengths.add(trips[k])
  if mem.get(CUR + 8) != A:
   mism += 1
 return {"encoding_violations": len(viol), "first_violations": viol[:3],
   "trip_instructions_certified_in_window": sorted(trip_lengths),
   "body_length": len(body), "prologue": len(pro),
   "exit_path": max(len(e) for e in exits.values()),
   "entries": nent, "runs": len(starts), "trips_checked": checked, "exit_statuses": statuses,
   "mismatches_with_definition": mism,
   "body_opcodes": {op: sum(1 for i in body if i[0] == op) for op in sorted({i[0] for i in body})}}



# =============================================================== Part I: integrated full-width walk
def part_I(gm, steps, ell=8):
 """The production walk at full BN254 width through BOTH generated programs: each step is one trip of the
    generated loop (part S) on a real table entry (F words from the row of t(Y T_e) in this layout, child
    pointers), each window end is the generated product (part L) on real data. Window length ell (production
    256), 8 real multipliers g^c z^d. Checked at every step against the definition s = floor(8 t(Y)/p)."""
 P = Tower(p)
 g = profile_to_tower(G_PROFILE, p)
 x = rng.randrange(1, r)
 z = P.pow(g, x)
 cs = [rng.randrange(r) for _ in range(8)]
 ds = [rng.randrange(r) for _ in range(8)]
 Ms = [P.mul(P.pow(g, cs[i]), P.pow(z, ds[i])) for i in range(8)]
 Rinv = pow(RMONT, -1, p)
 pro, body, exits = build_ttloop()
 mem = {}
 entries = {}            # key (multiset counts) -> (address P, T_e, alpha, beta)
 POISON = 0x7000000000 + 512
 nxt = [0x10000000]

 def make(key, T, al, be):
  addr = nxt[0] + 512
  nxt[0] += 4096
  for j, (lo, hi) in enumerate(fwords(row(P, T))):
   mem[addr + f_off(j)] = lo
   mem[addr + f_off(j) + 8] = hi
  for s in range(8):
   mem[addr + 8 * s] = 0 if sum(key) == ell else POISON
  entries[key] = (addr, T, al, be)
  return entries[key]

 def child(key, s):
  k2 = key[:s] + (key[s] + 1,) + key[s + 1:]
  if k2 not in entries:
   a0, T, al, be = entries[key]
   make(k2, P.mul(T, Ms[s]), (al + cs[s]) % r, (be + ds[s]) % r)
  mem[entries[key][0] + 8 * s] = entries[k2][0]
  return k2
 root = (0,) * 8
 make(root, P.one, 0, 0)
 a, b = rng.randrange(r), rng.randrange(r)
 Y0 = P.mul(P.pow(g, a), P.pow(z, b))
 YB = 0x20000000
 CUR = 0x30000000

 def put_y(Y):
  for m in range(6):
   for bb in range(2):
    v = Y[m][bb] * RMONT % p
    for i in range(4):
     mem[YB + 8 * (8 * m + 4 * bb + i)] = (v >> (64 * i)) & M64
 put_y(Y0)
 key = root
 trips, rebases, bad, label_checks, prod_fail = set(), 0, 0, 0, 0
 s_hist = []
 for step in range(steps):
  Ycur = P.mul(Y0, entries[key][1])
  tval = (Ycur[0][0] - 9 * Ycur[0][1]) % p
  s_true = (8 * tval) // p
  s_hist.append(s_true)
  depth_end = sum(key) == ell
  if not depth_end:
   k2 = child(key, s_true)
  regs = [0] * 33
  regs[0], regs[1], regs[2], regs[3], regs[SP] = entries[key][0], YB, 1, CUR, 0x40000000
  status, tr, total = run_tt(pro, body, exits, mem, regs)
  A = mem[CUR + 8]
  cert = ((A + 96) & M64) & MASK61 > 96
  bad += cert and (A >> 61) != s_true
  if not cert:
   bad += status != "EXACT"          # exact evaluation: the driver uses s_true from t directly
   k2 = None
  if depth_end and cert:
   bad += status != "REBASE"
  if (not depth_end) and cert:
   bad += status != "BUDGET" or mem[CUR] != entries[k2][0]
   trips.add(tr[-1])
  if depth_end or not cert:
   # re-base by the generated product: Y0 <- Y0 * T_e, labels added; then take the step from the root
   T_e = entries[key][1]
   mem2, regs2 = lazy_env(gm, Y0, T_e)
   run(gm.code, mem2, regs2)
   raw = lazy_out(mem2, regs2)
   newY = tuple(tuple(c_ * Rinv % p for c_ in cc) for cc in raw)
   prod_fail += newY != P.mul(Y0, T_e) or any(c_ >= p for cc in raw for c_ in cc)
   a, b = (a + entries[key][2]) % r, (b + entries[key][3]) % r
   Y0 = newY
   put_y(Y0)
   rebases += 1
   key = child(root, s_true)
  else:
   key = k2
  if step % 100 == 99:
   al, be = entries[key][2], entries[key][3]
   Yrep = P.mul(Y0, entries[key][1])
   label_checks += 1
   bad += Yrep != P.mul(P.pow(g, (a + al) % r), P.pow(z, (b + be) % r))
   bad += Yrep != P.mul(Ycur, Ms[s_true])
 return {"steps": steps, "window_length_in_test": ell, "window_rebases_by_generated_product": rebases,
   "generated_product_failures": prod_fail, "mismatches_with_definition": bad,
   "certified_in_window_trip_instructions": sorted(trips), "label_checks": label_checks,
   "entries_built": len(entries), "selector_histogram": [s_hist.count(i) for i in range(8)]}

