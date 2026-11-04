import os

SP = 0x7000000

SITE_KERNELS = """
\t.macro ENTER
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstp xzr,x30,[sp,#-16]!
\t.endm
\t.macro LEAVE
\tldp xzr,x30,[sp],#16
\tret
\t.endm
t_shx:
\tENTER
\tSHX x5,x25,x28,x26,x27
\tLEAVE
t_shxt:
\tENTER
\tSHXT x25,x28,x26,x27
\tLEAVE
t_fdx:
\tENTER
\tFDX x10,x11,x7,x24,x10,x11,x5,x6,x12,x13
\tLEAVE
t_cpf:
\tENTER
\tAD4 x13,x14,x15,x16,x17,x18,x9,x10,x11,x12
\tFDX x13,x14,x17,x18,x13,x14,x19,x17,x15,x16
\tLEAVE
t_cpk:
\tENTER
\tCPK 0,16,32,48,64,80,96,112,128,144,x3,32,48
\tLEAVE
t_squ:
\tENTER
\tSQ13 x17,x18,x19,x20,x21,x22,x23,x24,x5,x6,x7,x8,x9
\tF5 x21,x22,x23,x24,x8,x21,x22,x9,x8
\tLEAVE
t_sqv:
\tENTER
\tmov x30,x1
\tSQ13 x27,x28,x29,x30,x25,x26,x0,x2,x5,x6,x7,x8,x9
\tF5 x25,x26,x0,x2,x8,x25,x26,x9,x8
\tLEAVE
t_fs:
\tENTER
\tFS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
\tLEAVE
t_faddw:
\tENTER
\tFADDW x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x5
\tLEAVE
t_blk:
\tENTER
\tmov x15,#0xfff8000000000000
\tBLK0 x21,x22,x23,x24,x29,x27,x7,x1
\tLEAVE
"""

COMMON = dict(t_shx=4, t_shxt=3, t_fdx=6, t_cpf=6 + 6, t_squ=58 + 5, t_sqv=1 + 58 + 5, t_fs=7, t_faddw=7, t_blk=1 + 51)
EXTRA = dict(t_shx=3, t_shxt=3, t_fdx=3, t_cpf=3, t_squ=3, t_sqv=3, t_fs=4, t_faddw=4, t_blk=5)

def run_site(E, prog, labels, entry, regs, mem=None):
 mem = {} if mem is None else mem
 out = {}
 n = E.a64.a64_run(prog, labels, entry, dict(regs, sp=SP, x30=E.SENTINEL), mem, out_regs=out)
 return out, n, mem[SP - 16]

def site_checks(E, rng, fprog, flabels):
 p, C0, M64, put, get = E.p, E.C0, E.M64, E.put, E.get
 reg_int, as_regs = E.reg_int, E.as_regs
 prog, labels = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read() + SITE_KERNELS)
 res, fix = {}, {}

 def site(entry, regs, mem=None):
  return run_site(E, prog, labels, entry, regs, mem)

 def tally(name, cases):
  res[name] = dict(cases=len(cases), predicate_holds=sum(c[0] for c in cases),
       deviates_iff_predicate=all(c[0] == c[1] for c in cases), model=all(c[2] for c in cases))

 def fixed(entry, f, n, word, in_set):
  d = fix.setdefault(entry, dict(cases=0, fixups=0, word_ok=True, count_ok=True, word_set_ok=True))
  d["cases"] += 1
  d["fixups"] += f
  d["word_ok"] &= word == C0 * f
  d["count_ok"] &= n == 5 + COMMON[entry] + EXTRA[entry] * f
  d["word_set_ok"] &= in_set or not f

 def fold(v):
  return (v % 2 ** 256) + (v >> 256) * C0

 def low_carry(v, k):
  return (v % 2 ** (64 * k)) + (v >> 256) * C0 >= 2 ** (64 * k)

 cases = []
 for kind in (0, 1, 2, 3):
  for _ in range(48):
   b3 = rng.randrange(2 ** 40, 2 ** 64)
   hi = (b3 * C0) >> 64
   low = (2 ** 192 - hi + rng.randrange(hi), 2 ** 192 - hi - 1 - rng.randrange(2 ** 20), rng.getrandbits(192),
       rng.getrandbits(126) << 64 | (2 ** 64 - hi + rng.randrange(hi)))[kind]
   B = low + (b3 << 192)
   o, n, word = site("t_shx", as_regs(B, ("x25", "x26", "x27", "x28")))
   o = reg_int(o, ("x5", "x25", "x26", "x27"))
   want = B * 2 ** 64 - b3 * p
   cases.append((low + hi >= 2 ** 192, o != want, o == want % 2 ** 256))
   fixed("t_shx", (B % 2 ** 64) + hi >= 2 ** 64, n, word, B % 2 ** 64 > 2 ** 64 - C0)
   o, n, word = site("t_shxt", as_regs(B, ("x25", "x26", "x27", "x28")))
   o = reg_int(o, ("x25", "x26", "x27"))
   cases.append((low + hi >= 2 ** 192, o != want >> 64, o == want % 2 ** 256 >> 64))
   fixed("t_shxt", (B % 2 ** 64) + hi >= 2 ** 64, n, word, B % 2 ** 64 > 2 ** 64 - C0)
 tally("row_shift", cases)
 # CPK (E, F): sum_k a_k * B_k plus K added after the first row, bit-exact on random and extreme rows, words and K
 i0 = labels["t_cpk"][0]
 i1 = next(i for i in range(i0, len(prog)) if prog[i][0] == "ret")
 straight = all(not (m.startswith("b") or m.startswith("cb") or m.startswith("tb")) for m, _ in prog[i0:i1])
 cpk_ok, cpk_n = straight, set()
 X3 = 0x200000
 for c in range(64):
  ext = c % 4
  A = [rng.getrandbits(64) if ext == 0 else 2 ** 64 - 1 for _ in range(4)]
  B = [rng.randrange(p) if ext in (0, 1) else (p - 1 if ext == 2 else 2 ** 256 - 1) for _ in range(4)]
  K = rng.randrange(p) if ext == 0 else (p - 1 if ext in (1, 2) else 2 ** 256 - 1)
  mem = {}
  for w in range(4):
   mem[X3 + 8 * w] = A[w]
   mem[X3 + 32 + 8 * w] = K >> 64 * w & M64
  for k in range(4):
   for w in range(4):
    mem[SP - 16 + 32 + 32 * k + 8 * w] = B[k] >> 64 * w & M64
  o, n, word = site("t_cpk", dict(x3=X3), mem)
  got = reg_int(o, ("x13", "x14", "x15", "x16", "x17", "x18"))
  cpk_ok &= got == sum(A[k] * B[k] for k in range(4)) + K and word == 0
  cpk_n.add(n)
 res["cpk_exact"] = dict(cases=64, exact=cpk_ok, straight_line=straight, counts=sorted(cpk_n))
 for name, entry, top in (("fold", "t_fdx", 2 ** 66), ("fold_of_a_constant_plus_tabled_rows", "t_cpf", 2 ** 66 + 1)):
  cases = []
  for kind in (0, 1, 2, 3):
   for _ in range(48):
    u = rng.randrange(1 + (entry == "t_cpf"), top + 1)
    m = u * C0
    low = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), rng.getrandbits(256),
        rng.getrandbits(127) << 128 | (2 ** 128 - m + rng.randrange(m)))[kind]
    Z = low + (u << 256)
    if entry == "t_fdx":
     regs = as_regs(low, ("x10", "x11", "x12", "x13"))
     regs.update(x7=u & M64, x24=u >> 64)
     outs = ("x10", "x11", "x12", "x13")
    else:
     K = rng.randrange(p)
     regs = as_regs(Z - K, ("x13", "x14", "x15", "x16", "x17", "x18"))
     regs.update(as_regs(K, ("x9", "x10", "x11", "x12")))
     outs = ("x13", "x14", "x15", "x16")
    o, n, word = site(entry, regs)
    o = reg_int(o, outs)
    t = fold(Z)
    cases.append((t >= 2 ** 256, o != t, o == t % 2 ** 256))
    fixed(entry, low_carry(Z, 2), n, word, t % 2 ** 128 < u * C0 and t >= 2 ** 128)
  tally(name, cases)
 for name, entry, a_regs, o_regs in (
   ("square_fold_u", "t_squ", ("x17", "x18", "x19", "x20"), ("x21", "x22", "x23", "x24")),
   ("square_fold_v", "t_sqv", ("x27", "x28", "x29", "x1"), ("x25", "x26", "x0", "x2"))):
  cases = []
  for kind in (0, 1, 2):
   for _ in range(32):
    while True:
     T = (rng.randrange(2 ** 64), rng.randrange(2 ** 65, 2 ** 80), rng.randrange(p))[kind]
     a = E.sqrt_mod_p(T)
     if a is not None:
      break
    a += rng.choice((0, p)) if a < 2 ** 256 - p else 0
    o, n, word = site(entry, as_regs(a, a_regs))
    o = reg_int(o, o_regs)
    V = fold(a * a)
    t = fold(V)
    cases.append((t >= p, o != T, o == t % 2 ** 256))
    fixed(entry, low_carry(V, 2), n, word, t % 2 ** 128 < 2 ** 65)
  tally(name, cases)
 edge = [0, 1, p - 1, C0, 2 ** 255]
 pairs = [(a, b) for a in edge for b in edge] + [(rng.randrange(p), rng.randrange(p)) for _ in range(48)]
 for _ in range(24):
  a = rng.randrange(p >> 1)
  pairs.append((a, a + (rng.randrange(1, 2 ** 190) << 64) - rng.randrange(C0)))
  b = rng.randrange(p >> 1, p)
  pairs.append((2 ** 256 - b + (rng.randrange(1, 2 ** 120) << 64) + 2 ** 64 - 1 - rng.randrange(C0), b))
 fs = fa = True
 for a, b in pairs:
  if a >= p or b >= p:
   continue
  regs = as_regs(a, ("x5", "x6", "x7", "x8"))
  regs.update(as_regs(b, ("x19", "x20", "x21", "x22")))
  o, n, word = site("t_fs", regs)
  o = reg_int(o, ("x9", "x10", "x11", "x12"))
  fs &= o == (a - b) % p
  fixed("t_fs", a < b and (a - b) % 2 ** 64 < C0, n, word, o % 2 ** 64 >= 2 ** 64 - C0)
  regs = as_regs(a, ("x13", "x14", "x15", "x16"))
  regs.update(as_regs(b, ("x17", "x18", "x19", "x20")))
  o, n, word = site("t_faddw", regs)
  o = reg_int(o, ("x13", "x14", "x15", "x16"))
  fa &= o % p == (a + b) % p and o < 2 ** 256
  fixed("t_faddw", a + b >= 2 ** 256 and (a + b) % 2 ** 64 + C0 >= 2 ** 64, n, word, o % 2 ** 64 < C0)
 # beta window sums (BLK0): beta tables (no fixup at random w) and crafted entries summing to 5*2^256 - 1 (fixup)
 cases = []
 for kind in (0, 1, 2):
  for _ in range(24):
   w = rng.getrandbits(256)
   mem = E.LazyTable([])
   if kind == 2:  # low 128 bits all ones, t = 4: the fixup runs (carry out of limb 1), no carry past 2^256
    es = [rng.randrange(p - p // 7, p) for _ in range(4)]
    rest = 4 * 2 ** 256 + 2 ** 200 + 2 ** 128 - 1 - sum(es)
    es += [rest // 2, rest - rest // 2]
    assert all(0 <= e < p for e in es)
    for i in range(6):
     put(mem, E.bt_addr(w, i), es[i])
    tot = sum(es)
   elif kind:
    es = [rng.randrange(p - p // 7, p) for _ in range(5)]
    es.append(5 * 2 ** 256 - 1 - sum(es))
    assert 0 <= es[5] < p
    for i in range(6):
     put(mem, E.bt_addr(w, i), es[i])
    tot = sum(es)
   else:
    tot = sum(get(mem, E.bt_addr(w, i)) for i in range(6))
   o, n, word = site("t_blk", as_regs(w, ("x21", "x22", "x23", "x24")), mem)
   o = reg_int(o, ("x29", "x27", "x7", "x1"))
   f = int(fold(tot) >> 256)
   fixed("t_blk", int(tot % 2 ** 128 + (tot >> 256) * C0 >= 2 ** 128), n, word, kind >= 1)
   ref = fold(tot) % 2 ** 256 + C0 * f  # the full carry chain: fold, and C0 once more on its carry past 2^256
   cases.append((f == 1, o != fold(tot), o == ref and ref < 2 ** 256 and (kind == 0) == (o % p == E.BETA * w % p)))
 tally("window_sum_fold", cases)
 # beta tables (route.s tbld) on short windows: every entry of the 7 slots and the count 4 + 7*12 + 20*(6n + n5)
 rprog, rlabels = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "route.s")).read())
 tb_ok = True
 for n, n5 in ((7, 5), (12, 3)):
  mem = {}
  for k in range(7):
   put(mem, 0x1000 + 32 * k, E.BETA * pow(2, 43 * E.SLOTW[k], p) % p)
  steps = E.a64.a64_run(rprog, rlabels, "tbld", dict(x0=0x100000, x1=0x1000, x2=n, x3=n5, sp=SP, x30=E.SENTINEL), mem)
  tb_ok &= steps == 4 + 7 * 12 + 20 * (6 * n + n5) and all(
   get(mem, 0x100000 + 256 * t + 32 * k) == E.BETA * t * pow(2, 43 * E.SLOTW[k], p) % p
   for k in range(7) for t in range(n5 if k == 5 else n))
 res["beta_tables"] = dict(cases=2, exact=tb_ok)
 res["fs_canonical_exact"] = dict(cases=len(pairs), exact=fs)
 res["faddw_canonical_exact"] = dict(cases=len(pairs), exact=fa)
 res["fixups"] = fix
 res["fixups_ok"] = all(d["word_ok"] and d["count_ok"] and d["word_set_ok"] and d["fixups"] > 0 for d in fix.values())
 weak = [0, 1, C0, p - 1, p, p + 1, 2 ** 256 - 1] + [rng.getrandbits(256) for _ in range(12)]
 mul_in = [(a, b) for a in weak[:7] for b in weak[:7]] + [(a, rng.getrandbits(256)) for a in weak[7:]]
 for _ in range(12):
  a = rng.randrange(1, p)
  mul_in.append((a, rng.randrange(1, 2 ** 64) * pow(a, -1, p) % p))
 sqr_in = weak + [p - rng.randrange(1, 2 ** 100) for _ in range(12)] + [2 ** 256 - rng.randrange(1, 2 ** 40)
                   for _ in range(6)]
 while len(sqr_in) < len(weak) + 42:
  s_ = rng.randrange(1, 2 ** 64)
  if pow(s_, (p - 1) // 2, p) == 1:
   sqr_in.append(pow(s_, (p + 1) // 4, p))
 ok_m = ok_s = True
 for a, b in mul_in:
  mem = {}
  put(mem, 0x1000, a)
  put(mem, 0x1100, b)
  E.a64.a64_run(fprog, flabels, "fe_mul_f", dict(x0=0x2000, x1=0x1000, x2=0x1100, sp=0x100000, x30=E.SENTINEL), mem)
  ok_m &= get(mem, 0x2000) % p == a * b % p
 for a in sqr_in:
  mem = {}
  put(mem, 0x1000, a)
  E.a64.a64_run(fprog, flabels, "fe_sqr_f", dict(x0=0x2000, x1=0x1000, sp=0x100000, x30=E.SENTINEL), mem)
  ok_s &= get(mem, 0x2000) % p == a * a % p

 def final_carry(v):
  t = (v % 2 ** 256) + (v >> 256) * C0
  return (t % 2 ** 256) + (t >> 256) * C0 >= 2 ** 256
 res["fe_mul_f_exact"] = dict(cases=len(mul_in), exact=ok_m, final_carry_cases=sum(final_carry(a * b) for a, b in mul_in))
 res["fe_sqr_f_exact"] = dict(cases=len(sqr_in), exact=ok_s, final_carry_cases=sum(final_carry(a * a) for a in sqr_in))
 res["all_ok"] = res["fixups_ok"] and all(v.get("deviates_iff_predicate", True) and v.get("model", True) and v.get("exact", True)
         and v.get("predicate_holds", 1) > 0 for k, v in res.items()
         if isinstance(v, dict) and k != "fixups")
 return res
