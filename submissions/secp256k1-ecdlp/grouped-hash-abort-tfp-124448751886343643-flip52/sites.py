import os

SP = 0x7000000

SITE_KERNELS = """
\t.macro ENTER
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstr x30,[sp,#520]
\tstr xzr,[sp]
\t.endm
\t.macro LEAVE
\tldr x30,[sp,#520]
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
t_lke:
\tENTER
\tstp x0,x2,[sp,#464]
\tstr x1,[sp,#480]
\tLKE8
\tLEAVE
t_lkf:
\tENTER
\tstp x0,x2,[sp,#464]
\tstr x1,[sp,#480]
\tLKFF8
\tLEAVE
t_lkp:
\tENTER
\tLKP8 352
\tLEAVE
t_c2uv:
\tENTER
\tC2UV8
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
\tmov x5,#0xfff8000000000000
\torr x6,x5,#64
\torr x9,x5,#128
\tBLK6 x21,x22,x23,x24,x29,x27,x7,x1
\tLEAVE
"""

COMMON = dict(t_shx=4, t_shxt=3, t_fdx=6, t_lke=2 + 70, t_lkf=2 + 63,
       t_squ=58 + 5, t_sqv=1 + 58 + 5, t_fs=7, t_faddw=7, t_blk=3 + 46)
EXTRA = dict(t_shx=3, t_shxt=3, t_fdx=3, t_lke=5, t_lkf=5, t_squ=3, t_sqv=3, t_fs=4, t_faddw=4, t_blk=5)

def site_checks(E, rng, fprog, flabels):
 p, C0, M64, put, get = E.p, E.C0, E.M64, E.put, E.get
 reg_int, as_regs = E.reg_int, E.as_regs
 prog, labels = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read() + SITE_KERNELS)
 res, fix = {}, {}
 LS = ([(E.XT0, 1 << 51), (E.BT0, 1 << 64)], [])

 def site(entry, regs, mem=None):
  mem = {} if mem is None else mem
  out = {}
  n = E.a64.a64_run(prog, labels, entry, dict(regs, sp=SP, x30=E.SENTINEL), mem, out_regs=out, ls64=LS)
  return out, n, mem[SP]

 def tally(name, cases):
  res[name] = dict(cases=len(cases), predicate_holds=sum(c[0] for c in cases),
       deviates_iff_predicate=all(c[0] == c[1] for c in cases), model=all(c[2] for c in cases))

 def fixed(entry, f, n, word, in_set):
  d = fix.setdefault(entry, dict(cases=0, fixups=0, word_ok=True, count_ok=True, word_set_ok=True))
  d["cases"] += 1; d["fixups"] += f
  d["word_ok"] &= word == C0 * f
  d["count_ok"] &= n == 6 + COMMON[entry] + EXTRA[entry] * f
  d["word_set_ok"] &= in_set or not f

 fold = lambda v: (v % 2 ** 256) + (v >> 256) * C0
 low_carry = lambda v, k: (v % 2 ** (64 * k)) + (v >> 256) * C0 >= 2 ** (64 * k)

 # 1. SHX / SHXT
 cases = []
 for kind in (0, 1, 2, 3):
  for _ in range(48):
   b3 = rng.randrange(2 ** 40, 2 ** 64); hi = (b3 * C0) >> 64
   low = (2 ** 192 - hi + rng.randrange(hi), 2 ** 192 - hi - 1 - rng.randrange(2 ** 20), rng.getrandbits(192),
       rng.getrandbits(126) << 64 | (2 ** 64 - hi + rng.randrange(hi)))[kind]
   B = low + (b3 << 192)
   o, n, word = site("t_shx", as_regs(B, ("x25", "x26", "x27", "x28")))
   ov = reg_int(o, ("x5", "x25", "x26", "x27")); want = B * 2 ** 64 - b3 * p
   cases.append((low + hi >= 2 ** 192, ov != want, ov == want % 2 ** 256))
   fixed("t_shx", (B % 2 ** 64) + hi >= 2 ** 64, n, word, B % 2 ** 64 > 2 ** 64 - C0)
   o, n, word = site("t_shxt", as_regs(B, ("x25", "x26", "x27", "x28")))
   ov = reg_int(o, ("x25", "x26", "x27"))
   cases.append((low + hi >= 2 ** 192, ov != want >> 64, ov == want % 2 ** 256 >> 64))
   fixed("t_shxt", (B % 2 ** 64) + hi >= 2 ** 64, n, word, B % 2 ** 64 > 2 ** 64 - C0)
 tally("row_shift", cases)

 # 2. 8x32b LD64B lookup folds (LKE8: 9 terms, LKFF8: 8 terms, LKP8: 9 terms)
 cases = []
 for entry, o0, tbo, ci in (("t_lke", 272, 224, 0), ("t_lkf", 336, 232, 1)):
  for kind in (0, 0, 1, 2):
   for _ in range(12):
    a, tb, X3 = rng.randrange(1, p), E.XT0 + (ci << 40), 0x200000
    mem = E.LazyTable([], [(a, rng.randrange(1, p))])
    mem[X3 + tbo] = tb
    w = rng.getrandbits(256)
    offs = E.offs8(w)
    for k in range(8): mem[SP + o0 + 8 * k] = offs[k]
    b4 = rng.getrandbits(256) if entry == "t_lke" else 0
    if kind:
     b4 = (p - 1 - rng.randrange(100) if kind == 1 else rng.randrange(2 ** 200)) if entry == "t_lke" else 0
     T = ((8 if entry == "t_lke" else 7) * 2 ** 256 - 1 if kind == 1
          else 6 * 2 ** 256 + (rng.getrandbits(100) << 128) + 2 ** 128 - 1 - rng.randrange(C0))
     es = [p - 1 - rng.randrange(100) for _ in range(6 if kind == 1 else 5)]
     rest = T - b4 - sum(es)
     es += [rest // 2, rest - rest // 2] if kind == 1 else [rest // 3, rest // 3, rest - 2 * (rest // 3)]
     assert len(es) == 8 and all(0 <= e < p for e in es)
     for k in range(8): put(mem, tb + ((w >> 32 * k & 0xffffffff) << 8) + 32 * k, es[k])
    put(mem, SP + 16, b4)
    o, n, word = site(entry, dict(x3=X3, x0=11, x2=22, x1=33), mem)
    ov = reg_int(o, ("x13", "x14", "x9", "x10"))
    tot = sum(get(mem, tb + ((w >> 32 * k & 0xffffffff) << 8) + 32 * k) for k in range(8)) + b4
    t = tot >> 256
    want = (E.x_consts(a)[ci] * w + (b4 if ci == 0 else E.x_consts(a)[3])) % p
    cases.append((fold(tot) >= 2 ** 256, ov != fold(tot) % 2 ** 256,
        ov == E.lkfold(tot) and ov % p == tot % p and (kind > 0 or ov % p == want) and (o["x0"], o["x2"], o["x1"]) == (11, 22, 33)))
    fixed(entry, tot % 2 ** 128 + t * C0 >= 2 ** 128, n, word, tot % 2 ** 128 >= 2 ** 128 - t * C0)
 tally("lookup_sum_fold", cases)

 # 3. FDX
 cases = []
 for kind in (0, 1, 2, 3):
  for _ in range(48):
   u = rng.randrange(1, 2 ** 66 + 1); m = u * C0
   low = (2 ** 256 - m + rng.randrange(m), 2 ** 256 - m - 1 - rng.randrange(2 ** 20), rng.getrandbits(256),
       rng.getrandbits(127) << 128 | (2 ** 128 - m + rng.randrange(m)))[kind]
   Z = low + (u << 256)
   regs = as_regs(low, ("x10", "x11", "x12", "x13")); regs.update(x7=u & M64, x24=u >> 64)
   o, n, word = site("t_fdx", regs); ov = reg_int(o, ("x10", "x11", "x12", "x13")); t = fold(Z)
   cases.append((t >= 2 ** 256, ov != t, ov == t % 2 ** 256))
   fixed("t_fdx", low_carry(Z, 2), n, word, t % 2 ** 128 < u * C0 and t >= 2 ** 128)
 tally("fold", cases)

 # 4. F5 square folds
 for name, entry, a_regs, o_regs in (
   ("square_fold_u", "t_squ", ("x17", "x18", "x19", "x20"), ("x21", "x22", "x23", "x24")),
   ("square_fold_v", "t_sqv", ("x27", "x28", "x29", "x1"), ("x25", "x26", "x0", "x2"))):
  cases = []
  for kind in (0, 1, 2):
   for _ in range(24):
    while True:
     T = (rng.randrange(2 ** 64), rng.randrange(2 ** 65, 2 ** 80), rng.randrange(p))[kind]
     a = E.sqrt_mod_p(T)
     if a is not None: break
    a += rng.choice((0, p)) if a < 2 ** 256 - p else 0
    o, n, word = site(entry, as_regs(a, a_regs)); ov = reg_int(o, o_regs); V = fold(a * a); t = fold(V)
    cases.append((t >= p, ov != T, ov == t % 2 ** 256))
    fixed(entry, low_carry(V, 2), n, word, t % 2 ** 128 < 2 ** 65)
  tally(name, cases)

 # 5. FS and FADDW (including C2 weak operands < 2^256 testing both exact and second-borrow/second-carry predicate!)
 fs_cases, fa_cases = [], []
 pairs = [(rng.randrange(p), rng.randrange(p)) for _ in range(48)]
 for _ in range(24):
  a = rng.randrange(p >> 1)
  pairs.append((a, a + (rng.randrange(1, 2 ** 190) << 64) - rng.randrange(C0)))
  b = rng.randrange(p >> 1, p)
  pairs.append((2 ** 256 - b + (rng.randrange(1, 2 ** 120) << 64) + 2 ** 64 - 1 - rng.randrange(C0), b))
  # C2 weak error-site cases: L in [p+1, 2^256-1], A < L - p (FS fails) and A, L in [p+1, 2^256-1] with A+L >= 2^256+p (FADDW fails)
  L_bad = p + 1 + rng.randrange(C0 - 1)
  A_bad = rng.randrange(L_bad - p)
  pairs.append((A_bad, L_bad))
  pairs.append((p + C0 - 1 - rng.randrange(10), p + 10 + rng.randrange(C0 - 10)))
 for a, b in pairs:
  if a >= 2 ** 256 or b >= 2 ** 256: continue
  regs = as_regs(a, ("x5", "x6", "x7", "x8")); regs.update(as_regs(b, ("x19", "x20", "x21", "x22")))
  o, n, word = site("t_fs", regs); ov = reg_int(o, ("x9", "x10", "x11", "x12"))
  pred_fs = a - b + p < 0
  want_fs = a - b if a >= b else a - b + p
  fs_cases.append((pred_fs, ov != want_fs % 2 ** 256 or pred_fs, (ov == want_fs) if not pred_fs else ((a - b) % p >= p - C0 + 1)))
  if not pred_fs:
   fixed("t_fs", a < b and (a - b) % 2 ** 64 < C0, n, word, ov % 2 ** 64 >= 2 ** 64 - C0)
  regs = as_regs(a, ("x13", "x14", "x15", "x16")); regs.update(as_regs(b, ("x17", "x18", "x19", "x20")))
  o, n, word = site("t_faddw", regs); ov = reg_int(o, ("x13", "x14", "x15", "x16"))
  pred_fa = a + b >= 2 ** 256 + p
  want_fa = a + b if a + b < 2 ** 256 else a + b - p
  fa_cases.append((pred_fa, ov != want_fa % 2 ** 256 or pred_fa, (ov == want_fa) if not pred_fa else (C0 <= (a + b) % p <= 2 * C0 - 2)))
  if not pred_fa:
   fixed("t_faddw", a + b >= 2 ** 256 and (a + b) % 2 ** 64 + C0 >= 2 ** 64, n, word, ov % 2 ** 64 < C0)
 tally("fs_c2", fs_cases)
 tally("faddw_c2", fa_cases)

 # 6. BLK6 (C4a LD64B beta window sum)
 cases = []
 for kind in (0, 1, 2):
  for _ in range(20):
   w = rng.getrandbits(256); mem = E.LazyTable([])
   if kind == 2:
    es = [rng.randrange(p - p // 7, p) for _ in range(4)]
    rest = 4 * 2 ** 256 + 2 ** 200 + 2 ** 128 - 1 - sum(es)
    es += [rest // 2, rest - rest // 2]
   elif kind == 1:
    es = [p - 1 - rng.randrange(100) for _ in range(5)]
    es.append(5 * 2 ** 256 - 1 - sum(es))
    assert 0 <= es[5] < p
   if kind:
    w = sum((i + 1 + 8 * rng.randrange(1000)) << (43 * i) for i in range(6))
    for i in range(6): put(mem, E.bt_addr(w, i), es[i])
   o, n, word = site("t_blk", as_regs(w, ("x21", "x22", "x23", "x24")), mem)
   ov = reg_int(o, ("x29", "x27", "x7", "x1"))
   tot = sum(get(mem, E.bt_addr(w, i)) for i in range(6)); t = tot >> 256
   want = E.blook(w) if kind == 0 else fold(fold(tot)) % 2 ** 256
   cases.append((fold(tot) >= 2 ** 256, ov != fold(tot) % 2 ** 256, ov == want and ov % p == tot % p))
   fixed("t_blk", tot % 2 ** 128 + t * C0 >= 2 ** 128, n, word, tot % 2 ** 128 >= 2 ** 128 - t * C0)
 tally("beta_ls64_fold", cases)

 # 7. Unified tbld check (7-slot beta and 8-slot X-table)
 rprog, rlab = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "route.s")).read())
 bmem, bbase, bdesc = {}, E.BT0, 0x200000
 for idx, wi in enumerate(E.SLOTW):
  put(bmem, bdesc + 64 * idx, E.BETA * pow(2, 43 * wi, p) % p)
  put(bmem, bdesc + 64 * idx + 32, 0)
 nb = E.a64.a64_run(rprog, rlab, "tbld", dict(x0=bbase, x1=bdesc, x2=8, x3=4, x4=7, x30=E.SENTINEL), bmem)
 bok = nb == 4 + 7 * 10 + 20 * (6 * 8 + 4) and all(
  get(bmem, bbase + 256 * t + 32 * k) == E.BETA * t * pow(2, 43 * E.SLOTW[k], p) % p
  for k in range(7) for t in range(4 if k == 5 else 8))

 xmem, xbase, xdesc, cval, add0 = {}, E.XT0, 0x200000, rng.randrange(1, p), rng.randrange(1, p)
 for k in range(8):
  put(xmem, xdesc + 64 * k, cval * pow(2, 32 * k, p) % p)
  put(xmem, xdesc + 64 * k + 32, add0 if k == 0 else 0)
 nx = E.a64.a64_run(rprog, rlab, "tbld", dict(x0=xbase, x1=xdesc, x2=8, x3=8, x4=8, x30=E.SENTINEL), xmem)
 xok = nx == 4 + 8 * 10 + 20 * 64 and all(
  get(xmem, xbase + 256 * t + 32 * k) == (cval * t * pow(2, 32 * k, p) + (add0 if k == 0 else 0)) % p
  for k in range(8) for t in range(8))
 res["tbld_unified"] = dict(beta_ok=bok, xtable_ok=xok)

 res["fixup_branches"] = fix
 res["all_true"] = (bok and xok and all(d["predicate_holds"] > 0 and d["deviates_iff_predicate"] and d["model"]
         for k, d in res.items() if isinstance(d, dict) and "model" in d)
         and all(d["fixups"] > 0 and d["word_ok"] and d["count_ok"] and d["word_set_ok"] for d in fix.values()))
 return res
