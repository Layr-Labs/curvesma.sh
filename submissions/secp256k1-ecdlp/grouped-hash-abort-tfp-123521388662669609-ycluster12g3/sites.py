import os
SP = (34 << 40) + 0x7000000
SITE_KERNELS = """
\t.macro ENTER
\tmov x4,#977
\tmovk x4,#1,lsl #32
\tstr x30,[sp,#504]
\tstr xzr,[sp]
\t.endm
\t.macro LEAVE
\tldr x30,[sp,#504]
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
\tadd x30,sp,#256
\tldp x5,x6,[sp,#384]
\tldp x7,x26,[sp,#400]
\tldp x27,x28,[sp,#416]
\tldr x29,[sp,#432]
\tEPAIR_ONE NOSTEP
\tLEAVE
t_lkf:
\tENTER
\tldp x5,x6,[sp,#96]
\tldp x7,x26,[sp,#112]
\tldp x27,x28,[sp,#352]
\tldr x29,[sp,#368]
\tFGPAIR_ONE 208,NOSTEP
\tLEAVE
t_p1:
\tENTER
\tP1_7WIN
\tstr x30,[sp,#496]
\tLEAVE
t_sqinv:
	ENTER
	SQ13 x8,x9,x10,x11,x13,x14,x15,x16,x5,x6,x7,x12,x17
	F5 x13,x14,x15,x16,x12,x13,x14,x17,x12
	LEAVE
t_fs:
	ENTER
	FS x9,x10,x11,x12,x5,x6,x7,x8,x19,x20,x21,x22,x13
	LEAVE
t_faddw:
	ENTER
	FADDW x13,x14,x15,x16,x13,x14,x15,x16,x17,x18,x19,x20,x5
	LEAVE
t_blk:
	ENTER
	mov x30,x11
	PAIRB6
	LEAVE
"""
COMMON = dict(t_shx=4, t_shxt=3, t_fdx=6, t_lke=4 + 96, t_lkf=4 + 84,
       t_p1=253 + 1, t_sqinv=58 + 5, t_fs=7, t_faddw=7, t_blk=1 + 3 + 45 + 45)
EXTRA = dict(t_shx=3, t_shxt=3, t_fdx=3, t_lke=5, t_lkf=5, t_sqinv=3, t_fs=4, t_faddw=4, t_blk=5)
def site_checks(E, rng, fprog, flabels):
 p, C0, M64, put, get = E.p, E.C0, E.M64, E.put, E.get
 reg_int, as_regs = E.reg_int, E.as_regs
 prog, labels = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read() + SITE_KERNELS)
 res, fix = {}, {}
 LS = ([(E.XT0, 1 << 51), (E.BT0, 1 << 64), (SP + 128, SP + 320)], [(0x200000, 0x200100)])
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
   cases.append((low + hi >= 2 ** 192, ov != want >> 64, ov == want % 2 ** 256 >> 48 >> 16))
   fixed("t_shxt", (B % 2 ** 64) + hi >= 2 ** 64, n, word, B % 2 ** 64 > 2 ** 64 - C0)
 tally("row_shift", cases)
 cases = []
 for entry, offs_pos, tbo, ci in (("t_lke", (384, 392, 400, 408, 416, 424, 432), 264, 0),
                                   ("t_lkf", (96, 104, 112, 120, 352, 360, 368), 208, 1)):
  for kind in (0, 0, 1, 2):
   for _ in range(12):
    a1, a2 = rng.randrange(1, p), rng.randrange(1, p)
    tb, X2 = E.XT0 + (ci * 49 * (1 << 40)), 0x200040
    mem = E.LazyTable([], [(a1, rng.randrange(1, p)), (a2, rng.randrange(1, p))] * 6)
    mem[X2 + tbo] = tb
    w = rng.getrandbits(256)
    offs = E.offs7(w)
    for k, pos in enumerate(offs_pos): mem[SP + pos] = offs[k]
    b4 = rng.getrandbits(256) if entry == "t_lke" else 0
    if kind:
     b4 = (p - 1 - rng.randrange(100) if kind == 1 else rng.randrange(2 ** 200)) if entry == "t_lke" else 0
     T = ((7 if entry == "t_lke" else 6) * 2 ** 256 - 1 if kind == 1
          else 5 * 2 ** 256 + (rng.getrandbits(100) << 128) + 2 ** 128 - 1 - rng.randrange(C0))
     es = [p - 1 - rng.randrange(100) for _ in range(5 if kind == 1 else 4)]
     rest = T - b4 - sum(es)
     es += [rest // 2, rest - rest // 2] if kind == 1 else [rest // 3, rest // 3, rest - 2 * (rest // 3)]
     assert len(es) == 7 and all(0 <= e < p for e in es)
     es2 = [rng.randrange(p // 8) for _ in range(7)]
     for k in range(7):
      tk = (w >> E.W7_SHIFTS[k]) & E.W7_MASKS[k]
      put(mem, tb + (k << 43) + (tk << 6), es[k])
      put(mem, tb + (k << 43) + (tk << 6) + 32, es2[k])
    put(mem, SP + 256, b4); put(mem, SP + 288, 0)
    o, n, word = site(entry, dict(x2=X2), mem)
    ov1, ov2 = get(mem, X2), get(mem, X2 + 32)
    tot1 = sum(get(mem, tb + (k << 43) + (((w >> E.W7_SHIFTS[k]) & E.W7_MASKS[k]) << 6)) for k in range(7)) + b4
    tot2 = sum(get(mem, tb + (k << 43) + (((w >> E.W7_SHIFTS[k]) & E.W7_MASKS[k]) << 6) + 32) for k in range(7)) + b4
    t1, t2 = tot1 >> 256, tot2 >> 256
    f1 = tot1 % 2 ** 128 + t1 * C0 >= 2 ** 128
    f2 = tot2 % 2 ** 128 + t2 * C0 >= 2 ** 128
    want1 = (E.x_consts(a1)[ci] * w + (b4 if ci == 0 else E.x_consts(a1)[3])) % p
    want2 = (E.x_consts(a2)[ci] * w + (b4 if ci == 0 else E.x_consts(a2)[3])) % p
    cases.append((fold(tot1) >= 2 ** 256, ov1 != fold(tot1) % 2 ** 256,
        ov1 == E.lkfold(tot1) and ov2 == E.lkfold(tot2) and ov1 % p == tot1 % p and ov2 % p == tot2 % p
        and (kind > 0 or (ov1 % p == want1 and ov2 % p == want2))))
    if not f2:
     fixed(entry, f1, n, word, tot1 % 2 ** 128 >= 2 ** 128 - t1 * C0)
 tally("lookup_sum_fold", cases)
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
 cases = []
 for kind in (0, 1, 2):
  for _ in range(24):
   while True:
    T = (rng.randrange(2 ** 64), rng.randrange(2 ** 65, 2 ** 80), rng.randrange(p))[kind]
    a = E.sqrt_mod_p(T)
    if a is not None: break
   a += rng.choice((0, p)) if a < 2 ** 256 - p else 0
   o, n, word = site("t_sqinv", as_regs(a, ("x8", "x9", "x10", "x11")))
   ov = reg_int(o, ("x13", "x14", "x15", "x16"))
   kept = reg_int(o, ("x8", "x9", "x10", "x11")) == a
   V = fold(a * a); t = fold(V)
   cases.append((t >= p, ov != T, ov == t % 2 ** 256 and kept))
   fixed("t_sqinv", low_carry(V, 2), n, word, t % 2 ** 128 < 2 ** 65)
 tally("square_fold_inv2", cases)
 p1_ok = True
 p1_clean_cnt = 0
 for _ in range(32):
  X = E.PROD.mul(rng.randrange(1, E.r), E.GP)
  Y = E.PROD.mul(rng.randrange(1, E.r), E.GP)
  if X[0] == Y[0]: continue
  a, yx = X; b, yy = Y
  inv = pow((yy * yy - yx * yx) % p, p - 2, p)
  mem = E.LazyTable([], [X] * 12); X3 = 0x200000
  mem[X3 + 256] = E.AXYX0
  for k in range(4): put(mem, SP + 128 + 32 * k, yy * pow(2, 64 * k, p) % p)
  o, n, word = site("t_p1", dict(as_regs(inv, ("x8", "x9", "x10", "x11")), x3=X3), mem)
  u2_got = o["x28"] | (o["x29"] << 64) | (mem[SP + 496] << 128) | (o["x1"] << 192)
  v2_got = reg_int(o, ("x0", "x2", "x5", "x6"))
  inv2 = fold(fold(inv * inv)) % 2 ** 256
  Mval = E.lkfold(inv + sum((2 * yx * yx % p) * ((inv2 >> sh) & m) * pow(2, sh, p) % p for sh, m in zip(E.W7_SHIFTS, E.W7_MASKS)))
  Pval = E.lookup7(2 * yx % p, inv2)
  Nval = fold(sum(((Pval >> 64 * k) & M64) * (yy * pow(2, 64 * k, p) % p) for k in range(4))) % 2 ** 256
  p1_ok &= (u2_got == E.fs_op(Mval, Nval) and v2_got == E.faddw_op(Mval, Nval)
            and u2_got % p == pow(yy + yx, -2, p) and v2_got % p == pow(yy - yx, -2, p))
  if word == 0: p1_ok &= (n == 6 + COMMON["t_p1"]); p1_clean_cnt += 1
 res["p1_direct_u2v2"] = dict(ok=p1_ok and p1_clean_cnt > 0, clean_runs=p1_clean_cnt)
 fs_cases, fa_cases = [], []
 pairs = [(rng.randrange(p), rng.randrange(p)) for _ in range(48)]
 for _ in range(24):
  a = rng.randrange(p >> 1)
  pairs.append((a, a + (rng.randrange(1, 2 ** 190) << 64) - rng.randrange(C0)))
  b = rng.randrange(p >> 1, p)
  pairs.append((2 ** 256 - b + (rng.randrange(1, 2 ** 120) << 64) + 2 ** 64 - 1 - rng.randrange(C0), b))
  N_bad = p + 1 + rng.randrange(C0 - 1)
  M_bad = rng.randrange(N_bad - p)
  pairs.append((M_bad, N_bad))
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
 tally("fs_p1", fs_cases)
 tally("faddw_p1", fa_cases)
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
   w2 = rng.randrange(1 << 40)
   regs_blk = dict(as_regs(w, ("x28", "x29", "x11", "x1")))
   regs_blk.update(as_regs(w2, ("x0", "x2", "x5", "x6")))
   o, n, word = site("t_blk", regs_blk, mem)
   ov = reg_int(o, ("x24", "x25", "x26", "x27"))
   ov2 = reg_int(o, ("x20", "x21", "x22", "x23"))
   tot = sum(get(mem, E.bt_addr(w, i)) for i in range(6)); t = tot >> 256
   want = E.blook(w) if kind == 0 else fold(fold(tot)) % 2 ** 256
   cases.append((fold(tot) >= 2 ** 256, ov != fold(tot) % 2 ** 256, ov == want and ov % p == tot % p and (kind > 0 or ov2 == E.blook(w2))))
   fixed("t_blk", tot % 2 ** 128 + t * C0 >= 2 ** 128, n, word, tot % 2 ** 128 >= 2 ** 128 - t * C0)
 tally("beta_ls64_fold", cases)
 rprog, rlab = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "route.s")).read())
 bmem, bbase, bdesc = {}, E.BT0, 0x200000
 for idx, wi in enumerate(E.SLOTW):
  put(bmem, bdesc + 64 * idx, E.BETA * pow(2, 43 * wi, p) % p); put(bmem, bdesc + 64 * idx + 32, 0)
 nb = E.a64.a64_run(rprog, rlab, "tbld", dict(x0=bbase, x1=bdesc, x2=8, x3=2, x4=7, x5=256, x6=1 << 5, x30=E.SENTINEL), bmem)
 bok = (nb == 4 + 7 * 11 + 20 * 50
  and all(get(bmem, bbase + 256 * t + 32 * k) == E.BETA * t * pow(2, 43 * E.SLOTW[k], p) % p for k in range(7) for t in range(2 if k == 5 else 8))
  and all((bbase + 256 * t + 32 * 5) not in bmem and (bbase + 256 * t + 32 * 6) in bmem for t in range(2, 8)))
 xmem, xbase, xdesc, cval, add0 = {}, E.DEL0, 0x200000, rng.randrange(1, p), rng.randrange(1, p)
 for k, sh in enumerate(E.W7_SHIFTS):
  put(xmem, xdesc + 64 * k, cval * pow(2, sh, p) % p); put(xmem, xdesc + 64 * k + 32, add0 if k == 0 else 0)
 nx = E.a64.a64_run(rprog, rlab, "tbld", dict(x0=xbase, x1=xdesc, x2=8, x3=2, x4=7, x5=256, x6=1 << 6, x30=E.SENTINEL), xmem)
 xok = (nx == 4 + 7 * 11 + 20 * 50
  and all(get(xmem, xbase + 256 * t + 32 * k) == (cval * t * pow(2, E.W7_SHIFTS[k], p) + (add0 if k == 0 else 0)) % p for k in range(7) for t in range(2 if k == 6 else 8))
  and all((xbase + 256 * t + 32 * 5) in xmem and (xbase + 256 * t + 32 * 6) not in xmem for t in range(2, 8)))
 pmem, pbase, pdesc, c_ax, c_yx = {}, E.AXYX0, 0x200000, rng.randrange(1, p), rng.randrange(1, p)
 for k, sh in enumerate(E.W7_SHIFTS):
  put(pmem, pdesc + 128 * k, c_ax * pow(2, sh, p) % p); put(pmem, pdesc + 128 * k + 32, 0)
  put(pmem, pdesc + 128 * k + 64, c_yx * pow(2, sh, p) % p); put(pmem, pdesc + 128 * k + 96, 0)
 np_ = E.a64.a64_run(rprog, rlab, "tbld", dict(x0=pbase, x1=pdesc, x2=8, x3=2, x4=14, x5=512, x6=3 << 12, x30=E.SENTINEL), pmem)
 pok = (np_ == 4 + 14 * 11 + 20 * 100
  and all(get(pmem, pbase + 512 * t + 64 * k) == c_ax * t * pow(2, E.W7_SHIFTS[k], p) % p
   and get(pmem, pbase + 512 * t + 64 * k + 32) == c_yx * t * pow(2, E.W7_SHIFTS[k], p) % p for k in range(7) for t in range(2 if k == 6 else 8))
  and all((pbase + 512 * t + 64 * 5 + o) in pmem and (pbase + 512 * t + 64 * 6 + o) not in pmem for t in range(2, 8) for o in (0, 32)))
 cmem, cbase, cdesc, c1v, a1v, c2v, a2v = {}, E.XT0, 0x200000, rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p), rng.randrange(1, p)
 nc = 0
 for k, sh in enumerate(E.W7_SHIFTS):
  dk = cdesc + 128 * k
  put(cmem, dk, c1v * pow(2, sh, p) % p); put(cmem, dk + 32, a1v if k == 0 else 0)
  put(cmem, dk + 64, c2v * pow(2, sh, p) % p); put(cmem, dk + 96, a2v if k == 0 else 0)
  nc += E.a64.a64_run(rprog, rlab, "tbld", dict(x0=cbase + (k << 43), x1=dk, x2=8, x3=2, x4=2, x5=64, x6=3 if k == 6 else 0, x30=E.SENTINEL), cmem)
 cok = (nc == 7 * (4 + 2 * 11) + 20 * 100
  and all(get(cmem, cbase + (k << 43) + 64 * t) == (c1v * t * pow(2, E.W7_SHIFTS[k], p) + (a1v if k == 0 else 0)) % p
   and get(cmem, cbase + (k << 43) + 64 * t + 32) == (c2v * t * pow(2, E.W7_SHIFTS[k], p) + (a2v if k == 0 else 0)) % p
   for k in range(7) for t in range(2 if k == 6 else 8))
  and all((cbase + (6 << 43) + 64 * t + o) not in cmem for t in range(2, 8) for o in (0, 32)))
 res["tbld_unified"] = dict(beta_ok=bok, xtable7_ok=xok, paired7_ok=pok, contig_pair_ok=cok)
 res["fixup_branches"] = fix
 res["all_true"] = (bok and xok and pok and cok and res["p1_direct_u2v2"]["ok"]
         and all(d["predicate_holds"] > 0 and d["deviates_iff_predicate"] and d["model"]
          for k, d in res.items() if isinstance(d, dict) and "model" in d)
         and all(d["fixups"] > 0 and d["word_ok"] and d["count_ok"] and d["word_set_ok"] for d in fix.values()))
 return res
