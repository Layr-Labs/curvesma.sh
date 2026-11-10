M64 = (1 << 64) - 1

def expand(text):
 macros, out, cur = {}, [], None
 for raw in text.split("\n"):
  s = raw.split("//")[0].strip()
  if s.startswith(".macro"):
   parts = s[len(".macro"):].replace(",", " ").split()
   cur = macros[parts[0]] = (parts[1:], [])
   continue
  if s.startswith(".endm"):
   cur = None
   continue
  (cur[1] if cur is not None else out).append(raw)

 def emit(lines, depth):
  assert depth < 16
  res = []
  for raw in lines:
   s = raw.split("//")[0].strip()
   name = s.split()[0] if s else ""
   if raw[:1].isspace() and name in macros:
    params, body = macros[name]
    args = [a.strip() for a in s[len(name):].split(",")] if s[len(name):].strip() else []
    sub = []
    for b in body:
     for p_, a in zip(params, args):
      b = b.replace("\\" + p_, a)
     sub.append(b)
    res += emit(sub, depth + 1)
   else:
    res.append(raw)
  return res
 return "\n".join(emit(out, 0))

def _a64_ops(s):
 ops, cur, depth = [], "", 0
 for ch in s:
  if ch == "[":
   depth += 1
  if ch == "]":
   depth -= 1
  if ch == "," and depth == 0:
   ops.append(cur.strip())
   cur = ""
  else:
   cur += ch
 if cur.strip():
  ops.append(cur.strip())
 return ops

def a64_load(text):
 prog, labels = [], {}
 for raw in expand(text).split("\n"):
  line = raw.split("//")[0].rstrip()
  if not line.strip():
   continue
  if not line[0].isspace():
   labels.setdefault(line.rstrip(":"), []).append(len(prog))
   continue
  ins = line.strip()
  if ins.startswith("."):
   continue
  parts = ins.split(None, 1)
  mn, rest = parts[0], parts[1] if len(parts) > 1 else ""
  prog.append((mn, _a64_ops(rest)))
 return prog, labels

def a64_run(prog, labels, entry, regs, mem, limit=10**7, out_regs=None, visits=None, dev=None):
 M = M64
 R = {("x%d" % i): 0 for i in range(31)}
 R["sp"] = 0
 R.update(regs)
 lr0 = R["x30"]
 flags = {"C": 0, "Z": 0}
 pc = labels[entry][0]
 steps = 0
 depth = 0

 def rd(o):
  if o in ("xzr", "wzr"):
   return 0
  if o[0] == "w":
   return R["x" + o[1:]] & 0xFFFFFFFF
  return R[o]

 def wr(o, v):
  if o in ("xzr", "wzr"):
   return
  if o[0] == "w":
   R["x" + o[1:]] = v & 0xFFFFFFFF
  else:
   R[o] = v & M

 def imm(o):
  return int(o.lstrip("#"), 0)

 def ea(op, post, regoff=False):
  inner = op.strip()
  pre = inner.endswith("!")
  inner = inner.rstrip("!")[1:-1]
  parts = [p.strip() for p in inner.split(",")]
  base = parts[0]
  if len(parts) > 1 and parts[1][:1] == "x":
   assert regoff and post is None and not pre, "register offset only for plain LDR/STR"
   sh = parts[2].split() if len(parts) > 2 else ["lsl", "#0"]
   assert len(parts) <= 3 and sh[0] == "lsl" and imm(sh[1]) in (0, 3)
   return (rd(base) + (rd(parts[1]) << imm(sh[1]))) & M, None, base
  off = imm(parts[1]) if len(parts) > 1 else 0
  if post is not None:
   return rd(base), (rd(base) + imm(post)) & M, base
  if pre:
   return (rd(base) + off) & M, (rd(base) + off) & M, base
  return (rd(base) + off) & M, None, base

 def shifted(v, sh):
  kind, k = sh.split()
  assert kind in ("lsl", "lsr") and 0 <= imm(k) < 64
  return (v << imm(k)) & M if kind == "lsl" else v >> imm(k)

 def cond(c):
  return {"eq": flags["Z"] == 1, "ne": flags["Z"] == 0, "hs": flags["C"] == 1, "cs": flags["C"] == 1,
    "lo": flags["C"] == 0, "cc": flags["C"] == 0}[c]

 def store_dev(addr, v):
  if dev[0] <= addr < dev[0] + dev[1]:
   dev[2](addr - dev[0], v)

 def tgt(t):
  if t[:-1].isdigit():
   c = labels[t[:-1]]
   return [p for p in c if p < pc][-1] if t[-1] == "b" else [p for p in c if p >= pc][0]
  return labels[t][0]

 while True:
  mn, o = prog[pc]
  steps += 1
  assert steps < limit
  if visits is not None and pc in visits:
   visits[pc] += 1
  pc += 1
  if mn in ("mul", "umulh"):
   a, b = rd(o[1]), rd(o[2])
   wr(o[0], (a * b) & M if mn == "mul" else (a * b) >> 64)
  elif mn in ("adds", "adc", "adcs", "add") and not (mn == "add" and (o[2].startswith("#") or len(o) == 4)):
   a, b = rd(o[1]), rd(o[2])
   if len(o) == 4:
    sh = o[3].split()
    assert mn == "adds" and sh[0] == "lsl"
    b = (b << imm(sh[1])) & M
   cin = flags["C"] if mn in ("adc", "adcs") else 0
   t = a + b + cin
   wr(o[0], t)
   if mn in ("adds", "adcs"):
    flags["C"] = t >> 64
    flags["Z"] = 1 if (t & M) == 0 else 0
  elif mn == "add":
   if o[2].startswith("#"):
    wr(o[0], rd(o[1]) + imm(o[2]))
   else:
    wr(o[0], rd(o[1]) + shifted(rd(o[2]), o[3]))
  elif mn in ("sub", "subs", "sbc", "sbcs"):
   a = rd(o[1])
   b = imm(o[2]) if o[2].startswith("#") else rd(o[2])
   cin = flags["C"] if mn in ("sbc", "sbcs") else 1
   t = a + ((~b) & M) + cin
   wr(o[0], t)
   if mn in ("subs", "sbcs"):
    flags["C"] = t >> 64
    flags["Z"] = 1 if (t & M) == 0 else 0
  elif mn == "csel":
   wr(o[0], rd(o[1]) if cond(o[3]) else rd(o[2]))
  elif mn == "madd":
   wr(o[0], rd(o[1]) * rd(o[2]) + rd(o[3]))
  elif mn in ("cmp", "cmn", "ccmp", "ccmn"):
   if mn[1] == "c" and not cond(o[3]):
    flags["C"], flags["Z"] = imm(o[2]) >> 1 & 1, imm(o[2]) >> 2 & 1
   else:
    a, b = rd(o[0]), imm(o[1]) if o[1].startswith("#") else rd(o[1])
    if len(o) == 3:
     assert mn == "cmp" and o[2][:4] == "asr "
     b = ((b - (b >> 63 << 64)) >> imm(o[2][4:])) & M
    t = a + b if mn[-1] == "n" else a + ((~b) & M) + 1
    flags["C"], flags["Z"] = t >> 64, 1 if (t & M) == 0 else 0
  elif mn == "tst":
   t = rd(o[0]) & (imm(o[1]) if o[1].startswith("#") else rd(o[1]))
   flags["C"], flags["Z"] = 0, 1 if t == 0 else 0
  elif mn in ("orr", "and", "eor", "bic"):
   assert mn != "bic" or not o[2].startswith("#")
   b = imm(o[2]) if o[2].startswith("#") else rd(o[2])
   if len(o) == 4:
    assert not o[2].startswith("#")
    b = shifted(b, o[3])
   assert len(o) <= 4
   if mn == "bic":
    b = ~b & M
   wr(o[0], rd(o[1]) | b if mn == "orr" else rd(o[1]) ^ b if mn == "eor" else rd(o[1]) & b)
  elif mn == "lsl":
   wr(o[0], rd(o[1]) << imm(o[2]))
  elif mn == "cinc":
   wr(o[0], rd(o[1]) + (1 if cond(o[2]) else 0))
  elif mn == "cset":
   wr(o[0], 1 if cond(o[1]) else 0)
  elif mn == "bfi":
   mask = ((1 << imm(o[3])) - 1) << imm(o[2])
   wr(o[0], (rd(o[0]) & ~mask) | ((rd(o[1]) << imm(o[2])) & mask))
  elif mn == "extr":
   wr(o[0], ((rd(o[1]) << 64 | rd(o[2])) >> imm(o[3])))
  elif mn == "lsr":
   wr(o[0], rd(o[1]) >> imm(o[2]))
  elif mn == "ubfx":
   assert imm(o[2]) + imm(o[3]) <= 64
   wr(o[0], (rd(o[1]) >> imm(o[2])) & ((1 << imm(o[3])) - 1))
  elif mn == "mvn":
   assert len(o) == 2
   wr(o[0], ~rd(o[1]) & M)
  elif mn == "mov":
   wr(o[0], imm(o[1]) if o[1].startswith("#") else rd(o[1]))
  elif mn == "movk":
   sh = imm(o[2].split()[1]) if len(o) > 2 else 0
   wr(o[0], (rd(o[0]) & ~(0xFFFF << sh)) | ((imm(o[1]) & 0xFFFF) << sh))
  elif mn in ("ldp", "stp"):
   post = o[3] if len(o) > 3 else None
   addr, newb, base = ea(o[2], post)
   if mn == "ldp":
    wr(o[0], mem[addr])
    wr(o[1], mem[addr + 8])
   else:
    mem[addr] = rd(o[0])
    mem[addr + 8] = rd(o[1])
    if dev is not None:
     store_dev(addr, mem[addr])
     store_dev(addr + 8, mem[addr + 8])
   if newb is not None:
    R[base] = newb
  elif mn in ("ldr", "str"):
   post = o[2] if len(o) > 2 else None
   addr, newb, base = ea(o[1], post, regoff=True)
   if mn == "ldr":
    wr(o[0], mem[addr])
   else:
    mem[addr] = rd(o[0])
    if dev is not None:
     store_dev(addr, mem[addr])
   if newb is not None:
    R[base] = newb
  elif mn == "stadd":
   addr, newb, base = ea(o[1], None)
   mem[addr] = (mem[addr] + rd(o[0])) & M
  elif mn == "ldadd":
   addr, newb, base = ea(o[2], None)
   old = mem[addr]
   mem[addr] = (old + rd(o[0])) & M
   wr(o[1], old)
  elif mn in ("b", "b.eq", "b.ne", "b.lo", "b.hs", "b.cs", "b.cc"):
   if mn == "b" or cond(mn[2:]):
    pc = tgt(o[0])
  elif mn in ("cbz", "cbnz", "tbz", "tbnz"):
   v = rd(o[0]) if mn[0] == "c" else rd(o[0]) >> imm(o[1]) & 1
   if (v == 0) == (mn in ("cbz", "tbz")):
    pc = tgt(o[-1])
  elif mn == "bl":
   R["x30"] = pc
   depth += 1
   pc = tgt(o[0])
  elif mn == "ret":
   if depth:
    depth -= 1
    pc = R["x30"]
    continue
   assert R["x30"] == lr0, "outermost RET with a clobbered x30"
   if out_regs is not None:
    out_regs.update(R)
   return steps
  else:
   raise ValueError("unsupported instruction: " + mn)
