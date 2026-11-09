"""AArch64 interpreter for exercised kernel/audit opcodes"""
from functools import lru_cache
M64 = (1 << 64) - 1
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
@lru_cache(maxsize=16)
def _a64_parse(text):
 prog, labels = [], {}
 for raw in text.split("\n"):
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
def a64_load(text):
 p, l = _a64_parse(text)
 return [(mn, ops.copy()) for mn, ops in p], {k:v.copy() for k,v in l.items()}
def a64_run(prog, labels, entry, regs, mem, limit=10**7, out_regs=None, visits=None, dev=None):
 M = M64
 R = {("x%d" % i): 0 for i in range(31)}
 R["sp"] = 0
 R.update(regs)
 lr0 = R["x30"]
 flags = {"C": 0, "Z": 0}
 pc = labels[entry][0]
 steps = 0
 def rd(o):
  if o in R:
   return R[o]
  if o in ("xzr", "wzr"):
   return 0
  if o[0] == "w":
   return R["x" + o[1:]] & 0xFFFFFFFF
  return R[o]
 def wr(o, v):
  if o in R:
   R[o] = v & M
   return
  if o in ("xzr", "wzr"):
   return
  if o[0] == "w":
   R["x" + o[1:]] = v & 0xFFFFFFFF
  else:
   R[o] = v & M
 def imm(o):
  return int(o.lstrip("#"), 0)
 _ea_c = {}
 def ea(op, post, regoff=False):
  p = _ea_c.get(op)
  if p is None:
   inner = op.strip()
   pre = inner.endswith("!")
   inner = inner.rstrip("!")[1:-1]
   parts = [x.strip() for x in inner.split(",")]
   base = parts[0]
   if len(parts) > 1 and parts[1][:1] == "x":
    sh = parts[2].split() if len(parts) > 2 else ["lsl", "#0"]
    p = (True, base, parts[1], imm(sh[1]), pre)
   else:
    off = imm(parts[1]) if len(parts) > 1 else 0
    p = (False, base, off, pre)
   _ea_c[op] = p
  if p[0]:
   _, base, r2, sh, pre = p
   assert regoff and post is None and not pre
   return (rd(base) + (rd(r2) << sh)) & M, None, base
  _, base, off, pre = p
  if post is not None:
   return rd(base), (rd(base) + imm(post)) & M, base
  if pre:
   return (rd(base) + off) & M, (rd(base) + off) & M, base
  return (rd(base) + off) & M, None, base
 def cond(c):
  if c == "ne": return flags["Z"] == 0
  if c == "eq": return flags["Z"] == 1
  if c in ("hs", "cs"): return flags["C"] == 1
  if c in ("lo", "cc"): return flags["C"] == 0
  return False
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
   a = rd(o[1])
   b = imm(o[2]) if o[2].startswith("#") else rd(o[2])
   if len(o) == 4:
    sh = o[3].split()
    assert mn == "adds" and sh[0] in ("lsl", "lsr")
    b = ((b << imm(sh[1])) & M) if sh[0] == "lsl" else (b >> imm(sh[1]))
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
    wr(o[0], rd(o[1]) + (rd(o[2]) << imm(o[3].split()[1])))
  elif mn in ("sub", "subs", "sbc", "sbcs"):
   a = rd(o[1])
   b = imm(o[2]) if o[2].startswith("#") else rd(o[2])
   if len(o) == 4:
    sh = o[3].split()
    assert mn in ("sub", "subs") and sh[0] in ("lsl", "lsr")
    b = ((b << imm(sh[1])) & M) if sh[0] == "lsl" else (b >> imm(sh[1]))
   cin = flags["C"] if mn in ("sbc", "sbcs") else 1
   t = a + ((~b) & M) + cin
   wr(o[0], t)
   if mn in ("subs", "sbcs"):
    flags["C"] = t >> 64
    flags["Z"] = 1 if (t & M) == 0 else 0
  elif mn == "csel":
   wr(o[0], rd(o[1]) if cond(o[3]) else rd(o[2]))
  elif mn in ("cmp", "cmn", "ccmp", "ccmn"):
   if mn[1] == "c" and not cond(o[3]):
    flags["C"], flags["Z"] = imm(o[2]) >> 1 & 1, imm(o[2]) >> 2 & 1
   else:
    a, b = rd(o[0]), imm(o[1]) if o[1].startswith("#") else rd(o[1])
    if len(o) == 3:
     assert mn == "cmp" and o[2][:4] == "lsl "
     b = (b << imm(o[2][4:])) & M
    t = a + b if mn[-1] == "n" else a + ((~b) & M) + 1
    flags["C"], flags["Z"] = t >> 64, 1 if (t & M) == 0 else 0
  elif mn == "and":
   b = imm(o[2]) if o[2].startswith("#") else rd(o[2])
   wr(o[0], rd(o[1]) & b)
  elif mn == "lsr":
   wr(o[0], rd(o[1]) >> imm(o[2]))
  elif mn == "udiv":
   b = rd(o[2]); wr(o[0], rd(o[1]) // b if b else 0)
  elif mn == "extr":
   wr(o[0], ((rd(o[1]) << 64 | rd(o[2])) >> imm(o[3])))
  elif mn == "bfi":
   f = ((1 << imm(o[3])) - 1) << imm(o[2])
   wr(o[0], rd(o[0]) & ~f | rd(o[1]) << imm(o[2]) & f)
  elif mn == "mov":
   wr(o[0], imm(o[1]) if o[1].startswith("#") else rd(o[1]))
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
  elif mn in ("ldr", "str", "ldur"):
   post = o[2] if len(o) > 2 else None
   assert mn != "ldur" or (post is None and "!" not in o[1]), "ldur: plain offset only"
   addr, newb, base = ea(o[1], post, regoff=(mn != "ldur"))
   if mn in ("ldr", "ldur"):
    wr(o[0], mem[addr])
   else:
    mem[addr] = rd(o[0])
    if dev is not None:
     store_dev(addr, mem[addr])
   if newb is not None:
    R[base] = newb
  elif mn in ("b", "b.eq", "b.ne", "b.lo", "b.hs", "b.cs", "b.cc"):
   if mn == "b" or cond(mn[2:]):
    pc = tgt(o[0])
  elif mn in ("cbz", "cbnz", "tbz", "tbnz"):
   v = rd(o[0]) if mn[0] == "c" else rd(o[0]) >> imm(o[1]) & 1
   if (v == 0) == (mn in ("cbz", "tbz")):
    pc = tgt(o[-1])
  elif mn == "ret":
   assert R["x30"] == lr0, "outermost RET with a clobbered x30"
   if out_regs is not None:
    out_regs.update(R)
   return steps
  else:
   raise ValueError("unsupported instruction: " + mn)
