"""AArch64 interpreter and instruction counter of experiment.py (proof.md 2.6); macros are expanded first."""
M64 = (1 << 64) - 1

def expand(text):
 """Expand .macro definitions and their invocations (nested), returning plain assembly text."""
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

def _asm_counts(text):
 """Per function: instructions, and per loop the instructions of one iteration (iteration) and of its main
    path (main: every correction skip taken).  Directive lines (.text, .globl, ...) are skipped; labels are
    `name:` (functions), `.Lname:` or `N:` (local); a branch to `Nb` is a backward branch to the nearest preceding label N."""
 funcs = {}
 cur = None
 for line in expand(text).split("\n"):
  if not line.strip():
   continue
  if not line[0].isspace():
   name = line.split(":")[0]
   if name.startswith(".L") or name.isdigit():
    funcs[cur]["items"].append(("label", name))
   else:
    cur = name
    funcs[cur] = {"items": []}
   continue
  if line.strip().startswith("."):
   continue
  funcs[cur]["items"].append(("ins", line.strip()))
 res = {}
 for name, f in funcs.items():
  items = f["items"]
  pos = {}
  idx = 0
  for kind, val in items:
   if kind == "label":
    pos.setdefault(val, []).append(idx)
   else:
    idx += 1
  total = idx
  loops = []
  idx = 0
  for kind, val in items:
   if kind == "ins":
    mn = val.split()[0]
    tgt = val.split()[-1].split(",")[-1]
    key = tgt[:-1] if (tgt[:-1].isdigit() and tgt.endswith("b")) else tgt
    if (mn in ("b", "cbz", "cbnz", "tbz", "tbnz") or mn.startswith("b.")) and key in pos and pos[key][0] <= idx:
     start = max(p for p in pos[key] if p <= idx)
     body = [v for k, v in items if k == "ins"][start:idx + 1]
     inner = [(j, b) for j, b in enumerate(body[:-1]) if b.split()[0] in ("b", "bl", "cbz", "cbnz", "tbz", "tbnz", "ret") or b.split()[0].startswith("b.")]
     exits = [b for j, b in inner if b.split()[-1].startswith(".Lr")]
     # correction skips (proof.md 2.4; secp_loops.s): forward branches to a .Lb label in the body, taken on the main path
     skip = [pos[t_][0] - start - j - 1 for j, b in inner for t_ in [b.split()[-1]] if t_.startswith(".Lb") and start + j < pos[t_][0] <= idx]
     loops.append(dict(iteration=len(body), main=len(body) - sum(skip), skips=len(skip),
         inner_branches=len(inner) - len(exits) - len(skip)))
    idx += 1
  res[name] = dict(instructions=total, loops=loops,
      calls=sum(1 for k, v in items if k == "ins" and v.split()[0] == "bl"))
 return res

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
  mn, _, rest = ins.partition(" ")
  prog.append((mn, _a64_ops(rest)))
 return prog, labels

def a64_run(prog, labels, entry, regs, mem, limit=10**7, out_regs=None, visits=None):
 """Run from label `entry` until RET; returns the number of executed instructions. regs: dict x0..x30, sp.
    visits: optional {instruction index: count}, incremented each time one of those instructions runs."""
 M = M64
 R = {("x%d" % i): 0 for i in range(31)}
 R["sp"] = 0
 R.update(regs)
 flags = {"C": 0, "Z": 0}
 pc = labels[entry][0]
 steps = 0

 def rd(o):
  if o == "xzr":
   return 0
  return R[o]

 def wr(o, v):
  if o != "xzr":
   R[o] = v & M

 def imm(o):
  return int(o.lstrip("#"), 0)

 def ea(op, post):
  inner = op.strip()
  pre = inner.endswith("!")
  inner = inner.rstrip("!")[1:-1]
  parts = [p.strip() for p in inner.split(",")]
  base = parts[0]
  off = imm(parts[1]) if len(parts) > 1 else 0
  if post is not None:
   return rd(base), (rd(base) + imm(post)) & M, base
  if pre:
   return (rd(base) + off) & M, (rd(base) + off) & M, base
  return (rd(base) + off) & M, None, base

 def cond(c):
  return {"eq": flags["Z"] == 1, "ne": flags["Z"] == 0, "hs": flags["C"] == 1, "cs": flags["C"] == 1,
    "lo": flags["C"] == 0, "cc": flags["C"] == 0}[c]

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
    wr(o[0], rd(o[1]) + (rd(o[2]) << imm(o[3].split()[1])))
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
    t = a + b if mn[-1] == "n" else a + ((~b) & M) + 1
    flags["C"], flags["Z"] = t >> 64, 1 if (t & M) == 0 else 0
  elif mn in ("orr", "and"):
   wr(o[0], rd(o[1]) | imm(o[2]) if mn == "orr" else rd(o[1]) & imm(o[2]))
  elif mn == "bfi":
   mask = ((1 << imm(o[3])) - 1) << imm(o[2])
   wr(o[0], (rd(o[0]) & ~mask) | ((rd(o[1]) << imm(o[2])) & mask))
  elif mn == "extr":
   wr(o[0], ((rd(o[1]) << 64 | rd(o[2])) >> imm(o[3])))
  elif mn == "lsr":
   wr(o[0], rd(o[1]) >> imm(o[2]))
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
   if newb is not None:
    R[base] = newb
  elif mn in ("ldr", "str"):
   post = o[2] if len(o) > 2 else None
   addr, newb, base = ea(o[1], post)
   if mn == "ldr":
    wr(o[0], mem[addr])
   else:
    mem[addr] = rd(o[0])
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
   if out_regs is not None:
    out_regs.update(R)
   return steps
  else:
   raise ValueError("unsupported instruction: " + mn)

def exit_check(text, P):
 """Rare exits of kernels/secp_loops.s on constructed inputs against the limb ports P (evidence.md item 5)."""
 import random
 q, C0, rr = P["q"], P["C0"], random.Random(P["seed"] + 57)
 sub, L, V, rows, rmul = P["p_sub"], P["limbs"], P["val"], P["p_rows"], P["p_rmul"]
 prog, labels = a64_load(text)
 visits = {labels[k][0]: 0 for k in labels if k.startswith(".Lr")}
 WB, TB, PE, RB, SB, AB = 0x10000, 0x20000, 3 << 26, 0x50000, 0x90000, 1 << 26
 res = dict(runs=0, walks=0, exits=0, borrows=0, mismatches=0)

 def put(mem, a, v):
  for k in range(4):
   mem[a + 8 * k] = L(v)[k]

 def get(mem, a):
  return [mem[a + 8 * k] for k in range(4)]

 def run(fn, n, mem, want_exits=None, k=None, **regs):  # k: the borrows prod or pass 2 must return in x0
  for a in range(0, 160, 8):
   mem[SB - 112 + a] = 0
  pb = PE - 8 * n  # tp ends where bit 26 is first set; zero words beside it (proof.md 2.5)
  for i in range(-2, n + 2):
   mem[pb + 8 * i] = TB + 256 * i if 0 <= i < n else 0
  out = {}
  a64_run(prog, labels, fn, dict(regs, x0=n, x1=WB, x2=pb, x3=AB, x19=10 ** 9, sp=SB), mem, out_regs=out, visits=visits)
  drop = 10 ** 9 - out["x19"]
  res["runs"] += 1
  res["walks"] += n
  res["exits"] += drop // 20
  res["borrows"] += k or 0
  return drop % 20 != 0 or (want_exits is not None and drop != 20 * want_exits) or (k is not None and out["x0"] != k)

 def borrow_pair():  # a < b, (a - b) mod 2^64 < C0: fe_sub(a, b) exits; half the time its borrow reaches limb 3
  b, k = rr.randrange(2 ** 255, q), rr.choice([rr.getrandbits(160) + 2, rr.randrange(1, 2 ** 20) << 128])
  return b - (k << 64) + rr.randrange(1, C0), b

 def carry_operand():
  while True:
   H = (M64 << 64 | 2 ** 64 - rr.randrange(1, C0)) * pow(C0, -1, 2 ** 128) % 2 ** 128
   a = (H << 128) + 2 ** 128 - rr.randrange(1, max(2, H * C0 >> 128))
   if a < q and (a % 2 ** 128) * 2 ** 128 + H * C0 >= 2 ** 256:
    return a

 mi = lambda: (M64 - rr.randrange(2 ** 20)) * pow(C0, -1, 2 ** 64) % 2 ** 64

 def row_targets(w):
  # Rs: rows whose three shifts carry (w: the first wraps); a0*Rd in [T, T + a0): both folds carry and wrap
  b2 = mi()
  while b2 >> 53 == 0:
   b2 = mi()
  Rs = (mi() << 192) + (2 ** 192 - 1 if w else M64 + (rr.randrange(2 ** 63, 2 ** 64) << 64) + (b2 << 128))
  a0 = rr.randrange(2 ** 63, 2 ** 64)
  return Rs, -(-(2 ** 318 + 2 ** 256 - 2 ** 65) // a0), a0

 def row_walk(c, A, Rs, Rd, a0):  # (operand of the first product, dx) for the row kind c at rows of A
  if c == "F":
   return a0, a0
  return rr.randrange(q), (Rd if c == "S" else Rs) * pow(V(A), -1, q) % q

 bad = 0
 # pass 1's dx, prod's dy: borrows in every walk (t = 0) or the even walks (t = 1); 20 per predicted exit
 for fn in ("pass1", "prod"):
  for n in ((129, 12) if fn == "pass1" else (257, 12)):
   for t in range(2):
    mem, want, ex, A = {}, [], 0, L(1)
    for i in range(n):
     (at, aw) = borrow_pair() if not (t and i % 2) else (lambda b: (rr.randrange(b, q), b))(rr.randrange(q))
     ex += at < aw
     d = sub(L(at), L(aw))
     if fn == "pass1":
      put(mem, WB + 64 * i, aw), put(mem, TB + 256 * i, at)
      want.append(d)
     else:
      dx = rr.randrange(1, q)
      put(mem, WB + 64 * i + 32, aw), put(mem, TB + 256 * i + 32, at), put(mem, AB + 64 * i, dx)
      B = rows(A)
      want.append(rmul(B, d))
      A = rmul(B, L(dx))
    bad += run(fn, n, mem, ex, ex if fn == "prod" else None)
    off = 0 if fn == "pass1" else 32
    bad += sum(get(mem, AB + 64 * i + off) != want[i] for i in range(n)) + (fn == "prod" and get(mem, AB + 64 * n + 32) != A)
 # prod: acc's rows (Z sets acc = Rs, S shifts Rs and sets Rd, F's dy = dx = a0 carry both folds); y(W) = 0
 for n, t, w in [(n, t, 0) for n in (257, 20) for t in range(3)] + [(12, 0, 1)]:
  mem, A, want = {}, L(1), []
  Rs, Rd, a0 = row_targets(w)
  for i in range(n):
   dy, dx = row_walk("ZSF"[(t + i) % 3], A, Rs, Rd, a0)
   put(mem, AB + 64 * i, dx), put(mem, TB + 256 * i + 32, dy), put(mem, WB + 64 * i + 32, 0)
   B = rows(A)
   want.append(rmul(B, L(dy)))
   A = rmul(B, L(dx))
  bad += run("prod", n, mem, k=0)
  bad += sum(get(mem, AB + 64 * i + 32) != want[i] for i in range(n)) + (get(mem, AB + 64 * n + 32) != A)
 # pass 2: R's rows as in prod; then dx = 1, R = 1 (lambda = e): u's exit, y3's carry, x3's carry (evidence item 5)
 plan = [("rows", t, n, 0) for n in (257, 20) for t in range(3)] + [("rows", 0, 12, 1)] \
  + [("lam", t, n, 0) for n in (257, 20) for t in range(5)]
 for kind, t, n, w in plan:
  mem, ws, ts, slots = {}, [], [], [None] * n
  if kind == "rows":
   Rs, Rd, a0 = row_targets(w)
   R = L([rr.randrange(q), Rs, Rd][t % 3])
   rinv, A = V(R), R
   for i in range(n - 1, -1, -1):
    slots[i] = row_walk("ZSF"[(t + n - 1 - i) % 3], A, Rs, Rd, a0)
    A = rmul(rows(A), L(slots[i][1]))
   for i in range(n):
    ws.append((rr.randrange(q), rr.randrange(1, q)))
    ts.append(rr.randrange(q - 2 * C0))
  else:
   rinv = 1
   for i in range(n):
    c, xw, yw, lam = (i + t) % 5, rr.randrange(q), rr.randrange(1, q), 2 ** 128
    if c == 0:
     xt, lam = rr.randrange(q - 2 * C0), rr.randrange(q)
    elif c < 3:
     if c == 2:
      yw = q - rr.randrange(1, 2 ** 32)
     (x3, xw), lam = borrow_pair(), 1
     xt = (1 - xw - x3) % q
    elif c == 3:
     yw = 0
     xt = (C0 - 2 * xw - carry_operand()) % q
    else:
     m = (rr.getrandbits(120) << 128) + 2 ** 128 - rr.randrange(1, C0)
     xw = rr.randrange(q - m - C0)
     xt = q - m - xw
    ws.append((xw, yw))
    ts.append(xt)
    slots[i] = (lam, 1)
  for i in range(n):
   assert ts[i] + 2 * C0 < 2 ** 256
   put(mem, WB + 64 * i, ws[i][0]), put(mem, WB + 64 * i + 32, ws[i][1]), put(mem, TB + 256 * i + 192, ts[i] + 2 * C0)
   put(mem, AB + 64 * i + 32, slots[i][0]), put(mem, AB + 64 * i, slots[i][1])
  put(mem, RB, rinv)
  R, want = L(rinv), [None] * n
  for i in range(n - 1, -1, -1):
   B = rows(R)
   lam, R = rmul(B, L(slots[i][0])), rmul(B, L(slots[i][1]))
   xw, yw = L(ws[i][0]), L(ws[i][1])
   x3 = P["p_sqr"](lam, [xw, L(ts[i] + 2 * C0)], True)
   want[i] = (x3, P["p_mul"](lam, sub(x3, xw), add=yw))
  bad += run("pass2", n, mem, k=sum(V(x3) < x for (x3, _), (x, _) in zip(want, ws)), x4=RB)
  bad += sum((get(mem, WB + 64 * i), get(mem, WB + 64 * i + 32)) != want[i] for i in range(n))
 res["mismatches"] = bad
 res["stubs"], res["stubs_entered"] = len(visits), sum(1 for v in visits.values() if v)
 res["stub_entries"] = sum(visits.values())
 return res
