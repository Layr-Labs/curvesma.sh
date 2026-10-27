"""The AArch64 interpreter and instruction counter of experiment.py (proof.md section 2.6).

_asm_counts counts each function's instructions and the instructions per iteration of its single-block loops;
a64_load and a64_run execute the shipped kernels/secp_loops.s, kernels/secp_tests.s and kernels/arrive.s;
exit_check runs the loops' rare exits on constructed inputs (proof.md 2.5).  Both
first expand assembler macros (.macro NAME [params] ... .endm; a parameter p is referenced as backslash-p), as the assembler does.
"""
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
  """Per function: instructions, and for single-block loops the instructions per iteration.
    Directive lines (.text, .globl, ...) are skipped; labels are `name:` (functions), `.Lname:` or `N:`
    (local); a branch to `Nb` is a backward branch to the nearest preceding label N."""
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
        tgt = val.split()[-1]
        key = tgt[:-1] if (tgt[:-1].isdigit() and tgt.endswith("b")) else tgt
        if (mn in ("b", "cbz", "cbnz", "tbz", "tbnz") or mn.startswith("b.")) and key in pos and pos[key][0] <= idx:
          start = max(p for p in pos[key] if p <= idx)
          body = [v for k, v in items if k == "ins"][start:idx + 1]
          inner = [b for b in body[:-1] if b.split()[0] in ("b", "bl", "cbz", "cbnz", "tbz", "tbnz", "ret") or b.split()[0].startswith("b.")]
          exits = [b for b in inner if b.split()[-1].startswith(".Lr")]     # rare exits to stubs
          loops.append(dict(label=tgt, iteration=len(body), inner_branches=len(inner) - len(exits),
                   rare_exits=len(exits)))
        idx += 1
    res[name] = dict(instructions=total, loops=loops,
             loop_iteration_max=max([l["iteration"] for l in loops]) if loops else 0,
             calls=sum(1 for k, v in items if k == "ins" and v.split()[0] == "bl"))
  return res


# ---- a small AArch64 interpreter: executes the shipped kernels/secp_loops.s and kernels/secp_tests.s
# (proof.md section 2.6); flags C and Z only, which is all the kernels read
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
    # returns (address, new_base_value or None, base_register)
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

  def tgt(t):           # a named label, or N b / N f: the nearest label N before / after this instruction
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
    if mn in ("mul", "umulh", "smulh"):
      a, b = rd(o[1]), rd(o[2])
      if mn == "smulh":
        a, b = a - (a >> 63 << 64), b - (b >> 63 << 64)
      wr(o[0], (a * b) & M if mn == "mul" else (a * b) >> 64)
    elif mn == "asr":
      a = rd(o[1])
      wr(o[0], (a - (a >> 63 << 64)) >> imm(o[2]))
    elif mn in ("adds", "adc", "adcs", "add") and not (mn == "add" and (o[2].startswith("#") or len(o) == 4)):
      a, b = rd(o[1]), rd(o[2])
      if len(o) == 4:       # adds xd, xn, xm, lsl #k (shifted register; the bits shifted out are dropped)
        sh = o[3].split()
        assert mn == "adds" and sh[0] == "lsl"
        b = (b << imm(sh[1])) & M
      cin = flags["C"] if mn in ("adc", "adcs") else 0
      t = a + b + cin
      wr(o[0], t)
      if mn in ("adds", "adcs"):
        flags["C"] = t >> 64
        flags["Z"] = 1 if (t & M) == 0 else 0
    elif mn == "add":     # add xd, xn, #imm | add xd, xn, xm, lsl #k
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
    elif mn in ("cmp", "cmn", "ccmp", "ccmn"):     # flags C and Z; ccmp/ccmn set #nzcv when cond fails
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
    elif mn == "lsl":
      wr(o[0], rd(o[1]) << imm(o[2]))
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
    elif mn in ("cbz", "tbz", "tbnz"):
      bit = rd(o[0]) if mn == "cbz" else rd(o[0]) >> imm(o[1]) & 1
      if (bit == 0) == (mn != "tbnz"):
        pc = tgt(o[-1])
    elif mn == "ret":
      if out_regs is not None:
        out_regs.update(R)
      return steps
    else:
      raise ValueError("unsupported instruction: " + mn)


def exit_check(text, P):
  """Run pass 1 and pass 2 of kernels/secp_loops.s on constructed inputs that take the rare exits (proof.md 2.5),
    rotating each kind through every walk position, and compare every output word with the limb ports of the
    inline kernels (P: p_sub, p_mul, p_sqr, limbs, val, C0, q, seed); x19 must drop by 20 per exit."""
  import random
  q, C0, rr = P["q"], P["C0"], random.Random(P["seed"] + 57)
  p_sub, p_mul, p_sqr, L, _val = P["p_sub"], P["p_mul"], P["p_sqr"], P["limbs"], P["val"]
  prog, labels = a64_load(text)
  visits = {labels[k][0]: 0 for k in labels if k.startswith(".Lr")}
  WB, TB, PB, RB, SB, AB = 0x10000, 0x20000, 0x30000, 0x50000, 0x90000, (1 << 27) + (1 << 26) - 64

  def borrow_pair():                  # a - b borrows and (a - b) mod 2^64 < C0: fe_sub's second borrow chain
    b = rr.randrange(2 ** 255, q)
    a = b - (((rr.getrandbits(160) + 1) << 64) + 2 ** 64 - rr.randrange(1, C0))
    return a, b

  def carry_operand():                # a * 2^128 reduces with s4 = 1, s1 = 2^64 - 1: carry exit
    while True:
      H = (M64 << 64 | 2 ** 64 - rr.randrange(1, C0)) * pow(C0, -1, 2 ** 128) % 2 ** 128
      a = (H << 128) + 2 ** 128 - rr.randrange(1, max(2, H * C0 >> 128))
      if a < q and (a % 2 ** 128) * 2 ** 128 + H * C0 >= 2 ** 256:
        return a

  def fin(w, t, inv):                 # the walk's addition by the ports, as in _batch_port_check
    dy = p_sub(t[1], w[1])
    lam = p_mul(dy, inv)
    x3 = p_sqr(lam, [w[0], L(_val(t[0]) + 2 * C0)], True)
    return x3, p_mul(lam, p_sub(x3, w[0]), add=w[1])     # -(W + T)

  def run(ws, ts, slots, rinv):
    n, mem, regs = len(ws), {}, {}
    for i in range(n):
      for k in range(4):
        mem[WB + 64 * i + 8 * k], mem[WB + 64 * i + 32 + 8 * k] = L(ws[i][0])[k], L(ws[i][1])[k]
        mem[TB + 256 * i + 8 * k], mem[TB + 256 * i + 32 + 8 * k] = L(ts[i][0])[k], L(ts[i][1])[k]
        mem[TB + 256 * i + 192 + 8 * k] = L(ts[i][0] + 2 * C0)[k]     # the field x + 2C0
      assert ts[i][0] + 2 * C0 < 2 ** 256
      mem[PB + 8 * i] = TB + 256 * i
    for a in range(0, 160, 8):
      mem[SB - 112 + a] = 0
    W = [(L(x), L(y)) for x, y in ws]
    T = [(L(x), L(y)) for x, y in ts]
    x19 = 10 ** 9
    if slots is None:               # pass 1, then pass 2 on its slots with the true inverse
      a64_run(prog, labels, "pass1", {"x0": n, "x1": WB, "x2": PB, "x3": AB, "x19": x19, "sp": SB}, mem,
          out_regs=regs, visits=visits)
      x19 = regs["x19"]
      acc = [p_sub(T[0][0], W[0][0])]
      dxs = [None]
      for i in range(1, n):
        dxs.append(p_sub(T[i][0], W[i][0]))
        acc.append(p_mul(dxs[i], acc[i - 1]))
      bad = sum(L(sum(mem[AB + 64 * i + 32 + 8 * k] << (64 * k) for k in range(4))) != acc[i - 1]
           or [mem[AB + 64 * i + 8 * k] for k in range(4)] != dxs[i] for i in range(1, n))
      R = L(pow(_val(acc[n - 1]), -1, q))
    else:
      for i in range(1, n):
        for k in range(4):
          mem[AB + 64 * i + 8 * k], mem[AB + 64 * i + 32 + 8 * k] = L(slots[i][1])[k], L(slots[i][0])[k]
      acc, dxs, bad, R = [L(sl[0]) for sl in slots[1:]], [None] + [L(sl[1]) for sl in slots[1:]], 0, L(rinv)
    for k in range(4):
      mem[RB + 8 * k] = R[k]
    a64_run(prog, labels, "pass2", {"x0": n, "x1": WB, "x2": PB, "x3": AB, "x4": RB, "x19": x19, "sp": SB}, mem,
        out_regs=regs, visits=visits)
    out = [None] * n
    for i in range(n - 1, 0, -1):
      inv = p_mul(R, acc[i - 1])
      R = p_mul(dxs[i], R)
      out[i] = fin(W[i], T[i], inv)
    out[0] = fin(W[0], T[0], R)
    for i in range(n):
      got = ([mem[WB + 64 * i + 8 * k] for k in range(4)], [mem[WB + 64 * i + 32 + 8 * k] for k in range(4)])
      bad += got != (list(out[i][0]), list(out[i][1]))
    drop = 10 ** 9 - regs["x19"]
    return bad + (drop % 20 != 0), drop // 20

  res = dict(runs=0, walks=0, exits=0, mismatches=0)
  B = 2 ** 128
  plan = [("p1", t, n) for n in (129, 12) for t in range(3)] + [("carry", t, 20) for t in range(3)] \
    + [("lam1", t, 20) for t in range(5)] + [("lamcarry", t, 20) for t in range(2)]
  for kind, t, n in plan:
    ws, ts, slots, rinv = [], [], None, None
    if kind == "p1":                # pass 1, by (i + t) mod 3: second borrow, acc <- 2^128, carry
      xw = rr.randrange(q)
      if t % 3 == 1:              # walk 1 is a carry walk: acc_0 = dx_0 = 2^128
        xt, acc = (xw + B) % q, B
      else:                       # walk 0 takes the second borrow
        xt, xw = borrow_pair()
        acc = (xt - xw) % q
      ws.append((xw, rr.randrange(1, q)))
      ts.append((xt, rr.randrange(1, q)))
      for i in range(1, n):
        xw, c = rr.randrange(q), (i + t) % 3
        if c == 0:
          xt, xw = borrow_pair()
          dx = (xt - xw) % q
        else:
          dx = B * pow(acc, -1, q) % q if c == 1 else carry_operand()
          xt = (xw + dx) % q
        acc = acc * dx % q
        ws.append((xw, rr.randrange(1, q)))
        ts.append((xt, rr.randrange(1, q)))
    elif kind == "carry":           # pass 2, by (t - i) mod 3: inv and R carry, R <- 2^128, inv carry
      rinv, slots, R = B, [None] * n, B
      for i in range(n - 1, 0, -1):
        c = (t - i) % 3
        acc = carry_operand() if c != 1 else rr.randrange(1, q)
        dx = carry_operand() if c == 0 else B * pow(R, -1, q) % q if c == 1 else 1
        slots[i], R = (acc, dx), R * dx % q
      for i in range(n):
        ws.append((rr.randrange(q), rr.randrange(1, q)))
        ts.append((rr.randrange(q), rr.randrange(q)))
    else:                           # lambda = dy (R = 1) or dy * 2^128 (R = 2^128), every slot (1, 1)
      rinv, slots = (1 if kind == "lam1" else B), [None] + [(1, 1)] * (n - 1)
      for i in range(n):
        c, xw, yw = (i + t) % 5, rr.randrange(q), rr.randrange(1, q)
        if kind == "lamcarry":  # lambda = 2^128 * a: lambda's carry
          yt, xt = (yw + carry_operand()) % q, rr.randrange(q)
        elif c == 0:            # dy = y_T - y_W with a second borrow
          (yt, yw), xt = borrow_pair(), rr.randrange(q)
        elif c < 3:             # lambda = 1, x3 - x_W with a second borrow (c = 2: y_W near p)
          if c == 2:
            yw = q - rr.randrange(1, 2 ** 32)
          (x3, xw), yt = borrow_pair(), (yw + 1) % q
          xt = (1 - xw - x3) % q
        elif c == 3:            # lambda = 2^128, y_W = 0, x3 - x_W = a: y3's carry
          yw, yt = 0, B
          xt = (C0 - 2 * xw - carry_operand()) % q
        else:                   # lambda = 2^128: s' = 2^256 + m with m = p - x_W - x_T, m + C0 carries at limb 1
          m = (rr.getrandbits(120) << 128) + 2 ** 128 - rr.randrange(1, C0)
          xw, yt = rr.randrange(q - m - C0), (yw + B) % q
          xt = q - m - xw
        ws.append((xw, yw))
        ts.append((xt, yt))
    bad, ex = run(ws, ts, slots, rinv)
    res["runs"] += 1
    res["walks"] += n
    res["exits"] += ex
    res["mismatches"] += bad
  res["stubs"], res["stubs_entered"] = len(visits), sum(1 for v in visits.values() if v)
  res["stub_entries"] = sum(visits.values())
  return res
