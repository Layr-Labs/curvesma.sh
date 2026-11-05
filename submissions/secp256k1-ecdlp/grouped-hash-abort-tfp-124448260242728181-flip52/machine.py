import os, random
from collections import defaultdict

K1, K2B = 16, 32
B1, B2 = 64 * K1, 64 * K2B
UP = (1 << 64) - (1 << 52)
LO1, H1 = 1 << 51, 1 << 41
LC2, RC2 = 1 << 44, 1 << 45
MB2 = (1 << 48) - 128
PW = (1 << 48) - 168
CNT = (1 << 48) - 8192
CNTE = (1 << 48) - 6144
DONE = (1 << 48) - 4096
ST = (1 << 48) - 160
HT = (1 << 64) - (1 << 50)
LINK = 0xC000
WINDOW = (LINK, LINK + 4096)
STALE = 0x5A5A5A5A5A5A5A5A
M41 = 2 ** 41 - 1
DM = (2 ** 64 - 1) ^ M41
H2S_REGS = {}

def blk1(d): return LO1 + B1 * d if d < H1 else UP + B1 * (d - H1)
def blk2(d): return UP + B2 * d
LC1_TABLES = [(0, 8 * H1, LO1 + 32, B1), (8 * H1, 16 * H1, UP + 32, B1), (LO1, LO1 + B1 * H1, 0, 0), (UP, UP + B1 * H1, 0, 0)]
def dest(word): return word >> 3 & M41
def dest1(word): return word >> 3

class Port1:
 def __init__(self, mem):
  self.mem, self.blocks, self.pend = mem, [], []
  self.dev = (0, 1 << 64, self.store, self.read, lambda: not self.pend)
 def store(self, a, v):
  if a % B1 == 0 and (LO1 <= a < LO1 + B1 * H1 or UP <= a < UP + B1 * H1): self.pend.append((a, v))
 def read(self):
  for a, v in self.pend:
   self.blocks.append([v] + [self.mem.get(a + 8 * w, 0) for w in (1, 2, 3)] + [self.mem.get(a + 32 + 8 * w, STALE) for w in range(4 * K1)])
  self.pend = []

class Port2:
 def __init__(self, mem, mb=MB2):
  self.mem, self.msgs, self.blocks, self.mb, self.pend = mem, [], [], mb, []
  self.dev = (0, 1 << 64, self.store, self.read, lambda: not self.pend)
 def store(self, a, v):
  mb = self.mb
  if a == mb + 40 or (not mb <= a < mb + 48 and a >= UP and (a - UP) % B2 == 0): self.pend.append((a, v))
 def read(self):
  for a, v in self.pend:
   if a == self.mb + 40: self.msgs.append([self.mem.get(self.mb + 8 * w, 0) for w in range(5)] + [v])
   else: self.blocks.append([v, self.mem.get(a + 8, 0), self.mem.get(a + 16, 0), self.mem.get(a + 24, 0)] + [self.mem.get(a + 32 + 8 * w, STALE) for w in range(4 * K2B)])
  self.pend = []

class LinkView:
 def __init__(self, mem, msg, writes): self.mem, self.msg, self.writes = mem, msg, writes
 def __getitem__(self, a): return self.msg[(a - LINK) // 8] if LINK <= a < LINK + 8 * len(self.msg) else self.mem[a]
 def get(self, a, d=None): return self[a] if LINK <= a < LINK + 8 * len(self.msg) else self.mem.get(a, d)
 def __setitem__(self, a, v):
  assert not LINK <= a < LINK + 8 * len(self.msg)
  if a not in (CNT, CNTE): self.writes.append(a)
  self.mem[a] = v

def deliver(E, prog, labels, mem, msg, calls, cert, writes=None, fused=False, r2=None, formula=None, cnt=CNT):
 e, o = "dlv" if len(msg) == 6 else "dlvh", {}
 mem.setdefault(cnt, 0); c0 = mem[cnt]
 n = E.a64.a64_run(prog, labels, e, dict(x0=LINK, x28=1, x29=cnt, x30=E.SENTINEL), LinkView(mem, msg, [] if writes is None else writes),
       out_regs=o, ls64=([WINDOW] + ([r2] if r2 else []), [r2] if r2 else []))
 ok = o["x28"] == 1 and o["x29"] == cnt and mem[cnt] == c0 + 1
 calls.append((e, n if ok else -1, cert.DLV + cert.RCV if e == "dlv" else formula))

def h2s_model(L, cert):
 if L == K1: return 6 + 2 + 6 * L
 if L == 0: return 7
 return 14 + 1 + (L % 2 == 0) + 2 * ((L - 1) // 8) + 6 * L

def dlvh_model(n1, n2, K):
 if n1 == 0 and n2 == K: return 15 + 8 + 2 * (K - 2)
 def seg(n, m):
  if m == 0: return 1
  head = n % 2; r = m - head
  return 7 + head * 5 + 4 * (r // 2) + 4 * (r % 2) + 2 * (r // 8)
 return 12 - (n1 > 0) + 12 + seg(1, n1) + seg(1 + n1, n2) + seg(1 + n1 + n2, K - n1 - n2)

def h2s_state(o): return o["x14"] << 64 | o["x13"], o["x5"] << 64 | o["x15"]

def h2s0(E, prog, labels, nd, end, gens):
 r0 = {}
 E.a64.a64_run(prog, labels, "h2s0", dict(x12=end, x0=gens >> 64, x8=gens % 2 ** 64, x30=E.SENTINEL), nd.mem, out_regs=r0)
 assert (r0["x3"], r0["x10"], r0["x18"], r0["x12"], r0["x29"], r0["x6"], r0["x28"]) == (PW, LC2, 32 * K2B, end, DONE, 1 << 62, 0)
 assert h2s_state(r0) == (1 << 126, gens)
 return r0

def deliver_arr(E, prog, labels, nd, msg, end, extra, calls, cert, port, tag=0, gens=1 << 68):
 L, o = msg[1] // 32, {}
 r0 = H2S_REGS.get((id(nd), end)) or h2s0(E, prog, labels, nd, end, gens)
 nd.mem[PW] = tag; nd.mem.setdefault(DONE, 0); d0, (D0, f0) = nd.mem[DONE], h2s_state(r0)
 regs = dict(r0, x1=STALE, x19=STALE, sp=nd.SP, x30=E.SENTINEL)
 n = E.a64.a64_run(prog, labels, "h2s", regs, LinkView(nd.mem, msg, []), out_regs=o, dev=port and port.dev, ls64=([WINDOW], []))
 H2S_REGS[(id(nd), end)] = o
 assert h2s_model(L, cert) == cert.h2s_count(L)
 fin = L < K1
 D1_ = D0 + (msg[2] if fin else 0) - (1 << 126 if fin and f0 == 1 else 0) - 1
 done = D1_ == 0
 assert nd.mem[DONE] == d0 + done and h2s_state(o) == ((1 << 126, gens) if done else (D1_, f0 - fin))
 calls.append(("h2s", n, h2s_model(L, cert) + extra + (cert.FIN_CNT + (f0 == 1) if fin else 0) + cert.RCV1 + cert.H2S_DONE * done))
 return done

class Sender1:
 def __init__(self, cert): self.cert, self.n, self.slots, self.cnt, self.msgs, self.seq, self.extra = cert, {}, {}, {}, [], [], 0
 def block(self, d, n, c): return [8 * d, 32 * n, c, 0] + [w for r in self.slots.get(d, [None] * K1) for w in (r or [STALE] * 4)]
 def put(self, d, rec):
  s, n = self.slots.setdefault(d, [None] * K1), self.n.get(d, 0)
  s[n] = rec; self.seq.append((d, rec)); self.n[d] = (n + 1) % K1
  if n + 1 == K1:
   self.cnt[d] = self.cnt.get(d, 0) + 1; self.msgs.append(self.block(d, K1, self.cnt[d])); self.extra += self.cert.FLUSH1_X
 def finals(self, ds):
  out = []
  for d in ds: out.append(self.block(d, self.n.get(d, 0), self.cnt.get(d, 0) + 1)); self.cnt[d], self.n[d] = 0, 0
  return out

def state1_ok(mem, snd):
 for d, s in snd.slots.items():
  b = blk1(d)
  if mem[8 * d] != b + 32 + 32 * snd.n.get(d, 0) or mem[b + 16] != snd.cnt.get(d, 0): return False
  if any([mem[b + 32 * t + 32 + 8 * w] for w in range(4)] != s[t] for t in range(K1) if s[t] is not None): return False
 return True

class Sender:
 def __init__(self, rc0, end, cert, tag=0):
  self.rc0, self.end, self.cert, self.tag = rc0, end, cert, tag
  self.n, self.rc, self.slots, self.msgs, self.seq, self.sent, self.drops, self.extra = {}, {}, {}, [], [], {}, 0, 0
 def cursor(self, d): return self.rc.get(d, self.rc0(d))
 def put(self, d, rec):
  s, n = self.slots.setdefault(d, [None] * K2B), self.n.get(d, 0)
  s[n] = rec; self.seq.append((d, rec)); self.n[d] = (n + 1) % K2B
  if n + 1 < K2B: return
  c = self.cursor(d); self.rc[d] = c + 32 * K2B
  self.extra += self.cert.FLUSH2_X if c < self.end else self.cert.FLUSH2_D
  if c < self.end: self.msgs.append([LC2 + 8 * d, c, 32 * K2B, self.tag] + [w for r in s for w in r]); self.sent[d] = self.sent.get(d, 0) + 1
  else: self.drops += 1; self.rc[d] = self.end
 def final(self, d):
  n, c = self.n.get(d, 0), self.cursor(d)
  msg = [LC2 + 8 * d, c, 32 * n, self.tag] + [w for r in self.slots[d] for w in (r or [STALE] * 4)] if n and c < self.end else None
  return msg, min(c + 32 * n, self.end)
 def fin_count(self, ds):
  f = 1
  for d in ds:
   n = self.n.get(d, 0)
   f += self.cert.FINH_DEST + (0 if not n else self.cert.FINH_SENT if self.cursor(d) < self.end else self.cert.FIN_HELD)
  return f

def keys_reference(E, mem, sb, q2, W2, KN, snd, cert):
 meter, out = 0, []
 for l in range(q2 - 1, -1, -1):
  e = sb + E.ENTRY * l
  for s, ptr in ((0, e + 128), (1, e - 256)):
   f6 = E.fracs(mem, e, s); x18 = (KN + (ptr << 30)) & E.M64; dev = 0
   for t in range(3):
    f = E.record_fields(E.leaf_w(f6, t))
    meter += cert.JPATH if f[0] else 0; dev |= bool(f[0])
    snd.put(f[1], [f[3], f[4], W2, x18]); out.append(f)
   meter += cert.SJOIN if dev else 0
 return meter, out

def plant(E, mem, sb, sixes):
 for l, six in enumerate(sixes):
  for s in (0, 1):
   v = [E.t128(E.BT[m] * c % E.p) for c in E.inverse_dft(six[3 * s:3 * s + 3]) for m in (0, 1)]
   for i in (4, 5): v[i] = v[i] >> 64 << 64 | (v[i] & E.M64) >> 3
   for o, w in zip(E.FRAC_OFF[s], v): E.put(mem, sb + E.ENTRY * l + o, w, 2)

def state_ok(mem, snd):
 for d, s in snd.slots.items():
  b, e = blk2(d), LC2 + 8 * d
  if mem[e] != b + 32 + 32 * snd.n.get(d, 0) or mem[e | RC2] != snd.cursor(d): return False
  if any([mem[b + 32 * t + 32 + 8 * w] for w in range(4)] != s[t] for t in range(K2B) if s[t] is not None): return False
 return True

def unpack(dm, m, wr=None):
 if len(m) == 6:
  for w in range(4): dm[m[1] + 8 * w] = m[2 + w]
  return
 for t in range(m[2] // 32):
  for w in range(4): dm[m[1] + 32 * t + 8 * w] = m[4 + 4 * t + w]
  if wr is not None: wr.append(m[1] + 32 * t)

def rule_ok(dm, writes, snd, finals):
 ok, cnt = True, {}
 for d, rec in snd.seq:
  t = cnt[d] = cnt.get(d, -1) + 1; c = snd.rc0(d) + 32 * t
  if c < snd.end and (finals is not None or t < K2B * snd.sent.get(d, 0)):
   ok &= [dm[d].get(c + 8 * w) for w in range(4)] == rec
 if finals is not None: ok &= all(f == min(snd.rc0(d) + 32 * (cnt.get(d, -1) + 1), snd.end) for d, f in finals.items())
 return ok and all(snd.rc0(d) <= a < snd.end for d, w in writes.items() for a in w)

class Node:
 HB, SB, SP = 0x100000, 1 << 35, 0x7000000
 R2, E2, E2B = 0x3000020, 0x4000000, 0x4800000
 CL = 0x6200000
 def __init__(self, mem=None): self.mem = {} if mem is None else mem

def run_k(E, prog, labels, entry, regs, nd, calls, formula, port=None):
 o = {}
 regs = dict(regs, sp=nd.SP, x30=E.SENTINEL)
 calls.append((entry, E.a64.a64_run(prog, labels, entry, regs, nd.mem, out_regs=o, dev=port and port.dev), formula))
 return o

class HomeJoin:
 def __init__(self, cert, CC, ncl):
  self.cert, self.CC, self.ncl = cert, CC, ncl
  self.chains, self.heads, self.nextd, self.recs, self.cum, self.seen = {}, {}, {}, {}, [], []
  self.nls0, self.lists, self.paths, self.segs = 0, [0, 0], [0, 0], {}
 def block(self, m, F, NB):
  cert, c, cnt = self.cert, m[1], m[2]
  lo, hi = max(F, c), c + cnt
  n1, n2 = (lo - c) // 32, (hi - lo) // 32
  fast = F <= c and cnt == 32 * K2B
  self.paths[0 if fast else 1] += 1
  f = dlvh_model(n1, n2, K2B)
  ok = f == cert.dlvh_copy(n1, n2) and (f == cert.DLVHB if fast else f <= cert.DLVHS)
  k = (n1 % 2, n2 % 2, n1 > 0, n2 > 0); self.segs[k] = self.segs.get(k, 0) + 1
  for t in range(n1, n1 + n2):
   a, rec = c + 32 * t, m[4 + 4 * t:8 + 4 * t]
   late = not ((t + 1) % 2 == 0 and t + 2 <= n1 + n2)
   self.recs[a] = rec; self.seen.append((rec[0], rec[1], a))
   b = rec[1] & (NB - 1); ch = self.chains.setdefault(b, [])
   f += cert.HCHG_NEW
   if ch:
    match = None
    for j, e in enumerate(reversed(ch)):
     f += cert.HCHG_VISIT
     if e[1] == rec[1]:
      f += cert.HCHG_VER
      if (e[0] & DM) == (rec[0] & DM): match = (e, j); break
     if j == len(ch) - 1: f += cert.HCHG_END * (j % 2); break
    if match is not None:
     e, j = match; f -= 1
     if len(self.cum) < self.ncl:
      self.cum.append((e[2], a)); f += cert.HCHG_LS1 - cert.HCHG_LSA * (not late); self.lists[late] += 1
     else: self.nls0 += 1; f += cert.HCHG_LS0
   self.nextd[a] = self.heads.get(b, 0)
   self.heads[b] = (a >> 1) + self.CC & (2 ** 64 - 1)
   ch.append((rec[0], rec[1], a))
  return f + cert.RCV, ok

def home_state(nd, ncl, CC):
 nd.mem[ST], nd.mem[ST + 8], nd.mem[ST + 16], nd.mem[ST + 24] = nd.CL, nd.CL + 16 * ncl, CC, -2 * CC & (2 ** 64 - 1)

def part_b(E, cert):
 rg = random.Random(E.SEED + 32452843)
 rt = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "route.s")).read())
 lk = E.a64.a64_load(open(os.path.join(E.HERE, "kernels", "leaf.s")).read()) + ("keys",)
 D1S, D2S, G, q2 = (0, 1, H1, H1 + 1), (0, 1, H1 - 2, H1 - 1), 4, 12
 calls = []
 gens = {g: Node(defaultdict(int)) for g in range(G)}
 rels = {d1: Node() for d1 in D1S}
 homes = {(d1, d2): Node(defaultdict(int)) for d1 in D1S for d2 in D2S}
 pool = [E.planted_x(rg, rg.choice(D1S), rg.choice(D2S))[1] for _ in range(10)]
 batches, shared = {}, 0
 for g in gens:
  lst = []
  for _ in range(6 * q2):
   if rg.random() < 0.15: lst.append(E.planted_x(rg, 0, 0, rg.choice(pool))[0]); shared += 1
   else: lst.append(E.planted_x(rg, rg.choice(D1S), rg.choice(D2S))[0])
  batches[g] = [lst[6 * l:6 * l + 6] for l in range(q2)]
 msgs_ok, lab_ok, blocks = True, True, 0
 arrivals = {d1: [] for d1 in D1S}
 for g, nd in gens.items():
  for t0, n, first in ((0, 2, LO1 + 32), (8 * H1, 2, UP + 32)):
   run_k(E, rt[0], rt[1], "init3", dict(x0=n, x1=t0, x2=first, x3=B1), nd, calls, cert.INIT3_BKT * n + 1)
  plant(E, nd.mem, nd.SB, batches[g])
  KN = (g >> 64) | (5 + g) << 4
  snd = Sender1(cert)
  meter, fl = keys_reference(E, nd.mem, nd.SB, q2, g, KN, snd, cert)
  lab_ok &= [f[1:] for f in fl] == [E.class_label(x) for six in batches[g][::-1] for x in six]
  port = Port1(nd.mem)
  run_k(E, lk[0], lk[1], "keys", dict(x3=nd.SB + E.ENTRY * (q2 - 1) + E.KEYS0, x14=KN, x17=g, x18=0), nd, calls,
     cert.KEY_FIX + cert.KEY_INC * q2 + meter + snd.extra, port)
  ok = port.blocks == snd.msgs and state1_ok(nd.mem, snd)
  nb = len(port.blocks)
  for t0, n in ((0, 2), (8 * H1, 2)):
   run_k(E, rt[0], rt[1], "fin1", dict(x7=t0, x1=n), nd, calls, 1 + cert.FIN1_X * n, port)
  fins = snd.finals(D1S)
  ok &= port.blocks[nb:] == fins and all(nd.mem[8 * d] == blk1(d) + 32 and nd.mem[blk1(d) + 16] == 0 for d in D1S)
  msgs_ok &= ok; blocks += len(port.blocks)
  for m in port.blocks: arrivals[dest1(m[0])].append(m)
 R2END = Node.R2 + 3 * 32 * K2B
 ref = {h: [] for h in homes}; home_msgs = {h: [] for h in homes}; comp_ok = True
 for d1, nd in rels.items():
  for t0, n in ((0, 2), (H1 - 2, 2)):
   run_k(E, rt[0], rt[1], "init3", dict(x0=n, x1=LC2 + 8 * t0, x2=blk2(t0) + 32, x3=B2), nd, calls, cert.INIT3_BKT * n + 1)
   run_k(E, rt[0], rt[1], "tinit", dict(x0=n, x1=(LC2 + 8 * t0) | RC2, x2=Node.R2), nd, calls, cert.TINIT_ENTRY * n + 1)
  snd = Sender(lambda d: Node.R2, R2END, cert, Node.E2); port = Port2(nd.mem); arr = arrivals[d1]; rg.shuffle(arr); dn = []
  for m in arr:
   e0 = snd.extra
   for t in range(m[1] // 32):
    rec = list(m[4 + 4 * t:8 + 4 * t]); snd.put(rec[0] & M41, rec); ref[(d1, rec[0] & M41)].append(tuple(rec))
   dn.append(deliver_arr(E, rt[0], rt[1], nd, m, R2END, snd.extra - e0, calls, cert, port, snd.tag, G))
  comp_ok &= (dn == [False] * (len(arr) - 1) + [True] and nd.mem[DONE] == 1
       and h2s_state(H2S_REGS[(id(nd), R2END)]) == (1 << 126, G) and sum(m[2] for m in arr if m[1] < 32 * K1) == len(arr))
  ok = port.blocks == snd.msgs and not port.msgs and state_ok(nd.mem, snd); fport = Port2(nd.mem)
  for t0, n in ((0, 2), (H1 - 2, 2)):
   run_k(E, rt[0], rt[1], "fin", dict(x0=MB2, x7=LC2 + 8 * t0, x3=(LC2 + 8 * t0) | RC2, x8=blk2(t0) + 32, x1=n, x2=Node.E2B,
      x11=snd.tag, x6=R2END), nd, calls, snd.fin_count(range(t0, t0 + n)), fport)
  fins = [snd.final(d) for d in D2S]
  ok &= (fport.blocks == [m for m, _ in fins if m] and fport.msgs == [[LC2 + 8 * d, Node.E2B, c, 0, 0, 0] for d, (_, c) in zip(D2S, fins)])
  msgs_ok &= ok; blocks += len(port.blocks) + len(fport.blocks)
  for m in port.blocks + fport.blocks + fport.msgs: home_msgs[(d1, dest(m[0]))].append(m)
 slots_ok, tags_ok, join_ok, cands = True, True, True, 0
 CC = (UP - (Node.R2 >> 1)) % 2 ** 64
 for h, nd in homes.items():
  nd.mem[Node.E2] = Node.R2; home_state(nd, 64, CC); hj = HomeJoin(cert, CC, 64)
  for m in home_msgs[h]:
   if len(m) == 6: deliver(E, rt[0], rt[1], nd.mem, m, calls, cert, cnt=CNTE)
   else:
    f, okb = hj.block(m, nd.mem[Node.E2], 2 ** 47); join_ok &= okb
    deliver(E, rt[0], rt[1], nd.mem, m, calls, cert, fused=True, r2=(Node.R2, R2END), formula=f)
  F, F0 = nd.mem[Node.E2B], Node.R2
  start = F0 - (F0 - Node.R2) % (32 * K2B); expect = 1 + -(-(F - start) // (32 * K2B))
  comp_ok &= nd.mem[CNTE] == 1 and nd.mem[CNT] + 1 == len(home_msgs[h]) == expect
  got = [tuple(nd.mem[x + 8 * k] for k in range(4)) for x in range(Node.R2, F, 32)]
  slots_ok &= got == ref[h][:len(got)] and F == min(Node.R2 + 32 * len(ref[h]), R2END)
  tags_ok &= all(rec[0] & M41 == h[1] and rec[2] in gens and rec[3] >> 36 & 7 in (2, 4) for rec in got)
  join_ok &= ([(nd.mem[nd.CL + 16 * i], nd.mem[nd.CL + 16 * i + 8]) for i in range(len(hj.cum))] == hj.cum
     and all(nd.mem[HT + 8 * b] == v for b, v in hj.heads.items()))
  cands += len(hj.cum)
 groups = {}
 for g in gens:
  for six in batches[g]:
   for x in six: groups.setdefault(E.exact_key(x)[0], set()).add(g)
 cnt, form = {}, {}
 for e, c, f in calls: cnt[e], form[e] = cnt.get(e, 0) + c, form.get(e, 0) + f
 res = dict(generators=G, relays=len(rels), homes=len(homes), leaves=6 * q2 * G, shared_key_leaves=shared,
      cross_generator_key_groups=sum(len(v) >= 2 for v in groups.values()), home_candidates=cands,
      records_carry_class_labels=lab_ok, messages_equal_sender_model=msgs_ok, slots_equal_reference=slots_ok,
      source_tags_ok=tags_ok, joins_match_chain_model=join_ok, completion_counters_ok=comp_ok,
      kernel_calls=len(calls), counts=cnt, block_messages=blocks,
      formulas_ok=all(c == f for _, c, f in calls) and lab_ok and msgs_ok and slots_ok and tags_ok and join_ok and comp_ok and cands > 0)
 res["block_paths"] = block_paths(E, rt, lk, rg, cert)
 res["grouped_paths"] = grouped_paths(E, rt[0], rt[1], rg, cert)
 res["ls64"] = ls64_checks(E, rg)
 res["va52"] = va52_checks(E)
 res["phase4_round"] = phase4_round(E, rt, rg, cert)
 return res

def phase4_round(E, rt, rg, c):
 C, ok, N, B = [], True, Node, 0x8000000
 d1, d2, gs = rg.randrange(1 << 42), rg.randrange(1 << 41), [rg.randrange(15 << 64, 1 << 68), rg.randrange(1 << 64)]
 h, rl = N(defaultdict(int)), N(defaultdict(int))
 def run(nd, e, regs, f, base=None, msg=None):
  p = base and Port2(nd.mem, base)
  C.append((E.a64.a64_run(rt[0], rt[1], e, dict(regs, sp=N.SP, x30=E.SENTINEL, x28=1, x29=CNT),
        nd.mem if msg is None else LinkView(nd.mem, msg, []), dev=p and p.dev), f))
  return p and p.msgs
 def q(nd, m, cap=64):
  ok_ = nd.mem[CNT] + 1
  run(nd, "dlvq", dict(x0=LINK, x9=64, x10=B - 8, x11=B + 64 * cap), (c.DLVQ if nd.mem[B - 8] < B + 64 * cap else c.DLVQ_DROP) + c.RCV, msg=m)
  return nd.mem[CNT] == ok_
 R = [[rg.randrange(1 << 64), 0, g % 2 ** 64, g >> 64 | rg.randrange(1 << 60) << 4] for g in gs]
 for i in range(8): h.mem[N.R2 + 8 * i] = R[i // 4][i % 4]
 h.mem[N.CL], h.mem[N.CL + 8], rl.mem[B - 8], hd = N.R2, N.R2 + 32, B, 2 ** 44 + 8 * d2
 m1 = run(h, "reqh", dict(x0=MB2, x12=N.CL, x9=1, x10=hd, x13=B), c.REQH_C + 1, MB2)
 ok &= m1 == [[0, 0, r[2], r[3], hd, B + 256 * i] for i, r in enumerate(R)] and all(q(rl, m) for m in m1)
 m2 = run(rl, "fwdq", dict(x0=MB2 + 64, x1=B, x6=B + 128, x7=8 * d1), 2 * c.FWDQ + 3, MB2 + 64)
 rep = [[rg.randrange(1 << 64) for _ in range(32)] for _ in gs]
 for i, m in enumerate(m2):
  g = N(defaultdict(int)); g.mem[B - 8] = B
  ok &= (m[1] & 15) << 64 | m[0] == gs[i] and q(g, m) and [g.mem[B + 8 * k] for k in range(4)] == [gs[i] >> 64, R[i][3], B + 64 * i << 8, 8 * d1]
  for k in range(32): g.mem[B + 4096 + 8 * k] = rep[i][k]
  m3 = run(g, "sendr", dict(x0=MB2, x1=B + 4096, x9=8, x7=8 * d1, x8=g.mem[B + 16]), 8 * c.SENDR + 1, MB2)
  ok &= len(m3) == 8 and all(q(rl, m) for m in m3)
 m4 = run(rl, "fwdr", dict(x0=MB2, x1=B + 128, x6=B + 64 * 18), 16 * c.FWDR + 3, MB2)
 for m in m4: run(h, "dlv", dict(x0=LINK), c.DLV + c.RCV, msg=m)
 d = N(defaultdict(int)); d.mem[B - 8] = B + 4096
 ok &= ([h.mem[B + 8 * k] for k in range(64)] == rep[0] + rep[1] and len(m4) == 16 and q(d, [1] * 6)
     and d.mem[B] == 0 and all(x == y for x, y in C) and max(c.PH4_HOPS.values()) <= c.MSG)
 return dict(ok=ok, calls=len(C), hops=c.PH4_HOPS)

def va52_checks(E):
 a64 = E.a64
 prog, lab = a64.a64_load("ld:\n\tldr x2,[x1]\n\tret\n")
 out = {}
 for name, a in (("lower_top", 2 ** 52 - 8), ("upper_bottom", UP), ("top", 2 ** 64 - 8), ("lower_2^52", 2 ** 52),
     ("upper_minus_8", UP - 8), ("2^63", 2 ** 63)):
  try: a64.a64_run(prog, lab, "ld", dict(x1=a, x30=E.SENTINEL), defaultdict(int)); out[name] = "ran"
  except AssertionError as ex: out[name] = "rejected" if "52-bit" in str(ex) else "other: " + str(ex)
 out["ok"] = [out[k] for k in ("lower_top", "upper_bottom", "top")] == ["ran"] * 3 and all(out[k] == "rejected" for k in ("lower_2^52", "upper_minus_8", "2^63"))
 return out

def ls64_checks(E, rg):
 a64 = E.a64
 src = open(os.path.join(E.HERE, "kernels", "route.s")).read()
 t_prog, t_lab = a64.a64_load("ld:\n\tld64b x20,[x1]\n\tret\nst:\n\tst64b x16,[x2]\n\tret\nhi:\n\tld64b x24,[x1]\n\tret\nodd:\n\tld64b x21,[x1]\n\tret\novl:\n\tld64b x12,[x15]\n\tret\n")
 R2 = (Node.R2, Node.R2 + 4 * 32 * K2B)
 msg = [rg.getrandbits(64) for _ in range(4 + 4 * K2B)]
 def run(prog, labels, entry, regs, mem, lsr=([WINDOW, R2], [R2])):
  try: a64.a64_run(prog, labels, entry, dict(regs, x30=E.SENTINEL), mem, ls64=lsr); return "ran"
  except AssertionError as ex: return "LS64 rejected" if "LS64" in str(ex) else "other: " + str(ex)
 out, o = {}, {}
 a64.a64_run(t_prog, t_lab, "ld", dict(x1=LINK + 64, x30=E.SENTINEL), LinkView({}, msg, []), out_regs=o, ls64=([WINDOW], []))
 out["ld64b_window_slot_2_3"] = [o["x%d" % (20 + k)] for k in range(8)] == msg[8:16]
 mem = {}
 out["st64b_r2_aligned"] = run(t_prog, t_lab, "st", dict(x2=R2[0] + 32), mem)
 out["st64b_r2_written"] = all(mem.get(R2[0] + 32 + 8 * k) == 0 for k in range(8))
 out["ld64b_odd_slot"] = run(t_prog, t_lab, "ld", dict(x1=LINK + 32), LinkView({}, msg, []))
 out["ld64b_outside_regions"] = run(t_prog, t_lab, "ld", dict(x1=0x5000040), {})
 out["st64b_into_window"] = run(t_prog, t_lab, "st", dict(x2=LINK + 64), LinkView({}, msg, []))
 out["st64b_r2_misaligned"] = run(t_prog, t_lab, "st", dict(x2=R2[0]), {})
 out["ld64b_x24"] = run(t_prog, t_lab, "hi", dict(x1=LINK + 64), LinkView({}, msg, []))
 out["ld64b_x21"] = run(t_prog, t_lab, "odd", dict(x1=LINK + 64), LinkView({}, msg, []))
 out["ld64b_base_overlap"] = run(t_prog, t_lab, "ovl", dict(x15=LINK + 64), LinkView({}, msg, []))
 def fused(text, c, r2):
  prog, labels = a64.a64_load(text); mem = defaultdict(int); CC = (UP - (r2[0] >> 1)) % 2 ** 64
  mem[ST], mem[ST + 8], mem[ST + 16], mem[ST + 24] = 0x6200000, 0x6200000 + 1024, CC, -2 * CC % 2 ** 64
  mem[0x6600000] = c; m = [0, c, 32 * K2B, 0x6600000] + msg[4:]
  r = run(prog, labels, "dlvh", dict(x0=LINK, x28=1, x29=CNT), LinkView(mem, m, []), ([WINDOW, r2], [r2]))
  return r if r != "ran" else ("ran, R2 exact" if all(mem[c + 8 * w] == m[4 + w] for w in range(4 * K2B)) else "ran, R2 wrong")
 out["dlvh_fresh_block"] = fused(src, R2[0], R2)
 r2z = (R2[0] - 32, R2[1] - 32)
 out["dlvh_r2_at_0_mod_64"] = fused(src, r2z[0], r2z)
 mut = src.replace("6:\n\tadd x13,x14,x2,lsr #1\n\tDLVHS\n", "6:\n\tadd x13,x14,x2,lsr #1\n\tDLVHP\n")
 out["dlvh_mutant_ld64b_at_slot_1"] = fused(mut, R2[0], R2) if mut != src else "mutation not applied"
 hm = src.replace("\tH2R1\n\ttbz x16,#6,3f\n", "\tH2P\n\ttbz x16,#6,3f\n")
 hp, hl = a64.a64_load(hm)
 regs = dict(x10=LC2, x18=32 * K2B, x12=1 << 40, x3=PW, x29=DONE, x6=1 << 62, x14=1 << 62, x5=16)
 out["h2s_mutant_ld64b_at_slot_1"] = run(hp, hl, "h2s", regs, LinkView(defaultdict(int), [0, 32 * 9, 1, 0] + msg[4:4 + 4 * K1], []), ([WINDOW], [])) if hm != src else "mutation not applied"
 rej = ("ld64b_odd_slot", "ld64b_outside_regions", "st64b_into_window", "st64b_r2_misaligned", "ld64b_x24", "ld64b_x21",
    "ld64b_base_overlap", "dlvh_r2_at_0_mod_64", "dlvh_mutant_ld64b_at_slot_1", "h2s_mutant_ld64b_at_slot_1")
 out["ok"] = (out["ld64b_window_slot_2_3"] and out["st64b_r2_aligned"] == "ran" and out["st64b_r2_written"]
     and out["dlvh_fresh_block"] == "ran, R2 exact" and all(out[k] == "LS64 rejected" for k in rej))
 return out

def grouped_paths(E, prog, labels, rg, cert):
 GG, HB_T, NSLOT = 3, 6, 2
 NB = 1 << HB_T; SLOTB = 32 * K2B * 64; SIZES = ((300, 9), (520, 262), (250, 7))
 def run_case(stop, ncl=64):
  calls, nd = [], Node(defaultdict(int))
  HRB, AFB = 0x6300000, 0x6500000
  E2P = (nd.E2, nd.E2 + 0x800000)
  pf = lambda g, s: HRB + 8 * s if g == 0 else E2P[(g - 1) % 2] + 32 * s
  aff = {}
  for w in range(cert.ZW): nd.mem[HT + 8 * w] = 0x5A5A5A5A5A5A5A5A
  run_k(E, prog, labels, "zero16384", dict(x0=cert.ZW, x1=HT), nd, calls, cert.ZERO_PASS + 1)
  ok = all(nd.mem[HT + 8 * w] == 0 for w in range(cert.ZW))
  run_k(E, prog, labels, "init3", dict(x0=NSLOT, x1=HRB, x2=nd.R2, x3=SLOTB), nd, calls, cert.INIT3_BKT * NSLOT + 1)
  ok &= all(nd.mem[HRB + 8 * s] == nd.R2 + SLOTB * s for s in range(NSLOT))
  CC = (UP - (nd.R2 >> 1)) % 2 ** 64; home_state(nd, ncl, CC); hj = HomeJoin(cert, CC, ncl)
  ran, stopped = [], False
  rc = [nd.R2 + SLOTB * s for s in range(NSLOT)]; F = list(rc)
  blk, nb = [[None] * K2B for _ in range(NSLOT)], [0] * NSLOT
  for g in range(GG):
   if stopped: break
   ran.append(g); aff[g] = [(rg.getrandbits(64), rg.getrandbits(64)) for _ in range(3)]
   for k, (px, py) in enumerate(aff[g]): nd.mem[AFB + (g * 3 + k) * 16] = px; nd.mem[AFB + (g * 3 + k) * 16 + 8] = py
   streams = []
   for s in range(NSLOT):
    nrec = SIZES[g][s] + (0 if stop or s else 1); recs = []
    for _ in range(nrec): recs.append([rg.getrandbits(64), rg.getrandbits(7) << 57 | rg.getrandbits(HB_T), rg.getrandbits(64), rg.getrandbits(64)])
    if hj.seen:
     pk2, pk3 = hj.seen[rg.randrange(len(hj.seen))][:2]; recs[0][1] = pk3; recs[0][0] = (pk2 & DM) | rg.getrandbits(41)
    recs[1][1] = recs[0][1]; recs[1][0] = (recs[0][0] & DM) | rg.getrandbits(41); msgs = []
    for rec in recs:
     blk[s][nb[s]] = rec; nb[s] += 1
     if nb[s] == K2B: msgs.append([0, rc[s], 32 * K2B, pf(g, s)] + [w for r in blk[s] for w in r]); rc[s] += 32 * K2B; nb[s] = 0
    if nb[s]: msgs.append([0, rc[s], 32 * nb[s], pf(g, s)] + [w for r in blk[s] for w in r or [STALE] * 4])
    streams.append(msgs)
   order = [m for st in streams for m in st]; rg.shuffle(order)
   for m in order:
    s = [t for t in range(NSLOT) if pf(g, t) == m[3]][0]
    f, okb = hj.block(m, F[s], NB); ok &= okb
    deliver(E, prog, labels, nd.mem, m, calls, cert, fused=True, r2=(nd.R2, nd.R2 + NSLOT * SLOTB), formula=f)
    ok &= all(nd.mem[m[1] + 8 * w] == m[4 + w] for w in range(4 * K2B))
   ok &= all(nd.mem[HT + 8 * b] == hj.heads.get(b, 0) for b in range(NB))
   ok &= all(nd.mem[(a >> 1) + CC & (2 ** 64 - 1)] == hj.nextd[a] and nd.mem[((a >> 1) + CC & (2 ** 64 - 1)) + 8] == hj.recs[a][1] for a in hj.nextd)
   ok &= [(nd.mem[nd.CL + 16 * i], nd.mem[nd.CL + 16 * i + 8]) for i in range(len(hj.cum))] == hj.cum
   ok &= nd.mem[ST] == nd.CL + 16 * len(hj.cum) and len(hj.cum) <= ncl
   for s in range(NSLOT): F[s] = rc[s] + 32 * nb[s]; deliver(E, prog, labels, nd.mem, [0, E2P[g % 2] + 32 * s, F[s], 0, 0, 0], calls, cert)
   ok &= all(nd.mem[E2P[g % 2] + 32 * s] == F[s] for s in range(NSLOT))
   slot_addrs = [nd.R2 + SLOTB * s + 32 * k for s in range(NSLOT) for k in range((F[s] - (nd.R2 + SLOTB * s)) // 32)]
   ok &= sorted(hj.nextd) == sorted(slot_addrs)
   reached = []
   for b in range(NB):
    u = hj.heads.get(b, 0)
    while u: a = 2 * ((u - CC) % 2 ** 64); reached.append(a); u = hj.nextd[a]
   ok &= sorted(reached) == sorted(slot_addrs) and all([nd.mem[a + 8 * w] for w in range(4)] == hj.recs[a] for a in slot_addrs)
   if stop and g == 1: stopped = any(nd.mem[a + 8] == nd.mem[b + 8] and nd.mem[a] & DM == nd.mem[b] & DM for a, b in hj.cum)
  ok &= all(nd.mem[a + 8] == nd.mem[b + 8] and nd.mem[a] & DM == nd.mem[b] & DM for a, b in hj.cum)
  parent = {}
  def find(u):
   while parent.get(u, u) != u: u = parent[u]
   return u
  for a, b in hj.cum: parent[find(a)] = find(b)
  by_key = {}
  for k2, k3, a in hj.seen: by_key.setdefault((k3, k2 & DM), []).append(a)
  true_pairs = [(r[i], r[i + 1]) for r in by_key.values() for i in range(len(r) - 1)]
  if ncl > len(hj.cum): ok &= all(find(a) == find(b) for a, b in true_pairs) and len(true_pairs) > 0
  ok &= all(c == f_ for _, c, f_ in calls)
  old_x_ok = [(nd.mem[AFB + 16 * k], nd.mem[AFB + 16 * k + 8]) for k in range(3)] == aff[0]
  wide = sum(1 for _, k3, _ in hj.seen if k3 & 0x30)
  ok &= old_x_ok and wide > 0 and len({k3 & (NB - 1) for _, k3, _ in hj.seen}) > 16
  return dict(ok=ok, groups_run=ran, candidates=len(hj.cum), stopped=stopped,
     stop_pair_listed=any(nd.mem[a + 8] == nd.mem[b + 8] for a, b in hj.cum),
     true_pairs=len(true_pairs), list_skips=hj.nls0, old_x_served=old_x_ok,
     records_with_bucket_bits_4_5=wide, fast_slow_blocks=hj.paths, listings_direct_late=hj.lists)
 a, b, c = run_case(True), run_case(False), run_case(False, ncl=1)
 return dict(stop_case=a, full_case=b, full_list_case=c,
    ok=a["ok"] and b["ok"] and c["ok"] and a["groups_run"] == [0, 1] and b["groups_run"] == [0, 1, 2]
    and a["stopped"] and a["stop_pair_listed"] and a["candidates"] > 0 and b["candidates"] > 0 and c["candidates"] == 1 and c["list_skips"] > 0,
    old_x_served=a["old_x_served"] and b["old_x_served"] and c["old_x_served"])

def block_paths(E, rt, lk, rg, cert):
 calls, out = [], {}; nd = Node(defaultdict(int)); D1_ = (0, H1, 1, 2, H1 + 1)
 xs = [E.planted_x(rg, d, rg.randrange(4), j=t % 3)[0] for t, d in enumerate(D1_)]
 ds = [t for t, n in enumerate((534, 140, 128, 2, 0)) for _ in range(n)]; rg.shuffle(ds)
 q2 = len(ds) // 6; leaves = [[xs[t] for t in ds[6 * l:6 * l + 6]] for l in range(q2)]
 plant(E, nd.mem, nd.SB, leaves)
 for t0, n, first in ((0, 3, LO1 + 32), (8 * H1, 2, UP + 32)):
  run_k(E, rt[0], rt[1], "init3", dict(x0=n, x1=t0, x2=first, x3=B1), nd, calls, cert.INIT3_BKT * n + 1)
 snd = Sender1(cert); meter, fl = keys_reference(E, nd.mem, nd.SB, q2, 7, 0x70, snd, cert); port = Port1(nd.mem)
 run_k(E, lk[0], lk[1], "keys", dict(x3=nd.SB + E.ENTRY * (q2 - 1) + E.KEYS0, x14=0x70, x17=7, x18=0), nd, calls,
    cert.KEY_FIX + cert.KEY_INC * q2 + meter + snd.extra, port)
 ok = port.blocks == snd.msgs and state1_ok(nd.mem, snd); nb = len(port.blocks)
 for t0, n in ((0, 3), (8 * H1, 2)): run_k(E, rt[0], rt[1], "fin1", dict(x7=t0, x1=n), nd, calls, 1 + cert.FIN1_X * n, port)
 fins = snd.finals(sorted(D1_)); ok &= port.blocks[nb:] == fins
 ok &= [f[:2] for f in fl] == [(t % 3, D1_[t]) for l in range(q2 - 1, -1, -1) for t in ds[6 * l:6 * l + 6]]
 out["keys_fin1"] = dict(records=len(ds), blocks_sent=nb, finals=len(fins), empty_finals=sum(m[1] == 0 for m in fins),
        upper_half_blocks=sum(m[0] >= 8 * H1 for m in port.blocks), model_ok=ok)
 nd = Node(); D2_ = [rg.getrandbits(40) for _ in range(4)]
 sizes = [K1] * 12 + list(range(K1)) + [K1] * 20; rg.shuffle(sizes)
 dl = [D2_[t] for t, n in enumerate((560, 140, 50, 20)) for _ in range(n)]; rg.shuffle(dl)
 end, m = nd.R2 + 2 * 32 * K2B, {D2_[0]: 1, D2_[1]: 2, D2_[2]: 1, D2_[3]: 1}
 snd = Sender(lambda d: end - 32 * K2B * m[d], end, cert)
 for d in D2_: nd.mem[LC2 + 8 * d], nd.mem[(LC2 + 8 * d) | RC2] = blk2(d) + 32, end - 32 * K2B * m[d]
 k, arr = 0, []
 for n in sizes:
  recs = [[rg.getrandbits(23) << 41 | dl[(k + t) % len(dl)], rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64)] for t in range(n)]
  k += n; arr.append(recs)
 port = Port2(nd.mem); nfin = 0
 for chunk in arr:
  e0 = snd.extra
  for rec in chunk: snd.put(rec[0] & M41, rec)
  nfin += len(chunk) < K1
  deliver_arr(E, rt[0], rt[1], nd, [0, 32 * len(chunk), 1, 0] + [w for r in chunk for w in r] + [STALE] * 4 * (K1 - len(chunk)), end, snd.extra - e0, calls, cert, port)
 ok2 = port.blocks == snd.msgs and not port.msgs and state_ok(nd.mem, snd)
 ok2 &= h2s_state(H2S_REGS[(id(nd), end)]) == ((1 << 126) + nfin - len(arr), (1 << 68) - nfin) and nd.mem[DONE] == 0
 dm, wr = {d: {} for d in D2_}, {d: [] for d in D2_}
 for msg in port.blocks: unpack(dm[dest(msg[0])], msg, wr[dest(msg[0])])
 ok2 &= rule_ok(dm, wr, snd, None)
 out["h2s"] = dict(records=k, blocks_sent=len(port.blocks), blocks_dropped=snd.drops, end_stored_back=snd.cursor(D2_[0]) == end, model_and_rule_ok=ok2)
 out["two_group_finals"] = tg = two_group_finals(E, rt, rg, cert, calls)
 out["wide_counts"] = wide_counts(E, rt, cert, calls)
 out["counts_equal_formulas"] = all(c == f for _, c, f in calls)
 out["ok"] = (out["counts_equal_formulas"] and ok and ok2 and out["keys_fin1"]["empty_finals"] == 2
     and out["keys_fin1"]["upper_half_blocks"] > 2 and out["h2s"]["blocks_sent"] >= 1 and out["h2s"]["blocks_dropped"] >= 1
     and out["h2s"]["end_stored_back"] and tg["ok"] and out["wide_counts"])
 return out

def wide_counts(E, rt, cert, calls):
 nd, end = Node(), Node.R2 + 32 * K2B; key = (id(nd), end)
 run = lambda c: deliver_arr(E, rt[0], rt[1], nd, [0, 0, c, 0] + [STALE] * 4 * K1, end, 0, calls, cert, None, 0, 1 << 68)
 def put(D, f): H2S_REGS[key] = dict(H2S_REGS[key], x13=D % 2 ** 64, x14=D >> 64, x15=f % 2 ** 64, x5=f >> 64)
 ok = not run(1) and h2s_state(H2S_REGS[key]) == (1 << 126, (1 << 68) - 1)
 put((1 << 126) + (1 << 64), 1); ok &= not run(1) and h2s_state(H2S_REGS[key]) == (1 << 64, 0)
 put((1 << 126) - 1, 1)
 return ok and run(2) and nd.mem[DONE] == 1 and h2s_state(H2S_REGS[key]) == (1 << 126, 1 << 68)

def two_group_finals(E, rt, rg, cert, calls):
 nd = Node(); d, base = H1 - 1, nd.R2; end = base + 30 * 32 * K2B
 nd.mem[LC2 + 8 * d], nd.mem[(LC2 + 8 * d) | RC2] = blk2(d) + 32, base
 snd = Sender(lambda _d: base, end, cert); dm, wr, finals, recs = {}, [], [], []
 for n in (300, 300):
  group = [[rg.getrandbits(23) << 41 | d, rg.getrandbits(64), rg.getrandbits(64), rg.getrandbits(64)] for _ in range(n)]
  sent0 = len(snd.msgs); recs.append(group); port = Port2(nd.mem)
  for c in range(0, n, K1):
   e0 = snd.extra; chunk = group[c:c + K1]
   for rec in chunk: snd.put(d, rec)
   deliver_arr(E, rt[0], rt[1], nd, [0, 32 * len(chunk), 1, 0] + [w for r in chunk for w in r] + [STALE] * 4 * (K1 - len(chunk)), end, snd.extra - e0, calls, cert, port)
  okg = port.blocks == snd.msgs[sent0:] and not port.msgs and state_ok(nd.mem, snd); fport = Port2(nd.mem)
  run_k(E, rt[0], rt[1], "fin", dict(x0=MB2, x7=LC2 + 8 * d, x3=(LC2 + 8 * d) | RC2, x8=blk2(d) + 32, x1=1, x2=nd.E2, x6=end), nd, calls, snd.fin_count([d]), fport)
  fmsg, fcur = snd.final(d)
  okg &= fport.blocks == ([fmsg] if fmsg else []) and fport.msgs == [[LC2 + 8 * d, nd.E2, fcur, 0, 0, 0]]
  for m in port.blocks + fport.blocks + fport.msgs: unpack(dm, m, wr if len(m) > 6 else None)
  finals.append((fcur, okg, port.blocks))
 f1, f2, blk2_ = finals[0][0], finals[1][0], finals[1][2]; q = 300 // K2B * K2B
 resent = len(blk2_) >= 1 and blk2_[0][1] == base + 32 * q and blk2_[0][4:4 + 4 * (300 - q)] == [w for r in recs[0][q:300] for w in r]
 ok = (all(okg for _, okg, _ in finals) and rule_ok({d: dm}, {d: wr}, snd, {d: f2}) and resent
    and (f1, f2, dm[nd.E2]) == (base + 32 * 300, base + 32 * 600, f2)
    and snd.cursor(d) == nd.mem[(LC2 + 8 * d) | RC2] == base + 600 // K2B * 32 * K2B
    and all([dm.get(base + 32 * t + 8 * w) for w in range(4)] == recs[0][t] for t in range(300))
    and all([dm.get(f1 + 32 * t + 8 * w) for w in range(4)] == recs[1][t] for t in range(300)))
 return dict(ok=ok, first_final_rewritten_identically=resent)
