"""Tag-tracing step loop v2 (213 instructions per certified trip) and the ledger at that price.

Part S2 generates the new loop, validates the encodings of the scalar AArch64 subset it uses, and
        executes it trip by trip on a small simulator. At every trip it compares the tag A, the
        certificate, the selector, the next entry and the exit status with the definition of proof.md
        section 2 (A = sum over 48 limbs of UMULH(y, Flo) + y * Fhi mod 2^64).
        It runs on random synthetic entries with 128-bit F words, forced uncertified tags, budget-limited
        calls and window leaves (whose child pointers point to the leaf itself).
        Three changes from the parent's 217-instruction loop:
          (1) y is read SP-relative from the caller's frame, so x1 holds a resident limb;
          (2) x9/x10 are not used as scratch: the certificate and selector use x4 after the last
              accumulation, so x9 and x10 hold resident limbs. That gives 22 resident limbs and 26
              reloaded by 13 LDP (the parent had 19 resident and reloaded 29 with 14 LDP + 1 LDR);
          (3) the child is loaded straight into x0 and the null test is removed. A window leaf's eight
              child pointers point to the leaf itself, and the driver calls with
              budget = min(steps left, walk left, 256). A call that is not limited by a counter therefore
              ends with BUDGET after exactly 256 trips, the last of them at the leaf, and the driver
              re-bases with s = A >> 61 from the cursor.
Part G2 the parent's expected-work ledger, recomputed with 213 per step and the re-stated per-call
        driver and entry/exit counts. Same formula, same constants, same caps.
Standard library only. Optional JSON on stdin: {"seed": int, "entries": int, "starts": int}.
"""
import json
import random
import sys
from decimal import Decimal, ROUND_CEILING, getcontext
from pathlib import Path

try:
 cfg = json.loads(sys.stdin.read() or "{}")
except ValueError:
 cfg = {}
SEED = int(cfg.get("seed", 1))
NENT = int(cfg.get("entries", 64))
NSTART = int(cfg.get("starts", 64))
rng = random.Random(SEED)
base = Path("/snapshot/files")
if not (base / "target-profile.json").is_file():
 base = Path(__file__).resolve().parent
profile = json.loads((base / "target-profile.json").read_text())
r = int(profile["subgroup_order"], 16)

M64 = (1 << 64) - 1
MASK61 = (1 << 61) - 1
SP, ZR = 31, 32
FRAME = 112                       # kernel frame; y lives in the caller's frame at kernel SP + 112 + 8 j
YOFF = FRAME
RES = [1, 9, 10, 11, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]   # 22 resident limbs
NRES = len(RES)
WINDOW = 256


def f_off(j):                     # parent's entry layout (P = entry + 512)
 return -512 + 16 * j if j < 32 else 64 + 16 * (j - 32)


def build():
 pro = [("addi", SP, SP, -FRAME), ("str", 3, SP, 96)]
 for k in range(6):
  pro.append(("stp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k))
 for j in range(0, NRES, 2):
  pro.append(("ldp", RES[j], RES[j + 1], SP, YOFF + 8 * j))
 body = []
 first = True
 j = 0
 while j < 48:
  if j < NRES:
   ys = [RES[j]]
  else:
   body.append(("ldp", 7, 8, SP, YOFF + 8 * j))
   ys = [7, 8]
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
 body += [("addi", 4, 3, 96),
    ("andi", 4, 4, MASK61),
    ("subs", ZR, 4, 96),
    ("b.ls", "EXACT"),
    ("lsr", 4, 3, 61),
    ("ldrx", 0, 0, 4),            # P <- child[s]; a leaf's children are the leaf itself
    ("subsi", 2, 2, 1),
    ("b.ne", "LOOP")]

 def exit_path(code):
  return [("ldr", 4, SP, 96), ("stp", 0, 3, 4, 0), ("str", 2, 4, 16), ("movi", 0, code)] + \
    [("ldp", 19 + 2 * k, 20 + 2 * k, SP, 16 * k) for k in range(6)] + [("addi", SP, SP, FRAME), ("ret",)]
 return pro, body, {"BUDGET": exit_path(0), "EXACT": exit_path(1)}


GPR = set(range(31))


def validate(code):
 bad = []
 for ins in code:
  op, a = ins[0], ins[1:]
  if op in ("ldp", "stp"):
   ok = a[0] in GPR and a[1] in GPR and a[0] != a[1] and (a[2] in GPR or a[2] == SP) \
    and a[3] % 8 == 0 and -512 <= a[3] <= 504
  elif op in ("ldr", "str"):
   ok = a[0] in GPR and (a[1] in GPR or a[1] == SP) and a[2] % 8 == 0 and 0 <= a[2] <= 32760
  elif op == "ldrx":
   ok = all(x in GPR for x in a)
  elif op in ("umulh", "mul", "add"):
   ok = all(x in GPR for x in a)
  elif op == "madd":
   ok = all(x in GPR for x in a)
  elif op == "addi":
   ok = (a[0] in GPR or a[0] == SP) and (a[1] in GPR or a[1] == SP) and -4095 <= a[2] <= 4095
  elif op == "andi":
   ok = a[2] == MASK61
  elif op == "subs":
   ok = a[0] == ZR and a[1] in GPR and 0 <= a[2] <= 4095
  elif op == "subsi":
   ok = a[0] in GPR and a[1] in GPR and 0 <= a[2] <= 4095
  elif op == "lsr":
   ok = 0 <= a[2] <= 63
  elif op == "movi":
   ok = 0 <= a[1] < 65536
  else:
   ok = op in ("b.ls", "b.ne", "ret")
  regs_used = [x for x in a if isinstance(x, int) and x in GPR]
  if op not in ("ldp", "stp", "ldr", "str", "addi", "andi", "subs", "subsi", "lsr", "movi"):
   pass
  if 18 in [x for x in a[:2] if isinstance(x, int)] and op not in ("stp", "str", "ldp"):
   ok = False
  if not ok:
   bad.append(ins)
 return bad


def run(pro, body, exits, mem, regs):
 C = Z = 0
 trips = []

 def ex(ins):
  nonlocal C, Z
  op = ins[0]
  if op == "ldp":
   ad = regs[ins[3]] + ins[4]
   regs[ins[1]], regs[ins[2]] = mem[ad], mem[ad + 8]
  elif op == "stp":
   ad = regs[ins[3]] + ins[4]
   mem[ad], mem[ad + 8] = regs[ins[1]], regs[ins[2]]
  elif op == "ldr":
   regs[ins[1]] = mem[regs[ins[2]] + ins[3]]
  elif op == "str":
   mem[regs[ins[2]] + ins[3]] = regs[ins[1]]
  elif op == "ldrx":
   regs[ins[1]] = mem[regs[ins[2]] + 8 * regs[ins[3]]]
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
  elif op in ("subs", "subsi"):
   d = regs[ins[2]] - ins[3]
   C, Z = int(d >= 0), int(d == 0)
   if ins[1] != ZR:
    regs[ins[1]] = d & M64
  elif op == "lsr":
   regs[ins[1]] = regs[ins[2]] >> ins[3]
  elif op == "movi":
   regs[ins[1]] = ins[2]
  elif op == "ret":
   pass
  else:
   raise RuntimeError(op)
  regs[ZR] = 0

 total = 0
 for ins in pro:
  ex(ins)
  total += 1
 while True:
  n, target = 0, None
  for ins in body:
   n += 1
   if ins[0] == "b.ls":
    if (not C) or Z:
     target = ins[1]
     break
   elif ins[0] == "b.ne":
    if Z:
     target = "BUDGET"
     break
   else:
    ex(ins)
  trips.append(n)
  total += n
  if target is not None:
   for ins in exits[target]:
    ex(ins)
    total += 1
   return target, trips, total


def tag(y, F):
 A = 0
 for j in range(48):
  lo, hi = F[j]
  A = (A + ((y[j] * lo) >> 64) + y[j] * hi) & M64
 return A


def part_S2():
 pro, body, exits = build()
 viol = validate(pro + body + [i for e in exits.values() for i in e])
 # synthetic window: entries with depth; children of depth-256 entries point to themselves
 base_addr = 0x10000000
 P = [base_addr + 4096 * i + 512 for i in range(NENT)]
 depth = [1 + (i * 37) % 256 for i in range(NENT)]
 depth[NENT - 1] = 256
 depth[NENT - 2] = 256
 mem, Fw = {}, []
 for i in range(NENT):
  F = [(rng.getrandbits(64), rng.getrandbits(64)) for _ in range(48)]
  Fw.append(F)
 kids = {}
 for i in range(NENT):
  for s in range(8):
   if depth[i] == 256:
    kids[(i, s)] = i
   else:
    kids[(i, s)] = rng.randrange(NENT)
 y = [rng.getrandbits(64) for _ in range(48)]
 y[47] |= 1
 forced = set(i for i in range(NENT) if i % 9 == 4)
 for i in forced:                        # force (A + 96) mod 2^61 <= 96: uncertified
  F = Fw[i]
  A = tag(y, F)
  target = (rng.getrandbits(3) << 61) | ((1 << 61) - 96 + rng.randrange(0, 96))
  lo, hi = F[47]
  F[47] = (lo, (hi + (target - A) * pow(y[47], -1, 1 << 64)) & M64)
 for i in range(NENT):
  for j, (lo, hi) in enumerate(Fw[i]):
   mem[P[i] + f_off(j)] = lo
   mem[P[i] + f_off(j) + 8] = hi
  for s in range(8):
   mem[P[i] + 8 * s] = P[kids[(i, s)]]
 SPV = 0x40000000
 for j in range(48):
  mem[SPV + 8 * j] = y[j]                 # caller frame holds Y0 at caller SP; kernel SP + 112
 CUR = 0x30000000
 mism, checked, lens, statuses, leaf_ends, extra = 0, 0, set(), {}, 0, []
 for k in range(NSTART):
  start = rng.randrange(NENT)
  budget = rng.choice([256, 256, 256, rng.randrange(1, 257), 1, 2])
  regs = [rng.getrandbits(64) for _ in range(33)]
  regs[0], regs[2], regs[3], regs[SP] = P[start], budget, CUR, SPV
  snap = dict(mem)
  status, trips, total = run(pro, body, exits, mem, regs)
  statuses[status] = statuses.get(status, 0) + 1
  # replay with the definition
  e, n_done, exp_status = start, 0, "BUDGET"
  for t in range(len(trips)):
   A = tag(y, Fw[e])
   cert = (((A + 96) & M64) & MASK61) > 96
   checked += 1
   if not cert:
    exp_status = "EXACT"
    break
   e = kids[(e, A >> 61)]
   n_done += 1
   if t < len(trips) - 1:
    lens.add(trips[t])
   if n_done == budget:
    break
  mism += status != exp_status
  mism += mem[CUR] != P[e] or mem[CUR + 8] != A or mem[CUR + 16] != budget - n_done
  if status == "BUDGET" and depth[e] == 256 and kids[(e, A >> 61)] == e and Fw[e]:
   leaf_ends += 1
  # callee-saved registers and SP restored, y untouched
  mism += regs[SP] != SPV or any(mem[SPV + 8 * j] != y[j] for j in range(48))
  if status == "BUDGET" and trips:
   lens.add(trips[-1])
 return {"encoding_violations": len(viol), "first_violations": [str(v) for v in viol[:3]],
   "prologue": len(pro), "body_length": len(body), "exit_path": max(len(e) for e in exits.values()),
   "trip_instructions_certified": sorted(lens), "resident_limbs": NRES,
   "reload_loads_per_trip": sum(1 for i in body if i[0] == "ldp" and i[3] == SP),
   "entries": NENT, "calls": NSTART, "trips_checked": checked, "exit_statuses": statuses,
   "calls_ending_on_a_leaf_self_pointer": leaf_ends,
   "mismatches_with_definition": mism,
   "body_opcodes": {op: sum(1 for i in body if i[0] == op) for op in sorted({i[0] for i in body})}}


# ------------------------------------------------------------------ Part R: full driver vs counted steps
def part_R(W=4, steps_total=6000, calls_cap=10**6):
 """Window W (production 256) multiset tree over 8 multipliers, the generated kernel, and the driver of
    machine.md (budget = min(steps left, W + 1 - d0); BUDGET with no counter limit = leaf trip = re-base).
    Compared step by step with the parent's semantics (null leaf children, REBASE exit at the leaf trip,
    EXACT exit, every certified or exact trip counted as one mathematical step Y <- Y * M_s)."""
 pro, body, exits = build()
 from itertools import combinations_with_replacement
 ents = [()]
 for k in range(1, W + 1):
  ents += list(combinations_with_replacement(range(8), k))
 idx = {e: n for n, e in enumerate(ents)}
 base_addr = 0x10000000
 P = [base_addr + 4096 * n + 512 for n in range(len(ents))]
 Fw = [[(rng.getrandbits(64), rng.getrandbits(64)) for _ in range(48)] for _ in ents]
 for n in range(7, len(ents), 23):                # all-zero F words: A = 0, uncertified for every y (EXACT)
  Fw[n] = [(0, 0)] * 48
 mem = {}
 for n, e in enumerate(ents):
  for j, (lo, hi) in enumerate(Fw[n]):
   mem[P[n] + f_off(j)] = lo
   mem[P[n] + f_off(j) + 8] = hi
  for s in range(8):
   child = tuple(sorted(e + (s,))) if len(e) < W else None
   mem[P[n] + 8 * s] = P[idx[child]] if child is not None else P[n]      # leaf: self-pointer
 SPV, CUR = 0x40000000, 0x30000000

 def new_y(y, e, s):                     # stands for the re-base Y0 <- Y0 * T_e (value irrelevant here)
  rr = random.Random(hash((tuple(y[:4]), e, s)) & M64)
  return [rr.getrandbits(64) for _ in range(48)]

 def reference(y0, L):
  y, e, log = y0, (), []
  for _ in range(L):
   A = tag(y, Fw[idx[e]])
   if (((A + 96) & M64) & MASK61) <= 96:      # uncertified: exact evaluation, re-base, root child
    s = A >> 61
    y, e = new_y(y, e, ("x", s)), (s,)
    log.append(("X", s))
   elif len(e) < W:
    s = A >> 61
    e = tuple(sorted(e + (s,)))
    log.append(("C", s))
   else:                                       # leaf trip: REBASE, step into the root's child s
    s = A >> 61
    y, e = new_y(y, e, s), (s,)
    log.append(("R", s))
  return log

 def driver(y0, L):
  y, e, d0, left, log, calls = y0, (), 0, L, [], 0
  while left > 0 and calls < calls_cap:
   calls += 1
   for j in range(48):
    mem[SPV + 8 * j] = y[j]
   budget = min(left, W + 1 - d0)
   regs = [0] * 33
   regs[0], regs[2], regs[3], regs[SP] = P[idx[e]], budget, CUR, SPV
   status, trips, _ = run(pro, body, exits, mem, regs)
   n_cert = budget - mem[CUR + 16]             # certified trips taken in this call
   ent = ents[P.index(mem[CUR])]
   A = mem[CUR + 8]
   # replay the certified trips for the log (the kernel's own path)
   path_e = e
   for k in range(n_cert):
    At = tag(y, Fw[idx[path_e]])
    s = At >> 61
    if len(path_e) < W:
     path_e = tuple(sorted(path_e + (s,)))
     log.append(("C", s))
    else:
     log.append(("R", s))
   left -= n_cert
   if status == "EXACT":
    s = A >> 61
    y, e, d0 = new_y(y, ent, ("x", s)), (s,), 1
    log.append(("X", s))
    left -= 1
   elif n_cert == W + 1 - d0:                  # not counter-limited: the last trip was the leaf trip
    assert ent == path_e and len(ent) == W
    s = A >> 61
    y, e, d0 = new_y(y, ent, s), (s,), 1
   else:                                       # counter-limited: the run (here: the test) ends
    e = ent
    d0 = len(ent)
  return log

 y0 = [rng.getrandbits(64) for _ in range(48)]
 ref = reference(y0, steps_total)
 got = driver(y0, steps_total)
 return {"window_test": W, "entries": len(ents), "steps": steps_total,
   "reference_steps": len(ref), "driver_steps": len(got),
   "identical_step_sequence": ref == got[:len(ref)] and len(got) == len(ref),
   "rebases": sum(1 for x in ref if x[0] == "R"), "exact": sum(1 for x in ref if x[0] == "X")}


# ------------------------------------------------------------------ Part D: deterministic DP store
# Key = the FULL canonical distinguished element Y0 (after the exact evaluation's re-base): its 12 Fp
# coordinates as canonical Montgomery residues, 4 limbs each, 48 limbs = 3072 bits; limb (k, i) occupies
# key bits 256 k + 64 i .. +63. Montgomery form is a bijection of Fp, so key equality is element equality.
KEY_BITS = 3072
LEAF_BYTES, INODE_BYTES = 512, 32       # leaf: key 384 + labels a, b 64 + walk id 8 + flags 8 + pad 48
TRAV_I = 16          # per internal node: lsr (limb index), ldr key limb [reg, lsl 3], lsr var, and #1,
      # ldr child [reg, lsl 3], tag test, branch, loop: 10 itemised, 16 charged
CMP_I = 4 * 48       # 24 LDP pairs from each side, EOR, ORR-accumulate, final CBZ: <= 4 per limb
DIFF_I = 4 * 48 + 8  # scan to the first differing limb (EOR, CBNZ, index), CLZ, bit index arithmetic
LINK_I = 16          # allocate internal node (index bump), 3 stores, parent child store, tag
LEAF_W_I = 96        # 24 LDP + 24 STP of the key (materialisation of the DP element), labels, id, flags
CALL_I = 64
STORE_OP_CHARGE = 2 ** 17
OPS_CAP_LOG2 = 52    # one store operation per walk; walk cap 2^52
TABLE_ENTRIES, ENTRY_BYTES = 525783425977953, 2624


class CritBit:
 """Binary Patricia (crit-bit) trie keyed by the full 3072-bit element. Internal node: (bit, child0,
    child1); leaf: (key, record). Branch bit positions strictly decrease from the root, so a path has at
    most 3072 internal nodes. Equal key (the same group element): the existing record is returned and kept."""

 def __init__(self):
  self.root = None
  self.leaves, self.inodes = [], []
  self.max_depth = 0
  self.max_instr = 0

 def _walk(self, t):
  node, depth = self.root, 0
  while node[0] == "I":
   b, c0, c1 = self.inodes[node[1]]
   node = c1 if (t >> b) & 1 else c0
   depth += 1
  return node[1], depth

 def insert_or_find(self, t, rec):
  instr = CALL_I
  if self.root is None:
   self.leaves.append((t, rec))
   self.root = ("L", 0)
   self.max_instr = max(self.max_instr, instr + LEAF_W_I)
   return None
  li, d1 = self._walk(t)
  instr += TRAV_I * d1 + CMP_I
  lt, lrec = self.leaves[li]
  if lt == t:
   self.max_instr = max(self.max_instr, instr)
   self.max_depth = max(self.max_depth, d1)
   return lrec
  crit = (lt ^ t).bit_length() - 1
  instr += DIFF_I
  newleaf = ("L", len(self.leaves))
  self.leaves.append((t, rec))
  parent, side, node, d2 = None, None, self.root, 0
  while node[0] == "I" and self.inodes[node[1]][0] > crit:
   b, c0, c1 = self.inodes[node[1]]
   parent, side = node[1], (t >> b) & 1
   node = c1 if side else c0
   d2 += 1
  nb = (t >> crit) & 1
  self.inodes.append((crit, node, newleaf) if nb else (crit, newleaf, node))
  newi = ("I", len(self.inodes) - 1)
  if parent is None:
   self.root = newi
  else:
   b, c0, c1 = self.inodes[parent]
   self.inodes[parent] = (b, c0, newi) if side else (b, newi, c1)
  instr += TRAV_I * d2 + LINK_I + LEAF_W_I
  self.max_depth = max(self.max_depth, d1, d2)
  self.max_instr = max(self.max_instr, instr)
  return None


def element_key(coords):
 """12 canonical residues (< p) -> 3072-bit key (coordinate k at bits 256 k .. 256 k + 255)."""
 return sum(c << (256 * k) for k, c in enumerate(coords))


def part_D():
 st, ref = CritBit(), {}
 pmod = int(profile["field_modulus"], 16)
 rnd = lambda: [rng.randrange(pmod) for _ in range(12)]
 elems = [rnd() for _ in range(600)]
 base_e = rnd()
 # adversarial: single-bit neighbours of one element at every limb boundary of every coordinate
 flips = []
 for k in range(12):
  for b in (0, 1, 62, 63, 64, 65, 127, 128, 191, 192, 252, 253):
   e = list(base_e)
   e[k] ^= 1 << b
   if e[k] < pmod:
    flips.append(e)
 # nested shared prefixes: elements agreeing with base_e on the top coordinates and more and more bits
 nested = []
 for k in (11, 6, 0):
  for b in range(253, -1, -23):
   e = list(base_e)
   e[k] = (e[k] ^ (1 << b)) % pmod
   for kk in range(k):
    e[kk] = rng.randrange(pmod)
   nested.append(e)
 edge = [[0] * 12, [1] + [0] * 11, [pmod - 1] * 12, [0] * 11 + [pmod - 1]]
 keys = [element_key(e) for e in elems + [base_e] + flips + nested + edge]
 keys += keys[:150] + keys[600:650]                        # duplicates (exact element equality)
 rng.shuffle(keys)
 bad = 0
 for n, k in enumerate(keys):
  rec = (n, k % 1000003)
  got = st.insert_or_find(k, rec)
  exp = ref.get(k)
  if exp is None:
   ref[k] = rec
  bad += got != exp
 for k, rec in ref.items():
  li, _ = st._walk(k)
  bad += st.leaves[li] != (k, rec)
 distinct = len(ref)
 worst = CALL_I + 2 * KEY_BITS * TRAV_I + CMP_I + DIFF_I + LINK_I + LEAF_W_I
 N = 2 ** OPS_CAP_LOG2
 bytes_store = N * LEAF_BYTES + (N - 1) * INODE_BYTES
 table_bytes = TABLE_ENTRIES * ENTRY_BYTES
 peak = table_bytes + bytes_store + 2 ** 30
 lg = lambda v: float(Decimal(v).ln() / Decimal(2).ln())
 lgup = lambda v: float((Decimal(v).ln() / Decimal(2).ln() * 10000).to_integral_value(rounding=ROUND_CEILING) / 10000)
 return {"key_bits": KEY_BITS, "operations": len(keys), "distinct_keys": distinct, "mismatches": bad,
   "adversarial_flip_keys": len(flips), "nested_prefix_keys": len(nested),
   "leaves": len(st.leaves), "internal_nodes": len(st.inodes),
   "nodes_ok": len(st.leaves) == distinct and len(st.inodes) == distinct - 1,
   "max_depth_seen": st.max_depth, "max_instr_seen": st.max_instr,
   "depth_bound": KEY_BITS,
   "analytic_worst_instr_per_op": worst, "charged_per_op": STORE_OP_CHARGE,
   "charge_ok": worst <= STORE_OP_CHARGE and st.max_instr <= worst,
   "store_ops_cap_log2": OPS_CAP_LOG2, "store_total_log2_instr": OPS_CAP_LOG2 + 17,
   "store_bytes_log2": lg(bytes_store), "table_bytes_log2": lg(table_bytes),
   "peak_bytes_log2_rounded_up": lgup(peak), "peak_below_2^63": peak < 2 ** 63}


# ------------------------------------------------------------------ Part G2: ledger
STEP = 213
TT_ENTRY = 19 + 12          # prologue (2 + 6 STP + 11 LDP) + one exit path
DRIVER = 168                # parent's 160 plus the leaf test (compare budget with 256, branch) and slack
MUL_LEN = 8679              # parent's generated window product (lazy.py part L)
LABELS = 2 * 34 + 20
WEND = 9216
FIXED = 2 ** 110
DP_BITS = 84
CLAIM = Decimal("119.65")


def part_G2(allow_percent=5, cap_bits=133):
 getcontext().prec = 80
 ln2 = Decimal(2).ln()
 lg = lambda v: Decimal(v).ln() / ln2

 def arctan_inv(n):
  x = Decimal(1) / n
  total, term, k, n2 = x, x, 1, n * n
  while abs(term) > Decimal(10) ** -75:
   term = -term / n2
   k += 2
   total += term / k
  return total
 PI = 4 * (4 * arctan_inv(5) - arctan_inv(239))
 wend_exact = MUL_LEN + LABELS + TT_ENTRY + DRIVER
 kappa = Decimal(8) / 7
 model = (PI * kappa * Decimal(r) / 2).sqrt()
 E_touch = Decimal(100 + allow_percent) / 100 * model + 1          # discrete tail: + 1
 S_cap = Decimal(2) ** cap_bits
 p_self = Decimal(2) ** (52 + 2 * DP_BITS + 1) / Decimal(r)       # some walk closes a cycle before its DP
 p_no_dp = Decimal(2) ** 52 * Decimal(-64).exp()                  # some walk has no DP within 2^90 steps
 E_steps = E_touch + Decimal(2) ** (DP_BITS + 1) + (p_self + p_no_dp) * S_cap   # exceptions charged S_cap
 E_walks = E_steps / Decimal(2) ** DP_BITS * 2 + 1
 E_exact = E_steps * 97 / Decimal(2) ** 61
 losses = {
  "step_cap_markov": E_steps / S_cap,
  "walk_cap_2^52_markov": E_walks / Decimal(2) ** 52,
  "exact_meter_2^90_markov": E_exact / Decimal(2) ** 90,
  "walk_length_cap_and_self_cycles": p_self + p_no_dp,
  "degeneracy_and_sampling": Decimal(2) ** -200,
 }
 fail = sum(losses.values())
 succ = Decimal("0.97")
 K = int((Decimal("0.39") / succ * Decimal(2) ** 64).to_integral_value(rounding=ROUND_CEILING))
 gate = Decimal(K) / Decimal(2) ** 64
 R_run = STEP * E_steps + WEND * (E_steps / WINDOW + E_walks + 1) + FIXED
 W_run = R_run / 20000
 total = 1 + gate * W_run
 every_run = (STEP * S_cap + WEND * (S_cap / WINDOW + Decimal(2) ** 52 + 1) + FIXED) / 20000
 parent_total = 1 + gate * ((217 * E_steps + WEND * (E_steps / WINDOW + E_walks + 1) + FIXED) / 20000)
 up = lambda v, d: float((lg(v) * 10 ** d).to_integral_value(rounding=ROUND_CEILING) / 10 ** d)
 return {"window_end_instr": wend_exact, "window_end_charged": WEND, "window_ok": wend_exact <= WEND,
   "per_step_instr": str(STEP + Decimal(WEND) / WINDOW),
   "log2_charged_mean_steps": float(lg(E_steps)),
   "losses": {k: float(v) for k, v in losses.items()},
   "success_past_gate_lower": float(1 - fail), "budget_ok": 1 - fail >= succ,
   "gate_K": str(K), "success_ok": gate * succ >= Decimal("0.39"),
   "K_is_least": Decimal(K - 1) / Decimal(2) ** 64 * succ < Decimal("0.39"),
   "log2_expected_run_units": float(lg(W_run)), "log2_expected_total_units": float(lg(total)),
   "parent_log2_expected_total_units_same_formula": float(lg(parent_total)),
   "improvement_bits_same_formula": float(lg(parent_total) - lg(total)),
   "claim": str(CLAIM), "claim_ok": CLAIM >= lg(total),
   "claim_minus_computed_bits": float(CLAIM - lg(total)),
   "log2_every_run_bound_units": float(lg(every_run)),
   "displayed_upper_bounds_rounded_up": {"E_steps_log2": up(E_steps, 4), "E_W_run_log2": up(W_run, 4),
             "E_W_log2": up(total, 4), "every_run_log2": up(every_run, 4)}}


S2 = part_S2()
R2 = part_R(4, int(cfg.get("driver_steps", 6000)))
D2 = part_D()
G2 = part_G2()
ok = (R2["identical_step_sequence"] and R2["exact"] > 0 and R2["rebases"] > 0 and D2["mismatches"] == 0 and D2["nodes_ok"] and D2["charge_ok"]
  and D2["max_depth_seen"] <= KEY_BITS and D2["peak_below_2^63"] and S2["encoding_violations"] == 0 and S2["mismatches_with_definition"] == 0
  and S2["trip_instructions_certified"] == [STEP] and set(S2["exit_statuses"]) == {"BUDGET", "EXACT"}
  and S2["calls_ending_on_a_leaf_self_pointer"] > 0
  and G2["window_ok"] and G2["budget_ok"] and G2["success_ok"] and G2["K_is_least"] and G2["claim_ok"])
out = {"seed": SEED, "partS2_step_loop_v2": S2, "partR_driver_equivalence": R2, "partD_dp_store": D2,
  "partG2_ledger": G2, "all_ok": ok,
  "scope": "The step loop is executed on a simulator of the scalar AArch64 subset it uses, against the tag "
    "definition; the window product, toy collision times and base checks are the parent's (lazy.py, "
    "development-evidence*.json). No production DLP is run."}
print(json.dumps(out, indent=1, default=str))
sys.exit(0 if ok else 1)
