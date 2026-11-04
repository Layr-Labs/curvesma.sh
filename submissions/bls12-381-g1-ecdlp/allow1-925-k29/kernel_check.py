"""Audit of measure/side_lanes.S (proof.md s1/s24): 24 paths (worst 932, charged mean 925); word-level calls vs
exact arithmetic (edge/random/curve/extreme lanes, branches both ways, fixed == FIXED, X-counter exact every call,
inputs/callee-saved unchanged); stride-128 footprint under canaries (96-stride negative control)."""
import sys, json, re, random, time
from pathlib import Path
root = Path(__file__).resolve().parent
if (root / 'files').is_dir(): root = root / 'files'
sys.path.insert(0, str(root))
import a64
params = json.loads((root / 'params.json').read_text()); p = int(params['p'], 16)
R = pow(2, 384, p); RI = pow(R, -1, p); MASK = (1 << 64) - 1; BASE = 0x80000
CBASE = 0x90000   # test address of the node-local X counter (BSS in production)
FRSZ = 432      # frame bytes; run() re-parses and asserts
WORST = 932
FIXED = 42      # measured fixed part per call: 32 (s20 row) + 10 s24 (init 1 + aggregate 9)
CHARGED = 40
assert FIXED - 10 + 1 + 5 + 2 == CHARGED   # exit + staging 5 + round control 2 (s20); s24's 10 in its own row

def lane_paths(text):
 lines = text.split("\n"); st = lines.index(".Llane:"); code, labels = [], {}
 for l in lines[st:]:
  l = l.split("//")[0].rstrip()
  if l.startswith(".Lexit:"): break
  if re.match(r"^\.L\w+:", l): labels[l[:-1]] = len(code); continue
  if l.startswith("\t") and l.strip(): code.append(l.strip())
 paths = []
 def walk(i, n, tr):
  while True:
   ins = code[i]; n += 1; op = ins.split()[0]
   if op == "cbnz": paths.append((n, tr)); return
   m = re.match(r"(b\.\w+|tbz|tbnz|cbz|b)\b.*?(\.L\w+)$", ins)
   if m and op == "b": assert labels[m.group(2)] > i; i = labels[m.group(2)]; continue
   if m: assert labels[m.group(2)] > i; walk(labels[m.group(2)], n, tr + [op + "+"]); tr = tr + [op + "-"]
   i += 1
 walk(0, 0, [])
 return paths

def prepare(text):
 code, data = text.split('.section .rodata', 1)
 vals = []
 for line in data.splitlines():
  if '.quad' in line: vals += [int(x.strip(), 0) for x in line.split('.quad', 1)[1].split(',')]
 assert sum(vals[i] << (64 * i) for i in range(6)) == p and vals[6] == (-pow(p, -1, 1 << 64)) & MASK
 code = code.replace('__SIZEOF_POINTER__', '8')
 code = re.sub(r'(?<![A-Za-z_])(-?\d+)\*(\d+)', lambda m: str(int(m[1]) * int(m[2])), code)
 code, n = re.subn(r'\badrp\s+x0,\.Lsk_consts', 'mov x0,#' + str(BASE), code); assert n == 1
 code, n = re.subn(r'\badd\s+x0,x0,:lo12:\.Lsk_consts', 'add x0,x0,#0', code); assert n == 1
 code, n = re.subn(r'\badrp\s+x0,x_count_words', 'mov x0,#' + str(CBASE), code); assert n == 1
 code, n = re.subn(r'\badd\s+x0,x0,:lo12:x_count_words', 'add x0,x0,#0', code); assert n == 1
 code = '\n'.join(l for l in code.splitlines() if not l.lstrip().startswith('#'))
 pr, labels = a64.a64_load(code)
 return pr, labels, {BASE + 8 * i: v for i, v in enumerate(vals)}

def put(m, a, v):
 for i in range(6): m[a + 8 * i] = (v >> (64 * i)) & MASK
def get(m, a): return sum(m.get(a + 8 * i, 0) << (64 * i) for i in range(6))

def add(A, B):
 x, y = A; u, v = B
 t = (3 * x * x * pow(2 * y, -1, p) if x == u else (v - y) * pow(u - x, -1, p)) % p
 xx = (t * t - x - u) % p
 return xx, (t * (x - xx) - y) % p

def expect(x1, y1, x2, y2, iv):   # residues in, residues out (kernel formula)
 mu = (y1 + y2) * iv * RI % p; xx = (mu * mu * RI - x1 - x2) % p
 yy = (mu * ((xx - x2) % p) * RI + y2) % p
 e = yy if yy % 2 == 0 else p - yy
 return xx, yy, [(e >> (64 * j)) & MASK for j in (1, 2, 3)]

def stride_audit(prog, labels, cm, rng, stride, n=8, sp=0x7000000):
 """Production slot layout (proof.md section 8): n out slots contiguous at the given stride from a 16-aligned base,
    every slot word and a 128-byte trailing region pre-filled with random canaries. Asserts: every word the kernel writes
    outside its own stack frame lies inside the union of the slots; every word outside the kernel's output ranges
    [out, out+120) keeps its canary (pads, neighbours, trailing region, inputs); outputs are correct. At stride 128 this
    passes; at stride 96 the last lane's key lands in the trailing region, so the first assertion fails (negative control)."""
 base, m, want = 0x20000, dict(cm), []
 for i in range(n):
  x1, y1, x2, y2, iv = (rng.randrange(p) for _ in range(5)); y2 = y2 or 1
  out, pp, tt, ii = base + stride * i, 0x30000 + 256 * i, 0x40000 + 256 * i, 0x50000 + 256 * i
  for j, v in enumerate((out, pp, tt, ii)): m[0x10000 + 32 * i + 8 * j] = v
  put(m, pp, x1); put(m, pp + 48, y1); put(m, tt, x2); put(m, tt + 48, y2)
  put(m, tt + 96, p - x2); put(m, tt + 144, p - y2); put(m, ii, iv)
  want.append(expect(x1, y1, x2, y2, iv))
 for a_ in range(base, base + stride * n + 128, 8): m[a_] = rng.getrandbits(64)   # canaries: all slots + trailing region
 m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0
 m[CBASE] = m[CBASE + 8] = 0
 regs = {f'x{i}': rng.getrandbits(64) for i in range(31)}; regs.update(x0=0x10000, x1=n, sp=sp)
 mem = dict(m); ro = {}
 z = [i for i, (mn, o) in enumerate(prog) if (mn, tuple(o)) == ('b.lo', ('.Lx_done',))][0] + 4; vis = {z: 0}
 a64.a64_run(prog, labels, 'side_lanes', regs, mem, limit=50000, out_regs=ro, visits=vis)
 assert mem[CBASE] == vis[z] and mem[sp - FRSZ + 240] == vis[z]   # counter exact here too
 slots, outs = set(), set()
 for i in range(n):
  slots.update(range(base + stride * i, base + stride * (i + 1), 8)); outs.update(range(base + stride * i, base + stride * i + 120, 8))
 written = {a_ for a_ in mem if a_ not in m or mem[a_] != m[a_]} - {CBASE, CBASE + 8}   # BSS counter outside slots
 frame = {a_ for a_ in written if sp - 1024 <= a_ < sp}                 # the kernel's own stack frame
 assert written - frame <= slots, sorted(hex(a_) for a_ in (written - frame) - slots)   # footprint inside the slots
 for a_, v in m.items():
  if a_ in (CBASE, CBASE + 8): continue   # BSS asserted separately
  assert mem[a_] == v or a_ in outs, hex(a_)     # pads, trailing region, inputs untouched
 for i, (xx, yy, key) in enumerate(want):
  o = base + stride * i
  assert get(mem, o) == xx and get(mem, o + 48) == yy and [mem[o + 96 + 8 * j] for j in range(3)] == key, i
 assert ro['sp'] == sp
 return n

def run(cfg):
 t0 = time.time(); rng = random.Random(cfg.get('seed', 1))
 text = (root / 'measure/side_lanes.S').read_text()
 paths = lane_paths(text); worst = max(n for n, _ in paths)
 assert worst == WORST, worst
 assert sorted(n for n, _ in paths) == sorted([918] * 6 + [922] * 6 + [928] * 6 + [932] * 6)
 assert re.findall(r'\[sp,#-(\d+)\]!', text) == ['288'] and re.findall(r'sub\s+sp,sp,#(\d+)', text) == ['144']
 prog, labels, cm = prepare(text)
 flat = [(mn, tuple(o)) for mn, o in prog]
 init = [i for i, t in enumerate(flat) if t == ('str', ('xzr', '[sp,#240]'))]; assert len(init) == 1
 xarm = [i for i, t in enumerate(flat) if t == ('ldr', ('x0', '[sp,#240]'))]; assert len(xarm) == 1
 assert flat[xarm[0] + 1] == ('adc', ('x0', 'x0', 'xzr')) and flat[xarm[0] + 2] == ('str', ('x0', '[sp,#240]'))
 epi = [i for i, t in enumerate(flat) if t == ('mov', ('x0', '#' + str(CBASE)))]; assert len(epi) == 1
 assert [mn for mn, _ in flat[epi[0]:epi[0] + 9]] == 'mov add ldr ldr ldr adds adc str str'.split()
 assert flat[epi[0] + 2] == ('ldr', ('x1', '[sp,#240]'))   # epilogue shape (behavior below)
 zbr = [i for i, t in enumerate(flat) if t == ('b.lo', ('.Lx_done',))]; assert len(zbr) == 1
 assert xarm[0] == zbr[0] + 1; ZPC = zbr[0] + 4; assert prog[ZPC][0] == 'subs'   # first Z-sub runs iff X >= p
 FIXPCS = init + [epi[0] + k for k in range(9)]
 vals = [0, 1, p - 1, p - 2, R]; cases = []
 for i, a in enumerate(vals):
  for j, b in enumerate(vals):
   for iv in vals: cases.append((vals[(i + j) % 5], a, vals[(2 * i + j) % 5], b or 1, iv))
 cases += [tuple(rng.randrange(1 if k == 3 else 0, p) for k in range(5)) for _ in range(cfg.get('samples', 200))]
 cases += [(0, rng.randrange(p), 0, rng.randrange(1, p), rng.randrange(p)) for _ in range(60)]
 cases += [(rng.randrange(p), rng.randrange(p), rng.randrange(p), p - rng.randrange(1, 1 << 300), rng.randrange(p)) for _ in range(60)]
 for k in range(40):                                            # operand extremes (3)
  y2 = p - 1 if k % 2 else p - 1 - rng.randrange(1 << 64); iv = p - 1 if k % 4 < 2 else rng.randrange(p - (1 << 64), p)
  mu = (p - 1 + y2) * iv * RI % p; x1 = (mu * mu * RI - (p - 1)) % p     # x' = p - 1 with x2 = 0
  cases.append((x1, p - 1, 0, y2, iv))
 G = (int(params['gx'], 16), int(params['gy'], 16)); A, B = G, add(G, G); curve = []
 for _ in range(24):
  A = add(A, G); B = add(B, add(G, G))
  if A[0] == B[0]: continue
  D = add(A, (B[0], -B[1] % p)); curve.append(len(cases))
  cases.append((A[0] * R % p, A[1] * R % p, B[0] * R % p, B[1] * R % p, pow((B[0] - A[0]) % p, -1, p) * R % p))
  assert expect(*cases[-1])[:2] == (D[0] * R % p, D[1] * R % p)
 fixed, lanes, per_call = set(), 0, []
 lo = labels['.Llane'][0]; hi = labels['.Lexit'][0]
 body_visits = {i: 0 for i in range(lo, hi)}
 for begin in range(0, len(cases), 8):
  grp = cases[begin:begin + 8]; n = len(grp); m = dict(cm); want = []
  for i, (x1, y1, x2, y2, iv) in enumerate(grp):
   out, pp, tt, ii = 0x20000 + 256 * i, 0x30000 + 256 * i, 0x40000 + 256 * i, 0x50000 + 256 * i
   for j, v in enumerate((out, pp, tt, ii)): m[0x10000 + 32 * i + 8 * j] = v
   put(m, pp, x1); put(m, pp + 48, y1); put(m, tt, x2); put(m, tt + 48, y2)
   put(m, tt + 96, p - x2); put(m, tt + 144, p - y2); put(m, ii, iv)
   want.append(expect(x1, y1, x2, y2, iv))
  m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0         # terminator args[n] = {NULL, NULL}
  m[CBASE] = m[CBASE + 8] = 0                                  # BSS starts zero
  regs = {f'x{i}': rng.getrandbits(64) for i in range(31)}; regs.update(x0=0x10000, x1=n, sp=0x7000000)
  mem = dict(m); ro = {}; vis = {i: 0 for i in list(range(lo, hi)) + FIXPCS}
  c = a64.a64_run(prog, labels, 'side_lanes', regs, mem, limit=50000, out_regs=ro, visits=vis)
  body = sum(v for i, v in vis.items() if lo <= i < hi); assert body <= WORST * n + 1, (body, n)   # +1: the final lane's exit branch b .Lexit (fixed part)
  fixed.add(c - body)
  for i, v in vis.items():
   if lo <= i < hi: body_visits[i] += v
  for pc in FIXPCS: assert vis[pc] == 1, (pc, vis[pc])
  exp = vis[ZPC]
  assert mem[regs['sp'] - FRSZ + 240] == exp and mem[CBASE] == exp and mem[CBASE + 8] == 0   # counters == expensive lanes
  for a_, v in m.items():
   if a_ in (CBASE, CBASE + 8): continue
   assert mem[a_] == v or 0x20000 <= a_ < 0x30000, hex(a_)
  for i, (xx, yy, key) in enumerate(want):
   o = 0x20000 + 256 * i
   assert get(mem, o) == xx and get(mem, o + 48) == yy and [mem[o + 96 + 8 * j] for j in range(3)] == key, (begin + i)
  for reg in ['sp'] + [f'x{i}' for i in range(18, 30)]: assert ro[reg] == regs[reg], reg
  per_call.append((n, body)); lanes += n
 assert len(fixed) == 1 and max(fixed) == FIXED, fixed
 m0 = dict(cm); m0[CBASE] = m0[CBASE + 8] = m0[0x10000] = m0[0x10008] = 0
 r0 = {f'x{i}': 0 for i in range(31)}; r0.update(x0=0x10000, x1=0, sp=0x7000000); mem0, ro0 = dict(m0), {}
 assert a64.a64_run(prog, labels, 'side_lanes', r0, mem0, limit=5000, out_regs=ro0) == 3 and mem0[CBASE] == 0 and ro0['sp'] == 0x7000000   # empty call bypass
 stride_lanes = stride_audit(prog, labels, cm, rng, 128)      # production layout passes
 try: stride_audit(prog, labels, cm, rng, 96); negative = False   # the superseded 96-byte stride must be caught
 except AssertionError: negative = True
 assert negative
 branches = 0
 for i in range(lo, hi):
  mn, ops = prog[i]
  if mn.startswith('b.') or mn in ('tbz', 'tbnz', 'cbz'):   # (cbnz = loop control; its exit is the fixed part)
   tgt = labels[ops[-1]][0]
   assert body_visits[tgt] > 0 and body_visits.get(i + 1, 1) > 0, (i, mn)   # both successors reached (i+1=.Lexit: loop exit)
   branches += 1
 return {'status': 'passed', 'lane_paths': len(paths), 'worst_path_per_lane': worst,
   'best_path_per_lane': min(n for n, _ in paths), 'lanes': lanes, 'genuine_curve_cases': len(curve),
   'fixed_per_call': max(fixed), 'fixed_s20': 32, 'fixed_s24': 10, 'fixed_charged_per_call': CHARGED, 'stride128_footprint_lanes': stride_lanes, 'stride96_negative_control_caught': True, 'conditional_branches_both_ways': branches,
   'max_lane_body_per_lane': max(b / n for n, b in per_call), 'seconds': round(time.time() - t0, 1),
   'scope': 'Word-level a64.py interpretation vs exact arithmetic; path enumeration; X-counter exact every call (s24). Callgrind: measure/measurement.log.'}

if __name__ == '__main__': print(json.dumps(run(json.load(sys.stdin)), separators=(',', ':')))
