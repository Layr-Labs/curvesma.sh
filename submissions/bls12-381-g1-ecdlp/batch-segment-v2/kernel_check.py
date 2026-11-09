"""Audit of the s31 kernel esum_text(measure/side_lanes.S) (s1; plus: s26 variant): path enumeration (worst charged); a64.py word-level runs vs
exact arithmetic (fixed part, regs preserved); stride-192 footprint, 80 negative control; unused keys removed."""
import sys, json, re, random, time
from pathlib import Path
root = Path(__file__).resolve().parent
if (root / 'files').is_dir(): root = root / 'files'
sys.path.insert(0, str(root))
import a64
params = json.loads((root / 'params.json').read_text()); p = int(params['p'], 16)
R = pow(2, 384, p); RI = pow(R, -1, p); MASK = (1 << 64) - 1; BASE = 0x80000; IV = 117 * p // 100
WORST = 916
EOLD = ('\tldp\tx0,x29,[x30,#8]\t\t// P, T (T stays in x29 through the second product)\n\tldp\tx19,x20,[x0]\t\t// x1\n\tldp\tx21,x22,[x0,#16]\n\tldp\tx23,x24,[x0,#32]\n'
 '\tldp\tx25,x26,[x29,#96]\t// v2 = p - x2\n\tldp\tx27,x28,[x29,#112]\n\tldp\tx17,x3,[x29,#128]\n\tadds\tx11,x11,x25\n\tadcs\tx12,x12,x26\n\tadcs\tx13,x13,x27\n'
 '\tadcs\tx14,x14,x28\n\tadcs\tx15,x15,x17\n\tadc\tx16,x16,x3\n\tsubs\tx11,x11,x19\n\tsbcs\tx12,x12,x20\n\tsbcs\tx13,x13,x21\n\tsbcs\tx14,x14,x22\n\tsbcs\tx15,x15,x23\n'
 '\tsbcs\tx16,x16,x24\t\t// lo iff X < 0\n\tb.hs\t.Lx_nonneg')
ENEW = ('\tldp\tx29,x0,[x30,#16]\t// T, buf (s31: e = v2 - x1 at buf+96, 384-bit two\'s complement)\n\tldp\tx19,x20,[x0,#96]\n\tldp\tx21,x22,[x0,#112]\n\tldp\tx23,x24,[x0,#128]\n'
 '\tldp\tx25,x26,[x29,#96]\t// v2 = p - x2\n\tldp\tx27,x28,[x29,#112]\n\tldp\tx17,x3,[x29,#128]\n\tadds\tx11,x11,x19\n\tadcs\tx12,x12,x20\n\tadcs\tx13,x13,x21\n'
 '\tadcs\tx14,x14,x22\n\tadcs\tx15,x15,x23\n\tadc\tx16,x16,x24\t\t// X = V + e in (-p, 2.16p): bit 383 is the sign\n\ttbz\tx16,#63,.Lx_nonneg')
def esum_text(tx):
 """s31: X = V + e with e = v2 - x1 precomputed once per lane op by the batch forward pass (buf+96..143); -6 per lane."""
 assert tx.count(EOLD) == 1; return tx.replace(EOLD, ENEW)
def drop_parity_text(tx):
 """Remove four unused output-key blocks; coordinates are already stored."""
 lines=tx.splitlines(True); out=[]; i=0; blocks=0
 while i<len(lines):
  if re.match(r'\ttbz\t(?:x19|x11),#0,\.Lk_[ed]2?',lines[i]):
   j=i+1
   while not lines[j].startswith('\tldp\tx0,x1,[x30,#32]!'): j+=1
   ins=[v.split('//')[0].strip() for v in lines[i:j] if v.startswith('\t')]
   assert [v.split()[0] for v in ins]==['tbz','subs','sbcs','sbcs','sbc','stp','str']
   assert ins[-2].endswith('[x2,#96]') and ins[-1].endswith('[x2,#112]')
   blocks+=1; i=j
  else: out.append(lines[i]); i+=1
 assert blocks==4
 return ''.join(out)
def kernel_text(plus=False):
 tx = drop_parity_text(esum_text((root / 'measure/side_lanes.S').read_text()))
 if plus:
  from batch_audit import plus_text; tx = plus_text(tx)
 return tx
def put_e(m, ii, x1, x2): put(m, ii + 96, (p - x2 - x1) % (1 << 384))
FIXED = 32
CHARGED = 40
assert FIXED + 1 + 5 + 2 == CHARGED

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

def expect(x1, y1, x2, y2, iv):
 mu = (y1 + y2) * iv * RI % p; xx = (mu * mu * RI - x1 - x2) % p
 yy = (mu * ((xx - x2) % p) * RI + y2) % p
 return xx, yy

def stride_audit(prog, labels, cm, rng, stride, n=8, sp=0x7000000, plus=False):
 """s8 layout: n out slots at the stride under canaries; writes outside the frame lie in the slots, all else intact,
    outputs exact. Stride 80 must fail (coordinates cross a slot boundary)."""
 base, m, want = 0x20000, dict(cm), []
 for i in range(n):
  x1, y1, x2, y2, iv = (rng.randrange(p) for _ in range(5)); y2 = y2 or 1
  out, pp, tt, ii = base + stride * i, 0x30000 + 256 * i, 0x40000 + 256 * i, 0x50000 + 256 * i
  for j, v in enumerate((out, pp, tt, ii)): m[0x10000 + 32 * i + 8 * j] = v
  put(m, pp, x1); put(m, pp + 48, y1); put(m, tt, x2); put(m, tt + 48, y2)
  put(m, tt + 96, p - x2); put(m, tt + 144, p - y2); put(m, ii, iv); put_e(m, ii, x1, x2)
  if plus: put(m, tt + 48, p - y2); put(m, tt + 144, y2)
  want.append(expect(x1, y1, x2, y2, iv))
 for a_ in range(base, base + stride * n + 128, 8): m[a_] = rng.getrandbits(64)
 m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0
 regs = {f'x{i}': rng.getrandbits(64) for i in range(31)}; regs.update(x0=0x10000, x1=n, sp=sp)
 mem = dict(m); ro = {}
 a64.a64_run(prog, labels, 'side_lanes', regs, mem, limit=50000, out_regs=ro)
 slots, outs = set(), set()
 for i in range(n):
  slots.update(range(base + stride * i, base + stride * (i + 1), 8)); outs.update(range(base + stride * i, base + stride * i + 96, 8))
 written = {a_ for a_ in mem if a_ not in m or mem[a_] != m[a_]}
 frame = {a_ for a_ in written if sp - 1024 <= a_ < sp}
 assert written - frame <= slots, sorted(hex(a_) for a_ in (written - frame) - slots)
 for a_, v in m.items(): assert mem[a_] == v or a_ in outs, hex(a_)
 for i, (xx, yy) in enumerate(want):
  o = base + stride * i
  assert get(mem, o) == xx and get(mem, o + 48) == yy, i
 assert ro['sp'] == sp
 return n

def run(cfg):
 t0 = time.time(); rng = random.Random(cfg.get('seed', 1)); plus = cfg.get('plus', False)
 text = kernel_text(plus)
 paths = lane_paths(text); worst = max(n for n, _ in paths)
 assert worst == WORST, worst
 prog, labels, cm = prepare(text)
 vals = [0, 1, p - 1, p - 2, R]; cases = []
 for i, a in enumerate(vals):
  for j, b in enumerate(vals):
   for iv in vals: cases.append((vals[(i + j) % 5], a, vals[(2 * i + j) % 5], b or 1, iv))
 cases += [tuple(rng.randrange(1 if k == 3 else 0, IV if k == 4 else p) for k in range(5)) for _ in range(cfg.get('samples', 200))]
 cases += [(0, rng.randrange(p), 0, rng.randrange(1, p), rng.randrange(p)) for _ in range(60)]
 cases += [(rng.randrange(p), rng.randrange(p), rng.randrange(p), p - rng.randrange(1, 1 << 300), rng.randrange(p)) for _ in range(60)]
 for k in range(40):
  y2 = p - 1 if k % 2 else p - 1 - rng.randrange(1 << 64); iv = (p - 1, IV - 1)[k % 2] if k % 4 < 2 else rng.randrange(IV - (1 << 64), IV)
  mu = (p - 1 + y2) * iv * RI % p; x1 = (mu * mu * RI - (p - 1)) % p
  cases.append((x1, p - 1, 0, y2, iv))
 G = (int(params['gx'], 16), int(params['gy'], 16)); A, B = G, add(G, G); curve = []
 for _ in range(24):
  A = add(A, G); B = add(B, add(G, G))
  if A[0] == B[0]: continue
  D = add(A, (B[0], -B[1] % p)); curve.append(len(cases))
  cases.append((A[0] * R % p, A[1] * R % p, B[0] * R % p, B[1] * R % p, pow((B[0] - A[0]) % p, -1, p) * R % p))
  assert expect(*cases[-1])[:2] == (D[0] * R % p, D[1] * R % p)
 unreduced = sum(c[4] >= p for c in cases); assert unreduced >= 40
 fixed, lanes, per_call = set(), 0, []
 lo = labels['.Llane'][0]; hi = labels['.Lexit'][0]
 body_visits = {i: 0 for i in range(lo, hi)}
 for begin in range(0, len(cases), 8):
  grp = cases[begin:begin + 8]; n = len(grp); m = dict(cm); want = []
  for i, (x1, y1, x2, y2, iv) in enumerate(grp):
   out, pp, tt, ii = 0x20000 + 256 * i, 0x30000 + 256 * i, 0x40000 + 256 * i, 0x50000 + 256 * i
   for j, v in enumerate((out, pp, tt, ii)): m[0x10000 + 32 * i + 8 * j] = v
   put(m, pp, x1); put(m, pp + 48, y1); put(m, tt, x2); put(m, tt + 48, y2)
   put(m, tt + 96, p - x2); put(m, tt + 144, p - y2); put(m, ii, iv); put_e(m, ii, x1, x2)
   if plus: put(m, tt + 48, p - y2); put(m, tt + 144, y2)
   want.append(expect(x1, y1, x2, y2, iv))
  m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0
  regs = {f'x{i}': rng.getrandbits(64) for i in range(31)}; regs.update(x0=0x10000, x1=n, sp=0x7000000)
  mem = dict(m); ro = {}; vis = {i: 0 for i in range(lo, hi)}
  c = a64.a64_run(prog, labels, 'side_lanes', regs, mem, limit=50000, out_regs=ro, visits=vis)
  body = sum(vis.values()); assert body <= WORST * n + 1, (body, n)
  fixed.add(c - body)
  for i, v in vis.items(): body_visits[i] += v
  for a_, v in m.items(): assert mem[a_] == v or 0x20000 <= a_ < 0x30000, hex(a_)
  for i, (xx, yy) in enumerate(want):
   o = 0x20000 + 256 * i
   assert get(mem, o) == xx and get(mem, o + 48) == yy, (begin + i)
  for reg in ['sp'] + [f'x{i}' for i in range(18, 30)]: assert ro[reg] == regs[reg], reg
  per_call.append((n, body)); lanes += n
 assert len(fixed) == 1 and max(fixed) == FIXED, fixed
 stride_lanes = stride_audit(prog, labels, cm, rng, 192, plus=plus)
 try: stride_audit(prog, labels, cm, rng, 80, plus=plus); negative = False
 except AssertionError: negative = True
 assert negative
 branches = 0
 for i in range(lo, hi):
  mn, ops = prog[i]
  if mn.startswith('b.') or mn in ('tbz', 'tbnz', 'cbz'):
   tgt = labels[ops[-1]][0]
   assert body_visits[tgt] > 0 and body_visits.get(i + 1, 1) > 0, (i, mn)
   branches += 1
 return {'status': 'passed', 'lane_paths': len(paths), 'worst_path_per_lane': worst,
   'best_path_per_lane': min(n for n, _ in paths), 'lanes': lanes, 'genuine_curve_cases': len(curve), 'unreduced_inverse_lanes': unreduced,
   'fixed_per_call': max(fixed), 'fixed_charged_per_call': CHARGED, 'stride192_footprint_lanes': stride_lanes, 'stride80_negative_control_caught': True, 'conditional_branches_both_ways': branches,
   'max_lane_body_per_lane': max(b / n for n, b in per_call), 'seconds': round(time.time() - t0, 1),
   'scope': 'Word-level interpretation (a64.py) of the derived coordinate-only listing against exact arithmetic; '
     'static path enumeration of the lane body. Callgrind record: measure/measurement.log.'}

def leaf_text(tx, plus=False):
 """s30 leaf kernel: the shipped listing cut after the x' store (y' tail and mu spill deleted); plus: y2 -> n2 loads."""
 L = tx.split('\n'); i, j, e = L.index('.Lx_done:'), L.index('.Lx_y:'), L.index('.Lexit:')
 nx = ['\tldp\tx0,x1,[x30,#32]!', '\tcbnz\tx0,.Llane']
 assert L[i + 1].startswith('\tldr\tx2,[x30]') and L[j + 1].startswith('\tldr\tx2,[x30]') and all(L[k].startswith('\tstp') for k in (i + 2, i + 3, i + 4, j + 2, j + 3, j + 4))
 t = '\n'.join(L[:i + 5] + nx + ['\tb\t.Lexit'] + L[j:j + 5] + nx + L[e:])
 vl = '\tldp\tx25,x26,[x29,#96]\t// v2 = p - x2\n\tldp\tx27,x28,[x29,#112]\n\tldp\tx17,x3,[x29,#128]\n'
 assert t.count(vl) == 1; t = t.replace(vl, '')
 sp = '\tstp\tx11,x12,[sp,#96]\n\tstp\tx13,x14,[sp,#112]\n\tstp\tx15,x16,[sp,#128]\n'
 assert t.count(sp) == 1; t = t.replace(sp, '')
 for o in (48, 64, 80) if plus else ():
  a = '[x3,#%d]' % o; assert t.count(a) == 1; t = t.replace(a, '[x3,#%d]' % (o + 96))
 for a, b in (('b.lo\t.Lx_done\t\t// X < p\n', 'ldr\tx0,[sp,#96]\n\tadd\tx0,x0,#1\n\tstr\tx0,[sp,#96]\n'), ('sub\tsp,sp,#144\t\t// [sp,#0..95] 768-bit square, [sp,#96..143] mu\n', 'str\txzr,[sp,#96]\n')):
  assert t.count(a) == 1; t = t.replace(a, a + '\t' + b)
 for a, b in (('\tadd\tsp,sp,#144\n', '\tldr\tx0,[sp,#96]\n'), ('.Lnone:\n', '\tret\n')):
  assert t.count(a) == 1; t = t.replace(a, b + a)
 t = t.replace('.Lnone:\n', '.Lnone:\n\tmov\tx0,xzr\n')
 return t

def leaf_worst(): return max(n for n, _ in lane_paths(leaf_text(kernel_text())))
def xgep(x1, y1, x2, y2, iv):
 S = lambda a, b: (a * b + (-a * b * pow(p, -1, 1 << 384) % (1 << 384)) * p) >> 384
 mu = S(y1 + y2, iv); return S(mu, mu) + p - x2 - x1 >= p

def run_leaf(cfg):
 """s30: leaf kernel (plus: P+T variant) executed vs exact x'; writes exactly [out, out+48); regs kept; worst 596."""
 rng = random.Random(cfg.get('seed', 1)); plus = cfg.get('plus', False)
 text = leaf_text(kernel_text(), plus)
 assert text == leaf_text(esum_text((root/'measure/side_lanes.S').read_text()), plus)
 paths = lane_paths(text)
 worst = max(n for n, _ in paths); assert (worst, min(n for n, _ in paths)) == (590, 580); tc = 0
 prog, labels, cm = prepare(text); lo, hi = labels['.Llane'][0], labels['.Lexit'][0]; vis = {i: 0 for i in range(lo, hi)}
 v5 = [0, 1, p - 1, p - 2, R]
 cs = [(v5[(i + j) % 5], a, v5[(2 * i + j) % 5], b or 1, iv) for i, a in enumerate(v5) for j, b in enumerate(v5) for iv in v5]
 cs += [tuple(rng.randrange(1 if k == 3 else 0, IV if k == 4 else p) for k in range(5)) for _ in range(cfg.get('samples', 200))]
 for k in range(40):
  y2 = p - 1 - (k % 2) * rng.randrange(1 << 64); iv = IV - 1 if k % 2 else p - 1
  mu = (p - 1 + y2) * iv * RI % p; cs.append(((mu * mu * RI - (p - 1)) % p, p - 1, 0, y2, iv))
 for b0 in range(0, len(cs), 8):
  g = cs[b0:b0 + 8]; n = len(g); m = dict(cm)
  for i, (x1, y1, x2, y2, iv) in enumerate(g):
   for j, v in enumerate((0x20000 + 256 * i, 0x30000 + 256 * i, 0x40000 + 256 * i, 0x50000 + 256 * i)): m[0x10000 + 32 * i + 8 * j] = v
   for a_, v in ((0x30000, x1), (0x30030, y1), (0x40000, x2), (0x40030, p - y2 if plus else y2), (0x40060, p - x2), (0x40090, y2 if plus else p - y2), (0x50000, iv)): put(m, a_ + 256 * i, v)
   put_e(m, 0x50000 + 256 * i, x1, x2)
  m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0
  rg = {f'x{i}': rng.getrandbits(64) for i in range(31)}; rg.update(x0=0x10000, x1=n, sp=0x7000000); mem, ro, vv = dict(m), {}, {i: 0 for i in range(lo, hi)}
  c = a64.a64_run(prog, labels, 'side_lanes', rg, mem, limit=50000, out_regs=ro, visits=vv); body = sum(vv.values())
  cnt = sum(xgep(*e) for e in g); tc += cnt
  assert ro['x0'] == cnt and body - 580 * n - 10 * cnt in (0, 1) and c - body == FIXED + 2
  out = {0x20000 + 256 * i + 8 * k for i in range(n) for k in range(6)}
  assert all(mem[a_] == v or a_ in out for a_, v in m.items()) and {a_ for a_ in mem if a_ not in m} - out <= set(range(0x7000000 - 1024, 0x7000000, 8))
  assert all(get(mem, 0x20000 + 256 * i) == expect(*e)[0] for i, e in enumerate(g)) and all(ro[r_] == rg[r_] for r_ in ['sp'] + [f'x{i}' for i in range(18, 30)])
  for i in vv: vis[i] += vv[i]
 for i in range(lo, hi):
  mn, ops = prog[i]
  if mn.startswith('b.'): assert vis[labels[ops[-1]][0]] > 0 and vis[i + 1] > 0, i
 return {'status': 'passed', 'plus': plus, 'lane_paths': len(paths), 'worst_path_per_lane': worst, 'best_path_per_lane': min(n for n, _ in paths),
   'lanes': len(cs), 'unreduced_inverse_lanes': sum(e[4] >= p for e in cs), 'fixed_per_call': FIXED + 2, 'x_ge_p_lanes': tc, 'scope': "s30/s37 leaf kernel: x' exact, footprint [out, out+48), X >= p count returned"}

if __name__ == '__main__': print(json.dumps(run(json.load(sys.stdin)), separators=(',', ':')))
