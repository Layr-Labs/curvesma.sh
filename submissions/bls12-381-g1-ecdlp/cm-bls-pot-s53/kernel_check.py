"""Kernel audit (s1, s26, s40): path enumeration, a64.py runs against exact arithmetic, stride footprint."""
import sys, json, re, random, time
from pathlib import Path
root = Path(__file__).resolve().parent
if (root / 'files').is_dir(): root = root / 'files'
sys.path.insert(0, str(root))
import a64, rows as RW
params = json.loads((root / 'params.json').read_text()); p = int(params['p'], 16)
R = pow(2, 384, p); RI = pow(R, -1, p); MASK = (1 << 64) - 1; BASE = 0x80000; IV = 117 * p // 100
WORST = None
EOLD = ('\tldp\tx0,x29,[x30,#8]\t\t// P, T (T stays in x29 through the second product)\n\tldp\tx19,x20,[x0]\t\t// x1\n\tldp\tx21,x22,[x0,#16]\n\tldp\tx23,x24,[x0,#32]\n'
 '\tldp\tx25,x26,[x29,#96]\t// v2 = p - x2\n\tldp\tx27,x28,[x29,#112]\n\tldp\tx17,x3,[x29,#128]\n\tadds\tx11,x11,x25\n\tadcs\tx12,x12,x26\n\tadcs\tx13,x13,x27\n'
 '\tadcs\tx14,x14,x28\n\tadcs\tx15,x15,x17\n\tadc\tx16,x16,x3\n\tsubs\tx11,x11,x19\n\tsbcs\tx12,x12,x20\n\tsbcs\tx13,x13,x21\n\tsbcs\tx14,x14,x22\n\tsbcs\tx15,x15,x23\n'
 '\tsbcs\tx16,x16,x24\t\t// lo iff X < 0\n\tb.hs\t.Lx_nonneg')
ENEW = ('\tldp\tx29,x0,[x30,#16]\t// T, buf (s31: e = v2 - x1 at buf+96, 384-bit two\'s complement)\n\tldp\tx19,x20,[x0,#96]\n\tldp\tx21,x22,[x0,#112]\n\tldp\tx23,x24,[x0,#128]\n'
 '\tldp\tx25,x26,[x29,#96]\t// v2 = p - x2\n\tldp\tx27,x28,[x29,#112]\n\tldp\tx17,x3,[x29,#128]\n\tadds\tx11,x11,x19\n\tadcs\tx12,x12,x20\n\tadcs\tx13,x13,x21\n'
 '\tadcs\tx14,x14,x22\n\tadcs\tx15,x15,x23\n\tadc\tx16,x16,x24\t\t// X = V + e in (-p, 2.16p): bit 383 is the sign\n\ttbz\tx16,#63,.Lx_nonneg')
def esum_text(tx):
 """s32: X = V + e, e precomputed by the batch pass."""
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
MU0 = '\n.Llane:\n'
def mu_text(tx, plus=False):
 """s46: the batch round supplies mu; the lane start (t = y1 + y2 and the first product) becomes four loads."""
 i = tx.index(MU0) + len(MU0); j = tx.index('\tstp\tx11,x12,[sp,#96]'); o = RW.OMP if plus else RW.OMM
 assert tx.count(MU0) == 1 and tx.count('\tstp\tx11,x12,[sp,#96]') == 1 and tx[i:j].count('// ---- end MULNR') == 1
 tx = tx[:i] + '\tldr\tx2,[x30,#24]\t\t// buf (s46: mu from the batch round)\n' + ''.join('\tldp\tx%d,x%d,[x2,#%d]\n' % (11 + 2 * k, 12 + 2 * k, o + 16 * k) for k in range(3)) + tx[j:]
 for k, r in enumerate(('x19,x20', 'x21,x22', 'x23,x24')):
  a = '\tldp\t%s,[x0,#%d]' % (r, 96 + 16 * k); assert tx.count(a) == 1; tx = tx.replace(a, '\tldp\t%s,[x0,#%d]' % (r, RW.OE + 16 * k))
 return tx
def rows_text(tx):
 """s46: a full child stores rows 1..5 of y' (out-240) and of x' (out+192) before the next lane."""
 nxt = "\tldp\tx0,x1,[x30,#32]!\t// next lane's out (NULL at args[n]) and P\n"
 parts = tx.split(nxt); assert len(parts) == 5
 out = parts[0]
 for k in range(4):
  y19 = parts[k].rstrip().endswith('stp\tx23,x24,[x2,#80]'); assert y19 or parts[k].rstrip().endswith('stp\tx15,x16,[x2,#80]')
  R, V, Q, Tm = (['x19', 'x20', 'x21', 'x22', 'x23', 'x24'], ['x0', 'x1', 'x3', 'x4', 'x11', 'x12'], 'x13', 'x14') if y19 else (['x11', 'x12', 'x13', 'x14', 'x15', 'x16'], ['x0', 'x1', 'x3', 'x4', 'x19', 'x20'], 'x21', 'x22')
  L = RW.rowbuild(R, 'x2', RW.yoff, V, Q, Tm) + RW.ldrow(R, 'x2', 0) + RW.rowbuild(R, 'x2', RW.xoff, V, Q, Tm)
  out += '\t// s46 rows of y\' and x\'\n' + ''.join('\t' + l.replace(' ', '\t', 1) + '\n' for l in L) + nxt + parts[k + 1]
 return out
def base_text(): return drop_parity_text(esum_text((root / 'measure/side_lanes.S').read_text()))
def sqr_1t():
 """s53: 1-temp interleaved 384-bit squaring of x11..x16 -> (x17, x19..x28, x16); 91 instructions,
 uses only x0 as temporary, leaving x1, x2 (out), x3, x4, x5..x10 (p), x18 (TMRB), x29 (buf), x30 (args) live."""
 A = ['x11', 'x12', 'x13', 'x14', 'x15', 'x16']
 W = ['x19', 'x20', 'x21', 'x22', 'x23', 'x24', 'x25', 'x26', 'x27', 'x28']  # words 1..10
 o = ['mul %s,%s,%s' % (W[j - 1], A[j], A[0]) for j in range(1, 6)]
 for j in range(1, 6):
  op = 'adds' if j == 1 else ('adc' if j == 5 else 'adcs')
  o += ['umulh x0,%s,%s' % (A[j], A[0]), '%s %s,%s,x0' % (op, W[j], 'xzr' if j == 5 else W[j])]
 for i in range(1, 5):
  for j in range(i + 1, 6):
   op = 'adds' if j == i + 1 else 'adcs'
   o += ['mul x0,%s,%s' % (A[j], A[i]), '%s %s,%s,x0' % (op, W[i + j - 1], W[i + j - 1])]
  if i < 4:
   o += ['adc %s,xzr,xzr' % W[i + 5]]
   for j in range(i + 1, 6):
    op = 'adds' if j == i + 1 else ('adc' if j == 5 else 'adcs')
    o += ['umulh x0,%s,%s' % (A[j], A[i]), '%s %s,%s,x0' % (op, W[i + j], W[i + j])]
  else:
   o += ['umulh x0,x16,x15', 'adc x28,x0,xzr']
 o += ['adds x19,x19,x19'] + ['adcs %s,%s,%s' % (W[k], W[k], W[k]) for k in range(1, 10)]
 o += ['mul x17,x11,x11', 'umulh x0,x11,x11', 'adds x19,x19,x0']
 for i in range(1, 5):
  o += ['mul x0,%s,%s' % (A[i], A[i]), 'adcs %s,%s,x0' % (W[2 * i - 1], W[2 * i - 1]),
        'umulh x0,%s,%s' % (A[i], A[i]), 'adcs %s,%s,x0' % (W[2 * i], W[2 * i])]
 o += ['mul x0,x16,x16', 'adcs x28,x28,x0', 'umulh x16,x16,x16', 'adc x16,x16,xzr']
 assert len(o) == 91, len(o)
 return _fmt_ins(o)
def _fmt_ins(L): return ''.join('\t' + l.replace(' ', '\t', 1) + '\n' for l in L)
def tmr_leaf_block(branch=True):
 """s53: 7-chunk TMR with x1=384, x2=out, x29=buf pinned across SQR+TMR; reads e directly from [x29,#OE]."""
 o = ['and x3,x17,#0x7fffffffffffff', 'madd x0,x3,x1,x18',
  'ldp x11,x12,[x0,#-256]', 'ldp x13,x14,[x0,#-240]', 'ldp x15,x0,[x0,#-224]',
  'adds x11,x24,x11', 'adcs x12,x25,x12', 'adcs x13,x26,x13', 'adcs x14,x27,x14', 'adcs x15,x28,x15', 'adcs x16,x16,x0']
 L = ['x17', 'x19', 'x20', 'x21', 'x22', 'x23']
 for i in range(1, 7):
  off = 48 * i - 256
  o += [x.format(C='x3', **{('L%d' % k): L[k] for k in range(6)}) for x in RW.TMR_EXT[i]]
  o += ['madd x0,x3,x1,x18', 'ldp x24,x25,[x0,#%d]' % off, 'ldp x26,x27,[x0,#%d]' % (off + 16), 'ldp x28,x0,[x0,#%d]' % (off + 32),
   'adds x11,x11,x24', 'adcs x12,x12,x25', 'adcs x13,x13,x26', 'adcs x14,x14,x27', 'adcs x15,x15,x28', 'adcs x16,x16,x0']
 o += ['ldp x24,x25,[x29,#%d]' % RW.OE, 'ldp x26,x27,[x29,#%d]' % (RW.OE + 16), 'ldp x28,x0,[x29,#%d]' % (RW.OE + 32),
  'adds x11,x11,x24', 'adcs x12,x12,x25', 'adcs x13,x13,x26', 'adcs x14,x14,x27', 'adcs x15,x15,x28', 'adcs x16,x16,x0']
 o += ['lsr x3,x16,#56', 'madd x0,x3,x1,x18', 'ldp x24,x25,[x0,#80]', 'ldp x26,x27,[x0,#96]', 'ldp x28,x0,[x0,#112]',
  'subs x11,x11,x24', 'sbcs x12,x12,x25', 'sbcs x13,x13,x26', 'sbcs x14,x14,x27', 'sbcs x15,x15,x28', 'sbc x16,x16,x0']
 if branch:
  o += ['subs x19,x11,x5', 'sbcs x20,x12,x6', 'sbcs x21,x13,x7', 'sbcs x22,x14,x8', 'sbcs x23,x15,x9', 'sbcs x24,x16,x10', 'b.hs .Lx_y']
 return _fmt_ins(o)
def kernel_text(plus=False, rows=True):
 """s53: full-child kernel with 1-temp SQR keeping x2=out & x29=buf live, direct [x29,#omu] mu read (no sp spill),
 and branch canonicalization for x' and y' (worst 765, best 764, FIXED=32)."""
 raw = (root / 'measure/side_lanes.S').read_text()
 pro = raw[:raw.index('.Llane:\n')]
 epi = raw[raw.index('.Lexit:\n'):]
 pro = pro.replace('\tldr\tx1,[x30,#8]\t\t// P of lane 0 (later lanes: loaded with the next out)\n', '\tmov\tx1,#384\t\t\t// s52: pinned TMR stride\n', 1)
 old_p = '\tstp\tx5,x6,[x29,#128]\n\tstp\tx7,x8,[x29,#144]\n\tstp\tx9,x10,[x29,#160]\n'
 new_p = '\tldr\tx3,[x0,#56]\t\t// s53: TMR table base (keeps x2=out live)\n\tstp\tx18,x3,[sp,#0]\t\t// [sp,#0]=n0, [sp,#8]=TMRB\n\tmov\tx18,x3\n'
 assert pro.count(old_p) == 1 and '\tmov\tx1,#384' in pro
 pro = pro.replace(old_p, new_p, 1)
 omu = RW.OMP if plus else RW.OMM
 on2 = 48 if plus else 144
 b = ['ldr x29,[x30,#24]'] + ['ldp x%d,x%d,[x29,#%d]' % (11 + 2 * k, 12 + 2 * k, omu + 16 * k) for k in range(3)]
 body = '.Llane:\n' + _fmt_ins(b) + sqr_1t() + tmr_leaf_block(branch=False)
 W12 = ['x0', 'x2', 'x3', 'x4', 'x11', 'x12', 'x13', 'x14', 'x15', 'x16', 'x17', 'x25']
 Y6 = ['x19', 'x20', 'x21', 'x22', 'x23', 'x24']
 X6 = ['x11', 'x12', 'x13', 'x14', 'x15', 'x16']
 T6 = ['x0', 'x2', 'x3', 'x4', 'x11', 'x12']
 def _x_finish(src):
  return (RW.strow(src, 'x2', 0) + ['ldr x0,[x30,#16]', 'ldp x3,x4,[x0,#112]', 'ldp x25,x26,[x0,#128]', 'ldp x0,x2,[x0,#96]']
   + ['adds x19,%s,x0' % src[0], 'adcs x20,%s,x2' % src[1], 'adcs x21,%s,x3' % src[2], 'adcs x22,%s,x4' % src[3], 'adcs x23,%s,x25' % src[4], 'adc x24,%s,x26' % src[5]])
 c = ['subs x19,x11,x5', 'sbcs x20,x12,x6', 'sbcs x21,x13,x7', 'sbcs x22,x14,x8', 'sbcs x23,x15,x9', 'sbcs x24,x16,x10', 'b.hs .Lx_ge']
 c += _x_finish(X6) + ['b .Lx_join']
 c_ge = ['.Lx_ge:'] + _x_finish(Y6) + ['.Lx_join:']
 c2 = RW.mul6x6_1t(Y6, 'x29', omu, W12, 'x28', ['x26', 'x27'])
 c2 += RW.tmr_l(W12[:6], Y6, 'x28', 'x28', 'x27', W12[6:], 'x18', first_into=W12[6:], stride_reg='x1')
 c2 += ['ldr x29,[x30,#16]'] + RW.ldrow(T6, 'x29', on2)
 c2 += ['subs x19,x19,x0', 'sbcs x20,x20,x2', 'sbcs x21,x21,x3', 'sbcs x22,x22,x4', 'sbcs x23,x23,x11', 'sbc x24,x24,x12']
 c2 += RW.tmr_corr_l(Y6, 'x28', 'x28', 'x27', T6, 'x18', 'x1')
 c2 += ['subs x11,x19,x5', 'sbcs x12,x20,x6', 'sbcs x13,x21,x7', 'sbcs x14,x22,x8', 'sbcs x15,x23,x9', 'sbcs x16,x24,x10', 'b.hs .Ly_ge']
 def _y_tail(src, V6, Q, Tm, fall_exit=False):
  t = ['ldr x2,[x30]'] + RW.strow(src, 'x2', 48)
  if rows:
   t += ['ldr x18,[sp,#0]'] + RW.rowbuild(src, 'x2', RW.yoff, V6, Q, Tm) + RW.ldrow(src, 'x2', 0) + RW.rowbuild(src, 'x2', RW.xoff, V6, Q, Tm) + ['ldr x18,[sp,#8]']
  t += ['ldr x2,[x30,#32]!', 'cbnz x2,.Llane']
  if not fall_exit: t += ['b .Lexit']
  return t
 c2 += _y_tail(Y6, ['x0', 'x3', 'x4', 'x11', 'x12', 'x13'], 'x14', 'x15', False)
 c3 = _y_tail(X6, ['x0', 'x3', 'x4', 'x19', 'x20', 'x23'], 'x21', 'x22', True)
 return pro + body + _fmt_ins(c) + '\n'.join(l if l.endswith(':') else '\t' + l.replace(' ', '\t', 1) for l in c_ge) + '\n' + _fmt_ins(c2) + '.Ly_ge:\n' + _fmt_ins(c3) + epi
def put_e(m, ii, x1, x2): put(m, ii + RW.OE, (p - x2 - x1) % (1 << 384))
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

def prepare(text, base=None):
 code, data = text.split('.section .rodata', 1)
 vals = []
 for line in data.splitlines():
  if '.quad' in line: vals += [int(x.strip(), 0) for x in line.split('.quad', 1)[1].split(',')]
 assert sum(vals[i] << (64 * i) for i in range(6)) == p and vals[6] == (-pow(p, -1, 1 << 64)) & MASK
 if base is not None: vals = vals[:7] + [base]
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

def mu_of(y1, y2, iv, rng=None):
 """a congruent mu in the batch round's output range [0, p + 2^376)."""
 mu = (y1 + y2) * iv * RI % p
 return mu + p if rng and mu < (1 << 376) and rng.random() < 0.5 else mu
def expect_mu(x1, y1, x2, y2, mu):
 xx = (mu * mu * RI - x1 - x2) % p; yy = (mu * ((xx - x2) % p) * RI + y2) % p
 return xx, yy
def expect(x1, y1, x2, y2, iv):
 mu = (y1 + y2) * iv * RI % p; xx = (mu * mu * RI - x1 - x2) % p
 yy = (mu * ((xx - x2) % p) * RI + y2) % p
 return xx, yy

SLOT = RW.XR + 240 - RW.YR
def lane_mem(m, i, x1, y1, x2, y2, iv, plus, rng, out):
 pp, tt, ii = 0x30000 + 256 * i, 0x40000 + 256 * i, 0x50000 + 512 * i
 for j, v in enumerate((out, pp, tt, ii)): m[0x10000 + 32 * i + 8 * j] = v
 put(m, pp, x1); put(m, pp + 48, y1); put(m, tt, x2); put(m, tt + 48, p - y2 if plus else y2)
 put(m, tt + 96, p - x2); put(m, tt + 144, y2 if plus else p - y2); put_e(m, ii, x1, x2)
 mu = mu_of(y1, y2, iv, rng); put(m, ii + (RW.OMP if plus else RW.OMM), mu)
 return expect_mu(x1, y1, x2, y2, mu)
def rows_ok(mem, o, xx, yy):
 return all(get(mem, o + RW.yoff(j)) == RW.rows_of(yy)[j] and get(mem, o + RW.xoff(j)) == RW.rows_of(xx)[j] for j in range(6))
def stride_audit(prog, labels, cm, rng, stride, n=8, sp=0x7000000, plus=False):
 """s8/s46 layout: full-child slots [out-240, out+432) under canaries, keys 96..191 untouched; an overlapping stride must fail."""
 base, m, want = 0x20000 + 256, dict(cm), []
 for i in range(n):
  x1, y1, x2, y2, iv = (rng.randrange(p) for _ in range(5)); y2 = y2 or 1
  want.append(lane_mem(m, i, x1, y1, x2, y2, iv, plus, rng, base + stride * i))
  mu = get(m, 0x50000 + 512 * i + (RW.OMP if plus else RW.OMM))
  e = get(m, 0x50000 + 512 * i + RW.OE); e = e - (1 << 384) if e >> 383 else e
  RW.tmr_put(m, TB, mu * mu, e, put)
  xx_k = RW.tmr_val(mu * mu, e) % p
  RW.tmr_put(m, TB, (xx_k + p - x2) * mu, -(p - y2), put)
 for a_ in range(base + RW.YR, base + stride * n + RW.XR + 240 + 128, 8): m[a_] = rng.getrandbits(64)
 m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0
 regs = {f'x{i}': rng.getrandbits(64) for i in range(31)}; regs.update(x0=0x10000, x1=n, sp=sp)
 mem = dict(m); ro = {}
 a64.a64_run(prog, labels, 'side_lanes', regs, mem, limit=50000, out_regs=ro)
 outs = set()
 for i in range(n):
  o = base + stride * i; outs.update(range(o + RW.YR, o + 96, 8)); outs.update(range(o + RW.XR, o + RW.XR + 240, 8))
 written = {a_ for a_ in mem if a_ not in m or mem[a_] != m[a_]}
 frame = {a_ for a_ in written if sp - 1024 <= a_ < sp}
 assert written - frame <= outs, sorted(hex(a_) for a_ in (written - frame) - outs)
 for i, (xx, yy) in enumerate(want):
  o = base + stride * i
  assert get(mem, o) == xx and get(mem, o + 48) == yy and rows_ok(mem, o, xx, yy), i
 assert ro['sp'] == sp
 return n

def run(cfg):
 t0 = time.time(); rng = random.Random(cfg.get('seed', 1)); plus = cfg.get('plus', False)
 text = kernel_text(plus)
 paths = lane_paths(text); worst = max(n for n, _ in paths)
 assert worst == 765, worst
 prog, labels, cm = prepare(text, TB)
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
   want.append(lane_mem(m, i, x1, y1, x2, y2, iv, plus, rng, 0x20000 + 256 + 1024 * i))
   mu = get(m, 0x50000 + 512 * i + (RW.OMP if plus else RW.OMM))
   e = get(m, 0x50000 + 512 * i + RW.OE); e = e - (1 << 384) if e >> 383 else e
   RW.tmr_put(m, TB, mu * mu, e, put)
   xx_k = RW.tmr_val(mu * mu, e) % p
   RW.tmr_put(m, TB, (xx_k + p - x2) * mu, -(p - y2), put)
  m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0
  regs = {f'x{i}': rng.getrandbits(64) for i in range(31)}; regs.update(x0=0x10000, x1=n, sp=0x7000000)
  mem = dict(m); ro = {}; vis = {i: 0 for i in range(lo, hi)}
  c = a64.a64_run(prog, labels, 'side_lanes', regs, mem, limit=50000, out_regs=ro, visits=vis)
  body = sum(vis.values()); assert body <= worst * n + 1, (body, n)
  fixed.add(c - body)
  for i, v in vis.items(): body_visits[i] += v
  for a_, v in m.items(): assert mem[a_] == v or 0x20000 <= a_ < 0x30000, hex(a_)
  for i, (xx, yy) in enumerate(want):
   o = 0x20000 + 256 + 1024 * i
   assert get(mem, o) == xx and get(mem, o + 48) == yy and rows_ok(mem, o, xx, yy), (begin + i)
  for reg in ['sp'] + [f'x{i}' for i in range(18, 30)]: assert ro[reg] == regs[reg], reg
  per_call.append((n, body)); lanes += n
 assert len(fixed) == 1 and max(fixed) == FIXED, fixed
 stride_lanes = stride_audit(prog, labels, cm, rng, SLOT, plus=plus)
 try: stride_audit(prog, labels, cm, rng, SLOT - 32, plus=plus); negative = False
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
 'fixed_per_call': max(fixed), 'fixed_charged_per_call': CHARGED, 'stride672_footprint_lanes': stride_lanes, 'stride640_negative_control_caught': True, 'rows_checked': True, 'conditional_branches_both_ways': branches,
 'max_lane_body_per_lane': max(b / n for n, b in per_call), 'seconds': round(time.time() - t0, 1),
 'scope': 'Word-level interpretation (a64.py) of the s53 TMR full-child listing against exact arithmetic; '
 'static path enumeration of the lane body. Callgrind record: measure/measurement.log.'}

def leaf_text(tx=None, plus=False):
 """s53 leaf kernel: 1-temp SQR keeps x2=out & x29=buf live across SQR+TMR, x1=384, x18=TMRB (209 both paths)."""
 raw = (root / 'measure/side_lanes.S').read_text()
 pro = raw[:raw.index('.Llane:\n')]
 epi = raw[raw.index('.Lexit:\n'):]
 pro = pro.replace('\tldr\tx1,[x30,#8]\t\t// P of lane 0 (later lanes: loaded with the next out)\n', '\tmov\tx1,#384\t\t\t// s52: pinned TMR stride\n', 1)
 a = '\tldr\tx18,[x0,#48]\t\t// n0, pinned in x18'; assert pro.count(a) == 1
 pro = pro.replace(a, '\tldr\tx18,[x0,#56]\t\t// s52: TMR table base (n0 unused by leaf)', 1)
 omu = RW.OMP if plus else RW.OMM
 b = ['ldr x29,[x30,#24]'] + ['ldp x%d,x%d,[x29,#%d]' % (11 + 2 * k, 12 + 2 * k, omu + 16 * k) for k in range(3)]
 done = ['.Lx_done:', '\tstp\tx11,x12,[x2]', '\tstp\tx13,x14,[x2,#16]', '\tstp\tx15,x16,[x2,#32]',
  '\tldr\tx2,[x30,#32]!', '\tcbnz\tx2,.Llane', '\tb\t.Lexit',
  '.Lx_y:', '\tstp\tx19,x20,[x2]', '\tstp\tx21,x22,[x2,#16]', '\tstp\tx23,x24,[x2,#32]',
  '\tldr\tx2,[x30,#32]!', '\tcbnz\tx2,.Llane']
 return pro + '.Llane:\n' + _fmt_ins(b) + sqr_1t() + tmr_leaf_block(branch=True) + '\n'.join(done) + '\n' + epi

def leaf_worst(): return max(n for n, _ in lane_paths(leaf_text()))
TB = 1 << 40
def run_leaf(cfg):
 """s30/s47/s52/s53: leaf kernel executed against exact x' (table entries the lanes read written from the exact model)."""
 rng = random.Random(cfg.get('seed', 1)); plus = cfg.get('plus', False)
 text = leaf_text(plus=plus)
 paths = lane_paths(text)
 worst, best = max(n for n, _ in paths), min(n for n, _ in paths); assert (worst, best) == (209, 209), (worst, best)
 prog, labels, cm = prepare(text, TB); lo, hi = labels['.Llane'][0], labels['.Lexit'][0]; vis = {i: 0 for i in range(lo, hi)}
 v5 = [0, 1, p - 1, p - 2, R]
 cs = [(v5[(i + j) % 5], a, v5[(2 * i + j) % 5], b or 1, iv) for i, a in enumerate(v5) for j, b in enumerate(v5) for iv in v5]
 cs += [tuple(rng.randrange(1 if k == 3 else 0, IV if k == 4 else p) for k in range(5)) for _ in range(cfg.get('samples', 200))]
 for k in range(40):
  y2 = p - 1 - (k % 2) * rng.randrange(1 << 64); iv = IV - 1 if k % 2 else p - 1
  mu = (p - 1 + y2) * iv * RI % p; cs.append(((mu * mu * RI - (p - 1)) % p, p - 1, 0, y2, iv))
 for b0 in range(0, len(cs), 8):
  g = cs[b0:b0 + 8]; n = len(g); m = dict(cm)
  mus = []
  for i, (x1, y1, x2, y2, iv) in enumerate(g):
   lane_mem(m, i, x1, y1, x2, y2, iv, plus, rng, 0x20000 + 256 * i)
   mu = get(m, 0x50000 + 512 * i + (RW.OMP if plus else RW.OMM)); mus.append(mu)
   e = get(m, 0x50000 + 512 * i + RW.OE); e = e - (1 << 384) if e >> 383 else e
   RW.tmr_put(m, TB, mu * mu, e, put)
  m[0x10000 + 32 * n] = m[0x10000 + 32 * n + 8] = 0
  rg = {f'x{i}': rng.getrandbits(64) for i in range(31)}; rg.update(x0=0x10000, x1=n, sp=0x7000000); mem, ro, vv = dict(m), {}, {i: 0 for i in range(lo, hi)}
  c = a64.a64_run(prog, labels, 'side_lanes', rg, mem, limit=50000, out_regs=ro, visits=vv); body = sum(vv.values())
  assert 209 * n <= body <= 209 * n + 1 and c - body == FIXED, (body, n, c - body)
  out = {0x20000 + 256 * i + 8 * k for i in range(n) for k in range(6)}
  assert all(mem[a_] == v or a_ in out for a_, v in m.items()) and {a_ for a_ in mem if a_ not in m} - out <= set(range(0x7000000 - 1024, 0x7000000, 8))
  assert all(get(mem, 0x20000 + 256 * i) == expect_mu(e[0], e[1], e[2], e[3], mu)[0] for i, (e, mu) in enumerate(zip(g, mus))) and all(ro[r_] == rg[r_] for r_ in ['sp'] + [f'x{i}' for i in range(19, 30)])
  for i in vv: vis[i] += vv[i]
 for i in range(lo, hi):
  mn, ops = prog[i]
  if mn.startswith('b.'): assert vis[labels[ops[-1]][0]] > 0 and vis[i + 1] > 0, i
 return {'status': 'passed', 'plus': plus, 'lane_paths': len(paths), 'worst_path_per_lane': worst, 'best_path_per_lane': best,
 'lanes': len(cs), 'fixed_per_call': FIXED, 'scope': "s53 leaf kernel with 1-temp SQR, live out/buf, pinned stride and TMR: x' exact and canonical, footprint [out, out+48)"}

if __name__ == '__main__': print(json.dumps(run(json.load(sys.stdin)), separators=(',', ':')))
