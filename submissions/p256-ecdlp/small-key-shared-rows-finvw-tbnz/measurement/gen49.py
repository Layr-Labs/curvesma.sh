"""Shared rows: derive bloop_branchy.h and floop.h from the y3-fold revision's (proof.md, "Shared rows").
usage: gen49.py <y3-fold measurement dir> <output dir>; every edit must match once."""
import re, sys
C, M, POOL, RR = 'x29', 'x28', ['x%d' % i for i in range(19)] + ['x30'], ['x19', 'x20', 'x21', 'x22']
class P:
 def __init__(s):
  s.ins, s.n = [], 0
 def v(s):
  s.n += 1
  return '%%v%d' % s.n
 def __call__(s, op, *a, sh=None, imm=None, d=None):
  d = d or s.v()
  s.ins.append((op, [d], list(a), sh, imm))
  return d
 def m(s, op, mem, ds=(), src=()):
  b = re.match(r'\[(%v\d+)', mem)
  s.ins.append((op, list(ds), list(src) + ([b.group(1)] if b else []), None, mem))
  return list(ds)
 def ldp(s, mem):
  return s.m('ldp', mem, (s.v(), s.v()))
 def raw(s, t):
  s.ins.append(('raw', [], [], None, t))
def row(p, X):
 t = p('lsr', X[0], imm='#32')
 b0, b1 = p('adds', X[1], X[0], sh='lsl'), p('adcs', X[2], t)
 lo, hi = p('mul', X[0], C), p('umulh', X[0], C)
 return [b0, b1, p('adcs', X[3], lo), p('adc', hi, 'xzr')]
def product(p, R, B1, B2, a):
 z = {}
 for i, Y in ((0, B2), (1, B1), (2, R)):
  s = a(i)
  if not z:
   z[0] = p('mul', s, Y[0])
   for k in (1, 2, 3):
    z[k] = p('adcs' if k > 1 else 'adds', p('mul', s, Y[k]), p('umulh', s, Y[k - 1]))
   z[4] = p('adc', p('umulh', s, Y[3]), 'xzr')
   continue
  top = p('umulh', s, Y[3])
  for k in range(4):
   z[k] = p('adcs' if k else 'adds', z[k], p('mul', s, Y[k]))
  top = p('adc', top, 'xzr')
  for k in (1, 2, 3):
   z[k] = p('adcs' if k > 1 else 'adds', z[k], p('umulh', s, Y[k - 1]))
  z[4] = p('adcs', z[4], top)
  z[5] = p('adc', z[5] if 5 in z else 'xzr', 'xzr')
 m0, s = z[0], a(3)
 lo, hi = p('mul', m0, C), p('umulh', m0, C)
 top = p('umulh', s, R[3])
 for k in range(4):
  z[k + 1] = p('adcs' if k else 'adds', z[k + 1], p('mul', s, R[k]))
 z[5] = p('adc', z[5], top)
 z[1] = p('adds', z[1], m0, sh='lsl')
 for k, w in ((2, p('extr', z[1], m0, imm='#32')), (3, lo), (4, hi)):
  z[k] = p('adcs', z[k], w)
 z[5] = p('adc', z[5], 'xzr')
 for k in (2, 3, 4):
  z[k] = p('adcs' if k > 2 else 'adds', z[k], p('umulh', s, R[k - 2]))
 z[5] = p('adc', z[5], 'xzr')
 lo, hi = p('mul', z[1], C), p('umulh', z[1], C)
 z[3] = p('adds', z[3], z[1], sh='lsr')
 z[4], z[5] = p('adcs', z[4], lo), p('adcs', z[5], hi)
 return [z[2], z[3], z[4], z[5]]
def fix(p, w, ops, label, cond):
 p.raw('b.%s %s' % (cond, label))
 for k, (op, x) in enumerate(ops):
  p(op, w[k], *([x] if x[0] != '#' else []), imm=x if x[0] == '#' else None, d=w[k])
 p.raw('add x27,x27,#5')
 p.raw(label + ':')
ADDP = [('subs', '#1'), ('adcs', M), ('adcs', 'xzr'), ('adc', C)]
SUBP = [('adds', '#1'), ('sbcs', M), ('sbcs', 'xzr'), ('sbc', C)]
def sub(p, a, b, label):
 d = [p('subs', a[0], b[0])] + [p('sbcs', a[k], b[k]) for k in (1, 2, 3)]
 fix(p, d, ADDP, label, 'cs')
 return d
def allocate(p, pins):
 first, last, uses = {}, {}, {}
 for i, (op, ds, src, sh, imm) in enumerate(p.ins):
  for x in src + ds:
   if x.startswith('%'):
    first.setdefault(x, i)
    last[x] = i
   elif x != 'xzr':
    uses.setdefault(x, []).append(i)
 for v in pins:
  last[v] = len(p.ins)
 phys, busy = {}, {}
 for i, (op, ds, src, sh, imm) in enumerate(p.ins):
  for x in src:
   if x.startswith('%') and last[x] == i:
    busy.pop(phys[x], None)
  for x in ds:
   if not x.startswith('%') or x in phys:
    if x in phys and last[x] > i:
     busy[phys[x]] = x
    continue
   free = lambda c: (c not in busy and not any(i <= u <= last[x] for u in uses.get(c, []))
                              and not any(first[v] <= last[x] and i <= last[v] for v, r in pins.items() if r == c))
   r = pins.get(x) or next(c for c in POOL if free(c))
   assert r not in busy and not any(i < u <= last[x] for u in uses.get(r, [])), (x, r)
   phys[x] = r
   if last[x] > i:
    busy[r] = x
 R = lambda t: re.sub(r'%v\d+', lambda m: phys[m.group(0)], t)
 out = []
 for op, ds, src, sh, imm in p.ins:
  if op == 'raw':
   out.append(imm)
  elif op in ('ldp', 'stp', 'ldr'):
   out.append('%s %s,%s' % (op, ','.join(R(x) for x in ds + src[:2 if op == 'stp' else 0]), R(imm)))
  else:
   out.append('%s %s' % (op, ','.join(R(x) for x in ds + src + [imm] if x)) + (',%s #32' % sh if sh else ''))
 return '\\n'.join(re.sub(r'(L[pq]\d+)(:?)$', r'\1" S "_%=\2', l) for l in out) + '\\n'
def forward():
 p = P()
 ad = p.m('ldr', '[x23],#32', (p.v(),))[0]
 dy = sub(p, p.ldp('[%s,#32]' % ad) + p.ldp('[%s,#48]' % ad), p.ldp('[x24,#32]') + p.ldp('[x24,#48]'), 'Lq13')
 B1 = row(p, RR)
 B2 = row(p, B1)
 q = product(p, RR, B1, B2, lambda i: dy[i])
 fix(p, q, SUBP, 'Lq14', 'cc')
 p.m('stp', '[x25,#32]', src=q[:2])
 p.m('stp', '[x25,#48]', src=q[2:])
 xa = p.ldp('[%s]' % ad) + p.ldp('[%s,#16]' % ad)
 hi = p.ldp('[x24,#16]')
 D = sub(p, p.ldp('[x24],#64') + hi, xa, 'Lq11')
 p.m('stp', '[x25,#16]', src=D[2:])
 p.m('stp', '[x25],#64', src=D[:2])
 o = product(p, RR, B1, B2, lambda i: D[i])
 fix(p, o, SUBP, 'Lq12', 'cc')
 return allocate(p, dict(zip(o, RR)))
def backward():
 p = P()
 B1 = row(p, RR)
 B2 = row(p, B1)
 sc = {}
 def lazy(mems):
  def a(i):
   if i not in sc:
    sc[i & 2], sc[i | 1] = p.ldp(mems[i // 2])
   return sc[i]
  return a
 lam = product(p, RR, B1, B2, lazy(['[x23,#32]', '[x23,#48]']))
 fix(p, lam, SUBP, 'Lq4', 'cc')
 sc.clear()
 o = product(p, RR, B1, B2, lazy(['[x23],#-64', '[x23,#80]']))
 fix(p, o, SUBP, 'Lq2', 'cc')
 p.raw('ldr x7,[x26]')
 pins = dict(zip(lam, ['x0', 'x2', 'x18', 'x30']))
 pins.update(zip(o, RR))
 return allocate(p, pins)
def once(t, old, new):
 assert t.count(old) == 1, old[:50]
 return t.replace(old, new)
src, dst = sys.argv[1], sys.argv[2]
t = open(src + '/bloop_branchy.h').read()
i = t.index('#define BLOOP_FIELD_S(S) "') + 26
end = 'add x27,x27,#5\\nLq4" S "_%=:\\n'
assert t[i:].startswith('ldp x7,x8,[x23,#32]\\n') and t.count(end) == 1
t = once(t[:i] + backward() + t[t.index(end) + len(end):], '(y3 fold)', '(y3 fold), then gen49.py (shared rows)')
open(dst + '/bloop_branchy.h', 'w').write(t)
f = open(src + '/floop.h').read()
f = re.sub(r'#define FLOOP_S\(S\) ".*"\n', lambda m: '#define FLOOP_S(S) "' + forward() + '"\n', f)
f = once(f, '(paired reduction words): forward prefix products, PRE in x19..x22',
         '(paired reduction words), then gen49.py (shared rows): forward prefix products PRE in x19..x22 and '
         'Q = dy*PRE in the pre slot')
open(dst + '/floop.h', 'w').write(f)
print('ok')
