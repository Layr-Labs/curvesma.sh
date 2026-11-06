"""SMX audit (proof.md s15): restart routine [b]Q_a as one AArch64 listing from the shipped kernel's Montgomery
bodies (MULNR, SQR, BY1), executed end to end by a64.py against integer arithmetic; SMX is the worst-path count."""
import json, random
from pathlib import Path
P = int('1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab', 16)
RA = 1 << 384; N0 = (-pow(P, -1, RA)) & ((1 << 64) - 1); RM = RA % P
T = ['x%d' % i for i in (11, 12, 13, 14, 15, 16)]
U = ['x%d' % i for i in (19, 20, 21, 22, 23, 24)]
R = ['x%d' % i for i in (1, 3, 4, 25, 26, 27)]
PM = ['x%d' % i for i in (5, 6, 7, 8, 9, 10)]
NSLOT = 16   # 0 X, 1 Y, 2 Z, 3 Qx, 4 Qy, 5..15 temporaries
def slot(k):
 assert 0 <= k < NSLOT
 return ('sp', 48 * k) if k < 10 else ('x30', 48 * (k - 10))
SAVE, STATE = 288, 368   # x30-relative: 5 reg pairs; scalar words,
PTRS, LR, FRAME = 1024, 1048, 1056   # sp-relative: pointers; saved x30; frame size
def ld(regs, k):
 b, o = slot(k); return ['ldp %s,%s,[%s,#%d]' % (regs[2 * i], regs[2 * i + 1], b, o + 16 * i) for i in range(3)]
def st(regs, k):
 b, o = slot(k); return ['stp %s,%s,[%s,#%d]' % (regs[2 * i], regs[2 * i + 1], b, o + 16 * i) for i in range(3)]
S6 = ['subs %s,%s,%s' % (R[0], T[0], PM[0])] + ['sbcs %s,%s,%s' % (R[i], T[i], PM[i]) for i in range(1, 6)]
CL = ['csel %s,%s,%s,lo' % (T[i], T[i], R[i]) for i in range(6)]
C1 = S6 + ['sbcs xzr,x0,xzr'] + CL   # subtract p once when (carry, T) >= p
C2 = S6 + CL   # canonicalise a value below 2p
def bodies():
 o = Path(__file__).resolve().parent
 if (o / 'files').is_dir(): o = o / 'files'
 tx = (o / 'measure/side_lanes.S').read_text()
 def rg(a, b):
  s = tx[tx.index(a):tx.index(b, tx.index(a))]
  return [l.split('//')[0].strip() for l in s.splitlines() if l.split('//')[0].strip() and not l.split('//')[0].strip().startswith('.')]
 mn = rg('// ---- inlined MULNR', '// ---- end MULNR'); sq = rg('// ---- inlined SQR', '// ---- end SQR'); m1 = rg('// ---- inlined BY1', '// ---- end BY1')
 assert (len(mn), len(sq), len(m1)) == (295, 91, 150)
 pre = [l.split('//')[0].strip() for l in tx[:tx.index('// ---- inlined MULNR')].splitlines() if l.split('//')[0].strip()][-1]
 assert pre == 'ldp\tx17,x29,[x2]'
 ha = [' '.join(l.split('//')[0].split()) for l in tx[tx.index('// ---- end BY1'):].splitlines()[1:] if l.split('//')[0].strip()][:6]
 assert [x.split()[0] for x in ha] == ['adds', 'adcs', 'adcs', 'adcs', 'adcs', 'adc']
 ST, US = {'adds', 'adcs', 'subs', 'sbcs'}, {'adcs', 'adc', 'sbcs', 'sbc', 'csel'}
 for c in (mn, sq, m1):
  f = False
  for ln in c:
   w = ln.split()[0]
   if w in ST: f = True
   if w in US: assert f
   assert 'sp' not in ln.split()[1:] and 'x30' not in ln and 'x18' not in ln.split(None, 1)[1].split(',')[0] and not ln.startswith(('stp', 'str'))
 return mn, sq, m1, ha
MN, SQ, M1, HA = bodies()
def Mc(a, b, r):
 bb, ob = slot(b)
 return ['add x2,%s,#%d' % (bb, ob)] + ld(T, a) + ['ldp x17,x29,[x2]'] + MN + C2 + st(T, r)
def Sc(a, r): return ld(T, a) + SQ + M1 + HA + C2 + st(T, r)
def a1(a, b, r):   # canonical a + b mod p
 return ld(T, a) + ld(U, b) + ['adds %s,%s,%s' % (T[0], T[0], U[0])] + ['adcs %s,%s,%s' % (T[i], T[i], U[i]) for i in range(1, 6)] + ['adc x0,xzr,xzr'] + C1 + st(T, r)
def sb(a, b, r):   # canonical a - b mod p
 return ld(T, a) + ld(U, b) + ['subs %s,%s,%s' % (T[0], T[0], U[0])] + ['sbcs %s,%s,%s' % (T[i], T[i], U[i]) for i in range(1, 6)] + ['adc x0,xzr,xzr'] + \
  ['adds %s,%s,%s' % (R[0], T[0], PM[0])] + ['adcs %s,%s,%s' % (R[i], T[i], PM[i]) for i in range(1, 5)] + ['adc %s,%s,%s' % (R[5], T[5], PM[5]), 'cmp x0,#1'] + ['csel %s,%s,%s,eq' % (T[i], T[i], R[i]) for i in range(6)] + st(T, r)
def xc(a, r):   # canonical 2a mod p
 return ld(T, a) + ['adds %s,%s,%s' % (T[0], T[0], T[0])] + ['adcs %s,%s,%s' % (T[i], T[i], T[i]) for i in range(1, 6)] + ['adc x0,xzr,xzr'] + C1 + st(T, r)
def x3(a, r, t): return xc(a, t) + a1(t, a, r)
def x4(a, r, t): return xc(a, t) + xc(t, r)
def x8(a, r, t): return xc(a, t) + xc(t, r) + xc(r, r)
assert [len(x) for x in (Mc(0, 1, 2), Sc(0, 1), a1(0, 1, 2), sb(0, 1, 2), xc(0, 1), x3(0, 1, 2), x4(0, 1, 2), x8(0, 1, 2))] == [315, 265, 29, 29, 26, 55, 52, 78]
def DBL():   # dbl-2007-bl on (X, Y, Z) = slots 0, 1, 2;
 return (Sc(0, 5) + Sc(1, 6) + Sc(6, 7) +   # A = X^2, B = Y^2, C = B^2
   a1(0, 6, 8) + Sc(8, 9) + sb(9, 5, 9) + sb(9, 7, 9) + xc(9, 10) +   # D = 2((X+B)^2 - A - C)
   x3(5, 11, 12) + Sc(11, 12) +   # E = 3A, F = E^2
   Mc(1, 2, 15) + xc(15, 2) +   # Z' = 2 Y Z (old Y, old Z)
   xc(10, 13) + sb(12, 13, 0) +   # X' = F - 2D
   sb(10, 0, 13) + Mc(11, 13, 13) + x8(7, 14, 15) + sb(13, 14, 1))   # Y' = E (D - X') - 8C
def ADD():   # madd-2007-bl: (X, Y, Z) += (Qx, Qy) = slots 3, 4
 return (Sc(2, 5) + Mc(3, 5, 6) + Mc(2, 5, 7) + Mc(4, 7, 7) +  # Z1Z1, U2 = Qx Z1Z1, S2 = Qy Z Z1Z1
   sb(6, 0, 8) + Sc(8, 9) + x4(9, 10, 11) + Mc(8, 10, 11) +   # H, HH, I = 4HH, J = H I
   sb(7, 1, 12) + xc(12, 12) + Mc(0, 10, 13) +   # r = 2 (S2 - Y), V = X I
   a1(2, 8, 2) + Sc(2, 2) + sb(2, 5, 2) + sb(2, 9, 2) +   # Z' = (Z + H)^2 - Z1Z1 - HH
   Sc(12, 14) + sb(14, 11, 14) + xc(13, 15) + sb(14, 15, 0) +   # X' = r^2 - J - 2V
   sb(13, 0, 6) + Mc(12, 6, 6) + Mc(1, 11, 14) + xc(14, 14) + sb(6, 14, 1))   # Y' = r (V - X') - 2 Y J
DBL_OP, ADD_OP = len(DBL()), len(ADD())
E = P - 2
FERMAT = []
for i in range(E.bit_length() - 2, -1, -1):
 FERMAT += Sc(2, 2)
 if (E >> i) & 1: FERMAT += Mc(2, 5, 2)
FERMAT_OP = len(FERMAT)
CB = 0x40000   # constants: p (6 words), n0, R mod p (6 words)
def restart_listing():
 L = ['restart:', 'sub sp,sp,#%d' % FRAME, 'str x30,[sp,#%d]' % LR, 'add x30,sp,#512']
 L += ['stp x%d,x%d,[x30,#%d]' % (19 + 2 * i, 20 + 2 * i, SAVE + 16 * i) for i in range(5)]   # callee-saved x19..x28
 L += ['str x0,[sp,#%d]' % PTRS, 'str x1,[sp,#%d]' % (PTRS + 8), 'str x3,[sp,#%d]' % (PTRS + 16)]   # Q_a, scalar, restart slot
 L += ['mov x4,#%d' % CB, 'add x4,x4,#0', 'ldp x5,x6,[x4]', 'ldp x7,x8,[x4,#16]', 'ldp x9,x10,[x4,#32]', 'ldr x18,[x4,#48]']   # p->x5..x10, n0->x18 (adrp/add)
 L += ['ldp x11,x12,[x4,#56]', 'ldp x13,x14,[x4,#72]', 'ldp x15,x16,[x4,#88]'] + st(T, 2)   # Z = R mod p
 L += ['ldp %s,%s,[x0,#%d]' % (T[2 * i], T[2 * i + 1], 16 * i) for i in range(3)] + st(T, 0) + st(T, 3)   # X = Qx (Montgomery residues)
 L += ['ldp %s,%s,[x0,#%d]' % (T[2 * i], T[2 * i + 1], 48 + 16 * i) for i in range(3)] + st(T, 1) + st(T, 4)   # Y = Qy
 L += ['ldp x20,x19,[x1]', 'mov x22,#%d' % (1 << 51), 'umulh x23,x20,x22', 'mul x20,x20,x22', 'mul x19,x19,x22', 'add x19,x19,x23', 'mov x21,#76']   # b << 51: bit 76 at bit 63 of x19
 L += ['.Lscan:', 'tbnz x19,#63,.Lfound', 'adds x20,x20,x20', 'adc x19,x19,x19', 'subs x21,x21,#1', 'b.ne .Lscan']   # first set bit (b >= 1)
 L += ['.Lfound:', 'adds x20,x20,x20', 'adc x19,x19,x19', 'stp x19,x20,[x30,#%d]' % STATE, 'str x21,[x30,#%d]' % (STATE + 16), 'cbz x21,.Linv']   # drop leading bit; x21 = remaining
 L += ['.Lloop:'] + DBL()
 L += ['ldp x19,x20,[x30,#%d]' % STATE, 'ldr x21,[x30,#%d]' % (STATE + 16), 'mov x22,x19', 'adds x20,x20,x20', 'adc x19,x19,x19', 'subs x21,x21,#1',
   'stp x19,x20,[x30,#%d]' % STATE, 'str x21,[x30,#%d]' % (STATE + 16), 'tbz x22,#63,.Lskip']
 L += ADD() + ['.Lskip:', 'ldr x21,[x30,#%d]' % (STATE + 16), 'cbnz x21,.Lloop']
 L += ['.Linv:'] + ld(T, 2) + st(T, 5) + FERMAT   # slot 5 = Z; slot 2 <- Z^(p-2)
 L += Sc(2, 6) + Mc(0, 6, 0) + Mc(1, 6, 1) + Mc(1, 2, 1)   # x = X Z^-2, y = Y Z^-3
 L += ['ldr x0,[sp,#%d]' % (PTRS + 16)] + ld(T, 0) + ['stp %s,%s,[x0,#%d]' % (T[2 * i], T[2 * i + 1], 16 * i) for i in range(3)] + ld(T, 1) + ['stp %s,%s,[x0,#%d]' % (T[2 * i], T[2 * i + 1], 48 + 16 * i) for i in range(3)]
 L += ['ldp x%d,x%d,[x30,#%d]' % (19 + 2 * i, 20 + 2 * i, SAVE + 16 * i) for i in range(5)] + ['ldr x30,[sp,#%d]' % LR, 'add sp,sp,#%d' % FRAME, 'ret']
 return L
LISTING = restart_listing()
ITER_CTL = 11   # per ladder iteration outside DBL/ADD
SCAN_MAX = 76 * 5   # b = 1: 76 shifts of 5 instructions beyond
def _count(sec): return len([l for l in sec if not l.endswith(':')])
ENTRY_OP = _count(LISTING[:LISTING.index('.Lscan:')]); EXIT_OP = _count(LISTING[LISTING.index('.Linv:') + 1:]) - FERMAT_OP - 6
FOUND_OP, COPY_OP = 5, 6   # drop leading bit, store state, test zero
SMX = ENTRY_OP + SCAN_MAX + FOUND_OP + 76 * (DBL_OP + ITER_CTL + ADD_OP) + COPY_OP + FERMAT_OP + EXIT_OP   # scan maximum and ladder maximum summed
assert (DBL_OP, ADD_OP, FERMAT_OP, ENTRY_OP, EXIT_OP, SMX) == (2340, 3656, 172520, 48, 1231, 630722)
def aff_mul(k, pt):
 def add(p1, q1):
  if p1 is None: return q1
  if q1 is None: return p1
  (x1, y1), (x2, y2) = p1, q1
  if x1 == x2 and (y1 + y2) % P == 0: return None
  l = (3 * x1 * x1) * pow(2 * y1, -1, P) % P if x1 == x2 else ((y2 - y1) * pow(x2 - x1, -1, P)) % P
  x3_ = (l * l - x1 - x2) % P
  return (x3_, (l * (x1 - x3_) - y1) % P)
 o = None
 while k:
  if k & 1: o = add(o, pt)
  pt = add(pt, pt); k >>= 1
 return o
def audit_smx(trials=1):
 import a64
 o = Path(__file__).resolve().parent
 if (o / 'files').is_dir(): o = o / 'files'
 r = random.Random(20261103); M64 = (1 << 64) - 1
 prog, labels = a64.a64_load('\n'.join(l if l.endswith(':') else '\t' + l for l in LISTING))
 assert _count(LISTING) == len(prog)
 qm = json.loads((o / 'params.json').read_text()); G = (int(qm['gx'], 16), int(qm['gy'], 16))
 assert (G[1] ** 2 - G[0] ** 3 - 4) % P == 0
 def put(m, a, v):
  for i in range(6): m[a + 8 * i] = (v >> (64 * i)) & M64
 def get(m, a): return sum(m[a + 8 * i] << (64 * i) for i in range(6))
 QA, SC, OUT, SP = 0x30000, 0x31000, 0x32000, 0x20000
 counts = {}
 cases = [(1 << 77) - 1, 1, 2] + [r.getrandbits(77) or 1 for _ in range(trials)]
 for b in cases:
  Qa = aff_mul(r.randrange(1, 1 << 20), G)
  m = {}
  for i, v in enumerate([(P >> (64 * i)) & M64 for i in range(6)] + [N0] + [(RM >> (64 * i)) & M64 for i in range(6)]): m[CB + 8 * i] = v
  put(m, QA, Qa[0] * RA % P); put(m, QA + 48, Qa[1] * RA % P)   # residues, as in gap table
  m[SC] = b & M64; m[SC + 8] = b >> 64
  for i in range(12): m[OUT + 8 * i] = r.getrandbits(64)
  regs = {'x%d' % i: r.getrandbits(64) for i in range(31)}; regs.update(x0=QA, x1=SC, x3=OUT, sp=SP)
  mem = dict(m); ro = {}
  c = a64.a64_run(prog, labels, 'restart', regs, mem, limit=2 * SMX, out_regs=ro)
  want = aff_mul(b, Qa)
  assert get(mem, OUT) == want[0] * RA % P and get(mem, OUT + 48) == want[1] * RA % P, b   # canonical Montgomery residues
  for a_, v in m.items(): assert mem[a_] == v or OUT <= a_ < OUT + 96, hex(a_)   # inputs and constants untouched
  for reg in ['sp', 'x30'] + ['x%d' % i for i in range(19, 29)]: assert ro[reg] == regs[reg], reg
  assert c <= SMX, (b, c, SMX)
  counts[b] = c
 assert counts[(1 << 77) - 1] == SMX - SCAN_MAX + 1 and counts[1] == SMX - 76 * (DBL_OP + ITER_CTL + ADD_OP)   # all-ones: one scan test; b = 1: full scan,
 return {'SMX': SMX, 'DBL': DBL_OP, 'ADD': ADD_OP, 'Fermat': FERMAT_OP, 'entry': ENTRY_OP, 'exit': EXIT_OP, 'iteration_control': ITER_CTL,
   'scan_max': SCAN_MAX, 'listing_instructions': _count(LISTING), 'emulated_scalars': len(cases), 'counts': {str(k)[:8]: v for k, v in counts.items()},
   'scope': 'The generated restart listing executed end to end by a64.py on Montgomery-form inputs; worst path = all-ones scalar plus the longest scan.'}
if __name__ == '__main__': print(json.dumps(audit_smx()))
